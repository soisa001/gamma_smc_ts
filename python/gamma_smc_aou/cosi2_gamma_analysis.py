"""Unit-level Gamma-SMC inference for the focused CoSi2 comparison.

The inferential unit is one independently simulated CoSi2 locus. Pair,
position, and threshold measurements must be reduced upstream to one score per
unit, statistic, window, and threshold before this module is called. The
natural-neutral reference is matched by demography only; final and sampled
allele frequency are recorded but are not used to condition the neutral bank.

With one selected unit, the reported AUC is only its neutral-bank midrank. It
is deliberately labelled descriptive and is not accompanied by a bootstrap
confidence interval over the selected model.
"""

from __future__ import annotations

import os
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.stats import beta

matplotlib.use("Agg")


SCHEMA_VERSION = "gamma-smc.cosi2-gamma-analysis/v1"
PANEL_HAPLOTYPES = 200
TMRCA_THRESHOLDS_YEARS = (
    1_000,
    4_500,
    10_000,
    20_000,
    30_000,
    40_000,
    50_000,
)
SUPPORTED_STATISTICS = (
    "paper_frac_posterior_mean_below",
    "mean_p_tmrca_lt",
)
SUPPORTED_WINDOWS = ("local_100kb", "focal")
PRIMARY_WINDOW = "local_100kb"
SIMULATION_CLASSES = ("neutral", "selected")
INPUT_COLUMNS = (
    "unit_id",
    "demography",
    "simulation_class",
    "selection_coefficient",
    "seed",
    "generation_time_years",
    "statistic",
    "window",
    "threshold_years",
    "score",
    "final_population_af",
    "sample_alt_count",
    "sample_af",
)
UNIT_METADATA_COLUMNS = (
    "unit_id",
    "demography",
    "simulation_class",
    "selection_coefficient",
    "seed",
    "generation_time_years",
    "final_population_af",
    "sample_alt_count",
    "sample_af",
)
SCORE_KEY = ("unit_id", "statistic", "window", "threshold_years")
_EPSILON = 1e-14

_STATISTIC_LABELS = {
    "paper_frac_posterior_mean_below": (
        "Paper-compatible fraction\n(posterior mean TMRCA < x)"
    ),
    "mean_p_tmrca_lt": "Mean posterior P(TMRCA < x)",
}
_WINDOW_COLORS = {"local_100kb": "#2171b5", "focal": "#d95f0e"}
_WINDOW_LABELS = {"local_100kb": "+/-100 kb (primary)", "focal": "Focal"}


def _finite_numeric(frame: pd.DataFrame, column: str) -> pd.Series:
    values = pd.to_numeric(frame[column], errors="coerce")
    if not np.all(np.isfinite(values.to_numpy(dtype=float))):
        raise ValueError(f"{column} must contain only finite numeric values")
    return values


def _integer_numeric(frame: pd.DataFrame, column: str) -> pd.Series:
    values = _finite_numeric(frame, column)
    numeric = values.to_numpy(dtype=float)
    if not np.array_equal(numeric, np.floor(numeric)):
        raise ValueError(f"{column} must contain only integers")
    return values.astype(np.int64)


def validate_scores(scores: pd.DataFrame) -> pd.DataFrame:
    """Validate and deterministically order complete CoSi2 score curves."""

    if not isinstance(scores, pd.DataFrame):
        raise TypeError("scores must be a pandas DataFrame")
    missing = sorted(set(INPUT_COLUMNS).difference(scores.columns))
    if missing:
        raise ValueError(f"score table is missing columns: {', '.join(missing)}")
    if scores.empty:
        raise ValueError("score table is empty")

    result = scores.loc[:, INPUT_COLUMNS].copy()
    text_columns = (
        "unit_id",
        "demography",
        "simulation_class",
        "statistic",
        "window",
    )
    for column in text_columns:
        if result[column].isna().any():
            raise ValueError(f"{column} contains missing values")
        result[column] = result[column].astype(str).str.strip()
        if (result[column] == "").any():
            raise ValueError(f"{column} contains empty values")

    result["selection_coefficient"] = _finite_numeric(
        result, "selection_coefficient"
    ).astype(float)
    result["seed"] = _integer_numeric(result, "seed")
    result["generation_time_years"] = _finite_numeric(
        result, "generation_time_years"
    ).astype(float)
    result["threshold_years"] = _integer_numeric(result, "threshold_years")
    result["score"] = _finite_numeric(result, "score").astype(float)
    result["final_population_af"] = _finite_numeric(
        result, "final_population_af"
    ).astype(float)
    result["sample_alt_count"] = _integer_numeric(result, "sample_alt_count")
    result["sample_af"] = _finite_numeric(result, "sample_af").astype(float)

    if set(result["simulation_class"]) != set(SIMULATION_CLASSES):
        raise ValueError(
            "simulation_class must contain exactly neutral and selected units"
        )
    if set(result["statistic"]) != set(SUPPORTED_STATISTICS):
        raise ValueError(
            "statistic must contain exactly the two supported score definitions"
        )
    if set(result["window"]) != set(SUPPORTED_WINDOWS):
        raise ValueError("window must contain exactly local_100kb and focal")
    if set(result["threshold_years"].astype(int)) != set(TMRCA_THRESHOLDS_YEARS):
        raise ValueError("threshold_years does not match the fixed seven thresholds")

    bounded_columns = ("score", "final_population_af", "sample_af")
    for column in bounded_columns:
        if ((result[column] < 0.0) | (result[column] > 1.0)).any():
            raise ValueError(f"{column} must lie in [0, 1]")
    if (result["final_population_af"] <= 0.0).any():
        raise ValueError(
            "final_population_af must be positive for every survival-conditioned unit"
        )
    if (
        (result["sample_alt_count"] < 0)
        | (result["sample_alt_count"] > PANEL_HAPLOTYPES)
    ).any():
        raise ValueError(f"sample_alt_count must lie in [0, {PANEL_HAPLOTYPES}]")
    expected_sample_af = result["sample_alt_count"].to_numpy(dtype=float) / float(
        PANEL_HAPLOTYPES
    )
    if not np.allclose(
        result["sample_af"].to_numpy(dtype=float),
        expected_sample_af,
        rtol=0.0,
        atol=1e-12,
    ):
        raise ValueError("sample_af disagrees with sample_alt_count / 200")

    neutral_coefficients = result.loc[
        result["simulation_class"] == "neutral", "selection_coefficient"
    ].to_numpy(dtype=float)
    selected_coefficients = result.loc[
        result["simulation_class"] == "selected", "selection_coefficient"
    ].to_numpy(dtype=float)
    if not np.allclose(neutral_coefficients, 0.0, rtol=0.0, atol=1e-15):
        raise ValueError("neutral units must have selection_coefficient 0")
    if np.any(selected_coefficients <= 0.0):
        raise ValueError("selected units must have a positive selection_coefficient")
    if (result["seed"] < 0).any():
        raise ValueError("seed must be nonnegative")
    if (result["generation_time_years"] <= 0.0).any():
        raise ValueError("generation_time_years must be positive")

    if result.duplicated(list(SCORE_KEY), keep=False).any():
        raise ValueError("score table contains duplicate unit/statistic/window rows")

    metadata = result.loc[:, UNIT_METADATA_COLUMNS].drop_duplicates()
    if metadata["unit_id"].duplicated(keep=False).any():
        raise ValueError("unit metadata changes across score rows")
    if metadata.duplicated(["demography", "seed"], keep=False).any():
        raise ValueError("independent units must have unique seeds within demography")

    expected_curves = {
        (statistic, window, int(threshold))
        for statistic in SUPPORTED_STATISTICS
        for window in SUPPORTED_WINDOWS
        for threshold in TMRCA_THRESHOLDS_YEARS
    }
    for unit_id, unit in result.groupby("unit_id", sort=False):
        observed = set(
            zip(
                unit["statistic"],
                unit["window"],
                unit["threshold_years"].astype(int),
                strict=True,
            )
        )
        if observed != expected_curves:
            raise ValueError(f"unit {unit_id} lacks the exact complete score grid")
        for labels, curve in unit.groupby(["statistic", "window"], sort=False):
            ordered = curve.sort_values("threshold_years")
            if tuple(ordered["threshold_years"].astype(int)) != (
                TMRCA_THRESHOLDS_YEARS
            ):
                raise ValueError(f"unit {unit_id} {labels} lacks a complete curve")
            if np.any(np.diff(ordered["score"].to_numpy(dtype=float)) < -1e-12):
                raise ValueError(f"unit {unit_id} {labels} has a nonmonotone curve")

    for demography, group in metadata.groupby("demography", sort=False):
        if set(group["simulation_class"]) != set(SIMULATION_CLASSES):
            raise ValueError(
                f"demography {demography} lacks a neutral or selected unit"
            )
        if group["generation_time_years"].nunique(dropna=False) != 1:
            raise ValueError(
                f"generation_time_years changes within demography {demography}"
            )

    return result.sort_values(
        ["demography", "unit_id", "statistic", "window", "threshold_years"],
        kind="mergesort",
    ).reset_index(drop=True)


def clopper_pearson_interval(
    successes: int, trials: int, *, confidence_level: float = 0.95
) -> tuple[float, float]:
    """Return an exact equal-tailed binomial confidence interval."""

    if isinstance(successes, (bool, np.bool_)) or not isinstance(
        successes, (int, np.integer)
    ):
        raise TypeError("successes must be an integer")
    if isinstance(trials, (bool, np.bool_)) or not isinstance(
        trials, (int, np.integer)
    ):
        raise TypeError("trials must be an integer")
    successes = int(successes)
    trials = int(trials)
    if trials < 1 or successes < 0 or successes > trials:
        raise ValueError("require 0 <= successes <= trials and trials >= 1")
    if not np.isfinite(confidence_level) or not 0.0 < confidence_level < 1.0:
        raise ValueError("confidence_level must lie strictly between zero and one")

    alpha = 1.0 - float(confidence_level)
    low = (
        0.0
        if successes == 0
        else float(beta.ppf(alpha / 2.0, successes, trials - successes + 1))
    )
    high = (
        1.0
        if successes == trials
        else float(beta.ppf(1.0 - alpha / 2.0, successes + 1, trials - successes))
    )
    return low, high


def _tail_summary(observed: float, neutral: np.ndarray) -> dict[str, Any]:
    upper_count = int(np.count_nonzero(neutral >= observed - _EPSILON))
    lower_count = int(np.count_nonzero(neutral <= observed + _EPSILON))
    strictly_less = int(np.count_nonzero(neutral < observed - _EPSILON))
    strictly_greater = int(np.count_nonzero(neutral > observed + _EPSILON))
    tie_count = int(len(neutral) - strictly_less - strictly_greater)
    denominator = len(neutral) + 1
    p_upper = float((1 + upper_count) / denominator)
    p_lower = float((1 + lower_count) / denominator)
    p_two_sided = float(min(1.0, 2.0 * min(p_upper, p_lower)))
    upper_low, upper_high = clopper_pearson_interval(upper_count, len(neutral))
    lower_low, lower_high = clopper_pearson_interval(lower_count, len(neutral))
    return {
        "neutral_upper_exceedance_count": upper_count,
        "neutral_lower_exceedance_count": lower_count,
        "neutral_strictly_less_count": strictly_less,
        "neutral_tie_count": tie_count,
        "mc_p_upper": p_upper,
        "mc_p_lower": p_lower,
        "mc_p_two_sided": p_two_sided,
        "upper_tail_cp95_low": upper_low,
        "upper_tail_cp95_high": upper_high,
        "lower_tail_cp95_low": lower_low,
        "lower_tail_cp95_high": lower_high,
        "neutral_midrank_descriptive_auc": float(
            (strictly_less + 0.5 * tie_count) / len(neutral)
        ),
    }


def pointwise_pvalues(scores: pd.DataFrame) -> pd.DataFrame:
    """Compare every selected unit with its demography-matched neutral bank."""

    frame = validate_scores(scores)
    group_columns = ["demography", "statistic", "window", "threshold_years"]
    rows: list[dict[str, Any]] = []
    for labels, group in frame.groupby(group_columns, sort=True, dropna=False):
        label = dict(zip(group_columns, labels, strict=True))
        neutral = group.loc[group["simulation_class"] == "neutral", "score"].to_numpy(
            dtype=float
        )
        if not len(neutral):
            raise ValueError(f"matched neutral bank is empty for {label}")
        quantiles = np.quantile(neutral, [0.025, 0.05, 0.5, 0.95, 0.975])
        selected = group[group["simulation_class"] == "selected"]
        for record in selected.to_dict(orient="records"):
            observed = float(record["score"])
            rows.append(
                {
                    "schema": SCHEMA_VERSION,
                    **{column: record[column] for column in UNIT_METADATA_COLUMNS},
                    **label,
                    "threshold_generations": (
                        float(record["threshold_years"])
                        / float(record["generation_time_years"])
                    ),
                    "observed_score": observed,
                    "n_neutral": len(neutral),
                    "neutral_mean": float(np.mean(neutral)),
                    "neutral_sd": (
                        float(np.std(neutral, ddof=1)) if len(neutral) > 1 else np.nan
                    ),
                    "neutral_q025": float(quantiles[0]),
                    "neutral_q05": float(quantiles[1]),
                    "neutral_median": float(quantiles[2]),
                    "neutral_q95": float(quantiles[3]),
                    "neutral_q975": float(quantiles[4]),
                    **_tail_summary(observed, neutral),
                    "minimum_attainable_p": float(1 / (len(neutral) + 1)),
                    "confidence_interval_method": (
                        "exact_clopper_pearson_for_unadjusted_null_tail_probability"
                    ),
                    "pvalue_method": "neutral_tail_plus_one_conservative_ties",
                    "neutral_bank": (
                        "same_demography_statistic_window_threshold_"
                        "natural_survivors_no_af_conditioning"
                    ),
                    "primary_window": str(record["window"]) == PRIMARY_WINDOW,
                    "score_direction": "larger_is_more_selected_no_flipping",
                    "auc_interpretation": (
                        "descriptive_selected_unit_neutral_midrank_not_power"
                    ),
                }
            )
    return (
        pd.DataFrame(rows)
        .sort_values(
            [
                "demography",
                "unit_id",
                "statistic",
                "window",
                "threshold_years",
            ],
            kind="mergesort",
        )
        .reset_index(drop=True)
    )


def _exchangeable_minp(
    neutral_profiles: np.ndarray, observed_profile: np.ndarray
) -> dict[str, Any]:
    if neutral_profiles.ndim != 2:
        raise ValueError("neutral_profiles must be a two-dimensional matrix")
    if observed_profile.shape != (neutral_profiles.shape[1],):
        raise ValueError("observed profile does not match neutral thresholds")
    if not len(neutral_profiles):
        raise ValueError("neutral profile bank is empty")
    if not np.all(np.isfinite(neutral_profiles)) or not np.all(
        np.isfinite(observed_profile)
    ):
        raise ValueError("minP profiles must be finite")

    pooled = np.vstack((neutral_profiles, observed_profile[None, :]))
    denominator = len(pooled)
    upper = np.empty_like(pooled, dtype=float)
    lower = np.empty_like(pooled, dtype=float)
    for index in range(pooled.shape[1]):
        values = pooled[:, index]
        upper[:, index] = (
            np.count_nonzero(
                values[None, :] >= values[:, None] - _EPSILON,
                axis=1,
            )
            / denominator
        )
        lower[:, index] = (
            np.count_nonzero(
                values[None, :] <= values[:, None] + _EPSILON,
                axis=1,
            )
            / denominator
        )
    two_sided = np.minimum(1.0, 2.0 * np.minimum(upper, lower))

    result: dict[str, Any] = {}
    for tail, marginal in (
        ("upper", upper),
        ("lower", lower),
        ("two_sided", two_sided),
    ):
        minimum = np.min(marginal, axis=1)
        observed_minimum = float(minimum[-1])
        extreme_count = int(
            np.count_nonzero(minimum[:-1] <= observed_minimum + _EPSILON)
        )
        result[f"observed_min_pointwise_p_{tail}"] = observed_minimum
        result[f"threshold_index_at_min_pointwise_p_{tail}"] = int(
            np.argmin(marginal[-1])
        )
        result[f"neutral_minp_as_or_more_extreme_count_{tail}"] = extreme_count
        result[f"omnibus_mc_p_{tail}"] = float((1 + extreme_count) / denominator)
    return result


def minp_omnibus(scores: pd.DataFrame) -> pd.DataFrame:
    """Calibrate seven-threshold minP while retaining whole unit vectors."""

    frame = validate_scores(scores)
    rows: list[dict[str, Any]] = []
    group_columns = ["demography", "statistic", "window"]
    thresholds = tuple(int(value) for value in TMRCA_THRESHOLDS_YEARS)
    for labels, group in frame.groupby(group_columns, sort=True, dropna=False):
        label = dict(zip(group_columns, labels, strict=True))
        profiles = group.pivot(
            index="unit_id", columns="threshold_years", values="score"
        ).reindex(columns=thresholds)
        metadata = (
            group.loc[:, UNIT_METADATA_COLUMNS].drop_duplicates().set_index("unit_id")
        )
        neutral_ids = metadata.index[
            metadata["simulation_class"].astype(str) == "neutral"
        ]
        selected_ids = metadata.index[
            metadata["simulation_class"].astype(str) == "selected"
        ]
        neutral = profiles.loc[neutral_ids].sort_index()
        selected = profiles.loc[selected_ids].sort_index()
        if neutral.empty or selected.empty or neutral.isna().any().any():
            raise ValueError(f"incomplete minP profile bank for {label}")
        if selected.isna().any().any():
            raise ValueError(f"incomplete selected minP profiles for {label}")

        neutral_values = neutral.to_numpy(dtype=float)
        for unit_id, observed_row in selected.iterrows():
            observed = observed_row.to_numpy(dtype=float)
            inference = _exchangeable_minp(neutral_values, observed)
            record = metadata.loc[unit_id]
            output: dict[str, Any] = {
                "schema": SCHEMA_VERSION,
                "unit_id": str(unit_id),
                "demography": str(record["demography"]),
                "simulation_class": str(record["simulation_class"]),
                "selection_coefficient": float(record["selection_coefficient"]),
                "seed": int(record["seed"]),
                "generation_time_years": float(record["generation_time_years"]),
                "statistic": label["statistic"],
                "window": label["window"],
                "final_population_af": float(record["final_population_af"]),
                "sample_alt_count": int(record["sample_alt_count"]),
                "sample_af": float(record["sample_af"]),
                "n_neutral": len(neutral_values),
                "n_thresholds": len(thresholds),
                "thresholds_years": ",".join(map(str, thresholds)),
                "minimum_attainable_p": float(1 / (len(neutral_values) + 1)),
                "primary_window": label["window"] == PRIMARY_WINDOW,
                "omnibus_method": (
                    "exchangeable_whole_unit_minP_over_seven_thresholds"
                ),
                "correlation_handling": ("complete_neutral_threshold_vectors_retained"),
                "score_direction": "larger_is_more_selected_no_flipping",
            }
            for tail in ("upper", "lower", "two_sided"):
                index = int(inference[f"threshold_index_at_min_pointwise_p_{tail}"])
                output[f"observed_min_pointwise_p_{tail}"] = float(
                    inference[f"observed_min_pointwise_p_{tail}"]
                )
                output[f"threshold_years_at_min_pointwise_p_{tail}"] = thresholds[index]
                output[f"observed_score_at_min_pointwise_p_{tail}"] = float(
                    observed[index]
                )
                output[f"neutral_minp_as_or_more_extreme_count_{tail}"] = int(
                    inference[f"neutral_minp_as_or_more_extreme_count_{tail}"]
                )
                output[f"omnibus_mc_p_{tail}"] = float(
                    inference[f"omnibus_mc_p_{tail}"]
                )
            rows.append(output)
    return (
        pd.DataFrame(rows)
        .sort_values(["demography", "unit_id", "statistic", "window"], kind="mergesort")
        .reset_index(drop=True)
    )


def analyze_scores(scores: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Return validated inputs and both prespecified inference tables."""

    validated = validate_scores(scores)
    return {
        "scores": validated,
        "pointwise_pvalues": pointwise_pvalues(validated),
        "minp_omnibus": minp_omnibus(validated),
    }


def _atomic_save_figure_pair(
    figure: Figure, output_dir: Path, stem: str
) -> Mapping[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "png": output_dir / f"{stem}.png",
        "pdf": output_dir / f"{stem}.pdf",
    }
    for extension, path in outputs.items():
        temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
        try:
            metadata = {
                "Creator": "gamma_smc_aou.cosi2_gamma_analysis",
                "CreationDate": None,
                "ModDate": None,
            }
            figure.savefig(
                temporary,
                format=extension,
                dpi=300 if extension == "png" else None,
                facecolor="white",
                metadata=(
                    {"Software": "gamma_smc_aou.cosi2_gamma_analysis"}
                    if extension == "png"
                    else metadata
                ),
            )
            os.replace(temporary, path)
        finally:
            temporary.unlink(missing_ok=True)
    return outputs


def plot_score_profiles(
    scores: pd.DataFrame,
    output_dir: str | Path,
    *,
    stem: str = "cosi2_gamma_selected_vs_neutral_profiles",
) -> Mapping[str, Path]:
    """Write a publication-sized Letter-landscape PNG/PDF profile figure."""

    frame = validate_scores(scores)
    if not stem or Path(stem).name != stem:
        raise ValueError("stem must be a nonempty filename stem")
    destination = Path(output_dir).resolve()
    demographies = tuple(sorted(frame["demography"].unique()))
    n_columns = len(demographies)
    figure = Figure(figsize=(11.0, 8.5))
    axes = figure.subplots(
        len(SUPPORTED_STATISTICS),
        n_columns,
        sharex=True,
        sharey=True,
        squeeze=False,
    )
    x = np.asarray(TMRCA_THRESHOLDS_YEARS, dtype=float) / 1_000.0
    try:
        for row_index, statistic in enumerate(SUPPORTED_STATISTICS):
            for column_index, demography in enumerate(demographies):
                axis = axes[row_index, column_index]
                cell = frame[
                    (frame["demography"] == demography)
                    & (frame["statistic"] == statistic)
                ]
                unit_metadata = cell.loc[
                    :, ["unit_id", "simulation_class"]
                ].drop_duplicates()
                n_neutral = int((unit_metadata["simulation_class"] == "neutral").sum())
                n_selected = int(
                    (unit_metadata["simulation_class"] == "selected").sum()
                )
                for window in SUPPORTED_WINDOWS:
                    color = _WINDOW_COLORS[window]
                    window_frame = cell[cell["window"] == window]
                    profiles = window_frame.pivot(
                        index="unit_id", columns="threshold_years", values="score"
                    ).reindex(columns=TMRCA_THRESHOLDS_YEARS)
                    classes = (
                        window_frame.loc[:, ["unit_id", "simulation_class"]]
                        .drop_duplicates()
                        .set_index("unit_id")["simulation_class"]
                    )
                    neutral = profiles.loc[
                        classes.index[classes.astype(str) == "neutral"]
                    ].to_numpy(dtype=float)
                    selected = profiles.loc[
                        classes.index[classes.astype(str) == "selected"]
                    ].to_numpy(dtype=float)
                    neutral_low, neutral_median, neutral_high = np.quantile(
                        neutral, [0.025, 0.5, 0.975], axis=0
                    )
                    axis.fill_between(
                        x,
                        neutral_low,
                        neutral_high,
                        color=color,
                        alpha=0.13,
                        linewidth=0,
                    )
                    axis.plot(
                        x,
                        neutral_median,
                        color=color,
                        linestyle="--",
                        linewidth=2.2,
                    )
                    axis.plot(
                        x,
                        np.median(selected, axis=0),
                        color=color,
                        linestyle="-",
                        marker="o" if window == PRIMARY_WINDOW else "s",
                        markersize=5.0,
                        linewidth=2.6,
                    )
                axis.set_title(
                    f"{demography} - {_STATISTIC_LABELS[statistic]}",
                    fontsize=15,
                    pad=10,
                )
                axis.set_ylim(0.0, 1.0)
                axis.set_xticks(x)
                axis.set_xticklabels(
                    ["1", "4.5", "10", "20", "30", "40", "50"], fontsize=11
                )
                axis.tick_params(axis="y", labelsize=12)
                axis.grid(True, color="#d9d9d9", linewidth=0.7, alpha=0.8)
                axis.text(
                    0.03,
                    0.95,
                    f"neutral n={n_neutral}; selected n={n_selected}",
                    transform=axis.transAxes,
                    ha="left",
                    va="top",
                    fontsize=11,
                )
                if column_index == 0:
                    axis.set_ylabel("Unit-level score", fontsize=14)
                if row_index == len(SUPPORTED_STATISTICS) - 1:
                    axis.set_xlabel("TMRCA threshold (kya)", fontsize=14)

        legend_handles: list[Any] = [
            Line2D(
                [0],
                [0],
                color=_WINDOW_COLORS[window],
                linewidth=3,
                label=_WINDOW_LABELS[window],
            )
            for window in SUPPORTED_WINDOWS
        ]
        legend_handles.extend(
            [
                Line2D(
                    [0],
                    [0],
                    color="#333333",
                    linestyle="-",
                    linewidth=2.6,
                    label="Selected median",
                ),
                Line2D(
                    [0],
                    [0],
                    color="#333333",
                    linestyle="--",
                    linewidth=2.2,
                    label="Neutral median",
                ),
                Patch(
                    facecolor="#999999",
                    edgecolor="none",
                    alpha=0.2,
                    label="Neutral 95% range",
                ),
            ]
        )
        figure.legend(
            handles=legend_handles,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.985),
            ncol=5,
            frameon=False,
            fontsize=12.5,
        )
        figure.suptitle(
            "CoSi2 Gamma-SMC: selected versus natural-neutral profiles",
            fontsize=20,
            y=0.915,
        )
        figure.subplots_adjust(
            left=0.08,
            right=0.985,
            bottom=0.09,
            top=0.80,
            wspace=0.17,
            hspace=0.34,
        )
        return _atomic_save_figure_pair(figure, destination, stem)
    finally:
        figure.clear()


__all__ = [
    "INPUT_COLUMNS",
    "PANEL_HAPLOTYPES",
    "PRIMARY_WINDOW",
    "SCHEMA_VERSION",
    "SUPPORTED_STATISTICS",
    "SUPPORTED_WINDOWS",
    "TMRCA_THRESHOLDS_YEARS",
    "analyze_scores",
    "clopper_pearson_interval",
    "minp_omnibus",
    "plot_score_profiles",
    "pointwise_pvalues",
    "validate_scores",
]
