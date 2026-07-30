from __future__ import annotations

import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from tempfile import TemporaryDirectory
from time import perf_counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .calibration import bh_fdr
from .selection import run_slim_recent_sweep, within_individual_tmrca_grid


FONTS = {
    "suptitle": 23,
    "title": 20,
    "axis": 18,
    "tick": 15,
    "legend": 15,
    "annotation": 14,
}


def window_centers(sequence_length: int, window_size: int, center: int) -> np.ndarray:
    """Return an evenly spaced position grid containing the selected site."""
    if sequence_length < 1 or window_size < 1:
        raise ValueError("sequence length and window size must be positive")
    left = np.arange(center, -1, -window_size, dtype=float)[::-1]
    right = np.arange(center + window_size, sequence_length, window_size, dtype=float)
    return np.unique(np.r_[left, right, sequence_length - 1])


def calibrate_spatial_windows(
    observed: pd.DataFrame,
    neutral: pd.DataFrame,
    *,
    statistic_column: str = "mean_p_tmrca_lt_threshold",
) -> pd.DataFrame:
    """Calculate pointwise Monte Carlo p-values at matched spatial windows.

    ``statistic_column`` may be the soft mean posterior probability or one of
    the decoder's hard-call fractions (for example ``frac_recent_4500``).
    """
    required = {"position_0based", "n_pairs", statistic_column}
    missing_observed = required.difference(observed.columns)
    missing_neutral = (required | {"replicate"}).difference(neutral.columns)
    if missing_observed:
        raise ValueError(
            "observed profile is missing columns: "
            + ", ".join(sorted(missing_observed))
        )
    if missing_neutral:
        raise ValueError(
            "neutral profiles are missing columns: "
            + ", ".join(sorted(missing_neutral))
        )
    rows = []
    replicate_count = neutral["replicate"].nunique()
    if replicate_count < 1:
        raise ValueError("neutral profiles contain no replicates")
    coverage = neutral.groupby("position_0based")["replicate"].nunique()
    complete_positions = coverage.index[coverage.eq(replicate_count)].to_numpy(
        dtype=float
    )
    matched_observed = observed[
        observed["position_0based"].isin(complete_positions)
    ]
    if matched_observed.empty:
        raise ValueError("no output position is present in every neutral replicate")
    for _, row in matched_observed.sort_values("position_0based").iterrows():
        position = float(row["position_0based"])
        values = neutral.loc[
            np.isclose(neutral["position_0based"], position),
            statistic_column,
        ].to_numpy(dtype=float)
        if len(values) != replicate_count:
            raise ValueError(
                f"position {position:g} has {len(values)} null values; "
                f"expected {replicate_count}"
            )
        observed_value = float(row[statistic_column])
        exceedances = int(np.count_nonzero(values >= observed_value))
        rows.append({
            "position_0based": position,
            "n_pairs": int(row["n_pairs"]),
            "observed_fraction_recent": observed_value,
            "neutral_mean_fraction_recent": float(values.mean()),
            "neutral_median_fraction_recent": float(np.median(values)),
            "neutral_ci95_lower": float(np.quantile(values, 0.025)),
            "neutral_ci95_upper": float(np.quantile(values, 0.975)),
            "neutral_exceedances": exceedances,
            "n_neutral_replicates": int(replicate_count),
            "p_upper": float((1 + exceedances) / (1 + replicate_count)),
        })
    result = pd.DataFrame(rows)
    result["q_bh"] = bh_fdr(result["p_upper"].to_numpy(dtype=float))
    return result


def significant_regions(
    scan: pd.DataFrame,
    *,
    sequence_length: int,
    window_size: int,
    alpha: float = 0.05,
) -> pd.DataFrame:
    """Merge adjacent significant grid windows and report midpoint breakpoints."""
    significant = scan[scan["p_upper"] < alpha].sort_values("position_0based")
    columns = [
        "region_id",
        "start_0based",
        "end_0based_exclusive",
        "start_mb",
        "end_mb",
        "n_windows",
        "minimum_p",
        "maximum_observed_fraction_recent",
    ]
    if significant.empty:
        return pd.DataFrame(columns=columns)
    positions = significant["position_0based"].to_numpy(dtype=float)
    groups = np.cumsum(np.r_[False, np.diff(positions) > window_size * 1.01])
    rows = []
    half = window_size / 2
    for region_id, (_, group) in enumerate(significant.groupby(groups), start=1):
        start = max(0.0, float(group["position_0based"].min() - half))
        end = min(
            float(sequence_length),
            float(group["position_0based"].max() + half),
        )
        rows.append({
            "region_id": region_id,
            "start_0based": int(round(start)),
            "end_0based_exclusive": int(round(end)),
            "start_mb": start / 1e6,
            "end_mb": end / 1e6,
            "n_windows": int(len(group)),
            "minimum_p": float(group["p_upper"].min()),
            "maximum_observed_fraction_recent": float(
                group["observed_fraction_recent"].max()
            ),
        })
    return pd.DataFrame(rows, columns=columns)


def _format_axis(axis, *, ylabel: str | None = None) -> None:
    axis.set_xlabel("Position (Mb)", fontsize=FONTS["axis"])
    if ylabel:
        axis.set_ylabel(ylabel, fontsize=FONTS["axis"])
    axis.tick_params(axis="both", labelsize=FONTS["tick"])
    axis.grid(alpha=0.16, linewidth=0.7)


def _plot_observed_profile(
    scan: pd.DataFrame,
    *,
    center: int,
    sequence_length: int,
    zoom_half_width: int,
    threshold_years: float,
    window_size: int,
    output_path: Path,
    pair_count: int = 2_000,
    series_label: str = "selected pseudo-empirical decode",
    statistic_label: str = "Mean inferred P",
) -> None:
    fig, axes = plt.subplots(
        1, 2, figsize=(19, 7.2), sharey=True, constrained_layout=True
    )
    limits = [(0, sequence_length), (center - zoom_half_width, center + zoom_half_width)]
    titles = ["Full 10 Mb region", f"Selected-site zoom (+/-{zoom_half_width / 1e3:g} kb)"]
    for axis, (left, right), title in zip(axes, limits, titles):
        subset = scan[scan["position_0based"].between(left, right)]
        axis.plot(
            subset["position_0based"] / 1e6,
            subset["observed_fraction_recent"],
            color="#7c3aed",
            lw=2.2,
            marker="o",
            markersize=3.2,
            label=series_label,
        )
        axis.axvline(
            center / 1e6,
            color="#e69f00",
            ls="--",
            lw=1.6,
            label="selected site",
        )
        axis.set_xlim(left / 1e6, right / 1e6)
        axis.set_title(title, fontsize=FONTS["title"])
        _format_axis(
            axis,
            ylabel=(
                f"{statistic_label} across {pair_count:,} pairs\n"
                f"(TMRCA < {threshold_years:g} years)"
                if axis is axes[0]
                else None
            ),
        )
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper left",
        bbox_to_anchor=(1.005, 0.90),
        fontsize=FONTS["legend"],
    )
    fig.suptitle(
        f"Recent-coalescence probability profile ({window_size / 1e3:g} kb grid)",
        fontsize=FONTS["suptitle"],
    )
    fig.savefig(output_path, dpi=190, bbox_inches="tight")
    plt.close(fig)


def _shade_regions(axis, regions: pd.DataFrame, *, label: bool = False) -> None:
    for index, region in enumerate(regions.itertuples(index=False)):
        axis.axvspan(
            region.start_mb,
            region.end_mb,
            color="#d55e00",
            alpha=0.14,
            label="pointwise p<0.05" if label and index == 0 else None,
        )
        axis.axvline(region.start_mb, color="#d55e00", lw=1.0, ls=":")
        axis.axvline(region.end_mb, color="#d55e00", lw=1.0, ls=":")


def _plot_null_spatial_calibration(
    scan: pd.DataFrame,
    regions: pd.DataFrame,
    *,
    center: int,
    sequence_length: int,
    zoom_half_width: int,
    threshold_years: float,
    window_size: int,
    output_path: Path,
    pair_count: int = 2_000,
    series_label: str = "selected pseudo-empirical decode",
    statistic_label: str = "Mean inferred P",
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(20, 13), constrained_layout=True)
    limits = [(0, sequence_length), (center - zoom_half_width, center + zoom_half_width)]
    titles = ["Full 10 Mb region", f"Selected-site zoom (+/-{zoom_half_width / 1e3:g} kb)"]
    for column, ((left, right), title) in enumerate(zip(limits, titles)):
        subset = scan[scan["position_0based"].between(left, right)]
        local_regions = regions[
            (regions["end_0based_exclusive"] > left)
            & (regions["start_0based"] < right)
        ]
        x = subset["position_0based"].to_numpy(dtype=float) / 1e6
        top = axes[0, column]
        top.fill_between(
            x,
            subset["neutral_ci95_lower"].to_numpy(dtype=float),
            subset["neutral_ci95_upper"].to_numpy(dtype=float),
            color="0.72",
            alpha=0.48,
            linewidth=0,
            label="neutral 95% interval",
        )
        top.plot(
            x,
            subset["neutral_median_fraction_recent"],
            color="0.35",
            lw=1.4,
            label="neutral median",
        )
        top.plot(
            x,
            subset["observed_fraction_recent"],
            color="#7c3aed",
            lw=2.0,
            label=series_label,
        )
        _shade_regions(top, local_regions, label=True)
        top.axvline(center / 1e6, color="#e69f00", ls="--", lw=1.5)
        top.set_xlim(left / 1e6, right / 1e6)
        top.set_title(title, fontsize=FONTS["title"])
        _format_axis(
            top,
            ylabel=(
                f"{statistic_label} across {pair_count:,} pairs\n"
                f"(TMRCA < {threshold_years:g} years)"
                if column == 0
                else None
            ),
        )

        bottom = axes[1, column]
        bottom.plot(
            x,
            -np.log10(subset["p_upper"].to_numpy(dtype=float)),
            color="#0072b2",
            lw=1.8,
            marker="o",
            markersize=3,
        )
        bottom.axhline(-np.log10(0.05), color="#d55e00", ls="--", lw=1.5)
        _shade_regions(bottom, local_regions)
        bottom.axvline(center / 1e6, color="#e69f00", ls="--", lw=1.5)
        bottom.set_xlim(left / 1e6, right / 1e6)
        _format_axis(bottom, ylabel="-log10(pointwise p)" if column == 0 else None)
        if column == 1:
            center_row = scan.iloc[
                np.argmin(
                    np.abs(scan["position_0based"].to_numpy(dtype=float) - center)
                )
            ]
            bottom.text(
                0.97,
                0.94,
                f"selected base: p={center_row['p_upper']:.4f}\n"
                f"BH q={center_row['q_bh']:.4f}",
                transform=bottom.transAxes,
                ha="right",
                va="top",
                fontsize=FONTS["annotation"],
            )
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper left",
        bbox_to_anchor=(1.005, 0.93),
        fontsize=FONTS["legend"],
    )
    fig.suptitle(
        f"Selected decode versus neutral decoded simulations ({window_size / 1e3:g} kb grid)",
        fontsize=FONTS["suptitle"],
    )
    fig.savefig(output_path, dpi=190, bbox_inches="tight")
    plt.close(fig)


def analyze_two_epoch_spatial_truth(
    source_dir: str | Path,
    output_dir: str | Path | None = None,
    *,
    executable: str | Path | None = None,
    workers: int = 20,
    window_size: int = 20_000,
    zoom_half_width: int = 500_000,
) -> dict:
    """Reconstruct the accepted run and matched nulls for a spatial truth scan."""
    if workers < 1 or window_size < 1 or zoom_half_width < 1:
        raise ValueError("workers, window size, and zoom width must be positive")
    source_dir = Path(source_dir).resolve()
    output_dir = source_dir if output_dir is None else Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    with (source_dir / "metrics.json").open(encoding="utf-8") as handle:
        metrics = json.load(handle)
    neutral_expected = pd.read_csv(
        source_dir / "neutral_statistics.tsv", sep="\t"
    ).set_index("replicate")
    selected_expected = pd.read_csv(
        source_dir / "selected_observation.tsv", sep="\t"
    ).iloc[0]

    ancestral_size = int(metrics["demography"]["ancestral_population_size"])
    present_size = int(metrics["demography"]["present_population_size"])
    change_age = int(metrics["demography"]["size_change_generations_ago"])
    sample_diploids = int(metrics["sample_diploids"])
    sequence_length = int(metrics["sequence_length"])
    variant_age = int(metrics["variant_age_generations"])
    generation_time = float(metrics["generation_time_years"])
    selection_coefficient = float(metrics["selection_coefficient"])
    mutation_rate = float(metrics["mutation_rate"])
    recombination_rate = float(metrics["recombination_rate"])
    seed = int(metrics["base_seed"])
    center = sequence_length // 2
    positions = window_centers(sequence_length, window_size, center)
    started = perf_counter()

    with TemporaryDirectory(prefix="gamma_smc_spatial_truth_") as temporary:
        temporary = Path(temporary)

        def neutral_task(replicate: int) -> pd.DataFrame:
            ts, _ = run_slim_recent_sweep(
                temporary / f"neutral_{replicate}.trees",
                executable=executable,
                population_size=ancestral_size,
                present_population_size=present_size,
                size_change_generations_ago=change_age,
                sample_diploids=sample_diploids,
                sequence_length=sequence_length,
                sweep_position=center,
                selection_coefficient=0.0,
                age_generations=variant_age,
                mutation_rate=mutation_rate,
                recombination_rate=recombination_rate,
                seed=seed + 1_000_000 + replicate * 10,
            )
            profile = within_individual_tmrca_grid(ts, positions, variant_age)
            center_value = float(
                profile.loc[
                    np.isclose(profile["position_0based"], center),
                    "mean_p_tmrca_lt_threshold",
                ].iloc[0]
            )
            expected = float(neutral_expected.loc[replicate, "center_fraction_recent"])
            if not np.isclose(center_value, expected, atol=1e-12, rtol=0):
                raise RuntimeError(f"neutral replicate {replicate} did not reproduce")
            profile["replicate"] = int(replicate)
            return profile

        replicate_ids = neutral_expected.index.astype(int).tolist()
        if workers == 1:
            neutral_profiles = [neutral_task(value) for value in replicate_ids]
        else:
            with ThreadPoolExecutor(max_workers=workers) as executor:
                neutral_profiles = list(executor.map(neutral_task, replicate_ids))

        attempt = int(selected_expected["attempt"])
        selected_ts, selected_run = run_slim_recent_sweep(
            temporary / f"selected_{attempt}.trees",
            executable=executable,
            population_size=ancestral_size,
            present_population_size=present_size,
            size_change_generations_ago=change_age,
            sample_diploids=sample_diploids,
            sequence_length=sequence_length,
            sweep_position=center,
            selection_coefficient=selection_coefficient,
            age_generations=variant_age,
            mutation_rate=mutation_rate,
            recombination_rate=recombination_rate,
            seed=seed + 10_000_000 + attempt * 10,
            capture_focal_genotypes=True,
        )
        if not np.isclose(
            selected_run["realized_population_allele_frequency"],
            selected_expected["population_allele_frequency"],
            atol=1e-12,
            rtol=0,
        ):
            raise RuntimeError("selected population allele frequency did not reproduce")
        observed = within_individual_tmrca_grid(
            selected_ts, positions, variant_age
        )

    center_observed = float(
        observed.loc[
            np.isclose(observed["position_0based"], center),
            "mean_p_tmrca_lt_threshold",
        ].iloc[0]
    )
    if not np.isclose(
        center_observed,
        selected_expected["center_fraction_recent"],
        atol=1e-12,
        rtol=0,
    ):
        raise RuntimeError("selected center statistic did not reproduce")

    neutral = pd.concat(neutral_profiles, ignore_index=True)
    scan = calibrate_spatial_windows(observed, neutral)
    regions = significant_regions(
        scan,
        sequence_length=sequence_length,
        window_size=window_size,
    )
    observed.to_csv(
        output_dir / "selected_recent_probability_profile.tsv", sep="\t", index=False
    )
    neutral.to_csv(
        output_dir / "neutral_recent_probability_profiles.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    scan.to_csv(
        output_dir / "spatial_pointwise_calibration.tsv", sep="\t", index=False
    )
    regions.to_csv(
        output_dir / "spatial_significant_regions.tsv", sep="\t", index=False
    )
    threshold_years = variant_age * generation_time
    _plot_observed_profile(
        scan,
        center=center,
        sequence_length=sequence_length,
        zoom_half_width=zoom_half_width,
        threshold_years=threshold_years,
        window_size=window_size,
        output_path=output_dir / "selected_recent_probability_spatial.png",
        pair_count=sample_diploids,
        series_label="selected simulation truth",
        statistic_label="Fraction with true",
    )
    _plot_null_spatial_calibration(
        scan,
        regions,
        center=center,
        sequence_length=sequence_length,
        zoom_half_width=zoom_half_width,
        threshold_years=threshold_years,
        window_size=window_size,
        output_path=output_dir / "selected_vs_neutral_spatial_pvalues.png",
        pair_count=sample_diploids,
        series_label="selected simulation truth",
        statistic_label="Fraction with true",
    )

    center_row = scan.loc[np.isclose(scan["position_0based"], center)].iloc[0]
    result = {
        "source_dir": str(source_dir),
        "truth_not_decoder_posterior": True,
        "statistic": (
            f"fraction of {sample_diploids} within-individual pairs with true "
            f"TMRCA < {threshold_years:g} years"
        ),
        "window_grid_bp": int(window_size),
        "window_interpretation": "statistic evaluated at each grid-window center; breakpoints are halfway between adjacent centers",
        "sequence_length": sequence_length,
        "zoom_half_width": int(zoom_half_width),
        "sample_diploids": sample_diploids,
        "sample_haplotypes": int(2 * sample_diploids),
        "neutral_replicates": int(neutral["replicate"].nunique()),
        "selected_attempt": int(selected_expected["attempt"]),
        "center_observed_fraction_recent": float(center_row["observed_fraction_recent"]),
        "center_p_upper": float(center_row["p_upper"]),
        "pointwise_alpha": 0.05,
        "multiple_testing_note": "highlighted regions use raw pointwise Monte Carlo p<0.05; q_bh is reported separately in the calibration TSV",
        "n_significant_windows": int(np.count_nonzero(scan["p_upper"] < 0.05)),
        "n_significant_regions": int(len(regions)),
        "significant_regions": regions.to_dict("records"),
        "deterministic_reconstruction_verified": True,
        "workers": int(workers),
        "elapsed_seconds": float(perf_counter() - started),
    }
    with (output_dir / "spatial_scan_metrics.json").open("w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2)
    return result
