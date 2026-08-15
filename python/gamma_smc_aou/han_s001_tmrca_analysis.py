"""Additive inference for natural-AF Han ``s=0.001`` TMRCA scores.

The analysis unit is one independently simulated locus.  Pair, position, and
threshold observations remain nested inside that unit.  The primary score is
the overall-pair mean ``P(TMRCA < x)`` in a prespecified window; no focal
genotype-class contrast or final-frequency acceptance condition is introduced
here.
"""

from __future__ import annotations

import gzip
import hashlib
import io
import math
import os
import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.stats import rankdata

from .eas_sweep_models import GENERATION_TIME_YEARS, TMRCA_THRESHOLDS_YEARS
from .eas_sweep_study import BASE_SEED


SCHEMA_VERSION = "gamma-smc.han-s001-tmrca-analysis/v1"
EXPECTED_SELECTION_COEFFICIENT = 0.001
PANEL_DIPLOIDS = 100
PANEL_HAPLOTYPES = 2 * PANEL_DIPLOIDS
ALLOWED_CLASSES = ("neutral", "selected")
ALLOWED_SOURCES = ("gamma_smc", "tree_truth")
ALLOWED_WINDOWS = ("local_100kb", "focal")
THRESHOLDS_YEARS = tuple(int(value) for value in TMRCA_THRESHOLDS_YEARS)
DEFAULT_BOOTSTRAP_DRAWS = 2_000
DEFAULT_PERMUTATION_DRAWS = 9_999
DEFAULT_AF_STRATUM_WIDTH = 0.05
PRIMARY_EMPIRICAL_ESTIMAND = "primary_empirical_overall_p_tmrca"
PRIMARY_EMPIRICAL_NULL_MODEL = (
    "all_natural_neutral_population_survivors_no_af_conditioning"
)
SECONDARY_EMPIRICAL_ESTIMAND = "secondary_empirical_sample_detection_and_af_matched"
SECONDARY_EMPIRICAL_NULL_MODEL = (
    "natural_neutral_population_survivors_matched_on_sample_detection_and_af"
)

SCORE_COLUMNS = (
    "unit_id",
    "simulation_class",
    "selection_coefficient",
    "seed",
    "source",
    "window",
    "threshold_years",
    "score",
    "final_population_af",
    "sample_alt_count",
    "sample_af",
    "sample_detected",
)
SCORE_KEY = ("unit_id", "source", "window", "threshold_years")
UNIT_METADATA_COLUMNS = (
    "unit_id",
    "simulation_class",
    "selection_coefficient",
    "seed",
    "final_population_af",
    "sample_alt_count",
    "sample_af",
    "sample_detected",
)

EMPIRICAL_COLUMNS = (
    "empirical_id",
    "variant_id",
    "population",
    "source",
    "window",
    "threshold_years",
    "score",
    "sample_alt_count",
    "sample_af",
    "sample_detected",
    "source_data_sha256",
    "sample_manifest_sha256",
    "decoder_binary_sha256",
    "decoder_parameters_sha256",
    "ne_history_sha256",
    "pairing_scheme",
    "qc_status",
    "data_status",
)
EMPIRICAL_RESULT_COLUMNS = (
    "status",
    "estimand",
    "null_model",
    "empirical_id",
    "variant_id",
    "source",
    "window",
    "threshold_years",
    "threshold_generations",
    "observed_score",
    "observed_sample_af",
    "af_match_tolerance",
    "sample_detection_matched",
    "n_matched_neutral",
    "neutral_upper_exceedance_count",
    "mc_p_upper",
    "minimum_attainable_p",
    "null_definition",
)

EMPIRICAL_OMNIBUS_COLUMNS = (
    "status",
    "estimand",
    "null_model",
    "empirical_id",
    "variant_id",
    "source",
    "window",
    "n_thresholds",
    "thresholds_years",
    "observed_min_pointwise_p",
    "threshold_years_at_min_pointwise_p",
    "n_shared_neutral_profiles",
    "neutral_profile_exceedance_count",
    "omnibus_mc_p_upper",
    "minimum_attainable_p",
    "af_match_tolerance",
    "sample_detection_matched",
    "omnibus_method",
    "null_definition",
)

TABLE_FILENAMES = {
    "auc_by_threshold": "han_s001_auc_by_threshold.tsv",
    "auc_omnibus": "han_s001_auc_omnibus.tsv",
    "selected_unit_pvalues": "han_s001_selected_unit_pvalues.tsv.gz",
    "af_only_benchmarks": "han_s001_af_only_benchmarks.tsv",
    "sample_af_stratified_auc": "han_s001_sample_af_stratified_auc.tsv",
    "sample_af_stratified_omnibus": ("han_s001_sample_af_stratified_omnibus.tsv"),
    "empirical_pvalues": "han_s001_empirical_pvalues.tsv",
    "empirical_omnibus": "han_s001_empirical_omnibus.tsv",
    "empirical_input": "han_s001_empirical_input.tsv",
}


def _stable_seed(base_seed: int, label: str) -> int:
    validated_seed = _nonnegative_integer(base_seed, label="base_seed")
    payload = f"{validated_seed}:{label}".encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little")


def _positive_integer(value: Any, *, label: str) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise ValueError(f"{label} must be a positive integer")
    result = int(value)
    if result < 1:
        raise ValueError(f"{label} must be a positive integer")
    return result


def _nonnegative_integer(value: Any, *, label: str) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise ValueError(f"{label} must be a nonnegative integer")
    result = int(value)
    if result < 0:
        raise ValueError(f"{label} must be a nonnegative integer")
    return result


def _strict_boolean(series: pd.Series, *, label: str) -> np.ndarray:
    values: list[bool] = []
    for value in series.tolist():
        if isinstance(value, (bool, np.bool_)):
            values.append(bool(value))
            continue
        if isinstance(value, (int, np.integer)) and int(value) in {0, 1}:
            values.append(bool(value))
            continue
        text = str(value).strip().casefold()
        if text in {"true", "1"}:
            values.append(True)
        elif text in {"false", "0"}:
            values.append(False)
        else:
            raise ValueError(f"{label} contains a non-boolean value: {value!r}")
    return np.asarray(values, dtype=bool)


def _numeric_column(frame: pd.DataFrame, column: str) -> np.ndarray:
    values = pd.to_numeric(frame[column], errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError(f"{column} must contain only finite numeric values")
    return values


def validate_replicate_scores(scores: pd.DataFrame) -> pd.DataFrame:
    """Validate and normalize the exact natural-AF replicate-score contract."""

    if not isinstance(scores, pd.DataFrame):
        raise TypeError("scores must be a pandas DataFrame")
    missing = sorted(set(SCORE_COLUMNS).difference(scores.columns))
    unexpected = sorted(set(scores.columns).difference(SCORE_COLUMNS))
    if missing or unexpected:
        pieces = []
        if missing:
            pieces.append("missing: " + ", ".join(missing))
        if unexpected:
            pieces.append("unexpected: " + ", ".join(unexpected))
        raise ValueError(
            "replicate-score columns are incompatible (" + "; ".join(pieces) + ")"
        )
    if scores.empty:
        raise ValueError("replicate-score table is empty")

    result = scores.loc[:, SCORE_COLUMNS].copy()
    for column in ("unit_id", "simulation_class", "source", "window"):
        if result[column].isna().any():
            raise ValueError(f"{column} contains a missing value")
        result[column] = result[column].astype(str).str.strip()
        if result[column].eq("").any():
            raise ValueError(f"{column} contains an empty value")
    if set(result["simulation_class"]) != set(ALLOWED_CLASSES):
        raise ValueError("simulation_class must contain exactly neutral and selected")
    if set(result["source"]) != set(ALLOWED_SOURCES):
        raise ValueError("source must contain exactly gamma_smc and tree_truth")
    if set(result["window"]) != set(ALLOWED_WINDOWS):
        raise ValueError("window must contain exactly local_100kb and focal")

    coefficient = _numeric_column(result, "selection_coefficient")
    neutral = result["simulation_class"].eq("neutral").to_numpy()
    selected = ~neutral
    if not np.allclose(coefficient[neutral], 0.0, rtol=0.0, atol=1e-15):
        raise ValueError("neutral rows must have selection_coefficient 0")
    if not np.allclose(
        coefficient[selected],
        EXPECTED_SELECTION_COEFFICIENT,
        rtol=0.0,
        atol=1e-15,
    ):
        raise ValueError("selected rows must have selection_coefficient 0.001")
    result["selection_coefficient"] = coefficient

    seed = _numeric_column(result, "seed")
    if np.any(seed < 0) or not np.array_equal(seed, np.floor(seed)):
        raise ValueError("seed values must be nonnegative integers")
    result["seed"] = seed.astype(np.int64)

    thresholds = _numeric_column(result, "threshold_years")
    if not np.array_equal(thresholds, np.floor(thresholds)):
        raise ValueError("threshold_years values must be integers")
    result["threshold_years"] = thresholds.astype(np.int64)
    if set(result["threshold_years"].astype(int)) != set(THRESHOLDS_YEARS):
        raise ValueError("replicate-score table lacks the exact seven thresholds")

    score = _numeric_column(result, "score")
    if np.any((score < 0.0) | (score > 1.0)):
        raise ValueError("score must lie in [0, 1]")
    result["score"] = score

    final_af = _numeric_column(result, "final_population_af")
    if np.any((final_af <= 0.0) | (final_af > 1.0)):
        raise ValueError(
            "final_population_af must lie in (0, 1] under population survival"
        )
    result["final_population_af"] = final_af

    sample_count = _numeric_column(result, "sample_alt_count")
    if np.any(
        (sample_count < 0) | (sample_count > PANEL_HAPLOTYPES)
    ) or not np.array_equal(sample_count, np.floor(sample_count)):
        raise ValueError("sample_alt_count must be an integer in [0, 200]")
    result["sample_alt_count"] = sample_count.astype(np.int64)
    sample_af = _numeric_column(result, "sample_af")
    if not np.allclose(
        sample_af,
        sample_count / PANEL_HAPLOTYPES,
        rtol=0.0,
        atol=1e-12,
    ):
        raise ValueError("sample_af disagrees with sample_alt_count / 200")
    result["sample_af"] = sample_af
    detected = _strict_boolean(result["sample_detected"], label="sample_detected")
    if not np.array_equal(detected, sample_count > 0):
        raise ValueError("sample_detected disagrees with sample_alt_count > 0")
    result["sample_detected"] = detected

    if result.duplicated(list(SCORE_KEY), keep=False).any():
        raise ValueError("replicate-score table contains a duplicate key")

    expected_combinations = {
        (source, window) for source in ALLOWED_SOURCES for window in ALLOWED_WINDOWS
    }
    unit_metadata = result.loc[:, UNIT_METADATA_COLUMNS].drop_duplicates()
    if unit_metadata["unit_id"].duplicated(keep=False).any():
        raise ValueError("unit metadata changes across score rows")
    if unit_metadata["seed"].duplicated(keep=False).any():
        raise ValueError("independent simulation units must have unique seeds")
    if unit_metadata["unit_id"].nunique() < 2:
        raise ValueError("at least two independent units are required")

    for unit_id, unit in result.groupby("unit_id", sort=False):
        combinations = set(zip(unit["source"], unit["window"], strict=True))
        if combinations != expected_combinations:
            raise ValueError(f"unit {unit_id} lacks a source/window combination")
        for labels, curve in unit.groupby(["source", "window"], sort=False):
            ordered = curve.sort_values("threshold_years")
            if tuple(ordered["threshold_years"].astype(int)) != THRESHOLDS_YEARS:
                raise ValueError(
                    f"unit {unit_id} {labels} lacks the exact threshold curve"
                )
            if np.any(np.diff(ordered["score"].to_numpy(dtype=float)) < -1e-12):
                raise ValueError(
                    f"unit {unit_id} {labels} has a nonmonotone probability curve"
                )

    for labels, group in result.groupby(["source", "window"], sort=False):
        if set(group["simulation_class"]) != set(ALLOWED_CLASSES):
            raise ValueError(f"source/window cell {labels} lacks a simulation class")

    return result.sort_values(list(SCORE_KEY), kind="mergesort").reset_index(drop=True)


def load_replicate_scores(path: str | Path) -> pd.DataFrame:
    """Read and validate a plain or gzip-compressed TSV score table."""

    resolved = Path(path).resolve()
    if not resolved.is_file():
        raise ValueError(f"replicate-score table is absent: {resolved}")
    return validate_replicate_scores(pd.read_csv(resolved, sep="\t"))


def _roc_auc(selected: Sequence[float], neutral: Sequence[float]) -> float:
    selected_values = np.asarray(selected, dtype=float)
    neutral_values = np.asarray(neutral, dtype=float)
    if not len(selected_values) or not len(neutral_values):
        return np.nan
    greater = (selected_values[:, None] > neutral_values[None, :]).sum()
    tied = (selected_values[:, None] == neutral_values[None, :]).sum()
    return float((greater + 0.5 * tied) / (len(selected_values) * len(neutral_values)))


def _profile_matrices(
    frame: pd.DataFrame, *, source: str, window: str
) -> tuple[pd.DataFrame, pd.DataFrame]:
    subset = frame[(frame["source"] == source) & (frame["window"] == window)]
    pivot = subset.pivot(
        index="unit_id", columns="threshold_years", values="score"
    ).reindex(columns=THRESHOLDS_YEARS)
    metadata = frame.loc[:, ["unit_id", "simulation_class"]].drop_duplicates()
    pivot = pivot.join(metadata.set_index("unit_id"), validate="one_to_one")
    selected = pivot[pivot["simulation_class"] == "selected"].drop(
        columns="simulation_class"
    )
    neutral = pivot[pivot["simulation_class"] == "neutral"].drop(
        columns="simulation_class"
    )
    if (
        selected.empty
        or neutral.empty
        or selected.isna().any().any()
        or neutral.isna().any().any()
    ):
        raise ValueError(f"incomplete selected/neutral profile for {source}/{window}")
    return selected.sort_index(), neutral.sort_index()


def _bootstrap_auc_profile(
    selected: np.ndarray,
    neutral: np.ndarray,
    *,
    draws: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    draw_count = _positive_integer(draws, label="bootstrap_draws")
    rng = np.random.default_rng(seed)
    n_selected, n_thresholds = selected.shape
    n_neutral = neutral.shape[0]
    samples = np.empty((draw_count, n_thresholds), dtype=float)
    completed = 0
    chunk_size = 32
    while completed < draw_count:
        size = min(chunk_size, draw_count - completed)
        selected_draws = selected[rng.integers(0, n_selected, size=(size, n_selected))]
        neutral_draws = neutral[rng.integers(0, n_neutral, size=(size, n_neutral))]
        greater = (selected_draws[:, :, None, :] > neutral_draws[:, None, :, :]).sum(
            axis=(1, 2)
        )
        tied = (selected_draws[:, :, None, :] == neutral_draws[:, None, :, :]).sum(
            axis=(1, 2)
        )
        samples[completed : completed + size, :] = (greater + 0.5 * tied) / (
            n_selected * n_neutral
        )
        completed += size
    return (
        np.quantile(samples, 0.025, axis=0),
        np.quantile(samples, 0.975, axis=0),
    )


def _permutation_auc_profile(
    selected: np.ndarray,
    neutral: np.ndarray,
    *,
    draws: int,
    seed: int,
) -> dict[str, Any]:
    draw_count = _positive_integer(draws, label="permutation_draws")
    combined = np.vstack((selected, neutral))
    n_selected = len(selected)
    n_neutral = len(neutral)
    observed = np.asarray(
        [
            _roc_auc(selected[:, index], neutral[:, index])
            for index in range(selected.shape[1])
        ],
        dtype=float,
    )
    ranks = np.column_stack(
        [
            rankdata(combined[:, index], method="average")
            for index in range(combined.shape[1])
        ]
    )
    point_exceedances = np.zeros(combined.shape[1], dtype=np.int64)
    global_exceedances = 0
    observed_max = float(np.max(observed))
    rng = np.random.default_rng(seed)
    completed = 0
    chunk_size = 256
    while completed < draw_count:
        size = min(chunk_size, draw_count - completed)
        random_keys = rng.random((size, len(combined)))
        selected_indices = np.argpartition(random_keys, n_selected - 1, axis=1)[
            :, :n_selected
        ]
        rank_sums = ranks[selected_indices].sum(axis=1)
        permuted = (rank_sums - n_selected * (n_selected + 1) / 2) / (
            n_selected * n_neutral
        )
        point_exceedances += np.count_nonzero(
            permuted >= observed[None, :] - 1e-14, axis=0
        )
        global_exceedances += int(
            np.count_nonzero(np.max(permuted, axis=1) >= observed_max - 1e-14)
        )
        completed += size
    denominator = draw_count + 1
    return {
        "observed": observed,
        "point_exceedances": point_exceedances,
        "point_p": (1 + point_exceedances) / denominator,
        "observed_max": observed_max,
        "observed_max_index": int(np.argmax(observed)),
        "global_exceedances": int(global_exceedances),
        "global_p": float((1 + global_exceedances) / denominator),
    }


def _profile_inference(
    selected: np.ndarray,
    neutral: np.ndarray,
    *,
    bootstrap_draws: int,
    permutation_draws: int,
    bootstrap_seed: int,
    permutation_seed: int,
) -> dict[str, Any]:
    low, high = _bootstrap_auc_profile(
        selected, neutral, draws=bootstrap_draws, seed=bootstrap_seed
    )
    permutation = _permutation_auc_profile(
        selected, neutral, draws=permutation_draws, seed=permutation_seed
    )
    return {"ci_low": low, "ci_high": high, **permutation}


def auc_by_threshold(
    scores: pd.DataFrame,
    *,
    bootstrap_draws: int = DEFAULT_BOOTSTRAP_DRAWS,
    permutation_draws: int = DEFAULT_PERMUTATION_DRAWS,
    base_seed: int = BASE_SEED,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute primary, unflipped selected-versus-neutral ROC AUC profiles."""

    base_seed = _nonnegative_integer(base_seed, label="base_seed")
    frame = validate_replicate_scores(scores)
    rows: list[dict[str, Any]] = []
    omnibus_rows: list[dict[str, Any]] = []
    for source in ALLOWED_SOURCES:
        for window in ALLOWED_WINDOWS:
            selected, neutral = _profile_matrices(frame, source=source, window=window)
            selected_values = selected.to_numpy(dtype=float)
            neutral_values = neutral.to_numpy(dtype=float)
            seed_label = f"primary:{source}:{window}"
            bootstrap_seed = _stable_seed(base_seed, "bootstrap:" + seed_label)
            permutation_seed = _stable_seed(base_seed, "permutation:" + seed_label)
            inference = _profile_inference(
                selected_values,
                neutral_values,
                bootstrap_draws=bootstrap_draws,
                permutation_draws=permutation_draws,
                bootstrap_seed=bootstrap_seed,
                permutation_seed=permutation_seed,
            )
            for index, threshold in enumerate(THRESHOLDS_YEARS):
                selected_column = selected_values[:, index]
                neutral_column = neutral_values[:, index]
                rows.append(
                    {
                        "schema": SCHEMA_VERSION,
                        "estimand": ("primary_natural_survivor_overall_p_tmrca"),
                        "cohort": "population_survival_only_no_final_af_conditioning",
                        "source": source,
                        "window": window,
                        "threshold_years": threshold,
                        "threshold_generations": threshold / GENERATION_TIME_YEARS,
                        "n_selected": len(selected_column),
                        "n_neutral": len(neutral_column),
                        "mean_selected_score": float(np.mean(selected_column)),
                        "median_selected_score": float(np.median(selected_column)),
                        "mean_neutral_score": float(np.mean(neutral_column)),
                        "median_neutral_score": float(np.median(neutral_column)),
                        "mean_selected_minus_neutral": float(
                            np.mean(selected_column) - np.mean(neutral_column)
                        ),
                        "roc_auc": float(inference["observed"][index]),
                        "roc_auc_ci95_low": float(inference["ci_low"][index]),
                        "roc_auc_ci95_high": float(inference["ci_high"][index]),
                        "auc_permutation_exceedances": int(
                            inference["point_exceedances"][index]
                        ),
                        "auc_permutation_p_upper": float(inference["point_p"][index]),
                        "bootstrap_draws": int(bootstrap_draws),
                        "bootstrap_seed": int(bootstrap_seed),
                        "permutation_draws": int(permutation_draws),
                        "permutation_seed": int(permutation_seed),
                        "permutation_family": (
                            "shared_whole_unit_labels_across_seven_thresholds"
                        ),
                        "score_direction": "larger_is_more_selected_no_flipping",
                    }
                )
            best_index = int(inference["observed_max_index"])
            omnibus_rows.append(
                {
                    "schema": SCHEMA_VERSION,
                    "estimand": "primary_natural_survivor_overall_p_tmrca",
                    "cohort": "population_survival_only_no_final_af_conditioning",
                    "source": source,
                    "window": window,
                    "n_thresholds": len(THRESHOLDS_YEARS),
                    "thresholds_years": ",".join(map(str, THRESHOLDS_YEARS)),
                    "max_observed_auc": float(inference["observed_max"]),
                    "threshold_years_at_max_auc": THRESHOLDS_YEARS[best_index],
                    "max_auc_minus_half": float(inference["observed_max"] - 0.5),
                    "max_auc_permutation_exceedances": int(
                        inference["global_exceedances"]
                    ),
                    "max_auc_permutation_p_upper": float(inference["global_p"]),
                    "permutation_draws": int(permutation_draws),
                    "permutation_seed": int(permutation_seed),
                    "omnibus_method": (
                        "max_auc_shared_whole_unit_label_permutation_plus_one"
                    ),
                    "score_direction": "larger_is_more_selected_no_flipping",
                }
            )
    return pd.DataFrame(rows), pd.DataFrame(omnibus_rows)


def selected_unit_neutral_tail_pvalues(scores: pd.DataFrame) -> pd.DataFrame:
    """Score each selected unit against the matched independent neutral bank."""

    frame = validate_replicate_scores(scores)
    rows: list[dict[str, Any]] = []
    for (source, window, threshold), group in frame.groupby(
        ["source", "window", "threshold_years"], sort=True
    ):
        null = group.loc[group["simulation_class"] == "neutral", "score"].to_numpy(
            dtype=float
        )
        selected = group[group["simulation_class"] == "selected"]
        for record in selected.to_dict(orient="records"):
            observed = float(record["score"])
            exceedances = int(np.count_nonzero(null >= observed - 1e-14))
            rows.append(
                {
                    "schema": SCHEMA_VERSION,
                    "unit_id": str(record["unit_id"]),
                    "seed": int(record["seed"]),
                    "source": source,
                    "window": window,
                    "threshold_years": int(threshold),
                    "threshold_generations": (float(threshold) / GENERATION_TIME_YEARS),
                    "observed_score": observed,
                    "n_neutral": len(null),
                    "neutral_mean": float(np.mean(null)),
                    "neutral_median": float(np.median(null)),
                    "neutral_upper_exceedance_count": exceedances,
                    "mc_p_upper": float((1 + exceedances) / (len(null) + 1)),
                    "minimum_attainable_p": float(1 / (len(null) + 1)),
                    "pvalue_method": "neutral_upper_tail_plus_one",
                    "neutral_bank": "same_source_window_threshold_natural_survivors",
                }
            )
    return (
        pd.DataFrame(rows)
        .sort_values(
            ["source", "window", "unit_id", "threshold_years"], kind="mergesort"
        )
        .reset_index(drop=True)
    )


def _unit_table(frame: pd.DataFrame) -> pd.DataFrame:
    return (
        frame.loc[:, UNIT_METADATA_COLUMNS]
        .drop_duplicates()
        .sort_values("unit_id", kind="mergesort")
        .reset_index(drop=True)
    )


def af_only_benchmarks(
    scores: pd.DataFrame,
    *,
    bootstrap_draws: int = DEFAULT_BOOTSTRAP_DRAWS,
    permutation_draws: int = DEFAULT_PERMUTATION_DRAWS,
    base_seed: int = BASE_SEED,
) -> pd.DataFrame:
    """Return AF-only discrimination benchmarks on the same simulation units."""

    base_seed = _nonnegative_integer(base_seed, label="base_seed")
    frame = validate_replicate_scores(scores)
    units = _unit_table(frame)
    metrics = ("final_population_af", "sample_af")
    selected = units[units["simulation_class"] == "selected"]
    neutral = units[units["simulation_class"] == "neutral"]
    selected_values = selected.loc[:, metrics].to_numpy(dtype=float)
    neutral_values = neutral.loc[:, metrics].to_numpy(dtype=float)
    bootstrap_seed = _stable_seed(base_seed, "bootstrap:af-only")
    permutation_seed = _stable_seed(base_seed, "permutation:af-only")
    inference = _profile_inference(
        selected_values,
        neutral_values,
        bootstrap_draws=bootstrap_draws,
        permutation_draws=permutation_draws,
        bootstrap_seed=bootstrap_seed,
        permutation_seed=permutation_seed,
    )
    rows = []
    for index, metric in enumerate(metrics):
        rows.append(
            {
                "schema": SCHEMA_VERSION,
                "estimand": "af_only_natural_survivor_benchmark",
                "metric": metric,
                "n_selected": len(selected),
                "n_neutral": len(neutral),
                "mean_selected": float(np.mean(selected_values[:, index])),
                "median_selected": float(np.median(selected_values[:, index])),
                "mean_neutral": float(np.mean(neutral_values[:, index])),
                "median_neutral": float(np.median(neutral_values[:, index])),
                "roc_auc": float(inference["observed"][index]),
                "roc_auc_ci95_low": float(inference["ci_low"][index]),
                "roc_auc_ci95_high": float(inference["ci_high"][index]),
                "auc_permutation_exceedances": int(
                    inference["point_exceedances"][index]
                ),
                "auc_permutation_p_upper": float(inference["point_p"][index]),
                "bootstrap_draws": int(bootstrap_draws),
                "bootstrap_seed": int(bootstrap_seed),
                "permutation_draws": int(permutation_draws),
                "permutation_seed": int(permutation_seed),
                "interpretation": (
                    "benchmark_only; primary TMRCA AUC includes AF-mediated signal"
                ),
            }
        )
    return pd.DataFrame(rows)


def _af_strata(
    values: np.ndarray, width: float, sample_detected: np.ndarray
) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    detected = np.asarray(sample_detected)
    if values.shape != detected.shape:
        raise ValueError("sample-AF values and detection labels must have equal shape")
    if detected.dtype != np.bool_:
        raise ValueError("sample-detection labels must be boolean")
    width = float(width)
    if not np.isfinite(width) or width <= 0 or width > 1:
        raise ValueError("af_stratum_width must evenly divide [0, 1]")
    reciprocal = 1.0 / width
    n_strata = int(round(reciprocal))
    if not math.isclose(reciprocal, n_strata, rel_tol=0.0, abs_tol=1e-10):
        raise ValueError("af_stratum_width must evenly divide [0, 1]")
    strata = np.minimum(np.floor(values / width).astype(int), n_strata - 1)
    strata[~detected] = -1
    return strata


def _format_af_strata(strata: Sequence[int]) -> str:
    return ",".join(
        "sample_undetected_af_zero"
        if int(stratum) == -1
        else f"sample_detected_af_bin_{int(stratum)}"
        for stratum in strata
    )


def _stratified_auc_profile(
    selected: np.ndarray,
    neutral: np.ndarray,
    selected_strata: np.ndarray,
    neutral_strata: np.ndarray,
) -> tuple[np.ndarray, list[int], int]:
    common = sorted(set(selected_strata).intersection(neutral_strata))
    if not common:
        raise ValueError("sample-AF stratification has no common-support stratum")
    numerator = np.zeros(selected.shape[1], dtype=float)
    denominator = 0
    for stratum in common:
        selected_values = selected[selected_strata == stratum]
        neutral_values = neutral[neutral_strata == stratum]
        pairs = len(selected_values) * len(neutral_values)
        denominator += pairs
        for index in range(selected.shape[1]):
            numerator[index] += pairs * _roc_auc(
                selected_values[:, index], neutral_values[:, index]
            )
    return numerator / denominator, common, denominator


def _bootstrap_stratified_auc_profile(
    selected: np.ndarray,
    neutral: np.ndarray,
    selected_strata: np.ndarray,
    neutral_strata: np.ndarray,
    common: Sequence[int],
    *,
    draws: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    draw_count = _positive_integer(draws, label="bootstrap_draws")
    rng = np.random.default_rng(seed)
    output = np.zeros((draw_count, selected.shape[1]), dtype=float)
    denominator = sum(
        int(np.count_nonzero(selected_strata == stratum))
        * int(np.count_nonzero(neutral_strata == stratum))
        for stratum in common
    )
    for stratum in common:
        selected_values = selected[selected_strata == stratum]
        neutral_values = neutral[neutral_strata == stratum]
        n_selected = len(selected_values)
        n_neutral = len(neutral_values)
        selected_indices = rng.integers(0, n_selected, size=(draw_count, n_selected))
        neutral_indices = rng.integers(0, n_neutral, size=(draw_count, n_neutral))
        completed = 0
        chunk_size = 32
        while completed < draw_count:
            end = min(completed + chunk_size, draw_count)
            selected_draws = selected_values[selected_indices[completed:end]]
            neutral_draws = neutral_values[neutral_indices[completed:end]]
            greater = (
                selected_draws[:, :, None, :] > neutral_draws[:, None, :, :]
            ).sum(axis=(1, 2))
            tied = (selected_draws[:, :, None, :] == neutral_draws[:, None, :, :]).sum(
                axis=(1, 2)
            )
            output[completed:end] += greater + 0.5 * tied
            completed = end
    output /= denominator
    return np.quantile(output, 0.025, axis=0), np.quantile(output, 0.975, axis=0)


def _permutation_stratified_auc_profile(
    selected: np.ndarray,
    neutral: np.ndarray,
    selected_strata: np.ndarray,
    neutral_strata: np.ndarray,
    common: Sequence[int],
    observed: np.ndarray,
    *,
    draws: int,
    seed: int,
) -> dict[str, Any]:
    draw_count = _positive_integer(draws, label="permutation_draws")
    stratum_payloads = []
    denominator = 0
    for stratum in common:
        selected_values = selected[selected_strata == stratum]
        neutral_values = neutral[neutral_strata == stratum]
        combined = np.vstack((selected_values, neutral_values))
        n_selected = len(selected_values)
        n_neutral = len(neutral_values)
        ranks = np.column_stack(
            [
                rankdata(combined[:, index], method="average")
                for index in range(combined.shape[1])
            ]
        )
        denominator += n_selected * n_neutral
        stratum_payloads.append((ranks, n_selected, n_neutral))
    rng = np.random.default_rng(seed)
    point_exceedances = np.zeros(len(observed), dtype=np.int64)
    global_exceedances = 0
    observed_max = float(np.max(observed))
    completed = 0
    chunk_size = 128
    while completed < draw_count:
        size = min(chunk_size, draw_count - completed)
        numerator = np.zeros((size, len(observed)), dtype=float)
        for ranks, n_selected, n_neutral in stratum_payloads:
            random_keys = rng.random((size, len(ranks)))
            indices = np.argpartition(random_keys, n_selected - 1, axis=1)[
                :, :n_selected
            ]
            rank_sums = ranks[indices].sum(axis=1)
            numerator += rank_sums - n_selected * (n_selected + 1) / 2
        permuted = numerator / denominator
        point_exceedances += np.count_nonzero(
            permuted >= observed[None, :] - 1e-14, axis=0
        )
        global_exceedances += int(
            np.count_nonzero(np.max(permuted, axis=1) >= observed_max - 1e-14)
        )
        completed += size
    denominator_p = draw_count + 1
    return {
        "point_exceedances": point_exceedances,
        "point_p": (1 + point_exceedances) / denominator_p,
        "global_exceedances": global_exceedances,
        "global_p": float((1 + global_exceedances) / denominator_p),
        "observed_max": observed_max,
        "observed_max_index": int(np.argmax(observed)),
    }


def sample_af_stratified_auc(
    scores: pd.DataFrame,
    *,
    af_stratum_width: float = DEFAULT_AF_STRATUM_WIDTH,
    bootstrap_draws: int = DEFAULT_BOOTSTRAP_DRAWS,
    permutation_draws: int = DEFAULT_PERMUTATION_DRAWS,
    base_seed: int = BASE_SEED,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return a clearly secondary, common-support sample-AF-stratified AUC."""

    base_seed = _nonnegative_integer(base_seed, label="base_seed")
    frame = validate_replicate_scores(scores)
    units = _unit_table(frame).set_index("unit_id")
    rows: list[dict[str, Any]] = []
    omnibus_rows: list[dict[str, Any]] = []
    for source in ALLOWED_SOURCES:
        for window in ALLOWED_WINDOWS:
            selected_frame, neutral_frame = _profile_matrices(
                frame, source=source, window=window
            )
            selected = selected_frame.to_numpy(dtype=float)
            neutral = neutral_frame.to_numpy(dtype=float)
            selected_af = units.loc[selected_frame.index, "sample_af"].to_numpy(float)
            neutral_af = units.loc[neutral_frame.index, "sample_af"].to_numpy(float)
            selected_detected = units.loc[
                selected_frame.index, "sample_detected"
            ].to_numpy(bool)
            neutral_detected = units.loc[
                neutral_frame.index, "sample_detected"
            ].to_numpy(bool)
            selected_strata = _af_strata(
                selected_af, af_stratum_width, selected_detected
            )
            neutral_strata = _af_strata(neutral_af, af_stratum_width, neutral_detected)
            observed, common, comparable_pairs = _stratified_auc_profile(
                selected, neutral, selected_strata, neutral_strata
            )
            seed_label = f"af-stratified:{source}:{window}:{af_stratum_width:g}"
            bootstrap_seed = _stable_seed(base_seed, "bootstrap:" + seed_label)
            permutation_seed = _stable_seed(base_seed, "permutation:" + seed_label)
            low, high = _bootstrap_stratified_auc_profile(
                selected,
                neutral,
                selected_strata,
                neutral_strata,
                common,
                draws=bootstrap_draws,
                seed=bootstrap_seed,
            )
            permutation = _permutation_stratified_auc_profile(
                selected,
                neutral,
                selected_strata,
                neutral_strata,
                common,
                observed,
                draws=permutation_draws,
                seed=permutation_seed,
            )
            selected_common = int(np.count_nonzero(np.isin(selected_strata, common)))
            neutral_common = int(np.count_nonzero(np.isin(neutral_strata, common)))
            for index, threshold in enumerate(THRESHOLDS_YEARS):
                rows.append(
                    {
                        "schema": SCHEMA_VERSION,
                        "estimand": (
                            "secondary_sample_af_stratified_common_support_pair_weighted"
                        ),
                        "source": source,
                        "window": window,
                        "threshold_years": threshold,
                        "threshold_generations": threshold / GENERATION_TIME_YEARS,
                        "af_stratum_width": float(af_stratum_width),
                        "common_strata": _format_af_strata(common),
                        "af_stratum_definition": (
                            "sample-undetected AF=0 is separate; detected samples use "
                            "floor(sample_af/af_stratum_width)"
                        ),
                        "n_selected_total": len(selected),
                        "n_neutral_total": len(neutral),
                        "n_selected_common_support": selected_common,
                        "n_neutral_common_support": neutral_common,
                        "n_comparable_selected_neutral_pairs": comparable_pairs,
                        "roc_auc": float(observed[index]),
                        "roc_auc_ci95_low": float(low[index]),
                        "roc_auc_ci95_high": float(high[index]),
                        "auc_permutation_exceedances": int(
                            permutation["point_exceedances"][index]
                        ),
                        "auc_permutation_p_upper": float(permutation["point_p"][index]),
                        "bootstrap_draws": int(bootstrap_draws),
                        "bootstrap_seed": int(bootstrap_seed),
                        "permutation_draws": int(permutation_draws),
                        "permutation_seed": int(permutation_seed),
                        "permutation_method": (
                            "within_sample_af_stratum_whole_unit_labels_plus_one"
                        ),
                        "interpretation": (
                            "secondary incremental genealogy analysis; not the total "
                            "natural-selection effect"
                        ),
                        "score_direction": "larger_is_more_selected_no_flipping",
                    }
                )
            best_index = int(permutation["observed_max_index"])
            omnibus_rows.append(
                {
                    "schema": SCHEMA_VERSION,
                    "estimand": (
                        "secondary_sample_af_stratified_common_support_pair_weighted"
                    ),
                    "source": source,
                    "window": window,
                    "af_stratum_width": float(af_stratum_width),
                    "common_strata": _format_af_strata(common),
                    "af_stratum_definition": (
                        "sample-undetected AF=0 is separate; detected samples use "
                        "floor(sample_af/af_stratum_width)"
                    ),
                    "n_thresholds": len(THRESHOLDS_YEARS),
                    "max_observed_auc": float(permutation["observed_max"]),
                    "threshold_years_at_max_auc": THRESHOLDS_YEARS[best_index],
                    "max_auc_permutation_exceedances": int(
                        permutation["global_exceedances"]
                    ),
                    "max_auc_permutation_p_upper": float(permutation["global_p"]),
                    "permutation_draws": int(permutation_draws),
                    "permutation_seed": int(permutation_seed),
                    "omnibus_method": (
                        "max_auc_within_af_stratum_unit_label_permutation_plus_one"
                    ),
                }
            )
    return pd.DataFrame(rows), pd.DataFrame(omnibus_rows)


def empirical_input_template() -> pd.DataFrame:
    """Return an intentionally incomplete template that cannot yield a result."""

    rows = []
    for window in ALLOWED_WINDOWS:
        for threshold in THRESHOLDS_YEARS:
            rows.append(
                {
                    "empirical_id": "replace_with_empirical_id",
                    "variant_id": "replace_with_variant_id",
                    "population": "Han",
                    "source": "gamma_smc",
                    "window": window,
                    "threshold_years": threshold,
                    "score": np.nan,
                    "sample_alt_count": np.nan,
                    "sample_af": np.nan,
                    "sample_detected": np.nan,
                    "source_data_sha256": "replace_with_sha256",
                    "sample_manifest_sha256": "replace_with_sha256",
                    "decoder_binary_sha256": "replace_with_sha256",
                    "decoder_parameters_sha256": "replace_with_sha256",
                    "ne_history_sha256": "replace_with_sha256",
                    "pairing_scheme": "within_individual_diploid_overall",
                    "qc_status": "replace_with_pass_after_qc",
                    "data_status": "template_unfilled",
                }
            )
    return pd.DataFrame(rows, columns=EMPIRICAL_COLUMNS)


def validate_empirical_input(empirical: pd.DataFrame) -> pd.DataFrame:
    """Fail closed unless a complete, provenance-bound empirical table is given."""

    if not isinstance(empirical, pd.DataFrame):
        raise TypeError("empirical must be a pandas DataFrame")
    missing = sorted(set(EMPIRICAL_COLUMNS).difference(empirical.columns))
    unexpected = sorted(set(empirical.columns).difference(EMPIRICAL_COLUMNS))
    if missing or unexpected:
        raise ValueError("empirical columns are incompatible")
    if empirical.empty:
        raise ValueError("empirical table is empty")
    result = empirical.loc[:, EMPIRICAL_COLUMNS].copy()
    text_columns = (
        "empirical_id",
        "variant_id",
        "population",
        "source",
        "window",
        "source_data_sha256",
        "sample_manifest_sha256",
        "decoder_binary_sha256",
        "decoder_parameters_sha256",
        "ne_history_sha256",
        "pairing_scheme",
        "qc_status",
        "data_status",
    )
    for column in text_columns:
        if result[column].isna().any():
            raise ValueError(f"empirical provenance field is missing: {column}")
        result[column] = result[column].astype(str).str.strip()
        if (
            result[column].eq("").any()
            or result[column].str.startswith("replace_with_").any()
        ):
            raise ValueError(f"empirical provenance field is unset: {column}")
    if set(result["population"]) != {"Han"}:
        raise ValueError("empirical population must be Han")
    if set(result["source"]) != {"gamma_smc"}:
        raise ValueError(
            "empirical source must be gamma_smc; tree truth is unavailable"
        )
    if set(result["window"]) != set(ALLOWED_WINDOWS):
        raise ValueError("empirical input must contain both prespecified windows")
    if set(result["pairing_scheme"]) != {"within_individual_diploid_overall"}:
        raise ValueError("empirical pairing_scheme is incompatible")
    if set(result["qc_status"].str.casefold()) != {"pass"}:
        raise ValueError("empirical qc_status must be pass")
    if set(result["data_status"]) != {"empirical_observed"}:
        raise ValueError("empirical data_status must be empirical_observed")
    sha_pattern = re.compile(r"^[0-9a-f]{64}$")
    for column in (
        "source_data_sha256",
        "sample_manifest_sha256",
        "decoder_binary_sha256",
        "decoder_parameters_sha256",
        "ne_history_sha256",
    ):
        result[column] = result[column].str.casefold()
        if (
            not result[column]
            .map(lambda value: bool(sha_pattern.fullmatch(value)))
            .all()
        ):
            raise ValueError(f"empirical SHA-256 is malformed: {column}")

    threshold = _numeric_column(result, "threshold_years")
    if not np.array_equal(threshold, np.floor(threshold)):
        raise ValueError("empirical threshold_years must be integer")
    result["threshold_years"] = threshold.astype(np.int64)
    score = _numeric_column(result, "score")
    if np.any((score < 0) | (score > 1)):
        raise ValueError("empirical score must lie in [0, 1]")
    result["score"] = score
    sample_count = _numeric_column(result, "sample_alt_count")
    if np.any(
        (sample_count < 0) | (sample_count > PANEL_HAPLOTYPES)
    ) or not np.array_equal(sample_count, np.floor(sample_count)):
        raise ValueError("empirical sample_alt_count must be an integer in [0, 200]")
    result["sample_alt_count"] = sample_count.astype(np.int64)
    sample_af = _numeric_column(result, "sample_af")
    if not np.allclose(
        sample_af,
        sample_count / PANEL_HAPLOTYPES,
        rtol=0.0,
        atol=1e-12,
    ):
        raise ValueError("empirical sample_af disagrees with sample_alt_count / 200")
    result["sample_af"] = sample_af
    detected = _strict_boolean(result["sample_detected"], label="sample_detected")
    if not np.array_equal(detected, sample_count > 0):
        raise ValueError("empirical sample_detected is inconsistent")
    result["sample_detected"] = detected
    key = ["empirical_id", "source", "window", "threshold_years"]
    if result.duplicated(key, keep=False).any():
        raise ValueError("empirical table contains duplicate statistic rows")
    for labels, curve in result.groupby(
        ["empirical_id", "source", "window"], sort=False
    ):
        ordered = curve.sort_values("threshold_years")
        if tuple(ordered["threshold_years"].astype(int)) != THRESHOLDS_YEARS:
            raise ValueError(f"empirical curve {labels} lacks the seven thresholds")
        if np.any(np.diff(ordered["score"].to_numpy(float)) < -1e-12):
            raise ValueError(f"empirical curve {labels} is nonmonotone")
    provenance = [
        column
        for column in EMPIRICAL_COLUMNS
        if column not in {"window", "threshold_years", "score"}
    ]
    expected_profiles = {("gamma_smc", window) for window in ALLOWED_WINDOWS}
    for empirical_id, group in result.groupby("empirical_id", sort=False):
        profiles = set(zip(group["source"], group["window"], strict=True))
        if profiles != expected_profiles:
            raise ValueError(
                f"empirical_id {empirical_id} lacks both windows and seven "
                "thresholds per window"
            )
        if len(group.loc[:, provenance].drop_duplicates()) != 1:
            raise ValueError(f"empirical provenance changes within {empirical_id}")
    return result.sort_values(key, kind="mergesort").reset_index(drop=True)


def unavailable_empirical_results() -> pd.DataFrame:
    """Return explicit no-claim rows for both planned empirical estimands."""

    rows = []
    for estimand, null_model, _, null_definition in _empirical_null_models():
        row = {column: np.nan for column in EMPIRICAL_RESULT_COLUMNS}
        row.update(
            {
                "status": "unavailable_no_explicit_empirical_table",
                "estimand": estimand,
                "null_model": null_model,
                "null_definition": (
                    "no empirical claim; supply a validated empirical table; "
                    + null_definition
                ),
            }
        )
        rows.append(row)
    return pd.DataFrame(rows, columns=EMPIRICAL_RESULT_COLUMNS)


def unavailable_empirical_omnibus_results() -> pd.DataFrame:
    """Return explicit no-claim omnibus rows for both empirical estimands."""

    rows = []
    for (
        estimand,
        null_model,
        detection_matched,
        null_definition,
    ) in _empirical_null_models():
        row = {column: np.nan for column in EMPIRICAL_OMNIBUS_COLUMNS}
        row.update(
            {
                "status": "unavailable_no_explicit_empirical_table",
                "estimand": estimand,
                "null_model": null_model,
                "n_thresholds": len(THRESHOLDS_YEARS),
                "thresholds_years": ",".join(map(str, THRESHOLDS_YEARS)),
                "sample_detection_matched": detection_matched,
                "null_definition": (
                    "no empirical claim; supply a validated empirical table; "
                    + null_definition
                ),
            }
        )
        rows.append(row)
    return pd.DataFrame(rows, columns=EMPIRICAL_OMNIBUS_COLUMNS)


def _empirical_null_models() -> tuple[tuple[str, str, bool, str], ...]:
    return (
        (
            PRIMARY_EMPIRICAL_ESTIMAND,
            PRIMARY_EMPIRICAL_NULL_MODEL,
            False,
            "all natural neutral population survivors; no final-population-AF, "
            "sample-AF, or sample-detection conditioning",
        ),
        (
            SECONDARY_EMPIRICAL_ESTIMAND,
            SECONDARY_EMPIRICAL_NULL_MODEL,
            True,
            "secondary natural neutral survivors matched on sample detection and "
            "the empirical sample-AF caliper",
        ),
    )


def _validate_af_match_tolerance(value: float) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError("af_match_tolerance must lie in [0, 0.05]") from error
    if not np.isfinite(result) or result < 0 or result > 0.05:
        raise ValueError("af_match_tolerance must lie in [0, 0.05]")
    return result


def _empirical_neutral_subset(
    neutral: pd.DataFrame,
    record: Mapping[str, Any],
    *,
    secondary: bool,
    af_match_tolerance: float,
    include_threshold: bool,
) -> pd.DataFrame:
    mask = (neutral["source"] == record["source"]) & (
        neutral["window"] == record["window"]
    )
    if include_threshold:
        mask &= neutral["threshold_years"] == int(record["threshold_years"])
    if secondary:
        mask &= neutral["sample_detected"] == bool(record["sample_detected"])
        mask &= (
            np.abs(neutral["sample_af"] - float(record["sample_af"]))
            <= af_match_tolerance + 1e-12
        )
    return neutral[mask]


def empirical_neutral_pvalues(
    empirical: pd.DataFrame | None,
    scores: pd.DataFrame,
    *,
    af_match_tolerance: float = 0.025,
) -> pd.DataFrame:
    """Return primary natural-neutral and secondary AF-matched pointwise tests."""

    tolerance = _validate_af_match_tolerance(af_match_tolerance)
    if empirical is None:
        return unavailable_empirical_results()
    frame = validate_replicate_scores(scores)
    data = validate_empirical_input(empirical)
    neutral = frame[frame["simulation_class"] == "neutral"]
    rows = []
    for record in data.to_dict(orient="records"):
        for (
            estimand,
            null_model,
            secondary,
            null_definition,
        ) in _empirical_null_models():
            matched = _empirical_neutral_subset(
                neutral,
                record,
                secondary=secondary,
                af_match_tolerance=tolerance,
                include_threshold=True,
            )
            null = matched["score"].to_numpy(dtype=float)
            output = {
                "status": "available_explicit_empirical_input",
                "estimand": estimand,
                "null_model": null_model,
                "empirical_id": record["empirical_id"],
                "variant_id": record["variant_id"],
                "source": record["source"],
                "window": record["window"],
                "threshold_years": int(record["threshold_years"]),
                "threshold_generations": (
                    float(record["threshold_years"]) / GENERATION_TIME_YEARS
                ),
                "observed_score": float(record["score"]),
                "observed_sample_af": float(record["sample_af"]),
                "af_match_tolerance": tolerance if secondary else np.nan,
                "sample_detection_matched": secondary,
                "n_matched_neutral": len(null),
                "null_definition": null_definition,
            }
            if len(null):
                exceedances = int(
                    np.count_nonzero(null >= float(record["score"]) - 1e-14)
                )
                output.update(
                    {
                        "neutral_upper_exceedance_count": exceedances,
                        "mc_p_upper": float((1 + exceedances) / (len(null) + 1)),
                        "minimum_attainable_p": float(1 / (len(null) + 1)),
                    }
                )
            else:
                output.update(
                    {
                        "status": "unavailable_no_matched_neutral_profiles",
                        "neutral_upper_exceedance_count": np.nan,
                        "mc_p_upper": np.nan,
                        "minimum_attainable_p": np.nan,
                    }
                )
            rows.append(output)
    return pd.DataFrame(rows, columns=EMPIRICAL_RESULT_COLUMNS)


def _exchangeable_min_p_omnibus(
    neutral_profiles: np.ndarray, observed_profile: np.ndarray
) -> dict[str, Any]:
    """Calibrate min pointwise upper-tail p over correlated threshold profiles."""

    pooled = np.vstack((neutral_profiles, observed_profile[None, :]))
    marginal_p = np.empty_like(pooled, dtype=float)
    for index in range(pooled.shape[1]):
        values = pooled[:, index]
        marginal_p[:, index] = np.count_nonzero(
            values[None, :] >= values[:, None] - 1e-14,
            axis=1,
        ) / len(pooled)
    minimum_p = np.min(marginal_p, axis=1)
    observed_minimum = float(minimum_p[-1])
    exceedances = int(np.count_nonzero(minimum_p[:-1] <= observed_minimum + 1e-14))
    return {
        "observed_minimum": observed_minimum,
        "observed_minimum_index": int(np.argmin(marginal_p[-1])),
        "exceedances": exceedances,
        "pvalue": float((1 + exceedances) / len(pooled)),
    }


def empirical_neutral_omnibus_pvalues(
    empirical: pd.DataFrame | None,
    scores: pd.DataFrame,
    *,
    af_match_tolerance: float = 0.025,
) -> pd.DataFrame:
    """Test each empirical seven-threshold curve using shared neutral profiles."""

    tolerance = _validate_af_match_tolerance(af_match_tolerance)
    if empirical is None:
        return unavailable_empirical_omnibus_results()
    frame = validate_replicate_scores(scores)
    data = validate_empirical_input(empirical)
    neutral = frame[frame["simulation_class"] == "neutral"]
    rows = []
    for (empirical_id, source, window), curve in data.groupby(
        ["empirical_id", "source", "window"], sort=True
    ):
        ordered = curve.sort_values("threshold_years")
        record = ordered.iloc[0].to_dict()
        observed = ordered["score"].to_numpy(dtype=float)
        for (
            estimand,
            null_model,
            secondary,
            null_definition,
        ) in _empirical_null_models():
            matched = _empirical_neutral_subset(
                neutral,
                record,
                secondary=secondary,
                af_match_tolerance=tolerance,
                include_threshold=False,
            )
            profiles = matched.pivot(
                index="unit_id", columns="threshold_years", values="score"
            ).reindex(columns=THRESHOLDS_YEARS)
            output = {
                "status": "available_explicit_empirical_input",
                "estimand": estimand,
                "null_model": null_model,
                "empirical_id": empirical_id,
                "variant_id": record["variant_id"],
                "source": source,
                "window": window,
                "n_thresholds": len(THRESHOLDS_YEARS),
                "thresholds_years": ",".join(map(str, THRESHOLDS_YEARS)),
                "n_shared_neutral_profiles": len(profiles),
                "af_match_tolerance": tolerance if secondary else np.nan,
                "sample_detection_matched": secondary,
                "omnibus_method": (
                    "exchangeable_min_pointwise_upper_tail_p_over_shared_"
                    "seven_threshold_neutral_unit_profiles"
                ),
                "null_definition": null_definition,
            }
            if profiles.empty:
                output.update(
                    {
                        "status": "unavailable_no_matched_neutral_profiles",
                        "observed_min_pointwise_p": np.nan,
                        "threshold_years_at_min_pointwise_p": np.nan,
                        "neutral_profile_exceedance_count": np.nan,
                        "omnibus_mc_p_upper": np.nan,
                        "minimum_attainable_p": np.nan,
                    }
                )
            else:
                if profiles.isna().any().any():
                    raise ValueError(
                        "empirical omnibus neutral bank lacks shared complete profiles"
                    )
                inference = _exchangeable_min_p_omnibus(
                    profiles.to_numpy(dtype=float), observed
                )
                best_index = int(inference["observed_minimum_index"])
                output.update(
                    {
                        "observed_min_pointwise_p": float(
                            inference["observed_minimum"]
                        ),
                        "threshold_years_at_min_pointwise_p": (
                            THRESHOLDS_YEARS[best_index]
                        ),
                        "neutral_profile_exceedance_count": int(
                            inference["exceedances"]
                        ),
                        "omnibus_mc_p_upper": float(inference["pvalue"]),
                        "minimum_attainable_p": float(1 / (len(profiles) + 1)),
                    }
                )
            rows.append(output)
    return pd.DataFrame(rows, columns=EMPIRICAL_OMNIBUS_COLUMNS)


def analyze_replicate_scores(
    scores: pd.DataFrame,
    *,
    empirical: pd.DataFrame | None = None,
    bootstrap_draws: int = DEFAULT_BOOTSTRAP_DRAWS,
    permutation_draws: int = DEFAULT_PERMUTATION_DRAWS,
    base_seed: int = BASE_SEED,
    include_af_stratified: bool = True,
    af_stratum_width: float = DEFAULT_AF_STRATUM_WIDTH,
    empirical_af_match_tolerance: float = 0.025,
) -> dict[str, pd.DataFrame]:
    """Run all additive inference without simulating or decoding any unit."""

    base_seed = _nonnegative_integer(base_seed, label="base_seed")
    frame = validate_replicate_scores(scores)
    validated_empirical = (
        empirical_input_template()
        if empirical is None
        else validate_empirical_input(empirical)
    )
    empirical_for_inference = None if empirical is None else validated_empirical
    auc, omnibus = auc_by_threshold(
        frame,
        bootstrap_draws=bootstrap_draws,
        permutation_draws=permutation_draws,
        base_seed=base_seed,
    )
    if include_af_stratified:
        stratified, stratified_omnibus = sample_af_stratified_auc(
            frame,
            af_stratum_width=af_stratum_width,
            bootstrap_draws=bootstrap_draws,
            permutation_draws=permutation_draws,
            base_seed=base_seed,
        )
    else:
        stratified = pd.DataFrame()
        stratified_omnibus = pd.DataFrame()
    return {
        "auc_by_threshold": auc,
        "auc_omnibus": omnibus,
        "selected_unit_pvalues": selected_unit_neutral_tail_pvalues(frame),
        "af_only_benchmarks": af_only_benchmarks(
            frame,
            bootstrap_draws=bootstrap_draws,
            permutation_draws=permutation_draws,
            base_seed=base_seed,
        ),
        "sample_af_stratified_auc": stratified,
        "sample_af_stratified_omnibus": stratified_omnibus,
        "empirical_pvalues": empirical_neutral_pvalues(
            empirical_for_inference,
            frame,
            af_match_tolerance=empirical_af_match_tolerance,
        ),
        "empirical_omnibus": empirical_neutral_omnibus_pvalues(
            empirical_for_inference,
            frame,
            af_match_tolerance=empirical_af_match_tolerance,
        ),
        "empirical_input": validated_empirical,
    }


def _atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        if path.name.endswith(".gz"):
            with temporary.open("wb") as raw:
                with gzip.GzipFile(
                    filename="", fileobj=raw, mode="wb", mtime=0
                ) as compressed:
                    with io.TextIOWrapper(
                        compressed, encoding="utf-8", newline=""
                    ) as text:
                        frame.to_csv(text, sep="\t", index=False, lineterminator="\n")
        else:
            frame.to_csv(temporary, sep="\t", index=False, lineterminator="\n")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def write_analysis_tables(
    results: Mapping[str, pd.DataFrame], output_dir: str | Path
) -> dict[str, Path]:
    """Atomically write the complete additive analysis table set."""

    missing = sorted(set(TABLE_FILENAMES).difference(results))
    if missing:
        raise ValueError("analysis results are missing tables: " + ", ".join(missing))
    destination = Path(output_dir).resolve()
    paths = {}
    for key, filename in TABLE_FILENAMES.items():
        frame = results[key]
        if not isinstance(frame, pd.DataFrame):
            raise TypeError(f"analysis result {key} is not a DataFrame")
        path = destination / filename
        _atomic_tsv(frame, path)
        paths[key] = path
    return paths


def _bootstrap_mean_curve(
    values: np.ndarray, *, draws: int, seed: int
) -> tuple[np.ndarray, np.ndarray]:
    draw_count = _positive_integer(draws, label="bootstrap_draws")
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, len(values), size=(draw_count, len(values)))
    means = values[indices].mean(axis=1)
    return np.quantile(means, 0.025, axis=0), np.quantile(means, 0.975, axis=0)


def _atomic_figure(figure: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        save_options: dict[str, Any] = {
            "format": path.suffix.removeprefix("."),
            "dpi": 300 if path.suffix == ".png" else None,
            "bbox_inches": None,
        }
        if path.suffix == ".pdf":
            save_options["metadata"] = {"CreationDate": None, "ModDate": None}
        figure.savefig(temporary, **save_options)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _save_figure_pair(
    figure: plt.Figure, output_dir: str | Path, stem: str
) -> tuple[Path, Path]:
    destination = Path(output_dir).resolve()
    png = destination / f"{stem}.png"
    pdf = destination / f"{stem}.pdf"
    try:
        _atomic_figure(figure, png)
        _atomic_figure(figure, pdf)
    finally:
        plt.close(figure)
    return png, pdf


def plot_selected_neutral_curves(
    scores: pd.DataFrame,
    output_dir: str | Path,
    *,
    bootstrap_draws: int = DEFAULT_BOOTSTRAP_DRAWS,
    base_seed: int = BASE_SEED,
    stem: str = "han_s001_selected_vs_neutral_tmrca_curves",
) -> tuple[Path, Path]:
    """Plot selected and neutral score curves as a large-text Letter figure."""

    base_seed = _nonnegative_integer(base_seed, label="base_seed")
    frame = validate_replicate_scores(scores)
    colors = {"neutral": "#4C78A8", "selected": "#D1495B"}
    markers = {"neutral": "o", "selected": "s"}
    figure, axes = plt.subplots(
        2, 2, figsize=(11, 8.5), sharex=True, sharey=True, squeeze=False
    )
    x = np.asarray(THRESHOLDS_YEARS, dtype=float) / 1_000
    for row, source in enumerate(ALLOWED_SOURCES):
        for column, window in enumerate(ALLOWED_WINDOWS):
            axis = axes[row, column]
            selected, neutral = _profile_matrices(frame, source=source, window=window)
            for simulation_class, matrix in (
                ("neutral", neutral),
                ("selected", selected),
            ):
                values = matrix.to_numpy(dtype=float)
                low, high = _bootstrap_mean_curve(
                    values,
                    draws=bootstrap_draws,
                    seed=_stable_seed(
                        base_seed,
                        f"plot-mean:{source}:{window}:{simulation_class}",
                    ),
                )
                mean = values.mean(axis=0)
                axis.plot(
                    x,
                    mean,
                    color=colors[simulation_class],
                    marker=markers[simulation_class],
                    linewidth=2.7,
                    markersize=6,
                )
                axis.fill_between(
                    x, low, high, color=colors[simulation_class], alpha=0.18
                )
            axis.set_title(
                f"{source.replace('_', ' ').title()} | {window.replace('_', ' ')}",
                fontsize=15,
            )
            axis.set_ylim(0, 1)
            axis.grid(axis="y", color="#D0D0D0", linewidth=0.8, alpha=0.65)
            axis.tick_params(labelsize=13)
            if row == 1:
                axis.set_xlabel("TMRCA threshold x (kya)", fontsize=15)
            if column == 0:
                axis.set_ylabel("Mean P(TMRCA < x)", fontsize=15)
    figure.suptitle(
        "Han introgression (Q=5): s=0.001 versus natural neutral",
        fontsize=21,
        y=0.985,
    )
    figure.legend(
        handles=[
            Line2D(
                [0],
                [0],
                color=colors["neutral"],
                marker=markers["neutral"],
                linewidth=2.7,
                label="Neutral introgressed region",
            ),
            Line2D(
                [0],
                [0],
                color=colors["selected"],
                marker=markers["selected"],
                linewidth=2.7,
                label="Selected, s=0.001",
            ),
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, 0.94),
        ncol=2,
        frameon=False,
        fontsize=14,
    )
    figure.subplots_adjust(top=0.83, bottom=0.11, left=0.10, right=0.98, hspace=0.28)
    return _save_figure_pair(figure, output_dir, stem)


def plot_auc_profiles(
    auc: pd.DataFrame,
    af_benchmarks: pd.DataFrame,
    output_dir: str | Path,
    *,
    stratified_auc: pd.DataFrame | None = None,
    stem: str = "han_s001_tmrca_roc_auc_profiles",
) -> tuple[Path, Path]:
    """Plot unflipped ROC AUC profiles and the AF-only benchmark."""

    required = {
        "source",
        "window",
        "threshold_years",
        "roc_auc",
        "roc_auc_ci95_low",
        "roc_auc_ci95_high",
    }
    if required.difference(auc.columns):
        raise ValueError("primary AUC table lacks required plotting columns")
    benchmark = af_benchmarks[af_benchmarks["metric"].astype(str) == "sample_af"]
    if len(benchmark) != 1:
        raise ValueError("AF benchmark table must contain exactly one sample_af row")
    af_auc = float(benchmark.iloc[0]["roc_auc"])
    figure, axes = plt.subplots(
        2, 2, figsize=(11, 8.5), sharex=True, sharey=True, squeeze=False
    )
    x = np.asarray(THRESHOLDS_YEARS, dtype=float) / 1_000
    for row, source in enumerate(ALLOWED_SOURCES):
        for column, window in enumerate(ALLOWED_WINDOWS):
            axis = axes[row, column]
            cell = auc[
                (auc["source"] == source) & (auc["window"] == window)
            ].sort_values("threshold_years")
            if tuple(cell["threshold_years"].astype(int)) != THRESHOLDS_YEARS:
                raise ValueError("primary AUC plot lacks a complete threshold curve")
            y = cell["roc_auc"].to_numpy(float)
            low = cell["roc_auc_ci95_low"].to_numpy(float)
            high = cell["roc_auc_ci95_high"].to_numpy(float)
            axis.plot(x, y, color="#D1495B", marker="o", linewidth=2.8)
            axis.fill_between(x, low, high, color="#D1495B", alpha=0.18)
            if stratified_auc is not None and not stratified_auc.empty:
                secondary = stratified_auc[
                    (stratified_auc["source"] == source)
                    & (stratified_auc["window"] == window)
                ].sort_values("threshold_years")
                if tuple(secondary["threshold_years"].astype(int)) != THRESHOLDS_YEARS:
                    raise ValueError("stratified AUC plot lacks a complete curve")
                axis.plot(
                    x,
                    secondary["roc_auc"].to_numpy(float),
                    color="#59A14F",
                    marker="s",
                    linestyle="--",
                    linewidth=2.4,
                )
            axis.axhline(0.5, color="#666666", linewidth=1.5, linestyle=":")
            axis.axhline(af_auc, color="#4C78A8", linewidth=1.7, linestyle="-.")
            axis.set_title(
                f"{source.replace('_', ' ').title()} | {window.replace('_', ' ')}",
                fontsize=15,
            )
            axis.set_ylim(0, 1)
            axis.grid(axis="y", color="#D0D0D0", linewidth=0.8, alpha=0.65)
            axis.tick_params(labelsize=13)
            if row == 1:
                axis.set_xlabel("TMRCA threshold x (kya)", fontsize=15)
            if column == 0:
                axis.set_ylabel("ROC AUC", fontsize=15)
    handles = [
        Line2D(
            [0],
            [0],
            color="#D1495B",
            marker="o",
            linewidth=2.8,
            label="Natural-survivor TMRCA AUC",
        ),
        Line2D(
            [0],
            [0],
            color="#4C78A8",
            linestyle="-.",
            linewidth=1.7,
            label="Sample-AF-only AUC",
        ),
        Line2D(
            [0], [0], color="#666666", linestyle=":", linewidth=1.5, label="AUC = 0.5"
        ),
    ]
    if stratified_auc is not None and not stratified_auc.empty:
        handles.insert(
            1,
            Line2D(
                [0],
                [0],
                color="#59A14F",
                marker="s",
                linestyle="--",
                linewidth=2.4,
                label="Secondary sample-AF-stratified AUC",
            ),
        )
    figure.suptitle(
        "Han selected-versus-neutral TMRCA discrimination (Q=5)",
        fontsize=21,
        y=0.985,
    )
    figure.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.94),
        ncol=2,
        frameon=False,
        fontsize=13,
    )
    figure.subplots_adjust(top=0.81, bottom=0.11, left=0.10, right=0.98, hspace=0.28)
    return _save_figure_pair(figure, output_dir, stem)


__all__ = [
    "ALLOWED_SOURCES",
    "ALLOWED_WINDOWS",
    "DEFAULT_AF_STRATUM_WIDTH",
    "DEFAULT_BOOTSTRAP_DRAWS",
    "DEFAULT_PERMUTATION_DRAWS",
    "EMPIRICAL_COLUMNS",
    "EMPIRICAL_OMNIBUS_COLUMNS",
    "EMPIRICAL_RESULT_COLUMNS",
    "PRIMARY_EMPIRICAL_ESTIMAND",
    "PRIMARY_EMPIRICAL_NULL_MODEL",
    "SCORE_COLUMNS",
    "SCHEMA_VERSION",
    "SECONDARY_EMPIRICAL_ESTIMAND",
    "SECONDARY_EMPIRICAL_NULL_MODEL",
    "TABLE_FILENAMES",
    "THRESHOLDS_YEARS",
    "af_only_benchmarks",
    "analyze_replicate_scores",
    "auc_by_threshold",
    "empirical_input_template",
    "empirical_neutral_omnibus_pvalues",
    "empirical_neutral_pvalues",
    "load_replicate_scores",
    "plot_auc_profiles",
    "plot_selected_neutral_curves",
    "sample_af_stratified_auc",
    "selected_unit_neutral_tail_pvalues",
    "unavailable_empirical_results",
    "unavailable_empirical_omnibus_results",
    "validate_empirical_input",
    "validate_replicate_scores",
    "write_analysis_tables",
]
