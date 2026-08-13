"""Matched-null inference and figures for the focused EAS simulation campaign.

The primary selection score is the within-panel contrast
``P(T < x | hom_alt) - P(T < x | hom_ref)``.  Matching that contrast to the
same statistic in independently simulated neutral panels controls the focal
allele-frequency and demographic conditioning.  Direct genotype-class scores
are retained as secondary between-simulation comparisons.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
import os
import platform
import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.integrate import trapezoid
from scipy.stats import (
    binomtest,
    kstest,
    mannwhitneyu,
    pearsonr,
    rankdata,
    spearmanr,
)

from .calibration import bh_fdr, monte_carlo_pvalue
from .eas_sweep_models import (
    FOCAL_POSITION_BP,
    SEQUENCE_LENGTH_BP,
    TMRCA_THRESHOLDS_YEARS,
)
from .eas_sweep_study import BASE_SEED
from .focused_selection_decode import GENOTYPE_CLASSES, sha256_file

SCHEMA_VERSION = "gamma-smc.focused-selection-analysis/v3"
DEFAULT_BOOTSTRAP_DRAWS = 2_000
DEFAULT_LABEL_PERMUTATION_DRAWS = 9_999
DEFAULT_CROSSFIT_FOLDS = 5
UNIT_COLUMNS = (
    "unit_id",
    "demography_id",
    "simulation_class",
    "selection_coefficient",
    "target_allele_frequency",
    "replicate_index",
)
CELL_COLUMNS = ("demography_id", "target_allele_frequency")
SCOPES = {
    "region_mean": (
        "region_mean_p_tmrca_lt_threshold",
        "region_mean_tmrca_generations",
    ),
    "focal_nearest": (
        "focal_mean_p_tmrca_lt_threshold",
        "focal_mean_tmrca_generations",
    ),
}
PRIMARY_STATISTIC = "p_hom_alt_minus_hom_ref"
EMPIRICAL_STATISTIC_COLUMNS = (
    "demography_id",
    "target_allele_frequency",
    "statistic_scope",
    "metric",
    "statistic",
    "threshold_years",
    "threshold_generations",
    "observed_p_overall",
    "observed_p_hom_ref",
    "observed_p_heterozygous",
    "observed_p_hom_alt",
    "observed_value",
)
EMPIRICAL_PROVENANCE_COLUMNS = (
    "empirical_id",
    "variant_id",
    "source_label",
    "source_data_sha256",
    "sample_manifest_sha256",
    "pair_table_sha256",
    "demography_assignment",
    "population",
    "observed_allele_frequency",
    "af_match_tolerance",
    "n_diploid_samples",
    "n_overall_pairs",
    "n_hom_ref_pairs",
    "n_heterozygous_pairs",
    "n_hom_alt_pairs",
    "alt_allele_count",
    "pairing_scheme",
    "decoder_name",
    "decoder_version",
    "reference_decoder_binary_sha256",
    "empirical_decoder_binary_sha256",
    "reference_decoder_parameters_sha256",
    "empirical_decoder_parameters_sha256",
    "reference_ne_history_sha256",
    "empirical_ne_history_sha256",
    "generation_time_years",
    "mutation_rate_per_bp_per_generation",
    "sequence_length_bp",
    "genome_build",
    "contig",
    "position_1based",
    "ref_allele",
    "alt_allele",
    "region_start_0based",
    "region_end_0based",
    "focal_position_0based",
    "mask_manifest_sha256",
    "qc_status",
)
EMPIRICAL_COLUMNS = (*EMPIRICAL_STATISTIC_COLUMNS, *EMPIRICAL_PROVENANCE_COLUMNS)
TABLE_FILENAMES = {
    "replicate_statistics": "replicate_statistics.tsv.gz",
    "selected_pvalues": "selected_vs_neutral_pvalues.tsv.gz",
    "within_selected": "within_selected_summary.tsv",
    "cell_comparisons": "cell_roc_auc.tsv",
    "neutral_loo": "neutral_loo_pvalues.tsv.gz",
    "neutral_calibration": "neutral_loo_descriptive_calibration.tsv",
    "neutral_crossfit": "neutral_crossfit_pvalues.tsv.gz",
    "neutral_crossfit_calibration": "neutral_crossfit_calibration.tsv",
    "primary_omnibus": "primary_7_threshold_omnibus.tsv.gz",
    "empirical_pvalues": "empirical_pvalues.tsv",
    "empirical_template": "empirical_input_template.tsv",
    "unit_pair_tests": "within_selected_unit_pair_tests.tsv.gz",
    "conditional_power": "conditional_power_summary.tsv",
}
ANALYSIS_IMPLEMENTATION_SOURCES = tuple(
    Path(__file__).resolve().with_name(filename)
    for filename in (
        "focused_selection_analysis.py",
        "calibration.py",
        "eas_sweep_models.py",
        "eas_sweep_study.py",
        "focused_selection_decode.py",
    )
)


def _canonical_sha256(value: Any) -> str:
    payload = json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _stable_seed(base_seed: int, label: str) -> int:
    digest = hashlib.sha256(f"{base_seed}:{label}".encode()).digest()
    return int(1 + int.from_bytes(digest[:8], "big") % (2**32 - 2))


def _analysis_implementation_contract() -> dict[str, Any]:
    sources: dict[str, dict[str, str]] = {}
    for path in ANALYSIS_IMPLEMENTATION_SOURCES:
        if not path.is_file():
            raise ValueError(f"analysis implementation source is absent: {path}")
        sources[path.name] = {
            "path": str(path),
            "sha256": sha256_file(path),
        }
    return {
        "sources": sources,
        "software": {
            "python": platform.python_version(),
            "matplotlib": matplotlib.__version__,
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
        },
    }


def _atomic_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        temporary.write_text(
            json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_frame(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    compression: str | Mapping[str, Any] | None = None
    if path.suffix == ".gz":
        compression = {"method": "gzip", "compresslevel": 6, "mtime": 0}
    try:
        frame.to_csv(
            temporary,
            sep="\t",
            index=False,
            compression=compression,
            float_format="%.12g",
        )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_figure(figure: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        figure.savefig(
            temporary,
            format=path.suffix.removeprefix("."),
            dpi=300 if path.suffix == ".png" else None,
            bbox_inches="tight",
        )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _lower_tail_pvalue(observed: float, null: Sequence[float]) -> float:
    values = np.asarray(null, dtype=float)
    values = values[np.isfinite(values)]
    if not np.isfinite(observed) or not len(values):
        return np.nan
    return float((1 + np.count_nonzero(values <= observed)) / (len(values) + 1))


def _two_sided_pvalue(upper: float, lower: float) -> float:
    if not np.isfinite(upper) or not np.isfinite(lower):
        return np.nan
    return float(min(1.0, 2.0 * min(upper, lower)))


def _roc_auc(selected: Sequence[float], neutral: Sequence[float]) -> float:
    selected_values = np.asarray(selected, dtype=float)
    neutral_values = np.asarray(neutral, dtype=float)
    selected_values = selected_values[np.isfinite(selected_values)]
    neutral_values = neutral_values[np.isfinite(neutral_values)]
    if not len(selected_values) or not len(neutral_values):
        return np.nan
    greater = (selected_values[:, None] > neutral_values[None, :]).sum()
    tied = (selected_values[:, None] == neutral_values[None, :]).sum()
    return float((greater + 0.5 * tied) / (len(selected_values) * len(neutral_values)))


def _genotype_label_permutation_pvalue(
    alt_scores: Sequence[float],
    ref_scores: Sequence[float],
    *,
    seed: int,
    draws: int,
    exact_limit: int = 50_000,
) -> tuple[float, str, int, int]:
    """One-sided label-permutation p-value for mean alt-minus-ref score."""

    alt = np.asarray(alt_scores, dtype=float)
    ref = np.asarray(ref_scores, dtype=float)
    alt = alt[np.isfinite(alt)]
    ref = ref[np.isfinite(ref)]
    if not len(alt) or not len(ref):
        return np.nan, "unavailable", 0, 0
    combined = np.concatenate((alt, ref))
    n_alt = len(alt)
    n_total = len(combined)
    observed = float(np.mean(alt) - np.mean(ref))
    total_assignments = math.comb(n_total, n_alt)
    tolerance = 1e-14
    if total_assignments <= int(exact_limit):
        total_sum = float(combined.sum())
        exceedances = 0
        for indices in itertools.combinations(range(n_total), n_alt):
            alt_sum = float(combined[np.fromiter(indices, dtype=np.int64)].sum())
            statistic = alt_sum / n_alt - (total_sum - alt_sum) / (n_total - n_alt)
            exceedances += statistic >= observed - tolerance
        return (
            float(exceedances / total_assignments),
            "exact_all_label_assignments",
            int(total_assignments),
            int(exceedances),
        )
    if int(draws) < 1:
        raise ValueError("permutation_draws must be positive")
    rng = np.random.default_rng(seed)
    total_sum = float(combined.sum())
    exceedances = 0
    completed = 0
    chunk_size = 256
    while completed < int(draws):
        size = min(chunk_size, int(draws) - completed)
        random_keys = rng.random((size, n_total))
        alt_indices = np.argpartition(random_keys, n_alt - 1, axis=1)[:, :n_alt]
        alt_sum = combined[alt_indices].sum(axis=1)
        statistics = alt_sum / n_alt - (total_sum - alt_sum) / (n_total - n_alt)
        exceedances += int(np.count_nonzero(statistics >= observed - tolerance))
        completed += size
    return (
        float((1 + exceedances) / (int(draws) + 1)),
        "seeded_monte_carlo_plus_one",
        int(draws),
        int(exceedances),
    )


def _unit_label_permutation_pvalues(
    selected: Sequence[float],
    neutral: Sequence[float],
    *,
    selected_indices: np.ndarray,
) -> tuple[float, float, int, int]:
    """Return one-sided mean-difference and AUC label-permutation p-values."""

    selected_values = np.asarray(selected, dtype=float)
    neutral_values = np.asarray(neutral, dtype=float)
    if (
        not len(selected_values)
        or not len(neutral_values)
        or not np.isfinite(selected_values).all()
        or not np.isfinite(neutral_values).all()
    ):
        return np.nan, np.nan, 0, 0
    indices = np.asarray(selected_indices, dtype=np.int64)
    if indices.ndim != 2 or indices.shape[1] != len(selected_values):
        raise ValueError("selected permutation-index matrix has the wrong shape")
    combined = np.concatenate((selected_values, neutral_values))
    if indices.size and (indices.min() < 0 or indices.max() >= len(combined)):
        raise ValueError("selected permutation-index matrix is out of bounds")
    n_selected = len(selected_values)
    n_neutral = len(neutral_values)
    observed_mean = float(np.mean(selected_values) - np.mean(neutral_values))
    observed_auc = _roc_auc(selected_values, neutral_values)
    selected_sums = combined[indices].sum(axis=1)
    permuted_means = (
        selected_sums / n_selected - (combined.sum() - selected_sums) / n_neutral
    )
    ranks = rankdata(combined, method="average")
    permuted_rank_sums = ranks[indices].sum(axis=1)
    permuted_auc = (permuted_rank_sums - n_selected * (n_selected + 1) / 2) / (
        n_selected * n_neutral
    )
    tolerance = 1e-14
    mean_exceedances = int(
        np.count_nonzero(permuted_means >= observed_mean - tolerance)
    )
    auc_exceedances = int(np.count_nonzero(permuted_auc >= observed_auc - tolerance))
    denominator = len(indices) + 1
    return (
        float((1 + mean_exceedances) / denominator),
        float((1 + auc_exceedances) / denominator),
        mean_exceedances,
        auc_exceedances,
    )


def _wilson_interval(detected: int, total: int) -> tuple[float, float, float]:
    if int(total) < 1 or not 0 <= int(detected) <= int(total):
        return np.nan, np.nan, np.nan
    fraction = int(detected) / int(total)
    z = 1.959963984540054
    denominator = 1 + z**2 / total
    center = (fraction + z**2 / (2 * total)) / denominator
    half_width = (
        z
        * np.sqrt(fraction * (1 - fraction) / total + z**2 / (4 * total**2))
        / denominator
    )
    return (
        float(fraction),
        float(max(0.0, center - half_width)),
        float(min(1.0, center + half_width)),
    )


def _bootstrap_auc_interval(
    selected: Sequence[float],
    neutral: Sequence[float],
    *,
    draws: int,
    seed: int,
) -> tuple[float, float]:
    selected_values = np.asarray(selected, dtype=float)
    neutral_values = np.asarray(neutral, dtype=float)
    selected_values = selected_values[np.isfinite(selected_values)]
    neutral_values = neutral_values[np.isfinite(neutral_values)]
    if draws < 1 or not len(selected_values) or not len(neutral_values):
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    values = np.empty(draws, dtype=float)
    chunk = 100
    for start in range(0, draws, chunk):
        stop = min(draws, start + chunk)
        size = stop - start
        selected_draws = selected_values[
            rng.integers(0, len(selected_values), size=(size, len(selected_values)))
        ]
        neutral_draws = neutral_values[
            rng.integers(0, len(neutral_values), size=(size, len(neutral_values)))
        ]
        greater = (selected_draws[:, :, None] > neutral_draws[:, None, :]).sum(
            axis=(1, 2)
        )
        tied = (selected_draws[:, :, None] == neutral_draws[:, None, :]).sum(
            axis=(1, 2)
        )
        values[start:stop] = (greater + 0.5 * tied) / (
            len(selected_values) * len(neutral_values)
        )
    low, high = np.quantile(values, [0.025, 0.975])
    return float(low), float(high)


def _bootstrap_mean_interval(
    values: Sequence[float], *, draws: int, seed: int
) -> tuple[float, float]:
    array = np.asarray(values, dtype=float)
    array = array[np.isfinite(array)]
    if draws < 1 or not len(array):
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, len(array), size=(draws, len(array)))
    means = array[indices].mean(axis=1)
    low, high = np.quantile(means, [0.025, 0.975])
    return float(low), float(high)


def _bootstrap_difference_interval(
    selected: Sequence[float],
    neutral: Sequence[float],
    *,
    draws: int,
    seed: int,
) -> tuple[float, float]:
    selected_values = np.asarray(selected, dtype=float)
    neutral_values = np.asarray(neutral, dtype=float)
    selected_values = selected_values[np.isfinite(selected_values)]
    neutral_values = neutral_values[np.isfinite(neutral_values)]
    if draws < 1 or not len(selected_values) or not len(neutral_values):
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    selected_indices = rng.integers(
        0, len(selected_values), size=(draws, len(selected_values))
    )
    neutral_indices = rng.integers(
        0, len(neutral_values), size=(draws, len(neutral_values))
    )
    differences = selected_values[selected_indices].mean(axis=1) - neutral_values[
        neutral_indices
    ].mean(axis=1)
    low, high = np.quantile(differences, [0.025, 0.975])
    return float(low), float(high)


def validate_class_summaries(
    frame: pd.DataFrame,
    *,
    thresholds_years: Sequence[float] = TMRCA_THRESHOLDS_YEARS,
    minimum_selected_per_cell: int = 10,
    minimum_neutral_per_cell: int = 100,
) -> pd.DataFrame:
    """Validate the complete six-cell table before any inferential analysis."""

    required = {
        *UNIT_COLUMNS,
        "source",
        "genotype_class",
        "threshold_years",
        "n_pairs",
        *(column for columns in SCOPES.values() for column in columns),
    }
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"class summaries are missing columns: {', '.join(missing)}")
    if frame.empty:
        raise ValueError("class summaries are empty")
    result = frame.copy()
    if result[[*UNIT_COLUMNS, "source"]].isna().any().any():
        raise ValueError("class summaries contain missing unit/cell identifiers")
    result["simulation_class"] = result["simulation_class"].astype(str)
    if set(result["simulation_class"]) - {"selected", "neutral"}:
        raise ValueError("simulation_class must contain only selected and neutral")
    result["genotype_class"] = result["genotype_class"].astype(str)
    if set(result["genotype_class"]) != set(GENOTYPE_CLASSES):
        raise ValueError("class summaries must contain all four genotype classes")
    thresholds = tuple(float(value) for value in thresholds_years)
    if set(result["threshold_years"].astype(float)) != set(thresholds):
        raise ValueError("class summaries do not contain the required thresholds")
    numeric_columns = [
        "target_allele_frequency",
        "selection_coefficient",
        "replicate_index",
        "threshold_years",
        "n_pairs",
        *(column for columns in SCOPES.values() for column in columns),
    ]
    for column in numeric_columns:
        result[column] = pd.to_numeric(result[column], errors="coerce")
        if not np.isfinite(result[column].to_numpy(dtype=float)).all():
            raise ValueError(f"class-summary column is not finite: {column}")
    if (result["n_pairs"] != np.floor(result["n_pairs"])).any():
        raise ValueError("class-summary n_pairs values must be integers")
    probability_columns = [columns[0] for columns in SCOPES.values()]
    if any(
        ((result[column] < 0) | (result[column] > 1)).any()
        for column in probability_columns
    ):
        raise ValueError("P(TMRCA < x) values must lie in [0, 1]")
    tmrca_columns = [columns[1] for columns in SCOPES.values()]
    if any((result[column] <= 0).any() for column in tmrca_columns):
        raise ValueError("mean TMRCA values must be positive")
    row_key = [*UNIT_COLUMNS, "genotype_class", "threshold_years"]
    if result.duplicated(row_key, keep=False).any():
        raise ValueError("class summaries contain duplicate unit/class/threshold rows")
    expected_rows = len(GENOTYPE_CLASSES) * len(thresholds)
    unit_sizes = result.groupby("unit_id", sort=False).size()
    if (unit_sizes != expected_rows).any():
        raise ValueError(
            "every unit must contain each genotype class and threshold once"
        )
    for unit_id, group in result.groupby("unit_id", sort=False):
        metadata = group[list(UNIT_COLUMNS)].drop_duplicates()
        if len(metadata) != 1:
            raise ValueError(f"unit metadata changes within {unit_id}")
        if group["source"].astype(str).nunique() != 1:
            raise ValueError(f"analysis source changes within {unit_id}")
        pair_counts = (
            group.groupby("genotype_class", sort=False)["n_pairs"]
            .agg(["min", "max"])
            .astype(float)
        )
        if (pair_counts["min"] != pair_counts["max"]).any() or (
            pair_counts["min"] <= 0
        ).any():
            raise ValueError(f"genotype pair counts are invalid within {unit_id}")
        pair_count = pair_counts["min"]
        if not np.isclose(
            pair_count["overall"],
            pair_count[["hom_ref", "heterozygous", "hom_alt"]].sum(),
        ):
            raise ValueError(f"genotype pair counts do not sum within {unit_id}")
        for column in tmrca_columns:
            spread = group.groupby("genotype_class")[column].nunique(dropna=False)
            if (spread != 1).any():
                raise ValueError(
                    f"threshold-invariant TMRCA values change within {unit_id}"
                )
        for column in probability_columns:
            for _, curve in group.groupby("genotype_class", sort=False):
                values = curve.sort_values("threshold_years")[column].to_numpy(
                    dtype=float
                )
                if np.any(np.diff(values) < -1e-10):
                    raise ValueError(
                        f"P(TMRCA < x) is not nondecreasing within {unit_id}"
                    )
    unit_table = result[list(UNIT_COLUMNS)].drop_duplicates()
    selected_coefficients = set(
        unit_table.loc[
            unit_table["simulation_class"] == "selected", "selection_coefficient"
        ].astype(float)
    )
    if not selected_coefficients or any(value <= 0 for value in selected_coefficients):
        raise ValueError("selected units must have a positive selection coefficient")
    if set(
        unit_table.loc[
            unit_table["simulation_class"] == "neutral", "selection_coefficient"
        ].astype(float)
    ) != {0.0}:
        raise ValueError("shared neutral units must have selection_coefficient 0")
    neutral_counts = (
        unit_table[unit_table["simulation_class"] == "neutral"]
        .groupby(list(CELL_COLUMNS))
        .size()
    )
    selected_counts = (
        unit_table[unit_table["simulation_class"] == "selected"]
        .groupby([*CELL_COLUMNS, "selection_coefficient"])
        .size()
    )
    cells = unit_table[list(CELL_COLUMNS)].drop_duplicates()
    for cell in cells.itertuples(index=False, name=None):
        neutral_count = int(neutral_counts.get(cell, 0))
        if neutral_count < int(minimum_neutral_per_cell):
            raise ValueError(
                f"cell {cell} has {neutral_count} neutral units; at least "
                f"{minimum_neutral_per_cell} are required"
            )
        for coefficient in sorted(selected_coefficients):
            count = int(selected_counts.get((*cell, coefficient), 0))
            if count < int(minimum_selected_per_cell):
                raise ValueError(
                    f"cell {cell}, s={coefficient:g} has {count} selected units; "
                    f"at least {minimum_selected_per_cell} are required"
                )
    return result


def validate_spatial_summaries(
    frame: pd.DataFrame,
    class_summaries: pd.DataFrame,
) -> pd.DataFrame:
    """Validate cell-level genomic profiles against the complete unit table."""

    required = {
        "source",
        "demography_id",
        "simulation_class",
        "selection_coefficient",
        "target_allele_frequency",
        "genotype_class",
        "threshold_years",
        "position_0based",
        "n_units",
        "mean_p_tmrca_lt_threshold",
        "sd_p_tmrca_lt_threshold",
        "sem_p_tmrca_lt_threshold",
        "mean_tmrca_generations",
        "sd_tmrca_generations",
        "sem_tmrca_generations",
    }
    missing = sorted(required.difference(frame.columns))
    if missing or frame.empty:
        detail = ", ".join(missing) if missing else "empty table"
        raise ValueError(f"aggregated spatial summaries are invalid: {detail}")
    result = frame.copy()
    numeric = [
        "selection_coefficient",
        "target_allele_frequency",
        "threshold_years",
        "position_0based",
        "n_units",
        "mean_p_tmrca_lt_threshold",
        "sd_p_tmrca_lt_threshold",
        "sem_p_tmrca_lt_threshold",
        "mean_tmrca_generations",
        "sd_tmrca_generations",
        "sem_tmrca_generations",
    ]
    for column in numeric:
        result[column] = pd.to_numeric(result[column], errors="coerce")
        if not np.isfinite(result[column].to_numpy(dtype=float)).all():
            raise ValueError(f"spatial-summary column is not finite: {column}")
    if set(result["genotype_class"].astype(str)) != set(GENOTYPE_CLASSES):
        raise ValueError("spatial summaries must contain all four genotype classes")
    if set(result["threshold_years"].astype(float)) != set(TMRCA_THRESHOLDS_YEARS):
        raise ValueError("spatial summaries lack a prespecified threshold")
    positions = result["position_0based"].to_numpy(dtype=float)
    if (
        not np.array_equal(positions, positions.astype(np.int64))
        or np.any(positions < 0)
        or np.any(positions >= SEQUENCE_LENGTH_BP)
        or FOCAL_POSITION_BP not in set(positions.astype(np.int64))
    ):
        raise ValueError("spatial output positions are outside the 10-Mb contract")
    if (
        (result["n_units"] < 1).any()
        or (result["n_units"] != np.floor(result["n_units"])).any()
        or (result["mean_p_tmrca_lt_threshold"] < 0).any()
        or (result["mean_p_tmrca_lt_threshold"] > 1).any()
        or (result["mean_tmrca_generations"] <= 0).any()
        or (
            result[[column for column in numeric if column.startswith(("sd_", "sem_"))]]
            < 0
        )
        .any()
        .any()
    ):
        raise ValueError("spatial means, uncertainty, or unit counts are invalid")
    key = [
        "demography_id",
        "simulation_class",
        "selection_coefficient",
        "target_allele_frequency",
        "genotype_class",
        "threshold_years",
        "position_0based",
    ]
    if result.duplicated(key).any():
        raise ValueError("spatial summaries contain duplicate cell/profile rows")
    class_sources = set(class_summaries["source"].astype(str))
    if len(class_sources) != 1 or set(result["source"].astype(str)) != class_sources:
        raise ValueError("spatial and class-summary analysis sources disagree")
    unit_table = class_summaries[list(UNIT_COLUMNS)].drop_duplicates()
    expected_counts = (
        unit_table.groupby(
            [
                "demography_id",
                "simulation_class",
                "selection_coefficient",
                "target_allele_frequency",
            ],
            dropna=False,
        )
        .size()
        .to_dict()
    )
    observed_groups = result.groupby(key[:-1], sort=False, dropna=False)
    expected_positions: np.ndarray | None = None
    observed_keys: set[tuple[Any, ...]] = set()
    for labels, group in observed_groups:
        observed_keys.add(tuple(labels))
        cell = tuple(labels[:4])
        expected_n = expected_counts.get(cell)
        if expected_n is None or set(group["n_units"].astype(int)) != {int(expected_n)}:
            raise ValueError(
                "spatial profile unit count disagrees with class summaries"
            )
        group_positions = np.sort(group["position_0based"].to_numpy(dtype=np.int64))
        if expected_positions is None:
            expected_positions = group_positions
        elif not np.array_equal(expected_positions, group_positions):
            raise ValueError("spatial output grid changes between cells")
    expected_keys = {
        (
            str(row.demography_id),
            str(row.simulation_class),
            float(row.selection_coefficient),
            float(row.target_allele_frequency),
            genotype_class,
            float(threshold),
        )
        for row in unit_table[
            [
                "demography_id",
                "simulation_class",
                "selection_coefficient",
                "target_allele_frequency",
            ]
        ]
        .drop_duplicates()
        .itertuples(index=False)
        for genotype_class in GENOTYPE_CLASSES
        for threshold in TMRCA_THRESHOLDS_YEARS
    }
    normalized_observed = {
        (
            str(values[0]),
            str(values[1]),
            float(values[2]),
            float(values[3]),
            str(values[4]),
            float(values[5]),
        )
        for values in observed_keys
    }
    if normalized_observed != expected_keys:
        raise ValueError("aggregated spatial profile grid is incomplete")
    return result


def build_replicate_statistics(class_summaries: pd.DataFrame) -> pd.DataFrame:
    """Create selection-directed unit scores from compact class summaries."""

    records: list[dict[str, Any]] = []
    metadata_columns = list(UNIT_COLUMNS)
    for scope, (probability_column, tmrca_column) in SCOPES.items():
        for unit_id, unit in class_summaries.groupby("unit_id", sort=True):
            metadata = unit.iloc[0][metadata_columns].to_dict()
            overall_pairs = unit.loc[
                unit["genotype_class"] == "overall", "n_pairs"
            ].astype(int)
            metadata["panel_n_diploid_samples"] = int(overall_pairs.iloc[0])
            metadata["analysis_source"] = str(unit.iloc[0]["source"])
            probability = unit.pivot(
                index="threshold_years",
                columns="genotype_class",
                values=probability_column,
            ).sort_index()
            for threshold, row in probability.iterrows():
                for genotype_class in GENOTYPE_CLASSES:
                    records.append(
                        {
                            **metadata,
                            "statistic_scope": scope,
                            "metric": "p_tmrca_lt_threshold",
                            "statistic": f"p_{genotype_class}",
                            "threshold_years": float(threshold),
                            "selection_score": float(row[genotype_class]),
                            "raw_value": float(row[genotype_class]),
                            "selection_direction": "larger_is_more_recent",
                        }
                    )
                contrast = float(row["hom_alt"] - row["hom_ref"])
                records.append(
                    {
                        **metadata,
                        "statistic_scope": scope,
                        "metric": "p_tmrca_lt_threshold",
                        "statistic": PRIMARY_STATISTIC,
                        "threshold_years": float(threshold),
                        "selection_score": contrast,
                        "raw_value": contrast,
                        "selection_direction": "larger_is_younger_alt_vs_ref",
                    }
                )

            thresholds = probability.index.to_numpy(dtype=float)
            x = np.concatenate(([0.0], thresholds))
            areas: dict[str, float] = {}
            for genotype_class in GENOTYPE_CLASSES:
                y = np.concatenate(
                    ([0.0], probability[genotype_class].to_numpy(dtype=float))
                )
                areas[genotype_class] = float(trapezoid(y, x) / x[-1])
                records.append(
                    {
                        **metadata,
                        "statistic_scope": scope,
                        "metric": "normalized_cdf_area",
                        "statistic": f"cdf_area_{genotype_class}",
                        "threshold_years": 0.0,
                        "selection_score": areas[genotype_class],
                        "raw_value": areas[genotype_class],
                        "selection_direction": "larger_is_more_recent",
                    }
                )
            area_contrast = areas["hom_alt"] - areas["hom_ref"]
            records.append(
                {
                    **metadata,
                    "statistic_scope": scope,
                    "metric": "normalized_cdf_area",
                    "statistic": "cdf_area_hom_alt_minus_hom_ref",
                    "threshold_years": 0.0,
                    "selection_score": area_contrast,
                    "raw_value": area_contrast,
                    "selection_direction": "larger_is_younger_alt_vs_ref",
                }
            )

            first_threshold = float(probability.index[0])
            tmrca = (
                unit[unit["threshold_years"].astype(float) == first_threshold]
                .set_index("genotype_class")[tmrca_column]
                .astype(float)
            )
            tmrca_contrast = float(tmrca["hom_ref"] - tmrca["hom_alt"])
            records.append(
                {
                    **metadata,
                    "statistic_scope": scope,
                    "metric": "mean_tmrca_generations",
                    "statistic": "tmrca_hom_ref_minus_hom_alt",
                    "threshold_years": 0.0,
                    "selection_score": tmrca_contrast,
                    "raw_value": tmrca_contrast,
                    "selection_direction": "larger_is_younger_alt_vs_ref",
                }
            )
    result = pd.DataFrame(records)
    key = [
        *UNIT_COLUMNS,
        "statistic_scope",
        "metric",
        "statistic",
        "threshold_years",
    ]
    if result.duplicated(key, keep=False).any():
        raise RuntimeError("replicate-statistic construction produced duplicate rows")
    return result


def selected_vs_neutral_pvalues(scores: pd.DataFrame) -> pd.DataFrame:
    """Score every selected replicate against its independently seeded null."""

    group_columns = [
        *CELL_COLUMNS,
        "statistic_scope",
        "metric",
        "statistic",
        "threshold_years",
    ]
    rows: list[dict[str, Any]] = []
    for labels, group in scores.groupby(group_columns, sort=True, dropna=False):
        label = dict(zip(group_columns, labels, strict=True))
        null = group.loc[
            group["simulation_class"] == "neutral", "selection_score"
        ].to_numpy(dtype=float)
        selected = group[group["simulation_class"] == "selected"]
        for record in selected.to_dict(orient="records"):
            observed = float(record["selection_score"])
            upper = monte_carlo_pvalue(observed, null)
            lower = _lower_tail_pvalue(observed, null)
            rows.append(
                {
                    **{column: record[column] for column in UNIT_COLUMNS},
                    **label,
                    "observed_selection_score": observed,
                    "null_mean": float(np.mean(null)),
                    "null_median": float(np.median(null)),
                    "null_sd": float(np.std(null, ddof=1)) if len(null) > 1 else np.nan,
                    "null_q025": float(np.quantile(null, 0.025)),
                    "null_q05": float(np.quantile(null, 0.05)),
                    "null_q95": float(np.quantile(null, 0.95)),
                    "null_q975": float(np.quantile(null, 0.975)),
                    "n_null": len(null),
                    "null_upper_exceedance_count": int(
                        np.count_nonzero(null >= observed)
                    ),
                    "mc_p_upper": upper,
                    "mc_p_lower": lower,
                    "mc_p_two_sided": _two_sided_pvalue(upper, lower),
                    "null_percentile_midrank": float(
                        (
                            np.count_nonzero(null < observed)
                            + 0.5 * np.count_nonzero(null == observed)
                        )
                        / len(null)
                    ),
                    "neutral_bank_scope": (
                        "shared_by_demography_and_af_across_selection_coefficients"
                    ),
                }
            )
    result = pd.DataFrame(rows)
    if result.empty:
        return result
    family = [
        "unit_id",
        "selection_coefficient",
        "statistic_scope",
        "metric",
        "statistic",
    ]
    result["bh_q_upper"] = result.groupby(family, sort=False)["mc_p_upper"].transform(
        lambda values: bh_fdr(values.to_numpy(dtype=float))
    )
    result["bh_q_two_sided"] = result.groupby(family, sort=False)[
        "mc_p_two_sided"
    ].transform(lambda values: bh_fdr(values.to_numpy(dtype=float)))
    return result


def primary_threshold_omnibus(scores: pd.DataFrame) -> pd.DataFrame:
    """Correlation-aware one-sided minP test across the seven primary thresholds.

    For each selected unit, its seven-score vector is pooled with the matched
    neutral unit vectors. Every pooled row is left out in turn when its seven
    marginal upper-tail p-values are calculated. Ranking the observed minimum
    p-value among these exchangeable leave-one-out minima preserves the neutral
    cross-threshold dependence and yields an exact finite-bank permutation p-value.
    """

    thresholds = tuple(float(value) for value in TMRCA_THRESHOLDS_YEARS)
    primary = scores[
        (scores["metric"] == "p_tmrca_lt_threshold")
        & (scores["statistic"] == PRIMARY_STATISTIC)
        & scores["threshold_years"].astype(float).isin(thresholds)
    ].copy()
    group_columns = [*CELL_COLUMNS, "statistic_scope"]
    rows: list[dict[str, Any]] = []
    for labels, group in primary.groupby(group_columns, sort=True, dropna=False):
        label = dict(zip(group_columns, labels, strict=True))
        neutral = (
            group[group["simulation_class"] == "neutral"]
            .pivot(index="unit_id", columns="threshold_years", values="selection_score")
            .reindex(columns=thresholds)
            .sort_index()
        )
        if neutral.empty or neutral.isna().any().any():
            raise ValueError("primary omnibus neutral vectors are incomplete")
        neutral_values = neutral.to_numpy(dtype=float)
        selected = group[group["simulation_class"] == "selected"]
        for coefficient, coefficient_group in selected.groupby(
            "selection_coefficient", sort=True
        ):
            selected_wide = (
                coefficient_group.pivot(
                    index="unit_id",
                    columns="threshold_years",
                    values="selection_score",
                )
                .reindex(columns=thresholds)
                .sort_index()
            )
            if selected_wide.isna().any().any():
                raise ValueError("primary omnibus selected vectors are incomplete")
            for unit_id, observed_row in selected_wide.iterrows():
                unit_metadata = coefficient_group[
                    coefficient_group["unit_id"].astype(str) == str(unit_id)
                ].iloc[0]
                observed = observed_row.to_numpy(dtype=float)
                pooled = np.vstack((neutral_values, observed))
                marginal = np.empty_like(pooled, dtype=float)
                for threshold_index in range(len(thresholds)):
                    values = pooled[:, threshold_index]
                    upper_counts_including_self = (
                        values[None, :] >= values[:, None] - 1e-14
                    ).sum(axis=1)
                    marginal[:, threshold_index] = upper_counts_including_self / len(
                        pooled
                    )
                min_p = marginal.min(axis=1)
                observed_min_p = float(min_p[-1])
                neutral_exceedances = int(
                    np.count_nonzero(min_p[:-1] <= observed_min_p + 1e-14)
                )
                best_index = int(np.argmin(marginal[-1]))
                rows.append(
                    {
                        **label,
                        "simulation_class": "selected",
                        "selection_coefficient": float(coefficient),
                        "unit_id": str(unit_id),
                        "replicate_index": int(unit_metadata["replicate_index"]),
                        "metric": "p_tmrca_lt_threshold",
                        "statistic": PRIMARY_STATISTIC,
                        "n_null": len(neutral_values),
                        "n_thresholds": len(thresholds),
                        "thresholds_years": ",".join(
                            f"{value:g}" for value in thresholds
                        ),
                        "observed_min_marginal_p_upper": observed_min_p,
                        "most_extreme_threshold_years": thresholds[best_index],
                        "observed_score_at_most_extreme_threshold": float(
                            observed[best_index]
                        ),
                        "neutral_minp_exceedance_count": neutral_exceedances,
                        "omnibus_p_upper": float(
                            (1 + neutral_exceedances) / len(pooled)
                        ),
                        "omnibus_method": (
                            "pooled_leave_one_out_minP_exact_unit_exchangeability"
                        ),
                        "correlation_handling": (
                            "whole neutral seven-threshold vectors retained"
                        ),
                        "neutral_bank_scope": (
                            "shared_by_demography_and_af_across_selection_coefficients"
                        ),
                    }
                )
    return pd.DataFrame(rows)


def conditional_power_summary(
    selected_pvalues: pd.DataFrame, *, alpha: float = 0.05
) -> pd.DataFrame:
    """Estimate raw and within-seven-BH AF-conditioned detection power."""

    if not 0 < float(alpha) < 1:
        raise ValueError("alpha must lie strictly between zero and one")
    if selected_pvalues.empty:
        return pd.DataFrame()
    group_columns = [
        *CELL_COLUMNS,
        "selection_coefficient",
        "statistic_scope",
        "metric",
        "statistic",
        "threshold_years",
    ]
    rows: list[dict[str, Any]] = []
    for labels, group in selected_pvalues.groupby(
        group_columns, sort=True, dropna=False
    ):
        raw_pvalues = group["mc_p_upper"].to_numpy(dtype=float)
        adjusted_pvalues = group["bh_q_upper"].to_numpy(dtype=float)
        valid = np.isfinite(raw_pvalues) & np.isfinite(adjusted_pvalues)
        raw_pvalues = raw_pvalues[valid]
        adjusted_pvalues = adjusted_pvalues[valid]
        n = len(raw_pvalues)
        raw_detected = int(np.count_nonzero(raw_pvalues <= alpha))
        adjusted_detected = int(np.count_nonzero(adjusted_pvalues <= alpha))
        raw_fraction, raw_low, raw_high = _wilson_interval(raw_detected, n)
        adjusted_fraction, adjusted_low, adjusted_high = _wilson_interval(
            adjusted_detected, n
        )
        metric = str(labels[group_columns.index("metric")])
        family_size = (
            len(TMRCA_THRESHOLDS_YEARS) if metric == "p_tmrca_lt_threshold" else 1
        )
        rows.append(
            {
                **dict(zip(group_columns, labels, strict=True)),
                "alpha": float(alpha),
                "n_selected": n,
                # Backward-compatible aliases are explicitly the raw result.
                "n_detected": raw_detected,
                "conditional_power": raw_fraction,
                "wilson_ci95_low": raw_low,
                "wilson_ci95_high": raw_high,
                "n_detected_raw": raw_detected,
                "conditional_power_raw": raw_fraction,
                "raw_wilson_ci95_low": raw_low,
                "raw_wilson_ci95_high": raw_high,
                "n_detected_bh_within_7": adjusted_detected,
                "conditional_power_bh_within_7": adjusted_fraction,
                "bh_within_7_wilson_ci95_low": adjusted_low,
                "bh_within_7_wilson_ci95_high": adjusted_high,
                "bh_family_size": int(family_size),
                "bh_family_definition": (
                    "thresholds within unit, selection coefficient, scope, metric, statistic"
                ),
                "conditioning": ("conditional_on_accepted_final_af_and_genotype_qc"),
                "neutral_bank_scope": (
                    "shared_by_demography_and_af_across_selection_coefficients"
                ),
            }
        )
    return pd.DataFrame(rows)


def neutral_leave_one_out(scores: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return descriptive-only LOO neutral ranks and calibration summaries."""

    neutral = scores[scores["simulation_class"] == "neutral"].copy()
    group_columns = [
        *CELL_COLUMNS,
        "statistic_scope",
        "metric",
        "statistic",
        "threshold_years",
    ]
    rows: list[dict[str, Any]] = []
    for labels, group in neutral.groupby(group_columns, sort=True, dropna=False):
        label = dict(zip(group_columns, labels, strict=True))
        values = group["selection_score"].to_numpy(dtype=float)
        unit_ids = group["unit_id"].astype(str).to_numpy()
        for index, (unit_id, observed) in enumerate(zip(unit_ids, values, strict=True)):
            null = np.delete(values, index)
            rows.append(
                {
                    **label,
                    "unit_id": unit_id,
                    "observed_selection_score": float(observed),
                    "n_loo_null": len(null),
                    "loo_mc_p_upper": monte_carlo_pvalue(observed, null),
                    "diagnostic_role": "descriptive_only_not_out_of_sample",
                }
            )
    loo = pd.DataFrame(rows)
    calibration_rows: list[dict[str, Any]] = []
    for labels, group in loo.groupby(group_columns, sort=True, dropna=False):
        label = dict(zip(group_columns, labels, strict=True))
        pvalues = group["loo_mc_p_upper"].to_numpy(dtype=float)
        pvalues = pvalues[np.isfinite(pvalues)]
        calibration_rows.append(
            {
                **label,
                "n_neutral": len(pvalues),
                "mean_loo_p": float(np.mean(pvalues)),
                "fraction_loo_p_lt_0_05": float(np.mean(pvalues < 0.05)),
                "ks_uniform_p": float(kstest(pvalues, "uniform").pvalue),
                "calibration_note": (
                    "descriptive only; LOO reuses the neutral bank and finite-null "
                    "p-values are discrete"
                ),
            }
        )
    return loo, pd.DataFrame(calibration_rows)


def neutral_crossfit_calibration(
    scores: pd.DataFrame,
    *,
    n_folds: int = DEFAULT_CROSSFIT_FOLDS,
    base_seed: int = BASE_SEED,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Evaluate neutral calibration with deterministic held-out unit folds."""

    if int(n_folds) < 2:
        raise ValueError("n_folds must be at least two")
    neutral = scores[scores["simulation_class"] == "neutral"].copy()
    unit_table = neutral[[*CELL_COLUMNS, "unit_id"]].drop_duplicates()
    fold_rows: list[dict[str, Any]] = []
    for labels, group in unit_table.groupby(list(CELL_COLUMNS), sort=True):
        unit_ids = np.asarray(sorted(group["unit_id"].astype(str)), dtype=object)
        if len(unit_ids) < int(n_folds):
            raise ValueError(
                f"neutral cross-fit cell {labels} has fewer units than folds"
            )
        seed_label = ":".join(map(str, labels))
        seed = _stable_seed(base_seed, "neutral-crossfit-folds:" + seed_label)
        rng = np.random.default_rng(seed)
        shuffled = unit_ids[rng.permutation(len(unit_ids))]
        for order, unit_id in enumerate(shuffled):
            fold_rows.append(
                {
                    **dict(zip(CELL_COLUMNS, labels, strict=True)),
                    "unit_id": str(unit_id),
                    "crossfit_fold": int(order % int(n_folds)),
                    "crossfit_assignment_seed": int(seed),
                }
            )
    assignments = pd.DataFrame(fold_rows)
    neutral = neutral.merge(
        assignments,
        on=[*CELL_COLUMNS, "unit_id"],
        how="left",
        validate="many_to_one",
    )
    group_columns = [
        *CELL_COLUMNS,
        "statistic_scope",
        "metric",
        "statistic",
        "threshold_years",
    ]
    rows: list[dict[str, Any]] = []
    for labels, group in neutral.groupby(group_columns, sort=True, dropna=False):
        label = dict(zip(group_columns, labels, strict=True))
        for fold, held_out in group.groupby("crossfit_fold", sort=True):
            training = group[group["crossfit_fold"] != fold]
            null = training["selection_score"].to_numpy(dtype=float)
            if not len(null):
                raise ValueError("neutral cross-fit training bank is empty")
            for record in held_out.to_dict(orient="records"):
                observed = float(record["selection_score"])
                upper = monte_carlo_pvalue(observed, null)
                lower = _lower_tail_pvalue(observed, null)
                rows.append(
                    {
                        **label,
                        "unit_id": str(record["unit_id"]),
                        "crossfit_fold": int(fold),
                        "crossfit_assignment_seed": int(
                            record["crossfit_assignment_seed"]
                        ),
                        "observed_selection_score": observed,
                        "n_training_null": len(null),
                        "crossfit_mc_p_upper": upper,
                        "crossfit_mc_p_lower": lower,
                        "crossfit_mc_p_two_sided": _two_sided_pvalue(upper, lower),
                        "calibration_role": "primary_out_of_sample_crossfit",
                    }
                )
    crossfit = pd.DataFrame(rows)
    calibration_rows: list[dict[str, Any]] = []
    for labels, group in crossfit.groupby(group_columns, sort=True, dropna=False):
        label = dict(zip(group_columns, labels, strict=True))
        pvalues = group["crossfit_mc_p_upper"].to_numpy(dtype=float)
        pvalues = pvalues[np.isfinite(pvalues)]
        fold_sizes = group.groupby("crossfit_fold", sort=True).size().to_dict()
        calibration_rows.append(
            {
                **label,
                "n_neutral": len(pvalues),
                "n_folds": int(n_folds),
                "fold_sizes": json.dumps(fold_sizes, sort_keys=True),
                "minimum_training_null": int(group["n_training_null"].min()),
                "maximum_training_null": int(group["n_training_null"].max()),
                "mean_crossfit_p": float(np.mean(pvalues)),
                "fraction_crossfit_p_le_0_05": float(np.mean(pvalues <= 0.05)),
                "ks_uniform_p": float(kstest(pvalues, "uniform").pvalue),
                "calibration_note": (
                    "primary held-out diagnostic; deterministic folds are shared "
                    "across estimands and finite-null p-values are discrete"
                ),
            }
        )
    return crossfit, pd.DataFrame(calibration_rows)


def cell_roc_comparisons(
    scores: pd.DataFrame,
    *,
    bootstrap_draws: int = DEFAULT_BOOTSTRAP_DRAWS,
    permutation_draws: int = DEFAULT_LABEL_PERMUTATION_DRAWS,
    base_seed: int = BASE_SEED,
) -> pd.DataFrame:
    """Summarize selected-versus-neutral discrimination in each matched cell."""

    if int(permutation_draws) < 1:
        raise ValueError("permutation_draws must be positive")

    group_columns = [
        *CELL_COLUMNS,
        "statistic_scope",
        "metric",
        "statistic",
        "threshold_years",
    ]
    rows: list[dict[str, Any]] = []
    permutation_cache: dict[tuple[Any, ...], tuple[int, np.ndarray]] = {}
    for labels, group in scores.groupby(group_columns, sort=True, dropna=False):
        label = dict(zip(group_columns, labels, strict=True))
        neutral_frame = group[group["simulation_class"] == "neutral"].sort_values(
            "unit_id"
        )
        neutral = neutral_frame["selection_score"].to_numpy(dtype=float)
        selected_group = group[group["simulation_class"] == "selected"]
        for coefficient, coefficient_group in selected_group.groupby(
            "selection_coefficient", sort=True
        ):
            coefficient_group = coefficient_group.sort_values("unit_id")
            selected = coefficient_group["selection_score"].to_numpy(dtype=float)
            auc = _roc_auc(selected, neutral)
            seed_label = ":".join(map(str, (*labels, coefficient)))
            seed = _stable_seed(base_seed, "roc-bootstrap:" + seed_label)
            low, high = _bootstrap_auc_interval(
                selected, neutral, draws=bootstrap_draws, seed=seed
            )
            difference_seed = _stable_seed(
                base_seed, "difference-bootstrap:" + seed_label
            )
            difference_low, difference_high = _bootstrap_difference_interval(
                selected,
                neutral,
                draws=bootstrap_draws,
                seed=difference_seed,
            )
            test = mannwhitneyu(selected, neutral, alternative="greater", method="auto")
            permutation_key = (
                str(labels[0]),
                float(labels[1]),
                float(coefficient),
                tuple(coefficient_group["unit_id"].astype(str)),
                tuple(neutral_frame["unit_id"].astype(str)),
            )
            cached_permutation = permutation_cache.get(permutation_key)
            if cached_permutation is None:
                permutation_seed = _stable_seed(
                    base_seed,
                    "unit-label-permutation:" + ":".join(map(str, permutation_key[:3])),
                )
                rng = np.random.default_rng(permutation_seed)
                random_keys = rng.random(
                    (int(permutation_draws), len(selected) + len(neutral))
                )
                selected_indices = np.argpartition(
                    random_keys, len(selected) - 1, axis=1
                )[:, : len(selected)].copy()
                cached_permutation = (permutation_seed, selected_indices)
                permutation_cache[permutation_key] = cached_permutation
            permutation_seed, selected_indices = cached_permutation
            (
                mean_permutation_p,
                auc_permutation_p,
                mean_permutation_exceedances,
                auc_permutation_exceedances,
            ) = _unit_label_permutation_pvalues(
                selected,
                neutral,
                selected_indices=selected_indices,
            )
            rows.append(
                {
                    **label,
                    "selection_coefficient": float(coefficient),
                    "n_selected": len(selected),
                    "n_neutral": len(neutral),
                    "mean_selected_score": float(np.mean(selected)),
                    "mean_neutral_score": float(np.mean(neutral)),
                    "mean_selected_minus_neutral": float(
                        np.mean(selected) - np.mean(neutral)
                    ),
                    "mean_difference_ci95_low": difference_low,
                    "mean_difference_ci95_high": difference_high,
                    "roc_auc": auc,
                    "roc_auc_ci95_low": low,
                    "roc_auc_ci95_high": high,
                    "bootstrap_draws": int(bootstrap_draws),
                    "bootstrap_seed": seed,
                    "difference_bootstrap_seed": difference_seed,
                    "mann_whitney_u": float(test.statistic),
                    "mann_whitney_p_upper": float(test.pvalue),
                    "mean_unit_label_permutation_p_upper": mean_permutation_p,
                    "auc_unit_label_permutation_p_upper": auc_permutation_p,
                    "unit_label_permutation_draws": int(permutation_draws),
                    "unit_label_permutation_seed": int(permutation_seed),
                    "mean_unit_label_permutation_exceedances": int(
                        mean_permutation_exceedances
                    ),
                    "auc_unit_label_permutation_exceedances": int(
                        auc_permutation_exceedances
                    ),
                    "unit_label_permutation_method": (
                        "seeded_monte_carlo_preserve_group_sizes_plus_one"
                    ),
                    "unit_label_permutation_family": (
                        "same unit-label draws reused across estimands within cell and s"
                    ),
                    "neutral_bank_scope": (
                        "shared_by_demography_and_af_across_selection_coefficients"
                    ),
                }
            )
    result = pd.DataFrame(rows)
    family = [
        *CELL_COLUMNS,
        "selection_coefficient",
        "statistic_scope",
        "metric",
        "statistic",
    ]
    result["mann_whitney_bh_q_upper"] = result.groupby(family, sort=False)[
        "mann_whitney_p_upper"
    ].transform(lambda values: bh_fdr(values.to_numpy(dtype=float)))
    global_family = ["statistic_scope", "metric", "statistic"]
    result["mann_whitney_global_bh_q_upper"] = result.groupby(
        global_family, sort=False
    )["mann_whitney_p_upper"].transform(
        lambda values: bh_fdr(values.to_numpy(dtype=float))
    )
    for prefix in ("mean", "auc"):
        pvalue_column = f"{prefix}_unit_label_permutation_p_upper"
        result[f"{prefix}_unit_label_permutation_bh_q_upper"] = result.groupby(
            family, sort=False
        )[pvalue_column].transform(lambda values: bh_fdr(values.to_numpy(dtype=float)))
        result[f"{prefix}_unit_label_permutation_global_bh_q_upper"] = result.groupby(
            global_family, sort=False
        )[pvalue_column].transform(lambda values: bh_fdr(values.to_numpy(dtype=float)))
    return result


def within_selected_summary(
    scores: pd.DataFrame,
    *,
    bootstrap_draws: int = DEFAULT_BOOTSTRAP_DRAWS,
    base_seed: int = BASE_SEED,
) -> pd.DataFrame:
    """Summarize paired hom-alt versus hom-ref effects within selected panels."""

    selected = scores[
        (scores["simulation_class"] == "selected")
        & scores["statistic"]
        .astype(str)
        .str.contains("hom_alt_minus_hom_ref|hom_ref_minus_hom_alt")
    ]
    group_columns = [
        *CELL_COLUMNS,
        "selection_coefficient",
        "statistic_scope",
        "metric",
        "statistic",
        "threshold_years",
    ]
    rows: list[dict[str, Any]] = []
    for labels, group in selected.groupby(group_columns, sort=True, dropna=False):
        label = dict(zip(group_columns, labels, strict=True))
        values = group["selection_score"].to_numpy(dtype=float)
        nonzero = values[values != 0]
        positives = int(np.count_nonzero(nonzero > 0))
        sign_p = (
            float(binomtest(positives, len(nonzero), 0.5, alternative="greater").pvalue)
            if len(nonzero)
            else np.nan
        )
        seed = _stable_seed(base_seed, "within-bootstrap:" + ":".join(map(str, labels)))
        low, high = _bootstrap_mean_interval(values, draws=bootstrap_draws, seed=seed)
        rows.append(
            {
                **label,
                "n_selected": len(values),
                "mean_within_selection_score": float(np.mean(values)),
                "median_within_selection_score": float(np.median(values)),
                "sd_within_selection_score": (
                    float(np.std(values, ddof=1)) if len(values) > 1 else np.nan
                ),
                "mean_ci95_low": low,
                "mean_ci95_high": high,
                "n_nonzero": len(nonzero),
                "n_positive": positives,
                "sign_test_p_upper": sign_p,
                "bootstrap_draws": int(bootstrap_draws),
                "bootstrap_seed": seed,
            }
        )
    result = pd.DataFrame(rows)
    family = [
        *CELL_COLUMNS,
        "selection_coefficient",
        "statistic_scope",
        "metric",
        "statistic",
    ]
    if not result.empty:
        result["sign_test_bh_q_upper"] = result.groupby(family, sort=False)[
            "sign_test_p_upper"
        ].transform(lambda values: bh_fdr(values.to_numpy(dtype=float)))
    return result


def within_selected_unit_pair_tests(
    pair_summaries: pd.DataFrame,
    *,
    permutation_draws: int = 9_999,
    base_seed: int = BASE_SEED,
) -> pd.DataFrame:
    """Test alt/alt versus ref/ref pair posteriors within each selected unit.

    These tests diagnose genotype association inside an accepted realization;
    they are not demographic selection p-values.  Labels are permuted only
    among homozygous pairs, preserving the observed hom-alt/ref class counts.
    """

    required = {
        *UNIT_COLUMNS,
        "pair_index",
        "genotype_class",
        "threshold_years",
        "region_mean_p_tmrca_lt_threshold",
        "focal_p_tmrca_lt_threshold",
        "region_mean_tmrca_generations",
        "focal_tmrca_generations",
    }
    missing = sorted(required.difference(pair_summaries.columns))
    if missing:
        raise ValueError(f"pair summaries are missing columns: {', '.join(missing)}")
    if int(permutation_draws) < 1:
        raise ValueError("permutation_draws must be positive")
    frame = pair_summaries.copy()
    frame = frame[frame["simulation_class"].astype(str) == "selected"]
    if frame.empty:
        return pd.DataFrame()
    for column in (
        "selection_coefficient",
        "target_allele_frequency",
        "replicate_index",
        "pair_index",
        "threshold_years",
        "region_mean_p_tmrca_lt_threshold",
        "focal_p_tmrca_lt_threshold",
        "region_mean_tmrca_generations",
        "focal_tmrca_generations",
    ):
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
        if not np.isfinite(frame[column].to_numpy(dtype=float)).all():
            raise ValueError(f"pair-summary column is not finite: {column}")
    if (frame["selection_coefficient"] <= 0).any():
        raise ValueError(
            "selected pair summaries require positive selection coefficients"
        )
    frame["genotype_class"] = frame["genotype_class"].astype(str)
    if set(frame["genotype_class"]) - {"hom_ref", "heterozygous", "hom_alt"}:
        raise ValueError("pair summaries contain unexpected genotype classes")
    key = ["unit_id", "pair_index", "threshold_years"]
    if frame.duplicated(key, keep=False).any():
        raise ValueError("pair summaries contain duplicate unit/pair/threshold rows")
    if (
        (frame["region_mean_p_tmrca_lt_threshold"] < 0)
        | (frame["region_mean_p_tmrca_lt_threshold"] > 1)
        | (frame["focal_p_tmrca_lt_threshold"] < 0)
        | (frame["focal_p_tmrca_lt_threshold"] > 1)
    ).any():
        raise ValueError("pair P(TMRCA < x) values must lie in [0, 1]")
    for (unit_id, pair_index), curve in frame.groupby(
        ["unit_id", "pair_index"], sort=False
    ):
        ordered = curve.sort_values("threshold_years")
        for column in (
            "region_mean_p_tmrca_lt_threshold",
            "focal_p_tmrca_lt_threshold",
        ):
            if np.any(np.diff(ordered[column].to_numpy(dtype=float)) < -1e-10):
                raise ValueError(
                    f"pair P(TMRCA < x) is not nondecreasing for "
                    f"{unit_id} pair {pair_index}"
                )
    rows: list[dict[str, Any]] = []
    metadata_columns = list(UNIT_COLUMNS)

    def add_test(
        unit: pd.DataFrame,
        *,
        scope: str,
        metric: str,
        statistic: str,
        threshold_years: float,
        alt_scores: np.ndarray,
        ref_scores: np.ndarray,
    ) -> None:
        metadata = unit.iloc[0][metadata_columns].to_dict()
        seed_label = ":".join(
            map(
                str,
                (
                    metadata["unit_id"],
                    scope,
                    metric,
                    statistic,
                    threshold_years,
                ),
            )
        )
        seed = _stable_seed(base_seed, "within-label-permutation:" + seed_label)
        pvalue, method, effective_draws, exceedances = (
            _genotype_label_permutation_pvalue(
                alt_scores,
                ref_scores,
                seed=seed,
                draws=permutation_draws,
            )
        )
        rows.append(
            {
                **metadata,
                "statistic_scope": scope,
                "metric": metric,
                "statistic": statistic,
                "threshold_years": float(threshold_years),
                "n_hom_alt_pairs": len(alt_scores),
                "n_hom_ref_pairs": len(ref_scores),
                "mean_hom_alt_selection_score": float(np.mean(alt_scores)),
                "mean_hom_ref_selection_score": float(np.mean(ref_scores)),
                "hom_alt_minus_hom_ref_selection_score": float(
                    np.mean(alt_scores) - np.mean(ref_scores)
                ),
                "within_pair_roc_auc": _roc_auc(alt_scores, ref_scores),
                "genotype_label_permutation_p_upper": pvalue,
                "permutation_method": method,
                "permutation_draws_or_assignments": effective_draws,
                "permutation_exceedances": exceedances,
                "permutation_seed": seed,
                "interpretation_scope": (
                    "within_realization_genotype_association_not_demographic_selection_p"
                ),
            }
        )

    for _, unit in frame.groupby("unit_id", sort=True):
        unit_metadata = unit[list(UNIT_COLUMNS)].drop_duplicates()
        if len(unit_metadata) != 1:
            raise ValueError("unit metadata changes within pair summaries")
        counts = (
            unit[["pair_index", "genotype_class"]]
            .drop_duplicates()
            .groupby("genotype_class")
            .size()
        )
        if int(counts.get("hom_alt", 0)) < 1 or int(counts.get("hom_ref", 0)) < 1:
            raise ValueError("within-unit tests require hom-alt and hom-ref pairs")
        for scope, probability_column in (
            ("region_mean", "region_mean_p_tmrca_lt_threshold"),
            ("focal_nearest", "focal_p_tmrca_lt_threshold"),
        ):
            for threshold, threshold_group in unit.groupby(
                "threshold_years", sort=True
            ):
                alt = threshold_group.loc[
                    threshold_group["genotype_class"] == "hom_alt",
                    probability_column,
                ].to_numpy(dtype=float)
                ref = threshold_group.loc[
                    threshold_group["genotype_class"] == "hom_ref",
                    probability_column,
                ].to_numpy(dtype=float)
                add_test(
                    unit,
                    scope=scope,
                    metric="p_tmrca_lt_threshold",
                    statistic=PRIMARY_STATISTIC,
                    threshold_years=float(threshold),
                    alt_scores=alt,
                    ref_scores=ref,
                )

            homozygous = unit[unit["genotype_class"].isin(["hom_ref", "hom_alt"])]
            probability = homozygous.pivot(
                index=["pair_index", "genotype_class"],
                columns="threshold_years",
                values=probability_column,
            ).sort_index(axis=1)
            x = np.concatenate(([0.0], probability.columns.to_numpy(dtype=float)))
            areas = (
                trapezoid(
                    np.column_stack(
                        (np.zeros(len(probability)), probability.to_numpy(dtype=float))
                    ),
                    x=x,
                    axis=1,
                )
                / x[-1]
            )
            area_classes = probability.index.get_level_values("genotype_class")
            add_test(
                unit,
                scope=scope,
                metric="normalized_cdf_area",
                statistic="cdf_area_hom_alt_minus_hom_ref",
                threshold_years=0.0,
                alt_scores=areas[area_classes == "hom_alt"],
                ref_scores=areas[area_classes == "hom_ref"],
            )

            tmrca_column = (
                "region_mean_tmrca_generations"
                if scope == "region_mean"
                else "focal_tmrca_generations"
            )
            first = homozygous.sort_values("threshold_years").drop_duplicates(
                ["pair_index"], keep="first"
            )
            # Negate TMRCA so every reported selection score has the same
            # direction: larger means younger/more sweep-like.
            alt_tmrca_score = -first.loc[
                first["genotype_class"] == "hom_alt", tmrca_column
            ].to_numpy(dtype=float)
            ref_tmrca_score = -first.loc[
                first["genotype_class"] == "hom_ref", tmrca_column
            ].to_numpy(dtype=float)
            add_test(
                unit,
                scope=scope,
                metric="mean_tmrca_generations",
                statistic="negative_tmrca_hom_alt_minus_hom_ref",
                threshold_years=0.0,
                alt_scores=alt_tmrca_score,
                ref_scores=ref_tmrca_score,
            )
    result = pd.DataFrame(rows)
    family = [
        "unit_id",
        "statistic_scope",
        "metric",
        "statistic",
    ]
    result["genotype_label_permutation_bh_q_upper"] = result.groupby(
        family, sort=False
    )["genotype_label_permutation_p_upper"].transform(
        lambda values: bh_fdr(values.to_numpy(dtype=float))
    )
    return result


def empirical_input_template(class_summaries: pd.DataFrame) -> pd.DataFrame:
    """Build a deliberately incomplete, fail-closed empirical input template."""

    rows: list[dict[str, Any]] = []
    thresholds = sorted(class_summaries["threshold_years"].astype(float).unique())
    for labels, cell_frame in class_summaries.groupby(
        list(CELL_COLUMNS), sort=True, dropna=False
    ):
        overall_pairs = cell_frame.loc[
            cell_frame["genotype_class"] == "overall", "n_pairs"
        ].astype(int)
        if overall_pairs.nunique() != 1:
            raise ValueError("empirical template requires one panel size per cell")
        sources = cell_frame["source"].astype(str).unique()
        if len(sources) != 1:
            raise ValueError("empirical template requires one analysis source per cell")
        cell = dict(zip(CELL_COLUMNS, labels, strict=True))
        for scope in SCOPES:
            for statistic in (
                PRIMARY_STATISTIC,
                "p_hom_alt",
                "p_hom_ref",
                "p_heterozygous",
                "p_overall",
            ):
                for threshold in thresholds:
                    rows.append(
                        {
                            **cell,
                            "statistic_scope": scope,
                            "metric": "p_tmrca_lt_threshold",
                            "statistic": statistic,
                            "threshold_years": threshold,
                            "threshold_generations": threshold / 25.0,
                            "observed_p_overall": np.nan,
                            "observed_p_hom_ref": np.nan,
                            "observed_p_heterozygous": np.nan,
                            "observed_p_hom_alt": np.nan,
                            "observed_value": np.nan,
                            "empirical_id": "replace_with_empirical_locus_id",
                            "variant_id": "replace_with_variant_id",
                            "source_label": "replace_with_empirical_dataset_id",
                            "source_data_sha256": "replace_with_sha256",
                            "sample_manifest_sha256": "replace_with_sha256",
                            "pair_table_sha256": "replace_with_sha256",
                            "demography_assignment": (
                                "replace_with_population_and_demography_assignment"
                            ),
                            "population": "replace_with_population",
                            "observed_allele_frequency": np.nan,
                            "af_match_tolerance": 0.025,
                            "n_diploid_samples": int(overall_pairs.iloc[0]),
                            "n_overall_pairs": int(overall_pairs.iloc[0]),
                            "n_hom_ref_pairs": np.nan,
                            "n_heterozygous_pairs": np.nan,
                            "n_hom_alt_pairs": np.nan,
                            "alt_allele_count": np.nan,
                            "pairing_scheme": "within_individual_diploid",
                            "decoder_name": str(sources[0]),
                            "decoder_version": "replace_with_decoder_version",
                            "reference_decoder_binary_sha256": "replace_with_sha256",
                            "empirical_decoder_binary_sha256": "replace_with_sha256",
                            "reference_decoder_parameters_sha256": (
                                "replace_with_sha256"
                            ),
                            "empirical_decoder_parameters_sha256": (
                                "replace_with_sha256"
                            ),
                            "reference_ne_history_sha256": "replace_with_sha256",
                            "empirical_ne_history_sha256": "replace_with_sha256",
                            "generation_time_years": 25.0,
                            "mutation_rate_per_bp_per_generation": 1.25e-8,
                            "sequence_length_bp": 10_000_000,
                            "genome_build": "replace_with_genome_build",
                            "contig": "replace_with_contig",
                            "position_1based": np.nan,
                            "ref_allele": "replace_with_ref",
                            "alt_allele": "replace_with_alt",
                            "region_start_0based": np.nan,
                            "region_end_0based": np.nan,
                            "focal_position_0based": np.nan,
                            "mask_manifest_sha256": "replace_with_sha256",
                            "qc_status": "replace_with_pass_after_qc",
                        }
                    )
    return pd.DataFrame(rows, columns=EMPIRICAL_COLUMNS)


def _validate_empirical_provenance(
    data: pd.DataFrame, scores: pd.DataFrame
) -> pd.DataFrame:
    result = data.copy()
    numeric_columns = (
        "target_allele_frequency",
        "threshold_years",
        "threshold_generations",
        "observed_p_overall",
        "observed_p_hom_ref",
        "observed_p_heterozygous",
        "observed_p_hom_alt",
        "observed_value",
        "observed_allele_frequency",
        "af_match_tolerance",
        "n_diploid_samples",
        "n_overall_pairs",
        "n_hom_ref_pairs",
        "n_heterozygous_pairs",
        "n_hom_alt_pairs",
        "alt_allele_count",
        "generation_time_years",
        "mutation_rate_per_bp_per_generation",
        "sequence_length_bp",
        "position_1based",
        "region_start_0based",
        "region_end_0based",
        "focal_position_0based",
    )
    for column in numeric_columns:
        result[column] = pd.to_numeric(result[column], errors="coerce")
    if not np.isfinite(result[list(numeric_columns)].to_numpy(dtype=float)).all():
        raise ValueError("empirical TSV contains missing or nonnumeric required values")
    text_columns = [
        column
        for column in EMPIRICAL_PROVENANCE_COLUMNS
        if column not in numeric_columns
    ]
    for column in text_columns:
        values = result[column].astype(str).str.strip()
        if (
            values.eq("").any()
            or values.str.lower().isin({"nan", "none", "<na>"}).any()
            or values.str.startswith("replace_with_").any()
        ):
            raise ValueError(f"empirical provenance field is unset: {column}")
        result[column] = values
    sha_columns = (
        "source_data_sha256",
        "sample_manifest_sha256",
        "pair_table_sha256",
        "reference_decoder_binary_sha256",
        "empirical_decoder_binary_sha256",
        "reference_decoder_parameters_sha256",
        "empirical_decoder_parameters_sha256",
        "reference_ne_history_sha256",
        "empirical_ne_history_sha256",
        "mask_manifest_sha256",
    )
    sha_pattern = re.compile(r"^[0-9a-f]{64}$")
    for column in sha_columns:
        result[column] = result[column].str.lower()
        if (
            not result[column]
            .map(lambda value: bool(sha_pattern.fullmatch(value)))
            .all()
        ):
            raise ValueError(f"empirical provenance SHA-256 is malformed: {column}")
    for reference_column, empirical_column in (
        ("reference_decoder_binary_sha256", "empirical_decoder_binary_sha256"),
        (
            "reference_decoder_parameters_sha256",
            "empirical_decoder_parameters_sha256",
        ),
        ("reference_ne_history_sha256", "empirical_ne_history_sha256"),
    ):
        if not (
            result[reference_column].to_numpy() == result[empirical_column].to_numpy()
        ).all():
            raise ValueError(
                "empirical decoder/Ne setting does not match its neutral reference: "
                f"{empirical_column}"
            )
    integer_columns = (
        "n_diploid_samples",
        "n_overall_pairs",
        "n_hom_ref_pairs",
        "n_heterozygous_pairs",
        "n_hom_alt_pairs",
        "alt_allele_count",
        "sequence_length_bp",
        "position_1based",
        "region_start_0based",
        "region_end_0based",
        "focal_position_0based",
    )
    for column in integer_columns:
        if (result[column] != np.floor(result[column])).any():
            raise ValueError(
                f"empirical provenance count/coordinate is not integer: {column}"
            )
        result[column] = result[column].astype(np.int64)
    required_thresholds = {float(value) for value in TMRCA_THRESHOLDS_YEARS}
    if not set(result["threshold_years"].astype(float)).issubset(required_thresholds):
        raise ValueError("empirical threshold is not in the prespecified seven")
    if not np.allclose(
        result["threshold_generations"],
        result["threshold_years"] / result["generation_time_years"],
        rtol=0.0,
        atol=1e-12,
    ):
        raise ValueError(
            "empirical threshold-year/generation conversion is inconsistent"
        )
    probability_columns = (
        "observed_p_overall",
        "observed_p_hom_ref",
        "observed_p_heterozygous",
        "observed_p_hom_alt",
    )
    if any(
        ((result[column] < 0) | (result[column] > 1)).any()
        for column in probability_columns
    ):
        raise ValueError("empirical genotype-class probabilities must lie in [0, 1]")
    expected_statistic_values = {
        PRIMARY_STATISTIC: result["observed_p_hom_alt"] - result["observed_p_hom_ref"],
        "p_hom_alt": result["observed_p_hom_alt"],
        "p_hom_ref": result["observed_p_hom_ref"],
        "p_heterozygous": result["observed_p_heterozygous"],
        "p_overall": result["observed_p_overall"],
    }
    if set(result["statistic"]) - set(expected_statistic_values):
        raise ValueError("empirical statistic is not supported by the aggregate schema")
    expected_observed = np.empty(len(result), dtype=float)
    for statistic, values in expected_statistic_values.items():
        mask = result["statistic"] == statistic
        expected_observed[mask] = values[mask].to_numpy(dtype=float)
    if not np.allclose(
        result["observed_value"], expected_observed, rtol=0.0, atol=1e-12
    ):
        raise ValueError("empirical observed statistic disagrees with genotype classes")
    curve_key = ["empirical_id", *CELL_COLUMNS, "statistic_scope"]
    for labels, group in result.groupby(curve_key, sort=False, dropna=False):
        for threshold, threshold_group in group.groupby("threshold_years", sort=True):
            if any(
                threshold_group[column].nunique(dropna=False) != 1
                for column in probability_columns
            ):
                raise ValueError(
                    "empirical genotype probability changes across statistic rows at "
                    f"{':'.join(map(str, labels))}:{threshold:g}"
                )
        curve = group.drop_duplicates("threshold_years").sort_values("threshold_years")
        if any(
            np.any(np.diff(curve[column].to_numpy(dtype=float)) < -1e-12)
            for column in probability_columns
        ):
            raise ValueError(
                "empirical genotype probability is not nondecreasing for "
                + ":".join(map(str, labels))
            )
    pair_counts = result[["n_hom_ref_pairs", "n_heterozygous_pairs", "n_hom_alt_pairs"]]
    if (
        (result["n_hom_ref_pairs"] < 2)
        | (result["n_heterozygous_pairs"] < 1)
        | (result["n_hom_alt_pairs"] < 2)
    ).any():
        raise ValueError(
            "empirical genotype pair counts fail the prespecified QC minima"
        )
    if (result["n_diploid_samples"] <= 0).any():
        raise ValueError("empirical diploid sample count must be positive")
    if not (
        pair_counts.sum(axis=1).to_numpy() == result["n_diploid_samples"].to_numpy()
    ).all():
        raise ValueError(
            "empirical genotype pair counts do not sum to the sample count"
        )
    if not (
        result["n_overall_pairs"].to_numpy() == result["n_diploid_samples"].to_numpy()
    ).all():
        raise ValueError("empirical overall pair count does not equal the sample count")
    derived_alt_count = 2 * result["n_hom_alt_pairs"] + result["n_heterozygous_pairs"]
    if not (
        derived_alt_count.to_numpy() == result["alt_allele_count"].to_numpy()
    ).all():
        raise ValueError("empirical alt count disagrees with genotype pair counts")
    derived_af = derived_alt_count / (2 * result["n_diploid_samples"])
    if not np.allclose(
        derived_af,
        result["observed_allele_frequency"],
        rtol=0.0,
        atol=1e-12,
    ):
        raise ValueError("empirical allele frequency disagrees with allele counts")
    if (
        (result["target_allele_frequency"] <= 0)
        | (result["target_allele_frequency"] >= 1)
        | (result["region_start_0based"] < 0)
        | (result["region_end_0based"] <= result["region_start_0based"])
        | (result["observed_allele_frequency"] <= 0)
        | (result["observed_allele_frequency"] >= 1)
        | (result["af_match_tolerance"] < 0)
        | (result["af_match_tolerance"] > 0.05)
    ).any():
        raise ValueError(
            "empirical allele frequency, AF tolerance, or region is invalid"
        )
    af_difference = (
        result["observed_allele_frequency"] - result["target_allele_frequency"]
    ).abs()
    if (af_difference > result["af_match_tolerance"] + 1e-12).any():
        raise ValueError("empirical allele frequency is outside its matched AF stratum")
    if set(result["pairing_scheme"]) != {"within_individual_diploid"}:
        raise ValueError("empirical pairing_scheme must be within_individual_diploid")
    if set(result["qc_status"].str.lower()) != {"pass"}:
        raise ValueError("empirical qc_status must be pass")
    if not np.allclose(result["generation_time_years"], 25.0, rtol=0.0, atol=1e-12):
        raise ValueError("empirical generation time must match 25 years")
    if not np.allclose(
        result["mutation_rate_per_bp_per_generation"],
        1.25e-8,
        rtol=0.0,
        atol=1e-20,
    ):
        raise ValueError("empirical mutation rate does not match the neutral reference")
    if (result["sequence_length_bp"] != 10_000_000).any():
        raise ValueError("empirical sequence length must match 10 Mb")
    if not (
        result["region_end_0based"] - result["region_start_0based"]
        == result["sequence_length_bp"]
    ).all():
        raise ValueError("empirical region coordinates do not span 10 Mb")
    if (
        (result["focal_position_0based"] < result["region_start_0based"])
        | (result["focal_position_0based"] >= result["region_end_0based"])
    ).any():
        raise ValueError("empirical focal coordinate lies outside the decoded region")
    if not (
        result["position_1based"].to_numpy()
        == result["focal_position_0based"].to_numpy() + 1
    ).all():
        raise ValueError("empirical 1-based and 0-based focal coordinates disagree")
    valid_alleles = {"A", "C", "G", "T"}
    result["ref_allele"] = result["ref_allele"].str.upper()
    result["alt_allele"] = result["alt_allele"].str.upper()
    if (
        set(result["ref_allele"]) - valid_alleles
        or set(result["alt_allele"]) - valid_alleles
        or (result["ref_allele"] == result["alt_allele"]).any()
    ):
        raise ValueError("empirical ref/alt alleles must be distinct SNV bases")
    neutral = scores[scores["simulation_class"] == "neutral"]
    reference = neutral[
        [*CELL_COLUMNS, "panel_n_diploid_samples", "analysis_source"]
    ].drop_duplicates()
    if reference.duplicated(list(CELL_COLUMNS), keep=False).any():
        raise ValueError("neutral reference has inconsistent panel or decoder settings")
    checked = result.merge(
        reference,
        on=list(CELL_COLUMNS),
        how="left",
        validate="many_to_one",
    )
    if checked[["panel_n_diploid_samples", "analysis_source"]].isna().any().any():
        raise ValueError(
            "empirical demography/AF cell has no neutral reference settings"
        )
    if not (
        checked["n_diploid_samples"].to_numpy()
        == checked["panel_n_diploid_samples"].astype(int).to_numpy()
    ).all():
        raise ValueError("empirical sample size does not match the neutral reference")
    if not (
        checked["decoder_name"].astype(str).to_numpy()
        == checked["analysis_source"].astype(str).to_numpy()
    ).all():
        raise ValueError("empirical decoder name does not match the neutral reference")
    provenance_key = ["source_label", *CELL_COLUMNS]
    for labels, group in checked.groupby(provenance_key, sort=False, dropna=False):
        if len(group[list(EMPIRICAL_PROVENANCE_COLUMNS)].drop_duplicates()) != 1:
            raise ValueError(
                "empirical provenance changes across statistic rows for "
                + ":".join(map(str, labels))
            )
    checked["af_absolute_difference"] = (
        checked["observed_allele_frequency"] - checked["target_allele_frequency"]
    ).abs()
    checked["genotype_derived_alt_count"] = (
        2 * checked["n_hom_alt_pairs"] + checked["n_heterozygous_pairs"]
    ).astype(int)
    checked["provenance_qc_status"] = "pass"
    return checked


def empirical_null_pvalues(
    empirical: pd.DataFrame | None, scores: pd.DataFrame
) -> pd.DataFrame:
    """Match optional empirical selection-directed scores to the neutral null."""

    output_columns = [
        *EMPIRICAL_COLUMNS,
        "n_null",
        "null_mean",
        "null_sd",
        "mc_p_upper",
        "mc_p_lower",
        "mc_p_two_sided",
        "bh_q_upper",
        "bh_q_two_sided",
        "panel_n_diploid_samples",
        "analysis_source",
        "af_absolute_difference",
        "genotype_derived_alt_count",
        "provenance_qc_status",
    ]
    if empirical is None or empirical.empty:
        return pd.DataFrame(columns=output_columns)
    missing = sorted(set(EMPIRICAL_COLUMNS).difference(empirical.columns))
    if missing:
        raise ValueError(f"empirical TSV is missing columns: {', '.join(missing)}")
    data = _validate_empirical_provenance(
        empirical.loc[:, EMPIRICAL_COLUMNS].copy(), scores
    )
    key = [
        *CELL_COLUMNS,
        "statistic_scope",
        "metric",
        "statistic",
        "threshold_years",
    ]
    if data.duplicated(key, keep=False).any():
        raise ValueError("empirical TSV contains duplicate statistic rows")
    neutral = scores[scores["simulation_class"] == "neutral"]
    rows: list[dict[str, Any]] = []
    for record in data.to_dict(orient="records"):
        mask = np.ones(len(neutral), dtype=bool)
        for column in key:
            if column in {"target_allele_frequency", "threshold_years"}:
                mask &= np.isclose(
                    neutral[column].to_numpy(dtype=float),
                    float(record[column]),
                    rtol=0.0,
                    atol=1e-9,
                )
            else:
                mask &= neutral[column].astype(str).to_numpy() == str(record[column])
        null = neutral.loc[mask, "selection_score"].to_numpy(dtype=float)
        if not len(null):
            raise ValueError(
                "empirical statistic has no exactly matched neutral null: "
                + ":".join(str(record[column]) for column in key)
            )
        observed = float(record["observed_value"])
        upper = monte_carlo_pvalue(observed, null)
        lower = _lower_tail_pvalue(observed, null)
        rows.append(
            {
                **record,
                "n_null": len(null),
                "null_mean": float(np.mean(null)),
                "null_sd": float(np.std(null, ddof=1)) if len(null) > 1 else np.nan,
                "mc_p_upper": upper,
                "mc_p_lower": lower,
                "mc_p_two_sided": _two_sided_pvalue(upper, lower),
            }
        )
    result = pd.DataFrame(rows)
    family = [
        *CELL_COLUMNS,
        "statistic_scope",
        "metric",
        "statistic",
    ]
    result["bh_q_upper"] = result.groupby(family, sort=False)["mc_p_upper"].transform(
        lambda values: bh_fdr(values.to_numpy(dtype=float))
    )
    result["bh_q_two_sided"] = result.groupby(family, sort=False)[
        "mc_p_two_sided"
    ].transform(lambda values: bh_fdr(values.to_numpy(dtype=float)))
    return result.loc[:, output_columns]


def _cell_label(demography: str, af: float, coefficient: float | None = None) -> str:
    short = "EAS median" if "eas_phlash" in demography else "Han introgression"
    suffix = "" if coefficient is None else f", s={coefficient:g}"
    return f"{short}, AF {100 * af:.0f}%{suffix}"


def _source_display_name(frame: pd.DataFrame) -> str:
    sources = set(frame["source"].astype(str)) if "source" in frame else set()
    if sources == {"gamma_smc"}:
        return "Gamma-SMC"
    if sources == {"tree_truth"}:
        return "Tree truth"
    if len(sources) == 1:
        return next(iter(sources)).replace("_", " ").title()
    raise ValueError("plot input must contain exactly one analysis source")


def plot_probability_curves(
    class_summaries: pd.DataFrame,
    figure_dir: Path,
    *,
    source_label: str | None = None,
) -> list[Path]:
    """Plot selected/ref genotype profiles against matched neutral profiles."""

    outputs: list[Path] = []
    display_source = source_label or _source_display_name(class_summaries)
    colors = {
        "overall": "#000000",
        "hom_ref": "#595959",
        "heterozygous": "#E69F00",
        "hom_alt": "#0072B2",
    }
    for demography, demographic in class_summaries.groupby("demography_id", sort=True):
        coefficients = sorted(
            demographic.loc[
                demographic["simulation_class"] == "selected",
                "selection_coefficient",
            ]
            .astype(float)
            .unique()
        )
        for coefficient in coefficients:
            comparison = demographic[
                (demographic["simulation_class"] == "neutral")
                | (
                    (demographic["simulation_class"] == "selected")
                    & np.isclose(demographic["selection_coefficient"], coefficient)
                )
            ]
            afs = sorted(comparison["target_allele_frequency"].astype(float).unique())
            figure, axes = plt.subplots(1, len(afs), figsize=(11, 8.5), sharey=True)
            axes = np.atleast_1d(axes)
            for axis, af in zip(axes, afs, strict=True):
                cell = comparison[np.isclose(comparison["target_allele_frequency"], af)]
                for simulation_class, linestyle in (
                    ("selected", "-"),
                    ("neutral", "--"),
                ):
                    for genotype_class in GENOTYPE_CLASSES:
                        group = cell[
                            (cell["simulation_class"] == simulation_class)
                            & (cell["genotype_class"] == genotype_class)
                        ]
                        curve = group.groupby("threshold_years")[
                            "focal_mean_p_tmrca_lt_threshold"
                        ].mean()
                        axis.plot(
                            curve.index.to_numpy(dtype=float) / 1_000,
                            curve.to_numpy(dtype=float),
                            color=colors[genotype_class],
                            linestyle=linestyle,
                            linewidth=2.8,
                            label=f"{simulation_class.capitalize()} {genotype_class.replace('_', '/')}",
                        )
                axis.set_title(f"Final AF {100 * af:.0f}%", fontsize=18)
                axis.set_xlabel("TMRCA threshold (kya)", fontsize=17)
                axis.tick_params(labelsize=14)
                axis.grid(alpha=0.25)
            axes[0].set_ylabel("Mean P(TMRCA < threshold)", fontsize=17)
            handles, labels = axes[-1].get_legend_handles_labels()
            figure.legend(handles, labels, loc="lower center", ncol=2, fontsize=15)
            title = (
                "EAS median demography"
                if "eas_phlash" in str(demography)
                else "Han AncientEurasia introgression"
            )
            figure.suptitle(
                f"{display_source} focal TMRCA profiles: {title}, s={coefficient:g}",
                fontsize=21,
            )
            figure.subplots_adjust(bottom=0.22, top=0.88, wspace=0.18)
            safe = str(demography).replace("/", "_")
            coefficient_slug = str(coefficient).replace(".", "p")
            for extension in ("png", "pdf"):
                path = (
                    figure_dir
                    / f"probability_curves_{safe}_s{coefficient_slug}.{extension}"
                )
                _atomic_figure(figure, path)
                outputs.append(path)
            plt.close(figure)
    return outputs


def plot_spatial_probability_profiles(
    spatial_summaries: pd.DataFrame,
    figure_dir: Path,
    *,
    source_label: str | None = None,
) -> list[Path]:
    """Plot genomic P(TMRCA<x) profiles for every class, AF, s, and threshold."""

    display_source = source_label or _source_display_name(spatial_summaries)
    outputs: list[Path] = []
    colors = {
        "overall": "#000000",
        "hom_ref": "#595959",
        "heterozygous": "#E69F00",
        "hom_alt": "#0072B2",
    }
    for demography, demographic in spatial_summaries.groupby(
        "demography_id", sort=True
    ):
        coefficients = sorted(
            demographic.loc[
                demographic["simulation_class"] == "selected",
                "selection_coefficient",
            ]
            .astype(float)
            .unique()
        )
        for coefficient in coefficients:
            comparison = demographic[
                (demographic["simulation_class"] == "neutral")
                | (
                    (demographic["simulation_class"] == "selected")
                    & np.isclose(
                        demographic["selection_coefficient"], float(coefficient)
                    )
                )
            ]
            afs = sorted(comparison["target_allele_frequency"].astype(float).unique())
            for threshold in TMRCA_THRESHOLDS_YEARS:
                threshold_frame = comparison[
                    np.isclose(comparison["threshold_years"], float(threshold))
                ]
                figure, axes = plt.subplots(
                    1, len(afs), figsize=(11, 8.5), sharex=True, sharey=True
                )
                axes = np.atleast_1d(axes)
                for axis, af in zip(axes, afs, strict=True):
                    cell = threshold_frame[
                        np.isclose(threshold_frame["target_allele_frequency"], af)
                    ]
                    for simulation_class, linestyle in (
                        ("selected", "-"),
                        ("neutral", "--"),
                    ):
                        for genotype_class in GENOTYPE_CLASSES:
                            curve = cell[
                                (cell["simulation_class"] == simulation_class)
                                & (cell["genotype_class"] == genotype_class)
                            ].sort_values("position_0based")
                            if curve.empty:
                                raise ValueError(
                                    "spatial plot grid lacks a selected/neutral class"
                                )
                            x = curve["position_0based"].to_numpy(dtype=float) / 1e6
                            y = curve["mean_p_tmrca_lt_threshold"].to_numpy(dtype=float)
                            sem = curve["sem_p_tmrca_lt_threshold"].to_numpy(
                                dtype=float
                            )
                            color = colors[genotype_class]
                            axis.plot(
                                x,
                                y,
                                color=color,
                                linestyle=linestyle,
                                linewidth=2.2,
                                label=(
                                    f"{simulation_class.capitalize()} "
                                    f"{genotype_class.replace('_', '/')}"
                                ),
                            )
                            axis.fill_between(
                                x,
                                np.clip(y - 1.96 * sem, 0, 1),
                                np.clip(y + 1.96 * sem, 0, 1),
                                color=color,
                                alpha=0.07 if simulation_class == "selected" else 0.035,
                            )
                    axis.axvline(
                        FOCAL_POSITION_BP / 1e6,
                        color="#CC3311",
                        linewidth=1.4,
                        alpha=0.8,
                    )
                    axis.set_title(f"Final AF {100 * af:.0f}%", fontsize=18)
                    axis.set_xlabel("Position in 10-Mb region (Mb)", fontsize=16)
                    axis.tick_params(labelsize=13)
                    axis.grid(alpha=0.20)
                axes[0].set_ylabel(
                    f"Mean P(TMRCA < {threshold / 1_000:g} kya)", fontsize=16
                )
                handles, labels = axes[-1].get_legend_handles_labels()
                figure.legend(
                    handles,
                    labels,
                    loc="lower center",
                    ncol=2,
                    fontsize=14,
                )
                title = (
                    "EAS median demography"
                    if "eas_phlash" in str(demography)
                    else "Han AncientEurasia introgression"
                )
                figure.suptitle(
                    f"{display_source} genomic profiles: {title}, "
                    f"s={coefficient:g}, x={threshold / 1_000:g} kya",
                    fontsize=20,
                )
                figure.subplots_adjust(bottom=0.23, top=0.88, wspace=0.18)
                safe = str(demography).replace("/", "_")
                coefficient_slug = str(coefficient).replace(".", "p")
                threshold_slug = f"{float(threshold):g}".replace(".", "p")
                for extension in ("png", "pdf"):
                    path = figure_dir / (
                        f"spatial_probability_{safe}_s{coefficient_slug}_"
                        f"t{threshold_slug}.{extension}"
                    )
                    _atomic_figure(figure, path)
                    outputs.append(path)
                plt.close(figure)
    return outputs


def plot_primary_contrast_curves(
    scores: pd.DataFrame,
    figure_dir: Path,
    *,
    source_label: str = "Gamma-SMC",
) -> list[Path]:
    """Plot the prespecified alt-minus-ref score and replicate distributions."""

    primary = scores[
        (scores["statistic_scope"] == "focal_nearest")
        & (scores["metric"] == "p_tmrca_lt_threshold")
        & (scores["statistic"] == PRIMARY_STATISTIC)
    ]
    outputs: list[Path] = []
    for demography, demographic in primary.groupby("demography_id", sort=True):
        coefficients = sorted(
            demographic.loc[
                demographic["simulation_class"] == "selected",
                "selection_coefficient",
            ]
            .astype(float)
            .unique()
        )
        for coefficient in coefficients:
            comparison = demographic[
                (demographic["simulation_class"] == "neutral")
                | (
                    (demographic["simulation_class"] == "selected")
                    & np.isclose(demographic["selection_coefficient"], coefficient)
                )
            ]
            afs = sorted(comparison["target_allele_frequency"].astype(float).unique())
            figure, axes = plt.subplots(1, len(afs), figsize=(11, 8.5), sharey=True)
            axes = np.atleast_1d(axes)
            for axis, af in zip(axes, afs, strict=True):
                cell = comparison[np.isclose(comparison["target_allele_frequency"], af)]
                for simulation_class, color in (
                    ("selected", "#D55E00"),
                    ("neutral", "#666666"),
                ):
                    selected = cell[cell["simulation_class"] == simulation_class]
                    aggregate = selected.groupby("threshold_years")[
                        "selection_score"
                    ].agg(
                        mean="mean",
                        low=lambda values: np.quantile(values, 0.025),
                        high=lambda values: np.quantile(values, 0.975),
                    )
                    x = aggregate.index.to_numpy(dtype=float) / 1_000
                    axis.plot(
                        x,
                        aggregate["mean"],
                        color=color,
                        linewidth=3,
                        label=simulation_class.capitalize(),
                    )
                    axis.fill_between(
                        x,
                        aggregate["low"],
                        aggregate["high"],
                        color=color,
                        alpha=0.16,
                    )
                axis.axhline(0, color="black", linewidth=1, alpha=0.6)
                axis.set_title(f"Final AF {100 * af:.0f}%", fontsize=18)
                axis.set_xlabel("TMRCA threshold (kya)", fontsize=17)
                axis.tick_params(labelsize=14)
                axis.grid(alpha=0.25)
            axes[0].set_ylabel("P(T<x | alt/alt) - P(T<x | ref/ref)", fontsize=16)
            handles, labels = axes[-1].get_legend_handles_labels()
            figure.legend(handles, labels, loc="lower center", ncol=2, fontsize=15)
            title = (
                "EAS median demography"
                if "eas_phlash" in str(demography)
                else "Han AncientEurasia introgression"
            )
            figure.suptitle(
                f"Primary {source_label} selection contrast: {title}, s={coefficient:g}",
                fontsize=21,
            )
            figure.subplots_adjust(bottom=0.18, top=0.88, wspace=0.20)
            safe = str(demography).replace("/", "_")
            coefficient_slug = str(coefficient).replace(".", "p")
            for extension in ("png", "pdf"):
                path = (
                    figure_dir
                    / f"primary_contrast_curves_{safe}_s{coefficient_slug}.{extension}"
                )
                _atomic_figure(figure, path)
                outputs.append(path)
            plt.close(figure)
    return outputs


def plot_primary_auc(
    cell_comparisons: pd.DataFrame,
    figure_dir: Path,
    *,
    source_label: str = "Gamma-SMC",
) -> list[Path]:
    """Plot threshold-specific ROC AUC for the primary within-panel contrast."""

    selected = cell_comparisons[
        (cell_comparisons["metric"] == "p_tmrca_lt_threshold")
        & (cell_comparisons["statistic"] == PRIMARY_STATISTIC)
        & (cell_comparisons["statistic_scope"] == "focal_nearest")
    ].copy()
    selected["cell"] = [
        _cell_label(str(d), float(af), float(s))
        for d, af, s in zip(
            selected["demography_id"],
            selected["target_allele_frequency"],
            selected["selection_coefficient"],
            strict=True,
        )
    ]
    matrix = selected.pivot(index="cell", columns="threshold_years", values="roc_auc")
    matrix = matrix.sort_index()
    figure, axis = plt.subplots(figsize=(11, 8.5))
    image = axis.imshow(
        matrix.to_numpy(dtype=float), vmin=0, vmax=1, cmap="viridis", aspect="auto"
    )
    axis.set_xticks(
        np.arange(len(matrix.columns)),
        [f"{x / 1000:g}" for x in matrix.columns],
        fontsize=14,
    )
    axis.set_yticks(np.arange(len(matrix.index)), matrix.index, fontsize=14)
    axis.set_xlabel("TMRCA threshold (kya)", fontsize=18)
    axis.set_title(
        f"{source_label} ROC AUC: selected vs neutral primary contrast", fontsize=21
    )
    for row in range(len(matrix.index)):
        for column in range(len(matrix.columns)):
            value = float(matrix.iloc[row, column])
            axis.text(
                column,
                row,
                f"{value:.2f}",
                ha="center",
                va="center",
                color="white" if value < 0.75 else "black",
                fontsize=13,
            )
    colorbar = figure.colorbar(image, ax=axis, pad=0.02)
    colorbar.set_label("ROC AUC", fontsize=17)
    colorbar.ax.tick_params(labelsize=13)
    outputs: list[Path] = []
    for extension in ("png", "pdf"):
        path = figure_dir / f"primary_contrast_roc_auc.{extension}"
        _atomic_figure(figure, path)
        outputs.append(path)
    plt.close(figure)
    return outputs


def plot_primary_pvalues(
    cell_comparisons: pd.DataFrame,
    figure_dir: Path,
    *,
    source_label: str = "Gamma-SMC",
) -> list[Path]:
    """Plot one-sided Mann-Whitney evidence for the primary contrast."""

    selected = cell_comparisons[
        (cell_comparisons["metric"] == "p_tmrca_lt_threshold")
        & (cell_comparisons["statistic"] == PRIMARY_STATISTIC)
        & (cell_comparisons["statistic_scope"] == "focal_nearest")
    ].copy()
    selected["cell"] = [
        _cell_label(str(d), float(af), float(s))
        for d, af, s in zip(
            selected["demography_id"],
            selected["target_allele_frequency"],
            selected["selection_coefficient"],
            strict=True,
        )
    ]
    selected["minus_log10_p"] = -np.log10(
        np.maximum(selected["mann_whitney_p_upper"].astype(float), 1e-300)
    )
    matrix = selected.pivot(
        index="cell", columns="threshold_years", values="minus_log10_p"
    ).sort_index()
    figure, axis = plt.subplots(figsize=(11, 8.5))
    image = axis.imshow(
        matrix.to_numpy(dtype=float), vmin=0, cmap="magma", aspect="auto"
    )
    axis.set_xticks(
        np.arange(len(matrix.columns)),
        [f"{x / 1000:g}" for x in matrix.columns],
        fontsize=14,
    )
    axis.set_yticks(np.arange(len(matrix.index)), matrix.index, fontsize=14)
    axis.set_xlabel("TMRCA threshold (kya)", fontsize=18)
    axis.set_title(
        f"{source_label} Mann-Whitney evidence: selected > neutral", fontsize=21
    )
    for row in range(len(matrix.index)):
        for column in range(len(matrix.columns)):
            value = float(matrix.iloc[row, column])
            axis.text(
                column,
                row,
                f"{value:.2f}",
                ha="center",
                va="center",
                color="white",
                fontsize=13,
            )
    colorbar = figure.colorbar(image, ax=axis, pad=0.02)
    colorbar.set_label("-log10(one-sided Mann-Whitney p)", fontsize=17)
    colorbar.ax.tick_params(labelsize=13)
    outputs: list[Path] = []
    for extension in ("png", "pdf"):
        path = figure_dir / f"primary_contrast_pvalues.{extension}"
        _atomic_figure(figure, path)
        outputs.append(path)
    plt.close(figure)
    return outputs


def _output_record(
    path: Path, root: Path, frame: pd.DataFrame | None = None
) -> dict[str, Any]:
    record: dict[str, Any] = {
        "path": str(path.relative_to(root)).replace("\\", "/"),
        "sha256": sha256_file(path),
        "size_bytes": path.stat().st_size,
    }
    if frame is not None:
        record["rows"] = len(frame)
        record["columns"] = list(frame.columns)
    return record


def _validate_cached_analysis(
    output_dir: Path,
    completion: Mapping[str, Any],
    contract: Mapping[str, Any],
) -> dict[str, Path]:
    if completion.get("schema") != SCHEMA_VERSION:
        raise ValueError("analysis completion schema is incompatible")
    if completion.get("status") != "complete":
        raise ValueError("analysis completion status is not complete")
    if completion.get("contract_sha256") != _canonical_sha256(
        completion.get("contract")
    ):
        raise ValueError("analysis completion contract checksum is invalid")
    if _canonical_sha256(completion.get("contract")) != _canonical_sha256(contract):
        raise ValueError("existing analysis contract is incompatible")
    outputs = completion.get("outputs")
    if not isinstance(outputs, Mapping) or not outputs:
        raise ValueError("analysis completion output manifest is malformed")
    paths: dict[str, Path] = {}
    for label, record in outputs.items():
        path = output_dir / str(record.get("path", ""))
        paths[str(label)] = path
        if not path.is_file() or sha256_file(path) != record.get("sha256"):
            raise ValueError(f"analysis output checksum failed: {label}")
        if path.stat().st_size != int(record.get("size_bytes", -1)):
            raise ValueError(f"analysis output size failed: {label}")
        if "rows" in record:
            try:
                frame = pd.read_csv(path, sep="\t")
            except pd.errors.EmptyDataError:
                frame = pd.DataFrame()
            except Exception as error:
                raise ValueError(f"analysis output is unreadable: {label}") from error
            if len(frame) != int(record["rows"]) or list(frame.columns) != list(
                record.get("columns", [])
            ):
                raise ValueError(f"analysis output schema failed: {label}")
    paths["completion"] = output_dir / "analysis_completion.json"
    return paths


def run_focused_analysis(
    class_summaries_path: str | Path,
    output_dir: str | Path,
    *,
    pair_summaries_path: str | Path | None = None,
    spatial_summaries_path: str | Path | None = None,
    empirical_path: str | Path | None = None,
    minimum_selected_per_cell: int = 10,
    minimum_neutral_per_cell: int = 100,
    bootstrap_draws: int = DEFAULT_BOOTSTRAP_DRAWS,
    permutation_draws: int = DEFAULT_LABEL_PERMUTATION_DRAWS,
    crossfit_folds: int = DEFAULT_CROSSFIT_FOLDS,
    base_seed: int = BASE_SEED,
    make_plots: bool = True,
) -> dict[str, Path]:
    """Run and checksum-contract the full focused statistical analysis."""

    class_summaries_path = Path(class_summaries_path).resolve()
    output_dir = Path(output_dir).resolve()
    pair_summaries = (
        Path(pair_summaries_path).resolve() if pair_summaries_path is not None else None
    )
    spatial_summaries = (
        Path(spatial_summaries_path).resolve()
        if spatial_summaries_path is not None
        else None
    )
    empirical = Path(empirical_path).resolve() if empirical_path is not None else None
    if not class_summaries_path.is_file():
        raise ValueError("combined class summaries are absent")
    if empirical is not None and not empirical.is_file():
        raise ValueError("requested empirical TSV is absent")
    if pair_summaries is not None and not pair_summaries.is_file():
        raise ValueError("requested combined pair summaries are absent")
    if spatial_summaries is not None and not spatial_summaries.is_file():
        raise ValueError("requested aggregated spatial summaries are absent")
    if int(bootstrap_draws) < 1:
        raise ValueError("bootstrap_draws must be positive")
    if int(permutation_draws) < 1:
        raise ValueError("permutation_draws must be positive")
    if int(crossfit_folds) < 2:
        raise ValueError("crossfit_folds must be at least two")
    contract = {
        "schema": SCHEMA_VERSION,
        "implementation": _analysis_implementation_contract(),
        "inputs": {
            "class_summaries_path": str(class_summaries_path),
            "class_summaries_sha256": sha256_file(class_summaries_path),
            "pair_summaries_path": (
                str(pair_summaries) if pair_summaries is not None else None
            ),
            "pair_summaries_sha256": (
                sha256_file(pair_summaries) if pair_summaries is not None else None
            ),
            "spatial_summaries_path": (
                str(spatial_summaries) if spatial_summaries is not None else None
            ),
            "spatial_summaries_sha256": (
                sha256_file(spatial_summaries)
                if spatial_summaries is not None
                else None
            ),
            "empirical_path": str(empirical) if empirical is not None else None,
            "empirical_sha256": sha256_file(empirical)
            if empirical is not None
            else None,
        },
        "parameters": {
            "thresholds_years": list(TMRCA_THRESHOLDS_YEARS),
            "minimum_selected_per_cell": int(minimum_selected_per_cell),
            "minimum_neutral_per_cell": int(minimum_neutral_per_cell),
            "bootstrap_draws": int(bootstrap_draws),
            "unit_label_permutation_draws": int(permutation_draws),
            "neutral_crossfit_folds": int(crossfit_folds),
            "base_seed": int(base_seed),
            "make_plots": bool(make_plots),
            "primary_statistic": PRIMARY_STATISTIC,
            "monte_carlo_pvalue": "one-sided upper-tail with +1 correction",
            "bh_family": "seven thresholds within cell/scope/metric/statistic",
            "primary_omnibus": (
                "pooled leave-one-out minP retaining whole unit threshold vectors"
            ),
            "neutral_calibration": (
                f"deterministic {int(crossfit_folds)}-fold held-out units"
            ),
        },
    }
    completion_path = output_dir / "analysis_completion.json"
    if completion_path.is_file():
        try:
            completion = json.loads(completion_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            raise ValueError("analysis completion is unreadable") from error
        return _validate_cached_analysis(output_dir, completion, contract)
    if output_dir.exists() and any(output_dir.iterdir()):
        raise ValueError(
            "analysis directory is nonempty without completion; quarantine or remove "
            "the partial analysis before retrying"
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    try:
        raw = pd.read_csv(class_summaries_path, sep="\t")
    except Exception as error:
        raise ValueError("combined class summaries are unreadable") from error
    classes = validate_class_summaries(
        raw,
        minimum_selected_per_cell=minimum_selected_per_cell,
        minimum_neutral_per_cell=minimum_neutral_per_cell,
    )
    spatial_frame = (
        validate_spatial_summaries(pd.read_csv(spatial_summaries, sep="\t"), classes)
        if spatial_summaries is not None
        else None
    )
    scores = build_replicate_statistics(classes)
    selected_pvalues = selected_vs_neutral_pvalues(scores)
    power = conditional_power_summary(selected_pvalues)
    loo, calibration = neutral_leave_one_out(scores)
    crossfit, crossfit_calibration = neutral_crossfit_calibration(
        scores, n_folds=crossfit_folds, base_seed=base_seed
    )
    omnibus = primary_threshold_omnibus(scores)
    comparisons = cell_roc_comparisons(
        scores,
        bootstrap_draws=bootstrap_draws,
        permutation_draws=permutation_draws,
        base_seed=base_seed,
    )
    within = within_selected_summary(
        scores, bootstrap_draws=bootstrap_draws, base_seed=base_seed
    )
    pair_frame = (
        pd.read_csv(pair_summaries, sep="\t")
        if pair_summaries is not None
        else pd.DataFrame()
    )
    unit_pair_tests = (
        within_selected_unit_pair_tests(pair_frame, base_seed=base_seed)
        if pair_summaries is not None
        else pd.DataFrame()
    )
    template = empirical_input_template(classes)
    empirical_frame = (
        pd.read_csv(empirical, sep="\t") if empirical is not None else None
    )
    empirical_results = empirical_null_pvalues(empirical_frame, scores)
    frames = {
        "replicate_statistics": scores,
        "selected_pvalues": selected_pvalues,
        "within_selected": within,
        "cell_comparisons": comparisons,
        "neutral_loo": loo,
        "neutral_calibration": calibration,
        "neutral_crossfit": crossfit,
        "neutral_crossfit_calibration": crossfit_calibration,
        "primary_omnibus": omnibus,
        "empirical_pvalues": empirical_results,
        "empirical_template": template,
        "unit_pair_tests": unit_pair_tests,
        "conditional_power": power,
    }
    paths: dict[str, Path] = {}
    for label, filename in TABLE_FILENAMES.items():
        path = output_dir / filename
        _atomic_frame(path, frames[label])
        paths[label] = path
    figure_paths: list[Path] = []
    if make_plots:
        figure_dir = output_dir / "figures"
        source_label = _source_display_name(classes)
        figure_paths.extend(
            plot_probability_curves(classes, figure_dir, source_label=source_label)
        )
        if spatial_frame is not None:
            figure_paths.extend(
                plot_spatial_probability_profiles(
                    spatial_frame, figure_dir, source_label=source_label
                )
            )
        figure_paths.extend(
            plot_primary_contrast_curves(scores, figure_dir, source_label=source_label)
        )
        figure_paths.extend(
            plot_primary_auc(comparisons, figure_dir, source_label=source_label)
        )
        figure_paths.extend(
            plot_primary_pvalues(comparisons, figure_dir, source_label=source_label)
        )
    outputs: dict[str, Any] = {
        label: _output_record(paths[label], output_dir, frames[label])
        for label in paths
    }
    for index, path in enumerate(figure_paths):
        label = f"figure_{index:02d}_{path.stem}_{path.suffix.removeprefix('.')}"
        outputs[label] = _output_record(path, output_dir)
        paths[label] = path
    completion = {
        "schema": SCHEMA_VERSION,
        "status": "complete",
        "contract": contract,
        "contract_sha256": _canonical_sha256(contract),
        "outputs": outputs,
        "interpretation": {
            "primary_score": (
                "P(T<x | hom_alt) minus P(T<x | hom_ref), matched to the same "
                "contrast in neutral simulations"
            ),
            "roc_auc": (
                "probability that a randomly drawn selected replicate has a larger "
                "selection-directed score than a randomly drawn neutral replicate"
            ),
            "normalized_cdf_area_is_roc_auc": False,
            "empirical_status": "provided"
            if empirical is not None
            else "template_only",
            "within_unit_pair_tests": (
                "provided"
                if pair_summaries is not None
                else "not_run_pair_summaries_not_supplied"
            ),
            "spatial_profiles": (
                "provided" if spatial_summaries is not None else "not_provided"
            ),
            "minimum_attainable_empirical_p": float(
                1.0 / (minimum_neutral_per_cell + 1)
            ),
            "neutral_calibration_primary": (
                f"deterministic_{int(crossfit_folds)}_fold_crossfit"
            ),
            "neutral_loo_role": "descriptive_only",
            "threshold_omnibus": (
                "correlation-aware pooled leave-one-out minP across seven thresholds"
            ),
        },
    }
    _atomic_json(completion_path, completion)
    paths["completion"] = completion_path
    return paths


def run_truth_gamma_comparison(
    truth_class_summaries_path: str | Path,
    gamma_class_summaries_path: str | Path,
    output_dir: str | Path,
    *,
    minimum_selected_per_cell: int = 10,
    minimum_neutral_per_cell: int = 100,
) -> dict[str, Path]:
    """Compare identical tree-truth and Gamma-SMC unit-level estimands."""

    truth_path = Path(truth_class_summaries_path).resolve()
    gamma_path = Path(gamma_class_summaries_path).resolve()
    destination = Path(output_dir).resolve()
    if not truth_path.is_file() or not gamma_path.is_file():
        raise ValueError("truth and Gamma class-summary inputs are required")
    contract = {
        "schema": SCHEMA_VERSION,
        "implementation": _analysis_implementation_contract(),
        "inputs": {
            "truth_sha256": sha256_file(truth_path),
            "gamma_sha256": sha256_file(gamma_path),
        },
        "parameters": {
            "minimum_selected_per_cell": int(minimum_selected_per_cell),
            "minimum_neutral_per_cell": int(minimum_neutral_per_cell),
        },
    }
    completion_path = destination / "truth_gamma_completion.json"
    if completion_path.is_file():
        try:
            completion = json.loads(completion_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            raise ValueError("truth-Gamma completion is unreadable") from error
        if (
            completion.get("schema") != SCHEMA_VERSION
            or completion.get("status") != "complete"
            or completion.get("contract_sha256") != _canonical_sha256(contract)
        ):
            raise ValueError("truth-Gamma completion contract is incompatible")
        paths: dict[str, Path] = {"completion": completion_path}
        for label, record in completion.get("outputs", {}).items():
            records = record if isinstance(record, list) else [record]
            for index, item in enumerate(records):
                if not isinstance(item, Mapping):
                    raise TypeError(
                        f"truth-Gamma output manifest is malformed: {label}"
                    )
                path = destination / str(item.get("path", ""))
                if not path.is_file() or sha256_file(path) != item.get("sha256"):
                    raise ValueError(f"truth-Gamma output checksum failed: {label}")
                key = str(label) if len(records) == 1 else f"{label}_{index}"
                paths[key] = path
        return paths
    if destination.exists() and any(destination.iterdir()):
        raise ValueError(
            "truth-Gamma directory is nonempty without completion; quarantine or "
            "remove the partial analysis before retrying"
        )
    truth = validate_class_summaries(
        pd.read_csv(truth_path, sep="\t"),
        minimum_selected_per_cell=minimum_selected_per_cell,
        minimum_neutral_per_cell=minimum_neutral_per_cell,
    )
    gamma = validate_class_summaries(
        pd.read_csv(gamma_path, sep="\t"),
        minimum_selected_per_cell=minimum_selected_per_cell,
        minimum_neutral_per_cell=minimum_neutral_per_cell,
    )
    truth_scores = build_replicate_statistics(truth)
    gamma_scores = build_replicate_statistics(gamma)
    key = [
        *UNIT_COLUMNS,
        "statistic_scope",
        "metric",
        "statistic",
        "threshold_years",
    ]
    paired = truth_scores[key + ["selection_score"]].merge(
        gamma_scores[key + ["selection_score"]],
        on=key,
        how="outer",
        validate="one_to_one",
        suffixes=("_tree_truth", "_gamma_smc"),
        indicator=True,
    )
    if set(paired["_merge"]) != {"both"}:
        raise ValueError("tree truth and Gamma-SMC unit estimands do not align exactly")
    paired = paired.drop(columns="_merge")
    paired["gamma_minus_truth"] = (
        paired["selection_score_gamma_smc"] - paired["selection_score_tree_truth"]
    )
    paired["absolute_error"] = paired["gamma_minus_truth"].abs()
    paired["squared_error"] = paired["gamma_minus_truth"] ** 2
    paired["truth_positive"] = paired["selection_score_tree_truth"] > 0
    paired["gamma_positive"] = paired["selection_score_gamma_smc"] > 0
    paired["sign_agreement"] = paired["truth_positive"] == paired["gamma_positive"]
    group_columns = [
        "simulation_class",
        *CELL_COLUMNS,
        "selection_coefficient",
        "statistic_scope",
        "metric",
        "statistic",
        "threshold_years",
    ]
    rows: list[dict[str, Any]] = []
    for labels, group in paired.groupby(group_columns, sort=True, dropna=False):
        tree = group["selection_score_tree_truth"].to_numpy(dtype=float)
        gamma = group["selection_score_gamma_smc"].to_numpy(dtype=float)
        pearson = (
            pearsonr(tree, gamma).statistic
            if len(group) >= 2 and np.std(tree) > 0 and np.std(gamma) > 0
            else np.nan
        )
        spearman = (
            spearmanr(tree, gamma).statistic
            if len(group) >= 2 and np.std(tree) > 0 and np.std(gamma) > 0
            else np.nan
        )
        rows.append(
            {
                **dict(zip(group_columns, labels, strict=True)),
                "n_units": len(group),
                "mean_tree_truth": float(np.mean(tree)),
                "mean_gamma_smc": float(np.mean(gamma)),
                "bias_gamma_minus_truth": float(np.mean(gamma - tree)),
                "mae": float(np.mean(np.abs(gamma - tree))),
                "rmse": float(np.sqrt(np.mean((gamma - tree) ** 2))),
                "pearson_r": float(pearson),
                "spearman_rho": float(spearman),
                "sign_agreement_fraction": float(group["sign_agreement"].mean()),
                "truth_positive_fraction": float(group["truth_positive"].mean()),
                "gamma_positive_fraction": float(group["gamma_positive"].mean()),
            }
        )
    summary = pd.DataFrame(rows)
    destination.mkdir(parents=True, exist_ok=True)
    paired_path = destination / "truth_gamma_paired_statistics.tsv.gz"
    summary_path = destination / "truth_gamma_accuracy.tsv"
    _atomic_frame(paired_path, paired)
    _atomic_frame(summary_path, summary)
    primary = paired[
        (paired["simulation_class"] == "selected")
        & (paired["statistic"] == PRIMARY_STATISTIC)
        & (paired["metric"] == "p_tmrca_lt_threshold")
    ]
    figure_paths: list[Path] = []
    if not primary.empty:
        for scope, scoped in primary.groupby("statistic_scope", sort=True):
            figure, axis = plt.subplots(figsize=(11, 8.5))
            for (demography, coefficient), group in scoped.groupby(
                ["demography_id", "selection_coefficient"], sort=True
            ):
                axis.scatter(
                    group["selection_score_tree_truth"],
                    group["selection_score_gamma_smc"],
                    s=55,
                    alpha=0.65,
                    label=f"{demography}, s={coefficient:g}",
                )
            limits = np.asarray(
                [
                    scoped["selection_score_tree_truth"].min(),
                    scoped["selection_score_tree_truth"].max(),
                    scoped["selection_score_gamma_smc"].min(),
                    scoped["selection_score_gamma_smc"].max(),
                ]
            )
            low, high = float(limits.min()), float(limits.max())
            axis.plot([low, high], [low, high], color="black", linestyle="--")
            axis.axhline(0, color="grey", linewidth=1)
            axis.axvline(0, color="grey", linewidth=1)
            axis.set_xlabel("Tree-truth alt/alt minus ref/ref P(T<x)", fontsize=16)
            axis.set_ylabel("Gamma-SMC alt/alt minus ref/ref P(T<x)", fontsize=16)
            axis.set_title(f"Truth versus Gamma-SMC: {scope}", fontsize=19)
            axis.tick_params(labelsize=13)
            axis.legend(fontsize=14, loc="best")
            figure.tight_layout()
            for extension in ("png", "pdf"):
                path = destination / f"truth_gamma_primary_{scope}.{extension}"
                _atomic_figure(figure, path)
                figure_paths.append(path)
            plt.close(figure)
    completion = {
        "schema": SCHEMA_VERSION,
        "status": "complete",
        "contract": contract,
        "contract_sha256": _canonical_sha256(contract),
        "outputs": {
            "paired": {"path": paired_path.name, "sha256": sha256_file(paired_path)},
            "summary": {"path": summary_path.name, "sha256": sha256_file(summary_path)},
            "figures": [
                {"path": path.name, "sha256": sha256_file(path)}
                for path in figure_paths
            ],
        },
    }
    _atomic_json(completion_path, completion)
    return {
        "paired": paired_path,
        "summary": summary_path,
        "completion": completion_path,
        **{f"figure_{index}": path for index, path in enumerate(figure_paths)},
    }


__all__ = [
    "DEFAULT_BOOTSTRAP_DRAWS",
    "DEFAULT_CROSSFIT_FOLDS",
    "DEFAULT_LABEL_PERMUTATION_DRAWS",
    "EMPIRICAL_COLUMNS",
    "EMPIRICAL_PROVENANCE_COLUMNS",
    "EMPIRICAL_STATISTIC_COLUMNS",
    "PRIMARY_STATISTIC",
    "SCHEMA_VERSION",
    "build_replicate_statistics",
    "cell_roc_comparisons",
    "conditional_power_summary",
    "empirical_input_template",
    "empirical_null_pvalues",
    "neutral_crossfit_calibration",
    "neutral_leave_one_out",
    "plot_primary_auc",
    "plot_primary_contrast_curves",
    "plot_primary_pvalues",
    "plot_probability_curves",
    "plot_spatial_probability_profiles",
    "primary_threshold_omnibus",
    "run_focused_analysis",
    "run_truth_gamma_comparison",
    "selected_vs_neutral_pvalues",
    "validate_class_summaries",
    "validate_spatial_summaries",
    "within_selected_summary",
    "within_selected_unit_pair_tests",
]
