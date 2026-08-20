"""Statistics, null distribution and p-values for the run2 study.

The test statistic is deliberately the *overall* within-individual
``P(TMRCA < x)`` at the focal base.  The null contains no focal allele, so there
is no genotype class to condition on and any carrier-stratified statistic would
be uncomparable to it.  Carrier-stratified numbers are still carried through the
tables, clearly marked as descriptive.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

from .run2_config import (
    MODES,
    Run2Arm,
    SIGNIFICANCE_LEVEL,
    TMRCA_THRESHOLDS_YEARS,
)

TEST_CLASS = "overall"

RESULT_OUTPUTS = {
    "final_af": "final_af.tsv",
    "observed_statistics": "observed_statistics.tsv",
    "null_distribution": "null_distribution.tsv",
    "pvalues": "pvalues.tsv",
    "power": "power_summary.tsv",
    "null_calibration": "null_calibration.tsv",
    "spatial_profiles": "spatial_profiles.tsv.gz",
    "replicate_status": "replicate_status.tsv",
}

_NULL_QUANTILES = (0.5, 0.9, 0.95, 0.99)


# ---------------------------------------------------------------------------
# Collection
# ---------------------------------------------------------------------------


def _replicate_directories(study_root: str | Path, arm: Run2Arm, mode: str) -> list[Path]:
    base = Path(study_root) / arm.arm_id / mode / "replicates"
    if not base.is_dir():
        return []
    return sorted(path for path in base.iterdir() if path.is_dir())


def load_endpoints(study_root: str | Path, arm: Run2Arm, mode: str) -> list[dict[str, Any]]:
    endpoints = []
    for directory in _replicate_directories(study_root, arm, mode):
        path = directory / "endpoint.json"
        if path.is_file():
            endpoints.append(json.loads(path.read_text(encoding="utf-8")))
    return endpoints


def replicate_status(study_root: str | Path, arm: Run2Arm) -> pd.DataFrame:
    """Return one row per attempted replicate with its completion status."""

    rows = []
    for mode in MODES:
        for endpoint in load_endpoints(study_root, arm, mode):
            rows.append(
                {
                    "arm_id": arm.arm_id,
                    "mode": mode,
                    "unit_id": endpoint["unit_id"],
                    "replicate_index": endpoint["replicate_index"],
                    "seed": endpoint["seed"],
                    "status": endpoint["status"],
                    "elapsed_seconds": endpoint.get("elapsed_seconds"),
                    "error": endpoint.get("error", ""),
                }
            )
    return pd.DataFrame(rows).sort_values(["mode", "replicate_index"], ignore_index=True)


def final_allele_frequencies(study_root: str | Path, arm: Run2Arm) -> pd.DataFrame:
    """Return the per-replicate final allele frequency table for the selected mode."""

    rows = []
    for endpoint in load_endpoints(study_root, arm, "selected"):
        if endpoint["status"] != "completed":
            continue
        final = endpoint["final_allele_frequency"]
        standing = endpoint["standing_variation"]
        census_af = final.get("census_af")
        census_af = float(census_af) if census_af is not None else float("nan")
        rows.append(
            {
                "arm_id": arm.arm_id,
                "unit_id": endpoint["unit_id"],
                "replicate_index": endpoint["replicate_index"],
                "seed": endpoint["seed"],
                "standing_carrier_genomes": standing.get("carrier_genomes"),
                "standing_population_genomes": standing.get("population_genomes"),
                "standing_realized_frequency": standing.get("realized_frequency"),
                "census_alt_count": final.get("census_alt_count"),
                "census_total_count": final.get("census_total_count"),
                "census_af": census_af,
                "census_state": final.get("census_state"),
                "sample_af": final.get("sample_af"),
                "sample_carrier_diploids": final.get("sample_carrier_diploids"),
                "lost": bool(census_af == 0.0),
                "fixed": bool(census_af == 1.0),
                "segregating": bool(0.0 < census_af < 1.0),
            }
        )
    return pd.DataFrame(rows).sort_values("replicate_index", ignore_index=True)


def observed_statistics(study_root: str | Path, arm: Run2Arm) -> pd.DataFrame:
    """Return the long-form focal statistic for every completed replicate."""

    frames = []
    for mode in MODES:
        for directory in _replicate_directories(study_root, arm, mode):
            endpoint_path = directory / "endpoint.json"
            scores_path = directory / "focal_scores.tsv"
            if not (endpoint_path.is_file() and scores_path.is_file()):
                continue
            endpoint = json.loads(endpoint_path.read_text(encoding="utf-8"))
            if endpoint["status"] != "completed":
                continue
            frame = pd.read_csv(scores_path, sep="\t")
            frame.insert(0, "arm_id", arm.arm_id)
            frame.insert(1, "mode", mode)
            frame["seed"] = endpoint["seed"]
            frames.append(frame)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True).sort_values(
        ["mode", "replicate_index", "genotype_class", "threshold_years"],
        ignore_index=True,
    )


def aggregate_spatial_profiles(study_root: str | Path, arm: Run2Arm) -> pd.DataFrame:
    """Return per-position mean and quantile profiles for each mode.

    Per-replicate profiles stay in the replicate directories; only this
    aggregate is a study result.
    """
    frames = []
    for mode in MODES:
        stack = []
        for directory in _replicate_directories(study_root, arm, mode):
            path = directory / "spatial_profile.tsv.gz"
            if path.is_file():
                stack.append(pd.read_csv(path, sep="\t"))
        if not stack:
            continue
        positions = stack[0]["position_0based"].to_numpy()
        for threshold_years in TMRCA_THRESHOLDS_YEARS:
            column = f"p_lt_{int(threshold_years)}y"
            values = np.vstack([frame[column].to_numpy() for frame in stack])
            frames.append(
                pd.DataFrame(
                    {
                        "arm_id": arm.arm_id,
                        "mode": mode,
                        "threshold_years": float(threshold_years),
                        "position_0based": positions,
                        "n_replicates": values.shape[0],
                        "mean": values.mean(axis=0),
                        "median": np.median(values, axis=0),
                        "q05": np.quantile(values, 0.05, axis=0),
                        "q95": np.quantile(values, 0.95, axis=0),
                    }
                )
            )
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


# ---------------------------------------------------------------------------
# Null and p-values
# ---------------------------------------------------------------------------


def _test_values(statistics: pd.DataFrame, mode: str) -> pd.DataFrame:
    subset = statistics[
        (statistics["mode"] == mode) & (statistics["genotype_class"] == TEST_CLASS)
    ]
    return subset[
        ["replicate_index", "unit_id", "threshold_years", "p_tmrca_lt_threshold"]
    ].reset_index(drop=True)


def null_distribution(statistics: pd.DataFrame, arm: Run2Arm) -> pd.DataFrame:
    """Summarize the neutral statistic per threshold."""

    neutral = _test_values(statistics, "neutral")
    rows = []
    for threshold_years, group in neutral.groupby("threshold_years"):
        values = group["p_tmrca_lt_threshold"].to_numpy(dtype=float)
        row: dict[str, Any] = {
            "arm_id": arm.arm_id,
            "threshold_years": float(threshold_years),
            "n_null": int(values.size),
            "mean": float(values.mean()) if values.size else float("nan"),
            "sd": float(values.std(ddof=1)) if values.size > 1 else float("nan"),
            "min": float(values.min()) if values.size else float("nan"),
            "max": float(values.max()) if values.size else float("nan"),
        }
        for quantile in _NULL_QUANTILES:
            row[f"q{int(quantile * 100):02d}"] = (
                float(np.quantile(values, quantile)) if values.size else float("nan")
            )
        rows.append(row)
    return pd.DataFrame(rows).sort_values("threshold_years", ignore_index=True)


def upper_tail_pvalue(observed: float, null_values: np.ndarray) -> float:
    """Return ``(1 + #{null >= observed}) / (1 + n_null)``."""
    if not np.isfinite(observed) or null_values.size == 0:
        return float("nan")
    exceed = int(np.sum(null_values >= observed))
    return (1.0 + exceed) / (1.0 + null_values.size)


def pvalues(statistics: pd.DataFrame, arm: Run2Arm) -> pd.DataFrame:
    """Return one p-value per selected replicate per threshold."""

    selected = _test_values(statistics, "selected")
    neutral = _test_values(statistics, "neutral")
    rows = []
    for threshold_years, group in selected.groupby("threshold_years"):
        null_values = neutral.loc[
            neutral["threshold_years"] == threshold_years, "p_tmrca_lt_threshold"
        ].to_numpy(dtype=float)
        null_values = null_values[np.isfinite(null_values)]
        for record in group.itertuples(index=False):
            observed = float(record.p_tmrca_lt_threshold)
            p = upper_tail_pvalue(observed, null_values)
            rows.append(
                {
                    "arm_id": arm.arm_id,
                    "unit_id": record.unit_id,
                    "replicate_index": int(record.replicate_index),
                    "threshold_years": float(threshold_years),
                    "observed": observed,
                    "n_null": int(null_values.size),
                    "p_value": p,
                    "p_value_floor": 1.0 / (1.0 + null_values.size)
                    if null_values.size
                    else float("nan"),
                    "significant": bool(np.isfinite(p) and p <= SIGNIFICANCE_LEVEL),
                }
            )
    return pd.DataFrame(rows).sort_values(
        ["threshold_years", "replicate_index"], ignore_index=True
    )


def power_summary(pvalue_table: pd.DataFrame, arm: Run2Arm) -> pd.DataFrame:
    """Return the headline number: detection rate per threshold."""

    rows = []
    for threshold_years, group in pvalue_table.groupby("threshold_years"):
        finite = group[np.isfinite(group["p_value"])]
        rows.append(
            {
                "arm_id": arm.arm_id,
                "threshold_years": float(threshold_years),
                "n_selected": int(len(finite)),
                "n_significant": int(finite["significant"].sum()),
                "power": float(finite["significant"].mean()) if len(finite) else float("nan"),
                "median_p_value": float(finite["p_value"].median())
                if len(finite)
                else float("nan"),
                "median_observed": float(finite["observed"].median())
                if len(finite)
                else float("nan"),
                "significance_level": SIGNIFICANCE_LEVEL,
                "p_value_floor": float(finite["p_value_floor"].iloc[0])
                if len(finite)
                else float("nan"),
            }
        )
    return pd.DataFrame(rows).sort_values("threshold_years", ignore_index=True)


def null_calibration(statistics: pd.DataFrame, arm: Run2Arm) -> pd.DataFrame:
    """Leave-one-out false-positive rate of the null against itself.

    If the p-value machinery is sound this sits near the nominal level.  It
    costs nothing and it is the only check that distinguishes "the test works"
    from "the test always fires".
    """
    neutral = _test_values(statistics, "neutral")
    rows = []
    for threshold_years, group in neutral.groupby("threshold_years"):
        values = group["p_tmrca_lt_threshold"].to_numpy(dtype=float)
        values = values[np.isfinite(values)]
        if values.size < 2:
            continue
        leave_one_out = []
        for index in range(values.size):
            others = np.delete(values, index)
            leave_one_out.append(upper_tail_pvalue(values[index], others))
        array = np.asarray(leave_one_out, dtype=float)
        rows.append(
            {
                "arm_id": arm.arm_id,
                "threshold_years": float(threshold_years),
                "n_null": int(values.size),
                "false_positive_rate": float(np.mean(array <= SIGNIFICANCE_LEVEL)),
                "median_p_value": float(np.median(array)),
                "significance_level": SIGNIFICANCE_LEVEL,
            }
        )
    return pd.DataFrame(rows).sort_values("threshold_years", ignore_index=True)


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------


def analyse_arm(study_root: str | Path, arm: Run2Arm) -> dict[str, Path]:
    """Compute and write every result table for one arm."""

    results = Path(study_root) / arm.arm_id / "results"
    results.mkdir(parents=True, exist_ok=True)

    statistics = observed_statistics(study_root, arm)
    if statistics.empty:
        raise RuntimeError(f"no completed replicates found for arm {arm.arm_id}")

    tables = {
        "replicate_status": replicate_status(study_root, arm),
        "final_af": final_allele_frequencies(study_root, arm),
        "observed_statistics": statistics,
        "null_distribution": null_distribution(statistics, arm),
    }
    pvalue_table = pvalues(statistics, arm)
    tables["pvalues"] = pvalue_table
    tables["power"] = power_summary(pvalue_table, arm)
    tables["null_calibration"] = null_calibration(statistics, arm)
    tables["spatial_profiles"] = aggregate_spatial_profiles(study_root, arm)

    written: dict[str, Path] = {}
    for key, frame in tables.items():
        path = results / RESULT_OUTPUTS[key]
        frame.to_csv(path, sep="\t", index=False)
        written[key] = path
    return written


def cross_arm_comparison(study_root: str | Path, arms: Iterable[Run2Arm]) -> pd.DataFrame:
    """Stack the per-arm power curves for the cross-arm figure."""

    frames = []
    for arm in arms:
        path = Path(study_root) / arm.arm_id / "results" / RESULT_OUTPUTS["power"]
        if path.is_file():
            frame = pd.read_csv(path, sep="\t")
            frame["arm_label"] = arm.label
            frames.append(frame)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


__all__ = [
    "RESULT_OUTPUTS",
    "TEST_CLASS",
    "aggregate_spatial_profiles",
    "analyse_arm",
    "cross_arm_comparison",
    "final_allele_frequencies",
    "load_endpoints",
    "null_calibration",
    "null_distribution",
    "observed_statistics",
    "power_summary",
    "pvalues",
    "replicate_status",
    "upper_tail_pvalue",
]
