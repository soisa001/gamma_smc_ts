"""Analysis for run7: power, AUC, p-values and P(TMRCA < x) across the s grid.

The comparison follows what run6 established. The statistic is computed two ways
and both are carried through, because they answer different questions:

``unconditional``
    Every sampled individual, carrier or not. This is the statistic the
    empirical analysis can actually run against a simulated neutral null,
    because it has no ascertainment for the null to match.
``hom_carrier``
    Only individuals homozygous for the introgressed allele. More powerful, but
    ascertained, so an unconditional null scores it optimistically -- run6
    measured that inflation at 1.2-3x depending on frequency.

The null is the shared neutral arm, drawn at a *matched pair count* by
hypergeometric sampling so that a class holding 20 pairs is not compared with a
null built from 100.
"""

from __future__ import annotations

import glob
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd

from .run6_scan import build_matched_null
from .run7_config import (
    FOCAL_POSITION_BP,
    GENERATION_TIME_YEARS,
    SELECTION_COEFFICIENTS,
    THRESHOLDS_YEARS,
    arm_id_for,
)

__all__ = [
    "ALPHA",
    "ARM_SPECS",
    "MIN_CLASS_PAIRS",
    "Replicate",
    "focal_table",
    "load_arm",
    "load_mode",
    "scan_table",
    "summarise",
]

#: The de novo arms are a separate origin over their own coefficients. They share
#: the model, the onset and the neutral null with the introgressed arms, so the
#: only thing that differs is where the selected allele came from.
DENOVO_COEFFICIENTS: tuple[float, ...] = (
    0.001,
    0.002,
    0.003,
    0.005,
    0.01,
    0.02,
    0.05,
)
#: Analysis order is ascending for plotting; the runner keeps its own order
#: because seeds are derived from position in that tuple.


def denovo_arm_id_for(selection_coefficient: float) -> str:
    return f"denovo_s{selection_coefficient:g}".replace(".", "p")


#: (origin, s, arm_id) for every selected arm in the study.
ARM_SPECS: tuple[tuple[str, float, str], ...] = tuple(
    [("introgressed", s, arm_id_for(s)) for s in SELECTION_COEFFICIENTS]
    + [("de_novo", s, denovo_arm_id_for(s)) for s in DENOVO_COEFFICIENTS]
)

ALPHA = 0.05
MIN_CLASS_PAIRS = 10
CLASSES = ("unconditional", "hom_carrier")


@dataclass(frozen=True)
class Replicate:
    positions: np.ndarray
    tmrca: np.ndarray
    genotypes: np.ndarray | None
    replicate_index: int
    selection_coefficient: float
    sample_af: float
    hom_carrier_diploids: int
    elapsed_seconds: float


def _load_directory(pattern: str) -> list[Replicate]:
    records: list[Replicate] = []
    for path in sorted(glob.glob(pattern)):
        directory = Path(path)
        matrix = directory / "tmrca_matrix.npz"
        endpoint = directory / "endpoint.json"
        if not (matrix.is_file() and endpoint.is_file()):
            continue
        record = json.loads(endpoint.read_text(encoding="utf-8"))
        if record.get("status") != "completed":
            continue
        with np.load(matrix) as data:
            positions = np.asarray(data["positions"], dtype=float)
            tmrca = np.asarray(data["tmrca_generations"], dtype=np.float32)
        genotype_path = directory / "focal_genotypes.npy"
        genotypes = np.load(genotype_path) if genotype_path.is_file() else None
        frequency = record.get("final_allele_frequency", {})
        records.append(
            Replicate(
                positions=positions,
                tmrca=tmrca,
                genotypes=genotypes,
                replicate_index=int(record["replicate_index"]),
                selection_coefficient=float(record.get("selection_coefficient", 0.0)),
                sample_af=float(frequency.get("sample_af", 0.0)),
                hom_carrier_diploids=int(frequency.get("sample_hom_carrier_diploids", 0)),
                elapsed_seconds=float(record.get("elapsed_seconds", float("nan"))),
            )
        )
    return records


def load_mode(study_root: str | Path, selection_coefficient: float | None) -> list[Replicate]:
    """Load one introgressed arm, or the shared neutral arm when ``s`` is None."""
    root = Path(study_root)
    if selection_coefficient is None:
        return _load_directory(str(root / "neutral" / "replicates" / "*"))
    return load_arm(study_root, arm_id_for(selection_coefficient))


def load_arm(study_root: str | Path, arm_id: str) -> list[Replicate]:
    root = Path(study_root)
    return _load_directory(str(root / arm_id / "selected" / "replicates" / "*"))


def _class_mask(replicate: Replicate, name: str) -> np.ndarray:
    if name == "unconditional":
        return np.ones(replicate.tmrca.shape[1], dtype=bool)
    if replicate.genotypes is None:
        return np.zeros(replicate.tmrca.shape[1], dtype=bool)
    return replicate.genotypes == 2


def _tail_p(observed: float, null: np.ndarray) -> float:
    """Upper-tail p-value with the conventional +1 correction."""
    if null.size == 0:
        return float("nan")
    return float((np.sum(null >= observed) + 1) / (null.size + 1))


def focal_table(
    study_root: str | Path,
    *,
    thresholds_years: Sequence[float] = THRESHOLDS_YEARS,
    generation_time: float = GENERATION_TIME_YEARS,
    seed: int = 20261102,
) -> pd.DataFrame:
    """Per replicate, per class, per threshold: statistic, p-value, percentile."""
    neutral = load_mode(study_root, None)
    if not neutral:
        raise RuntimeError(f"no neutral replicates under {study_root}")
    rng = np.random.default_rng(seed)
    null_cache: dict[tuple[float, int], np.ndarray] = {}

    rows: list[dict[str, Any]] = []
    for origin, coefficient, arm in ARM_SPECS:
        selected = load_arm(study_root, arm)
        for replicate in selected:
            focal = int(np.argmin(np.abs(replicate.positions - FOCAL_POSITION_BP)))
            for name in CLASSES:
                mask = _class_mask(replicate, name)
                size = int(mask.sum())
                if size < MIN_CLASS_PAIRS:
                    continue
                for years in thresholds_years:
                    generations = float(years) / generation_time
                    value = float(
                        np.mean(replicate.tmrca[focal, mask] < generations)
                    )
                    key = (float(years), size)
                    if key not in null_cache:
                        null_cache[key] = build_matched_null(
                            neutral, generations, size, rng=rng
                        )
                    null = null_cache[key]
                    rows.append(
                        {
                            "origin": origin,
                            "arm_id": arm,
                            "selection_coefficient": coefficient,
                            "replicate_index": replicate.replicate_index,
                            "sample_af": replicate.sample_af,
                            "hom_carrier_diploids": replicate.hom_carrier_diploids,
                            "genotype_class": name,
                            "n_pairs": size,
                            "threshold_years": float(years),
                            "statistic": value,
                            "null_mean": float(null.mean()) if null.size else np.nan,
                            "p_value": _tail_p(value, null),
                            "percentile": float(np.mean(null < value))
                            + 0.5 * float(np.mean(null == value)),
                        }
                    )

    # The neutral arm scored against its own null, which calibrates the test:
    # a well-behaved p-value here is uniform, and power at alpha is the false
    # positive rate.
    for replicate in neutral:
        focal = int(np.argmin(np.abs(replicate.positions - FOCAL_POSITION_BP)))
        size = replicate.tmrca.shape[1]
        for years in thresholds_years:
            generations = float(years) / generation_time
            value = float(np.mean(replicate.tmrca[focal] < generations))
            key = (float(years), size)
            if key not in null_cache:
                null_cache[key] = build_matched_null(
                    neutral, generations, size, rng=rng
                )
            null = null_cache[key]
            rows.append(
                {
                    "origin": "neutral",
                    "arm_id": "neutral",
                    "selection_coefficient": 0.0,
                    "replicate_index": replicate.replicate_index,
                    "sample_af": 0.0,
                    "hom_carrier_diploids": 0,
                    "genotype_class": "neutral",
                    "n_pairs": size,
                    "threshold_years": float(years),
                    "statistic": value,
                    "null_mean": float(null.mean()) if null.size else np.nan,
                    "p_value": _tail_p(value, null),
                    "percentile": float(np.mean(null < value))
                    + 0.5 * float(np.mean(null == value)),
                }
            )
    return pd.DataFrame(rows)


def summarise(focal: pd.DataFrame, *, alpha: float = ALPHA) -> pd.DataFrame:
    """Power, mean AUC and mean statistic per arm, class and threshold."""
    return (
        focal.groupby(["origin", "selection_coefficient", "genotype_class", "threshold_years"])
        .agg(
            n_replicates=("p_value", "size"),
            mean_sample_af=("sample_af", "mean"),
            mean_statistic=("statistic", "mean"),
            mean_null=("null_mean", "mean"),
            auc=("percentile", "mean"),
            power=("p_value", lambda s: float(np.mean(s < alpha))),
            median_p=("p_value", "median"),
        )
        .reset_index()
    )


def scan_table(
    study_root: str | Path,
    *,
    thresholds_years: Sequence[float] = (20_000.0, 30_000.0, 50_000.0),
    generation_time: float = GENERATION_TIME_YEARS,
    seed: int = 20261103,
) -> pd.DataFrame:
    """The spatial profile: mean statistic and AUC per position and threshold."""
    neutral = load_mode(study_root, None)
    rng = np.random.default_rng(seed)
    positions = neutral[0].positions
    rows = []

    for origin, coefficient, arm in ARM_SPECS:
        selected = load_arm(study_root, arm)
        if not selected:
            continue
        for name in CLASSES:
            for years in thresholds_years:
                generations = float(years) / generation_time
                statistics = []
                percentiles = []
                for replicate in selected:
                    mask = _class_mask(replicate, name)
                    size = int(mask.sum())
                    if size < MIN_CLASS_PAIRS:
                        continue
                    values = (replicate.tmrca[:, mask] < generations).mean(axis=1)
                    null = build_matched_null(neutral, generations, size, rng=rng)
                    statistics.append(values)
                    sorted_null = null
                    left = np.searchsorted(sorted_null, values, side="left")
                    right = np.searchsorted(sorted_null, values, side="right")
                    percentiles.append(
                        (left + 0.5 * (right - left)) / max(sorted_null.size, 1)
                    )
                if not statistics:
                    continue
                rows.append(
                    pd.DataFrame(
                        {
                            "position_0based": positions,
                            "origin": origin,
                            "selection_coefficient": coefficient,
                            "genotype_class": name,
                            "threshold_years": float(years),
                            "n_replicates": len(statistics),
                            "mean_statistic": np.mean(statistics, axis=0),
                            "auc_vs_neutral": np.mean(percentiles, axis=0),
                        }
                    )
                )

    for years in thresholds_years:
        generations = float(years) / generation_time
        values = np.vstack([(r.tmrca < generations).mean(axis=1) for r in neutral])
        rows.append(
            pd.DataFrame(
                {
                    "position_0based": positions,
                    "origin": "neutral",
                    "selection_coefficient": 0.0,
                    "genotype_class": "neutral",
                    "threshold_years": float(years),
                    "n_replicates": len(neutral),
                    "mean_statistic": values.mean(axis=0),
                    "auc_vs_neutral": 0.5,
                }
            )
        )
    return pd.concat(rows, ignore_index=True)
