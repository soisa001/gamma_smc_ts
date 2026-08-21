"""Analysis for run8's three checks, on tree truth and Gamma-SMC decode alike.

The three questions, in the order they were posed:

1. **Unconditional.** Without conditioning on carrier status at all, does
   ``P(TMRCA < x)`` spike locally at the selected site? This is the only version
   the empirical pipeline can currently run, because hmmix enters as a post-hoc
   overlap rather than as an input.
2. **Conditioned.** Stratifying pairs by carrier status -- carrier-carrier,
   carrier-noncarrier, noncarrier-noncarrier -- is the spike stronger? This is an
   upper bound on what perfect ancestry calls would buy, not a proposed method.
3. **Calibration.** p-value distributions and AUC against a null that is
   unconditioned and built only from the overall statistic.

Two methods are carried side by side throughout:

``truth``
    The genealogical answer, computed on the same 10,000 haplotype pairs the
    decoder drew.
``decode``
    Gamma-SMC's ``frac_recent_<T>``, which with ``--recent_call mean`` is the
    proportion of posterior means below ``T`` -- the paper's statistic. Its
    ``mean_p_lt_<T>`` is carried as a secondary.

Only ``truth`` can be stratified by carrier class here: the decoder emits
per-position aggregates over its whole pair set, so a per-class decode would need
three separate runs against three pair files. Check 2 is therefore reported on
truth, which is sufficient for an upper bound.

The significance rule follows the study design: with 100 neutral replicates, a
selected replicate is called significant when **at most one** neutral replicate
matches or exceeds it. Both the raw exceedance count and the conventional
``(k+1)/(n+1)`` p-value are reported, because the two differ at exactly this
boundary -- ``k = 1`` gives ``p = 0.0198`` under the correction and ``p = 0.01``
without it.
"""

from __future__ import annotations

import glob
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd

from .run8_config import (
    ALPHA,
    FOCAL_POSITION_BP,
    GENERATION_TIME_YEARS,
    SCENARIOS,
    THRESHOLDS_YEARS,
)

__all__ = [
    "CLASSES",
    "Replicate",
    "exceedance_table",
    "focal_table",
    "load_arm",
    "scan_table",
]

CLASSES = ("overall", "carrier_carrier", "carrier_noncarrier", "noncarrier_noncarrier")
MIN_CLASS_PAIRS = 50


@dataclass(frozen=True)
class Replicate:
    positions: np.ndarray
    tmrca: np.ndarray  # (n_positions, n_pairs) generations
    pair_haplotypes: np.ndarray  # (n_pairs, 2) indices into the haplotype panel
    carriers: np.ndarray  # (n_haplotypes,) bool
    decoded: pd.DataFrame | None
    replicate_index: int
    scenario_id: str
    sample_af: float

    def class_mask(self, name: str) -> np.ndarray:
        if name == "overall":
            return np.ones(self.tmrca.shape[1], dtype=bool)
        left = self.carriers[self.pair_haplotypes[:, 0]]
        right = self.carriers[self.pair_haplotypes[:, 1]]
        if name == "carrier_carrier":
            return left & right
        if name == "carrier_noncarrier":
            return left ^ right
        if name == "noncarrier_noncarrier":
            return ~left & ~right
        raise ValueError(f"unknown class {name!r}")


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
            pair_haplotypes = np.asarray(data["pair_haplotypes"], dtype=np.int64)
        carriers_path = directory / "focal_haplotype_carriers.npy"
        carriers = (
            np.load(carriers_path)
            if carriers_path.is_file()
            else np.zeros(int(pair_haplotypes.max()) + 1, dtype=bool)
        )
        summary_path = directory / "gamma_smc_summary.tsv"
        decoded = pd.read_csv(summary_path, sep="\t") if summary_path.is_file() else None
        records.append(
            Replicate(
                positions=positions,
                tmrca=tmrca,
                pair_haplotypes=pair_haplotypes,
                carriers=carriers,
                decoded=decoded,
                replicate_index=int(record["replicate_index"]),
                scenario_id=record["scenario_id"],
                sample_af=float(
                    record.get("final_allele_frequency", {}).get("sample_af", 0.0)
                ),
            )
        )
    return records


def load_arm(study_root: str | Path, scenario_id: str | None) -> list[Replicate]:
    """Load one selected scenario, or the shared neutral arm when ``None``."""
    root = Path(study_root)
    if scenario_id is None:
        return _load_directory(str(root / "neutral" / "replicates" / "*"))
    return _load_directory(str(root / scenario_id / "selected" / "replicates" / "*"))


def _truth_profile(replicate: Replicate, years: float, mask: np.ndarray) -> np.ndarray:
    generations = float(years) / GENERATION_TIME_YEARS
    return (replicate.tmrca[:, mask] < generations).mean(axis=1)


def _decode_profile(replicate: Replicate, years: float, column: str) -> np.ndarray | None:
    if replicate.decoded is None:
        return None
    name = f"{column}_{int(years)}"
    if name not in replicate.decoded.columns:
        return None
    return replicate.decoded[name].to_numpy(dtype=float)


def focal_table(
    study_root: str | Path,
    *,
    thresholds_years: Sequence[float] = THRESHOLDS_YEARS,
) -> pd.DataFrame:
    """Per replicate, per class, per threshold, per method: the focal statistic."""
    rows: list[dict[str, Any]] = []
    arms: list[tuple[str, str | None]] = [("neutral", None)] + [
        (s.scenario_id, s.scenario_id) for s in SCENARIOS
    ]
    for arm_label, scenario_id in arms:
        for replicate in load_arm(study_root, scenario_id):
            focal = int(np.argmin(np.abs(replicate.positions - FOCAL_POSITION_BP)))
            for years in thresholds_years:
                for name in CLASSES:
                    mask = replicate.class_mask(name)
                    size = int(mask.sum())
                    if size < MIN_CLASS_PAIRS:
                        continue
                    rows.append(
                        {
                            "arm": arm_label,
                            "replicate_index": replicate.replicate_index,
                            "sample_af": replicate.sample_af,
                            "threshold_years": float(years),
                            "genotype_class": name,
                            "n_pairs": size,
                            "method": "truth",
                            "statistic": float(
                                _truth_profile(replicate, years, mask)[focal]
                            ),
                        }
                    )
                # The decoder aggregates over its whole pair set, so its output
                # exists only for the unstratified class.
                for method, column in (
                    ("decode", "frac_recent"),
                    ("decode_meanp", "mean_p_lt"),
                ):
                    profile = _decode_profile(replicate, years, column)
                    if profile is None:
                        continue
                    rows.append(
                        {
                            "arm": arm_label,
                            "replicate_index": replicate.replicate_index,
                            "sample_af": replicate.sample_af,
                            "threshold_years": float(years),
                            "genotype_class": "overall",
                            "n_pairs": int(replicate.pair_haplotypes.shape[0]),
                            "method": method,
                            "statistic": float(profile[focal]),
                        }
                    )
    return pd.DataFrame(rows)


def exceedance_table(focal: pd.DataFrame, *, alpha: float = ALPHA) -> pd.DataFrame:
    """Score every selected replicate against the unconditioned neutral null.

    The null is always the neutral arm's *overall* statistic under the same
    method and threshold -- never a class-matched null. That is the comparison
    the study asked for, and for the stratified classes it is deliberately
    optimistic: conditioning on carrier status shifts the statistic under
    neutrality too, so those rows measure headroom rather than calibrated power.
    """
    rows: list[dict[str, Any]] = []
    neutral = focal[focal["arm"] == "neutral"]
    for (years, method), block in focal.groupby(["threshold_years", "method"]):
        null = neutral[
            (neutral["threshold_years"] == years)
            & (neutral["method"] == method)
            & (neutral["genotype_class"] == "overall")
        ]["statistic"].to_numpy(dtype=float)
        if null.size == 0:
            continue
        for arm, arm_block in block.groupby("arm"):
            if arm == "neutral":
                continue
            for name, class_block in arm_block.groupby("genotype_class"):
                values = class_block["statistic"].to_numpy(dtype=float)
                exceed = np.array([int(np.sum(null >= v)) for v in values])
                p_corrected = (exceed + 1.0) / (null.size + 1.0)
                p_raw = exceed / float(null.size)
                percentile = np.array(
                    [
                        float(np.mean(null < v)) + 0.5 * float(np.mean(null == v))
                        for v in values
                    ]
                )
                rows.append(
                    {
                        "arm": arm,
                        "method": method,
                        "genotype_class": name,
                        "threshold_years": float(years),
                        "n_replicates": int(values.size),
                        "n_null": int(null.size),
                        "mean_statistic": float(values.mean()),
                        "null_mean": float(null.mean()),
                        "auc": float(percentile.mean()),
                        # The study rule: at most one neutral replicate may match
                        # or exceed the observation.
                        "power_rule": float(np.mean(exceed <= 1)),
                        "power_p_corrected": float(np.mean(p_corrected <= alpha)),
                        "median_exceedance": float(np.median(exceed)),
                        "median_p_corrected": float(np.median(p_corrected)),
                        "median_p_raw": float(np.median(p_raw)),
                    }
                )
    return pd.DataFrame(rows)


def scan_table(
    study_root: str | Path,
    *,
    thresholds_years: Sequence[float] = THRESHOLDS_YEARS,
) -> pd.DataFrame:
    """Spatial profile of the mean statistic along the contig, per arm and class."""
    rows = []
    arms: list[tuple[str, str | None]] = [("neutral", None)] + [
        (s.scenario_id, s.scenario_id) for s in SCENARIOS
    ]
    for arm_label, scenario_id in arms:
        replicates = load_arm(study_root, scenario_id)
        if not replicates:
            continue
        positions = replicates[0].positions
        for years in thresholds_years:
            for name in CLASSES:
                profiles = []
                for replicate in replicates:
                    mask = replicate.class_mask(name)
                    if int(mask.sum()) < MIN_CLASS_PAIRS:
                        continue
                    profiles.append(_truth_profile(replicate, years, mask))
                if profiles:
                    rows.append(
                        pd.DataFrame(
                            {
                                "position_0based": positions,
                                "arm": arm_label,
                                "threshold_years": float(years),
                                "genotype_class": name,
                                "method": "truth",
                                "n_replicates": len(profiles),
                                "mean_statistic": np.mean(profiles, axis=0),
                            }
                        )
                    )
            for method, column in (("decode", "frac_recent"), ("decode_meanp", "mean_p_lt")):
                profiles = [
                    p
                    for p in (_decode_profile(r, years, column) for r in replicates)
                    if p is not None
                ]
                if profiles:
                    rows.append(
                        pd.DataFrame(
                            {
                                "position_0based": positions,
                                "arm": arm_label,
                                "threshold_years": float(years),
                                "genotype_class": "overall",
                                "method": method,
                                "n_replicates": len(profiles),
                                "mean_statistic": np.mean(profiles, axis=0),
                            }
                        )
                    )
    return pd.concat(rows, ignore_index=True)
