"""Analysis for run9.

Reuses run8's loading and class machinery unchanged; what is new is how results
are summarised.

run8 reported one AF-median replicate per arm, which was its weakest part -- a
single draw says little, and picking it by frequency conflates "is this arm
detectable" with "was this replicate lucky". run9 replaces that with:

* **quantile curves** -- the p-value at the 10th, 25th, 50th, 75th and 90th
  percentile of the statistic. Reporting the whole curve answers "is the median
  simply weak?" without cherry-picking, and note that *"the top 10% are
  detectable"* is exactly the statement *power >= 0.10*, which the power column
  already carries.
* **allele-frequency stratification** -- pooling every arm gives a wide AF range,
  and since only about ``AF^2`` of random pairs are carrier-carrier, AF is the
  quantity the selection coefficient mostly acts *through*. Binning on it
  estimates the portable answer: the frequency above which unconditional
  detection starts to work.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd

from .run8_analysis import (  # noqa: F401  (re-exported: run9 changes neither)
    CLASSES,
    MIN_CLASS_PAIRS,
    Replicate,
    _decode_profile,
    _truth_profile,
    load_arm,
)
from .run9_config import (
    ALPHA,
    FOCAL_POSITION_BP,
    QUANTILES,
    SCENARIOS,
    THRESHOLDS_YEARS,
)

__all__ = [
    "af_stratified_table",
    "exceedance_table",
    "focal_table",
    "quantile_table",
    "scan_table",
]

#: Bins wide enough to hold tens of replicates each once all arms are pooled.
AF_EDGES = (0.0, 0.05, 0.15, 0.30, 0.45, 0.60, 0.75, 0.90, 1.01)


def _arms() -> list[tuple[str, str | None]]:
    return [("neutral", None)] + [(s.scenario_id, s.scenario_id) for s in SCENARIOS]


def focal_table(
    study_root: str | Path,
    *,
    thresholds_years: Sequence[float] = THRESHOLDS_YEARS,
) -> pd.DataFrame:
    """Per replicate, per class, per threshold, per method: the focal statistic."""
    rows: list[dict[str, Any]] = []
    for arm_label, scenario_id in _arms():
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


def _null_for(neutral: pd.DataFrame, years: float, method: str) -> np.ndarray:
    return neutral[
        (neutral["threshold_years"] == years)
        & (neutral["method"] == method)
        & (neutral["genotype_class"] == "overall")
    ]["statistic"].to_numpy(dtype=float)


def _p_values(values: np.ndarray, null: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    exceed = np.array([int(np.sum(null >= v)) for v in values])
    return exceed, (exceed + 1.0) / (null.size + 1.0)


def exceedance_table(focal: pd.DataFrame, *, alpha: float = ALPHA) -> pd.DataFrame:
    """Score every selected replicate against the unconditioned neutral null."""
    rows: list[dict[str, Any]] = []
    neutral = focal[focal["arm"] == "neutral"]
    for (years, method), block in focal.groupby(["threshold_years", "method"]):
        null = _null_for(neutral, years, method)
        if null.size == 0:
            continue
        for arm, arm_block in block.groupby("arm"):
            if arm == "neutral":
                continue
            for name, class_block in arm_block.groupby("genotype_class"):
                values = class_block["statistic"].to_numpy(dtype=float)
                exceed, p = _p_values(values, null)
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
                        "mean_af": float(class_block["sample_af"].mean()),
                        "mean_statistic": float(values.mean()),
                        "null_mean": float(null.mean()),
                        "auc": float(percentile.mean()),
                        "power_rule": float(np.mean(exceed <= 1)),
                        "power_p_corrected": float(np.mean(p <= alpha)),
                        "median_p": float(np.median(p)),
                    }
                )
    return pd.DataFrame(rows)


def quantile_table(
    focal: pd.DataFrame, *, quantiles: Sequence[float] = QUANTILES
) -> pd.DataFrame:
    """p-value at each quantile of the statistic, replacing run8's single replicate.

    The quantile is taken over the *statistic* and the p-value computed for that
    value, rather than taking a quantile of the p-values: the two coincide because
    the p-value is monotone decreasing in the statistic, but this way the reported
    statistic and p-value describe the same hypothetical replicate.
    """
    rows: list[dict[str, Any]] = []
    neutral = focal[focal["arm"] == "neutral"]
    for (arm, method, name, years), block in focal.groupby(
        ["arm", "method", "genotype_class", "threshold_years"]
    ):
        if arm == "neutral":
            continue
        null = _null_for(neutral, years, method)
        if null.size == 0:
            continue
        values = block["statistic"].to_numpy(dtype=float)
        for q in quantiles:
            value = float(np.quantile(values, q))
            exceed = int(np.sum(null >= value))
            rows.append(
                {
                    "arm": arm,
                    "method": method,
                    "genotype_class": name,
                    "threshold_years": float(years),
                    "quantile": float(q),
                    "statistic": value,
                    "null_mean": float(null.mean()),
                    "n_neutral_exceeding": exceed,
                    "p_corrected": (exceed + 1.0) / (null.size + 1.0),
                    "significant_rule": bool(exceed <= 1),
                }
            )
    return pd.DataFrame(rows)


def af_stratified_table(focal: pd.DataFrame) -> pd.DataFrame:
    """Power and AUC by final allele frequency, pooled across every arm.

    Selection coefficient acts on detection almost entirely through the frequency
    it produces, so pooling arms and binning on AF gives a relationship that
    transfers to any scenario reaching that frequency, rather than one tied to a
    particular ``s``.
    """
    rows: list[dict[str, Any]] = []
    neutral = focal[focal["arm"] == "neutral"]
    selected = focal[focal["arm"] != "neutral"].copy()
    selected["af_bin"] = pd.cut(selected["sample_af"], AF_EDGES, right=False)

    for (name, method, years, af_bin), block in selected.groupby(
        ["genotype_class", "method", "threshold_years", "af_bin"], observed=True
    ):
        null = _null_for(neutral, years, method)
        if null.size == 0 or block.empty:
            continue
        values = block["statistic"].to_numpy(dtype=float)
        exceed, p = _p_values(values, null)
        percentile = np.array(
            [float(np.mean(null < v)) + 0.5 * float(np.mean(null == v)) for v in values]
        )
        rows.append(
            {
                "genotype_class": name,
                "method": method,
                "threshold_years": float(years),
                "af_bin": str(af_bin),
                "af_low": float(af_bin.left),
                "af_mid": float(block["sample_af"].mean()),
                "n_replicates": int(values.size),
                "mean_statistic": float(values.mean()),
                "null_mean": float(null.mean()),
                "auc": float(percentile.mean()),
                "power_rule": float(np.mean(exceed <= 1)),
            }
        )
    return pd.DataFrame(rows)


def scan_table(
    study_root: str | Path,
    *,
    thresholds_years: Sequence[float] = THRESHOLDS_YEARS,
) -> pd.DataFrame:
    """Spatial profile of the mean statistic along the contig."""
    rows = []
    for arm_label, scenario_id in _arms():
        replicates = load_arm(study_root, scenario_id)
        if not replicates:
            continue
        positions = replicates[0].positions
        for years in thresholds_years:
            for name in CLASSES:
                profiles = [
                    _truth_profile(r, years, r.class_mask(name))
                    for r in replicates
                    if int(r.class_mask(name).sum()) >= MIN_CLASS_PAIRS
                ]
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
