"""Statistics, null distribution and p-values for the run5 study.

Everything here reads the *raw* per-pair TMRCA arrays written by
``run5_simulate`` and the raw decoder summaries written by ``run5_decode``, so
thresholds, genotype stratifications and conditioning are all post-hoc.  Nothing
requires re-simulating.

Three estimators are carried side by side:

* ``tree_truth`` -- the fraction of within-individual pairs whose true TMRCA at
  the focal base is below the threshold.
* ``gamma_smc_frac`` -- the same fraction using Gamma-SMC's *called* TMRCA.
* ``gamma_smc_meanp`` -- Gamma-SMC's mean posterior probability that the TMRCA is
  below the threshold, which uses the whole posterior instead of a point call.

Genotype stratification is available for tree truth only: the decoder summary is
already aggregated over pairs, so per-pair carrier classes cannot be recovered
from it without a separate per-pair decode.  The headline truth-versus-decoded
comparison therefore uses the ``overall`` class for both.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd

from .run5_config import (
    FOCAL_POSITION_BP,
    MODES,
    Run5Arm,
    SIGNIFICANCE_LEVEL,
    TMRCA_THRESHOLDS_YEARS,
)
from .run5_simulate import GENOTYPE_CLASSES, REPLICATE_OUTPUTS

ESTIMATORS = ("tree_truth", "gamma_smc_frac", "gamma_smc_meanp")

RESULT_OUTPUTS = {
    "replicate_status": "replicate_status.tsv",
    "final_af": "final_af.tsv",
    "statistics": "statistics.tsv",
    "null_distribution": "null_distribution.tsv",
    "pvalues": "pvalues.tsv",
    "power": "power_summary.tsv",
    "null_calibration": "null_calibration.tsv",
    "truth_vs_decoded": "truth_vs_decoded.tsv",
    "spatial_profiles": "spatial_profiles.tsv.gz",
}

_NULL_QUANTILES = (0.5, 0.9, 0.95, 0.975, 0.99)


def _replicate_directories(study_root: str | Path, arm: Run5Arm, mode: str) -> list[Path]:
    base = Path(study_root) / arm.arm_id / mode / "replicates"
    if not base.is_dir():
        return []
    return sorted(p for p in base.iterdir() if p.is_dir())


def _load_endpoint(directory: Path) -> dict[str, Any] | None:
    path = directory / REPLICATE_OUTPUTS["endpoint"]
    if not path.is_file():
        return None
    record = json.loads(path.read_text(encoding="utf-8"))
    return record if record.get("status") == "completed" else None


def _class_mask(counts: np.ndarray, genotype_class: str) -> np.ndarray:
    if genotype_class == "overall":
        return np.ones(counts.shape, dtype=bool)
    if genotype_class == "any_carrier":
        return counts > 0
    if genotype_class == "hom_carrier":
        return counts == 2
    if genotype_class == "het":
        return counts == 1
    if genotype_class == "hom_noncarrier":
        return counts == 0
    raise ValueError(f"unknown genotype class {genotype_class!r}")


# ---------------------------------------------------------------------------
# Per-replicate statistics
# ---------------------------------------------------------------------------


def replicate_statistics(
    study_root: str | Path,
    arm: Run5Arm,
    *,
    thresholds_years: Sequence[float] = TMRCA_THRESHOLDS_YEARS,
) -> pd.DataFrame:
    """Return the long-form focal statistic for every completed replicate.

    One row per (mode, replicate, estimator, genotype class, threshold).
    """
    rows: list[dict[str, Any]] = []
    for mode in MODES:
        for directory in _replicate_directories(study_root, arm, mode):
            endpoint = _load_endpoint(directory)
            if endpoint is None:
                continue
            index = int(endpoint["replicate_index"])
            present = bool(
                endpoint["final_allele_frequency"]["archaic_present_at_focal_base"]
            )
            base = {
                "arm_id": arm.arm_id,
                "mode": mode,
                "replicate_index": index,
                "unit_id": endpoint["unit_id"],
                "seed": endpoint["seed"],
                "allele_present": present,
            }

            focal = np.load(directory / REPLICATE_OUTPUTS["focal_tmrca"])
            counts = np.load(directory / REPLICATE_OUTPUTS["focal_genotypes"])
            for genotype_class in GENOTYPE_CLASSES:
                mask = _class_mask(counts, genotype_class)
                selected = focal[mask]
                for years in thresholds_years:
                    generations = float(years) / arm.generation_time
                    rows.append(
                        {
                            **base,
                            "estimator": "tree_truth",
                            "genotype_class": genotype_class,
                            "n_pairs": int(mask.sum()),
                            "threshold_years": float(years),
                            "value": (
                                float(np.mean(selected < generations))
                                if selected.size
                                else float("nan")
                            ),
                        }
                    )

            summary_path = directory / "gamma_smc_summary.tsv"
            if summary_path.is_file():
                frame = pd.read_csv(summary_path, sep="\t")
                row = frame.iloc[
                    (frame["position_0based"] - FOCAL_POSITION_BP).abs().idxmin()
                ]
                for years in thresholds_years:
                    tag = int(years)
                    for estimator, column in (
                        ("gamma_smc_frac", f"frac_recent_{tag}"),
                        ("gamma_smc_meanp", f"mean_p_lt_{tag}"),
                    ):
                        if column not in frame.columns:
                            continue
                        rows.append(
                            {
                                **base,
                                "estimator": estimator,
                                "genotype_class": "overall",
                                "n_pairs": int(row.get("n_pairs", len(counts))),
                                "threshold_years": float(years),
                                "value": float(row[column]),
                            }
                        )
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).sort_values(
        ["mode", "estimator", "genotype_class", "threshold_years", "replicate_index"],
        ignore_index=True,
    )


def replicate_status(study_root: str | Path, arm: Run5Arm) -> pd.DataFrame:
    rows = []
    for mode in MODES:
        for directory in _replicate_directories(study_root, arm, mode):
            path = directory / REPLICATE_OUTPUTS["endpoint"]
            if not path.is_file():
                continue
            record = json.loads(path.read_text(encoding="utf-8"))
            decoded = (directory / "gamma_smc_receipt.json").is_file()
            rows.append(
                {
                    "arm_id": arm.arm_id,
                    "mode": mode,
                    "unit_id": record["unit_id"],
                    "replicate_index": record["replicate_index"],
                    "seed": record["seed"],
                    "status": record["status"],
                    "elapsed_seconds": record.get("elapsed_seconds"),
                    "decoded": decoded,
                    "error": record.get("error", ""),
                }
            )
    return pd.DataFrame(rows).sort_values(["mode", "replicate_index"], ignore_index=True)


def final_allele_frequencies(study_root: str | Path, arm: Run5Arm) -> pd.DataFrame:
    rows = []
    for mode in MODES:
        for directory in _replicate_directories(study_root, arm, mode):
            endpoint = _load_endpoint(directory)
            if endpoint is None:
                continue
            final = endpoint["final_allele_frequency"]
            placement = endpoint["placement"]
            summary = endpoint["focal_tmrca_summary"]
            census = final.get("census_af")
            census = float(census) if census is not None else float("nan")
            rows.append(
                {
                    "arm_id": arm.arm_id,
                    "mode": mode,
                    "replicate_index": endpoint["replicate_index"],
                    "seed": endpoint["seed"],
                    "placement_frequency": placement.get("frequency"),
                    "placement_carrier_genomes": placement.get("carrier_genomes"),
                    "census_af": census,
                    "census_state": final.get("census_state"),
                    "sample_af": final.get("sample_af"),
                    "sample_carrier_diploids": final.get("sample_carrier_diploids"),
                    "sample_hom_carrier_diploids": final.get("sample_hom_carrier_diploids"),
                    "allele_present": final.get("archaic_present_at_focal_base"),
                    "focal_mean_tmrca_generations": summary.get("mean_generations"),
                    "focal_median_tmrca_generations": summary.get("median_generations"),
                }
            )
    return pd.DataFrame(rows).sort_values(
        ["mode", "replicate_index"], ignore_index=True
    )


# ---------------------------------------------------------------------------
# Null, p-values and power
# ---------------------------------------------------------------------------


def upper_tail_pvalue(observed: float, null_values: np.ndarray) -> float:
    if not np.isfinite(observed) or null_values.size == 0:
        return float("nan")
    return (1.0 + int(np.sum(null_values >= observed))) / (1.0 + null_values.size)


def _subset(statistics: pd.DataFrame, mode: str, estimator: str, genotype_class: str):
    return statistics[
        (statistics["mode"] == mode)
        & (statistics["estimator"] == estimator)
        & (statistics["genotype_class"] == genotype_class)
    ]


def null_distribution(statistics: pd.DataFrame, arm: Run5Arm) -> pd.DataFrame:
    rows = []
    for estimator in statistics["estimator"].unique():
        for genotype_class in statistics["genotype_class"].unique():
            frame = _subset(statistics, "neutral", estimator, genotype_class)
            if frame.empty:
                continue
            for years, group in frame.groupby("threshold_years"):
                values = group["value"].to_numpy(dtype=float)
                values = values[np.isfinite(values)]
                row: dict[str, Any] = {
                    "arm_id": arm.arm_id,
                    "estimator": estimator,
                    "genotype_class": genotype_class,
                    "threshold_years": float(years),
                    "n_null": int(values.size),
                    "mean": float(values.mean()) if values.size else float("nan"),
                    "sd": float(values.std(ddof=1)) if values.size > 1 else float("nan"),
                }
                for q in _NULL_QUANTILES:
                    row[f"q{int(q * 1000):03d}"] = (
                        float(np.quantile(values, q)) if values.size else float("nan")
                    )
                rows.append(row)
    return pd.DataFrame(rows)


def pvalues(
    statistics: pd.DataFrame, arm: Run5Arm, *, condition_on_presence: bool = False
) -> pd.DataFrame:
    """Return one p-value per selected replicate, estimator, class and threshold.

    ``condition_on_presence`` restricts *both* arms to replicates that carry the
    focal allele at all.  In run5 that is a substantive restriction: archaic
    ancestry is absent at any single base in most replicates, and a replicate
    with no allele is neutral by construction in both modes.
    """
    rows = []
    for estimator in statistics["estimator"].unique():
        for genotype_class in statistics["genotype_class"].unique():
            selected = _subset(statistics, "selected", estimator, genotype_class)
            neutral = _subset(statistics, "neutral", estimator, genotype_class)
            if selected.empty or neutral.empty:
                continue
            if condition_on_presence:
                selected = selected[selected["allele_present"]]
                neutral = neutral[neutral["allele_present"]]
            for years in sorted(selected["threshold_years"].unique()):
                null_values = neutral.loc[
                    neutral["threshold_years"] == years, "value"
                ].to_numpy(dtype=float)
                null_values = null_values[np.isfinite(null_values)]
                for record in selected[selected["threshold_years"] == years].itertuples(
                    index=False
                ):
                    p = upper_tail_pvalue(float(record.value), null_values)
                    rows.append(
                        {
                            "arm_id": arm.arm_id,
                            "estimator": estimator,
                            "genotype_class": genotype_class,
                            "conditioned_on_presence": condition_on_presence,
                            "replicate_index": int(record.replicate_index),
                            "allele_present": bool(record.allele_present),
                            "threshold_years": float(years),
                            "observed": float(record.value),
                            "n_null": int(null_values.size),
                            "p_value": p,
                            "significant": bool(
                                np.isfinite(p) and p <= SIGNIFICANCE_LEVEL
                            ),
                        }
                    )
    return pd.DataFrame(rows)


def power_summary(pvalue_table: pd.DataFrame, arm: Run5Arm) -> pd.DataFrame:
    if pvalue_table.empty:
        return pd.DataFrame()
    grouped = pvalue_table.groupby(
        ["estimator", "genotype_class", "conditioned_on_presence", "threshold_years"]
    )
    rows = []
    for (estimator, genotype_class, conditioned, years), group in grouped:
        finite = group[np.isfinite(group["p_value"])]
        rows.append(
            {
                "arm_id": arm.arm_id,
                "estimator": estimator,
                "genotype_class": genotype_class,
                "conditioned_on_presence": conditioned,
                "threshold_years": float(years),
                "n_selected": int(len(finite)),
                "n_significant": int(finite["significant"].sum()),
                "power": float(finite["significant"].mean()) if len(finite) else float("nan"),
                "median_p_value": float(finite["p_value"].median()) if len(finite) else float("nan"),
                "median_observed": float(finite["observed"].median()) if len(finite) else float("nan"),
                "n_null": int(finite["n_null"].iloc[0]) if len(finite) else 0,
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["estimator", "genotype_class", "conditioned_on_presence", "threshold_years"],
        ignore_index=True,
    )


def null_calibration(statistics: pd.DataFrame, arm: Run5Arm) -> pd.DataFrame:
    """Leave-one-out false-positive rate of the null against itself."""
    rows = []
    for estimator in statistics["estimator"].unique():
        frame = _subset(statistics, "neutral", estimator, "overall")
        if frame.empty:
            continue
        for years, group in frame.groupby("threshold_years"):
            values = group["value"].to_numpy(dtype=float)
            values = values[np.isfinite(values)]
            if values.size < 2:
                continue
            loo = np.array(
                [
                    upper_tail_pvalue(values[i], np.delete(values, i))
                    for i in range(values.size)
                ]
            )
            rows.append(
                {
                    "arm_id": arm.arm_id,
                    "estimator": estimator,
                    "threshold_years": float(years),
                    "n_null": int(values.size),
                    "false_positive_rate": float(np.mean(loo <= SIGNIFICANCE_LEVEL)),
                    "median_p_value": float(np.median(loo)),
                }
            )
    return pd.DataFrame(rows)


def truth_vs_decoded(statistics: pd.DataFrame, arm: Run5Arm) -> pd.DataFrame:
    """Per-replicate agreement between the truth statistic and each decoder one."""
    overall = statistics[statistics["genotype_class"] == "overall"]
    wide = overall.pivot_table(
        index=["mode", "replicate_index", "threshold_years"],
        columns="estimator",
        values="value",
    ).reset_index()
    rows = []
    for estimator in ("gamma_smc_frac", "gamma_smc_meanp"):
        if estimator not in wide.columns or "tree_truth" not in wide.columns:
            continue
        for (mode, years), group in wide.groupby(["mode", "threshold_years"]):
            a = group["tree_truth"].to_numpy(dtype=float)
            b = group[estimator].to_numpy(dtype=float)
            ok = np.isfinite(a) & np.isfinite(b)
            if ok.sum() < 3:
                continue
            rows.append(
                {
                    "arm_id": arm.arm_id,
                    "mode": mode,
                    "estimator": estimator,
                    "threshold_years": float(years),
                    "n": int(ok.sum()),
                    "truth_mean": float(a[ok].mean()),
                    "decoded_mean": float(b[ok].mean()),
                    "bias": float(b[ok].mean() - a[ok].mean()),
                    "rmse": float(np.sqrt(np.mean((b[ok] - a[ok]) ** 2))),
                    "pearson_r": (
                        float(np.corrcoef(a[ok], b[ok])[0, 1])
                        if a[ok].std() > 0 and b[ok].std() > 0
                        else float("nan")
                    ),
                }
            )
    return pd.DataFrame(rows)


def aggregate_spatial_profiles(
    study_root: str | Path,
    arm: Run5Arm,
    *,
    thresholds_years: Sequence[float] = TMRCA_THRESHOLDS_YEARS,
) -> pd.DataFrame:
    """Per-position truth profiles, summarised across replicates, per mode."""
    frames = []
    for mode in MODES:
        stack, positions = [], None
        for directory in _replicate_directories(study_root, arm, mode):
            path = directory / REPLICATE_OUTPUTS["tmrca_matrix"]
            if not path.is_file():
                continue
            with np.load(path) as data:
                positions = data["positions"]
                stack.append(data["tmrca_generations"].astype(np.float32))
        if not stack or positions is None:
            continue
        for years in thresholds_years:
            generations = float(years) / arm.generation_time
            values = np.vstack(
                [np.mean(m < generations, axis=1)[None, :] for m in stack]
            )
            frames.append(
                pd.DataFrame(
                    {
                        "arm_id": arm.arm_id,
                        "mode": mode,
                        "threshold_years": float(years),
                        "position_0based": positions,
                        "n_replicates": values.shape[0],
                        "mean": values.mean(axis=0),
                        "median": np.median(values, axis=0),
                        "q025": np.quantile(values, 0.025, axis=0),
                        "q975": np.quantile(values, 0.975, axis=0),
                    }
                )
            )
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------


def analyse_arm(study_root: str | Path, arm: Run5Arm) -> dict[str, Path]:
    results = Path(study_root) / arm.arm_id / "results"
    results.mkdir(parents=True, exist_ok=True)

    statistics = replicate_statistics(study_root, arm)
    if statistics.empty:
        raise RuntimeError(f"no completed replicates found for arm {arm.arm_id}")

    unconditional = pvalues(statistics, arm, condition_on_presence=False)
    conditional = pvalues(statistics, arm, condition_on_presence=True)
    both = pd.concat([unconditional, conditional], ignore_index=True)

    tables = {
        "replicate_status": replicate_status(study_root, arm),
        "final_af": final_allele_frequencies(study_root, arm),
        "statistics": statistics,
        "null_distribution": null_distribution(statistics, arm),
        "pvalues": both,
        "power": power_summary(both, arm),
        "null_calibration": null_calibration(statistics, arm),
        "truth_vs_decoded": truth_vs_decoded(statistics, arm),
        "spatial_profiles": aggregate_spatial_profiles(study_root, arm),
    }
    written: dict[str, Path] = {}
    for key, frame in tables.items():
        path = results / RESULT_OUTPUTS[key]
        frame.to_csv(path, sep="\t", index=False)
        written[key] = path
    return written


def cross_arm_power(study_root: str | Path, arms: Iterable[Run5Arm]) -> pd.DataFrame:
    frames = []
    for arm in arms:
        path = Path(study_root) / arm.arm_id / "results" / RESULT_OUTPUTS["power"]
        if path.is_file():
            frame = pd.read_csv(path, sep="\t")
            frame["arm_label"] = arm.label
            frames.append(frame)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


__all__ = [
    "ESTIMATORS",
    "RESULT_OUTPUTS",
    "aggregate_spatial_profiles",
    "analyse_arm",
    "cross_arm_power",
    "final_allele_frequencies",
    "null_calibration",
    "null_distribution",
    "power_summary",
    "pvalues",
    "replicate_statistics",
    "replicate_status",
    "truth_vs_decoded",
    "upper_tail_pvalue",
]
