"""Threshold-free power diagnostics for run5: AUC, and p-value distributions.

``power`` in the existing tables is a *per-replicate detection rate*: each
selected replicate is tested individually against the pooled neutral
distribution, and power is the fraction of them reaching ``p <= 0.05``.  That
answers "if I ran a scan on one such locus, how often would I call it", but it
depends on a threshold and on the size of the null.

Two complements, both computed post hoc from the stored statistics:

* **AUC** -- the probability that a randomly chosen selected replicate scores
  above a randomly chosen neutral one, with ties split.  Threshold-free, and
  unaffected by the p-value floor, so it separates "how well does the statistic
  discriminate" from "how often does it clear 0.05 with 100 null replicates".
  Averaging each selected replicate's percentile against the null gives exactly
  the AUC, which is what lets AUC be related to per-replicate quantities such as
  the final allele frequency.
* **p-value distributions** -- selected against the leave-one-out null. Under the
  null the p-values should be near-uniform; under selection they pile at the
  floor. Plotting both together shows calibration and signal in one panel.

Ties matter here: with 100 diploid pairs the statistic only takes values in
steps of 0.01, so the 0.5-weight tie convention is used throughout rather than a
strict inequality.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .run5_analysis import RESULT_OUTPUTS, upper_tail_pvalue

#: The statistic the test uses. Carrier stratification exists for tree truth
#: only, and the null carries no allele, so "overall" is the comparable class.
TEST_CLASS = "overall"
from .run5_config import ARMS, SIGNIFICANCE_LEVEL, Run5Arm, TMRCA_THRESHOLDS_YEARS

_COLOURS = {
    "tree_truth": "#111111",
    "gamma_smc_frac": "#b2182b",
    "gamma_smc_meanp": "#ef8a62",
}
POWER_OUTPUTS = {
    "auc": "auc_summary.tsv",
    "percentiles": "replicate_percentiles.tsv",
}


def _save(fig, stem: Path) -> dict[str, Path]:
    stem.parent.mkdir(parents=True, exist_ok=True)
    written = {}
    for suffix in ("png", "pdf"):
        path = stem.with_suffix(f".{suffix}")
        fig.savefig(path, dpi=170, bbox_inches="tight")
        written[suffix] = path
    plt.close(fig)
    return written


def _load(results: Path, key: str) -> pd.DataFrame:
    path = results / RESULT_OUTPUTS[key]
    if not path.is_file():
        raise FileNotFoundError(f"missing result table: {path}")
    return pd.read_csv(path, sep="\t")


def _split(statistics: pd.DataFrame, estimator: str, years: float):
    frame = statistics[
        (statistics["estimator"] == estimator)
        & (statistics["genotype_class"] == TEST_CLASS)
        & (statistics["threshold_years"] == years)
    ]
    selected = frame[frame["mode"] == "selected"]
    neutral = frame[frame["mode"] == "neutral"]["value"].to_numpy(dtype=float)
    return selected, neutral[np.isfinite(neutral)]


def percentile_against_null(value: float, null_values: np.ndarray) -> float:
    """Fraction of the null below ``value``, splitting ties.

    Averaged over selected replicates this is exactly the AUC, so it is the
    per-replicate quantity that AUC decomposes into.
    """
    if not np.isfinite(value) or null_values.size == 0:
        return float("nan")
    below = float(np.sum(null_values < value))
    tied = float(np.sum(null_values == value))
    return (below + 0.5 * tied) / null_values.size


def replicate_percentiles(
    study_root: str | Path,
    arm: Run5Arm,
    *,
    thresholds_years: Sequence[float] = TMRCA_THRESHOLDS_YEARS,
) -> pd.DataFrame:
    """Per selected replicate: percentile against the null, p-value, final AF."""
    results = Path(study_root) / arm.arm_id / "results"
    statistics = _load(results, "statistics")
    final_af = _load(results, "final_af")
    af_by_index = (
        final_af[final_af["mode"] == "selected"]
        .set_index("replicate_index")["census_af"]
        .to_dict()
    )

    rows = []
    for estimator in sorted(statistics["estimator"].unique()):
        for years in thresholds_years:
            selected, null_values = _split(statistics, estimator, float(years))
            if selected.empty or null_values.size == 0:
                continue
            for record in selected.itertuples(index=False):
                value = float(record.value)
                rows.append(
                    {
                        "arm_id": arm.arm_id,
                        "estimator": estimator,
                        "threshold_years": float(years),
                        "replicate_index": int(record.replicate_index),
                        "observed": value,
                        "percentile": percentile_against_null(value, null_values),
                        "p_value": upper_tail_pvalue(value, null_values),
                        "census_af": af_by_index.get(int(record.replicate_index), np.nan),
                        "n_null": int(null_values.size),
                    }
                )
    return pd.DataFrame(rows)


def null_percentiles(
    study_root: str | Path,
    arm: Run5Arm,
    *,
    thresholds_years: Sequence[float] = TMRCA_THRESHOLDS_YEARS,
) -> pd.DataFrame:
    """Leave-one-out p-values for the neutral replicates, for calibration."""
    results = Path(study_root) / arm.arm_id / "results"
    statistics = _load(results, "statistics")
    rows = []
    for estimator in sorted(statistics["estimator"].unique()):
        for years in thresholds_years:
            _, null_values = _split(statistics, estimator, float(years))
            if null_values.size < 2:
                continue
            for index in range(null_values.size):
                others = np.delete(null_values, index)
                rows.append(
                    {
                        "arm_id": arm.arm_id,
                        "estimator": estimator,
                        "threshold_years": float(years),
                        "p_value": upper_tail_pvalue(null_values[index], others),
                    }
                )
    return pd.DataFrame(rows)


def auc_summary(percentiles: pd.DataFrame, arm: Run5Arm) -> pd.DataFrame:
    """AUC per estimator and threshold, plus AUC within final-AF bins."""
    rows = []
    for (estimator, years), group in percentiles.groupby(
        ["estimator", "threshold_years"]
    ):
        finite = group[np.isfinite(group["percentile"])]
        if finite.empty:
            continue
        row: dict[str, Any] = {
            "arm_id": arm.arm_id,
            "estimator": estimator,
            "threshold_years": float(years),
            "n_selected": int(len(finite)),
            "n_null": int(finite["n_null"].iloc[0]),
            "auc": float(finite["percentile"].mean()),
            "power_at_alpha": float((finite["p_value"] <= SIGNIFICANCE_LEVEL).mean()),
            "median_p_value": float(finite["p_value"].median()),
        }
        af = finite["census_af"]
        if af.notna().any():
            for label, mask in (
                ("af_lt_0p4", af < 0.4),
                ("af_0p4_0p7", (af >= 0.4) & (af < 0.7)),
                ("af_ge_0p7", af >= 0.7),
            ):
                subset = finite[mask.fillna(False)]
                row[f"auc_{label}"] = (
                    float(subset["percentile"].mean()) if len(subset) else float("nan")
                )
                row[f"n_{label}"] = int(len(subset))
        rows.append(row)
    return pd.DataFrame(rows).sort_values(
        ["estimator", "threshold_years"], ignore_index=True
    )


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------


def plot_pvalue_distribution(
    arm: Run5Arm,
    percentiles: pd.DataFrame,
    null_p: pd.DataFrame,
    figures: Path,
    *,
    estimator: str = "tree_truth",
) -> dict[str, Path]:
    """Selected p-values against the leave-one-out null, one panel per cutoff."""
    thresholds = [
        t for t in TMRCA_THRESHOLDS_YEARS if t in set(percentiles["threshold_years"])
    ]
    columns = 4
    rows = int(np.ceil(len(thresholds) / columns))
    fig, axes = plt.subplots(rows, columns, figsize=(3.9 * columns, 3.0 * rows))
    axes = np.atleast_1d(axes).ravel()
    bins = np.linspace(0, 1, 21)
    for ax, years in zip(axes, thresholds):
        sel = percentiles[
            (percentiles["estimator"] == estimator)
            & (percentiles["threshold_years"] == years)
        ]["p_value"].to_numpy(dtype=float)
        nul = null_p[
            (null_p["estimator"] == estimator) & (null_p["threshold_years"] == years)
        ]["p_value"].to_numpy(dtype=float)
        ax.hist(nul, bins=bins, color="#2166ac", alpha=0.55, density=True, label="null (LOO)")
        ax.hist(sel, bins=bins, color="#b2182b", alpha=0.55, density=True, label="selected")
        ax.axvline(SIGNIFICANCE_LEVEL, color="black", ls="--", lw=1.1)
        ax.axhline(1.0, color="grey", ls=":", lw=1)
        power = float(np.mean(sel <= SIGNIFICANCE_LEVEL)) if sel.size else float("nan")
        ax.set_title(f"x = {int(years):,} y   power {power:.2f}", fontsize=9)
        ax.set_xlabel("p-value", fontsize=8)
        ax.set_ylabel("density", fontsize=8)
        ax.tick_params(labelsize=8)
    for ax in axes[len(thresholds) :]:
        ax.set_visible(False)
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle(
        f"{arm.label}\np-value distribution, {estimator}; dashed line alpha = "
        f"{SIGNIFICANCE_LEVEL:g}, dotted line = uniform",
        y=1.0,
    )
    fig.tight_layout()
    return _save(fig, figures / f"pvalue_distribution_{estimator}")


def plot_null_vs_selected(
    arm: Run5Arm, study_root: str | Path, figures: Path, *, estimator: str = "tree_truth"
) -> dict[str, Path]:
    """The raw statistic: null distribution with selected values overlaid."""
    statistics = _load(Path(study_root) / arm.arm_id / "results", "statistics")
    thresholds = [
        t for t in TMRCA_THRESHOLDS_YEARS if t in set(statistics["threshold_years"])
    ]
    columns = 4
    rows = int(np.ceil(len(thresholds) / columns))
    fig, axes = plt.subplots(rows, columns, figsize=(3.9 * columns, 3.0 * rows))
    axes = np.atleast_1d(axes).ravel()
    for ax, years in zip(axes, thresholds):
        selected, null_values = _split(statistics, estimator, float(years))
        obs = selected["value"].to_numpy(dtype=float)
        combined = np.concatenate([null_values, obs])
        combined = combined[np.isfinite(combined)]
        top = float(combined.max()) if combined.size else 1.0
        bins = np.linspace(0, max(top * 1.05, 1e-6), 26)
        ax.hist(null_values, bins=bins, color="#2166ac", alpha=0.7, label="null")
        ax.hist(obs, bins=bins, color="#b2182b", alpha=0.55, label="selected")
        ax.set_title(f"x = {int(years):,} y", fontsize=9)
        ax.set_xlabel("P(TMRCA < x)", fontsize=8)
        ax.set_ylabel("replicates", fontsize=8)
        ax.tick_params(labelsize=8)
    for ax in axes[len(thresholds) :]:
        ax.set_visible(False)
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle(f"{arm.label}\nnull versus selected statistic, {estimator}", y=1.0)
    fig.tight_layout()
    return _save(fig, figures / f"null_vs_selected_{estimator}")


def plot_auc(arm: Run5Arm, auc: pd.DataFrame, figures: Path) -> dict[str, Path]:
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6))
    for estimator, group in auc.groupby("estimator"):
        group = group.sort_values("threshold_years")
        axes[0].plot(group["threshold_years"], group["auc"], marker="o", lw=2,
                     color=_COLOURS.get(estimator), label=estimator)
        axes[1].plot(group["threshold_years"], group["power_at_alpha"], marker="o", lw=2,
                     color=_COLOURS.get(estimator), label=estimator)
    axes[0].axhline(0.5, color="grey", ls="--", lw=1)
    axes[0].set(xlabel="Threshold x (years)", ylabel="AUC", ylim=(0.4, 1.02),
                title="Threshold-free discrimination (AUC)")
    axes[1].axhline(SIGNIFICANCE_LEVEL, color="grey", ls="--", lw=1)
    axes[1].set(xlabel="Threshold x (years)",
                ylabel=f"Fraction with p <= {SIGNIFICANCE_LEVEL:g}", ylim=(0, 1.02),
                title="Per-replicate detection rate (power)")
    for ax in axes:
        ax.legend(frameon=False, fontsize=9)
    fig.suptitle(f"{arm.label} - AUC versus power", y=1.02)
    return _save(fig, figures / "auc_vs_power")


def plot_auc_vs_final_af(
    arm: Run5Arm,
    percentiles: pd.DataFrame,
    figures: Path,
    *,
    estimator: str = "tree_truth",
    thresholds: Sequence[float] = (20_000.0, 30_000.0, 40_000.0),
) -> dict[str, Path]:
    """Per-replicate percentile against the null, versus final allele frequency.

    Each point's y value is that replicate's contribution to the AUC, so the
    binned means are AUC restricted to that final-frequency band.
    """
    available = [t for t in thresholds if t in set(percentiles["threshold_years"])]
    fig, axes = plt.subplots(1, len(available), figsize=(4.6 * len(available), 4.4),
                             sharey=True)
    axes = np.atleast_1d(axes)
    edges = np.array([0.0, 0.4, 0.7, 1.01])
    for ax, years in zip(axes, available):
        frame = percentiles[
            (percentiles["estimator"] == estimator)
            & (percentiles["threshold_years"] == years)
        ]
        frame = frame[np.isfinite(frame["census_af"]) & np.isfinite(frame["percentile"])]
        if frame.empty:
            ax.set_visible(False)
            continue
        ax.scatter(frame["census_af"], frame["percentile"], s=18, alpha=0.55,
                   color="#b2182b")
        for lo, hi in zip(edges[:-1], edges[1:]):
            band = frame[(frame["census_af"] >= lo) & (frame["census_af"] < hi)]
            if band.empty:
                continue
            ax.plot([lo, min(hi, 1.0)], [band["percentile"].mean()] * 2,
                    color="black", lw=2.4)
            ax.text((lo + min(hi, 1.0)) / 2, band["percentile"].mean() + 0.03,
                    f"AUC {band['percentile'].mean():.2f}\n(n={len(band)})",
                    ha="center", fontsize=8)
        ax.axhline(0.5, color="grey", ls="--", lw=1)
        ax.set(xlabel="final census allele frequency", xlim=(0, 1), ylim=(0, 1.08),
               title=f"x = {int(years):,} y")
    axes[0].set_ylabel("percentile against the null (AUC contribution)")
    fig.suptitle(
        f"{arm.label}\ndiscrimination versus final allele frequency, {estimator}",
        y=1.02,
    )
    fig.tight_layout()
    return _save(fig, figures / f"auc_vs_final_af_{estimator}")


def analyse_arm_power(study_root: str | Path, arm: Run5Arm) -> dict[str, Any]:
    root = Path(study_root)
    results = root / arm.arm_id / "results"
    figures = root / arm.arm_id / "figures"
    percentiles = replicate_percentiles(root, arm)
    null_p = null_percentiles(root, arm)
    auc = auc_summary(percentiles, arm)

    percentiles.to_csv(results / POWER_OUTPUTS["percentiles"], sep="\t", index=False)
    auc.to_csv(results / POWER_OUTPUTS["auc"], sep="\t", index=False)

    written: dict[str, Any] = {
        "auc_table": str(results / POWER_OUTPUTS["auc"]),
        "percentiles_table": str(results / POWER_OUTPUTS["percentiles"]),
        "auc_vs_power": plot_auc(arm, auc, figures),
    }
    for estimator in ("tree_truth", "gamma_smc_frac"):
        written[f"pvalue_distribution_{estimator}"] = plot_pvalue_distribution(
            arm, percentiles, null_p, figures, estimator=estimator
        )
        written[f"null_vs_selected_{estimator}"] = plot_null_vs_selected(
            arm, root, figures, estimator=estimator
        )
    if percentiles["census_af"].notna().any():
        written["auc_vs_final_af"] = plot_auc_vs_final_af(arm, percentiles, figures)
    return written


def analyse_all(study_root: str | Path) -> dict[str, Any]:
    return {arm.arm_id: analyse_arm_power(study_root, arm) for arm in ARMS}


__all__ = [
    "POWER_OUTPUTS",
    "analyse_all",
    "analyse_arm_power",
    "auc_summary",
    "null_percentiles",
    "percentile_against_null",
    "plot_auc",
    "plot_auc_vs_final_af",
    "plot_null_vs_selected",
    "plot_pvalue_distribution",
    "replicate_percentiles",
]
