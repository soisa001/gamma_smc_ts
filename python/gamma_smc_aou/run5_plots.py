"""Figures for the run5 study.

The two additions relative to run2 are that every power figure carries one line
per estimator (tree truth versus the two Gamma-SMC statistics), and that power
is shown both unconditionally and conditional on the focal allele being present
at all -- which in run5 is the difference between including and excluding the
majority of replicates that carry no archaic ancestry at the focal base.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .run5_analysis import RESULT_OUTPUTS, cross_arm_power
from .run5_config import ARMS, SIGNIFICANCE_LEVEL, Run5Arm, TMRCA_THRESHOLDS_YEARS

_COLOURS = {
    "tree_truth": "#111111",
    "gamma_smc_frac": "#b2182b",
    "gamma_smc_meanp": "#ef8a62",
}
_MODE_COLOURS = {"selected": "#b2182b", "neutral": "#2166ac"}


def _save(fig, stem: Path) -> dict[str, Path]:
    stem.parent.mkdir(parents=True, exist_ok=True)
    written = {}
    for suffix in ("png", "pdf"):
        path = stem.with_suffix(f".{suffix}")
        fig.savefig(path, dpi=170, bbox_inches="tight")
        written[suffix] = path
    plt.close(fig)
    return written


def _read(results: Path, key: str) -> pd.DataFrame:
    path = results / RESULT_OUTPUTS[key]
    if not path.is_file():
        raise FileNotFoundError(f"missing result table: {path}")
    return pd.read_csv(path, sep="\t")


# ---------------------------------------------------------------------------


def plot_final_af(arm: Run5Arm, results: Path, figures: Path) -> dict[str, Path]:
    frame = _read(results, "final_af")
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), width_ratios=[1, 2, 2])

    absent = frame.groupby("mode")["allele_present"].apply(lambda s: 1 - s.mean())
    axes[0].bar(absent.index, absent.to_numpy(), color=[_MODE_COLOURS[m] for m in absent.index])
    axes[0].set(ylabel="Fraction with no allele at the focal base", title="Allele absent", ylim=(0, 1))
    for i, v in enumerate(absent.to_numpy()):
        axes[0].text(i, v, f"{v:.2f}", ha="center", va="bottom", fontsize=9)

    for mode, colour in _MODE_COLOURS.items():
        values = frame.loc[frame["mode"] == mode, "census_af"].to_numpy(dtype=float)
        values = values[np.isfinite(values)]
        present = values[values > 0]
        axes[1].hist(present, bins=np.linspace(0, max(0.05, present.max() if present.size else 0.05), 30),
                     color=colour, alpha=0.55, label=f"{mode} (n={present.size} present)")
        axes[2].hist(values, bins=np.linspace(0, 1, 41), color=colour, alpha=0.55, label=mode)
    axes[1].set(xlabel="Census allele frequency", ylabel="Replicates",
                title="Present replicates, low-frequency detail")
    axes[2].set(xlabel="Census allele frequency", ylabel="Replicates",
                title="Full range", xlim=(0, 1))
    for ax in axes[1:]:
        ax.legend(frameon=False, fontsize=8)
    fig.suptitle(f"{arm.label} — focal allele outcome", y=1.03)
    return _save(fig, figures / "final_af")


def plot_power_by_estimator(arm: Run5Arm, results: Path, figures: Path) -> dict[str, Path]:
    power = _read(results, "power")
    power = power[power["genotype_class"] == "overall"]
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6), sharey=True)
    for ax, conditioned in zip(axes, (False, True)):
        subset = power[power["conditioned_on_presence"] == conditioned]
        for estimator, group in subset.groupby("estimator"):
            group = group.sort_values("threshold_years")
            if group.empty:
                continue
            ax.plot(group["threshold_years"], group["power"], marker="o", lw=2,
                    color=_COLOURS.get(estimator, None), label=estimator)
            n = int(group["n_selected"].iloc[0])
        ax.axhline(SIGNIFICANCE_LEVEL, color="grey", ls="--", lw=1)
        ax.set(xlabel="Threshold x (years)", xscale="log", ylim=(0, 1.02),
               title=("conditional on the allele being present" if conditioned
                      else "all replicates"))
        ax.legend(frameon=False, fontsize=9)
    axes[0].set_ylabel(f"Fraction of replicates with p <= {SIGNIFICANCE_LEVEL:g}")
    fig.suptitle(f"{arm.label} — detection power, tree truth vs Gamma-SMC", y=1.02)
    return _save(fig, figures / "power_by_estimator")


def plot_truth_vs_decoded(arm: Run5Arm, results: Path, figures: Path) -> dict[str, Path]:
    statistics = _read(results, "statistics")
    overall = statistics[statistics["genotype_class"] == "overall"]
    wide = overall.pivot_table(
        index=["mode", "replicate_index", "threshold_years"],
        columns="estimator", values="value",
    ).reset_index()
    if "gamma_smc_frac" not in wide.columns:
        raise RuntimeError("no decoded statistics are available")
    thresholds = [t for t in TMRCA_THRESHOLDS_YEARS if t in set(wide["threshold_years"])]
    thresholds = thresholds[-6:] if len(thresholds) > 6 else thresholds
    columns = 3
    rows = int(np.ceil(len(thresholds) / columns))
    fig, axes = plt.subplots(rows, columns, figsize=(4.2 * columns, 3.7 * rows))
    axes = np.atleast_1d(axes).ravel()
    for ax, years in zip(axes, thresholds):
        sub = wide[wide["threshold_years"] == years]
        for mode, colour in _MODE_COLOURS.items():
            s = sub[sub["mode"] == mode]
            ax.scatter(s["tree_truth"], s["gamma_smc_frac"], s=14, alpha=0.6,
                       color=colour, label=mode)
        lo = 0.0
        hi = float(np.nanmax([sub["tree_truth"].max(), sub["gamma_smc_frac"].max(), 1e-3])) * 1.05
        ax.plot([lo, hi], [lo, hi], color="grey", ls="--", lw=1)
        ax.set(xlabel="tree truth", ylabel="Gamma-SMC frac_recent",
               title=f"x = {int(years):,} y", xlim=(lo, hi), ylim=(lo, hi))
        ax.tick_params(labelsize=8)
    for ax in axes[len(thresholds):]:
        ax.set_visible(False)
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle(f"{arm.label} — decoded vs true statistic, per replicate", y=1.0)
    fig.tight_layout()
    return _save(fig, figures / "truth_vs_decoded")


def plot_spatial_profile(
    arm: Run5Arm, results: Path, figures: Path,
    thresholds: Sequence[float] = (20_000.0, 50_000.0, 100_000.0, 200_000.0),
) -> dict[str, Path]:
    frame = _read(results, "spatial_profiles")
    available = [t for t in thresholds if t in set(frame["threshold_years"])]
    fig, axes = plt.subplots(len(available), 1, figsize=(11, 2.4 * len(available)), sharex=True)
    axes = np.atleast_1d(axes)
    for ax, years in zip(axes, available):
        for mode, colour in _MODE_COLOURS.items():
            sub = frame[(frame["mode"] == mode) & (frame["threshold_years"] == years)]
            sub = sub.sort_values("position_0based")
            if sub.empty:
                continue
            x = sub["position_0based"].to_numpy() / 1e6
            ax.fill_between(x, sub["q025"], sub["q975"], color=colour, alpha=0.2, lw=0)
            ax.plot(x, sub["mean"], color=colour, lw=1.7, label=mode)
        ax.axvline(5.0, color="black", ls=":", lw=1.1)
        ax.set_ylabel(f"P(TMRCA <\n{int(years):,} y)", fontsize=9)
        ax.tick_params(labelsize=8)
    axes[0].legend(frameon=False, fontsize=9, ncol=2)
    axes[-1].set_xlabel("Position (Mb); dotted line marks the focal base", fontsize=10)
    fig.suptitle(f"{arm.label} — tree-truth spatial profile (mean, 95% band)", y=0.998)
    fig.tight_layout()
    return _save(fig, figures / "spatial_profile")


def plot_arm(study_root: str | Path, arm: Run5Arm) -> dict[str, Any]:
    root = Path(study_root)
    results = root / arm.arm_id / "results"
    figures = root / arm.arm_id / "figures"
    written: dict[str, Any] = {"final_af": plot_final_af(arm, results, figures)}
    for name, fn in (
        ("power_by_estimator", plot_power_by_estimator),
        ("truth_vs_decoded", plot_truth_vs_decoded),
        ("spatial_profile", plot_spatial_profile),
    ):
        try:
            written[name] = fn(arm, results, figures)
        except (RuntimeError, FileNotFoundError, ValueError) as error:
            written[name] = {"skipped": str(error)}
    return written


def plot_cross_arm(study_root: str | Path) -> dict[str, Path]:
    root = Path(study_root)
    frame = cross_arm_power(root, ARMS)
    if frame.empty:
        raise RuntimeError("no per-arm power tables were found")
    frame = frame[(frame["genotype_class"] == "overall") & (frame["estimator"] == "tree_truth")]
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.4), sharey=True)
    for ax, conditioned in zip(axes, (False, True)):
        subset = frame[frame["conditioned_on_presence"] == conditioned]
        for arm_id, group in subset.groupby("arm_id"):
            group = group.sort_values("threshold_years")
            ax.plot(group["threshold_years"], group["power"], marker="o", lw=2,
                    label=group["arm_label"].iloc[0])
        ax.axhline(SIGNIFICANCE_LEVEL, color="grey", ls="--", lw=1)
        ax.set(xlabel="Threshold x (years)", xscale="log", ylim=(0, 1.02),
               title=("conditional on presence" if conditioned else "all replicates"))
        ax.legend(frameon=False, fontsize=9)
    axes[0].set_ylabel(f"Fraction with p <= {SIGNIFICANCE_LEVEL:g}")
    fig.suptitle("run5 — tree-truth power by arm", y=1.02)
    return _save(fig, root / "cross_arm" / "figures" / "power_by_arm")


__all__ = [
    "plot_arm",
    "plot_cross_arm",
    "plot_final_af",
    "plot_power_by_estimator",
    "plot_spatial_profile",
    "plot_truth_vs_decoded",
]
