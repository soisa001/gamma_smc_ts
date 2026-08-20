"""Figures for the run2 study.

Every figure is written as both PNG and PDF.  Figures 2 and 3 deliberately share
their y-axis so the selected and neutral panels can be compared by eye without
reading the tick labels.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .run2_analysis import (
    RESULT_OUTPUTS,
    TEST_CLASS,
    cross_arm_comparison,
    representative_replicates,
)
from .run2_config import (
    ARMS,
    SIGNIFICANCE_LEVEL,
    STANDING_FREQUENCY,
    Run2Arm,
    TMRCA_THRESHOLDS_YEARS,
)

_SELECTED_COLOUR = "#b2182b"
_NEUTRAL_COLOUR = "#2166ac"
_SPATIAL_THRESHOLDS = (4_500.0, 20_000.0)


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


def _test_frame(statistics: pd.DataFrame, mode: str) -> pd.DataFrame:
    return statistics[
        (statistics["mode"] == mode) & (statistics["genotype_class"] == TEST_CLASS)
    ]


def _threshold_matrix(frame: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """Return the threshold axis and a ``(n_replicates, n_thresholds)`` matrix."""
    wide = frame.pivot_table(
        index="replicate_index", columns="threshold_years", values="p_tmrca_lt_threshold"
    ).sort_index(axis=1)
    return wide.columns.to_numpy(dtype=float), wide.to_numpy(dtype=float)


# ---------------------------------------------------------------------------
# Figure 1: final allele frequency
# ---------------------------------------------------------------------------


def plot_final_af(arm: Run2Arm, results: Path, figures: Path) -> dict[str, Path]:
    frame = _read(results, "final_af")
    values = frame["census_af"].to_numpy(dtype=float)
    values = values[np.isfinite(values)]
    lost = int(np.sum(values == 0.0))
    fixed = int(np.sum(values == 1.0))

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), width_ratios=[1, 3])
    axes[0].bar(
        ["lost", "segregating", "fixed"],
        [lost, len(values) - lost - fixed, fixed],
        color=["#999999", _SELECTED_COLOUR, "#4d004b"],
    )
    axes[0].set(ylabel="Replicates", title="Outcome")
    for index, count in enumerate([lost, len(values) - lost - fixed, fixed]):
        axes[0].text(index, count, str(count), ha="center", va="bottom", fontsize=9)

    surviving = values[values > 0.0]
    axes[1].hist(
        surviving, bins=np.linspace(0, 1, 41), color=_SELECTED_COLOUR, alpha=0.85
    )
    axes[1].axvline(
        STANDING_FREQUENCY,
        color="black",
        ls="--",
        lw=1,
        label=f"starting frequency {STANDING_FREQUENCY:g}",
    )
    axes[1].set(
        xlabel="Present-day census allele frequency",
        ylabel="Replicates",
        title=f"Surviving replicates (n = {len(surviving)}); {lost} lost, not shown",
        xlim=(0, 1),
    )
    axes[1].legend(frameon=False, fontsize=9)
    fig.suptitle(f"{arm.label} — final frequency of the selected allele", y=1.02)
    return _save(fig, figures / "final_af")


# ---------------------------------------------------------------------------
# Figures 2 and 3: P(TMRCA < x) against threshold
# ---------------------------------------------------------------------------


def _plot_threshold_panel(ax, thresholds, matrix, colour, title) -> None:
    for row in matrix:
        ax.plot(thresholds, row, color=colour, alpha=0.08, lw=0.8)
    median = np.nanmedian(matrix, axis=0)
    lower = np.nanquantile(matrix, 0.25, axis=0)
    upper = np.nanquantile(matrix, 0.75, axis=0)
    ax.fill_between(thresholds, lower, upper, color=colour, alpha=0.25, lw=0)
    ax.plot(thresholds, median, color=colour, lw=2.4, label="median")
    ax.set(
        xlabel="Threshold x (years)",
        ylabel="P(TMRCA < x) at the focal base",
        title=title,
        xscale="log",
    )
    ax.legend(frameon=False, fontsize=9, loc="upper left")


def plot_threshold_curves(arm: Run2Arm, results: Path, figures: Path) -> dict[str, Any]:
    statistics = _read(results, "observed_statistics")
    written: dict[str, Any] = {}
    panels = {}
    for mode, colour in (("selected", _SELECTED_COLOUR), ("neutral", _NEUTRAL_COLOUR)):
        thresholds, matrix = _threshold_matrix(_test_frame(statistics, mode))
        panels[mode] = (thresholds, matrix, colour)

    top = max(
        float(np.nanmax(matrix)) for _, matrix, _ in panels.values() if matrix.size
    )
    ylim = (0.0, min(1.0, max(0.05, top * 1.15)))

    for mode, (thresholds, matrix, colour) in panels.items():
        fig, ax = plt.subplots(figsize=(7.2, 4.6))
        _plot_threshold_panel(
            ax,
            thresholds,
            matrix,
            colour,
            f"{arm.label} — {mode} (n = {matrix.shape[0]})",
        )
        ax.set_ylim(*ylim)
        written[mode] = _save(fig, figures / f"ptmrca_{mode}")
    return written


# ---------------------------------------------------------------------------
# Figure 4: null distributions and p-values
# ---------------------------------------------------------------------------


def plot_pvalue_panels(arm: Run2Arm, results: Path, figures: Path) -> dict[str, Path]:
    statistics = _read(results, "observed_statistics")
    power = _read(results, "power")
    neutral = _test_frame(statistics, "neutral")
    selected = _test_frame(statistics, "selected")

    thresholds = list(TMRCA_THRESHOLDS_YEARS)
    columns = 4
    rows = int(np.ceil(len(thresholds) / columns))
    fig, axes = plt.subplots(rows, columns, figsize=(4 * columns, 3.1 * rows))
    axes = np.atleast_1d(axes).ravel()

    for index, threshold in enumerate(thresholds):
        ax = axes[index]
        null_values = neutral.loc[
            neutral["threshold_years"] == threshold, "p_tmrca_lt_threshold"
        ].to_numpy(dtype=float)
        observed = selected.loc[
            selected["threshold_years"] == threshold, "p_tmrca_lt_threshold"
        ].to_numpy(dtype=float)
        combined = np.concatenate([null_values, observed])
        combined = combined[np.isfinite(combined)]
        bins = np.linspace(0, max(combined.max(), 1e-6) * 1.05, 26) if combined.size else 10
        ax.hist(null_values, bins=bins, color=_NEUTRAL_COLOUR, alpha=0.75, label="null")
        ax.hist(
            observed, bins=bins, color=_SELECTED_COLOUR, alpha=0.55, label="selected"
        )
        row = power[power["threshold_years"] == threshold]
        if len(row):
            record = row.iloc[0]
            ax.axvline(
                float(record["median_observed"]), color="black", ls="--", lw=1.2
            )
            ax.set_title(
                f"x = {int(threshold):,} y\n"
                f"power {record['power']:.2f}, median p {record['median_p_value']:.3f}",
                fontsize=9,
            )
        else:
            ax.set_title(f"x = {int(threshold):,} y", fontsize=9)
        ax.set_xlabel("P(TMRCA < x)", fontsize=8)
        ax.set_ylabel("Replicates", fontsize=8)
        ax.tick_params(labelsize=8)
        if index == 0:
            ax.legend(frameon=False, fontsize=8)

    for ax in axes[len(thresholds) :]:
        ax.set_visible(False)
    fig.suptitle(
        f"{arm.label} — null vs selected, one-sided upper tail "
        f"(alpha = {SIGNIFICANCE_LEVEL:g}; dashed line = median selected replicate)",
        y=1.0,
    )
    fig.tight_layout()
    return _save(fig, figures / "pvalue_panels")


# ---------------------------------------------------------------------------
# Figure 5: spatial profile
# ---------------------------------------------------------------------------


def plot_spatial_profile(
    arm: Run2Arm,
    results: Path,
    figures: Path,
    thresholds: Sequence[float] = _SPATIAL_THRESHOLDS,
) -> dict[str, Path]:
    frame = _read(results, "spatial_profiles")
    fig, axes = plt.subplots(
        len(thresholds), 1, figsize=(9, 3.2 * len(thresholds)), sharex=True
    )
    axes = np.atleast_1d(axes)
    for ax, threshold in zip(axes, thresholds):
        for mode, colour in (
            ("neutral", _NEUTRAL_COLOUR),
            ("selected", _SELECTED_COLOUR),
        ):
            subset = frame[
                (frame["mode"] == mode) & (frame["threshold_years"] == threshold)
            ].sort_values("position_0based")
            if subset.empty:
                continue
            positions = subset["position_0based"].to_numpy() / 1e6
            ax.fill_between(
                positions,
                subset["q05"].to_numpy(),
                subset["q95"].to_numpy(),
                color=colour,
                alpha=0.2,
                lw=0,
            )
            ax.plot(
                positions, subset["mean"].to_numpy(), color=colour, lw=1.6, label=mode
            )
        ax.axvline(5.0, color="black", ls=":", lw=1)
        ax.set(ylabel=f"P(TMRCA < {int(threshold):,} y)")
        ax.legend(frameon=False, fontsize=9, loc="upper right")
    axes[-1].set_xlabel("Position (Mb); dotted line marks the focal base")
    fig.suptitle(f"{arm.label} — spatial profile (mean, 5-95% band)", y=0.99)
    fig.tight_layout()
    return _save(fig, figures / "spatial_profile")


# ---------------------------------------------------------------------------
# Figure 6: descriptive carrier stratification
# ---------------------------------------------------------------------------


def plot_by_genotype(arm: Run2Arm, results: Path, figures: Path) -> dict[str, Path]:
    statistics = _read(results, "observed_statistics")
    selected = statistics[statistics["mode"] == "selected"]
    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    colours = {
        "overall": "black",
        "hom_carrier": "#b2182b",
        "het": "#ef8a62",
        "hom_noncarrier": "#67a9cf",
    }
    for genotype_class, colour in colours.items():
        subset = selected[selected["genotype_class"] == genotype_class]
        if subset.empty:
            continue
        summary = subset.groupby("threshold_years")["p_tmrca_lt_threshold"].median()
        counts = subset.groupby("threshold_years")["n_pairs"].mean()
        if summary.isna().all():
            continue
        ax.plot(
            summary.index.to_numpy(dtype=float),
            summary.to_numpy(dtype=float),
            color=colour,
            lw=2,
            marker="o",
            ms=4,
            label=f"{genotype_class} (mean n = {counts.mean():.0f})",
        )
    ax.set(
        xlabel="Threshold x (years)",
        ylabel="Median P(TMRCA < x) at the focal base",
        title=f"{arm.label} — descriptive carrier stratification (not the test)",
        xscale="log",
    )
    ax.legend(frameon=False, fontsize=9)
    return _save(fig, figures / "ptmrca_by_genotype")


# ---------------------------------------------------------------------------
# Figure 7: cross-arm comparison
# ---------------------------------------------------------------------------


def plot_arm_comparison(study_root: Path, figures: Path) -> dict[str, Path]:
    frame = cross_arm_comparison(study_root, ARMS)
    if frame.empty:
        raise RuntimeError("no per-arm power tables were found")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4))
    markers = {"chb_ancient_eurasia": "o", "eas_phlash": "s"}
    for arm_id, group in frame.groupby("arm_id"):
        group = group.sort_values("threshold_years")
        axes[0].plot(
            group["threshold_years"],
            group["power"],
            marker=markers.get(arm_id, "o"),
            lw=2,
            label=group["arm_label"].iloc[0],
        )
        axes[1].plot(
            group["threshold_years"],
            group["median_observed"],
            marker=markers.get(arm_id, "o"),
            lw=2,
            label=group["arm_label"].iloc[0],
        )
    axes[0].axhline(SIGNIFICANCE_LEVEL, color="grey", ls="--", lw=1)
    axes[0].set(
        xlabel="Threshold x (years)",
        ylabel=f"Fraction of replicates with p <= {SIGNIFICANCE_LEVEL:g}",
        title="Detection power",
        xscale="log",
        ylim=(0, 1.02),
    )
    axes[1].set(
        xlabel="Threshold x (years)",
        ylabel="Median P(TMRCA < x)",
        title="Median selected statistic",
        xscale="log",
    )
    for ax in axes:
        ax.legend(frameon=False, fontsize=9)
    fig.suptitle("run2 — introgressed vs non-introgressed standing variation", y=1.02)
    return _save(fig, figures / "arm_comparison")


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------


def plot_arm(study_root: str | Path, arm: Run2Arm) -> dict[str, Any]:
    root = Path(study_root)
    results = root / arm.arm_id / "results"
    figures = root / arm.arm_id / "figures"
    return {
        "final_af": plot_final_af(arm, results, figures),
        "threshold_curves": plot_threshold_curves(arm, results, figures),
        "pvalue_panels": plot_pvalue_panels(arm, results, figures),
        "spatial_profile": plot_spatial_profile(arm, results, figures),
        "by_genotype": plot_by_genotype(arm, results, figures),
        "selected_vs_neutral_band": plot_selected_vs_neutral_band(root, arm, figures),
    }


def plot_cross_arm(study_root: str | Path) -> dict[str, Path]:
    root = Path(study_root)
    return plot_arm_comparison(root, root / "cross_arm" / "figures")


__all__ = [
    "plot_arm",
    "plot_arm_comparison",
    "plot_by_genotype",
    "plot_cross_arm",
    "plot_final_af",
    "plot_pvalue_panels",
    "plot_selected_vs_neutral_band",
    "plot_spatial_profile",
    "plot_threshold_curves",
]


def plot_selected_vs_neutral_band(
    study_root: str | Path,
    arm: Run2Arm,
    figures: Path | None = None,
    *,
    thresholds: Sequence[float] = TMRCA_THRESHOLDS_YEARS,
    smooth_bp: int = 250_000,
) -> dict[str, Path]:
    """One panel per cutoff: a single selected replicate against the whole null.

    A single neutral replicate is far too noisy to read against -- with 100
    diploid pairs its trace moves in steps of 0.01 -- so the neutral side is all
    100 replicates summarised as a median and a 2.5-97.5% band.  The selected
    side stays a single replicate on purpose: averaging the selected arm would
    blend the sweep with the ~20% of replicates that lost the allele and are
    neutral by construction.

    Every curve, including the band edges, is smoothed with the same rolling
    mean so the comparison is like for like.
    """
    root = Path(study_root)
    figures = figures if figures is not None else root / arm.arm_id / "figures"
    chosen = representative_replicates(root, arm)
    if "selected" not in chosen:
        raise RuntimeError("no representative selected replicate is available")
    record = chosen["selected"]

    path = Path(record["spatial_profile"])
    if not path.is_file():
        raise FileNotFoundError(
            f"per-replicate profile is missing: {path}. These live under "
            "<arm>/<mode>/replicates/ and are gitignored, so this plot must be "
            "produced where the study was run."
        )
    selected = pd.read_csv(path, sep="\t").sort_values("position_0based")
    aggregate = _read(root / arm.arm_id / "results", "spatial_profiles")
    neutral = aggregate[aggregate["mode"] == "neutral"]
    if neutral.empty:
        raise RuntimeError("the aggregate spatial profile has no neutral rows")
    for column in ("q025", "q975"):
        if column not in neutral.columns:
            raise RuntimeError(
                "the aggregate spatial profile predates the 95% band; re-run "
                "the analyze phase"
            )

    positions = selected["position_0based"].to_numpy(dtype=float)
    step = float(np.median(np.diff(positions))) if len(positions) > 1 else 1.0
    window = max(1, int(round(smooth_bp / step)))

    def smooth(values):
        return (
            pd.Series(np.asarray(values, dtype=float))
            .rolling(window, center=True, min_periods=1)
            .mean()
            .to_numpy()
        )

    n_neutral = int(neutral["n_replicates"].iloc[0])
    fig, axes = plt.subplots(
        len(thresholds), 1, figsize=(11, 2.05 * len(thresholds)), sharex=True
    )
    axes = np.atleast_1d(axes)
    for ax, threshold in zip(axes, thresholds):
        band = neutral[neutral["threshold_years"] == threshold].sort_values(
            "position_0based"
        )
        x = band["position_0based"].to_numpy(dtype=float) / 1e6
        ax.fill_between(
            x,
            smooth(band["q025"]),
            smooth(band["q975"]),
            color=_NEUTRAL_COLOUR,
            alpha=0.25,
            lw=0,
            label=f"neutral 95% band (n = {n_neutral})",
        )
        ax.plot(
            x, smooth(band["median"]), color=_NEUTRAL_COLOUR, lw=1.7, label="neutral median"
        )

        column = f"p_lt_{int(threshold)}y"
        values = selected[column].to_numpy(dtype=float)
        ax.plot(positions / 1e6, values, color=_SELECTED_COLOUR, lw=0.5, alpha=0.25)
        ax.plot(
            positions / 1e6,
            smooth(values),
            color=_SELECTED_COLOUR,
            lw=2.0,
            label="selected replicate",
        )
        ax.axvline(5.0, color="black", ls=":", lw=1.1)
        ax.set_ylabel(f"P(TMRCA <\n{int(threshold):,} y)", fontsize=9)
        ax.tick_params(labelsize=8)
        ax.margins(x=0.005)
    axes[0].legend(frameon=False, fontsize=8.5, ncol=3, loc="upper left")
    axes[-1].set_xlabel(
        "Position (Mb); dotted line marks the focal base at 5 Mb", fontsize=10
    )
    fig.suptitle(
        f"{arm.label} — selected rep{record['replicate_index']:03d}"
        f"{' (fixed)' if record['restricted_to_fixed'] else ''} "
        f"against the neutral distribution\n"
        f"all curves smoothed with a {smooth_bp // 1000} kb rolling mean; "
        f"faint red = raw 10 kb grid",
        fontsize=11,
        y=0.997,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    return _save(fig, figures / "selected_vs_neutral_band_by_threshold")
