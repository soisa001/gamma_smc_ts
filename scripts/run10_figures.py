#!/usr/bin/env python
"""Figures for run10: what survives a frequency-matched null."""

from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "python"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

RUN9 = REPO / "sim_results_run9"
RUN10 = REPO / "sim_results_run10"
OUT = RUN10 / "figures"
AF_EDGES = (0.10, 0.15, 0.25, 0.35, 0.45, 0.60, 0.71)
PULSE_ARMS = ("pulse_s0p003", "pulse_s0p004", "pulse_s0p005", "pulse_delayed_s0p005")


def figure_distributions() -> None:
    """The whole argument in one panel: which null you use decides the answer."""
    null = pd.read_csv(RUN10 / "tables" / "matched_null_raw.tsv", sep="\t")
    null["af_bin"] = pd.cut(null["archaic_frequency"], AF_EDGES, right=False)
    selected = pd.read_csv(RUN9 / "tables" / "regional_raw.tsv", sep="\t")
    selected = selected[
        (selected["window_bp"] == 10_000)
        & (selected["genotype_class"] == "carrier_carrier")
        & (selected["arm"].isin(PULSE_ARMS))
    ].copy()
    selected["af_bin"] = pd.cut(selected["sample_af"], AF_EDGES, right=False)
    # The carrier-free null used by every earlier run, for contrast.
    carrier_free = pd.read_csv(RUN9 / "tables" / "regional_raw.tsv", sep="\t")
    carrier_free = carrier_free[
        (carrier_free["window_bp"] == 10_000)
        & (carrier_free["arm"] == "neutral")
        & (carrier_free["genotype_class"] == "overall")
    ]

    bands = ["[0.25, 0.35)", "[0.35, 0.45)", "[0.45, 0.6)"]
    thresholds = [10_000.0, 30_000.0]
    fig, axes = plt.subplots(len(thresholds), len(bands), figsize=(15, 7.5))
    for row, years in enumerate(thresholds):
        for column, band in enumerate(bands):
            axis = axes[row, column]
            ref = null[
                (null["af_bin"].astype(str) == band)
                & (null["threshold_years"] == years)
            ]["carrier_carrier"].to_numpy(dtype=float)
            sel = selected[
                (selected["af_bin"].astype(str) == band)
                & (selected["threshold_years"] == years)
            ]["statistic"].to_numpy(dtype=float)
            free = carrier_free[carrier_free["threshold_years"] == years][
                "statistic"
            ].to_numpy(dtype=float)
            if ref.size == 0 or sel.size == 0:
                axis.axis("off")
                continue
            hi = max(ref.max(), sel.max(), free.max()) * 1.02
            bins = np.linspace(0, hi, 40)
            axis.hist(free, bins=bins, density=True, color="#bbbbbb", alpha=0.65,
                      label="carrier-free null\n(used in runs 6–9)")
            axis.hist(ref, bins=bins, density=True, histtype="step", lw=1.8,
                      color="#023047", label="matched null\n(neutral introgression)")
            axis.hist(sel, bins=bins, density=True, histtype="step", lw=1.8,
                      color="#c1121f", label="selected")
            axis.axvline(np.quantile(ref, 0.95), color="#023047", ls="--", lw=1.0)
            if row == 0:
                axis.set_title(f"AF {band}", fontsize=10)
            if column == 0:
                axis.set_ylabel(f"P(TMRCA < {years:,.0f} y)\ndensity")
            if row == len(thresholds) - 1:
                axis.set_xlabel("carrier–carrier statistic")
    axes[0, 0].legend(fontsize=7, loc="upper right")
    fig.suptitle(
        "The null decides the answer: selected regions separate cleanly from a "
        "carrier-free null, far less from neutral introgression at the same frequency",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(OUT / "run10_distributions.png", dpi=150)
    plt.close(fig)


def figure_auc_power() -> None:
    table = pd.read_csv(RUN10 / "tables" / "headtohead.tsv", sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))
    for years, colour in ((10_000.0, "#219ebc"), (30_000.0, "#c1121f"),
                          (50_000.0, "#8ecae6")):
        part = table[table["threshold_years"] == years].sort_values("af_bin")
        if part.empty:
            continue
        axes[0].plot(part["af_bin"], part["auc"], marker="o", ms=5, lw=1.6,
                     color=colour, label=f"{years:,.0f} y")
        axes[1].plot(part["af_bin"], part["power_p05"], marker="o", ms=5, lw=1.6,
                     color=colour, label=f"{years:,.0f} y")
    axes[0].axhline(0.5, color="black", ls=":", lw=0.9)
    axes[0].set_ylabel("AUC vs frequency-matched null")
    axes[0].set_ylim(0.4, 1.02)
    axes[1].axhline(0.05, color="black", ls=":", lw=0.9)
    axes[1].set_ylabel("power at p ≤ 0.05")
    axes[1].set_ylim(0, 1.02)
    for axis in axes:
        axis.set_xlabel("archaic frequency band")
        axis.tick_params(axis="x", rotation=30, labelsize=8)
        axis.legend(fontsize=8)
    fig.suptitle(
        "Selected introgressed regions against neutrally drifted ones at the same "
        "frequency — real but weak separation",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(OUT / "run10_auc_power.png", dpi=150)
    plt.close(fig)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    figure_distributions()
    figure_auc_power()
    print("figures ->", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
