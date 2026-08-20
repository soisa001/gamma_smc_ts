#!/usr/bin/env python
"""Figures for run6."""

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

OUT = REPO / "sim_results_run6" / "figures"
FOCAL = 5_000_000
CLASS_STYLE = {
    "hom_carrier": ("#1b5e9e", "hom carrier"),
    "het": ("#b8860b", "heterozygote"),
    "hom_noncarrier": ("#777777", "hom non-carrier"),
    "overall": ("#2e7d32", "all individuals"),
}


def figure_scan() -> None:
    scan = pd.read_csv("sim_results_run6/eas_standing/scan.tsv", sep="\t")
    thresholds = [10_000.0, 20_000.0, 30_000.0, 50_000.0]
    fig, axes = plt.subplots(1, 4, figsize=(17, 4.0), sharey=True)
    for axis, years in zip(axes, thresholds):
        block = scan[scan["threshold_years"] == years]
        for name, (colour, label) in CLASS_STYLE.items():
            part = block[block["genotype_class"] == name].sort_values("position_0based")
            if part.empty:
                continue
            axis.plot(
                (part["position_0based"] - FOCAL) / 1e6,
                part["auc_vs_neutral"],
                color=colour,
                lw=1.2,
                label=label,
            )
        axis.axhline(0.5, color="black", lw=0.8, ls=":")
        axis.axvline(0.0, color="red", lw=0.8, ls="--", alpha=0.6)
        axis.set_title(f"P(TMRCA < {years:,.0f} y)")
        axis.set_xlabel("distance from focal site (Mb)")
        axis.set_xlim(-2.5, 2.5)
    axes[0].set_ylabel("AUC vs neutral")
    axes[0].set_ylim(0, 1)
    axes[0].legend(fontsize=8, loc="lower left")
    fig.suptitle(
        "run6 EAS arm (SLiM, full PHLASH history): genotype-stratified scan",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(OUT / "run6_eas_scan.png", dpi=150)
    plt.close(fig)


def figure_conditioning() -> None:
    """How much does the *ascertainment alone* shift the null?"""
    eas = pd.read_csv(
        "sim_results_run6/conditioned_null/conditioned_null_summary.tsv", sep="\t"
    )
    chb = pd.read_csv(
        "sim_results_run6/conditioned_null/conditioned_null_chb_summary.tsv", sep="\t"
    )
    tract = pd.read_csv(
        "sim_results_run6/introgressed_null/introgressed_null_summary.tsv", sep="\t"
    )
    tract["band_low"] = (
        tract["frequency_bin"].str.extract(r"\[([0-9.]+),")[0].astype(float)
    )
    tract["ratio"] = tract["mean_hom"] / tract["mean_uncond"]

    thresholds = [20_000.0, 30_000.0, 50_000.0]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), sharey=True)
    for axis, years in zip(axes, thresholds):
        for frame, colour, label, marker in (
            (eas, "#1b5e9e", "hom-derived SNP, EAS history", "o"),
            (chb, "#7b1fa2", "hom-derived SNP, 5R19 history", "s"),
            (tract, "#c62828", "hom-introgressed tract, 5R19 history", "^"),
        ):
            part = frame[frame["threshold_years"] == years].sort_values("band_low")
            axis.plot(
                part["band_low"],
                part["ratio"],
                marker=marker,
                color=colour,
                lw=1.4,
                ms=5,
                label=label,
            )
        axis.axhline(1.0, color="black", lw=0.9, ls=":")
        axis.set_title(f"P(TMRCA < {years:,.0f} y)")
        axis.set_xlabel("frequency of the conditioning class")
    axes[0].set_ylabel("null inflation  (conditioned / unconditional)")
    axes[0].legend(fontsize=8)
    fig.suptitle(
        "Conditioning on the ascertainment shifts the null under strict neutrality",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(OUT / "run6_conditioning_inflation.png", dpi=150)
    plt.close(fig)


def figure_power() -> None:
    summary = pd.read_csv(
        "sim_results_run6/null_comparison/null_comparison_summary.tsv", sep="\t"
    )
    order = [
        "chb_unconditional",
        "eas_unconditional",
        "snp_hom_matched",
        "introgressed_hom_matched",
    ]
    colours = ["#999999", "#1b5e9e", "#7b1fa2", "#c62828"]
    thresholds = sorted(summary["threshold_years"].unique())
    x = np.arange(len(thresholds))
    width = 0.2

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.4))
    for offset, (name, colour) in enumerate(zip(order, colours)):
        part = summary[summary["null"] == name].sort_values("threshold_years")
        axes[0].bar(
            x + (offset - 1.5) * width, part["mean_null"], width, color=colour, label=name
        )
        axes[1].bar(
            x + (offset - 1.5) * width, part["power"], width, color=colour, label=name
        )
    observed = (
        summary[summary["null"] == "chb_unconditional"]
        .sort_values("threshold_years")["mean_observed"]
        .to_numpy()
    )
    axes[0].plot(x, observed, "k*--", ms=12, lw=1.2, label="observed (selected)")
    axes[0].set_ylabel("mean P(TMRCA < x)")
    axes[0].set_title("Where each candidate null sits")
    axes[1].axhline(0.05, color="black", lw=0.9, ls=":")
    axes[1].set_ylabel(r"power at $\alpha$ = 0.05")
    axes[1].set_title("Power for the CHB hom-carrier statistic")
    for axis in axes:
        axis.set_xticks(x)
        axis.set_xticklabels([f"{t:,.0f} y" for t in thresholds])
        axis.legend(fontsize=8)
    fig.suptitle(
        "The same observation scored against four nulls "
        "(introgression pulse ≈ 36,000 y)",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(OUT / "run6_null_comparison.png", dpi=150)
    plt.close(fig)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    figure_scan()
    figure_conditioning()
    figure_power()
    print("figures written to", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
