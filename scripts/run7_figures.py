#!/usr/bin/env python
"""The run7 plot suite: final AF, P(TMRCA<x), AUC, p-values, power and the scan."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "python"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from gamma_smc_aou.run7_analysis import (  # noqa: E402
    ALPHA,
    DENOVO_COEFFICIENTS,
    focal_table,
    scan_table,
    summarise,
)
from gamma_smc_aou.run7_config import (  # noqa: E402
    FOCAL_POSITION_BP,
    GENERATION_TIME_YEARS,
    PULSE_GENERATIONS,
    SELECTION_COEFFICIENTS,
    THRESHOLDS_YEARS,
)

STUDY = REPO / "sim_results_run7"
OUT = STUDY / "figures"
PULSE_YEARS = PULSE_GENERATIONS * GENERATION_TIME_YEARS
COLOURS = {0.002: "#8ecae6", 0.003: "#219ebc", 0.005: "#fb8500", 0.01: "#c1121f"}
DENOVO_COLOURS = {0.005: "#b5e48c", 0.01: "#52b788", 0.02: "#2d6a4f", 0.05: "#1b4332"}
CLASS_LABEL = {"unconditional": "all individuals", "hom_carrier": "hom carriers"}


def _pulse_marker(axis, orientation="v"):
    if orientation == "v":
        axis.axvline(
            PULSE_YEARS, color="grey", ls="--", lw=1.0, alpha=0.8, zorder=0
        )


def figure_final_af(focal: pd.DataFrame) -> None:
    per_replicate = (
        focal[(focal["genotype_class"] == "unconditional")
              & (focal["origin"] == "introgressed")]
        .groupby(["selection_coefficient", "replicate_index"])["sample_af"]
        .first()
        .reset_index()
    )
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
    data = [
        per_replicate.loc[
            per_replicate["selection_coefficient"] == s, "sample_af"
        ].to_numpy()
        for s in SELECTION_COEFFICIENTS
    ]
    parts = axes[0].violinplot(data, showmedians=True, widths=0.8)
    for body, s in zip(parts["bodies"], SELECTION_COEFFICIENTS):
        body.set_facecolor(COLOURS[s])
        body.set_alpha(0.75)
    axes[0].set_xticks(range(1, len(SELECTION_COEFFICIENTS) + 1))
    axes[0].set_xticklabels([f"s = {s:g}" for s in SELECTION_COEFFICIENTS])
    axes[0].set_ylabel("final allele frequency")
    axes[0].set_ylim(0, 1.02)
    axes[0].axhline(2 / 3, color="black", ls=":", lw=1.0)
    axes[0].text(
        0.02, 2 / 3 + 0.02, "f = 2/3 crossover", fontsize=8, transform=axes[0].get_yaxis_transform()
    )
    axes[0].set_title("Where the sweep ends up")

    for s in SELECTION_COEFFICIENTS:
        values = per_replicate.loc[
            per_replicate["selection_coefficient"] == s, "sample_af"
        ].to_numpy()
        axes[1].hist(
            values, bins=np.linspace(0, 1, 26), histtype="step", lw=1.6,
            color=COLOURS[s], label=f"s = {s:g} (mean {values.mean():.2f})",
        )
    axes[1].set_xlabel("final allele frequency")
    axes[1].set_ylabel("replicates")
    axes[1].legend(fontsize=8)
    axes[1].set_title("Distribution across replicates")
    fig.suptitle(
        "run7: introgressed allele from 2.5% at the 55 kya pulse, conditioned on survival",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(OUT / "run7_final_af.png", dpi=150)
    plt.close(fig)


def figure_ptmrca(summary: pd.DataFrame, focal: pd.DataFrame) -> None:
    """P(TMRCA<x) against the neutral band, per class."""
    neutral = focal[focal["genotype_class"] == "neutral"]
    bands = neutral.groupby("threshold_years")["statistic"].agg(
        median="median",
        low=lambda s: float(np.quantile(s, 0.025)),
        high=lambda s: float(np.quantile(s, 0.975)),
    ).reset_index()

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6), sharey=True)
    for axis, name in zip(axes, ("unconditional", "hom_carrier")):
        axis.fill_between(
            bands["threshold_years"], bands["low"], bands["high"],
            color="grey", alpha=0.25, label="neutral 95% band",
        )
        axis.plot(
            bands["threshold_years"], bands["median"], color="black", lw=1.6,
            ls="--", label="neutral median",
        )
        for s in SELECTION_COEFFICIENTS:
            block = summary[
                (summary["origin"] == "introgressed")
                & (summary["selection_coefficient"] == s)
                & (summary["genotype_class"] == name)
            ].sort_values("threshold_years")
            if block.empty:
                continue
            axis.plot(
                block["threshold_years"], block["mean_statistic"],
                marker="o", ms=4, lw=1.6, color=COLOURS[s], label=f"s = {s:g}",
            )
        _pulse_marker(axis)
        axis.set_xscale("log")
        axis.set_xlabel("TMRCA cutoff x (years)")
        axis.set_title(CLASS_LABEL[name])
    axes[0].set_ylabel("P(TMRCA < x)")
    axes[0].legend(fontsize=8, loc="upper left")
    axes[1].text(
        PULSE_YEARS * 1.05, 0.05, "pulse", fontsize=8, color="grey", rotation=90
    )
    fig.suptitle("P(TMRCA < x) at the focal site, against the neutral band", fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT / "run7_ptmrca.png", dpi=150)
    plt.close(fig)


def figure_auc_and_power(summary: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13, 8.4), sharex=True)
    for column, name in enumerate(("unconditional", "hom_carrier")):
        for s in SELECTION_COEFFICIENTS:
            block = summary[
                (summary["origin"] == "introgressed")
                & (summary["selection_coefficient"] == s)
                & (summary["genotype_class"] == name)
            ].sort_values("threshold_years")
            if block.empty:
                continue
            axes[0, column].plot(
                block["threshold_years"], block["auc"], marker="o", ms=4,
                lw=1.6, color=COLOURS[s], label=f"s = {s:g}",
            )
            axes[1, column].plot(
                block["threshold_years"], block["power"], marker="o", ms=4,
                lw=1.6, color=COLOURS[s], label=f"s = {s:g}",
            )
        axes[0, column].axhline(0.5, color="black", ls=":", lw=1.0)
        axes[1, column].axhline(ALPHA, color="black", ls=":", lw=1.0)
        for row in (0, 1):
            _pulse_marker(axes[row, column])
            axes[row, column].set_xscale("log")
        axes[0, column].set_title(CLASS_LABEL[name])
        axes[1, column].set_xlabel("TMRCA cutoff x (years)")
    axes[0, 0].set_ylabel("AUC vs neutral")
    axes[1, 0].set_ylabel(rf"power at $\alpha$ = {ALPHA}")
    axes[0, 0].set_ylim(0, 1.02)
    axes[1, 0].set_ylim(0, 1.02)
    axes[0, 0].legend(fontsize=8, loc="upper left")
    fig.suptitle(
        "Detection is confined to cutoffs below the pulse (dashed grey, 55 kya)",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(OUT / "run7_auc_power.png", dpi=150)
    plt.close(fig)


def figure_auc_vs_af(focal: pd.DataFrame) -> None:
    """The relationship the whole design turns on."""
    thresholds = [20_000.0, 30_000.0, 50_000.0]
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharex=True, sharey=True)
    for row, name in enumerate(("unconditional", "hom_carrier")):
        for column, years in enumerate(thresholds):
            axis = axes[row, column]
            block = focal[
                (focal["genotype_class"] == name)
                & (focal["threshold_years"] == years)
                & (focal["origin"] == "introgressed")
            ]
            for s in SELECTION_COEFFICIENTS:
                part = block[(block["selection_coefficient"] == s)
                             & (block["origin"] == "introgressed")]
                axis.scatter(
                    part["sample_af"], part["percentile"], s=16, alpha=0.65,
                    color=COLOURS[s], label=f"s = {s:g}" if (row == 0 and column == 0) else None,
                )
            if len(block) > 10:
                order = np.argsort(block["sample_af"].to_numpy())
                af = block["sample_af"].to_numpy()[order]
                auc = block["percentile"].to_numpy()[order]
                width = max(len(af) // 10, 5)
                smooth = pd.Series(auc).rolling(width, center=True, min_periods=3).mean()
                axis.plot(af, smooth, color="black", lw=1.8)
            axis.axhline(0.5, color="black", ls=":", lw=1.0)
            axis.axvline(2 / 3, color="grey", ls="--", lw=1.0)
            if row == 0:
                axis.set_title(f"P(TMRCA < {years:,.0f} y)")
            if row == 1:
                axis.set_xlabel("final allele frequency")
        axes[row, 0].set_ylabel(f"AUC vs neutral\n({CLASS_LABEL[name]})")
    axes[0, 0].legend(fontsize=8, loc="lower right")
    fig.suptitle(
        "AUC against final allele frequency; grey line marks the f = 2/3 crossover",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(OUT / "run7_auc_vs_af.png", dpi=150)
    plt.close(fig)


def figure_pvalues(focal: pd.DataFrame) -> None:
    thresholds = [20_000.0, 30_000.0, 50_000.0]
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharex=True)
    edges = np.linspace(0, 1, 21)
    for row, name in enumerate(("unconditional", "hom_carrier")):
        for column, years in enumerate(thresholds):
            axis = axes[row, column]
            neutral = focal[
                (focal["genotype_class"] == "neutral")
                & (focal["threshold_years"] == years)
            ]["p_value"]
            axis.hist(
                neutral, bins=edges, density=True, color="grey", alpha=0.45,
                label="neutral (null)",
            )
            for s in SELECTION_COEFFICIENTS:
                part = focal[
                    (focal["genotype_class"] == name)
                    & (focal["threshold_years"] == years)
                    & (focal["origin"] == "introgressed")
                    & (focal["selection_coefficient"] == s)
                ]["p_value"]
                if part.empty:
                    continue
                axis.hist(
                    part, bins=edges, density=True, histtype="step", lw=1.6,
                    color=COLOURS[s], label=f"s = {s:g}",
                )
            axis.axvline(ALPHA, color="red", ls="--", lw=1.2)
            axis.axhline(1.0, color="black", ls=":", lw=1.0)
            if row == 0:
                axis.set_title(f"P(TMRCA < {years:,.0f} y)")
            if row == 1:
                axis.set_xlabel("p-value")
        axes[row, 0].set_ylabel(f"density\n({CLASS_LABEL[name]})")
    axes[0, 0].legend(fontsize=8)
    fig.suptitle(
        "p-value distributions; the neutral arm should be flat at 1 if the null is calibrated",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(OUT / "run7_pvalues.png", dpi=150)
    plt.close(fig)


def figure_scan(scan: pd.DataFrame) -> None:
    thresholds = sorted(scan["threshold_years"].unique())
    fig, axes = plt.subplots(2, len(thresholds), figsize=(15, 8), sharey="row")
    for column, years in enumerate(thresholds):
        for row, name in enumerate(("unconditional", "hom_carrier")):
            axis = axes[row, column]
            for s in SELECTION_COEFFICIENTS:
                part = scan[
                    (scan["threshold_years"] == years)
                    & (scan["genotype_class"] == name)
                    & (scan["origin"] == "introgressed")
                    & (scan["selection_coefficient"] == s)
                ].sort_values("position_0based")
                if part.empty:
                    continue
                axis.plot(
                    (part["position_0based"] - FOCAL_POSITION_BP) / 1e6,
                    part["auc_vs_neutral"], lw=1.3, color=COLOURS[s],
                    label=f"s = {s:g}" if (row == 0 and column == 0) else None,
                )
            axis.axhline(0.5, color="black", ls=":", lw=0.9)
            axis.axvline(0.0, color="red", ls="--", lw=0.9, alpha=0.6)
            axis.set_xlim(-2.5, 2.5)
            if row == 0:
                axis.set_title(f"P(TMRCA < {years:,.0f} y)")
            if row == 1:
                axis.set_xlabel("distance from focal site (Mb)")
        axes[row, 0].set_ylabel(f"AUC vs neutral\n({CLASS_LABEL[name]})")
    axes[0, 0].set_ylim(0, 1.02)
    axes[0, 0].legend(fontsize=8, loc="upper left")
    fig.suptitle("Spatial extent of the signal along the 10 Mb contig", fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT / "run7_scan.png", dpi=150)
    plt.close(fig)



def figure_origin_contrast(summary: pd.DataFrame) -> None:
    """De novo against introgressed: the same sweep from a different origin.

    The introgressed arms are floored by the pulse -- their carriers descend from
    many archaic haplotypes whose mutual coalescence sits in the archaic branch,
    so nothing about them can be recent. A de novo sweep has no such floor, and
    the gap at recent cutoffs is the whole difference between the two.
    """
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), sharey=True)
    for axis, name in zip(axes, ("unconditional", "hom_carrier")):
        for s in SELECTION_COEFFICIENTS:
            block = summary[
                (summary["origin"] == "introgressed")
                & (summary["selection_coefficient"] == s)
                & (summary["genotype_class"] == name)
            ].sort_values("threshold_years")
            if block.empty:
                continue
            axis.plot(
                block["threshold_years"], block["auc"], marker="o", ms=4, lw=1.6,
                color=COLOURS[s], label=f"introgressed, s = {s:g}",
            )
        for s in DENOVO_COEFFICIENTS:
            block = summary[
                (summary["origin"] == "de_novo")
                & (summary["selection_coefficient"] == s)
                & (summary["genotype_class"] == name)
            ].sort_values("threshold_years")
            if block.empty:
                continue
            axis.plot(
                block["threshold_years"], block["auc"], marker="s", ms=4, lw=1.6,
                ls="--", color=DENOVO_COLOURS[s], label=f"de novo, s = {s:g}",
            )
        axis.axhline(0.5, color="black", ls=":", lw=1.0)
        _pulse_marker(axis)
        axis.set_xscale("log")
        axis.set_xlabel("TMRCA cutoff x (years)")
        axis.set_title(CLASS_LABEL[name])
    axes[0].set_ylabel("AUC vs neutral")
    axes[0].set_ylim(0, 1.02)
    axes[0].legend(fontsize=7, loc="upper left", ncol=2)
    fig.suptitle(
        "A de novo sweep is detectable at cutoffs adaptive introgression cannot reach",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(OUT / "run7_origin_contrast.png", dpi=150)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study-root", default=str(STUDY))
    args = parser.parse_args(argv)

    OUT.mkdir(parents=True, exist_ok=True)
    tables = STUDY / "tables"
    tables.mkdir(parents=True, exist_ok=True)

    focal = focal_table(args.study_root)
    focal.to_csv(tables / "focal_table.tsv", sep="\t", index=False)
    summary = summarise(focal)
    summary.to_csv(tables / "summary.tsv", sep="\t", index=False)
    scan = scan_table(args.study_root)
    scan.to_csv(tables / "scan.tsv", sep="\t", index=False)

    figure_final_af(focal)
    figure_ptmrca(summary, focal)
    figure_auc_and_power(summary)
    figure_auc_vs_af(focal)
    figure_pvalues(focal)
    figure_scan(scan)
    figure_origin_contrast(summary)

    pd.set_option("display.width", 250)
    view = summary[
        summary["threshold_years"].isin(
            [4_500.0, 10_000.0, 20_000.0, 30_000.0, 50_000.0]
        )
    ]
    for origin in ("introgressed", "de_novo"):
        print(f"=== {origin} ===")
        print(
            view[view["origin"] == origin]
            .pivot_table(
                index=["genotype_class", "threshold_years"],
                columns="selection_coefficient",
                values=["auc", "power"],
            )
            .to_string(float_format=lambda x: f"{x:.3f}")
        )
        print()
    print()
    print("figures ->", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
