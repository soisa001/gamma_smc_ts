#!/usr/bin/env python
"""Figures for run9."""

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

from gamma_smc_aou.run9_config import (  # noqa: E402
    ALPHA,
    FOCAL_POSITION_BP,
    N_NEUTRAL,
    SCENARIOS,
)

STUDY = REPO / "sim_results_run9"
OUT = STUDY / "figures"
TABLES = STUDY / "tables"

PULSE = [s.scenario_id for s in SCENARIOS if s.origin == "introgressed_pulse"]
DENOVO = [s.scenario_id for s in SCENARIOS if s.origin == "de_novo"]
ORDER = [s.scenario_id for s in SCENARIOS]
FLOOR = 1.0 / (N_NEUTRAL + 1.0)

COLOUR = {
    "pulse_s0p003": "#8ecae6",
    "pulse_s0p004": "#219ebc",
    "pulse_s0p005": "#023047",
    "pulse_delayed_s0p005": "#c1121f",
    "denovo_s0p003": "#b5e48c",
    "denovo_s0p004": "#52b788",
    "denovo_s0p005": "#1b4332",
    "denovo_delayed_s0p01": "#ffb703",
    "denovo_delayed_s0p05": "#fb8500",
}
STYLE = {s.scenario_id: ("--" if s.is_delayed else "-") for s in SCENARIOS}
CLASS_COLOUR = {
    "carrier_carrier": "#1b5e9e",
    "carrier_noncarrier": "#b8860b",
    "noncarrier_noncarrier": "#777777",
    "overall": "#2e7d32",
}


def figure_power(exceed: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(14, 8.5), sharex=True)
    overall = exceed[exceed["genotype_class"] == "overall"]
    for column, method in enumerate(("truth", "decode")):
        block = overall[overall["method"] == method]
        for arm in ORDER:
            part = block[block["arm"] == arm].sort_values("threshold_years")
            if part.empty:
                continue
            for row, metric in enumerate(("auc", "power_rule")):
                axes[row, column].plot(
                    part["threshold_years"], part[metric], STYLE[arm],
                    marker="o", ms=3.5, lw=1.5, color=COLOUR[arm],
                    label=arm if row == 0 else None,
                )
        axes[0, column].axhline(0.5, color="black", ls=":", lw=0.9)
        axes[1, column].axhline(ALPHA, color="black", ls=":", lw=0.9)
        axes[0, column].set_title(
            "tree truth" if method == "truth" else "Gamma-SMC decode", fontsize=11
        )
        axes[1, column].set_xlabel("TMRCA cutoff x (years)")
        for row in (0, 1):
            axes[row, column].set_xscale("log")
    axes[0, 0].set_ylabel("AUC vs neutral")
    axes[1, 0].set_ylabel(f"power (≤1 of {N_NEUTRAL} neutral exceeds)")
    axes[0, 0].set_ylim(0, 1.02)
    axes[1, 0].set_ylim(0, 1.02)
    axes[0, 0].legend(fontsize=6.5, loc="lower left", ncol=2)
    fig.suptitle(
        "run9: unconditional AUC and power. Solid = onset at the 55 kya pulse, "
        "dashed = onset at 10 kya",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(OUT / "run9_power.png", dpi=150)
    plt.close(fig)


def figure_af_threshold(strat: pd.DataFrame) -> None:
    """The portable result: detection as a function of final allele frequency."""
    thresholds = [10_000.0, 30_000.0, 50_000.0]
    fig, axes = plt.subplots(2, len(thresholds), figsize=(15, 8), sharey="row")
    for column, years in enumerate(thresholds):
        for row, method in enumerate(("truth", "decode")):
            axis = axes[row, column]
            part = strat[
                (strat["genotype_class"] == "overall")
                & (strat["method"] == method)
                & (strat["threshold_years"] == years)
            ].sort_values("af_mid")
            if part.empty:
                continue
            axis.plot(part["af_mid"], part["power_rule"], marker="o", ms=5,
                      lw=1.8, color="#c1121f", label="power")
            axis.plot(part["af_mid"], part["auc"], marker="s", ms=5, lw=1.8,
                      ls="--", color="#023047", label="AUC")
            for _, r in part.iterrows():
                axis.annotate(f"{int(r['n_replicates'])}", (r["af_mid"], r["power_rule"]),
                              textcoords="offset points", xytext=(0, 7),
                              ha="center", fontsize=6, color="grey")
            axis.axhline(0.5, color="black", ls=":", lw=0.8)
            axis.axhline(ALPHA, color="grey", ls=":", lw=0.8)
            axis.set_ylim(0, 1.02)
            if row == 0:
                axis.set_title(f"P(TMRCA < {years:,.0f} y)", fontsize=10)
            else:
                axis.set_xlabel("final allele frequency")
        axes[0, column].legend(fontsize=8)
    axes[0, 0].set_ylabel("tree truth")
    axes[1, 0].set_ylabel("Gamma-SMC decode")
    fig.suptitle(
        "Detection against final allele frequency, pooled over all nine arms "
        "(grey numbers = replicates per bin)",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(OUT / "run9_af_threshold.png", dpi=150)
    plt.close(fig)


def figure_quantiles(quant: pd.DataFrame) -> None:
    """Quantile curves: is the median replicate simply weak?"""
    fig, axes = plt.subplots(2, 5, figsize=(19, 7.5), sharey=True)
    shown = PULSE + DENOVO
    for index, arm in enumerate(shown[:10]):
        axis = axes[index // 5, index % 5]
        part = quant[
            (quant["arm"] == arm)
            & (quant["method"] == "truth")
            & (quant["genotype_class"] == "overall")
        ]
        if part.empty:
            continue
        for q, group in part.groupby("quantile"):
            group = group.sort_values("threshold_years")
            axis.plot(
                group["threshold_years"], group["p_corrected"], marker="o", ms=3,
                lw=1.3, label=f"q{int(q * 100)}",
            )
        axis.axhline(ALPHA, color="red", ls="--", lw=1.0)
        axis.axhline(FLOOR, color="grey", ls=":", lw=1.0)
        axis.set_yscale("log")
        axis.set_xscale("log")
        axis.set_title(arm, fontsize=8)
        if index // 5 == 1:
            axis.set_xlabel("cutoff (years)")
    for index in range(len(shown), 10):
        axes[index // 5, index % 5].axis("off")
    axes[0, 0].set_ylabel("p-value (tree truth)")
    axes[1, 0].set_ylabel("p-value (tree truth)")
    axes[0, 0].legend(fontsize=7)
    fig.suptitle(
        f"Quantile curves: p at each percentile of the statistic. Red = α, "
        f"grey = the 1/{N_NEUTRAL + 1} floor",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(OUT / "run9_quantiles.png", dpi=150)
    plt.close(fig)


def figure_classes(exceed: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 5, figsize=(19, 7.5), sharey=True)
    shown = PULSE + DENOVO
    for index, arm in enumerate(shown[:10]):
        axis = axes[index // 5, index % 5]
        part = exceed[(exceed["arm"] == arm) & (exceed["method"] == "truth")]
        for name, group in part.groupby("genotype_class"):
            group = group.sort_values("threshold_years")
            axis.plot(
                group["threshold_years"], group["auc"], marker="o", ms=3, lw=1.3,
                color=CLASS_COLOUR.get(name, "black"),
                label=name.replace("_", "–") if index == 0 else None,
            )
        axis.axhline(0.5, color="black", ls=":", lw=0.9)
        axis.set_xscale("log")
        axis.set_ylim(0, 1.02)
        axis.set_title(arm, fontsize=8)
        if index // 5 == 1:
            axis.set_xlabel("cutoff (years)")
    for index in range(len(shown), 10):
        axes[index // 5, index % 5].axis("off")
    axes[0, 0].set_ylabel("AUC vs neutral")
    axes[1, 0].set_ylabel("AUC vs neutral")
    axes[0, 0].legend(fontsize=7, loc="lower right")
    fig.suptitle(
        "Carrier-pair stratification by arm (tree truth), scored against the "
        "unconditioned null — an upper bound, not calibrated power",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(OUT / "run9_classes.png", dpi=150)
    plt.close(fig)


def figure_scan(scan: pd.DataFrame) -> None:
    arms = ["pulse_s0p005", "pulse_delayed_s0p005", "denovo_s0p005", "denovo_delayed_s0p05"]
    cutoffs = {"pulse_s0p005": 40_000.0, "pulse_delayed_s0p005": 10_000.0,
               "denovo_s0p005": 30_000.0, "denovo_delayed_s0p05": 5_000.0}
    fig, axes = plt.subplots(2, len(arms), figsize=(17, 7.5))
    for column, arm in enumerate(arms):
        years = cutoffs[arm]
        for row, method in enumerate(("truth", "decode")):
            axis = axes[row, column]
            for source, style, colour in ((arm, "-", COLOUR[arm]), ("neutral", ":", "black")):
                part = scan[
                    (scan["arm"] == source)
                    & (scan["method"] == method)
                    & (scan["genotype_class"] == "overall")
                    & (scan["threshold_years"] == years)
                ].sort_values("position_0based")
                if part.empty:
                    continue
                axis.plot(
                    (part["position_0based"] - FOCAL_POSITION_BP) / 1e6,
                    part["mean_statistic"], style, lw=1.4, color=colour,
                    label="selected" if source == arm else "neutral",
                )
            axis.axvline(0.0, color="red", ls="--", lw=0.8, alpha=0.5)
            axis.set_xlim(-2.5, 2.5)
            if row == 0:
                axis.set_title(f"{arm}\ncutoff {years:,.0f} y", fontsize=9)
            else:
                axis.set_xlabel("distance from focal site (Mb)")
        axes[0, column].legend(fontsize=7)
    axes[0, 0].set_ylabel("tree truth")
    axes[1, 0].set_ylabel("Gamma-SMC decode")
    fig.suptitle("Unconditional statistic along the contig", fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT / "run9_scan.png", dpi=150)
    plt.close(fig)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    exceed = pd.read_csv(TABLES / "exceedance.tsv", sep="\t")
    quant = pd.read_csv(TABLES / "quantiles.tsv", sep="\t")
    strat = pd.read_csv(TABLES / "af_stratified.tsv", sep="\t")
    scan = pd.read_csv(TABLES / "scan.tsv", sep="\t")

    figure_power(exceed)
    figure_af_threshold(strat)
    figure_quantiles(quant)
    figure_classes(exceed)
    figure_scan(scan)

    pd.set_option("display.width", 250)
    print("=== AF-stratified detection (overall class) ===")
    for method in ("truth", "decode"):
        part = strat[
            (strat["genotype_class"] == "overall")
            & (strat["method"] == method)
            & (strat["threshold_years"].isin([10_000.0, 30_000.0, 50_000.0]))
        ]
        if part.empty:
            continue
        print(f"\n-- {method}: power --")
        print(
            part.pivot_table(index="af_bin", columns="threshold_years",
                             values="power_rule").to_string(float_format=lambda x: f"{x:.3f}")
        )
    print("\nfigures ->", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
