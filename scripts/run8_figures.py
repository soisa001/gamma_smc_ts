#!/usr/bin/env python
"""Figures for run8's three checks."""

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

from gamma_smc_aou.run8_analysis import CLASSES, load_arm  # noqa: E402
from gamma_smc_aou.run8_config import (  # noqa: E402
    ALPHA,
    FOCAL_POSITION_BP,
    SCENARIOS,
    THRESHOLDS_YEARS,
)

STUDY = REPO / "sim_results_run8"
OUT = STUDY / "figures"
TABLES = STUDY / "tables"
ORDER = [s.scenario_id for s in SCENARIOS]
LABEL = {
    "recent_denovo": "recent de novo\n5 kya, s=0.1",
    "postintro_denovo": "de novo at pulse\n55 kya, s=0.002",
    "introgressed_pulse": "introgressed pulse\n55 kya, s=0.002",
    "introgressed_pulse_weak": "introgressed pulse\n55 kya, s=0.001",
}
COLOUR = {
    "recent_denovo": "#c1121f",
    "postintro_denovo": "#fb8500",
    "introgressed_pulse": "#219ebc",
    "introgressed_pulse_weak": "#8ecae6",
}
CLASS_COLOUR = {
    "carrier_carrier": "#1b5e9e",
    "carrier_noncarrier": "#b8860b",
    "noncarrier_noncarrier": "#777777",
    "overall": "#2e7d32",
}


def figure_scan(scan: pd.DataFrame) -> None:
    """Check 1: does the unconditional statistic spike locally?"""
    fig, axes = plt.subplots(2, len(ORDER), figsize=(18, 7.5))
    for column, arm in enumerate(ORDER):
        cutoff = next(s.expected_cutoff_years for s in SCENARIOS if s.scenario_id == arm)
        for row, method in enumerate(("truth", "decode")):
            axis = axes[row, column]
            for source, style in ((arm, "-"), ("neutral", ":")):
                part = scan[
                    (scan["arm"] == source)
                    & (scan["method"] == method)
                    & (scan["genotype_class"] == "overall")
                    & (scan["threshold_years"] == cutoff)
                ].sort_values("position_0based")
                if part.empty:
                    continue
                axis.plot(
                    (part["position_0based"] - FOCAL_POSITION_BP) / 1e6,
                    part["mean_statistic"],
                    style,
                    lw=1.4,
                    color=COLOUR[arm] if source == arm else "black",
                    label="selected" if source == arm else "neutral",
                )
            axis.axvline(0.0, color="red", ls="--", lw=0.8, alpha=0.5)
            axis.set_xlim(-2.5, 2.5)
            if row == 0:
                axis.set_title(f"{LABEL[arm]}\ncutoff {cutoff:,.0f} y", fontsize=9)
            else:
                axis.set_xlabel("distance from focal site (Mb)")
        axes[0, column].legend(fontsize=7)
    axes[0, 0].set_ylabel("tree truth\nP(TMRCA < x)")
    axes[1, 0].set_ylabel("Gamma-SMC decode\nfrac. posterior means < x")
    fig.suptitle(
        "Check 1: unconditional statistic along the contig, no carrier conditioning",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(OUT / "run8_check1_scan.png", dpi=150)
    plt.close(fig)


def figure_classes(scan: pd.DataFrame) -> None:
    """Check 2: stratifying pairs by carrier status."""
    fig, axes = plt.subplots(1, len(ORDER), figsize=(18, 4.2), sharey=True)
    for column, arm in enumerate(ORDER):
        cutoff = next(s.expected_cutoff_years for s in SCENARIOS if s.scenario_id == arm)
        axis = axes[column]
        for name in CLASSES:
            part = scan[
                (scan["arm"] == arm)
                & (scan["method"] == "truth")
                & (scan["genotype_class"] == name)
                & (scan["threshold_years"] == cutoff)
            ].sort_values("position_0based")
            if part.empty:
                continue
            axis.plot(
                (part["position_0based"] - FOCAL_POSITION_BP) / 1e6,
                part["mean_statistic"],
                lw=1.4,
                color=CLASS_COLOUR[name],
                label=name.replace("_", "–"),
            )
        neutral = scan[
            (scan["arm"] == "neutral")
            & (scan["method"] == "truth")
            & (scan["genotype_class"] == "overall")
            & (scan["threshold_years"] == cutoff)
        ].sort_values("position_0based")
        if not neutral.empty:
            axis.plot(
                (neutral["position_0based"] - FOCAL_POSITION_BP) / 1e6,
                neutral["mean_statistic"],
                ":", lw=1.4, color="black", label="neutral",
            )
        axis.axvline(0.0, color="red", ls="--", lw=0.8, alpha=0.5)
        axis.set_xlim(-2.5, 2.5)
        axis.set_title(f"{LABEL[arm]}\ncutoff {cutoff:,.0f} y", fontsize=9)
        axis.set_xlabel("distance from focal site (Mb)")
    axes[0].set_ylabel("P(TMRCA < x), tree truth")
    axes[0].legend(fontsize=7)
    fig.suptitle(
        "Check 2: carrier-pair stratification — the spike is carried entirely by "
        "carrier–carrier pairs",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(OUT / "run8_check2_classes.png", dpi=150)
    plt.close(fig)


def figure_power(exceed: pd.DataFrame) -> None:
    """Check 3: AUC and power against the unconditioned null."""
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharex=True)
    overall = exceed[exceed["genotype_class"] == "overall"]
    for column, method in enumerate(("truth", "decode")):
        block = overall[overall["method"] == method]
        for arm in ORDER:
            part = block[block["arm"] == arm].sort_values("threshold_years")
            if part.empty:
                continue
            axes[0, column].plot(
                part["threshold_years"], part["auc"], marker="o", ms=4, lw=1.5,
                color=COLOUR[arm], label=LABEL[arm].replace("\n", " "),
            )
            axes[1, column].plot(
                part["threshold_years"], part["power_rule"], marker="o", ms=4,
                lw=1.5, color=COLOUR[arm],
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
    axes[1, 0].set_ylabel("power (≤1 of 100 neutral exceeds)")
    axes[0, 0].set_ylim(0, 1.02)
    axes[1, 0].set_ylim(0, 1.02)
    axes[0, 0].legend(fontsize=8, loc="lower left")
    fig.suptitle(
        "Check 3: AUC and power against an unconditioned null built from the "
        "overall statistic",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(OUT / "run8_check3_power.png", dpi=150)
    plt.close(fig)


def figure_pvalues(focal: pd.DataFrame) -> None:
    """Check 3: p-value distributions, selected against the neutral null."""
    neutral = focal[focal["arm"] == "neutral"]
    fig, axes = plt.subplots(2, len(ORDER), figsize=(18, 7.5), sharex=True)
    edges = np.linspace(0, 1, 21)
    for column, arm in enumerate(ORDER):
        cutoff = next(s.expected_cutoff_years for s in SCENARIOS if s.scenario_id == arm)
        for row, method in enumerate(("truth", "decode")):
            axis = axes[row, column]
            null = neutral[
                (neutral["method"] == method)
                & (neutral["genotype_class"] == "overall")
                & (neutral["threshold_years"] == cutoff)
            ]["statistic"].to_numpy(dtype=float)
            for source, colour, label in (
                (arm, COLOUR[arm], "selected"),
                ("neutral", "grey", "neutral (null)"),
            ):
                values = focal[
                    (focal["arm"] == source)
                    & (focal["method"] == method)
                    & (focal["genotype_class"] == "overall")
                    & (focal["threshold_years"] == cutoff)
                ]["statistic"].to_numpy(dtype=float)
                if values.size == 0 or null.size == 0:
                    continue
                p = np.array([(np.sum(null >= v) + 1.0) / (null.size + 1.0) for v in values])
                axis.hist(
                    p, bins=edges, density=True, histtype="stepfilled" if source == "neutral" else "step",
                    lw=1.6, alpha=0.45 if source == "neutral" else 1.0,
                    color=colour, label=label,
                )
            axis.axvline(ALPHA, color="red", ls="--", lw=1.2)
            axis.axhline(1.0, color="black", ls=":", lw=0.9)
            if row == 0:
                axis.set_title(f"{LABEL[arm]}\ncutoff {cutoff:,.0f} y", fontsize=9)
            else:
                axis.set_xlabel("p-value")
        axes[0, column].legend(fontsize=7)
    axes[0, 0].set_ylabel("density (tree truth)")
    axes[1, 0].set_ylabel("density (decode)")
    fig.suptitle(
        "Check 3: p-value distributions at each scenario's expected cutoff; the "
        "neutral arm should be flat",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(OUT / "run8_check3_pvalues.png", dpi=150)
    plt.close(fig)


def median_replicate_report(focal: pd.DataFrame) -> pd.DataFrame:
    """The p-value for the AF-median selected replicate, as the study asked."""
    neutral = focal[focal["arm"] == "neutral"]
    rows = []
    for scenario in SCENARIOS:
        arm = scenario.scenario_id
        block = focal[(focal["arm"] == arm) & (focal["genotype_class"] == "overall")]
        if block.empty:
            continue
        # One replicate per arm: the one whose final frequency is the median.
        per_replicate = block.groupby("replicate_index")["sample_af"].first()
        target = per_replicate.iloc[
            int(np.argmin(np.abs(per_replicate.to_numpy() - per_replicate.median())))
        ]
        chosen = int(per_replicate.index[
            int(np.argmin(np.abs(per_replicate.to_numpy() - per_replicate.median())))
        ])
        for method in ("truth", "decode"):
            for years in THRESHOLDS_YEARS:
                null = neutral[
                    (neutral["method"] == method)
                    & (neutral["genotype_class"] == "overall")
                    & (neutral["threshold_years"] == years)
                ]["statistic"].to_numpy(dtype=float)
                value = block[
                    (block["replicate_index"] == chosen)
                    & (block["method"] == method)
                    & (block["threshold_years"] == years)
                ]["statistic"]
                if value.empty or null.size == 0:
                    continue
                observed = float(value.iloc[0])
                exceed = int(np.sum(null >= observed))
                rows.append(
                    {
                        "arm": arm,
                        "replicate_index": chosen,
                        "sample_af": float(target),
                        "method": method,
                        "threshold_years": float(years),
                        "statistic": observed,
                        "null_mean": float(null.mean()),
                        "n_neutral_exceeding": exceed,
                        "p_corrected": (exceed + 1.0) / (null.size + 1.0),
                        "significant_rule": bool(exceed <= 1),
                    }
                )
    return pd.DataFrame(rows)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    focal = pd.read_csv(TABLES / "focal_table.tsv", sep="\t")
    exceed = pd.read_csv(TABLES / "exceedance.tsv", sep="\t")
    scan = pd.read_csv(TABLES / "scan.tsv", sep="\t")

    figure_scan(scan)
    figure_classes(scan)
    figure_power(exceed)
    figure_pvalues(focal)

    median = median_replicate_report(focal)
    median.to_csv(TABLES / "median_replicate.tsv", sep="\t", index=False)

    pd.set_option("display.width", 250)
    print("=== AF-median selected replicate, p-values vs the neutral null ===")
    for method in ("truth", "decode"):
        print(f"\n-- {method} --")
        block = median[median["method"] == method]
        print(
            block.pivot_table(
                index="threshold_years", columns="arm", values="p_corrected"
            ).reindex(columns=ORDER).to_string(float_format=lambda x: f"{x:.4f}")
        )
    print("\nfigures ->", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
