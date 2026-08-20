#!/usr/bin/env python
"""Decompose the run5 within-individual TMRCA statistic by genotype class.

An introgressed swept locus produces a *trimodal* within-individual TMRCA
distribution, and the three modes pull in opposite directions:

* **hom-carrier** pairs are two archaic haplotypes, which under a sweep coalesce
  after introgression -- so they are far *younger* than neutral.
* **heterozygous** pairs are one archaic and one modern haplotype, which cannot
  coalesce until before the human-Neanderthal split -- so they are far *older*
  than neutral, and can never be recent.
* **hom-non-carrier** pairs are ordinary CHB pairs, so neutral-like.

``P(TMRCA < x)`` sums over all three. It therefore only shows net signal once the
hom-carrier mass outweighs the heterozygous mass, which at Hardy-Weinberg
proportions happens at ``f^2 > 2f(1-f)``, i.e. **f > 2/3**. Below that the
heterozygous class actively pushes the statistic below neutral, which is why AUC
can fall under 0.5 at low final frequency.

Writes a table and a figure to the arm's results/figures directories.
"""

from __future__ import annotations

import glob
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "python"))

GENERATION_TIME = 29.0
CLASSES = ("hom_carrier", "het", "hom_noncarrier")
BANDS = ((0.0, 0.4), (0.4, 0.7), (0.7, 1.01))


def _load(arm_dir: Path, mode: str):
    records = []
    for path in sorted(glob.glob(str(arm_dir / mode / "replicates" / "*"))):
        directory = Path(path)
        try:
            tmrca = np.load(directory / "focal_tmrca.npy")
            endpoint = json.loads((directory / "endpoint.json").read_text())
        except (OSError, ValueError):
            continue
        if endpoint.get("status") != "completed":
            continue
        genotypes = None
        if (directory / "focal_genotypes.npy").is_file():
            genotypes = np.load(directory / "focal_genotypes.npy")
        census = endpoint["final_allele_frequency"].get("census_af")
        records.append((tmrca, genotypes, float(census) if census is not None else np.nan))
    return records


def main(study_root: str = "sim_results_run5", arm_id: str = "chb_from_introgression") -> int:
    arm_dir = Path(study_root) / arm_id
    selected = _load(arm_dir, "selected")
    neutral = _load(arm_dir, "neutral")
    if not selected or not neutral:
        raise SystemExit(f"no replicates found under {arm_dir}")

    neutral_tmrca = np.concatenate([t for t, _, _ in neutral])
    rows = []
    for tmrca, genotypes, census in selected:
        if genotypes is None:
            continue
        for name, mask in (
            ("hom_carrier", genotypes == 2),
            ("het", genotypes == 1),
            ("hom_noncarrier", genotypes == 0),
        ):
            if not mask.any():
                continue
            rows.append(
                {
                    "census_af": census,
                    "genotype_class": name,
                    "n_pairs": int(mask.sum()),
                    "median_tmrca_generations": float(np.median(tmrca[mask])),
                    "p_lt_30ky": float(np.mean(tmrca[mask] < 30_000 / GENERATION_TIME)),
                }
            )
    frame = pd.DataFrame(rows)

    results = arm_dir / "results"
    results.mkdir(parents=True, exist_ok=True)
    frame.to_csv(results / "genotype_decomposition.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    bins = np.logspace(np.log10(50), np.log10(2e5), 60)
    axes[0].hist(neutral_tmrca, bins=bins, density=True, color="#2166ac", alpha=0.45,
                 label="neutral (all pairs)")
    colours = {"hom_carrier": "#b2182b", "het": "#4d004b", "hom_noncarrier": "#999999"}
    for name in CLASSES:
        pooled = np.concatenate(
            [t[(g == {"hom_carrier": 2, "het": 1, "hom_noncarrier": 0}[name])]
             for t, g, _ in selected if g is not None]
        )
        if pooled.size:
            axes[0].hist(pooled, bins=bins, density=True, histtype="step", lw=2,
                         color=colours[name], label=f"selected: {name}")
    axes[0].axvline(30_000 / GENERATION_TIME, color="black", ls="--", lw=1.1)
    axes[0].set(xscale="log", xlabel="within-individual TMRCA (generations)",
                ylabel="density",
                title="Trimodal by genotype class\n(dashed line = 30 ky cutoff)")
    axes[0].legend(frameon=False, fontsize=8)

    for name in CLASSES:
        xs, ys = [], []
        for lo, hi in BANDS:
            band = frame[(frame["census_af"] >= lo) & (frame["census_af"] < hi)
                         & (frame["genotype_class"] == name)]
            if band.empty:
                continue
            xs.append((lo + min(hi, 1.0)) / 2)
            ys.append(np.average(band["p_lt_30ky"], weights=band["n_pairs"]))
        if xs:
            axes[1].plot(xs, ys, marker="o", lw=2, color=colours[name], label=name)
    axes[1].axhline(float(np.mean(neutral_tmrca < 30_000 / GENERATION_TIME)),
                    color="#2166ac", ls="--", lw=1.2, label="neutral baseline")
    axes[1].axvline(2 / 3, color="black", ls=":", lw=1.2)
    axes[1].set(xlabel="final census allele frequency", ylabel="P(TMRCA < 30 ky)",
                ylim=(-0.03, 1.03),
                title="Dotted line: f = 2/3, where hom-carrier\npairs start to outnumber heterozygous ones")
    axes[1].legend(frameon=False, fontsize=8)
    fig.suptitle(f"{arm_id} - why P(TMRCA < x) depends on final frequency", y=1.02)

    figures = arm_dir / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "pdf"):
        fig.savefig(figures / f"genotype_decomposition.{suffix}", dpi=170,
                    bbox_inches="tight")
    plt.close(fig)

    print(f"neutral: median {np.median(neutral_tmrca):.0f} gen, "
          f"P(<30ky) {np.mean(neutral_tmrca < 30_000 / GENERATION_TIME):.3f}")
    for name in CLASSES:
        subset = frame[frame["genotype_class"] == name]
        weights = subset["n_pairs"]
        print(f"{name:15s} pairs={int(weights.sum()):6d}  "
              f"median {np.average(subset['median_tmrca_generations'], weights=weights):8.0f} gen  "
              f"P(<30ky) {np.average(subset['p_lt_30ky'], weights=weights):.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(*sys.argv[1:]))
