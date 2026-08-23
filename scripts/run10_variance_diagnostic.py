#!/usr/bin/env python
"""Is the null's width sampling noise or genealogical variance?

This decides whether a bigger panel can rescue detection at realistic
introgression frequencies. The observed per-window statistic varies for two
reasons:

    Var(observed) = Var_between_replicates(true local P)  +  E[sampling variance]

The second term shrinks as more haplotypes are sampled at a locus; the first does
not, because it is real variation in the local genealogy from one replicate to
the next. If the genealogical term dominates, then going from 200 haplotypes to
4,000 buys almost nothing and the AF ~ 0.6 detection threshold is a property of
the statistic rather than of our sample size.

The test subsamples run9's *neutral* replicates two ways and watches the variance:

* by **pair count** -- 1,000 / 3,000 / 10,000 of the drawn pairs;
* by **haplotype count** -- 50 / 100 / 200 haplotypes, keeping only pairs whose
  both members survive, which is the lever a larger panel actually pulls.

Pairs are heavily dependent within a small haplotype set, so subsampling pairs
alone understates the benefit of a larger panel; the haplotype axis is the honest
one.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "python"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from gamma_smc_aou.run8_analysis import load_arm  # noqa: E402
from gamma_smc_aou.run9_config import (  # noqa: E402
    FOCAL_POSITION_BP,
    GENERATION_TIME_YEARS,
)

THRESHOLDS = (10_000.0, 30_000.0, 50_000.0)
PAIR_COUNTS = (1_000, 3_000, 10_000)
HAPLOTYPE_COUNTS = (50, 100, 200)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study-root", default=str(REPO / "sim_results_run9"))
    parser.add_argument("--seed", type=int, default=20290101)
    args = parser.parse_args(argv)

    rng = np.random.default_rng(args.seed)
    neutral = load_arm(args.study_root, None)
    if not neutral:
        raise RuntimeError("no neutral replicates found")
    print(f"neutral replicates: {len(neutral)}")

    rows = []
    for years in THRESHOLDS:
        generations = years / GENERATION_TIME_YEARS

        # Per replicate, the full-data statistic and the two subsample ladders.
        full, by_pairs, by_haps = [], {k: [] for k in PAIR_COUNTS}, {k: [] for k in HAPLOTYPE_COUNTS}
        for replicate in neutral:
            focal = int(np.argmin(np.abs(replicate.positions - FOCAL_POSITION_BP)))
            below = replicate.tmrca[focal] < generations
            full.append(float(below.mean()))

            for k in PAIR_COUNTS:
                if k > below.size:
                    continue
                pick = rng.choice(below.size, k, replace=False)
                by_pairs[k].append(float(below[pick].mean()))

            haplotypes = replicate.pair_haplotypes
            n_haps = int(haplotypes.max()) + 1
            for k in HAPLOTYPE_COUNTS:
                if k > n_haps:
                    continue
                keep = rng.choice(n_haps, k, replace=False)
                member = np.zeros(n_haps, dtype=bool)
                member[keep] = True
                mask = member[haplotypes[:, 0]] & member[haplotypes[:, 1]]
                if mask.sum() < 50:
                    continue
                by_haps[k].append(float(below[mask].mean()))

        full = np.asarray(full)
        # Binomial variance if the 10,000 pairs were independent, which they are
        # not -- included only as the floor sampling noise could contribute.
        p = full.mean()
        binomial = p * (1 - p) / 10_000.0
        rows.append(
            {
                "threshold_years": years,
                "axis": "full",
                "k": 10_000,
                "mean": p,
                "var_between_replicates": float(full.var(ddof=1)),
                "binomial_var_if_independent": binomial,
            }
        )
        for label, table in (("pairs", by_pairs), ("haplotypes", by_haps)):
            for k, values in table.items():
                if not values:
                    continue
                v = np.asarray(values)
                rows.append(
                    {
                        "threshold_years": years,
                        "axis": label,
                        "k": k,
                        "mean": float(v.mean()),
                        "var_between_replicates": float(v.var(ddof=1)),
                        "binomial_var_if_independent": np.nan,
                    }
                )

    table = pd.DataFrame(rows)
    out = Path(args.study_root) / "tables"
    out.mkdir(parents=True, exist_ok=True)
    table.to_csv(out / "variance_decomposition.tsv", sep="\t", index=False)

    pd.set_option("display.width", 220)
    for years in THRESHOLDS:
        block = table[table["threshold_years"] == years]
        full_var = float(block[block["axis"] == "full"]["var_between_replicates"].iloc[0])
        binom = float(block[block["axis"] == "full"]["binomial_var_if_independent"].iloc[0])
        print("=" * 76)
        print(f"P(TMRCA < {years:,.0f} y)   mean = {block['mean'].iloc[0]:.4f}")
        print("=" * 76)
        print(f"  variance across 300 neutral replicates : {full_var:.3e}  (sd {np.sqrt(full_var):.4f})")
        print(f"  binomial variance if pairs independent : {binom:.3e}")
        print(f"  ratio                                  : {full_var / binom:.1f}x")
        print()
        for axis in ("pairs", "haplotypes"):
            part = block[block["axis"] == axis].sort_values("k")
            if part.empty:
                continue
            print(f"  subsampling by {axis}:")
            for _, r in part.iterrows():
                ratio = r["var_between_replicates"] / full_var
                print(
                    f"    k = {int(r['k']):5d}  sd = {np.sqrt(r['var_between_replicates']):.4f}"
                    f"   variance x{ratio:.2f} of full"
                )
            print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
