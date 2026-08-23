#!/usr/bin/env python
"""The carrier-matched null: is conditioning detecting selection, or introgression?

Every carrier-carrier result so far is scored against a neutral arm that contains
no carriers at all, because the neutral simulations carry no allele. That
comparison cannot separate two very different claims:

    (a) this locus was under selection, and
    (b) this locus is introgressed.

Two haplotypes that both carry an archaic tract are already unusual -- they must
either share a post-pulse ancestor or meet back in the archaic branch -- so their
``P(TMRCA < x)`` differs from a random pair's under strict neutrality. run6
measured that inflation at 1.4-4.8x depending on frequency. Against a null with no
carriers, (b) alone would look like a detection.

This builds the honest comparator. Neutral replicates of the *same* grafted model
are simulated with a census at the pulse, archaic ancestry is called per
haplotype, and every stride whose archaic frequency drifted into the range real
introgressed regions reach is scored exactly as the selected arms are. Binning on
that frequency gives a null matched on the thing being conditioned upon.

msprime rather than SLiM: the null is neutral, so there is no sweep to forward
simulate, and a census event gives exact ancestry rather than an hmmix-style
reconstruction.

Cost control: ancestry calling is cheap relative to pairwise TMRCA, so strides
are filtered on frequency first and only survivors pay for the TMRCA pass.
"""

from __future__ import annotations

import argparse
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "python"))

import msprime  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from gamma_smc_aou.run7_models import (  # noqa: E402
    build_eas_with_archaic,
    call_archaic_ancestry,
)
from gamma_smc_aou.run7_config import (  # noqa: E402
    GENERATION_TIME_YEARS,
    RECOMBINATION_RATE,
    SEQUENCE_LENGTH_BP,
)
from gamma_smc_aou.run9_config import SAMPLE_DIPLOIDS  # noqa: E402

THRESHOLDS = (10_000.0, 30_000.0, 50_000.0)
STRIDE_BP = 50_000
#: Only strides in this band are worth the TMRCA pass; below it there are too few
#: carrier pairs to score, above it drift essentially never reaches.
FREQUENCY_BAND = (0.10, 0.70)
MIN_CARRIER_PAIRS = 50
N_PAIRS = 10_000
PAIRS_SEED = 1729


def one_replicate(job) -> pd.DataFrame:
    seed, repo_root = job
    model = build_eas_with_archaic(repo_root, add_census=True)
    ts = msprime.sim_ancestry(
        samples={"EAS": SAMPLE_DIPLOIDS},
        demography=model.demography,
        sequence_length=SEQUENCE_LENGTH_BP,
        recombination_rate=RECOMBINATION_RATE,
        random_seed=seed,
    )
    samples = np.asarray(ts.samples(), dtype=np.int64)
    n_haps = samples.size

    rng = np.random.default_rng(PAIRS_SEED)
    all_pairs = np.array(
        [(i, j) for i in range(n_haps) for j in range(i + 1, n_haps)], dtype=np.int64
    )
    pick = rng.choice(all_pairs.shape[0], min(N_PAIRS, all_pairs.shape[0]), replace=False)
    pair_index = all_pairs[pick]
    node_pairs = samples[pair_index]

    positions = np.arange(STRIDE_BP // 2, SEQUENCE_LENGTH_BP, STRIDE_BP, dtype=float)
    rows = []
    for position in positions:
        carriers = call_archaic_ancestry(ts, model, float(position), samples)
        frequency = float(carriers.mean())
        if not (FREQUENCY_BAND[0] <= frequency <= FREQUENCY_BAND[1]):
            continue
        left = carriers[pair_index[:, 0]]
        right = carriers[pair_index[:, 1]]
        cc = left & right
        if int(cc.sum()) < MIN_CARRIER_PAIRS:
            continue

        tree = ts.at(float(position))
        tmrca = np.array(
            [tree.tmrca(int(a), int(b)) for a, b in node_pairs], dtype=float
        )
        for years in THRESHOLDS:
            generations = years / GENERATION_TIME_YEARS
            below = tmrca < generations
            rows.append(
                {
                    "seed": seed,
                    "position": float(position),
                    "archaic_frequency": frequency,
                    "n_carrier_pairs": int(cc.sum()),
                    "threshold_years": years,
                    "carrier_carrier": float(below[cc].mean()),
                    "overall": float(below.mean()),
                }
            )
    return pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=600)
    parser.add_argument("--workers", type=int, default=20)
    parser.add_argument("--repo-root", default=str(REPO))
    parser.add_argument("--out", default=str(REPO / "sim_results_run10"))
    args = parser.parse_args(argv)

    jobs = [(20300101 + i, args.repo_root) for i in range(args.replicates)]
    if args.workers == 1:
        frames = [one_replicate(job) for job in jobs]
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            frames = list(executor.map(one_replicate, jobs))
    table = pd.concat([f for f in frames if len(f)], ignore_index=True)

    out = Path(args.out) / "tables"
    out.mkdir(parents=True, exist_ok=True)
    table.to_csv(out / "matched_null_raw.tsv", sep="\t", index=False)

    edges = [0.10, 0.15, 0.25, 0.35, 0.45, 0.60, 0.71]
    table["af_bin"] = pd.cut(table["archaic_frequency"], edges, right=False)
    summary = (
        table.groupby(["af_bin", "threshold_years"], observed=True)
        .agg(
            n_loci=("carrier_carrier", "size"),
            mean_frequency=("archaic_frequency", "mean"),
            mean_carrier_pairs=("n_carrier_pairs", "mean"),
            overall_mean=("overall", "mean"),
            cc_mean=("carrier_carrier", "mean"),
            cc_q95=("carrier_carrier", lambda s: float(np.quantile(s, 0.95))),
            cc_q99=("carrier_carrier", lambda s: float(np.quantile(s, 0.99))),
        )
        .reset_index()
    )
    summary["inflation_vs_overall"] = summary["cc_mean"] / summary["overall_mean"]
    summary.to_csv(out / "matched_null_summary.tsv", sep="\t", index=False)

    pd.set_option("display.width", 250)
    print(f"replicates: {args.replicates}   scored loci: {table['position'].size // len(THRESHOLDS)}")
    print()
    print(summary.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    print()
    print("  cc_mean is the carrier-carrier statistic under strict NEUTRALITY at")
    print("  the matching archaic frequency. Selected replicates must beat this,")
    print("  not the carrier-free neutral arm used so far.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
