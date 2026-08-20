#!/usr/bin/env python
"""Check that the grafted EAS model produces realistic introgression.

A grafted model is only useful for power simulation if the haplotypes it
delivers look like real archaic tracts. Four things are checked against
independent expectations rather than against each other:

1. **Ancestry fraction** should land on the admixture proportion (~2.5%).
2. **Tract length** should be near ``1 / (r t)`` for the pulse time -- about
   44 kb at 2,272 generations -- because that is how far recombination has had
   time to erode an intact segment.
3. **Archaic vs modern depth** should exceed the split time, since those two
   haplotypes cannot meet more recently than the branch they sit on.
4. **The dissociation** -- tract length divided by the length a panmictic
   population would give at the same coalescent depth -- should be far above 1.
   This is the quantity no standing variant can reproduce, and it is the reason
   the graft exists.

Check 4 is the load-bearing one. ``modern_vs_modern`` pairs supply the calibration:
being ordinary within-population pairs, their ratio must come out at 1, and it
does, so a ratio far above 1 for archaic tracts is a real dissociation rather than
an artefact of how the span is measured.
"""

from __future__ import annotations

import argparse
import json
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

SEQUENCE_LENGTH_BP = 10_000_000
RECOMBINATION_RATE = 1.0e-8
SAMPLE_DIPLOIDS = 100
PROBE_POSITIONS = 25
MAX_PAIRS_PER_CLASS = 60


def shared_span(ts, a: int, b: int, position: float) -> float:
    tree = ts.at(position)
    reference = tree.tmrca(a, b)
    left, right = tree.interval.left, tree.interval.right

    probe = tree.copy()
    while probe.index > 0:
        probe.prev()
        if abs(probe.tmrca(a, b) - reference) > 1e-9:
            break
        left = probe.interval.left

    probe = tree.copy()
    while probe.index < ts.num_trees - 1:
        probe.next()
        if abs(probe.tmrca(a, b) - reference) > 1e-9:
            break
        right = probe.interval.right
    return right - left


TRACT_WINDOW_BP = 2_000_000
TRACT_STEP_BP = 5_000


def archaic_tract_lengths(ts, model, samples) -> list[float]:
    """Contiguous archaic run lengths, which is what an hmmix call reports.

    This is *not* the pairwise IBD span: two archaic haplotypes from different
    founding lineages share no intact segment, yet each still carries a long
    archaic tract of its own. Tract length is the quantity set by the pulse time,
    and it is the one that has to look realistic.
    """
    centre = ts.sequence_length / 2
    left = centre - TRACT_WINDOW_BP / 2
    grid = np.arange(left, left + TRACT_WINDOW_BP, TRACT_STEP_BP)
    calls = np.vstack(
        [call_archaic_ancestry(ts, model, float(p), samples) for p in grid]
    )  # (n_positions, n_samples)

    lengths = []
    for column in range(calls.shape[1]):
        series = calls[:, column]
        if not series.any():
            continue
        padded = np.concatenate(([False], series, [False]))
        edges = np.flatnonzero(padded[1:] != padded[:-1])
        for start, stop in zip(edges[0::2], edges[1::2]):
            # Runs touching a boundary are censored, so they would bias the mean.
            if start == 0 or stop == series.size:
                continue
            lengths.append(float((stop - start) * TRACT_STEP_BP))
    return lengths


def one_replicate(job) -> pd.DataFrame:
    seed, repo_root = job
    model = build_eas_with_archaic(repo_root)
    ts = msprime.sim_ancestry(
        samples={"EAS": SAMPLE_DIPLOIDS},
        demography=model.demography,
        sequence_length=SEQUENCE_LENGTH_BP,
        recombination_rate=RECOMBINATION_RATE,
        random_seed=seed,
    )
    rng = np.random.default_rng(seed)
    samples = np.asarray(ts.samples(), dtype=np.int64)
    positions = np.linspace(
        0.05 * SEQUENCE_LENGTH_BP, 0.95 * SEQUENCE_LENGTH_BP, PROBE_POSITIONS
    )

    rows = []
    for position in positions:
        flags = call_archaic_ancestry(ts, model, float(position), samples)
        archaic = samples[flags]
        modern = samples[~flags]
        tree = ts.at(float(position))

        def sample_pairs(left_pool, right_pool, label):
            if left_pool.size < 1 or right_pool.size < 1:
                return
            if label == "archaic_vs_archaic":
                if left_pool.size < 2:
                    return
                candidates = [
                    (int(left_pool[i]), int(left_pool[j]))
                    for i in range(left_pool.size)
                    for j in range(i + 1, left_pool.size)
                ]
            else:
                candidates = [
                    (int(a), int(b))
                    for a in left_pool[: min(left_pool.size, 30)]
                    for b in right_pool[: min(right_pool.size, 30)]
                ]
            if len(candidates) > MAX_PAIRS_PER_CLASS:
                picked = rng.choice(
                    len(candidates), MAX_PAIRS_PER_CLASS, replace=False
                )
                candidates = [candidates[k] for k in picked]
            for a, b in candidates:
                depth = float(tree.tmrca(a, b))
                rows.append(
                    {
                        "seed": seed,
                        "position": float(position),
                        "archaic_fraction": float(flags.mean()),
                        "pair_class": label,
                        "tmrca_generations": depth,
                        "shared_span_bp": shared_span(ts, a, b, float(position)),
                    }
                )

        sample_pairs(archaic, archaic, "archaic_vs_archaic")
        sample_pairs(archaic, modern, "archaic_vs_modern")
        sample_pairs(modern, modern, "modern_vs_modern")

    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame, []
    # Two archaic haplotypes are either recent copies of one introgressing
    # lineage or descendants of different ones. Only the second kind carries the
    # dissociation, so pooling them hides exactly what this is measuring.
    deep = (frame["pair_class"] == "archaic_vs_archaic") & (
        frame["tmrca_generations"] >= model.pulse_generations
    )
    recent = (frame["pair_class"] == "archaic_vs_archaic") & ~deep
    frame.loc[deep, "pair_class"] = "archaic_pair_distinct_founders"
    frame.loc[recent, "pair_class"] = "archaic_pair_same_founder"
    return frame, archaic_tract_lengths(ts, model, samples)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=30)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--repo-root", default=str(REPO))
    parser.add_argument("--out", default=str(REPO / "sim_results_run7" / "graft_validation"))
    args = parser.parse_args(argv)

    model = build_eas_with_archaic(args.repo_root)
    jobs = [(20261001 + i, args.repo_root) for i in range(args.replicates)]
    if args.workers == 1:
        results = [one_replicate(job) for job in jobs]
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            results = list(executor.map(one_replicate, jobs))
    frames = [f for f, _ in results if len(f)]
    tract_lengths = np.concatenate([np.asarray(t) for _, t in results if len(t)])
    table = pd.concat(frames, ignore_index=True)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    table.to_csv(out / "graft_validation_raw.tsv", sep="\t", index=False)
    (out / "grafted_model.json").write_text(
        json.dumps(model.record, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    summary = (
        table.groupby("pair_class")
        .agg(
            n_pairs=("tmrca_generations", "size"),
            median_tmrca_gen=("tmrca_generations", "median"),
            median_tmrca_ky=(
                "tmrca_generations",
                lambda s: float(np.median(s)) * model.generation_time / 1000,
            ),
            q25_tmrca_gen=("tmrca_generations", lambda s: float(np.quantile(s, 0.25))),
            q75_tmrca_gen=("tmrca_generations", lambda s: float(np.quantile(s, 0.75))),
            median_span_kb=("shared_span_bp", lambda s: float(np.median(s)) / 1000),
        )
        .reset_index()
    )
    summary["panmictic_span_kb"] = (
        1.0 / (RECOMBINATION_RATE * summary["median_tmrca_gen"]) / 1000
    )
    summary["dissociation"] = summary["median_span_kb"] / summary["panmictic_span_kb"]
    summary.to_csv(out / "graft_validation_summary.tsv", sep="\t", index=False)

    fraction = table.groupby(["seed", "position"])["archaic_fraction"].first()
    modern_depth = float(
        table.loc[table["pair_class"] == "archaic_vs_modern", "tmrca_generations"].median()
    )
    checks = {
        "archaic_fraction_mean": float(fraction.mean()),
        "archaic_fraction_target": model.admixture_proportion,
        "expected_tract_length_kb": model.expected_tract_length_bp(RECOMBINATION_RATE)
        / 1000,
        "observed_tract_length_kb_mean": float(np.mean(tract_lengths)) / 1000,
        "observed_tract_length_kb_median": float(np.median(tract_lengths)) / 1000,
        "n_tracts": int(tract_lengths.size),
        "archaic_vs_modern_depth_generations": modern_depth,
        "panmictic_span_at_that_depth_kb": 1.0
        / (RECOMBINATION_RATE * modern_depth)
        / 1000,
        "dissociation_tract_vs_panmictic": (
            float(np.mean(tract_lengths))
            / (1.0 / (RECOMBINATION_RATE * modern_depth))
        ),
        "split_generations": model.split_generations,
        "pulse_generations": model.pulse_generations,
    }

    pd.set_option("display.width", 250)
    print(json.dumps(checks, indent=2))
    print()
    print(summary.to_string(index=False, float_format=lambda x: f"{x:.2f}"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
