#!/usr/bin/env python
"""The frequency-matched *neutral introgression* null for the hom-carrier statistic.

``run6_conditioned_null.py`` conditions on being homozygous for a single derived
SNP.  That is the wrong ascertainment, and wrong in a knowable direction: every
copy of a neutral SNP descends from one mutation, so two hom-derived haplotypes
are *forced* to coalesce after it.  An introgressed tract carries no such
guarantee -- the pulse delivered many archaic lineages, so two haplotypes that
are both introgressed may coalesce recently (same post-pulse founder) or deep in
the archaic branch (different founders).  The SNP null is therefore biased
*young* relative to the real ascertainment, which makes it conservative but not
correct.

This script builds the right null.  It simulates the CHB arm's own demography
(``OutOfAfricaArchaicAdmixture_5R19``) with **no selection at all**, calls archaic
ancestry per haplotype from the recorded migrations, and reports, for every
stride where enough homozygous-introgressed individuals exist, their
``P(TMRCA < x)`` alongside the local introgression frequency.

Pooling those strides genome-wide gives the null as a *function of introgression
frequency*, which is what the observation has to be matched against: a neutral
locus that has drifted to 60% introgression is the honest comparator for a
selected locus at 60%, and neither is comparable to a locus at 3%.
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
import stdpopsim  # noqa: E402

MODEL_ID = "OutOfAfricaArchaicAdmixture_5R19"
TARGET_POPULATION = "CHB"
ARCHAIC_POPULATION = "Neanderthal"
GENERATION_TIME = 29.0
SEQUENCE_LENGTH_BP = 10_000_000
RECOMBINATION_RATE = 1.0e-8
SAMPLE_DIPLOIDS = 100
STRIDE_BP = 50_000
THRESHOLDS_YEARS = (10_000.0, 20_000.0, 30_000.0, 50_000.0, 100_000.0)
MIN_HOM_INDIVIDUALS = 5


def archaic_population_ids(demography) -> set[int]:
    ids = set()
    for index, population in enumerate(demography.populations):
        if ARCHAIC_POPULATION.lower() in population.name.lower():
            ids.add(index)
    if not ids:
        raise RuntimeError(f"no population matching {ARCHAIC_POPULATION!r}")
    return ids


def one_replicate(seed: int) -> pd.DataFrame:
    species = stdpopsim.get_species("HomSap")
    model = species.get_demographic_model(MODEL_ID)
    demography = model.model
    archaic = archaic_population_ids(demography)

    ts = msprime.sim_ancestry(
        samples={TARGET_POPULATION: SAMPLE_DIPLOIDS},
        demography=demography,
        sequence_length=SEQUENCE_LENGTH_BP,
        recombination_rate=RECOMBINATION_RATE,
        record_migrations=True,
        random_seed=seed,
    )

    # A lineage is archaic over an interval if, going backwards, it migrated into
    # an archaic deme there. Collect those (node, interval) records once.
    into_archaic: dict[int, list[tuple[float, float]]] = {}
    migrations = ts.tables.migrations
    for node, left, right, dest in zip(
        migrations.node, migrations.left, migrations.right, migrations.dest
    ):
        if int(dest) in archaic:
            into_archaic.setdefault(int(node), []).append((float(left), float(right)))

    pairs = np.array([i.nodes for i in ts.individuals()], dtype=np.int64)
    samples = pairs.ravel()
    positions = np.arange(STRIDE_BP // 2, SEQUENCE_LENGTH_BP, STRIDE_BP, dtype=float)

    rows = []
    cursor = 0
    for tree in ts.trees():
        left, right = tree.interval
        while cursor < positions.size and positions[cursor] < right:
            position = positions[cursor]
            cursor += 1
            if position < left:
                continue

            introgressed = np.zeros(samples.size, dtype=bool)
            for index, sample in enumerate(samples):
                node = int(sample)
                while node != -1:
                    for lo, hi in into_archaic.get(node, ()):
                        if lo <= position < hi:
                            introgressed[index] = True
                            break
                    if introgressed[index]:
                        break
                    node = tree.parent(node)

            per_individual = introgressed.reshape(-1, 2)
            hom = per_individual.all(axis=1)
            frequency = float(introgressed.mean())
            if hom.sum() < MIN_HOM_INDIVIDUALS:
                continue
            tmrca = np.array(
                [tree.tmrca(int(a), int(b)) for a, b in pairs], dtype=float
            )
            for years in THRESHOLDS_YEARS:
                generations = years / GENERATION_TIME
                rows.append(
                    {
                        "seed": seed,
                        "position": position,
                        "introgression_frequency": frequency,
                        "n_hom_introgressed": int(hom.sum()),
                        "threshold_years": float(years),
                        "p_hom_introgressed": float(np.mean(tmrca[hom] < generations)),
                        "p_unconditional": float(np.mean(tmrca < generations)),
                    }
                )
        if cursor >= positions.size:
            break
    return pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=60)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument(
        "--out", default=str(REPO / "sim_results_run6" / "introgressed_null")
    )
    args = parser.parse_args(argv)

    seeds = [20260825 + i for i in range(args.replicates)]
    if args.workers == 1:
        frames = [one_replicate(s) for s in seeds]
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            frames = list(executor.map(one_replicate, seeds))
    table = pd.concat([f for f in frames if len(f)], ignore_index=True)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    table.to_csv(out / "introgressed_null_raw.tsv", sep="\t", index=False)

    edges = [0.0, 0.05, 0.10, 0.20, 0.35, 0.50, 0.70, 1.01]
    table["frequency_bin"] = pd.cut(table["introgression_frequency"], edges, right=False)
    summary = (
        table.groupby(["frequency_bin", "threshold_years"], observed=True)
        .agg(
            n_loci=("p_hom_introgressed", "size"),
            mean_n_hom=("n_hom_introgressed", "mean"),
            mean_uncond=("p_unconditional", "mean"),
            mean_hom=("p_hom_introgressed", "mean"),
            q95_hom=("p_hom_introgressed", lambda s: float(np.quantile(s, 0.95))),
            q99_hom=("p_hom_introgressed", lambda s: float(np.quantile(s, 0.99))),
        )
        .reset_index()
    )
    summary["inflation"] = summary["mean_hom"] - summary["mean_uncond"]
    summary.to_csv(out / "introgressed_null_summary.tsv", sep="\t", index=False)

    pd.set_option("display.width", 250)
    print(json.dumps({"replicates": args.replicates, "loci_rows": int(len(table))}, indent=2))
    print()
    print(summary.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
