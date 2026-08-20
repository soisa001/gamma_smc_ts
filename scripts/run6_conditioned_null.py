#!/usr/bin/env python
"""Measure the ascertainment bias in a hom-carrier TMRCA null.

The hom-carrier statistic is computed on pairs selected for carrying the *same
derived allele on both haplotypes*. That is not a random pair: two haplotypes
carrying the same derived allele must coalesce at or after that allele's origin,
so they share a more recent ancestor than two random haplotypes **even under
complete neutrality**. Comparing an ascertained observation against an
unconditional null therefore manufactures signal.

This script quantifies the bias directly. It simulates the EAS demography
neutrally, finds segregating sites near the focal base at a target derived
frequency, takes the individuals homozygous for the derived allele, and compares
their within-individual ``P(TMRCA < x)`` against the unconditional statistic at
the same position.

The size of the bias depends on frequency: a common derived allele is old, so
conditioning on it buys little, while a rare one is young and buys a lot. That
is why the null has to be matched on frequency, not merely conditioned.

msprime is used rather than SLiM because the null is neutral -- there is no sweep,
so the restriction that forced run5 to flatten the PHLASH history does not apply
and the full ten-thousand-epoch curve is carried exactly.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "python"))

import msprime  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from gamma_smc_aou.eas_sweep_models import (  # noqa: E402
    build_eas_demography_models,
    load_phlash_eas_npz,
)
from gamma_smc_aou.run6_config import (  # noqa: E402
    FOCAL_POSITION_BP,
    MUTATION_RATE,
    PHLASH_EAS_SHA256,
    PHLASH_EAS_RELATIVE_PATH,
    RECOMBINATION_RATE,
    SAMPLE_DIPLOIDS,
    SEQUENCE_LENGTH_BP,
)

DEMOGRAPHIES = {
    "eas": {"generation_time": 25.0, "population": "EAS"},
    "chb": {"generation_time": 29.0, "population": "CHB"},
}
THRESHOLDS_YEARS = (10_000.0, 20_000.0, 30_000.0, 50_000.0, 100_000.0)
FREQUENCY_BANDS = ((0.05, 0.15), (0.25, 0.35), (0.45, 0.55), (0.65, 0.75), (0.85, 0.95))
WINDOW_BP = 250_000
MIN_HOM_INDIVIDUALS = 10


def build_demography(repo_root: Path, kind: str) -> msprime.Demography:
    """EAS is the inferred history; CHB is 5R19, which controls for demography.

    The tract-based null can only be built under a model that *has* archaic
    admixture, so isolating the effect of the ascertainment from the effect of
    the demography needs this SNP-based null run under 5R19 as well.
    """
    if kind == "chb":
        import stdpopsim

        return stdpopsim.get_species("HomSap").get_demographic_model(
            "OutOfAfricaArchaicAdmixture_5R19"
        ).model
    artifact = load_phlash_eas_npz(
        repo_root / PHLASH_EAS_RELATIVE_PATH, expected_sha256=PHLASH_EAS_SHA256
    )
    return build_eas_demography_models(artifact)["median"].msprime_demography


def one_replicate(demography, seed: int, kind: str = "eas") -> pd.DataFrame:
    settings = DEMOGRAPHIES[kind]
    generation_time = settings["generation_time"]
    ts = msprime.sim_ancestry(
        samples=[
            msprime.SampleSet(SAMPLE_DIPLOIDS, population=settings["population"], ploidy=2)
        ],
        demography=demography,
        sequence_length=SEQUENCE_LENGTH_BP,
        recombination_rate=RECOMBINATION_RATE,
        random_seed=seed,
    )
    ts = msprime.sim_mutations(ts, rate=MUTATION_RATE, random_seed=seed + 1)

    pairs = np.array(
        [(i.nodes[0], i.nodes[1]) for i in ts.individuals()], dtype=np.int64
    )
    ordered = pairs.ravel()
    low = FOCAL_POSITION_BP - WINDOW_BP
    high = FOCAL_POSITION_BP + WINDOW_BP

    # One candidate site per frequency band, nearest the focal base.
    wanted = {band: None for band in FREQUENCY_BANDS}
    for variant in ts.variants(samples=ordered, left=low, right=high):
        genotypes = np.asarray(variant.genotypes)
        if genotypes.max() < 1:
            continue
        derived = (genotypes != 0).astype(np.int8)
        frequency = derived.mean()
        for band in FREQUENCY_BANDS:
            if not (band[0] <= frequency < band[1]):
                continue
            distance = abs(variant.site.position - FOCAL_POSITION_BP)
            current = wanted[band]
            if current is None or distance < current[0]:
                wanted[band] = (distance, variant.site.position, derived.reshape(-1, 2))

    rows = []
    for band, hit in wanted.items():
        if hit is None:
            continue
        _, position, per_individual = hit
        counts = per_individual.sum(axis=1)
        tree = ts.at(position)
        tmrca = np.array(
            [tree.tmrca(int(a), int(b)) for a, b in pairs], dtype=float
        )
        hom = counts == 2
        if hom.sum() < MIN_HOM_INDIVIDUALS:
            continue
        for years in THRESHOLDS_YEARS:
            generations = years / generation_time
            rows.append(
                {
                    "seed": seed,
                    "band_low": band[0],
                    "band_high": band[1],
                    "position": float(position),
                    "derived_frequency": float(counts.sum() / (2 * len(pairs))),
                    "n_hom": int(hom.sum()),
                    "threshold_years": float(years),
                    "p_hom_derived": float(np.mean(tmrca[hom] < generations)),
                    "p_unconditional": float(np.mean(tmrca < generations)),
                }
            )
    return pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=200)
    parser.add_argument("--demography", choices=sorted(DEMOGRAPHIES), default="eas")
    parser.add_argument("--repo-root", default=str(REPO))
    parser.add_argument(
        "--out", default=str(REPO / "sim_results_run6" / "conditioned_null")
    )
    args = parser.parse_args(argv)

    repo_root = Path(args.repo_root)
    demography = build_demography(repo_root, args.demography)
    frames = [
        one_replicate(demography, seed=20260824 + 2 * i, kind=args.demography)
        for i in range(args.replicates)
    ]
    table = pd.concat([f for f in frames if len(f)], ignore_index=True)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    suffix = "" if args.demography == "eas" else f"_{args.demography}"
    table.to_csv(out / f"conditioned_null{suffix}_raw.tsv", sep="\t", index=False)

    summary = (
        table.groupby(["band_low", "band_high", "threshold_years"])
        .agg(
            n=("p_hom_derived", "size"),
            mean_hom=("p_hom_derived", "mean"),
            mean_unconditional=("p_unconditional", "mean"),
            q95_hom=("p_hom_derived", lambda s: float(np.quantile(s, 0.95))),
            q95_unconditional=("p_unconditional", lambda s: float(np.quantile(s, 0.95))),
            mean_n_hom=("n_hom", "mean"),
        )
        .reset_index()
    )
    summary["inflation"] = summary["mean_hom"] - summary["mean_unconditional"]
    summary["ratio"] = summary["mean_hom"] / summary["mean_unconditional"].replace(0, np.nan)
    summary.to_csv(out / f"conditioned_null{suffix}_summary.tsv", sep="\t", index=False)

    pd.set_option("display.width", 220)
    print(json.dumps({"replicates": args.replicates, "rows": int(len(table))}, indent=2))
    print()
    print(
        summary[
            ["band_low", "threshold_years", "n", "mean_n_hom",
             "mean_unconditional", "mean_hom", "inflation", "ratio"]
        ].to_string(index=False, float_format=lambda x: f"{x:.3f}")
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
