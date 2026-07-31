#!/usr/bin/env python3
"""Derive a selected-study source at a new mutation rate on fixed ancestry."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path
from time import perf_counter

import msprime
import tskit

from gamma_smc_aou.tree_sequence import tree_sequence_to_vcf


TREE_NAME = "selected_s0p05_af30.trees"
VCF_NAME = "selected_s0p05_af30.vcf.gz"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _mutation_free(ts: tskit.TreeSequence) -> tskit.TreeSequence:
    tables = ts.dump_tables()
    tables.mutations.clear()
    tables.sites.clear()
    return tables.tree_sequence()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--template-source", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--mutation-rate", type=float, required=True)
    parser.add_argument("--mutation-seed", type=int, required=True)
    args = parser.parse_args()

    if args.mutation_rate <= 0 or args.mutation_seed <= 0:
        parser.error("mutation rate and seed must be positive")

    source = args.template_source.resolve()
    output = args.output_dir.resolve()
    if output.exists() and any(output.iterdir()):
        parser.error("--output-dir must be empty to avoid mixing study versions")
    output.mkdir(parents=True, exist_ok=True)
    started = perf_counter()

    excluded = {
        TREE_NAME,
        VCF_NAME,
        "metrics.json",
        "mutation_overlay_metadata.json",
    }
    for path in source.iterdir():
        if path.is_file() and path.name not in excluded:
            shutil.copy2(path, output / path.name)

    template = tskit.load(source / TREE_NAME)
    ancestry = _mutation_free(template)
    selected = msprime.sim_mutations(
        ancestry,
        rate=args.mutation_rate,
        model=msprime.BinaryMutationModel(),
        random_seed=args.mutation_seed,
    )
    selected.dump(output / TREE_NAME)
    tree_sequence_to_vcf(output / TREE_NAME, output / VCF_NAME)

    with (source / "metrics.json").open(encoding="utf-8") as handle:
        metrics = json.load(handle)
    metrics["mutation_rate"] = float(args.mutation_rate)
    metrics["mutation_overlay_derivation"] = {
        "genealogy_source": str(
            Path("..") / source.name / TREE_NAME
        ).replace("\\", "/"),
        "sites_and_mutations_cleared": True,
        "mutation_free_tables_equal_ignoring_provenance": bool(
            ancestry.tables.equals(
                _mutation_free(selected).tables,
                ignore_provenance=True,
            )
        ),
        "model": "msprime.BinaryMutationModel",
        "random_seed": int(args.mutation_seed),
        "retained_sites": int(selected.num_sites),
    }
    metrics["matched_overlay_elapsed_seconds"] = float(perf_counter() - started)
    with (output / "metrics.json").open("w", encoding="utf-8") as handle:
        json.dump(metrics, handle, indent=2)
        handle.write("\n")

    metadata = {
        **metrics["mutation_overlay_derivation"],
        "mutation_rate": float(args.mutation_rate),
        "tree_sha256": _sha256(output / TREE_NAME),
        "vcf_sha256": _sha256(output / VCF_NAME),
    }
    with (output / "mutation_overlay_metadata.json").open(
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(metadata, handle, indent=2)
        handle.write("\n")

    print(json.dumps({
        "output_dir": str(output),
        "mutation_rate": args.mutation_rate,
        "mutation_seed": args.mutation_seed,
        "retained_sites": selected.num_sites,
        "retained_mutations": selected.num_mutations,
        "elapsed_seconds": perf_counter() - started,
    }, indent=2))


if __name__ == "__main__":
    main()
