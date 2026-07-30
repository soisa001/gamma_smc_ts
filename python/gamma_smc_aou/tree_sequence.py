from __future__ import annotations

import gzip
from pathlib import Path

import tskit


def load_tree_sequence(path: str | Path):
    path = Path(path)
    if path.suffix.lower() == ".tsz":
        try:
            import tszip
        except ImportError as error:
            raise RuntimeError("Reading .tsz requires: pip install tszip") from error
        return tszip.load(path)
    return tskit.load(path)


def diploid_individuals(ts) -> list[int]:
    sample_set = set(ts.samples())
    individuals = []
    assigned = set()
    for individual in ts.individuals():
        nodes = [node for node in individual.nodes if node in sample_set]
        if nodes:
            if len(nodes) != 2:
                raise ValueError(f"individual {individual.id} has {len(nodes)} sample nodes; expected 2")
            individuals.append(individual.id)
            assigned.update(nodes)
    if assigned and assigned != sample_set:
        raise ValueError("some sample nodes are not assigned to diploid individuals")
    if not assigned and ts.num_samples % 2:
        raise ValueError("unassigned samples cannot be grouped into diploids: odd sample count")
    return individuals


def tree_sequence_to_vcf(source: str | Path, destination: str | Path) -> Path:
    ts = load_tree_sequence(source)
    individuals = diploid_individuals(ts)
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.suffix.lower() == ".gz":
        output_handle = gzip.open(destination, "wt", encoding="utf-8")
    else:
        output_handle = destination.open("w", encoding="utf-8")
    # Tskit coordinates are 0-based, while VCF POS must be at least 1. The
    # legacy transform preserves the usual rounded positions, maps coordinate
    # zero to one, and disambiguates any resulting duplicate integer positions.
    write_options = {"position_transform": "legacy"}
    with output_handle as output:
        if individuals:
            ts.write_vcf(output, individuals=individuals, **write_options)
        else:
            ts.write_vcf(output, ploidy=2, **write_options)
    return destination
