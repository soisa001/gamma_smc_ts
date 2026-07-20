from __future__ import annotations

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
    with destination.open("w", encoding="utf-8") as output:
        if individuals:
            ts.write_vcf(output, individuals=individuals)
        else:
            ts.write_vcf(output, ploidy=2)
    return destination
