"""The EAS history with an archaic branch grafted on, for power simulation only.

The study uses two models with two different jobs, and it matters that they stay
separate:

* **The base EAS history** calibrates p-values on empirical windows. It was
  inferred on the target samples themselves, so it is the null of record. It
  contains no archaic admixture, and it does not need any.
* **This grafted model** exists only to ask whether the statistic *can* detect
  selection on an introgressed haplotype -- power, and the assumptions behind it.
  Nothing computed from it is ever used as a null.

The graft is deliberately minimal. The inferred EAS trajectory is left exactly as
it is and carries the whole human lineage; a Neanderthal branch is attached to it
going backwards, and a single pulse delivers archaic haplotypes into EAS. Making
Neanderthal a *derived* population that merges into EAS at the split -- rather
than splitting EAS into an ancestral population -- is what keeps the PHLASH
epochs valid all the way back, and it mirrors how ``OutOfAfricaArchaicAdmixture_5R19``
attaches its archaic branches.

Why a graft is needed at all: in a single panmictic population an allele of age
``T`` carries a shared haplotype of order ``1 / (r T)``, so depth and tract
length are locked together and older necessarily means shorter. Introgression
breaks that link -- the tract entered at the pulse, so it is long, while the two
archaic haplotypes it sits on meet in the archaic branch, so they are deep. No
standing variant in a one-population model can be both, which is why an
artificial "introgression from standing variation" arm cannot stand in for the
real thing: conditioning it to be young enough to have long tracts also makes its
carriers recently coalescing, which manufactures the very signal under test.

Parameter provenance
--------------------
``HUMAN_NEANDERTHAL_SPLIT_GENERATIONS`` (27,840 gen = 696 ky at 25 y/gen) and
``INTROGRESSION_PULSE_GENERATIONS`` (2,272 gen = 56.8 ky) are the values this
repository already used for the ``AncientEurasia_9K19`` work, and both sit inside
the published ranges -- a 550-765 ky split (Prufer et al. 2014) and a 47-65 ky
admixture date (Sankararaman et al. 2012; Fu et al. 2014, Ust'-Ishim).  The
archaic effective size of 3,600 and the option of continuous rather than pulsed
admixture are taken from ``OutOfAfricaArchaicAdmixture_5R19`` (Ragsdale & Gravel
2019).  The default admixture proportion of 2.5% is the East Asian Neanderthal
ancestry estimate (Prufer et al. 2014; Vernot & Akey 2015); the CHB arms in this
repository used 2.96%, which is within the same range.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import msprime
import numpy as np

from .eas_sweep_models import (
    ARCHAIC_POPULATION,
    EAS_POPULATION,
    GENERATION_TIME_YEARS,
    HUMAN_NEANDERTHAL_SPLIT_GENERATIONS,
    INTROGRESSION_PULSE_GENERATIONS,
    build_eas_demography_models,
    load_phlash_eas_npz,
)
from .run6_config import PHLASH_EAS_RELATIVE_PATH, PHLASH_EAS_SHA256

__all__ = [
    "ARCHAIC_EFFECTIVE_SIZE",
    "ADMIXTURE_PROPORTION",
    "CENSUS_OFFSET_GENERATIONS",
    "GraftedModel",
    "build_eas_with_archaic",
    "call_archaic_ancestry",
]

#: Ragsdale & Gravel (2019), as implemented in OutOfAfricaArchaicAdmixture_5R19.
ARCHAIC_EFFECTIVE_SIZE = 3_600.0
#: East Asian Neanderthal ancestry (Prufer et al. 2014; Vernot & Akey 2015).
ADMIXTURE_PROPORTION = 0.025
#: The census sits just *older* than the pulse, so a lineage sampled there is
#: already in whichever population it came from.
CENSUS_OFFSET_GENERATIONS = 1.0


@dataclass(frozen=True)
class GraftedModel:
    demography: msprime.Demography
    split_generations: float
    pulse_generations: float
    census_generations: float
    archaic_effective_size: float
    admixture_proportion: float
    generation_time: float
    record: dict[str, Any] = field(default_factory=dict)

    def expected_tract_length_bp(self, recombination_rate: float) -> float:
        """Mean intact span of a pulse tract, ``1 / (r t)``."""
        return 1.0 / (recombination_rate * self.pulse_generations)


def build_eas_with_archaic(
    repo_root: str | Path,
    *,
    quantile: str = "median",
    split_generations: float = HUMAN_NEANDERTHAL_SPLIT_GENERATIONS,
    pulse_generations: float = INTROGRESSION_PULSE_GENERATIONS,
    archaic_effective_size: float = ARCHAIC_EFFECTIVE_SIZE,
    admixture_proportion: float = ADMIXTURE_PROPORTION,
    add_census: bool = True,
) -> GraftedModel:
    """The inferred EAS history plus a Neanderthal branch and one pulse."""
    if not 0.0 < admixture_proportion < 1.0:
        raise ValueError("admixture_proportion must lie in (0, 1)")
    if not 0.0 < pulse_generations < split_generations:
        raise ValueError("the pulse must fall between the present and the split")

    artifact = load_phlash_eas_npz(
        Path(repo_root) / PHLASH_EAS_RELATIVE_PATH, expected_sha256=PHLASH_EAS_SHA256
    )
    base = build_eas_demography_models(artifact)[quantile]

    # Rebuilt rather than mutated: the shared model object is reused elsewhere.
    demography = msprime.Demography()
    demography.add_population(
        name=EAS_POPULATION, initial_size=float(base.ne[0])
    )
    demography.add_population(
        name=ARCHAIC_POPULATION, initial_size=float(archaic_effective_size)
    )
    for time, size in zip(base.time_generations[1:], base.ne[1:]):
        # The PHLASH epochs describe the human lineage, which here runs the whole
        # way back; the archaic branch keeps its own constant size until it joins.
        demography.add_population_parameters_change(
            time=float(time), initial_size=float(size), population=EAS_POPULATION
        )

    census_generations = float(pulse_generations) + CENSUS_OFFSET_GENERATIONS
    # Backwards in time a fraction of EAS lineages step into the archaic branch.
    demography.add_mass_migration(
        time=float(pulse_generations),
        source=EAS_POPULATION,
        dest=ARCHAIC_POPULATION,
        proportion=float(admixture_proportion),
    )
    if add_census:
        demography.add_census(time=census_generations)
    # A mass migration rather than ``add_population_split``: a split marks the
    # ancestral population inactive until its split time, which would leave EAS
    # unable to coalesce for the whole 27,840 generations before it and push
    # every pair past the split. 5R19 joins its archaic branches the same way.
    demography.add_mass_migration(
        time=float(split_generations),
        source=ARCHAIC_POPULATION,
        dest=EAS_POPULATION,
        proportion=1.0,
    )
    demography.sort_events()
    demography.debug()

    record = {
        "schema": "gamma-smc.run7-grafted-model/v1",
        "purpose": "power simulation only; never used as a null",
        "base_model": base.label,
        "generation_time_years": GENERATION_TIME_YEARS,
        "split_generations": float(split_generations),
        "split_years": float(split_generations) * GENERATION_TIME_YEARS,
        "pulse_generations": float(pulse_generations),
        "pulse_years": float(pulse_generations) * GENERATION_TIME_YEARS,
        "census_generations": census_generations,
        "archaic_effective_size": float(archaic_effective_size),
        "admixture_proportion": float(admixture_proportion),
        "provenance": {
            "split": "Prufer et al. 2014 (550-765 ky)",
            "pulse": "Sankararaman et al. 2012; Fu et al. 2014 (47-65 ky)",
            "archaic_ne": "Ragsdale & Gravel 2019 (OutOfAfricaArchaicAdmixture_5R19)",
            "proportion": "Prufer et al. 2014; Vernot & Akey 2015 (~2-3% in East Asia)",
        },
    }
    return GraftedModel(
        demography=demography,
        split_generations=float(split_generations),
        pulse_generations=float(pulse_generations),
        census_generations=census_generations,
        archaic_effective_size=float(archaic_effective_size),
        admixture_proportion=float(admixture_proportion),
        generation_time=GENERATION_TIME_YEARS,
        record=record,
    )


def call_archaic_ancestry(
    ts, model: GraftedModel, position: float, samples: np.ndarray | None = None
) -> np.ndarray:
    """Per-haplotype archaic ancestry at ``position``, from the census nodes.

    Migration records cannot be used for this: ``sim_ancestry`` simplifies the
    output, and the nodes those records name are frequently gone from the tree.
    A census inserts a node on *every* lineage at its time and those nodes are
    retained, so the ancestor of a sample at the census is always reachable by
    walking the local tree, and its population is the answer.
    """
    if samples is None:
        samples = np.asarray(ts.samples(), dtype=np.int64)

    archaic_id = [
        index
        for index, population in enumerate(ts.tables.populations)
        if (population.metadata or {}).get("name") == ARCHAIC_POPULATION
    ]
    if not archaic_id:
        raise RuntimeError(f"no {ARCHAIC_POPULATION!r} population in the tree sequence")
    archaic = archaic_id[0]

    node_time = ts.tables.nodes.time
    node_population = ts.tables.nodes.population
    tree = ts.at(position)

    flags = np.zeros(samples.size, dtype=bool)
    for index, sample in enumerate(samples):
        node = int(sample)
        while node != -1 and node_time[node] < model.census_generations:
            node = tree.parent(node)
        if node != -1:
            flags[index] = node_population[node] == archaic
    return flags
