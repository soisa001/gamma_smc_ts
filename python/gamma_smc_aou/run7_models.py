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
Defaults come from :mod:`run7_config` and are the rounded literature values: a
700 kya split (Prufer et al. 2014 give 550-765 ky), a 55 kya pulse (Sankararaman
et al. 2012 and Fu et al. 2014 give 47-65 ky), an archaic effective size of 3,600
(Ragsdale & Gravel 2019), and 2.5% admixture (Prufer et al. 2014; Vernot & Akey
2015).  At 25 years per generation that is 28,000 and 2,200 generations.
"""

from __future__ import annotations

import hashlib
import re
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterator

import msprime
import numpy as np
import stdpopsim
from stdpopsim import slim_engine

from .eas_sweep_models import (
    build_eas_demography_models,
    load_phlash_eas_npz,
)
from .run7_config import (
    ADMIXTURE_PROPORTION as _ADMIXTURE_PROPORTION,
    ARCHAIC_EFFECTIVE_SIZE as _ARCHAIC_EFFECTIVE_SIZE,
    ARCHAIC_POPULATION,
    DOMINANCE_COEFFICIENT,
    EAS_POPULATION,
    FOCAL_POSITION_BP,
    FOCAL_SITE_ID,
    GENERATION_TIME_YEARS,
    MUTATION_RATE,
    PHLASH_EAS_RELATIVE_PATH,
    PHLASH_EAS_SHA256,
    PULSE_GENERATIONS,
    RECOMBINATION_RATE,
    REQUIRED_STDPOPSIM_VERSION,
    CONDITION_ON_SURVIVAL,
    SEQUENCE_LENGTH_BP,
    SLIM_SCALING_FACTOR,
    SPLIT_GENERATIONS,
    TARGET_FREQUENCY_BAND,
)

__all__ = [
    "ARCHAIC_EFFECTIVE_SIZE",
    "ADMIXTURE_PROPORTION",
    "CENSUS_OFFSET_GENERATIONS",
    "GraftedModel",
    "build_contig",
    "build_eas_with_archaic",
    "build_denovo_events",
    "build_extended_events",
    "build_stdpopsim_model",
    "call_archaic_ancestry",
    "scoped_slim_patch",
]

#: Ragsdale & Gravel (2019), as implemented in OutOfAfricaArchaicAdmixture_5R19.
ARCHAIC_EFFECTIVE_SIZE = _ARCHAIC_EFFECTIVE_SIZE
#: East Asian Neanderthal ancestry (Prufer et al. 2014; Vernot & Akey 2015).
ADMIXTURE_PROPORTION = _ADMIXTURE_PROPORTION
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
    split_generations: float = SPLIT_GENERATIONS,
    pulse_generations: float = PULSE_GENERATIONS,
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
    # unable to coalesce for the whole 28,000 generations before it and push
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


# ---------------------------------------------------------------------------
# SLiM building blocks
# ---------------------------------------------------------------------------


def _require_stdpopsim() -> None:
    if stdpopsim.__version__ != REQUIRED_STDPOPSIM_VERSION:
        raise RuntimeError(
            f"run7 requires stdpopsim {REQUIRED_STDPOPSIM_VERSION}; "
            f"found {stdpopsim.__version__}"
        )


def build_stdpopsim_model(repo_root: str | Path, **kwargs) -> stdpopsim.DemographicModel:
    """The grafted demography wrapped for the SLiM engine.

    The census is dropped here: it exists so that ancestry can be read off the
    tree, and the SLiM engine has no use for it -- SLiM marks the focal allele
    itself, which is what the genotype stratification uses.
    """
    _require_stdpopsim()
    kwargs.setdefault("add_census", False)
    grafted = build_eas_with_archaic(repo_root, **kwargs)
    model = stdpopsim.DemographicModel(
        id="PhlashEASArchaicGraft",
        description="PHLASH EAS median history with a grafted Neanderthal branch",
        long_description=(
            "The inferred EAS history carrying the whole human lineage, with a "
            "Neanderthal branch of constant size joining it at the split and a "
            "single admixture pulse into EAS."
        ),
        generation_time=GENERATION_TIME_YEARS,
        mutation_rate=MUTATION_RATE,
        recombination_rate=RECOMBINATION_RATE,
        model=grafted.demography,
    )
    return model


def build_contig() -> stdpopsim.Contig:
    _require_stdpopsim()
    contig = stdpopsim.get_species("HomSap").get_contig(
        length=SEQUENCE_LENGTH_BP,
        mutation_rate=MUTATION_RATE,
        recombination_rate=RECOMBINATION_RATE,
    )
    contig.add_single_site(
        id=FOCAL_SITE_ID,
        coordinate=FOCAL_POSITION_BP,
        description="run7 focal selected site",
    )
    return contig


def snap_generations(value: float) -> float:
    """Round a time onto the tick grid the SLiM engine actually uses."""
    return round(float(value) / SLIM_SCALING_FACTOR) * SLIM_SCALING_FACTOR


def build_extended_events(mode: str, selection_coefficient: float) -> tuple[Any, ...]:
    """Fix the allele in the archaic branch, then select it in EAS after the pulse.

    The two frequency conditions together are the ascertainment: the allele has
    to arrive in EAS at roughly the admixture proportion. Fixing it in the
    archaic branch delivers 2.5% in expectation, but the pulse is a single
    binomial draw, so the realised frequency is required rather than assumed.
    """
    _require_stdpopsim()
    if mode == "neutral":
        return ()
    if mode != "selected":
        raise ValueError(f"mode must be 'selected' or 'neutral', got {mode!r}")

    placement = stdpopsim.GenerationAfter(snap_generations(SPLIT_GENERATIONS))
    onset = snap_generations(PULSE_GENERATIONS)
    # One tick *after* the pulse, as a plain time rather than a GenerationAfter:
    # the condition has to be evaluated once the migration has actually moved
    # archaic genomes into EAS, and a GenerationAfter on both ends of the window
    # collapses it so that the check never runs at all.
    check = onset - SLIM_SCALING_FACTOR
    low, high = TARGET_FREQUENCY_BAND
    events = (
        stdpopsim.DrawMutation(
            time=placement,
            single_site_id=FOCAL_SITE_ID,
            population=ARCHAIC_POPULATION,
        ),
        stdpopsim.ChangeMutationFitness(
            start_time=onset,
            end_time=0.0,
            single_site_id=FOCAL_SITE_ID,
            population=EAS_POPULATION,
            selection_coeff=float(selection_coefficient),
            dominance_coeff=DOMINANCE_COEFFICIENT,
        ),
        stdpopsim.ConditionOnAlleleFrequency(
            start_time=check,
            end_time=check,
            single_site_id=FOCAL_SITE_ID,
            population=EAS_POPULATION,
            op=">=",
            allele_frequency=float(low),
        ),
        stdpopsim.ConditionOnAlleleFrequency(
            start_time=check,
            end_time=check,
            single_site_id=FOCAL_SITE_ID,
            population=EAS_POPULATION,
            op="<=",
            allele_frequency=float(high),
        ),
    )
    if CONDITION_ON_SURVIVAL:
        # Only ever scan tracts that exist: the empirical analysis sees hmmix
        # calls, never the loci where introgression failed to establish. The
        # establishment probability is measured separately rather than lost here.
        events = events + (
            stdpopsim.ConditionOnAlleleFrequency(
                start_time=check,
                end_time=0.0,
                single_site_id=FOCAL_SITE_ID,
                population=EAS_POPULATION,
                op=">",
                allele_frequency=0.0,
            ),
        )
    return events


_ADD_MUT_PATTERN = re.compile(
    r"// Add `mut_type` mutation at `pos`, to a single individual in `pop`\.\n"
    r"function \(void\)add_mut\(object\$ mut_type, object\$ pop, integer\$ pos\) \{"
    r".*?\n\}",
    re.DOTALL,
)
_END_PATTERN = re.compile(
    r"// Output tree sequence file and end the simulation\.\n"
    r"function \(void\)end\(void\) \{.*?\n\}",
    re.DOTALL,
)

#: Fixing the allele in the archaic branch is what makes the pulse deliver it on
#: genuine archaic haplotypes; drawing a single copy instead would hand it a
#: random modern background and defeat the whole point of the graft.
_PLACEMENT = """// run7: fix the focal allele in the archaic branch.
function (void)add_mut(object$ mut_type, object$ pop, integer$ pos) {
    if (size(pop.individuals) == 0)
        err("run7: the archaic population is empty at the placement tick");
    carrier_genomes = pop.genomes;
    carrier_genomes.addNewDrawnMutation(mut_type, pos);
    muts = sim.mutationsOfType(mut_type);
    if (size(muts) != 1)
        err("run7: placement did not create exactly one focal mutation");
    metadata.setValue("run7_focal_mutation_type_id", mut_type.id);
    metadata.setValue("run7_placement_tick", community.tick);
    metadata.setValue("run7_placement_mode", "fixed_in_archaic_population");
    metadata.setValue("run7_placement_carrier_genomes", size(carrier_genomes));
    metadata.setValue("run7_placement_frequency", sim.mutationFrequencies(pop, muts[0]));
}"""


def _census_end_function(target_population_id: int) -> str:
    return f"""// run7: record the present-day census frequency, then finish.
function (void)end(void) {{
    mt_id = metadata.getValue("run7_focal_mutation_type_id");
    pop = sim.subpopulations[sim.subpopulations.id == {target_population_id}];
    if (size(pop) != 1)
        err("run7: the target subpopulation is not present at the end");
    pop = pop[0];
    total_count = size(pop.genomes);
    alt_count = 0;
    state = "absent";
    if (!isNULL(mt_id)) {{
        mt_matches = sim.mutationTypes[sim.mutationTypes.id == mt_id];
        if (size(mt_matches) == 1) {{
            mt = mt_matches[0];
            muts = sim.mutationsOfType(mt);
            subs = sim.substitutions[sim.substitutions.mutationType == mt];
            if (size(subs) == 1) {{
                alt_count = total_count;
                state = "fixed_as_substitution";
            }} else if (size(muts) == 1) {{
                alt_count = sum(pop.genomes.containsMutations(muts[0]));
                state = "segregating_or_fixed";
            }} else {{
                alt_count = 0;
                state = "lost";
            }}
        }}
    }}
    metadata.setValue("run7_final_census_alt_count", alt_count);
    metadata.setValue("run7_final_census_total_count", total_count);
    metadata.setValue("run7_final_census_af", alt_count / total_count);
    metadata.setValue("run7_final_census_state", state);
    sim.treeSeqOutput(trees_file, metadata=metadata);
    sim.simulationFinished();
}}"""


def _replace_once(text: str, pattern: re.Pattern[str], replacement: str, label: str) -> str:
    matches = pattern.findall(text)
    if len(matches) != 1:
        raise RuntimeError(
            f"stdpopsim {label} contract changed: expected one match, got {len(matches)}"
        )
    return pattern.sub(lambda _: replacement, text, count=1)


@contextmanager
def scoped_slim_patch(
    mode: str, target_population_id: int = 0, *, patch_placement: bool = True
) -> Iterator[dict[str, Any]]:
    """Patch the SLiM preamble for the duration of one simulation.

    ``patch_placement`` is False for the de novo arms: stdpopsim already draws a
    single copy, which is exactly what a de novo origin means, so only the census
    at the end needs replacing.
    """
    original = slim_engine._slim_functions
    record: dict[str, Any] = {"mode": mode, "patch_placement": patch_placement}
    if mode == "neutral":
        record["patched"] = False
        yield record
        return
    patched = str(original)
    if patch_placement:
        patched = _replace_once(patched, _ADD_MUT_PATTERN, _PLACEMENT, "add_mut")
    patched = _replace_once(
        patched, _END_PATTERN, _census_end_function(target_population_id), "end"
    )
    record.update(
        {
            "patched": True,
            "target_population_id": target_population_id,
            "patched_sha256": hashlib.sha256(patched.encode()).hexdigest(),
        }
    )
    slim_engine._slim_functions = patched
    try:
        yield record
    finally:
        slim_engine._slim_functions = original


def build_denovo_events(selection_coefficient: float) -> tuple[Any, ...]:
    """A single de novo copy in EAS at the pulse time, conditioned on survival.

    Everything except the origin is held identical to the introgressed arms --
    same model, same onset, same null -- so the contrast isolates one variable:
    a sweep starting from one copy on a random modern background against one
    starting from ~39 copies carried on archaic haplotypes.

    Survival conditioning is doing far more work here. A single copy establishes
    with probability about ``2 h s``, so most draws are lost and SLiM restarts;
    the restart is cheap only because the checkpoint sits at the draw itself
    rather than back at the archaic split.
    """
    _require_stdpopsim()
    onset = snap_generations(PULSE_GENERATIONS)
    check = stdpopsim.GenerationAfter(onset)
    return (
        stdpopsim.DrawMutation(
            time=onset, single_site_id=FOCAL_SITE_ID, population=EAS_POPULATION
        ),
        stdpopsim.ChangeMutationFitness(
            start_time=onset,
            end_time=0.0,
            single_site_id=FOCAL_SITE_ID,
            population=EAS_POPULATION,
            selection_coeff=float(selection_coefficient),
            dominance_coeff=DOMINANCE_COEFFICIENT,
        ),
        stdpopsim.ConditionOnAlleleFrequency(
            start_time=check,
            end_time=0.0,
            single_site_id=FOCAL_SITE_ID,
            population=EAS_POPULATION,
            op=">",
            allele_frequency=0.0,
        ),
    )
