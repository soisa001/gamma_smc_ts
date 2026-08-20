"""SLiM model for run6's EAS arm.

The full PHLASH history is used, unflattened.  msprime's sweep model rejects
demographic events during a sweep, which forced run5 to hold the EAS history
constant across the sweep window; SLiM has no such restriction, so the
ten-thousand-epoch curve is carried as-is.  SLiM also marks the selected allele,
which msprime's structured-coalescent sweep does not, and the mark is what makes
the genotype stratification possible at all.
"""

from __future__ import annotations

import hashlib
import re
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator

import numpy as np
import stdpopsim
from stdpopsim import slim_engine

from .eas_sweep_models import build_eas_demography_models, load_phlash_eas_npz
from .run6_config import (
    DOMINANCE_COEFFICIENT,
    SELECTION_COEFFICIENT,
    EAS_ONSET_GENERATIONS,
    EAS_PLACEMENT_FREQUENCY,
    EAS_POPULATION,
    FOCAL_POSITION_BP,
    FOCAL_SITE_ID,
    MUTATION_RATE,
    PHLASH_EAS_RELATIVE_PATH,
    PHLASH_EAS_SHA256,
    RECOMBINATION_RATE,
    REQUIRED_STDPOPSIM_VERSION,
    SEQUENCE_LENGTH_BP,
    SLIM_SCALING_FACTOR,
)

__all__ = [
    "build_contig",
    "build_eas_model",
    "build_extended_events",
    "onset_generations",
    "scoped_slim_patch",
]


def _require_stdpopsim() -> None:
    if stdpopsim.__version__ != REQUIRED_STDPOPSIM_VERSION:
        raise RuntimeError(
            f"run6 requires stdpopsim {REQUIRED_STDPOPSIM_VERSION}; "
            f"found {stdpopsim.__version__}"
        )


def onset_generations() -> float:
    return round(EAS_ONSET_GENERATIONS / SLIM_SCALING_FACTOR) * SLIM_SCALING_FACTOR


def build_eas_model(repo_root: str | Path):
    """The PHLASH EAS median history, full resolution."""
    _require_stdpopsim()
    artifact = load_phlash_eas_npz(
        Path(repo_root) / PHLASH_EAS_RELATIVE_PATH, expected_sha256=PHLASH_EAS_SHA256
    )
    return build_eas_demography_models(artifact)["median"].stdpopsim_model


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
        description="run6 focal selected site",
    )
    return contig


def build_extended_events(mode: str) -> tuple[Any, ...]:
    """Selected: place at 2% and select. Neutral: no allele at all."""
    _require_stdpopsim()
    if mode == "neutral":
        return ()
    if mode != "selected":
        raise ValueError(f"mode must be 'selected' or 'neutral', got {mode!r}")
    onset = onset_generations()
    return (
        stdpopsim.DrawMutation(
            time=onset, single_site_id=FOCAL_SITE_ID, population=EAS_POPULATION
        ),
        stdpopsim.ChangeMutationFitness(
            start_time=onset,
            end_time=0.0,
            single_site_id=FOCAL_SITE_ID,
            population=EAS_POPULATION,
            selection_coeff=SELECTION_COEFFICIENT,
            dominance_coeff=DOMINANCE_COEFFICIENT,
        ),
        # Standing variation at 2% is very unlikely to be lost at s = 0.01, but
        # conditioning costs nothing here because the checkpoint sits at the
        # placement tick rather than deep in the past.
        stdpopsim.ConditionOnAlleleFrequency(
            start_time=stdpopsim.GenerationAfter(onset),
            end_time=0.0,
            single_site_id=FOCAL_SITE_ID,
            population=EAS_POPULATION,
            op=">",
            allele_frequency=0.0,
        ),
    )


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

_PLACEMENT = f"""// run6: place the allele as standing variation.
function (void)add_mut(object$ mut_type, object$ pop, integer$ pos) {{
    all_genomes = pop.genomes;
    n_target = asInteger(round({EAS_PLACEMENT_FREQUENCY:.10g} * size(all_genomes)));
    if (n_target < 1)
        err("run6: the placement target is fewer than one genome");
    carrier_genomes = sample(all_genomes, n_target);
    carrier_genomes.addNewDrawnMutation(mut_type, pos);
    muts = sim.mutationsOfType(mut_type);
    if (size(muts) != 1)
        err("run6: placement did not create exactly one focal mutation");
    metadata.setValue("run6_focal_mutation_type_id", mut_type.id);
    metadata.setValue("run6_placement_tick", community.tick);
    metadata.setValue("run6_placement_mode", "uniform_random_genomes");
    metadata.setValue("run6_placement_frequency", sim.mutationFrequencies(pop, muts[0]));
}}"""

_CENSUS_END = """// run6: record the present-day census frequency, then finish.
function (void)end(void) {
    mt_id = metadata.getValue("run6_focal_mutation_type_id");
    pop = sim.subpopulations[0];
    total_count = size(pop.genomes);
    alt_count = 0;
    state = "absent";
    if (!isNULL(mt_id)) {
        mt_matches = sim.mutationTypes[sim.mutationTypes.id == mt_id];
        if (size(mt_matches) == 1) {
            mt = mt_matches[0];
            muts = sim.mutationsOfType(mt);
            subs = sim.substitutions[sim.substitutions.mutationType == mt];
            if (size(subs) == 1) {
                alt_count = total_count;
                state = "fixed_as_substitution";
            } else if (size(muts) == 1) {
                alt_count = sum(pop.genomes.containsMutations(muts[0]));
                state = "segregating_or_fixed";
            } else {
                alt_count = 0;
                state = "lost";
            }
        }
    }
    metadata.setValue("run6_final_census_alt_count", alt_count);
    metadata.setValue("run6_final_census_total_count", total_count);
    metadata.setValue("run6_final_census_af", alt_count / total_count);
    metadata.setValue("run6_final_census_state", state);
    sim.treeSeqOutput(trees_file, metadata=metadata);
    sim.simulationFinished();
}"""


def _replace_once(text: str, pattern: re.Pattern[str], replacement: str, label: str) -> str:
    matches = pattern.findall(text)
    if len(matches) != 1:
        raise RuntimeError(
            f"stdpopsim {label} contract changed: expected one match, got {len(matches)}"
        )
    return pattern.sub(lambda _: replacement, text, count=1)


@contextmanager
def scoped_slim_patch(mode: str) -> Iterator[dict[str, Any]]:
    original = slim_engine._slim_functions
    record: dict[str, Any] = {"mode": mode}
    if mode == "neutral":
        record["patched"] = False
        yield record
        return
    patched = _replace_once(str(original), _ADD_MUT_PATTERN, _PLACEMENT, "add_mut")
    patched = _replace_once(patched, _END_PATTERN, _CENSUS_END, "end")
    record.update(
        {"patched": True, "patched_sha256": hashlib.sha256(patched.encode()).hexdigest()}
    )
    slim_engine._slim_functions = patched
    try:
        yield record
    finally:
        slim_engine._slim_functions = original
