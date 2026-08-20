"""Demography, contig, events and SLiM patches for the run3 study.

Two things differ substantially from run2.

**The allele marks archaic ancestry.**  ``add_mut`` is replaced so that the
focal mutation is placed on *every* genome of the Neanderthal population one
tick after the archaic split.  From then on it is carried into CHB by the
catalog's own continuous Neanderthal migration, so its frequency in CHB is
exactly the archaic ancestry fraction at the focal base.  This happens in both
modes; only the fitness callback distinguishes them.  The neutral arm therefore
contains introgression *without* selection, which is the null the question
actually calls for.

**Timing needs no ancestral fudge.**  In ``OutOfAfricaArchaicAdmixture_5R19``
CHB is founded at the same moment its Neanderthal migration switches on, so
"selection starts when introgression starts" is expressible exactly, in CHB
alone, with no ancestral population carrying the allele under selection first.

Tick arithmetic mirrors stdpopsim's own, as in run2, but the generation time is
now per-arm: 29 years for the catalog model, 25 for the PHLASH EAS history.
"""

from __future__ import annotations

import hashlib
import io
import re
from contextlib import contextmanager, redirect_stdout
from pathlib import Path
from typing import Any, Iterator

import numpy as np
import stdpopsim
from stdpopsim import slim_engine

from .eas_sweep_models import (
    build_eas_demography_models,
    load_phlash_eas_npz,
    serialize_extended_events,
)
from .run2_models import (  # generic, arm-independent guards
    mask_focal_rate_map,
    scoped_focal_overlay_patch,
)
from .run3_config import (
    ARCHAIC_POPULATION,
    CHB_ARCHAIC_MIGRATION_END_GENERATIONS,
    CHB_SPLIT_GENERATIONS,
    DOMINANCE_COEFFICIENT,
    EAS_PLACEMENT_FREQUENCY,
    FOCAL_POSITION_BP,
    FOCAL_SITE_ID,
    MUTATION_RATE,
    NEANDERTHAL_SPLIT_GENERATIONS,
    OOA_GENERATION_TIME,
    OOA_MODEL_ID,
    PHLASH_EAS_RELATIVE_PATH,
    PHLASH_EAS_SHA256,
    RECOMBINATION_RATE,
    REQUIRED_STDPOPSIM_VERSION,
    Run3Arm,
    SAMPLE_DIPLOIDS,
    SELECTION_COEFFICIENT,
    SEQUENCE_LENGTH_BP,
    SLIM_BURN_IN,
    SLIM_SCALING_FACTOR,
)

__all__ = [
    "build_contig",
    "build_demographic_model",
    "build_extended_events",
    "eas_selection_onset_generations",
    "epoch_tick_offset",
    "event_tick_offset",
    "generate_slim_script",
    "mask_focal_rate_map",
    "model_record",
    "population_id",
    "scoped_focal_overlay_patch",
    "scoped_slim_patch",
    "selection_onset_generations",
    "snap_generations",
    "tick_schedule",
]


# ---------------------------------------------------------------------------
# Tick arithmetic
# ---------------------------------------------------------------------------


def snap_generations(time: Any, scaling_factor: float = SLIM_SCALING_FACTOR) -> float:
    """Snap an extended-event time onto the SLiM tick grid (stdpopsim fix_time)."""
    snapped = round(float(time) / scaling_factor) * scaling_factor
    if isinstance(time, stdpopsim.GenerationAfter):
        snapped -= scaling_factor
    if snapped < 0:
        raise ValueError(f"snapped time {snapped} is negative")
    return float(snapped)


def event_tick_offset(
    time: Any,
    generation_time: float,
    *,
    scaling_factor: float = SLIM_SCALING_FACTOR,
) -> int:
    """Ticks before the final tick for an extended-event time."""
    years = snap_generations(time, scaling_factor) * generation_time
    return int(years / generation_time / scaling_factor)


def epoch_tick_offset(
    generations_ago: float,
    generation_time: float,
    *,
    scaling_factor: float = SLIM_SCALING_FACTOR,
) -> int:
    """Ticks before the final tick for a demographic epoch boundary."""
    years = round(float(generations_ago) * generation_time)
    return int(years / generation_time / scaling_factor)


def selection_onset_generations() -> float:
    """Realized CHB selection onset, in generations, after tick snapping."""
    return snap_generations(CHB_SPLIT_GENERATIONS)


def eas_selection_onset_generations(arm: Run3Arm) -> float:
    """Realized EAS onset, matched to CHB in *years* then snapped.

    The two models declare different generation times (29 vs 25), so matching in
    years keeps the year-denominated TMRCA thresholds comparable across arms.
    Matching in generations instead would put the EAS onset 5,000 years later.
    """
    years = selection_onset_generations() * OOA_GENERATION_TIME
    return snap_generations(years / arm.generation_time)


def tick_schedule(arm: Run3Arm) -> dict[str, Any]:
    """Return the nominal and realized tick schedule, asserting the invariants."""

    q = SLIM_SCALING_FACTOR
    gt = arm.generation_time
    schedule: dict[str, Any] = {
        "scaling_factor": q,
        "generation_time_years": gt,
        "tick_offsets_before_final_tick": {"present": 0},
        "nominal_generations_ago": {"present": 0.0},
        "realized_generations_ago": {"present": 0.0},
        "realized_years_ago": {"present": 0.0},
    }

    if arm.archaic:
        split_offset = epoch_tick_offset(NEANDERTHAL_SPLIT_GENERATIONS, gt)
        fixation_time = stdpopsim.GenerationAfter(NEANDERTHAL_SPLIT_GENERATIONS)
        fixation_offset = event_tick_offset(fixation_time, gt)
        if fixation_offset != split_offset - 1:
            raise RuntimeError(
                "archaic fixation must fall exactly one tick after the split: "
                f"split {split_offset}, fixation {fixation_offset}"
            )
        chb_split_offset = epoch_tick_offset(CHB_SPLIT_GENERATIONS, gt)
        onset_offset = event_tick_offset(CHB_SPLIT_GENERATIONS, gt)
        if onset_offset != chb_split_offset:
            raise RuntimeError(
                "selection must start on the tick CHB is founded and archaic "
                f"migration begins: founding {chb_split_offset}, onset {onset_offset}"
            )
        migration_end_offset = epoch_tick_offset(
            CHB_ARCHAIC_MIGRATION_END_GENERATIONS, gt
        )
        if not split_offset > chb_split_offset > migration_end_offset > 0:
            raise RuntimeError("the archaic event schedule is out of order")

        onset = snap_generations(CHB_SPLIT_GENERATIONS)
        schedule["tick_offsets_before_final_tick"].update(
            {
                "archaic_split": split_offset,
                "archaic_fixation": fixation_offset,
                "chb_founding_and_selection_onset": chb_split_offset,
                "archaic_migration_end": migration_end_offset,
            }
        )
        schedule["nominal_generations_ago"].update(
            {
                "archaic_split": NEANDERTHAL_SPLIT_GENERATIONS,
                "archaic_fixation": NEANDERTHAL_SPLIT_GENERATIONS,
                "chb_founding_and_selection_onset": CHB_SPLIT_GENERATIONS,
                "archaic_migration_end": CHB_ARCHAIC_MIGRATION_END_GENERATIONS,
            }
        )
        schedule["realized_generations_ago"].update(
            {
                "archaic_split": split_offset * q,
                "archaic_fixation": fixation_offset * q,
                "chb_founding_and_selection_onset": onset,
                "archaic_migration_end": migration_end_offset * q,
            }
        )
        schedule["selection_intervals_generations_ago"] = [
            {"population": arm.target_population, "start": onset, "end": 0.0}
        ]
    else:
        onset = eas_selection_onset_generations(arm)
        onset_offset = event_tick_offset(onset, gt)
        schedule["tick_offsets_before_final_tick"]["selection_onset"] = onset_offset
        schedule["nominal_generations_ago"]["selection_onset"] = onset
        schedule["realized_generations_ago"]["selection_onset"] = onset_offset * q
        schedule["selection_intervals_generations_ago"] = [
            {"population": arm.target_population, "start": onset, "end": 0.0}
        ]

    schedule["realized_years_ago"] = {
        key: value * gt for key, value in schedule["realized_generations_ago"].items()
    }
    return schedule


# ---------------------------------------------------------------------------
# Demography and contig
# ---------------------------------------------------------------------------


def _require_stdpopsim() -> None:
    if stdpopsim.__version__ != REQUIRED_STDPOPSIM_VERSION:
        raise RuntimeError(
            f"run3 requires stdpopsim {REQUIRED_STDPOPSIM_VERSION}; "
            f"found {stdpopsim.__version__}"
        )


def _load_ooa_archaic_model() -> stdpopsim.DemographicModel:
    """Load ``OutOfAfricaArchaicAdmixture_5R19`` and contract-check its schedule."""

    _require_stdpopsim()
    model = stdpopsim.get_species("HomSap").get_demographic_model(OOA_MODEL_ID)
    if not np.isclose(model.generation_time, OOA_GENERATION_TIME):
        raise RuntimeError(
            f"{OOA_MODEL_ID} generation time changed from {OOA_GENERATION_TIME}"
        )
    names = {population.name for population in model.model.populations}
    for required in (ARCHAIC_POPULATION, "CHB", "CEU", "YRI"):
        if required not in names:
            raise RuntimeError(f"{OOA_MODEL_ID} is missing population {required!r}")

    events = model.model.events
    archaic_id = population_id(model, ARCHAIC_POPULATION)
    chb_id = population_id(model, "CHB")

    splits = [
        e
        for e in events
        if type(e).__name__ == "MassMigration"
        and e.source == archaic_id
        and float(e.proportion) == 1.0
    ]
    if len(splits) != 1 or not np.isclose(
        float(splits[0].time), NEANDERTHAL_SPLIT_GENERATIONS
    ):
        raise RuntimeError(f"{OOA_MODEL_ID} archaic split time changed")

    chb_splits = [
        e
        for e in events
        if type(e).__name__ == "MassMigration"
        and e.source == chb_id
        and float(e.proportion) == 1.0
    ]
    if len(chb_splits) != 1 or not np.isclose(
        float(chb_splits[0].time), CHB_SPLIT_GENERATIONS
    ):
        raise RuntimeError(f"{OOA_MODEL_ID} CHB split time changed")

    archaic_into_chb = [
        e
        for e in events
        if type(e).__name__ == "MigrationRateChange"
        and getattr(e, "matrix_index", None) == (chb_id, archaic_id)
        and float(e.rate) > 0
    ]
    if not archaic_into_chb:
        raise RuntimeError(f"{OOA_MODEL_ID} has no CHB-to-Neanderthal migration")

    model.mutation_rate = MUTATION_RATE
    model.recombination_rate = RECOMBINATION_RATE
    return model


def build_demographic_model(
    arm: Run3Arm, repo_root: str | Path
) -> stdpopsim.DemographicModel:
    _require_stdpopsim()
    if arm.archaic:
        return _load_ooa_archaic_model()
    artifact = load_phlash_eas_npz(
        Path(repo_root) / PHLASH_EAS_RELATIVE_PATH, expected_sha256=PHLASH_EAS_SHA256
    )
    return build_eas_demography_models(artifact)["median"].stdpopsim_model


def build_contig() -> stdpopsim.Contig:
    """Return the 10 Mb contig with the focal base reserved as a single site."""
    _require_stdpopsim()
    species = stdpopsim.get_species("HomSap")
    contig = species.get_contig(
        length=SEQUENCE_LENGTH_BP,
        mutation_rate=MUTATION_RATE,
        recombination_rate=RECOMBINATION_RATE,
    )
    contig.add_single_site(
        id=FOCAL_SITE_ID,
        coordinate=FOCAL_POSITION_BP,
        description="run3 focal selected site",
        long_description=(
            "Reserved base. In the archaic arm the allele is fixed in "
            "Neanderthal after the split and marks archaic ancestry exactly."
        ),
    )
    return contig


def population_id(model: stdpopsim.DemographicModel, name: str) -> int:
    for index, population in enumerate(model.model.populations):
        if population.name == name:
            return index
    raise ValueError(f"population {name!r} is absent from the model")


def build_extended_events(arm: Run3Arm, mode: str) -> tuple[Any, ...]:
    """Return the extended events for one arm and mode.

    The mutation is drawn in **both** modes, so the neutral arm carries the same
    introgressed material without any fitness effect.  Only the
    ``ChangeMutationFitness`` event distinguishes the modes, and there are no
    conditioning events at all.
    """
    _require_stdpopsim()
    if mode not in ("selected", "neutral"):
        raise ValueError(f"mode must be 'selected' or 'neutral', got {mode!r}")
    tick_schedule(arm)

    if arm.archaic:
        draw_time: Any = stdpopsim.GenerationAfter(NEANDERTHAL_SPLIT_GENERATIONS)
        onset: Any = CHB_SPLIT_GENERATIONS
    else:
        draw_time = eas_selection_onset_generations(arm)
        onset = draw_time

    events: list[Any] = [
        stdpopsim.DrawMutation(
            time=draw_time,
            single_site_id=FOCAL_SITE_ID,
            population=arm.origin_population,
        )
    ]
    if mode == "selected":
        events.append(
            stdpopsim.ChangeMutationFitness(
                start_time=onset,
                end_time=0.0,
                single_site_id=FOCAL_SITE_ID,
                population=arm.target_population,
                selection_coeff=SELECTION_COEFFICIENT,
                dominance_coeff=DOMINANCE_COEFFICIENT,
            )
        )
    return tuple(events)


# ---------------------------------------------------------------------------
# SLiM patches
# ---------------------------------------------------------------------------

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


def _placement_function(arm: Run3Arm) -> str:
    if arm.archaic:
        selection = """    carrier_genomes = pop.genomes;
    if (size(carrier_genomes) == 0)
        err("run3: the archaic population is empty at the fixation tick");"""
        placement = "fixed_in_archaic_population"
    else:
        selection = f"""    all_genomes = pop.genomes;
    n_target = asInteger(round({EAS_PLACEMENT_FREQUENCY:.10g} * size(all_genomes)));
    if (n_target < 1)
        err("run3: placement target is fewer than one genome");
    carrier_genomes = sample(all_genomes, n_target);"""
        placement = "uniform_random_genomes"

    return f"""// run3: introduce the focal allele (fixed in the archaic population, or at a
// set frequency in the single-population arm).
function (void)add_mut(object$ mut_type, object$ pop, integer$ pos) {{
    inds = pop.individuals;
    if (size(inds) == 0)
        err("run3: the origin population is empty at the placement tick");
{selection}
    carrier_genomes.addNewDrawnMutation(mut_type, pos);
    muts = sim.mutationsOfType(mut_type);
    if (size(muts) != 1)
        err("run3: placement did not create exactly one focal mutation");
    frequency = sim.mutationFrequencies(pop, muts[0]);
    metadata.setValue("run3_focal_mutation_type_id", mut_type.id);
    metadata.setValue("run3_placement_tick", community.tick);
    metadata.setValue("run3_placement_population_id", pop.id);
    metadata.setValue("run3_placement_mode", "{placement}");
    metadata.setValue("run3_placement_carrier_genomes", size(carrier_genomes));
    metadata.setValue("run3_placement_population_genomes", size(pop.genomes));
    metadata.setValue("run3_placement_frequency", frequency);
}}"""


def _census_end_function(target_population_id: int) -> str:
    return f"""// run3: record the present-day census frequency, then finish.
function (void)end(void) {{
    mt_id = metadata.getValue("run3_focal_mutation_type_id");
    mt_matches = sim.mutationTypes[sim.mutationTypes.id == mt_id];
    if (size(mt_matches) != 1)
        err("run3: the focal mutation type is not unique at the end");
    mt = mt_matches[0];
    pop_matches = sim.subpopulations[sim.subpopulations.id == {target_population_id}];
    if (size(pop_matches) != 1)
        err("run3: the target population is not unique at the end");
    pop = pop_matches[0];
    total_count = size(pop.genomes);
    if (total_count <= 0)
        err("run3: the target population is empty at the end");
    muts = sim.mutationsOfType(mt);
    subs = sim.substitutions[sim.substitutions.mutationType == mt];
    if (size(muts) + size(subs) > 1)
        err("run3: more than one focal allele is present at the end");
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
    metadata.setValue("run3_final_census_alt_count", alt_count);
    metadata.setValue("run3_final_census_total_count", total_count);
    metadata.setValue("run3_final_census_af", alt_count / total_count);
    metadata.setValue("run3_final_census_state", state);
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
    arm: Run3Arm, mode: str, target_population_id: int
) -> Iterator[dict[str, Any]]:
    """Patch stdpopsim's SLiM helpers for one engine call, in both modes."""

    original_object = slim_engine._slim_functions
    original = str(original_object)
    placement = _placement_function(arm)
    end = _census_end_function(target_population_id)
    patched = _replace_once(original, _ADD_MUT_PATTERN, placement, "add_mut")
    patched = _replace_once(patched, _END_PATTERN, end, "end")
    record = {
        "arm_id": arm.arm_id,
        "mode": mode,
        "target_population_id": target_population_id,
        "original_sha256": hashlib.sha256(original.encode()).hexdigest(),
        "placement_sha256": hashlib.sha256(placement.encode()).hexdigest(),
        "end_sha256": hashlib.sha256(end.encode()).hexdigest(),
        "patched_sha256": hashlib.sha256(patched.encode()).hexdigest(),
    }
    slim_engine._slim_functions = patched
    try:
        yield record
    finally:
        slim_engine._slim_functions = original_object


# ---------------------------------------------------------------------------
# Script generation
# ---------------------------------------------------------------------------


def generate_slim_script(arm: Run3Arm, mode: str, repo_root: str | Path) -> str:
    """Return the SLiM script stdpopsim would run, with run3's patches applied."""
    model = build_demographic_model(arm, repo_root)
    contig = build_contig()
    events = build_extended_events(arm, mode)
    target_id = population_id(model, arm.target_population)
    engine = stdpopsim.get_engine("slim")
    buffer = io.StringIO()
    with scoped_slim_patch(arm, mode, target_id):
        with redirect_stdout(buffer):
            engine.simulate(
                model,
                contig,
                {arm.target_population: SAMPLE_DIPLOIDS},
                seed=1,
                extended_events=list(events),
                slim_script=True,
                slim_scaling_factor=SLIM_SCALING_FACTOR,
                slim_burn_in=SLIM_BURN_IN,
            )
    return buffer.getvalue()


def model_record(arm: Run3Arm, mode: str, repo_root: str | Path) -> dict[str, Any]:
    model = build_demographic_model(arm, repo_root)
    events = build_extended_events(arm, mode)
    return {
        "schema": "gamma-smc.run3-model/v1",
        "arm_id": arm.arm_id,
        "mode": mode,
        "demographic_model_id": model.id,
        "generation_time_years": arm.generation_time,
        "target_population": arm.target_population,
        "target_population_id": population_id(model, arm.target_population),
        "origin_population": arm.origin_population,
        "placement": arm.placement,
        "allele_present": True,
        "selection_coefficient": SELECTION_COEFFICIENT if mode == "selected" else 0.0,
        "conditioning_events": 0,
        "tick_schedule": tick_schedule(arm),
        "extended_events": serialize_extended_events(events),
        "stdpopsim_version": stdpopsim.__version__,
    }
