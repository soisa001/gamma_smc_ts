"""Demography, contig, events and SLiM patches for the run2 study.

Two things here are worth reading carefully.

**Standing variation without rejection sampling.**  stdpopsim's ``DrawMutation``
puts a single copy of the focal allele on one random genome.  run2 needs the
allele to start at exactly the introgression fraction, so this module replaces
stdpopsim's ``add_mut`` SLiM function for the duration of one engine call:

* CHB arm -- the allele goes on *every* introgressing genome at the tick after
  the pulse.  That is exactly equivalent to the allele being fixed in the
  archaic source, so the post-pulse recipient frequency is the pulse proportion
  and the carriers are exactly the archaic haplotypes.
* EAS arm -- there is no archaic source, so the allele goes on a uniformly
  random 2.96% of genomes at the same tick.

Neither placement rejects anything, so the accept rate is 1.

**Tick arithmetic.**  With ``Q = 5`` SLiM only has ticks every five
generations, and stdpopsim snaps event times onto that grid.  The functions
below mirror stdpopsim's own arithmetic exactly so that the realized times can
be asserted rather than assumed:

* demographic epoch times use ``int(round(t * generation_time) / generation_time / Q)``
* extended-event times are first snapped with ``round(t / Q) * Q``, and
  ``GenerationAfter`` subtracts one further ``Q``

which puts the 2272-generation pulse on tick ``G0 - 454`` (2270 generations
ago), the standing-variation placement one tick later on ``G0 - 453`` (2265
generations ago), and the 2016-generation Han split on ``G0 - 403`` (2015
generations ago).
"""

from __future__ import annotations

import hashlib
import inspect
import io
import re
from contextlib import contextmanager, redirect_stdout
from pathlib import Path
from typing import Any, Iterator

import msprime
import numpy as np
import stdpopsim
from stdpopsim import slim_engine

from .eas_sweep_models import (
    build_eas_demography_models,
    load_ancient_eurasia_model,
    load_phlash_eas_npz,
    serialize_extended_events,
)
from .run2_config import (
    DOMINANCE_COEFFICIENT,
    FOCAL_POSITION_BP,
    FOCAL_SITE_ID,
    GENERATION_TIME_YEARS,
    HAN_SPLIT_GENERATIONS,
    MUTATION_RATE,
    PHLASH_EAS_RELATIVE_PATH,
    PHLASH_EAS_SHA256,
    PULSE_GENERATIONS,
    RECOMBINATION_RATE,
    REQUIRED_STDPOPSIM_VERSION,
    Run2Arm,
    SAMPLE_DIPLOIDS,
    SELECTION_COEFFICIENT,
    SEQUENCE_LENGTH_BP,
    SLIM_BURN_IN,
    SLIM_SCALING_FACTOR,
    STANDING_FREQUENCY,
)

# ---------------------------------------------------------------------------
# Tick arithmetic (mirrors stdpopsim.slim_engine)
# ---------------------------------------------------------------------------


def snap_generations(time: Any, scaling_factor: float = SLIM_SCALING_FACTOR) -> float:
    """Snap an extended-event time onto the SLiM tick grid.

    Mirrors ``slim_makescript.fix_time``: ``round(t / Q) * Q``, with a further
    ``- Q`` for :class:`stdpopsim.GenerationAfter`.
    """
    snapped = round(float(time) / scaling_factor) * scaling_factor
    if isinstance(time, stdpopsim.GenerationAfter):
        snapped -= scaling_factor
    if snapped < 0:
        raise ValueError(f"snapped time {snapped} is negative")
    return float(snapped)


def event_tick_offset(
    time: Any,
    *,
    scaling_factor: float = SLIM_SCALING_FACTOR,
    generation_time: float = GENERATION_TIME_YEARS,
) -> int:
    """Ticks before the final tick for an extended-event time."""
    years = snap_generations(time, scaling_factor) * generation_time
    return int(years / generation_time / scaling_factor)


def epoch_tick_offset(
    generations_ago: float,
    *,
    scaling_factor: float = SLIM_SCALING_FACTOR,
    generation_time: float = GENERATION_TIME_YEARS,
) -> int:
    """Ticks before the final tick for a demographic epoch boundary.

    Epoch times are *not* snapped by ``fix_time``; stdpopsim converts them to
    years with ``round(t * generation_time)`` and truncates in ``time_to_tick``.
    """
    years = round(float(generations_ago) * generation_time)
    return int(years / generation_time / scaling_factor)


def realized_generations(tick_offset: int, scaling_factor: float = SLIM_SCALING_FACTOR) -> float:
    """Return the generations-ago actually represented by a tick offset."""
    return float(tick_offset) * float(scaling_factor)


def tick_schedule(arm: Run2Arm) -> dict[str, Any]:
    """Return the nominal and realized tick schedule for one arm.

    Raises if the standing-variation placement does not land exactly one tick
    after the pulse, which is the invariant the whole design rests on.
    """
    standing_time = stdpopsim.GenerationAfter(PULSE_GENERATIONS)
    pulse_offset = epoch_tick_offset(PULSE_GENERATIONS)
    standing_offset = event_tick_offset(standing_time)

    if standing_offset != pulse_offset - 1:
        raise RuntimeError(
            "standing-variation placement must fall exactly one tick after the "
            f"pulse: pulse offset {pulse_offset}, placement offset {standing_offset}"
        )

    schedule: dict[str, Any] = {
        "scaling_factor": SLIM_SCALING_FACTOR,
        "generation_time_years": GENERATION_TIME_YEARS,
        "tick_offsets_before_final_tick": {
            "pulse": pulse_offset,
            "standing_variation": standing_offset,
            "present": 0,
        },
        "nominal_generations_ago": {
            "pulse": PULSE_GENERATIONS,
            "standing_variation": PULSE_GENERATIONS,
            "present": 0.0,
        },
        "realized_generations_ago": {
            "pulse": realized_generations(pulse_offset),
            "standing_variation": realized_generations(standing_offset),
            "present": 0.0,
        },
        "selection_intervals_generations_ago": [],
    }

    if arm.introgression:
        split_offset = epoch_tick_offset(HAN_SPLIT_GENERATIONS)
        split_event_offset = event_tick_offset(HAN_SPLIT_GENERATIONS)
        if split_event_offset != split_offset:
            raise RuntimeError(
                "the Han fitness callback does not start on the split tick: "
                f"split offset {split_offset}, callback offset {split_event_offset}"
            )
        if not standing_offset > split_offset > 0:
            raise RuntimeError("Han split must fall between the pulse and the present")
        schedule["tick_offsets_before_final_tick"]["han_split"] = split_offset
        schedule["nominal_generations_ago"]["han_split"] = HAN_SPLIT_GENERATIONS
        schedule["realized_generations_ago"]["han_split"] = realized_generations(
            split_offset
        )
        schedule["selection_intervals_generations_ago"] = [
            {
                "population": arm.origin_population,
                "start": realized_generations(standing_offset),
                "end": realized_generations(split_offset),
            },
            {
                "population": arm.target_population,
                "start": realized_generations(split_offset),
                "end": 0.0,
            },
        ]
    else:
        schedule["selection_intervals_generations_ago"] = [
            {
                "population": arm.target_population,
                "start": realized_generations(standing_offset),
                "end": 0.0,
            }
        ]
    return schedule


# ---------------------------------------------------------------------------
# Demography and contig
# ---------------------------------------------------------------------------


def _require_stdpopsim() -> None:
    if stdpopsim.__version__ != REQUIRED_STDPOPSIM_VERSION:
        raise RuntimeError(
            f"run2 requires stdpopsim {REQUIRED_STDPOPSIM_VERSION}; "
            f"found {stdpopsim.__version__}"
        )


def phlash_artifact_path(repo_root: str | Path) -> Path:
    return Path(repo_root) / PHLASH_EAS_RELATIVE_PATH


def build_demographic_model(
    arm: Run2Arm, repo_root: str | Path
) -> stdpopsim.DemographicModel:
    """Return the stdpopsim demographic model backing one arm."""

    _require_stdpopsim()
    if arm.introgression:
        return load_ancient_eurasia_model()
    artifact = load_phlash_eas_npz(
        phlash_artifact_path(repo_root), expected_sha256=PHLASH_EAS_SHA256
    )
    return build_eas_demography_models(artifact)["median"].stdpopsim_model


def build_contig() -> stdpopsim.Contig:
    """Return the 10 Mb contig with the focal base reserved as a single site.

    The same contig is used for both modes so that the neutral background is
    identical between selected and neutral replicates.
    """
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
        description="run2 focal selected site",
        long_description=(
            "The focal base is reserved: the selected allele is placed here by "
            "the run2 standing-variation patch, and neutral mutations are "
            "masked out of this base in every overlay."
        ),
    )
    return contig


def build_extended_events(arm: Run2Arm, mode: str) -> tuple[Any, ...]:
    """Return the stdpopsim extended events for one arm and mode.

    The neutral mode has none at all: no mutation is drawn, no fitness event is
    registered and nothing is conditioned on, so the focal base may well be
    monomorphic.  The selected mode has a draw and one fitness event per
    population that carries the allele -- and, deliberately, no
    ``ConditionOnAlleleFrequency`` events anywhere.
    """
    _require_stdpopsim()
    if mode == "neutral":
        return ()
    if mode != "selected":
        raise ValueError(f"mode must be 'selected' or 'neutral', got {mode!r}")

    tick_schedule(arm)  # assert the tick invariants before building events
    standing_time = stdpopsim.GenerationAfter(PULSE_GENERATIONS)

    events: list[Any] = [
        stdpopsim.DrawMutation(
            time=standing_time,
            single_site_id=FOCAL_SITE_ID,
            population=arm.origin_population,
        )
    ]
    if arm.introgression:
        events.append(
            stdpopsim.ChangeMutationFitness(
                start_time=standing_time,
                end_time=HAN_SPLIT_GENERATIONS,
                single_site_id=FOCAL_SITE_ID,
                population=arm.origin_population,
                selection_coeff=SELECTION_COEFFICIENT,
                dominance_coeff=DOMINANCE_COEFFICIENT,
            )
        )
        events.append(
            stdpopsim.ChangeMutationFitness(
                start_time=HAN_SPLIT_GENERATIONS,
                end_time=0.0,
                single_site_id=FOCAL_SITE_ID,
                population=arm.target_population,
                selection_coeff=SELECTION_COEFFICIENT,
                dominance_coeff=DOMINANCE_COEFFICIENT,
            )
        )
    else:
        events.append(
            stdpopsim.ChangeMutationFitness(
                start_time=standing_time,
                end_time=0.0,
                single_site_id=FOCAL_SITE_ID,
                population=arm.target_population,
                selection_coeff=SELECTION_COEFFICIENT,
                dominance_coeff=DOMINANCE_COEFFICIENT,
            )
        )
    return tuple(events)


def population_id(model: stdpopsim.DemographicModel, name: str) -> int:
    for index, population in enumerate(model.model.populations):
        if population.name == name:
            return index
    raise ValueError(f"population {name!r} is absent from the model")


# ---------------------------------------------------------------------------
# SLiM patch: standing variation
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


def _standing_variation_function(arm: Run2Arm) -> str:
    """Return the SLiM ``add_mut`` replacement for one arm."""

    if arm.introgression:
        selection = f"""    carriers = inds[inds.migrant];
    n_carriers = size(carriers);
    if (n_carriers == 0)
        err("run2: no introgressing individuals at the standing-variation tick");
    carrier_genomes = carriers.genomes;"""
        placement = "all_introgressing_genomes"
    else:
        selection = f"""    all_genomes = pop.genomes;
    n_target = asInteger(round({STANDING_FREQUENCY:.10g} * size(all_genomes)));
    if (n_target < 1)
        err("run2: standing-variation target is fewer than one genome");
    carrier_genomes = sample(all_genomes, n_target);
    n_carriers = size(carrier_genomes);"""
        placement = "uniform_random_genomes"

    return f"""// run2: place the focal allele as standing variation at the pulse tick.
function (void)add_mut(object$ mut_type, object$ pop, integer$ pos) {{
    inds = pop.individuals;
    n_individuals = size(inds);
    if (n_individuals == 0)
        err("run2: the origin population is empty at the standing-variation tick");
{selection}
    carrier_genomes.addNewDrawnMutation(mut_type, pos);
    muts = sim.mutationsOfType(mut_type);
    if (size(muts) != 1)
        err("run2: standing variation did not create exactly one focal mutation");
    frequency = sim.mutationFrequencies(pop, muts[0]);
    metadata.setValue("run2_focal_mutation_type_id", mut_type.id);
    metadata.setValue("run2_standing_tick", community.tick);
    metadata.setValue("run2_standing_population_id", pop.id);
    metadata.setValue("run2_standing_placement", "{placement}");
    metadata.setValue("run2_standing_carrier_genomes", size(carrier_genomes));
    metadata.setValue("run2_standing_population_genomes", size(pop.genomes));
    metadata.setValue("run2_standing_population_individuals", n_individuals);
    metadata.setValue("run2_standing_frequency", frequency);
}}"""


def _census_end_function(target_population_id: int) -> str:
    """Return an ``end()`` that records the true present-day census frequency.

    The focal mutation may be segregating, lost, or fixed.  Lost and fixed are
    both legitimate outcomes -- run2 conditions on nothing -- so this counts all
    three cases rather than treating any of them as an error.  Single-site
    mutation types are created with ``convertToSubstitution = F``, but the
    substitution branch is kept so a fixed allele can never be miscounted as a
    lost one.
    """
    return f"""// run2: record the present-day census frequency, then finish.
function (void)end(void) {{
    mt_id = metadata.getValue("run2_focal_mutation_type_id");
    mt_matches = sim.mutationTypes[sim.mutationTypes.id == mt_id];
    if (size(mt_matches) != 1)
        err("run2: the focal mutation type is not uniquely present at the end");
    mt = mt_matches[0];
    pop_matches = sim.subpopulations[sim.subpopulations.id == {target_population_id}];
    if (size(pop_matches) != 1)
        err("run2: the target population is not uniquely present at the end");
    pop = pop_matches[0];
    total_count = size(pop.genomes);
    if (total_count <= 0)
        err("run2: the target population is empty at the end");
    muts = sim.mutationsOfType(mt);
    subs = sim.substitutions[sim.substitutions.mutationType == mt];
    if (size(muts) + size(subs) > 1)
        err("run2: more than one focal allele is present at the end");
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
    metadata.setValue("run2_final_census_alt_count", alt_count);
    metadata.setValue("run2_final_census_total_count", total_count);
    metadata.setValue("run2_final_census_af", alt_count / total_count);
    metadata.setValue("run2_final_census_state", state);
    metadata.setValue("run2_final_census_population_id", {target_population_id});
    sim.treeSeqOutput(trees_file, metadata=metadata);
    sim.simulationFinished();
}}"""


_ENTRY_FUNCTION = """
// run2: record the focal frequency at the tick the target population is founded.
function (void)run2_record_entry(integer$ target_id) {
    mt_id = metadata.getValue("run2_focal_mutation_type_id");
    mt_matches = sim.mutationTypes[sim.mutationTypes.id == mt_id];
    if (size(mt_matches) != 1)
        err("run2: the focal mutation type is not unique at the entry tick");
    pop_matches = sim.subpopulations[sim.subpopulations.id == target_id];
    if (size(pop_matches) != 1)
        err("run2: the target population does not exist at the entry tick");
    pop = pop_matches[0];
    metadata.setValue("run2_entry_tick", community.tick);
    metadata.setValue("run2_entry_population_id", target_id);
    metadata.setValue("run2_entry_total_genomes", size(pop.genomes));
    metadata.setValue("run2_entry_af", af(mt_matches[0], pop));
}
"""

_END_REGISTRATION = (
    '    community.registerLateEvent(NULL, "{dbg(self.source); end();}", G_end, G_end);'
)


def _entry_registration(target_population_id: int, entry_generations: float) -> str:
    """Register the entry recorder one step before stdpopsim's own end hook.

    Registration order is execution order for same-tick late events, and the
    split that founds the target population is registered earlier in the block,
    so this observes the population as founded.
    """
    years = float(entry_generations) * GENERATION_TIME_YEARS
    return (
        "    // run2: record the focal frequency as the target population is founded.\n"
        f'    community.registerLateEvent(NULL, "{{run2_record_entry({target_population_id});}}",\n'
        f"        time_to_tick({years:.1f}), time_to_tick({years:.1f}));\n\n"
        f"{_END_REGISTRATION}"
    )


def _replace_once(text: str, pattern: re.Pattern[str], replacement: str, label: str) -> str:
    matches = pattern.findall(text)
    if len(matches) != 1:
        raise RuntimeError(
            f"stdpopsim {label} contract changed: expected exactly one match, "
            f"found {len(matches)}"
        )
    return pattern.sub(lambda _: replacement, text, count=1)


@contextmanager
def scoped_slim_patch(
    arm: Run2Arm, mode: str, target_population_id: int
) -> Iterator[dict[str, Any]]:
    """Patch stdpopsim's SLiM helper functions for the duration of one call.

    In ``selected`` mode this swaps ``add_mut`` for the standing-variation
    placement and ``end`` for the census recorder.  In ``neutral`` mode nothing
    is patched -- there is no focal mutation to place or count.
    """
    original_object = slim_engine._slim_functions
    original = str(original_object)
    record: dict[str, Any] = {
        "mode": mode,
        "arm_id": arm.arm_id,
        "original_sha256": hashlib.sha256(original.encode()).hexdigest(),
        "add_mut_replaced": False,
        "end_replaced": False,
    }
    if mode == "neutral":
        record["patched_sha256"] = record["original_sha256"]
        yield record
        return

    add_mut = _standing_variation_function(arm)
    end = _census_end_function(target_population_id)
    patched = _replace_once(original, _ADD_MUT_PATTERN, add_mut, "add_mut")
    patched = _replace_once(patched, _END_PATTERN, end, "end")
    record.update(
        {
            "add_mut_replaced": True,
            "end_replaced": True,
            "add_mut_sha256": hashlib.sha256(add_mut.encode()).hexdigest(),
            "end_sha256": hashlib.sha256(end.encode()).hexdigest(),
            "target_population_id": target_population_id,
        }
    )

    # Only the introgression arm has a founding event to observe: the EAS arm's
    # entry frequency is the placement itself, already recorded by add_mut.
    original_main_object = slim_engine._slim_main
    original_main = str(original_main_object)
    patched_main = original_main
    if arm.introgression:
        entry_generations = tick_schedule(arm)["realized_generations_ago"]["han_split"]
        patched += _ENTRY_FUNCTION
        if original_main.count(_END_REGISTRATION) != 1:
            raise RuntimeError(
                "stdpopsim end-hook registration contract changed; cannot place "
                "the run2 entry recorder"
            )
        patched_main = original_main.replace(
            _END_REGISTRATION,
            _entry_registration(target_population_id, entry_generations),
            1,
        )
        record.update(
            {
                "entry_recorder_installed": True,
                "entry_generations_ago": entry_generations,
                "patched_main_sha256": hashlib.sha256(patched_main.encode()).hexdigest(),
            }
        )
    else:
        record["entry_recorder_installed"] = False

    record["patched_sha256"] = hashlib.sha256(patched.encode()).hexdigest()
    slim_engine._slim_functions = patched
    slim_engine._slim_main = patched_main
    try:
        yield record
    finally:
        slim_engine._slim_functions = original_object
        slim_engine._slim_main = original_main_object


# ---------------------------------------------------------------------------
# msprime patch: reserve the focal base from the neutral overlay
# ---------------------------------------------------------------------------


def mask_focal_rate_map(
    rate_map: msprime.RateMap, focal_position: int = FOCAL_POSITION_BP
) -> tuple[msprime.RateMap, bool]:
    """Zero the mutation rate on ``[focal_position, focal_position + 1)``.

    Rates on every other interval, and therefore every other integrated rate,
    are preserved exactly.  Returns the input unchanged when the focal base
    already carries no mutational mass.
    """
    positions = np.asarray(rate_map.position, dtype=float)
    rates = np.asarray(rate_map.rate, dtype=float)
    low = float(focal_position)
    high = low + 1.0
    if low < positions[0] or high > positions[-1]:
        return rate_map, False

    refined = np.unique(np.concatenate([positions, [low, high]]))
    midpoints = 0.5 * (refined[:-1] + refined[1:])
    source = np.searchsorted(positions, midpoints, side="right") - 1
    source = np.clip(source, 0, len(rates) - 1)
    refined_rates = rates[source].copy()

    target = (refined[:-1] >= low) & (refined[1:] <= high)
    if not np.any(refined_rates[target] > 0):
        return rate_map, False
    refined_rates[target] = 0.0
    return msprime.RateMap(position=refined, rate=refined_rates), True


@contextmanager
def scoped_focal_overlay_patch(
    focal_position: int = FOCAL_POSITION_BP,
) -> Iterator[dict[str, Any]]:
    """Mask focal mutational mass in every neutral overlay of one engine call.

    stdpopsim already excludes a single-site DFE referenced by an extended event
    from the background overlay, but the recapitation DFE spans the whole contig
    and would otherwise drop neutral mutations onto the reserved base.  This
    guard masks *every* intercepted rate map that carries positive rate there,
    so it behaves identically in selected and neutral mode.
    """
    original = msprime.sim_mutations
    signature = inspect.signature(original)
    receipt: dict[str, Any] = {
        "intercepted_calls": 0,
        "masked_calls": 0,
        "focal_position_0based": int(focal_position),
        "target_interval_half_open": [int(focal_position), int(focal_position) + 1],
        "restored": False,
    }

    def guarded(*args: Any, **kwargs: Any):
        bound = signature.bind(*args, **kwargs)
        bound.apply_defaults()
        receipt["intercepted_calls"] += 1
        rate = bound.arguments.get("rate")
        if isinstance(rate, msprime.RateMap):
            replacement, masked = mask_focal_rate_map(rate, focal_position)
            if masked:
                bound.arguments["rate"] = replacement
                receipt["masked_calls"] += 1
        return original(*bound.args, **bound.kwargs)

    msprime.sim_mutations = guarded
    try:
        yield receipt
    finally:
        msprime.sim_mutations = original
        receipt["restored"] = True


# ---------------------------------------------------------------------------
# Script generation (validation without a SLiM binary)
# ---------------------------------------------------------------------------


def generate_slim_script(arm: Run2Arm, mode: str, repo_root: str | Path) -> str:
    """Return the SLiM script stdpopsim would run, with run2's patches applied.

    This needs stdpopsim but not the SLiM binary, so the whole event schedule
    can be inspected and asserted before any simulation is launched.
    """
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
                extended_events=list(events) or None,
                slim_script=True,
                slim_scaling_factor=SLIM_SCALING_FACTOR,
                slim_burn_in=SLIM_BURN_IN,
            )
    return buffer.getvalue()


def model_record(arm: Run2Arm, mode: str, repo_root: str | Path) -> dict[str, Any]:
    """Return a JSON-serializable description of one arm/mode simulation law."""

    model = build_demographic_model(arm, repo_root)
    events = build_extended_events(arm, mode)
    return {
        "schema": "gamma-smc.run2-model/v1",
        "arm_id": arm.arm_id,
        "mode": mode,
        "demographic_model_id": model.id,
        "target_population": arm.target_population,
        "target_population_id": population_id(model, arm.target_population),
        "origin_population": arm.origin_population,
        "standing_variation_placement": arm.placement if mode == "selected" else None,
        "standing_frequency": STANDING_FREQUENCY if mode == "selected" else None,
        "selection_coefficient": SELECTION_COEFFICIENT if mode == "selected" else 0.0,
        "dominance_coefficient": DOMINANCE_COEFFICIENT if mode == "selected" else None,
        "conditioning_events": 0,
        "tick_schedule": tick_schedule(arm),
        "extended_events": serialize_extended_events(events),
        "stdpopsim_version": stdpopsim.__version__,
    }


__all__ = [
    "build_contig",
    "build_demographic_model",
    "build_extended_events",
    "epoch_tick_offset",
    "event_tick_offset",
    "generate_slim_script",
    "mask_focal_rate_map",
    "model_record",
    "phlash_artifact_path",
    "population_id",
    "realized_generations",
    "scoped_focal_overlay_patch",
    "scoped_slim_patch",
    "snap_generations",
    "tick_schedule",
]
