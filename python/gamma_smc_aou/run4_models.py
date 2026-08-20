"""Demography, contig, events and SLiM patches for the run4 study.

Three placements, one per arm:

* **archaic single founder** -- the allele is fixed in Neanderthal one tick after
  the archaic split, introgresses into CHB by continuous migration, and is then
  reduced at the selection onset to a *single* archaic founding genome.  Every
  present-day carrier therefore descends from one introgressed haplotype, which
  is what EPAS1 is.  If no CHB genome carries archaic ancestry at the focal base
  at that tick, the event calls stdpopsim's own ``restore()`` and the replicate
  is re-drawn from the checkpoint -- the ascertainment that makes such a locus
  worth studying at all.
* **de novo** -- a single copy on one random genome of the target population at
  the selection onset.  No archaic background.

Both selected placements are hard sweeps from one founder, so both need
conditioning on the allele surviving: a single beneficial copy is lost with
probability roughly ``1 - 2hs``.  That conditioning is expressed as a stdpopsim
``ConditionOnAlleleFrequency``, which SLiM implements with its own
save/restore, so for the de novo arms a rejection replays only the ~41 ticks
since the draw rather than the whole simulation.

The neutral mode has **no extended events at all**: no allele is drawn, nothing
is conditioned, and the focal base may be empty.
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
from .run4_config import (
    MARKER_POSITION_BP,
    MARKER_SITE_ID,
    RESERVED_POSITIONS,
    ARCHAIC_POPULATION,
    CHB_ARCHAIC_MIGRATION_END_GENERATIONS,
    CHB_SPLIT_GENERATIONS,
    DOMINANCE_COEFFICIENT,
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
    Run4Arm,
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
    "epoch_tick_offset",
    "event_tick_offset",
    "generate_slim_script",
    "mask_reserved_rate_map",
    "model_record",
    "population_id",
    "reserved_positions",
    "scoped_reserved_overlay_patch",
    "scoped_slim_patch",
    "snap_generations",
    "tick_schedule",
]


def snap_generations(time: Any, scaling_factor: float = SLIM_SCALING_FACTOR) -> float:
    snapped = round(float(time) / scaling_factor) * scaling_factor
    if isinstance(time, stdpopsim.GenerationAfter):
        snapped -= scaling_factor
    if snapped < 0:
        raise ValueError(f"snapped time {snapped} is negative")
    return float(snapped)


def event_tick_offset(
    time: Any, generation_time: float, *, scaling_factor: float = SLIM_SCALING_FACTOR
) -> int:
    years = snap_generations(time, scaling_factor) * generation_time
    return int(years / generation_time / scaling_factor)


def epoch_tick_offset(
    generations_ago: float,
    generation_time: float,
    *,
    scaling_factor: float = SLIM_SCALING_FACTOR,
) -> int:
    years = round(float(generations_ago) * generation_time)
    return int(years / generation_time / scaling_factor)


def tick_schedule(arm: Run4Arm) -> dict[str, Any]:
    q = SLIM_SCALING_FACTOR
    gt = arm.generation_time
    onset = arm.onset_generations()
    onset_offset = event_tick_offset(onset, gt)
    schedule: dict[str, Any] = {
        "scaling_factor": q,
        "generation_time_years": gt,
        "tick_offsets_before_final_tick": {"present": 0, "selection_onset": onset_offset},
        "realized_generations_ago": {"present": 0.0, "selection_onset": onset_offset * q},
        "selection_intervals_generations_ago": [
            {"population": arm.target_population, "start": onset, "end": 0.0}
        ],
    }
    if arm.archaic:
        split_offset = epoch_tick_offset(NEANDERTHAL_SPLIT_GENERATIONS, gt)
        fixation_offset = event_tick_offset(
            stdpopsim.GenerationAfter(NEANDERTHAL_SPLIT_GENERATIONS), gt
        )
        if fixation_offset != split_offset - 1:
            raise RuntimeError("archaic fixation must be one tick after the split")
        migration_end = epoch_tick_offset(CHB_ARCHAIC_MIGRATION_END_GENERATIONS, gt)
        chb_split = epoch_tick_offset(CHB_SPLIT_GENERATIONS, gt)
        if not split_offset > chb_split > migration_end > onset_offset > 0:
            raise RuntimeError(
                "the archaic schedule must run split > CHB founding > migration "
                "end > selection onset > present"
            )
        schedule["tick_offsets_before_final_tick"].update(
            {
                "archaic_split": split_offset,
                "archaic_fixation": fixation_offset,
                "chb_founding": chb_split,
                "archaic_migration_end": migration_end,
            }
        )
        schedule["realized_generations_ago"].update(
            {
                "archaic_split": split_offset * q,
                "archaic_fixation": fixation_offset * q,
                "chb_founding": chb_split * q,
                "archaic_migration_end": migration_end * q,
            }
        )
    schedule["realized_years_ago"] = {
        k: v * gt for k, v in schedule["realized_generations_ago"].items()
    }
    return schedule


# ---------------------------------------------------------------------------
# Demography and contig
# ---------------------------------------------------------------------------


def _require_stdpopsim() -> None:
    if stdpopsim.__version__ != REQUIRED_STDPOPSIM_VERSION:
        raise RuntimeError(
            f"run4 requires stdpopsim {REQUIRED_STDPOPSIM_VERSION}; "
            f"found {stdpopsim.__version__}"
        )


def _load_ooa_archaic_model() -> stdpopsim.DemographicModel:
    _require_stdpopsim()
    model = stdpopsim.get_species("HomSap").get_demographic_model(OOA_MODEL_ID)
    if not np.isclose(model.generation_time, OOA_GENERATION_TIME):
        raise RuntimeError(f"{OOA_MODEL_ID} generation time changed")
    names = {p.name for p in model.model.populations}
    for required in (ARCHAIC_POPULATION, "CHB"):
        if required not in names:
            raise RuntimeError(f"{OOA_MODEL_ID} is missing population {required!r}")
    model.mutation_rate = MUTATION_RATE
    model.recombination_rate = RECOMBINATION_RATE
    return model


def build_demographic_model(arm: Run4Arm, repo_root: str | Path):
    _require_stdpopsim()
    if arm.demography_id == OOA_MODEL_ID:
        return _load_ooa_archaic_model()
    artifact = load_phlash_eas_npz(
        Path(repo_root) / PHLASH_EAS_RELATIVE_PATH, expected_sha256=PHLASH_EAS_SHA256
    )
    return build_eas_demography_models(artifact)["median"].stdpopsim_model


def build_contig(arm: Run4Arm | None = None) -> stdpopsim.Contig:
    """Return the contig. The archaic arm reserves a second base for the marker."""
    _require_stdpopsim()
    contig = stdpopsim.get_species("HomSap").get_contig(
        length=SEQUENCE_LENGTH_BP,
        mutation_rate=MUTATION_RATE,
        recombination_rate=RECOMBINATION_RATE,
    )
    contig.add_single_site(
        id=FOCAL_SITE_ID,
        coordinate=FOCAL_POSITION_BP,
        description="run4 focal selected site",
    )
    return contig


def reserved_positions(arm: Run4Arm) -> tuple[int, ...]:
    return (FOCAL_POSITION_BP,)


def mask_reserved_rate_map(rate_map, positions):
    """Zero the mutation rate on each reserved base, preserving all other rates."""
    import msprime

    result, changed_any = rate_map, False
    for position in positions:
        pos = np.asarray(result.position, dtype=float)
        rates = np.asarray(result.rate, dtype=float)
        low, high = float(position), float(position) + 1.0
        if low < pos[0] or high > pos[-1]:
            continue
        refined = np.unique(np.concatenate([pos, [low, high]]))
        mid = 0.5 * (refined[:-1] + refined[1:])
        source = np.clip(np.searchsorted(pos, mid, side="right") - 1, 0, len(rates) - 1)
        refined_rates = rates[source].copy()
        target = (refined[:-1] >= low) & (refined[1:] <= high)
        if not np.any(refined_rates[target] > 0):
            continue
        refined_rates[target] = 0.0
        result = msprime.RateMap(position=refined, rate=refined_rates)
        changed_any = True
    return result, changed_any


@contextmanager
def scoped_reserved_overlay_patch(positions) -> Iterator[dict[str, Any]]:
    """Mask every reserved base in each neutral overlay of one engine call."""
    import inspect

    import msprime

    original = msprime.sim_mutations
    signature = inspect.signature(original)
    receipt: dict[str, Any] = {
        "intercepted_calls": 0,
        "masked_calls": 0,
        "reserved_positions": [int(p) for p in positions],
        "restored": False,
    }

    def guarded(*args: Any, **kwargs: Any):
        bound = signature.bind(*args, **kwargs)
        bound.apply_defaults()
        receipt["intercepted_calls"] += 1
        rate = bound.arguments.get("rate")
        if isinstance(rate, msprime.RateMap):
            replacement, masked = mask_reserved_rate_map(rate, positions)
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


def population_id(model, name: str) -> int:
    for index, population in enumerate(model.model.populations):
        if population.name == name:
            return index
    raise ValueError(f"population {name!r} is absent from the model")


def build_extended_events(arm: Run4Arm, mode: str) -> tuple[Any, ...]:
    """Return the extended events. The neutral mode has none at all."""
    _require_stdpopsim()
    if mode == "neutral":
        return ()
    if mode != "selected":
        raise ValueError(f"mode must be 'selected' or 'neutral', got {mode!r}")
    tick_schedule(arm)

    onset = arm.onset_generations()
    events: list[Any] = []

    if arm.archaic:
        # The allele is fixed in Neanderthal and introgresses, so at the onset it
        # is the whole archaic haplotype class at ~1.5%. That head start is what
        # lets EPAS1 reach 63-87% in 6,000 years; a single copy cannot.
        events.append(
            stdpopsim.DrawMutation(
                time=stdpopsim.GenerationAfter(NEANDERTHAL_SPLIT_GENERATIONS),
                single_site_id=FOCAL_SITE_ID,
                population=ARCHAIC_POPULATION,
            )
        )
        # Ascertainment: archaic ancestry must actually be present at the focal
        # base when selection starts. SLiM handles the re-draw with restore().
        events.append(
            stdpopsim.ConditionOnAlleleFrequency(
                start_time=onset,
                end_time=onset,
                single_site_id=FOCAL_SITE_ID,
                population=arm.target_population,
                op=">",
                allele_frequency=0.0,
            )
        )
    else:
        events.append(
            stdpopsim.DrawMutation(
                time=onset,
                single_site_id=FOCAL_SITE_ID,
                population=arm.target_population,
            )
        )

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
    events.append(
        stdpopsim.ConditionOnAlleleFrequency(
            start_time=stdpopsim.GenerationAfter(onset),
            end_time=0.0,
            single_site_id=FOCAL_SITE_ID,
            population=arm.target_population,
            op=">",
            allele_frequency=0.0,
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
_END_REGISTRATION = (
    '    community.registerLateEvent(NULL, "{dbg(self.source); end();}", G_end, G_end);'
)


def _placement_function(arm: Run4Arm) -> str:
    """Return the ``add_mut`` replacement for one arm.

    The archaic arm fixes the allele across the whole Neanderthal population, so
    it marks archaic ancestry exactly and arrives in CHB as a haplotype class at
    roughly the archaic ancestry proportion.  The de novo arms place a single
    copy.  Nothing is ever removed: SLiM's ``removeMutations`` leaves empty
    derived states that stdpopsim's recapitation cannot parse.
    """
    if arm.archaic:
        body = """    carrier_genomes = pop.genomes;
    if (size(carrier_genomes) == 0)
        err("run4: the archaic population is empty at the fixation tick");"""
        label = "fixed_in_archaic_population"
    else:
        body = """    carrier_genomes = sample(pop.genomes, 1);"""
        label = "single_de_novo_copy"

    return f"""// run4: introduce the focal allele.
function (void)add_mut(object$ mut_type, object$ pop, integer$ pos) {{
    if (size(pop.individuals) == 0)
        err("run4: the origin population is empty at the placement tick");
{body}
    carrier_genomes.addNewDrawnMutation(mut_type, pos);
    muts = sim.mutationsOfType(mut_type);
    if (size(muts) != 1)
        err("run4: placement did not create exactly one focal mutation");
    metadata.setValue("run4_focal_mutation_type_id", mut_type.id);
    metadata.setValue("run4_focal_position", pos);
    metadata.setValue("run4_placement_tick", community.tick);
    metadata.setValue("run4_placement_mode", "{label}");
    metadata.setValue("run4_placement_carrier_genomes", size(carrier_genomes));
    metadata.setValue("run4_placement_frequency", sim.mutationFrequencies(pop, muts[0]));
}}"""


def _census_end_function(target_population_id: int) -> str:
    return f"""// run4: record the present-day census frequency, then finish.
function (void)end(void) {{
    mt_id = metadata.getValue("run4_focal_mutation_type_id");
    total_count = 0;
    alt_count = 0;
    state = "absent";
    pop_matches = sim.subpopulations[sim.subpopulations.id == {target_population_id}];
    if (size(pop_matches) != 1)
        err("run4: the target population is not unique at the end");
    pop = pop_matches[0];
    total_count = size(pop.genomes);
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
    metadata.setValue("run4_final_census_alt_count", alt_count);
    metadata.setValue("run4_final_census_total_count", total_count);
    metadata.setValue("run4_final_census_af", alt_count / total_count);
    metadata.setValue("run4_final_census_state", state);
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


def _reduction_registration(target_population_id: int, onset_years: float) -> str:
    return (
        "    // run4: reduce to a single archaic founder at the selection onset.\n"
        f'    community.registerLateEvent(NULL, "{{run4_single_founder({target_population_id});}}",\n'
        f"        time_to_tick({onset_years:.1f}), time_to_tick({onset_years:.1f}));\n\n"
        f"{_END_REGISTRATION}"
    )


@contextmanager
def scoped_slim_patch(
    arm: Run4Arm, mode: str, target_population_id: int
) -> Iterator[dict[str, Any]]:
    """Patch stdpopsim's SLiM helpers for one engine call."""
    original_functions = slim_engine._slim_functions
    original_main = slim_engine._slim_main
    functions = str(original_functions)
    main = str(original_main)
    record: dict[str, Any] = {"arm_id": arm.arm_id, "mode": mode}

    if mode == "neutral":
        # Nothing is drawn, so nothing needs patching; the stock end() is used.
        record["patched"] = False
        yield record
        return

    patched_functions = _replace_once(
        functions, _ADD_MUT_PATTERN, _placement_function(arm), "add_mut"
    )
    patched_functions = _replace_once(
        patched_functions, _END_PATTERN, _census_end_function(target_population_id), "end"
    )
    # No main-block patch is needed any more: the single-founder seeding happens
    # inside add_mut, which stdpopsim already schedules at the draw tick.
    patched_main = main
    record["single_founder_seeding"] = bool(arm.archaic)

    record.update(
        {
            "patched": True,
            "functions_sha256": hashlib.sha256(patched_functions.encode()).hexdigest(),
            "main_sha256": hashlib.sha256(patched_main.encode()).hexdigest(),
            "target_population_id": target_population_id,
        }
    )
    slim_engine._slim_functions = patched_functions
    slim_engine._slim_main = patched_main
    try:
        yield record
    finally:
        slim_engine._slim_functions = original_functions
        slim_engine._slim_main = original_main


def generate_slim_script(arm: Run4Arm, mode: str, repo_root: str | Path) -> str:
    model = build_demographic_model(arm, repo_root)
    contig = build_contig(arm)
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


def model_record(arm: Run4Arm, mode: str, repo_root: str | Path) -> dict[str, Any]:
    model = build_demographic_model(arm, repo_root)
    events = build_extended_events(arm, mode)
    return {
        "schema": "gamma-smc.run4-model/v1",
        "arm_id": arm.arm_id,
        "mode": mode,
        "demographic_model_id": model.id,
        "generation_time_years": arm.generation_time,
        "target_population": arm.target_population,
        "target_population_id": population_id(model, arm.target_population),
        "origin": arm.origin,
        "selection_coefficient": SELECTION_COEFFICIENT if mode == "selected" else 0.0,
        "selection_onset_generations": arm.onset_generations(),
        "conditioned_on_survival": mode == "selected",
        "tick_schedule": tick_schedule(arm),
        "extended_events": serialize_extended_events(events),
        "stdpopsim_version": stdpopsim.__version__,
    }
