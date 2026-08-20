"""Demography, contig, events and engine dispatch for the run5 study.

Two engines:

* **SLiM**, for the two CHB arms, where archaic ancestry has to be simulated
  forwards so that introgressed material is genuinely archaic in descent.
* **msprime**, for the EAS arm, which is an idealized structured-coalescent
  sweep used purely as a calibration reference.  It runs in under a second per
  replicate.

msprime rejects demographic events *during* a sweep, and the PHLASH curve has
ten thousand of them.  Rather than flattening the whole history, the EAS arm
keeps it as **discrete constant-size epochs** everywhere older than the sweep
window and holds the size flat only across the window itself, at the harmonic
mean.  The only constraint msprime imposes is that no epoch boundary falls
inside the sweep, so this satisfies it while preserving the deep EAS history --
which matters, because that history sets the neutral TMRCA background and hence
the null distribution.
"""

from __future__ import annotations

import hashlib
import io
import re
from contextlib import contextmanager, redirect_stdout
from pathlib import Path
from typing import Any, Iterator

import msprime
import numpy as np
import stdpopsim
from stdpopsim import slim_engine

from .eas_sweep_models import serialize_extended_events
from .run5_config import (
    ARCHAIC_POPULATION,
    CHB_ARCHAIC_MIGRATION_END_GENERATIONS,
    CHB_SPLIT_GENERATIONS,
    DOMINANCE_COEFFICIENT,
    EAS_END_FREQUENCY,
    EAS_START_FREQUENCY,
    EAS_EPOCH_COUNT,
    EAS_FLAT_WINDOW_GENERATIONS,
    EAS_OLDEST_EPOCH_GENERATIONS,
    EAS_SWEEP_DT,
    FOCAL_POSITION_BP,
    FOCAL_SITE_ID,
    MIN_INTROGRESSED_FREQUENCY,
    MUTATION_RATE,
    NEANDERTHAL_SPLIT_GENERATIONS,
    NEA_STARTING_FREQUENCY,
    OOA_GENERATION_TIME,
    OOA_MODEL_ID,
    PHLASH_EAS_RELATIVE_PATH,
    RECOMBINATION_RATE,
    REQUIRED_STDPOPSIM_VERSION,
    Run5Arm,
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
    "build_eas_demography",
    "eas_effective_size",
    "epoch_tick_offset",
    "event_tick_offset",
    "generate_slim_script",
    "model_record",
    "population_id",
    "simulate_msprime",
    "scoped_slim_patch",
    "snap_generations",
    "tick_schedule",
]


# ---------------------------------------------------------------------------
# Tick arithmetic
# ---------------------------------------------------------------------------


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


def tick_schedule(arm: Run5Arm) -> dict[str, Any]:
    q = SLIM_SCALING_FACTOR
    gt = arm.generation_time
    if arm.engine == "msprime":
        return {
            "engine": "msprime",
            "generation_time_years": gt,
            "sweep_start_frequency": EAS_START_FREQUENCY,
            "sweep_end_frequency": EAS_END_FREQUENCY,
            "selection_coefficient": SELECTION_COEFFICIENT,
        }

    onset = arm.onset_generations()
    onset_offset = event_tick_offset(onset, gt)
    split_offset = epoch_tick_offset(NEANDERTHAL_SPLIT_GENERATIONS, gt)
    placement_offset = event_tick_offset(
        stdpopsim.GenerationAfter(NEANDERTHAL_SPLIT_GENERATIONS), gt
    )
    migration_end = epoch_tick_offset(CHB_ARCHAIC_MIGRATION_END_GENERATIONS, gt)
    chb_split = epoch_tick_offset(CHB_SPLIT_GENERATIONS, gt)

    if placement_offset != split_offset - 1:
        raise RuntimeError("archaic placement must be one tick after the split")
    if onset_offset != chb_split:
        raise RuntimeError(
            "selection must start on the tick CHB is founded and archaic "
            f"migration begins: founding {chb_split}, onset {onset_offset}"
        )
    if not split_offset > chb_split > migration_end > 0:
        raise RuntimeError("the archaic event schedule is out of order")

    schedule = {
        "engine": "slim",
        "scaling_factor": q,
        "generation_time_years": gt,
        "tick_offsets_before_final_tick": {
            "archaic_split": split_offset,
            "archaic_placement": placement_offset,
            "chb_founding_and_selection_onset": chb_split,
            "archaic_migration_end": migration_end,
            "present": 0,
        },
        "realized_generations_ago": {
            "archaic_split": split_offset * q,
            "archaic_placement": placement_offset * q,
            "chb_founding_and_selection_onset": chb_split * q,
            "archaic_migration_end": migration_end * q,
            "present": 0.0,
        },
        "selection_intervals_generations_ago": (
            [{"population": ARCHAIC_POPULATION, "start": placement_offset * q, "end": 0.0}]
            if arm.nea_selection
            else []
        )
        + [{"population": arm.target_population, "start": chb_split * q, "end": 0.0}],
        "ascertainment": {
            "population": arm.target_population,
            "generations_ago": migration_end * q,
            "minimum_frequency": MIN_INTROGRESSED_FREQUENCY,
        },
    }
    schedule["realized_years_ago"] = {
        k: v * gt for k, v in schedule["realized_generations_ago"].items()
    }
    return schedule


# ---------------------------------------------------------------------------
# Demography
# ---------------------------------------------------------------------------


def _require_stdpopsim() -> None:
    if stdpopsim.__version__ != REQUIRED_STDPOPSIM_VERSION:
        raise RuntimeError(
            f"run5 requires stdpopsim {REQUIRED_STDPOPSIM_VERSION}; "
            f"found {stdpopsim.__version__}"
        )


def _phlash_median(repo_root: str | Path):
    path = Path(repo_root) / PHLASH_EAS_RELATIVE_PATH
    with np.load(path, allow_pickle=False) as data:
        times = np.asarray(data["time"], dtype=float)
        median = np.median(np.asarray(data["bootstrap_ne"], dtype=float), axis=0)
    return times, median


def eas_effective_size(
    repo_root: str | Path, generations: float = EAS_FLAT_WINDOW_GENERATIONS
) -> float:
    """Harmonic-mean PHLASH median Ne across the sweep window.

    The harmonic mean is the size that governs coalescence over a period of
    changing size, so holding the window flat at this value distorts the sweep
    phase as little as a single number can.
    """
    times, median = _phlash_median(repo_root)

    def size_at(generation: float) -> float:
        if generation <= times[0]:
            return float(median[0])
        return float(median[int(np.argmin(np.abs(times - generation)))])

    grid = np.arange(1, int(generations) + 1)
    return float(len(grid) / np.sum([1.0 / size_at(g) for g in grid]))


def build_eas_demography(repo_root: str | Path) -> msprime.Demography:
    """PHLASH EAS as discrete epochs, held flat across the sweep window.

    msprime rejects demographic events during a sweep, so the only constraint is
    that no epoch boundary falls inside it. Everything older keeps the PHLASH
    shape, which is what the neutral TMRCA background -- and hence the null --
    depends on.
    """
    times, median = _phlash_median(repo_root)

    def size_at(generation: float) -> float:
        if generation <= times[0]:
            return float(median[0])
        return float(median[int(np.argmin(np.abs(times - generation)))])

    demography = msprime.Demography()
    demography.add_population(
        name="EAS", initial_size=eas_effective_size(repo_root)
    )
    edges = np.unique(
        np.round(
            np.geomspace(
                EAS_FLAT_WINDOW_GENERATIONS,
                EAS_OLDEST_EPOCH_GENERATIONS,
                EAS_EPOCH_COUNT,
            )
        )
    )
    for edge in edges:
        demography.add_population_parameters_change(
            time=float(edge), initial_size=size_at(float(edge)), population="EAS"
        )
    return demography


def build_demographic_model(arm: Run5Arm, repo_root: str | Path):
    """Return a stdpopsim model (SLiM arms) or an msprime Demography (EAS)."""
    if arm.engine == "msprime":
        return build_eas_demography(repo_root)
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
        description="run5 focal selected site",
    )
    return contig


def population_id(model, name: str) -> int:
    for index, population in enumerate(model.model.populations):
        if population.name == name:
            return index
    raise ValueError(f"population {name!r} is absent from the model")


# ---------------------------------------------------------------------------
# Extended events (SLiM arms)
# ---------------------------------------------------------------------------


def build_extended_events(arm: Run5Arm, mode: str) -> tuple[Any, ...]:
    """Return the extended events. The neutral mode has none at all."""
    if arm.engine == "msprime":
        return ()
    _require_stdpopsim()
    if mode == "neutral":
        return ()
    if mode != "selected":
        raise ValueError(f"mode must be 'selected' or 'neutral', got {mode!r}")
    tick_schedule(arm)

    placement = stdpopsim.GenerationAfter(NEANDERTHAL_SPLIT_GENERATIONS)
    onset = arm.onset_generations()
    events: list[Any] = [
        stdpopsim.DrawMutation(
            time=placement,
            single_site_id=FOCAL_SITE_ID,
            population=ARCHAIC_POPULATION,
        )
    ]

    if arm.nea_selection:
        # A sweep inside the archaic source compresses archaic haplotype
        # diversity, so what introgresses descends from fewer, longer haplotypes.
        events.append(
            stdpopsim.ChangeMutationFitness(
                start_time=placement,
                end_time=0.0,
                single_site_id=FOCAL_SITE_ID,
                population=ARCHAIC_POPULATION,
                selection_coeff=SELECTION_COEFFICIENT,
                dominance_coeff=DOMINANCE_COEFFICIENT,
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
    # Ascertainment: the allele must actually have introgressed to an
    # appreciable frequency by the time archaic migration stops.
    events.append(
        stdpopsim.ConditionOnAlleleFrequency(
            start_time=snap_generations(CHB_ARCHAIC_MIGRATION_END_GENERATIONS),
            end_time=snap_generations(CHB_ARCHAIC_MIGRATION_END_GENERATIONS),
            single_site_id=FOCAL_SITE_ID,
            population=arm.target_population,
            op=">=",
            allele_frequency=MIN_INTROGRESSED_FREQUENCY,
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


def _placement_function(arm: Run5Arm) -> str:
    if arm.nea_selection:
        body = f"""    all_genomes = pop.genomes;
    n_target = asInteger(round({NEA_STARTING_FREQUENCY:.10g} * size(all_genomes)));
    if (n_target < 1)
        err("run5: the archaic standing-variation target is fewer than one genome");
    carrier_genomes = sample(all_genomes, n_target);"""
        label = "archaic_standing_variation"
    else:
        body = """    carrier_genomes = pop.genomes;
    if (size(carrier_genomes) == 0)
        err("run5: the archaic population is empty at the placement tick");"""
        label = "fixed_in_archaic_population"

    return f"""// run5: introduce the focal allele into the archaic population.
function (void)add_mut(object$ mut_type, object$ pop, integer$ pos) {{
    if (size(pop.individuals) == 0)
        err("run5: the origin population is empty at the placement tick");
{body}
    carrier_genomes.addNewDrawnMutation(mut_type, pos);
    muts = sim.mutationsOfType(mut_type);
    if (size(muts) != 1)
        err("run5: placement did not create exactly one focal mutation");
    metadata.setValue("run5_focal_mutation_type_id", mut_type.id);
    metadata.setValue("run5_placement_tick", community.tick);
    metadata.setValue("run5_placement_mode", "{label}");
    metadata.setValue("run5_placement_carrier_genomes", size(carrier_genomes));
    metadata.setValue("run5_placement_frequency", sim.mutationFrequencies(pop, muts[0]));
}}"""


def _census_end_function(target_population_id: int) -> str:
    return f"""// run5: record the present-day census frequency, then finish.
function (void)end(void) {{
    mt_id = metadata.getValue("run5_focal_mutation_type_id");
    pop_matches = sim.subpopulations[sim.subpopulations.id == {target_population_id}];
    if (size(pop_matches) != 1)
        err("run5: the target population is not unique at the end");
    pop = pop_matches[0];
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
    metadata.setValue("run5_final_census_alt_count", alt_count);
    metadata.setValue("run5_final_census_total_count", total_count);
    metadata.setValue("run5_final_census_af", alt_count / total_count);
    metadata.setValue("run5_final_census_state", state);
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
    arm: Run5Arm, mode: str, target_population_id: int
) -> Iterator[dict[str, Any]]:
    original = slim_engine._slim_functions
    text = str(original)
    record: dict[str, Any] = {"arm_id": arm.arm_id, "mode": mode}
    if mode == "neutral" or arm.engine == "msprime":
        record["patched"] = False
        yield record
        return
    patched = _replace_once(text, _ADD_MUT_PATTERN, _placement_function(arm), "add_mut")
    patched = _replace_once(
        patched, _END_PATTERN, _census_end_function(target_population_id), "end"
    )
    record.update(
        {
            "patched": True,
            "patched_sha256": hashlib.sha256(patched.encode()).hexdigest(),
            "target_population_id": target_population_id,
        }
    )
    slim_engine._slim_functions = patched
    try:
        yield record
    finally:
        slim_engine._slim_functions = original


# ---------------------------------------------------------------------------
# msprime arm
# ---------------------------------------------------------------------------


def simulate_msprime(arm: Run5Arm, mode: str, seed: int, repo_root: str | Path):
    """Simulate one EAS replicate under the structured-coalescent sweep model.

    The sweep is a partial sweep from ``EAS_START_FREQUENCY`` to
    ``EAS_END_FREQUENCY``, which is why the resulting TMRCA distribution is
    bimodal: swept pairs share a recent ancestor, unswept and mixed pairs do not.
    """
    demography = build_demographic_model(arm, repo_root)
    if mode == "selected":
        models: list[Any] = [
            msprime.SweepGenicSelection(
                position=FOCAL_POSITION_BP,
                start_frequency=EAS_START_FREQUENCY,
                end_frequency=EAS_END_FREQUENCY,
                s=SELECTION_COEFFICIENT,
                dt=EAS_SWEEP_DT,
            ),
            msprime.StandardCoalescent(),
        ]
    elif mode == "neutral":
        models = [msprime.StandardCoalescent()]
    else:
        raise ValueError(f"mode must be 'selected' or 'neutral', got {mode!r}")

    ts = msprime.sim_ancestry(
        samples=[
            msprime.SampleSet(SAMPLE_DIPLOIDS, population=arm.target_population, ploidy=2)
        ],
        demography=demography,
        sequence_length=SEQUENCE_LENGTH_BP,
        recombination_rate=RECOMBINATION_RATE,
        model=models,
        random_seed=seed,
    )
    # Neutral mutations everywhere except the reserved focal base, matching the
    # SLiM arms so the decoder sees the same mutational background.
    rate_map = msprime.RateMap(
        position=np.array(
            [0.0, float(FOCAL_POSITION_BP), float(FOCAL_POSITION_BP) + 1.0,
             float(SEQUENCE_LENGTH_BP)]
        ),
        rate=np.array([MUTATION_RATE, 0.0, MUTATION_RATE]),
    )
    return msprime.sim_mutations(ts, rate=rate_map, random_seed=seed + 1)


# ---------------------------------------------------------------------------
# Script generation and records
# ---------------------------------------------------------------------------


def generate_slim_script(arm: Run5Arm, mode: str, repo_root: str | Path) -> str:
    if arm.engine == "msprime":
        raise ValueError(f"{arm.arm_id} does not use SLiM")
    model = build_demographic_model(arm, repo_root)
    contig = build_contig()
    events = build_extended_events(arm, mode)
    target_id = population_id(model, arm.target_population)
    buffer = io.StringIO()
    with scoped_slim_patch(arm, mode, target_id):
        with redirect_stdout(buffer):
            stdpopsim.get_engine("slim").simulate(
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


def model_record(arm: Run5Arm, mode: str, repo_root: str | Path) -> dict[str, Any]:
    record: dict[str, Any] = {
        "schema": "gamma-smc.run5-model/v1",
        "arm_id": arm.arm_id,
        "mode": mode,
        "engine": arm.engine,
        "generation_time_years": arm.generation_time,
        "target_population": arm.target_population,
        "selection_coefficient": SELECTION_COEFFICIENT if mode == "selected" else 0.0,
        "tick_schedule": tick_schedule(arm),
    }
    if arm.engine == "msprime":
        record["flat_window_effective_size"] = eas_effective_size(repo_root)
        record["flat_window_generations"] = EAS_FLAT_WINDOW_GENERATIONS
        record["epoch_count"] = EAS_EPOCH_COUNT
        record["sweep"] = {
            "start_frequency": EAS_START_FREQUENCY,
            "end_frequency": EAS_END_FREQUENCY,
            "dt": EAS_SWEEP_DT,
        }
    else:
        model = build_demographic_model(arm, repo_root)
        record["demographic_model_id"] = model.id
        record["target_population_id"] = population_id(model, arm.target_population)
        record["nea_selection"] = arm.nea_selection
        record["extended_events"] = serialize_extended_events(
            build_extended_events(arm, mode)
        )
        record["stdpopsim_version"] = stdpopsim.__version__
    return record
