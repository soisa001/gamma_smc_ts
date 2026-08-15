"""Natural-neutral versus ``s=0.001`` Han introgression study.

The focal allele is fixed in the Neanderthal donor immediately before the
2.96% pulse.  Both classes are conditioned only on population-level survival
to the present; terminal population or sample allele frequency is never used
as an acceptance gate.  A uniformly drawn 100-diploid panel is retained from
each 500-diploid sample even when the focal allele is not observed in it.

Every phase is restartable and checksum contracted.  Outputs are isolated
from the older frequency-conditioned focused-selection campaign.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import platform
from collections.abc import Mapping, Sequence
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import contextmanager
from pathlib import Path
from time import perf_counter
from typing import Any

import numpy as np
import pandas as pd
import msprime
import pyslim
import scipy
import stdpopsim
import stdpopsim.slim_engine as slim_engine
import tskit

from .defaults import DEFAULT_CACHE_SIZE
from .decoder import run_within_decoder
from .eas_sweep_analysis import (
    build_diploid_pair_table,
    gamma_summary_to_long,
    tree_truth_profiles,
)
from .eas_sweep_models import (
    ARCHAIC_POPULATION,
    FOCAL_POSITION_BP,
    FOCAL_SITE_ID,
    GENERATION_TIME_YEARS,
    HAN_SPLIT_GENERATIONS,
    INTROGRESSION_PULSE_GENERATIONS,
    INTROGRESSION_RECIPIENT_POPULATION,
    INTROGRESSION_TARGET_POPULATION,
    MUTATION_RATE,
    RECOMBINATION_RATE,
    SEQUENCE_LENGTH_BP,
    build_focal_contig,
    load_ancient_eurasia_model,
)
from .eas_sweep_study import BASE_SEED


SCHEMA_VERSION = "gamma-smc.han-s001-tmrca-study/v1"
SIMULATION_SCHEMA = "gamma-smc.han-s001-tmrca-simulation/v1"
DECODE_SCHEMA = "gamma-smc.han-s001-tmrca-decode/v1"
AGGREGATE_SCHEMA = "gamma-smc.han-s001-tmrca-aggregate/v1"

DEFAULT_RESULTS_DIR = Path("focused_selection_EAS_sim/results/han_s001_tmrca")
DEFAULT_STUDY_DIR = DEFAULT_RESULTS_DIR / "study"
DEFAULT_ANALYSIS_DIR = DEFAULT_RESULTS_DIR / "analysis"
PLAN_RELATIVE_PATH = DEFAULT_STUDY_DIR / "study_plan.json"
WORK_RELATIVE_DIR = Path("focused_selection_EAS_sim/work/han_s001_tmrca")
SIMULATION_STATUS_RELATIVE_PATH = DEFAULT_STUDY_DIR / "simulation_status.tsv"
DECODE_STATUS_RELATIVE_PATH = DEFAULT_STUDY_DIR / "decode_status.tsv"
SCORES_RELATIVE_PATH = DEFAULT_STUDY_DIR / "replicate_scores.tsv"
STUDY_COMPLETION_RELATIVE_PATH = DEFAULT_STUDY_DIR / "study_completion.json"

TMRCA_THRESHOLDS_YEARS = (1_000, 4_500, 10_000, 20_000, 30_000, 40_000, 50_000)
DEFAULT_SELECTION_COEFFICIENTS = {"neutral": 0.0, "selected": 0.001}
DEFAULT_NEUTRAL_REPLICATES = 100
DEFAULT_SELECTED_REPLICATES = 100
DEFAULT_WORKERS = 20
DEFAULT_DECODE_THREADS = 1
DEFAULT_POOL_DIPLOIDS = 500
DEFAULT_PANEL_DIPLOIDS = 100
DEFAULT_SLIM_SCALING_FACTOR = 5.0
DEFAULT_SLIM_BURN_IN = 0.1
DEFAULT_DOMINANCE_COEFFICIENT = 0.5
DEFAULT_OUTPUT_STRIDE_BP = 10_000
LOCAL_WINDOW_HALF_WIDTH_BP = 100_000
INTROGRESSION_PULSE_PROPORTION = 0.0296

SOURCE_PATHS = (
    "python/gamma_smc_aou/han_s001_tmrca_simulation.py",
    "python/gamma_smc_aou/eas_sweep_models.py",
    "python/gamma_smc_aou/eas_sweep_analysis.py",
    "python/gamma_smc_aou/decoder.py",
    "python/gamma_smc_aou/tree_sequence.py",
)
UNIT_COLUMNS = (
    "unit_id",
    "simulation_class",
    "selection_coefficient",
    "replicate_index",
    "seed",
)
SCORE_COLUMNS = (
    "unit_id",
    "simulation_class",
    "selection_coefficient",
    "seed",
    "source",
    "window",
    "threshold_years",
    "score",
    "final_population_af",
    "sample_alt_count",
    "sample_af",
    "sample_detected",
)

_ORIGINAL_ADD_MUT = """// Add `mut_type` mutation at `pos`, to a single individual in `pop`.
function (void)add_mut(object$ mut_type, object$ pop, integer$ pos) {
   targets = sample(pop.genomes, 1);
   targets.addNewDrawnMutation(mut_type, pos);
}"""
_DONOR_FIXED_ADD_MUT = """// Add the focal mutation to every donor genome (donor fixed by design).
function (void)add_mut(object$ mut_type, object$ pop, integer$ pos) {
   targets = sample(pop.genomes, 1);
   if (size(targets) != 1)
       err(\"donor-fixed draw did not select exactly one genome\");
   mut = targets.addNewDrawnMutation(mut_type, pos);
   if (size(mut) != 1)
       err(\"donor-fixed draw did not create exactly one mutation object\");
   remaining = setDifference(pop.genomes, targets);
   remaining.addMutations(mut);
   if (size(sim.mutationsOfType(mut_type)) != 1)
       err(\"donor-fixed draw has more than one extant focal mutation object\");
   source_af = af(mut_type, pop);
   if (source_af != 1.0)
       err(\"donor-fixed draw did not reach source AF 1\");
   metadata.setValue(\"donor_fixed_source_af_at_draw\", source_af);
}"""
_ORIGINAL_END = """// Output tree sequence file and end the simulation.
function (void)end(void) {
    sim.treeSeqOutput(trees_file, metadata=metadata);
    sim.simulationFinished();
}"""
_EXPECTED_SLIM_USER_METADATA_KEYS = {
    "donor_fixed_source_af_at_draw",
    "donor_fixed_final_population_alt_count",
    "donor_fixed_final_population_total_count",
    "donor_fixed_final_population_af",
    "donor_fixed_final_population_survived",
    "donor_fixed_final_population_fixed",
}


def _sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _sha256_source(path: str | Path) -> str:
    text = Path(path).read_text(encoding="utf-8")
    canonical_lf = text.replace("\r\n", "\n").replace("\r", "\n")
    return hashlib.sha256(canonical_lf.encode("utf-8")).hexdigest()


def _canonical_sha256(value: Any) -> str:
    payload = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _atomic_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        temporary.write_text(
            json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_frame(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        frame.to_csv(
            temporary,
            sep="\t",
            index=False,
            lineterminator="\n",
            float_format="%.12g",
        )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _positive_integer(value: Any, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{label} must be a positive integer")
    return value


def _stable_seed(unit_id: str) -> int:
    token = f"{BASE_SEED}:han-s001-tmrca:{unit_id}".encode("utf-8")
    return 1 + int.from_bytes(hashlib.sha256(token).digest()[:8], "big") % (2**31 - 2)


def _source_hashes(repo_root: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    for relative in SOURCE_PATHS:
        path = repo_root / relative
        if not path.is_file():
            raise ValueError(f"study implementation source is absent: {path}")
        records[relative] = _sha256_source(path)
    return records


def _runtime_versions() -> dict[str, str]:
    return {
        "python": platform.python_version(),
        "numpy": str(np.__version__),
        "pandas": str(pd.__version__),
        "scipy": str(scipy.__version__),
        "stdpopsim": str(stdpopsim.__version__),
        "pyslim": str(pyslim.__version__),
        "msprime": str(msprime.__version__),
        "tskit": str(tskit.__version__),
    }


def _slim_patch_contract() -> dict[str, Any]:
    functions = str(slim_engine._slim_functions)
    add_count = functions.count(_ORIGINAL_ADD_MUT)
    end_count = functions.count(_ORIGINAL_END)
    if add_count != 1 or end_count != 1:
        raise RuntimeError(
            "stdpopsim 0.3.0 SLiM helpers no longer match the donor-fixed "
            f"patch contract (add_mut={add_count}, end={end_count})"
        )
    return {
        "stdpopsim_version": str(stdpopsim.__version__),
        "slim_engine_module": "stdpopsim.slim_engine",
        "slim_functions_sha256": hashlib.sha256(functions.encode()).hexdigest(),
        "original_add_mut_sha256": hashlib.sha256(
            _ORIGINAL_ADD_MUT.encode()
        ).hexdigest(),
        "replacement_add_mut_sha256": hashlib.sha256(
            _DONOR_FIXED_ADD_MUT.encode()
        ).hexdigest(),
        "original_end_sha256": hashlib.sha256(_ORIGINAL_END.encode()).hexdigest(),
        "replacement_count": {"add_mut": add_count, "end": end_count},
    }


def plan_path(repo_root: str | Path) -> Path:
    """Return the study plan path without creating or changing anything."""

    return Path(repo_root).resolve() / PLAN_RELATIVE_PATH


def _build_units(n_neutral: int, n_selected: int) -> list[dict[str, Any]]:
    units: list[dict[str, Any]] = []
    for simulation_class, count in (
        ("neutral", n_neutral),
        ("selected", n_selected),
    ):
        for replicate_index in range(count):
            unit_id = f"{simulation_class}_{replicate_index:04d}"
            units.append(
                {
                    "unit_id": unit_id,
                    "simulation_class": simulation_class,
                    "selection_coefficient": DEFAULT_SELECTION_COEFFICIENTS[
                        simulation_class
                    ],
                    "replicate_index": replicate_index,
                    "seed": _stable_seed(unit_id),
                }
            )
    return units


def _plan_payload(repo_root: Path, n_neutral: int, n_selected: int) -> dict[str, Any]:
    units = _build_units(n_neutral, n_selected)
    payload: dict[str, Any] = {
        "schema": SCHEMA_VERSION,
        "n_neutral": n_neutral,
        "n_selected": n_selected,
        "base_seed": BASE_SEED,
        "units": units,
        "design": {
            "demographic_model": "AncientEurasia_9K19",
            "population": INTROGRESSION_TARGET_POPULATION,
            "donor_population": ARCHAIC_POPULATION,
            "introgression_recipient_population": (INTROGRESSION_RECIPIENT_POPULATION),
            "requested_pulse_generations_ago": INTROGRESSION_PULSE_GENERATIONS,
            "realized_pulse_generations_ago_q5": (
                round(INTROGRESSION_PULSE_GENERATIONS / DEFAULT_SLIM_SCALING_FACTOR)
                * DEFAULT_SLIM_SCALING_FACTOR
            ),
            "pulse_proportion": INTROGRESSION_PULSE_PROPORTION,
            "donor_frequency_at_pulse": 1.0,
            "selection_coefficients": DEFAULT_SELECTION_COEFFICIENTS,
            "effective_q_times_selection_coefficients": {
                key: DEFAULT_SLIM_SCALING_FACTOR * value
                for key, value in DEFAULT_SELECTION_COEFFICIENTS.items()
            },
            "dominance_coefficient": DEFAULT_DOMINANCE_COEFFICIENT,
            "selection_interval": "postpulse_through_present",
            "requested_han_split_generations_ago": HAN_SPLIT_GENERATIONS,
            "realized_han_split_generations_ago_q5": (
                round(HAN_SPLIT_GENERATIONS / DEFAULT_SLIM_SCALING_FACTOR)
                * DEFAULT_SLIM_SCALING_FACTOR
            ),
            "conditioning": "population_survival_only",
            "terminal_frequency_conditioning": False,
            "pool_diploids": DEFAULT_POOL_DIPLOIDS,
            "panel_diploids": DEFAULT_PANEL_DIPLOIDS,
            "panel_sampling": "uniform_without_frequency_conditioning",
            "pool_or_panel_undetected_is_accepted": True,
            "sequence_length_bp": SEQUENCE_LENGTH_BP,
            "focal_position_0based": FOCAL_POSITION_BP,
            "local_window_half_width_bp": LOCAL_WINDOW_HALF_WIDTH_BP,
            "output_stride_bp": DEFAULT_OUTPUT_STRIDE_BP,
            "thresholds_years": list(TMRCA_THRESHOLDS_YEARS),
            "generation_time_years": GENERATION_TIME_YEARS,
            "slim_scaling_factor": DEFAULT_SLIM_SCALING_FACTOR,
            "slim_burn_in": DEFAULT_SLIM_BURN_IN,
        },
        "implementation_sha256s": _source_hashes(repo_root),
        "runtime_versions": _runtime_versions(),
        "stdpopsim_patch": _slim_patch_contract(),
    }
    payload["contract_sha256"] = _canonical_sha256(payload)
    return payload


def _validate_plan(payload: Mapping[str, Any], repo_root: Path) -> dict[str, Any]:
    result = dict(payload)
    if result.get("schema") != SCHEMA_VERSION:
        raise ValueError("study plan schema is incompatible")
    observed_hash = result.pop("contract_sha256", None)
    if observed_hash != _canonical_sha256(result):
        raise ValueError("study plan contract checksum failed")
    result["contract_sha256"] = observed_hash
    n_neutral = _positive_integer(result.get("n_neutral"), "n_neutral")
    n_selected = _positive_integer(result.get("n_selected"), "n_selected")
    expected_units = _build_units(n_neutral, n_selected)
    if result.get("units") != expected_units:
        raise ValueError("study plan execution units are incompatible")
    if result.get("implementation_sha256s") != _source_hashes(repo_root):
        raise ValueError("study plan implementation source checksums changed")
    if result.get("stdpopsim_patch") != _slim_patch_contract():
        raise ValueError("study plan stdpopsim patch contract changed")
    design = result.get("design")
    if not isinstance(design, Mapping):
        raise ValueError("study plan design is absent")
    expected = _plan_payload(repo_root, n_neutral, n_selected)
    if result != expected:
        raise ValueError("study plan design differs from the frozen contract")
    return result


def write_plan(
    repo_root: str | Path,
    *,
    n_neutral: int = DEFAULT_NEUTRAL_REPLICATES,
    n_selected: int = DEFAULT_SELECTED_REPLICATES,
    force: bool = False,
) -> Path:
    """Create the immutable execution plan, or validate and reuse it."""

    root = Path(repo_root).resolve()
    n_neutral = _positive_integer(n_neutral, "n_neutral")
    n_selected = _positive_integer(n_selected, "n_selected")
    path = plan_path(root)
    requested = _plan_payload(root, n_neutral, n_selected)
    if path.is_file() and not force:
        existing = load_plan(root)
        if existing != requested:
            raise ValueError(
                "existing study plan differs; use a new study root instead of "
                "reinterpreting execution units"
            )
        return path
    if force and (root / WORK_RELATIVE_DIR).exists():
        if any((root / WORK_RELATIVE_DIR).iterdir()):
            raise ValueError("cannot force-replace a plan after unit work exists")
    _atomic_json(path, requested)
    return path


def load_plan(repo_root: str | Path) -> dict[str, Any]:
    """Read and strictly validate the existing study plan."""

    root = Path(repo_root).resolve()
    path = plan_path(root)
    if not path.is_file():
        raise ValueError(f"study plan is absent: {path}")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError("study plan is unreadable") from error
    if not isinstance(payload, Mapping):
        raise ValueError("study plan must be a JSON object")
    return _validate_plan(payload, root)


def _population_id(model: stdpopsim.DemographicModel, name: str) -> int:
    matches = [
        index
        for index, population in enumerate(model.model.populations)
        if population.name == name
    ]
    if len(matches) != 1:
        raise ValueError(f"demographic population {name!r} is not unique")
    return matches[0]


def _donor_fixed_end(han_population_id: int) -> str:
    return f"""// Record exact population AF, output the tree sequence, and finish.
function (void)end(void) {{
    mut_type = m1;
    muts = sim.mutationsOfType(mut_type);
    if (size(muts) != 1)
        err(\"present-day state does not contain exactly one focal mutation object\");
    han_matches = sim.subpopulations[sim.subpopulations.id == {han_population_id}];
    if (size(han_matches) != 1)
        err(\"present-day state does not contain exactly one Han population\");
    han = han_matches[0];
    total_count = size(han.genomes);
    alt_count = sum(han.genomes.containsMutations(muts[0]));
    if ((total_count <= 0) | (alt_count <= 0) | (alt_count > total_count))
        err(\"present-day focal population counts violate survival\");
    final_af = alt_count / total_count;
    metadata.setValue(\"donor_fixed_final_population_alt_count\", alt_count);
    metadata.setValue(\"donor_fixed_final_population_total_count\", total_count);
    metadata.setValue(\"donor_fixed_final_population_af\", final_af);
    metadata.setValue(\"donor_fixed_final_population_survived\", alt_count > 0);
    metadata.setValue(\"donor_fixed_final_population_fixed\", alt_count == total_count);
    sim.treeSeqOutput(trees_file, metadata=metadata);
    sim.simulationFinished();
}}"""


@contextmanager
def _scoped_donor_fixed_patch(han_population_id: int):
    """Patch only this process's stdpopsim script helpers, then restore them."""

    original_object = slim_engine._slim_functions
    original = str(original_object)
    _slim_patch_contract()
    patched = original.replace(_ORIGINAL_ADD_MUT, _DONOR_FIXED_ADD_MUT, 1)
    patched = patched.replace(_ORIGINAL_END, _donor_fixed_end(han_population_id), 1)
    if patched == original:
        raise RuntimeError("donor-fixed stdpopsim patch made no change")
    slim_engine._slim_functions = patched
    try:
        yield {
            "original_sha256": hashlib.sha256(original.encode()).hexdigest(),
            "patched_sha256": hashlib.sha256(patched.encode()).hexdigest(),
            "han_population_id": han_population_id,
        }
    finally:
        slim_engine._slim_functions = original_object


def _build_extended_events(selection_coefficient: float) -> list[Any]:
    coefficient = float(selection_coefficient)
    if not math.isfinite(coefficient) or coefficient < 0:
        raise ValueError("selection coefficient must be finite and nonnegative")
    q = DEFAULT_SLIM_SCALING_FACTOR
    realized_split = round(HAN_SPLIT_GENERATIONS / q) * q
    events: list[Any] = [
        stdpopsim.DrawMutation(
            time=INTROGRESSION_PULSE_GENERATIONS,
            single_site_id=FOCAL_SITE_ID,
            population=ARCHAIC_POPULATION,
        ),
        stdpopsim.ConditionOnAlleleFrequency(
            start_time=stdpopsim.GenerationAfter(INTROGRESSION_PULSE_GENERATIONS),
            end_time=realized_split,
            single_site_id=FOCAL_SITE_ID,
            population=INTROGRESSION_RECIPIENT_POPULATION,
            op=">",
            allele_frequency=0.0,
        ),
        stdpopsim.ConditionOnAlleleFrequency(
            start_time=stdpopsim.GenerationAfter(realized_split),
            end_time=0.0,
            single_site_id=FOCAL_SITE_ID,
            population=INTROGRESSION_TARGET_POPULATION,
            op=">",
            allele_frequency=0.0,
        ),
    ]
    if coefficient > 0:
        events.extend(
            (
                stdpopsim.ChangeMutationFitness(
                    start_time=INTROGRESSION_PULSE_GENERATIONS,
                    end_time=realized_split,
                    single_site_id=FOCAL_SITE_ID,
                    population=INTROGRESSION_RECIPIENT_POPULATION,
                    selection_coeff=coefficient,
                    dominance_coeff=DEFAULT_DOMINANCE_COEFFICIENT,
                ),
                stdpopsim.ChangeMutationFitness(
                    start_time=realized_split,
                    end_time=0.0,
                    single_site_id=FOCAL_SITE_ID,
                    population=INTROGRESSION_TARGET_POPULATION,
                    selection_coeff=coefficient,
                    dominance_coeff=DEFAULT_DOMINANCE_COEFFICIENT,
                ),
            )
        )
    return events


def _sampled_diploids(ts: tskit.TreeSequence) -> list[tuple[int, tuple[int, int]]]:
    sample_nodes = {int(node) for node in ts.samples()}
    records: list[tuple[int, tuple[int, int]]] = []
    assigned: list[int] = []
    for individual in ts.individuals():
        nodes = tuple(int(node) for node in individual.nodes if node in sample_nodes)
        if not nodes:
            continue
        if len(nodes) != 2:
            raise ValueError("sampled Han individuals must be diploid")
        records.append((int(individual.id), (nodes[0], nodes[1])))
        assigned.extend(nodes)
    if len(assigned) != len(sample_nodes) or set(assigned) != sample_nodes:
        raise ValueError("every sampled Han node must belong to one diploid")
    return records


def _focal_counts_or_zero(ts: tskit.TreeSequence) -> np.ndarray:
    records = _sampled_diploids(ts)
    sites = [
        site
        for site in ts.sites()
        if math.isclose(float(site.position), FOCAL_POSITION_BP, abs_tol=1e-9)
    ]
    if not sites:
        return np.zeros(len(records), dtype=np.int8)
    if len(sites) != 1:
        raise ValueError("focal position contains more than one tree-sequence site")
    site = sites[0]
    derived_states = {str(mutation.derived_state) for mutation in site.mutations}
    if len(derived_states) != 1 or site.ancestral_state in derived_states:
        raise ValueError("focal site does not contain one unique causal derived allele")
    for mutation in site.mutations:
        raw_types: list[Any] = []

        def collect_types(value: Any) -> None:
            if isinstance(value, Mapping):
                if "mutation_type" in value:
                    raw_types.append(value["mutation_type"])
                for nested in value.values():
                    collect_types(nested)
            elif isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
                for nested in value:
                    collect_types(nested)

        collect_types(mutation.metadata)
        if len(raw_types) != 1:
            raise ValueError(
                "each focal mutation row must identify one causal mutation type"
            )
        for raw_type in raw_types:
            while isinstance(raw_type, Sequence) and not isinstance(
                raw_type, (str, bytes)
            ):
                if len(raw_type) != 1:
                    raise ValueError("focal mutation type metadata is not scalar")
                raw_type = raw_type[0]
            if int(raw_type) != 1:
                raise ValueError(
                    "focal site contains a noncausal/background mutation type"
                )
    ordered_nodes = [node for _, nodes in records for node in nodes]
    variant = tskit.Variant(ts, samples=ordered_nodes)
    variant.decode(site.id)
    genotypes = np.asarray(variant.genotypes, dtype=np.int64)
    if np.any(genotypes < 0):
        raise ValueError("focal genotypes contain missing values")
    causal_state = next(iter(derived_states))
    causal_indices = [
        index for index, allele in enumerate(variant.alleles) if allele == causal_state
    ]
    if len(causal_indices) != 1 or causal_indices[0] == 0:
        raise ValueError("focal causal allele is ambiguous in decoded genotypes")
    noncausal_derived = set(genotypes.tolist()).difference({0, causal_indices[0]})
    if noncausal_derived:
        raise ValueError("focal genotypes contain a background derived allele")
    counts = (genotypes.reshape(len(records), 2) == causal_indices[0]).sum(axis=1)
    return counts.astype(np.int8, copy=False)


def _uniform_panel(
    pool_ts: tskit.TreeSequence,
    pool_counts: np.ndarray,
    *,
    seed: int,
) -> tuple[tskit.TreeSequence, np.ndarray, list[int]]:
    pool = _sampled_diploids(pool_ts)
    if len(pool) != DEFAULT_POOL_DIPLOIDS or len(pool_counts) != len(pool):
        raise ValueError("SLiM candidate pool does not contain 500 diploids")
    rng = np.random.default_rng(seed)
    chosen_rows = sorted(
        int(value)
        for value in rng.choice(
            len(pool), size=DEFAULT_PANEL_DIPLOIDS, replace=False
        ).tolist()
    )
    chosen_nodes = np.asarray(
        [node for row in chosen_rows for node in pool[row][1]], dtype=np.int32
    )
    tables = pool_ts.dump_tables()
    flags = np.asarray(tables.nodes.flags, dtype=np.uint32).copy()
    flags &= np.bitwise_not(np.uint32(tskit.NODE_IS_SAMPLE))
    flags[chosen_nodes] |= np.uint32(tskit.NODE_IS_SAMPLE)
    nodes = tables.nodes
    nodes.set_columns(
        flags=flags,
        time=nodes.time,
        population=nodes.population,
        individual=nodes.individual,
        metadata=nodes.metadata,
        metadata_offset=nodes.metadata_offset,
    )
    panel_ts = tables.tree_sequence()
    panel_counts = np.asarray(pool_counts[chosen_rows], dtype=np.int8)
    if len(_sampled_diploids(panel_ts)) != DEFAULT_PANEL_DIPLOIDS:
        raise RuntimeError("uniform panel construction changed sample cardinality")
    return panel_ts, panel_counts, chosen_rows


def _slim_user_metadata(metadata: Any) -> Mapping[str, Any]:
    if not isinstance(metadata, Mapping):
        raise ValueError("tree-sequence metadata must be a mapping")
    slim = metadata.get("SLiM")
    if not isinstance(slim, Mapping):
        raise ValueError("tree metadata lacks exact SLiM mapping")
    user = slim.get("user_metadata")
    if not isinstance(user, Mapping):
        raise ValueError("tree metadata lacks exact SLiM.user_metadata mapping")
    missing = sorted(_EXPECTED_SLIM_USER_METADATA_KEYS.difference(user))
    if missing:
        raise ValueError(
            "SLiM.user_metadata lacks donor/final keys: " + ", ".join(missing)
        )
    for key in _EXPECTED_SLIM_USER_METADATA_KEYS:
        occurrences = 0

        def count(candidate: Any) -> None:
            nonlocal occurrences
            if isinstance(candidate, Mapping):
                occurrences += int(key in candidate)
                for nested in candidate.values():
                    count(nested)
            elif isinstance(candidate, Sequence) and not isinstance(
                candidate, (str, bytes)
            ):
                for nested in candidate:
                    count(nested)

        count(metadata)
        if occurrences != 1:
            raise ValueError(
                f"tree metadata key {key} occurs outside its one exact path"
            )
    return user


def _metadata_scalar(metadata: Any, key: str) -> float:
    if key not in _EXPECTED_SLIM_USER_METADATA_KEYS:
        raise ValueError(f"unexpected contracted SLiM metadata key: {key}")
    value = _slim_user_metadata(metadata)[key]
    while isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        if len(value) != 1:
            raise ValueError(f"tree metadata {key} is not scalar")
        value = value[0]
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"tree metadata {key} is nonfinite")
    return result


def _metadata_integer(metadata: Any, key: str) -> int:
    value = _metadata_scalar(metadata, key)
    if value != math.floor(value):
        raise ValueError(f"tree metadata {key} is not an integer")
    return int(value)


def _portable_binary_identity(path: Path, repo_root: Path) -> str:
    try:
        return path.resolve().relative_to(repo_root.resolve()).as_posix()
    except ValueError as error:
        raise ValueError(
            "production binaries must be inside repo_root for a resolvable contract"
        ) from error


def _resolve_binary_identity(identity: Any, repo_root: Path, label: str) -> Path:
    relative = Path(str(identity))
    if relative.is_absolute() or ".." in relative.parts or not relative.parts:
        raise ValueError(f"{label} binary identity is not repo relative")
    resolved = (repo_root / relative).resolve()
    try:
        resolved.relative_to(repo_root.resolve())
    except ValueError as error:
        raise ValueError(f"{label} binary identity escapes repo_root") from error
    if not resolved.is_file():
        raise ValueError(f"{label} contracted binary is absent: {resolved}")
    return resolved


def _output_record(path: Path, *, rows: int | None = None) -> dict[str, Any]:
    record: dict[str, Any] = {
        "path": path.name,
        "sha256": _sha256_file(path),
        "size_bytes": path.stat().st_size,
    }
    if rows is not None:
        record["rows"] = int(rows)
    return record


def _validate_output_record(
    record: Any, path: Path, *, expected_rows: int | None = None
) -> None:
    expected_fields = {"path", "sha256", "size_bytes"}
    if expected_rows is not None:
        expected_fields.add("rows")
    if not isinstance(record, Mapping) or set(record) != expected_fields:
        raise ValueError(f"output record fields changed: {path.name}")
    digest = str(record.get("sha256", ""))
    if len(digest) != 64 or any(char not in "0123456789abcdef" for char in digest):
        raise ValueError(f"output record SHA-256 is invalid: {path.name}")
    if (
        not path.is_file()
        or record.get("path") != path.name
        or _sha256_file(path) != digest
        or path.stat().st_size != int(record.get("size_bytes", -1))
    ):
        raise ValueError(f"output checksum failed: {path.name}")
    if expected_rows is not None and int(record.get("rows", -1)) != expected_rows:
        raise ValueError(f"output row count record changed: {path.name}")


def _unit_dir(repo_root: Path, unit_id: str) -> Path:
    return repo_root / WORK_RELATIVE_DIR / unit_id


def _simulation_contract(
    plan: Mapping[str, Any],
    unit: Mapping[str, Any],
    slim_bin: Path,
    repo_root: Path,
) -> dict[str, Any]:
    return {
        "schema": SIMULATION_SCHEMA,
        "plan_contract_sha256": plan["contract_sha256"],
        "unit": dict(unit),
        "design": dict(plan["design"]),
        "implementation_sha256s": dict(plan["implementation_sha256s"]),
        "runtime_versions": dict(plan["runtime_versions"]),
        "stdpopsim_patch": dict(plan["stdpopsim_patch"]),
        "slim_binary": {
            "identity": _portable_binary_identity(slim_bin, repo_root),
            "sha256": _sha256_file(slim_bin),
            "size_bytes": slim_bin.stat().st_size,
        },
    }


def _validate_simulation_completion(
    unit_dir: Path, expected_contract: Mapping[str, Any]
) -> dict[str, Any]:
    path = unit_dir / "simulation_complete.json"
    try:
        completion = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError("simulation completion is absent or unreadable") from error
    expected_completion_fields = {
        "schema",
        "status",
        "unit_id",
        "contract",
        "contract_sha256",
        "metadata",
        "elapsed_seconds",
        "outputs",
    }
    if set(completion) != expected_completion_fields:
        raise ValueError("simulation completion top-level fields changed")
    if completion.get("schema") != SIMULATION_SCHEMA:
        raise ValueError("simulation completion schema is incompatible")
    if completion.get("status") != "complete":
        raise ValueError("simulation completion status is not complete")
    if completion.get("unit_id") != expected_contract["unit"]["unit_id"]:
        raise ValueError("simulation completion unit_id changed")
    elapsed = float(completion.get("elapsed_seconds", math.nan))
    if not math.isfinite(elapsed) or elapsed < 0:
        raise ValueError("simulation completion elapsed time is invalid")
    if completion.get("contract") != expected_contract:
        raise ValueError("simulation completion contract changed")
    if completion.get("contract_sha256") != _canonical_sha256(expected_contract):
        raise ValueError("simulation completion contract checksum failed")
    expected_names = {
        "simulation.trees",
        "sample_manifest.tsv",
        "truth_overall.tsv",
        "simulation_contract.json",
        "slim.tsv",
    }
    outputs = completion.get("outputs")
    if not isinstance(outputs, Mapping) or set(outputs) != expected_names:
        raise ValueError("simulation completion output inventory changed")
    manifest = pd.read_csv(unit_dir / "sample_manifest.tsv", sep="\t")
    truth = pd.read_csv(unit_dir / "truth_overall.tsv", sep="\t")
    if len(manifest) != DEFAULT_PANEL_DIPLOIDS:
        raise ValueError("simulation sample manifest is not 100 diploids")
    if len(truth) != (SEQUENCE_LENGTH_BP // DEFAULT_OUTPUT_STRIDE_BP) * len(
        TMRCA_THRESHOLDS_YEARS
    ):
        raise ValueError("simulation truth profile row count changed")
    for name, record in outputs.items():
        rows = (
            len(manifest)
            if name == "sample_manifest.tsv"
            else len(truth)
            if name == "truth_overall.tsv"
            else None
        )
        _validate_output_record(record, unit_dir / name, expected_rows=rows)
    _validate_exact_profile_grid(truth, "tree-truth")
    metadata = completion.get("metadata", {})
    expected_metadata_fields = {
        "donor_source_af_at_draw",
        "final_population_af",
        "final_population_alt_count",
        "final_population_total_count",
        "final_population_survived",
        "final_population_fixed",
        "pool_alt_count",
        "pool_af",
        "pool_detected",
        "sample_alt_count",
        "sample_af",
        "sample_detected",
        "panel_seed",
        "chosen_pool_rows_zero_based",
        "panel_frequency_used_as_acceptance_gate",
        "population_survival_only",
        "requested_pulse_generations_ago",
        "realized_pulse_generations_ago_q5",
        "requested_han_split_generations_ago",
        "realized_han_split_generations_ago_q5",
        "effective_q_times_selection_coefficient",
        "patch",
    }
    if not isinstance(metadata, Mapping) or set(metadata) != expected_metadata_fields:
        raise ValueError("simulation completion metadata fields changed")
    patch = metadata["patch"]
    if not isinstance(patch, Mapping) or set(patch) != {
        "original_sha256",
        "patched_sha256",
        "han_population_id",
    }:
        raise ValueError("simulation patch metadata fields changed")
    final_af = float(metadata.get("final_population_af", math.nan))
    final_alt = int(metadata.get("final_population_alt_count", -1))
    final_total = int(metadata.get("final_population_total_count", -1))
    if not 0.0 < final_af <= 1.0:
        raise ValueError("population-survival simulation has invalid final AF")
    if (
        final_alt < 1
        or final_total < final_alt
        or not math.isclose(final_af, final_alt / final_total, abs_tol=1e-12)
    ):
        raise ValueError("final population AF disagrees with exact allele counts")
    for boolean_key in (
        "final_population_survived",
        "final_population_fixed",
        "pool_detected",
        "sample_detected",
        "panel_frequency_used_as_acceptance_gate",
        "population_survival_only",
    ):
        if not isinstance(metadata[boolean_key], bool):
            raise ValueError(
                f"simulation metadata {boolean_key} must be a strict boolean"
            )
    if (
        not metadata["final_population_survived"]
        or metadata["final_population_fixed"] != (final_alt == final_total)
        or not math.isclose(
            float(metadata["donor_source_af_at_draw"]), 1.0, abs_tol=1e-12
        )
    ):
        raise ValueError("donor fixation/final survival metadata changed")
    pool_alt = int(metadata["pool_alt_count"])
    pool_af = float(metadata["pool_af"])
    if (
        not 0 <= pool_alt <= 2 * DEFAULT_POOL_DIPLOIDS
        or not math.isclose(
            pool_af, pool_alt / (2 * DEFAULT_POOL_DIPLOIDS), abs_tol=1e-12
        )
        or metadata["pool_detected"] != (pool_alt > 0)
    ):
        raise ValueError("pool AF/count/detection metadata is inconsistent")
    sample_alt = int(metadata["sample_alt_count"])
    sample_af = float(metadata["sample_af"])
    manifest_alt = int(manifest["focal_selected_allele_count"].sum())
    if (
        sample_alt != manifest_alt
        or not 0 <= sample_alt <= 2 * DEFAULT_PANEL_DIPLOIDS
        or not math.isclose(
            sample_af, sample_alt / (2 * DEFAULT_PANEL_DIPLOIDS), abs_tol=1e-12
        )
        or bool(metadata["sample_detected"]) != (sample_alt > 0)
    ):
        raise ValueError("sample AF metadata disagrees with the sample manifest")
    chosen = metadata["chosen_pool_rows_zero_based"]
    if (
        not isinstance(chosen, list)
        or len(chosen) != DEFAULT_PANEL_DIPLOIDS
        or chosen != sorted(set(int(value) for value in chosen))
        or min(chosen) < 0
        or max(chosen) >= DEFAULT_POOL_DIPLOIDS
    ):
        raise ValueError("uniform panel row selection metadata is invalid")
    if bool(metadata["panel_frequency_used_as_acceptance_gate"]) or not bool(
        metadata["population_survival_only"]
    ):
        raise ValueError("simulation ascertainment semantics changed")
    design = expected_contract["design"]
    unit = expected_contract["unit"]
    expected_semantics = {
        "requested_pulse_generations_ago": design["requested_pulse_generations_ago"],
        "realized_pulse_generations_ago_q5": design[
            "realized_pulse_generations_ago_q5"
        ],
        "requested_han_split_generations_ago": design[
            "requested_han_split_generations_ago"
        ],
        "realized_han_split_generations_ago_q5": design[
            "realized_han_split_generations_ago_q5"
        ],
        "effective_q_times_selection_coefficient": (
            design["slim_scaling_factor"] * unit["selection_coefficient"]
        ),
    }
    for key, expected_value in expected_semantics.items():
        if not math.isclose(float(metadata[key]), float(expected_value), abs_tol=1e-12):
            raise ValueError(f"simulation metadata {key} changed")

    tree = tskit.load(unit_dir / "simulation.trees")
    if tree.num_samples != 2 * DEFAULT_PANEL_DIPLOIDS or not math.isclose(
        float(tree.sequence_length), SEQUENCE_LENGTH_BP, abs_tol=1e-9
    ):
        raise ValueError("simulation tree sample/sequence dimensions changed")
    tree_source_af = _metadata_scalar(tree.metadata, "donor_fixed_source_af_at_draw")
    tree_alt = _metadata_integer(
        tree.metadata, "donor_fixed_final_population_alt_count"
    )
    tree_total = _metadata_integer(
        tree.metadata, "donor_fixed_final_population_total_count"
    )
    tree_af = _metadata_scalar(tree.metadata, "donor_fixed_final_population_af")
    tree_survived = bool(
        _metadata_integer(tree.metadata, "donor_fixed_final_population_survived")
    )
    tree_fixed = bool(
        _metadata_integer(tree.metadata, "donor_fixed_final_population_fixed")
    )
    if (
        not math.isclose(tree_source_af, 1.0, abs_tol=1e-12)
        or tree_alt != final_alt
        or tree_total != final_total
        or not math.isclose(tree_af, final_af, abs_tol=1e-12)
        or tree_survived is not True
        or tree_fixed != metadata["final_population_fixed"]
    ):
        raise ValueError("completion metadata disagrees with exact SLiM tree metadata")
    tree_counts = _focal_counts_or_zero(tree)
    manifest_counts = manifest["focal_selected_allele_count"].to_numpy(dtype=int)
    if not np.array_equal(tree_counts.astype(int), manifest_counts):
        raise ValueError("sample manifest focal counts disagree with simulation tree")
    return completion


def _simulate_unit(payload: Mapping[str, Any]) -> dict[str, Any]:
    from .focused_selection_simulation import UnitLockHeld, exclusive_unit_lock

    repo_root = Path(payload["repo_root"])
    plan = dict(payload["plan"])
    unit = dict(payload["unit"])
    slim_bin = Path(payload["slim_bin"])
    directory = _unit_dir(repo_root, str(unit["unit_id"]))
    contract = _simulation_contract(plan, unit, slim_bin, repo_root)
    directory.mkdir(parents=True, exist_ok=True)
    try:
        with exclusive_unit_lock(
            directory / ".simulation.lock",
            unit_id=str(unit["unit_id"]),
            phase="han_s001_simulate",
        ):
            if (directory / "simulation_complete.json").is_file():
                completion = _validate_simulation_completion(directory, contract)
                return {
                    "unit_id": unit["unit_id"],
                    "status": "cached",
                    "error": "",
                    "final_population_af": completion["metadata"][
                        "final_population_af"
                    ],
                }
            _atomic_json(directory / "simulation_contract.json", contract)
            started = perf_counter()
            model = load_ancient_eurasia_model()
            han_population_id = _population_id(model, INTROGRESSION_TARGET_POPULATION)
            contig = build_focal_contig()
            events = _build_extended_events(float(unit["selection_coefficient"]))
            engine = stdpopsim.get_engine("slim")
            logfile = directory / "slim.tsv"
            with _scoped_donor_fixed_patch(han_population_id) as patch_record:
                pool_ts = engine.simulate(
                    model,
                    contig,
                    {INTROGRESSION_TARGET_POPULATION: DEFAULT_POOL_DIPLOIDS},
                    seed=int(unit["seed"]),
                    extended_events=events,
                    slim_path=str(slim_bin),
                    slim_scaling_factor=DEFAULT_SLIM_SCALING_FACTOR,
                    slim_burn_in=DEFAULT_SLIM_BURN_IN,
                    logfile=str(logfile),
                    logfile_interval=100,
                    keep_mutation_ids_as_alleles=False,
                )
            source_af = _metadata_scalar(
                pool_ts.metadata, "donor_fixed_source_af_at_draw"
            )
            final_population_af = _metadata_scalar(
                pool_ts.metadata, "donor_fixed_final_population_af"
            )
            final_population_alt_count = _metadata_integer(
                pool_ts.metadata, "donor_fixed_final_population_alt_count"
            )
            final_population_total_count = _metadata_integer(
                pool_ts.metadata, "donor_fixed_final_population_total_count"
            )
            final_population_survived = bool(
                _metadata_integer(
                    pool_ts.metadata, "donor_fixed_final_population_survived"
                )
            )
            final_population_fixed = bool(
                _metadata_integer(
                    pool_ts.metadata, "donor_fixed_final_population_fixed"
                )
            )
            if not math.isclose(source_af, 1.0, abs_tol=1e-12):
                raise ValueError("patched focal mutation was not donor fixed")
            if not 0.0 < final_population_af <= 1.0:
                raise ValueError("recipient survival condition failed at present")
            if (
                not final_population_survived
                or final_population_alt_count < 1
                or final_population_total_count < final_population_alt_count
                or not math.isclose(
                    final_population_af,
                    final_population_alt_count / final_population_total_count,
                    abs_tol=1e-12,
                )
                or final_population_fixed
                != (final_population_alt_count == final_population_total_count)
            ):
                raise ValueError("exact final population counts disagree with AF flags")
            pool_counts = _focal_counts_or_zero(pool_ts)
            panel_seed = _stable_seed(str(unit["unit_id"]) + ":panel")
            panel_ts, panel_counts, chosen_rows = _uniform_panel(
                pool_ts, pool_counts, seed=panel_seed
            )
            pair_table = build_diploid_pair_table(
                panel_ts,
                FOCAL_POSITION_BP,
                focal_genotype_counts=panel_counts,
            )
            tree_path = directory / "simulation.trees"
            temporary_tree = tree_path.with_name(tree_path.name + f".tmp.{os.getpid()}")
            try:
                panel_ts.dump(temporary_tree)
                os.replace(temporary_tree, tree_path)
            finally:
                temporary_tree.unlink(missing_ok=True)
            _atomic_frame(directory / "sample_manifest.tsv", pair_table)
            positions = np.arange(
                0,
                SEQUENCE_LENGTH_BP,
                DEFAULT_OUTPUT_STRIDE_BP,
                dtype=float,
            )
            truth = tree_truth_profiles(
                panel_ts,
                positions,
                pair_table,
                thresholds_years=TMRCA_THRESHOLDS_YEARS,
                generation_time_years=GENERATION_TIME_YEARS,
                require_nonempty=False,
            )
            truth = truth[truth["genotype_class"] == "overall"].reset_index(drop=True)
            _atomic_frame(directory / "truth_overall.tsv", truth)
            sample_alt_count = int(panel_counts.sum())
            metadata = {
                "donor_source_af_at_draw": source_af,
                "final_population_af": final_population_af,
                "final_population_alt_count": final_population_alt_count,
                "final_population_total_count": final_population_total_count,
                "final_population_survived": final_population_survived,
                "final_population_fixed": final_population_fixed,
                "pool_alt_count": int(pool_counts.sum()),
                "pool_af": float(pool_counts.sum() / (2 * len(pool_counts))),
                "pool_detected": bool(pool_counts.sum() > 0),
                "sample_alt_count": sample_alt_count,
                "sample_af": float(sample_alt_count / (2 * DEFAULT_PANEL_DIPLOIDS)),
                "sample_detected": bool(sample_alt_count > 0),
                "panel_seed": panel_seed,
                "chosen_pool_rows_zero_based": chosen_rows,
                "panel_frequency_used_as_acceptance_gate": False,
                "population_survival_only": True,
                "requested_pulse_generations_ago": (INTROGRESSION_PULSE_GENERATIONS),
                "realized_pulse_generations_ago_q5": (
                    round(INTROGRESSION_PULSE_GENERATIONS / DEFAULT_SLIM_SCALING_FACTOR)
                    * DEFAULT_SLIM_SCALING_FACTOR
                ),
                "requested_han_split_generations_ago": HAN_SPLIT_GENERATIONS,
                "realized_han_split_generations_ago_q5": (
                    round(HAN_SPLIT_GENERATIONS / DEFAULT_SLIM_SCALING_FACTOR)
                    * DEFAULT_SLIM_SCALING_FACTOR
                ),
                "effective_q_times_selection_coefficient": (
                    DEFAULT_SLIM_SCALING_FACTOR * float(unit["selection_coefficient"])
                ),
                "patch": patch_record,
            }
            outputs = {
                name: _output_record(
                    directory / name,
                    rows=(
                        len(pair_table)
                        if name == "sample_manifest.tsv"
                        else len(truth)
                        if name == "truth_overall.tsv"
                        else None
                    ),
                )
                for name in (
                    "simulation.trees",
                    "sample_manifest.tsv",
                    "truth_overall.tsv",
                    "simulation_contract.json",
                    "slim.tsv",
                )
            }
            completion = {
                "schema": SIMULATION_SCHEMA,
                "status": "complete",
                "unit_id": unit["unit_id"],
                "contract": contract,
                "contract_sha256": _canonical_sha256(contract),
                "metadata": metadata,
                "elapsed_seconds": perf_counter() - started,
                "outputs": outputs,
            }
            _atomic_json(directory / "simulation_complete.json", completion)
            (directory / "simulation_failed.json").unlink(missing_ok=True)
            return {
                "unit_id": unit["unit_id"],
                "status": "complete",
                "error": "",
                "final_population_af": final_population_af,
            }
    except UnitLockHeld as error:
        return {
            "unit_id": unit["unit_id"],
            "status": "locked",
            "error": str(error),
            "final_population_af": math.nan,
        }
    except Exception as error:
        failure = {
            "schema": SIMULATION_SCHEMA,
            "status": "failed",
            "unit_id": unit["unit_id"],
            "error_type": type(error).__name__,
            "error": str(error),
            "contract_sha256": _canonical_sha256(contract),
        }
        _atomic_json(directory / "simulation_failed.json", failure)
        return {
            "unit_id": unit["unit_id"],
            "status": "failed",
            "error": f"{type(error).__name__}: {error}",
            "final_population_af": math.nan,
        }


def simulate_study(
    repo_root: str | Path,
    slim_bin: str | Path,
    *,
    workers: int = DEFAULT_WORKERS,
) -> pd.DataFrame:
    """Run or resume all population-survival simulation units."""

    root = Path(repo_root).resolve()
    plan = load_plan(root)
    binary = Path(slim_bin).resolve()
    if not binary.is_file():
        raise ValueError(f"SLiM binary is absent: {binary}")
    n_workers = _positive_integer(workers, "workers")
    payloads = [
        {
            "repo_root": str(root),
            "plan": plan,
            "unit": unit,
            "slim_bin": str(binary),
        }
        for unit in plan["units"]
    ]
    results: list[dict[str, Any]] = []
    if n_workers == 1:
        results = [_simulate_unit(payload) for payload in payloads]
    else:
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = [executor.submit(_simulate_unit, payload) for payload in payloads]
            for future in as_completed(futures):
                results.append(future.result())
    order = {unit["unit_id"]: index for index, unit in enumerate(plan["units"])}
    frame = pd.DataFrame(results).sort_values(
        "unit_id", key=lambda values: values.map(order)
    )
    _atomic_frame(root / SIMULATION_STATUS_RELATIVE_PATH, frame)
    return frame.reset_index(drop=True)


def simulation_status(repo_root: str | Path) -> pd.DataFrame:
    """Return read-only checksum-validated per-unit simulation status."""

    root = Path(repo_root).resolve()
    plan = load_plan(root)
    rows: list[dict[str, Any]] = []
    for unit in plan["units"]:
        directory = _unit_dir(root, str(unit["unit_id"]))
        contract_path = directory / "simulation_contract.json"
        if (
            directory / "simulation_complete.json"
        ).is_file() and contract_path.is_file():
            try:
                contract = json.loads(contract_path.read_text(encoding="utf-8"))
                binary = _resolve_binary_identity(
                    contract["slim_binary"]["identity"], root, "SLiM"
                )
                expected_contract = _simulation_contract(plan, unit, binary, root)
                if contract != expected_contract:
                    raise ValueError("simulation contract differs from current plan")
                completion = _validate_simulation_completion(
                    directory, expected_contract
                )
                status = "complete"
                error = ""
                final_af = completion["metadata"]["final_population_af"]
            except (
                KeyError,
                OSError,
                TypeError,
                ValueError,
                json.JSONDecodeError,
            ) as validation_error:
                status = "invalid"
                error = str(validation_error)
                final_af = math.nan
        elif (directory / "simulation_failed.json").is_file():
            status = "failed"
            failure = json.loads(
                (directory / "simulation_failed.json").read_text(encoding="utf-8")
            )
            error = str(failure.get("error", ""))
            final_af = math.nan
        elif directory.exists():
            status = "partial"
            error = ""
            final_af = math.nan
        else:
            status = "not_started"
            error = ""
            final_af = math.nan
        rows.append(
            {
                **unit,
                "status": status,
                "error": error,
                "final_population_af": final_af,
            }
        )
    return pd.DataFrame(rows)


def _present_han_ne() -> float:
    model = load_ancient_eurasia_model()
    values = [
        float(population.initial_size)
        for population in model.model.populations
        if population.name == INTROGRESSION_TARGET_POPULATION
    ]
    if len(values) != 1 or not math.isfinite(values[0]) or values[0] <= 0:
        raise ValueError("AncientEurasia Han present-day Ne is invalid")
    return values[0]


def _validate_exact_profile_grid(frame: pd.DataFrame, label: str) -> None:
    required = {"position_0based", "threshold_years"}
    if not required.issubset(frame.columns):
        raise ValueError(f"{label} profile lacks position/threshold columns")
    positions = frame["position_0based"].to_numpy(dtype=float)
    thresholds = frame["threshold_years"].to_numpy(dtype=float)
    if (
        not np.all(np.isfinite(positions))
        or not np.all(np.isfinite(thresholds))
        or not np.all(positions == np.floor(positions))
        or not np.all(thresholds == np.floor(thresholds))
    ):
        raise ValueError(f"{label} profile grid contains noninteger values")
    observed = pd.MultiIndex.from_arrays(
        [positions.astype(np.int64), thresholds.astype(np.int64)]
    )
    expected_positions = np.arange(
        0, SEQUENCE_LENGTH_BP, DEFAULT_OUTPUT_STRIDE_BP, dtype=np.int64
    )
    expected = pd.MultiIndex.from_product([expected_positions, TMRCA_THRESHOLDS_YEARS])
    if (
        len(observed) != len(expected)
        or observed.has_duplicates
        or set(observed.tolist()) != set(expected.tolist())
    ):
        raise ValueError(
            f"{label} profile must contain exactly one row for every 10 kb "
            "position and TMRCA threshold"
        )


def _validate_profile_estimand(frame: pd.DataFrame, label: str) -> None:
    required = {
        "position_0based",
        "position_1based",
        "threshold_years",
        "threshold_generations",
        "generation_time_years",
        "n_pairs",
    }
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(
            f"{label} profile lacks estimand columns: {', '.join(missing)}"
        )
    n_pairs = pd.to_numeric(frame["n_pairs"], errors="coerce").to_numpy(dtype=float)
    if (
        not np.all(np.isfinite(n_pairs))
        or not np.all(n_pairs == np.floor(n_pairs))
        or not np.all(n_pairs == DEFAULT_PANEL_DIPLOIDS)
    ):
        raise ValueError(
            f"{label} profile must use exactly {DEFAULT_PANEL_DIPLOIDS} pairs"
        )
    positions_0 = frame["position_0based"].to_numpy(dtype=float)
    positions_1 = frame["position_1based"].to_numpy(dtype=float)
    generation_time = frame["generation_time_years"].to_numpy(dtype=float)
    thresholds_years = frame["threshold_years"].to_numpy(dtype=float)
    thresholds_generations = frame["threshold_generations"].to_numpy(dtype=float)
    if not np.allclose(positions_1, positions_0 + 1, rtol=0.0, atol=1e-9):
        raise ValueError(f"{label} one-based positions are inconsistent")
    if not np.allclose(generation_time, GENERATION_TIME_YEARS, rtol=0.0, atol=1e-12):
        raise ValueError(f"{label} generation-time semantics changed")
    if not np.allclose(
        thresholds_generations,
        thresholds_years / GENERATION_TIME_YEARS,
        rtol=0.0,
        atol=1e-12,
    ):
        raise ValueError(f"{label} threshold-generation semantics changed")


def _decode_contract(
    plan: Mapping[str, Any],
    unit: Mapping[str, Any],
    unit_dir: Path,
    decoder_bin: Path,
    repo_root: Path,
    threads: int,
) -> dict[str, Any]:
    if threads != DEFAULT_DECODE_THREADS:
        raise ValueError(
            f"this study requires exactly {DEFAULT_DECODE_THREADS} decode thread"
        )
    simulation_completion = unit_dir / "simulation_complete.json"
    tree = unit_dir / "simulation.trees"
    manifest = unit_dir / "sample_manifest.tsv"
    return {
        "schema": DECODE_SCHEMA,
        "plan_contract_sha256": plan["contract_sha256"],
        "unit": dict(unit),
        "simulation_completion_sha256": _sha256_file(simulation_completion),
        "tree_sha256": _sha256_file(tree),
        "sample_manifest_sha256": _sha256_file(manifest),
        "decoder_binary": {
            "identity": _portable_binary_identity(decoder_bin, repo_root),
            "sha256": _sha256_file(decoder_bin),
            "size_bytes": decoder_bin.stat().st_size,
        },
        "runtime_versions": dict(plan["runtime_versions"]),
        "settings": {
            "input_format": "trees_via_streamed_vcf",
            "pair_selector": "within_diploid_overall",
            "n_pairs": DEFAULT_PANEL_DIPLOIDS,
            "scaled_mutation_rate": 4 * _present_han_ne() * MUTATION_RATE,
            "unscaled_mutation_rate": MUTATION_RATE,
            "recombination_to_mutation_ratio": (RECOMBINATION_RATE / MUTATION_RATE),
            "thresholds_years": list(TMRCA_THRESHOLDS_YEARS),
            "generation_time_years": GENERATION_TIME_YEARS,
            "output_at_stride": DEFAULT_OUTPUT_STRIDE_BP,
            "output_at_hets": False,
            "cache_size": DEFAULT_CACHE_SIZE,
            "pair_block": 256,
            "exp10": "accurate",
            "backward_alignment": "fixed",
            "recent_call": "median",
            "recent_call_probability": 0.5,
            "threads": threads,
        },
    }


def _gamma_scores(summary: pd.DataFrame) -> pd.DataFrame:
    long = gamma_summary_to_long(
        summary,
        genotype_class="overall",
        generation_time_years=GENERATION_TIME_YEARS,
        source="gamma_smc",
    )
    _validate_exact_profile_grid(long, "Gamma")
    _validate_profile_estimand(long, "Gamma")
    positions = long["position_0based"].to_numpy(dtype=float)
    focal = np.isclose(positions, FOCAL_POSITION_BP, rtol=0.0, atol=1e-9)
    if focal.sum() != len(TMRCA_THRESHOLDS_YEARS):
        raise ValueError("Gamma output does not contain one focal row per threshold")
    local = (positions >= FOCAL_POSITION_BP - LOCAL_WINDOW_HALF_WIDTH_BP) & (
        positions <= FOCAL_POSITION_BP + LOCAL_WINDOW_HALF_WIDTH_BP
    )
    rows: list[dict[str, Any]] = []
    for window, selected in (("focal", focal), ("local_100kb", local)):
        frame = long.loc[selected]
        for threshold in TMRCA_THRESHOLDS_YEARS:
            values = frame.loc[
                np.isclose(
                    frame["threshold_years"].to_numpy(dtype=float),
                    threshold,
                    rtol=0.0,
                    atol=1e-9,
                ),
                "p_tmrca_lt_threshold",
            ].to_numpy(dtype=float)
            if not len(values) or not np.all(np.isfinite(values)):
                raise ValueError("Gamma score window has missing posterior values")
            expected_count = 1 if window == "focal" else 21
            if len(values) != expected_count:
                raise ValueError(
                    f"Gamma {window} window must contain {expected_count} positions"
                )
            score = float(values.mean())
            if not 0.0 <= score <= 1.0:
                raise ValueError("Gamma posterior probability lies outside [0, 1]")
            rows.append(
                {
                    "source": "gamma_smc",
                    "window": window,
                    "threshold_years": int(threshold),
                    "score": score,
                    "n_output_positions": len(values),
                }
            )
    return pd.DataFrame(rows)


def _validate_decode_completion(
    unit_dir: Path, expected_contract: Mapping[str, Any]
) -> dict[str, Any]:
    decoded = unit_dir / "decoded"
    path = decoded / "decode_complete.json"
    try:
        completion = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError("decode completion is absent or unreadable") from error
    expected_completion_fields = {
        "schema",
        "status",
        "unit_id",
        "contract",
        "contract_sha256",
        "elapsed_seconds",
        "outputs",
    }
    if set(completion) != expected_completion_fields:
        raise ValueError("decode completion top-level fields changed")
    if completion.get("schema") != DECODE_SCHEMA:
        raise ValueError("decode completion schema is incompatible")
    if completion.get("status") != "complete":
        raise ValueError("decode completion status is not complete")
    if completion.get("unit_id") != expected_contract["unit"]["unit_id"]:
        raise ValueError("decode completion unit_id changed")
    elapsed = float(completion.get("elapsed_seconds", math.nan))
    if not math.isfinite(elapsed) or elapsed < 0:
        raise ValueError("decode completion elapsed time is invalid")
    if completion.get("contract") != expected_contract:
        raise ValueError("decode completion contract changed")
    if completion.get("contract_sha256") != _canonical_sha256(expected_contract):
        raise ValueError("decode completion contract checksum failed")
    expected_names = {
        "decode_contract.json",
        "decoder_run.json",
        "overall.summary.tsv",
        "posterior.zst",
        "posterior.zst.meta",
        "gamma_unit_scores.tsv",
    }
    outputs = completion.get("outputs")
    if not isinstance(outputs, Mapping) or set(outputs) != expected_names:
        raise ValueError("decode completion output inventory changed")
    scores = pd.read_csv(decoded / "gamma_unit_scores.tsv", sep="\t")
    if len(scores) != 2 * len(TMRCA_THRESHOLDS_YEARS):
        raise ValueError("decode score row count changed")
    if set(scores["window"]) != {"focal", "local_100kb"}:
        raise ValueError("decode score windows changed")
    if set(scores["threshold_years"].astype(int)) != set(TMRCA_THRESHOLDS_YEARS):
        raise ValueError("decode score thresholds changed")
    for name, record in outputs.items():
        rows = len(scores) if name == "gamma_unit_scores.tsv" else None
        _validate_output_record(record, decoded / name, expected_rows=rows)
    return completion


def _decode_unit(payload: Mapping[str, Any]) -> dict[str, Any]:
    from .focused_selection_simulation import UnitLockHeld, exclusive_unit_lock

    root = Path(payload["repo_root"])
    plan = dict(payload["plan"])
    unit = dict(payload["unit"])
    decoder_bin = Path(payload["decoder_bin"])
    threads = int(payload["threads"])
    directory = _unit_dir(root, str(unit["unit_id"]))
    decoded = directory / "decoded"
    if not (directory / "simulation_complete.json").is_file():
        return {
            "unit_id": unit["unit_id"],
            "status": "not_simulated",
            "error": "",
        }
    simulation_contract = json.loads(
        (directory / "simulation_contract.json").read_text(encoding="utf-8")
    )
    try:
        contracted_slim = _resolve_binary_identity(
            simulation_contract["slim_binary"]["identity"], root, "SLiM"
        )
    except (KeyError, TypeError) as error:
        raise ValueError("simulation contract lacks a binary identity") from error
    expected_simulation_contract = _simulation_contract(
        plan, unit, contracted_slim, root
    )
    if simulation_contract != expected_simulation_contract:
        raise ValueError("simulation contract differs from the frozen plan")
    _validate_simulation_completion(directory, expected_simulation_contract)
    contract = _decode_contract(plan, unit, directory, decoder_bin, root, threads)
    decoded.mkdir(parents=True, exist_ok=True)
    try:
        with exclusive_unit_lock(
            decoded / ".decode.lock",
            unit_id=str(unit["unit_id"]),
            phase="han_s001_decode",
        ):
            if (decoded / "decode_complete.json").is_file():
                _validate_decode_completion(directory, contract)
                return {
                    "unit_id": unit["unit_id"],
                    "status": "cached",
                    "error": "",
                }
            _atomic_json(decoded / "decode_contract.json", contract)
            started = perf_counter()
            summary_path = decoded / "overall.summary.tsv"
            raw_path = decoded / "posterior.zst"
            run = run_within_decoder(
                decoder_bin,
                directory / "simulation.trees",
                summary_path,
                scaled_mutation_rate=contract["settings"]["scaled_mutation_rate"],
                recombination_to_mutation_ratio=(
                    contract["settings"]["recombination_to_mutation_ratio"]
                ),
                mutation_rate=MUTATION_RATE,
                threshold_years=TMRCA_THRESHOLDS_YEARS,
                generation_time=GENERATION_TIME_YEARS,
                input_format="trees",
                raw_output=raw_path,
                output_at_stride=DEFAULT_OUTPUT_STRIDE_BP,
                output_at_hets=False,
                only_within=True,
                threads=threads,
                cache_size=DEFAULT_CACHE_SIZE,
                pair_block=256,
                exp10="accurate",
                backward_alignment="fixed",
            )
            _atomic_json(decoded / "decoder_run.json", run)
            summary = pd.read_csv(summary_path, sep="\t")
            scores = _gamma_scores(summary)
            _atomic_frame(decoded / "gamma_unit_scores.tsv", scores)
            output_names = (
                "decode_contract.json",
                "decoder_run.json",
                "overall.summary.tsv",
                "posterior.zst",
                "posterior.zst.meta",
                "gamma_unit_scores.tsv",
            )
            outputs = {
                name: _output_record(
                    decoded / name,
                    rows=(len(scores) if name == "gamma_unit_scores.tsv" else None),
                )
                for name in output_names
            }
            completion = {
                "schema": DECODE_SCHEMA,
                "status": "complete",
                "unit_id": unit["unit_id"],
                "contract": contract,
                "contract_sha256": _canonical_sha256(contract),
                "elapsed_seconds": perf_counter() - started,
                "outputs": outputs,
            }
            _atomic_json(decoded / "decode_complete.json", completion)
            (decoded / "decode_failed.json").unlink(missing_ok=True)
            return {
                "unit_id": unit["unit_id"],
                "status": "complete",
                "error": "",
            }
    except UnitLockHeld as error:
        return {
            "unit_id": unit["unit_id"],
            "status": "locked",
            "error": str(error),
        }
    except Exception as error:
        _atomic_json(
            decoded / "decode_failed.json",
            {
                "schema": DECODE_SCHEMA,
                "status": "failed",
                "unit_id": unit["unit_id"],
                "error_type": type(error).__name__,
                "error": str(error),
                "contract_sha256": _canonical_sha256(contract),
            },
        )
        return {
            "unit_id": unit["unit_id"],
            "status": "failed",
            "error": f"{type(error).__name__}: {error}",
        }


def decode_study(
    repo_root: str | Path,
    decoder_bin: str | Path,
    *,
    workers: int = DEFAULT_WORKERS,
    threads: int = DEFAULT_DECODE_THREADS,
) -> pd.DataFrame:
    """Run or resume Gamma-SMC overall-pair decoding for every unit."""

    root = Path(repo_root).resolve()
    plan = load_plan(root)
    binary = Path(decoder_bin).resolve()
    if not binary.is_file():
        raise ValueError(f"Gamma-SMC binary is absent: {binary}")
    n_workers = _positive_integer(workers, "workers")
    n_threads = _positive_integer(threads, "threads")
    if n_threads != DEFAULT_DECODE_THREADS:
        raise ValueError(
            f"this study requires exactly {DEFAULT_DECODE_THREADS} decode thread"
        )
    payloads = [
        {
            "repo_root": str(root),
            "plan": plan,
            "unit": unit,
            "decoder_bin": str(binary),
            "threads": n_threads,
        }
        for unit in plan["units"]
    ]
    if n_workers == 1:
        results = [_decode_unit(payload) for payload in payloads]
    else:
        results = []
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = [executor.submit(_decode_unit, payload) for payload in payloads]
            for future in as_completed(futures):
                results.append(future.result())
    order = {unit["unit_id"]: index for index, unit in enumerate(plan["units"])}
    frame = pd.DataFrame(results).sort_values(
        "unit_id", key=lambda values: values.map(order)
    )
    frame = frame.reset_index(drop=True)
    _atomic_frame(root / DECODE_STATUS_RELATIVE_PATH, frame)
    if set(frame["status"]).issubset({"complete", "cached"}):
        collect_replicate_scores(root)
    return frame


def decode_status(repo_root: str | Path) -> pd.DataFrame:
    """Return read-only checksum-validated per-unit decode status."""

    root = Path(repo_root).resolve()
    plan = load_plan(root)
    rows: list[dict[str, Any]] = []
    for unit in plan["units"]:
        directory = _unit_dir(root, str(unit["unit_id"]))
        decoded = directory / "decoded"
        contract_path = decoded / "decode_contract.json"
        if (decoded / "decode_complete.json").is_file() and contract_path.is_file():
            try:
                _validated_unit_artifacts(root, plan, unit)
                status = "complete"
                error = ""
            except (
                KeyError,
                OSError,
                TypeError,
                ValueError,
                json.JSONDecodeError,
            ) as validation_error:
                status = "invalid"
                error = str(validation_error)
        elif (decoded / "decode_failed.json").is_file():
            status = "failed"
            failure = json.loads(
                (decoded / "decode_failed.json").read_text(encoding="utf-8")
            )
            error = str(failure.get("error", ""))
        elif not (directory / "simulation_complete.json").is_file():
            status = "not_simulated"
            error = ""
        elif decoded.exists():
            status = "partial"
            error = ""
        else:
            status = "not_started"
            error = ""
        rows.append({**unit, "status": status, "error": error})
    return pd.DataFrame(rows)


def _truth_scores(truth: pd.DataFrame) -> pd.DataFrame:
    required = {
        "position_0based",
        "threshold_years",
        "p_tmrca_lt_threshold",
        "genotype_class",
    }
    if not required.issubset(truth.columns) or set(truth["genotype_class"]) != {
        "overall"
    }:
        raise ValueError("tree-truth profile schema or overall class changed")
    _validate_exact_profile_grid(truth, "tree-truth")
    _validate_profile_estimand(truth, "tree-truth")
    positions = truth["position_0based"].to_numpy(dtype=float)
    focal = np.isclose(positions, FOCAL_POSITION_BP, rtol=0.0, atol=1e-9)
    if focal.sum() != len(TMRCA_THRESHOLDS_YEARS):
        raise ValueError("tree truth does not contain one focal row per threshold")
    local = (positions >= FOCAL_POSITION_BP - LOCAL_WINDOW_HALF_WIDTH_BP) & (
        positions <= FOCAL_POSITION_BP + LOCAL_WINDOW_HALF_WIDTH_BP
    )
    rows: list[dict[str, Any]] = []
    for window, selected in (("focal", focal), ("local_100kb", local)):
        frame = truth.loc[selected]
        for threshold in TMRCA_THRESHOLDS_YEARS:
            values = frame.loc[
                np.isclose(
                    frame["threshold_years"].to_numpy(dtype=float),
                    threshold,
                    rtol=0.0,
                    atol=1e-9,
                ),
                "p_tmrca_lt_threshold",
            ].to_numpy(dtype=float)
            if not len(values) or not np.all(np.isfinite(values)):
                raise ValueError("tree-truth score window has missing values")
            expected_count = 1 if window == "focal" else 21
            if len(values) != expected_count:
                raise ValueError(
                    f"tree-truth {window} window must contain "
                    f"{expected_count} positions"
                )
            score = float(values.mean())
            if not 0.0 <= score <= 1.0:
                raise ValueError("tree-truth probability lies outside [0, 1]")
            rows.append(
                {
                    "source": "tree_truth",
                    "window": window,
                    "threshold_years": int(threshold),
                    "score": score,
                }
            )
    return pd.DataFrame(rows)


def _validated_unit_artifacts(
    root: Path, plan: Mapping[str, Any], unit: Mapping[str, Any]
) -> tuple[dict[str, Any], dict[str, Any], Path]:
    directory = _unit_dir(root, str(unit["unit_id"]))
    try:
        simulation_contract = json.loads(
            (directory / "simulation_contract.json").read_text(encoding="utf-8")
        )
        decode_contract = json.loads(
            (directory / "decoded/decode_contract.json").read_text(encoding="utf-8")
        )
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"unit {unit['unit_id']} contracts are unreadable") from error
    try:
        slim_identity = simulation_contract["slim_binary"]["identity"]
    except (KeyError, TypeError) as error:
        raise ValueError("simulation contract lacks a binary identity") from error
    slim_bin = _resolve_binary_identity(slim_identity, root, "SLiM")
    expected_simulation_contract = _simulation_contract(plan, unit, slim_bin, root)
    if simulation_contract != expected_simulation_contract:
        raise ValueError(f"unit {unit['unit_id']} simulation contract is stale")
    simulation_completion = _validate_simulation_completion(
        directory, expected_simulation_contract
    )
    try:
        decoder_identity = decode_contract["decoder_binary"]["identity"]
        threads = _positive_integer(
            decode_contract["settings"]["threads"], "contracted decode threads"
        )
    except (KeyError, TypeError) as error:
        raise ValueError("decode contract lacks binary/settings identity") from error
    decoder_bin = _resolve_binary_identity(decoder_identity, root, "Gamma-SMC")
    expected_decode_contract = _decode_contract(
        plan, unit, directory, decoder_bin, root, threads
    )
    if decode_contract != expected_decode_contract:
        raise ValueError(f"unit {unit['unit_id']} decode contract is stale")
    decode_completion = _validate_decode_completion(directory, expected_decode_contract)
    return simulation_completion, decode_completion, directory


def _collect_payload(
    root: Path, plan: Mapping[str, Any]
) -> tuple[pd.DataFrame, dict[str, Any]]:
    rows: list[pd.DataFrame] = []
    simulation_hashes: dict[str, str] = {}
    decode_hashes: dict[str, str] = {}
    shared_slim_binary: dict[str, Any] | None = None
    shared_decoder_binary: dict[str, Any] | None = None
    shared_decode_settings: dict[str, Any] | None = None
    for unit in plan["units"]:
        simulation, decode, directory = _validated_unit_artifacts(root, plan, unit)
        simulation_binary = dict(simulation["contract"]["slim_binary"])
        decoder_binary = dict(decode["contract"]["decoder_binary"])
        decode_settings = dict(decode["contract"]["settings"])
        if shared_slim_binary is None:
            shared_slim_binary = simulation_binary
            shared_decoder_binary = decoder_binary
            shared_decode_settings = decode_settings
        elif (
            simulation_binary != shared_slim_binary
            or decoder_binary != shared_decoder_binary
            or decode_settings != shared_decode_settings
        ):
            raise ValueError("study units mix SLiM/Gamma binaries or decode settings")
        simulation_path = directory / "simulation_complete.json"
        decode_path = directory / "decoded/decode_complete.json"
        simulation_hashes[str(unit["unit_id"])] = _sha256_file(simulation_path)
        decode_hashes[str(unit["unit_id"])] = _sha256_file(decode_path)
        truth = _truth_scores(pd.read_csv(directory / "truth_overall.tsv", sep="\t"))
        gamma = pd.read_csv(directory / "decoded/gamma_unit_scores.tsv", sep="\t")[
            ["source", "window", "threshold_years", "score"]
        ]
        scores = pd.concat([gamma, truth], ignore_index=True)
        metadata = simulation["metadata"]
        scores.insert(0, "seed", int(unit["seed"]))
        scores.insert(0, "selection_coefficient", float(unit["selection_coefficient"]))
        scores.insert(0, "simulation_class", str(unit["simulation_class"]))
        scores.insert(0, "unit_id", str(unit["unit_id"]))
        scores["final_population_af"] = float(metadata["final_population_af"])
        scores["sample_alt_count"] = int(metadata["sample_alt_count"])
        scores["sample_af"] = float(metadata["sample_af"])
        scores["sample_detected"] = bool(metadata["sample_detected"])
        rows.append(scores.loc[:, SCORE_COLUMNS])
        if decode.get("status") != "complete":
            raise ValueError(f"unit {unit['unit_id']} decode is not complete")
    combined = pd.concat(rows, ignore_index=True)
    expected_rows = len(plan["units"]) * 4 * len(TMRCA_THRESHOLDS_YEARS)
    if len(combined) != expected_rows:
        raise ValueError("collected replicate score row count changed")
    from .han_s001_tmrca_analysis import validate_replicate_scores

    combined = validate_replicate_scores(combined)
    provenance = {
        "schema": AGGREGATE_SCHEMA,
        "status": "complete",
        "plan_contract_sha256": plan["contract_sha256"],
        "n_units": len(plan["units"]),
        "n_neutral": int(plan["n_neutral"]),
        "n_selected": int(plan["n_selected"]),
        "score_rows": len(combined),
        "score_columns": list(SCORE_COLUMNS),
        "simulation_completion_sha256s": simulation_hashes,
        "decode_completion_sha256s": decode_hashes,
        "simulation_completion_sha256": _canonical_sha256(simulation_hashes),
        "decode_completion_sha256": _canonical_sha256(decode_hashes),
        "shared_slim_binary": shared_slim_binary,
        "shared_decoder_binary": shared_decoder_binary,
        "shared_decode_settings": shared_decode_settings,
    }
    return combined, provenance


def collect_replicate_scores(repo_root: str | Path) -> pd.DataFrame:
    """Collect exact focal and +/-100 kb truth/Gamma scores for all units."""

    root = Path(repo_root).resolve()
    plan = load_plan(root)
    scores, provenance = _collect_payload(root, plan)
    score_path = root / SCORES_RELATIVE_PATH
    _atomic_frame(score_path, scores)
    completion = {
        **provenance,
        "scores": _output_record(score_path, rows=len(scores)),
    }
    completion["aggregate_contract_sha256"] = _canonical_sha256(completion)
    _atomic_json(root / STUDY_COMPLETION_RELATIVE_PATH, completion)
    return scores


def verify_study(repo_root: str | Path) -> dict[str, Any]:
    """Strictly verify the complete plan, unit artifacts, and aggregate scores."""

    root = Path(repo_root).resolve()
    plan = load_plan(root)
    expected_scores, expected = _collect_payload(root, plan)
    score_path = root / SCORES_RELATIVE_PATH
    completion_path = root / STUDY_COMPLETION_RELATIVE_PATH
    if not score_path.is_file() or not completion_path.is_file():
        raise ValueError("aggregate scores/completion are absent; finish decode first")
    try:
        completion = json.loads(completion_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError("study completion is unreadable") from error
    observed_contract_sha256 = completion.pop("aggregate_contract_sha256", None)
    if observed_contract_sha256 != _canonical_sha256(completion):
        raise ValueError("study aggregate contract checksum failed")
    expected_completion = {
        **expected,
        "scores": _output_record(score_path, rows=len(expected_scores)),
    }
    if completion != expected_completion:
        raise ValueError("study completion differs from verified unit artifacts")
    observed_scores = pd.read_csv(score_path, sep="\t")
    from .han_s001_tmrca_analysis import validate_replicate_scores

    observed_scores = validate_replicate_scores(observed_scores)
    pd.testing.assert_frame_equal(
        observed_scores.reset_index(drop=True),
        expected_scores.reset_index(drop=True),
        check_dtype=False,
        check_exact=False,
        rtol=1e-11,
        atol=1e-12,
    )
    return {
        **completion,
        "aggregate_contract_sha256": observed_contract_sha256,
    }


__all__ = [
    "DEFAULT_ANALYSIS_DIR",
    "DEFAULT_NEUTRAL_REPLICATES",
    "DEFAULT_SELECTED_REPLICATES",
    "DEFAULT_SELECTION_COEFFICIENTS",
    "DEFAULT_WORKERS",
    "PLAN_RELATIVE_PATH",
    "SCORE_COLUMNS",
    "TMRCA_THRESHOLDS_YEARS",
    "collect_replicate_scores",
    "decode_status",
    "decode_study",
    "load_plan",
    "plan_path",
    "simulate_study",
    "simulation_status",
    "verify_study",
    "write_plan",
]
