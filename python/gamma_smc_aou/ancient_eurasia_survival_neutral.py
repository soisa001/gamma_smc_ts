"""Restartable single-origin neutral bank under ``AncientEurasia_9K19``.

The workflow is additive and deliberately narrow: every replicate begins with
the same one-copy Neanderthal mutation used by the focused selected simulation,
follows the Neanderthal -> Loschbour -> Han route, has ``s = 0`` everywhere,
and is conditioned only on strict segregation in the complete present-day Han
census.  The returned 500-diploid observational pool and uniformly sampled
100-diploid panel are never used as acceptance gates.
"""

from __future__ import annotations

import argparse
import hashlib
import inspect
import json
import math
import os
import platform
import sys
from collections.abc import Mapping, Sequence
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import contextmanager
from dataclasses import replace
from datetime import UTC, datetime
from itertools import pairwise
from pathlib import Path
from time import perf_counter
from typing import Any

import msprime
import numpy as np
import pandas as pd
import pyslim
import scipy
import stdpopsim
import tskit
from stdpopsim import slim_engine

from .decoder import run_within_decoder
from .defaults import DEFAULT_CACHE_SIZE
from .eas_sweep_analysis import (
    build_diploid_pair_table,
    gamma_pair_tables,
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
    serialize_extended_events,
)
from .eas_sweep_study import BASE_SEED
from .focused_selection_simulation import (
    UnitLockHeld,
    _selected_focal_expectation,
    build_raw_han_selected_sweep,
    exclusive_unit_lock,
    inspect_selected_focal_identity,
    validate_selected_focal_identity,
)
from .han_s001_tmrca_simulation import (
    _focal_counts_or_zero,
    _present_han_ne,
    _sampled_diploids,
    _uniform_panel,
)

STUDY_ID = "ancient_eurasia_single_origin_survival"
SCENARIO = "ancient_eurasia_single_origin_survival_neutral"
PLAN_SCHEMA = "gamma-smc.ancient-eurasia-survival-neutral-plan/v2"
SIMULATION_SCHEMA = "gamma-smc.ancient-eurasia-survival-neutral-simulation/v2"
DECODE_SCHEMA = "gamma-smc.ancient-eurasia-survival-neutral-decode/v2"
RESULTS_SCHEMA = "gamma-smc.ancient-eurasia-survival-neutral-results/v2"
DECODER_RUN_SCHEMA = "gamma-smc.ancient-eurasia-survival-neutral-decoder-run/v2"
FOCAL_OVERLAY_PATCH_SCHEMA = "gamma-smc.ancient-eurasia-reserved-focal-overlay-patch/v1"
INDIVIDUAL_LOCATION_SCHEMA = (
    "gamma-smc.ancient-eurasia-individual-location-canonicalization/v1"
)

DEFAULT_STUDY_ROOT = Path(
    "focused_selection_EAS_sim/ancient_eurasia_single_origin_survival_v2"
)
INSTRUMENTED_SMOKE_ROOT_BASENAME = (
    "ancient_eurasia_single_origin_survival_v2_instrumented_rep051_smoke"
)
INSTRUMENTED_SMOKE_UNIT_ID = f"{STUDY_ID}__neutral__rep051"
DEFAULT_REPLICATES = 100
DEFAULT_WORKERS = 4
MAX_WORKERS = 24
DEFAULT_DECODE_THREADS = 1
DEFAULT_POOL_DIPLOIDS = 500
DEFAULT_PANEL_DIPLOIDS = 100
DEFAULT_SLIM_SCALING_FACTOR = 5.0
DEFAULT_SLIM_BURN_IN = 0.1
DEFAULT_OUTPUT_STRIDE_BP = 10_000
TMRCA_THRESHOLDS_YEARS = (1_000, 4_500, 10_000, 20_000, 30_000, 40_000, 50_000)
# Frozen from the isolated instrumented rep051 smoke under stdpopsim 0.3.0.
EXPECTED_OVERLAY_CALL_COUNT: int | None = 2
EXPECTED_OVERLAY_PATCHED_COUNT: int | None = 1
MUTATION_ORIGIN_GENERATIONS = 2_400.0
PRE_PULSE_SOURCE_CHECK_GENERATIONS = 2_273.0
REALIZED_HAN_SPLIT_GENERATIONS_Q5 = 2_015.0
REALIZED_INTROGRESSION_PULSE_GENERATIONS_Q5 = 2_270.0
EXPECTED_SLIM_SHA256 = (
    "67409ff190808c2967c949cae6697bf235010dbe832389cb4976d9089b5f5e1d"
)
EXPECTED_DECODER_SHA256 = (
    "057774590af12efbbedc9d1ffb736747c0dba2604ef0028dcf9b65e275e70a22"
)

SOURCE_PATHS = (
    "python/gamma_smc_aou/ancient_eurasia_survival_neutral.py",
    "python/gamma_smc_aou/focused_selection_simulation.py",
    "python/gamma_smc_aou/focused_selection_campaign.py",
    "python/gamma_smc_aou/focused_selection_decode.py",
    "python/gamma_smc_aou/eas_sweep_models.py",
    "python/gamma_smc_aou/eas_sweep_analysis.py",
    "python/gamma_smc_aou/decoder.py",
    "python/gamma_smc_aou/tree_sequence.py",
    "python/gamma_smc_aou/han_s001_tmrca_simulation.py",
)

EXECUTION_COLUMNS = (
    "unit_id",
    "replicate_index",
    "seed",
    "panel_seed",
    "scenario",
    "selection_coefficient",
    "mutation_origin_population",
    "mutation_age_generations",
    "pulse_generations",
    "present_conditioning",
    "pool_frequency_conditioning",
    "terminal_af_target_or_band_conditioning",
    "sample_detection_conditioning",
    "pool_diploids",
    "sample_diploids",
)

SCORE_COLUMNS = (
    "unit_id",
    "replicate_index",
    "seed",
    "source",
    "threshold_years",
    "overall_focal_p_tmrca_lt_threshold",
    "final_pool_alt_count",
    "final_pool_total_count",
    "final_pool_af",
    "sample_alt_count",
    "sample_total_count",
    "sample_af",
    "sample_detected",
)

ENDPOINT_COLUMNS = (
    "unit_id",
    "replicate_index",
    "seed",
    "simulation_completion_sha256",
    "final_census_alt_count",
    "final_census_total_count",
    "final_census_af",
    "final_census_segregating",
    "pool_diploids",
    "pool_alt_count",
    "pool_total_count",
    "pool_af",
    "pool_detected",
    "pool_fixed",
    "sample_diploids",
    "sample_alt_count",
    "sample_total_count",
    "sample_af",
    "sample_detected",
    "sample_fixed",
    "panel_seed",
)

SIMULATION_OUTPUT_PATHS = {
    "simulation_contract": "simulation_contract.json",
    "tree": "simulation.trees",
    "sample_manifest": "sample_manifest.tsv",
    "slim_log": "slim.tsv",
    "truth_profiles": "truth_profiles.tsv.gz",
    "truth_class_summaries": "truth_class_summaries.tsv",
    "truth_unit_scores": "truth_unit_scores.tsv",
}

DECODE_OUTPUT_PATHS = {
    "decode_contract": "decode_contract.json",
    "decoder_run": "decoder_run.json",
    "overall_pairs": "overall.pairs.tsv",
    "overall_summary": "overall.summary.tsv",
    "posterior": "posterior.zst",
    "posterior_metadata": "posterior.zst.meta",
    "gamma_unit_scores": "gamma_unit_scores.tsv",
}

RESULT_OUTPUT_PATHS = {
    "endpoint_summary": "endpoint_summary.tsv",
    "replicate_scores": "replicate_scores.tsv",
    "simulation_status": "simulation_status.tsv",
    "decode_status": "decode_status.tsv",
}

SIMULATION_PARAMETERS = {
    "sequence_length_bp": SEQUENCE_LENGTH_BP,
    "focal_position_0based": FOCAL_POSITION_BP,
    "mutation_rate": MUTATION_RATE,
    "recombination_rate": RECOMBINATION_RATE,
    "slim_scaling_factor": DEFAULT_SLIM_SCALING_FACTOR,
    "slim_burn_in": DEFAULT_SLIM_BURN_IN,
    "pool_diploids": DEFAULT_POOL_DIPLOIDS,
    "sample_diploids": DEFAULT_PANEL_DIPLOIDS,
    "generation_time_years": GENERATION_TIME_YEARS,
    "output_stride_bp": DEFAULT_OUTPUT_STRIDE_BP,
    "present_han_unscaled_ne": 6_300.0,
    "present_han_scaled_census_diploids_q5": 1_260,
    "present_han_scaled_census_genomes_q5": 2_520,
    "census_semantics": "complete_simulated_q5_present_han_census",
}

DECODER_SETTINGS = {
    "input_format": "trees_via_streamed_vcf",
    "pair_selector": "explicit_within_individual_pairs",
    "n_pairs": DEFAULT_PANEL_DIPLOIDS,
    "scaled_mutation_rate": 0.000315,
    "unscaled_mutation_rate": MUTATION_RATE,
    "recombination_to_mutation_ratio": RECOMBINATION_RATE / MUTATION_RATE,
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
    "threads": DEFAULT_DECODE_THREADS,
}

ESTIMAND = {
    "source_model": "stdpopsim AncientEurasia_9K19",
    "mutation_origin_population": ARCHAIC_POPULATION,
    "mutation_origin_mode": "single_copy",
    "mutation_age_generations": MUTATION_ORIGIN_GENERATIONS,
    "pulse_generations": INTROGRESSION_PULSE_GENERATIONS,
    "selection_coefficient": 0.0,
    "present_conditioning": "strict_segregating_in_present_han_census",
    "pool_frequency_conditioning": False,
    "terminal_af_target_or_band_conditioning": False,
    "sample_detection_conditioning": False,
    "panel_sampling": "uniform_100_of_500",
    "sample_diploids": DEFAULT_PANEL_DIPLOIDS,
    "candidate_pool_diploids": DEFAULT_POOL_DIPLOIDS,
    "thresholds_years": list(TMRCA_THRESHOLDS_YEARS),
}

SIMULATION_LAW = {
    "engine": "stdpopsim_slim_forward",
    "demographic_model": "AncientEurasia_9K19",
    "sequence_scope": "direct_full_10mb",
    "mutation_introduction": "single_copy_drawmutation_neanderthal_at_2400_generations",
    "conditioning": "internal_slim_extended_event_rejection",
    "branch_enumeration": False,
    "ancestry_scan": False,
    "branch_span_weighting": False,
    "cropping": False,
    "posthoc_causal_mutation_placement": False,
    "background_mutation_overlay": "stdpopsim_post_slim_neutral_msprime_overlay",
    "background_overlay_reserved_interval": (
        "every_layer_masks_0based_half_open_[5000000,5000001)"
    ),
    "reserved_focal_base_overlay_guard": (
        "scoped_msprime_sim_mutations_rate_map_mask_during_engine_simulate"
    ),
    "nonspatial_individual_location_canonicalization": (
        "post_uniform_panel_pre_analysis_exact_500x3_float64_zeros"
    ),
    "external_pool_frequency_rejection": False,
    "external_sample_detection_rejection": False,
    "external_terminal_target_af_rejection": False,
}


def _expected_overlay_inventory() -> list[dict[str, Any]] | None:
    if EXPECTED_OVERLAY_CALL_COUNT is None or EXPECTED_OVERLAY_PATCHED_COUNT is None:
        return None
    focal = float(FOCAL_POSITION_BP)
    length = float(SEQUENCE_LENGTH_BP)
    mu = float(MUTATION_RATE)
    return [
        {
            "layer": "ordinary_neutral_dfe",
            "model_class": "msprime.mutations.SLiMMutationModel",
            "model_type": 0,
            "start_time": None,
            "end_time": "finite_shared_slim_tick",
            "original_positions": [0.0, focal, focal + 1.0, length],
            "original_rates": [mu, 0.0, mu],
            "passed_positions": [0.0, focal, focal + 1.0, length],
            "passed_rates": [mu, 0.0, mu],
            "patched": False,
        },
        {
            "layer": "recapitation",
            "model_class": "msprime.mutations.SLiMMutationModel",
            "model_type": 2,
            "start_time": "finite_shared_slim_tick",
            "end_time": None,
            "original_positions": [0.0, length],
            "original_rates": [mu],
            "passed_positions": [0.0, focal, focal + 1.0, length],
            "passed_rates": [mu, 0.0, mu],
            "patched": True,
        },
    ]


def _focal_overlay_patch_contract() -> dict[str, Any]:
    """Return the immutable contract for the v2 reserved-base guard."""

    original = msprime.sim_mutations
    try:
        source = inspect.getsource(original)
    except (OSError, TypeError) as error:
        raise RuntimeError("msprime.sim_mutations source is unavailable") from error
    engine_method = type(stdpopsim.get_engine("slim"))._recap_and_rescale
    try:
        engine_source = inspect.getsource(engine_method)
    except (OSError, TypeError) as error:
        raise RuntimeError(
            "SlimEngine._recap_and_rescale source is unavailable"
        ) from error
    return {
        "schema": FOCAL_OVERLAY_PATCH_SCHEMA,
        "scope": "only_during_one_stdpopsim_slim_engine_simulate_call",
        "patched_callable": f"{original.__module__}.{original.__qualname__}",
        "patched_callable_signature": str(inspect.signature(original)),
        "patched_callable_source_sha256": hashlib.sha256(
            source.encode("utf-8")
        ).hexdigest(),
        "upstream_overlay_method": (
            f"{engine_method.__module__}.{engine_method.__qualname__}"
        ),
        "upstream_overlay_method_signature": str(inspect.signature(engine_method)),
        "upstream_overlay_method_source_sha256": hashlib.sha256(
            engine_source.encode("utf-8")
        ).hexdigest(),
        "target_interval_0based_half_open": [
            FOCAL_POSITION_BP,
            FOCAL_POSITION_BP + 1,
        ],
        "detection": (
            "every intercepted msprime.RateMap with positive integrated rate "
            "on the target interval"
        ),
        "replacement": (
            "insert target breakpoints and set only the target interval rate to zero"
        ),
        "non_rate_arguments": "preserved_by_identity",
        "flank_rates_and_integrals": "exactly_preserved",
        "minimum_patched_call_count": 1,
        "expected_call_count": EXPECTED_OVERLAY_CALL_COUNT,
        "expected_patched_call_count": EXPECTED_OVERLAY_PATCHED_COUNT,
        "expected_call_inventory": _expected_overlay_inventory(),
        "inventory_freeze_state": (
            "frozen_after_rep051_smoke"
            if EXPECTED_OVERLAY_CALL_COUNT is not None
            and EXPECTED_OVERLAY_PATCHED_COUNT is not None
            else "pending_instrumented_rep051_smoke"
        ),
        "restoration": "original_callable_restored_on_success_and_exception",
    }


def _individual_location_contract() -> dict[str, Any]:
    zeros = np.zeros(3 * DEFAULT_POOL_DIPLOIDS, dtype=np.float64)
    offsets = np.arange(0, 3 * DEFAULT_POOL_DIPLOIDS + 1, 3, dtype=np.uint32)
    method = tskit.IndividualTable.packset_location
    try:
        source = inspect.getsource(method)
    except (OSError, TypeError) as error:
        raise RuntimeError(
            "tskit IndividualTable.packset_location source is unavailable"
        ) from error
    return {
        "schema": INDIVIDUAL_LOCATION_SCHEMA,
        "scope": "immediately_after_uniform_panel_before_analysis_or_persistence",
        "expected_individual_count": DEFAULT_POOL_DIPLOIDS,
        "required_location_width_per_individual": 3,
        "precanonicalization_values": "arbitrary_float64_bits_semantically_discarded",
        "canonical_location": [0.0, 0.0, 0.0],
        "canonical_flat_location_record": _array_record(zeros),
        "canonical_location_offset_record": _array_record(offsets),
        "mutation": "IndividualTable.packset_location_only",
        "mutation_api": f"{method.__module__}.{method.__qualname__}",
        "mutation_api_signature": str(inspect.signature(method)),
        "mutation_api_source_sha256": hashlib.sha256(
            source.encode("utf-8")
        ).hexdigest(),
        "preserve_location_offsets_and_widths": True,
        "preserve_all_nonlocation_individual_columns": True,
        "preserve_all_other_tables_ts_metadata_and_reference": True,
    }


def _selected_template_unit() -> dict[str, Any]:
    """Return the minimal valid selected-unit record used only as an event template."""

    return {
        "unit_id": "ancient_eurasia_single_origin_survival__event_template",
        "demography_id": "ancient_eurasia_han_introgression",
        "demography_kind": "introgression",
        "population": INTROGRESSION_TARGET_POPULATION,
        "simulation_class": "selected",
        "selection_coefficient": 0.01,
        "target_allele_frequency": 0.20,
        "population_af_lower": 0.175,
        "population_af_upper": 0.225,
        "exact_sample_alt_count": 40,
        "sample_diploids": DEFAULT_PANEL_DIPLOIDS,
        "replicate_index": 0,
        "seed": 1,
        "sequence_length_bp": SEQUENCE_LENGTH_BP,
        "focal_position_bp": FOCAL_POSITION_BP,
    }


def build_neutral_event_schedule() -> tuple[tuple[Any, ...], dict[str, Any]]:
    """Derive and validate the exact six-event neutral schedule.

    The transformation intentionally starts from the raw focused selected
    constructor.  It removes exactly its two fitness events and exactly its two
    terminal Han AF-band events, then adds the strict present-day ``AF < 1``
    condition.  The retained Han ``AF > 0`` condition plus this added event is
    strict census segregation; neither the pool nor panel participates.
    """

    sweep = build_raw_han_selected_sweep(
        _selected_template_unit(), slim_scaling_factor=DEFAULT_SLIM_SCALING_FACTOR
    )
    retained: list[Any] = []
    removed_fitness: list[Any] = []
    removed_terminal: list[Any] = []
    for event in sweep.extended_events:
        if isinstance(event, stdpopsim.ChangeMutationFitness):
            removed_fitness.append(event)
            continue
        terminal_han_band = (
            isinstance(event, stdpopsim.ConditionOnAlleleFrequency)
            and str(event.population) == INTROGRESSION_TARGET_POPULATION
            and math.isclose(float(event.start_time), 0.0, abs_tol=1e-12)
            and math.isclose(float(event.end_time), 0.0, abs_tol=1e-12)
            and str(event.op) in {">=", "<="}
        )
        if terminal_han_band:
            removed_terminal.append(event)
            continue
        retained.append(event)
    if len(removed_fitness) != 2:
        raise RuntimeError(
            f"raw selected template must contain exactly two fitness events; found {len(removed_fitness)}"
        )
    if len(removed_terminal) != 2 or {str(event.op) for event in removed_terminal} != {
        ">=",
        "<=",
    }:
        raise RuntimeError(
            "raw selected template must contain exactly two terminal Han AF-band events"
        )
    retained.append(
        stdpopsim.ConditionOnAlleleFrequency(
            start_time=0.0,
            end_time=0.0,
            single_site_id=FOCAL_SITE_ID,
            population=INTROGRESSION_TARGET_POPULATION,
            op="<",
            allele_frequency=1.0,
        )
    )
    serialized = serialize_extended_events(retained)
    exact = "exact_generation"
    after = "generation_after_in_forward_time"
    expected_schedule = [
        {
            "event_type": "DrawMutation",
            "single_site_id": FOCAL_SITE_ID,
            "population": ARCHAIC_POPULATION,
            "time": {"generations_ago": 2_400.0, "semantics": exact},
        },
        {
            "event_type": "ConditionOnAlleleFrequency",
            "single_site_id": FOCAL_SITE_ID,
            "population": ARCHAIC_POPULATION,
            "start_time": {"generations_ago": 2_400.0, "semantics": after},
            "end_time": {"generations_ago": 2_273.0, "semantics": exact},
            "operator": ">",
            "allele_frequency": 0.0,
        },
        {
            "event_type": "ConditionOnAlleleFrequency",
            "single_site_id": FOCAL_SITE_ID,
            "population": ARCHAIC_POPULATION,
            "start_time": {"generations_ago": 2_273.0, "semantics": exact},
            "end_time": {"generations_ago": 2_273.0, "semantics": exact},
            "operator": ">=",
            "allele_frequency": 1e-9,
        },
        {
            "event_type": "ConditionOnAlleleFrequency",
            "single_site_id": FOCAL_SITE_ID,
            "population": INTROGRESSION_RECIPIENT_POPULATION,
            "start_time": {"generations_ago": 2_272.0, "semantics": after},
            "end_time": {"generations_ago": 2_015.0, "semantics": exact},
            "operator": ">",
            "allele_frequency": 0.0,
        },
        {
            "event_type": "ConditionOnAlleleFrequency",
            "single_site_id": FOCAL_SITE_ID,
            "population": INTROGRESSION_TARGET_POPULATION,
            "start_time": {"generations_ago": 2_015.0, "semantics": after},
            "end_time": {"generations_ago": 0.0, "semantics": exact},
            "operator": ">",
            "allele_frequency": 0.0,
        },
        {
            "event_type": "ConditionOnAlleleFrequency",
            "single_site_id": FOCAL_SITE_ID,
            "population": INTROGRESSION_TARGET_POPULATION,
            "start_time": {"generations_ago": 0.0, "semantics": exact},
            "end_time": {"generations_ago": 0.0, "semantics": exact},
            "operator": "<",
            "allele_frequency": 1.0,
        },
    ]
    if len(retained) != 6 or serialized != expected_schedule:
        raise RuntimeError(
            "neutral transformation did not yield the exact six-event ancestry schedule"
        )
    transform = {
        "source_constructor": "focused_selection_simulation.build_raw_han_selected_sweep",
        "source_constructor_slim_scaling_factor": DEFAULT_SLIM_SCALING_FACTOR,
        "removed_change_mutation_fitness_count": len(removed_fitness),
        "removed_terminal_han_af_band_count": len(removed_terminal),
        "added_terminal_han_strict_upper_count": 1,
        "selection_coefficient_after_transform": 0.0,
        "nominal_pulse_generations_ago": INTROGRESSION_PULSE_GENERATIONS,
        "realized_pulse_generations_ago_q5": REALIZED_INTROGRESSION_PULSE_GENERATIONS_Q5,
        "nominal_han_split_generations_ago": HAN_SPLIT_GENERATIONS,
        "realized_han_split_generations_ago_q5": REALIZED_HAN_SPLIT_GENERATIONS_Q5,
        "event_schedule": serialized,
    }
    return tuple(retained), transform


_ORIGINAL_END = """// Output tree sequence file and end the simulation.
function (void)end(void) {
    sim.treeSeqOutput(trees_file, metadata=metadata);
    sim.simulationFinished();
}"""
_CENSUS_METADATA_KEYS = {
    "survival_neutral_final_census_alt_count",
    "survival_neutral_final_census_total_count",
    "survival_neutral_final_census_af",
    "survival_neutral_final_census_segregating",
    "survival_neutral_final_census_fixed",
}


def _sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_sha256(value: Any) -> str:
    payload = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _canonical_source_sha256(path: str | Path) -> str:
    text = Path(path).read_text(encoding="utf-8")
    return hashlib.sha256(
        text.replace("\r\n", "\n").replace("\r", "\n").encode()
    ).hexdigest()


def _atomic_json(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        temporary.write_text(
            json.dumps(
                value, indent=2, sort_keys=True, ensure_ascii=True, allow_nan=False
            )
            + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_frame(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    compression: str | dict[str, Any] | None = None
    if path.name.endswith(".gz"):
        compression = {"method": "gzip", "mtime": 0}
    try:
        frame.to_csv(
            temporary,
            sep="\t",
            index=False,
            lineterminator="\n",
            float_format="%.12g",
            compression=compression,
        )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_tree(path: Path, tree: tskit.TreeSequence) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        tree.dump(temporary)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _output_record(path: Path, *, rows: int | None = None) -> dict[str, Any]:
    if not path.is_file() or path.stat().st_size <= 0:
        raise ValueError(f"output is absent or empty: {path}")
    return {
        "path": path.name,
        "sha256": _sha256_file(path),
        "size_bytes": path.stat().st_size,
        "rows": rows,
    }


def _headered_text_row_count(path: Path) -> int:
    lines = [
        line
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith("#")
    ]
    return max(0, len(lines) - 1)


def _stable_seed(label: str) -> int:
    token = f"{BASE_SEED}:{STUDY_ID}:{label}".encode()
    return 1 + int.from_bytes(hashlib.sha256(token).digest()[:8], "big") % (2**31 - 2)


def _validate_workers(workers: int) -> int:
    if (
        isinstance(workers, bool)
        or not isinstance(workers, int)
        or not 1 <= workers <= MAX_WORKERS
    ):
        raise ValueError(f"workers must be an integer in [1, {MAX_WORKERS}]")
    return workers


def _path_identity(path: Path, repo_root: Path) -> str:
    try:
        return path.resolve().relative_to(repo_root.resolve()).as_posix()
    except ValueError as error:
        raise ValueError("contracted binaries must be beneath repo_root") from error


def _binary_record(
    path: str | Path, repo_root: Path, *, expected_sha256: str, label: str
) -> dict[str, Any]:
    binary = Path(path).resolve()
    if not binary.is_file():
        raise ValueError(f"{label} binary is absent: {binary}")
    observed = _sha256_file(binary)
    if observed != expected_sha256:
        raise ValueError(
            f"{label} binary SHA-256 differs: expected {expected_sha256}, observed {observed}"
        )
    return {
        "identity": _path_identity(binary, repo_root),
        "sha256": observed,
        "size_bytes": binary.stat().st_size,
    }


def _source_hashes(repo_root: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for identity in SOURCE_PATHS:
        path = repo_root / identity
        if not path.is_file():
            raise ValueError(f"implementation source is absent: {identity}")
        result[identity] = _canonical_source_sha256(path)
    return result


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


def _array_record(values: np.ndarray) -> dict[str, Any]:
    array = np.ascontiguousarray(values)
    return {
        "dtype": array.dtype.str,
        "shape": list(array.shape),
        "sha256": hashlib.sha256(array.tobytes(order="C")).hexdigest(),
    }


def _binary_structure_sha256(value: Any) -> str:
    """Hash nested tskit ``asdict`` data without expanding large arrays."""

    digest = hashlib.sha256()

    def add(blob: bytes) -> None:
        digest.update(len(blob).to_bytes(8, "big"))
        digest.update(blob)

    def visit(candidate: Any) -> None:
        if isinstance(candidate, Mapping):
            add(b"mapping")
            for key in sorted(candidate, key=str):
                add(str(key).encode("utf-8"))
                visit(candidate[key])
        elif isinstance(candidate, np.ndarray):
            add(b"ndarray")
            add(candidate.dtype.str.encode("ascii"))
            add(json.dumps(list(candidate.shape)).encode("ascii"))
            add(np.ascontiguousarray(candidate).tobytes(order="C"))
        elif isinstance(candidate, bytes):
            add(b"bytes")
            add(candidate)
        elif isinstance(candidate, (list, tuple)):
            add(type(candidate).__name__.encode("ascii"))
            for item in candidate:
                visit(item)
        elif isinstance(candidate, np.generic):
            visit(candidate.item())
        elif candidate is None or isinstance(candidate, (str, int, float, bool)):
            add(type(candidate).__name__.encode("ascii"))
            add(
                json.dumps(
                    candidate, sort_keys=True, separators=(",", ":"), allow_nan=False
                ).encode("utf-8")
            )
        else:
            raise TypeError(
                "unsupported value in tskit table digest: "
                f"{type(candidate).__module__}.{type(candidate).__qualname__}"
            )

    visit(value)
    return digest.hexdigest()


def _individual_location_state(tables: tskit.TableCollection) -> dict[str, Any]:
    payload = tables.asdict()
    individual_payload = dict(payload["individuals"])
    location = np.asarray(individual_payload.pop("location"))
    offsets = np.asarray(individual_payload["location_offset"])
    widths = np.diff(offsets)
    other_payload = dict(payload)
    other_payload.pop("individuals")
    without_location = dict(payload)
    without_location["individuals"] = individual_payload
    return {
        "individual_count": int(tables.individuals.num_rows),
        "location_value_count": int(location.size),
        "unique_location_widths": sorted({int(value) for value in widths}),
        "location_all_finite": bool(np.all(np.isfinite(location))),
        "location_finite_count": int(np.count_nonzero(np.isfinite(location))),
        "location_nan_count": int(np.count_nonzero(np.isnan(location))),
        "location_positive_infinity_count": int(
            np.count_nonzero(np.isposinf(location))
        ),
        "location_negative_infinity_count": int(
            np.count_nonzero(np.isneginf(location))
        ),
        "location_all_zero": bool(np.all(location == 0.0)),
        "location_nonzero_count": int(np.count_nonzero(location)),
        "location_record": _array_record(location),
        "location_offset_record": _array_record(offsets),
        "individual_nonlocation_sha256": _binary_structure_sha256(individual_payload),
        "other_tables_ts_metadata_reference_sha256": _binary_structure_sha256(
            other_payload
        ),
        "all_except_individual_location_sha256": _binary_structure_sha256(
            without_location
        ),
    }


def _validate_array_record(
    record: Any,
    *,
    expected_dtype: str | None = None,
    expected_shape: Sequence[int] | None = None,
) -> None:
    if not isinstance(record, Mapping) or set(record) != {"dtype", "shape", "sha256"}:
        raise ValueError("array record fields differ")
    try:
        int(str(record["sha256"]), 16)
    except (TypeError, ValueError) as error:
        raise ValueError("array record SHA-256 is malformed") from error
    if len(str(record["sha256"])) != 64:
        raise ValueError("array record SHA-256 is malformed")
    if expected_dtype is not None and record["dtype"] != expected_dtype:
        raise ValueError("array record dtype differs")
    if expected_shape is not None and record["shape"] != list(expected_shape):
        raise ValueError("array record shape differs")


def _validate_individual_location_receipt(
    receipt: Mapping[str, Any],
    contract: Mapping[str, Any],
    *,
    tree: tskit.TreeSequence | None = None,
) -> None:
    if (
        set(receipt)
        != {
            "schema",
            "contract",
            "contract_sha256",
            "before",
            "after",
            "stored_tree",
            "proof",
        }
        or receipt.get("schema") != INDIVIDUAL_LOCATION_SCHEMA
    ):
        raise ValueError("individual-location receipt fields/schema differ")
    if receipt.get("contract") != contract or receipt.get(
        "contract_sha256"
    ) != _canonical_sha256(contract):
        raise ValueError("individual-location canonicalization contract differs")
    before = receipt.get("before")
    after = receipt.get("after")
    stored_tree = receipt.get("stored_tree")
    state_fields = {
        "individual_count",
        "location_value_count",
        "unique_location_widths",
        "location_all_finite",
        "location_finite_count",
        "location_nan_count",
        "location_positive_infinity_count",
        "location_negative_infinity_count",
        "location_all_zero",
        "location_nonzero_count",
        "location_record",
        "location_offset_record",
        "individual_nonlocation_sha256",
        "other_tables_ts_metadata_reference_sha256",
        "all_except_individual_location_sha256",
    }
    if (
        not isinstance(before, Mapping)
        or not isinstance(after, Mapping)
        or not isinstance(stored_tree, Mapping)
        or set(before) != state_fields
        or set(after) != state_fields
        or set(stored_tree) != state_fields
    ):
        raise ValueError("individual-location before/after state fields differ")
    expected_count = int(contract["expected_individual_count"])
    expected_width = int(contract["required_location_width_per_individual"])
    location_dtype = str(contract["canonical_flat_location_record"]["dtype"])
    location_shape = contract["canonical_flat_location_record"]["shape"]
    offset_dtype = str(contract["canonical_location_offset_record"]["dtype"])
    offset_shape = contract["canonical_location_offset_record"]["shape"]
    for state in (before, after, stored_tree):
        if (
            state["individual_count"] != expected_count
            or state["location_value_count"] != expected_count * expected_width
            or state["unique_location_widths"] != [expected_width]
            or state["location_finite_count"]
            + state["location_nan_count"]
            + state["location_positive_infinity_count"]
            + state["location_negative_infinity_count"]
            != state["location_value_count"]
            or bool(state["location_all_finite"])
            != (state["location_finite_count"] == state["location_value_count"])
        ):
            raise ValueError("individual-location dimensions differ")
        for key in (
            "individual_nonlocation_sha256",
            "other_tables_ts_metadata_reference_sha256",
            "all_except_individual_location_sha256",
        ):
            if not isinstance(state[key], str) or len(state[key]) != 64:
                raise ValueError("individual-location preservation hash is malformed")
        _validate_array_record(
            state["location_record"],
            expected_dtype=location_dtype,
            expected_shape=location_shape,
        )
        _validate_array_record(
            state["location_offset_record"],
            expected_dtype=offset_dtype,
            expected_shape=offset_shape,
        )
    if (
        before["location_offset_record"] != contract["canonical_location_offset_record"]
        or after["location_offset_record"]
        != contract["canonical_location_offset_record"]
        or after["location_record"] != contract["canonical_flat_location_record"]
        or not after["location_all_finite"]
        or after["location_finite_count"] != expected_count * expected_width
        or after["location_nan_count"] != 0
        or after["location_positive_infinity_count"] != 0
        or after["location_negative_infinity_count"] != 0
        or not after["location_all_zero"]
        or after["location_nonzero_count"] != 0
        or stored_tree["location_record"] != after["location_record"]
        or stored_tree["location_offset_record"] != after["location_offset_record"]
        or stored_tree["individual_nonlocation_sha256"]
        != after["individual_nonlocation_sha256"]
        or not stored_tree["location_all_finite"]
        or stored_tree["location_finite_count"] != expected_count * expected_width
        or stored_tree["location_nan_count"] != 0
        or stored_tree["location_positive_infinity_count"] != 0
        or stored_tree["location_negative_infinity_count"] != 0
        or not stored_tree["location_all_zero"]
        or stored_tree["location_nonzero_count"] != 0
    ):
        raise ValueError("canonical individual locations are not exact zeros")
    unchanged_keys = (
        "location_offset_record",
        "individual_nonlocation_sha256",
        "other_tables_ts_metadata_reference_sha256",
        "all_except_individual_location_sha256",
    )
    if any(before[key] != after[key] for key in unchanged_keys):
        raise ValueError("individual-location canonicalization changed other data")
    if dict(after) != dict(stored_tree):
        raise ValueError("canonical after state differs from stored tree binding")
    expected_proof = {
        "location_offsets_and_widths_unchanged": True,
        "all_nonlocation_individual_columns_unchanged": True,
        "all_other_tables_ts_metadata_reference_unchanged": True,
        "only_individual_location_values_canonicalized": True,
        "stored_tree_state_bound_at_canonicalization": True,
    }
    if receipt.get("proof") != expected_proof:
        raise ValueError("individual-location preservation proof differs")
    if tree is not None:
        observed = _individual_location_state(tree.dump_tables())
        if observed != stored_tree:
            raise ValueError("stored tree state differs from location receipt binding")


def _canonicalize_individual_locations(
    tree: tskit.TreeSequence,
) -> tuple[tskit.TreeSequence, dict[str, Any]]:
    """Replace nondeterministic nonspatial SLiM locations with exact zeros."""

    contract = _individual_location_contract()
    tables = tree.dump_tables()
    before = _individual_location_state(tables)
    if (
        before["individual_count"] != DEFAULT_POOL_DIPLOIDS
        or before["unique_location_widths"] != [3]
        or before["location_value_count"] != 3 * DEFAULT_POOL_DIPLOIDS
    ):
        raise ValueError("engine-returned individual location layout differs")
    tables.individuals.packset_location(
        [np.zeros(3, dtype=np.float64) for _ in range(DEFAULT_POOL_DIPLOIDS)]
    )
    after = _individual_location_state(tables)
    result = tables.tree_sequence()
    roundtrip = _individual_location_state(result.dump_tables())
    if after != roundtrip:
        raise RuntimeError("individual-location canonicalization changed on round trip")
    proof = {
        "location_offsets_and_widths_unchanged": (
            before["location_offset_record"] == after["location_offset_record"]
            and before["unique_location_widths"] == after["unique_location_widths"]
        ),
        "all_nonlocation_individual_columns_unchanged": (
            before["individual_nonlocation_sha256"]
            == after["individual_nonlocation_sha256"]
        ),
        "all_other_tables_ts_metadata_reference_unchanged": (
            before["other_tables_ts_metadata_reference_sha256"]
            == after["other_tables_ts_metadata_reference_sha256"]
        ),
        "only_individual_location_values_canonicalized": (
            before["all_except_individual_location_sha256"]
            == after["all_except_individual_location_sha256"]
        ),
        "stored_tree_state_bound_at_canonicalization": True,
    }
    receipt = {
        "schema": INDIVIDUAL_LOCATION_SCHEMA,
        "contract": contract,
        "contract_sha256": _canonical_sha256(contract),
        "before": before,
        "after": after,
        "stored_tree": roundtrip,
        "proof": proof,
    }
    _validate_individual_location_receipt(receipt, contract, tree=result)
    return result, receipt


def _rate_map_mass(
    positions: np.ndarray, rates: np.ndarray, start: float, end: float
) -> float:
    overlaps = np.minimum(positions[1:], end) - np.maximum(positions[:-1], start)
    keep = overlaps > 0
    if not np.any(keep):
        return 0.0
    return float(np.sum(rates[keep] * overlaps[keep]))


def _rate_map_arrays(rate_map: msprime.RateMap) -> tuple[np.ndarray, np.ndarray]:
    positions = np.asarray(rate_map.position, dtype=float)
    rates = np.asarray(rate_map.rate, dtype=float)
    if (
        len(positions) != len(rates) + 1
        or len(rates) == 0
        or not np.all(np.isfinite(positions))
        or not np.all(np.isfinite(rates))
        or np.any(rates < 0)
        or np.any(np.diff(positions) <= 0)
    ):
        raise ValueError("intercepted mutation RateMap is not finite and nonnegative")
    return positions, rates


def _rate_map_metrics(rate_map: msprime.RateMap) -> dict[str, Any]:
    positions, rates = _rate_map_arrays(rate_map)
    focal_start = float(FOCAL_POSITION_BP)
    focal_end = focal_start + 1.0
    return {
        "positions": positions.tolist(),
        "rates": rates.tolist(),
        "sequence_length": float(rate_map.sequence_length),
        "total_integral": _rate_map_mass(
            positions, rates, float(positions[0]), float(positions[-1])
        ),
        "left_flank_integral": _rate_map_mass(
            positions, rates, float(positions[0]), focal_start
        ),
        "focal_integral": _rate_map_mass(positions, rates, focal_start, focal_end),
        "right_flank_integral": _rate_map_mass(
            positions, rates, focal_end, float(positions[-1])
        ),
    }


def _mask_focal_rate_map(
    rate_map: msprime.RateMap,
) -> tuple[msprime.RateMap, dict[str, Any], dict[str, Any], bool]:
    """Zero only ``[focal, focal + 1)`` and preserve every flank segment."""

    positions, rates = _rate_map_arrays(rate_map)
    before = _rate_map_metrics(rate_map)
    focal_start = float(FOCAL_POSITION_BP)
    focal_end = focal_start + 1.0
    positive = before["focal_integral"] > 0.0
    if not positive:
        return rate_map, before, dict(before), False
    if positions[0] > focal_start or positions[-1] < focal_end:
        raise ValueError("positive focal mutation rate is not on a full focal base")
    new_positions = np.unique(
        np.concatenate((positions, np.asarray([focal_start, focal_end])))
    )
    new_rates: list[float] = []
    for left, right in pairwise(new_positions):
        midpoint = float(left + (right - left) / 2)
        old_index = int(np.searchsorted(positions, midpoint, side="right") - 1)
        old_rate = float(rates[old_index])
        new_rates.append(
            0.0 if left >= focal_start and right <= focal_end else old_rate
        )
    replacement = msprime.RateMap(position=new_positions, rate=new_rates)
    after = _rate_map_metrics(replacement)
    if after["focal_integral"] != 0.0:
        raise RuntimeError("reserved focal-base mutation rate remains positive")
    for key in ("left_flank_integral", "right_flank_integral"):
        if not math.isclose(
            float(before[key]), float(after[key]), rel_tol=1e-15, abs_tol=1e-18
        ):
            raise RuntimeError("focal masking changed a flank integral")
    comparison_breaks = np.unique(np.concatenate((positions, new_positions)))
    for left, right in pairwise(comparison_breaks):
        if left >= focal_start and right <= focal_end:
            continue
        midpoint = float(left + (right - left) / 2)
        if rate_map.get_rate(midpoint) != replacement.get_rate(midpoint):
            raise RuntimeError("focal masking changed a flank rate segment")
    if not math.isclose(
        float(after["total_integral"]),
        float(before["total_integral"]) - float(before["focal_integral"]),
        rel_tol=1e-15,
        abs_tol=1e-18,
    ):
        raise RuntimeError("focal masking changed mutation mass outside the focal base")
    return replacement, before, after, True


def _model_class(value: Any) -> str:
    return f"{type(value).__module__}.{type(value).__qualname__}"


def _model_integer(value: Any, attribute: str) -> int | None:
    observed = getattr(value, attribute, None)
    if observed is None:
        return None
    if isinstance(observed, bool) or int(observed) != observed:
        raise ValueError(f"mutation model {attribute} is not an integer")
    return int(observed)


def _overlay_layer(
    start_time: float | None,
    end_time: float | None,
    rate_map: Mapping[str, Any] | None,
) -> str:
    if rate_map is None:
        return "non_rate_map"
    if start_time is not None and end_time is None:
        return "recapitation"
    if start_time is None and end_time is not None:
        ordinary = _expected_overlay_inventory()
        if ordinary is not None and (
            rate_map["positions"] == ordinary[0]["original_positions"]
            and rate_map["rates"] == ordinary[0]["original_rates"]
        ):
            return "ordinary_neutral_dfe"
        if rate_map["focal_integral"] > 0 and math.isclose(
            rate_map["total_integral"], rate_map["focal_integral"]
        ):
            return "single_site_neutral_dfe"
        return "other_post_slim_dfe"
    return "other_rate_map"


def _optional_number(value: Any) -> int | float | None:
    if value is None:
        return None
    if isinstance(value, (int, np.integer)) and not isinstance(value, bool):
        return int(value)
    result = float(value)
    if not math.isfinite(result):
        raise ValueError("mutation-overlay call contains a nonfinite scalar")
    return result


def _validate_focal_overlay_patch_receipt(
    receipt: Mapping[str, Any], contract: Mapping[str, Any]
) -> None:
    required = {
        "schema",
        "contract",
        "contract_sha256",
        "intercepted_call_count",
        "patched_call_count",
        "all_positive_focal_rates_masked",
        "callable_restored",
        "calls",
    }
    if set(receipt) != required or receipt.get("schema") != FOCAL_OVERLAY_PATCH_SCHEMA:
        raise ValueError("focal-overlay patch receipt fields/schema differ")
    if receipt.get("contract") != contract or receipt.get(
        "contract_sha256"
    ) != _canonical_sha256(contract):
        raise ValueError("focal-overlay patch contract differs")
    calls = receipt.get("calls")
    if not isinstance(calls, list) or int(receipt["intercepted_call_count"]) != len(
        calls
    ):
        raise ValueError("focal-overlay intercepted call inventory differs")
    patched = 0
    call_fields = {
        "call_index",
        "layer",
        "rate_kind",
        "model_class",
        "model_type",
        "model_next_id",
        "start_time",
        "end_time",
        "random_seed",
        "keep",
        "discrete_genome",
        "record_provenance",
        "non_rate_arguments_preserved_by_identity",
        "model_identity_preserved",
        "rate_argument_replaced",
        "positive_focal_rate_before",
        "positive_focal_rate_after",
        "original_focal_rate",
        "passed_focal_rate",
        "original_rate_map",
        "passed_rate_map",
        "call_succeeded",
    }
    for index, call in enumerate(calls, start=1):
        if not isinstance(call, Mapping) or set(call) != call_fields:
            raise ValueError("focal-overlay call receipt fields differ")
        if int(call["call_index"]) != index or not call["call_succeeded"]:
            raise ValueError("focal-overlay call order/completion differs")
        if (
            not call["non_rate_arguments_preserved_by_identity"]
            or not call["model_identity_preserved"]
        ):
            raise ValueError("focal-overlay patch changed a non-rate argument")
        is_rate_map = call["rate_kind"] == "msprime.RateMap"
        if is_rate_map:
            original = call["original_rate_map"]
            passed = call["passed_rate_map"]
            try:
                original_map = msprime.RateMap(
                    position=original["positions"], rate=original["rates"]
                )
                expected_map, before, after, should_patch = _mask_focal_rate_map(
                    original_map
                )
            except (KeyError, TypeError, ValueError) as error:
                raise ValueError(
                    "focal-overlay rate-map receipt is malformed"
                ) from error
            expected = _rate_map_metrics(expected_map)
            if original != before or passed != expected or expected != after:
                raise ValueError("focal-overlay rate-map receipt differs")
            if call["layer"] != _overlay_layer(
                call["start_time"], call["end_time"], before
            ):
                raise ValueError("focal-overlay layer receipt differs")
            if bool(call["rate_argument_replaced"]) != should_patch:
                raise ValueError("focal-overlay replacement flag differs")
            if bool(call["positive_focal_rate_before"]) != should_patch:
                raise ValueError("focal-overlay original focal-rate flag differs")
            if bool(call["positive_focal_rate_after"]):
                raise ValueError("focal-overlay passed a positive focal rate")
            if (
                call["original_focal_rate"] != before["focal_integral"]
                or call["passed_focal_rate"] != after["focal_integral"]
            ):
                raise ValueError("focal-overlay focal-rate receipt differs")
            patched += int(should_patch)
        elif call["layer"] != "non_rate_map" or any(
            (
                call["original_rate_map"] is not None,
                call["passed_rate_map"] is not None,
                call["rate_argument_replaced"],
                call["positive_focal_rate_before"],
                call["positive_focal_rate_after"],
                call["original_focal_rate"] is not None,
                call["passed_focal_rate"] is not None,
            )
        ):
            raise ValueError("non-RateMap mutation call has a mask receipt")
    if patched != int(receipt["patched_call_count"]):
        raise ValueError("focal-overlay patched call count differs")
    minimum = int(contract["minimum_patched_call_count"])
    if patched < minimum or not receipt["all_positive_focal_rates_masked"]:
        raise ValueError("focal-overlay guard did not mask every positive focal rate")
    expected_calls = contract["expected_call_count"]
    expected_patched = contract["expected_patched_call_count"]
    if expected_calls is not None and len(calls) != int(expected_calls):
        raise ValueError("focal-overlay exact call count differs")
    if expected_patched is not None and patched != int(expected_patched):
        raise ValueError("focal-overlay exact patched call count differs")
    expected_inventory = contract["expected_call_inventory"]
    if expected_inventory is not None:
        if len(expected_inventory) != len(calls):
            raise ValueError("focal-overlay exact call inventory length differs")
        shared_tick = calls[0]["end_time"]
        if (
            not isinstance(shared_tick, (int, float))
            or isinstance(shared_tick, bool)
            or not math.isfinite(float(shared_tick))
            or float(shared_tick) <= 0
            or calls[1]["start_time"] != shared_tick
        ):
            raise ValueError("focal-overlay calls do not share the finite SLiM tick")
        random_seeds: list[int] = []
        for call, expected in zip(calls, expected_inventory, strict=True):
            original = call["original_rate_map"]
            passed = call["passed_rate_map"]
            expected_start = (
                shared_tick
                if expected["start_time"] == "finite_shared_slim_tick"
                else None
            )
            expected_end = (
                shared_tick
                if expected["end_time"] == "finite_shared_slim_tick"
                else None
            )
            if (
                call["model_class"] != expected["model_class"]
                or call["layer"] != expected["layer"]
                or call["model_type"] != expected["model_type"]
                or not isinstance(call["model_next_id"], int)
                or isinstance(call["model_next_id"], bool)
                or call["model_next_id"] < 0
                or call["start_time"] != expected_start
                or call["end_time"] != expected_end
                or original["positions"] != expected["original_positions"]
                or original["rates"] != expected["original_rates"]
                or passed["positions"] != expected["passed_positions"]
                or passed["rates"] != expected["passed_rates"]
                or bool(call["rate_argument_replaced"]) != expected["patched"]
                or call["keep"] is not True
                or call["record_provenance"] is not True
                or call["discrete_genome"] is not None
            ):
                raise ValueError("focal-overlay exact call inventory shape differs")
            seed = call["random_seed"]
            if (
                not isinstance(seed, int)
                or isinstance(seed, bool)
                or not 1 <= seed < 2**32
            ):
                raise ValueError("focal-overlay mutation random seed differs")
            random_seeds.append(seed)
        if len(set(random_seeds)) != len(random_seeds):
            raise ValueError("focal-overlay mutation random seeds are not unique")
        if calls[0]["start_time"] is not None or calls[1]["end_time"] is not None:
            raise ValueError("focal-overlay start/end layer ordering differs")
    if not receipt["callable_restored"]:
        raise ValueError("msprime.sim_mutations was not restored")


@contextmanager
def _scoped_focal_overlay_patch():
    """Mask focal mutation mass for all overlays in one engine call."""

    original = msprime.sim_mutations
    signature = inspect.signature(original)
    contract = _focal_overlay_patch_contract()
    receipt: dict[str, Any] = {
        "schema": FOCAL_OVERLAY_PATCH_SCHEMA,
        "contract": contract,
        "contract_sha256": _canonical_sha256(contract),
        "intercepted_call_count": 0,
        "patched_call_count": 0,
        "all_positive_focal_rates_masked": True,
        "callable_restored": False,
        "calls": [],
    }

    def guarded(*args: Any, **kwargs: Any):
        bound = signature.bind(*args, **kwargs)
        bound.apply_defaults()
        rate = bound.arguments["rate"]
        model = bound.arguments["model"]
        replacement = rate
        before: dict[str, Any] | None = None
        after: dict[str, Any] | None = None
        patched = False
        if isinstance(rate, msprime.RateMap):
            replacement, before, after, patched = _mask_focal_rate_map(rate)
        call_args = list(args)
        call_kwargs = dict(kwargs)
        if patched:
            if len(call_args) >= 2:
                call_args[1] = replacement
            else:
                call_kwargs["rate"] = replacement
        passed = signature.bind(*call_args, **call_kwargs)
        passed.apply_defaults()
        non_rate_preserved = all(
            bound.arguments[name] is passed.arguments[name]
            for name in bound.arguments
            if name != "rate"
        )
        start_time = _optional_number(bound.arguments["start_time"])
        end_time = _optional_number(bound.arguments["end_time"])
        call = {
            "call_index": len(receipt["calls"]) + 1,
            "layer": _overlay_layer(start_time, end_time, before),
            "rate_kind": (
                "msprime.RateMap"
                if isinstance(rate, msprime.RateMap)
                else _model_class(rate)
            ),
            "model_class": _model_class(model),
            "model_type": _model_integer(model, "type"),
            "model_next_id": _model_integer(model, "next_id"),
            "start_time": start_time,
            "end_time": end_time,
            "random_seed": _optional_number(bound.arguments["random_seed"]),
            "keep": bound.arguments["keep"],
            "discrete_genome": bound.arguments["discrete_genome"],
            "record_provenance": bound.arguments["record_provenance"],
            "non_rate_arguments_preserved_by_identity": non_rate_preserved,
            "model_identity_preserved": model is passed.arguments["model"],
            "rate_argument_replaced": patched,
            "positive_focal_rate_before": bool(
                before is not None and before["focal_integral"] > 0.0
            ),
            "positive_focal_rate_after": bool(
                after is not None and after["focal_integral"] > 0.0
            ),
            "original_focal_rate": (
                None if before is None else before["focal_integral"]
            ),
            "passed_focal_rate": None if after is None else after["focal_integral"],
            "original_rate_map": before,
            "passed_rate_map": after,
            "call_succeeded": False,
        }
        receipt["calls"].append(call)
        receipt["intercepted_call_count"] += 1
        receipt["patched_call_count"] += int(patched)
        if call["positive_focal_rate_after"]:
            receipt["all_positive_focal_rates_masked"] = False
            raise RuntimeError("positive focal mutation rate escaped overlay guard")
        result = original(*call_args, **call_kwargs)
        call["call_succeeded"] = True
        return result

    msprime.sim_mutations = guarded
    failed = False
    try:
        yield receipt
    except BaseException:
        failed = True
        raise
    finally:
        msprime.sim_mutations = original
        receipt["callable_restored"] = msprime.sim_mutations is original
        if not failed:
            _validate_focal_overlay_patch_receipt(receipt, contract)


def _utc_now() -> str:
    return datetime.now(UTC).isoformat(timespec="seconds").replace("+00:00", "Z")


def _execution_units(replicates: int = DEFAULT_REPLICATES) -> pd.DataFrame:
    if replicates != DEFAULT_REPLICATES:
        raise ValueError(
            f"the frozen study requires exactly {DEFAULT_REPLICATES} replicates"
        )
    rows = []
    for replicate_index in range(DEFAULT_REPLICATES):
        unit_id = f"{STUDY_ID}__neutral__rep{replicate_index:03d}"
        rows.append(
            {
                "unit_id": unit_id,
                "replicate_index": replicate_index,
                "seed": _stable_seed(unit_id),
                "panel_seed": _stable_seed(f"{unit_id}:panel"),
                "scenario": SCENARIO,
                "selection_coefficient": 0.0,
                "mutation_origin_population": ARCHAIC_POPULATION,
                "mutation_age_generations": MUTATION_ORIGIN_GENERATIONS,
                "pulse_generations": INTROGRESSION_PULSE_GENERATIONS,
                "present_conditioning": "strict_segregating_in_present_han_census",
                "pool_frequency_conditioning": False,
                "terminal_af_target_or_band_conditioning": False,
                "sample_detection_conditioning": False,
                "pool_diploids": DEFAULT_POOL_DIPLOIDS,
                "sample_diploids": DEFAULT_PANEL_DIPLOIDS,
            }
        )
    frame = pd.DataFrame(rows, columns=EXECUTION_COLUMNS)
    if frame["unit_id"].duplicated().any() or frame["seed"].duplicated().any():
        raise RuntimeError("execution unit IDs/seeds are not unique")
    return frame


def _census_end(han_population_id: int) -> str:
    return f"""// Record exact present-Han census segregation and finish.
function (void)end(void) {{
    muts = sim.mutationsOfType(m1);
    if (size(muts) != 1)
        err("present-day state does not contain exactly one focal mutation object");
    han_matches = sim.subpopulations[sim.subpopulations.id == {han_population_id}];
    if (size(han_matches) != 1)
        err("present-day state does not contain exactly one Han population");
    han = han_matches[0];
    total_count = size(han.genomes);
    alt_count = sum(han.genomes.containsMutations(muts[0]));
    if ((total_count <= 0) | (alt_count <= 0) | (alt_count >= total_count))
        err("present-day Han census is not strictly segregating");
    final_af = alt_count / total_count;
    metadata.setValue("survival_neutral_final_census_alt_count", alt_count);
    metadata.setValue("survival_neutral_final_census_total_count", total_count);
    metadata.setValue("survival_neutral_final_census_af", final_af);
    metadata.setValue("survival_neutral_final_census_segregating", (alt_count > 0) & (alt_count < total_count));
    metadata.setValue("survival_neutral_final_census_fixed", alt_count == total_count);
    sim.treeSeqOutput(trees_file, metadata=metadata);
    sim.simulationFinished();
}}"""


@contextmanager
def _scoped_census_end_patch(han_population_id: int):
    """Patch only stdpopsim's end hook in this worker and restore it on exit."""

    original_object = slim_engine._slim_functions
    original = str(original_object)
    count = original.count(_ORIGINAL_END)
    if count != 1:
        raise RuntimeError(f"stdpopsim end-hook contract changed (matches={count})")
    replacement = _census_end(han_population_id)
    patched = original.replace(_ORIGINAL_END, replacement, 1)
    slim_engine._slim_functions = patched
    try:
        yield {
            "original_sha256": hashlib.sha256(original.encode()).hexdigest(),
            "patched_sha256": hashlib.sha256(patched.encode()).hexdigest(),
            "replacement_sha256": hashlib.sha256(replacement.encode()).hexdigest(),
            "replacement_count": 1,
            "han_population_id": han_population_id,
            "mutation_draw_helper_changed": False,
        }
    finally:
        slim_engine._slim_functions = original_object


def _census_end_patch_contract(han_population_id: int) -> dict[str, Any]:
    original = str(slim_engine._slim_functions)
    if original.count(_ORIGINAL_END) != 1:
        raise RuntimeError("stdpopsim end-hook contract changed")
    replacement = _census_end(han_population_id)
    patched = original.replace(_ORIGINAL_END, replacement, 1)
    return {
        "original_sha256": hashlib.sha256(original.encode()).hexdigest(),
        "patched_sha256": hashlib.sha256(patched.encode()).hexdigest(),
        "replacement_sha256": hashlib.sha256(replacement.encode()).hexdigest(),
        "replacement_count": 1,
        "han_population_id": han_population_id,
        "mutation_draw_helper_changed": False,
    }


def _slim_user_metadata(metadata: Any) -> Mapping[str, Any]:
    if not isinstance(metadata, Mapping):
        raise TypeError("tree-sequence metadata must be a mapping")
    slim = metadata.get("SLiM")
    user = slim.get("user_metadata") if isinstance(slim, Mapping) else None
    if not isinstance(user, Mapping):
        raise TypeError("tree metadata lacks exact SLiM.user_metadata mapping")
    if not _CENSUS_METADATA_KEYS.issubset(user):
        raise ValueError("tree metadata lacks exact present-Han census keys")
    for key in _CENSUS_METADATA_KEYS:
        occurrences = 0

        def count(candidate: Any, target_key: str = key) -> None:
            nonlocal occurrences
            if isinstance(candidate, Mapping):
                occurrences += int(target_key in candidate)
                for nested in candidate.values():
                    count(nested, target_key)
            elif isinstance(candidate, Sequence) and not isinstance(
                candidate, (str, bytes)
            ):
                for nested in candidate:
                    count(nested, target_key)

        count(metadata)
        if occurrences != 1:
            raise ValueError(f"tree metadata key {key} is not uniquely scoped")
    return user


def _metadata_scalar(metadata: Any, key: str) -> float:
    if key not in _CENSUS_METADATA_KEYS:
        raise ValueError(f"uncontracted census metadata key: {key}")
    value: Any = _slim_user_metadata(metadata)[key]
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


def _population_id(model: stdpopsim.DemographicModel, name: str) -> int:
    matches = [
        index
        for index, population in enumerate(model.model.populations)
        if population.name == name
    ]
    if len(matches) != 1:
        raise ValueError(f"demographic population {name!r} is not unique")
    return matches[0]


def study_root(repo_root: str | Path, root: str | Path = DEFAULT_STUDY_ROOT) -> Path:
    base = Path(repo_root).resolve()
    candidate = Path(root)
    return (
        (base / candidate).resolve()
        if not candidate.is_absolute()
        else candidate.resolve()
    )


def _plan_path(root: Path) -> Path:
    return root / "plan/study_plan.json"


def _plan_contract(
    repo: Path, decoder: Path, execution_units_sha256: str
) -> dict[str, Any]:
    _events, transform = build_neutral_event_schedule()
    return {
        "study_id": STUDY_ID,
        "expected_units": DEFAULT_REPLICATES,
        "estimand": ESTIMAND,
        "simulation_parameters": SIMULATION_PARAMETERS,
        "decoder_settings": DECODER_SETTINGS,
        "decoder_binary": _binary_record(
            decoder, repo, expected_sha256=EXPECTED_DECODER_SHA256, label="Gamma-SMC"
        ),
        "execution_units_sha256": execution_units_sha256,
        "source_hashes": _source_hashes(repo),
        "event_transform": transform,
        "simulation_law": SIMULATION_LAW,
        "focal_overlay_patch": _focal_overlay_patch_contract(),
        "individual_location_canonicalization": _individual_location_contract(),
        "runtime_versions": _runtime_versions(),
        "base_seed": BASE_SEED,
        "seed_derivation": (
            "1 + int.from_bytes(sha256(f'{base_seed}:{study_id}:{label}'.encode('utf-8'))[:8], "
            "'big') % (2**31 - 2)"
        ),
    }


def write_plan(
    repo_root: str | Path,
    *,
    root: str | Path = DEFAULT_STUDY_ROOT,
    decoder_bin: str | Path | None = None,
    force: bool = False,
) -> Path:
    """Create the immutable 100-unit execution plan and its exact receipt."""

    repo = Path(repo_root).resolve()
    destination = study_root(repo, root)
    plan_dir = destination / "plan"
    plan_path = _plan_path(destination)
    decoder = Path(decoder_bin or repo / "bin/gamma_smc").resolve()
    units = _execution_units()
    plan_dir.mkdir(parents=True, exist_ok=True)
    execution_path = plan_dir / "execution_units.tsv"
    if plan_path.exists() and not force:
        load_plan(repo, root=root)
        return plan_path
    if force:
        work = destination / "work"
        results = destination / "results"
        if (work.exists() and any(work.iterdir())) or (
            results.exists() and any(results.iterdir())
        ):
            raise ValueError("refusing to overwrite a plan after work/results exist")
    _atomic_frame(execution_path, units)
    contract = _plan_contract(repo, decoder, _sha256_file(execution_path))
    completion = {
        "schema": PLAN_SCHEMA,
        "status": "complete",
        "created_utc": _utc_now(),
        "contract": contract,
        "contract_sha256": _canonical_sha256(contract),
        "outputs": {"execution_units": _output_record(execution_path, rows=len(units))},
    }
    _atomic_json(plan_path, completion)
    return plan_path


def load_plan(
    repo_root: str | Path, *, root: str | Path = DEFAULT_STUDY_ROOT
) -> dict[str, Any]:
    destination = study_root(repo_root, root)
    path = _plan_path(destination)
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError("study plan is absent or unreadable") from error
    if set(payload) != {
        "schema",
        "status",
        "created_utc",
        "contract",
        "contract_sha256",
        "outputs",
    }:
        raise ValueError("study plan top-level fields differ")
    if payload.get("schema") != PLAN_SCHEMA or payload.get("status") != "complete":
        raise ValueError("study plan schema/status differs")
    if not isinstance(payload.get("created_utc"), str) or not payload[
        "created_utc"
    ].endswith("Z"):
        raise ValueError("study plan creation timestamp differs")
    contract = payload.get("contract")
    if not isinstance(contract, Mapping) or payload.get(
        "contract_sha256"
    ) != _canonical_sha256(contract):
        raise ValueError("study plan contract checksum failed")
    execution_path = destination / "plan/execution_units.tsv"
    units = pd.read_csv(execution_path, sep="\t")
    expected = _execution_units()
    pd.testing.assert_frame_equal(units, expected, check_dtype=False)
    execution_sha = _sha256_file(execution_path)
    try:
        decoder_identity = str(contract["decoder_binary"]["identity"])
    except (KeyError, TypeError) as error:
        raise ValueError("study plan decoder identity is absent") from error
    repo = Path(repo_root).resolve()
    decoder = (repo / decoder_identity).resolve()
    expected_contract = _plan_contract(repo, decoder, execution_sha)
    if dict(contract) != expected_contract:
        raise ValueError(
            "study plan differs from current binaries/runtime/source/event contract"
        )
    outputs = payload.get("outputs")
    if not isinstance(outputs, Mapping) or set(outputs) != {"execution_units"}:
        raise ValueError("study plan output inventory differs")
    _validate_file_record(
        outputs["execution_units"], execution_path, expected_rows=DEFAULT_REPLICATES
    )
    result = dict(payload)
    result["units"] = units.to_dict(orient="records")
    return result


def _unit_dir(root: Path, unit_id: str) -> Path:
    return root / "work" / unit_id


def _focal_expectation(
    model: stdpopsim.DemographicModel, events: Sequence[Any]
) -> dict[str, Any]:
    template = build_raw_han_selected_sweep(
        _selected_template_unit(), slim_scaling_factor=DEFAULT_SLIM_SCALING_FACTOR
    )
    neutral = replace(template, extended_events=tuple(events))
    expectation = _selected_focal_expectation(
        model, neutral, DEFAULT_SLIM_SCALING_FACTOR
    )
    if (
        expectation["source_population"] != ARCHAIC_POPULATION
        or expectation["requested_origin_time_generations"]
        != MUTATION_ORIGIN_GENERATIONS
        or expectation["realized_origin_time_generations"]
        != MUTATION_ORIGIN_GENERATIONS
        or expectation["expected_selection_coeff"] != 0.0
    ):
        raise RuntimeError(
            "focal identity expectation differs from the frozen neutral origin"
        )
    return expectation


def _simulation_contract(
    repo_root: Path,
    root: Path,
    plan: Mapping[str, Any],
    unit: Mapping[str, Any],
    slim_bin: Path,
) -> tuple[dict[str, Any], tuple[Any, ...]]:
    events, transform = build_neutral_event_schedule()
    model = load_ancient_eurasia_model()
    contract = {
        "unit_id": str(unit["unit_id"]),
        "replicate_index": int(unit["replicate_index"]),
        "seed": int(unit["seed"]),
        "plan_sha256": str(plan["contract_sha256"]),
        "estimand": ESTIMAND,
        "event_schedule": transform["event_schedule"],
        "event_schedule_sha256": _canonical_sha256(transform["event_schedule"]),
        "simulation_parameters": SIMULATION_PARAMETERS,
        "source_hashes": dict(plan["contract"]["source_hashes"]),
        "panel_seed": int(unit["panel_seed"]),
        "slim_binary": _binary_record(
            slim_bin,
            repo_root,
            expected_sha256=EXPECTED_SLIM_SHA256,
            label="SLiM",
        ),
        "focal_identity_expectation": _focal_expectation(model, events),
        "focal_overlay_patch": _focal_overlay_patch_contract(),
        "individual_location_canonicalization": _individual_location_contract(),
        "slim_end_hook_patch": _census_end_patch_contract(
            _population_id(model, INTROGRESSION_TARGET_POPULATION)
        ),
    }
    return contract, events


def _validate_file_record(
    record: Any, path: Path, *, expected_rows: int | None = None
) -> None:
    if not isinstance(record, Mapping) or set(record) != {
        "path",
        "sha256",
        "size_bytes",
        "rows",
    }:
        raise ValueError(f"output record is malformed: {path.name}")
    if record["path"] != path.name or not path.is_file():
        raise ValueError(f"output path is absent or changed: {path.name}")
    if (
        record["sha256"] != _sha256_file(path)
        or record["size_bytes"] != path.stat().st_size
    ):
        raise ValueError(f"output checksum/size differs: {path.name}")
    tabular = path.name.endswith((".tsv", ".tsv.gz"))
    if tabular:
        rows = record["rows"]
        if isinstance(rows, bool) or not isinstance(rows, int) or rows < 0:
            raise ValueError(f"tabular output row declaration is invalid: {path.name}")
        if expected_rows is not None and rows != expected_rows:
            raise ValueError(f"output row declaration differs: {path.name}")
    elif record["rows"] is not None:
        raise ValueError(f"non-tabular output declares rows: {path.name}")


def _truth_class_summaries(truth: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (genotype_class, threshold), frame in truth.groupby(
        ["genotype_class", "threshold_years"], sort=True
    ):
        focal = frame[np.isclose(frame["position_0based"], FOCAL_POSITION_BP)]
        if len(focal) != 1:
            raise ValueError(
                "truth profile lacks one focal row for a represented class"
            )
        rows.append(
            {
                "source": "tree_truth",
                "genotype_class": str(genotype_class),
                "threshold_years": int(threshold),
                "n_pairs": int(focal.iloc[0]["n_pairs"]),
                "focal_p_tmrca_lt_threshold": float(
                    focal.iloc[0]["p_tmrca_lt_threshold"]
                ),
                "region_mean_p_tmrca_lt_threshold": float(
                    frame["p_tmrca_lt_threshold"].mean()
                ),
            }
        )
    return pd.DataFrame(rows)


def _score_rows(
    unit: Mapping[str, Any],
    *,
    source: str,
    probabilities: Mapping[int, float],
    endpoint: Mapping[str, Any],
) -> pd.DataFrame:
    rows = []
    for threshold in TMRCA_THRESHOLDS_YEARS:
        value = float(probabilities[int(threshold)])
        if not 0.0 <= value <= 1.0:
            raise ValueError("focal TMRCA probability is outside [0,1]")
        rows.append(
            {
                "unit_id": str(unit["unit_id"]),
                "replicate_index": int(unit["replicate_index"]),
                "seed": int(unit["seed"]),
                "source": source,
                "threshold_years": int(threshold),
                "overall_focal_p_tmrca_lt_threshold": value,
                "final_pool_alt_count": int(endpoint["pool_alt_count"]),
                "final_pool_total_count": int(endpoint["pool_total_count"]),
                "final_pool_af": float(endpoint["pool_af"]),
                "sample_alt_count": int(endpoint["sample_alt_count"]),
                "sample_total_count": int(endpoint["sample_total_count"]),
                "sample_af": float(endpoint["sample_af"]),
                "sample_detected": bool(endpoint["sample_detected"]),
            }
        )
    return pd.DataFrame(rows, columns=SCORE_COLUMNS)


def _truth_scores(
    unit: Mapping[str, Any], truth: pd.DataFrame, endpoint: Mapping[str, Any]
) -> pd.DataFrame:
    focal = truth[
        (truth["genotype_class"].astype(str) == "overall")
        & np.isclose(truth["position_0based"], FOCAL_POSITION_BP)
    ].sort_values("threshold_years")
    if len(focal) != len(TMRCA_THRESHOLDS_YEARS) or set(
        focal["threshold_years"].astype(int)
    ) != set(TMRCA_THRESHOLDS_YEARS):
        raise ValueError("truth profile lacks the exact focal threshold grid")
    if set(focal["n_pairs"].astype(int)) != {DEFAULT_PANEL_DIPLOIDS}:
        raise ValueError("truth focal reducer does not use exactly 100 pairs")
    probabilities = dict(
        zip(
            focal["threshold_years"].astype(int),
            focal["p_tmrca_lt_threshold"].astype(float),
            strict=True,
        )
    )
    return _score_rows(
        unit, source="tree_truth", probabilities=probabilities, endpoint=endpoint
    )


def _validate_endpoint(endpoint: Mapping[str, Any], unit: Mapping[str, Any]) -> None:
    expected = {
        "final_census_alt_count",
        "final_census_total_count",
        "final_census_af",
        "final_census_segregating",
        "pool_diploids",
        "pool_alt_count",
        "pool_total_count",
        "pool_af",
        "pool_detected",
        "pool_fixed",
        "sample_diploids",
        "sample_alt_count",
        "sample_total_count",
        "sample_af",
        "sample_detected",
        "sample_fixed",
        "panel_seed",
        "chosen_panel_rows_zero_based",
    }
    if set(endpoint) != expected:
        raise ValueError("simulation endpoint fields differ")
    census_alt = int(endpoint["final_census_alt_count"])
    census_total = int(endpoint["final_census_total_count"])
    if census_total != 2_520 or not 0 < census_alt < census_total:
        raise ValueError("Q=5 present-Han census is not strictly segregating")
    if not math.isclose(float(endpoint["final_census_af"]), census_alt / census_total):
        raise ValueError("present-Han census AF differs from exact counts")
    if endpoint["final_census_segregating"] is not True:
        raise ValueError("present-Han census segregation flag is false")
    chosen = endpoint["chosen_panel_rows_zero_based"]
    expected_chosen = sorted(
        int(value)
        for value in np.random.default_rng(int(unit["panel_seed"])).choice(
            DEFAULT_POOL_DIPLOIDS, size=DEFAULT_PANEL_DIPLOIDS, replace=False
        )
    )
    if chosen != expected_chosen:
        raise ValueError("chosen panel rows differ from deterministic uniform sampling")
    for prefix, n_diploids in (
        ("pool", DEFAULT_POOL_DIPLOIDS),
        ("sample", DEFAULT_PANEL_DIPLOIDS),
    ):
        alt = int(endpoint[f"{prefix}_alt_count"])
        total = int(endpoint[f"{prefix}_total_count"])
        if endpoint[f"{prefix}_diploids"] != n_diploids or total != 2 * n_diploids:
            raise ValueError(f"{prefix} denominator differs")
        if not 0 <= alt <= total or not math.isclose(
            float(endpoint[f"{prefix}_af"]), alt / total
        ):
            raise ValueError(f"{prefix} AF/count differs")
        if bool(endpoint[f"{prefix}_detected"]) != (alt > 0) or bool(
            endpoint[f"{prefix}_fixed"]
        ) != (alt == total):
            raise ValueError(f"{prefix} flags differ")


def _validate_simulation_completion(
    root: Path,
    repo_root: Path,
    plan: Mapping[str, Any],
    unit: Mapping[str, Any],
    slim_bin: Path,
) -> dict[str, Any]:
    directory = _unit_dir(root, str(unit["unit_id"]))
    path = directory / "simulation_complete.json"
    try:
        completion = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError("simulation completion is absent or unreadable") from error
    contract, _events = _simulation_contract(repo_root, root, plan, unit, slim_bin)
    if (
        completion.get("schema") != SIMULATION_SCHEMA
        or completion.get("status") != "complete"
        or completion.get("contract") != contract
        or completion.get("contract_sha256") != _canonical_sha256(contract)
    ):
        raise ValueError("simulation completion contract differs")
    standalone_contract = json.loads(
        (directory / "simulation_contract.json").read_text(encoding="utf-8")
    )
    if standalone_contract != contract:
        raise ValueError("standalone simulation contract differs from completion")
    if completion.get("slim_end_hook_patch") != contract["slim_end_hook_patch"]:
        raise ValueError("simulation end-hook patch receipt differs from contract")
    overlay_receipt = completion.get("focal_overlay_patch")
    if not isinstance(overlay_receipt, Mapping):
        raise TypeError("simulation focal-overlay patch receipt is absent")
    _validate_focal_overlay_patch_receipt(
        overlay_receipt, contract["focal_overlay_patch"]
    )
    location_receipt = completion.get("individual_location_canonicalization")
    if not isinstance(location_receipt, Mapping):
        raise TypeError("simulation individual-location receipt is absent")
    endpoint = completion.get("endpoint")
    if not isinstance(endpoint, Mapping):
        raise TypeError("simulation endpoint is absent")
    _validate_endpoint(endpoint, unit)
    outputs = completion.get("outputs")
    if not isinstance(outputs, Mapping) or set(outputs) != set(SIMULATION_OUTPUT_PATHS):
        raise ValueError("simulation output inventory differs")
    tabular_rows = {
        "sample_manifest": DEFAULT_PANEL_DIPLOIDS,
        "truth_unit_scores": len(TMRCA_THRESHOLDS_YEARS),
    }
    for key, name in SIMULATION_OUTPUT_PATHS.items():
        _validate_file_record(
            outputs[key], directory / name, expected_rows=tabular_rows.get(key)
        )
    tree = tskit.load(directory / "simulation.trees")
    if (
        tree.num_samples != 2 * DEFAULT_PANEL_DIPLOIDS
        or tree.sequence_length != SEQUENCE_LENGTH_BP
    ):
        raise ValueError("simulation tree dimensions differ")
    _validate_individual_location_receipt(
        location_receipt,
        contract["individual_location_canonicalization"],
        tree=tree,
    )
    tree_census = {
        "final_census_alt_count": _metadata_integer(
            tree.metadata, "survival_neutral_final_census_alt_count"
        ),
        "final_census_total_count": _metadata_integer(
            tree.metadata, "survival_neutral_final_census_total_count"
        ),
        "final_census_af": _metadata_scalar(
            tree.metadata, "survival_neutral_final_census_af"
        ),
        "final_census_segregating": bool(
            _metadata_integer(
                tree.metadata, "survival_neutral_final_census_segregating"
            )
        ),
    }
    tree_fixed = bool(
        _metadata_integer(tree.metadata, "survival_neutral_final_census_fixed")
    )
    for key, observed in tree_census.items():
        expected = endpoint[key]
        equal = (
            math.isclose(float(observed), float(expected), abs_tol=1e-12)
            if isinstance(observed, float)
            else observed == expected
        )
        if not equal:
            raise ValueError(f"tree census metadata differs from endpoint: {key}")
    if tree_fixed or tree_census["final_census_total_count"] != 2_520:
        raise ValueError("tree census metadata is fixed or has the wrong Q=5 size")
    identity = _validated_focal_identity(
        tree,
        contract["focal_identity_expectation"],
        pool_alt_count=int(endpoint["pool_alt_count"]),
    )
    if completion.get("focal_identity") != identity:
        raise ValueError("simulation focal identity receipt differs")
    manifest = pd.read_csv(directory / "sample_manifest.tsv", sep="\t")
    observed_counts = _focal_counts_or_zero(tree)
    if not np.array_equal(
        manifest["focal_selected_allele_count"].to_numpy(dtype=int), observed_counts
    ):
        raise ValueError("simulation tree and sample manifest focal counts differ")
    if int(observed_counts.sum()) != int(endpoint["sample_alt_count"]):
        raise ValueError("sample manifest/tree counts differ from endpoint")
    pool_records = _retained_present_han_pool(tree)
    pool_counts = _focal_counts_for_records(tree, pool_records)
    if int(pool_counts.sum()) != int(endpoint["pool_alt_count"]):
        raise ValueError("retained 500-diploid pool counts differ from endpoint")
    chosen = endpoint["chosen_panel_rows_zero_based"]
    expected_ids = [pool_records[int(row)][0] for row in chosen]
    if manifest["tree_sequence_individual_id"].astype(int).tolist() != expected_ids:
        raise ValueError("sample manifest IDs differ from chosen retained pool rows")
    if not np.array_equal(observed_counts, pool_counts[np.asarray(chosen, dtype=int)]):
        raise ValueError("sample focal counts differ from chosen retained pool rows")
    return completion


def _validated_focal_identity(
    tree: tskit.TreeSequence,
    expectation: Mapping[str, Any],
    *,
    pool_alt_count: int,
) -> dict[str, Any]:
    inspected = inspect_selected_focal_identity(tree, expectation)
    if inspected["status"] == "valid":
        return validate_selected_focal_identity(tree, expectation)
    if (
        int(pool_alt_count) == 0
        and inspected["status"] == "selected_focal_mutation_absent"
    ):
        return {
            **inspected,
            "validation": "causal_site_absent_from_observational_pool; exact census survival retained in SLiM metadata",
        }
    raise ValueError(
        "focal identity is invalid for the observed pool: " + str(inspected["status"])
    )


def _retained_present_han_pool(
    tree: tskit.TreeSequence,
) -> list[tuple[int, tuple[int, int]]]:
    han_id = _population_id(
        load_ancient_eurasia_model(), INTROGRESSION_TARGET_POPULATION
    )
    records: list[tuple[int, tuple[int, int]]] = []
    for individual in tree.individuals():
        nodes = tuple(int(node) for node in individual.nodes)
        if len(nodes) != 2:
            continue
        node_records = [tree.node(node) for node in nodes]
        if all(
            math.isclose(float(node.time), 0.0, abs_tol=1e-12)
            and int(node.population) == han_id
            for node in node_records
        ):
            records.append((int(individual.id), nodes))
    records.sort(key=lambda record: record[0])
    if len(records) != DEFAULT_POOL_DIPLOIDS:
        raise ValueError(
            "panel tree does not retain exactly 500 present-Han pool individuals"
        )
    return records


def _focal_counts_for_records(
    tree: tskit.TreeSequence, records: Sequence[tuple[int, tuple[int, int]]]
) -> np.ndarray:
    sites = [
        site
        for site in tree.sites()
        if math.isclose(float(site.position), FOCAL_POSITION_BP, abs_tol=1e-9)
    ]
    if not sites:
        return np.zeros(len(records), dtype=np.int8)
    if len(sites) != 1 or len(sites[0].mutations) != 1:
        raise ValueError("retained focal site is not unique and nonrecurrent")
    site = sites[0]
    derived = str(site.mutations[0].derived_state)
    ordered_nodes = [node for _, nodes in records for node in nodes]
    variant = tskit.Variant(tree, samples=ordered_nodes)
    variant.decode(site.id)
    derived_indices = [
        index for index, allele in enumerate(variant.alleles) if allele == derived
    ]
    if len(derived_indices) != 1 or derived_indices[0] == 0:
        raise ValueError("retained focal allele is ambiguous")
    genotypes = np.asarray(variant.genotypes, dtype=int)
    if np.any(genotypes < 0) or set(genotypes.tolist()).difference(
        {0, derived_indices[0]}
    ):
        raise ValueError("retained focal genotypes are missing or noncausal")
    return (
        (genotypes.reshape(len(records), 2) == derived_indices[0])
        .sum(axis=1)
        .astype(np.int8)
    )


def _overlay_inventory_launch_mode(
    repo: Path,
    destination: Path,
    *,
    selected: Sequence[Mapping[str, Any]],
    unit_ids: Sequence[str] | None,
) -> str:
    """Fail closed until the instrumented rep051 inventory is frozen."""

    counts = (EXPECTED_OVERLAY_CALL_COUNT, EXPECTED_OVERLAY_PATCHED_COUNT)
    if all(count is not None for count in counts):
        return "frozen"
    if any(count is not None for count in counts):
        raise RuntimeError("focal-overlay expected counts are only partially frozen")
    explicit = [] if unit_ids is None else list(dict.fromkeys(map(str, unit_ids)))
    selected_ids = [str(unit["unit_id"]) for unit in selected]
    default_root = study_root(repo, DEFAULT_STUDY_ROOT)
    allowed = (
        explicit == [INSTRUMENTED_SMOKE_UNIT_ID]
        and selected_ids == [INSTRUMENTED_SMOKE_UNIT_ID]
        and destination != default_root
        and destination.name == INSTRUMENTED_SMOKE_ROOT_BASENAME
    )
    work = destination / "work"
    foreign_work = (
        [
            path.name
            for path in work.iterdir()
            if path.name != INSTRUMENTED_SMOKE_UNIT_ID
        ]
        if work.is_dir()
        else []
    )
    if not allowed or foreign_work:
        raise RuntimeError(
            "provisional focal-overlay inventory allows only an explicit rep051 "
            f"run in an isolated {INSTRUMENTED_SMOKE_ROOT_BASENAME!r} root"
        )
    return "instrumented_rep051_smoke"


def _simulate_unit(payload: Mapping[str, Any]) -> dict[str, Any]:
    repo = Path(payload["repo_root"]).resolve()
    root = Path(payload["study_root"]).resolve()
    plan = dict(payload["plan"])
    unit = dict(payload["unit"])
    slim_bin = Path(payload["slim_bin"]).resolve()
    directory = _unit_dir(root, str(unit["unit_id"]))
    directory.mkdir(parents=True, exist_ok=True)
    contract, events = _simulation_contract(repo, root, plan, unit, slim_bin)
    overlay_record: dict[str, Any] | None = None
    location_record: dict[str, Any] | None = None
    inventory_mode = str(payload.get("overlay_inventory_mode", ""))
    if inventory_mode not in {"frozen", "instrumented_rep051_smoke"}:
        raise ValueError("simulation lacks a valid focal-overlay inventory launch mode")
    if (
        inventory_mode == "instrumented_rep051_smoke"
        and str(unit["unit_id"]) != INSTRUMENTED_SMOKE_UNIT_ID
    ):
        raise ValueError("provisional focal-overlay inventory permits only rep051")
    try:
        with exclusive_unit_lock(
            directory / ".simulation.lock",
            unit_id=str(unit["unit_id"]),
            phase="ancient_eurasia_survival_neutral_simulate",
        ):
            if (directory / "simulation_complete.json").is_file():
                completion = _validate_simulation_completion(
                    root, repo, plan, unit, slim_bin
                )
                return {
                    **unit,
                    "status": "cached",
                    "error": "",
                    "final_census_af": float(completion["endpoint"]["final_census_af"]),
                    "simulation_completion_sha256": _sha256_file(
                        directory / "simulation_complete.json"
                    ),
                }
            _atomic_json(directory / "simulation_contract.json", contract)
            started = perf_counter()
            model = load_ancient_eurasia_model()
            if not math.isclose(_present_han_ne(), 6_300.0, abs_tol=1e-12):
                raise ValueError("AncientEurasia present Han Ne differs from 6300")
            han_population_id = _population_id(model, INTROGRESSION_TARGET_POPULATION)
            engine = stdpopsim.get_engine("slim")
            slim_log = directory / "slim.tsv"
            with (
                _scoped_census_end_patch(han_population_id) as patch_record,
                _scoped_focal_overlay_patch() as overlay_record,
            ):
                pool_ts = engine.simulate(
                    model,
                    build_focal_contig(),
                    {INTROGRESSION_TARGET_POPULATION: DEFAULT_POOL_DIPLOIDS},
                    seed=int(unit["seed"]),
                    extended_events=events,
                    slim_path=str(slim_bin),
                    slim_scaling_factor=DEFAULT_SLIM_SCALING_FACTOR,
                    slim_burn_in=DEFAULT_SLIM_BURN_IN,
                    logfile=str(slim_log),
                    logfile_interval=100,
                    keep_mutation_ids_as_alleles=False,
                )
            census_alt = _metadata_integer(
                pool_ts.metadata, "survival_neutral_final_census_alt_count"
            )
            census_total = _metadata_integer(
                pool_ts.metadata, "survival_neutral_final_census_total_count"
            )
            census_af = _metadata_scalar(
                pool_ts.metadata, "survival_neutral_final_census_af"
            )
            census_segregating = bool(
                _metadata_integer(
                    pool_ts.metadata, "survival_neutral_final_census_segregating"
                )
            )
            census_fixed = bool(
                _metadata_integer(
                    pool_ts.metadata, "survival_neutral_final_census_fixed"
                )
            )
            if (
                census_total != 2_520
                or not 0 < census_alt < census_total
                or not math.isclose(census_af, census_alt / census_total)
                or not census_segregating
                or census_fixed
            ):
                raise ValueError(
                    "exact Q=5 present-Han census metadata is inconsistent"
                )
            pool_records = _sampled_diploids(pool_ts)
            if len(pool_records) != DEFAULT_POOL_DIPLOIDS:
                raise ValueError("SLiM observational pool is not exactly 500 diploids")
            pool_counts = _focal_counts_or_zero(pool_ts)
            if len(pool_counts) != DEFAULT_POOL_DIPLOIDS:
                raise ValueError("pool focal counts do not cover exactly 500 diploids")
            panel_ts, panel_counts, chosen_rows = _uniform_panel(
                pool_ts, pool_counts, seed=int(unit["panel_seed"])
            )
            panel_ts, location_record = _canonicalize_individual_locations(panel_ts)
            expected_chosen = sorted(
                int(value)
                for value in np.random.default_rng(int(unit["panel_seed"])).choice(
                    DEFAULT_POOL_DIPLOIDS,
                    size=DEFAULT_PANEL_DIPLOIDS,
                    replace=False,
                )
            )
            if chosen_rows != expected_chosen:
                raise RuntimeError(
                    "uniform panel helper changed its deterministic draw"
                )
            pair_table = build_diploid_pair_table(
                panel_ts,
                FOCAL_POSITION_BP,
                focal_genotype_counts=panel_counts,
            )
            expected_individual_ids = [pool_records[row][0] for row in chosen_rows]
            if (
                pair_table["tree_sequence_individual_id"].astype(int).tolist()
                != expected_individual_ids
            ):
                raise ValueError(
                    "panel manifest individual IDs differ from chosen pool rows"
                )
            if pair_table["vcf_diploid_index"].astype(int).tolist() != list(
                range(DEFAULT_PANEL_DIPLOIDS)
            ):
                raise ValueError("panel manifest VCF ordering differs")
            pool_alt = int(pool_counts.sum())
            sample_alt = int(panel_counts.sum())
            endpoint = {
                "final_census_alt_count": census_alt,
                "final_census_total_count": census_total,
                "final_census_af": census_af,
                "final_census_segregating": True,
                "pool_diploids": DEFAULT_POOL_DIPLOIDS,
                "pool_alt_count": pool_alt,
                "pool_total_count": 2 * DEFAULT_POOL_DIPLOIDS,
                "pool_af": pool_alt / (2 * DEFAULT_POOL_DIPLOIDS),
                "pool_detected": bool(pool_alt > 0),
                "pool_fixed": bool(pool_alt == 2 * DEFAULT_POOL_DIPLOIDS),
                "sample_diploids": DEFAULT_PANEL_DIPLOIDS,
                "sample_alt_count": sample_alt,
                "sample_total_count": 2 * DEFAULT_PANEL_DIPLOIDS,
                "sample_af": sample_alt / (2 * DEFAULT_PANEL_DIPLOIDS),
                "sample_detected": bool(sample_alt > 0),
                "sample_fixed": bool(sample_alt == 2 * DEFAULT_PANEL_DIPLOIDS),
                "panel_seed": int(unit["panel_seed"]),
                "chosen_panel_rows_zero_based": chosen_rows,
            }
            _validate_endpoint(endpoint, unit)
            focal_identity = _validated_focal_identity(
                panel_ts,
                contract["focal_identity_expectation"],
                pool_alt_count=pool_alt,
            )
            tree_path = directory / "simulation.trees"
            _atomic_tree(tree_path, panel_ts)
            _atomic_frame(directory / "sample_manifest.tsv", pair_table)
            positions = np.arange(
                0, SEQUENCE_LENGTH_BP, DEFAULT_OUTPUT_STRIDE_BP, dtype=float
            )
            truth = tree_truth_profiles(
                panel_ts,
                positions,
                pair_table,
                thresholds_years=TMRCA_THRESHOLDS_YEARS,
                generation_time_years=GENERATION_TIME_YEARS,
                require_nonempty=False,
            )
            truth_path = directory / "truth_profiles.tsv.gz"
            _atomic_frame(truth_path, truth)
            truth_classes = _truth_class_summaries(truth)
            _atomic_frame(directory / "truth_class_summaries.tsv", truth_classes)
            truth_scores = _truth_scores(unit, truth, endpoint)
            _atomic_frame(directory / "truth_unit_scores.tsv", truth_scores)
            output_rows = {
                "sample_manifest": len(pair_table),
                "slim_log": _headered_text_row_count(slim_log),
                "truth_profiles": len(truth),
                "truth_class_summaries": len(truth_classes),
                "truth_unit_scores": len(truth_scores),
            }
            outputs = {
                key: _output_record(directory / name, rows=output_rows.get(key))
                for key, name in SIMULATION_OUTPUT_PATHS.items()
            }
            completion = {
                "schema": SIMULATION_SCHEMA,
                "status": "complete",
                "contract": contract,
                "contract_sha256": _canonical_sha256(contract),
                "elapsed_seconds": perf_counter() - started,
                "endpoint": endpoint,
                "focal_identity": focal_identity,
                "focal_overlay_patch": overlay_record,
                "individual_location_canonicalization": location_record,
                "slim_end_hook_patch": patch_record,
                "outputs": outputs,
            }
            _atomic_json(directory / "simulation_complete.json", completion)
            (directory / "simulation_failed.json").unlink(missing_ok=True)
            _validate_simulation_completion(root, repo, plan, unit, slim_bin)
            return {
                **unit,
                "status": "complete",
                "error": "",
                "final_census_af": census_af,
                "simulation_completion_sha256": _sha256_file(
                    directory / "simulation_complete.json"
                ),
            }
    except UnitLockHeld as error:
        status = "locked"
        message = str(error)
    except Exception as error:  # noqa: BLE001 - persist arbitrary worker failures
        status = "failed"
        message = f"{type(error).__name__}: {error}"
        _atomic_json(
            directory / "simulation_failed.json",
            {
                "schema": SIMULATION_SCHEMA,
                "status": "failed",
                "unit_id": str(unit["unit_id"]),
                "error_type": type(error).__name__,
                "error": str(error),
                "contract_sha256": _canonical_sha256(contract),
                "focal_overlay_patch": overlay_record,
                "individual_location_canonicalization": location_record,
            },
        )
    return {
        **unit,
        "status": status,
        "error": message,
        "final_census_af": math.nan,
        "simulation_completion_sha256": "",
    }


def simulate_study(
    repo_root: str | Path,
    *,
    root: str | Path = DEFAULT_STUDY_ROOT,
    slim_bin: str | Path | None = None,
    workers: int = DEFAULT_WORKERS,
    unit_ids: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Run or resume all 100 direct-SLiM neutral simulations."""

    repo = Path(repo_root).resolve()
    destination = study_root(repo, root)
    plan = load_plan(repo, root=root)
    binary = Path(slim_bin or repo / ".native-stdpopsim/bin/slim").resolve()
    _binary_record(binary, repo, expected_sha256=EXPECTED_SLIM_SHA256, label="SLiM")
    n_workers = _validate_workers(workers)
    selected = _select_units(plan["units"], unit_ids)
    inventory_mode = _overlay_inventory_launch_mode(
        repo, destination, selected=selected, unit_ids=unit_ids
    )
    payloads = [
        {
            "repo_root": str(repo),
            "study_root": str(destination),
            "plan": plan,
            "unit": unit,
            "slim_bin": str(binary),
            "overlay_inventory_mode": inventory_mode,
        }
        for unit in selected
    ]
    if n_workers == 1:
        rows = [_simulate_unit(payload) for payload in payloads]
    else:
        rows = []
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = [executor.submit(_simulate_unit, payload) for payload in payloads]
            for future in as_completed(futures):
                rows.append(future.result())
    order = {str(unit["unit_id"]): index for index, unit in enumerate(plan["units"])}
    frame = (
        pd.DataFrame(rows)
        .sort_values("unit_id", key=lambda values: values.map(order), kind="stable")
        .reset_index(drop=True)
    )
    columns = (
        *EXECUTION_COLUMNS,
        "status",
        "error",
        "final_census_af",
        "simulation_completion_sha256",
    )
    frame = frame.loc[:, columns]
    _atomic_frame(destination / "results/simulation_status.tsv", frame)
    return frame


def simulation_status(
    repo_root: str | Path,
    *,
    root: str | Path = DEFAULT_STUDY_ROOT,
) -> pd.DataFrame:
    """Deep-check every existing simulation receipt without changing state."""

    repo = Path(repo_root).resolve()
    destination = study_root(repo, root)
    plan = load_plan(repo, root=root)
    rows = []
    for unit in plan["units"]:
        directory = _unit_dir(destination, str(unit["unit_id"]))
        contract_path = directory / "simulation_contract.json"
        if (
            directory / "simulation_complete.json"
        ).is_file() and contract_path.is_file():
            try:
                contract = json.loads(contract_path.read_text(encoding="utf-8"))
                identity = contract["slim_binary"]["identity"]
                slim_bin = (repo / str(identity)).resolve()
                completion = _validate_simulation_completion(
                    destination, repo, plan, unit, slim_bin
                )
                status = "complete"
                error = ""
                final_af = float(completion["endpoint"]["final_census_af"])
                receipt = _sha256_file(directory / "simulation_complete.json")
            except Exception as exc:  # noqa: BLE001 - status reports invalid receipts
                status = "invalid"
                error = f"{type(exc).__name__}: {exc}"
                final_af = math.nan
                receipt = ""
        elif (directory / "simulation_failed.json").is_file():
            failure = json.loads(
                (directory / "simulation_failed.json").read_text(encoding="utf-8")
            )
            status = "failed"
            error = f"{failure.get('error_type', 'Error')}: {failure.get('error', '')}"
            final_af = math.nan
            receipt = ""
        elif directory.exists() and any(directory.iterdir()):
            status = "partial"
            error = ""
            final_af = math.nan
            receipt = ""
        else:
            status = "not_started"
            error = ""
            final_af = math.nan
            receipt = ""
        rows.append(
            {
                **unit,
                "status": status,
                "error": error,
                "final_census_af": final_af,
                "simulation_completion_sha256": receipt,
            }
        )
    return pd.DataFrame(rows)


def _atomic_pair_file(path: Path, pair_table: pd.DataFrame) -> None:
    pairs = gamma_pair_tables(pair_table, require_nonempty=False)["overall"]
    if len(pairs) != DEFAULT_PANEL_DIPLOIDS:
        raise ValueError("explicit overall pair table is not exactly 100 pairs")
    expected = pd.DataFrame(
        {
            "haplotype_0": np.arange(0, 2 * DEFAULT_PANEL_DIPLOIDS, 2),
            "haplotype_1": np.arange(1, 2 * DEFAULT_PANEL_DIPLOIDS, 2),
        }
    )
    pd.testing.assert_frame_equal(pairs, expected, check_dtype=False)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        with temporary.open("w", encoding="utf-8", newline="\n") as stream:
            stream.write("# Gamma-SMC explicit haplotype pairs\n")
            stream.write("# genotype_class\toverall\n")
            stream.write(f"# n_pairs\t{DEFAULT_PANEL_DIPLOIDS}\n")
            pairs.to_csv(
                stream, sep="\t", header=False, index=False, lineterminator="\n"
            )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _read_pairs(path: Path) -> pd.DataFrame:
    pairs = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=None,
        names=["haplotype_0", "haplotype_1"],
    )
    if len(pairs) != DEFAULT_PANEL_DIPLOIDS:
        raise ValueError("overall pair file does not contain exactly 100 data rows")
    expected = np.arange(2 * DEFAULT_PANEL_DIPLOIDS, dtype=int).reshape(-1, 2)
    if not np.array_equal(pairs.to_numpy(dtype=int), expected):
        raise ValueError("overall pair file differs from within-diploid VCF ordering")
    return pairs


def _raw_input_record(
    path: Path, relative_path: str, *, rows: int | None
) -> dict[str, Any]:
    if not path.is_file() or path.stat().st_size <= 0:
        raise ValueError(f"raw decode input is absent: {path}")
    return {
        "path": relative_path,
        "sha256": _sha256_file(path),
        "size_bytes": path.stat().st_size,
        "rows": rows,
    }


def _decode_contract(
    repo: Path,
    root: Path,
    plan: Mapping[str, Any],
    unit: Mapping[str, Any],
    decoder_bin: Path,
    overall_pairs: Path,
) -> dict[str, Any]:
    directory = _unit_dir(root, str(unit["unit_id"]))
    manifest = directory / "sample_manifest.tsv"
    tree = directory / "simulation.trees"
    return {
        "unit_id": str(unit["unit_id"]),
        "replicate_index": int(unit["replicate_index"]),
        "seed": int(unit["seed"]),
        "plan_sha256": str(plan["contract_sha256"]),
        "simulation_completion_sha256": _sha256_file(
            directory / "simulation_complete.json"
        ),
        "estimand": ESTIMAND,
        "decoder_settings": DECODER_SETTINGS,
        "decoder_binary": _binary_record(
            decoder_bin,
            repo,
            expected_sha256=EXPECTED_DECODER_SHA256,
            label="Gamma-SMC",
        ),
        "source_hashes": dict(plan["contract"]["source_hashes"]),
        "raw_inputs": {
            "simulation_tree": _raw_input_record(
                tree, "../simulation.trees", rows=None
            ),
            "sample_manifest": _raw_input_record(
                manifest, "../sample_manifest.tsv", rows=DEFAULT_PANEL_DIPLOIDS
            ),
            "overall_pairs": _raw_input_record(
                overall_pairs, "overall.pairs.tsv", rows=DEFAULT_PANEL_DIPLOIDS
            ),
        },
    }


def _gamma_scores(
    unit: Mapping[str, Any], summary: pd.DataFrame, endpoint: Mapping[str, Any]
) -> pd.DataFrame:
    required = {
        "position_0based",
        "n_pairs",
        *(f"mean_p_lt_{threshold}" for threshold in TMRCA_THRESHOLDS_YEARS),
    }
    if not required.issubset(summary.columns):
        raise ValueError("Gamma summary lacks focal reducer columns")
    focal = summary[np.isclose(summary["position_0based"], FOCAL_POSITION_BP)]
    if len(focal) != 1 or int(focal.iloc[0]["n_pairs"]) != DEFAULT_PANEL_DIPLOIDS:
        raise ValueError("Gamma summary lacks one 100-pair focal output")
    probabilities = {
        int(threshold): float(focal.iloc[0][f"mean_p_lt_{threshold}"])
        for threshold in TMRCA_THRESHOLDS_YEARS
    }
    return _score_rows(
        unit, source="gamma_smc", probabilities=probabilities, endpoint=endpoint
    )


def _expected_decoder_command(decoder_bin: Path, directory: Path) -> list[str]:
    summary = directory / "decoded/overall.summary.tsv"
    posterior = directory / "decoded/posterior.zst"
    pair_path = directory / "decoded/overall.pairs.tsv"
    return [
        str(decoder_bin),
        "--input",
        "/dev/stdin",
        "--input_format",
        "vcf",
        "--scaled_mutation_rate",
        str(DECODER_SETTINGS["scaled_mutation_rate"]),
        "--recombination_to_mutation_ratio",
        str(DECODER_SETTINGS["recombination_to_mutation_ratio"]),
        "--unscaled_mutation_rate",
        str(MUTATION_RATE),
        "--recent_threshold_years",
        ",".join(str(value) for value in TMRCA_THRESHOLDS_YEARS),
        "--generation_time",
        str(GENERATION_TIME_YEARS),
        "--recent_summary",
        str(summary),
        "--recent_call",
        "median",
        "--recent_call_probability",
        "0.5",
        "--output_at_hets=false",
        "--output_at_stride",
        str(DEFAULT_OUTPUT_STRIDE_BP),
        "--cache_size",
        str(DEFAULT_CACHE_SIZE),
        "--threads",
        str(DEFAULT_DECODE_THREADS),
        "--pair_block",
        "256",
        "--backward_alignment",
        "fixed",
        "--exp10",
        "accurate",
        "--pairs_file",
        str(pair_path),
        "--output",
        str(posterior),
    ]


def _expected_producer_command(tree: Path) -> list[str]:
    script = (
        "import sys; "
        "from gamma_smc_aou.tree_sequence import stream_tree_sequence_vcf; "
        "stream_tree_sequence_vcf(sys.argv[1], sys.stdout, input_format=sys.argv[2])"
    )
    return [sys.executable, "-c", script, str(tree), "trees"]


def _validate_decoder_run(
    record: Mapping[str, Any],
    *,
    unit: Mapping[str, Any],
    contract: Mapping[str, Any],
    directory: Path,
    decoder_bin: Path,
) -> None:
    if set(record) != {"schema", "unit_id", "run", "raw_inputs"}:
        raise ValueError("decoder run receipt fields differ")
    if record["schema"] != DECODER_RUN_SCHEMA:
        raise ValueError("decoder run schema differs")
    if (
        record["unit_id"] != unit["unit_id"]
        or record["raw_inputs"] != contract["raw_inputs"]
    ):
        raise ValueError("decoder run raw-input binding differs")
    run = record["run"]
    if not isinstance(run, Mapping):
        raise TypeError("decoder run command record is absent")
    expected_command = _expected_decoder_command(decoder_bin, directory)
    command = [str(value) for value in run.get("command", [])]
    if command != expected_command:
        raise ValueError("decoder run command differs from every frozen setting/path")
    tree = directory / "simulation.trees"
    expected_producer = _expected_producer_command(tree)
    if run.get("tree_sequence_vcf_producer_command") != expected_producer:
        raise ValueError("decoder run tree-to-VCF producer command differs")
    if (
        run.get("input_format") != "trees"
        or Path(str(run.get("input_path"))).resolve() != tree.resolve()
    ):
        raise ValueError("decoder run tree input binding differs")
    exact_scalars = {
        "stride_bp": DEFAULT_OUTPUT_STRIDE_BP,
        "cache_size_bp": DEFAULT_CACHE_SIZE,
        "n_output_positions": SEQUENCE_LENGTH_BP // DEFAULT_OUTPUT_STRIDE_BP,
        "pairs_manifest": None,
        "n_pairs_recorded": None,
    }
    for key, expected in exact_scalars.items():
        if run.get(key) != expected:
            raise ValueError(f"decoder run scalar differs: {key}")
    seconds = run.get("decode_seconds")
    if (
        isinstance(seconds, bool)
        or not isinstance(seconds, (int, float))
        or not math.isfinite(seconds)
        or seconds < 0
    ):
        raise ValueError("decoder run elapsed time is invalid")
    for key in ("stdout", "stderr"):
        if not isinstance(run.get(key), str):
            raise TypeError(f"decoder run {key} is not text")


def _validate_decode_completion(
    root: Path,
    repo: Path,
    plan: Mapping[str, Any],
    unit: Mapping[str, Any],
    decoder_bin: Path,
) -> dict[str, Any]:
    directory = _unit_dir(root, str(unit["unit_id"]))
    decoded = directory / "decoded"
    overall_pairs = decoded / "overall.pairs.tsv"
    contract = _decode_contract(repo, root, plan, unit, decoder_bin, overall_pairs)
    path = decoded / "decode_complete.json"
    try:
        completion = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError("decode completion is absent or unreadable") from error
    if (
        completion.get("schema") != DECODE_SCHEMA
        or completion.get("status") != "complete"
        or completion.get("contract") != contract
        or completion.get("contract_sha256") != _canonical_sha256(contract)
    ):
        raise ValueError("decode completion contract differs")
    standalone_contract = json.loads(
        (decoded / "decode_contract.json").read_text(encoding="utf-8")
    )
    if standalone_contract != contract:
        raise ValueError("standalone decode contract differs from completion")
    outputs = completion.get("outputs")
    if not isinstance(outputs, Mapping) or set(outputs) != set(DECODE_OUTPUT_PATHS):
        raise ValueError("decode output inventory differs")
    expected_rows = {
        "overall_pairs": DEFAULT_PANEL_DIPLOIDS,
        "gamma_unit_scores": len(TMRCA_THRESHOLDS_YEARS),
    }
    for key, name in DECODE_OUTPUT_PATHS.items():
        _validate_file_record(
            outputs[key], decoded / name, expected_rows=expected_rows.get(key)
        )
    _read_pairs(overall_pairs)
    pair_table = pd.read_csv(directory / "sample_manifest.tsv", sep="\t")
    explicit = gamma_pair_tables(pair_table, require_nonempty=False)["overall"]
    persisted = _read_pairs(overall_pairs)
    pd.testing.assert_frame_equal(persisted, explicit, check_dtype=False)
    run = json.loads((decoded / "decoder_run.json").read_text(encoding="utf-8"))
    _validate_decoder_run(
        run,
        unit=unit,
        contract=contract,
        directory=directory,
        decoder_bin=decoder_bin,
    )
    endpoint = json.loads(
        (directory / "simulation_complete.json").read_text(encoding="utf-8")
    )["endpoint"]
    summary = pd.read_csv(decoded / "overall.summary.tsv", sep="\t")
    expected_scores = _gamma_scores(unit, summary, endpoint)
    observed_scores = pd.read_csv(decoded / "gamma_unit_scores.tsv", sep="\t")
    pd.testing.assert_frame_equal(
        observed_scores, expected_scores, check_dtype=False, rtol=1e-12, atol=1e-15
    )
    return completion


def _decode_unit(payload: Mapping[str, Any]) -> dict[str, Any]:
    repo = Path(payload["repo_root"]).resolve()
    root = Path(payload["study_root"]).resolve()
    plan = dict(payload["plan"])
    unit = dict(payload["unit"])
    decoder_bin = Path(payload["decoder_bin"]).resolve()
    directory = _unit_dir(root, str(unit["unit_id"]))
    decoded = directory / "decoded"
    if not (directory / "simulation_complete.json").is_file():
        return {
            **unit,
            "status": "not_simulated",
            "error": "",
            "decode_completion_sha256": "",
        }
    try:
        simulation_contract = json.loads(
            (directory / "simulation_contract.json").read_text(encoding="utf-8")
        )
        slim_bin = (repo / simulation_contract["slim_binary"]["identity"]).resolve()
        simulation = _validate_simulation_completion(root, repo, plan, unit, slim_bin)
        decoded.mkdir(parents=True, exist_ok=True)
        with exclusive_unit_lock(
            decoded / ".decode.lock",
            unit_id=str(unit["unit_id"]),
            phase="ancient_eurasia_survival_neutral_decode",
        ):
            if (decoded / "decode_complete.json").is_file():
                _validate_decode_completion(root, repo, plan, unit, decoder_bin)
                return {
                    **unit,
                    "status": "cached",
                    "error": "",
                    "decode_completion_sha256": _sha256_file(
                        decoded / "decode_complete.json"
                    ),
                }
            pair_table = pd.read_csv(directory / "sample_manifest.tsv", sep="\t")
            overall_pairs = decoded / "overall.pairs.tsv"
            _atomic_pair_file(overall_pairs, pair_table)
            contract = _decode_contract(
                repo, root, plan, unit, decoder_bin, overall_pairs
            )
            _atomic_json(decoded / "decode_contract.json", contract)
            started = perf_counter()
            summary_path = decoded / "overall.summary.tsv"
            posterior_path = decoded / "posterior.zst"
            run = run_within_decoder(
                decoder_bin,
                directory / "simulation.trees",
                summary_path,
                scaled_mutation_rate=DECODER_SETTINGS["scaled_mutation_rate"],
                recombination_to_mutation_ratio=DECODER_SETTINGS[
                    "recombination_to_mutation_ratio"
                ],
                mutation_rate=MUTATION_RATE,
                threshold_years=TMRCA_THRESHOLDS_YEARS,
                generation_time=GENERATION_TIME_YEARS,
                input_format="trees",
                raw_output=posterior_path,
                output_at_stride=DEFAULT_OUTPUT_STRIDE_BP,
                output_at_hets=False,
                only_within=False,
                pairs_file=overall_pairs,
                recent_call="median",
                recent_call_probability=0.5,
                threads=DEFAULT_DECODE_THREADS,
                cache_size=DEFAULT_CACHE_SIZE,
                pair_block=256,
                exp10="accurate",
                backward_alignment="fixed",
            )
            run_receipt = {
                "schema": DECODER_RUN_SCHEMA,
                "unit_id": str(unit["unit_id"]),
                "run": run,
                "raw_inputs": contract["raw_inputs"],
            }
            _atomic_json(decoded / "decoder_run.json", run_receipt)
            summary_run = summary_path.with_suffix(summary_path.suffix + ".run.json")
            summary_run.unlink(missing_ok=True)
            summary = pd.read_csv(summary_path, sep="\t")
            gamma_scores = _gamma_scores(unit, summary, simulation["endpoint"])
            _atomic_frame(decoded / "gamma_unit_scores.tsv", gamma_scores)
            output_rows = {
                "overall_pairs": DEFAULT_PANEL_DIPLOIDS,
                "overall_summary": len(summary),
                "gamma_unit_scores": len(gamma_scores),
            }
            outputs = {
                key: _output_record(decoded / name, rows=output_rows.get(key))
                for key, name in DECODE_OUTPUT_PATHS.items()
            }
            completion = {
                "schema": DECODE_SCHEMA,
                "status": "complete",
                "contract": contract,
                "contract_sha256": _canonical_sha256(contract),
                "elapsed_seconds": perf_counter() - started,
                "outputs": outputs,
            }
            _atomic_json(decoded / "decode_complete.json", completion)
            (decoded / "decode_failed.json").unlink(missing_ok=True)
            _validate_decode_completion(root, repo, plan, unit, decoder_bin)
            return {
                **unit,
                "status": "complete",
                "error": "",
                "decode_completion_sha256": _sha256_file(
                    decoded / "decode_complete.json"
                ),
            }
    except UnitLockHeld as error:
        status = "locked"
        message = str(error)
    except Exception as error:  # noqa: BLE001 - persist arbitrary worker failures
        status = "failed"
        message = f"{type(error).__name__}: {error}"
        decoded.mkdir(parents=True, exist_ok=True)
        _atomic_json(
            decoded / "decode_failed.json",
            {
                "schema": DECODE_SCHEMA,
                "status": "failed",
                "unit_id": str(unit["unit_id"]),
                "error_type": type(error).__name__,
                "error": str(error),
            },
        )
    return {
        **unit,
        "status": status,
        "error": message,
        "decode_completion_sha256": "",
    }


def _select_units(
    units: Sequence[Mapping[str, Any]], unit_ids: Sequence[str] | None
) -> list[dict[str, Any]]:
    records = [dict(unit) for unit in units]
    if not unit_ids:
        return records
    requested = list(dict.fromkeys(str(unit_id) for unit_id in unit_ids))
    by_id = {str(unit["unit_id"]): unit for unit in records}
    missing = sorted(set(requested).difference(by_id))
    if missing:
        raise ValueError(
            "requested unit IDs are absent from the plan: " + ", ".join(missing)
        )
    return [by_id[unit_id] for unit_id in requested]


def decode_study(
    repo_root: str | Path,
    *,
    root: str | Path = DEFAULT_STUDY_ROOT,
    decoder_bin: str | Path | None = None,
    workers: int = DEFAULT_WORKERS,
    threads: int = DEFAULT_DECODE_THREADS,
    unit_ids: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Run or resume exact-pair Gamma-SMC decoding for planned units."""

    if threads != DEFAULT_DECODE_THREADS:
        raise ValueError(
            "the frozen workflow requires exactly one Gamma-SMC thread per unit"
        )
    repo = Path(repo_root).resolve()
    destination = study_root(repo, root)
    plan = load_plan(repo, root=root)
    binary = Path(decoder_bin or repo / "bin/gamma_smc").resolve()
    _binary_record(
        binary, repo, expected_sha256=EXPECTED_DECODER_SHA256, label="Gamma-SMC"
    )
    selected = _select_units(plan["units"], unit_ids)
    n_workers = _validate_workers(workers)
    payloads = [
        {
            "repo_root": str(repo),
            "study_root": str(destination),
            "plan": plan,
            "unit": unit,
            "decoder_bin": str(binary),
        }
        for unit in selected
    ]
    if n_workers == 1:
        rows = [_decode_unit(payload) for payload in payloads]
    else:
        rows = []
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = [executor.submit(_decode_unit, payload) for payload in payloads]
            for future in as_completed(futures):
                rows.append(future.result())
    order = {str(unit["unit_id"]): index for index, unit in enumerate(plan["units"])}
    frame = (
        pd.DataFrame(rows)
        .sort_values("unit_id", key=lambda values: values.map(order), kind="stable")
        .reset_index(drop=True)
    )
    columns = (*EXECUTION_COLUMNS, "status", "error", "decode_completion_sha256")
    frame = frame.loc[:, columns]
    _atomic_frame(destination / "results/decode_status.tsv", frame)
    full_status = decode_status(repo, root=root)
    if set(full_status["status"].astype(str)) == {"complete"}:
        collect_results(repo, root=root)
    return frame


def decode_status(
    repo_root: str | Path,
    *,
    root: str | Path = DEFAULT_STUDY_ROOT,
) -> pd.DataFrame:
    """Deep-check every existing decode receipt without changing state."""

    repo = Path(repo_root).resolve()
    destination = study_root(repo, root)
    plan = load_plan(repo, root=root)
    rows = []
    for unit in plan["units"]:
        directory = _unit_dir(destination, str(unit["unit_id"]))
        decoded = directory / "decoded"
        contract_path = decoded / "decode_contract.json"
        if (decoded / "decode_complete.json").is_file() and contract_path.is_file():
            try:
                contract = json.loads(contract_path.read_text(encoding="utf-8"))
                decoder_bin = (repo / contract["decoder_binary"]["identity"]).resolve()
                _validate_decode_completion(destination, repo, plan, unit, decoder_bin)
                status = "complete"
                error = ""
                receipt = _sha256_file(decoded / "decode_complete.json")
            except Exception as exc:  # noqa: BLE001 - status reports invalid receipts
                status = "invalid"
                error = f"{type(exc).__name__}: {exc}"
                receipt = ""
        elif (decoded / "decode_failed.json").is_file():
            failure = json.loads(
                (decoded / "decode_failed.json").read_text(encoding="utf-8")
            )
            status = "failed"
            error = f"{failure.get('error_type', 'Error')}: {failure.get('error', '')}"
            receipt = ""
        elif decoded.exists() and any(decoded.iterdir()):
            status = "partial"
            error = ""
            receipt = ""
        elif (directory / "simulation_complete.json").is_file():
            status = "not_started"
            error = ""
            receipt = ""
        else:
            status = "not_simulated"
            error = ""
            receipt = ""
        rows.append(
            {
                **unit,
                "status": status,
                "error": error,
                "decode_completion_sha256": receipt,
            }
        )
    return pd.DataFrame(rows)


def _aggregate_payload(
    repo: Path, root: Path, plan: Mapping[str, Any]
) -> tuple[
    pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, list[dict[str, Any]]
]:
    simulation_rows = simulation_status(repo, root=root)
    decode_rows = decode_status(repo, root=root)
    if set(simulation_rows["status"].astype(str)) != {"complete"}:
        raise ValueError("all 100 simulations must be complete before aggregation")
    if set(decode_rows["status"].astype(str)) != {"complete"}:
        raise ValueError("all 100 decodes must be complete before aggregation")
    scores: list[pd.DataFrame] = []
    endpoints: list[dict[str, Any]] = []
    inventory: list[dict[str, Any]] = []
    for unit in plan["units"]:
        unit_id = str(unit["unit_id"])
        directory = _unit_dir(root, unit_id)
        decoded = directory / "decoded"
        sim_path = directory / "simulation_complete.json"
        decode_path = decoded / "decode_complete.json"
        simulation = json.loads(sim_path.read_text(encoding="utf-8"))
        decode = json.loads(decode_path.read_text(encoding="utf-8"))
        truth_scores = pd.read_csv(directory / "truth_unit_scores.tsv", sep="\t")
        gamma_scores = pd.read_csv(decoded / "gamma_unit_scores.tsv", sep="\t")
        scores.extend((truth_scores, gamma_scores))
        endpoint = dict(simulation["endpoint"])
        endpoint.pop("chosen_panel_rows_zero_based")
        endpoints.append(
            {
                "unit_id": unit_id,
                "replicate_index": int(unit["replicate_index"]),
                "seed": int(unit["seed"]),
                "simulation_completion_sha256": _sha256_file(sim_path),
                **endpoint,
            }
        )
        inventory.append(
            {
                "unit_id": unit_id,
                "replicate_index": int(unit["replicate_index"]),
                "seed": int(unit["seed"]),
                "simulation_completion_path": f"work/{unit_id}/simulation_complete.json",
                "simulation_completion_sha256": _sha256_file(sim_path),
                "decode_completion_path": f"work/{unit_id}/decoded/decode_complete.json",
                "decode_completion_sha256": _sha256_file(decode_path),
            }
        )
        if simulation.get("status") != "complete" or decode.get("status") != "complete":
            raise ValueError(f"unit completion status differs: {unit_id}")
    combined = pd.concat(scores, ignore_index=True).loc[:, SCORE_COLUMNS]
    combined = combined.sort_values(
        ["source", "unit_id", "threshold_years"], kind="stable"
    ).reset_index(drop=True)
    if len(combined) != 2 * DEFAULT_REPLICATES * len(TMRCA_THRESHOLDS_YEARS):
        raise ValueError("aggregate score row count differs")
    endpoint_frame = (
        pd.DataFrame(endpoints, columns=ENDPOINT_COLUMNS)
        .sort_values("unit_id", kind="stable")
        .reset_index(drop=True)
    )
    simulation_rows = simulation_rows.copy()
    decode_rows = decode_rows.copy()
    simulation_rows["status"] = "complete"
    decode_rows["status"] = "complete"
    return combined, endpoint_frame, simulation_rows, decode_rows, inventory


def collect_results(
    repo_root: str | Path,
    *,
    root: str | Path = DEFAULT_STUDY_ROOT,
) -> Path:
    """Atomically publish the compact 100-unit result receipts and scores."""

    repo = Path(repo_root).resolve()
    destination = study_root(repo, root)
    plan = load_plan(repo, root=root)
    scores, endpoints, simulation_rows, decode_rows, inventory = _aggregate_payload(
        repo, destination, plan
    )
    results = destination / "results"
    results.mkdir(parents=True, exist_ok=True)
    frames = {
        "endpoint_summary": endpoints,
        "replicate_scores": scores,
        "simulation_status": simulation_rows,
        "decode_status": decode_rows,
    }
    for key, name in RESULT_OUTPUT_PATHS.items():
        _atomic_frame(results / name, frames[key])
    execution_path = destination / "plan/execution_units.tsv"
    contract = {
        "study_id": STUDY_ID,
        "expected_units": DEFAULT_REPLICATES,
        "plan_sha256": str(plan["contract_sha256"]),
        "execution_units_sha256": _sha256_file(execution_path),
        "estimand": ESTIMAND,
        "simulation_law": SIMULATION_LAW,
        "decoder_settings": DECODER_SETTINGS,
        "decoder_binary": dict(plan["contract"]["decoder_binary"]),
        "source_hashes": dict(plan["contract"]["source_hashes"]),
        "schemas": {
            "plan": PLAN_SCHEMA,
            "simulation": SIMULATION_SCHEMA,
            "decode": DECODE_SCHEMA,
            "results": RESULTS_SCHEMA,
        },
    }
    outputs = {
        key: _output_record(results / name, rows=len(frames[key]))
        for key, name in RESULT_OUTPUT_PATHS.items()
    }
    completion = {
        "schema": RESULTS_SCHEMA,
        "status": "complete",
        "created_utc": _utc_now(),
        "contract": contract,
        "contract_sha256": _canonical_sha256(contract),
        "outputs": outputs,
        "unit_inventory": inventory,
    }
    _atomic_json(results / "RESULTS_COMPLETION.json", completion)
    verify_study(repo, root=root)
    return results / "RESULTS_COMPLETION.json"


def verify_study(
    repo_root: str | Path,
    *,
    root: str | Path = DEFAULT_STUDY_ROOT,
) -> dict[str, Any]:
    """Deep-verify the immutable plan, every unit, and compact result bundle."""

    repo = Path(repo_root).resolve()
    destination = study_root(repo, root)
    plan = load_plan(repo, root=root)
    results = destination / "results"
    completion_path = results / "RESULTS_COMPLETION.json"
    try:
        completion = json.loads(completion_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError("results completion is absent or unreadable") from error
    contract = completion.get("contract")
    if (
        completion.get("schema") != RESULTS_SCHEMA
        or completion.get("status") != "complete"
        or not isinstance(contract, Mapping)
        or completion.get("contract_sha256") != _canonical_sha256(contract)
    ):
        raise ValueError("results completion contract differs")
    expected_contract = {
        "study_id": STUDY_ID,
        "expected_units": DEFAULT_REPLICATES,
        "plan_sha256": str(plan["contract_sha256"]),
        "execution_units_sha256": _sha256_file(
            destination / "plan/execution_units.tsv"
        ),
        "estimand": ESTIMAND,
        "simulation_law": SIMULATION_LAW,
        "decoder_settings": DECODER_SETTINGS,
        "decoder_binary": dict(plan["contract"]["decoder_binary"]),
        "source_hashes": dict(plan["contract"]["source_hashes"]),
        "schemas": {
            "plan": PLAN_SCHEMA,
            "simulation": SIMULATION_SCHEMA,
            "decode": DECODE_SCHEMA,
            "results": RESULTS_SCHEMA,
        },
    }
    if contract != expected_contract:
        raise ValueError("results root contract differs from the frozen plan")
    (
        expected_scores,
        expected_endpoints,
        expected_sim,
        expected_decode,
        expected_inventory,
    ) = _aggregate_payload(repo, destination, plan)
    if completion.get("unit_inventory") != expected_inventory:
        raise ValueError("results unit receipt inventory differs")
    outputs = completion.get("outputs")
    if not isinstance(outputs, Mapping) or set(outputs) != set(RESULT_OUTPUT_PATHS):
        raise ValueError("results output inventory differs")
    expected_frames = {
        "endpoint_summary": expected_endpoints,
        "replicate_scores": expected_scores,
        "simulation_status": expected_sim,
        "decode_status": expected_decode,
    }
    for key, name in RESULT_OUTPUT_PATHS.items():
        path = results / name
        _validate_file_record(
            outputs[key], path, expected_rows=len(expected_frames[key])
        )
        observed = pd.read_csv(path, sep="\t", keep_default_na=False)
        pd.testing.assert_frame_equal(
            observed,
            expected_frames[key],
            check_dtype=False,
            check_exact=False,
            rtol=1e-11,
            atol=1e-12,
        )
    actual = {path.name for path in results.iterdir() if path.is_file()}
    expected_files = {"RESULTS_COMPLETION.json", *RESULT_OUTPUT_PATHS.values()}
    if actual != expected_files:
        raise ValueError("results directory contains missing or extra files")
    return dict(completion)


def main(argv: list[str] | None = None) -> int:
    """Dispatch an explicitly selected restartable workflow phase."""

    parser = argparse.ArgumentParser(
        description=(
            "AncientEurasia_9K19 single-origin neutral simulations conditioned "
            "only on strict present-Han census segregation"
        )
    )
    parser.add_argument(
        "phase",
        choices=("plan", "simulate", "decode", "status", "collect", "verify", "all"),
    )
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--study-root", type=Path, default=DEFAULT_STUDY_ROOT)
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument("--threads", type=int, default=DEFAULT_DECODE_THREADS)
    parser.add_argument("--slim-bin", type=Path)
    parser.add_argument("--decoder-bin", type=Path)
    parser.add_argument(
        "--unit-id",
        action="append",
        dest="unit_ids",
        help="run only this planned unit (repeatable; the plan remains exactly 100 units)",
    )
    parser.add_argument("--force-plan", action="store_true")
    args = parser.parse_args(argv)
    repo = args.repo_root.resolve()
    root = args.study_root
    _validate_workers(args.workers)
    if args.threads != DEFAULT_DECODE_THREADS:
        raise ValueError("--threads must remain 1 for the frozen decoder contract")
    if args.force_plan and args.phase != "plan":
        raise ValueError("--force-plan is accepted only for the plan phase")
    plan_file = _plan_path(study_root(repo, root))
    if args.phase in {"plan", "all"}:
        path = write_plan(
            repo,
            root=root,
            decoder_bin=args.decoder_bin,
            force=bool(args.force_plan),
        )
        print(path)
        if args.phase == "plan":
            return 0
    elif not plan_file.is_file():
        raise ValueError("study plan is absent; run the plan phase first")

    def report(label: str, frame: pd.DataFrame) -> None:
        counts = frame.groupby("status", dropna=False).size().to_dict()
        print(json.dumps({"phase": label, "status_counts": counts}, sort_keys=True))

    if args.phase in {"simulate", "all"}:
        frame = simulate_study(
            repo,
            root=root,
            slim_bin=args.slim_bin,
            workers=args.workers,
            unit_ids=args.unit_ids,
        )
        report("simulate", frame)
        incomplete = set(frame["status"].astype(str)).difference({"complete", "cached"})
        if incomplete:
            raise RuntimeError(
                f"simulation phase has incomplete statuses: {sorted(incomplete)}"
            )
    if args.phase in {"decode", "all"}:
        frame = decode_study(
            repo,
            root=root,
            decoder_bin=args.decoder_bin,
            workers=args.workers,
            threads=args.threads,
            unit_ids=args.unit_ids,
        )
        report("decode", frame)
        incomplete = set(frame["status"].astype(str)).difference({"complete", "cached"})
        if incomplete:
            raise RuntimeError(
                f"decode phase has incomplete statuses: {sorted(incomplete)}"
            )
    if args.phase == "status":
        report("simulate", simulation_status(repo, root=root))
        report("decode", decode_status(repo, root=root))
        return 0
    if args.phase == "collect":
        print(collect_results(repo, root=root))
        return 0
    if args.phase in {"verify", "all"}:
        verification = verify_study(repo, root=root)
        print(
            json.dumps(
                {
                    "status": verification["status"],
                    "contract_sha256": verification["contract_sha256"],
                    "units": len(verification["unit_inventory"]),
                },
                sort_keys=True,
            )
        )
    return 0


__all__ = [
    "BASE_SEED",
    "DECODER_SETTINGS",
    "DECODE_SCHEMA",
    "DEFAULT_REPLICATES",
    "DEFAULT_STUDY_ROOT",
    "ENDPOINT_COLUMNS",
    "ESTIMAND",
    "EXECUTION_COLUMNS",
    "PLAN_SCHEMA",
    "RESULTS_SCHEMA",
    "SCORE_COLUMNS",
    "SIMULATION_LAW",
    "SIMULATION_PARAMETERS",
    "SIMULATION_SCHEMA",
    "TMRCA_THRESHOLDS_YEARS",
    "build_neutral_event_schedule",
    "collect_results",
    "decode_status",
    "decode_study",
    "load_plan",
    "main",
    "simulate_study",
    "simulation_status",
    "verify_study",
    "write_plan",
]
