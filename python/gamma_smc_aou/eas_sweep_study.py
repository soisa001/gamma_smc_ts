"""Phased orchestration for the EAS selected-sweep simulation studies.

The low-level demographic and selected-site contracts live in
``eas_sweep_models``.  This module owns execution, integrity-checked caching,
Gamma-SMC decoding, aggregation, and provenance.  Raw tree sequences and
decoder intermediates live below ``work/``; compact tables and figures live
below ``results/``.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import multiprocessing
import os
import platform
import queue
import shutil
import signal
import subprocess
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path
from time import perf_counter, sleep
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
import stdpopsim
import tskit

from .decoder import run_within_decoder
from .eas_sweep_analysis import (
    DEFAULT_THRESHOLDS_YEARS,
    build_diploid_pair_table,
    gamma_summaries_to_long,
    internal_genotype_contrasts,
    plot_ancient_eurasia_introgression_timeline,
    plot_eas_ne_origin_design,
    plot_simulation_status_grid,
    tree_truth_profiles,
    write_cross_spec_study_plots,
    write_gamma_pair_files,
    write_publication_plots,
)
from .eas_sweep_models import (
    DE_NOVO_DEFAULT_LOWER_FREQUENCY,
    DE_NOVO_DEFAULT_UPPER_FREQUENCY,
    DE_NOVO_TARGET_FREQUENCY,
    DEFAULT_ARCHAIC_SOURCE_MINIMUM_FREQUENCY,
    DOMINANCE_COEFFICIENT,
    EAS_POPULATION,
    ARCHAIC_POPULATION,
    FOCAL_POSITION_BP,
    GENERATION_TIME_YEARS,
    HAN_SPLIT_GENERATIONS,
    INTROGRESSION_PULSE_GENERATIONS,
    INTROGRESSION_TARGET_POPULATION,
    INTROGRESSION_TARGET_FREQUENCIES,
    MUTATION_RATE,
    RECOMBINATION_RATE,
    SELECTION_COEFFICIENTS,
    SEQUENCE_LENGTH_BP,
    build_de_novo_origin_age_grid,
    build_eas_demography_models,
    build_introgression_sweep_spec,
    build_no_introgression_sweep_spec,
    load_phlash_eas_npz,
)


EAS_RESOURCE_SHA256 = "a5bdfd843629d48050cd56a84da4d70e356d0339017d3c1f67be646d727db928"
NO_INTROGRESSION_DIR = "no_introgression_EAS_sim"
INTROGRESSION_DIR = "introgression_EAS_sim"
BASE_SEED = 20240523
DEFAULT_SAMPLE_DIPLOIDS = 100
DEFAULT_THREADS = 4
DEFAULT_SIMULATION_WORKERS = 4
DEFAULT_OUTPUT_STRIDE_BP = 10_000
DEFAULT_CACHE_SIZE_BP = 10_000
DEFAULT_SLIM_SCALING_FACTOR = 10.0
DEFAULT_SLIM_BURN_IN = 0.1
DEFAULT_MAX_EXTERNAL_DRAWS = 20
DEFAULT_MIN_HOMOZYGOUS_DIPLOIDS = 2
DEFAULT_INTROGRESSION_MAXIMUM_AF = 0.95
INTROGRESSION_ORIGIN_GENERATIONS = 2_400.0
SCHEMA_VERSION = "gamma-smc.eas-sweep-study/v1"
EXTERNAL_DRAW_STATE_SCHEMA = "gamma-smc.eas-sweep-external-draw-state/v1"
EXTERNAL_DRAW_COLUMNS = (
    "external_draw_zero_based",
    "seed",
    "elapsed_seconds",
    "accepted",
    "achieved_allele_frequency",
    "achieved_allele_frequency_scope",
    "population_frequency_condition_lower",
    "population_frequency_condition_upper",
    "n_sample_diploids",
    "n_hom_ref",
    "n_heterozygous",
    "n_hom_alt",
    "focal_site_id",
    "focal_mutation_count",
    "focal_mutation_node_population",
    "focal_mutation_declared_origin_population",
    "focal_mutation_node_population_is_origin_evidence",
    "focal_mutation_origin_evidence",
    "n_segregating_nucleotide_snps",
    "vcf_allele_encoding",
)
REPRESENTATIVE_SPECIFICATIONS = {
    "no_introgression": "no_intro_s0p010_median_age1289g",
    "introgression": "intro_s0p005_af50",
}


@dataclass(frozen=True)
class StudyRuntime:
    """Execution parameters shared by all requested specifications."""

    repo_root: str
    sample_diploids: int = DEFAULT_SAMPLE_DIPLOIDS
    sequence_length_bp: int = SEQUENCE_LENGTH_BP
    focal_position_bp: int = FOCAL_POSITION_BP
    mutation_rate: float = MUTATION_RATE
    recombination_rate: float = RECOMBINATION_RATE
    generation_time_years: float = GENERATION_TIME_YEARS
    output_stride_bp: int = DEFAULT_OUTPUT_STRIDE_BP
    cache_size_bp: int = DEFAULT_CACHE_SIZE_BP
    threads: int = DEFAULT_THREADS
    simulation_workers: int = DEFAULT_SIMULATION_WORKERS
    slim_scaling_factor: float = DEFAULT_SLIM_SCALING_FACTOR
    slim_burn_in: float = DEFAULT_SLIM_BURN_IN
    max_external_draws: int = DEFAULT_MAX_EXTERNAL_DRAWS
    max_specification_seconds: float = 0.0
    minimum_homozygous_diploids: int = DEFAULT_MIN_HOMOZYGOUS_DIPLOIDS
    base_seed: int = BASE_SEED

    def validate(self) -> None:
        if self.sequence_length_bp != SEQUENCE_LENGTH_BP:
            raise ValueError("the current focal-contig contract is exactly 10 Mb")
        if self.focal_position_bp != FOCAL_POSITION_BP:
            raise ValueError("the current selected-site contract is fixed at 5 Mb")
        for name in (
            "sample_diploids",
            "output_stride_bp",
            "cache_size_bp",
            "threads",
            "simulation_workers",
            "max_external_draws",
            "minimum_homozygous_diploids",
        ):
            if int(getattr(self, name)) < 1:
                raise ValueError(f"{name} must be positive")
        if not 0 < self.slim_scaling_factor:
            raise ValueError("slim_scaling_factor must be positive")
        if self.slim_burn_in < 0:
            raise ValueError("slim_burn_in must be non-negative")
        if self.max_specification_seconds < 0:
            raise ValueError("max_specification_seconds must be non-negative")


def _canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _contract_sha256(value: Any) -> str:
    return hashlib.sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def _sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _atomic_json(path: str | Path, payload: Any) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + f".tmp.{os.getpid()}")
    try:
        temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)
    return destination


def _atomic_frame(path: str | Path, frame: pd.DataFrame) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + f".tmp.{os.getpid()}")
    try:
        frame.to_csv(
            temporary,
            sep="\t",
            index=False,
            compression="gzip" if destination.suffix == ".gz" else None,
        )
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)
    return destination


def _status_column_last(frame: pd.DataFrame) -> pd.DataFrame:
    """Put the required nonempty status field last in tracked status TSVs."""
    if "status" not in frame.columns:
        raise ValueError("status table is missing the status column")
    columns = [column for column in frame.columns if column != "status"]
    return frame.loc[:, [*columns, "status"]]


def _atomic_copy(source: str | Path, destination: str | Path) -> Path:
    source = Path(source)
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + f".tmp.{os.getpid()}")
    try:
        shutil.copy2(source, temporary)
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)
    return destination


def _stable_seed(base_seed: int, specification_id: str, draw: int) -> int:
    token = f"{base_seed}:{specification_id}:{draw}".encode("utf-8")
    offset = int.from_bytes(hashlib.sha256(token).digest()[:8], "big")
    return int(1 + offset % (2**31 - 2))


def _slug_s(selection_coefficient: float) -> str:
    return f"{selection_coefficient:.3f}".replace(".", "p")


def _scaled_slim_time(time_generations: float, scaling_factor: float) -> float:
    """Mirror stdpopsim 0.3's integral-tick rounding for extended events."""
    return round(float(time_generations) / float(scaling_factor)) * float(
        scaling_factor
    )


def _scaled_demographic_time(time_generations: float, scaling_factor: float) -> float:
    """Return the positive-time tick used by generated demographic events.

    stdpopsim 0.3 writes catalog demographic times to ``_T`` without applying
    its Python extended-event rounding. The generated SLiM ``time_to_tick``
    function converts the positive quotient with Eidos ``asInteger()``, which
    truncates it. Keep this distinct from :func:`_scaled_slim_time` so a split
    and a fitness-population transition cannot silently land on adjacent ticks.
    """

    return math.floor(float(time_generations) / float(scaling_factor)) * float(
        scaling_factor
    )


def _portable_runtime_record(runtime: StudyRuntime) -> dict[str, Any]:
    """Return the tracked design record without a workstation-specific path."""
    record = asdict(runtime)
    record["repo_root"] = "."
    return record


def _simulation_runtime_record(runtime: StudyRuntime) -> dict[str, Any]:
    """Return only result-affecting simulation and truth-profile settings.

    Worker counts and the external-search ceiling do not change an already
    accepted seed-ordered draw, so they deliberately do not invalidate it.
    """
    return {
        "sample_diploids": runtime.sample_diploids,
        "sequence_length_bp": runtime.sequence_length_bp,
        "focal_position_bp": runtime.focal_position_bp,
        "mutation_rate": runtime.mutation_rate,
        "recombination_rate": runtime.recombination_rate,
        "generation_time_years": runtime.generation_time_years,
        "output_stride_bp": runtime.output_stride_bp,
        "slim_scaling_factor": runtime.slim_scaling_factor,
        "slim_burn_in": runtime.slim_burn_in,
        "minimum_homozygous_diploids": runtime.minimum_homozygous_diploids,
        "base_seed": runtime.base_seed,
    }


def _software_record(*, slim_path: str | Path | None = None) -> dict[str, Any]:
    record: dict[str, Any] = {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "stdpopsim": stdpopsim.__version__,
        "tskit": tskit.__version__,
    }
    if slim_path is not None:
        executable = Path(slim_path).resolve()
        record["slim_path"] = str(executable)
        record["slim_sha256"] = _sha256(executable)
        completed = subprocess.run(
            [str(executable), "-v"], text=True, capture_output=True, check=False
        )
        record["slim_version_output"] = (completed.stdout + completed.stderr).strip()
    return record


def _study_dir(repo_root: str | Path, scenario: str) -> Path:
    if scenario == "no_introgression":
        return Path(repo_root) / NO_INTROGRESSION_DIR
    if scenario == "introgression":
        return Path(repo_root) / INTROGRESSION_DIR
    raise ValueError(f"unknown scenario {scenario!r}")


def _no_introgression_specifications(runtime: StudyRuntime) -> list[dict[str, Any]]:
    study_dir = _study_dir(runtime.repo_root, "no_introgression")
    artifact = load_phlash_eas_npz(
        study_dir / "resources" / "EAS.npz",
        expected_sha256=EAS_RESOURCE_SHA256,
    )
    demographies = build_eas_demography_models(artifact)
    ages = build_de_novo_origin_age_grid(
        demographies, slim_scaling_factor=runtime.slim_scaling_factor
    )
    records = []
    for estimate in ages:
        specification_id = (
            f"no_intro_s{_slug_s(estimate.selection_coefficient)}_"
            f"{estimate.trajectory_label}_age{int(estimate.age_generations)}g"
        )
        demographic = demographies[estimate.trajectory_label]
        records.append(
            {
                "scenario": "no_introgression_de_novo",
                "study_type": "no_introgression",
                "specification_id": specification_id,
                "selection_coefficient": estimate.selection_coefficient,
                "demography_curve": estimate.trajectory_label,
                "demography_quantile": demographic.quantile,
                "target_allele_frequency": DE_NOVO_TARGET_FREQUENCY,
                "lower_allele_frequency": DE_NOVO_DEFAULT_LOWER_FREQUENCY,
                "upper_allele_frequency": DE_NOVO_DEFAULT_UPPER_FREQUENCY,
                "origin_age_generations": estimate.age_generations,
                "origin_age_years": estimate.age_years,
                "origin_age_method": estimate.to_record()["method"],
                "origin_age_is_exact": False,
                "slim_scaled_origin_age_generations": _scaled_slim_time(
                    estimate.age_generations, runtime.slim_scaling_factor
                ),
                "slim_scaled_post_origin_survival_check_generations": (
                    _scaled_slim_time(
                        estimate.age_generations, runtime.slim_scaling_factor
                    )
                    - runtime.slim_scaling_factor
                ),
                "ne_at_origin": estimate.ne_at_origin,
                "present_ne": float(demographic.ne[0]),
                "theta_for_gamma_smc": 4
                * float(demographic.ne[0])
                * runtime.mutation_rate,
                "deterministic_initial_frequency": estimate.initial_frequency,
                "deterministic_initial_frequency_slim_scaling_factor": (
                    estimate.slim_scaling_factor
                ),
                "deterministic_initial_frequency_semantics": (
                    "one copy in the Q-scaled diploid population, approximated "
                    "as Q/(2*Ne) before integral population-size rounding"
                ),
                "deterministic_fixed_point_residual_generations": (
                    estimate.fixed_point_residual_generations
                ),
                "deterministic_boundary_limited": estimate.boundary_limited,
                "eligible_for_simulation": True,
                "eligibility_reason": "conditioned_de_novo_frequency_interval",
                "model_record": demographic.provenance_record(),
                "origin_record": estimate.to_record(),
            }
        )
    return records


def _required_introgressed_frequency(
    final_frequency: float,
    selection_coefficient: float,
    *,
    dominance_coefficient: float = DOMINANCE_COEFFICIENT,
) -> float:
    """Deterministic rare-allele approximation for the pulse-time frequency."""
    duration = INTROGRESSION_PULSE_GENERATIONS
    logit_final = math.log(final_frequency / (1 - final_frequency))
    logit_initial = logit_final - (
        dominance_coefficient * selection_coefficient * duration
    )
    return 1 / (1 + math.exp(-logit_initial))


def _introgression_specifications(runtime: StudyRuntime) -> list[dict[str, Any]]:
    records = []
    realized_han_split = _scaled_demographic_time(
        HAN_SPLIT_GENERATIONS, runtime.slim_scaling_factor
    )
    for selection_coefficient in SELECTION_COEFFICIENTS:
        for target in INTROGRESSION_TARGET_FREQUENCIES:
            required = _required_introgressed_frequency(target, selection_coefficient)
            deterministic_feasible = required <= 0.0296
            specification_id = (
                f"intro_s{_slug_s(selection_coefficient)}_af{round(target * 100):02d}"
            )
            records.append(
                {
                    "scenario": "neanderthal_introgression_standing_variation",
                    "study_type": "introgression",
                    "specification_id": specification_id,
                    "selection_coefficient": selection_coefficient,
                    "demography_curve": None,
                    "target_allele_frequency": target,
                    "lower_allele_frequency": target,
                    "upper_allele_frequency": DEFAULT_INTROGRESSION_MAXIMUM_AF,
                    "origin_age_generations": INTROGRESSION_ORIGIN_GENERATIONS,
                    "origin_age_years": (
                        INTROGRESSION_ORIGIN_GENERATIONS * runtime.generation_time_years
                    ),
                    "source_minimum_frequency": (
                        DEFAULT_ARCHAIC_SOURCE_MINIMUM_FREQUENCY
                    ),
                    "pulse_proportion": 0.0296,
                    "pulse_generations_ago": INTROGRESSION_PULSE_GENERATIONS,
                    "han_split_generations_ago": HAN_SPLIT_GENERATIONS,
                    "requested_han_split_generations_ago": HAN_SPLIT_GENERATIONS,
                    "realized_han_split_generations_ago": realized_han_split,
                    "slim_scaled_origin_generations_ago": _scaled_slim_time(
                        INTROGRESSION_ORIGIN_GENERATIONS,
                        runtime.slim_scaling_factor,
                    ),
                    "slim_scaled_source_check_generations_ago": _scaled_slim_time(
                        2_273.0, runtime.slim_scaling_factor
                    ),
                    "slim_scaled_pulse_generations_ago": _scaled_slim_time(
                        INTROGRESSION_PULSE_GENERATIONS,
                        runtime.slim_scaling_factor,
                    ),
                    "slim_scaled_han_split_generations_ago": _scaled_slim_time(
                        realized_han_split, runtime.slim_scaling_factor
                    ),
                    "slim_demographic_han_split_generations_ago": realized_han_split,
                    "slim_selection_transition_generations_ago": realized_han_split,
                    "slim_demographic_time_rule": (
                        "floor(time/Q)*Q for positive catalog demographic events"
                    ),
                    "deterministic_required_pulse_frequency": required,
                    "deterministic_feasible_from_fixed_archaic_source": (
                        deterministic_feasible
                    ),
                    "feasibility_interpretation": (
                        "screening approximation only; drift and conditioning remain stochastic"
                    ),
                    "eligible_for_simulation": True,
                    "eligibility_reason": (
                        "conditioned_stochastic_tail"
                        if not deterministic_feasible
                        else "deterministically_reachable_from_maximum_pulse"
                    ),
                    "present_ne": 6300.0,
                    "theta_for_gamma_smc": 4 * 6300.0 * runtime.mutation_rate,
                }
            )
    return records


def build_specifications(
    runtime: StudyRuntime, scenarios: Sequence[str]
) -> dict[str, list[dict[str, Any]]]:
    runtime.validate()
    result: dict[str, list[dict[str, Any]]] = {}
    for scenario in scenarios:
        if scenario == "no_introgression":
            result[scenario] = _no_introgression_specifications(runtime)
        elif scenario == "introgression":
            result[scenario] = _introgression_specifications(runtime)
        else:
            raise ValueError(f"unknown scenario {scenario!r}")
    return result


def write_study_plans(
    runtime: StudyRuntime, scenarios: Sequence[str]
) -> dict[str, list[dict[str, Any]]]:
    specifications = build_specifications(runtime, scenarios)
    for scenario, records in specifications.items():
        study_dir = _study_dir(runtime.repo_root, scenario)
        results_dir = study_dir / "results"
        results_dir.mkdir(parents=True, exist_ok=True)
        serializable = [
            {k: v for k, v in row.items() if k not in {"model_record", "origin_record"}}
            for row in records
        ]
        _atomic_frame(study_dir / "specifications.tsv", pd.DataFrame(serializable))
        design = {
            "schema": SCHEMA_VERSION,
            "phase": "plan",
            "scenario": scenario,
            "runtime": _portable_runtime_record(runtime),
            "tmrca_thresholds_years": list(DEFAULT_THRESHOLDS_YEARS),
            "selection_coefficients": list(SELECTION_COEFFICIENTS),
            "rng": {
                "base_seed": runtime.base_seed,
                "derivation": "sha256(base_seed:specification_id:external_draw)",
            },
            "frequency_semantics": (
                "SLiM conditions population AF; accepted sample AF and genotype counts "
                "are separately validated and reported"
            ),
            "auc_semantics": (
                "normalized within-genotype TMRCA-CDF area, not ROC AUC; no neutral null"
            ),
            "scaling_caveat": (
                "SLiM scaling is an exploratory approximation and is recorded per run"
            ),
            "slim_time_rounding": (
                "stdpopsim rounds extended events to the nearest scaled tick, "
                "while its generated demographic time_to_tick truncates positive "
                "catalog times; specifications.tsv records both requested and "
                "realized study-critical times"
            ),
            "specification_count": len(records),
        }
        if scenario == "no_introgression":
            design["eas_demography_resource"] = {
                "path": "resources/EAS.npz",
                "sha256": EAS_RESOURCE_SHA256,
                "quantiles": [0.025, 0.5, 0.975],
                "interval_semantics": "empirical_pointwise_bootstrap_quantiles",
            }
        else:
            realized_han_split = _scaled_demographic_time(
                HAN_SPLIT_GENERATIONS, runtime.slim_scaling_factor
            )
            design["scaled_introgression_event_times_generations_ago"] = {
                "archaic_mutation_origin": _scaled_slim_time(
                    INTROGRESSION_ORIGIN_GENERATIONS, runtime.slim_scaling_factor
                ),
                "requested_pre_pulse_source_check": _scaled_slim_time(
                    2_273.0, runtime.slim_scaling_factor
                ),
                "introgression_pulse_and_selection_onset": _scaled_slim_time(
                    INTROGRESSION_PULSE_GENERATIONS, runtime.slim_scaling_factor
                ),
                "requested_catalog_han_split": HAN_SPLIT_GENERATIONS,
                "realized_demographic_han_split": realized_han_split,
                "selection_population_transition": realized_han_split,
            }
            design["scaled_time_rules"] = {
                "extended_events": "round(time/Q)*Q using Python round",
                "positive_catalog_demographic_events": (
                    "floor(time/Q)*Q through generated Eidos asInteger"
                ),
                "han_split_alignment_validated": True,
            }
            design["scaled_time_collision"] = (
                "the requested 2,273-generation source check and 2,272-generation "
                "pulse both execute at 2,270 generations under Q=10"
            )
        _atomic_json(study_dir / "study_design.json", design)
        if scenario == "no_introgression":
            detailed = {
                row["specification_id"]: {
                    "demography": copy.deepcopy(row["model_record"]),
                    "origin": copy.deepcopy(row["origin_record"]),
                }
                for row in records
            }
            for payload in detailed.values():
                payload["demography"]["source"]["source_path"] = "resources/EAS.npz"
            _atomic_json(study_dir / "origin_age_contracts.json", detailed)
            artifact = load_phlash_eas_npz(
                study_dir / "resources" / "EAS.npz",
                expected_sha256=EAS_RESOURCE_SHA256,
            )
            demographies = build_eas_demography_models(artifact)
            trajectories = pd.concat(
                [
                    pd.DataFrame(
                        {
                            "demography_curve": label,
                            "time_generations": model.time_generations,
                            "effective_population_size": model.ne,
                        }
                    )
                    for label, model in demographies.items()
                ],
                ignore_index=True,
            )
            origins = pd.DataFrame(
                [
                    {
                        "demography_curve": row["demography_curve"],
                        "selection_coefficient": row["selection_coefficient"],
                        "origin_age_generations": row["origin_age_generations"],
                        "ne_at_origin": row["ne_at_origin"],
                    }
                    for row in records
                ]
            )
            plot_eas_ne_origin_design(
                trajectories,
                origins,
                results_dir / "figures" / "design" / "eas_demography_and_origins",
                generation_time_years=runtime.generation_time_years,
            )
        else:
            plot_ancient_eurasia_introgression_timeline(
                results_dir / "figures" / "design" / "ancient_eurasia_timeline",
                introgression_pulse_generations=INTROGRESSION_PULSE_GENERATIONS,
                realized_introgression_pulse_generations=_scaled_slim_time(
                    INTROGRESSION_PULSE_GENERATIONS, runtime.slim_scaling_factor
                ),
                han_split_generations=HAN_SPLIT_GENERATIONS,
                realized_han_split_generations=_scaled_demographic_time(
                    HAN_SPLIT_GENERATIONS, runtime.slim_scaling_factor
                ),
                generation_time_years=runtime.generation_time_years,
            )
    return specifications


def _selected_site(ts: tskit.TreeSequence, focal_position: float):
    sites = [
        site
        for site in ts.sites()
        if math.isclose(site.position, focal_position, rel_tol=0.0, abs_tol=1e-9)
    ]
    if len(sites) != 1:
        raise ValueError(
            f"expected one selected site at {focal_position:g}; found {len(sites)}"
        )
    if len(sites[0].mutations) != 1:
        raise ValueError("the selected site must contain one nonrecurrent mutation")
    return sites[0]


def _mutation_population_name(ts: tskit.TreeSequence, mutation) -> str | None:
    population_id = ts.node(mutation.node).population
    if population_id == tskit.NULL:
        return None
    metadata = ts.population(population_id).metadata
    if isinstance(metadata, Mapping):
        return metadata.get("name")
    return None


def _validate_accepted_tree(
    ts: tskit.TreeSequence,
    specification: Mapping[str, Any],
    runtime: StudyRuntime,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    if not math.isclose(ts.sequence_length, runtime.sequence_length_bp):
        raise ValueError("tree sequence length does not match the 10-Mb contract")
    segregating_snps = 0
    for variant in ts.variants():
        observed = np.unique(variant.genotypes[variant.genotypes >= 0])
        alleles = [allele for allele in variant.alleles if allele is not None]
        nucleotide_alleles = all(
            len(allele) == 1 and allele in {"A", "C", "G", "T"} for allele in alleles
        )
        if len(observed) > 1 and nucleotide_alleles:
            segregating_snps += 1
        elif len(observed) > 1:
            raise ValueError(
                "segregating variants must use single-nucleotide A/C/G/T alleles "
                "before Gamma-SMC VCF streaming"
            )
    if segregating_snps == 0:
        raise ValueError("accepted tree contains no segregating nucleotide SNPs")
    site = _selected_site(ts, runtime.focal_position_bp)
    pair_table = build_diploid_pair_table(ts, runtime.focal_position_bp)
    if len(pair_table) != runtime.sample_diploids:
        raise ValueError(
            f"expected {runtime.sample_diploids} sampled diploids; found {len(pair_table)}"
        )
    counts = pair_table["focal_selected_allele_count"].to_numpy(dtype=int)
    sample_af = float(counts.sum() / (2 * len(counts)))
    class_counts = pair_table["genotype_class"].value_counts().to_dict()
    minimum = runtime.minimum_homozygous_diploids
    accepted = (
        class_counts.get("hom_ref", 0) >= minimum
        and class_counts.get("hom_alt", 0) >= minimum
        and class_counts.get("heterozygous", 0) >= 1
        and sample_af >= float(specification["lower_allele_frequency"])
        and sample_af <= float(specification["upper_allele_frequency"])
    )
    mutation = site.mutations[0]
    mutation_population = _mutation_population_name(ts, mutation)
    expected_origin_population = (
        EAS_POPULATION
        if specification["study_type"] == "no_introgression"
        else ARCHAIC_POPULATION
    )
    validation = {
        "accepted": bool(accepted),
        "achieved_allele_frequency": sample_af,
        "achieved_allele_frequency_scope": "sample_estimate",
        "population_frequency_condition_lower": float(
            specification["lower_allele_frequency"]
        ),
        "population_frequency_condition_upper": float(
            specification["upper_allele_frequency"]
        ),
        "n_sample_diploids": len(pair_table),
        "n_hom_ref": int(class_counts.get("hom_ref", 0)),
        "n_heterozygous": int(class_counts.get("heterozygous", 0)),
        "n_hom_alt": int(class_counts.get("hom_alt", 0)),
        "focal_site_id": int(site.id),
        "focal_mutation_count": len(site.mutations),
        "focal_mutation_node_population": mutation_population,
        "focal_mutation_declared_origin_population": expected_origin_population,
        "focal_mutation_node_population_is_origin_evidence": False,
        "focal_mutation_origin_evidence": "serialized_forward_event_contract",
        "n_segregating_nucleotide_snps": segregating_snps,
        "vcf_allele_encoding": "single_nucleotide_ACGT",
    }
    return pair_table, validation


def _simulation_contract(
    specification: Mapping[str, Any], runtime: StudyRuntime, slim_path: str | Path
) -> dict[str, Any]:
    return {
        "schema": SCHEMA_VERSION,
        "phase": "simulate",
        "specification": {
            key: value
            for key, value in specification.items()
            if key not in {"model_record", "origin_record"}
        },
        "runtime": _simulation_runtime_record(runtime),
        "software": _software_record(slim_path=slim_path),
        "sample_manifest_contract": {
            "population": (
                EAS_POPULATION
                if specification["study_type"] == "no_introgression"
                else INTROGRESSION_TARGET_POPULATION
            ),
            "diploid_individuals": runtime.sample_diploids,
            "ploidy": 2,
        },
        "vcf_allele_encoding": "pyslim_generated_nucleotides",
    }


def _validated_simulation_contract(
    specification: Mapping[str, Any],
    runtime: StudyRuntime,
    slim_path: str | Path,
    model_record: Mapping[str, Any],
) -> dict[str, Any]:
    expected_origin_population = (
        EAS_POPULATION
        if specification["study_type"] == "no_introgression"
        else ARCHAIC_POPULATION
    )
    declared_origin_population = (
        model_record.get("population")
        if specification["study_type"] == "no_introgression"
        else model_record.get("mutation_origin", {}).get("population")
    )
    if declared_origin_population != expected_origin_population:
        raise ValueError(
            "selected-site event contract has the wrong mutation origin: "
            f"{declared_origin_population!r} != {expected_origin_population!r}"
        )
    contract = _simulation_contract(specification, runtime, slim_path)
    contract["selected_site_model"] = dict(model_record)
    contract["mutation_origin_validation"] = {
        "declared_population": declared_origin_population,
        "expected_population": expected_origin_population,
        "validated": True,
        "evidence": "serialized_forward_event_contract",
        "tree_node_population_is_origin_evidence": False,
    }
    return contract


def _valid_simulation_cache(
    specification_dir: Path,
    contract: Mapping[str, Any],
    specification: Mapping[str, Any],
    runtime: StudyRuntime,
) -> dict[str, Any] | None:
    completion_path = specification_dir / "simulation_complete.json"
    tree_path = specification_dir / "selected.trees"
    pair_path = specification_dir / "sample_manifest.tsv"
    truth_path = specification_dir / "truth_profiles.tsv.gz"
    contrast_path = specification_dir / "truth_internal_contrasts.tsv.gz"
    required_paths = (
        completion_path,
        tree_path,
        pair_path,
        truth_path,
        contrast_path,
    )
    if not all(path.is_file() for path in required_paths):
        return None
    try:
        completion = json.loads(completion_path.read_text(encoding="utf-8"))
        if completion["contract_sha256"] != _contract_sha256(contract):
            return None
        if completion["tree_sha256"] != _sha256(tree_path):
            raise ValueError("cached selected.trees checksum mismatch")
        if completion["sample_manifest_sha256"] != _sha256(pair_path):
            raise ValueError("cached sample manifest checksum mismatch")
        if completion["truth_profiles_sha256"] != _sha256(truth_path):
            raise ValueError("cached truth profile checksum mismatch")
        if completion["truth_internal_contrasts_sha256"] != _sha256(contrast_path):
            raise ValueError("cached truth contrast checksum mismatch")
        for pair_record in completion["pair_files"].values():
            # Completion records may be validated from Linux/WSL after a
            # Windows SLiM run. Treat the stored value as a portable relative
            # path; older records used ``os.path.relpath`` and can therefore
            # contain Windows separators.
            relative_pair_path = str(pair_record["path"]).replace("\\", "/")
            cached_pair_path = specification_dir.joinpath(
                *Path(relative_pair_path).parts
            )
            if not cached_pair_path.is_file():
                raise ValueError("cached Gamma-SMC pair file is missing")
            if pair_record["sha256"] != _sha256(cached_pair_path):
                raise ValueError("cached Gamma-SMC pair checksum mismatch")
        if pd.read_csv(truth_path, sep="\t").empty:
            raise ValueError("cached truth profile is empty")
        if pd.read_csv(contrast_path, sep="\t").empty:
            raise ValueError("cached truth contrast is empty")
        ts = tskit.load(tree_path)
        pair_table, validation = _validate_accepted_tree(ts, specification, runtime)
        if not validation["accepted"]:
            raise ValueError("cached tree no longer passes acceptance checks")
        if len(pair_table) != completion["validation"]["n_sample_diploids"]:
            raise ValueError("cached sample manifest count mismatch")
        return completion
    except (KeyError, OSError, ValueError, json.JSONDecodeError, tskit.FileFormatError):
        return None


def _build_simulation_objects(
    specification: Mapping[str, Any], runtime: StudyRuntime
) -> tuple[Any, Any, dict[str, int], dict[str, Any]]:
    if specification["study_type"] == "no_introgression":
        study_dir = _study_dir(runtime.repo_root, "no_introgression")
        artifact = load_phlash_eas_npz(
            study_dir / "resources" / "EAS.npz",
            expected_sha256=EAS_RESOURCE_SHA256,
        )
        demographies = build_eas_demography_models(artifact)
        demographic = demographies[str(specification["demography_curve"])]
        sweep = build_no_introgression_sweep_spec(
            origin_age_generations=float(specification["origin_age_generations"]),
            selection_coefficient=float(specification["selection_coefficient"]),
            target_frequency=float(specification["target_allele_frequency"]),
            lower_frequency=float(specification["lower_allele_frequency"]),
            upper_frequency=float(specification["upper_allele_frequency"]),
        )
        return (
            demographic.stdpopsim_model,
            sweep,
            {EAS_POPULATION: runtime.sample_diploids},
            sweep.provenance_record(),
        )
    sweep = build_introgression_sweep_spec(
        selection_coefficient=float(specification["selection_coefficient"]),
        minimum_han_frequency=float(specification["lower_allele_frequency"]),
        maximum_han_frequency=float(specification["upper_allele_frequency"]),
        mutation_age_generations=float(specification["origin_age_generations"]),
        source_minimum_frequency=float(specification["source_minimum_frequency"]),
        realized_han_split_generations=float(
            specification["realized_han_split_generations_ago"]
        ),
    )
    return (
        sweep.demographic_model,
        sweep,
        {INTROGRESSION_TARGET_POPULATION: runtime.sample_diploids},
        sweep.provenance_record(),
    )


def _current_selected_site_model_record(
    specification: Mapping[str, Any],
) -> dict[str, Any]:
    """Rebuild the result-affecting selected-site event record without SLiM."""
    if specification["study_type"] == "no_introgression":
        sweep = build_no_introgression_sweep_spec(
            origin_age_generations=float(specification["origin_age_generations"]),
            selection_coefficient=float(specification["selection_coefficient"]),
            target_frequency=float(specification["target_allele_frequency"]),
            lower_frequency=float(specification["lower_allele_frequency"]),
            upper_frequency=float(specification["upper_allele_frequency"]),
        )
    else:
        sweep = build_introgression_sweep_spec(
            selection_coefficient=float(specification["selection_coefficient"]),
            minimum_han_frequency=float(specification["lower_allele_frequency"]),
            maximum_han_frequency=float(specification["upper_allele_frequency"]),
            mutation_age_generations=float(specification["origin_age_generations"]),
            source_minimum_frequency=float(specification["source_minimum_frequency"]),
            realized_han_split_generations=float(
                specification["realized_han_split_generations_ago"]
            ),
        )
    return dict(sweep.provenance_record())


def _recorded_simulation_contract_if_current(
    specification_dir: Path,
    specification: Mapping[str, Any],
    runtime: StudyRuntime,
) -> dict[str, Any] | None:
    """Return a recorded contract only when all current result inputs agree."""
    contract_path = specification_dir / "simulation_contract.json"
    if not contract_path.is_file():
        return None
    try:
        contract = json.loads(contract_path.read_text(encoding="utf-8"))
        expected_specification = {
            key: value
            for key, value in specification.items()
            if key not in {"model_record", "origin_record"}
        }
        expected_population = (
            EAS_POPULATION
            if specification["study_type"] == "no_introgression"
            else INTROGRESSION_TARGET_POPULATION
        )
        expected_origin = (
            EAS_POPULATION
            if specification["study_type"] == "no_introgression"
            else ARCHAIC_POPULATION
        )
        expected_manifest = {
            "population": expected_population,
            "diploid_individuals": runtime.sample_diploids,
            "ploidy": 2,
        }
        origin_validation = contract["mutation_origin_validation"]
        if (
            contract.get("schema") != SCHEMA_VERSION
            or contract.get("phase") != "simulate"
            or _canonical_json(contract.get("specification"))
            != _canonical_json(expected_specification)
            or _canonical_json(contract.get("runtime"))
            != _canonical_json(_simulation_runtime_record(runtime))
            or _canonical_json(contract.get("sample_manifest_contract"))
            != _canonical_json(expected_manifest)
            or contract.get("vcf_allele_encoding") != "pyslim_generated_nucleotides"
            or _canonical_json(contract.get("selected_site_model"))
            != _canonical_json(_current_selected_site_model_record(specification))
            or origin_validation.get("declared_population") != expected_origin
            or origin_validation.get("expected_population") != expected_origin
            or origin_validation.get("validated") is not True
            or origin_validation.get("evidence") != "serialized_forward_event_contract"
            or origin_validation.get("tree_node_population_is_origin_evidence")
            is not False
            or not isinstance(contract.get("software"), dict)
            or not contract["software"].get("slim_sha256")
        ):
            return None
        return contract
    except (KeyError, OSError, TypeError, ValueError, json.JSONDecodeError):
        return None


def _valid_recorded_simulation_cache(
    specification_dir: Path,
    specification: Mapping[str, Any],
    runtime: StudyRuntime,
) -> dict[str, Any] | None:
    """Validate an on-disk simulation against current inputs and artifact hashes."""
    contract = _recorded_simulation_contract_if_current(
        specification_dir, specification, runtime
    )
    if contract is None:
        return None
    return _valid_simulation_cache(specification_dir, contract, specification, runtime)


def _new_external_draw_state(contract_sha256: str) -> dict[str, Any]:
    return {
        "schema": EXTERNAL_DRAW_STATE_SCHEMA,
        "simulation_contract_sha256": contract_sha256,
        "next_external_draw_zero_based": 0,
        "completed_external_sample_panel_draws": 0,
        "active_internal_conditioned_trajectory": None,
        "interrupted_internal_conditioned_trajectories": [],
    }


def _validated_external_draw_progress(
    specification_dir: Path,
    contract: Mapping[str, Any],
    specification: Mapping[str, Any],
    runtime: StudyRuntime,
) -> tuple[dict[str, Any], list[dict[str, Any]]] | None:
    """Read a contract-bound ledger without changing interruption state."""
    state_path = specification_dir / "external_draw_state.json"
    attempts_path = specification_dir / "external_draws.tsv"
    if not state_path.is_file():
        return None
    try:
        state = json.loads(state_path.read_text(encoding="utf-8"))
        if state.get("schema") != EXTERNAL_DRAW_STATE_SCHEMA or state.get(
            "simulation_contract_sha256"
        ) != _contract_sha256(contract):
            return None
        if attempts_path.is_file():
            attempts_frame = pd.read_csv(attempts_path, sep="\t")
            attempts = attempts_frame.to_dict(orient="records")
        else:
            attempts = []
        if int(state["completed_external_sample_panel_draws"]) != len(attempts):
            return None
        draw_indices: list[int] = []
        for attempt in attempts:
            draw = int(attempt["external_draw_zero_based"])
            seed = int(attempt["seed"])
            if seed != _stable_seed(
                runtime.base_seed, specification["specification_id"], draw
            ):
                return None
            draw_indices.append(draw)
        if len(draw_indices) != len(set(draw_indices)):
            return None
        next_draw = int(state["next_external_draw_zero_based"])
        if draw_indices and next_draw <= max(draw_indices):
            return None
        active = state.get("active_internal_conditioned_trajectory")
        if active is not None:
            active_draw = int(active["external_draw_zero_based"])
            if active_draw != next_draw or int(active["seed"]) != _stable_seed(
                runtime.base_seed,
                specification["specification_id"],
                active_draw,
            ):
                return None
        interruptions = state.get("interrupted_internal_conditioned_trajectories")
        if not isinstance(interruptions, list):
            return None
        return state, attempts
    except (
        KeyError,
        OSError,
        TypeError,
        ValueError,
        json.JSONDecodeError,
        pd.errors.EmptyDataError,
        pd.errors.ParserError,
    ):
        return None


def _resume_external_draw_progress(
    specification_dir: Path,
    contract: Mapping[str, Any],
    specification: Mapping[str, Any],
    runtime: StudyRuntime,
    *,
    allow_legacy_progress: bool = False,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Resume completed draws and advance past an interrupted internal run."""
    state_path = specification_dir / "external_draw_state.json"
    attempts_path = specification_dir / "external_draws.tsv"
    progress = _validated_external_draw_progress(
        specification_dir, contract, specification, runtime
    )
    if progress is None:
        if allow_legacy_progress and state_path.is_file():
            raise ValueError(
                "external draw state exists but is corrupt or incompatible with "
                "its unchanged simulation contract"
            )
        state = _new_external_draw_state(_contract_sha256(contract))
        attempts: list[dict[str, Any]] = []
        if allow_legacy_progress:
            if attempts_path.is_file():
                try:
                    attempts = pd.read_csv(attempts_path, sep="\t").to_dict(
                        orient="records"
                    )
                except (
                    OSError,
                    pd.errors.EmptyDataError,
                    pd.errors.ParserError,
                ) as error:
                    raise ValueError(
                        "legacy external draw ledger is unreadable"
                    ) from error
            completed_indices: set[int] = set()
            for attempt in attempts:
                try:
                    draw = int(attempt["external_draw_zero_based"])
                    seed = int(attempt["seed"])
                except (KeyError, TypeError, ValueError) as error:
                    raise ValueError(
                        "legacy external draw ledger lacks a valid draw/seed"
                    ) from error
                if draw in completed_indices or seed != _stable_seed(
                    runtime.base_seed, specification["specification_id"], draw
                ):
                    raise ValueError(
                        "legacy external draw ledger has duplicate draws or a seed "
                        "that disagrees with the current contract"
                    )
                completed_indices.add(draw)
            launched_indices: set[int] = set()
            for logfile in specification_dir.glob("slim_draw_*.csv"):
                try:
                    launched_indices.add(int(logfile.stem.rsplit("_", 1)[-1]))
                except ValueError as error:
                    raise ValueError(
                        f"unparseable legacy SLiM draw logfile: {logfile.name}"
                    ) from error
            interrupted_indices = sorted(launched_indices - completed_indices)
            for draw in interrupted_indices:
                state["interrupted_internal_conditioned_trajectories"].append(
                    {
                        "external_draw_zero_based": draw,
                        "seed": _stable_seed(
                            runtime.base_seed,
                            specification["specification_id"],
                            draw,
                        ),
                        "stage": "legacy_uncheckpointed_stdpopsim_conditioning",
                        "outcome": (
                            "legacy_uncheckpointed_stdpopsim_internal_"
                            "conditioned_trajectory"
                        ),
                        "external_sample_panel_completed": False,
                    }
                )
            observed_indices = completed_indices | launched_indices
            if observed_indices:
                state["next_external_draw_zero_based"] = max(observed_indices) + 1
        _atomic_frame(
            attempts_path,
            pd.DataFrame(attempts, columns=EXTERNAL_DRAW_COLUMNS),
        )
    else:
        state, attempts = progress
        active = state.get("active_internal_conditioned_trajectory")
        if active is not None:
            outcome = (
                "interrupted_stdpopsim_internal_conditioned_trajectory"
                if active.get("stage") == "stdpopsim_internal_conditioning"
                else "interrupted_before_external_sample_panel_checkpoint"
            )
            interrupted = {
                **active,
                "outcome": outcome,
                "external_sample_panel_completed": False,
            }
            state["interrupted_internal_conditioned_trajectories"].append(interrupted)
            state["next_external_draw_zero_based"] = (
                int(active["external_draw_zero_based"]) + 1
            )
            state["active_internal_conditioned_trajectory"] = None
    state["completed_external_sample_panel_draws"] = len(attempts)
    _atomic_json(state_path, state)
    return state, attempts


def _simulate_specification_task(task: Mapping[str, Any]) -> dict[str, Any]:
    specification = dict(task["specification"])
    runtime = StudyRuntime(**task["runtime"])
    slim_path = Path(task["slim_path"]).resolve()
    specification_dir = Path(task["specification_dir"])
    specification_dir.mkdir(parents=True, exist_ok=True)
    demographic_model, sweep, samples, model_record = _build_simulation_objects(
        specification, runtime
    )
    contract = _validated_simulation_contract(
        specification, runtime, slim_path, model_record
    )
    prior_contract_matches = False
    prior_contract_path = specification_dir / "simulation_contract.json"
    if prior_contract_path.is_file():
        try:
            prior_contract = json.loads(prior_contract_path.read_text(encoding="utf-8"))
            prior_contract_matches = _canonical_json(prior_contract) == _canonical_json(
                contract
            )
        except (OSError, TypeError, ValueError, json.JSONDecodeError):
            prior_contract_matches = False
    cached = _valid_simulation_cache(
        specification_dir, contract, specification, runtime
    )
    if cached is not None:
        return {"status": "cached", "completion": cached}

    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    _atomic_json(specification_dir / "simulation_contract.json", contract)
    attempts_path = specification_dir / "external_draws.tsv"
    state_path = specification_dir / "external_draw_state.json"
    state, attempts = _resume_external_draw_progress(
        specification_dir,
        contract,
        specification,
        runtime,
        allow_legacy_progress=prior_contract_matches,
    )
    engine = stdpopsim.get_engine("slim")
    started_all = perf_counter()
    accepted_payload: tuple[tskit.TreeSequence, pd.DataFrame, dict[str, Any]] | None = (
        None
    )
    while len(attempts) < runtime.max_external_draws:
        draw = int(state["next_external_draw_zero_based"])
        seed = _stable_seed(runtime.base_seed, specification["specification_id"], draw)
        logfile = specification_dir / f"slim_draw_{draw:03d}.csv"
        state["active_internal_conditioned_trajectory"] = {
            "external_draw_zero_based": draw,
            "seed": seed,
            "stage": "stdpopsim_internal_conditioning",
        }
        _atomic_json(state_path, state)
        started = perf_counter()
        ts = engine.simulate(
            demographic_model,
            sweep.contig,
            samples,
            seed=seed,
            extended_events=list(sweep.extended_events),
            slim_path=str(slim_path),
            slim_scaling_factor=runtime.slim_scaling_factor,
            slim_burn_in=runtime.slim_burn_in,
            logfile=str(logfile),
            logfile_interval=100,
            keep_mutation_ids_as_alleles=False,
        )
        elapsed = perf_counter() - started
        state["active_internal_conditioned_trajectory"]["stage"] = (
            "sample_panel_validation"
        )
        _atomic_json(state_path, state)
        pair_table, validation = _validate_accepted_tree(ts, specification, runtime)
        attempts.append(
            {
                "external_draw_zero_based": draw,
                "seed": seed,
                "elapsed_seconds": elapsed,
                **validation,
            }
        )
        _atomic_frame(attempts_path, pd.DataFrame(attempts))
        state["active_internal_conditioned_trajectory"] = None
        state["next_external_draw_zero_based"] = draw + 1
        state["completed_external_sample_panel_draws"] = len(attempts)
        _atomic_json(state_path, state)
        if validation["accepted"]:
            accepted_payload = (ts, pair_table, validation)
            break
    if accepted_payload is None:
        raise RuntimeError(
            f"{specification['specification_id']} produced no acceptable sampled "
            f"genotype panel in {runtime.max_external_draws} external SLiM draws"
        )

    ts, pair_table, validation = accepted_payload
    tree_path = specification_dir / "selected.trees"
    temporary_tree = tree_path.with_name(tree_path.name + f".tmp.{os.getpid()}")
    ts.dump(temporary_tree)
    os.replace(temporary_tree, tree_path)
    pair_path = _atomic_frame(specification_dir / "sample_manifest.tsv", pair_table)
    pair_files = write_gamma_pair_files(
        pair_table, specification_dir / "pairs", require_nonempty=True
    )
    positions = np.unique(
        np.r_[
            np.arange(
                0,
                runtime.sequence_length_bp,
                runtime.output_stride_bp,
                dtype=float,
            ),
            float(runtime.focal_position_bp),
        ]
    )
    truth = tree_truth_profiles(
        ts,
        positions,
        pair_table,
        thresholds_years=DEFAULT_THRESHOLDS_YEARS,
        generation_time_years=runtime.generation_time_years,
    )
    truth_path = _atomic_frame(specification_dir / "truth_profiles.tsv.gz", truth)
    truth_contrasts = internal_genotype_contrasts(truth)
    truth_contrast_path = _atomic_frame(
        specification_dir / "truth_internal_contrasts.tsv.gz", truth_contrasts
    )
    completion = {
        "schema": SCHEMA_VERSION,
        "phase": "simulate",
        "status": "complete",
        "specification_id": specification["specification_id"],
        "contract_sha256": _contract_sha256(contract),
        "tree_sha256": _sha256(tree_path),
        "sample_manifest_sha256": _sha256(pair_path),
        "truth_profiles_sha256": _sha256(truth_path),
        "truth_internal_contrasts_sha256": _sha256(truth_contrast_path),
        "pair_files": {
            key: {
                "path": Path(path).relative_to(specification_dir).as_posix(),
                "sha256": _sha256(path),
            }
            for key, path in pair_files.items()
        },
        "accepted_external_draw_zero_based": int(
            attempts[-1]["external_draw_zero_based"]
        ),
        "attempts_to_accept": len(attempts),
        "completed_external_sample_panel_draws": len(attempts),
        "interrupted_internal_conditioned_trajectories": len(
            state["interrupted_internal_conditioned_trajectories"]
        ),
        "accepted_seed": int(attempts[-1]["seed"]),
        "validation": validation,
        "elapsed_seconds": perf_counter() - started_all,
    }
    _atomic_json(specification_dir / "simulation_complete.json", completion)
    (specification_dir / "simulation_failed.json").unlink(missing_ok=True)
    return {"status": "complete", "completion": completion}


def _filtered_records(
    records: Iterable[dict[str, Any]], specification_filters: Sequence[str]
) -> list[dict[str, Any]]:
    records = list(records)
    if not specification_filters:
        return records
    selected = [
        row
        for row in records
        if any(token in row["specification_id"] for token in specification_filters)
    ]
    if not selected:
        raise ValueError(
            "no specification IDs matched: " + ", ".join(specification_filters)
        )
    return selected


def _simulation_process_entry(task: Mapping[str, Any], result_queue) -> None:
    """Run one specification in an independently terminable process group."""
    if os.name != "nt":
        os.setsid()
    try:
        result_queue.put(
            {
                "kind": "result",
                "value": _simulate_specification_task(task),
            }
        )
    except BaseException as error:  # pragma: no cover - exercised by integration runs
        result_queue.put(
            {
                "kind": "error",
                "error_type": type(error).__name__,
                "error": str(error),
                "traceback": traceback.format_exc(),
            }
        )


def _terminate_process_tree(process: multiprocessing.Process) -> None:
    """Terminate the exact worker and its SLiM descendant after a declared timeout."""
    if not process.is_alive():
        process.join(timeout=1)
        return
    if os.name == "nt":
        subprocess.run(
            ["taskkill", "/PID", str(process.pid), "/T", "/F"],
            check=False,
            capture_output=True,
            text=True,
        )
    else:
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
        process.join(timeout=5)
        if process.is_alive():
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
    process.join(timeout=5)


def _failed_simulation_result(
    task: Mapping[str, Any],
    *,
    status: str,
    elapsed_seconds: float,
    error_type: str,
    error: str,
    traceback_text: str | None = None,
) -> dict[str, Any]:
    specification = task["specification"]
    runtime = StudyRuntime(**task["runtime"])
    specification_dir = Path(task["specification_dir"])
    contract = None
    contract_path = specification_dir / "simulation_contract.json"
    try:
        if contract_path.is_file():
            contract = json.loads(contract_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        contract = None
    progress = (
        _validated_external_draw_progress(
            specification_dir, contract, specification, runtime
        )
        if contract is not None
        else None
    )
    state = progress[0] if progress is not None else {}
    active = state.get("active_internal_conditioned_trajectory")
    failure = {
        "schema": SCHEMA_VERSION,
        "phase": "simulate",
        "status": status,
        "specification_id": specification["specification_id"],
        "elapsed_seconds": float(elapsed_seconds),
        "declared_timeout_seconds": float(
            task.get("runtime", {}).get("max_specification_seconds", 0.0)
        ),
        "error_type": error_type,
        "error": error,
        "traceback": traceback_text,
        "retryable": True,
        "simulation_contract_sha256": (
            _contract_sha256(contract) if contract is not None else None
        ),
        "completed_external_sample_panel_draws": int(
            state.get("completed_external_sample_panel_draws", 0)
        ),
        "active_internal_conditioned_trajectory": active,
        "interrupted_during_stdpopsim_internal_conditioning": bool(
            active is not None
            and active.get("stage") == "stdpopsim_internal_conditioning"
        ),
    }
    _atomic_json(specification_dir / "simulation_failed.json", failure)
    return {"status": status, "failure": failure}


def _simulate_tasks_with_timeouts(
    tasks: Sequence[Mapping[str, Any]], runtime: StudyRuntime
) -> dict[str, dict[str, Any]]:
    context = multiprocessing.get_context("spawn")
    pending = list(tasks)
    active: dict[int, dict[str, Any]] = {}
    results: dict[str, dict[str, Any]] = {}
    while pending or active:
        while pending and len(active) < runtime.simulation_workers:
            task = pending.pop(0)
            result_queue = context.Queue()
            process = context.Process(
                target=_simulation_process_entry,
                args=(task, result_queue),
            )
            process.start()
            active[int(process.pid)] = {
                "process": process,
                "queue": result_queue,
                "task": task,
                "started": perf_counter(),
            }
        for pid, state in list(active.items()):
            process = state["process"]
            task = state["task"]
            specification_id = task["specification"]["specification_id"]
            elapsed = perf_counter() - state["started"]
            message = None
            try:
                message = state["queue"].get_nowait()
            except queue.Empty:
                pass
            if message is not None:
                process.join(timeout=5)
                if message["kind"] == "result":
                    results[specification_id] = message["value"]
                else:
                    results[specification_id] = _failed_simulation_result(
                        task,
                        status="failed",
                        elapsed_seconds=elapsed,
                        error_type=message["error_type"],
                        error=message["error"],
                        traceback_text=message["traceback"],
                    )
                del active[pid]
                continue
            if (
                runtime.max_specification_seconds > 0
                and elapsed >= runtime.max_specification_seconds
            ):
                _terminate_process_tree(process)
                results[specification_id] = _failed_simulation_result(
                    task,
                    status="timed_out",
                    elapsed_seconds=elapsed,
                    error_type="SpecificationTimeout",
                    error=(
                        "SLiM conditioning exceeded the declared per-specification "
                        f"limit of {runtime.max_specification_seconds:g} seconds"
                    ),
                )
                del active[pid]
                continue
            if not process.is_alive():
                process.join(timeout=1)
                results[specification_id] = _failed_simulation_result(
                    task,
                    status="failed",
                    elapsed_seconds=elapsed,
                    error_type="WorkerExit",
                    error=f"simulation worker exited with code {process.exitcode}",
                )
                del active[pid]
        if active:
            sleep(0.25)
    return results


def simulate_studies(
    runtime: StudyRuntime,
    specifications: Mapping[str, Sequence[dict[str, Any]]],
    *,
    slim_path: str | Path,
    specification_filters: Sequence[str] = (),
) -> dict[str, dict[str, Any]]:
    slim_path = Path(slim_path).resolve()
    if not slim_path.is_file():
        raise FileNotFoundError(f"SLiM executable not found: {slim_path}")
    tasks = []
    for scenario, records in specifications.items():
        study_dir = _study_dir(runtime.repo_root, scenario)
        for record in _filtered_records(records, specification_filters):
            tasks.append(
                {
                    "specification": record,
                    "runtime": asdict(runtime),
                    "slim_path": str(slim_path),
                    "specification_dir": str(
                        study_dir / "work" / record["specification_id"]
                    ),
                }
            )
    results: dict[str, dict[str, Any]] = {}
    if runtime.max_specification_seconds > 0:
        results = _simulate_tasks_with_timeouts(tasks, runtime)
        _write_simulation_status_tables(runtime, specifications, results)
        return results
    if runtime.simulation_workers == 1:
        for task in tasks:
            specification_id = task["specification"]["specification_id"]
            started = perf_counter()
            try:
                results[specification_id] = _simulate_specification_task(task)
            except Exception as error:
                results[specification_id] = _failed_simulation_result(
                    task,
                    status="failed",
                    elapsed_seconds=perf_counter() - started,
                    error_type=type(error).__name__,
                    error=str(error),
                    traceback_text=traceback.format_exc(),
                )
        _write_simulation_status_tables(runtime, specifications, results)
        return results
    with ProcessPoolExecutor(max_workers=runtime.simulation_workers) as executor:
        futures = {
            executor.submit(_simulate_specification_task, task): task for task in tasks
        }
        for future in as_completed(futures):
            task = futures[future]
            specification_id = task["specification"]["specification_id"]
            try:
                results[specification_id] = future.result()
            except Exception as error:
                results[specification_id] = _failed_simulation_result(
                    task,
                    status="failed",
                    elapsed_seconds=0.0,
                    error_type=type(error).__name__,
                    error=str(error),
                    traceback_text=traceback.format_exc(),
                )
    _write_simulation_status_tables(runtime, specifications, results)
    return results


def _write_simulation_status_tables(
    runtime: StudyRuntime,
    specifications: Mapping[str, Sequence[dict[str, Any]]],
    results: Mapping[str, Mapping[str, Any]],
) -> None:
    for scenario, records in specifications.items():
        rows = []
        for specification in records:
            specification_id = specification["specification_id"]
            specification_dir = (
                _study_dir(runtime.repo_root, scenario) / "work" / specification_id
            )
            recorded_contract = _recorded_simulation_contract_if_current(
                specification_dir, specification, runtime
            )
            valid_completion = (
                _valid_simulation_cache(
                    specification_dir,
                    recorded_contract,
                    specification,
                    runtime,
                )
                if recorded_contract is not None
                else None
            )
            supplied_result = results.get(specification_id)
            result = supplied_result
            stale_cache = False
            if result is not None and "completion" in result:
                if valid_completion is None:
                    result = None
                    stale_cache = True
                else:
                    result = {**result, "completion": valid_completion}
            if result is None and valid_completion is not None:
                result = {"status": "complete", "completion": valid_completion}
            if result is None:
                completion_path = specification_dir / "simulation_complete.json"
                failure_path = specification_dir / "simulation_failed.json"
                stale_cache = stale_cache or completion_path.is_file()
                if failure_path.is_file():
                    try:
                        failure = json.loads(failure_path.read_text(encoding="utf-8"))
                        failure_is_current = (
                            recorded_contract is not None
                            and failure.get("simulation_contract_sha256")
                            == _contract_sha256(recorded_contract)
                        )
                        if failure_is_current:
                            result = {
                                "status": failure["status"],
                                "failure": failure,
                            }
                        else:
                            stale_cache = True
                    except (KeyError, OSError, json.JSONDecodeError):
                        stale_cache = True
            if result is None:
                status = "not_requested"
                detail = (
                    "stale or incompatible simulation cache ignored for the "
                    "current specification/runtime contract"
                    if stale_cache
                    else None
                )
            else:
                status = str(result["status"])
                detail = result.get("failure", {}).get("error")
            payload = (
                result.get("failure", result.get("completion", {}))
                if result is not None
                else {}
            )
            progress = (
                _validated_external_draw_progress(
                    specification_dir,
                    recorded_contract,
                    specification,
                    runtime,
                )
                if recorded_contract is not None
                else None
            )
            attempts = progress[1] if progress is not None else []
            draw_state = progress[0] if progress is not None else {}
            last_attempt = attempts[-1] if attempts else {}
            completed_draws = len(attempts)
            if not attempts and valid_completion is not None:
                completed_draws = int(
                    valid_completion.get(
                        "completed_external_sample_panel_draws",
                        valid_completion.get("attempts_to_accept", 0),
                    )
                )
                last_attempt = {
                    "achieved_allele_frequency": valid_completion["validation"].get(
                        "achieved_allele_frequency"
                    ),
                    "accepted": True,
                }
            rows.append(
                {
                    "specification_id": specification_id,
                    "selection_coefficient": specification["selection_coefficient"],
                    "demography_curve": specification.get("demography_curve"),
                    "target_allele_frequency": specification["target_allele_frequency"],
                    "deterministic_feasible_from_fixed_archaic_source": (
                        specification.get(
                            "deterministic_feasible_from_fixed_archaic_source"
                        )
                    ),
                    "declared_timeout_seconds": payload.get(
                        "declared_timeout_seconds", runtime.max_specification_seconds
                    ),
                    "elapsed_seconds": payload.get("elapsed_seconds"),
                    "completed_external_draws": completed_draws,
                    "interrupted_internal_conditioned_trajectories": (
                        len(
                            draw_state.get(
                                "interrupted_internal_conditioned_trajectories", []
                            )
                        )
                        + int(
                            draw_state.get("active_internal_conditioned_trajectory")
                            is not None
                            and status in {"timed_out", "failed"}
                        )
                    ),
                    "last_completed_draw_sample_af": last_attempt.get(
                        "achieved_allele_frequency"
                    ),
                    "last_completed_draw_accepted": last_attempt.get("accepted"),
                    "status": status,
                    "detail": detail,
                }
            )
        _atomic_frame(
            _study_dir(runtime.repo_root, scenario)
            / "results"
            / "simulation_status.tsv",
            _status_column_last(pd.DataFrame(rows)),
        )


def _decode_contract(
    specification: Mapping[str, Any],
    runtime: StudyRuntime,
    *,
    tree_path: Path,
    pairs_path: Path,
    decoder_path: Path,
    genotype_class: str,
) -> dict[str, Any]:
    return _decode_contract_with_sha256(
        specification,
        runtime,
        tree_path=tree_path,
        pairs_path=pairs_path,
        decoder_sha256=_sha256(decoder_path),
        genotype_class=genotype_class,
    )


def _decode_contract_with_sha256(
    specification: Mapping[str, Any],
    runtime: StudyRuntime,
    *,
    tree_path: Path,
    pairs_path: Path,
    decoder_sha256: str,
    genotype_class: str,
) -> dict[str, Any]:
    return {
        "schema": SCHEMA_VERSION,
        "phase": "decode",
        "specification_id": specification["specification_id"],
        "genotype_class": genotype_class,
        "selected_tree_sha256": _sha256(tree_path),
        "pairs_sha256": _sha256(pairs_path),
        "decoder_sha256": decoder_sha256,
        "scaled_mutation_rate": float(specification["theta_for_gamma_smc"]),
        "unscaled_mutation_rate": runtime.mutation_rate,
        "recombination_to_mutation_ratio": (
            runtime.recombination_rate / runtime.mutation_rate
        ),
        "thresholds_years": list(DEFAULT_THRESHOLDS_YEARS),
        "generation_time_years": runtime.generation_time_years,
        "output_stride_bp": runtime.output_stride_bp,
        "cache_size_bp": runtime.cache_size_bp,
        "recent_call": "median",
        "pair_scope": "explicit_within_diploid_genotype_class",
    }


def _valid_decode_cache(output: Path, contract: Mapping[str, Any]) -> bool:
    completion_path = output.with_suffix(output.suffix + ".complete.json")
    run_path = output.with_suffix(output.suffix + ".run.json")
    if not (output.is_file() and completion_path.is_file() and run_path.is_file()):
        return False
    try:
        completion = json.loads(completion_path.read_text(encoding="utf-8"))
        return (
            completion["contract_sha256"] == _contract_sha256(contract)
            and completion["summary_sha256"] == _sha256(output)
            and completion["run_sha256"] == _sha256(run_path)
            and len(pd.read_csv(output, sep="\t")) > 0
        )
    except (KeyError, OSError, ValueError, json.JSONDecodeError, pd.errors.ParserError):
        return False


def _decode_bundle_contract(
    specification: Mapping[str, Any],
    specification_dir: Path,
    simulation_completion: Mapping[str, Any],
    class_contract_sha256s: Mapping[str, str],
    class_completion_sha256s: Mapping[str, str],
) -> dict[str, Any]:
    return {
        "schema": SCHEMA_VERSION,
        "phase": "decode_bundle",
        "specification_id": specification["specification_id"],
        "simulation_contract_sha256": simulation_completion["contract_sha256"],
        "simulation_completion_sha256": _sha256(
            specification_dir / "simulation_complete.json"
        ),
        "class_decode_contract_sha256s": dict(class_contract_sha256s),
        "class_decode_completion_sha256s": dict(class_completion_sha256s),
    }


def _valid_recorded_decode_cache(
    specification_dir: Path,
    specification: Mapping[str, Any],
    runtime: StudyRuntime,
    *,
    decoder_path: str | Path | None = None,
) -> dict[str, Any] | None:
    """Validate decoded derivatives against the current simulation and Q/runtime."""
    simulation_completion = _valid_recorded_simulation_cache(
        specification_dir, specification, runtime
    )
    if simulation_completion is None:
        return None
    tree_path = specification_dir / "selected.trees"
    current_decoder_sha256 = None
    if decoder_path is not None:
        current_decoder = Path(decoder_path).resolve()
        if not current_decoder.is_file():
            return None
        current_decoder_sha256 = _sha256(current_decoder)
    class_contract_sha256s: dict[str, str] = {}
    class_completion_sha256s: dict[str, str] = {}
    try:
        for genotype_class in ("overall", "hom_ref", "heterozygous", "hom_alt"):
            pairs_path = specification_dir / "pairs" / f"{genotype_class}.pairs.tsv"
            output = specification_dir / "decoded" / f"{genotype_class}.summary.tsv"
            contract_path = output.with_suffix(output.suffix + ".contract.json")
            completion_path = output.with_suffix(output.suffix + ".complete.json")
            if not contract_path.is_file():
                return None
            recorded_contract = json.loads(contract_path.read_text(encoding="utf-8"))
            decoder_sha256 = (
                current_decoder_sha256
                if current_decoder_sha256 is not None
                else str(recorded_contract["decoder_sha256"])
            )
            expected_contract = _decode_contract_with_sha256(
                specification,
                runtime,
                tree_path=tree_path,
                pairs_path=pairs_path,
                decoder_sha256=decoder_sha256,
                genotype_class=genotype_class,
            )
            if _canonical_json(recorded_contract) != _canonical_json(
                expected_contract
            ) or not _valid_decode_cache(output, recorded_contract):
                return None
            class_contract_sha256s[genotype_class] = _contract_sha256(recorded_contract)
            class_completion_sha256s[genotype_class] = _sha256(completion_path)
        expected_bundle = _decode_bundle_contract(
            specification,
            specification_dir,
            simulation_completion,
            class_contract_sha256s,
            class_completion_sha256s,
        )
        bundle_path = specification_dir / "decode_contract.json"
        completion_path = specification_dir / "decode_complete.json"
        decoded_path = specification_dir / "decoded_profiles.tsv.gz"
        contrast_path = specification_dir / "decoded_internal_contrasts.tsv.gz"
        if not all(
            path.is_file()
            for path in (bundle_path, completion_path, decoded_path, contrast_path)
        ):
            return None
        recorded_bundle = json.loads(bundle_path.read_text(encoding="utf-8"))
        completion = json.loads(completion_path.read_text(encoding="utf-8"))
        if (
            _canonical_json(recorded_bundle) != _canonical_json(expected_bundle)
            or completion.get("schema") != SCHEMA_VERSION
            or completion.get("phase") != "decode"
            or completion.get("status") != "complete"
            or completion.get("specification_id") != specification["specification_id"]
            or completion.get("contract_sha256") != _contract_sha256(recorded_bundle)
            or completion.get("decoded_profiles_sha256") != _sha256(decoded_path)
            or completion.get("decoded_internal_contrasts_sha256")
            != _sha256(contrast_path)
            or pd.read_csv(decoded_path, sep="\t").empty
            or pd.read_csv(contrast_path, sep="\t").empty
        ):
            return None
        return completion
    except (
        KeyError,
        OSError,
        TypeError,
        ValueError,
        json.JSONDecodeError,
        pd.errors.EmptyDataError,
        pd.errors.ParserError,
    ):
        return None


def decode_specification(
    specification: Mapping[str, Any],
    runtime: StudyRuntime,
    *,
    decoder_path: str | Path,
) -> dict[str, Any]:
    scenario = str(specification["study_type"])
    study_dir = _study_dir(runtime.repo_root, scenario)
    specification_dir = study_dir / "work" / specification["specification_id"]
    tree_path = specification_dir / "selected.trees"
    simulation_completion = _valid_recorded_simulation_cache(
        specification_dir, specification, runtime
    )
    if simulation_completion is None:
        raise ValueError(
            "simulation is incomplete, corrupt, or incompatible with the current "
            f"specification/runtime contract: {specification['specification_id']}"
        )
    _validate_accepted_tree(tskit.load(tree_path), specification, runtime)
    decoder_path = Path(decoder_path).resolve()
    if not decoder_path.is_file():
        raise FileNotFoundError(f"Gamma-SMC executable not found: {decoder_path}")
    summaries: dict[str, pd.DataFrame] = {}
    class_contract_sha256s: dict[str, str] = {}
    class_completion_sha256s: dict[str, str] = {}
    decode_dir = specification_dir / "decoded"
    decode_dir.mkdir(parents=True, exist_ok=True)
    for genotype_class in ("overall", "hom_ref", "heterozygous", "hom_alt"):
        pairs_path = specification_dir / "pairs" / f"{genotype_class}.pairs.tsv"
        output = decode_dir / f"{genotype_class}.summary.tsv"
        contract = _decode_contract(
            specification,
            runtime,
            tree_path=tree_path,
            pairs_path=pairs_path,
            decoder_path=decoder_path,
            genotype_class=genotype_class,
        )
        contract_path = output.with_suffix(output.suffix + ".contract.json")
        _atomic_json(contract_path, contract)
        if not _valid_decode_cache(output, contract):
            run_within_decoder(
                decoder_path,
                tree_path,
                output,
                scaled_mutation_rate=float(specification["theta_for_gamma_smc"]),
                recombination_to_mutation_ratio=(
                    runtime.recombination_rate / runtime.mutation_rate
                ),
                mutation_rate=runtime.mutation_rate,
                threshold_years=DEFAULT_THRESHOLDS_YEARS,
                generation_time=runtime.generation_time_years,
                input_format="trees",
                output_at_stride=runtime.output_stride_bp,
                output_at_hets=False,
                only_within=False,
                pairs_file=pairs_path,
                recent_call="median",
                threads=runtime.threads,
                cache_size=runtime.cache_size_bp,
            )
            run_path = output.with_suffix(output.suffix + ".run.json")
            _atomic_json(
                output.with_suffix(output.suffix + ".complete.json"),
                {
                    "contract_sha256": _contract_sha256(contract),
                    "summary_sha256": _sha256(output),
                    "run_sha256": _sha256(run_path),
                },
            )
        summaries[genotype_class] = pd.read_csv(output, sep="\t")
        class_contract_sha256s[genotype_class] = _contract_sha256(contract)
        class_completion_sha256s[genotype_class] = _sha256(
            output.with_suffix(output.suffix + ".complete.json")
        )
    decoded = gamma_summaries_to_long(
        summaries,
        generation_time_years=runtime.generation_time_years,
        source="gamma_smc",
    )
    decoded_path = _atomic_frame(specification_dir / "decoded_profiles.tsv.gz", decoded)
    decoded_contrasts = internal_genotype_contrasts(decoded)
    contrast_path = _atomic_frame(
        specification_dir / "decoded_internal_contrasts.tsv.gz", decoded_contrasts
    )
    bundle_contract = _decode_bundle_contract(
        specification,
        specification_dir,
        simulation_completion,
        class_contract_sha256s,
        class_completion_sha256s,
    )
    _atomic_json(specification_dir / "decode_contract.json", bundle_contract)
    completion = {
        "schema": SCHEMA_VERSION,
        "phase": "decode",
        "status": "complete",
        "specification_id": specification["specification_id"],
        "contract_sha256": _contract_sha256(bundle_contract),
        "decoded_profiles_sha256": _sha256(decoded_path),
        "decoded_internal_contrasts_sha256": _sha256(contrast_path),
    }
    _atomic_json(specification_dir / "decode_complete.json", completion)
    return completion


def _decode_specification_task(
    specification: Mapping[str, Any],
    runtime_record: Mapping[str, Any],
    decoder_path: str,
) -> dict[str, Any]:
    return decode_specification(
        specification,
        StudyRuntime(**runtime_record),
        decoder_path=decoder_path,
    )


def decode_studies(
    runtime: StudyRuntime,
    specifications: Mapping[str, Sequence[dict[str, Any]]],
    *,
    decoder_path: str | Path,
    specification_filters: Sequence[str] = (),
) -> dict[str, dict[str, Any]]:
    if os.name == "nt":
        raise RuntimeError(
            "the bundled Gamma-SMC executable is Linux/AVX2; run the decode phase "
            "inside WSL2 or Linux"
        )
    decoder_resolved = Path(decoder_path).resolve()
    if not decoder_resolved.is_file():
        raise FileNotFoundError(f"Gamma-SMC executable not found: {decoder_resolved}")
    decoder = str(decoder_resolved)
    results: dict[str, dict[str, Any]] = {}
    for scenario, records in specifications.items():
        selected = _filtered_records(records, specification_filters)
        selected_ids = {record["specification_id"] for record in selected}
        study_dir = _study_dir(runtime.repo_root, scenario)
        ready: list[dict[str, Any]] = []
        selected_status: dict[str, dict[str, Any]] = {}
        for record in selected:
            specification_id = record["specification_id"]
            specification_dir = study_dir / "work" / specification_id
            if (
                _valid_recorded_simulation_cache(specification_dir, record, runtime)
                is not None
            ):
                ready.append(record)
            else:
                selected_status[specification_id] = {
                    "status": "not_simulated",
                    "error_type": "ValueError",
                    "error": (
                        "simulation is incomplete, corrupt, or incompatible with "
                        "the current specification/runtime contract"
                    ),
                }
        if ready:
            with ProcessPoolExecutor(
                max_workers=min(runtime.simulation_workers, len(ready))
            ) as executor:
                futures = {
                    executor.submit(
                        _decode_specification_task,
                        record,
                        asdict(runtime),
                        decoder,
                    ): record
                    for record in ready
                }
                for future in as_completed(futures):
                    record = futures[future]
                    specification_id = record["specification_id"]
                    try:
                        completion = future.result()
                    except Exception as error:  # pragma: no cover - integration path
                        selected_status[specification_id] = {
                            "status": "failed",
                            "error_type": type(error).__name__,
                            "error": str(error),
                        }
                    else:
                        results[specification_id] = completion
                        selected_status[specification_id] = {
                            "status": "complete",
                            "error_type": None,
                            "error": None,
                        }
        status_rows: list[dict[str, Any]] = []
        for record in records:
            specification_id = record["specification_id"]
            specification_dir = study_dir / "work" / specification_id
            status = selected_status.get(specification_id)
            if status is None and (
                _valid_recorded_decode_cache(
                    specification_dir,
                    record,
                    runtime,
                    decoder_path=decoder_resolved,
                )
                is not None
            ):
                status = {
                    "status": "complete",
                    "error_type": None,
                    "error": None,
                }
            elif status is None and specification_id not in selected_ids:
                simulation_ready = (
                    _valid_recorded_simulation_cache(specification_dir, record, runtime)
                    is not None
                )
                status = {
                    "status": "not_requested" if simulation_ready else "not_simulated",
                    "error_type": None,
                    "error": None,
                }
            if status is None:
                status = {
                    "status": "failed",
                    "error_type": "DecodeStateError",
                    "error": "decode task produced no current-contract completion",
                }
            status_rows.append({"specification_id": specification_id, **status})
        _atomic_frame(
            study_dir / "results" / "decode_status.tsv",
            _status_column_last(
                pd.DataFrame(status_rows).sort_values("specification_id")
            ),
        )
    return results


def _center_rows(frame: pd.DataFrame, focal_position: float) -> pd.DataFrame:
    mask = np.isclose(
        frame["position_0based"].to_numpy(dtype=float),
        float(focal_position),
        rtol=0.0,
        atol=1e-9,
    )
    center = frame.loc[mask].copy()
    if center.empty:
        raise ValueError(f"profile contains no focal row at {focal_position:g}")
    return center


def _validated_plot_inputs(
    runtime: StudyRuntime,
    scenario: str,
    specifications: Sequence[Mapping[str, Any]],
    *,
    decoder_path: str | Path | None,
) -> tuple[dict[str, dict[str, Any]], dict[str, dict[str, Any]], Path]:
    """Validate plot inputs and rewrite decode status from the same decision.

    With no explicit decoder, recorded decoder contracts remain portable and
    are validated by their recorded executable and output hashes. Supplying a
    decoder makes its SHA-256 an additional requirement. The returned cache
    maps are reused by aggregation so status cannot claim ``complete`` for a
    decode that plotting then silently excludes.
    """

    expected_decoder: Path | None = None
    expected_decoder_sha256: str | None = None
    if decoder_path is not None:
        expected_decoder = Path(decoder_path).resolve()
        if not expected_decoder.is_file():
            raise FileNotFoundError(
                f"Gamma-SMC executable not found for plot validation: {expected_decoder}"
            )
        expected_decoder_sha256 = _sha256(expected_decoder)

    study_dir = _study_dir(runtime.repo_root, scenario)
    simulations: dict[str, dict[str, Any]] = {}
    decodes: dict[str, dict[str, Any]] = {}
    status_rows: list[dict[str, Any]] = []
    for specification in specifications:
        specification_id = str(specification["specification_id"])
        specification_dir = study_dir / "work" / specification_id
        simulation = _valid_recorded_simulation_cache(
            specification_dir, specification, runtime
        )
        if simulation is not None:
            simulations[specification_id] = simulation
        decode = (
            _valid_recorded_decode_cache(
                specification_dir,
                specification,
                runtime,
                decoder_path=expected_decoder,
            )
            if simulation is not None
            else None
        )
        if decode is not None:
            decodes[specification_id] = decode
        status_rows.append(
            {
                "specification_id": specification_id,
                "status": (
                    "complete"
                    if decode is not None
                    else "not_requested"
                    if simulation is not None
                    else "not_simulated"
                ),
                "error_type": None,
                "error": None,
                "decoder_validation_mode": (
                    "explicit_decoder_sha256"
                    if expected_decoder is not None
                    else "recorded_decode_contract"
                ),
                "expected_decoder_sha256": expected_decoder_sha256,
            }
        )
    status_path = _atomic_frame(
        study_dir / "results" / "decode_status.tsv",
        _status_column_last(pd.DataFrame(status_rows).sort_values("specification_id")),
    )
    return simulations, decodes, status_path


def plot_and_aggregate_studies(
    runtime: StudyRuntime,
    specifications: Mapping[str, Sequence[dict[str, Any]]],
    *,
    specification_filters: Sequence[str] = (),
    decoder_path: str | Path | None = None,
) -> dict[str, dict[str, str]]:
    _write_simulation_status_tables(runtime, specifications, {})
    outputs: dict[str, dict[str, str]] = {}
    for scenario, all_records in specifications.items():
        plot_records = _filtered_records(all_records, specification_filters)
        plot_specification_ids = {
            specification["specification_id"] for specification in plot_records
        }
        study_dir = _study_dir(runtime.repo_root, scenario)
        results_dir = study_dir / "results"
        valid_simulations, valid_decodes, decode_status_path = _validated_plot_inputs(
            runtime,
            scenario,
            all_records,
            decoder_path=decoder_path,
        )
        status_path = results_dir / "simulation_status.tsv"
        simulation_status = pd.read_csv(status_path, sep="\t")
        simulation_status_counts = {
            str(status): int(count)
            for status, count in simulation_status["status"].value_counts().items()
        }
        status_plot_outputs = plot_simulation_status_grid(
            simulation_status,
            results_dir
            / "figures"
            / "cross_specification"
            / f"{scenario}_simulation_status",
        )
        design_and_contract_artifacts: dict[str, Path] = {
            "specifications_tsv": study_dir / "specifications.tsv",
            "study_design_json": study_dir / "study_design.json",
        }
        design_and_contract_artifacts["decode_status_tsv"] = decode_status_path
        decode_status = pd.read_csv(decode_status_path, sep="\t")
        decode_status_counts = {
            str(status): int(count)
            for status, count in decode_status["status"].value_counts().items()
        }
        if scenario == "no_introgression":
            design_and_contract_artifacts.update(
                {
                    "origin_age_contracts_json": study_dir
                    / "origin_age_contracts.json",
                    "eas_resource_npz": study_dir / "resources" / "EAS.npz",
                    "eas_resource_provenance_json": study_dir
                    / "resources"
                    / "EAS.provenance.json",
                }
            )
        for design_path in sorted((results_dir / "figures" / "design").glob("*")):
            if design_path.is_file() and design_path.suffix in {".png", ".pdf"}:
                design_and_contract_artifacts[
                    f"design_{design_path.stem}_{design_path.suffix[1:]}"
                ] = design_path
        acceptance_rows: list[dict[str, Any]] = []
        center_contrast_rows: list[pd.DataFrame] = []
        center_profile_rows: list[pd.DataFrame] = []
        representative_paths: dict[str, Path] = {}
        simulation_software: dict[str, dict[str, Any]] = {}
        decoder_sha256s: set[str] = set()
        for specification in all_records:
            specification_id = specification["specification_id"]
            specification_dir = study_dir / "work" / specification_id
            simulation = valid_simulations.get(specification_id)
            if simulation is None:
                continue
            simulation_contract = _recorded_simulation_contract_if_current(
                specification_dir, specification, runtime
            )
            if simulation_contract is not None:
                software = dict(simulation_contract.get("software", {}))
                software.pop("slim_path", None)
                simulation_software[_canonical_json(software)] = software
            truth = pd.read_csv(specification_dir / "truth_profiles.tsv.gz", sep="\t")
            decoded_path = specification_dir / "decoded_profiles.tsv.gz"
            profile = truth
            valid_decode = valid_decodes.get(specification_id)
            if valid_decode is not None:
                decoded = pd.read_csv(decoded_path, sep="\t")
                profile = pd.concat([truth, decoded], ignore_index=True)
                for decoder_contract_path in (specification_dir / "decoded").glob(
                    "*.contract.json"
                ):
                    decoder_contract = json.loads(
                        decoder_contract_path.read_text(encoding="utf-8")
                    )
                    decoder_sha256s.add(decoder_contract["decoder_sha256"])
            if specification_id in plot_specification_ids:
                figure_dir = results_dir / "figures" / specification_id
                plot_outputs = write_publication_plots(
                    profile,
                    runtime.focal_position_bp,
                    figure_dir,
                    prefix=specification_id,
                    generation_time_years=runtime.generation_time_years,
                )
                if specification_id == REPRESENTATIVE_SPECIFICATIONS[scenario]:
                    representative_dir = results_dir / "figures" / "representative"
                    for plot_name, formats in plot_outputs.items():
                        for suffix, source_path in formats.items():
                            destination = representative_dir / Path(source_path).name
                            _atomic_copy(source_path, destination)
                            representative_paths[f"{plot_name}_{suffix}"] = destination
            acceptance_rows.append(
                {
                    "scenario": specification["scenario"],
                    "study_type": specification["study_type"],
                    "specification_id": specification_id,
                    "selection_coefficient": specification["selection_coefficient"],
                    "demography_curve": specification.get("demography_curve"),
                    "target_allele_frequency": specification["target_allele_frequency"],
                    "achieved_allele_frequency": simulation["validation"][
                        "achieved_allele_frequency"
                    ],
                    "achieved_allele_frequency_scope": "sample_estimate",
                    "attempts_to_accept": simulation["attempts_to_accept"],
                    "completed_external_sample_panel_draws": simulation.get(
                        "completed_external_sample_panel_draws",
                        simulation["attempts_to_accept"],
                    ),
                    "interrupted_internal_conditioned_trajectories": simulation.get(
                        "interrupted_internal_conditioned_trajectories", 0
                    ),
                    "accepted_seed": simulation["accepted_seed"],
                    "n_hom_ref": simulation["validation"]["n_hom_ref"],
                    "n_heterozygous": simulation["validation"]["n_heterozygous"],
                    "n_hom_alt": simulation["validation"]["n_hom_alt"],
                    "slim_scaling_factor": runtime.slim_scaling_factor,
                }
            )
            contrasts = internal_genotype_contrasts(profile)
            contrasts = _center_rows(contrasts, runtime.focal_position_bp)
            for key in (
                "scenario",
                "study_type",
                "specification_id",
                "selection_coefficient",
                "demography_curve",
                "target_allele_frequency",
            ):
                contrasts[key] = specification.get(key)
            center_contrast_rows.append(contrasts)
            center = _center_rows(profile, runtime.focal_position_bp)
            for key in (
                "scenario",
                "study_type",
                "specification_id",
                "selection_coefficient",
                "demography_curve",
                "target_allele_frequency",
                "origin_age_generations",
            ):
                center[key] = specification.get(key)
            center_profile_rows.append(center)
        status_artifacts = {
            "simulation_status_tsv": status_path,
            **design_and_contract_artifacts,
            **{
                f"simulation_status_{suffix}": path
                for suffix, path in status_plot_outputs.items()
            },
        }
        if not acceptance_rows:
            _atomic_json(
                results_dir / "results_manifest.json",
                {
                    "schema": SCHEMA_VERSION,
                    "phase": "plot",
                    "specifications_plotted": 0,
                    "simulation_status_counts": simulation_status_counts,
                    "decode_status_counts": decode_status_counts,
                    "artifacts": {
                        key: {
                            "path": value.relative_to(study_dir).as_posix(),
                            "sha256": _sha256(value),
                        }
                        for key, value in sorted(status_artifacts.items())
                    },
                    "interpretation": (
                        "no valid current-contract simulations available; status only"
                    ),
                    "runtime": _portable_runtime_record(runtime),
                    "simulation_software": [],
                    "gamma_smc_decoder_sha256s": [],
                },
            )
            continue
        acceptance = pd.DataFrame(acceptance_rows)
        center_contrasts = pd.concat(center_contrast_rows, ignore_index=True)
        center_profiles = pd.concat(center_profile_rows, ignore_index=True)
        acceptance_path = _atomic_frame(
            results_dir / "acceptance_summary.tsv", acceptance
        )
        contrast_path = _atomic_frame(
            results_dir / "center_internal_contrasts.tsv", center_contrasts
        )
        profile_path = _atomic_frame(
            results_dir / "center_tmrca_profiles.tsv", center_profiles
        )
        cross_plot_outputs = write_cross_spec_study_plots(
            center_contrasts,
            acceptance,
            results_dir / "figures" / "cross_specification",
            prefix=scenario,
        )
        result = {
            "acceptance_summary": str(acceptance_path),
            "center_internal_contrasts": str(contrast_path),
            "center_tmrca_profiles": str(profile_path),
        }
        artifact_paths = {
            **{key: Path(value) for key, value in result.items()},
            **status_artifacts,
            **{
                f"representative_{key}": value
                for key, value in representative_paths.items()
            },
        }
        for scenario_outputs in cross_plot_outputs.values():
            for plot_name, formats in scenario_outputs.items():
                for suffix, path in formats.items():
                    artifact_paths[f"cross_{plot_name}_{suffix}"] = path
        _atomic_json(
            results_dir / "results_manifest.json",
            {
                "schema": SCHEMA_VERSION,
                "phase": "plot",
                "specifications_plotted": len(acceptance),
                "simulation_status_counts": simulation_status_counts,
                "decode_status_counts": decode_status_counts,
                "artifacts": {
                    key: {
                        "path": Path(value).relative_to(study_dir).as_posix(),
                        "sha256": _sha256(value),
                    }
                    for key, value in sorted(artifact_paths.items())
                },
                "interpretation": (
                    "descriptive internal genotype contrasts; no neutral-null p-values"
                ),
                "runtime": _portable_runtime_record(runtime),
                "simulation_software": list(simulation_software.values()),
                "gamma_smc_decoder_sha256s": sorted(decoder_sha256s),
            },
        )
        outputs[scenario] = result
    return outputs


def _resolve_scenarios(value: str) -> tuple[str, ...]:
    if value == "both":
        return ("no_introgression", "introgression")
    return (value,)


def _default_slim(repo_root: Path) -> Path:
    if os.name == "nt":
        return repo_root / ".native-stdpopsim" / "Library" / "bin" / "slim.exe"
    return repo_root / ".native-stdpopsim" / "bin" / "slim"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run phased EAS no-introgression and introgression sweep studies"
    )
    parser.add_argument("phase", choices=("plan", "simulate", "decode", "plot", "all"))
    parser.add_argument(
        "--scenario",
        choices=("no_introgression", "introgression", "both"),
        default="both",
    )
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--slim-bin", type=Path)
    parser.add_argument("--decoder-bin", type=Path)
    parser.add_argument("--sample-diploids", type=int, default=DEFAULT_SAMPLE_DIPLOIDS)
    parser.add_argument("--workers", type=int, default=DEFAULT_SIMULATION_WORKERS)
    parser.add_argument("--threads", type=int, default=DEFAULT_THREADS)
    parser.add_argument(
        "--slim-scaling-factor", type=float, default=DEFAULT_SLIM_SCALING_FACTOR
    )
    parser.add_argument("--slim-burn-in", type=float, default=DEFAULT_SLIM_BURN_IN)
    parser.add_argument(
        "--max-external-draws", type=int, default=DEFAULT_MAX_EXTERNAL_DRAWS
    )
    parser.add_argument(
        "--spec-timeout-minutes",
        type=float,
        default=0.0,
        help="terminate and provenance-mark a conditioning cell after this wall time; 0 is unlimited",
    )
    parser.add_argument(
        "--spec",
        action="append",
        default=[],
        help="run only specification IDs containing this text (repeatable)",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    repo_root = args.repo_root.resolve()
    runtime = StudyRuntime(
        repo_root=str(repo_root),
        sample_diploids=args.sample_diploids,
        threads=args.threads,
        simulation_workers=args.workers,
        slim_scaling_factor=args.slim_scaling_factor,
        slim_burn_in=args.slim_burn_in,
        max_external_draws=args.max_external_draws,
        max_specification_seconds=args.spec_timeout_minutes * 60.0,
    )
    scenarios = _resolve_scenarios(args.scenario)
    specifications = write_study_plans(runtime, scenarios)
    slim_path = args.slim_bin or _default_slim(repo_root)
    decoder_path = args.decoder_bin or (repo_root / "bin" / "gamma_smc")
    if args.phase in {"simulate", "all"}:
        simulate_studies(
            runtime,
            specifications,
            slim_path=slim_path,
            specification_filters=args.spec,
        )
    if args.phase in {"decode", "all"}:
        decode_studies(
            runtime,
            specifications,
            decoder_path=decoder_path,
            specification_filters=args.spec,
        )
    if args.phase in {"plot", "all"}:
        plot_decoder_path = decoder_path if args.phase == "all" else args.decoder_bin
        plot_and_aggregate_studies(
            runtime,
            specifications,
            specification_filters=args.spec,
            decoder_path=plot_decoder_path,
        )
    return 0


__all__ = [
    "BASE_SEED",
    "EAS_RESOURCE_SHA256",
    "StudyRuntime",
    "build_parser",
    "build_specifications",
    "decode_specification",
    "decode_studies",
    "main",
    "plot_and_aggregate_studies",
    "simulate_studies",
    "write_study_plans",
]
