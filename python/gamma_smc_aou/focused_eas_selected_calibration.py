"""Frozen EAS selected-trajectory acceptance calibration.

This module is intentionally separate from the focused campaign executor.  It
calibrates the acceptance rate of the *existing* Q=5 EAS de-novo sweep cells;
it does not tune selection, mutation age, Gamma-SMC, or any demographic
parameter.  The selected mutation retains continuous selection from its
prespecified origin to the present and restores the original lower and upper
terminal population-frequency conditions inside SLiM.  A completed trajectory
must still pass the external 500-diploid AF gate and admit an exact-k,
genotype-QC 100-diploid panel.

Three seed-disjoint phases are mandatory before a production authorization can
be frozen:

``screen20``
    Twenty 2-bp forward-only draws per s-by-AF cell.
``confirm100``
    One hundred fresh 2-bp forward-only draws per cell.
``sensitivity10mb20``
    Twenty fresh 10-Mb draws per cell using normal stdpopsim recapitation and
    the production exact-panel construction.

The 2-bp phases are speed proxies only and can never be used as production tree
sequences.  The frozen artifact is an authorization record that a later,
explicit production integration can consume; importing this module does not
change the campaign simulation contract.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
from collections.abc import Mapping, Sequence
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, replace
from functools import lru_cache
from pathlib import Path
from time import perf_counter
from typing import Any

import numpy as np
import pandas as pd
import pyslim
import stdpopsim
import tskit

from . import focused_selection_simulation as focused_simulation
from .eas_sweep_analysis import build_diploid_pair_table
from .eas_sweep_models import (
    DOMINANCE_COEFFICIENT,
    EAS_POPULATION,
    FOCAL_POSITION_BP,
    FOCAL_SITE_ID,
    SEQUENCE_LENGTH_BP,
    build_de_novo_origin_age_grid,
    build_eas_demography_models,
    build_focal_contig,
    build_no_introgression_sweep_spec,
    load_phlash_eas_npz,
    serialize_extended_events,
)
from .eas_sweep_study import BASE_SEED, EAS_RESOURCE_SHA256
from .focused_selection_simulation import (
    DEFAULT_SELECTED_POOL_DIPLOIDS,
    DEFAULT_SLIM_BURN_IN,
    DEFAULT_SLIM_SCALING_FACTOR,
    SimulationDrawTimeout,
    _canonical_sha256,
    _selected_focal_expectation,
    _stable_seed,
    _wall_clock_timeout,
    candidate_pool_frequency_gate,
    exclusive_unit_lock,
    feasible_exact_panel_genotype_triples,
    inspect_selected_focal_identity,
    select_exact_sample_panel,
    sha256_file,
    validate_selected_focal_identity,
)

SCHEMA_VERSION = "gamma-smc.eas-selected-calibration/v2"
PLAN_MANIFEST_SCHEMA = "gamma-smc.eas-selected-calibration-plan/v2"
FROZEN_SCHEMA_VERSION = "gamma-smc.eas-selected-production-authorization/v2"
DRAW_CONTRACT_SCHEMA = "gamma-smc.eas-selected-calibration-draw/v2"

PHASE_SCREEN = "screen20"
PHASE_CONFIRM = "confirm100"
PHASE_SENSITIVITY = "sensitivity10mb20"
PHASE_ORDER = (PHASE_SCREEN, PHASE_CONFIRM, PHASE_SENSITIVITY)
PHASE_DRAWS = {
    PHASE_SCREEN: 20,
    PHASE_CONFIRM: 100,
    PHASE_SENSITIVITY: 20,
}
PHASE_MODES = {
    PHASE_SCREEN: "focal_2bp_forward_only_no_recap_no_neutral_overlay",
    PHASE_CONFIRM: "focal_2bp_forward_only_no_recap_no_neutral_overlay",
    PHASE_SENSITIVITY: "full_10mb_normal_stdpopsim_processing_exact_panel",
}

DEFAULT_SELECTION_COEFFICIENTS = (0.01, 0.005)
DEFAULT_TARGET_FREQUENCIES = (0.10, 0.20, 0.30)
DEFAULT_AF_HALF_WIDTH = 0.025
DEFAULT_SAMPLE_DIPLOIDS = 100
DEFAULT_POOL_DIPLOIDS = DEFAULT_SELECTED_POOL_DIPLOIDS
DEFAULT_MINIMUM_HIT_RATE = 0.20
DEFAULT_DRAW_TIMEOUT_SECONDS = 5 * 60.0
FOCAL_ONLY_LENGTH_BP = 2
FOCAL_ONLY_POSITION_BP = 1
EAS_RESOURCE_PATH = "no_introgression_EAS_sim/resources/EAS.npz"
MODULE_SOURCE_PATH = "python/gamma_smc_aou/focused_eas_selected_calibration.py"
IMPLEMENTATION_SOURCE_PATHS = (
    MODULE_SOURCE_PATH,
    "python/gamma_smc_aou/focused_selection_simulation.py",
    "python/gamma_smc_aou/eas_sweep_models.py",
    "python/gamma_smc_aou/eas_sweep_analysis.py",
    "python/gamma_smc_aou/eas_sweep_study.py",
)
HAN_SEED_PLAN_FILENAMES = (
    "han_selection_end_plan.tsv",
    "han_selection_end_reconciliation_plan.tsv",
    "han_selection_end_extension_plan.tsv",
    "han_frequency_conditioned_all_phase_plan.tsv",
)
HAN_FALLBACK_PLAN_FILENAME = "han_frequency_conditioned_all_phase_plan.tsv"
HAN_FALLBACK_MANIFEST_RELATIVE_PATH = (
    "calibration/han_frequency_conditioned/han_frequency_conditioned_manifest.json"
)
HAN_FALLBACK_MANIFEST_SCHEMA = "gamma-smc.han-frequency-conditioned-manifest/v1"
PRODUCTION_INTEGRATION_MUTABLE_SOURCE_PATHS = (
    "python/gamma_smc_aou/focused_selection_simulation.py",
)

PLAN_COLUMNS = (
    "calibration_id",
    "phase",
    "processing_mode",
    "selection_coefficient",
    "target_allele_frequency",
    "population_af_lower",
    "population_af_upper",
    "candidate_pool_diploids",
    "sample_diploids",
    "exact_sample_alt_count",
    "draw_index",
    "seed",
    "panel_seed",
    "slim_scaling_factor",
    "slim_burn_in",
    "sequence_length_bp",
    "focal_position_bp",
    "origin_age_generations",
    "origin_contract_sha256",
    "event_contract_sha256",
    "cell_contract_sha256",
    "slim_sha256",
)

RESULT_COLUMNS = PLAN_COLUMNS + (
    "status",
    "evaluable",
    "terminal_class",
    "candidate_pool_alt_count",
    "candidate_pool_af",
    "candidate_pool_frequency_gate_passed",
    "exact_panel_feasible",
    "pool_band_and_panel_hit",
    "sample_alt_count",
    "sample_af",
    "n_hom_ref",
    "n_heterozygous",
    "n_hom_alt",
    "stdpopsim_recapitation_applied",
    "post_slim_neutral_mutation_overlay_applied",
    "elapsed_seconds",
    "trajectory_sha256",
    "panel_tree_sha256",
    "panel_manifest_sha256",
    "error",
)


@dataclass(frozen=True)
class EasSelectedCalibrationDesign:
    """Immutable scientific and execution grid for EAS acceptance calibration."""

    selection_coefficients: tuple[float, ...] = DEFAULT_SELECTION_COEFFICIENTS
    target_frequencies: tuple[float, ...] = DEFAULT_TARGET_FREQUENCIES
    slim_scaling_factor: float = DEFAULT_SLIM_SCALING_FACTOR
    slim_burn_in: float = DEFAULT_SLIM_BURN_IN
    candidate_pool_diploids: int = DEFAULT_POOL_DIPLOIDS
    sample_diploids: int = DEFAULT_SAMPLE_DIPLOIDS
    af_half_width: float = DEFAULT_AF_HALF_WIDTH
    minimum_hit_rate: float = DEFAULT_MINIMUM_HIT_RATE
    base_seed: int = BASE_SEED

    def validate(self) -> None:
        if tuple(self.selection_coefficients) != DEFAULT_SELECTION_COEFFICIENTS:
            raise ValueError("EAS calibration requires the frozen s=0.01,0.005 grid")
        if tuple(self.target_frequencies) != DEFAULT_TARGET_FREQUENCIES:
            raise ValueError(
                "EAS calibration requires the frozen AF=0.10,0.20,0.30 grid"
            )
        if not math.isclose(
            float(self.slim_scaling_factor),
            DEFAULT_SLIM_SCALING_FACTOR,
            rel_tol=0.0,
            abs_tol=1e-12,
        ):
            raise ValueError("EAS calibration requires SLiM scaling Q=5")
        if not math.isclose(
            float(self.slim_burn_in),
            DEFAULT_SLIM_BURN_IN,
            rel_tol=0.0,
            abs_tol=1e-12,
        ):
            raise ValueError("EAS calibration requires the frozen SLiM burn-in")
        if int(self.candidate_pool_diploids) != DEFAULT_POOL_DIPLOIDS:
            raise ValueError("EAS calibration requires 500 candidate diploids")
        if int(self.sample_diploids) != DEFAULT_SAMPLE_DIPLOIDS:
            raise ValueError("EAS calibration requires 100 analysis diploids")
        if not math.isclose(
            float(self.af_half_width),
            DEFAULT_AF_HALF_WIDTH,
            rel_tol=0.0,
            abs_tol=1e-12,
        ):
            raise ValueError("EAS calibration requires the target +/- 0.025 pool band")
        if not 0.0 < float(self.minimum_hit_rate) <= 1.0:
            raise ValueError("minimum hit rate must lie in (0, 1]")
        if int(self.base_seed) <= 0:
            raise ValueError("base seed must be positive")

    def to_record(self) -> dict[str, Any]:
        self.validate()
        return {
            "selection_coefficients": [float(x) for x in self.selection_coefficients],
            "target_allele_frequencies": [float(x) for x in self.target_frequencies],
            "slim_scaling_factor": float(self.slim_scaling_factor),
            "slim_burn_in": float(self.slim_burn_in),
            "candidate_pool_diploids": int(self.candidate_pool_diploids),
            "sample_diploids": int(self.sample_diploids),
            "af_half_width": float(self.af_half_width),
            "minimum_hit_rate": float(self.minimum_hit_rate),
            "base_seed": int(self.base_seed),
            "phase_draws_per_cell": dict(PHASE_DRAWS),
        }


def _atomic_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f"{path.name}.tmp.{os.getpid()}")
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
    temporary = path.with_name(f"{path.name}.tmp.{os.getpid()}")
    try:
        frame.to_csv(temporary, sep="\t", index=False, lineterminator="\n")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _coerce_boolean(value: Any) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, str):
        normalized = value.strip().casefold()
        if normalized in {"true", "1"}:
            return True
        if normalized in {"false", "0"}:
            return False
    if isinstance(value, (int, np.integer)) and int(value) in {0, 1}:
        return bool(value)
    raise ValueError(f"value is not an explicit boolean: {value!r}")


def _inside_repo(path: str | Path, repo_root: str | Path, *, label: str) -> Path:
    root = Path(repo_root).resolve()
    resolved = Path(path).resolve()
    try:
        resolved.relative_to(root)
    except ValueError as error:
        raise ValueError(f"{label} must remain inside repo_root") from error
    if any("onedrive" in part.casefold() for part in resolved.parts):
        raise ValueError(f"{label} must not be inside OneDrive")
    return resolved


def _source_hashes(repo_root: Path) -> dict[str, str]:
    hashes: dict[str, str] = {}
    for relative in IMPLEMENTATION_SOURCE_PATHS:
        path = repo_root / relative
        if not path.is_file():
            raise ValueError(f"calibration implementation source is absent: {path}")
        hashes[relative] = sha256_file(path)
    return hashes


def _seed_digest(values: Sequence[int]) -> str:
    normalized = sorted({int(value) for value in values})
    return hashlib.sha256(
        json.dumps(normalized, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def seed_reservation_record(
    reservations: Mapping[str, Sequence[int]],
) -> dict[str, Any]:
    """Return a non-secret, exact digest record for forbidden RNG seeds."""

    if "production" not in reservations or "han" not in reservations:
        raise ValueError("seed reservations require production and Han sources")
    normalized = {
        str(label): sorted({int(value) for value in values})
        for label, values in reservations.items()
    }
    for label, values in normalized.items():
        if any(value <= 0 or value >= 2**31 - 1 for value in values):
            raise ValueError(f"{label} contains a seed outside the SLiM range")
    return {
        "sources": {
            label: {"n_unique": len(values), "sha256": _seed_digest(values)}
            for label, values in sorted(normalized.items())
        },
        "union_n_unique": len(set().union(*(set(x) for x in normalized.values()))),
        "union_sha256": _seed_digest(
            list(set().union(*(set(x) for x in normalized.values())))
        ),
    }


def collect_seed_reservations(campaign_dir: str | Path) -> dict[str, list[int]]:
    """Collect production and all fixed/fallback Han plan RNG streams.

    The named Han plan projections are required so an EAS plan cannot be
    generated before every intended Han simulation and panel stream has been
    atomically published. Result/ledger TSVs are deliberately ignored.
    """

    campaign = Path(campaign_dir).resolve()
    units_path = campaign / "execution_units.tsv"
    if not units_path.is_file():
        raise ValueError(f"execution-unit table is absent: {units_path}")
    units = pd.read_csv(units_path, sep="\t")
    if "seed" not in units:
        raise ValueError("execution-unit table has no seed column")
    production = [int(value) for value in units["seed"]]
    han: list[int] = []
    calibration_dir = campaign / "calibration"
    for name in HAN_SEED_PLAN_FILENAMES:
        path = calibration_dir / name
        if not path.is_file():
            raise ValueError(f"required Han calibration seed plan is absent: {path}")
        frame = pd.read_csv(path, sep="\t")
        if "seed" not in frame:
            raise ValueError(
                f"Han calibration plan has no simulation seed column: {path}"
            )
        if name == HAN_FALLBACK_PLAN_FILENAME and "panel_seed" not in frame:
            raise ValueError(f"Han fallback plan has no panel seed column: {path}")
        if name == HAN_FALLBACK_PLAN_FILENAME:
            manifest_path = campaign / HAN_FALLBACK_MANIFEST_RELATIVE_PATH
            try:
                manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError) as error:
                raise ValueError(
                    "authoritative Han fallback manifest is absent or unreadable"
                ) from error
            reservation = manifest.get("all_phase_seed_reservation")
            if (
                manifest.get("schema") != HAN_FALLBACK_MANIFEST_SCHEMA
                or manifest.get("status") != "planned"
                or not isinstance(reservation, Mapping)
            ):
                raise ValueError("authoritative Han fallback manifest is incompatible")
            relative_path = Path(str(reservation.get("path", "")))
            if relative_path.is_absolute():
                raise ValueError("Han fallback reservation path must be repo-relative")
            bound_path = (campaign.parent / relative_path).resolve()
            if (
                bound_path != path.resolve()
                or sha256_file(path) != str(reservation.get("sha256", ""))
                or int(reservation.get("rows", -1)) != len(frame)
                or reservation.get("simulation_seed_column") != "seed"
                or reservation.get("panel_seed_column") != "panel_seed"
                or reservation.get(
                    "authoritative_only_when_manifest_exists_and_binds_this_sha256"
                )
                is not True
            ):
                raise ValueError(
                    "Han fallback all-phase plan is not bound by its authoritative manifest"
                )
            if str(reservation.get("simulation_seed_sha256", "")) != _seed_digest(
                frame["seed"]
            ) or str(reservation.get("panel_seed_sha256", "")) != _seed_digest(
                frame["panel_seed"]
            ):
                raise ValueError("Han fallback seed-stream digests differ")
        for column in ("seed", "panel_seed"):
            if column not in frame:
                continue
            numeric = pd.to_numeric(frame[column], errors="raise").dropna()
            if not np.all(numeric.to_numpy(dtype=float) == np.floor(numeric)):
                raise ValueError(
                    f"Han calibration {column} values are not integers: {path}"
                )
            han.extend(int(value) for value in numeric)
    if not han:
        raise ValueError(
            "no Han calibration seed table is available; exact cross-phase "
            "seed-disjointness cannot be established"
        )
    return {"production": production, "han": han}


def _unit_for_row(row: Mapping[str, Any]) -> dict[str, Any]:
    target = float(row["target_allele_frequency"])
    return {
        "unit_id": str(row["calibration_id"]),
        "demography_id": "eas_phlash_median",
        "demography_kind": "no_introgression",
        "population": EAS_POPULATION,
        "source_model": "PHLASH EAS pointwise median trajectory",
        "selection_origin": "de_novo",
        "simulation_class": "selected",
        "selection_coefficient": float(row["selection_coefficient"]),
        "target_allele_frequency": target,
        "population_af_lower": float(row["population_af_lower"]),
        "population_af_upper": float(row["population_af_upper"]),
        "exact_sample_alt_count": int(row["exact_sample_alt_count"]),
        "sample_diploids": int(row["sample_diploids"]),
        "replicate_index": int(row["draw_index"]) + 1,
        "seed": int(row["seed"]),
        "sequence_length_bp": SEQUENCE_LENGTH_BP,
        "focal_position_bp": FOCAL_POSITION_BP,
    }


@lru_cache(maxsize=8)
def _cached_eas_median_demography(repo_root: str, resource_sha256: str) -> Any:
    """Load/build the immutable EAS median once per worker and resource digest."""

    root = Path(repo_root).resolve()
    if resource_sha256 != EAS_RESOURCE_SHA256:
        raise ValueError("EAS demographic resource checksum differs")
    artifact = load_phlash_eas_npz(
        root / EAS_RESOURCE_PATH,
        expected_sha256=resource_sha256,
    )
    return build_eas_demography_models(artifact)["median"]


def _load_eas_cell(
    repo_root: Path,
    *,
    selection_coefficient: float,
    target_frequency: float,
    af_half_width: float,
    slim_scaling_factor: float,
) -> tuple[Any, Any, Any, dict[str, Any], dict[str, Any]]:
    resource_digest = sha256_file(repo_root / EAS_RESOURCE_PATH)
    if resource_digest != EAS_RESOURCE_SHA256:
        raise ValueError("EAS demographic resource checksum differs")
    demographic = _cached_eas_median_demography(
        str(repo_root.resolve()), resource_digest
    )
    estimates = build_de_novo_origin_age_grid(
        {"q025": demographic, "median": demographic, "q975": demographic},
        selection_coefficients=(float(selection_coefficient),),
        target_frequency=float(target_frequency),
        slim_scaling_factor=float(slim_scaling_factor),
    )
    estimate = next(item for item in estimates if item.trajectory_label == "median")
    sweep = build_no_introgression_sweep_spec(
        origin_age_generations=float(estimate.age_generations),
        selection_coefficient=float(selection_coefficient),
        target_frequency=float(target_frequency),
        lower_frequency=float(target_frequency - af_half_width),
        upper_frequency=float(target_frequency + af_half_width),
    )
    origin = estimate.to_record()
    events = serialize_extended_events(sweep.extended_events)
    _validate_event_contract(
        events,
        coefficient=float(selection_coefficient),
        target=float(target_frequency),
        lower=float(target_frequency - af_half_width),
        upper=float(target_frequency + af_half_width),
        origin_age=float(estimate.age_generations),
    )
    return (
        demographic.stdpopsim_model,
        sweep,
        demographic,
        origin,
        {
            "events": events,
            "event_contract_sha256": _canonical_sha256(events),
        },
    )


def _event_time(event: Mapping[str, Any], key: str) -> float:
    value = event.get(key)
    if not isinstance(value, Mapping) or "generations_ago" not in value:
        raise ValueError(f"event {event.get('event_type')} lacks {key}")
    return float(value["generations_ago"])


def _validate_event_contract(
    events: Sequence[Mapping[str, Any]],
    *,
    coefficient: float,
    target: float,
    lower: float,
    upper: float,
    origin_age: float,
) -> None:
    """Fail closed if origin, continuous selection, or terminal bounds drift."""

    if len(events) != 5:
        raise ValueError("EAS sweep must contain exactly five focal events")
    draws = [event for event in events if event.get("event_type") == "DrawMutation"]
    fitness = [
        event for event in events if event.get("event_type") == "ChangeMutationFitness"
    ]
    conditions = [
        event
        for event in events
        if event.get("event_type") == "ConditionOnAlleleFrequency"
    ]
    if len(draws) != 1 or len(fitness) != 1 or len(conditions) != 3:
        raise ValueError("EAS sweep event topology differs from the frozen contract")
    draw = draws[0]
    if (
        str(draw.get("population")) != EAS_POPULATION
        or str(draw.get("single_site_id")) != FOCAL_SITE_ID
        or not math.isclose(_event_time(draw, "time"), origin_age, abs_tol=1e-12)
    ):
        raise ValueError("EAS de-novo mutation origin event differs")
    selected = fitness[0]
    if (
        str(selected.get("population")) != EAS_POPULATION
        or not math.isclose(
            _event_time(selected, "start_time"), origin_age, abs_tol=1e-12
        )
        or not math.isclose(_event_time(selected, "end_time"), 0.0, abs_tol=1e-12)
        or not math.isclose(
            float(selected.get("selection_coefficient", math.nan)),
            coefficient,
            abs_tol=1e-12,
        )
        or not math.isclose(
            float(selected.get("dominance_coefficient", math.nan)),
            DOMINANCE_COEFFICIENT,
            abs_tol=1e-12,
        )
    ):
        raise ValueError(
            "EAS selection is not the original continuous-to-present event"
        )
    nonloss = [
        event
        for event in conditions
        if event.get("operator") == ">"
        and math.isclose(float(event.get("allele_frequency", math.nan)), 0.0)
    ]
    terminal = [
        event
        for event in conditions
        if math.isclose(_event_time(event, "start_time"), 0.0, abs_tol=1e-12)
        and math.isclose(_event_time(event, "end_time"), 0.0, abs_tol=1e-12)
    ]
    terminal_values = {
        str(event.get("operator")): float(event.get("allele_frequency", math.nan))
        for event in terminal
    }
    if len(nonloss) != 1 or len(terminal) != 2:
        raise ValueError("EAS survival or terminal frequency conditions differ")
    if not (
        math.isclose(terminal_values.get(">=", math.nan), lower, abs_tol=1e-12)
        and math.isclose(terminal_values.get("<=", math.nan), upper, abs_tol=1e-12)
        and lower <= target <= upper
    ):
        raise ValueError("EAS lower/upper terminal AF bounds differ")


def _cell_contracts(
    design: EasSelectedCalibrationDesign, repo_root: Path
) -> dict[str, dict[str, Any]]:
    cells: dict[str, dict[str, Any]] = {}
    for coefficient in design.selection_coefficients:
        for target in design.target_frequencies:
            _, _, _, origin, event = _load_eas_cell(
                repo_root,
                selection_coefficient=coefficient,
                target_frequency=target,
                af_half_width=design.af_half_width,
                slim_scaling_factor=design.slim_scaling_factor,
            )
            origin_hash = _canonical_sha256(origin)
            record = {
                "demography_id": "eas_phlash_median",
                "population": EAS_POPULATION,
                "selection_origin": "de_novo",
                "selection_coefficient": float(coefficient),
                "target_allele_frequency": float(target),
                "population_af_lower_inclusive": float(target - design.af_half_width),
                "population_af_upper_inclusive": float(target + design.af_half_width),
                "origin_age_generations": float(origin["age_generations"]),
                "origin_contract": origin,
                "origin_contract_sha256": origin_hash,
                "extended_events": event["events"],
                "event_contract_sha256": event["event_contract_sha256"],
                "continuous_selection_to_present": True,
                "terminal_population_af_conditioning_in_slim": True,
                "gamma_smc_statistics_used": False,
            }
            record["cell_contract_sha256"] = _canonical_sha256(record)
            key = f"s={coefficient:.6f}|af={target:.6f}"
            cells[key] = record
    if len(cells) != 6:
        raise RuntimeError("EAS calibration must contain exactly six cells")
    return cells


def build_calibration_plans(
    design: EasSelectedCalibrationDesign,
    *,
    repo_root: str | Path,
    slim_path: str | Path,
    seed_reservations: Mapping[str, Sequence[int]],
) -> tuple[dict[str, pd.DataFrame], dict[str, Any]]:
    """Build all three immutable plans and their source/binary/event binding."""

    design.validate()
    root = Path(repo_root).resolve()
    slim = Path(slim_path).resolve()
    if not slim.is_file():
        raise ValueError(f"SLiM binary is absent: {slim}")
    resource = root / EAS_RESOURCE_PATH
    if sha256_file(resource) != EAS_RESOURCE_SHA256:
        raise ValueError("EAS demographic resource checksum differs")
    reservation_record = seed_reservation_record(seed_reservations)
    forbidden = set().union(
        *({int(value) for value in values} for values in seed_reservations.values())
    )
    cells = _cell_contracts(design, root)
    slim_digest = sha256_file(slim)
    plans: dict[str, pd.DataFrame] = {}
    all_stream_seeds: set[int] = set()
    for phase in PHASE_ORDER:
        rows: list[dict[str, Any]] = []
        for coefficient in design.selection_coefficients:
            for target in design.target_frequencies:
                cell = cells[f"s={coefficient:.6f}|af={target:.6f}"]
                for draw_index in range(PHASE_DRAWS[phase]):
                    identity = (
                        f"eas_selected__{phase}__s{coefficient:.6f}__"
                        f"af{target:.6f}__q{design.slim_scaling_factor:g}__"
                        f"draw{draw_index:03d}"
                    )
                    seed = _stable_seed(
                        int(design.base_seed),
                        f"{SCHEMA_VERSION}:{identity}:simulation",
                    )
                    panel_seed = _stable_seed(
                        int(design.base_seed),
                        f"{SCHEMA_VERSION}:{identity}:panel",
                    )
                    for stream, value in (
                        ("simulation", seed),
                        ("panel", panel_seed),
                    ):
                        if value in forbidden:
                            raise ValueError(
                                f"EAS calibration {stream} seed collides with a "
                                f"reserved seed: {value}"
                            )
                        if value in all_stream_seeds:
                            raise ValueError(
                                f"EAS calibration {stream} seed is duplicated: {value}"
                            )
                        all_stream_seeds.add(value)
                    sequence_length = (
                        SEQUENCE_LENGTH_BP
                        if phase == PHASE_SENSITIVITY
                        else FOCAL_ONLY_LENGTH_BP
                    )
                    focal_position = (
                        FOCAL_POSITION_BP
                        if phase == PHASE_SENSITIVITY
                        else FOCAL_ONLY_POSITION_BP
                    )
                    rows.append(
                        {
                            "calibration_id": identity,
                            "phase": phase,
                            "processing_mode": PHASE_MODES[phase],
                            "selection_coefficient": float(coefficient),
                            "target_allele_frequency": float(target),
                            "population_af_lower": float(target - design.af_half_width),
                            "population_af_upper": float(target + design.af_half_width),
                            "candidate_pool_diploids": int(
                                design.candidate_pool_diploids
                            ),
                            "sample_diploids": int(design.sample_diploids),
                            "exact_sample_alt_count": round(
                                2 * design.sample_diploids * target
                            ),
                            "draw_index": int(draw_index),
                            "seed": int(seed),
                            "panel_seed": int(panel_seed),
                            "slim_scaling_factor": float(design.slim_scaling_factor),
                            "slim_burn_in": float(design.slim_burn_in),
                            "sequence_length_bp": int(sequence_length),
                            "focal_position_bp": int(focal_position),
                            "origin_age_generations": float(
                                cell["origin_age_generations"]
                            ),
                            "origin_contract_sha256": str(
                                cell["origin_contract_sha256"]
                            ),
                            "event_contract_sha256": str(cell["event_contract_sha256"]),
                            "cell_contract_sha256": str(cell["cell_contract_sha256"]),
                            "slim_sha256": slim_digest,
                        }
                    )
        frame = pd.DataFrame(rows, columns=PLAN_COLUMNS)
        _validate_phase_plan(frame, phase)
        plans[phase] = frame
    manifest = {
        "schema": PLAN_MANIFEST_SCHEMA,
        "status": "planned",
        "scientific_purpose": (
            "acceptance-rate validation only; no selection, origin-age, "
            "demography, endpoint, or Gamma-SMC tuning"
        ),
        "design": design.to_record(),
        "phase_order": list(PHASE_ORDER),
        "phase_modes": dict(PHASE_MODES),
        "cells": cells,
        "source_sha256": _source_hashes(root),
        "eas_resource": {
            "path": EAS_RESOURCE_PATH,
            "sha256": EAS_RESOURCE_SHA256,
        },
        "slim": {
            "path_at_plan_time": str(slim),
            "sha256": slim_digest,
        },
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "pyslim": pyslim.__version__,
            "stdpopsim": stdpopsim.__version__,
            "tskit": tskit.__version__,
        },
        "seed_reservations_at_plan_time": reservation_record,
        "seed_derivation": {
            "algorithm": "sha256 label mapped to 1..2^31-2",
            "simulation_label": f"{SCHEMA_VERSION}:identity:simulation",
            "panel_label": f"{SCHEMA_VERSION}:identity:panel",
            "streams_per_draw": 2,
            "all_streams_reserved_at_plan_time": True,
        },
        "plans": {
            phase: {
                "rows": len(frame),
                "cells": int(
                    frame.groupby(
                        ["selection_coefficient", "target_allele_frequency"]
                    ).ngroups
                ),
                "content_sha256": _canonical_sha256(frame.to_dict(orient="records")),
            }
            for phase, frame in plans.items()
        },
    }
    manifest["contract_sha256"] = _canonical_sha256(manifest)
    return plans, manifest


def _validate_phase_plan(plan: pd.DataFrame, phase: str) -> None:
    if phase not in PHASE_ORDER:
        raise ValueError(f"unknown EAS calibration phase: {phase}")
    if list(plan.columns) != list(PLAN_COLUMNS):
        raise ValueError("EAS calibration plan columns are incompatible")
    expected_rows = 6 * PHASE_DRAWS[phase]
    if len(plan) != expected_rows:
        raise ValueError(f"{phase} requires exactly {expected_rows} draws")
    if (
        plan["calibration_id"].duplicated().any()
        or plan["seed"].duplicated().any()
        or plan["panel_seed"].duplicated().any()
    ):
        raise ValueError(f"{phase} plan contains duplicate IDs or RNG seeds")
    simulation_seeds = {int(value) for value in plan["seed"]}
    panel_seeds = {int(value) for value in plan["panel_seed"]}
    if simulation_seeds.intersection(panel_seeds):
        raise ValueError(f"{phase} simulation and panel RNG streams overlap")
    if any(
        value <= 0 or value >= 2**31 - 1
        for value in simulation_seeds.union(panel_seeds)
    ):
        raise ValueError(f"{phase} contains an RNG seed outside the SLiM range")
    if set(plan["phase"].astype(str)) != {phase}:
        raise ValueError(f"{phase} plan contains another phase")
    if set(plan["processing_mode"].astype(str)) != {PHASE_MODES[phase]}:
        raise ValueError(f"{phase} processing mode differs")
    groups = plan.groupby(["selection_coefficient", "target_allele_frequency"])
    if groups.ngroups != 6 or set(groups.size()) != {PHASE_DRAWS[phase]}:
        raise ValueError(f"{phase} does not have exact six-cell draw coverage")
    for _, group in groups:
        if sorted(int(value) for value in group["draw_index"]) != list(
            range(PHASE_DRAWS[phase])
        ):
            raise ValueError(f"{phase} draw indexes are not contiguous from zero")
    expected_length = (
        SEQUENCE_LENGTH_BP if phase == PHASE_SENSITIVITY else FOCAL_ONLY_LENGTH_BP
    )
    expected_focal = (
        FOCAL_POSITION_BP if phase == PHASE_SENSITIVITY else FOCAL_ONLY_POSITION_BP
    )
    if set(plan["sequence_length_bp"].astype(int)) != {expected_length} or set(
        plan["focal_position_bp"].astype(int)
    ) != {expected_focal}:
        raise ValueError(f"{phase} contig contract differs")


def write_calibration_plans(
    output_dir: str | Path,
    plans: Mapping[str, pd.DataFrame],
    manifest: Mapping[str, Any],
) -> dict[str, Path]:
    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}
    updated = json.loads(json.dumps(manifest))
    for phase in PHASE_ORDER:
        plan = plans[phase]
        _validate_phase_plan(plan, phase)
        path = output / f"{phase}_plan.tsv"
        _atomic_frame(path, plan)
        serialized = pd.read_csv(path, sep="\t")
        _validate_phase_plan(serialized, phase)
        updated["plans"][phase]["path"] = path.name
        updated["plans"][phase]["sha256"] = sha256_file(path)
        updated["plans"][phase]["content_sha256"] = _canonical_sha256(
            serialized.to_dict(orient="records")
        )
        written[phase] = path
    updated.pop("contract_sha256", None)
    updated["contract_sha256"] = _canonical_sha256(updated)
    manifest_path = output / "eas_selected_calibration_plan.json"
    _atomic_json(manifest_path, updated)
    written["manifest"] = manifest_path
    return written


def read_plan_bundle(
    output_dir: str | Path,
) -> tuple[dict[str, pd.DataFrame], dict[str, Any]]:
    output = Path(output_dir).resolve()
    manifest_path = output / "eas_selected_calibration_plan.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    expected_hash = manifest.get("contract_sha256")
    unhashed = dict(manifest)
    unhashed.pop("contract_sha256", None)
    if manifest.get(
        "schema"
    ) != PLAN_MANIFEST_SCHEMA or expected_hash != _canonical_sha256(unhashed):
        raise ValueError("EAS calibration plan manifest is incompatible or corrupt")
    plans: dict[str, pd.DataFrame] = {}
    all_stream_seeds: set[int] = set()
    for phase in PHASE_ORDER:
        entry = manifest["plans"][phase]
        path = output / str(entry["path"])
        if sha256_file(path) != str(entry["sha256"]):
            raise ValueError(f"{phase} plan checksum differs")
        plan = pd.read_csv(path, sep="\t")
        _validate_phase_plan(plan, phase)
        content_hash = _canonical_sha256(plan.to_dict(orient="records"))
        if content_hash != str(entry["content_sha256"]):
            raise ValueError(f"{phase} plan content binding differs")
        seeds = {
            int(value) for column in ("seed", "panel_seed") for value in plan[column]
        }
        if seeds.intersection(all_stream_seeds):
            raise ValueError("EAS calibration phase RNG streams overlap")
        all_stream_seeds.update(seeds)
        plans[phase] = plan
    return plans, manifest


def _validate_manifest_integrity(manifest: Mapping[str, Any]) -> None:
    unhashed = dict(manifest)
    expected = unhashed.pop("contract_sha256", None)
    if manifest.get("schema") != PLAN_MANIFEST_SCHEMA or expected != _canonical_sha256(
        unhashed
    ):
        raise ValueError("EAS calibration plan manifest is incompatible or corrupt")


def _validate_plan_manifest_binding(
    plan: pd.DataFrame, phase: str, manifest: Mapping[str, Any]
) -> None:
    _validate_manifest_integrity(manifest)
    _validate_phase_plan(plan, phase)
    entry = manifest.get("plans", {}).get(phase)
    if not isinstance(entry, Mapping):
        raise TypeError(f"{phase} plan manifest entry is absent")
    if _canonical_sha256(plan.to_dict(orient="records")) != str(
        entry.get("content_sha256")
    ):
        raise ValueError(f"{phase} plan differs from its manifest binding")
    cells = manifest.get("cells")
    if not isinstance(cells, Mapping):
        raise TypeError("EAS plan manifest cell contracts are absent")
    for row in plan.itertuples(index=False):
        key = (
            f"s={float(row.selection_coefficient):.6f}|"
            f"af={float(row.target_allele_frequency):.6f}"
        )
        cell = cells.get(key)
        if not isinstance(cell, Mapping) or any(
            str(getattr(row, field)) != str(cell[field])
            for field in (
                "origin_contract_sha256",
                "event_contract_sha256",
                "cell_contract_sha256",
            )
        ):
            raise ValueError(f"{phase} row differs from scientific cell contract")


def _manifest_current_binding(
    manifest: Mapping[str, Any], repo_root: Path, slim_path: Path
) -> None:
    _validate_manifest_integrity(manifest)
    if manifest.get("source_sha256") != _source_hashes(repo_root):
        raise ValueError("EAS calibration implementation sources differ from plan")
    resource = repo_root / EAS_RESOURCE_PATH
    if sha256_file(resource) != str(manifest["eas_resource"]["sha256"]):
        raise ValueError("EAS demographic resource differs from plan")
    if sha256_file(slim_path) != str(manifest["slim"]["sha256"]):
        raise ValueError("SLiM binary differs from the EAS calibration plan")
    design = EasSelectedCalibrationDesign(
        selection_coefficients=tuple(
            float(value) for value in manifest["design"]["selection_coefficients"]
        ),
        target_frequencies=tuple(
            float(value) for value in manifest["design"]["target_allele_frequencies"]
        ),
        slim_scaling_factor=float(manifest["design"]["slim_scaling_factor"]),
        slim_burn_in=float(manifest["design"]["slim_burn_in"]),
        candidate_pool_diploids=int(manifest["design"]["candidate_pool_diploids"]),
        sample_diploids=int(manifest["design"]["sample_diploids"]),
        af_half_width=float(manifest["design"]["af_half_width"]),
        minimum_hit_rate=float(manifest["design"]["minimum_hit_rate"]),
        base_seed=int(manifest["design"]["base_seed"]),
    )
    current_cells = _cell_contracts(design, repo_root)
    if current_cells != manifest.get("cells"):
        raise ValueError("EAS origin-age or extended-event contracts differ from plan")


def validate_seed_disjointness(
    plans: Mapping[str, pd.DataFrame],
    reservations: Mapping[str, Sequence[int]],
) -> dict[str, Any]:
    """Audit all EAS phase seeds against current production and Han seeds."""

    reservation_record = seed_reservation_record(reservations)
    forbidden = set().union(
        *({int(value) for value in values} for values in reservations.values())
    )
    observed: set[int] = set()
    simulation: set[int] = set()
    panel: set[int] = set()
    for phase in PHASE_ORDER:
        _validate_phase_plan(plans[phase], phase)
        phase_simulation = {int(value) for value in plans[phase]["seed"]}
        phase_panel = {int(value) for value in plans[phase]["panel_seed"]}
        seeds = phase_simulation.union(phase_panel)
        overlap = seeds.intersection(forbidden)
        if overlap:
            preview = ", ".join(str(value) for value in sorted(overlap)[:3])
            raise ValueError(
                f"{phase} collides with production/Han seed reservations: {preview}"
            )
        if seeds.intersection(observed):
            raise ValueError("EAS calibration phases have overlapping seeds")
        observed.update(seeds)
        simulation.update(phase_simulation)
        panel.update(phase_panel)
    return {
        **reservation_record,
        "eas_calibration_n_unique": len(observed),
        "eas_calibration_sha256": _seed_digest(list(observed)),
        "eas_simulation_n_unique": len(simulation),
        "eas_simulation_sha256": _seed_digest(list(simulation)),
        "eas_panel_n_unique": len(panel),
        "eas_panel_sha256": _seed_digest(list(panel)),
        "disjoint": True,
    }


def _normalized_row(row: Mapping[str, Any]) -> dict[str, Any]:
    return {
        key: (value.item() if isinstance(value, np.generic) else value)
        for key, value in row.items()
    }


def _draw_contract(
    row: Mapping[str, Any],
    *,
    repo_root: Path,
    slim_path: Path,
    manifest: Mapping[str, Any],
) -> dict[str, Any]:
    phase = str(row["phase"])
    return {
        "schema": DRAW_CONTRACT_SCHEMA,
        "plan_manifest_contract_sha256": str(manifest["contract_sha256"]),
        "plan_row": {key: _normalized_row(row)[key] for key in PLAN_COLUMNS},
        "scientific_contract": {
            "mutation_origin": "prespecified_Q5_EAS_de_novo_origin_age",
            "selection": "continuous_from_origin_to_present",
            "terminal_population_af_conditions": "lower_and_upper_inclusive_in_slim",
            "external_ascertainment": (
                "500_diploid_pool_band_then_exact_k_N100_genotype_QC"
            ),
            "gamma_smc_statistics_used": False,
        },
        "processing": {
            "mode": PHASE_MODES[phase],
            "sequence_length_bp": int(row["sequence_length_bp"]),
            "focal_position_bp": int(row["focal_position_bp"]),
            "stdpopsim_recapitation_applied": phase == PHASE_SENSITIVITY,
            "post_slim_neutral_mutation_overlay_applied": (phase == PHASE_SENSITIVITY),
            "production_tree_use_permitted": phase == PHASE_SENSITIVITY,
        },
        # Current files are checked once at phase entry/audit by
        # _manifest_current_binding.  Reusing the immutable manifest digests
        # here avoids re-hashing multi-megabyte sources for every one of the
        # 840 draws while retaining an exact per-draw source binding.
        "sources": dict(manifest["source_sha256"]),
        "eas_resource_sha256": str(manifest["eas_resource"]["sha256"]),
        "slim": {
            "path": str(slim_path),
            "sha256": str(manifest["slim"]["sha256"]),
        },
        "software": {
            "python": platform.python_version(),
            "pyslim": pyslim.__version__,
            "stdpopsim": stdpopsim.__version__,
            "tskit": tskit.__version__,
        },
    }


def _raw_focal_identity_and_counts(
    ts: tskit.TreeSequence, expectation: Mapping[str, Any]
) -> tuple[dict[str, Any], np.ndarray]:
    """Validate and count the type-1 focal mutation before post-processing."""

    coordinate = float(expectation["focal_position_0based"])
    focal_sites = [
        site
        for site in ts.sites()
        if math.isclose(float(site.position), coordinate, rel_tol=0.0, abs_tol=1e-9)
    ]
    entries: list[tuple[tskit.Mutation, Mapping[str, Any]]] = []
    for site in focal_sites:
        for mutation in site.mutations:
            metadata = mutation.metadata
            if not isinstance(metadata, Mapping):
                continue
            mutation_list = metadata.get("mutation_list")
            if not isinstance(mutation_list, list):
                continue
            entries.extend(
                (mutation, entry)
                for entry in mutation_list
                if isinstance(entry, Mapping)
                and int(entry.get("mutation_type", -1))
                == int(expectation["expected_mutation_type"])
            )
    if len(focal_sites) != 1 or len(entries) != 1:
        raise focused_simulation.SelectedFocalIdentityError(
            "raw EAS focal identity requires one type-1 nonrecurrent mutation"
        )
    mutation, entry = entries[0]
    if mutation.parent != tskit.NULL:
        raise focused_simulation.SelectedFocalIdentityError(
            "raw EAS focal mutation must be a root mutation"
        )
    if int(entry.get("subpopulation", -1)) != int(expectation["source_subpopulation"]):
        raise focused_simulation.SelectedFocalIdentityError(
            "raw EAS focal source population differs"
        )
    expected_coefficient = float(expectation["expected_selection_coeff"])
    if not math.isclose(
        float(entry.get("selection_coeff", math.nan)),
        expected_coefficient,
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        raise focused_simulation.SelectedFocalIdentityError(
            "raw EAS focal mutation type coefficient differs"
        )
    q = float(expectation["slim_scaling_factor"])
    observed_age = float(mutation.time) * q
    if not math.isclose(
        observed_age,
        float(expectation["realized_origin_time_generations"]),
        rel_tol=0.0,
        abs_tol=1e-9,
    ):
        raise focused_simulation.SelectedFocalIdentityError(
            "raw EAS focal origin age differs"
        )
    try:
        observed_q = float(ts.metadata["SLiM"]["user_metadata"]["Q"][0])
    except (KeyError, TypeError, ValueError, IndexError) as error:
        raise focused_simulation.SelectedFocalIdentityError(
            "raw EAS focal tree lacks Q metadata"
        ) from error
    if not math.isclose(observed_q, q, rel_tol=0.0, abs_tol=1e-12):
        raise focused_simulation.SelectedFocalIdentityError(
            "raw EAS focal tree Q differs"
        )
    mutation_node = int(mutation.node)
    tree = ts.at(coordinate)
    carriers = {int(node) for node in tree.samples(mutation_node)}
    samples = {int(node) for node in ts.samples()}
    counts: list[int] = []
    assigned: list[int] = []
    for individual in ts.individuals():
        nodes = [int(node) for node in individual.nodes if int(node) in samples]
        if not nodes:
            continue
        if len(nodes) != 2:
            raise ValueError("remembered EAS calibration individual is not diploid")
        assigned.extend(nodes)
        counts.append(sum(node in carriers for node in nodes))
    if not counts or set(assigned) != samples or len(assigned) != len(set(assigned)):
        raise ValueError("remembered EAS calibration samples are not exact pairs")
    identity = {
        "status": "valid",
        "focal_position_0based": coordinate,
        "expected_mutation_type": int(expectation["expected_mutation_type"]),
        "observed_mutation_node": mutation_node,
        "observed_source_subpopulation": int(entry["subpopulation"]),
        "observed_origin_time_generations": observed_age,
        "observed_slim_scaling_factor": observed_q,
    }
    return identity, np.asarray(counts, dtype=np.int8)


def _result_base(row: Mapping[str, Any]) -> dict[str, Any]:
    normalized = _normalized_row(row)
    return {key: normalized[key] for key in PLAN_COLUMNS}


def _failed_result(
    row: Mapping[str, Any],
    *,
    status: str,
    elapsed_seconds: float,
    error: str,
) -> dict[str, Any]:
    phase = str(row["phase"])
    return {
        **_result_base(row),
        "status": status,
        "evaluable": False,
        "terminal_class": status,
        "candidate_pool_alt_count": None,
        "candidate_pool_af": None,
        "candidate_pool_frequency_gate_passed": False,
        "exact_panel_feasible": False,
        "pool_band_and_panel_hit": False,
        "sample_alt_count": None,
        "sample_af": None,
        "n_hom_ref": None,
        "n_heterozygous": None,
        "n_hom_alt": None,
        "stdpopsim_recapitation_applied": phase == PHASE_SENSITIVITY,
        "post_slim_neutral_mutation_overlay_applied": phase == PHASE_SENSITIVITY,
        "elapsed_seconds": float(elapsed_seconds),
        "trajectory_sha256": "",
        "panel_tree_sha256": "",
        "panel_manifest_sha256": "",
        "error": str(error),
    }


def _run_draw(task: Mapping[str, Any]) -> dict[str, Any]:
    row = _normalized_row(task["row"])
    root = Path(task["repo_root"]).resolve()
    work_root = Path(task["work_root"]).resolve()
    slim = Path(task["slim_path"]).resolve()
    timeout_seconds = float(task["timeout_seconds"])
    manifest = dict(task["manifest"])
    phase = str(row["phase"])
    directory = work_root / phase / str(row["calibration_id"])
    directory.mkdir(parents=True, exist_ok=True)
    completion_path = directory / "completion.json"
    contract = _draw_contract(row, repo_root=root, slim_path=slim, manifest=manifest)
    contract_hash = _canonical_sha256(contract)
    with exclusive_unit_lock(
        directory / ".calibration.lock",
        unit_id=str(row["calibration_id"]),
        phase=f"eas_selected_{phase}",
    ):
        if completion_path.is_file():
            payload = json.loads(completion_path.read_text(encoding="utf-8"))
            if (
                payload.get("schema") != SCHEMA_VERSION
                or payload.get("contract_sha256") != contract_hash
                or payload.get("contract") != contract
            ):
                raise ValueError("cached EAS calibration draw is incompatible")
            _validate_result_artifact_checksums(directory, payload["result"])
            return dict(payload["result"])
        started = perf_counter()
        trajectory_path = directory / "trajectory.csv"
        try:
            model, sweep, _, origin, event = _load_eas_cell(
                root,
                selection_coefficient=float(row["selection_coefficient"]),
                target_frequency=float(row["target_allele_frequency"]),
                af_half_width=(
                    float(row["target_allele_frequency"])
                    - float(row["population_af_lower"])
                ),
                slim_scaling_factor=float(row["slim_scaling_factor"]),
            )
            if _canonical_sha256(origin) != str(row["origin_contract_sha256"]):
                raise ValueError("draw origin-age contract differs from plan")
            if event["event_contract_sha256"] != str(row["event_contract_sha256"]):
                raise ValueError("draw extended-event contract differs from plan")
            if phase != PHASE_SENSITIVITY:
                sweep = replace(
                    sweep,
                    contig=build_focal_contig(
                        sequence_length=FOCAL_ONLY_LENGTH_BP,
                        focal_position=FOCAL_ONLY_POSITION_BP,
                    ),
                )
            samples = {EAS_POPULATION: int(row["candidate_pool_diploids"])}
            expectation = _selected_focal_expectation(
                model, sweep, float(row["slim_scaling_factor"])
            )
            expectation["focal_position_0based"] = int(row["focal_position_bp"])
            engine = stdpopsim.get_engine("slim")
            simulate_kwargs = {
                "seed": int(row["seed"]),
                "extended_events": list(sweep.extended_events),
                "slim_path": str(slim),
                "slim_scaling_factor": float(row["slim_scaling_factor"]),
                "slim_burn_in": float(row["slim_burn_in"]),
                "logfile": str(trajectory_path),
                "logfile_interval": 100,
            }
            with _wall_clock_timeout(timeout_seconds):
                if phase == PHASE_SENSITIVITY:
                    ts = engine.simulate(
                        model,
                        sweep.contig,
                        samples,
                        keep_mutation_ids_as_alleles=False,
                        **simulate_kwargs,
                    )
                    identity = inspect_selected_focal_identity(ts, expectation)
                    if identity["status"] != "valid":
                        raise focused_simulation.SelectedFocalIdentityError(
                            f"10-Mb selected focal identity is {identity['status']}"
                        )
                    pool_table = build_diploid_pair_table(ts, FOCAL_POSITION_BP)
                    counts = pool_table["focal_selected_allele_count"].to_numpy(
                        dtype=int
                    )
                else:
                    raw = engine.simulate(
                        model,
                        sweep.contig,
                        samples,
                        keep_mutation_ids_as_alleles=True,
                        _recap_and_rescale=False,
                        **simulate_kwargs,
                    )
                    ts = engine._simplify_remembered(raw)
                    identity, counts = _raw_focal_identity_and_counts(ts, expectation)
            if len(counts) != int(row["candidate_pool_diploids"]):
                raise ValueError("EAS calibration candidate pool has wrong size")
            unit = _unit_for_row(row)
            gate = candidate_pool_frequency_gate(unit, counts)
            triples = (
                feasible_exact_panel_genotype_triples(
                    counts,
                    sample_diploids=int(row["sample_diploids"]),
                    exact_alt_count=int(row["exact_sample_alt_count"]),
                )
                if gate["passed"]
                else []
            )
            hit = bool(gate["passed"] and triples)
            sample_counts: np.ndarray | None = None
            panel_tree_hash = ""
            panel_manifest_hash = ""
            if phase == PHASE_SENSITIVITY and hit:
                panel_ts, _ = select_exact_sample_panel(
                    ts,
                    sample_diploids=int(row["sample_diploids"]),
                    exact_alt_count=int(row["exact_sample_alt_count"]),
                    seed=int(row["panel_seed"]),
                )
                validate_selected_focal_identity(panel_ts, expectation)
                panel_table = build_diploid_pair_table(panel_ts, FOCAL_POSITION_BP)
                sample_counts = panel_table["focal_selected_allele_count"].to_numpy(
                    dtype=int
                )
                panel_tree = directory / "accepted_panel.trees"
                panel_manifest = directory / "accepted_panel_manifest.tsv"
                panel_ts.dump(panel_tree)
                _atomic_frame(panel_manifest, panel_table)
                panel_tree_hash = sha256_file(panel_tree)
                panel_manifest_hash = sha256_file(panel_manifest)
            alt_count = int(np.sum(counts))
            terminal_class = (
                "loss"
                if alt_count == 0
                else "fixed"
                if alt_count == 2 * len(counts)
                else "segregating"
            )
            trajectory_hash = (
                sha256_file(trajectory_path) if trajectory_path.is_file() else ""
            )
            result = {
                **_result_base(row),
                "status": "complete",
                "evaluable": True,
                "terminal_class": terminal_class,
                "candidate_pool_alt_count": alt_count,
                "candidate_pool_af": float(alt_count / (2 * len(counts))),
                "candidate_pool_frequency_gate_passed": bool(gate["passed"]),
                "exact_panel_feasible": bool(triples),
                "pool_band_and_panel_hit": hit,
                "sample_alt_count": (
                    int(sample_counts.sum()) if sample_counts is not None else None
                ),
                "sample_af": (
                    float(sample_counts.mean() / 2)
                    if sample_counts is not None
                    else None
                ),
                "n_hom_ref": (
                    int(np.count_nonzero(sample_counts == 0))
                    if sample_counts is not None
                    else None
                ),
                "n_heterozygous": (
                    int(np.count_nonzero(sample_counts == 1))
                    if sample_counts is not None
                    else None
                ),
                "n_hom_alt": (
                    int(np.count_nonzero(sample_counts == 2))
                    if sample_counts is not None
                    else None
                ),
                "stdpopsim_recapitation_applied": phase == PHASE_SENSITIVITY,
                "post_slim_neutral_mutation_overlay_applied": (
                    phase == PHASE_SENSITIVITY
                ),
                "elapsed_seconds": float(perf_counter() - started),
                "trajectory_sha256": trajectory_hash,
                "panel_tree_sha256": panel_tree_hash,
                "panel_manifest_sha256": panel_manifest_hash,
                "error": "",
            }
        except SimulationDrawTimeout as error:
            result = _failed_result(
                row,
                status="timed_out",
                elapsed_seconds=perf_counter() - started,
                error=str(error),
            )
        except Exception as error:  # noqa: BLE001 - preserve exact per-draw failure
            result = _failed_result(
                row,
                status="failed",
                elapsed_seconds=perf_counter() - started,
                error=f"{type(error).__name__}: {error}",
            )
        payload = {
            "schema": SCHEMA_VERSION,
            "contract": contract,
            "contract_sha256": contract_hash,
            "result": result,
        }
        if result["evaluable"]:
            _atomic_json(completion_path, payload)
            (directory / "last_failure.json").unlink(missing_ok=True)
        else:
            # Failures are deliberately not cached as completion: rerunning uses
            # the same immutable seed and can resume after a larger time budget.
            _atomic_json(directory / "last_failure.json", payload)
        return result


def _validate_result_artifact_checksums(
    directory: Path, result: Mapping[str, Any]
) -> None:
    expected = {
        "trajectory.csv": str(result.get("trajectory_sha256", "")),
        "accepted_panel.trees": str(result.get("panel_tree_sha256", "")),
        "accepted_panel_manifest.tsv": str(result.get("panel_manifest_sha256", "")),
    }
    for name, digest in expected.items():
        path = directory / name
        if digest:
            if not path.is_file() or sha256_file(path) != digest:
                raise ValueError(f"cached EAS calibration artifact differs: {path}")
        elif path.is_file() and name != "trajectory.csv":
            raise ValueError(f"unexpected cached EAS panel artifact exists: {path}")


def run_phase(
    plan: pd.DataFrame,
    *,
    repo_root: str | Path,
    work_root: str | Path,
    slim_path: str | Path,
    manifest: Mapping[str, Any],
    workers: int = 4,
    timeout_seconds: float = DEFAULT_DRAW_TIMEOUT_SECONDS,
) -> pd.DataFrame:
    """Run/resume one immutable phase; each draw has an independent timeout."""

    phase_values = set(plan["phase"].astype(str))
    if len(phase_values) != 1:
        raise ValueError("run_phase requires one phase")
    phase = next(iter(phase_values))
    _validate_plan_manifest_binding(plan, phase, manifest)
    if int(workers) < 1:
        raise ValueError("workers must be positive")
    if not math.isfinite(timeout_seconds) or timeout_seconds <= 0:
        raise ValueError("draw timeout must be finite and positive")
    root = Path(repo_root).resolve()
    work = _inside_repo(work_root, root, label="EAS calibration work root")
    slim = Path(slim_path).resolve()
    _manifest_current_binding(manifest, root, slim)
    tasks = [
        {
            "row": row,
            "repo_root": str(root),
            "work_root": str(work),
            "slim_path": str(slim),
            "timeout_seconds": float(timeout_seconds),
            "manifest": dict(manifest),
        }
        for row in plan.to_dict(orient="records")
    ]
    rows: list[dict[str, Any]] = []
    if int(workers) == 1:
        rows = [_run_draw(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=int(workers)) as pool:
            futures = {pool.submit(_run_draw, task): task for task in tasks}
            for future in as_completed(futures):
                task = futures[future]
                try:
                    rows.append(future.result())
                except Exception as error:  # noqa: BLE001 - preserve worker failure
                    rows.append(
                        _failed_result(
                            task["row"],
                            status="failed",
                            elapsed_seconds=math.nan,
                            error=f"{type(error).__name__}: {error}",
                        )
                    )
    return (
        pd.DataFrame(rows, columns=RESULT_COLUMNS)
        .sort_values("calibration_id")
        .reset_index(drop=True)
    )


def _values_equal(left: Any, right: Any) -> bool:
    if pd.isna(left) and pd.isna(right):
        return True
    if isinstance(left, (bool, np.bool_)) or isinstance(right, (bool, np.bool_)):
        try:
            return _coerce_boolean(left) == _coerce_boolean(right)
        except ValueError:
            return False
    if isinstance(left, (int, float, np.number)) and isinstance(
        right, (int, float, np.number)
    ):
        return math.isclose(float(left), float(right), rel_tol=0.0, abs_tol=1e-12)
    return str(left) == str(right)


def validate_phase_ledger(
    plan: pd.DataFrame, ledger: pd.DataFrame, phase: str
) -> pd.DataFrame:
    """Validate exact coverage and biological/result invariants for one phase."""

    _validate_phase_plan(plan, phase)
    missing_columns = [column for column in RESULT_COLUMNS if column not in ledger]
    if missing_columns:
        raise ValueError(
            "EAS calibration ledger is missing columns: " + ", ".join(missing_columns)
        )
    if ledger["calibration_id"].duplicated().any():
        raise ValueError(f"{phase} ledger contains duplicate calibration IDs")
    if set(ledger["calibration_id"].astype(str)) != set(
        plan["calibration_id"].astype(str)
    ):
        raise ValueError(f"{phase} ledger does not cover its plan exactly")
    canonical = (
        ledger.loc[:, RESULT_COLUMNS]
        .sort_values("calibration_id")
        .reset_index(drop=True)
    )
    indexed_plan = plan.set_index("calibration_id", drop=False)
    for row in canonical.to_dict(orient="records"):
        identity = str(row["calibration_id"])
        planned = indexed_plan.loc[identity]
        for column in PLAN_COLUMNS:
            if not _values_equal(row[column], planned[column]):
                raise ValueError(
                    f"{phase} ledger differs from plan for {identity}: {column}"
                )
        if str(row["status"]) != "complete" or not _coerce_boolean(row["evaluable"]):
            raise ValueError(f"{phase} has an unevaluable planned draw: {identity}")
        gate = _coerce_boolean(row["candidate_pool_frequency_gate_passed"])
        feasible = _coerce_boolean(row["exact_panel_feasible"])
        hit = _coerce_boolean(row["pool_band_and_panel_hit"])
        if hit != (gate and feasible):
            raise ValueError(f"{phase} hit indicator is inconsistent: {identity}")
        pool_diploids = int(row["candidate_pool_diploids"])
        alt_count = int(row["candidate_pool_alt_count"])
        if not 0 <= alt_count <= 2 * pool_diploids or not math.isclose(
            float(row["candidate_pool_af"]),
            alt_count / (2 * pool_diploids),
            rel_tol=0.0,
            abs_tol=1e-12,
        ):
            raise ValueError(f"{phase} candidate-pool count/AF differs: {identity}")
        expected_recap = phase == PHASE_SENSITIVITY
        if _coerce_boolean(row["stdpopsim_recapitation_applied"]) != expected_recap:
            raise ValueError(f"{phase} recapitation flag differs: {identity}")
        if (
            _coerce_boolean(row["post_slim_neutral_mutation_overlay_applied"])
            != expected_recap
        ):
            raise ValueError(f"{phase} mutation-overlay flag differs: {identity}")
        panel_hash = str(row["panel_tree_sha256"])
        manifest_hash = str(row["panel_manifest_sha256"])
        if expected_recap and hit:
            counts = [
                int(row["n_hom_ref"]),
                int(row["n_heterozygous"]),
                int(row["n_hom_alt"]),
            ]
            if (
                sum(counts) != int(row["sample_diploids"])
                or counts[1] + 2 * counts[2] != int(row["exact_sample_alt_count"])
                or counts[0] < 2
                or counts[1] < 1
                or counts[2] < 2
                or int(row["sample_alt_count"]) != int(row["exact_sample_alt_count"])
                or not math.isclose(
                    float(row["sample_af"]),
                    int(row["exact_sample_alt_count"])
                    / (2 * int(row["sample_diploids"])),
                    rel_tol=0.0,
                    abs_tol=1e-12,
                )
            ):
                raise ValueError(
                    f"{phase} accepted panel does not have exact AF/QC: {identity}"
                )
            if len(panel_hash) != 64 or len(manifest_hash) != 64:
                raise ValueError(
                    f"{phase} accepted panel checksums are absent: {identity}"
                )
        elif panel_hash or manifest_hash:
            raise ValueError(f"{phase} has unexpected accepted-panel files: {identity}")
    return canonical


def audit_phase_completions(
    plan: pd.DataFrame,
    ledger: pd.DataFrame,
    *,
    repo_root: str | Path,
    work_root: str | Path,
    slim_path: str | Path,
    manifest: Mapping[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Rebuild a phase ledger from checksum-validated per-draw completions."""

    phase = str(plan["phase"].iloc[0])
    _validate_plan_manifest_binding(plan, phase, manifest)
    canonical_ledger = validate_phase_ledger(plan, ledger, phase)
    root = Path(repo_root).resolve()
    work = _inside_repo(work_root, root, label="EAS calibration work root")
    slim = Path(slim_path).resolve()
    _manifest_current_binding(manifest, root, slim)
    supplied = canonical_ledger.set_index("calibration_id", drop=False)
    completion_results: list[dict[str, Any]] = []
    audit_rows: list[dict[str, Any]] = []
    for planned in plan.to_dict(orient="records"):
        identity = str(planned["calibration_id"])
        directory = work / phase / identity
        completion = directory / "completion.json"
        if not completion.is_file():
            raise ValueError(f"{phase} completion is absent: {completion}")
        payload = json.loads(completion.read_text(encoding="utf-8"))
        contract = _draw_contract(
            planned, repo_root=root, slim_path=slim, manifest=manifest
        )
        if (
            payload.get("schema") != SCHEMA_VERSION
            or payload.get("contract") != contract
            or payload.get("contract_sha256") != _canonical_sha256(contract)
        ):
            raise ValueError(f"{phase} completion contract differs: {identity}")
        result = payload.get("result")
        if not isinstance(result, Mapping):
            raise TypeError(f"{phase} completion result is malformed: {identity}")
        missing = [column for column in RESULT_COLUMNS if column not in result]
        if missing:
            raise ValueError(f"{phase} completion result is incomplete: {identity}")
        _validate_result_artifact_checksums(directory, result)
        for column in RESULT_COLUMNS:
            if not _values_equal(result[column], supplied.loc[identity, column]):
                raise ValueError(
                    f"{phase} ledger differs from completion for {identity}: {column}"
                )
        completion_results.append({column: result[column] for column in RESULT_COLUMNS})
        audit_rows.append(
            {
                "calibration_id": identity,
                "completion_path": completion.relative_to(root).as_posix(),
                "completion_sha256": sha256_file(completion),
                "trajectory_sha256": str(result["trajectory_sha256"]),
                "panel_tree_sha256": str(result["panel_tree_sha256"]),
                "panel_manifest_sha256": str(result["panel_manifest_sha256"]),
            }
        )
    rebuilt = (
        pd.DataFrame(completion_results, columns=RESULT_COLUMNS)
        .sort_values("calibration_id")
        .reset_index(drop=True)
    )
    validate_phase_ledger(plan, rebuilt, phase)
    audit = (
        pd.DataFrame(audit_rows).sort_values("calibration_id").reset_index(drop=True)
    )
    return rebuilt, audit


def summarize_phase(
    plan: pd.DataFrame, ledger: pd.DataFrame, phase: str
) -> pd.DataFrame:
    canonical = validate_phase_ledger(plan, ledger, phase)
    rows: list[dict[str, Any]] = []
    for keys, group in canonical.groupby(
        ["selection_coefficient", "target_allele_frequency"], sort=True
    ):
        coefficient, target = (float(value) for value in keys)
        hits = group["pool_band_and_panel_hit"].map(_coerce_boolean)
        rows.append(
            {
                "phase": phase,
                "selection_coefficient": coefficient,
                "target_allele_frequency": target,
                "n_planned": len(group),
                "n_evaluable": int(group["evaluable"].map(_coerce_boolean).sum()),
                "n_pool_band_and_panel_hits": int(hits.sum()),
                "pool_band_and_panel_hit_rate": float(hits.mean()),
                "minimum_required_hit_rate": DEFAULT_MINIMUM_HIT_RATE,
                "passed": bool(hits.mean() >= DEFAULT_MINIMUM_HIT_RATE),
                "median_candidate_pool_af": float(
                    pd.to_numeric(group["candidate_pool_af"]).median()
                ),
                "total_elapsed_seconds": float(
                    pd.to_numeric(group["elapsed_seconds"]).sum()
                ),
                "origin_contract_sha256": str(group["origin_contract_sha256"].iloc[0]),
                "event_contract_sha256": str(group["event_contract_sha256"].iloc[0]),
                "cell_contract_sha256": str(group["cell_contract_sha256"].iloc[0]),
            }
        )
    summary = (
        pd.DataFrame(rows)
        .sort_values(
            ["selection_coefficient", "target_allele_frequency"],
            ascending=[False, True],
        )
        .reset_index(drop=True)
    )
    if len(summary) != 6 or not summary["passed"].map(_coerce_boolean).all():
        # A caller can still persist this useful summary, but freezing will
        # fail closed.  Do not raise here solely for an inadequate biological
        # hit rate.
        pass
    return summary


def _production_seed_digest(execution_units: pd.DataFrame) -> str:
    if "seed" not in execution_units:
        raise ValueError("execution-unit table has no seed column")
    return _seed_digest([int(value) for value in execution_units["seed"]])


def _relative_artifact_record(paths: Sequence[Path], repo_root: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    for path in paths:
        resolved = path.resolve()
        try:
            relative = resolved.relative_to(repo_root).as_posix()
        except ValueError as error:
            raise ValueError(
                "frozen calibration artifact lies outside repo_root"
            ) from error
        records[relative] = sha256_file(resolved)
    return records


def freeze_calibration(
    plans: Mapping[str, pd.DataFrame],
    ledgers: Mapping[str, pd.DataFrame],
    manifest: Mapping[str, Any],
    *,
    repo_root: str | Path,
    output_dir: str | Path,
    work_root: str | Path,
    slim_path: str | Path,
    execution_units_path: str | Path,
    seed_reservations: Mapping[str, Sequence[int]],
) -> dict[str, Any]:
    """Audit all phases and write an EAS production authorization if adequate."""

    root = Path(repo_root).resolve()
    output = _inside_repo(output_dir, root, label="EAS calibration output")
    work = _inside_repo(work_root, root, label="EAS calibration work root")
    units_path = _inside_repo(execution_units_path, root, label="execution-unit table")
    slim = Path(slim_path).resolve()
    output.mkdir(parents=True, exist_ok=True)
    _manifest_current_binding(manifest, root, slim)
    seed_audit = validate_seed_disjointness(plans, seed_reservations)
    execution_units = pd.read_csv(units_path, sep="\t")
    selected_eas = execution_units[
        execution_units["demography_id"].astype(str).eq("eas_phlash_median")
        & execution_units["simulation_class"].astype(str).eq("selected")
    ]
    production_cells = {
        (float(row.selection_coefficient), float(row.target_allele_frequency))
        for row in selected_eas.itertuples(index=False)
    }
    expected_cells = {
        (float(coefficient), float(target))
        for coefficient in DEFAULT_SELECTION_COEFFICIENTS
        for target in DEFAULT_TARGET_FREQUENCIES
    }
    if production_cells != expected_cells:
        raise ValueError("execution-unit table lacks the exact six EAS selected cells")
    artifact_paths: list[Path] = [units_path]
    plan_manifest_path = output / "eas_selected_calibration_plan.json"
    if not plan_manifest_path.is_file():
        raise ValueError("EAS calibration plan manifest is absent from output_dir")
    artifact_paths.append(plan_manifest_path)
    summaries: dict[str, pd.DataFrame] = {}
    phase_records: dict[str, Any] = {}
    frozen_path = output / "eas_selected_calibration_frozen.json"
    not_frozen_path = output / "eas_selected_calibration_not_frozen.json"
    all_adequate = True
    cross_phase_cells: dict[tuple[float, float], tuple[str, str, str]] = {}
    for phase in PHASE_ORDER:
        plan = plans[phase]
        supplied_ledger = ledgers[phase]
        rebuilt, audit = audit_phase_completions(
            plan,
            supplied_ledger,
            repo_root=root,
            work_root=work,
            slim_path=slim,
            manifest=manifest,
        )
        summary = summarize_phase(plan, rebuilt, phase)
        summaries[phase] = summary
        plan_path = output / f"{phase}_plan.tsv"
        ledger_path = output / f"{phase}_ledger.tsv"
        summary_path = output / f"{phase}_summary.tsv"
        audit_path = output / f"{phase}_completion_checksums.tsv"
        if not plan_path.is_file() or sha256_file(plan_path) != str(
            manifest["plans"][phase]["sha256"]
        ):
            raise ValueError(f"{phase} frozen plan file differs from manifest")
        _atomic_frame(ledger_path, rebuilt)
        _atomic_frame(summary_path, summary)
        _atomic_frame(audit_path, audit)
        artifact_paths.extend((plan_path, ledger_path, summary_path, audit_path))
        phase_passed = bool(
            len(summary) == 6
            and summary["passed"].map(_coerce_boolean).all()
            and (summary["n_evaluable"] == summary["n_planned"]).all()
            and set(summary["n_planned"].astype(int)) == {PHASE_DRAWS[phase]}
        )
        all_adequate = all_adequate and phase_passed
        phase_records[phase] = {
            "processing_mode": PHASE_MODES[phase],
            "draws_per_cell": PHASE_DRAWS[phase],
            "total_draws": len(rebuilt),
            "minimum_hit_rate": DEFAULT_MINIMUM_HIT_RATE,
            "all_six_cells_passed": phase_passed,
            "summary": summary.to_dict(orient="records"),
        }
        for row in summary.itertuples(index=False):
            key = (float(row.selection_coefficient), float(row.target_allele_frequency))
            hashes = (
                str(row.origin_contract_sha256),
                str(row.event_contract_sha256),
                str(row.cell_contract_sha256),
            )
            previous = cross_phase_cells.setdefault(key, hashes)
            if previous != hashes:
                raise ValueError("origin/event cell contract differs across phases")
    if set(cross_phase_cells) != expected_cells:
        raise ValueError("frozen EAS calibration lacks exact six-cell coverage")
    authorizations: dict[str, Any] = {}
    for key, cell in manifest["cells"].items():
        coefficient = float(cell["selection_coefficient"])
        target = float(cell["target_allele_frequency"])
        authorizations[key] = {
            "selection_coefficient": coefficient,
            "target_allele_frequency": target,
            "origin_age_generations": float(cell["origin_age_generations"]),
            "origin_contract_sha256": str(cell["origin_contract_sha256"]),
            "event_contract_sha256": str(cell["event_contract_sha256"]),
            "cell_contract_sha256": str(cell["cell_contract_sha256"]),
            "phase_hit_rates": {
                phase: float(
                    summaries[phase]
                    .loc[
                        np.isclose(
                            summaries[phase]["selection_coefficient"], coefficient
                        )
                        & np.isclose(
                            summaries[phase]["target_allele_frequency"], target
                        ),
                        "pool_band_and_panel_hit_rate",
                    ]
                    .iloc[0]
                )
                for phase in PHASE_ORDER
            },
        }
    payload = {
        "schema": FROZEN_SCHEMA_VERSION,
        "status": "frozen" if all_adequate else "not_frozen",
        "gamma_smc_statistics_used": False,
        "minimum_pool_band_and_panel_hit_rate": DEFAULT_MINIMUM_HIT_RATE,
        "scientific_caveat": (
            "Terminal lower and upper population-AF conditions make these "
            "explicitly conditional selected trajectories; focal-only phases "
            "are acceptance proxies, not production genealogies."
        ),
        "plan_manifest_contract_sha256": str(manifest["contract_sha256"]),
        "source_sha256": dict(manifest["source_sha256"]),
        "slim": {
            "sha256": sha256_file(slim),
            "slim_scaling_factor": DEFAULT_SLIM_SCALING_FACTOR,
            "slim_burn_in": DEFAULT_SLIM_BURN_IN,
        },
        "eas_resource": dict(manifest["eas_resource"]),
        "seed_disjointness_at_freeze": seed_audit,
        "execution_units": {
            "path": units_path.relative_to(root).as_posix(),
            "sha256": sha256_file(units_path),
            "production_seed_digest": _production_seed_digest(execution_units),
            "eas_selected_cells": 6,
        },
        "production_contract": {
            "authorization_scope": "eas_phlash_median_selected_only",
            "integration_status": "not_wired_into_campaign_executor",
            "source_paths_permitted_to_change_only_for_explicit_integration": list(
                PRODUCTION_INTEGRATION_MUTABLE_SOURCE_PATHS
            ),
            "selection_origin": "de_novo_at_prespecified_Q5_origin_age",
            "selection_episode": "continuous_from_origin_to_present",
            "terminal_population_af_conditioning_in_slim": True,
            "terminal_conditions": "lower_and_upper_inclusive_target_plus_minus_0.025",
            "candidate_pool_diploids": DEFAULT_POOL_DIPLOIDS,
            "external_candidate_pool_af_gate_retained": True,
            "exact_k_sample_panel_retained": True,
            "sample_diploids": DEFAULT_SAMPLE_DIPLOIDS,
            "genotype_qc": {
                "minimum_hom_ref": 2,
                "minimum_het": 1,
                "minimum_hom_alt": 2,
            },
            "focal_only_tree_production_use_permitted": False,
            "ten_mb_sensitivity_required": True,
            "ten_mb_sensitivity_passed": bool(
                phase_records[PHASE_SENSITIVITY]["all_six_cells_passed"]
            ),
        },
        "phases": phase_records,
        "cell_authorizations": authorizations,
        "artifacts": _relative_artifact_record(artifact_paths, root),
    }
    payload["payload_sha256"] = _canonical_sha256(payload)
    if payload["status"] == "frozen":
        _atomic_json(frozen_path, payload)
        not_frozen_path.unlink(missing_ok=True)
    else:
        frozen_path.unlink(missing_ok=True)
        _atomic_json(not_frozen_path, payload)
    return payload


def load_frozen_eas_selected_calibration(
    path: str | Path,
    *,
    repo_root: str | Path,
    execution_units_path: str | Path,
    slim_path: str | Path,
) -> dict[str, Any]:
    """Fail-closed loader for a future explicit production integration."""

    root = Path(repo_root).resolve()
    frozen_path = _inside_repo(path, root, label="frozen EAS calibration")
    units_path = _inside_repo(execution_units_path, root, label="execution-unit table")
    payload = json.loads(frozen_path.read_text(encoding="utf-8"))
    digest = payload.get("payload_sha256")
    unhashed = dict(payload)
    unhashed.pop("payload_sha256", None)
    if (
        payload.get("schema") != FROZEN_SCHEMA_VERSION
        or payload.get("status") != "frozen"
        or digest != _canonical_sha256(unhashed)
    ):
        raise ValueError("frozen EAS calibration is incompatible or corrupt")
    stored_sources = payload.get("source_sha256")
    if not isinstance(stored_sources, Mapping):
        raise TypeError("frozen EAS calibration source binding is absent")
    current_sources = _source_hashes(root)
    immutable_sources = set(IMPLEMENTATION_SOURCE_PATHS).difference(
        PRODUCTION_INTEGRATION_MUTABLE_SOURCE_PATHS
    )
    if any(
        str(stored_sources.get(relative)) != current_sources[relative]
        for relative in immutable_sources
    ):
        raise ValueError("frozen EAS immutable calibration source binding differs")
    if sha256_file(slim_path) != str(payload["slim"]["sha256"]):
        raise ValueError("frozen EAS calibration SLiM binary differs")
    if sha256_file(root / EAS_RESOURCE_PATH) != str(payload["eas_resource"]["sha256"]):
        raise ValueError("frozen EAS demographic resource differs")
    if sha256_file(units_path) != str(payload["execution_units"]["sha256"]):
        raise ValueError("frozen EAS execution-unit table differs")
    units = pd.read_csv(units_path, sep="\t")
    if _production_seed_digest(units) != str(
        payload["execution_units"]["production_seed_digest"]
    ):
        raise ValueError("frozen EAS production seed digest differs")
    artifacts = payload.get("artifacts")
    if not isinstance(artifacts, Mapping) or not artifacts:
        raise ValueError("frozen EAS artifact manifest is absent")
    for relative, expected in artifacts.items():
        artifact = _inside_repo(root / str(relative), root, label="EAS artifact")
        if not artifact.is_file() or sha256_file(artifact) != str(expected):
            raise ValueError(f"frozen EAS artifact checksum differs: {relative}")
    contract = payload.get("production_contract", {})
    required_true = (
        "terminal_population_af_conditioning_in_slim",
        "external_candidate_pool_af_gate_retained",
        "exact_k_sample_panel_retained",
        "ten_mb_sensitivity_required",
        "ten_mb_sensitivity_passed",
    )
    if any(contract.get(key) is not True for key in required_true):
        raise ValueError("frozen EAS production authorization is incomplete")
    if contract.get("integration_status") != "not_wired_into_campaign_executor":
        raise ValueError("unexpected EAS production integration state")
    cells = payload.get("cell_authorizations")
    if not isinstance(cells, Mapping) or len(cells) != 6:
        raise ValueError("frozen EAS calibration does not authorize six cells")
    # Production wiring necessarily edits focused_selection_simulation.py.  Its
    # calibration-time digest remains in the frozen provenance, while the
    # scientific origin/event hashes are recomputed from the still-immutable
    # model/calibration sources.  This permits explicit wiring without allowing
    # mutation ages, selection intervals, or terminal conditions to drift.
    current_cells = _cell_contracts(EasSelectedCalibrationDesign(), root)
    for key, current in current_cells.items():
        authorization = cells.get(key)
        if not isinstance(authorization, Mapping) or any(
            str(authorization.get(field)) != str(current[field])
            for field in (
                "origin_contract_sha256",
                "event_contract_sha256",
                "cell_contract_sha256",
            )
        ):
            raise ValueError(f"frozen EAS scientific cell contract differs: {key}")
    for phase in PHASE_ORDER:
        record = payload["phases"].get(phase, {})
        if (
            record.get("all_six_cells_passed") is not True
            or int(record.get("draws_per_cell", -1)) != PHASE_DRAWS[phase]
        ):
            raise ValueError(f"frozen EAS phase authorization differs: {phase}")
    return payload


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("plan", "run", "validate", "freeze"))
    parser.add_argument("--phase", choices=PHASE_ORDER)
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--campaign-dir", type=Path)
    parser.add_argument("--slim-bin", type=Path)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--timeout-minutes", type=float, default=5.0)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    root = args.repo_root.resolve()
    campaign = _inside_repo(
        args.campaign_dir or root / "focused_selection_EAS_sim",
        root,
        label="focused campaign directory",
    )
    output = campaign / "calibration" / "eas_selected"
    work = campaign / "work" / "eas_selected_calibration"
    units_path = campaign / "execution_units.tsv"
    reservations = collect_seed_reservations(campaign)
    if args.slim_bin is None:
        raise ValueError(f"{args.action} requires --slim-bin for binary binding")
    slim = args.slim_bin.resolve()
    if args.action == "plan":
        plans, manifest = build_calibration_plans(
            EasSelectedCalibrationDesign(),
            repo_root=root,
            slim_path=slim,
            seed_reservations=reservations,
        )
        write_calibration_plans(output, plans, manifest)
        return 0
    plans, manifest = read_plan_bundle(output)
    _manifest_current_binding(manifest, root, slim)
    validate_seed_disjointness(plans, reservations)
    if args.action in {"run", "validate"} and args.phase is None:
        raise ValueError(f"{args.action} requires --phase")
    if args.action == "run":
        ledger = run_phase(
            plans[args.phase],
            repo_root=root,
            work_root=work,
            slim_path=slim,
            manifest=manifest,
            workers=args.workers,
            timeout_seconds=args.timeout_minutes * 60.0,
        )
        _atomic_frame(output / f"{args.phase}_ledger.tsv", ledger)
        return 0 if ledger["evaluable"].map(_coerce_boolean).all() else 1
    if args.action == "validate":
        ledger = pd.read_csv(output / f"{args.phase}_ledger.tsv", sep="\t")
        rebuilt, audit = audit_phase_completions(
            plans[args.phase],
            ledger,
            repo_root=root,
            work_root=work,
            slim_path=slim,
            manifest=manifest,
        )
        summary = summarize_phase(plans[args.phase], rebuilt, args.phase)
        _atomic_frame(output / f"{args.phase}_ledger.tsv", rebuilt)
        _atomic_frame(output / f"{args.phase}_summary.tsv", summary)
        _atomic_frame(output / f"{args.phase}_completion_checksums.tsv", audit)
        return 0 if summary["passed"].map(_coerce_boolean).all() else 1
    ledgers = {
        phase: pd.read_csv(output / f"{phase}_ledger.tsv", sep="\t")
        for phase in PHASE_ORDER
    }
    payload = freeze_calibration(
        plans,
        ledgers,
        manifest,
        repo_root=root,
        output_dir=output,
        work_root=work,
        slim_path=slim,
        execution_units_path=units_path,
        seed_reservations=reservations,
    )
    return 0 if payload["status"] == "frozen" else 1


__all__ = [
    "PHASE_CONFIRM",
    "PHASE_ORDER",
    "PHASE_SCREEN",
    "PHASE_SENSITIVITY",
    "EasSelectedCalibrationDesign",
    "audit_phase_completions",
    "build_calibration_plans",
    "collect_seed_reservations",
    "freeze_calibration",
    "load_frozen_eas_selected_calibration",
    "main",
    "read_plan_bundle",
    "run_phase",
    "seed_reservation_record",
    "summarize_phase",
    "validate_phase_ledger",
    "validate_seed_disjointness",
    "write_calibration_plans",
]
