"""Independent Han calibration with terminal AF conditioning inside SLiM.

This fallback deliberately does not reuse observations from the failed
fixed-cessation calibration. The only contact with those artifacts is seed
exclusion and immutable failed-model provenance. Selection begins at the
introgression pulse and remains active through the present; no cessation time
is fitted and no Gamma-SMC statistic is available to any gate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import shutil
import subprocess
from collections.abc import Mapping, Sequence
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, replace
from pathlib import Path
from time import perf_counter
from typing import Any

import numpy as np
import pandas as pd
import stdpopsim
import tskit

from . import focused_selection_simulation as focused_simulation
from .eas_sweep_analysis import build_diploid_pair_table
from .eas_sweep_models import (
    DOMINANCE_COEFFICIENT,
    FOCAL_SITE_ID,
    HAN_SPLIT_GENERATIONS,
    INTROGRESSION_PULSE_GENERATIONS,
    INTROGRESSION_TARGET_POPULATION,
    PRE_PULSE_SOURCE_CHECK_GENERATIONS,
    build_focal_contig,
    serialize_extended_events,
)
from .eas_sweep_study import BASE_SEED
from .focused_selection_simulation import (
    DEFAULT_SLIM_BURN_IN,
    DEFAULT_NEUTRAL_ANCESTRY_BATCH_SIZE,
    DEFAULT_NEUTRAL_MAX_ATTEMPTS,
    DEFAULT_SELECTED_MAX_DRAWS,
    INTROGRESSION_MUTATION_AGE_GENERATIONS,
    INTROGRESSION_SOURCE_PRESENCE_EPSILON,
    SimulationDrawTimeout,
    _atomic_tree,
    _canonical_sha256,
    _selected_focal_expectation,
    _stable_seed,
    _wall_clock_timeout,
    build_raw_han_selected_sweep,
    candidate_pool_frequency_gate,
    exclusive_unit_lock,
    feasible_exact_panel_genotype_triples,
    select_exact_sample_panel,
    sha256_file,
    validate_selected_focal_identity,
)

MANIFEST_SCHEMA = "gamma-smc.han-frequency-conditioned-manifest/v1"
DRAW_SCHEMA = "gamma-smc.han-frequency-conditioned-draw/v1"
AUDIT_SCHEMA = "gamma-smc.han-frequency-conditioned-phase-audit/v1"
FROZEN_SCHEMA = "gamma-smc.han-frequency-conditioned-calibration-frozen/v1"
SEED_REGISTRY_SCHEMA = "gamma-smc.han-frequency-conditioned-seeds/v1"

SELECTION_COEFFICIENTS = (0.01, 0.005)
TARGET_FREQUENCIES = (0.10, 0.20, 0.30)
AF_HALF_WIDTH = 0.025
SLIM_SCALING_FACTOR = 5.0
CANDIDATE_POOL_DIPLOIDS = 500
SAMPLE_DIPLOIDS = 100
MINIMUM_HIT_RATE = 0.20
DEFAULT_DRAW_TIMEOUT_SECONDS = 1_800.0
EXPECTED_STDPOPSIM_SLIM_VERSION = "4.2.2"
FOCAL_ONLY_LENGTH_BP = 2
FOCAL_ONLY_POSITION_BP = 1
PRODUCTION_LENGTH_BP = 10_000_000
PRODUCTION_POSITION_BP = 5_000_000
PHASE_DRAW_COUNTS = {
    "screen": 20,
    "confirmation": 100,
    "sensitivity_10mb": 20,
}
FOCAL_ONLY_MODE = "focal_2bp_raw_forward_no_recap_no_neutral_overlay"
PRODUCTION_MODE = (
    "production_10mb_normal_recap_neutral_overlay_strict_identity_exact_panel"
)
MODEL_MODE = "selection_through_present_with_internal_terminal_han_af_conditioning"

ROW_COLUMNS = (
    "draw_id",
    "phase",
    "selection_coefficient",
    "target_allele_frequency",
    "population_af_lower",
    "population_af_upper",
    "exact_sample_alt_count",
    "sample_diploids",
    "candidate_pool_diploids",
    "slim_scaling_factor",
    "mutation_age_generations",
    "selection_start_generations_ago",
    "selection_end_generations_ago",
    "draw_index",
    "seed",
    "seed_nonce",
    "panel_seed",
    "panel_seed_nonce",
    "timeout_seconds",
    "sequence_length_bp",
    "focal_position_bp",
    "processing_mode",
)


@dataclass(frozen=True)
class FrequencyConditionedDesign:
    selection_coefficients: tuple[float, ...] = SELECTION_COEFFICIENTS
    target_frequencies: tuple[float, ...] = TARGET_FREQUENCIES
    af_half_width: float = AF_HALF_WIDTH
    slim_scaling_factor: float = SLIM_SCALING_FACTOR
    candidate_pool_diploids: int = CANDIDATE_POOL_DIPLOIDS
    sample_diploids: int = SAMPLE_DIPLOIDS
    minimum_hit_rate: float = MINIMUM_HIT_RATE
    draw_timeout_seconds: float = DEFAULT_DRAW_TIMEOUT_SECONDS
    base_seed: int = BASE_SEED

    def validate(self) -> None:
        if tuple(self.selection_coefficients) != SELECTION_COEFFICIENTS:
            raise ValueError("fallback selection-coefficient grid is frozen")
        if tuple(self.target_frequencies) != TARGET_FREQUENCIES:
            raise ValueError("fallback target-frequency grid is frozen")
        if not math.isclose(self.af_half_width, AF_HALF_WIDTH, abs_tol=1e-12):
            raise ValueError("fallback AF half width is frozen")
        if not math.isclose(
            self.slim_scaling_factor, SLIM_SCALING_FACTOR, abs_tol=1e-12
        ):
            raise ValueError("fallback Q is frozen at 5")
        if int(self.candidate_pool_diploids) != CANDIDATE_POOL_DIPLOIDS:
            raise ValueError("fallback candidate pool is frozen at 500 diploids")
        if int(self.sample_diploids) != SAMPLE_DIPLOIDS:
            raise ValueError("fallback panel is frozen at 100 diploids")
        if not math.isclose(self.minimum_hit_rate, MINIMUM_HIT_RATE, abs_tol=1e-12):
            raise ValueError("fallback minimum hit rate is frozen at 20%")
        if not math.isclose(
            self.draw_timeout_seconds,
            DEFAULT_DRAW_TIMEOUT_SECONDS,
            abs_tol=1e-12,
        ):
            raise ValueError("fallback per-draw timeout is frozen at 1800 seconds")


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


def _require_sha256(value: str, *, label: str) -> str:
    normalized = str(value).strip().casefold()
    if len(normalized) != 64 or any(
        character not in "0123456789abcdef" for character in normalized
    ):
        raise ValueError(f"{label} is not a SHA-256 digest")
    return normalized


def _seed_digest(values: Sequence[int]) -> str:
    seeds = sorted(int(value) for value in values)
    return hashlib.sha256(
        json.dumps(seeds, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


_ROW_STRING_COLUMNS = {"draw_id", "phase", "processing_mode"}
_ROW_INTEGER_COLUMNS = {
    "exact_sample_alt_count",
    "sample_diploids",
    "candidate_pool_diploids",
    "draw_index",
    "seed",
    "seed_nonce",
    "panel_seed",
    "panel_seed_nonce",
    "sequence_length_bp",
    "focal_position_bp",
}


def _normalized_manifest_row(row: Mapping[str, Any]) -> dict[str, Any]:
    if set(ROW_COLUMNS).difference(row):
        raise ValueError("fallback draw record lacks manifest fields")
    normalized = {}
    for column in ROW_COLUMNS:
        value = row[column]
        if column in _ROW_STRING_COLUMNS:
            normalized[column] = str(value)
        elif column in _ROW_INTEGER_COLUMNS:
            normalized[column] = int(value)
        else:
            normalized[column] = float(value)
    return normalized


def _implementation_paths(repo_root: Path) -> tuple[Path, ...]:
    return (
        Path(__file__).resolve(),
        repo_root / "python/gamma_smc_aou/focused_selection_simulation.py",
        repo_root / "python/gamma_smc_aou/eas_sweep_models.py",
        repo_root / "python/gamma_smc_aou/eas_sweep_analysis.py",
        repo_root / "scripts/calibrate_han_frequency_conditioned.py",
        repo_root / "pyproject.toml",
        repo_root / "uv.lock",
    )


def _implementation_sources(repo_root: Path) -> dict[str, str]:
    sources = {}
    for path in _implementation_paths(repo_root):
        resolved = path.resolve()
        if not resolved.is_file():
            raise ValueError(f"required implementation source is absent: {resolved}")
        try:
            key = resolved.relative_to(repo_root.resolve()).as_posix()
        except ValueError as error:
            raise ValueError("implementation source lies outside repo_root") from error
        sources[key] = sha256_file(resolved)
    return sources


def _require_repo_path(path: str | Path, repo_root: str | Path, *, label: str) -> Path:
    resolved = Path(path).resolve()
    root = Path(repo_root).resolve()
    try:
        resolved.relative_to(root)
    except ValueError as error:
        raise ValueError(f"fallback {label} must remain inside repo_root") from error
    if any("onedrive" in part.casefold() for part in resolved.parts):
        raise ValueError(f"fallback {label} must not be inside OneDrive")
    return resolved


def _slim_version(slim_path: Path) -> str:
    if not slim_path.is_file():
        raise ValueError(f"fallback SLiM executable is absent: {slim_path}")
    try:
        completed = subprocess.run(
            [str(slim_path), "-v"],
            check=True,
            capture_output=True,
            text=True,
            timeout=30,
        )
    except (OSError, subprocess.CalledProcessError, subprocess.TimeoutExpired) as error:
        raise ValueError("fallback SLiM version probe failed") from error
    version = " ".join((completed.stdout or completed.stderr).split())
    if "SLiM" not in version or EXPECTED_STDPOPSIM_SLIM_VERSION not in version:
        raise ValueError(
            "fallback requires the stdpopsim-compatible SLiM "
            f"{EXPECTED_STDPOPSIM_SLIM_VERSION} runtime"
        )
    return version


def _runtime_record(repo_root: Path, slim_path: str | Path) -> dict[str, Any]:
    slim = _require_repo_path(slim_path, repo_root, label="SLiM executable")
    if not slim.is_file():
        raise ValueError(f"fallback SLiM executable is absent: {slim}")
    return {
        "python": platform.python_version(),
        "stdpopsim": stdpopsim.__version__,
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "tskit": tskit.__version__,
        "slim": {
            "path": slim.relative_to(repo_root.resolve()).as_posix(),
            "sha256": sha256_file(slim),
            "version": _slim_version(slim),
            "expected_version": EXPECTED_STDPOPSIM_SLIM_VERSION,
        },
    }


def _artifact_evidence(
    paths: Sequence[str | Path], repo_root: Path, *, label: str
) -> list[dict[str, Any]]:
    if not paths:
        raise ValueError(f"{label} artifacts are required")
    records = []
    observed = set()
    for raw in paths:
        path = _require_repo_path(raw, repo_root, label=label)
        relative = path.relative_to(repo_root.resolve()).as_posix()
        if relative in observed:
            raise ValueError(f"duplicate {label} artifact: {relative}")
        observed.add(relative)
        if not path.is_file():
            raise ValueError(f"{label} artifact is absent: {path}")
        records.append(
            {
                "path": relative,
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
            }
        )
    return sorted(records, key=lambda item: item["path"])


def _excluded_seed_evidence(
    paths: Sequence[str | Path], repo_root: Path
) -> tuple[set[int], list[dict[str, Any]]]:
    if not paths:
        raise ValueError("excluded seed artifacts are required")
    excluded: set[int] = set()
    evidence = []
    observed_paths = set()
    for raw_path in paths:
        path = Path(raw_path).resolve()
        try:
            relative = path.relative_to(repo_root.resolve()).as_posix()
        except ValueError as error:
            raise ValueError("excluded seed artifact lies outside repo_root") from error
        if relative in observed_paths:
            raise ValueError(f"duplicate excluded seed artifact: {relative}")
        observed_paths.add(relative)
        if not path.is_file():
            raise ValueError(f"excluded seed artifact is absent: {path}")
        try:
            frame = pd.read_csv(path, sep="\t")
        except Exception as error:
            raise ValueError(f"excluded seed artifact is unreadable: {path}") from error
        if "seed" not in frame.columns or frame.empty:
            raise ValueError(f"excluded seed artifact lacks seeds: {path}")
        values = [int(value) for value in frame["seed"]]
        base_unique = sorted(set(values))
        expanded = set(base_unique)
        expansion_mode = "direct_rng_seed"
        if {
            "unit_id",
            "simulation_class",
            "seed",
        }.issubset(frame.columns) and path.name == "execution_units.tsv":
            # Production unit seeds are namespaces, not usually direct RNG
            # invocations. Reserve every bounded derived seed used by the
            # production selected and neutral executors as well as the base
            # namespace itself. This makes "disjoint" an RNG-level claim.
            expansion_mode = "focused_campaign_bounded_rng_namespace_v1"
            for row in frame.to_dict(orient="records"):
                base = int(row["seed"])
                simulation_class = str(row["simulation_class"])
                if simulation_class == "selected":
                    for attempt in range(DEFAULT_SELECTED_MAX_DRAWS):
                        expanded.add(_stable_seed(base, f"selected:{attempt}"))
                        expanded.add(_stable_seed(base, f"selected-panel:{attempt}"))
                elif simulation_class == "neutral":
                    for attempt in range(DEFAULT_NEUTRAL_MAX_ATTEMPTS):
                        for component in ("ancestry", "mutations", "panel"):
                            expanded.add(
                                _stable_seed(base, f"neutral:{attempt}:{component}")
                            )
                    batches = math.ceil(
                        DEFAULT_NEUTRAL_MAX_ATTEMPTS
                        / DEFAULT_NEUTRAL_ANCESTRY_BATCH_SIZE
                    )
                    for batch in range(batches):
                        expanded.add(
                            _stable_seed(base, f"neutral-batch:{batch}:candidate")
                        )
                else:
                    raise ValueError(
                        "execution-unit seed evidence has an unknown class"
                    )
        unique = sorted(expanded)
        excluded.update(unique)
        evidence.append(
            {
                "path": relative,
                "sha256": sha256_file(path),
                "rows": int(len(frame)),
                "seed_column": "seed",
                "base_unique_seeds": len(base_unique),
                "base_seed_sha256": _seed_digest(base_unique),
                "rng_expansion_mode": expansion_mode,
                "reserved_rng_seeds": len(unique),
                "reserved_rng_seed_sha256": _seed_digest(unique),
            }
        )
    return excluded, sorted(evidence, key=lambda item: item["path"])


def _allocate_seed(
    base_seed: int, identity: str, occupied: set[int]
) -> tuple[int, int]:
    for nonce in range(10_000):
        label = identity if nonce == 0 else f"{identity}:nonce={nonce}"
        seed = _stable_seed(int(base_seed), label)
        if seed not in occupied:
            occupied.add(seed)
            return seed, nonce
    raise RuntimeError("could not allocate a disjoint fallback seed")


def build_atomic_manifest(
    design: FrequencyConditionedDesign,
    *,
    repo_root: str | Path,
    slim_path: str | Path,
    excluded_seed_artifacts: Sequence[str | Path],
    failed_fixed_cessation_artifacts: Sequence[str | Path],
) -> dict[str, Any]:
    """Allocate all screen, confirmation, and sensitivity seeds at once."""

    design.validate()
    root = Path(repo_root).resolve()
    excluded, excluded_evidence = _excluded_seed_evidence(excluded_seed_artifacts, root)
    failed_evidence = _artifact_evidence(
        failed_fixed_cessation_artifacts,
        root,
        label="failed fixed-cessation evidence",
    )
    occupied = set(excluded)
    rows = []
    phase_seeds: dict[str, dict[str, list[int]]] = {
        phase: {"simulation": [], "panel": []} for phase in PHASE_DRAW_COUNTS
    }
    for phase, draws in PHASE_DRAW_COUNTS.items():
        focal_only = phase in {"screen", "confirmation"}
        length = FOCAL_ONLY_LENGTH_BP if focal_only else PRODUCTION_LENGTH_BP
        position = FOCAL_ONLY_POSITION_BP if focal_only else PRODUCTION_POSITION_BP
        mode = FOCAL_ONLY_MODE if focal_only else PRODUCTION_MODE
        for coefficient in design.selection_coefficients:
            for target in design.target_frequencies:
                lower = float(target - design.af_half_width)
                upper = float(target + design.af_half_width)
                for draw_index in range(draws):
                    identity = (
                        f"han_frequency_conditioned__{phase}__s{coefficient:.6f}__"
                        f"af{target:.6f}__q{design.slim_scaling_factor:g}__"
                        f"draw{draw_index:03d}"
                    )
                    seed, nonce = _allocate_seed(
                        int(design.base_seed), identity, occupied
                    )
                    panel_seed, panel_nonce = _allocate_seed(
                        int(design.base_seed),
                        f"{identity}:uniform-exact-k-panel",
                        occupied,
                    )
                    phase_seeds[phase]["simulation"].append(seed)
                    phase_seeds[phase]["panel"].append(panel_seed)
                    rows.append(
                        {
                            "draw_id": identity,
                            "phase": phase,
                            "selection_coefficient": float(coefficient),
                            "target_allele_frequency": float(target),
                            "population_af_lower": lower,
                            "population_af_upper": upper,
                            "exact_sample_alt_count": round(
                                2 * int(design.sample_diploids) * float(target)
                            ),
                            "sample_diploids": int(design.sample_diploids),
                            "candidate_pool_diploids": int(
                                design.candidate_pool_diploids
                            ),
                            "slim_scaling_factor": float(design.slim_scaling_factor),
                            "mutation_age_generations": (
                                INTROGRESSION_MUTATION_AGE_GENERATIONS
                            ),
                            "selection_start_generations_ago": (
                                INTROGRESSION_PULSE_GENERATIONS
                            ),
                            "selection_end_generations_ago": 0.0,
                            "draw_index": draw_index,
                            "seed": seed,
                            "seed_nonce": nonce,
                            "panel_seed": panel_seed,
                            "panel_seed_nonce": panel_nonce,
                            "timeout_seconds": float(design.draw_timeout_seconds),
                            "sequence_length_bp": length,
                            "focal_position_bp": position,
                            "processing_mode": mode,
                        }
                    )
    if len(rows) != 840:
        raise RuntimeError("fallback manifest must contain exactly 840 draws")
    frame = pd.DataFrame(rows, columns=ROW_COLUMNS)
    fallback_rng_seeds = pd.concat(
        [frame["seed"], frame["panel_seed"]], ignore_index=True
    )
    if (
        frame["draw_id"].duplicated().any()
        or frame["seed"].duplicated().any()
        or frame["panel_seed"].duplicated().any()
        or fallback_rng_seeds.duplicated().any()
    ):
        raise RuntimeError("fallback draw identities and seeds must be unique")
    if set(fallback_rng_seeds).intersection(excluded):
        raise RuntimeError("fallback seeds collide with excluded seeds")
    seed_registry = {
        "schema": SEED_REGISTRY_SCHEMA,
        "excluded_union_count": len(excluded),
        "excluded_union_sha256": _seed_digest(excluded),
        "fallback": {
            phase: {
                component: {
                    "count": len(values),
                    "seed_sha256": _seed_digest(values),
                }
                for component, values in components.items()
            }
            for phase, components in phase_seeds.items()
        },
        "fallback_union_count": len(fallback_rng_seeds),
        "fallback_union_sha256": _seed_digest(fallback_rng_seeds),
        "simulation_seed_count": len(frame),
        "simulation_seed_sha256": _seed_digest(frame["seed"]),
        "panel_seed_count": len(frame),
        "panel_seed_sha256": _seed_digest(frame["panel_seed"]),
        "pairwise_disjoint": True,
    }
    manifest = {
        "schema": MANIFEST_SCHEMA,
        "status": "planned",
        "model": {
            "mode": MODEL_MODE,
            "archaic_specific_no_ils_by_construction": True,
            "mutation_population": "Neanderthal",
            "mutation_age_generations": INTROGRESSION_MUTATION_AGE_GENERATIONS,
            "introgression_pulse_generations_ago": INTROGRESSION_PULSE_GENERATIONS,
            "selection_start_generations_ago": INTROGRESSION_PULSE_GENERATIONS,
            "selection_end_generations_ago": 0.0,
            "selection_piecewise_populations": ["Loschbour", "Han"],
            "post_pulse_recipient_nonloss_conditioning": True,
            "terminal_han_frequency_conditioning_in_slim": True,
            "external_candidate_pool_frequency_gate": True,
            "uniform_exact_k_panel": True,
            "gamma_smc_statistics_used": False,
            "failed_fixed_cessation_observations_used": False,
        },
        "design": {
            "selection_coefficients": list(design.selection_coefficients),
            "target_allele_frequencies": list(design.target_frequencies),
            "af_half_width": design.af_half_width,
            "slim_scaling_factor": design.slim_scaling_factor,
            "candidate_pool_diploids": design.candidate_pool_diploids,
            "sample_diploids": design.sample_diploids,
            "minimum_hit_rate": design.minimum_hit_rate,
            "draw_timeout_seconds": design.draw_timeout_seconds,
            "phase_draws_per_cell": dict(PHASE_DRAW_COUNTS),
        },
        "scientific_caveats": [
            (
                "Internal terminal population-AF rejection and the external "
                "500-diploid AF gate are two distinct conditioning layers."
            ),
            (
                "Selection remains active through the present; this model does "
                "not estimate a selection-cessation time."
            ),
            (
                "Recipient nonloss and target-AF attainment condition on successful "
                "introgression and cannot estimate unconditional sweep prevalence."
            ),
            "SLiM Q=5 is an approximation and is frozen for every phase.",
            "No Gamma-SMC or TMRCA statistic enters planning, gating, or freezing.",
        ],
        "implementation_sources": _implementation_sources(root),
        "runtime": _runtime_record(root, slim_path),
        "excluded_seed_artifacts": excluded_evidence,
        "failed_fixed_cessation_evidence": failed_evidence,
        "seed_registry": seed_registry,
        "rows_sha256": _canonical_sha256(rows),
        "rows": rows,
    }
    return manifest


def _snapshot_sources(
    sources: Mapping[str, str], repo_root: Path, snapshot_dir: Path
) -> dict[str, dict[str, Any]]:
    snapshot_dir.mkdir(parents=True, exist_ok=True)
    records = {}
    for relative, expected_hash in sorted(sources.items()):
        source = (repo_root / relative).resolve()
        digest = _require_sha256(expected_hash, label=f"source {relative}")
        if sha256_file(source) != digest:
            raise ValueError(
                f"implementation source changed before snapshot: {relative}"
            )
        safe_name = relative.replace("/", "__").replace("\\", "__")
        destination = snapshot_dir / f"{safe_name}.sha256_{digest}"
        if destination.is_file():
            if sha256_file(destination) != digest:
                raise ValueError(f"source snapshot differs: {destination}")
        else:
            temporary = destination.with_name(f"{destination.name}.tmp.{os.getpid()}")
            try:
                shutil.copyfile(source, temporary)
                if sha256_file(temporary) != digest:
                    raise ValueError("temporary source snapshot differs")
                os.replace(temporary, destination)
            finally:
                temporary.unlink(missing_ok=True)
        try:
            destination.chmod(0o444)
        except OSError:
            pass
        records[relative] = {
            "path": destination.name,
            "sha256": digest,
            "bytes": destination.stat().st_size,
        }
    return records


def write_atomic_manifest(
    manifest: Mapping[str, Any],
    path: str | Path,
    *,
    repo_root: str | Path,
) -> dict[str, Any]:
    """Snapshot sources and write the one authoritative all-phase manifest."""

    root = Path(repo_root).resolve()
    destination = _require_repo_path(path, root, label="manifest")
    payload = json.loads(json.dumps(manifest, allow_nan=False))
    if payload.get("schema") != MANIFEST_SCHEMA:
        raise ValueError("fallback manifest schema is incompatible")
    reservation_path = _require_repo_path(
        destination.parent.parent / "han_frequency_conditioned_all_phase_plan.tsv",
        root,
        label="all-phase seed-reservation projection",
    )
    reservation_frame = pd.DataFrame(payload["rows"]).loc[:, list(ROW_COLUMNS)]
    if reservation_path.is_file():
        current_reservation = pd.read_csv(reservation_path, sep="\t")
        if set(current_reservation.columns) != set(ROW_COLUMNS):
            raise ValueError("existing fallback reservation projection differs")
        current_reservation = current_reservation.loc[:, list(ROW_COLUMNS)]
        pd.testing.assert_frame_equal(
            reservation_frame.reset_index(drop=True),
            current_reservation.reset_index(drop=True),
            check_dtype=False,
            check_exact=False,
            rtol=0.0,
            atol=1e-12,
        )
    if (
        reservation_frame[["seed", "panel_seed"]].isna().any().any()
        or not reservation_frame["seed"]
        .map(lambda value: int(value) == float(value))
        .all()
        or not reservation_frame["panel_seed"]
        .map(lambda value: int(value) == float(value))
        .all()
    ):
        raise ValueError("fallback reservation RNG columns are not non-null integers")
    staged_reservation: Path | None = None
    reservation_evidence_path = reservation_path
    if not reservation_path.is_file():
        staged_reservation = destination.parent / (
            f".han_frequency_conditioned_all_phase_plan.staged.{os.getpid()}.tsv"
        )
        _atomic_frame(staged_reservation, reservation_frame)
        reservation_evidence_path = staged_reservation
    payload["all_phase_seed_reservation"] = {
        "path": reservation_path.relative_to(root).as_posix(),
        "sha256": sha256_file(reservation_evidence_path),
        "rows": len(reservation_frame),
        "simulation_seed_column": "seed",
        "simulation_seed_sha256": _seed_digest(reservation_frame["seed"]),
        "panel_seed_column": "panel_seed",
        "panel_seed_sha256": _seed_digest(reservation_frame["panel_seed"]),
        "published_after_authoritative_manifest_commit": True,
        "authoritative_only_when_manifest_exists_and_binds_this_sha256": True,
    }
    snapshots = _snapshot_sources(
        payload["implementation_sources"],
        root,
        destination.parent / "source_snapshots",
    )
    payload["source_snapshots"] = snapshots
    payload["source_snapshots_sha256"] = _canonical_sha256(snapshots)
    try:
        if destination.exists():
            current = json.loads(destination.read_text(encoding="utf-8"))
            if current != payload:
                raise ValueError(
                    "existing fallback manifest differs; refusing overwrite"
                )
        else:
            _atomic_json(destination, payload)
        if staged_reservation is not None:
            os.replace(staged_reservation, reservation_path)
            staged_reservation = None
        return payload
    finally:
        if staged_reservation is not None:
            staged_reservation.unlink(missing_ok=True)


def _validate_manifest_rows(manifest: Mapping[str, Any]) -> pd.DataFrame:
    rows = manifest.get("rows")
    if not isinstance(rows, list) or manifest.get("rows_sha256") != _canonical_sha256(
        rows
    ):
        raise ValueError("fallback manifest row self-check failed")
    frame = pd.DataFrame(rows)
    if set(frame.columns) != set(ROW_COLUMNS) or len(frame.columns) != len(ROW_COLUMNS):
        raise ValueError("fallback manifest row columns are incompatible")
    frame = frame.loc[:, list(ROW_COLUMNS)]
    if len(frame) != 840 or frame["draw_id"].duplicated().any():
        raise ValueError("fallback manifest draw identities are incomplete")
    fallback_rng_seeds = pd.concat(
        [frame["seed"], frame["panel_seed"]], ignore_index=True
    )
    if fallback_rng_seeds.duplicated().any():
        raise ValueError("fallback manifest simulation/panel seeds are not unique")
    for phase, count in PHASE_DRAW_COUNTS.items():
        subset = frame[frame["phase"].astype(str) == phase]
        expected = count * len(SELECTION_COEFFICIENTS) * len(TARGET_FREQUENCIES)
        if len(subset) != expected:
            raise ValueError(f"fallback {phase} draw count is incompatible")
        per_cell = subset.groupby(
            ["selection_coefficient", "target_allele_frequency"]
        ).size()
        if len(per_cell) != 6 or set(per_cell.astype(int)) != {count}:
            raise ValueError(f"fallback {phase} per-cell plan is incompatible")
    if set(frame["phase"].astype(str)) != set(PHASE_DRAW_COUNTS):
        raise ValueError("fallback manifest contains an unknown phase")
    if set(frame["selection_coefficient"].astype(float)) != set(
        SELECTION_COEFFICIENTS
    ) or set(frame["target_allele_frequency"].astype(float)) != set(TARGET_FREQUENCIES):
        raise ValueError("fallback manifest does not contain the frozen six-cell grid")
    if set(frame["slim_scaling_factor"].astype(float)) != {SLIM_SCALING_FACTOR}:
        raise ValueError("fallback manifest Q differs from 5")
    if set(frame["candidate_pool_diploids"].astype(int)) != {CANDIDATE_POOL_DIPLOIDS}:
        raise ValueError("fallback manifest pool size differs")
    if set(frame["selection_end_generations_ago"].astype(float)) != {0.0}:
        raise ValueError("fallback selection must continue through the present")
    if set(frame["mutation_age_generations"].astype(float)) != {
        INTROGRESSION_MUTATION_AGE_GENERATIONS
    }:
        raise ValueError("fallback mutation age differs from 2400 generations")
    if set(frame["selection_start_generations_ago"].astype(float)) != {
        INTROGRESSION_PULSE_GENERATIONS
    }:
        raise ValueError("fallback selection onset differs from the pulse")
    if set(frame["sample_diploids"].astype(int)) != {SAMPLE_DIPLOIDS}:
        raise ValueError("fallback sample panel differs from N=100")
    if set(frame["timeout_seconds"].astype(float)) != {DEFAULT_DRAW_TIMEOUT_SECONDS}:
        raise ValueError("fallback per-draw timeout differs from 1800 seconds")
    for row in frame.to_dict(orient="records"):
        target = float(row["target_allele_frequency"])
        if not math.isclose(
            float(row["population_af_lower"]), target - AF_HALF_WIDTH, abs_tol=1e-12
        ) or not math.isclose(
            float(row["population_af_upper"]), target + AF_HALF_WIDTH, abs_tol=1e-12
        ):
            raise ValueError("fallback AF bounds differ from target +/- 0.025")
        if int(row["exact_sample_alt_count"]) != round(2 * SAMPLE_DIPLOIDS * target):
            raise ValueError("fallback exact-k count differs from target")
        phase = str(row["phase"])
        expected_length = (
            FOCAL_ONLY_LENGTH_BP
            if phase in {"screen", "confirmation"}
            else PRODUCTION_LENGTH_BP
        )
        expected_position = (
            FOCAL_ONLY_POSITION_BP
            if phase in {"screen", "confirmation"}
            else PRODUCTION_POSITION_BP
        )
        expected_mode = (
            FOCAL_ONLY_MODE if phase in {"screen", "confirmation"} else PRODUCTION_MODE
        )
        if (
            int(row["sequence_length_bp"]) != expected_length
            or int(row["focal_position_bp"]) != expected_position
            or str(row["processing_mode"]) != expected_mode
        ):
            raise ValueError("fallback phase geometry/processing differs")
        identity = (
            f"han_frequency_conditioned__{phase}__"
            f"s{float(row['selection_coefficient']):.6f}__"
            f"af{target:.6f}__q{SLIM_SCALING_FACTOR:g}__"
            f"draw{int(row['draw_index']):03d}"
        )
        if str(row["draw_id"]) != identity:
            raise ValueError("fallback draw identity differs from its fields")
        if (
            int(row["seed"]) < 1
            or int(row["panel_seed"]) < 1
            or int(row["seed_nonce"]) < 0
            or int(row["panel_seed_nonce"]) < 0
        ):
            raise ValueError("fallback RNG seed fields are invalid")
    for (phase, coefficient, target), group in frame.groupby(
        ["phase", "selection_coefficient", "target_allele_frequency"]
    ):
        if sorted(group["draw_index"].astype(int)) != list(
            range(PHASE_DRAW_COUNTS[str(phase)])
        ):
            raise ValueError(
                f"fallback draw indices are incomplete for {coefficient}/{target}"
            )
    screen_modes = set(
        frame.loc[
            frame["phase"].isin(["screen", "confirmation"]), "processing_mode"
        ].astype(str)
    )
    sensitivity_modes = set(
        frame.loc[frame["phase"].eq("sensitivity_10mb"), "processing_mode"].astype(str)
    )
    if screen_modes != {FOCAL_ONLY_MODE} or sensitivity_modes != {PRODUCTION_MODE}:
        raise ValueError("fallback processing modes are incompatible")
    registry = manifest.get("seed_registry")
    if (
        not isinstance(registry, Mapping)
        or registry.get("schema") != SEED_REGISTRY_SCHEMA
    ):
        raise ValueError("fallback seed registry is absent")
    if registry.get("fallback_union_sha256") != _seed_digest(fallback_rng_seeds):
        raise ValueError("fallback seed registry differs from rows")
    if int(registry.get("fallback_union_count", -1)) != len(fallback_rng_seeds):
        raise ValueError("fallback seed registry count differs")
    if registry.get("simulation_seed_sha256") != _seed_digest(frame["seed"]):
        raise ValueError("fallback simulation-seed registry differs")
    if registry.get("panel_seed_sha256") != _seed_digest(frame["panel_seed"]):
        raise ValueError("fallback panel-seed registry differs")
    fallback_registry = registry.get("fallback")
    if not isinstance(fallback_registry, Mapping):
        raise ValueError("fallback per-phase seed registry is absent")
    for phase in PHASE_DRAW_COUNTS:
        phase_rows = frame[frame["phase"].astype(str).eq(phase)]
        components = fallback_registry.get(phase)
        if not isinstance(components, Mapping):
            raise ValueError(f"fallback {phase} seed registry is absent")
        for component, column in (("simulation", "seed"), ("panel", "panel_seed")):
            values = components.get(component)
            if (
                not isinstance(values, Mapping)
                or int(values.get("count", -1)) != len(phase_rows)
                or values.get("seed_sha256") != _seed_digest(phase_rows[column])
            ):
                raise ValueError(f"fallback {phase} {component} seed registry differs")
    design = manifest.get("design")
    model = manifest.get("model")
    expected_design = {
        "selection_coefficients": list(SELECTION_COEFFICIENTS),
        "target_allele_frequencies": list(TARGET_FREQUENCIES),
        "af_half_width": AF_HALF_WIDTH,
        "slim_scaling_factor": SLIM_SCALING_FACTOR,
        "candidate_pool_diploids": CANDIDATE_POOL_DIPLOIDS,
        "sample_diploids": SAMPLE_DIPLOIDS,
        "minimum_hit_rate": MINIMUM_HIT_RATE,
        "draw_timeout_seconds": DEFAULT_DRAW_TIMEOUT_SECONDS,
        "phase_draws_per_cell": dict(PHASE_DRAW_COUNTS),
    }
    if design != expected_design:
        raise ValueError("fallback manifest design record differs")
    required_model = {
        "mode": MODEL_MODE,
        "archaic_specific_no_ils_by_construction": True,
        "mutation_population": "Neanderthal",
        "mutation_age_generations": INTROGRESSION_MUTATION_AGE_GENERATIONS,
        "introgression_pulse_generations_ago": INTROGRESSION_PULSE_GENERATIONS,
        "selection_start_generations_ago": INTROGRESSION_PULSE_GENERATIONS,
        "selection_end_generations_ago": 0.0,
        "selection_piecewise_populations": ["Loschbour", "Han"],
        "post_pulse_recipient_nonloss_conditioning": True,
        "terminal_han_frequency_conditioning_in_slim": True,
        "external_candidate_pool_frequency_gate": True,
        "uniform_exact_k_panel": True,
        "gamma_smc_statistics_used": False,
        "failed_fixed_cessation_observations_used": False,
    }
    if model != required_model:
        raise ValueError("fallback manifest model record differs")
    return frame


def load_manifest(
    path: str | Path,
    *,
    repo_root: str | Path,
    require_current_sources: bool = True,
) -> tuple[dict[str, Any], pd.DataFrame, str]:
    root = Path(repo_root).resolve()
    manifest_path = _require_repo_path(path, root, label="manifest")
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError("fallback manifest is unreadable") from error
    if manifest.get("schema") != MANIFEST_SCHEMA or manifest.get("status") != "planned":
        raise ValueError("fallback manifest schema/status is incompatible")
    frame = _validate_manifest_rows(manifest)
    reservation = manifest.get("all_phase_seed_reservation")
    if not isinstance(reservation, Mapping):
        raise ValueError("fallback all-phase seed reservation is absent")
    reservation_path = _require_repo_path(
        root / str(reservation.get("path", "")),
        root,
        label="all-phase seed-reservation projection",
    )
    if (
        not reservation_path.is_file()
        or sha256_file(reservation_path) != reservation.get("sha256")
        or int(reservation.get("rows", -1)) != len(frame)
        or reservation.get("simulation_seed_column") != "seed"
        or reservation.get("panel_seed_column") != "panel_seed"
        or reservation.get("simulation_seed_sha256") != _seed_digest(frame["seed"])
        or reservation.get("panel_seed_sha256") != _seed_digest(frame["panel_seed"])
    ):
        raise ValueError("fallback all-phase seed reservation differs")
    reservation_frame = pd.read_csv(reservation_path, sep="\t")
    if set(reservation_frame.columns) != set(ROW_COLUMNS):
        raise ValueError("fallback all-phase reservation columns differ")
    reservation_frame = reservation_frame.loc[:, list(ROW_COLUMNS)]
    if (
        reservation_frame[["seed", "panel_seed"]].isna().any().any()
        or not reservation_frame["seed"]
        .map(lambda value: int(value) == float(value))
        .all()
        or not reservation_frame["panel_seed"]
        .map(lambda value: int(value) == float(value))
        .all()
    ):
        raise ValueError("fallback reservation RNG columns are not non-null integers")
    pd.testing.assert_frame_equal(
        frame.reset_index(drop=True),
        reservation_frame.reset_index(drop=True),
        check_dtype=False,
        check_exact=False,
        rtol=0.0,
        atol=1e-12,
    )
    sources = manifest.get("implementation_sources")
    snapshots = manifest.get("source_snapshots")
    if not isinstance(sources, Mapping) or not isinstance(snapshots, Mapping):
        raise ValueError("fallback source provenance is incomplete")
    if manifest.get("source_snapshots_sha256") != _canonical_sha256(snapshots):
        raise ValueError("fallback source snapshot self-check failed")
    if set(sources) != set(snapshots):
        raise ValueError("fallback source snapshots do not cover implementation")
    for relative, expected_hash in sources.items():
        digest = _require_sha256(expected_hash, label=f"source {relative}")
        snapshot = (
            manifest_path.parent / "source_snapshots" / str(snapshots[relative]["path"])
        )
        if not snapshot.is_file() or sha256_file(snapshot) != digest:
            raise ValueError(f"fallback source snapshot checksum failed: {relative}")
        if require_current_sources:
            current = root / relative
            if not current.is_file() or sha256_file(current) != digest:
                raise ValueError(f"fallback implementation source drifted: {relative}")
    excluded_evidence = manifest.get("excluded_seed_artifacts", [])
    if not isinstance(excluded_evidence, list) or not excluded_evidence:
        raise ValueError("fallback excluded-seed evidence is absent")
    evidence_paths = []
    for evidence in excluded_evidence:
        path = root / str(evidence["path"])
        if not path.is_file() or sha256_file(path) != evidence["sha256"]:
            raise ValueError(
                f"excluded seed evidence changed after planning: {evidence['path']}"
            )
        evidence_paths.append(path)
    excluded, recomputed_evidence = _excluded_seed_evidence(evidence_paths, root)
    if recomputed_evidence != excluded_evidence:
        raise ValueError("fallback excluded RNG-seed expansion drifted")
    if manifest["seed_registry"].get("excluded_union_count") != len(excluded):
        raise ValueError("fallback excluded seed-registry count differs")
    if manifest["seed_registry"].get("excluded_union_sha256") != _seed_digest(excluded):
        raise ValueError("fallback excluded seed-registry digest differs")
    fallback_rng = set(frame["seed"].astype(int)) | set(frame["panel_seed"].astype(int))
    if fallback_rng.intersection(excluded):
        raise ValueError("fallback RNG seeds collide with reserved seeds")
    runtime = manifest.get("runtime")
    if not isinstance(runtime, Mapping) or not isinstance(runtime.get("slim"), Mapping):
        raise ValueError("fallback runtime binding is absent")
    slim_relative = str(runtime["slim"].get("path", ""))
    slim = _require_repo_path(root / slim_relative, root, label="SLiM executable")
    if require_current_sources and _runtime_record(root, slim) != runtime:
        raise ValueError("fallback runtime binding drifted")
    failed_evidence = manifest.get("failed_fixed_cessation_evidence")
    if not isinstance(failed_evidence, list) or not failed_evidence:
        raise ValueError("fallback failed fixed-cessation evidence is absent")
    failed_paths = [root / str(record["path"]) for record in failed_evidence]
    if (
        _artifact_evidence(failed_paths, root, label="failed fixed-cessation evidence")
        != failed_evidence
    ):
        raise ValueError("fallback failed fixed-cessation evidence drifted")
    return manifest, frame, sha256_file(manifest_path)


def _unit_for_row(row: Mapping[str, Any]) -> dict[str, Any]:
    target = float(row["target_allele_frequency"])
    return {
        "unit_id": str(row["draw_id"]),
        "demography_id": "ancient_eurasia_han_introgression",
        "demography_kind": "introgression",
        "population": INTROGRESSION_TARGET_POPULATION,
        "source_model": "stdpopsim AncientEurasia_9K19",
        "selection_origin": "archaic_specific_introgressed_standing_variation",
        "simulation_class": "selected",
        "selection_coefficient": float(row["selection_coefficient"]),
        "target_allele_frequency": target,
        "population_af_lower": float(row["population_af_lower"]),
        "population_af_upper": float(row["population_af_upper"]),
        "exact_sample_alt_count": int(row["exact_sample_alt_count"]),
        "sample_diploids": int(row["sample_diploids"]),
        "replicate_index": int(row["draw_index"]) + 1,
        "seed": int(row["seed"]),
        # Shared production validators require the canonical 10-Mb/5-Mb
        # semantic unit. Focal-only calibration changes only the sweep contig
        # after the raw production-shaped model has been constructed.
        "sequence_length_bp": PRODUCTION_LENGTH_BP,
        "focal_position_bp": PRODUCTION_POSITION_BP,
    }


def _audit_frequency_conditioned_events(
    sweep: Any, row: Mapping[str, Any]
) -> dict[str, Any]:
    events = list(sweep.extended_events)
    focal = [
        event
        for event in events
        if str(getattr(event, "single_site_id", "")) == FOCAL_SITE_ID
    ]
    draws = [event for event in focal if isinstance(event, stdpopsim.DrawMutation)]
    if len(focal) != 9 or len(events) != 9 or len(draws) != 1:
        raise ValueError("fallback sweep must contain one focal DrawMutation")
    draw = draws[0]
    if str(draw.population) != "Neanderthal" or not math.isclose(
        float(draw.time),
        INTROGRESSION_MUTATION_AGE_GENERATIONS,
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        raise ValueError("fallback mutation origin differs from Neanderthal age 2400")
    fitness = [
        event for event in focal if isinstance(event, stdpopsim.ChangeMutationFitness)
    ]
    if len(fitness) != 2:
        raise ValueError("fallback sweep must contain two fitness intervals")
    by_population = {str(event.population): event for event in fitness}
    if set(by_population) != {"Loschbour", INTROGRESSION_TARGET_POPULATION}:
        raise ValueError("fallback fitness populations are incompatible")
    q = float(row["slim_scaling_factor"])
    realized_split = math.floor(HAN_SPLIT_GENERATIONS / q) * q
    coefficient = float(row["selection_coefficient"])
    expected_intervals = {
        "Loschbour": (INTROGRESSION_PULSE_GENERATIONS, realized_split),
        INTROGRESSION_TARGET_POPULATION: (realized_split, 0.0),
    }
    for population, event in by_population.items():
        expected_start, expected_end = expected_intervals[population]
        if not math.isclose(
            float(event.start_time), expected_start, abs_tol=1e-12
        ) or not math.isclose(float(event.end_time), expected_end, abs_tol=1e-12):
            raise ValueError(
                "fallback selection interval differs from pulse-to-present"
            )
        if not math.isclose(float(event.selection_coeff), coefficient, abs_tol=1e-12):
            raise ValueError("fallback selection coefficient differs")
        if not math.isclose(
            float(event.dominance_coeff), DOMINANCE_COEFFICIENT, abs_tol=1e-12
        ):
            raise ValueError("fallback dominance coefficient differs")
    conditions = [
        event
        for event in focal
        if isinstance(event, stdpopsim.ConditionOnAlleleFrequency)
    ]
    terminal = [
        event
        for event in conditions
        if str(event.population) == INTROGRESSION_TARGET_POPULATION
        and math.isclose(float(event.start_time), 0.0, abs_tol=1e-12)
        and math.isclose(float(event.end_time), 0.0, abs_tol=1e-12)
    ]
    if len(terminal) != 2 or {str(event.op) for event in terminal} != {">=", "<="}:
        raise ValueError("fallback must retain exactly two terminal Han AF conditions")
    observed_bounds = {
        str(event.op): float(event.allele_frequency) for event in terminal
    }
    expected_bounds = {
        ">=": float(row["population_af_lower"]),
        "<=": float(row["population_af_upper"]),
    }
    if any(
        not math.isclose(observed_bounds[op], value, abs_tol=1e-12)
        for op, value in expected_bounds.items()
    ):
        raise ValueError("fallback terminal Han AF bounds differ from the plan")
    recipient_nonloss = [
        event
        for event in conditions
        if str(event.op) == ">"
        and math.isclose(float(event.allele_frequency), 0.0, abs_tol=1e-12)
        and str(event.population) in {"Loschbour", INTROGRESSION_TARGET_POPULATION}
    ]
    if len(recipient_nonloss) != 2 or {
        str(event.population) for event in recipient_nonloss
    } != {"Loschbour", INTROGRESSION_TARGET_POPULATION}:
        raise ValueError("fallback recipient nonloss conditioning is incomplete")
    source_nonloss = [
        event
        for event in conditions
        if str(event.population) == "Neanderthal"
        and str(event.op) == ">"
        and math.isclose(float(event.allele_frequency), 0.0, abs_tol=1e-12)
        and isinstance(event.start_time, stdpopsim.GenerationAfter)
        and math.isclose(
            float(event.start_time),
            INTROGRESSION_MUTATION_AGE_GENERATIONS,
            abs_tol=1e-12,
        )
        and math.isclose(
            float(event.end_time), PRE_PULSE_SOURCE_CHECK_GENERATIONS, abs_tol=1e-12
        )
    ]
    source_presence = [
        event
        for event in conditions
        if str(event.population) == "Neanderthal"
        and str(event.op) == ">="
        and math.isclose(
            float(event.allele_frequency),
            INTROGRESSION_SOURCE_PRESENCE_EPSILON,
            abs_tol=1e-15,
        )
        and math.isclose(
            float(event.start_time), PRE_PULSE_SOURCE_CHECK_GENERATIONS, abs_tol=1e-12
        )
        and math.isclose(
            float(event.end_time), PRE_PULSE_SOURCE_CHECK_GENERATIONS, abs_tol=1e-12
        )
    ]
    if len(source_nonloss) != 1 or len(source_presence) != 1:
        raise ValueError("fallback Neanderthal source conditioning differs")
    recipient_by_population = {
        str(event.population): event for event in recipient_nonloss
    }
    if not isinstance(
        recipient_by_population["Loschbour"].start_time,
        stdpopsim.GenerationAfter,
    ) or not isinstance(
        recipient_by_population[INTROGRESSION_TARGET_POPULATION].start_time,
        stdpopsim.GenerationAfter,
    ):
        raise ValueError("fallback recipient nonloss generation semantics differ")
    if not math.isclose(
        float(recipient_by_population["Loschbour"].start_time),
        INTROGRESSION_PULSE_GENERATIONS,
        abs_tol=1e-12,
    ) or not math.isclose(
        float(recipient_by_population[INTROGRESSION_TARGET_POPULATION].start_time),
        realized_split,
        abs_tol=1e-12,
    ):
        raise ValueError("fallback recipient nonloss intervals differ")
    if not math.isclose(
        float(recipient_by_population["Loschbour"].end_time),
        realized_split,
        abs_tol=1e-12,
    ) or not math.isclose(
        float(recipient_by_population[INTROGRESSION_TARGET_POPULATION].end_time),
        0.0,
        abs_tol=1e-12,
    ):
        raise ValueError("fallback recipient nonloss endpoints differ")
    record = sweep.provenance_record()
    if "cessation" in record.get("selection", {}):
        raise ValueError("fallback sweep must not contain fixed-cessation provenance")
    if record.get("focused_sample_endpoint_ascertainment") is not None:
        raise ValueError("fallback sweep must not remove internal terminal conditions")
    single_site_dfes = [dfe for dfe in sweep.contig.dfe_list if dfe.id == FOCAL_SITE_ID]
    if len(single_site_dfes) != 1 or len(single_site_dfes[0].mutation_types) != 1:
        raise ValueError("fallback focal mutation-type definition differs")
    mutation_type = single_site_dfes[0].mutation_types[0]
    if mutation_type.distribution_type != "f" or not math.isclose(
        float(mutation_type.distribution_args[0]), 0.0, abs_tol=1e-12
    ):
        raise ValueError("fallback focal type must be neutral outside callbacks")
    realized = {
        "mutation_age_generations": round(float(draw.time) / q) * q,
        "source_nonloss_generation_after_start": (
            round(float(source_nonloss[0].start_time) / q) * q - q
        ),
        "source_check_generations": round(PRE_PULSE_SOURCE_CHECK_GENERATIONS / q) * q,
        "introgression_and_selection_onset_generations": round(
            INTROGRESSION_PULSE_GENERATIONS / q
        )
        * q,
        "loschbour_nonloss_generation_after_start": (
            round(float(recipient_by_population["Loschbour"].start_time) / q) * q - q
        ),
        "demographic_han_split_generations": math.floor(HAN_SPLIT_GENERATIONS / q) * q,
        "selection_population_transition_generations": round(realized_split / q) * q,
        "han_nonloss_generation_after_start": (
            round(
                float(
                    recipient_by_population[INTROGRESSION_TARGET_POPULATION].start_time
                )
                / q
            )
            * q
            - q
        ),
    }
    if realized != {
        "mutation_age_generations": 2400.0,
        "source_nonloss_generation_after_start": 2395.0,
        "source_check_generations": 2275.0,
        "introgression_and_selection_onset_generations": 2270.0,
        "loschbour_nonloss_generation_after_start": 2265.0,
        "demographic_han_split_generations": 2015.0,
        "selection_population_transition_generations": 2015.0,
        "han_nonloss_generation_after_start": 2010.0,
    }:
        raise ValueError("fallback Q5 realized event schedule differs")
    serialized = serialize_extended_events(events)
    return {
        "schema": "gamma-smc.han-frequency-conditioned-events/v1",
        "mutation": {
            "population": str(draw.population),
            "age_generations": float(draw.time),
        },
        "realized_han_split_generations_ago": realized_split,
        "selection_intervals": [
            {
                "population": population,
                "start_time": expected_intervals[population][0],
                "end_time": expected_intervals[population][1],
                "selection_coefficient": coefficient,
            }
            for population in ("Loschbour", INTROGRESSION_TARGET_POPULATION)
        ],
        "recipient_nonloss_condition_count": len(recipient_nonloss),
        "source_condition_count": len(source_nonloss) + len(source_presence),
        "terminal_han_conditions": observed_bounds,
        "requested_event_times_generations": {
            "mutation_age": INTROGRESSION_MUTATION_AGE_GENERATIONS,
            "source_check": PRE_PULSE_SOURCE_CHECK_GENERATIONS,
            "introgression_and_selection_onset": INTROGRESSION_PULSE_GENERATIONS,
            "catalog_han_split": HAN_SPLIT_GENERATIONS,
            "selection_population_transition": realized_split,
        },
        "realized_q5_event_times_generations": realized,
        "realized_event_schedule_sha256": _canonical_sha256(realized),
        "focal_event_count": len(focal),
        "focal_mutation_type_distribution": "fixed_zero_callback_selected",
        "dominance_coefficient": DOMINANCE_COEFFICIENT,
        "selection_cessation_applied": False,
        "external_terminal_condition_removal_applied": False,
        "serialized_extended_events": serialized,
        "serialized_extended_events_sha256": _canonical_sha256(serialized),
    }


def build_frequency_conditioned_objects(
    row: Mapping[str, Any],
) -> tuple[stdpopsim.DemographicModel, Any, dict[str, int], dict[str, Any]]:
    """Build the original pulse-to-present internally conditioned Han model."""

    unit = _unit_for_row(row)
    q = float(row["slim_scaling_factor"])
    sweep = build_raw_han_selected_sweep(unit, slim_scaling_factor=q)
    event_record = _audit_frequency_conditioned_events(sweep, row)
    sweep = replace(
        sweep,
        contig=build_focal_contig(
            sequence_length=int(row["sequence_length_bp"]),
            focal_position=int(row["focal_position_bp"]),
        ),
    )
    return (
        sweep.demographic_model,
        sweep,
        {INTROGRESSION_TARGET_POPULATION: int(row["candidate_pool_diploids"])},
        {
            "model_mode": MODEL_MODE,
            "event_audit": event_record,
            "sweep": sweep.provenance_record(),
        },
    )


def _inspect_raw_type1(
    ts: tskit.TreeSequence,
    expectation: Mapping[str, Any],
) -> dict[str, Any]:
    coordinate = float(expectation["focal_position_0based"])
    sites = [
        site
        for site in ts.sites()
        if math.isclose(float(site.position), coordinate, abs_tol=1e-9)
    ]
    entries: list[tuple[tskit.Mutation, Mapping[str, Any]]] = []
    for site in sites:
        for mutation in site.mutations:
            metadata = mutation.metadata
            if not isinstance(metadata, Mapping):
                raise focused_simulation.SelectedFocalIdentityError(
                    "raw focal mutation metadata is malformed"
                )
            mutation_list = metadata.get("mutation_list")
            if not isinstance(mutation_list, list):
                raise focused_simulation.SelectedFocalIdentityError(
                    "raw focal mutation list is malformed"
                )
            for entry in mutation_list:
                if not isinstance(entry, Mapping):
                    raise focused_simulation.SelectedFocalIdentityError(
                        "raw focal mutation entry is malformed"
                    )
                entries.append((mutation, entry))
    selected = [
        (mutation, entry)
        for mutation, entry in entries
        if int(entry.get("mutation_type", -1))
        == int(expectation["expected_mutation_type"])
    ]
    if len(sites) != 1 or len(sites[0].mutations) != 1 or len(entries) != 1:
        raise focused_simulation.SelectedFocalIdentityError(
            "raw focal identity requires one nonrecurrent mutation and metadata entry"
        )
    if len(selected) != 1:
        raise focused_simulation.SelectedFocalIdentityError(
            "raw focal identity lacks exactly one type-1 mutation"
        )
    mutation, entry = selected[0]
    required = {"selection_coeff", "subpopulation", "slim_time", "nucleotide"}
    missing = sorted(required.difference(entry))
    if missing:
        raise focused_simulation.SelectedFocalIdentityError(
            "raw focal metadata is missing: " + ", ".join(missing)
        )
    if mutation.parent != tskit.NULL:
        raise focused_simulation.SelectedFocalIdentityError(
            "raw focal type-1 mutation is recurrent"
        )
    q = float(expectation["slim_scaling_factor"])
    observed_age = float(mutation.time) * q
    expected_age = float(expectation["realized_origin_time_generations"])
    if not math.isclose(observed_age, expected_age, abs_tol=1e-9):
        raise focused_simulation.SelectedFocalIdentityError(
            "raw focal mutation age differs"
        )
    if int(entry["subpopulation"]) != int(expectation["source_subpopulation"]):
        raise focused_simulation.SelectedFocalIdentityError(
            "raw focal source subpopulation differs"
        )
    if not math.isclose(
        float(entry["selection_coeff"]),
        float(expectation["expected_selection_coeff"]),
        abs_tol=1e-12,
    ):
        raise focused_simulation.SelectedFocalIdentityError(
            "raw focal mutation coefficient differs"
        )
    try:
        cycle = int(ts.metadata["SLiM"]["cycle"])
        observed_q = float(ts.metadata["SLiM"]["user_metadata"]["Q"][0])
    except (KeyError, TypeError, ValueError, IndexError) as error:
        raise focused_simulation.SelectedFocalIdentityError(
            "raw focal SLiM Q/cycle provenance is absent"
        ) from error
    if not math.isclose(observed_q, q, abs_tol=1e-12) or int(
        entry["slim_time"]
    ) != round(cycle - float(mutation.time)):
        raise focused_simulation.SelectedFocalIdentityError(
            "raw focal Q/time provenance differs"
        )
    return {
        "schema": "gamma-smc.han-frequency-conditioned-raw-type1/v1",
        "status": "valid",
        "focal_position_0based": coordinate,
        "mutation_node": int(mutation.node),
        "mutation_type": int(entry["mutation_type"]),
        "source_subpopulation": int(entry["subpopulation"]),
        "origin_time_raw_q_units": float(mutation.time),
        "origin_time_generations": observed_age,
        "slim_scaling_factor": observed_q,
        "slim_cycle": cycle,
        "slim_time": int(entry["slim_time"]),
    }


def _raw_type1_counts(
    ts: tskit.TreeSequence, identity: Mapping[str, Any]
) -> np.ndarray:
    samples = {int(node) for node in ts.samples()}
    tree = ts.at(float(identity["focal_position_0based"]))
    carriers = {int(node) for node in tree.samples(int(identity["mutation_node"]))}
    counts = []
    assigned = []
    for individual in ts.individuals():
        nodes = [int(node) for node in individual.nodes if int(node) in samples]
        if not nodes:
            continue
        if len(nodes) != 2:
            raise ValueError("fallback sampled individual is not diploid")
        assigned.extend(nodes)
        counts.append(sum(node in carriers for node in nodes))
    if set(assigned) != samples or len(assigned) != len(set(assigned)) or not counts:
        raise ValueError("fallback sample nodes are not exact diploid pairs")
    return np.asarray(counts, dtype=np.int8)


def _select_exact_sample_panel_at_position(
    ts: tskit.TreeSequence,
    *,
    focal_position: int,
    focal_genotype_counts: Sequence[int] | None,
    sample_diploids: int,
    exact_alt_count: int,
    seed: int,
) -> tuple[tskit.TreeSequence, dict[str, Any]]:
    """Calibration-local, position-aware copy of the production sampler.

    Genotype-count triples are weighted by the exact number of underlying
    subsets and individuals are sampled uniformly within each class. This is
    uniform over all eligible exact-k/QC panels while avoiding mutation of the
    production module's 5-Mb global focal coordinate.
    """

    pool = build_diploid_pair_table(
        ts, int(focal_position), focal_genotype_counts=focal_genotype_counts
    )
    pool_counts = pool["focal_selected_allele_count"].to_numpy(dtype=int)
    available = {
        value: int(np.count_nonzero(pool_counts == value)) for value in (0, 1, 2)
    }
    triples = feasible_exact_panel_genotype_triples(
        pool_counts,
        sample_diploids=int(sample_diploids),
        exact_alt_count=int(exact_alt_count),
    )
    if not triples:
        raise focused_simulation.SimulationExhausted(
            "candidate pool has no exact-AF subset satisfying genotype QC"
        )
    log_weights = np.asarray(
        [record["log_number_of_subsets"] for record in triples], dtype=float
    )
    weights = np.exp(log_weights - log_weights.max())
    rng = np.random.default_rng(int(seed))
    triple_index = int(rng.choice(len(triples), p=weights / weights.sum()))
    chosen = triples[triple_index]
    chosen_numbers = (
        int(chosen["hom_ref"]),
        int(chosen["heterozygous"]),
        int(chosen["hom_alt"]),
    )
    chosen_rows: list[int] = []
    for genotype_count, number in enumerate(chosen_numbers):
        choices = np.flatnonzero(pool_counts == genotype_count)
        chosen_rows.extend(
            int(value)
            for value in rng.choice(choices, size=number, replace=False).tolist()
        )
    chosen_rows.sort()
    nodes = (
        pool.iloc[chosen_rows][["sample_node_0", "sample_node_1"]]
        .to_numpy(dtype=np.int32)
        .reshape(-1)
    )
    if ts.num_migrations:
        tables = ts.dump_tables()
        node_flags = np.asarray(tables.nodes.flags, dtype=np.uint32).copy()
        node_flags &= np.bitwise_not(np.uint32(tskit.NODE_IS_SAMPLE))
        node_flags[nodes] |= np.uint32(tskit.NODE_IS_SAMPLE)
        table_nodes = tables.nodes
        table_nodes.set_columns(
            flags=node_flags,
            time=table_nodes.time,
            population=table_nodes.population,
            individual=table_nodes.individual,
            metadata=table_nodes.metadata,
            metadata_offset=table_nodes.metadata_offset,
        )
        panel = tables.tree_sequence()
    else:
        panel = ts.simplify(
            samples=nodes,
            filter_populations=False,
            filter_individuals=True,
            filter_sites=True,
            keep_unary=False,
        )
    return panel, {
        "schema": "gamma-smc.position-aware-uniform-exact-k-panel/v1",
        "candidate_pool_diploids": len(pool),
        "candidate_pool_genotype_counts": {
            "hom_ref": available[0],
            "heterozygous": available[1],
            "hom_alt": available[2],
        },
        "eligible_genotype_count_triples": len(triples),
        "chosen_sample_genotype_counts": {
            "hom_ref": chosen_numbers[0],
            "heterozygous": chosen_numbers[1],
            "hom_alt": chosen_numbers[2],
        },
        "selection_method": (
            "uniform_over_all_size_N_exact_alt_count_subsets_passing_genotype_qc"
        ),
        "panel_seed": int(seed),
        "focal_position_0based": int(focal_position),
        "post_simplification_manifest_rebuilt_by_caller": True,
    }


def _draw_contract(
    row: Mapping[str, Any],
    *,
    manifest_sha256: str,
    implementation_sources: Mapping[str, str],
    runtime: Mapping[str, Any],
) -> dict[str, Any]:
    normal = str(row["processing_mode"]) == PRODUCTION_MODE
    return {
        "schema": DRAW_SCHEMA,
        "manifest_sha256": _require_sha256(manifest_sha256, label="fallback manifest"),
        "manifest_row": _normalized_manifest_row(row),
        "model": {
            "mode": MODEL_MODE,
            "terminal_han_af_conditioning_in_slim": True,
            "post_pulse_recipient_nonloss_conditioning": True,
            "selection_end_generations_ago": 0.0,
            "failed_fixed_cessation_observations_used": False,
            "gamma_smc_statistics_used": False,
        },
        "processing": {
            "mode": str(row["processing_mode"]),
            "stdpopsim_recapitation_applied": normal,
            "post_slim_neutral_mutation_overlay_applied": normal,
            "strict_type1_identity": True,
            "external_candidate_pool_gate": True,
            "uniform_exact_k_panel_sampler": True,
            "panel_sampler_source": (
                "focused_selection_simulation.select_exact_sample_panel"
                if normal
                else "calibration_local_position_aware_equivalent"
            ),
        },
        "slim_burn_in": DEFAULT_SLIM_BURN_IN,
        "draw_timeout_seconds": float(row["timeout_seconds"]),
        "implementation_sources": dict(implementation_sources),
        "runtime": dict(runtime),
    }


def _select_validate_and_save_panel(
    ts: tskit.TreeSequence,
    row: Mapping[str, Any],
    expectation: Mapping[str, Any],
    directory: Path,
    *,
    raw_forward: bool,
) -> dict[str, Any]:
    coordinate = int(row["focal_position_bp"])
    panel_seed = int(row["panel_seed"])
    raw_counts = None
    if raw_forward:
        raw_identity = _inspect_raw_type1(ts, expectation)
        raw_counts = _raw_type1_counts(ts, raw_identity)
        panel_ts, panel_record = _select_exact_sample_panel_at_position(
            ts,
            focal_position=coordinate,
            focal_genotype_counts=raw_counts,
            sample_diploids=int(row["sample_diploids"]),
            exact_alt_count=int(row["exact_sample_alt_count"]),
            seed=panel_seed,
        )
    else:
        if coordinate != PRODUCTION_POSITION_BP:
            raise ValueError(
                "normal sensitivity focal position differs from production"
            )
        panel_ts, panel_record = select_exact_sample_panel(
            ts,
            sample_diploids=int(row["sample_diploids"]),
            exact_alt_count=int(row["exact_sample_alt_count"]),
            seed=panel_seed,
        )
    if raw_forward:
        identity = _inspect_raw_type1(panel_ts, expectation)
        panel_counts = _raw_type1_counts(panel_ts, identity)
        pair_table = build_diploid_pair_table(
            panel_ts, coordinate, focal_genotype_counts=panel_counts
        )
    else:
        identity = validate_selected_focal_identity(panel_ts, expectation)
        pair_table = build_diploid_pair_table(panel_ts, coordinate)
        panel_counts = pair_table["focal_selected_allele_count"].to_numpy(dtype=int)
    if (
        len(pair_table) != int(row["sample_diploids"])
        or int(panel_counts.sum()) != int(row["exact_sample_alt_count"])
        or np.count_nonzero(panel_counts == 0) < 2
        or np.count_nonzero(panel_counts == 1) < 1
        or np.count_nonzero(panel_counts == 2) < 2
    ):
        raise ValueError("fallback exact panel violates size, exact-k, or genotype QC")
    tree_path = directory / "accepted_panel.trees"
    manifest_path = directory / "accepted_panel_manifest.tsv"
    _atomic_tree(tree_path, panel_ts)
    _atomic_frame(manifest_path, pair_table)
    reloaded = tskit.load(tree_path)
    if raw_forward:
        reloaded_identity = _inspect_raw_type1(reloaded, expectation)
        reloaded_counts = _raw_type1_counts(reloaded, reloaded_identity)
        rebuilt = build_diploid_pair_table(
            reloaded, coordinate, focal_genotype_counts=reloaded_counts
        )
    else:
        reloaded_identity = validate_selected_focal_identity(reloaded, expectation)
        rebuilt = build_diploid_pair_table(reloaded, coordinate)
    persisted = pd.read_csv(manifest_path, sep="\t")
    pd.testing.assert_frame_equal(
        rebuilt.reset_index(drop=True),
        persisted.reset_index(drop=True),
        check_dtype=False,
        check_exact=False,
        rtol=0.0,
        atol=1e-12,
    )
    return {
        "panel_seed": panel_seed,
        "panel_sampler_source": (
            "calibration_local_position_aware_equivalent"
            if raw_forward
            else "focused_selection_simulation.select_exact_sample_panel"
        ),
        "panel_record_json": json.dumps(
            panel_record, sort_keys=True, separators=(",", ":")
        ),
        "panel_identity_json": json.dumps(
            reloaded_identity, sort_keys=True, separators=(",", ":")
        ),
        "panel_tree_path": tree_path.name,
        "panel_tree_sha256": sha256_file(tree_path),
        "panel_manifest_path": manifest_path.name,
        "panel_manifest_sha256": sha256_file(manifest_path),
        "panel_alt_count": int(panel_counts.sum()),
        "panel_af": float(panel_counts.mean() / 2),
        "panel_hom_ref": int(np.count_nonzero(panel_counts == 0)),
        "panel_heterozygous": int(np.count_nonzero(panel_counts == 1)),
        "panel_hom_alt": int(np.count_nonzero(panel_counts == 2)),
        "panel_roundtrip_validated": True,
    }


def _empty_panel_result() -> dict[str, Any]:
    return {
        "panel_seed": None,
        "panel_sampler_source": "",
        "panel_record_json": "",
        "panel_identity_json": "",
        "panel_tree_path": "",
        "panel_tree_sha256": "",
        "panel_manifest_path": "",
        "panel_manifest_sha256": "",
        "panel_alt_count": None,
        "panel_af": None,
        "panel_hom_ref": None,
        "panel_heterozygous": None,
        "panel_hom_alt": None,
        "panel_roundtrip_validated": False,
    }


def _evaluate_completed_tree(
    ts: tskit.TreeSequence,
    row: Mapping[str, Any],
    expectation: Mapping[str, Any],
    unit: Mapping[str, Any],
    directory: Path,
) -> dict[str, Any]:
    raw_forward = str(row["processing_mode"]) == FOCAL_ONLY_MODE
    coordinate = int(row["focal_position_bp"])
    for name in ("accepted_panel.trees", "accepted_panel_manifest.tsv"):
        (directory / name).unlink(missing_ok=True)
    try:
        if raw_forward:
            identity = _inspect_raw_type1(ts, expectation)
            counts = _raw_type1_counts(ts, identity)
            build_diploid_pair_table(ts, coordinate, focal_genotype_counts=counts)
        else:
            identity = focused_simulation.inspect_selected_focal_identity(
                ts, expectation
            )
            if identity["status"] != "valid":
                return {
                    "terminal_class": "identity_collision",
                    "strict_type1_identity_passed": False,
                    "identity_error": str(identity["status"]),
                    "identity_json": json.dumps(
                        identity, sort_keys=True, separators=(",", ":")
                    ),
                    "candidate_pool_alt_count": None,
                    "candidate_pool_af": None,
                    "candidate_pool_frequency_gate_passed": False,
                    "exact_panel_feasible": False,
                    "hit": False,
                    **_empty_panel_result(),
                }
            pool_table = build_diploid_pair_table(ts, coordinate)
            counts = pool_table["focal_selected_allele_count"].to_numpy(dtype=int)
    except focused_simulation.SelectedFocalIdentityError as error:
        if raw_forward:
            raise
        # Normal recap/background recurrence is a prespecified failed
        # sensitivity observation, not a seed retry.
        return {
            "terminal_class": "identity_collision",
            "strict_type1_identity_passed": False,
            "identity_error": f"{type(error).__name__}: {error}",
            "identity_json": "",
            "candidate_pool_alt_count": None,
            "candidate_pool_af": None,
            "candidate_pool_frequency_gate_passed": False,
            "exact_panel_feasible": False,
            "hit": False,
            **_empty_panel_result(),
        }
    if len(counts) != int(row["candidate_pool_diploids"]):
        raise ValueError("fallback candidate pool has the wrong diploid count")
    gate = candidate_pool_frequency_gate(unit, counts)
    feasible = (
        feasible_exact_panel_genotype_triples(
            counts,
            sample_diploids=int(row["sample_diploids"]),
            exact_alt_count=int(row["exact_sample_alt_count"]),
        )
        if gate["passed"]
        else []
    )
    alt_count = int(counts.sum())
    terminal_class = (
        "loss"
        if alt_count == 0
        else "fixed"
        if alt_count == 2 * len(counts)
        else "segregating"
    )
    panel_result = _empty_panel_result()
    if gate["passed"] and feasible:
        panel_result = _select_validate_and_save_panel(
            ts,
            row,
            expectation,
            directory,
            raw_forward=raw_forward,
        )
    else:
        (directory / "accepted_panel.trees").unlink(missing_ok=True)
        (directory / "accepted_panel_manifest.tsv").unlink(missing_ok=True)
    hit = bool(
        gate["passed"] and feasible and panel_result["panel_roundtrip_validated"]
    )
    return {
        "terminal_class": terminal_class,
        "strict_type1_identity_passed": True,
        "identity_error": "",
        "identity_json": json.dumps(identity, sort_keys=True, separators=(",", ":")),
        "candidate_pool_alt_count": alt_count,
        "candidate_pool_af": float(alt_count / (2 * len(counts))),
        "candidate_pool_frequency_gate_passed": bool(gate["passed"]),
        "candidate_pool_frequency_gate_json": json.dumps(
            gate, sort_keys=True, separators=(",", ":")
        ),
        "exact_panel_feasible": bool(feasible),
        "hit": hit,
        **panel_result,
    }


def _run_draw(task: Mapping[str, Any]) -> dict[str, Any]:
    row = _normalized_manifest_row(task["row"])
    work_root = Path(task["work_root"]).resolve()
    slim_path = Path(task["slim_path"]).resolve()
    timeout_seconds = float(task["timeout_seconds"])
    if not math.isclose(timeout_seconds, float(row["timeout_seconds"]), abs_tol=1e-12):
        raise ValueError("fallback task timeout differs from its manifest row")
    manifest_sha256 = str(task["manifest_sha256"])
    sources = dict(task["implementation_sources"])
    runtime = dict(task["runtime"])
    directory = work_root / str(row["draw_id"])
    directory.mkdir(parents=True, exist_ok=True)
    completion_path = directory / "completion.json"
    contract = _draw_contract(
        row,
        manifest_sha256=manifest_sha256,
        implementation_sources=sources,
        runtime=runtime,
    )
    with exclusive_unit_lock(
        directory / ".draw.lock",
        unit_id=str(row["draw_id"]),
        phase=f"han_frequency_conditioned_{row['phase']}",
    ):
        if completion_path.is_file():
            payload = json.loads(completion_path.read_text(encoding="utf-8"))
            if (
                payload.get("schema") != DRAW_SCHEMA
                or payload.get("contract_sha256") != _canonical_sha256(contract)
                or payload.get("contract") != contract
            ):
                raise ValueError("cached fallback draw is incompatible")
            result = dict(payload["result"])
            for path_key, hash_key in (
                ("panel_tree_path", "panel_tree_sha256"),
                ("panel_manifest_path", "panel_manifest_sha256"),
            ):
                if result.get(path_key):
                    path = directory / str(result[path_key])
                    if not path.is_file() or sha256_file(path) != result[hash_key]:
                        raise ValueError("cached fallback panel artifact differs")
            return result
        started = perf_counter()
        trajectory_path = directory / "trajectory.csv"
        trajectory_path.unlink(missing_ok=True)
        try:
            unit = _unit_for_row(row)
            model, sweep, samples, model_record = build_frequency_conditioned_objects(
                row
            )
            expectation = _selected_focal_expectation(
                model, sweep, float(row["slim_scaling_factor"])
            )
            expectation["focal_position_0based"] = int(row["focal_position_bp"])
            engine = stdpopsim.get_engine("slim")
            normal = str(row["processing_mode"]) == PRODUCTION_MODE
            with _wall_clock_timeout(timeout_seconds):
                if normal:
                    ts = engine.simulate(
                        model,
                        sweep.contig,
                        samples,
                        seed=int(row["seed"]),
                        extended_events=list(sweep.extended_events),
                        slim_path=str(slim_path),
                        slim_scaling_factor=float(row["slim_scaling_factor"]),
                        slim_burn_in=DEFAULT_SLIM_BURN_IN,
                        logfile=str(trajectory_path),
                        logfile_interval=100,
                        keep_mutation_ids_as_alleles=False,
                    )
                else:
                    raw = engine.simulate(
                        model,
                        sweep.contig,
                        samples,
                        seed=int(row["seed"]),
                        extended_events=list(sweep.extended_events),
                        slim_path=str(slim_path),
                        slim_scaling_factor=float(row["slim_scaling_factor"]),
                        slim_burn_in=DEFAULT_SLIM_BURN_IN,
                        logfile=str(trajectory_path),
                        logfile_interval=100,
                        keep_mutation_ids_as_alleles=True,
                        _recap_and_rescale=False,
                    )
                    ts = engine._simplify_remembered(raw)
            evaluation = _evaluate_completed_tree(ts, row, expectation, unit, directory)
            if not trajectory_path.is_file():
                raise ValueError("fallback SLiM trajectory log is absent")
            result = {
                **row,
                "status": "complete",
                "evaluable": True,
                "model_mode": MODEL_MODE,
                "internal_terminal_han_af_conditioning": True,
                "selection_through_present": True,
                "recipient_nonloss_conditioning": True,
                "stdpopsim_recapitation_applied": normal,
                "post_slim_neutral_mutation_overlay_applied": normal,
                **evaluation,
                "elapsed_seconds": perf_counter() - started,
                "error": evaluation.get("identity_error", ""),
                "event_audit_json": json.dumps(
                    model_record["event_audit"],
                    sort_keys=True,
                    separators=(",", ":"),
                ),
                "trajectory_path": trajectory_path.name,
                "trajectory_sha256": sha256_file(trajectory_path),
            }
        except SimulationDrawTimeout as error:
            result = {
                **row,
                "status": "timed_out",
                "evaluable": False,
                "terminal_class": "timed_out",
                "candidate_pool_alt_count": None,
                "candidate_pool_af": None,
                "candidate_pool_frequency_gate_passed": False,
                "exact_panel_feasible": False,
                "hit": False,
                "elapsed_seconds": perf_counter() - started,
                "error": str(error),
            }
        except Exception as error:  # noqa: BLE001 - preserve exact draw failure
            result = {
                **row,
                "status": "failed",
                "evaluable": False,
                "terminal_class": "failed",
                "candidate_pool_alt_count": None,
                "candidate_pool_af": None,
                "candidate_pool_frequency_gate_passed": False,
                "exact_panel_feasible": False,
                "hit": False,
                "elapsed_seconds": perf_counter() - started,
                "error": f"{type(error).__name__}: {error}",
            }
        payload = {
            "schema": DRAW_SCHEMA,
            "contract": contract,
            "contract_sha256": _canonical_sha256(contract),
            "result": result,
        }
        if result["evaluable"]:
            _atomic_json(completion_path, payload)
            (directory / "last_failure.json").unlink(missing_ok=True)
        else:
            _atomic_json(directory / "last_failure.json", payload)
        return result


def run_phase(
    phase_plan: pd.DataFrame,
    *,
    repo_root: str | Path,
    work_root: str | Path,
    slim_path: str | Path,
    manifest_sha256: str,
    implementation_sources: Mapping[str, str],
    runtime: Mapping[str, Any],
    workers: int = 4,
    timeout_seconds: float = 1_800.0,
) -> pd.DataFrame:
    if list(phase_plan.columns) != list(ROW_COLUMNS):
        raise ValueError("fallback phase plan columns are incompatible")
    phases = set(phase_plan["phase"].astype(str))
    if len(phases) != 1 or not phases.issubset(PHASE_DRAW_COUNTS):
        raise ValueError("fallback phase plan must contain one known phase")
    phase = next(iter(phases))
    expected = PHASE_DRAW_COUNTS[phase] * 6
    if len(phase_plan) != expected:
        raise ValueError("fallback phase plan draw count is incompatible")
    if (
        phase_plan["draw_id"].duplicated().any()
        or phase_plan["seed"].duplicated().any()
        or phase_plan["panel_seed"].duplicated().any()
    ):
        raise ValueError("fallback phase plan identities/seeds are not unique")
    if not 1 <= int(workers) <= 24:
        raise ValueError("fallback workers must be in [1, 24]")
    if not math.isfinite(float(timeout_seconds)) or float(timeout_seconds) <= 0:
        raise ValueError("fallback draw timeout must be finite and positive")
    if set(phase_plan["timeout_seconds"].astype(float)) != {float(timeout_seconds)}:
        raise ValueError("fallback run timeout differs from its atomic manifest")
    root = Path(repo_root).resolve()
    work = _require_repo_path(work_root, root, label="phase work directory")
    slim = _require_repo_path(slim_path, root, label="SLiM executable")
    current_sources = _implementation_sources(root)
    if dict(implementation_sources) != current_sources:
        raise ValueError("fallback implementation source binding drifted")
    current_runtime = _runtime_record(root, slim)
    if dict(runtime) != current_runtime:
        raise ValueError("fallback runtime binding drifted")
    tasks = [
        {
            "row": _normalized_manifest_row(row),
            "work_root": str(work),
            "slim_path": str(slim),
            "timeout_seconds": float(timeout_seconds),
            "manifest_sha256": str(manifest_sha256),
            "implementation_sources": current_sources,
            "runtime": current_runtime,
        }
        for row in phase_plan.to_dict(orient="records")
    ]
    rows = []
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
                        {
                            **task["row"],
                            "status": "failed",
                            "evaluable": False,
                            "terminal_class": "failed",
                            "candidate_pool_alt_count": None,
                            "candidate_pool_af": None,
                            "candidate_pool_frequency_gate_passed": False,
                            "exact_panel_feasible": False,
                            "hit": False,
                            "elapsed_seconds": np.nan,
                            "error": f"{type(error).__name__}: {error}",
                        }
                    )
    return pd.DataFrame(rows).sort_values("draw_id").reset_index(drop=True)


def _validate_result_row_fields(
    result: Mapping[str, Any], row: Mapping[str, Any]
) -> None:
    expected = _normalized_manifest_row(row)
    observed = _normalized_manifest_row(result)
    if observed != expected:
        raise ValueError("fallback completion manifest row differs from its plan")


def _validate_passing_panel_artifacts(
    result: Mapping[str, Any],
    row: Mapping[str, Any],
    directory: Path,
    expectation_cache: dict[tuple[Any, ...], dict[str, Any]],
) -> None:
    tree_path = directory / str(result.get("panel_tree_path", ""))
    manifest_path = directory / str(result.get("panel_manifest_path", ""))
    if (
        tree_path.name != "accepted_panel.trees"
        or manifest_path.name != "accepted_panel_manifest.tsv"
        or not tree_path.is_file()
        or not manifest_path.is_file()
        or sha256_file(tree_path) != str(result.get("panel_tree_sha256", ""))
        or sha256_file(manifest_path) != str(result.get("panel_manifest_sha256", ""))
    ):
        raise ValueError("fallback passing panel artifact checksum failed")
    key = (
        float(row["selection_coefficient"]),
        float(row["target_allele_frequency"]),
        int(row["sequence_length_bp"]),
        int(row["focal_position_bp"]),
    )
    if key not in expectation_cache:
        model, sweep, _, _ = build_frequency_conditioned_objects(row)
        expectation = _selected_focal_expectation(
            model, sweep, float(row["slim_scaling_factor"])
        )
        expectation["focal_position_0based"] = int(row["focal_position_bp"])
        expectation_cache[key] = expectation
    expectation = expectation_cache[key]
    panel = tskit.load(tree_path)
    coordinate = int(row["focal_position_bp"])
    raw = str(row["processing_mode"]) == FOCAL_ONLY_MODE
    if raw:
        identity = _inspect_raw_type1(panel, expectation)
        counts = _raw_type1_counts(panel, identity)
        rebuilt = build_diploid_pair_table(
            panel, coordinate, focal_genotype_counts=counts
        )
    else:
        validate_selected_focal_identity(panel, expectation)
        rebuilt = build_diploid_pair_table(panel, coordinate)
        counts = rebuilt["focal_selected_allele_count"].to_numpy(dtype=int)
    persisted = pd.read_csv(manifest_path, sep="\t")
    pd.testing.assert_frame_equal(
        rebuilt.reset_index(drop=True),
        persisted.reset_index(drop=True),
        check_dtype=False,
        check_exact=False,
        rtol=0.0,
        atol=1e-12,
    )
    n0, n1, n2 = (int(np.count_nonzero(counts == value)) for value in (0, 1, 2))
    if (
        len(counts) != int(row["sample_diploids"])
        or int(counts.sum()) != int(row["exact_sample_alt_count"])
        or n0 < 2
        or n1 < 1
        or n2 < 2
        or int(result.get("panel_seed", -1)) != int(row["panel_seed"])
        or str(result.get("panel_sampler_source"))
        != (
            "calibration_local_position_aware_equivalent"
            if raw
            else "focused_selection_simulation.select_exact_sample_panel"
        )
        or int(result.get("panel_alt_count", -1)) != int(counts.sum())
        or not math.isclose(
            float(result.get("panel_af", math.nan)),
            float(counts.mean() / 2),
            abs_tol=1e-12,
        )
        or (n0, n1, n2)
        != (
            int(result.get("panel_hom_ref", -1)),
            int(result.get("panel_heterozygous", -1)),
            int(result.get("panel_hom_alt", -1)),
        )
    ):
        raise ValueError("fallback passing panel semantics differ")


def _validate_evaluable_result(
    result: Mapping[str, Any],
    row: Mapping[str, Any],
    directory: Path,
    expectation_cache: dict[tuple[Any, ...], dict[str, Any]],
) -> None:
    if str(result.get("status")) != "complete" or not _coerce_boolean(
        result.get("evaluable")
    ):
        raise ValueError("fallback completion is not evaluable")
    for column in (
        "internal_terminal_han_af_conditioning",
        "selection_through_present",
        "recipient_nonloss_conditioning",
    ):
        if not _coerce_boolean(result.get(column)):
            raise ValueError(f"fallback completion lacks {column}")
    normal = str(row["phase"]) == "sensitivity_10mb"
    if _coerce_boolean(result.get("stdpopsim_recapitation_applied")) != normal or (
        _coerce_boolean(result.get("post_slim_neutral_mutation_overlay_applied"))
        != normal
    ):
        raise ValueError("fallback completion processing mode differs")
    trajectory = directory / str(result.get("trajectory_path", ""))
    if (
        trajectory.name != "trajectory.csv"
        or not trajectory.is_file()
        or sha256_file(trajectory) != str(result.get("trajectory_sha256", ""))
    ):
        raise ValueError("fallback trajectory checksum failed")
    event_audit = json.loads(str(result.get("event_audit_json", "")))
    event_key = (
        "event_audit",
        float(row["selection_coefficient"]),
        float(row["target_allele_frequency"]),
        int(row["sequence_length_bp"]),
        int(row["focal_position_bp"]),
    )
    if event_key not in expectation_cache:
        _, _, _, expected_model_record = build_frequency_conditioned_objects(row)
        expectation_cache[event_key] = expected_model_record["event_audit"]
    if (
        event_audit != expectation_cache[event_key]
        or event_audit.get("schema") != "gamma-smc.han-frequency-conditioned-events/v1"
        or event_audit.get("focal_event_count") != 9
        or event_audit.get("selection_cessation_applied") is not False
        or event_audit.get("external_terminal_condition_removal_applied") is not False
        or event_audit.get("realized_event_schedule_sha256")
        != _canonical_sha256(event_audit.get("realized_q5_event_times_generations"))
        or event_audit.get("serialized_extended_events_sha256")
        != _canonical_sha256(event_audit.get("serialized_extended_events"))
    ):
        raise ValueError("fallback completion event audit differs")
    strict_identity = _coerce_boolean(result.get("strict_type1_identity_passed"))
    gate = _coerce_boolean(result.get("candidate_pool_frequency_gate_passed"))
    feasible = _coerce_boolean(result.get("exact_panel_feasible"))
    roundtrip = _coerce_boolean(result.get("panel_roundtrip_validated"))
    hit = _coerce_boolean(result.get("hit"))
    if hit != bool(strict_identity and gate and feasible and roundtrip):
        raise ValueError("fallback hit is inconsistent with its four acceptance gates")
    tree_path = directory / "accepted_panel.trees"
    panel_manifest = directory / "accepted_panel_manifest.tsv"
    if not strict_identity:
        if not normal or str(result.get("terminal_class")) != "identity_collision":
            raise ValueError("only sensitivity identity collisions may be evaluable")
        if gate or feasible or roundtrip or hit:
            raise ValueError("fallback identity collision passed a later gate")
        if tree_path.exists() or panel_manifest.exists():
            raise ValueError("fallback rejected identity retained panel artifacts")
        return
    pool_alt = int(result.get("candidate_pool_alt_count", -1))
    pool_af = float(result.get("candidate_pool_af", math.nan))
    haplotypes = 2 * int(row["candidate_pool_diploids"])
    if not 0 <= pool_alt <= haplotypes or not math.isclose(
        pool_af, pool_alt / haplotypes, abs_tol=1e-12
    ):
        raise ValueError("fallback candidate-pool AF/count differ")
    expected_class = (
        "loss"
        if pool_alt == 0
        else "fixed"
        if pool_alt == haplotypes
        else "segregating"
    )
    if str(result.get("terminal_class")) != expected_class:
        raise ValueError("fallback terminal class differs from pool count")
    gate_record = json.loads(str(result.get("candidate_pool_frequency_gate_json", "")))
    if (
        _coerce_boolean(gate_record.get("passed")) != gate
        or int(gate_record.get("observed_alt_count", -1)) != pool_alt
        or not math.isclose(
            float(gate_record.get("observed_af", math.nan)), pool_af, abs_tol=1e-12
        )
    ):
        raise ValueError("fallback candidate-pool gate record differs")
    classes = gate_record.get("candidate_pool_genotype_counts")
    if not isinstance(classes, Mapping):
        raise ValueError("fallback pool genotype counts are absent")
    n0 = int(classes.get("hom_ref", -1))
    n1 = int(classes.get("heterozygous", -1))
    n2 = int(classes.get("hom_alt", -1))
    if n0 + n1 + n2 != int(row["candidate_pool_diploids"]) or n1 + 2 * n2 != pool_alt:
        raise ValueError("fallback pool genotype counts differ")
    reconstructed = np.repeat(np.asarray([0, 1, 2]), [n0, n1, n2])
    expected_feasible = (
        bool(
            feasible_exact_panel_genotype_triples(
                reconstructed,
                sample_diploids=int(row["sample_diploids"]),
                exact_alt_count=int(row["exact_sample_alt_count"]),
            )
        )
        if gate
        else False
    )
    if feasible != expected_feasible:
        raise ValueError("fallback exact-panel feasibility differs")
    if hit:
        _validate_passing_panel_artifacts(result, row, directory, expectation_cache)
    elif (
        roundtrip
        or result.get("panel_tree_path")
        or result.get("panel_manifest_path")
        or tree_path.exists()
        or panel_manifest.exists()
    ):
        raise ValueError("fallback rejected draw retained a panel artifact")


def rebuild_phase_ledger(
    phase_plan: pd.DataFrame,
    *,
    repo_root: str | Path,
    work_root: str | Path,
    manifest_sha256: str,
    implementation_sources: Mapping[str, str],
    runtime: Mapping[str, Any],
) -> pd.DataFrame:
    """Rebuild authoritative phase evidence from per-draw contract files."""

    root = Path(repo_root).resolve()
    work = _require_repo_path(work_root, root, label="phase work directory")
    expectation_cache: dict[tuple[Any, ...], dict[str, Any]] = {}
    results = []
    for raw_row in phase_plan.to_dict(orient="records"):
        row = _normalized_manifest_row(raw_row)
        directory = work / row["draw_id"]
        completion = directory / "completion.json"
        failure = directory / "last_failure.json"
        present = [path for path in (completion, failure) if path.is_file()]
        if len(present) != 1:
            raise ValueError(
                f"fallback draw must have exactly one terminal record: {row['draw_id']}"
            )
        try:
            payload = json.loads(present[0].read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            raise ValueError(
                f"fallback draw record is unreadable: {row['draw_id']}"
            ) from error
        expected_contract = _draw_contract(
            row,
            manifest_sha256=manifest_sha256,
            implementation_sources=implementation_sources,
            runtime=runtime,
        )
        if (
            payload.get("schema") != DRAW_SCHEMA
            or payload.get("contract") != expected_contract
            or payload.get("contract_sha256") != _canonical_sha256(expected_contract)
        ):
            raise ValueError(f"fallback draw contract differs: {row['draw_id']}")
        result = payload.get("result")
        if not isinstance(result, Mapping):
            raise ValueError(f"fallback draw result is absent: {row['draw_id']}")
        _validate_result_row_fields(result, row)
        evaluable = _coerce_boolean(result.get("evaluable"))
        if present[0] == completion:
            if not evaluable:
                raise ValueError(
                    "fallback completion file contains an unevaluable draw"
                )
            _validate_evaluable_result(result, row, directory, expectation_cache)
        elif evaluable or str(result.get("status")) not in {"failed", "timed_out"}:
            raise ValueError("fallback failure file contains an evaluable draw")
        results.append(dict(result))
    return pd.DataFrame(results).sort_values("draw_id").reset_index(drop=True)


def summarize_phase(
    phase_plan: pd.DataFrame,
    ledger: pd.DataFrame,
    *,
    work_root: str | Path,
    minimum_hit_rate: float = MINIMUM_HIT_RATE,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    if list(phase_plan.columns) != list(ROW_COLUMNS):
        raise ValueError("fallback phase plan columns are incompatible")
    phases = set(phase_plan["phase"].astype(str))
    if len(phases) != 1:
        raise ValueError("fallback phase plan contains multiple phases")
    phase = next(iter(phases))
    if phase not in PHASE_DRAW_COUNTS:
        raise ValueError("fallback phase is unknown")
    expected_rows = (
        PHASE_DRAW_COUNTS[phase] * len(SELECTION_COEFFICIENTS) * len(TARGET_FREQUENCIES)
    )
    if len(phase_plan) != expected_rows:
        raise ValueError("fallback phase plan does not cover all six cells")
    if set(phase_plan["draw_id"]) != set(ledger["draw_id"]):
        raise ValueError("fallback phase ledger does not cover its plan exactly")
    if ledger["draw_id"].duplicated().any():
        raise ValueError("fallback phase ledger contains duplicate draw IDs")
    plan_indexed = phase_plan.set_index("draw_id")
    for row in ledger.to_dict(orient="records"):
        identity = str(row["draw_id"])
        planned = plan_indexed.loc[identity].to_dict()
        planned["draw_id"] = identity
        _validate_result_row_fields(row, planned)
    evaluable = ledger["evaluable"].map(_coerce_boolean)
    hits = ledger["hit"].map(_coerce_boolean)
    if (hits & ~evaluable).any():
        raise ValueError("unevaluable fallback draw cannot be a hit")
    for row in ledger.loc[evaluable].to_dict(orient="records"):
        expected_hit = bool(
            _coerce_boolean(row["strict_type1_identity_passed"])
            and _coerce_boolean(row["candidate_pool_frequency_gate_passed"])
            and _coerce_boolean(row["exact_panel_feasible"])
            and _coerce_boolean(row["panel_roundtrip_validated"])
        )
        if _coerce_boolean(row["hit"]) != expected_hit:
            raise ValueError("fallback hit differs from its four acceptance gates")
    completed = ledger.loc[evaluable]
    if len(completed):
        required_true = (
            "internal_terminal_han_af_conditioning",
            "selection_through_present",
            "recipient_nonloss_conditioning",
        )
        for column in required_true:
            if (
                column not in completed
                or not completed[column].map(_coerce_boolean).all()
            ):
                raise ValueError(f"fallback completed draw lacks {column}")
        normal = completed["phase"].astype(str).eq("sensitivity_10mb")
        if (
            not (
                completed["stdpopsim_recapitation_applied"].map(_coerce_boolean)
                == normal
            ).all()
            or not (
                completed["post_slim_neutral_mutation_overlay_applied"].map(
                    _coerce_boolean
                )
                == normal
            ).all()
        ):
            raise ValueError("fallback phase postprocessing differs from its contract")
    work = Path(work_root).resolve()
    for row in ledger.loc[hits].to_dict(orient="records"):
        directory = work / str(row["draw_id"])
        for path_key, hash_key in (
            ("panel_tree_path", "panel_tree_sha256"),
            ("panel_manifest_path", "panel_manifest_sha256"),
        ):
            path = directory / str(row[path_key])
            if not path.is_file() or sha256_file(path) != str(row[hash_key]):
                raise ValueError("fallback passing panel artifact failed checksum")
        if not _coerce_boolean(row["panel_roundtrip_validated"]):
            raise ValueError("fallback passing panel lacks roundtrip validation")
    merged = phase_plan.merge(
        ledger[
            [
                "draw_id",
                "evaluable",
                "terminal_class",
                "candidate_pool_af",
                "candidate_pool_frequency_gate_passed",
                "exact_panel_feasible",
                "hit",
                "elapsed_seconds",
                "error",
            ]
        ],
        on="draw_id",
        how="left",
        validate="one_to_one",
    )
    rows = []
    for keys, group in merged.groupby(
        [
            "phase",
            "selection_coefficient",
            "target_allele_frequency",
            "population_af_lower",
            "population_af_upper",
        ],
        sort=True,
    ):
        phase_value, coefficient, target, lower, upper = keys
        group_evaluable = group["evaluable"].map(_coerce_boolean)
        group_hits = group["hit"].map(_coerce_boolean) & group_evaluable
        af = pd.to_numeric(
            group.loc[group_evaluable, "candidate_pool_af"], errors="coerce"
        )
        rows.append(
            {
                "phase": str(phase_value),
                "selection_coefficient": float(coefficient),
                "target_allele_frequency": float(target),
                "population_af_lower": float(lower),
                "population_af_upper": float(upper),
                "n_planned": int(len(group)),
                "n_evaluable": int(group_evaluable.sum()),
                "n_hits": int(group_hits.sum()),
                # Denominator is the fixed plan, so execution failures cannot
                # inflate a cell while the separate complete gate fails closed.
                "hit_rate": float(group_hits.sum() / len(group)),
                "n_identity_collision": int(
                    group["terminal_class"].astype(str).eq("identity_collision").sum()
                ),
                "n_loss": int(group["terminal_class"].astype(str).eq("loss").sum()),
                "n_fixed": int(group["terminal_class"].astype(str).eq("fixed").sum()),
                "median_candidate_pool_af": (
                    float(af.median()) if af.notna().any() else np.nan
                ),
                "total_elapsed_seconds": float(
                    pd.to_numeric(group["elapsed_seconds"], errors="coerce").sum()
                ),
            }
        )
    summary = pd.DataFrame(rows)
    expected_per_cell = PHASE_DRAW_COUNTS[phase]
    complete = bool(
        (summary["n_planned"] == expected_per_cell).all()
        and (summary["n_evaluable"] == expected_per_cell).all()
    )
    adequate = bool((summary["hit_rate"] >= float(minimum_hit_rate)).all())
    audit = {
        "schema": AUDIT_SCHEMA,
        "phase": phase,
        "status": "passed" if complete and adequate else "failed",
        "model_mode": MODEL_MODE,
        "draws_per_cell": expected_per_cell,
        "minimum_hit_rate": float(minimum_hit_rate),
        "all_draws_evaluable": complete,
        "all_six_cells_passed": adequate,
        "gamma_smc_statistics_used": False,
        "failed_fixed_cessation_observations_used": False,
        "identity_collisions_count_as_hits": False,
    }
    return summary, audit


def write_phase_audit(
    phase: str,
    *,
    repo_root: str | Path,
    output_dir: str | Path,
    work_root: str | Path,
    manifest_path: str | Path,
) -> dict[str, Any]:
    if phase not in PHASE_DRAW_COUNTS:
        raise ValueError("fallback phase is unknown")
    root = Path(repo_root).resolve()
    output = _require_repo_path(output_dir, root, label="phase-audit output")
    manifest, frame, manifest_sha256 = load_manifest(
        manifest_path, repo_root=root, require_current_sources=True
    )
    phase_plan = frame[frame["phase"].astype(str).eq(phase)].reset_index(drop=True)
    ledger = rebuild_phase_ledger(
        phase_plan,
        repo_root=root,
        work_root=work_root,
        manifest_sha256=manifest_sha256,
        implementation_sources=manifest["implementation_sources"],
        runtime=manifest["runtime"],
    )
    summary, audit = summarize_phase(phase_plan, ledger, work_root=work_root)
    ledger_path = output / f"{phase}_draws.tsv"
    summary_path = output / f"{phase}_summary.tsv"
    audit_path = output / f"{phase}_audit.json"
    _atomic_frame(ledger_path, ledger.sort_values("draw_id"))
    _atomic_frame(summary_path, summary)
    payload = {
        **audit,
        "manifest_sha256": manifest_sha256,
        "runtime": manifest["runtime"],
        "implementation_sources": manifest["implementation_sources"],
        "artifacts": {
            ledger_path.name: sha256_file(ledger_path),
            summary_path.name: sha256_file(summary_path),
        },
    }
    _atomic_json(audit_path, payload)
    return payload


def _require_passed_audit(
    path: Path,
    *,
    phase: str,
    manifest_sha256: str,
    runtime: Mapping[str, Any],
    implementation_sources: Mapping[str, str],
) -> dict[str, Any]:
    try:
        audit = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"fallback {phase} audit is unreadable") from error
    if (
        audit.get("schema") != AUDIT_SCHEMA
        or audit.get("phase") != phase
        or audit.get("status") != "passed"
        or audit.get("manifest_sha256") != manifest_sha256
        or audit.get("runtime") != runtime
        or audit.get("implementation_sources") != implementation_sources
    ):
        raise ValueError(f"fallback {phase} audit has not passed")
    artifacts = audit.get("artifacts")
    if not isinstance(artifacts, Mapping) or set(artifacts) != {
        f"{phase}_draws.tsv",
        f"{phase}_summary.tsv",
    }:
        raise ValueError(f"fallback {phase} audit artifacts are incompatible")
    for name, digest in artifacts.items():
        child = path.parent / str(name)
        if not child.is_file() or sha256_file(child) != str(digest):
            raise ValueError(f"fallback {phase} audited child checksum failed")
    return audit


def finalize_frequency_conditioned(
    *,
    manifest_path: str | Path,
    repo_root: str | Path,
    output_dir: str | Path,
    phase_work_roots: Mapping[str, str | Path],
) -> dict[str, Any]:
    root = Path(repo_root).resolve()
    output = _require_repo_path(output_dir, root, label="frozen output")
    manifest, frame, manifest_sha256 = load_manifest(
        manifest_path, repo_root=root, require_current_sources=True
    )
    summaries = {}
    artifacts = {Path(manifest_path).name: manifest_sha256}
    for phase in PHASE_DRAW_COUNTS:
        phase_plan = frame[frame["phase"].eq(phase)].reset_index(drop=True)
        ledger_path = output / f"{phase}_draws.tsv"
        summary_path = output / f"{phase}_summary.tsv"
        audit_path = output / f"{phase}_audit.json"
        audit = _require_passed_audit(
            audit_path,
            phase=phase,
            manifest_sha256=manifest_sha256,
            runtime=manifest["runtime"],
            implementation_sources=manifest["implementation_sources"],
        )
        ledger = rebuild_phase_ledger(
            phase_plan,
            repo_root=root,
            work_root=phase_work_roots[phase],
            manifest_sha256=manifest_sha256,
            implementation_sources=manifest["implementation_sources"],
            runtime=manifest["runtime"],
        )
        persisted_ledger = pd.read_csv(ledger_path, sep="\t")
        pd.testing.assert_frame_equal(
            ledger.reset_index(drop=True),
            persisted_ledger.sort_values("draw_id").reset_index(drop=True),
            check_dtype=False,
            check_exact=False,
            rtol=0.0,
            atol=1e-12,
        )
        summary, recomputed = summarize_phase(
            phase_plan, ledger, work_root=phase_work_roots[phase]
        )
        if recomputed["status"] != "passed":
            raise ValueError(f"fallback {phase} evidence no longer passes")
        persisted = pd.read_csv(summary_path, sep="\t")
        pd.testing.assert_frame_equal(
            summary.reset_index(drop=True),
            persisted.reset_index(drop=True),
            check_dtype=False,
            check_exact=False,
            rtol=0.0,
            atol=1e-12,
        )
        expected_artifacts = audit.get("artifacts")
        if not isinstance(expected_artifacts, Mapping):
            raise ValueError(f"fallback {phase} artifact manifest is absent")
        for path in (ledger_path, summary_path):
            if sha256_file(path) != expected_artifacts.get(path.name):
                raise ValueError(f"fallback {phase} child checksum failed")
            artifacts[path.name] = sha256_file(path)
        artifacts[audit_path.name] = sha256_file(audit_path)
        summaries[phase] = summary
    cells = {}
    for coefficient in SELECTION_COEFFICIENTS:
        for target in TARGET_FREQUENCIES:
            key = f"s={coefficient:.6f}|af={target:.6f}"
            entry = {
                "selection_coefficient": coefficient,
                "target_allele_frequency": target,
                "population_af_lower": target - AF_HALF_WIDTH,
                "population_af_upper": target + AF_HALF_WIDTH,
                "selection_start_generations_ago": INTROGRESSION_PULSE_GENERATIONS,
                "selection_end_generations_ago": 0.0,
                "terminal_han_frequency_conditioning_in_slim": True,
            }
            for phase, summary in summaries.items():
                row = summary[
                    np.isclose(summary["selection_coefficient"], coefficient)
                    & np.isclose(summary["target_allele_frequency"], target)
                ]
                if len(row) != 1:
                    raise ValueError("fallback phase summary lacks a unique cell")
                record = row.iloc[0]
                entry[f"{phase}_hit_rate"] = float(record["hit_rate"])
                entry[f"{phase}_n_hits"] = int(record["n_hits"])
                entry[f"{phase}_n_evaluable"] = int(record["n_evaluable"])
            cells[key] = entry
    payload = {
        "schema": FROZEN_SCHEMA,
        "status": "frozen",
        "manifest_sha256": manifest_sha256,
        "model_mode": MODEL_MODE,
        "production_integration_status": "not_wired_into_campaign_executor",
        "production_integration_requirement": (
            "future production code must explicitly support endpoint=0 and retain "
            "the two terminal Han AF conditions inside SLiM"
        ),
        "selection_endpoint_generations_ago": 0.0,
        "terminal_han_frequency_conditioning_in_slim": True,
        "post_pulse_recipient_nonloss_conditioning": True,
        "mutation_age_generations": INTROGRESSION_MUTATION_AGE_GENERATIONS,
        "slim_scaling_factor": SLIM_SCALING_FACTOR,
        "candidate_pool_diploids": CANDIDATE_POOL_DIPLOIDS,
        "sample_diploids": SAMPLE_DIPLOIDS,
        "af_half_width": AF_HALF_WIDTH,
        "minimum_hit_rate": MINIMUM_HIT_RATE,
        "all_three_phase_gates_passed": True,
        "gamma_smc_statistics_used": False,
        "failed_fixed_cessation_observations_used": False,
        "seed_exclusion_evidence": manifest["excluded_seed_artifacts"],
        "failed_fixed_cessation_evidence": manifest["failed_fixed_cessation_evidence"],
        "seed_registry": manifest["seed_registry"],
        "runtime": manifest["runtime"],
        "implementation_sources": manifest["implementation_sources"],
        "source_snapshots": manifest["source_snapshots"],
        "cells": cells,
        "artifacts": dict(sorted(artifacts.items())),
    }
    frozen_path = output / "han_frequency_conditioned_calibration_frozen.json"
    if frozen_path.is_file():
        existing = json.loads(frozen_path.read_text(encoding="utf-8"))
        if existing != payload:
            raise ValueError("existing fallback frozen contract differs")
        return existing
    _atomic_json(frozen_path, payload)
    return payload


def load_frozen_frequency_conditioned(
    path: str | Path,
    *,
    repo_root: str | Path,
    manifest_path: str | Path,
) -> dict[str, Any]:
    """Load the separate fallback authorization artifact fail-closed."""

    root = Path(repo_root).resolve()
    frozen_path = _require_repo_path(path, root, label="frozen contract")
    try:
        payload = json.loads(frozen_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError("fallback frozen contract is unreadable") from error
    manifest, _, manifest_sha256 = load_manifest(
        manifest_path, repo_root=root, require_current_sources=True
    )
    if (
        payload.get("schema") != FROZEN_SCHEMA
        or payload.get("status") != "frozen"
        or payload.get("manifest_sha256") != manifest_sha256
        or payload.get("model_mode") != MODEL_MODE
        or payload.get("production_integration_status")
        != "not_wired_into_campaign_executor"
        or payload.get("selection_endpoint_generations_ago") != 0.0
        or payload.get("terminal_han_frequency_conditioning_in_slim") is not True
        or payload.get("all_three_phase_gates_passed") is not True
        or payload.get("gamma_smc_statistics_used") is not False
        or payload.get("failed_fixed_cessation_observations_used") is not False
        or "mapping" in payload
    ):
        raise ValueError("fallback frozen contract semantics differ")
    cells = payload.get("cells")
    if not isinstance(cells, Mapping) or len(cells) != 6:
        raise ValueError("fallback frozen contract does not contain six cells")
    for record in cells.values():
        if (
            float(record.get("selection_end_generations_ago", math.nan)) != 0.0
            or record.get("terminal_han_frequency_conditioning_in_slim") is not True
            or any(
                float(record.get(f"{phase}_hit_rate", -1.0)) < MINIMUM_HIT_RATE
                for phase in PHASE_DRAW_COUNTS
            )
        ):
            raise ValueError("fallback frozen cell failed its phase gates")
    if payload.get("runtime") != manifest.get("runtime"):
        raise ValueError("fallback frozen runtime differs from its manifest")
    if payload.get("implementation_sources") != manifest.get("implementation_sources"):
        raise ValueError("fallback frozen source binding differs")
    artifacts = payload.get("artifacts")
    if not isinstance(artifacts, Mapping) or not artifacts:
        raise ValueError("fallback frozen child artifacts are absent")
    for name, expected in artifacts.items():
        child = frozen_path.parent / str(name)
        if not child.is_file() or sha256_file(child) != str(expected):
            raise ValueError(f"fallback frozen child checksum failed: {name}")
    return payload


def _campaign_paths(
    repo_root: Path, campaign_dir: str | Path
) -> tuple[Path, Path, Path]:
    campaign = Path(campaign_dir)
    if not campaign.is_absolute():
        campaign = repo_root / campaign
    campaign = _require_repo_path(campaign, repo_root, label="campaign directory")
    output = _require_repo_path(
        campaign / "calibration/han_frequency_conditioned",
        repo_root,
        label="calibration output",
    )
    work = _require_repo_path(
        campaign / "work/han_frequency_conditioned",
        repo_root,
        label="calibration work",
    )
    return campaign, output, work


def _default_seed_exclusion_paths(campaign: Path) -> list[Path]:
    return [
        campaign / "execution_units.tsv",
        campaign / "calibration/han_selection_end_plan.tsv",
        campaign / "calibration/han_selection_end_extension_plan.tsv",
        campaign / "calibration/han_selection_end_reconciliation_plan.tsv",
    ]


def _default_failed_fixed_cessation_paths(campaign: Path) -> list[Path]:
    calibration = campaign / "calibration"
    return [
        calibration / "han_selection_end_terminal_draws_reconciled.tsv",
        calibration / "han_selection_end_extension_draws.tsv",
        calibration / "han_selection_end_extension_select.log",
        calibration / "han_selection_end_reconciliation_audit.json",
        calibration / "fixed_cessation_q5_not_frozen/NOT_FROZEN.json",
    ]


def _phase_prerequisites(phase: str) -> tuple[str, ...]:
    return {
        "screen": (),
        "confirmation": ("screen",),
        "sensitivity_10mb": ("screen", "confirmation"),
    }[phase]


def _require_phase_prerequisites(
    phase: str,
    *,
    output: Path,
    manifest_sha256: str,
    manifest: Mapping[str, Any],
) -> None:
    for prerequisite in _phase_prerequisites(phase):
        _require_passed_audit(
            output / f"{prerequisite}_audit.json",
            phase=prerequisite,
            manifest_sha256=manifest_sha256,
            runtime=manifest["runtime"],
            implementation_sources=manifest["implementation_sources"],
        )


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Plan, run, and audit the isolated Han terminal-AF-conditioned "
            "fallback calibration."
        )
    )
    parser.add_argument(
        "command",
        choices=(
            "plan",
            "screen-run",
            "screen-audit",
            "confirmation-run",
            "confirmation-audit",
            "sensitivity-run",
            "sensitivity-audit",
            "finalize",
        ),
    )
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument(
        "--campaign-dir", type=Path, default=Path("focused_selection_EAS_sim")
    )
    parser.add_argument(
        "--slim-bin",
        type=Path,
        default=Path(".native-stdpopsim/bin/slim"),
        help="stdpopsim-compatible SLiM executable; hash/version are frozen at plan",
    )
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument(
        "--timeout-seconds", type=float, default=DEFAULT_DRAW_TIMEOUT_SECONDS
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    if not 1 <= int(args.workers) <= 24:
        raise ValueError("fallback workers must be in [1, 24]")
    if not math.isclose(
        float(args.timeout_seconds), DEFAULT_DRAW_TIMEOUT_SECONDS, abs_tol=1e-12
    ):
        raise ValueError("fallback timeout is frozen at 1800 seconds per draw")
    root = Path(args.repo_root).resolve()
    campaign, output, work = _campaign_paths(root, args.campaign_dir)
    slim = Path(args.slim_bin)
    if not slim.is_absolute():
        slim = root / slim
    slim = _require_repo_path(slim, root, label="SLiM executable")
    manifest_path = output / "han_frequency_conditioned_manifest.json"
    if args.command == "plan":
        manifest = build_atomic_manifest(
            FrequencyConditionedDesign(),
            repo_root=root,
            slim_path=slim,
            excluded_seed_artifacts=_default_seed_exclusion_paths(campaign),
            failed_fixed_cessation_artifacts=(
                _default_failed_fixed_cessation_paths(campaign)
            ),
        )
        written = write_atomic_manifest(manifest, manifest_path, repo_root=root)
        print(
            json.dumps(
                {
                    "status": "planned",
                    "manifest": str(manifest_path),
                    "manifest_sha256": sha256_file(manifest_path),
                    "draws": len(written["rows"]),
                    "phase_draws_per_cell": PHASE_DRAW_COUNTS,
                    "all_phase_rng_seeds_reserved_atomically": True,
                },
                indent=2,
                sort_keys=True,
            )
        )
        return 0
    manifest, frame, manifest_sha256 = load_manifest(
        manifest_path, repo_root=root, require_current_sources=True
    )
    command_to_phase = {
        "screen-run": "screen",
        "screen-audit": "screen",
        "confirmation-run": "confirmation",
        "confirmation-audit": "confirmation",
        "sensitivity-run": "sensitivity_10mb",
        "sensitivity-audit": "sensitivity_10mb",
    }
    if args.command in command_to_phase:
        phase = command_to_phase[args.command]
        _require_phase_prerequisites(
            phase,
            output=output,
            manifest_sha256=manifest_sha256,
            manifest=manifest,
        )
        phase_plan = frame[frame["phase"].astype(str).eq(phase)].reset_index(drop=True)
        phase_work = work / phase
        if args.command.endswith("-run"):
            ledger = run_phase(
                phase_plan,
                repo_root=root,
                work_root=phase_work,
                slim_path=slim,
                manifest_sha256=manifest_sha256,
                implementation_sources=manifest["implementation_sources"],
                runtime=manifest["runtime"],
                workers=args.workers,
                timeout_seconds=args.timeout_seconds,
            )
            output.mkdir(parents=True, exist_ok=True)
            snapshot = output / f"{phase}_run_snapshot.tsv"
            _atomic_frame(snapshot, ledger)
            complete = bool(ledger["evaluable"].map(_coerce_boolean).all())
            print(
                json.dumps(
                    {
                        "phase": phase,
                        "status": "complete" if complete else "retry_required",
                        "planned": len(phase_plan),
                        "evaluable": int(
                            ledger["evaluable"].map(_coerce_boolean).sum()
                        ),
                        "snapshot": str(snapshot),
                        "snapshot_sha256": sha256_file(snapshot),
                    },
                    indent=2,
                    sort_keys=True,
                )
            )
            return 0 if complete else 2
        audit = write_phase_audit(
            phase,
            repo_root=root,
            output_dir=output,
            work_root=phase_work,
            manifest_path=manifest_path,
        )
        print(json.dumps(audit, indent=2, sort_keys=True))
        return 0 if audit["status"] == "passed" else 2
    phase_roots = {phase: work / phase for phase in PHASE_DRAW_COUNTS}
    payload = finalize_frequency_conditioned(
        manifest_path=manifest_path,
        repo_root=root,
        output_dir=output,
        phase_work_roots=phase_roots,
    )
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
