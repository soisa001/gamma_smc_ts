"""Simulation executors for the focused frequency-matched EAS campaign.

Selected trajectories reuse the validated stdpopsim/SLiM event models. Selected
and neutral candidates share a 500-diploid AF-band gate and the same uniform
exact-k 100-diploid panel sampler. Neutral units use real mutations on matched
genealogies generated with 200 recent generations of DTWF followed by SMC-prime.
The no-introgression null conditions a branch in the 5-Mb marginal tree of the
requested 10-Mb sequence. The AncientEurasia null uses a 20-Mb ascertainment
ancestry and indexes only branches reached through the Neanderthal pulse before
cropping a 10-Mb window. Thus its focal allele is archaic-specific and cannot
be ILS.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import platform
import signal
from collections.abc import Mapping, Sequence
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import contextmanager
from dataclasses import dataclass, replace
from itertools import pairwise
from pathlib import Path
from time import perf_counter
from typing import Any

import msprime
import numpy as np
import pandas as pd
import pyslim
import stdpopsim
import tskit

from .eas_sweep_analysis import (
    build_diploid_pair_table,
    tree_truth_profiles,
    write_gamma_pair_files,
)
from .eas_sweep_models import (
    ARCHAIC_POPULATION,
    EAS_POPULATION,
    FOCAL_POSITION_BP,
    FOCAL_SITE_ID,
    GENERATION_TIME_YEARS,
    HAN_SPLIT_GENERATIONS,
    HUMAN_NEANDERTHAL_SPLIT_GENERATIONS,
    INTROGRESSION_PULSE_GENERATIONS,
    INTROGRESSION_TARGET_POPULATION,
    MUTATION_RATE,
    RECOMBINATION_RATE,
    SEQUENCE_LENGTH_BP,
    TMRCA_THRESHOLDS_YEARS,
    build_de_novo_origin_age_grid,
    build_eas_demography_models,
    build_introgression_sweep_spec,
    build_no_introgression_sweep_spec,
    load_ancient_eurasia_model,
    load_phlash_eas_npz,
    serialize_extended_events,
)
from .eas_sweep_study import EAS_RESOURCE_SHA256

SCHEMA_VERSION = "gamma-smc.focused-selection-simulation/v5"
FOCAL_SEMANTIC_VALIDATION_SCHEMA = "gamma-smc.focused-focal-semantics/v1"
SELECTED_FOCAL_IDENTITY_SCHEMA = "gamma-smc.selected-focal-identity/v1"
CANDIDATE_POOL_FREQUENCY_GATE_SCHEMA = "gamma-smc.candidate-pool-frequency-gate/v1"
HAN_SELECTION_CALIBRATION_SCHEMA = "gamma-smc.han-selection-end-calibration/v1"
HAN_SELECTION_PRODUCTION_CONTRACT_SCHEMA = (
    "gamma-smc.han-selection-production-contract/v1"
)
SIMULATION_SOURCE_PATH = "python/gamma_smc_aou/focused_selection_simulation.py"
DEFAULT_HAN_SELECTION_CALIBRATION_PATH = (
    "focused_selection_EAS_sim/calibration/han_selection_end_frozen.json"
)
INTROGRESSION_ASCERTAINMENT_LENGTH_BP = 20_000_000
ASCERTAINMENT_INTERIOR_LEFT_BP = SEQUENCE_LENGTH_BP // 2
ASCERTAINMENT_INTERIOR_RIGHT_BP = (
    INTROGRESSION_ASCERTAINMENT_LENGTH_BP - SEQUENCE_LENGTH_BP // 2
)
NEUTRAL_RECENT_DTWF_DURATION_GENERATIONS = 200.0
NEUTRAL_ANCESTRY_MODEL_KEYWORDS = ("dtwf", "smc_prime")
NEUTRAL_ANCESTRY_MODEL_NAME = "DTWF200ThenSmcPrimeApproxCoalescent"
DEFAULT_SLIM_SCALING_FACTOR = 5.0
DEFAULT_SLIM_BURN_IN = 0.1
DEFAULT_NEUTRAL_MAX_ATTEMPTS = 100
DEFAULT_NEUTRAL_ANCESTRY_BATCH_SIZE = 4
DEFAULT_SELECTED_MAX_DRAWS = 50
DEFAULT_SELECTED_POOL_DIPLOIDS = 500
DEFAULT_NEUTRAL_POOL_DIPLOIDS = DEFAULT_SELECTED_POOL_DIPLOIDS
DEFAULT_SELECTED_DRAW_TIMEOUT_SECONDS = 5 * 60.0
DEFAULT_SELECTED_CUMULATIVE_TIMEOUT_SECONDS = 45 * 60.0
INTROGRESSION_SOURCE_PRESENCE_EPSILON = 1e-9
INTROGRESSION_MUTATION_AGE_GENERATIONS = 2_400.0
INTROGRESSION_PULSE_PROPORTION = 0.0296
EXPECTED_SCALED_ARCHAIC_DIPLOIDS_AT_PULSE = 8
IMPLEMENTATION_PATHS = (
    SIMULATION_SOURCE_PATH,
    "python/gamma_smc_aou/eas_sweep_models.py",
    "python/gamma_smc_aou/eas_sweep_analysis.py",
)


class SimulationExhausted(RuntimeError):
    """Raised when a unit consumes its frozen bounded search budget."""


class SimulationDrawTimeout(TimeoutError):
    """Raised when one stdpopsim/SLiM draw exceeds its wall-clock budget."""


class UnitLockHeld(RuntimeError):
    """Raised when another process already owns a unit-phase lock."""


class SelectedFocalIdentityError(ValueError):
    """Raised when a selected tree lacks the prescribed focal mutation identity."""


@dataclass(frozen=True)
class SimulationArtifacts:
    unit_dir: Path
    tree_path: Path
    pair_table_path: Path
    overall_pairs_path: Path
    truth_profiles_path: Path
    truth_class_summaries_path: Path
    completion_path: Path
    cache_hit: bool


@dataclass(frozen=True)
class NeutralCandidate:
    left: float
    right: float
    node: int
    mutation_time_lower: float
    mutation_time_upper: float
    mutation_origin_population: str
    candidate_pool_alt_count: int
    candidate_pool_genotype_counts: tuple[int, int, int]
    eligible_genotype_count_triples: int
    n_introgressed_descendant_lineages: int
    weight: float


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_sha256(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(
            value, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode()
    ).hexdigest()


def _stable_seed(base_seed: int, label: str) -> int:
    digest = hashlib.sha256(f"{base_seed}:{label}".encode()).digest()
    return int(1 + int.from_bytes(digest[:8], "big") % (2**31 - 2))


def _portable_path(path: Path, repo_root: Path) -> str:
    try:
        return path.resolve().relative_to(repo_root.resolve()).as_posix()
    except ValueError:
        return str(path.resolve())


def _han_endpoint_key(selection_coefficient: float, target_frequency: float) -> str:
    return f"s={float(selection_coefficient):.6f}|af={float(target_frequency):.6f}"


def _production_seed_digest(units: pd.DataFrame) -> str:
    seeds = sorted(int(value) for value in units["seed"])
    return hashlib.sha256(
        json.dumps(seeds, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def load_frozen_han_selection_calibration(
    path: str | Path,
    *,
    repo_root: str | Path,
    execution_units_path: str | Path,
    slim_scaling_factor: float,
    candidate_pool_diploids: int,
) -> dict[str, Any]:
    """Validate and bind the frozen Han endpoint calibration for production.

    The JSON is not treated as a loose parameter file. Its design contract,
    production-seed digest, selected endpoint table, and every calibration
    artifact checksum must agree before a Han selected draw can start.
    """

    root = Path(repo_root).resolve()
    frozen_path = Path(path).resolve()
    units_path = Path(execution_units_path).resolve()
    try:
        frozen_path.relative_to(root)
        units_path.relative_to(root)
    except ValueError as error:
        raise ValueError(
            "Han calibration and execution units must remain inside repo_root"
        ) from error
    if any("onedrive" in part.casefold() for part in frozen_path.parts):
        raise ValueError("Han calibration artifact must not be inside OneDrive")
    if not frozen_path.is_file():
        raise ValueError(f"frozen Han selection calibration is absent: {frozen_path}")
    try:
        payload = json.loads(frozen_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError("frozen Han selection calibration is unreadable") from error
    if (
        payload.get("schema") != HAN_SELECTION_CALIBRATION_SCHEMA
        or payload.get("status") != "frozen"
    ):
        raise ValueError(
            "Han selection calibration is not a frozen compatible artifact"
        )

    production = payload.get("production_contract")
    if (
        not isinstance(production, Mapping)
        or production.get("schema") != HAN_SELECTION_PRODUCTION_CONTRACT_SCHEMA
    ):
        raise ValueError("frozen Han selection production contract is absent")
    q = float(slim_scaling_factor)
    pool_size = int(candidate_pool_diploids)
    if not math.isclose(
        float(production.get("slim_scaling_factor", math.nan)),
        q,
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        raise ValueError(
            "frozen Han calibration Q differs from production SLiM scaling factor"
        )
    if int(production.get("candidate_pool_diploids", -1)) != pool_size:
        raise ValueError(
            "frozen Han calibration candidate-pool size differs from production"
        )
    if production.get("selection_endpoint_mode") != (
        "fixed_endpoint_from_disjoint_terminal_af_calibration"
    ):
        raise ValueError("frozen Han calibration endpoint mode is incompatible")
    if production.get("post_pulse_recipient_nonloss_conditioning") is not True:
        raise ValueError("frozen Han calibration lacks recipient nonloss conditioning")

    try:
        units = pd.read_csv(units_path, sep="\t")
    except Exception as error:
        raise ValueError("focused execution-unit table is unreadable") from error
    required_columns = {
        "seed",
        "demography_id",
        "simulation_class",
        "selection_coefficient",
        "target_allele_frequency",
        "population_af_lower",
        "population_af_upper",
    }
    if required_columns.difference(units.columns) or units.empty:
        raise ValueError("focused execution-unit table is incompatible")
    if units["seed"].duplicated().any():
        raise ValueError("focused execution-unit seeds are not unique")
    expected_seed_digest = _production_seed_digest(units)
    if payload.get("production_seed_digest") != expected_seed_digest:
        raise ValueError("frozen Han calibration production-seed digest is stale")
    han = units[
        (units["demography_id"].astype(str) == "ancient_eurasia_han_introgression")
        & (units["simulation_class"].astype(str) == "selected")
    ].copy()
    if han.empty:
        raise ValueError("execution-unit table lacks Han selected cells")
    coefficients = sorted(set(han["selection_coefficient"].astype(float)))
    targets = sorted(set(han["target_allele_frequency"].astype(float)))
    if coefficients != sorted(
        float(value) for value in production.get("selection_coefficients", [])
    ) or targets != sorted(
        float(value) for value in production.get("target_allele_frequencies", [])
    ):
        raise ValueError("frozen Han calibration grid differs from production")
    half_widths = np.r_[
        han["target_allele_frequency"].astype(float)
        - han["population_af_lower"].astype(float),
        han["population_af_upper"].astype(float)
        - han["target_allele_frequency"].astype(float),
    ]
    if not np.allclose(
        half_widths,
        float(production.get("af_half_width", math.nan)),
        rtol=0.0,
        atol=1e-12,
    ):
        raise ValueError("frozen Han calibration AF band differs from production")

    expected_keys = {
        _han_endpoint_key(coefficient, target)
        for coefficient in coefficients
        for target in targets
    }
    mapping = payload.get("mapping")
    if not isinstance(mapping, Mapping) or set(mapping) != expected_keys:
        raise ValueError("frozen Han calibration endpoint mapping is incomplete")
    normalized_mapping: dict[str, dict[str, Any]] = {}
    for key in sorted(expected_keys):
        entry = mapping[key]
        if not isinstance(entry, Mapping):
            raise TypeError(f"frozen Han endpoint entry is malformed: {key}")
        coefficient = float(entry.get("selection_coefficient", math.nan))
        target = float(entry.get("target_allele_frequency", math.nan))
        endpoint = float(entry.get("realized_selection_end_generations_ago", math.nan))
        if key != _han_endpoint_key(coefficient, target):
            raise ValueError(f"frozen Han endpoint key disagrees with its entry: {key}")
        if not (0.0 < endpoint < INTROGRESSION_PULSE_GENERATIONS) or not math.isclose(
            endpoint / q, round(endpoint / q), rel_tol=0.0, abs_tol=1e-9
        ):
            raise ValueError(f"frozen Han endpoint is outside the Q grid: {key}")
        normalized_mapping[key] = dict(entry)

    artifact_hashes = payload.get("artifacts")
    required_artifacts = {
        "han_selection_end_plan.tsv",
        "han_selection_end_terminal_draws.tsv",
        "han_selection_end_summary.tsv",
        "han_selection_end_selected.tsv",
    }
    if not isinstance(artifact_hashes, Mapping) or not required_artifacts.issubset(
        artifact_hashes
    ):
        raise ValueError("frozen Han calibration artifact manifest is incomplete")
    verified_artifacts: dict[str, str] = {}
    for name, expected_sha256 in sorted(artifact_hashes.items()):
        if Path(str(name)).name != str(name):
            raise ValueError("frozen Han calibration artifact name is unsafe")
        artifact_path = frozen_path.parent / str(name)
        if not artifact_path.is_file() or sha256_file(artifact_path) != expected_sha256:
            raise ValueError(f"frozen Han calibration artifact checksum failed: {name}")
        verified_artifacts[str(name)] = str(expected_sha256)
    selected_path = frozen_path.parent / "han_selection_end_selected.tsv"
    try:
        selected = pd.read_csv(selected_path, sep="\t")
    except Exception as error:
        raise ValueError("frozen Han selected-endpoint table is unreadable") from error
    selected_required = {
        "selection_coefficient",
        "target_allele_frequency",
        "realized_selection_end_generations_ago",
    }
    if selected_required.difference(selected.columns) or len(selected) != len(
        expected_keys
    ):
        raise ValueError("frozen Han selected-endpoint table is incompatible")
    table_mapping: dict[str, float] = {}
    for row in selected.itertuples(index=False):
        key = _han_endpoint_key(
            float(row.selection_coefficient), float(row.target_allele_frequency)
        )
        if key in table_mapping:
            raise ValueError("frozen Han selected-endpoint table has duplicate cells")
        table_mapping[key] = float(row.realized_selection_end_generations_ago)
    if set(table_mapping) != expected_keys or any(
        not math.isclose(
            table_mapping[key],
            float(normalized_mapping[key]["realized_selection_end_generations_ago"]),
            rel_tol=0.0,
            abs_tol=1e-12,
        )
        for key in expected_keys
    ):
        raise ValueError("frozen Han JSON mapping differs from selected-endpoint table")

    return {
        "schema": HAN_SELECTION_PRODUCTION_CONTRACT_SCHEMA,
        "calibration_schema": HAN_SELECTION_CALIBRATION_SCHEMA,
        "path": _portable_path(frozen_path, root),
        "sha256": sha256_file(frozen_path),
        "production_contract": dict(production),
        "production_seed_digest": expected_seed_digest,
        "artifact_sha256": verified_artifacts,
        "mapping_sha256": _canonical_sha256(normalized_mapping),
        "mapping": normalized_mapping,
    }


def resolve_han_selection_endpoint(
    calibration: Mapping[str, Any],
    *,
    selection_coefficient: float,
    target_frequency: float,
) -> dict[str, Any]:
    """Resolve one s-by-AF endpoint from an already validated frozen artifact."""

    if calibration.get("schema") != HAN_SELECTION_PRODUCTION_CONTRACT_SCHEMA:
        raise ValueError("validated Han calibration binding is incompatible")
    key = _han_endpoint_key(selection_coefficient, target_frequency)
    mapping = calibration.get("mapping")
    if not isinstance(mapping, Mapping) or key not in mapping:
        raise ValueError(f"frozen Han calibration lacks production cell {key}")
    entry = dict(mapping[key])
    return {
        "mapping_key": key,
        "selection_coefficient": float(selection_coefficient),
        "target_allele_frequency": float(target_frequency),
        "selection_end_generations_ago": float(
            entry["realized_selection_end_generations_ago"]
        ),
        "calibration_json_path": str(calibration["path"]),
        "calibration_json_sha256": str(calibration["sha256"]),
        "selected_endpoint_table_sha256": str(
            calibration["artifact_sha256"]["han_selection_end_selected.tsv"]
        ),
        "mapping_sha256": str(calibration["mapping_sha256"]),
        "calibration_entry": entry,
    }


def _implementation_contract(
    repo_root: Path,
    *,
    simulation_class: str,
    slim_path: str | Path | None,
) -> dict[str, Any]:
    sources: dict[str, str] = {}
    for relative in IMPLEMENTATION_PATHS:
        path = repo_root / relative
        if not path.is_file():
            raise ValueError(f"required implementation input is absent: {path}")
        sources[relative] = sha256_file(path)
    slim: dict[str, Any] | None = None
    if simulation_class == "selected":
        if slim_path is None:
            raise ValueError("selected simulation requires --slim-bin")
        executable = Path(slim_path).resolve()
        if not executable.is_file():
            raise ValueError(f"SLiM executable is absent: {executable}")
        slim = {
            "path": _portable_path(executable, repo_root),
            "sha256": sha256_file(executable),
        }
    return {
        "sources": sources,
        "slim": slim,
        "eas_resource_sha256": EAS_RESOURCE_SHA256,
    }


@contextmanager
def _wall_clock_timeout(seconds: float):
    """Interrupt a blocking SLiM subprocess wait on POSIX worker processes."""

    timeout = float(seconds)
    if not math.isfinite(timeout) or timeout <= 0:
        raise ValueError("wall-clock timeout must be finite and positive")
    if os.name != "posix":
        raise RuntimeError("selected SLiM timeout enforcement requires POSIX/WSL")
    previous_handler = signal.getsignal(signal.SIGALRM)

    def _raise_timeout(_signum, _frame):
        raise SimulationDrawTimeout(
            f"SLiM draw exceeded {timeout:.3f} wall-clock seconds"
        )

    signal.signal(signal.SIGALRM, _raise_timeout)
    previous_timer = signal.setitimer(signal.ITIMER_REAL, timeout)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, *previous_timer)
        signal.signal(signal.SIGALRM, previous_handler)


def _atomic_json(path: Path, payload: Any) -> None:
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
    compression: str | Mapping[str, Any] | None = None
    if path.suffix == ".gz":
        compression = {"method": "gzip", "compresslevel": 6, "mtime": 0}
    try:
        frame.to_csv(temporary, sep="\t", index=False, compression=compression)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_tree(path: Path, ts: tskit.TreeSequence) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        ts.dump(temporary)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _lock_conflict(error: OSError) -> bool:
    """Return whether an OS locking error means another owner holds the lock."""

    return isinstance(error, BlockingIOError) or error.errno in {
        11,  # EAGAIN on POSIX
        13,  # EACCES from msvcrt.locking on Windows
        35,  # EAGAIN on macOS
    }


@contextmanager
def exclusive_unit_lock(
    lock_path: str | Path,
    *,
    unit_id: str,
    phase: str,
):
    """Hold a nonblocking OS lock for one unit phase.

    The lock file is deliberately retained after release. Kernel-backed locks,
    rather than file existence, define ownership, so an interrupted worker is
    immediately restartable without stale-lock deletion.
    """

    path = Path(lock_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = path.open("a+b")
    locked = False
    backend = ""
    try:
        # Windows byte-range locks require the locked byte to exist. This is
        # harmless on POSIX and is done before either process claims ownership.
        handle.seek(0, os.SEEK_END)
        if handle.tell() == 0:
            handle.write(b"\0")
            handle.flush()
        handle.seek(0)
        try:
            if os.name == "posix":
                import fcntl

                fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                backend = "fcntl.flock"
            elif os.name == "nt":
                import msvcrt

                msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
                backend = "msvcrt.locking"
            else:  # pragma: no cover - supported production/test OSes are above
                raise RuntimeError(
                    f"unit locking is unsupported on os.name={os.name!r}"
                )
        except OSError as error:
            if _lock_conflict(error):
                raise UnitLockHeld(
                    f"{phase} lock is already held for unit {unit_id}: {path}"
                ) from error
            raise
        locked = True
        owner = {
            "schema": "gamma-smc.unit-lock/v1",
            "unit_id": str(unit_id),
            "phase": str(phase),
            "pid": os.getpid(),
            "hostname": platform.node(),
            "backend": backend,
        }
        handle.seek(0)
        handle.truncate()
        handle.write(
            (json.dumps(owner, sort_keys=True, allow_nan=False) + "\n").encode("utf-8")
        )
        handle.flush()
        os.fsync(handle.fileno())
        yield path
    finally:
        if locked:
            handle.seek(0)
            if os.name == "posix":
                import fcntl

                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
            elif os.name == "nt":
                import msvcrt

                msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
        handle.close()


def unit_lock_is_held(lock_path: str | Path) -> bool:
    """Return whether another process currently owns a unit lock."""

    path = Path(lock_path)
    if not path.is_file():
        return False
    with path.open("r+b") as handle:
        handle.seek(0, os.SEEK_END)
        if handle.tell() == 0:
            return False
        handle.seek(0)
        try:
            if os.name == "posix":
                import fcntl

                fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
            elif os.name == "nt":
                import msvcrt

                msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
                handle.seek(0)
                msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
            else:  # pragma: no cover - supported production/test OSes are above
                raise RuntimeError(
                    f"unit locking is unsupported on os.name={os.name!r}"
                )
        except OSError as error:
            if _lock_conflict(error):
                return True
            raise
    return False


def _validate_unit(unit: Mapping[str, Any]) -> dict[str, Any]:
    required = {
        "unit_id",
        "demography_id",
        "demography_kind",
        "population",
        "simulation_class",
        "selection_coefficient",
        "target_allele_frequency",
        "population_af_lower",
        "population_af_upper",
        "exact_sample_alt_count",
        "sample_diploids",
        "replicate_index",
        "seed",
        "sequence_length_bp",
        "focal_position_bp",
    }
    missing = sorted(required.difference(unit))
    if missing:
        raise ValueError(f"execution unit is missing fields: {', '.join(missing)}")
    record = {
        str(key): (value.item() if isinstance(value, np.generic) else value)
        for key, value in unit.items()
    }
    if int(record["sequence_length_bp"]) != SEQUENCE_LENGTH_BP:
        raise ValueError("focused simulation requires exactly 10 Mb")
    if int(record["focal_position_bp"]) != FOCAL_POSITION_BP:
        raise ValueError("focused simulation requires the 5-Mb focal position")
    sample_diploids = int(record["sample_diploids"])
    target = float(record["target_allele_frequency"])
    expected = round(2 * sample_diploids * target)
    if int(record["exact_sample_alt_count"]) != expected:
        raise ValueError("execution unit exact sample alternate count is inconsistent")
    simulation_class = str(record["simulation_class"])
    coefficient = float(record["selection_coefficient"])
    if simulation_class == "neutral" and coefficient != 0:
        raise ValueError("neutral units must have s=0")
    if simulation_class == "selected" and coefficient <= 0:
        raise ValueError("selected units must have positive s")
    if simulation_class not in {"neutral", "selected"}:
        raise ValueError("simulation_class must be selected or neutral")
    return record


def _population_name(demography: msprime.Demography, population_id: int) -> str:
    return str(demography.populations[int(population_id)].name)


def _tree_sequence_population_name(ts: tskit.TreeSequence, population_id: int) -> str:
    """Return a population name from a stdpopsim/msprime tree sequence."""

    population = ts.population(int(population_id))
    metadata = population.metadata
    if isinstance(metadata, Mapping):
        for key in ("name", "id"):
            value = metadata.get(key)
            if value:
                return str(value)
    raise ValueError(f"tree-sequence population {population_id} lacks name/id metadata")


def load_cell_demography(
    unit: Mapping[str, Any], repo_root: str | Path
) -> tuple[msprime.Demography, str, float]:
    """Load the exact cell demography and its Gamma-SMC present-day theta."""

    record = _validate_unit(unit)
    root = Path(repo_root).resolve()
    if record["demography_id"] == "eas_phlash_median":
        artifact = load_phlash_eas_npz(
            root / "no_introgression_EAS_sim" / "resources" / "EAS.npz",
            expected_sha256=EAS_RESOURCE_SHA256,
        )
        model = build_eas_demography_models(artifact)["median"]
        return model.msprime_demography, EAS_POPULATION, float(model.ne[0])
    if record["demography_id"] == "ancient_eurasia_han_introgression":
        model = load_ancient_eurasia_model()
        present_ne = next(
            float(population.initial_size)
            for population in model.model.populations
            if population.name == INTROGRESSION_TARGET_POPULATION
        )
        return model.model, INTROGRESSION_TARGET_POPULATION, present_ne
    raise ValueError(f"unknown focused demography {record['demography_id']!r}")


def _genotype_counts_from_samples(
    sample_nodes: np.ndarray, descendants: set[int]
) -> tuple[int, ...]:
    reshaped = np.asarray(sample_nodes, dtype=np.int64).reshape(-1, 2)
    return tuple(
        int(int(left) in descendants) + int(int(right) in descendants)
        for left, right in reshaped
    )


def _candidate_pool_alt_count_bounds(
    record: Mapping[str, Any], pool_diploids: int
) -> tuple[int, int]:
    n_diploids = int(pool_diploids)
    if n_diploids < 1:
        raise ValueError("candidate pool must contain at least one diploid")
    lower = float(record["population_af_lower"])
    upper = float(record["population_af_upper"])
    if not (
        math.isfinite(lower) and math.isfinite(upper) and 0.0 <= lower <= upper <= 1.0
    ):
        raise ValueError("candidate-pool AF interval must be finite within [0, 1]")
    haplotypes = 2 * n_diploids
    # The execution table stores decimal AF boundaries as binary floats. The
    # tolerance prevents an exact integer boundary such as 0.075 * 1000 from
    # becoming 75.00000000000001 and incorrectly excluding 75 copies.
    tolerance = 1e-9
    minimum = math.ceil(haplotypes * lower - tolerance)
    maximum = math.floor(haplotypes * upper + tolerance)
    if minimum > maximum:
        raise ValueError("candidate-pool AF interval contains no integer allele count")
    return minimum, maximum


def candidate_pool_frequency_gate(
    unit: Mapping[str, Any], genotype_counts: Sequence[int]
) -> dict[str, Any]:
    """Evaluate the shared, inclusive AF gate on a diploid candidate pool.

    This is the outer ascertainment used for both SLiM-selected trajectories
    and msprime neutral focal branches. Exact-k sampling of 100 diploids is a
    separate second stage and cannot rescue a candidate pool outside this band.
    """

    record = _validate_unit(unit)
    raw_values = np.asarray(genotype_counts)
    if raw_values.ndim != 1 or len(raw_values) < 1:
        raise ValueError("candidate-pool genotypes must be a nonempty vector")
    if not np.issubdtype(raw_values.dtype, np.number) or not np.all(
        np.isfinite(raw_values)
    ):
        raise ValueError("candidate-pool diploid allele counts must be finite numbers")
    if not np.all(raw_values == np.floor(raw_values)):
        raise ValueError("candidate-pool diploid allele counts must be integers")
    values = raw_values.astype(np.int64, copy=False)
    if np.any((values < 0) | (values > 2)):
        raise ValueError("candidate-pool diploid allele counts must be 0, 1, or 2")
    minimum, maximum = _candidate_pool_alt_count_bounds(record, len(values))
    observed = int(values.sum())
    passed = bool(minimum <= observed <= maximum)
    available = {value: int(np.count_nonzero(values == value)) for value in (0, 1, 2)}
    return {
        "schema": CANDIDATE_POOL_FREQUENCY_GATE_SCHEMA,
        "status": "pass" if passed else "fail",
        "passed": passed,
        "candidate_pool_diploids": len(values),
        "candidate_pool_haplotypes": 2 * len(values),
        "population_af_lower_inclusive": float(record["population_af_lower"]),
        "population_af_upper_inclusive": float(record["population_af_upper"]),
        "minimum_alt_count_inclusive": minimum,
        "maximum_alt_count_inclusive": maximum,
        "observed_alt_count": observed,
        "observed_af": float(observed / (2 * len(values))),
        "candidate_pool_genotype_counts": {
            "hom_ref": available[0],
            "heterozygous": available[1],
            "hom_alt": available[2],
        },
        "gate_order": "before_exact_sample_panel",
    }


def _validate_candidate_pool_frequency_gate_record(
    gate: Mapping[str, Any],
    unit: Mapping[str, Any],
    *,
    expected_pool_diploids: int,
) -> dict[str, Any]:
    record = _validate_unit(unit)
    pool_diploids = int(expected_pool_diploids)
    if gate.get("schema") != CANDIDATE_POOL_FREQUENCY_GATE_SCHEMA:
        raise ValueError("candidate-pool frequency-gate provenance is incompatible")
    if gate.get("passed") is not True or gate.get("status") != "pass":
        raise ValueError("accepted simulation did not pass its candidate-pool AF gate")
    if int(gate.get("candidate_pool_diploids", -1)) != pool_diploids:
        raise ValueError("candidate-pool frequency-gate size is incompatible")
    minimum, maximum = _candidate_pool_alt_count_bounds(record, pool_diploids)
    expected_scalars = {
        "candidate_pool_haplotypes": 2 * pool_diploids,
        "minimum_alt_count_inclusive": minimum,
        "maximum_alt_count_inclusive": maximum,
    }
    if any(int(gate.get(key, -1)) != value for key, value in expected_scalars.items()):
        raise ValueError("candidate-pool frequency-gate integer bounds are invalid")
    if not math.isclose(
        float(gate.get("population_af_lower_inclusive", math.nan)),
        float(record["population_af_lower"]),
        rel_tol=0.0,
        abs_tol=1e-12,
    ) or not math.isclose(
        float(gate.get("population_af_upper_inclusive", math.nan)),
        float(record["population_af_upper"]),
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        raise ValueError("candidate-pool frequency-gate AF bounds are invalid")
    observed = int(gate.get("observed_alt_count", -1))
    if not minimum <= observed <= maximum:
        raise ValueError(
            "candidate-pool frequency-gate observed count is outside bounds"
        )
    if not math.isclose(
        float(gate.get("observed_af", math.nan)),
        observed / (2 * pool_diploids),
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        raise ValueError("candidate-pool frequency-gate observed AF is inconsistent")
    classes = gate.get("candidate_pool_genotype_counts")
    if not isinstance(classes, Mapping):
        raise TypeError("candidate-pool frequency-gate genotype counts are absent")
    hom_ref = int(classes.get("hom_ref", -1))
    heterozygous = int(classes.get("heterozygous", -1))
    hom_alt = int(classes.get("hom_alt", -1))
    if hom_ref + heterozygous + hom_alt != pool_diploids or (
        heterozygous + 2 * hom_alt != observed
    ):
        raise ValueError(
            "candidate-pool frequency-gate genotype counts are inconsistent"
        )
    if gate.get("gate_order") != "before_exact_sample_panel":
        raise ValueError("candidate-pool frequency gate was applied in the wrong order")
    return dict(gate)


def _validate_exact_sample_panel_record(
    panel: Mapping[str, Any],
    unit: Mapping[str, Any],
    pool_gate: Mapping[str, Any],
) -> dict[str, Any]:
    record = _validate_unit(unit)
    if int(panel.get("candidate_pool_diploids", -1)) != int(
        pool_gate["candidate_pool_diploids"]
    ):
        raise ValueError("exact-panel candidate-pool size is inconsistent")
    if panel.get("candidate_pool_genotype_counts") != pool_gate.get(
        "candidate_pool_genotype_counts"
    ):
        raise ValueError("exact-panel candidate-pool genotype counts are inconsistent")
    if int(panel.get("eligible_genotype_count_triples", 0)) < 1:
        raise ValueError("accepted exact-panel provenance has no eligible triple")
    chosen = panel.get("chosen_sample_genotype_counts")
    if not isinstance(chosen, Mapping):
        raise TypeError("accepted exact-panel genotype counts are absent")
    n_ref = int(chosen.get("hom_ref", -1))
    n_het = int(chosen.get("heterozygous", -1))
    n_alt = int(chosen.get("hom_alt", -1))
    if (
        n_ref + n_het + n_alt != int(record["sample_diploids"])
        or n_het + 2 * n_alt != int(record["exact_sample_alt_count"])
        or n_ref < 2
        or n_het < 1
        or n_alt < 2
    ):
        raise ValueError("accepted exact-panel genotype counts violate AF/QC")
    if panel.get("selection_method") != (
        "uniform_over_all_size_N_exact_alt_count_subsets_passing_genotype_qc"
    ):
        raise ValueError("accepted exact-panel selection method is incompatible")
    if int(panel.get("panel_seed", 0)) < 1:
        raise ValueError("accepted exact-panel seed is invalid")
    return dict(panel)


def _passes_genotype_qc(counts: Sequence[int], exact_alt_count: int) -> bool:
    values = np.asarray(counts, dtype=np.int8)
    return bool(
        int(values.sum()) == int(exact_alt_count)
        and np.count_nonzero(values == 0) >= 2
        and np.count_nonzero(values == 1) >= 1
        and np.count_nonzero(values == 2) >= 2
    )


def feasible_exact_panel_genotype_triples(
    genotype_counts: Sequence[int],
    *,
    sample_diploids: int,
    exact_alt_count: int,
) -> list[dict[str, Any]]:
    """Enumerate exact-k, genotype-QC-feasible subsets of a diploid pool."""

    raw_values = np.asarray(genotype_counts)
    n_sample = int(sample_diploids)
    exact = int(exact_alt_count)
    if raw_values.ndim != 1 or len(raw_values) < 1:
        raise ValueError("candidate-pool genotypes must be a nonempty vector")
    if not np.issubdtype(raw_values.dtype, np.number) or not np.all(
        np.isfinite(raw_values)
    ):
        raise ValueError("candidate-pool diploid allele counts must be finite numbers")
    if not np.all(raw_values == np.floor(raw_values)):
        raise ValueError("candidate-pool diploid allele counts must be integers")
    values = raw_values.astype(np.int64, copy=False)
    if np.any((values < 0) | (values > 2)):
        raise ValueError("candidate-pool diploid allele counts must be 0, 1, or 2")
    if n_sample < 1 or n_sample > len(values):
        raise ValueError("sample diploid count must be within the candidate pool")
    if exact < 0 or exact > 2 * n_sample:
        raise ValueError("exact sample alternate count is outside the sample range")
    available = {value: int(np.count_nonzero(values == value)) for value in (0, 1, 2)}
    triples: list[dict[str, Any]] = []
    for n_alt in range(2, min(available[2], n_sample) + 1):
        n_het = int(exact - 2 * n_alt)
        n_ref = int(n_sample - n_alt - n_het)
        if n_het < 1 or n_ref < 2 or n_het > available[1] or n_ref > available[0]:
            continue
        triples.append(
            {
                "hom_ref": n_ref,
                "heterozygous": n_het,
                "hom_alt": n_alt,
                "log_number_of_subsets": (
                    _log_choose(available[0], n_ref)
                    + _log_choose(available[1], n_het)
                    + _log_choose(available[2], n_alt)
                ),
            }
        )
    return triples


def _introgression_migrations_by_tree(
    ts: tskit.TreeSequence, demography: msprime.Demography
) -> dict[int, list[tskit.Migration]]:
    result: dict[int, list[tskit.Migration]] = {}
    breakpoints = np.fromiter(ts.breakpoints(), dtype=float)
    for migration in ts.migrations():
        if not math.isclose(
            float(migration.time), INTROGRESSION_PULSE_GENERATIONS, abs_tol=1e-9
        ):
            continue
        if (
            _population_name(demography, migration.source) != "Loschbour"
            or _population_name(demography, migration.dest) != ARCHAIC_POPULATION
        ):
            continue
        first = max(
            0, int(np.searchsorted(breakpoints, migration.left, side="right") - 1)
        )
        last = min(
            ts.num_trees - 1,
            int(np.searchsorted(breakpoints, migration.right, side="left") - 1),
        )
        for tree_index in range(first, last + 1):
            left = max(float(migration.left), float(breakpoints[tree_index]))
            right = min(float(migration.right), float(breakpoints[tree_index + 1]))
            if right > left:
                result.setdefault(tree_index, []).append(migration)
    return result


def find_neutral_candidates(
    ts: tskit.TreeSequence,
    unit: Mapping[str, Any],
    demography: msprime.Demography,
    *,
    candidate_pool_diploids: int = DEFAULT_NEUTRAL_POOL_DIPLOIDS,
) -> list[NeutralCandidate]:
    """Enumerate pool-band and exact-panel-feasible neutral branches.

    The two demographies intentionally use different spatial ascertainment.
    EAS is conditioned at the requested 5-Mb site of an already 10-Mb
    ancestry, with eligible branches weighted only by mutation-time length.
    Ancient Eurasia is first simulated over 20 Mb because high-frequency Han
    archaic branches are rare at one point. Its search is driven by the pulse
    migration table, rather than scanning every node in every marginal tree,
    and eligible fixed-age branches are weighted by genomic span.
    """

    record = _validate_unit(unit)
    exact_count = int(record["exact_sample_alt_count"])
    target_population = str(record["population"])
    target_population_id = next(
        index
        for index, population in enumerate(demography.populations)
        if population.name == target_population
    )
    sample_nodes = np.asarray(
        ts.samples(population=target_population_id), dtype=np.int64
    )
    pool_diploids = int(candidate_pool_diploids)
    if pool_diploids < int(record["sample_diploids"]):
        raise ValueError("neutral candidate pool is smaller than the sample panel")
    if len(sample_nodes) != 2 * pool_diploids:
        raise ValueError("neutral ancestry sample count does not match candidate pool")
    minimum_pool_alt_count, maximum_pool_alt_count = _candidate_pool_alt_count_bounds(
        record, pool_diploids
    )

    def candidate_details(
        tree: tskit.Tree, node: int
    ) -> tuple[dict[str, Any], int] | None:
        branch_sample_count = int(tree.num_samples(node))
        if not (
            minimum_pool_alt_count <= branch_sample_count <= maximum_pool_alt_count
        ):
            return None
        descendants = {int(value) for value in tree.samples(node)}
        genotype_counts = _genotype_counts_from_samples(sample_nodes, descendants)
        pool_gate = candidate_pool_frequency_gate(record, genotype_counts)
        if not pool_gate["passed"]:
            return None
        triples = feasible_exact_panel_genotype_triples(
            genotype_counts,
            sample_diploids=int(record["sample_diploids"]),
            exact_alt_count=exact_count,
        )
        if not triples:
            return None
        return pool_gate, len(triples)

    if record["demography_kind"] != "introgression":
        if float(ts.sequence_length) < float(SEQUENCE_LENGTH_BP):
            raise ValueError("EAS neutral ancestry is shorter than the requested 10 Mb")
        tree = ts.at(float(FOCAL_POSITION_BP))
        candidates: list[NeutralCandidate] = []
        for raw_node in tree.nodes():
            node = int(raw_node)
            parent = int(tree.parent(node))
            if parent == tskit.NULL:
                continue
            details = candidate_details(tree, node)
            if details is None:
                continue
            pool_gate, n_triples = details
            lower = float(ts.node(node).time)
            upper = float(ts.node(parent).time)
            if upper <= lower:
                continue
            # The one-base interval forces the constructed mutation to the
            # fixed 5-Mb site after floor-alignment. Its span is deliberately
            # absent from the probability measure; at one site only branch
            # time determines the conditional neutral mutation probability.
            candidates.append(
                NeutralCandidate(
                    left=float(FOCAL_POSITION_BP),
                    right=float(FOCAL_POSITION_BP + 1),
                    node=node,
                    mutation_time_lower=lower,
                    mutation_time_upper=upper,
                    mutation_origin_population=target_population,
                    candidate_pool_alt_count=int(pool_gate["observed_alt_count"]),
                    candidate_pool_genotype_counts=(
                        int(pool_gate["candidate_pool_genotype_counts"]["hom_ref"]),
                        int(
                            pool_gate["candidate_pool_genotype_counts"]["heterozygous"]
                        ),
                        int(pool_gate["candidate_pool_genotype_counts"]["hom_alt"]),
                    ),
                    eligible_genotype_count_triples=n_triples,
                    n_introgressed_descendant_lineages=0,
                    weight=upper - lower,
                )
            )
        return candidates

    if not math.isclose(
        float(ts.sequence_length),
        float(INTROGRESSION_ASCERTAINMENT_LENGTH_BP),
        rel_tol=0.0,
        abs_tol=1e-9,
    ):
        raise ValueError("Han neutral ascertainment ancestry must be exactly 20 Mb")
    pulse_migrations = _introgression_migrations_by_tree(ts, demography)
    candidates = []
    for tree_index, migrations in pulse_migrations.items():
        tree = ts.at_index(int(tree_index))
        tree_left = max(
            float(tree.interval.left), float(ASCERTAINMENT_INTERIOR_LEFT_BP)
        )
        tree_right = min(
            float(tree.interval.right), float(ASCERTAINMENT_INTERIOR_RIGHT_BP)
        )
        if tree_right <= tree_left:
            continue
        routes_by_origin_node: dict[int, list[tuple[float, float, int]]] = {}
        for migration in migrations:
            left = max(tree_left, float(migration.left))
            right = min(tree_right, float(migration.right))
            if right <= left:
                continue
            node = int(migration.node)
            parent = int(tree.parent(node))
            while (
                parent != tskit.NULL
                and float(ts.node(parent).time)
                <= INTROGRESSION_MUTATION_AGE_GENERATIONS
            ):
                node = parent
                parent = int(tree.parent(node))
            if parent == tskit.NULL or not (
                float(ts.node(node).time)
                <= INTROGRESSION_MUTATION_AGE_GENERATIONS
                < float(ts.node(parent).time)
            ):
                continue
            routes_by_origin_node.setdefault(node, []).append(
                (left, right, int(migration.node))
            )

        for node, routes in routes_by_origin_node.items():
            details = candidate_details(tree, node)
            if details is None:
                continue
            pool_gate, n_triples = details
            # Split the union of migration intervals into nonoverlapping pieces
            # so shared ancestral branches are not double-weighted when several
            # pulse lineages descend from the same 2,400-generation branch.
            endpoints = sorted(
                {value for left, right, _ in routes for value in (left, right)}
            )
            for left, right in pairwise(endpoints):
                if right <= left:
                    continue
                midpoint = (left + right) / 2
                lineage_nodes = {
                    migration_node
                    for route_left, route_right, migration_node in routes
                    if route_left <= midpoint < route_right
                }
                if not lineage_nodes:
                    continue
                candidates.append(
                    NeutralCandidate(
                        left=left,
                        right=right,
                        node=node,
                        mutation_time_lower=(INTROGRESSION_MUTATION_AGE_GENERATIONS),
                        mutation_time_upper=(INTROGRESSION_MUTATION_AGE_GENERATIONS),
                        mutation_origin_population=ARCHAIC_POPULATION,
                        candidate_pool_alt_count=int(pool_gate["observed_alt_count"]),
                        candidate_pool_genotype_counts=(
                            int(pool_gate["candidate_pool_genotype_counts"]["hom_ref"]),
                            int(
                                pool_gate["candidate_pool_genotype_counts"][
                                    "heterozygous"
                                ]
                            ),
                            int(pool_gate["candidate_pool_genotype_counts"]["hom_alt"]),
                        ),
                        eligible_genotype_count_triples=n_triples,
                        n_introgressed_descendant_lineages=len(lineage_nodes),
                        weight=right - left,
                    )
                )
    return candidates


def _crop_and_add_mutations(
    ancestry: tskit.TreeSequence,
    *,
    candidate: NeutralCandidate,
    position: float,
    mutation_time: float,
    mutation_seed: int,
) -> tskit.TreeSequence:
    # Integer-align the focal draw before cropping. Subtracting arbitrary
    # floating coordinates can otherwise leave a 10,000,000.000000002 span,
    # which msprime correctly rejects as a mutation-map length mismatch.
    position = float(math.floor(position))
    window_left = float(position - FOCAL_POSITION_BP)
    window_right = float(window_left + SEQUENCE_LENGTH_BP)
    cropped = ancestry.keep_intervals(
        [[window_left, window_right]], simplify=False
    ).trim()
    if float(cropped.sequence_length) != float(SEQUENCE_LENGTH_BP):
        raise ValueError("ascertainment crop did not produce exactly 10 Mb")
    focal = float(position - window_left)
    rate_map = msprime.RateMap(
        position=[0.0, focal, focal + 1.0, float(SEQUENCE_LENGTH_BP)],
        rate=[MUTATION_RATE, 0.0, MUTATION_RATE],
    )
    mutated = msprime.sim_mutations(
        cropped,
        rate=rate_map,
        model=msprime.JC69(),
        random_seed=mutation_seed,
        discrete_genome=True,
    )
    tree = mutated.at(focal)
    # simplify=False preserves node IDs through interval cropping.
    if candidate.node not in set(tree.nodes()):
        raise ValueError("chosen neutral focal node was lost during cropping")
    tables = mutated.dump_tables()
    site_id = tables.sites.add_row(position=focal, ancestral_state="A")
    tables.mutations.add_row(
        site=site_id,
        node=int(candidate.node),
        derived_state="G",
        time=float(mutation_time),
    )
    tables.sort()
    tables.build_index()
    tables.compute_mutation_parents()
    return tables.tree_sequence()


def _ensure_sample_individuals(ts: tskit.TreeSequence) -> tskit.TreeSequence:
    """Attach consecutive sample pairs to diploid individuals when absent.

    msprime's mapping-form samples currently carry individuals, but retaining
    this explicit repair makes the decoder pairing contract independent of that
    convenience behavior and fails closed on partially assigned samples.
    """

    samples = np.asarray(ts.samples(), dtype=np.int32)
    individual_ids = ts.tables.nodes.individual[samples]
    if np.all(individual_ids != tskit.NULL):
        return ts
    if np.any(individual_ids != tskit.NULL) or len(samples) % 2:
        raise ValueError("sample nodes have incomplete diploid individual assignment")
    tables = ts.dump_tables()
    nodes = tables.nodes
    individual_column = np.array(nodes.individual, copy=True)
    for left, right in samples.reshape(-1, 2):
        individual_id = tables.individuals.add_row()
        individual_column[int(left)] = individual_id
        individual_column[int(right)] = individual_id
    nodes.set_columns(
        flags=nodes.flags,
        time=nodes.time,
        population=nodes.population,
        individual=individual_column,
        metadata=nodes.metadata,
        metadata_offset=nodes.metadata_offset,
    )
    return tables.tree_sequence()


def _neutral_attempt_seed(unit_seed: int, attempt: int, component: str) -> int:
    return _stable_seed(unit_seed, f"neutral:{attempt}:{component}")


def _neutral_ancestry_models() -> list[msprime.AncestryModel]:
    """Return a fresh recent-DTWF then SMC-prime ancestry model chain."""

    return [
        msprime.DiscreteTimeWrightFisher(
            duration=NEUTRAL_RECENT_DTWF_DURATION_GENERATIONS
        ),
        msprime.SmcPrimeApproxCoalescent(),
    ]


def simulate_neutral_unit(
    unit: Mapping[str, Any],
    repo_root: str | Path,
    unit_dir: str | Path,
    *,
    max_attempts: int = DEFAULT_NEUTRAL_MAX_ATTEMPTS,
    ancestry_batch_size: int = DEFAULT_NEUTRAL_ANCESTRY_BATCH_SIZE,
    candidate_pool_diploids: int = DEFAULT_NEUTRAL_POOL_DIPLOIDS,
) -> tuple[tskit.TreeSequence, pd.DataFrame, dict[str, Any], list[dict[str, Any]]]:
    record = _validate_unit(unit)
    if record["simulation_class"] != "neutral":
        raise ValueError("simulate_neutral_unit received a selected unit")
    pool_diploids = int(candidate_pool_diploids)
    if pool_diploids < int(record["sample_diploids"]):
        raise ValueError("neutral candidate pool is smaller than the sample panel")
    minimum_pool_alt_count, maximum_pool_alt_count = _candidate_pool_alt_count_bounds(
        record, pool_diploids
    )
    demography, population, _ = load_cell_demography(record, repo_root)
    introgression = record["demography_kind"] == "introgression"
    ancestry_length_bp = (
        INTROGRESSION_ASCERTAINMENT_LENGTH_BP if introgression else SEQUENCE_LENGTH_BP
    )
    directory = Path(unit_dir)
    directory.mkdir(parents=True, exist_ok=True)
    ledger_path = directory / "attempts.tsv"
    attempts: list[dict[str, Any]] = []
    if ledger_path.is_file():
        attempts = pd.read_csv(ledger_path, sep="\t").to_dict(orient="records")
        observed = [int(row["attempt_zero_based"]) for row in attempts]
        if observed != list(range(len(attempts))):
            raise ValueError("neutral attempt ledger is not contiguous from zero")
        # A killed process may leave a partial fixed-size ancestry batch. Redo
        # that complete deterministic batch so the candidate weighting stays
        # invariant across restarts.
        keep = (len(attempts) // int(ancestry_batch_size)) * int(ancestry_batch_size)
        if any(bool(row.get("accepted", False)) for row in attempts):
            accepted_index = next(
                index
                for index, row in enumerate(attempts)
                if bool(row.get("accepted", False))
            )
            keep = (accepted_index // int(ancestry_batch_size)) * int(
                ancestry_batch_size
            )
        attempts = attempts[:keep]
        _atomic_frame(ledger_path, pd.DataFrame(attempts))
    if int(max_attempts) < 1 or int(ancestry_batch_size) < 1:
        raise ValueError("neutral search and ancestry-batch sizes must be positive")
    for batch_start in range(
        len(attempts), int(max_attempts), int(ancestry_batch_size)
    ):
        batch_stop = min(int(max_attempts), batch_start + int(ancestry_batch_size))
        entries: list[tuple[tskit.TreeSequence, NeutralCandidate, int]] = []
        for attempt in range(batch_start, batch_stop):
            ancestry_seed = _neutral_attempt_seed(
                int(record["seed"]), attempt, "ancestry"
            )
            mutation_seed = _neutral_attempt_seed(
                int(record["seed"]), attempt, "mutations"
            )
            started = perf_counter()
            ancestry = msprime.sim_ancestry(
                samples={population: pool_diploids},
                ploidy=2,
                demography=demography,
                sequence_length=ancestry_length_bp,
                recombination_rate=RECOMBINATION_RATE,
                random_seed=ancestry_seed,
                model=_neutral_ancestry_models(),
                record_migrations=introgression,
            )
            candidates = find_neutral_candidates(
                ancestry,
                record,
                demography,
                candidate_pool_diploids=pool_diploids,
            )
            attempts.append(
                {
                    "attempt_zero_based": attempt,
                    "ancestry_batch_zero_based": (
                        batch_start // int(ancestry_batch_size)
                    ),
                    "ancestry_seed": ancestry_seed,
                    "mutation_seed": mutation_seed,
                    "candidate_pool_diploids": pool_diploids,
                    "candidate_pool_minimum_alt_count": minimum_pool_alt_count,
                    "candidate_pool_maximum_alt_count": maximum_pool_alt_count,
                    "ancestry_model": NEUTRAL_ANCESTRY_MODEL_NAME,
                    "recent_dtwf_duration_generations": (
                        NEUTRAL_RECENT_DTWF_DURATION_GENERATIONS
                    ),
                    "ancestry_length_bp": ancestry_length_bp,
                    "n_eligible_branch_rectangles": len(candidates),
                    "total_eligible_branch_measure": float(
                        sum(candidate.weight for candidate in candidates)
                    ),
                    "n_trees": ancestry.num_trees,
                    "n_migrations": ancestry.num_migrations,
                    "elapsed_seconds": perf_counter() - started,
                    "eligible": bool(candidates),
                    "accepted": False,
                }
            )
            entries.extend((ancestry, candidate, attempt) for candidate in candidates)
            _atomic_frame(ledger_path, pd.DataFrame(attempts))
        if not entries:
            continue
        draw_seed = _stable_seed(
            int(record["seed"]),
            f"neutral-batch:{batch_start // int(ancestry_batch_size)}:candidate",
        )
        weights = np.asarray([entry[1].weight for entry in entries], dtype=float)
        if not np.all(np.isfinite(weights)) or np.any(weights <= 0):
            raise ValueError("neutral candidate weights must be finite and positive")
        rng = np.random.default_rng(draw_seed)
        entry_index = int(rng.choice(len(entries), p=weights / weights.sum()))
        ancestry, candidate, accepted_attempt = entries[entry_index]
        position = (
            float(rng.uniform(candidate.left, candidate.right))
            if introgression
            else float(FOCAL_POSITION_BP)
        )
        mutation_time = float(
            rng.uniform(candidate.mutation_time_lower, candidate.mutation_time_upper)
        )
        mutation_seed = _neutral_attempt_seed(
            int(record["seed"]), accepted_attempt, "mutations"
        )
        pool_ts = _crop_and_add_mutations(
            ancestry,
            candidate=candidate,
            position=position,
            mutation_time=mutation_time,
            mutation_seed=mutation_seed,
        )
        pool_ts = _ensure_sample_individuals(pool_ts)
        pool_table = build_diploid_pair_table(pool_ts, FOCAL_POSITION_BP)
        pool_counts = pool_table["focal_selected_allele_count"].to_numpy(dtype=int)
        if len(pool_counts) != pool_diploids:
            raise ValueError("constructed neutral candidate pool has the wrong size")
        pool_gate = candidate_pool_frequency_gate(record, pool_counts)
        if not pool_gate["passed"]:
            raise ValueError(
                "constructed neutral candidate pool failed its prespecified AF gate"
            )
        observed_pool_genotype_counts = tuple(
            int(np.count_nonzero(pool_counts == value)) for value in (0, 1, 2)
        )
        if int(pool_counts.sum()) != candidate.candidate_pool_alt_count or (
            observed_pool_genotype_counts != candidate.candidate_pool_genotype_counts
        ):
            raise ValueError(
                "constructed neutral candidate-pool genotypes changed during cropping"
            )
        panel_seed = _neutral_attempt_seed(
            int(record["seed"]), accepted_attempt, "panel"
        )
        ts, panel_record = select_exact_sample_panel(
            pool_ts,
            sample_diploids=int(record["sample_diploids"]),
            exact_alt_count=int(record["exact_sample_alt_count"]),
            seed=panel_seed,
        )
        pair_table = build_diploid_pair_table(ts, FOCAL_POSITION_BP)
        counts = pair_table["focal_selected_allele_count"].to_numpy(dtype=int)
        if not _passes_genotype_qc(counts, int(record["exact_sample_alt_count"])):
            raise ValueError(
                "constructed neutral focal variant failed exact AF/genotype QC"
            )
        attempts[accepted_attempt].update(
            {
                "accepted": True,
                "candidate_seed": draw_seed,
                "panel_seed": panel_seed,
                "candidate_pool_alt_count": int(pool_gate["observed_alt_count"]),
                "candidate_pool_af": float(pool_gate["observed_af"]),
                "eligible_genotype_count_triples": int(
                    candidate.eligible_genotype_count_triples
                ),
                "sample_alt_count": int(counts.sum()),
                "sample_af": float(counts.mean() / 2),
                "n_hom_ref": int(np.count_nonzero(counts == 0)),
                "n_heterozygous": int(np.count_nonzero(counts == 1)),
                "n_hom_alt": int(np.count_nonzero(counts == 2)),
            }
        )
        _atomic_frame(ledger_path, pd.DataFrame(attempts))
        provenance = {
            "engine": "msprime",
            "msprime_version": msprime.__version__,
            "model": NEUTRAL_ANCESTRY_MODEL_NAME,
            "model_keywords": list(NEUTRAL_ANCESTRY_MODEL_KEYWORDS),
            "recent_dtwf_duration_generations": (
                NEUTRAL_RECENT_DTWF_DURATION_GENERATIONS
            ),
            "coalescent_approximation": "SMC-prime after 200 generations",
            "ascertainment_length_bp": ancestry_length_bp,
            "window_length_bp": SEQUENCE_LENGTH_BP,
            "candidate_weight": (
                "genomic_span_bp_at_fixed_2400_generation_origin"
                if introgression
                else "eligible_mutation_time_generations_at_fixed_5Mb_site"
            ),
            "candidate_index": (
                "Loschbour_to_Neanderthal_pulse_migrations_then_ascend_to_"
                "branch_crossing_2400_generations"
                if introgression
                else "marginal_tree_at_fixed_5Mb_site"
            ),
            "ancestry_batch_size": int(ancestry_batch_size),
            "candidate_pool_diploids": pool_diploids,
            "candidate_pool_frequency_gate": pool_gate,
            "sample_panel": panel_record,
            "matched_selected_neutral_ascertainment": (
                "same_candidate_pool_size_and_AF_band_then_same_exact_k_panel_sampler"
            ),
            "accepted_attempt_zero_based": accepted_attempt,
            "accepted_batch_zero_based": (batch_start // int(ancestry_batch_size)),
            "ancestry_seed": int(attempts[accepted_attempt]["ancestry_seed"]),
            "candidate_seed": draw_seed,
            "mutation_seed": mutation_seed,
            "panel_seed": panel_seed,
            "pre_crop_focal_position": position,
            "mutation_time_generations_ago": mutation_time,
            "mutation_origin_population": candidate.mutation_origin_population,
            "archaic_specific_no_ils": introgression,
            "entry_route": (
                "Loschbour_to_Neanderthal_backward_mass_migration_at_2272_generations"
                if introgression
                else "single_population_EAS_lineage"
            ),
            "n_eligible_branch_rectangles_in_batch": len(entries),
            "neutral_origin_matching": (
                "archaic_population_branch_with_one_or_more_pulse_lineages; "
                "nonzero standing variation at pulse"
                if introgression
                else "single EAS branch"
            ),
            "n_introgressed_descendant_lineages": (
                candidate.n_introgressed_descendant_lineages
            ),
        }
        return ts, pair_table, provenance, attempts
    raise SimulationExhausted(
        f"{record['unit_id']} found no pool-band, panel-feasible neutral branch "
        f"in {max_attempts} "
        f"independent {ancestry_length_bp / 1_000_000:g}-Mb ancestries"
    )


def with_post_pulse_recipient_nonloss(
    sweep: Any,
    *,
    realized_han_split_generations: float,
) -> Any:
    """Condition the introgressed allele on survival along the Han recipient path."""

    boundary = float(realized_han_split_generations)
    if not 0.0 < boundary < INTROGRESSION_PULSE_GENERATIONS:
        raise ValueError("realized Han split lies outside the post-pulse interval")
    for event in sweep.extended_events:
        if (
            isinstance(event, stdpopsim.ConditionOnAlleleFrequency)
            and str(event.single_site_id) == FOCAL_SITE_ID
            and str(event.population) in {"Loschbour", INTROGRESSION_TARGET_POPULATION}
            and str(event.op) == ">"
            and math.isclose(float(event.allele_frequency), 0.0, abs_tol=1e-12)
        ):
            raise ValueError("recipient nonloss conditioning is already present")
    recipient_conditions = (
        stdpopsim.ConditionOnAlleleFrequency(
            start_time=stdpopsim.GenerationAfter(INTROGRESSION_PULSE_GENERATIONS),
            end_time=boundary,
            single_site_id=FOCAL_SITE_ID,
            population="Loschbour",
            op=">",
            allele_frequency=0.0,
        ),
        stdpopsim.ConditionOnAlleleFrequency(
            start_time=stdpopsim.GenerationAfter(boundary),
            end_time=0.0,
            single_site_id=FOCAL_SITE_ID,
            population=INTROGRESSION_TARGET_POPULATION,
            op=">",
            allele_frequency=0.0,
        ),
    )
    events = list(sweep.extended_events)
    insertion = next(
        (
            index
            for index, event in enumerate(events)
            if isinstance(event, stdpopsim.ChangeMutationFitness)
        ),
        len(events),
    )
    events[insertion:insertion] = recipient_conditions
    record = sweep.provenance_record()
    record["post_pulse_recipient_nonloss"] = {
        "mode": "piecewise_loschbour_then_han_nonloss",
        "conditioned_on_introgressed_copy_reaching_recipient": True,
        "population_intervals": [
            {
                "population": "Loschbour",
                "start_time": INTROGRESSION_PULSE_GENERATIONS,
                "start_time_semantics": "generation_after_in_forward_time",
                "end_time": boundary,
            },
            {
                "population": INTROGRESSION_TARGET_POPULATION,
                "start_time": boundary,
                "start_time_semantics": "generation_after_in_forward_time",
                "end_time": 0.0,
            },
        ],
    }
    record["extended_events"] = serialize_extended_events(events)
    return replace(sweep, extended_events=tuple(events), record=record)


def build_raw_han_selected_sweep(
    unit: Mapping[str, Any],
    *,
    slim_scaling_factor: float,
) -> Any:
    """Build the endpoint-unmodified Han sweep shared by calibration/production."""

    record = _validate_unit(unit)
    if (
        record["simulation_class"] != "selected"
        or record["demography_id"] != "ancient_eurasia_han_introgression"
    ):
        raise ValueError("raw Han sweep builder requires a Han selected unit")
    q = float(slim_scaling_factor)
    realized_split = math.floor(HAN_SPLIT_GENERATIONS / q) * q
    sweep = build_introgression_sweep_spec(
        selection_coefficient=float(record["selection_coefficient"]),
        minimum_han_frequency=float(record["population_af_lower"]),
        maximum_han_frequency=float(record["population_af_upper"]),
        mutation_age_generations=INTROGRESSION_MUTATION_AGE_GENERATIONS,
        source_minimum_frequency=INTROGRESSION_SOURCE_PRESENCE_EPSILON,
        realized_han_split_generations=realized_split,
    )
    return with_post_pulse_recipient_nonloss(
        sweep, realized_han_split_generations=realized_split
    )


def apply_han_selection_end(sweep: Any, endpoint: float) -> tuple[Any, dict[str, Any]]:
    """Apply one fixed endpoint across Loschbour/Han fitness intervals."""

    endpoint = float(endpoint)
    events: list[stdpopsim.ExtendedEvent] = []
    loschbour: stdpopsim.ChangeMutationFitness | None = None
    han: stdpopsim.ChangeMutationFitness | None = None
    for event in sweep.extended_events:
        if (
            isinstance(event, stdpopsim.ChangeMutationFitness)
            and str(event.single_site_id) == FOCAL_SITE_ID
        ):
            if str(event.population) == "Loschbour":
                loschbour = event
                continue
            if str(event.population) == INTROGRESSION_TARGET_POPULATION:
                han = event
                continue
        events.append(event)
    if loschbour is None or han is None:
        raise ValueError("expected one Loschbour and one Han fitness interval")
    boundary = float(han.start_time)
    coefficient = float(han.selection_coeff)
    if not 0.0 < endpoint < float(loschbour.start_time):
        raise ValueError("calibrated endpoint lies outside selected history")
    intervals: list[dict[str, Any]] = []
    if endpoint >= boundary:
        events.append(
            stdpopsim.ChangeMutationFitness(
                start_time=float(loschbour.start_time),
                end_time=endpoint,
                single_site_id=FOCAL_SITE_ID,
                population="Loschbour",
                selection_coeff=coefficient,
                dominance_coeff=float(loschbour.dominance_coeff),
            )
        )
        intervals.append(
            {
                "population": "Loschbour",
                "start_time": float(loschbour.start_time),
                "end_time": endpoint,
            }
        )
    else:
        events.append(loschbour)
        events.append(
            stdpopsim.ChangeMutationFitness(
                start_time=boundary,
                end_time=endpoint,
                single_site_id=FOCAL_SITE_ID,
                population=INTROGRESSION_TARGET_POPULATION,
                selection_coeff=coefficient,
                dominance_coeff=float(han.dominance_coeff),
            )
        )
        intervals.extend(
            (
                {
                    "population": "Loschbour",
                    "start_time": float(loschbour.start_time),
                    "end_time": float(loschbour.end_time),
                },
                {
                    "population": INTROGRESSION_TARGET_POPULATION,
                    "start_time": boundary,
                    "end_time": endpoint,
                },
            )
        )
    record = sweep.provenance_record()
    cessation = {
        "mode": "fixed_endpoint_from_disjoint_terminal_af_calibration",
        "selection_end_generations_ago": endpoint,
        "population_intervals": intervals,
        "frequency_after_cessation": "neutral_drift_to_present",
        "gamma_smc_statistics_used_for_calibration": False,
    }
    record["selection"]["mode"] = (
        "standing_variation_at_introgression_onset_frozen_finite_episode"
    )
    record["selection"]["cessation"] = cessation
    record["selection"]["han_interval_generations_ago"] = (
        [boundary, endpoint] if endpoint < boundary else None
    )
    record["extended_events"] = serialize_extended_events(events)
    return replace(sweep, extended_events=tuple(events), record=record), cessation


def with_external_sample_endpoint_ascertainment(
    sweep: Any,
    *,
    population: str,
    target_frequency: float,
    population_af_lower: float,
    population_af_upper: float,
    sample_diploids: int,
    exact_sample_alt_count: int,
    selected_pool_diploids: int,
) -> Any:
    """Replace terminal population-AF rejection with exact-panel ascertainment.

    The generic introgression sweep adds lower and upper Han frequency
    conditions at generation zero.  Those conditions make SLiM rejection-sample
    a narrow population-frequency band internally and can dominate runtime.
    The focused campaign instead applies the prespecified AF band to the
    completed 500-diploid candidate pool, then accepts only a uniformly drawn
    diploid panel with the exact alternate-copy count and genotype-class QC.

    Mutation origin, nonloss/presence conditions, and selected fitness intervals
    are retained unchanged. Exactly two terminal conditions must be removed.
    """

    target = float(target_frequency)
    lower = float(population_af_lower)
    upper = float(population_af_upper)
    n_diploids = int(sample_diploids)
    exact_count = int(exact_sample_alt_count)
    pool_diploids = int(selected_pool_diploids)
    if n_diploids <= 0 or pool_diploids < n_diploids:
        raise ValueError("selected panel and candidate-pool sizes are inconsistent")
    if not (0.0 <= lower <= target <= upper <= 1.0):
        raise ValueError("candidate-pool AF interval does not contain the target")
    if exact_count != round(2 * n_diploids * target):
        raise ValueError("exact sample alternate count is inconsistent with target AF")

    retained: list[stdpopsim.ExtendedEvent] = []
    removed: list[stdpopsim.ConditionOnAlleleFrequency] = []
    for event in sweep.extended_events:
        terminal_population_condition = (
            isinstance(event, stdpopsim.ConditionOnAlleleFrequency)
            and str(event.population) == str(population)
            and math.isclose(float(event.start_time), 0.0, abs_tol=1e-12)
            and math.isclose(float(event.end_time), 0.0, abs_tol=1e-12)
        )
        if terminal_population_condition:
            removed.append(event)
        else:
            retained.append(event)
    if len(removed) != 2 or {str(event.op) for event in removed} != {">=", "<="}:
        raise RuntimeError(
            "expected exactly the lower and upper terminal population AF conditions"
        )

    record = sweep.provenance_record()
    frequency_key = (
        "present_han_frequency"
        if str(population) == INTROGRESSION_TARGET_POPULATION
        else "present_frequency_interval"
    )
    requested = dict(record.get(frequency_key, {}))
    record[frequency_key] = {
        "population_conditioning_in_slim": False,
        "removed_terminal_condition_count": len(removed),
        "requested_population_interval_not_applied": requested,
        "endpoint_ascertainment": (
            "external_candidate_pool_af_band_then_exact_sample_panel"
        ),
    }
    record["focused_sample_endpoint_ascertainment"] = {
        "population": str(population),
        "population_frequency_conditioning_in_slim": False,
        "target_sample_allele_frequency": target,
        "sample_diploids": n_diploids,
        "exact_sample_alt_count": exact_count,
        "candidate_pool_diploids": pool_diploids,
        "candidate_pool_frequency_gate": {
            "schema": CANDIDATE_POOL_FREQUENCY_GATE_SCHEMA,
            "lower_inclusive": lower,
            "upper_inclusive": upper,
            "applied_after_completed_slim_trajectory": True,
            "applied_before_exact_sample_panel": True,
        },
        "minimum_hom_ref_diploids": 2,
        "minimum_heterozygous_diploids": 1,
        "minimum_hom_alt_diploids": 2,
        "panel_sampling": (
            "uniform_over_all_size_N_exact_alt_count_subsets_passing_genotype_qc"
        ),
        "inference_conditioning": (
            "selected results are conditional on the completed candidate-pool AF "
            "band, exact sampled AF, and genotype QC; no terminal population-AF "
            "rejection is performed inside SLiM"
        ),
    }
    record["extended_events"] = serialize_extended_events(retained)
    return replace(sweep, extended_events=tuple(retained), record=record)


def _selected_objects(
    unit: Mapping[str, Any],
    repo_root: str | Path,
    slim_scaling_factor: float,
    selected_pool_diploids: int,
    *,
    han_selection_calibration: Mapping[str, Any] | None = None,
    han_selection_end_generations_ago: float | None = None,
) -> tuple[stdpopsim.DemographicModel, Any, dict[str, int], dict[str, Any]]:
    record = _validate_unit(unit)
    target = float(record["target_allele_frequency"])
    lower = float(record["population_af_lower"])
    upper = float(record["population_af_upper"])
    coefficient = float(record["selection_coefficient"])
    if record["demography_id"] == "eas_phlash_median":
        artifact = load_phlash_eas_npz(
            Path(repo_root) / "no_introgression_EAS_sim" / "resources" / "EAS.npz",
            expected_sha256=EAS_RESOURCE_SHA256,
        )
        demographic = build_eas_demography_models(artifact)["median"]
        estimates = build_de_novo_origin_age_grid(
            {"q025": demographic, "median": demographic, "q975": demographic},
            selection_coefficients=(coefficient,),
            target_frequency=target,
            slim_scaling_factor=slim_scaling_factor,
        )
        estimate = next(item for item in estimates if item.trajectory_label == "median")
        sweep = build_no_introgression_sweep_spec(
            origin_age_generations=estimate.age_generations,
            selection_coefficient=coefficient,
            target_frequency=target,
            lower_frequency=lower,
            upper_frequency=upper,
        )
        sweep = with_external_sample_endpoint_ascertainment(
            sweep,
            population=EAS_POPULATION,
            target_frequency=target,
            population_af_lower=lower,
            population_af_upper=upper,
            sample_diploids=int(record["sample_diploids"]),
            exact_sample_alt_count=int(record["exact_sample_alt_count"]),
            selected_pool_diploids=int(selected_pool_diploids),
        )
        return (
            demographic.stdpopsim_model,
            sweep,
            {EAS_POPULATION: int(selected_pool_diploids)},
            {"sweep": sweep.provenance_record(), "origin_age": estimate.to_record()},
        )
    if record["demography_id"] == "ancient_eurasia_han_introgression":
        if (han_selection_calibration is None) == (
            han_selection_end_generations_ago is None
        ):
            raise ValueError(
                "Han selected construction requires exactly one frozen production "
                "calibration or explicit calibration-candidate endpoint"
            )
        sweep = build_raw_han_selected_sweep(
            record, slim_scaling_factor=slim_scaling_factor
        )
        if han_selection_calibration is not None:
            endpoint_binding = resolve_han_selection_endpoint(
                han_selection_calibration,
                selection_coefficient=coefficient,
                target_frequency=target,
            )
            endpoint = float(endpoint_binding["selection_end_generations_ago"])
            endpoint_source = "frozen_production_calibration"
        else:
            endpoint = float(han_selection_end_generations_ago)
            endpoint_binding = {
                "selection_coefficient": coefficient,
                "target_allele_frequency": target,
                "selection_end_generations_ago": endpoint,
            }
            endpoint_source = "explicit_calibration_candidate"
        sweep, cessation = apply_han_selection_end(sweep, endpoint)
        sweep = with_external_sample_endpoint_ascertainment(
            sweep,
            population=INTROGRESSION_TARGET_POPULATION,
            target_frequency=target,
            population_af_lower=lower,
            population_af_upper=upper,
            sample_diploids=int(record["sample_diploids"]),
            exact_sample_alt_count=int(record["exact_sample_alt_count"]),
            selected_pool_diploids=int(selected_pool_diploids),
        )
        return (
            sweep.demographic_model,
            sweep,
            {INTROGRESSION_TARGET_POPULATION: int(selected_pool_diploids)},
            {
                "sweep": sweep.provenance_record(),
                "selection_endpoint": {
                    "source": endpoint_source,
                    "binding": endpoint_binding,
                    "cessation": cessation,
                },
            },
        )
    raise ValueError(f"unknown selected demography {record['demography_id']!r}")


def _selected_focal_expectation(
    demographic_model: stdpopsim.DemographicModel,
    sweep: Any,
    slim_scaling_factor: float,
) -> dict[str, Any]:
    """Derive the exact selected-mutation identity expected in a SLiM tree."""

    draws = [
        event
        for event in sweep.extended_events
        if isinstance(event, stdpopsim.DrawMutation)
        and str(event.single_site_id) == FOCAL_SITE_ID
    ]
    if len(draws) != 1:
        raise ValueError("selected sweep must contain exactly one focal DrawMutation")
    draw = draws[0]
    source_name = str(draw.population)
    source_ids = [
        int(population.id)
        for population in demographic_model.populations
        if str(population.name) == source_name
    ]
    if len(source_ids) != 1:
        raise ValueError(
            f"selected mutation source population is not unique: {source_name}"
        )
    scaling = float(slim_scaling_factor)
    requested_time = float(draw.time)
    if not math.isfinite(scaling) or scaling <= 0:
        raise ValueError("SLiM scaling factor must be finite and positive")
    # stdpopsim 0.3 maps event times onto the nearest integer SLiM tick.
    realized_time = float(round(requested_time / scaling) * scaling)
    return {
        "schema": SELECTED_FOCAL_IDENTITY_SCHEMA,
        "single_site_id": FOCAL_SITE_ID,
        "focal_position_0based": int(FOCAL_POSITION_BP),
        "expected_mutation_type": 1,
        "expected_selection_coeff": 0.0,
        "source_population": source_name,
        "source_subpopulation": source_ids[0],
        "requested_origin_time_generations": requested_time,
        "realized_origin_time_generations": realized_time,
        "slim_scaling_factor": scaling,
        "slim_time_rule": "cycle - realized_origin_time / Q",
    }


def _focal_sites(ts: tskit.TreeSequence) -> list[tskit.Site]:
    return [
        site
        for site in ts.sites()
        if math.isclose(
            float(site.position),
            float(FOCAL_POSITION_BP),
            rel_tol=0.0,
            abs_tol=1e-9,
        )
    ]


def inspect_selected_focal_identity(
    ts: tskit.TreeSequence,
    expectation: Mapping[str, Any],
) -> dict[str, Any]:
    """Inspect selected-site identity, distinguishing stochastic loss from corruption."""

    if expectation.get("schema") != SELECTED_FOCAL_IDENTITY_SCHEMA:
        raise ValueError("selected focal identity expectation schema is incompatible")
    focal_sites = _focal_sites(ts)
    observed_entries: list[tuple[tskit.Mutation, Mapping[str, Any]]] = []
    malformed_metadata = False
    for site in focal_sites:
        for mutation in site.mutations:
            metadata = mutation.metadata
            if not isinstance(metadata, Mapping):
                malformed_metadata = True
                continue
            mutation_list = metadata.get("mutation_list")
            if not isinstance(mutation_list, list) or not mutation_list:
                malformed_metadata = True
                continue
            for entry in mutation_list:
                if not isinstance(entry, Mapping):
                    malformed_metadata = True
                    continue
                observed_entries.append((mutation, entry))
    observed_types = sorted(
        {
            int(entry["mutation_type"])
            for _, entry in observed_entries
            if "mutation_type" in entry
        }
    )
    expected_type = int(expectation["expected_mutation_type"])
    selected_entries = [
        (mutation, entry)
        for mutation, entry in observed_entries
        if int(entry.get("mutation_type", -1)) == expected_type
    ]
    base = {
        "schema": SELECTED_FOCAL_IDENTITY_SCHEMA,
        "focal_position_0based": int(FOCAL_POSITION_BP),
        "expected_mutation_type": expected_type,
        "expected_source_population": str(expectation["source_population"]),
        "expected_source_subpopulation": int(expectation["source_subpopulation"]),
        "expected_origin_time_generations": float(
            expectation["realized_origin_time_generations"]
        ),
        "observed_focal_site_count": len(focal_sites),
        "observed_tskit_mutation_count": sum(
            len(site.mutations) for site in focal_sites
        ),
        "observed_metadata_entry_count": len(observed_entries),
        "observed_mutation_types": observed_types,
    }
    if not selected_entries:
        if malformed_metadata:
            raise SelectedFocalIdentityError(
                "selected focal identity metadata is absent or malformed"
            )
        status = (
            "selected_focal_mutation_absent"
            if not focal_sites
            else "selected_focal_mutation_absent_background_only"
        )
        return {**base, "status": status}
    if len(focal_sites) != 1:
        raise SelectedFocalIdentityError(
            "selected focal identity requires exactly one focal site"
        )
    site = focal_sites[0]
    if len(site.mutations) != 1 or len(observed_entries) != 1:
        raise SelectedFocalIdentityError(
            "selected focal identity requires one nonrecurrent mutation and metadata entry"
        )
    if len(selected_entries) != 1:
        raise SelectedFocalIdentityError(
            "selected focal identity contains multiple type-1 mutation entries"
        )
    mutation, entry = selected_entries[0]
    required_metadata = {
        "mutation_type",
        "selection_coeff",
        "subpopulation",
        "slim_time",
        "nucleotide",
    }
    missing = sorted(required_metadata.difference(entry))
    if missing:
        raise SelectedFocalIdentityError(
            "selected focal identity metadata is missing: " + ", ".join(missing)
        )
    if mutation.parent != tskit.NULL:
        raise SelectedFocalIdentityError(
            "selected focal identity must be a single nonrecurrent root mutation"
        )
    expected_source = int(expectation["source_subpopulation"])
    observed_source = int(entry["subpopulation"])
    if observed_source != expected_source:
        raise SelectedFocalIdentityError(
            "selected focal identity source subpopulation differs from expectation: "
            f"observed {observed_source}, expected {expected_source}"
        )
    expected_time = float(expectation["realized_origin_time_generations"])
    observed_time = float(mutation.time)
    if not math.isclose(observed_time, expected_time, rel_tol=0.0, abs_tol=1e-9):
        raise SelectedFocalIdentityError(
            "selected focal identity origin time differs from expectation: "
            f"observed {observed_time:g}, expected {expected_time:g} generations"
        )
    expected_coefficient = float(expectation["expected_selection_coeff"])
    observed_coefficient = float(entry["selection_coeff"])
    if not math.isclose(
        observed_coefficient, expected_coefficient, rel_tol=0.0, abs_tol=1e-12
    ):
        raise SelectedFocalIdentityError(
            "selected focal identity mutation selection coefficient is invalid"
        )
    nucleotide = int(entry["nucleotide"])
    if nucleotide not in range(4) or mutation.derived_state != "ACGT"[nucleotide]:
        raise SelectedFocalIdentityError(
            "selected focal identity nucleotide metadata does not match derived state"
        )
    if site.ancestral_state == mutation.derived_state:
        raise SelectedFocalIdentityError(
            "selected focal identity ancestral and derived states are identical"
        )
    try:
        slim = ts.metadata["SLiM"]
        cycle = int(slim["cycle"])
        q_values = slim["user_metadata"]["Q"]
        q = float(q_values[0])
    except (KeyError, TypeError, ValueError, IndexError) as error:
        raise SelectedFocalIdentityError(
            "selected focal identity lacks SLiM cycle/Q provenance"
        ) from error
    expected_q = float(expectation["slim_scaling_factor"])
    if not math.isclose(q, expected_q, rel_tol=0.0, abs_tol=1e-12):
        raise SelectedFocalIdentityError(
            f"selected focal identity Q differs: observed {q:g}, expected {expected_q:g}"
        )
    slim_time = int(entry["slim_time"])
    expected_slim_time = round(cycle - observed_time / q)
    if slim_time != expected_slim_time:
        raise SelectedFocalIdentityError(
            "selected focal identity slim_time is inconsistent with cycle, age, and Q"
        )
    return {
        **base,
        "status": "valid",
        "observed_mutation_type": int(entry["mutation_type"]),
        "observed_source_subpopulation": observed_source,
        "observed_origin_time_generations": observed_time,
        "observed_selection_coeff": observed_coefficient,
        "observed_slim_time": slim_time,
        "observed_slim_cycle": cycle,
        "observed_slim_scaling_factor": q,
        "observed_ancestral_state": site.ancestral_state,
        "observed_derived_state": mutation.derived_state,
        "observed_nucleotide": nucleotide,
    }


def validate_selected_focal_identity(
    ts: tskit.TreeSequence,
    expectation: Mapping[str, Any],
) -> dict[str, Any]:
    """Require the focal allele to be the prescribed selected SLiM mutation."""

    record = inspect_selected_focal_identity(ts, expectation)
    if record["status"] != "valid":
        descriptor = (
            "background-only focal site"
            if record["status"].endswith("background_only")
            else "selected mutation absent at the focal site"
        )
        raise SelectedFocalIdentityError(
            f"selected focal identity failed: {descriptor}; "
            f"observed mutation types={record['observed_mutation_types']}"
        )
    return record


def _neutral_focal_expectation(record: Mapping[str, Any]) -> dict[str, Any]:
    """Build the immutable focal identity contract for one neutral unit."""

    expectation: dict[str, Any] = {
        "schema": "gamma-smc.neutral-focal-identity/v1",
        "focal_position_0based": int(FOCAL_POSITION_BP),
        "ancestral_state": "A",
        "derived_state": "G",
        "single_nonrecurrent_mutation": True,
    }
    if str(record.get("demography_kind", "no_introgression")) == "introgression":
        expectation.update(
            {
                "expected_origin_time_generations": float(
                    INTROGRESSION_MUTATION_AGE_GENERATIONS
                ),
                "expected_source_population": ARCHAIC_POPULATION,
                "pulse_route": {
                    "time_generations": float(INTROGRESSION_PULSE_GENERATIONS),
                    "source_population": "Loschbour",
                    "destination_population": ARCHAIC_POPULATION,
                    "minimum_descendant_lineages": 1,
                },
                "human_neanderthal_split_generations": float(
                    HUMAN_NEANDERTHAL_SPLIT_GENERATIONS
                ),
                "archaic_specific_no_ils_by_construction": True,
            }
        )
    return expectation


def _focal_semantic_validation_contract(
    record: Mapping[str, Any],
    repo_root: Path,
    *,
    slim_scaling_factor: float,
    selected_pool_diploids: int,
    han_selection_calibration: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    if str(record["simulation_class"]) == "selected":
        model, sweep, _, _ = _selected_objects(
            record,
            repo_root,
            slim_scaling_factor,
            selected_pool_diploids,
            han_selection_calibration=han_selection_calibration,
        )
        identity = _selected_focal_expectation(model, sweep, slim_scaling_factor)
    else:
        identity = _neutral_focal_expectation(record)
    return {
        "schema": FOCAL_SEMANTIC_VALIDATION_SCHEMA,
        "simulation_class": str(record["simulation_class"]),
        "sample_diploids": int(record["sample_diploids"]),
        "exact_sample_alt_count": int(record["exact_sample_alt_count"]),
        "minimum_hom_ref_diploids": 2,
        "minimum_heterozygous_diploids": 1,
        "minimum_hom_alt_diploids": 2,
        "identity": identity,
        "tree_manifest_exact_match": True,
    }


def _lineage_population_at_time(
    ts: tskit.TreeSequence,
    *,
    node: int,
    position: float,
    time: float,
) -> str:
    """Reconstruct one lineage's population along its focal branch."""

    tree = ts.at(float(position))
    parent = int(tree.parent(int(node)))
    if parent == tskit.NULL or not (
        float(ts.node(int(node)).time) <= float(time) < float(ts.node(parent).time)
    ):
        raise ValueError("neutral focal mutation time is outside its tree branch")
    current = _tree_sequence_population_name(ts, ts.node(int(node)).population)
    migrations = sorted(
        (
            migration
            for migration in ts.migrations()
            if int(migration.node) == int(node)
            and float(migration.left) <= float(position) < float(migration.right)
            and float(ts.node(int(node)).time) <= float(migration.time) <= float(time)
        ),
        key=lambda migration: float(migration.time),
    )
    for migration in migrations:
        source = _tree_sequence_population_name(ts, migration.source)
        destination = _tree_sequence_population_name(ts, migration.dest)
        if source != current:
            raise ValueError(
                "neutral focal lineage migration source differs from its "
                "reconstructed population"
            )
        current = destination
    return current


def _validate_neutral_focal_identity(
    ts: tskit.TreeSequence, expectation: Mapping[str, Any]
) -> dict[str, Any]:
    sites = _focal_sites(ts)
    if len(sites) != 1 or len(sites[0].mutations) != 1:
        raise ValueError("neutral focal identity requires one site and one mutation")
    site = sites[0]
    mutation = site.mutations[0]
    if (
        site.ancestral_state != "A"
        or mutation.derived_state != "G"
        or mutation.parent != tskit.NULL
        or not math.isfinite(float(mutation.time))
    ):
        raise ValueError(
            "neutral focal identity is not the constructed A-to-G mutation"
        )
    result: dict[str, Any] = {
        "schema": "gamma-smc.neutral-focal-identity/v1",
        "status": "valid",
        "focal_position_0based": int(FOCAL_POSITION_BP),
        "ancestral_state": site.ancestral_state,
        "derived_state": mutation.derived_state,
        "origin_time_generations": float(mutation.time),
        "mutation_metadata_empty": mutation.metadata in (b"", None),
    }
    if "expected_origin_time_generations" not in expectation:
        return result

    expected_time = float(expectation["expected_origin_time_generations"])
    observed_time = float(mutation.time)
    if not math.isclose(observed_time, expected_time, rel_tol=0.0, abs_tol=1e-9):
        raise ValueError(
            "introgressed neutral focal mutation age differs from its fixed contract"
        )
    source_population = _lineage_population_at_time(
        ts,
        node=int(mutation.node),
        position=float(site.position),
        time=observed_time,
    )
    expected_source = str(expectation["expected_source_population"])
    if source_population != expected_source:
        raise ValueError(
            "introgressed neutral focal branch is not resident in the expected "
            f"source population: observed {source_population}, expected "
            f"{expected_source}"
        )
    route = expectation.get("pulse_route")
    if not isinstance(route, Mapping):
        raise TypeError("introgressed neutral focal pulse-route contract is absent")
    tree = ts.at(float(site.position))
    matching_migrations = [
        migration
        for migration in ts.migrations()
        if float(migration.left) <= float(site.position) < float(migration.right)
        and math.isclose(
            float(migration.time),
            float(route["time_generations"]),
            rel_tol=0.0,
            abs_tol=1e-9,
        )
        and _tree_sequence_population_name(ts, migration.source)
        == str(route["source_population"])
        and _tree_sequence_population_name(ts, migration.dest)
        == str(route["destination_population"])
        and tree.is_descendant(int(migration.node), int(mutation.node))
    ]
    minimum_lineages = int(route["minimum_descendant_lineages"])
    lineage_nodes = sorted({int(migration.node) for migration in matching_migrations})
    if len(lineage_nodes) < minimum_lineages:
        raise ValueError(
            "introgressed neutral focal allele lacks reconstructable descendant "
            "pulse-route evidence"
        )
    split_time = float(expectation["human_neanderthal_split_generations"])
    if not float(route["time_generations"]) < observed_time < split_time:
        raise ValueError("introgressed neutral focal age does not exclude ILS")
    result.update(
        {
            "observed_source_population": source_population,
            "pulse_route_evidence": {
                "time_generations": float(route["time_generations"]),
                "source_population": str(route["source_population"]),
                "destination_population": str(route["destination_population"]),
                "descendant_lineage_nodes": lineage_nodes,
                "n_descendant_lineages": len(lineage_nodes),
            },
            "archaic_specific_no_ils_by_construction": True,
        }
    )
    return result


def _validate_unit_semantics(
    ts: tskit.TreeSequence,
    pair_table: pd.DataFrame,
    record: Mapping[str, Any],
    validation_contract: Mapping[str, Any],
) -> dict[str, Any]:
    """Validate focal identity, exact AF/QC, and tree-manifest equivalence."""

    if validation_contract.get("schema") != FOCAL_SEMANTIC_VALIDATION_SCHEMA:
        raise ValueError("focal semantic validation contract is incompatible")
    simulation_class = str(record["simulation_class"])
    if str(validation_contract.get("simulation_class")) != simulation_class:
        raise ValueError("focal validation simulation class differs from the unit")
    for field in ("sample_diploids", "exact_sample_alt_count"):
        if int(validation_contract.get(field, -1)) != int(record[field]):
            raise ValueError(f"focal validation {field} differs from the unit")
    expected_qc = {
        "minimum_hom_ref_diploids": 2,
        "minimum_heterozygous_diploids": 1,
        "minimum_hom_alt_diploids": 2,
    }
    if (
        any(
            int(validation_contract.get(field, -1)) != minimum
            for field, minimum in expected_qc.items()
        )
        or validation_contract.get("tree_manifest_exact_match") is not True
    ):
        raise ValueError("focal validation genotype/manifest contract is incompatible")
    if simulation_class == "selected":
        identity = validate_selected_focal_identity(ts, validation_contract["identity"])
    else:
        expected_identity = _neutral_focal_expectation(record)
        if _canonical_sha256(validation_contract.get("identity")) != (
            _canonical_sha256(expected_identity)
        ):
            raise ValueError("neutral focal identity contract is incompatible")
        identity = _validate_neutral_focal_identity(ts, expected_identity)
    expected = build_diploid_pair_table(ts, FOCAL_POSITION_BP)
    columns = list(expected.columns)
    if list(pair_table.columns) != columns:
        raise ValueError("sample manifest columns/order do not exactly match the tree")
    try:
        pd.testing.assert_frame_equal(
            pair_table[columns].reset_index(drop=True),
            expected[columns].reset_index(drop=True),
            check_dtype=False,
            check_exact=True,
        )
    except AssertionError as error:
        raise ValueError(
            "sample manifest does not exactly match the focal tree"
        ) from error
    counts = expected["focal_selected_allele_count"].to_numpy(dtype=int)
    exact = int(record["exact_sample_alt_count"])
    if len(expected) != int(record["sample_diploids"]):
        raise ValueError("tree/manifest diploid count differs from the unit")
    if not _passes_genotype_qc(counts, exact):
        raise ValueError("tree/manifest exact AF or genotype-class QC failed")
    genotype_counts = expected["genotype_class"].value_counts().to_dict()
    return {
        "schema": FOCAL_SEMANTIC_VALIDATION_SCHEMA,
        "status": "valid",
        "identity": identity,
        "tree_manifest_exact_match": True,
        "sample_diploids": len(expected),
        "sample_alt_count": int(counts.sum()),
        "sample_af": float(counts.mean() / 2),
        "genotype_counts": {
            str(key): int(value) for key, value in genotype_counts.items()
        },
    }


def _log_choose(n: int, k: int) -> float:
    if k < 0 or k > n:
        return -math.inf
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def select_exact_sample_panel(
    ts: tskit.TreeSequence,
    *,
    sample_diploids: int,
    exact_alt_count: int,
    seed: int,
) -> tuple[tskit.TreeSequence, dict[str, Any]]:
    """Uniformly draw a diploid subset conditional on size, AF, and QC.

    Feasible genotype-count triples are weighted by their exact number of
    underlying subsets, ``C(N0,n0) C(N1,n1) C(N2,n2)``. Individuals are then
    sampled uniformly within genotype classes. The result is therefore uniform
    over all eligible subsets rather than the first panel found by rejection.
    """

    pool = build_diploid_pair_table(ts, FOCAL_POSITION_BP)
    pool_counts = pool["focal_selected_allele_count"].to_numpy(dtype=int)
    available = {
        value: int(np.count_nonzero(pool_counts == value)) for value in (0, 1, 2)
    }
    triples = feasible_exact_panel_genotype_triples(
        pool_counts,
        sample_diploids=sample_diploids,
        exact_alt_count=exact_alt_count,
    )
    if not triples:
        raise SimulationExhausted(
            "candidate pool has no exact-AF subset satisfying genotype QC"
        )
    values = np.asarray(
        [record["log_number_of_subsets"] for record in triples], dtype=float
    )
    weights = np.exp(values - values.max())
    rng = np.random.default_rng(seed)
    triple_index = int(rng.choice(len(triples), p=weights / weights.sum()))
    chosen_triple = triples[triple_index]
    chosen_rows: list[int] = []
    chosen_numbers = (
        int(chosen_triple["hom_ref"]),
        int(chosen_triple["heterozygous"]),
        int(chosen_triple["hom_alt"]),
    )
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
        # tskit cannot simplify tables containing migrations. Preserve the
        # complete ancestry/migration record needed to prove the Han pulse
        # route, while making the chosen diploid nodes the only samples seen by
        # VCF/Gamma-SMC and downstream pair construction.
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
    manifest = build_diploid_pair_table(panel, FOCAL_POSITION_BP)
    observed = manifest["focal_selected_allele_count"].to_numpy(dtype=int)
    if len(manifest) != sample_diploids or not _passes_genotype_qc(
        observed, exact_alt_count
    ):
        raise RuntimeError("exact sample-panel construction failed validation")
    return panel, {
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
    }


def simulate_selected_unit(
    unit: Mapping[str, Any],
    repo_root: str | Path,
    unit_dir: str | Path,
    *,
    slim_path: str | Path,
    slim_scaling_factor: float = DEFAULT_SLIM_SCALING_FACTOR,
    slim_burn_in: float = DEFAULT_SLIM_BURN_IN,
    max_draws: int = DEFAULT_SELECTED_MAX_DRAWS,
    selected_pool_diploids: int = DEFAULT_SELECTED_POOL_DIPLOIDS,
    draw_timeout_seconds: float = DEFAULT_SELECTED_DRAW_TIMEOUT_SECONDS,
    cumulative_timeout_seconds: float = DEFAULT_SELECTED_CUMULATIVE_TIMEOUT_SECONDS,
    han_selection_calibration: Mapping[str, Any] | None = None,
) -> tuple[tskit.TreeSequence, pd.DataFrame, dict[str, Any], list[dict[str, Any]]]:
    record = _validate_unit(unit)
    if record["simulation_class"] != "selected":
        raise ValueError("simulate_selected_unit received a neutral unit")
    directory = Path(unit_dir)
    directory.mkdir(parents=True, exist_ok=True)
    model, sweep, samples, model_record = _selected_objects(
        record,
        repo_root,
        slim_scaling_factor,
        selected_pool_diploids,
        han_selection_calibration=han_selection_calibration,
    )
    focal_expectation = _selected_focal_expectation(model, sweep, slim_scaling_factor)
    engine = stdpopsim.get_engine("slim")
    ledger_path = directory / "attempts.tsv"
    attempts: list[dict[str, Any]] = []
    if ledger_path.is_file():
        attempts = pd.read_csv(ledger_path, sep="\t").to_dict(orient="records")
        observed = [int(row["attempt_zero_based"]) for row in attempts]
        if observed != list(range(len(attempts))):
            raise ValueError("selected attempt ledger is not contiguous from zero")
    unit_started = perf_counter()
    previous_elapsed = float(
        sum(float(row.get("elapsed_seconds", 0.0)) for row in attempts)
    )
    for attempt in range(len(attempts), int(max_draws)):
        elapsed_before_draw = previous_elapsed + perf_counter() - unit_started
        remaining = float(cumulative_timeout_seconds) - elapsed_before_draw
        if remaining <= 0:
            raise SimulationExhausted(
                f"{record['unit_id']} exhausted its cumulative selected-simulation "
                f"budget of {cumulative_timeout_seconds:g} seconds"
            )
        effective_timeout = min(float(draw_timeout_seconds), remaining)
        seed = _stable_seed(int(record["seed"]), f"selected:{attempt}")
        logfile = directory / f"slim_draw_{attempt:03d}.csv"
        started = perf_counter()
        try:
            with _wall_clock_timeout(effective_timeout):
                ts = engine.simulate(
                    model,
                    sweep.contig,
                    samples,
                    seed=seed,
                    extended_events=list(sweep.extended_events),
                    slim_path=str(Path(slim_path).resolve()),
                    slim_scaling_factor=slim_scaling_factor,
                    slim_burn_in=slim_burn_in,
                    logfile=str(logfile),
                    logfile_interval=100,
                    keep_mutation_ids_as_alleles=False,
                )
        except SimulationDrawTimeout as error:
            attempts.append(
                {
                    "attempt_zero_based": attempt,
                    "seed": seed,
                    "status": "timed_out",
                    "timeout_seconds": effective_timeout,
                    "elapsed_seconds": perf_counter() - started,
                    "error": str(error),
                    "accepted": False,
                }
            )
            _atomic_frame(ledger_path, pd.DataFrame(attempts))
            continue
        focal_identity = inspect_selected_focal_identity(ts, focal_expectation)
        if focal_identity["status"] != "valid":
            # Once terminal population-frequency rejection is removed, loss of
            # the selected type-1 mutation is an expected rejected outer draw.
            # A coincident recapitation/background mutation at 5 Mb must not be
            # mistaken for the selected allele.
            attempts.append(
                {
                    "attempt_zero_based": attempt,
                    "seed": seed,
                    "status": "complete",
                    "timeout_seconds": effective_timeout,
                    "candidate_pool_alt_count": 0,
                    "candidate_pool_af": 0.0,
                    "sample_alt_count": None,
                    "sample_af": None,
                    "n_hom_ref": None,
                    "n_heterozygous": None,
                    "n_hom_alt": None,
                    "elapsed_seconds": perf_counter() - started,
                    "accepted": False,
                    "rejection_reason": focal_identity["status"],
                    "observed_focal_site_count": focal_identity[
                        "observed_focal_site_count"
                    ],
                    "observed_focal_mutation_types": json.dumps(
                        focal_identity["observed_mutation_types"],
                        separators=(",", ":"),
                    ),
                    "selected_focal_identity_schema": (SELECTED_FOCAL_IDENTITY_SCHEMA),
                }
            )
            _atomic_frame(ledger_path, pd.DataFrame(attempts))
            continue
        pool_table = build_diploid_pair_table(ts, FOCAL_POSITION_BP)
        pool_counts = pool_table["focal_selected_allele_count"].to_numpy(dtype=int)
        if len(pool_counts) != int(selected_pool_diploids):
            raise ValueError("selected candidate pool has the wrong diploid count")
        pool_gate = candidate_pool_frequency_gate(record, pool_counts)
        exact = int(record["exact_sample_alt_count"])
        panel_seed = _stable_seed(int(record["seed"]), f"selected-panel:{attempt}")
        feasible_triples = (
            feasible_exact_panel_genotype_triples(
                pool_counts,
                sample_diploids=int(record["sample_diploids"]),
                exact_alt_count=exact,
            )
            if pool_gate["passed"]
            else []
        )
        panel_ts = None
        panel_record = None
        panel_identity = None
        rejection_reason: str | None = None
        if not pool_gate["passed"]:
            pair_table = pool_table
            counts = pool_counts
            accepted = False
            rejection_reason = "candidate_pool_af_outside_prespecified_band"
        elif not feasible_triples:
            pair_table = pool_table
            counts = pool_counts
            accepted = False
            rejection_reason = "candidate_pool_has_no_exact_af_genotype_qc_panel"
        else:
            panel_ts, panel_record = select_exact_sample_panel(
                ts,
                sample_diploids=int(record["sample_diploids"]),
                exact_alt_count=exact,
                seed=panel_seed,
            )
            panel_identity = validate_selected_focal_identity(
                panel_ts, focal_expectation
            )
            pair_table = build_diploid_pair_table(panel_ts, FOCAL_POSITION_BP)
            counts = pair_table["focal_selected_allele_count"].to_numpy(dtype=int)
            accepted = True
        attempts.append(
            {
                "attempt_zero_based": attempt,
                "seed": seed,
                "status": "complete",
                "timeout_seconds": effective_timeout,
                "candidate_pool_alt_count": int(pool_counts.sum()),
                "candidate_pool_af": float(pool_counts.mean() / 2),
                "candidate_pool_minimum_alt_count": int(
                    pool_gate["minimum_alt_count_inclusive"]
                ),
                "candidate_pool_maximum_alt_count": int(
                    pool_gate["maximum_alt_count_inclusive"]
                ),
                "candidate_pool_frequency_gate_passed": bool(pool_gate["passed"]),
                "eligible_genotype_count_triples": len(feasible_triples),
                "sample_alt_count": int(counts.sum()) if accepted else None,
                "sample_af": float(counts.mean() / 2) if accepted else None,
                "n_hom_ref": int(np.count_nonzero(counts == 0)) if accepted else None,
                "n_heterozygous": int(np.count_nonzero(counts == 1))
                if accepted
                else None,
                "n_hom_alt": int(np.count_nonzero(counts == 2)) if accepted else None,
                "elapsed_seconds": perf_counter() - started,
                "accepted": accepted,
                "rejection_reason": rejection_reason,
            }
        )
        _atomic_frame(ledger_path, pd.DataFrame(attempts))
        if accepted:
            provenance = {
                "engine": "stdpopsim_slim",
                "stdpopsim_version": stdpopsim.__version__,
                "pyslim_version": pyslim.__version__,
                "slim_path": str(Path(slim_path).resolve()),
                "slim_sha256": sha256_file(slim_path),
                "slim_scaling_factor": slim_scaling_factor,
                "slim_burn_in": slim_burn_in,
                "draw_timeout_seconds": float(draw_timeout_seconds),
                "cumulative_timeout_seconds": float(cumulative_timeout_seconds),
                "accepted_draw_zero_based": attempt,
                "accepted_seed": seed,
                "candidate_pool_frequency_gate": pool_gate,
                "sample_panel": panel_record,
                "matched_selected_neutral_ascertainment": (
                    "same_candidate_pool_size_and_AF_band_then_same_exact_k_"
                    "panel_sampler"
                ),
                "selected_focal_identity": panel_identity,
                **model_record,
            }
            return panel_ts, pair_table, provenance, attempts
    raise SimulationExhausted(
        f"{record['unit_id']} produced no exact sampled AF in {max_draws} "
        "completed SLiM draws"
    )


def _completion_artifacts(unit_dir: Path) -> SimulationArtifacts:
    return SimulationArtifacts(
        unit_dir=unit_dir,
        tree_path=unit_dir / "simulation.trees",
        pair_table_path=unit_dir / "sample_manifest.tsv",
        overall_pairs_path=unit_dir / "pairs" / "overall.pairs.tsv",
        truth_profiles_path=unit_dir / "truth_profiles.tsv.gz",
        truth_class_summaries_path=unit_dir / "truth_class_summaries.tsv",
        completion_path=unit_dir / "simulation_complete.json",
        cache_hit=False,
    )


def _validate_completion(
    artifacts: SimulationArtifacts, contract: Mapping[str, Any]
) -> SimulationArtifacts:
    try:
        completion = json.loads(artifacts.completion_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError("simulation completion is unreadable") from error
    if (
        completion.get("schema") != SCHEMA_VERSION
        or completion.get("status") != "complete"
    ):
        raise ValueError("simulation completion schema/status is incompatible")
    recorded_contract = completion.get("contract")
    if not isinstance(recorded_contract, Mapping) or completion.get(
        "contract_sha256"
    ) != _canonical_sha256(recorded_contract):
        raise ValueError("simulation completion contract self-check is invalid")
    if _canonical_sha256(recorded_contract) != _canonical_sha256(contract):
        raise ValueError("simulation completion contract is incompatible")
    contract_path = artifacts.unit_dir / "simulation_contract.json"
    try:
        disk_contract = json.loads(contract_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError("completed simulation contract file is unreadable") from error
    if _canonical_sha256(disk_contract) != _canonical_sha256(contract):
        raise ValueError("completed simulation contract file is incompatible")
    output_paths = {
        "tree": artifacts.tree_path,
        "pair_table": artifacts.pair_table_path,
        "overall_pairs": artifacts.overall_pairs_path,
        "truth_profiles": artifacts.truth_profiles_path,
        "truth_class_summaries": artifacts.truth_class_summaries_path,
    }
    for label, path in output_paths.items():
        output_record = completion.get("outputs", {}).get(label, {})
        expected_relative = path.relative_to(artifacts.unit_dir).as_posix()
        if (
            not isinstance(output_record, Mapping)
            or str(output_record.get("path", "")).replace("\\", "/")
            != expected_relative
            or not path.is_file()
            or sha256_file(path) != output_record.get("sha256")
        ):
            raise ValueError(f"completed simulation {label} failed checksum validation")
    ts = tskit.load(artifacts.tree_path)
    if not math.isclose(ts.sequence_length, SEQUENCE_LENGTH_BP):
        raise ValueError("cached simulation tree length is invalid")
    try:
        pair_table = pd.read_csv(artifacts.pair_table_path, sep="\t")
    except Exception as error:
        raise ValueError("cached simulation sample manifest is unreadable") from error
    semantic = _validate_unit_semantics(
        ts,
        pair_table,
        contract["unit"],
        contract["parameters"]["focal_semantic_validation"],
    )
    provenance = completion.get("provenance")
    if not isinstance(provenance, Mapping):
        raise TypeError("simulation completion provenance is malformed")
    recorded_semantic = provenance.get("focal_semantic_validation")
    if _canonical_sha256(recorded_semantic) != _canonical_sha256(semantic):
        raise ValueError(
            "simulation completion focal semantic validation is absent or stale"
        )
    validated_pool_gate = _validate_candidate_pool_frequency_gate_record(
        provenance.get("candidate_pool_frequency_gate", {}),
        contract["unit"],
        expected_pool_diploids=int(contract["parameters"]["candidate_pool_diploids"]),
    )
    _validate_exact_sample_panel_record(
        provenance.get("sample_panel", {}),
        contract["unit"],
        validated_pool_gate,
    )
    if str(contract["unit"]["simulation_class"]) == "selected" and (
        _canonical_sha256(provenance.get("selected_focal_identity"))
        != _canonical_sha256(semantic["identity"])
    ):
        raise ValueError(
            "simulation completion selected focal identity provenance is invalid"
        )
    if str(contract["unit"]["simulation_class"]) == "neutral" and (
        _canonical_sha256(provenance.get("neutral_focal_identity"))
        != _canonical_sha256(semantic["identity"])
    ):
        raise ValueError(
            "simulation completion neutral focal identity provenance is invalid"
        )
    return SimulationArtifacts(**{**artifacts.__dict__, "cache_hit": True})


def _simulate_unit_unlocked(
    unit: Mapping[str, Any],
    repo_root: str | Path,
    campaign_dir: str | Path,
    *,
    slim_path: str | Path | None = None,
    neutral_max_attempts: int = DEFAULT_NEUTRAL_MAX_ATTEMPTS,
    selected_max_draws: int = DEFAULT_SELECTED_MAX_DRAWS,
    slim_scaling_factor: float = DEFAULT_SLIM_SCALING_FACTOR,
    slim_burn_in: float = DEFAULT_SLIM_BURN_IN,
    selected_draw_timeout_seconds: float = DEFAULT_SELECTED_DRAW_TIMEOUT_SECONDS,
    selected_cumulative_timeout_seconds: float = (
        DEFAULT_SELECTED_CUMULATIVE_TIMEOUT_SECONDS
    ),
    han_selection_calibration: Mapping[str, Any] | None = None,
) -> SimulationArtifacts:
    """Simulate one unit with checksum-valid, idempotent terminal artifacts."""

    record = _validate_unit(unit)
    root = Path(repo_root).resolve()
    campaign = Path(campaign_dir).resolve()
    if any("onedrive" in part.casefold() for part in campaign.parts):
        raise ValueError("focused campaign output must not be inside OneDrive")
    try:
        campaign.relative_to(root)
    except ValueError as error:
        raise ValueError(
            "focused campaign output must remain inside repo_root"
        ) from error
    unit_dir = campaign / "work" / str(record["unit_id"])
    unit_dir.mkdir(parents=True, exist_ok=True)
    artifacts = _completion_artifacts(unit_dir)
    focal_validation_contract = _focal_semantic_validation_contract(
        record,
        root,
        slim_scaling_factor=slim_scaling_factor,
        selected_pool_diploids=DEFAULT_SELECTED_POOL_DIPLOIDS,
        han_selection_calibration=han_selection_calibration,
    )
    contract = {
        "schema": SCHEMA_VERSION,
        "unit": record,
        "parameters": {
            "slim_scaling_factor": float(slim_scaling_factor),
            "slim_burn_in": float(slim_burn_in),
            "candidate_pool_diploids": DEFAULT_SELECTED_POOL_DIPLOIDS,
            "selected_pool_diploids": DEFAULT_SELECTED_POOL_DIPLOIDS,
            "neutral_pool_diploids": DEFAULT_NEUTRAL_POOL_DIPLOIDS,
            "candidate_pool_frequency_gate": {
                "schema": CANDIDATE_POOL_FREQUENCY_GATE_SCHEMA,
                "scope": "selected_and_neutral",
                "interval_source": (
                    "unit.population_af_lower_to_population_af_upper_inclusive"
                ),
                "order": "before_exact_sample_panel",
            },
            "neutral_ancestry_batch_size": DEFAULT_NEUTRAL_ANCESTRY_BATCH_SIZE,
            "neutral_ancestry_model": NEUTRAL_ANCESTRY_MODEL_NAME,
            "neutral_ancestry_model_keywords": list(NEUTRAL_ANCESTRY_MODEL_KEYWORDS),
            "neutral_recent_dtwf_duration_generations": (
                NEUTRAL_RECENT_DTWF_DURATION_GENERATIONS
            ),
            "neutral_eas_ascertainment_length_bp": SEQUENCE_LENGTH_BP,
            "neutral_han_ascertainment_length_bp": (
                INTROGRESSION_ASCERTAINMENT_LENGTH_BP
            ),
            "neutral_candidate_index": {
                "no_introgression": "marginal_tree_at_fixed_5Mb_site",
                "introgression": (
                    "pulse_migrations_then_branch_crossing_2400_generations"
                ),
            },
            "introgression_source_presence_epsilon": (
                INTROGRESSION_SOURCE_PRESENCE_EPSILON
            ),
            "introgression_mutation_age_generations": (
                INTROGRESSION_MUTATION_AGE_GENERATIONS
            ),
            "han_selection_calibration": (
                dict(han_selection_calibration)
                if record["simulation_class"] == "selected"
                and record["demography_id"] == "ancient_eurasia_han_introgression"
                else None
            ),
            "focal_semantic_validation": focal_validation_contract,
        },
        "software": {
            "python": platform.python_version(),
            "msprime": msprime.__version__,
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "pyslim": pyslim.__version__,
            "stdpopsim": stdpopsim.__version__,
            "tskit": tskit.__version__,
        },
        "implementation": _implementation_contract(
            root,
            simulation_class=str(record["simulation_class"]),
            slim_path=slim_path,
        ),
    }
    if artifacts.completion_path.is_file():
        return _validate_completion(artifacts, contract)
    contract_path = unit_dir / "simulation_contract.json"
    if contract_path.is_file():
        try:
            recorded_contract = json.loads(contract_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            raise ValueError("partial simulation contract is unreadable") from error
        if _canonical_sha256(recorded_contract) != _canonical_sha256(contract):
            raise ValueError(
                "partial simulation was created under an incompatible contract"
            )
    else:
        if (unit_dir / "attempts.tsv").exists() or (
            unit_dir / "simulation_failed.json"
        ).exists():
            raise ValueError(
                "partial simulation lacks a contract; quarantine before retrying"
            )
        _atomic_json(contract_path, contract)
    if any(path.exists() for path in (artifacts.tree_path, artifacts.pair_table_path)):
        raise ValueError("partial simulation artifacts exist without completion")
    started = perf_counter()
    try:
        if record["simulation_class"] == "neutral":
            ts, pair_table, provenance, attempts = simulate_neutral_unit(
                record,
                root,
                unit_dir,
                max_attempts=neutral_max_attempts,
                ancestry_batch_size=DEFAULT_NEUTRAL_ANCESTRY_BATCH_SIZE,
                candidate_pool_diploids=DEFAULT_NEUTRAL_POOL_DIPLOIDS,
            )
        else:
            if slim_path is None:
                raise ValueError("selected simulation requires --slim-bin")
            ts, pair_table, provenance, attempts = simulate_selected_unit(
                record,
                root,
                unit_dir,
                slim_path=slim_path,
                slim_scaling_factor=slim_scaling_factor,
                slim_burn_in=slim_burn_in,
                max_draws=selected_max_draws,
                selected_pool_diploids=DEFAULT_SELECTED_POOL_DIPLOIDS,
                draw_timeout_seconds=selected_draw_timeout_seconds,
                cumulative_timeout_seconds=selected_cumulative_timeout_seconds,
                han_selection_calibration=han_selection_calibration,
            )
    except Exception as error:
        _atomic_json(
            unit_dir / "simulation_failed.json",
            {
                "schema": SCHEMA_VERSION,
                "status": "failed",
                "unit_id": record["unit_id"],
                "error_type": type(error).__name__,
                "error": str(error),
                "contract_sha256": _canonical_sha256(contract),
                "search_budget": {
                    "neutral_max_attempts": int(neutral_max_attempts),
                    "selected_max_draws": int(selected_max_draws),
                    "selected_draw_timeout_seconds": float(
                        selected_draw_timeout_seconds
                    ),
                    "selected_cumulative_timeout_seconds": float(
                        selected_cumulative_timeout_seconds
                    ),
                },
            },
        )
        raise
    semantic_validation = _validate_unit_semantics(
        ts,
        pair_table,
        record,
        focal_validation_contract,
    )
    provenance = dict(provenance)
    provenance["candidate_pool_frequency_gate"] = (
        _validate_candidate_pool_frequency_gate_record(
            provenance.get("candidate_pool_frequency_gate", {}),
            record,
            expected_pool_diploids=DEFAULT_SELECTED_POOL_DIPLOIDS,
        )
    )
    provenance["sample_panel"] = _validate_exact_sample_panel_record(
        provenance.get("sample_panel", {}),
        record,
        provenance["candidate_pool_frequency_gate"],
    )
    if record["simulation_class"] == "selected":
        if _canonical_sha256(
            provenance.get("selected_focal_identity")
        ) != _canonical_sha256(semantic_validation["identity"]):
            raise SelectedFocalIdentityError(
                "selected focal identity provenance differs from the accepted tree"
            )
        provenance["selected_focal_identity"] = semantic_validation["identity"]
    else:
        provenance["neutral_focal_identity"] = semantic_validation["identity"]
    provenance["focal_semantic_validation"] = semantic_validation
    _atomic_tree(artifacts.tree_path, ts)
    _atomic_frame(artifacts.pair_table_path, pair_table)
    pair_paths = write_gamma_pair_files(pair_table, unit_dir / "pairs")
    positions = np.unique(
        np.r_[np.arange(0, SEQUENCE_LENGTH_BP, 10_000), FOCAL_POSITION_BP]
    )
    truth = tree_truth_profiles(
        ts,
        positions,
        pair_table,
        thresholds_years=TMRCA_THRESHOLDS_YEARS,
        generation_time_years=GENERATION_TIME_YEARS,
    )
    prefix = {
        key: record[key]
        for key in (
            "unit_id",
            "demography_id",
            "demography_kind",
            "population",
            "simulation_class",
            "selection_coefficient",
            "target_allele_frequency",
            "replicate_index",
            "seed",
        )
    }
    truth = truth.assign(**prefix)
    _atomic_frame(artifacts.truth_profiles_path, truth)
    focal = truth[np.isclose(truth["position_0based"], FOCAL_POSITION_BP)]
    if len(focal) != len(TMRCA_THRESHOLDS_YEARS) * 4:
        raise ValueError("tree-truth focal profile is incomplete")
    truth_class_rows: list[dict[str, Any]] = []
    for (genotype_class, threshold), group in truth.groupby(
        ["genotype_class", "threshold_years"], sort=True
    ):
        center = focal[
            (focal["genotype_class"] == genotype_class)
            & np.isclose(focal["threshold_years"], threshold)
        ].iloc[0]
        truth_class_rows.append(
            {
                **prefix,
                "source": "tree_truth",
                "genotype_class": genotype_class,
                "threshold_years": float(threshold),
                "region_mean_p_tmrca_lt_threshold": float(
                    group["mean_p_tmrca_lt_threshold"].mean()
                ),
                "focal_mean_p_tmrca_lt_threshold": float(
                    center["mean_p_tmrca_lt_threshold"]
                ),
                "region_mean_tmrca_generations": float(
                    group["mean_tmrca_generations"].mean()
                ),
                "focal_mean_tmrca_generations": float(center["mean_tmrca_generations"]),
                "n_pairs": int(center["n_pairs"]),
                "focal_output_position_0based": int(FOCAL_POSITION_BP),
                "focal_offset_bp": 0,
            }
        )
    truth_class = pd.DataFrame(truth_class_rows)
    _atomic_frame(artifacts.truth_class_summaries_path, truth_class)
    outputs = {
        "tree": {
            "path": artifacts.tree_path.name,
            "sha256": sha256_file(artifacts.tree_path),
        },
        "pair_table": {
            "path": artifacts.pair_table_path.name,
            "sha256": sha256_file(artifacts.pair_table_path),
        },
        "overall_pairs": {
            "path": str(pair_paths["overall"].relative_to(unit_dir)),
            "sha256": sha256_file(pair_paths["overall"]),
        },
        "truth_profiles": {
            "path": artifacts.truth_profiles_path.name,
            "sha256": sha256_file(artifacts.truth_profiles_path),
        },
        "truth_class_summaries": {
            "path": artifacts.truth_class_summaries_path.name,
            "sha256": sha256_file(artifacts.truth_class_summaries_path),
        },
    }
    completion = {
        "schema": SCHEMA_VERSION,
        "status": "complete",
        "contract": contract,
        "contract_sha256": _canonical_sha256(contract),
        "provenance": provenance,
        "attempts_completed": len(attempts),
        "search_budget": {
            "neutral_max_attempts": int(neutral_max_attempts),
            "selected_max_draws": int(selected_max_draws),
            "selected_draw_timeout_seconds": float(selected_draw_timeout_seconds),
            "selected_cumulative_timeout_seconds": float(
                selected_cumulative_timeout_seconds
            ),
        },
        "elapsed_seconds": perf_counter() - started,
        "genotype_counts": pair_table["genotype_class"].value_counts().to_dict(),
        "sample_alt_count": int(pair_table["focal_selected_allele_count"].sum()),
        "sample_af": float(pair_table["focal_selected_allele_count"].mean() / 2),
        "outputs": outputs,
    }
    _atomic_json(artifacts.completion_path, completion)
    (unit_dir / "simulation_failed.json").unlink(missing_ok=True)
    return artifacts


def simulate_unit(
    unit: Mapping[str, Any],
    repo_root: str | Path,
    campaign_dir: str | Path,
    *,
    slim_path: str | Path | None = None,
    neutral_max_attempts: int = DEFAULT_NEUTRAL_MAX_ATTEMPTS,
    selected_max_draws: int = DEFAULT_SELECTED_MAX_DRAWS,
    slim_scaling_factor: float = DEFAULT_SLIM_SCALING_FACTOR,
    slim_burn_in: float = DEFAULT_SLIM_BURN_IN,
    selected_draw_timeout_seconds: float = DEFAULT_SELECTED_DRAW_TIMEOUT_SECONDS,
    selected_cumulative_timeout_seconds: float = (
        DEFAULT_SELECTED_CUMULATIVE_TIMEOUT_SECONDS
    ),
    han_selection_calibration: Mapping[str, Any] | None = None,
) -> SimulationArtifacts:
    """Simulate one unit while holding its restart-safe process lock."""

    record = _validate_unit(unit)
    root = Path(repo_root).resolve()
    campaign = Path(campaign_dir).resolve()
    if any("onedrive" in part.casefold() for part in campaign.parts):
        raise ValueError("focused campaign output must not be inside OneDrive")
    try:
        campaign.relative_to(root)
    except ValueError as error:
        raise ValueError(
            "focused campaign output must remain inside repo_root"
        ) from error
    unit_dir = campaign / "work" / str(record["unit_id"])
    with exclusive_unit_lock(
        unit_dir / ".simulation.lock",
        unit_id=str(record["unit_id"]),
        phase="simulation",
    ):
        return _simulate_unit_unlocked(
            record,
            root,
            campaign,
            slim_path=slim_path,
            neutral_max_attempts=neutral_max_attempts,
            selected_max_draws=selected_max_draws,
            slim_scaling_factor=slim_scaling_factor,
            slim_burn_in=slim_burn_in,
            selected_draw_timeout_seconds=selected_draw_timeout_seconds,
            selected_cumulative_timeout_seconds=selected_cumulative_timeout_seconds,
            han_selection_calibration=han_selection_calibration,
        )


def _simulate_task(task: Mapping[str, Any]) -> dict[str, Any]:
    try:
        artifacts = simulate_unit(**task)
        return {
            "unit_id": task["unit"]["unit_id"],
            "status": "cached" if artifacts.cache_hit else "complete",
            "completion_path": str(artifacts.completion_path),
            "error_type": "",
            "error": "",
        }
    except Exception as error:  # noqa: BLE001 - worker must return all task failures
        return {
            "unit_id": task["unit"]["unit_id"],
            "status": "failed",
            "completion_path": "",
            "error_type": type(error).__name__,
            "error": str(error),
        }


def simulate_units(
    units: pd.DataFrame,
    repo_root: str | Path,
    campaign_dir: str | Path,
    *,
    slim_path: str | Path | None = None,
    workers: int = 4,
    neutral_max_attempts: int = DEFAULT_NEUTRAL_MAX_ATTEMPTS,
    selected_max_draws: int = DEFAULT_SELECTED_MAX_DRAWS,
    selected_draw_timeout_seconds: float = DEFAULT_SELECTED_DRAW_TIMEOUT_SECONDS,
    selected_cumulative_timeout_seconds: float = (
        DEFAULT_SELECTED_CUMULATIVE_TIMEOUT_SECONDS
    ),
    han_selection_calibration_path: str | Path | None = None,
    slim_scaling_factor: float = DEFAULT_SLIM_SCALING_FACTOR,
) -> pd.DataFrame:
    """Execute a filtered unit frame in parallel and return terminal statuses."""

    root = Path(repo_root).resolve()
    campaign = Path(campaign_dir).resolve()
    han_selected = units[
        (units["simulation_class"].astype(str) == "selected")
        & (units["demography_id"].astype(str) == "ancient_eurasia_han_introgression")
    ]
    han_selection_calibration = None
    if not han_selected.empty:
        if han_selection_calibration_path is None:
            raise ValueError(
                "Han selected production requires --han-selection-calibration"
            )
        han_selection_calibration = load_frozen_han_selection_calibration(
            han_selection_calibration_path,
            repo_root=root,
            execution_units_path=campaign / "execution_units.tsv",
            slim_scaling_factor=slim_scaling_factor,
            candidate_pool_diploids=DEFAULT_SELECTED_POOL_DIPLOIDS,
        )
    tasks = [
        {
            "unit": row,
            "repo_root": str(root),
            "campaign_dir": str(campaign),
            "slim_path": str(Path(slim_path).resolve()) if slim_path else None,
            "neutral_max_attempts": neutral_max_attempts,
            "selected_max_draws": selected_max_draws,
            "selected_draw_timeout_seconds": selected_draw_timeout_seconds,
            "selected_cumulative_timeout_seconds": (
                selected_cumulative_timeout_seconds
            ),
            "slim_scaling_factor": float(slim_scaling_factor),
            "han_selection_calibration": han_selection_calibration,
        }
        for row in units.to_dict(orient="records")
    ]
    rows: list[dict[str, Any]] = []
    if int(workers) == 1:
        rows = [_simulate_task(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=int(workers)) as pool:
            futures = {pool.submit(_simulate_task, task): task for task in tasks}
            for future in as_completed(futures):
                rows.append(future.result())
    return pd.DataFrame(rows).sort_values("unit_id").reset_index(drop=True)


__all__ = [
    "CANDIDATE_POOL_FREQUENCY_GATE_SCHEMA",
    "DEFAULT_NEUTRAL_ANCESTRY_BATCH_SIZE",
    "DEFAULT_NEUTRAL_MAX_ATTEMPTS",
    "DEFAULT_NEUTRAL_POOL_DIPLOIDS",
    "DEFAULT_SELECTED_CUMULATIVE_TIMEOUT_SECONDS",
    "DEFAULT_SELECTED_DRAW_TIMEOUT_SECONDS",
    "DEFAULT_SELECTED_MAX_DRAWS",
    "DEFAULT_SELECTED_POOL_DIPLOIDS",
    "FOCAL_SEMANTIC_VALIDATION_SCHEMA",
    "INTROGRESSION_ASCERTAINMENT_LENGTH_BP",
    "SCHEMA_VERSION",
    "SELECTED_FOCAL_IDENTITY_SCHEMA",
    "NeutralCandidate",
    "SelectedFocalIdentityError",
    "SimulationArtifacts",
    "SimulationDrawTimeout",
    "SimulationExhausted",
    "UnitLockHeld",
    "candidate_pool_frequency_gate",
    "exclusive_unit_lock",
    "feasible_exact_panel_genotype_triples",
    "find_neutral_candidates",
    "inspect_selected_focal_identity",
    "load_cell_demography",
    "select_exact_sample_panel",
    "simulate_neutral_unit",
    "simulate_selected_unit",
    "simulate_unit",
    "simulate_units",
    "unit_lock_is_held",
    "validate_selected_focal_identity",
]
