"""Direct, bounded Han selected-production executor (v3).

This additive module replaces the failed fixed-v2 *rate gate* with the direct
operational estimand: ten accepted, production-length selected simulations in
each prespecified s-by-AF cell.  It does not edit or authorize any neutral
simulation and does not use Gamma-SMC output for planning, acceptance, or
completion.

The first production replicate in every cell is a six-unit pilot under the
exact final contract.  Those accepted units count toward the final sixty.  A
separate ``simulate`` stage runs replicates 2--10 only after the pilot audit
passes.  Every unit has exactly 300 prospectively allocated SLiM seeds and 300
distinct panel seeds, a five-minute per-draw timeout, and a 25-hour cumulative
budget.  Exhaustion is a terminal failure, never a silently reduced sample.
"""

from __future__ import annotations

import argparse
import ast
import ctypes
import hashlib
import json
import math
import multiprocessing
import os
import platform
import re
import select
import shutil
import signal
import subprocess
import sys
import tempfile
import threading
import traceback
from collections.abc import Callable, Iterator, Mapping, Sequence
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from time import monotonic, perf_counter, sleep, time_ns
from typing import Any

import msprime
import numpy as np
import pandas as pd
import pyslim
import stdpopsim
import tskit

from . import focused_han_fixed_v2_calibration as fixed_v2
from . import focused_selection_simulation as simulation
from .eas_sweep_study import BASE_SEED

MANIFEST_SCHEMA = "gamma-smc.han-direct-production-v3-plan/v1"
AUTHORIZATION_SCHEMA = "gamma-smc.han-direct-production-v3-authorization/v1"
STAGE_AUDIT_SCHEMA = "gamma-smc.han-direct-production-v3-stage-audit/v1"
QUARANTINE_PLAN_SCHEMA = "gamma-smc.han-direct-production-v3-quarantine-plan/v1"
ATTEMPT_LEDGER_SCHEMA = "gamma-smc.han-direct-production-v3-attempt-ledger/v1"
SUPERVISOR_LEDGER_SCHEMA = "gamma-smc.han-direct-production-v3-process-supervisor/v1"
SUPERVISOR_ACTIVE_SCHEMA = "gamma-smc.han-direct-production-v3-active-supervision/v1"
SUPERVISOR_RECONCILIATION_SCHEMA = (
    "gamma-smc.han-direct-production-v3-supervision-reconciliation/v1"
)
SUPERVISOR_STARTUP_CLEANUP_SCHEMA = (
    "gamma-smc.han-direct-production-v3-startup-cleanup/v1"
)
UNIT_ATOMIC_TARGET_NAMES = (
    "attempts.tsv",
    "supervision.tsv",
    "supervision_reconciliation.tsv",
    "supervision_startup_cleanup.tsv",
    "simulation_contract.json",
    "simulation_failed.json",
    "simulation.trees",
    "sample_manifest.tsv",
    "truth_profiles.tsv.gz",
    "truth_class_summaries.tsv",
    "simulation_complete.json",
)
EXCLUDED_SEED_REGISTRY_SCHEMA = (
    "gamma-smc.han-direct-production-v3-excluded-seed-registry/v1"
)
EXCLUDED_SEED_INVENTORY_SCHEMA = (
    "gamma-smc.han-direct-production-v3-excluded-seed-inventory/v1"
)

MODULE_SOURCE_PATH = "python/gamma_smc_aou/focused_han_direct_v3.py"
WRAPPER_SOURCE_PATH = "scripts/run_focused_han_direct_v3.py"
BOOTSTRAP_SOURCE_PATH = "scripts/bootstrap_focused_han_direct_v3.sh"
DEFAULT_CAMPAIGN_RELATIVE_PATH = "focused_selection_EAS_sim"
DEFAULT_OUTPUT_RELATIVE_PATH = (
    "focused_selection_EAS_sim/calibration/han_direct_production_v3"
)
DEFAULT_PLAN_MANIFEST_NAME = "han_direct_v3_manifest.json"
DEFAULT_SEED_PLAN_NAME = "han_direct_v3_attempt_seed_plan.tsv"
DEFAULT_EXCLUDED_SEED_REGISTRY_NAME = "han_direct_v3_excluded_seed_registry.tsv"
DEFAULT_EXCLUDED_SEED_INVENTORY_NAME = "han_direct_v3_excluded_seed_inventory.json"
FAILED_V2_QUARANTINE_ROOT = (
    "focused_selection_EAS_sim/calibration/han_fixed_v2_retry_quarantine/"
    "20260813T113900-0500_panel_seed_result_contract"
)
FAILED_V2_PILOT_RELATIVE_PATH = (
    f"{FAILED_V2_QUARANTINE_ROOT}/planning_bundle/confirmation_run_snapshot.tsv"
)
FAILED_V2_COMPACT_SUMMARY_RELATIVE_PATH = (
    f"{FAILED_V2_QUARANTINE_ROOT}/quarantine_inventory.json"
)

# This source-controlled allowlist is the sole planning-time seed namespace.
# New local files are never recursively discovered.  Optional local evidence
# may be checked, but its seeds must already be a subset of this registry.
CANONICAL_SEED_ARTIFACT_ALLOWLIST = (
    "focused_selection_EAS_sim/execution_units.tsv",
    "focused_selection_EAS_sim/calibration/eas_selected/screen20_plan.tsv",
    "focused_selection_EAS_sim/calibration/eas_selected/confirm100_plan.tsv",
    "focused_selection_EAS_sim/calibration/eas_selected/sensitivity10mb20_plan.tsv",
    (
        "focused_selection_EAS_sim/calibration/exploratory_seed_provenance/"
        "exploratory_seed_reservation.tsv"
    ),
    (
        "focused_selection_EAS_sim/calibration/exploratory_seed_provenance/"
        "previously_omitted_seed_reservation.tsv"
    ),
    (
        "focused_selection_EAS_sim/calibration/han_fixed_cessation_v2/"
        "han_fixed_v2_all_phase_plan.tsv"
    ),
    "focused_selection_EAS_sim/calibration/han_frequency_conditioned_all_phase_plan.tsv",
    "focused_selection_EAS_sim/calibration/han_selection_end_plan.tsv",
    "focused_selection_EAS_sim/calibration/han_selection_end_extension_plan.tsv",
    "focused_selection_EAS_sim/calibration/han_selection_end_reconciliation_plan.tsv",
    "focused_selection_EAS_sim/calibration/han_selection_end_terminal_draws_reconciled.tsv",
    (f"{FAILED_V2_QUARANTINE_ROOT}/planning_bundle/han_fixed_v2_all_phase_plan.tsv"),
)

SELECTION_COEFFICIENTS = (0.005, 0.01)
TARGET_FREQUENCIES = (0.10, 0.20, 0.30)
REPLICATES_PER_CELL = 10
ATTEMPTS_PER_UNIT = 300
CANDIDATE_POOL_DIPLOIDS = 500
SAMPLE_DIPLOIDS = 100
AF_HALF_WIDTH = 0.025
SLIM_SCALING_FACTOR = 5.0
SLIM_BURN_IN = 0.1
DRAW_TIMEOUT_SECONDS = 5.0 * 60.0
CUMULATIVE_TIMEOUT_SECONDS = ATTEMPTS_PER_UNIT * DRAW_TIMEOUT_SECONDS
SUPERVISOR_TERM_GRACE_SECONDS = 5.0
SUPERVISOR_POLL_SECONDS = 0.02
EXPECTED_SLIM_VERSION = "4.2.2"
EXPECTED_SLIM_BUILD = "h5888daf_1"
EXPECTED_SLIM_PACKAGE_SHA256 = (
    "d9ba79f22fe916bb5f76eb66284e6ec45dcb0e8694319a9f0112e07f9145aec3"
)
EXPECTED_SLIM_BINARY_SHA256 = (
    "67409ff190808c2967c949cae6697bf235010dbe832389cb4976d9089b5f5e1d"
)
EXPECTED_RUNTIME_VERSIONS = {
    "python": "3.11.15",
    "msprime": "1.4.2",
    "numpy": "2.4.6",
    "pandas": "3.0.3",
    "pyslim": "1.1.1",
    "stdpopsim": "0.3.0",
    "tskit": "1.0.3",
}
PILOT_REPLICATE_INDEX = 1
STAGE_PILOT = "pilot"
STAGE_SIMULATE = "simulate"
STAGE_ALL = "all"

PLAN_COLUMNS = (
    "unit_id",
    "stage",
    "selection_coefficient",
    "target_allele_frequency",
    "replicate_index",
    "production_unit_seed",
    "attempt_zero_based",
    "simulation_seed",
    "simulation_seed_nonce",
    "panel_seed",
    "panel_seed_nonce",
    "realized_selection_end_generations_ago",
    "duration_multiplier",
    "sequence_length_bp",
    "focal_position_bp",
    "candidate_pool_diploids",
    "sample_diploids",
    "draw_timeout_seconds",
    "cumulative_timeout_seconds",
)

SOURCE_PATHS = (
    ".gitattributes",
    MODULE_SOURCE_PATH,
    WRAPPER_SOURCE_PATH,
    BOOTSTRAP_SOURCE_PATH,
    "scripts/bootstrap_eas_sweep_study.sh",
    "scripts/bootstrap_uv.sh",
    "tests/test_focused_han_direct_v3.py",
    simulation.SIMULATION_SOURCE_PATH,
    "python/gamma_smc_aou/focused_han_fixed_v2_calibration.py",
    "python/gamma_smc_aou/eas_sweep_models.py",
    "python/gamma_smc_aou/eas_sweep_analysis.py",
    "python/gamma_smc_aou/eas_sweep_study.py",
    "pyproject.toml",
    "uv.lock",
)

_PATCH_LOCK = threading.RLock()
_WORKER_BUNDLE: DirectV3Bundle | None = None
# Cache identity must not depend on another additive executor's temporary
# process-local replacement of simulation.IMPLEMENTATION_PATHS.
_BASE_SIMULATION_IMPLEMENTATION_PATHS = tuple(simulation.IMPLEMENTATION_PATHS)


@dataclass(frozen=True)
class DirectV3Bundle:
    repo_root: Path
    campaign_dir: Path
    output_dir: Path
    manifest_path: Path
    seed_plan_path: Path
    excluded_seed_registry_path: Path
    excluded_seed_inventory_path: Path
    slim_path: Path
    manifest: Mapping[str, Any]
    seed_plan: pd.DataFrame
    excluded_seed_registry: pd.DataFrame
    excluded_seed_inventory: Mapping[str, Any]
    units: pd.DataFrame


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
        ).encode("utf-8")
    ).hexdigest()


def _seed_digest(values: Sequence[int] | set[int] | pd.Series) -> str:
    normalized = sorted(int(value) for value in values)
    return hashlib.sha256(
        json.dumps(normalized, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _inside_repo(path: str | Path, repo_root: Path, *, label: str) -> Path:
    candidate = Path(path)
    if not candidate.is_absolute():
        candidate = repo_root / candidate
    resolved = candidate.resolve()
    try:
        resolved.relative_to(repo_root)
    except ValueError as error:
        raise ValueError(f"Han v3 {label} must remain inside repo_root") from error
    if any("onedrive" in part.casefold() for part in resolved.parts):
        raise ValueError(f"Han v3 {label} must not be inside OneDrive")
    return resolved


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


def _is_unit_atomic_temp_name(name: str) -> bool:
    return (
        any(
            re.fullmatch(re.escape(target) + r"\.tmp\.[0-9]+", name)
            for target in UNIT_ATOMIC_TARGET_NAMES
        )
        or re.fullmatch(r"\.han_v3_active_[0-9]{3}\.json\.tmp\.[0-9]+", name)
        is not None
    )


def _unit_atomic_temp_paths(unit_dir: Path) -> list[Path]:
    paths: list[Path] = []
    for target in UNIT_ATOMIC_TARGET_NAMES:
        for candidate in unit_dir.glob(f"{target}.tmp.*"):
            suffix = candidate.name.removeprefix(f"{target}.tmp.")
            if (
                candidate.is_file()
                and re.fullmatch(r"[0-9]+", suffix)
                and _is_unit_atomic_temp_name(candidate.name)
            ):
                paths.append(candidate)
    for candidate in unit_dir.glob(".han_v3_active_*.json.tmp.*"):
        if candidate.is_file() and _is_unit_atomic_temp_name(candidate.name):
            paths.append(candidate)
    return sorted(paths, key=lambda path: path.name)


def _immutable_json(path: Path, payload: Any) -> None:
    if path.is_file():
        current = json.loads(path.read_text(encoding="utf-8"))
        if current != payload:
            raise ValueError(f"immutable Han v3 JSON differs: {path}")
        return
    _atomic_json(path, payload)


def _source_records(repo_root: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    for relative in SOURCE_PATHS:
        path = _inside_repo(repo_root / relative, repo_root, label="source")
        if not path.is_file():
            raise ValueError(f"Han v3 source is absent: {relative}")
        records[relative] = sha256_file(path)
    return records


def _slim_version(slim_path: Path) -> str:
    if not slim_path.is_file():
        raise ValueError(f"Han v3 SLiM executable is absent: {slim_path}")
    try:
        completed = subprocess.run(
            [str(slim_path), "-v"],
            check=True,
            capture_output=True,
            text=True,
            timeout=30,
        )
    except (OSError, subprocess.SubprocessError) as error:
        raise ValueError("Han v3 SLiM version probe failed") from error
    observed = " ".join((completed.stdout or completed.stderr).split())
    if "SLiM" not in observed or EXPECTED_SLIM_VERSION not in observed:
        raise ValueError(f"Han v3 requires SLiM {EXPECTED_SLIM_VERSION}")
    return observed


def _runtime_record(repo_root: Path, slim_path: Path) -> dict[str, Any]:
    if os.name != "posix" or platform.system() != "Linux":
        raise ValueError(
            "Han v3 planning/simulation requires the pinned Linux/WSL runtime "
            "used for process-group SLiM draw supervision"
        )
    if "fork" not in multiprocessing.get_all_start_methods():
        raise ValueError("Han v3 requires multiprocessing fork supervision")
    observed_versions = {
        "python": platform.python_version(),
        "msprime": msprime.__version__,
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "pyslim": pyslim.__version__,
        "stdpopsim": stdpopsim.__version__,
        "tskit": tskit.__version__,
    }
    if observed_versions != EXPECTED_RUNTIME_VERSIONS:
        raise ValueError(
            "Han v3 runtime versions differ from exact pins: "
            f"observed={observed_versions!r}"
        )
    slim = _inside_repo(slim_path, repo_root, label="SLiM executable")
    package_metadata = (
        slim.parent.parent
        / "conda-meta"
        / (f"slim-{EXPECTED_SLIM_VERSION}-{EXPECTED_SLIM_BUILD}.json")
    )
    if not package_metadata.is_file():
        raise ValueError("Han v3 pinned SLiM package metadata is absent")
    try:
        package = json.loads(package_metadata.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError("Han v3 pinned SLiM package metadata is unreadable") from error
    if (
        package.get("name") != "slim"
        or package.get("version") != EXPECTED_SLIM_VERSION
        or package.get("build") != EXPECTED_SLIM_BUILD
        or package.get("sha256") != EXPECTED_SLIM_PACKAGE_SHA256
        or sha256_file(slim) != EXPECTED_SLIM_BINARY_SHA256
    ):
        raise ValueError("Han v3 SLiM package/build/binary identity differs")
    executable = Path(sys.executable).resolve()
    if not executable.is_file():
        raise ValueError("Han v3 Python executable is absent")
    portable = {
        "os_name": os.name,
        "platform_system": platform.system(),
        "platform_machine": platform.machine(),
        "required_versions": EXPECTED_RUNTIME_VERSIONS,
        "observed_versions": observed_versions,
        "process_supervisor": {
            "schema": SUPERVISOR_LEDGER_SCHEMA,
            "multiprocessing_start_method": "fork",
            "child_session": "setsid",
            "linux_child_subreaper": True,
            "timeout_signal_sequence": "SIGTERM_then_grace_then_SIGKILL",
            "term_grace_seconds": SUPERVISOR_TERM_GRACE_SECONDS,
            "direct_and_adopted_children_reaped": True,
            "attempt_temp_tree_removed_before_return": True,
            "parent_death_watchdog": (
                "control_pipe_eof_keeps_session_leader_alive_to_terminate_and_"
                "reap_all_descendants_before_exit"
            ),
        },
        "slim": {
            "relative_path": slim.relative_to(repo_root).as_posix(),
            "binary_sha256": sha256_file(slim),
            "version": _slim_version(slim),
            "required_version": EXPECTED_SLIM_VERSION,
            "required_build": EXPECTED_SLIM_BUILD,
            "package_sha256": package["sha256"],
        },
    }
    return {
        "portable_contract": portable,
        "host_diagnostics_not_compared": {
            "platform_release": platform.release(),
            "libc": list(platform.libc_ver()),
            "python_full_version": sys.version,
            "python_executable_path": str(executable),
            "python_executable_sha256": sha256_file(executable),
            "slim_absolute_path": str(slim),
            "slim_package_metadata_path": str(package_metadata),
            "slim_package_url": package.get("url"),
        },
    }


def _endpoint_records() -> dict[tuple[float, float], dict[str, float]]:
    table = fixed_v2.endpoint_table()
    records: dict[tuple[float, float], dict[str, float]] = {}
    for raw in table.to_dict(orient="records"):
        key = (
            float(raw["selection_coefficient"]),
            float(raw["target_allele_frequency"]),
        )
        records[key] = {
            "duration_multiplier": float(raw["duration_multiplier"]),
            "requested_selection_end_generations_ago": float(
                raw["requested_selection_end_generations_ago"]
            ),
            "realized_selection_end_generations_ago": float(
                raw["realized_selection_end_generations_ago"]
            ),
        }
    if set(records) != {
        (coefficient, target)
        for coefficient in SELECTION_COEFFICIENTS
        for target in TARGET_FREQUENCIES
    }:
        raise ValueError("Han v3 fixed endpoint table differs from six-cell contract")
    return records


def _exact_integer_series(
    frame: pd.DataFrame,
    column: str,
    *,
    label: str,
    minimum: int = 0,
) -> pd.Series:
    numeric = pd.to_numeric(frame[column], errors="coerce")
    floating = numeric.astype(float)
    if (
        numeric.isna().any()
        or not np.isfinite(floating).all()
        or not np.equal(floating, np.rint(floating)).all()
        or (floating < minimum).any()
    ):
        raise ValueError(f"Han v3 {label} {column} must contain exact integers")
    return numeric.astype(np.int64)


def _validate_production_units(frame: pd.DataFrame) -> pd.DataFrame:
    required = {
        "unit_id",
        "demography_id",
        "demography_kind",
        "population",
        "source_model",
        "selection_origin",
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
    if required.difference(frame.columns):
        raise ValueError("Han v3 execution-unit table lacks required columns")

    integer_columns = {
        column: _exact_integer_series(frame, column, label="execution", minimum=minimum)
        for column, minimum in {
            "replicate_index": 1,
            "seed": 1,
            "exact_sample_alt_count": 0,
            "sample_diploids": 1,
            "sequence_length_bp": 1,
            "focal_position_bp": 0,
        }.items()
    }
    frame = frame.copy()
    for column, values in integer_columns.items():
        frame[column] = values
    units = frame[
        frame["demography_id"].astype(str).eq("ancient_eurasia_han_introgression")
        & frame["simulation_class"].astype(str).eq("selected")
    ].copy()
    if len(units) != 60 or units["unit_id"].duplicated().any():
        raise ValueError("Han v3 requires exactly 60 unique Han selected units")
    if units["seed"].astype(int).duplicated().any():
        raise ValueError("Han v3 production unit seeds are not unique")
    cell_counts = units.groupby(
        ["selection_coefficient", "target_allele_frequency"]
    ).size()
    replicate_contract = units.groupby(
        ["selection_coefficient", "target_allele_frequency"]
    )["replicate_index"].apply(lambda values: sorted(values.astype(int)))
    if (
        len(cell_counts) != 6
        or set(cell_counts.astype(int)) != {REPLICATES_PER_CELL}
        or set(units["selection_coefficient"].astype(float))
        != set(SELECTION_COEFFICIENTS)
        or set(units["target_allele_frequency"].astype(float))
        != set(TARGET_FREQUENCIES)
        or set(units["replicate_index"].astype(int))
        != set(range(1, REPLICATES_PER_CELL + 1))
        or any(
            values != list(range(1, REPLICATES_PER_CELL + 1))
            for values in replicate_contract
        )
    ):
        raise ValueError("Han v3 production six-cell/replicate grid differs")
    exact_source_fields = {
        "demography_kind": "introgression",
        "population": "Han",
        "source_model": "stdpopsim AncientEurasia_9K19",
        "selection_origin": "archaic_specific_introgressed_standing_variation",
    }
    for column, expected in exact_source_fields.items():
        if set(units[column].astype(str)) != {expected}:
            raise ValueError(f"Han v3 production {column} field differs")
    expected_unit_ids = []
    for raw in units.to_dict(orient="records"):
        coefficient = float(raw["selection_coefficient"])
        target = float(raw["target_allele_frequency"])
        coefficient_label = f"{coefficient:.3f}".replace(".", "p")
        target_label = round(100 * target)
        expected_unit_ids.append(
            "ancient_eurasia_han_introgression"
            f"__af{target_label}__selected_s{coefficient_label}"
            f"__rep{int(raw['replicate_index']):03d}"
        )
    if list(units["unit_id"].astype(str)) != expected_unit_ids:
        raise ValueError("Han v3 unit IDs do not correspond to AF/s/rep fields")
    expected_alt = np.rint(
        2
        * units["sample_diploids"].astype(int)
        * units["target_allele_frequency"].astype(float)
    ).astype(int)
    if (
        set(units["sample_diploids"].astype(int)) != {SAMPLE_DIPLOIDS}
        or not np.array_equal(
            units["exact_sample_alt_count"].astype(int).to_numpy(),
            expected_alt.to_numpy(),
        )
        or set(units["sequence_length_bp"].astype(int))
        != {simulation.SEQUENCE_LENGTH_BP}
        or set(units["focal_position_bp"].astype(int)) != {simulation.FOCAL_POSITION_BP}
    ):
        raise ValueError("Han v3 exact panel or 10-Mb geometry differs")
    if not np.allclose(
        units["population_af_lower"].astype(float),
        units["target_allele_frequency"].astype(float) - AF_HALF_WIDTH,
        rtol=0.0,
        atol=1e-12,
    ) or not np.allclose(
        units["population_af_upper"].astype(float),
        units["target_allele_frequency"].astype(float) + AF_HALF_WIDTH,
        rtol=0.0,
        atol=1e-12,
    ):
        raise ValueError("Han v3 candidate-pool AF gate differs from target +/- 0.025")
    return units.sort_values(
        ["selection_coefficient", "target_allele_frequency", "replicate_index"]
    ).reset_index(drop=True)


def _seed_columns(frame: pd.DataFrame) -> list[str]:
    result = []
    for column in frame.columns:
        normalized = str(column).casefold()
        if "seed" not in normalized or any(
            token in normalized for token in ("sha", "digest", "nonce", "count")
        ):
            continue
        if (
            normalized == "seed"
            or normalized.endswith("_seed")
            or normalized.startswith("seed_")
        ):
            result.append(str(column))
    return result


def _numeric_seeds(frame: pd.DataFrame) -> set[int]:
    values: set[int] = set()
    for column in _seed_columns(frame):
        numeric = pd.to_numeric(frame[column], errors="coerce").dropna()
        floating = numeric.astype(float)
        if (
            not np.isfinite(floating).all()
            or not np.equal(floating, np.rint(floating)).all()
        ):
            raise ValueError(f"Han v3 seed artifact {column} is not integral")
        for raw in numeric:
            value = int(raw)
            if value >= 1:
                values.add(value)
    return values


def _classify_seed_artifact(relative: str) -> str:
    lowered = relative.casefold()
    if lowered.endswith("execution_units.tsv"):
        return "production"
    if "exploratory_seed_provenance" in lowered:
        return "exploratory"
    if "eas_selected" in lowered:
        return "eas"
    if "han_fixed" in lowered or "fixed_cessation_v2" in lowered:
        return "v2"
    if "han_selection_end" in lowered or "han_frequency_conditioned" in lowered:
        return "v1"
    return "fallback"


def _require_git_tracked(paths: Sequence[Path], repo_root: Path) -> None:
    """Fail unless every canonical seed artifact is tracked by this checkout."""

    relatives = [path.relative_to(repo_root).as_posix() for path in paths]
    completed = subprocess.run(
        ["git", "-C", str(repo_root), "ls-files", "--error-unmatch", "--", *relatives],
        check=False,
        capture_output=True,
        text=True,
        timeout=60,
    )
    observed = {line.strip() for line in completed.stdout.splitlines() if line.strip()}
    if completed.returncode != 0 or observed != set(relatives):
        missing = sorted(set(relatives).difference(observed))
        raise ValueError(
            "Han v3 canonical seed artifacts must all be tracked: " + ", ".join(missing)
        )


def canonical_seed_artifacts(*, repo_root: str | Path) -> list[Path]:
    """Return the fixed, source-controlled seed-exclusion allowlist."""

    root = Path(repo_root).resolve()
    paths = [
        _inside_repo(root / relative, root, label="canonical seed artifact")
        for relative in CANONICAL_SEED_ARTIFACT_ALLOWLIST
    ]
    if len(paths) != len(set(paths)) or any(not path.is_file() for path in paths):
        raise ValueError(
            "Han v3 canonical seed-artifact allowlist is absent/duplicated"
        )
    _require_git_tracked(paths, root)
    return paths


def _validate_local_seed_subset(
    paths: Sequence[str | Path],
    *,
    repo_root: Path,
    excluded: set[int],
) -> list[dict[str, Any]]:
    """Allow local evidence only when it introduces no seed outside the registry."""

    records: list[dict[str, Any]] = []
    observed: set[Path] = set()
    for raw in paths:
        path = _inside_repo(raw, repo_root, label="local seed evidence")
        if path in observed or not path.is_file():
            raise ValueError("Han v3 local seed evidence is absent or duplicated")
        observed.add(path)
        frame = pd.read_csv(path, sep="\t")
        seeds = _numeric_seeds(frame)
        if not seeds or not seeds.issubset(excluded):
            raise ValueError(
                "Han v3 local seed evidence must be a nonempty subset of the "
                "canonical excluded-seed registry"
            )
        records.append(
            {
                "path": path.relative_to(repo_root).as_posix(),
                "sha256": sha256_file(path),
                "rows": len(frame),
                "seed_columns": _seed_columns(frame),
                "direct_unique_count": len(seeds),
                "direct_seed_sha256": _seed_digest(seeds),
                "subset_of_frozen_registry": True,
                "planning_authority": False,
            }
        )
    return sorted(records, key=lambda item: item["path"])


def _expand_production_fallback_namespace(frame: pd.DataFrame) -> set[int]:
    """Reserve existing unit seeds and all bounded fallback derivations.

    Selected bases are conservatively expanded through 300 attempts, which
    covers both the old 50/68/100-draw paths and any direct retry that might
    otherwise reuse a familiar ``selected:*`` label.  Neutral expansions match
    their unchanged bounded executor.
    """

    occupied = {int(value) for value in frame["seed"]}
    for row in frame.to_dict(orient="records"):
        base = int(row["seed"])
        if str(row["simulation_class"]) == "selected":
            for attempt in range(ATTEMPTS_PER_UNIT):
                occupied.add(simulation._stable_seed(base, f"selected:{attempt}"))
                occupied.add(simulation._stable_seed(base, f"selected-panel:{attempt}"))
        elif str(row["simulation_class"]) == "neutral":
            for attempt in range(simulation.DEFAULT_NEUTRAL_MAX_ATTEMPTS):
                for component in ("ancestry", "mutations", "panel"):
                    occupied.add(
                        simulation._stable_seed(base, f"neutral:{attempt}:{component}")
                    )
            batches = math.ceil(
                simulation.DEFAULT_NEUTRAL_MAX_ATTEMPTS
                / simulation.DEFAULT_NEUTRAL_ANCESTRY_BATCH_SIZE
            )
            for batch in range(batches):
                occupied.add(
                    simulation._stable_seed(base, f"neutral-batch:{batch}:candidate")
                )
        else:
            raise ValueError("Han v3 production table contains unknown class")
    return occupied


def collect_seed_exclusions(
    paths: Sequence[str | Path], *, repo_root: str | Path
) -> tuple[set[int], list[dict[str, Any]], dict[str, set[int]]]:
    root = Path(repo_root).resolve()
    if not paths:
        raise ValueError("Han v3 seed-exclusion artifacts are required")
    occupied: set[int] = set()
    categories: dict[str, set[int]] = {
        key: set()
        for key in ("production", "fallback", "v1", "v2", "eas", "exploratory")
    }
    records: list[dict[str, Any]] = []
    observed: set[str] = set()
    for raw in paths:
        path = _inside_repo(raw, root, label="seed artifact")
        relative = path.relative_to(root).as_posix()
        if relative in observed or not path.is_file():
            raise ValueError(f"Han v3 seed artifact is invalid: {relative}")
        observed.add(relative)
        frame = pd.read_csv(path, sep="\t")
        direct = _numeric_seeds(frame)
        if not direct:
            raise ValueError(f"Han v3 seed artifact has no numeric seeds: {relative}")
        category = _classify_seed_artifact(relative)
        expanded = set(direct)
        expansion = "explicit_seed_columns"
        if path.name == "execution_units.tsv":
            expanded = _expand_production_fallback_namespace(frame)
            categories["production"].update(direct)
            categories["fallback"].update(expanded.difference(direct))
            expansion = "production_bases_plus_300_draw_selected_fallback_namespace"
        else:
            categories[category].update(expanded)
            if "han_selection_end" in relative.casefold():
                derived = {
                    simulation._stable_seed(seed, "sensitivity-exact-panel")
                    for seed in direct
                }
                expanded.update(derived)
                categories[category].update(derived)
                expansion += "+legacy_sensitivity_panel_derivation"
        occupied.update(expanded)
        records.append(
            {
                "path": relative,
                "sha256": sha256_file(path),
                "category": category,
                "rows": len(frame),
                "seed_columns": _seed_columns(frame),
                "direct_unique_count": len(direct),
                "direct_seed_sha256": _seed_digest(direct),
                "expansion": expansion,
                "reserved_unique_count": len(expanded),
                "reserved_seed_sha256": _seed_digest(expanded),
            }
        )
    absent = [name for name, seeds in categories.items() if not seeds]
    if absent:
        raise ValueError(
            "Han v3 seed exclusions lack required namespaces: " + ", ".join(absent)
        )
    return occupied, sorted(records, key=lambda item: item["path"]), categories


def _excluded_seed_registry_frame(
    excluded: set[int], categories: Mapping[str, set[int]]
) -> pd.DataFrame:
    rows = []
    for seed in sorted(excluded):
        namespaces = sorted(
            name for name, values in categories.items() if int(seed) in values
        )
        if not namespaces:
            raise ValueError("Han v3 excluded seed lacks a namespace")
        rows.append(
            {
                "schema": EXCLUDED_SEED_REGISTRY_SCHEMA,
                "seed": int(seed),
                "namespaces": ",".join(namespaces),
            }
        )
    frame = pd.DataFrame(rows, columns=("schema", "seed", "namespaces"))
    if frame.empty or frame["seed"].duplicated().any():
        raise ValueError("Han v3 excluded-seed registry is empty or duplicated")
    return frame


def _validate_frozen_seed_registry(
    frame: pd.DataFrame, expected: Mapping[str, Any]
) -> set[int]:
    if list(frame.columns) != ["schema", "seed", "namespaces"] or frame.empty:
        raise ValueError("Han v3 frozen excluded-seed registry shape differs")
    if set(frame["schema"].astype(str)) != {EXCLUDED_SEED_REGISTRY_SCHEMA}:
        raise ValueError("Han v3 frozen excluded-seed registry schema differs")
    seeds = _exact_integer_series(
        frame, "seed", label="frozen excluded-seed registry", minimum=1
    )
    if seeds.duplicated().any() or list(seeds) != sorted(seeds):
        raise ValueError("Han v3 frozen excluded-seed registry is not canonical")
    namespace_sets = (
        frame["namespaces"]
        .astype(str)
        .map(lambda value: tuple(sorted(filter(None, value.split(",")))))
    )
    allowed = {"production", "fallback", "v1", "v2", "eas", "exploratory"}
    if any(
        not values or not set(values).issubset(allowed) for values in namespace_sets
    ):
        raise ValueError("Han v3 frozen excluded-seed namespaces differ")
    result = set(seeds)
    if len(result) != int(expected["excluded_union_count"]) or _seed_digest(
        result
    ) != str(expected["excluded_union_sha256"]):
        raise ValueError("Han v3 frozen excluded-seed registry digest differs")
    observed_namespaces = {}
    for name in sorted(allowed):
        values = {
            int(seed)
            for seed, namespaces in zip(seeds, namespace_sets, strict=True)
            if name in namespaces
        }
        observed_namespaces[name] = {
            "count": len(values),
            "seed_sha256": _seed_digest(values),
        }
    if observed_namespaces != expected["excluded_namespaces"]:
        raise ValueError("Han v3 frozen excluded-seed namespace digest differs")
    return result


def _allocate_seed(identity: str, occupied: set[int]) -> tuple[int, int]:
    for nonce in range(100_000):
        label = identity if nonce == 0 else f"{identity}:nonce={nonce}"
        value = simulation._stable_seed(BASE_SEED, label)
        if value not in occupied:
            occupied.add(value)
            return value, nonce
    raise RuntimeError("Han v3 could not allocate a disjoint RNG seed")


def _failed_v2_budget_evidence(
    snapshot_path: Path, summary_path: Path, repo_root: Path
) -> dict[str, Any]:
    artifact = _inside_repo(
        snapshot_path, repo_root, label="failed-v2 pilot snapshot evidence"
    )
    summary_artifact = _inside_repo(
        summary_path, repo_root, label="failed-v2 compact summary evidence"
    )
    if not artifact.is_file() or not summary_artifact.is_file():
        raise ValueError("Han v3 failed-v2 budget evidence is absent")
    _require_git_tracked((artifact, summary_artifact), repo_root)
    frame = pd.read_csv(artifact, sep="\t")
    if frame.empty or not {"status", "elapsed_seconds"}.issubset(frame.columns):
        raise ValueError("Han v3 failed-v2 budget evidence is malformed")
    elapsed = pd.to_numeric(frame["elapsed_seconds"], errors="coerce").dropna()
    statuses = frame["status"].fillna("").astype(str).value_counts().to_dict()
    if int(statuses.get("failed", 0)) < 1:
        raise ValueError("Han v3 pilot evidence is not from the failed v2 run")
    try:
        summary = json.loads(summary_artifact.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError("Han v3 failed-v2 compact summary is malformed") from error
    work = summary.get("work_evidence", {})
    completed_rows = int(statuses.get("complete", 0))
    direct_failure_rows = int(
        frame["error"]
        .fillna("")
        .astype(str)
        .str.contains("exited with code -15", regex=False)
        .sum()
    )
    if (
        summary.get("schema") != "gamma-smc.han-fixed-cessation-v2-retry-quarantine/v1"
        or summary.get("status") != "canonical_noninferential_retry_provenance"
        or summary.get("inference_use") is not False
        or completed_rows != int(work.get("valid_completions", -1))
        or direct_failure_rows != int(work.get("sigterm_failure_records", -1))
        or int(work.get("started_directories", -1))
        != completed_rows
        + int(work.get("failure_records", -1))
        + int(work.get("partial_directories", -1))
    ):
        raise ValueError("Han v3 failed-v2 compact summary semantics differ")
    return {
        "snapshot": {
            "path": artifact.relative_to(repo_root).as_posix(),
            "sha256": sha256_file(artifact),
            "rows": len(frame),
        },
        "compact_summary": {
            "path": summary_artifact.relative_to(repo_root).as_posix(),
            "sha256": sha256_file(summary_artifact),
            "schema": summary["schema"],
            "status": summary["status"],
            "work_evidence": work,
        },
        "status_counts": {str(key): int(value) for key, value in statuses.items()},
        "direct_sigterm_failure_rows": direct_failure_rows,
        "completed_elapsed_seconds": {
            "n": len(elapsed),
            "mean": float(elapsed.mean()) if len(elapsed) else None,
            "maximum": float(elapsed.max()) if len(elapsed) else None,
        },
        "use": "pilot_compute_budget_evidence_only",
        "inferential_use": False,
        "endpoint_authorization": False,
        "production_unit_acceptance_use": False,
        "gamma_smc_statistics_used": False,
    }


def build_plan(
    *,
    repo_root: str | Path,
    campaign_dir: str | Path,
    slim_path: str | Path,
    local_seed_artifacts: Sequence[str | Path] = (),
    failed_v2_evidence_path: str | Path = FAILED_V2_PILOT_RELATIVE_PATH,
    failed_v2_summary_path: str | Path = FAILED_V2_COMPACT_SUMMARY_RELATIVE_PATH,
) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    root = Path(repo_root).resolve()
    if any("onedrive" in part.casefold() for part in root.parts):
        raise ValueError("Han v3 repo_root must not be inside OneDrive")
    campaign = _inside_repo(campaign_dir, root, label="campaign directory")
    slim = _inside_repo(slim_path, root, label="SLiM executable")
    execution_path = campaign / "execution_units.tsv"
    execution = pd.read_csv(execution_path, sep="\t")
    units = _validate_production_units(execution)
    paths = canonical_seed_artifacts(repo_root=root)
    excluded, evidence, categories = collect_seed_exclusions(paths, repo_root=root)
    local_seed_evidence = _validate_local_seed_subset(
        local_seed_artifacts, repo_root=root, excluded=excluded
    )
    excluded_registry = _excluded_seed_registry_frame(excluded, categories)
    excluded_inventory: dict[str, Any] = {
        "schema": EXCLUDED_SEED_INVENTORY_SCHEMA,
        "status": "frozen_at_atomic_plan_publication",
        "canonical_allowlist_source": MODULE_SOURCE_PATH,
        "canonical_allowlist": list(CANONICAL_SEED_ARTIFACT_ALLOWLIST),
        "canonical_artifacts": evidence,
        "local_subset_evidence": local_seed_evidence,
        "local_evidence_planning_authority": False,
        "excluded_union_count": len(excluded),
        "excluded_union_sha256": _seed_digest(excluded),
    }
    excluded_inventory["payload_sha256"] = _canonical_sha256(excluded_inventory)
    occupied = set(excluded)
    endpoints = _endpoint_records()
    rows: list[dict[str, Any]] = []
    for unit in units.to_dict(orient="records"):
        coefficient = float(unit["selection_coefficient"])
        target = float(unit["target_allele_frequency"])
        endpoint = endpoints[(coefficient, target)]
        stage = (
            STAGE_PILOT
            if int(unit["replicate_index"]) == PILOT_REPLICATE_INDEX
            else STAGE_SIMULATE
        )
        for attempt in range(ATTEMPTS_PER_UNIT):
            identity = (
                f"han_direct_v3:{unit['unit_id']}:attempt={attempt:03d}:simulation"
            )
            sim_seed, sim_nonce = _allocate_seed(identity, occupied)
            panel_seed, panel_nonce = _allocate_seed(
                f"han_direct_v3:{unit['unit_id']}:attempt={attempt:03d}:panel",
                occupied,
            )
            rows.append(
                {
                    "unit_id": str(unit["unit_id"]),
                    "stage": stage,
                    "selection_coefficient": coefficient,
                    "target_allele_frequency": target,
                    "replicate_index": int(unit["replicate_index"]),
                    "production_unit_seed": int(unit["seed"]),
                    "attempt_zero_based": attempt,
                    "simulation_seed": sim_seed,
                    "simulation_seed_nonce": sim_nonce,
                    "panel_seed": panel_seed,
                    "panel_seed_nonce": panel_nonce,
                    "realized_selection_end_generations_ago": endpoint[
                        "realized_selection_end_generations_ago"
                    ],
                    "duration_multiplier": endpoint["duration_multiplier"],
                    "sequence_length_bp": simulation.SEQUENCE_LENGTH_BP,
                    "focal_position_bp": simulation.FOCAL_POSITION_BP,
                    "candidate_pool_diploids": CANDIDATE_POOL_DIPLOIDS,
                    "sample_diploids": SAMPLE_DIPLOIDS,
                    "draw_timeout_seconds": DRAW_TIMEOUT_SECONDS,
                    "cumulative_timeout_seconds": CUMULATIVE_TIMEOUT_SECONDS,
                }
            )
    plan = pd.DataFrame(rows, columns=PLAN_COLUMNS)
    _validate_seed_plan(plan, units=units, excluded=excluded)
    all_v3_seeds = pd.concat(
        [plan["simulation_seed"], plan["panel_seed"]], ignore_index=True
    )
    namespace_records = {
        name: {"count": len(values), "seed_sha256": _seed_digest(values)}
        for name, values in categories.items()
    }
    manifest: dict[str, Any] = {
        "schema": MANIFEST_SCHEMA,
        "status": "planned",
        "model": {
            "demography_id": "ancient_eurasia_han_introgression",
            "source_model": "stdpopsim AncientEurasia_9K19",
            "mode": "direct_operational_estimand_fixed_cessation_q5",
            "selection_coefficients": list(SELECTION_COEFFICIENTS),
            "target_allele_frequencies": list(TARGET_FREQUENCIES),
            "fixed_endpoints": [
                {
                    "selection_coefficient": coefficient,
                    "target_allele_frequency": target,
                    **record,
                }
                for (coefficient, target), record in sorted(endpoints.items())
            ],
            "mutation_origin_population": "Neanderthal",
            "mutation_age_generations": 2400.0,
            "introgression_pulse_generations": 2272.0,
            "archaic_specific_no_ils_by_construction": True,
            "post_pulse_recipient_nonloss_conditioning": True,
            "terminal_population_af_conditioning_in_slim": False,
            "external_candidate_pool_af_gate": "target_plus_or_minus_0.025_inclusive",
            "uniform_exact_k_panel": True,
            "strict_selected_type1_identity": True,
            "normal_stdpopsim_processing": (
                "10Mb_forward_then_recapitation_and_neutral_mutation_overlay"
            ),
            "gamma_smc_statistics_used": False,
        },
        "design": {
            "units": 60,
            "cells": 6,
            "replicates_per_cell": REPLICATES_PER_CELL,
            "pilot_units": 6,
            "remaining_simulation_units": 54,
            "pilot_replicate_index": PILOT_REPLICATE_INDEX,
            "pilot_units_count_toward_final": True,
            "attempts_per_unit": ATTEMPTS_PER_UNIT,
            "simulation_seeds_per_unit": ATTEMPTS_PER_UNIT,
            "panel_seeds_per_unit": ATTEMPTS_PER_UNIT,
            "draw_timeout_seconds": DRAW_TIMEOUT_SECONDS,
            "cumulative_timeout_seconds": CUMULATIVE_TIMEOUT_SECONDS,
            "fail_closed_on_any_unit_exhaustion": True,
            "sequence_length_bp": simulation.SEQUENCE_LENGTH_BP,
            "focal_position_bp": simulation.FOCAL_POSITION_BP,
            "candidate_pool_diploids": CANDIDATE_POOL_DIPLOIDS,
            "sample_diploids": SAMPLE_DIPLOIDS,
            "slim_scaling_factor": SLIM_SCALING_FACTOR,
            "slim_burn_in": SLIM_BURN_IN,
        },
        "production_execution_units": {
            "path": execution_path.relative_to(root).as_posix(),
            "sha256": sha256_file(execution_path),
            "all_rows": len(execution),
            "han_selected_rows": len(units),
            "unit_id_sha256": _canonical_sha256(sorted(units["unit_id"].astype(str))),
            "production_seed_sha256": _seed_digest(units["seed"]),
        },
        "failed_v2_pilot_evidence": _failed_v2_budget_evidence(
            _inside_repo(failed_v2_evidence_path, root, label="v2 evidence"),
            _inside_repo(failed_v2_summary_path, root, label="v2 compact summary"),
            root,
        ),
        "seed_registry": {
            "base_seed": BASE_SEED,
            "all_attempt_and_panel_seeds_allocated_atomically": True,
            "attempt_seed_count": len(plan),
            "attempt_seed_sha256": _seed_digest(plan["simulation_seed"]),
            "panel_seed_count": len(plan),
            "panel_seed_sha256": _seed_digest(plan["panel_seed"]),
            "v3_union_count": len(all_v3_seeds),
            "v3_union_sha256": _seed_digest(all_v3_seeds),
            "v3_pairwise_distinct": True,
            "excluded_union_count": len(excluded),
            "excluded_union_sha256": _seed_digest(excluded),
            "excluded_namespaces": namespace_records,
            "disjoint_from_all_excluded_namespaces": True,
        },
        "seed_exclusion_policy": {
            "canonical_allowlist_is_source_controlled": True,
            "recursive_discovery_used": False,
            "local_seed_artifacts_must_be_subset": True,
            "load_validates_frozen_bundle_only": True,
            "canonical_artifact_count": len(evidence),
            "local_subset_artifact_count": len(local_seed_evidence),
        },
        "implementation_sources": _source_records(root),
        "runtime": _runtime_record(root, slim),
        "seed_plan_rows_sha256": _canonical_sha256(plan.to_dict(orient="records")),
        "inference_contract": {
            "estimand": (
                "Gamma-SMC performance conditional on a valid type-1 selected "
                "trajectory reaching the prespecified 500-diploid AF band and "
                "then an exact-k 100-diploid genotype-QC panel"
            ),
            "endpoint_rule_selected_after_exploratory_screen": True,
            "exploratory_and_failed_v2_runs_inferential_use": False,
            "pilot_is_final_contract_and_counts_toward_sixty": True,
            "all_sixty_accepted_units_required": True,
            "neutral_rarity_or_unconditional_selection_probability_estimated": False,
            "gamma_smc_statistics_used_for_simulation_acceptance": False,
        },
        "scientific_caveats": [
            (
                "Inference is conditional on survival of the introgressed type-1 "
                "allele, the completed 500-diploid AF-band gate, exact sampled "
                "allele frequency, and genotype QC."
            ),
            (
                "The six fixed cessation endpoints were selected after an "
                "exploratory screen; v3 estimates method performance at those "
                "operational endpoints and does not validate their selection rule."
            ),
            (
                "SLiM Q=5 is a scaled approximation; all v3 pilot and final units "
                "use the identical Q=5, 10-Mb production contract."
            ),
            (
                "The failed v2 snapshot informs compute budgeting only and supplies "
                "no endpoint, trajectory, TMRCA, or Gamma-SMC evidence."
            ),
            (
                "This selected design does not estimate the probability that a "
                "neutral or unconditional selected trajectory attains the AF gate."
            ),
        ],
    }
    return manifest, plan, excluded_registry, excluded_inventory


def _validate_seed_plan(
    plan: pd.DataFrame, *, units: pd.DataFrame, excluded: set[int]
) -> None:
    if list(plan.columns) != list(PLAN_COLUMNS) or len(plan) != 60 * ATTEMPTS_PER_UNIT:
        raise ValueError("Han v3 seed plan shape/columns differ")
    plan = plan.copy()
    for column, minimum in {
        "replicate_index": 1,
        "production_unit_seed": 1,
        "attempt_zero_based": 0,
        "simulation_seed": 1,
        "simulation_seed_nonce": 0,
        "panel_seed": 1,
        "panel_seed_nonce": 0,
        "sequence_length_bp": 1,
        "focal_position_bp": 0,
        "candidate_pool_diploids": 1,
        "sample_diploids": 1,
    }.items():
        plan[column] = _exact_integer_series(
            plan, column, label="seed plan", minimum=minimum
        )
    if set(plan["unit_id"].astype(str)) != set(units["unit_id"].astype(str)):
        raise ValueError("Han v3 seed plan unit IDs differ from production")
    counts = plan.groupby("unit_id").size()
    if set(counts.astype(int)) != {ATTEMPTS_PER_UNIT}:
        raise ValueError("Han v3 seed plan does not allocate 300 attempts per unit")
    endpoints = _endpoint_records()
    unit_records = {str(row["unit_id"]): row for row in units.to_dict(orient="records")}
    for unit_id, group in plan.groupby("unit_id", sort=False):
        unit = unit_records[str(unit_id)]
        coefficient = float(unit["selection_coefficient"])
        target = float(unit["target_allele_frequency"])
        endpoint = endpoints[(coefficient, target)]
        expected_stage = (
            STAGE_PILOT
            if int(unit["replicate_index"]) == PILOT_REPLICATE_INDEX
            else STAGE_SIMULATE
        )
        if list(group["attempt_zero_based"].astype(int)) != list(
            range(ATTEMPTS_PER_UNIT)
        ):
            raise ValueError("Han v3 per-unit attempts are not contiguous from zero")
        exact_values: dict[str, int | float | str] = {
            "unit_id": str(unit_id),
            "stage": expected_stage,
            "selection_coefficient": coefficient,
            "target_allele_frequency": target,
            "replicate_index": int(unit["replicate_index"]),
            "production_unit_seed": int(unit["seed"]),
            "realized_selection_end_generations_ago": endpoint[
                "realized_selection_end_generations_ago"
            ],
            "duration_multiplier": endpoint["duration_multiplier"],
        }
        for column, expected in exact_values.items():
            if isinstance(expected, float):
                valid = np.allclose(
                    group[column].astype(float), expected, rtol=0.0, atol=1e-12
                )
            elif isinstance(expected, int):
                valid = set(group[column].astype(int)) == {expected}
            else:
                valid = set(group[column].astype(str)) == {expected}
            if not valid:
                raise ValueError(f"Han v3 schedule {column} differs for unit {unit_id}")
    all_seeds = pd.concat(
        [plan["simulation_seed"], plan["panel_seed"]], ignore_index=True
    ).astype(int)
    if all_seeds.duplicated().any() or set(all_seeds).intersection(excluded):
        raise ValueError("Han v3 attempt/panel seeds are not globally disjoint")
    if set(plan["stage"].astype(str).value_counts().to_dict().items()) != {
        (STAGE_PILOT, 6 * ATTEMPTS_PER_UNIT),
        (STAGE_SIMULATE, 54 * ATTEMPTS_PER_UNIT),
    }:
        raise ValueError("Han v3 pilot/simulate allocation differs")
    fixed = {
        "sequence_length_bp": simulation.SEQUENCE_LENGTH_BP,
        "focal_position_bp": simulation.FOCAL_POSITION_BP,
        "candidate_pool_diploids": CANDIDATE_POOL_DIPLOIDS,
        "sample_diploids": SAMPLE_DIPLOIDS,
        "draw_timeout_seconds": DRAW_TIMEOUT_SECONDS,
        "cumulative_timeout_seconds": CUMULATIVE_TIMEOUT_SECONDS,
    }
    for column, expected in fixed.items():
        if set(plan[column].astype(float)) != {float(expected)}:
            raise ValueError(f"Han v3 seed plan {column} differs")


def publish_plan_bundle(
    output_dir: str | Path,
    *,
    repo_root: str | Path,
    manifest: Mapping[str, Any],
    seed_plan: pd.DataFrame,
    excluded_seed_registry: pd.DataFrame,
    excluded_seed_inventory: Mapping[str, Any],
) -> dict[str, Any]:
    """Publish plan, frozen seed evidence, and registry in one directory rename."""

    root = Path(repo_root).resolve()
    output = _inside_repo(output_dir, root, label="plan output")
    if output.exists():
        bundle = load_bundle(output / DEFAULT_PLAN_MANIFEST_NAME, repo_root=root)
        return dict(bundle.manifest)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.tmp.{os.getpid()}")
    if temporary.exists():
        raise ValueError(f"Han v3 temporary plan directory already exists: {temporary}")
    temporary.mkdir()
    try:
        seed_path = temporary / DEFAULT_SEED_PLAN_NAME
        seed_plan.to_csv(seed_path, sep="\t", index=False, lineterminator="\n")
        registry_path = temporary / DEFAULT_EXCLUDED_SEED_REGISTRY_NAME
        excluded_seed_registry.to_csv(
            registry_path, sep="\t", index=False, lineterminator="\n"
        )
        inventory_path = temporary / DEFAULT_EXCLUDED_SEED_INVENTORY_NAME
        inventory_path.write_text(
            json.dumps(
                dict(excluded_seed_inventory),
                indent=2,
                sort_keys=True,
                allow_nan=False,
            )
            + "\n",
            encoding="utf-8",
        )
        payload = dict(manifest)
        payload["seed_plan_artifact"] = {
            "path": f"{output.relative_to(root).as_posix()}/{DEFAULT_SEED_PLAN_NAME}",
            "sha256": sha256_file(seed_path),
            "rows": len(seed_plan),
        }
        payload["excluded_seed_registry_artifact"] = {
            "schema": EXCLUDED_SEED_REGISTRY_SCHEMA,
            "path": (
                f"{output.relative_to(root).as_posix()}/"
                f"{DEFAULT_EXCLUDED_SEED_REGISTRY_NAME}"
            ),
            "sha256": sha256_file(registry_path),
            "rows": len(excluded_seed_registry),
        }
        payload["excluded_seed_inventory_artifact"] = {
            "schema": EXCLUDED_SEED_INVENTORY_SCHEMA,
            "path": (
                f"{output.relative_to(root).as_posix()}/"
                f"{DEFAULT_EXCLUDED_SEED_INVENTORY_NAME}"
            ),
            "sha256": sha256_file(inventory_path),
            "canonical_artifacts": len(excluded_seed_inventory["canonical_artifacts"]),
            "local_subset_artifacts": len(
                excluded_seed_inventory["local_subset_evidence"]
            ),
        }
        payload["payload_sha256"] = _canonical_sha256(payload)
        (temporary / DEFAULT_PLAN_MANIFEST_NAME).write_text(
            json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, output)
        return payload
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)


def _read_manifest(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError("Han v3 manifest is unreadable") from error
    unhashed = dict(payload)
    digest = unhashed.pop("payload_sha256", None)
    if (
        payload.get("schema") != MANIFEST_SCHEMA
        or payload.get("status") != "planned"
        or digest != _canonical_sha256(unhashed)
    ):
        raise ValueError("Han v3 manifest is corrupt or incompatible")
    return payload


def load_bundle(
    manifest_path: str | Path,
    *,
    repo_root: str | Path,
    slim_path: str | Path | None = None,
) -> DirectV3Bundle:
    root = Path(repo_root).resolve()
    manifest_file = _inside_repo(manifest_path, root, label="manifest")
    manifest = _read_manifest(manifest_file)
    execution_record = manifest["production_execution_units"]
    execution_path = _inside_repo(
        root / execution_record["path"], root, label="execution units"
    )
    campaign = _inside_repo(execution_path.parent, root, label="campaign directory")
    if sha256_file(execution_path) != execution_record["sha256"]:
        raise ValueError("Han v3 execution-unit binding drifted")
    execution = pd.read_csv(execution_path, sep="\t")
    units = _validate_production_units(execution)
    if (
        _canonical_sha256(sorted(units["unit_id"].astype(str)))
        != execution_record["unit_id_sha256"]
    ):
        raise ValueError("Han v3 immutable production unit IDs drifted")
    sources = manifest.get("implementation_sources", {})
    if sources != _source_records(root):
        raise ValueError("Han v3 implementation source binding drifted")
    runtime = manifest.get("runtime")
    if not isinstance(runtime, Mapping) or not isinstance(
        runtime.get("portable_contract"), Mapping
    ):
        raise TypeError("Han v3 portable runtime contract is absent")
    portable_runtime = runtime["portable_contract"]
    recorded_slim = _inside_repo(
        root / portable_runtime["slim"]["relative_path"],
        root,
        label="recorded SLiM",
    )
    requested_slim = (
        recorded_slim
        if slim_path is None
        else _inside_repo(slim_path, root, label="requested SLiM")
    )
    current_runtime = _runtime_record(root, requested_slim)
    if current_runtime["portable_contract"] != portable_runtime:
        raise ValueError("Han v3 runtime/SLiM binding drifted")
    pilot = manifest["failed_v2_pilot_evidence"]
    if any(
        pilot.get(key) is not False
        for key in (
            "inferential_use",
            "endpoint_authorization",
            "production_unit_acceptance_use",
            "gamma_smc_statistics_used",
        )
    ):
        raise ValueError("Han v3 failed-v2 evidence semantics drifted")
    artifact = manifest["seed_plan_artifact"]
    seed_path = _inside_repo(root / artifact["path"], root, label="seed plan")
    if sha256_file(seed_path) != artifact["sha256"]:
        raise ValueError("Han v3 seed-plan checksum drifted")
    plan = pd.read_csv(seed_path, sep="\t")
    registry_artifact = manifest["excluded_seed_registry_artifact"]
    registry_path = _inside_repo(
        root / registry_artifact["path"], root, label="frozen excluded-seed registry"
    )
    if sha256_file(registry_path) != registry_artifact["sha256"]:
        raise ValueError("Han v3 frozen excluded-seed registry checksum drifted")
    registry = pd.read_csv(registry_path, sep="\t")
    excluded = _validate_frozen_seed_registry(registry, manifest["seed_registry"])
    if len(registry) != int(registry_artifact["rows"]):
        raise ValueError("Han v3 frozen excluded-seed registry row count drifted")
    inventory_artifact = manifest["excluded_seed_inventory_artifact"]
    inventory_path = _inside_repo(
        root / inventory_artifact["path"],
        root,
        label="frozen excluded-seed inventory",
    )
    if sha256_file(inventory_path) != inventory_artifact["sha256"]:
        raise ValueError("Han v3 frozen excluded-seed inventory checksum drifted")
    try:
        inventory = json.loads(inventory_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(
            "Han v3 frozen excluded-seed inventory is unreadable"
        ) from error
    unhashed_inventory = dict(inventory)
    inventory_digest = unhashed_inventory.pop("payload_sha256", None)
    if (
        inventory.get("schema") != EXCLUDED_SEED_INVENTORY_SCHEMA
        or inventory.get("status") != "frozen_at_atomic_plan_publication"
        or inventory_digest != _canonical_sha256(unhashed_inventory)
        or int(inventory.get("excluded_union_count", -1)) != len(excluded)
        or inventory.get("excluded_union_sha256") != _seed_digest(excluded)
    ):
        raise ValueError("Han v3 frozen excluded-seed inventory semantics drifted")
    _validate_seed_plan(plan, units=units, excluded=excluded)
    if (
        _canonical_sha256(plan.to_dict(orient="records"))
        != manifest["seed_plan_rows_sha256"]
    ):
        raise ValueError("Han v3 seed-plan semantic digest drifted")
    return DirectV3Bundle(
        repo_root=root,
        campaign_dir=campaign,
        output_dir=manifest_file.parent,
        manifest_path=manifest_file,
        seed_plan_path=seed_path,
        excluded_seed_registry_path=registry_path,
        excluded_seed_inventory_path=inventory_path,
        slim_path=requested_slim,
        manifest=manifest,
        seed_plan=plan,
        excluded_seed_registry=registry,
        excluded_seed_inventory=inventory,
        units=units,
    )


def _unit_seed_schedule(bundle: DirectV3Bundle, unit_id: str) -> pd.DataFrame:
    schedule = bundle.seed_plan[
        bundle.seed_plan["unit_id"].astype(str).eq(str(unit_id))
    ].copy()
    if len(schedule) != ATTEMPTS_PER_UNIT:
        raise ValueError(f"Han v3 seed plan lacks 300 rows for {unit_id}")
    return schedule.sort_values("attempt_zero_based").reset_index(drop=True)


def _v3_implementation_paths(bundle: DirectV3Bundle) -> tuple[str, ...]:
    return tuple(
        dict.fromkeys(
            (
                *_BASE_SIMULATION_IMPLEMENTATION_PATHS,
                MODULE_SOURCE_PATH,
                WRAPPER_SOURCE_PATH,
                "tests/test_focused_han_direct_v3.py",
                "python/gamma_smc_aou/focused_han_fixed_v2_calibration.py",
                "python/gamma_smc_aou/eas_sweep_study.py",
                bundle.manifest_path.relative_to(bundle.repo_root).as_posix(),
                bundle.seed_plan_path.relative_to(bundle.repo_root).as_posix(),
                bundle.excluded_seed_registry_path.relative_to(
                    bundle.repo_root
                ).as_posix(),
                bundle.excluded_seed_inventory_path.relative_to(
                    bundle.repo_root
                ).as_posix(),
            )
        )
    )


def _unit_binding(bundle: DirectV3Bundle, unit: Mapping[str, Any]) -> dict[str, Any]:
    schedule = _unit_seed_schedule(bundle, str(unit["unit_id"]))
    return {
        "schema": AUTHORIZATION_SCHEMA,
        "status": "planned_direct_production_unit",
        "manifest_path": bundle.manifest_path.relative_to(bundle.repo_root).as_posix(),
        "manifest_sha256": sha256_file(bundle.manifest_path),
        "seed_plan_path": bundle.seed_plan_path.relative_to(
            bundle.repo_root
        ).as_posix(),
        "seed_plan_sha256": sha256_file(bundle.seed_plan_path),
        "unit_id": str(unit["unit_id"]),
        "attempt_schedule_sha256": _canonical_sha256(
            schedule.to_dict(orient="records")
        ),
        "realized_selection_end_generations_ago": float(
            schedule.iloc[0]["realized_selection_end_generations_ago"]
        ),
        "gamma_smc_statistics_used": False,
        "failed_v2_inferential_use": False,
    }


def _validate_recorded_event_audit(
    bundle: DirectV3Bundle,
    unit: Mapping[str, Any],
    provenance: Mapping[str, Any],
) -> None:
    schedule = _unit_seed_schedule(bundle, str(unit["unit_id"]))
    endpoint = float(schedule.iloc[0]["realized_selection_end_generations_ago"])
    coefficient = float(unit["selection_coefficient"])
    endpoint_manifest = next(
        item
        for item in bundle.manifest["model"]["fixed_endpoints"]
        if math.isclose(
            float(item["selection_coefficient"]), coefficient, abs_tol=1e-12
        )
        and math.isclose(
            float(item["target_allele_frequency"]),
            float(unit["target_allele_frequency"]),
            abs_tol=1e-12,
        )
    )
    audit = provenance.get("event_audit")
    if not isinstance(audit, Mapping):
        raise TypeError("Han v3 cached completion event audit is absent")
    realized = audit.get("realized_q5_event_times_generations")
    serialized = audit.get("serialized_extended_events")
    boundary = (
        math.floor(fixed_v2.HAN_SPLIT_GENERATIONS / SLIM_SCALING_FACTOR)
        * SLIM_SCALING_FACTOR
    )
    selection_start = float(fixed_v2.INTROGRESSION_PULSE_GENERATIONS)
    expected_intervals = (
        [
            {
                "population": "Loschbour",
                "start_time": selection_start,
                "end_time": endpoint,
            }
        ]
        if endpoint >= boundary
        else [
            {
                "population": "Loschbour",
                "start_time": selection_start,
                "end_time": boundary,
            },
            {
                "population": fixed_v2.INTROGRESSION_TARGET_POPULATION,
                "start_time": boundary,
                "end_time": endpoint,
            },
        ]
    )
    expected_realized = {
        "mutation_age_generations": 2400.0,
        "source_check_generations": 2275.0,
        "introgression_and_selection_onset_generations": 2270.0,
        "demographic_han_split_generations": 2015.0,
        "selection_end_generations": endpoint,
    }

    def event_time(event: Mapping[str, Any], key: str) -> float:
        value = event.get(key)
        if not isinstance(value, Mapping):
            return math.nan
        return float(value.get("generations_ago", math.nan))

    serialized_valid = isinstance(serialized, list) and all(
        isinstance(event, Mapping) for event in serialized
    )
    serialized_events = list(serialized) if serialized_valid else []
    draws = [
        event
        for event in serialized_events
        if event.get("event_type") == "DrawMutation"
    ]
    fitness = [
        event
        for event in serialized_events
        if event.get("event_type") == "ChangeMutationFitness"
    ]
    conditions = [
        event
        for event in serialized_events
        if event.get("event_type") == "ConditionOnAlleleFrequency"
    ]
    observed_intervals = sorted(
        (
            str(event.get("population")),
            event_time(event, "start_time"),
            event_time(event, "end_time"),
        )
        for event in fitness
    )
    expected_interval_tuples = sorted(
        (item["population"], item["start_time"], item["end_time"])
        for item in expected_intervals
    )
    recipient_nonloss = [
        event
        for event in conditions
        if event.get("operator") == ">"
        and math.isclose(
            float(event.get("allele_frequency", math.nan)), 0.0, abs_tol=1e-12
        )
        and event.get("population")
        in {"Loschbour", fixed_v2.INTROGRESSION_TARGET_POPULATION}
    ]
    source_conditions = [
        event for event in conditions if event.get("population") == "Neanderthal"
    ]
    terminal_han = [
        event
        for event in conditions
        if event.get("population") == fixed_v2.INTROGRESSION_TARGET_POPULATION
        and math.isclose(event_time(event, "start_time"), 0.0, abs_tol=1e-12)
        and math.isclose(event_time(event, "end_time"), 0.0, abs_tol=1e-12)
    ]
    serialized_semantics_valid = (
        serialized_valid
        and len(draws) == 1
        and draws[0].get("population") == "Neanderthal"
        and math.isclose(event_time(draws[0], "time"), 2400.0, abs_tol=1e-12)
        and observed_intervals == expected_interval_tuples
        and all(
            math.isclose(
                float(event.get("selection_coefficient", math.nan)),
                coefficient,
                abs_tol=1e-12,
            )
            and math.isclose(
                float(event.get("dominance_coefficient", math.nan)),
                0.5,
                abs_tol=1e-12,
            )
            for event in fitness
        )
        and len(recipient_nonloss) == 2
        and len(source_conditions) == 2
        and not terminal_han
        and all(
            event.get("single_site_id") == fixed_v2.FOCAL_SITE_ID
            for event in serialized_events
        )
    )
    if (
        provenance.get("conditional_inference") is not True
        or not math.isclose(
            float(provenance.get("fixed_cessation_endpoint", math.nan)),
            endpoint,
            abs_tol=1e-12,
        )
        or audit.get("schema") != fixed_v2.EVENT_AUDIT_SCHEMA
        or audit.get("mutation_population") != "Neanderthal"
        or not math.isclose(
            float(audit.get("mutation_age_generations", math.nan)),
            2400.0,
            abs_tol=1e-12,
        )
        or not math.isclose(
            float(audit.get("selection_coefficient", math.nan)),
            coefficient,
            abs_tol=1e-12,
        )
        or not math.isclose(
            float(audit.get("duration_multiplier", math.nan)),
            float(endpoint_manifest["duration_multiplier"]),
            abs_tol=1e-12,
        )
        or not math.isclose(
            float(audit.get("requested_selection_end_generations_ago", math.nan)),
            float(endpoint_manifest["requested_selection_end_generations_ago"]),
            abs_tol=1e-12,
        )
        or not math.isclose(
            float(audit.get("realized_selection_end_generations_ago", math.nan)),
            endpoint,
            abs_tol=1e-12,
        )
        or audit.get("selection_intervals") != expected_intervals
        or int(audit.get("selection_interval_count", -1)) != len(expected_intervals)
        or int(audit.get("recipient_nonloss_condition_count", -1)) != 2
        or int(audit.get("source_condition_count", -1)) != 2
        or int(audit.get("terminal_han_condition_count", -1)) != 0
        or audit.get("external_terminal_ascertainment") is not True
        or realized != expected_realized
        or audit.get("realized_event_schedule_sha256")
        != _canonical_sha256(expected_realized)
        or not serialized_semantics_valid
        or audit.get("serialized_extended_events_sha256")
        != _canonical_sha256(serialized)
        or int(audit.get("focal_event_count", -1)) != len(serialized)
        or not math.isclose(
            endpoint / SLIM_SCALING_FACTOR, round(endpoint / 5), abs_tol=1e-9
        )
    ):
        raise ValueError("Han v3 cached completion event audit differs")


_DRAW_TIMEOUT_LOCAL = threading.local()
_PR_SET_CHILD_SUBREAPER = 36


@contextmanager
def _supervised_draw_timeout(seconds: float) -> Iterator[None]:
    timeout = float(seconds)
    if not math.isfinite(timeout) or timeout <= 0:
        raise ValueError("Han v3 supervised timeout must be finite and positive")
    previous = getattr(_DRAW_TIMEOUT_LOCAL, "seconds", None)
    _DRAW_TIMEOUT_LOCAL.seconds = timeout
    try:
        yield
    finally:
        if previous is None:
            try:
                del _DRAW_TIMEOUT_LOCAL.seconds
            except AttributeError:
                pass
        else:
            _DRAW_TIMEOUT_LOCAL.seconds = previous


def _enable_linux_child_subreaper() -> None:
    if os.name != "posix" or platform.system() != "Linux":
        raise RuntimeError("Han v3 process supervision requires Linux/WSL")
    libc = ctypes.CDLL(None, use_errno=True)
    result = libc.prctl(_PR_SET_CHILD_SUBREAPER, 1, 0, 0, 0)
    if result != 0:
        error_number = ctypes.get_errno()
        raise OSError(error_number, "PR_SET_CHILD_SUBREAPER failed")


def _process_group_exists(pgid: int | None) -> bool:
    if pgid is None:
        return False
    try:
        os.killpg(int(pgid), 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def _linux_process_identity(pid: int) -> dict[str, Any] | None:
    try:
        boot_id = (
            Path("/proc/sys/kernel/random/boot_id").read_text(encoding="utf-8").strip()
        )
        stat = Path(f"/proc/{int(pid)}/stat").read_text(encoding="utf-8")
        fields = stat[stat.rfind(")") + 2 :].split()
        start_ticks = int(fields[19])
    except (FileNotFoundError, IndexError, OSError, ValueError):
        return None
    return {"pid": int(pid), "boot_id": boot_id, "start_ticks": start_ticks}


def _linux_child_processes(parent_pid: int) -> list[int]:
    try:
        text = Path(
            f"/proc/{int(parent_pid)}/task/{int(parent_pid)}/children"
        ).read_text(encoding="utf-8")
    except OSError:
        return []
    return [int(value) for value in text.split()]


def _linux_descendant_processes(parent_pid: int) -> list[int]:
    """Snapshot every descendant, including children adopted by a subreaper."""

    descendants: list[int] = []
    pending = list(_linux_child_processes(parent_pid))
    seen: set[int] = set()
    while pending:
        child_pid = pending.pop()
        if child_pid in seen:
            continue
        seen.add(child_pid)
        descendants.append(child_pid)
        pending.extend(_linux_child_processes(child_pid))
    return descendants


def _reap_available_children() -> None:
    while True:
        try:
            child_pid, _status = os.waitpid(-1, os.WNOHANG)
        except ChildProcessError:
            return
        if child_pid <= 0:
            return


def _terminate_descendants_and_reap(
    parent_pid: int, *, term_grace_seconds: float
) -> None:
    """Watchdog stays alive until its complete descendant tree is dead and reaped."""

    grace = max(0.0, float(term_grace_seconds))
    known: dict[int, dict[str, Any]] = {}

    def refresh_known() -> None:
        for child_pid in _linux_descendant_processes(parent_pid):
            identity = _linux_process_identity(child_pid)
            if identity is not None:
                known[child_pid] = identity

    def live_known() -> list[int]:
        return [pid for pid, identity in known.items() if _same_linux_process(identity)]

    refresh_known()
    for child_pid in live_known():
        try:
            os.kill(child_pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
    deadline = monotonic() + grace
    while monotonic() < deadline:
        _reap_available_children()
        refresh_known()
        if not live_known() and not _linux_descendant_processes(parent_pid):
            return
        sleep(SUPERVISOR_POLL_SECONDS)

    # Re-enumerate on every pass, but retain prior identities as well. After an
    # intermediate parent dies, reparenting can transiently hide a descendant
    # from /proc/.../children even though its original PID is still live.
    deadline = monotonic() + max(1.0, grace)
    while monotonic() < deadline:
        refresh_known()
        for child_pid in live_known():
            try:
                os.kill(child_pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        _reap_available_children()
        refresh_known()
        if not live_known() and not _linux_descendant_processes(parent_pid):
            return
        sleep(SUPERVISOR_POLL_SECONDS)
    raise RuntimeError("Han v3 watchdog could not reap every descendant")


def _same_linux_process(identity: Mapping[str, Any] | None) -> bool:
    if not isinstance(identity, Mapping) or "pid" not in identity:
        return False
    return _linux_process_identity(int(identity["pid"])) == dict(identity)


def _parse_linux_process_identity(value: Any) -> Mapping[str, Any] | None:
    if isinstance(value, Mapping):
        return value
    if not isinstance(value, str) or not value.strip():
        return None
    try:
        parsed = ast.literal_eval(value)
    except (SyntaxError, ValueError):
        return None
    return parsed if isinstance(parsed, Mapping) else None


def _validate_linux_process_identity(
    value: Any, *, expected_pid: int
) -> Mapping[str, Any]:
    identity = _parse_linux_process_identity(value)
    if (
        identity is None
        or set(identity) != {"pid", "boot_id", "start_ticks"}
        or isinstance(identity["pid"], bool)
        or isinstance(identity["start_ticks"], bool)
        or not isinstance(identity["pid"], int)
        or not isinstance(identity["start_ticks"], int)
        or not isinstance(identity["boot_id"], str)
        or not identity["boot_id"].strip()
        or int(identity["pid"]) != int(expected_pid)
        or int(identity["pid"]) < 1
        or int(identity["start_ticks"]) < 0
    ):
        raise ValueError("Han v3 Linux process identity is malformed")
    return identity


def _process_state(pid: int) -> str | None:
    try:
        stat = Path(f"/proc/{int(pid)}/stat").read_text(encoding="utf-8")
    except OSError:
        return None
    fields = stat[stat.rfind(")") + 2 :].split()
    return fields[0] if fields else None


def _process_absent_or_zombie(pid: int) -> bool:
    return _process_state(int(pid)) in {None, "Z", "X"}


def _process_group_has_live_members(pgid: int | None) -> bool:
    if pgid is None:
        return False
    for stat_path in Path("/proc").glob("[0-9]*/stat"):
        try:
            stat = stat_path.read_text(encoding="utf-8")
            fields = stat[stat.rfind(")") + 2 :].split()
            state = fields[0]
            process_group = int(fields[2])
        except (IndexError, OSError, ValueError):
            continue
        if process_group == int(pgid) and state not in {"Z", "X"}:
            return True
    return False


def _reap_if_our_child(pid: int) -> bool:
    try:
        waited_pid, _status = os.waitpid(int(pid), os.WNOHANG)
    except ChildProcessError:
        return False
    return waited_pid == int(pid)


def _write_active_supervision(path: Path, record: Mapping[str, Any]) -> None:
    payload = dict(record)
    payload["payload_sha256"] = _canonical_sha256(payload)
    _atomic_json(path, payload)


def _read_active_supervision(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError("Han v3 active-supervision record is unreadable") from error
    unhashed = dict(payload)
    digest = unhashed.pop("payload_sha256", None)
    if (
        payload.get("schema") != SUPERVISOR_ACTIVE_SCHEMA
        or payload.get("status") != "active"
        or digest != _canonical_sha256(unhashed)
    ):
        raise ValueError("Han v3 active-supervision record is corrupt")
    return payload


def _reap_adopted_group_children(pgid: int, *, deadline: float) -> int:
    reaped = 0
    while monotonic() < deadline:
        try:
            child_pid, _status = os.waitpid(-int(pgid), os.WNOHANG)
        except ChildProcessError:
            break
        if child_pid > 0:
            reaped += 1
            continue
        if not _process_group_exists(pgid):
            break
        sleep(SUPERVISOR_POLL_SECONDS)
    return reaped


def _slim_engine_child(
    sender: Any,
    engine: Any,
    args: tuple[Any, ...],
    kwargs: dict[str, Any],
    tree_path: str,
    watchdog_read_fd: int,
    watchdog_write_fd: int,
    expected_parent_pid: int,
    term_grace_seconds: float,
) -> None:
    """Supervise an engine wrapper inside one fresh session/process group."""

    try:
        os.close(watchdog_write_fd)
        os.setsid()
        _enable_linux_child_subreaper()
        if os.getppid() != int(expected_parent_pid):
            _terminate_descendants_and_reap(
                os.getpid(), term_grace_seconds=term_grace_seconds
            )
            os._exit(70)
        attempt_temp_dir = str(Path(tree_path).parent)
        os.environ["TMPDIR"] = attempt_temp_dir
        os.environ["TEMP"] = attempt_temp_dir
        os.environ["TMP"] = attempt_temp_dir
        tempfile.tempdir = attempt_temp_dir
        session_identity_path = Path(tree_path).with_name(
            "session_process_identity.json"
        )
        _atomic_json(
            session_identity_path,
            {
                "session_pid": os.getpid(),
                "session_pgid": os.getpgid(0),
                "session_sid": os.getsid(0),
                "session_process_identity": _linux_process_identity(os.getpid()),
                "parent_death_watchdog_armed": True,
            },
        )
        sender.send(
            {
                "event": "started",
                "pid": os.getpid(),
                "pgid": os.getpgid(0),
                "sid": os.getsid(0),
                "parent_death_watchdog_armed": True,
                "parent_death_watchdog_kind": "session_leader_control_pipe_loop",
            }
        )
        status_path = Path(tree_path).with_name("engine_status.json")
        engine_pid = os.fork()
        if engine_pid == 0:
            try:
                tree_sequence = engine.simulate(*args, **kwargs)
                tree_sequence.dump(tree_path)
                _atomic_json(status_path, {"event": "complete"})
                os._exit(0)
            except BaseException as error:  # noqa: BLE001
                _atomic_json(
                    status_path,
                    {
                        "event": "error",
                        "error_type": type(error).__name__,
                        "error": str(error),
                        "traceback": traceback.format_exc(),
                    },
                )
                os._exit(1)
        engine_identity = _linux_process_identity(engine_pid)
        _atomic_json(
            Path(tree_path).with_name("engine_process_identity.json"),
            {
                "engine_pid": engine_pid,
                "engine_process_identity": engine_identity,
            },
        )
        sender.send(
            {
                "event": "engine_started",
                "engine_pid": engine_pid,
                "engine_process_identity": engine_identity,
            }
        )
        wait_status: int | None = None
        while wait_status is None:
            waited_pid, candidate_status = os.waitpid(engine_pid, os.WNOHANG)
            if waited_pid == engine_pid:
                wait_status = candidate_status
                break
            try:
                readable, _, _ = select.select(
                    [watchdog_read_fd], [], [], SUPERVISOR_POLL_SECONDS
                )
                parent_alive = not readable or os.read(watchdog_read_fd, 1) != b""
            except OSError:
                parent_alive = False
            if not parent_alive or os.getppid() != int(expected_parent_pid):
                _terminate_descendants_and_reap(
                    os.getpid(), term_grace_seconds=term_grace_seconds
                )
                os._exit(70)
        if wait_status is None:
            raise RuntimeError("Han v3 engine wrapper was not reaped")
        try:
            status = json.loads(status_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            raise RuntimeError(
                f"Han v3 engine wrapper exited with wait status {wait_status}"
            ) from error
        sender.send(status)
    except BaseException as error:  # noqa: BLE001 - marshal the child failure
        try:
            sender.send(
                {
                    "event": "error",
                    "error_type": type(error).__name__,
                    "error": str(error),
                    "traceback": traceback.format_exc(),
                }
            )
        except (BrokenPipeError, EOFError, OSError):
            pass
    finally:
        try:
            os.close(watchdog_read_fd)
        except OSError:
            pass
        sender.close()


def _stop_supervised_process(
    process: multiprocessing.Process,
    *,
    pgid: int | None,
    term_grace_seconds: float,
    record: dict[str, Any],
) -> None:
    grace = float(term_grace_seconds)
    if pgid is None and process.pid is not None:
        # setsid() is the child contract; this remains the only exact candidate
        # group even if the startup handshake loses a scheduling race.
        pgid = int(process.pid)
        record["child_pgid"] = pgid
        record["child_sid"] = pgid
    group_established = False
    if pgid is not None and process.pid is not None:
        try:
            group_established = os.getpgid(int(process.pid)) == int(pgid)
        except ProcessLookupError:
            group_established = _process_group_exists(pgid)
    if group_established:
        try:
            os.killpg(int(pgid), signal.SIGTERM)
            record["sigterm_sent"] = True
            record["parent_death_watchdog_armed"] = True
        except ProcessLookupError:
            pass
    if not group_established and process.is_alive():
        process.terminate()
        record["sigterm_sent"] = True
    grace_deadline = monotonic() + grace
    while monotonic() < grace_deadline:
        process.join(timeout=SUPERVISOR_POLL_SECONDS)
        if not process.is_alive() and not _process_group_exists(pgid):
            break
        sleep(SUPERVISOR_POLL_SECONDS)
    if _process_group_exists(pgid):
        try:
            os.killpg(int(pgid), signal.SIGKILL)
            record["sigkill_sent"] = True
        except ProcessLookupError:
            pass
    elif process.is_alive():
        process.kill()
        record["sigkill_sent"] = True
    process.join(timeout=max(1.0, grace))
    if process.is_alive():
        raise RuntimeError("Han v3 supervisor could not reap its direct child")
    record["direct_child_reaped"] = process.exitcode is not None
    if pgid is not None:
        record["adopted_children_reaped"] = _reap_adopted_group_children(
            pgid, deadline=monotonic() + max(1.0, grace)
        )
    record["process_group_gone"] = pgid is not None and not _process_group_exists(pgid)
    if not record["process_group_gone"]:
        raise RuntimeError("Han v3 supervisor left a live process-group descendant")


def _supervise_engine_call(
    engine: Any,
    args: tuple[Any, ...],
    kwargs: dict[str, Any],
    *,
    tree_path: Path,
    timeout_seconds: float,
    term_grace_seconds: float,
    record: dict[str, Any],
    progress_callback: Callable[[], None],
) -> tskit.TreeSequence:
    """Return one tree sequence or reap the entire attempt process group."""

    _enable_linux_child_subreaper()
    context = multiprocessing.get_context("fork")
    receiver, sender = context.Pipe(duplex=False)
    watchdog_read_fd, watchdog_write_fd = os.pipe()
    expected_parent_pid = os.getpid()
    process = context.Process(
        target=_slim_engine_child,
        args=(
            sender,
            engine,
            args,
            kwargs,
            str(tree_path),
            watchdog_read_fd,
            watchdog_write_fd,
            expected_parent_pid,
            term_grace_seconds,
        ),
        name="han-v3-slim-draw",
    )
    started = monotonic()
    try:
        process.start()
    except BaseException:
        os.close(watchdog_read_fd)
        os.close(watchdog_write_fd)
        receiver.close()
        sender.close()
        raise
    os.close(watchdog_read_fd)
    sender.close()
    record.update(
        {
            "supervisor_pid": os.getpid(),
            "supervisor_process_identity": _linux_process_identity(os.getpid()),
            "child_pid": int(process.pid),
            "child_process_identity": _linux_process_identity(int(process.pid)),
            "child_identity_observed": True,
            "child_pgid": int(process.pid),
            "child_sid": int(process.pid),
            "parent_death_watchdog_armed": False,
            "engine_wrapper_pid": None,
            "engine_wrapper_process_identity": None,
            "sigterm_sent": False,
            "sigkill_sent": False,
            "direct_child_reaped": False,
            "adopted_children_reaped": 0,
            "process_group_gone": False,
        }
    )
    progress_callback()
    deadline = started + float(timeout_seconds)
    pgid: int | None = None
    terminal_message: dict[str, Any] | None = None
    try:
        while monotonic() < deadline:
            wait_seconds = min(SUPERVISOR_POLL_SECONDS, deadline - monotonic())
            if receiver.poll(max(0.0, wait_seconds)):
                message = receiver.recv()
                if message.get("event") == "started":
                    child_pid = int(message["pid"])
                    candidate_pgid = int(message["pgid"])
                    candidate_sid = int(message["sid"])
                    if (
                        child_pid != int(process.pid)
                        or candidate_pgid != child_pid
                        or candidate_sid != child_pid
                        or message.get("parent_death_watchdog_armed") is not True
                    ):
                        terminal_message = {
                            "event": "error",
                            "error_type": "SupervisorProtocolError",
                            "error": "child did not enter its own session/process group",
                        }
                        break
                    pgid = candidate_pgid
                    record["child_pgid"] = pgid
                    record["child_sid"] = candidate_sid
                    record["parent_death_watchdog_armed"] = True
                    record["child_process_identity"] = _linux_process_identity(
                        child_pid
                    )
                    progress_callback()
                    continue
                if message.get("event") == "engine_started":
                    record["engine_wrapper_pid"] = int(message["engine_pid"])
                    record["engine_wrapper_process_identity"] = message.get(
                        "engine_process_identity"
                    )
                    progress_callback()
                    continue
                terminal_message = message
                break
            if not process.is_alive() and not receiver.poll():
                terminal_message = {
                    "event": "error",
                    "error_type": "ChildExitError",
                    "error": "supervised child exited without a terminal message",
                }
                break
        if terminal_message is None:
            record["supervisor_status"] = "timed_out"
            _stop_supervised_process(
                process,
                pgid=pgid,
                term_grace_seconds=term_grace_seconds,
                record=record,
            )
            raise simulation.SimulationDrawTimeout(
                f"SLiM draw exceeded {float(timeout_seconds):.3f} wall-clock seconds"
            )
        remaining = max(0.0, deadline - monotonic())
        process.join(timeout=remaining)
        if process.is_alive():
            record["supervisor_status"] = "timed_out_after_terminal_message"
            _stop_supervised_process(
                process,
                pgid=pgid,
                term_grace_seconds=term_grace_seconds,
                record=record,
            )
            raise simulation.SimulationDrawTimeout(
                f"SLiM draw exceeded {float(timeout_seconds):.3f} wall-clock seconds"
            )
        record["direct_child_reaped"] = process.exitcode is not None
        if _process_group_exists(pgid):
            record["supervisor_status"] = "descendant_leak"
            _stop_supervised_process(
                process,
                pgid=pgid,
                term_grace_seconds=term_grace_seconds,
                record=record,
            )
            raise RuntimeError("Han v3 supervised draw left a descendant process")
        record["process_group_gone"] = True
        if terminal_message.get("event") == "error":
            record["supervisor_status"] = "child_error"
            record["child_error_type"] = str(terminal_message.get("error_type", ""))
            record["child_error"] = str(terminal_message.get("error", ""))
            raise RuntimeError(
                "Han v3 supervised SLiM child failed: "
                f"{record['child_error_type']}: {record['child_error']}"
            )
        if terminal_message.get("event") != "complete" or not tree_path.is_file():
            record["supervisor_status"] = "protocol_error"
            raise RuntimeError("Han v3 supervised SLiM result tree is absent")
        record["temp_tree_created"] = True
        tree_sequence = tskit.load(tree_path)
        record["temp_tree_loaded"] = True
        record["supervisor_status"] = "complete"
        return tree_sequence
    finally:
        os.close(watchdog_write_fd)
        receiver.close()
        record["supervisor_elapsed_seconds"] = monotonic() - started


class _SupervisedSlimEngine:
    """Worker-local engine proxy enforcing one exact process group per draw."""

    def __init__(
        self,
        engine: Any,
        *,
        unit_dir: Path,
        schedule: pd.DataFrame,
        term_grace_seconds: float = SUPERVISOR_TERM_GRACE_SECONDS,
    ) -> None:
        self._engine = engine
        self._unit_dir = unit_dir
        self._schedule = schedule.set_index("attempt_zero_based")
        self._term_grace_seconds = float(term_grace_seconds)
        self._ledger_path = unit_dir / "supervision.tsv"
        self._reconciliation_path = unit_dir / "supervision_reconciliation.tsv"
        self._startup_cleanup_path = unit_dir / "supervision_startup_cleanup.tsv"
        self.records: dict[int, dict[str, Any]] = {}
        self._reconcile_pre_marker_orphans()
        if self._ledger_path.is_file():
            prior = pd.read_csv(self._ledger_path, sep="\t")
            if prior["attempt_zero_based"].astype(int).duplicated().any():
                raise ValueError("Han v3 supervisor ledger contains duplicate attempts")
            self.records = {
                int(row["attempt_zero_based"]): row
                for row in prior.to_dict(orient="records")
            }
        self._reconcile_interrupted_attempts()
        self._reconcile_uncommitted_supervision()

    @staticmethod
    def _startup_orphan_descriptor(path: Path) -> dict[str, Any]:
        if path.is_symlink():
            raise ValueError("Han v3 startup orphan must not be a symlink")
        if path.is_file():
            artifact_kind = (
                "active_atomic_temp"
                if path.name.startswith(".han_v3_active_")
                else "known_atomic_temp"
            )
            inventory = [
                {
                    "path": path.name,
                    "kind": "file",
                    "size": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
            ]
        elif path.is_dir():
            artifact_kind = "attempt_scratch_directory"
            inventory = []
            for candidate in sorted(path.rglob("*"), key=lambda item: item.as_posix()):
                if candidate.is_symlink():
                    raise ValueError("Han v3 startup orphan contains a symlink")
                relative = candidate.relative_to(path).as_posix()
                if candidate.is_dir():
                    inventory.append({"path": relative, "kind": "directory"})
                elif candidate.is_file():
                    inventory.append(
                        {
                            "path": relative,
                            "kind": "file",
                            "size": candidate.stat().st_size,
                            "sha256": sha256_file(candidate),
                        }
                    )
                else:
                    raise ValueError("Han v3 startup orphan has a special file")
        else:
            raise ValueError("Han v3 startup orphan is not a file/directory")
        descriptor = {
            "artifact_name": path.name,
            "artifact_kind": artifact_kind,
            "inventory_sha256": _canonical_sha256(inventory),
        }
        descriptor["cleanup_id"] = _canonical_sha256(descriptor)
        return descriptor

    def _reconcile_pre_marker_orphans(self) -> None:
        records: list[dict[str, Any]] = []
        if self._startup_cleanup_path.is_file():
            frame = pd.read_csv(self._startup_cleanup_path, sep="\t")
            required = {
                "schema",
                "cleanup_id",
                "artifact_name",
                "artifact_kind",
                "inventory_sha256",
                "status",
                "discovered_unix_ns",
                "removed_unix_ns",
            }
            if required.difference(frame.columns):
                raise ValueError("Han v3 startup-cleanup ledger differs")
            records = frame.to_dict(orient="records")
        referenced_scratch: set[str] = set()
        canonical_active_attempts: set[str] = set()
        for active_path in sorted(self._unit_dir.glob(".han_v3_active_*.json")):
            match = re.fullmatch(r"\.han_v3_active_([0-9]{3})\.json", active_path.name)
            if match is None:
                continue
            canonical_active_attempts.add(match.group(1))
            active = _read_active_supervision(active_path)
            temp_name = str(active.get("attempt_record", {}).get("temp_directory", ""))
            if temp_name:
                referenced_scratch.add(temp_name)

        def eligible(path: Path) -> bool:
            if path.is_dir():
                return (
                    re.fullmatch(
                        r"\.han_v3_attempt_[0-9]{3}_seed_[0-9]+_[A-Za-z0-9_-]+",
                        path.name,
                    )
                    is not None
                    and path.name not in referenced_scratch
                )
            match = re.fullmatch(
                r"\.han_v3_active_([0-9]{3})\.json\.tmp\.[0-9]+", path.name
            )
            if match is not None:
                return match.group(1) not in canonical_active_attempts
            return (
                _is_unit_atomic_temp_name(path.name) and not canonical_active_attempts
            )

        candidates = sorted(
            {
                *self._unit_dir.glob(".han_v3_attempt_*"),
                *self._unit_dir.glob(".han_v3_active_*.json.tmp.*"),
                *_unit_atomic_temp_paths(self._unit_dir),
            },
            key=lambda path: path.name,
        )
        for path in candidates:
            if not eligible(path):
                continue
            descriptor = self._startup_orphan_descriptor(path)
            matches = [
                row
                for row in records
                if str(row.get("cleanup_id", "")) == descriptor["cleanup_id"]
            ]
            if len(matches) > 1:
                raise ValueError("Han v3 startup-cleanup identity duplicates")
            if matches:
                row = matches[0]
                if any(
                    str(row.get(field, "")) != str(descriptor[field])
                    for field in (
                        "artifact_name",
                        "artifact_kind",
                        "inventory_sha256",
                    )
                ):
                    raise ValueError("Han v3 startup-cleanup identity conflicts")
            else:
                row = {
                    "schema": SUPERVISOR_STARTUP_CLEANUP_SCHEMA,
                    **descriptor,
                    "status": "pending",
                    "discovered_unix_ns": time_ns(),
                    "removed_unix_ns": 0,
                }
                records.append(row)
                _atomic_frame(self._startup_cleanup_path, pd.DataFrame(records))
            if path.is_dir():
                shutil.rmtree(path)
            else:
                path.unlink()
            row["status"] = "removed"
            row["removed_unix_ns"] = time_ns()
            _atomic_frame(self._startup_cleanup_path, pd.DataFrame(records))

        changed = False
        by_name = {path.name: path for path in candidates}
        for row in records:
            if str(row.get("status")) != "pending":
                continue
            path = by_name.get(str(row["artifact_name"]), self._unit_dir / "__absent__")
            if path.exists():
                raise ValueError("Han v3 pending startup cleanup remains present")
            row["status"] = "removed"
            row["removed_unix_ns"] = time_ns()
            changed = True
        if changed:
            _atomic_frame(self._startup_cleanup_path, pd.DataFrame(records))

    def _append_reconciliation(self, row: Mapping[str, Any]) -> dict[str, Any]:
        records: list[dict[str, Any]] = []
        if self._reconciliation_path.is_file():
            records = pd.read_csv(self._reconciliation_path, sep="\t").to_dict(
                orient="records"
            )
        execution_id = str(row.get("execution_id", ""))
        matching = [
            record
            for record in records
            if str(record.get("execution_id", "")) == execution_id
        ]
        if execution_id and matching:
            if len(matching) != 1:
                raise ValueError("Han v3 reconciliation execution identity duplicates")
            durable = dict(matching[0])

            numeric_invariants = {
                "attempt_zero_based",
                "simulation_seed",
                "prior_supervisor_pid",
                "prior_child_pid",
                "prior_child_pgid",
            }

            def comparable(field: str, value: Any) -> str:
                if value is None or pd.isna(value):
                    return ""
                if field in numeric_invariants:
                    numeric = float(value)
                    if not math.isfinite(numeric) or numeric != round(numeric):
                        return "invalid-numeric"
                    return str(int(numeric))
                if field == "process_group_gone":
                    return str(value).casefold()
                return str(value)

            invariant_fields = (
                "schema",
                "reason",
                "execution_id",
                "attempt_zero_based",
                "simulation_seed",
                "active_record_path",
                "active_record_sha256",
                "prior_supervisor_pid",
                "prior_child_pid",
                "prior_child_pgid",
                "process_group_gone",
            )
            if any(
                comparable(field, durable.get(field))
                != comparable(field, row.get(field))
                for field in invariant_fields
            ):
                raise ValueError("Han v3 reconciliation execution identity conflicts")
            try:
                prior_removed = json.loads(
                    str(durable.get("atomic_temp_files_removed", "[]"))
                )
                replay_removed = json.loads(
                    str(row.get("atomic_temp_files_removed", "[]"))
                )
            except json.JSONDecodeError as error:
                raise ValueError(
                    "Han v3 reconciliation cleanup journal is corrupt"
                ) from error
            if not isinstance(prior_removed, list) or not isinstance(
                replay_removed, list
            ):
                raise ValueError("Han v3 reconciliation cleanup journal is corrupt")
            combined_removed = sorted(set(prior_removed) | set(replay_removed))
            if combined_removed != prior_removed:
                durable["atomic_temp_files_removed"] = json.dumps(
                    combined_removed, separators=(",", ":")
                )
                durable["atomic_temp_file_count_removed"] = len(combined_removed)
                records = [
                    durable
                    if str(record.get("execution_id", "")) == execution_id
                    else record
                    for record in records
                ]
                _atomic_frame(self._reconciliation_path, pd.DataFrame(records))
            return durable
        durable = dict(row)
        records.append(durable)
        _atomic_frame(self._reconciliation_path, pd.DataFrame(records))
        return durable

    def _reconciliation_for_execution(self, execution_id: str) -> dict[str, Any] | None:
        if not self._reconciliation_path.is_file():
            return None
        rows = pd.read_csv(self._reconciliation_path, sep="\t")
        if "execution_id" not in rows:
            raise ValueError("Han v3 reconciliation lacks execution identity")
        matches = rows[rows["execution_id"].astype(str).eq(str(execution_id))]
        if len(matches) > 1:
            raise ValueError("Han v3 reconciliation execution identity duplicates")
        return None if matches.empty else matches.iloc[0].to_dict()

    @staticmethod
    def _execution_id(record: Mapping[str, Any]) -> str:
        raw_existing = record.get("execution_id")
        existing = (
            ""
            if raw_existing is None or pd.isna(raw_existing)
            else str(raw_existing).strip()
        )
        if existing:
            if re.fullmatch(r"[0-9a-f]{64}", existing) is None:
                raise ValueError("Han v3 execution identity is malformed")
            return existing
        return _canonical_sha256(
            {
                "attempt_zero_based": int(record["attempt_zero_based"]),
                "simulation_seed": int(record["simulation_seed"]),
                "supervisor_pid": int(record.get("supervisor_pid", 0)),
                "supervisor_started_unix_ns": int(
                    record.get("supervisor_started_unix_ns", 0)
                ),
            }
        )

    def _read_base_attempts(self) -> pd.DataFrame:
        path = self._unit_dir / "attempts.tsv"
        if not path.is_file():
            return pd.DataFrame()
        attempts = pd.read_csv(path, sep="\t")
        if "attempt_zero_based" not in attempts:
            raise ValueError("Han v3 base attempt ledger lacks attempt identity")
        observed = list(
            _exact_integer_series(
                attempts,
                "attempt_zero_based",
                label="base attempt ledger",
                minimum=0,
            )
        )
        if observed != list(range(len(observed))):
            raise ValueError("Han v3 base attempt ledger is noncontiguous")
        return attempts

    def _interrupted_attempt_row(
        self, record: Mapping[str, Any], *, elapsed_seconds: float
    ) -> dict[str, Any]:
        attempt = int(record["attempt_zero_based"])
        schedule = self._schedule.loc[attempt]
        return {
            "attempt_zero_based": attempt,
            "seed": int(record["simulation_seed"]),
            "status": "timed_out",
            "timeout_seconds": float(record.get("timeout_seconds", 0.0)),
            "elapsed_seconds": float(elapsed_seconds),
            "error": "worker interrupted before base attempt commit",
            "accepted": False,
            "rejection_reason": "worker_interrupted_before_base_attempt_commit",
            "planned_simulation_seed": int(record["simulation_seed"]),
            "panel_seed": int(schedule["panel_seed"]),
            "panel_seed_used": False,
        }

    def _normalize_interrupted_supervision(
        self,
        record: Mapping[str, Any],
        *,
        elapsed_seconds: float,
        temp_removed: bool,
    ) -> dict[str, Any]:
        normalized = dict(record)
        child_pid = int(normalized.get("child_pid") or 0)
        startup_interrupted = child_pid <= 0
        if startup_interrupted:
            normalized["child_pid"] = None
            normalized["child_pgid"] = None
            normalized["child_sid"] = None
            normalized["child_process_identity"] = None
        persisted_supervisor_elapsed = pd.to_numeric(
            pd.Series([normalized.get("supervisor_elapsed_seconds")]),
            errors="coerce",
        ).iloc[0]
        if pd.isna(persisted_supervisor_elapsed) or persisted_supervisor_elapsed <= 0:
            persisted_supervisor_elapsed = min(
                float(elapsed_seconds), float(normalized["timeout_seconds"])
            )
        normalized.update(
            {
                "schema": SUPERVISOR_LEDGER_SCHEMA,
                "execution_id": self._execution_id(normalized),
                "supervisor_status": (
                    "startup_interrupted_before_child_observed"
                    if startup_interrupted
                    else "interrupted_by_worker_death"
                ),
                "child_pgid": (
                    None
                    if startup_interrupted
                    else int(normalized.get("child_pgid") or child_pid)
                ),
                "child_sid": (
                    None
                    if startup_interrupted
                    else int(normalized.get("child_sid") or child_pid)
                ),
                "child_identity_observed": not startup_interrupted,
                "direct_child_reaped": not startup_interrupted,
                "process_group_gone": True,
                "temp_tree_removed": bool(temp_removed),
                "supervisor_elapsed_seconds": max(
                    float(persisted_supervisor_elapsed), 1e-9
                ),
                "adopted_children_reaped": int(
                    normalized.get("adopted_children_reaped", 0) or 0
                ),
                "sigterm_sent": (
                    False
                    if startup_interrupted
                    else bool(normalized.get("sigterm_sent", False))
                ),
                "sigkill_sent": (
                    False
                    if startup_interrupted
                    else bool(normalized.get("sigkill_sent", False))
                ),
                "parent_death_watchdog_armed": (
                    False
                    if startup_interrupted
                    else bool(normalized.get("parent_death_watchdog_armed", True))
                ),
            }
        )
        return normalized

    def _clean_incomplete_unit_artifacts(self) -> None:
        for relative in (
            "simulation.trees",
            "sample_manifest.tsv",
            "truth_profiles.tsv.gz",
            "truth_class_summaries.tsv",
            "simulation_failed.json",
        ):
            (self._unit_dir / relative).unlink(missing_ok=True)
        pairs = self._unit_dir / "pairs"
        if pairs.is_dir():
            shutil.rmtree(pairs)

    def _unit_atomic_temp_paths(self) -> list[Path]:
        return _unit_atomic_temp_paths(self._unit_dir)

    def _remove_unit_atomic_temps(self) -> list[str]:
        removed: list[str] = []
        for path in self._unit_atomic_temp_paths():
            path.unlink()
            removed.append(path.name)
        return removed

    def _write_base_attempts(self, attempts: pd.DataFrame) -> None:
        _atomic_frame(self._unit_dir / "attempts.tsv", attempts)

    def _reconcile_uncommitted_supervision(self) -> None:
        attempts = self._read_base_attempts()
        committed_count = len(attempts)
        observed_supervision = sorted(self.records)
        if observed_supervision != list(range(len(observed_supervision))):
            raise ValueError("Han v3 supervisor ledger is noncontiguous")
        if len(observed_supervision) < committed_count:
            raise ValueError("Han v3 supervisor ledger lacks a committed attempt")
        dangling = [
            attempt for attempt in observed_supervision if attempt >= committed_count
        ]
        for attempt in dangling:
            if attempt != len(attempts):
                raise ValueError("Han v3 dangling supervisor attempt is noncontiguous")
            record = dict(self.records[attempt])
            record["execution_id"] = self._execution_id(record)
            record["supervisor_status"] = "engine_complete_without_base_attempt_commit"
            record["lifecycle_phase"] = "engine_complete_without_base_attempt_commit"
            self.records[attempt] = record
            removed_atomic_temps = self._remove_unit_atomic_temps()
            elapsed = max(1e-9, float(record.get("supervisor_elapsed_seconds", 0.0)))
            durable_reconciliation = self._append_reconciliation(
                {
                    "schema": SUPERVISOR_RECONCILIATION_SCHEMA,
                    "reason": "engine_complete_without_base_attempt_commit",
                    "execution_id": self._execution_id(record),
                    "attempt_zero_based": attempt,
                    "simulation_seed": int(record["simulation_seed"]),
                    "prior_supervisor_pid": int(record["supervisor_pid"]),
                    "prior_child_pid": record.get("child_pid"),
                    "prior_child_pgid": record.get("child_pgid"),
                    "process_group_gone": bool(record.get("process_group_gone")),
                    "temp_directory_removed": bool(record.get("temp_tree_removed")),
                    "charged_elapsed_seconds": 0.0,
                    "transferred_to_attempt_elapsed_seconds": elapsed,
                    "atomic_temp_files_removed": json.dumps(
                        removed_atomic_temps, separators=(",", ":")
                    ),
                    "atomic_temp_file_count_removed": len(removed_atomic_temps),
                    "reconciled_unix_ns": time_ns(),
                }
            )
            durable_elapsed = float(
                durable_reconciliation["transferred_to_attempt_elapsed_seconds"]
            )
            if (
                float(durable_reconciliation["charged_elapsed_seconds"]) != 0
                or not math.isfinite(durable_elapsed)
                or durable_elapsed <= 0
            ):
                raise ValueError(
                    "Han v3 uncommitted-supervision reconciliation differs"
                )
            attempts = pd.concat(
                [
                    attempts,
                    pd.DataFrame(
                        [
                            self._interrupted_attempt_row(
                                record, elapsed_seconds=durable_elapsed
                            )
                        ]
                    ),
                ],
                ignore_index=True,
            )
        if dangling:
            self._persist()
            self._write_base_attempts(attempts)

    @property
    def interrupted_elapsed_seconds(self) -> float:
        if not self._reconciliation_path.is_file():
            return 0.0
        frame = pd.read_csv(self._reconciliation_path, sep="\t")
        if "charged_elapsed_seconds" not in frame:
            return 0.0
        values = pd.to_numeric(
            frame["charged_elapsed_seconds"], errors="coerce"
        ).fillna(0.0)
        if (values < 0).any():
            raise ValueError("Han v3 interrupted elapsed charge is negative")
        return float(values.sum())

    def _reconcile_interrupted_attempts(self) -> None:
        attempts = self._read_base_attempts()
        attempts_changed = False
        records_changed = False
        for active_path in sorted(self._unit_dir.glob(".han_v3_active_*.json")):
            active_sha256 = sha256_file(active_path)
            active = _read_active_supervision(active_path)
            attempt_record = dict(active.get("attempt_record", {}))
            attempt = int(attempt_record["attempt_zero_based"])
            seed = int(attempt_record["simulation_seed"])
            if attempt not in self._schedule.index or seed != int(
                self._schedule.loc[attempt, "simulation_seed"]
            ):
                raise ValueError("Han v3 active attempt differs from seed plan")
            persisted = self.records.get(attempt)
            if persisted is not None:
                if int(persisted["simulation_seed"]) != seed:
                    raise ValueError("Han v3 active/supervision seed differs")
                if self._execution_id(persisted) != self._execution_id(attempt_record):
                    raise ValueError("Han v3 active/supervision execution differs")
                attempt_record = {**attempt_record, **persisted}
                persisted_status = str(persisted.get("supervisor_status", ""))
                if persisted_status in {
                    "child_error",
                    "protocol_error",
                    "descendant_leak",
                }:
                    error_type = str(persisted.get("child_error_type", ""))
                    error_text = str(persisted.get("child_error", ""))
                    detail = ": ".join(
                        value for value in (error_type, error_text) if value
                    )
                    suffix = f" ({detail})" if detail else ""
                    raise RuntimeError(
                        "Han v3 persisted nonretryable supervisor failure: "
                        f"{persisted_status}{suffix}"
                    )
            completion_path = self._unit_dir / "simulation_complete.json"
            if completion_path.is_file() and attempt < len(attempts):
                # The atomic completion is the terminal unit transaction commit.
                # A worker can die after that commit but before the context
                # removes its pending marker.  The durable supervisor and
                # attempt rows are canonical, so inspect but never rewrite them.
                if persisted is None:
                    raise ValueError(
                        "Han v3 completed stale marker lacks persisted supervision"
                    )
                ledger_seed = int(attempts.iloc[attempt]["seed"])
                accepted_value = attempts.iloc[attempt].get("accepted", False)
                ledger_accepted = (
                    bool(accepted_value)
                    if isinstance(accepted_value, (bool, np.bool_))
                    else str(accepted_value).casefold() == "true"
                )
                expected_attempt_sha256 = str(
                    attempt_record.get("base_attempt_ledger_sha256", "")
                )
                if (
                    ledger_seed != seed
                    or not ledger_accepted
                    or (
                        expected_attempt_sha256
                        and expected_attempt_sha256
                        != sha256_file(self._unit_dir / "attempts.tsv")
                    )
                ):
                    raise ValueError(
                        "Han v3 completed stale marker differs from terminal attempt"
                    )
                if self._unit_atomic_temp_paths():
                    raise ValueError(
                        "Han v3 completed stale marker has atomic temp artifacts"
                    )
                active_path.unlink()
                continue
            temp_name = str(attempt_record.get("temp_directory", ""))
            temp_dir = (self._unit_dir / temp_name).resolve()
            if (
                temp_dir.parent != self._unit_dir.resolve()
                or not temp_dir.name.startswith(".han_v3_attempt_")
            ):
                raise ValueError("Han v3 interrupted temp-directory identity differs")
            session_identity_path = temp_dir / "session_process_identity.json"
            session_identity_deadline = monotonic() + max(
                1.0, 2.0 * self._term_grace_seconds
            )
            while (
                not session_identity_path.is_file()
                and monotonic() < session_identity_deadline
            ):
                sleep(SUPERVISOR_POLL_SECONDS)
            if session_identity_path.is_file():
                try:
                    session_identity = json.loads(
                        session_identity_path.read_text(encoding="utf-8")
                    )
                    session_pid = int(session_identity["session_pid"])
                    session_pgid = int(session_identity["session_pgid"])
                    session_sid = int(session_identity["session_sid"])
                    session_process_identity = _validate_linux_process_identity(
                        session_identity["session_process_identity"],
                        expected_pid=session_pid,
                    )
                except (
                    KeyError,
                    OSError,
                    TypeError,
                    ValueError,
                    json.JSONDecodeError,
                ) as error:
                    raise ValueError(
                        "Han v3 interrupted session identity is unreadable"
                    ) from error
                if (
                    session_pid != session_pgid
                    or session_pid != session_sid
                    or session_identity.get("parent_death_watchdog_armed") is not True
                ):
                    raise ValueError("Han v3 interrupted session identity differs")
                attempt_record["child_pid"] = session_pid
                attempt_record["child_pgid"] = session_pgid
                attempt_record["child_sid"] = session_sid
                attempt_record["child_process_identity"] = session_process_identity
                attempt_record["child_identity_observed"] = True
                attempt_record["parent_death_watchdog_armed"] = True
            engine_identity_path = temp_dir / "engine_process_identity.json"
            engine_identity_deadline = monotonic() + max(
                1.0, 2.0 * self._term_grace_seconds
            )
            while (
                not engine_identity_path.is_file()
                and attempt_record.get("child_pid") is not None
                and not _process_absent_or_zombie(int(attempt_record["child_pid"]))
                and monotonic() < engine_identity_deadline
            ):
                sleep(SUPERVISOR_POLL_SECONDS)
            if engine_identity_path.is_file():
                try:
                    engine_identity = json.loads(
                        engine_identity_path.read_text(encoding="utf-8")
                    )
                except (OSError, json.JSONDecodeError) as error:
                    raise ValueError(
                        "Han v3 interrupted engine identity is unreadable"
                    ) from error
                attempt_record["engine_wrapper_pid"] = int(
                    engine_identity["engine_pid"]
                )
                attempt_record["engine_wrapper_process_identity"] = engine_identity[
                    "engine_process_identity"
                ]
            if _same_linux_process(attempt_record.get("supervisor_process_identity")):
                raise ValueError("Han v3 active supervisor process is still alive")
            deadline = monotonic() + max(1.0, self._term_grace_seconds)
            child_identity = attempt_record.get("child_process_identity")
            child_pid = attempt_record.get("child_pid")
            while monotonic() < deadline:
                if child_pid is not None and not pd.isna(child_pid):
                    _reap_if_our_child(int(child_pid))
                child_live = (
                    child_pid is not None
                    and not pd.isna(child_pid)
                    and _same_linux_process(child_identity)
                    and not _process_absent_or_zombie(int(child_pid))
                )
                if not child_live:
                    engine_identity = attempt_record.get(
                        "engine_wrapper_process_identity"
                    )
                    engine_pid = attempt_record.get("engine_wrapper_pid")
                    engine_live = (
                        engine_pid is not None
                        and not pd.isna(engine_pid)
                        and _same_linux_process(engine_identity)
                        and not _process_absent_or_zombie(int(engine_pid))
                    )
                    group_live = (
                        attempt_record.get("child_pgid") is not None
                        and not pd.isna(attempt_record.get("child_pgid"))
                        and _process_group_has_live_members(
                            int(attempt_record["child_pgid"])
                        )
                    )
                    if not engine_live and not group_live:
                        break
                sleep(SUPERVISOR_POLL_SECONDS)
            if (
                child_pid is not None
                and not pd.isna(child_pid)
                and _same_linux_process(child_identity)
                and not _process_absent_or_zombie(int(child_pid))
            ):
                raise ValueError("Han v3 interrupted child cleanup is still pending")
            if (
                child_pid is not None
                and not pd.isna(child_pid)
                and not (_process_absent_or_zombie(int(child_pid)))
            ):
                raise ValueError("Han v3 interrupted child process is still alive")
            engine_pid = attempt_record.get("engine_wrapper_pid")
            engine_identity = attempt_record.get("engine_wrapper_process_identity")
            if (
                engine_pid is not None
                and not pd.isna(engine_pid)
                and _same_linux_process(engine_identity)
                and not _process_absent_or_zombie(int(engine_pid))
            ):
                raise ValueError("Han v3 interrupted engine wrapper is still alive")
            child_pgid = attempt_record.get("child_pgid")
            if (
                child_pgid is not None
                and not pd.isna(child_pgid)
                and _process_group_has_live_members(int(child_pgid))
            ):
                raise ValueError("Han v3 interrupted process group is still live")
            temp_existed = temp_dir.exists()
            if temp_existed:
                shutil.rmtree(temp_dir)
            if temp_dir.exists():
                raise RuntimeError("Han v3 interrupted temp-directory cleanup failed")
            started_ns = int(attempt_record.get("supervisor_started_unix_ns", 0))
            wall_elapsed = min(
                CUMULATIVE_TIMEOUT_SECONDS,
                max(1e-9, (time_ns() - started_ns) / 1e9)
                if started_ns > 0
                else max(
                    1e-9,
                    float(attempt_record.get("supervisor_elapsed_seconds", 0.0)),
                ),
            )
            attempt_record = self._normalize_interrupted_supervision(
                attempt_record,
                elapsed_seconds=wall_elapsed,
                temp_removed=True,
            )
            self.records[attempt] = attempt_record
            records_changed = True
            if attempt > len(attempts):
                raise ValueError("Han v3 interrupted attempt is noncontiguous")
            removed_atomic_temps = self._remove_unit_atomic_temps()
            execution_id = self._execution_id(attempt_record)
            prior_reconciliation = self._reconciliation_for_execution(execution_id)
            current_accepted = False
            if attempt < len(attempts):
                accepted_value = attempts.iloc[attempt].get("accepted", False)
                current_accepted = (
                    bool(accepted_value)
                    if isinstance(accepted_value, (bool, np.bool_))
                    else str(accepted_value).casefold() == "true"
                )
            if prior_reconciliation is not None:
                reason = str(prior_reconciliation["reason"])
            elif attempt == len(attempts):
                reason = "worker_interrupted_active_engine_call"
            elif current_accepted:
                reason = "accepted_attempt_interrupted_before_unit_completion"
            else:
                reason = "attempt_interrupted_after_base_commit"

            precommit_reason = reason == "worker_interrupted_active_engine_call"
            accepted_interruption = (
                reason == "accepted_attempt_interrupted_before_unit_completion"
            )
            supplemental_charge = 0.0
            transferred_elapsed = wall_elapsed if precommit_reason else 0.0
            if not precommit_reason:
                if attempt >= len(attempts):
                    raise ValueError(
                        "Han v3 post-commit reconciliation lacks base attempt"
                    )
                committed_elapsed = float(
                    pd.to_numeric(
                        pd.Series([attempts.iloc[attempt].get("elapsed_seconds")]),
                        errors="coerce",
                    ).iloc[0]
                )
                if not math.isfinite(committed_elapsed) or committed_elapsed < 0:
                    raise ValueError("Han v3 interrupted committed elapsed differs")
                commit_ns = int(attempt_record.get("base_attempt_committed_unix_ns", 0))
                supplemental_charge = (
                    max(1e-9, (time_ns() - commit_ns) / 1e9)
                    if commit_ns > 0
                    else max(1e-9, wall_elapsed - committed_elapsed)
                )
            if accepted_interruption:
                self._clean_incomplete_unit_artifacts()
            durable_reconciliation = self._append_reconciliation(
                {
                    "schema": SUPERVISOR_RECONCILIATION_SCHEMA,
                    "reason": reason,
                    "execution_id": execution_id,
                    "attempt_zero_based": attempt,
                    "simulation_seed": seed,
                    "active_record_path": active_path.name,
                    "active_record_sha256": active_sha256,
                    "prior_supervisor_pid": int(attempt_record["supervisor_pid"]),
                    "prior_child_pid": attempt_record.get("child_pid"),
                    "prior_child_pgid": attempt_record.get("child_pgid"),
                    "process_group_gone": True,
                    "temp_directory_existed": temp_existed,
                    "temp_directory_removed": True,
                    "charged_elapsed_seconds": supplemental_charge,
                    "transferred_to_attempt_elapsed_seconds": transferred_elapsed,
                    "atomic_temp_files_removed": json.dumps(
                        removed_atomic_temps, separators=(",", ":")
                    ),
                    "atomic_temp_file_count_removed": len(removed_atomic_temps),
                    "reconciled_unix_ns": time_ns(),
                }
            )
            durable_charge = float(
                pd.to_numeric(
                    pd.Series([durable_reconciliation["charged_elapsed_seconds"]]),
                    errors="coerce",
                ).iloc[0]
            )
            durable_transfer = float(
                pd.to_numeric(
                    pd.Series(
                        [
                            durable_reconciliation[
                                "transferred_to_attempt_elapsed_seconds"
                            ]
                        ]
                    ),
                    errors="coerce",
                ).iloc[0]
            )
            if (
                not math.isfinite(durable_charge)
                or durable_charge < 0
                or not math.isfinite(durable_transfer)
                or durable_transfer < 0
            ):
                raise ValueError("Han v3 reconciliation elapsed journal is corrupt")
            reason = str(durable_reconciliation["reason"])
            if reason == "worker_interrupted_active_engine_call":
                if durable_charge != 0 or durable_transfer <= 0:
                    raise ValueError("Han v3 pre-commit reconciliation journal differs")
                interrupted_row = self._interrupted_attempt_row(
                    attempt_record, elapsed_seconds=durable_transfer
                )
                if attempt == len(attempts):
                    attempts = pd.concat(
                        [attempts, pd.DataFrame([interrupted_row])],
                        ignore_index=True,
                    )
                else:
                    for column, value in interrupted_row.items():
                        attempts.loc[attempt, column] = value
                attempts_changed = True
            elif reason in {
                "attempt_interrupted_after_base_commit",
                "accepted_attempt_interrupted_before_unit_completion",
            }:
                if durable_transfer != 0 or durable_charge <= 0:
                    raise ValueError(
                        "Han v3 post-commit reconciliation journal differs"
                    )
                attempts.loc[attempt, "status"] = "timed_out"
                attempts.loc[attempt, "accepted"] = False
                attempts.loc[attempt, "panel_seed_used"] = False
                attempts.loc[attempt, "rejection_reason"] = (
                    "worker_interrupted_after_acceptance_before_unit_completion"
                    if reason == "accepted_attempt_interrupted_before_unit_completion"
                    else "worker_interrupted_after_base_attempt_commit"
                )
                attempts_changed = True
            else:
                raise ValueError("Han v3 active reconciliation reason differs")
            # All durable state needed for a restart is published before the
            # active transaction marker is removed. Each operation is atomic;
            # replay is coalesced by execution_id.
            self._persist()
            records_changed = False
            self._write_base_attempts(attempts)
            attempts_changed = False
            active_path.unlink()
        if attempts_changed:
            self._write_base_attempts(attempts)
        if records_changed:
            self._persist()

    def _persist(self) -> None:
        frame = pd.DataFrame([self.records[key] for key in sorted(self.records)])
        if frame.empty:
            self._ledger_path.unlink(missing_ok=True)
        else:
            _atomic_frame(self._ledger_path, frame)

    def mark_base_attempt_ledger_committed(self, frame: pd.DataFrame) -> None:
        """Advance pending engine calls only after base attempts.tsv is durable."""

        attempts_path = self._unit_dir / "attempts.tsv"
        if not attempts_path.is_file():
            raise RuntimeError("Han v3 base attempt ledger commit is absent")
        for raw in frame.to_dict(orient="records"):
            attempt = int(raw["attempt_zero_based"])
            active_path = self._unit_dir / f".han_v3_active_{attempt:03d}.json"
            if not active_path.is_file():
                continue
            if attempt not in self.records:
                raise RuntimeError("Han v3 committed attempt lacks supervision")
            accepted_value = raw.get("accepted", False)
            accepted = (
                bool(accepted_value)
                if isinstance(accepted_value, (bool, np.bool_))
                else str(accepted_value).casefold() == "true"
            )
            if not accepted:
                active_path.unlink()
                continue
            record = self.records[attempt]
            record["lifecycle_phase"] = (
                "accepted_attempt_committed_pending_unit_completion"
            )
            record["base_attempt_committed_unix_ns"] = time_ns()
            record["base_attempt_ledger_sha256"] = sha256_file(attempts_path)
            self.records[attempt] = record
            _write_active_supervision(
                active_path,
                {
                    "schema": SUPERVISOR_ACTIVE_SCHEMA,
                    "status": "active",
                    "attempt_record": record,
                },
            )
        self._persist()

    def finalize_unit_completion(self) -> None:
        active_paths = list(self._unit_dir.glob(".han_v3_active_*.json"))
        if not active_paths:
            return
        completion = self._unit_dir / "simulation_complete.json"
        if not completion.is_file():
            raise RuntimeError("Han v3 unit completion commit is absent")
        for active_path in active_paths:
            active_path.unlink()

    def simulate(self, *args: Any, **kwargs: Any) -> tskit.TreeSequence:
        timeout = getattr(_DRAW_TIMEOUT_LOCAL, "seconds", None)
        if timeout is None:
            raise RuntimeError("Han v3 engine call lacks supervised timeout context")
        logfile = Path(kwargs.get("logfile", ""))
        match = re.fullmatch(r"slim_draw_([0-9]+)\.csv", logfile.name)
        if match is None or logfile.parent.resolve() != self._unit_dir.resolve():
            raise ValueError("Han v3 engine call has an unexpected attempt logfile")
        attempt = int(match.group(1))
        if attempt not in self._schedule.index or attempt in self.records:
            raise ValueError("Han v3 engine call attempt is unallocated or duplicated")
        seed = int(kwargs.get("seed", -1))
        if seed != int(self._schedule.loc[attempt, "simulation_seed"]):
            raise ValueError("Han v3 supervised engine seed differs from its plan")
        temp_dir = Path(
            tempfile.mkdtemp(
                prefix=f".han_v3_attempt_{attempt:03d}_seed_{seed}_",
                dir=self._unit_dir,
            )
        )
        tree_path = temp_dir / "result.trees"
        active_path = self._unit_dir / f".han_v3_active_{attempt:03d}.json"
        if active_path.exists():
            raise ValueError("Han v3 active attempt record already exists")
        record: dict[str, Any] = {
            "schema": SUPERVISOR_LEDGER_SCHEMA,
            "attempt_zero_based": attempt,
            "simulation_seed": seed,
            "supervisor_started_unix_ns": time_ns(),
            "supervisor_pid": os.getpid(),
            "supervisor_process_identity": _linux_process_identity(os.getpid()),
            "timeout_seconds": float(timeout),
            "term_grace_seconds": self._term_grace_seconds,
            "temp_directory": temp_dir.name,
            "temp_tree_created": False,
            "temp_tree_loaded": False,
            "temp_tree_removed": False,
        }
        record["execution_id"] = self._execution_id(record)

        def persist_active() -> None:
            _write_active_supervision(
                active_path,
                {
                    "schema": SUPERVISOR_ACTIVE_SCHEMA,
                    "status": "active",
                    "attempt_record": record,
                },
            )

        persist_active()
        try:
            return _supervise_engine_call(
                self._engine,
                tuple(args),
                dict(kwargs),
                tree_path=tree_path,
                timeout_seconds=float(timeout),
                term_grace_seconds=self._term_grace_seconds,
                record=record,
                progress_callback=persist_active,
            )
        finally:
            shutil.rmtree(temp_dir, ignore_errors=False)
            record["temp_tree_removed"] = not temp_dir.exists()
            record["lifecycle_phase"] = "engine_returned_pending_base_attempt_commit"
            self.records[attempt] = record
            self._persist()
            persist_active()


@contextmanager
def patched_v3_executor(
    bundle: DirectV3Bundle, unit: Mapping[str, Any]
) -> Iterator[_SupervisedSlimEngine]:
    """Patch only one worker process and one immutable Han production unit."""

    record = simulation._validate_unit(unit)
    if (
        record["simulation_class"] != "selected"
        or record["demography_id"] != "ancient_eurasia_han_introgression"
    ):
        raise ValueError("Han v3 process-local patch accepts Han selected units only")
    schedule = _unit_seed_schedule(bundle, str(record["unit_id"]))
    expected_base = int(record["seed"])
    if set(schedule["production_unit_seed"].astype(int)) != {expected_base}:
        raise ValueError("Han v3 attempt schedule production seed differs")
    by_attempt = schedule.set_index("attempt_zero_based")
    binding = _unit_binding(bundle, record)
    unit_dir = bundle.campaign_dir / "work" / str(record["unit_id"])
    unit_dir.mkdir(parents=True, exist_ok=True)
    endpoint = float(schedule.iloc[0]["realized_selection_end_generations_ago"])
    endpoint_row = {
        "selection_coefficient": float(record["selection_coefficient"]),
        "target_allele_frequency": float(record["target_allele_frequency"]),
        "duration_multiplier": float(schedule.iloc[0]["duration_multiplier"]),
        "requested_selection_end_generations_ago": next(
            item["requested_selection_end_generations_ago"]
            for item in bundle.manifest["model"]["fixed_endpoints"]
            if math.isclose(
                float(item["selection_coefficient"]),
                float(record["selection_coefficient"]),
                abs_tol=1e-12,
            )
            and math.isclose(
                float(item["target_allele_frequency"]),
                float(record["target_allele_frequency"]),
                abs_tol=1e-12,
            )
        ),
        "realized_selection_end_generations_ago": endpoint,
        "slim_scaling_factor": SLIM_SCALING_FACTOR,
    }
    with _PATCH_LOCK:
        original_objects = simulation._selected_objects
        original_stable_seed = simulation._stable_seed
        original_atomic_frame = simulation._atomic_frame
        original_paths = simulation.IMPLEMENTATION_PATHS
        original_wall_clock_timeout = simulation._wall_clock_timeout
        original_get_engine = simulation.stdpopsim.get_engine
        supervised_engine = _SupervisedSlimEngine(
            original_get_engine("slim"), unit_dir=unit_dir, schedule=schedule
        )

        def selected_objects(
            selected_unit: Mapping[str, Any],
            repo_root: str | Path,
            slim_scaling_factor: float,
            selected_pool_diploids: int,
            *,
            han_selection_calibration: Mapping[str, Any] | None = None,
            han_selection_end_generations_ago: float | None = None,
        ):
            selected = simulation._validate_unit(selected_unit)
            if selected["unit_id"] != record["unit_id"]:
                raise ValueError("Han v3 worker received an unexpected production unit")
            if (
                not math.isclose(
                    float(slim_scaling_factor), SLIM_SCALING_FACTOR, abs_tol=1e-12
                )
                or int(selected_pool_diploids) != CANDIDATE_POOL_DIPLOIDS
                or han_selection_end_generations_ago is not None
                or dict(han_selection_calibration or {}) != binding
            ):
                raise ValueError("Han v3 selected-object authorization differs")
            model, sweep, samples, provenance = original_objects(
                selected,
                repo_root,
                slim_scaling_factor,
                selected_pool_diploids,
                han_selection_end_generations_ago=endpoint,
            )
            event_audit = fixed_v2._audit_fixed_events(sweep, endpoint_row)
            provenance = dict(provenance)
            provenance["han_direct_production_v3"] = {
                **binding,
                "direct_operational_estimand": True,
                "fixed_cessation_endpoint": endpoint,
                "event_audit": event_audit,
                "candidate_pool_af_gate": "target_plus_or_minus_0.025_inclusive",
                "uniform_exact_k_panel": True,
                "strict_selected_type1_identity": True,
                "normal_10mb_stdpopsim_processing": True,
                "conditional_inference": True,
            }
            return model, sweep, samples, provenance

        def stable_seed(base: int, label: str) -> int:
            if int(base) != expected_base:
                return original_stable_seed(base, label)
            match = re.fullmatch(r"selected(-panel)?:([0-9]+)", str(label))
            if match is None:
                return original_stable_seed(base, label)
            attempt = int(match.group(2))
            if attempt not in by_attempt.index:
                raise ValueError("Han v3 requested an unallocated attempt seed")
            column = "panel_seed" if match.group(1) else "simulation_seed"
            return int(by_attempt.loc[attempt, column])

        def atomic_frame(path: Path, frame: pd.DataFrame) -> None:
            if Path(path).name != "attempts.tsv" or "attempt_zero_based" not in frame:
                original_atomic_frame(path, frame)
                return
            enriched = frame.copy()
            simulation_seeds: list[int] = []
            panel_seeds: list[int] = []
            for index, raw_attempt in enumerate(enriched["attempt_zero_based"]):
                attempt = int(raw_attempt)
                if attempt not in by_attempt.index:
                    raise ValueError(
                        "Han v3 attempt ledger contains unallocated attempt"
                    )
                planned_sim = int(by_attempt.loc[attempt, "simulation_seed"])
                planned_panel = int(by_attempt.loc[attempt, "panel_seed"])
                if (
                    "seed" in enriched
                    and int(enriched.iloc[index]["seed"]) != planned_sim
                ):
                    raise ValueError("Han v3 attempt ledger simulation seed differs")
                simulation_seeds.append(planned_sim)
                panel_seeds.append(planned_panel)
            enriched["planned_simulation_seed"] = simulation_seeds
            enriched["panel_seed"] = panel_seeds
            enriched["panel_seed_used"] = (
                enriched.get("accepted", pd.Series(False, index=enriched.index))
                .fillna(False)
                .astype(bool)
            )
            supervisor_columns = (
                ("supervisor_schema", "schema"),
                ("supervisor_execution_id", "execution_id"),
                ("supervisor_status", "supervisor_status"),
                ("supervisor_pid", "supervisor_pid"),
                (
                    "supervisor_process_identity",
                    "supervisor_process_identity",
                ),
                ("supervisor_started_unix_ns", "supervisor_started_unix_ns"),
                ("child_pid", "child_pid"),
                ("child_pgid", "child_pgid"),
                ("child_sid", "child_sid"),
                ("child_process_identity", "child_process_identity"),
                ("child_identity_observed", "child_identity_observed"),
                ("engine_wrapper_pid", "engine_wrapper_pid"),
                (
                    "engine_wrapper_process_identity",
                    "engine_wrapper_process_identity",
                ),
                (
                    "parent_death_watchdog_armed",
                    "parent_death_watchdog_armed",
                ),
                ("supervisor_timeout_seconds", "timeout_seconds"),
                ("supervisor_term_grace_seconds", "term_grace_seconds"),
                ("supervisor_sigterm_sent", "sigterm_sent"),
                ("supervisor_sigkill_sent", "sigkill_sent"),
                ("direct_child_reaped", "direct_child_reaped"),
                ("adopted_children_reaped", "adopted_children_reaped"),
                ("process_group_gone", "process_group_gone"),
                ("attempt_temp_directory", "temp_directory"),
                ("attempt_temp_tree_created", "temp_tree_created"),
                ("attempt_temp_tree_loaded", "temp_tree_loaded"),
                ("attempt_temp_tree_removed", "temp_tree_removed"),
                ("supervisor_elapsed_seconds", "supervisor_elapsed_seconds"),
            )
            for output_column, record_column in supervisor_columns:
                values = []
                for raw_attempt in enriched["attempt_zero_based"]:
                    attempt = int(raw_attempt)
                    if attempt not in supervised_engine.records:
                        raise ValueError(
                            "Han v3 attempt lacks process-supervisor provenance"
                        )
                    values.append(supervised_engine.records[attempt].get(record_column))
                enriched[output_column] = values
            if "rejection_reason" not in enriched:
                enriched["rejection_reason"] = ""
            timed_out = (
                enriched.get("status", pd.Series("", index=enriched.index))
                .astype(str)
                .eq("timed_out")
            )
            blank_reason = enriched["rejection_reason"].fillna("").astype(str).eq("")
            enriched.loc[timed_out & blank_reason, "rejection_reason"] = "draw_timeout"
            enriched["attempt_ledger_schema"] = ATTEMPT_LEDGER_SCHEMA
            enriched["seed_plan_sha256"] = sha256_file(bundle.seed_plan_path)
            original_atomic_frame(path, enriched)
            supervised_engine.mark_base_attempt_ledger_committed(enriched)

        implementation_paths = _v3_implementation_paths(bundle)

        def get_engine(name: str, *args: Any, **kwargs: Any) -> Any:
            if name == "slim":
                if args or kwargs:
                    raise ValueError("Han v3 SLiM engine lookup arguments differ")
                return supervised_engine
            return original_get_engine(name, *args, **kwargs)

        simulation._selected_objects = selected_objects
        simulation._stable_seed = stable_seed
        simulation._atomic_frame = atomic_frame
        simulation.IMPLEMENTATION_PATHS = implementation_paths
        simulation._wall_clock_timeout = _supervised_draw_timeout
        simulation.stdpopsim.get_engine = get_engine
        unit_completed = False
        try:
            yield supervised_engine
            unit_completed = True
        finally:
            try:
                if unit_completed:
                    supervised_engine.finalize_unit_completion()
            finally:
                simulation._selected_objects = original_objects
                simulation._stable_seed = original_stable_seed
                simulation._atomic_frame = original_atomic_frame
                simulation.IMPLEMENTATION_PATHS = original_paths
                simulation._wall_clock_timeout = original_wall_clock_timeout
                simulation.stdpopsim.get_engine = original_get_engine


def _unit_row(bundle: DirectV3Bundle, unit_id: str) -> dict[str, Any]:
    rows = bundle.units[bundle.units["unit_id"].astype(str).eq(str(unit_id))]
    if len(rows) != 1:
        raise ValueError(f"Han v3 production unit lookup differs: {unit_id}")
    return rows.iloc[0].to_dict()


def validate_v3_cached_unit(
    bundle: DirectV3Bundle, unit: Mapping[str, Any]
) -> dict[str, Any]:
    """Read-only validation seam for one current v3 completion.

    The recorded focal-validation contract is checked directly, so this path
    never reconstructs or consults any mutable Han calibration object.
    """

    record = simulation._validate_unit(unit)
    if (
        record["simulation_class"] != "selected"
        or record["demography_id"] != "ancient_eurasia_han_introgression"
    ):
        raise ValueError("Han v3 cache seam accepts Han selected units only")
    unit_dir = bundle.campaign_dir / "work" / str(record["unit_id"])
    artifacts = simulation._completion_artifacts(unit_dir)
    try:
        completion = json.loads(artifacts.completion_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError("Han v3 cached completion is unreadable") from error
    contract = completion.get("contract")
    if (
        not isinstance(contract, Mapping)
        or completion.get("schema") != simulation.SCHEMA_VERSION
        or completion.get("status") != "complete"
        or completion.get("contract_sha256") != _canonical_sha256(contract)
        or _canonical_sha256(contract.get("unit")) != _canonical_sha256(record)
    ):
        raise ValueError("Han v3 cached completion contract/unit differs")
    parameters = contract.get("parameters", {})
    binding = _unit_binding(bundle, record)
    if (
        not math.isclose(
            float(parameters.get("slim_scaling_factor", math.nan)),
            SLIM_SCALING_FACTOR,
            abs_tol=1e-12,
        )
        or not math.isclose(
            float(parameters.get("slim_burn_in", math.nan)),
            SLIM_BURN_IN,
            abs_tol=1e-12,
        )
        or int(parameters.get("candidate_pool_diploids", -1)) != CANDIDATE_POOL_DIPLOIDS
        or _canonical_sha256(parameters.get("han_selection_calibration"))
        != _canonical_sha256(binding)
    ):
        raise ValueError("Han v3 cached completion parameter/authorization differs")
    implementation = contract.get("implementation", {})
    expected_paths = _v3_implementation_paths(bundle)
    expected_sources = {
        relative: sha256_file(bundle.repo_root / relative)
        for relative in expected_paths
    }
    slim_record = implementation.get("slim", {})
    if (
        implementation.get("sources") != expected_sources
        or implementation.get("eas_resource_sha256") != simulation.EAS_RESOURCE_SHA256
        or slim_record.get("path")
        != simulation._portable_path(bundle.slim_path, bundle.repo_root)
        or slim_record.get("sha256") != sha256_file(bundle.slim_path)
    ):
        raise ValueError("Han v3 cached completion implementation binding differs")
    validated = simulation._validate_completion(artifacts, contract)
    if not validated.cache_hit:
        raise ValueError("Han v3 base cache validation did not report a cache hit")
    provenance = completion.get("provenance", {})
    v3 = provenance.get("han_direct_production_v3", {})
    for key, value in binding.items():
        if _canonical_sha256(v3.get(key)) != _canonical_sha256(value):
            raise ValueError(
                "Han v3 cached completion provenance authorization differs"
            )
    if (
        v3.get("direct_operational_estimand") is not True
        or v3.get("normal_10mb_stdpopsim_processing") is not True
        or v3.get("strict_selected_type1_identity") is not True
    ):
        raise ValueError("Han v3 cached completion scientific semantics differ")
    _validate_recorded_event_audit(bundle, record, v3)
    attempt = _validate_attempt_ledger(bundle, record)
    search_budget = completion.get("search_budget", {})
    expected_remaining_cumulative = float(
        attempt["remaining_cumulative_timeout_seconds"]
    )
    if (
        int(search_budget.get("selected_max_draws", -1)) != ATTEMPTS_PER_UNIT
        or not math.isclose(
            float(search_budget.get("selected_draw_timeout_seconds", math.nan)),
            DRAW_TIMEOUT_SECONDS,
            abs_tol=1e-12,
        )
        or not math.isclose(
            float(search_budget.get("selected_cumulative_timeout_seconds", math.nan)),
            expected_remaining_cumulative,
            abs_tol=1e-12,
        )
    ):
        raise ValueError("Han v3 cached completion search budget differs")
    identity = provenance.get("selected_focal_identity", {})
    pool_gate = provenance.get("candidate_pool_frequency_gate", {})
    panel = provenance.get("sample_panel", {})
    if (
        identity.get("status") != "valid"
        or int(identity.get("observed_mutation_type", -1)) != 1
        or pool_gate.get("passed") is not True
        or int(panel.get("panel_seed", -1)) != int(attempt["accepted_panel_seed"])
    ):
        raise ValueError("Han v3 cached completion acceptance semantics differ")
    if int(completion.get("attempts_completed", -1)) != int(
        attempt["attempts_completed"]
    ):
        raise ValueError("Han v3 completion/attempt-ledger counts differ")
    return {
        "unit_id": str(record["unit_id"]),
        "selection_coefficient": float(record["selection_coefficient"]),
        "target_allele_frequency": float(record["target_allele_frequency"]),
        "replicate_index": int(record["replicate_index"]),
        "completion_path": artifacts.completion_path.relative_to(
            bundle.repo_root
        ).as_posix(),
        "completion_sha256": sha256_file(artifacts.completion_path),
        "contract_sha256": completion["contract_sha256"],
        "cache_validated_read_only": True,
        "authorization_validated": True,
        "current_sources_validated": True,
        "strict_type1_identity": True,
        "pool_gate_passed": True,
        "exact_k_panel_passed": True,
        **attempt,
    }


def aggregate_v3_cached_completions(
    bundle: DirectV3Bundle, stage: str = STAGE_ALL
) -> pd.DataFrame:
    """Aggregate only checksum- and semantics-valid current v3 completions."""

    units = _stage_units(bundle, stage)
    rows = [
        validate_v3_cached_unit(bundle, unit)
        for unit in units.to_dict(orient="records")
    ]
    return (
        pd.DataFrame(rows)
        .sort_values(
            ["selection_coefficient", "target_allele_frequency", "replicate_index"]
        )
        .reset_index(drop=True)
    )


@contextmanager
def patched_v3_cache_validation_dispatch(
    bundle: DirectV3Bundle,
    *,
    selected_bundle: Any | None = None,
) -> Iterator[None]:
    """Dispatch standard aggregation's Han cache check through frozen v3.

    Non-Han units are delegated byte-for-byte to the prior validator, so this
    can be nested inside the EAS selected-production context while neutral
    validation remains unchanged.
    """

    with _PATCH_LOCK:
        original = simulation._simulate_unit_unlocked

        def dispatch(
            unit: Mapping[str, Any],
            repo_root: str | Path,
            campaign_dir: str | Path,
            *args: Any,
            **kwargs: Any,
        ) -> simulation.SimulationArtifacts:
            record = simulation._validate_unit(unit)
            is_v3_han = (
                record["simulation_class"] == "selected"
                and record["demography_id"] == "ancient_eurasia_han_introgression"
            )
            if not is_v3_han:
                is_eas_selected = (
                    record["simulation_class"] == "selected"
                    and record["demography_id"] == "eas_phlash_median"
                )
                if selected_bundle is not None and is_eas_selected:
                    from . import focused_selected_production as selected_production

                    with selected_production.patched_selected_executor(selected_bundle):
                        return original(unit, repo_root, campaign_dir, *args, **kwargs)
                return original(unit, repo_root, campaign_dir, *args, **kwargs)
            if (
                Path(repo_root).resolve() != bundle.repo_root
                or Path(campaign_dir).resolve() != bundle.campaign_dir
            ):
                raise ValueError("Han v3 aggregate dispatch root/campaign differs")
            validate_v3_cached_unit(bundle, record)
            artifacts = simulation._completion_artifacts(
                bundle.campaign_dir / "work" / str(record["unit_id"])
            )
            return simulation.SimulationArtifacts(
                **{**artifacts.__dict__, "cache_hit": True}
            )

        simulation._simulate_unit_unlocked = dispatch
        try:
            yield
        finally:
            simulation._simulate_unit_unlocked = original


def aggregate_decoded_outputs_with_v3_dispatch(
    bundle: DirectV3Bundle,
    units: pd.DataFrame,
    *,
    selected_integration_manifest: str | Path | None = None,
) -> tuple[Path, Path, Path, Path, Path]:
    """Run the standard decoded-output aggregator with real v3 cache dispatch."""

    from . import focused_selection_campaign as focused_campaign

    selected_bundle: Any | None = None
    if selected_integration_manifest is not None:
        from . import focused_selected_production as selected_production

        selected_bundle = selected_production.load_selected_production_bundle(
            selected_integration_manifest,
            repo_root=bundle.repo_root,
            slim_path=bundle.slim_path,
        )
    with patched_v3_cache_validation_dispatch(bundle, selected_bundle=selected_bundle):
        return focused_campaign._aggregate_compact_decode(bundle.campaign_dir, units)


def _simulate_one(bundle: DirectV3Bundle, unit: Mapping[str, Any]) -> dict[str, Any]:
    started = perf_counter()
    try:
        binding = _unit_binding(bundle, unit)
        unit_dir = bundle.campaign_dir / "work" / str(unit["unit_id"])
        with (
            simulation.exclusive_unit_lock(
                unit_dir / ".simulation.lock",
                unit_id=str(unit["unit_id"]),
                phase="han_direct_production_v3",
            ),
            patched_v3_executor(bundle, unit) as supervised_engine,
        ):
            remaining_cumulative_seconds = (
                CUMULATIVE_TIMEOUT_SECONDS
                - supervised_engine.interrupted_elapsed_seconds
            )
            if remaining_cumulative_seconds <= 0:
                raise simulation.SimulationExhausted(
                    f"{unit['unit_id']} exhausted its cumulative selected-"
                    "simulation budget including interrupted supervised draws"
                )
            artifacts = simulation._simulate_unit_unlocked(
                unit,
                bundle.repo_root,
                bundle.campaign_dir,
                slim_path=bundle.slim_path,
                selected_max_draws=ATTEMPTS_PER_UNIT,
                slim_scaling_factor=SLIM_SCALING_FACTOR,
                slim_burn_in=SLIM_BURN_IN,
                selected_draw_timeout_seconds=DRAW_TIMEOUT_SECONDS,
                selected_cumulative_timeout_seconds=remaining_cumulative_seconds,
                han_selection_calibration=binding,
            )
        return {
            "unit_id": str(unit["unit_id"]),
            "status": "cached" if artifacts.cache_hit else "complete",
            "completion_path": artifacts.completion_path.relative_to(
                bundle.repo_root
            ).as_posix(),
            "elapsed_seconds": perf_counter() - started,
            "error_type": "",
            "error": "",
        }
    except Exception as error:  # noqa: BLE001 - worker must report all failures
        return {
            "unit_id": str(unit["unit_id"]),
            "status": "failed",
            "completion_path": "",
            "elapsed_seconds": perf_counter() - started,
            "error_type": type(error).__name__,
            "error": str(error),
        }


def _worker_initialize(manifest_path: str, repo_root: str, slim_path: str) -> None:
    global _WORKER_BUNDLE
    _WORKER_BUNDLE = load_bundle(
        manifest_path, repo_root=repo_root, slim_path=slim_path
    )


def _worker_simulate(unit: Mapping[str, Any]) -> dict[str, Any]:
    if _WORKER_BUNDLE is None:
        raise RuntimeError("Han v3 worker bundle is not initialized")
    return _simulate_one(_WORKER_BUNDLE, unit)


def _stage_units(bundle: DirectV3Bundle, stage: str) -> pd.DataFrame:
    if stage == STAGE_ALL:
        return bundle.units.copy()
    replicate_is_pilot = (
        bundle.units["replicate_index"].astype(int).eq(PILOT_REPLICATE_INDEX)
    )
    if stage == STAGE_PILOT:
        units = bundle.units[replicate_is_pilot].copy()
    elif stage == STAGE_SIMULATE:
        units = bundle.units[~replicate_is_pilot].copy()
    else:
        raise ValueError(f"Han v3 stage is unknown: {stage}")
    expected = 6 if stage == STAGE_PILOT else 54
    if len(units) != expected:
        raise ValueError(f"Han v3 {stage} unit count differs")
    return units


def _seed_plan_audit(bundle: DirectV3Bundle) -> dict[str, Any]:
    excluded = _validate_frozen_seed_registry(
        bundle.excluded_seed_registry, bundle.manifest["seed_registry"]
    )
    _validate_seed_plan(bundle.seed_plan, units=bundle.units, excluded=excluded)
    plan = bundle.seed_plan
    cell_counts = plan.groupby(
        ["selection_coefficient", "target_allele_frequency", "replicate_index"]
    ).size()
    if len(cell_counts) != 60 or set(cell_counts.astype(int)) != {ATTEMPTS_PER_UNIT}:
        raise ValueError("Han v3 audited schedule cell/unit attempts differ")
    fields = (
        "unit_id",
        "selection_coefficient",
        "target_allele_frequency",
        "replicate_index",
        "production_unit_seed",
        "stage",
        "attempt_zero_based",
        "realized_selection_end_generations_ago",
        "duration_multiplier",
    )
    return {
        "rows": len(plan),
        "units": int(plan["unit_id"].nunique()),
        "cells": int(
            plan[["selection_coefficient", "target_allele_frequency"]]
            .drop_duplicates()
            .shape[0]
        ),
        "attempts_per_unit": ATTEMPTS_PER_UNIT,
        "stage_row_counts": {
            str(key): int(value)
            for key, value in plan["stage"].value_counts().sort_index().items()
        },
        "schedule_contract_fields": list(fields),
        "schedule_contract_rows_sha256": _canonical_sha256(
            plan[list(fields)].to_dict(orient="records")
        ),
        "every_schedule_row_unit_cell_endpoint_stage_attempt_and_production_seed_validated": True,
        "attempt_and_panel_seeds_disjoint_from_frozen_registry": True,
    }


def _required_pilot_audit(bundle: DirectV3Bundle) -> dict[str, Any]:
    path = bundle.output_dir / "pilot_audit.json"
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(
            "Han v3 simulate stage requires a passed pilot audit"
        ) from error
    unhashed = dict(payload)
    digest = unhashed.pop("payload_sha256", None)
    if (
        payload.get("schema") != STAGE_AUDIT_SCHEMA
        or payload.get("stage") != STAGE_PILOT
        or payload.get("status") != "passed"
        or payload.get("manifest_sha256") != sha256_file(bundle.manifest_path)
        or payload.get("seed_plan_sha256") != sha256_file(bundle.seed_plan_path)
        or int(payload.get("accepted_units", -1)) != 6
        or payload.get("all_current_completions_cache_validated_read_only") is not True
        or payload.get("all_current_completion_budgets_match_authorization") is not True
        or payload.get("seed_plan_audit") != _seed_plan_audit(bundle)
        or digest != _canonical_sha256(unhashed)
    ):
        raise ValueError("Han v3 pilot audit is corrupt, stale, or not passed")
    current_rows = aggregate_v3_cached_completions(bundle, STAGE_PILOT).to_dict(
        orient="records"
    )
    if (
        _canonical_sha256(current_rows) != _canonical_sha256(payload.get("rows"))
        or len(current_rows) != 6
    ):
        raise ValueError(
            "Han v3 pilot audit no longer matches current validated completions"
        )
    return payload


def run_stage(bundle: DirectV3Bundle, stage: str, *, workers: int) -> pd.DataFrame:
    if stage == STAGE_SIMULATE:
        _required_pilot_audit(bundle)
    if not 1 <= int(workers) <= 24:
        raise ValueError("Han v3 workers must be in [1, 24]")
    units = _stage_units(bundle, stage)
    records = units.to_dict(orient="records")
    if int(workers) == 1:
        results = [_simulate_one(bundle, unit) for unit in records]
    else:
        results = []
        with ProcessPoolExecutor(
            max_workers=int(workers),
            initializer=_worker_initialize,
            initargs=(
                str(bundle.manifest_path),
                str(bundle.repo_root),
                str(bundle.slim_path),
            ),
        ) as pool:
            futures = {pool.submit(_worker_simulate, unit): unit for unit in records}
            for future in as_completed(futures):
                try:
                    results.append(future.result())
                except Exception as error:  # noqa: BLE001
                    unit = futures[future]
                    results.append(
                        {
                            "unit_id": str(unit["unit_id"]),
                            "status": "failed",
                            "completion_path": "",
                            "elapsed_seconds": 0.0,
                            "error_type": type(error).__name__,
                            "error": str(error),
                        }
                    )
    frame = pd.DataFrame(results).sort_values("unit_id").reset_index(drop=True)
    frame.insert(0, "stage", stage)
    frame["manifest_sha256"] = sha256_file(bundle.manifest_path)
    frame["seed_plan_sha256"] = sha256_file(bundle.seed_plan_path)
    frame["attempts_per_unit"] = ATTEMPTS_PER_UNIT
    frame["draw_timeout_seconds"] = DRAW_TIMEOUT_SECONDS
    frame["cumulative_timeout_seconds"] = CUMULATIVE_TIMEOUT_SECONDS
    frame["slim_scaling_factor"] = SLIM_SCALING_FACTOR
    frame["slim_burn_in"] = SLIM_BURN_IN
    snapshot = bundle.output_dir / f"{stage}_run_snapshot.tsv"
    _atomic_frame(snapshot, frame)
    return frame


def _validate_attempt_ledger(
    bundle: DirectV3Bundle, unit: Mapping[str, Any]
) -> dict[str, Any]:
    unit_id = str(unit["unit_id"])
    path = bundle.campaign_dir / "work" / unit_id / "attempts.tsv"
    if not path.is_file():
        raise ValueError(f"Han v3 attempt ledger is absent: {unit_id}")
    ledger = pd.read_csv(path, sep="\t")
    if not 1 <= len(ledger) <= ATTEMPTS_PER_UNIT:
        raise ValueError(f"Han v3 attempt ledger size differs: {unit_id}")
    required = {
        "attempt_zero_based",
        "seed",
        "planned_simulation_seed",
        "panel_seed",
        "panel_seed_used",
        "accepted",
        "rejection_reason",
        "status",
        "timeout_seconds",
        "elapsed_seconds",
        "attempt_ledger_schema",
        "seed_plan_sha256",
        "supervisor_schema",
        "supervisor_execution_id",
        "supervisor_status",
        "supervisor_pid",
        "supervisor_process_identity",
        "supervisor_started_unix_ns",
        "child_pid",
        "child_pgid",
        "child_sid",
        "child_process_identity",
        "child_identity_observed",
        "engine_wrapper_pid",
        "engine_wrapper_process_identity",
        "parent_death_watchdog_armed",
        "supervisor_timeout_seconds",
        "supervisor_term_grace_seconds",
        "supervisor_sigterm_sent",
        "supervisor_sigkill_sent",
        "direct_child_reaped",
        "adopted_children_reaped",
        "process_group_gone",
        "attempt_temp_directory",
        "attempt_temp_tree_created",
        "attempt_temp_tree_loaded",
        "attempt_temp_tree_removed",
        "supervisor_elapsed_seconds",
    }
    if required.difference(ledger.columns):
        raise ValueError(f"Han v3 attempt ledger lacks v3 provenance: {unit_id}")
    for column, minimum in {
        "attempt_zero_based": 0,
        "seed": 1,
        "planned_simulation_seed": 1,
        "panel_seed": 1,
        "supervisor_pid": 1,
        "supervisor_started_unix_ns": 1,
        "adopted_children_reaped": 0,
    }.items():
        ledger[column] = _exact_integer_series(
            ledger, column, label="attempt ledger", minimum=minimum
        )
    attempts = list(ledger["attempt_zero_based"])
    if attempts != list(range(len(ledger))):
        raise ValueError(f"Han v3 attempt ledger is noncontiguous: {unit_id}")
    schedule = _unit_seed_schedule(bundle, unit_id).iloc[: len(ledger)]
    if (
        not np.array_equal(
            ledger["seed"].astype(int).to_numpy(),
            schedule["simulation_seed"].astype(int).to_numpy(),
        )
        or not np.array_equal(
            ledger["planned_simulation_seed"].astype(int).to_numpy(),
            schedule["simulation_seed"].astype(int).to_numpy(),
        )
        or not np.array_equal(
            ledger["panel_seed"].astype(int).to_numpy(),
            schedule["panel_seed"].astype(int).to_numpy(),
        )
        or set(ledger["attempt_ledger_schema"].astype(str)) != {ATTEMPT_LEDGER_SCHEMA}
        or set(ledger["seed_plan_sha256"].astype(str))
        != {sha256_file(bundle.seed_plan_path)}
    ):
        raise ValueError(f"Han v3 attempt ledger seed provenance differs: {unit_id}")

    def boolean_series(column: str) -> pd.Series:
        values = ledger[column]
        if values.dtype == bool:
            return values
        normalized = values.fillna("").astype(str).str.casefold()
        if not normalized.isin({"true", "false"}).all():
            raise ValueError(
                f"Han v3 attempt ledger boolean {column} differs: {unit_id}"
            )
        return normalized.eq("true")

    accepted = boolean_series("accepted")
    panel_seed_used = boolean_series("panel_seed_used")
    direct_reaped = boolean_series("direct_child_reaped")
    group_gone = boolean_series("process_group_gone")
    temp_removed = boolean_series("attempt_temp_tree_removed")
    watchdog_armed = boolean_series("parent_death_watchdog_armed")
    sigterm_sent = boolean_series("supervisor_sigterm_sent")
    sigkill_sent = boolean_series("supervisor_sigkill_sent")
    child_identity_observed = boolean_series("child_identity_observed")
    if not panel_seed_used.equals(accepted):
        raise ValueError(f"Han v3 panel-seed-used flags differ: {unit_id}")
    if int(accepted.sum()) != 1 or not bool(accepted.iloc[-1]):
        raise ValueError(
            f"Han v3 attempt ledger lacks one terminal acceptance: {unit_id}"
        )
    if (
        not ledger.loc[~accepted, "rejection_reason"]
        .fillna("")
        .astype(str)
        .str.len()
        .gt(0)
        .all()
    ):
        raise ValueError(f"Han v3 rejection reason is not persisted: {unit_id}")
    child_pid = pd.to_numeric(ledger["child_pid"], errors="coerce")
    child_pgid = pd.to_numeric(ledger["child_pgid"], errors="coerce")
    child_sid = pd.to_numeric(ledger["child_sid"], errors="coerce")
    observed_child_values = pd.concat([child_pid, child_pgid, child_sid], axis=1)[
        child_identity_observed
    ]
    unobserved_child_values = pd.concat([child_pid, child_pgid, child_sid], axis=1)[
        ~child_identity_observed
    ]
    raw_unobserved_child_values = ledger.loc[
        ~child_identity_observed, ["child_pid", "child_pgid", "child_sid"]
    ]
    if (
        observed_child_values.isna().any().any()
        or not np.isfinite(observed_child_values.astype(float)).all().all()
        or (observed_child_values < 1).any().any()
        or not np.equal(
            observed_child_values.astype(float),
            np.rint(observed_child_values.astype(float)),
        )
        .all()
        .all()
        or unobserved_child_values.notna().any().any()
        or any(
            not bool(pd.isna(value))
            for value in raw_unobserved_child_values.to_numpy().ravel()
        )
    ):
        raise ValueError(f"Han v3 child PID provenance differs: {unit_id}")
    if (
        set(ledger["supervisor_schema"].astype(str)) != {SUPERVISOR_LEDGER_SCHEMA}
        or not np.array_equal(
            child_pid[child_identity_observed].astype(int),
            child_pgid[child_identity_observed].astype(int),
        )
        or not np.array_equal(
            child_pid[child_identity_observed].astype(int),
            child_sid[child_identity_observed].astype(int),
        )
        or not direct_reaped[child_identity_observed].all()
        or not group_gone.all()
        or not temp_removed.all()
        or not watchdog_armed[child_identity_observed].all()
        or not np.allclose(
            ledger["supervisor_timeout_seconds"].astype(float),
            ledger["timeout_seconds"].astype(float),
            rtol=0.0,
            atol=1e-9,
        )
        or set(ledger["supervisor_term_grace_seconds"].astype(float))
        != {SUPERVISOR_TERM_GRACE_SECONDS}
        or (ledger["adopted_children_reaped"].astype(int) < 0).any()
        or (ledger["supervisor_elapsed_seconds"].astype(float) <= 0).any()
        or (ledger["supervisor_started_unix_ns"].astype(int) <= 0).any()
    ):
        raise ValueError(f"Han v3 process-supervisor provenance differs: {unit_id}")
    if (
        direct_reaped[~child_identity_observed].any()
        or watchdog_armed[~child_identity_observed].any()
        or sigterm_sent[~child_identity_observed].any()
        or sigkill_sent[~child_identity_observed].any()
    ):
        raise ValueError(f"Han v3 unobserved child lifecycle was fabricated: {unit_id}")
    wrapper_pid = pd.to_numeric(ledger["engine_wrapper_pid"], errors="coerce")
    wrapper_observed = wrapper_pid.notna()
    if (
        (wrapper_pid[wrapper_observed] < 1).any()
        or not np.isfinite(wrapper_pid[wrapper_observed].astype(float)).all()
        or not np.equal(
            wrapper_pid[wrapper_observed].astype(float),
            np.rint(wrapper_pid[wrapper_observed].astype(float)),
        ).all()
    ):
        raise ValueError(f"Han v3 engine-wrapper PID differs: {unit_id}")
    supervisor_status = ledger["supervisor_status"].astype(str)
    completed_supervisor = supervisor_status.eq("complete")
    status = ledger["status"].astype(str)
    timed_out = status.eq("timed_out")
    complete = status.eq("complete")
    rejection_reason = ledger["rejection_reason"].fillna("").astype(str)
    scientific_rejections = {
        "selected_focal_mutation_absent",
        "selected_focal_mutation_absent_background_only",
        "candidate_pool_af_outside_prespecified_band",
        "candidate_pool_has_no_exact_af_genotype_qc_panel",
    }
    timeout_rejections = {
        "draw_timeout",
        "worker_interrupted_before_base_attempt_commit",
        "worker_interrupted_after_base_attempt_commit",
        "worker_interrupted_after_acceptance_before_unit_completion",
    }
    if (
        not status.isin({"complete", "timed_out"}).all()
        or not (complete[accepted] & rejection_reason[accepted].eq("")).all()
        or not rejection_reason[complete & ~accepted].isin(scientific_rejections).all()
        or not rejection_reason[timed_out].isin(timeout_rejections).all()
        or accepted[timed_out].any()
    ):
        raise ValueError(f"Han v3 base attempt status matrix differs: {unit_id}")
    allowed_timeout_supervisor = supervisor_status.isin(
        {
            "timed_out",
            "timed_out_after_terminal_message",
            "interrupted_by_worker_death",
            "startup_interrupted_before_child_observed",
            "engine_complete_without_base_attempt_commit",
        }
    )
    if (
        not (
            completed_supervisor.eq(~timed_out)
            & (completed_supervisor | allowed_timeout_supervisor)
        ).all()
        or not ledger.loc[completed_supervisor, "attempt_temp_tree_created"]
        .pipe(
            lambda values: (
                values
                if values.dtype == bool
                else values.astype(str).str.casefold().eq("true")
            )
        )
        .all()
        or not ledger.loc[completed_supervisor, "attempt_temp_tree_loaded"]
        .pipe(
            lambda values: (
                values
                if values.dtype == bool
                else values.astype(str).str.casefold().eq("true")
            )
        )
        .all()
    ):
        raise ValueError(f"Han v3 supervised attempt lifecycle differs: {unit_id}")
    supervision_path = path.with_name("supervision.tsv")
    if not supervision_path.is_file():
        raise ValueError(f"Han v3 supervisor ledger is absent: {unit_id}")
    supervision = pd.read_csv(supervision_path, sep="\t")
    if {
        "schema",
        "execution_id",
        "attempt_zero_based",
        "simulation_seed",
    }.difference(supervision.columns):
        raise ValueError(f"Han v3 supervisor ledger schema differs: {unit_id}")
    supervision["attempt_zero_based"] = _exact_integer_series(
        supervision,
        "attempt_zero_based",
        label="supervisor ledger",
        minimum=0,
    )
    supervision["simulation_seed"] = _exact_integer_series(
        supervision, "simulation_seed", label="supervisor ledger", minimum=1
    )
    if (
        len(supervision) != len(ledger)
        or list(supervision["attempt_zero_based"].astype(int)) != attempts
        or set(supervision["schema"].astype(str)) != {SUPERVISOR_LEDGER_SCHEMA}
        or not np.array_equal(
            supervision["simulation_seed"].astype(int).to_numpy(),
            schedule["simulation_seed"].astype(int).to_numpy(),
        )
        or supervision["execution_id"].fillna("").astype(str).eq("").any()
        or supervision["execution_id"].astype(str).duplicated().any()
        or not supervision["execution_id"]
        .astype(str)
        .str.fullmatch(r"[0-9a-f]{64}")
        .all()
    ):
        raise ValueError(f"Han v3 supervisor/attempt ledgers differ: {unit_id}")
    base_elapsed = pd.to_numeric(ledger["elapsed_seconds"], errors="coerce")
    supervisor_elapsed = pd.to_numeric(
        ledger["supervisor_elapsed_seconds"], errors="coerce"
    )
    if (
        base_elapsed.isna().any()
        or supervisor_elapsed.isna().any()
        or not np.isfinite(base_elapsed.astype(float)).all()
        or not np.isfinite(supervisor_elapsed.astype(float)).all()
        or (base_elapsed < 0).any()
        or (supervisor_elapsed <= 0).any()
        or (base_elapsed + 1e-9 < supervisor_elapsed).any()
    ):
        raise ValueError(f"Han v3 committed elapsed accounting differs: {unit_id}")
    copied_supervisor_fields = {
        "schema": "supervisor_schema",
        "execution_id": "supervisor_execution_id",
        "supervisor_status": "supervisor_status",
        "supervisor_pid": "supervisor_pid",
        "supervisor_process_identity": "supervisor_process_identity",
        "supervisor_started_unix_ns": "supervisor_started_unix_ns",
        "child_pid": "child_pid",
        "child_pgid": "child_pgid",
        "child_sid": "child_sid",
        "child_process_identity": "child_process_identity",
        "child_identity_observed": "child_identity_observed",
        "engine_wrapper_pid": "engine_wrapper_pid",
        "engine_wrapper_process_identity": "engine_wrapper_process_identity",
        "parent_death_watchdog_armed": "parent_death_watchdog_armed",
        "timeout_seconds": "supervisor_timeout_seconds",
        "term_grace_seconds": "supervisor_term_grace_seconds",
        "sigterm_sent": "supervisor_sigterm_sent",
        "sigkill_sent": "supervisor_sigkill_sent",
        "direct_child_reaped": "direct_child_reaped",
        "adopted_children_reaped": "adopted_children_reaped",
        "process_group_gone": "process_group_gone",
        "temp_directory": "attempt_temp_directory",
        "temp_tree_created": "attempt_temp_tree_created",
        "temp_tree_loaded": "attempt_temp_tree_loaded",
        "temp_tree_removed": "attempt_temp_tree_removed",
        "supervisor_elapsed_seconds": "supervisor_elapsed_seconds",
    }

    def normalized_values(values: pd.Series) -> list[str]:
        return [
            "" if pd.isna(value) else str(value).strip().casefold() for value in values
        ]

    for source_column, attempt_column in copied_supervisor_fields.items():
        if source_column not in supervision.columns:
            raise ValueError(
                f"Han v3 supervisor field {source_column} is absent: {unit_id}"
            )
        if normalized_values(supervision[source_column]) != normalized_values(
            ledger[attempt_column]
        ):
            raise ValueError(
                f"Han v3 supervisor field {source_column} differs: {unit_id}"
            )
    reconciliation_path = path.with_name("supervision_reconciliation.tsv")
    interrupted_charge = 0.0
    reconciliation_sha256 = None
    if reconciliation_path.is_file():
        reconciliation = pd.read_csv(reconciliation_path, sep="\t")
        required_reconciliation = {
            "schema",
            "reason",
            "execution_id",
            "attempt_zero_based",
            "simulation_seed",
            "process_group_gone",
            "temp_directory_removed",
            "charged_elapsed_seconds",
            "transferred_to_attempt_elapsed_seconds",
            "atomic_temp_files_removed",
            "atomic_temp_file_count_removed",
            "reconciled_unix_ns",
        }
        if required_reconciliation.difference(reconciliation.columns):
            raise ValueError(f"Han v3 reconciliation fields differ: {unit_id}")
        if set(reconciliation["schema"].astype(str)) != {
            SUPERVISOR_RECONCILIATION_SCHEMA
        }:
            raise ValueError(f"Han v3 reconciliation schema differs: {unit_id}")
        charges = pd.to_numeric(
            reconciliation["charged_elapsed_seconds"], errors="coerce"
        )
        transferred = pd.to_numeric(
            reconciliation["transferred_to_attempt_elapsed_seconds"],
            errors="coerce",
        )
        reconciliation["attempt_zero_based"] = _exact_integer_series(
            reconciliation,
            "attempt_zero_based",
            label="reconciliation ledger",
            minimum=0,
        )
        reconciliation["simulation_seed"] = _exact_integer_series(
            reconciliation,
            "simulation_seed",
            label="reconciliation ledger",
            minimum=1,
        )
        reconciliation["reconciled_unix_ns"] = _exact_integer_series(
            reconciliation,
            "reconciled_unix_ns",
            label="reconciliation ledger",
            minimum=1,
        )
        reconciliation["atomic_temp_file_count_removed"] = _exact_integer_series(
            reconciliation,
            "atomic_temp_file_count_removed",
            label="reconciliation ledger",
            minimum=0,
        )
        reasons = set(reconciliation["reason"].astype(str))
        execution_ids = reconciliation["execution_id"].fillna("").astype(str)
        reconciliation_attempts = reconciliation["attempt_zero_based"]
        expected_reconciliation_seeds = [
            int(_unit_seed_schedule(bundle, unit_id).iloc[attempt]["simulation_seed"])
            for attempt in reconciliation_attempts
        ]

        def reconciliation_bool(column: str) -> pd.Series:
            values = reconciliation[column]
            if values.dtype == bool:
                return values
            normalized = values.fillna("").astype(str).str.casefold()
            if not normalized.isin({"true", "false"}).all():
                raise ValueError(
                    f"Han v3 reconciliation boolean {column} differs: {unit_id}"
                )
            return normalized.eq("true")

        if (
            charges.isna().any()
            or not np.isfinite(charges.astype(float)).all()
            or (charges < 0).any()
            or transferred.isna().any()
            or not np.isfinite(transferred.astype(float)).all()
            or (transferred < 0).any()
            or not reasons.issubset(
                {
                    "worker_interrupted_active_engine_call",
                    "engine_complete_without_base_attempt_commit",
                    "attempt_interrupted_after_base_commit",
                    "accepted_attempt_interrupted_before_unit_completion",
                }
            )
            or (reconciliation_attempts < 0).any()
            or (reconciliation_attempts >= len(ledger)).any()
            or reconciliation["reconciled_unix_ns"].duplicated().any()
            or execution_ids.eq("").any()
            or execution_ids.duplicated().any()
            or not execution_ids.str.fullmatch(r"[0-9a-f]{64}").all()
            or reconciliation_attempts.duplicated().any()
            or list(reconciliation["simulation_seed"]) != expected_reconciliation_seeds
            or not reconciliation_bool("process_group_gone").all()
            or not reconciliation_bool("temp_directory_removed").all()
        ):
            raise ValueError(f"Han v3 interrupted charge differs: {unit_id}")
        supervisor_execution_by_attempt = {
            int(row["attempt_zero_based"]): str(row["execution_id"])
            for row in supervision.to_dict(orient="records")
        }
        special_dispositions: set[int] = set()
        for row_index, row in reconciliation.iterrows():
            attempt = int(row["attempt_zero_based"])
            reason = str(row["reason"])
            attempt_row = ledger.iloc[attempt]
            supervisor_row = supervision.iloc[attempt]
            rejection_reason = str(attempt_row.get("rejection_reason", ""))
            if (
                str(row["execution_id"]) != supervisor_execution_by_attempt.get(attempt)
                or str(attempt_row["supervisor_execution_id"])
                != str(row["execution_id"])
                or str(attempt_row["status"]) != "timed_out"
                or bool(accepted.iloc[attempt])
                or bool(panel_seed_used.iloc[attempt])
            ):
                raise ValueError(
                    f"Han v3 reconciliation disposition differs: {unit_id}"
                )
            charge = float(charges.iloc[row_index])
            transfer = float(transferred.iloc[row_index])
            supervisor_lifecycle = str(supervisor_row["supervisor_status"])
            try:
                removed_atomic_temps = json.loads(str(row["atomic_temp_files_removed"]))
            except json.JSONDecodeError as error:
                raise ValueError(
                    f"Han v3 reconciliation atomic-temp record differs: {unit_id}"
                ) from error
            if (
                not isinstance(removed_atomic_temps, list)
                or any(not isinstance(name, str) for name in removed_atomic_temps)
                or len(removed_atomic_temps) != len(set(removed_atomic_temps))
                or len(removed_atomic_temps)
                != int(row["atomic_temp_file_count_removed"])
                or any(
                    not _is_unit_atomic_temp_name(name) for name in removed_atomic_temps
                )
            ):
                raise ValueError(
                    f"Han v3 reconciliation atomic-temp record differs: {unit_id}"
                )
            if reason == "worker_interrupted_active_engine_call":
                valid = (
                    rejection_reason == "worker_interrupted_before_base_attempt_commit"
                    and supervisor_lifecycle
                    in {
                        "interrupted_by_worker_death",
                        "startup_interrupted_before_child_observed",
                    }
                    and transfer > 0
                    and charge == 0
                    and math.isclose(
                        float(base_elapsed.iloc[attempt]),
                        transfer,
                        rel_tol=0.0,
                        abs_tol=1e-9,
                    )
                )
            elif reason == "engine_complete_without_base_attempt_commit":
                valid = (
                    rejection_reason == "worker_interrupted_before_base_attempt_commit"
                    and supervisor_lifecycle
                    == "engine_complete_without_base_attempt_commit"
                    and transfer > 0
                    and charge == 0
                    and math.isclose(
                        float(base_elapsed.iloc[attempt]),
                        transfer,
                        rel_tol=0.0,
                        abs_tol=1e-9,
                    )
                )
            elif reason == "attempt_interrupted_after_base_commit":
                valid = (
                    rejection_reason == "worker_interrupted_after_base_attempt_commit"
                    and supervisor_lifecycle == "interrupted_by_worker_death"
                    and transfer == 0
                    and charge > 0
                )
            else:
                valid = (
                    rejection_reason
                    == "worker_interrupted_after_acceptance_before_unit_completion"
                    and supervisor_lifecycle == "interrupted_by_worker_death"
                    and transfer == 0
                    and charge > 0
                )
            if not valid:
                raise ValueError(
                    f"Han v3 reconciliation reason matrix differs: {unit_id}"
                )
            special_dispositions.add(attempt)
        expected_special_dispositions = {
            index
            for index, value in enumerate(
                ledger["rejection_reason"].fillna("").astype(str)
            )
            if value
            in {
                "worker_interrupted_before_base_attempt_commit",
                "worker_interrupted_after_base_attempt_commit",
                "worker_interrupted_after_acceptance_before_unit_completion",
            }
        }
        if special_dispositions != expected_special_dispositions:
            raise ValueError(
                f"Han v3 reconciliation terminal dispositions differ: {unit_id}"
            )
        interrupted_charge = float(charges.sum())
        reconciliation_sha256 = sha256_file(reconciliation_path)
    elif (
        ledger["rejection_reason"]
        .fillna("")
        .astype(str)
        .isin(
            {
                "worker_interrupted_before_base_attempt_commit",
                "worker_interrupted_after_base_attempt_commit",
                "worker_interrupted_after_acceptance_before_unit_completion",
            }
        )
        .any()
    ):
        raise ValueError(
            f"Han v3 interrupted disposition lacks reconciliation: {unit_id}"
        )
    active_markers = list(path.parent.glob(".han_v3_active_*.json"))
    attempt_scratch = list(path.parent.glob(".han_v3_attempt_*"))
    atomic_temps = _unit_atomic_temp_paths(path.parent)
    if active_markers or attempt_scratch or atomic_temps:
        raise ValueError(f"Han v3 active/scratch attempt state remains: {unit_id}")
    startup_cleanup_path = path.with_name("supervision_startup_cleanup.tsv")
    startup_cleanup_count = 0
    startup_cleanup_sha256 = None
    if startup_cleanup_path.is_file():
        startup_cleanup = pd.read_csv(startup_cleanup_path, sep="\t")
        required_startup_cleanup = {
            "schema",
            "cleanup_id",
            "artifact_name",
            "artifact_kind",
            "inventory_sha256",
            "status",
            "discovered_unix_ns",
            "removed_unix_ns",
        }
        if required_startup_cleanup.difference(startup_cleanup.columns):
            raise ValueError(f"Han v3 startup-cleanup fields differ: {unit_id}")
        discovered = _exact_integer_series(
            startup_cleanup,
            "discovered_unix_ns",
            label="startup-cleanup ledger",
            minimum=1,
        )
        removed = _exact_integer_series(
            startup_cleanup,
            "removed_unix_ns",
            label="startup-cleanup ledger",
            minimum=1,
        )
        expected_cleanup_ids = [
            _canonical_sha256(
                {
                    "artifact_name": str(row["artifact_name"]),
                    "artifact_kind": str(row["artifact_kind"]),
                    "inventory_sha256": str(row["inventory_sha256"]),
                }
            )
            for row in startup_cleanup.to_dict(orient="records")
        ]
        if (
            set(startup_cleanup["schema"].astype(str))
            != {SUPERVISOR_STARTUP_CLEANUP_SCHEMA}
            or set(startup_cleanup["status"].astype(str)) != {"removed"}
            or not startup_cleanup["artifact_kind"]
            .astype(str)
            .isin(
                {
                    "attempt_scratch_directory",
                    "active_atomic_temp",
                    "known_atomic_temp",
                }
            )
            .all()
            or startup_cleanup["cleanup_id"].astype(str).duplicated().any()
            or not startup_cleanup["cleanup_id"]
            .astype(str)
            .str.fullmatch(r"[0-9a-f]{64}")
            .all()
            or list(startup_cleanup["cleanup_id"].astype(str)) != expected_cleanup_ids
            or not startup_cleanup["inventory_sha256"]
            .astype(str)
            .str.fullmatch(r"[0-9a-f]{64}")
            .all()
            or (removed < discovered).any()
            or any(
                (path.parent / str(name)).exists()
                for name in startup_cleanup["artifact_name"]
            )
        ):
            raise ValueError(f"Han v3 startup-cleanup provenance differs: {unit_id}")
        startup_cleanup_count = len(startup_cleanup)
        startup_cleanup_sha256 = sha256_file(startup_cleanup_path)
    if simulation.unit_lock_is_held(path.parent / ".simulation.lock"):
        raise ValueError(f"Han v3 completed unit lock is still held: {unit_id}")
    for index, value in enumerate(ledger["supervisor_process_identity"]):
        _validate_linux_process_identity(
            value, expected_pid=int(ledger.iloc[index]["supervisor_pid"])
        )
    for index, value in enumerate(ledger["child_process_identity"]):
        if bool(child_identity_observed.iloc[index]):
            identity = _validate_linux_process_identity(
                value, expected_pid=int(child_pid.iloc[index])
            )
            if _same_linux_process(identity):
                raise ValueError(
                    f"Han v3 supervised descendant identity is still live: {unit_id}"
                )
        elif _parse_linux_process_identity(value) is not None or (
            isinstance(value, str) and value.strip()
        ):
            raise ValueError(f"Han v3 startup-interrupted identity differs: {unit_id}")
    for index, value in enumerate(ledger["engine_wrapper_process_identity"]):
        if bool(wrapper_observed.iloc[index]):
            identity = _validate_linux_process_identity(
                value, expected_pid=int(wrapper_pid.iloc[index])
            )
            if _same_linux_process(identity):
                raise ValueError(
                    f"Han v3 supervised wrapper identity is still live: {unit_id}"
                )
        elif _parse_linux_process_identity(value) is not None or (
            isinstance(value, str) and value.strip()
        ):
            raise ValueError(
                f"Han v3 unobserved wrapper has process identity: {unit_id}"
            )
    if interrupted_charge >= CUMULATIVE_TIMEOUT_SECONDS:
        raise ValueError(
            f"Han v3 interrupted charge exhausted authorization: {unit_id}"
        )
    committed_elapsed = float(base_elapsed.sum())
    cumulative_elapsed = committed_elapsed + interrupted_charge
    # A terminal accepted draw may cross the budget by only its authorized
    # effective per-draw timeout; no subsequent attempt may start afterward.
    if cumulative_elapsed > CUMULATIVE_TIMEOUT_SECONDS + DRAW_TIMEOUT_SECONDS + 1e-6:
        raise ValueError(f"Han v3 cumulative elapsed accounting differs: {unit_id}")
    return {
        "path": path.relative_to(bundle.repo_root).as_posix(),
        "sha256": sha256_file(path),
        "attempts_completed": len(ledger),
        "rejections": len(ledger) - 1,
        "accepted_attempt_zero_based": int(ledger.iloc[-1]["attempt_zero_based"]),
        "accepted_simulation_seed": int(ledger.iloc[-1]["seed"]),
        "accepted_panel_seed": int(ledger.iloc[-1]["panel_seed"]),
        "panel_seed_used_equals_accepted_for_every_attempt": True,
        "all_rejections_persisted": True,
        "all_attempt_and_panel_seeds_match_plan": True,
        "supervision_path": supervision_path.relative_to(bundle.repo_root).as_posix(),
        "supervision_sha256": sha256_file(supervision_path),
        "all_process_groups_gone_and_children_reaped": True,
        "all_attempt_temp_trees_removed": True,
        "active_supervision_markers": 0,
        "attempt_scratch_directories": 0,
        "unit_atomic_temp_files": 0,
        "pre_marker_startup_cleanups": startup_cleanup_count,
        "startup_cleanup_path": (
            startup_cleanup_path.relative_to(bundle.repo_root).as_posix()
            if startup_cleanup_path.is_file()
            else None
        ),
        "startup_cleanup_sha256": startup_cleanup_sha256,
        "held_simulation_locks": 0,
        "live_supervised_descendant_identities": 0,
        "terminal_attempt_seed_dispositions": len(ledger),
        "interrupted_nonterminal_dispositions": (
            len(pd.read_csv(reconciliation_path, sep="\t"))
            if reconciliation_path.is_file()
            else 0
        ),
        "one_terminal_disposition_per_attempted_seed": True,
        "supervisor_fields_equal_attempt_ledger": True,
        "interrupted_elapsed_seconds_charged": interrupted_charge,
        "committed_elapsed_seconds": committed_elapsed,
        "accounted_cumulative_elapsed_seconds": cumulative_elapsed,
        "authorized_cumulative_timeout_seconds": CUMULATIVE_TIMEOUT_SECONDS,
        "remaining_cumulative_timeout_seconds": (
            CUMULATIVE_TIMEOUT_SECONDS - interrupted_charge
        ),
        "reconciliation_path": (
            reconciliation_path.relative_to(bundle.repo_root).as_posix()
            if reconciliation_path.is_file()
            else None
        ),
        "reconciliation_sha256": reconciliation_sha256,
    }


def audit_stage(
    bundle: DirectV3Bundle, stage: str, *, write: bool = True
) -> dict[str, Any]:
    if stage == STAGE_SIMULATE:
        _required_pilot_audit(bundle)
        audit_label = STAGE_ALL
    elif stage == STAGE_PILOT:
        audit_label = STAGE_PILOT
    else:
        raise ValueError("Han v3 audit stage must be pilot or simulate")
    aggregate_stage = STAGE_PILOT if audit_label == STAGE_PILOT else STAGE_ALL
    records = aggregate_v3_cached_completions(bundle, aggregate_stage).to_dict(
        orient="records"
    )
    expected = 6 if audit_label == STAGE_PILOT else 60
    if len(records) != expected:
        raise ValueError("Han v3 audit accepted-unit count differs")
    expected_replicates = (
        {PILOT_REPLICATE_INDEX}
        if audit_label == STAGE_PILOT
        else set(range(1, REPLICATES_PER_CELL + 1))
    )
    audit_frame = pd.DataFrame(records)
    per_cell_replicates = audit_frame.groupby(
        ["selection_coefficient", "target_allele_frequency"]
    )["replicate_index"].apply(lambda values: set(values.astype(int)))
    if len(per_cell_replicates) != 6 or any(
        values != expected_replicates for values in per_cell_replicates
    ):
        raise ValueError("Han v3 audited per-cell replicate grid differs")
    schedule_audit = _seed_plan_audit(bundle)
    payload: dict[str, Any] = {
        "schema": STAGE_AUDIT_SCHEMA,
        "stage": audit_label,
        "status": "passed",
        "manifest_sha256": sha256_file(bundle.manifest_path),
        "seed_plan_sha256": sha256_file(bundle.seed_plan_path),
        "accepted_units": len(records),
        "expected_units": expected,
        "pilot_units_count_toward_final": True,
        "all_units_use_exact_final_contract": True,
        "all_current_completions_cache_validated_read_only": True,
        "all_current_completion_counts_match_attempt_ledgers": True,
        "all_current_completion_budgets_match_authorization": True,
        "all_current_completion_sources_match_pinned_runtime": True,
        "all_unit_authorizations_match_manifest_and_seed_plan": True,
        "per_cell_replicates": sorted(expected_replicates),
        "seed_plan_audit": schedule_audit,
        "all_units_strict_type1_identity": True,
        "all_rejections_persisted": True,
        "gamma_smc_statistics_used": False,
        "failed_v2_inferential_use": False,
        "rows": records,
    }
    payload["payload_sha256"] = _canonical_sha256(payload)
    if write:
        name = (
            "pilot_audit.json" if audit_label == STAGE_PILOT else "all_units_audit.json"
        )
        _immutable_json(bundle.output_dir / name, payload)
        if audit_label == STAGE_ALL:
            authorization = {
                "schema": AUTHORIZATION_SCHEMA,
                "status": "authorized",
                "manifest_sha256": sha256_file(bundle.manifest_path),
                "seed_plan_sha256": sha256_file(bundle.seed_plan_path),
                "all_units_audit_payload_sha256": payload["payload_sha256"],
                "accepted_units": 60,
                "cells": 6,
                "replicates_per_cell": 10,
                "pilot_units_included": 6,
                "validated_completion_count": 60,
                "validated_attempt_count": int(
                    sum(int(row["attempts_completed"]) for row in records)
                ),
                "attempts_per_unit_budget": ATTEMPTS_PER_UNIT,
                "draw_timeout_seconds": DRAW_TIMEOUT_SECONDS,
                "cumulative_timeout_seconds": CUMULATIVE_TIMEOUT_SECONDS,
                "all_completion_counts_budgets_and_authorizations_validated": True,
                "all_process_groups_gone_and_children_reaped": True,
                "pinned_runtime": bundle.manifest["runtime"],
                "direct_operational_estimand": True,
                "conditional_inference": True,
                "gamma_smc_statistics_used": False,
                "failed_v2_inferential_use": False,
                "scientific_caveats": bundle.manifest["scientific_caveats"],
            }
            authorization["payload_sha256"] = _canonical_sha256(authorization)
            _immutable_json(
                bundle.output_dir / "han_direct_v3_authorization.json",
                authorization,
            )
    return payload


def _directory_digest(path: Path) -> tuple[int, int, str]:
    records = []
    total = 0
    for child in sorted(item for item in path.rglob("*") if item.is_file()):
        relative = child.relative_to(path).as_posix()
        size = child.stat().st_size
        total += size
        records.append(
            {"path": relative, "size_bytes": size, "sha256": sha256_file(child)}
        )
    return len(records), total, _canonical_sha256(records)


def build_quarantine_plan(
    *,
    repo_root: str | Path,
    campaign_dir: str | Path,
    quarantine_id: str = "han_direct_v3_stale_selected",
    expected_stale_count: int = 0,
    manifest_path: str | Path | None = None,
) -> dict[str, Any]:
    """Inventory stale Han selected directories; never move or delete them."""

    root = Path(repo_root).resolve()
    try:
        stale_count = int(expected_stale_count)
        stale_count_float = float(expected_stale_count)
    except (TypeError, ValueError) as error:
        raise ValueError(
            "Han v3 expected stale count must be an exact integer"
        ) from error
    if (
        not math.isfinite(stale_count_float)
        or stale_count_float != stale_count
        or stale_count < 0
    ):
        raise ValueError(
            "Han v3 expected stale count must be a nonnegative exact integer"
        )
    campaign = _inside_repo(campaign_dir, root, label="campaign directory")
    units = _validate_production_units(
        pd.read_csv(campaign / "execution_units.tsv", sep="\t")
    )
    bundle = (
        load_bundle(manifest_path, repo_root=root)
        if manifest_path is not None
        else None
    )
    destination_root = _inside_repo(
        campaign / "calibration/quarantine" / quarantine_id,
        root,
        label="quarantine proposal",
    )
    rows = []
    for unit_id in sorted(units["unit_id"].astype(str)):
        source = campaign / "work" / unit_id
        if not source.is_dir():
            continue
        completion = source / "simulation_complete.json"
        compatible_v3 = False
        if completion.is_file() and bundle is not None:
            try:
                unit = _unit_row(bundle, unit_id)
                compatible_v3 = bool(
                    validate_v3_cached_unit(bundle, unit).get(
                        "cache_validated_read_only"
                    )
                )
            except (OSError, ValueError, TypeError, KeyError):
                compatible_v3 = False
        if compatible_v3:
            continue
        lock_path = source / ".simulation.lock"
        if simulation.unit_lock_is_held(lock_path):
            raise ValueError(f"Han v3 stale unit lock is held: {unit_id}")
        count, size, digest = _directory_digest(source)
        rows.append(
            {
                "unit_id": unit_id,
                "source": source.relative_to(root).as_posix(),
                "proposed_destination": (destination_root / unit_id)
                .relative_to(root)
                .as_posix(),
                "file_count": count,
                "size_bytes": size,
                "payload_sha256": digest,
                "source_exists": True,
                "destination_exists": (destination_root / unit_id).exists(),
            }
        )
    if len(rows) != stale_count:
        raise ValueError(
            f"Han v3 observed {len(rows)} stale Han selected directories, "
            f"expected {stale_count}"
        )
    if any(row["destination_exists"] for row in rows):
        raise ValueError("Han v3 quarantine proposal destination already exists")
    return {
        "schema": QUARANTINE_PLAN_SCHEMA,
        "status": "planning_only",
        "scope": "stale_han_selected_production_directories_only",
        "observed_stale_selected_directories": len(rows),
        "expected_stale_selected_directories": stale_count,
        "proposed_destination_root": destination_root.relative_to(root).as_posix(),
        "moves_performed": 0,
        "deletions_performed": 0,
        "neutral_directories_touched": 0,
        "rows": rows,
    }


def _default_paths(args: argparse.Namespace) -> tuple[Path, Path, Path, Path]:
    root = Path(args.repo_root).resolve()
    campaign = _inside_repo(args.campaign_dir, root, label="campaign directory")
    output = _inside_repo(args.output_dir, root, label="output directory")
    slim = _inside_repo(args.slim_bin, root, label="SLiM executable")
    return root, campaign, output, slim


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Plan, pilot, simulate, and audit Han direct production v3."
    )
    parser.add_argument(
        "command",
        choices=(
            "plan",
            "pilot",
            "pilot-audit",
            "simulate",
            "audit",
            "aggregate",
            "quarantine-plan",
        ),
    )
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument(
        "--campaign-dir", type=Path, default=Path(DEFAULT_CAMPAIGN_RELATIVE_PATH)
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path(DEFAULT_OUTPUT_RELATIVE_PATH)
    )
    parser.add_argument(
        "--slim-bin", type=Path, default=Path(".native-stdpopsim/bin/slim")
    )
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument(
        "--selected-integration-manifest",
        type=Path,
        help=(
            "Optional EAS selected-production integration manifest; when set, "
            "the v3 dispatch is nested with its EAS validator before standard "
            "decoded-output aggregation."
        ),
    )
    parser.add_argument(
        "--local-seed-artifact",
        type=Path,
        action="append",
        default=[],
        help=(
            "Optional local seed evidence; every seed must already be a subset "
            "of the canonical tracked registry."
        ),
    )
    parser.add_argument(
        "--expected-stale-count",
        type=int,
        default=0,
        help="Explicit stale Han selected directory count; defaults fail-closed to 0.",
    )
    parser.add_argument("--quarantine-id", default="han_direct_v3_stale_selected")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    root, campaign, output, slim = _default_paths(args)
    if args.command == "plan":
        manifest, seed_plan, excluded_registry, excluded_inventory = build_plan(
            repo_root=root,
            campaign_dir=campaign,
            slim_path=slim,
            local_seed_artifacts=args.local_seed_artifact,
        )
        payload = publish_plan_bundle(
            output,
            repo_root=root,
            manifest=manifest,
            seed_plan=seed_plan,
            excluded_seed_registry=excluded_registry,
            excluded_seed_inventory=excluded_inventory,
        )
        print(
            json.dumps(
                {
                    "status": "planned",
                    "manifest": str(output / DEFAULT_PLAN_MANIFEST_NAME),
                    "manifest_sha256": sha256_file(output / DEFAULT_PLAN_MANIFEST_NAME),
                    "units": payload["design"]["units"],
                    "attempt_seed_rows": payload["seed_registry"]["attempt_seed_count"],
                    "excluded_seed_rows": payload["excluded_seed_registry_artifact"][
                        "rows"
                    ],
                    "pilot_units": payload["design"]["pilot_units"],
                    "remaining_units": payload["design"]["remaining_simulation_units"],
                    "gamma_smc_statistics_used": False,
                },
                indent=2,
                sort_keys=True,
            )
        )
        return 0
    if args.command == "quarantine-plan":
        manifest_path = output / DEFAULT_PLAN_MANIFEST_NAME
        payload = build_quarantine_plan(
            repo_root=root,
            campaign_dir=campaign,
            quarantine_id=args.quarantine_id,
            expected_stale_count=args.expected_stale_count,
            manifest_path=manifest_path if manifest_path.is_file() else None,
        )
        payload["payload_sha256"] = _canonical_sha256(payload)
        _immutable_json(output.parent / "han_direct_v3_quarantine_plan.json", payload)
        print(json.dumps(payload, indent=2, sort_keys=True))
        return 0
    bundle = load_bundle(
        output / DEFAULT_PLAN_MANIFEST_NAME,
        repo_root=root,
        slim_path=slim,
    )
    if args.command == "pilot-audit":
        payload = audit_stage(bundle, STAGE_PILOT)
        print(json.dumps(payload, indent=2, sort_keys=True))
        return 0
    if args.command == "audit":
        payload = audit_stage(bundle, STAGE_SIMULATE)
        print(json.dumps(payload, indent=2, sort_keys=True))
        return 0
    if args.command == "aggregate":
        units = pd.read_csv(bundle.campaign_dir / "execution_units.tsv", sep="\t")
        outputs = aggregate_decoded_outputs_with_v3_dispatch(
            bundle,
            units,
            selected_integration_manifest=args.selected_integration_manifest,
        )
        print(
            json.dumps(
                {
                    "status": "aggregated",
                    "all_execution_units": len(units),
                    "han_v3_units_validated": 60,
                    "standard_downstream_inputs": [
                        path.relative_to(root).as_posix() for path in outputs
                    ],
                },
                indent=2,
                sort_keys=True,
            )
        )
        return 0
    stage = STAGE_PILOT if args.command == "pilot" else STAGE_SIMULATE
    results = run_stage(bundle, stage, workers=args.workers)
    failed = results[results["status"].astype(str).eq("failed")]
    response: dict[str, Any] = {
        "stage": stage,
        "status": "failed" if len(failed) else "complete",
        "planned_units": len(results),
        "completed_or_cached_units": len(results) - len(failed),
        "failed_units": len(failed),
        "fail_closed_on_any_unit_exhaustion": True,
        "run_snapshot": str(output / f"{stage}_run_snapshot.tsv"),
        "run_snapshot_sha256": sha256_file(output / f"{stage}_run_snapshot.tsv"),
    }
    if not len(failed):
        response["audit"] = audit_stage(bundle, stage)
    print(json.dumps(response, indent=2, sort_keys=True))
    return 2 if len(failed) else 0


if __name__ == "__main__":
    raise SystemExit(main())
