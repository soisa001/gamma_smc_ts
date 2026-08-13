"""Additive, provenance-preserving recovery for the Han v3 dtype failure.

The published Han direct-production v3 plan and implementation are immutable.
This module therefore does not replace either one.  It strict-loads the
published bundle, delegates all seed, timeout, reconciliation, and simulation
semantics to the frozen v3 executor, and changes only the in-memory dtype of
the ``rejection_reason`` column before the frozen attempt-ledger callback runs.

``dry-plan`` is deliberately read-only.  ``resume`` requires the exact digest
printed by a current dry plan, snapshots every incomplete unit byte-for-byte
into checksum-bound history, publishes an additive recovery authorization,
and only then starts worker-local recovery contexts.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import shutil
from collections.abc import Iterator, Mapping, Sequence
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter, time_ns
from typing import Any

import numpy as np
import pandas as pd

from . import focused_han_direct_v3 as v3
from . import focused_selection_simulation as simulation

RECOVERY_PLAN_SCHEMA = "gamma-smc.han-direct-production-v3-dtype-recovery-plan/v1"
RECOVERY_HISTORY_SCHEMA = (
    "gamma-smc.han-direct-production-v3-dtype-recovery-history/v1"
)
RECOVERY_AUTHORIZATION_SCHEMA = (
    "gamma-smc.han-direct-production-v3-dtype-recovery-authorization/v1"
)
RECOVERY_COMPLETION_AUDIT_SCHEMA = (
    "gamma-smc.han-direct-production-v3-dtype-recovery-completion-audit/v1"
)

MODULE_SOURCE_PATH = "python/gamma_smc_aou/focused_han_direct_v3_recovery.py"
WRAPPER_SOURCE_PATH = "scripts/run_focused_han_direct_v3_recovery.py"
TEST_SOURCE_PATH = "tests/test_focused_han_direct_v3_recovery.py"
RECOVERY_SOURCE_PATHS = (MODULE_SOURCE_PATH, WRAPPER_SOURCE_PATH)
CANONICAL_AUTHORIZATION_NAME = "han_direct_v3_dtype_recovery_authorization.json"
CANONICAL_AUDIT_NAME = "han_direct_v3_dtype_recovery_audit.json"

KNOWN_DTYPE_ERROR_TYPE = "TypeError"
KNOWN_DTYPE_ERROR_PATTERN = re.compile(
    r"Invalid value ['\"]draw_timeout['\"] for dtype ['\"]float64['\"]"
)
EXPECTED_STOPPED_COMPLETED_SIMULATE_UNITS = 7
EXPECTED_STOPPED_ACTIVE_UNITS = 25
EXPECTED_STOPPED_UNSTARTED_UNITS = 22
EXPECTED_STOPPED_DTYPE_FAILURE_UNITS = 1
EXPECTED_STOPPED_RECOVERABLE_UNITS = 47
_FROZEN_PATCHED_EXECUTOR = v3.patched_v3_executor
_REAL_PANDAS = pd
_WORKER_RECOVERY: RecoveryBundle | None = None


@dataclass(frozen=True)
class RecoveryBundle:
    direct: v3.DirectV3Bundle
    recovery_sources: Mapping[str, str]

    @property
    def repo_root(self) -> Path:
        return self.direct.repo_root

    @property
    def campaign_dir(self) -> Path:
        return self.direct.campaign_dir

    @property
    def output_dir(self) -> Path:
        return self.direct.output_dir


class _RecoveryPandasProxy:
    """Delegate pandas while dtype-fixing only reads of attempts.tsv."""

    def __getattr__(self, name: str) -> Any:
        return getattr(_REAL_PANDAS, name)

    @staticmethod
    def read_csv(*args: Any, **kwargs: Any) -> pd.DataFrame:
        frame = _REAL_PANDAS.read_csv(*args, **kwargs)
        source = args[0] if args else kwargs.get("filepath_or_buffer")
        try:
            name = Path(source).name
        except TypeError:
            name = ""
        if name == "attempts.tsv" and "rejection_reason" in frame:
            frame = _coerce_rejection_reason_object(frame)
        return frame


def _canonical_sha256(value: Any) -> str:
    return v3._canonical_sha256(value)


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


def _immutable_json(path: Path, payload: Any) -> None:
    if path.is_file():
        try:
            current = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            raise ValueError(f"recovery JSON is unreadable: {path}") from error
        if current != payload:
            raise ValueError(f"immutable recovery JSON differs: {path}")
        return
    _atomic_json(path, payload)


def _recovery_source_records(repo_root: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    for relative in RECOVERY_SOURCE_PATHS:
        path = repo_root / relative
        if not path.is_file():
            raise ValueError(f"Han v3 recovery source is absent: {relative}")
        records[relative] = v3.sha256_file(path)
    return records


def load_recovery_bundle(
    manifest_path: str | Path,
    *,
    repo_root: str | Path,
    slim_path: str | Path | None = None,
) -> RecoveryBundle:
    """Strict-load the original bundle before authorizing the additive seam."""

    direct = v3.load_bundle(
        manifest_path,
        repo_root=repo_root,
        slim_path=slim_path,
    )
    bound_sources = direct.manifest.get("implementation_sources", {})
    current_v3_sha = v3.sha256_file(direct.repo_root / v3.MODULE_SOURCE_PATH)
    if bound_sources.get(v3.MODULE_SOURCE_PATH) != current_v3_sha:
        raise ValueError("Han v3 recovery original source binding differs")
    current_simulation_sha = v3.sha256_file(
        direct.repo_root / simulation.SIMULATION_SOURCE_PATH
    )
    if bound_sources.get(simulation.SIMULATION_SOURCE_PATH) != current_simulation_sha:
        raise ValueError("Han v3 recovery simulation source binding differs")
    return RecoveryBundle(
        direct=direct,
        recovery_sources=_recovery_source_records(direct.repo_root),
    )


def _coerce_rejection_reason_object(frame: pd.DataFrame) -> pd.DataFrame:
    """Return a copy whose sole semantic change is one object-capable column."""

    normalized = frame.copy()
    if "rejection_reason" not in normalized:
        normalized["rejection_reason"] = pd.Series(
            [""] * len(normalized), index=normalized.index, dtype="object"
        )
    else:
        normalized["rejection_reason"] = normalized["rejection_reason"].astype(
            "object"
        )
    return normalized


@contextmanager
def patched_recovery_executor(
    bundle: v3.DirectV3Bundle,
    unit: Mapping[str, Any],
) -> Iterator[Any]:
    """Apply only the dtype seam around the frozen process-local executor."""

    prior_v3_pandas = v3.pd
    v3.pd = _RecoveryPandasProxy()
    try:
        with _FROZEN_PATCHED_EXECUTOR(bundle, unit) as supervised_engine:
            frozen_atomic_frame = simulation._atomic_frame

            def dtype_safe_atomic_frame(path: Path, frame: pd.DataFrame) -> None:
                if Path(path).name == "attempts.tsv" and (
                    "attempt_zero_based" in frame
                ):
                    frame = _coerce_rejection_reason_object(frame)
                frozen_atomic_frame(path, frame)

            simulation._atomic_frame = dtype_safe_atomic_frame
            try:
                yield supervised_engine
            finally:
                simulation._atomic_frame = frozen_atomic_frame
    finally:
        v3.pd = prior_v3_pandas


@contextmanager
def _install_recovery_executor() -> Iterator[None]:
    prior = v3.patched_v3_executor
    v3.patched_v3_executor = patched_recovery_executor
    try:
        yield
    finally:
        v3.patched_v3_executor = prior


def _unit_inventory(unit_dir: Path) -> tuple[list[dict[str, Any]], int]:
    if unit_dir.is_symlink():
        raise ValueError("Han v3 recovery unit directory must not be a symlink")
    rows: list[dict[str, Any]] = []
    total = 0
    if not unit_dir.exists():
        return rows, total
    for path in sorted(unit_dir.rglob("*"), key=lambda item: item.as_posix()):
        if path.is_symlink():
            raise ValueError("Han v3 recovery evidence must not contain symlinks")
        if path.is_dir():
            continue
        if not path.is_file():
            raise ValueError("Han v3 recovery evidence contains a special file")
        relative = path.relative_to(unit_dir).as_posix()
        size = path.stat().st_size
        total += size
        rows.append(
            {
                "path": relative,
                "size_bytes": size,
                "sha256": v3.sha256_file(path),
            }
        )
    return rows, total


def _read_optional_frame(path: Path) -> pd.DataFrame:
    if not path.is_file():
        return pd.DataFrame()
    frame = pd.read_csv(path, sep="\t")
    if path.name == "attempts.tsv" and "rejection_reason" in frame:
        frame = _coerce_rejection_reason_object(frame)
    return frame


def _exact_attempts(
    bundle: v3.DirectV3Bundle,
    unit: Mapping[str, Any],
    unit_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    schedule = v3._unit_seed_schedule(bundle, str(unit["unit_id"]))
    attempts = _read_optional_frame(unit_dir / "attempts.tsv")
    if not attempts.empty:
        observed = list(
            v3._exact_integer_series(
                attempts,
                "attempt_zero_based",
                label="recovery base attempt ledger",
                minimum=0,
            )
        )
        if observed != list(range(len(observed))):
            raise ValueError("Han v3 recovery attempt ledger is noncontiguous")
        if len(attempts) >= v3.ATTEMPTS_PER_UNIT:
            raise ValueError("Han v3 recovery incomplete unit exhausted its attempts")
        seeds = v3._exact_integer_series(
            attempts,
            "seed",
            label="recovery base attempt ledger",
            minimum=1,
        )
        expected = schedule.iloc[: len(attempts)]["simulation_seed"].astype(int)
        if list(seeds) != list(expected):
            raise ValueError("Han v3 recovery committed attempt seeds differ")
        elapsed = pd.to_numeric(attempts.get("elapsed_seconds"), errors="coerce")
        if elapsed.isna().any() or (elapsed < 0).any():
            raise ValueError("Han v3 recovery committed elapsed time differs")
        if float(elapsed.sum()) > v3.CUMULATIVE_TIMEOUT_SECONDS + 1e-6:
            raise ValueError("Han v3 recovery committed cumulative budget differs")
    supervision = _read_optional_frame(unit_dir / "supervision.tsv")
    if not supervision.empty:
        observed = list(
            v3._exact_integer_series(
                supervision,
                "attempt_zero_based",
                label="recovery supervision ledger",
                minimum=0,
            )
        )
        if observed != list(range(len(observed))):
            raise ValueError("Han v3 recovery supervision ledger is noncontiguous")
        seeds = v3._exact_integer_series(
            supervision,
            "simulation_seed",
            label="recovery supervision ledger",
            minimum=1,
        )
        expected = schedule.iloc[: len(supervision)]["simulation_seed"].astype(int)
        if list(seeds) != list(expected):
            raise ValueError("Han v3 recovery supervision seeds differ")
        if len(supervision) < len(attempts) or len(supervision) > len(attempts) + 1:
            raise ValueError("Han v3 recovery supervision/attempt counts differ")
    return attempts, supervision


def _strict_partial_contract(
    bundle: v3.DirectV3Bundle,
    unit: Mapping[str, Any],
    unit_dir: Path,
) -> tuple[dict[str, Any], str]:
    path = unit_dir / "simulation_contract.json"
    try:
        contract = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError("Han v3 recovery partial contract is unreadable") from error
    record = simulation._validate_unit(unit)
    if _canonical_sha256(contract.get("unit")) != _canonical_sha256(record):
        raise ValueError("Han v3 recovery partial contract unit differs")
    expected_sources = {
        relative: v3.sha256_file(bundle.repo_root / relative)
        for relative in v3._v3_implementation_paths(bundle)
    }
    implementation = contract.get("implementation", {})
    if implementation.get("sources") != expected_sources:
        raise ValueError("Han v3 recovery partial implementation binding differs")
    binding = v3._unit_binding(bundle, record)
    parameters = contract.get("parameters", {})
    if _canonical_sha256(parameters.get("han_selection_calibration")) != (
        _canonical_sha256(binding)
    ):
        raise ValueError("Han v3 recovery partial authorization differs")
    return contract, _canonical_sha256(contract)


def _known_dtype_failure(
    path: Path,
    *,
    unit_id: str,
    contract_sha256: str,
) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError("Han v3 recovery failure record is unreadable") from error
    if (
        payload.get("schema") != simulation.SCHEMA_VERSION
        or payload.get("status") != "failed"
        or payload.get("unit_id") != unit_id
        or payload.get("contract_sha256") != contract_sha256
        or payload.get("error_type") != KNOWN_DTYPE_ERROR_TYPE
        or KNOWN_DTYPE_ERROR_PATTERN.fullmatch(str(payload.get("error", ""))) is None
    ):
        raise ValueError("Han v3 recovery encountered a non-dtype failure")
    budget = payload.get("search_budget", {})
    if (
        int(budget.get("selected_max_draws", -1)) != v3.ATTEMPTS_PER_UNIT
        or not math.isclose(
            float(budget.get("selected_draw_timeout_seconds", math.nan)),
            v3.DRAW_TIMEOUT_SECONDS,
            abs_tol=1e-12,
        )
        or not math.isclose(
            float(budget.get("selected_cumulative_timeout_seconds", math.nan)),
            v3.CUMULATIVE_TIMEOUT_SECONDS,
            abs_tol=1e-12,
        )
    ):
        raise ValueError("Han v3 recovery failure budget differs")
    return payload


def _active_marker_record(
    bundle: v3.DirectV3Bundle,
    unit: Mapping[str, Any],
    unit_dir: Path,
    attempts: pd.DataFrame,
    supervision: pd.DataFrame,
) -> tuple[Path, dict[str, Any]] | None:
    active_paths = sorted(unit_dir.glob(".han_v3_active_*.json"))
    if len(active_paths) > 1:
        raise ValueError("Han v3 recovery unit has multiple active markers")
    if not active_paths:
        return None
    active_path = active_paths[0]
    match = re.fullmatch(r"\.han_v3_active_([0-9]{3})\.json", active_path.name)
    if match is None:
        raise ValueError("Han v3 recovery active marker name differs")
    payload = v3._read_active_supervision(active_path)
    record = dict(payload["attempt_record"])
    attempt = int(record["attempt_zero_based"])
    if attempt != int(match.group(1)):
        raise ValueError("Han v3 recovery active attempt/name differs")
    schedule = v3._unit_seed_schedule(bundle, str(unit["unit_id"]))
    if attempt not in schedule.index or int(record["simulation_seed"]) != int(
        schedule.loc[attempt, "simulation_seed"]
    ):
        raise ValueError("Han v3 recovery active attempt seed differs")
    if attempt not in {len(attempts), max(0, len(attempts) - 1)}:
        raise ValueError("Han v3 recovery active/committed attempt position differs")
    if not supervision.empty and attempt < len(supervision):
        persisted = supervision.iloc[attempt]
        if (
            int(persisted["simulation_seed"]) != int(record["simulation_seed"])
            or str(persisted.get("execution_id", ""))
            != v3._SupervisedSlimEngine._execution_id(record)
        ):
            raise ValueError("Han v3 recovery active/supervision identity differs")

    if v3._same_linux_process(record.get("supervisor_process_identity")):
        raise ValueError("Han v3 recovery supervisor process is still alive")
    for label, pid_field, identity_field in (
        ("child", "child_pid", "child_process_identity"),
        ("engine wrapper", "engine_wrapper_pid", "engine_wrapper_process_identity"),
    ):
        raw_pid = record.get(pid_field)
        if raw_pid is None or pd.isna(raw_pid):
            continue
        if v3._same_linux_process(record.get(identity_field)) and not (
            v3._process_absent_or_zombie(int(raw_pid))
        ):
            raise ValueError(f"Han v3 recovery {label} process is still alive")
    raw_pgid = record.get("child_pgid")
    if raw_pgid is not None and not pd.isna(raw_pgid) and (
        v3._process_group_has_live_members(int(raw_pgid))
    ):
        raise ValueError("Han v3 recovery child process group is still alive")
    return active_path, record


def _inspect_incomplete_unit(
    bundle: v3.DirectV3Bundle,
    unit: Mapping[str, Any],
) -> dict[str, Any]:
    unit_id = str(unit["unit_id"])
    unit_dir = bundle.campaign_dir / "work" / unit_id
    if not unit_dir.exists():
        return {
            "unit_id": unit_id,
            "state": "unstarted",
            "replicate_index": int(unit["replicate_index"]),
            "selection_coefficient": float(unit["selection_coefficient"]),
            "target_allele_frequency": float(unit["target_allele_frequency"]),
            "committed_attempts": 0,
            "supervision_records": 0,
            "active_attempt_zero_based": None,
            "active_simulation_seed": None,
            "active_execution_id": None,
            "active_marker_path": None,
            "active_marker_sha256": None,
            "failure_record_sha256": None,
            "partial_contract_sha256": None,
            "partial_contract_matches_original_frozen_sources": None,
            "all_recorded_processes_and_groups_absent": True,
            "artifact_file_count": 0,
            "artifact_total_bytes": 0,
            "artifact_inventory_sha256": _canonical_sha256([]),
        }
    if not unit_dir.is_dir() or unit_dir.is_symlink():
        raise ValueError("Han v3 recovery unit path is not a real directory")
    if (unit_dir / "simulation_complete.json").exists():
        raise ValueError("Han v3 recovery incomplete inspector received a completion")
    contract, contract_sha256 = _strict_partial_contract(bundle, unit, unit_dir)
    attempts, supervision = _exact_attempts(bundle, unit, unit_dir)
    active = _active_marker_record(bundle, unit, unit_dir, attempts, supervision)
    failure = _known_dtype_failure(
        unit_dir / "simulation_failed.json",
        unit_id=unit_id,
        contract_sha256=contract_sha256,
    )
    quiescent_between_attempts = False
    if active is None and failure is None:
        accepted = (
            pd.Series(False, index=attempts.index)
            if "accepted" not in attempts
            else attempts["accepted"].map(
                lambda value: (
                    bool(value)
                    if isinstance(value, (bool, np.bool_))
                    else str(value).casefold() == "true"
                )
            )
        )
        forbidden_outputs = (
            "simulation.trees",
            "sample_manifest.tsv",
            "truth_profiles.tsv.gz",
            "truth_class_summaries.tsv",
            "simulation_complete.json",
        )
        if (
            accepted.any()
            or len(supervision) != len(attempts)
            or list(unit_dir.glob(".han_v3_attempt_*"))
            or v3._unit_atomic_temp_paths(unit_dir)
            or any((unit_dir / name).exists() for name in forbidden_outputs)
            or (unit_dir / "pairs").exists()
        ):
            raise ValueError(
                "Han v3 recovery incomplete unit has no recoverable disposition"
            )
        quiescent_between_attempts = True
    inventory, total = _unit_inventory(unit_dir)
    active_path, active_record = active if active is not None else (None, None)
    return {
        "unit_id": unit_id,
        "state": (
            "known_dtype_failure_with_active_marker"
            if failure is not None and active is not None
            else "known_dtype_failure"
            if failure is not None
            else "quiescent_between_attempts"
            if quiescent_between_attempts
            else "interrupted_active_attempt"
        ),
        "replicate_index": int(unit["replicate_index"]),
        "selection_coefficient": float(unit["selection_coefficient"]),
        "target_allele_frequency": float(unit["target_allele_frequency"]),
        "committed_attempts": len(attempts),
        "supervision_records": len(supervision),
        "active_attempt_zero_based": (
            None if active_record is None else int(active_record["attempt_zero_based"])
        ),
        "active_simulation_seed": (
            None if active_record is None else int(active_record["simulation_seed"])
        ),
        "active_execution_id": (
            None
            if active_record is None
            else v3._SupervisedSlimEngine._execution_id(active_record)
        ),
        "active_marker_path": None if active_path is None else active_path.name,
        "active_marker_sha256": (
            None if active_path is None else v3.sha256_file(active_path)
        ),
        "failure_record_sha256": (
            None
            if failure is None
            else v3.sha256_file(unit_dir / "simulation_failed.json")
        ),
        "partial_contract_sha256": contract_sha256,
        "partial_contract_matches_original_frozen_sources": bool(contract),
        "all_recorded_processes_and_groups_absent": True,
        "artifact_file_count": len(inventory),
        "artifact_total_bytes": total,
        "artifact_inventory_sha256": _canonical_sha256(inventory),
    }


def build_recovery_plan(bundle: RecoveryBundle) -> dict[str, Any]:
    """Build a deterministic, read-only plan over the exact stopped state."""

    v3._required_pilot_audit(bundle.direct)
    stage_units = v3._stage_units(bundle.direct, v3.STAGE_SIMULATE)
    completed: list[dict[str, Any]] = []
    recover: list[dict[str, Any]] = []
    for unit in stage_units.to_dict(orient="records"):
        unit_dir = bundle.campaign_dir / "work" / str(unit["unit_id"])
        if (unit_dir / "simulation_complete.json").is_file():
            validated = v3.validate_v3_cached_unit(bundle.direct, unit)
            completed.append(
                {
                    "unit_id": str(unit["unit_id"]),
                    "completion_sha256": validated["completion_sha256"],
                    "attempts_completed": int(validated["attempts_completed"]),
                }
            )
        else:
            recover.append(_inspect_incomplete_unit(bundle.direct, unit))
    if len(completed) + len(recover) != 54 or not recover:
        raise ValueError("Han v3 recovery stage partition differs")
    active_count = sum(row["active_marker_path"] is not None for row in recover)
    dtype_failure_count = sum(
        row["failure_record_sha256"] is not None for row in recover
    )
    unstarted_count = sum(row["state"] == "unstarted" for row in recover)
    expected_stopped_counts = {
        "completed": EXPECTED_STOPPED_COMPLETED_SIMULATE_UNITS,
        "recoverable": EXPECTED_STOPPED_RECOVERABLE_UNITS,
        "active": EXPECTED_STOPPED_ACTIVE_UNITS,
        "unstarted": EXPECTED_STOPPED_UNSTARTED_UNITS,
        "dtype_failure": EXPECTED_STOPPED_DTYPE_FAILURE_UNITS,
    }
    observed_stopped_counts = {
        "completed": len(completed),
        "recoverable": len(recover),
        "active": active_count,
        "unstarted": unstarted_count,
        "dtype_failure": dtype_failure_count,
    }
    if observed_stopped_counts != expected_stopped_counts:
        raise ValueError(
            "Han v3 recovery stopped-state counts differ: "
            f"observed={observed_stopped_counts!r}, "
            f"expected={expected_stopped_counts!r}"
        )
    payload: dict[str, Any] = {
        "schema": RECOVERY_PLAN_SCHEMA,
        "status": "ready_read_only",
        "manifest_path": bundle.direct.manifest_path.relative_to(
            bundle.repo_root
        ).as_posix(),
        "manifest_sha256": v3.sha256_file(bundle.direct.manifest_path),
        "seed_plan_sha256": v3.sha256_file(bundle.direct.seed_plan_path),
        "original_implementation_sources": dict(
            bundle.direct.manifest["implementation_sources"]
        ),
        "recovery_sources": dict(bundle.recovery_sources),
        "stage": v3.STAGE_SIMULATE,
        "expected_stage_units": 54,
        "validated_completed_units": len(completed),
        "recoverable_units": len(recover),
        "active_markers": active_count,
        "unstarted_units": unstarted_count,
        "known_dtype_failure_units": dtype_failure_count,
        "expected_stopped_counts": expected_stopped_counts,
        "completed": sorted(completed, key=lambda row: row["unit_id"]),
        "recover": sorted(recover, key=lambda row: row["unit_id"]),
        "mutation_scope": "process_local_rejection_reason_frame_dtype_only",
        "original_plan_manifest_seed_schedule_and_scientific_contract_unchanged": True,
        "original_attempt_supervision_and_active_artifacts_unchanged_by_dry_plan": True,
        "every_active_attempt_uses_its_original_planned_seed": True,
        "frozen_reconciliation_owns_elapsed_time_and_restart_disposition": True,
        "moves_performed": 0,
        "deletions_performed": 0,
        "writes_performed": 0,
    }
    payload["payload_sha256"] = _canonical_sha256(payload)
    return payload


def _validate_plan_payload(payload: Mapping[str, Any]) -> None:
    unhashed = dict(payload)
    digest = unhashed.pop("payload_sha256", None)
    if (
        payload.get("schema") != RECOVERY_PLAN_SCHEMA
        or payload.get("status") != "ready_read_only"
        or digest != _canonical_sha256(unhashed)
    ):
        raise ValueError("Han v3 recovery dry plan is corrupt")


def snapshot_recovery_history(
    bundle: RecoveryBundle,
    plan: Mapping[str, Any],
) -> tuple[Path, dict[str, Any]]:
    """Copy stopped-state evidence before the frozen reconciler may change it."""

    _validate_plan_payload(plan)
    plan_sha256 = str(plan["payload_sha256"])
    history_root = bundle.output_dir / "recovery_history"
    destination = history_root / plan_sha256
    manifest_path = destination / "history_manifest.json"
    if destination.exists():
        try:
            payload = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            raise ValueError("Han v3 recovery history is unreadable") from error
        unhashed = dict(payload)
        digest = unhashed.pop("payload_sha256", None)
        if (
            payload.get("schema") != RECOVERY_HISTORY_SCHEMA
            or payload.get("dry_plan_payload_sha256") != plan_sha256
            or digest != _canonical_sha256(unhashed)
        ):
            raise ValueError("Han v3 recovery history differs")
        for record in payload["files"]:
            snapshot = destination / str(record["snapshot_path"])
            if (
                not snapshot.is_file()
                or snapshot.stat().st_size != int(record["size_bytes"])
                or v3.sha256_file(snapshot) != record["sha256"]
            ):
                raise ValueError("Han v3 recovery history payload drifted")
        return manifest_path, payload

    history_root.mkdir(parents=True, exist_ok=True)
    temporary = history_root / f".{plan_sha256}.tmp.{os.getpid()}"
    if temporary.exists():
        raise ValueError("Han v3 recovery temporary history already exists")
    temporary.mkdir()
    copied: list[dict[str, Any]] = []
    try:
        for unit_record in plan["recover"]:
            unit_id = str(unit_record["unit_id"])
            unit_dir = bundle.campaign_dir / "work" / unit_id
            inventory, total = _unit_inventory(unit_dir)
            if (
                len(inventory) != int(unit_record["artifact_file_count"])
                or total != int(unit_record["artifact_total_bytes"])
                or _canonical_sha256(inventory)
                != unit_record["artifact_inventory_sha256"]
            ):
                raise ValueError("Han v3 recovery evidence changed after dry plan")
            for record in inventory:
                source = unit_dir / str(record["path"])
                snapshot_relative = Path("units") / unit_id / str(record["path"])
                snapshot = temporary / snapshot_relative
                snapshot.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, snapshot)
                if v3.sha256_file(snapshot) != record["sha256"]:
                    raise ValueError("Han v3 recovery snapshot checksum differs")
                copied.append(
                    {
                        "source_path": source.relative_to(
                            bundle.repo_root
                        ).as_posix(),
                        "snapshot_path": snapshot_relative.as_posix(),
                        **record,
                    }
                )
        history: dict[str, Any] = {
            "schema": RECOVERY_HISTORY_SCHEMA,
            "status": "immutable_stopped_state_snapshot",
            "dry_plan_payload_sha256": plan_sha256,
            "manifest_sha256": plan["manifest_sha256"],
            "seed_plan_sha256": plan["seed_plan_sha256"],
            "original_implementation_sources": plan[
                "original_implementation_sources"
            ],
            "recovery_sources": plan["recovery_sources"],
            "recoverable_units": int(plan["recoverable_units"]),
            "file_count": len(copied),
            "total_bytes": sum(int(row["size_bytes"]) for row in copied),
            "files": copied,
            "original_evidence_moved": False,
            "original_evidence_deleted": False,
        }
        history["payload_sha256"] = _canonical_sha256(history)
        (temporary / "history_manifest.json").write_text(
            json.dumps(history, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, destination)
        return manifest_path, history
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)


def publish_recovery_authorization(
    bundle: RecoveryBundle,
    plan: Mapping[str, Any],
    history_manifest_path: Path,
    history: Mapping[str, Any],
) -> Path:
    _validate_plan_payload(plan)
    relative_history = history_manifest_path.relative_to(bundle.repo_root).as_posix()
    payload: dict[str, Any] = {
        "schema": RECOVERY_AUTHORIZATION_SCHEMA,
        "status": "authorized_after_immutable_history_snapshot",
        "dry_plan_payload_sha256": plan["payload_sha256"],
        "manifest_sha256": plan["manifest_sha256"],
        "seed_plan_sha256": plan["seed_plan_sha256"],
        "history_manifest_path": relative_history,
        "history_manifest_sha256": v3.sha256_file(history_manifest_path),
        "history_payload_sha256": history["payload_sha256"],
        "original_implementation_sources": plan["original_implementation_sources"],
        "recovery_sources": plan["recovery_sources"],
        "recoverable_unit_ids": [row["unit_id"] for row in plan["recover"]],
        "process_local_patch": {
            "input": "attempts.tsv frame passed to frozen v3 callback",
            "column": "rejection_reason",
            "change": "cast existing or new column to object dtype",
            "planned_seed_change": False,
            "retry_change": False,
            "timeout_change": False,
            "budget_change": False,
            "reconciliation_change": False,
            "scientific_model_change": False,
        },
    }
    payload["payload_sha256"] = _canonical_sha256(payload)
    path = (
        bundle.output_dir
        / "recovery_authorizations"
        / f"recovery_{plan['payload_sha256']}.json"
    )
    _immutable_json(path, payload)
    _immutable_json(bundle.output_dir / CANONICAL_AUTHORIZATION_NAME, payload)
    return path


def _load_authorization(
    bundle: RecoveryBundle,
    path: str | Path,
) -> dict[str, Any]:
    authorization_path = Path(path).resolve()
    try:
        authorization_path.relative_to(bundle.repo_root)
    except ValueError as error:
        raise ValueError("Han v3 recovery authorization is outside repo") from error
    try:
        payload = json.loads(authorization_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError("Han v3 recovery authorization is unreadable") from error
    unhashed = dict(payload)
    digest = unhashed.pop("payload_sha256", None)
    if (
        payload.get("schema") != RECOVERY_AUTHORIZATION_SCHEMA
        or payload.get("status") != "authorized_after_immutable_history_snapshot"
        or payload.get("manifest_sha256")
        != v3.sha256_file(bundle.direct.manifest_path)
        or payload.get("seed_plan_sha256")
        != v3.sha256_file(bundle.direct.seed_plan_path)
        or payload.get("original_implementation_sources")
        != bundle.direct.manifest["implementation_sources"]
        or payload.get("recovery_sources") != dict(bundle.recovery_sources)
        or digest != _canonical_sha256(unhashed)
    ):
        raise ValueError("Han v3 recovery authorization differs")
    history_path = bundle.repo_root / str(payload["history_manifest_path"])
    if (
        not history_path.is_file()
        or v3.sha256_file(history_path) != payload["history_manifest_sha256"]
    ):
        raise ValueError("Han v3 recovery history binding differs")
    return payload


def _simulate_one_recovery(
    bundle: RecoveryBundle,
    unit: Mapping[str, Any],
    authorization_path: Path,
) -> dict[str, Any]:
    started = perf_counter()
    try:
        authorization = _load_authorization(bundle, authorization_path)
        if str(unit["unit_id"]) not in authorization["recoverable_unit_ids"]:
            raise ValueError("Han v3 recovery unit is not authorized")
        with _install_recovery_executor():
            result = v3._simulate_one(bundle.direct, unit)
        result["recovery_authorization_sha256"] = v3.sha256_file(
            authorization_path
        )
        return result
    except Exception as error:  # noqa: BLE001 - worker result must fail closed
        return {
            "unit_id": str(unit["unit_id"]),
            "status": "failed",
            "completion_path": "",
            "elapsed_seconds": perf_counter() - started,
            "error_type": type(error).__name__,
            "error": str(error),
            "recovery_authorization_sha256": "",
        }


def _worker_initialize(
    manifest_path: str,
    repo_root: str,
    slim_path: str,
    authorization_path: str,
) -> None:
    global _WORKER_RECOVERY
    _WORKER_RECOVERY = load_recovery_bundle(
        manifest_path,
        repo_root=repo_root,
        slim_path=slim_path,
    )
    _load_authorization(_WORKER_RECOVERY, authorization_path)


def _worker_simulate(
    unit: Mapping[str, Any], authorization_path: str
) -> dict[str, Any]:
    if _WORKER_RECOVERY is None:
        raise RuntimeError("Han v3 recovery worker is not initialized")
    return _simulate_one_recovery(
        _WORKER_RECOVERY,
        unit,
        Path(authorization_path),
    )


def _authorized_current_state(
    bundle: RecoveryBundle,
    authorization: Mapping[str, Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Validate caches and inspect only the immutable authorized recovery set."""

    authorized_ids = sorted(
        str(value) for value in authorization["recoverable_unit_ids"]
    )
    stage_units = v3._stage_units(bundle.direct, v3.STAGE_SIMULATE)
    stage_by_id = {
        str(row["unit_id"]): row for row in stage_units.to_dict(orient="records")
    }
    if len(authorized_ids) != EXPECTED_STOPPED_RECOVERABLE_UNITS or any(
        unit_id not in stage_by_id for unit_id in authorized_ids
    ):
        raise ValueError("Han v3 recovery authorized unit grid differs")
    completed: list[dict[str, Any]] = []
    incomplete: list[dict[str, Any]] = []
    for unit_id in authorized_ids:
        unit = stage_by_id[unit_id]
        completion_path = (
            bundle.campaign_dir / "work" / unit_id / "simulation_complete.json"
        )
        if completion_path.is_file():
            validated = v3.validate_v3_cached_unit(bundle.direct, unit)
            completed.append(
                {
                    "unit_id": validated["unit_id"],
                    "completion_sha256": validated["completion_sha256"],
                    "attempts_completed": int(validated["attempts_completed"]),
                }
            )
        else:
            incomplete.append(_inspect_incomplete_unit(bundle.direct, unit))
    if len(completed) + len(incomplete) != len(authorized_ids):
        raise ValueError("Han v3 recovery current authorized partition differs")
    return completed, incomplete


def run_recovery(
    bundle: RecoveryBundle,
    *,
    expected_plan_sha256: str,
    workers: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    if not 1 <= int(workers) <= 24:
        raise ValueError("Han v3 recovery workers must be in [1, 24]")
    canonical_authorization_path = bundle.output_dir / CANONICAL_AUTHORIZATION_NAME
    if canonical_authorization_path.is_file():
        authorization_path = canonical_authorization_path
        authorization = _load_authorization(bundle, authorization_path)
        if authorization["dry_plan_payload_sha256"] != str(expected_plan_sha256):
            raise ValueError("Han v3 recovery requested initial dry-plan SHA differs")
        history_path = bundle.repo_root / str(authorization["history_manifest_path"])
        history = json.loads(history_path.read_text(encoding="utf-8"))
        if history.get("payload_sha256") != authorization["history_payload_sha256"]:
            raise ValueError("Han v3 recovery canonical history payload differs")
    else:
        plan = build_recovery_plan(bundle)
        if plan["payload_sha256"] != str(expected_plan_sha256):
            raise ValueError("Han v3 recovery state changed after the read-only dry plan")
        history_path, history = snapshot_recovery_history(bundle, plan)
        authorization_path = publish_recovery_authorization(
            bundle, plan, history_path, history
        )
        authorization = _load_authorization(bundle, authorization_path)

    precompleted, preincomplete = _authorized_current_state(bundle, authorization)
    unit_ids = {str(row["unit_id"]) for row in preincomplete}
    units = bundle.direct.units[
        bundle.direct.units["unit_id"].astype(str).isin(unit_ids)
    ].copy()
    records = units.sort_values("unit_id").to_dict(orient="records")
    if len(records) != len(preincomplete):
        raise ValueError("Han v3 recovery remaining unit set differs")
    invocation_started_ns = time_ns()
    invocation_prestate = {
        "authorization_payload_sha256": authorization["payload_sha256"],
        "already_completed": precompleted,
        "still_incomplete": preincomplete,
        "invocation_started_unix_ns": invocation_started_ns,
        "worker_count": int(workers),
    }
    invocation_id = _canonical_sha256(invocation_prestate)
    if not records:
        results: list[dict[str, Any]] = []
    elif int(workers) == 1:
        results = [
            _simulate_one_recovery(bundle, unit, authorization_path)
            for unit in records
        ]
    else:
        results = []
        with ProcessPoolExecutor(
            max_workers=int(workers),
            initializer=_worker_initialize,
            initargs=(
                str(bundle.direct.manifest_path),
                str(bundle.repo_root),
                str(bundle.direct.slim_path),
                str(authorization_path),
            ),
        ) as pool:
            futures = {
                pool.submit(_worker_simulate, unit, str(authorization_path)): unit
                for unit in records
            }
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
                            "recovery_authorization_sha256": "",
                        }
                    )
    result_columns = (
        "unit_id",
        "status",
        "completion_path",
        "elapsed_seconds",
        "error_type",
        "error",
        "recovery_authorization_sha256",
    )
    frame = pd.DataFrame(results, columns=result_columns)
    if not frame.empty:
        frame = frame.sort_values("unit_id").reset_index(drop=True)
    snapshot_path = (
        bundle.output_dir
        / "recovery_run_snapshots"
        / f"invocation_{invocation_id}.tsv"
    )
    v3._atomic_frame(snapshot_path, frame)
    failed = frame[frame["status"].astype(str).eq("failed")]
    completed_rows, incomplete_rows = _authorized_current_state(bundle, authorization)
    base_audit: dict[str, Any] | None = None
    if not incomplete_rows:
        base_audit = v3.audit_stage(
            bundle.direct, v3.STAGE_SIMULATE, write=False
        )
        if (
            base_audit.get("status") != "passed"
            or int(base_audit.get("accepted_units", -1)) != 60
        ):
            raise ValueError("Han v3 recovery post-run base audit differs")
    recovered_ids = sorted(row["unit_id"] for row in completed_rows)
    incomplete_ids = sorted(row["unit_id"] for row in incomplete_rows)
    invocation_failed_ids = sorted(failed["unit_id"].astype(str))
    authorized_ids = sorted(str(value) for value in authorization["recoverable_unit_ids"])
    if not incomplete_ids and recovered_ids != authorized_ids:
        raise ValueError("Han v3 recovery completion/authorization IDs differ")
    audit: dict[str, Any] = {
        "schema": RECOVERY_COMPLETION_AUDIT_SCHEMA,
        "status": "passed" if not incomplete_ids else "incomplete",
        "invocation_id": invocation_id,
        "invocation_started_unix_ns": invocation_started_ns,
        "dry_plan_payload_sha256": authorization["dry_plan_payload_sha256"],
        "authorization_path": authorization_path.relative_to(
            bundle.repo_root
        ).as_posix(),
        "authorization_sha256": v3.sha256_file(authorization_path),
        "canonical_authorization_path": canonical_authorization_path.relative_to(
            bundle.repo_root
        ).as_posix(),
        "canonical_authorization_sha256": v3.sha256_file(
            canonical_authorization_path
        ),
        "history_manifest_path": history_path.relative_to(
            bundle.repo_root
        ).as_posix(),
        "history_manifest_sha256": v3.sha256_file(history_path),
        "history_payload_sha256": history["payload_sha256"],
        "history_file_count": int(history["file_count"]),
        "history_total_bytes": int(history["total_bytes"]),
        "manifest_sha256": authorization["manifest_sha256"],
        "seed_plan_sha256": authorization["seed_plan_sha256"],
        "original_implementation_sources": authorization[
            "original_implementation_sources"
        ],
        "recovery_sources": authorization["recovery_sources"],
        "run_snapshot_path": snapshot_path.relative_to(bundle.repo_root).as_posix(),
        "run_snapshot_sha256": v3.sha256_file(snapshot_path),
        "initial_authorized_recovery_units": len(authorized_ids),
        "invocation_scheduled_units": len(records),
        "pre_invocation_already_recovered_unit_ids": sorted(
            row["unit_id"] for row in precompleted
        ),
        "pre_invocation_incomplete_unit_ids": sorted(
            row["unit_id"] for row in preincomplete
        ),
        "pre_recovery_incomplete_unit_ids": authorized_ids,
        "post_recovery_recovered_unit_ids": recovered_ids,
        "post_recovery_incomplete_unit_ids": incomplete_ids,
        "invocation_failed_unit_ids": invocation_failed_ids,
        "validated_recovered_completions": len(completed_rows),
        "remaining_incomplete_units": len(incomplete_ids),
        "invocation_failed_units": len(failed),
        "rows": completed_rows,
        "post_run_base_audit_status": (
            None if base_audit is None else base_audit["status"]
        ),
        "post_run_base_audit_accepted_units": (
            None if base_audit is None else int(base_audit["accepted_units"])
        ),
        "post_run_base_audit_payload_sha256": (
            None if base_audit is None else base_audit["payload_sha256"]
        ),
        "post_run_base_audit": base_audit,
        "original_attempt_seed_timeout_budget_and_reconciliation_semantics_preserved": (
            True
        ),
    }
    audit["payload_sha256"] = _canonical_sha256(audit)
    audit_path = (
        bundle.output_dir
        / "recovery_completion_audits"
        / f"invocation_{invocation_id}.json"
    )
    _immutable_json(audit_path, audit)
    if not incomplete_ids:
        _immutable_json(bundle.output_dir / CANONICAL_AUDIT_NAME, audit)
    return frame, audit


def _default_paths(args: argparse.Namespace) -> tuple[Path, Path, Path]:
    root = Path(args.repo_root).resolve()
    output = Path(args.output_dir)
    if not output.is_absolute():
        output = root / output
    slim = Path(args.slim_bin)
    if not slim.is_absolute():
        slim = root / slim
    return root, output.resolve(), slim.resolve()


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Dry-plan or run the additive Han v3 dtype recovery."
    )
    parser.add_argument("command", choices=("dry-plan", "resume", "audit"))
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument(
        "--output-dir", type=Path, default=Path(v3.DEFAULT_OUTPUT_RELATIVE_PATH)
    )
    parser.add_argument(
        "--slim-bin", type=Path, default=Path(".native-stdpopsim/bin/slim")
    )
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--expected-plan-sha256")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    root, output, slim = _default_paths(args)
    bundle = load_recovery_bundle(
        output / v3.DEFAULT_PLAN_MANIFEST_NAME,
        repo_root=root,
        slim_path=slim,
    )
    if args.command == "dry-plan":
        print(json.dumps(build_recovery_plan(bundle), indent=2, sort_keys=True))
        return 0
    if args.command == "audit":
        payload = v3.audit_stage(bundle.direct, v3.STAGE_SIMULATE)
        print(json.dumps(payload, indent=2, sort_keys=True))
        return 0
    if not args.expected_plan_sha256:
        raise ValueError("resume requires --expected-plan-sha256 from dry-plan")
    frame, audit = run_recovery(
        bundle,
        expected_plan_sha256=args.expected_plan_sha256,
        workers=args.workers,
    )
    response = {
        "status": audit["status"],
        "planned_recovery_units": len(frame),
        "completed_or_cached_units": int(
            (~frame["status"].astype(str).eq("failed")).sum()
        ),
        "failed_units": int(frame["status"].astype(str).eq("failed").sum()),
        "completion_audit_payload_sha256": audit["payload_sha256"],
    }
    print(json.dumps(response, indent=2, sort_keys=True))
    return 0 if audit["status"] == "passed" else 2


if __name__ == "__main__":
    raise SystemExit(main())
