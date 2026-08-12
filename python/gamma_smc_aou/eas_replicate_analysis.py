"""Aggregation for replicated AncientEurasia introgressed-sweep campaigns.

The unit of replication in this module is one independently seeded, accepted
SLiM trajectory.  Diploid pairs within a trajectory are never treated as
independent replicates.  All summaries are descriptive and conditional on
successful frequency conditioning; they do not provide neutral-null p-values
or estimates of statistical power.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import shutil
import socket
import subprocess
import sys
import uuid
from contextlib import contextmanager
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .eas_sweep_analysis import internal_genotype_contrasts
from .eas_sweep_study import (
    StudyRuntime,
    _atomic_frame,
    _atomic_json,
    _canonical_json,
    _center_rows,
    _contract_sha256,
    _pid_is_alive,
    _recorded_simulation_contract_if_current,
    _replace_with_access_retry,
    _sha256,
    _simulation_invocation_history,
    _stable_seed,
    _valid_recorded_decode_cache,
    _valid_recorded_simulation_cache,
    _validated_external_draw_progress,
)


REPLICATE_RESULTS_SCHEMA = "gamma-smc.eas-introgression-replicates/v1"
INPUT_SNAPSHOT_SCHEMA = "gamma-smc.eas-introgression-input-snapshot/v1"
BOOTSTRAP_DRAWS = 5_000
MIN_CELL_PAIRED_REPLICATES = 2
_METRICS = {
    "alt_minus_ref_normalized_cdf_auc": "positive",
    "alt_minus_ref_mean_tmrca_generations": "negative",
}
_SOURCES = ("tree_truth", "gamma_smc")
_FONT = {"title": 21, "axis": 18, "tick": 15, "legend": 13, "annotation": 12}
_COLORS = {0.01: "#0072B2", 0.005: "#D55E00", 0.001: "#009E73"}
_LETTER_LANDSCAPE = (11.0, 8.5)
_PROFILE_CLASSES = ("overall", "hom_ref", "hom_alt")
_REPORT_METRIC_LABELS = {
    "alt_minus_ref_normalized_cdf_auc": (
        "Normalized integrated TMRCA-CDF area (alt/alt - ref/ref; expected > 0)"
    ),
    "alt_minus_ref_mean_tmrca_generations": (
        "Mean TMRCA in generations (alt/alt - ref/ref; expected < 0)"
    ),
}
_REPORT_CLASS_LABELS = {
    "overall": "Overall",
    "hom_ref": "Ref/ref",
    "hom_alt": "Alt/alt",
    "heterozygous": "Heterozygous",
}
_AGGREGATION_LOCK_NAME = ".replicate_aggregation.lock"
_COMPACT_SPECIFICATION_INPUTS = (
    "simulation_contract.json",
    "simulation_complete.json",
    "simulation_failed.json",
    "external_draw_state.json",
    "external_draws.tsv",
    "simulation_invocations.tsv",
    "sample_manifest.tsv",
    "truth_profiles.tsv.gz",
    "truth_internal_contrasts.tsv.gz",
    "decode_contract.json",
    "decode_complete.json",
    "decoded_profiles.tsv.gz",
    "decoded_internal_contrasts.tsv.gz",
)


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


@contextmanager
def _campaign_aggregation_lock(campaign_dir: Path):
    """Serialize aggregations and recover a demonstrably stale local lock."""
    campaign_dir.mkdir(parents=True, exist_ok=True)
    lock_path = campaign_dir / _AGGREGATION_LOCK_NAME
    hostname = socket.gethostname()
    payload = {
        "pid": os.getpid(),
        "hostname": hostname,
        "created_at_utc": _utc_now(),
        "token": uuid.uuid4().hex,
    }
    for _ in range(2):
        try:
            descriptor = os.open(
                lock_path,
                os.O_CREAT | os.O_EXCL | os.O_WRONLY,
                0o600,
            )
        except FileExistsError:
            try:
                owner = json.loads(lock_path.read_text(encoding="utf-8"))
                owner_pid = int(owner["pid"])
                owner_host = str(owner["hostname"])
            except (
                KeyError,
                OSError,
                TypeError,
                ValueError,
                json.JSONDecodeError,
            ) as error:
                raise RuntimeError(
                    f"aggregation lock is unreadable and must be inspected: {lock_path}"
                ) from error
            if owner_host != hostname or _pid_is_alive(owner_pid):
                raise RuntimeError(
                    "replicate aggregation is already active at "
                    f"pid={owner_pid} host={owner_host}: {lock_path}"
                )
            stale_path = lock_path.with_name(
                f"{lock_path.name}.stale.{owner_pid}.{os.getpid()}"
            )
            try:
                _replace_with_access_retry(lock_path, stale_path)
            except FileNotFoundError:
                continue
            continue
        else:
            try:
                with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
                    json.dump(payload, handle, sort_keys=True)
                    handle.write("\n")
                yield lock_path
            finally:
                try:
                    current = json.loads(lock_path.read_text(encoding="utf-8"))
                except (OSError, ValueError, json.JSONDecodeError):
                    current = None
                if current == payload:
                    lock_path.unlink(missing_ok=True)
            return
    raise RuntimeError(f"could not acquire aggregation lock: {lock_path}")


def _assert_campaign_quiescent(campaign_dir: Path) -> None:
    """Fail closed when a simulation/decode/plan task can still mutate inputs."""
    hostname = socket.gethostname()
    lock_paths = [campaign_dir / ".study_plan_creation.lock"]
    lock_paths.extend((campaign_dir / "work").glob("*/.simulation.lock"))
    lock_paths.extend((campaign_dir / "work").glob("*/.decode.lock"))
    for lock_path in lock_paths:
        if not lock_path.is_file():
            continue
        try:
            owner = json.loads(lock_path.read_text(encoding="utf-8"))
            owner_pid = int(owner["pid"])
            owner_host = owner.get("hostname")
        except (
            KeyError,
            OSError,
            TypeError,
            ValueError,
            json.JSONDecodeError,
        ) as error:
            raise RuntimeError(
                f"campaign task lock is unreadable; snapshot refused: {lock_path}"
            ) from error
        owner_state = (
            "remote_or_unknown"
            if owner_host not in (None, hostname)
            else "active"
            if _pid_is_alive(owner_pid)
            else "stale_unreconciled"
        )
        raise RuntimeError(
            "campaign task lock must be reconciled before snapshot "
            f"(state={owner_state}, pid={owner_pid}): {lock_path}"
        )

    events: dict[str, dict[str, Any]] = {}
    for event_path in sorted((campaign_dir / "phase_invocations").glob("*.json")):
        try:
            event = json.loads(event_path.read_text(encoding="utf-8"))
            invocation_id = str(event["invocation_id"])
            event_name = str(event["event"])
        except (
            KeyError,
            OSError,
            TypeError,
            ValueError,
            json.JSONDecodeError,
        ) as error:
            raise RuntimeError(
                f"phase invocation event is unreadable; snapshot refused: {event_path}"
            ) from error
        state = events.setdefault(invocation_id, {"started": None, "terminal": False})
        if event_name == "started":
            if state["started"] is not None:
                raise RuntimeError(
                    "phase invocation contains duplicate started events; snapshot "
                    f"refused: invocation_id={invocation_id}"
                )
            state["started"] = event
        elif event_name in {"completed", "failed"}:
            if state["terminal"]:
                raise RuntimeError(
                    "phase invocation contains multiple terminal events; snapshot "
                    f"refused: invocation_id={invocation_id}"
                )
            state["terminal"] = True
        else:
            raise RuntimeError(
                "phase invocation contains an unsupported event; snapshot refused: "
                f"invocation_id={invocation_id} event={event_name}"
            )
    for invocation_id, state in events.items():
        started = state["started"]
        if started is None:
            raise RuntimeError(
                "phase invocation has a terminal event without a start; snapshot "
                f"refused: invocation_id={invocation_id}"
            )
        if state["terminal"]:
            continue
        owner_pid = int(started.get("process_id", -1))
        if owner_pid == os.getpid():
            # The orchestrator records its plot/all event before entering this
            # function; its preceding actions are synchronous and have returned.
            continue
        raise RuntimeError(
            "unterminated campaign phase must be reconciled before snapshot: "
            f"invocation_id={invocation_id} pid={owner_pid}"
        )


def _snapshot_input_paths(
    campaign_dir: Path,
    expanded: Sequence[Mapping[str, Any]],
) -> dict[Path, str]:
    """Return every compact analysis input plus hash-only large tree inputs."""
    paths: dict[Path, str] = {}

    def include(path: Path, mode: str = "copied") -> None:
        if path.is_file():
            paths[path.resolve()] = mode

    for relative in (
        "study_design.json",
        "specifications.tsv",
        "results/simulation_status.tsv",
        "results/decode_status.tsv",
    ):
        include(campaign_dir / relative)
    for path in sorted((campaign_dir / "phase_invocations").glob("*.json")):
        include(path)
    for path in sorted(campaign_dir.glob(f"{_AGGREGATION_LOCK_NAME}.stale.*")):
        include(path)
    for specification in expanded:
        specification_dir = (
            campaign_dir / "work" / str(specification["specification_id"])
        )
        for filename in _COMPACT_SPECIFICATION_INPUTS:
            include(specification_dir / filename)
        for subdirectory in ("pairs", "decoded"):
            for path in sorted((specification_dir / subdirectory).glob("*")):
                include(path)
        for pattern in (".simulation.lock.stale.*", ".decode.lock.stale.*"):
            for path in sorted(specification_dir.glob(pattern)):
                include(path)
        include(specification_dir / "selected.trees", mode="hash_only")
    return paths


def _portable_relative(path: Path, root: Path) -> str:
    return path.resolve().relative_to(root.resolve()).as_posix()


def _build_input_snapshot(
    campaign_dir: Path,
    expanded: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Materialize a content-addressed snapshot before any aggregation."""
    provenance_dir = campaign_dir / "results" / "provenance"
    snapshots_dir = provenance_dir / "input_snapshots"
    provenance_dir.mkdir(parents=True, exist_ok=True)
    snapshots_dir.mkdir(parents=True, exist_ok=True)
    staging = provenance_dir / f".input_snapshot.{os.getpid()}.{uuid.uuid4().hex}"
    staging.mkdir()
    records: list[dict[str, Any]] = []
    try:
        for source, mode in sorted(
            _snapshot_input_paths(campaign_dir, expanded).items(),
            key=lambda item: _portable_relative(item[0], campaign_dir),
        ):
            relative = _portable_relative(source, campaign_dir)
            before_stat = source.stat()
            source_sha256 = _sha256(source)
            snapshot_path: str | None = None
            if mode == "copied":
                destination = staging / relative
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(source, destination)
                if _sha256(destination) != source_sha256:
                    raise RuntimeError(f"snapshot copy checksum mismatch: {source}")
                snapshot_path = relative
            after_stat = source.stat()
            if (
                before_stat.st_size != after_stat.st_size
                or before_stat.st_mtime_ns != after_stat.st_mtime_ns
                or _sha256(source) != source_sha256
            ):
                raise RuntimeError(f"campaign input changed during snapshot: {source}")
            records.append(
                {
                    "original_path": relative,
                    "snapshot_path": snapshot_path,
                    "snapshot_mode": mode,
                    "sha256": source_sha256,
                    "size_bytes": int(after_stat.st_size),
                    "mtime_ns_at_snapshot": int(after_stat.st_mtime_ns),
                }
            )
        fingerprint_records = [
            {
                key: record[key]
                for key in (
                    "original_path",
                    "snapshot_path",
                    "snapshot_mode",
                    "sha256",
                    "size_bytes",
                    "mtime_ns_at_snapshot",
                )
            }
            for record in records
        ]
        snapshot_id = _contract_sha256(fingerprint_records)
        payload = {
            "schema": INPUT_SNAPSHOT_SCHEMA,
            "snapshot_id": snapshot_id,
            "created_at_utc": _utc_now(),
            "campaign_root": ".",
            "files": records,
        }
        inventory_path = staging / "input_snapshot_manifest.json"
        inventory_path.write_text(
            json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        destination = snapshots_dir / snapshot_id
        if destination.exists():
            existing = json.loads(
                (destination / "input_snapshot_manifest.json").read_text(
                    encoding="utf-8"
                )
            )
            if existing.get("snapshot_id") != snapshot_id:
                raise RuntimeError(f"snapshot ID collision: {destination}")
            shutil.rmtree(staging)
            payload = existing
        else:
            _replace_with_access_retry(staging, destination)
        return {
            "snapshot_id": snapshot_id,
            "root": destination,
            "manifest": payload,
            "manifest_path": destination / "input_snapshot_manifest.json",
        }
    except BaseException:
        if staging.exists():
            shutil.rmtree(staging)
        raise


def _verify_snapshot_sources_unchanged(
    campaign_dir: Path,
    expanded: Sequence[Mapping[str, Any]],
    snapshot: Mapping[str, Any],
) -> None:
    """Prove the snapshotted source set stayed fixed through aggregation."""
    current = _snapshot_input_paths(campaign_dir, expanded)
    current_relatives = {
        _portable_relative(path, campaign_dir): mode for path, mode in current.items()
    }
    records = snapshot["manifest"]["files"]
    expected_relatives = {
        str(record["original_path"]): str(record["snapshot_mode"]) for record in records
    }
    if current_relatives != expected_relatives:
        raise RuntimeError("campaign input file set changed after snapshot")
    for record in records:
        source = campaign_dir / str(record["original_path"])
        stat = source.stat()
        if stat.st_size != int(record["size_bytes"]) or stat.st_mtime_ns != int(
            record["mtime_ns_at_snapshot"]
        ):
            raise RuntimeError(f"campaign input changed after snapshot: {source}")
        if _sha256(source) != record["sha256"]:
            raise RuntimeError(
                f"campaign input checksum changed after snapshot: {source}"
            )


def _same_biological_value(left: Any, right: Any) -> bool:
    if left is None or right is None:
        return left is right
    if isinstance(left, (int, float, np.number)) and isinstance(
        right, (int, float, np.number)
    ):
        return bool(np.isclose(float(left), float(right), rtol=0, atol=1e-12))
    return left == right


def _normalise_specifications(
    base_specs: Sequence[Mapping[str, Any]],
    expanded_specs: Sequence[Mapping[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    bases = [dict(record) for record in base_specs]
    if not bases:
        raise ValueError("base_specs must contain at least one biological cell")
    required_base = {
        "specification_id",
        "study_type",
        "selection_coefficient",
        "target_allele_frequency",
    }
    for base in bases:
        missing = sorted(required_base.difference(base))
        if missing:
            raise ValueError(f"base specification is missing: {', '.join(missing)}")
        if base["study_type"] != "introgression":
            raise ValueError("replicate aggregation only accepts introgression specs")
        selection = float(base["selection_coefficient"])
        frequency = float(base["target_allele_frequency"])
        if not np.isfinite(selection) or selection <= 0:
            raise ValueError("selection coefficients must be finite and positive")
        if not np.isfinite(frequency) or not 0 < frequency < 1:
            raise ValueError(
                "target allele frequencies must lie strictly within (0, 1)"
            )
    base_lookup = {str(base["specification_id"]): base for base in bases}
    if len(base_lookup) != len(bases):
        raise ValueError("base specification IDs must be unique")

    biological_fields = (
        "study_type",
        "selection_coefficient",
        "target_allele_frequency",
        "lower_allele_frequency",
        "upper_allele_frequency",
        "origin_age_generations",
        "source_minimum_frequency",
        "slim_scaled_source_check_generations_ago",
        "slim_scaled_pulse_generations_ago",
    )
    expanded: list[dict[str, Any]] = []
    seen_ids: set[str] = set()
    seen_replicate_ids: set[str] = set()
    seen_indices: set[tuple[str, int]] = set()
    for supplied in expanded_specs:
        supplied = dict(supplied)
        if "base_specification_id" not in supplied:
            raise ValueError("expanded specification lacks base_specification_id")
        base_id = str(supplied["base_specification_id"])
        if base_id not in base_lookup:
            raise ValueError(f"unknown base_specification_id: {base_id}")
        if "replicate_index" not in supplied:
            raise ValueError("expanded specification lacks replicate_index")
        if isinstance(supplied["replicate_index"], bool):
            raise ValueError("replicate_index must be a non-negative integer")
        replicate_index = int(supplied["replicate_index"])
        if replicate_index < 0 or replicate_index != float(supplied["replicate_index"]):
            raise ValueError("replicate_index must be a non-negative integer")
        specification_id = str(
            supplied.get("specification_id", supplied.get("replicate_id", ""))
        )
        if not specification_id:
            raise ValueError("expanded specification lacks specification_id")
        replicate_id = str(supplied.get("replicate_id", specification_id))
        if specification_id in seen_ids:
            raise ValueError(f"duplicate expanded specification ID: {specification_id}")
        if replicate_id in seen_replicate_ids:
            raise ValueError(f"duplicate replicate_id: {replicate_id}")
        index_key = (base_id, replicate_index)
        if index_key in seen_indices:
            raise ValueError(
                f"duplicate replicate_index {replicate_index} for {base_id}"
            )

        base = base_lookup[base_id]
        for field in biological_fields:
            if (
                field in supplied
                and field in base
                and not _same_biological_value(supplied[field], base[field])
            ):
                raise ValueError(
                    f"expanded {specification_id} changes biological field {field}"
                )
        record = {**base, **supplied}
        record.update(
            {
                "base_specification_id": base_id,
                "replicate_index": replicate_index,
                "replicate_id": replicate_id,
                "specification_id": specification_id,
            }
        )
        expanded.append(record)
        seen_ids.add(specification_id)
        seen_replicate_ids.add(replicate_id)
        seen_indices.add(index_key)
    return bases, sorted(
        expanded,
        key=lambda row: (row["base_specification_id"], row["replicate_index"]),
    )


def _failure_status(
    specification_dir: Path,
    specification: Mapping[str, Any],
    runtime: StudyRuntime,
) -> tuple[str, str | None, str | None]:
    """Return only contract-current failure information."""
    failure_path = specification_dir / "simulation_failed.json"
    contract = _recorded_simulation_contract_if_current(
        specification_dir, specification, runtime
    )
    if failure_path.is_file() and contract is not None:
        try:
            failure = json.loads(failure_path.read_text(encoding="utf-8"))
            if failure.get("simulation_contract_sha256") == _contract_sha256(contract):
                status = str(failure.get("status", "failed"))
                return status, failure.get("error_type"), failure.get("error")
        except (OSError, ValueError, json.JSONDecodeError):
            pass
    evidence_names = (
        "simulation_complete.json",
        "simulation_failed.json",
        "selected.trees",
        "truth_profiles.tsv.gz",
    )
    if any((specification_dir / name).exists() for name in evidence_names):
        return (
            "stale_or_invalid",
            "CacheValidationError",
            "on-disk simulation artifacts do not satisfy the current contract",
        )
    return "not_requested", None, None


def _decode_status_evidence(
    specification_dir: Path,
    simulation_complete: bool,
    decode: Mapping[str, Any] | None,
    recorded_status: Mapping[str, Any] | None,
) -> tuple[str, str | None, str | None]:
    """Classify decode state without treating every invalid cache as unrequested."""
    if not simulation_complete:
        return "not_simulated", None, None
    if decode is not None:
        return "complete", None, None

    recorded_status = dict(recorded_status or {})
    recorded = str(recorded_status.get("status", "")).strip().lower()
    recorded_error_type = recorded_status.get("error_type")
    recorded_detail = recorded_status.get("error")
    decode_complete = (specification_dir / "decode_complete.json").is_file()
    decode_contract = (specification_dir / "decode_contract.json").is_file()
    derivative_paths = (
        specification_dir / "decoded_profiles.tsv.gz",
        specification_dir / "decoded_internal_contrasts.tsv.gz",
    )
    class_artifacts = any(
        path.is_file() for path in (specification_dir / "decoded").glob("*")
    )
    has_any_evidence = (
        decode_complete
        or decode_contract
        or any(path.is_file() for path in derivative_paths)
        or class_artifacts
    )
    if recorded == "failed":
        return "failed", recorded_error_type, recorded_detail
    if decode_complete or recorded == "complete":
        return (
            "stale_or_incompatible",
            "DecodeCacheValidationError",
            "decode completion exists but does not satisfy the requested contract",
        )
    if has_any_evidence:
        return (
            "incomplete_or_corrupt",
            "DecodeCacheValidationError",
            "partial decode artifacts exist without a valid completion",
        )
    if recorded in {"pending", "queued", "running"}:
        return "pending", recorded_error_type, recorded_detail
    if recorded == "not_requested":
        return "not_requested", None, None
    if recorded == "not_simulated":
        return (
            "pending",
            None,
            "simulation became available after the recorded decode-status snapshot",
        )
    return "not_requested", None, None


def _recorded_decode_status(source_study_dir: Path) -> dict[str, dict[str, Any]]:
    path = source_study_dir / "results" / "decode_status.tsv"
    if not path.is_file():
        return {}
    try:
        frame = pd.read_csv(path, sep="\t")
    except (OSError, pd.errors.EmptyDataError, pd.errors.ParserError) as error:
        raise RuntimeError(f"decode status table is unreadable: {path}") from error
    if "specification_id" not in frame or frame["specification_id"].duplicated().any():
        raise RuntimeError(
            "decode status table must contain one row per specification_id"
        )
    return frame.set_index("specification_id").to_dict(orient="index")


def _safe_number(mapping: Mapping[str, Any], key: str) -> Any:
    value = mapping.get(key)
    return None if value is None else value


def _metadata(record: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "base_specification_id": record["base_specification_id"],
        "replicate_id": record["replicate_id"],
        "replicate_index": int(record["replicate_index"]),
        "specification_id": record["specification_id"],
        "selection_coefficient": float(record["selection_coefficient"]),
        "target_allele_frequency": float(record["target_allele_frequency"]),
    }


def _read_center_profile(path: Path, runtime: StudyRuntime) -> pd.DataFrame:
    profile = pd.read_csv(path, sep="\t")
    if profile.empty:
        raise ValueError(f"validated profile is empty: {path}")
    return _center_rows(profile, runtime.focal_position_bp)


def _external_draw_counts(
    draw_state: Mapping[str, Any], attempts: Sequence[Mapping[str, Any]]
) -> tuple[int, int]:
    active_draw = draw_state.get("active_internal_conditioned_trajectory")
    launched_draws = int(
        draw_state.get("next_external_draw_zero_based", len(attempts))
    ) + int(active_draw is not None)
    # An active trajectory is launched, but it is not interrupted until the
    # simulation state records it in the terminal interruption ledger.
    interrupted_draws = len(
        draw_state.get("interrupted_internal_conditioned_trajectories", [])
    )
    return launched_draws, interrupted_draws


def _collect_replicates(
    runtime: StudyRuntime,
    expanded: Sequence[Mapping[str, Any]],
    study_dir: Path,
    decoder_path: str | Path | None,
    *,
    source_study_dir: Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    source_study_dir = source_study_dir or study_dir
    decode_status_lookup = _recorded_decode_status(source_study_dir)
    status_rows: list[dict[str, Any]] = []
    acceptance_rows: list[dict[str, Any]] = []
    contrast_frames: list[pd.DataFrame] = []
    profile_frames: list[pd.DataFrame] = []
    for specification in expanded:
        meta = _metadata(specification)
        specification_dir = study_dir / "work" / str(specification["specification_id"])
        source_specification_dir = (
            source_study_dir / "work" / str(specification["specification_id"])
        )
        recorded_contract = _recorded_simulation_contract_if_current(
            specification_dir, specification, runtime
        )
        progress = (
            _validated_external_draw_progress(
                source_specification_dir,
                recorded_contract,
                specification,
                runtime,
            )
            if recorded_contract is not None
            else None
        )
        draw_state, attempts = progress if progress is not None else ({}, [])
        launched_draws, interrupted_draws = _external_draw_counts(draw_state, attempts)
        invocation_history = (
            _simulation_invocation_history(
                source_specification_dir,
                expected_contract_sha256=_contract_sha256(recorded_contract),
                allow_legacy_migration=True,
            )
            if recorded_contract is not None
            else pd.DataFrame()
        )
        cumulative_elapsed = float(
            pd.to_numeric(
                invocation_history.get("elapsed_seconds", pd.Series(dtype=float)),
                errors="coerce",
            )
            .fillna(0.0)
            .sum()
        )
        simulation = _valid_recorded_simulation_cache(
            specification_dir, specification, runtime
        )
        if simulation is None:
            sim_status, error_type, detail = _failure_status(
                source_specification_dir, specification, runtime
            )
            decode = None
        else:
            sim_status, error_type, detail = "complete", None, None
            decode = _valid_recorded_decode_cache(
                specification_dir,
                specification,
                runtime,
                decoder_path=decoder_path,
            )
        decode_status, decode_error_type, decode_detail = _decode_status_evidence(
            source_specification_dir,
            simulation is not None,
            decode,
            decode_status_lookup.get(str(specification["specification_id"])),
        )
        status_rows.append(
            {
                **meta,
                "simulation_status": sim_status,
                "decode_status": decode_status,
                "error_type": error_type,
                "detail": detail,
                "decode_error_type": decode_error_type,
                "decode_detail": decode_detail,
                "invocation_count": int(len(invocation_history)),
                "cumulative_elapsed_seconds": cumulative_elapsed,
                "completed_external_draws": int(len(attempts)),
                "launched_external_draws": launched_draws,
                "interrupted_internal_conditioned_trajectories": interrupted_draws,
                "last_completed_draw_sample_af": (
                    attempts[-1].get("achieved_allele_frequency") if attempts else None
                ),
                "last_completed_draw_n_hom_ref": (
                    attempts[-1].get("n_hom_ref") if attempts else None
                ),
                "conditioning_missingness": "nonrandom",
                "status": sim_status,
            }
        )
        if simulation is None:
            continue

        validation = simulation["validation"]
        acceptance_rows.append(
            {
                **meta,
                "achieved_allele_frequency": float(
                    validation["achieved_allele_frequency"]
                ),
                "achieved_allele_frequency_scope": validation.get(
                    "achieved_allele_frequency_scope", "sample_estimate"
                ),
                "accepted_seed": int(simulation["accepted_seed"]),
                "attempts_to_accept": int(simulation["attempts_to_accept"]),
                "completed_external_sample_panel_draws": int(
                    simulation.get(
                        "completed_external_sample_panel_draws",
                        simulation["attempts_to_accept"],
                    )
                ),
                "interrupted_internal_conditioned_trajectories": int(
                    simulation.get("interrupted_internal_conditioned_trajectories", 0)
                ),
                "launched_external_draws": launched_draws,
                "invocation_count": int(len(invocation_history)),
                "cumulative_elapsed_seconds": cumulative_elapsed,
                "n_hom_ref": int(validation["n_hom_ref"]),
                "n_heterozygous": int(validation["n_heterozygous"]),
                "n_hom_alt": int(validation["n_hom_alt"]),
                "decoded": decode is not None,
                "conditioning_missingness": "nonrandom",
            }
        )

        truth = _read_center_profile(
            source_specification_dir / "truth_profiles.tsv.gz", runtime
        )
        source_profiles = [truth]
        if decode is not None:
            source_profiles.append(
                _read_center_profile(
                    source_specification_dir / "decoded_profiles.tsv.gz",
                    runtime,
                )
            )
        for profile in source_profiles:
            if profile["source"].nunique() != 1:
                raise ValueError(
                    f"one replicate profile must contain one source: {specification_dir}"
                )
            enriched = profile.copy()
            for key, value in meta.items():
                enriched[key] = value
            profile_frames.append(enriched)
            contrasts = internal_genotype_contrasts(enriched)
            if len(contrasts) != 1:
                raise ValueError(
                    "each replicate/source must yield exactly one center contrast"
                )
            for key, value in meta.items():
                contrasts[key] = value
            contrast_frames.append(contrasts)

    status_columns = [
        "base_specification_id",
        "replicate_id",
        "replicate_index",
        "specification_id",
        "selection_coefficient",
        "target_allele_frequency",
        "simulation_status",
        "decode_status",
        "error_type",
        "detail",
        "decode_error_type",
        "decode_detail",
        "invocation_count",
        "cumulative_elapsed_seconds",
        "completed_external_draws",
        "launched_external_draws",
        "interrupted_internal_conditioned_trajectories",
        "last_completed_draw_sample_af",
        "last_completed_draw_n_hom_ref",
        "conditioning_missingness",
        "status",
    ]
    metadata_columns = list(_metadata(expanded[0]).keys()) if expanded else []
    acceptance_columns = [
        *metadata_columns,
        "achieved_allele_frequency",
        "achieved_allele_frequency_scope",
        "accepted_seed",
        "attempts_to_accept",
        "completed_external_sample_panel_draws",
        "interrupted_internal_conditioned_trajectories",
        "launched_external_draws",
        "invocation_count",
        "cumulative_elapsed_seconds",
        "n_hom_ref",
        "n_heterozygous",
        "n_hom_alt",
        "decoded",
        "conditioning_missingness",
    ]
    status = pd.DataFrame(status_rows, columns=status_columns)
    acceptance = pd.DataFrame(acceptance_rows, columns=acceptance_columns)
    contrasts = (
        pd.concat(contrast_frames, ignore_index=True)
        if contrast_frames
        else pd.DataFrame()
    )
    profiles = (
        pd.concat(profile_frames, ignore_index=True)
        if profile_frames
        else pd.DataFrame()
    )
    return status, acceptance, contrasts, profiles


def _cell_coverage(
    bases: Sequence[Mapping[str, Any]],
    expanded: Sequence[Mapping[str, Any]],
    status: pd.DataFrame,
) -> pd.DataFrame:
    target_counts: dict[str, int] = {}
    for record in expanded:
        base_id = str(record["base_specification_id"])
        target_counts[base_id] = target_counts.get(base_id, 0) + 1
    rows = []
    for base in bases:
        base_id = str(base["specification_id"])
        selected = status[status["base_specification_id"] == base_id]
        accepted = int((selected["simulation_status"] == "complete").sum())
        decoded = int((selected["decode_status"] == "complete").sum())
        target = int(target_counts.get(base_id, 0))
        counts = {
            str(key): int(value)
            for key, value in selected["simulation_status"].value_counts().items()
        }
        completed_draws = int(selected["completed_external_draws"].sum())
        launched_draws = int(selected["launched_external_draws"].sum())
        interrupted_draws = int(
            selected["interrupted_internal_conditioned_trajectories"].sum()
        )
        cumulative_elapsed = float(selected["cumulative_elapsed_seconds"].sum())
        attempted = int(
            (
                (selected["invocation_count"] > 0)
                | (selected["simulation_status"] == "complete")
            ).sum()
        )
        rows.append(
            {
                "base_specification_id": base_id,
                "selection_coefficient": float(base["selection_coefficient"]),
                "target_allele_frequency": float(base["target_allele_frequency"]),
                "target_replicates": target,
                "accepted_replicates": accepted,
                "decoded_replicates": decoded,
                "accepted_fraction_of_target": accepted / target if target else np.nan,
                "decoded_fraction_of_accepted": decoded / accepted
                if accepted
                else np.nan,
                "attempted_replicates": attempted,
                "timed_out_replicates": counts.get("timed_out", 0),
                "exhausted_replicates": counts.get("exhausted", 0),
                "failed_replicates": counts.get("failed", 0),
                "completed_external_draws": completed_draws,
                "launched_external_draws": launched_draws,
                "interrupted_internal_conditioned_trajectories": interrupted_draws,
                "cumulative_elapsed_seconds": cumulative_elapsed,
                "simulation_status_counts_json": _canonical_json(counts),
                "conditioning_missingness": "nonrandom",
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["selection_coefficient", "target_allele_frequency"],
        ascending=[False, True],
    )


def _wilson_interval(successes: int, trials: int) -> tuple[float, float]:
    if trials == 0:
        return np.nan, np.nan
    z = 1.959963984540054
    p = successes / trials
    denominator = 1 + z * z / trials
    center = (p + z * z / (2 * trials)) / denominator
    half = z * math.sqrt(p * (1 - p) / trials + z * z / (4 * trials**2))
    half /= denominator
    return center - half, center + half


def _bootstrap_mean_interval(
    values: np.ndarray, *, seed: int, draws: int = BOOTSTRAP_DRAWS
) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return np.nan, np.nan
    if len(values) == 1:
        return float(values[0]), float(values[0])
    generator = np.random.default_rng(seed)
    means = np.empty(draws, dtype=float)
    # Chunking avoids an unnecessarily large n_replicates x draws matrix.
    for start in range(0, draws, 1_000):
        stop = min(start + 1_000, draws)
        indices = generator.integers(0, len(values), size=(stop - start, len(values)))
        means[start:stop] = values[indices].mean(axis=1)
    low, high = np.quantile(means, [0.025, 0.975])
    return float(low), float(high)


def _effect_summary(
    runtime: StudyRuntime,
    bases: Sequence[Mapping[str, Any]],
    coverage: pd.DataFrame,
    contrasts: pd.DataFrame,
) -> pd.DataFrame:
    lookup = coverage.set_index("base_specification_id").to_dict(orient="index")
    rows: list[dict[str, Any]] = []
    for base in bases:
        base_id = str(base["specification_id"])
        for source in _SOURCES:
            for metric, expected_sign in _METRICS.items():
                if contrasts.empty:
                    values = np.array([], dtype=float)
                else:
                    mask = (contrasts["base_specification_id"] == base_id) & (
                        contrasts["source"] == source
                    )
                    values = contrasts.loc[mask, metric].to_numpy(dtype=float)
                    values = values[np.isfinite(values)]
                successes = int(
                    np.sum(values > 0)
                    if expected_sign == "positive"
                    else np.sum(values < 0)
                )
                sign_low, sign_high = _wilson_interval(successes, len(values))
                seed = _stable_seed(
                    runtime.base_seed,
                    f"replicate-bootstrap:{base_id}:{source}:{metric}",
                    0,
                )
                boot_low, boot_high = _bootstrap_mean_interval(values, seed=seed)
                q25, q75 = (
                    tuple(map(float, np.quantile(values, [0.25, 0.75])))
                    if len(values)
                    else (np.nan, np.nan)
                )
                cell = lookup[base_id]
                rows.append(
                    {
                        "base_specification_id": base_id,
                        "selection_coefficient": float(base["selection_coefficient"]),
                        "target_allele_frequency": float(
                            base["target_allele_frequency"]
                        ),
                        "source": source,
                        "metric": metric,
                        "expected_sweep_sign": expected_sign,
                        "target_replicates": int(cell["target_replicates"]),
                        "accepted_replicates": int(cell["accepted_replicates"]),
                        "decoded_replicates": int(cell["decoded_replicates"]),
                        "n_observed_replicates": int(len(values)),
                        "mean": float(np.mean(values)) if len(values) else np.nan,
                        "median": float(np.median(values)) if len(values) else np.nan,
                        "sample_sd": (
                            float(np.std(values, ddof=1)) if len(values) > 1 else np.nan
                        ),
                        "q25": q25,
                        "q75": q75,
                        "iqr": q75 - q25 if len(values) else np.nan,
                        "bootstrap_mean_ci95_low": boot_low,
                        "bootstrap_mean_ci95_high": boot_high,
                        "bootstrap_draws": BOOTSTRAP_DRAWS,
                        "bootstrap_seed": seed,
                        "expected_sign_count": successes,
                        "expected_sign_rate": successes / len(values)
                        if len(values)
                        else np.nan,
                        "expected_sign_wilson95_low": sign_low,
                        "expected_sign_wilson95_high": sign_high,
                        "conditioning_missingness": "nonrandom",
                    }
                )
    return pd.DataFrame(rows)


def _correlation(left: np.ndarray, right: np.ndarray, *, rank: bool) -> float:
    if len(left) < 2:
        return np.nan
    if rank:
        left = pd.Series(left).rank(method="average").to_numpy(dtype=float)
        right = pd.Series(right).rank(method="average").to_numpy(dtype=float)
    if np.isclose(np.std(left), 0) or np.isclose(np.std(right), 0):
        return np.nan
    return float(np.corrcoef(left, right)[0, 1])


def _paired_contrasts(contrasts: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "base_specification_id",
        "replicate_id",
        "replicate_index",
        "specification_id",
        "selection_coefficient",
        "target_allele_frequency",
    ]
    if contrasts.empty:
        return pd.DataFrame(
            columns=[
                *columns,
                "metric",
                "tree_truth",
                "gamma_smc",
                "tree_truth_expected_sign",
                "gamma_smc_expected_sign",
                "raw_truth_gamma_sign_agreement",
                "gamma_detected_given_truth_expected_sign",
            ]
        )
    long = contrasts.melt(
        id_vars=[*columns, "source"],
        value_vars=list(_METRICS),
        var_name="metric",
        value_name="effect",
    )
    duplicates = long.duplicated([*columns, "metric", "source"], keep=False)
    if duplicates.any():
        raise ValueError("replicate contrasts are not one observation per source")
    paired = long.pivot(index=[*columns, "metric"], columns="source", values="effect")
    paired = paired.reset_index().rename_axis(columns=None)
    for source in _SOURCES:
        if source not in paired:
            paired[source] = np.nan
    paired = paired.replace([np.inf, -np.inf], np.nan)
    paired = paired.dropna(subset=list(_SOURCES)).reset_index(drop=True)
    expected_positive = paired["metric"].map(_METRICS).eq("positive")
    paired["tree_truth_expected_sign"] = np.where(
        expected_positive,
        paired["tree_truth"] > 0,
        paired["tree_truth"] < 0,
    )
    paired["gamma_smc_expected_sign"] = np.where(
        expected_positive,
        paired["gamma_smc"] > 0,
        paired["gamma_smc"] < 0,
    )
    paired["raw_truth_gamma_sign_agreement"] = np.sign(paired["tree_truth"]) == np.sign(
        paired["gamma_smc"]
    )
    paired["gamma_detected_given_truth_expected_sign"] = (
        paired["tree_truth_expected_sign"] & paired["gamma_smc_expected_sign"]
    )
    return paired


def _trajectory_accuracy_record(
    *,
    scope: str,
    base: Mapping[str, Any] | None,
    metric: str,
    expected_sign: str,
    selected: pd.DataFrame,
    aggregation_weighting: str,
    is_headline: bool,
    total_biological_cells: int,
) -> dict[str, Any]:
    truth = selected.get("tree_truth", pd.Series(dtype=float)).to_numpy(dtype=float)
    gamma = selected.get("gamma_smc", pd.Series(dtype=float)).to_numpy(dtype=float)
    errors = gamma - truth
    truth_expected = truth > 0 if expected_sign == "positive" else truth < 0
    gamma_expected = gamma > 0 if expected_sign == "positive" else gamma < 0
    truth_count = int(np.sum(truth_expected))
    gamma_count = int(np.sum(gamma_expected))
    raw_agreement = int(np.sum(np.sign(truth) == np.sign(gamma)))
    detected = int(np.sum(truth_expected & gamma_expected))
    truth_low, truth_high = _wilson_interval(truth_count, len(truth))
    gamma_low, gamma_high = _wilson_interval(gamma_count, len(gamma))
    agreement_low, agreement_high = _wilson_interval(raw_agreement, len(truth))
    detection_low, detection_high = _wilson_interval(detected, truth_count)
    base_id = None if base is None else str(base["specification_id"])
    cells_with_pairs = (
        int(selected["base_specification_id"].nunique()) if not selected.empty else 0
    )
    return {
        "summary_scope": scope,
        "aggregation_weighting": aggregation_weighting,
        "is_headline": is_headline,
        "minimum_paired_replicates_per_cell": (
            MIN_CELL_PAIRED_REPLICATES if is_headline else np.nan
        ),
        "n_total_biological_cells": total_biological_cells,
        "n_biological_cells_with_pairs": cells_with_pairs,
        "n_biological_cells_meeting_minimum": (
            cells_with_pairs if is_headline else np.nan
        ),
        "n_biological_cells_in_detection_rate": (
            cells_with_pairs if truth_count else 0
        ),
        "base_specification_id": base_id,
        "selection_coefficient": (
            np.nan if base is None else float(base["selection_coefficient"])
        ),
        "target_allele_frequency": (
            np.nan if base is None else float(base["target_allele_frequency"])
        ),
        "metric": metric,
        "n_paired_replicates": int(len(truth)),
        "tree_truth_mean": float(np.mean(truth)) if len(truth) else np.nan,
        "gamma_smc_mean": float(np.mean(gamma)) if len(gamma) else np.nan,
        "gamma_minus_truth_bias": (float(np.mean(errors)) if len(errors) else np.nan),
        "mae": float(np.mean(np.abs(errors))) if len(errors) else np.nan,
        "rmse": float(np.sqrt(np.mean(errors**2))) if len(errors) else np.nan,
        "pearson": _correlation(truth, gamma, rank=False),
        "spearman": _correlation(truth, gamma, rank=True),
        "truth_expected_sign_count": truth_count,
        "truth_expected_sign_rate": (
            truth_count / len(truth) if len(truth) else np.nan
        ),
        "truth_expected_sign_wilson95_low": truth_low,
        "truth_expected_sign_wilson95_high": truth_high,
        "gamma_expected_sign_count": gamma_count,
        "gamma_expected_sign_rate": (
            gamma_count / len(gamma) if len(gamma) else np.nan
        ),
        "gamma_expected_sign_wilson95_low": gamma_low,
        "gamma_expected_sign_wilson95_high": gamma_high,
        "raw_truth_gamma_sign_agreement_count": raw_agreement,
        "raw_truth_gamma_sign_agreement_rate": (
            raw_agreement / len(truth) if len(truth) else np.nan
        ),
        "raw_truth_gamma_sign_agreement_wilson95_low": agreement_low,
        "raw_truth_gamma_sign_agreement_wilson95_high": agreement_high,
        "truth_expected_sign_denominator": truth_count,
        "gamma_detection_given_truth_expected_sign_count": detected,
        "gamma_detection_given_truth_expected_sign_rate": (
            detected / truth_count if truth_count else np.nan
        ),
        "gamma_detection_given_truth_expected_sign_wilson95_low": detection_low,
        "gamma_detection_given_truth_expected_sign_wilson95_high": detection_high,
        "interval_method": "wilson_trajectory_binomial",
        "raw_sign_agreement_rule": "sign(tree_truth)==sign(gamma_smc)",
        "detection_rule": (
            "Gamma has expected sweep sign among replicates whose tree truth has "
            "the expected sweep sign; jointly wrong signs never enter the numerator "
            "or denominator"
        ),
        "pairing_rule": "one_truth_and_one_gamma_value_per_replicate",
    }


def _macro_cell_equal_accuracy(
    metric: str,
    cell_rows: Sequence[Mapping[str, Any]],
    *,
    total_biological_cells: int,
) -> dict[str, Any]:
    eligible = [
        dict(row)
        for row in cell_rows
        if int(row["n_paired_replicates"]) >= MIN_CELL_PAIRED_REPLICATES
    ]

    def mean_field(field: str, rows: Sequence[Mapping[str, Any]] = eligible) -> float:
        values = np.asarray([row[field] for row in rows], dtype=float)
        values = values[np.isfinite(values)]
        return float(np.mean(values)) if len(values) else np.nan

    detection_rows = [
        row
        for row in eligible
        if np.isfinite(float(row["gamma_detection_given_truth_expected_sign_rate"]))
    ]
    paired_count = int(sum(int(row["n_paired_replicates"]) for row in eligible))
    return {
        "summary_scope": "overall_macro_cell_equal",
        "aggregation_weighting": (
            "biological_cell_equal_among_cells_with_"
            f"n_paired>={MIN_CELL_PAIRED_REPLICATES}"
        ),
        "is_headline": True,
        "minimum_paired_replicates_per_cell": MIN_CELL_PAIRED_REPLICATES,
        "n_total_biological_cells": total_biological_cells,
        "n_biological_cells_with_pairs": int(
            sum(int(row["n_paired_replicates"]) > 0 for row in cell_rows)
        ),
        "n_biological_cells_meeting_minimum": len(eligible),
        "n_biological_cells_in_detection_rate": len(detection_rows),
        "base_specification_id": None,
        "selection_coefficient": np.nan,
        "target_allele_frequency": np.nan,
        "metric": metric,
        "n_paired_replicates": paired_count,
        "tree_truth_mean": mean_field("tree_truth_mean"),
        "gamma_smc_mean": mean_field("gamma_smc_mean"),
        "gamma_minus_truth_bias": mean_field("gamma_minus_truth_bias"),
        "mae": mean_field("mae"),
        "rmse": mean_field("rmse"),
        "pearson": np.nan,
        "spearman": np.nan,
        "truth_expected_sign_count": np.nan,
        "truth_expected_sign_rate": mean_field("truth_expected_sign_rate"),
        "truth_expected_sign_wilson95_low": np.nan,
        "truth_expected_sign_wilson95_high": np.nan,
        "gamma_expected_sign_count": np.nan,
        "gamma_expected_sign_rate": mean_field("gamma_expected_sign_rate"),
        "gamma_expected_sign_wilson95_low": np.nan,
        "gamma_expected_sign_wilson95_high": np.nan,
        "raw_truth_gamma_sign_agreement_count": np.nan,
        "raw_truth_gamma_sign_agreement_rate": mean_field(
            "raw_truth_gamma_sign_agreement_rate"
        ),
        "raw_truth_gamma_sign_agreement_wilson95_low": np.nan,
        "raw_truth_gamma_sign_agreement_wilson95_high": np.nan,
        "truth_expected_sign_denominator": np.nan,
        "gamma_detection_given_truth_expected_sign_count": np.nan,
        "gamma_detection_given_truth_expected_sign_rate": mean_field(
            "gamma_detection_given_truth_expected_sign_rate", detection_rows
        ),
        "gamma_detection_given_truth_expected_sign_wilson95_low": np.nan,
        "gamma_detection_given_truth_expected_sign_wilson95_high": np.nan,
        "interval_method": "not_computed_for_macro_cell_equal",
        "raw_sign_agreement_rule": "sign(tree_truth)==sign(gamma_smc)",
        "detection_rule": (
            "Equal-cell mean of within-cell Gamma expected-sign rates among cells "
            "whose tree-truth expected-sign denominator is nonzero"
        ),
        "pairing_rule": (
            "one_truth_and_one_gamma_value_per_replicate_then_equal_cell_macro"
        ),
    }


def _paired_accuracy(
    bases: Sequence[Mapping[str, Any]], paired: pd.DataFrame
) -> pd.DataFrame:
    cell_rows_by_metric: dict[str, list[dict[str, Any]]] = {
        metric: [] for metric in _METRICS
    }
    micro_rows: list[dict[str, Any]] = []
    for metric, expected_sign in _METRICS.items():
        metric_rows = paired[paired["metric"] == metric]
        micro_rows.append(
            _trajectory_accuracy_record(
                scope="overall_micro_trajectory_weighted",
                base=None,
                metric=metric,
                expected_sign=expected_sign,
                selected=metric_rows,
                aggregation_weighting="accepted_trajectory_equal_micro_secondary",
                is_headline=False,
                total_biological_cells=len(bases),
            )
        )
        for base in bases:
            base_id = str(base["specification_id"])
            selected = metric_rows[metric_rows["base_specification_id"] == base_id]
            cell_rows_by_metric[metric].append(
                _trajectory_accuracy_record(
                    scope="cell",
                    base=base,
                    metric=metric,
                    expected_sign=expected_sign,
                    selected=selected,
                    aggregation_weighting="trajectory_equal_within_biological_cell",
                    is_headline=False,
                    total_biological_cells=1,
                )
            )
    macro_rows = [
        _macro_cell_equal_accuracy(
            metric,
            cell_rows_by_metric[metric],
            total_biological_cells=len(bases),
        )
        for metric in _METRICS
    ]
    cell_rows = [row for metric in _METRICS for row in cell_rows_by_metric[metric]]
    return pd.DataFrame([*macro_rows, *micro_rows, *cell_rows])


def _aggregate_profiles(profiles: pd.DataFrame, runtime: StudyRuntime) -> pd.DataFrame:
    columns = [
        "base_specification_id",
        "selection_coefficient",
        "target_allele_frequency",
        "source",
        "genotype_class",
        "threshold_years",
        "n_observed_replicates",
        "mean_p_tmrca_lt_threshold",
        "median_p_tmrca_lt_threshold",
        "sample_sd_p_tmrca_lt_threshold",
        "q25_p_tmrca_lt_threshold",
        "q75_p_tmrca_lt_threshold",
        "bootstrap_mean_ci95_low",
        "bootstrap_mean_ci95_high",
        "bootstrap_draws",
        "bootstrap_seed",
        "conditioning_missingness",
    ]
    if profiles.empty:
        return pd.DataFrame(columns=columns)
    keys = [
        "base_specification_id",
        "selection_coefficient",
        "target_allele_frequency",
        "source",
        "genotype_class",
        "threshold_years",
    ]
    unit_keys = ["specification_id", *keys]
    if profiles.duplicated(unit_keys, keep=False).any():
        raise ValueError("profiles are not one threshold value per replicate/source")
    unique = profiles
    rows = []
    for key, group in unique.groupby(keys, dropna=False, sort=True):
        values = group["p_tmrca_lt_threshold"].to_numpy(dtype=float)
        values = values[np.isfinite(values)]
        seed = _stable_seed(
            runtime.base_seed,
            "replicate-profile-bootstrap:" + ":".join(map(str, key)),
            0,
        )
        low, high = _bootstrap_mean_interval(values, seed=seed)
        q25, q75 = (
            tuple(map(float, np.quantile(values, [0.25, 0.75])))
            if len(values)
            else (np.nan, np.nan)
        )
        rows.append(
            {
                **dict(zip(keys, key, strict=True)),
                "n_observed_replicates": int(len(values)),
                "mean_p_tmrca_lt_threshold": (
                    float(np.mean(values)) if len(values) else np.nan
                ),
                "median_p_tmrca_lt_threshold": (
                    float(np.median(values)) if len(values) else np.nan
                ),
                "sample_sd_p_tmrca_lt_threshold": (
                    float(np.std(values, ddof=1)) if len(values) > 1 else np.nan
                ),
                "q25_p_tmrca_lt_threshold": q25,
                "q75_p_tmrca_lt_threshold": q75,
                "bootstrap_mean_ci95_low": low,
                "bootstrap_mean_ci95_high": high,
                "bootstrap_draws": BOOTSTRAP_DRAWS,
                "bootstrap_seed": seed,
                "conditioning_missingness": "nonrandom",
            }
        )
    return pd.DataFrame(rows, columns=columns)


def _threshold_accuracy(profiles: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "summary_scope",
        "base_specification_id",
        "selection_coefficient",
        "target_allele_frequency",
        "genotype_class",
        "threshold_years",
        "n_paired_replicates",
        "gamma_minus_truth_bias",
        "mae",
        "rmse",
        "aggregation_weighting",
        "pairing_rule",
    ]
    if profiles.empty:
        return pd.DataFrame(columns=columns)
    ids = [
        "base_specification_id",
        "replicate_id",
        "replicate_index",
        "specification_id",
        "selection_coefficient",
        "target_allele_frequency",
        "genotype_class",
        "threshold_years",
    ]
    values = profiles[ids + ["source", "p_tmrca_lt_threshold"]]
    duplicates = values.duplicated(ids + ["source"], keep=False)
    if duplicates.any():
        raise ValueError("profiles are not one threshold value per replicate/source")
    paired = values.pivot(
        index=ids, columns="source", values="p_tmrca_lt_threshold"
    ).reset_index()
    if not set(_SOURCES).issubset(paired.columns):
        return pd.DataFrame(columns=columns)
    paired = paired.replace([np.inf, -np.inf], np.nan).dropna(subset=list(_SOURCES))
    paired["error"] = paired["gamma_smc"] - paired["tree_truth"]
    rows: list[dict[str, Any]] = []
    groupings = [
        (
            "overall_micro_trajectory_weighted",
            ["genotype_class", "threshold_years"],
        ),
        (
            "cell",
            [
                "base_specification_id",
                "selection_coefficient",
                "target_allele_frequency",
                "genotype_class",
                "threshold_years",
            ],
        ),
    ]
    for scope, keys in groupings:
        for key, group in paired.groupby(keys, dropna=False, sort=True):
            if not isinstance(key, tuple):
                key = (key,)
            labels = dict(zip(keys, key, strict=True))
            errors = group["error"].to_numpy(dtype=float)
            rows.append(
                {
                    "summary_scope": scope,
                    "base_specification_id": labels.get("base_specification_id"),
                    "selection_coefficient": labels.get(
                        "selection_coefficient", np.nan
                    ),
                    "target_allele_frequency": labels.get(
                        "target_allele_frequency", np.nan
                    ),
                    "genotype_class": labels["genotype_class"],
                    "threshold_years": labels["threshold_years"],
                    "n_paired_replicates": int(len(group)),
                    "gamma_minus_truth_bias": float(np.mean(errors)),
                    "mae": float(np.mean(np.abs(errors))),
                    "rmse": float(np.sqrt(np.mean(errors**2))),
                    "aggregation_weighting": (
                        "accepted_trajectory_equal_micro_across_biological_cells"
                        if scope == "overall_micro_trajectory_weighted"
                        else "trajectory_equal_within_biological_cell"
                    ),
                    "pairing_rule": "one_truth_and_one_gamma_value_per_replicate",
                }
            )
    return pd.DataFrame(rows, columns=columns)


def _report_float(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _report_int(value: Any) -> str:
    number = _report_float(value)
    return "—" if number is None else f"{int(round(number)):,}"


def _report_number(value: Any) -> str:
    number = _report_float(value)
    if number is None:
        return "—"
    if number == 0:
        return "0"
    absolute = abs(number)
    if absolute >= 1_000:
        return f"{number:,.1f}"
    if absolute >= 1:
        return f"{number:.3f}"
    return f"{number:.4f}"


def _report_rate(value: Any) -> str:
    number = _report_float(value)
    return "—" if number is None else f"{100.0 * number:.1f}%"


def _report_rate_with_interval(
    record: Mapping[str, Any],
    field: str,
    low_field: str,
    high_field: str,
) -> str:
    rate = _report_float(record.get(field))
    if rate is None:
        return "—"
    low = _report_float(record.get(low_field))
    high = _report_float(record.get(high_field))
    rendered = _report_rate(rate)
    if low is not None and high is not None:
        rendered += f" ({_report_rate(low)}–{_report_rate(high)})"
    return rendered


def _report_markdown_cell(value: Any) -> str:
    return str(value).replace("|", "\\|").replace("\n", " ")


def _append_report_table(
    lines: list[str], headers: Sequence[str], rows: Sequence[Sequence[Any]]
) -> None:
    lines.append("| " + " | ".join(map(_report_markdown_cell, headers)) + " |")
    lines.append("|" + "|".join("---" for _ in headers) + "|")
    for row in rows:
        lines.append("| " + " | ".join(map(_report_markdown_cell, row)) + " |")
    lines.append("")


def _report_frame(frame: pd.DataFrame, columns: Sequence[str]) -> pd.DataFrame:
    prepared = frame.copy()
    for column in columns:
        if column not in prepared:
            prepared[column] = np.nan
    return prepared


def _single_specification_number(
    specifications: Sequence[Mapping[str, Any]], key: str, label: str
) -> float:
    if not specifications:
        raise ValueError(f"RUN_RESULTS requires at least one specification for {label}")
    values = []
    for specification in specifications:
        value = _report_float(specification.get(key))
        if value is None:
            raise ValueError(
                f"RUN_RESULTS specification lacks finite {label}: "
                f"{specification.get('specification_id', '<unknown>')}"
            )
        values.append(value)
    first = values[0]
    if any(not _same_biological_value(value, first) for value in values[1:]):
        raise ValueError(f"RUN_RESULTS specifications disagree on {label}")
    return first


def _report_campaign_facts(
    runtime: StudyRuntime,
    specifications: Sequence[Mapping[str, Any]],
    campaign_design: Mapping[str, Any],
) -> dict[str, Any]:
    """Extract model-specific report statements from validated campaign inputs."""
    sample_diploids = _report_float(runtime.sample_diploids)
    if (
        sample_diploids is None
        or sample_diploids <= 0
        or not sample_diploids.is_integer()
    ):
        raise ValueError("RUN_RESULTS requires a positive integer sample_diploids")
    af_ceiling = _single_specification_number(
        specifications, "upper_allele_frequency", "AF ceiling"
    )
    if not 0 < af_ceiling <= 1:
        raise ValueError("RUN_RESULTS AF ceiling must be in (0, 1]")
    scaling = _report_float(runtime.slim_scaling_factor)
    burn_in = _report_float(runtime.slim_burn_in)
    if scaling is None or scaling <= 0:
        raise ValueError("RUN_RESULTS requires a positive SLiM scaling factor")
    if burn_in is None or burn_in < 0:
        raise ValueError("RUN_RESULTS requires a nonnegative SLiM burn-in")
    source_check = _single_specification_number(
        specifications,
        "slim_scaled_source_check_generations_ago",
        "scaled source-check time",
    )
    pulse = _single_specification_number(
        specifications,
        "slim_scaled_pulse_generations_ago",
        "scaled introgression-pulse time",
    )
    if source_check < 0 or pulse < 0:
        raise ValueError("RUN_RESULTS scaled event times must be nonnegative")
    fixed_model = campaign_design.get("fixed_model_contract")
    campaign_contract = campaign_design.get("campaign_contract")
    if not isinstance(fixed_model, Mapping) or not isinstance(
        campaign_contract, Mapping
    ):
        raise ValueError("RUN_RESULTS campaign design lacks fixed model contracts")
    acceptance_contract = campaign_contract.get("acceptance")
    archaic_contract = (
        acceptance_contract.get("archaic_specific_no_ils")
        if isinstance(acceptance_contract, Mapping)
        else None
    )
    if not isinstance(archaic_contract, Mapping):
        raise ValueError("RUN_RESULTS campaign design lacks archaic/no-ILS contract")
    mutation_mode = fixed_model.get("mutation_mode")
    origin_population = archaic_contract.get("mutation_declared_origin_population")
    after_split = archaic_contract.get("origin_after_human_neanderthal_split")
    before_pulse = archaic_contract.get("origin_before_introgression_pulse")
    entry_route = archaic_contract.get("han_entry_route")
    node_is_evidence = archaic_contract.get("tree_node_population_is_origin_evidence")
    origin_evidence = archaic_contract.get("origin_evidence")
    expected_archaic_contract = {
        "mutation_mode": "single archaic-specific de novo origin",
        "mutation_declared_origin_population": "Neanderthal",
        "origin_after_human_neanderthal_split": True,
        "origin_before_introgression_pulse": True,
        "han_entry_route": "introgression_pulse_only",
        "tree_node_population_is_origin_evidence": False,
        "origin_evidence": "serialized_forward_event_contract",
    }
    observed_archaic_contract = {
        "mutation_mode": mutation_mode,
        "mutation_declared_origin_population": origin_population,
        "origin_after_human_neanderthal_split": after_split,
        "origin_before_introgression_pulse": before_pulse,
        "han_entry_route": entry_route,
        "tree_node_population_is_origin_evidence": node_is_evidence,
        "origin_evidence": origin_evidence,
    }
    if observed_archaic_contract != expected_archaic_contract:
        raise ValueError(
            "RUN_RESULTS archaic/no-ILS contract is absent or incompatible: "
            + _canonical_json(observed_archaic_contract)
        )
    return {
        "sample_diploids": int(sample_diploids),
        "af_ceiling": af_ceiling,
        "slim_scaling_factor": scaling,
        "slim_burn_in": burn_in,
        "scaled_source_check_generations_ago": source_check,
        "scaled_introgression_pulse_generations_ago": pulse,
        "scaled_source_pulse_collision": _same_biological_value(source_check, pulse),
        **observed_archaic_contract,
    }


def _validated_report_campaign_facts(facts: Mapping[str, Any]) -> dict[str, Any]:
    required = {
        "sample_diploids",
        "af_ceiling",
        "slim_scaling_factor",
        "slim_burn_in",
        "scaled_source_check_generations_ago",
        "scaled_introgression_pulse_generations_ago",
        "scaled_source_pulse_collision",
        "mutation_mode",
        "mutation_declared_origin_population",
        "origin_after_human_neanderthal_split",
        "origin_before_introgression_pulse",
        "han_entry_route",
        "tree_node_population_is_origin_evidence",
        "origin_evidence",
    }
    missing = sorted(required - set(facts))
    if missing:
        raise ValueError(
            "RUN_RESULTS campaign facts are incomplete: " + ", ".join(missing)
        )
    numeric_keys = {
        "sample_diploids",
        "af_ceiling",
        "slim_scaling_factor",
        "slim_burn_in",
        "scaled_source_check_generations_ago",
        "scaled_introgression_pulse_generations_ago",
    }
    numeric = {key: _report_float(facts[key]) for key in numeric_keys}
    if any(value is None for value in numeric.values()):
        raise ValueError("RUN_RESULTS campaign facts contain nonfinite values")
    sample = numeric["sample_diploids"]
    af_ceiling = numeric["af_ceiling"]
    scaling = numeric["slim_scaling_factor"]
    burn_in = numeric["slim_burn_in"]
    source_check = numeric["scaled_source_check_generations_ago"]
    pulse = numeric["scaled_introgression_pulse_generations_ago"]
    if sample <= 0 or not sample.is_integer():
        raise ValueError("RUN_RESULTS sample_diploids must be a positive integer")
    if not 0 < af_ceiling <= 1 or scaling <= 0 or burn_in < 0:
        raise ValueError("RUN_RESULTS campaign facts have invalid model parameters")
    if source_check < 0 or pulse < 0:
        raise ValueError("RUN_RESULTS campaign facts have invalid event times")
    collision = facts["scaled_source_pulse_collision"]
    if not isinstance(collision, (bool, np.bool_)):
        raise ValueError("RUN_RESULTS source/pulse collision flag must be boolean")
    observed_collision = _same_biological_value(source_check, pulse)
    if bool(collision) != observed_collision:
        raise ValueError("RUN_RESULTS source/pulse collision flag disagrees with times")
    expected_archaic_contract = {
        "mutation_mode": "single archaic-specific de novo origin",
        "mutation_declared_origin_population": "Neanderthal",
        "origin_after_human_neanderthal_split": True,
        "origin_before_introgression_pulse": True,
        "han_entry_route": "introgression_pulse_only",
        "tree_node_population_is_origin_evidence": False,
        "origin_evidence": "serialized_forward_event_contract",
    }
    observed_archaic_contract = {key: facts[key] for key in expected_archaic_contract}
    if observed_archaic_contract != expected_archaic_contract:
        raise ValueError("RUN_RESULTS campaign facts violate archaic/no-ILS contract")
    return {
        **numeric,
        "sample_diploids": int(sample),
        "scaled_source_pulse_collision": observed_collision,
        **observed_archaic_contract,
    }


def _report_accuracy_rows(accuracy: pd.DataFrame, scope: str) -> list[dict[str, Any]]:
    prepared = _report_frame(
        accuracy,
        [
            "summary_scope",
            "metric",
            "n_paired_replicates",
            "n_biological_cells_meeting_minimum",
            "n_biological_cells_in_detection_rate",
        ],
    )
    selected = prepared[prepared["summary_scope"].astype(str) == scope].copy()
    paired = pd.to_numeric(selected["n_paired_replicates"], errors="coerce").fillna(0)
    if scope == "overall_macro_cell_equal":
        eligible = pd.to_numeric(
            selected["n_biological_cells_meeting_minimum"], errors="coerce"
        ).fillna(0)
        selected = selected[(paired > 0) & (eligible > 0)].copy()
    elif scope == "overall_micro_trajectory_weighted":
        selected = selected[paired > 0].copy()
    if selected.duplicated(["metric"], keep=False).any():
        raise ValueError(f"RUN_RESULTS has duplicate {scope!r} accuracy rows")
    metric_order = {metric: index for index, metric in enumerate(_METRICS)}
    selected["_metric_order"] = (
        selected["metric"].map(metric_order).fillna(len(metric_order))
    )
    selected = selected.sort_values(
        ["_metric_order", "metric"], kind="mergesort", na_position="last"
    )
    return selected.drop(columns="_metric_order").to_dict(orient="records")


def _render_run_results_markdown(
    summary: Mapping[str, Any],
    coverage: pd.DataFrame,
    acceptance: pd.DataFrame,
    accuracy: pd.DataFrame,
    effects: pd.DataFrame,
    threshold_accuracy: pd.DataFrame,
    campaign_facts: Mapping[str, Any],
) -> str:
    """Render a byte-stable campaign report from finalized aggregate values.

    The report deliberately contains no wall-clock generation timestamp.  Its
    bytes depend only on the supplied result objects, so repeated aggregation
    of the same validated inputs produces the same Markdown and checksum.
    """
    del effects  # Cell-stratified effects remain authoritative in their TSV.
    facts = _validated_report_campaign_facts(campaign_facts)
    coverage = _report_frame(
        coverage,
        [
            "base_specification_id",
            "selection_coefficient",
            "target_allele_frequency",
            "target_replicates",
            "attempted_replicates",
            "accepted_replicates",
            "decoded_replicates",
            "timed_out_replicates",
            "exhausted_replicates",
            "failed_replicates",
        ],
    )
    coverage_keys = [
        "base_specification_id"
        if coverage["base_specification_id"].notna().any()
        else "selection_coefficient",
        "target_allele_frequency",
    ]
    if not coverage.empty and coverage.duplicated(coverage_keys, keep=False).any():
        raise ValueError("RUN_RESULTS coverage contains duplicate biological cells")
    coverage = coverage.sort_values(
        ["selection_coefficient", "target_allele_frequency"],
        ascending=[False, True],
        kind="mergesort",
        na_position="last",
    ).reset_index(drop=True)
    acceptance = _report_frame(
        acceptance,
        [
            "specification_id",
            "achieved_allele_frequency",
            "n_hom_ref",
            "n_heterozygous",
            "n_hom_alt",
        ],
    ).sort_values(["specification_id"], kind="mergesort", na_position="last")

    requested_cells = summary.get("requested_biological_cells", len(coverage))
    requested_trajectories = summary.get(
        "requested_trajectories", coverage["target_replicates"].sum()
    )
    attempted_trajectories = summary.get(
        "attempted_trajectories", coverage["attempted_replicates"].sum()
    )
    accepted_trajectories = summary.get("accepted_trajectories", len(acceptance))
    decoded_trajectories = summary.get(
        "decoded_trajectories", coverage["decoded_replicates"].sum()
    )
    represented_cells = int(
        (pd.to_numeric(coverage["accepted_replicates"], errors="coerce") > 0).sum()
    )
    requested_cell_label = "cell" if _report_float(requested_cells) == 1 else "cells"
    represented_cell_label = "cell" if represented_cells == 1 else "cells"
    sample_diploids = int(facts["sample_diploids"])
    af_ceiling = float(facts["af_ceiling"])
    scaling = float(facts["slim_scaling_factor"])
    burn_in = float(facts["slim_burn_in"])
    source_check = float(facts["scaled_source_check_generations_ago"])
    pulse = float(facts["scaled_introgression_pulse_generations_ago"])
    if facts["scaled_source_pulse_collision"]:
        scaled_event_caveat = (
            f"the scaled source check and introgression pulse both occur at "
            f"{_report_number(source_check)} generations ago, so their requested "
            "separation is not retained"
        )
    else:
        scaled_event_caveat = (
            f"the scaled source check occurs at {_report_number(source_check)} "
            f"generations ago and the introgression pulse at "
            f"{_report_number(pulse)} generations ago"
        )

    lines = [
        "# Replicated EAS introgressed-sweep results",
        "",
        "## Bottom line",
        "",
        (
            f"This bounded campaign requested {_report_int(requested_trajectories)} "
            f"independently seeded trajectory slots across "
            f"{_report_int(requested_cells)} biological selection-by-AF-floor "
            f"{requested_cell_label}. {_report_int(attempted_trajectories)} slots "
            "were attempted, "
            f"{_report_int(accepted_trajectories)} produced accepted "
            f"{sample_diploids}-diploid Han panels, and "
            f"{_report_int(decoded_trajectories)} were decoded with Gamma-SMC. "
            f"{_report_int(represented_cells)} biological "
            f"{represented_cell_label} had at least one accepted trajectory."
        ),
        "",
    ]

    macro_rows = _report_accuracy_rows(accuracy, "overall_macro_cell_equal")
    recovery = []
    for record in macro_rows:
        rate = _report_float(
            record.get("gamma_detection_given_truth_expected_sign_rate")
        )
        cells = _report_float(record.get("n_biological_cells_in_detection_rate"))
        if rate is not None and cells is not None and cells > 0:
            metric = str(record.get("metric"))
            label = _REPORT_METRIC_LABELS.get(metric, metric)
            recovery.append(
                f"{label}: {_report_rate(rate)} across "
                f"{_report_int(cells)} contributing cells"
            )
    if recovery:
        lines.extend(
            [
                "The headline cell-equal conditional recovery estimates were "
                + "; ".join(recovery)
                + ". These are qualitative recovery checks in accepted, "
                "tree-truth-positive trajectories, not estimates of power.",
                "",
            ]
        )
    elif macro_rows:
        lines.extend(
            [
                "Eligible cells were available for the headline cell-equal "
                "comparison, but none contributed a tree-truth-positive "
                "trajectory to conditional recovery.",
                "",
            ]
        )
    else:
        lines.extend(
            [
                "No biological cell met the prespecified minimum for a headline "
                "cell-equal truth-versus-Gamma summary.",
                "",
            ]
        )

    lines.extend(["## Execution and coverage", ""])
    selection_values = sorted(
        {
            value
            for value in (
                _report_float(item)
                for item in coverage["selection_coefficient"].tolist()
            )
            if value is not None
        },
        reverse=True,
    )
    if selection_values:
        coverage_rows = []
        for selection in selection_values:
            group = coverage[
                np.isclose(
                    pd.to_numeric(
                        coverage["selection_coefficient"], errors="coerce"
                    ).to_numpy(dtype=float),
                    selection,
                    equal_nan=False,
                )
            ]
            accepted_values = pd.to_numeric(
                group["accepted_replicates"], errors="coerce"
            ).fillna(0)
            coverage_rows.append(
                [
                    f"{selection:g}",
                    _report_int(group["target_replicates"].sum()),
                    _report_int(group["attempted_replicates"].sum()),
                    _report_int(accepted_values.sum()),
                    _report_int(group["decoded_replicates"].sum()),
                    f"{int((accepted_values > 0).sum())}/{len(group)}",
                ]
            )
        _append_report_table(
            lines,
            [
                "Selection coefficient",
                "Requested slots",
                "Attempted",
                "Accepted",
                "Decoded",
                "Cells represented",
            ],
            coverage_rows,
        )

        af_values = sorted(
            {
                value
                for value in (
                    _report_float(item)
                    for item in coverage["target_allele_frequency"].tolist()
                )
                if value is not None
            }
        )
        lines.extend(
            [
                "### Accepted/requested slots by nested AF floor",
                "",
            ]
        )
        grid_rows = []
        for selection in selection_values:
            row = [f"s={selection:g}"]
            for af in af_values:
                selected = coverage[
                    np.isclose(
                        pd.to_numeric(
                            coverage["selection_coefficient"], errors="coerce"
                        ).to_numpy(dtype=float),
                        selection,
                        equal_nan=False,
                    )
                    & np.isclose(
                        pd.to_numeric(
                            coverage["target_allele_frequency"], errors="coerce"
                        ).to_numpy(dtype=float),
                        af,
                        equal_nan=False,
                    )
                ]
                if selected.empty:
                    row.append("—")
                else:
                    record = selected.iloc[0]
                    row.append(
                        f"{_report_int(record['accepted_replicates'])}/"
                        f"{_report_int(record['target_replicates'])} "
                        f"(d={_report_int(record['decoded_replicates'])})"
                    )
            grid_rows.append(row)
        _append_report_table(
            lines,
            ["Selection"] + [f"AF≥{100 * af:g}%" for af in af_values],
            grid_rows,
        )
    else:
        lines.extend(
            [
                "No biological-cell coverage rows were available at aggregation.",
                "",
            ]
        )

    if acceptance.empty:
        lines.extend(["No trajectories satisfied the acceptance contract.", ""])
    else:
        af = pd.to_numeric(
            acceptance["achieved_allele_frequency"], errors="coerce"
        ).dropna()
        lines.append(
            "Among accepted sample panels, achieved selected-allele AF was "
            + (
                f"{100 * af.min():.1f}%–{100 * af.max():.1f}% "
                f"(median {100 * af.median():.1f}%)."
                if not af.empty
                else "not available."
            )
        )
        genotype_fragments = []
        for column, label in (
            ("n_hom_ref", "ref/ref"),
            ("n_heterozygous", "heterozygous"),
            ("n_hom_alt", "alt/alt"),
        ):
            values = pd.to_numeric(acceptance[column], errors="coerce").dropna()
            if not values.empty:
                genotype_fragments.append(
                    f"{label} {_report_int(values.min())}–{_report_int(values.max())}"
                )
        if genotype_fragments:
            lines.append(
                "Observed per-panel genotype-count ranges were "
                + ", ".join(genotype_fragments)
                + "."
            )
        lines.append("")

    lines.extend(
        [
            f"The AF targets are nested acceptance floors `[f, {af_ceiling:g}]`, "
            "not disjoint bins. An achieved AF can satisfy every lower floor "
            "beneath it, although each cell has its own seed stream. Cross-floor "
            "differences therefore do not estimate a response across mutually "
            "exclusive AF strata.",
            "",
            "All effects are conditional on trajectories that returned within "
            "the bounded time/draw budgets and passed population AF, sample AF, "
            "and genotype-comparator checks. Failures are nonrandom. Intervals "
            "describe variation among accepted trajectories; they do not account "
            "for conditioning selection or missing cells. Tree-truth summaries "
            "use accepted trajectories, whereas Gamma and paired summaries "
            "additionally require a valid decode.",
            "",
            "## Truth-versus-Gamma recovery",
            "",
            "### Headline: biological-cell-equal macro summary",
            "",
        ]
    )
    if macro_rows:
        macro_table = []
        for record in macro_rows:
            metric = str(record.get("metric"))
            macro_table.append(
                [
                    _REPORT_METRIC_LABELS.get(metric, metric),
                    _report_int(record.get("n_biological_cells_meeting_minimum")),
                    _report_int(record.get("n_paired_replicates")),
                    _report_rate(record.get("truth_expected_sign_rate")),
                    _report_rate(record.get("gamma_expected_sign_rate")),
                    _report_rate(record.get("raw_truth_gamma_sign_agreement_rate")),
                    _report_rate(
                        record.get("gamma_detection_given_truth_expected_sign_rate")
                    ),
                    _report_number(record.get("gamma_minus_truth_bias")),
                    _report_number(record.get("mae")),
                    _report_number(record.get("rmse")),
                ]
            )
        _append_report_table(
            lines,
            [
                "Metric",
                "Eligible cells",
                "Paired trajectories",
                "Truth expected sign",
                "Gamma expected sign",
                "Raw sign agreement",
                "Gamma detection given truth",
                "Bias",
                "MAE",
                "RMSE",
            ],
            macro_table,
        )
    else:
        lines.extend(["No cell-equal macro estimate was available.", ""])
    lines.extend(
        [
            "The macro summary gives each eligible biological cell equal weight. "
            "Eligibility requires at least two paired trajectories; conditional "
            "detection also requires at least one truth-positive trajectory in "
            "the cell. No confidence interval is computed for this macro mean.",
            "",
            "### Secondary: accepted-trajectory-weighted micro summary",
            "",
        ]
    )
    micro_rows = _report_accuracy_rows(accuracy, "overall_micro_trajectory_weighted")
    if micro_rows:
        micro_table = []
        for record in micro_rows:
            metric = str(record.get("metric"))
            micro_table.append(
                [
                    _REPORT_METRIC_LABELS.get(metric, metric),
                    _report_int(record.get("n_paired_replicates")),
                    _report_rate_with_interval(
                        record,
                        "truth_expected_sign_rate",
                        "truth_expected_sign_wilson95_low",
                        "truth_expected_sign_wilson95_high",
                    ),
                    _report_rate_with_interval(
                        record,
                        "gamma_expected_sign_rate",
                        "gamma_expected_sign_wilson95_low",
                        "gamma_expected_sign_wilson95_high",
                    ),
                    _report_rate_with_interval(
                        record,
                        "raw_truth_gamma_sign_agreement_rate",
                        "raw_truth_gamma_sign_agreement_wilson95_low",
                        "raw_truth_gamma_sign_agreement_wilson95_high",
                    ),
                    _report_rate_with_interval(
                        record,
                        "gamma_detection_given_truth_expected_sign_rate",
                        "gamma_detection_given_truth_expected_sign_wilson95_low",
                        "gamma_detection_given_truth_expected_sign_wilson95_high",
                    ),
                    _report_number(record.get("gamma_minus_truth_bias")),
                    _report_number(record.get("mae")),
                    _report_number(record.get("rmse")),
                ]
            )
        _append_report_table(
            lines,
            [
                "Metric",
                "Paired trajectories",
                "Truth expected sign (95% CI)",
                "Gamma expected sign (95% CI)",
                "Raw sign agreement (95% CI)",
                "Gamma detection given truth (95% CI)",
                "Bias",
                "MAE",
                "RMSE",
            ],
            micro_table,
        )
    else:
        lines.extend(["No paired trajectory-level estimate was available.", ""])
    lines.extend(
        [
            "The micro summary gives each accepted paired trajectory equal "
            "weight and therefore overweights easier, higher-coverage cells. It "
            "is a different estimand from the headline cell-equal macro summary.",
            "",
            "Raw sign agreement can count a jointly contrary truth/Gamma pair as "
            "agreement. Conditional detection does not: it includes only "
            "trajectories with the expected tree-truth sign and asks whether "
            "Gamma-SMC also has that sign.",
            "",
            "No matched neutral/null trajectories were simulated. These rates "
            "are descriptive checks, not p-values, sensitivity, specificity, "
            "false-positive rate, power, or ROC AUC. The normalized integrated "
            "CDF area is an area under a TMRCA CDF, not classifier AUC.",
            "",
            "## TMRCA probability-profile accuracy: overall micro summary",
            "",
            "This threshold-error summary pools accepted, decoded paired "
            "trajectories across biological cells with equal weight per "
            "trajectory. It is therefore a secondary overall "
            "micro/trajectory-weighted estimand: easier, higher-coverage cells "
            "contribute more heavily and can dominate the pooled error.",
            "",
        ]
    )
    thresholds = _report_frame(
        threshold_accuracy,
        [
            "summary_scope",
            "genotype_class",
            "threshold_years",
            "n_paired_replicates",
            "gamma_minus_truth_bias",
            "mae",
            "rmse",
        ],
    )
    thresholds = thresholds[
        thresholds["summary_scope"].astype(str) == "overall_micro_trajectory_weighted"
    ].copy()
    threshold_pairs = pd.to_numeric(
        thresholds["n_paired_replicates"], errors="coerce"
    ).fillna(0)
    thresholds = thresholds[threshold_pairs > 0].copy()
    threshold_values = pd.to_numeric(thresholds["threshold_years"], errors="coerce")
    if not thresholds.empty and threshold_values.notna().any():
        maximum = float(threshold_values.max())
        selected = thresholds[np.isclose(threshold_values, maximum)].copy()
        class_order = {name: index for index, name in enumerate(_PROFILE_CLASSES)}
        selected["_class_order"] = (
            selected["genotype_class"].map(class_order).fillna(len(class_order))
        )
        selected = selected.sort_values(
            ["_class_order", "genotype_class"], kind="mergesort"
        )
        threshold_rows = []
        for record in selected.to_dict(orient="records"):
            genotype_class = str(record.get("genotype_class"))
            threshold_rows.append(
                [
                    _REPORT_CLASS_LABELS.get(genotype_class, genotype_class),
                    _report_int(record.get("n_paired_replicates")),
                    _report_number(record.get("gamma_minus_truth_bias")),
                    _report_number(record.get("mae")),
                    _report_number(record.get("rmse")),
                ]
            )
        lines.append(f"At the largest decoded cutoff ({maximum / 1_000:g} kya):")
        lines.append("")
        _append_report_table(
            lines,
            ["Genotype class", "Paired", "Bias", "MAE", "RMSE"],
            threshold_rows,
        )
    else:
        lines.extend(["No paired truth/Gamma TMRCA-threshold rows were available.", ""])
    lines.extend(
        [
            "Cell-stratified overall, ref/ref, and alt/alt curves and their "
            "trajectory-bootstrap bands are in "
            "`results/aggregate_cdf_profiles.tsv` and the PNG/PDF figure set. "
            "They are not pooled into an AF-response curve.",
            "",
            "## Interpretation and provenance limits",
            "",
            f"- These results characterize the configured `Q={scaling:g}`, "
            f"`slim_burn_in={burn_in:g}` process. Q rescaling can materially "
            "reduce small simulated populations and collapse short demographic "
            f"phases; in this campaign, {scaled_event_caveat}. Q/burn-in "
            "sensitivity runs are required before biological calibration.",
            "- Archaic specificity/no ILS is bound to the validated "
            f"`{facts['origin_evidence']}`: mutation mode "
            f"`{facts['mutation_mode']}`, declared origin population "
            f"`{facts['mutation_declared_origin_population']}`, origin after the "
            "human-Neanderthal split and before the pulse, and Han entry route "
            f"`{facts['han_entry_route']}`. The population attached to the "
            "retained mutation node is contractually not origin evidence; donor "
            "AF, migrant-copy count, and the full trajectory are not recovered "
            "from the returned tree.",
            "- Artifact validity does not establish semantic equivalence across "
            "implementation source epochs. Any format-only transition that is "
            "not cryptographically bound per task remains an operator "
            "attestation; a frozen-source rerun is required to remove that "
            "provenance limitation.",
            "",
            "## Machine-readable outputs",
            "",
            "- `results/cell_coverage.tsv`: requested and observed coverage for "
            "all biological cells.",
            "- `results/replicate_acceptance.tsv`: achieved sample AF, genotype "
            "counts, seeds, and effort for accepted panels.",
            "- `results/cell_effect_summary.tsv`: exact cell/source focal-effect "
            "summaries and trajectory-bootstrap intervals.",
            "- `results/truth_gamma_accuracy.tsv`: explicitly labeled macro and "
            "micro recovery estimands.",
            "- `results/aggregate_cdf_profiles.tsv` and "
            "`results/threshold_accuracy.tsv`: TMRCA-CDF profiles and paired "
            "truth/Gamma errors; the overall threshold rows are explicitly "
            "trajectory-weighted micro summaries.",
            "- `results/replicate_results_manifest.json`: artifact-relative "
            "paths, SHA-256 checksums, contracts, decoder hashes, and environment "
            "provenance. This `RUN_RESULTS.md` file is itself checksum-covered by "
            "that manifest.",
            "",
        ]
    )
    return "\n".join(lines)


def _atomic_text(path: Path, value: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f"{path.name}.tmp.{os.getpid()}.{uuid.uuid4().hex}")
    try:
        with temporary.open("w", encoding="utf-8", newline="\n") as handle:
            handle.write(value)
        _replace_with_access_retry(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)
    return path


def _save_figure(fig, stem: Path) -> dict[str, Path]:
    stem.parent.mkdir(parents=True, exist_ok=True)
    outputs: dict[str, Path] = {}
    try:
        for suffix, kwargs in (
            ("png", {"dpi": 220}),
            ("pdf", {}),
        ):
            output = stem.with_suffix(f".{suffix}")
            temporary = output.with_name(
                f"{output.stem}.tmp.{os.getpid()}{output.suffix}"
            )
            try:
                # Do not use bbox_inches="tight": preserving the declared
                # 11 x 8.5 inch canvas is part of the publication contract.
                fig.savefig(temporary, format=suffix, **kwargs)
                _replace_with_access_retry(temporary, output)
            finally:
                temporary.unlink(missing_ok=True)
            outputs[suffix] = output
    finally:
        plt.close(fig)
    return outputs


def _heatmap(
    ax,
    frame: pd.DataFrame,
    value: str,
    text,
    *,
    title: str,
    cmap: str = "Blues",
    vmin: float | None = 0,
    vmax: float | None = 1,
) -> None:
    selections = sorted(frame["selection_coefficient"].unique(), reverse=True)
    frequencies = sorted(frame["target_allele_frequency"].unique())
    matrix = np.full((len(selections), len(frequencies)), np.nan)
    labels: dict[tuple[int, int], str] = {}
    for _, row in frame.iterrows():
        i = selections.index(row["selection_coefficient"])
        j = frequencies.index(row["target_allele_frequency"])
        matrix[i, j] = row[value]
        labels[(i, j)] = text(row)
    masked = np.ma.masked_invalid(matrix)
    image = ax.imshow(masked, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xticks(range(len(frequencies)), [f"{x:.0%}" for x in frequencies])
    ax.set_yticks(range(len(selections)), [f"{x:g}" for x in selections])
    ax.set_xlabel("Required final Han AF floor")
    ax.set_ylabel("Selection coefficient")
    ax.set_title(title)
    for (i, j), label in labels.items():
        ax.text(j, i, label, ha="center", va="center", fontsize=_FONT["annotation"])
    return image


def _plot_coverage(coverage: pd.DataFrame, stem: Path) -> dict[str, Path]:
    with plt.rc_context({"font.size": _FONT["tick"]}):
        fig, axes = plt.subplots(
            1, 2, figsize=_LETTER_LANDSCAPE, constrained_layout=True
        )
        accepted_image = _heatmap(
            axes[0],
            coverage,
            "accepted_fraction_of_target",
            lambda row: (
                f"{int(row['accepted_replicates'])}/{int(row['target_replicates'])}"
            ),
            title="Accepted trajectories / requested",
        )
        decoded_image = _heatmap(
            axes[1],
            coverage,
            "decoded_fraction_of_accepted",
            lambda row: (
                f"{int(row['decoded_replicates'])}/{int(row['accepted_replicates'])}"
            ),
            title="Decoded / accepted trajectories",
        )
        fig.colorbar(accepted_image, ax=axes[0], label="Fraction")
        fig.colorbar(decoded_image, ax=axes[1], label="Fraction")
        fig.suptitle(
            "Replicate availability (conditioning failures are nonrandom)",
            fontsize=_FONT["title"],
        )
    return _save_figure(fig, stem)


def _plot_achieved(acceptance: pd.DataFrame, stem: Path) -> dict[str, Path]:
    with plt.rc_context({"font.size": _FONT["tick"]}):
        fig, axes = plt.subplots(
            1, 2, figsize=_LETTER_LANDSCAPE, constrained_layout=True
        )
        if acceptance.empty:
            for ax in axes:
                ax.text(0.5, 0.5, "No accepted trajectories", ha="center", va="center")
                ax.set_axis_off()
        else:
            for selection, group in acceptance.groupby("selection_coefficient"):
                color = _COLORS.get(float(selection), "#666666")
                jitter = (
                    (group["replicate_index"].to_numpy(dtype=float) % 7) - 3
                ) * 0.002
                axes[0].scatter(
                    group["target_allele_frequency"] + jitter,
                    group["achieved_allele_frequency"],
                    s=55,
                    alpha=0.8,
                    color=color,
                    label=f"s={selection:g}",
                )
                axes[1].scatter(
                    group["achieved_allele_frequency"],
                    group["n_hom_ref"],
                    s=48,
                    marker="o",
                    color=color,
                    alpha=0.7,
                )
                axes[1].scatter(
                    group["achieved_allele_frequency"],
                    group["n_heterozygous"],
                    s=48,
                    marker="^",
                    color=color,
                    alpha=0.7,
                )
                axes[1].scatter(
                    group["achieved_allele_frequency"],
                    group["n_hom_alt"],
                    s=48,
                    marker="s",
                    color=color,
                    alpha=0.7,
                )
            axes[0].plot([0, 1], [0, 1], "k--", linewidth=1.5)
            axes[0].set(xlim=(0, 1), ylim=(0, 1))
            axes[0].set_xlabel("Required Han AF floor")
            axes[0].set_ylabel("Achieved sample AF")
            axes[0].set_title("Frequency-conditioned accepted draws")
            axes[0].legend(fontsize=_FONT["legend"])
            axes[1].set_xlabel("Achieved sample AF")
            axes[1].set_ylabel("Diploid count")
            axes[1].set_title(
                "Genotype counts (circle ref/ref, triangle ref/alt, square alt/alt)"
            )
            for ax in axes:
                ax.grid(alpha=0.25)
        fig.suptitle(
            "Accepted introgressed-sweep replicate panels", fontsize=_FONT["title"]
        )
    return _save_figure(fig, stem)


def _plot_truth_vs_gamma(paired: pd.DataFrame, stem: Path) -> dict[str, Path]:
    with plt.rc_context({"font.size": _FONT["tick"]}):
        fig, axes = plt.subplots(
            1, 2, figsize=_LETTER_LANDSCAPE, constrained_layout=True
        )
        for ax, (metric, expected) in zip(axes, _METRICS.items(), strict=True):
            selected = paired[paired["metric"] == metric]
            if selected.empty:
                ax.text(
                    0.5, 0.5, "No paired decoded replicates", ha="center", va="center"
                )
                ax.set_axis_off()
                continue
            for selection, group in selected.groupby("selection_coefficient"):
                ax.scatter(
                    group["tree_truth"],
                    group["gamma_smc"],
                    s=65,
                    alpha=0.8,
                    color=_COLORS.get(float(selection), "#666666"),
                    label=f"s={selection:g}",
                )
            extrema = np.r_[selected["tree_truth"], selected["gamma_smc"]]
            low, high = float(np.min(extrema)), float(np.max(extrema))
            padding = max((high - low) * 0.08, 1e-12)
            ax.plot(
                [low - padding, high + padding], [low - padding, high + padding], "k--"
            )
            ax.set_xlabel("Tree truth")
            ax.set_ylabel("Gamma-SMC")
            title = (
                "Alt-ref normalized CDF area"
                if expected == "positive"
                else "Alt-ref mean TMRCA (generations)"
            )
            ax.set_title(title)
            ax.grid(alpha=0.25)
            ax.legend(fontsize=_FONT["legend"])
        fig.suptitle(
            "Per-trajectory truth versus Gamma-SMC effects", fontsize=_FONT["title"]
        )
    return _save_figure(fig, stem)


def _plot_effect_sign(summary: pd.DataFrame, stem: Path) -> dict[str, Path]:
    with plt.rc_context({"font.size": _FONT["tick"]}):
        fig, axes = plt.subplots(
            2, 2, figsize=_LETTER_LANDSCAPE, constrained_layout=True
        )
        for row_index, source in enumerate(_SOURCES):
            for column_index, metric in enumerate(_METRICS):
                ax = axes[row_index, column_index]
                selected = summary[
                    (summary["source"] == source) & (summary["metric"] == metric)
                ].copy()
                values = selected["mean"].to_numpy(dtype=float)
                finite = np.abs(values[np.isfinite(values)])
                bound = float(finite.max()) if len(finite) else 1.0
                if not np.isfinite(bound) or bound <= 0:
                    bound = 1.0
                title_metric = (
                    "Delta normalized CDF area"
                    if "cdf" in metric
                    else "Delta mean TMRCA (generations)"
                )
                image = _heatmap(
                    ax,
                    selected,
                    "mean",
                    lambda record: (
                        "NA\n0/0"
                        if pd.isna(record["mean"])
                        else f"{record['mean']:.3g}\n{int(record['expected_sign_count'])}/{int(record['n_observed_replicates'])} sign"
                    ),
                    title=f"{source}: {title_metric}",
                    cmap="RdBu_r",
                    vmin=-bound,
                    vmax=bound,
                )
                fig.colorbar(image, ax=ax, shrink=0.8, label="Mean alt - ref effect")
        fig.suptitle(
            "Observed replicate effects and expected-sign counts",
            fontsize=_FONT["title"],
        )
    return _save_figure(fig, stem)


def _selection_slug(value: float) -> str:
    return f"{float(value):.3f}".replace(".", "p")


def _plot_cdf_profiles(
    aggregated: pd.DataFrame,
    bases: Sequence[Mapping[str, Any]],
    stem: Path,
) -> dict[str, Path]:
    """Plot exact s-by-AF cells on readable letter-landscape pages.

    A page covers one selection coefficient and at most three AF floors. Rows
    are overall, ref/ref, and alt/alt. Tree-truth and Gamma-SMC means each have
    a trajectory-bootstrap 95% confidence band. Empty panels are retained
    because conditioning failures are nonrandom.
    """
    cells = pd.DataFrame(
        [
            {
                "base_specification_id": str(base["specification_id"]),
                "selection_coefficient": float(base["selection_coefficient"]),
                "target_allele_frequency": float(base["target_allele_frequency"]),
            }
            for base in bases
        ]
    ).drop_duplicates()
    outputs: dict[str, Path] = {}
    source_styles = {
        "tree_truth": ("#1B1B1B", "-", "Tree truth"),
        "gamma_smc": ("#0072B2", "--", "Gamma-SMC"),
    }
    class_labels = {
        "overall": "Overall",
        "hom_ref": "Ref/ref",
        "hom_alt": "Alt/alt",
    }
    for selection in sorted(cells["selection_coefficient"].unique(), reverse=True):
        selection_cells = cells[
            np.isclose(cells["selection_coefficient"], selection)
        ].sort_values("target_allele_frequency")
        for chunk_start in range(0, len(selection_cells), 3):
            chunk = selection_cells.iloc[chunk_start : chunk_start + 3]
            with plt.rc_context({"font.size": 11}):
                fig, axes = plt.subplots(
                    len(_PROFILE_CLASSES),
                    len(chunk),
                    figsize=_LETTER_LANDSCAPE,
                    sharex=True,
                    sharey=True,
                    squeeze=False,
                    constrained_layout=True,
                )
                for column, (_, cell) in enumerate(chunk.iterrows()):
                    cell_profile = (
                        aggregated[
                            (
                                aggregated["base_specification_id"]
                                == cell["base_specification_id"]
                            )
                            & np.isclose(aggregated["selection_coefficient"], selection)
                            & np.isclose(
                                aggregated["target_allele_frequency"],
                                cell["target_allele_frequency"],
                            )
                        ]
                        if not aggregated.empty
                        else pd.DataFrame()
                    )
                    for row, genotype_class in enumerate(_PROFILE_CLASSES):
                        ax = axes[row, column]
                        plotted = False
                        for source in _SOURCES:
                            group = (
                                cell_profile[
                                    (cell_profile["genotype_class"] == genotype_class)
                                    & (cell_profile["source"] == source)
                                ].sort_values("threshold_years")
                                if not cell_profile.empty
                                else pd.DataFrame()
                            )
                            if group.empty:
                                continue
                            x = group["threshold_years"].to_numpy(dtype=float) / 1000
                            mean = group["mean_p_tmrca_lt_threshold"].to_numpy(
                                dtype=float
                            )
                            low = group["bootstrap_mean_ci95_low"].to_numpy(dtype=float)
                            high = group["bootstrap_mean_ci95_high"].to_numpy(
                                dtype=float
                            )
                            finite_line = np.isfinite(x) & np.isfinite(mean)
                            if not finite_line.any():
                                continue
                            color, linestyle, label = source_styles[source]
                            ax.plot(
                                x[finite_line],
                                mean[finite_line],
                                color=color,
                                linestyle=linestyle,
                                marker="o",
                                markersize=3.5,
                                linewidth=2,
                                label=label,
                            )
                            finite_band = (
                                np.isfinite(x) & np.isfinite(low) & np.isfinite(high)
                            )
                            if finite_band.any():
                                ax.fill_between(
                                    x[finite_band],
                                    low[finite_band],
                                    high[finite_band],
                                    color=color,
                                    alpha=0.14,
                                )
                            plotted = True
                        if not plotted:
                            ax.text(
                                0.5,
                                0.5,
                                "No observed\ntrajectories",
                                ha="center",
                                va="center",
                                transform=ax.transAxes,
                                fontsize=10,
                            )
                        if row == 0:
                            ax.set_title(
                                f"Han AF floor {cell['target_allele_frequency']:.0%}",
                                fontsize=14,
                            )
                        if column == 0:
                            ax.set_ylabel(
                                f"{class_labels[genotype_class]}\nP(TMRCA < x)",
                                fontsize=12,
                            )
                        if row == len(_PROFILE_CLASSES) - 1:
                            ax.set_xlabel("Threshold (kya)", fontsize=12)
                        ax.set_ylim(-0.03, 1.03)
                        ax.grid(alpha=0.22)
                        ax.tick_params(labelsize=10)
                handles = [
                    plt.Line2D(
                        [],
                        [],
                        color=color,
                        linestyle=linestyle,
                        marker="o",
                        label=label,
                    )
                    for color, linestyle, label in source_styles.values()
                ]
                fig.legend(
                    handles=handles,
                    loc="lower center",
                    ncol=2,
                    fontsize=11,
                    bbox_to_anchor=(0.5, -0.015),
                )
                fig.suptitle(
                    f"Focal-site TMRCA CDFs: s={selection:g} "
                    "(trajectory bootstrap 95% CI)",
                    fontsize=17,
                )
            low_af = int(round(float(chunk.iloc[0]["target_allele_frequency"]) * 100))
            high_af = int(round(float(chunk.iloc[-1]["target_allele_frequency"]) * 100))
            page_key = f"s{_selection_slug(selection)}_af{low_af:02d}_{high_af:02d}"
            page_paths = _save_figure(fig, stem.with_name(f"{stem.name}_{page_key}"))
            for suffix, path in page_paths.items():
                outputs[f"{page_key}_{suffix}"] = path
    return outputs


def _plot_threshold_mae(summary: pd.DataFrame, stem: Path) -> dict[str, Path]:
    overall = (
        summary[summary["summary_scope"] == "overall_micro_trajectory_weighted"]
        if not summary.empty
        else pd.DataFrame()
    )
    with plt.rc_context({"font.size": _FONT["tick"]}):
        fig, ax = plt.subplots(figsize=_LETTER_LANDSCAPE, constrained_layout=True)
        colors = {
            "overall": "#000000",
            "hom_ref": "#666666",
            "heterozygous": "#E69F00",
            "hom_alt": "#0072B2",
        }
        for genotype_class, group in (
            overall.groupby("genotype_class") if not overall.empty else []
        ):
            group = group.sort_values("threshold_years")
            ax.plot(
                group["threshold_years"] / 1000,
                group["mae"],
                marker="o",
                linewidth=2.5,
                color=colors.get(str(genotype_class), "#666666"),
                label=str(genotype_class).replace("_", "/"),
            )
        if overall.empty:
            ax.text(
                0.5,
                0.5,
                "No paired decoded profiles",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
        ax.set_xlabel("TMRCA threshold (kya)")
        ax.set_ylabel("Mean absolute error in P(TMRCA < threshold)")
        ax.set_title(
            "Gamma-SMC threshold accuracy\n"
            "Overall micro (trajectory-weighted across accepted pairs)"
        )
        ax.grid(alpha=0.25)
        if not overall.empty:
            ax.legend(fontsize=_FONT["legend"])
    return _save_figure(fig, stem)


def _json_clean(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_clean(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_json_clean(item) for item in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return None if not np.isfinite(value) else float(value)
    return value


def _portable_specification_records(
    specifications: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    records = [
        {
            key: value
            for key, value in specification.items()
            if key not in {"model_record", "origin_record"}
        }
        for specification in specifications
    ]
    return sorted(records, key=lambda record: str(record["specification_id"]))


def _file_provenance(path: Path, campaign_dir: Path) -> dict[str, Any]:
    record: dict[str, Any] = {
        "path": path.relative_to(campaign_dir).as_posix(),
        "present": path.is_file(),
        "sha256": None,
    }
    if path.is_file():
        record["sha256"] = _sha256(path)
    return record


def _portable_runtime_record(
    runtime: StudyRuntime, campaign_dir: Path
) -> dict[str, Any]:
    record = asdict(runtime)
    relative_repo = os.path.relpath(
        Path(runtime.repo_root).resolve(), campaign_dir.resolve()
    ).replace("\\", "/")
    record["repo_root"] = relative_repo
    record["repo_root_path_semantics"] = "relative_to_campaign_directory"
    return record


def _command_stdout(command: Sequence[str], cwd: Path) -> str | None:
    try:
        completed = subprocess.run(
            list(command),
            cwd=cwd,
            text=True,
            capture_output=True,
            check=False,
            timeout=30,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return completed.stdout.strip() if completed.returncode == 0 else None


def _environment_provenance(repo_root: Path) -> dict[str, Any]:
    repository = repo_root.resolve()
    head = _command_stdout(["git", "rev-parse", "HEAD"], repository)
    branch = _command_stdout(["git", "rev-parse", "--abbrev-ref", "HEAD"], repository)
    dirty_patch = _command_stdout(
        ["git", "diff", "--binary", "HEAD", "--", "."], repository
    )
    untracked = _command_stdout(
        ["git", "ls-files", "--others", "--exclude-standard"], repository
    )
    tracked_patch_sha256 = (
        hashlib.sha256(dirty_patch.encode("utf-8")).hexdigest()
        if dirty_patch is not None
        else None
    )
    lockfiles = {}
    for filename in ("pyproject.toml", "uv.lock"):
        path = repository / filename
        lockfiles[filename] = {
            "path": filename,
            "present": path.is_file(),
            "sha256": _sha256(path) if path.is_file() else None,
        }
    return {
        "git": {
            "head_sha": head,
            "branch": branch,
            "tracked_dirty_patch_sha256": tracked_patch_sha256,
            "tracked_dirty_patch_present": bool(dirty_patch),
            "untracked_path_count": (
                len(untracked.splitlines()) if untracked is not None else None
            ),
            "untracked_content_hashed": False,
            "provenance_status": (
                "available"
                if head is not None and branch is not None
                else "unavailable"
            ),
        },
        "python": {
            "version": sys.version,
            "implementation": sys.implementation.name,
            "executable_basename": Path(sys.executable).name,
        },
        "dependency_manifests": lockfiles,
    }


def _campaign_input_provenance(
    campaign_dir: Path,
    bases: Sequence[Mapping[str, Any]],
    expanded: Sequence[Mapping[str, Any]],
    runtime: StudyRuntime,
    snapshot: Mapping[str, Any],
) -> dict[str, Any]:
    portable_bases = _portable_specification_records(bases)
    portable_expanded = _portable_specification_records(expanded)
    snapshot_root = Path(snapshot["root"])
    source_path = Path(__file__).resolve()
    repo_root = Path(runtime.repo_root).resolve()
    try:
        portable_source_path = source_path.relative_to(repo_root).as_posix()
    except ValueError:
        portable_source_path = source_path.name
    return {
        "study_design_json": _file_provenance(
            snapshot_root / "study_design.json", campaign_dir
        ),
        "specifications_tsv": _file_provenance(
            snapshot_root / "specifications.tsv", campaign_dir
        ),
        "input_snapshot": {
            "snapshot_id": snapshot["snapshot_id"],
            "manifest": _file_provenance(Path(snapshot["manifest_path"]), campaign_dir),
            "file_count": len(snapshot["manifest"]["files"]),
            "files": snapshot["manifest"]["files"],
        },
        "base_specifications": {
            "count": len(portable_bases),
            "canonical_sha256": _contract_sha256(portable_bases),
        },
        "expanded_specifications": {
            "count": len(portable_expanded),
            "canonical_sha256": _contract_sha256(portable_expanded),
        },
        "analysis_source": {
            "path": portable_source_path,
            "present": True,
            "sha256": _sha256(source_path),
        },
    }


def _observed_decoder_sha256s(
    source_study_dir: Path,
    expanded: Sequence[Mapping[str, Any]],
    explicit_decoder_sha256: str | None,
) -> list[str]:
    """Fail closed on mixed decoders or an explicit/recorded SHA mismatch."""
    observed: set[str] = set()
    decoded_bundles = 0
    for specification in expanded:
        specification_dir = (
            source_study_dir / "work" / str(specification["specification_id"])
        )
        if not (specification_dir / "decode_complete.json").is_file():
            continue
        decoded_bundles += 1
        bundle_shas: set[str] = set()
        for contract_path in sorted(
            (specification_dir / "decoded").glob("*.contract.json")
        ):
            try:
                payload = json.loads(contract_path.read_text(encoding="utf-8"))
                decoder_sha256 = str(payload["decoder_sha256"])
            except (
                KeyError,
                OSError,
                TypeError,
                ValueError,
                json.JSONDecodeError,
            ) as error:
                raise RuntimeError(
                    f"decoded class contract lacks decoder SHA: {contract_path}"
                ) from error
            valid_sha256 = len(decoder_sha256) == 64 and all(
                character in "0123456789abcdefABCDEF" for character in decoder_sha256
            )
            if not valid_sha256:
                raise RuntimeError(
                    f"decoded class contract has malformed decoder SHA: {contract_path}"
                )
            bundle_shas.add(decoder_sha256)
        if len(bundle_shas) != 1:
            raise RuntimeError(
                "one decoded bundle must contain exactly one Gamma-SMC binary SHA: "
                f"{specification['specification_id']} observed={sorted(bundle_shas)}"
            )
        observed.update(bundle_shas)
    if decoded_bundles and len(observed) != 1:
        raise RuntimeError(
            "mixed Gamma-SMC binary SHA sets require explicit stratified analysis; "
            f"observed={sorted(observed)}"
        )
    if (
        explicit_decoder_sha256 is not None
        and observed
        and observed != {explicit_decoder_sha256}
    ):
        raise RuntimeError(
            "explicit Gamma-SMC decoder SHA does not match recorded decoded outputs: "
            f"explicit={explicit_decoder_sha256} observed={sorted(observed)}"
        )
    return sorted(observed)


def _completion_provenance(
    campaign_dir: Path,
    expanded: Sequence[Mapping[str, Any]],
    status: pd.DataFrame,
    explicit_decoder_sha256: str | None,
    *,
    source_study_dir: Path | None = None,
    observed_decoder_sha256s: Sequence[str] = (),
) -> dict[str, Any]:
    source_study_dir = source_study_dir or campaign_dir
    status_lookup = status.set_index("specification_id").to_dict(orient="index")
    simulation_completions: list[dict[str, Any]] = []
    simulation_contracts: list[dict[str, Any]] = []
    decode_completions: list[dict[str, Any]] = []
    decode_contracts: list[dict[str, Any]] = []
    class_decoder_contracts: list[dict[str, Any]] = []
    contract_decoder_sha256s: set[str] = set()
    for specification in expanded:
        specification_id = str(specification["specification_id"])
        row = status_lookup.get(specification_id, {})
        specification_dir = source_study_dir / "work" / specification_id
        if row.get("simulation_status") == "complete":
            for collection, filename in (
                (simulation_contracts, "simulation_contract.json"),
                (simulation_completions, "simulation_complete.json"),
            ):
                record = _file_provenance(specification_dir / filename, campaign_dir)
                record["specification_id"] = specification_id
                collection.append(record)
        if row.get("decode_status") != "complete":
            continue
        for collection, filename in (
            (decode_contracts, "decode_contract.json"),
            (decode_completions, "decode_complete.json"),
        ):
            record = _file_provenance(specification_dir / filename, campaign_dir)
            record["specification_id"] = specification_id
            collection.append(record)
        for contract_path in sorted(
            (specification_dir / "decoded").glob("*.contract.json")
        ):
            record = _file_provenance(contract_path, campaign_dir)
            record["specification_id"] = specification_id
            try:
                payload = json.loads(contract_path.read_text(encoding="utf-8"))
                decoder_sha256 = str(payload["decoder_sha256"])
                record["decoder_sha256"] = decoder_sha256
                contract_decoder_sha256s.add(decoder_sha256)
            except (KeyError, OSError, TypeError, ValueError, json.JSONDecodeError):
                record["decoder_sha256"] = None
            class_decoder_contracts.append(record)
    if not contract_decoder_sha256s.issubset(set(observed_decoder_sha256s)):
        raise RuntimeError(
            "validated decode provenance disagrees with the decoder SHA snapshot"
        )
    return {
        "validated_simulation_contracts": simulation_contracts,
        "validated_simulation_completions": simulation_completions,
        "validated_decode_contracts": decode_contracts,
        "validated_decode_completions": decode_completions,
        "validated_class_decoder_contracts": class_decoder_contracts,
        "observed_decoder_sha256s": sorted(observed_decoder_sha256s),
        "validated_decoder_sha256s": sorted(contract_decoder_sha256s),
        "explicit_decoder_sha256": explicit_decoder_sha256,
        "decoder_sha256s": sorted(observed_decoder_sha256s),
    }


def _aggregate_introgression_replicates_locked(
    runtime: StudyRuntime,
    base_specs: Sequence[Mapping[str, Any]],
    expanded_specs: Sequence[Mapping[str, Any]],
    study_dir: str | Path,
    decoder_path: str | Path | None = None,
) -> dict[str, str]:
    """Validate, aggregate, and plot an introgressed-sweep replicate campaign.

    Parameters are deliberately explicit so a campaign may live outside the
    core single-realization study directories.  Only checksum- and
    contract-valid simulation/decode caches contribute scientific values.
    """
    runtime.validate()
    bases, expanded = _normalise_specifications(base_specs, expanded_specs)
    campaign_dir = Path(study_dir).resolve()
    results_dir = campaign_dir / "results"
    figure_dir = results_dir / "figures"
    results_dir.mkdir(parents=True, exist_ok=True)

    resolved_decoder: Path | None = None
    decoder_sha256 = None
    if decoder_path is not None:
        resolved_decoder = Path(decoder_path).resolve()
        if not resolved_decoder.is_file():
            raise FileNotFoundError(f"Gamma-SMC decoder not found: {resolved_decoder}")
        decoder_sha256 = _sha256(resolved_decoder)

    snapshot = _build_input_snapshot(campaign_dir, expanded)
    snapshot_root = Path(snapshot["root"])
    observed_decoder_sha256s = _observed_decoder_sha256s(
        snapshot_root, expanded, decoder_sha256
    )
    status, acceptance, contrasts, profiles = _collect_replicates(
        runtime,
        expanded,
        campaign_dir,
        resolved_decoder,
        source_study_dir=snapshot_root,
    )
    _verify_snapshot_sources_unchanged(campaign_dir, expanded, snapshot)
    coverage = _cell_coverage(bases, expanded, status)
    effect_summary = _effect_summary(runtime, bases, coverage, contrasts)
    paired = _paired_contrasts(contrasts)
    paired_accuracy = _paired_accuracy(bases, paired)
    aggregate_profiles = _aggregate_profiles(profiles, runtime)
    threshold_accuracy = _threshold_accuracy(profiles)

    artifact_paths: dict[str, Path] = {
        "input_snapshot_manifest_json": Path(snapshot["manifest_path"]),
        "coverage_input_inventory_tsv": _atomic_frame(
            results_dir / "provenance" / "coverage_input_inventory.tsv",
            pd.DataFrame(snapshot["manifest"]["files"]),
        ),
        "replicate_simulation_status_tsv": _atomic_frame(
            results_dir / "replicate_simulation_status.tsv", status
        ),
        "cell_coverage_tsv": _atomic_frame(results_dir / "cell_coverage.tsv", coverage),
        "replicate_acceptance_tsv": _atomic_frame(
            results_dir / "replicate_acceptance.tsv", acceptance
        ),
        "replicate_center_contrasts_tsv": _atomic_frame(
            results_dir / "replicate_center_contrasts.tsv", contrasts
        ),
        "replicate_center_profiles_tsv_gz": _atomic_frame(
            results_dir / "replicate_center_profiles.tsv.gz", profiles
        ),
        "cell_effect_summary_tsv": _atomic_frame(
            results_dir / "cell_effect_summary.tsv", effect_summary
        ),
        "truth_gamma_paired_replicates_tsv": _atomic_frame(
            results_dir / "truth_gamma_paired_replicates.tsv", paired
        ),
        "truth_gamma_accuracy_tsv": _atomic_frame(
            results_dir / "truth_gamma_accuracy.tsv", paired_accuracy
        ),
        "aggregate_cdf_profiles_tsv": _atomic_frame(
            results_dir / "aggregate_cdf_profiles.tsv", aggregate_profiles
        ),
        "threshold_accuracy_tsv": _atomic_frame(
            results_dir / "threshold_accuracy.tsv", threshold_accuracy
        ),
    }

    plots = {
        "coverage": _plot_coverage(coverage, figure_dir / "replicate_coverage"),
        "achieved_af_genotypes": _plot_achieved(
            acceptance, figure_dir / "replicate_achieved_af_genotypes"
        ),
        "truth_vs_gamma": _plot_truth_vs_gamma(
            paired, figure_dir / "replicate_truth_vs_gamma"
        ),
        "effect_sign": _plot_effect_sign(
            effect_summary, figure_dir / "replicate_effect_sign"
        ),
        "tmrca_cdf": _plot_cdf_profiles(
            aggregate_profiles, bases, figure_dir / "replicate_tmrca_cdf"
        ),
        "threshold_mae": _plot_threshold_mae(
            threshold_accuracy, figure_dir / "replicate_threshold_mae"
        ),
    }
    for plot_name, formats in plots.items():
        for suffix, path in formats.items():
            artifact_paths[f"figure_{plot_name}_{suffix}"] = path

    status_counts = {
        str(key): int(value)
        for key, value in status["simulation_status"].value_counts().items()
    }
    headline_accuracy = paired_accuracy[
        paired_accuracy["summary_scope"] == "overall_macro_cell_equal"
    ].to_dict(orient="records")
    secondary_micro_accuracy = paired_accuracy[
        paired_accuracy["summary_scope"] == "overall_micro_trajectory_weighted"
    ].to_dict(orient="records")
    summary_payload = _json_clean(
        {
            "schema": REPLICATE_RESULTS_SCHEMA,
            "requested_biological_cells": len(bases),
            "requested_trajectories": len(expanded),
            "accepted_trajectories": len(acceptance),
            "decoded_trajectories": int(
                acceptance.get("decoded", pd.Series(dtype=bool)).sum()
            ),
            "attempted_trajectories": int(coverage["attempted_replicates"].sum()),
            "completed_external_draws": int(coverage["completed_external_draws"].sum()),
            "launched_external_draws": int(coverage["launched_external_draws"].sum()),
            "interrupted_internal_conditioned_trajectories": int(
                coverage["interrupted_internal_conditioned_trajectories"].sum()
            ),
            "cumulative_simulation_elapsed_seconds": float(
                coverage["cumulative_elapsed_seconds"].sum()
            ),
            "simulation_status_counts": status_counts,
            "headline_macro_cell_equal_truth_gamma_accuracy": headline_accuracy,
            "secondary_micro_trajectory_weighted_truth_gamma_accuracy": (
                secondary_micro_accuracy
            ),
            "headline_minimum_paired_replicates_per_cell": (MIN_CELL_PAIRED_REPLICATES),
            "interpretation": [
                "The trajectory, not a within-individual pair, is the replicate unit.",
                "Conditioning failures are nonrandom; effect summaries describe accepted trajectories only.",
                "Headline overall accuracy is an equal-cell macro mean over biological cells meeting the explicit paired-replicate minimum.",
                "The secondary micro summary pools accepted trajectories and is therefore trajectory-weighted across nonrandom cell coverage.",
                "Overall TMRCA-threshold accuracy is a secondary accepted-trajectory-weighted micro summary across biological cells.",
                "Truth and Gamma expected-sign rates are reported separately.",
                "Raw truth/Gamma sign agreement is descriptive and may include jointly wrong signs.",
                "Gamma detection uses only the denominator of replicates with the expected tree-truth sign; jointly wrong signs are never successes.",
                "Normalized CDF area is a within-genotype TMRCA-CDF area, not ROC AUC.",
                "No neutral-null p-values or power claims are produced.",
            ],
        }
    )
    summary_path = _atomic_json(
        results_dir / "replicate_results_summary.json", summary_payload
    )
    artifact_paths["replicate_results_summary_json"] = summary_path
    try:
        campaign_design = json.loads(
            (snapshot_root / "study_design.json").read_text(encoding="utf-8")
        )
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError("RUN_RESULTS campaign design is unreadable") from error
    if not isinstance(campaign_design, Mapping):
        raise ValueError("RUN_RESULTS campaign design must be a JSON object")
    report_path = _atomic_text(
        campaign_dir / "RUN_RESULTS.md",
        _render_run_results_markdown(
            summary_payload,
            coverage,
            acceptance,
            paired_accuracy,
            effect_summary,
            threshold_accuracy,
            _report_campaign_facts(runtime, bases, campaign_design),
        ),
    )
    artifact_paths["run_results_markdown"] = report_path

    input_provenance = _campaign_input_provenance(
        campaign_dir, bases, expanded, runtime, snapshot
    )
    completion_provenance = _completion_provenance(
        campaign_dir,
        expanded,
        status,
        decoder_sha256,
        source_study_dir=snapshot_root,
        observed_decoder_sha256s=observed_decoder_sha256s,
    )
    manifest_payload = {
        "schema": REPLICATE_RESULTS_SCHEMA,
        "phase": "aggregate",
        "runtime": _portable_runtime_record(runtime, campaign_dir),
        "base_specification_count": len(bases),
        "expanded_specification_count": len(expanded),
        "decoder_validation": (
            "explicit_decoder_sha256"
            if decoder_path is not None
            else "recorded_contract"
        ),
        "decoder_sha256": decoder_sha256,
        "decoder_sha256s": observed_decoder_sha256s,
        "decoder_stratification": "single_sha_required_fail_closed",
        "replicate_unit": "independently_seeded_accepted_slim_trajectory",
        "conditioning_missingness": "nonrandom",
        "neutral_null": False,
        "environment": _environment_provenance(Path(runtime.repo_root)),
        "inputs": input_provenance,
        "validated_contracts_and_completions": completion_provenance,
        "artifacts": {
            name: {
                "path": path.relative_to(campaign_dir).as_posix(),
                "sha256": _sha256(path),
            }
            for name, path in sorted(artifact_paths.items())
        },
    }
    manifest_path = _atomic_json(
        results_dir / "replicate_results_manifest.json", manifest_payload
    )
    artifact_paths["replicate_results_manifest_json"] = manifest_path
    _verify_snapshot_sources_unchanged(campaign_dir, expanded, snapshot)
    return {name: str(path) for name, path in sorted(artifact_paths.items())}


def aggregate_introgression_replicates(
    runtime: StudyRuntime,
    base_specs: Sequence[Mapping[str, Any]],
    expanded_specs: Sequence[Mapping[str, Any]],
    study_dir: str | Path,
    decoder_path: str | Path | None = None,
) -> dict[str, str]:
    """Aggregate one quiescent campaign from an immutable compact snapshot."""
    campaign_dir = Path(study_dir).resolve()
    with _campaign_aggregation_lock(campaign_dir):
        _assert_campaign_quiescent(campaign_dir)
        artifacts = _aggregate_introgression_replicates_locked(
            runtime,
            base_specs,
            expanded_specs,
            campaign_dir,
            decoder_path,
        )
        _assert_campaign_quiescent(campaign_dir)
        return artifacts


__all__ = [
    "BOOTSTRAP_DRAWS",
    "INPUT_SNAPSHOT_SCHEMA",
    "MIN_CELL_PAIRED_REPLICATES",
    "REPLICATE_RESULTS_SCHEMA",
    "aggregate_introgression_replicates",
]
