from __future__ import annotations

import json
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from gamma_smc_aou import focused_han_direct_v3 as v3
from gamma_smc_aou import focused_han_direct_v3_recovery as recovery
from gamma_smc_aou import focused_selection_simulation as simulation


def _fake_frozen_atomic_capture(captured: list[pd.DataFrame]):
    def atomic(_path: Path, frame: pd.DataFrame) -> None:
        enriched = frame.copy()
        if "rejection_reason" not in enriched:
            enriched["rejection_reason"] = ""
        timed_out = enriched["status"].astype(str).eq("timed_out")
        blank = enriched["rejection_reason"].fillna("").astype(str).eq("")
        enriched.loc[timed_out & blank, "rejection_reason"] = "draw_timeout"
        captured.append(enriched)

    return atomic


def _exercise_dtype_context(
    monkeypatch: pytest.MonkeyPatch, frame: pd.DataFrame
) -> pd.DataFrame:
    captured: list[pd.DataFrame] = []
    prior_atomic = simulation._atomic_frame
    prior_pandas = v3.pd

    @contextmanager
    def fake_frozen_executor(_bundle, _unit):
        original = simulation._atomic_frame
        simulation._atomic_frame = _fake_frozen_atomic_capture(captured)
        try:
            yield SimpleNamespace()
        finally:
            simulation._atomic_frame = original

    monkeypatch.setattr(recovery, "_FROZEN_PATCHED_EXECUTOR", fake_frozen_executor)
    with recovery.patched_recovery_executor(SimpleNamespace(), {}):
        simulation._atomic_frame(Path("attempts.tsv"), frame)
    assert simulation._atomic_frame is prior_atomic
    assert v3.pd is prior_pandas
    assert len(captured) == 1
    return captured[0]


def test_timeout_then_scientific_reject_uses_object_reason_dtype(monkeypatch):
    first = pd.DataFrame(
        {
            "attempt_zero_based": [0],
            "status": ["timed_out"],
            "accepted": [False],
            "rejection_reason": pd.Series([np.nan], dtype="float64"),
        }
    )
    first_written = _exercise_dtype_context(monkeypatch, first)
    assert first_written["rejection_reason"].dtype == object
    assert first_written["rejection_reason"].tolist() == ["draw_timeout"]

    second = pd.DataFrame(
        {
            "attempt_zero_based": [0, 1],
            "status": ["timed_out", "complete"],
            "accepted": [False, False],
            "rejection_reason": [
                "draw_timeout",
                "candidate_pool_af_outside_prespecified_band",
            ],
        }
    )
    second_written = _exercise_dtype_context(monkeypatch, second)
    assert second_written["rejection_reason"].tolist() == [
        "draw_timeout",
        "candidate_pool_af_outside_prespecified_band",
    ]


def test_timeout_then_accept_retains_blank_accept_reason(monkeypatch):
    frame = pd.DataFrame(
        {
            "attempt_zero_based": [0, 1],
            "status": ["timed_out", "complete"],
            "accepted": [False, True],
            "rejection_reason": pd.Series([np.nan, np.nan], dtype="float64"),
        }
    )
    written = _exercise_dtype_context(monkeypatch, frame)
    assert written["rejection_reason"].dtype == object
    assert written.iloc[0]["rejection_reason"] == "draw_timeout"
    assert pd.isna(written.iloc[1]["rejection_reason"])
    assert bool(written.iloc[1]["accepted"])


def test_missing_reason_column_is_added_as_object_without_other_changes():
    frame = pd.DataFrame(
        {
            "attempt_zero_based": [0],
            "seed": [123],
            "status": ["complete"],
            "accepted": [True],
        }
    )
    normalized = recovery._coerce_rejection_reason_object(frame)
    pd.testing.assert_frame_equal(
        normalized.drop(columns="rejection_reason"), frame, check_dtype=True
    )
    assert normalized["rejection_reason"].dtype == object
    assert normalized["rejection_reason"].tolist() == [""]


def test_restart_read_proxy_prevents_blank_reason_float_inference(tmp_path: Path):
    path = tmp_path / "attempts.tsv"
    path.write_text(
        "attempt_zero_based\tstatus\trejection_reason\n0\tcomplete\t\n",
        encoding="utf-8",
    )
    ordinary = pd.read_csv(path, sep="\t")
    assert ordinary["rejection_reason"].dtype == np.dtype("float64")
    recovered = recovery._RecoveryPandasProxy.read_csv(path, sep="\t")
    assert recovered["rejection_reason"].dtype == object
    recovered.loc[0, "rejection_reason"] = (
        "worker_interrupted_after_base_attempt_commit"
    )
    assert recovered.loc[0, "rejection_reason"].startswith("worker_interrupted")


def test_postcommit_active_restart_uses_frozen_reconciliation_with_dtype_proxy(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    unit_dir = tmp_path / "unit"
    unit_dir.mkdir()
    schedule = pd.DataFrame(
        [{"attempt_zero_based": 0, "simulation_seed": 12345, "panel_seed": 54321}]
    )
    attempts = pd.DataFrame(
        [
            {
                "attempt_zero_based": 0,
                "seed": 12345,
                "status": "complete",
                "timeout_seconds": 300.0,
                "elapsed_seconds": 1.0,
                "error": "",
                "accepted": False,
                "rejection_reason": np.nan,
                "planned_simulation_seed": 12345,
                "panel_seed": 54321,
                "panel_seed_used": False,
            }
        ]
    )
    attempts.to_csv(unit_dir / "attempts.tsv", sep="\t", index=False)
    record = {
        "schema": v3.SUPERVISOR_LEDGER_SCHEMA,
        "attempt_zero_based": 0,
        "simulation_seed": 12345,
        "supervisor_started_unix_ns": 1,
        "supervisor_pid": 999_001,
        "supervisor_process_identity": {
            "pid": 999_001,
            "boot_id": "dead",
            "start_ticks": 1,
        },
        "timeout_seconds": 300.0,
        "term_grace_seconds": 0.01,
        "temp_directory": ".han_v3_attempt_000_seed_12345_test",
        "temp_tree_created": True,
        "temp_tree_loaded": True,
        "temp_tree_removed": True,
        "child_pid": 999_002,
        "child_pgid": 999_002,
        "child_sid": 999_002,
        "child_process_identity": {
            "pid": 999_002,
            "boot_id": "dead",
            "start_ticks": 2,
        },
        "child_identity_observed": True,
        "engine_wrapper_pid": 999_003,
        "engine_wrapper_process_identity": {
            "pid": 999_003,
            "boot_id": "dead",
            "start_ticks": 3,
        },
        "direct_child_reaped": True,
        "adopted_children_reaped": 0,
        "process_group_gone": True,
        "supervisor_status": "complete",
        "supervisor_elapsed_seconds": 1.0,
        "lifecycle_phase": "engine_returned_pending_base_attempt_commit",
    }
    record["execution_id"] = v3._SupervisedSlimEngine._execution_id(record)
    pd.DataFrame([record]).to_csv(
        unit_dir / "supervision.tsv", sep="\t", index=False
    )
    v3._write_active_supervision(
        unit_dir / ".han_v3_active_000.json",
        {
            "schema": v3.SUPERVISOR_ACTIVE_SCHEMA,
            "status": "active",
            "attempt_record": record,
        },
    )
    monkeypatch.setattr(v3, "_same_linux_process", lambda _identity: False)
    monkeypatch.setattr(v3, "_process_group_has_live_members", lambda _pgid: False)
    monkeypatch.setattr(v3, "_process_absent_or_zombie", lambda _pid: True)
    prior = v3.pd
    v3.pd = recovery._RecoveryPandasProxy()
    try:
        v3._SupervisedSlimEngine(object(), unit_dir=unit_dir, schedule=schedule)
    finally:
        v3.pd = prior
    reconciled = pd.read_csv(unit_dir / "attempts.tsv", sep="\t")
    assert reconciled.loc[0, "status"] == "timed_out"
    assert reconciled.loc[0, "rejection_reason"] == (
        "worker_interrupted_after_base_attempt_commit"
    )
    assert int(reconciled.loc[0, "seed"]) == 12345
    assert not list(unit_dir.glob(".han_v3_active_*.json"))
    journal = pd.read_csv(unit_dir / "supervision_reconciliation.tsv", sep="\t")
    assert journal.loc[0, "reason"] == "attempt_interrupted_after_base_commit"


def test_active_preflight_rejects_live_or_ambiguous_marker(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    unit_dir = tmp_path / "unit"
    unit_dir.mkdir()
    unit = {"unit_id": "unit"}
    schedule = pd.DataFrame(
        [{"attempt_zero_based": 0, "simulation_seed": 123, "panel_seed": 456}]
    )
    monkeypatch.setattr(v3, "_unit_seed_schedule", lambda _bundle, _unit_id: schedule)
    record = {
        "schema": v3.SUPERVISOR_LEDGER_SCHEMA,
        "attempt_zero_based": 0,
        "simulation_seed": 123,
        "supervisor_started_unix_ns": 1,
        "supervisor_pid": 12,
        "supervisor_process_identity": {
            "pid": 12,
            "boot_id": "live",
            "start_ticks": 2,
        },
        "timeout_seconds": 300.0,
        "temp_directory": ".han_v3_attempt_000_seed_123_test",
    }
    record["execution_id"] = v3._SupervisedSlimEngine._execution_id(record)
    for attempt in (0,):
        v3._write_active_supervision(
            unit_dir / f".han_v3_active_{attempt:03d}.json",
            {
                "schema": v3.SUPERVISOR_ACTIVE_SCHEMA,
                "status": "active",
                "attempt_record": record,
            },
        )
    monkeypatch.setattr(v3, "_same_linux_process", lambda _identity: True)
    with pytest.raises(ValueError, match="supervisor process is still alive"):
        recovery._active_marker_record(
            SimpleNamespace(), unit, unit_dir, pd.DataFrame(), pd.DataFrame()
        )
    (unit_dir / ".han_v3_active_001.json").write_text(
        (unit_dir / ".han_v3_active_000.json").read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="multiple active markers"):
        recovery._active_marker_record(
            SimpleNamespace(), unit, unit_dir, pd.DataFrame(), pd.DataFrame()
        )


def test_unstarted_unit_is_authorized_from_attempt_zero_without_writes(tmp_path: Path):
    campaign = tmp_path / "campaign"
    campaign.mkdir()
    unit = {
        "unit_id": "unstarted",
        "replicate_index": 2,
        "selection_coefficient": 0.005,
        "target_allele_frequency": 0.3,
    }
    row = recovery._inspect_incomplete_unit(
        SimpleNamespace(campaign_dir=campaign), unit
    )
    assert row["state"] == "unstarted"
    assert row["committed_attempts"] == 0
    assert row["active_simulation_seed"] is None
    assert row["artifact_file_count"] == 0
    assert not (campaign / "work" / "unstarted").exists()


def test_quiescent_crash_boundary_between_attempts_is_restartable(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    campaign = tmp_path / "campaign"
    unit_dir = campaign / "work" / "unit"
    unit_dir.mkdir(parents=True)
    (unit_dir / "simulation_contract.json").write_text("{}", encoding="utf-8")
    unit = {
        "unit_id": "unit",
        "replicate_index": 2,
        "selection_coefficient": 0.005,
        "target_allele_frequency": 0.3,
    }
    attempts = pd.DataFrame(
        [{"attempt_zero_based": 0, "accepted": False, "seed": 123}]
    )
    supervision = pd.DataFrame(
        [{"attempt_zero_based": 0, "simulation_seed": 123}]
    )
    monkeypatch.setattr(
        recovery,
        "_strict_partial_contract",
        lambda _bundle, _unit, _dir: ({"bound": True}, "c" * 64),
    )
    monkeypatch.setattr(
        recovery,
        "_exact_attempts",
        lambda _bundle, _unit, _dir: (attempts, supervision),
    )
    monkeypatch.setattr(
        recovery, "_active_marker_record", lambda *_args: None
    )
    monkeypatch.setattr(recovery, "_known_dtype_failure", lambda *_args, **_kw: None)
    row = recovery._inspect_incomplete_unit(
        SimpleNamespace(campaign_dir=campaign), unit
    )
    assert row["state"] == "quiescent_between_attempts"
    assert row["committed_attempts"] == 1

    attempts.loc[0, "accepted"] = True
    with pytest.raises(ValueError, match="no recoverable disposition"):
        recovery._inspect_incomplete_unit(
            SimpleNamespace(campaign_dir=campaign), unit
        )


def test_history_snapshot_preserves_failure_and_active_bytes(tmp_path: Path):
    root = tmp_path / "repo"
    campaign = root / "campaign"
    output = root / "output"
    unit_dir = campaign / "work" / "unit"
    unit_dir.mkdir(parents=True)
    output.mkdir(parents=True)
    (unit_dir / "simulation_failed.json").write_bytes(b"failure-bytes\n")
    (unit_dir / ".han_v3_active_000.json").write_bytes(b"active-bytes\n")
    inventory, total = recovery._unit_inventory(unit_dir)
    plan = {
        "schema": recovery.RECOVERY_PLAN_SCHEMA,
        "status": "ready_read_only",
        "manifest_sha256": "m" * 64,
        "seed_plan_sha256": "s" * 64,
        "original_implementation_sources": {"base": "b" * 64},
        "recovery_sources": {"adapter": "a" * 64},
        "recoverable_units": 1,
        "recover": [
            {
                "unit_id": "unit",
                "artifact_file_count": len(inventory),
                "artifact_total_bytes": total,
                "artifact_inventory_sha256": recovery._canonical_sha256(inventory),
            }
        ],
    }
    plan["payload_sha256"] = recovery._canonical_sha256(plan)
    bundle = recovery.RecoveryBundle(
        direct=SimpleNamespace(
            repo_root=root,
            campaign_dir=campaign,
            output_dir=output,
        ),
        recovery_sources={},
    )
    manifest_path, history = recovery.snapshot_recovery_history(bundle, plan)
    assert history["file_count"] == 2
    snapshot_root = manifest_path.parent / "units" / "unit"
    assert (snapshot_root / "simulation_failed.json").read_bytes() == b"failure-bytes\n"
    assert (snapshot_root / ".han_v3_active_000.json").read_bytes() == b"active-bytes\n"
    (unit_dir / "simulation_failed.json").write_bytes(b"changed\n")
    assert (snapshot_root / "simulation_failed.json").read_bytes() == b"failure-bytes\n"
    same_path, same_history = recovery.snapshot_recovery_history(bundle, plan)
    assert same_path == manifest_path
    assert same_history == history


def test_failure_preflight_accepts_only_exact_dtype_bug(tmp_path: Path):
    path = tmp_path / "simulation_failed.json"
    payload = {
        "schema": simulation.SCHEMA_VERSION,
        "status": "failed",
        "unit_id": "unit",
        "error_type": "TypeError",
        "error": "Invalid value 'draw_timeout' for dtype 'float64'",
        "contract_sha256": "c" * 64,
        "search_budget": {
            "selected_max_draws": v3.ATTEMPTS_PER_UNIT,
            "selected_draw_timeout_seconds": v3.DRAW_TIMEOUT_SECONDS,
            "selected_cumulative_timeout_seconds": v3.CUMULATIVE_TIMEOUT_SECONDS,
        },
    }
    path.write_text(json.dumps(payload), encoding="utf-8")
    assert recovery._known_dtype_failure(
        path, unit_id="unit", contract_sha256="c" * 64
    ) == payload
    payload["error"] = "some other TypeError"
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="non-dtype failure"):
        recovery._known_dtype_failure(
            path, unit_id="unit", contract_sha256="c" * 64
        )


def test_authorized_rerun_schedules_only_still_incomplete_units(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    campaign = tmp_path / "campaign"
    units = pd.DataFrame(
        [
            {"unit_id": "unit-a", "replicate_index": 2},
            {"unit_id": "unit-b", "replicate_index": 3},
        ]
    )
    bundle = recovery.RecoveryBundle(
        direct=SimpleNamespace(campaign_dir=campaign), recovery_sources={}
    )
    authorization = {"recoverable_unit_ids": ["unit-a", "unit-b"]}
    monkeypatch.setattr(recovery, "EXPECTED_STOPPED_RECOVERABLE_UNITS", 2)
    monkeypatch.setattr(v3, "_stage_units", lambda _bundle, _stage: units)
    monkeypatch.setattr(
        v3,
        "validate_v3_cached_unit",
        lambda _bundle, unit: {
            "unit_id": unit["unit_id"],
            "completion_sha256": unit["unit_id"] + "-sha",
            "attempts_completed": 2,
        },
    )
    monkeypatch.setattr(
        recovery,
        "_inspect_incomplete_unit",
        lambda _bundle, unit: {"unit_id": unit["unit_id"], "state": "unstarted"},
    )
    first_completion = campaign / "work" / "unit-a" / "simulation_complete.json"
    first_completion.parent.mkdir(parents=True)
    first_completion.write_text("{}", encoding="utf-8")
    completed, incomplete = recovery._authorized_current_state(
        bundle, authorization
    )
    assert [row["unit_id"] for row in completed] == ["unit-a"]
    assert [row["unit_id"] for row in incomplete] == ["unit-b"]

    second_completion = campaign / "work" / "unit-b" / "simulation_complete.json"
    second_completion.parent.mkdir(parents=True)
    second_completion.write_text("{}", encoding="utf-8")
    completed, incomplete = recovery._authorized_current_state(
        bundle, authorization
    )
    assert [row["unit_id"] for row in completed] == ["unit-a", "unit-b"]
    assert incomplete == []
