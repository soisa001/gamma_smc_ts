from __future__ import annotations

import json
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from gamma_smc_aou import focused_han_direct_v3 as v3
from gamma_smc_aou import focused_han_direct_v3_recovery_audit as audit_adapter


def _real_shape_frame() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "unit_id": "no-interruption",
                "interrupted_nonterminal_dispositions": 0,
                "reconciliation_path": None,
                "reconciliation_sha256": None,
                "committed_elapsed_seconds": 10.0,
            },
            {
                "unit_id": "has-interruption",
                "interrupted_nonterminal_dispositions": 1,
                "reconciliation_path": (
                    "focused_selection_EAS_sim/work/has-interruption/"
                    "supervision_reconciliation.tsv"
                ),
                "reconciliation_sha256": "a" * 64,
                "committed_elapsed_seconds": 11.0,
            },
        ]
    )


def test_mixed_real_shape_normalizes_only_structural_absence():
    frame = _real_shape_frame()
    assert np.isnan(frame.loc[0, "reconciliation_path"])
    assert np.isnan(frame.loc[0, "reconciliation_sha256"])

    normalized = audit_adapter._normalize_reconciliation_absence(frame)
    records = normalized.to_dict(orient="records")
    assert records[0]["reconciliation_path"] is None
    assert records[0]["reconciliation_sha256"] is None
    assert records[1] == frame.to_dict(orient="records")[1]
    assert records[0]["committed_elapsed_seconds"] == 10.0
    assert isinstance(audit_adapter._canonical_sha256(records), str)


@pytest.mark.parametrize(
    ("column", "value", "message"),
    [
        ("committed_elapsed_seconds", np.nan, "unexpected nonfinite value"),
        ("committed_elapsed_seconds", np.inf, "unexpected nonfinite value"),
        ("reconciliation_path", np.inf, "unexpected nonfinite value"),
    ],
)
def test_unexpected_nonfinite_values_fail_closed(column, value, message):
    frame = _real_shape_frame()
    frame[column] = frame[column].astype("object")
    frame.at[0, column] = value
    with pytest.raises(ValueError, match=message):
        audit_adapter._normalize_reconciliation_absence(frame)


def test_reconciliation_reference_cannot_be_missing_when_count_is_positive():
    frame = _real_shape_frame()
    frame.at[1, "reconciliation_sha256"] = np.nan
    with pytest.raises(ValueError, match="unexpected nonfinite value"):
        audit_adapter._normalize_reconciliation_absence(frame)


def test_structural_absence_cannot_mask_an_unexpected_string():
    frame = _real_shape_frame()
    frame.at[0, "reconciliation_path"] = "unexpected.tsv"
    with pytest.raises(ValueError, match="structural absence differs"):
        audit_adapter._normalize_reconciliation_absence(frame)


def test_patch_is_scoped_and_restores_frozen_aggregator(monkeypatch):
    frame = _real_shape_frame()

    def frozen_aggregate(_bundle, _stage=v3.STAGE_ALL):
        return frame.copy()

    monkeypatch.setattr(v3, "aggregate_v3_cached_completions", frozen_aggregate)
    with audit_adapter.patched_recovery_audit_aggregation():
        normalized = v3.aggregate_v3_cached_completions(SimpleNamespace(), v3.STAGE_ALL)
        assert normalized.loc[0, "reconciliation_path"] is None
    assert v3.aggregate_v3_cached_completions is frozen_aggregate


def test_finalizer_refuses_incomplete_state_before_any_audit_write(
    tmp_path, monkeypatch
):
    authorization_path = tmp_path / "han_direct_v3_dtype_recovery_authorization.json"
    authorization_path.write_text("{}", encoding="utf-8")
    bundle = SimpleNamespace(output_dir=tmp_path)
    authorization = {
        "dry_plan_payload_sha256": "p" * 64,
        "recoverable_unit_ids": ["unit"],
    }
    monkeypatch.setattr(
        audit_adapter.recovery,
        "_load_authorization",
        lambda _bundle, _path: authorization,
    )
    monkeypatch.setattr(
        audit_adapter.recovery,
        "_authorized_current_state",
        lambda _bundle, _authorization: ([], [{"unit_id": "unit"}]),
    )
    with pytest.raises(ValueError, match="refuses to launch simulations"):
        audit_adapter.finalize_recovery_audit(bundle, expected_plan_sha256="p" * 64)


def test_adapter_has_no_simulation_capable_recovery_call():
    assert "run_recovery" not in audit_adapter.finalize_recovery_audit.__code__.co_names
    assert "_simulate_one_recovery" not in audit_adapter.__dict__


def test_completed_audit_publisher_records_zero_scheduled_units(tmp_path, monkeypatch):
    output = tmp_path / "output"
    output.mkdir()
    authorization_path = output / "han_direct_v3_dtype_recovery_authorization.json"
    authorization_path.write_text("{}\n", encoding="utf-8")
    history_path = output / "recovery_history" / "history.json"
    history = {
        "schema": audit_adapter.recovery.RECOVERY_HISTORY_SCHEMA,
        "status": "immutable_stopped_state_snapshot",
        "file_count": 2,
        "total_bytes": 10,
    }
    history["payload_sha256"] = audit_adapter._canonical_sha256(history)
    history_path.parent.mkdir()
    history_path.write_text(json.dumps(history) + "\n", encoding="utf-8")
    authorization = {
        "payload_sha256": "a" * 64,
        "dry_plan_payload_sha256": "p" * 64,
        "history_manifest_path": history_path.relative_to(tmp_path).as_posix(),
        "history_manifest_sha256": v3.sha256_file(history_path),
        "history_payload_sha256": history["payload_sha256"],
        "manifest_sha256": "m" * 64,
        "seed_plan_sha256": "s" * 64,
        "original_implementation_sources": {"base": "b" * 64},
        "recovery_sources": {"recovery": "r" * 64},
        "recoverable_unit_ids": ["unit"],
    }
    completed = [
        {
            "unit_id": "unit",
            "completion_sha256": "c" * 64,
            "attempts_completed": 2,
        }
    ]
    base_audit = {
        "status": "passed",
        "accepted_units": 60,
        "payload_sha256": "f" * 64,
    }
    adapter_binding = {
        "sources": {"adapter": "z" * 64},
        "audit_only_no_simulation_route": True,
    }
    monkeypatch.setattr(
        audit_adapter,
        "_audit_adapter_binding",
        lambda _bundle: adapter_binding,
    )
    payload = audit_adapter._publish_completed_recovery_audit(
        SimpleNamespace(repo_root=tmp_path, output_dir=output),
        authorization,
        completed,
        base_audit,
    )
    assert payload["invocation_scheduled_units"] == 0
    assert payload["pre_invocation_incomplete_unit_ids"] == []
    assert payload["post_recovery_recovered_unit_ids"] == ["unit"]
    assert payload["post_run_base_audit"] == base_audit
    assert payload["audit_adapter"] == adapter_binding
    assert (output / audit_adapter.recovery.CANONICAL_AUDIT_NAME).is_file()
