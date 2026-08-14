from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
from contextlib import contextmanager
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from gamma_smc_aou import focused_han_direct_v3 as han_v3
from gamma_smc_aou import focused_selected_postproduction as postproduction


def _signed(payload: dict) -> dict:
    result = dict(payload)
    result["payload_sha256"] = postproduction._canonical_sha256(result)
    return result


def test_canonical_text_hash_accepts_lf_and_crlf_without_unicode_normalization(
    tmp_path,
):
    lf = tmp_path / "lf.py"
    crlf = tmp_path / "crlf.py"
    decomposed = "value = 'e\N{COMBINING ACUTE ACCENT}'\n".encode()
    lf.write_bytes(decomposed)
    crlf.write_bytes(decomposed.replace(b"\n", b"\r\n"))

    expected = hashlib.sha256(decomposed).hexdigest()
    assert postproduction.canonical_text_sha256(lf) == expected
    assert postproduction.canonical_text_sha256(crlf) == expected
    assert (
        expected
        != hashlib.sha256(
            "value = '\N{LATIN SMALL LETTER E WITH ACUTE}'\n".encode()
        ).hexdigest()
    )


@pytest.mark.parametrize(
    ("content", "message"),
    [
        (b"\xef\xbb\xbfvalue = 1\n", "UTF-8 BOM"),
        (b"value = 1\rvalue = 2\n", "bare carriage return"),
        (b"value = '\xff'\n", "not strict UTF-8"),
    ],
)
def test_canonical_text_hash_rejects_unsafe_encodings(tmp_path, content, message):
    source = tmp_path / "source.py"
    source.write_bytes(content)
    with pytest.raises(ValueError, match=message):
        postproduction.canonical_text_sha256(source)


def test_canonical_source_binding_accepts_core_autocrlf_checkout(tmp_path):
    if shutil.which("git") is None:
        pytest.skip("git is required for the core.autocrlf checkout test")
    relative = "adapter.py"
    source = tmp_path / relative
    source.write_bytes(b"line_one = 1\nline_two = 2\n")
    recorded = {relative: postproduction.canonical_text_sha256(source)}
    subprocess.run(["git", "init", "-q"], cwd=tmp_path, check=True)
    subprocess.run(
        ["git", "-c", "core.autocrlf=false", "add", relative],
        cwd=tmp_path,
        check=True,
    )
    source.unlink()
    subprocess.run(
        ["git", "-c", "core.autocrlf=true", "checkout-index", "-f", "--", relative],
        cwd=tmp_path,
        check=True,
    )
    assert source.read_bytes() == b"line_one = 1\r\nline_two = 2\r\n"

    observed = postproduction._validate_canonical_source_records(
        recorded,
        repo_root=tmp_path,
        relative_paths=(relative,),
        label="fresh checkout",
    )
    assert observed == recorded


def test_raw_nan_fails_without_context_and_succeeds_inside_context(monkeypatch):
    frame = pd.DataFrame(
        [
            {
                "unit_id": "no-interruption",
                "interrupted_nonterminal_dispositions": 0,
                "reconciliation_path": None,
                "reconciliation_sha256": None,
            },
            {
                "unit_id": "interrupted",
                "interrupted_nonterminal_dispositions": 1,
                "reconciliation_path": "work/interrupted/supervision_reconciliation.tsv",
                "reconciliation_sha256": "a" * 64,
            },
        ]
    )
    assert np.isnan(frame.loc[0, "reconciliation_path"])

    def aggregate(_bundle, _stage=han_v3.STAGE_ALL):
        return frame.copy()

    monkeypatch.setattr(han_v3, "aggregate_v3_cached_completions", aggregate)
    with pytest.raises(ValueError, match="Out of range float values"):
        postproduction._canonical_sha256(frame.to_dict(orient="records"))
    with postproduction.recovery_audit.patched_recovery_audit_aggregation():
        normalized = han_v3.aggregate_v3_cached_completions(
            SimpleNamespace(), han_v3.STAGE_SIMULATE
        )
        digest = postproduction._canonical_sha256(normalized.to_dict(orient="records"))
    assert len(digest) == 64
    assert normalized.loc[0, "reconciliation_path"] is None


def test_hashed_normalization_corruption_fails_closed():
    payload = _signed(
        {
            "schema": postproduction.recovery_audit.AUDIT_NORMALIZATION_SCHEMA,
            "status": "passed",
            "validated_rows": 60,
        }
    )
    payload["validated_rows"] = 59
    with pytest.raises(ValueError, match="corrupt or incompatible"):
        postproduction._validate_hashed_payload(
            payload,
            schema=postproduction.recovery_audit.AUDIT_NORMALIZATION_SCHEMA,
            status="passed",
            label="normalization",
        )


@pytest.mark.parametrize(
    ("absent", "present"),
    [
        (34, 26),
        (36, 24),
        (True, 25),
        (35, False),
    ],
)
def test_normalization_partition_requires_exact_strict_integer_counts(absent, present):
    with pytest.raises(ValueError, match="row partition differs"):
        postproduction._validate_normalization_partition(
            {
                "structurally_absent_rows": absent,
                "reconciliation_present_rows": present,
            }
        )


def test_base_audit_recomputes_exact_reconciliation_partition():
    rows = [
        {
            "interrupted_nonterminal_dispositions": 0,
            "reconciliation_path": None,
            "reconciliation_sha256": None,
        }
        for _ in range(35)
    ]
    rows.extend(
        {
            "interrupted_nonterminal_dispositions": 1,
            "reconciliation_path": f"work/unit-{index}/supervision_reconciliation.tsv",
            "reconciliation_sha256": "a" * 64,
        }
        for index in range(25)
    )
    assert postproduction._base_reconciliation_partition({"rows": rows}) == (35, 25)
    rows[0] = {
        "interrupted_nonterminal_dispositions": 1,
        "reconciliation_path": "work/changed/supervision_reconciliation.tsv",
        "reconciliation_sha256": "b" * 64,
    }
    with pytest.raises(ValueError, match="partition differs"):
        postproduction._base_reconciliation_partition({"rows": rows})


def test_inside_repo_rejects_an_in_repo_source_symlink(tmp_path):
    root = (tmp_path / "root").resolve()
    root.mkdir()
    target = root / "target.py"
    target.write_text("value = 1\n", encoding="utf-8")
    link = root / "link.py"
    try:
        link.symlink_to(target)
    except OSError:
        pytest.skip("filesystem does not permit symlink creation")
    with pytest.raises(ValueError, match="contains a symlink"):
        postproduction._inside_repo(link, root, label="source")


def test_immutable_authorization_rejects_different_payload(tmp_path):
    destination = tmp_path / "authorization.json"
    first = _signed(
        {
            "schema": postproduction.POSTPRODUCTION_AUTHORIZATION_SCHEMA,
            "status": "authorized",
            "value": 1,
        }
    )
    second = _signed({**first, "value": 2, "payload_sha256": "removed"})
    second.pop("payload_sha256", None)
    second["payload_sha256"] = postproduction._canonical_sha256(second)
    postproduction.write_postproduction_authorization(destination, first)
    with pytest.raises(FileExistsError, match="immutable"):
        postproduction.write_postproduction_authorization(destination, second)


def test_authorization_load_rejects_current_evidence_drift(tmp_path, monkeypatch):
    root = tmp_path.resolve()
    destination = root / postproduction.DEFAULT_AUTHORIZATION_RELATIVE_PATH
    destination.parent.mkdir(parents=True)
    payload = _signed(
        {
            "schema": postproduction.POSTPRODUCTION_AUTHORIZATION_SCHEMA,
            "status": "authorized",
            "selected_integration": {
                "path": postproduction.selected.DEFAULT_INTEGRATION_RELATIVE_PATH
            },
            "han_normalization_evidence": {
                "normalization": {
                    "path": postproduction.DEFAULT_NORMALIZATION_RELATIVE_PATH
                }
            },
            "value": 1,
        }
    )
    destination.write_text(json.dumps(payload), encoding="utf-8")
    current = _signed({key: value for key, value in payload.items() if key != "value"})
    monkeypatch.setattr(
        postproduction, "build_postproduction_authorization", lambda **_kwargs: current
    )
    with pytest.raises(ValueError, match="differs from current evidence"):
        postproduction.load_postproduction_authorization(
            destination,
            repo_root=root,
            slim_path=root / "slim",
        )


def test_parser_and_module_have_no_simulation_capable_action():
    assert postproduction.POSTPRODUCTION_ACTIONS == (
        "authorize",
        "validate",
        "aggregate-validate",
        "aggregate",
        "analyze",
        "report",
    )
    assert "simulate" not in postproduction._build_parser().format_help()
    assert "simulate_selected_units" not in postproduction.__dict__
    assert "run_recovery" not in postproduction.__dict__


@pytest.mark.parametrize(
    "action", ["aggregate-validate", "aggregate", "analyze", "report"]
)
def test_every_postproduction_route_delegates_inside_authorized_context(
    tmp_path, monkeypatch, action
):
    root = tmp_path.resolve()
    slim = root / ".native-stdpopsim/bin/slim"
    slim.parent.mkdir(parents=True)
    slim.write_bytes(b"slim")
    active = {"value": False}
    calls: list[str] = []

    @contextmanager
    def authorized_context(**_kwargs):
        active["value"] = True
        try:
            yield {"payload_sha256": "a" * 64}
        finally:
            active["value"] = False

    monkeypatch.setattr(
        postproduction, "authorized_postproduction_context", authorized_context
    )

    def require_context(label):
        assert active["value"] is True
        calls.append(label)

    def build_readiness(**_kwargs):
        require_context("aggregate-validate")
        return _signed(
            {
                "schema": postproduction.selected.AGGREGATE_READINESS_SCHEMA,
                "status": "ready_for_aggregate",
            }
        )

    def write_readiness(_path, _payload):
        assert active["value"] is True

    aggregate_paths = tuple(root / f"result-{index}.tsv" for index in range(5))
    monkeypatch.setattr(
        postproduction.selected, "build_aggregate_readiness", build_readiness
    )
    monkeypatch.setattr(
        postproduction.selected, "write_aggregate_readiness", write_readiness
    )
    monkeypatch.setattr(
        postproduction.selected,
        "aggregate_with_readiness",
        lambda **_kwargs: (require_context("aggregate"), aggregate_paths)[1],
    )
    monkeypatch.setattr(
        postproduction.selected,
        "analyze_with_readiness",
        lambda **_kwargs: (
            require_context("analyze"),
            {"status": "analyzed_and_reported"},
        )[1],
    )
    monkeypatch.setattr(
        postproduction.selected,
        "report_with_readiness",
        lambda **_kwargs: (require_context("report"), {"status": "reported"})[1],
    )

    assert (
        postproduction.main(
            [
                action,
                "--repo-root",
                str(root),
                "--slim-bin",
                str(slim),
            ]
        )
        == 0
    )
    assert calls == [action]
    assert active["value"] is False


def test_authorization_build_validates_evidence_before_entering_patch(
    tmp_path, monkeypatch
):
    root = tmp_path.resolve()
    slim = root / "slim"
    slim.write_bytes(b"slim")
    integration = root / postproduction.selected.DEFAULT_INTEGRATION_RELATIVE_PATH
    integration.parent.mkdir(parents=True)
    integration_payload = _signed({"schema": "integration", "status": "eas_ready"})
    integration.write_text(json.dumps(integration_payload), encoding="utf-8")
    normalization = root / postproduction.DEFAULT_NORMALIZATION_RELATIVE_PATH
    normalization.parent.mkdir(parents=True)
    normalization.write_text("{}", encoding="utf-8")
    for relative in (
        *postproduction.POSTPRODUCTION_SOURCE_PATHS,
        *postproduction.DELEGATED_SOURCE_PATHS,
    ):
        source = root / relative
        source.parent.mkdir(parents=True, exist_ok=True)
        source.write_text(f"# {relative}\n", encoding="utf-8")

    order: list[str] = []
    evidence = {
        "normalization": {"path": normalization.relative_to(root).as_posix()},
        "canonical_han_artifacts": {
            "recovery_authorization": {"sha256": "r" * 64},
            "recovery_completion_audit": {"sha256": "c" * 64},
            "base_all_units_audit": {"sha256": "b" * 64},
            "base_authorization": {"sha256": "a" * 64},
        },
    }

    def validate_evidence(**_kwargs):
        order.append("evidence")
        return evidence

    @contextmanager
    def patched():
        order.append("enter")
        yield
        order.append("exit")

    preflight = {
        "payload_sha256": "p" * 64,
        "counts": {"han_selected": 60},
        "han_v3": {
            "all_units_audit": {"sha256": "b" * 64},
            "final_operational_authorization": {"sha256": "a" * 64},
            "dtype_recovery": {
                "evidence_present": True,
                "binding": {
                    "authorization": {"sha256": "r" * 64},
                    "audit": {"sha256": "c" * 64},
                },
            },
        },
    }

    def build_readiness(**_kwargs):
        order.append("preflight")
        return preflight

    monkeypatch.setattr(
        postproduction,
        "_validate_normalization_evidence_for_bundle",
        validate_evidence,
    )
    monkeypatch.setattr(
        postproduction.recovery_audit,
        "patched_recovery_audit_aggregation",
        patched,
    )
    monkeypatch.setattr(
        postproduction.selected, "build_aggregate_readiness", build_readiness
    )
    payload = postproduction.build_postproduction_authorization(
        repo_root=root,
        slim_path=slim,
        integration_manifest_path=integration,
        normalization_path=normalization,
    )
    assert order == ["evidence", "enter", "preflight", "exit"]
    assert payload["simulation_route_present"] is False
