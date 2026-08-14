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
        "attest-eas",
        "supersede-legacy",
        "authorize",
        "validate",
        "aggregate-validate",
        "aggregate",
        "decode",
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
    eas_attestation_path = root / postproduction.DEFAULT_EAS_ATTESTATION_RELATIVE_PATH
    eas_attestation_path.parent.mkdir(parents=True, exist_ok=True)
    eas_attestation_path.write_text("{}", encoding="utf-8")
    legacy = _signed(
        {
            "schema": postproduction.LEGACY_POSTPRODUCTION_AUTHORIZATION_SCHEMA,
            "status": "authorized",
        }
    )
    legacy_path = (
        root
        / postproduction.DEFAULT_SUPERSEDED_AUTHORIZATION_RELATIVE_DIR
        / f"postproduction_authorization_v1__{legacy['payload_sha256']}.json"
    )
    legacy_path.parent.mkdir(parents=True)
    legacy_path.write_text(json.dumps(legacy), encoding="utf-8")
    for relative in (
        *postproduction.POSTPRODUCTION_SOURCE_PATHS,
        *postproduction.DELEGATED_SOURCE_PATHS,
    ):
        source = root / relative
        source.parent.mkdir(parents=True, exist_ok=True)
        source.write_text(f"# {relative}\n", encoding="utf-8")

    order: list[str] = []
    eas_attestation = {
        "payload_sha256": "e" * 64,
        "cache_audit": {
            "units": 60,
            "cache_validated_read_only": True,
            "variant_counts": {"canonical": 20, "transient_audit_text": 40},
        },
    }
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
        order.append("han-enter")
        yield
        order.append("han-exit")

    @contextmanager
    def patched_eas(*_args, **_kwargs):
        order.append("eas-enter")
        yield
        order.append("eas-exit")

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
    monkeypatch.setattr(
        postproduction,
        "load_eas_contract_attestation",
        lambda *_args, **_kwargs: (order.append("eas-audit"), eas_attestation)[1],
    )
    monkeypatch.setattr(
        postproduction,
        "patched_eas_contract_attestation_validation",
        patched_eas,
    )
    monkeypatch.setattr(
        postproduction,
        "_validate_eas_preflight_binding",
        lambda *_args, **_kwargs: order.append("eas-binding"),
    )
    payload = postproduction.build_postproduction_authorization(
        repo_root=root,
        slim_path=slim,
        integration_manifest_path=integration,
        normalization_path=normalization,
        eas_attestation_path=eas_attestation_path,
    )
    assert order == [
        "evidence",
        "eas-audit",
        "han-enter",
        "eas-enter",
        "preflight",
        "eas-exit",
        "han-exit",
        "eas-binding",
    ]
    assert payload["simulation_route_present"] is False


def _contract_with_attribute_pair(pair, *, science_value=1):
    return {
        "schema": "simulation-contract",
        "science": {"fixed_value": science_value},
        "implementation": {
            "sources": {
                postproduction.EAS_ATTRIBUTE_SOURCE_PATHS[0]: pair[0],
                postproduction.EAS_ATTRIBUTE_SOURCE_PATHS[1]: pair[1],
                "python/gamma_smc_aou/focused_selection_simulation.py": "s" * 64,
            }
        },
    }


def test_legacy_authorization_is_atomically_preserved_and_bound(tmp_path):
    root = tmp_path.resolve()
    canonical = root / postproduction.DEFAULT_AUTHORIZATION_RELATIVE_PATH
    canonical.parent.mkdir(parents=True)
    legacy = _signed(
        {
            "schema": postproduction.LEGACY_POSTPRODUCTION_AUTHORIZATION_SCHEMA,
            "status": "authorized",
            "evidence": "legacy-v1",
        }
    )
    content = (json.dumps(legacy, indent=2, sort_keys=True) + "\n").encode()
    canonical.write_bytes(content)
    record = postproduction.supersede_legacy_postproduction_authorization(
        repo_root=root
    )
    archive = root / record["path"]
    assert not canonical.exists()
    assert archive.read_bytes() == content
    assert record["payload_sha256"] == legacy["payload_sha256"]
    assert record["preservation_operation"] == "same_volume_atomic_rename"
    assert record["bytes_deleted"] is False
    assert (
        postproduction.supersede_legacy_postproduction_authorization(repo_root=root)
        == record
    )


def test_unknown_legacy_authorization_is_not_moved(tmp_path):
    root = tmp_path.resolve()
    canonical = root / postproduction.DEFAULT_AUTHORIZATION_RELATIVE_PATH
    canonical.parent.mkdir(parents=True)
    content = b'{"schema":"unknown","status":"authorized"}\n'
    canonical.write_bytes(content)
    with pytest.raises(ValueError, match="schema is unknown"):
        postproduction.supersede_legacy_postproduction_authorization(
            repo_root=root
        )
    assert canonical.read_bytes() == content
    archive_dir = root / postproduction.DEFAULT_SUPERSEDED_AUTHORIZATION_RELATIVE_DIR
    assert not archive_dir.exists()


def test_transient_attribute_evidence_is_exact_and_hash_bound(tmp_path):
    records = postproduction._transient_attribute_evidence_records(tmp_path)
    assert tuple(record["sha256"] for record in records.values()) == (
        postproduction.EAS_TRANSIENT_ATTRIBUTE_HASH_PAIR
    )
    for source_path, record in records.items():
        fixed = postproduction.EAS_TRANSIENT_ATTRIBUTE_EVIDENCE[source_path]
        content = postproduction.base64.b64decode(fixed["content_base64"])
        assert len(content) == record["size_bytes"]
        assert hashlib.sha256(content).hexdigest() == record["sha256"]
        assert content.endswith(b"focused_han_direct_v3_recovery_audit.py -text\n")


def test_contract_normalization_allows_only_exact_two_attribute_differences():
    expected = _contract_with_attribute_pair(
        postproduction.EAS_CANONICAL_ATTRIBUTE_HASH_PAIR
    )
    transient = _contract_with_attribute_pair(
        postproduction.EAS_TRANSIENT_ATTRIBUTE_HASH_PAIR
    )
    differences = postproduction._deep_difference_paths(transient, expected)
    assert set(differences) == set(postproduction.EAS_ATTRIBUTE_CONTRACT_PATHS)
    assert (
        postproduction._classify_eas_contract_variant(
            postproduction.EAS_TRANSIENT_ATTRIBUTE_HASH_PAIR, differences
        )
        == "transient_audit_text"
    )
    assert postproduction._normalized_eas_contract(transient, expected) == expected

    science_drift = _contract_with_attribute_pair(
        postproduction.EAS_TRANSIENT_ATTRIBUTE_HASH_PAIR, science_value=2
    )
    drift = postproduction._deep_difference_paths(science_drift, expected)
    with pytest.raises(ValueError, match="outside two attributes"):
        postproduction._classify_eas_contract_variant(
            postproduction.EAS_TRANSIENT_ATTRIBUTE_HASH_PAIR, drift
        )
    assert postproduction._normalized_eas_contract(science_drift, expected) != expected


@pytest.mark.parametrize(
    "pair",
    [
        ("f" * 64, "e" * 64),
        (
            postproduction.EAS_CANONICAL_ATTRIBUTE_HASH_PAIR[0],
            postproduction.EAS_TRANSIENT_ATTRIBUTE_HASH_PAIR[1],
        ),
        (
            postproduction.EAS_TRANSIENT_ATTRIBUTE_HASH_PAIR[0],
            postproduction.EAS_CANONICAL_ATTRIBUTE_HASH_PAIR[1],
        ),
    ],
)
def test_unknown_or_mixed_attribute_pairs_fail_closed(pair):
    with pytest.raises(ValueError, match="unknown or mixed"):
        postproduction._classify_eas_contract_variant(
            pair, postproduction.EAS_ATTRIBUTE_CONTRACT_PATHS
        )


def test_eas_attestation_payload_corruption_fails_closed():
    payload = _signed(
        {
            "schema": postproduction.EAS_CONTRACT_ATTESTATION_SCHEMA,
            "status": "passed",
            "cache_audit": {"units": 60},
        }
    )
    payload["cache_audit"]["units"] = 59
    with pytest.raises(ValueError, match="corrupt or incompatible"):
        postproduction._validate_hashed_payload(
            payload,
            schema=postproduction.EAS_CONTRACT_ATTESTATION_SCHEMA,
            status="passed",
            label="EAS attestation",
        )


def test_output_artifact_tamper_is_detected(tmp_path):
    unit = tmp_path / "work" / "unit"
    artifacts = postproduction.simulation._completion_artifacts(unit)
    paths = {
        "tree": artifacts.tree_path,
        "pair_table": artifacts.pair_table_path,
        "overall_pairs": artifacts.overall_pairs_path,
        "truth_profiles": artifacts.truth_profiles_path,
        "truth_class_summaries": artifacts.truth_class_summaries_path,
    }
    outputs = {}
    for index, (label, path) in enumerate(paths.items()):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(f"artifact-{index}\n".encode())
        outputs[label] = {
            "path": path.relative_to(unit).as_posix(),
            "sha256": han_v3.sha256_file(path),
        }
    completion = {"outputs": outputs}
    records = postproduction._attested_output_records(
        unit_dir=unit, completion=completion, repo_root=tmp_path
    )
    assert set(records) == set(paths)
    artifacts.tree_path.write_bytes(b"tampered\n")
    with pytest.raises(ValueError, match="checksum differs: tree"):
        postproduction._attested_output_records(
            unit_dir=unit, completion=completion, repo_root=tmp_path
        )


def test_standard_raw_validator_fails_but_scoped_attestation_patch_succeeds_and_restores(
    tmp_path, monkeypatch
):
    root = tmp_path.resolve()
    integration = root / "integration.json"
    units = root / "execution_units.tsv"
    attestation_path = root / postproduction.DEFAULT_EAS_ATTESTATION_RELATIVE_PATH
    attestation_path.parent.mkdir(parents=True)
    integration.write_text("integration\n", encoding="utf-8")
    units.write_text("units\n", encoding="utf-8")
    attestation_path.write_text("attestation\n", encoding="utf-8")
    payload = {
        "payload_sha256": "a" * 64,
        "selected_integration": {
            **postproduction._file_record(integration, root),
            "payload_sha256": "i" * 64,
        },
        "execution_units": postproduction._file_record(units, root),
        "cache_audit": {"units": 60, "cache_validated_read_only": True},
    }
    bundle = SimpleNamespace(
        integration_manifest_path=integration,
        execution_units_path=units,
    )

    def raw_failure(_bundle):
        raise ValueError("raw current-contract mismatch")

    monkeypatch.setattr(
        postproduction.selected,
        "validate_eas_cached_units_read_only",
        raw_failure,
    )
    with pytest.raises(ValueError, match="raw current-contract mismatch"):
        postproduction.selected.validate_eas_cached_units_read_only(bundle)
    with postproduction.patched_eas_contract_attestation_validation(
        payload, path=attestation_path, repo_root=root
    ):
        patched = postproduction.selected.validate_eas_cached_units_read_only(bundle)
        assert patched["units"] == 60
        assert patched["contract_attestation"]["payload_sha256"] == "a" * 64
    assert postproduction.selected.validate_eas_cached_units_read_only is raw_failure


def test_decode_delegates_exact_720_defaults_inside_authorized_context(
    tmp_path, monkeypatch
):
    root = tmp_path.resolve()
    slim = root / ".native-stdpopsim/bin/slim"
    decoder = root / postproduction.DEFAULT_DECODER_RELATIVE_PATH
    slim.parent.mkdir(parents=True)
    decoder.parent.mkdir(parents=True)
    slim.write_bytes(b"slim")
    decoder.write_bytes(b"decoder")
    active = {"value": False}
    observed = {}

    @contextmanager
    def authorized(**_kwargs):
        active["value"] = True
        try:
            yield {"payload_sha256": "p" * 64}
        finally:
            active["value"] = False

    def decode_main(argv):
        assert active["value"] is True
        assert argv[0] == "decode"
        assert "simulate" not in argv
        assert "--remove-raw-after-success" not in argv
        observed["argv"] = argv
        results = (
            root
            / postproduction.selected.DEFAULT_CAMPAIGN_RELATIVE_PATH
            / "results"
        )
        results.mkdir(parents=True)
        pd.DataFrame(
            {"unit_id": [f"unit-{index:03d}" for index in range(720)], "status": "cached"}
        ).to_csv(results / "decode_status.tsv", sep="\t", index=False)
        return 0

    monkeypatch.setattr(postproduction, "authorized_postproduction_context", authorized)
    monkeypatch.setattr(postproduction.campaign, "main", decode_main)
    assert (
        postproduction.main(
            ["decode", "--repo-root", str(root), "--slim-bin", str(slim)]
        )
        == 0
    )
    assert active["value"] is False
    argv = observed["argv"]
    assert argv[argv.index("--workers") + 1] == "20"
    assert argv[argv.index("--threads") + 1] == "1"
