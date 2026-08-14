from __future__ import annotations

import json
import shutil
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
from gamma_smc_aou import focused_han_direct_v3 as han_v3
from gamma_smc_aou import focused_selected_postdecode as postdecode
from gamma_smc_aou import focused_selection_campaign as campaign
from gamma_smc_aou import focused_selection_simulation as simulation


def _signed(payload: dict) -> dict:
    result = dict(payload)
    result["payload_sha256"] = postdecode.postproduction._canonical_sha256(result)
    return result


def _unit(unit_id: str, *, demography: str = "eas_phlash_median") -> dict:
    if demography == "eas_phlash_median":
        population = "EAS"
        kind = "no_introgression"
        source = "PHLASH EAS pointwise median trajectory"
        origin = "de_novo"
    else:
        population = "Han"
        kind = "introgression"
        source = "stdpopsim AncientEurasia_9K19"
        origin = "archaic_specific_introgressed_standing_variation"
    return {
        "unit_id": unit_id,
        "demography_id": demography,
        "demography_kind": kind,
        "population": population,
        "source_model": source,
        "selection_origin": origin,
        "simulation_class": "selected",
        "selection_coefficient": 0.005,
        "target_allele_frequency": 0.1,
        "population_af_lower": 0.075,
        "population_af_upper": 0.125,
        "exact_sample_alt_count": 20,
        "sample_diploids": 100,
        "replicate_index": 1,
        "seed": 12345,
        "sequence_length_bp": 10_000_000,
        "focal_position_bp": 5_000_000,
    }


def _contract(unit: dict, pair: tuple[str, str], *, science: int = 1) -> dict:
    return {
        "schema": simulation.SCHEMA_VERSION,
        "unit": dict(unit),
        "parameters": {"science": science},
        "implementation": {
            "sources": {
                postdecode.postproduction.EAS_ATTRIBUTE_SOURCE_PATHS[0]: pair[0],
                postdecode.postproduction.EAS_ATTRIBUTE_SOURCE_PATHS[1]: pair[1],
                "python/gamma_smc_aou/focused_selection_simulation.py": "f" * 64,
            }
        },
    }


def _dummy_file_record(path: str) -> dict:
    return {"path": path, "sha256": "d" * 64, "size_bytes": 1}


def _row(
    unit_id: str,
    *,
    variant: str,
    stored_sha256: str = "1" * 64,
    current_sha256: str = "2" * 64,
) -> dict:
    pair = (
        postdecode.postproduction.EAS_CANONICAL_ATTRIBUTE_HASH_PAIR
        if variant == "canonical"
        else postdecode.postproduction.EAS_TRANSIENT_ATTRIBUTE_HASH_PAIR
    )
    differences = (
        []
        if variant == "canonical"
        else [
            list(path)
            for path in sorted(postdecode.postproduction.EAS_ATTRIBUTE_CONTRACT_PATHS)
        ]
    )
    unit_root = f"{postdecode.selected.DEFAULT_CAMPAIGN_RELATIVE_PATH}/work/{unit_id}"
    outputs = {
        "tree": _dummy_file_record(f"{unit_root}/simulation.trees"),
        "pair_table": _dummy_file_record(f"{unit_root}/sample_manifest.tsv"),
        "overall_pairs": _dummy_file_record(f"{unit_root}/pairs/overall.pairs.tsv"),
        "truth_profiles": _dummy_file_record(f"{unit_root}/truth_profiles.tsv.gz"),
        "truth_class_summaries": _dummy_file_record(
            f"{unit_root}/truth_class_summaries.tsv"
        ),
    }
    return {
        "unit_id": unit_id,
        "contract_variant": variant,
        "completion": _dummy_file_record(f"{unit_root}/simulation_complete.json"),
        "contract_file": _dummy_file_record(f"{unit_root}/simulation_contract.json"),
        "stored_contract_sha256": stored_sha256,
        "normalized_contract_sha256": current_sha256,
        "current_expected_contract_sha256": current_sha256,
        "attribute_source_hashes": dict(
            zip(
                postdecode.postproduction.EAS_ATTRIBUTE_SOURCE_PATHS,
                pair,
                strict=True,
            )
        ),
        "contract_difference_paths": differences,
        "outputs": outputs,
        "outputs_canonical_sha256": postdecode.postproduction._canonical_sha256(
            outputs
        ),
        "base_validation_against_stored_contract": True,
        "current_contract_rebuilt_without_simulation": True,
        "normalized_equality_except_allowed_attribute_hashes": True,
        "attempt_provenance": {"validated": True},
    }


def _attestation(rows: list[dict]) -> dict:
    audit = {
        "units": 60,
        "cells": 6,
        "replicates_per_cell": 10,
        "cache_validated_read_only": True,
        "base_validated_against_each_stored_contract": True,
        "current_contract_independently_rebuilt_for_each_unit": True,
        "normalized_differences_limited_to_two_attribute_hashes": True,
        "simulation_launched": False,
        "artifacts_changed": False,
        "variant_counts": {"canonical": 20, "transient_audit_text": 40},
        "rows_canonical_sha256": postdecode.postproduction._canonical_sha256(rows),
        "rows": rows,
    }
    return _signed(
        {
            "schema": postdecode.postproduction.EAS_CONTRACT_ATTESTATION_SCHEMA,
            "status": "passed",
            "cache_audit": audit,
        }
    )


def _write_attestation(root: Path, payload: dict) -> Path:
    path = root / postdecode.postproduction.DEFAULT_EAS_ATTESTATION_RELATIVE_PATH
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path


def _attested_fixture(
    root: Path,
    *,
    stored_pair: tuple[str, str] | None = None,
) -> SimpleNamespace:
    unit_id = "eas_phlash_median__af10__selected_s0p005__rep001"
    unit = _unit(unit_id)
    current = _contract(
        unit, postdecode.postproduction.EAS_CANONICAL_ATTRIBUTE_HASH_PAIR
    )
    stored = _contract(
        unit,
        stored_pair or postdecode.postproduction.EAS_TRANSIENT_ATTRIBUTE_HASH_PAIR,
    )
    unit_dir = (
        root / postdecode.selected.DEFAULT_CAMPAIGN_RELATIVE_PATH / "work" / unit_id
    )
    artifacts = simulation._completion_artifacts(unit_dir)
    output_paths = {
        "tree": artifacts.tree_path,
        "pair_table": artifacts.pair_table_path,
        "overall_pairs": artifacts.overall_pairs_path,
        "truth_profiles": artifacts.truth_profiles_path,
        "truth_class_summaries": artifacts.truth_class_summaries_path,
    }
    outputs = {}
    for index, (label, path) in enumerate(output_paths.items()):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(f"output-{index}\n".encode())
        outputs[label] = postdecode.postproduction._file_record(path, root)
    contract_path = unit_dir / "simulation_contract.json"
    contract_path.write_text(json.dumps(stored), encoding="utf-8")
    completion = {
        "schema": simulation.SCHEMA_VERSION,
        "status": "complete",
        "contract": stored,
        "contract_sha256": postdecode.postproduction._canonical_sha256(stored),
        "outputs": {
            label: {
                "path": path.relative_to(unit_dir).as_posix(),
                "sha256": outputs[label]["sha256"],
            }
            for label, path in output_paths.items()
        },
    }
    artifacts.completion_path.write_text(json.dumps(completion), encoding="utf-8")
    active = _row(
        unit_id,
        variant="transient_audit_text",
        stored_sha256=postdecode.postproduction._canonical_sha256(stored),
        current_sha256=postdecode.postproduction._canonical_sha256(current),
    )
    active["completion"] = postdecode.postproduction._file_record(
        artifacts.completion_path, root
    )
    active["contract_file"] = postdecode.postproduction._file_record(
        contract_path, root
    )
    active["outputs"] = outputs
    active["outputs_canonical_sha256"] = postdecode.postproduction._canonical_sha256(
        outputs
    )
    rows = [active]
    rows.extend(
        _row(
            f"eas_phlash_median__dummy_canonical__rep{index:03d}",
            variant="canonical",
        )
        for index in range(1, 21)
    )
    rows.extend(
        _row(
            f"eas_phlash_median__dummy_transient__rep{index:03d}",
            variant="transient_audit_text",
        )
        for index in range(1, 40)
    )
    payload = _attestation(rows)
    path = _write_attestation(root, payload)
    return SimpleNamespace(
        unit=unit,
        current=current,
        stored=stored,
        artifacts=artifacts,
        payload=payload,
        path=path,
        rows=rows,
    )


def test_action_surface_has_only_authorize_and_postdecode_routes():
    assert postdecode.POSTDECODE_ACTIONS == (
        "authorize",
        "aggregate",
        "analyze",
        "report",
    )
    for forbidden in ("simulate", "decode", "recovery"):
        assert forbidden not in postdecode.POSTDECODE_ACTIONS
        assert forbidden not in postdecode.__dict__


def test_attested_eas_validator_uses_stored_contract_and_restores(
    tmp_path, monkeypatch
):
    root = tmp_path.resolve()
    fixture = _attested_fixture(root)
    calls = []

    def base_validator(artifacts, contract):
        calls.append(contract)
        assert artifacts is fixture.artifacts
        assert contract == fixture.stored
        return simulation.SimulationArtifacts(
            **{**artifacts.__dict__, "cache_hit": True}
        )

    monkeypatch.setattr(simulation, "_validate_completion", base_validator)
    with postdecode.patched_attested_eas_aggregate_validation(
        fixture.payload,
        path=fixture.path,
        repo_root=root,
        require_complete_coverage=False,
    ) as coverage:
        result = simulation._validate_completion(fixture.artifacts, fixture.current)
        assert result.cache_hit is True
    assert coverage["validated_unit_ids"] == [fixture.unit["unit_id"]]
    assert calls == [fixture.stored]
    assert simulation._validate_completion is base_validator


def test_non_eas_validation_delegates_byte_for_byte(tmp_path, monkeypatch):
    root = tmp_path.resolve()
    fixture = _attested_fixture(root)
    sentinel = object()
    neutral = {
        "unit": {"simulation_class": "neutral", "demography_id": "eas_phlash_median"}
    }

    def base_validator(artifacts, contract):
        assert artifacts is sentinel
        assert contract is neutral
        return sentinel

    monkeypatch.setattr(simulation, "_validate_completion", base_validator)
    with postdecode.patched_attested_eas_aggregate_validation(
        fixture.payload,
        path=fixture.path,
        repo_root=root,
        require_complete_coverage=False,
    ):
        assert simulation._validate_completion(sentinel, neutral) is sentinel
    assert simulation._validate_completion is base_validator


def test_current_science_drift_fails_closed_and_restores(tmp_path, monkeypatch):
    root = tmp_path.resolve()
    fixture = _attested_fixture(root)

    def base_validator(*_args):
        raise AssertionError("base validator must not receive unattested science drift")

    monkeypatch.setattr(simulation, "_validate_completion", base_validator)
    drift = _contract(
        fixture.unit,
        postdecode.postproduction.EAS_CANONICAL_ATTRIBUTE_HASH_PAIR,
        science=2,
    )
    with (
        pytest.raises(ValueError, match="current EAS aggregate contract differs"),
        postdecode.patched_attested_eas_aggregate_validation(
            fixture.payload,
            path=fixture.path,
            repo_root=root,
            require_complete_coverage=False,
        ),
    ):
        simulation._validate_completion(fixture.artifacts, drift)
    assert simulation._validate_completion is base_validator


def test_output_tamper_fails_before_base_validator(tmp_path, monkeypatch):
    root = tmp_path.resolve()
    fixture = _attested_fixture(root)

    def base_validator(*_args):
        raise AssertionError("base validator must not receive tampered output")

    monkeypatch.setattr(simulation, "_validate_completion", base_validator)
    fixture.artifacts.tree_path.write_bytes(b"tampered\n")
    with (
        pytest.raises(ValueError, match="checksum differs: tree"),
        postdecode.patched_attested_eas_aggregate_validation(
            fixture.payload,
            path=fixture.path,
            repo_root=root,
            require_complete_coverage=False,
        ),
    ):
        simulation._validate_completion(fixture.artifacts, fixture.current)


def test_redirected_unit_directory_fails_before_base_validator(tmp_path, monkeypatch):
    root = tmp_path.resolve()
    fixture = _attested_fixture(root)
    relocated = root / "redirected-unit"
    shutil.move(str(fixture.artifacts.unit_dir), relocated)
    try:
        fixture.artifacts.unit_dir.symlink_to(relocated, target_is_directory=True)
    except OSError:
        pytest.skip("filesystem does not permit directory symlinks")
    monkeypatch.setattr(
        simulation,
        "_validate_completion",
        lambda *_args: pytest.fail("redirected unit reached base validator"),
    )
    with (
        pytest.raises(ValueError, match="path contains a symlink"),
        postdecode.patched_attested_eas_aggregate_validation(
            fixture.payload,
            path=fixture.path,
            repo_root=root,
            require_complete_coverage=False,
        ),
    ):
        simulation._validate_completion(fixture.artifacts, fixture.current)


def test_mixed_attribute_pair_fails_closed(tmp_path, monkeypatch):
    root = tmp_path.resolve()
    mixed = (
        postdecode.postproduction.EAS_CANONICAL_ATTRIBUTE_HASH_PAIR[0],
        postdecode.postproduction.EAS_TRANSIENT_ATTRIBUTE_HASH_PAIR[1],
    )
    fixture = _attested_fixture(root, stored_pair=mixed)
    monkeypatch.setattr(
        simulation,
        "_validate_completion",
        lambda *_args: pytest.fail("mixed pair reached base validator"),
    )
    with (
        pytest.raises(ValueError, match="unknown or mixed"),
        postdecode.patched_attested_eas_aggregate_validation(
            fixture.payload,
            path=fixture.path,
            repo_root=root,
            require_complete_coverage=False,
        ),
    ):
        simulation._validate_completion(fixture.artifacts, fixture.current)


def test_corrupt_attestation_fails_before_patch(tmp_path, monkeypatch):
    root = tmp_path.resolve()
    fixture = _attested_fixture(root)
    fixture.payload["cache_audit"]["units"] = 59
    _write_attestation(root, fixture.payload)
    original = simulation._validate_completion
    with (
        pytest.raises(ValueError, match="corrupt or incompatible"),
        postdecode.patched_attested_eas_aggregate_validation(
            fixture.payload,
            path=fixture.path,
            repo_root=root,
            require_complete_coverage=False,
        ),
    ):
        pytest.fail("corrupt attestation entered patch")
    assert simulation._validate_completion is original


def test_exact_60_coverage_is_required_only_on_success(tmp_path, monkeypatch):
    root = tmp_path.resolve()
    fixture = _attested_fixture(root)
    original = simulation._validate_completion
    monkeypatch.setattr(
        postdecode,
        "_validate_attested_eas_completion",
        lambda artifacts, *_args, **_kwargs: artifacts,
    )
    with postdecode.patched_attested_eas_aggregate_validation(
        fixture.payload,
        path=fixture.path,
        repo_root=root,
        require_complete_coverage=True,
    ):
        for row in fixture.rows:
            artifacts = SimpleNamespace(cache_hit=True)
            simulation._validate_completion(
                artifacts,
                {
                    "unit": {
                        "unit_id": row["unit_id"],
                        "simulation_class": "selected",
                        "demography_id": "eas_phlash_median",
                    }
                },
            )
    assert simulation._validate_completion is original

    with (
        pytest.raises(ValueError, match="every EAS attested unit exactly once"),
        postdecode.patched_attested_eas_aggregate_validation(
            fixture.payload,
            path=fixture.path,
            repo_root=root,
            require_complete_coverage=True,
        ),
    ):
        for row in fixture.rows[:-1]:
            simulation._validate_completion(
                SimpleNamespace(cache_hit=True),
                {
                    "unit": {
                        "unit_id": row["unit_id"],
                        "simulation_class": "selected",
                        "demography_id": "eas_phlash_median",
                    }
                },
            )
    assert simulation._validate_completion is original

    with postdecode.patched_attested_eas_aggregate_validation(
        fixture.payload,
        path=fixture.path,
        repo_root=root,
        require_complete_coverage=False,
    ):
        pass


def _route_evidence(root: Path, fixture: SimpleNamespace) -> SimpleNamespace:
    for relative in postdecode.ADAPTER_SOURCE_PATHS:
        source = root / relative
        source.parent.mkdir(parents=True, exist_ok=True)
        source.write_text(f"# {relative}\n", encoding="utf-8")
    integration = root / postdecode.selected.DEFAULT_INTEGRATION_RELATIVE_PATH
    integration.parent.mkdir(parents=True, exist_ok=True)
    integration.write_text("integration\n", encoding="utf-8")
    integration_record = postdecode.postproduction._file_record(integration, root)
    attestation_record = {
        **postdecode.postproduction._file_record(fixture.path, root),
        "payload_sha256": fixture.payload["payload_sha256"],
    }
    readiness = _signed(
        {
            "schema": postdecode.selected.AGGREGATE_READINESS_SCHEMA,
            "status": "ready_for_aggregate",
            "eas": {
                "integration_manifest": integration_record,
                "integration_payload_sha256": "i" * 64,
                "cache_audit": {
                    **fixture.payload["cache_audit"],
                    "contract_attestation": attestation_record,
                },
            },
        }
    )
    readiness_path = (
        root / postdecode.selected.DEFAULT_AGGREGATE_READINESS_RELATIVE_PATH
    )
    readiness_path.write_text(json.dumps(readiness), encoding="utf-8")
    upstream = _signed(
        {
            "schema": postdecode.postproduction.POSTPRODUCTION_AUTHORIZATION_SCHEMA,
            "status": "authorized",
            "selected_integration": {
                **integration_record,
                "payload_sha256": "i" * 64,
                "strict_loaded": True,
            },
            "eas_contract_attestation": {
                **attestation_record,
                "validated_units": 60,
                "variant_counts": {
                    "canonical": 20,
                    "transient_audit_text": 40,
                },
            },
            "preflight_payload_sha256": readiness["payload_sha256"],
            "every_postproduction_route_inside_both_contexts": True,
            "eas_cache_validator_patched_only_inside_postproduction_context": True,
            "actions": ["aggregate", "analyze", "report"],
        }
    )
    upstream_path = root / postdecode.postproduction.DEFAULT_AUTHORIZATION_RELATIVE_PATH
    upstream_path.write_text(json.dumps(upstream), encoding="utf-8")
    return SimpleNamespace(
        upstream=upstream,
        upstream_path=upstream_path,
        readiness=readiness,
        readiness_path=readiness_path,
        integration=integration,
        attestation_record=attestation_record,
    )


def test_immutable_authorization_binds_sources_and_upstream(tmp_path):
    root = tmp_path.resolve()
    fixture = _attested_fixture(root)
    evidence = _route_evidence(root, fixture)
    payload = postdecode.build_postdecode_authorization(
        repo_root=root,
        upstream_authorization=evidence.upstream,
        readiness=evidence.readiness,
        attestation=fixture.payload,
        upstream_authorization_path=evidence.upstream_path,
        readiness_path=evidence.readiness_path,
        attestation_path=fixture.path,
    )
    destination = root / postdecode.DEFAULT_POSTDECODE_AUTHORIZATION_RELATIVE_PATH
    postdecode.write_postdecode_authorization(destination, payload, repo_root=root)
    assert payload["adapter_sources"] == (
        postdecode.postproduction._canonical_source_records(
            root, postdecode.ADAPTER_SOURCE_PATHS
        )
    )
    assert payload["upstream_artifacts"]["aggregate_readiness"] == {
        **postdecode.postproduction._file_record(evidence.readiness_path, root),
        "payload_sha256": evidence.readiness["payload_sha256"],
    }
    assert (
        postdecode.load_postdecode_authorization(
            destination, repo_root=root, expected=payload
        )
        == payload
    )
    source = root / postdecode.MODULE_SOURCE_PATH
    source.write_text("# tampered\n", encoding="utf-8")
    drift = postdecode.build_postdecode_authorization(
        repo_root=root,
        upstream_authorization=evidence.upstream,
        readiness=evidence.readiness,
        attestation=fixture.payload,
        upstream_authorization_path=evidence.upstream_path,
        readiness_path=evidence.readiness_path,
        attestation_path=fixture.path,
    )
    with pytest.raises(ValueError, match="differs from current sources/evidence"):
        postdecode.load_postdecode_authorization(
            destination, repo_root=root, expected=drift
        )
    with pytest.raises(FileExistsError, match="immutable"):
        postdecode.write_postdecode_authorization(destination, drift, repo_root=root)


def test_actual_han_and_eas_dispatch_reaches_nested_aggregate_with_both_contexts(
    tmp_path, monkeypatch
):
    root = tmp_path.resolve()
    fixture = _attested_fixture(root)
    campaign_dir = root / postdecode.selected.DEFAULT_CAMPAIGN_RELATIVE_PATH
    bundle = SimpleNamespace(
        repo_root=root,
        campaign_dir=campaign_dir,
        slim_path=root / "slim",
    )
    han = _unit(
        "ancient_eurasia_han_introgression__af10__selected_s0p005__rep001",
        demography="ancient_eurasia_han_introgression",
    )
    units = pd.DataFrame([han, fixture.unit])
    raw_han = pd.DataFrame(
        [
            {
                "unit_id": han["unit_id"],
                "interrupted_nonterminal_dispositions": 0,
                "reconciliation_path": None,
                "reconciliation_sha256": None,
            }
        ]
    )
    raw_aggregate = lambda *_args, **_kwargs: raw_han.copy()
    monkeypatch.setattr(han_v3, "aggregate_v3_cached_completions", raw_aggregate)
    monkeypatch.setattr(
        han_v3,
        "validate_v3_cached_unit",
        lambda *_args, **_kwargs: {"cache_validated_read_only": True},
    )

    @contextmanager
    def selected_context(_bundle):
        yield

    monkeypatch.setattr(
        postdecode.selected,
        "load_selected_production_bundle",
        lambda *_args, **_kwargs: object(),
    )
    monkeypatch.setattr(
        postdecode.selected, "patched_selected_executor", selected_context
    )
    base_calls = []

    def base_validator(artifacts, contract):
        base_calls.append(contract)
        assert contract == fixture.stored
        return simulation.SimulationArtifacts(
            **{**artifacts.__dict__, "cache_hit": True}
        )

    monkeypatch.setattr(simulation, "_validate_completion", base_validator)

    def delegated_simulation(unit, *_args, **_kwargs):
        assert unit["unit_id"] == fixture.unit["unit_id"]
        return simulation._validate_completion(fixture.artifacts, fixture.current)

    monkeypatch.setattr(simulation, "_simulate_unit_unlocked", delegated_simulation)
    outputs = tuple(root / f"aggregate-{index}.tsv" for index in range(5))

    def nested_aggregate(observed_campaign, observed_units):
        assert observed_campaign == campaign_dir
        normalized = han_v3.aggregate_v3_cached_completions(bundle)
        assert normalized.loc[0, "reconciliation_path"] is None
        artifacts = [
            simulation._simulate_unit_unlocked(unit, root, campaign_dir)
            for unit in observed_units.to_dict(orient="records")
        ]
        assert all(result.cache_hit is True for result in artifacts)
        return outputs

    monkeypatch.setattr(campaign, "_aggregate_compact_decode", nested_aggregate)
    with (
        postdecode.postproduction.recovery_audit.patched_recovery_audit_aggregation(),
        postdecode.patched_attested_eas_aggregate_validation(
            fixture.payload,
            path=fixture.path,
            repo_root=root,
            require_complete_coverage=False,
        ),
    ):
        observed = han_v3.aggregate_decoded_outputs_with_v3_dispatch(
            bundle,
            units,
            selected_integration_manifest=root / "integration.json",
        )
    assert observed == outputs
    assert base_calls == [fixture.stored]
    assert simulation._validate_completion is base_validator
    assert han_v3.aggregate_v3_cached_completions is raw_aggregate


def test_both_global_contexts_restore_after_nested_error(tmp_path, monkeypatch):
    root = tmp_path.resolve()
    fixture = _attested_fixture(root)
    raw_han = pd.DataFrame(
        [
            {
                "unit_id": "han",
                "interrupted_nonterminal_dispositions": 0,
                "reconciliation_path": None,
                "reconciliation_sha256": None,
            }
        ]
    )
    raw_aggregate = lambda *_args, **_kwargs: raw_han.copy()
    base_validator = simulation._validate_completion
    monkeypatch.setattr(han_v3, "aggregate_v3_cached_completions", raw_aggregate)
    with (
        pytest.raises(RuntimeError, match="injected nested aggregate failure"),
        postdecode.postproduction.recovery_audit.patched_recovery_audit_aggregation(),
        postdecode.patched_attested_eas_aggregate_validation(
            fixture.payload,
            path=fixture.path,
            repo_root=root,
            require_complete_coverage=True,
        ),
    ):
        assert simulation._validate_completion is not base_validator
        assert han_v3.aggregate_v3_cached_completions is not raw_aggregate
        raise RuntimeError("injected nested aggregate failure")
    assert simulation._validate_completion is base_validator
    assert han_v3.aggregate_v3_cached_completions is raw_aggregate


@pytest.mark.parametrize("action", ["aggregate", "analyze", "report"])
def test_every_execution_route_keeps_contexts_and_restores_on_error(
    tmp_path, monkeypatch, action
):
    root = tmp_path.resolve()
    fixture = _attested_fixture(root)
    evidence = _route_evidence(root, fixture)
    slim = root / ".native-stdpopsim/bin/slim"
    slim.parent.mkdir(parents=True)
    slim.write_bytes(b"slim")
    own_path = root / postdecode.DEFAULT_POSTDECODE_AUTHORIZATION_RELATIVE_PATH
    own_path.write_text("{}\n", encoding="utf-8")
    own = _signed(
        {
            "schema": postdecode.POSTDECODE_AUTHORIZATION_SCHEMA,
            "status": "authorized",
        }
    )
    raw_han = pd.DataFrame(
        [
            {
                "unit_id": "han",
                "interrupted_nonterminal_dispositions": 0,
                "reconciliation_path": None,
                "reconciliation_sha256": None,
            }
        ]
    )
    raw_aggregate = lambda *_args, **_kwargs: raw_han.copy()
    monkeypatch.setattr(han_v3, "aggregate_v3_cached_completions", raw_aggregate)
    base_validator = simulation._validate_completion

    @contextmanager
    def upstream_context(**_kwargs):
        with postdecode.postproduction.recovery_audit.patched_recovery_audit_aggregation():
            yield evidence.upstream

    monkeypatch.setattr(
        postdecode.postproduction,
        "authorized_postproduction_context",
        upstream_context,
    )
    monkeypatch.setattr(
        postdecode.postproduction,
        "load_eas_contract_attestation",
        lambda *_args, **_kwargs: fixture.payload,
    )
    monkeypatch.setattr(
        postdecode.selected,
        "load_aggregate_readiness",
        lambda *_args, **_kwargs: evidence.readiness,
    )
    monkeypatch.setattr(postdecode, "_route_binding", lambda *_args, **_kwargs: {})
    monkeypatch.setattr(
        postdecode,
        "build_postdecode_authorization",
        lambda **_kwargs: own,
    )
    monkeypatch.setattr(
        postdecode,
        "load_postdecode_authorization",
        lambda *_args, **_kwargs: own,
    )
    monkeypatch.setattr(
        postdecode,
        "_validate_attested_eas_completion",
        lambda artifacts, *_args, **_kwargs: artifacts,
    )

    def route(**_kwargs):
        assert simulation._validate_completion is not base_validator
        normalized = han_v3.aggregate_v3_cached_completions(object())
        assert normalized.loc[0, "reconciliation_path"] is None
        if action != "report":
            for row in fixture.rows:
                simulation._validate_completion(
                    SimpleNamespace(cache_hit=True),
                    {
                        "unit": {
                            "unit_id": row["unit_id"],
                            "simulation_class": "selected",
                            "demography_id": "eas_phlash_median",
                        }
                    },
                )
        if action == "aggregate":
            return tuple(root / f"result-{index}.tsv" for index in range(5))
        return {"status": action}

    monkeypatch.setattr(postdecode.selected, "aggregate_with_readiness", route)
    monkeypatch.setattr(postdecode.selected, "analyze_with_readiness", route)
    monkeypatch.setattr(postdecode.selected, "report_with_readiness", route)
    assert (
        postdecode.main([action, "--repo-root", str(root), "--slim-bin", str(slim)])
        == 0
    )
    assert simulation._validate_completion is base_validator
    assert han_v3.aggregate_v3_cached_completions is raw_aggregate
