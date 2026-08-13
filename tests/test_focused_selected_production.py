from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
from gamma_smc_aou import focused_selected_production as production
from gamma_smc_aou import focused_selection_campaign as campaign
from gamma_smc_aou import focused_selection_simulation as simulation


def _execution_units(root: Path) -> pd.DataFrame:
    campaign_dir = root / "focused_selection_EAS_sim"
    plan = campaign.FocusedCampaignPlan(
        repo_root=root,
        campaign_dir=campaign_dir,
    )
    units = campaign.build_execution_units(plan)
    campaign_dir.mkdir(parents=True, exist_ok=True)
    units.to_csv(
        campaign_dir / "execution_units.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    return units


def _selected_unit(demography_id: str = "eas_phlash_median") -> dict[str, object]:
    population = "EAS" if demography_id == "eas_phlash_median" else "Han"
    return {
        "unit_id": f"test-{demography_id}",
        "demography_id": demography_id,
        "demography_kind": (
            "no_introgression"
            if demography_id == "eas_phlash_median"
            else "introgression"
        ),
        "population": population,
        "source_model": "test",
        "selection_origin": "test",
        "simulation_class": "selected",
        "selection_coefficient": 0.01,
        "target_allele_frequency": 0.10,
        "population_af_lower": 0.075,
        "population_af_upper": 0.125,
        "exact_sample_alt_count": 20,
        "sample_diploids": 100,
        "replicate_index": 1,
        "seed": 41,
        "sequence_length_bp": 10_000_000,
        "focal_position_bp": 5_000_000,
    }


def _bundle(root: Path) -> production.SelectedProductionBundle:
    manifest_path = root / production.DEFAULT_INTEGRATION_RELATIVE_PATH
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text("{}\n", encoding="utf-8")
    slim = root / "bin/slim"
    slim.parent.mkdir(parents=True, exist_ok=True)
    slim.write_bytes(b"slim")
    origin = {"age_generations": 999.0}
    origin_hash = production.eas_calibration._canonical_sha256(origin)
    return production.SelectedProductionBundle(
        repo_root=root,
        campaign_dir=root / "focused_selection_EAS_sim",
        integration_manifest_path=manifest_path,
        slim_path=slim,
        execution_units_path=root / "focused_selection_EAS_sim/execution_units.tsv",
        manifest={"payload_sha256": "integration"},
        eas_authorization={"payload_sha256": "adapter"},
        eas_frozen={
            "payload_sha256": "eas-frozen",
            "cell_authorizations": {
                "s=0.010000|af=0.100000": {
                    "origin_age_generations": 999.0,
                    "origin_contract_sha256": origin_hash,
                    "event_contract_sha256": "events",
                }
            },
        },
        han_auxiliary={"manifest_sha256": "han-manifest"},
        han_binding={"schema": "han-binding", "mapping": {}},
        implementation_paths=production.INTEGRATION_IMPLEMENTATION_PATHS,
    )


def test_byte_bound_sources_are_preserved_across_platforms():
    root = Path(__file__).resolve().parents[1]
    rules = {
        line.split()[0]: line.split()[1:]
        for line in (root / ".gitattributes").read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    }
    for relative in (production.MODULE_SOURCE_PATH, production.WRAPPER_SOURCE_PATH):
        assert "-text" in rules.get(relative, [])


def _fake_authorizations():
    eas_authorization = {"payload_sha256": "eas-adapter"}
    eas_frozen = {
        "payload_sha256": "eas-frozen",
        "cell_authorizations": {f"cell-{index}": {} for index in range(6)},
    }
    han_auxiliary = {
        "manifest_sha256": "han-manifest",
        "cells": {f"cell-{index}": {} for index in range(6)},
    }
    han_binding = {"schema": "han-binding", "mapping": {}}
    return eas_authorization, eas_frozen, han_auxiliary, han_binding


def _fake_integration_repo(tmp_path: Path):
    root = tmp_path / "repo"
    root.mkdir()
    units = _execution_units(root)
    for relative in production.INTEGRATION_IMPLEMENTATION_PATHS[:-1]:
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f"source: {relative}\n", encoding="utf-8")
    slim = root / "bin/slim"
    slim.parent.mkdir(parents=True, exist_ok=True)
    slim.write_bytes(b"slim-4.2.2")
    paths = {}
    for key, relative in {
        "eas_adapter": production.DEFAULT_EAS_ADAPTER_AUTHORIZATION_RELATIVE_PATH,
        "eas_frozen": production.DEFAULT_EAS_FROZEN_RELATIVE_PATH,
        "han_aux": production.DEFAULT_HAN_AUXILIARY_RELATIVE_PATH,
        "han_frozen": production.DEFAULT_HAN_FROZEN_RELATIVE_PATH,
        "han_manifest": production.DEFAULT_HAN_MANIFEST_RELATIVE_PATH,
    }.items():
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f"{key}\n", encoding="utf-8")
        paths[key] = path
    return root, units, slim, paths


def test_execution_unit_contract_is_exact(tmp_path):
    root = tmp_path / "repo"
    root.mkdir()
    units = _execution_units(root)
    counts = production._validate_execution_units(units)
    assert counts == {
        "total_units": 720,
        "selected_units": 120,
        "neutral_units": 600,
        "selected_cells": 12,
        "selected_replicates_per_cell": 10,
        "neutral_cells": 6,
        "neutral_replicates_per_cell": 100,
    }
    broken = units.iloc[:-1].copy()
    with pytest.raises(ValueError, match="720-unit"):
        production._validate_execution_units(broken)


def test_integration_manifest_is_immutable_and_revalidates_sources(
    tmp_path, monkeypatch
):
    root, _, slim, paths = _fake_integration_repo(tmp_path)
    eas_authorization, eas_frozen, han_auxiliary, han_binding = _fake_authorizations()
    monkeypatch.setattr(
        production,
        "_load_eas_authorization_bundle",
        lambda **_kwargs: (eas_authorization, eas_frozen),
    )
    monkeypatch.setattr(
        production,
        "_load_han_authorization_bundle",
        lambda **_kwargs: (han_auxiliary, han_binding),
    )
    manifest_path = root / production.DEFAULT_INTEGRATION_RELATIVE_PATH
    payload = production.build_integration_manifest(
        repo_root=root,
        campaign_dir=root / "focused_selection_EAS_sim",
        slim_path=slim,
        eas_adapter_authorization_path=paths["eas_adapter"],
        eas_frozen_path=paths["eas_frozen"],
        han_auxiliary_path=paths["han_aux"],
        han_frozen_path=paths["han_frozen"],
        han_manifest_path=paths["han_manifest"],
        integration_manifest_path=manifest_path,
    )
    production.write_integration_manifest(manifest_path, payload)
    production.write_integration_manifest(manifest_path, payload)
    bundle = production.load_selected_production_bundle(
        manifest_path, repo_root=root, slim_path=slim
    )
    assert bundle.manifest["payload_sha256"] == payload["payload_sha256"]
    assert (
        bundle.implementation_paths[-1] == production.DEFAULT_INTEGRATION_RELATIVE_PATH
    )
    changed_source = root / production.MODULE_SOURCE_PATH
    changed_source.write_text("changed\n", encoding="utf-8")
    with pytest.raises(ValueError, match="implementation"):
        production.load_selected_production_bundle(
            manifest_path, repo_root=root, slim_path=slim
        )


def test_integration_manifest_payload_tamper_fails_before_loader(tmp_path, monkeypatch):
    root, _, slim, paths = _fake_integration_repo(tmp_path)
    eas_authorization, eas_frozen, han_auxiliary, han_binding = _fake_authorizations()
    monkeypatch.setattr(
        production,
        "_load_eas_authorization_bundle",
        lambda **_kwargs: (eas_authorization, eas_frozen),
    )
    monkeypatch.setattr(
        production,
        "_load_han_authorization_bundle",
        lambda **_kwargs: (han_auxiliary, han_binding),
    )
    manifest_path = root / production.DEFAULT_INTEGRATION_RELATIVE_PATH
    payload = production.build_integration_manifest(
        repo_root=root,
        campaign_dir=root / "focused_selection_EAS_sim",
        slim_path=slim,
        eas_adapter_authorization_path=paths["eas_adapter"],
        eas_frozen_path=paths["eas_frozen"],
        han_auxiliary_path=paths["han_aux"],
        han_frozen_path=paths["han_frozen"],
        han_manifest_path=paths["han_manifest"],
        integration_manifest_path=manifest_path,
    )
    production.write_integration_manifest(manifest_path, payload)
    tampered = json.loads(manifest_path.read_text(encoding="utf-8"))
    tampered["status"] = "tampered"
    manifest_path.write_text(json.dumps(tampered), encoding="utf-8")
    with pytest.raises(ValueError, match="corrupt"):
        production.load_selected_production_bundle(
            manifest_path, repo_root=root, slim_path=slim
        )


def test_eas_adapter_authorization_checks_every_bound_sidecar(tmp_path, monkeypatch):
    root = tmp_path / "repo"
    root.mkdir()
    campaign_dir = root / "focused_selection_EAS_sim"
    output = campaign_dir / "calibration/eas_selected"
    output.mkdir(parents=True)
    units = campaign_dir / "execution_units.tsv"
    units.write_text("seed\n1\n", encoding="utf-8")
    slim = root / "bin/slim"
    slim.parent.mkdir(parents=True)
    slim.write_bytes(b"slim")
    adapter_source = root / production.eas_audit.MODULE_SOURCE_PATH
    adapter_source.parent.mkdir(parents=True, exist_ok=True)
    adapter_source.write_text("adapter\n", encoding="utf-8")
    frozen_path = output / "eas_selected_calibration_frozen.json"
    frozen_path.write_text("frozen\n", encoding="utf-8")
    frozen = {
        "payload_sha256": "frozen-payload",
        "slim": {"slim_scaling_factor": 5.0, "slim_burn_in": 0.1},
        "production_contract": {
            "selection_origin": "de_novo_at_prespecified_Q5_origin_age",
            "selection_episode": "continuous_from_origin_to_present",
            "terminal_conditions": (
                "lower_and_upper_inclusive_target_plus_minus_0.025"
            ),
            "candidate_pool_diploids": 500,
            "sample_diploids": 100,
        },
    }
    plan_manifest = {
        "source_sha256": {"legacy": "source"},
        "contract_sha256": "plan-contract",
    }
    sidecars = {}
    for phase in production.eas_calibration.PHASE_ORDER:
        path = output / f"{phase}_adapter_audit.json"
        path.write_text(f"{phase}\n", encoding="utf-8")
        sidecars[phase] = {"payload_sha256": f"payload-{phase}"}
    monkeypatch.setattr(
        production.eas_calibration,
        "load_frozen_eas_selected_calibration",
        lambda *_args, **_kwargs: frozen,
    )
    monkeypatch.setattr(
        production.eas_audit,
        "_load_bound_inputs",
        lambda **_kwargs: (output, campaign_dir / "work", {}, plan_manifest),
    )
    monkeypatch.setattr(
        production.eas_audit,
        "_load_phase_sidecar",
        lambda phase, **_kwargs: sidecars[phase],
    )
    authorization = {
        "schema": production.eas_audit.AUTHORIZATION_SCHEMA_VERSION,
        "status": "authorized",
        "adapter": {
            "path": production.eas_audit.MODULE_SOURCE_PATH,
            "sha256": production.sha256_file(adapter_source),
        },
        "legacy_source_sha256": plan_manifest["source_sha256"],
        "plan_manifest_contract_sha256": plan_manifest["contract_sha256"],
        "phase_audits": {
            phase: {
                **production._file_record(output / f"{phase}_adapter_audit.json", root),
                "payload_sha256": sidecars[phase]["payload_sha256"],
            }
            for phase in production.eas_calibration.PHASE_ORDER
        },
        "frozen_authorization": {
            **production._file_record(frozen_path, root),
            "payload_sha256": frozen["payload_sha256"],
        },
    }
    authorization["payload_sha256"] = production._canonical_sha256(authorization)
    authorization_path = output / "eas_selected_calibration_frozen.adapter.json"
    authorization_path.write_text(json.dumps(authorization), encoding="utf-8")
    production._load_eas_authorization_bundle(
        repo_root=root,
        campaign_dir=campaign_dir,
        slim_path=slim,
        execution_units_path=units,
        adapter_authorization_path=authorization_path,
        frozen_path=frozen_path,
    )
    (output / "screen20_adapter_audit.json").write_text("tampered\n", encoding="utf-8")
    with pytest.raises(ValueError, match="screen20 audit"):
        production._load_eas_authorization_bundle(
            repo_root=root,
            campaign_dir=campaign_dir,
            slim_path=slim,
            execution_units_path=units,
            adapter_authorization_path=authorization_path,
            frozen_path=frozen_path,
        )


def test_han_auxiliary_must_match_loader_binding(tmp_path, monkeypatch):
    root = tmp_path / "repo"
    root.mkdir()
    auxiliary_path = root / "aux.json"
    frozen_path = root / "frozen.json"
    manifest_path = root / "manifest.json"
    units_path = root / "units.tsv"
    for path in (auxiliary_path, frozen_path, manifest_path, units_path):
        path.write_text(path.name, encoding="utf-8")
    binding = {"schema": "binding", "mapping": {}}
    auxiliary = {
        "loader_compatible_frozen_path": frozen_path.relative_to(root).as_posix(),
        "loader_compatible_frozen_sha256": production.sha256_file(frozen_path),
        "production_loader_binding_sha256": production._canonical_sha256(binding),
    }
    monkeypatch.setattr(
        production.han_fixed_v2,
        "load_frozen_fixed_v2",
        lambda *_args, **_kwargs: auxiliary,
    )
    monkeypatch.setattr(
        production.simulation,
        "load_frozen_han_selection_calibration",
        lambda *_args, **_kwargs: binding,
    )
    loaded_auxiliary, loaded_binding = production._load_han_authorization_bundle(
        repo_root=root,
        execution_units_path=units_path,
        auxiliary_path=auxiliary_path,
        frozen_path=frozen_path,
        manifest_path=manifest_path,
        slim_scaling_factor=5.0,
    )
    assert loaded_auxiliary is auxiliary
    assert loaded_binding is binding
    auxiliary["production_loader_binding_sha256"] = "wrong"
    with pytest.raises(ValueError, match="loader binding"):
        production._load_han_authorization_bundle(
            repo_root=root,
            execution_units_path=units_path,
            auxiliary_path=auxiliary_path,
            frozen_path=frozen_path,
            manifest_path=manifest_path,
            slim_scaling_factor=5.0,
        )


def test_process_local_patch_uses_authorized_eas_cell_and_restores(
    tmp_path, monkeypatch
):
    root = tmp_path / "repo"
    root.mkdir()
    bundle = _bundle(root)
    original_objects = simulation._selected_objects
    original_paths = simulation.IMPLEMENTATION_PATHS
    origin = {"age_generations": 999.0}
    sweep = SimpleNamespace(
        extended_events=(), provenance_record=lambda: {"events": "calibrated"}
    )
    serialized = [
        {
            "event_type": "ConditionOnAlleleFrequency",
            "operator": operator,
            "start_time": {"generations_ago": 0.0},
            "end_time": {"generations_ago": 0.0},
        }
        for operator in (">=", "<=")
    ]
    monkeypatch.setattr(
        production.eas_calibration,
        "_load_eas_cell",
        lambda *_args, **_kwargs: (
            "model",
            sweep,
            "demographic",
            origin,
            {"event_contract_sha256": "events"},
        ),
    )
    monkeypatch.setattr(
        production.eas_calibration,
        "serialize_extended_events",
        lambda _events: serialized,
    )
    with production.patched_selected_executor(bundle):
        model, loaded_sweep, samples, provenance = simulation._selected_objects(
            _selected_unit(),
            root,
            5.0,
            500,
            han_selection_calibration=bundle.han_binding,
        )
        assert model == "model"
        assert loaded_sweep is sweep
        assert samples == {"EAS": 500}
        assert (
            provenance["selected_production_authorization"][
                "terminal_population_af_conditioning_in_slim"
            ]
            is True
        )
        assert simulation.IMPLEMENTATION_PATHS == bundle.implementation_paths
    assert simulation._selected_objects is original_objects
    assert simulation.IMPLEMENTATION_PATHS == original_paths


def test_process_local_patch_requires_exact_han_binding(tmp_path, monkeypatch):
    root = tmp_path / "repo"
    root.mkdir()
    bundle = _bundle(root)

    def fake_original(*_args, **_kwargs):
        return "model", "sweep", {"Han": 500}, {"selection_endpoint": "fixed"}

    monkeypatch.setattr(simulation, "_selected_objects", fake_original)
    with production.patched_selected_executor(bundle):
        _, _, _, provenance = simulation._selected_objects(
            _selected_unit("ancient_eurasia_han_introgression"),
            root,
            5.0,
            500,
            han_selection_calibration=bundle.han_binding,
        )
        assert (
            provenance["selected_production_authorization"][
                "post_pulse_recipient_nonloss_conditioning"
            ]
            is True
        )
        with pytest.raises(ValueError, match="fixed-v2"):
            simulation._selected_objects(
                _selected_unit("ancient_eurasia_han_introgression"),
                root,
                5.0,
                500,
                han_selection_calibration={"wrong": True},
            )


def test_quarantine_plan_inventories_exactly_24_selected_dirs_without_moving(
    tmp_path, monkeypatch
):
    root = tmp_path / "repo"
    root.mkdir()
    units = _execution_units(root)
    campaign_dir = root / "focused_selection_EAS_sim"
    work = campaign_dir / "work"
    selected = units[units["simulation_class"].eq("selected")].head(24)
    for index, unit_id in enumerate(selected["unit_id"]):
        directory = work / unit_id
        directory.mkdir(parents=True)
        (directory / "attempts.tsv").write_text(f"attempt\n{index}\n", encoding="utf-8")
    neutral_id = units[units["simulation_class"].eq("neutral")].iloc[0]["unit_id"]
    neutral = work / str(neutral_id)
    neutral.mkdir(parents=True)
    (neutral / "simulation_complete.json").write_text("neutral\n", encoding="utf-8")
    monkeypatch.setattr(production.simulation, "unit_lock_is_held", lambda _path: False)
    payload = production.build_quarantine_plan(
        repo_root=root,
        campaign_dir=campaign_dir,
        expected_stale_count=24,
    )
    assert payload["observed_stale_selected_directories"] == 24
    assert payload["moves_performed"] == payload["deletions_performed"] == 0
    assert all(row["simulation_class"] == "selected" for row in payload["rows"])
    assert all(Path(root / row["source"]).is_dir() for row in payload["rows"])
    assert all(not Path(root / row["destination"]).exists() for row in payload["rows"])
    assert neutral.is_dir()
    with pytest.raises(ValueError, match="observed 24, expected 23"):
        production.build_quarantine_plan(
            repo_root=root,
            campaign_dir=campaign_dir,
            expected_stale_count=23,
        )


def test_quarantine_plan_excludes_already_integrated_selected_directory(
    tmp_path, monkeypatch
):
    root = tmp_path / "repo"
    root.mkdir()
    units = _execution_units(root)
    campaign_dir = root / "focused_selection_EAS_sim"
    selected = units[units["simulation_class"].eq("selected")].head(2)
    for unit_id in selected["unit_id"]:
        directory = campaign_dir / "work" / unit_id
        directory.mkdir(parents=True)
        (directory / "data").write_text("data\n", encoding="utf-8")
    module_source = root / production.MODULE_SOURCE_PATH
    module_source.parent.mkdir(parents=True, exist_ok=True)
    module_source.write_text("adapter\n", encoding="utf-8")
    integration_manifest = root / production.DEFAULT_INTEGRATION_RELATIVE_PATH
    integration_manifest.parent.mkdir(parents=True, exist_ok=True)
    integration_manifest.write_text("integration\n", encoding="utf-8")
    integrated = campaign_dir / "work" / selected.iloc[0]["unit_id"]
    (integrated / "simulation_contract.json").write_text(
        json.dumps(
            {
                "implementation": {
                    "sources": {
                        production.MODULE_SOURCE_PATH: production.sha256_file(
                            module_source
                        ),
                        production.DEFAULT_INTEGRATION_RELATIVE_PATH: (
                            production.sha256_file(integration_manifest)
                        ),
                    }
                }
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr(production.simulation, "unit_lock_is_held", lambda _path: False)
    payload = production.build_quarantine_plan(
        repo_root=root,
        campaign_dir=campaign_dir,
        expected_stale_count=1,
    )
    assert payload["rows"][0]["unit_id"] == selected.iloc[1]["unit_id"]


def test_selected_request_refuses_neutral_unit(tmp_path):
    root = tmp_path / "repo"
    root.mkdir()
    units = _execution_units(root)
    neutral = units[units["simulation_class"].eq("neutral")].head(1)
    with pytest.raises(ValueError, match="refuses neutral"):
        production._validate_requested_units(neutral, units)


def test_selected_request_rejects_duplicate_unit_and_retry_budget_drift(
    tmp_path, monkeypatch
):
    root = tmp_path / "repo"
    root.mkdir()
    units = _execution_units(root)
    selected = units[units["simulation_class"].eq("selected")].head(1)
    duplicated = pd.concat([selected, selected], ignore_index=True)
    with pytest.raises(ValueError, match="duplicate unit"):
        production._validate_requested_units(duplicated, units)

    bundle = _bundle(root)
    bundle = production.SelectedProductionBundle(
        **{
            **bundle.__dict__,
            "manifest": {
                "payload_sha256": "integration",
                "runtime": {"eas_max_draws": 68, "han_max_draws": 100},
            },
        }
    )
    monkeypatch.setattr(
        production, "load_selected_production_bundle", lambda *_args, **_kwargs: bundle
    )
    with pytest.raises(ValueError, match="retry budgets"):
        production.simulate_selected_units(
            selected,
            repo_root=root,
            integration_manifest_path=bundle.integration_manifest_path,
            slim_path=bundle.slim_path,
            eas_max_draws=69,
        )
