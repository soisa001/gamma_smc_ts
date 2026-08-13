from __future__ import annotations

import json
import multiprocessing
import os
import subprocess
import time
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
from gamma_smc_aou import focused_han_direct_v3 as han_v3
from gamma_smc_aou import focused_han_direct_v3_recovery as han_recovery
from gamma_smc_aou import focused_selected_production as production
from gamma_smc_aou import focused_selection_analysis as focused_analysis
from gamma_smc_aou import focused_selection_campaign as campaign
from gamma_smc_aou import focused_selection_report as focused_report
from gamma_smc_aou import focused_selection_simulation as simulation


class _HangingEngineWithGrandchild:
    def simulate(self, *_args, **kwargs):
        child = subprocess.Popen(["/bin/sh", "-c", "trap '' TERM; sleep 60"])
        Path(str(kwargs["logfile"]) + ".grandchild.pid").write_text(
            str(child.pid), encoding="utf-8"
        )
        while True:
            time.sleep(1)


class _FastTreeEngine:
    def simulate(self, *_args, **_kwargs):
        import tskit

        return tskit.TableCollection(1).tree_sequence()


def _write_accepted_partial_worker(unit_dir_raw: str) -> None:
    unit_dir = Path(unit_dir_raw)
    schedule = pd.DataFrame(
        [{"attempt_zero_based": 0, "simulation_seed": 12345, "panel_seed": 54321}]
    )
    proxy = han_v3._SupervisedSlimEngine(
        _FastTreeEngine(), unit_dir=unit_dir, schedule=schedule, term_grace_seconds=0.1
    )
    with han_v3._supervised_draw_timeout(5.0):
        proxy.simulate(
            seed=12345,
            logfile=str(unit_dir / "slim_draw_000.csv"),
        )
    attempts = pd.DataFrame(
        [
            {
                "attempt_zero_based": 0,
                "seed": 12345,
                "status": "complete",
                "timeout_seconds": 5.0,
                "elapsed_seconds": 0.1,
                "accepted": True,
            }
        ]
    )
    simulation._atomic_frame(unit_dir / "attempts.tsv", attempts)
    proxy.mark_base_attempt_ledger_committed(attempts)
    for name in ("simulation.trees", "sample_manifest.tsv", "truth_profiles.tsv.gz"):
        (unit_dir / name).write_bytes(b"partial")
    pairs = unit_dir / "pairs"
    pairs.mkdir()
    (pairs / "overall.pairs.tsv").write_bytes(b"partial")
    os._exit(0)


def _execution_units(root: Path) -> pd.DataFrame:
    campaign_dir = root / production.DEFAULT_CAMPAIGN_RELATIVE_PATH
    units = campaign.build_execution_units(
        campaign.FocusedCampaignPlan(repo_root=root, campaign_dir=campaign_dir)
    )
    campaign_dir.mkdir(parents=True, exist_ok=True)
    units.to_csv(
        campaign_dir / "execution_units.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    return units


def _selected_unit(demography_id: str = "eas_phlash_median") -> dict[str, object]:
    return {
        "unit_id": f"test-{demography_id}",
        "demography_id": demography_id,
        "demography_kind": (
            "no_introgression"
            if demography_id == "eas_phlash_median"
            else "introgression"
        ),
        "population": "EAS" if demography_id == "eas_phlash_median" else "Han",
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
    manifest = root / production.DEFAULT_INTEGRATION_RELATIVE_PATH
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_text("{}\n", encoding="utf-8")
    slim = root / "bin/slim"
    slim.parent.mkdir(parents=True, exist_ok=True)
    slim.write_bytes(b"slim")
    origin = {"age_generations": 999.0}
    return production.SelectedProductionBundle(
        repo_root=root,
        campaign_dir=root / production.DEFAULT_CAMPAIGN_RELATIVE_PATH,
        integration_manifest_path=manifest,
        slim_path=slim,
        execution_units_path=(
            root / production.DEFAULT_CAMPAIGN_RELATIVE_PATH / "execution_units.tsv"
        ),
        manifest={"payload_sha256": "integration"},
        eas_authorization={"payload_sha256": "adapter"},
        eas_frozen={
            "payload_sha256": "eas-frozen",
            "cell_authorizations": {
                "s=0.010000|af=0.100000": {
                    "origin_age_generations": 999.0,
                    "origin_contract_sha256": production._canonical_sha256(origin),
                    "event_contract_sha256": "events",
                }
            },
        },
        implementation_paths=production.INTEGRATION_IMPLEMENTATION_PATHS,
    )


def _write_valid_eas_restart_ledgers(
    bundle: production.SelectedProductionBundle,
    unit: dict[str, object],
    monkeypatch,
) -> tuple[Path, dict[str, object]]:
    supervisor_source = bundle.repo_root / han_v3.MODULE_SOURCE_PATH
    supervisor_source.parent.mkdir(parents=True, exist_ok=True)
    supervisor_source.write_bytes(b"hardened supervisor source\n")
    monkeypatch.setattr(
        simulation.stdpopsim,
        "get_engine",
        lambda name, *_args, **_kwargs: _FastTreeEngine(),
    )
    unit_dir = bundle.campaign_dir / "work" / str(unit["unit_id"])
    schedule = production._eas_supervisor_schedule(unit).set_index("attempt_zero_based")
    with production._patched_eas_process_supervisor(bundle, unit) as supervisor:
        engine = simulation.stdpopsim.get_engine("slim")
        for attempt in (0, 1):
            with simulation._wall_clock_timeout(
                production.DEFAULT_DRAW_TIMEOUT_SECONDS
            ):
                engine.simulate(
                    seed=int(schedule.loc[attempt, "simulation_seed"]),
                    logfile=str(unit_dir / f"slim_draw_{attempt:03d}.csv"),
                )
        supervisor.records[0]["supervisor_status"] = "interrupted_by_worker_death"
        supervisor.records[0]["lifecycle_phase"] = (
            "accepted_attempt_interrupted_before_unit_completion"
        )
        attempts = pd.DataFrame(
            [
                {
                    "attempt_zero_based": 0,
                    "seed": int(schedule.loc[0, "simulation_seed"]),
                    "status": "timed_out",
                    "timeout_seconds": production.DEFAULT_DRAW_TIMEOUT_SECONDS,
                    "elapsed_seconds": float(
                        supervisor.records[0]["supervisor_elapsed_seconds"]
                    )
                    + 0.01,
                    "accepted": False,
                    "rejection_reason": (
                        "worker_interrupted_after_acceptance_before_unit_completion"
                    ),
                },
                {
                    "attempt_zero_based": 1,
                    "seed": int(schedule.loc[1, "simulation_seed"]),
                    "status": "complete",
                    "timeout_seconds": production.DEFAULT_DRAW_TIMEOUT_SECONDS,
                    "elapsed_seconds": float(
                        supervisor.records[1]["supervisor_elapsed_seconds"]
                    )
                    + 0.01,
                    "accepted": True,
                    "rejection_reason": "",
                },
            ]
        )
        simulation._atomic_frame(unit_dir / "attempts.tsv", attempts)
        first = supervisor.records[0]
        reconciliation = pd.DataFrame(
            [
                {
                    "schema": han_v3.SUPERVISOR_RECONCILIATION_SCHEMA,
                    "reason": ("accepted_attempt_interrupted_before_unit_completion"),
                    "execution_id": first["execution_id"],
                    "attempt_zero_based": 0,
                    "simulation_seed": int(schedule.loc[0, "simulation_seed"]),
                    "active_record_path": ".han_v3_active_000.json",
                    "active_record_sha256": "a" * 64,
                    "prior_supervisor_pid": int(first["supervisor_pid"]),
                    "prior_child_pid": int(first["child_pid"]),
                    "prior_child_pgid": int(first["child_pgid"]),
                    "process_group_gone": True,
                    "temp_directory_existed": False,
                    "temp_directory_removed": True,
                    "charged_elapsed_seconds": 0.01,
                    "transferred_to_attempt_elapsed_seconds": 0.0,
                    "atomic_temp_files_removed": "[]",
                    "atomic_temp_file_count_removed": 0,
                    "reconciled_unix_ns": time.time_ns(),
                }
            ]
        )
        simulation._atomic_frame(
            unit_dir / "supervision_reconciliation.tsv", reconciliation
        )
        (unit_dir / "simulation_complete.json").write_text("{}\n", encoding="utf-8")
    return unit_dir, {"attempts_completed": 2}


def _write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _signed(payload: dict[str, object]) -> dict[str, object]:
    result = dict(payload)
    result["payload_sha256"] = production._canonical_sha256(result)
    return result


def _compact_phase(
    root: Path,
    phase: str,
    *,
    adapter_sha256: str,
    source_sha256: dict[str, str],
    plan_contract_sha256: str,
    slim: Path,
) -> Path:
    output = (
        root / production.DEFAULT_CAMPAIGN_RELATIVE_PATH / "calibration/eas_selected"
    )
    output.mkdir(parents=True, exist_ok=True)
    resource = root / production.eas_calibration.EAS_RESOURCE_PATH
    resource.parent.mkdir(parents=True, exist_ok=True)
    resource.write_bytes(b"EAS")
    plan_manifest = output / "eas_selected_calibration_plan.json"
    if not plan_manifest.exists():
        plan_manifest.write_text("{}\n", encoding="utf-8")
    plan = output / f"{phase}_plan.tsv"
    ledger = output / f"{phase}_ledger.tsv"
    summary = output / f"{phase}_summary.tsv"
    for path in (plan, ledger, summary):
        path.write_text(f"phase\n{phase}\n", encoding="utf-8")

    rows = []
    expected_rows = 6 * production.eas_calibration.PHASE_DRAWS[phase]
    for index in range(expected_rows):
        rows.append(
            {
                "calibration_id": f"{phase}-{index:04d}",
                "completion_path": (
                    f"{production.DEFAULT_CAMPAIGN_RELATIVE_PATH}/work/"
                    f"eas_selected_calibration/{phase}/{phase}-{index:04d}/"
                    "completion.json"
                ),
                "completion_sha256": f"{index + 1:064x}",
                "trajectory_sha256": f"{index + 10_001:064x}",
                "panel_tree_sha256": "" if index % 2 else f"{index + 20_001:064x}",
                "panel_manifest_sha256": (
                    "" if index % 2 else f"{index + 30_001:064x}"
                ),
            }
        )
    checksums = output / f"{phase}_completion_checksums.tsv"
    checksum_frame = pd.DataFrame(rows)
    checksum_frame.to_csv(checksums, sep="\t", index=False, lineterminator="\n")
    mapping = production._compact_completion_mapping(
        checksum_frame, phase=phase, repo_root=root
    )
    payload = _signed(
        {
            "schema": production.eas_audit.SCHEMA_VERSION,
            "status": "validated",
            "phase": phase,
            "adapter": {
                "path": production.eas_audit.MODULE_SOURCE_PATH,
                "sha256": adapter_sha256,
            },
            "legacy_source_sha256": source_sha256,
            "plan_manifest": {
                **production._file_record(plan_manifest, root),
                "contract_sha256": plan_contract_sha256,
            },
            "plan": production._file_record(plan, root),
            "outputs": {
                "rebuilt_ledger": production._file_record(ledger, root),
                "summary": production._file_record(summary, root),
                "completion_checksums": production._file_record(checksums, root),
            },
            "normalization": {
                "input_ledger": production._file_record(ledger, root),
                "policy": "tracked compact evidence",
            },
            "slim": production._file_record(slim, root),
            "eas_resource": production._file_record(resource, root),
            "counts": {
                "planned": expected_rows,
                "validated_completions": expected_rows,
                "completion_checksum_rows": expected_rows,
                "cells": 6,
                "passing_cells": 6,
            },
            "completion_mapping": mapping,
        }
    )
    sidecar = output / f"{phase}_adapter_audit.json"
    _write_json(sidecar, payload)
    return sidecar


def _fake_manifest_repo(tmp_path: Path) -> tuple[Path, Path, dict[str, Path]]:
    root = tmp_path / "repo"
    root.mkdir()
    _execution_units(root)
    for relative in production.INTEGRATION_IMPLEMENTATION_PATHS[:-1]:
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f"source: {relative}\n", encoding="utf-8")
    slim = root / "bin/slim"
    slim.parent.mkdir(parents=True, exist_ok=True)
    slim.write_bytes(b"slim-4.2.2")
    paths: dict[str, Path] = {}
    for key, relative in {
        "adapter": production.DEFAULT_EAS_ADAPTER_AUTHORIZATION_RELATIVE_PATH,
        "frozen": production.DEFAULT_EAS_FROZEN_RELATIVE_PATH,
    }.items():
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f"{key}\n", encoding="utf-8")
        paths[key] = path
    return root, slim, paths


def _fake_eas_authorizations() -> tuple[dict[str, object], dict[str, object]]:
    return (
        {"payload_sha256": "eas-adapter"},
        {
            "payload_sha256": "eas-frozen",
            "cell_authorizations": {f"cell-{index}": {} for index in range(6)},
        },
    )


def _han_audit(
    *, stage: str, units: int, replicates: int, manifest: Path, seed_plan: Path
) -> dict[str, object]:
    return _signed(
        {
            "schema": han_v3.STAGE_AUDIT_SCHEMA,
            "status": "passed",
            "stage": stage,
            "expected_units": units,
            "accepted_units": units,
            "per_cell_replicates": list(range(1, replicates + 1)),
            "manifest_sha256": production.sha256_file(manifest),
            "seed_plan_sha256": production.sha256_file(seed_plan),
            "failed_v2_inferential_use": False,
            "gamma_smc_statistics_used": False,
            "all_current_completion_budgets_match_authorization": True,
            "all_current_completion_counts_match_attempt_ledgers": True,
            "all_current_completion_sources_match_pinned_runtime": True,
            "all_current_completions_cache_validated_read_only": True,
            "all_rejections_persisted": True,
            "all_unit_authorizations_match_manifest_and_seed_plan": True,
            "all_units_strict_type1_identity": True,
            "all_units_use_exact_final_contract": True,
            "pilot_units_count_toward_final": True,
            "seed_plan_audit": {},
            "rows": [
                {"unit_id": f"unit-{index:03d}", "attempts_completed": 1}
                for index in range(units)
            ],
        }
    )


def _han_authorization(
    *,
    manifest: Path,
    seed_plan: Path,
    manifest_payload: dict[str, object],
    final_audit: dict[str, object],
) -> dict[str, object]:
    return _signed(
        {
            "schema": han_v3.AUTHORIZATION_SCHEMA,
            "status": "authorized",
            "manifest_sha256": production.sha256_file(manifest),
            "seed_plan_sha256": production.sha256_file(seed_plan),
            "all_units_audit_payload_sha256": final_audit["payload_sha256"],
            "accepted_units": 60,
            "cells": 6,
            "replicates_per_cell": 10,
            "pilot_units_included": 6,
            "validated_completion_count": 60,
            "validated_attempt_count": 60,
            "attempts_per_unit_budget": han_v3.ATTEMPTS_PER_UNIT,
            "draw_timeout_seconds": han_v3.DRAW_TIMEOUT_SECONDS,
            "cumulative_timeout_seconds": han_v3.CUMULATIVE_TIMEOUT_SECONDS,
            "all_completion_counts_budgets_and_authorizations_validated": True,
            "all_process_groups_gone_and_children_reaped": True,
            "pinned_runtime": manifest_payload["runtime"],
            "direct_operational_estimand": True,
            "conditional_inference": True,
            "gamma_smc_statistics_used": False,
            "failed_v2_inferential_use": False,
            "scientific_caveats": manifest_payload["scientific_caveats"],
        }
    )


def test_byte_bound_adapter_sources_are_preserved_across_platforms():
    root = Path(__file__).resolve().parents[1]
    for relative in (production.MODULE_SOURCE_PATH, production.WRAPPER_SOURCE_PATH):
        source = root / relative
        rules = {
            line.split()[0]: line.split()[1:]
            for line in (source.parent / ".gitattributes")
            .read_text(encoding="utf-8")
            .splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        }
        assert "-text" in rules[source.name]


def test_execution_contract_and_eas_grid_are_exact(tmp_path: Path):
    units = _execution_units(tmp_path)
    counts = production._validate_execution_units(units)
    eas = production._eas_production_units(units)
    assert counts["total_units"] == 720
    assert counts["selected_units"] == 120
    assert counts["neutral_units"] == 600
    assert len(eas) == 60
    assert set(eas["demography_id"]) == {"eas_phlash_median"}
    assert set(
        eas.groupby(["selection_coefficient", "target_allele_frequency"]).size()
    ) == {10}


def test_compact_sidecar_validates_with_ignored_raw_work_tree_absent(tmp_path: Path):
    root = tmp_path / "repo"
    root.mkdir()
    slim = root / "bin/slim"
    slim.parent.mkdir(parents=True)
    slim.write_bytes(b"slim")
    source = {"immutable.py": "a" * 64}
    sidecar = _compact_phase(
        root,
        "screen20",
        adapter_sha256="b" * 64,
        source_sha256=source,
        plan_contract_sha256="c" * 64,
        slim=slim,
    )
    assert not (root / production.DEFAULT_CAMPAIGN_RELATIVE_PATH / "work").exists()
    loaded = production._load_compact_phase_sidecar(
        "screen20",
        repo_root=root,
        output=sidecar.parent,
        adapter_sha256="b" * 64,
        expected_source_sha256=source,
        expected_plan_contract_sha256="c" * 64,
        slim_path=slim,
    )
    assert loaded["completion_mapping"]["rows"] == 120
    assert not (root / production.DEFAULT_CAMPAIGN_RELATIVE_PATH / "work").exists()


def test_compact_sidecar_rejects_checksum_table_and_mapping_corruption(tmp_path: Path):
    root = tmp_path / "repo"
    root.mkdir()
    slim = root / "slim"
    slim.write_bytes(b"slim")
    source = {"immutable.py": "a" * 64}
    sidecar = _compact_phase(
        root,
        "sensitivity10mb20",
        adapter_sha256="b" * 64,
        source_sha256=source,
        plan_contract_sha256="c" * 64,
        slim=slim,
    )
    checksums = sidecar.parent / "sensitivity10mb20_completion_checksums.tsv"
    checksums.write_text("corrupt\n", encoding="utf-8")
    with pytest.raises(ValueError, match="completion_checksums file binding"):
        production._load_compact_phase_sidecar(
            "sensitivity10mb20",
            repo_root=root,
            output=sidecar.parent,
            adapter_sha256="b" * 64,
            expected_source_sha256=source,
            expected_plan_contract_sha256="c" * 64,
            slim_path=slim,
        )


def test_eas_authorization_uses_only_compact_tracked_chain(tmp_path: Path, monkeypatch):
    root = tmp_path / "repo"
    root.mkdir()
    campaign_dir = root / production.DEFAULT_CAMPAIGN_RELATIVE_PATH
    campaign_dir.mkdir()
    units = campaign_dir / "execution_units.tsv"
    units.write_text("seed\n1\n", encoding="utf-8")
    slim = root / "bin/slim"
    slim.parent.mkdir(parents=True)
    slim.write_bytes(b"slim")
    adapter_source = root / production.eas_audit.MODULE_SOURCE_PATH
    adapter_source.parent.mkdir(parents=True, exist_ok=True)
    adapter_source.write_text("adapter\n", encoding="utf-8")
    source = {"immutable.py": "a" * 64}
    contract = "c" * 64
    sidecars = {
        phase: _compact_phase(
            root,
            phase,
            adapter_sha256=production.sha256_file(adapter_source),
            source_sha256=source,
            plan_contract_sha256=contract,
            slim=slim,
        )
        for phase in production.eas_calibration.PHASE_ORDER
    }
    output = next(iter(sidecars.values())).parent
    frozen_path = output / "eas_selected_calibration_frozen.json"
    frozen_path.write_text("frozen\n", encoding="utf-8")
    frozen = {
        "payload_sha256": "frozen-payload",
        "source_sha256": source,
        "plan_manifest_contract_sha256": contract,
        "slim": {"slim_scaling_factor": 5.0, "slim_burn_in": 0.1},
        "production_contract": {
            "selection_origin": "de_novo_at_prespecified_Q5_origin_age",
            "selection_episode": "continuous_from_origin_to_present",
            "terminal_conditions": "lower_and_upper_inclusive_target_plus_minus_0.025",
            "candidate_pool_diploids": 500,
            "sample_diploids": 100,
        },
    }
    monkeypatch.setattr(
        production.eas_calibration,
        "load_frozen_eas_selected_calibration",
        lambda *_args, **_kwargs: frozen,
    )
    authorization = _signed(
        {
            "schema": production.eas_audit.AUTHORIZATION_SCHEMA_VERSION,
            "status": "authorized",
            "adapter": {
                "path": production.eas_audit.MODULE_SOURCE_PATH,
                "sha256": production.sha256_file(adapter_source),
            },
            "legacy_source_sha256": source,
            "plan_manifest_contract_sha256": contract,
            "phase_audits": {
                phase: {
                    **production._file_record(path, root),
                    "payload_sha256": json.loads(path.read_text())["payload_sha256"],
                }
                for phase, path in sidecars.items()
            },
            "frozen_authorization": {
                **production._file_record(frozen_path, root),
                "payload_sha256": frozen["payload_sha256"],
            },
        }
    )
    authorization_path = output / "eas_selected_calibration_frozen.adapter.json"
    _write_json(authorization_path, authorization)
    loaded, loaded_frozen = production._load_eas_authorization_bundle(
        repo_root=root,
        campaign_dir=campaign_dir,
        slim_path=slim,
        execution_units_path=units,
        adapter_authorization_path=authorization_path,
        frozen_path=frozen_path,
    )
    assert loaded["payload_sha256"] == authorization["payload_sha256"]
    assert loaded_frozen is frozen
    assert not (campaign_dir / "work").exists()


def test_eas_manifest_is_independent_of_han_and_revalidates_sources(
    tmp_path: Path, monkeypatch
):
    root, slim, paths = _fake_manifest_repo(tmp_path)
    authorizations = _fake_eas_authorizations()
    monkeypatch.setattr(
        production, "_load_eas_authorization_bundle", lambda **_kwargs: authorizations
    )
    manifest_path = root / production.DEFAULT_INTEGRATION_RELATIVE_PATH
    payload = production.build_integration_manifest(
        repo_root=root,
        campaign_dir=root / production.DEFAULT_CAMPAIGN_RELATIVE_PATH,
        slim_path=slim,
        eas_adapter_authorization_path=paths["adapter"],
        eas_frozen_path=paths["frozen"],
        integration_manifest_path=manifest_path,
    )
    assert set(payload["authorization_artifacts"]) == {"eas_adapter", "eas_frozen"}
    assert payload["scientific_contract"]["han"].startswith("separate")
    production.write_integration_manifest(manifest_path, payload)
    production.write_integration_manifest(manifest_path, payload)
    bundle = production.load_selected_production_bundle(
        manifest_path, repo_root=root, slim_path=slim
    )
    assert bundle.manifest["status"] == "eas_ready"
    (root / production.MODULE_SOURCE_PATH).write_text("changed\n", encoding="utf-8")
    with pytest.raises(ValueError, match="implementation"):
        production.load_selected_production_bundle(
            manifest_path, repo_root=root, slim_path=slim
        )


def test_process_patch_authorizes_eas_and_refuses_han(tmp_path: Path, monkeypatch):
    root = tmp_path / "repo"
    root.mkdir()
    bundle = _bundle(root)
    original_objects = simulation._selected_objects
    original_paths = simulation.IMPLEMENTATION_PATHS
    origin = {"age_generations": 999.0}
    sweep = SimpleNamespace(
        extended_events=(), provenance_record=lambda: {"events": "calibrated"}
    )
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
        lambda _events: [
            {
                "event_type": "ConditionOnAlleleFrequency",
                "operator": operator,
                "start_time": {"generations_ago": 0.0},
                "end_time": {"generations_ago": 0.0},
            }
            for operator in (">=", "<=")
        ],
    )
    with production.patched_selected_executor(bundle):
        model, _, samples, provenance = simulation._selected_objects(
            _selected_unit(), root, 5.0, 500
        )
        assert model == "model"
        assert samples == {"EAS": 500}
        assert provenance["selected_production_authorization"][
            "terminal_population_af_conditioning_in_slim"
        ]
        with pytest.raises(ValueError, match="focused_han_direct_v3"):
            simulation._selected_objects(
                _selected_unit("ancient_eurasia_han_introgression"), root, 5.0, 500
            )
    assert simulation._selected_objects is original_objects
    assert simulation.IMPLEMENTATION_PATHS == original_paths


def test_requested_units_require_exact_60_eas_set(tmp_path: Path):
    units = _execution_units(tmp_path)
    eas = production._eas_production_units(units)
    assert len(production._validate_requested_units(eas, units)) == 60
    with pytest.raises(ValueError, match="exactly 60"):
        production._validate_requested_units(eas.iloc[:-1], units)
    mixed = pd.concat(
        [
            eas.iloc[:-1],
            units[units["demography_id"].eq("ancient_eurasia_han_introgression")].iloc[
                :1
            ],
        ]
    )
    with pytest.raises(ValueError, match="neutral or Han"):
        production._validate_requested_units(mixed, units)


def test_retry_and_timeout_bounds_are_not_user_tunable(tmp_path: Path):
    units = _execution_units(tmp_path)
    eas = production._eas_production_units(units)
    with pytest.raises(ValueError, match="frozen at 68"):
        production.simulate_selected_units(
            eas,
            repo_root=tmp_path,
            integration_manifest_path=tmp_path / "manifest.json",
            slim_path=tmp_path / "slim",
            eas_max_draws=69,
        )
    with pytest.raises(ValueError, match="frozen at 30"):
        production.simulate_selected_units(
            eas,
            repo_root=tmp_path,
            integration_manifest_path=tmp_path / "manifest.json",
            slim_path=tmp_path / "slim",
            draw_timeout_seconds=1799,
        )


def test_quarantine_plan_is_read_only_and_excludes_han(tmp_path: Path, monkeypatch):
    root = tmp_path / "repo"
    root.mkdir()
    units = _execution_units(root)
    campaign_dir = root / production.DEFAULT_CAMPAIGN_RELATIVE_PATH
    eas = production._eas_production_units(units).head(2)
    han = units[
        units["simulation_class"].eq("selected")
        & units["demography_id"].eq("ancient_eurasia_han_introgression")
    ].head(1)
    for row in pd.concat([eas, han]).to_dict(orient="records"):
        directory = campaign_dir / "work" / str(row["unit_id"])
        directory.mkdir(parents=True)
        (directory / "attempts.tsv").write_text("attempt\n", encoding="utf-8")

    def lock_check(path: Path) -> bool:
        assert str(han.iloc[0]["unit_id"]) not in str(path)
        return False

    monkeypatch.setattr(production.simulation, "unit_lock_is_held", lock_check)
    payload = production.build_quarantine_plan(
        repo_root=root, campaign_dir=campaign_dir, expected_stale_count=2
    )
    assert payload["scope"] == "stale_eas_selected_directories_only"
    assert payload["moves_performed"] == payload["deletions_performed"] == 0
    assert {row["unit_id"] for row in payload["rows"]} == set(eas["unit_id"])
    assert (campaign_dir / "work" / str(han.iloc[0]["unit_id"])).is_dir()


def test_han_stage_audit_rejects_semantic_corruption(tmp_path: Path):
    root = tmp_path / "repo"
    root.mkdir()
    manifest = root / "manifest.json"
    seed_plan = root / "seed.tsv"
    manifest.write_text("manifest\n", encoding="utf-8")
    seed_plan.write_text("seed\n", encoding="utf-8")
    audit_path = root / "pilot.json"
    audit = _han_audit(
        stage=han_v3.STAGE_PILOT,
        units=6,
        replicates=1,
        manifest=manifest,
        seed_plan=seed_plan,
    )
    _write_json(audit_path, audit)
    loaded = production._load_han_v3_stage_audit(
        audit_path,
        repo_root=root,
        stage=han_v3.STAGE_PILOT,
        expected_units=6,
        expected_replicates=[1],
        manifest_path=manifest,
        seed_plan_path=seed_plan,
        expected_unit_ids={f"unit-{index:03d}" for index in range(6)},
        expected_seed_plan_audit={},
    )
    assert loaded["accepted_units"] == 6
    corrupt = dict(audit)
    corrupt["all_units_strict_type1_identity"] = False
    corrupt.pop("payload_sha256")
    _write_json(audit_path, _signed(corrupt))
    with pytest.raises(ValueError, match="contract differs"):
        production._load_han_v3_stage_audit(
            audit_path,
            repo_root=root,
            stage=han_v3.STAGE_PILOT,
            expected_units=6,
            expected_replicates=[1],
            manifest_path=manifest,
            seed_plan_path=seed_plan,
            expected_unit_ids={f"unit-{index:03d}" for index in range(6)},
            expected_seed_plan_audit={},
        )


def test_aggregate_readiness_binds_pilot_then_requires_final(
    tmp_path: Path, monkeypatch
):
    root = tmp_path / "repo"
    root.mkdir()
    slim = root / "slim"
    slim.write_bytes(b"slim")
    eas_manifest = root / production.DEFAULT_INTEGRATION_RELATIVE_PATH
    _write_json(eas_manifest, {"eas": True})
    eas_bundle = _bundle(root)
    monkeypatch.setattr(
        production,
        "load_selected_production_bundle",
        lambda *_args, **_kwargs: eas_bundle,
    )
    monkeypatch.setattr(
        production,
        "validate_eas_cached_units_read_only",
        lambda _bundle: {"units": 60, "cache_validated_read_only": True},
    )
    han_manifest = root / production.DEFAULT_HAN_V3_MANIFEST_RELATIVE_PATH
    _write_json(han_manifest, {"han": True})
    seed_plan = han_manifest.parent / han_v3.DEFAULT_SEED_PLAN_NAME
    seed_plan.write_text("seed\n", encoding="utf-8")
    han_manifest_payload = {
        "payload_sha256": "han-payload",
        "runtime": {"slim": "pinned"},
        "scientific_caveats": ["conditional inference"],
    }
    han_bundle = SimpleNamespace(
        manifest=han_manifest_payload,
        seed_plan_path=seed_plan,
        output_dir=han_manifest.parent,
        units=pd.DataFrame(
            {
                "unit_id": [f"unit-{index:03d}" for index in range(60)],
                "replicate_index": [1] * 6 + [2] * 54,
            }
        ),
    )
    monkeypatch.setattr(han_v3, "load_bundle", lambda *_args, **_kwargs: han_bundle)
    monkeypatch.setattr(han_v3, "_seed_plan_audit", lambda _bundle: {})
    pilot_path = root / production.DEFAULT_HAN_V3_PILOT_AUDIT_RELATIVE_PATH
    _write_json(
        pilot_path,
        _han_audit(
            stage=han_v3.STAGE_PILOT,
            units=6,
            replicates=1,
            manifest=han_manifest,
            seed_plan=seed_plan,
        ),
    )
    monkeypatch.setattr(
        han_v3,
        "audit_stage",
        lambda _bundle, stage, *, write: json.loads(
            (
                pilot_path
                if stage == han_v3.STAGE_PILOT
                else root / production.DEFAULT_HAN_V3_FINAL_AUDIT_RELATIVE_PATH
            ).read_text(encoding="utf-8")
        ),
    )
    status = production.build_aggregate_readiness(
        repo_root=root,
        eas_integration_manifest_path=eas_manifest,
        slim_path=slim,
        require_final=False,
    )
    assert status["status"] == "pilot_passed_final_pending"
    with pytest.raises(ValueError, match="all-unit audit is required"):
        production.build_aggregate_readiness(
            repo_root=root,
            eas_integration_manifest_path=eas_manifest,
            slim_path=slim,
        )

    final_path = root / production.DEFAULT_HAN_V3_FINAL_AUDIT_RELATIVE_PATH
    final_audit = _han_audit(
        stage=han_v3.STAGE_ALL,
        units=60,
        replicates=10,
        manifest=han_manifest,
        seed_plan=seed_plan,
    )
    _write_json(final_path, final_audit)
    with pytest.raises(ValueError, match="operational authorization is required"):
        production.build_aggregate_readiness(
            repo_root=root,
            eas_integration_manifest_path=eas_manifest,
            slim_path=slim,
        )
    authorization_path = root / production.DEFAULT_HAN_V3_AUTHORIZATION_RELATIVE_PATH
    _write_json(
        authorization_path,
        _han_authorization(
            manifest=han_manifest,
            seed_plan=seed_plan,
            manifest_payload=han_manifest_payload,
            final_audit=final_audit,
        ),
    )
    ready = production.build_aggregate_readiness(
        repo_root=root,
        eas_integration_manifest_path=eas_manifest,
        slim_path=slim,
    )
    assert ready["status"] == "ready_for_aggregate"
    assert ready["counts"] == {
        "eas_selected": 60,
        "han_selected": 60,
        "selected_ready": 120,
        "neutral_preserved": 600,
        "campaign_total_when_ready": 720,
    }
    (han_manifest.parent / "recovery_history").mkdir()
    with pytest.raises(ValueError, match="requires canonical authorization and audit"):
        production.build_aggregate_readiness(
            repo_root=root,
            eas_integration_manifest_path=eas_manifest,
            slim_path=slim,
        )
    (han_manifest.parent / "recovery_history").rmdir()
    output = root / production.DEFAULT_AGGREGATE_READINESS_RELATIVE_PATH
    production.write_aggregate_readiness(output, ready)
    production.write_aggregate_readiness(output, ready)
    monkeypatch.setattr(
        production, "build_aggregate_readiness", lambda **_kwargs: ready
    )
    assert (
        production.load_aggregate_readiness(output, repo_root=root, slim_path=slim)[
            "payload_sha256"
        ]
        == ready["payload_sha256"]
    )
    changed = dict(ready)
    changed["counts"] = {**ready["counts"], "selected_ready": 119}
    changed.pop("payload_sha256")
    with pytest.raises(FileExistsError):
        production.write_aggregate_readiness(output, _signed(changed))


def test_han_recovery_binding_requires_exact_additive_patch_and_history(
    tmp_path: Path, monkeypatch
):
    root = tmp_path / "repo"
    root.mkdir()
    manifest_path = root / production.DEFAULT_HAN_V3_MANIFEST_RELATIVE_PATH
    manifest_path.parent.mkdir(parents=True)
    manifest_path.write_bytes(b"manifest\n")
    seed_plan_path = manifest_path.parent / "seed_plan.tsv"
    seed_plan_path.write_bytes(b"seed-plan\n")
    original_sources = {"python/original.py": "a" * 64}
    recovery_sources = {
        han_recovery.MODULE_SOURCE_PATH: "b" * 64,
        han_recovery.WRAPPER_SOURCE_PATH: "c" * 64,
    }
    unit_rows = [
        {
            "unit_id": f"han-r{replicate:02d}-c{cell:02d}",
            "replicate_index": replicate,
        }
        for replicate in range(1, 11)
        for cell in range(6)
    ]
    direct = SimpleNamespace(
        manifest={"implementation_sources": original_sources},
        manifest_path=manifest_path,
        seed_plan_path=seed_plan_path,
        campaign_dir=root / production.DEFAULT_CAMPAIGN_RELATIVE_PATH,
        units=pd.DataFrame(unit_rows),
    )
    bundle = SimpleNamespace(direct=direct, recovery_sources=recovery_sources)
    authorized_ids = sorted(
        row["unit_id"] for row in unit_rows if row["replicate_index"] != 1
    )[: han_recovery.EXPECTED_STOPPED_RECOVERABLE_UNITS]
    history_path = (
        manifest_path.parent / "recovery_history" / "dry-plan" / "history_manifest.json"
    )
    history = _signed(
        {
            "schema": han_recovery.RECOVERY_HISTORY_SCHEMA,
            "status": "immutable_stopped_state_snapshot",
            "dry_plan_payload_sha256": "d" * 64,
            "manifest_sha256": production.sha256_file(manifest_path),
            "seed_plan_sha256": production.sha256_file(seed_plan_path),
            "original_implementation_sources": original_sources,
            "recovery_sources": recovery_sources,
            "recoverable_units": len(authorized_ids),
            "file_count": 0,
            "total_bytes": 0,
            "files": [],
            "original_evidence_moved": False,
            "original_evidence_deleted": False,
        }
    )
    _write_json(history_path, history)
    process_local_patch = {
        "input": "attempts.tsv frame passed to frozen v3 callback",
        "column": "rejection_reason",
        "change": "cast existing or new column to object dtype",
        "planned_seed_change": False,
        "retry_change": False,
        "timeout_change": False,
        "budget_change": False,
        "reconciliation_change": False,
        "scientific_model_change": False,
    }
    authorization_path = (
        root / production.DEFAULT_HAN_V3_RECOVERY_AUTHORIZATION_RELATIVE_PATH
    )
    authorization = _signed(
        {
            "schema": han_recovery.RECOVERY_AUTHORIZATION_SCHEMA,
            "status": "authorized_after_immutable_history_snapshot",
            "dry_plan_payload_sha256": history["dry_plan_payload_sha256"],
            "manifest_sha256": history["manifest_sha256"],
            "seed_plan_sha256": history["seed_plan_sha256"],
            "history_manifest_path": history_path.relative_to(root).as_posix(),
            "history_manifest_sha256": production.sha256_file(history_path),
            "history_payload_sha256": history["payload_sha256"],
            "original_implementation_sources": original_sources,
            "recovery_sources": recovery_sources,
            "recoverable_unit_ids": authorized_ids,
            "process_local_patch": process_local_patch,
        }
    )
    _write_json(authorization_path, authorization)
    completion_rows = []
    for unit_id in authorized_ids:
        completion = direct.campaign_dir / "work" / unit_id / "simulation_complete.json"
        completion.parent.mkdir(parents=True)
        completion.write_text(f"{unit_id}\n", encoding="utf-8")
        completion_rows.append(
            {
                "unit_id": unit_id,
                "completion_sha256": production.sha256_file(completion),
                "attempts_completed": 1,
            }
        )
    final_audit = {"status": "passed", "payload_sha256": "e" * 64}
    recovery_audit = _signed(
        {
            "schema": han_recovery.RECOVERY_COMPLETION_AUDIT_SCHEMA,
            "status": "passed",
            "dry_plan_payload_sha256": authorization["dry_plan_payload_sha256"],
            "manifest_sha256": authorization["manifest_sha256"],
            "seed_plan_sha256": authorization["seed_plan_sha256"],
            "history_manifest_path": authorization["history_manifest_path"],
            "history_manifest_sha256": authorization["history_manifest_sha256"],
            "history_payload_sha256": authorization["history_payload_sha256"],
            "original_implementation_sources": original_sources,
            "recovery_sources": recovery_sources,
            "canonical_authorization_path": authorization_path.relative_to(
                root
            ).as_posix(),
            "canonical_authorization_sha256": production.sha256_file(
                authorization_path
            ),
            "initial_authorized_recovery_units": len(authorized_ids),
            "pre_recovery_incomplete_unit_ids": authorized_ids,
            "post_recovery_recovered_unit_ids": authorized_ids,
            "post_recovery_incomplete_unit_ids": [],
            "validated_recovered_completions": len(authorized_ids),
            "remaining_incomplete_units": 0,
            "rows": completion_rows,
            "post_run_base_audit_status": "passed",
            "post_run_base_audit_accepted_units": 60,
            "post_run_base_audit_payload_sha256": final_audit["payload_sha256"],
            "post_run_base_audit": final_audit,
            "original_attempt_seed_timeout_budget_and_reconciliation_semantics_preserved": (
                True
            ),
        }
    )
    recovery_audit_path = root / production.DEFAULT_HAN_V3_RECOVERY_AUDIT_RELATIVE_PATH
    _write_json(recovery_audit_path, recovery_audit)
    monkeypatch.setattr(
        han_recovery, "load_recovery_bundle", lambda *_args, **_kwargs: bundle
    )
    monkeypatch.setattr(
        han_recovery, "_load_authorization", lambda *_args, **_kwargs: authorization
    )
    binding = production._load_han_v3_recovery_binding(
        repo_root=root,
        slim_path=root / "slim",
        han_manifest_path=manifest_path,
        final_audit=final_audit,
        authorization_path=authorization_path,
        audit_path=recovery_audit_path,
    )
    assert binding["recovered_units"] == len(authorized_ids) == 47
    assert binding["process_local_patch_only"] is True

    corrupted_authorization = {
        **authorization,
        "process_local_patch": {
            **process_local_patch,
            "planned_seed_change": True,
        },
    }
    monkeypatch.setattr(
        han_recovery,
        "_load_authorization",
        lambda *_args, **_kwargs: corrupted_authorization,
    )
    with pytest.raises(ValueError, match="dtype-recovery audit binding differs"):
        production._load_han_v3_recovery_binding(
            repo_root=root,
            slim_path=root / "slim",
            han_manifest_path=manifest_path,
            final_audit=final_audit,
            authorization_path=authorization_path,
            audit_path=recovery_audit_path,
        )


def test_parser_exposes_only_exact_eas_runner_not_cell_filters():
    parser = production.build_parser()
    help_text = parser.format_help()
    assert "simulate-eas" in help_text
    assert "analyze" in help_text
    assert "report" in help_text
    assert "--af" not in help_text
    assert "--unit" not in help_text
    assert "--selection-coefficient" not in help_text
    assert "fixed-v2" not in help_text


def test_analysis_and_report_use_only_readiness_safe_routes(
    tmp_path: Path, monkeypatch
):
    root = tmp_path / "repo"
    root.mkdir()
    aggregate_outputs = []
    for name in (
        "classes.tsv",
        "pairs.tsv",
        "truth.tsv",
        "spatial.tsv",
        "truth_spatial.tsv",
    ):
        path = root / name
        path.write_text("test\n", encoding="utf-8")
        aggregate_outputs.append(path)
    analysis_output = root / "analysis.json"
    analysis_output.write_text("{}\n", encoding="utf-8")
    report_path = root / "report.md"

    def unsafe_standard_aggregate(*_args, **_kwargs):
        raise AssertionError("standard unsafe aggregator must never be called")

    def safe_aggregate(**_kwargs):
        return tuple(aggregate_outputs)

    def fake_report(_campaign_dir):
        report_path.write_text("report\n", encoding="utf-8")
        return report_path

    monkeypatch.setattr(
        campaign, "_aggregate_compact_decode", unsafe_standard_aggregate
    )
    monkeypatch.setattr(production, "aggregate_with_readiness", safe_aggregate)
    monkeypatch.setattr(
        focused_analysis,
        "run_focused_analysis",
        lambda *_args, **_kwargs: {"completion": analysis_output},
    )
    monkeypatch.setattr(
        focused_analysis,
        "run_truth_gamma_comparison",
        lambda *_args, **_kwargs: {"completion": analysis_output},
    )
    monkeypatch.setattr(focused_report, "generate_run_report", fake_report)
    analyzed = production.analyze_with_readiness(
        repo_root=root,
        slim_path=root / "slim",
    )
    assert analyzed["status"] == "analyzed_and_reported"
    assert analyzed["all_simulation_cache_validation_dispatched_through_readiness"]

    monkeypatch.setattr(
        production,
        "load_aggregate_readiness",
        lambda *_args, **_kwargs: {"payload_sha256": "ready"},
    )
    reported = production.report_with_readiness(
        repo_root=root,
        slim_path=root / "slim",
    )
    assert reported["status"] == "reported"
    assert reported["unsafe_standard_reaggregation_used"] is False
    assert reported["report_sha256"] == production.sha256_file(report_path)


def test_read_only_eas_cache_validator_detects_missing_completion(tmp_path: Path):
    root = tmp_path / "repo"
    root.mkdir()
    _execution_units(root)
    bundle = _bundle(root)
    with pytest.raises(ValueError, match="completion is absent"):
        production.validate_eas_cached_units_read_only(bundle)


def test_eas_supervisor_schedule_is_exact_and_deterministic():
    unit = _selected_unit()
    schedule = production._eas_supervisor_schedule(unit)
    assert len(schedule) == production.DEFAULT_EAS_MAX_DRAWS == 68
    assert list(schedule["attempt_zero_based"]) == list(range(68))
    assert schedule["simulation_seed"].is_unique
    assert schedule["panel_seed"].is_unique
    assert set(schedule["simulation_seed"]).isdisjoint(schedule["panel_seed"])
    assert schedule.equals(production._eas_supervisor_schedule(unit))


@pytest.mark.parametrize(
    ("terminal_reason", "terminal_accepted"),
    [
        ("candidate_pool_af_outside_prespecified_band", False),
        (None, True),
    ],
    ids=("timeout_then_reject", "timeout_then_accept"),
)
def test_eas_timeout_reason_cast_handles_pandas3_terminal_sequences(
    tmp_path: Path,
    monkeypatch,
    terminal_reason: str | None,
    terminal_accepted: bool,
):
    root = tmp_path / "repo"
    root.mkdir()
    bundle = _bundle(root)
    supervisor_source = root / han_v3.MODULE_SOURCE_PATH
    supervisor_source.parent.mkdir(parents=True, exist_ok=True)
    supervisor_source.write_bytes(b"hardened supervisor source\n")
    captured: list[pd.DataFrame] = []

    class FakeSupervisor:
        def __init__(self, _engine, *, unit_dir, schedule):
            self.records = {
                attempt: {
                    "schema": han_v3.SUPERVISOR_LEDGER_SCHEMA,
                    "execution_id": f"{attempt + 1:064x}",
                    "supervisor_status": ("timed_out" if attempt == 0 else "complete"),
                }
                for attempt in schedule["attempt_zero_based"].astype(int)
            }

        def mark_base_attempt_ledger_committed(self, frame):
            captured.append(frame.copy())

        def finalize_unit_completion(self):
            return None

    monkeypatch.setattr(han_v3, "_SupervisedSlimEngine", FakeSupervisor)
    monkeypatch.setattr(
        simulation.stdpopsim,
        "get_engine",
        lambda name, *_args, **_kwargs: object(),
    )
    unit = _selected_unit()
    unit_dir = bundle.campaign_dir / "work" / str(unit["unit_id"])
    timeout = pd.DataFrame(
        {
            "attempt_zero_based": [0],
            "seed": [
                int(
                    production._eas_supervisor_schedule(unit).iloc[0]["simulation_seed"]
                )
            ],
            "status": ["timed_out"],
            "accepted": [False],
            "rejection_reason": pd.Series([float("nan")], dtype="float64"),
        }
    )
    terminal = pd.DataFrame(
        {
            "attempt_zero_based": [0, 1],
            "seed": production._eas_supervisor_schedule(unit)
            .iloc[:2]["simulation_seed"]
            .astype(int)
            .tolist(),
            "status": ["timed_out", "complete"],
            "accepted": [False, terminal_accepted],
        }
    )
    terminal["rejection_reason"] = (
        pd.Series([float("nan"), terminal_reason], dtype="object")
        if terminal_reason is not None
        else pd.Series([float("nan"), float("nan")], dtype="float64")
    )
    with production._patched_eas_process_supervisor(bundle, unit):
        simulation._atomic_frame(unit_dir / "attempts.tsv", timeout)
        simulation._atomic_frame(unit_dir / "attempts.tsv", terminal)

    assert len(captured) == 2
    assert captured[0].iloc[0]["rejection_reason"] == "draw_timeout"
    assert captured[1].iloc[0]["rejection_reason"] == "draw_timeout"
    observed_terminal = captured[1].iloc[1]["rejection_reason"]
    if terminal_reason is None:
        assert pd.isna(observed_terminal)
    else:
        assert observed_terminal == terminal_reason


@pytest.mark.skipif(os.name != "posix", reason="requires Linux/WSL process groups")
def test_eas_attempt_audit_strictly_binds_supervision_and_reconciliation(
    tmp_path: Path, monkeypatch
):
    root = tmp_path / "repo"
    root.mkdir()
    bundle = _bundle(root)
    unit = _selected_unit()
    unit_dir, completion = _write_valid_eas_restart_ledgers(bundle, unit, monkeypatch)
    audit = production._validate_eas_attempt_provenance(bundle, unit, completion)
    assert audit["attempts_completed"] == 2
    assert audit["crash_reconciliation_rows"] == 1
    assert audit["reconciliation_sha256"] == production.sha256_file(
        unit_dir / "supervision_reconciliation.tsv"
    )
    assert audit["supervisor_fields_equal_attempt_ledger"] is True

    attempts_path = unit_dir / "attempts.tsv"
    original_attempts = attempts_path.read_bytes()
    attempts = pd.read_csv(attempts_path, sep="\t")
    attempts["eas_process_group_gone"] = attempts["eas_process_group_gone"].astype(
        object
    )
    attempts.loc[0, "eas_process_group_gone"] = "malformed"
    attempts.to_csv(attempts_path, sep="\t", index=False, lineterminator="\n")
    with pytest.raises(ValueError, match="boolean eas_process_group_gone"):
        production._validate_eas_attempt_provenance(bundle, unit, completion)
    attempts_path.write_bytes(original_attempts)

    supervision_path = unit_dir / "supervision.tsv"
    original_supervision = supervision_path.read_bytes()
    supervision = pd.read_csv(supervision_path, sep="\t")
    supervision.loc[0, "execution_id"] = "f" * 64
    supervision.to_csv(supervision_path, sep="\t", index=False, lineterminator="\n")
    with pytest.raises(ValueError, match="supervisor field execution_id differs"):
        production._validate_eas_attempt_provenance(bundle, unit, completion)
    supervision_path.write_bytes(original_supervision)

    reconciliation_path = unit_dir / "supervision_reconciliation.tsv"
    reconciliation = pd.read_csv(reconciliation_path, sep="\t")
    reconciliation["process_group_gone"] = reconciliation["process_group_gone"].astype(
        object
    )
    reconciliation.loc[0, "process_group_gone"] = "malformed"
    reconciliation.to_csv(
        reconciliation_path, sep="\t", index=False, lineterminator="\n"
    )
    with pytest.raises(ValueError, match="boolean process_group_gone"):
        production._validate_eas_attempt_provenance(bundle, unit, completion)

    attempts_path.unlink()
    with pytest.raises(ValueError, match="attempt provenance is absent"):
        production._validate_eas_attempt_provenance(bundle, unit, completion)


@pytest.mark.skipif(os.name != "posix", reason="requires Linux/WSL process groups")
def test_eas_adapter_supervisor_reaps_real_grandchild_on_timeout(
    tmp_path: Path, monkeypatch
):
    root = tmp_path / "repo"
    root.mkdir()
    bundle = _bundle(root)
    unit = _selected_unit()
    monkeypatch.setattr(
        simulation.stdpopsim,
        "get_engine",
        lambda name, *_args, **_kwargs: _HangingEngineWithGrandchild(),
    )
    schedule = production._eas_supervisor_schedule(unit).set_index("attempt_zero_based")
    started = time.monotonic()
    with (
        pytest.raises(simulation.SimulationDrawTimeout),
        production._patched_eas_process_supervisor(bundle, unit) as supervisor,
        simulation._wall_clock_timeout(0.4),
    ):
        supervisor.simulate(
            seed=int(schedule.loc[0, "simulation_seed"]),
            logfile=str(
                bundle.campaign_dir
                / "work"
                / str(unit["unit_id"])
                / "slim_draw_000.csv"
            ),
        )
    assert time.monotonic() - started < 7.0
    record = supervisor.records[0]
    assert record["supervisor_status"] == "timed_out"
    assert record["direct_child_reaped"] is True
    assert record["process_group_gone"] is True
    assert record["sigterm_sent"] is True
    assert record["sigkill_sent"] is True
    grandchild_pid = int(
        (
            bundle.campaign_dir
            / "work"
            / str(unit["unit_id"])
            / "slim_draw_000.csv.grandchild.pid"
        ).read_text(encoding="utf-8")
    )
    with pytest.raises(ProcessLookupError):
        os.kill(grandchild_pid, 0)


@pytest.mark.skipif(os.name != "posix", reason="requires Linux/WSL process groups")
def test_accepted_precompletion_crash_reconciles_only_known_partial_artifacts(
    tmp_path: Path,
):
    unit_dir = tmp_path / "unit"
    unit_dir.mkdir()
    context = multiprocessing.get_context("fork")
    worker = context.Process(
        target=_write_accepted_partial_worker, args=(str(unit_dir),)
    )
    worker.start()
    worker.join(timeout=10.0)
    assert not worker.is_alive() and worker.exitcode == 0
    assert list(unit_dir.glob(".han_v3_active_*.json"))
    schedule = pd.DataFrame(
        [{"attempt_zero_based": 0, "simulation_seed": 12345, "panel_seed": 54321}]
    )
    han_v3._SupervisedSlimEngine(
        _FastTreeEngine(),
        unit_dir=unit_dir,
        schedule=schedule,
        term_grace_seconds=0.1,
    )
    attempts = pd.read_csv(unit_dir / "attempts.tsv", sep="\t")
    assert str(attempts.iloc[0]["accepted"]).casefold() == "false"
    assert attempts.iloc[0]["status"] == "timed_out"
    assert (
        attempts.iloc[0]["rejection_reason"]
        == "worker_interrupted_after_acceptance_before_unit_completion"
    )
    for name in ("simulation.trees", "sample_manifest.tsv", "truth_profiles.tsv.gz"):
        assert not (unit_dir / name).exists()
    assert not (unit_dir / "pairs").exists()
    assert not list(unit_dir.glob(".han_v3_active_*.json"))
    reconciliation = pd.read_csv(unit_dir / "supervision_reconciliation.tsv", sep="\t")
    assert (
        reconciliation.iloc[0]["reason"]
        == "accepted_attempt_interrupted_before_unit_completion"
    )
    assert bool(reconciliation.iloc[0]["process_group_gone"])


@contextmanager
def _no_op_context(_bundle):
    yield
