from __future__ import annotations

import json
import multiprocessing
import os
import subprocess
import time
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
from gamma_smc_aou import focused_han_direct_v3 as v3
from gamma_smc_aou import focused_selection_campaign as focused_campaign
from gamma_smc_aou import focused_selection_simulation as simulation


class _HangingEngineWithGrandchild:
    def simulate(self, *_args, **kwargs):
        Path(str(kwargs["logfile"]) + ".wrapper.pid").write_text(
            str(os.getpid()), encoding="utf-8"
        )
        child = subprocess.Popen(["/bin/sh", "-c", "trap '' TERM; sleep 60"])
        Path(str(kwargs["logfile"]) + ".grandchild.pid").write_text(
            str(child.pid), encoding="utf-8"
        )
        while True:
            time.sleep(1)


def _run_hanging_proxy_worker(unit_dir_raw: str) -> None:
    unit_dir = Path(unit_dir_raw)
    schedule = pd.DataFrame(
        [{"attempt_zero_based": 0, "simulation_seed": 12345, "panel_seed": 54321}]
    )
    proxy = v3._SupervisedSlimEngine(
        _HangingEngineWithGrandchild(),
        unit_dir=unit_dir,
        schedule=schedule,
        term_grace_seconds=0.2,
    )
    with v3._supervised_draw_timeout(60.0):
        proxy.simulate(
            seed=12345,
            logfile=str(unit_dir / "slim_draw_000.csv"),
        )


def _execution_units(root: Path) -> pd.DataFrame:
    rows = []
    for coefficient in v3.SELECTION_COEFFICIENTS:
        for target in v3.TARGET_FREQUENCIES:
            for replicate in range(1, v3.REPLICATES_PER_CELL + 1):
                rows.append(
                    {
                        "unit_id": (
                            f"ancient_eurasia_han_introgression__af{round(100 * target)}"
                            + f"__selected_s{coefficient:.3f}".replace(".", "p")
                            + f"__rep{replicate:03d}"
                        ),
                        "demography_id": "ancient_eurasia_han_introgression",
                        "demography_kind": "introgression",
                        "population": "Han",
                        "source_model": "stdpopsim AncientEurasia_9K19",
                        "selection_origin": "archaic_specific_introgressed_standing_variation",
                        "simulation_class": "selected",
                        "selection_coefficient": coefficient,
                        "target_allele_frequency": target,
                        "population_af_lower": target - v3.AF_HALF_WIDTH,
                        "population_af_upper": target + v3.AF_HALF_WIDTH,
                        "exact_sample_alt_count": round(2 * 100 * target),
                        "sample_diploids": 100,
                        "replicate_index": replicate,
                        "seed": 1000 + len(rows),
                        "sequence_length_bp": 10_000_000,
                        "focal_position_bp": 5_000_000,
                    }
                )
    frame = pd.DataFrame(rows)
    path = root / "focused_selection_EAS_sim/execution_units.tsv"
    path.parent.mkdir(parents=True)
    frame.to_csv(path, sep="\t", index=False)
    return frame


def _write_seed(path: Path, seed: int, *, panel_seed: int | None = None) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame({"seed": [seed]})
    if panel_seed is not None:
        frame["panel_seed"] = [panel_seed]
    frame.to_csv(path, sep="\t", index=False)
    return path


def _fake_plan_inputs(root: Path):
    units = _execution_units(root)
    campaign = root / "focused_selection_EAS_sim"
    seed_paths = [campaign / "execution_units.tsv"]
    seed_paths.append(
        _write_seed(campaign / "calibration/han_selection_end_plan.tsv", 80_001)
    )
    seed_paths.append(
        _write_seed(
            campaign
            / "calibration/han_fixed_cessation_v2/han_fixed_v2_all_phase_plan.tsv",
            80_002,
            panel_seed=80_003,
        )
    )
    seed_paths.append(
        _write_seed(campaign / "calibration/eas_selected/screen20_plan.tsv", 80_004)
    )
    seed_paths.append(
        _write_seed(
            campaign / "calibration/exploratory_seed_provenance/"
            "exploratory_seed_reservation.tsv",
            80_005,
        )
    )
    failed = (
        campaign / "calibration/han_fixed_cessation_v2/confirmation_run_snapshot.tsv"
    )
    failed.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        {
            "status": ["complete", "failed"],
            "elapsed_seconds": [91.0, None],
            "error": ["", "SLiM exited with code -15"],
        }
    ).to_csv(failed, sep="\t", index=False)
    summary = failed.parent.parent / "quarantine_inventory.json"
    summary.write_text(
        json.dumps(
            {
                "schema": "gamma-smc.han-fixed-cessation-v2-retry-quarantine/v1",
                "status": "canonical_noninferential_retry_provenance",
                "inference_use": False,
                "work_evidence": {
                    "started_directories": 2,
                    "valid_completions": 1,
                    "sigterm_failure_records": 1,
                    "failure_records": 1,
                    "partial_directories": 0,
                },
            }
        ),
        encoding="utf-8",
    )
    slim = root / "bin/slim"
    slim.parent.mkdir(parents=True)
    slim.write_bytes(b"SLiM 4.2.2")
    return units, campaign, seed_paths, failed, summary, slim


def _runtime(_root: Path, slim: Path):
    return {
        "portable_contract": {
            "test": True,
            "slim": {"relative_path": slim.relative_to(_root).as_posix()},
        },
        "host_diagnostics_not_compared": {"root": str(_root)},
    }


@pytest.fixture
def planned_bundle(tmp_path: Path, monkeypatch):
    root = tmp_path / "repo"
    root.mkdir()
    units, campaign, seed_paths, failed, summary, slim = _fake_plan_inputs(root)
    monkeypatch.setattr(v3, "_source_records", lambda _root: {"source": "bound"})
    monkeypatch.setattr(v3, "_runtime_record", _runtime)
    monkeypatch.setattr(
        v3,
        "CANONICAL_SEED_ARTIFACT_ALLOWLIST",
        tuple(path.relative_to(root).as_posix() for path in seed_paths),
    )
    monkeypatch.setattr(v3, "_require_git_tracked", lambda _paths, _root: None)
    manifest, plan, registry, inventory = v3.build_plan(
        repo_root=root,
        campaign_dir=campaign,
        slim_path=slim,
        failed_v2_evidence_path=failed,
        failed_v2_summary_path=summary,
    )
    output = root / v3.DEFAULT_OUTPUT_RELATIVE_PATH
    v3.publish_plan_bundle(
        output,
        repo_root=root,
        manifest=manifest,
        seed_plan=plan,
        excluded_seed_registry=registry,
        excluded_seed_inventory=inventory,
    )
    bundle = v3.load_bundle(
        output / v3.DEFAULT_PLAN_MANIFEST_NAME,
        repo_root=root,
        slim_path=slim,
    )
    return root, units, campaign, output, bundle


def test_plan_allocates_exact_final_sixty_unit_contract(planned_bundle):
    _, _, _, _, bundle = planned_bundle
    plan = bundle.seed_plan
    assert len(plan) == 60 * 300
    assert plan.groupby("unit_id").size().eq(300).all()
    assert plan["stage"].value_counts().to_dict() == {
        "simulate": 54 * 300,
        "pilot": 6 * 300,
    }
    all_seeds = pd.concat([plan["simulation_seed"], plan["panel_seed"]])
    assert len(all_seeds) == 36_000
    assert not all_seeds.duplicated().any()
    assert bundle.manifest["design"] == {
        **bundle.manifest["design"],
        "pilot_units": 6,
        "remaining_simulation_units": 54,
        "attempts_per_unit": 300,
        "draw_timeout_seconds": 300.0,
        "cumulative_timeout_seconds": 90_000.0,
    }
    assert bundle.manifest["model"]["strict_selected_type1_identity"] is True
    assert bundle.manifest["model"]["gamma_smc_statistics_used"] is False
    assert bundle.manifest["seed_exclusion_policy"] == {
        **bundle.manifest["seed_exclusion_policy"],
        "canonical_allowlist_is_source_controlled": True,
        "recursive_discovery_used": False,
        "load_validates_frozen_bundle_only": True,
    }
    assert (
        len(bundle.excluded_seed_registry)
        == bundle.manifest["seed_registry"]["excluded_union_count"]
    )


def test_runtime_source_binding_includes_transitive_bootstraps():
    repo_root = Path(__file__).resolve().parents[1]
    sources = v3._source_records(repo_root)
    assert "scripts/bootstrap_focused_han_direct_v3.sh" in sources
    assert "scripts/bootstrap_eas_sweep_study.sh" in sources
    assert "scripts/bootstrap_uv.sh" in sources
    assert all(len(sources[path]) == 64 for path in sources)


def test_seed_registry_excludes_every_named_namespace(planned_bundle):
    _, _, _, _, bundle = planned_bundle
    namespaces = bundle.manifest["seed_registry"]["excluded_namespaces"]
    assert set(namespaces) == {
        "production",
        "fallback",
        "v1",
        "v2",
        "eas",
        "exploratory",
    }
    assert all(record["count"] > 0 for record in namespaces.values())
    assert bundle.manifest["seed_registry"][
        "all_attempt_and_panel_seeds_allocated_atomically"
    ]
    assert bundle.manifest["seed_registry"]["disjoint_from_all_excluded_namespaces"]


def test_failed_v2_is_compute_budget_evidence_only(planned_bundle):
    _, _, _, _, bundle = planned_bundle
    evidence = bundle.manifest["failed_v2_pilot_evidence"]
    assert evidence["status_counts"] == {"complete": 1, "failed": 1}
    assert evidence["use"] == "pilot_compute_budget_evidence_only"
    assert evidence["inferential_use"] is False
    assert evidence["endpoint_authorization"] is False
    assert evidence["production_unit_acceptance_use"] is False
    assert evidence["gamma_smc_statistics_used"] is False
    assert evidence["snapshot"]["rows"] == 2
    assert evidence["compact_summary"]["work_evidence"]["started_directories"] == 2


def test_process_local_patch_uses_preallocated_seeds_and_restores(
    planned_bundle, monkeypatch
):
    _, _, _, _, bundle = planned_bundle
    unit = bundle.units.iloc[0].to_dict()
    schedule = v3._unit_seed_schedule(bundle, unit["unit_id"])
    original_objects = simulation._selected_objects
    original_seed = simulation._stable_seed
    original_frame = simulation._atomic_frame
    original_paths = simulation.IMPLEMENTATION_PATHS
    original_timeout = simulation._wall_clock_timeout
    original_get_engine = simulation.stdpopsim.get_engine

    fake_sweep = SimpleNamespace(provenance_record=lambda: {"fixed": True})

    def fake_objects(*_args, **kwargs):
        assert (
            kwargs["han_selection_end_generations_ago"]
            == schedule.iloc[0]["realized_selection_end_generations_ago"]
        )
        return "model", fake_sweep, {"Han": 500}, {"selection_endpoint": "fixed"}

    monkeypatch.setattr(simulation, "_selected_objects", fake_objects)
    monkeypatch.setattr(v3.fixed_v2, "_audit_fixed_events", lambda *_a, **_k: {})

    class FakeSupervisor:
        def __init__(self, *_args, **_kwargs):
            self.records = {
                0: {
                    "schema": v3.SUPERVISOR_LEDGER_SCHEMA,
                    "execution_id": "a" * 64,
                    "supervisor_status": "complete",
                    "supervisor_pid": 10,
                    "supervisor_process_identity": {
                        "pid": 10,
                        "boot_id": "test",
                        "start_ticks": 1,
                    },
                    "supervisor_started_unix_ns": 1,
                    "child_pid": 11,
                    "child_pgid": 11,
                    "child_sid": 11,
                    "child_process_identity": {
                        "pid": 11,
                        "boot_id": "test",
                        "start_ticks": 1,
                    },
                    "child_identity_observed": True,
                    "engine_wrapper_pid": 12,
                    "engine_wrapper_process_identity": {
                        "pid": 12,
                        "boot_id": "test",
                        "start_ticks": 1,
                    },
                    "parent_death_watchdog_armed": True,
                    "timeout_seconds": 300.0,
                    "term_grace_seconds": 5.0,
                    "sigterm_sent": False,
                    "sigkill_sent": False,
                    "direct_child_reaped": True,
                    "adopted_children_reaped": 0,
                    "process_group_gone": True,
                    "temp_directory": ".temp",
                    "temp_tree_created": True,
                    "temp_tree_loaded": True,
                    "temp_tree_removed": True,
                    "supervisor_elapsed_seconds": 1.0,
                }
            }

        def mark_base_attempt_ledger_committed(self, _frame):
            return None

        def finalize_unit_completion(self):
            return None

    monkeypatch.setattr(v3, "_SupervisedSlimEngine", FakeSupervisor)
    fake_original = simulation._selected_objects
    ledger_path = bundle.repo_root / "attempts.tsv"
    binding = v3._unit_binding(bundle, unit)
    with v3.patched_v3_executor(bundle, unit):
        assert simulation._stable_seed(unit["seed"], "selected:0") == int(
            schedule.iloc[0]["simulation_seed"]
        )
        assert simulation._stable_seed(unit["seed"], "selected-panel:0") == int(
            schedule.iloc[0]["panel_seed"]
        )
        _, _, _, provenance = simulation._selected_objects(
            unit,
            bundle.repo_root,
            5.0,
            500,
            han_selection_calibration=binding,
        )
        assert provenance["han_direct_production_v3"]["direct_operational_estimand"]
        simulation._atomic_frame(
            ledger_path,
            pd.DataFrame(
                [
                    {
                        "attempt_zero_based": 0,
                        "seed": int(schedule.iloc[0]["simulation_seed"]),
                        "status": "complete",
                        "accepted": False,
                        "rejection_reason": "outside_band",
                    }
                ]
            ),
        )
        ledger = pd.read_csv(ledger_path, sep="\t")
        assert int(ledger.iloc[0]["panel_seed"]) == int(schedule.iloc[0]["panel_seed"])
        assert ledger.iloc[0]["attempt_ledger_schema"] == v3.ATTEMPT_LEDGER_SCHEMA
    assert simulation._selected_objects is fake_original
    assert simulation._stable_seed is original_seed
    assert simulation._atomic_frame is original_frame
    assert simulation.IMPLEMENTATION_PATHS is original_paths
    assert simulation._wall_clock_timeout is original_timeout
    assert simulation.stdpopsim.get_engine is original_get_engine
    monkeypatch.setattr(simulation, "_selected_objects", original_objects)


def test_simulate_stage_refuses_to_start_without_passed_pilot(planned_bundle):
    _, _, _, _, bundle = planned_bundle
    with pytest.raises(ValueError, match="passed pilot audit"):
        v3.run_stage(bundle, v3.STAGE_SIMULATE, workers=1)


def test_quarantine_plan_is_read_only_and_han_selected_only(
    planned_bundle, monkeypatch
):
    root, _, campaign, _, _ = planned_bundle
    unit_ids = pd.read_csv(campaign / "execution_units.tsv", sep="\t")["unit_id"].head(
        2
    )
    for index, unit_id in enumerate(unit_ids):
        directory = campaign / "work" / unit_id
        directory.mkdir(parents=True)
        (directory / "attempts.tsv").write_text(f"attempt\n{index}\n", encoding="utf-8")
        if index == 0:
            (directory / "simulation_complete.json").write_text(
                json.dumps(
                    {
                        "provenance": {
                            "han_direct_production_v3": {
                                "direct_operational_estimand": True
                            }
                        }
                    }
                ),
                encoding="utf-8",
            )
    neutral = campaign / "work/eas_phlash_median__af10__neutral__rep001"
    neutral.mkdir(parents=True)
    (neutral / "keep.txt").write_text("neutral\n", encoding="utf-8")
    monkeypatch.setattr(simulation, "unit_lock_is_held", lambda _path: False)
    with pytest.raises(ValueError, match="nonnegative"):
        v3.build_quarantine_plan(
            repo_root=root,
            campaign_dir=campaign,
            expected_stale_count=2.5,
        )
    payload = v3.build_quarantine_plan(
        repo_root=root,
        campaign_dir=campaign,
        expected_stale_count=2,
    )
    assert payload["moves_performed"] == payload["deletions_performed"] == 0
    assert payload["neutral_directories_touched"] == 0
    assert len(payload["rows"]) == 2
    assert all(Path(root / row["source"]).is_dir() for row in payload["rows"])
    assert all(
        not Path(root / row["proposed_destination"]).exists() for row in payload["rows"]
    )
    assert neutral.is_dir()


def test_manifest_and_seed_plan_tampering_fail_closed(planned_bundle):
    root, _, _, output, bundle = planned_bundle
    tampered = pd.read_csv(bundle.seed_plan_path, sep="\t")
    tampered.loc[0, "panel_seed"] = tampered.loc[0, "simulation_seed"]
    tampered.to_csv(bundle.seed_plan_path, sep="\t", index=False)
    with pytest.raises(ValueError, match="checksum"):
        v3.load_bundle(
            output / v3.DEFAULT_PLAN_MANIFEST_NAME,
            repo_root=root,
            slim_path=bundle.slim_path,
        )


def test_load_uses_frozen_seed_registry_not_mutable_recursive_discovery(
    planned_bundle, monkeypatch
):
    root, _, campaign, output, bundle = planned_bundle
    local = _write_seed(campaign / "calibration/new_local_extra.tsv", 999_999_999)
    assert local.is_file()
    monkeypatch.setattr(
        v3,
        "canonical_seed_artifacts",
        lambda **_kwargs: (_ for _ in ()).throw(AssertionError("must not discover")),
    )
    reloaded = v3.load_bundle(
        output / v3.DEFAULT_PLAN_MANIFEST_NAME,
        repo_root=root,
        slim_path=bundle.slim_path,
    )
    assert len(reloaded.excluded_seed_registry) == len(bundle.excluded_seed_registry)


def test_load_compares_portable_runtime_not_host_diagnostics(
    planned_bundle, monkeypatch
):
    root, _, _, output, bundle = planned_bundle
    portable = bundle.manifest["runtime"]["portable_contract"]
    monkeypatch.setattr(
        v3,
        "_runtime_record",
        lambda _root, _slim: {
            "portable_contract": portable,
            "host_diagnostics_not_compared": {
                "python_executable_path": "/different/fresh/clone/.venv/bin/python",
                "platform_release": "different-kernel",
            },
        },
    )
    reloaded = v3.load_bundle(
        output / v3.DEFAULT_PLAN_MANIFEST_NAME,
        repo_root=root,
        slim_path=bundle.slim_path,
    )
    assert reloaded.manifest["runtime"]["portable_contract"] == portable


def test_local_seed_evidence_must_be_subset_of_canonical_registry(planned_bundle):
    root, _, campaign, _, bundle = planned_bundle
    included_seed = int(bundle.excluded_seed_registry.iloc[0]["seed"])
    subset = _write_seed(campaign / "calibration/local_subset.tsv", included_seed)
    records = v3._validate_local_seed_subset(
        [subset],
        repo_root=root,
        excluded=set(bundle.excluded_seed_registry["seed"].astype(int)),
    )
    assert records[0]["subset_of_frozen_registry"] is True
    outside = _write_seed(campaign / "calibration/local_outside.tsv", 2_147_483_646)
    with pytest.raises(ValueError, match="subset"):
        v3._validate_local_seed_subset(
            [outside],
            repo_root=root,
            excluded=set(bundle.excluded_seed_registry["seed"].astype(int)),
        )


def test_schedule_validation_covers_endpoint_stage_attempt_and_production_seed(
    planned_bundle,
):
    _, _, _, _, bundle = planned_bundle
    tampered = bundle.seed_plan.copy()
    tampered.loc[0, "production_unit_seed"] += 1
    with pytest.raises(ValueError, match="production_unit_seed"):
        v3._validate_seed_plan(
            tampered,
            units=bundle.units,
            excluded=set(bundle.excluded_seed_registry["seed"].astype(int)),
        )
    nonintegral = bundle.seed_plan.copy()
    nonintegral["simulation_seed"] = nonintegral["simulation_seed"].astype(float)
    nonintegral.loc[0, "simulation_seed"] += 0.5
    with pytest.raises(ValueError, match="exact integers"):
        v3._validate_seed_plan(
            nonintegral,
            units=bundle.units,
            excluded=set(bundle.excluded_seed_registry["seed"].astype(int)),
        )


def test_execution_unit_integrality_and_unit_id_correspondence(tmp_path: Path):
    root = tmp_path / "repo"
    root.mkdir()
    units = _execution_units(root)
    v3._validate_production_units(units)
    nonintegral = units.copy()
    nonintegral["replicate_index"] = nonintegral["replicate_index"].astype(float)
    nonintegral.loc[0, "replicate_index"] = 1.5
    with pytest.raises(ValueError, match="exact integers"):
        v3._validate_production_units(nonintegral)
    wrong_id = units.copy()
    wrong_id.loc[0, "unit_id"] = wrong_id.loc[0, "unit_id"].replace("af10", "af11")
    with pytest.raises(ValueError, match="unit IDs"):
        v3._validate_production_units(wrong_id)


def _write_valid_attempt_ledgers(bundle, unit) -> None:
    schedule = v3._unit_seed_schedule(bundle, unit["unit_id"])
    row = schedule.iloc[0]
    unit_dir = bundle.campaign_dir / "work" / unit["unit_id"]
    unit_dir.mkdir(parents=True, exist_ok=True)
    attempt_row = {
        "attempt_zero_based": 0,
        "seed": int(row["simulation_seed"]),
        "planned_simulation_seed": int(row["simulation_seed"]),
        "panel_seed": int(row["panel_seed"]),
        "panel_seed_used": True,
        "accepted": True,
        "rejection_reason": "",
        "status": "complete",
        "timeout_seconds": 300.0,
        "elapsed_seconds": 1.0,
        "attempt_ledger_schema": v3.ATTEMPT_LEDGER_SCHEMA,
        "seed_plan_sha256": v3.sha256_file(bundle.seed_plan_path),
        "supervisor_schema": v3.SUPERVISOR_LEDGER_SCHEMA,
        "supervisor_execution_id": "a" * 64,
        "supervisor_status": "complete",
        "supervisor_pid": 100,
        "supervisor_process_identity": (
            "{'pid': 100, 'boot_id': 'test', 'start_ticks': 1}"
        ),
        "supervisor_started_unix_ns": 1,
        "child_pid": 101,
        "child_pgid": 101,
        "child_sid": 101,
        "child_process_identity": ("{'pid': 101, 'boot_id': 'test', 'start_ticks': 1}"),
        "child_identity_observed": True,
        "engine_wrapper_pid": 102,
        "engine_wrapper_process_identity": (
            "{'pid': 102, 'boot_id': 'test', 'start_ticks': 1}"
        ),
        "parent_death_watchdog_armed": True,
        "supervisor_timeout_seconds": 300.0,
        "supervisor_term_grace_seconds": 5.0,
        "supervisor_sigterm_sent": False,
        "supervisor_sigkill_sent": False,
        "direct_child_reaped": True,
        "adopted_children_reaped": 0,
        "process_group_gone": True,
        "attempt_temp_directory": ".temp",
        "attempt_temp_tree_created": True,
        "attempt_temp_tree_loaded": True,
        "attempt_temp_tree_removed": True,
        "supervisor_elapsed_seconds": 1.0,
    }
    pd.DataFrame([attempt_row]).to_csv(unit_dir / "attempts.tsv", sep="\t", index=False)
    supervisor_row = {
        "schema": v3.SUPERVISOR_LEDGER_SCHEMA,
        "execution_id": "a" * 64,
        "attempt_zero_based": 0,
        "simulation_seed": int(row["simulation_seed"]),
        "supervisor_status": "complete",
        "supervisor_pid": 100,
        "supervisor_process_identity": (
            "{'pid': 100, 'boot_id': 'test', 'start_ticks': 1}"
        ),
        "supervisor_started_unix_ns": 1,
        "child_pid": 101,
        "child_pgid": 101,
        "child_sid": 101,
        "child_process_identity": ("{'pid': 101, 'boot_id': 'test', 'start_ticks': 1}"),
        "child_identity_observed": True,
        "engine_wrapper_pid": 102,
        "engine_wrapper_process_identity": (
            "{'pid': 102, 'boot_id': 'test', 'start_ticks': 1}"
        ),
        "parent_death_watchdog_armed": True,
        "timeout_seconds": 300.0,
        "term_grace_seconds": 5.0,
        "sigterm_sent": False,
        "sigkill_sent": False,
        "direct_child_reaped": True,
        "adopted_children_reaped": 0,
        "process_group_gone": True,
        "temp_directory": ".temp",
        "temp_tree_created": True,
        "temp_tree_loaded": True,
        "temp_tree_removed": True,
        "supervisor_elapsed_seconds": 1.0,
    }
    pd.DataFrame([supervisor_row]).to_csv(
        unit_dir / "supervision.tsv", sep="\t", index=False
    )


def _write_recovered_then_accepted_ledgers(bundle, unit, *, startup: bool) -> Path:
    _write_valid_attempt_ledgers(bundle, unit)
    schedule = v3._unit_seed_schedule(bundle, unit["unit_id"])
    unit_dir = bundle.campaign_dir / "work" / unit["unit_id"]
    attempts = pd.read_csv(unit_dir / "attempts.tsv", sep="\t")
    supervision = pd.read_csv(unit_dir / "supervision.tsv", sep="\t")
    first_attempt = attempts.iloc[0].to_dict()
    first_supervision = supervision.iloc[0].to_dict()
    first_attempt.update(
        {
            "panel_seed_used": False,
            "accepted": False,
            "rejection_reason": "worker_interrupted_before_base_attempt_commit",
            "status": "timed_out",
            "elapsed_seconds": 1.0,
            "supervisor_status": (
                "startup_interrupted_before_child_observed"
                if startup
                else "engine_complete_without_base_attempt_commit"
            ),
        }
    )
    first_supervision["supervisor_status"] = first_attempt["supervisor_status"]
    reason = (
        "worker_interrupted_active_engine_call"
        if startup
        else "engine_complete_without_base_attempt_commit"
    )
    if startup:
        first_attempt.update(
            {
                "child_pid": None,
                "child_pgid": None,
                "child_sid": None,
                "child_process_identity": None,
                "child_identity_observed": False,
                "engine_wrapper_pid": None,
                "engine_wrapper_process_identity": None,
                "parent_death_watchdog_armed": False,
                "supervisor_sigterm_sent": False,
                "supervisor_sigkill_sent": False,
                "direct_child_reaped": False,
                "attempt_temp_tree_created": False,
                "attempt_temp_tree_loaded": False,
            }
        )
        first_supervision.update(
            {
                "child_pid": None,
                "child_pgid": None,
                "child_sid": None,
                "child_process_identity": None,
                "child_identity_observed": False,
                "engine_wrapper_pid": None,
                "engine_wrapper_process_identity": None,
                "parent_death_watchdog_armed": False,
                "sigterm_sent": False,
                "sigkill_sent": False,
                "direct_child_reaped": False,
                "temp_tree_created": False,
                "temp_tree_loaded": False,
            }
        )
    second_attempt = attempts.iloc[0].to_dict()
    second_supervision = supervision.iloc[0].to_dict()
    second_schedule = schedule.iloc[1]
    second_attempt.update(
        {
            "attempt_zero_based": 1,
            "seed": int(second_schedule["simulation_seed"]),
            "planned_simulation_seed": int(second_schedule["simulation_seed"]),
            "panel_seed": int(second_schedule["panel_seed"]),
            "supervisor_execution_id": "b" * 64,
            "supervisor_pid": 200,
            "supervisor_process_identity": (
                "{'pid': 200, 'boot_id': 'test', 'start_ticks': 1}"
            ),
            "child_pid": 201,
            "child_pgid": 201,
            "child_sid": 201,
            "child_process_identity": (
                "{'pid': 201, 'boot_id': 'test', 'start_ticks': 1}"
            ),
            "engine_wrapper_pid": 202,
            "engine_wrapper_process_identity": (
                "{'pid': 202, 'boot_id': 'test', 'start_ticks': 1}"
            ),
        }
    )
    second_supervision.update(
        {
            "execution_id": "b" * 64,
            "attempt_zero_based": 1,
            "simulation_seed": int(second_schedule["simulation_seed"]),
            "supervisor_pid": 200,
            "supervisor_process_identity": (
                "{'pid': 200, 'boot_id': 'test', 'start_ticks': 1}"
            ),
            "child_pid": 201,
            "child_pgid": 201,
            "child_sid": 201,
            "child_process_identity": (
                "{'pid': 201, 'boot_id': 'test', 'start_ticks': 1}"
            ),
            "engine_wrapper_pid": 202,
            "engine_wrapper_process_identity": (
                "{'pid': 202, 'boot_id': 'test', 'start_ticks': 1}"
            ),
        }
    )
    pd.DataFrame([first_attempt, second_attempt]).to_csv(
        unit_dir / "attempts.tsv", sep="\t", index=False
    )
    pd.DataFrame([first_supervision, second_supervision]).to_csv(
        unit_dir / "supervision.tsv", sep="\t", index=False
    )
    reconciliation_path = unit_dir / "supervision_reconciliation.tsv"
    pd.DataFrame(
        [
            {
                "schema": v3.SUPERVISOR_RECONCILIATION_SCHEMA,
                "reason": reason,
                "execution_id": "a" * 64,
                "attempt_zero_based": 0,
                "simulation_seed": int(schedule.iloc[0]["simulation_seed"]),
                "prior_supervisor_pid": 100,
                "prior_child_pid": None if startup else 101,
                "prior_child_pgid": None if startup else 101,
                "process_group_gone": True,
                "temp_directory_removed": True,
                "charged_elapsed_seconds": 0.0,
                "transferred_to_attempt_elapsed_seconds": 1.0,
                "atomic_temp_files_removed": "[]",
                "atomic_temp_file_count_removed": 0,
                "reconciled_unix_ns": 1,
            }
        ]
    ).to_csv(reconciliation_path, sep="\t", index=False)
    return reconciliation_path


def test_attempt_audit_requires_panel_seed_used_equals_accepted(planned_bundle):
    _, _, _, _, bundle = planned_bundle
    unit = bundle.units.iloc[0].to_dict()
    _write_valid_attempt_ledgers(bundle, unit)
    result = v3._validate_attempt_ledger(bundle, unit)
    assert result["panel_seed_used_equals_accepted_for_every_attempt"] is True
    path = bundle.campaign_dir / "work" / unit["unit_id"] / "attempts.tsv"
    ledger = pd.read_csv(path, sep="\t")
    ledger["panel_seed_used"] = False
    ledger.to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match="panel-seed-used"):
        v3._validate_attempt_ledger(bundle, unit)


def test_uncommitted_supervision_is_reconciled_and_charged(tmp_path: Path):
    unit_dir = tmp_path / "unit"
    unit_dir.mkdir()
    pd.DataFrame(
        [
            {
                "schema": v3.SUPERVISOR_LEDGER_SCHEMA,
                "attempt_zero_based": 0,
                "simulation_seed": 12345,
                "supervisor_pid": 100,
                "child_pid": 101,
                "child_pgid": 101,
                "process_group_gone": True,
                "temp_tree_removed": True,
                "supervisor_elapsed_seconds": 2.5,
            }
        ]
    ).to_csv(unit_dir / "supervision.tsv", sep="\t", index=False)
    proxy = v3._SupervisedSlimEngine(
        object(),
        unit_dir=unit_dir,
        schedule=pd.DataFrame(
            [
                {
                    "attempt_zero_based": 0,
                    "simulation_seed": 12345,
                    "panel_seed": 54321,
                }
            ]
        ),
    )
    assert list(proxy.records) == [0]
    assert proxy.interrupted_elapsed_seconds == pytest.approx(0.0)
    attempts = pd.read_csv(unit_dir / "attempts.tsv", sep="\t")
    assert float(attempts.iloc[0]["elapsed_seconds"]) == pytest.approx(2.5)
    assert not bool(attempts.iloc[0]["accepted"])
    reconciliation = pd.read_csv(unit_dir / "supervision_reconciliation.tsv", sep="\t")
    assert reconciliation.iloc[0]["reason"] == (
        "engine_complete_without_base_attempt_commit"
    )
    assert float(reconciliation.iloc[0]["charged_elapsed_seconds"]) == 0.0
    assert float(reconciliation.iloc[0]["transferred_to_attempt_elapsed_seconds"]) == (
        pytest.approx(2.5)
    )


def test_active_only_startup_interruption_becomes_one_durable_disposition(
    tmp_path: Path,
):
    unit_dir = tmp_path / "unit"
    unit_dir.mkdir()
    scratch = unit_dir / ".han_v3_attempt_000_seed_12345_startup"
    scratch.mkdir()
    record = {
        "schema": v3.SUPERVISOR_LEDGER_SCHEMA,
        "attempt_zero_based": 0,
        "simulation_seed": 12345,
        "supervisor_started_unix_ns": time.time_ns() - 10_000_000,
        "supervisor_pid": 999_999_999,
        "supervisor_process_identity": {
            "pid": 999_999_999,
            "boot_id": "absent",
            "start_ticks": 1,
        },
        "timeout_seconds": 300.0,
        "term_grace_seconds": 0.2,
        "temp_directory": scratch.name,
        "temp_tree_created": False,
        "temp_tree_loaded": False,
        "temp_tree_removed": False,
    }
    v3._write_active_supervision(
        unit_dir / ".han_v3_active_000.json",
        {
            "schema": v3.SUPERVISOR_ACTIVE_SCHEMA,
            "status": "active",
            "attempt_record": record,
        },
    )
    schedule = pd.DataFrame(
        [{"attempt_zero_based": 0, "simulation_seed": 12345, "panel_seed": 54321}]
    )
    proxy = v3._SupervisedSlimEngine(
        object(), unit_dir=unit_dir, schedule=schedule, term_grace_seconds=0.2
    )
    assert proxy.records[0]["supervisor_status"] == (
        "startup_interrupted_before_child_observed"
    )
    assert proxy.records[0]["child_identity_observed"] is False
    assert not list(unit_dir.glob(".han_v3_active_*.json"))
    assert not list(unit_dir.glob(".han_v3_attempt_*"))
    attempts = pd.read_csv(unit_dir / "attempts.tsv", sep="\t")
    assert len(attempts) == 1
    assert int(attempts.iloc[0]["seed"]) == 12345
    reconciliation = pd.read_csv(unit_dir / "supervision_reconciliation.tsv", sep="\t")
    assert len(reconciliation) == 1


def test_pre_marker_orphan_cleanup_is_journaled_and_replayable(
    tmp_path: Path, monkeypatch
):
    unit_dir = tmp_path / "unit"
    unit_dir.mkdir()
    scratch = unit_dir / ".han_v3_attempt_000_seed_12345_orphaned"
    scratch.mkdir()
    (scratch / "partial.txt").write_text("partial", encoding="utf-8")
    active_temp = unit_dir / ".han_v3_active_000.json.tmp.999999"
    active_temp.write_text("partial", encoding="utf-8")
    schedule = pd.DataFrame(
        [{"attempt_zero_based": 0, "simulation_seed": 12345, "panel_seed": 54321}]
    )
    original_rmtree = v3.shutil.rmtree
    injected = {"raised": False}

    def fail_once(path, *args, **kwargs):
        if Path(path) == scratch and not injected["raised"]:
            injected["raised"] = True
            raise RuntimeError("injected cleanup crash")
        return original_rmtree(path, *args, **kwargs)

    monkeypatch.setattr(v3.shutil, "rmtree", fail_once)
    with pytest.raises(RuntimeError, match="injected cleanup crash"):
        v3._SupervisedSlimEngine(object(), unit_dir=unit_dir, schedule=schedule)
    cleanup_path = unit_dir / "supervision_startup_cleanup.tsv"
    pending = pd.read_csv(cleanup_path, sep="\t")
    assert "pending" in set(pending["status"])
    monkeypatch.setattr(v3.shutil, "rmtree", original_rmtree)
    v3._SupervisedSlimEngine(object(), unit_dir=unit_dir, schedule=schedule)
    assert not scratch.exists()
    assert not active_temp.exists()
    cleanup = pd.read_csv(cleanup_path, sep="\t")
    assert len(cleanup) == 2
    assert set(cleanup["status"]) == {"removed"}
    assert set(cleanup["artifact_kind"]) == {
        "attempt_scratch_directory",
        "active_atomic_temp",
    }


def test_reconciliation_journal_replays_after_crash_before_ledger_commit(
    tmp_path: Path, monkeypatch
):
    unit_dir = tmp_path / "unit"
    unit_dir.mkdir()
    scratch = unit_dir / ".han_v3_attempt_000_seed_12345_replay"
    scratch.mkdir()
    record = {
        "schema": v3.SUPERVISOR_LEDGER_SCHEMA,
        "execution_id": "a" * 64,
        "attempt_zero_based": 0,
        "simulation_seed": 12345,
        "supervisor_started_unix_ns": time.time_ns() - 10_000_000,
        "supervisor_pid": 999_999_999,
        "supervisor_process_identity": {
            "pid": 999_999_999,
            "boot_id": "absent",
            "start_ticks": 1,
        },
        "timeout_seconds": 300.0,
        "term_grace_seconds": 0.2,
        "temp_directory": scratch.name,
        "temp_tree_created": False,
        "temp_tree_loaded": False,
        "temp_tree_removed": False,
    }
    v3._write_active_supervision(
        unit_dir / ".han_v3_active_000.json",
        {
            "schema": v3.SUPERVISOR_ACTIVE_SCHEMA,
            "status": "active",
            "attempt_record": record,
        },
    )
    schedule = pd.DataFrame(
        [{"attempt_zero_based": 0, "simulation_seed": 12345, "panel_seed": 54321}]
    )
    original_persist = v3._SupervisedSlimEngine._persist
    injected = {"raised": False}

    def fail_once(self):
        if not injected["raised"]:
            injected["raised"] = True
            raise RuntimeError("injected post-journal crash")
        return original_persist(self)

    monkeypatch.setattr(v3._SupervisedSlimEngine, "_persist", fail_once)
    with pytest.raises(RuntimeError, match="injected post-journal crash"):
        v3._SupervisedSlimEngine(object(), unit_dir=unit_dir, schedule=schedule)
    reconciliation_path = unit_dir / "supervision_reconciliation.tsv"
    durable_before = reconciliation_path.read_bytes()
    durable = pd.read_csv(reconciliation_path, sep="\t").iloc[0]
    monkeypatch.setattr(v3._SupervisedSlimEngine, "_persist", original_persist)
    v3._SupervisedSlimEngine(object(), unit_dir=unit_dir, schedule=schedule)
    assert reconciliation_path.read_bytes() == durable_before
    attempts = pd.read_csv(unit_dir / "attempts.tsv", sep="\t")
    assert float(attempts.iloc[0]["elapsed_seconds"]) == pytest.approx(
        float(durable["transferred_to_attempt_elapsed_seconds"])
    )
    assert not list(unit_dir.glob(".han_v3_active_*.json"))


def test_reconciliation_replay_normalizes_nullable_integer_columns(tmp_path: Path):
    unit_dir = tmp_path / "unit"
    unit_dir.mkdir()
    schedule = pd.DataFrame(
        [
            {"attempt_zero_based": 0, "simulation_seed": 12345, "panel_seed": 54321},
            {"attempt_zero_based": 1, "simulation_seed": 12346, "panel_seed": 54322},
        ]
    )
    proxy = v3._SupervisedSlimEngine(object(), unit_dir=unit_dir, schedule=schedule)
    base = {
        "schema": v3.SUPERVISOR_RECONCILIATION_SCHEMA,
        "reason": "worker_interrupted_active_engine_call",
        "active_record_path": ".han_v3_active_000.json",
        "active_record_sha256": "c" * 64,
        "prior_supervisor_pid": 100,
        "process_group_gone": True,
        "temp_directory_removed": True,
        "charged_elapsed_seconds": 0.0,
        "transferred_to_attempt_elapsed_seconds": 1.0,
        "atomic_temp_files_removed": "[]",
        "atomic_temp_file_count_removed": 0,
        "reconciled_unix_ns": 1,
    }
    proxy._append_reconciliation(
        {
            **base,
            "execution_id": "a" * 64,
            "attempt_zero_based": 0,
            "simulation_seed": 12345,
            "prior_child_pid": None,
            "prior_child_pgid": None,
        }
    )
    observed = {
        **base,
        "execution_id": "b" * 64,
        "attempt_zero_based": 1,
        "simulation_seed": 12346,
        "active_record_path": ".han_v3_active_001.json",
        "prior_child_pid": 101,
        "prior_child_pgid": 101,
        "reconciled_unix_ns": 2,
    }
    proxy._append_reconciliation(observed)
    replayed = proxy._append_reconciliation(observed)
    assert int(float(replayed["prior_child_pid"])) == 101
    assert len(pd.read_csv(unit_dir / "supervision_reconciliation.tsv", sep="\t")) == 2


def test_attempt_audit_validates_reconciliation_disposition_and_budget(
    planned_bundle,
):
    _, _, _, _, bundle = planned_bundle
    unit = bundle.units.iloc[0].to_dict()
    reconciliation_path = _write_recovered_then_accepted_ledgers(
        bundle, unit, startup=False
    )
    result = v3._validate_attempt_ledger(bundle, unit)
    assert result["interrupted_elapsed_seconds_charged"] == pytest.approx(0.0)
    assert result["one_terminal_disposition_per_attempted_seed"] is True
    tampered = pd.read_csv(reconciliation_path, sep="\t")
    tampered.loc[0, "execution_id"] = "c" * 64
    tampered.to_csv(reconciliation_path, sep="\t", index=False)
    with pytest.raises(ValueError, match="reconciliation disposition"):
        v3._validate_attempt_ledger(bundle, unit)


def test_startup_null_identity_recovery_passes_eventual_full_attempt_audit(
    planned_bundle,
):
    _, _, _, _, bundle = planned_bundle
    unit = bundle.units.iloc[0].to_dict()
    _write_recovered_then_accepted_ledgers(bundle, unit, startup=True)
    result = v3._validate_attempt_ledger(bundle, unit)
    assert result["attempts_completed"] == 2
    assert result["all_process_groups_gone_and_children_reaped"] is True


def test_completion_committed_stale_marker_does_not_rewrite_ledgers(
    planned_bundle,
):
    _, _, _, _, bundle = planned_bundle
    unit = bundle.units.iloc[0].to_dict()
    _write_valid_attempt_ledgers(bundle, unit)
    unit_dir = bundle.campaign_dir / "work" / unit["unit_id"]
    attempts_path = unit_dir / "attempts.tsv"
    supervision_path = unit_dir / "supervision.tsv"
    attempts_before = attempts_path.read_bytes()
    supervision_before = supervision_path.read_bytes()
    record = pd.read_csv(supervision_path, sep="\t").iloc[0].to_dict()
    record.update(
        {
            "lifecycle_phase": ("accepted_attempt_committed_pending_unit_completion"),
            "base_attempt_committed_unix_ns": time.time_ns(),
            "base_attempt_ledger_sha256": v3.sha256_file(attempts_path),
        }
    )
    v3._write_active_supervision(
        unit_dir / ".han_v3_active_000.json",
        {
            "schema": v3.SUPERVISOR_ACTIVE_SCHEMA,
            "status": "active",
            "attempt_record": record,
        },
    )
    (unit_dir / "simulation_complete.json").write_text("{}", encoding="utf-8")
    v3._SupervisedSlimEngine(
        object(),
        unit_dir=unit_dir,
        schedule=v3._unit_seed_schedule(bundle, unit["unit_id"]),
    )
    assert not list(unit_dir.glob(".han_v3_active_*.json"))
    assert attempts_path.read_bytes() == attempts_before
    assert supervision_path.read_bytes() == supervision_before
    assert not (unit_dir / "supervision_reconciliation.tsv").exists()


def test_persisted_engine_error_is_fail_closed_and_never_retried(tmp_path: Path):
    unit_dir = tmp_path / "unit"
    unit_dir.mkdir()
    scratch = unit_dir / ".han_v3_attempt_000_seed_12345_systemic"
    scratch.mkdir()
    record = {
        "schema": v3.SUPERVISOR_LEDGER_SCHEMA,
        "execution_id": "a" * 64,
        "attempt_zero_based": 0,
        "simulation_seed": 12345,
        "supervisor_started_unix_ns": time.time_ns() - 10_000_000,
        "supervisor_pid": 999_999_999,
        "supervisor_process_identity": {
            "pid": 999_999_999,
            "boot_id": "absent",
            "start_ticks": 1,
        },
        "supervisor_status": "child_error",
        "child_error_type": "ValueError",
        "child_error": "systemic engine failure",
        "timeout_seconds": 300.0,
        "term_grace_seconds": 0.2,
        "temp_directory": scratch.name,
        "temp_tree_created": False,
        "temp_tree_loaded": False,
        "temp_tree_removed": False,
    }
    supervision_path = unit_dir / "supervision.tsv"
    pd.DataFrame([record]).to_csv(supervision_path, sep="\t", index=False)
    v3._write_active_supervision(
        unit_dir / ".han_v3_active_000.json",
        {
            "schema": v3.SUPERVISOR_ACTIVE_SCHEMA,
            "status": "active",
            "attempt_record": record,
        },
    )
    supervision_before = supervision_path.read_bytes()
    schedule = pd.DataFrame(
        [{"attempt_zero_based": 0, "simulation_seed": 12345, "panel_seed": 54321}]
    )
    with pytest.raises(
        RuntimeError,
        match=("persisted nonretryable supervisor failure: child_error.*ValueError"),
    ):
        v3._SupervisedSlimEngine(object(), unit_dir=unit_dir, schedule=schedule)
    assert supervision_path.read_bytes() == supervision_before
    assert (unit_dir / ".han_v3_active_000.json").is_file()
    assert scratch.is_dir()
    assert not (unit_dir / "attempts.tsv").exists()
    assert not (unit_dir / "supervision_reconciliation.tsv").exists()


def test_read_only_cache_validation_does_not_rebuild_han_calibration(
    planned_bundle, monkeypatch
):
    _, _, _, _, bundle = planned_bundle
    unit = bundle.units.iloc[0].to_dict()
    record = simulation._validate_unit(unit)
    binding = v3._unit_binding(bundle, unit)
    monkeypatch.setattr(v3, "_v3_implementation_paths", lambda _bundle: ())
    monkeypatch.setattr(
        simulation,
        "_selected_objects",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("cache validation must not rebuild Han calibration")
        ),
    )
    monkeypatch.setattr(
        simulation,
        "_validate_completion",
        lambda _artifacts, _contract: SimpleNamespace(cache_hit=True),
    )
    monkeypatch.setattr(
        v3,
        "_validate_attempt_ledger",
        lambda _bundle, _unit: {
            "attempts_completed": 1,
            "remaining_cumulative_timeout_seconds": 90_000.0,
            "accepted_panel_seed": int(
                v3._unit_seed_schedule(bundle, unit["unit_id"]).iloc[0]["panel_seed"]
            ),
        },
    )
    contract = {
        "schema": simulation.SCHEMA_VERSION,
        "unit": record,
        "parameters": {
            "slim_scaling_factor": 5.0,
            "slim_burn_in": 0.1,
            "candidate_pool_diploids": 500,
            "han_selection_calibration": binding,
        },
        "implementation": {
            "sources": {},
            "eas_resource_sha256": simulation.EAS_RESOURCE_SHA256,
            "slim": {
                "path": simulation._portable_path(bundle.slim_path, bundle.repo_root),
                "sha256": v3.sha256_file(bundle.slim_path),
            },
        },
    }
    completion = {
        "schema": simulation.SCHEMA_VERSION,
        "status": "complete",
        "contract": contract,
        "contract_sha256": v3._canonical_sha256(contract),
        "provenance": {
            "han_direct_production_v3": {
                **binding,
                "direct_operational_estimand": True,
                "normal_10mb_stdpopsim_processing": True,
                "strict_selected_type1_identity": True,
                "fixed_cessation_endpoint": float(
                    v3._unit_seed_schedule(bundle, unit["unit_id"]).iloc[0][
                        "realized_selection_end_generations_ago"
                    ]
                ),
                "conditional_inference": True,
            },
            "selected_focal_identity": {
                "status": "valid",
                "observed_mutation_type": 1,
            },
            "candidate_pool_frequency_gate": {"passed": True},
            "sample_panel": {
                "panel_seed": int(
                    v3._unit_seed_schedule(bundle, unit["unit_id"]).iloc[0][
                        "panel_seed"
                    ]
                )
            },
        },
        "search_budget": {
            "selected_max_draws": 300,
            "selected_draw_timeout_seconds": 300.0,
            "selected_cumulative_timeout_seconds": 90_000.0,
        },
        "attempts_completed": 1,
    }
    unit_dir = bundle.campaign_dir / "work" / unit["unit_id"]
    endpoint = float(
        v3._unit_seed_schedule(bundle, unit["unit_id"]).iloc[0][
            "realized_selection_end_generations_ago"
        ]
    )
    endpoint_manifest = next(
        item
        for item in bundle.manifest["model"]["fixed_endpoints"]
        if item["selection_coefficient"] == unit["selection_coefficient"]
        and item["target_allele_frequency"] == unit["target_allele_frequency"]
    )
    boundary = 2015.0
    intervals = (
        [{"population": "Loschbour", "start_time": 2272.0, "end_time": endpoint}]
        if endpoint >= boundary
        else [
            {"population": "Loschbour", "start_time": 2272.0, "end_time": boundary},
            {"population": "Han", "start_time": boundary, "end_time": endpoint},
        ]
    )
    serialized = [
        {
            "event_type": "DrawMutation",
            "population": "Neanderthal",
            "single_site_id": v3.fixed_v2.FOCAL_SITE_ID,
            "time": {"generations_ago": 2400.0},
        },
        *[
            {
                "event_type": "ConditionOnAlleleFrequency",
                "population": population,
                "operator": operator,
                "allele_frequency": frequency,
                "single_site_id": v3.fixed_v2.FOCAL_SITE_ID,
                "start_time": {"generations_ago": start},
                "end_time": {"generations_ago": end},
            }
            for population, operator, frequency, start, end in (
                ("Neanderthal", ">", 0.0, 2400.0, 2273.0),
                ("Neanderthal", ">=", 1e-9, 2273.0, 2273.0),
                ("Loschbour", ">", 0.0, 2272.0, 2015.0),
                ("Han", ">", 0.0, 2015.0, 0.0),
            )
        ],
        *[
            {
                "event_type": "ChangeMutationFitness",
                "population": interval["population"],
                "selection_coefficient": unit["selection_coefficient"],
                "dominance_coefficient": 0.5,
                "single_site_id": v3.fixed_v2.FOCAL_SITE_ID,
                "start_time": {"generations_ago": interval["start_time"]},
                "end_time": {"generations_ago": interval["end_time"]},
            }
            for interval in intervals
        ],
    ]
    realized = {
        "mutation_age_generations": 2400.0,
        "source_check_generations": 2275.0,
        "introgression_and_selection_onset_generations": 2270.0,
        "demographic_han_split_generations": 2015.0,
        "selection_end_generations": endpoint,
    }
    completion["provenance"]["han_direct_production_v3"]["event_audit"] = {
        "schema": v3.fixed_v2.EVENT_AUDIT_SCHEMA,
        "mutation_population": "Neanderthal",
        "mutation_age_generations": 2400.0,
        "selection_coefficient": unit["selection_coefficient"],
        "duration_multiplier": endpoint_manifest["duration_multiplier"],
        "requested_selection_end_generations_ago": endpoint_manifest[
            "requested_selection_end_generations_ago"
        ],
        "realized_selection_end_generations_ago": endpoint,
        "selection_intervals": intervals,
        "selection_interval_count": len(intervals),
        "recipient_nonloss_condition_count": 2,
        "source_condition_count": 2,
        "terminal_han_condition_count": 0,
        "external_terminal_ascertainment": True,
        "realized_q5_event_times_generations": realized,
        "realized_event_schedule_sha256": v3._canonical_sha256(realized),
        "focal_event_count": len(serialized),
        "serialized_extended_events": serialized,
        "serialized_extended_events_sha256": v3._canonical_sha256(serialized),
    }
    unit_dir.mkdir(parents=True, exist_ok=True)
    (unit_dir / "simulation_complete.json").write_text(
        json.dumps(completion), encoding="utf-8"
    )
    result = v3.validate_v3_cached_unit(bundle, unit)
    assert result["cache_validated_read_only"] is True
    assert result["authorization_validated"] is True


def test_additive_dispatch_routes_han_cache_and_delegates_non_han(
    planned_bundle, monkeypatch
):
    _, _, _, _, bundle = planned_bundle
    han = bundle.units.iloc[0].to_dict()
    eas = dict(han)
    eas.update(
        {
            "unit_id": "eas_phlash_median__af10__selected_s0p005__rep001",
            "demography_id": "eas_phlash_median",
            "demography_kind": "no_introgression",
            "population": "EAS",
            "source_model": "PHLASH EAS pointwise median trajectory",
            "selection_origin": "de_novo",
        }
    )
    delegated = object()
    delegate_calls = []

    def delegate(unit, *_args, **_kwargs):
        delegate_calls.append(unit["unit_id"])
        return delegated

    monkeypatch.setattr(simulation, "_simulate_unit_unlocked", delegate)
    monkeypatch.setattr(
        v3,
        "validate_v3_cached_unit",
        lambda _bundle, _unit: {"cache_validated_read_only": True},
    )
    with v3.patched_v3_cache_validation_dispatch(bundle):
        han_result = simulation._simulate_unit_unlocked(
            han, bundle.repo_root, bundle.campaign_dir
        )
        eas_result = simulation._simulate_unit_unlocked(
            eas, bundle.repo_root, bundle.campaign_dir
        )
    assert han_result.cache_hit is True
    assert eas_result is delegated
    assert delegate_calls == [eas["unit_id"]]
    assert simulation._simulate_unit_unlocked is delegate


def test_standard_aggregate_route_executes_under_v3_dispatch(
    planned_bundle, monkeypatch
):
    _, _, _, _, bundle = planned_bundle
    han = bundle.units.iloc[[0]].copy()
    monkeypatch.setattr(
        v3,
        "validate_v3_cached_unit",
        lambda _bundle, _unit: {"cache_validated_read_only": True},
    )
    outputs = tuple(bundle.repo_root / f"result-{index}.tsv" for index in range(5))

    def fake_standard_aggregate(campaign_dir, units):
        assert campaign_dir == bundle.campaign_dir
        artifacts = simulation._simulate_unit_unlocked(
            units.iloc[0].to_dict(), bundle.repo_root, bundle.campaign_dir
        )
        assert artifacts.cache_hit is True
        return outputs

    monkeypatch.setattr(
        focused_campaign, "_aggregate_compact_decode", fake_standard_aggregate
    )
    observed = v3.aggregate_decoded_outputs_with_v3_dispatch(bundle, han)
    assert observed == outputs


def test_v3_cache_identity_is_stable_inside_eas_selected_context(
    planned_bundle, monkeypatch
):
    _, _, _, _, bundle = planned_bundle
    expected = v3._v3_implementation_paths(bundle)
    mutated = ("temporary/eas/integration/path.py",)

    @v3.contextmanager
    def fake_selected_context(_selected_bundle):
        original = simulation.IMPLEMENTATION_PATHS
        simulation.IMPLEMENTATION_PATHS = mutated
        try:
            assert v3._v3_implementation_paths(bundle) == expected
            yield
        finally:
            simulation.IMPLEMENTATION_PATHS = original

    from gamma_smc_aou import focused_selected_production as selected_production

    monkeypatch.setattr(
        selected_production,
        "load_selected_production_bundle",
        lambda *_args, **_kwargs: object(),
    )
    monkeypatch.setattr(
        selected_production, "patched_selected_executor", fake_selected_context
    )
    outputs = tuple(bundle.repo_root / f"nested-{index}.tsv" for index in range(5))
    delegated = []

    def delegate(unit, *_args, **_kwargs):
        delegated.append((unit["simulation_class"], simulation.IMPLEMENTATION_PATHS))
        return SimpleNamespace(cache_hit=True)

    monkeypatch.setattr(simulation, "_simulate_unit_unlocked", delegate)
    monkeypatch.setattr(
        v3,
        "validate_v3_cached_unit",
        lambda _bundle, _unit: {"cache_validated_read_only": True},
    )

    def fake_aggregate(_campaign, units):
        assert (
            simulation.IMPLEMENTATION_PATHS == v3._BASE_SIMULATION_IMPLEMENTATION_PATHS
        )
        assert v3._v3_implementation_paths(bundle) == expected
        for unit in units.to_dict(orient="records"):
            simulation._simulate_unit_unlocked(
                unit, bundle.repo_root, bundle.campaign_dir
            )
        return outputs

    monkeypatch.setattr(focused_campaign, "_aggregate_compact_decode", fake_aggregate)
    observed = v3.aggregate_decoded_outputs_with_v3_dispatch(
        bundle,
        pd.DataFrame(
            [
                bundle.units.iloc[0].to_dict(),
                {
                    **bundle.units.iloc[0].to_dict(),
                    "unit_id": "eas_phlash_median__af10__selected_s0p005__rep001",
                    "demography_id": "eas_phlash_median",
                    "demography_kind": "no_introgression",
                    "population": "EAS",
                    "source_model": "PHLASH EAS pointwise median trajectory",
                    "selection_origin": "de_novo",
                },
                {
                    **bundle.units.iloc[0].to_dict(),
                    "unit_id": "eas_phlash_median__af10__neutral__rep001",
                    "demography_id": "eas_phlash_median",
                    "demography_kind": "no_introgression",
                    "population": "EAS",
                    "source_model": "PHLASH EAS pointwise median trajectory",
                    "selection_origin": "neutral",
                    "simulation_class": "neutral",
                    "selection_coefficient": 0.0,
                },
            ]
        ),
        selected_integration_manifest=bundle.manifest_path,
    )
    assert observed == outputs
    assert delegated == [
        ("selected", mutated),
        ("neutral", v3._BASE_SIMULATION_IMPLEMENTATION_PATHS),
    ]
    assert simulation.IMPLEMENTATION_PATHS == v3._BASE_SIMULATION_IMPLEMENTATION_PATHS


def test_pilot_gate_revalidates_current_completion_rows(planned_bundle, monkeypatch):
    _, _, _, output, bundle = planned_bundle
    rows = [{"unit_id": f"pilot-{index}"} for index in range(6)]
    payload = {
        "schema": v3.STAGE_AUDIT_SCHEMA,
        "stage": v3.STAGE_PILOT,
        "status": "passed",
        "manifest_sha256": v3.sha256_file(bundle.manifest_path),
        "seed_plan_sha256": v3.sha256_file(bundle.seed_plan_path),
        "accepted_units": 6,
        "all_current_completions_cache_validated_read_only": True,
        "all_current_completion_budgets_match_authorization": True,
        "seed_plan_audit": v3._seed_plan_audit(bundle),
        "rows": rows,
    }
    payload["payload_sha256"] = v3._canonical_sha256(payload)
    (output / "pilot_audit.json").write_text(json.dumps(payload), encoding="utf-8")
    monkeypatch.setattr(
        v3,
        "aggregate_v3_cached_completions",
        lambda *_args, **_kwargs: pd.DataFrame(
            [{"unit_id": f"changed-{index}"} for index in range(6)]
        ),
    )
    with pytest.raises(ValueError, match="current validated completions"):
        v3._required_pilot_audit(bundle)


@pytest.mark.skipif(os.name != "posix", reason="requires Linux/WSL process groups")
def test_supervisor_reaps_real_child_and_grandchild_on_timeout(tmp_path: Path):
    unit_dir = tmp_path / "unit"
    unit_dir.mkdir()
    schedule = pd.DataFrame(
        [{"attempt_zero_based": 0, "simulation_seed": 12345, "panel_seed": 54321}]
    )
    proxy = v3._SupervisedSlimEngine(
        _HangingEngineWithGrandchild(),
        unit_dir=unit_dir,
        schedule=schedule,
        term_grace_seconds=0.2,
    )
    started = time.monotonic()
    with (
        pytest.raises(simulation.SimulationDrawTimeout),
        v3._supervised_draw_timeout(0.4),
    ):
        proxy.simulate(
            seed=12345,
            logfile=str(unit_dir / "slim_draw_000.csv"),
        )
    assert time.monotonic() - started < 5.0
    record = proxy.records[0]
    assert record["supervisor_status"] == "timed_out"
    assert record["child_pid"] == record["child_pgid"] == record["child_sid"]
    assert record["parent_death_watchdog_armed"] is True
    assert record["direct_child_reaped"] is True
    assert record["process_group_gone"] is True
    assert record["sigterm_sent"] is True
    assert record["sigkill_sent"] is True
    assert record["temp_tree_removed"] is True
    grandchild_pid = int(
        (unit_dir / "slim_draw_000.csv.grandchild.pid").read_text(encoding="utf-8")
    )
    with pytest.raises(ProcessLookupError):
        os.kill(grandchild_pid, 0)
    assert not list(unit_dir.glob(".han_v3_attempt_*"))


@pytest.mark.skipif(os.name != "posix", reason="requires Linux/WSL process groups")
def test_parent_death_watchdog_kills_group_and_restart_reconciles_temp(
    tmp_path: Path,
):
    unit_dir = tmp_path / "unit"
    unit_dir.mkdir()
    context = multiprocessing.get_context("fork")
    worker = context.Process(target=_run_hanging_proxy_worker, args=(str(unit_dir),))
    worker.start()
    wrapper_path = unit_dir / "slim_draw_000.csv.wrapper.pid"
    grandchild_path = unit_dir / "slim_draw_000.csv.grandchild.pid"
    deadline = time.monotonic() + 10.0
    while (
        not wrapper_path.is_file() or not grandchild_path.is_file()
    ) and time.monotonic() < deadline:
        time.sleep(0.02)
    assert wrapper_path.is_file() and grandchild_path.is_file()
    wrapper_pid = int(wrapper_path.read_text(encoding="utf-8"))
    grandchild_pid = int(grandchild_path.read_text(encoding="utf-8"))
    worker.kill()
    worker.join(timeout=5.0)
    assert not worker.is_alive()
    active = list(unit_dir.glob(".han_v3_active_*.json"))
    scratch = list(unit_dir.glob(".han_v3_attempt_*"))
    assert len(active) == len(scratch) == 1
    schedule = pd.DataFrame(
        [{"attempt_zero_based": 0, "simulation_seed": 12345, "panel_seed": 54321}]
    )
    reconciled_proxy = v3._SupervisedSlimEngine(
        _HangingEngineWithGrandchild(),
        unit_dir=unit_dir,
        schedule=schedule,
        term_grace_seconds=0.2,
    )
    living = [
        pid
        for pid in (wrapper_pid, grandchild_pid)
        if not v3._process_absent_or_zombie(pid)
    ]
    assert not living
    assert not list(unit_dir.glob(".han_v3_active_*.json"))
    assert not list(unit_dir.glob(".han_v3_attempt_*"))
    assert reconciled_proxy.interrupted_elapsed_seconds == 0
    interrupted_attempt = pd.read_csv(unit_dir / "attempts.tsv", sep="\t").iloc[0]
    assert float(interrupted_attempt["elapsed_seconds"]) > 0
    assert str(interrupted_attempt["accepted"]).casefold() == "false"
    reconciliation = pd.read_csv(unit_dir / "supervision_reconciliation.tsv", sep="\t")
    assert reconciliation.iloc[0]["schema"] == v3.SUPERVISOR_RECONCILIATION_SCHEMA
    assert bool(reconciliation.iloc[0]["process_group_gone"])
    assert bool(reconciliation.iloc[0]["temp_directory_removed"])
    assert float(reconciliation.iloc[0]["charged_elapsed_seconds"]) == 0
    assert float(reconciliation.iloc[0]["transferred_to_attempt_elapsed_seconds"]) > 0


def test_conditional_inference_caveats_are_explicit(planned_bundle):
    _, _, _, _, bundle = planned_bundle
    contract = bundle.manifest["inference_contract"]
    assert contract["endpoint_rule_selected_after_exploratory_screen"] is True
    assert contract["exploratory_and_failed_v2_runs_inferential_use"] is False
    assert (
        contract["neutral_rarity_or_unconditional_selection_probability_estimated"]
        is False
    )
    assert contract["pilot_is_final_contract_and_counts_toward_sixty"] is True
    assert len(bundle.manifest["scientific_caveats"]) >= 5
