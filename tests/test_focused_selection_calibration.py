from __future__ import annotations

import json
from dataclasses import dataclass

import numpy as np
import pandas as pd
import pytest
import stdpopsim
import tskit
from gamma_smc_aou import focused_selection_calibration as calibration
from gamma_smc_aou import focused_selection_campaign as campaign
from gamma_smc_aou import focused_selection_simulation as simulation
from gamma_smc_aou.eas_sweep_models import FOCAL_SITE_ID


def _plan(tmp_path, *, draws=20):
    design = calibration.HanSelectionEndCalibrationDesign(draws_per_candidate=draws)
    production = campaign.build_execution_units(
        campaign.FocusedCampaignPlan(
            repo_root=tmp_path,
            campaign_dir=tmp_path / "focused_selection_EAS_sim",
            selected_replicates=1,
            neutral_replicates=1,
        )
    )
    return calibration.build_calibration_plan(
        design, production_seeds=production["seed"]
    )


def test_default_plan_is_840_seed_disjoint_q5_draws(tmp_path):
    production = campaign.build_execution_units(
        campaign.FocusedCampaignPlan(
            repo_root=tmp_path,
            campaign_dir=tmp_path / "focused_selection_EAS_sim",
        )
    )
    plan = calibration.build_calibration_plan(
        calibration.HanSelectionEndCalibrationDesign(),
        production_seeds=production["seed"],
    )
    assert len(plan) == 840
    assert plan["calibration_id"].is_unique
    assert plan["seed"].is_unique
    assert set(plan["seed"]).isdisjoint(set(production["seed"]))
    assert set(plan["slim_scaling_factor"]) == {5.0}
    assert set(
        plan.groupby(["selection_coefficient", "target_allele_frequency"]).size()
    ) == {140}
    for _, group in plan.groupby(["selection_coefficient", "target_allele_frequency"]):
        assert group["realized_selection_end_generations_ago"].nunique() == 7


def _toy_sweep():
    events = (
        stdpopsim.ChangeMutationFitness(
            start_time=2272,
            end_time=2015,
            single_site_id=FOCAL_SITE_ID,
            population="Loschbour",
            selection_coeff=0.01,
            dominance_coeff=0.5,
        ),
        stdpopsim.ChangeMutationFitness(
            start_time=2015,
            end_time=1200,
            single_site_id=FOCAL_SITE_ID,
            population="Han",
            selection_coeff=0.01,
            dominance_coeff=0.5,
        ),
    )

    @dataclass(frozen=True)
    class Sweep:
        extended_events: tuple
        record: dict

        def provenance_record(self):
            return self.record

    return Sweep(events, {"selection": {}})


def test_piecewise_endpoint_shortens_loschbour_or_han():
    early, early_record = calibration.apply_han_selection_end(_toy_sweep(), 2150)
    fitness = [
        event
        for event in early.extended_events
        if isinstance(event, stdpopsim.ChangeMutationFitness)
    ]
    assert [
        (event.population, event.start_time, event.end_time) for event in fitness
    ] == [("Loschbour", 2272, 2150)]
    assert early_record["selection_end_generations_ago"] == 2150

    recent, _ = calibration.apply_han_selection_end(_toy_sweep(), 1500)
    fitness = [
        event
        for event in recent.extended_events
        if isinstance(event, stdpopsim.ChangeMutationFitness)
    ]
    assert [
        (event.population, event.start_time, event.end_time) for event in fitness
    ] == [
        ("Loschbour", 2272, 2015),
        ("Han", 2015, 1500),
    ]


def _synthetic_ledger(plan: pd.DataFrame, *, good=True) -> pd.DataFrame:
    rows = []
    for row in plan.to_dict(orient="records"):
        multiplier = float(row["duration_multiplier"])
        hit = np.isclose(multiplier, 0.25) or (
            np.isclose(multiplier, 0.40) and int(row["draw_index"]) == 0
        )
        af = (
            float(row["target_allele_frequency"])
            if hit
            else min(0.99, float(row["target_allele_frequency"]) + 0.30)
        )
        if not good and np.isclose(row["selection_coefficient"], 0.005):
            hit = False
            af = 0.8
        rows.append(
            {
                **row,
                "status": "complete",
                "evaluable": True,
                "terminal_class": "segregating",
                "candidate_pool_alt_count": round(1000 * af),
                "candidate_pool_af": af,
                "candidate_pool_frequency_gate_passed": hit,
                "exact_panel_feasible": hit,
                "elapsed_seconds": 1.0,
                "error": "",
            }
        )
    return pd.DataFrame(rows)


def test_summary_prioritizes_hit_rate_and_finalize_is_fail_closed(tmp_path):
    plan = _plan(tmp_path, draws=5)
    ledger = _synthetic_ledger(plan)
    summary = calibration.summarize_terminal_draws(plan, ledger)
    selected = calibration.select_calibration_endpoints(summary)
    assert set(selected["duration_multiplier"]) == {0.25}
    sensitivity = tmp_path / "calibration" / "sensitivity.tsv"
    sensitivity.parent.mkdir()
    sensitivity.write_text("cell\tpassed\nall\ttrue\n", encoding="utf-8")
    payload = calibration.finalize_calibration(
        plan,
        ledger,
        tmp_path / "calibration",
        production_seed_digest="abc",
        minimum_band_hit_rate=0.2,
        sensitivity_artifacts={sensitivity.name: simulation.sha256_file(sensitivity)},
        sensitivity_passed=True,
    )
    assert payload["status"] == "frozen"
    frozen = tmp_path / "calibration" / "han_selection_end_frozen.json"
    assert frozen.is_file()
    assert len(json.loads(frozen.read_text())["mapping"]) == 6
    assert (tmp_path / "calibration" / "han_selection_end_response.png").is_file()
    assert (tmp_path / "calibration" / "han_selection_end_response.pdf").is_file()

    failed = calibration.finalize_calibration(
        plan,
        _synthetic_ledger(plan, good=False),
        tmp_path / "calibration",
        production_seed_digest="abc",
        minimum_band_hit_rate=0.2,
        sensitivity_artifacts={sensitivity.name: simulation.sha256_file(sensitivity)},
        sensitivity_passed=True,
    )
    assert failed["status"] == "not_frozen"
    assert not frozen.exists()


def test_finalize_requires_ten_mb_sensitivity_before_freeze(tmp_path):
    plan = _plan(tmp_path, draws=5)
    payload = calibration.finalize_calibration(
        plan,
        _synthetic_ledger(plan),
        tmp_path / "calibration",
        production_seed_digest="abc",
        minimum_band_hit_rate=0.2,
    )
    assert payload["status"] == "not_frozen"
    contract = payload["production_contract"]
    assert contract["slim_scaling_factor"] == 5
    assert contract["candidate_pool_diploids"] == 500
    assert contract["af_half_width"] == pytest.approx(0.025)
    assert contract["post_pulse_recipient_nonloss_conditioning"] is True
    assert contract["screen_contig_production_use_permitted"] is False
    assert contract["ten_mb_sensitivity_passed"] is False


def test_summary_rejects_missing_or_duplicate_terminal_draws(tmp_path):
    plan = _plan(tmp_path, draws=1)
    ledger = _synthetic_ledger(plan)
    with pytest.raises(ValueError, match="cover"):
        calibration.summarize_terminal_draws(plan, ledger.iloc[:-1])
    with pytest.raises(ValueError, match="duplicate"):
        calibration.summarize_terminal_draws(
            plan, pd.concat([ledger, ledger.iloc[[0]]], ignore_index=True)
        )


def test_cached_draw_contract_fails_closed(tmp_path, monkeypatch):
    plan = _plan(tmp_path, draws=1).iloc[[0]]
    row = plan.iloc[0].to_dict()
    slim = tmp_path / "slim"
    slim.write_bytes(b"slim")
    source = tmp_path / "python/gamma_smc_aou"
    source.mkdir(parents=True)
    (source / "focused_selection_simulation.py").write_bytes(b"simulation")
    monkeypatch.setattr(calibration, "__file__", str(tmp_path / "calibration.py"))
    (tmp_path / "calibration.py").write_bytes(b"calibration")
    work = tmp_path / "work"
    directory = work / row["calibration_id"]
    directory.mkdir(parents=True)
    contract = calibration._draw_contract(row, tmp_path, slim)
    result = {**row, "status": "complete", "evaluable": True}
    (directory / "completion.json").write_text(
        json.dumps(
            {
                "schema": calibration.SCHEMA_VERSION,
                "contract": contract,
                "contract_sha256": simulation._canonical_sha256(contract),
                "result": result,
            }
        )
    )
    assert (
        calibration._run_draw(
            {
                "row": row,
                "repo_root": str(tmp_path),
                "work_root": str(work),
                "slim_path": str(slim),
                "timeout_seconds": 10,
            }
        )["status"]
        == "complete"
    )
    payload = json.loads((directory / "completion.json").read_text())
    payload["contract"]["slim_burn_in"] = 0.2
    (directory / "completion.json").write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="incompatible"):
        calibration._run_draw(
            {
                "row": row,
                "repo_root": str(tmp_path),
                "work_root": str(work),
                "slim_path": str(slim),
                "timeout_seconds": 10,
            }
        )


def _screen_with_identity_failures(plan: pd.DataFrame, count: int = 3) -> pd.DataFrame:
    ledger = _synthetic_ledger(plan)
    indexes = ledger.index[:count]
    ledger.loc[indexes, "status"] = "failed"
    ledger.loc[indexes, "evaluable"] = False
    ledger.loc[indexes, "terminal_class"] = "failed"
    ledger.loc[indexes, "candidate_pool_alt_count"] = np.nan
    ledger.loc[indexes, "candidate_pool_af"] = np.nan
    ledger.loc[indexes, "candidate_pool_frequency_gate_passed"] = False
    ledger.loc[indexes, "exact_panel_feasible"] = False
    ledger.loc[indexes, "error"] = calibration.RECONCILABLE_IDENTITY_ERROR
    return ledger


def _reconciliation_binding() -> dict:
    return {
        "original_calibration_source_sha256": "1" * 64,
        "original_simulation_source_sha256": "2" * 64,
        "original_slim_sha256": "3" * 64,
        "original_screen_contracts_audited": 39,
    }


def test_reconciliation_plan_and_merge_count_each_screen_draw_once(tmp_path):
    plan = _plan(tmp_path, draws=1)
    original = _screen_with_identity_failures(plan)
    reconciliation_plan = calibration.build_reconciliation_plan(
        plan,
        original,
        original_plan_sha256="a" * 64,
        original_ledger_sha256="b" * 64,
        screen_contract_binding=_reconciliation_binding(),
    )
    assert len(reconciliation_plan) == 3
    assert reconciliation_plan["calibration_id"].is_unique
    assert set(reconciliation_plan["seed"]).issubset(set(plan["seed"]))
    assert set(reconciliation_plan["original_error"]) == {
        calibration.RECONCILABLE_IDENTITY_ERROR
    }

    rows = []
    for row in reconciliation_plan.to_dict(orient="records"):
        rows.append(
            {
                **{key: row[key] for key in calibration.PLAN_COLUMNS},
                "reconciliation_id": row["reconciliation_id"],
                "reconciliation_mode": calibration.RECONCILIATION_MODE,
                "original_status": row["original_status"],
                "original_terminal_class": row["original_terminal_class"],
                "original_error": row["original_error"],
                "original_plan_sha256": row["original_plan_sha256"],
                "original_ledger_sha256": row["original_ledger_sha256"],
                "stdpopsim_recapitation_applied": False,
                "post_slim_neutral_mutation_overlay_applied": False,
                "status": "complete",
                "evaluable": True,
                "terminal_class": "segregating",
                "candidate_pool_alt_count": 100,
                "candidate_pool_af": 0.1,
                "candidate_pool_frequency_gate_passed": True,
                "exact_panel_feasible": True,
                "elapsed_seconds": 1.0,
                "error": "",
            }
        )
    canonical = calibration.merge_reconciliation_results(
        plan,
        original,
        reconciliation_plan,
        pd.DataFrame(rows),
    )
    assert len(canonical) == len(plan)
    assert canonical["calibration_id"].is_unique
    assert canonical["evaluable"].map(bool).all()
    assert canonical["canonical_result_source"].value_counts().to_dict() == {
        "original_screen": len(plan) - 3,
        "reconciled_forward_only_no_overlay": 3,
    }


def test_reconciliation_rejects_unrelated_screen_failure(tmp_path):
    plan = _plan(tmp_path, draws=1)
    original = _screen_with_identity_failures(plan, count=1)
    original.loc[original.index[0], "error"] = "RuntimeError: unrelated"
    with pytest.raises(ValueError, match="outside the approved"):
        calibration.build_reconciliation_plan(
            plan,
            original,
            original_plan_sha256="a" * 64,
            original_ledger_sha256="b" * 64,
            screen_contract_binding=_reconciliation_binding(),
        )


def _raw_selected_tree() -> tskit.TreeSequence:
    tables = tskit.TableCollection(sequence_length=2)
    tables.metadata_schema = tskit.MetadataSchema.permissive_json()
    tables.metadata = {
        "SLiM": {"cycle": 5000, "tick": 5000, "user_metadata": {"Q": [5.0]}}
    }
    tables.mutations.metadata_schema = tskit.MetadataSchema.permissive_json()
    population = tables.populations.add_row()
    samples = []
    for _ in range(4):
        individual = tables.individuals.add_row()
        samples.extend(
            tables.nodes.add_row(
                flags=tskit.NODE_IS_SAMPLE,
                time=0,
                population=population,
                individual=individual,
            )
            for _ in range(2)
        )
    selected_ancestor = tables.nodes.add_row(time=1, population=population)
    root = tables.nodes.add_row(time=600, population=population)
    for sample in samples[:4]:
        tables.edges.add_row(0, 2, selected_ancestor, sample)
    tables.edges.add_row(0, 2, root, selected_ancestor)
    for sample in samples[4:]:
        tables.edges.add_row(0, 2, root, sample)
    site = tables.sites.add_row(position=1, ancestral_state="")
    tables.mutations.add_row(
        site=site,
        node=selected_ancestor,
        time=480,
        derived_state="0",
        metadata={
            "mutation_list": [
                {
                    "mutation_type": 1,
                    "selection_coeff": 0.0,
                    "subpopulation": 7,
                    "slim_time": 4520,
                    "nucleotide": -1,
                }
            ]
        },
    )
    tables.sort()
    return tables.tree_sequence()


def test_raw_identity_and_carriers_are_validated_before_overlay():
    expectation = {
        "focal_position_0based": 1,
        "expected_mutation_type": 1,
        "expected_selection_coeff": 0.0,
        "source_subpopulation": 7,
        "realized_origin_time_generations": 2400.0,
        "slim_scaling_factor": 5.0,
    }
    identity = calibration._inspect_raw_selected_focal_identity(
        _raw_selected_tree(), expectation
    )
    assert identity["status"] == "valid"
    assert identity["observed_origin_time_raw_scaled_generations"] == 480
    assert identity["observed_origin_time_generations"] == 2400
    assert calibration._selected_counts_from_raw_tree(
        _raw_selected_tree(), identity
    ).tolist() == [2, 2, 0, 0]


def _completed_followup_ledger(plan: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for row in plan.to_dict(orient="records"):
        hit = True
        if row["calibration_phase"] == "extension":
            hit = (
                np.isclose(row["duration_multiplier"], 0.60)
                and int(row["draw_index"]) < 4
            )
        target = float(row["target_allele_frequency"])
        af = target if hit else min(0.95, target + 0.40)
        result = {
            **row,
            "stdpopsim_recapitation_applied": row["calibration_phase"]
            == "sensitivity_10mb",
            "post_slim_neutral_mutation_overlay_applied": row["calibration_phase"]
            == "sensitivity_10mb",
            "status": "complete",
            "evaluable": True,
            "terminal_class": "segregating",
            "candidate_pool_alt_count": round(1000 * af),
            "candidate_pool_af": af,
            "candidate_pool_frequency_gate_passed": hit,
            "exact_panel_feasible": hit,
            "elapsed_seconds": 1.0,
            "error": "",
        }
        if row["calibration_phase"] == "sensitivity_10mb":
            result.update(
                {
                    "sensitivity_draw_passed": hit,
                    "normal_production_processing_validated": True,
                    "production_exact_panel_artifact_validated": hit,
                }
            )
        rows.append(result)
    return pd.DataFrame(rows)


def test_followup_plans_extend_only_failures_and_use_disjoint_seed_namespaces(
    tmp_path,
):
    screen_plan = _plan(tmp_path, draws=5)
    canonical = _synthetic_ledger(screen_plan)
    failed = np.isclose(canonical["selection_coefficient"], 0.005) & np.isclose(
        canonical["target_allele_frequency"], 0.30
    )
    canonical.loc[failed, "candidate_pool_frequency_gate_passed"] = False
    canonical.loc[failed, "exact_panel_feasible"] = False
    canonical.loc[failed, "candidate_pool_alt_count"] = 800
    canonical.loc[failed, "candidate_pool_af"] = 0.8
    screen_seeds = set(screen_plan["seed"].astype(int))
    production = [
        seed for seed in range(2**31 - 1, 2**31 - 100, -1) if seed not in screen_seeds
    ][:2]
    reserved = {"production": production, "screen": screen_plan["seed"]}
    extension = calibration.build_extension_plan(
        screen_plan,
        canonical,
        upstream_artifact_sha256="c" * 64,
        reserved_seed_groups=reserved,
    )
    assert len(extension) == 80
    assert calibration._calibration_cells(extension) == {(0.005, 0.30)}
    assert set(extension["duration_multiplier"]) == {0.60, 0.70, 0.85, 1.00}
    assert set(extension["simulation_mode"]) == {calibration.FOCAL_ONLY_FORWARD_MODE}
    assert set(extension["sequence_length_bp"]) == {2}
    assert set(extension["focal_position_bp"]) == {1}

    extension_ledger = _completed_followup_ledger(extension)
    _, _, _, selected = calibration.combine_screen_and_extension(
        screen_plan,
        canonical,
        extension,
        extension_ledger,
    )
    chosen = selected[
        np.isclose(selected["selection_coefficient"], 0.005)
        & np.isclose(selected["target_allele_frequency"], 0.30)
    ].iloc[0]
    assert chosen["duration_multiplier"] == pytest.approx(0.60)
    assert chosen["band_hit_rate"] == pytest.approx(0.20)

    confirmation_reserved = {
        **reserved,
        "extension": extension["seed"],
    }
    confirmation = calibration.build_confirmation_plan(
        selected,
        upstream_artifact_sha256="d" * 64,
        reserved_seed_groups=confirmation_reserved,
    )
    assert len(confirmation) == 600
    assert set(confirmation["simulation_mode"]) == {calibration.FOCAL_ONLY_FORWARD_MODE}
    sensitivity_reserved = {
        **confirmation_reserved,
        "confirmation": confirmation["seed"],
    }
    sensitivity = calibration.build_sensitivity_plan(
        selected,
        upstream_artifact_sha256="d" * 64,
        reserved_seed_groups=sensitivity_reserved,
    )
    assert len(sensitivity) == 120
    assert set(sensitivity["simulation_mode"]) == {
        calibration.PRODUCTION_10MB_SENSITIVITY_MODE
    }
    assert set(sensitivity["sequence_length_bp"]) == {10_000_000}
    assert set(sensitivity["focal_position_bp"]) == {5_000_000}
    audit = calibration.audit_disjoint_seed_groups(
        {**sensitivity_reserved, "sensitivity_10mb": sensitivity["seed"]}
    )
    assert audit["status"] == "pairwise_disjoint"
    assert set(audit["groups"]) == {
        "production",
        "screen",
        "extension",
        "confirmation",
        "sensitivity_10mb",
    }

    confirmation_summary = calibration.summarize_endpoint_validation_phase(
        confirmation, _completed_followup_ledger(confirmation)
    )
    sensitivity_summary = calibration.summarize_endpoint_validation_phase(
        sensitivity, _completed_followup_ledger(sensitivity)
    )
    assert confirmation_summary["phase_gate_passed"].all()
    assert sensitivity_summary["phase_gate_passed"].all()


def test_seed_audit_rejects_cross_phase_reuse():
    with pytest.raises(ValueError, match="not disjoint"):
        calibration.audit_disjoint_seed_groups(
            {"screen": [101, 102], "confirmation": [102, 103]}
        )


def test_followup_finalizer_requires_and_binds_all_phase_evidence(
    tmp_path, monkeypatch
):
    output = tmp_path / "calibration"
    output.mkdir()
    screen_plan = _plan(tmp_path, draws=5)
    original = _synthetic_ledger(screen_plan)
    failure_index = original.index[0]
    original.loc[failure_index, "status"] = "failed"
    original.loc[failure_index, "evaluable"] = False
    original.loc[failure_index, "terminal_class"] = "failed"
    original.loc[failure_index, "candidate_pool_alt_count"] = np.nan
    original.loc[failure_index, "candidate_pool_af"] = np.nan
    original.loc[failure_index, "candidate_pool_frequency_gate_passed"] = False
    original.loc[failure_index, "exact_panel_feasible"] = False
    original.loc[failure_index, "error"] = calibration.RECONCILABLE_IDENTITY_ERROR
    canonical = original.copy()
    canonical.loc[failure_index, "status"] = "complete"
    canonical.loc[failure_index, "evaluable"] = True
    canonical.loc[failure_index, "terminal_class"] = "segregating"
    canonical.loc[failure_index, "candidate_pool_alt_count"] = 100
    canonical.loc[failure_index, "candidate_pool_af"] = 0.1
    canonical.loc[failure_index, "candidate_pool_frequency_gate_passed"] = True
    canonical.loc[failure_index, "exact_panel_feasible"] = True
    canonical.loc[failure_index, "error"] = ""
    failed_cell = np.isclose(canonical["selection_coefficient"], 0.005) & np.isclose(
        canonical["target_allele_frequency"], 0.30
    )
    canonical.loc[failed_cell, "candidate_pool_frequency_gate_passed"] = False
    canonical.loc[failed_cell, "exact_panel_feasible"] = False
    canonical.loc[failed_cell, "candidate_pool_alt_count"] = 800
    canonical.loc[failed_cell, "candidate_pool_af"] = 0.8

    source_dir = tmp_path / "python/gamma_smc_aou"
    source_dir.mkdir(parents=True)
    calibration_source = source_dir / "focused_selection_calibration.py"
    simulation_source = source_dir / "focused_selection_simulation.py"
    calibration_source.write_bytes(b"finalizer calibration source\n")
    simulation_source.write_bytes(b"simulation source\n")
    monkeypatch.setattr(calibration, "__file__", str(calibration_source))
    source_hash = simulation.sha256_file(calibration_source)
    simulation_hash = simulation.sha256_file(simulation_source)
    snapshots = output / "source_snapshots"
    snapshots.mkdir()
    (snapshots / f"focused_selection_calibration.py.sha256_{source_hash}").write_bytes(
        calibration_source.read_bytes()
    )
    reconciliation_plan = calibration.build_reconciliation_plan(
        screen_plan,
        original,
        original_plan_sha256="a" * 64,
        original_ledger_sha256="b" * 64,
        screen_contract_binding={
            "original_calibration_source_sha256": source_hash,
            "original_simulation_source_sha256": simulation_hash,
            "original_slim_sha256": "c" * 64,
            "original_screen_contracts_audited": len(screen_plan) - 1,
        },
    )
    screen_seeds = set(screen_plan["seed"].astype(int))
    production = [
        seed for seed in range(2**31 - 1, 2**31 - 100, -1) if seed not in screen_seeds
    ][:2]
    reserved = {"production": production, "screen": screen_plan["seed"]}
    extension = calibration.build_extension_plan(
        screen_plan,
        canonical,
        upstream_artifact_sha256="d" * 64,
        reserved_seed_groups=reserved,
    )
    extension_ledger = _completed_followup_ledger(extension)
    _, _, _, selected = calibration.combine_screen_and_extension(
        screen_plan, canonical, extension, extension_ledger
    )
    confirmation = calibration.build_confirmation_plan(
        selected,
        upstream_artifact_sha256="e" * 64,
        reserved_seed_groups={**reserved, "extension": extension["seed"]},
    )
    confirmation_ledger = _completed_followup_ledger(confirmation)
    sensitivity = calibration.build_sensitivity_plan(
        selected,
        upstream_artifact_sha256="e" * 64,
        reserved_seed_groups={
            **reserved,
            "extension": extension["seed"],
            "confirmation": confirmation["seed"],
        },
    )
    sensitivity_ledger = _completed_followup_ledger(sensitivity)
    sensitivity_work = tmp_path / "sensitivity-work"
    for index, row in sensitivity_ledger.iterrows():
        directory = sensitivity_work / str(row["calibration_id"])
        directory.mkdir(parents=True)
        tree = directory / "accepted_panel.trees"
        manifest = directory / "accepted_panel_manifest.tsv"
        tree.write_bytes(f"tree:{row['calibration_id']}\n".encode())
        manifest.write_bytes(f"manifest:{row['calibration_id']}\n".encode())
        sensitivity_ledger.loc[index, "panel_tree_path"] = tree.name
        sensitivity_ledger.loc[index, "panel_tree_sha256"] = simulation.sha256_file(
            tree
        )
        sensitivity_ledger.loc[index, "panel_manifest_path"] = manifest.name
        sensitivity_ledger.loc[index, "panel_manifest_sha256"] = simulation.sha256_file(
            manifest
        )
    screen_plan.to_csv(output / "han_selection_end_plan.tsv", sep="\t", index=False)
    original.to_csv(
        output / "han_selection_end_terminal_draws.tsv", sep="\t", index=False
    )
    audit_paths = []
    for name, status in {
        "han_selection_end_reconciliation_audit.json": "complete",
        "han_selection_end_endpoint_selection_audit.json": "complete",
        "han_selection_end_confirmation_audit.json": "passed",
        "han_selection_end_sensitivity_10mb_audit.json": "passed",
    }.items():
        path = output / name
        path.write_text(json.dumps({"status": status}) + "\n")
        audit_paths.append(path)
    payload = calibration.finalize_followup_calibration(
        screen_plan=screen_plan,
        original_screen_ledger=original,
        canonical_screen_ledger=canonical,
        reconciliation_plan=reconciliation_plan,
        extension_plan=extension,
        extension_ledger=extension_ledger,
        selected_artifact=selected,
        confirmation_plan=confirmation,
        confirmation_ledger=confirmation_ledger,
        sensitivity_plan=sensitivity,
        sensitivity_ledger=sensitivity_ledger,
        output_dir=output,
        sensitivity_work_root=sensitivity_work,
        production_seeds=production,
        artifact_paths=audit_paths,
    )
    assert payload["status"] == "frozen"
    assert payload["production_contract"]["confirmation_passed"] is True
    assert payload["production_contract"]["ten_mb_sensitivity_passed"] is True
    assert len(payload["mapping"]) == 6
    assert (output / "han_selection_end_frozen.json").is_file()
    assert (output / "han_selection_end_sensitivity_10mb_panel_artifacts.tsv").is_file()
