from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import stdpopsim
import tskit
from gamma_smc_aou import focused_han_frequency_conditioned_calibration as calibration
from gamma_smc_aou.focused_selection_simulation import sha256_file


def _fake_repo(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    root = tmp_path / "repo"
    root.mkdir()
    sources = []
    for name in ("module.py", "wrapper.py", "pyproject.toml", "uv.lock"):
        path = root / name
        path.write_text(f"source={name}\n", encoding="utf-8")
        sources.append(path)
    monkeypatch.setattr(calibration, "_implementation_paths", lambda _root: sources)
    monkeypatch.setattr(
        calibration, "_slim_version", lambda _path: "SLiM version 4.2.2"
    )
    slim = root / "slim"
    slim.write_bytes(b"test slim binary")
    campaign = root / "focused_selection_EAS_sim"
    fixed = campaign / "calibration"
    fixed.mkdir(parents=True)
    execution = campaign / "execution_units.tsv"
    pd.DataFrame(
        [
            {"unit_id": "selected", "simulation_class": "selected", "seed": 11},
            {"unit_id": "neutral", "simulation_class": "neutral", "seed": 12},
        ]
    ).to_csv(execution, sep="\t", index=False)
    exclusions = [execution]
    for index, name in enumerate(
        (
            "han_selection_end_plan.tsv",
            "han_selection_end_extension_plan.tsv",
            "han_selection_end_reconciliation_plan.tsv",
        ),
        start=20,
    ):
        path = fixed / name
        pd.DataFrame({"seed": [index]}).to_csv(path, sep="\t", index=False)
        exclusions.append(path)
    failed = []
    for name in ("fixed_draws.tsv", "NOT_FROZEN.json"):
        path = fixed / name
        path.write_text(f"failed fixed cessation: {name}\n", encoding="utf-8")
        failed.append(path)
    manifest_path = (
        fixed / "han_frequency_conditioned" / "han_frequency_conditioned_manifest.json"
    )
    return root, slim, exclusions, failed, manifest_path


def _build_manifest(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    root, slim, exclusions, failed, manifest_path = _fake_repo(tmp_path, monkeypatch)
    manifest = calibration.build_atomic_manifest(
        calibration.FrequencyConditionedDesign(),
        repo_root=root,
        slim_path=slim,
        excluded_seed_artifacts=exclusions,
        failed_fixed_cessation_artifacts=failed,
    )
    calibration.write_atomic_manifest(manifest, manifest_path, repo_root=root)
    loaded, frame, digest = calibration.load_manifest(manifest_path, repo_root=root)
    return root, manifest_path, loaded, frame, digest


def _row(*, phase="screen", coefficient=0.01, target=0.10, draw=0):
    focal = phase in {"screen", "confirmation"}
    identity = (
        f"han_frequency_conditioned__{phase}__s{coefficient:.6f}__"
        f"af{target:.6f}__q5__draw{draw:03d}"
    )
    return {
        "draw_id": identity,
        "phase": phase,
        "selection_coefficient": coefficient,
        "target_allele_frequency": target,
        "population_af_lower": target - 0.025,
        "population_af_upper": target + 0.025,
        "exact_sample_alt_count": round(200 * target),
        "sample_diploids": 100,
        "candidate_pool_diploids": 500,
        "slim_scaling_factor": 5.0,
        "mutation_age_generations": 2400.0,
        "selection_start_generations_ago": 2272.0,
        "selection_end_generations_ago": 0.0,
        "draw_index": draw,
        "seed": 1000 + draw,
        "seed_nonce": 0,
        "panel_seed": 2000 + draw,
        "panel_seed_nonce": 0,
        "timeout_seconds": 1800.0,
        "sequence_length_bp": 2 if focal else 10_000_000,
        "focal_position_bp": 1 if focal else 5_000_000,
        "processing_mode": (
            calibration.FOCAL_ONLY_MODE if focal else calibration.PRODUCTION_MODE
        ),
    }


def test_manifest_roundtrip_atomically_reserves_all_phase_rng_seeds(
    tmp_path, monkeypatch
):
    root, manifest_path, manifest, frame, digest = _build_manifest(
        tmp_path, monkeypatch
    )
    assert len(frame) == 840
    assert frame.groupby("phase").size().to_dict() == {
        "screen": 120,
        "confirmation": 600,
        "sensitivity_10mb": 120,
    }
    rng_seeds = pd.concat([frame["seed"], frame["panel_seed"]])
    assert len(rng_seeds) == 1680
    assert rng_seeds.is_unique
    assert manifest["seed_registry"]["fallback_union_count"] == 1680
    reservation = (
        root
        / "focused_selection_EAS_sim/calibration/han_frequency_conditioned_all_phase_plan.tsv"
    )
    assert reservation.is_file()
    assert manifest["all_phase_seed_reservation"]["sha256"] == sha256_file(reservation)
    projected = pd.read_csv(reservation, sep="\t")
    assert len(projected) == 840
    assert set(projected["seed"]) == set(frame["seed"])
    assert set(projected["panel_seed"]) == set(frame["panel_seed"])
    assert digest == sha256_file(manifest_path)
    original = json.loads(manifest_path.read_text(encoding="utf-8"))
    calibration.write_atomic_manifest(
        calibration.build_atomic_manifest(
            calibration.FrequencyConditionedDesign(),
            repo_root=root,
            slim_path=root / "slim",
            excluded_seed_artifacts=[
                root / "focused_selection_EAS_sim/execution_units.tsv",
                root
                / "focused_selection_EAS_sim/calibration/han_selection_end_plan.tsv",
                root
                / "focused_selection_EAS_sim/calibration/han_selection_end_extension_plan.tsv",
                root
                / "focused_selection_EAS_sim/calibration/han_selection_end_reconciliation_plan.tsv",
            ],
            failed_fixed_cessation_artifacts=[
                root / "focused_selection_EAS_sim/calibration/fixed_draws.tsv",
                root / "focused_selection_EAS_sim/calibration/NOT_FROZEN.json",
            ],
        ),
        manifest_path,
        repo_root=root,
    )
    assert json.loads(manifest_path.read_text(encoding="utf-8")) == original
    assert (
        calibration._require_repo_path(root / "safe", root, label="test")
        == (root / "safe").resolve()
    )
    with pytest.raises(ValueError, match="OneDrive"):
        calibration._require_repo_path(root / "OneDrive/output", root, label="test")


def test_manifest_detects_reservation_and_source_drift(tmp_path, monkeypatch):
    root, manifest_path, manifest, _, _ = _build_manifest(tmp_path, monkeypatch)
    reservation = root / manifest["all_phase_seed_reservation"]["path"]
    reservation.write_text(
        reservation.read_text(encoding="utf-8") + "tamper\n", encoding="utf-8"
    )
    with pytest.raises(ValueError, match="reservation"):
        calibration.load_manifest(manifest_path, repo_root=root)


def test_projection_is_published_only_after_manifest_commit(tmp_path, monkeypatch):
    root, slim, exclusions, failed, manifest_path = _fake_repo(tmp_path, monkeypatch)
    manifest = calibration.build_atomic_manifest(
        calibration.FrequencyConditionedDesign(),
        repo_root=root,
        slim_path=slim,
        excluded_seed_artifacts=exclusions,
        failed_fixed_cessation_artifacts=failed,
    )
    reservation = (
        root
        / "focused_selection_EAS_sim/calibration/han_frequency_conditioned_all_phase_plan.tsv"
    )
    original_atomic_json = calibration._atomic_json

    def checked_atomic_json(path, payload):
        if Path(path) == manifest_path:
            assert not reservation.exists()
        original_atomic_json(path, payload)

    monkeypatch.setattr(calibration, "_atomic_json", checked_atomic_json)
    calibration.write_atomic_manifest(manifest, manifest_path, repo_root=root)
    assert manifest_path.is_file()
    assert reservation.is_file()


def test_frequency_conditioned_event_topology_and_q5_realized_schedule():
    row = _row()
    _, sweep, samples, record = calibration.build_frequency_conditioned_objects(row)
    audit = record["event_audit"]
    assert len(sweep.extended_events) == audit["focal_event_count"] == 9
    assert samples == {"Han": 500}
    assert audit["mutation"] == {
        "population": "Neanderthal",
        "age_generations": 2400.0,
    }
    assert audit["recipient_nonloss_condition_count"] == 2
    assert audit["source_condition_count"] == 2
    assert audit["terminal_han_conditions"][">="] == pytest.approx(0.075)
    assert audit["terminal_han_conditions"]["<="] == pytest.approx(0.125)
    assert audit["realized_q5_event_times_generations"] == {
        "mutation_age_generations": 2400.0,
        "source_nonloss_generation_after_start": 2395.0,
        "source_check_generations": 2275.0,
        "introgression_and_selection_onset_generations": 2270.0,
        "loschbour_nonloss_generation_after_start": 2265.0,
        "demographic_han_split_generations": 2015.0,
        "selection_population_transition_generations": 2015.0,
        "han_nonloss_generation_after_start": 2010.0,
    }
    fitness = [
        event
        for event in sweep.extended_events
        if isinstance(event, stdpopsim.ChangeMutationFitness)
    ]
    assert [
        (event.population, event.start_time, event.end_time, event.dominance_coeff)
        for event in fitness
    ] == [
        ("Loschbour", 2272.0, 2015.0, 0.5),
        ("Han", 2015.0, 0.0, 0.5),
    ]
    invalid = replace(
        sweep,
        extended_events=sweep.extended_events + (sweep.extended_events[-1],),
    )
    with pytest.raises(ValueError, match="one focal DrawMutation"):
        calibration._audit_frequency_conditioned_events(invalid, row)


def _raw_pool_tree() -> tskit.TreeSequence:
    tables = tskit.TableCollection(sequence_length=2)
    tables.metadata_schema = tskit.MetadataSchema.permissive_json()
    tables.metadata = {
        "SLiM": {"cycle": 5000, "tick": 5000, "user_metadata": {"Q": [5.0]}}
    }
    tables.mutations.metadata_schema = tskit.MetadataSchema.permissive_json()
    population = tables.populations.add_row()
    individuals = []
    for _ in range(500):
        individual = tables.individuals.add_row()
        individuals.append(
            tuple(
                tables.nodes.add_row(
                    flags=tskit.NODE_IS_SAMPLE,
                    time=0,
                    population=population,
                    individual=individual,
                )
                for _ in range(2)
            )
        )
    carriers = {nodes[0] for nodes in individuals[:80]}
    carriers.update(node for nodes in individuals[80:90] for node in nodes)
    selected_ancestor = tables.nodes.add_row(time=1, population=population)
    root = tables.nodes.add_row(time=600, population=population)
    for nodes in individuals:
        for sample in nodes:
            parent = selected_ancestor if sample in carriers else root
            tables.edges.add_row(0, 2, parent, sample)
    tables.edges.add_row(0, 2, root, selected_ancestor)
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


def test_focal_only_uniform_exact_k_panel_uses_reserved_seed():
    pool = _raw_pool_tree()
    expectation = {
        "focal_position_0based": 1,
        "expected_mutation_type": 1,
        "expected_selection_coeff": 0.0,
        "source_subpopulation": 7,
        "realized_origin_time_generations": 2400.0,
        "slim_scaling_factor": 5.0,
    }
    identity = calibration._inspect_raw_type1(pool, expectation)
    counts = calibration._raw_type1_counts(pool, identity)
    assert len(counts) == 500
    assert int(counts.sum()) == 100
    panel, record = calibration._select_exact_sample_panel_at_position(
        pool,
        focal_position=1,
        focal_genotype_counts=counts,
        sample_diploids=100,
        exact_alt_count=20,
        seed=987654,
    )
    panel_identity = calibration._inspect_raw_type1(panel, expectation)
    observed = calibration._raw_type1_counts(panel, panel_identity)
    assert len(observed) == 100
    assert int(observed.sum()) == 20
    assert np.count_nonzero(observed == 0) >= 2
    assert np.count_nonzero(observed == 1) >= 1
    assert np.count_nonzero(observed == 2) >= 2
    assert record["panel_seed"] == 987654
    assert record["selection_method"] == (
        "uniform_over_all_size_N_exact_alt_count_subsets_passing_genotype_qc"
    )


def _synthetic_screen_ledger(plan: pd.DataFrame, work: Path) -> pd.DataFrame:
    rows = []
    for row in plan.to_dict(orient="records"):
        hit = int(row["draw_index"]) < 4
        directory = work / str(row["draw_id"])
        directory.mkdir(parents=True)
        panel_tree = directory / "accepted_panel.trees"
        panel_manifest = directory / "accepted_panel_manifest.tsv"
        if hit:
            panel_tree.write_bytes(b"synthetic tree")
            panel_manifest.write_text("synthetic\n", encoding="utf-8")
        rows.append(
            {
                **row,
                "status": "complete",
                "evaluable": True,
                "terminal_class": "segregating",
                "strict_type1_identity_passed": True,
                "candidate_pool_alt_count": int(
                    1000 * float(row["target_allele_frequency"])
                ),
                "candidate_pool_af": float(row["target_allele_frequency"]),
                "candidate_pool_frequency_gate_passed": hit,
                "exact_panel_feasible": hit,
                "panel_roundtrip_validated": hit,
                "hit": hit,
                "panel_tree_path": panel_tree.name if hit else "",
                "panel_tree_sha256": sha256_file(panel_tree) if hit else "",
                "panel_manifest_path": panel_manifest.name if hit else "",
                "panel_manifest_sha256": (sha256_file(panel_manifest) if hit else ""),
                "internal_terminal_han_af_conditioning": True,
                "selection_through_present": True,
                "recipient_nonloss_conditioning": True,
                "stdpopsim_recapitation_applied": False,
                "post_slim_neutral_mutation_overlay_applied": False,
                "elapsed_seconds": 1.0,
                "error": "",
            }
        )
    return pd.DataFrame(rows)


def test_phase_gate_uses_fixed_plan_denominator(tmp_path, monkeypatch):
    _, _, _, frame, _ = _build_manifest(tmp_path, monkeypatch)
    plan = frame[frame["phase"].eq("screen")].reset_index(drop=True)
    work = tmp_path / "repo/work"
    ledger = _synthetic_screen_ledger(plan, work)
    summary, audit = calibration.summarize_phase(plan, ledger, work_root=work)
    assert set(summary["n_hits"]) == {4}
    assert set(summary["hit_rate"]) == {0.20}
    assert audit["status"] == "passed"
    first_cell = ledger[
        np.isclose(ledger["selection_coefficient"], 0.005)
        & np.isclose(ledger["target_allele_frequency"], 0.10)
        & ledger["hit"]
    ].index[0]
    ledger.loc[first_cell, "hit"] = False
    ledger.loc[first_cell, "candidate_pool_frequency_gate_passed"] = False
    ledger.loc[first_cell, "exact_panel_feasible"] = False
    ledger.loc[first_cell, "panel_roundtrip_validated"] = False
    summary, audit = calibration.summarize_phase(plan, ledger, work_root=work)
    failed = summary[
        np.isclose(summary["selection_coefficient"], 0.005)
        & np.isclose(summary["target_allele_frequency"], 0.10)
    ].iloc[0]
    assert failed["n_hits"] == 3
    assert failed["hit_rate"] == 0.15
    assert audit["status"] == "failed"


def test_authoritative_ledger_rebuild_rejects_contract_tampering(tmp_path):
    root = tmp_path / "repo"
    root.mkdir()
    work = root / "work"
    row = _row(phase="sensitivity_10mb")
    directory = work / row["draw_id"]
    directory.mkdir(parents=True)
    trajectory = directory / "trajectory.csv"
    trajectory.write_text("tick,p1\n1,0.1\n", encoding="utf-8")
    _, _, _, model_record = calibration.build_frequency_conditioned_objects(row)
    event_audit = model_record["event_audit"]
    result = {
        **row,
        "status": "complete",
        "evaluable": True,
        "terminal_class": "identity_collision",
        "strict_type1_identity_passed": False,
        "candidate_pool_alt_count": None,
        "candidate_pool_af": None,
        "candidate_pool_frequency_gate_passed": False,
        "exact_panel_feasible": False,
        "panel_roundtrip_validated": False,
        "hit": False,
        "panel_tree_path": "",
        "panel_manifest_path": "",
        "internal_terminal_han_af_conditioning": True,
        "selection_through_present": True,
        "recipient_nonloss_conditioning": True,
        "stdpopsim_recapitation_applied": True,
        "post_slim_neutral_mutation_overlay_applied": True,
        "trajectory_path": trajectory.name,
        "trajectory_sha256": sha256_file(trajectory),
        "event_audit_json": json.dumps(event_audit),
    }
    sources = {"module.py": "a" * 64}
    runtime = {"slim": {"sha256": "b" * 64}}
    contract = calibration._draw_contract(
        row,
        manifest_sha256="c" * 64,
        implementation_sources=sources,
        runtime=runtime,
    )
    completion = directory / "completion.json"
    completion.write_text(
        json.dumps(
            {
                "schema": calibration.DRAW_SCHEMA,
                "contract": contract,
                "contract_sha256": calibration._canonical_sha256(contract),
                "result": result,
            }
        ),
        encoding="utf-8",
    )
    plan = pd.DataFrame([row], columns=calibration.ROW_COLUMNS)
    ledger = calibration.rebuild_phase_ledger(
        plan,
        repo_root=root,
        work_root=work,
        manifest_sha256="c" * 64,
        implementation_sources=sources,
        runtime=runtime,
    )
    assert ledger.loc[0, "terminal_class"] == "identity_collision"
    payload = json.loads(completion.read_text(encoding="utf-8"))
    payload["contract"]["model"]["selection_end_generations_ago"] = 10
    completion.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="contract differs"):
        calibration.rebuild_phase_ledger(
            plan,
            repo_root=root,
            work_root=work,
            manifest_sha256="c" * 64,
            implementation_sources=sources,
            runtime=runtime,
        )
