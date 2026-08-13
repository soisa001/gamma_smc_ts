from __future__ import annotations

import copy
import json
import shutil
from pathlib import Path

import msprime
import pandas as pd
import pyslim
import pytest
import tskit
from gamma_smc_aou import focused_han_fixed_v2_calibration as calibration
from gamma_smc_aou.focused_selection_simulation import sha256_file

SCREEN_SHA = calibration.EXPLORATORY_SCREEN_ARTIFACTS[
    calibration.EXPLORATORY_SCREEN_LEDGER_RELATIVE_PATH
]


def _screen_record():
    return {
        "status": "bound_exploratory_only",
        "canonical_rows": 840,
        "canonical_evaluable_rows": 840,
        "canonical_ledger_path": calibration.EXPLORATORY_SCREEN_LEDGER_RELATIVE_PATH,
        "canonical_ledger_sha256": SCREEN_SHA,
        "endpoint_table_path": calibration.EXPLORATORY_SELECTED_RELATIVE_PATH,
        "endpoint_table_sha256": calibration.EXPLORATORY_SCREEN_ARTIFACTS[
            calibration.EXPLORATORY_SELECTED_RELATIVE_PATH
        ],
        "endpoint_rule_changed_after_screen": True,
        "screen_draws_used_for_inferential_validation": False,
        "fresh_confirmation_and_sensitivity_only": True,
        "artifacts": {},
    }


def _seed_provenance_record(path: Path):
    return {
        "schema": calibration.EXPLORATORY_SEED_PROVENANCE_SCHEMA,
        "status": "bound_noninferential_seed_evidence",
        "canonical_reservation_path": path.as_posix(),
        "canonical_seed_count": 1,
        "canonical_seed_sha256": calibration._seed_digest([41]),
    }


def _fake_repo(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    root = tmp_path / "repo"
    root.mkdir()
    source = root / "source.py"
    source.write_text("source\n", encoding="utf-8")
    slim = root / "slim"
    slim.write_bytes(b"slim 4.2.2")
    campaign = root / "focused_selection_EAS_sim"
    calibration_dir = campaign / "calibration"
    calibration_dir.mkdir(parents=True)
    execution = campaign / "execution_units.tsv"
    pd.DataFrame(
        [
            {"unit_id": "selected", "simulation_class": "selected", "seed": 11},
            {"unit_id": "neutral", "simulation_class": "neutral", "seed": 12},
        ]
    ).to_csv(execution, sep="\t", index=False)
    fixed = calibration_dir / "fixed.tsv"
    pd.DataFrame({"seed": [21, 22], "panel_seed": [31, 32]}).to_csv(
        fixed, sep="\t", index=False
    )
    seed_reservation = (
        calibration_dir / "exploratory_seed_provenance/exploratory_seed_reservation.tsv"
    )
    seed_reservation.parent.mkdir(parents=True)
    pd.DataFrame({"seed": [41]}).to_csv(seed_reservation, sep="\t", index=False)
    monkeypatch.setattr(calibration, "EXPLORATORY_SEED_COUNT", 1)
    monkeypatch.setattr(
        calibration, "EXPLORATORY_SEED_SHA256", calibration._seed_digest([41])
    )
    monkeypatch.setattr(calibration, "EXPLORATORY_PREVIOUSLY_OMITTED_SEED_COUNT", 0)
    monkeypatch.setattr(
        calibration,
        "EXPLORATORY_PREVIOUSLY_OMITTED_SEED_SHA256",
        calibration._seed_digest([]),
    )
    monkeypatch.setattr(
        calibration, "_bind_exploratory_screen", lambda _root: _screen_record()
    )
    monkeypatch.setattr(
        calibration,
        "_bind_exploratory_seed_provenance",
        lambda _root: _seed_provenance_record(seed_reservation.relative_to(root)),
    )
    monkeypatch.setattr(
        calibration,
        "_implementation_sources",
        lambda _root: {"source.py": sha256_file(source)},
    )
    runtime = {
        "python": "test",
        "stdpopsim": "test",
        "msprime": "test",
        "pyslim": "test",
        "numpy": "test",
        "pandas": "test",
        "tskit": "test",
        "slim": {
            "path": "slim",
            "sha256": sha256_file(slim),
            "version": "SLiM version 4.2.2",
            "expected_version": "4.2.2",
        },
    }
    monkeypatch.setattr(calibration, "_runtime_record", lambda _root, _slim: runtime)
    manifest_path = (
        calibration_dir / "han_fixed_cessation_v2/han_fixed_v2_manifest.json"
    )
    return root, slim, [execution, fixed, seed_reservation], manifest_path, runtime


def _built_manifest(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    root, slim, exclusions, manifest_path, runtime = _fake_repo(tmp_path, monkeypatch)
    manifest = calibration.build_atomic_manifest(
        calibration.HanFixedV2Design(),
        repo_root=root,
        slim_path=slim,
        excluded_seed_artifacts=exclusions,
    )
    calibration.write_atomic_manifest(manifest, manifest_path, repo_root=root)
    loaded, frame, digest = calibration.load_manifest(manifest_path, repo_root=root)
    return root, manifest_path, loaded, frame, digest, runtime


def test_exploratory_screen_binding_and_post_screen_endpoints_are_exact():
    root = Path(__file__).resolve().parents[1]
    record = calibration._bind_exploratory_screen(root)
    assert record["canonical_rows"] == record["canonical_evaluable_rows"] == 840
    assert record["canonical_ledger_sha256"] == SCREEN_SHA
    assert record["endpoint_rule_changed_after_screen"] is True
    assert record["screen_draws_used_for_inferential_validation"] is False
    endpoints = calibration.endpoint_table().set_index(
        ["selection_coefficient", "target_allele_frequency"]
    )
    for cell, (multiplier, endpoint) in calibration.FIXED_ENDPOINTS.items():
        assert endpoints.loc[cell, "duration_multiplier"] == pytest.approx(multiplier)
        assert endpoints.loc[
            cell, "realized_selection_end_generations_ago"
        ] == pytest.approx(endpoint)


def test_exploratory_seed_archive_is_byte_bound_and_reserves_all_344_seeds():
    root = Path(__file__).resolve().parents[1]
    record = calibration._bind_exploratory_seed_provenance(root)
    assert record["canonical_seed_count"] == 344
    assert record["canonical_seed_sha256"] == (
        "22e28351535828559093f24461439c5ce406855657a904aca86599a764e45eb1"
    )
    assert record["previously_omitted_seed_count"] == 212
    assert record["previously_omitted_seed_sha256"] == (
        "3373a7db6c2ef5dc5d8d6e1036a8039f6b65d256db75129d66703ef926296f9d"
    )
    assert len(record["archive"]) == 8
    reservation = pd.read_csv(root / record["canonical_reservation_path"], sep="\t")
    assert reservation["seed"].is_unique
    assert set(reservation["prearchive_v2_exclusion_status"]) == {
        "already_reserved",
        "omitted",
    }


def test_exploratory_seed_archive_rejects_absence_and_tamper(tmp_path):
    source = (
        Path(__file__).resolve().parents[1]
        / calibration.EXPLORATORY_SEED_PROVENANCE_RELATIVE_DIR
    )
    root = tmp_path / "repo"
    destination = root / calibration.EXPLORATORY_SEED_PROVENANCE_RELATIVE_DIR
    shutil.copytree(source, destination)
    calibration._bind_exploratory_seed_provenance(root)
    reservation = root / calibration.EXPLORATORY_SEED_RESERVATION_RELATIVE_PATH
    original = reservation.read_bytes()
    reservation.write_bytes(original + b"tamper\n")
    with pytest.raises(ValueError, match="reservation"):
        calibration._bind_exploratory_seed_provenance(root)
    reservation.write_bytes(original)
    inventory = root / calibration.EXPLORATORY_SEED_INVENTORY_RELATIVE_PATH
    inventory.unlink()
    with pytest.raises(ValueError, match="inventory"):
        calibration._bind_exploratory_seed_provenance(root)


def test_manifest_source_and_runtime_bindings_cover_transitive_study_dependencies(
    tmp_path, monkeypatch
):
    relative = "python/gamma_smc_aou/eas_sweep_study.py"
    assert relative in calibration.IMPLEMENTATION_SOURCE_PATHS
    root = tmp_path / "repo"
    root.mkdir()
    study = root / relative
    study.parent.mkdir(parents=True)
    study.write_text("BASE_SEED = 1\n", encoding="utf-8")
    monkeypatch.setattr(calibration, "IMPLEMENTATION_SOURCE_PATHS", (relative,))
    sources = calibration._implementation_sources(root)
    study.write_text("BASE_SEED = 2\n", encoding="utf-8")
    assert calibration._implementation_sources(root) != sources
    slim = root / "slim"
    slim.write_bytes(b"slim")
    monkeypatch.setattr(calibration, "_slim_version", lambda _path: "SLiM 4.2.2")
    runtime = calibration._runtime_record(root, slim)
    assert runtime["msprime"] == msprime.__version__
    assert runtime["pyslim"] == pyslim.__version__
    monkeypatch.setattr(calibration.msprime, "__version__", "drifted")
    assert calibration._runtime_record(root, slim) != runtime


def test_gitattributes_preserve_every_byte_bound_source():
    root = Path(__file__).resolve().parents[1]
    rules = {
        line.split()[0]: line.split()[1:]
        for line in (root / ".gitattributes").read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    }
    missing = {
        path
        for path in calibration.IMPLEMENTATION_SOURCE_PATHS
        if "-text" not in rules.get(path, [])
    }
    assert not missing


def test_manifest_atomically_reserves_fresh_explicit_sim_and_panel_seeds(
    tmp_path, monkeypatch
):
    root, manifest_path, manifest, frame, digest, _ = _built_manifest(
        tmp_path, monkeypatch
    )
    assert len(frame) == 720
    assert frame.groupby("phase").size().to_dict() == {
        "confirmation": 600,
        "sensitivity_10mb": 120,
    }
    all_rng = pd.concat([frame["seed"], frame["panel_seed"]])
    assert len(all_rng) == 1440
    assert all_rng.is_unique
    assert manifest["seed_registry"]["fixed_v2_union_count"] == 1440
    assert manifest["inference_contract"] == {
        "endpoint_rule_prespecified_before_exploratory_screen": False,
        "endpoint_rule_changed_after_exploratory_screen": True,
        "exploratory_screen_is_descriptive_only": True,
        "fresh_confirmation_and_sensitivity_are_inferential": True,
        "confirmation_all_evaluable_and_minimum_hits_per_cell": 10,
        "sensitivity_all_evaluable_and_minimum_hits_per_cell": 2,
    }
    reservation = root / manifest["all_phase_seed_reservation"]["path"]
    assert reservation.is_file()
    assert sha256_file(reservation) == manifest["all_phase_seed_reservation"]["sha256"]
    assert digest == sha256_file(manifest_path)
    before = manifest_path.read_bytes()
    calibration.write_atomic_manifest(manifest, manifest_path, repo_root=root)
    assert manifest_path.read_bytes() == before


def test_manifest_rejects_reservation_drift_and_onedrive(tmp_path, monkeypatch):
    root, manifest_path, manifest, _, _, _ = _built_manifest(tmp_path, monkeypatch)
    reservation = root / manifest["all_phase_seed_reservation"]["path"]
    reservation.write_text(reservation.read_text() + "tamper\n", encoding="utf-8")
    with pytest.raises(ValueError, match="reservation"):
        calibration.load_manifest(manifest_path, repo_root=root)
    with pytest.raises(ValueError, match="OneDrive"):
        calibration._require_repo_path(
            root / "OneDrive/output", root, label="test output"
        )


def test_manifest_rejects_canonical_exploratory_seed_tamper(tmp_path, monkeypatch):
    root, manifest_path, _, _, _, _ = _built_manifest(tmp_path, monkeypatch)
    reservation = root / calibration.EXPLORATORY_SEED_RESERVATION_RELATIVE_PATH
    reservation.write_bytes(reservation.read_bytes() + b"42\n")
    with pytest.raises(ValueError, match="excluded seed evidence drifted"):
        calibration.load_manifest(manifest_path, repo_root=root)


def _row(phase="confirmation", coefficient=0.01, target=0.10, draw=0):
    endpoint = (
        calibration.endpoint_table()
        .set_index(["selection_coefficient", "target_allele_frequency"])
        .loc[(coefficient, target)]
    )
    focal = phase == "confirmation"
    return {
        "draw_id": (
            f"han_fixed_v2__{phase}__s{coefficient:.6f}__af{target:.6f}__"
            f"q5__draw{draw:03d}"
        ),
        "phase": phase,
        "selection_coefficient": coefficient,
        "target_allele_frequency": target,
        "population_af_lower": target - 0.025,
        "population_af_upper": target + 0.025,
        "duration_multiplier": float(endpoint["duration_multiplier"]),
        "requested_selection_end_generations_ago": float(
            endpoint["requested_selection_end_generations_ago"]
        ),
        "realized_selection_end_generations_ago": float(
            endpoint["realized_selection_end_generations_ago"]
        ),
        "exact_sample_alt_count": round(200 * target),
        "sample_diploids": 100,
        "candidate_pool_diploids": 500,
        "slim_scaling_factor": 5.0,
        "mutation_age_generations": 2400.0,
        "selection_start_generations_ago": 2272.0,
        "draw_index": draw,
        "seed": 100000 + draw,
        "seed_nonce": 0,
        "panel_seed": 200000 + draw,
        "panel_seed_nonce": 0,
        "timeout_seconds": 1800.0,
        "sequence_length_bp": 2 if focal else 10_000_000,
        "focal_position_bp": 1 if focal else 5_000_000,
        "processing_mode": (
            calibration.FOCAL_ONLY_MODE if focal else calibration.PRODUCTION_MODE
        ),
        "exploratory_screen_ledger_sha256": SCREEN_SHA,
    }


def _sensitivity_identity_expectation():
    return {
        "schema": calibration.focused_simulation.SELECTED_FOCAL_IDENTITY_SCHEMA,
        "focal_position_0based": 5_000_000,
        "expected_mutation_type": 1,
        "expected_selection_coeff": 0.0,
        "source_population": "Neanderthal",
        "source_subpopulation": 7,
        "requested_origin_time_generations": 2_400.0,
        "realized_origin_time_generations": 2_400.0,
        "slim_scaling_factor": 5.0,
        "slim_time_rule": "cycle - realized_origin_time / Q",
    }


def _sensitivity_identity_tree(
    *,
    selected: bool = True,
    overlay_type: int | None = 0,
    selected_source: int = 7,
    malformed_overlay: bool = False,
):
    ancestry = msprime.sim_ancestry(
        samples=2,
        ploidy=2,
        sequence_length=10_000_000,
        recombination_rate=0,
        population_size=100_000,
        random_seed=73,
    )
    annotated = pyslim.annotate(ancestry, model_type="WF", tick=5_000, stage="late")
    tables = annotated.dump_tables()
    top_metadata = copy.deepcopy(tables.metadata)
    top_metadata["SLiM"]["user_metadata"] = {"Q": [5.0]}
    tables.metadata = top_metadata
    if malformed_overlay:
        tables.mutations.metadata_schema = tskit.MetadataSchema.permissive_json()
    site = tables.sites.add_row(position=5_000_000, ancestral_state="A")
    node = int(annotated.samples()[0])
    if selected:
        tables.mutations.add_row(
            site=site,
            node=node,
            time=2_400.0,
            derived_state="G",
            metadata={
                "mutation_list": [
                    {
                        "mutation_type": 1,
                        "selection_coeff": 0.0,
                        "subpopulation": selected_source,
                        "slim_time": 4_520,
                        "nucleotide": 2,
                    }
                ]
            },
        )
    if overlay_type is not None:
        metadata = (
            {}
            if malformed_overlay
            else {
                "mutation_list": [
                    {
                        "mutation_type": overlay_type,
                        "selection_coeff": 0.0,
                        "subpopulation": -1,
                        "slim_time": 4_999,
                        "nucleotide": 3,
                    }
                ]
            }
        )
        tables.mutations.add_row(
            site=site,
            node=node,
            time=5.0,
            derived_state="T",
            metadata=metadata,
        )
    tables.sort()
    tables.build_index()
    tables.compute_mutation_parents()
    return tables.tree_sequence()


def test_both_phases_require_retained_valid_type1_and_rethrow_absence(
    tmp_path, monkeypatch
):
    sensitivity = _row("sensitivity_10mb")
    with pytest.raises(
        calibration.focused_simulation.SelectedFocalIdentityError,
        match="type-1 identity is required",
    ):
        calibration._evaluate_tree(
            _sensitivity_identity_tree(selected=False),
            sensitivity,
            _sensitivity_identity_expectation(),
            {},
            tmp_path,
        )

    monkeypatch.setattr(
        calibration,
        "_inspect_raw_selected_focal_identity",
        lambda _ts, _expectation: {"status": "selected_focal_mutation_absent"},
    )
    with pytest.raises(
        calibration.focused_simulation.SelectedFocalIdentityError,
        match="confirmation requires",
    ):
        calibration._evaluate_tree(object(), _row("confirmation"), {}, {}, tmp_path)


def test_sensitivity_only_accepts_structured_type0_overlay_recurrence(tmp_path):
    tree = _sensitivity_identity_tree()
    evidence = calibration._inspect_neutral_overlay_recurrence(
        tree, _sensitivity_identity_expectation()
    )
    calibration._validate_identity_collision_evidence(evidence)
    assert evidence["selected_type1_valid"] is True
    assert evidence["neutral_overlay_only"] is True
    assert evidence["selected_identity"]["status"] == "valid"
    result = calibration._evaluate_tree(
        tree,
        _row("sensitivity_10mb"),
        _sensitivity_identity_expectation(),
        {},
        tmp_path,
    )
    assert result["terminal_class"] == "identity_collision"
    assert result["strict_type1_identity_passed"] is False
    assert result["hit"] is False
    assert json.loads(result["identity_json"])["schema"] == (
        calibration.IDENTITY_COLLISION_SCHEMA
    )


@pytest.mark.parametrize(
    "tree",
    [
        _sensitivity_identity_tree(overlay_type=2),
        _sensitivity_identity_tree(selected_source=6),
        _sensitivity_identity_tree(malformed_overlay=True),
    ],
)
def test_wrong_or_malformed_sensitivity_identity_is_operational_failure(tmp_path, tree):
    with pytest.raises(calibration.focused_simulation.SelectedFocalIdentityError):
        calibration._evaluate_tree(
            tree,
            _row("sensitivity_10mb"),
            _sensitivity_identity_expectation(),
            {},
            tmp_path,
        )


def test_collision_validator_rejects_unstructured_evidence():
    evidence = calibration._inspect_neutral_overlay_recurrence(
        _sensitivity_identity_tree(), _sensitivity_identity_expectation()
    )
    evidence.pop("selected_identity")
    with pytest.raises(ValueError, match="not structured"):
        calibration._validate_identity_collision_evidence(evidence)


def test_fixed_event_topology_uses_external_af_gate_and_prescribed_endpoint():
    root = Path(__file__).resolve().parents[1]
    for coefficient, target in calibration.FIXED_ENDPOINTS:
        row = _row(coefficient=coefficient, target=target)
        _, sweep, samples, record = calibration.build_fixed_v2_objects(
            row, repo_root=root
        )
        audit = record["event_audit"]
        assert samples == {"Han": 500}
        assert audit["mutation_population"] == "Neanderthal"
        assert audit["mutation_age_generations"] == 2400.0
        assert audit["terminal_han_condition_count"] == 0
        assert audit["external_terminal_ascertainment"] is True
        assert audit["recipient_nonloss_condition_count"] == 2
        assert audit["source_condition_count"] == 2
        assert audit["realized_selection_end_generations_ago"] == pytest.approx(
            calibration.FIXED_ENDPOINTS[(coefficient, target)][1]
        )
        fitness = [
            event
            for event in sweep.extended_events
            if isinstance(event, __import__("stdpopsim").ChangeMutationFitness)
        ]
        assert len(fitness) == audit["selection_interval_count"]


def _synthetic_plan_and_ledger(phase: str, hits_per_cell: int):
    draws = calibration.PHASE_DRAW_COUNTS[phase]
    plan_rows = []
    ledger_rows = []
    for coefficient in calibration.SELECTION_COEFFICIENTS:
        for target in calibration.TARGET_FREQUENCIES:
            for draw in range(draws):
                row = _row(phase, coefficient, target, draw)
                # Make cell seeds unique across s/AF for synthetic validation.
                offset = int(coefficient * 1_000_000 + target * 100_000)
                row["seed"] += offset
                row["panel_seed"] += offset
                plan_rows.append(row)
                ledger_rows.append(
                    {
                        **row,
                        "evaluable": True,
                        "terminal_class": "segregating",
                        "candidate_pool_af": target,
                        "candidate_pool_frequency_gate_passed": draw < hits_per_cell,
                        "exact_panel_feasible": draw < hits_per_cell,
                        "hit": draw < hits_per_cell,
                        "elapsed_seconds": 1.0,
                        "error": "",
                    }
                )
    return (
        pd.DataFrame(plan_rows, columns=calibration.ROW_COLUMNS),
        pd.DataFrame(ledger_rows),
    )


@pytest.mark.parametrize(
    ("phase", "passing", "failing"),
    [("confirmation", 10, 9), ("sensitivity_10mb", 2, 1)],
)
def test_phase_gate_requires_all_evaluable_and_exact_hit_floor(phase, passing, failing):
    plan, ledger = _synthetic_plan_and_ledger(phase, passing)
    summary, audit = calibration.summarize_phase(plan, ledger)
    assert audit["status"] == "passed"
    assert set(summary["n_evaluable"]) == {calibration.PHASE_DRAW_COUNTS[phase]}
    assert set(summary["n_hits"]) == {passing}
    _, failing_ledger = _synthetic_plan_and_ledger(phase, failing)
    _, failed = calibration.summarize_phase(plan, failing_ledger)
    assert failed["status"] == "failed"
    incomplete = ledger.copy()
    incomplete.loc[incomplete.index[0], "evaluable"] = False
    incomplete.loc[incomplete.index[0], "hit"] = False
    _, failed = calibration.summarize_phase(plan, incomplete)
    assert failed["status"] == "failed"


def test_draw_contract_binds_explicit_rng_and_fresh_validation_only():
    row = _row("sensitivity_10mb")
    contract = calibration._draw_contract(
        row,
        manifest_sha256="a" * 64,
        implementation_sources={"source": "b" * 64},
        runtime={"slim": {"sha256": "c" * 64}},
    )
    assert contract["processing"]["explicit_simulation_seed"] == row["seed"]
    assert contract["processing"]["explicit_panel_seed"] == row["panel_seed"]
    assert contract["model"]["exploratory_screen_inferential_use"] is False
    assert contract["model"]["post_screen_endpoint_selection"] is True
    assert contract["processing"]["stdpopsim_recapitation_applied"] is True


def test_default_exclusions_cover_production_fixed_fallback_and_eas_streams():
    root = Path(__file__).resolve().parents[1]
    campaign = root / "focused_selection_EAS_sim"
    paths = calibration._default_seed_exclusion_paths(campaign)
    relative = {path.relative_to(campaign).as_posix() for path in paths}
    assert {
        "execution_units.tsv",
        "calibration/exploratory_seed_provenance/exploratory_seed_reservation.tsv",
        "calibration/han_selection_end_plan.tsv",
        "calibration/han_selection_end_extension_plan.tsv",
        "calibration/han_selection_end_reconciliation_plan.tsv",
        "calibration/han_frequency_conditioned_all_phase_plan.tsv",
        "calibration/eas_selected/screen20_plan.tsv",
        "calibration/eas_selected/confirm100_plan.tsv",
        "calibration/eas_selected/sensitivity10mb20_plan.tsv",
    }.issubset(relative)


def test_default_exclusions_fail_closed_when_canonical_reservation_is_absent(
    tmp_path,
):
    campaign = tmp_path / "focused_selection_EAS_sim"
    for relative in (
        "execution_units.tsv",
        "calibration/han_selection_end_plan.tsv",
        "calibration/han_selection_end_extension_plan.tsv",
        "calibration/han_selection_end_reconciliation_plan.tsv",
        "calibration/han_frequency_conditioned_all_phase_plan.tsv",
        "calibration/eas_selected/screen20_plan.tsv",
        "calibration/eas_selected/confirm100_plan.tsv",
        "calibration/eas_selected/sensitivity10mb20_plan.tsv",
    ):
        path = campaign / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"seed": [1]}).to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match="exploratory_seed_reservation.tsv"):
        calibration._default_seed_exclusion_paths(campaign)


def test_default_exclusions_ignore_own_published_projection_on_plan_rerun(tmp_path):
    campaign = tmp_path / "focused_selection_EAS_sim"
    required = (
        "execution_units.tsv",
        "calibration/exploratory_seed_provenance/exploratory_seed_reservation.tsv",
        "calibration/han_selection_end_plan.tsv",
        "calibration/han_selection_end_extension_plan.tsv",
        "calibration/han_selection_end_reconciliation_plan.tsv",
        "calibration/han_frequency_conditioned_all_phase_plan.tsv",
        "calibration/eas_selected/screen20_plan.tsv",
        "calibration/eas_selected/confirm100_plan.tsv",
        "calibration/eas_selected/sensitivity10mb20_plan.tsv",
        "calibration/han_fixed_cessation_v2/han_fixed_v2_all_phase_plan.tsv",
    )
    for relative in required:
        path = campaign / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"seed": [1]}).to_csv(path, sep="\t", index=False)
    paths = calibration._default_seed_exclusion_paths(campaign)
    own = (
        campaign / "calibration/han_fixed_cessation_v2/han_fixed_v2_all_phase_plan.tsv"
    ).resolve()
    assert own not in paths
    assert len(paths) == 9


def test_finalize_is_idempotent_and_existing_production_loader_compatible(
    tmp_path, monkeypatch
):
    root = tmp_path / "repo"
    root.mkdir()
    campaign = root / "focused_selection_EAS_sim"
    output = campaign / "calibration/han_fixed_cessation_v2"
    output.mkdir(parents=True)
    execution = campaign / "execution_units.tsv"
    unit_rows = []
    for index, (coefficient, target) in enumerate(calibration.FIXED_ENDPOINTS):
        unit_rows.append(
            {
                "seed": index + 1,
                "demography_id": "ancient_eurasia_han_introgression",
                "simulation_class": "selected",
                "selection_coefficient": coefficient,
                "target_allele_frequency": target,
                "population_af_lower": target - 0.025,
                "population_af_upper": target + 0.025,
            }
        )
    pd.DataFrame(unit_rows).to_csv(execution, sep="\t", index=False)
    plans = {}
    ledgers = {}
    combined_plans = []
    for phase, hits in (("confirmation", 10), ("sensitivity_10mb", 2)):
        plan, ledger = _synthetic_plan_and_ledger(phase, hits)
        plans[phase] = plan
        ledgers[phase] = ledger
        combined_plans.append(plan)
    frame = pd.concat(combined_plans, ignore_index=True)
    source = root / "source.py"
    source.write_text("source\n", encoding="utf-8")
    source_hash = sha256_file(source)
    runtime = {"slim": {"path": "slim", "sha256": "c" * 64}}
    manifest = {
        "runtime": runtime,
        "implementation_sources": {"source.py": source_hash},
        "source_snapshots": {},
        "seed_registry": {"schema": calibration.SEED_REGISTRY_SCHEMA},
        "production_execution_units": {
            "path": execution.relative_to(root).as_posix(),
            "sha256": sha256_file(execution),
            "production_seed_digest": calibration._seed_digest(
                pd.DataFrame(unit_rows)["seed"]
            ),
        },
    }
    manifest_path = output / "han_fixed_v2_manifest.json"
    manifest_path.write_text("{}\n", encoding="utf-8")
    manifest_sha = sha256_file(manifest_path)

    for phase in calibration.PHASE_DRAW_COUNTS:
        summary, audit = calibration.summarize_phase(plans[phase], ledgers[phase])
        draws_path = output / f"{phase}_draws.tsv"
        summary_path = output / f"{phase}_summary.tsv"
        ledgers[phase].to_csv(draws_path, sep="\t", index=False)
        summary.to_csv(summary_path, sep="\t", index=False)
        payload = {
            **audit,
            "manifest_sha256": manifest_sha,
            "runtime": runtime,
            "implementation_sources": manifest["implementation_sources"],
            "exploratory_screen_ledger_sha256": SCREEN_SHA,
            "artifacts": {
                draws_path.name: sha256_file(draws_path),
                summary_path.name: sha256_file(summary_path),
            },
        }
        (output / f"{phase}_audit.json").write_text(
            json.dumps(payload, sort_keys=True), encoding="utf-8"
        )

    monkeypatch.setattr(
        calibration,
        "load_manifest",
        lambda *_args, **_kwargs: (manifest, frame, manifest_sha),
    )
    monkeypatch.setattr(
        calibration,
        "rebuild_phase_ledger",
        lambda phase_plan, **_kwargs: ledgers[str(phase_plan["phase"].iloc[0])],
    )
    phase_roots = {
        phase: campaign / f"work/han_fixed_cessation_v2/{phase}"
        for phase in calibration.PHASE_DRAW_COUNTS
    }
    first = calibration.finalize_fixed_v2(
        manifest_path=manifest_path,
        repo_root=root,
        output_dir=output,
        phase_work_roots=phase_roots,
    )
    frozen = output / "han_selection_end_frozen.json"
    auxiliary = output / "han_fixed_v2_authorization.json"
    assert frozen.is_file() and auxiliary.is_file()
    assert first["production_integration_status"] == "not_wired_into_campaign_executor"
    assert first["exploratory_screen_inferential_use"] is False
    before = {
        path.name: path.read_bytes() for path in output.iterdir() if path.is_file()
    }
    second = calibration.finalize_fixed_v2(
        manifest_path=manifest_path,
        repo_root=root,
        output_dir=output,
        phase_work_roots=phase_roots,
    )
    assert second == first
    assert before == {
        path.name: path.read_bytes() for path in output.iterdir() if path.is_file()
    }
    loaded = calibration.load_frozen_fixed_v2(
        auxiliary, repo_root=root, manifest_path=manifest_path
    )
    assert loaded == first
    frozen_payload = json.loads(frozen.read_text())
    assert frozen_payload["schema"] == calibration.HAN_SELECTION_CALIBRATION_SCHEMA
    assert len(frozen_payload["mapping"]) == 6
    frozen_payload["mapping"]["s=0.010000|af=0.100000"][
        "realized_selection_end_generations_ago"
    ] = 2060.0
    frozen.write_text(json.dumps(frozen_payload), encoding="utf-8")
    with pytest.raises(ValueError, match="frozen artifact differs"):
        calibration.load_frozen_fixed_v2(
            auxiliary, repo_root=root, manifest_path=manifest_path
        )
