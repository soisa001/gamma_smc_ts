from __future__ import annotations

import json
import math
import shutil
from pathlib import Path

import pandas as pd
import pytest
from gamma_smc_aou import focused_eas_selected_calibration as calibration

SOURCE_REPO = Path(__file__).resolve().parents[1]


def _isolated_repo(tmp_path: Path) -> tuple[Path, Path]:
    root = tmp_path / "repo"
    for relative in calibration.IMPLEMENTATION_SOURCE_PATHS:
        destination = root / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(SOURCE_REPO / relative, destination)
    resource = root / calibration.EAS_RESOURCE_PATH
    resource.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(SOURCE_REPO / calibration.EAS_RESOURCE_PATH, resource)
    slim = root / "bin" / "slim"
    slim.parent.mkdir(parents=True)
    slim.write_bytes(b"test-slim-binary")
    return root, slim


def _production_units() -> pd.DataFrame:
    rows = []
    seed = 1_000_000
    for coefficient in calibration.DEFAULT_SELECTION_COEFFICIENTS:
        for target in calibration.DEFAULT_TARGET_FREQUENCIES:
            rows.append(
                {
                    "unit_id": f"production-{coefficient}-{target}",
                    "demography_id": "eas_phlash_median",
                    "simulation_class": "selected",
                    "selection_coefficient": coefficient,
                    "target_allele_frequency": target,
                    "seed": seed,
                }
            )
            seed += 1
    return pd.DataFrame(rows)


def _bundle(tmp_path: Path):
    root, slim = _isolated_repo(tmp_path)
    units = _production_units()
    reservations = {
        "production": units["seed"].tolist(),
        "han": list(range(2_000_000, 2_001_000)),
    }
    plans, manifest = calibration.build_calibration_plans(
        calibration.EasSelectedCalibrationDesign(),
        repo_root=root,
        slim_path=slim,
        seed_reservations=reservations,
    )
    output = root / "focused_selection_EAS_sim" / "calibration" / "eas_selected"
    calibration.write_calibration_plans(output, plans, manifest)
    plans, manifest = calibration.read_plan_bundle(output)
    units_path = root / "focused_selection_EAS_sim" / "execution_units.tsv"
    units_path.parent.mkdir(parents=True, exist_ok=True)
    units.to_csv(units_path, sep="\t", index=False, lineterminator="\n")
    work = root / "focused_selection_EAS_sim" / "work" / "eas_selected_calibration"
    return root, slim, output, work, units_path, reservations, plans, manifest


def _clone_bundle(planned_bundle, tmp_path: Path):
    source_root, _, _, _, _, reservations, _, _ = planned_bundle
    root = tmp_path / "repo"
    shutil.copytree(source_root, root)
    slim = root / "bin" / "slim"
    output = root / "focused_selection_EAS_sim" / "calibration" / "eas_selected"
    plans, manifest = calibration.read_plan_bundle(output)
    units_path = root / "focused_selection_EAS_sim" / "execution_units.tsv"
    work = root / "focused_selection_EAS_sim" / "work" / "eas_selected_calibration"
    return root, slim, output, work, units_path, reservations, plans, manifest


@pytest.fixture(scope="module")
def planned_bundle(tmp_path_factory):
    return _bundle(tmp_path_factory.mktemp("eas-calibration-plan"))


def test_plans_have_exact_three_phase_grid_and_disjoint_seeds(planned_bundle):
    _, _, _, _, _, reservations, plans, manifest = planned_bundle
    assert {phase: len(plan) for phase, plan in plans.items()} == {
        calibration.PHASE_SCREEN: 120,
        calibration.PHASE_CONFIRM: 600,
        calibration.PHASE_SENSITIVITY: 120,
    }
    all_seeds = []
    for phase in calibration.PHASE_ORDER:
        plan = plans[phase]
        assert plan["seed"].is_unique
        assert set(
            plan.groupby(["selection_coefficient", "target_allele_frequency"]).size()
        ) == {calibration.PHASE_DRAWS[phase]}
        all_seeds.extend(int(value) for value in plan["seed"])
    assert len(all_seeds) == len(set(all_seeds)) == 840
    assert set(all_seeds).isdisjoint(
        set(reservations["production"]) | set(reservations["han"])
    )
    assert len(manifest["cells"]) == 6
    assert {
        key: value["origin_age_generations"] for key, value in manifest["cells"].items()
    } == {
        "s=0.010000|af=0.100000": 999.0,
        "s=0.010000|af=0.200000": 1152.0,
        "s=0.010000|af=0.300000": 1258.0,
        "s=0.005000|af=0.100000": 2038.0,
        "s=0.005000|af=0.200000": 2409.0,
        "s=0.005000|af=0.300000": 2657.0,
    }


def test_cell_contract_retains_continuous_selection_and_both_terminal_bounds(
    planned_bundle,
):
    *_, manifest = planned_bundle
    for cell in manifest["cells"].values():
        assert cell["continuous_selection_to_present"] is True
        assert cell["terminal_population_af_conditioning_in_slim"] is True
        events = cell["extended_events"]
        fitness = [
            event for event in events if event["event_type"] == "ChangeMutationFitness"
        ]
        assert len(fitness) == 1
        assert fitness[0]["end_time"]["generations_ago"] == 0
        terminal = {
            event["operator"]: event["allele_frequency"]
            for event in events
            if event["event_type"] == "ConditionOnAlleleFrequency"
            and event["start_time"]["generations_ago"] == 0
            and event["end_time"]["generations_ago"] == 0
        }
        assert terminal == {
            ">=": pytest.approx(cell["population_af_lower_inclusive"]),
            "<=": pytest.approx(cell["population_af_upper_inclusive"]),
        }


def test_plan_manifest_and_plan_checksum_tampering_is_rejected(
    planned_bundle, tmp_path
):
    _, _, output, _, _, _, _, _ = _clone_bundle(planned_bundle, tmp_path)
    path = output / f"{calibration.PHASE_SCREEN}_plan.tsv"
    original = path.read_bytes()
    path.write_bytes(original + b"\n")
    with pytest.raises(ValueError, match="checksum"):
        calibration.read_plan_bundle(output)


def test_seed_collision_with_later_han_phase_is_rejected(planned_bundle):
    *_, reservations, plans, _ = planned_bundle
    collided = {key: list(value) for key, value in reservations.items()}
    collided["han"].append(int(plans[calibration.PHASE_CONFIRM]["seed"].iloc[0]))
    with pytest.raises(ValueError, match="collides"):
        calibration.validate_seed_disjointness(plans, collided)


def test_design_rejects_any_tunable_scientific_grid():
    with pytest.raises(ValueError, match="Q=5"):
        calibration.EasSelectedCalibrationDesign(slim_scaling_factor=10).validate()
    with pytest.raises(ValueError, match="AF=0.10"):
        calibration.EasSelectedCalibrationDesign(
            target_frequencies=(0.1, 0.2)
        ).validate()


def _synthetic_completions(
    root: Path,
    slim: Path,
    work: Path,
    plans: dict[str, pd.DataFrame],
    manifest: dict,
    *,
    hit_rate: float = 0.20,
    phases: tuple[str, ...] = calibration.PHASE_ORDER,
    write_completions: bool = True,
) -> dict[str, pd.DataFrame]:
    ledgers = {}
    for phase in phases:
        plan = plans[phase]
        rows = []
        for row in plan.to_dict(orient="records"):
            hits_required = math.ceil(calibration.PHASE_DRAWS[phase] * hit_rate)
            hit = int(row["draw_index"]) < hits_required
            target = float(row["target_allele_frequency"])
            pool_diploids = int(row["candidate_pool_diploids"])
            alt_count = (
                round(2 * pool_diploids * target)
                if hit
                else round(2 * pool_diploids * (target + 0.05))
            )
            directory = work / phase / str(row["calibration_id"])
            if write_completions:
                directory.mkdir(parents=True, exist_ok=True)
            panel_tree_hash = ""
            panel_manifest_hash = ""
            exact = int(row["exact_sample_alt_count"])
            if phase == calibration.PHASE_SENSITIVITY and hit:
                if write_completions:
                    panel_tree = directory / "accepted_panel.trees"
                    panel_manifest = directory / "accepted_panel_manifest.tsv"
                    panel_tree.write_bytes(f"tree:{row['calibration_id']}".encode())
                    panel_manifest.write_text("sample\tcount\n0\t0\n", encoding="utf-8")
                    panel_tree_hash = calibration.sha256_file(panel_tree)
                    panel_manifest_hash = calibration.sha256_file(panel_manifest)
                else:
                    panel_tree_hash = "a" * 64
                    panel_manifest_hash = "b" * 64
                n_alt = exact // 4
                n_het = exact - 2 * n_alt
                n_ref = int(row["sample_diploids"]) - n_alt - n_het
                sample_alt_count = exact
                sample_af = exact / (2 * int(row["sample_diploids"]))
            else:
                n_ref = n_het = n_alt = None
                sample_alt_count = sample_af = None
            result = {
                **{column: row[column] for column in calibration.PLAN_COLUMNS},
                "status": "complete",
                "evaluable": True,
                "terminal_class": "segregating",
                "candidate_pool_alt_count": alt_count,
                "candidate_pool_af": alt_count / (2 * pool_diploids),
                "candidate_pool_frequency_gate_passed": hit,
                "exact_panel_feasible": hit,
                "pool_band_and_panel_hit": hit,
                "sample_alt_count": sample_alt_count,
                "sample_af": sample_af,
                "n_hom_ref": n_ref,
                "n_heterozygous": n_het,
                "n_hom_alt": n_alt,
                "stdpopsim_recapitation_applied": (
                    phase == calibration.PHASE_SENSITIVITY
                ),
                "post_slim_neutral_mutation_overlay_applied": (
                    phase == calibration.PHASE_SENSITIVITY
                ),
                "elapsed_seconds": 1.0,
                "trajectory_sha256": "",
                "panel_tree_sha256": panel_tree_hash,
                "panel_manifest_sha256": panel_manifest_hash,
                "error": "",
            }
            if write_completions:
                contract = calibration._draw_contract(
                    row, repo_root=root, slim_path=slim, manifest=manifest
                )
                payload = {
                    "schema": calibration.SCHEMA_VERSION,
                    "contract": contract,
                    "contract_sha256": calibration._canonical_sha256(contract),
                    "result": result,
                }
                (directory / "completion.json").write_text(
                    json.dumps(payload, sort_keys=True), encoding="utf-8"
                )
            rows.append(result)
        ledgers[phase] = pd.DataFrame(rows, columns=calibration.RESULT_COLUMNS)
    return ledgers


def test_summary_requires_at_least_twenty_percent_in_every_cell(
    planned_bundle, tmp_path
):
    root, slim, _, work, _, _, plans, manifest = _clone_bundle(planned_bundle, tmp_path)
    ledgers = _synthetic_completions(
        root, slim, work, plans, manifest, write_completions=False
    )
    for phase in calibration.PHASE_ORDER:
        summary = calibration.summarize_phase(plans[phase], ledgers[phase], phase)
        assert len(summary) == 6
        assert summary["passed"].all()
        assert set(summary["pool_band_and_panel_hit_rate"]) == {0.2}
    screen = ledgers[calibration.PHASE_SCREEN].copy()
    first_cell = screen[
        (screen["selection_coefficient"] == 0.01)
        & (screen["target_allele_frequency"] == 0.1)
        & screen["pool_band_and_panel_hit"]
    ].index[0]
    for column in (
        "candidate_pool_frequency_gate_passed",
        "exact_panel_feasible",
        "pool_band_and_panel_hit",
    ):
        screen.loc[first_cell, column] = False
    summary = calibration.summarize_phase(
        plans[calibration.PHASE_SCREEN], screen, calibration.PHASE_SCREEN
    )
    failed = summary[
        (summary["selection_coefficient"] == 0.01)
        & (summary["target_allele_frequency"] == 0.1)
    ].iloc[0]
    assert failed["pool_band_and_panel_hit_rate"] == pytest.approx(0.15)
    assert not bool(failed["passed"])


def test_freeze_loader_and_artifact_tamper_are_fail_closed(planned_bundle, tmp_path):
    root, slim, output, work, units_path, reservations, plans, manifest = _clone_bundle(
        planned_bundle, tmp_path
    )
    ledgers = _synthetic_completions(root, slim, work, plans, manifest)
    payload = calibration.freeze_calibration(
        plans,
        ledgers,
        manifest,
        repo_root=root,
        output_dir=output,
        work_root=work,
        slim_path=slim,
        execution_units_path=units_path,
        seed_reservations=reservations,
    )
    assert payload["status"] == "frozen"
    assert len(payload["cell_authorizations"]) == 6
    assert set(payload["phases"]) == set(calibration.PHASE_ORDER)
    assert all(phase["all_six_cells_passed"] for phase in payload["phases"].values())
    frozen = output / "eas_selected_calibration_frozen.json"
    loaded = calibration.load_frozen_eas_selected_calibration(
        frozen,
        repo_root=root,
        execution_units_path=units_path,
        slim_path=slim,
    )
    assert loaded["payload_sha256"] == payload["payload_sha256"]

    # The later, explicitly authorized production wiring must edit this one
    # executor file without invalidating the independently frozen origin/event
    # contracts.
    simulation_source = root / "python/gamma_smc_aou/focused_selection_simulation.py"
    simulation_source.write_text(
        simulation_source.read_text(encoding="utf-8") + "\n# production wiring\n",
        encoding="utf-8",
    )
    calibration.load_frozen_eas_selected_calibration(
        frozen,
        repo_root=root,
        execution_units_path=units_path,
        slim_path=slim,
    )

    summary = output / f"{calibration.PHASE_SCREEN}_summary.tsv"
    original = summary.read_bytes()
    summary.write_bytes(original + b"\n")
    with pytest.raises(ValueError, match="artifact checksum"):
        calibration.load_frozen_eas_selected_calibration(
            frozen,
            repo_root=root,
            execution_units_path=units_path,
            slim_path=slim,
        )


def test_completion_contract_tampering_is_rejected(planned_bundle, tmp_path):
    root, slim, _, work, _, _, plans, manifest = _clone_bundle(planned_bundle, tmp_path)
    ledgers = _synthetic_completions(
        root,
        slim,
        work,
        plans,
        manifest,
        phases=(calibration.PHASE_SCREEN,),
    )
    phase = calibration.PHASE_SCREEN
    identity = str(plans[phase]["calibration_id"].iloc[0])
    completion = work / phase / identity / "completion.json"
    payload = json.loads(completion.read_text(encoding="utf-8"))
    payload["contract"]["scientific_contract"]["selection"] = "changed"
    completion.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="contract differs"):
        calibration.audit_phase_completions(
            plans[phase],
            ledgers[phase],
            repo_root=root,
            work_root=work,
            slim_path=slim,
            manifest=manifest,
        )
