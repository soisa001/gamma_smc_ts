from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from gamma_smc_aou.cosi2_framework import (
    COSI2_COMMIT,
    HAN_CHB_SPLIT_GENERATIONS,
    HAN_DEMOGRAPHY_ID,
    HAN_MIGRATION_END_GENERATIONS,
    HAN_NEANDERTHAL_MIGRATION_RATE,
    HAN_OOA_SPLIT_GENERATIONS,
    HAN_SYMMETRIC_CONTACT_ANCESTRY_EXPECTATION,
    FrameworkPlan,
    _han_bar_offsets,
    _plot_framework,
    _present_af_token,
    _readme,
    _recombination_map_text,
    _run_commands,
    _selection_coefficient_index,
    _selection_plot_style,
    build_age_response,
    build_han_event_ledger,
    build_model_cells,
    derive_eas_birth_age,
    derive_han_selection_onset,
    lint_parameter_text,
    load_eas_history,
    render_eas_parameter_file,
    render_han_parameter_file,
    rle_eas_history,
    validate_framework,
    validate_saved_trajectory,
)

REPO_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def history() -> pd.DataFrame:
    return load_eas_history(REPO_ROOT)


@pytest.fixture(scope="module")
def plan() -> FrameworkPlan:
    return FrameworkPlan()


@pytest.fixture(scope="module")
def cells(history: pd.DataFrame, plan: FrameworkPlan) -> pd.DataFrame:
    return build_model_cells(history, plan)


def test_eas_history_is_exact_half_up_integer_rle(history: pd.DataFrame):
    rle = rle_eas_history(history)
    assert len(history) == 10_001
    assert len(rle) == 4_667
    assert int(rle.iloc[0]["cosi2_integer_ne"]) == 37_637
    assert int(rle.iloc[-1]["cosi2_integer_ne"]) == 13_138
    assert history["time_generations"].is_monotonic_increasing


@pytest.mark.parametrize(
    ("coefficient", "expected_age", "expected_ne"),
    [(0.01, 1_475, 3_197), (0.005, 3_146, 5_214)],
)
def test_eas_fixed_point_centers(
    history: pd.DataFrame,
    coefficient: float,
    expected_age: int,
    expected_ne: int,
):
    record = derive_eas_birth_age(
        history, selection_coefficient=coefficient, target_af=0.20
    )
    assert record["age_generations"] == expected_age
    assert record["diploid_ne_at_birth"] == expected_ne
    assert record["deterministic_final_af"] == pytest.approx(0.20, abs=0.001)


@pytest.mark.parametrize(("coefficient", "expected_age"), [(0.01, 605), (0.005, 1_210)])
def test_han_target_bridge_duration_is_diagnostic(
    coefficient: float, expected_age: int
):
    record = derive_han_selection_onset(
        selection_coefficient=coefficient, initial_af=0.012, target_af=0.20
    )
    assert record["age_generations"] == expected_age
    assert record["age_years"] == expected_age * 29
    assert record["age_kya"] == pytest.approx(expected_age * 29 / 1_000)
    assert record["deterministic_final_af"] == pytest.approx(0.20, abs=0.001)


def test_model_cells_separate_genomewide_ancestry_from_allele_af(
    cells: pd.DataFrame, plan: FrameworkPlan
):
    assert len(cells) == 2
    assert set(cells["selection_coefficient"]) == {0.01}
    np.testing.assert_allclose(cells["present_af_lower"], 0.195)
    np.testing.assert_allclose(cells["present_af_upper"], 0.205)
    han = cells[cells["sample_population"] == "Han"]
    assert set(han["demography_id"]) == {HAN_DEMOGRAPHY_ID}
    assert set(han["selection_population"]) == {"CHB"}
    assert set(han["selection_onset_generations_ago"]) == {
        HAN_MIGRATION_END_GENERATIONS
    }
    assert set(han["source_generation_time_years"]) == {29.0}
    assert set(han["introgression_ancestry_fraction"]) == {0.012}
    assert set(han["introgressed_allele_af_target"]) == {0.012}
    np.testing.assert_allclose(han["introgressed_allele_af_lower"], 0.011)
    np.testing.assert_allclose(han["introgressed_allele_af_upper"], 0.013)
    assert set(han["archaic_migration_rate_per_generation"]) == {
        HAN_NEANDERTHAL_MIGRATION_RATE
    }
    assert han["introgression_enforcement"].str.contains("posthoc").all()
    assert set(cells.loc[cells["sample_population"] == "EAS", "seed"]) == {20240524}
    assert set(han["seed"]) == {20240624}
    assert _selection_coefficient_index(0.01) == 1
    assert _selection_plot_style(0.01) == ("#2171b5", "s")
    endpoints = dict(
        zip(
            han["selection_coefficient"],
            han["deterministic_present_af_from_gate_center"],
            strict=True,
        )
    )
    assert endpoints[0.01] == pytest.approx(0.234015, abs=1e-6)
    assert _present_af_token(plan) == "0.195-0.205"
    assert _recombination_map_text(plan) == "1\t1e-08\n"


def test_han_event_ledger_is_exact_5r19_continuous_migration(plan: FrameworkPlan):
    ledger = build_han_event_ledger(plan)
    assert "admix" not in set(ledger["event_type"])
    assert (ledger["event_type"] == "migration_rate").sum() == 14
    contact = ledger[
        ledger["label"].isin(["CEU_Neanderthal_contact_a", "CHB_Neanderthal_contact_a"])
    ]
    assert len(contact) == 2
    np.testing.assert_allclose(
        contact["migration_rate"], HAN_NEANDERTHAL_MIGRATION_RATE
    )
    assert set(contact["cosi2_time_generations"]) == {HAN_MIGRATION_END_GENERATIONS}
    assert ledger["cosi2_time_generations"].is_monotonic_increasing
    assert HAN_CHB_SPLIT_GENERATIONS in set(ledger["cosi2_time_generations"])
    assert HAN_OOA_SPLIT_GENERATIONS in set(ledger["cosi2_time_generations"])
    boundary = ledger[
        ledger["label"].isin(["OOA_bottleneck", "YRI_OOA_a", "YRI_OOA_b"])
    ]
    assert set(boundary["cosi2_time_generations"]) == {HAN_CHB_SPLIT_GENERATIONS}
    assert ledger["rounding_delta_years"].abs().max() < 12
    assert HAN_SYMMETRIC_CONTACT_ANCESTRY_EXPECTATION == pytest.approx(
        0.011806645, abs=1e-10
    )


def test_rendered_parameter_files_have_modern_sweep_contract(
    history: pd.DataFrame, cells: pd.DataFrame, plan: FrameworkPlan
):
    rle = rle_eas_history(history)
    events = build_han_event_ledger(plan)
    for row in cells.to_dict(orient="records"):
        if row["sample_population"] == "EAS":
            text = render_eas_parameter_file(rle, row, plan)
            assert text.count("pop_event change_size") == 4_666
            assert "pop_define 1 EAS" in text
        else:
            text = render_han_parameter_file(events, row, plan)
            assert "pop_define 4 Neanderthal" in text
            assert "pop_define 3 Han_CHB" in text
            assert {
                line for line in text.splitlines() if line.startswith("sample_size ")
            } == {
                "sample_size 1 0",
                "sample_size 2 0",
                "sample_size 3 200",
                "sample_size 4 0",
                "sample_size 5 0",
            }
            assert "pop_event admix" not in text
            assert "0.011-0.013 and is applied post hoc" in text
            assert (
                'pop_event migration_rate "CHB_Neanderthal_contact_a" 3 4 645 8.25e-06'
            ) in text
            expected_sweep = (
                f'pop_event sweep_mult_standing "{row["cell_id"]}" 4 2400 '
                f"{row['selection_coefficient']:g} 0.5 0.195-0.205 3 645"
            )
            assert expected_sweep in text
        audit = lint_parameter_text(text, expected_cell=row)
        assert audit["endpoint_token"] == "0.195-0.205"
        assert audit["n_sample_requests"] == (
            1 if row["sample_population"] == "EAS" else 5
        )
        assert audit["n_positive_sample_requests"] == 1
        assert text.count("pop_event sweep_mult_standing") == 1
        assert COSI2_COMMIT in text


def test_linter_rejects_wrong_han_selection_population(
    cells: pd.DataFrame, plan: FrameworkPlan
):
    row = cells[cells["sample_population"] == "Han"].iloc[0].to_dict()
    text = render_han_parameter_file(build_han_event_ledger(plan), row, plan)
    broken = text.replace("0.195-0.205 3 645", "0.195-0.205 2 645")
    with pytest.raises(ValueError, match="birth/selection population IDs"):
        lint_parameter_text(broken, expected_cell=row)


def test_linter_rejects_wrong_present_af_interval(
    cells: pd.DataFrame, plan: FrameworkPlan
):
    row = cells[cells["sample_population"] == "Han"].iloc[0].to_dict()
    text = render_han_parameter_file(build_han_event_ledger(plan), row, plan)
    broken = text.replace("0.195-0.205 3 645", "0.1-0.9 3 645")
    with pytest.raises(ValueError, match="differs from the model cell"):
        lint_parameter_text(broken, expected_cell=row)


def test_linter_requires_zero_sample_requests_for_unsampled_5r19_populations(
    cells: pd.DataFrame, plan: FrameworkPlan
):
    row = cells[cells["sample_population"] == "Han"].iloc[0].to_dict()
    text = render_han_parameter_file(build_han_event_ledger(plan), row, plan)
    broken = text.replace("sample_size 1 0\n", "")
    with pytest.raises(ValueError, match="sample_size request for every"):
        lint_parameter_text(broken, expected_cell=row)


def test_age_response_contains_both_models_and_coefficients(
    history: pd.DataFrame, plan: FrameworkPlan
):
    response = build_age_response(history, plan)
    assert set(response["selection_coefficient"]) == {0.01}
    assert set(response["demography_id"]) == {
        "eas_phlash_median",
        HAN_DEMOGRAPHY_ID,
    }
    assert response["deterministic_final_af"].between(0, 1).all()
    assert np.isfinite(response.select_dtypes(include=["number"])).all(axis=None)


def test_smoke_commands_use_top_level_cosi2_binary(cells: pd.DataFrame):
    commands = _run_commands(cells)
    assert "../../external/cosi2/coalescent" in commands
    assert "../../external/cosi2/cosi/coalescent" not in commands
    assert "MAXATTEMPTS=${COSI_MAXATTEMPTS:-10000}" in commands
    assert 'COSI_MAXATTEMPTS="$MAXATTEMPTS"' in commands
    assert commands.count("COSI_SAVE_TRAJ=") == 2
    assert commands.count("validate-trajectory") == 2
    assert "Syntax smoke only" in commands
    assert "does not implement the required outer rejection loop" in commands
    assert "RUN_ROOT=../cosi2_framework_runs" in commands
    assert 'mkdir -p "$RUN_ROOT/smoke"' in commands
    assert commands.count('COSI_SAVE_TRAJ="$RUN_ROOT/smoke/') == 2
    assert commands.count('--trajectory "$RUN_ROOT/smoke/') == 2
    assert "mkdir -p work/" not in commands
    assert "COSI_SAVE_TRAJ=work/" not in commands


def test_sibling_run_root_preserves_exact_bundle_inventory(tmp_path: Path):
    source = REPO_ROOT / "focused_selection_EAS_sim" / "cosi2_framework"
    bundle = tmp_path / "cosi2_framework"
    shutil.copytree(source, bundle)
    before = {
        path.relative_to(bundle).as_posix()
        for path in bundle.rglob("*")
        if path.is_file()
    }
    run_output = tmp_path / "cosi2_framework_runs" / "smoke" / "eas_s0p01"
    run_output.mkdir(parents=True)
    (run_output / "simulation.ms").write_text("pilot\n", encoding="utf-8")
    result = validate_framework(REPO_ROOT, output_dir=bundle)
    after = {
        path.relative_to(bundle).as_posix()
        for path in bundle.rglob("*")
        if path.is_file()
    }
    assert result["status"] == "complete"
    assert before == after
    assert (run_output / "simulation.ms").is_file()


def test_readme_does_not_promise_one_shot_han_acceptance(
    cells: pd.DataFrame, plan: FrameworkPlan
):
    readme = _readme(plan, cells)
    assert "one joint acceptance per 28 million" in readme
    assert "does not promise that its one Han draw" in readme
    assert "explicit `sample_size 0` requests" in readme
    assert "zero requests are required CoSi2" in readme
    assert "bookkeeping; they do not add sampled chromosomes" in readme


def test_saved_han_trajectory_requires_endpoint_and_migration_end_gate(
    tmp_path: Path, cells: pd.DataFrame
):
    cell = cells[cells["sample_population"] == "Han"].iloc[0].to_dict()
    generations = np.asarray([2_400, 2_093, 1_241, 645, 0], dtype=float)
    frame = pd.DataFrame(
        {
            "sim": 1,
            "gen": generations,
            "selfreq_2": [0.0, 0.005, 0.01, 0.01, 0.03],
            "selfreq_3": [0.0, 0.0, 0.005, 0.012, 0.20],
            "selfreq_4": [1 / 7_200, 0.08, 0.05, 0.02, 0.0],
        }
    )
    path = tmp_path / "trajectory.tsv"
    frame.to_csv(path, sep="\t", index=False)
    result = validate_saved_trajectory(path, cell=cell)
    assert result["status"] == "valid"
    assert result["selection_start_gate_passed"] is True

    broken = frame.copy()
    broken.loc[broken["gen"] == 645, "selfreq_3"] = 0.04
    broken.to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match="migration-end allele-AF gate"):
        validate_saved_trajectory(path, cell=cell)


def test_saved_eas_trajectory_uses_numeric_population_id(
    tmp_path: Path, cells: pd.DataFrame
):
    cell = cells[cells["sample_population"] == "EAS"].iloc[0].to_dict()
    frame = pd.DataFrame(
        {
            "sim": 1,
            "gen": [3_146, 1_000, 0],
            "selfreq_1": [1 / (2 * 5_214), 0.05, 0.20],
        }
    )
    path = tmp_path / "eas_trajectory.tsv"
    frame.to_csv(path, sep="\t", index=False)
    result = validate_saved_trajectory(path, cell=cell)
    assert result["status"] == "valid"
    assert result["selection_start_gate_required"] is False


def test_each_han_simulation_requires_its_own_boundary_row(
    tmp_path: Path, cells: pd.DataFrame
):
    cell = cells[cells["sample_population"] == "Han"].iloc[0].to_dict()
    first = pd.DataFrame(
        {
            "sim": 1,
            "gen": [2_400, 645, 0],
            "selfreq_3": [0.0, 0.012, 0.20],
            "selfreq_4": [1 / 7_200, 0.01, 0.0],
        }
    )
    second = pd.DataFrame(
        {
            "sim": 2,
            "gen": [2_400, 0],
            "selfreq_3": [0.0, 0.20],
            "selfreq_4": [1 / 7_200, 0.0],
        }
    )
    path = tmp_path / "multi_sim.tsv"
    pd.concat([first, second], ignore_index=True).to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match="exactly one migration-end row"):
        validate_saved_trajectory(path, cell=cell)


def test_plan_rejects_exact_zero_tolerance_and_invalid_han_birth():
    with pytest.raises(ValueError, match="present_af_tolerance"):
        FrameworkPlan(present_af_tolerance=0.0).validate()
    with pytest.raises(ValueError, match="Han mutation birth"):
        FrameworkPlan(han_mutation_birth_generations=2_000).validate()
    with pytest.raises(ValueError, match="integer generation"):
        FrameworkPlan(han_mutation_birth_generations=2_400.5).validate()
    with pytest.raises(ValueError, match="selection-start allele-AF gate"):
        FrameworkPlan(han_selection_start_af_tolerance=0.02).validate()


def test_model_cells_preserve_legacy_seed_by_selection_coefficient(
    history: pd.DataFrame,
):
    cells = build_model_cells(history, FrameworkPlan(selection_coefficients=(0.005,)))
    assert len(cells) == 2
    assert set(cells["sample_population"]) == {"EAS", "Han"}
    assert set(cells["selection_coefficient"]) == {0.005}
    assert set(cells.loc[cells["sample_population"] == "EAS", "seed"]) == {20240523}
    assert set(cells.loc[cells["sample_population"] == "Han", "seed"]) == {20240623}


def test_custom_af_contract_and_additional_coefficient_render_dynamically(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    history: pd.DataFrame,
):
    custom = FrameworkPlan(
        selection_coefficients=(0.001,),
        present_af=0.18,
        present_af_tolerance=0.01,
        han_selection_start_af_target=0.02,
        han_selection_start_af_tolerance=0.005,
    )
    custom_cells = build_model_cells(history, custom)
    custom_response = build_age_response(history, custom)
    eas = custom_cells[custom_cells["sample_population"] == "EAS"].iloc[0]
    eas_curve = custom_response[custom_response["demography_id"] == "eas_phlash_median"]
    marker = eas_curve[
        eas_curve["time_generations"] == eas["mutation_birth_generations_ago"]
    ]
    assert len(marker) == 1
    assert float(marker.iloc[0]["deterministic_final_af"]) == pytest.approx(
        custom.present_af, abs=0.001
    )
    han = custom_cells[custom_cells["sample_population"] == "Han"].iloc[0].to_dict()
    text = render_han_parameter_file(build_han_event_ledger(custom), han, custom)
    assert "0.17-0.19 3 645" in text
    assert "0.015-0.025 and is applied post hoc" in text

    snapshots: list[list[tuple[float, float]]] = []

    def fake_savefig(figure, *_args, **_kwargs):
        snapshots.append([axis.get_xlim() for axis in figure.axes])

    monkeypatch.setattr("matplotlib.figure.Figure.savefig", fake_savefig)
    _plot_framework(
        tmp_path / "unused.png",
        tmp_path / "unused.pdf",
        history=history,
        age_response=custom_response,
        cells=custom_cells,
        plan=custom,
    )
    assert len(snapshots) == 2
    assert all(lower == pytest.approx(0.0) for lower, _upper in snapshots[0])


def test_plot_time_axes_put_present_on_left(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    history: pd.DataFrame,
    cells: pd.DataFrame,
    plan: FrameworkPlan,
):
    np.testing.assert_array_equal(_han_bar_offsets(1), [0.0])
    np.testing.assert_allclose(_han_bar_offsets(2), [-0.11, 0.11])
    snapshots: list[list[tuple[float, float]]] = []

    def fake_savefig(figure, *_args, **_kwargs):
        snapshots.append([axis.get_xlim() for axis in figure.axes])

    monkeypatch.setattr("matplotlib.figure.Figure.savefig", fake_savefig)
    _plot_framework(
        tmp_path / "unused.png",
        tmp_path / "unused.pdf",
        history=history,
        age_response=build_age_response(history, plan),
        cells=cells,
        plan=plan,
    )
    assert len(snapshots) == 2
    for lower, upper in snapshots[0]:
        assert lower == pytest.approx(0.0)
        assert upper > lower
