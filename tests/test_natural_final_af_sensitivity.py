from __future__ import annotations

import copy
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from gamma_smc_aou import natural_final_af as base
from gamma_smc_aou import natural_final_af_sensitivity as sensitivity
from gamma_smc_aou.natural_final_af_sensitivity import (
    COMBINED_SELECTION_COEFFICIENTS,
    DEFAULT_EAS_S001_ATTEMPTS,
    DEFAULT_EAS_S005_ATTEMPTS,
    DEFAULT_HAN_S001_ATTEMPTS,
    DEFAULT_HAN_S005_ATTEMPTS,
    SensitivityPlan,
    build_cells,
    build_comparisons,
    combine_analysis_tables,
)


ROOT = Path(__file__).resolve().parents[1]
CAMPAIGN = ROOT / base.DEFAULT_CAMPAIGN_DIR


def _production_plan() -> SensitivityPlan:
    return SensitivityPlan(repo_root=ROOT, campaign_dir=CAMPAIGN)


def _tiny_cells_and_attempts() -> tuple[list[dict], pd.DataFrame]:
    cells = []
    frames = []
    schedules = {
        "eas_phlash_median": base.build_eas_schedule(ROOT),
        "ancient_eurasia_han_introgression": base.build_han_schedule(
            validate_catalog=False
        ),
    }
    for index, original in enumerate(build_cells(_production_plan())):
        cell = {**original, "attempts": 2, "batch_size": 2, "n_batches": 1}
        cells.append(cell)
        frames.append(
            base.simulate_trajectory_batch(
                cell,
                schedules[cell["demography_id"]],
                attempt_start=0,
                attempt_stop=2,
                batch_index=0,
                batch_seed=91_000 + index,
            )
        )
    return cells, pd.concat(frames, ignore_index=True)


def _base_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    result_dir = CAMPAIGN / base.DEFAULT_RESULTS_SUBDIR
    summary = pd.read_csv(result_dir / base.ANALYSIS_OUTPUT_NAMES["summary"], sep="\t")
    survivors = pd.read_csv(
        result_dir / base.ANALYSIS_OUTPUT_NAMES["survivors"], sep="\t"
    )
    return summary, survivors


def _combined_test_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    cells, attempts = _tiny_cells_and_attempts()
    base_summary, base_survivors = _base_tables()
    summary, survivors, _ = combine_analysis_tables(
        attempts, base_summary, base_survivors, cells
    )
    return summary, survivors


def test_default_cells_are_exactly_the_four_new_fixed_attempt_cells():
    cells = build_cells(_production_plan())
    observed = {
        (
            cell["demography_id"],
            cell["selection_coefficient"],
            cell["attempts"],
            cell["batch_size"],
        )
        for cell in cells
    }
    assert observed == {
        (
            "eas_phlash_median",
            0.001,
            DEFAULT_EAS_S001_ATTEMPTS,
            50_000,
        ),
        (
            "eas_phlash_median",
            0.005,
            DEFAULT_EAS_S005_ATTEMPTS,
            50_000,
        ),
        (
            "ancient_eurasia_han_introgression",
            0.001,
            DEFAULT_HAN_S001_ATTEMPTS,
            1_000,
        ),
        (
            "ancient_eurasia_han_introgression",
            0.005,
            DEFAULT_HAN_S005_ATTEMPTS,
            1_000,
        ),
    }
    assert len(cells) == 4
    assert {cell["simulation_class"] for cell in cells} == {"selected"}
    assert sum(cell["n_batches"] for cell in cells) == 80
    assert all(
        "s0p001" in cell["cell_id"] or "s0p005" in cell["cell_id"] for cell in cells
    )


def test_plan_payload_binds_base_without_changing_base_files():
    result_dir = CAMPAIGN / base.DEFAULT_RESULTS_SUBDIR
    paths = sorted(path for path in result_dir.iterdir() if path.is_file())
    before = {path.name: base.sha256_file(path) for path in paths}
    payload = sensitivity._build_plan_payload(_production_plan())
    after = {path.name: base.sha256_file(path) for path in paths}

    assert before == after
    assert payload["scientific_contract"]["new_selection_coefficients"] == [
        0.001,
        0.005,
    ]
    assert payload["scientific_contract"]["combined_selection_coefficients"] == [
        0.0,
        0.001,
        0.005,
        0.01,
    ]
    assert payload["paths"]["work_subdir"] != base.DEFAULT_WORK_SUBDIR
    assert payload["paths"]["results_subdir"] != base.DEFAULT_RESULTS_SUBDIR
    assert len(payload["base_bundle"]["cells"]) == 4
    assert len(payload["cells"]) == 4


def test_bound_base_bundle_rejects_even_manifest_only_tampering():
    plan = _production_plan()
    contract = sensitivity._base_bundle_contract(plan)
    bad = copy.deepcopy(contract)
    bad["summary"]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="changed after sensitivity planning"):
        sensitivity._load_bound_base_tables(plan, bad)


def test_comparisons_are_separate_for_every_positive_coefficient():
    rows = []
    values = {
        0.0: (0.1, 0.2),
        0.001: (0.15,),
        0.005: (0.2,),
        0.01: (0.3,),
    }
    for demography in sensitivity.DEMOGRAPHY_LABELS:
        for coefficient, frequencies in values.items():
            rows.extend(
                {
                    "demography_id": demography,
                    "selection_coefficient": coefficient,
                    "final_population_af": frequency,
                    "population_survived": True,
                }
                for frequency in frequencies
            )
    comparisons = build_comparisons(pd.DataFrame(rows))

    assert len(comparisons) == 6
    assert set(comparisons["selected_selection_coefficient"]) == {
        0.001,
        0.005,
        0.01,
    }
    for demography in sensitivity.DEMOGRAPHY_LABELS:
        cell = comparisons[comparisons["demography_id"] == demography].set_index(
            "selected_selection_coefficient"
        )
        assert cell.loc[0.001, "population_af_roc_auc_survivors"] == pytest.approx(0.5)
        assert cell.loc[0.005, "population_af_roc_auc_survivors"] == pytest.approx(0.75)
        assert cell.loc[0.01, "population_af_roc_auc_survivors"] == pytest.approx(1.0)


def test_combination_has_eight_summary_cells_and_survival_is_only_filter():
    cells, attempts = _tiny_cells_and_attempts()
    base_summary, base_survivors = _base_tables()
    summary, survivors, comparisons = combine_analysis_tables(
        attempts, base_summary, base_survivors, cells
    )

    assert len(summary) == 8
    assert len(comparisons) == 6
    for demography in sensitivity.DEMOGRAPHY_LABELS:
        coefficients = summary.loc[
            summary["demography_id"] == demography, "selection_coefficient"
        ].to_numpy(float)
        assert np.allclose(coefficients, COMBINED_SELECTION_COEFFICIENTS)
    assert survivors["population_survived"].astype(bool).all()
    assert (survivors["final_alt_count"].astype(int) > 0).all()
    assert len(survivors) == int(summary["n_population_survived"].sum())
    assert int(summary["n_attempted"].sum()) == int(
        base_summary["n_attempted"].sum()
    ) + len(attempts)

    tampered = attempts.copy()
    tampered.loc[0, "selection_coefficient"] = 0.01
    with pytest.raises(ValueError, match="unexpected coefficient"):
        combine_analysis_tables(tampered, base_summary, base_survivors, cells)


@pytest.mark.parametrize(
    ("plotter", "stem"),
    (
        (sensitivity._plot_unconditional_fates, "unconditional_fate_fractions"),
        (sensitivity._plot_ecdf, "final_population_af_survivor_ecdf"),
        (
            sensitivity._plot_population_vs_sample,
            "population_vs_sample_af_survivor_ecdf",
        ),
    ),
)
def test_plot_layouts_include_all_four_coefficients_without_rendering(
    tmp_path, plotter, stem
):
    summary, survivors = _combined_test_tables()
    seen: dict[str, object] = {}

    def inspect_only(figure, result_dir, observed_stem):
        assert result_dir == tmp_path
        seen["stem"] = observed_stem
        seen["size"] = tuple(figure.get_size_inches())
        text = [item.get_text() for item in figure.texts]
        for axis in figure.axes:
            seen.setdefault("axis_text", []).extend(
                item.get_text() for item in axis.texts
            )
            text.extend(item.get_text() for item in axis.texts)
            text.extend(item.get_text() for item in axis.get_xticklabels())
        for legend in figure.legends:
            text.extend(item.get_text() for item in legend.get_texts())
        seen["text"] = "\n".join(text)
        return []

    data = summary if plotter is sensitivity._plot_unconditional_fates else survivors
    assert plotter(data, tmp_path, save_pair=inspect_only) == []
    assert seen["stem"] == stem
    assert seen["size"] == pytest.approx(base.FIGURE_SIZE_INCHES)
    for label in ("s=0", "s=0.001", "s=0.005", "s=0.01"):
        assert label in seen["text"]
    if plotter is sensitivity._plot_unconditional_fates:
        attempt_blocks = [
            value for value in seen["axis_text"] if value.count("n=") == 4
        ]
        assert len(attempt_blocks) == 2
        assert all(value.count("\n") == 1 for value in attempt_blocks)
    assert not list(tmp_path.iterdir())
    plt.close("all")


def test_cli_defaults_match_the_frozen_production_attempt_budget():
    args = sensitivity.build_parser().parse_args(["plan"])
    assert args.eas_s001_attempts == 1_500_000
    assert args.eas_s005_attempts == 500_000
    assert args.han_s001_attempts == 20_000
    assert args.han_s005_attempts == 20_000
    assert args.eas_batch_size == 50_000
    assert args.han_batch_size == 1_000
