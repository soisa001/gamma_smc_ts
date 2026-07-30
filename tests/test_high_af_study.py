import pandas as pd
import pytest

from gamma_smc_aou.container_study import (
    _assert_matched_soft_profiles,
    _plot_decoded_center_calibration,
    _write_recent_call_comparison,
)
from gamma_smc_aou.high_af_study import _check_null_compatibility


def _design():
    return {
        "ancestral_population_size": 10_000,
        "present_population_size": 20_000,
        "size_change_generations_ago": 100,
        "sample_diploids": 2_000,
        "sequence_length": 10_000_000,
        "variant_age_generations": 180,
        "generation_time_years": 25,
        "mutation_rate": 1.25e-8,
        "recombination_rate": 1e-8,
    }


def test_reused_truth_null_must_match_selected_design():
    design = _design()
    metrics = {
        **{key: value for key, value in design.items() if "population_size" not in key and key != "size_change_generations_ago"},
        "demography": {
            "ancestral_population_size": 10_000,
            "present_population_size": 20_000,
            "size_change_generations_ago": 100,
        },
    }
    _check_null_compatibility(metrics, design)
    metrics["mutation_rate"] = 2e-8
    with pytest.raises(ValueError, match="mutation_rate"):
        _check_null_compatibility(metrics, design)


def test_decoded_center_plot_reports_monte_carlo_distribution(tmp_path):
    neutral = pd.DataFrame({
        "replicate": list(range(100)),
        "position_0based": [5_000_000] * 100,
        "mean_p_tmrca_lt_threshold": [index / 1000 for index in range(100)],
    })
    output = tmp_path / "decoded_center.png"
    _plot_decoded_center_calibration(
        neutral,
        center_position=5_000_000,
        observed_value=0.075,
        pvalue=26 / 101,
        output_path=output,
    )
    assert output.exists()
    assert output.stat().st_size > 10_000


def test_matched_recent_call_comparison_reports_each_statistic(tmp_path):
    positions = [0.0, 10_000.0, 20_000.0]
    base = {
        "position_0based": positions,
        "position_1based": [1.0, 10_001.0, 20_001.0],
        "n_pairs": [2_000, 2_000, 2_000],
        "mean_p_tmrca_lt_threshold": [0.1, 0.2, 0.1],
        "mean_tmrca_generations": [500.0, 400.0, 500.0],
        "mean_p_lt_4500": [0.1, 0.2, 0.1],
    }
    mean_observed = pd.DataFrame({
        **base,
        "n_recent_4500": [100, 400, 100],
        "frac_recent_4500": [0.05, 0.20, 0.05],
    })
    median_observed = pd.DataFrame({
        **base,
        "n_recent_4500": [200, 600, 200],
        "frac_recent_4500": [0.10, 0.30, 0.10],
    })
    _assert_matched_soft_profiles(mean_observed, median_observed)

    mean_null = []
    median_null = []
    for replicate in range(20):
        mean_null.append(pd.DataFrame({
            **base,
            "mean_p_tmrca_lt_threshold": [0.05, 0.05, 0.05],
            "mean_p_lt_4500": [0.05, 0.05, 0.05],
            "n_recent_4500": [20, 20, 20],
            "frac_recent_4500": [0.01, 0.01, 0.01],
            "replicate": replicate,
        }))
        median_null.append(pd.DataFrame({
            **base,
            "mean_p_tmrca_lt_threshold": [0.05, 0.05, 0.05],
            "mean_p_lt_4500": [0.05, 0.05, 0.05],
            "n_recent_4500": [40, 40, 40],
            "frac_recent_4500": [0.02, 0.02, 0.02],
            "replicate": replicate,
        }))
    truth = pd.DataFrame({
        "position_0based": positions,
        "truth_fraction_recent": [0.05, 0.25, 0.05],
    })
    metrics = _write_recent_call_comparison(
        {"mean": mean_observed, "median": median_observed},
        {
            "mean": pd.concat(mean_null, ignore_index=True),
            "median": pd.concat(median_null, ignore_index=True),
        },
        truth,
        threshold_years=4_500,
        center=10_000,
        sequence_length=30_000,
        stride=10_000,
        output_dir=tmp_path,
    )
    assert set(metrics) == {
        "soft_probability",
        "posterior_mean_call",
        "posterior_median_call",
        "truth",
    }
    assert metrics["posterior_mean_call"]["center_selected"] == 0.20
    assert metrics["posterior_median_call"]["center_selected"] == 0.30
    assert metrics["soft_probability"]["center_monte_carlo_p_upper"] == 1 / 21
    assert (tmp_path / "posterior_mean_median_rule_comparison.png").exists()
    assert (
        tmp_path / "posterior_summary_rule_comparison_metrics.json"
    ).exists()
