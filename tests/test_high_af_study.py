import json

import pandas as pd
import pytest

from gamma_smc_aou.container_study import _plot_decoded_center_calibration
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
