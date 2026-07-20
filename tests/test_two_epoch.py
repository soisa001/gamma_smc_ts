import os
import shutil

import numpy as np
import pandas as pd
import pytest

from gamma_smc_aou.selection import run_slim_recent_sweep
from gamma_smc_aou.two_epoch import (
    _plot_demography_and_variant,
    _plot_null_calibration,
    two_epoch_recent_probability,
)


def test_two_epoch_probability_integrates_both_population_sizes():
    observed = two_epoch_recent_probability(
        ancestral_population_size=10_000,
        present_population_size=20_000,
        size_change_generations_ago=100,
        threshold_generations=180,
    )
    expected = 1 - np.exp(-(100 / 40_000 + 80 / 20_000))
    assert observed == pytest.approx(expected)


def test_two_epoch_plots_show_timing_and_upper_tail_pvalue(tmp_path):
    demography_path = tmp_path / "demography.png"
    _plot_demography_and_variant(
        demography_path,
        ancestral_population_size=10_000,
        present_population_size=20_000,
        size_change_generations_ago=100,
        variant_age_generations=180,
        generation_time_years=25,
        selection_coefficient=0.05,
    )
    assert demography_path.exists()

    neutral = pd.DataFrame({"center_fraction_recent": [0.005, 0.006, 0.007]})
    calibration_path = tmp_path / "calibration.png"
    exceedances, pvalue = _plot_null_calibration(
        neutral,
        0.05,
        theoretical_fraction_recent=0.0065,
        threshold_generations=180,
        generation_time_years=25,
        selection_coefficient=0.05,
        output_path=calibration_path,
    )
    assert exceedances == 0
    assert pvalue == 0.25
    assert calibration_path.exists()


@pytest.mark.skipif(
    not (os.environ.get("SLIM_BIN") or shutil.which("slim")),
    reason="SLiM executable not available",
)
def test_recent_sweep_changes_population_size_at_requested_time(tmp_path):
    _, run = run_slim_recent_sweep(
        tmp_path / "two_epoch.trees",
        population_size=100,
        present_population_size=200,
        size_change_generations_ago=10,
        sample_diploids=50,
        sequence_length=200_000,
        selection_coefficient=0.0,
        age_generations=30,
        recombination_rate=2e-7,
        seed=12345,
    )
    assert run["ancestral_population_size"] == 100
    assert run["present_population_size"] == 200
    assert "SIZE_CHANGE population_size=200 tick=21" in run["slim_stdout"]
