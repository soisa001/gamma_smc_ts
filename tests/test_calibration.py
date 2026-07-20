import numpy as np
import pandas as pd

from gamma_smc_aou.calibration import bh_fdr, calibrate_sites, monte_carlo_pvalue, randomized_rank_pvalue


def frame(positions, values):
    return pd.DataFrame({
        "position_0based": positions,
        "position_1based": np.asarray(positions) + 1,
        "mean_p_tmrca_lt_threshold": values,
    })


def test_monte_carlo_upper_tail_and_plus_one():
    assert monte_carlo_pvalue(0.8, [0.1, 0.2, 0.9]) == 0.5
    assert monte_carlo_pvalue(1.0, np.zeros(999)) == 0.001
    assert randomized_rank_pvalue(0.5, [0.5, 0.5, 0.2], 0.5) == 0.375


def test_nearest_uses_one_value_per_replicate():
    observed = frame([50], [0.8])
    simulations = [frame([49, 90], [0.7, 1.0]), frame([10, 52], [0.2, 0.9])]
    result, null = calibrate_sites(
        observed, simulations, observed_length=100, simulation_lengths=[100, 100], match="nearest"
    )
    np.testing.assert_allclose(null, [[0.7, 0.9]])
    assert result.loc[0, "n_null"] == 2
    assert result.loc[0, "mc_p_upper"] == 2 / 3


def test_region_max_is_conservative_and_bh_monotone():
    observed = frame([20, 80], [0.8, 0.6])
    simulations = [frame([10, 90], [0.2, 0.9]), frame([40], [0.7])]
    result, null = calibrate_sites(
        observed, simulations, observed_length=100, match="region_max"
    )
    np.testing.assert_allclose(null, [[0.9, 0.7], [0.9, 0.7]])
    assert np.all((result["bh_q"] >= result["mc_p_upper"]) | np.isclose(result["bh_q"], result["mc_p_upper"]))
    np.testing.assert_allclose(bh_fdr([0.01, 0.04, 0.03]), [0.03, 0.04, 0.04])
