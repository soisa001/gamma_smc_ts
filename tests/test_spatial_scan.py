import numpy as np
import pandas as pd

from gamma_smc_aou.spatial_scan import (
    _plot_null_spatial_calibration,
    _plot_observed_profile,
    calibrate_spatial_windows,
    significant_regions,
    window_centers,
)


def test_window_grid_contains_selected_center():
    positions = window_centers(10_000_000, 20_000, 5_000_000)
    assert len(positions) == 501
    assert 5_000_000 in positions
    assert positions[0] == 0
    assert positions[-1] == 9_999_999


def test_spatial_pvalues_and_breakpoints_are_replicate_based(tmp_path):
    positions = np.asarray([0, 20_000, 40_000, 100_000], dtype=float)
    observed = pd.DataFrame({
        "position_0based": positions,
        "n_pairs": 2_000,
        "mean_p_tmrca_lt_threshold": [0.5, 0.5, 0.0, 0.5],
    })
    neutral = pd.concat(
        [
            pd.DataFrame({
                "position_0based": positions,
                "n_pairs": 2_000,
                "mean_p_tmrca_lt_threshold": [0.1, 0.1, 0.1, 0.1],
                "replicate": replicate,
            })
            for replicate in range(20)
        ],
        ignore_index=True,
    )
    scan = calibrate_spatial_windows(observed, neutral)
    assert scan.loc[0, "p_upper"] == 1 / 21
    assert scan.loc[2, "p_upper"] == 1.0
    regions = significant_regions(
        scan, sequence_length=120_000, window_size=20_000
    )
    assert regions[["start_0based", "end_0based_exclusive"]].values.tolist() == [
        [0, 30_000],
        [90_000, 110_000],
    ]

    observed_plot = tmp_path / "observed.png"
    null_plot = tmp_path / "null.png"
    _plot_observed_profile(
        scan,
        center=100_000,
        sequence_length=120_000,
        zoom_half_width=20_000,
        threshold_years=4_500,
        window_size=20_000,
        output_path=observed_plot,
    )
    _plot_null_spatial_calibration(
        scan,
        regions,
        center=100_000,
        sequence_length=120_000,
        zoom_half_width=20_000,
        threshold_years=4_500,
        window_size=20_000,
        output_path=null_plot,
    )
    assert observed_plot.exists()
    assert null_plot.exists()


def test_spatial_calibration_drops_positions_missing_from_any_null():
    observed = pd.DataFrame({
        "position_0based": [0.0, 1000.0],
        "n_pairs": [2, 2],
        "mean_p_tmrca_lt_threshold": [0.2, 0.3],
    })
    neutral = pd.DataFrame({
        "position_0based": [0.0, 1000.0, 0.0],
        "n_pairs": [2, 2, 2],
        "mean_p_tmrca_lt_threshold": [0.1, 0.1, 0.1],
        "replicate": [0, 0, 1],
    })
    scan = calibrate_spatial_windows(observed, neutral)
    assert scan["position_0based"].tolist() == [0.0]
