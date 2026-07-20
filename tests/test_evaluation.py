import pandas as pd

from gamma_smc_aou.evaluation import align_nearest, error_metrics


def test_decoder_error_metrics_and_nearest_alignment():
    truth = pd.DataFrame({
        "position_0based": [10, 50],
        "mean_p_tmrca_lt_threshold": [0.2, 0.8],
        "mean_tmrca_generations": [1000, 200],
    })
    decoded = pd.DataFrame({
        "position_0based": [12, 48],
        "mean_p_tmrca_lt_threshold": [0.3, 0.7],
        "mean_tmrca_generations": [900, 250],
    })
    aligned = align_nearest(truth, decoded)
    metrics = error_metrics(aligned)
    assert metrics["n"] == 2
    assert abs(metrics["mae_recent_probability"] - 0.1) < 1e-12
    assert abs(metrics["rmse_recent_probability"] - 0.1) < 1e-12
    assert abs(metrics["bias_recent_probability"]) < 1e-12
