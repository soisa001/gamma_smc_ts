import importlib.util
from pathlib import Path

import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "analyze_mean_median_strides.py"
SPEC = importlib.util.spec_from_file_location("mean_median_analysis", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def profile(positions, called):
    return pd.DataFrame(
        {
            "position_0based": positions,
            "position_1based": [position + 1 for position in positions],
            "n_pairs": 2,
            "mean_p_tmrca_lt_threshold": [0.1] * len(positions),
            "mean_tmrca_generations": [200.0] * len(positions),
            "mean_p_lt_4500": [0.1] * len(positions),
            "frac_recent_4500": called,
        }
    )


def test_soft_identity_allows_hard_calls_to_differ():
    mean = profile([0, 10_000], [0.0, 0.5])
    median = profile([0, 10_000], [0.5, 1.0])
    MODULE._assert_soft_identity(mean, median)

    median.loc[0, "mean_tmrca_generations"] += 1
    with pytest.raises(RuntimeError, match="soft column"):
        MODULE._assert_soft_identity(mean, median)


def test_shared_stride_metrics_uses_only_common_grid():
    one_kb = profile(
        [0, 1_000, 10_000, 11_000, 20_000],
        [0.0, 0.25, 0.5, 0.75, 1.0],
    )
    ten_kb = profile([0, 10_000, 20_000], [0.0, 0.5, 0.5])
    metrics = MODULE._shared_stride_metrics(
        one_kb,
        ten_kb,
        keys=["position_0based"],
        value="frac_recent_4500",
    )

    assert metrics["shared_rows"] == 3
    assert metrics["identical_rows"] == 2
    assert metrics["maximum_absolute_difference"] == 0.5
