import numpy as np
import pandas as pd

from gamma_smc_aou.partial_sweep_cv import (
    make_folds,
    normalize,
    select_threshold,
    summit_candidates,
    threshold_scores,
)


def test_folds_hold_out_whole_replicates_and_are_deterministic():
    meta = pd.DataFrame(
        {
            "mode": ["neutral"] * 10 + ["selected"] * 10,
            "s": [0.0] * 10 + [0.005] * 10,
            "replicate": list(range(10)) * 2,
        }
    )
    folds = make_folds(meta, 20380101, 5)
    np.testing.assert_equal(folds, make_folds(meta, 20380101, 5))
    assert all(np.count_nonzero(folds[:10] == fold) == 2 for fold in range(5))
    assert all(np.count_nonzero(folds[10:] == fold) == 2 for fold in range(5))


def test_normalization_never_uses_test_values():
    x = np.array([[1.0, 3.0], [3.0, 5.0], [20.0, 30.0]])
    first = normalize(x, [0, 1])
    x[2] *= 100
    second = normalize(x, [0, 1])
    np.testing.assert_equal(first[:2], second[:2])


def test_wide_peak_cannot_become_a_true_summit_by_overlap():
    peaks = np.array([[5.0, 2, 0, 10]])
    score = threshold_scores(peaks, 4, 0, 4, 6)
    assert score["detected"] == 0 and score["fp"] == 1
    assert score["overlap_detected"] == 1
    empty = threshold_scores(peaks, 6, 0, 4, 6)
    assert empty["calls"] == empty["fdp"] == empty["f1"] == 0


def test_boundaries_and_peak_separation():
    peaks = summit_candidates(np.array([[4.0, 0, 2, 0, 3.0]]), 3)[0]
    np.testing.assert_equal(peaks[:, 1], [0, 4])


def test_threshold_choice_does_not_look_at_test_peak():
    peaks = [
        np.array([[5.0, 5, 5, 5]]),
        np.array([[4.0, 5, 5, 5]]),
        np.array([[3.0, 1, 1, 1]]),
        np.array([[100.0, 5, 5, 5]]),
    ]
    chosen = select_threshold(peaks, [0, 1], [2], 0.7, [0, 1], 4, 6)
    peaks[3] = np.array([[0.1, 1, 1, 1]])
    assert chosen == select_threshold(peaks, [0, 1], [2], 0.7, [0, 1], 4, 6)
    assert chosen["training_power"] == 1 and chosen["training_pooled_fdr"] == 0


def test_neutral_error_objective_takes_priority_over_localization_fdr():
    peaks = [
        np.array([[100.0, 0, 0, 0], [strength, 5, 5, 5]])
        for strength in (10.0, 9.0, 8.0, 7.0)
    ]
    peaks.extend([np.array([[9.0, 0, 0, 0]]), np.array([[8.0, 0, 0, 0]])])
    chosen = select_threshold(peaks, [0, 1, 2, 3], [4, 5], 0.5, [0, 1], 4, 6)
    assert chosen["threshold"] == 9
    assert chosen["training_power"] == 0.5
    assert chosen["training_neutral_any"] == 0.5


def test_region_sensitivity_does_not_require_call_near_selected_site():
    peaks = [
        np.array([[5.0, 0, 0, 0]]),
        np.array([[4.0, 9, 9, 9]]),
        np.array([[3.0, 0, 0, 0]]),
    ]
    assert select_threshold(peaks, [0, 1], [2], 0.7, [1], 4, 6) is None
    chosen = select_threshold(
        peaks, [0, 1], [2], 0.7, [1], 4, 6, endpoint="region_any_call"
    )
    assert chosen["training_power"] == 1
    assert chosen["training_neutral_any"] == 0
