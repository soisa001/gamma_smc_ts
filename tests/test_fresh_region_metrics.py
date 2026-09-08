import numpy as np

from gamma_smc_aou.fresh_region_metrics import (
    aligned_pvalues,
    runs,
    score_peaks,
    score_strides,
)


def test_upper_tail_with_ties_and_neutral_leave_one_out():
    neutral = np.array([0, 2, 2, 5]).reshape(4, 1, 1)
    selected = np.array([1, 2, 6]).reshape(3, 1, 1)
    null_p, selected_p = aligned_pvalues(neutral, selected)
    np.testing.assert_allclose(null_p.ravel(), [1, 0.75, 0.75, 0.25])
    np.testing.assert_allclose(selected_p.ravel(), [0.8, 0.8, 0.2])


def test_neutral_calibration_is_specific_to_position():
    neutral = np.array([[1, 8], [2, 9], [3, 10]]).reshape(3, 2, 1)
    _, selected_p = aligned_pvalues(neutral, np.array([4, 4]).reshape(1, 2, 1))
    np.testing.assert_allclose(selected_p.ravel(), [0.25, 1])


def test_stride_counts_and_no_call_convention():
    target = np.array([False, True, True, True, False])
    result = score_strides(np.array([True, True, False, False, False]), target)
    assert (result["tp"], result["fp"], result["fn"]) == (1, 1, 2)
    assert result["fdp"] == 0.5 and result["f1"] == 0.4
    no_calls = score_strides(np.zeros(5, dtype=bool), target)
    assert no_calls["fdp"] == no_calls["f1"] == 0


def test_peaks_handle_edges_width_and_one_to_one_matching():
    mask = np.array([1, 1, 0, 1, 0, 1, 1], dtype=bool)
    assert runs(mask) == [(0, 2), (3, 4), (5, 7)]
    assert runs(mask, 2) == [(0, 2), (5, 7)]
    assert runs(np.zeros(7, dtype=bool)) == []
    # Two target-overlapping fragments can match only one true event.
    result = score_peaks(runs(mask), np.array([0, 1, 1, 1, 0, 0, 0], dtype=bool))
    assert (result["tp"], result["fp"], result["fn"]) == (1, 2, 0)
    assert result["f1"] == 0.5
