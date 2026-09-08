import numpy as np

from gamma_smc_aou.paired_peak_grid import peak_groups
from gamma_smc_aou.fresh_region_metrics import aligned_pvalues


def test_inclusive_point001_rank_boundary_has_valid_nonzero_calls():
    neutral = np.arange(1000).reshape(1000, 1, 1)
    selected = np.array([1000, 999, 998]).reshape(3, 1, 1)
    null_p, selected_p = aligned_pvalues(neutral, selected)
    assert np.count_nonzero(null_p <= 0.001) == 1
    np.testing.assert_equal((selected_p <= 0.001).ravel(), [True, False, False])
    assert np.all(null_p > 0) and np.all(selected_p > 0)


def test_gap_merging_counts_evidence_not_total_span():
    mask = np.array([1, 1, 0, 1, 0, 0, 1], dtype=bool)
    values = np.array([1, 3, 0, 2, 0, 0, 9])
    assert peak_groups(mask, values, 0) == [(2, 1, 0, 1), (1, 3, 3, 3), (1, 6, 6, 6)]
    assert peak_groups(mask, values, 1) == [(3, 1, 0, 3), (1, 6, 6, 6)]
