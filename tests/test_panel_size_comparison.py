import numpy as np
import pandas as pd

from gamma_smc_aou.panel_size_comparison import (
    empirical_upper_p,
    evaluate,
    fold_roles,
    make_folds,
    target_threshold,
)


def test_region_rank_includes_ties_and_cannot_return_zero():
    np.testing.assert_array_equal(
        empirical_upper_p(np.array([2, 4, 4]), np.array([1, 4, 5])), [1, 0.75, 0.25]
    )
    assert empirical_upper_p(np.ones(100), np.ones(1))[0] == 1


def test_selected_threshold_retains_at_least_requested_fraction_with_ties():
    selected = np.array([1, 2, 2, 3, 3, 3, 4, 5, 6, 7])
    threshold = target_threshold(selected, 0.7)
    assert threshold == 3 and np.mean(selected >= threshold) == 0.7


def test_held_out_values_do_not_change_their_training_thresholds():
    meta = pd.DataFrame(
        [
            dict(
                task_id=f"{mode}/{i}",
                mode=mode,
                s=s,
                replicate=i,
                sample_af=0.5,
                fixed=False,
            )
            for mode, s, n in [("neutral", 0.0, 20), ("selected", 0.005, 10)]
            for i in range(n)
        ]
    )
    seed = 20380101
    folds = make_folds(meta, seed, 5)
    test, fit, calibration, selected_train = fold_roles(meta, folds, 0)
    assert set(test).isdisjoint(set(fit) | set(calibration) | set(selected_train))
    assert set(fit).isdisjoint(calibration)
    assert len(set(test) | set(fit) | set(calibration) | set(selected_train)) == len(
        meta
    )
    values = np.random.default_rng(seed).uniform(size=(len(meta), 4, 1))
    _, before = evaluate(values, meta, [10000], seed, 400, "decoded")
    values[test] = 1e6
    _, after = evaluate(values, meta, [10000], seed, 400, "decoded")
    assert [r for r in before if r["fold"] == 0] == [r for r in after if r["fold"] == 0]
