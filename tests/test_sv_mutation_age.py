from dataclasses import dataclass

import numpy as np
import pandas as pd
import pytest
import tskit
from gamma_smc_aou.sv_mutation_age import (
    cached,
    complete,
    composite_curves,
    exact_bounds,
    fit_tempered_composite,
    introgressed_samples,
    mask_intervals,
    pair_classes,
    split_name,
)


def example_tree():
    tables = tskit.TableCollection(10)
    for _ in range(3):
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
    tables.nodes.add_row(time=10)
    tables.nodes.add_row(time=150)
    for parent, child in [(3, 0), (3, 1), (4, 3), (4, 2)]:
        tables.edges.add_row(0, 10, parent, child)
    tables.sort()
    return tables.tree_sequence().first()


def test_branch_spanning_introgression_does_not_identify_mutation_age():
    tree = example_tree()
    low, high = exact_bounds(tree, [0, 1, 2], np.array([True, True, False]))
    assert (low, high) == (10, 150)
    assert low < 20 < 50 < 80 < high  # both new and old SV fit identical observed tree


def test_migration_mapping_uses_interval_and_descendants():
    @dataclass
    class Migration:
        left: int = 0
        right: int = 5
        node: int = 3
        source: int = 0
        dest: int = 1
        time: int = 50

    tree = example_tree()
    assert introgressed_samples(tree, [0, 1, 2], [Migration()], 2, 50).tolist() == [
        True,
        True,
        False,
    ]
    assert not introgressed_samples(tree, [0, 1, 2], [Migration()], 5, 50).any()


def test_composite_support_converges_to_true_branch_not_carrier_crown():
    grid = np.array([5.0, 20.0, 80.0, 200.0])
    # Concentrated pair posteriors: carrier MRCA=10, opposite-state MRCA=150.
    curves = composite_curves(
        np.array([1e5, 1e5]),
        np.array([1e4, 1e5 / 150]),
        {"CC": [0], "CN": [1], "ICN": []},
        grid,
        1,
    )
    score = curves["CC"].sum(0) + curves["CN"].sum(0)
    assert score[1] == pytest.approx(score[2])
    assert score[0] < score[1] - 100 and score[3] < score[2] - 100


def test_pair_sampling_never_uses_wrong_allele_class():
    carrier = np.array([True, True, True, False, False, False])
    ancestry = np.array([True, True, True, True, False, False])
    pairs = pair_classes(carrier, ancestry, 2, np.random.default_rng(1729))
    assert all(carrier[a] and carrier[b] for a, b in pairs["CC"])
    assert all(carrier[a] != carrier[b] for a, b in pairs["CN"])
    assert all(
        ancestry[a] and ancestry[b] and carrier[a] != carrier[b]
        for a, b in pairs["ICN"]
    )
    assert all(len(p) <= 2 for p in pairs.values())


def test_ancestry_composite_does_not_count_the_same_pair_twice():
    curves = composite_curves(
        np.ones(2),
        np.ones(2),
        {"CC": [0], "CN": [1], "ICN": [1]},
        np.array([1.0, 2.0]),
        1.0,
    )
    assert curves["ICN"].shape == (0, 2)


def test_masks_merge_overlaps_and_preserve_half_open_coordinates():
    assert mask_intervals(20, [5, 7], 2) == [(0, 3), (10, 20)]
    assert mask_intervals(20, [0, 19], 0) == [(1, 19)]


def test_cache_detects_corruption_without_commit_dependency(tmp_path):
    (tmp_path / "data.txt").write_text("complete")
    spec = {"sample_manifest": "abc"}
    complete(tmp_path, "receipt", spec, ["data.txt"])
    assert cached(tmp_path, "receipt", spec)
    (tmp_path / "data.txt").write_text("corrupt")
    with pytest.raises(ValueError, match="Corrupt"):
        cached(tmp_path, "receipt", spec)


def test_region_split_keeps_all_loci_and_mask_conditions_together():
    assert [split_name(i) for i in range(8)] == [
        "train",
        "train",
        "calibration",
        "test",
    ] * 2


def test_tempered_fit_does_not_use_test_mutation_ages():
    grid = np.geomspace(100, 1000000, 100)
    curve = -0.5 * ((np.log(grid) - np.log(20000)) / 0.2) ** 2
    frame = pd.DataFrame(
        {
            "id": [str(i) for i in range(50)],
            "replicate": np.arange(50),
            "mask_radius": 0,
            "split": ["train"] * 30 + ["calibration"] * 10 + ["test"] * 10,
            "age_years": 20000.0,
            "log_compatibility": [curve] * 50,
        }
    )
    cfg = {"pulse_years": 50000, "interval_alpha": 0.1}
    indices = (np.arange(30), np.arange(30, 40), np.arange(40, 50))
    first, first_fit = fit_tempered_composite(
        frame, grid, np.gradient(grid), cfg, *indices
    )
    frame.loc[40:, "age_years"] = 1000000.0
    second, second_fit = fit_tempered_composite(
        frame, grid, np.gradient(grid), cfg, *indices
    )
    assert first_fit == second_fit
    assert first == second
