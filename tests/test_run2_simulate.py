"""Tree-truth extraction tests for the run2 study.

A small msprime simulation stands in for a SLiM replicate: the statistic code
only ever sees a tree sequence, so it can be exercised without SLiM.
"""

from __future__ import annotations

import msprime
import numpy as np
import pytest

from gamma_smc_aou import run2_simulate as simulate
from gamma_smc_aou.run2_config import FOCAL_POSITION_BP, TMRCA_THRESHOLDS_YEARS


@pytest.fixture(scope="module")
def small_ts():
    return msprime.sim_ancestry(
        samples=20,
        population_size=1_000,
        sequence_length=200_000,
        recombination_rate=1e-8,
        ploidy=2,
        random_seed=42,
    )


# ---------------------------------------------------------------------------
# Pairs and positions
# ---------------------------------------------------------------------------


def test_sampled_diploid_pairs_partition_the_sample_nodes(small_ts):
    pairs = simulate.sampled_diploid_pairs(small_ts)
    assert pairs.shape == (20, 2)
    flat = pairs.ravel()
    assert sorted(flat.tolist()) == sorted(int(n) for n in small_ts.samples())


def test_profile_positions_always_include_the_focal_base():
    positions = simulate.profile_positions()
    assert float(FOCAL_POSITION_BP) in set(positions.tolist())
    assert np.all(np.diff(positions) > 0)
    assert positions[0] == 0.0
    assert positions[-1] < 10_000_000


def test_profile_positions_handle_a_focal_base_off_the_stride_grid():
    positions = simulate.profile_positions(
        sequence_length=1_000, stride=100, focal_position=555
    )
    assert 555.0 in set(positions.tolist())
    assert len(positions) == 11


# ---------------------------------------------------------------------------
# TMRCA matrix
# ---------------------------------------------------------------------------


def test_pairwise_tmrca_matrix_matches_direct_lookup(small_ts):
    pairs = simulate.sampled_diploid_pairs(small_ts)
    positions = np.array([0.0, 50_000.0, 123_456.0, 199_999.0])
    matrix = simulate.pairwise_tmrca_matrix(small_ts, pairs, positions)
    assert matrix.shape == (4, 20)
    assert np.all(matrix > 0)
    for row, position in enumerate(positions):
        tree = small_ts.at(position)
        expected = [tree.tmrca(int(a), int(b)) for a, b in pairs]
        assert matrix[row] == pytest.approx(expected)


def test_pairwise_tmrca_matrix_covers_the_whole_stride_grid(small_ts):
    pairs = simulate.sampled_diploid_pairs(small_ts)
    positions = simulate.profile_positions(
        sequence_length=200_000, stride=10_000, focal_position=100_000
    )
    matrix = simulate.pairwise_tmrca_matrix(small_ts, pairs, positions)
    assert matrix.shape == (len(positions), 20)
    assert np.all(np.isfinite(matrix))


# ---------------------------------------------------------------------------
# Scores
# ---------------------------------------------------------------------------


def _scores_for(counts, n_pairs=20):
    positions = np.array([0.0, 100.0])
    # Two pairs' worth of contrived TMRCAs: everything coalesces at 10 generations.
    tmrca = np.full((2, n_pairs), 10.0)
    return simulate.focal_scores(
        tmrca, positions, counts, focal_position=100, thresholds_years=(1000.0,)
    )


def test_focal_scores_cover_every_genotype_class():
    frame = _scores_for(np.array([2] * 10 + [1] * 5 + [0] * 5, dtype=np.int8))
    assert set(frame["genotype_class"]) == set(simulate.GENOTYPE_CLASSES)
    by_class = frame.set_index("genotype_class")
    assert by_class.loc["overall", "n_pairs"] == 20
    assert by_class.loc["hom_carrier", "n_pairs"] == 10
    assert by_class.loc["het", "n_pairs"] == 5
    assert by_class.loc["hom_noncarrier", "n_pairs"] == 5
    # 1000 years is 40 generations, so every pair coalesces below it.
    assert by_class.loc["overall", "p_tmrca_lt_threshold"] == pytest.approx(1.0)


def test_focal_scores_treat_a_missing_focal_site_as_all_noncarrier():
    frame = _scores_for(None)
    by_class = frame.set_index("genotype_class")
    assert by_class.loc["overall", "n_pairs"] == 20
    assert by_class.loc["hom_noncarrier", "n_pairs"] == 20
    assert by_class.loc["hom_carrier", "n_pairs"] == 0
    assert np.isnan(by_class.loc["hom_carrier", "p_tmrca_lt_threshold"])
    # The overall class -- the one the test uses -- is always defined.
    assert np.isfinite(by_class.loc["overall", "p_tmrca_lt_threshold"])


def test_focal_scores_require_the_focal_position_on_the_grid():
    with pytest.raises(ValueError, match="focal position is missing"):
        simulate.focal_scores(
            np.full((1, 4), 10.0), np.array([0.0]), None, focal_position=999
        )


def test_focal_scores_use_a_strict_inequality():
    positions = np.array([0.0])
    tmrca = np.full((1, 4), 40.0)  # exactly 1000 years at 25 y/generation
    frame = simulate.focal_scores(
        tmrca, positions, None, focal_position=0, thresholds_years=(1000.0,)
    )
    overall = frame.set_index("genotype_class").loc["overall"]
    assert overall["p_tmrca_lt_threshold"] == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# Spatial profile
# ---------------------------------------------------------------------------


def test_spatial_profile_has_one_column_per_threshold(small_ts):
    pairs = simulate.sampled_diploid_pairs(small_ts)
    positions = simulate.profile_positions(
        sequence_length=200_000, stride=20_000, focal_position=100_000
    )
    tmrca = simulate.pairwise_tmrca_matrix(small_ts, pairs, positions)
    frame = simulate.spatial_profile(tmrca, positions)
    for threshold in TMRCA_THRESHOLDS_YEARS:
        column = f"p_lt_{int(threshold)}y"
        assert column in frame.columns
        assert frame[column].between(0.0, 1.0).all()
    assert len(frame) == len(positions)
    assert (frame["mean_tmrca_generations"] > 0).all()


def test_spatial_profile_probabilities_increase_with_the_threshold(small_ts):
    pairs = simulate.sampled_diploid_pairs(small_ts)
    positions = np.array([0.0, 100_000.0])
    tmrca = simulate.pairwise_tmrca_matrix(small_ts, pairs, positions)
    frame = simulate.spatial_profile(tmrca, positions)
    columns = [f"p_lt_{int(t)}y" for t in TMRCA_THRESHOLDS_YEARS]
    values = frame[columns].to_numpy()
    assert np.all(np.diff(values, axis=1) >= 0)


# ---------------------------------------------------------------------------
# Task construction
# ---------------------------------------------------------------------------


def test_build_tasks_is_deterministic_and_indexable(tmp_path):
    from gamma_smc_aou.run2_config import CHB_ARM

    tasks = simulate.build_tasks(
        CHB_ARM,
        "selected",
        study_root=tmp_path,
        repo_root=tmp_path,
        slim_path="slim",
        replicates=3,
    )
    assert [task.replicate_index for task in tasks] == [0, 1, 2]
    assert [task.seed for task in tasks] == [CHB_ARM.seed("selected", i) for i in range(3)]
    assert tasks[0].unit_id.endswith("selected__rep000")

    subset = simulate.build_tasks(
        CHB_ARM,
        "selected",
        study_root=tmp_path,
        repo_root=tmp_path,
        slim_path="slim",
        replicates=100,
        indices=[7],
    )
    assert len(subset) == 1
    assert subset[0].seed == CHB_ARM.seed("selected", 7)


def test_build_tasks_rejects_an_unknown_mode(tmp_path):
    from gamma_smc_aou.run2_config import CHB_ARM

    with pytest.raises(ValueError):
        simulate.build_tasks(
            CHB_ARM,
            "conditioned",
            study_root=tmp_path,
            repo_root=tmp_path,
            slim_path="slim",
            replicates=1,
        )
