import os
import shutil

import msprime
import numpy as np
import pandas as pd
import pytest

from gamma_smc_aou.selection import (
    validate_recent_sweep_grid,
    validate_slim_hard_sweep,
    within_individual_tmrca_grid,
)


def test_tmrca_grid_uses_one_pair_per_diploid():
    ts = msprime.sim_ancestry(
        samples=[msprime.SampleSet(6, ploidy=2)], population_size=100,
        sequence_length=1_000, recombination_rate=1e-7,
        model=msprime.StandardCoalescent(), random_seed=12,
    )
    result = within_individual_tmrca_grid(ts, np.asarray([100, 500, 900]), 150)
    assert len(result) == 3
    assert result["n_pairs"].eq(6).all()
    assert result["mean_p_tmrca_lt_threshold"].between(0, 1).all()


@pytest.mark.skipif(
    not (os.environ.get("SLIM_BIN") or shutil.which("slim")),
    reason="SLiM executable not available",
)
def test_slim_hard_sweep_has_recent_local_genealogy(tmp_path):
    metrics = validate_slim_hard_sweep(
        tmp_path,
        population_size=100,
        sequence_length=50_000,
        selection_coefficient=0.75,
        # Elevated recombination localizes the sweep in this short, fast fixture.
        recombination_rate=2e-6,
        threshold_generations=100,
        neutral_replicates=39,
        seed=97531,
    )
    assert metrics["center_mean_tmrca_generations"] < metrics["flank_mean_tmrca_generations"]
    assert metrics["center_to_flank_tmrca_ratio"] < 0.5
    assert metrics["center_fraction_recent"] > metrics["neutral_mean_fraction_recent"]
    assert metrics["mc_p_upper"] <= 0.05
    assert (tmp_path / "hard_sweep_validation.png").exists()


@pytest.mark.skipif(
    not (os.environ.get("SLIM_BIN") or shutil.which("slim")),
    reason="SLiM executable not available",
)
def test_recent_sweep_grid_records_full_region_and_mc_test(tmp_path):
    metrics = validate_recent_sweep_grid(
        tmp_path,
        population_size=100,
        sample_diploids=50,
        sequence_length=200_000,
        selection_coefficients=(0.0, 0.5),
        age_generations=30,
        mutation_rate=1.25e-8,
        recombination_rate=2e-7,
        neutral_replicates=9,
        selected_replicates=1,
        workers=2,
        seed=314159,
    )
    assert metrics["threshold_years"] == 750
    assert metrics["neutral_replicates"] == 9
    assert (tmp_path / "recent_sweep_10mb_validation.png").exists()
    assert (tmp_path / "replicate_statistics.tsv").exists()
    assert (tmp_path / "truth_profiles.tsv").exists()
    assert (tmp_path / "statistical_summary.json").exists()
    assert (tmp_path / "power_summary.tsv").exists()
    rows = pd.read_csv(tmp_path / "replicate_statistics.tsv", sep="\t")
    neutral = rows[rows["selection_coefficient"] == 0]
    assert len(neutral) == 9
    assert neutral["focal_allele_outcome"].eq("lost").any()
    assert neutral["center_fraction_recent"].notna().all()
