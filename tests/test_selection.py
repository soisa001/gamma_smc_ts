import os
import shutil

import msprime
import numpy as np
import pandas as pd
import pytest

from gamma_smc_aou.selection import (
    retained_sweep_calibration_table,
    summarize_carrier_tmrca,
    validate_recent_sweep_grid,
    validate_slim_hard_sweep,
    within_individual_tmrca_details,
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


def test_pair_details_attach_focal_carrier_counts():
    ts = msprime.sim_ancestry(
        samples=[msprime.SampleSet(4, ploidy=2)], population_size=100,
        sequence_length=1_000, recombination_rate=1e-7,
        model=msprime.StandardCoalescent(), random_seed=121,
    )
    counts = np.asarray([0, 1, 2, 0], dtype=np.int8)
    result = within_individual_tmrca_details(ts, 500, 150, counts)
    assert len(result) == 4
    assert result["focal_carrier_copies"].tolist() == counts.tolist()
    assert result["tmrca_lt_threshold"].dtype == bool


def test_retained_calibration_uses_full_null_and_first_retained_replicates():
    stats = pd.DataFrame({
        "selection_coefficient": [0.0, 0.0, 0.0, 0.1, 0.1, 0.1],
        "replicate": [0, 1, 2, 8, 2, 5],
        "focal_allele_outcome": [
            "lost", "lost", "lost", "segregating", "lost", "segregating"
        ],
        "realized_population_allele_frequency": [0, 0, 0, 0.8, 0, 0.5],
        "center_fraction_recent": [0.01, 0.02, 0.03, 0.8, 0.02, 0.5],
    })
    neutral, retained = retained_sweep_calibration_table(
        stats,
        selection_coefficient=0.1,
        retained_replicates=2,
        sample_diploids=2_000,
    )
    assert len(neutral) == 3
    assert retained["replicate"].tolist() == [5, 8]
    assert retained["neutral_exceedances"].eq(0).all()
    assert retained["mc_p_upper"].eq(0.25).all()


def test_carrier_summary_compares_two_copy_with_zero_copy_pairs(tmp_path):
    details = pd.DataFrame({
        "replicate": [7] * 6,
        "focal_carrier_copies": [0, 0, 1, 1, 2, 2],
        "tmrca_generations": [1_000, 2_000, 800, 900, 10, 20],
        "tmrca_lt_threshold": [False, False, False, False, True, True],
    })
    selected = pd.DataFrame({
        "replicate": [7],
        "realized_population_allele_frequency": [0.5],
    })
    copies, direct, effects = summarize_carrier_tmrca(
        details, selected, threshold_generations=180, output_dir=tmp_path
    )
    assert copies["focal_carrier_copies"].tolist() == [0, 1, 2]
    assert direct["focal_carrier_copies"].tolist() == [0, 2]
    assert effects.iloc[0]["carrier_fraction_tmrca_lt_threshold"] == 1
    assert effects.iloc[0]["noncarrier_fraction_tmrca_lt_threshold"] == 0
    assert (tmp_path / "carrier_vs_noncarrier_tmrca_ecdf.png").exists()


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

    reused = validate_recent_sweep_grid(
        tmp_path / "reused",
        population_size=100,
        sample_diploids=50,
        sequence_length=200_000,
        selection_coefficients=(0.75,),
        age_generations=30,
        mutation_rate=1.25e-8,
        recombination_rate=2e-7,
        selected_replicates=1,
        workers=2,
        reuse_null_from=tmp_path,
        save_trees=False,
        seed=271828,
    )
    assert reused["neutral_replicates"] == 9
    assert reused["neutral_reused_from"] == str(tmp_path.resolve())
    power = pd.read_csv(tmp_path / "reused" / "power_summary.tsv", sep="\t").iloc[0]
    assert power["n_focal_allele_lost"] + power["n_focal_allele_present"] == 1
    assert "center_fraction_recent_power_p_le_0_05_given_focal_allele_present" in power
