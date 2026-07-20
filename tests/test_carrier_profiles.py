import msprime
import numpy as np

from gamma_smc_aou.carrier_profiles import (
    pair_tmrca_profile_by_focal_copy,
    plot_carrier_profile_figure,
)


def test_pair_tmrca_profile_reports_mean_ci_for_hom_ref_and_hom_alt(tmp_path):
    ts = msprime.sim_ancestry(
        samples=[msprime.SampleSet(4, ploidy=2)],
        population_size=100,
        sequence_length=1_000,
        recombination_rate=1e-7,
        model=msprime.StandardCoalescent(),
        random_seed=818,
    )
    carrier_counts = np.asarray([0, 0, 2, 2], dtype=np.int8)
    profile = pair_tmrca_profile_by_focal_copy(
        ts, np.asarray([0, 400, 500, 600, 999]), carrier_counts
    )
    assert set(profile["pair_class"]) == {"hom_ref", "hom_alt"}
    assert profile["n_pairs"].eq(2).all()
    assert (
        profile["ci95_lower_generations"] <= profile["mean_tmrca_generations"]
    ).all()
    assert (
        profile["mean_tmrca_generations"] <= profile["ci95_upper_generations"]
    ).all()

    output = tmp_path / "profile.png"
    plot_carrier_profile_figure(
        profile,
        replicate=7,
        population_allele_frequency=0.5,
        sequence_length=1_000,
        center=500,
        zoom_half_width=100,
        output_path=output,
    )
    assert output.exists()
