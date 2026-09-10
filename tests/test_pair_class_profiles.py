import subprocess

import msprime
import numpy as np
import pandas as pd
import pytest
import zstandard

from gamma_smc_aou.fresh_power import REPO, ordered_nodes
from gamma_smc_aou.pair_class_profiles import (
    build_frame,
    pair_classes,
    posterior_arrays,
    summarize_decoded,
    summarize_truth,
)


CFG = {
    "decoder_scaled_mutation_rate": 0.00075,
    "mutation_rate": 1.25e-8,
    "generation_time_years": 25,
    "tmrca_cutoffs_years": [5000, 10000, 20000, 30000, 40000, 50000],
}


def test_classes_include_cross_individual_pairs_of_heterozygotes():
    # Every diploid is heterozygous. There are still alt/alt and ref/ref pairs
    # between their haplotypes, which a homozygous-individual filter would lose.
    carriers = np.array([1, 0, 1, 0, 1, 0, 1, 0])
    pairs = np.array([[0, 2], [1, 3], [0, 1], [2, 3], [4, 6], [5, 7]])
    np.testing.assert_array_equal(pair_classes(carriers, pairs), [2, 0, 1, 1, 2, 0])
    with pytest.raises(ValueError, match="indices"):
        pair_classes(carriers, np.array([[0, 0]]))


def test_truth_class_counts_match_explicit_pair_ages_and_retain_fixed():
    ts = msprime.sim_ancestry(
        samples=4,
        population_size=700,
        sequence_length=1000,
        recombination_rate=1e-6,
        random_seed=20380101,
    )
    pairs = np.array([[0, 2], [1, 3], [0, 1], [2, 3], [4, 6], [5, 7]])
    positions = np.array([0, 200, 500, 999])
    classes = pair_classes(np.ones(8, dtype=bool), pairs)
    counts, sums = summarize_truth(ts, pairs, classes, positions, CFG)
    nodes = ordered_nodes(ts)
    expected_ages = np.array(
        [[ts.at(p).tmrca(nodes[a], nodes[b]) * 25 for a, b in pairs] for p in positions]
    )
    np.testing.assert_array_equal(
        counts[2],
        (expected_ages[:, :, None] < np.array(CFG["tmrca_cutoffs_years"])).sum(axis=1),
    )
    np.testing.assert_allclose(sums[3], expected_ages.sum(axis=1))
    frame = build_frame(counts, sums, classes, positions, CFG, "truth")
    assert len(frame) == len(positions) * 4
    assert frame[frame.pair_class == "ref/ref"].n_pairs.eq(0).all()
    assert frame[frame.pair_class == "ref/ref"].mean_tmrca_years.isna().all()
    assert frame[frame.pair_class == "ref/ref"].frac_recent_50000.isna().all()
    np.testing.assert_array_equal(counts[2], counts[3])


def test_padded_posterior_classes_reproduce_native_mean_counts(tmp_path):
    n_pairs, n_positions = 11, 3
    rng = np.random.default_rng(20380101)
    data = np.empty((2, 2, n_positions, 8), dtype=np.float32)
    data[:, 0] = rng.uniform(1, 10, (2, n_positions, 8))
    data[:, 1] = rng.uniform(20, 500, (2, n_positions, 8))
    data[1, :, :, 3:] = np.nan  # padding must never enter any class
    raw = tmp_path / "posteriors.zst"
    raw.write_bytes(
        zstandard.ZstdCompressor(write_checksum=True).compress(data.tobytes())
    )
    alpha, beta = posterior_arrays(raw, n_pairs, n_positions)
    classes = np.arange(n_pairs) % 3
    counts, sums = summarize_decoded(alpha, beta, classes, CFG)
    np.testing.assert_array_equal(counts[:3].sum(axis=0), counts[3])
    units = float(np.float32(0.00075)) / (2 * 1.25e-8) * 25
    np.testing.assert_allclose(
        sums[3], (alpha.astype(float) / beta.astype(float) * units).sum(axis=1)
    )
    positions = tmp_path / "positions.txt"
    np.savetxt(positions, [0, 10, 20], fmt="%d")
    helper = REPO / "bin/summarize_recent_rules"
    subprocess.run(
        [
            str(helper),
            str(raw),
            str(n_pairs),
            str(positions),
            "0.00075",
            "1.25e-8",
            "25",
            ",".join(map(str, CFG["tmrca_cutoffs_years"])),
            str(tmp_path / "native"),
        ],
        check=True,
        capture_output=True,
    )
    original = (
        pd.read_csv(tmp_path / "native/mean.tsv", sep="\t").iloc[:, 1:].to_numpy()
    )
    np.testing.assert_array_equal(np.rint(original * n_pairs), counts[3])
    raw.write_bytes(zstandard.ZstdCompressor().compress(data.tobytes()[:-4]))
    with pytest.raises(ValueError, match="Truncated"):
        posterior_arrays(raw, n_pairs, n_positions)
