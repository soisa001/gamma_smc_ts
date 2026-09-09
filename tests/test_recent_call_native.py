"""Small native-decoder comparison against direct incomplete-gamma evaluation."""

import sys
import json
import subprocess

import msprime
import numpy as np
import pandas as pd
import pytest
import zstandard
from scipy.special import gammainc

from gamma_smc_aou.decoder import run_within_decoder
from gamma_smc_aou.fresh_power import REPO


@pytest.mark.skipif(sys.platform != "linux", reason="Native Linux replay helper")
def test_saved_posterior_replay_matches_cdf_and_rejects_truncation(tmp_path):
    decoder = REPO / "bin/gamma_smc"
    helper = REPO / "bin/summarize_recent_rules"
    if not decoder.exists() or not helper.exists():
        pytest.skip("Native tools not built")
    ts = msprime.sim_ancestry(
        samples=4,
        population_size=2000,
        sequence_length=50000,
        recombination_rate=1e-8,
        random_seed=20380101,
    )
    ts = msprime.sim_mutations(ts, rate=1.25e-8, random_seed=20380102)
    source = tmp_path / "input.trees"
    ts.dump(source)
    pairs = tmp_path / "pairs.tsv"
    np.savetxt(
        pairs, np.column_stack(np.triu_indices(8, k=1)), fmt="%d", delimiter="\t"
    )
    positions = tmp_path / "positions.txt"
    np.savetxt(positions, np.arange(0, 50000, 1000), fmt="%d")
    cutoffs = [5000, 10000, 20000, 30000, 40000, 50000]
    raw = tmp_path / "posterior.zst"
    summary = tmp_path / "native.tsv"
    run_within_decoder(
        decoder,
        source,
        summary,
        raw_output=raw,
        mutation_rate=1.25e-8,
        threshold_years=cutoffs,
        generation_time=25,
        output_at_stride=-1,
        output_positions_file=positions,
        only_within=False,
        pairs_file=pairs,
        threads=1,
        recent_call="mean",
        vcf_position_transform="one_based",
        extra_args=["--no_recent_probability"],
    )
    meta = json.loads(raw.with_name(raw.name + ".meta").read_text())
    assert (
        meta["num_pairs"] == 28
        and meta["chunk_size"] == 8
        and meta["sequence_length"] == 50
    )
    command = [
        str(helper),
        str(raw),
        "28",
        str(positions),
        "0.00075",
        "1.25e-8",
        "25",
        ",".join(map(str, cutoffs)),
        str(tmp_path / "replayed"),
    ]
    subprocess.run(command, check=True, capture_output=True)
    expected = {
        name: np.zeros((50, 6), dtype=int) for name in ("median", "prob80", "prob90")
    }
    with raw.open("rb") as stream:
        with zstandard.ZstdDecompressor().stream_reader(stream) as reader:
            data = np.frombuffer(reader.read(), dtype=np.float32).reshape(4, 2, 50, 8)
    two_ne = float(np.float32(0.00075)) / (2 * 1.25e-8)
    for i, chunk in enumerate(data):
        n = min(8, 28 - i * 8)
        alpha, beta = chunk[:, :, :n].astype(float)
        probabilities = gammainc(
            alpha[:, :, None], beta[:, :, None] * np.array(cutoffs) / 25 / two_ne
        )
        for name, q in (("median", 0.5), ("prob80", 0.8), ("prob90", 0.9)):
            expected[name] += (probabilities >= q).sum(axis=1)
    columns = [f"frac_recent_{t}" for t in cutoffs]
    original = pd.read_csv(summary, sep="\t")[columns].to_numpy()
    replayed = pd.read_csv(tmp_path / "replayed/mean.tsv", sep="\t")[columns].to_numpy()
    np.testing.assert_array_equal(np.rint(original * 28), np.rint(replayed * 28))
    for name in expected:
        actual = pd.read_csv(tmp_path / f"replayed/{name}.tsv", sep="\t")[
            columns
        ].to_numpy()
        np.testing.assert_array_equal(np.rint(actual * 28), expected[name])
    raw.write_bytes(raw.read_bytes()[:-8])
    assert subprocess.run(command, capture_output=True).returncode != 0


@pytest.mark.skipif(
    sys.platform != "linux", reason="Native decoder uses the pinned Linux environment"
)
def test_native_hard_calls_match_exact_gamma_evaluation(tmp_path):
    decoder = REPO / "bin/gamma_smc"
    if not decoder.is_file():
        pytest.skip("Native decoder not built")
    ts = msprime.sim_ancestry(
        samples=4,
        population_size=2000,
        sequence_length=50000,
        recombination_rate=1e-8,
        random_seed=20380101,
    )
    ts = msprime.sim_mutations(ts, rate=1.25e-8, random_seed=20380102)
    tree_file = tmp_path / "input.trees"
    ts.dump(tree_file)
    pairs_file = tmp_path / "pairs.tsv"
    np.savetxt(
        pairs_file, np.column_stack(np.triu_indices(8, k=1)), fmt="%d", delimiter="\t"
    )
    positions_file = tmp_path / "positions.txt"
    np.savetxt(positions_file, np.arange(0, 50000, 1000), fmt="%d")
    cutoffs = [5000, 10000, 20000, 30000, 40000, 50000]
    columns = [f"n_recent_{t}" for t in cutoffs]
    profiles = {}
    for name, call, probability in [
        ("mean", "mean", 0.5),
        ("median", "median", 0.5),
        ("prob80", "prob", 0.8),
        ("prob90", "prob", 0.9),
    ]:
        for exact in (False, True):
            path = tmp_path / f"{name}_{exact}.tsv"
            run_within_decoder(
                decoder,
                tree_file,
                path,
                mutation_rate=1.25e-8,
                threshold_years=cutoffs,
                generation_time=25,
                output_at_stride=-1,
                output_positions_file=positions_file,
                only_within=False,
                pairs_file=pairs_file,
                threads=1,
                recent_call=call,
                recent_call_probability=probability,
                vcf_position_transform="one_based",
                extra_args=["--no_recent_probability"]
                + (["--exact_recent_stats"] if exact else []),
            )
            frame = pd.read_csv(path, sep="\t")
            assert (frame.n_pairs == 28).all()
            profiles[name, exact] = frame[columns].to_numpy()
        np.testing.assert_array_equal(profiles[name, False], profiles[name, True])
    assert np.any(profiles["median", False] > profiles["mean", False])
    assert np.any(profiles["prob80", False] > profiles["prob90", False])
    for lower, upper in (
        ("mean", "median"),
        ("prob80", "median"),
        ("prob90", "prob80"),
    ):
        assert np.all(profiles[lower, False] <= profiles[upper, False])
