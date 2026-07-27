"""End-to-end tests of the compiled decoder's new outputs.

These need GAMMA_SMC_BIN and are skipped without it. What they check that the
pure-Python tests cannot: that the bit matrix the C++ writes agrees with the
summary TSV the same run produced, that threading does not change the answer,
and that random pair sampling is reproducible.
"""

import json
import os
import subprocess

import numpy as np
import pandas as pd
import pytest

from gamma_smc_aou import bitmatrix

from test_simulation import diploid_ts


pytestmark = pytest.mark.skipif(
    "GAMMA_SMC_BIN" not in os.environ, reason="compiled Linux binary not available"
)

THRESHOLDS = (4500.0, 10000.0)
STRIDE = 1000


def decode(tmp_path, ts, *, name, extra=(), n_random_pairs=0, threads=1, bits=True):
    source = tmp_path / f"{name}.trees"
    ts.dump(source)
    summary = tmp_path / f"{name}.tsv"
    command = [
        os.environ["GAMMA_SMC_BIN"],
        "--input", str(source), "--input_format", "trees",
        "--scaled_mutation_rate", "0.0008",
        "--recombination_to_mutation_ratio", "0.05",
        "--unscaled_mutation_rate", "2e-7",
        "--recent_threshold_years", ",".join(str(value) for value in THRESHOLDS),
        "--generation_time", "25",
        "--recent_summary", str(summary),
        "--output_at_hets=false",
        "--output_at_stride", str(STRIDE),
        "--threads", str(threads),
    ]
    if n_random_pairs:
        command += ["--n_random_pairs", str(n_random_pairs), "--pairs_seed", "20240727"]
    else:
        command += ["--only_within"]
    if bits:
        command += ["--recent_bitmatrix", str(tmp_path / f"{name}.bits")]
    command += list(extra)

    completed = subprocess.run(command, check=True, text=True, capture_output=True)
    return summary, tmp_path / f"{name}.bits", completed


def test_bitmatrix_counts_match_the_summary_tsv(tmp_path):
    ts = diploid_ts(length=60_000)
    summary, bits_path, _ = decode(tmp_path, ts, name="within", n_random_pairs=0)

    frame = pd.read_csv(summary, sep="\t")
    counts = bitmatrix.position_counts(bits_path)
    meta = bitmatrix.read_meta(bits_path)

    assert list(meta["thresholds_years"]) == list(THRESHOLDS)
    assert len(frame) == meta["n_positions"]
    np.testing.assert_array_equal(
        bitmatrix.output_positions(meta), frame["position_0based"].to_numpy()
    )
    np.testing.assert_array_equal(counts[0], frame["n_recent_4500"].to_numpy())
    np.testing.assert_array_equal(counts[1], frame["n_recent_10000"].to_numpy())


def test_summary_keeps_the_legacy_schema(tmp_path):
    ts = diploid_ts(length=60_000)
    summary, _, _ = decode(tmp_path, ts, name="legacy")
    frame = pd.read_csv(summary, sep="\t")

    for column in (
        "position_0based", "position_1based", "n_pairs",
        "mean_p_tmrca_lt_threshold", "mean_tmrca_generations",
    ):
        assert column in frame.columns

    # The legacy alias is the first threshold's mean probability.
    np.testing.assert_allclose(
        frame["mean_p_tmrca_lt_threshold"].to_numpy(),
        frame["mean_p_lt_4500"].to_numpy(),
        rtol=0, atol=0,
    )
    assert frame["mean_p_tmrca_lt_threshold"].between(0, 1).all()
    assert frame["frac_recent_4500"].between(0, 1).all()
    # A longer horizon can only catch more pairs.
    assert (frame["n_recent_10000"] >= frame["n_recent_4500"]).all()
    assert (frame["mean_p_lt_10000"] >= frame["mean_p_lt_4500"] - 1e-6).all()


def test_thread_count_does_not_change_the_result(tmp_path):
    ts = diploid_ts(length=60_000)
    one, bits_one, _ = decode(tmp_path, ts, name="t1", n_random_pairs=64, threads=1)
    many, bits_many, _ = decode(tmp_path, ts, name="t8", n_random_pairs=64, threads=8)

    frame_one = pd.read_csv(one, sep="\t")
    frame_many = pd.read_csv(many, sep="\t")
    pd.testing.assert_frame_equal(frame_one, frame_many)

    np.testing.assert_array_equal(
        bitmatrix.position_counts(bits_one), bitmatrix.position_counts(bits_many)
    )


def test_pair_block_does_not_change_the_result(tmp_path):
    ts = diploid_ts(length=60_000)
    small, bits_small, _ = decode(
        tmp_path, ts, name="b8", n_random_pairs=50, threads=4, extra=("--pair_block", "8")
    )
    large, bits_large, _ = decode(
        tmp_path, ts, name="b512", n_random_pairs=50, threads=4, extra=("--pair_block", "512")
    )
    pd.testing.assert_frame_equal(pd.read_csv(small, sep="\t"), pd.read_csv(large, sep="\t"))
    np.testing.assert_array_equal(
        bitmatrix.position_counts(bits_small), bitmatrix.position_counts(bits_large)
    )


def test_random_pairs_are_reproducible_and_well_formed(tmp_path):
    ts = diploid_ts(length=60_000)
    _, bits_a, _ = decode(tmp_path, ts, name="ra", n_random_pairs=20)
    _, bits_b, _ = decode(tmp_path, ts, name="rb", n_random_pairs=20)

    meta_a = bitmatrix.read_meta(bits_a)
    meta_b = bitmatrix.read_meta(bits_b)
    assert meta_a["pairs"] == meta_b["pairs"]
    assert meta_a["n_pairs"] == 20

    pairs = bitmatrix.pairs(meta_a)
    n_haplotypes = 2 * len(ts.individuals())
    assert pairs.min() >= 0 and pairs.max() < n_haplotypes
    assert (pairs[:, 0] < pairs[:, 1]).all()
    assert len({tuple(pair) for pair in pairs}) == 20


def test_exclude_within_drops_same_individual_pairs(tmp_path):
    ts = diploid_ts(length=60_000)
    _, bits_path, _ = decode(
        tmp_path, ts, name="xw", n_random_pairs=20, extra=("--exclude_within",)
    )
    pairs = bitmatrix.pairs(bitmatrix.read_meta(bits_path))
    assert ((pairs[:, 0] >> 1) != (pairs[:, 1] >> 1)).all()


def test_pairs_file_is_honoured(tmp_path):
    ts = diploid_ts(length=60_000)
    wanted = [(0, 3), (1, 6), (2, 7)]
    pairs_file = tmp_path / "pairs.txt"
    pairs_file.write_text(
        "# haplotype pairs\n" + "\n".join(f"{i}\t{j}" for i, j in wanted) + "\n"
    )
    _, bits_path, _ = decode(
        tmp_path, ts, name="pf", n_random_pairs=0, extra=("--pairs_file", str(pairs_file))
    )
    meta = bitmatrix.read_meta(bits_path)
    assert [tuple(pair) for pair in meta["pairs"]] == wanted


def test_asking_for_more_pairs_than_exist_fails_cleanly(tmp_path):
    ts = diploid_ts(length=60_000)
    with pytest.raises(subprocess.CalledProcessError) as excinfo:
        decode(tmp_path, ts, name="toomany", n_random_pairs=10_000)
    assert "exceeds" in excinfo.value.stdout


def test_lookup_tables_self_check_within_tolerance(tmp_path):
    ts = diploid_ts(length=60_000)
    _, _, completed = decode(tmp_path, ts, name="selfcheck")
    log = completed.stdout

    error = float(log.split("max abs error vs Boost")[1].split()[0])
    assert error < 1e-3, f"P(T<t) table error {error:g} is worse than expected"

    disagreement = float(log.split("disagreement with Boost")[1].split()[0])
    assert disagreement < 1e-4, f"call table disagreement {disagreement:g} is too high"


def test_mean_call_rule_is_a_threshold_on_the_posterior_mean(tmp_path):
    ts = diploid_ts(length=60_000)
    summary, _, _ = decode(
        tmp_path, ts, name="meanrule", extra=("--recent_call", "mean"), bits=False
    )
    frame = pd.read_csv(summary, sep="\t")
    # 10000 years at 25 years/generation is 400 generations.
    assert frame["frac_recent_10000"].between(0, 1).all()
    assert (frame["n_recent_10000"] >= frame["n_recent_4500"]).all()


def test_raw_posteriors_still_round_trip(tmp_path):
    ts = diploid_ts(length=60_000)
    raw = tmp_path / "posteriors.zst"
    summary, _, _ = decode(
        tmp_path, ts, name="raw", n_random_pairs=20, threads=4,
        extra=("--output", str(raw)), bits=False,
    )
    meta = json.loads((tmp_path / "posteriors.zst.meta").read_text())
    assert meta["num_pairs"] == 20
    assert meta["chunk_size"] == 8
    assert meta["sequence_length"] == len(pd.read_csv(summary, sep="\t"))

    import zstandard

    raw_floats = np.frombuffer(
        zstandard.ZstdDecompressor().stream_reader(raw.open("rb")).read(), dtype=np.float32
    )
    n_chunks = -(-meta["num_pairs"] // meta["chunk_size"])
    assert raw_floats.size == n_chunks * 2 * meta["sequence_length"] * meta["chunk_size"]


def test_pair_subset_counts_agree_with_the_whole(tmp_path):
    ts = diploid_ts(length=60_000)
    _, bits_path, _ = decode(tmp_path, ts, name="subset", n_random_pairs=40, threads=4)
    meta = bitmatrix.read_meta(bits_path)
    total = bitmatrix.position_counts(bits_path, meta)
    first = bitmatrix.position_counts(bits_path, meta, pair_indices=range(0, 20))
    second = bitmatrix.position_counts(bits_path, meta, pair_indices=range(20, 40))
    np.testing.assert_array_equal(total, first + second)
