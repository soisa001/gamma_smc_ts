"""Bit-matrix reader tests.

These build the file the way the C++ writer does -- independent zstd frames in
arbitrary order, indexed by the sidecar metadata -- so they pin the on-disk
contract without needing the compiled binary.
"""

import json

import numpy as np
import pytest
import zstandard

from gamma_smc_aou import bitmatrix


PAIRS_PER_BYTE = 8


def build_bitmatrix(tmp_path, calls, *, pair_block, thresholds=(4500.0, 10000.0),
                    stride=1000, shuffle_frames=False):
    """Write a bit matrix from a (n_thresholds, n_positions, n_pairs) bool array.

    Mirrors the C++ writer exactly: payload[k][chunk][position], bit p of each
    byte is pair block_first_pair + 8*chunk + p, frames concatenated in an order
    the metadata does not have to agree with.
    """
    calls = np.asarray(calls, dtype=bool)
    n_thresholds, n_positions, n_pairs = calls.shape
    assert n_thresholds == len(thresholds)
    assert pair_block % PAIRS_PER_BYTE == 0

    n_blocks = (n_pairs + pair_block - 1) // pair_block
    frames = []
    for block in range(n_blocks):
        first_pair = block * pair_block
        block_n_pairs = min(pair_block, n_pairs - first_pair)
        n_chunks = (block_n_pairs + PAIRS_PER_BYTE - 1) // PAIRS_PER_BYTE

        payload = np.zeros((n_thresholds, n_chunks, n_positions), dtype=np.uint8)
        for k in range(n_thresholds):
            for chunk in range(n_chunks):
                for lane in range(PAIRS_PER_BYTE):
                    pair = first_pair + chunk * PAIRS_PER_BYTE + lane
                    if pair >= n_pairs:
                        continue  # padding bits stay 0
                    payload[k, chunk] |= calls[k, :, pair].astype(np.uint8) << lane

        raw = payload.tobytes()
        frames.append({
            "block": block,
            "first_pair": first_pair,
            "n_pairs": block_n_pairs,
            "raw_size": len(raw),
            "data": zstandard.ZstdCompressor(level=1).compress(raw),
        })

    order = list(range(n_blocks))
    if shuffle_frames:
        order = order[::-1]

    path = tmp_path / "calls.bits"
    offsets = [0] * n_blocks
    with path.open("wb") as handle:
        for block in order:
            offsets[block] = handle.tell()
            handle.write(frames[block]["data"])

    meta = {
        "format": "gamma_smc_bitmatrix_v1",
        "n_pairs": n_pairs,
        "n_positions": n_positions,
        "pair_block": pair_block,
        "pairs_per_byte": PAIRS_PER_BYTE,
        "n_blocks": n_blocks,
        "thresholds_years": list(thresholds),
        "positions_are_strided": True,
        "stride": stride,
        "first_position": 0,
        "pairs": [[2 * i, 2 * i + 1] for i in range(n_pairs)],
        "block_first_pair": [frame["first_pair"] for frame in frames],
        "block_n_pairs": [frame["n_pairs"] for frame in frames],
        "frame_offsets": offsets,
        "frame_sizes": [len(frame["data"]) for frame in frames],
        "frame_uncompressed_sizes": [frame["raw_size"] for frame in frames],
    }
    with path.with_name(path.name + ".meta").open("w", encoding="utf-8") as handle:
        json.dump(meta, handle)
    return path, calls


@pytest.fixture
def sample_matrix(tmp_path):
    rng = np.random.default_rng(20240727)
    # 21 pairs over a 16-pair block: exercises a partial final block and a
    # partial final chunk, where the padding bits live.
    calls = rng.random((2, 37, 21)) < 0.3
    path, truth = build_bitmatrix(tmp_path, calls, pair_block=16)
    return path, truth


def test_round_trip_recovers_every_bit(sample_matrix):
    path, truth = sample_matrix
    meta = bitmatrix.read_meta(path)
    recovered = np.concatenate(
        [block for _, _, block in bitmatrix.iter_blocks(path, meta)], axis=2
    )
    assert recovered.shape == truth.shape
    assert np.array_equal(recovered, truth)


def test_position_counts_match_a_direct_sum(sample_matrix):
    path, truth = sample_matrix
    assert np.array_equal(bitmatrix.position_counts(path), truth.sum(axis=2))


def test_padding_bits_are_not_counted(tmp_path):
    # 5 pairs in an 8-pair chunk: the three padding lanes must not show up as
    # recent calls even though every real pair is set at every position.
    calls = np.ones((1, 4, 5), dtype=bool)
    path, _ = build_bitmatrix(tmp_path, calls, pair_block=8, thresholds=(4500.0,))
    assert np.array_equal(bitmatrix.position_counts(path), np.full((1, 4), 5))


def test_frames_may_be_stored_out_of_order(tmp_path):
    rng = np.random.default_rng(7)
    calls = rng.random((2, 11, 40)) < 0.5
    path, truth = build_bitmatrix(tmp_path, calls, pair_block=8, shuffle_frames=True)
    assert np.array_equal(bitmatrix.position_counts(path), truth.sum(axis=2))
    recovered = np.concatenate(
        [block for _, _, block in bitmatrix.iter_blocks(path)], axis=2
    )
    assert np.array_equal(recovered, truth)


def test_pair_profiles_select_the_requested_pairs(sample_matrix):
    path, truth = sample_matrix
    wanted = [20, 0, 17, 3]
    profiles = bitmatrix.pair_profiles(path, wanted)
    assert profiles.shape == (2, 37, len(wanted))
    assert np.array_equal(profiles, truth[:, :, wanted])


def test_counts_can_be_restricted_to_a_pair_subset(sample_matrix):
    path, truth = sample_matrix
    subset = [1, 2, 3, 16, 17, 20]
    assert np.array_equal(
        bitmatrix.position_counts(path, pair_indices=subset),
        truth[:, :, subset].sum(axis=2),
    )


def test_out_of_range_pairs_are_rejected(sample_matrix):
    path, truth = sample_matrix
    with pytest.raises(IndexError):
        bitmatrix.pair_profiles(path, [truth.shape[2]])


def test_summary_frame_has_the_expected_schema(sample_matrix):
    path, truth = sample_matrix
    frame = bitmatrix.to_frame(path)
    assert list(frame["position_0based"]) == list(range(0, 37 * 1000, 1000))
    assert (frame["position_1based"] == frame["position_0based"] + 1).all()
    assert frame["n_pairs"].eq(21).all()
    assert np.array_equal(frame["n_recent_4500"].to_numpy(), truth[0].sum(axis=1))
    assert np.array_equal(frame["n_recent_10000"].to_numpy(), truth[1].sum(axis=1))
    np.testing.assert_allclose(
        frame["frac_recent_4500"].to_numpy(), truth[0].sum(axis=1) / 21
    )


def test_strided_positions_are_expanded(sample_matrix):
    path, _ = sample_matrix
    meta = bitmatrix.read_meta(path)
    positions = bitmatrix.output_positions(meta)
    assert positions[0] == 0
    assert positions[1] - positions[0] == 1000
    assert len(positions) == meta["n_positions"]


def test_unknown_format_is_rejected(tmp_path):
    path = tmp_path / "bad.bits"
    path.write_bytes(b"")
    with path.with_name(path.name + ".meta").open("w", encoding="utf-8") as handle:
        json.dump({"format": "something_else"}, handle)
    with pytest.raises(ValueError):
        bitmatrix.read_meta(path)
