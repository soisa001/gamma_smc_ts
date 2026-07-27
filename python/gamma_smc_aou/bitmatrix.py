"""Reader for the packed per-pair recent-coalescence bit matrix.

The decoder writes one bit per (haplotype pair, output position, time
threshold): set when that pair's posterior TMRCA is called below the threshold.
At 100,000 pairs and stride 1000 that is ~3 GB per chromosome per threshold,
against ~200 GB for the raw alpha/beta posteriors the same run would produce.

FILE LAYOUT
-----------
The ``.bits`` file is a bare concatenation of independent zstd frames, one per
block of pairs, in whatever order the worker threads finished. The companion
``.bits.meta`` JSON carries ``frame_offsets``/``frame_sizes`` indexed by block,
so any single block decompresses on its own without reading the rest of the
file -- that is what makes ``pair_profiles`` cheap for a handful of pairs.

Uncompressed payload of one block, with ``K`` thresholds, ``C`` eight-pair
chunks and ``P`` output positions, is C-contiguous ``payload[k][c][pos]``, one
byte per entry. Bit ``p`` of that byte (bit 0 = least significant) is the call
for pair ``block_first_pair + 8 * c + p``. Padding bits past ``n_pairs`` are 0.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import zstandard


PAIRS_PER_BYTE = 8

# byte value -> number of set bits, for counting without unpacking to bool
_POPCOUNT = np.unpackbits(
    np.arange(256, dtype=np.uint8)[:, None], axis=1
).sum(axis=1).astype(np.int64)


def meta_path(path: str | Path) -> Path:
    path = Path(path)
    return path.with_name(path.name + ".meta")


def read_meta(path: str | Path) -> dict:
    """Load the sidecar metadata for a bit-matrix file."""
    with meta_path(path).open(encoding="utf-8") as handle:
        meta = json.load(handle)
    if meta.get("format") != "gamma_smc_bitmatrix_v1":
        raise ValueError(f"unexpected bit-matrix format: {meta.get('format')!r}")
    if int(meta["pairs_per_byte"]) != PAIRS_PER_BYTE:
        raise ValueError("only 8 pairs per byte is supported")
    return meta


def output_positions(meta: dict) -> np.ndarray:
    """0-based genomic positions of the columns, expanding a plain stride."""
    if meta.get("positions_are_strided"):
        stride = int(meta["stride"])
        first = int(meta["first_position"])
        return first + stride * np.arange(int(meta["n_positions"]), dtype=np.int64)
    return np.asarray(meta["output_positions"], dtype=np.int64)


def pairs(meta: dict) -> np.ndarray:
    """(n_pairs, 2) array of 0-based haplotype indices, in output order."""
    return np.asarray(meta["pairs"], dtype=np.int64)


def _block_bytes(path: Path, meta: dict, block_index: int) -> np.ndarray:
    """Decompress one block into its (K, C, P) byte array."""
    offset = int(meta["frame_offsets"][block_index])
    size = int(meta["frame_sizes"][block_index])
    expected = int(meta["frame_uncompressed_sizes"][block_index])

    with path.open("rb") as handle:
        handle.seek(offset)
        frame = handle.read(size)
    if len(frame) != size:
        raise EOFError(f"block {block_index} is truncated in {path}")

    raw = zstandard.ZstdDecompressor().decompress(frame, max_output_size=expected)
    if len(raw) != expected:
        raise ValueError(
            f"block {block_index} decompressed to {len(raw)} bytes, expected {expected}"
        )

    n_thresholds = len(meta["thresholds_years"])
    n_positions = int(meta["n_positions"])
    n_chunks = expected // (n_thresholds * n_positions) if n_positions else 0
    return np.frombuffer(raw, dtype=np.uint8).reshape(n_thresholds, n_chunks, n_positions)


def read_block(path: str | Path, block_index: int, meta: dict | None = None) -> np.ndarray:
    """One block as a bool array of shape (n_thresholds, n_positions, block_n_pairs)."""
    path = Path(path)
    meta = meta if meta is not None else read_meta(path)
    packed = _block_bytes(path, meta, block_index)
    block_n_pairs = int(meta["block_n_pairs"][block_index])
    # (K, C, P) -> (K, P, C) so unpacking the byte axis yields pair order.
    bits = np.unpackbits(np.ascontiguousarray(packed.transpose(0, 2, 1)), axis=-1, bitorder="little")
    return bits[:, :, :block_n_pairs].astype(bool)


def iter_blocks(path: str | Path, meta: dict | None = None):
    """Yield (block_index, first_pair, bool array) for every block in order."""
    path = Path(path)
    meta = meta if meta is not None else read_meta(path)
    for block_index in range(int(meta["n_blocks"])):
        yield (
            block_index,
            int(meta["block_first_pair"][block_index]),
            read_block(path, block_index, meta),
        )


def position_counts(
    path: str | Path,
    meta: dict | None = None,
    pair_indices: np.ndarray | list[int] | None = None,
) -> np.ndarray:
    """Number of pairs called recent at each position, per threshold.

    Shape (n_thresholds, n_positions). With ``pair_indices`` the count is
    restricted to those pairs, which is the carrier-versus-noncarrier split;
    only the blocks containing them are read.
    """
    path = Path(path)
    meta = meta if meta is not None else read_meta(path)
    n_thresholds = len(meta["thresholds_years"])
    n_positions = int(meta["n_positions"])
    counts = np.zeros((n_thresholds, n_positions), dtype=np.int64)

    if pair_indices is None:
        for block_index in range(int(meta["n_blocks"])):
            packed = _block_bytes(path, meta, block_index)
            counts += _POPCOUNT[packed].sum(axis=1)
        return counts

    wanted = np.unique(np.asarray(pair_indices, dtype=np.int64))
    if wanted.size and (wanted[0] < 0 or wanted[-1] >= int(meta["n_pairs"])):
        raise IndexError("pair index out of range")
    pair_block = int(meta["pair_block"])
    for block_index in np.unique(wanted // pair_block):
        block_index = int(block_index)
        first_pair = int(meta["block_first_pair"][block_index])
        block_n_pairs = int(meta["block_n_pairs"][block_index])
        local = wanted[(wanted >= first_pair) & (wanted < first_pair + block_n_pairs)] - first_pair
        packed = _block_bytes(path, meta, block_index)
        bits = np.unpackbits(
            np.ascontiguousarray(packed.transpose(0, 2, 1)), axis=-1, bitorder="little"
        )
        counts += bits[:, :, local].sum(axis=2).astype(np.int64)
    return counts


def pair_profiles(
    path: str | Path,
    pair_indices: np.ndarray | list[int],
    meta: dict | None = None,
) -> np.ndarray:
    """Calls for specific pairs: (n_thresholds, n_positions, len(pair_indices)) bool.

    Pairs come back in the order requested. Only the blocks holding them are
    decompressed.
    """
    path = Path(path)
    meta = meta if meta is not None else read_meta(path)
    requested = np.asarray(pair_indices, dtype=np.int64)
    if requested.size and (requested.min() < 0 or requested.max() >= int(meta["n_pairs"])):
        raise IndexError("pair index out of range")

    n_thresholds = len(meta["thresholds_years"])
    n_positions = int(meta["n_positions"])
    out = np.zeros((n_thresholds, n_positions, requested.size), dtype=bool)

    pair_block = int(meta["pair_block"])
    for block_index in np.unique(requested // pair_block):
        block_index = int(block_index)
        first_pair = int(meta["block_first_pair"][block_index])
        block_n_pairs = int(meta["block_n_pairs"][block_index])
        selector = (requested >= first_pair) & (requested < first_pair + block_n_pairs)
        packed = _block_bytes(path, meta, block_index)
        bits = np.unpackbits(
            np.ascontiguousarray(packed.transpose(0, 2, 1)), axis=-1, bitorder="little"
        )
        out[:, :, selector] = bits[:, :, requested[selector] - first_pair].astype(bool)
    return out


def to_frame(
    path: str | Path,
    meta: dict | None = None,
    pair_indices: np.ndarray | list[int] | None = None,
) -> pd.DataFrame:
    """Per-position counts and fractions, matching the decoder's summary TSV.

    Columns: ``position_0based``, ``position_1based``, ``n_pairs``, and
    ``n_recent_<years>`` / ``frac_recent_<years>`` per threshold. Recomputing
    this from the bit matrix and comparing against the decoder's own
    ``--recent_summary`` is a direct end-to-end check of both.
    """
    path = Path(path)
    meta = meta if meta is not None else read_meta(path)
    counts = position_counts(path, meta, pair_indices)
    positions = output_positions(meta)
    n_pairs = (
        int(meta["n_pairs"])
        if pair_indices is None
        else int(np.unique(np.asarray(pair_indices)).size)
    )

    frame = pd.DataFrame({
        "position_0based": positions,
        "position_1based": positions + 1,
        "n_pairs": n_pairs,
    })
    for index, years in enumerate(meta["thresholds_years"]):
        suffix = _format_threshold(years)
        frame[f"n_recent_{suffix}"] = counts[index]
        frame[f"frac_recent_{suffix}"] = counts[index] / n_pairs if n_pairs else np.nan
    return frame


def _format_threshold(years: float) -> str:
    """Match the C++ column suffix: integral thresholds print without a point."""
    return str(int(years)) if float(years).is_integer() else str(years)
