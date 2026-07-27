#pragma once

// Packed per-pair recent-coalescence bit matrix.
//
// Writing the raw Gamma posterior (alpha, beta as float32) for every pair at
// every output position costs 8 bytes per cell: 100,000 pairs x 250,000 stride-1kb
// positions on chr1 is 200 GB. The scan only needs, per pair and position, whether
// that pair's posterior places the TMRCA below each of a small number of time
// thresholds, so one bit per (pair, position, threshold) is enough -- 64x smaller,
// and the per-position counts fall out of a popcount.
//
// FILE LAYOUT
// -----------
// The file is a bare concatenation of independent zstd frames, one per pair
// block, written in whatever order the worker threads finish. Nothing about the
// order is implied: `.meta` carries frame_offsets/frame_sizes indexed by block,
// so any single block can be decompressed on its own without touching the rest.
//
// Uncompressed payload of the frame for block b, with K thresholds,
// n_chunks = ceil(block_n_pairs / 8) and P output positions:
//
//     payload[k][c][pos]        k in [0,K), c in [0,n_chunks), pos in [0,P)
//
// laid out C-contiguous, one byte per (k, c, pos). Bit p of that byte (bit 0 =
// least significant) is the call for pair
//
//     block_first_pair + 8 * c + p
//
// under threshold k at output position pos. Padding bits, where the last chunk
// of the last block runs past n_pairs, are always 0.
//
// Position is the fastest-varying axis on purpose: each 8-pair chunk writes one
// contiguous P-byte run, so packing is a sequential store, and one pair's whole
// genome-wide profile is a strided read inside a single contiguous run.

#include "common.h"

struct BitMatrixFrame {
    long block_index = -1;
    long first_pair = 0;
    long n_pairs = 0;
    long offset = 0;
    long compressed_size = 0;
    long uncompressed_size = 0;
};

// Serialises frame writes and hands back the byte offset each frame landed at.
class BitMatrixWriter {
  public:
    ostream* _output;
    long _offset = 0;
    vector<BitMatrixFrame> _frames;

    explicit BitMatrixWriter(ostream* output) : _output(output) {}

    // Thread-safe: callers hold no lock, this one does. Returns the record so
    // the caller can keep it if it wants; the writer also retains a copy.
    BitMatrixFrame write_frame(
        long block_index,
        long first_pair,
        long n_pairs_in_block,
        const char* compressed,
        size_t compressed_size,
        size_t uncompressed_size
    ) {
        BitMatrixFrame frame;
        {
#ifdef _OPENMP
#pragma omp critical(gamma_smc_bitmatrix_write)
#endif
            {
                frame.block_index = block_index;
                frame.first_pair = first_pair;
                frame.n_pairs = n_pairs_in_block;
                frame.offset = _offset;
                frame.compressed_size = (long) compressed_size;
                frame.uncompressed_size = (long) uncompressed_size;
                _output->write(compressed, (std::streamsize) compressed_size);
                _offset += (long) compressed_size;
                _frames.push_back(frame);
            }
        }
        return frame;
    }

    // Frames are appended in completion order; sort so `.meta` can index by block.
    void sort_frames() {
        std::sort(
            _frames.begin(), _frames.end(),
            [](const BitMatrixFrame& a, const BitMatrixFrame& b) {
                return a.block_index < b.block_index;
            }
        );
    }
};
