#pragma once

// Construction of the haplotype-pair work list.
//
// Haplotypes are numbered 0..2*n_samples-1; haplotypes 2*s and 2*s+1 belong to
// diploid sample s. The scans in this branch need ~1e5 pairs drawn at random
// from the whole panel rather than the O(n^2) exhaustive list, so this header
// centralises every way the list can be built:
//
//   * within-individual only          (--only_within)
//   * first list against second list  (--samples_against)
//   * uniform random sample           (--n_random_pairs)
//   * explicit list from a file       (--pairs_file)
//   * exhaustive                      (default)

#include "common.h"

#include <boost/algorithm/string.hpp>

#include <random>
#include <sstream>
#include <unordered_set>

// Pack an unordered pair into one key so a hash set can deduplicate draws.
inline uint64_t pair_key(int i, int j) {
    return (((uint64_t) i) << 32) | ((uint64_t) (uint32_t) j);
}

// Uniform sample of `n_requested` distinct unordered haplotype pairs.
//
// Draws both haplotypes uniformly and rejects i == j, which is uniform over
// unordered pairs; `exclude_within` additionally rejects the two haplotypes of
// the same diploid. Pairs are returned sorted, which is deterministic for a
// given seed and keeps the per-site allele reads of a pair block local.
inline void sample_random_pairs(
    int n_haplotypes,
    long n_requested,
    uint64_t seed,
    bool exclude_within,
    vector<pair<int, int>>& out
) {
    if (n_haplotypes < 2) {
        cout << "Error: need at least two haplotypes to form a pair." << endl;
        exit(-1);
    }

    // Total number of distinct unordered pairs available under the current
    // exclusion rule; sampling more than that can never terminate.
    const double n_hap = (double) n_haplotypes;
    double available = n_hap * (n_hap - 1) / 2;
    if (exclude_within) {
        available -= (double) (n_haplotypes / 2);
    }
    if ((double) n_requested > available) {
        cout << boost::format(
            "Error: --n_random_pairs %ld exceeds the %.0f distinct haplotype pairs available "
            "from %d haplotypes%s.\n"
        ) % n_requested % available % n_haplotypes
          % (exclude_within ? " (within-individual pairs excluded)" : "");
        exit(-1);
    }

    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<int> pick(0, n_haplotypes - 1);

    std::unordered_set<uint64_t> seen;
    seen.reserve((size_t) (n_requested * 2));

    out.clear();
    out.reserve((size_t) n_requested);

    while ((long) out.size() < n_requested) {
        int i = pick(rng);
        int j = pick(rng);
        if (i == j) {
            continue;
        }
        if (i > j) {
            std::swap(i, j);
        }
        if (exclude_within && ((i >> 1) == (j >> 1))) {
            continue;
        }
        if (!seen.insert(pair_key(i, j)).second) {
            continue;
        }
        out.push_back(make_pair(i, j));
    }

    std::sort(out.begin(), out.end());
}

// Read an explicit pair list: two whitespace-separated 0-based haplotype
// indices per line. Blank lines and lines starting with '#' are ignored.
inline void read_pairs_file(
    const string& filename,
    int n_haplotypes,
    vector<pair<int, int>>& out
) {
    ifstream input(filename);
    if (!input.is_open() || !input.good()) {
        cout << boost::format("Error: Cannot open --pairs_file %s\n") % filename;
        exit(-1);
    }

    out.clear();
    string line;
    long line_number = 0;
    while (getline(input, line)) {
        line_number++;
        boost::algorithm::trim(line);
        if (line.empty() || line[0] == '#') {
            continue;
        }
        std::istringstream stream(line);
        long i = -1, j = -1;
        if (!(stream >> i >> j)) {
            cout << boost::format("Error: --pairs_file %s line %ld is not two integers.\n")
                    % filename % line_number;
            exit(-1);
        }
        if (i < 0 || j < 0 || i >= n_haplotypes || j >= n_haplotypes) {
            cout << boost::format(
                "Error: --pairs_file %s line %ld references haplotype outside [0, %d).\n"
            ) % filename % line_number % n_haplotypes;
            exit(-1);
        }
        if (i == j) {
            cout << boost::format("Error: --pairs_file %s line %ld pairs a haplotype with itself.\n")
                    % filename % line_number;
            exit(-1);
        }
        out.push_back(make_pair((int) std::min(i, j), (int) std::max(i, j)));
    }

    if (out.empty()) {
        cout << boost::format("Error: --pairs_file %s contains no pairs.\n") % filename;
        exit(-1);
    }
}
