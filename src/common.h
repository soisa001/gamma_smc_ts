#pragma once

#include <iostream>
#include <fstream>
#include <utility>
#include <stdexcept>
#include <exception>
#include <memory>
#include <chrono>
#include <algorithm>
#include <unordered_map>
#include <filesystem>
#include <cstring>
#include <limits>
#include <vector>
#include <cmath>

#include <boost/assert.hpp>
#include <boost/format.hpp>

#include <math.h>
#include <stdlib.h>
#include <stdint.h>

#include <zstd.h>

// TODO: AVX DEFINES
#include <immintrin.h>
#include <emmintrin.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#if defined(__linux__)
#include <sys/mman.h>
#endif

#include "indicators.h"
using namespace indicators;

using namespace std;

typedef long position_t;


// TODO: Unclear if this is needed anymore so much. Now that each segment is
// a stretch of missing and of hom, the difference is just the emission.
enum segment_type : int8_t {
    MISSING_STRETCH_MISSING_SITE, 
    HOM_STRETCH_HOM_SITE, 
    HOM_STRETCH_HET_SITE,
    HOM_STRETCH                 // with missing site in the end; TODO: Change name to HOM_STRETCH_MISSING_SITE?
};

// TODO: Maybe unite those two, they are not so different
struct SegSite_t {
    position_t pos;
    segment_type type;  
    vector<int32_t> alleles;
};

struct SegregatingSite {
    position_t pos;
    vector<int8_t> alleles;
};

// Bit-packed genotype matrix for the retained biallelic segregating sites.
//
// One byte per haplotype per site (the old vector<int8_t> alleles) is 60 GB for
// 3M chr1 sites x 20,000 haplotypes, which does not fit the 64 GB target. Two
// bit planes bring the allele plane to 7.5 GB, and the missing plane is stored
// only for the sites that actually carry a missing genotype -- usually none in
// phased panel data.
struct SiteMatrix {
    // Allele rows live in fixed-size blocks rather than one vector, because a
    // single growing vector would have to copy a multi-GB buffer on each
    // reallocation and would peak at ~2.5x the final footprint.
    static const long sites_per_block = 4096;

    long n_sites = 0;
    int n_haplotypes = 0;
    long words_per_site = 0;               // ceil(n_haplotypes / 64)
    position_t sequence_length = 0;         // bases, from the VCF contig header
    string contig_name;

    vector<position_t> pos;                // n_sites
    vector<vector<uint64_t>> alt_blocks;   // bit set = ALT allele
    vector<int32_t> missing_row;           // n_sites; -1 when the site has no missing call
    vector<uint64_t> missing;              // n_missing_rows * words_per_site; bit set = missing

    inline bool empty() const { return n_sites == 0; }
    inline long size() const { return n_sites; }

    void init(int n_haplotypes_) {
        n_haplotypes = n_haplotypes_;
        words_per_site = (n_haplotypes_ + 63) / 64;
        n_sites = 0;
        sequence_length = 0;
        contig_name.clear();
        pos.clear();
        alt_blocks.clear();
        missing_row.clear();
        missing.clear();
    }

    // Appends a zeroed allele row and returns it for the caller to fill.
    uint64_t* append_site(position_t position) {
        const long block = n_sites / sites_per_block;
        if ((long) alt_blocks.size() <= block) {
            alt_blocks.emplace_back();
            alt_blocks.back().assign((size_t) (sites_per_block * words_per_site), 0ULL);
        }
        uint64_t* row = alt_blocks[block].data() + (n_sites % sites_per_block) * words_per_site;
        pos.push_back(position);
        missing_row.push_back(-1);
        n_sites++;
        return row;
    }

    // Attaches a zeroed missing plane to the most recently appended site.
    uint64_t* append_missing_row() {
        const long row_index = (long) (missing.size() / (size_t) words_per_site);
        missing.resize(missing.size() + (size_t) words_per_site, 0ULL);
        missing_row[n_sites - 1] = (int32_t) row_index;
        return missing.data() + row_index * words_per_site;
    }

    inline const uint64_t* alt_row(long site) const {
        return alt_blocks[site / sites_per_block].data()
             + (site % sites_per_block) * words_per_site;
    }

    // NULL when this site has no missing genotype, which lets the emission loop
    // skip the missing test entirely for the overwhelmingly common case.
    inline const uint64_t* missing_row_ptr(long site) const {
        const int32_t row = missing_row[site];
        return (row < 0) ? NULL : (missing.data() + ((long) row) * words_per_site);
    }

    static inline bool test_bit(const uint64_t* row, int haplotype) {
        return (row[haplotype >> 6] >> (haplotype & 63)) & 1ULL;
    }

    static inline void set_bit(uint64_t* row, int haplotype) {
        row[haplotype >> 6] |= (1ULL << (haplotype & 63));
    }

    // 0/1 allele; only meaningful when the haplotype is not missing.
    inline int allele(long site, int haplotype) const {
        return (int) test_bit(alt_row(site), haplotype);
    }

    inline bool is_missing(long site, int haplotype) const {
        const uint64_t* row = missing_row_ptr(site);
        return (row != NULL) && test_bit(row, haplotype);
    }

    double bytes() const {
        double total = (double) (
            pos.size() * sizeof(position_t) +
            missing_row.size() * sizeof(int32_t) +
            missing.size() * sizeof(uint64_t)
        );
        for (const auto& block : alt_blocks) {
            total += (double) (block.size() * sizeof(uint64_t));
        }
        return total;
    }
};


struct Segment_t {
    position_t pos;   // TODO: Make this start_pos and end_pos, or at least clarify which is it
    int length;    
    segment_type type;
    bool output_at_start;
    bool output_at_end;
    int seg_site_index;
};

// Memory and AVX

const int parallel_vector_size = 8;
const size_t parallel_vector_size_in_bytes_float = parallel_vector_size * sizeof(float);
const size_t memory_alignment = 64;

// aligned_alloc requires the size to be a multiple of the alignment, so round up.
inline void* aligned_alloc_bytes(size_t n_bytes, bool reset) {
    size_t alloc_size = ((n_bytes + memory_alignment - 1) / memory_alignment) * memory_alignment;
    if (alloc_size == 0) {
        alloc_size = memory_alignment;
    }
    void* ptr = aligned_alloc(memory_alignment, alloc_size);
    if (ptr == NULL) {
        cerr << boost::format("Error: failed to allocate %.3f GB.\n")
                % (alloc_size / 1024.0 / 1024.0 / 1024.0);
        exit(-1);
    }
    if (reset) memset(ptr, 0, alloc_size);
    return ptr;
}

// The flow-field cache is a few hundred MB read with an effectively random
// access pattern, so 4 KB pages make it TLB-bound once every thread is hitting
// it. Ask for transparent huge pages; silently ignore a refusal.
inline void advise_huge_pages(void* ptr, size_t n_bytes) {
#if defined(__linux__) && defined(MADV_HUGEPAGE)
    const uintptr_t page = 2u * 1024u * 1024u;
    uintptr_t start = ((uintptr_t) ptr + page - 1) & ~(page - 1);
    uintptr_t end = ((uintptr_t) ptr + n_bytes) & ~(page - 1);
    if (end > start) {
        madvise((void*) start, (size_t) (end - start), MADV_HUGEPAGE);
    }
#else
    (void) ptr;
    (void) n_bytes;
#endif
}

inline float* aligned_alloc_float(size_t n_elements, bool reset = true) {
    return (float*) aligned_alloc_bytes(n_elements * sizeof(float), reset);
}


inline int32_t* aligned_alloc_int32(size_t n_elements, bool reset = true) {
    return (int32_t*) aligned_alloc_bytes(n_elements * sizeof(int32_t), reset);
}


inline int8_t* aligned_alloc_int8(size_t n_elements, bool reset = true) {
    return (int8_t*) aligned_alloc_bytes(n_elements * sizeof(int8_t), reset);
}
