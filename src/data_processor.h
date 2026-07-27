#pragma once

#include "common.h"

class DataProcessor {
  public:

    const SiteMatrix& _sites;
    const vector<string>& _sample_names;
    const vector<pair<int, int>>& _global_mask;
    unordered_map<string, vector<pair<int, int>>>& _mask_map;
    int _posterior_every;
    int _flow_field_cache_n_steps;
    bool _output_at_every;
    bool _output_at_hets;

    bool _is_global_mask;
    vector<position_t> _output_positions;
    vector<Segment_t> _segments;

    // Which sites fall inside a mask, i.e. which sites are *called*. This is
    // the inverse of the old _is_seg_site_missing.
    //
    // With one global mask the answer does not depend on the sample, so a single
    // bit per site replaces the old n_samples x n_sites bool matrix: 375 KB
    // instead of 3.75 GB for 10,000 samples over 3M chr1 sites, and the emission
    // loop stops chasing one cache line per lane.
    vector<uint64_t> _site_called_global;

    // Per-sample masks keep the sample dimension, but site-major so that the
    // eight samples touched by one AVX chunk sit in adjacent bits.
    vector<uint64_t> _site_sample_called;
    long _sample_words_per_site = 0;

    // Index of the segment that ends at each segregating site. Every site ends
    // exactly one segment, so this lets the emission pass visit only the
    // segments whose type actually varies by pair.
    vector<int32_t> _segment_index_of_site;

    long _n_segments;
    position_t _seq_length;

    DataProcessor(
        const SiteMatrix& sites,
        const vector<string>& sample_names,
        const vector<pair<int, int>>& global_mask,
        unordered_map<string, vector<pair<int, int>>>& mask_map,
        int posterior_every,
        int flow_field_cache_n_steps,
        bool output_at_hets
    ) :
        _sites(sites),
        _sample_names(sample_names),
        _global_mask(global_mask),
        _mask_map(mask_map),
        _posterior_every(posterior_every),
        _flow_field_cache_n_steps(flow_field_cache_n_steps),
        _output_at_every(_posterior_every > 0),
        _output_at_hets(output_at_hets),
        _is_global_mask(_mask_map.size() == 0)
    {
        prepare_segments();
        precalculate_seg_sites_in_masks();
    }

    bool is_global_mask() {
        return _is_global_mask;
    }

    // Segment layout is unchanged from upstream; the only addition is
    // _segment_index_of_site, which lets the emission pass touch only the
    // segments whose type depends on the pair instead of walking all of them.
    void prepare_segments() {
        bool prev_to_output = false;
        bool to_output;
        position_t n_pos = -1;
        _seq_length = 0;

        position_t last = _sites.pos.back();
        position_t jump_to_pos;
        position_t next_output = 0;
        long next_site_index = 0;

        _segment_index_of_site.assign((size_t) _sites.size(), -1);

        while (n_pos < last) {
            jump_to_pos = _sites.pos[next_site_index];
            if (_output_at_every) {
                jump_to_pos = std::min(jump_to_pos, next_output);
            }

            // While we can't jump to the next site within the cache
            while (jump_to_pos - n_pos > _flow_field_cache_n_steps) {
                _segments.push_back({
                    n_pos + _flow_field_cache_n_steps,      // pos
                    _flow_field_cache_n_steps,              // segment length
                    HOM_STRETCH,                            // segment type - meaningless default here; TODO: Remove field?
                    prev_to_output,                         // output at start
                    false,                                  // output at end
                    -1,                                     // No corresponding site
                });

                n_pos += _flow_field_cache_n_steps;
                prev_to_output = false;
            }

            // output only if: (i) we output hets; or (ii) it's an output site
            to_output = (_output_at_every && (jump_to_pos == next_output)) || _output_at_hets;

            const bool ends_at_site = (jump_to_pos == _sites.pos[next_site_index]);

            // Now we can jump
            _segments.push_back({
                jump_to_pos,                       // pos
                (int) (jump_to_pos - n_pos),  // segment length
                HOM_STRETCH,                     // segment type - meaningless default here; TODO: Remove field?
                prev_to_output,     // output at start
                to_output,          // output at end
                (ends_at_site ? (int) next_site_index : -1)
            });

            if (to_output) { _output_positions.push_back(jump_to_pos); }

            // Update pointers
            n_pos = jump_to_pos;
            prev_to_output = to_output;

            if (ends_at_site) {
                _segment_index_of_site[next_site_index] = (int32_t) (_segments.size() - 1);
                next_site_index++;
            }
            if (_output_at_every && (jump_to_pos == next_output)) {
                next_output += _posterior_every;
            }
        }

        _n_segments = (long) _segments.size();
        _seq_length = (position_t) _output_positions.size();

        if (_n_segments > (long) std::numeric_limits<int32_t>::max()) {
            cout << boost::format(
                "Error: %ld segments exceeds the 2^31 limit; raise --cache_size.\n"
            ) % _n_segments;
            exit(-1);
        }
    }

    // Marks which sites are callable, from the global mask or per-sample masks.
    void precalculate_seg_sites_in_masks() {
        const long n_sites = _sites.size();
        const long site_words = (n_sites + 63) / 64;

        if (is_global_mask()) {
            _site_called_global.assign((size_t) site_words, 0ULL);
            uint64_t* called = _site_called_global.data();
            size_t n_interval = 0;
            for (long n_seg_site = 0; n_seg_site < n_sites; n_seg_site++) {
                const position_t pos = _sites.pos[n_seg_site];
                while ((n_interval < _global_mask.size()) && (_global_mask[n_interval].first <= pos)) {
                    if (_global_mask[n_interval].second > pos) {
                        called[n_seg_site >> 6] |= (1ULL << (n_seg_site & 63));
                        break;
                    }
                    n_interval++;
                }
            }
            return;
        }

        _sample_words_per_site = (long) ((_sample_names.size() + 63) / 64);
        _site_sample_called.assign((size_t) (n_sites * _sample_words_per_site), 0ULL);

        for (size_t n_sample = 0; n_sample < _sample_names.size(); n_sample++) {
            const vector<pair<int, int>>& mask = _mask_map.at(_sample_names[n_sample]);
            const uint64_t sample_bit = 1ULL << (n_sample & 63);
            const long sample_word = (long) (n_sample >> 6);
            size_t n_interval = 0;
            for (long n_seg_site = 0; n_seg_site < n_sites; n_seg_site++) {
                const position_t pos = _sites.pos[n_seg_site];
                while ((n_interval < mask.size()) && (mask[n_interval].first <= pos)) {
                    if (mask[n_interval].second > pos) {
                        _site_sample_called[(size_t) (n_seg_site * _sample_words_per_site + sample_word)] |= sample_bit;
                        break;
                    }
                    n_interval++;
                }
            }
        }
    }

    inline bool site_is_called_global(long site) const {
        return (_site_called_global[(size_t) (site >> 6)] >> (site & 63)) & 1ULL;
    }

    inline bool site_is_called_for_sample(long site, int sample) const {
        return (_site_sample_called[(size_t) (site * _sample_words_per_site + (sample >> 6))]
                >> (sample & 63)) & 1ULL;
    }

    float calculate_heterozygosity_for_sample(int n_sample) {
        const vector<pair<int, int>>* mask;
            if (is_global_mask()) {
                mask = &_global_mask;
            } else {
                mask = &(_mask_map.at(_sample_names[n_sample]));
            };

        // Count hets not missing and not in mask
        long n_hets = 0;
        const int hap_a = 2 * n_sample;
        const int hap_b = 2 * n_sample + 1;
        for (long n_seg_site = 0; n_seg_site < _sites.size(); n_seg_site++) {
            // If missing, skip
            const uint64_t* missing = _sites.missing_row_ptr(n_seg_site);
            if ((missing != NULL)
                && (SiteMatrix::test_bit(missing, hap_a) || SiteMatrix::test_bit(missing, hap_b))) {
                continue;
            }

            // If hom, skip
            const uint64_t* alt = _sites.alt_row(n_seg_site);
            if (SiteMatrix::test_bit(alt, hap_a) == SiteMatrix::test_bit(alt, hap_b)) {
                continue;
            }

            // If masked, skip
            const bool called = is_global_mask()
                ? site_is_called_global(n_seg_site)
                : site_is_called_for_sample(n_seg_site, n_sample);
            if (!called) {
                continue;
            }

            // Otherwise, count
            n_hets++;
        }

        // Count number of sites in mask
        double n_total_sites = 0;
        for (auto& interval : *mask) {
            n_total_sites += (double) (interval.second-interval.first);
        }

        return (float) (((double) n_hets) / n_total_sites);
    }

    float calculate_heterozygosity() {
        double theta = 0.0;
        const long n_samples = (long) _sample_names.size();
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : theta) schedule(static)
#endif
        for (long n_sample = 0; n_sample < n_samples; n_sample++) {
            theta += (double) calculate_heterozygosity_for_sample((int) n_sample);
        }
        return (float) (theta / n_samples);
    }

    void intersect_masks(
        const vector<pair<int, int>>& mask1,
        const vector<pair<int, int>>& mask2,
        int32_t* n_called_array,
        long ptr_jump = 1            // stride into n_called_array
        ) const {
        size_t i = 0, j = 0;
        long k = 0;
        size_t n = mask1.size(), m = mask2.size();

        position_t cur_segment_start = 0;
        position_t cur_segment_end = _segments[0].pos;

        int32_t n_called = 0;

        // Loop through all intervals unless one of the interval gets exhausted
        while ((i < n) && (j < m) && (k < _n_segments)) {
            // Left bound for intersecting interval
            position_t l = max(max((position_t) mask1[i].first, (position_t) mask2[j].first), cur_segment_start);

            // Right bound for intersecting interval
            position_t r = min(min((position_t) mask1[i].second, (position_t) mask2[j].second), cur_segment_end);

            // If interval is valid, count it
            if (l < r) {
                n_called += (int32_t) (r - l);
            }

            if (cur_segment_end < mask1[i].second) {
                if (cur_segment_end < mask2[j].second) {
                    // Belt and braces against a mask that still manages to
                    // double-count: n_called above the segment length would make
                    // n_missing negative in at_flat_vectorized and index the
                    // cached tables from before their base.
                    const int32_t span = (int32_t) (cur_segment_end - cur_segment_start);
                    *n_called_array = min(n_called, span);
                    n_called_array += ptr_jump;
                    n_called = 0;

                    k++;
                    if (k >= _n_segments) {
                        break;
                    }
                    cur_segment_start = cur_segment_end;
                    cur_segment_end = _segments[k].pos;
                } else {
                    j++;
                }
            } else {
                if (mask1[i].second < mask2[j].second) {
                    i++;
                } else {
                    j++;
                }
            }
        }

        while (k < _n_segments) {
            *n_called_array = 0;
            n_called_array += ptr_jump;
            k++;
        }

    }

    void intersect_global_mask(
        int32_t* n_called_array,
        long ptr_jump = 1            // stride into n_called_array
        ) const {
        intersect_masks(
            _global_mask,
            _global_mask,
            n_called_array,
            ptr_jump
        );
    }

    void intersect_masks_at_ids(
        int sample_id_1,
        int sample_id_2,
        int32_t* n_called_array,
        long ptr_jump = 1            // stride into n_called_array
        ) const {
        auto& mask1 = _mask_map.at(_sample_names[sample_id_1]);
        auto& mask2 = _mask_map.at(_sample_names[sample_id_2]);

        intersect_masks(
            mask1,
            mask2,
            n_called_array,
            ptr_jump
        );
    }


};
