#pragma once

#include "common.h"
#include "flow_field.h"
#include "data_processor.h"
#include "recent_stats.h"
#include "bitmatrix.h"

#include <atomic>
#include <iomanip>
#include <sstream>

// https://gist.github.com/andersx/8057b2a6fd3d715d35eb

// Approximation for EXP(x) -- very fast, but not super accurate
static inline __m256 _mm256_expfaster_ps(const __m256 &q) {

     // REGEV: Added this to make it 10**x instead of exp(x)
    const __m256 ln10 =  _mm256_set1_ps (2.30258509299f);

    const __m256 C1 = _mm256_set1_ps(1064872507.1541044f);
    const __m256 C2 = _mm256_set1_ps(12102203.161561485f);

    return _mm256_castsi256_ps(_mm256_cvttps_epi32(_mm256_fmadd_ps(C2, _mm256_mul_ps(q, ln10), C1)));
}

// Accurate 10^q. This is the default; --exp10 fast selects the approximation
// above.
//
// The Schraudolph trick replaces 2^f by the straight line 1+f inside each
// binade, giving -3.89% .. +2.01% relative error. Measured over the reachable
// range that is a near zero-mean sawtooth (mean +0.03%), not a systematic bias
// -- but it is a deterministic function of the value rather than noise, so it
// does not cancel across pairs whose posteriors land in the same part of the
// sawtooth, and alpha and beta are perturbed independently of one another. It
// lands directly on the posterior shape and rate and therefore on P(T < t).
//
// This version splits 10^q = 2^n * 2^f with f in [-0.5, 0.5] and evaluates 2^f
// with its degree-6 Taylor series, whose truncation term at |f| = 0.5 is 1.19e-7
// against a float32 epsilon of 1.19e-7 -- one ULP, as good as single precision
// allows. End to end the routine is accurate to 1.19e-6, about 10 ULP; the
// residual is the float32 range reduction q * log2(10), not the polynomial.
// That is 40x tighter than the P(T<t) lookup table consuming it, so splitting
// log2(10) Cody-Waite style to recover the last few ULP would buy nothing.
//
// It is called once per output position, not once per segment, so the extra
// arithmetic over the fast path measures +2.4% of decode.
//
// PRECONDITION: neither this nor the fast path above clamps the exponent. Here,
// (n + 127) << 23 walks into the sign bit once n <= -127, returning large
// wrong-signed values rather than underflowing to zero -- multiplying by the
// constructed power of two does not rescue that. What keeps both safe is the
// caller: clip_mean_cv_log10 holds the message state inside the flow-field grid
// (mean in [1e-5, 100], cv in [0.01, 1]), which bounds q to [-2, 9] and n to
// [-7, 30]. Widen that grid and this needs an explicit clamp.
static inline __m256 _mm256_exp10_accurate_ps(const __m256 &q) {
    const __m256 log2_10 = _mm256_set1_ps(3.32192809488736f);
    const __m256 y = _mm256_mul_ps(q, log2_10);

    const __m256 n = _mm256_round_ps(y, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
    const __m256 f = _mm256_sub_ps(y, n);

    // 2^f = sum_k (ln2)^k f^k / k!
    __m256 poly = _mm256_set1_ps(1.5403530393381609e-4f);
    poly = _mm256_fmadd_ps(poly, f, _mm256_set1_ps(1.3333558146428441e-3f));
    poly = _mm256_fmadd_ps(poly, f, _mm256_set1_ps(9.6181291076284772e-3f));
    poly = _mm256_fmadd_ps(poly, f, _mm256_set1_ps(5.5504108664821579e-2f));
    poly = _mm256_fmadd_ps(poly, f, _mm256_set1_ps(2.4022650695910071e-1f));
    poly = _mm256_fmadd_ps(poly, f, _mm256_set1_ps(6.9314718055994531e-1f));
    poly = _mm256_fmadd_ps(poly, f, _mm256_set1_ps(1.0f));

    const __m256i exponent = _mm256_slli_epi32(
        _mm256_add_epi32(_mm256_cvtps_epi32(n), _mm256_set1_epi32(127)), 23
    );
    return _mm256_mul_ps(poly, _mm256_castsi256_ps(exponent));
}

// --------------------------------------------------------------------------------------------------

// Everything one worker thread mutates while it decodes a block of pairs.
//
// The forward/backward passes are embarrassingly parallel across pairs, so the
// only thing standing between the original single-threaded loop and 32 cores
// was that this scratch used to live in class members. Pulling it into a
// per-thread struct is the whole parallelisation.
//
// Buffers touched by AVX aligned loads must come from aligned_alloc_*; a
// std::vector<float> only guarantees alignof(float) and would segfault.
struct PairWorkspace {
    float* posteriors_alpha = NULL;      // seq_length * 8
    float* posteriors_beta = NULL;
    float* mean_log10 = NULL;            // 8
    float* cv_log10 = NULL;              // 8
    int8_t* segment_types = NULL;        // n_segments * 8
    int32_t* n_called = NULL;            // n_segments * 8; shared unless per-sample masks
    bool owns_n_called = false;

    // Per-position partial sums, merged once at the end of the run.
    vector<double> sum_probability;      // n_thresholds * seq_length
    vector<double> sum_tmrca;            // seq_length
    vector<int64_t> n_valid;             // seq_length
    vector<int64_t> n_recent;            // n_thresholds * seq_length

    // Packed calls for the block in flight, plus its compression scratch.
    vector<uint8_t> bits;
    vector<char> compressed;
    ZSTD_CCtx* cctx = NULL;

    double timer_emissions = 0.0;
    double timer_forward = 0.0;
    double timer_backward = 0.0;
    double timer_output = 0.0;

    PairWorkspace() = default;
    PairWorkspace(const PairWorkspace&) = delete;
    PairWorkspace& operator=(const PairWorkspace&) = delete;

    ~PairWorkspace() {
        free(posteriors_alpha);
        free(posteriors_beta);
        free(mean_log10);
        free(cv_log10);
        free(segment_types);
        if (owns_n_called) {
            free(n_called);
        }
        if (cctx != NULL) {
            ZSTD_freeCCtx(cctx);
        }
    }
};

class CachedPairwiseGammaSMC {
  public:
    const SiteMatrix& _sites;
    const vector<pair<int, int>>& _haplotype_pairs;
    float _scaled_recombination_rate;
    float _scaled_mutation_rate;
    unique_ptr<FlowFieldCache> _flow_field_cache;
    DataProcessor& _data_processor;
    int _posterior_every;
    int _flow_field_cache_n_steps;

    long _n_pairs;
    long _n_pairs_in_chunk;

    bool _output_at_hets;
    bool _only_forward;
    bool _only_backward;

    const vector<position_t>& _output_positions;
    const vector<Segment_t>& _segments;

    long _n_segments;
    position_t _seq_length;

    ofstream* _output_file_raw_header;
    ostream* _output_file_raw;
    ostream* _recent_summary_output;

    // Recent-coalescence configuration and outputs
    vector<RecentThreshold> _thresholds;
    const RecentProbabilityTable* _probability_table = NULL;
    recent_call_t _recent_call = RECENT_CALL_MEDIAN;
    bool _accumulate_probability = true;
    double _two_ne_generations = -1.0;
    vector<double> _sum_probability;
    vector<double> _sum_tmrca_generations;
    vector<int64_t> _n_summary_pairs;
    vector<int64_t> _n_recent_pairs;

    // Packed bit-matrix output
    ostream* _bitmatrix_output = NULL;
    string _bitmatrix_filename;
    unique_ptr<BitMatrixWriter> _bitmatrix_writer;

    int _zstd_compression_level;

    // Numerical/algorithmic switches. Both default to the CORRECTED behaviour;
    // --exp10 fast and --backward_alignment legacy reproduce upstream, for
    // comparing against results generated with the old binary. Set from main()
    // after construction.
    bool _accurate_exp10 = true;
    bool _fix_backward_alignment = true;
    bool _exact_recent_stats = false;
    double _recent_call_probability = 0.5;

    // Blocking and threading
    long _pair_block = 256;
    long _chunks_per_block = 32;
    long _n_blocks = 0;
    int _n_threads = 1;

    // Shared, read-only for the whole run when a single global mask applies.
    int32_t* _global_n_called = NULL;

    double _timer_emissions = 0.0;
    double _timer_forward = 0.0;
    double _timer_backward = 0.0;
    double _timer_output = 0.0;
    double _timer_wall = 0.0;

    CachedPairwiseGammaSMC(
        const SiteMatrix& sites,
        const vector<pair<int, int>>& haplotype_pairs,
        float scaled_recombination_rate,
        float scaled_mutation_rate,
        unique_ptr<FlowFieldCache> flow_field_cache,  // TODO: Do we really need unique_ptr rather than const&
        DataProcessor& data_processor,
        int posterior_every,
        bool output_at_hets,
        bool only_forward,
        bool only_backward,
        ofstream* output_file_raw_header,
        ostream* output_file_raw,
        int zstd_compression_level,
        ostream* recent_summary_output,
        const vector<RecentThreshold>& thresholds,
        const RecentProbabilityTable* probability_table,
        recent_call_t recent_call,
        bool accumulate_probability,
        double two_ne_generations,
        ostream* bitmatrix_output,
        const string& bitmatrix_filename,
        long pair_block,
        int n_threads
    ) :
        _sites(sites),
        _haplotype_pairs(haplotype_pairs),
        _scaled_recombination_rate(scaled_recombination_rate),
        _scaled_mutation_rate(scaled_mutation_rate),
        _flow_field_cache(move(flow_field_cache)),
        _data_processor(data_processor),
        _posterior_every(posterior_every),
        _flow_field_cache_n_steps(_flow_field_cache->_n_steps),
        _n_pairs(_haplotype_pairs.size()),
        _n_pairs_in_chunk(parallel_vector_size),
        _output_at_hets(output_at_hets),
        _only_forward(only_forward),
        _only_backward(only_backward),
        _output_positions(data_processor._output_positions),
        _segments(data_processor._segments),
        _n_segments(data_processor._n_segments),
        _seq_length(data_processor._seq_length),
        _output_file_raw_header(output_file_raw_header),
        _output_file_raw(output_file_raw),
        _recent_summary_output(recent_summary_output),
        _thresholds(thresholds),
        _probability_table(probability_table),
        _recent_call(recent_call),
        _accumulate_probability(accumulate_probability),
        _two_ne_generations(two_ne_generations),
        _bitmatrix_output(bitmatrix_output),
        _bitmatrix_filename(bitmatrix_filename),
        _zstd_compression_level(zstd_compression_level),
        _n_threads(n_threads)
    {
        // The raw posterior stream is one continuous zstd frame whose reader
        // reconstructs chunk k from byte offset k, so it can only be written in
        // order. Keeping blocks at one chunk lets an ordered section emit it
        // without buffering a block's worth of posteriors per thread.
        _pair_block = (_output_file_raw != NULL) ? _n_pairs_in_chunk : pair_block;
        if (_pair_block < _n_pairs_in_chunk) {
            _pair_block = _n_pairs_in_chunk;
        }
        _pair_block -= (_pair_block % _n_pairs_in_chunk);
        _chunks_per_block = _pair_block / _n_pairs_in_chunk;
        _n_blocks = (_n_pairs + _pair_block - 1) / _pair_block;

        if (_data_processor.is_global_mask()) {
            // Identical for every pair and every chunk, so one shared copy.
            _global_n_called = aligned_alloc_int32((size_t) (_n_segments * _n_pairs_in_chunk), true);
            _data_processor.intersect_global_mask(_global_n_called, _n_pairs_in_chunk);
            // Broadcast lane 0 across the other seven instead of intersecting
            // the same two masks eight times.
            for (long segment = 0; segment < _n_segments; segment++) {
                int32_t* row = _global_n_called + segment * _n_pairs_in_chunk;
                for (long lane = 1; lane < _n_pairs_in_chunk; lane++) {
                    row[lane] = row[0];
                }
            }
        }

        if (_recent_summary_output != NULL || _bitmatrix_output != NULL) {
            const size_t per_threshold = (size_t) (_thresholds.size() * _seq_length);
            _sum_probability.assign(per_threshold, 0.0);
            _n_recent_pairs.assign(per_threshold, 0);
            _sum_tmrca_generations.assign((size_t) _seq_length, 0.0);
            _n_summary_pairs.assign((size_t) _seq_length, 0);
        }

        if (_bitmatrix_output != NULL) {
            _bitmatrix_writer.reset(new BitMatrixWriter(_bitmatrix_output));
        }
    }

    // ----------------------------------------------------------------------
    // Per-thread workspace
    // ----------------------------------------------------------------------

    void allocate_workspace(PairWorkspace& ws) const {
        const size_t posterior_elements = (size_t) (_seq_length * _n_pairs_in_chunk);
        ws.posteriors_alpha = aligned_alloc_float(posterior_elements, true);
        ws.posteriors_beta = aligned_alloc_float(posterior_elements, true);
        ws.mean_log10 = aligned_alloc_float((size_t) _n_pairs_in_chunk, true);
        ws.cv_log10 = aligned_alloc_float((size_t) _n_pairs_in_chunk, true);

        // Segments that do not end at a segregating site always emit
        // HOM_STRETCH_HOM_SITE, whatever the pair. Filling them once here
        // instead of once per chunk removes the dominant term of the emission
        // pass: there are millions of such segments and only one pass over the
        // segregating sites is actually pair-dependent.
        ws.segment_types = aligned_alloc_int8((size_t) (_n_segments * _n_pairs_in_chunk), false);
        memset((void*) ws.segment_types, HOM_STRETCH_HOM_SITE,
               (size_t) (_n_segments * _n_pairs_in_chunk));

        if (_data_processor.is_global_mask()) {
            ws.n_called = _global_n_called;
            ws.owns_n_called = false;
        } else {
            ws.n_called = aligned_alloc_int32((size_t) (_n_segments * _n_pairs_in_chunk), true);
            ws.owns_n_called = true;
        }

        if (needs_statistics()) {
            const size_t per_threshold = (size_t) (_thresholds.size() * _seq_length);
            ws.sum_probability.assign(per_threshold, 0.0);
            ws.n_recent.assign(per_threshold, 0);
            ws.sum_tmrca.assign((size_t) _seq_length, 0.0);
            ws.n_valid.assign((size_t) _seq_length, 0);
        }

        if (_bitmatrix_output != NULL) {
            const size_t payload = bitmatrix_payload_bytes(_chunks_per_block);
            ws.bits.assign(payload, 0);
            ws.compressed.resize(ZSTD_compressBound(payload));
            ws.cctx = ZSTD_createCCtx();
            ZSTD_CCtx_setParameter(ws.cctx, ZSTD_c_compressionLevel, _zstd_compression_level);
        }
    }

    inline bool needs_statistics() const {
        return (_recent_summary_output != NULL) || (_bitmatrix_output != NULL);
    }

    inline size_t bitmatrix_payload_bytes(long chunks_in_block) const {
        return (size_t) _thresholds.size() * (size_t) chunks_in_block * (size_t) _seq_length;
    }

    // ----------------------------------------------------------------------
    // Emissions
    // ----------------------------------------------------------------------

    void prepare_pairwise_n_called(PairWorkspace& ws, long starting_n_pair) const {
        int i, j;
        for (long n_pair = starting_n_pair;
             n_pair < min(_n_pairs, starting_n_pair + _n_pairs_in_chunk); n_pair++) {
            tie(i, j) = _haplotype_pairs[n_pair];
            _data_processor.intersect_masks_at_ids(
                i >> 1,         // First sample ID
                j >> 1,         // Second sample ID
                ws.n_called + (n_pair - starting_n_pair),      // Pointer to place in big array
                _n_pairs_in_chunk             // Stride into array
            );
        }
    }

    // Only the segments that end at a segregating site depend on the pair; the
    // rest were filled once in allocate_workspace().
    void prepare_pairwise_emissions(PairWorkspace& ws, long starting_n_pair) const {
        const bool global_mask = _data_processor._is_global_mask;
        const long n_sites = _sites.size();
        const long last_pair = min(_n_pairs, starting_n_pair + _n_pairs_in_chunk);
        const long real_lanes = last_pair - starting_n_pair;

        for (long site = 0; site < n_sites; site++) {
            const int32_t segment_index = _data_processor._segment_index_of_site[site];
            if (segment_index < 0) {
                continue;
            }
            int8_t* cur_ptr = ws.segment_types + _n_pairs_in_chunk * (long) segment_index;

            const uint64_t* alt = _sites.alt_row(site);
            const uint64_t* missing = _sites.missing_row_ptr(site);
            const bool site_called_globally = global_mask
                ? _data_processor.site_is_called_global(site)
                : true;

            if (global_mask && !site_called_globally) {
                // Outside the mask for every sample: the site contributes no
                // emission for any pair in this chunk.
                memset((void*) cur_ptr, HOM_STRETCH, (size_t) real_lanes);
            } else {
                for (long lane = 0; lane < real_lanes; lane++) {
                    const pair<int, int>& haplotypes = _haplotype_pairs[starting_n_pair + lane];
                    const int i = haplotypes.first;
                    const int j = haplotypes.second;

                    // Missing means either haplotype is uncalled at this site,
                    // or the site is outside a sample's mask. Bit-packed, so
                    // this is the same test the old
                    //   is_missing_i || is_missing_j || alleles[i] < 0 || alleles[j] < 0
                    // performed, one cache line instead of four.
                    bool is_missing =
                        (missing != NULL)
                        && (SiteMatrix::test_bit(missing, i) || SiteMatrix::test_bit(missing, j));
                    if (!global_mask) {
                        is_missing = is_missing
                            || !_data_processor.site_is_called_for_sample(site, i >> 1)
                            || !_data_processor.site_is_called_for_sample(site, j >> 1);
                    }

                    if (is_missing) {
                        cur_ptr[lane] = HOM_STRETCH;
                    } else {
                        cur_ptr[lane] = (SiteMatrix::test_bit(alt, i) != SiteMatrix::test_bit(alt, j))
                            ? HOM_STRETCH_HET_SITE
                            : HOM_STRETCH_HOM_SITE;
                    }
                }
            }

            if (real_lanes < _n_pairs_in_chunk) {
                // Padding lanes past the end of the pair list; value is
                // irrelevant but must not be left stale.
                memset((void*) (cur_ptr + real_lanes), MISSING_STRETCH_MISSING_SITE,
                       (size_t) (_n_pairs_in_chunk - real_lanes));
            }
        }
    }

    // ----------------------------------------------------------------------
    // Forward / backward
    // ----------------------------------------------------------------------

    inline __m256 exp10_vec(__m256 value) const {
        return _accurate_exp10 ? _mm256_exp10_accurate_ps(value) : _mm256_expfaster_ps(value);
    }

    void mean_cv_to_alpha_beta_log10_vec_forward(
        const float* mean_log10,
        const float* cv_log10,
        float* alpha_log10_output,
        float* beta_log10_output
        ) const {

        const __m256 minus_two = _mm256_set1_ps(-2);
        const __m256 alpha_log10_vec = _mm256_mul_ps(minus_two, _mm256_load_ps(cv_log10));
        const __m256 beta_log10_vec = _mm256_sub_ps(alpha_log10_vec, _mm256_load_ps(mean_log10));

        _mm256_store_ps(alpha_log10_output, exp10_vec(alpha_log10_vec));
        _mm256_store_ps(beta_log10_output, exp10_vec(beta_log10_vec));
    }

    void mean_cv_to_alpha_beta_log10_vec_backward(
        const float* mean_log10,
        const float* cv_log10,
        float* alpha_log10_output,
        float* beta_log10_output
        ) const {

        const __m256 minus_two = _mm256_set1_ps(-2);
        const __m256 minus_one = _mm256_set1_ps(-1);

        const __m256 alpha_log10_vec = _mm256_mul_ps(minus_two, _mm256_load_ps(cv_log10));
        const __m256 beta_log10_vec = _mm256_sub_ps(alpha_log10_vec, _mm256_load_ps(mean_log10));

        _mm256_store_ps(
            alpha_log10_output,
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_load_ps(alpha_log10_output),
                    exp10_vec(alpha_log10_vec)
                ),
                minus_one)
        );

        _mm256_store_ps(
            beta_log10_output,
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_load_ps(beta_log10_output),
                    exp10_vec(beta_log10_vec)
                ),
                minus_one)
        );
    }

    void forward_vectorized(PairWorkspace& ws) const {
        memset((void *) ws.mean_log10, 0, _n_pairs_in_chunk * sizeof(float));
        memset((void *) ws.cv_log10, 0, _n_pairs_in_chunk * sizeof(float));

        const segment_type* segment_types_ptr = (const segment_type*) ws.segment_types;
        const int32_t* n_called_ptr = ws.n_called;

        float* posteriors_alpha_ptr = ws.posteriors_alpha;
        float* posteriors_beta_ptr = ws.posteriors_beta;

        for (long next_processed_segment_index = 0;
             next_processed_segment_index < _n_segments; next_processed_segment_index++) {
            const Segment_t& next_segment = _segments[next_processed_segment_index];

            _flow_field_cache->at_flat_vectorized(
                ws.mean_log10,
                ws.cv_log10,
                (int32_t) next_segment.length,
                n_called_ptr,
                segment_types_ptr,
                true           // forward
            );

            segment_types_ptr += _n_pairs_in_chunk;
            n_called_ptr += _n_pairs_in_chunk;

            // If we need to output, do so
            if (next_segment.output_at_end) {
                mean_cv_to_alpha_beta_log10_vec_forward(
                    ws.mean_log10,
                    ws.cv_log10,
                    posteriors_alpha_ptr,
                    posteriors_beta_ptr
                );

                posteriors_alpha_ptr += _n_pairs_in_chunk;
                posteriors_beta_ptr += _n_pairs_in_chunk;
            }
        }
    }

    void backward_vectorized(PairWorkspace& ws) const {
        memset((void *) ws.mean_log10, 0, _n_pairs_in_chunk * sizeof(float));
        memset((void *) ws.cv_log10, 0, _n_pairs_in_chunk * sizeof(float));

        // Segment k writes the backward message at _segments[k-1].pos, because
        // output_at_start[k] == output_at_end[k-1]. There are therefore
        // _seq_length - 1 backward writes whenever the final segment is itself
        // an output position, one fewer than the forward pass makes. Starting
        // the pointer at _seq_length-1 then pairs every backward message with
        // the forward message one output position to its right and leaves
        // position 0 with no backward message at all.
        //
        // Whether that happens depends on the output mode. With
        // --output_at_hets every segment ending at a site is an output
        // position, including the last, so the shift always occurs. With
        // stride-only output the last segment ends at the last segregating
        // site, which is an output position only by coincidence, so the two
        // settings normally agree -- measured identical over 5,000 positions.
        // Leaving the last position with only a forward contribution is
        // correct: the backward message there is the Exp(1) prior, and the
        // alpha_f + alpha_b - 1 convention makes that a no-op.
        const position_t start_index = _fix_backward_alignment
            ? max((position_t) 0, _seq_length - 1 - (_segments.back().output_at_end ? 1 : 0))
            : (_seq_length - 1);

        float* posteriors_alpha_ptr =
            ws.posteriors_alpha + _n_pairs_in_chunk * start_index;
        float* posteriors_beta_ptr =
            ws.posteriors_beta + _n_pairs_in_chunk * start_index;

        for (long next_processed_segment_index = _n_segments-1;
             next_processed_segment_index >= 0; next_processed_segment_index--) {
            const Segment_t& next_segment = _segments[next_processed_segment_index];

            _flow_field_cache->at_flat_vectorized(
                ws.mean_log10,
                ws.cv_log10,
                (int32_t) next_segment.length,
                ws.n_called + _n_pairs_in_chunk * next_processed_segment_index,
                (const segment_type*) (ws.segment_types + _n_pairs_in_chunk * next_processed_segment_index),
                false           // backward
            );

            if (next_segment.output_at_start) {
                mean_cv_to_alpha_beta_log10_vec_backward(
                    ws.mean_log10,
                    ws.cv_log10,
                    posteriors_alpha_ptr,
                    posteriors_beta_ptr
                );

                posteriors_alpha_ptr -= _n_pairs_in_chunk;
                posteriors_beta_ptr -= _n_pairs_in_chunk;
            }
        }
    }

    // ----------------------------------------------------------------------
    // Recent-coalescence calls, packed one bit per pair
    // ----------------------------------------------------------------------

    // Reference implementation: boost::math::gamma_p per element, exactly what
    // upstream did. Far too slow for a real run -- it is the reason the lookup
    // tables exist -- but it isolates the tables' contribution to both runtime
    // and output when the two are run over the same decoded posteriors.
    void accumulate_and_pack_exact(
        PairWorkspace& ws,
        long chunk_first_pair,
        long chunk_in_block,
        long chunks_in_block
    ) const {
        const long real_lanes = min(_n_pairs_in_chunk, _n_pairs - chunk_first_pair);
        const int n_thresholds = (int) _thresholds.size();
        const bool write_bits = (_bitmatrix_output != NULL);
        const double call_probability =
            (_recent_call == RECENT_CALL_MEDIAN) ? 0.5 : _recent_call_probability;

        for (position_t pos = 0; pos < _seq_length; pos++) {
            const long offset = pos * _n_pairs_in_chunk;
            for (long lane = 0; lane < real_lanes; lane++) {
                const double alpha = (double) ws.posteriors_alpha[offset + lane];
                const double beta = (double) ws.posteriors_beta[offset + lane];
                if (!(alpha > 0.0) || !(beta > 0.0)
                    || !std::isfinite(alpha) || !std::isfinite(beta)) {
                    continue;
                }
                ws.n_valid[(size_t) pos] += 1;
                ws.sum_tmrca[(size_t) pos] += (alpha / beta) * _two_ne_generations;

                for (int k = 0; k < n_thresholds; k++) {
                    const double x = beta * _thresholds[k].scaled;
                    const double probability = boost::math::gamma_p(alpha, x);
                    const bool called = (_recent_call == RECENT_CALL_MEAN)
                        ? ((alpha / beta) < _thresholds[k].scaled)
                        : (probability >= call_probability);
                    if (called) {
                        ws.n_recent[(size_t) (k * _seq_length + pos)] += 1;
                        if (write_bits) {
                            ws.bits[(size_t) (((k * chunks_in_block) + chunk_in_block)
                                              * _seq_length + pos)] |= (uint8_t) (1u << lane);
                        }
                    }
                    if (_accumulate_probability) {
                        ws.sum_probability[(size_t) (k * _seq_length + pos)] += probability;
                    }
                }
            }
        }
    }

    void accumulate_and_pack(
        PairWorkspace& ws,
        long chunk_first_pair,
        long chunk_in_block,
        long chunks_in_block
    ) const {
        if (_exact_recent_stats) {
            accumulate_and_pack_exact(ws, chunk_first_pair, chunk_in_block, chunks_in_block);
            return;
        }

        const long real_lanes = min(_n_pairs_in_chunk, _n_pairs - chunk_first_pair);
        // Lanes past the end of the pair list must never contribute a count or
        // set a bit, so mask them out here rather than trusting their contents.
        alignas(32) int32_t lane_bits[parallel_vector_size];
        for (long lane = 0; lane < _n_pairs_in_chunk; lane++) {
            lane_bits[lane] = (lane < real_lanes) ? -1 : 0;
        }
        const __m256i lane_mask = _mm256_load_si256((const __m256i*) lane_bits);

        const int n_thresholds = (int) _thresholds.size();
        const __m256 two_ne = _mm256_set1_ps((float) _two_ne_generations);
        const bool write_bits = (_bitmatrix_output != NULL);

        for (position_t pos = 0; pos < _seq_length; pos++) {
            const long offset = pos * _n_pairs_in_chunk;
            const __m256 alpha = _mm256_load_ps(ws.posteriors_alpha + offset);
            const __m256 beta = _mm256_load_ps(ws.posteriors_beta + offset);

            const __m256i valid_i = _mm256_and_si256(
                lane_mask,
                _mm256_and_si256(positive_finite_mask(alpha), positive_finite_mask(beta))
            );
            const int valid_bits = _mm256_movemask_ps(_mm256_castsi256_ps(valid_i));
            if (valid_bits == 0) {
                continue;
            }
            const __m256 valid = _mm256_castsi256_ps(valid_i);

            ws.n_valid[(size_t) pos] += (int64_t) __builtin_popcount((unsigned) valid_bits);

            // Posterior mean TMRCA in generations; invalid lanes zeroed so they
            // add nothing to the running sum.
            const __m256 tmrca = _mm256_and_ps(
                valid,
                _mm256_mul_ps(_mm256_div_ps(alpha, beta), two_ne)
            );
            ws.sum_tmrca[(size_t) pos] += (double) horizontal_sum(tmrca);

            for (int k = 0; k < n_thresholds; k++) {
                const RecentThreshold& threshold = _thresholds[k];
                const __m256 x = _mm256_mul_ps(beta, _mm256_set1_ps(threshold.scaled_float));

                const __m256 called = _mm256_and_ps(valid, recent_call_mask(threshold, alpha, x));
                const int call_bits = _mm256_movemask_ps(called);
                ws.n_recent[(size_t) (k * _seq_length + pos)] +=
                    (int64_t) __builtin_popcount((unsigned) call_bits);

                if (write_bits) {
                    ws.bits[(size_t) (((k * chunks_in_block) + chunk_in_block) * _seq_length + pos)] =
                        (uint8_t) call_bits;
                }

                if (_accumulate_probability) {
                    const __m256 probability = _mm256_and_ps(
                        valid, recent_probability(*_probability_table, alpha, x)
                    );
                    ws.sum_probability[(size_t) (k * _seq_length + pos)] +=
                        (double) horizontal_sum(probability);
                }
            }
        }
    }

    void merge_workspace(const PairWorkspace& ws) {
        const size_t per_threshold = (size_t) (_thresholds.size() * _seq_length);
        for (size_t i = 0; i < per_threshold; i++) {
            _sum_probability[i] += ws.sum_probability[i];
            _n_recent_pairs[i] += ws.n_recent[i];
        }
        for (position_t i = 0; i < _seq_length; i++) {
            _sum_tmrca_generations[(size_t) i] += ws.sum_tmrca[(size_t) i];
            _n_summary_pairs[(size_t) i] += ws.n_valid[(size_t) i];
        }
        _timer_emissions += ws.timer_emissions;
        _timer_forward += ws.timer_forward;
        _timer_backward += ws.timer_backward;
        _timer_output += ws.timer_output;
    }

    // ----------------------------------------------------------------------
    // Output
    // ----------------------------------------------------------------------

    static string format_threshold(double years) {
        std::ostringstream stream;
        if (years == std::floor(years) && std::fabs(years) < 1e15) {
            stream << (long long) years;
        } else {
            stream << years;
        }
        return stream.str();
    }

    void output_raw_header() {
        ostream& out = *_output_file_raw_header;
        out << "{\n";

        out << boost::format("\t\"scaled_mutation_rate\": %.10f,\n") % _scaled_mutation_rate;
        out << boost::format("\t\"scaled_recombination_rate\": %.10f,\n") % _scaled_recombination_rate;
        out << boost::format("\t\"sequence_length\": %ld,\n") % _seq_length;
        out << boost::format("\t\"chunk_size\": %ld,\n") % _n_pairs_in_chunk;
        out << boost::format("\t\"num_pairs\": %ld,\n") % _n_pairs;

        out << "\t\"sample_names\": {";
        for (size_t i = 0; i < _data_processor._sample_names.size(); i++) {
            const string& sample_name = _data_processor._sample_names[i];
            out << '"' << (i * 2) << "\": \"" << sample_name << ".0\", \""
                << (i * 2 + 1) << "\": \"" << sample_name << ".1\"";
            if (i + 1 != _data_processor._sample_names.size()) {
                out << ", ";
            }
        }
        out << "},\n";

        // boost::format per element costs seconds once there are 1e5 pairs and
        // 2.5e5 positions, so these two lists go straight to the stream.
        out << "\t\"output_positions\": [\n\t\t";
        for (position_t i = 0; i < _seq_length; i++) {
            out << _output_positions[(size_t) i];
            if (i + 1 != _seq_length) {
                out << ", ";
            }
        }
        out << "\n\t],\n";

        out << "\t\"pairs\": [\n\t\t";
        for (long i = 0; i < _n_pairs; i++) {
            out << '[' << _haplotype_pairs[(size_t) i].first << ", "
                << _haplotype_pairs[(size_t) i].second << ']';
            if (i + 1 != _n_pairs) {
                out << ", ";
            }
        }
        out << "\n\t]\n";
        out << "}\n";
    }

    void output_recent_summary() {
        ostream& out = *_recent_summary_output;

        // The first five columns are the historical schema, byte-identical, so
        // every existing consumer keeps working; mean_p_tmrca_lt_threshold is
        // the first threshold's mean probability. Per-threshold columns follow.
        out << "position_0based\tposition_1based\tn_pairs\tmean_p_tmrca_lt_threshold"
            << "\tmean_tmrca_generations";
        for (size_t k = 0; k < _thresholds.size(); k++) {
            const string suffix = format_threshold(_thresholds[k].years);
            out << "\tn_recent_" << suffix
                << "\tfrac_recent_" << suffix
                << "\tmean_p_lt_" << suffix;
        }
        out << "\n";

        out << std::setprecision(10);
        for (position_t i = 0; i < _seq_length; ++i) {
            const int64_t n_valid = _n_summary_pairs[(size_t) i];
            const double denom = (double) n_valid;
            out << _output_positions[(size_t) i] << '\t' << (_output_positions[(size_t) i] + 1) << '\t'
                << n_valid << '\t';
            if (n_valid > 0 && _accumulate_probability) {
                out << (_sum_probability[(size_t) i] / denom);
            } else {
                out << NAN;
            }
            out << '\t' << (n_valid > 0 ? (_sum_tmrca_generations[(size_t) i] / denom) : NAN);

            for (size_t k = 0; k < _thresholds.size(); k++) {
                const int64_t n_recent = _n_recent_pairs[(size_t) (k * _seq_length + i)];
                out << '\t' << n_recent
                    << '\t' << (n_valid > 0 ? ((double) n_recent / denom) : NAN)
                    << '\t';
                if (n_valid > 0 && _accumulate_probability) {
                    out << (_sum_probability[(size_t) (k * _seq_length + i)] / denom);
                } else {
                    out << NAN;
                }
            }
            out << '\n';
        }
    }

    // This just dumps the memory, so it later needs to be loaded in a particular way
    void output_raw_chunk(PairWorkspace& ws, bool last_chunk) {
        int finished;
        ZSTD_inBuffer input;
        ZSTD_EndDirective const mode = last_chunk ? ZSTD_e_end : ZSTD_e_continue;
        ZSTD_outBuffer output;

        input = {reinterpret_cast<void *>(ws.posteriors_alpha),
                 (size_t) (_seq_length * _n_pairs_in_chunk) * sizeof(float), 0};
        do {
            output = {_raw_compressed_buffer, _raw_compressed_buffer_size, 0};
            size_t const remaining = ZSTD_compressStream2(_raw_cctx, &output, &input, ZSTD_e_continue);
            _output_file_raw->write(
                reinterpret_cast<const char*>(_raw_compressed_buffer),
                (std::streamsize) output.pos
            );
            finished = last_chunk ? (remaining == 0) : (input.pos == input.size);
        } while (!finished);

        input = {reinterpret_cast<void *>(ws.posteriors_beta),
                 (size_t) (_seq_length * _n_pairs_in_chunk) * sizeof(float), 0};
        do {
            output = {_raw_compressed_buffer, _raw_compressed_buffer_size, 0};
            size_t const remaining = ZSTD_compressStream2(_raw_cctx, &output, &input, mode);
            _output_file_raw->write(
                reinterpret_cast<const char*>(_raw_compressed_buffer),
                (std::streamsize) output.pos
            );
            finished = last_chunk ? (remaining == 0) : (input.pos == input.size);
        } while (!finished);
    }

    void write_bitmatrix_meta(const string& filename) const {
        ofstream out(filename + ".meta", ios_base::out);
        out << "{\n";
        out << "\t\"format\": \"gamma_smc_bitmatrix_v1\",\n";
        out << boost::format("\t\"scaled_mutation_rate\": %.10f,\n") % _scaled_mutation_rate;
        out << boost::format("\t\"scaled_recombination_rate\": %.10f,\n") % _scaled_recombination_rate;
        out << boost::format("\t\"two_ne_generations\": %.10f,\n") % _two_ne_generations;
        out << "\t\"recent_call\": \"" << recent_call_name(_recent_call) << "\",\n";
        out << boost::format("\t\"n_pairs\": %ld,\n") % _n_pairs;
        out << boost::format("\t\"n_positions\": %ld,\n") % _seq_length;
        out << boost::format("\t\"pair_block\": %ld,\n") % _pair_block;
        out << boost::format("\t\"pairs_per_byte\": %ld,\n") % _n_pairs_in_chunk;
        out << boost::format("\t\"n_blocks\": %ld,\n") % _n_blocks;
        out << boost::format("\t\"zstd_compression_level\": %d,\n") % _zstd_compression_level;

        out << "\t\"thresholds_years\": [";
        for (size_t k = 0; k < _thresholds.size(); k++) {
            out << (k ? ", " : "") << format_threshold(_thresholds[k].years);
        }
        out << "],\n";
        out << "\t\"thresholds_generations\": [";
        for (size_t k = 0; k < _thresholds.size(); k++) {
            out << (k ? ", " : "") << _thresholds[k].generations;
        }
        out << "],\n";
        out << "\t\"thresholds_scaled\": [";
        for (size_t k = 0; k < _thresholds.size(); k++) {
            out << (k ? ", " : "") << _thresholds[k].scaled;
        }
        out << "],\n";

        // Positions on a plain stride are implied rather than listed; a
        // genome-wide list would otherwise dominate the metadata file.
        const bool strided = (!_output_at_hets) && (_posterior_every > 0);
        out << "\t\"positions_are_strided\": " << (strided ? "true" : "false") << ",\n";
        out << boost::format("\t\"stride\": %d,\n") % _posterior_every;
        if (strided) {
            out << "\t\"first_position\": " << (_seq_length > 0 ? _output_positions[0] : 0) << ",\n";
        } else {
            out << "\t\"output_positions\": [";
            for (position_t i = 0; i < _seq_length; i++) {
                out << (i ? ", " : "") << _output_positions[(size_t) i];
            }
            out << "],\n";
        }

        out << "\t\"pairs\": [";
        for (long i = 0; i < _n_pairs; i++) {
            out << (i ? ", " : "") << '[' << _haplotype_pairs[(size_t) i].first << ", "
                << _haplotype_pairs[(size_t) i].second << ']';
        }
        out << "],\n";

        out << "\t\"sample_names\": {";
        for (size_t i = 0; i < _data_processor._sample_names.size(); i++) {
            const string& sample_name = _data_processor._sample_names[i];
            out << (i ? ", " : "") << '"' << (i * 2) << "\": \"" << sample_name << ".0\", \""
                << (i * 2 + 1) << "\": \"" << sample_name << ".1\"";
        }
        out << "},\n";

        const vector<BitMatrixFrame>& frames = _bitmatrix_writer->_frames;
        out << "\t\"block_first_pair\": [";
        for (size_t i = 0; i < frames.size(); i++) {
            out << (i ? ", " : "") << frames[i].first_pair;
        }
        out << "],\n";
        out << "\t\"block_n_pairs\": [";
        for (size_t i = 0; i < frames.size(); i++) {
            out << (i ? ", " : "") << frames[i].n_pairs;
        }
        out << "],\n";
        out << "\t\"frame_offsets\": [";
        for (size_t i = 0; i < frames.size(); i++) {
            out << (i ? ", " : "") << frames[i].offset;
        }
        out << "],\n";
        out << "\t\"frame_sizes\": [";
        for (size_t i = 0; i < frames.size(); i++) {
            out << (i ? ", " : "") << frames[i].compressed_size;
        }
        out << "],\n";
        out << "\t\"frame_uncompressed_sizes\": [";
        for (size_t i = 0; i < frames.size(); i++) {
            out << (i ? ", " : "") << frames[i].uncompressed_size;
        }
        out << "]\n";
        out << "}\n";
        out.close();
    }

    // ----------------------------------------------------------------------
    // Driver
    // ----------------------------------------------------------------------

    void process_block(PairWorkspace& ws, long block_index) {
        const long block_first_pair = block_index * _pair_block;
        const long block_pairs = min(_pair_block, _n_pairs - block_first_pair);
        const long chunks_in_block = (block_pairs + _n_pairs_in_chunk - 1) / _n_pairs_in_chunk;

        if (_bitmatrix_output != NULL) {
            std::fill(ws.bits.begin(), ws.bits.begin() + bitmatrix_payload_bytes(chunks_in_block), 0);
        }

        for (long chunk = 0; chunk < chunks_in_block; chunk++) {
            const long chunk_first_pair = block_first_pair + chunk * _n_pairs_in_chunk;

            auto t1 = std::chrono::high_resolution_clock::now();
            if (!_data_processor.is_global_mask()) {
                prepare_pairwise_n_called(ws, chunk_first_pair);
            }
            prepare_pairwise_emissions(ws, chunk_first_pair);
            auto t2 = std::chrono::high_resolution_clock::now();
            ws.timer_emissions += std::chrono::duration<double>(t2 - t1).count();

            if (_only_backward) {
                // The backward pass accumulates onto the buffer, so it has to
                // start from 1 for every chunk, not once per allocation.
                const size_t elements = (size_t) (_seq_length * _n_pairs_in_chunk);
                std::fill(ws.posteriors_alpha, ws.posteriors_alpha + elements, 1.0f);
                std::fill(ws.posteriors_beta, ws.posteriors_beta + elements, 1.0f);
            }

            if (!_only_backward) {
                t1 = std::chrono::high_resolution_clock::now();
                forward_vectorized(ws);
                t2 = std::chrono::high_resolution_clock::now();
                ws.timer_forward += std::chrono::duration<double>(t2 - t1).count();
            }

            if (!_only_forward) {
                t1 = std::chrono::high_resolution_clock::now();
                backward_vectorized(ws);
                t2 = std::chrono::high_resolution_clock::now();
                ws.timer_backward += std::chrono::duration<double>(t2 - t1).count();
            }

            if (needs_statistics()) {
                t1 = std::chrono::high_resolution_clock::now();
                accumulate_and_pack(ws, chunk_first_pair, chunk, chunks_in_block);
                t2 = std::chrono::high_resolution_clock::now();
                ws.timer_output += std::chrono::duration<double>(t2 - t1).count();
            }
        }

        if (_bitmatrix_output != NULL) {
            const auto t1 = std::chrono::high_resolution_clock::now();
            const size_t payload = bitmatrix_payload_bytes(chunks_in_block);
            const size_t compressed_size = ZSTD_compress2(
                ws.cctx,
                ws.compressed.data(), ws.compressed.size(),
                ws.bits.data(), payload
            );
            if (ZSTD_isError(compressed_size)) {
                cerr << "Error: zstd failed on a bit-matrix block: "
                     << ZSTD_getErrorName(compressed_size) << endl;
                std::abort();
            }
            _bitmatrix_writer->write_frame(
                block_index, block_first_pair, block_pairs,
                ws.compressed.data(), compressed_size, payload
            );
            const auto t2 = std::chrono::high_resolution_clock::now();
            ws.timer_output += std::chrono::duration<double>(t2 - t1).count();
        }
    }

    void calculate_posteriors() {
        const auto wall_start = std::chrono::high_resolution_clock::now();

        if (_output_file_raw != NULL) {
            output_raw_header();
            _raw_cctx = ZSTD_createCCtx();
            _raw_compressed_buffer_size = ZSTD_CStreamOutSize();
            _raw_compressed_buffer = malloc(_raw_compressed_buffer_size);
            ZSTD_CCtx_setParameter(_raw_cctx, ZSTD_c_compressionLevel, _zstd_compression_level);
            ZSTD_CCtx_setParameter(_raw_cctx, ZSTD_c_checksumFlag, 1);
        }

        BlockProgressBar bar{
            option::BarWidth{70},
            option::FontStyles{
                std::vector<FontStyle>{FontStyle::bold}},
            option::ShowElapsedTime{true},
            option::ShowRemainingTime{true},
            option::MaxProgress{_n_pairs}
        };
        // Touch the bar once before any thread exists: termcolor's first call
        // writes the stream's iword storage, which is not safe to race on.
        bar.set_progress(0);

        std::atomic<long> pairs_done(0);

        // -ffast-math sets FTZ/DAZ in MXCSR for the main thread only, and MXCSR
        // is per-thread on x86-64. Copy it so a pair's result cannot depend on
        // which worker happened to pick up its block.
        const unsigned int main_mxcsr = _mm_getcsr();

        vector<PairWorkspace> workspaces((size_t) _n_threads);

#ifdef _OPENMP
#pragma omp parallel num_threads(_n_threads)
#endif
        {
#ifdef _OPENMP
            const int thread_id = omp_get_thread_num();
#else
            const int thread_id = 0;
#endif
            _mm_setcsr(main_mxcsr);
            PairWorkspace& ws = workspaces[(size_t) thread_id];
            allocate_workspace(ws);

            if (_output_file_raw != NULL) {
                // Raw posteriors must reach the stream in block order.
#ifdef _OPENMP
#pragma omp for schedule(static, 1) ordered
#endif
                for (long block = 0; block < _n_blocks; block++) {
                    process_block(ws, block);
#ifdef _OPENMP
#pragma omp ordered
#endif
                    {
                        const auto t1 = std::chrono::high_resolution_clock::now();
                        output_raw_chunk(ws, block == (_n_blocks - 1));
                        const auto t2 = std::chrono::high_resolution_clock::now();
                        ws.timer_output += std::chrono::duration<double>(t2 - t1).count();
                    }
                    report_progress(bar, pairs_done, block);
                }
            } else {
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
                for (long block = 0; block < _n_blocks; block++) {
                    process_block(ws, block);
                    report_progress(bar, pairs_done, block);
                }
            }

#ifdef _OPENMP
#pragma omp critical(gamma_smc_merge)
#endif
            {
                if (needs_statistics()) {
                    merge_workspace(ws);
                } else {
                    _timer_emissions += ws.timer_emissions;
                    _timer_forward += ws.timer_forward;
                    _timer_backward += ws.timer_backward;
                    _timer_output += ws.timer_output;
                }
            }
        }

        bar.set_option(option::PostfixText{
            std::to_string(_n_pairs) + "/" + std::to_string(_n_pairs)
        });
        bar.set_progress(_n_pairs);
        bar.mark_as_completed();
        indicators::show_console_cursor(true);

        if (_raw_cctx != NULL) {
            ZSTD_freeCCtx(_raw_cctx);
            _raw_cctx = NULL;
            free(_raw_compressed_buffer);
            _raw_compressed_buffer = NULL;
        }

        if (_recent_summary_output != NULL) {
            output_recent_summary();
        }

        if (_bitmatrix_output != NULL) {
            _bitmatrix_writer->sort_frames();
        }

        _timer_wall = std::chrono::duration<double>(
            std::chrono::high_resolution_clock::now() - wall_start).count();
    }

    void report_progress(BlockProgressBar& bar, std::atomic<long>& pairs_done, long block) {
        const long block_first_pair = block * _pair_block;
        const long block_pairs = min(_pair_block, _n_pairs - block_first_pair);
        const long done = pairs_done.fetch_add(block_pairs) + block_pairs;
#ifdef _OPENMP
#pragma omp critical(gamma_smc_progress)
#endif
        {
            bar.set_option(option::PostfixText{
                std::to_string(done) + "/" + std::to_string(_n_pairs)
            });
            bar.set_progress(done);
        }
    }

    ZSTD_CCtx* _raw_cctx = NULL;
    void* _raw_compressed_buffer = NULL;
    size_t _raw_compressed_buffer_size = 0;

    virtual ~CachedPairwiseGammaSMC() {
        free(_global_n_called);
        if (_raw_cctx != NULL) {
            ZSTD_freeCCtx(_raw_cctx);
        }
        free(_raw_compressed_buffer);
    }
};
