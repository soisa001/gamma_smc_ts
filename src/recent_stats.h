#pragma once

// Recent-coalescence statistics from a Gamma(alpha, beta) TMRCA posterior.
//
// For each (pair, position) and each time threshold t (in coalescent units)
// we need two things:
//
//   * a hard call  -- did this pair coalesce within t?      -> one bit
//   * P(T < t)     -- so the across-pair mean keeps its meaning
//
// with both evaluated ~2.5e10 times per chromosome at 100k pairs and stride
// 1000. boost::math::gamma_p per cell costs a few hundred nanoseconds, which on
// its own would take longer than the rest of the run, so both go through lookup
// tables built once at startup from the exact Boost functions.
//
// THE HARD CALL IS A ONE-DIMENSIONAL PROBLEM
// ------------------------------------------
// With x = beta * t,
//
//     P(T < t) = P(alpha, x)          (regularized lower incomplete gamma)
//
// and P is strictly increasing in x, so for a fixed probability p
//
//     P(alpha, x) >= p    <=>    x >= Q_p(alpha)
//
// where Q_p(alpha) = gamma_p_inv(alpha, p) is the p-quantile of Gamma(alpha, 1).
// That is a function of alpha alone: one table, one compare, no 2-D anything.
// The posterior-mean rule falls out of the same machinery, since
// alpha/beta < t <=> x > alpha, i.e. Q(alpha) = alpha.
//
// NO LOGARITHMS ANYWHERE
// ----------------------
// Both tables are indexed by
//
//     alpha_coord(alpha) = exponent(alpha) + mantissa_fraction(alpha)
//
// read straight off the float32 bit pattern. It is the usual piecewise-linear
// approximation of log2: strictly increasing, exact at powers of two, and
// obtainable with an integer shift and a mask. Accuracy does not matter,
// because the tables are *built* at the alpha values this coordinate decodes
// to, so build and lookup agree by construction and interpolation is exact at
// grid nodes.

#include "common.h"

#include <boost/math/special_functions/gamma.hpp>

// ---------------------------------------------------------------------------
// alpha <-> table coordinate
// ---------------------------------------------------------------------------

// Covered range of alpha: [2^-16, 2^16). Posteriors outside are clamped.
static const int recent_alpha_min_exponent = -16;
static const int recent_alpha_max_exponent = 16;

inline uint32_t float_bits(float value) {
    uint32_t bits;
    memcpy(&bits, &value, sizeof(bits));
    return bits;
}

inline float bits_to_float(uint32_t bits) {
    float value;
    memcpy(&value, &bits, sizeof(value));
    return value;
}

// exponent + mantissa fraction, i.e. the piecewise-linear log2 described above.
inline float alpha_coord(float alpha) {
    const uint32_t bits = float_bits(alpha);
    return (float) ((int32_t) (bits >> 23) - 127)
         + (float) (bits & 0x7FFFFFu) * (1.0f / 8388608.0f);
}

// Inverse of alpha_coord, used to build the tables.
inline float alpha_from_coord(double coord) {
    double exponent_part = std::floor(coord);
    double fraction = coord - exponent_part;
    if (exponent_part < recent_alpha_min_exponent) {
        exponent_part = recent_alpha_min_exponent;
        fraction = 0.0;
    }
    if (exponent_part > 127.0) {
        exponent_part = 127.0;
        fraction = 0.0;
    }
    uint32_t mantissa = (uint32_t) (fraction * 8388608.0);
    if (mantissa > 0x7FFFFFu) {
        mantissa = 0x7FFFFFu;
    }
    return bits_to_float((((uint32_t) ((int32_t) exponent_part + 127)) << 23) | mantissa);
}

// ---------------------------------------------------------------------------
// Decision rule
// ---------------------------------------------------------------------------

enum recent_call_t : int8_t {
    RECENT_CALL_MEDIAN = 0,   // P(T < t) >= 0.5, i.e. the posterior median is below t
    RECENT_CALL_PROB = 1,     // P(T < t) >= call_probability
    RECENT_CALL_MEAN = 2,     // posterior mean alpha/beta < t
};

inline const char* recent_call_name(recent_call_t rule) {
    switch (rule) {
        case RECENT_CALL_MEDIAN: return "median";
        case RECENT_CALL_PROB: return "prob";
        case RECENT_CALL_MEAN: return "mean";
    }
    return "unknown";
}

// ---------------------------------------------------------------------------
// Quantile table: the hard call
// ---------------------------------------------------------------------------

// 2^8 = 256 table points per octave of alpha.
static const int recent_quantile_mantissa_bits = 8;
static const int recent_quantile_shift = 23 - recent_quantile_mantissa_bits;   // 15
static const uint32_t recent_quantile_frac_mask = (1u << recent_quantile_shift) - 1u;
static const float recent_quantile_frac_scale = 1.0f / (float) (1u << recent_quantile_shift);
static const int recent_quantile_n_octaves = recent_alpha_max_exponent - recent_alpha_min_exponent;
static const int recent_quantile_table_size =
    recent_quantile_n_octaves * (1 << recent_quantile_mantissa_bits) + 2;

inline uint32_t recent_quantile_base_index() {
    return ((uint32_t) (recent_alpha_min_exponent + 127)) << recent_quantile_mantissa_bits;
}

// The alpha value that quantile-table slot `index` stands for.
inline float recent_quantile_alpha_at(int index) {
    return bits_to_float((recent_quantile_base_index() + (uint32_t) index) << recent_quantile_shift);
}

class RecentThreshold {
  public:
    double years = 0.0;
    double generations = 0.0;
    double scaled = 0.0;                 // threshold in coalescent (2Ne generation) units

    // Bit is set when beta * scaled >= quantile[idx(alpha)].
    vector<float> quantile;

    // Kept as a plain float rather than a broadcast __m256: an over-aligned
    // member would make vector<RecentThreshold> depend on aligned-new support.
    float scaled_float = 0.0f;
    float alpha_min = 0.0f;
    float alpha_max = 0.0f;

    void build(recent_call_t rule, double call_probability) {
        quantile.assign((size_t) recent_quantile_table_size, 0.0f);
        const bool use_mean_rule = (rule == RECENT_CALL_MEAN);
        const double p = (rule == RECENT_CALL_MEDIAN) ? 0.5 : call_probability;

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int index = 0; index < recent_quantile_table_size; index++) {
            const double alpha = (double) recent_quantile_alpha_at(index);
            quantile[(size_t) index] = use_mean_rule
                ? (float) alpha
                : (float) boost::math::gamma_p_inv(alpha, p);
        }

        alpha_min = recent_quantile_alpha_at(0);
        alpha_max = recent_quantile_alpha_at(recent_quantile_table_size - 2);
        scaled_float = (float) scaled;
    }

    // Scalar reference used by the self-check and by the low-volume paths.
    bool calls_recent(float alpha, float beta) const {
        if (!(alpha > 0.0f) || !(beta > 0.0f)) {
            return false;
        }
        const float clamped = std::min(std::max(alpha, alpha_min), alpha_max);
        const uint32_t bits = float_bits(clamped);
        const int index = (int) ((bits >> recent_quantile_shift) - recent_quantile_base_index());
        const float frac = (float) (bits & recent_quantile_frac_mask) * recent_quantile_frac_scale;
        const float q0 = quantile[(size_t) index];
        const float q1 = quantile[(size_t) index + 1];
        return (beta * (float) scaled) >= (q0 + frac * (q1 - q0));
    }
};

// ---------------------------------------------------------------------------
// P(alpha, x) table: the continuous statistic
// ---------------------------------------------------------------------------
//
// Second axis is the standardized deviate
//
//     z = (x / alpha - 1) * sqrt(alpha)
//
// under which the transition of P sits at z ~ 0 with O(1) width for every
// alpha, instead of collapsing to a spike of width 1/sqrt(alpha) as it would in
// x or log x. The lower end is safe because x >= 0 forces z >= -sqrt(alpha), so
// clamping to P = 0 below z_min is exact whenever alpha <= z_min^2; the upper
// end runs out to z = 20 so that the heavy small-alpha tail has converged
// before the clamp to 1.

static const int recent_prob_alpha_steps = 512;    // over alpha_coord in [-16, 16]
static const int recent_prob_z_steps = 2048;       // over z in [-8, 20]
static const float recent_prob_z_min = -8.0f;
static const float recent_prob_z_max = 20.0f;

class RecentProbabilityTable {
  public:
    vector<float> _table;              // recent_prob_alpha_steps x recent_prob_z_steps
    float _alpha_scale = 0.0f;         // grid = alpha_coord * scale + offset
    float _alpha_offset = 0.0f;
    float _z_scale = 0.0f;
    float _z_offset = 0.0f;
    double _max_abs_error = 0.0;

    void build() {
        _table.assign((size_t) recent_prob_alpha_steps * recent_prob_z_steps, 0.0f);

        const double alpha_step =
            (double) (recent_alpha_max_exponent - recent_alpha_min_exponent)
            / (recent_prob_alpha_steps - 1);
        const double z_step =
            (double) (recent_prob_z_max - recent_prob_z_min) / (recent_prob_z_steps - 1);

        _alpha_scale = (float) (1.0 / alpha_step);
        _alpha_offset = (float) (-(double) recent_alpha_min_exponent / alpha_step);
        _z_scale = (float) (1.0 / z_step);
        _z_offset = (float) (-(double) recent_prob_z_min / z_step);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int a = 0; a < recent_prob_alpha_steps; a++) {
            // Build at exactly the alpha the lookup coordinate decodes to.
            const double alpha =
                (double) alpha_from_coord(recent_alpha_min_exponent + a * alpha_step);
            const double sqrt_alpha = std::sqrt(alpha);
            for (int k = 0; k < recent_prob_z_steps; k++) {
                const double z = recent_prob_z_min + k * z_step;
                const double x = alpha + z * sqrt_alpha;
                _table[(size_t) a * recent_prob_z_steps + k] =
                    (x <= 0.0) ? 0.0f : (float) boost::math::gamma_p(alpha, x);
            }
        }
    }

    // Scalar reference; mirrors the SIMD path so the self-check measures what
    // the hot loop will actually compute.
    float lookup(float alpha, float x) const {
        if (!(alpha > 0.0f) || !(x > 0.0f)) {
            return 0.0f;
        }
        const float sqrt_alpha = std::sqrt(alpha);
        const float z = (x / alpha - 1.0f) * sqrt_alpha;
        if (z <= recent_prob_z_min) {
            return 0.0f;
        }
        if (z >= recent_prob_z_max) {
            return 1.0f;
        }

        float ga = alpha_coord(alpha) * _alpha_scale + _alpha_offset;
        ga = std::min(std::max(ga, 0.0f), (float) (recent_prob_alpha_steps - 1) - 1.0f / 512.0f);
        const float gz = z * _z_scale + _z_offset;

        const int a0 = (int) ga;
        const int k0 = (int) gz;
        const float wa = ga - a0;
        const float wk = gz - k0;

        const size_t base = (size_t) a0 * recent_prob_z_steps + k0;
        const float v00 = _table[base];
        const float v01 = _table[base + 1];
        const float v10 = _table[base + recent_prob_z_steps];
        const float v11 = _table[base + recent_prob_z_steps + 1];

        const float lower = v00 + wk * (v01 - v00);
        const float upper = v10 + wk * (v11 - v10);
        return lower + wa * (upper - lower);
    }

    // Compare against Boost over the parameter region the decoder actually
    // visits, deliberately landing between nodes where bilinear interpolation is
    // at its worst, so a bad grid shows up as a printed number rather than as a
    // silently biased statistic.
    void self_check() {
        double worst = 0.0;
        const int alpha_samples = 512;
        const int z_samples = 512;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(max : worst)
#endif
        for (int i = 0; i < alpha_samples; i++) {
            const double coord = -4.0 + (i + 0.5) * (18.0 / alpha_samples);
            const double alpha = (double) alpha_from_coord(coord);
            const double sqrt_alpha = std::sqrt(alpha);
            for (int k = 0; k < z_samples; k++) {
                const double z = -7.5 + (k + 0.5) * (27.0 / z_samples);
                const double x = alpha + z * sqrt_alpha;
                if (x <= 0.0) {
                    continue;
                }
                const double exact = boost::math::gamma_p(alpha, x);
                const double approx = (double) lookup((float) alpha, (float) x);
                worst = std::max(worst, std::fabs(exact - approx));
            }
        }
        _max_abs_error = worst;
    }
};

// ---------------------------------------------------------------------------
// SIMD kernels
// ---------------------------------------------------------------------------

// A float32 lane is a usable posterior parameter iff its bit pattern, read as a
// signed int32, lies in [1, 0x7F7FFFFF]: that excludes zero, every negative
// value (sign bit set => negative as int32), and inf/nan (>= 0x7F800000).
// Doing it on the integer side survives -ffast-math, which is free to delete
// isfinite() checks.
inline __m256i positive_finite_mask(__m256 value) {
    const __m256i bits = _mm256_castps_si256(value);
    return _mm256_and_si256(
        _mm256_cmpgt_epi32(bits, _mm256_setzero_si256()),
        _mm256_cmpgt_epi32(_mm256_set1_epi32(0x7F800000), bits)
    );
}

// Hard call for eight pairs at once: does x = beta*t reach the p-quantile of
// Gamma(alpha, 1)?
inline __m256 recent_call_mask(
    const RecentThreshold& threshold,
    __m256 alpha,
    __m256 x
) {
    const __m256 clamped = _mm256_min_ps(
        _mm256_max_ps(alpha, _mm256_set1_ps(threshold.alpha_min)),
        _mm256_set1_ps(threshold.alpha_max)
    );
    const __m256i bits = _mm256_castps_si256(clamped);
    const __m256i index = _mm256_sub_epi32(
        _mm256_srli_epi32(bits, recent_quantile_shift),
        _mm256_set1_epi32((int32_t) recent_quantile_base_index())
    );
    const __m256 frac = _mm256_mul_ps(
        _mm256_cvtepi32_ps(_mm256_and_si256(bits, _mm256_set1_epi32((int32_t) recent_quantile_frac_mask))),
        _mm256_set1_ps(recent_quantile_frac_scale)
    );

    const float* table = threshold.quantile.data();
    const __m256 q0 = _mm256_i32gather_ps(table, index, 4);
    const __m256 q1 = _mm256_i32gather_ps(table + 1, index, 4);
    const __m256 q = _mm256_fmadd_ps(frac, _mm256_sub_ps(q1, q0), q0);

    return _mm256_cmp_ps(x, q, _CMP_GE_OQ);
}

// Bilinear P(alpha, x) for eight pairs at once.
inline __m256 recent_probability(
    const RecentProbabilityTable& table,
    __m256 alpha,
    __m256 x
) {
    const __m256 one = _mm256_set1_ps(1.0f);
    const __m256 sqrt_alpha = _mm256_sqrt_ps(alpha);
    const __m256 z = _mm256_mul_ps(
        _mm256_sub_ps(_mm256_div_ps(x, alpha), one),
        sqrt_alpha
    );

    const __m256 below = _mm256_cmp_ps(z, _mm256_set1_ps(recent_prob_z_min), _CMP_LE_OQ);
    const __m256 above = _mm256_cmp_ps(z, _mm256_set1_ps(recent_prob_z_max), _CMP_GE_OQ);

    // alpha_coord without a log: exponent + mantissa fraction.
    const __m256i bits = _mm256_castps_si256(alpha);
    const __m256 coord = _mm256_add_ps(
        _mm256_cvtepi32_ps(_mm256_sub_epi32(
            _mm256_srli_epi32(bits, 23), _mm256_set1_epi32(127)
        )),
        _mm256_mul_ps(
            _mm256_cvtepi32_ps(_mm256_and_si256(bits, _mm256_set1_epi32(0x7FFFFF))),
            _mm256_set1_ps(1.0f / 8388608.0f)
        )
    );

    __m256 ga = _mm256_fmadd_ps(coord, _mm256_set1_ps(table._alpha_scale), _mm256_set1_ps(table._alpha_offset));
    ga = _mm256_min_ps(
        _mm256_max_ps(ga, _mm256_setzero_ps()),
        _mm256_set1_ps((float) (recent_prob_alpha_steps - 1) - 1.0f / 512.0f)
    );

    __m256 gz = _mm256_fmadd_ps(z, _mm256_set1_ps(table._z_scale), _mm256_set1_ps(table._z_offset));
    gz = _mm256_min_ps(
        _mm256_max_ps(gz, _mm256_setzero_ps()),
        _mm256_set1_ps((float) (recent_prob_z_steps - 2))
    );

    const __m256 ga_floor = _mm256_floor_ps(ga);
    const __m256 gz_floor = _mm256_floor_ps(gz);
    const __m256 wa = _mm256_sub_ps(ga, ga_floor);
    const __m256 wk = _mm256_sub_ps(gz, gz_floor);

    const __m256i base = _mm256_add_epi32(
        _mm256_mullo_epi32(_mm256_cvttps_epi32(ga_floor), _mm256_set1_epi32(recent_prob_z_steps)),
        _mm256_cvttps_epi32(gz_floor)
    );

    const float* data = table._table.data();
    const __m256 v00 = _mm256_i32gather_ps(data, base, 4);
    const __m256 v01 = _mm256_i32gather_ps(data + 1, base, 4);
    const __m256 v10 = _mm256_i32gather_ps(data + recent_prob_z_steps, base, 4);
    const __m256 v11 = _mm256_i32gather_ps(data + recent_prob_z_steps + 1, base, 4);

    const __m256 lower = _mm256_fmadd_ps(wk, _mm256_sub_ps(v01, v00), v00);
    const __m256 upper = _mm256_fmadd_ps(wk, _mm256_sub_ps(v11, v10), v10);
    __m256 result = _mm256_fmadd_ps(wa, _mm256_sub_ps(upper, lower), lower);

    result = _mm256_blendv_ps(result, _mm256_setzero_ps(), below);
    result = _mm256_blendv_ps(result, one, above);
    return result;
}

inline float horizontal_sum(__m256 v) {
    __m128 low = _mm256_castps256_ps128(v);
    const __m128 high = _mm256_extractf128_ps(v, 1);
    low = _mm_add_ps(low, high);
    low = _mm_add_ps(low, _mm_movehl_ps(low, low));
    low = _mm_add_ss(low, _mm_shuffle_ps(low, low, 1));
    return _mm_cvtss_f32(low);
}
