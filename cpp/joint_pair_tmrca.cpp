// Vectorized exact MRCA ages on one tskit marginal tree. No approximations.
// The decoder is unchanged; this accelerates extraction of simulation truth.
#include <cstdint>
#include <limits>

extern "C" int joint_pair_tmrca(const int32_t* parent, const double* times,
                                const int32_t* pairs, int64_t n_pairs,
                                int64_t n_nodes, double* output) {
    for (int64_t i = 0; i < n_pairs; ++i) {
        int32_t a = pairs[2*i], b = pairs[2*i+1];
        if (a < 0 || b < 0 || a >= n_nodes || b >= n_nodes) return 1;
        while (a != b) {
            if (times[a] < times[b]) a = parent[a];
            else b = parent[b];
            if (a < 0 || b < 0) return 2;
        }
        output[i] = times[a];
    }
    return 0;
}
