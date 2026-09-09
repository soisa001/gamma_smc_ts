// Replay the decoder's native hard-call kernel on saved gamma posteriors.
// No HMM, flow-field, demographic, or pair-selection changes are made here.
#include "recent_stats.h"
#include <sstream>
#include <iomanip>

int main(int argc, char** argv) {
    try {
        if (argc != 9) throw runtime_error("Usage: summarize_recent_rules RAW N_PAIRS POSITIONS THETA MU GEN YEARS_CSV OUTDIR");
        const long pairs = stol(argv[2]);
        vector<long> positions;
        ifstream position_file(argv[3]);
        long position;
        while (position_file >> position) positions.push_back(position);
        if (pairs <= 0 || positions.empty()) throw runtime_error("Empty pair/position manifest");
        const float theta = stof(argv[4]);
        const double two_ne = theta / (2.0 * stod(argv[5]));
        const double generation = stod(argv[6]);
        if (!(two_ne > 0 && generation > 0)) throw runtime_error("Invalid time units");
        vector<long> years;
        string value;
        stringstream input_years(argv[7]);
        while (getline(input_years, value, ',')) years.push_back(stol(value));
        if (years.empty()) throw runtime_error("No thresholds");
        vector<float> thresholds;
        for (long year : years) thresholds.push_back((float) ((year / generation) / two_ne));

        ifstream raw(argv[1], ios::binary | ios::ate);
        if (!raw) throw runtime_error("Cannot open raw posterior file");
        const size_t compressed_size = (size_t) raw.tellg();
        vector<char> compressed(compressed_size);
        raw.seekg(0);
        raw.read(compressed.data(), compressed.size());
        if (!raw) throw runtime_error("Truncated compressed input read");
        const size_t cells = (size_t) ((pairs + 7) / 8) * positions.size() * 8 * 2;
        vector<float> posterior(cells);
        const size_t written = ZSTD_decompress(posterior.data(), cells * sizeof(float),
                                                compressed.data(), compressed.size());
        if (ZSTD_isError(written) || written != cells * sizeof(float))
            throw runtime_error("Corrupt posterior stream or unexpected decompressed length");
        const vector<string> rules = {"mean", "median", "prob80", "prob90"};
        vector<RecentThreshold> tables(4);
        tables[0].build(RECENT_CALL_MEAN, 0.5);
        tables[1].build(RECENT_CALL_MEDIAN, 0.5);
        tables[2].build(RECENT_CALL_PROB, 0.8);
        tables[3].build(RECENT_CALL_PROB, 0.9);
        const size_t npos = positions.size();
        vector<vector<long>> counts(4, vector<long>(npos * years.size(), 0));
        const size_t block = npos * 8;
        for (long first = 0; first < pairs; first += 8) {
            const int real = min(8L, pairs - first);
            const unsigned lane_mask = (1u << real) - 1;
            const float* alpha = posterior.data() + (size_t) (first / 8) * block * 2;
            const float* beta = alpha + block;
            for (size_t p = 0; p < npos; p++) {
                for (int lane = 0; lane < real; lane++) {
                    for (float x : {alpha[p*8+lane], beta[p*8+lane]}) {
                        // Integer finite check is preserved under -ffast-math.
                        if (!(x > 0) || (float_bits(x) & 0x7fffffffu) >= 0x7f800000u)
                            throw runtime_error("Invalid posterior parameter");
                    }
                }
                const __m256 a = _mm256_loadu_ps(alpha + p*8);
                const __m256 b = _mm256_loadu_ps(beta + p*8);
                for (size_t t = 0; t < years.size(); t++) {
                    const __m256 x = _mm256_mul_ps(b, _mm256_set1_ps(thresholds[t]));
                    for (size_t rule = 0; rule < rules.size(); rule++) {
                        const unsigned calls = (unsigned) _mm256_movemask_ps(recent_call_mask(tables[rule], a, x));
                        counts[rule][p*years.size()+t] += __builtin_popcount(calls & lane_mask);
                    }
                }
            }
        }
        filesystem::create_directories(argv[8]);
        for (size_t rule = 0; rule < rules.size(); rule++) {
            const auto destination = filesystem::path(argv[8]) / (rules[rule] + ".tsv");
            ofstream out(destination);
            out << "position_0based";
            for (long year : years) out << "\tfrac_recent_" << year;
            out << '\n' << setprecision(10);
            for (size_t p = 0; p < npos; p++) {
                out << positions[p];
                for (size_t t = 0; t < years.size(); t++)
                    out << '\t' << (double) counts[rule][p*years.size()+t] / pairs;
                out << '\n';
            }
            out.close();
            if (!out) throw runtime_error("Failed to write summary");
        }
        cout << "Validated " << pairs << " pairs x " << npos << " positions; wrote four native hard-call summaries\n";
        return 0;
    } catch (const exception& error) {
        cerr << error.what() << '\n';
        return 1;
    }
}
