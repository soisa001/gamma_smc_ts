#pragma once

// A durable record of which haplotype pairs a run actually decoded.
//
// --n_random_pairs is seeded, so in principle the draw is reproducible from
// --pairs_seed alone. In practice that is a thin guarantee: the indices it
// produces mean nothing without the haplotype ordering they index into, which
// depends on the input file, on --samples, and on which records htslib kept.
// Re-deriving a draw months later from a seed, against a panel that has since
// been re-exported or re-subset, silently decodes different pairs.
//
// So the draw is written out literally. The manifest is a TSV whose header
// lines all start with '#', which read_pairs_file already skips, and whose
// first two columns are the haplotype indices, which read_pairs_file already
// parses while ignoring the rest. A manifest is therefore directly usable as
// --pairs_file with no conversion:
//
//     gamma_smc ... --n_random_pairs 100000 --pairs_manifest chr1.pairs.tsv
//     gamma_smc ... --pairs_file chr1.pairs.tsv          # same pairs, chr2
//
// The header also carries a digest of the ordered sample names. Reusing a
// manifest against a panel whose digest differs is refused rather than
// silently producing pairs that point somewhere else.

#include "common.h"

#include <boost/algorithm/string.hpp>

#include <ctime>
#include <sstream>

const char* const PAIR_MANIFEST_MAGIC = "gamma_smc_pair_manifest_v1";

// FNV-1a. Not cryptographic -- this guards against using a manifest with the
// wrong panel by accident, not against a forged one.
inline uint64_t fnv1a(const void* data, size_t length, uint64_t hash = 1469598103934665603ULL) {
    const unsigned char* bytes = (const unsigned char*) data;
    for (size_t i = 0; i < length; i++) {
        hash ^= (uint64_t) bytes[i];
        hash *= 1099511628211ULL;
    }
    return hash;
}

// Digest of the haplotype universe: the ordered sample names, which is exactly
// what the pair indices are relative to.
inline uint64_t panel_digest(const vector<string>& sample_names) {
    uint64_t hash = fnv1a("", 0);
    for (const string& name : sample_names) {
        hash = fnv1a(name.data(), name.size(), hash);
        hash = fnv1a("\0", 1, hash);
    }
    return hash;
}

inline uint64_t pairs_digest(const vector<pair<int, int>>& pairs) {
    uint64_t hash = fnv1a("", 0);
    for (const pair<int, int>& p : pairs) {
        const int32_t values[2] = {(int32_t) p.first, (int32_t) p.second};
        hash = fnv1a(values, sizeof(values), hash);
    }
    return hash;
}

inline string utc_timestamp() {
    const std::time_t now = std::time(NULL);
    std::tm parts;
#ifdef _WIN32
    gmtime_s(&parts, &now);
#else
    gmtime_r(&now, &parts);
#endif
    char buffer[32];
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%dT%H:%M:%SZ", &parts);
    return string(buffer);
}

inline string hex64(uint64_t value) {
    std::ostringstream stream;
    stream << "0x" << std::hex << std::setw(16) << std::setfill('0') << value;
    return stream.str();
}

struct PairManifestHeader {
    bool present = false;
    long n_pairs = -1;
    int n_haplotypes = -1;
    uint64_t panel_digest = 0;
    bool has_panel_digest = false;
    string source;
};

// Haplotype h is haplotype (h & 1) of diploid sample (h >> 1).
inline string haplotype_label(int haplotype, const vector<string>& sample_names) {
    const size_t sample = (size_t) (haplotype >> 1);
    if (sample >= sample_names.size()) {
        return string("?");
    }
    return sample_names[sample] + ((haplotype & 1) ? ".1" : ".0");
}

inline void write_pair_manifest(
    const string& filename,
    const vector<pair<int, int>>& pairs,
    const vector<string>& sample_names,
    const string& input_filename,
    const string& mode,
    long requested,
    uint64_t seed,
    bool exclude_within,
    bool seed_is_meaningful
) {
    ofstream out(filename, ios_base::out);
    if (!out.good()) {
        cout << boost::format("Error: Cannot write --pairs_manifest %s\n") % filename;
        exit(-1);
    }

    const int n_haplotypes = (int) (sample_names.size() * 2);
    out << "# " << PAIR_MANIFEST_MAGIC << "\n";
    out << "# created_utc\t" << utc_timestamp() << "\n";
    out << "# input\t" << input_filename << "\n";
    out << "# mode\t" << mode << "\n";
    out << boost::format("# n_pairs\t%ld\n") % (long) pairs.size();
    out << boost::format("# n_haplotypes\t%d\n") % n_haplotypes;
    out << boost::format("# n_samples\t%ld\n") % (long) sample_names.size();
    if (seed_is_meaningful) {
        out << boost::format("# pairs_seed\t%llu\n") % (unsigned long long) seed;
        out << boost::format("# n_random_pairs\t%ld\n") % requested;
        out << "# exclude_within\t" << (exclude_within ? "true" : "false") << "\n";
        out << "# rng\tmt19937_64; draw both haplotypes uniformly, reject i==j, "
               "deduplicate, sort ascending\n";
    }
    out << "# panel_digest\t" << hex64(panel_digest(sample_names)) << "\n";
    out << "# pairs_digest\t" << hex64(pairs_digest(pairs)) << "\n";
    out << "# Reuse with --pairs_file to decode exactly these pairs again.\n";
    out << "# hap_i\thap_j\thaplotype_i\thaplotype_j\n";

    for (const pair<int, int>& p : pairs) {
        out << p.first << '\t' << p.second << '\t'
            << haplotype_label(p.first, sample_names) << '\t'
            << haplotype_label(p.second, sample_names) << '\n';
    }
    out.close();
}

// Pull our header back out of a pairs file, if it has one. Files hand-written
// by the user have no header; that is fine and simply skips the checks.
inline PairManifestHeader read_pair_manifest_header(const string& filename) {
    PairManifestHeader header;
    ifstream input(filename);
    if (!input.is_open()) {
        return header;
    }

    string line;
    bool first = true;
    while (getline(input, line)) {
        boost::algorithm::trim(line);
        if (line.empty()) {
            continue;
        }
        if (line[0] != '#') {
            break;      // header ends at the first pair row
        }
        string body = line.substr(1);
        boost::algorithm::trim(body);
        if (first) {
            if (body != PAIR_MANIFEST_MAGIC) {
                return header;      // some other file's comments
            }
            header.present = true;
            first = false;
            continue;
        }

        std::istringstream stream(body);
        string key, value;
        if (!(stream >> key >> value)) {
            continue;
        }
        try {
            if (key == "n_pairs") {
                header.n_pairs = std::stol(value);
            } else if (key == "n_haplotypes") {
                header.n_haplotypes = std::stoi(value);
            } else if (key == "panel_digest") {
                header.panel_digest = std::stoull(value, NULL, 16);
                header.has_panel_digest = true;
            } else if (key == "input") {
                header.source = value;
            }
        } catch (const std::exception&) {
            // A malformed header field is not worth aborting a run over; the
            // checks below simply do not fire for it.
        }
    }
    return header;
}

// Refuse a manifest drawn against a different panel: the indices would point at
// different haplotypes and the run would look perfectly healthy.
inline void verify_pair_manifest(
    const PairManifestHeader& header,
    const vector<string>& sample_names,
    const string& filename,
    bool allow_mismatch
) {
    if (!header.present) {
        return;
    }
    const int n_haplotypes = (int) (sample_names.size() * 2);
    const uint64_t digest = panel_digest(sample_names);

    const bool haplotypes_differ =
        (header.n_haplotypes >= 0) && (header.n_haplotypes != n_haplotypes);
    const bool digest_differs =
        header.has_panel_digest && (header.panel_digest != digest);

    if (!haplotypes_differ && !digest_differs) {
        return;
    }

    cout << boost::format(
        "%s: --pairs_file %s was drawn against a different panel.\n"
    ) % (allow_mismatch ? "Warning" : "Error") % filename;
    if (haplotypes_differ) {
        cout << boost::format("  haplotypes: manifest %d, this run %d\n")
                % header.n_haplotypes % n_haplotypes;
    }
    if (digest_differs) {
        cout << boost::format("  panel digest: manifest %s, this run %s\n")
                % hex64(header.panel_digest) % hex64(digest);
        if (!header.source.empty()) {
            cout << boost::format("  manifest was drawn from: %s\n") % header.source;
        }
    }
    cout << "  Haplotype indices are relative to the sample order, so these pairs\n"
            "  would decode different haplotypes than the ones recorded.\n";
    if (!allow_mismatch) {
        cout << "  Pass --allow_panel_mismatch to proceed anyway.\n";
        exit(-1);
    }
}
