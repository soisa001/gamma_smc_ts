#include "common.h"

#include "sys.h"
#include "io.h"
#include "flow_field.h"
#include "gamma_smc.h"
#include "data_processor.h"
#include "pair_sampling.h"
#include "screenoutput.h"
#include "cxxopts.hpp"

int main(int argc, char** argv) {
    ScreenOutput screen;
    screen.print_header();

    //
    // Print command line
    //
    screen.print_subtitle("Command line:");
    for (int i = 0; i < argc; ++i) {
        cout << argv[i] << ' ';
    }
    cout << endl << endl;


    //
    // Parse flags
    //
    cxxopts::Options options("Gamma-SMC", "Fast inference of pairwise coalescence times");

    options.add_options()
        ("i,input", "Input file", cxxopts::value<std::string>())
        ("input_format", "Input format: auto, vcf, trees, or tsz (required for tree sequence stdin)", cxxopts::value<std::string>()->default_value("auto"))
        ("allow_unphased", "Allow unphased heterozygotes (not recommended for haplotype scans)")
        ("o,output", "Output file", cxxopts::value<std::string>())
        ("recent_summary", "Write an across-pair recent-coalescence TSV", cxxopts::value<std::string>())
        ("recent_bitmatrix", "Write the per-pair recent-coalescence calls as a packed bit matrix", cxxopts::value<std::string>())
        ("recent_threshold_years", "Recent-coalescence threshold in years; repeat for several", cxxopts::value<std::vector<double>>()->default_value("4500"))
        ("recent_call", "Per-pair call rule: median, mean, or prob", cxxopts::value<std::string>()->default_value("median"))
        ("recent_call_probability", "Probability used by --recent_call prob", cxxopts::value<double>()->default_value("0.5"))
        ("no_recent_probability", "Skip the across-pair mean of P(T<t); counts only")
        ("generation_time", "Generation time in years", cxxopts::value<double>()->default_value("30"))
        ("unscaled_mutation_rate", "Per-base per-generation mutation rate used to unscale time", cxxopts::value<double>())
        ("m,scaled_mutation_rate", "Scaled mutation rate", cxxopts::value<float>())
        ("r,scaled_recombination_rate", "Scaled recombination rate", cxxopts::value<float>())
        ("t,recombination_to_mutation_ratio", "Recombination to mutation rates ratio", cxxopts::value<float>())
        ("f,flow_field", "Flow field file", cxxopts::value<std::string>())
        ("a,mask", "File of a global mask (empty for no mask)", cxxopts::value<std::string>())
        ("b,masks_per_sample", "File of masks filenames per sample (empty for no masks)", cxxopts::value<std::string>())
        ("S,samples", "Filename of a list of subset of samples to take", cxxopts::value<std::string>())
        ("T,samples_against", "Filename of a second list of subset of samples to take, to infer against first list", cxxopts::value<std::string>())
        ("w,only_within", "Apply only to haplotype pairs within each diploid")
        ("n_random_pairs", "Sample this many haplotype pairs uniformly at random", cxxopts::value<long>()->default_value("0"))
        ("pairs_seed", "Seed for --n_random_pairs", cxxopts::value<unsigned long long>()->default_value("1729"))
        ("pairs_file", "File of explicit haplotype pairs, two 0-based indices per line", cxxopts::value<std::string>())
        ("exclude_within", "Exclude within-individual pairs when sampling at random")
        ("s,output_at_stride", "Output at positions which are multiples of this number", cxxopts::value<int>()->default_value("-1"))
        ("h,output_at_hets", "Output at segregating sites", cxxopts::value<bool>()->default_value("true"))
        ("z,cache_size", "Maximum cache size in basepairs", cxxopts::value<int>()->default_value("1000"))
        ("j,threads", "Worker threads (0 = all available)", cxxopts::value<int>()->default_value("0"))
        ("pair_block", "Pairs decoded per work unit and per bit-matrix frame", cxxopts::value<long>()->default_value("256"))
        ("exp10", "accurate (default) or fast (upstream's approximation, ~1% systematic bias)", cxxopts::value<std::string>()->default_value("accurate"))
        ("exact_recent_stats", "Evaluate P(T<t) with boost::math::gamma_p per element instead of the lookup tables (validation only; very slow)")
        ("backward_alignment", "fixed (default) or legacy (upstream's one-output-position shift of the backward message)", cxxopts::value<std::string>()->default_value("fixed"))
        ("y,only_forward", "Calculate only forward pass", cxxopts::value<bool>()->default_value("false"))
        ("d,only_backward", "Calculate only backward pass", cxxopts::value<bool>()->default_value("false"))
        ("zstd_compression_level", "zstd compression level", cxxopts::value<int>()->default_value("1"))
        ("help", "Produce help message", cxxopts::value<bool>()->default_value("false"))
    ;

    auto vm = options.parse(argc, argv);

    // Check for --help or --version here, before notify
    if (vm.count("help")) {
        std::cout << options.help() << "\n";
        exit(-1);
    }

    //
    // Validate flags
    //
    float scaled_mutation_rate = -1;
    if (vm.count("scaled_mutation_rate")) {
        scaled_mutation_rate = vm["scaled_mutation_rate"].as<float>();
        if (scaled_mutation_rate < 0) {
            cout << "Error: --scaled_mutation_rate must be positive." << endl;
            exit(-1);
        }
    }

    if ((vm.count("scaled_recombination_rate") > 0) && (vm.count("recombination_to_mutation_ratio") > 0)) {
        cout << "Error: --scaled_recombination_rate and --recombination_to_mutation_ratio are mutually exclusive." << endl;
        exit(-1);
    }

    if ((vm.count("scaled_recombination_rate") == 0) && (vm.count("recombination_to_mutation_ratio") == 0)) {
        cout << "Error: Either --scaled_recombination_rate or --recombination_to_mutation_ratio must be specified." << endl;
        exit(-1);
    }

    float scaled_recombination_rate;
    if (vm.count("scaled_recombination_rate")) {
        scaled_recombination_rate = vm["scaled_recombination_rate"].as<float>();
        if (scaled_recombination_rate < 0) {
            cout << "Error: --scaled_recombination_rate must be positive." << endl;
            exit(-1);
        }
    }

    float recombination_to_mutation_ratio = -1;
    if (vm.count("recombination_to_mutation_ratio")) {
        recombination_to_mutation_ratio = vm["recombination_to_mutation_ratio"].as<float>();
        if (recombination_to_mutation_ratio < 0) {
            cout << "Error: --recombination_to_mutation_ratio must be positive." << endl;
            exit(-1);
        }
    }

    string flow_field_filename;
    if (vm.count("flow_field")) {
        flow_field_filename = vm["flow_field"].as<string>();
        if (!std::filesystem::exists(flow_field_filename)) {
            cout << boost::format("Error: Cannot open --flow_field file: %s\n") % flow_field_filename;
            exit(-1);
        }
    }

    if (vm.count("input") == 0) {
        cout << boost::format("Error: --input required.\n");
        exit(-1);
    }

    string input_filename = vm["input"].as<string>();
    if (!std::filesystem::exists(input_filename)) {
        cout << boost::format("Error: Cannot open --input file: %s\n") % input_filename;
        exit(-1);
    }

    string input_format = vm["input_format"].as<string>();
    const vector<string> valid_input_formats{"auto", "vcf", "trees", "tsz"};
    if (std::find(valid_input_formats.begin(), valid_input_formats.end(), input_format) == valid_input_formats.end()) {
        cout << "Error: --input_format must be auto, vcf, trees, or tsz." << endl;
        exit(-1);
    }
    if (input_filename == "/dev/stdin" && input_format == "auto") {
        cout << "Error: --input_format is required when reading /dev/stdin." << endl;
        exit(-1);
    }

    if (vm.count("output") == 0 && vm.count("recent_summary") == 0 && vm.count("recent_bitmatrix") == 0) {
        cout << boost::format("Error: At least one of --output, --recent_summary or --recent_bitmatrix is required.\n");
        exit(-1);
    }
    string output_filename;
    if (vm.count("output")) {
        output_filename = vm["output"].as<string>();
        auto output_directory = std::filesystem::path(output_filename).parent_path();
        if (!output_directory.empty()) {
            std::filesystem::create_directories(output_directory);
        }
    }

    //
    // Recent-coalescence configuration
    //
    const bool wants_recent_summary = vm.count("recent_summary") > 0;
    const bool wants_bitmatrix = vm.count("recent_bitmatrix") > 0;
    const bool wants_recent = wants_recent_summary || wants_bitmatrix;

    string recent_summary_filename;
    if (wants_recent_summary) {
        recent_summary_filename = vm["recent_summary"].as<string>();
        auto summary_directory = std::filesystem::path(recent_summary_filename).parent_path();
        if (!summary_directory.empty()) {
            std::filesystem::create_directories(summary_directory);
        }
    }
    string bitmatrix_filename;
    if (wants_bitmatrix) {
        bitmatrix_filename = vm["recent_bitmatrix"].as<string>();
        auto bitmatrix_directory = std::filesystem::path(bitmatrix_filename).parent_path();
        if (!bitmatrix_directory.empty()) {
            std::filesystem::create_directories(bitmatrix_directory);
        }
    }

    vector<double> threshold_years = vm["recent_threshold_years"].as<std::vector<double>>();
    recent_call_t recent_call = RECENT_CALL_MEDIAN;
    double recent_call_probability = vm["recent_call_probability"].as<double>();
    const bool accumulate_probability = (vm.count("no_recent_probability") == 0);

    if (wants_recent) {
        if (!vm.count("unscaled_mutation_rate") || vm["unscaled_mutation_rate"].as<double>() <= 0.0) {
            cout << "Error: --recent_summary/--recent_bitmatrix require a positive --unscaled_mutation_rate." << endl;
            exit(-1);
        }
        if (vm["generation_time"].as<double>() <= 0.0) {
            cout << "Error: --generation_time must be positive." << endl;
            exit(-1);
        }
        if (threshold_years.empty()) {
            cout << "Error: --recent_threshold_years needs at least one value." << endl;
            exit(-1);
        }
        for (double years : threshold_years) {
            if (years <= 0.0) {
                cout << "Error: every --recent_threshold_years value must be positive." << endl;
                exit(-1);
            }
        }

        const string call_name = vm["recent_call"].as<string>();
        if (call_name == "median") {
            recent_call = RECENT_CALL_MEDIAN;
        } else if (call_name == "mean") {
            recent_call = RECENT_CALL_MEAN;
        } else if (call_name == "prob") {
            recent_call = RECENT_CALL_PROB;
            if (recent_call_probability <= 0.0 || recent_call_probability >= 1.0) {
                cout << "Error: --recent_call_probability must be strictly between 0 and 1." << endl;
                exit(-1);
            }
        } else {
            cout << "Error: --recent_call must be median, mean or prob." << endl;
            exit(-1);
        }
    }

    if (vm.count("only_within") && vm.count("samples_against")) {
        cout << "Error: --only_within and --samples_file_against are mutually exclusive." << endl;
        exit(-1);
    }
    if (vm.count("samples_against") && !vm.count("samples")) {
        cout << "Error: --samples_file_against requires --samples_file." << endl;
        exit(-1);
    }
    if (vm.count("mask") && vm.count("masks_per_sample")) {
        cout << "Error: --mask and --masks_per_sample are mutually exclusive." << endl;
        exit(-1);
    }

    const long n_random_pairs = vm["n_random_pairs"].as<long>();
    const bool wants_pairs_file = vm.count("pairs_file") > 0;
    if (n_random_pairs < 0) {
        cout << "Error: --n_random_pairs must not be negative." << endl;
        exit(-1);
    }
    {
        int n_pair_selectors = (n_random_pairs > 0) + (wants_pairs_file ? 1 : 0)
                             + (vm.count("only_within") > 0) + (vm.count("samples_against") > 0);
        if (n_pair_selectors > 1) {
            cout << "Error: --n_random_pairs, --pairs_file, --only_within and --samples_against "
                    "are mutually exclusive." << endl;
            exit(-1);
        }
    }
    string pairs_filename;
    if (wants_pairs_file) {
        pairs_filename = vm["pairs_file"].as<string>();
        if (!std::filesystem::exists(pairs_filename)) {
            cout << boost::format("Error: Cannot open --pairs_file: %s\n") % pairs_filename;
            exit(-1);
        }
    }

    string mask_filename;
    if (vm.count("mask")) {
        mask_filename = vm["mask"].as<string>();
        if (!std::filesystem::exists(mask_filename)) {
            cout << boost::format("Error: Cannot open --mask file: %s\n") % mask_filename;
            exit(-1);
        }
    }

    string masks_per_sample_filename;
    if (vm.count("masks_per_sample")) {
        masks_per_sample_filename = vm["masks_per_sample"].as<string>();
        if (!std::filesystem::exists(masks_per_sample_filename)) {
            cout << boost::format("Error: Cannot open --masks_per_sample file: %s\n") % masks_per_sample_filename;
            exit(-1);
        }
    }

    string samples_filename;
    if (vm.count("samples")) {
        samples_filename = vm["samples"].as<string>();
        if (!std::filesystem::exists(samples_filename)) {
            cout << boost::format("Error: Cannot open --samples file: %s\n") % samples_filename;
            exit(-1);
        }
    }

    string samples_against_filename;
    if (vm.count("samples_against")) {
        samples_against_filename = vm["samples_against"].as<string>();
        if (!std::filesystem::exists(samples_against_filename)) {
            cout << boost::format("Error: Cannot open --samples_against file: %s\n") % samples_against_filename;
            exit(-1);
        }
    }

    int output_at_stride = vm["output_at_stride"].as<int>();
    bool output_at_hets = vm["output_at_hets"].as<bool>();
    if (!output_at_hets && (output_at_stride == -1)) {
        cout << "Warning: No output flags provided.\n";
    }

    bool only_forward = (vm.count("only_forward") > 0);
    bool only_backward = (vm.count("only_backward") > 0);

    int cache_size = vm["cache_size"].as<int>();
    if (cache_size <= 0) {
        cout << boost::format("Error: --cache_size must be positive.\n");
        exit(-1);
    }

    long pair_block = vm["pair_block"].as<long>();
    if (pair_block < parallel_vector_size) {
        pair_block = parallel_vector_size;
    }

    const string backward_alignment = vm["backward_alignment"].as<string>();
    if (backward_alignment != "legacy" && backward_alignment != "fixed") {
        cout << "Error: --backward_alignment must be fixed or legacy." << endl;
        exit(-1);
    }

    const string exp10_mode = vm["exp10"].as<string>();
    if (exp10_mode != "accurate" && exp10_mode != "fast") {
        cout << "Error: --exp10 must be accurate or fast." << endl;
        exit(-1);
    }

    int n_threads = vm["threads"].as<int>();
#ifdef _OPENMP
    if (n_threads <= 0) {
        n_threads = omp_get_max_threads();
    }
    omp_set_num_threads(n_threads);
#else
    if (n_threads > 1) {
        cout << "Warning: built without OpenMP; running on one thread.\n";
    }
    n_threads = 1;
#endif

    //
    // Load input file
    //
    screen.print_subtitle("Reading input file...");

    SiteMatrix input_sites;
    vector<string> sample_names;
    vector<int> samples_indices;
    vector<int> samples_against_indices;

    readVcf(
        input_filename,
        input_sites,
        sample_names,
        samples_filename,
        samples_indices,
        samples_against_filename,
        samples_against_indices,
        input_format,
        vm.count("allow_unphased") > 0
    );

    screen.print_item(boost::str(boost::format("Read %d samples.") % sample_names.size()));
    screen.print_item(boost::str(boost::format("Read %d segregating sites.") % input_sites.size()));
    screen.print_item(boost::str(
        boost::format("Genotype matrix: %.3f GB (bit-packed).") % (input_sites.bytes() / 1073741824.0)
    ));
    if (input_sites.empty()) {
        cout << "Error: Input contains no segregating SNP sites after filtering." << endl;
        exit(-1);
    }
    screen.print_done();

    //
    // Load masks
    //
    screen.print_subtitle("Reading mask(s) file...");

    vector<pair<int, int>> global_mask;
    if (vm.count("mask")) {
        readMask(mask_filename, global_mask);
    } else {
        // If no global mask is given, assume no mask
        global_mask.push_back(make_pair(0, input_sites.pos.back()+1));
    }

    unordered_map<string, vector<pair<int, int>>> mask_map;
    if (vm.count("masks_per_sample")) {
        readMasks(masks_per_sample_filename, mask_map, sample_names);
        if (mask_map.size() != sample_names.size()) {
            cout << boost::format("Error: Expected one mask for each of %d samples, but read %d.\n")
                    % sample_names.size() % mask_map.size();
            exit(-1);
        }
        screen.print_item(boost::str(boost::format("Read %d masks.") % mask_map.size()));
    }

    screen.print_done();

    //
    // Process segments
    //
    screen.print_subtitle("Creating segments...");

    DataProcessor data_processor(
        input_sites,
        sample_names,
        global_mask,
        mask_map,
        output_at_stride,
        cache_size,
        output_at_hets
    );

    screen.print_item(boost::str(
        boost::format("Created %d segments, with %d output positions.") % data_processor._n_segments % data_processor._seq_length
    ));

    screen.print_done();

    //
    // Estimated scaled mutation rate if needed
    //
    screen.print_subtitle("Calculating rates...");
    if (vm.count("scaled_mutation_rate") == 0) {
        screen.print_item("Estimating mutation rate...");
        scaled_mutation_rate = data_processor.calculate_heterozygosity();
    }

    if (vm.count("recombination_to_mutation_ratio")) {
        scaled_recombination_rate = scaled_mutation_rate * recombination_to_mutation_ratio;
    }

    screen.print_item(boost::str(boost::format("Scaled mutation rate: %f") % scaled_mutation_rate));
    screen.print_item(boost::str(boost::format("Scaled recombination rate: %f") % scaled_recombination_rate));
    screen.print_done();

    //
    // Create a list of pairs to work on
    //
    vector<pair<int, int>> haplotype_pairs;
    uint n_samples = sample_names.size();
    const int n_haplotypes = (int) (2 * n_samples);

    if (n_random_pairs > 0) {
        sample_random_pairs(
            n_haplotypes,
            n_random_pairs,
            vm["pairs_seed"].as<unsigned long long>(),
            vm.count("exclude_within") > 0,
            haplotype_pairs
        );
    } else if (wants_pairs_file) {
        read_pairs_file(pairs_filename, n_haplotypes, haplotype_pairs);
    } else if (vm.count("only_within")) {
        for (uint i = 0; i < n_samples; i++) {
            haplotype_pairs.push_back(make_pair(2*i, 2*i+1));
        }
    } else {
        if (samples_against_indices.size()) {
            for (int i : samples_indices) {
                for (int j : samples_against_indices) {
                    haplotype_pairs.push_back(make_pair(2 * i, 2 * j));
                    haplotype_pairs.push_back(make_pair(2 * i, 2 * j + 1));
                    haplotype_pairs.push_back(make_pair(2 * i + 1, 2 * j));
                    haplotype_pairs.push_back(make_pair(2 * i + 1, 2 * j + 1));
                }
            }
        } else {
            for (uint i = 0; i < n_samples; i++) {
                haplotype_pairs.push_back(make_pair(2 * i, 2 * i + 1));
                for (uint j = i+1; j < n_samples; j++) {
                    haplotype_pairs.push_back(make_pair(2 * i, 2 * j));
                    haplotype_pairs.push_back(make_pair(2 * i, 2 * j + 1));
                    haplotype_pairs.push_back(make_pair(2 * i + 1, 2 * j));
                    haplotype_pairs.push_back(make_pair(2 * i + 1, 2 * j + 1));
                }
            }
        }
    }

    if (haplotype_pairs.empty()) {
        cout << "Error: no haplotype pairs selected." << endl;
        exit(-1);
    }

    screen.print_item(boost::str(
        boost::format("Applying to %d haplotype pairs") % haplotype_pairs.size()
    ));
    screen.print_item(boost::str(
        boost::format("Using %d thread(s), %d pairs per work unit") % n_threads % pair_block
    ));

    //
    // Memory budget. Per-thread scratch scales with the number of output
    // positions and with the number of segments, both of which explode when
    // --output_at_hets is left on for a large panel, so say so up front rather
    // than letting the run die in the allocator an hour later.
    //
    {
        const double bytes_per_thread_posteriors =
            2.0 * (double) data_processor._seq_length * parallel_vector_size * sizeof(float);
        const double bytes_per_thread_segment_types =
            (double) data_processor._n_segments * parallel_vector_size;
        const double bytes_per_thread_n_called = mask_map.empty()
            ? 0.0
            : (double) data_processor._n_segments * parallel_vector_size * sizeof(int32_t);
        const long effective_block = (output_filename.size() > 0) ? parallel_vector_size : pair_block;
        const double bytes_per_thread_bits = bitmatrix_filename.empty()
            ? 0.0
            : 2.0 * (double) threshold_years.size()
              * (double) (effective_block / parallel_vector_size)
              * (double) data_processor._seq_length;
        const double bytes_per_thread_accumulators = wants_recent
            ? (double) data_processor._seq_length * (threshold_years.size() * 16.0 + 16.0)
            : 0.0;

        const double per_thread = bytes_per_thread_posteriors + bytes_per_thread_segment_types
            + bytes_per_thread_n_called + bytes_per_thread_bits + bytes_per_thread_accumulators;
        const double shared = input_sites.bytes()
            + (mask_map.empty()
               ? (double) data_processor._n_segments * parallel_vector_size * sizeof(int32_t)
               : 0.0)
            + 489.6e6;   // flow-field cache, six flat tables

        screen.print_item(boost::str(
            boost::format("Memory estimate: %.2f GB shared + %.2f GB x %d threads = %.2f GB")
            % (shared / 1073741824.0)
            % (per_thread / 1073741824.0)
            % n_threads
            % ((shared + per_thread * n_threads) / 1073741824.0)
        ));
        if (output_at_hets && data_processor._seq_length > 1000000) {
            cout << boost::format(
                "Warning: --output_at_hets is on with %ld output positions. For a whole-genome "
                "scan use --output_at_hets=false --output_at_stride 1000 to cut per-thread "
                "memory and output volume by an order of magnitude.\n"
            ) % (long) data_processor._seq_length;
        }
        if (output_filename.size() > 0) {
            const double raw_bytes = 2.0 * (double) data_processor._seq_length
                                   * (double) haplotype_pairs.size() * sizeof(float);
            screen.print_item(boost::str(
                boost::format("Raw posterior output before compression: %.2f GB")
                % (raw_bytes / 1073741824.0)
            ));
            if (raw_bytes > 50e9) {
                cout << "Warning: --output at this scale writes a very large file; "
                        "--recent_bitmatrix stores one bit per pair and position instead.\n";
            }
        }
    }
    screen.print_done();


    //
    // Load flow field file
    //
    vector<float> mean_grid_def;
    vector<float> cv_grid_def;
    vector<float> flow_field_unravelled;

    if (vm.count("flow_field")) {
        read_flow_field_raw(
            flow_field_filename,
            mean_grid_def,
            cv_grid_def,
            flow_field_unravelled
        );
    } else {
        read_flow_field_default(
            mean_grid_def,
            cv_grid_def,
            flow_field_unravelled
        );
    }

    //
    // Prepare output files
    //
    // TODO: Check errors
    ofstream* output_file_raw_meta = NULL;
    ofstream* output_file_raw = NULL;
    if (output_filename.size() > 0) {
        output_file_raw = new ofstream(output_filename, ios_base::out | ios_base::binary);
        output_file_raw_meta = new ofstream(output_filename + ".meta", ios_base::out);
    }
    ofstream* recent_summary_file = NULL;
    if (recent_summary_filename.size() > 0) {
        recent_summary_file = new ofstream(recent_summary_filename, ios_base::out);
    }
    ofstream* bitmatrix_file = NULL;
    if (bitmatrix_filename.size() > 0) {
        bitmatrix_file = new ofstream(bitmatrix_filename, ios_base::out | ios_base::binary);
    }

    //
    // Build the recent-coalescence thresholds and lookup tables
    //
    vector<RecentThreshold> thresholds;
    RecentProbabilityTable probability_table;
    double two_ne_generations = -1.0;

    if (wants_recent) {
        screen.print_subtitle("Building recent-coalescence tables...");
        const double mutation_rate = vm["unscaled_mutation_rate"].as<double>();
        const double generation_time = vm["generation_time"].as<double>();
        two_ne_generations = scaled_mutation_rate / (2.0 * mutation_rate);

        for (double years : threshold_years) {
            RecentThreshold threshold;
            threshold.years = years;
            threshold.generations = years / generation_time;
            threshold.scaled = threshold.generations / two_ne_generations;
            threshold.build(recent_call, recent_call_probability);
            thresholds.push_back(std::move(threshold));
        }

        screen.print_item(boost::str(boost::format("2Ne: %.1f generations") % two_ne_generations));
        for (const auto& threshold : thresholds) {
            screen.print_item(boost::str(
                boost::format("Threshold %.0f years = %.1f generations = %.6g coalescent units")
                % threshold.years % threshold.generations % threshold.scaled
            ));
        }
        screen.print_item(boost::str(
            boost::format("Call rule: %s%s") % recent_call_name(recent_call)
            % (recent_call == RECENT_CALL_PROB
               ? boost::str(boost::format(" (p = %.4g)") % recent_call_probability)
               : string(""))
        ));

        if (accumulate_probability) {
            probability_table.build();
        }

        // Runs the SIMD kernels themselves against Boost, so a bad grid or a
        // mistake in the deviance shows up here rather than as a quietly wrong
        // statistic in the output.
        const RecentTableAccuracy accuracy = self_check_recent_tables(
            probability_table, thresholds.front(), accumulate_probability,
            recent_call, recent_call_probability
        );
        if (accumulate_probability) {
            screen.print_item(boost::str(
                boost::format("P(T<t) table: %.1f MB, max abs error vs Boost %.2e")
                % (probability_table._table.size() * sizeof(float) / 1048576.0)
                % accuracy.max_probability_error
            ));
            if (accuracy.max_probability_error > 1e-2) {
                cout << "Warning: the P(T<t) table is less accurate than expected; "
                        "treat mean_p_* columns with care.\n";
            }
        }
        if (recent_call != RECENT_CALL_MEAN) {
            screen.print_item(boost::str(
                boost::format("Call table: %ld probes, disagreement with Boost %.2e")
                % accuracy.n_samples % accuracy.call_disagreement
            ));
        }
        screen.print_done();
    }

    //
    // Construct flow field
    //
    unique_ptr<FlowField> FF(new FlowField(mean_grid_def, cv_grid_def, flow_field_unravelled));

    //
    // Construct flow field cache
    //
    screen.print_subtitle("Building flow field cache...");

    unique_ptr<FlowFieldCache> FFC(new FlowFieldCache(
        scaled_recombination_rate,
        scaled_mutation_rate,
        move(FF),
        cache_size,
        true                // entropy clipping
    ));

    screen.print_done();

    //
    // Run
    //
    screen.print_subtitle("Running...");
    CachedPairwiseGammaSMC PPC(
        input_sites,
        haplotype_pairs,
        scaled_recombination_rate,
        scaled_mutation_rate,
        move(FFC),
        data_processor,
        output_at_stride,
        output_at_hets,
        only_forward,
        only_backward,
        output_file_raw_meta,
        output_file_raw,
        vm["zstd_compression_level"].as<int>(),
        recent_summary_file,
        thresholds,
        accumulate_probability ? &probability_table : NULL,
        recent_call,
        accumulate_probability,
        two_ne_generations,
        bitmatrix_file,
        bitmatrix_filename,
        pair_block,
        n_threads
    );

    PPC._accurate_exp10 = (exp10_mode == "accurate");
    PPC._fix_backward_alignment = (backward_alignment == "fixed");
    PPC._exact_recent_stats = (vm.count("exact_recent_stats") > 0);
    PPC._recent_call_probability = recent_call_probability;
    if (!PPC._accurate_exp10 || !PPC._fix_backward_alignment) {
        screen.print_item(boost::str(boost::format(
            "Reproducing upstream numerics: exp10=%s, backward_alignment=%s"
        ) % exp10_mode % backward_alignment));
    }
    if (PPC._exact_recent_stats) {
        screen.print_item("Using exact boost::math::gamma_p instead of the lookup tables.");
    }

    PPC.calculate_posteriors();

    //
    // Close output files
    //
    if (output_file_raw != NULL) {
        output_file_raw->close();
        output_file_raw_meta->close();
    }
    if (recent_summary_file != NULL) {
        recent_summary_file->close();
    }
    if (bitmatrix_file != NULL) {
        bitmatrix_file->close();
        PPC.write_bitmatrix_meta(bitmatrix_filename);
    }

    double total_processing_time = PPC._timer_emissions + PPC._timer_forward + PPC._timer_backward;  // Excludes output time
    double total_basepairs = data_processor._segments.back().pos + data_processor._segments.back().length - data_processor._segments.front().pos;
    double time_per_bp_per_pair = total_processing_time / total_basepairs / haplotype_pairs.size();

    // The per-pass timers sum CPU time across workers, so they no longer add up
    // to elapsed time; both are reported to keep the distinction visible.
    screen.print_item(boost::str(boost::format("Emissions preparation time:\t%.3f secs (CPU)") % PPC._timer_emissions));
    screen.print_item(boost::str(boost::format("Forward pass time:\t\t%.3f secs (CPU)") % PPC._timer_forward));
    screen.print_item(boost::str(boost::format("Backward pass time:\t\t%.3f secs (CPU)") % PPC._timer_backward));
    screen.print_item(boost::str(boost::format("Output time:\t\t%.3f secs (CPU)") % PPC._timer_output));
    screen.print_item(boost::str(boost::format("Decoding wall time:\t\t%.3f secs") % PPC._timer_wall));
    screen.print_done();

    screen.print_subtitle("Summary:");

    screen.print_item(boost::str(boost::format("Time per Gbp per pair\t%.3f secs (CPU)") % (time_per_bp_per_pair * 1e9)));
    screen.print_item(boost::str(boost::format("Overall time:\t\t%.3f secs") % mp_cputime()));
    screen.print_item(boost::str(boost::format("Peak memory:\t\t%.3f GB") % (mp_peakrss() / 1024.0 / 1024.0 / 1024.0)));
    screen.print_done();

    return 0;
}
