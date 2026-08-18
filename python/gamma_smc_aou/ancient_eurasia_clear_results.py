"""Clear, additive results for the AncientEurasia 9K19 Han experiment.

This module is deliberately post-processing only.  It verifies and reuses the
frozen selected simulations/decodes, a separately completed survival-only
neutral bank, and the published final-AF-matched analysis. It exposes no
generation route.

The primary estimand is the focal, within-diploid overall
``P(TMRCA < x)`` in selected panels versus a neutral bank that has the same
single-copy archaic origin and pulse-survival event but no selection. It is
conditioned on strict present Han census segregation, with no terminal AF target,
band, or selected-AF matching beyond 0 < AF < 1. The selected
alt/alt-versus-ref/ref contrast is descriptive; the final-AF-matched neutral
contrast is a labeled sensitivity analysis rather than the primary null.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import re
import shutil
import struct
from collections.abc import Mapping, Sequence
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Any

os.environ["MPLBACKEND"] = "Agg"
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tskit
from matplotlib.lines import Line2D

SCHEMA_VERSION = "gamma-smc.ancient-eurasia-clear-results/v1"
ANALYSIS_SCHEMA_VERSION = "gamma-smc.focused-selection-analysis/v3"
FOCUSED_SIMULATION_SCHEMA_VERSION = "gamma-smc.focused-selection-simulation/v5"
FOCUSED_DECODE_SCHEMA_VERSION = "gamma-smc.focused-selection-decode/v1"
SURVIVAL_PLAN_SCHEMA_VERSION = "gamma-smc.ancient-eurasia-survival-neutral-plan/v2"
SURVIVAL_SIMULATION_SCHEMA_VERSION = (
    "gamma-smc.ancient-eurasia-survival-neutral-simulation/v2"
)
SURVIVAL_DECODE_SCHEMA_VERSION = "gamma-smc.ancient-eurasia-survival-neutral-decode/v2"
SURVIVAL_RESULTS_SCHEMA_VERSION = (
    "gamma-smc.ancient-eurasia-survival-neutral-results/v2"
)
SURVIVAL_DECODER_RUN_SCHEMA_VERSION = (
    "gamma-smc.ancient-eurasia-survival-neutral-decoder-run/v2"
)
SURVIVAL_FOCAL_IDENTITY_SCHEMA = "gamma-smc.selected-focal-identity/v1"
SURVIVAL_FOCAL_OVERLAY_PATCH_SCHEMA = (
    "gamma-smc.ancient-eurasia-reserved-focal-overlay-patch/v1"
)
SURVIVAL_INDIVIDUAL_LOCATION_SCHEMA = (
    "gamma-smc.ancient-eurasia-individual-location-canonicalization/v1"
)
FOCUSED_SIMULATION_RECEIPT_FIELDS = frozenset(
    {
        "schema",
        "status",
        "contract",
        "contract_sha256",
        "elapsed_seconds",
        "attempts_completed",
        "search_budget",
        "sample_alt_count",
        "sample_af",
        "genotype_counts",
        "provenance",
        "outputs",
    }
)
FOCUSED_DECODE_RECEIPT_FIELDS = frozenset(
    {
        "schema",
        "status",
        "contract",
        "contract_sha256",
        "software",
        "posterior_scale",
        "posterior_dimensions",
        "focal_center",
        "outputs",
    }
)
FOCUSED_SIMULATION_CONTRACT_FIELDS = frozenset(
    {"schema", "unit", "parameters", "software", "implementation"}
)
FOCUSED_DECODE_CONTRACT_FIELDS = frozenset({"static", "raw_inputs"})
FOCUSED_SIMULATION_OUTPUT_KEYS = frozenset(
    {"tree", "pair_table", "overall_pairs", "truth_profiles", "truth_class_summaries"}
)
FOCUSED_DECODE_OUTPUT_KEYS = frozenset(
    {
        "overall_summary",
        "class_summaries",
        "pair_summaries",
        "spatial_profiles",
        "run_manifest",
        "decoder_stdout",
        "decoder_stderr",
    }
)
FOCUSED_DECODE_PROVENANCE_FIELDS = frozenset(
    {
        "decoder_path",
        "decoder_sha256",
        "implementation_sources",
        "overall_pairs_sha256",
        "pair_table_sha256",
        "settings",
        "simulation_completion_sha256",
        "tree_sha256",
    }
)
FOCUSED_DECODER_SETTINGS = {
    "backward_alignment": "fixed",
    "cache_size": 1_000,
    "exp10": "accurate",
    "generation_time_years": 25.0,
    "input_format": "trees_via_streamed_vcf",
    "output_at_hets": False,
    "output_at_stride": 10_000,
    "pair_block": 256,
    "pair_selector": "explicit_within_individual_pairs",
    "recent_call": "median",
    "recent_call_probability": 0.5,
    "recombination_to_mutation_ratio": 0.8,
    "scaled_mutation_rate": 0.00031499999999999996,
    "threads": 1,
    "thresholds_years": [1_000, 4_500, 10_000, 20_000, 30_000, 40_000, 50_000],
    "unscaled_mutation_rate": 1.25e-8,
}
FOCUSED_SIMULATION_PARAMETER_FIELDS = frozenset(
    {
        "candidate_pool_diploids",
        "candidate_pool_frequency_gate",
        "focal_semantic_validation",
        "han_selection_calibration",
        "introgression_mutation_age_generations",
        "introgression_source_presence_epsilon",
        "neutral_ancestry_batch_size",
        "neutral_ancestry_model",
        "neutral_ancestry_model_keywords",
        "neutral_candidate_index",
        "neutral_eas_ascertainment_length_bp",
        "neutral_han_ascertainment_length_bp",
        "neutral_pool_diploids",
        "neutral_recent_dtwf_duration_generations",
        "selected_pool_diploids",
        "slim_burn_in",
        "slim_scaling_factor",
    }
)
FOCUSED_NEUTRAL_PROVENANCE_FIELDS = frozenset(
    {
        "accepted_attempt_zero_based",
        "accepted_batch_zero_based",
        "ancestry_batch_size",
        "ancestry_seed",
        "archaic_specific_no_ils",
        "ascertainment_length_bp",
        "candidate_index",
        "candidate_pool_diploids",
        "candidate_pool_frequency_gate",
        "candidate_seed",
        "candidate_weight",
        "coalescent_approximation",
        "engine",
        "entry_route",
        "focal_semantic_validation",
        "matched_selected_neutral_ascertainment",
        "model",
        "model_keywords",
        "msprime_version",
        "mutation_origin_population",
        "mutation_seed",
        "mutation_time_generations_ago",
        "n_eligible_branch_rectangles_in_batch",
        "n_introgressed_descendant_lineages",
        "neutral_focal_identity",
        "neutral_origin_matching",
        "panel_seed",
        "pre_crop_focal_position",
        "recent_dtwf_duration_generations",
        "sample_panel",
        "window_length_bp",
    }
)
FOCUSED_SELECTED_PROVENANCE_FIELDS = frozenset(
    {
        "accepted_draw_zero_based",
        "accepted_seed",
        "candidate_pool_frequency_gate",
        "cumulative_timeout_seconds",
        "draw_timeout_seconds",
        "engine",
        "focal_semantic_validation",
        "han_direct_production_v3",
        "matched_selected_neutral_ascertainment",
        "pyslim_version",
        "sample_panel",
        "selected_focal_identity",
        "selection_endpoint",
        "slim_burn_in",
        "slim_path",
        "slim_scaling_factor",
        "slim_sha256",
        "stdpopsim_version",
        "sweep",
    }
)
FOCUSED_SIMULATION_OUTPUT_PATHS = {
    "tree": "simulation.trees",
    "pair_table": "sample_manifest.tsv",
    "overall_pairs": "pairs/overall.pairs.tsv",
    "truth_profiles": "truth_profiles.tsv.gz",
    "truth_class_summaries": "truth_class_summaries.tsv",
}
FOCUSED_DECODE_OUTPUT_PATHS = {
    "overall_summary": "overall.summary.tsv",
    "class_summaries": "class_summaries.tsv",
    "pair_summaries": "pair_summaries.tsv.gz",
    "spatial_profiles": "spatial_profiles.tsv.gz",
    "run_manifest": "overall.summary.tsv.run.json",
    "decoder_stdout": "decoder.stdout.zst",
    "decoder_stderr": "decoder.stderr.zst",
}
FOCUSED_DECODE_TABULAR_OUTPUT_ROWS = {
    "class_summaries": 28,
    "pair_summaries": 700,
    "spatial_profiles": 28_000,
}
FOCUSED_CLASS_SUMMARY_COMMON_COLUMNS = (
    "unit_id",
    "demography_id",
    "demography_kind",
    "population",
    "simulation_class",
    "selection_coefficient",
    "target_allele_frequency",
    "replicate_index",
    "seed",
    "source",
    "genotype_class",
    "threshold_years",
    "region_mean_p_tmrca_lt_threshold",
    "focal_mean_p_tmrca_lt_threshold",
    "region_mean_tmrca_generations",
    "focal_mean_tmrca_generations",
    "n_pairs",
)
FOCUSED_TRUTH_CLASS_SUMMARY_COLUMNS = (
    *FOCUSED_CLASS_SUMMARY_COMMON_COLUMNS,
    "focal_output_position_0based",
    "focal_offset_bp",
)
FOCUSED_GAMMA_CLASS_SUMMARY_COLUMNS = (
    *FOCUSED_CLASS_SUMMARY_COMMON_COLUMNS,
    "n_valid_posterior_cells_region",
    "n_valid_pairs_focal",
    "focal_output_position_0based",
    "focal_offset_bp",
)
FOCUSED_DECODE_TABULAR_OUTPUT_COLUMNS = {
    "class_summaries": list(FOCUSED_GAMMA_CLASS_SUMMARY_COLUMNS),
    "pair_summaries": [
        "unit_id",
        "demography_id",
        "demography_kind",
        "population",
        "simulation_class",
        "selection_coefficient",
        "target_allele_frequency",
        "replicate_index",
        "seed",
        "source",
        "pair_index",
        "gamma_smc_haplotype_0",
        "gamma_smc_haplotype_1",
        "genotype_class",
        "threshold_years",
        "region_mean_p_tmrca_lt_threshold",
        "focal_p_tmrca_lt_threshold",
        "region_mean_tmrca_generations",
        "focal_tmrca_generations",
        "n_valid_positions",
        "focal_output_position_0based",
        "focal_offset_bp",
    ],
    "spatial_profiles": [
        "unit_id",
        "demography_id",
        "demography_kind",
        "population",
        "simulation_class",
        "selection_coefficient",
        "target_allele_frequency",
        "replicate_index",
        "seed",
        "source",
        "genotype_class",
        "threshold_years",
        "position_0based",
        "mean_p_tmrca_lt_threshold",
        "mean_tmrca_generations",
        "n_valid_pairs",
    ],
}
FOCUSED_SIMULATION_SOFTWARE = {
    "msprime": "1.4.2",
    "numpy": "2.4.6",
    "pandas": "3.0.3",
    "pyslim": "1.1.1",
    "python": "3.11.15",
    "stdpopsim": "0.3.0",
    "tskit": "1.0.3",
}
FOCUSED_DECODE_SOFTWARE = {
    "numpy": "2.4.6",
    "pandas": "3.0.3",
    "scipy": "1.17.1",
    "zstandard": "0.25.0",
}
FOCUSED_DECODE_PARAMETERS = {
    "focal_position_0based": 5_000_000.0,
    "generation_time_years": 25.0,
    "thresholds_years": [
        1_000.0,
        4_500.0,
        10_000.0,
        20_000.0,
        30_000.0,
        40_000.0,
        50_000.0,
    ],
    "unscaled_mutation_rate_per_bp_per_generation": 1.25e-8,
}
FOCUSED_DECODE_POSTERIOR_SCALE = {
    "scaled_mutation_rate": 0.000315,
    "two_ne_generations": 12_600.000000000002,
    "unscaled_mutation_rate": 1.25e-8,
}
FOCUSED_DECODE_POSTERIOR_DIMENSIONS = {
    "chunk_size": 8,
    "num_output_positions": 1_000,
    "num_pairs": 100,
    "total_cells": 100_000,
    "valid_cells": 100_000,
}
FOCUSED_DECODE_FOCAL_CENTER = {
    "offset_bp": 0.0,
    "output_position_0based": 5_000_000,
    "output_position_index": 500,
    "requested_focal_position_0based": 5_000_000.0,
}
FOCUSED_POSTPROCESSOR_SHA256 = (
    "7d5e64f3c56e7d18277ade91cbac9183044b72593681cec225da4ec6dad36b50"
)
FOCUSED_DECODE_SOURCE_HASHES = {
    "python/gamma_smc_aou/decoder.py": (
        "c6c7fd339a3545f328bed8e4c7e311eccf46e71459cd3bde9180e861cd23e41e"
    ),
    "python/gamma_smc_aou/eas_sweep_models.py": (
        "365854b146d4512b6443d8dc28abaae4e53d7ee0ad5403d03536e06e8b016ae6"
    ),
    "python/gamma_smc_aou/focused_selection_campaign.py": (
        "20d7a2123e1096bbb9cf4d84a9dff66e2060733152f9bdbcd307fdb7b3aead5e"
    ),
    "python/gamma_smc_aou/focused_selection_decode.py": FOCUSED_POSTPROCESSOR_SHA256,
    "python/gamma_smc_aou/tree_sequence.py": (
        "9604ba949220388c642ea06f98d8500c2773e1a170cec43ec33f1301b981923a"
    ),
}
FOCUSED_NEUTRAL_SIMULATION_SOURCE_HASHES = {
    "python/gamma_smc_aou/eas_sweep_analysis.py": (
        "e8f30519c986416daf4d4bb825aef459d226b8ade90606286d96bed603d4abe6"
    ),
    "python/gamma_smc_aou/eas_sweep_models.py": (
        "365854b146d4512b6443d8dc28abaae4e53d7ee0ad5403d03536e06e8b016ae6"
    ),
    "python/gamma_smc_aou/focused_selection_simulation.py": (
        "29b52d52adb76d2aa195fb0900baef829f92f1ed0b7b3edbff93cadce0f26448"
    ),
}
FOCUSED_SELECTED_SIMULATION_SOURCE_HASHES = {
    "focused_selection_EAS_sim/calibration/han_direct_production_v3/han_direct_v3_attempt_seed_plan.tsv": "04017a12d2a39709777be7fc9767fcf6834bf714519296bf3a7f0878e536f376",
    "focused_selection_EAS_sim/calibration/han_direct_production_v3/han_direct_v3_excluded_seed_inventory.json": "5ffd2085d4e9598855ff000a7486b6e9bee45e5b9076bd70357ca346095149b5",
    "focused_selection_EAS_sim/calibration/han_direct_production_v3/han_direct_v3_excluded_seed_registry.tsv": "85035c70242f16b614a43dedf8d04020f5b185feabf9ade85370d0867ee067be",
    "focused_selection_EAS_sim/calibration/han_direct_production_v3/han_direct_v3_manifest.json": "7ba2be1875ea804c1cdbec4ae3b9fa5dffb8307799d6af0c810ccc9df5c40e79",
    "python/gamma_smc_aou/eas_sweep_analysis.py": "e8f30519c986416daf4d4bb825aef459d226b8ade90606286d96bed603d4abe6",
    "python/gamma_smc_aou/eas_sweep_models.py": "365854b146d4512b6443d8dc28abaae4e53d7ee0ad5403d03536e06e8b016ae6",
    "python/gamma_smc_aou/eas_sweep_study.py": "caab1f7d682d5d70b39fd39196330bb2aff025b6a974f23d8b1f31b2157e7530",
    "python/gamma_smc_aou/focused_han_direct_v3.py": "4a471605fc5e95e993acc94902fd7d9514eedc74b786767e5d69e872245f1d00",
    "python/gamma_smc_aou/focused_han_fixed_v2_calibration.py": "dfa1abc5887b6ecbd1aa3c9a698bc9db16d3e9aebbb7edb0981584ccfc909c18",
    "python/gamma_smc_aou/focused_selection_simulation.py": "29b52d52adb76d2aa195fb0900baef829f92f1ed0b7b3edbff93cadce0f26448",
    "scripts/run_focused_han_direct_v3.py": "a0346b4c9eb1d80f229e0758c3046df6a60e14b0d713868220f03c5fe048703e",
    "tests/test_focused_han_direct_v3.py": "c903e689414f3e69a59b6363b3e96b5aa5bd2b8eda20bba53d2f526ef1dc1b51",
}
DEMOGRAPHY_ID = "ancient_eurasia_han_introgression"
SOURCE_MODEL = "stdpopsim AncientEurasia_9K19"
TARGET_AF = 0.20
SELECTION_COEFFICIENT = 0.01
SELECTED_REPLICATES = 10
SURVIVAL_NEUTRAL_REPLICATES = 100
MATCHED_AF_NEUTRAL_REPLICATES = 100
SURVIVAL_BASE_SEED = 20_240_523
SURVIVAL_SEED_DERIVATION = (
    "1 + int.from_bytes(sha256(f'{base_seed}:{study_id}:{label}'.encode('utf-8'))[:8], "
    "'big') % (2**31 - 2)"
)
SAMPLE_DIPLOIDS = 100
CANDIDATE_POOL_DIPLOIDS = 500
PRESENT_HAN_UNSCALED_NE = 6_300
PRESENT_HAN_SCALED_CENSUS_DIPLOIDS = 1_260
PRESENT_HAN_SCALED_CENSUS_GENOMES = 2_520
FOCAL_POSITION_BP = 5_000_000
FOCAL_SITE_ID = "eas_selected_site"
MUTATION_ORIGIN_GENERATIONS = 2_400.0
INTROGRESSION_PULSE_GENERATIONS = 2_272.0
SELECTION_END_GENERATIONS = 1_780.0
GENERATION_TIME_YEARS = 25.0
THRESHOLDS_YEARS = (1_000, 4_500, 10_000, 20_000, 30_000, 40_000, 50_000)
DISPLAY_THRESHOLD_YEARS = 50_000
PRIMARY_ESTIMAND = "survival_only_overall_vs_overall"
SENSITIVITY_ESTIMAND = "final_af_matched_alt_alt_minus_ref_ref"
DEFAULT_WORKERS = 4
MAX_WORKERS = 24
EXPECTED_GAMMA_BINARY_IDENTITY = "bin/gamma_smc"
EXPECTED_SURVIVAL_PRODUCER_PYTHON_IDENTITY = ".venv-wsl/bin/python"
EXPECTED_GAMMA_BINARY_SHA256 = (
    "057774590af12efbbedc9d1ffb736747c0dba2604ef0028dcf9b65e275e70a22"
)
EXPECTED_SLIM_BINARY_IDENTITY = ".native-stdpopsim/bin/slim"
EXPECTED_SLIM_BINARY_SHA256 = (
    "67409ff190808c2967c949cae6697bf235010dbe832389cb4976d9089b5f5e1d"
)
EXPECTED_SURVIVAL_GENERATOR_SOURCE_SHA256 = (
    "6b30517b768656f0de0ce7aab72e7e12ec0e92cdb00bbf65a51c01186fed3c33"
)

DEFAULT_CAMPAIGN_ROOT = Path("focused_selection_EAS_sim")
DEFAULT_OUTPUT_DIR = DEFAULT_CAMPAIGN_ROOT / "results/ancient_eurasia_af20_s0p01_clear"
DEFAULT_GAMMA_ANALYSIS_DIR = DEFAULT_CAMPAIGN_ROOT / "results/analysis"
DEFAULT_TRUTH_ANALYSIS_DIR = DEFAULT_CAMPAIGN_ROOT / "results/analysis_tree_truth"
DEFAULT_SURVIVAL_RESULTS_DIR = (
    DEFAULT_CAMPAIGN_ROOT / "ancient_eurasia_single_origin_survival_v2/results"
)
DEFAULT_TWO_EPOCH_TRUTH_METRICS = Path(
    "sim_results/two_epoch_growth_s0p05_n2000_mu1p29e9/metrics.json"
)
DEFAULT_TWO_EPOCH_GAMMA_METRICS = Path(
    "sim_results/gamma_smc_container_stride1000/decoded_study_metrics.json"
)
DEFAULT_TWO_EPOCH_GAMMA_RESULTS = Path(
    "sim_results/gamma_smc_container_stride1000/RESULTS.md"
)

TABLE_FILENAMES = {
    "replicate_statistics": "replicate_statistics.tsv.gz",
    "selected_pvalues": "selected_vs_neutral_pvalues.tsv.gz",
    "primary_omnibus": "primary_7_threshold_omnibus.tsv.gz",
    "conditional_power": "conditional_power_summary.tsv",
}
PLOT_STEMS = (
    "tree_truth_pulse_aligned_null",
    "gamma_smc_pulse_aligned_null",
    "threshold_scores_pvalues_power",
    "timing_endpoints_two_epoch_diagnostic",
)
TEXT_OUTPUTS = (
    "unit_scores.tsv",
    "selected_pointwise_pvalues.tsv",
    "selected_minp.tsv",
    "power_summary.tsv",
    "endpoint_af.tsv",
    "input_cache_audit.tsv",
    "comparison_design.tsv",
    "RUN_RESULTS.md",
)

SURVIVAL_SCORE_COLUMNS = (
    "unit_id",
    "replicate_index",
    "seed",
    "source",
    "threshold_years",
    "overall_focal_p_tmrca_lt_threshold",
    "final_pool_alt_count",
    "final_pool_total_count",
    "final_pool_af",
    "sample_alt_count",
    "sample_total_count",
    "sample_af",
    "sample_detected",
)
SURVIVAL_ENDPOINT_COLUMNS = (
    "unit_id",
    "replicate_index",
    "seed",
    "simulation_completion_sha256",
    "final_census_alt_count",
    "final_census_total_count",
    "final_census_af",
    "final_census_segregating",
    "pool_diploids",
    "pool_alt_count",
    "pool_total_count",
    "pool_af",
    "pool_detected",
    "pool_fixed",
    "sample_diploids",
    "sample_alt_count",
    "sample_total_count",
    "sample_af",
    "sample_detected",
    "sample_fixed",
    "panel_seed",
)
SURVIVAL_EXECUTION_COLUMNS = (
    "unit_id",
    "replicate_index",
    "seed",
    "panel_seed",
    "scenario",
    "selection_coefficient",
    "mutation_origin_population",
    "mutation_age_generations",
    "pulse_generations",
    "present_conditioning",
    "pool_frequency_conditioning",
    "terminal_af_target_or_band_conditioning",
    "sample_detection_conditioning",
    "pool_diploids",
    "sample_diploids",
)
SURVIVAL_SIMULATION_OUTPUT_PATHS = {
    "simulation_contract": "simulation_contract.json",
    "tree": "simulation.trees",
    "sample_manifest": "sample_manifest.tsv",
    "slim_log": "slim.tsv",
    "truth_profiles": "truth_profiles.tsv.gz",
    "truth_class_summaries": "truth_class_summaries.tsv",
    "truth_unit_scores": "truth_unit_scores.tsv",
}
SURVIVAL_DECODE_OUTPUT_PATHS = {
    "decode_contract": "decode_contract.json",
    "decoder_run": "decoder_run.json",
    "overall_pairs": "overall.pairs.tsv",
    "overall_summary": "overall.summary.tsv",
    "posterior": "posterior.zst",
    "posterior_metadata": "posterior.zst.meta",
    "gamma_unit_scores": "gamma_unit_scores.tsv",
}
SURVIVAL_RESULTS_OUTPUT_PATHS = {
    "endpoint_summary": "endpoint_summary.tsv",
    "replicate_scores": "replicate_scores.tsv",
    "simulation_status": "simulation_status.tsv",
    "decode_status": "decode_status.tsv",
}
SURVIVAL_SOURCE_PATHS = (
    "python/gamma_smc_aou/ancient_eurasia_survival_neutral.py",
    "python/gamma_smc_aou/focused_selection_simulation.py",
    "python/gamma_smc_aou/focused_selection_campaign.py",
    "python/gamma_smc_aou/focused_selection_decode.py",
    "python/gamma_smc_aou/eas_sweep_models.py",
    "python/gamma_smc_aou/eas_sweep_analysis.py",
    "python/gamma_smc_aou/decoder.py",
    "python/gamma_smc_aou/tree_sequence.py",
    "python/gamma_smc_aou/han_s001_tmrca_simulation.py",
)
SURVIVAL_SAMPLE_MANIFEST_COLUMNS = (
    "vcf_diploid_index",
    "tree_sequence_individual_id",
    "sample_node_0",
    "sample_node_1",
    "gamma_smc_haplotype_0",
    "gamma_smc_haplotype_1",
    "focal_position_0based",
    "focal_selected_allele_index",
    "focal_selected_allele_count",
    "genotype_class",
)
SURVIVAL_SIMULATION_PARAMETERS = {
    "sequence_length_bp": 10_000_000,
    "focal_position_0based": 5_000_000,
    "mutation_rate": 1.25e-8,
    "recombination_rate": 1e-8,
    "slim_scaling_factor": 5.0,
    "slim_burn_in": 0.1,
    "pool_diploids": 500,
    "sample_diploids": 100,
    "generation_time_years": 25.0,
    "output_stride_bp": 10_000,
    "present_han_unscaled_ne": 6_300.0,
    "present_han_scaled_census_diploids_q5": 1_260,
    "present_han_scaled_census_genomes_q5": 2_520,
    "census_semantics": "complete_simulated_q5_present_han_census",
}
SURVIVAL_DECODER_SETTINGS = {
    "input_format": "trees_via_streamed_vcf",
    "pair_selector": "explicit_within_individual_pairs",
    "n_pairs": 100,
    "scaled_mutation_rate": 0.000315,
    "unscaled_mutation_rate": 1.25e-8,
    "recombination_to_mutation_ratio": 0.8,
    "thresholds_years": list(THRESHOLDS_YEARS),
    "generation_time_years": 25.0,
    "output_at_stride": 10_000,
    "output_at_hets": False,
    "cache_size": 1_000,
    "pair_block": 256,
    "exp10": "accurate",
    "backward_alignment": "fixed",
    "recent_call": "median",
    "recent_call_probability": 0.5,
    "threads": 1,
}
SURVIVAL_SIMULATION_LAW = {
    "engine": "stdpopsim_slim_forward",
    "demographic_model": "AncientEurasia_9K19",
    "sequence_scope": "direct_full_10mb",
    "mutation_introduction": (
        "single_copy_drawmutation_neanderthal_at_2400_generations"
    ),
    "conditioning": "internal_slim_extended_event_rejection",
    "branch_enumeration": False,
    "ancestry_scan": False,
    "branch_span_weighting": False,
    "cropping": False,
    "posthoc_causal_mutation_placement": False,
    "background_mutation_overlay": "stdpopsim_post_slim_neutral_msprime_overlay",
    "background_overlay_reserved_interval": (
        "every_layer_masks_0based_half_open_[5000000,5000001)"
    ),
    "reserved_focal_base_overlay_guard": (
        "scoped_msprime_sim_mutations_rate_map_mask_during_engine_simulate"
    ),
    "nonspatial_individual_location_canonicalization": (
        "post_uniform_panel_pre_analysis_exact_500x3_float64_zeros"
    ),
    "external_pool_frequency_rejection": False,
    "external_sample_detection_rejection": False,
    "external_terminal_target_af_rejection": False,
}
SURVIVAL_FOCAL_OVERLAY_PATCH = {
    "schema": SURVIVAL_FOCAL_OVERLAY_PATCH_SCHEMA,
    "scope": "only_during_one_stdpopsim_slim_engine_simulate_call",
    "patched_callable": "msprime.mutations.sim_mutations",
    "patched_callable_signature": (
        "(tree_sequence, rate=None, *, random_seed=None, model=None, "
        "start_time=None, end_time=None, discrete_genome=None, keep=None, "
        "record_provenance=True)"
    ),
    "patched_callable_source_sha256": (
        "c4b787534d05f07b09d17e47a77115234377d987618e5ee046729be4ad30423f"
    ),
    "upstream_overlay_method": ("stdpopsim.slim_engine._SLiMEngine._recap_and_rescale"),
    "upstream_overlay_method_signature": (
        "(self, ts, seed, recap_epoch, contig, slim_scaling_factor, "
        "keep_mutation_ids_as_alleles, extended_events=None)"
    ),
    "upstream_overlay_method_source_sha256": (
        "dd6cb4f6373b73638483c5053cdd5e28c42820cb69b54641e9f9a1537820cda2"
    ),
    "target_interval_0based_half_open": [5_000_000, 5_000_001],
    "detection": (
        "every intercepted msprime.RateMap with positive integrated rate on the "
        "target interval"
    ),
    "replacement": (
        "insert target breakpoints and set only the target interval rate to zero"
    ),
    "non_rate_arguments": "preserved_by_identity",
    "flank_rates_and_integrals": "exactly_preserved",
    "minimum_patched_call_count": 1,
    "expected_call_count": 2,
    "expected_patched_call_count": 1,
    "expected_call_inventory": [
        {
            "layer": "ordinary_neutral_dfe",
            "model_class": "msprime.mutations.SLiMMutationModel",
            "model_type": 0,
            "start_time": None,
            "end_time": "finite_shared_slim_tick",
            "original_positions": [0.0, 5_000_000.0, 5_000_001.0, 10_000_000.0],
            "original_rates": [1.25e-8, 0.0, 1.25e-8],
            "passed_positions": [0.0, 5_000_000.0, 5_000_001.0, 10_000_000.0],
            "passed_rates": [1.25e-8, 0.0, 1.25e-8],
            "patched": False,
        },
        {
            "layer": "recapitation",
            "model_class": "msprime.mutations.SLiMMutationModel",
            "model_type": 2,
            "start_time": "finite_shared_slim_tick",
            "end_time": None,
            "original_positions": [0.0, 10_000_000.0],
            "original_rates": [1.25e-8],
            "passed_positions": [0.0, 5_000_000.0, 5_000_001.0, 10_000_000.0],
            "passed_rates": [1.25e-8, 0.0, 1.25e-8],
            "patched": True,
        },
    ],
    "inventory_freeze_state": "frozen_after_rep051_smoke",
    "restoration": "original_callable_restored_on_success_and_exception",
}
SURVIVAL_INDIVIDUAL_LOCATION_CANONICALIZATION = {
    "schema": SURVIVAL_INDIVIDUAL_LOCATION_SCHEMA,
    "scope": "immediately_after_uniform_panel_before_analysis_or_persistence",
    "expected_individual_count": 500,
    "required_location_width_per_individual": 3,
    "precanonicalization_values": "arbitrary_float64_bits_semantically_discarded",
    "canonical_location": [0.0, 0.0, 0.0],
    "canonical_flat_location_record": {
        "dtype": "<f8",
        "shape": [1_500],
        "sha256": "ff6698a6e831ffcf47af2fed388ffc262f319e72b26cd140929d1e19b1246ad4",
    },
    "canonical_location_offset_record": {
        "dtype": "<u4",
        "shape": [501],
        "sha256": "5d348c3a0a5849e32e7c21925d7532df9762cf2fa2435d1ff297b9c0f1602a9d",
    },
    "mutation": "IndividualTable.packset_location_only",
    "mutation_api": "tskit.tables.IndividualTable.packset_location",
    "mutation_api_signature": "(self, locations)",
    "mutation_api_source_sha256": (
        "a41ee8aecaeb0dc3829403bba3df77caf9e4748ee48ed595874bab1fa2ffa5f3"
    ),
    "preserve_location_offsets_and_widths": True,
    "preserve_all_nonlocation_individual_columns": True,
    "preserve_all_other_tables_ts_metadata_and_reference": True,
}

_SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")


@dataclass(frozen=True)
class AnalysisSourceSpec:
    key: str
    label: str
    analysis_source: str


ANALYSIS_SOURCES = (
    AnalysisSourceSpec("gamma_smc", "Gamma-SMC", "gamma_smc"),
    AnalysisSourceSpec("tree_truth", "Tree truth", "tree_truth"),
)


@dataclass(frozen=True)
class FileCheck:
    scope: str
    source: str
    unit_id: str
    stage: str
    record_key: str
    path: Path
    path_identity: str
    declared_sha256: str | None
    declared_size_bytes: int | None


@dataclass
class AnalysisBundle:
    spec: AnalysisSourceSpec
    completion: dict[str, Any]
    completion_path: Path
    tables: dict[str, pd.DataFrame]
    input_paths: dict[str, Path]
    input_contract: dict[str, Any]
    audit: pd.DataFrame


@dataclass
class ExistingUnitBundle:
    units: pd.DataFrame
    endpoints: pd.DataFrame
    unit_scores: pd.DataFrame
    class_summaries: dict[str, pd.DataFrame]
    input_contract: dict[str, Any]
    audit: pd.DataFrame


@dataclass
class SurvivalBundle:
    scores: pd.DataFrame
    endpoints: pd.DataFrame
    completion: dict[str, Any]
    completion_path: Path
    input_contract: dict[str, Any]
    audit: pd.DataFrame


@dataclass
class PreparedInputs:
    unit_scores: pd.DataFrame
    selected_pointwise: pd.DataFrame
    selected_minp: pd.DataFrame
    power: pd.DataFrame
    endpoints: pd.DataFrame
    cache_audit: pd.DataFrame
    comparison: pd.DataFrame
    contract: dict[str, Any]


def _canonical_sha256(value: Any) -> str:
    payload = json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _array_record(values: np.ndarray) -> dict[str, Any]:
    array = np.ascontiguousarray(values)
    return {
        "dtype": array.dtype.str,
        "shape": list(array.shape),
        "sha256": hashlib.sha256(array.tobytes(order="C")).hexdigest(),
    }


def _binary_structure_sha256(value: Any) -> str:
    """Hash nested tskit ``asdict`` data without expanding large arrays."""

    digest = hashlib.sha256()

    def add(blob: bytes) -> None:
        digest.update(len(blob).to_bytes(8, "big"))
        digest.update(blob)

    def visit(candidate: Any) -> None:
        if isinstance(candidate, Mapping):
            add(b"mapping")
            for key in sorted(candidate, key=str):
                add(str(key).encode("utf-8"))
                visit(candidate[key])
        elif isinstance(candidate, np.ndarray):
            add(b"ndarray")
            add(candidate.dtype.str.encode("ascii"))
            add(json.dumps(list(candidate.shape)).encode("ascii"))
            add(np.ascontiguousarray(candidate).tobytes(order="C"))
        elif isinstance(candidate, bytes):
            add(b"bytes")
            add(candidate)
        elif isinstance(candidate, (list, tuple)):
            add(type(candidate).__name__.encode("ascii"))
            for item in candidate:
                visit(item)
        elif isinstance(candidate, np.generic):
            visit(candidate.item())
        elif candidate is None or isinstance(candidate, (str, int, float, bool)):
            add(type(candidate).__name__.encode("ascii"))
            add(
                json.dumps(
                    candidate, sort_keys=True, separators=(",", ":"), allow_nan=False
                ).encode("utf-8")
            )
        else:
            raise TypeError(
                "unsupported value in tskit table digest: "
                f"{type(candidate).__module__}.{type(candidate).__qualname__}"
            )

    visit(value)
    return digest.hexdigest()


def _individual_location_state(tables: tskit.TableCollection) -> dict[str, Any]:
    payload = tables.asdict()
    individual_payload = dict(payload["individuals"])
    location = np.asarray(individual_payload.pop("location"))
    offsets = np.asarray(individual_payload["location_offset"])
    widths = np.diff(offsets)
    other_payload = dict(payload)
    other_payload.pop("individuals")
    without_location = dict(payload)
    without_location["individuals"] = individual_payload
    return {
        "individual_count": int(tables.individuals.num_rows),
        "location_value_count": int(location.size),
        "unique_location_widths": sorted({int(value) for value in widths}),
        "location_all_finite": bool(np.all(np.isfinite(location))),
        "location_finite_count": int(np.count_nonzero(np.isfinite(location))),
        "location_nan_count": int(np.count_nonzero(np.isnan(location))),
        "location_positive_infinity_count": int(
            np.count_nonzero(np.isposinf(location))
        ),
        "location_negative_infinity_count": int(
            np.count_nonzero(np.isneginf(location))
        ),
        "location_all_zero": bool(np.all(location == 0.0)),
        "location_nonzero_count": int(np.count_nonzero(location)),
        "location_record": _array_record(location),
        "location_offset_record": _array_record(offsets),
        "individual_nonlocation_sha256": _binary_structure_sha256(individual_payload),
        "other_tables_ts_metadata_reference_sha256": _binary_structure_sha256(
            other_payload
        ),
        "all_except_individual_location_sha256": _binary_structure_sha256(
            without_location
        ),
    }


def _sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_text_sha256(path: str | Path) -> str:
    source_path = Path(path)
    try:
        raw = source_path.read_bytes()
    except OSError as error:
        raise ValueError(f"source is unreadable: {source_path}") from error
    if raw.startswith(b"\xef\xbb\xbf"):
        raise ValueError(f"source has a UTF-8 BOM: {source_path}")
    try:
        text = raw.decode("utf-8", errors="strict")
    except UnicodeDecodeError as error:
        raise ValueError(f"source is not strict UTF-8: {source_path}") from error
    if "\r" in text.replace("\r\n", ""):
        raise ValueError(f"source contains a bare carriage return: {source_path}")
    return hashlib.sha256(text.replace("\r\n", "\n").encode("utf-8")).hexdigest()


def _strict_sha256(value: Any, *, label: str) -> str:
    if not isinstance(value, str) or _SHA256_PATTERN.fullmatch(value) is None:
        raise ValueError(f"{label} SHA-256 is malformed")
    return value


def _strict_positive_integer(value: Any, *, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{label} must be a positive integer")
    return value


def _validate_workers(value: Any) -> int:
    workers = _strict_positive_integer(value, label="workers")
    if workers > MAX_WORKERS:
        raise ValueError(f"workers must not exceed {MAX_WORKERS}")
    return workers


def _resolve_path(repo_root: Path, value: str | Path) -> Path:
    path = Path(value)
    if not path.is_absolute():
        path = repo_root / path
    return path.resolve()


def _path_identity(path: Path, repo_root: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(repo_root.resolve()).as_posix()
    except ValueError:
        return str(resolved)


def _require_regular_file(
    path: Path, *, boundary: Path | None = None, label: str = "required input"
) -> None:
    lexical = Path(os.path.abspath(path))
    if boundary is not None:
        root = boundary.resolve()
        try:
            relative = lexical.relative_to(root)
        except ValueError as error:
            raise ValueError(f"{label} is outside its declared root: {path}") from error
        cursor = root
        for component in relative.parts:
            cursor = cursor / component
            if cursor.is_symlink():
                raise ValueError(f"{label} traverses a symlink: {path}")
    elif lexical.is_symlink():
        raise ValueError(f"{label} is a symlink: {path}")
    if not lexical.is_file():
        raise ValueError(f"{label} is absent: {path}")


def _read_json(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"JSON input is unreadable: {path}") from error
    if not isinstance(payload, dict):
        raise TypeError(f"JSON input is not an object: {path}")
    return payload


def _validate_completion_contract(
    payload: Mapping[str, Any], *, expected_schema: str | None, label: str
) -> tuple[dict[str, Any], str]:
    if payload.get("status") != "complete":
        raise ValueError(f"{label} status is not complete")
    if expected_schema is not None and payload.get("schema") != expected_schema:
        raise ValueError(f"{label} schema is incompatible")
    contract = payload.get("contract")
    if not isinstance(contract, dict):
        raise TypeError(f"{label} contract is absent")
    observed = _canonical_sha256(contract)
    if payload.get("contract_sha256") != observed:
        raise ValueError(f"{label} contract hash is incompatible")
    return contract, observed


def _require_exact_receipt_fields(
    payload: Mapping[str, Any], expected: set[str], *, label: str
) -> None:
    if set(payload) != expected:
        raise ValueError(f"{label} receipt fields differ")
    if "created_utc" in expected:
        created = payload.get("created_utc")
        if (
            not isinstance(created, str)
            or re.fullmatch(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(?:\.\d+)?Z", created)
            is None
        ):
            raise ValueError(f"{label} creation timestamp differs")
    if "elapsed_seconds" in expected:
        elapsed = payload.get("elapsed_seconds")
        if (
            isinstance(elapsed, bool)
            or not isinstance(elapsed, (int, float))
            or not np.isfinite(float(elapsed))
            or float(elapsed) < 0
        ):
            raise ValueError(f"{label} elapsed time differs")


def _resolve_declared_child(base: Path, relative: Any, *, label: str) -> Path:
    if not isinstance(relative, str) or not relative:
        raise ValueError(f"{label} path is malformed")
    pure = PurePosixPath(relative.replace("\\", "/"))
    if pure.is_absolute() or ".." in pure.parts:
        raise ValueError(f"{label} path escapes its bundle")
    path = base.joinpath(*pure.parts)
    _require_regular_file(path, boundary=base, label=label)
    return path


def _file_record(
    path: Path,
    repo_root: Path,
    *,
    include_text_hash: bool = False,
) -> dict[str, Any]:
    record = {
        "path": _path_identity(path, repo_root),
        "sha256": _sha256_file(path),
        "size_bytes": path.stat().st_size,
    }
    if include_text_hash:
        record["canonical_text_sha256"] = _canonical_text_sha256(path)
    return record


def _check_from_record(
    *,
    scope: str,
    source: str,
    unit_id: str,
    stage: str,
    record_key: str,
    base: Path,
    record: Any,
    repo_root: Path,
) -> FileCheck:
    if not isinstance(record, Mapping):
        raise TypeError(f"{stage} output record is malformed: {record_key}")
    path = _resolve_declared_child(
        base, record.get("path"), label=f"{stage} output {record_key}"
    )
    sha = _strict_sha256(record.get("sha256"), label=f"{stage} {record_key}")
    size_value = record.get("size_bytes")
    size = None
    if size_value is not None:
        size = _strict_positive_integer(size_value, label=f"{stage} {record_key} size")
    return FileCheck(
        scope=scope,
        source=source,
        unit_id=unit_id,
        stage=stage,
        record_key=record_key,
        path=path,
        path_identity=_path_identity(path, repo_root),
        declared_sha256=sha,
        declared_size_bytes=size,
    )


def _audit_file_checks(checks: Sequence[FileCheck], *, workers: int) -> pd.DataFrame:
    def inspect(check: FileCheck) -> dict[str, Any]:
        observed_size = check.path.stat().st_size
        observed_sha = _sha256_file(check.path)
        if (
            check.declared_size_bytes is not None
            and observed_size != check.declared_size_bytes
        ):
            raise ValueError(f"input size differs from receipt: {check.path_identity}")
        if check.declared_sha256 is not None and observed_sha != check.declared_sha256:
            raise ValueError(f"input hash differs from receipt: {check.path_identity}")
        return {
            "scope": check.scope,
            "source": check.source,
            "unit_id": check.unit_id,
            "stage": check.stage,
            "record_key": check.record_key,
            "path": check.path_identity,
            "declared_sha256": check.declared_sha256 or observed_sha,
            "observed_sha256": observed_sha,
            "declared_size_bytes": (
                check.declared_size_bytes
                if check.declared_size_bytes is not None
                else observed_size
            ),
            "observed_size_bytes": observed_size,
            "status": "verified",
        }

    with ThreadPoolExecutor(max_workers=workers) as executor:
        rows = list(executor.map(inspect, checks))
    return pd.DataFrame(rows).sort_values(
        ["scope", "source", "unit_id", "stage", "record_key", "path"],
        kind="stable",
    )


def _atomic_json(path: Path, value: Any) -> None:
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        temporary.write_text(
            json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_frame(path: Path, frame: pd.DataFrame) -> None:
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        frame.to_csv(temporary, sep="\t", index=False, float_format="%.12g")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_figure(figure: plt.Figure, path: Path) -> None:
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        figure.savefig(
            temporary,
            format=path.suffix.removeprefix("."),
            dpi=300 if path.suffix == ".png" else None,
            bbox_inches=None,
        )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _save_pair(figure: plt.Figure, output_dir: Path, stem: str) -> list[Path]:
    paths = []
    for suffix in ("png", "pdf"):
        path = output_dir / f"{stem}.{suffix}"
        _atomic_figure(figure, path)
        paths.append(path)
    return paths


def _expected_output_paths() -> set[str]:
    paths = set(TEXT_OUTPUTS)
    for stem in PLOT_STEMS:
        paths.add(f"{stem}.png")
        paths.add(f"{stem}.pdf")
    return paths


def _exact_numeric(frame: pd.DataFrame, column: str) -> np.ndarray:
    try:
        values = pd.to_numeric(frame[column], errors="raise").to_numpy(dtype=float)
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"numeric column is invalid: {column}") from error
    if not np.isfinite(values).all():
        raise ValueError(f"numeric column is nonfinite: {column}")
    return values


def _exact_integer(frame: pd.DataFrame, column: str) -> np.ndarray:
    values = _exact_numeric(frame, column)
    rounded = np.rint(values)
    if not np.array_equal(values, rounded):
        raise ValueError(f"integer column is fractional: {column}")
    return rounded.astype(np.int64)


def _read_tsv(path: Path, *, label: str) -> pd.DataFrame:
    try:
        return pd.read_csv(path, sep="\t")
    except Exception as error:
        raise ValueError(f"{label} is unreadable: {path}") from error


def _analysis_output_checks(
    *,
    repo_root: Path,
    analysis_dir: Path,
    completion: Mapping[str, Any],
    spec: AnalysisSourceSpec,
) -> tuple[list[FileCheck], dict[str, dict[str, Any]]]:
    outputs = completion.get("outputs")
    if not isinstance(outputs, Mapping) or not outputs:
        raise ValueError(f"{spec.key} analysis completion lacks outputs")
    checks: list[FileCheck] = []
    records: dict[str, dict[str, Any]] = {}
    declared_paths: set[str] = set()
    for key, raw_record in sorted(outputs.items()):
        if not isinstance(key, str) or not isinstance(raw_record, Mapping):
            raise TypeError(f"{spec.key} analysis output record is malformed")
        check = _check_from_record(
            scope="analysis_bundle",
            source=spec.key,
            unit_id="",
            stage="analysis_output",
            record_key=key,
            base=analysis_dir,
            record=raw_record,
            repo_root=repo_root,
        )
        relative = check.path.relative_to(analysis_dir).as_posix()
        if relative in declared_paths:
            raise ValueError(f"{spec.key} analysis declares a duplicate output path")
        declared_paths.add(relative)
        checks.append(check)
        records[key] = dict(raw_record)
    actual_paths = {
        path.relative_to(analysis_dir).as_posix()
        for path in analysis_dir.rglob("*")
        if path.is_file() and path.name != "analysis_completion.json"
    }
    if actual_paths != declared_paths:
        raise ValueError(f"{spec.key} analysis disk inventory differs from completion")
    for path in analysis_dir.rglob("*"):
        if path.is_symlink():
            raise ValueError(f"{spec.key} analysis inventory traverses a symlink")
    return checks, records


def _analysis_source_audit(
    completion: Mapping[str, Any], *, repo_root: Path, spec: AnalysisSourceSpec
) -> tuple[pd.DataFrame, dict[str, dict[str, Any]]]:
    sources = completion.get("contract", {}).get("implementation", {}).get("sources")
    if not isinstance(sources, Mapping) or not sources:
        raise ValueError(f"{spec.key} analysis source inventory is absent")
    rows = []
    normalized: dict[str, dict[str, Any]] = {}
    for key, record in sorted(sources.items()):
        if not isinstance(record, Mapping):
            raise TypeError(f"{spec.key} analysis source record is malformed")
        relative = record.get("path")
        expected = _strict_sha256(
            record.get("sha256"), label=f"{spec.key} analysis source {key}"
        )
        if not isinstance(relative, str):
            raise TypeError(f"{spec.key} analysis source path is malformed")
        path = _resolve_path(repo_root, relative)
        _require_regular_file(path, boundary=repo_root, label="analysis source")
        observed = _canonical_text_sha256(path)
        if observed != expected:
            raise ValueError(f"{spec.key} analysis source hash differs: {relative}")
        size = path.stat().st_size
        identity = _path_identity(path, repo_root)
        rows.append(
            {
                "scope": "analysis_source",
                "source": spec.key,
                "unit_id": "",
                "stage": "implementation",
                "record_key": str(key),
                "path": identity,
                "declared_sha256": expected,
                "observed_sha256": observed,
                "declared_size_bytes": size,
                "observed_size_bytes": size,
                "status": "verified",
            }
        )
        normalized[str(key)] = {"path": identity, "sha256": observed}
    return pd.DataFrame(rows), normalized


def _analysis_input_checks(
    contract: Mapping[str, Any], *, repo_root: Path, spec: AnalysisSourceSpec
) -> tuple[list[FileCheck], dict[str, Path], dict[str, Any]]:
    inputs = contract.get("inputs")
    expected_fields = {
        "class_summaries_path",
        "class_summaries_sha256",
        "empirical_path",
        "empirical_sha256",
        "pair_summaries_path",
        "pair_summaries_sha256",
        "spatial_summaries_path",
        "spatial_summaries_sha256",
    }
    if not isinstance(inputs, Mapping) or set(inputs) != expected_fields:
        raise ValueError(f"{spec.key} analysis input fields differ")
    expected_paths: dict[str, str | None]
    if spec.key == "gamma_smc":
        expected_paths = {
            "class_summaries": (
                "focused_selection_EAS_sim/results/combined_class_summaries.tsv.gz"
            ),
            "pair_summaries": (
                "focused_selection_EAS_sim/results/combined_pair_summaries.tsv.gz"
            ),
            "spatial_summaries": (
                "focused_selection_EAS_sim/results/aggregated_spatial_profiles.tsv.gz"
            ),
            "empirical": None,
        }
    elif spec.key == "tree_truth":
        expected_paths = {
            "class_summaries": (
                "focused_selection_EAS_sim/results/"
                "combined_truth_class_summaries.tsv.gz"
            ),
            "pair_summaries": None,
            "spatial_summaries": (
                "focused_selection_EAS_sim/results/"
                "aggregated_truth_spatial_profiles.tsv.gz"
            ),
            "empirical": None,
        }
    else:
        raise ValueError(f"analysis source is incompatible: {spec.key}")
    checks: list[FileCheck] = []
    paths: dict[str, Path] = {}
    records: dict[str, Any] = {}
    for key, expected_relative in expected_paths.items():
        declared_path = inputs.get(f"{key}_path")
        declared_sha = inputs.get(f"{key}_sha256")
        if expected_relative is None:
            if declared_path is not None or declared_sha is not None:
                raise ValueError(f"{spec.key} null analysis input differs: {key}")
            records[key] = None
            continue
        if declared_path != expected_relative:
            raise ValueError(f"{spec.key} analysis input path differs: {key}")
        sha = _strict_sha256(declared_sha, label=f"{spec.key} analysis input {key}")
        path = _resolve_path(repo_root, expected_relative)
        _require_regular_file(path, boundary=repo_root, label="analysis input")
        paths[key] = path
        records[key] = {
            "path": _path_identity(path, repo_root),
            "sha256": sha,
        }
        checks.append(
            FileCheck(
                "analysis_input",
                spec.key,
                "",
                "analysis_input",
                key,
                path,
                _path_identity(path, repo_root),
                sha,
                None,
            )
        )
    return checks, paths, records


def _load_analysis_bundle(
    repo_root: Path,
    analysis_dir: Path,
    spec: AnalysisSourceSpec,
    *,
    workers: int,
) -> AnalysisBundle:
    completion_path = analysis_dir / "analysis_completion.json"
    _require_regular_file(
        completion_path, boundary=analysis_dir, label=f"{spec.key} completion"
    )
    completion = _read_json(completion_path)
    contract, contract_sha = _validate_completion_contract(
        completion,
        expected_schema=ANALYSIS_SCHEMA_VERSION,
        label=f"{spec.key} analysis",
    )
    input_checks, input_paths, input_records = _analysis_input_checks(
        contract, repo_root=repo_root, spec=spec
    )
    checks, output_records = _analysis_output_checks(
        repo_root=repo_root,
        analysis_dir=analysis_dir,
        completion=completion,
        spec=spec,
    )
    completion_check = FileCheck(
        scope="analysis_bundle",
        source=spec.key,
        unit_id="",
        stage="analysis_completion",
        record_key="analysis_completion",
        path=completion_path,
        path_identity=_path_identity(completion_path, repo_root),
        declared_sha256=None,
        declared_size_bytes=None,
    )
    audit = _audit_file_checks(
        [completion_check, *checks, *input_checks], workers=workers
    )
    source_audit, source_records = _analysis_source_audit(
        completion, repo_root=repo_root, spec=spec
    )
    audit = pd.concat([audit, source_audit], ignore_index=True).sort_values(
        ["scope", "source", "unit_id", "stage", "record_key", "path"],
        kind="stable",
    )
    tables: dict[str, pd.DataFrame] = {}
    table_contract: dict[str, Any] = {}
    for key, filename in TABLE_FILENAMES.items():
        record = output_records.get(key)
        if not isinstance(record, Mapping) or record.get("path") != filename:
            raise ValueError(f"{spec.key} analysis lacks required table: {key}")
        path = analysis_dir / filename
        frame = _read_tsv(path, label=f"{spec.key} {key}")
        declared_rows = record.get("rows")
        declared_columns = record.get("columns")
        if declared_rows != len(frame) or declared_columns != list(frame.columns):
            raise ValueError(f"{spec.key} {key} schema differs from completion")
        tables[key] = frame
        table_contract[key] = {
            "path": _path_identity(path, repo_root),
            "sha256": _sha256_file(path),
            "size_bytes": path.stat().st_size,
            "rows": len(frame),
            "columns": list(frame.columns),
        }
    completion_record = _file_record(completion_path, repo_root)
    return AnalysisBundle(
        spec=spec,
        completion=completion,
        completion_path=completion_path,
        tables=tables,
        input_paths=input_paths,
        audit=audit,
        input_contract={
            "completion": completion_record,
            "contract_sha256": contract_sha,
            "declared_output_count": len(output_records),
            "declared_outputs_sha256": _canonical_sha256(output_records),
            "required_tables": table_contract,
            "upstream_sources": source_records,
            "analysis_inputs": input_records,
            "analysis_inputs_sha256": _canonical_sha256(input_records),
        },
    )


def _validate_study_design(design: Mapping[str, Any]) -> None:
    if design.get("schema") != "gamma-smc.focused-selection-campaign/v1":
        raise ValueError("focused campaign study design schema is incompatible")
    grid = design.get("scientific_grid", {})
    demographies = grid.get("demographies")
    if not isinstance(demographies, list):
        raise TypeError("focused campaign demography grid is absent")
    matches = [
        item for item in demographies if item.get("demography_id") == DEMOGRAPHY_ID
    ]
    if len(matches) != 1 or matches[0].get("source_model") != SOURCE_MODEL:
        raise ValueError("focused campaign does not bind AncientEurasia_9K19")
    sequence = design.get("sequence", {})
    time = design.get("time", {})
    if sequence.get("focal_position_0based") != FOCAL_POSITION_BP:
        raise ValueError("focused campaign focal position is incompatible")
    if tuple(time.get("thresholds_years", ())) != THRESHOLDS_YEARS:
        raise ValueError("focused campaign threshold grid is incompatible")
    if not np.isclose(time.get("generation_time_years"), GENERATION_TIME_YEARS):
        raise ValueError("focused campaign generation time is incompatible")


def _select_existing_units(frame: pd.DataFrame) -> pd.DataFrame:
    required = {
        "unit_id",
        "demography_id",
        "source_model",
        "simulation_class",
        "selection_coefficient",
        "target_allele_frequency",
        "replicate_index",
        "seed",
        "sample_diploids",
        "exact_sample_alt_count",
    }
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"execution inventory lacks columns: {', '.join(missing)}")
    coefficient = _exact_numeric(frame, "selection_coefficient")
    target = _exact_numeric(frame, "target_allele_frequency")
    selected = frame[
        (frame["demography_id"].astype(str) == DEMOGRAPHY_ID)
        & (frame["source_model"].astype(str) == SOURCE_MODEL)
        & np.isclose(target, TARGET_AF)
        & (
            (
                (frame["simulation_class"].astype(str) == "neutral")
                & np.isclose(coefficient, 0.0)
            )
            | (
                (frame["simulation_class"].astype(str) == "selected")
                & np.isclose(coefficient, SELECTION_COEFFICIENT)
            )
        )
    ].copy()
    counts = selected.groupby("simulation_class").size().to_dict()
    if counts != {
        "neutral": MATCHED_AF_NEUTRAL_REPLICATES,
        "selected": SELECTED_REPLICATES,
    }:
        raise ValueError("focused execution inventory is not the exact 10+100 cell")
    expected_ids = {
        *(
            f"{DEMOGRAPHY_ID}__af20__selected_s0p010__rep{index:03d}"
            for index in range(1, SELECTED_REPLICATES + 1)
        ),
        *(
            f"{DEMOGRAPHY_ID}__af20__neutral__rep{index:03d}"
            for index in range(1, MATCHED_AF_NEUTRAL_REPLICATES + 1)
        ),
    }
    if set(selected["unit_id"].astype(str)) != expected_ids:
        raise ValueError("focused execution unit IDs are incompatible")
    if set(_exact_integer(selected, "sample_diploids")) != {SAMPLE_DIPLOIDS}:
        raise ValueError("focused execution sample size is incompatible")
    if set(_exact_integer(selected, "exact_sample_alt_count")) != {40}:
        raise ValueError("focused execution exact sample AF is incompatible")
    if selected["unit_id"].duplicated().any():
        raise ValueError("focused execution inventory contains duplicate unit IDs")
    exact_common = {
        "demography_id": DEMOGRAPHY_ID,
        "demography_kind": "introgression",
        "population": "Han",
        "source_model": SOURCE_MODEL,
        "selection_origin": "archaic_specific_introgressed_standing_variation",
        "target_allele_frequency": 0.2,
        "population_af_lower": 0.175,
        "population_af_upper": 0.225,
        "exact_sample_alt_count": 40,
        "sample_diploids": 100,
        "sequence_length_bp": 10_000_000,
        "focal_position_bp": 5_000_000,
    }
    for record in selected.to_dict(orient="records"):
        if any(record.get(key) != value for key, value in exact_common.items()):
            raise ValueError(f"focused execution semantics differ: {record['unit_id']}")
        simulation_class = str(record["simulation_class"])
        expected_selection = 0.01 if simulation_class == "selected" else 0.0
        suffix = re.search(r"__rep(\d{3})$", str(record["unit_id"]))
        seed = record.get("seed")
        replicate = record.get("replicate_index")
        if (
            simulation_class not in {"selected", "neutral"}
            or record.get("selection_coefficient") != expected_selection
            or suffix is None
            or isinstance(replicate, (bool, np.bool_))
            or not isinstance(replicate, (int, np.integer))
            or replicate != int(suffix.group(1))
            or isinstance(seed, (bool, np.bool_))
            or not isinstance(seed, (int, np.integer))
            or seed < 1
        ):
            raise ValueError(f"focused execution unit differs: {record['unit_id']}")
    return selected.sort_values("unit_id", kind="stable").reset_index(drop=True)


def _validate_focused_simulation_contract(
    contract: Mapping[str, Any], *, unit_id: str, simulation_class: str
) -> None:
    if contract.get("software") != FOCUSED_SIMULATION_SOFTWARE:
        raise ValueError(f"focused simulation software differs: {unit_id}")
    parameters = contract.get("parameters")
    if (
        not isinstance(parameters, Mapping)
        or set(parameters) != FOCUSED_SIMULATION_PARAMETER_FIELDS
    ):
        raise ValueError(f"focused simulation parameter fields differ: {unit_id}")
    common = {
        "candidate_pool_diploids": 500,
        "candidate_pool_frequency_gate": {
            "interval_source": (
                "unit.population_af_lower_to_population_af_upper_inclusive"
            ),
            "order": "before_exact_sample_panel",
            "schema": "gamma-smc.candidate-pool-frequency-gate/v1",
            "scope": "selected_and_neutral",
        },
        "introgression_mutation_age_generations": 2_400.0,
        "introgression_source_presence_epsilon": 1e-9,
        "neutral_ancestry_batch_size": 4,
        "neutral_ancestry_model": "DTWF200ThenSmcPrimeApproxCoalescent",
        "neutral_ancestry_model_keywords": ["dtwf", "smc_prime"],
        "neutral_candidate_index": {
            "introgression": "pulse_migrations_then_branch_crossing_2400_generations",
            "no_introgression": "marginal_tree_at_fixed_5Mb_site",
        },
        "neutral_eas_ascertainment_length_bp": 10_000_000,
        "neutral_han_ascertainment_length_bp": 20_000_000,
        "neutral_pool_diploids": 500,
        "neutral_recent_dtwf_duration_generations": 200.0,
        "selected_pool_diploids": 500,
        "slim_burn_in": 0.1,
        "slim_scaling_factor": 5.0,
    }
    if any(parameters.get(key) != value for key, value in common.items()):
        raise ValueError(f"focused simulation parameters differ: {unit_id}")
    focal = parameters.get("focal_semantic_validation")
    if (
        not isinstance(focal, Mapping)
        or set(focal)
        != {
            "exact_sample_alt_count",
            "identity",
            "minimum_heterozygous_diploids",
            "minimum_hom_alt_diploids",
            "minimum_hom_ref_diploids",
            "sample_diploids",
            "schema",
            "simulation_class",
            "tree_manifest_exact_match",
        }
        or focal.get("exact_sample_alt_count") != 40
        or focal.get("minimum_heterozygous_diploids") != 1
        or focal.get("minimum_hom_alt_diploids") != 2
        or focal.get("minimum_hom_ref_diploids") != 2
        or focal.get("sample_diploids") != 100
        or focal.get("schema") != "gamma-smc.focused-focal-semantics/v1"
        or focal.get("simulation_class") != simulation_class
        or focal.get("tree_manifest_exact_match") is not True
    ):
        raise ValueError(f"focused focal contract differs: {unit_id}")
    identity = focal.get("identity")
    if not isinstance(identity, Mapping):
        raise TypeError(f"focused focal identity contract differs: {unit_id}")
    if simulation_class == "neutral":
        identity_fixed = {
            "schema": "gamma-smc.neutral-focal-identity/v1",
            "focal_position_0based": FOCAL_POSITION_BP,
            "expected_origin_time_generations": 2_400.0,
            "expected_source_population": "Neanderthal",
            "ancestral_state": "A",
            "derived_state": "G",
            "single_nonrecurrent_mutation": True,
            "archaic_specific_no_ils_by_construction": True,
            "human_neanderthal_split_generations": 27_840.0,
        }
        pulse = identity.get("pulse_route")
        if (
            set(identity) != {*identity_fixed, "pulse_route"}
            or any(identity.get(key) != value for key, value in identity_fixed.items())
            or pulse
            != {
                "destination_population": "Neanderthal",
                "minimum_descendant_lineages": 1,
                "source_population": "Loschbour",
                "time_generations": 2_272.0,
            }
            or parameters.get("han_selection_calibration") is not None
        ):
            raise ValueError(f"focused neutral focal contract differs: {unit_id}")
    else:
        identity_expected = {
            "expected_mutation_type": 1,
            "expected_selection_coeff": 0.0,
            "focal_position_0based": FOCAL_POSITION_BP,
            "realized_origin_time_generations": 2_400.0,
            "requested_origin_time_generations": 2_400.0,
            "schema": "gamma-smc.selected-focal-identity/v1",
            "single_site_id": FOCAL_SITE_ID,
            "slim_scaling_factor": 5.0,
            "slim_time_rule": "cycle - realized_origin_time / Q",
            "source_population": "Neanderthal",
            "source_subpopulation": 7,
        }
        calibration = parameters.get("han_selection_calibration")
        if identity != identity_expected or not isinstance(calibration, Mapping):
            raise ValueError(f"focused selected focal contract differs: {unit_id}")
        calibration_fixed = {
            "failed_v2_inferential_use": False,
            "gamma_smc_statistics_used": False,
            "manifest_path": (
                "focused_selection_EAS_sim/calibration/han_direct_production_v3/"
                "han_direct_v3_manifest.json"
            ),
            "manifest_sha256": (
                "7ba2be1875ea804c1cdbec4ae3b9fa5dffb8307799d6af0c810ccc9df5c40e79"
            ),
            "realized_selection_end_generations_ago": 1_780.0,
            "schema": "gamma-smc.han-direct-production-v3-authorization/v1",
            "seed_plan_path": (
                "focused_selection_EAS_sim/calibration/han_direct_production_v3/"
                "han_direct_v3_attempt_seed_plan.tsv"
            ),
            "seed_plan_sha256": (
                "04017a12d2a39709777be7fc9767fcf6834bf714519296bf3a7f0878e536f376"
            ),
            "status": "planned_direct_production_unit",
        }
        if (
            set(calibration)
            != {*calibration_fixed, "attempt_schedule_sha256", "unit_id"}
            or any(
                calibration.get(key) != value
                for key, value in calibration_fixed.items()
            )
            or calibration.get("unit_id") != unit_id
        ):
            raise ValueError(f"focused selected calibration differs: {unit_id}")
        _strict_sha256(
            calibration.get("attempt_schedule_sha256"),
            label=f"{unit_id} selected attempt schedule",
        )
    expected_implementation = {
        "eas_resource_sha256": (
            "a5bdfd843629d48050cd56a84da4d70e356d0339017d3c1f67be646d727db928"
        ),
        "slim": (
            {
                "path": EXPECTED_SLIM_BINARY_IDENTITY,
                "sha256": EXPECTED_SLIM_BINARY_SHA256,
            }
            if simulation_class == "selected"
            else None
        ),
        "sources": (
            FOCUSED_SELECTED_SIMULATION_SOURCE_HASHES
            if simulation_class == "selected"
            else FOCUSED_NEUTRAL_SIMULATION_SOURCE_HASHES
        ),
    }
    if contract.get("implementation") != expected_implementation:
        raise ValueError(f"focused simulation implementation differs: {unit_id}")


def _focused_genotype_counts(value: Any, *, label: str) -> dict[str, int]:
    if not isinstance(value, Mapping) or set(value) != {
        "hom_ref",
        "heterozygous",
        "hom_alt",
    }:
        raise ValueError(f"{label} genotype fields differ")
    if any(
        isinstance(value[key], bool) or not isinstance(value[key], int) for key in value
    ):
        raise ValueError(f"{label} genotype counts differ")
    counts = {key: value[key] for key in value}
    if any(count < 0 for count in counts.values()):
        raise ValueError(f"{label} genotype counts differ")
    return counts


def _validate_focused_simulation_provenance(
    payload: Mapping[str, Any],
    *,
    unit_id: str,
    simulation_class: str,
    repo_root: Path,
) -> None:
    sample_counts = _focused_genotype_counts(
        payload.get("genotype_counts"), label=f"{unit_id} sample"
    )
    if (
        sum(sample_counts.values()) != 100
        or 2 * sample_counts["hom_alt"] + sample_counts["heterozygous"] != 40
        or payload.get("sample_alt_count") != 40
        or payload.get("sample_af") != 0.2
    ):
        raise ValueError(f"focused sample endpoint differs: {unit_id}")
    provenance = payload.get("provenance")
    expected_fields = (
        FOCUSED_SELECTED_PROVENANCE_FIELDS
        if simulation_class == "selected"
        else FOCUSED_NEUTRAL_PROVENANCE_FIELDS
    )
    if not isinstance(provenance, Mapping) or set(provenance) != expected_fields:
        raise ValueError(f"focused simulation provenance fields differ: {unit_id}")
    gate = provenance.get("candidate_pool_frequency_gate")
    gate_fields = {
        "candidate_pool_diploids",
        "candidate_pool_genotype_counts",
        "candidate_pool_haplotypes",
        "gate_order",
        "maximum_alt_count_inclusive",
        "minimum_alt_count_inclusive",
        "observed_af",
        "observed_alt_count",
        "passed",
        "population_af_lower_inclusive",
        "population_af_upper_inclusive",
        "schema",
        "status",
    }
    gate_fixed = {
        "candidate_pool_diploids": 500,
        "candidate_pool_haplotypes": 1_000,
        "gate_order": "before_exact_sample_panel",
        "maximum_alt_count_inclusive": 225,
        "minimum_alt_count_inclusive": 175,
        "passed": True,
        "population_af_lower_inclusive": 0.175,
        "population_af_upper_inclusive": 0.225,
        "schema": "gamma-smc.candidate-pool-frequency-gate/v1",
        "status": "pass",
    }
    if (
        not isinstance(gate, Mapping)
        or set(gate) != gate_fields
        or any(gate.get(key) != value for key, value in gate_fixed.items())
    ):
        raise ValueError(f"focused candidate-pool gate differs: {unit_id}")
    pool_counts = _focused_genotype_counts(
        gate.get("candidate_pool_genotype_counts"), label=f"{unit_id} pool"
    )
    pool_alt_raw = gate.get("observed_alt_count")
    if isinstance(pool_alt_raw, bool) or not isinstance(pool_alt_raw, int):
        raise TypeError(f"focused candidate-pool count differs: {unit_id}")
    pool_alt = pool_alt_raw
    if (
        sum(pool_counts.values()) != 500
        or 2 * pool_counts["hom_alt"] + pool_counts["heterozygous"] != pool_alt
        or not 175 <= pool_alt <= 225
        or gate.get("observed_af") != pool_alt / 1_000
    ):
        raise ValueError(f"focused candidate-pool counts differ: {unit_id}")
    focal = provenance.get("focal_semantic_validation")
    focal_fields = {
        "genotype_counts",
        "identity",
        "sample_af",
        "sample_alt_count",
        "sample_diploids",
        "schema",
        "status",
        "tree_manifest_exact_match",
    }
    if (
        not isinstance(focal, Mapping)
        or set(focal) != focal_fields
        or focal.get("genotype_counts") != sample_counts
        or focal.get("sample_af") != 0.2
        or focal.get("sample_alt_count") != 40
        or focal.get("sample_diploids") != 100
        or focal.get("schema") != "gamma-smc.focused-focal-semantics/v1"
        or focal.get("status") != "valid"
        or focal.get("tree_manifest_exact_match") is not True
    ):
        raise ValueError(f"focused focal provenance differs: {unit_id}")
    sample_panel = provenance.get("sample_panel")
    if (
        not isinstance(sample_panel, Mapping)
        or set(sample_panel)
        != {
            "candidate_pool_diploids",
            "candidate_pool_genotype_counts",
            "chosen_sample_genotype_counts",
            "eligible_genotype_count_triples",
            "panel_seed",
            "selection_method",
        }
        or sample_panel.get("candidate_pool_diploids") != 500
        or sample_panel.get("candidate_pool_genotype_counts") != pool_counts
        or sample_panel.get("chosen_sample_genotype_counts") != sample_counts
        or sample_panel.get("selection_method")
        != "uniform_over_all_size_N_exact_alt_count_subsets_passing_genotype_qc"
        or isinstance(sample_panel.get("panel_seed"), bool)
        or not isinstance(sample_panel.get("panel_seed"), int)
        or sample_panel.get("panel_seed") < 1
        or isinstance(sample_panel.get("eligible_genotype_count_triples"), bool)
        or not isinstance(sample_panel.get("eligible_genotype_count_triples"), int)
        or sample_panel.get("eligible_genotype_count_triples") < 1
    ):
        raise ValueError(f"focused sample-panel provenance differs: {unit_id}")
    if provenance.get("matched_selected_neutral_ascertainment") != (
        "same_candidate_pool_size_and_AF_band_then_same_exact_k_panel_sampler"
    ):
        raise ValueError(f"focused ascertainment provenance differs: {unit_id}")
    identity_key = (
        "selected_focal_identity"
        if simulation_class == "selected"
        else "neutral_focal_identity"
    )
    if provenance.get(identity_key) != focal.get("identity"):
        raise ValueError(f"focused identity provenance differs: {unit_id}")
    identity = focal["identity"]
    if simulation_class == "neutral":
        fixed = {
            "ancestry_batch_size": 4,
            "archaic_specific_no_ils": True,
            "ascertainment_length_bp": 20_000_000,
            "candidate_pool_diploids": 500,
            "candidate_weight": "genomic_span_bp_at_fixed_2400_generation_origin",
            "coalescent_approximation": "SMC-prime after 200 generations",
            "engine": "msprime",
            "entry_route": (
                "Loschbour_to_Neanderthal_backward_mass_migration_at_2272_generations"
            ),
            "model": "DTWF200ThenSmcPrimeApproxCoalescent",
            "model_keywords": ["dtwf", "smc_prime"],
            "msprime_version": "1.4.2",
            "mutation_origin_population": "Neanderthal",
            "mutation_time_generations_ago": 2_400.0,
            "neutral_origin_matching": (
                "archaic_population_branch_with_one_or_more_pulse_lineages; "
                "nonzero standing variation at pulse"
            ),
            "recent_dtwf_duration_generations": 200.0,
            "window_length_bp": 10_000_000,
        }
        if any(provenance.get(key) != value for key, value in fixed.items()):
            raise ValueError(f"focused neutral provenance differs: {unit_id}")
        pulse = (
            identity.get("pulse_route_evidence")
            if isinstance(identity, Mapping)
            else None
        )
        if (
            not isinstance(identity, Mapping)
            or set(identity)
            != {
                "ancestral_state",
                "archaic_specific_no_ils_by_construction",
                "derived_state",
                "focal_position_0based",
                "mutation_metadata_empty",
                "observed_source_population",
                "origin_time_generations",
                "pulse_route_evidence",
                "schema",
                "status",
            }
            or identity.get("ancestral_state") != "A"
            or identity.get("derived_state") != "G"
            or identity.get("archaic_specific_no_ils_by_construction") is not True
            or identity.get("focal_position_0based") != FOCAL_POSITION_BP
            or identity.get("mutation_metadata_empty") is not True
            or identity.get("observed_source_population") != "Neanderthal"
            or identity.get("origin_time_generations") != 2_400.0
            or identity.get("schema") != "gamma-smc.neutral-focal-identity/v1"
            or identity.get("status") != "valid"
            or not isinstance(pulse, Mapping)
            or set(pulse)
            != {
                "descendant_lineage_nodes",
                "destination_population",
                "n_descendant_lineages",
                "source_population",
                "time_generations",
            }
            or pulse.get("destination_population") != "Neanderthal"
            or pulse.get("source_population") != "Loschbour"
            or pulse.get("time_generations") != 2_272.0
            or isinstance(pulse.get("n_descendant_lineages"), bool)
            or not isinstance(pulse.get("n_descendant_lineages"), int)
            or pulse.get("n_descendant_lineages") < 1
        ):
            raise ValueError(f"focused neutral identity differs: {unit_id}")
        if provenance.get("panel_seed") != sample_panel.get("panel_seed"):
            raise ValueError(f"focused neutral panel seed differs: {unit_id}")
        strict_positive = (
            "ancestry_seed",
            "candidate_seed",
            "mutation_seed",
            "panel_seed",
            "n_eligible_branch_rectangles_in_batch",
            "n_introgressed_descendant_lineages",
        )
        if any(
            isinstance(provenance.get(key), bool)
            or not isinstance(provenance.get(key), int)
            or provenance.get(key) < 1
            for key in strict_positive
        ):
            raise ValueError(f"focused neutral integer provenance differs: {unit_id}")
        accepted_attempt = provenance.get("accepted_attempt_zero_based")
        pre_crop = provenance.get("pre_crop_focal_position")
        lineage_nodes = pulse.get("descendant_lineage_nodes")
        if (
            isinstance(provenance.get("accepted_batch_zero_based"), bool)
            or provenance.get("accepted_batch_zero_based") != 0
            or isinstance(accepted_attempt, bool)
            or not isinstance(accepted_attempt, int)
            or not 0 <= accepted_attempt <= 3
            or provenance.get("candidate_index")
            != (
                "Loschbour_to_Neanderthal_pulse_migrations_then_ascend_to_"
                "branch_crossing_2400_generations"
            )
            or isinstance(pre_crop, bool)
            or not isinstance(pre_crop, (int, float))
            or not np.isfinite(float(pre_crop))
            or not 5_000_000 <= float(pre_crop) < 15_000_000
            or not isinstance(lineage_nodes, list)
            or any(
                isinstance(node, bool) or not isinstance(node, int) or node < 0
                for node in lineage_nodes
            )
            or len(set(lineage_nodes)) != len(lineage_nodes)
            or len(lineage_nodes) != pulse.get("n_descendant_lineages")
            or len(lineage_nodes)
            != provenance.get("n_introgressed_descendant_lineages")
        ):
            raise ValueError(f"focused neutral route provenance differs: {unit_id}")
    else:
        fixed = {
            "engine": "stdpopsim_slim",
            "cumulative_timeout_seconds": 90_000.0,
            "draw_timeout_seconds": 300.0,
            "pyslim_version": "1.1.1",
            "slim_burn_in": 0.1,
            "slim_scaling_factor": 5.0,
            "slim_sha256": EXPECTED_SLIM_BINARY_SHA256,
            "stdpopsim_version": "0.3.0",
        }
        endpoint = provenance.get("selection_endpoint")
        sweep = provenance.get("sweep")
        accepted_draw = provenance.get("accepted_draw_zero_based")
        accepted_seed = provenance.get("accepted_seed")
        if (
            any(provenance.get(key) != value for key, value in fixed.items())
            or not _paths_equivalent(
                str(provenance.get("slim_path")),
                repo_root / EXPECTED_SLIM_BINARY_IDENTITY,
            )
            or isinstance(accepted_draw, bool)
            or not isinstance(accepted_draw, int)
            or not 0 <= accepted_draw < 300
            or isinstance(accepted_seed, bool)
            or not isinstance(accepted_seed, int)
            or accepted_seed < 1
            or not isinstance(endpoint, Mapping)
            or _canonical_sha256(endpoint)
            != "97b5143f8b452abbf144bd3d30fd76e40030019005f831abba082fcb1c0bf17c"
            or endpoint.get("binding")
            != {
                "selection_coefficient": 0.01,
                "selection_end_generations_ago": 1_780.0,
                "target_allele_frequency": 0.2,
            }
            or endpoint.get("cessation", {}).get("selection_end_generations_ago")
            != 1_780.0
            or not isinstance(sweep, Mapping)
            or _canonical_sha256(sweep)
            != "f4072a7185d2a145102f7f1eae5ceadb5fd773e6a066a1d43c40ff1a27881240"
            or sweep.get("schema") != "gamma-smc.introgression-sweep-spec/v1"
            or sweep.get("selection", {}).get("selection_coefficient") != 0.01
            or sweep.get("selection", {}).get("dominance_coefficient") != 0.5
            or sweep.get("mutation_origin", {}).get("age_generations") != 2_400.0
            or sweep.get("introgression", {}).get("pulse_generations_ago") != 2_272.0
        ):
            raise ValueError(f"focused selected provenance differs: {unit_id}")
        nucleotide = {"A": 0, "C": 1, "G": 2, "T": 3}
        if (
            not isinstance(identity, Mapping)
            or set(identity)
            != {
                "expected_mutation_type",
                "expected_origin_time_generations",
                "expected_source_population",
                "expected_source_subpopulation",
                "focal_position_0based",
                "observed_ancestral_state",
                "observed_derived_state",
                "observed_focal_site_count",
                "observed_metadata_entry_count",
                "observed_mutation_type",
                "observed_mutation_types",
                "observed_nucleotide",
                "observed_origin_time_generations",
                "observed_selection_coeff",
                "observed_slim_cycle",
                "observed_slim_scaling_factor",
                "observed_slim_time",
                "observed_source_subpopulation",
                "observed_tskit_mutation_count",
                "schema",
                "status",
            }
            or identity.get("expected_mutation_type") != 1
            or identity.get("expected_origin_time_generations") != 2_400.0
            or identity.get("expected_source_population") != "Neanderthal"
            or identity.get("expected_source_subpopulation") != 7
            or identity.get("focal_position_0based") != FOCAL_POSITION_BP
            or identity.get("observed_ancestral_state") not in nucleotide
            or identity.get("observed_derived_state") not in nucleotide
            or identity.get("observed_ancestral_state")
            == identity.get("observed_derived_state")
            or identity.get("observed_nucleotide")
            != nucleotide[identity.get("observed_derived_state")]
            or identity.get("observed_focal_site_count") != 1
            or identity.get("observed_metadata_entry_count") != 1
            or identity.get("observed_mutation_type") != 1
            or identity.get("observed_mutation_types") != [1]
            or identity.get("observed_origin_time_generations") != 2_400.0
            or identity.get("observed_selection_coeff") != 0.0
            or identity.get("observed_slim_cycle") != 5_933
            or identity.get("observed_slim_scaling_factor") != 5.0
            or identity.get("observed_slim_time") != 5_453
            or identity.get("observed_source_subpopulation") != 7
            or identity.get("observed_tskit_mutation_count") != 1
            or identity.get("schema") != "gamma-smc.selected-focal-identity/v1"
            or identity.get("status") != "valid"
        ):
            raise ValueError(f"focused selected identity differs: {unit_id}")
    search_budget = payload.get("search_budget")
    expected_budget = {
        "neutral_max_attempts": 100,
        "selected_cumulative_timeout_seconds": (
            90_000.0 if simulation_class == "selected" else 2_700.0
        ),
        "selected_draw_timeout_seconds": 300.0,
        "selected_max_draws": 300 if simulation_class == "selected" else 50,
    }
    attempts = payload.get("attempts_completed")
    if (
        search_budget != expected_budget
        or isinstance(attempts, bool)
        or not isinstance(attempts, int)
        or attempts < 1
        or (
            simulation_class == "selected"
            and attempts != provenance.get("accepted_draw_zero_based") + 1
        )
        or (simulation_class == "neutral" and attempts != 4)
    ):
        raise ValueError(f"focused simulation search budget differs: {unit_id}")


def _validate_focused_manifest(
    path: Path, pairs_path: Path, payload: Mapping[str, Any], *, unit_id: str
) -> None:
    manifest = _read_tsv(path, label=f"{unit_id} sample manifest")
    if (
        tuple(manifest.columns) != SURVIVAL_SAMPLE_MANIFEST_COLUMNS
        or len(manifest) != 100
    ):
        raise ValueError(f"focused sample manifest schema differs: {unit_id}")
    allele_count = _exact_integer(manifest, "focal_selected_allele_count")
    if np.any((allele_count < 0) | (allele_count > 2)):
        raise ValueError(f"focused sample manifest allele counts differ: {unit_id}")
    expected_classes = np.asarray(["hom_ref", "heterozygous", "hom_alt"], dtype=object)[
        allele_count
    ]
    diploid_index = _exact_integer(manifest, "vcf_diploid_index")
    haplotype_0 = _exact_integer(manifest, "gamma_smc_haplotype_0")
    haplotype_1 = _exact_integer(manifest, "gamma_smc_haplotype_1")
    individual_ids = _exact_integer(manifest, "tree_sequence_individual_id")
    sample_node_0 = _exact_integer(manifest, "sample_node_0")
    sample_node_1 = _exact_integer(manifest, "sample_node_1")
    if (
        not np.array_equal(diploid_index, np.arange(100))
        or not np.array_equal(haplotype_0, 2 * diploid_index)
        or not np.array_equal(haplotype_1, 2 * diploid_index + 1)
        or np.any(individual_ids < 0)
        or len(set(individual_ids)) != 100
        or np.any(sample_node_0 < 0)
        or np.any(sample_node_1 < 0)
        or len(set(np.r_[sample_node_0, sample_node_1])) != 200
        or not np.all(
            pd.to_numeric(manifest["focal_position_0based"]).to_numpy(dtype=float)
            == FOCAL_POSITION_BP
        )
        or set(_exact_integer(manifest, "focal_selected_allele_index")) != {1}
        or not np.array_equal(manifest["genotype_class"].astype(str), expected_classes)
        or int(allele_count.sum()) != 40
        or manifest["genotype_class"].value_counts().to_dict()
        != payload.get("genotype_counts")
    ):
        raise ValueError(f"focused sample manifest genotypes differ: {unit_id}")
    lines = pairs_path.read_text(encoding="utf-8").splitlines()
    expected_header = [
        "# Gamma-SMC explicit haplotype pairs",
        "# genotype_class\toverall",
        "# n_pairs\t100",
    ]
    if lines[:3] != expected_header or len(lines) != 103:
        raise ValueError(f"focused overall pair manifest differs: {unit_id}")
    observed_pairs = [tuple(map(int, line.split("\t"))) for line in lines[3:]]
    expected_pairs = list(
        zip(
            haplotype_0,
            haplotype_1,
            strict=True,
        )
    )
    if observed_pairs != expected_pairs:
        raise ValueError(f"focused overall pairs differ from manifest: {unit_id}")


def _unit_receipt_checks(
    repo_root: Path,
    campaign_root: Path,
    unit: Mapping[str, Any],
) -> tuple[list[FileCheck], dict[str, Any], dict[str, Any]]:
    unit_id = str(unit["unit_id"])
    unit_dir = campaign_root / "work" / unit_id
    sim_path = unit_dir / "simulation_complete.json"
    decode_dir = unit_dir / "decoded"
    decode_path = decode_dir / "completion.json"
    _require_regular_file(
        sim_path, boundary=campaign_root, label="simulation completion"
    )
    _require_regular_file(
        decode_path, boundary=campaign_root, label="decode completion"
    )
    sim = _read_json(sim_path)
    decode = _read_json(decode_path)
    _require_exact_receipt_fields(
        sim,
        set(FOCUSED_SIMULATION_RECEIPT_FIELDS),
        label=f"{unit_id} focused simulation",
    )
    _require_exact_receipt_fields(
        decode,
        set(FOCUSED_DECODE_RECEIPT_FIELDS),
        label=f"{unit_id} focused decode",
    )
    sim_contract, sim_contract_sha = _validate_completion_contract(
        sim,
        expected_schema=FOCUSED_SIMULATION_SCHEMA_VERSION,
        label=f"{unit_id} simulation",
    )
    decode_contract, decode_contract_sha = _validate_completion_contract(
        decode,
        expected_schema=FOCUSED_DECODE_SCHEMA_VERSION,
        label=f"{unit_id} decode",
    )
    if (
        set(sim_contract) != FOCUSED_SIMULATION_CONTRACT_FIELDS
        or sim_contract.get("schema") != FOCUSED_SIMULATION_SCHEMA_VERSION
        or set(decode_contract) != FOCUSED_DECODE_CONTRACT_FIELDS
    ):
        raise ValueError(f"focused unit contract fields differ: {unit_id}")
    simulation_contract_path = unit_dir / "simulation_contract.json"
    _require_regular_file(
        simulation_contract_path,
        boundary=unit_dir,
        label="focused simulation standalone contract",
    )
    if _read_json(simulation_contract_path) != sim_contract:
        raise ValueError(f"focused standalone simulation contract differs: {unit_id}")
    sim_unit = sim_contract.get("unit")
    if not isinstance(sim_unit, Mapping) or sim_unit.get("unit_id") != unit_id:
        raise ValueError(f"simulation unit binding differs: {unit_id}")
    if set(sim_unit) != set(unit):
        raise ValueError(f"simulation unit fields differ: {unit_id}")
    for key in unit:
        if sim_unit.get(key) != unit[key]:
            raise ValueError(f"simulation unit field differs for {unit_id}: {key}")
    simulation_class = str(unit["simulation_class"])
    _validate_focused_simulation_contract(
        sim_contract, unit_id=unit_id, simulation_class=simulation_class
    )
    _validate_focused_simulation_provenance(
        sim,
        unit_id=unit_id,
        simulation_class=simulation_class,
        repo_root=repo_root,
    )
    if simulation_class == "selected":
        calibration = sim_contract["parameters"]["han_selection_calibration"]
        authorization = sim["provenance"].get("han_direct_production_v3")
        authorization_fixed = {
            "candidate_pool_af_gate": "target_plus_or_minus_0.025_inclusive",
            "conditional_inference": True,
            "direct_operational_estimand": True,
            "failed_v2_inferential_use": False,
            "fixed_cessation_endpoint": 1_780.0,
            "gamma_smc_statistics_used": False,
            "manifest_path": calibration.get("manifest_path"),
            "manifest_sha256": calibration.get("manifest_sha256"),
            "normal_10mb_stdpopsim_processing": True,
            "realized_selection_end_generations_ago": 1_780.0,
            "schema": calibration.get("schema"),
            "seed_plan_path": calibration.get("seed_plan_path"),
            "seed_plan_sha256": calibration.get("seed_plan_sha256"),
            "status": calibration.get("status"),
            "strict_selected_type1_identity": True,
            "uniform_exact_k_panel": True,
            "unit_id": unit_id,
        }
        if (
            not isinstance(authorization, Mapping)
            or set(authorization)
            != {*authorization_fixed, "attempt_schedule_sha256", "event_audit"}
            or any(
                authorization.get(key) != value
                for key, value in authorization_fixed.items()
            )
            or authorization.get("unit_id") != unit_id
            or authorization.get("attempt_schedule_sha256")
            != calibration.get("attempt_schedule_sha256")
            or authorization.get("schema") != calibration.get("schema")
            or authorization.get("manifest_sha256")
            != calibration.get("manifest_sha256")
            or authorization.get("seed_plan_sha256")
            != calibration.get("seed_plan_sha256")
            or authorization.get("realized_selection_end_generations_ago") != 1_780.0
            or authorization.get("gamma_smc_statistics_used") is not False
            or authorization.get("failed_v2_inferential_use") is not False
        ):
            raise ValueError(f"focused selected authorization differs: {unit_id}")
        event_audit = authorization.get("event_audit")
        if (
            not isinstance(event_audit, Mapping)
            or event_audit.get("schema") != "gamma-smc.han-fixed-cessation-v2-events/v1"
            or event_audit.get("focal_event_count") != 7
            or event_audit.get("mutation_age_generations") != 2_400.0
            or event_audit.get("mutation_population") != "Neanderthal"
            or event_audit.get("selection_coefficient") != 0.01
            or event_audit.get("realized_selection_end_generations_ago") != 1_780.0
            or event_audit.get("serialized_extended_events_sha256")
            != "101733b4c5bebca901132824483ccb286027f556eda2f1c1bebfd165fda85042"
        ):
            raise ValueError(f"focused selected event audit differs: {unit_id}")
        if _canonical_sha256(event_audit) != (
            "2f269e4323ba42391605cf00b5972c75aa458c89a01f6b7d35d5039149ae4dd2"
        ):
            raise ValueError(f"focused selected event audit hash differs: {unit_id}")
    sim_outputs = sim.get("outputs")
    decode_outputs = decode.get("outputs")
    if (
        not isinstance(sim_outputs, Mapping)
        or set(sim_outputs) != FOCUSED_SIMULATION_OUTPUT_KEYS
        or not isinstance(decode_outputs, Mapping)
        or set(decode_outputs) != FOCUSED_DECODE_OUTPUT_KEYS
    ):
        raise ValueError(f"focused unit output inventory differs: {unit_id}")
    sim_output_records = {
        key: sim_outputs[key] for key in ("tree", "pair_table", "overall_pairs")
    }
    if not all(isinstance(record, Mapping) for record in sim_output_records.values()):
        raise TypeError(f"focused simulation output binding differs: {unit_id}")
    for key, expected_path in FOCUSED_SIMULATION_OUTPUT_PATHS.items():
        record = sim_outputs[key]
        if set(record) != {"path", "sha256"} or record.get("path") != expected_path:
            raise ValueError(
                f"focused simulation output record differs: {unit_id} {key}"
            )
        _strict_sha256(record.get("sha256"), label=f"{unit_id} simulation {key}")
    decode_static = decode_contract.get("static")
    raw_inputs = decode_contract.get("raw_inputs")
    if (
        not isinstance(decode_static, Mapping)
        or set(decode_static)
        != {
            "schema",
            "unit_record",
            "parameters",
            "pair_table",
            "implementation",
            "auxiliary_output_paths",
        }
        or decode_static.get("schema") != FOCUSED_DECODE_SCHEMA_VERSION
        or not isinstance(raw_inputs, Mapping)
        or set(raw_inputs)
        != {"posterior_path", "posterior_sha256", "metadata_path", "metadata_sha256"}
    ):
        raise ValueError(f"focused decode contract fields differ: {unit_id}")
    if (
        decode.get("software") != FOCUSED_DECODE_SOFTWARE
        or decode.get("posterior_scale") != FOCUSED_DECODE_POSTERIOR_SCALE
        or decode.get("posterior_dimensions") != FOCUSED_DECODE_POSTERIOR_DIMENSIONS
        or decode.get("focal_center") != FOCUSED_DECODE_FOCAL_CENTER
        or decode_static.get("parameters") != FOCUSED_DECODE_PARAMETERS
    ):
        raise ValueError(f"focused decode scientific metadata differs: {unit_id}")
    decode_unit = decode_static.get("unit_record")
    if not isinstance(decode_unit, Mapping) or decode_unit.get("unit_id") != unit_id:
        raise ValueError(f"decode unit binding differs: {unit_id}")
    if set(decode_unit) != {*unit, "present_ne", "decode_provenance"}:
        raise ValueError(f"decode unit fields differ: {unit_id}")
    for key in unit:
        if decode_unit.get(key) != unit[key]:
            raise ValueError(f"decode unit field differs for {unit_id}: {key}")
    if decode_unit.get("present_ne") != PRESENT_HAN_UNSCALED_NE:
        raise ValueError(f"decode present Han Ne differs: {unit_id}")
    decode_implementation = decode_static.get("implementation")
    if (
        not isinstance(decode_implementation, Mapping)
        or set(decode_implementation) != {"postprocessor_path", "postprocessor_sha256"}
        or not _paths_equivalent(
            str(decode_implementation.get("postprocessor_path")),
            repo_root / "python/gamma_smc_aou/focused_selection_decode.py",
        )
        or decode_implementation.get("postprocessor_sha256")
        != FOCUSED_POSTPROCESSOR_SHA256
    ):
        raise ValueError(f"focused decode implementation differs: {unit_id}")
    auxiliary = decode_static.get("auxiliary_output_paths")
    expected_auxiliary = {
        "decoder_stderr": decode_dir / "decoder.stderr.zst",
        "decoder_stdout": decode_dir / "decoder.stdout.zst",
        "overall_summary": decode_dir / "overall.summary.tsv",
        "run_manifest": decode_dir / "overall.summary.tsv.run.json",
    }
    if not isinstance(auxiliary, Mapping) or set(auxiliary) != set(expected_auxiliary):
        raise ValueError(f"focused decode auxiliary fields differ: {unit_id}")
    for key, expected_path in expected_auxiliary.items():
        if not _paths_equivalent(str(auxiliary.get(key)), expected_path):
            raise ValueError(f"focused decode auxiliary path differs: {unit_id} {key}")
    pair_table = decode_static.get("pair_table")
    expected_pair_table = unit_dir / "sample_manifest.tsv"
    if (
        not isinstance(pair_table, Mapping)
        or set(pair_table) != {"path", "sha256"}
        or not _paths_equivalent(str(pair_table.get("path")), expected_pair_table)
        or pair_table.get("sha256") != sim_output_records["pair_table"].get("sha256")
    ):
        raise ValueError(f"focused decode pair-table binding differs: {unit_id}")
    decode_provenance = decode_unit.get("decode_provenance")
    if (
        not isinstance(decode_provenance, Mapping)
        or set(decode_provenance) != FOCUSED_DECODE_PROVENANCE_FIELDS
        or decode_provenance.get("decoder_path") != EXPECTED_GAMMA_BINARY_IDENTITY
        or decode_provenance.get("decoder_sha256") != EXPECTED_GAMMA_BINARY_SHA256
        or decode_provenance.get("settings") != FOCUSED_DECODER_SETTINGS
        or decode_provenance.get("simulation_completion_sha256")
        != _sha256_file(sim_path)
        or decode_provenance.get("tree_sha256")
        != sim_output_records["tree"].get("sha256")
        or decode_provenance.get("pair_table_sha256")
        != sim_output_records["pair_table"].get("sha256")
        or decode_provenance.get("overall_pairs_sha256")
        != sim_output_records["overall_pairs"].get("sha256")
    ):
        raise ValueError(f"focused decode provenance binding differs: {unit_id}")
    if decode_provenance.get("implementation_sources") != FOCUSED_DECODE_SOURCE_HASHES:
        raise ValueError(f"focused decode implementation sources differ: {unit_id}")
    checks = [
        FileCheck(
            "focused_unit",
            "matched_af_or_selected",
            unit_id,
            "simulation_completion",
            "simulation_completion",
            sim_path,
            _path_identity(sim_path, repo_root),
            None,
            None,
        ),
        FileCheck(
            "focused_unit",
            "matched_af_or_selected",
            unit_id,
            "decode_completion",
            "decode_completion",
            decode_path,
            _path_identity(decode_path, repo_root),
            None,
            None,
        ),
        FileCheck(
            "focused_unit",
            "matched_af_or_selected",
            unit_id,
            "simulation_contract",
            "simulation_contract",
            simulation_contract_path,
            _path_identity(simulation_contract_path, repo_root),
            None,
            None,
        ),
    ]
    for key, record in sorted(sim_outputs.items()):
        checks.append(
            _check_from_record(
                scope="focused_unit",
                source="matched_af_or_selected",
                unit_id=unit_id,
                stage="simulation_output",
                record_key=str(key),
                base=unit_dir,
                record=record,
                repo_root=repo_root,
            )
        )
    for key, record in sorted(decode_outputs.items()):
        expected_fields = {"path", "sha256", "size_bytes"}
        if key in FOCUSED_DECODE_TABULAR_OUTPUT_ROWS:
            expected_fields |= {"rows", "columns"}
        if (
            not isinstance(record, Mapping)
            or set(record) != expected_fields
            or record.get("path") != FOCUSED_DECODE_OUTPUT_PATHS[key]
            or (
                key in FOCUSED_DECODE_TABULAR_OUTPUT_ROWS
                and (
                    record.get("rows") != FOCUSED_DECODE_TABULAR_OUTPUT_ROWS[key]
                    or record.get("columns")
                    != FOCUSED_DECODE_TABULAR_OUTPUT_COLUMNS[key]
                )
            )
        ):
            raise ValueError(f"focused decode output record differs: {unit_id} {key}")
        checks.append(
            _check_from_record(
                scope="focused_unit",
                source="matched_af_or_selected",
                unit_id=unit_id,
                stage="decode_output",
                record_key=str(key),
                base=decode_dir,
                record=record,
                repo_root=repo_root,
            )
        )
    _validate_focused_manifest(
        unit_dir / "sample_manifest.tsv",
        unit_dir / "pairs/overall.pairs.tsv",
        sim,
        unit_id=unit_id,
    )
    for key in ("posterior", "metadata"):
        declared_path = raw_inputs.get(f"{key}_path")
        sha = _strict_sha256(
            raw_inputs.get(f"{key}_sha256"), label=f"{unit_id} decode raw {key}"
        )
        if not isinstance(declared_path, str):
            raise TypeError(f"decode raw path is malformed: {unit_id} {key}")
        raw_path = decode_dir / (
            "posterior.zst" if key == "posterior" else "posterior.zst.meta"
        )
        if not _paths_equivalent(declared_path, raw_path):
            raise ValueError(f"decode raw path binding differs: {unit_id} {key}")
        _require_regular_file(raw_path, boundary=decode_dir, label="decode raw input")
        checks.append(
            FileCheck(
                "focused_unit",
                "matched_af_or_selected",
                unit_id,
                "decode_raw_input",
                key,
                raw_path,
                _path_identity(raw_path, repo_root),
                sha,
                None,
            )
        )
    binding = {
        "unit_id": unit_id,
        "simulation_completion_sha256": _sha256_file(sim_path),
        "simulation_contract_sha256": sim_contract_sha,
        "simulation_outputs_sha256": _canonical_sha256(sim_outputs),
        "decode_completion_sha256": _sha256_file(decode_path),
        "decode_contract_sha256": decode_contract_sha,
        "decode_outputs_sha256": _canonical_sha256(decode_outputs),
        "decode_raw_inputs_sha256": _canonical_sha256(raw_inputs),
    }
    return checks, binding, sim


def _endpoint_from_existing_unit(
    unit: Mapping[str, Any], completion: Mapping[str, Any]
) -> dict[str, Any]:
    provenance = completion.get("provenance", {})
    gate = provenance.get("candidate_pool_frequency_gate", {})
    genotype = completion.get("genotype_counts", {})
    sample_af = float(completion.get("sample_af"))
    sample_alt = int(completion.get("sample_alt_count"))
    if not np.isclose(sample_af, TARGET_AF) or sample_alt != 40:
        raise ValueError(f"focused endpoint AF differs: {unit['unit_id']}")
    pool_af = float(gate.get("observed_af"))
    pool_alt = int(gate.get("observed_alt_count"))
    if not (0.175 <= pool_af <= 0.225) or pool_alt != round(pool_af * 1_000):
        raise ValueError(f"focused candidate-pool AF differs: {unit['unit_id']}")
    simulation_class = str(unit["simulation_class"])
    selection_end = np.nan
    if simulation_class == "selected":
        selection_end = float(
            provenance.get("selection_endpoint", {})
            .get("cessation", {})
            .get("selection_end_generations_ago")
        )
        if not np.isclose(selection_end, SELECTION_END_GENERATIONS):
            raise ValueError(f"selected endpoint timing differs: {unit['unit_id']}")
    return {
        "unit_id": str(unit["unit_id"]),
        "bank_role": (
            "selected" if simulation_class == "selected" else "matched_af_sensitivity"
        ),
        "simulation_class": simulation_class,
        "replicate_index": int(unit["replicate_index"]),
        "seed": int(unit["seed"]),
        "selection_coefficient": float(unit["selection_coefficient"]),
        "present_conditioning": "candidate_pool_af_0.175_to_0.225",
        "simulation_completion_sha256": np.nan,
        "final_census_alt_count": np.nan,
        "final_census_total_count": np.nan,
        "final_census_af": np.nan,
        "final_census_segregating": np.nan,
        "pool_diploids": CANDIDATE_POOL_DIPLOIDS,
        "final_pool_alt_count": pool_alt,
        "final_pool_total_count": 1_000,
        "final_pool_af": pool_af,
        "pool_detected": pool_alt > 0,
        "pool_fixed": pool_alt == 1_000,
        "sample_diploids": SAMPLE_DIPLOIDS,
        "sample_alt_count": sample_alt,
        "sample_total_count": 200,
        "sample_af": sample_af,
        "sample_detected": True,
        "sample_fixed": sample_alt == 200,
        "panel_seed": np.nan,
        "hom_ref_diploids": int(genotype.get("hom_ref")),
        "heterozygous_diploids": int(genotype.get("heterozygous")),
        "hom_alt_diploids": int(genotype.get("hom_alt")),
        "mutation_origin_generations": MUTATION_ORIGIN_GENERATIONS,
        "introgression_pulse_generations": INTROGRESSION_PULSE_GENERATIONS,
        "selection_end_generations": selection_end,
        "realized_frequency_trajectory_persisted": False,
    }


def _focused_class_summaries_and_scores(
    campaign_root: Path,
    units: pd.DataFrame,
    genotype_counts: Mapping[str, Mapping[str, Any]],
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    frames: dict[str, list[pd.DataFrame]] = {"gamma_smc": [], "tree_truth": []}
    unit_lookup = {
        str(record["unit_id"]): record for record in units.to_dict(orient="records")
    }
    for unit_id, unit in unit_lookup.items():
        counts = _focused_genotype_counts(
            genotype_counts[unit_id], label=f"{unit_id} class-summary"
        )
        expected_n = {"overall": 100, **counts}
        paths = {
            "gamma_smc": campaign_root
            / "work"
            / unit_id
            / "decoded/class_summaries.tsv",
            "tree_truth": campaign_root
            / "work"
            / unit_id
            / "truth_class_summaries.tsv",
        }
        for source, path in paths.items():
            frame = _read_tsv(path, label=f"{unit_id} {source} class summaries")
            expected_columns = (
                FOCUSED_GAMMA_CLASS_SUMMARY_COLUMNS
                if source == "gamma_smc"
                else FOCUSED_TRUTH_CLASS_SUMMARY_COLUMNS
            )
            if tuple(frame.columns) != expected_columns or len(frame) != 28:
                raise ValueError(f"focused class-summary schema differs: {unit_id}")
            expected_grid = {
                (genotype, threshold)
                for genotype in ("overall", "hom_ref", "heterozygous", "hom_alt")
                for threshold in THRESHOLDS_YEARS
            }
            observed_grid = set(
                zip(
                    frame["genotype_class"].astype(str),
                    _exact_integer(frame, "threshold_years"),
                    strict=True,
                )
            )
            if observed_grid != expected_grid:
                raise ValueError(f"focused class-summary grid differs: {unit_id}")
            exact_metadata = {
                "unit_id": unit_id,
                "demography_id": unit["demography_id"],
                "demography_kind": unit["demography_kind"],
                "population": unit["population"],
                "simulation_class": unit["simulation_class"],
                "replicate_index": int(unit["replicate_index"]),
                "seed": int(unit["seed"]),
                "source": source,
                "focal_output_position_0based": FOCAL_POSITION_BP,
                "focal_offset_bp": 0,
            }
            for key, expected in exact_metadata.items():
                observed = frame[key]
                if not (observed == expected).all():
                    raise ValueError(
                        f"focused class-summary metadata differs: {unit_id} {key}"
                    )
            for key in ("selection_coefficient", "target_allele_frequency"):
                if not np.all(
                    pd.to_numeric(frame[key]).to_numpy(dtype=float) == float(unit[key])
                ):
                    raise ValueError(
                        f"focused class-summary metadata differs: {unit_id} {key}"
                    )
            for genotype, expected in expected_n.items():
                observed = _exact_integer(
                    frame[frame["genotype_class"] == genotype], "n_pairs"
                )
                if set(observed) != {expected}:
                    raise ValueError(
                        f"focused class-summary pair count differs: {unit_id}"
                    )
            probability = _exact_numeric(frame, "focal_mean_p_tmrca_lt_threshold")
            if np.any((probability < 0) | (probability > 1)):
                raise ValueError(f"focused focal probability differs: {unit_id}")
            for _, curve in frame.groupby("genotype_class", sort=False):
                ordered = curve.sort_values("threshold_years")
                if np.any(
                    np.diff(
                        ordered["focal_mean_p_tmrca_lt_threshold"].to_numpy(dtype=float)
                    )
                    < -1e-12
                ):
                    raise ValueError(f"focused focal CDF differs: {unit_id}")
            frames[source].append(frame)
    combined = {
        source: pd.concat(parts, ignore_index=True)
        .sort_values(["unit_id", "genotype_class", "threshold_years"], kind="stable")
        .reset_index(drop=True)
        for source, parts in frames.items()
    }
    score_rows: list[dict[str, Any]] = []
    for source, frame in combined.items():
        focal = frame.pivot(
            index=["unit_id", "threshold_years"],
            columns="genotype_class",
            values="focal_mean_p_tmrca_lt_threshold",
        ).reset_index()
        if len(focal) != 110 * len(THRESHOLDS_YEARS):
            raise ValueError(f"focused direct score inventory differs: {source}")
        for record in focal.to_dict(orient="records"):
            unit = unit_lookup[str(record["unit_id"])]
            hom_ref = float(record["hom_ref"])
            hom_alt = float(record["hom_alt"])
            score_rows.append(
                {
                    "source": source,
                    "unit_id": record["unit_id"],
                    "simulation_class": unit["simulation_class"],
                    "selection_coefficient": float(unit["selection_coefficient"]),
                    "target_allele_frequency": float(unit["target_allele_frequency"]),
                    "replicate_index": int(unit["replicate_index"]),
                    "panel_n_diploid_samples": 100,
                    "threshold_years": int(record["threshold_years"]),
                    "p_overall": float(record["overall"]),
                    "p_hom_ref": hom_ref,
                    "p_hom_alt": hom_alt,
                    "p_hom_alt_minus_hom_ref": hom_alt - hom_ref,
                    "bank_role": (
                        "selected"
                        if unit["simulation_class"] == "selected"
                        else "matched_af_sensitivity"
                    ),
                    "score_scope": "focal_nearest_within_diploid",
                }
            )
    scores = pd.DataFrame(score_rows).sort_values(
        ["source", "bank_role", "unit_id", "threshold_years"], kind="stable"
    )
    return combined, scores.reset_index(drop=True)


def _load_existing_unit_bundle(
    repo_root: Path, campaign_root: Path, *, workers: int
) -> ExistingUnitBundle:
    design_path = campaign_root / "study_design.json"
    execution_path = campaign_root / "execution_units.tsv"
    _require_regular_file(design_path, boundary=campaign_root, label="study design")
    _require_regular_file(
        execution_path, boundary=campaign_root, label="execution inventory"
    )
    design = _read_json(design_path)
    _validate_study_design(design)
    execution = _read_tsv(execution_path, label="execution inventory")
    selected = _select_existing_units(execution)
    expected_execution_sha = design.get("hashes", {}).get("execution_units_tsv_sha256")
    if expected_execution_sha != _sha256_file(execution_path):
        raise ValueError("execution inventory hash differs from study design")
    checks = [
        FileCheck(
            "focused_campaign",
            "study_design",
            "",
            "campaign_input",
            "study_design",
            design_path,
            _path_identity(design_path, repo_root),
            None,
            None,
        ),
        FileCheck(
            "focused_campaign",
            "execution_units",
            "",
            "campaign_input",
            "execution_units",
            execution_path,
            _path_identity(execution_path, repo_root),
            expected_execution_sha,
            None,
        ),
    ]
    bindings = []
    endpoints = []
    genotype_counts: dict[str, Mapping[str, Any]] = {}
    for unit in selected.to_dict(orient="records"):
        unit_checks, binding, completion = _unit_receipt_checks(
            repo_root, campaign_root, unit
        )
        checks.extend(unit_checks)
        bindings.append(binding)
        endpoint = _endpoint_from_existing_unit(unit, completion)
        endpoint["simulation_completion_sha256"] = binding[
            "simulation_completion_sha256"
        ]
        endpoints.append(endpoint)
        genotype_counts[str(unit["unit_id"])] = completion["genotype_counts"]
    audit = _audit_file_checks(checks, workers=workers)
    class_summaries, unit_scores = _focused_class_summaries_and_scores(
        campaign_root, selected, genotype_counts
    )
    return ExistingUnitBundle(
        units=selected,
        endpoints=pd.DataFrame(endpoints).sort_values(
            ["bank_role", "replicate_index"], kind="stable"
        ),
        unit_scores=unit_scores,
        class_summaries=class_summaries,
        audit=audit,
        input_contract={
            "study_design": _file_record(design_path, repo_root),
            "execution_units": _file_record(execution_path, repo_root),
            "filtered_execution_rows_sha256": _canonical_sha256(
                selected.to_dict(orient="records")
            ),
            "selected_unit_count": SELECTED_REPLICATES,
            "matched_af_neutral_unit_count": MATCHED_AF_NEUTRAL_REPLICATES,
            "unit_ids": selected["unit_id"].astype(str).tolist(),
            "unit_receipt_inventory_sha256": _canonical_sha256(bindings),
        },
    )


def _record_for_path(
    outputs: Mapping[str, Any], relative_path: str, *, label: str
) -> tuple[str, Mapping[str, Any]]:
    matches = [
        (str(key), record)
        for key, record in outputs.items()
        if isinstance(record, Mapping) and record.get("path") == relative_path
    ]
    if len(matches) != 1:
        raise ValueError(f"{label} completion does not uniquely bind {relative_path}")
    return matches[0]


def _validate_survival_estimand(contract: Mapping[str, Any]) -> None:
    estimand = contract.get("estimand")
    if not isinstance(estimand, Mapping):
        raise TypeError("survival-neutral completion lacks an estimand contract")
    exact = {
        "source_model": SOURCE_MODEL,
        "mutation_origin_population": "Neanderthal",
        "mutation_origin_mode": "single_copy",
        "present_conditioning": "strict_segregating_in_present_han_census",
        "pool_frequency_conditioning": False,
        "terminal_af_target_or_band_conditioning": False,
        "sample_detection_conditioning": False,
        "panel_sampling": "uniform_100_of_500",
        "sample_diploids": SAMPLE_DIPLOIDS,
        "candidate_pool_diploids": CANDIDATE_POOL_DIPLOIDS,
    }
    for key, expected in exact.items():
        if estimand.get(key) != expected:
            raise ValueError(f"survival-neutral estimand differs: {key}")
    numeric = {
        "mutation_age_generations": MUTATION_ORIGIN_GENERATIONS,
        "pulse_generations": INTROGRESSION_PULSE_GENERATIONS,
        "selection_coefficient": 0.0,
    }
    expected_keys = {*exact, *numeric, "thresholds_years"}
    if set(estimand) != expected_keys:
        raise ValueError("survival-neutral estimand fields differ from frozen schema")
    for key, expected in numeric.items():
        try:
            observed = float(estimand.get(key))
        except (TypeError, ValueError) as error:
            raise ValueError(f"survival-neutral estimand differs: {key}") from error
        if not np.isclose(observed, expected):
            raise ValueError(f"survival-neutral estimand differs: {key}")
    if tuple(estimand.get("thresholds_years", ())) != THRESHOLDS_YEARS:
        raise ValueError("survival-neutral threshold grid differs")


def _expected_survival_unit_ids() -> tuple[str, ...]:
    return tuple(
        f"ancient_eurasia_single_origin_survival__neutral__rep{index:03d}"
        for index in range(SURVIVAL_NEUTRAL_REPLICATES)
    )


def _survival_seed(label: str) -> int:
    token = f"{SURVIVAL_BASE_SEED}:ancient_eurasia_single_origin_survival:{label}"
    digest = hashlib.sha256(token.encode("utf-8")).digest()
    return 1 + int.from_bytes(digest[:8], "big") % (2**31 - 2)


def _validate_survival_seed_schedule(frame: pd.DataFrame) -> None:
    seeds = _exact_integer(frame, "seed")
    panel_seeds = _exact_integer(frame, "panel_seed")
    unit_ids = _expected_survival_unit_ids()
    expected_seeds = np.asarray([_survival_seed(unit_id) for unit_id in unit_ids])
    expected_panel_seeds = np.asarray(
        [_survival_seed(f"{unit_id}:panel") for unit_id in unit_ids]
    )
    if not np.array_equal(seeds, expected_seeds) or not np.array_equal(
        panel_seeds, expected_panel_seeds
    ):
        raise ValueError("survival-neutral deterministic seed schedule differs")
    if len(set(seeds)) != len(seeds) or len(set(panel_seeds)) != len(panel_seeds):
        raise ValueError("survival-neutral deterministic seeds are not unique")
    if set(seeds).intersection(panel_seeds):
        raise ValueError("survival-neutral unit and panel seed streams overlap")


def _validate_survival_scores(frame: pd.DataFrame) -> pd.DataFrame:
    if tuple(frame.columns) != SURVIVAL_SCORE_COLUMNS:
        raise ValueError("survival-neutral score columns differ from the frozen schema")
    result = frame.copy()
    if set(result["source"].astype(str)) != {"tree_truth", "gamma_smc"}:
        raise ValueError("survival-neutral scores have incompatible sources")
    expected_id_order = _expected_survival_unit_ids()
    expected_ids = set(expected_id_order)
    if set(result["unit_id"].astype(str)) != expected_ids:
        raise ValueError("survival-neutral unit IDs are incompatible")
    if result.duplicated(["unit_id", "source", "threshold_years"]).any():
        raise ValueError("survival-neutral scores contain duplicate rows")
    if len(result) != 2 * SURVIVAL_NEUTRAL_REPLICATES * len(THRESHOLDS_YEARS):
        raise ValueError("survival-neutral score row count is incompatible")
    if set(_exact_integer(result, "threshold_years")) != set(THRESHOLDS_YEARS):
        raise ValueError("survival-neutral thresholds are incompatible")
    expected_grid = {
        (unit_id, source, threshold)
        for unit_id in expected_ids
        for source in ("tree_truth", "gamma_smc")
        for threshold in THRESHOLDS_YEARS
    }
    observed_grid = set(
        zip(
            result["unit_id"].astype(str),
            result["source"].astype(str),
            _exact_integer(result, "threshold_years"),
            strict=True,
        )
    )
    if observed_grid != expected_grid:
        raise ValueError("survival-neutral score grid is incomplete")
    expected_replicates = {
        unit_id: index for index, unit_id in enumerate(expected_id_order)
    }
    observed_replicates = _exact_integer(result, "replicate_index")
    if any(
        observed != expected_replicates[unit_id]
        for unit_id, observed in zip(
            result["unit_id"].astype(str), observed_replicates, strict=True
        )
    ):
        raise ValueError("survival-neutral replicate index differs from unit ID")
    if result[["unit_id", "seed"]].drop_duplicates()["unit_id"].duplicated().any():
        raise ValueError("survival-neutral seed changes across score rows")
    _exact_integer(result, "seed")
    probability = _exact_numeric(result, "overall_focal_p_tmrca_lt_threshold")
    if np.any((probability < 0) | (probability > 1)):
        raise ValueError("survival-neutral focal probabilities lie outside [0, 1]")
    for _, curve in result.groupby(["unit_id", "source"], sort=False):
        ordered = curve.sort_values("threshold_years")
        values = ordered["overall_focal_p_tmrca_lt_threshold"].to_numpy(dtype=float)
        if np.any(np.diff(values) < -1e-10):
            raise ValueError("survival-neutral focal CDF is not nondecreasing")
    pool_alt = _exact_integer(result, "final_pool_alt_count")
    pool_total = _exact_integer(result, "final_pool_total_count")
    sample_alt = _exact_integer(result, "sample_alt_count")
    sample_total = _exact_integer(result, "sample_total_count")
    if set(pool_total) != {2 * CANDIDATE_POOL_DIPLOIDS}:
        raise ValueError("survival-neutral pool denominator is incompatible")
    if np.any((pool_alt < 0) | (pool_alt > pool_total)):
        raise ValueError("survival-neutral observed pool allele count is invalid")
    if set(sample_total) != {2 * SAMPLE_DIPLOIDS}:
        raise ValueError("survival-neutral sample denominator is incompatible")
    if np.any((sample_alt < 0) | (sample_alt > sample_total)):
        raise ValueError("survival-neutral sample allele count is invalid")
    if not np.allclose(_exact_numeric(result, "final_pool_af"), pool_alt / pool_total):
        raise ValueError("survival-neutral pool AF is inconsistent")
    if not np.allclose(_exact_numeric(result, "sample_af"), sample_alt / sample_total):
        raise ValueError("survival-neutral sample AF is inconsistent")
    detected = result["sample_detected"]
    if not all(isinstance(value, (bool, np.bool_)) for value in detected):
        raise ValueError("survival-neutral sample-detection field is not boolean")
    if not np.array_equal(detected.to_numpy(dtype=bool), sample_alt > 0):
        raise ValueError("survival-neutral sample-detection field is inconsistent")
    endpoint_columns = [
        "unit_id",
        "replicate_index",
        "seed",
        "final_pool_alt_count",
        "final_pool_total_count",
        "final_pool_af",
        "sample_alt_count",
        "sample_total_count",
        "sample_af",
        "sample_detected",
    ]
    endpoint_rows = result[endpoint_columns].drop_duplicates()
    if len(endpoint_rows) != SURVIVAL_NEUTRAL_REPLICATES:
        raise ValueError("survival-neutral endpoint fields change across score rows")
    return result.sort_values(
        ["source", "unit_id", "threshold_years"], kind="stable"
    ).reset_index(drop=True)


def _strict_boolean_column(frame: pd.DataFrame, column: str) -> np.ndarray:
    values = frame[column]
    if not all(isinstance(value, (bool, np.bool_)) for value in values):
        raise ValueError(f"boolean column is invalid: {column}")
    return values.to_numpy(dtype=bool)


def _validate_survival_endpoints(frame: pd.DataFrame) -> pd.DataFrame:
    if tuple(frame.columns) != SURVIVAL_ENDPOINT_COLUMNS:
        raise ValueError("survival-neutral endpoint columns differ from frozen schema")
    result = frame.copy()
    if (
        len(result) != SURVIVAL_NEUTRAL_REPLICATES
        or result["unit_id"].duplicated().any()
    ):
        raise ValueError("survival-neutral endpoints are not exactly one row per unit")
    expected_ids = _expected_survival_unit_ids()
    if set(result["unit_id"].astype(str)) != set(expected_ids):
        raise ValueError("survival-neutral endpoint unit IDs are incompatible")
    result = result.sort_values("unit_id", kind="stable").reset_index(drop=True)
    if not np.array_equal(
        _exact_integer(result, "replicate_index"),
        np.arange(SURVIVAL_NEUTRAL_REPLICATES),
    ):
        raise ValueError("survival-neutral endpoint replicate indices are incompatible")
    _validate_survival_seed_schedule(result)
    for value in result["simulation_completion_sha256"]:
        _strict_sha256(value, label="survival-neutral simulation completion")

    census_alt = _exact_integer(result, "final_census_alt_count")
    census_total = _exact_integer(result, "final_census_total_count")
    if set(census_total) != {PRESENT_HAN_SCALED_CENSUS_GENOMES} or np.any(
        (census_alt <= 0) | (census_alt >= census_total)
    ):
        raise ValueError(
            "survival-neutral complete simulated Q=5 present Han census differs or "
            "is not strictly segregating"
        )
    if not np.allclose(
        _exact_numeric(result, "final_census_af"), census_alt / census_total
    ):
        raise ValueError("survival-neutral census AF is inconsistent")
    if not _strict_boolean_column(result, "final_census_segregating").all():
        raise ValueError("survival-neutral census-segregating flag is false")

    pool_diploids = _exact_integer(result, "pool_diploids")
    pool_alt = _exact_integer(result, "pool_alt_count")
    pool_total = _exact_integer(result, "pool_total_count")
    if set(pool_diploids) != {CANDIDATE_POOL_DIPLOIDS} or not np.array_equal(
        pool_total, 2 * pool_diploids
    ):
        raise ValueError("survival-neutral observational pool size is incompatible")
    if np.any((pool_alt < 0) | (pool_alt > pool_total)):
        raise ValueError("survival-neutral observational pool count is invalid")
    if not np.allclose(_exact_numeric(result, "pool_af"), pool_alt / pool_total):
        raise ValueError("survival-neutral observational pool AF is inconsistent")
    if not np.array_equal(
        _strict_boolean_column(result, "pool_detected"), pool_alt > 0
    ) or not np.array_equal(
        _strict_boolean_column(result, "pool_fixed"), pool_alt == pool_total
    ):
        raise ValueError("survival-neutral observational pool flags are inconsistent")

    sample_diploids = _exact_integer(result, "sample_diploids")
    sample_alt = _exact_integer(result, "sample_alt_count")
    sample_total = _exact_integer(result, "sample_total_count")
    if set(sample_diploids) != {SAMPLE_DIPLOIDS} or not np.array_equal(
        sample_total, 2 * sample_diploids
    ):
        raise ValueError("survival-neutral sample size is incompatible")
    if np.any((sample_alt < 0) | (sample_alt > sample_total)):
        raise ValueError("survival-neutral sample count is invalid")
    if not np.allclose(_exact_numeric(result, "sample_af"), sample_alt / sample_total):
        raise ValueError("survival-neutral sample AF is inconsistent")
    if not np.array_equal(
        _strict_boolean_column(result, "sample_detected"), sample_alt > 0
    ) or not np.array_equal(
        _strict_boolean_column(result, "sample_fixed"), sample_alt == sample_total
    ):
        raise ValueError("survival-neutral sample flags are inconsistent")
    return result


def _validate_survival_unit_endpoint(
    endpoint: Any, unit: Mapping[str, Any]
) -> dict[str, Any]:
    expected_keys = set(SURVIVAL_ENDPOINT_COLUMNS).difference(
        {"unit_id", "replicate_index", "seed", "simulation_completion_sha256"}
    ) | {"chosen_panel_rows_zero_based"}
    if not isinstance(endpoint, Mapping) or set(endpoint) != expected_keys:
        raise ValueError("survival-neutral simulation endpoint fields differ")
    chosen = endpoint.get("chosen_panel_rows_zero_based")
    if (
        not isinstance(chosen, list)
        or len(chosen) != SAMPLE_DIPLOIDS
        or len(set(chosen)) != SAMPLE_DIPLOIDS
        or chosen != sorted(chosen)
        or any(
            isinstance(value, bool)
            or not isinstance(value, int)
            or value < 0
            or value >= CANDIDATE_POOL_DIPLOIDS
            for value in chosen
        )
    ):
        raise ValueError("survival-neutral chosen panel rows are invalid")
    expected_chosen = sorted(
        np.random.default_rng(int(endpoint["panel_seed"]))
        .choice(CANDIDATE_POOL_DIPLOIDS, size=SAMPLE_DIPLOIDS, replace=False)
        .tolist()
    )
    if chosen != expected_chosen:
        raise ValueError(
            "survival-neutral chosen panel rows do not reproduce from panel seed"
        )
    row = {
        "unit_id": str(unit["unit_id"]),
        "replicate_index": int(unit["replicate_index"]),
        "seed": int(unit["seed"]),
        "simulation_completion_sha256": "0" * 64,
        **{key: endpoint[key] for key in SURVIVAL_ENDPOINT_COLUMNS[4:]},
    }
    frame = pd.DataFrame([row], columns=SURVIVAL_ENDPOINT_COLUMNS)
    census_alt = _exact_integer(frame, "final_census_alt_count")
    census_total = _exact_integer(frame, "final_census_total_count")
    if census_total[0] != PRESENT_HAN_SCALED_CENSUS_GENOMES or not (
        0 < census_alt[0] < census_total[0]
    ):
        raise ValueError(
            "survival-neutral unit complete simulated Q=5 census differs or is not "
            "strictly segregating"
        )
    if not np.isclose(
        float(frame.loc[0, "final_census_af"]), census_alt[0] / census_total[0]
    ):
        raise ValueError("survival-neutral unit census AF is inconsistent")
    if not _strict_boolean_column(frame, "final_census_segregating")[0]:
        raise ValueError("survival-neutral unit census flag is false")
    for prefix, diploids in (
        ("pool", CANDIDATE_POOL_DIPLOIDS),
        ("sample", SAMPLE_DIPLOIDS),
    ):
        if int(frame.loc[0, f"{prefix}_diploids"]) != diploids:
            raise ValueError(f"survival-neutral unit {prefix} size differs")
        alt = int(frame.loc[0, f"{prefix}_alt_count"])
        total = int(frame.loc[0, f"{prefix}_total_count"])
        if total != 2 * diploids or not 0 <= alt <= total:
            raise ValueError(f"survival-neutral unit {prefix} count differs")
        if not np.isclose(float(frame.loc[0, f"{prefix}_af"]), alt / total):
            raise ValueError(f"survival-neutral unit {prefix} AF differs")
        if bool(frame.loc[0, f"{prefix}_detected"]) != (alt > 0) or bool(
            frame.loc[0, f"{prefix}_fixed"]
        ) != (alt == total):
            raise ValueError(f"survival-neutral unit {prefix} flags differ")
    if int(endpoint["panel_seed"]) != int(unit["panel_seed"]):
        raise ValueError("survival-neutral unit endpoint panel seed differs")
    return dict(endpoint)


def _verify_survival_score_endpoints(
    scores: pd.DataFrame, endpoints: pd.DataFrame
) -> None:
    score_columns = [
        "unit_id",
        "replicate_index",
        "seed",
        "final_pool_alt_count",
        "final_pool_total_count",
        "final_pool_af",
        "sample_alt_count",
        "sample_total_count",
        "sample_af",
        "sample_detected",
    ]
    score_rows = scores[score_columns].drop_duplicates().sort_values("unit_id")
    if len(score_rows) != SURVIVAL_NEUTRAL_REPLICATES:
        raise ValueError("survival-neutral score endpoints are not one row per unit")
    endpoint_rows = endpoints[
        [
            "unit_id",
            "replicate_index",
            "seed",
            "pool_alt_count",
            "pool_total_count",
            "pool_af",
            "sample_alt_count",
            "sample_total_count",
            "sample_af",
            "sample_detected",
        ]
    ].rename(
        columns={
            "pool_alt_count": "final_pool_alt_count",
            "pool_total_count": "final_pool_total_count",
            "pool_af": "final_pool_af",
        }
    )
    endpoint_rows = endpoint_rows.sort_values("unit_id")
    for column in (
        "unit_id",
        "replicate_index",
        "seed",
        "final_pool_alt_count",
        "final_pool_total_count",
        "sample_alt_count",
        "sample_total_count",
        "sample_detected",
    ):
        if not np.array_equal(
            score_rows[column].to_numpy(), endpoint_rows[column].to_numpy()
        ):
            raise ValueError(f"survival-neutral endpoint differs from scores: {column}")
    for column in ("final_pool_af", "sample_af"):
        if not np.allclose(
            score_rows[column].to_numpy(dtype=float),
            endpoint_rows[column].to_numpy(dtype=float),
        ):
            raise ValueError(f"survival-neutral endpoint differs from scores: {column}")


def _validate_frozen_output_mapping(
    outputs: Any,
    expected_paths: Mapping[str, str],
    *,
    scope: str,
    source: str,
    unit_id: str,
    stage: str,
    base: Path,
    repo_root: Path,
    expected_rows: Mapping[str, int] | None = None,
) -> list[FileCheck]:
    if not isinstance(outputs, Mapping) or set(outputs) != set(expected_paths):
        raise ValueError(f"{stage} output labels differ from frozen schema")
    checks = []
    expected_rows = expected_rows or {}
    for key, expected_path in expected_paths.items():
        record = outputs[key]
        if not isinstance(record, Mapping) or set(record) != {
            "path",
            "sha256",
            "size_bytes",
            "rows",
        }:
            raise ValueError(f"{stage} output record differs: {key}")
        if record.get("path") != expected_path:
            raise ValueError(f"{stage} output path differs: {key}")
        row_count = record.get("rows")
        tabular = expected_path.endswith((".tsv", ".tsv.gz"))
        if tabular:
            if (
                isinstance(row_count, bool)
                or not isinstance(row_count, int)
                or row_count < 0
            ):
                raise ValueError(f"{stage} row declaration is invalid: {key}")
        elif row_count is not None:
            raise ValueError(f"{stage} non-tabular row declaration is not null: {key}")
        if key in expected_rows and row_count != expected_rows[key]:
            raise ValueError(f"{stage} row declaration differs: {key}")
        checks.append(
            _check_from_record(
                scope=scope,
                source=source,
                unit_id=unit_id,
                stage=stage,
                record_key=key,
                base=base,
                record=record,
                repo_root=repo_root,
            )
        )
    return checks


def _portable_path_identity(value: str | Path) -> str:
    token = str(value).replace("\\", "/")
    match = re.fullmatch(r"/mnt/([A-Za-z])/(.*)", token)
    if match:
        return f"{match.group(1).lower()}:/{match.group(2)}".rstrip("/")
    drive_match = re.fullmatch(r"([A-Za-z]):/(.*)", token)
    if drive_match:
        return f"{drive_match.group(1).lower()}:/{drive_match.group(2)}".rstrip("/")
    return Path(token).resolve().as_posix().rstrip("/")


def _paths_equivalent(left: str | Path, right: str | Path) -> bool:
    return _portable_path_identity(left) == _portable_path_identity(right)


def _validate_survival_decode_raw_inputs(
    raw_inputs: Any,
    *,
    unit_id: str,
    unit_dir: Path,
    decode_dir: Path,
    repo_root: Path,
) -> list[FileCheck]:
    expected = {
        "simulation_tree": ("../simulation.trees", unit_dir / "simulation.trees", None),
        "sample_manifest": (
            "../sample_manifest.tsv",
            unit_dir / "sample_manifest.tsv",
            SAMPLE_DIPLOIDS,
        ),
        "overall_pairs": (
            "overall.pairs.tsv",
            decode_dir / "overall.pairs.tsv",
            SAMPLE_DIPLOIDS,
        ),
    }
    if not isinstance(raw_inputs, Mapping) or set(raw_inputs) != set(expected):
        raise ValueError(f"survival-neutral decode raw inputs differ: {unit_id}")
    checks = []
    for key, (relative, path, rows) in expected.items():
        record = raw_inputs[key]
        if not isinstance(record, Mapping) or set(record) != {
            "path",
            "sha256",
            "size_bytes",
            "rows",
        }:
            raise ValueError(
                f"survival-neutral decode raw-input record differs: {unit_id}:{key}"
            )
        if record.get("path") != relative or record.get("rows") != rows:
            raise ValueError(
                f"survival-neutral decode raw-input declaration differs: {unit_id}:{key}"
            )
        checks.append(
            FileCheck(
                "survival_unit",
                "single_origin_census_survival",
                unit_id,
                "decode_raw_input",
                key,
                path,
                _path_identity(path, repo_root),
                _strict_sha256(
                    record.get("sha256"),
                    label=f"survival decode raw input {unit_id}:{key}",
                ),
                _strict_positive_integer(
                    record.get("size_bytes"),
                    label=f"survival decode raw input size {unit_id}:{key}",
                ),
            )
        )
    return checks


def _expected_survival_decoder_command(repo_root: Path, decode_dir: Path) -> list[str]:
    settings = SURVIVAL_DECODER_SETTINGS
    return [
        str((repo_root / EXPECTED_GAMMA_BINARY_IDENTITY).resolve()),
        "--input",
        "/dev/stdin",
        "--input_format",
        "vcf",
        "--scaled_mutation_rate",
        str(settings["scaled_mutation_rate"]),
        "--recombination_to_mutation_ratio",
        str(settings["recombination_to_mutation_ratio"]),
        "--unscaled_mutation_rate",
        str(settings["unscaled_mutation_rate"]),
        "--recent_threshold_years",
        ",".join(str(value) for value in THRESHOLDS_YEARS),
        "--generation_time",
        str(settings["generation_time_years"]),
        "--recent_summary",
        str((decode_dir / "overall.summary.tsv").resolve()),
        "--recent_call",
        str(settings["recent_call"]),
        "--recent_call_probability",
        str(settings["recent_call_probability"]),
        "--output_at_hets=false",
        "--output_at_stride",
        str(settings["output_at_stride"]),
        "--cache_size",
        str(settings["cache_size"]),
        "--threads",
        str(settings["threads"]),
        "--pair_block",
        str(settings["pair_block"]),
        "--backward_alignment",
        str(settings["backward_alignment"]),
        "--exp10",
        str(settings["exp10"]),
        "--pairs_file",
        str((decode_dir / "overall.pairs.tsv").resolve()),
        "--output",
        str((decode_dir / "posterior.zst").resolve()),
    ]


def _require_exact_portable_command(
    observed: Sequence[str],
    expected: Sequence[str],
    *,
    path_indices: set[int],
    label: str,
) -> None:
    if len(observed) != len(expected):
        raise ValueError(f"{label} token count differs")
    for index, (observed_token, expected_token) in enumerate(zip(observed, expected)):
        matches = (
            _paths_equivalent(observed_token, expected_token)
            if index in path_indices
            else observed_token == expected_token
        )
        if not matches:
            raise ValueError(f"{label} token differs at index {index}")


def _validate_survival_decoder_run(
    path: Path,
    *,
    unit_id: str,
    raw_inputs: Mapping[str, Any],
    unit_dir: Path,
    decode_dir: Path,
    repo_root: Path,
) -> None:
    record = _read_json(path)
    if set(record) != {"schema", "unit_id", "run", "raw_inputs"}:
        raise ValueError(f"survival-neutral decoder-run fields differ: {unit_id}")
    if (
        record.get("schema") != SURVIVAL_DECODER_RUN_SCHEMA_VERSION
        or record.get("unit_id") != unit_id
        or record.get("raw_inputs") != raw_inputs
    ):
        raise ValueError(f"survival-neutral decoder-run binding differs: {unit_id}")
    run = record.get("run")
    expected_run_keys = {
        "command",
        "input_path",
        "input_format",
        "tree_sequence_vcf_producer_command",
        "stdout",
        "stderr",
        "decode_seconds",
        "stride_bp",
        "cache_size_bp",
        "n_output_positions",
        "pairs_manifest",
        "n_pairs_recorded",
    }
    if not isinstance(run, Mapping) or set(run) != expected_run_keys:
        raise ValueError(f"survival-neutral decoder run is absent: {unit_id}")
    command_raw = run.get("command")
    if not isinstance(command_raw, list) or not command_raw:
        raise ValueError(f"survival-neutral decoder command is absent: {unit_id}")
    command = [str(value) for value in command_raw]
    expected_command = _expected_survival_decoder_command(repo_root, decode_dir)
    _require_exact_portable_command(
        command,
        expected_command,
        path_indices={
            0,
            expected_command.index("--recent_summary") + 1,
            expected_command.index("--pairs_file") + 1,
            expected_command.index("--output") + 1,
        },
        label=f"survival-neutral decoder command for {unit_id}",
    )
    producer_raw = run.get("tree_sequence_vcf_producer_command")
    if not isinstance(producer_raw, list):
        raise TypeError(f"survival-neutral streamed tree input differs: {unit_id}")
    producer = [str(value) for value in producer_raw]
    producer_script = (
        "import sys; "
        "from gamma_smc_aou.tree_sequence import stream_tree_sequence_vcf; "
        "stream_tree_sequence_vcf(sys.argv[1], sys.stdout, input_format=sys.argv[2])"
    )
    expected_producer = [
        str(repo_root / EXPECTED_SURVIVAL_PRODUCER_PYTHON_IDENTITY),
        "-c",
        producer_script,
        str((unit_dir / "simulation.trees").resolve()),
        "trees",
    ]
    _require_exact_portable_command(
        producer,
        expected_producer,
        path_indices={0, 3},
        label=f"survival-neutral tree-to-VCF command for {unit_id}",
    )
    if (
        run.get("input_format") != "trees"
        or not _paths_equivalent(
            str(run.get("input_path")), unit_dir / "simulation.trees"
        )
        or run.get("stride_bp") != SURVIVAL_DECODER_SETTINGS["output_at_stride"]
        or run.get("cache_size_bp") != SURVIVAL_DECODER_SETTINGS["cache_size"]
        or run.get("n_output_positions")
        != SURVIVAL_SIMULATION_PARAMETERS["sequence_length_bp"]
        // SURVIVAL_DECODER_SETTINGS["output_at_stride"]
        or run.get("pairs_manifest") is not None
        or run.get("n_pairs_recorded") is not None
    ):
        raise ValueError(f"survival-neutral streamed tree input differs: {unit_id}")
    seconds = run.get("decode_seconds")
    if (
        isinstance(seconds, bool)
        or not isinstance(seconds, (int, float))
        or not np.isfinite(seconds)
        or seconds < 0
    ):
        raise ValueError(f"survival-neutral decoder elapsed time differs: {unit_id}")
    if not isinstance(run.get("stdout"), str) or not isinstance(run.get("stderr"), str):
        raise TypeError(f"survival-neutral decoder text capture differs: {unit_id}")


def _validate_survival_execution(frame: pd.DataFrame) -> pd.DataFrame:
    if tuple(frame.columns) != SURVIVAL_EXECUTION_COLUMNS:
        raise ValueError("survival-neutral execution columns differ from frozen schema")
    if len(frame) != SURVIVAL_NEUTRAL_REPLICATES or frame["unit_id"].duplicated().any():
        raise ValueError(
            "survival-neutral execution inventory is not exactly 100 units"
        )
    result = frame.sort_values("unit_id", kind="stable").reset_index(drop=True)
    if tuple(result["unit_id"].astype(str)) != _expected_survival_unit_ids():
        raise ValueError("survival-neutral execution unit IDs are incompatible")
    if not np.array_equal(
        _exact_integer(result, "replicate_index"),
        np.arange(SURVIVAL_NEUTRAL_REPLICATES),
    ):
        raise ValueError(
            "survival-neutral execution replicate indices are incompatible"
        )
    _validate_survival_seed_schedule(result)
    if set(result["scenario"].astype(str)) != {
        "ancient_eurasia_single_origin_survival_neutral"
    }:
        raise ValueError("survival-neutral execution scenario is incompatible")
    if not np.allclose(_exact_numeric(result, "selection_coefficient"), 0.0):
        raise ValueError("survival-neutral execution selection differs")
    if set(result["mutation_origin_population"].astype(str)) != {"Neanderthal"}:
        raise ValueError("survival-neutral execution origin population differs")
    for column, expected in (
        ("mutation_age_generations", MUTATION_ORIGIN_GENERATIONS),
        ("pulse_generations", INTROGRESSION_PULSE_GENERATIONS),
    ):
        if not np.allclose(_exact_numeric(result, column), expected):
            raise ValueError(f"survival-neutral execution timing differs: {column}")
    if set(result["present_conditioning"].astype(str)) != {
        "strict_segregating_in_present_han_census"
    }:
        raise ValueError("survival-neutral execution census conditioning differs")
    for column in (
        "pool_frequency_conditioning",
        "terminal_af_target_or_band_conditioning",
        "sample_detection_conditioning",
    ):
        if _strict_boolean_column(result, column).any():
            raise ValueError(f"survival-neutral execution gate is enabled: {column}")
    if set(_exact_integer(result, "pool_diploids")) != {CANDIDATE_POOL_DIPLOIDS}:
        raise ValueError("survival-neutral execution pool size differs")
    if set(_exact_integer(result, "sample_diploids")) != {SAMPLE_DIPLOIDS}:
        raise ValueError("survival-neutral execution sample size differs")
    return result


def _survival_source_audit(
    source_hashes: Any, *, repo_root: Path
) -> tuple[pd.DataFrame, dict[str, str]]:
    if not isinstance(source_hashes, Mapping) or not source_hashes:
        raise ValueError("survival-neutral source hash inventory is absent")
    if set(source_hashes) != set(SURVIVAL_SOURCE_PATHS):
        raise ValueError("survival-neutral source hash inventory differs")
    generator_identity = "python/gamma_smc_aou/ancient_eurasia_survival_neutral.py"
    if (
        source_hashes.get(generator_identity)
        != EXPECTED_SURVIVAL_GENERATOR_SOURCE_SHA256
    ):
        raise ValueError("survival-neutral frozen generator source hash differs")
    normalized: dict[str, str] = {}
    rows = []
    for identity, expected_value in sorted(source_hashes.items()):
        if not isinstance(identity, str):
            raise TypeError("survival-neutral source identity is malformed")
        expected = _strict_sha256(expected_value, label=f"survival source {identity}")
        path = _resolve_path(repo_root, identity)
        _require_regular_file(path, boundary=repo_root, label="survival source")
        observed = _canonical_text_sha256(path)
        if observed != expected:
            raise ValueError(f"survival-neutral source hash differs: {identity}")
        size = path.stat().st_size
        rows.append(
            {
                "scope": "survival_source",
                "source": "single_origin_census_survival",
                "unit_id": "",
                "stage": "implementation",
                "record_key": identity,
                "path": _path_identity(path, repo_root),
                "declared_sha256": expected,
                "observed_sha256": observed,
                "declared_size_bytes": size,
                "observed_size_bytes": size,
                "status": "verified",
            }
        )
        normalized[identity] = expected
    return pd.DataFrame(rows), normalized


def _binary_file_check(
    record: Any,
    *,
    repo_root: Path,
    unit_id: str,
    record_key: str,
) -> FileCheck:
    if not isinstance(record, Mapping) or set(record) != {
        "identity",
        "sha256",
        "size_bytes",
    }:
        raise ValueError(f"survival-neutral binary record differs: {record_key}")
    identity = record.get("identity")
    if not isinstance(identity, str):
        raise TypeError(f"survival-neutral binary identity differs: {record_key}")
    path = _resolve_path(repo_root, identity)
    _require_regular_file(path, boundary=repo_root, label="survival binary")
    return FileCheck(
        scope="survival_binary",
        source="single_origin_census_survival",
        unit_id=unit_id,
        stage="implementation",
        record_key=record_key,
        path=path,
        path_identity=_path_identity(path, repo_root),
        declared_sha256=_strict_sha256(
            record.get("sha256"), label=f"survival binary {record_key}"
        ),
        declared_size_bytes=_strict_positive_integer(
            record.get("size_bytes"), label=f"survival binary {record_key} size"
        ),
    )


def _validate_survival_event_schedule(
    schedule: Any, declared_sha256: Any
) -> list[dict[str, Any]]:
    if not isinstance(schedule, list):
        raise TypeError("survival-neutral event schedule is absent")
    if _strict_sha256(
        declared_sha256, label="survival-neutral event schedule"
    ) != _canonical_sha256(schedule):
        raise ValueError("survival-neutral event schedule hash differs")
    condition = "ConditionOnAlleleFrequency"
    site = FOCAL_SITE_ID
    exact = "exact_generation"
    after = "generation_after_in_forward_time"
    expected = [
        {
            "event_type": "DrawMutation",
            "single_site_id": site,
            "population": "Neanderthal",
            "time": {"generations_ago": 2_400.0, "semantics": exact},
        },
        {
            "event_type": condition,
            "single_site_id": site,
            "population": "Neanderthal",
            "start_time": {"generations_ago": 2_400.0, "semantics": after},
            "end_time": {"generations_ago": 2_273.0, "semantics": exact},
            "operator": ">",
            "allele_frequency": 0.0,
        },
        {
            "event_type": condition,
            "single_site_id": site,
            "population": "Neanderthal",
            "start_time": {"generations_ago": 2_273.0, "semantics": exact},
            "end_time": {"generations_ago": 2_273.0, "semantics": exact},
            "operator": ">=",
            "allele_frequency": 1e-9,
        },
        {
            "event_type": condition,
            "single_site_id": site,
            "population": "Loschbour",
            "start_time": {"generations_ago": 2_272.0, "semantics": after},
            "end_time": {"generations_ago": 2_015.0, "semantics": exact},
            "operator": ">",
            "allele_frequency": 0.0,
        },
        {
            "event_type": condition,
            "single_site_id": site,
            "population": "Han",
            "start_time": {"generations_ago": 2_015.0, "semantics": after},
            "end_time": {"generations_ago": 0.0, "semantics": exact},
            "operator": ">",
            "allele_frequency": 0.0,
        },
        {
            "event_type": condition,
            "single_site_id": site,
            "population": "Han",
            "start_time": {"generations_ago": 0.0, "semantics": exact},
            "end_time": {"generations_ago": 0.0, "semantics": exact},
            "operator": "<",
            "allele_frequency": 1.0,
        },
    ]
    if schedule != expected:
        raise ValueError(
            "survival-neutral event schedule differs from exact six-event route"
        )
    return schedule


def _validate_survival_event_transform(transform: Any) -> None:
    expected_keys = {
        "source_constructor",
        "source_constructor_slim_scaling_factor",
        "removed_change_mutation_fitness_count",
        "removed_terminal_han_af_band_count",
        "added_terminal_han_strict_upper_count",
        "selection_coefficient_after_transform",
        "nominal_pulse_generations_ago",
        "realized_pulse_generations_ago_q5",
        "nominal_han_split_generations_ago",
        "realized_han_split_generations_ago_q5",
        "event_schedule",
    }
    if not isinstance(transform, Mapping) or set(transform) != expected_keys:
        raise ValueError("survival-neutral event-transform fields differ")
    exact = {
        "source_constructor": (
            "focused_selection_simulation.build_raw_han_selected_sweep"
        ),
        "removed_change_mutation_fitness_count": 2,
        "removed_terminal_han_af_band_count": 2,
        "added_terminal_han_strict_upper_count": 1,
    }
    for key, expected in exact.items():
        if transform.get(key) != expected:
            raise ValueError(f"survival-neutral event-transform differs: {key}")
    numeric = {
        "source_constructor_slim_scaling_factor": 5.0,
        "selection_coefficient_after_transform": 0.0,
        "nominal_pulse_generations_ago": 2_272.0,
        "realized_pulse_generations_ago_q5": 2_270.0,
        "nominal_han_split_generations_ago": 2_016.0,
        "realized_han_split_generations_ago_q5": 2_015.0,
    }
    for key, expected in numeric.items():
        try:
            observed = float(transform.get(key))
        except (TypeError, ValueError) as error:
            raise ValueError(
                f"survival-neutral event-transform differs: {key}"
            ) from error
        if not np.isclose(observed, expected):
            raise ValueError(f"survival-neutral event-transform differs: {key}")
    schedule = transform.get("event_schedule")
    _validate_survival_event_schedule(schedule, _canonical_sha256(schedule))


def _validate_survival_focal_expectation(expectation: Any) -> dict[str, Any]:
    expected = {
        "schema": SURVIVAL_FOCAL_IDENTITY_SCHEMA,
        "single_site_id": FOCAL_SITE_ID,
        "focal_position_0based": FOCAL_POSITION_BP,
        "expected_mutation_type": 1,
        "expected_selection_coeff": 0.0,
        "source_population": "Neanderthal",
        "source_subpopulation": 7,
        "requested_origin_time_generations": MUTATION_ORIGIN_GENERATIONS,
        "realized_origin_time_generations": MUTATION_ORIGIN_GENERATIONS,
        "slim_scaling_factor": 5.0,
        "slim_time_rule": "cycle - realized_origin_time / Q",
    }
    if not isinstance(expectation, Mapping) or dict(expectation) != expected:
        raise ValueError("survival-neutral focal identity expectation differs")
    return dict(expectation)


def _validate_survival_slim_end_hook_patch(record: Any) -> dict[str, Any]:
    expected_keys = {
        "original_sha256",
        "patched_sha256",
        "replacement_sha256",
        "replacement_count",
        "han_population_id",
        "mutation_draw_helper_changed",
    }
    if not isinstance(record, Mapping) or set(record) != expected_keys:
        raise ValueError("survival-neutral SLiM end-hook patch fields differ")
    for key in ("original_sha256", "patched_sha256", "replacement_sha256"):
        _strict_sha256(record.get(key), label=f"survival-neutral end-hook {key}")
    if (
        record.get("replacement_count") != 1
        or record.get("han_population_id") != 5
        or record.get("mutation_draw_helper_changed") is not False
    ):
        raise ValueError("survival-neutral SLiM end-hook patch differs")
    if (
        len(
            {
                record["original_sha256"],
                record["patched_sha256"],
                record["replacement_sha256"],
            }
        )
        != 3
    ):
        raise ValueError("survival-neutral SLiM end-hook hashes are not distinct")
    return dict(record)


def _validate_survival_focal_overlay_contract(record: Any) -> dict[str, Any]:
    if (
        not isinstance(record, Mapping)
        or dict(record) != SURVIVAL_FOCAL_OVERLAY_PATCH
        or _canonical_sha256(record) != _canonical_sha256(SURVIVAL_FOCAL_OVERLAY_PATCH)
    ):
        raise ValueError("survival-neutral focal-overlay contract differs")
    return dict(record)


def _validate_survival_individual_location_contract(record: Any) -> dict[str, Any]:
    if (
        not isinstance(record, Mapping)
        or dict(record) != SURVIVAL_INDIVIDUAL_LOCATION_CANONICALIZATION
        or _canonical_sha256(record)
        != _canonical_sha256(SURVIVAL_INDIVIDUAL_LOCATION_CANONICALIZATION)
    ):
        raise ValueError("survival-neutral individual-location contract differs")
    return dict(record)


def _survival_overlay_map_record(
    positions: Sequence[float], rates: Sequence[float]
) -> dict[str, Any]:
    position_values = [float(value) for value in positions]
    rate_values = [float(value) for value in rates]
    if len(position_values) != len(rate_values) + 1:
        raise ValueError("survival-neutral focal-overlay map dimensions differ")

    def integrated_rate(left: float, right: float) -> float:
        result = 0.0
        for interval_left, interval_right, rate in zip(
            position_values[:-1], position_values[1:], rate_values, strict=True
        ):
            width = max(0.0, min(interval_right, right) - max(interval_left, left))
            result += width * rate
        return result

    length = position_values[-1]
    return {
        "positions": position_values,
        "rates": rate_values,
        "sequence_length": length,
        "total_integral": integrated_rate(0.0, length),
        "left_flank_integral": integrated_rate(0.0, float(FOCAL_POSITION_BP)),
        "focal_integral": integrated_rate(
            float(FOCAL_POSITION_BP), float(FOCAL_POSITION_BP + 1)
        ),
        "right_flank_integral": integrated_rate(float(FOCAL_POSITION_BP + 1), length),
    }


def _validate_survival_focal_overlay_receipt(record: Any) -> dict[str, Any]:
    required = {
        "schema",
        "contract",
        "contract_sha256",
        "intercepted_call_count",
        "patched_call_count",
        "all_positive_focal_rates_masked",
        "callable_restored",
        "calls",
    }
    if (
        not isinstance(record, Mapping)
        or set(record) != required
        or record.get("schema") != SURVIVAL_FOCAL_OVERLAY_PATCH_SCHEMA
    ):
        raise ValueError("survival-neutral focal-overlay receipt fields differ")
    contract = _validate_survival_focal_overlay_contract(record.get("contract"))
    if record.get("contract_sha256") != _canonical_sha256(contract):
        raise ValueError("survival-neutral focal-overlay contract hash differs")
    calls = record.get("calls")
    if (
        not isinstance(calls, list)
        or len(calls) != 2
        or isinstance(record.get("intercepted_call_count"), bool)
        or not isinstance(record.get("intercepted_call_count"), int)
        or record.get("intercepted_call_count") != 2
        or isinstance(record.get("patched_call_count"), bool)
        or not isinstance(record.get("patched_call_count"), int)
        or record.get("patched_call_count") != 1
        or record.get("all_positive_focal_rates_masked") is not True
        or record.get("callable_restored") is not True
    ):
        raise ValueError("survival-neutral focal-overlay call inventory differs")
    call_fields = {
        "call_index",
        "layer",
        "rate_kind",
        "model_class",
        "model_type",
        "model_next_id",
        "start_time",
        "end_time",
        "random_seed",
        "keep",
        "discrete_genome",
        "record_provenance",
        "non_rate_arguments_preserved_by_identity",
        "model_identity_preserved",
        "rate_argument_replaced",
        "positive_focal_rate_before",
        "positive_focal_rate_after",
        "original_focal_rate",
        "passed_focal_rate",
        "original_rate_map",
        "passed_rate_map",
        "call_succeeded",
    }
    shared_tick = calls[0].get("end_time") if isinstance(calls[0], Mapping) else None
    if (
        isinstance(shared_tick, bool)
        or not isinstance(shared_tick, (int, float))
        or not np.isfinite(shared_tick)
        or shared_tick <= 0
    ):
        raise ValueError("survival-neutral focal-overlay shared tick differs")
    random_seeds = []
    for index, (call, expected) in enumerate(
        zip(calls, contract["expected_call_inventory"], strict=True), start=1
    ):
        if not isinstance(call, Mapping) or set(call) != call_fields:
            raise ValueError("survival-neutral focal-overlay call fields differ")
        expected_start = (
            shared_tick if expected["start_time"] == "finite_shared_slim_tick" else None
        )
        expected_end = (
            shared_tick if expected["end_time"] == "finite_shared_slim_tick" else None
        )
        original = _survival_overlay_map_record(
            expected["original_positions"], expected["original_rates"]
        )
        passed = _survival_overlay_map_record(
            expected["passed_positions"], expected["passed_rates"]
        )
        if (
            isinstance(call.get("call_index"), bool)
            or not isinstance(call.get("call_index"), int)
            or call.get("call_index") != index
            or call.get("layer") != expected["layer"]
            or call.get("rate_kind") != "msprime.RateMap"
            or call.get("model_class") != expected["model_class"]
            or isinstance(call.get("model_type"), bool)
            or not isinstance(call.get("model_type"), int)
            or call.get("model_type") != expected["model_type"]
            or isinstance(call.get("model_next_id"), bool)
            or not isinstance(call.get("model_next_id"), int)
            or call["model_next_id"] < 0
            or call.get("start_time") != expected_start
            or call.get("end_time") != expected_end
            or call.get("keep") is not True
            or call.get("discrete_genome") is not None
            or call.get("record_provenance") is not True
            or call.get("non_rate_arguments_preserved_by_identity") is not True
            or call.get("model_identity_preserved") is not True
            or call.get("rate_argument_replaced") is not expected["patched"]
            or call.get("positive_focal_rate_before") is not expected["patched"]
            or call.get("positive_focal_rate_after") is not False
            or isinstance(call.get("original_focal_rate"), bool)
            or not isinstance(call.get("original_focal_rate"), (int, float))
            or call.get("original_focal_rate") != original["focal_integral"]
            or isinstance(call.get("passed_focal_rate"), bool)
            or not isinstance(call.get("passed_focal_rate"), (int, float))
            or call.get("passed_focal_rate") != passed["focal_integral"]
            or _canonical_sha256(call.get("original_rate_map"))
            != _canonical_sha256(original)
            or _canonical_sha256(call.get("passed_rate_map"))
            != _canonical_sha256(passed)
            or call.get("call_succeeded") is not True
        ):
            raise ValueError("survival-neutral focal-overlay call content differs")
        seed = call.get("random_seed")
        if isinstance(seed, bool) or not isinstance(seed, int) or not 1 <= seed < 2**32:
            raise ValueError("survival-neutral focal-overlay random seed differs")
        random_seeds.append(seed)
    if len(set(random_seeds)) != len(random_seeds):
        raise ValueError("survival-neutral focal-overlay random seeds are not unique")
    return dict(record)


def _validate_survival_array_record(
    record: Any, *, expected_dtype: str, expected_shape: Sequence[int]
) -> None:
    if not isinstance(record, Mapping) or set(record) != {"dtype", "shape", "sha256"}:
        raise ValueError("survival-neutral location array record fields differ")
    if record.get("dtype") != expected_dtype or record.get("shape") != list(
        expected_shape
    ):
        raise ValueError("survival-neutral location array layout differs")
    _strict_sha256(record.get("sha256"), label="survival-neutral location array")


def _validate_survival_individual_location_receipt(
    record: Any, *, tree: tskit.TreeSequence | None = None
) -> dict[str, Any]:
    required = {
        "schema",
        "contract",
        "contract_sha256",
        "before",
        "after",
        "stored_tree",
        "proof",
    }
    if (
        not isinstance(record, Mapping)
        or set(record) != required
        or record.get("schema") != SURVIVAL_INDIVIDUAL_LOCATION_SCHEMA
    ):
        raise ValueError("survival-neutral individual-location receipt fields differ")
    contract = _validate_survival_individual_location_contract(record.get("contract"))
    if record.get("contract_sha256") != _canonical_sha256(contract):
        raise ValueError("survival-neutral individual-location contract hash differs")
    state_fields = {
        "individual_count",
        "location_value_count",
        "unique_location_widths",
        "location_all_finite",
        "location_finite_count",
        "location_nan_count",
        "location_positive_infinity_count",
        "location_negative_infinity_count",
        "location_all_zero",
        "location_nonzero_count",
        "location_record",
        "location_offset_record",
        "individual_nonlocation_sha256",
        "other_tables_ts_metadata_reference_sha256",
        "all_except_individual_location_sha256",
    }
    states = [record.get(key) for key in ("before", "after", "stored_tree")]
    if any(
        not isinstance(state, Mapping) or set(state) != state_fields for state in states
    ):
        raise ValueError("survival-neutral individual-location state fields differ")
    before, after, stored_tree = states
    for state in states:
        category_counts = []
        for key in (
            "location_finite_count",
            "location_nan_count",
            "location_positive_infinity_count",
            "location_negative_infinity_count",
            "location_nonzero_count",
        ):
            value = state[key]
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ValueError("survival-neutral individual-location counts differ")
            category_counts.append(value)
        if (
            isinstance(state["individual_count"], bool)
            or not isinstance(state["individual_count"], int)
            or state["individual_count"] != CANDIDATE_POOL_DIPLOIDS
            or isinstance(state["location_value_count"], bool)
            or not isinstance(state["location_value_count"], int)
            or state["location_value_count"] != 3 * CANDIDATE_POOL_DIPLOIDS
            or state["unique_location_widths"] != [3]
            or sum(category_counts[:4]) != state["location_value_count"]
            or state["location_nonzero_count"] > state["location_value_count"]
            or state["location_all_finite"]
            is not (state["location_finite_count"] == state["location_value_count"])
            or not isinstance(state["location_all_zero"], bool)
            or state["location_all_zero"] is not (state["location_nonzero_count"] == 0)
        ):
            raise ValueError("survival-neutral individual-location dimensions differ")
        _validate_survival_array_record(
            state["location_record"], expected_dtype="<f8", expected_shape=[1_500]
        )
        _validate_survival_array_record(
            state["location_offset_record"],
            expected_dtype="<u4",
            expected_shape=[501],
        )
        for key in (
            "individual_nonlocation_sha256",
            "other_tables_ts_metadata_reference_sha256",
            "all_except_individual_location_sha256",
        ):
            _strict_sha256(state[key], label=f"survival-neutral location {key}")
    canonical_location = contract["canonical_flat_location_record"]
    canonical_offsets = contract["canonical_location_offset_record"]
    if (
        before["location_offset_record"] != canonical_offsets
        or after["location_offset_record"] != canonical_offsets
        or after["location_record"] != canonical_location
        or after["location_all_finite"] is not True
        or after["location_finite_count"] != 1_500
        or after["location_nan_count"] != 0
        or after["location_positive_infinity_count"] != 0
        or after["location_negative_infinity_count"] != 0
        or after["location_all_zero"] is not True
        or after["location_nonzero_count"] != 0
        or _canonical_sha256(after) != _canonical_sha256(stored_tree)
    ):
        raise ValueError("survival-neutral canonical individual locations differ")
    unchanged_keys = (
        "location_offset_record",
        "individual_nonlocation_sha256",
        "other_tables_ts_metadata_reference_sha256",
        "all_except_individual_location_sha256",
    )
    if any(before[key] != after[key] for key in unchanged_keys):
        raise ValueError(
            "survival-neutral location canonicalization changed other data"
        )
    expected_proof = {
        "location_offsets_and_widths_unchanged": True,
        "all_nonlocation_individual_columns_unchanged": True,
        "all_other_tables_ts_metadata_reference_unchanged": True,
        "only_individual_location_values_canonicalized": True,
        "stored_tree_state_bound_at_canonicalization": True,
    }
    if record.get("proof") != expected_proof:
        raise ValueError("survival-neutral individual-location proof differs")
    if tree is not None and _canonical_sha256(
        _individual_location_state(tree.dump_tables())
    ) != _canonical_sha256(stored_tree):
        raise ValueError("survival-neutral stored tree differs from location receipt")
    return dict(record)


def _validate_survival_plan_contract(contract: Mapping[str, Any]) -> None:
    expected_keys = {
        "study_id",
        "expected_units",
        "estimand",
        "simulation_parameters",
        "decoder_settings",
        "decoder_binary",
        "execution_units_sha256",
        "source_hashes",
        "event_transform",
        "runtime_versions",
        "base_seed",
        "seed_derivation",
        "simulation_law",
        "focal_overlay_patch",
        "individual_location_canonicalization",
    }
    if set(contract) != expected_keys:
        raise ValueError("survival-neutral plan contract fields differ")
    if contract.get("study_id") != "ancient_eurasia_single_origin_survival":
        raise ValueError("survival-neutral study ID differs")
    if contract.get("expected_units") != SURVIVAL_NEUTRAL_REPLICATES:
        raise ValueError("survival-neutral planned unit count differs")
    if contract.get("base_seed") != SURVIVAL_BASE_SEED:
        raise ValueError("survival-neutral plan base seed differs")
    if contract.get("seed_derivation") != SURVIVAL_SEED_DERIVATION:
        raise ValueError("survival-neutral plan seed derivation differs")
    if contract.get("simulation_law") != SURVIVAL_SIMULATION_LAW or _canonical_sha256(
        contract.get("simulation_law")
    ) != _canonical_sha256(SURVIVAL_SIMULATION_LAW):
        raise ValueError("survival-neutral direct-forward simulation law differs")
    _validate_survival_focal_overlay_contract(contract.get("focal_overlay_patch"))
    _validate_survival_individual_location_contract(
        contract.get("individual_location_canonicalization")
    )
    _validate_survival_estimand({"estimand": contract.get("estimand")})
    if contract.get("simulation_parameters") != SURVIVAL_SIMULATION_PARAMETERS:
        raise ValueError("survival-neutral simulation parameters differ")
    if contract.get("decoder_settings") != SURVIVAL_DECODER_SETTINGS:
        raise ValueError("survival-neutral decoder settings differ")
    decoder_binary = contract.get("decoder_binary")
    if not isinstance(decoder_binary, Mapping) or (
        decoder_binary.get("identity") != EXPECTED_GAMMA_BINARY_IDENTITY
        or decoder_binary.get("sha256") != EXPECTED_GAMMA_BINARY_SHA256
    ):
        raise ValueError("survival-neutral decoder binary differs from selected bank")
    runtime_versions = contract.get("runtime_versions")
    expected_runtime_keys = {
        "python",
        "numpy",
        "pandas",
        "scipy",
        "stdpopsim",
        "pyslim",
        "msprime",
        "tskit",
    }
    if (
        not isinstance(runtime_versions, Mapping)
        or set(runtime_versions) != expected_runtime_keys
        or not all(
            isinstance(value, str) and value for value in runtime_versions.values()
        )
    ):
        raise ValueError("survival-neutral runtime version inventory differs")
    _validate_survival_event_transform(contract.get("event_transform"))
    _strict_sha256(
        contract.get("execution_units_sha256"),
        label="survival-neutral execution inventory",
    )


def _validate_survival_root_contract(
    contract: Mapping[str, Any], *, plan_contract_sha256: str, execution_sha256: str
) -> None:
    expected_keys = {
        "study_id",
        "plan_sha256",
        "execution_units_sha256",
        "estimand",
        "expected_units",
        "schemas",
        "decoder_settings",
        "decoder_binary",
        "source_hashes",
        "simulation_law",
    }
    if set(contract) != expected_keys:
        raise ValueError("survival-neutral results contract fields differ")
    if contract.get("study_id") != "ancient_eurasia_single_origin_survival":
        raise ValueError("survival-neutral results study ID differs")
    if contract.get("expected_units") != SURVIVAL_NEUTRAL_REPLICATES:
        raise ValueError("survival-neutral results unit count differs")
    if contract.get("simulation_law") != SURVIVAL_SIMULATION_LAW or _canonical_sha256(
        contract.get("simulation_law")
    ) != _canonical_sha256(SURVIVAL_SIMULATION_LAW):
        raise ValueError("survival-neutral results simulation law differs")
    if contract.get("plan_sha256") != plan_contract_sha256:
        raise ValueError("survival-neutral results plan hash differs")
    if contract.get("execution_units_sha256") != execution_sha256:
        raise ValueError("survival-neutral results execution hash differs")
    _validate_survival_estimand({"estimand": contract.get("estimand")})
    if contract.get("decoder_settings") != SURVIVAL_DECODER_SETTINGS:
        raise ValueError("survival-neutral results decoder settings differ")
    decoder_binary = contract.get("decoder_binary")
    if not isinstance(decoder_binary, Mapping) or (
        decoder_binary.get("identity") != EXPECTED_GAMMA_BINARY_IDENTITY
        or decoder_binary.get("sha256") != EXPECTED_GAMMA_BINARY_SHA256
    ):
        raise ValueError("survival-neutral results decoder binary differs")
    schemas = contract.get("schemas")
    expected_schemas = {
        "plan": SURVIVAL_PLAN_SCHEMA_VERSION,
        "simulation": SURVIVAL_SIMULATION_SCHEMA_VERSION,
        "decode": SURVIVAL_DECODE_SCHEMA_VERSION,
        "results": SURVIVAL_RESULTS_SCHEMA_VERSION,
    }
    if not isinstance(schemas, Mapping) or dict(schemas) != expected_schemas:
        raise ValueError("survival-neutral results schema inventory differs")


def _load_survival_plan(
    repo_root: Path, study_root: Path
) -> tuple[dict[str, Any], str, pd.DataFrame, list[FileCheck], pd.DataFrame]:
    plan_dir = study_root / "plan"
    plan_path = plan_dir / "study_plan.json"
    execution_path = plan_dir / "execution_units.tsv"
    _require_regular_file(plan_path, boundary=study_root, label="survival-neutral plan")
    plan = _read_json(plan_path)
    _require_exact_receipt_fields(
        plan,
        {
            "schema",
            "status",
            "created_utc",
            "contract",
            "contract_sha256",
            "outputs",
        },
        label="survival-neutral plan",
    )
    plan_contract, plan_contract_sha = _validate_completion_contract(
        plan,
        expected_schema=SURVIVAL_PLAN_SCHEMA_VERSION,
        label="survival-neutral plan",
    )
    _validate_survival_plan_contract(plan_contract)
    checks = [
        FileCheck(
            "survival_plan",
            "single_origin_census_survival",
            "",
            "plan_completion",
            "study_plan",
            plan_path,
            _path_identity(plan_path, repo_root),
            None,
            None,
        )
    ]
    checks.extend(
        _validate_frozen_output_mapping(
            plan.get("outputs"),
            {"execution_units": "execution_units.tsv"},
            scope="survival_plan",
            source="single_origin_census_survival",
            unit_id="",
            stage="plan_output",
            base=plan_dir,
            repo_root=repo_root,
            expected_rows={"execution_units": SURVIVAL_NEUTRAL_REPLICATES},
        )
    )
    execution_sha = _sha256_file(execution_path)
    if plan_contract.get("execution_units_sha256") != execution_sha:
        raise ValueError("survival-neutral plan execution hash differs")
    execution = _validate_survival_execution(
        _read_tsv(execution_path, label="survival-neutral execution inventory")
    )
    source_audit, source_hashes = _survival_source_audit(
        plan_contract.get("source_hashes"), repo_root=repo_root
    )
    plan_contract["source_hashes"] = source_hashes
    checks.append(
        _binary_file_check(
            plan_contract.get("decoder_binary"),
            repo_root=repo_root,
            unit_id="",
            record_key="decoder_binary",
        )
    )
    return plan_contract, plan_contract_sha, execution, checks, source_audit


def _survival_unit_receipts(
    repo_root: Path,
    study_root: Path,
    unit: Mapping[str, Any],
    inventory: Mapping[str, Any],
    *,
    plan_contract: Mapping[str, Any],
    plan_contract_sha256: str,
) -> tuple[list[FileCheck], dict[str, Any], pd.DataFrame]:
    unit_id = str(unit["unit_id"])
    unit_dir = study_root / "work" / unit_id
    decode_dir = unit_dir / "decoded"
    sim_path = unit_dir / "simulation_complete.json"
    decode_path = decode_dir / "decode_complete.json"
    _require_regular_file(
        sim_path, boundary=study_root, label="survival simulation completion"
    )
    _require_regular_file(
        decode_path, boundary=study_root, label="survival decode completion"
    )
    expected_inventory = {
        "unit_id": unit_id,
        "replicate_index": int(unit["replicate_index"]),
        "seed": int(unit["seed"]),
        "simulation_completion_path": f"work/{unit_id}/simulation_complete.json",
        "simulation_completion_sha256": _sha256_file(sim_path),
        "decode_completion_path": f"work/{unit_id}/decoded/decode_complete.json",
        "decode_completion_sha256": _sha256_file(decode_path),
    }
    if dict(inventory) != expected_inventory:
        raise ValueError(f"survival-neutral root unit inventory differs: {unit_id}")
    sim = _read_json(sim_path)
    decode = _read_json(decode_path)
    _require_exact_receipt_fields(
        sim,
        {
            "schema",
            "status",
            "contract",
            "contract_sha256",
            "elapsed_seconds",
            "endpoint",
            "focal_identity",
            "focal_overlay_patch",
            "individual_location_canonicalization",
            "slim_end_hook_patch",
            "outputs",
        },
        label=f"{unit_id} survival simulation",
    )
    _require_exact_receipt_fields(
        decode,
        {
            "schema",
            "status",
            "contract",
            "contract_sha256",
            "elapsed_seconds",
            "outputs",
        },
        label=f"{unit_id} survival decode",
    )
    sim_contract, sim_contract_sha = _validate_completion_contract(
        sim,
        expected_schema=SURVIVAL_SIMULATION_SCHEMA_VERSION,
        label=f"{unit_id} survival simulation",
    )
    decode_contract, decode_contract_sha = _validate_completion_contract(
        decode,
        expected_schema=SURVIVAL_DECODE_SCHEMA_VERSION,
        label=f"{unit_id} survival decode",
    )
    expected_sim_keys = {
        "unit_id",
        "replicate_index",
        "seed",
        "plan_sha256",
        "estimand",
        "event_schedule",
        "event_schedule_sha256",
        "simulation_parameters",
        "source_hashes",
        "panel_seed",
        "slim_binary",
        "focal_identity_expectation",
        "focal_overlay_patch",
        "individual_location_canonicalization",
        "slim_end_hook_patch",
    }
    expected_decode_keys = {
        "unit_id",
        "replicate_index",
        "seed",
        "plan_sha256",
        "simulation_completion_sha256",
        "estimand",
        "decoder_settings",
        "decoder_binary",
        "source_hashes",
        "raw_inputs",
    }
    if (
        set(sim_contract) != expected_sim_keys
        or set(decode_contract) != expected_decode_keys
    ):
        raise ValueError(f"survival-neutral unit contract fields differ: {unit_id}")
    for contract, stage in ((sim_contract, "simulation"), (decode_contract, "decode")):
        for key in ("unit_id", "replicate_index", "seed"):
            if contract.get(key) != unit[key]:
                raise ValueError(
                    f"survival-neutral {stage} unit binding differs: {unit_id}"
                )
        if contract.get("plan_sha256") != plan_contract_sha256:
            raise ValueError(f"survival-neutral {stage} plan hash differs: {unit_id}")
        _validate_survival_estimand({"estimand": contract.get("estimand")})
        if contract.get("source_hashes") != plan_contract.get("source_hashes"):
            raise ValueError(
                f"survival-neutral {stage} source hashes differ: {unit_id}"
            )
    if sim_contract.get("panel_seed") != unit["panel_seed"]:
        raise ValueError(f"survival-neutral simulation panel seed differs: {unit_id}")
    if sim_contract.get("simulation_parameters") != SURVIVAL_SIMULATION_PARAMETERS:
        raise ValueError(f"survival-neutral simulation parameters differ: {unit_id}")
    slim_binary = sim_contract.get("slim_binary")
    if not isinstance(slim_binary, Mapping) or (
        slim_binary.get("identity") != EXPECTED_SLIM_BINARY_IDENTITY
        or slim_binary.get("sha256") != EXPECTED_SLIM_BINARY_SHA256
    ):
        raise ValueError(f"survival-neutral SLiM binary differs: {unit_id}")
    _validate_survival_event_schedule(
        sim_contract.get("event_schedule"), sim_contract.get("event_schedule_sha256")
    )
    focal_expectation = _validate_survival_focal_expectation(
        sim_contract.get("focal_identity_expectation")
    )
    focal_overlay_contract = _validate_survival_focal_overlay_contract(
        sim_contract.get("focal_overlay_patch")
    )
    individual_location_contract = _validate_survival_individual_location_contract(
        sim_contract.get("individual_location_canonicalization")
    )
    if focal_overlay_contract != plan_contract.get(
        "focal_overlay_patch"
    ) or individual_location_contract != plan_contract.get(
        "individual_location_canonicalization"
    ):
        raise ValueError(
            f"survival-neutral v2 plan/simulation contract differs: {unit_id}"
        )
    focal_overlay_receipt = _validate_survival_focal_overlay_receipt(
        sim.get("focal_overlay_patch")
    )
    individual_location_receipt = _validate_survival_individual_location_receipt(
        sim.get("individual_location_canonicalization")
    )
    if (
        focal_overlay_receipt["contract"] != focal_overlay_contract
        or individual_location_receipt["contract"] != individual_location_contract
    ):
        raise ValueError(f"survival-neutral v2 simulation receipt differs: {unit_id}")
    slim_end_hook_patch = _validate_survival_slim_end_hook_patch(
        sim_contract.get("slim_end_hook_patch")
    )
    if sim.get("slim_end_hook_patch") != slim_end_hook_patch:
        raise ValueError(f"survival-neutral end-hook receipt differs: {unit_id}")
    unit_endpoint = _validate_survival_unit_endpoint(sim.get("endpoint"), unit)
    if decode_contract.get("simulation_completion_sha256") != _sha256_file(sim_path):
        raise ValueError(f"survival-neutral decode simulation hash differs: {unit_id}")
    if decode_contract.get("decoder_settings") != SURVIVAL_DECODER_SETTINGS:
        raise ValueError(f"survival-neutral decode settings differ: {unit_id}")
    if decode_contract.get("decoder_binary") != plan_contract.get("decoder_binary"):
        raise ValueError(f"survival-neutral decoder binary differs: {unit_id}")
    decode_raw_checks = _validate_survival_decode_raw_inputs(
        decode_contract.get("raw_inputs"),
        unit_id=unit_id,
        unit_dir=unit_dir,
        decode_dir=decode_dir,
        repo_root=repo_root,
    )

    sim_contract_path = unit_dir / "simulation_contract.json"
    decode_contract_path = decode_dir / "decode_contract.json"
    if (
        _read_json(sim_contract_path) != sim_contract
        or _read_json(decode_contract_path) != decode_contract
    ):
        raise ValueError(f"survival-neutral standalone contract differs: {unit_id}")
    checks = [
        FileCheck(
            "survival_unit",
            "single_origin_census_survival",
            unit_id,
            "simulation_completion",
            "simulation_complete",
            sim_path,
            _path_identity(sim_path, repo_root),
            expected_inventory["simulation_completion_sha256"],
            sim_path.stat().st_size,
        ),
        FileCheck(
            "survival_unit",
            "single_origin_census_survival",
            unit_id,
            "decode_completion",
            "decode_complete",
            decode_path,
            _path_identity(decode_path, repo_root),
            expected_inventory["decode_completion_sha256"],
            decode_path.stat().st_size,
        ),
        _binary_file_check(
            sim_contract.get("slim_binary"),
            repo_root=repo_root,
            unit_id=unit_id,
            record_key="slim_binary",
        ),
    ]
    checks.extend(
        _validate_frozen_output_mapping(
            sim.get("outputs"),
            SURVIVAL_SIMULATION_OUTPUT_PATHS,
            scope="survival_unit",
            source="single_origin_census_survival",
            unit_id=unit_id,
            stage="simulation_output",
            base=unit_dir,
            repo_root=repo_root,
            expected_rows={
                "sample_manifest": SAMPLE_DIPLOIDS,
                "truth_unit_scores": len(THRESHOLDS_YEARS),
            },
        )
    )
    checks.extend(
        _validate_frozen_output_mapping(
            decode.get("outputs"),
            SURVIVAL_DECODE_OUTPUT_PATHS,
            scope="survival_unit",
            source="single_origin_census_survival",
            unit_id=unit_id,
            stage="decode_output",
            base=decode_dir,
            repo_root=repo_root,
            expected_rows={
                "overall_pairs": SAMPLE_DIPLOIDS,
                "gamma_unit_scores": len(THRESHOLDS_YEARS),
            },
        )
    )
    checks.extend(decode_raw_checks)
    raw_inputs = decode_contract["raw_inputs"]
    if (
        raw_inputs["simulation_tree"]["sha256"] != sim["outputs"]["tree"]["sha256"]
        or raw_inputs["sample_manifest"]["sha256"]
        != sim["outputs"]["sample_manifest"]["sha256"]
        or raw_inputs["overall_pairs"]["sha256"]
        != decode["outputs"]["overall_pairs"]["sha256"]
    ):
        raise ValueError(
            f"survival-neutral decode raw/output binding differs: {unit_id}"
        )
    unit_scores = pd.concat(
        [
            _read_tsv(unit_dir / "truth_unit_scores.tsv", label="truth unit scores"),
            _read_tsv(decode_dir / "gamma_unit_scores.tsv", label="Gamma unit scores"),
        ],
        ignore_index=True,
    )
    unit_scores = _validate_survival_scores_for_one_unit(unit_scores, unit)
    _validate_survival_panel_artifacts(
        unit_dir,
        decode_dir,
        unit_endpoint,
        focal_expectation,
        sim.get("focal_identity"),
        individual_location_receipt,
    )
    _validate_survival_decoder_run(
        decode_dir / "decoder_run.json",
        unit_id=unit_id,
        raw_inputs=raw_inputs,
        unit_dir=unit_dir,
        decode_dir=decode_dir,
        repo_root=repo_root,
    )
    _verify_survival_focal_reducers(unit_dir, decode_dir, unit_scores)
    return (
        checks,
        {
            **expected_inventory,
            "simulation_contract_sha256": sim_contract_sha,
            "simulation_outputs_sha256": _canonical_sha256(sim["outputs"]),
            "decode_contract_sha256": decode_contract_sha,
            "decode_outputs_sha256": _canonical_sha256(decode["outputs"]),
            "endpoint": unit_endpoint,
            "endpoint_sha256": _canonical_sha256(unit_endpoint),
            "focal_identity_sha256": _canonical_sha256(sim.get("focal_identity")),
            "focal_overlay_patch_sha256": _canonical_sha256(focal_overlay_receipt),
            "individual_location_canonicalization_sha256": _canonical_sha256(
                individual_location_receipt
            ),
            "slim_end_hook_patch_sha256": _canonical_sha256(slim_end_hook_patch),
        },
        unit_scores,
    )


def _validate_survival_scores_for_one_unit(
    frame: pd.DataFrame, unit: Mapping[str, Any]
) -> pd.DataFrame:
    if tuple(frame.columns) != SURVIVAL_SCORE_COLUMNS:
        raise ValueError(
            "survival-neutral unit score columns differ from frozen schema"
        )
    unit_id = str(unit["unit_id"])
    if len(frame) != 2 * len(THRESHOLDS_YEARS):
        raise ValueError(f"survival-neutral unit score row count differs: {unit_id}")
    if set(frame["unit_id"].astype(str)) != {unit_id}:
        raise ValueError(f"survival-neutral unit score ID differs: {unit_id}")
    if set(_exact_integer(frame, "replicate_index")) != {
        int(unit["replicate_index"])
    } or set(_exact_integer(frame, "seed")) != {int(unit["seed"])}:
        raise ValueError(f"survival-neutral unit score key differs: {unit_id}")
    if set(frame["source"].astype(str)) != {"tree_truth", "gamma_smc"}:
        raise ValueError(f"survival-neutral unit score source differs: {unit_id}")
    expected = {
        (source, threshold)
        for source in ("tree_truth", "gamma_smc")
        for threshold in THRESHOLDS_YEARS
    }
    observed = set(
        zip(
            frame["source"].astype(str),
            _exact_integer(frame, "threshold_years"),
            strict=True,
        )
    )
    if observed != expected:
        raise ValueError(f"survival-neutral unit score grid differs: {unit_id}")
    return frame.sort_values(["source", "threshold_years"], kind="stable")


def _singleton_metadata_value(value: Any, *, label: str) -> Any:
    while isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        if len(value) != 1:
            raise ValueError(f"survival-neutral tree metadata is not scalar: {label}")
        value = value[0]
    return value


def _integer_metadata_value(value: Any, *, label: str) -> int:
    raw = _singleton_metadata_value(value, label=label)
    if isinstance(raw, bool):
        raise TypeError(f"survival-neutral tree metadata is not integer: {label}")
    try:
        numeric = float(raw)
    except (TypeError, ValueError) as error:
        raise ValueError(
            f"survival-neutral tree metadata is not integer: {label}"
        ) from error
    if not np.isfinite(numeric) or numeric != np.floor(numeric):
        raise ValueError(f"survival-neutral tree metadata is not integer: {label}")
    return int(numeric)


def _boolean_metadata_value(value: Any, *, label: str) -> bool:
    raw = _singleton_metadata_value(value, label=label)
    if isinstance(raw, (bool, np.bool_)):
        return bool(raw)
    if isinstance(raw, (int, float)) and not isinstance(raw, bool) and raw in (0, 1):
        return bool(raw)
    raise ValueError(f"survival-neutral tree metadata is not boolean: {label}")


def _validate_survival_tree_census(
    tree_sequence: tskit.TreeSequence, endpoint: Mapping[str, Any]
) -> Mapping[str, Any]:
    metadata = tree_sequence.metadata
    if not isinstance(metadata, Mapping):
        raise TypeError("survival-neutral tree lacks SLiM metadata")
    slim = metadata.get("SLiM")
    if not isinstance(slim, Mapping):
        raise TypeError("survival-neutral tree lacks exact SLiM metadata")
    user = slim.get("user_metadata")
    if not isinstance(user, Mapping):
        raise TypeError("survival-neutral tree lacks SLiM user metadata")
    key_map = {
        "final_census_alt_count": "survival_neutral_final_census_alt_count",
        "final_census_total_count": "survival_neutral_final_census_total_count",
        "final_census_af": "survival_neutral_final_census_af",
        "final_census_segregating": "survival_neutral_final_census_segregating",
        "final_census_fixed": "survival_neutral_final_census_fixed",
    }
    for metadata_key in key_map.values():
        occurrences = 0

        def count(candidate: Any, target: str = metadata_key) -> None:
            nonlocal occurrences
            if isinstance(candidate, Mapping):
                occurrences += int(target in candidate)
                for nested in candidate.values():
                    count(nested)
            elif isinstance(candidate, Sequence) and not isinstance(
                candidate, (str, bytes)
            ):
                for nested in candidate:
                    count(nested)

        count(metadata)
        if occurrences != 1 or metadata_key not in user:
            raise ValueError(
                f"survival-neutral census metadata is not uniquely scoped: {metadata_key}"
            )
    observed_alt = _integer_metadata_value(
        user[key_map["final_census_alt_count"]],
        label=key_map["final_census_alt_count"],
    )
    observed_total = _integer_metadata_value(
        user[key_map["final_census_total_count"]],
        label=key_map["final_census_total_count"],
    )
    observed_af_raw = _singleton_metadata_value(
        user[key_map["final_census_af"]], label=key_map["final_census_af"]
    )
    try:
        observed_af = float(observed_af_raw)
    except (TypeError, ValueError) as error:
        raise ValueError("survival-neutral tree census AF is malformed") from error
    if not np.isfinite(observed_af):
        raise ValueError("survival-neutral tree census AF is nonfinite")
    observed_segregating = _boolean_metadata_value(
        user[key_map["final_census_segregating"]],
        label=key_map["final_census_segregating"],
    )
    observed_fixed = _boolean_metadata_value(
        user[key_map["final_census_fixed"]],
        label=key_map["final_census_fixed"],
    )
    if (
        observed_alt != int(endpoint["final_census_alt_count"])
        or observed_total != int(endpoint["final_census_total_count"])
        or not np.isclose(observed_af, float(endpoint["final_census_af"]))
        or observed_segregating is not True
        or observed_fixed is not False
        or bool(endpoint["final_census_segregating"]) is not True
    ):
        raise ValueError(
            "survival-neutral census endpoint differs from exact tree end-hook metadata"
        )
    q_value = _singleton_metadata_value(user.get("Q"), label="Q")
    if not np.isclose(float(q_value), 5.0):
        raise ValueError("survival-neutral tree SLiM scaling factor differs")
    return slim


def _validate_survival_tree_focal_identity(
    tree_sequence: tskit.TreeSequence,
    expectation: Mapping[str, Any],
    declared_identity: Any,
    *,
    pool_alt_count: int,
) -> tuple[tskit.Site | None, str | None]:
    expectation = _validate_survival_focal_expectation(expectation)
    focal_sites = [
        site
        for site in tree_sequence.sites()
        if np.isclose(float(site.position), FOCAL_POSITION_BP, rtol=0.0, atol=1e-9)
    ]
    entries: list[tuple[tskit.Mutation, Mapping[str, Any]]] = []
    malformed_metadata = False
    for site in focal_sites:
        for mutation in site.mutations:
            metadata = mutation.metadata
            if not isinstance(metadata, Mapping):
                malformed_metadata = True
                continue
            mutation_list = metadata.get("mutation_list")
            if not isinstance(mutation_list, list) or not mutation_list:
                malformed_metadata = True
                continue
            for entry in mutation_list:
                if not isinstance(entry, Mapping):
                    malformed_metadata = True
                else:
                    entries.append((mutation, entry))
    mutation_types = sorted(
        {
            int(entry["mutation_type"])
            for _, entry in entries
            if "mutation_type" in entry
        }
    )
    expected_type = int(expectation["expected_mutation_type"])
    selected_entries = [
        (mutation, entry)
        for mutation, entry in entries
        if int(entry.get("mutation_type", -1)) == expected_type
    ]
    base = {
        "schema": SURVIVAL_FOCAL_IDENTITY_SCHEMA,
        "focal_position_0based": FOCAL_POSITION_BP,
        "expected_mutation_type": expected_type,
        "expected_source_population": str(expectation["source_population"]),
        "expected_source_subpopulation": int(expectation["source_subpopulation"]),
        "expected_origin_time_generations": float(
            expectation["realized_origin_time_generations"]
        ),
        "observed_focal_site_count": len(focal_sites),
        "observed_tskit_mutation_count": sum(
            len(site.mutations) for site in focal_sites
        ),
        "observed_metadata_entry_count": len(entries),
        "observed_mutation_types": mutation_types,
    }
    if not selected_entries:
        if malformed_metadata:
            raise ValueError("survival-neutral focal identity metadata is malformed")
        status = (
            "selected_focal_mutation_absent"
            if not focal_sites
            else "selected_focal_mutation_absent_background_only"
        )
        observed = {
            **base,
            "status": status,
            "validation": (
                "causal_site_absent_from_observational_pool; exact census survival "
                "retained in SLiM metadata"
            ),
        }
        if pool_alt_count != 0 or status != "selected_focal_mutation_absent":
            raise ValueError(
                "survival-neutral focal mutation is absent despite pool detection"
            )
        if declared_identity != observed:
            raise ValueError("survival-neutral focal identity receipt differs")
        return None, None
    if (
        len(focal_sites) != 1
        or len(focal_sites[0].mutations) != 1
        or len(entries) != 1
        or len(selected_entries) != 1
    ):
        raise ValueError("survival-neutral focal identity is recurrent or ambiguous")
    site = focal_sites[0]
    mutation, entry = selected_entries[0]
    required = {
        "mutation_type",
        "selection_coeff",
        "subpopulation",
        "slim_time",
        "nucleotide",
    }
    if not required.issubset(entry) or mutation.parent != tskit.NULL:
        raise ValueError("survival-neutral focal identity metadata fields differ")
    expected_source = int(expectation["source_subpopulation"])
    expected_time = float(expectation["realized_origin_time_generations"])
    observed_source = int(entry["subpopulation"])
    observed_time = float(mutation.time)
    observed_coefficient = float(entry["selection_coeff"])
    nucleotide = int(entry["nucleotide"])
    if (
        observed_source != expected_source
        or not np.isclose(observed_time, expected_time, atol=1e-9)
        or not np.isclose(observed_coefficient, 0.0, atol=1e-12)
        or nucleotide not in range(4)
        or mutation.derived_state != "ACGT"[nucleotide]
        or site.ancestral_state == mutation.derived_state
    ):
        raise ValueError("survival-neutral focal identity values differ")
    try:
        slim = tree_sequence.metadata["SLiM"]
        cycle = int(slim["cycle"])
        q = float(_singleton_metadata_value(slim["user_metadata"]["Q"], label="Q"))
    except (KeyError, TypeError, ValueError, IndexError) as error:
        raise ValueError("survival-neutral focal identity lacks cycle/Q") from error
    slim_time = int(entry["slim_time"])
    if not np.isclose(q, 5.0) or slim_time != round(cycle - observed_time / q):
        raise ValueError("survival-neutral focal identity time scaling differs")
    observed = {
        **base,
        "status": "valid",
        "observed_mutation_type": int(entry["mutation_type"]),
        "observed_source_subpopulation": observed_source,
        "observed_origin_time_generations": observed_time,
        "observed_selection_coeff": observed_coefficient,
        "observed_slim_time": slim_time,
        "observed_slim_cycle": cycle,
        "observed_slim_scaling_factor": q,
        "observed_ancestral_state": site.ancestral_state,
        "observed_derived_state": mutation.derived_state,
        "observed_nucleotide": nucleotide,
    }
    if declared_identity != observed:
        raise ValueError("survival-neutral focal identity receipt differs")
    return site, str(mutation.derived_state)


def _validate_survival_panel_tree(
    tree_path: Path,
    manifest: pd.DataFrame,
    endpoint: Mapping[str, Any],
    focal_expectation: Mapping[str, Any],
    focal_identity: Any,
    individual_location_receipt: Mapping[str, Any],
) -> None:
    try:
        tree_sequence = tskit.load(tree_path)
    except (OSError, tskit.FileFormatError) as error:
        raise ValueError("survival-neutral panel tree is unreadable") from error
    _validate_survival_individual_location_receipt(
        individual_location_receipt, tree=tree_sequence
    )
    _validate_survival_tree_census(tree_sequence, endpoint)
    focal_site, focal_derived_state = _validate_survival_tree_focal_identity(
        tree_sequence,
        focal_expectation,
        focal_identity,
        pool_alt_count=int(endpoint["pool_alt_count"]),
    )
    han_population_ids = []
    for population in tree_sequence.populations():
        metadata = population.metadata
        if isinstance(metadata, Mapping) and (
            metadata.get("name") == "Han" or metadata.get("id") == "Han"
        ):
            han_population_ids.append(int(population.id))
    if len(han_population_ids) != 1:
        raise ValueError("survival-neutral panel tree has no unique Han population")
    han_population_id = han_population_ids[0]
    neanderthal_population_ids = [
        int(population.id)
        for population in tree_sequence.populations()
        if isinstance(population.metadata, Mapping)
        and (
            population.metadata.get("name") == "Neanderthal"
            or population.metadata.get("id") == "Neanderthal"
        )
    ]
    if len(neanderthal_population_ids) != 1:
        raise ValueError(
            "survival-neutral panel tree has no unique Neanderthal population"
        )
    if neanderthal_population_ids[0] != int(focal_expectation["source_subpopulation"]):
        raise ValueError("survival-neutral Neanderthal population ID differs")
    pool = []
    for individual in tree_sequence.individuals():
        nodes = tuple(int(node_id) for node_id in individual.nodes)
        if len(nodes) != 2:
            continue
        node_rows = [tree_sequence.node(node_id) for node_id in nodes]
        if all(
            np.isclose(float(node.time), 0.0)
            and int(node.population) == han_population_id
            for node in node_rows
        ):
            pool.append((int(individual.id), nodes))
    if len(pool) != CANDIDATE_POOL_DIPLOIDS:
        raise ValueError("survival-neutral panel tree lacks the exact 500-diploid pool")
    chosen_rows = endpoint["chosen_panel_rows_zero_based"]
    expected = [pool[int(row)] for row in chosen_rows]
    expected_individuals = np.asarray([individual_id for individual_id, _ in expected])
    expected_nodes_0 = np.asarray([nodes[0] for _, nodes in expected])
    expected_nodes_1 = np.asarray([nodes[1] for _, nodes in expected])
    if (
        not np.array_equal(
            _exact_integer(manifest, "tree_sequence_individual_id"),
            expected_individuals,
        )
        or not np.array_equal(
            _exact_integer(manifest, "sample_node_0"), expected_nodes_0
        )
        or not np.array_equal(
            _exact_integer(manifest, "sample_node_1"), expected_nodes_1
        )
    ):
        raise ValueError(
            "survival-neutral manifest does not match the seeded 100-of-500 panel"
        )
    expected_samples = {node_id for _, nodes in expected for node_id in nodes}
    observed_samples = {int(node_id) for node_id in tree_sequence.samples()}
    if observed_samples != expected_samples:
        raise ValueError(
            "survival-neutral panel tree sample flags differ from manifest"
        )

    if focal_site is None:
        pool_counts = np.zeros(CANDIDATE_POOL_DIPLOIDS, dtype=np.int8)
    else:
        ordered_pool_nodes = [node_id for _, nodes in pool for node_id in nodes]
        variant = tskit.Variant(tree_sequence, samples=ordered_pool_nodes)
        variant.decode(focal_site.id)
        causal_indices = [
            index
            for index, allele in enumerate(variant.alleles)
            if allele == focal_derived_state
        ]
        if len(causal_indices) != 1 or causal_indices[0] == 0:
            raise ValueError("survival-neutral focal causal allele is ambiguous")
        genotypes = np.asarray(variant.genotypes, dtype=np.int64)
        if np.any(genotypes < 0) or not set(genotypes).issubset({0, causal_indices[0]}):
            raise ValueError("survival-neutral focal tree genotypes differ")
        pool_counts = (
            genotypes.reshape(CANDIDATE_POOL_DIPLOIDS, 2) == causal_indices[0]
        ).sum(axis=1)
    if int(pool_counts.sum()) != int(endpoint["pool_alt_count"]):
        raise ValueError("survival-neutral tree pool AF differs from endpoint")
    panel_counts = np.asarray(pool_counts)[np.asarray(chosen_rows, dtype=int)]
    if not np.array_equal(
        panel_counts,
        _exact_integer(manifest, "focal_selected_allele_count"),
    ) or int(panel_counts.sum()) != int(endpoint["sample_alt_count"]):
        raise ValueError(
            "survival-neutral tree panel genotypes differ from manifest or endpoint"
        )


def _validate_survival_panel_artifacts(
    unit_dir: Path,
    decode_dir: Path,
    endpoint: Mapping[str, Any],
    focal_expectation: Mapping[str, Any],
    focal_identity: Any,
    individual_location_receipt: Mapping[str, Any],
) -> None:
    manifest = _read_tsv(
        unit_dir / "sample_manifest.tsv", label="survival sample manifest"
    )
    if (
        tuple(manifest.columns) != SURVIVAL_SAMPLE_MANIFEST_COLUMNS
        or len(manifest) != SAMPLE_DIPLOIDS
    ):
        raise ValueError("survival-neutral sample manifest schema differs")
    vcf_index = _exact_integer(manifest, "vcf_diploid_index")
    if not np.array_equal(vcf_index, np.arange(SAMPLE_DIPLOIDS)):
        raise ValueError("survival-neutral manifest diploid order differs")
    haplotype_0 = _exact_integer(manifest, "gamma_smc_haplotype_0")
    haplotype_1 = _exact_integer(manifest, "gamma_smc_haplotype_1")
    if not np.array_equal(haplotype_0, 2 * vcf_index) or not np.array_equal(
        haplotype_1, 2 * vcf_index + 1
    ):
        raise ValueError("survival-neutral manifest Gamma haplotypes differ")
    if not np.allclose(
        _exact_numeric(manifest, "focal_position_0based"), FOCAL_POSITION_BP
    ) or set(_exact_integer(manifest, "focal_selected_allele_index")) != {1}:
        raise ValueError("survival-neutral manifest focal marker differs")
    allele_counts = _exact_integer(manifest, "focal_selected_allele_count")
    if not np.isin(allele_counts, (0, 1, 2)).all():
        raise ValueError("survival-neutral manifest focal genotypes differ")
    expected_classes = np.asarray(
        [{0: "hom_ref", 1: "heterozygous", 2: "hom_alt"}[int(x)] for x in allele_counts]
    )
    if not np.array_equal(manifest["genotype_class"].astype(str), expected_classes):
        raise ValueError("survival-neutral manifest genotype classes differ")
    for column in ("tree_sequence_individual_id", "sample_node_0", "sample_node_1"):
        values = _exact_integer(manifest, column)
        if np.any(values < 0) or len(set(values)) != len(values):
            raise ValueError(f"survival-neutral manifest identifiers differ: {column}")
    all_nodes = np.concatenate(
        [
            _exact_integer(manifest, "sample_node_0"),
            _exact_integer(manifest, "sample_node_1"),
        ]
    )
    if len(set(all_nodes)) != 2 * SAMPLE_DIPLOIDS:
        raise ValueError("survival-neutral manifest sample nodes are not unique")
    if int(allele_counts.sum()) != int(endpoint["sample_alt_count"]):
        raise ValueError("survival-neutral manifest sample AF differs from endpoint")

    _validate_survival_panel_tree(
        unit_dir / "simulation.trees",
        manifest,
        endpoint,
        focal_expectation,
        focal_identity,
        individual_location_receipt,
    )

    pair_path = decode_dir / "overall.pairs.tsv"
    pair_lines = pair_path.read_text(encoding="utf-8").splitlines()
    expected_header = [
        "# Gamma-SMC explicit haplotype pairs",
        "# genotype_class\toverall",
        f"# n_pairs\t{SAMPLE_DIPLOIDS}",
    ]
    expected_data = [
        f"{left}\t{right}" for left, right in zip(haplotype_0, haplotype_1, strict=True)
    ]
    if pair_lines != [*expected_header, *expected_data]:
        raise ValueError(
            "survival-neutral explicit pair file differs from sample manifest"
        )


def _verify_survival_focal_reducers(
    unit_dir: Path, decode_dir: Path, unit_scores: pd.DataFrame
) -> None:
    truth = _read_tsv(
        unit_dir / "truth_profiles.tsv.gz", label="survival truth profiles"
    )
    truth_required = {
        "source",
        "position_0based",
        "genotype_class",
        "n_pairs",
        "threshold_years",
        "p_tmrca_lt_threshold",
    }
    if not truth_required.issubset(truth.columns):
        raise ValueError("survival truth profiles lack focal reducer columns")
    truth_focal = truth[
        (truth["source"].astype(str) == "tree_truth")
        & (truth["genotype_class"].astype(str) == "overall")
        & np.isclose(pd.to_numeric(truth["position_0based"]), FOCAL_POSITION_BP)
    ].sort_values("threshold_years", kind="stable")
    if len(truth_focal) != len(THRESHOLDS_YEARS) or not np.array_equal(
        _exact_integer(truth_focal, "threshold_years"),
        np.asarray(THRESHOLDS_YEARS),
    ):
        raise ValueError("survival truth focal threshold grid differs")
    if set(_exact_integer(truth_focal, "n_pairs")) != {SAMPLE_DIPLOIDS}:
        raise ValueError("survival truth focal pair count differs")
    truth_scores = unit_scores[unit_scores["source"] == "tree_truth"].sort_values(
        "threshold_years", kind="stable"
    )
    if not np.allclose(
        _exact_numeric(truth_scores, "overall_focal_p_tmrca_lt_threshold"),
        _exact_numeric(truth_focal, "p_tmrca_lt_threshold"),
        rtol=1e-12,
        atol=1e-15,
    ):
        raise ValueError("survival truth unit scores differ from focal truth profiles")

    gamma = _read_tsv(
        decode_dir / "overall.summary.tsv", label="survival Gamma overall summary"
    )
    gamma_required = {
        "position_0based",
        "n_pairs",
        *(f"mean_p_lt_{threshold}" for threshold in THRESHOLDS_YEARS),
    }
    if not gamma_required.issubset(gamma.columns):
        raise ValueError("survival Gamma summary lacks focal reducer columns")
    gamma_focal = gamma[
        np.isclose(pd.to_numeric(gamma["position_0based"]), FOCAL_POSITION_BP)
    ]
    if len(gamma_focal) != 1 or int(gamma_focal.iloc[0]["n_pairs"]) != SAMPLE_DIPLOIDS:
        raise ValueError("survival Gamma focal row differs")
    gamma_scores = unit_scores[unit_scores["source"] == "gamma_smc"].sort_values(
        "threshold_years", kind="stable"
    )
    gamma_values = np.asarray(
        [
            float(gamma_focal.iloc[0][f"mean_p_lt_{threshold}"])
            for threshold in THRESHOLDS_YEARS
        ]
    )
    if not np.allclose(
        _exact_numeric(gamma_scores, "overall_focal_p_tmrca_lt_threshold"),
        gamma_values,
        rtol=1e-12,
        atol=1e-15,
    ):
        raise ValueError("survival Gamma unit scores differ from focal decoder row")


def _load_survival_bundle(
    repo_root: Path, results_dir: Path, *, workers: int
) -> SurvivalBundle:
    study_root = results_dir.parent
    plan_contract, plan_contract_sha, execution, plan_checks, source_audit = (
        _load_survival_plan(repo_root, study_root)
    )
    execution_sha = _sha256_file(study_root / "plan" / "execution_units.tsv")
    completion_path = results_dir / "RESULTS_COMPLETION.json"
    _require_regular_file(
        completion_path, boundary=results_dir, label="survival-neutral completion"
    )
    completion = _read_json(completion_path)
    _require_exact_receipt_fields(
        completion,
        {
            "schema",
            "status",
            "created_utc",
            "contract",
            "contract_sha256",
            "outputs",
            "unit_inventory",
        },
        label="survival-neutral results",
    )
    contract, contract_sha = _validate_completion_contract(
        completion,
        expected_schema=SURVIVAL_RESULTS_SCHEMA_VERSION,
        label="survival-neutral results",
    )
    _validate_survival_root_contract(
        contract,
        plan_contract_sha256=plan_contract_sha,
        execution_sha256=execution_sha,
    )
    if contract.get("source_hashes") != plan_contract.get("source_hashes"):
        raise ValueError("survival-neutral results source hashes differ from plan")
    if contract.get("decoder_binary") != plan_contract.get("decoder_binary"):
        raise ValueError("survival-neutral results decoder binary differs from plan")
    outputs = completion.get("outputs")
    checks = [
        *plan_checks,
        FileCheck(
            "survival_neutral",
            "single_origin_census_survival",
            "",
            "results_completion",
            "RESULTS_COMPLETION",
            completion_path,
            _path_identity(completion_path, repo_root),
            None,
            None,
        ),
    ]
    checks.extend(
        _validate_frozen_output_mapping(
            outputs,
            SURVIVAL_RESULTS_OUTPUT_PATHS,
            scope="survival_neutral",
            source="single_origin_census_survival",
            unit_id="",
            stage="results_output",
            base=results_dir,
            repo_root=repo_root,
            expected_rows={
                "endpoint_summary": SURVIVAL_NEUTRAL_REPLICATES,
                "replicate_scores": 2
                * SURVIVAL_NEUTRAL_REPLICATES
                * len(THRESHOLDS_YEARS),
                "simulation_status": SURVIVAL_NEUTRAL_REPLICATES,
                "decode_status": SURVIVAL_NEUTRAL_REPLICATES,
            },
        )
    )
    declared_paths = set(SURVIVAL_RESULTS_OUTPUT_PATHS.values())
    actual_paths = {
        path.relative_to(results_dir).as_posix()
        for path in results_dir.rglob("*")
        if path.is_file() and path.name != "RESULTS_COMPLETION.json"
    }
    if actual_paths != declared_paths:
        raise ValueError("survival-neutral result inventory differs from completion")
    unit_inventory = completion.get("unit_inventory")
    if not isinstance(unit_inventory, list) or len(unit_inventory) != (
        SURVIVAL_NEUTRAL_REPLICATES
    ):
        raise ValueError("survival-neutral root unit inventory is incomplete")
    inventory_by_id = {}
    frozen_inventory_keys = {
        "unit_id",
        "replicate_index",
        "seed",
        "simulation_completion_path",
        "simulation_completion_sha256",
        "decode_completion_path",
        "decode_completion_sha256",
    }
    for record in unit_inventory:
        if not isinstance(record, Mapping) or set(record) != frozen_inventory_keys:
            raise ValueError("survival-neutral root unit record differs")
        unit_id = str(record["unit_id"])
        if unit_id in inventory_by_id:
            raise ValueError("survival-neutral root unit inventory has duplicates")
        inventory_by_id[unit_id] = record
    if set(inventory_by_id) != set(_expected_survival_unit_ids()):
        raise ValueError("survival-neutral root unit inventory IDs differ")

    unit_bindings = []
    per_unit_scores = []
    for unit in execution.to_dict(orient="records"):
        unit_checks, binding, unit_score = _survival_unit_receipts(
            repo_root,
            study_root,
            unit,
            inventory_by_id[str(unit["unit_id"])],
            plan_contract=plan_contract,
            plan_contract_sha256=plan_contract_sha,
        )
        checks.extend(unit_checks)
        unit_bindings.append(binding)
        per_unit_scores.append(unit_score)

    scores_path = results_dir / "replicate_scores.tsv"
    endpoints_path = results_dir / "endpoint_summary.tsv"
    scores = _validate_survival_scores(
        _read_tsv(scores_path, label="survival-neutral replicate scores")
    )
    endpoints_input = _validate_survival_endpoints(
        _read_tsv(endpoints_path, label="survival-neutral endpoints")
    )
    _verify_survival_score_endpoints(scores, endpoints_input)
    deep_scores = _validate_survival_scores(
        pd.concat(per_unit_scores, ignore_index=True)
    )
    try:
        pd.testing.assert_frame_equal(
            scores,
            deep_scores,
            check_dtype=False,
            check_exact=False,
            rtol=1e-12,
            atol=1e-15,
        )
    except AssertionError as error:
        raise ValueError(
            "survival-neutral aggregate scores differ from per-unit source outputs"
        ) from error
    completion_sha_by_id = {
        str(binding["unit_id"]): binding["simulation_completion_sha256"]
        for binding in unit_bindings
    }
    if any(
        completion_sha_by_id[str(row["unit_id"])]
        != str(row["simulation_completion_sha256"])
        for row in endpoints_input.to_dict(orient="records")
    ):
        raise ValueError("survival-neutral endpoint simulation receipt hash differs")
    endpoint_by_id = {
        str(binding["unit_id"]): binding["endpoint"] for binding in unit_bindings
    }
    for row in endpoints_input.to_dict(orient="records"):
        unit_endpoint = endpoint_by_id[str(row["unit_id"])]
        for key in SURVIVAL_ENDPOINT_COLUMNS[4:]:
            observed = row[key]
            expected = unit_endpoint[key]
            if isinstance(expected, float):
                equal = np.isclose(float(observed), expected)
            else:
                equal = observed == expected
            if not equal:
                raise ValueError(
                    f"survival-neutral aggregate endpoint differs from unit receipt: {key}"
                )

    simulation_status = _read_tsv(
        results_dir / "simulation_status.tsv", label="survival simulation status"
    )
    decode_status = _read_tsv(
        results_dir / "decode_status.tsv", label="survival decode status"
    )
    expected_simulation_status_columns = (
        *SURVIVAL_EXECUTION_COLUMNS,
        "status",
        "error",
        "final_census_af",
        "simulation_completion_sha256",
    )
    expected_decode_status_columns = (
        *SURVIVAL_EXECUTION_COLUMNS,
        "status",
        "error",
        "decode_completion_sha256",
    )
    if (
        tuple(simulation_status.columns) != expected_simulation_status_columns
        or tuple(decode_status.columns) != expected_decode_status_columns
    ):
        raise ValueError("survival-neutral status columns differ from frozen schema")
    simulation_execution = _validate_survival_execution(
        simulation_status[list(SURVIVAL_EXECUTION_COLUMNS)]
    )
    decode_execution = _validate_survival_execution(
        decode_status[list(SURVIVAL_EXECUTION_COLUMNS)]
    )
    try:
        pd.testing.assert_frame_equal(
            execution, simulation_execution, check_dtype=False
        )
        pd.testing.assert_frame_equal(execution, decode_execution, check_dtype=False)
    except AssertionError as error:
        raise ValueError("survival-neutral status execution fields differ") from error
    for status, stage in (
        (simulation_status, "simulation"),
        (decode_status, "decode"),
    ):
        if set(status["status"].astype(str)) != {"complete"}:
            raise ValueError(f"survival-neutral final {stage} status is not complete")
        if not (status["error"].fillna("").astype(str) == "").all():
            raise ValueError(f"survival-neutral final {stage} status has an error")
    endpoint_af_by_id = dict(
        zip(
            endpoints_input["unit_id"].astype(str),
            endpoints_input["final_census_af"].to_numpy(dtype=float),
            strict=True,
        )
    )
    binding_by_id = {str(binding["unit_id"]): binding for binding in unit_bindings}
    for row in simulation_status.to_dict(orient="records"):
        unit_id = str(row["unit_id"])
        if not np.isclose(float(row["final_census_af"]), endpoint_af_by_id[unit_id]):
            raise ValueError("survival-neutral simulation status census AF differs")
        if (
            row["simulation_completion_sha256"]
            != binding_by_id[unit_id]["simulation_completion_sha256"]
        ):
            raise ValueError("survival-neutral simulation status receipt hash differs")
    for row in decode_status.to_dict(orient="records"):
        unit_id = str(row["unit_id"])
        if (
            row["decode_completion_sha256"]
            != binding_by_id[unit_id]["decode_completion_sha256"]
        ):
            raise ValueError("survival-neutral decode status receipt hash differs")

    audit = _audit_file_checks(checks, workers=workers)
    audit = pd.concat([audit, source_audit], ignore_index=True).sort_values(
        ["scope", "source", "unit_id", "stage", "record_key", "path"],
        kind="stable",
    )
    endpoints = endpoints_input.rename(
        columns={
            "pool_alt_count": "final_pool_alt_count",
            "pool_total_count": "final_pool_total_count",
            "pool_af": "final_pool_af",
        }
    ).copy()
    endpoints.insert(1, "bank_role", "survival_only_primary")
    endpoints.insert(2, "simulation_class", "neutral")
    endpoints.insert(4, "selection_coefficient", 0.0)
    endpoints.insert(
        5,
        "present_conditioning",
        (
            "strict_complete_simulated_q5_present_han_census_segregation; "
            "pool_and_panel_observational"
        ),
    )
    endpoints["hom_ref_diploids"] = np.nan
    endpoints["heterozygous_diploids"] = np.nan
    endpoints["hom_alt_diploids"] = np.nan
    endpoints["mutation_origin_generations"] = MUTATION_ORIGIN_GENERATIONS
    endpoints["introgression_pulse_generations"] = INTROGRESSION_PULSE_GENERATIONS
    endpoints["selection_end_generations"] = np.nan
    endpoints["realized_frequency_trajectory_persisted"] = False
    return SurvivalBundle(
        scores=scores,
        endpoints=endpoints.sort_values("replicate_index", kind="stable"),
        completion=completion,
        completion_path=completion_path,
        audit=audit,
        input_contract={
            "completion": _file_record(completion_path, repo_root),
            "contract_sha256": contract_sha,
            "study_plan": _file_record(
                study_root / "plan" / "study_plan.json", repo_root
            ),
            "plan_contract_sha256": plan_contract_sha,
            "execution_units": _file_record(
                study_root / "plan" / "execution_units.tsv", repo_root
            ),
            "declared_output_count": len(outputs),
            "declared_outputs_sha256": _canonical_sha256(outputs),
            "replicate_scores": _file_record(scores_path, repo_root),
            "endpoint_summary": _file_record(endpoints_path, repo_root),
            "unit_count": SURVIVAL_NEUTRAL_REPLICATES,
            "unit_receipt_inventory_sha256": _canonical_sha256(unit_bindings),
            "unit_ids": list(_expected_survival_unit_ids()),
            "conditioning": (
                "strict complete simulated Q=5 present Han census segregation "
                "(unscaled Ne=6300; scaled census=1260 diploids/2520 genomes); no "
                "terminal AF target/band, pool-frequency, selected-AF matching, or "
                "sample-detection gate; uniform pool and panel endpoints are observational"
            ),
        },
    )


def _existing_unit_scores(
    bundle: AnalysisBundle, existing_units: pd.DataFrame
) -> pd.DataFrame:
    frame = bundle.tables["replicate_statistics"].copy()
    required = {
        "unit_id",
        "demography_id",
        "simulation_class",
        "selection_coefficient",
        "target_allele_frequency",
        "replicate_index",
        "panel_n_diploid_samples",
        "analysis_source",
        "statistic_scope",
        "metric",
        "statistic",
        "threshold_years",
        "selection_score",
    }
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"{bundle.spec.key} replicate statistics lack columns")
    expected_ids = set(existing_units["unit_id"].astype(str))
    stats = ("p_overall", "p_hom_ref", "p_hom_alt", "p_hom_alt_minus_hom_ref")
    frame = frame[
        (frame["unit_id"].astype(str).isin(expected_ids))
        & (frame["demography_id"].astype(str) == DEMOGRAPHY_ID)
        & np.isclose(pd.to_numeric(frame["target_allele_frequency"]), TARGET_AF)
        & (frame["analysis_source"].astype(str) == bundle.spec.analysis_source)
        & (frame["statistic_scope"].astype(str) == "focal_nearest")
        & (frame["metric"].astype(str) == "p_tmrca_lt_threshold")
        & (frame["statistic"].astype(str).isin(stats))
        & pd.to_numeric(frame["threshold_years"]).isin(THRESHOLDS_YEARS)
    ].copy()
    expected_rows = (
        (SELECTED_REPLICATES + MATCHED_AF_NEUTRAL_REPLICATES)
        * len(stats)
        * len(THRESHOLDS_YEARS)
    )
    if len(frame) != expected_rows:
        raise ValueError(f"{bundle.spec.key} focal score grid is incomplete")
    keys = [
        "unit_id",
        "simulation_class",
        "selection_coefficient",
        "target_allele_frequency",
        "replicate_index",
        "panel_n_diploid_samples",
        "threshold_years",
    ]
    if frame.duplicated([*keys, "statistic"]).any():
        raise ValueError(f"{bundle.spec.key} focal scores contain duplicates")
    wide = (
        frame.pivot(index=keys, columns="statistic", values="selection_score")
        .reset_index()
        .rename_axis(columns=None)
    )
    if set(wide["unit_id"].astype(str)) != expected_ids:
        raise ValueError(f"{bundle.spec.key} focal score units are incomplete")
    for column in ("p_overall", "p_hom_ref", "p_hom_alt"):
        values = _exact_numeric(wide, column)
        if np.any((values < 0) | (values > 1)):
            raise ValueError(f"{bundle.spec.key} {column} lies outside [0, 1]")
    delta = _exact_numeric(wide, "p_hom_alt_minus_hom_ref")
    if not np.allclose(
        delta,
        wide["p_hom_alt"].to_numpy(dtype=float)
        - wide["p_hom_ref"].to_numpy(dtype=float),
        atol=1e-10,
        rtol=1e-9,
    ):
        raise ValueError(f"{bundle.spec.key} alt/ref contrast is inconsistent")
    for _, curve in wide.groupby("unit_id", sort=False):
        ordered = curve.sort_values("threshold_years")
        for column in ("p_overall", "p_hom_ref", "p_hom_alt"):
            if np.any(np.diff(ordered[column].to_numpy(dtype=float)) < -1e-10):
                raise ValueError(f"{bundle.spec.key} focal CDF is not nondecreasing")
    wide.insert(0, "source", bundle.spec.key)
    wide["bank_role"] = np.where(
        wide["simulation_class"].astype(str) == "selected",
        "selected",
        "matched_af_sensitivity",
    )
    wide["score_scope"] = "focal_nearest_within_diploid"
    return wide.sort_values(
        ["source", "bank_role", "unit_id", "threshold_years"], kind="stable"
    ).reset_index(drop=True)


def _survival_unit_scores(bundle: SurvivalBundle) -> pd.DataFrame:
    result = bundle.scores[
        [
            "unit_id",
            "replicate_index",
            "source",
            "threshold_years",
            "overall_focal_p_tmrca_lt_threshold",
        ]
    ].rename(columns={"overall_focal_p_tmrca_lt_threshold": "p_overall"})
    result["simulation_class"] = "neutral"
    result["selection_coefficient"] = 0.0
    result["target_allele_frequency"] = np.nan
    result["panel_n_diploid_samples"] = SAMPLE_DIPLOIDS
    result["p_hom_ref"] = np.nan
    result["p_hom_alt"] = np.nan
    result["p_hom_alt_minus_hom_ref"] = np.nan
    result["bank_role"] = "survival_only_primary"
    result["score_scope"] = "focal_nearest_within_diploid"
    order = [
        "source",
        "unit_id",
        "simulation_class",
        "selection_coefficient",
        "target_allele_frequency",
        "replicate_index",
        "panel_n_diploid_samples",
        "threshold_years",
        "p_overall",
        "p_hom_ref",
        "p_hom_alt",
        "p_hom_alt_minus_hom_ref",
        "bank_role",
        "score_scope",
    ]
    return result[order].sort_values(
        ["source", "unit_id", "threshold_years"], kind="stable"
    )


def _upper_reference_pvalue(observed: float, null: np.ndarray) -> tuple[int, float]:
    if len(null) != SURVIVAL_NEUTRAL_REPLICATES:
        raise ValueError("reference p-value requires exactly 100 neutral units")
    exceedances = int(np.count_nonzero(null >= observed))
    return exceedances, float((exceedances + 1) / (len(null) + 1))


def _bh_fdr(values: Sequence[float]) -> np.ndarray:
    pvalues = np.asarray(values, dtype=float)
    if np.any(~np.isfinite(pvalues)) or np.any((pvalues < 0) | (pvalues > 1)):
        raise ValueError("BH inputs must be finite probabilities")
    order = np.argsort(pvalues, kind="stable")
    ranked = pvalues[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.clip(adjusted, 0, 1)
    return result


def _primary_pointwise(unit_scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for source in ("tree_truth", "gamma_smc"):
        source_scores = unit_scores[unit_scores["source"] == source]
        for threshold in THRESHOLDS_YEARS:
            null = source_scores[
                (source_scores["bank_role"] == "survival_only_primary")
                & (source_scores["threshold_years"] == threshold)
            ]["p_overall"].to_numpy(dtype=float)
            selected = source_scores[
                (source_scores["bank_role"] == "selected")
                & (source_scores["threshold_years"] == threshold)
            ]
            if len(null) != SURVIVAL_NEUTRAL_REPLICATES or len(selected) != 10:
                raise ValueError("primary pointwise score cell is incomplete")
            for record in selected.to_dict(orient="records"):
                observed = float(record["p_overall"])
                exceedances, pvalue = _upper_reference_pvalue(observed, null)
                rows.append(
                    {
                        "source": source,
                        "estimand": PRIMARY_ESTIMAND,
                        "unit_id": record["unit_id"],
                        "replicate_index": int(record["replicate_index"]),
                        "threshold_years": threshold,
                        "statistic": "p_overall",
                        "observed_selection_score": observed,
                        "n_null": len(null),
                        "null_mean": float(np.mean(null)),
                        "null_median": float(np.median(null)),
                        "null_sd": float(np.std(null, ddof=1)),
                        "null_q025": float(np.quantile(null, 0.025)),
                        "null_q05": float(np.quantile(null, 0.05)),
                        "null_q95": float(np.quantile(null, 0.95)),
                        "null_q975": float(np.quantile(null, 0.975)),
                        "null_upper_exceedance_count": exceedances,
                        "mc_p_upper": pvalue,
                        "neutral_bank_scope": (
                            "single_origin_present_census_survival_without_terminal_af_matching"
                        ),
                        "inference_scope": (
                            "simulation_reference_tail_combines_selection_and_endpoint_af_enrichment"
                        ),
                    }
                )
    result = pd.DataFrame(rows)
    result["bh_q_upper"] = np.nan
    for indices in result.groupby(["source", "unit_id"], sort=False).groups.values():
        result.loc[indices, "bh_q_upper"] = _bh_fdr(
            result.loc[indices, "mc_p_upper"].to_numpy(dtype=float)
        )
    return result


def _sensitivity_pointwise(bundle: AnalysisBundle) -> pd.DataFrame:
    frame = bundle.tables["selected_pvalues"].copy()
    frame = frame[
        (frame["demography_id"].astype(str) == DEMOGRAPHY_ID)
        & np.isclose(pd.to_numeric(frame["target_allele_frequency"]), TARGET_AF)
        & np.isclose(
            pd.to_numeric(frame["selection_coefficient"]), SELECTION_COEFFICIENT
        )
        & (frame["statistic_scope"].astype(str) == "focal_nearest")
        & (frame["metric"].astype(str) == "p_tmrca_lt_threshold")
        & (frame["statistic"].astype(str) == "p_hom_alt_minus_hom_ref")
        & pd.to_numeric(frame["threshold_years"]).isin(THRESHOLDS_YEARS)
    ].copy()
    if len(frame) != SELECTED_REPLICATES * len(THRESHOLDS_YEARS):
        raise ValueError(f"{bundle.spec.key} matched-AF sensitivity grid is incomplete")
    n_null = _exact_integer(frame, "n_null")
    exceed = _exact_integer(frame, "null_upper_exceedance_count")
    if set(n_null) != {MATCHED_AF_NEUTRAL_REPLICATES}:
        raise ValueError(f"{bundle.spec.key} matched-AF null size differs")
    expected_p = (exceed + 1) / (n_null + 1)
    if not np.allclose(expected_p, _exact_numeric(frame, "mc_p_upper")):
        raise ValueError(f"{bundle.spec.key} matched-AF +1 p-values differ")
    columns = [
        "unit_id",
        "replicate_index",
        "threshold_years",
        "statistic",
        "observed_selection_score",
        "n_null",
        "null_mean",
        "null_median",
        "null_sd",
        "null_q025",
        "null_q05",
        "null_q95",
        "null_q975",
        "null_upper_exceedance_count",
        "mc_p_upper",
        "bh_q_upper",
        "neutral_bank_scope",
    ]
    result = frame[columns].copy()
    result.insert(0, "estimand", SENSITIVITY_ESTIMAND)
    result.insert(0, "source", bundle.spec.key)
    result["inference_scope"] = (
        "selection_isolating_final_af_matched_genotype_contrast_sensitivity"
    )
    return result


def _sensitivity_pointwise_from_scores(unit_scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for source in ("tree_truth", "gamma_smc"):
        source_scores = unit_scores[unit_scores["source"] == source]
        for threshold in THRESHOLDS_YEARS:
            neutral = source_scores[
                (source_scores["bank_role"] == "matched_af_sensitivity")
                & (source_scores["threshold_years"] == threshold)
            ]["p_hom_alt_minus_hom_ref"].to_numpy(dtype=float)
            selected = source_scores[
                (source_scores["bank_role"] == "selected")
                & (source_scores["threshold_years"] == threshold)
            ]
            if len(neutral) != 100 or len(selected) != 10:
                raise ValueError("matched-AF direct score cell is incomplete")
            for record in selected.to_dict(orient="records"):
                observed = float(record["p_hom_alt_minus_hom_ref"])
                exceedances, pvalue = _upper_reference_pvalue(observed, neutral)
                rows.append(
                    {
                        "source": source,
                        "estimand": SENSITIVITY_ESTIMAND,
                        "unit_id": record["unit_id"],
                        "replicate_index": int(record["replicate_index"]),
                        "threshold_years": threshold,
                        "statistic": "p_hom_alt_minus_hom_ref",
                        "observed_selection_score": observed,
                        "n_null": 100,
                        "null_mean": float(np.mean(neutral)),
                        "null_median": float(np.median(neutral)),
                        "null_sd": float(np.std(neutral, ddof=1)),
                        "null_q025": float(np.quantile(neutral, 0.025)),
                        "null_q05": float(np.quantile(neutral, 0.05)),
                        "null_q95": float(np.quantile(neutral, 0.95)),
                        "null_q975": float(np.quantile(neutral, 0.975)),
                        "null_upper_exceedance_count": exceedances,
                        "mc_p_upper": pvalue,
                        "neutral_bank_scope": (
                            "same_demography_same_target_af_same_sample_size"
                        ),
                        "inference_scope": (
                            "selection_isolating_final_af_matched_genotype_contrast_sensitivity"
                        ),
                    }
                )
    result = pd.DataFrame(rows)
    result["bh_q_upper"] = np.nan
    for indices in result.groupby(["source", "unit_id"], sort=False).groups.values():
        result.loc[indices, "bh_q_upper"] = _bh_fdr(
            result.loc[indices, "mc_p_upper"].to_numpy(dtype=float)
        )
    columns = [
        "source",
        "estimand",
        "unit_id",
        "replicate_index",
        "threshold_years",
        "statistic",
        "observed_selection_score",
        "n_null",
        "null_mean",
        "null_median",
        "null_sd",
        "null_q025",
        "null_q05",
        "null_q95",
        "null_q975",
        "null_upper_exceedance_count",
        "mc_p_upper",
        "bh_q_upper",
        "neutral_bank_scope",
        "inference_scope",
    ]
    return result[columns]


def _inclusive_pooled_minp(
    neutral_values: np.ndarray, observed: np.ndarray
) -> tuple[np.ndarray, float, int, int, float]:
    neutral = np.asarray(neutral_values, dtype=float)
    selected = np.asarray(observed, dtype=float)
    if neutral.ndim != 2 or selected.shape != (neutral.shape[1],):
        raise ValueError("inclusive pooled minP inputs have incompatible shapes")
    if not np.isfinite(neutral).all() or not np.isfinite(selected).all():
        raise ValueError("inclusive pooled minP inputs are nonfinite")
    pooled = np.vstack((neutral, selected))
    marginal = np.empty_like(pooled)
    for index in range(pooled.shape[1]):
        values = pooled[:, index]
        marginal[:, index] = (values[None, :] >= values[:, None] - 1e-14).sum(
            axis=1
        ) / len(pooled)
    minp = marginal.min(axis=1)
    observed_min = float(minp[-1])
    exceedances = int(np.count_nonzero(minp[:-1] <= observed_min + 1e-14))
    best = int(np.argmin(marginal[-1]))
    reference_p = float((exceedances + 1) / len(pooled))
    return marginal[-1], observed_min, exceedances, best, reference_p


def _primary_minp(unit_scores: pd.DataFrame, pointwise: pd.DataFrame) -> pd.DataFrame:
    rows = []
    thresholds = np.asarray(THRESHOLDS_YEARS, dtype=float)
    for source in ("tree_truth", "gamma_smc"):
        source_scores = unit_scores[unit_scores["source"] == source]
        neutral = (
            source_scores[source_scores["bank_role"] == "survival_only_primary"]
            .pivot(index="unit_id", columns="threshold_years", values="p_overall")
            .reindex(columns=THRESHOLDS_YEARS)
            .sort_index()
        )
        selected = (
            source_scores[source_scores["bank_role"] == "selected"]
            .pivot(index="unit_id", columns="threshold_years", values="p_overall")
            .reindex(columns=THRESHOLDS_YEARS)
            .sort_index()
        )
        if neutral.shape != (100, 7) or selected.shape != (10, 7):
            raise ValueError("primary minP score matrix is incomplete")
        neutral_values = neutral.to_numpy(dtype=float)
        for unit_id, series in selected.iterrows():
            observed = series.to_numpy(dtype=float)
            (
                _observed_marginal,
                observed_min,
                exceedances,
                best,
                omnibus_p,
            ) = _inclusive_pooled_minp(neutral_values, observed)
            replicate = int(
                pointwise.loc[
                    (pointwise["source"] == source)
                    & (pointwise["estimand"] == PRIMARY_ESTIMAND)
                    & (pointwise["unit_id"] == unit_id),
                    "replicate_index",
                ].iloc[0]
            )
            rows.append(
                {
                    "source": source,
                    "estimand": PRIMARY_ESTIMAND,
                    "unit_id": unit_id,
                    "replicate_index": replicate,
                    "statistic": "p_overall",
                    "n_null": 100,
                    "n_thresholds": 7,
                    "thresholds_years": ",".join(map(str, THRESHOLDS_YEARS)),
                    "observed_min_marginal_p_upper": observed_min,
                    "most_extreme_threshold_years": int(thresholds[best]),
                    "observed_score_at_most_extreme_threshold": float(observed[best]),
                    "neutral_minp_exceedance_count": exceedances,
                    "omnibus_p_upper": omnibus_p,
                    "omnibus_method": (
                        "inclusive_pooled_101_marginal_ranks_minP_plus_one_reference_tail"
                    ),
                    "correlation_handling": (
                        "whole neutral seven-threshold vectors retained"
                    ),
                    "inference_scope": (
                        "nonexchangeable_simulation_reference_due_to_different_endpoint_af_conditioning"
                    ),
                }
            )
    return pd.DataFrame(rows)


def _sensitivity_minp(bundle: AnalysisBundle) -> pd.DataFrame:
    frame = bundle.tables["primary_omnibus"].copy()
    frame = frame[
        (frame["demography_id"].astype(str) == DEMOGRAPHY_ID)
        & np.isclose(pd.to_numeric(frame["target_allele_frequency"]), TARGET_AF)
        & np.isclose(
            pd.to_numeric(frame["selection_coefficient"]), SELECTION_COEFFICIENT
        )
        & (frame["statistic_scope"].astype(str) == "focal_nearest")
        & (frame["metric"].astype(str) == "p_tmrca_lt_threshold")
        & (frame["statistic"].astype(str) == "p_hom_alt_minus_hom_ref")
    ].copy()
    if len(frame) != SELECTED_REPLICATES:
        raise ValueError(f"{bundle.spec.key} sensitivity minP grid is incomplete")
    n_null = _exact_integer(frame, "n_null")
    exceed = _exact_integer(frame, "neutral_minp_exceedance_count")
    if not np.allclose(
        (exceed + 1) / (n_null + 1), _exact_numeric(frame, "omnibus_p_upper")
    ):
        raise ValueError(f"{bundle.spec.key} sensitivity minP +1 values differ")
    columns = [
        "unit_id",
        "replicate_index",
        "statistic",
        "n_null",
        "n_thresholds",
        "thresholds_years",
        "observed_min_marginal_p_upper",
        "most_extreme_threshold_years",
        "observed_score_at_most_extreme_threshold",
        "neutral_minp_exceedance_count",
        "omnibus_p_upper",
        "omnibus_method",
        "correlation_handling",
    ]
    result = frame[columns].copy()
    result.insert(0, "estimand", SENSITIVITY_ESTIMAND)
    result.insert(0, "source", bundle.spec.key)
    result["inference_scope"] = (
        "selection_isolating_final_af_matched_genotype_contrast_sensitivity"
    )
    return result


def _sensitivity_minp_from_scores(
    unit_scores: pd.DataFrame, pointwise: pd.DataFrame
) -> pd.DataFrame:
    rows = []
    thresholds = np.asarray(THRESHOLDS_YEARS, dtype=float)
    for source in ("tree_truth", "gamma_smc"):
        source_scores = unit_scores[unit_scores["source"] == source]
        neutral = (
            source_scores[source_scores["bank_role"] == "matched_af_sensitivity"]
            .pivot(
                index="unit_id",
                columns="threshold_years",
                values="p_hom_alt_minus_hom_ref",
            )
            .reindex(columns=THRESHOLDS_YEARS)
            .sort_index()
        )
        selected = (
            source_scores[source_scores["bank_role"] == "selected"]
            .pivot(
                index="unit_id",
                columns="threshold_years",
                values="p_hom_alt_minus_hom_ref",
            )
            .reindex(columns=THRESHOLDS_YEARS)
            .sort_index()
        )
        if neutral.shape != (100, 7) or selected.shape != (10, 7):
            raise ValueError("matched-AF direct minP matrix is incomplete")
        neutral_values = neutral.to_numpy(dtype=float)
        for unit_id, series in selected.iterrows():
            observed = series.to_numpy(dtype=float)
            _, observed_min, exceedances, best, omnibus_p = _inclusive_pooled_minp(
                neutral_values, observed
            )
            replicate = int(
                pointwise.loc[
                    (pointwise["source"] == source)
                    & (pointwise["estimand"] == SENSITIVITY_ESTIMAND)
                    & (pointwise["unit_id"] == unit_id),
                    "replicate_index",
                ].iloc[0]
            )
            rows.append(
                {
                    "source": source,
                    "estimand": SENSITIVITY_ESTIMAND,
                    "unit_id": unit_id,
                    "replicate_index": replicate,
                    "statistic": "p_hom_alt_minus_hom_ref",
                    "n_null": 100,
                    "n_thresholds": 7,
                    "thresholds_years": ",".join(map(str, THRESHOLDS_YEARS)),
                    "observed_min_marginal_p_upper": observed_min,
                    "most_extreme_threshold_years": int(thresholds[best]),
                    "observed_score_at_most_extreme_threshold": float(observed[best]),
                    "neutral_minp_exceedance_count": exceedances,
                    "omnibus_p_upper": omnibus_p,
                    "omnibus_method": (
                        "inclusive_pooled_101_marginal_ranks_minP_plus_one_reference_tail"
                    ),
                    "correlation_handling": (
                        "whole neutral seven-threshold vectors retained"
                    ),
                    "inference_scope": (
                        "selection_isolating_final_af_matched_genotype_contrast_sensitivity"
                    ),
                }
            )
    return pd.DataFrame(rows)


def _wilson_interval(successes: int, total: int) -> tuple[float, float, float]:
    if total < 1 or not 0 <= successes <= total:
        raise ValueError("Wilson interval counts are invalid")
    z = 1.959963984540054
    fraction = successes / total
    denominator = 1 + z * z / total
    center = (fraction + z * z / (2 * total)) / denominator
    half = (
        z
        * np.sqrt(fraction * (1 - fraction) / total + z * z / (4 * total**2))
        / denominator
    )
    return fraction, max(0.0, center - half), min(1.0, center + half)


def _power_from_pointwise(pointwise: pd.DataFrame) -> pd.DataFrame:
    rows = []
    group_columns = ["source", "estimand", "statistic", "threshold_years"]
    for labels, group in pointwise.groupby(group_columns, sort=True):
        raw = group["mc_p_upper"].to_numpy(dtype=float)
        adjusted = group["bh_q_upper"].to_numpy(dtype=float)
        if len(raw) != SELECTED_REPLICATES:
            raise ValueError("power summary does not contain ten selected replicates")
        raw_n = int(np.count_nonzero(raw <= 0.05))
        adjusted_n = int(np.count_nonzero(adjusted <= 0.05))
        raw_fraction, raw_low, raw_high = _wilson_interval(raw_n, len(raw))
        adjusted_fraction, adjusted_low, adjusted_high = _wilson_interval(
            adjusted_n, len(raw)
        )
        rows.append(
            {
                **dict(zip(group_columns, labels, strict=True)),
                "alpha": 0.05,
                "n_selected": len(raw),
                "n_detected_raw": raw_n,
                "conditional_power_raw": raw_fraction,
                "raw_wilson_ci95_low": raw_low,
                "raw_wilson_ci95_high": raw_high,
                "n_detected_bh_within_7": adjusted_n,
                "conditional_power_bh_within_7": adjusted_fraction,
                "bh_within_7_wilson_ci95_low": adjusted_low,
                "bh_within_7_wilson_ci95_high": adjusted_high,
                "bh_family_size": 7,
                "interpretation": (
                    "simulation_reference_detection_fraction"
                    if labels[1] == PRIMARY_ESTIMAND
                    else "final_af_matched_conditional_power_sensitivity"
                ),
            }
        )
    return pd.DataFrame(rows)


def _verify_sensitivity_power(
    bundle: AnalysisBundle, derived_power: pd.DataFrame
) -> None:
    published = bundle.tables["conditional_power"].copy()
    published = published[
        (published["demography_id"].astype(str) == DEMOGRAPHY_ID)
        & np.isclose(pd.to_numeric(published["target_allele_frequency"]), TARGET_AF)
        & np.isclose(
            pd.to_numeric(published["selection_coefficient"]), SELECTION_COEFFICIENT
        )
        & (published["statistic_scope"].astype(str) == "focal_nearest")
        & (published["metric"].astype(str) == "p_tmrca_lt_threshold")
        & (published["statistic"].astype(str) == "p_hom_alt_minus_hom_ref")
        & pd.to_numeric(published["threshold_years"]).isin(THRESHOLDS_YEARS)
    ].copy()
    derived = derived_power[
        (derived_power["source"] == bundle.spec.key)
        & (derived_power["estimand"] == SENSITIVITY_ESTIMAND)
    ].sort_values("threshold_years")
    published = published.sort_values("threshold_years")
    if len(published) != 7 or len(derived) != 7:
        raise ValueError(f"{bundle.spec.key} sensitivity power grid is incomplete")
    comparisons = {
        "n_detected_raw": "n_detected_raw",
        "conditional_power_raw": "conditional_power_raw",
        "n_detected_bh_within_7": "n_detected_bh_within_7",
        "conditional_power_bh_within_7": "conditional_power_bh_within_7",
    }
    for left, right in comparisons.items():
        if not np.allclose(
            pd.to_numeric(published[left]).to_numpy(dtype=float),
            pd.to_numeric(derived[right]).to_numpy(dtype=float),
        ):
            raise ValueError(f"{bundle.spec.key} sensitivity power differs: {left}")


def _comparator_inputs(
    repo_root: Path,
    *,
    truth_metrics_path: Path,
    gamma_metrics_path: Path,
    gamma_results_path: Path,
    workers: int,
) -> tuple[dict[str, Any], dict[str, Any], pd.DataFrame, dict[str, Any]]:
    paths = (truth_metrics_path, gamma_metrics_path, gamma_results_path)
    for path in paths:
        _require_regular_file(path, label="two-epoch comparator input")
    checks = [
        FileCheck(
            "two_epoch_comparator",
            "linked_low_af",
            "",
            "comparator_input",
            key,
            path,
            _path_identity(path, repo_root),
            None,
            None,
        )
        for key, path in zip(
            ("truth_metrics", "gamma_metrics", "gamma_results"), paths, strict=True
        )
    ]
    audit = _audit_file_checks(checks, workers=workers)
    truth = _read_json(truth_metrics_path)
    gamma = _read_json(gamma_metrics_path)
    gamma_results_text = gamma_results_path.read_text(encoding="utf-8")
    truth_selected = truth.get("selected_rejection", {})
    if not (
        np.isclose(float(truth.get("selection_coefficient")), 0.05)
        and np.isclose(float(truth.get("mutation_rate")), 1.29e-9)
        and int(truth.get("sample_diploids")) == 2_000
        and int(truth.get("neutral_replicates")) == 100
        and np.isclose(float(truth.get("variant_age_years")), 4_500)
        and int(truth_selected.get("accepted_attempt_zero_based")) == 8
        and np.isclose(
            float(truth_selected.get("population_allele_frequency")), 0.102125
        )
        and np.isclose(float(truth_selected.get("sample_allele_frequency")), 0.10275)
        and np.isclose(float(truth_selected.get("observed_fraction_recent")), 0.0185)
        and int(truth_selected.get("neutral_exceedances")) == 0
        and np.isclose(float(truth_selected.get("monte_carlo_p_upper")), 1 / 101)
    ):
        raise ValueError("two-epoch truth comparator is incompatible")
    if not (
        gamma.get("data_interpretation")
        == "retained selected simulation treated as pseudo-empirical data"
        and gamma.get("container_image") == "docker.io/regevsch/gamma_smc:v0.2"
        and gamma.get("container_runtime") == "docker"
        and int(gamma.get("sample_diploids")) == 2_000
        and int(gamma.get("within_individual_pairs")) == 2_000
        and int(gamma.get("neutral_replicates")) == 100
        and int(gamma.get("stride_bp")) == 1_000
        and np.isclose(float(gamma.get("threshold_years")), 4_500)
        and int(gamma.get("center_neutral_exceedances")) == 13
        and np.isclose(
            float(gamma.get("center_observed_mean_p_recent")), 0.0305023913780933
        )
        and np.isclose(
            float(gamma.get("center_null_mean_p_recent")), 0.023461553284794868
        )
        and np.isclose(float(gamma.get("center_monte_carlo_p_upper")), 14 / 101)
        and "mu=1.25e-8" in gamma_results_text
        and "official `regevsch/gamma_smc:v0.2` container" in gamma_results_text
    ):
        raise ValueError("two-epoch Gamma-SMC comparator is incompatible")
    return (
        truth,
        gamma,
        audit,
        {
            "truth_metrics": _file_record(truth_metrics_path, repo_root),
            "gamma_metrics": _file_record(gamma_metrics_path, repo_root),
            "gamma_results": _file_record(gamma_results_path, repo_root),
            "relationship": (
                "same linked low-AF selected genealogy, but not the same mutation "
                "layer: truth folder mu=1.29e-9 versus corresponding official v0.2 "
                "container decode mu=1.25e-8"
            ),
        },
    )


def _comparison_design(
    pointwise: pd.DataFrame,
    power: pd.DataFrame,
    truth_metrics: Mapping[str, Any],
    gamma_metrics: Mapping[str, Any],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for source, label in (("tree_truth", "Tree truth"), ("gamma_smc", "Gamma-SMC")):
        selected = pointwise[
            (pointwise["source"] == source)
            & (pointwise["estimand"] == PRIMARY_ESTIMAND)
            & (pointwise["threshold_years"] == DISPLAY_THRESHOLD_YEARS)
        ]
        power_row = power[
            (power["source"] == source)
            & (power["estimand"] == PRIMARY_ESTIMAND)
            & (power["threshold_years"] == DISPLAY_THRESHOLD_YEARS)
        ].iloc[0]
        rows.append(
            {
                "comparison_id": f"ancient_eurasia_{source}_primary",
                "model": "AncientEurasia_9K19 Han introgression",
                "source": label,
                "relationship_group": "ancient_eurasia_primary",
                "same_selected_genealogy_within_group": True,
                "selection_coefficient": SELECTION_COEFFICIENT,
                "selected_origin_years": 60_000.0,
                "introgression_pulse_years": 56_800.0,
                "selection_end_years": 44_500.0,
                "selected_replicates": SELECTED_REPLICATES,
                "neutral_replicates": SURVIVAL_NEUTRAL_REPLICATES,
                "sample_diploids": SAMPLE_DIPLOIDS,
                "within_diploid_pairs": SAMPLE_DIPLOIDS,
                "selected_endpoint_conditioning": "exact sampled AF 20%",
                "neutral_endpoint_conditioning": "present Han census segregating only",
                "mutation_rate": 1.25e-8,
                "output_stride_bp": 10_000 if source == "gamma_smc" else np.nan,
                "display_threshold_years": DISPLAY_THRESHOLD_YEARS,
                "statistic": "focal overall P(TMRCA<x)",
                "selected_score": float(
                    np.median(selected["observed_selection_score"])
                ),
                "neutral_score": float(np.median(selected["null_median"])),
                "neutral_exceedances": np.nan,
                "reference_p_upper": float(np.median(selected["mc_p_upper"])),
                "reference_detection_fraction": float(
                    power_row["conditional_power_raw"]
                ),
                "interpretation": (
                    "10-replicate simulation-reference summary; not a single paired p-value"
                ),
                "mutation_layer_relationship": (
                    "source-specific AncientEurasia score bank; no cross-source pairing claim"
                ),
                "decoder_backend": (
                    "not applicable (tree truth)"
                    if source == "tree_truth"
                    else "frozen focused-study Gamma-SMC decode"
                ),
            }
        )
    truth_selected = truth_metrics["selected_rejection"]
    rows.append(
        {
            "comparison_id": "two_epoch_low_af_tree_truth_focal",
            "model": "two-epoch growth",
            "source": "Tree truth",
            "relationship_group": "two_epoch_low_af_same_selected_base",
            "same_selected_genealogy_within_group": True,
            "selection_coefficient": float(truth_metrics["selection_coefficient"]),
            "selected_origin_years": float(truth_metrics["variant_age_years"]),
            "introgression_pulse_years": np.nan,
            "selection_end_years": 0.0,
            "selected_replicates": 1,
            "neutral_replicates": int(truth_metrics["neutral_replicates"]),
            "sample_diploids": int(truth_metrics["sample_diploids"]),
            "within_diploid_pairs": int(truth_metrics["sample_diploids"]),
            "selected_endpoint_conditioning": (
                f"population AF {truth_selected['population_allele_frequency']:.6g}; "
                f"sample AF {truth_selected['sample_allele_frequency']:.6g}"
            ),
            "neutral_endpoint_conditioning": "unconditioned neutral",
            "mutation_rate": float(truth_metrics["mutation_rate"]),
            "output_stride_bp": np.nan,
            "display_threshold_years": float(truth_metrics["variant_age_years"]),
            "statistic": "truth fraction of all pairs with TMRCA<x at focal site",
            "selected_score": float(truth_selected["observed_fraction_recent"]),
            "neutral_score": float(
                truth_metrics["neutral_observed_mean_fraction_recent"]
            ),
            "neutral_exceedances": int(truth_selected["neutral_exceedances"]),
            "reference_p_upper": float(truth_selected["monte_carlo_p_upper"]),
            "reference_detection_fraction": np.nan,
            "interpretation": "attached histogram is tree truth, not Gamma-SMC",
            "mutation_layer_relationship": (
                "truth layer mu=1.29e-9; differs from v0.2 decode layer"
            ),
            "decoder_backend": "not applicable (tree truth)",
        }
    )
    rows.append(
        {
            "comparison_id": "two_epoch_low_af_gamma_smc_focal",
            "model": "two-epoch growth",
            "source": "Gamma-SMC",
            "relationship_group": "two_epoch_low_af_same_selected_base",
            "same_selected_genealogy_within_group": True,
            "selection_coefficient": float(truth_metrics["selection_coefficient"]),
            "selected_origin_years": float(truth_metrics["variant_age_years"]),
            "introgression_pulse_years": np.nan,
            "selection_end_years": 0.0,
            "selected_replicates": 1,
            "neutral_replicates": int(gamma_metrics["neutral_replicates"]),
            "sample_diploids": int(gamma_metrics["sample_diploids"]),
            "within_diploid_pairs": int(gamma_metrics["within_individual_pairs"]),
            "selected_endpoint_conditioning": (
                f"same selected base; sample AF {truth_selected['sample_allele_frequency']:.6g}"
            ),
            "neutral_endpoint_conditioning": "matched demography, decoded identically",
            "mutation_rate": 1.25e-8,
            "output_stride_bp": int(gamma_metrics["stride_bp"]),
            "display_threshold_years": float(gamma_metrics["threshold_years"]),
            "statistic": "Gamma-SMC mean P(TMRCA<x) across all pairs at focal site",
            "selected_score": float(gamma_metrics["center_observed_mean_p_recent"]),
            "neutral_score": float(gamma_metrics["center_null_mean_p_recent"]),
            "neutral_exceedances": int(gamma_metrics["center_neutral_exceedances"]),
            "reference_p_upper": float(gamma_metrics["center_monte_carlo_p_upper"]),
            "reference_detection_fraction": np.nan,
            "interpretation": (
                "focal Gamma-SMC was not significant; regional density, not focal p, "
                "was significant; official v0.2 container, not native-optimized"
            ),
            "mutation_layer_relationship": (
                "decode layer mu=1.25e-8; same genealogy but differs from truth layer"
            ),
            "decoder_backend": "official regevsch/gamma_smc:v0.2 container",
        }
    )
    local = gamma_metrics["scan_density_calibration"]["local"]
    rows.append(
        {
            "comparison_id": "two_epoch_low_af_gamma_smc_local_density",
            "model": "two-epoch growth",
            "source": "Gamma-SMC regional",
            "relationship_group": "two_epoch_low_af_same_selected_base",
            "same_selected_genealogy_within_group": True,
            "selection_coefficient": float(truth_metrics["selection_coefficient"]),
            "selected_origin_years": float(truth_metrics["variant_age_years"]),
            "introgression_pulse_years": np.nan,
            "selection_end_years": 0.0,
            "selected_replicates": 1,
            "neutral_replicates": int(gamma_metrics["neutral_replicates"]),
            "sample_diploids": int(gamma_metrics["sample_diploids"]),
            "within_diploid_pairs": int(gamma_metrics["within_individual_pairs"]),
            "selected_endpoint_conditioning": "same selected base",
            "neutral_endpoint_conditioning": "leave-one-out scan-density reference",
            "mutation_rate": 1.25e-8,
            "output_stride_bp": int(gamma_metrics["stride_bp"]),
            "display_threshold_years": float(gamma_metrics["threshold_years"]),
            "statistic": "raw-p-positive windows within +/-500 kb",
            "selected_score": float(local["selected_positive_windows"]),
            "neutral_score": float(local["null_mean_positive_windows"]),
            "neutral_exceedances": np.nan,
            "reference_p_upper": float(local["monte_carlo_p_upper"]),
            "reference_detection_fraction": np.nan,
            "interpretation": "regional density is a different statistic from focal p",
            "mutation_layer_relationship": (
                "decode layer mu=1.25e-8; same genealogy but differs from truth layer"
            ),
            "decoder_backend": "official regevsch/gamma_smc:v0.2 container",
        }
    )
    return pd.DataFrame(rows)


def _rug_selected(
    axis: plt.Axes, values: np.ndarray, *, color: str, label: str
) -> None:
    for index, value in enumerate(values):
        axis.axvline(
            value,
            ymin=0.79 + 0.018 * (index % 3),
            ymax=0.97,
            color=color,
            linewidth=1.8,
            alpha=0.72,
            label=label if index == 0 else None,
        )
    axis.axvline(
        float(np.median(values)),
        color=color,
        linewidth=4,
        label="Selected median",
    )


def plot_pulse_aligned_null(
    unit_scores: pd.DataFrame,
    endpoints: pd.DataFrame,
    *,
    source: str,
    output_dir: Path,
) -> list[Path]:
    source_label = "Tree truth" if source == "tree_truth" else "Gamma-SMC"
    at_threshold = unit_scores[
        (unit_scores["source"] == source)
        & (unit_scores["threshold_years"] == DISPLAY_THRESHOLD_YEARS)
    ]
    selected = at_threshold[at_threshold["bank_role"] == "selected"].sort_values(
        "replicate_index"
    )
    survival = at_threshold[at_threshold["bank_role"] == "survival_only_primary"]
    matched = at_threshold[at_threshold["bank_role"] == "matched_af_sensitivity"]
    if (len(selected), len(survival), len(matched)) != (10, 100, 100):
        raise ValueError(f"{source} pulse-aligned plot inputs are incomplete")
    figure, axes = plt.subplots(1, 3, figsize=(11, 8.5), squeeze=False)
    axes = axes[0]
    neutral_color = "#999999"
    selected_color = "#6F3FE1"
    mean_color = "#198754"

    overall_null = survival["p_overall"].to_numpy(dtype=float)
    overall_selected = selected["p_overall"].to_numpy(dtype=float)
    axes[0].hist(overall_null, bins=18, color=neutral_color, alpha=0.9)
    axes[0].axvline(
        np.mean(overall_null),
        color=mean_color,
        linestyle="--",
        linewidth=3,
        label="Neutral mean",
    )
    _rug_selected(
        axes[0], overall_selected, color=selected_color, label="Selected replicates"
    )
    overall_p = np.asarray(
        [_upper_reference_pvalue(value, overall_null)[1] for value in overall_selected]
    )
    axes[0].text(
        0.03,
        0.97,
        "Reference tails\n"
        f"median {np.median(overall_p):.3f}\n"
        f"range {overall_p.min():.3f}–{overall_p.max():.3f}\n"
        f"{np.count_nonzero(overall_p <= 0.05)}/10 ≤ 0.05",
        transform=axes[0].transAxes,
        ha="left",
        va="top",
        fontsize=10.8,
    )
    axes[0].set_title("Primary simulation reference\nOverall vs overall", fontsize=18)
    axes[0].set_xlabel("P(TMRCA < 50 kya)", fontsize=15)
    axes[0].set_ylabel("Neutral simulations", fontsize=16)

    x = np.arange(1, SELECTED_REPLICATES + 1)
    rr = selected["p_hom_ref"].to_numpy(dtype=float)
    aa = selected["p_hom_alt"].to_numpy(dtype=float)
    for index in range(SELECTED_REPLICATES):
        axes[1].plot(
            [x[index], x[index]],
            [rr[index], aa[index]],
            color="#BBBBBB",
            linewidth=1.5,
            zorder=1,
        )
    axes[1].scatter(x, rr, marker="o", s=45, color="#555555", label="Ref/ref")
    axes[1].scatter(x, aa, marker="^", s=58, color="#0072B2", label="Alt/alt")
    axes[1].axhline(np.median(rr), color="#555555", linestyle="--", linewidth=2)
    axes[1].axhline(np.median(aa), color="#0072B2", linestyle="--", linewidth=2)
    axes[1].set_title("Selected panels\nGenotype description", fontsize=18)
    axes[1].set_xlabel("Selected replicate", fontsize=15)
    axes[1].set_ylabel("P(TMRCA < 50 kya)", fontsize=16)
    axes[1].set_xticks([1, 5, 10])
    selected_endpoints = endpoints[endpoints["bank_role"] == "selected"].sort_values(
        "replicate_index"
    )
    if len(selected_endpoints) != SELECTED_REPLICATES:
        raise ValueError("selected endpoint counts are incomplete for genotype plot")
    aa_n = _exact_integer(selected_endpoints, "hom_alt_diploids")
    rr_n = _exact_integer(selected_endpoints, "hom_ref_diploids")
    axes[1].text(
        0.03,
        0.97,
        f"AA n={aa_n.min()}-{aa_n.max()} "
        f"(median {np.median(aa_n):.1f})\n"
        f"RR n={rr_n.min()}-{rr_n.max()} "
        f"(median {np.median(rr_n):.1f})",
        transform=axes[1].transAxes,
        ha="left",
        va="top",
        fontsize=11.5,
        wrap=True,
    )

    delta_null = matched["p_hom_alt_minus_hom_ref"].to_numpy(dtype=float)
    delta_selected = selected["p_hom_alt_minus_hom_ref"].to_numpy(dtype=float)
    axes[2].hist(delta_null, bins=18, color=neutral_color, alpha=0.9)
    axes[2].axvline(np.mean(delta_null), color=mean_color, linestyle="--", linewidth=3)
    _rug_selected(
        axes[2], delta_selected, color=selected_color, label="Selected replicates"
    )
    sensitivity_p = np.asarray(
        [_upper_reference_pvalue(value, delta_null)[1] for value in delta_selected]
    )
    axes[2].text(
        0.03,
        0.97,
        "AF-matched tails\n"
        f"median {np.median(sensitivity_p):.3f}\n"
        f"range {sensitivity_p.min():.3f}–{sensitivity_p.max():.3f}\n"
        f"{np.count_nonzero(sensitivity_p <= 0.05)}/10 ≤ 0.05",
        transform=axes[2].transAxes,
        ha="left",
        va="top",
        fontsize=10.8,
    )
    axes[2].set_title("AF-matched sensitivity\nAlt/alt minus ref/ref", fontsize=18)
    axes[2].set_xlabel("Alt/alt - ref/ref probability", fontsize=15)
    axes[2].set_ylabel("Matched-neutral simulations", fontsize=16)

    for axis in axes:
        axis.tick_params(labelsize=13)
        axis.grid(axis="y", alpha=0.2)
    handles = [
        Line2D([0], [0], color=neutral_color, linewidth=8, label="Neutral null"),
        Line2D(
            [0], [0], color=mean_color, linestyle="--", linewidth=3, label="Null mean"
        ),
        Line2D(
            [0], [0], color=selected_color, linewidth=2, label="10 selected replicates"
        ),
        Line2D([0], [0], color=selected_color, linewidth=4, label="Selected median"),
        Line2D(
            [0],
            [0],
            color="#555555",
            marker="o",
            linestyle="",
            label="Selected ref/ref",
        ),
        Line2D(
            [0],
            [0],
            color="#0072B2",
            marker="^",
            linestyle="",
            label="Selected alt/alt",
        ),
    ]
    figure.suptitle(
        f"{source_label}: AncientEurasia 9K19 at 50 kya\n"
        "Pulse-aligned display cutoff\n"
        "Seven-threshold simulation-reference minP remains primary",
        fontsize=18.5,
        y=0.99,
    )
    figure.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.80),
        ncol=3,
        fontsize=11,
        frameon=False,
    )
    figure.subplots_adjust(left=0.09, right=0.975, bottom=0.15, top=0.62, wspace=0.36)
    stem = f"{source}_pulse_aligned_null"
    paths = _save_pair(figure, output_dir, stem)
    plt.close(figure)
    return paths


def plot_threshold_results(
    unit_scores: pd.DataFrame,
    pointwise: pd.DataFrame,
    power: pd.DataFrame,
    *,
    output_dir: Path,
) -> list[Path]:
    figure, axes = plt.subplots(2, 2, figsize=(11, 8.5), sharex="col", squeeze=False)
    x = np.asarray(THRESHOLDS_YEARS, dtype=float) / 1_000
    selected_color = "#6F3FE1"
    neutral_color = "#555555"
    power_color = "#0072B2"
    for column, (source, label) in enumerate(
        (("tree_truth", "Tree truth"), ("gamma_smc", "Gamma-SMC"))
    ):
        scores = unit_scores[unit_scores["source"] == source]
        selected = (
            scores[scores["bank_role"] == "selected"]
            .pivot(index="unit_id", columns="threshold_years", values="p_overall")
            .reindex(columns=THRESHOLDS_YEARS)
        )
        neutral = (
            scores[scores["bank_role"] == "survival_only_primary"]
            .pivot(index="unit_id", columns="threshold_years", values="p_overall")
            .reindex(columns=THRESHOLDS_YEARS)
        )
        top = axes[0, column]
        for values in selected.to_numpy(dtype=float):
            top.plot(x, values, color=selected_color, alpha=0.22, linewidth=1.3)
        top.fill_between(
            x,
            np.quantile(neutral.to_numpy(dtype=float), 0.05, axis=0),
            np.quantile(neutral.to_numpy(dtype=float), 0.95, axis=0),
            color=neutral_color,
            alpha=0.18,
        )
        top.plot(
            x,
            np.median(neutral.to_numpy(dtype=float), axis=0),
            color=neutral_color,
            linestyle="--",
            linewidth=2.7,
        )
        top.plot(
            x,
            np.median(selected.to_numpy(dtype=float), axis=0),
            color=selected_color,
            marker="o",
            linewidth=3.2,
        )
        top.set_title(label, fontsize=19)
        top.set_ylabel("Focal P(TMRCA < x)", fontsize=15)

        bottom = axes[1, column]
        pvalues = pointwise[
            (pointwise["source"] == source)
            & (pointwise["estimand"] == PRIMARY_ESTIMAND)
        ].pivot(index="unit_id", columns="threshold_years", values="mc_p_upper")
        pvalues = pvalues.reindex(columns=THRESHOLDS_YEARS)
        for values in pvalues.to_numpy(dtype=float):
            bottom.plot(x, values, color=selected_color, alpha=0.30, linewidth=1.3)
        power_curve = power[
            (power["source"] == source) & (power["estimand"] == PRIMARY_ESTIMAND)
        ].sort_values("threshold_years")
        bottom.plot(
            x,
            power_curve["conditional_power_raw"],
            color=power_color,
            marker="s",
            linewidth=3.0,
        )
        bottom.axhline(0.05, color="#BB2222", linestyle=":", linewidth=2.2)
        bottom.set_ylim(-0.02, 1.02)
        bottom.set_ylabel("Reference p / p≤0.05 fraction", fontsize=15)
        bottom.set_xlabel("TMRCA threshold (kya)", fontsize=16)
        bottom.set_xticks(x, [f"{value:g}" for value in x])
    for axis in axes.ravel():
        axis.tick_params(labelsize=12.5)
        axis.grid(alpha=0.2)
    handles = [
        Line2D(
            [0],
            [0],
            color=selected_color,
            alpha=0.35,
            linewidth=2,
            label="Each selected replicate",
        ),
        Line2D(
            [0],
            [0],
            color=selected_color,
            marker="o",
            linewidth=3,
            label="Selected median",
        ),
        Line2D(
            [0],
            [0],
            color=neutral_color,
            linestyle="--",
            linewidth=3,
            label="Neutral median",
        ),
        Line2D(
            [0],
            [0],
            color=neutral_color,
            linewidth=8,
            alpha=0.18,
            label="Neutral 5–95% band",
        ),
        Line2D(
            [0],
            [0],
            color=power_color,
            marker="s",
            linewidth=3,
            label="p≤0.05 fraction",
        ),
        Line2D(
            [0],
            [0],
            color="#BB2222",
            linestyle=":",
            linewidth=2.2,
            label="α = 0.05",
        ),
    ]
    figure.suptitle(
        "AncientEurasia: all seven TMRCA thresholds\n"
        "Overall selected vs survival-neutral; tails include endpoint-AF enrichment",
        fontsize=18.5,
        y=0.985,
    )
    figure.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.86),
        ncol=3,
        fontsize=11.5,
        frameon=False,
    )
    figure.subplots_adjust(
        left=0.11, right=0.98, bottom=0.11, top=0.72, hspace=0.34, wspace=0.25
    )
    paths = _save_pair(figure, output_dir, "threshold_scores_pvalues_power")
    plt.close(figure)
    return paths


def plot_timing_endpoints_diagnostic(
    endpoints: pd.DataFrame,
    comparison: pd.DataFrame,
    *,
    output_dir: Path,
) -> list[Path]:
    figure, axes = plt.subplots(2, 2, figsize=(11, 8.5), squeeze=False)
    timing = axes[0, 0]
    timing.set_xlim(65, -2)
    timing.set_ylim(0, 1)
    timing.hlines(0.48, 0, 60, color="#555555", linewidth=2)
    events = [
        (60.0, "Origin\n60 kya", 63.5, 0.68),
        (56.8, "Pulse\n56.8 kya", 60.0, 0.88),
        (50.375, "Han split\n50.4 kya", 49.0, 0.68),
        (44.5, "Selection ends\n44.5 kya", 41.5, 0.88),
        (0.0, "Present", 0.0, 0.80),
    ]
    for value, label, label_x, label_y in events:
        timing.vlines(value, 0.34, 0.63, color="#333333", linewidth=2)
        timing.annotate(
            label,
            xy=(value, 0.63),
            xytext=(label_x, label_y),
            ha="center",
            va="center",
            fontsize=10.8,
            arrowprops={"arrowstyle": "-", "color": "#555555", "linewidth": 1.0},
        )
    timing.axvspan(44.5, 56.8, ymin=0.43, ymax=0.54, color="#6F3FE1", alpha=0.35)
    timing.text(50.65, 0.18, "s=0.01 finite episode", ha="center", fontsize=14)
    timing.set_xlabel("Thousands of years before present", fontsize=15)
    timing.set_yticks([])
    timing.set_title("Timing schematic (not a trajectory)", fontsize=17)

    endpoint_axis = axes[0, 1]
    groups = [
        ("selected", "Selected", "#6F3FE1"),
        ("survival_only_primary", "Survival neutral", "#198754"),
        ("matched_af_sensitivity", "AF-matched neutral", "#777777"),
    ]
    for x_index, (role, label, color) in enumerate(groups):
        group = endpoints[endpoints["bank_role"] == role].sort_values("replicate_index")
        offset = (
            np.linspace(-0.12, 0.12, len(group)) if len(group) > 1 else np.array([0.0])
        )
        endpoint_axis.scatter(
            x_index + offset,
            group["final_pool_af"],
            s=22,
            alpha=0.55,
            color=color,
            label=f"{label}: observed 500-pool AF",
        )
        endpoint_axis.scatter(
            [x_index],
            [np.median(group["sample_af"])],
            s=95,
            facecolors="none",
            edgecolors=color,
            linewidths=2.5,
            marker="o",
        )
    endpoint_axis.axhline(0.2, color="#BB2222", linestyle="--", linewidth=2)
    endpoint_axis.set_xticks(
        range(3), ["Selected", "Survival\nneutral", "AF-matched\nneutral"]
    )
    endpoint_axis.set_ylim(-0.02, 1.02)
    endpoint_axis.set_ylabel("Allele frequency", fontsize=15)
    endpoint_axis.set_title(
        "Observed endpoint AF\n(no stored trajectories)", fontsize=17
    )

    diagnostic = axes[1, 0]
    diagnostic_rows = comparison[
        comparison["comparison_id"].isin(
            [
                "two_epoch_low_af_tree_truth_focal",
                "two_epoch_low_af_gamma_smc_focal",
                "two_epoch_low_af_gamma_smc_local_density",
            ]
        )
    ]
    order = [
        "two_epoch_low_af_tree_truth_focal",
        "two_epoch_low_af_gamma_smc_focal",
        "two_epoch_low_af_gamma_smc_local_density",
    ]
    diagnostic_rows = diagnostic_rows.set_index("comparison_id").loc[order]
    pvalues = diagnostic_rows["reference_p_upper"].to_numpy(dtype=float)
    diagnostic.plot([0, 1], pvalues[:2], color="#777777", linewidth=2, zorder=1)
    diagnostic.scatter([0], [pvalues[0]], color="#6F3FE1", s=100, marker="o", zorder=2)
    diagnostic.scatter([1], [pvalues[1]], color="#0072B2", s=100, marker="s", zorder=2)
    diagnostic.scatter([2], [pvalues[2]], color="#E69F00", s=110, marker="^", zorder=2)
    for index, value in enumerate(pvalues):
        diagnostic.text(
            index, value + 0.018, f"p={value:.3f}", ha="center", fontsize=13
        )
    diagnostic.axhline(0.05, color="#BB2222", linestyle="--", linewidth=2)
    diagnostic.set_xticks(
        [0, 1, 2],
        ["Tree truth\nfocal", "Gamma-SMC\nfocal", "Gamma-SMC\n±500 kb density"],
    )
    diagnostic.set_ylim(0, max(0.22, float(np.max(pvalues)) + 0.07))
    diagnostic.set_ylabel("Upper-tail reference p", fontsize=15)
    diagnostic.set_xlabel(
        "Linked genealogy; different mutation layers\n"
        "truth mu=1.29e-9; Gamma v0.2 mu=1.25e-8",
        fontsize=10.5,
    )
    diagnostic.set_title("Two-epoch focal diagnostic", fontsize=17)

    design = axes[1, 1]
    design.axis("off")
    table = design.table(
        cellText=[
            ["Selection", "s=0.01", "s=0.05"],
            ["Origin", "60 kya, archaic", "4.5 kya, de novo"],
            ["Selected pairs", "100", "2,000"],
            ["Gamma stride", "10 kb", "1 kb"],
            ["Mutation layer", "mu=1.25e-8", "truth 1.29e-9\nGamma 1.25e-8"],
            ["Primary null", "survival only", "unconditioned"],
        ],
        colLabels=["Feature", "AncientEurasia", "Two epoch"],
        colWidths=[0.30, 0.30, 0.40],
        loc="center",
        cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10.2)
    table.scale(1.0, 1.90)
    for column in range(3):
        mutation_cell = table[(5, column)]
        mutation_cell.set_height(mutation_cell.get_height() * 1.35)
    design.set_title("Why the designs differ", fontsize=17, pad=10)

    for axis in (timing, endpoint_axis, diagnostic):
        axis.tick_params(labelsize=12.5)
        axis.grid(axis="y", alpha=0.18)
    handles = [
        Line2D(
            [0],
            [0],
            color="#6F3FE1",
            marker="o",
            linestyle="",
            label="Selected / tree truth",
        ),
        Line2D(
            [0],
            [0],
            color="#198754",
            marker="o",
            linestyle="",
            label="Survival neutral",
        ),
        Line2D(
            [0], [0], color="#0072B2", marker="s", linestyle="", label="Gamma-SMC focal"
        ),
        Line2D(
            [0],
            [0],
            color="#E69F00",
            marker="^",
            linestyle="",
            label="Regional density (different statistic)",
        ),
        Line2D(
            [0],
            [0],
            color="#555555",
            marker="o",
            markerfacecolor="#555555",
            linestyle="",
            label="Each observed 500-pool AF (filled)",
        ),
        Line2D(
            [0],
            [0],
            color="#555555",
            marker="o",
            markerfacecolor="none",
            markeredgewidth=2,
            linestyle="",
            label="Median observed 100-panel AF (open)",
        ),
    ]
    figure.suptitle(
        "AncientEurasia timing, endpoints, and two-epoch diagnostic",
        fontsize=19,
        y=0.985,
    )
    figure.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.91),
        ncol=3,
        fontsize=10.8,
        frameon=False,
    )
    figure.subplots_adjust(
        left=0.09, right=0.975, bottom=0.14, top=0.76, hspace=0.52, wspace=0.30
    )
    paths = _save_pair(figure, output_dir, "timing_endpoints_two_epoch_diagnostic")
    plt.close(figure)
    return paths


def _implementation_contract(repo_root: Path) -> tuple[dict[str, Any], pd.DataFrame]:
    paths = (
        Path(__file__).resolve(),
        repo_root / "scripts/run_ancient_eurasia_clear_results.py",
    )
    records: dict[str, Any] = {}
    rows = []
    for path in paths:
        _require_regular_file(path, boundary=repo_root, label="clear-results source")
        sha = _canonical_text_sha256(path)
        size = path.stat().st_size
        identity = _path_identity(path, repo_root)
        records[path.name] = {"path": identity, "sha256": sha}
        rows.append(
            {
                "scope": "clear_results_source",
                "source": "implementation",
                "unit_id": "",
                "stage": "implementation",
                "record_key": path.name,
                "path": identity,
                "declared_sha256": sha,
                "observed_sha256": sha,
                "declared_size_bytes": size,
                "observed_size_bytes": size,
                "status": "verified",
            }
        )
    return (
        {
            "sources": records,
            "source_hashing": {
                "encoding": "strict UTF-8",
                "newline_canonicalization": "CRLF to LF",
                "unicode_normalization": "none",
                "reject_bom": True,
                "reject_bare_carriage_return": True,
            },
            "software": {
                "python": platform.python_version(),
                "matplotlib": matplotlib.__version__,
                "numpy": np.__version__,
                "pandas": pd.__version__,
            },
        },
        pd.DataFrame(rows),
    )


def _verify_sensitivity_minp_against_pointwise(
    pointwise: pd.DataFrame, minp: pd.DataFrame
) -> None:
    sensitivity_p = pointwise[pointwise["estimand"] == SENSITIVITY_ESTIMAND]
    sensitivity_minp = minp[minp["estimand"] == SENSITIVITY_ESTIMAND]
    for record in sensitivity_minp.to_dict(orient="records"):
        values = sensitivity_p[
            (sensitivity_p["source"] == record["source"])
            & (sensitivity_p["unit_id"] == record["unit_id"])
        ]["mc_p_upper"].to_numpy(dtype=float)
        if len(values) != 7 or not np.isclose(
            min(values), float(record["observed_min_marginal_p_upper"])
        ):
            raise ValueError("matched-AF sensitivity pointwise and minP tables differ")


def _assert_focused_class_lineage(
    aggregate: pd.DataFrame, direct: pd.DataFrame, *, source: str
) -> None:
    keys = ["unit_id", "genotype_class", "threshold_years"]
    aggregate = aggregate.sort_values(keys, kind="stable").reset_index(drop=True)
    direct = direct.sort_values(keys, kind="stable").reset_index(drop=True)
    if (
        len(aggregate) != len(direct)
        or aggregate.duplicated(keys).any()
        or direct.duplicated(keys).any()
        or tuple(aggregate.columns) != tuple(direct.columns)
    ):
        raise ValueError(f"{source} focused aggregate lineage grid differs")
    try:
        pd.testing.assert_frame_equal(
            aggregate,
            direct,
            check_dtype=False,
            check_exact=source == "gamma_smc",
            rtol=1e-11,
            atol=5e-12 if source == "tree_truth" else 0.0,
        )
    except AssertionError as error:
        raise ValueError(
            f"{source} focused aggregate differs from unit receipts"
        ) from error


def _validate_focused_analysis_lineage(
    bundles: Mapping[str, AnalysisBundle], existing: ExistingUnitBundle
) -> dict[str, Any]:
    unit_ids = set(existing.units["unit_id"].astype(str))
    records: dict[str, Any] = {}
    for source, bundle in bundles.items():
        aggregate = _read_tsv(
            bundle.input_paths["class_summaries"],
            label=f"{source} aggregate class summaries",
        )
        aggregate = aggregate[aggregate["unit_id"].astype(str).isin(unit_ids)].copy()
        direct = existing.class_summaries[source].copy()
        if len(aggregate) != 3_080 or len(direct) != 3_080:
            raise ValueError(f"{source} focused aggregate lineage grid differs")
        _assert_focused_class_lineage(aggregate, direct, source=source)
        legacy = _existing_unit_scores(bundle, existing.units).sort_values(
            ["unit_id", "threshold_years"], kind="stable"
        )
        derived = existing.unit_scores[
            existing.unit_scores["source"] == source
        ].sort_values(["unit_id", "threshold_years"], kind="stable")
        identity_columns = [
            "unit_id",
            "simulation_class",
            "replicate_index",
            "threshold_years",
            "bank_role",
            "score_scope",
        ]
        score_columns = [
            "selection_coefficient",
            "target_allele_frequency",
            "panel_n_diploid_samples",
            "p_overall",
            "p_hom_ref",
            "p_hom_alt",
            "p_hom_alt_minus_hom_ref",
        ]
        if len(legacy) != 770 or len(derived) != 770:
            raise ValueError(f"{source} focused direct score grid differs")
        if (
            not legacy[identity_columns]
            .reset_index(drop=True)
            .equals(derived[identity_columns].reset_index(drop=True))
        ):
            raise ValueError(f"{source} focused score identities differ")
        for column in score_columns:
            if not np.allclose(
                legacy[column].to_numpy(dtype=float),
                derived[column].to_numpy(dtype=float),
                rtol=1e-11,
                atol=5e-12,
            ):
                raise ValueError(f"{source} focused score differs: {column}")
        records[source] = {
            "aggregate_class_summaries_sha256": _sha256_file(
                bundle.input_paths["class_summaries"]
            ),
            "direct_unit_class_rows": len(direct),
            "direct_unit_scores_sha256": _canonical_sha256(
                derived.to_dict(orient="records")
            ),
            "legacy_replicate_statistics_crosscheck": "matched",
        }
    return records


def _prepare_inputs(
    repo_root: Path,
    *,
    campaign_root: Path,
    gamma_analysis_dir: Path,
    truth_analysis_dir: Path,
    survival_results_dir: Path,
    two_epoch_truth_metrics: Path,
    two_epoch_gamma_metrics: Path,
    two_epoch_gamma_results: Path,
    workers: int,
) -> PreparedInputs:
    analysis_dirs = {
        "gamma_smc": gamma_analysis_dir,
        "tree_truth": truth_analysis_dir,
    }
    bundles = {
        spec.key: _load_analysis_bundle(
            repo_root, analysis_dirs[spec.key], spec, workers=workers
        )
        for spec in ANALYSIS_SOURCES
    }
    existing = _load_existing_unit_bundle(repo_root, campaign_root, workers=workers)
    survival = _load_survival_bundle(repo_root, survival_results_dir, workers=workers)
    focused_lineage = _validate_focused_analysis_lineage(bundles, existing)

    unit_scores = pd.concat(
        [existing.unit_scores, _survival_unit_scores(survival)], ignore_index=True
    )
    unit_scores = unit_scores.sort_values(
        ["source", "bank_role", "unit_id", "threshold_years"], kind="stable"
    ).reset_index(drop=True)
    if len(unit_scores) != (
        2
        * (
            SELECTED_REPLICATES
            + MATCHED_AF_NEUTRAL_REPLICATES
            + SURVIVAL_NEUTRAL_REPLICATES
        )
        * len(THRESHOLDS_YEARS)
    ):
        raise ValueError("combined unit score inventory is incomplete")

    pointwise = pd.concat(
        [
            _primary_pointwise(unit_scores),
            _sensitivity_pointwise_from_scores(unit_scores),
        ],
        ignore_index=True,
    ).sort_values(["source", "estimand", "unit_id", "threshold_years"], kind="stable")
    minp = pd.concat(
        [
            _primary_minp(unit_scores, pointwise),
            _sensitivity_minp_from_scores(unit_scores, pointwise),
        ],
        ignore_index=True,
    ).sort_values(["source", "estimand", "unit_id"], kind="stable")
    _verify_sensitivity_minp_against_pointwise(pointwise, minp)
    power = _power_from_pointwise(pointwise).sort_values(
        ["source", "estimand", "threshold_years"], kind="stable"
    )
    for bundle in bundles.values():
        _verify_sensitivity_power(bundle, power)

    endpoint_columns = [
        "unit_id",
        "bank_role",
        "simulation_class",
        "replicate_index",
        "seed",
        "selection_coefficient",
        "present_conditioning",
        "simulation_completion_sha256",
        "final_census_alt_count",
        "final_census_total_count",
        "final_census_af",
        "final_census_segregating",
        "pool_diploids",
        "final_pool_alt_count",
        "final_pool_total_count",
        "final_pool_af",
        "pool_detected",
        "pool_fixed",
        "sample_diploids",
        "sample_alt_count",
        "sample_total_count",
        "sample_af",
        "sample_detected",
        "sample_fixed",
        "panel_seed",
        "hom_ref_diploids",
        "heterozygous_diploids",
        "hom_alt_diploids",
        "mutation_origin_generations",
        "introgression_pulse_generations",
        "selection_end_generations",
        "realized_frequency_trajectory_persisted",
    ]
    endpoints = pd.concat([existing.endpoints, survival.endpoints], ignore_index=True)[
        endpoint_columns
    ].sort_values(["bank_role", "replicate_index"], kind="stable")
    if len(endpoints) != 210:
        raise ValueError("combined endpoint inventory is incomplete")

    truth_metrics, gamma_metrics, comparator_audit, comparator_contract = (
        _comparator_inputs(
            repo_root,
            truth_metrics_path=two_epoch_truth_metrics,
            gamma_metrics_path=two_epoch_gamma_metrics,
            gamma_results_path=two_epoch_gamma_results,
            workers=workers,
        )
    )
    comparison = _comparison_design(pointwise, power, truth_metrics, gamma_metrics)
    implementation, implementation_audit = _implementation_contract(repo_root)
    cache_audit = pd.concat(
        [
            *(bundle.audit for bundle in bundles.values()),
            existing.audit,
            survival.audit,
            comparator_audit,
            implementation_audit,
        ],
        ignore_index=True,
    ).sort_values(
        ["scope", "source", "unit_id", "stage", "record_key", "path"],
        kind="stable",
    )
    if set(cache_audit["status"].astype(str)) != {"verified"}:
        raise ValueError("input cache audit contains an unverified record")
    contract = {
        "schema": SCHEMA_VERSION,
        "implementation": implementation,
        "inputs": {
            "campaign_root": _path_identity(campaign_root, repo_root),
            "existing_selected_and_matched_af": existing.input_contract,
            "focused_unit_to_analysis_lineage": focused_lineage,
            "survival_neutral": survival.input_contract,
            "analysis": {
                key: bundle.input_contract for key, bundle in sorted(bundles.items())
            },
            "two_epoch_comparator": comparator_contract,
        },
        "parameters": {
            "demography_id": DEMOGRAPHY_ID,
            "source_model": SOURCE_MODEL,
            "selected_replicates": SELECTED_REPLICATES,
            "survival_neutral_replicates": SURVIVAL_NEUTRAL_REPLICATES,
            "matched_af_sensitivity_replicates": MATCHED_AF_NEUTRAL_REPLICATES,
            "sample_diploids": SAMPLE_DIPLOIDS,
            "present_han_unscaled_ne": PRESENT_HAN_UNSCALED_NE,
            "present_han_scaled_census_diploids": (PRESENT_HAN_SCALED_CENSUS_DIPLOIDS),
            "present_han_scaled_census_genomes": PRESENT_HAN_SCALED_CENSUS_GENOMES,
            "target_selected_sample_af": TARGET_AF,
            "selection_coefficient": SELECTION_COEFFICIENT,
            "primary_estimand": PRIMARY_ESTIMAND,
            "primary_reference_interpretation": (
                "selection plus endpoint-AF enrichment relative to survival-only neutral; "
                "not an exchangeability-valid selection-only p-value"
            ),
            "sensitivity_estimand": SENSITIVITY_ESTIMAND,
            "thresholds_years": list(THRESHOLDS_YEARS),
            "display_threshold_years": DISPLAY_THRESHOLD_YEARS,
            "display_threshold_role": "pulse-aligned display only",
            "primary_simulation_reference_omnibus": (
                "seven-threshold dependence-preserving simulation-reference minP"
            ),
            "p_value": "one-sided upper-tail simulation-reference with +1 correction",
            "mutation_origin_generations": MUTATION_ORIGIN_GENERATIONS,
            "introgression_pulse_generations": INTROGRESSION_PULSE_GENERATIONS,
            "selection_end_generations": SELECTION_END_GENERATIONS,
            "generation_time_years": GENERATION_TIME_YEARS,
            "workers": workers,
            "figure_size_inches": [11, 8.5],
            "png_dpi": 300,
            "pdf_media_box_points": [0, 0, 792, 612],
            "realized_selected_frequency_trajectory_persisted": False,
        },
    }
    return PreparedInputs(
        unit_scores=unit_scores,
        selected_pointwise=pointwise,
        selected_minp=minp,
        power=power,
        endpoints=endpoints,
        cache_audit=cache_audit,
        comparison=comparison,
        contract=contract,
    )


def _run_results_text(prepared: PreparedInputs) -> str:
    selected_endpoints = prepared.endpoints[
        prepared.endpoints["bank_role"] == "selected"
    ].sort_values("replicate_index")
    if len(selected_endpoints) != SELECTED_REPLICATES:
        raise ValueError("selected endpoint counts are incomplete for RUN_RESULTS")
    aa_n = _exact_integer(selected_endpoints, "hom_alt_diploids")
    rr_n = _exact_integer(selected_endpoints, "hom_ref_diploids")
    lines = [
        "# AncientEurasia 9K19: clear focal results",
        "",
        (
            "This additive bundle reused completed simulations and Gamma-SMC decodes; "
            "it did not simulate or decode any unit."
        ),
        "",
        "## Statistical contract",
        "",
        (
            "- Primary: selected overall focal P(TMRCA < x) versus the new single-copy, "
            "complete-simulated-Q5-present-Han-census-segregating neutral bank (100 units, "
            "s=0; unscaled Ne=6300, scaled census=1260 diploids/2520 genomes; no terminal "
            "AF target/band, selected-AF matching, pool-frequency gate, or panel-detection "
            "gate beyond the strict census 0<AF<1 condition). Because selected panels "
            "are exact-AF20 "
            "conditioned and this neutral bank is not, these are simulation-reference "
            "tails combining selection and endpoint-AF enrichment—not exchangeability-valid "
            "selection-only p-values."
        ),
        (
            "- Selection-isolating sensitivity: selected alt/alt-minus-ref/ref versus the "
            "existing final-AF20-matched neutral contrast."
        ),
        (
            f"- Genotype-class limitation: each selected panel has only {aa_n.min()}-"
            f"{aa_n.max()} alt/alt pairs (median {np.median(aa_n):.1f}), versus "
            f"{rr_n.min()}-{rr_n.max()} ref/ref pairs (median {np.median(rr_n):.1f}). "
            "The alt/alt curves are therefore noisy descriptive summaries, not the primary test."
        ),
        (
            "- All seven thresholds were prespecified. 50 kya is a pulse-aligned display "
            "cutoff; the primary simulation-reference minP retains whole "
            "seven-threshold vectors."
        ),
        "",
        "## Pulse-aligned display summary (50 kya)",
        "",
        (
            "| Source | Median selected score | Median neutral score | Median reference p | "
            "p<=0.05 fraction |"
        ),
        "|---|---:|---:|---:|---:|",
    ]
    for source, label in (("tree_truth", "Tree truth"), ("gamma_smc", "Gamma-SMC")):
        rows = prepared.selected_pointwise[
            (prepared.selected_pointwise["source"] == source)
            & (prepared.selected_pointwise["estimand"] == PRIMARY_ESTIMAND)
            & (
                prepared.selected_pointwise["threshold_years"]
                == DISPLAY_THRESHOLD_YEARS
            )
        ]
        power = prepared.power[
            (prepared.power["source"] == source)
            & (prepared.power["estimand"] == PRIMARY_ESTIMAND)
            & (prepared.power["threshold_years"] == DISPLAY_THRESHOLD_YEARS)
        ].iloc[0]
        lines.append(
            f"| {label} | {np.median(rows['observed_selection_score']):.6g} | "
            f"{np.median(rows['null_median']):.6g} | "
            f"{np.median(rows['mc_p_upper']):.6g} | "
            f"{power['conditional_power_raw']:.3f} |"
        )
    lines.extend(
        [
            "",
            "## Corrected two-epoch comparison",
            "",
            (
                "The attached two-epoch neutral-null histogram is tree truth: its focal "
                "upper-tail result was 1/101. The corresponding Gamma-SMC focal decode of "
                "that linked low-AF selected base was 0.0305024 versus a neutral mean of "
                "0.0234616, with 13/100 neutral exceedances and p=14/101=0.1386. It was "
                "the +/-500 kb positive-window density statistic—not the focal Gamma-SMC "
                "p-value—that reached 1/101."
            ),
            "",
            (
                "This is a genealogy-linked diagnostic, not a fully identical mutation-layer "
                "comparison: the truth folder used mu=1.29e-9, while the corresponding "
                "official Gamma-SMC v0.2 container decode used mu=1.25e-8. It is not a "
                "native-optimized decode."
            ),
            "",
            (
                "The two designs are also materially different: s=0.05 versus s=0.01, "
                "a 4.5-kya de novo allele versus 60-kya archaic standing variation, 2,000 "
                "within-diploid pairs versus 100, and 1-kb versus 10-kb Gamma output stride."
            ),
            "",
            "## Frequency-path limitation",
            "",
            (
                "No historical selected population-frequency trajectory was persisted. "
                "The timing figure is a schematic and the frequency panel shows observed "
                "uniform-pool/sample endpoints only; neither is a reconstructed trajectory."
            ),
            "",
            "## Output tables",
            "",
            (
                "`unit_scores.tsv` contains the complete seven-threshold truth/Gamma score "
                "grid. `selected_pointwise_pvalues.tsv` and `selected_minp.tsv` retain "
                "primary and sensitivity estimands separately. `input_cache_audit.tsv` "
                "lists every verified input receipt/artifact hash."
            ),
            "",
        ]
    )
    return "\n".join(lines)


def _png_dimensions(path: Path) -> tuple[int, int]:
    with path.open("rb") as stream:
        signature = stream.read(24)
    if len(signature) != 24 or signature[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError(f"PNG output is invalid: {path}")
    return struct.unpack(">II", signature[16:24])


def _pdf_geometry(path: Path) -> tuple[int, float, float]:
    payload = path.read_bytes()
    if not payload.startswith(b"%PDF-"):
        raise ValueError(f"PDF output is invalid: {path}")
    pages = len(re.findall(rb"/Type\s*/Page(?!s)\b", payload))
    boxes = re.findall(
        rb"/MediaBox\s*\[\s*([-+0-9.]+)\s+([-+0-9.]+)\s+"
        rb"([-+0-9.]+)\s+([-+0-9.]+)\s*\]",
        payload,
    )
    if pages != 1 or len(boxes) != 1:
        raise ValueError(f"PDF page contract is incompatible: {path}")
    x0, y0, x1, y1 = (float(value) for value in boxes[0])
    width = x1 - x0
    height = y1 - y0
    if not np.isclose(width, 792) or not np.isclose(height, 612):
        raise ValueError(f"PDF is not exact landscape letter size: {path}")
    return pages, width, height


def _output_record(path: Path, output_root: Path) -> dict[str, Any]:
    record: dict[str, Any] = {
        "path": path.relative_to(output_root).as_posix(),
        "sha256": _sha256_file(path),
        "size_bytes": path.stat().st_size,
    }
    if path.suffix == ".tsv":
        frame = _read_tsv(path, label="generated output")
        record["rows"] = len(frame)
        record["columns"] = list(frame.columns)
    elif path.suffix == ".png":
        width, height = _png_dimensions(path)
        record["width_px"] = width
        record["height_px"] = height
    elif path.suffix == ".pdf":
        pages, width, height = _pdf_geometry(path)
        record["page_count"] = pages
        record["width_points"] = width
        record["height_points"] = height
    return record


def _completion_interpretation() -> dict[str, str]:
    return {
        "primary": (
            "overall selected versus overall single-origin census-survival "
            "neutral simulation-reference tails; combines selection and "
            "endpoint-AF enrichment"
        ),
        "sensitivity": (
            "final-AF-matched alt/alt-minus-ref/ref selection-isolating contrast"
        ),
        "display_threshold": "50 kya is pulse-aligned display only",
        "primary_omnibus": (
            "dependence-preserving seven-threshold simulation-reference minP"
        ),
        "trajectory": "timing schematic and endpoints only; no stored path",
    }


def _validate_cached_output(
    output_dir: Path, completion: Mapping[str, Any], contract: Mapping[str, Any]
) -> dict[str, Any]:
    if set(completion) != {
        "schema",
        "status",
        "contract",
        "contract_sha256",
        "outputs",
        "interpretation",
    }:
        raise ValueError("clear-results completion fields differ")
    if (
        completion.get("schema") != SCHEMA_VERSION
        or completion.get("status") != "complete"
    ):
        raise ValueError("clear-results completion is incompatible")
    if completion.get("contract") != contract:
        raise ValueError("clear-results cache contract differs")
    contract_sha = _canonical_sha256(contract)
    if completion.get("contract_sha256") != contract_sha:
        raise ValueError("clear-results contract hash differs")
    if completion.get("interpretation") != _completion_interpretation():
        raise ValueError("clear-results interpretation contract differs")
    outputs = completion.get("outputs")
    expected = _expected_output_paths()
    if not isinstance(outputs, Mapping) or set(outputs) != expected:
        raise ValueError("clear-results output set is incompatible")
    completion_path = output_dir / "completion.json"
    if completion_path.is_symlink() or not completion_path.is_file():
        raise ValueError("clear-results completion path is incompatible")
    actual = {
        path.relative_to(output_dir).as_posix()
        for path in output_dir.rglob("*")
        if path.is_file() and path.name != "completion.json"
    }
    if actual != expected:
        raise ValueError("clear-results disk inventory is incompatible")
    for path in output_dir.rglob("*"):
        if path.is_symlink():
            raise ValueError("clear-results output inventory traverses a symlink")
        if path.is_dir() and path != output_dir:
            raise ValueError("clear-results output contains an unexpected directory")
    for relative, raw_record in outputs.items():
        if not isinstance(raw_record, Mapping) or raw_record.get("path") != relative:
            raise ValueError(f"clear-results output record differs: {relative}")
        path = output_dir / relative
        expected_record_keys = {"path", "sha256", "size_bytes"}
        if path.suffix == ".tsv":
            expected_record_keys.update({"rows", "columns"})
        elif path.suffix == ".png":
            expected_record_keys.update({"width_px", "height_px"})
        elif path.suffix == ".pdf":
            expected_record_keys.update({"page_count", "width_points", "height_points"})
        if set(raw_record) != expected_record_keys:
            raise ValueError(f"clear-results output record fields differ: {relative}")
        if raw_record.get("size_bytes") != path.stat().st_size:
            raise ValueError(f"clear-results output size differs: {relative}")
        if raw_record.get("sha256") != _sha256_file(path):
            raise ValueError(f"clear-results output hash differs: {relative}")
        if path.suffix == ".tsv":
            frame = _read_tsv(path, label="cached clear-results output")
            if raw_record.get("rows") != len(frame):
                raise ValueError(f"clear-results output rows differ: {relative}")
            if raw_record.get("columns") != list(frame.columns):
                raise ValueError(f"clear-results output columns differ: {relative}")
        elif path.suffix == ".png":
            dimensions = _png_dimensions(path)
            if dimensions != (3_300, 2_550):
                raise ValueError(f"clear-results PNG geometry differs: {relative}")
            if (
                raw_record.get("width_px"),
                raw_record.get("height_px"),
            ) != dimensions:
                raise ValueError(f"clear-results PNG record differs: {relative}")
        elif path.suffix == ".pdf":
            pages, width, height = _pdf_geometry(path)
            if (
                raw_record.get("page_count") != pages
                or raw_record.get("width_points") != width
                or raw_record.get("height_points") != height
            ):
                raise ValueError(f"clear-results PDF record differs: {relative}")
    return {
        "status": "verified_cached",
        "output_dir": str(output_dir),
        "contract_sha256": contract_sha,
        "output_count": len(outputs),
        "png_count": 4,
        "pdf_count": 4,
        "selected_replicates": SELECTED_REPLICATES,
        "survival_neutral_replicates": SURVIVAL_NEUTRAL_REPLICATES,
        "matched_af_sensitivity_replicates": MATCHED_AF_NEUTRAL_REPLICATES,
    }


def _write_results_bundle(prepared: PreparedInputs, destination: Path) -> None:
    frames = {
        "unit_scores.tsv": prepared.unit_scores,
        "selected_pointwise_pvalues.tsv": prepared.selected_pointwise,
        "selected_minp.tsv": prepared.selected_minp,
        "power_summary.tsv": prepared.power,
        "endpoint_af.tsv": prepared.endpoints,
        "input_cache_audit.tsv": prepared.cache_audit,
        "comparison_design.tsv": prepared.comparison,
    }
    for filename, frame in frames.items():
        _atomic_frame(destination / filename, frame)
    (destination / "RUN_RESULTS.md").write_text(
        _run_results_text(prepared), encoding="utf-8", newline="\n"
    )
    plot_paths: list[Path] = []
    for source in ("tree_truth", "gamma_smc"):
        plot_paths.extend(
            plot_pulse_aligned_null(
                prepared.unit_scores,
                prepared.endpoints,
                source=source,
                output_dir=destination,
            )
        )
    plot_paths.extend(
        plot_threshold_results(
            prepared.unit_scores,
            prepared.selected_pointwise,
            prepared.power,
            output_dir=destination,
        )
    )
    plot_paths.extend(
        plot_timing_endpoints_diagnostic(
            prepared.endpoints, prepared.comparison, output_dir=destination
        )
    )
    expected_plots = {
        f"{stem}.{suffix}" for stem in PLOT_STEMS for suffix in ("png", "pdf")
    }
    if {path.name for path in plot_paths} != expected_plots:
        raise ValueError("clear-results plot inventory is incomplete")


def _analyze_clear_results_impl(
    repo_root: str | Path,
    *,
    campaign_root: str | Path,
    output_dir: str | Path,
    gamma_analysis_dir: str | Path,
    truth_analysis_dir: str | Path,
    survival_results_dir: str | Path,
    two_epoch_truth_metrics: str | Path,
    two_epoch_gamma_metrics: str | Path,
    two_epoch_gamma_results: str | Path,
    workers: int,
    verify_only: bool,
) -> dict[str, Any]:
    root = Path(repo_root).resolve()
    if not root.is_dir():
        raise ValueError(f"repository root is absent: {root}")
    worker_count = _validate_workers(workers)
    paths = {
        "campaign_root": _resolve_path(root, campaign_root),
        "output_dir": _resolve_path(root, output_dir),
        "gamma_analysis_dir": _resolve_path(root, gamma_analysis_dir),
        "truth_analysis_dir": _resolve_path(root, truth_analysis_dir),
        "survival_results_dir": _resolve_path(root, survival_results_dir),
        "two_epoch_truth_metrics": _resolve_path(root, two_epoch_truth_metrics),
        "two_epoch_gamma_metrics": _resolve_path(root, two_epoch_gamma_metrics),
        "two_epoch_gamma_results": _resolve_path(root, two_epoch_gamma_results),
    }
    prepared = _prepare_inputs(
        root,
        campaign_root=paths["campaign_root"],
        gamma_analysis_dir=paths["gamma_analysis_dir"],
        truth_analysis_dir=paths["truth_analysis_dir"],
        survival_results_dir=paths["survival_results_dir"],
        two_epoch_truth_metrics=paths["two_epoch_truth_metrics"],
        two_epoch_gamma_metrics=paths["two_epoch_gamma_metrics"],
        two_epoch_gamma_results=paths["two_epoch_gamma_results"],
        workers=worker_count,
    )
    destination = paths["output_dir"]
    completion_path = destination / "completion.json"
    if completion_path.exists() or completion_path.is_symlink():
        return _validate_cached_output(
            destination, _read_json(completion_path), prepared.contract
        )
    if verify_only:
        raise ValueError("clear-results completion is absent")
    if destination.exists() and any(destination.iterdir()):
        raise ValueError(
            "clear-results output directory is nonempty without completion; "
            "preserve and inspect it before retrying"
        )
    temporary = destination.with_name(destination.name + f".tmp.{os.getpid()}")
    if temporary.exists():
        raise ValueError(
            f"clear-results temporary directory already exists: {temporary}"
        )
    temporary.parent.mkdir(parents=True, exist_ok=True)
    temporary.mkdir()
    try:
        _write_results_bundle(prepared, temporary)
        actual = {
            path.relative_to(temporary).as_posix()
            for path in temporary.rglob("*")
            if path.is_file()
        }
        expected = _expected_output_paths()
        if actual != expected:
            raise ValueError("generated clear-results inventory is incomplete")
        records = {
            relative: _output_record(temporary / relative, temporary)
            for relative in sorted(expected)
        }
        completion = {
            "schema": SCHEMA_VERSION,
            "status": "complete",
            "contract": prepared.contract,
            "contract_sha256": _canonical_sha256(prepared.contract),
            "outputs": records,
            "interpretation": _completion_interpretation(),
        }
        _atomic_json(temporary / "completion.json", completion)
        if destination.exists():
            destination.rmdir()
        os.replace(temporary, destination)
    except Exception:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    result = _validate_cached_output(
        destination, _read_json(destination / "completion.json"), prepared.contract
    )
    result["status"] = "generated"
    return result


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Analyze or verify the additive AncientEurasia 9K19 clear-results "
            "bundle; only post-processing routes are available."
        )
    )
    parser.add_argument("action", choices=("analyze", "verify"))
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--campaign-root", type=Path, default=DEFAULT_CAMPAIGN_ROOT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--gamma-analysis-dir", type=Path, default=DEFAULT_GAMMA_ANALYSIS_DIR
    )
    parser.add_argument(
        "--truth-analysis-dir", type=Path, default=DEFAULT_TRUTH_ANALYSIS_DIR
    )
    parser.add_argument(
        "--survival-results-dir", type=Path, default=DEFAULT_SURVIVAL_RESULTS_DIR
    )
    parser.add_argument(
        "--two-epoch-truth-metrics",
        type=Path,
        default=DEFAULT_TWO_EPOCH_TRUTH_METRICS,
    )
    parser.add_argument(
        "--two-epoch-gamma-metrics",
        type=Path,
        default=DEFAULT_TWO_EPOCH_GAMMA_METRICS,
    )
    parser.add_argument(
        "--two-epoch-gamma-results",
        type=Path,
        default=DEFAULT_TWO_EPOCH_GAMMA_RESULTS,
    )
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    result = analyze_clear_results(
        args.repo_root,
        campaign_root=args.campaign_root,
        output_dir=args.output_dir,
        gamma_analysis_dir=args.gamma_analysis_dir,
        truth_analysis_dir=args.truth_analysis_dir,
        survival_results_dir=args.survival_results_dir,
        two_epoch_truth_metrics=args.two_epoch_truth_metrics,
        two_epoch_gamma_metrics=args.two_epoch_gamma_metrics,
        two_epoch_gamma_results=args.two_epoch_gamma_results,
        workers=args.workers,
        verify_only=args.action == "verify",
    )
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


def analyze_clear_results(
    repo_root: str | Path,
    *,
    campaign_root: str | Path = DEFAULT_CAMPAIGN_ROOT,
    output_dir: str | Path = DEFAULT_OUTPUT_DIR,
    gamma_analysis_dir: str | Path = DEFAULT_GAMMA_ANALYSIS_DIR,
    truth_analysis_dir: str | Path = DEFAULT_TRUTH_ANALYSIS_DIR,
    survival_results_dir: str | Path = DEFAULT_SURVIVAL_RESULTS_DIR,
    two_epoch_truth_metrics: str | Path = DEFAULT_TWO_EPOCH_TRUTH_METRICS,
    two_epoch_gamma_metrics: str | Path = DEFAULT_TWO_EPOCH_GAMMA_METRICS,
    two_epoch_gamma_results: str | Path = DEFAULT_TWO_EPOCH_GAMMA_RESULTS,
    workers: int = DEFAULT_WORKERS,
    verify_only: bool = False,
) -> dict[str, Any]:
    """Build or verify the clear-results bundle.

    The full implementation is below the input/plot helpers so this public API
    remains the single analyze/verify entry point.
    """

    return _analyze_clear_results_impl(
        repo_root,
        campaign_root=campaign_root,
        output_dir=output_dir,
        gamma_analysis_dir=gamma_analysis_dir,
        truth_analysis_dir=truth_analysis_dir,
        survival_results_dir=survival_results_dir,
        two_epoch_truth_metrics=two_epoch_truth_metrics,
        two_epoch_gamma_metrics=two_epoch_gamma_metrics,
        two_epoch_gamma_results=two_epoch_gamma_results,
        workers=workers,
        verify_only=verify_only,
    )


__all__ = [
    "DISPLAY_THRESHOLD_YEARS",
    "PLOT_STEMS",
    "THRESHOLDS_YEARS",
    "analyze_clear_results",
    "main",
]


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
