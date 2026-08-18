from __future__ import annotations

import copy
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import tskit
from gamma_smc_aou import ancient_eurasia_clear_results as clear
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.text import Text


def _survival_estimand() -> dict[str, object]:
    return {
        "source_model": clear.SOURCE_MODEL,
        "mutation_origin_population": "Neanderthal",
        "mutation_origin_mode": "single_copy",
        "mutation_age_generations": 2_400,
        "pulse_generations": 2_272,
        "selection_coefficient": 0,
        "present_conditioning": "strict_segregating_in_present_han_census",
        "pool_frequency_conditioning": False,
        "terminal_af_target_or_band_conditioning": False,
        "sample_detection_conditioning": False,
        "panel_sampling": "uniform_100_of_500",
        "sample_diploids": 100,
        "candidate_pool_diploids": 500,
        "thresholds_years": list(clear.THRESHOLDS_YEARS),
    }


def _survival_focal_expectation() -> dict[str, object]:
    return {
        "schema": clear.SURVIVAL_FOCAL_IDENTITY_SCHEMA,
        "single_site_id": clear.FOCAL_SITE_ID,
        "focal_position_0based": clear.FOCAL_POSITION_BP,
        "expected_mutation_type": 1,
        "expected_selection_coeff": 0.0,
        "source_population": "Neanderthal",
        "source_subpopulation": 7,
        "requested_origin_time_generations": 2_400.0,
        "realized_origin_time_generations": 2_400.0,
        "slim_scaling_factor": 5.0,
        "slim_time_rule": "cycle - realized_origin_time / Q",
    }


def _absent_focal_identity() -> dict[str, object]:
    return {
        "schema": clear.SURVIVAL_FOCAL_IDENTITY_SCHEMA,
        "focal_position_0based": clear.FOCAL_POSITION_BP,
        "expected_mutation_type": 1,
        "expected_source_population": "Neanderthal",
        "expected_source_subpopulation": 7,
        "expected_origin_time_generations": 2_400.0,
        "observed_focal_site_count": 0,
        "observed_tskit_mutation_count": 0,
        "observed_metadata_entry_count": 0,
        "observed_mutation_types": [],
        "status": "selected_focal_mutation_absent",
        "validation": (
            "causal_site_absent_from_observational_pool; exact census survival "
            "retained in SLiM metadata"
        ),
    }


def _valid_focal_identity(*, cycle: int = 500) -> dict[str, object]:
    return {
        "schema": clear.SURVIVAL_FOCAL_IDENTITY_SCHEMA,
        "focal_position_0based": clear.FOCAL_POSITION_BP,
        "expected_mutation_type": 1,
        "expected_source_population": "Neanderthal",
        "expected_source_subpopulation": 7,
        "expected_origin_time_generations": 2_400.0,
        "observed_focal_site_count": 1,
        "observed_tskit_mutation_count": 1,
        "observed_metadata_entry_count": 1,
        "observed_mutation_types": [1],
        "status": "valid",
        "observed_mutation_type": 1,
        "observed_source_subpopulation": 7,
        "observed_origin_time_generations": 2_400.0,
        "observed_selection_coeff": 0.0,
        "observed_slim_time": cycle - 480,
        "observed_slim_cycle": cycle,
        "observed_slim_scaling_factor": 5.0,
        "observed_ancestral_state": "A",
        "observed_derived_state": "C",
        "observed_nucleotide": 1,
    }


def _focal_overlay_receipt() -> dict[str, object]:
    contract = copy.deepcopy(clear.SURVIVAL_FOCAL_OVERLAY_PATCH)
    tick = 29_665.0
    calls = []
    for index, expected in enumerate(contract["expected_call_inventory"], start=1):
        start = tick if expected["start_time"] == "finite_shared_slim_tick" else None
        end = tick if expected["end_time"] == "finite_shared_slim_tick" else None
        original = clear._survival_overlay_map_record(
            expected["original_positions"], expected["original_rates"]
        )
        passed = clear._survival_overlay_map_record(
            expected["passed_positions"], expected["passed_rates"]
        )
        calls.append(
            {
                "call_index": index,
                "layer": expected["layer"],
                "rate_kind": "msprime.RateMap",
                "model_class": expected["model_class"],
                "model_type": expected["model_type"],
                "model_next_id": index * 100,
                "start_time": start,
                "end_time": end,
                "random_seed": index * 1_000,
                "keep": True,
                "discrete_genome": None,
                "record_provenance": True,
                "non_rate_arguments_preserved_by_identity": True,
                "model_identity_preserved": True,
                "rate_argument_replaced": expected["patched"],
                "positive_focal_rate_before": expected["patched"],
                "positive_focal_rate_after": False,
                "original_focal_rate": original["focal_integral"],
                "passed_focal_rate": passed["focal_integral"],
                "original_rate_map": original,
                "passed_rate_map": passed,
                "call_succeeded": True,
            }
        )
    return {
        "schema": clear.SURVIVAL_FOCAL_OVERLAY_PATCH_SCHEMA,
        "contract": contract,
        "contract_sha256": clear._canonical_sha256(contract),
        "intercepted_call_count": 2,
        "patched_call_count": 1,
        "all_positive_focal_rates_masked": True,
        "callable_restored": True,
        "calls": calls,
    }


def _individual_location_receipt(tree: tskit.TreeSequence) -> dict[str, object]:
    contract = copy.deepcopy(clear.SURVIVAL_INDIVIDUAL_LOCATION_CANONICALIZATION)
    state = clear._individual_location_state(tree.dump_tables())
    return {
        "schema": clear.SURVIVAL_INDIVIDUAL_LOCATION_SCHEMA,
        "contract": contract,
        "contract_sha256": clear._canonical_sha256(contract),
        "before": copy.deepcopy(state),
        "after": copy.deepcopy(state),
        "stored_tree": copy.deepcopy(state),
        "proof": {
            "location_offsets_and_widths_unchanged": True,
            "all_nonlocation_individual_columns_unchanged": True,
            "all_other_tables_ts_metadata_reference_unchanged": True,
            "only_individual_location_values_canonicalized": True,
            "stored_tree_state_bound_at_canonicalization": True,
        },
    }


def _survival_scores() -> pd.DataFrame:
    rows = []
    for index in range(100):
        unit_id = f"ancient_eurasia_single_origin_survival__neutral__rep{index:03d}"
        pool_alt = 0 if index == 0 else 1_000 if index == 1 else index + 3
        sample_alt = 0 if index % 29 == 0 else min(200, index // 2 + 1)
        for source_index, source in enumerate(("tree_truth", "gamma_smc")):
            for threshold_index, threshold in enumerate(clear.THRESHOLDS_YEARS):
                probability = min(
                    0.99,
                    0.002
                    + 0.0003 * index
                    + 0.006 * threshold_index
                    + 0.001 * source_index,
                )
                rows.append(
                    {
                        "unit_id": unit_id,
                        "replicate_index": index,
                        "seed": clear._survival_seed(unit_id),
                        "source": source,
                        "threshold_years": threshold,
                        "overall_focal_p_tmrca_lt_threshold": probability,
                        "final_pool_alt_count": pool_alt,
                        "final_pool_total_count": 1_000,
                        "final_pool_af": pool_alt / 1_000,
                        "sample_alt_count": sample_alt,
                        "sample_total_count": 200,
                        "sample_af": sample_alt / 200,
                        "sample_detected": sample_alt > 0,
                    }
                )
    return pd.DataFrame(rows)


def _survival_endpoints() -> pd.DataFrame:
    score_rows = (
        _survival_scores()
        .drop_duplicates(["unit_id"])
        .sort_values("unit_id", kind="stable")
    )
    rows = []
    for row in score_rows.to_dict(orient="records"):
        pool_alt = int(row["final_pool_alt_count"])
        sample_alt = int(row["sample_alt_count"])
        rows.append(
            {
                "unit_id": row["unit_id"],
                "replicate_index": int(row["replicate_index"]),
                "seed": int(row["seed"]),
                "simulation_completion_sha256": "a" * 64,
                "final_census_alt_count": int(row["replicate_index"]) + 1,
                "final_census_total_count": 2_520,
                "final_census_af": (int(row["replicate_index"]) + 1) / 2_520,
                "final_census_segregating": True,
                "pool_diploids": 500,
                "pool_alt_count": pool_alt,
                "pool_total_count": 1_000,
                "pool_af": pool_alt / 1_000,
                "pool_detected": pool_alt > 0,
                "pool_fixed": pool_alt == 1_000,
                "sample_diploids": 100,
                "sample_alt_count": sample_alt,
                "sample_total_count": 200,
                "sample_af": sample_alt / 200,
                "sample_detected": sample_alt > 0,
                "sample_fixed": sample_alt == 200,
                "panel_seed": clear._survival_seed(f"{row['unit_id']}:panel"),
            }
        )
    return pd.DataFrame(rows, columns=clear.SURVIVAL_ENDPOINT_COLUMNS)


def _survival_execution() -> pd.DataFrame:
    endpoints = _survival_endpoints()
    rows = []
    for row in endpoints.to_dict(orient="records"):
        rows.append(
            {
                "unit_id": row["unit_id"],
                "replicate_index": row["replicate_index"],
                "seed": row["seed"],
                "panel_seed": row["panel_seed"],
                "scenario": "ancient_eurasia_single_origin_survival_neutral",
                "selection_coefficient": 0.0,
                "mutation_origin_population": "Neanderthal",
                "mutation_age_generations": 2_400,
                "pulse_generations": 2_272,
                "present_conditioning": "strict_segregating_in_present_han_census",
                "pool_frequency_conditioning": False,
                "terminal_af_target_or_band_conditioning": False,
                "sample_detection_conditioning": False,
                "pool_diploids": 500,
                "sample_diploids": 100,
            }
        )
    return pd.DataFrame(rows, columns=clear.SURVIVAL_EXECUTION_COLUMNS)


def _output_record(path: Path, *, relative: str, rows: int | None) -> dict[str, object]:
    return {
        "path": relative,
        "sha256": clear._sha256_file(path),
        "size_bytes": path.stat().st_size,
        "rows": rows,
    }


def _focused_receipt_fixture(
    repo: Path,
) -> tuple[Path, dict[str, object], Path, Path]:
    campaign = repo / "focused_selection_EAS_sim"
    unit_id = "ancient_eurasia_han_introgression__af20__neutral__rep001"
    unit_dir = campaign / "work" / unit_id
    decode_dir = unit_dir / "decoded"
    (unit_dir / "pairs").mkdir(parents=True, exist_ok=True)
    decode_dir.mkdir(exist_ok=True)
    unit = {
        "unit_id": unit_id,
        "demography_id": clear.DEMOGRAPHY_ID,
        "demography_kind": "introgression",
        "population": "Han",
        "source_model": clear.SOURCE_MODEL,
        "selection_origin": "archaic_specific_introgressed_standing_variation",
        "simulation_class": "neutral",
        "selection_coefficient": 0.0,
        "target_allele_frequency": clear.TARGET_AF,
        "population_af_lower": 0.175,
        "population_af_upper": 0.225,
        "exact_sample_alt_count": 40,
        "sample_diploids": clear.SAMPLE_DIPLOIDS,
        "replicate_index": 1,
        "seed": 1_607_718_386,
        "sequence_length_bp": 10_000_000,
        "focal_position_bp": clear.FOCAL_POSITION_BP,
    }
    simulation_paths = {
        "tree": unit_dir / "simulation.trees",
        "pair_table": unit_dir / "sample_manifest.tsv",
        "overall_pairs": unit_dir / "pairs/overall.pairs.tsv",
        "truth_profiles": unit_dir / "truth_profiles.tsv.gz",
        "truth_class_summaries": unit_dir / "truth_class_summaries.tsv",
    }
    allele_count = np.r_[np.zeros(65, dtype=int), np.ones(30, dtype=int), np.full(5, 2)]
    diploid_index = np.arange(100)
    manifest = pd.DataFrame(
        {
            "vcf_diploid_index": diploid_index,
            "tree_sequence_individual_id": diploid_index,
            "sample_node_0": 2 * diploid_index,
            "sample_node_1": 2 * diploid_index + 1,
            "gamma_smc_haplotype_0": 2 * diploid_index,
            "gamma_smc_haplotype_1": 2 * diploid_index + 1,
            "focal_position_0based": float(clear.FOCAL_POSITION_BP),
            "focal_selected_allele_index": 1,
            "focal_selected_allele_count": allele_count,
            "genotype_class": np.asarray(
                ["hom_ref", "heterozygous", "hom_alt"], dtype=object
            )[allele_count],
        },
        columns=clear.SURVIVAL_SAMPLE_MANIFEST_COLUMNS,
    )
    manifest.to_csv(simulation_paths["pair_table"], sep="\t", index=False)
    pair_lines = [
        "# Gamma-SMC explicit haplotype pairs",
        "# genotype_class\toverall",
        "# n_pairs\t100",
        *(
            f"{left}\t{right}"
            for left, right in zip(
                2 * diploid_index, 2 * diploid_index + 1, strict=True
            )
        ),
    ]
    simulation_paths["overall_pairs"].write_text(
        "\n".join(pair_lines) + "\n", encoding="utf-8"
    )
    for key in ("tree", "truth_profiles", "truth_class_summaries"):
        simulation_paths[key].write_text(f"{key}\n", encoding="utf-8")
    simulation_outputs = {
        key: {
            "path": path.relative_to(unit_dir).as_posix(),
            "sha256": clear._sha256_file(path),
        }
        for key, path in simulation_paths.items()
    }
    focal_contract = {
        "exact_sample_alt_count": 40,
        "identity": {
            "ancestral_state": "A",
            "archaic_specific_no_ils_by_construction": True,
            "derived_state": "G",
            "expected_origin_time_generations": 2_400.0,
            "expected_source_population": "Neanderthal",
            "focal_position_0based": clear.FOCAL_POSITION_BP,
            "human_neanderthal_split_generations": 27_840.0,
            "pulse_route": {
                "destination_population": "Neanderthal",
                "minimum_descendant_lineages": 1,
                "source_population": "Loschbour",
                "time_generations": 2_272.0,
            },
            "schema": "gamma-smc.neutral-focal-identity/v1",
            "single_nonrecurrent_mutation": True,
        },
        "minimum_heterozygous_diploids": 1,
        "minimum_hom_alt_diploids": 2,
        "minimum_hom_ref_diploids": 2,
        "sample_diploids": 100,
        "schema": "gamma-smc.focused-focal-semantics/v1",
        "simulation_class": "neutral",
        "tree_manifest_exact_match": True,
    }
    parameters = {
        "candidate_pool_diploids": 500,
        "candidate_pool_frequency_gate": {
            "interval_source": (
                "unit.population_af_lower_to_population_af_upper_inclusive"
            ),
            "order": "before_exact_sample_panel",
            "schema": "gamma-smc.candidate-pool-frequency-gate/v1",
            "scope": "selected_and_neutral",
        },
        "focal_semantic_validation": focal_contract,
        "han_selection_calibration": None,
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
    simulation_contract = {
        "schema": clear.FOCUSED_SIMULATION_SCHEMA_VERSION,
        "unit": unit,
        "parameters": parameters,
        "software": copy.deepcopy(clear.FOCUSED_SIMULATION_SOFTWARE),
        "implementation": {
            "eas_resource_sha256": (
                "a5bdfd843629d48050cd56a84da4d70e356d0339017d3c1f67be646d727db928"
            ),
            "slim": None,
            "sources": copy.deepcopy(clear.FOCUSED_NEUTRAL_SIMULATION_SOURCE_HASHES),
        },
    }
    observed_identity = {
        "ancestral_state": "A",
        "archaic_specific_no_ils_by_construction": True,
        "derived_state": "G",
        "focal_position_0based": clear.FOCAL_POSITION_BP,
        "mutation_metadata_empty": True,
        "observed_source_population": "Neanderthal",
        "origin_time_generations": 2_400.0,
        "pulse_route_evidence": {
            "descendant_lineage_nodes": [1],
            "destination_population": "Neanderthal",
            "n_descendant_lineages": 1,
            "source_population": "Loschbour",
            "time_generations": 2_272.0,
        },
        "schema": "gamma-smc.neutral-focal-identity/v1",
        "status": "valid",
    }
    pool_counts = {"hom_ref": 324, "heterozygous": 150, "hom_alt": 26}
    sample_counts = {"hom_ref": 65, "heterozygous": 30, "hom_alt": 5}
    gate = {
        "candidate_pool_diploids": 500,
        "candidate_pool_genotype_counts": pool_counts,
        "candidate_pool_haplotypes": 1_000,
        "gate_order": "before_exact_sample_panel",
        "maximum_alt_count_inclusive": 225,
        "minimum_alt_count_inclusive": 175,
        "observed_af": 0.202,
        "observed_alt_count": 202,
        "passed": True,
        "population_af_lower_inclusive": 0.175,
        "population_af_upper_inclusive": 0.225,
        "schema": "gamma-smc.candidate-pool-frequency-gate/v1",
        "status": "pass",
    }
    focal_provenance = {
        "genotype_counts": sample_counts,
        "identity": observed_identity,
        "sample_af": 0.2,
        "sample_alt_count": 40,
        "sample_diploids": 100,
        "schema": "gamma-smc.focused-focal-semantics/v1",
        "status": "valid",
        "tree_manifest_exact_match": True,
    }
    simulation_provenance = {
        "accepted_attempt_zero_based": 2,
        "accepted_batch_zero_based": 0,
        "ancestry_batch_size": 4,
        "ancestry_seed": 1,
        "archaic_specific_no_ils": True,
        "ascertainment_length_bp": 20_000_000,
        "candidate_index": (
            "Loschbour_to_Neanderthal_pulse_migrations_then_ascend_to_"
            "branch_crossing_2400_generations"
        ),
        "candidate_pool_diploids": 500,
        "candidate_pool_frequency_gate": gate,
        "candidate_seed": 2,
        "candidate_weight": "genomic_span_bp_at_fixed_2400_generation_origin",
        "coalescent_approximation": "SMC-prime after 200 generations",
        "engine": "msprime",
        "entry_route": (
            "Loschbour_to_Neanderthal_backward_mass_migration_at_2272_generations"
        ),
        "focal_semantic_validation": focal_provenance,
        "matched_selected_neutral_ascertainment": (
            "same_candidate_pool_size_and_AF_band_then_same_exact_k_panel_sampler"
        ),
        "model": "DTWF200ThenSmcPrimeApproxCoalescent",
        "model_keywords": ["dtwf", "smc_prime"],
        "msprime_version": "1.4.2",
        "mutation_origin_population": "Neanderthal",
        "mutation_seed": 3,
        "mutation_time_generations_ago": 2_400.0,
        "n_eligible_branch_rectangles_in_batch": 1,
        "n_introgressed_descendant_lineages": 1,
        "neutral_focal_identity": observed_identity,
        "neutral_origin_matching": (
            "archaic_population_branch_with_one_or_more_pulse_lineages; "
            "nonzero standing variation at pulse"
        ),
        "panel_seed": 4,
        "pre_crop_focal_position": 5_000_000.0,
        "recent_dtwf_duration_generations": 200.0,
        "sample_panel": {
            "candidate_pool_diploids": 500,
            "candidate_pool_genotype_counts": pool_counts,
            "chosen_sample_genotype_counts": sample_counts,
            "eligible_genotype_count_triples": 1,
            "panel_seed": 4,
            "selection_method": (
                "uniform_over_all_size_N_exact_alt_count_subsets_passing_genotype_qc"
            ),
        },
        "window_length_bp": 10_000_000,
    }
    simulation = {
        "schema": clear.FOCUSED_SIMULATION_SCHEMA_VERSION,
        "status": "complete",
        "contract": simulation_contract,
        "contract_sha256": clear._canonical_sha256(simulation_contract),
        "elapsed_seconds": 1.0,
        "attempts_completed": 4,
        "search_budget": {
            "neutral_max_attempts": 100,
            "selected_cumulative_timeout_seconds": 2_700.0,
            "selected_draw_timeout_seconds": 300.0,
            "selected_max_draws": 50,
        },
        "sample_alt_count": 40,
        "sample_af": 0.2,
        "genotype_counts": sample_counts,
        "provenance": simulation_provenance,
        "outputs": simulation_outputs,
    }
    simulation_contract_path = unit_dir / "simulation_contract.json"
    simulation_path = unit_dir / "simulation_complete.json"
    simulation_contract_path.write_text(
        json.dumps(simulation_contract, sort_keys=True) + "\n", encoding="utf-8"
    )
    simulation_path.write_text(
        json.dumps(simulation, sort_keys=True) + "\n", encoding="utf-8"
    )
    posterior = decode_dir / "posterior.zst"
    metadata = decode_dir / "posterior.zst.meta"
    posterior.write_text("posterior\n", encoding="utf-8")
    metadata.write_text("metadata\n", encoding="utf-8")
    decode_paths = {
        "overall_summary": decode_dir / "overall.summary.tsv",
        "class_summaries": decode_dir / "class_summaries.tsv",
        "pair_summaries": decode_dir / "pair_summaries.tsv.gz",
        "spatial_profiles": decode_dir / "spatial_profiles.tsv.gz",
        "run_manifest": decode_dir / "overall.summary.tsv.run.json",
        "decoder_stdout": decode_dir / "decoder.stdout.zst",
        "decoder_stderr": decode_dir / "decoder.stderr.zst",
    }
    for key, path in decode_paths.items():
        path.write_text(f"{key}\n", encoding="utf-8")
    decode_outputs = {}
    for key, path in decode_paths.items():
        record = {
            "path": path.relative_to(decode_dir).as_posix(),
            "sha256": clear._sha256_file(path),
            "size_bytes": path.stat().st_size,
        }
        if key in clear.FOCUSED_DECODE_TABULAR_OUTPUT_ROWS:
            record.update(
                {
                    "rows": clear.FOCUSED_DECODE_TABULAR_OUTPUT_ROWS[key],
                    "columns": clear.FOCUSED_DECODE_TABULAR_OUTPUT_COLUMNS[key],
                }
            )
        decode_outputs[key] = record
    decode_provenance = {
        "decoder_path": clear.EXPECTED_GAMMA_BINARY_IDENTITY,
        "decoder_sha256": clear.EXPECTED_GAMMA_BINARY_SHA256,
        "implementation_sources": copy.deepcopy(clear.FOCUSED_DECODE_SOURCE_HASHES),
        "overall_pairs_sha256": simulation_outputs["overall_pairs"]["sha256"],
        "pair_table_sha256": simulation_outputs["pair_table"]["sha256"],
        "settings": copy.deepcopy(clear.FOCUSED_DECODER_SETTINGS),
        "simulation_completion_sha256": clear._sha256_file(simulation_path),
        "tree_sha256": simulation_outputs["tree"]["sha256"],
    }
    decode_unit = {
        **unit,
        "present_ne": clear.PRESENT_HAN_UNSCALED_NE,
        "decode_provenance": decode_provenance,
    }
    decode_contract = {
        "static": {
            "schema": clear.FOCUSED_DECODE_SCHEMA_VERSION,
            "unit_record": decode_unit,
            "parameters": copy.deepcopy(clear.FOCUSED_DECODE_PARAMETERS),
            "pair_table": {
                "path": str(simulation_paths["pair_table"].resolve()),
                "sha256": simulation_outputs["pair_table"]["sha256"],
            },
            "implementation": {
                "postprocessor_path": str(
                    (
                        repo / "python/gamma_smc_aou/focused_selection_decode.py"
                    ).resolve()
                ),
                "postprocessor_sha256": clear.FOCUSED_POSTPROCESSOR_SHA256,
            },
            "auxiliary_output_paths": {
                "decoder_stderr": str(decode_paths["decoder_stderr"].resolve()),
                "decoder_stdout": str(decode_paths["decoder_stdout"].resolve()),
                "overall_summary": str(decode_paths["overall_summary"].resolve()),
                "run_manifest": str(decode_paths["run_manifest"].resolve()),
            },
        },
        "raw_inputs": {
            "posterior_path": str(posterior.resolve()),
            "posterior_sha256": clear._sha256_file(posterior),
            "metadata_path": str(metadata.resolve()),
            "metadata_sha256": clear._sha256_file(metadata),
        },
    }
    decode = {
        "schema": clear.FOCUSED_DECODE_SCHEMA_VERSION,
        "status": "complete",
        "contract": decode_contract,
        "contract_sha256": clear._canonical_sha256(decode_contract),
        "software": copy.deepcopy(clear.FOCUSED_DECODE_SOFTWARE),
        "posterior_scale": copy.deepcopy(clear.FOCUSED_DECODE_POSTERIOR_SCALE),
        "posterior_dimensions": copy.deepcopy(
            clear.FOCUSED_DECODE_POSTERIOR_DIMENSIONS
        ),
        "focal_center": copy.deepcopy(clear.FOCUSED_DECODE_FOCAL_CENTER),
        "outputs": decode_outputs,
    }
    decode_path = decode_dir / "completion.json"
    decode_path.write_text(json.dumps(decode, sort_keys=True) + "\n", encoding="utf-8")
    return campaign, unit, simulation_path, decode_path


def _survival_event_schedule() -> list[dict[str, object]]:
    condition = "ConditionOnAlleleFrequency"
    exact = "exact_generation"
    after = "generation_after_in_forward_time"
    common = {"single_site_id": clear.FOCAL_SITE_ID}
    return [
        {
            **common,
            "event_type": "DrawMutation",
            "population": "Neanderthal",
            "time": {"generations_ago": 2_400.0, "semantics": exact},
        },
        {
            **common,
            "event_type": condition,
            "population": "Neanderthal",
            "start_time": {"generations_ago": 2_400.0, "semantics": after},
            "end_time": {"generations_ago": 2_273.0, "semantics": exact},
            "operator": ">",
            "allele_frequency": 0.0,
        },
        {
            **common,
            "event_type": condition,
            "population": "Neanderthal",
            "start_time": {"generations_ago": 2_273.0, "semantics": exact},
            "end_time": {"generations_ago": 2_273.0, "semantics": exact},
            "operator": ">=",
            "allele_frequency": 1e-9,
        },
        {
            **common,
            "event_type": condition,
            "population": "Loschbour",
            "start_time": {"generations_ago": 2_272.0, "semantics": after},
            "end_time": {"generations_ago": 2_015.0, "semantics": exact},
            "operator": ">",
            "allele_frequency": 0.0,
        },
        {
            **common,
            "event_type": condition,
            "population": "Han",
            "start_time": {"generations_ago": 2_015.0, "semantics": after},
            "end_time": {"generations_ago": 0.0, "semantics": exact},
            "operator": ">",
            "allele_frequency": 0.0,
        },
        {
            **common,
            "event_type": condition,
            "population": "Han",
            "start_time": {"generations_ago": 0.0, "semantics": exact},
            "end_time": {"generations_ago": 0.0, "semantics": exact},
            "operator": "<",
            "allele_frequency": 1.0,
        },
    ]


def _combined_scores() -> pd.DataFrame:
    rows = []
    for source_index, source in enumerate(("tree_truth", "gamma_smc")):
        for index in range(100):
            for threshold_index, threshold in enumerate(clear.THRESHOLDS_YEARS):
                base = 0.01 + 0.004 * threshold_index + 0.0002 * index
                rows.append(
                    {
                        "source": source,
                        "unit_id": f"survival_{index:03d}",
                        "simulation_class": "neutral",
                        "selection_coefficient": 0.0,
                        "target_allele_frequency": np.nan,
                        "replicate_index": index,
                        "panel_n_diploid_samples": 100,
                        "threshold_years": threshold,
                        "p_overall": base + 0.0005 * source_index,
                        "p_hom_ref": np.nan,
                        "p_hom_alt": np.nan,
                        "p_hom_alt_minus_hom_ref": np.nan,
                        "bank_role": "survival_only_primary",
                        "score_scope": "focal_nearest_within_diploid",
                    }
                )
        for index in range(100):
            for threshold_index, threshold in enumerate(clear.THRESHOLDS_YEARS):
                rr = 0.02 + 0.003 * threshold_index + 0.0001 * index
                delta = -0.03 + 0.0006 * index + 0.004 * threshold_index
                rows.append(
                    {
                        "source": source,
                        "unit_id": f"matched_{index:03d}",
                        "simulation_class": "neutral",
                        "selection_coefficient": 0.0,
                        "target_allele_frequency": 0.2,
                        "replicate_index": index,
                        "panel_n_diploid_samples": 100,
                        "threshold_years": threshold,
                        "p_overall": rr + 0.01,
                        "p_hom_ref": rr,
                        "p_hom_alt": rr + delta,
                        "p_hom_alt_minus_hom_ref": delta,
                        "bank_role": "matched_af_sensitivity",
                        "score_scope": "focal_nearest_within_diploid",
                    }
                )
        for index in range(10):
            for threshold_index, threshold in enumerate(clear.THRESHOLDS_YEARS):
                rr = 0.025 + 0.004 * threshold_index + 0.0005 * index
                delta = 0.02 + 0.009 * threshold_index + 0.004 * index
                rows.append(
                    {
                        "source": source,
                        "unit_id": f"selected_{index:03d}",
                        "simulation_class": "selected",
                        "selection_coefficient": 0.01,
                        "target_allele_frequency": 0.2,
                        "replicate_index": index,
                        "panel_n_diploid_samples": 100,
                        "threshold_years": threshold,
                        "p_overall": 0.04
                        + 0.008 * threshold_index
                        + 0.003 * index
                        + 0.002 * source_index,
                        "p_hom_ref": rr,
                        "p_hom_alt": rr + delta,
                        "p_hom_alt_minus_hom_ref": delta,
                        "bank_role": "selected",
                        "score_scope": "focal_nearest_within_diploid",
                    }
                )
    return pd.DataFrame(rows)


def _endpoints() -> pd.DataFrame:
    rows = []
    for role, count, af in (
        ("selected", 10, 0.2),
        ("survival_only_primary", 100, 0.08),
        ("matched_af_sensitivity", 100, 0.2),
    ):
        for index in range(count):
            rows.append(
                {
                    "unit_id": f"{role}_{index:03d}",
                    "bank_role": role,
                    "replicate_index": index,
                    "final_pool_af": min(1, af + 0.001 * (index % 7)),
                    "sample_af": af,
                    "hom_ref_diploids": 63 + index % 4
                    if role == "selected"
                    else np.nan,
                    "hom_alt_diploids": 3 + index % 4 if role == "selected" else np.nan,
                }
            )
    return pd.DataFrame(rows)


def _comparison(pointwise: pd.DataFrame, power: pd.DataFrame) -> pd.DataFrame:
    truth = {
        "selection_coefficient": 0.05,
        "variant_age_years": 4_500,
        "neutral_replicates": 100,
        "sample_diploids": 2_000,
        "mutation_rate": 1.29e-9,
        "neutral_observed_mean_fraction_recent": 0.006555,
        "selected_rejection": {
            "population_allele_frequency": 0.102125,
            "sample_allele_frequency": 0.10275,
            "observed_fraction_recent": 0.0185,
            "neutral_exceedances": 0,
            "monte_carlo_p_upper": 1 / 101,
        },
    }
    gamma = {
        "neutral_replicates": 100,
        "sample_diploids": 2_000,
        "within_individual_pairs": 2_000,
        "stride_bp": 1_000,
        "threshold_years": 4_500,
        "center_observed_mean_p_recent": 0.0305023914,
        "center_null_mean_p_recent": 0.0234615533,
        "center_neutral_exceedances": 13,
        "center_monte_carlo_p_upper": 14 / 101,
        "scan_density_calibration": {
            "local": {
                "selected_positive_windows": 309,
                "null_mean_positive_windows": 40.04,
                "monte_carlo_p_upper": 1 / 101,
            }
        },
    }
    return clear._comparison_design(pointwise, power, truth, gamma)


def test_survival_score_schema_accepts_observational_pool_extremes_and_nondetection():
    frame = clear._validate_survival_scores(_survival_scores())

    assert len(frame) == 1_400
    assert frame["final_pool_alt_count"].min() == 0
    assert frame["final_pool_alt_count"].max() == 1_000
    assert (~frame["sample_detected"]).any()

    broken = _survival_scores()
    mask = (
        (broken["replicate_index"] == 3)
        & (broken["source"] == "tree_truth")
        & (broken["threshold_years"] == 4_500)
    )
    broken.loc[mask, "overall_focal_p_tmrca_lt_threshold"] = 0
    with pytest.raises(ValueError, match="nondecreasing"):
        clear._validate_survival_scores(broken)


def test_focused_v5_v1_receipts_remain_separate_and_cross_bound(tmp_path):
    repo = tmp_path / "repo"
    campaign, unit, simulation_path, decode_path = _focused_receipt_fixture(repo)

    checks, binding, _ = clear._unit_receipt_checks(repo, campaign, unit)
    audit = clear._audit_file_checks(checks, workers=2)
    assert len(checks) == 17
    assert set(audit["status"]) == {"verified"}
    assert binding["simulation_completion_sha256"] == clear._sha256_file(
        simulation_path
    )

    decode = json.loads(decode_path.read_text(encoding="utf-8"))
    decode["contract"]["static"]["unit_record"]["decode_provenance"]["tree_sha256"] = (
        "b" * 64
    )
    decode["contract_sha256"] = clear._canonical_sha256(decode["contract"])
    decode_path.write_text(json.dumps(decode) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="decode provenance binding"):
        clear._unit_receipt_checks(repo, campaign, unit)

    _, _, simulation_path, _ = _focused_receipt_fixture(repo)
    standalone_path = simulation_path.with_name("simulation_contract.json")
    standalone = json.loads(standalone_path.read_text(encoding="utf-8"))
    standalone["parameters"] = {"tampered": True}
    standalone_path.write_text(json.dumps(standalone) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="standalone simulation contract differs"):
        clear._unit_receipt_checks(repo, campaign, unit)

    _, _, simulation_path, _ = _focused_receipt_fixture(repo)
    simulation = json.loads(simulation_path.read_text(encoding="utf-8"))
    simulation["individual_location_canonicalization"] = {}
    simulation_path.write_text(json.dumps(simulation) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="focused simulation receipt fields differ"):
        clear._unit_receipt_checks(repo, campaign, unit)


@pytest.mark.parametrize(
    ("tamper", "error"),
    [
        ("simulation_completion_sha256", "decode provenance binding"),
        ("tree_sha256", "decode provenance binding"),
        ("pair_table_sha256", "decode provenance binding"),
        ("overall_pairs_sha256", "decode provenance binding"),
        ("pair_table_static_path", "pair-table binding"),
        ("pair_table_static_sha256", "pair-table binding"),
        ("posterior_path", "decode raw path binding"),
        ("metadata_path", "decode raw path binding"),
        ("posterior_sha256", "input hash differs"),
        ("metadata_sha256", "input hash differs"),
    ],
)
def test_focused_decode_input_bindings_reject_tamper(tmp_path, tamper, error):
    repo = tmp_path / "repo"
    campaign, unit, _, decode_path = _focused_receipt_fixture(repo)
    decode = json.loads(decode_path.read_text(encoding="utf-8"))
    static = decode["contract"]["static"]
    raw = decode["contract"]["raw_inputs"]
    provenance = static["unit_record"]["decode_provenance"]
    if tamper in {
        "simulation_completion_sha256",
        "tree_sha256",
        "pair_table_sha256",
        "overall_pairs_sha256",
    }:
        provenance[tamper] = "b" * 64
    elif tamper == "pair_table_static_path":
        static["pair_table"]["path"] = str(
            (decode_path.parent / "decoder.stdout.zst").resolve()
        )
    elif tamper == "pair_table_static_sha256":
        static["pair_table"]["sha256"] = "b" * 64
    elif tamper.endswith("_path"):
        raw[tamper] = str((decode_path.parent / "decoder.stdout.zst").resolve())
    else:
        raw[tamper] = "b" * 64
    decode["contract_sha256"] = clear._canonical_sha256(decode["contract"])
    decode_path.write_text(json.dumps(decode) + "\n", encoding="utf-8")

    if tamper.endswith("_sha256") and tamper.startswith(("posterior", "metadata")):
        checks, _, _ = clear._unit_receipt_checks(repo, campaign, unit)
        with pytest.raises(ValueError, match=error):
            clear._audit_file_checks(checks, workers=2)
    else:
        with pytest.raises(ValueError, match=error):
            clear._unit_receipt_checks(repo, campaign, unit)


def test_focused_manifest_rejects_arbitrary_pair_selector(tmp_path):
    repo = tmp_path / "repo"
    campaign, unit, _, _ = _focused_receipt_fixture(repo)
    pairs = campaign / "work" / unit["unit_id"] / "pairs/overall.pairs.tsv"
    lines = pairs.read_text(encoding="utf-8").splitlines()
    lines[3] = "0\t2"
    pairs.write_text("\n".join(lines) + "\n", encoding="utf-8")

    with pytest.raises(ValueError, match="overall pairs differ from manifest"):
        clear._unit_receipt_checks(repo, campaign, unit)


def test_analysis_inputs_are_live_hash_checked_and_null_paired(tmp_path):
    repo = tmp_path / "repo"
    relative_paths = {
        "class_summaries": (
            "focused_selection_EAS_sim/results/combined_class_summaries.tsv.gz"
        ),
        "pair_summaries": (
            "focused_selection_EAS_sim/results/combined_pair_summaries.tsv.gz"
        ),
        "spatial_summaries": (
            "focused_selection_EAS_sim/results/aggregated_spatial_profiles.tsv.gz"
        ),
    }
    inputs = {"empirical_path": None, "empirical_sha256": None}
    for key, relative in relative_paths.items():
        path = repo / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f"{key}\n", encoding="utf-8")
        inputs[f"{key}_path"] = relative
        inputs[f"{key}_sha256"] = clear._sha256_file(path)
    contract = {"inputs": inputs}
    spec = clear.AnalysisSourceSpec("gamma_smc", "Gamma-SMC", "gamma_smc")

    checks, paths, records = clear._analysis_input_checks(
        contract, repo_root=repo, spec=spec
    )
    assert set(paths) == set(relative_paths)
    assert records["empirical"] is None
    assert set(clear._audit_file_checks(checks, workers=2)["status"]) == {"verified"}

    broken = copy.deepcopy(contract)
    broken["inputs"]["class_summaries_path"] = relative_paths["pair_summaries"]
    with pytest.raises(ValueError, match="input path differs"):
        clear._analysis_input_checks(broken, repo_root=repo, spec=spec)
    broken = copy.deepcopy(contract)
    broken["inputs"]["empirical_path"] = relative_paths["class_summaries"]
    with pytest.raises(ValueError, match="null analysis input differs"):
        clear._analysis_input_checks(broken, repo_root=repo, spec=spec)
    broken = copy.deepcopy(contract)
    broken["inputs"]["spatial_summaries_sha256"] = "b" * 64
    checks, _, _ = clear._analysis_input_checks(broken, repo_root=repo, spec=spec)
    with pytest.raises(ValueError, match="input hash differs"):
        clear._audit_file_checks(checks, workers=2)


def test_focused_aggregate_lineage_rejects_row_tamper():
    direct = pd.DataFrame(
        {
            "unit_id": ["u", "u"],
            "genotype_class": ["overall", "overall"],
            "threshold_years": [1_000, 4_500],
            "focal_mean_p_tmrca_lt_threshold": [0.01, 0.02],
            "focal_mean_tmrca_generations": [10_000.0, 10_000.0],
        }
    )
    clear._assert_focused_class_lineage(direct.copy(), direct, source="gamma_smc")
    truth_roundtrip = direct.copy()
    truth_roundtrip["focal_mean_p_tmrca_lt_threshold"] += 4e-13
    clear._assert_focused_class_lineage(truth_roundtrip, direct, source="tree_truth")
    truth_roundtrip.loc[0, "focal_mean_p_tmrca_lt_threshold"] += 1e-6
    with pytest.raises(ValueError, match="differs from unit receipts"):
        clear._assert_focused_class_lineage(
            truth_roundtrip, direct, source="tree_truth"
        )


def test_direct_sensitivity_pointwise_and_minp_use_matched_score_vectors():
    scores = _combined_scores()
    pointwise = clear._sensitivity_pointwise_from_scores(scores)
    minp = clear._sensitivity_minp_from_scores(scores, pointwise)

    assert len(pointwise) == 140
    assert len(minp) == 20
    row = pointwise[
        (pointwise["source"] == "gamma_smc")
        & (pointwise["replicate_index"] == 0)
        & (pointwise["threshold_years"] == 10_000)
    ].iloc[0]
    neutral = scores[
        (scores["source"] == "gamma_smc")
        & (scores["bank_role"] == "matched_af_sensitivity")
        & (scores["threshold_years"] == 10_000)
    ]["p_hom_alt_minus_hom_ref"].to_numpy(dtype=float)
    observed = float(row["observed_selection_score"])
    exceedances = int(np.count_nonzero(neutral >= observed))
    assert int(row["null_upper_exceedance_count"]) == exceedances
    assert float(row["mc_p_upper"]) == pytest.approx((exceedances + 1) / 101)
    family = pointwise[
        (pointwise["source"] == "gamma_smc") & (pointwise["unit_id"] == row["unit_id"])
    ].sort_values("threshold_years")
    assert np.allclose(
        family["bh_q_upper"], clear._bh_fdr(family["mc_p_upper"].to_numpy())
    )

    neutral_matrix = (
        scores[
            (scores["source"] == "gamma_smc")
            & (scores["bank_role"] == "matched_af_sensitivity")
        ]
        .pivot(
            index="unit_id",
            columns="threshold_years",
            values="p_hom_alt_minus_hom_ref",
        )
        .reindex(columns=clear.THRESHOLDS_YEARS)
        .sort_index()
        .to_numpy()
    )
    selected_vector = (
        scores[
            (scores["source"] == "gamma_smc")
            & (scores["bank_role"] == "selected")
            & (scores["unit_id"] == row["unit_id"])
        ]
        .set_index("threshold_years")
        .reindex(clear.THRESHOLDS_YEARS)["p_hom_alt_minus_hom_ref"]
        .to_numpy()
    )
    _, observed_min, minp_exceedances, _, omnibus = clear._inclusive_pooled_minp(
        neutral_matrix, selected_vector
    )
    omnibus_row = minp[
        (minp["source"] == "gamma_smc") & (minp["unit_id"] == row["unit_id"])
    ].iloc[0]
    assert omnibus_row["observed_min_marginal_p_upper"] == pytest.approx(observed_min)
    assert int(omnibus_row["neutral_minp_exceedance_count"]) == minp_exceedances
    assert float(omnibus_row["omnibus_p_upper"]) == pytest.approx(omnibus)
    legacy_pointwise = pointwise[pointwise["source"] == "gamma_smc"].drop(
        columns=["source", "estimand", "inference_scope"]
    )
    legacy_pointwise = legacy_pointwise.assign(
        demography_id=clear.DEMOGRAPHY_ID,
        target_allele_frequency=clear.TARGET_AF,
        selection_coefficient=clear.SELECTION_COEFFICIENT,
        statistic_scope="focal_nearest",
        metric="p_tmrca_lt_threshold",
    )
    legacy_minp = minp[minp["source"] == "gamma_smc"].drop(
        columns=["source", "estimand", "inference_scope"]
    )
    legacy_minp = legacy_minp.assign(
        demography_id=clear.DEMOGRAPHY_ID,
        target_allele_frequency=clear.TARGET_AF,
        selection_coefficient=clear.SELECTION_COEFFICIENT,
        statistic_scope="focal_nearest",
        metric="p_tmrca_lt_threshold",
    )
    bundle = clear.AnalysisBundle(
        spec=clear.AnalysisSourceSpec("gamma_smc", "Gamma-SMC", "gamma_smc"),
        completion={},
        completion_path=Path("completion.json"),
        tables={
            "selected_pvalues": legacy_pointwise,
            "primary_omnibus": legacy_minp,
        },
        input_paths={},
        input_contract={},
        audit=pd.DataFrame(),
    )
    published_pointwise = clear._sensitivity_pointwise(bundle).sort_values(
        ["unit_id", "threshold_years"]
    )
    derived_pointwise = pointwise[pointwise["source"] == "gamma_smc"].sort_values(
        ["unit_id", "threshold_years"]
    )
    pd.testing.assert_frame_equal(
        published_pointwise.reset_index(drop=True),
        derived_pointwise.reset_index(drop=True),
        check_dtype=False,
    )
    published_minp = clear._sensitivity_minp(bundle).sort_values("unit_id")
    derived_minp = minp[minp["source"] == "gamma_smc"].sort_values("unit_id")
    pd.testing.assert_frame_equal(
        published_minp.reset_index(drop=True),
        derived_minp.reset_index(drop=True),
        check_dtype=False,
    )


def test_real_focused_cache_receipt_if_available():
    repo = Path(__file__).resolve().parents[1]
    campaign = repo / clear.DEFAULT_CAMPAIGN_ROOT
    execution_path = campaign / "execution_units.tsv"
    unit_id = "ancient_eurasia_han_introgression__af20__neutral__rep001"
    completion_path = campaign / "work" / unit_id / "simulation_complete.json"
    if not execution_path.is_file() or not completion_path.is_file():
        pytest.skip("frozen focused-selection cache is not present")
    execution = clear._select_existing_units(
        pd.read_csv(execution_path, sep="\t", low_memory=False)
    )
    unit = execution.loc[execution["unit_id"] == unit_id].iloc[0].to_dict()
    checks, binding, completion = clear._unit_receipt_checks(repo, campaign, unit)
    audit = clear._audit_file_checks(checks, workers=4)

    assert set(audit["status"]) == {"verified"}
    assert completion["schema"] == clear.FOCUSED_SIMULATION_SCHEMA_VERSION
    assert "individual_location_canonicalization" not in completion
    assert binding["decode_contract_sha256"]


def test_clear_loader_defaults_only_to_frozen_survival_v2():
    assert clear.DEFAULT_SURVIVAL_RESULTS_DIR.as_posix().endswith(
        "ancient_eurasia_single_origin_survival_v2/results"
    )
    assert {
        clear.SURVIVAL_PLAN_SCHEMA_VERSION,
        clear.SURVIVAL_SIMULATION_SCHEMA_VERSION,
        clear.SURVIVAL_DECODE_SCHEMA_VERSION,
        clear.SURVIVAL_RESULTS_SCHEMA_VERSION,
        clear.SURVIVAL_DECODER_RUN_SCHEMA_VERSION,
    } == {
        "gamma-smc.ancient-eurasia-survival-neutral-plan/v2",
        "gamma-smc.ancient-eurasia-survival-neutral-simulation/v2",
        "gamma-smc.ancient-eurasia-survival-neutral-decode/v2",
        "gamma-smc.ancient-eurasia-survival-neutral-results/v2",
        "gamma-smc.ancient-eurasia-survival-neutral-decoder-run/v2",
    }


def test_survival_estimand_is_exact_and_census_conditioned():
    clear._validate_survival_estimand({"estimand": _survival_estimand()})

    broken = _survival_estimand()
    broken["present_conditioning"] = "strict_segregating_in_500_diploid_pool"
    with pytest.raises(ValueError, match="present_conditioning"):
        clear._validate_survival_estimand({"estimand": broken})


def test_survival_event_schedule_is_exact_six_event_pulse_route():
    schedule = _survival_event_schedule()
    clear._validate_survival_event_schedule(schedule, clear._canonical_sha256(schedule))

    lost_route = json.loads(json.dumps(schedule))
    lost_route[1]["allele_frequency"] = 1e-9
    with pytest.raises(ValueError, match="six-event route"):
        clear._validate_survival_event_schedule(
            lost_route, clear._canonical_sha256(lost_route)
        )


def test_survival_plan_contract_binds_seed_sources_and_event_transform():
    schedule = _survival_event_schedule()
    transform = {
        "source_constructor": (
            "focused_selection_simulation.build_raw_han_selected_sweep"
        ),
        "source_constructor_slim_scaling_factor": 5.0,
        "removed_change_mutation_fitness_count": 2,
        "removed_terminal_han_af_band_count": 2,
        "added_terminal_han_strict_upper_count": 1,
        "selection_coefficient_after_transform": 0.0,
        "nominal_pulse_generations_ago": 2_272.0,
        "realized_pulse_generations_ago_q5": 2_270.0,
        "nominal_han_split_generations_ago": 2_016.0,
        "realized_han_split_generations_ago_q5": 2_015.0,
        "event_schedule": schedule,
    }
    contract = {
        "study_id": "ancient_eurasia_single_origin_survival",
        "expected_units": 100,
        "estimand": _survival_estimand(),
        "simulation_parameters": clear.SURVIVAL_SIMULATION_PARAMETERS,
        "decoder_settings": clear.SURVIVAL_DECODER_SETTINGS,
        "decoder_binary": {
            "identity": clear.EXPECTED_GAMMA_BINARY_IDENTITY,
            "sha256": clear.EXPECTED_GAMMA_BINARY_SHA256,
            "size_bytes": 1,
        },
        "execution_units_sha256": "a" * 64,
        "source_hashes": {path: "b" * 64 for path in clear.SURVIVAL_SOURCE_PATHS},
        "event_transform": transform,
        "runtime_versions": {
            key: "1"
            for key in (
                "python",
                "numpy",
                "pandas",
                "scipy",
                "stdpopsim",
                "pyslim",
                "msprime",
                "tskit",
            )
        },
        "base_seed": clear.SURVIVAL_BASE_SEED,
        "seed_derivation": clear.SURVIVAL_SEED_DERIVATION,
        "simulation_law": clear.SURVIVAL_SIMULATION_LAW,
        "focal_overlay_patch": clear.SURVIVAL_FOCAL_OVERLAY_PATCH,
        "individual_location_canonicalization": (
            clear.SURVIVAL_INDIVIDUAL_LOCATION_CANONICALIZATION
        ),
    }
    clear._validate_survival_plan_contract(contract)
    broken = json.loads(json.dumps(contract))
    broken["event_transform"]["realized_pulse_generations_ago_q5"] = 2_272
    with pytest.raises(ValueError, match="event-transform"):
        clear._validate_survival_plan_contract(broken)

    broken = copy.deepcopy(contract)
    broken["focal_overlay_patch"]["expected_patched_call_count"] = 2
    with pytest.raises(ValueError, match="focal-overlay contract"):
        clear._validate_survival_plan_contract(broken)

    type_tamper = copy.deepcopy(contract)
    type_tamper["focal_overlay_patch"]["expected_patched_call_count"] = True
    with pytest.raises(ValueError, match="focal-overlay contract"):
        clear._validate_survival_plan_contract(type_tamper)

    broken = copy.deepcopy(contract)
    broken["individual_location_canonicalization"]["expected_individual_count"] = 499
    with pytest.raises(ValueError, match="individual-location contract"):
        clear._validate_survival_plan_contract(broken)


def test_survival_endpoint_schema_is_exact_and_census_conditioned():
    endpoints = clear._validate_survival_endpoints(_survival_endpoints())
    clear._verify_survival_score_endpoints(_survival_scores(), endpoints)
    assert endpoints["final_census_segregating"].all()
    assert endpoints["pool_alt_count"].min() == 0
    assert endpoints["pool_alt_count"].max() == 1_000

    missing = endpoints[["unit_id"]]
    with pytest.raises(ValueError, match="columns"):
        clear._validate_survival_endpoints(missing)

    failed_census = endpoints.copy()
    failed_census.loc[0, ["final_census_alt_count", "final_census_af"]] = 0
    with pytest.raises(ValueError, match="strictly segregating"):
        clear._validate_survival_endpoints(failed_census)

    failed_seed = _survival_execution()
    failed_seed.loc[0, "panel_seed"] += 1
    with pytest.raises(ValueError, match="deterministic seed schedule"):
        clear._validate_survival_execution(failed_seed)


def test_v2_overlay_receipt_is_exact_and_tamper_evident():
    receipt = _focal_overlay_receipt()
    clear._validate_survival_focal_overlay_receipt(receipt)

    wrong_schema = copy.deepcopy(receipt)
    wrong_schema["schema"] = "gamma-smc.ancient-eurasia-reserved-focal-overlay-patch/v0"
    with pytest.raises(ValueError, match="receipt fields"):
        clear._validate_survival_focal_overlay_receipt(wrong_schema)

    escaped = copy.deepcopy(receipt)
    escaped["calls"][1]["passed_rate_map"]["rates"][1] = 1.25e-8
    with pytest.raises(ValueError, match="call content"):
        clear._validate_survival_focal_overlay_receipt(escaped)

    unbound = copy.deepcopy(receipt)
    unbound["calls"][0]["model_next_id"] = -1
    with pytest.raises(ValueError, match="call content"):
        clear._validate_survival_focal_overlay_receipt(unbound)

    type_tamper = copy.deepcopy(receipt)
    type_tamper["intercepted_call_count"] = True
    with pytest.raises(ValueError, match="call inventory"):
        clear._validate_survival_focal_overlay_receipt(type_tamper)


def test_v2_frozen_generator_source_hash_is_required():
    repo = Path(__file__).resolve().parents[1]
    source_hashes = {
        identity: clear._canonical_text_sha256(repo / identity)
        for identity in clear.SURVIVAL_SOURCE_PATHS
    }
    assert (
        source_hashes["python/gamma_smc_aou/ancient_eurasia_survival_neutral.py"]
        == clear.EXPECTED_SURVIVAL_GENERATOR_SOURCE_SHA256
    )
    clear._survival_source_audit(source_hashes, repo_root=repo)

    source_hashes["python/gamma_smc_aou/ancient_eurasia_survival_neutral.py"] = "0" * 64
    with pytest.raises(ValueError, match="frozen generator source hash"):
        clear._survival_source_audit(source_hashes, repo_root=repo)


def test_survival_absent_focal_identity_remains_valid_only_for_zero_pool_count():
    tables = tskit.TableCollection(sequence_length=10_000_000)
    tables.metadata_schema = tskit.MetadataSchema.permissive_json()
    tables.metadata = {"SLiM": {"cycle": 500, "user_metadata": {"Q": [5.0]}}}
    tree_sequence = tables.tree_sequence()

    assert clear._validate_survival_tree_focal_identity(
        tree_sequence,
        _survival_focal_expectation(),
        _absent_focal_identity(),
        pool_alt_count=0,
    ) == (None, None)
    with pytest.raises(ValueError, match="absent despite pool detection"):
        clear._validate_survival_tree_focal_identity(
            tree_sequence,
            _survival_focal_expectation(),
            _absent_focal_identity(),
            pool_alt_count=1,
        )


def test_survival_panel_tree_manifest_and_explicit_pairs_are_bound(tmp_path):
    unit_dir = tmp_path / "unit"
    decode_dir = unit_dir / "decoded"
    decode_dir.mkdir(parents=True)
    endpoint_row = _survival_endpoints().iloc[0]
    endpoint = {
        key: endpoint_row[key].item()
        if isinstance(endpoint_row[key], np.generic)
        else endpoint_row[key]
        for key in clear.SURVIVAL_ENDPOINT_COLUMNS[4:]
    }
    chosen = sorted(
        np.random.default_rng(int(endpoint["panel_seed"]))
        .choice(500, size=100, replace=False)
        .tolist()
    )
    chosen_set = set(chosen)
    carrier_row = next(index for index in range(500) if index not in chosen_set)
    endpoint["chosen_panel_rows_zero_based"] = chosen
    endpoint.update(
        {
            "pool_alt_count": 1,
            "pool_af": 0.001,
            "pool_detected": True,
            "pool_fixed": False,
        }
    )

    tables = tskit.TableCollection(sequence_length=10_000_000)
    tables.metadata_schema = tskit.MetadataSchema.permissive_json()
    tables.metadata = {
        "SLiM": {
            "cycle": 500,
            "user_metadata": {
                "Q": [5.0],
                "survival_neutral_final_census_alt_count": [
                    int(endpoint["final_census_alt_count"])
                ],
                "survival_neutral_final_census_total_count": [
                    int(endpoint["final_census_total_count"])
                ],
                "survival_neutral_final_census_af": [
                    float(endpoint["final_census_af"])
                ],
                "survival_neutral_final_census_segregating": [True],
                "survival_neutral_final_census_fixed": [False],
            },
        }
    }
    tables.populations.metadata_schema = tskit.MetadataSchema.permissive_json()
    tables.mutations.metadata_schema = tskit.MetadataSchema.permissive_json()
    for population_id in range(8):
        name = (
            "Han"
            if population_id == 5
            else "Neanderthal"
            if population_id == 7
            else f"population_{population_id}"
        )
        tables.populations.add_row(metadata={"id": name, "name": name})
    han = 5
    for index in range(500):
        individual = tables.individuals.add_row(location=[0.0, 0.0, 0.0])
        flags = tskit.NODE_IS_SAMPLE if index in chosen_set else 0
        tables.nodes.add_row(flags=flags, time=0, population=han, individual=individual)
        tables.nodes.add_row(flags=flags, time=0, population=han, individual=individual)
    root = tables.nodes.add_row(time=3_000, population=7)
    for child in range(1_000):
        tables.edges.add_row(0, tables.sequence_length, parent=root, child=child)
    site = tables.sites.add_row(position=clear.FOCAL_POSITION_BP, ancestral_state="A")
    mutation = tables.mutations.add_row(
        site=site,
        node=2 * carrier_row,
        time=2_400,
        derived_state="C",
        metadata={
            "mutation_list": [
                {
                    "mutation_type": 1,
                    "selection_coeff": 0.0,
                    "subpopulation": 7,
                    "slim_time": 20,
                    "nucleotide": 1,
                }
            ]
        },
    )
    tables.sort()
    tree_sequence = tables.tree_sequence()
    assert tree_sequence.node(tree_sequence.mutation(mutation).node).population == han
    tree_sequence.dump(unit_dir / "simulation.trees")
    location_receipt = _individual_location_receipt(tree_sequence)

    manifest = pd.DataFrame(
        {
            "vcf_diploid_index": np.arange(100),
            "tree_sequence_individual_id": chosen,
            "sample_node_0": np.asarray(chosen) * 2,
            "sample_node_1": np.asarray(chosen) * 2 + 1,
            "gamma_smc_haplotype_0": np.arange(100) * 2,
            "gamma_smc_haplotype_1": np.arange(100) * 2 + 1,
            "focal_position_0based": clear.FOCAL_POSITION_BP,
            "focal_selected_allele_index": 1,
            "focal_selected_allele_count": 0,
            "genotype_class": "hom_ref",
        },
        columns=clear.SURVIVAL_SAMPLE_MANIFEST_COLUMNS,
    )
    manifest.to_csv(unit_dir / "sample_manifest.tsv", sep="\t", index=False)
    pair_lines = [
        "# Gamma-SMC explicit haplotype pairs",
        "# genotype_class\toverall",
        "# n_pairs\t100",
        *(f"{2 * index}\t{2 * index + 1}" for index in range(100)),
    ]
    pair_path = decode_dir / "overall.pairs.tsv"
    pair_path.write_text("\n".join(pair_lines) + "\n", encoding="utf-8")

    clear._validate_survival_panel_artifacts(
        unit_dir,
        decode_dir,
        endpoint,
        _survival_focal_expectation(),
        _valid_focal_identity(),
        location_receipt,
    )
    tampered_endpoint = dict(endpoint)
    tampered_endpoint.update(
        {"sample_alt_count": 1, "sample_af": 0.005, "sample_detected": True}
    )
    manifest.loc[0, ["focal_selected_allele_count", "genotype_class"]] = [
        1,
        "heterozygous",
    ]
    manifest.to_csv(unit_dir / "sample_manifest.tsv", sep="\t", index=False)
    with pytest.raises(ValueError, match="tree panel genotypes"):
        clear._validate_survival_panel_artifacts(
            unit_dir,
            decode_dir,
            tampered_endpoint,
            _survival_focal_expectation(),
            _valid_focal_identity(),
            location_receipt,
        )
    manifest.loc[0, ["focal_selected_allele_count", "genotype_class"]] = [0, "hom_ref"]
    manifest.to_csv(unit_dir / "sample_manifest.tsv", sep="\t", index=False)
    pair_lines[-1] = "0\t199"
    pair_path.write_text("\n".join(pair_lines) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="explicit pair file"):
        clear._validate_survival_panel_artifacts(
            unit_dir,
            decode_dir,
            endpoint,
            _survival_focal_expectation(),
            _valid_focal_identity(),
            location_receipt,
        )

    arbitrary_before = copy.deepcopy(location_receipt)
    arbitrary_before["before"].update(
        {
            "location_all_finite": False,
            "location_finite_count": 1_497,
            "location_nan_count": 1,
            "location_positive_infinity_count": 1,
            "location_negative_infinity_count": 1,
            "location_all_zero": False,
            "location_nonzero_count": 3,
            "location_record": {
                "dtype": "<f8",
                "shape": [1_500],
                "sha256": "a" * 64,
            },
        }
    )
    clear._validate_survival_individual_location_receipt(
        arbitrary_before, tree=tree_sequence
    )

    paired_tamper = copy.deepcopy(location_receipt)
    for stage in ("before", "after", "stored_tree"):
        paired_tamper[stage]["all_except_individual_location_sha256"] = "b" * 64
    with pytest.raises(ValueError, match="stored tree"):
        clear._validate_survival_individual_location_receipt(
            paired_tamper, tree=tree_sequence
        )


def test_survival_producer_python_command_is_cross_host_portable():
    observed = [
        "/mnt/c/Users/Lenovo/codex/gamma_smc_ts/.venv-wsl/bin/python",
        "-c",
        "producer",
    ]
    expected = [
        r"C:\Users\Lenovo\codex\gamma_smc_ts\.venv-wsl\bin\python",
        "-c",
        "producer",
    ]
    clear._require_exact_portable_command(
        observed,
        expected,
        path_indices={0},
        label="cross-host producer command",
    )


def test_survival_decode_raw_inputs_and_explicit_command_are_bound(tmp_path):
    repo = tmp_path / "repo"
    unit_dir = repo / "study/work/unit"
    decode_dir = unit_dir / "decoded"
    binary = repo / clear.EXPECTED_GAMMA_BINARY_IDENTITY
    producer_python = repo / clear.EXPECTED_SURVIVAL_PRODUCER_PYTHON_IDENTITY
    decode_dir.mkdir(parents=True)
    binary.parent.mkdir(parents=True)
    binary.write_bytes(b"gamma")
    producer_python.parent.mkdir(parents=True)
    producer_python.symlink_to("/usr/bin/python3")
    assert producer_python.resolve() != producer_python
    tree = unit_dir / "simulation.trees"
    manifest = unit_dir / "sample_manifest.tsv"
    pairs = decode_dir / "overall.pairs.tsv"
    tree.write_bytes(b"tree")
    manifest.write_text("x\n" + "1\n" * 100, encoding="utf-8")
    pairs.write_text("0\t1\n" * 100, encoding="utf-8")
    raw_inputs = {
        "simulation_tree": _output_record(
            tree, relative="../simulation.trees", rows=None
        ),
        "sample_manifest": _output_record(
            manifest, relative="../sample_manifest.tsv", rows=100
        ),
        "overall_pairs": _output_record(pairs, relative="overall.pairs.tsv", rows=100),
    }
    checks = clear._validate_survival_decode_raw_inputs(
        raw_inputs,
        unit_id="unit",
        unit_dir=unit_dir,
        decode_dir=decode_dir,
        repo_root=repo,
    )
    assert set(clear._audit_file_checks(checks, workers=1)["status"]) == {"verified"}

    command = [
        str(binary),
        "--input",
        "/dev/stdin",
        "--input_format",
        "vcf",
        "--scaled_mutation_rate",
        "0.000315",
        "--recombination_to_mutation_ratio",
        "0.8",
        "--unscaled_mutation_rate",
        "1.25e-08",
        "--recent_threshold_years",
        "1000,4500,10000,20000,30000,40000,50000",
        "--generation_time",
        "25.0",
        "--recent_summary",
        str(decode_dir / "overall.summary.tsv"),
        "--recent_call",
        "median",
        "--recent_call_probability",
        "0.5",
        "--output_at_hets=false",
        "--output_at_stride",
        "10000",
        "--cache_size",
        "1000",
        "--threads",
        "1",
        "--pair_block",
        "256",
        "--backward_alignment",
        "fixed",
        "--exp10",
        "accurate",
        "--pairs_file",
        str(pairs),
        "--output",
        str(decode_dir / "posterior.zst"),
    ]
    run_receipt = {
        "schema": clear.SURVIVAL_DECODER_RUN_SCHEMA_VERSION,
        "unit_id": "unit",
        "raw_inputs": raw_inputs,
        "run": {
            "command": command,
            "input_format": "trees",
            "input_path": str(tree),
            "stdout": "",
            "stderr": "",
            "decode_seconds": 1.25,
            "stride_bp": 10_000,
            "cache_size_bp": 1_000,
            "n_output_positions": 1_000,
            "pairs_manifest": None,
            "n_pairs_recorded": None,
            "tree_sequence_vcf_producer_command": [
                str(producer_python),
                "-c",
                (
                    "import sys; from gamma_smc_aou.tree_sequence import "
                    "stream_tree_sequence_vcf; stream_tree_sequence_vcf(sys.argv[1], "
                    "sys.stdout, input_format=sys.argv[2])"
                ),
                str(tree),
                "trees",
            ],
        },
    }
    run_path = decode_dir / "decoder_run.json"
    run_path.write_text(json.dumps(run_receipt) + "\n", encoding="utf-8")
    clear._validate_survival_decoder_run(
        run_path,
        unit_id="unit",
        raw_inputs=raw_inputs,
        unit_dir=unit_dir,
        decode_dir=decode_dir,
        repo_root=repo,
    )
    run_receipt["run"]["command"].append("--extra")
    run_path.write_text(json.dumps(run_receipt) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="token count differs"):
        clear._validate_survival_decoder_run(
            run_path,
            unit_id="unit",
            raw_inputs=raw_inputs,
            unit_dir=unit_dir,
            decode_dir=decode_dir,
            repo_root=repo,
        )


def test_survival_focal_reducers_match_truth_and_gamma_source_rows(tmp_path):
    unit_dir = tmp_path / "unit"
    decode_dir = unit_dir / "decoded"
    decode_dir.mkdir(parents=True)
    unit_scores = _survival_scores()
    unit_scores = unit_scores[unit_scores["replicate_index"] == 0].copy()
    truth_scores = unit_scores[unit_scores["source"] == "tree_truth"].sort_values(
        "threshold_years"
    )
    truth = pd.DataFrame(
        {
            "source": "tree_truth",
            "position_0based": clear.FOCAL_POSITION_BP,
            "genotype_class": "overall",
            "n_pairs": 100,
            "threshold_years": truth_scores["threshold_years"].to_numpy(),
            "p_tmrca_lt_threshold": truth_scores[
                "overall_focal_p_tmrca_lt_threshold"
            ].to_numpy(),
        }
    )
    truth.to_csv(unit_dir / "truth_profiles.tsv.gz", sep="\t", index=False)
    gamma_scores = unit_scores[unit_scores["source"] == "gamma_smc"].sort_values(
        "threshold_years"
    )
    gamma_row = {"position_0based": clear.FOCAL_POSITION_BP, "n_pairs": 100}
    for threshold, value in zip(
        gamma_scores["threshold_years"],
        gamma_scores["overall_focal_p_tmrca_lt_threshold"],
        strict=True,
    ):
        gamma_row[f"mean_p_lt_{int(threshold)}"] = value
    pd.DataFrame([gamma_row]).to_csv(
        decode_dir / "overall.summary.tsv", sep="\t", index=False
    )

    clear._verify_survival_focal_reducers(unit_dir, decode_dir, unit_scores)
    gamma_row["mean_p_lt_50000"] += 0.1
    pd.DataFrame([gamma_row]).to_csv(
        decode_dir / "overall.summary.tsv", sep="\t", index=False
    )
    with pytest.raises(ValueError, match="focal decoder row"):
        clear._verify_survival_focal_reducers(unit_dir, decode_dir, unit_scores)


def test_survival_root_schema_endpoint_cartesian_and_per_unit_score_binding(
    tmp_path, monkeypatch
):
    repo = tmp_path / "repo"
    study = repo / "focused_selection_EAS_sim/ancient_eurasia_single_origin_survival_v2"
    plan_dir = study / "plan"
    results = study / "results"
    plan_dir.mkdir(parents=True)
    results.mkdir()
    execution = _survival_execution()
    execution.to_csv(plan_dir / "execution_units.tsv", sep="\t", index=False)
    (plan_dir / "study_plan.json").write_text("{}\n", encoding="utf-8")
    scores = _survival_scores()
    endpoints = _survival_endpoints()
    scores.to_csv(results / "replicate_scores.tsv", sep="\t", index=False)
    endpoints.to_csv(results / "endpoint_summary.tsv", sep="\t", index=False)
    simulation_status = execution.copy()
    simulation_status["status"] = "complete"
    simulation_status["error"] = ""
    simulation_status["final_census_af"] = endpoints["final_census_af"]
    simulation_status["simulation_completion_sha256"] = "a" * 64
    decode_status = execution.copy()
    decode_status["status"] = "complete"
    decode_status["error"] = ""
    decode_status["decode_completion_sha256"] = "b" * 64
    simulation_status.to_csv(results / "simulation_status.tsv", sep="\t", index=False)
    decode_status.to_csv(results / "decode_status.tsv", sep="\t", index=False)

    decoder = {
        "identity": clear.EXPECTED_GAMMA_BINARY_IDENTITY,
        "sha256": clear.EXPECTED_GAMMA_BINARY_SHA256,
        "size_bytes": 1,
    }
    plan_contract = {
        "estimand": _survival_estimand(),
        "source_hashes": {},
        "decoder_binary": decoder,
    }
    plan_sha = "c" * 64
    execution_sha = clear._sha256_file(plan_dir / "execution_units.tsv")
    root_contract = {
        "study_id": "ancient_eurasia_single_origin_survival",
        "plan_sha256": plan_sha,
        "execution_units_sha256": execution_sha,
        "estimand": _survival_estimand(),
        "expected_units": 100,
        "schemas": {
            "plan": clear.SURVIVAL_PLAN_SCHEMA_VERSION,
            "simulation": clear.SURVIVAL_SIMULATION_SCHEMA_VERSION,
            "decode": clear.SURVIVAL_DECODE_SCHEMA_VERSION,
            "results": clear.SURVIVAL_RESULTS_SCHEMA_VERSION,
        },
        "decoder_settings": clear.SURVIVAL_DECODER_SETTINGS,
        "decoder_binary": decoder,
        "source_hashes": {},
        "simulation_law": clear.SURVIVAL_SIMULATION_LAW,
    }
    outputs = {
        key: _output_record(
            results / relative,
            relative=relative,
            rows=100 if key != "replicate_scores" else 1_400,
        )
        for key, relative in clear.SURVIVAL_RESULTS_OUTPUT_PATHS.items()
    }
    inventory = [
        {
            "unit_id": unit_id,
            "replicate_index": index,
            "seed": clear._survival_seed(unit_id),
            "simulation_completion_path": f"work/{unit_id}/simulation_complete.json",
            "simulation_completion_sha256": "a" * 64,
            "decode_completion_path": f"work/{unit_id}/decoded/decode_complete.json",
            "decode_completion_sha256": "b" * 64,
        }
        for index, unit_id in enumerate(clear._expected_survival_unit_ids())
    ]
    completion = {
        "schema": clear.SURVIVAL_RESULTS_SCHEMA_VERSION,
        "status": "complete",
        "created_utc": "2026-08-17T00:00:00Z",
        "contract": root_contract,
        "contract_sha256": clear._canonical_sha256(root_contract),
        "outputs": outputs,
        "unit_inventory": inventory,
    }
    completion_path = results / "RESULTS_COMPLETION.json"
    completion_path.write_text(
        json.dumps(completion, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    monkeypatch.setattr(
        clear,
        "_load_survival_plan",
        lambda *_args, **_kwargs: (
            plan_contract,
            plan_sha,
            clear._validate_survival_execution(execution),
            [],
            pd.DataFrame(),
        ),
    )

    def fake_unit(_repo, _study, unit, _inventory, **_kwargs):
        subset = scores[scores["unit_id"] == unit["unit_id"]].copy()
        endpoint_row = endpoints[endpoints["unit_id"] == unit["unit_id"]].iloc[0]
        endpoint = {
            key: (
                endpoint_row[key].item()
                if isinstance(endpoint_row[key], np.generic)
                else endpoint_row[key]
            )
            for key in clear.SURVIVAL_ENDPOINT_COLUMNS[4:]
        }
        endpoint["chosen_panel_rows_zero_based"] = sorted(
            np.random.default_rng(int(endpoint["panel_seed"]))
            .choice(500, size=100, replace=False)
            .tolist()
        )
        return (
            [],
            {
                "unit_id": unit["unit_id"],
                "simulation_completion_sha256": "a" * 64,
                "decode_completion_sha256": "b" * 64,
                "endpoint": endpoint,
            },
            subset,
        )

    monkeypatch.setattr(clear, "_survival_unit_receipts", fake_unit)
    bundle = clear._load_survival_bundle(repo, results, workers=4)
    assert len(bundle.scores) == 1_400
    assert len(bundle.endpoints) == 100
    assert bundle.input_contract["unit_receipt_inventory_sha256"]

    completion["schema"] = "gamma-smc.ancient-eurasia-survival-neutral-results/v0"
    completion_path.write_text(json.dumps(completion) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="schema"):
        clear._load_survival_bundle(repo, results, workers=4)
    completion["schema"] = clear.SURVIVAL_RESULTS_SCHEMA_VERSION

    endpoints[["unit_id"]].to_csv(
        results / "endpoint_summary.tsv", sep="\t", index=False
    )
    completion["outputs"]["endpoint_summary"] = _output_record(
        results / "endpoint_summary.tsv", relative="endpoint_summary.tsv", rows=100
    )
    completion_path.write_text(json.dumps(completion) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="endpoint columns"):
        clear._load_survival_bundle(repo, results, workers=4)


def test_plus_one_pointwise_minp_and_power_math():
    scores = _combined_scores()
    pointwise = clear._primary_pointwise(scores)
    minp = clear._primary_minp(scores, pointwise)
    power = clear._power_from_pointwise(pointwise)

    assert len(pointwise) == 2 * 10 * 7
    assert np.allclose(
        pointwise["mc_p_upper"],
        (pointwise["null_upper_exceedance_count"] + 1) / (pointwise["n_null"] + 1),
    )
    assert len(minp) == 20
    assert np.allclose(
        minp["omnibus_p_upper"],
        (minp["neutral_minp_exceedance_count"] + 1) / (minp["n_null"] + 1),
    )
    assert len(power) == 14
    assert set(power["n_selected"]) == {10}
    assert set(pointwise["inference_scope"]) == {
        "simulation_reference_tail_combines_selection_and_endpoint_af_enrichment"
    }
    assert "exchangeability" not in " ".join(minp["omnibus_method"])
    assert set(minp["omnibus_method"]) == {
        "inclusive_pooled_101_marginal_ranks_minP_plus_one_reference_tail"
    }


def test_inclusive_pooled_minp_ties_are_hand_calculable():
    neutral = np.array([[1.0, 3.0], [2.0, 2.0], [3.0, 1.0]])
    observed = np.array([2.0, 2.0])
    marginal, observed_min, exceedances, best, reference_p = (
        clear._inclusive_pooled_minp(neutral, observed)
    )

    assert marginal == pytest.approx([0.75, 0.75])
    assert observed_min == pytest.approx(0.75)
    assert exceedances == 3
    assert best == 0
    assert reference_p == pytest.approx(1.0)


def test_bh_is_monotone_after_unordering_and_rejects_invalid_input():
    pvalues = np.array([0.04, 0.001, 0.03, 0.7])
    adjusted = clear._bh_fdr(pvalues)
    order = np.argsort(pvalues)
    assert np.all(np.diff(adjusted[order]) >= -1e-15)
    assert np.all(adjusted >= pvalues)
    with pytest.raises(ValueError, match="finite probabilities"):
        clear._bh_fdr([0.2, np.nan])


def test_plot_layouts_are_letter_sized_and_do_not_author_images(tmp_path, monkeypatch):
    scores = _combined_scores()
    pointwise = clear._primary_pointwise(scores)
    power = clear._power_from_pointwise(pointwise)
    comparison = _comparison(pointwise, power)
    captured = []

    def capture(figure, directory, stem):
        captured.append((figure, stem))
        return [directory / f"{stem}.png", directory / f"{stem}.pdf"]

    monkeypatch.setattr(clear, "_save_pair", capture)
    truth_paths = clear.plot_pulse_aligned_null(
        scores, _endpoints(), source="tree_truth", output_dir=tmp_path
    )
    threshold_paths = clear.plot_threshold_results(
        scores, pointwise, power, output_dir=tmp_path
    )
    timing_paths = clear.plot_timing_endpoints_diagnostic(
        _endpoints(), comparison, output_dir=tmp_path
    )

    assert len(truth_paths) == len(threshold_paths) == len(timing_paths) == 2
    assert not list(tmp_path.iterdir())
    assert [item[1] for item in captured] == [
        "tree_truth_pulse_aligned_null",
        "threshold_scores_pvalues_power",
        "timing_endpoints_two_epoch_diagnostic",
    ]
    for figure, stem in captured:
        assert tuple(figure.get_size_inches()) == pytest.approx((11, 8.5))
        assert figure.legends[0]._ncols >= 3
        renderer_canvas = FigureCanvasAgg(figure)
        renderer_canvas.draw()
        canvas = figure.bbox
        for artist in figure.findobj(match=Text):
            if (
                not artist.get_visible()
                or not artist.get_text().strip()
                or artist.get_clip_on()
            ):
                continue
            bounds = artist.get_window_extent(renderer_canvas.get_renderer())
            label = f"{stem}: {artist.get_text()!r}"
            assert bounds.x0 >= canvas.x0 - 2, label
            assert bounds.y0 >= canvas.y0 - 2, label
            assert bounds.x1 <= canvas.x1 + 2, label
            assert bounds.y1 <= canvas.y1 + 2, label
    histogram = captured[0][0]
    annotation_text = " ".join(
        text.get_text() for axis in histogram.axes for text in axis.texts
    )
    assert "Reference tails" in annotation_text
    assert "AF-matched tails" in annotation_text
    assert "AA n=3" in annotation_text
    assert "Pulse-aligned display cutoff" in histogram._suptitle.get_text()
    histogram_legend = {text.get_text() for text in histogram.legends[0].texts}
    assert "Selected ref/ref" in histogram_legend
    assert "Selected alt/alt" in histogram_legend
    timing_figure = captured[2][0]
    assert "not a trajectory" in timing_figure.axes[0].get_title().lower()
    timing_canvas = FigureCanvasAgg(timing_figure)
    timing_canvas.draw()
    event_labels = {
        "Origin\n60 kya",
        "Pulse\n56.8 kya",
        "Han split\n50.4 kya",
        "Selection ends\n44.5 kya",
        "Present",
    }
    event_boxes = [
        Text.get_window_extent(text, timing_canvas.get_renderer())
        for text in timing_figure.axes[0].texts
        if text.get_text() in event_labels
    ]
    assert len(event_boxes) == len(event_labels)
    for left_index, left_box in enumerate(event_boxes):
        for right_box in event_boxes[left_index + 1 :]:
            assert not left_box.overlaps(right_box)
    timing_legend = {text.get_text() for text in timing_figure.legends[0].texts}
    assert "Each observed 500-pool AF (filled)" in timing_legend
    assert "Median observed 100-panel AF (open)" in timing_legend


def test_expected_inventory_has_exactly_four_png_pdf_pairs():
    expected = clear._expected_output_paths()
    assert {path for path in expected if path.endswith(".png")} == {
        f"{stem}.png" for stem in clear.PLOT_STEMS
    }
    assert {path for path in expected if path.endswith(".pdf")} == {
        f"{stem}.pdf" for stem in clear.PLOT_STEMS
    }
    assert set(clear.TEXT_OUTPUTS).issubset(expected)


def test_cached_output_exact_inventory_and_tamper_detection(tmp_path, monkeypatch):
    output = tmp_path / "output"
    output.mkdir()
    table = output / "unit_scores.tsv"
    table.write_text("value\n1\n", encoding="utf-8")
    report = output / "RUN_RESULTS.md"
    report.write_text("results\n", encoding="utf-8")
    monkeypatch.setattr(
        clear, "_expected_output_paths", lambda: {"unit_scores.tsv", "RUN_RESULTS.md"}
    )
    contract = {"schema": "test", "parameters": {"custom": True}}
    outputs = {
        path.name: clear._output_record(path, output) for path in (table, report)
    }
    completion = {
        "schema": clear.SCHEMA_VERSION,
        "status": "complete",
        "contract": contract,
        "contract_sha256": clear._canonical_sha256(contract),
        "outputs": outputs,
        "interpretation": clear._completion_interpretation(),
    }
    (output / "completion.json").write_text(
        json.dumps(completion, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    result = clear._validate_cached_output(output, completion, contract)
    assert result["status"] == "verified_cached"

    table.write_text("value\n2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="hash"):
        clear._validate_cached_output(output, completion, contract)


def test_custom_paths_resolve_and_failed_write_removes_atomic_temporary(
    tmp_path, monkeypatch
):
    root = tmp_path / "repo"
    root.mkdir()
    custom_output = tmp_path / "custom" / "clear"
    captured = {}
    prepared = clear.PreparedInputs(
        unit_scores=pd.DataFrame(),
        selected_pointwise=pd.DataFrame(),
        selected_minp=pd.DataFrame(),
        power=pd.DataFrame(),
        endpoints=pd.DataFrame(),
        cache_audit=pd.DataFrame(),
        comparison=pd.DataFrame(),
        contract={"schema": "synthetic"},
    )

    def fake_prepare(repo_root, **kwargs):
        captured["root"] = repo_root
        captured.update(kwargs)
        return prepared

    monkeypatch.setattr(clear, "_prepare_inputs", fake_prepare)
    monkeypatch.setattr(
        clear,
        "_write_results_bundle",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("injected")),
    )
    with pytest.raises(RuntimeError, match="injected"):
        clear._analyze_clear_results_impl(
            root,
            campaign_root="custom_campaign",
            output_dir=custom_output,
            gamma_analysis_dir="custom_gamma",
            truth_analysis_dir="custom_truth",
            survival_results_dir="custom_survival",
            two_epoch_truth_metrics="custom_truth_metrics.json",
            two_epoch_gamma_metrics="custom_gamma_metrics.json",
            two_epoch_gamma_results="custom_results.md",
            workers=4,
            verify_only=False,
        )

    assert captured["campaign_root"] == (root / "custom_campaign").resolve()
    assert captured["survival_results_dir"] == (root / "custom_survival").resolve()
    assert not custom_output.exists()
    assert not custom_output.with_name(
        custom_output.name + f".tmp.{clear.os.getpid()}"
    ).exists()


def test_cli_exposes_analyze_verify_only_and_worker_bounds():
    parser = clear._parser()
    action = next(item for item in parser._actions if item.dest == "action")
    assert tuple(action.choices) == ("analyze", "verify")
    assert "simulate" not in parser.format_help().lower()
    assert "decode" not in parser.format_help().lower()
    assert clear._validate_workers(4) == 4
    for bad in (0, 25, True, 4.0):
        with pytest.raises(ValueError, match="workers"):
            clear._validate_workers(bad)
