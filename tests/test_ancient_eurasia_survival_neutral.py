"""Focused contracts for the AncientEurasia survival-only neutral bank."""

from __future__ import annotations

import copy
import json
from pathlib import Path

import msprime
import numpy as np
import pandas as pd
import pytest
import stdpopsim
import tskit
from gamma_smc_aou import ancient_eurasia_survival_neutral as workflow

REPO_ROOT = Path(__file__).resolve().parents[1]
_MUTATION_CALLS: list[dict] = []


@pytest.fixture(autouse=True)
def _bind_decoder_contract_to_current_test_build(monkeypatch) -> None:
    """Keep test plans strict without assuming compiler-identical binaries."""

    decoder = REPO_ROOT / "bin/gamma_smc"
    if decoder.is_file():
        monkeypatch.setattr(
            workflow, "EXPECTED_DECODER_SHA256", workflow._sha256_file(decoder)
        )


def test_binary_record_requires_the_exact_supplied_digest(tmp_path: Path) -> None:
    binary = tmp_path / "gamma_smc"
    binary.write_bytes(b"test binary")
    observed = workflow._sha256_file(binary)
    with pytest.raises(ValueError, match="binary SHA-256 differs"):
        workflow._binary_record(
            binary, tmp_path, expected_sha256="0" * 64, label="test"
        )
    assert workflow._binary_record(
        binary, tmp_path, expected_sha256=observed, label="test"
    )["sha256"] == observed


def _fake_sim_mutations(
    tree_sequence,
    rate=None,
    *,
    random_seed=None,
    model=None,
    start_time=None,
    end_time=None,
    discrete_genome=None,
    keep=None,
    record_provenance=True,
):
    _MUTATION_CALLS.append(
        {
            "tree_sequence": tree_sequence,
            "rate": rate,
            "random_seed": random_seed,
            "model": model,
            "start_time": start_time,
            "end_time": end_time,
            "discrete_genome": discrete_genome,
            "keep": keep,
            "record_provenance": record_provenance,
        }
    )
    if tree_sequence == "raise":
        raise RuntimeError("synthetic mutation failure")
    return tree_sequence


def _pair_table() -> pd.DataFrame:
    rows = []
    for index in range(workflow.DEFAULT_PANEL_DIPLOIDS):
        rows.append(
            {
                "vcf_diploid_index": index,
                "tree_sequence_individual_id": 1000 + index,
                "sample_node_0": 2 * index,
                "sample_node_1": 2 * index + 1,
                "gamma_smc_haplotype_0": 2 * index,
                "gamma_smc_haplotype_1": 2 * index + 1,
                "focal_position_0based": workflow.FOCAL_POSITION_BP,
                "focal_selected_allele_index": 1,
                "focal_selected_allele_count": 0,
                "genotype_class": "hom_ref",
            }
        )
    return pd.DataFrame(rows)


def _endpoint(unit: dict, *, pool_alt: int = 0, sample_alt: int = 0) -> dict:
    chosen = sorted(
        int(value)
        for value in np.random.default_rng(int(unit["panel_seed"])).choice(
            workflow.DEFAULT_POOL_DIPLOIDS,
            size=workflow.DEFAULT_PANEL_DIPLOIDS,
            replace=False,
        )
    )
    return {
        "final_census_alt_count": 1,
        "final_census_total_count": 2520,
        "final_census_af": 1 / 2520,
        "final_census_segregating": True,
        "pool_diploids": 500,
        "pool_alt_count": pool_alt,
        "pool_total_count": 1000,
        "pool_af": pool_alt / 1000,
        "pool_detected": pool_alt > 0,
        "pool_fixed": pool_alt == 1000,
        "sample_diploids": 100,
        "sample_alt_count": sample_alt,
        "sample_total_count": 200,
        "sample_af": sample_alt / 200,
        "sample_detected": sample_alt > 0,
        "sample_fixed": sample_alt == 200,
        "panel_seed": int(unit["panel_seed"]),
        "chosen_panel_rows_zero_based": chosen,
    }


def _minimal_plan() -> dict:
    return {
        "contract_sha256": "a" * 64,
        "contract": {"source_hashes": {"producer.py": "b" * 64}},
    }


def _decoder_run_record(
    directory: Path,
    decoder: Path,
    unit: dict,
    raw_inputs: dict,
) -> dict:
    tree = directory / "simulation.trees"
    return {
        "schema": workflow.DECODER_RUN_SCHEMA,
        "unit_id": unit["unit_id"],
        "raw_inputs": raw_inputs,
        "run": {
            "command": workflow._expected_decoder_command(decoder, directory),
            "input_path": str(tree),
            "input_format": "trees",
            "tree_sequence_vcf_producer_command": workflow._expected_producer_command(
                tree
            ),
            "stdout": "",
            "stderr": "",
            "decode_seconds": 1.25,
            "stride_bp": 10_000,
            "cache_size_bp": 1_000,
            "n_output_positions": 1_000,
            "pairs_manifest": None,
            "n_pairs_recorded": None,
        },
    }


def test_frozen_schema_names_are_additive() -> None:
    assert workflow.STUDY_ID == "ancient_eurasia_single_origin_survival"
    assert workflow.PLAN_SCHEMA.endswith("plan/v2")
    assert workflow.SIMULATION_SCHEMA.endswith("simulation/v2")
    assert workflow.DECODE_SCHEMA.endswith("decode/v2")
    assert workflow.RESULTS_SCHEMA.endswith("results/v2")
    assert str(workflow.DEFAULT_STUDY_ROOT).endswith(
        "ancient_eurasia_single_origin_survival_v2"
    )
    assert "python/gamma_smc_aou/han_s001_tmrca_simulation.py" in workflow.SOURCE_PATHS


def _overlay_maps() -> tuple[msprime.RateMap, msprime.RateMap, msprime.RateMap]:
    focal = workflow.FOCAL_POSITION_BP
    length = workflow.SEQUENCE_LENGTH_BP
    mu = workflow.MUTATION_RATE
    ordinary = msprime.RateMap(position=[0, focal, focal + 1, length], rate=[mu, 0, mu])
    focal_dfe = msprime.RateMap(position=[0, focal, focal + 1, length], rate=[0, mu, 0])
    recap = msprime.RateMap(position=[0, length], rate=[mu])
    return ordinary, focal_dfe, recap


def _tree_with_500_individual_locations(
    *, wrong_width: bool = False, nonfinite: bool = False
):
    tables = tskit.TableCollection(sequence_length=workflow.SEQUENCE_LENGTH_BP)
    tables.metadata_schema = tskit.MetadataSchema.permissive_json()
    tables.metadata = {"location_test": "preserve"}
    population = tables.populations.add_row()
    for index in range(workflow.DEFAULT_POOL_DIPLOIDS):
        width = 2 if wrong_width and index == 0 else 3
        location = np.arange(width, dtype=float) + index + 0.25
        if nonfinite and index == 0:
            location[:] = [np.nan, np.inf, -np.inf]
        tables.individuals.add_row(flags=index % 2, location=location)
    tables.nodes.add_row(
        flags=tskit.NODE_IS_SAMPLE, time=0, population=population, individual=0
    )
    tables.provenances.add_row(record='{"location_test":true}', timestamp="test")
    return tables.tree_sequence()


def test_individual_locations_are_canonical_and_only_location_changes() -> None:
    original = _tree_with_500_individual_locations()
    normalized, receipt = workflow._canonicalize_individual_locations(original)
    contract = workflow._individual_location_contract()
    workflow._validate_individual_location_receipt(receipt, contract, tree=normalized)
    assert receipt["before"]["location_all_finite"] is True
    assert receipt["before"]["location_all_zero"] is False
    assert receipt["after"]["location_all_zero"] is True
    assert receipt["after"]["location_record"]["sha256"] == (
        "ff6698a6e831ffcf47af2fed388ffc262f319e72b26cd140929d1e19b1246ad4"
    )
    assert set(receipt["proof"].values()) == {True}
    original_tables = original.dump_tables()
    normalized_tables = normalized.dump_tables()
    original_tables.nodes.assert_equals(normalized_tables.nodes)
    original_tables.provenances.assert_equals(normalized_tables.provenances)
    assert original_tables.metadata == normalized_tables.metadata
    assert np.array_equal(
        original_tables.individuals.flags, normalized_tables.individuals.flags
    )
    assert np.array_equal(
        original_tables.individuals.location_offset,
        normalized_tables.individuals.location_offset,
    )
    renormalized, second = workflow._canonicalize_individual_locations(normalized)
    assert np.array_equal(
        normalized_tables.individuals.location,
        renormalized.dump_tables().individuals.location,
    )
    assert second["before"]["location_all_zero"] is True


def test_individual_location_layout_and_receipt_tamper_are_rejected() -> None:
    with pytest.raises(ValueError, match="location layout differs"):
        workflow._canonicalize_individual_locations(
            _tree_with_500_individual_locations(wrong_width=True)
        )
    normalized, receipt = workflow._canonicalize_individual_locations(
        _tree_with_500_individual_locations()
    )
    tampered = copy.deepcopy(receipt)
    tampered["after"]["location_record"]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="not exact zeros"):
        workflow._validate_individual_location_receipt(
            tampered, receipt["contract"], tree=normalized
        )
    stored_tamper = copy.deepcopy(receipt)
    stored_tamper["stored_tree"]["other_tables_ts_metadata_reference_sha256"] = "1" * 64
    with pytest.raises(ValueError, match="after state differs"):
        workflow._validate_individual_location_receipt(
            stored_tamper, receipt["contract"], tree=normalized
        )
    for field in (
        "other_tables_ts_metadata_reference_sha256",
        "all_except_individual_location_sha256",
    ):
        paired_tamper = copy.deepcopy(receipt)
        paired_tamper["before"][field] = "a" * 64
        paired_tamper["after"][field] = "a" * 64
        with pytest.raises(ValueError, match="after state differs"):
            workflow._validate_individual_location_receipt(
                paired_tamper, receipt["contract"], tree=normalized
            )


def test_individual_location_nonfinite_garbage_is_discarded() -> None:
    normalized, receipt = workflow._canonicalize_individual_locations(
        _tree_with_500_individual_locations(nonfinite=True)
    )
    assert receipt["before"]["location_finite_count"] == 1497
    assert receipt["before"]["location_nan_count"] == 1
    assert receipt["before"]["location_positive_infinity_count"] == 1
    assert receipt["before"]["location_negative_infinity_count"] == 1
    assert receipt["after"]["location_finite_count"] == 1500
    assert receipt["after"]["location_nan_count"] == 0
    workflow._validate_individual_location_receipt(
        receipt, receipt["contract"], tree=normalized
    )


def test_individual_location_receipt_binds_final_sample_flags() -> None:
    raw = _tree_with_500_individual_locations()
    tables = raw.dump_tables()
    nodes = tables.nodes
    nodes.set_columns(
        flags=np.zeros_like(nodes.flags),
        time=nodes.time,
        population=nodes.population,
        individual=nodes.individual,
        metadata=nodes.metadata,
        metadata_offset=nodes.metadata_offset,
    )
    panel_state = tables.tree_sequence()
    normalized, receipt = workflow._canonicalize_individual_locations(panel_state)
    workflow._validate_individual_location_receipt(
        receipt, receipt["contract"], tree=normalized
    )
    altered_tables = normalized.dump_tables()
    nodes = altered_tables.nodes
    nodes.set_columns(
        flags=np.full_like(nodes.flags, tskit.NODE_IS_SAMPLE),
        time=nodes.time,
        population=nodes.population,
        individual=nodes.individual,
        metadata=nodes.metadata,
        metadata_offset=nodes.metadata_offset,
    )
    with pytest.raises(ValueError, match="stored tree state differs"):
        workflow._validate_individual_location_receipt(
            receipt,
            receipt["contract"],
            tree=altered_tables.tree_sequence(),
        )


def test_scoped_overlay_guard_masks_all_positive_layers_and_preserves_args(
    monkeypatch,
) -> None:
    monkeypatch.setattr(workflow, "EXPECTED_OVERLAY_CALL_COUNT", None)
    monkeypatch.setattr(workflow, "EXPECTED_OVERLAY_PATCHED_COUNT", None)
    original = workflow.msprime.sim_mutations
    monkeypatch.setattr(workflow.msprime, "sim_mutations", _fake_sim_mutations)
    _MUTATION_CALLS.clear()
    ordinary, focal_dfe, recap = _overlay_maps()
    models = [msprime.SLiMMutationModel(type=index) for index in (1, 2, 3)]
    with workflow._scoped_focal_overlay_patch() as receipt:
        workflow.msprime.sim_mutations(
            "ordinary",
            rate=ordinary,
            model=models[0],
            end_time=12.0,
            keep=True,
            random_seed=101,
        )
        workflow.msprime.sim_mutations(
            "focal",
            rate=focal_dfe,
            model=models[1],
            end_time=12.0,
            keep=True,
            random_seed=102,
        )
        workflow.msprime.sim_mutations(
            "recap",
            rate=recap,
            model=models[2],
            start_time=12.0,
            keep=True,
            random_seed=103,
        )
    assert workflow.msprime.sim_mutations is _fake_sim_mutations
    assert receipt["intercepted_call_count"] == 3
    assert receipt["patched_call_count"] == 2
    assert [call["rate_argument_replaced"] for call in receipt["calls"]] == [
        False,
        True,
        True,
    ]
    assert _MUTATION_CALLS[0]["rate"] is ordinary
    assert _MUTATION_CALLS[0]["model"] is models[0]
    assert _MUTATION_CALLS[0]["random_seed"] == 101
    assert _MUTATION_CALLS[0]["keep"] is True
    for observed, model in zip(_MUTATION_CALLS, models, strict=True):
        assert observed["model"] is model
        assert observed["rate"].get_rate(workflow.FOCAL_POSITION_BP) == 0.0
    focal_passed = workflow._rate_map_metrics(_MUTATION_CALLS[1]["rate"])
    assert focal_passed["total_integral"] == 0.0
    recap_passed = workflow._rate_map_metrics(_MUTATION_CALLS[2]["rate"])
    recap_original = workflow._rate_map_metrics(recap)
    assert recap_passed["left_flank_integral"] == recap_original["left_flank_integral"]
    assert (
        recap_passed["right_flank_integral"] == recap_original["right_flank_integral"]
    )
    monkeypatch.setattr(workflow.msprime, "sim_mutations", original)


def test_scoped_overlay_guard_does_not_patch_wrong_or_zero_calls(monkeypatch) -> None:
    monkeypatch.setattr(workflow, "EXPECTED_OVERLAY_CALL_COUNT", None)
    monkeypatch.setattr(workflow, "EXPECTED_OVERLAY_PATCHED_COUNT", None)
    monkeypatch.setattr(workflow.msprime, "sim_mutations", _fake_sim_mutations)
    _MUTATION_CALLS.clear()
    ordinary, _focal_dfe, recap = _overlay_maps()
    model = msprime.SLiMMutationModel(type=1)
    with workflow._scoped_focal_overlay_patch() as receipt:
        workflow.msprime.sim_mutations("scalar", rate=1e-8, model=model)
        workflow.msprime.sim_mutations("zero", rate=ordinary, model=model)
        workflow.msprime.sim_mutations("eligible", rate=recap, model=model)
    assert receipt["patched_call_count"] == 1
    assert receipt["calls"][0]["rate_argument_replaced"] is False
    assert receipt["calls"][1]["rate_argument_replaced"] is False
    assert _MUTATION_CALLS[0]["rate"] == 1e-8
    assert _MUTATION_CALLS[1]["rate"] is ordinary


def test_scoped_overlay_guard_restores_after_exception(monkeypatch) -> None:
    monkeypatch.setattr(workflow.msprime, "sim_mutations", _fake_sim_mutations)
    _MUTATION_CALLS.clear()
    _ordinary, _focal_dfe, recap = _overlay_maps()
    with (
        pytest.raises(RuntimeError, match="synthetic mutation failure"),
        workflow._scoped_focal_overlay_patch() as receipt,
    ):
        workflow.msprime.sim_mutations(
            "raise",
            rate=recap,
            model=msprime.SLiMMutationModel(type=1),
        )
    assert workflow.msprime.sim_mutations is _fake_sim_mutations
    assert receipt["callable_restored"] is True
    assert receipt["calls"][0]["call_succeeded"] is False


def test_overlay_receipt_tamper_is_rejected(monkeypatch) -> None:
    monkeypatch.setattr(workflow, "EXPECTED_OVERLAY_CALL_COUNT", None)
    monkeypatch.setattr(workflow, "EXPECTED_OVERLAY_PATCHED_COUNT", None)
    monkeypatch.setattr(workflow.msprime, "sim_mutations", _fake_sim_mutations)
    _MUTATION_CALLS.clear()
    _ordinary, _focal_dfe, recap = _overlay_maps()
    with workflow._scoped_focal_overlay_patch() as receipt:
        workflow.msprime.sim_mutations(
            "eligible",
            rate=recap,
            model=msprime.SLiMMutationModel(type=1),
        )
    tampered = copy.deepcopy(receipt)
    tampered["calls"][0]["passed_rate_map"]["rates"][1] = workflow.MUTATION_RATE
    with pytest.raises(ValueError, match="rate-map receipt differs"):
        workflow._validate_focal_overlay_patch_receipt(tampered, receipt["contract"])


def test_frozen_overlay_inventory_is_exact_two_calls_one_patch(monkeypatch) -> None:
    monkeypatch.setattr(workflow.msprime, "sim_mutations", _fake_sim_mutations)
    _MUTATION_CALLS.clear()
    ordinary, _focal_dfe, recap = _overlay_maps()
    tick = 29_665.0
    models = [
        msprime.SLiMMutationModel(type=0, next_id=7),
        msprime.SLiMMutationModel(type=2, next_id=9),
    ]
    with workflow._scoped_focal_overlay_patch() as receipt:
        workflow.msprime.sim_mutations(
            "ordinary",
            rate=ordinary,
            model=models[0],
            end_time=tick,
            keep=True,
            random_seed=2_460_052_835,
        )
        workflow.msprime.sim_mutations(
            "recap",
            rate=recap,
            model=models[1],
            start_time=tick,
            keep=True,
            random_seed=2_048_104_378,
        )
    assert receipt["contract"]["inventory_freeze_state"] == (
        "frozen_after_rep051_smoke"
    )
    assert receipt["intercepted_call_count"] == 2
    assert receipt["patched_call_count"] == 1
    assert [call["model_type"] for call in receipt["calls"]] == [0, 2]
    assert [call["model_next_id"] for call in receipt["calls"]] == [7, 9]
    assert receipt["calls"][0]["end_time"] == receipt["calls"][1]["start_time"]
    assert receipt["calls"][0]["rate_argument_replaced"] is False
    assert receipt["calls"][1]["rate_argument_replaced"] is True
    tampered = copy.deepcopy(receipt)
    tampered["calls"][1]["model_type"] = 1
    with pytest.raises(ValueError, match="exact call inventory shape"):
        workflow._validate_focal_overlay_patch_receipt(tampered, receipt["contract"])


def test_exact_six_event_schedule_and_transform_counts() -> None:
    events, transform = workflow.build_neutral_event_schedule()
    assert len(events) == 6
    assert not any(
        isinstance(event, stdpopsim.ChangeMutationFitness) for event in events
    )
    assert transform["removed_change_mutation_fitness_count"] == 2
    assert transform["removed_terminal_han_af_band_count"] == 2
    assert transform["added_terminal_han_strict_upper_count"] == 1
    assert transform["nominal_pulse_generations_ago"] == 2272
    assert transform["realized_pulse_generations_ago_q5"] == 2270
    assert transform["nominal_han_split_generations_ago"] == 2016
    assert transform["realized_han_split_generations_ago_q5"] == 2015
    schedule = transform["event_schedule"]
    assert [record["population"] for record in schedule] == [
        "Neanderthal",
        "Neanderthal",
        "Neanderthal",
        "Loschbour",
        "Han",
        "Han",
    ]
    assert schedule[0]["time"] == {
        "generations_ago": 2400.0,
        "semantics": "exact_generation",
    }
    assert schedule[1]["start_time"]["semantics"] == (
        "generation_after_in_forward_time"
    )
    assert schedule[1]["end_time"]["generations_ago"] == 2273
    assert (schedule[2]["operator"], schedule[2]["allele_frequency"]) == (
        ">=",
        1e-9,
    )
    assert schedule[3]["start_time"]["generations_ago"] == 2272
    assert schedule[3]["end_time"]["generations_ago"] == 2015
    assert schedule[4]["start_time"]["generations_ago"] == 2015
    assert (schedule[5]["operator"], schedule[5]["allele_frequency"]) == ("<", 1.0)


def test_units_seeds_and_simulation_law_have_no_observational_gate() -> None:
    units = workflow._execution_units()
    assert len(units) == 100
    assert tuple(units.columns) == workflow.EXECUTION_COLUMNS
    assert units.iloc[0]["unit_id"].endswith("rep000")
    assert units.iloc[-1]["unit_id"].endswith("rep099")
    first = units.iloc[0].to_dict()
    assert first["seed"] == workflow._stable_seed(first["unit_id"])
    assert first["panel_seed"] == workflow._stable_seed(f"{first['unit_id']}:panel")
    assert (
        not units[
            [
                "pool_frequency_conditioning",
                "terminal_af_target_or_band_conditioning",
                "sample_detection_conditioning",
            ]
        ]
        .to_numpy(dtype=bool)
        .any()
    )
    assert workflow.SIMULATION_LAW["engine"] == "stdpopsim_slim_forward"
    assert workflow.SIMULATION_LAW["sequence_scope"] == "direct_full_10mb"
    assert not any(
        workflow.SIMULATION_LAW[key]
        for key in (
            "branch_enumeration",
            "ancestry_scan",
            "branch_span_weighting",
            "cropping",
            "posthoc_causal_mutation_placement",
            "external_pool_frequency_rejection",
            "external_sample_detection_rejection",
            "external_terminal_target_af_rejection",
        )
    )


def test_scoped_end_hook_changes_only_end_and_restores() -> None:
    original_object = workflow.slim_engine._slim_functions
    original = str(original_object)
    expected = workflow._census_end_patch_contract(13)
    with workflow._scoped_census_end_patch(13) as observed:
        patched = str(workflow.slim_engine._slim_functions)
        assert observed == expected
        assert patched != original
        assert patched.count("survival_neutral_final_census_alt_count") == 1
        assert "targets.addNewDrawnMutation(mut_type, pos);" in patched
        assert observed["mutation_draw_helper_changed"] is False
    assert workflow.slim_engine._slim_functions is original_object


def test_census_metadata_reader_requires_one_exact_path() -> None:
    user = {
        "survival_neutral_final_census_alt_count": [4],
        "survival_neutral_final_census_total_count": [2520],
        "survival_neutral_final_census_af": [4 / 2520],
        "survival_neutral_final_census_segregating": [1],
        "survival_neutral_final_census_fixed": [0],
    }
    metadata = {"SLiM": {"user_metadata": user}}
    assert (
        workflow._metadata_integer(metadata, "survival_neutral_final_census_alt_count")
        == 4
    )
    duplicate = {"SLiM": {"user_metadata": user}, "duplicate": copy.deepcopy(user)}
    with pytest.raises(ValueError, match="uniquely scoped"):
        workflow._slim_user_metadata(duplicate)


def test_endpoint_accepts_zero_or_fixed_observational_samples_but_not_bad_panel() -> (
    None
):
    unit = workflow._execution_units().iloc[0].to_dict()
    endpoint = _endpoint(unit, pool_alt=0, sample_alt=200)
    workflow._validate_endpoint(endpoint, unit)
    assert endpoint["pool_detected"] is False
    assert endpoint["sample_fixed"] is True
    tampered = copy.deepcopy(endpoint)
    tampered["chosen_panel_rows_zero_based"][0] = tampered[
        "chosen_panel_rows_zero_based"
    ][1]
    with pytest.raises(ValueError, match="deterministic uniform"):
        workflow._validate_endpoint(tampered, unit)


def test_plan_is_exact_cached_and_tamper_detected(tmp_path: Path) -> None:
    study = tmp_path / "study"
    path = workflow.write_plan(REPO_ROOT, root=study)
    payload = workflow.load_plan(REPO_ROOT, root=study)
    assert path == study / "plan/study_plan.json"
    assert set(json.loads(path.read_text())) == {
        "schema",
        "status",
        "created_utc",
        "contract",
        "contract_sha256",
        "outputs",
    }
    assert len(payload["units"]) == 100
    assert payload["contract"]["simulation_law"] == workflow.SIMULATION_LAW
    overlay = payload["contract"]["focal_overlay_patch"]
    assert overlay["inventory_freeze_state"] == "frozen_after_rep051_smoke"
    assert overlay["expected_call_count"] == 2
    assert overlay["expected_patched_call_count"] == 1
    assert [row["model_type"] for row in overlay["expected_call_inventory"]] == [
        0,
        2,
    ]
    assert overlay["upstream_overlay_method"].endswith("._recap_and_rescale")
    assert len(overlay["upstream_overlay_method_source_sha256"]) == 64
    location = payload["contract"]["individual_location_canonicalization"]
    assert location["expected_individual_count"] == 500
    assert location["required_location_width_per_individual"] == 3
    assert location["canonical_flat_location_record"]["sha256"] == (
        "ff6698a6e831ffcf47af2fed388ffc262f319e72b26cd140929d1e19b1246ad4"
    )
    assert location["canonical_location_offset_record"]["sha256"] == (
        "5d348c3a0a5849e32e7c21925d7532df9762cf2fa2435d1ff297b9c0f1602a9d"
    )
    assert location["mutation_api_source_sha256"] == (
        "a41ee8aecaeb0dc3829403bba3df77caf9e4748ee48ed595874bab1fa2ffa5f3"
    )
    assert workflow.write_plan(REPO_ROOT, root=study) == path
    execution = study / "plan/execution_units.tsv"
    execution.write_text(execution.read_text() + "\n", encoding="utf-8")
    with pytest.raises((ValueError, AssertionError)):
        workflow.load_plan(REPO_ROOT, root=study)


def test_force_plan_refuses_after_work_exists(tmp_path: Path) -> None:
    study = tmp_path / "study"
    workflow.write_plan(REPO_ROOT, root=study)
    work = study / "work/unit"
    work.mkdir(parents=True)
    (work / "receipt.txt").write_text("owned work", encoding="utf-8")
    with pytest.raises(ValueError, match="work/results exist"):
        workflow.write_plan(REPO_ROOT, root=study, force=True)


def test_pair_file_is_exact_100_within_diploid_pairs(tmp_path: Path) -> None:
    path = tmp_path / "overall.pairs.tsv"
    workflow._atomic_pair_file(path, _pair_table())
    observed = workflow._read_pairs(path)
    assert len(observed) == 100
    assert np.array_equal(observed.to_numpy(), np.arange(200).reshape(100, 2))


def test_simulate_unit_selector_keeps_100_unit_plan(
    monkeypatch, tmp_path: Path
) -> None:
    study = tmp_path / workflow.INSTRUMENTED_SMOKE_ROOT_BASENAME
    plan_path = workflow.write_plan(REPO_ROOT, root=study)
    plan = workflow.load_plan(REPO_ROOT, root=study)
    selected_id = plan["units"][51]["unit_id"]

    def fake_simulate(payload):
        return {
            **payload["unit"],
            "status": "complete",
            "error": "",
            "final_census_af": 0.01,
            "simulation_completion_sha256": "c" * 64,
        }

    monkeypatch.setattr(workflow, "_simulate_unit", fake_simulate)
    status = workflow.simulate_study(
        REPO_ROOT, root=study, workers=1, unit_ids=[selected_id]
    )
    assert len(status) == 1
    assert status.iloc[0]["unit_id"] == selected_id
    assert len(pd.read_csv(plan_path.parent / "execution_units.tsv", sep="\t")) == 100


def test_provisional_overlay_inventory_fails_closed_outside_rep051_smoke(
    monkeypatch, tmp_path: Path
) -> None:
    monkeypatch.setattr(workflow, "EXPECTED_OVERLAY_CALL_COUNT", None)
    monkeypatch.setattr(workflow, "EXPECTED_OVERLAY_PATCHED_COUNT", None)

    def forbidden(_payload):
        raise AssertionError("worker must not launch")

    monkeypatch.setattr(workflow, "_simulate_unit", forbidden)
    wrong_root = tmp_path / "not_the_instrumented_smoke"
    workflow.write_plan(REPO_ROOT, root=wrong_root)
    with pytest.raises(RuntimeError, match="allows only an explicit rep051"):
        workflow.simulate_study(
            REPO_ROOT,
            root=wrong_root,
            workers=1,
            unit_ids=[workflow.INSTRUMENTED_SMOKE_UNIT_ID],
        )
    smoke_root = tmp_path / workflow.INSTRUMENTED_SMOKE_ROOT_BASENAME
    workflow.write_plan(REPO_ROOT, root=smoke_root)
    with pytest.raises(RuntimeError, match="allows only an explicit rep051"):
        workflow.simulate_study(REPO_ROOT, root=smoke_root, workers=1)


def test_decode_raw_inputs_command_and_standalone_tamper(tmp_path: Path) -> None:
    repo = REPO_ROOT
    root = tmp_path / "study"
    unit = workflow._execution_units().iloc[0].to_dict()
    directory = root / "work" / unit["unit_id"]
    decoded = directory / "decoded"
    decoded.mkdir(parents=True)
    pair_table = _pair_table()
    pair_table.to_csv(directory / "sample_manifest.tsv", sep="\t", index=False)
    (directory / "simulation.trees").write_bytes(b"mock-tree-input")
    endpoint = _endpoint(unit)
    (directory / "simulation_complete.json").write_text(
        json.dumps({"endpoint": endpoint}) + "\n", encoding="utf-8"
    )
    workflow._atomic_pair_file(decoded / "overall.pairs.tsv", pair_table)
    plan = _minimal_plan()
    decoder = (repo / "bin/gamma_smc").resolve()
    contract = workflow._decode_contract(
        repo, root, plan, unit, decoder, decoded / "overall.pairs.tsv"
    )
    assert set(contract["raw_inputs"]) == {
        "simulation_tree",
        "sample_manifest",
        "overall_pairs",
    }
    assert contract["decoder_settings"]["pair_selector"] == (
        "explicit_within_individual_pairs"
    )
    workflow._atomic_json(decoded / "decode_contract.json", contract)
    probabilities = {
        f"mean_p_lt_{threshold}": index / 10
        for index, threshold in enumerate(workflow.TMRCA_THRESHOLDS_YEARS, start=1)
    }
    summary = pd.DataFrame(
        [
            {
                "position_0based": workflow.FOCAL_POSITION_BP,
                "n_pairs": 100,
                **probabilities,
            }
        ]
    )
    summary.to_csv(decoded / "overall.summary.tsv", sep="\t", index=False)
    (decoded / "posterior.zst").write_bytes(b"posterior")
    (decoded / "posterior.zst.meta").write_bytes(b"metadata")
    run_record = _decoder_run_record(directory, decoder, unit, contract["raw_inputs"])
    workflow._atomic_json(decoded / "decoder_run.json", run_record)
    scores = workflow._gamma_scores(unit, summary, endpoint)
    scores.to_csv(decoded / "gamma_unit_scores.tsv", sep="\t", index=False)
    row_counts = {
        "overall_pairs": 100,
        "overall_summary": 1,
        "gamma_unit_scores": 7,
    }
    outputs = {
        key: workflow._output_record(decoded / name, rows=row_counts.get(key))
        for key, name in workflow.DECODE_OUTPUT_PATHS.items()
    }
    completion = {
        "schema": workflow.DECODE_SCHEMA,
        "status": "complete",
        "contract": contract,
        "contract_sha256": workflow._canonical_sha256(contract),
        "elapsed_seconds": 1.0,
        "outputs": outputs,
    }
    workflow._atomic_json(decoded / "decode_complete.json", completion)
    validated = workflow._validate_decode_completion(root, repo, plan, unit, decoder)
    assert validated["status"] == "complete"
    tampered = copy.deepcopy(contract)
    tampered["decoder_settings"]["cache_size"] = 999
    workflow._atomic_json(decoded / "decode_contract.json", tampered)
    with pytest.raises(ValueError, match="standalone decode contract"):
        workflow._validate_decode_completion(root, repo, plan, unit, decoder)


def test_decoder_run_rejects_any_command_setting_tamper(tmp_path: Path) -> None:
    unit = workflow._execution_units().iloc[0].to_dict()
    directory = tmp_path / "work" / unit["unit_id"]
    (directory / "decoded").mkdir(parents=True)
    decoder = (REPO_ROOT / "bin/gamma_smc").resolve()
    raw_inputs = {"simulation_tree": {}, "sample_manifest": {}, "overall_pairs": {}}
    record = _decoder_run_record(directory, decoder, unit, raw_inputs)
    workflow._validate_decoder_run(
        record,
        unit=unit,
        contract={"raw_inputs": raw_inputs},
        directory=directory,
        decoder_bin=decoder,
    )
    tampered = copy.deepcopy(record)
    command = tampered["run"]["command"]
    command[command.index("--cache_size") + 1] = "999"
    with pytest.raises(ValueError, match="every frozen setting"):
        workflow._validate_decoder_run(
            tampered,
            unit=unit,
            contract={"raw_inputs": raw_inputs},
            directory=directory,
            decoder_bin=decoder,
        )
    extra_selector = copy.deepcopy(record)
    extra_selector["run"]["command"].append("--only_within")
    with pytest.raises(ValueError, match="every frozen setting"):
        workflow._validate_decoder_run(
            extra_selector,
            unit=unit,
            contract={"raw_inputs": raw_inputs},
            directory=directory,
            decoder_bin=decoder,
        )
    producer_tamper = copy.deepcopy(record)
    producer_tamper["run"]["tree_sequence_vcf_producer_command"][1] = "--version"
    with pytest.raises(ValueError, match="producer command"):
        workflow._validate_decoder_run(
            producer_tamper,
            unit=unit,
            contract={"raw_inputs": raw_inputs},
            directory=directory,
            decoder_bin=decoder,
        )


def test_cli_help(capsys) -> None:
    with pytest.raises(SystemExit) as error:
        workflow.main(["--help"])
    assert error.value.code == 0
    assert "strict present-Han census segregation" in capsys.readouterr().out
