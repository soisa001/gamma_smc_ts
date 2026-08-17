from __future__ import annotations

import copy
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from gamma_smc_aou import cosi2_conditioned_null as conditioned


REPO_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def eas_sizes() -> np.ndarray:
    return conditioned.eas_integer_ne_by_generation(REPO_ROOT)


@pytest.fixture(scope="module")
def small_eas_blocks(eas_sizes: np.ndarray) -> list[dict[str, object]]:
    return [
        conditioned.simulate_eas_proposal_block(
            eas_sizes,
            block_index=index,
            proposals_per_block=64,
        )
        for index in range(4)
    ]


def _synthetic_base_manifest() -> dict[str, object]:
    source_binding = conditioned._frozen_eas_base_source_binding(REPO_ROOT)
    entries = [
        {
            "block_index": index,
            "block_seed": conditioned.EAS_PROPOSAL_BLOCK_SEED_START + index,
            "proposal_count": conditioned.EAS_PROPOSALS_PER_BLOCK,
            "completion_size_bytes": 1,
            "completion_sha256": "a" * 64,
            "proposal_block_size_bytes": 1,
            "proposal_block_sha256": "b" * 64,
            "summary_size_bytes": 1,
            "summary_file_sha256": "c" * 64,
            "summary_payload_sha256": "d" * 64,
        }
        for index in range(conditioned.EAS_BASE_PROPOSAL_BLOCK_COUNT)
    ]
    failure = conditioned._base_bank_failure_record(
        {
            "total_proposals": conditioned.EAS_BASE_TOTAL_PROPOSALS,
            "valid_proposals": 484_118,
            "global_log_max_weight": 10.069685602283542,
            "global_weight_sum_over_max": 49.095072451783,
            "global_effective_sample_size": 2207.3496608931655,
            "global_max_normalized_weight": 1.0 / 49.095072451783,
            "cap_hit_count": 0,
            "proposal_bank_binding_sha256": "e" * 64,
        }
    )
    attestation: dict[str, object] = {
        "status": "exact_match",
        "proposal_kernel_id": conditioned.EAS_PROPOSAL_KERNEL_ID,
        "attested_block_index": 0,
        "attested_block_seed": conditioned.EAS_PROPOSAL_BLOCK_SEED_START,
        "attested_proposal_count": conditioned.EAS_PROPOSALS_PER_BLOCK,
        "frozen_implementation_sha256": conditioned.EAS_BASE_IMPLEMENTATION_SHA256,
        "current_implementation_sha256": conditioned.sha256_file(
            conditioned.MODULE_PATH
        ),
        "numpy_version": np.__version__,
        "array_sha256": {
            name: "f" * 64
            for name in (
                "outcomes",
                "birth_generations",
                "log_forward",
                "log_reverse",
                "log_mutation",
                "log_weights",
            )
        },
    }
    attestation["attestation_payload_sha256"] = conditioned._canonical_sha256(
        attestation
    )
    manifest: dict[str, object] = {
        "schema": conditioned.EAS_BASE_BANK_MANIFEST_SCHEMA,
        "status": "complete",
        "contract": conditioned._eas_base_bank_contract_record(),
        "source_binding": source_binding,
        "source_binding_sha256": conditioned._canonical_sha256(source_binding),
        "blocks": entries,
        "base_failure_diagnostics": failure,
        "kernel_equivalence_attestation": attestation,
    }
    manifest["manifest_payload_sha256"] = conditioned._canonical_sha256(manifest)
    return manifest


def _write_synthetic_production_block(
    destination: Path,
    *,
    block_index: int,
    manifest: dict[str, object],
    plan: dict[str, object],
) -> None:
    count = conditioned.EAS_PROPOSALS_PER_BLOCK
    arrays = {
        "outcomes": np.full(count, "invalid_loss", dtype="U12"),
        "birth_generations": np.full(count, -1, dtype=np.int64),
        "log_forward": np.full(count, np.nan),
        "log_reverse": np.full(count, np.nan),
        "log_mutation": np.full(count, np.nan),
        "log_weights": np.full(count, np.nan),
    }
    block = {
        "schema": conditioned.EAS_PROPOSAL_BLOCK_SCHEMA,
        "block_index": block_index,
        "block_seed": conditioned.EAS_PROPOSAL_BLOCK_SEED_START + block_index,
        "proposal_count": count,
        "present_count": conditioned.EAS_SELECTED_PRESENT_COUNT,
        "max_birth_generation": conditioned.EAS_BRIDGE_MAX_GENERATION,
        **arrays,
    }
    summary = conditioned.validate_eas_proposal_block(block)
    destination.mkdir()
    np.savez_compressed(destination / "proposal_block.npz", **arrays)
    conditioned._atomic_json(destination / "summary.json", summary)
    outputs = {
        name: {
            "size_bytes": (destination / name).stat().st_size,
            "sha256": conditioned.sha256_file(destination / name),
        }
        for name in ("proposal_block.npz", "summary.json")
    }
    completion = {
        "schema": conditioned.EAS_PROPOSAL_BLOCK_SCHEMA,
        "status": "complete",
        "created_utc": "2026-08-17T00:00:00Z",
        "block_index": block_index,
        "block_seed": conditioned.EAS_PROPOSAL_BLOCK_SEED_START + block_index,
        "proposal_count": count,
        "present_count": conditioned.EAS_SELECTED_PRESENT_COUNT,
        "present_chromosomes": conditioned.EAS_PRESENT_CHROMOSOMES,
        "max_birth_generation": conditioned.EAS_BRIDGE_MAX_GENERATION,
        "absorption_draw_generation": conditioned.EAS_BRIDGE_ABSORPTION_GENERATION,
        "source_binding": conditioned._eas_proposal_source_binding(REPO_ROOT),
        "summary_sha256": conditioned._canonical_sha256(summary),
        "outputs": outputs,
    }
    if block_index >= conditioned.EAS_BASE_PROPOSAL_BLOCK_COUNT:
        completion["extension_bank_binding"] = conditioned._eas_extension_block_binding(
            manifest, plan
        )
    conditioned._atomic_json(destination / "completion.json", completion)


def test_exact_endpoint_and_fixed_ledgers_are_disjoint():
    assert conditioned.EAS_SELECTED_PRESENT_COUNT == 14_704
    assert conditioned.EAS_SELECTED_PRESENT_AF == pytest.approx(0.195339692324)
    eas = conditioned.build_eas_unit_ledger()
    han = conditioned.build_han_candidate_ledger()
    assert len(eas) == 100
    assert len(han) == 5_000
    assert eas["seed"].tolist() == list(range(2_026_083_600, 2_026_083_700))
    assert han.iloc[0]["source_kind"] == ("existing_complete_natural_neutral_candidate")
    assert han.iloc[99]["seed"] == 2_026_082_699
    assert han.iloc[100]["seed"] == 2_026_084_600
    assert not set(eas["seed"]).intersection(han["seed"])


def test_fixed_v2_extension_contract_records_the_post_pilot_gate_revision():
    contract = conditioned.EasProposalBankContract()
    contract.validate()
    assert contract.base_block_count == 20
    assert contract.target_block_count == 200
    assert contract.extension_block_count == 180
    assert contract.total_proposals == 10_000_000
    assert contract.adaptive_stopping_allowed is False

    manifest = _synthetic_base_manifest()
    conditioned.validate_eas_base_bank_manifest(REPO_ROOT, manifest)
    provenance = manifest["contract"]
    assert provenance["quality_gates_revised_post_pilot"] is True
    assert "quality_gates_unchanged" not in provenance
    assert provenance["pilot_quality_gate_contract"]["min_weight_sum_over_max"] == 1000
    assert provenance["v2_quality_gate_contract"]["min_weight_sum_over_max"] == 100
    assert (
        "max_normalized_weight_exclusive" in provenance["pilot_quality_gate_contract"]
    )
    assert "max_normalized_weight_inclusive" in provenance["v2_quality_gate_contract"]

    plan = conditioned.build_eas_extension_plan(REPO_ROOT, manifest)
    validation = conditioned.validate_eas_extension_plan(
        REPO_ROOT,
        plan,
        base_manifest=manifest,
    )
    assert validation["target_block_count"] == 200
    assert plan["block_ledger"][19]["required_action"] == (
        "reuse_and_verify_exact_bytes"
    )
    assert plan["block_ledger"][20]["required_action"] == "generate_or_verify"
    assert plan["block_ledger"][-1]["block_index"] == 199
    assert [row["block_seed"] for row in plan["block_ledger"]] == list(
        range(
            conditioned.EAS_PROPOSAL_BLOCK_SEED_START,
            conditioned.EAS_PROPOSAL_BLOCK_SEED_START + 200,
        )
    )
    assert plan["post_pilot_decision"]["adaptive_stopping_allowed"] is False
    assert plan["invalidation"]["reuse_v1_selected_particles"] is False

    with pytest.raises(ValueError, match="fixed v2 200-block"):
        conditioned.build_eas_extension_plan(
            REPO_ROOT,
            manifest,
            target_block_count=201,
        )
    tampered = copy.deepcopy(plan)
    tampered["post_pilot_decision"]["adaptive_stopping_allowed"] = True
    tampered["plan_payload_sha256"] = conditioned._canonical_sha256(
        {key: value for key, value in tampered.items() if key != "plan_payload_sha256"}
    )
    with pytest.raises(ValueError, match="fixed contract"):
        conditioned.validate_eas_extension_plan(
            REPO_ROOT,
            tampered,
            base_manifest=manifest,
        )


def test_v2_quality_boundary_is_inclusive_and_redundant_weight_views_agree():
    boundary = conditioned.evaluate_eas_v2_bank_quality(
        {
            "global_effective_sample_size": 10_000,
            "global_weight_sum_over_max": 100.0,
            "global_max_normalized_weight": 0.01,
            "cap_hit_count": 0,
        }
    )
    assert boundary["quality_gates_passed"] is True
    assert boundary["failed_quality_gates"] == []

    ess_only = conditioned.evaluate_eas_v2_bank_quality(
        {
            "global_effective_sample_size": 9_999,
            "global_weight_sum_over_max": 100.0,
            "global_max_normalized_weight": 0.01,
            "cap_hit_count": 0,
        }
    )
    assert ess_only["failed_quality_gates"] == ["global_effective_sample_size"]

    cap_only = conditioned.evaluate_eas_v2_bank_quality(
        {
            "global_effective_sample_size": 10_000,
            "global_weight_sum_over_max": 100.0,
            "global_max_normalized_weight": 0.01,
            "cap_hit_count": 1,
        }
    )
    assert cap_only["failed_quality_gates"] == ["cap_hit_count"]

    below = conditioned.evaluate_eas_v2_bank_quality(
        {
            "global_effective_sample_size": 9_999,
            "global_weight_sum_over_max": 99.0,
            "global_max_normalized_weight": 1.0 / 99.0,
            "cap_hit_count": 1,
        }
    )
    assert set(below["failed_quality_gates"]) == {
        "global_effective_sample_size",
        "global_weight_sum_over_max",
        "global_max_normalized_weight",
        "cap_hit_count",
    }


def test_extended_loader_rejects_partial_fixed_inventory_before_loading(tmp_path):
    manifest = _synthetic_base_manifest()
    plan = conditioned.build_eas_extension_plan(REPO_ROOT, manifest)
    proposal_root = tmp_path / "proposals"
    proposal_root.mkdir()
    for index in range(conditioned.EAS_V2_TARGET_PROPOSAL_BLOCK_COUNT - 1):
        (proposal_root / f"block_{index:02d}").mkdir()
    with pytest.raises(ValueError, match="exactly fixed blocks"):
        conditioned.load_eas_extended_proposal_bank(
            REPO_ROOT,
            proposal_root,
            base_manifest=manifest,
            extension_plan=plan,
        )


def test_v2_blocks_reject_wrong_base_extension_provenance_and_base_regeneration(
    tmp_path,
):
    manifest = _synthetic_base_manifest()
    plan = conditioned.build_eas_extension_plan(REPO_ROOT, manifest)

    extension_dir = tmp_path / "block_20"
    _write_synthetic_production_block(
        extension_dir,
        block_index=20,
        manifest=manifest,
        plan=plan,
    )
    observed = conditioned.verify_eas_proposal_block(
        REPO_ROOT,
        extension_dir,
        expected_block_index=20,
        base_manifest=manifest,
        extension_plan=plan,
    )
    assert observed["source_binding_kind"] == "current_extension_implementation"
    assert (
        observed["base_manifest_payload_sha256"] == manifest["manifest_payload_sha256"]
    )
    assert observed["extension_plan_payload_sha256"] == plan["plan_payload_sha256"]

    completion_path = extension_dir / "completion.json"
    completion = json.loads(completion_path.read_text(encoding="utf-8"))
    completion["extension_bank_binding"]["extension_plan_payload_sha256"] = "0" * 64
    conditioned._atomic_json(completion_path, completion)
    with pytest.raises(ValueError, match="fixed v2 plan"):
        conditioned.verify_eas_proposal_block(
            REPO_ROOT,
            extension_dir,
            expected_block_index=20,
            base_manifest=manifest,
            extension_plan=plan,
        )

    fake_base_dir = tmp_path / "block_00_current_source"
    _write_synthetic_production_block(
        fake_base_dir,
        block_index=0,
        manifest=manifest,
        plan=plan,
    )
    with pytest.raises(ValueError, match="frozen EAS base block source"):
        conditioned.verify_eas_proposal_block(
            REPO_ROOT,
            fake_base_dir,
            expected_block_index=0,
            base_manifest=manifest,
            extension_plan=plan,
        )

    with pytest.raises(ValueError, match="cannot be regenerated"):
        conditioned.write_eas_proposal_block(
            REPO_ROOT,
            tmp_path / "missing_base",
            block_index=0,
            base_manifest=manifest,
            extension_plan=plan,
        )


def test_base_manifest_byte_validation_uses_stored_kernel_attestation_once(
    tmp_path, monkeypatch
):
    manifest = _synthetic_base_manifest()
    arrays = {
        "outcomes": np.asarray(["invalid_loss"]),
        "birth_generations": np.asarray([-1], dtype=np.int64),
        "log_forward": np.asarray([np.nan]),
        "log_reverse": np.asarray([np.nan]),
        "log_mutation": np.asarray([np.nan]),
        "log_weights": np.asarray([np.nan]),
    }
    attestation = manifest["kernel_equivalence_attestation"]
    attestation["array_sha256"] = {
        name: conditioned._proposal_array_sha256(value)
        for name, value in arrays.items()
    }
    attestation["attestation_payload_sha256"] = conditioned._canonical_sha256(
        {
            key: value
            for key, value in attestation.items()
            if key != "attestation_payload_sha256"
        }
    )
    manifest["manifest_payload_sha256"] = conditioned._canonical_sha256(
        {
            key: value
            for key, value in manifest.items()
            if key != "manifest_payload_sha256"
        }
    )
    observed_indexes: list[int] = []

    def fake_verify(*_args, **kwargs):
        observed_indexes.append(int(kwargs["expected_block_index"]))
        return arrays

    def forbidden_replay(*_args, **_kwargs):
        raise AssertionError("ordinary manifest validation reran the proposal kernel")

    monkeypatch.setattr(conditioned, "verify_eas_proposal_block", fake_verify)
    monkeypatch.setattr(conditioned, "simulate_eas_proposal_block", forbidden_replay)
    conditioned.validate_eas_base_bank_manifest(
        REPO_ROOT,
        manifest,
        proposal_dir=tmp_path,
    )
    assert observed_indexes == list(range(20))


def test_integer_eas_history_holds_ancient_ne_constant(eas_sizes):
    assert len(eas_sizes) == 500_002
    assert eas_sizes[0] == 37_637
    assert eas_sizes[40_000] == 13_138
    assert np.unique(eas_sizes[40_001:]).tolist() == [13_138]


def test_terminal_loss_weight_and_birth_cap_use_the_explicit_older_generation():
    sizes = np.asarray([10, 11, 13], dtype=np.int64)
    counts = np.asarray([3, 1], dtype=np.int64)
    observed = conditioned._reverse_path_log_components(sizes, counts)

    from scipy.stats import binom

    expected_forward = binom.logpmf(3, 20, 1 / 22)
    expected_reverse = binom.logpmf(1, 22, 3 / 20) + binom.logpmf(0, 26, 1 / 22)
    assert observed["log_forward"] == pytest.approx(expected_forward)
    assert observed["log_reverse"] == pytest.approx(expected_reverse)
    # Slatkin's mutation-opportunity factor belongs to the first zero/loss
    # generation (Ne=13), not to the one-copy birth row (Ne=11).
    assert observed["log_mutation"] == pytest.approx(np.log(13))
    assert observed["log_weight"] == pytest.approx(
        expected_forward - expected_reverse + np.log(13)
    )


def test_tiny_constant_ne_weighted_age_matches_forward_transition_oracle():
    from scipy.stats import binom

    ne = 6
    chromosomes = 2 * ne
    endpoint_count = 3
    birth_cap = 200
    transition = np.asarray(
        [
            [
                binom.pmf(younger, chromosomes, older / chromosomes)
                for younger in range(chromosomes + 1)
            ]
            for older in range(chromosomes + 1)
        ]
    )
    distribution = np.zeros(chromosomes + 1)
    distribution[1] = 1.0
    exact_age_mass = []
    for _age in range(birth_cap + 1):
        exact_age_mass.append(ne * distribution[endpoint_count])
        distribution = distribution @ transition
    exact_age_mass = np.asarray(exact_age_mass)
    exact_age_mass /= exact_age_mass.sum()
    exact_mean = float(np.dot(np.arange(birth_cap + 1), exact_age_mass))

    block = conditioned._simulate_reverse_proposal_block(
        np.full(birth_cap + 2, ne, dtype=np.int64),
        block_index=0,
        block_seed=12_345,
        proposals_per_block=20_000,
        present_count=endpoint_count,
        max_birth_generation=birth_cap,
    )
    valid = block["outcomes"] == "valid"
    log_weights = block["log_weights"][valid]
    weights = np.exp(log_weights - np.max(log_weights))
    estimated_mean = float(
        np.average(block["birth_generations"][valid], weights=weights)
    )
    assert np.count_nonzero(block["outcomes"] == "cap") == 0
    assert estimated_mean == pytest.approx(exact_mean, abs=0.5)


def test_small_proposal_block_replays_selected_path_exactly(
    eas_sizes, small_eas_blocks
):
    first = small_eas_blocks[0]
    valid = np.flatnonzero(first["outcomes"] == "valid")
    assert len(valid) > 0
    selected = int(valid[0])
    replay = conditioned.simulate_eas_proposal_block(
        eas_sizes,
        block_index=0,
        proposals_per_block=64,
        capture_indices=(selected,),
    )
    for key in (
        "outcomes",
        "birth_generations",
        "log_forward",
        "log_reverse",
        "log_mutation",
        "log_weights",
    ):
        np.testing.assert_array_equal(first[key], replay[key])
    path = replay["captured_paths"][selected]
    assert len(path) == int(first["birth_generations"][selected]) + 1
    assert path[0] == 14_704
    assert path[-1] == 1
    expected_weight = (
        first["log_forward"][selected]
        - first["log_reverse"][selected]
        + first["log_mutation"][selected]
    )
    assert first["log_weights"][selected] == pytest.approx(expected_weight)


def test_kernel_equivalence_attestation_exactly_replays_arrays(
    eas_sizes, small_eas_blocks
):
    frozen = small_eas_blocks[0]
    attestation = conditioned.attest_eas_proposal_kernel_equivalence(eas_sizes, frozen)
    assert attestation["status"] == "exact_match"
    assert attestation["attested_proposal_count"] == 64
    assert set(attestation["array_sha256"]) == {
        "outcomes",
        "birth_generations",
        "log_forward",
        "log_reverse",
        "log_mutation",
        "log_weights",
    }

    tampered = dict(frozen)
    for name in attestation["array_sha256"]:
        tampered[name] = np.asarray(frozen[name]).copy()
    valid_index = int(np.flatnonzero(tampered["outcomes"] == "valid")[0])
    tampered["log_forward"][valid_index] += 0.125
    tampered["log_weights"][valid_index] += 0.125
    with pytest.raises(RuntimeError, match="proposal kernel changed"):
        conditioned.attest_eas_proposal_kernel_equivalence(eas_sizes, tampered)


def test_global_selection_and_pass_b_with_small_blocks(eas_sizes, small_eas_blocks):
    selected, diagnostics = conditioned.select_eas_global_proposals(
        small_eas_blocks,
        accepted_count=1,
        enforce_production_contract=False,
    )
    assert len(selected) == 1
    assert diagnostics["proposal_blocks"] == 4
    assert diagnostics["total_proposals"] == 256
    paths = conditioned.replay_eas_global_selection(
        eas_sizes, small_eas_blocks, selected
    )
    path = paths[selected.iloc[0]["unit_id"]]
    assert path[0] == 14_704
    assert path[-1] == 1


def test_global_categorical_selection_retains_duplicate_occurrences(
    small_eas_blocks,
):
    selected, diagnostics = conditioned.select_eas_global_proposals(
        small_eas_blocks,
        accepted_count=500,
        enforce_production_contract=False,
    )
    assert len(selected) == 500
    assert diagnostics["resampling_with_replacement"] is True
    assert diagnostics["selected_duplicate_occurrences"] > 0
    assert selected["unit_id"].is_unique
    duplicated = selected["global_proposal_index"].duplicated(keep=False)
    assert duplicated.any()
    assert (selected.loc[duplicated, "source_particle_draw_multiplicity"] > 1).all()


def test_restartable_eas_block_and_selected_path_bundles_are_hash_bound(
    tmp_path, eas_sizes
):
    block_dir = tmp_path / "block_000"
    block = conditioned.write_eas_proposal_block(
        REPO_ROOT,
        block_dir,
        block_index=0,
        proposals_per_block=64,
        enforce_production_contract=False,
    )
    assert len(block["proposal_block_sha256"]) == 64
    selected, diagnostics = conditioned.select_eas_global_proposals(
        [block],
        accepted_count=1,
        enforce_production_contract=False,
    )
    paths = conditioned.replay_eas_global_selection(eas_sizes, [block], selected)
    row = selected.iloc[0].to_dict()
    bundle_dir = tmp_path / row["unit_id"]
    completion = conditioned.write_eas_selected_path_bundle(
        REPO_ROOT,
        bundle_dir,
        selected_row=row,
        derived_counts_ascending=paths[row["unit_id"]],
        global_diagnostics=diagnostics,
    )
    assert completion["status"] == "complete"
    assert Path(completion["trajectory_path"]).is_file()

    unexpected = block_dir / "unexpected.txt"
    unexpected.write_text("tamper\n", encoding="utf-8")
    with pytest.raises(ValueError, match="directory inventory"):
        conditioned.verify_eas_proposal_block(
            REPO_ROOT,
            block_dir,
            expected_block_index=0,
            expected_proposals=64,
            enforce_production_contract=False,
        )


def test_replay_configs_are_s0_and_keep_population_not_sample_conditioning():
    eas = conditioned.render_eas_replay_parameter_text(
        REPO_ROOT,
        birth_generation=10_000,
    )
    eas_audit = conditioned.lint_conditioned_parameter_text(
        eas,
        model_id=conditioned.EAS_MODEL_ID,
        birth_generation=10_000,
        screen=False,
    )
    assert eas_audit["selection_coefficient"] == 0
    assert eas_audit["endpoint_token"] == (
        f"{conditioned.EAS_SELECTED_PRESENT_AF:.17g}"
    )
    assert eas_audit["sample_af_conditioned"] is False

    screen = conditioned.render_han_parameter_text(REPO_ROOT, screen=True)
    full = conditioned.render_han_parameter_text(REPO_ROOT, screen=False)
    assert "length 2" in screen
    assert "sample_size 3 2" in screen
    assert "length 10000000" in full
    assert "sample_size 3 200" in full
    for text, is_screen in ((screen, True), (full, False)):
        audit = conditioned.lint_conditioned_parameter_text(
            text,
            model_id=conditioned.HAN_MODEL_ID,
            birth_generation=conditioned.HAN_BIRTH_GENERATION,
            screen=is_screen,
        )
        assert audit["selection_coefficient"] == 0
        assert audit["endpoint_token"] == conditioned.HAN_PRESENT_TOKEN
        assert audit["sample_af_conditioned"] is False


def _han_trajectory(*, boundary_af: float, endpoint_af: float) -> pd.DataFrame:
    generations = np.arange(conditioned.HAN_BIRTH_GENERATION, -1, -1)
    rows = len(generations)
    data: dict[str, np.ndarray] = {
        "sim": np.ones(rows, dtype=int),
        "gen": generations,
    }
    for population in range(1, 6):
        data[f"selfreq_{population}"] = np.zeros(rows, dtype=float)
        data[f"popsize_{population}"] = np.full(rows, 10_000, dtype=int)
    data["selfreq_4"][0] = 1 / 20_000
    boundary_index = int(
        conditioned.HAN_BIRTH_GENERATION - conditioned.HAN_GATE_GENERATION
    )
    data["selfreq_3"][boundary_index:] = np.linspace(
        boundary_af, endpoint_af, rows - boundary_index
    )
    return pd.DataFrame(data)


def test_han_trajectory_requires_g645_gate_and_strict_present_segregation(tmp_path):
    valid_path = tmp_path / "valid.tsv"
    _han_trajectory(boundary_af=0.012, endpoint_af=0.02).to_csv(
        valid_path, sep="\t", index=False
    )
    audit = conditioned.validate_han_trajectory(valid_path)
    assert audit["han_generation_645_gate_passed"] is True
    assert audit["present_strictly_segregating"] is True

    outside = tmp_path / "outside.tsv"
    _han_trajectory(boundary_af=0.02, endpoint_af=0.03).to_csv(
        outside, sep="\t", index=False
    )
    with pytest.raises(ValueError, match="generation-645"):
        conditioned.validate_han_trajectory(outside)

    fixed = tmp_path / "fixed.tsv"
    _han_trajectory(boundary_af=0.012, endpoint_af=1.0).to_csv(
        fixed, sep="\t", index=False
    )
    with pytest.raises(ValueError, match="strictly segregating"):
        conditioned.validate_han_trajectory(fixed)


def test_han_first_100_pass_rule_is_candidate_order_not_input_order():
    ledger = conditioned.build_han_candidate_ledger().iloc[:110].copy()
    results = pd.DataFrame(
        {
            "candidate_id": ledger["candidate_id"].sample(frac=1, random_state=17),
        }
    )
    results["trajectory_sha256"] = "a" * 64
    results["present_chb_af"] = 0.02
    results["han_generation_645_chb_af"] = 0.012
    results["gate_passed"] = True
    results["strictly_segregating"] = True
    results["screen_completion_sha256"] = "b" * 64
    accepted = conditioned.select_han_acceptances(ledger, results)
    assert accepted["candidate_order"].tolist() == list(range(100))
    assert accepted["accepted_index"].tolist() == list(range(100))


def test_han_acceptance_boolean_strings_are_parsed_strictly():
    ledger = conditioned.build_han_candidate_ledger().iloc[:101].copy()
    results = pd.DataFrame({"candidate_id": ledger["candidate_id"]})
    results["trajectory_sha256"] = "a" * 64
    results["present_chb_af"] = 0.02
    results["han_generation_645_chb_af"] = 0.012
    results["gate_passed"] = "true"
    results["strictly_segregating"] = "true"
    results.loc[0, "gate_passed"] = "False"
    results["screen_completion_sha256"] = "b" * 64
    accepted = conditioned.select_han_acceptances(ledger, results)
    assert 0 not in set(accepted["candidate_order"])
    assert accepted["candidate_order"].tolist() == list(range(1, 101))

    results.loc[0, "gate_passed"] = "not-a-bool"
    with pytest.raises(ValueError, match="non-boolean"):
        conditioned.select_han_acceptances(ledger, results)
