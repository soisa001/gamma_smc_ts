from __future__ import annotations

import hashlib
import json
import os
import queue
import socket
import sys
from dataclasses import asdict
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
import tskit

import gamma_smc_aou.eas_sweep_study as study
from gamma_smc_aou.eas_sweep_models import (
    INTROGRESSION_TARGET_FREQUENCIES,
    PHLASH_MVN_SCHEMA,
    SELECTION_COEFFICIENTS,
)


def _write_tiny_eas_resource(repo_root):
    resource_dir = repo_root / study.NO_INTROGRESSION_DIR / "resources"
    resource_dir.mkdir(parents=True)
    time = np.array([100.0, 1_000.0, 10_000.0, 40_000.0])
    bootstrap_ne = np.array(
        [
            [8_000.0, 9_000.0, 10_000.0, 11_000.0],
            [10_000.0, 11_000.0, 12_000.0, 13_000.0],
            [12_000.0, 13_000.0, 14_000.0, 15_000.0],
            [14_000.0, 15_000.0, 16_000.0, 17_000.0],
        ]
    )
    path = resource_dir / "EAS.npz"
    np.savez_compressed(
        path,
        schema=np.asarray(PHLASH_MVN_SCHEMA),
        population=np.asarray("EAS"),
        time=time,
        mean_log_ne=np.log(np.median(bootstrap_ne, axis=0)),
        covariance_factor=np.full((2, time.size), 0.01),
        bootstrap_ne=bootstrap_ne,
        jitter=np.asarray(0.001),
    )
    return path, hashlib.sha256(path.read_bytes()).hexdigest()


def _accepted_tiny_tree_sequence(
    selected_sample_indices=(4, 6, 7, 8, 9),
):
    """Five diploids at 50% AF: 2 ref/ref, 1 ref/alt, and 2 alt/alt."""
    tables = tskit.TableCollection(sequence_length=10_000_000)
    tables.populations.metadata_schema = tskit.MetadataSchema.permissive_json()
    population = tables.populations.add_row(metadata={"name": "EAS"})
    individuals = [tables.individuals.add_row() for _ in range(5)]
    samples = [
        tables.nodes.add_row(
            flags=tskit.NODE_IS_SAMPLE,
            time=0,
            population=population,
            individual=individuals[index // 2],
        )
        for index in range(10)
    ]
    selected_parent = tables.nodes.add_row(time=1, population=population)
    root = tables.nodes.add_row(time=2, population=population)
    # Selected haplotypes yield per-diploid counts [0, 0, 1, 2, 2].
    selected_samples = tuple(samples[index] for index in selected_sample_indices)
    for child in selected_samples:
        tables.edges.add_row(0, 10_000_000, selected_parent, child)
    tables.edges.add_row(0, 10_000_000, root, selected_parent)
    for child in samples:
        if child not in selected_samples:
            tables.edges.add_row(0, 10_000_000, root, child)
    site = tables.sites.add_row(position=5_000_000, ancestral_state="A")
    tables.mutations.add_row(site=site, node=selected_parent, derived_state="G")
    tables.sort()
    return tables.tree_sequence()


def _runtime(repo_root, **overrides):
    values = {
        "repo_root": str(repo_root),
        "sample_diploids": 5,
        "simulation_workers": 1,
        "max_external_draws": 1,
    }
    values.update(overrides)
    return study.StudyRuntime(**values)


def _accepted_specification():
    return {
        "scenario": "no_introgression_de_novo",
        "study_type": "no_introgression",
        "specification_id": "tiny_no_intro_s0p010_median",
        "selection_coefficient": 0.01,
        "demography_curve": "median",
        "target_allele_frequency": 0.5,
        "lower_allele_frequency": 0.45,
        "upper_allele_frequency": 0.55,
        "origin_age_generations": 1_000.0,
        "theta_for_gamma_smc": 5e-4,
    }


def test_plan_grids_cover_nine_de_novo_and_twenty_seven_introgression_specs(
    tmp_path, monkeypatch
):
    _, digest = _write_tiny_eas_resource(tmp_path)
    monkeypatch.setattr(study, "EAS_RESOURCE_SHA256", digest)
    runtime = _runtime(tmp_path)

    specifications = study.build_specifications(
        runtime, ("no_introgression", "introgression")
    )

    no_intro = specifications["no_introgression"]
    intro = specifications["introgression"]
    assert len(no_intro) == 9
    assert len(intro) == 27
    assert len({row["specification_id"] for row in no_intro + intro}) == 36
    assert {
        (row["selection_coefficient"], row["demography_curve"]) for row in no_intro
    } == {
        (selection, curve)
        for selection in SELECTION_COEFFICIENTS
        for curve in ("q025", "median", "q975")
    }
    assert {
        (row["selection_coefficient"], row["target_allele_frequency"]) for row in intro
    } == {
        (selection, target)
        for selection in SELECTION_COEFFICIENTS
        for target in INTROGRESSION_TARGET_FREQUENCIES
    }
    assert all(row["target_allele_frequency"] == 0.5 for row in no_intro)
    assert all(
        row["deterministic_initial_frequency"]
        == pytest.approx(runtime.slim_scaling_factor / (2 * row["ne_at_origin"]))
        for row in no_intro
    )
    assert all(
        row["deterministic_initial_frequency_slim_scaling_factor"]
        == runtime.slim_scaling_factor
        for row in no_intro
    )
    assert all(row["upper_allele_frequency"] == 0.95 for row in intro)
    assert all(
        row["slim_scaled_post_origin_survival_check_generations"]
        == row["slim_scaled_origin_age_generations"] - 10
        for row in no_intro
    )
    assert all(
        row["slim_scaled_source_check_generations_ago"]
        == row["slim_scaled_pulse_generations_ago"]
        == 2_270
        for row in intro
    )
    assert all(
        row["requested_han_split_generations_ago"] == 2_016
        and row["realized_han_split_generations_ago"] == 2_010
        and row["slim_demographic_han_split_generations_ago"] == 2_010
        and row["slim_selection_transition_generations_ago"] == 2_010
        for row in intro
    )
    assert study._scaled_demographic_time(2_016, 1) == 2_016  # noqa: SLF001


def test_introgression_feasibility_is_screening_only_and_seeds_are_stable(tmp_path):
    runtime = _runtime(tmp_path)
    intro = study.build_specifications(runtime, ("introgression",))["introgression"]

    feasible = {
        (row["selection_coefficient"], row["target_allele_frequency"]): row[
            "deterministic_feasible_from_fixed_archaic_source"
        ]
        for row in intro
    }
    assert feasible[(0.01, 0.9)] is True
    assert feasible[(0.005, 0.8)] is True
    assert feasible[(0.005, 0.9)] is False
    assert feasible[(0.001, 0.1)] is False
    assert all(row["eligible_for_simulation"] for row in intro)
    assert all(
        row["eligibility_reason"] == "conditioned_stochastic_tail"
        for row in intro
        if not row["deterministic_feasible_from_fixed_archaic_source"]
    )

    assert (
        study._stable_seed(  # noqa: SLF001
            study.BASE_SEED, "no_intro_s0p010_median_age2000g", 0
        )
        == 495_446_250
    )
    assert (
        study._stable_seed(  # noqa: SLF001
            study.BASE_SEED, "no_intro_s0p010_median_age2000g", 1
        )
        == 275_164_609
    )
    assert (
        study._stable_seed(  # noqa: SLF001
            study.BASE_SEED, "intro_s0p005_af50", 0
        )
        == 804_816_346
    )


def test_write_study_plans_records_resource_identity_and_design(tmp_path, monkeypatch):
    resource, digest = _write_tiny_eas_resource(tmp_path)
    monkeypatch.setattr(study, "EAS_RESOURCE_SHA256", digest)
    runtime = _runtime(tmp_path)

    specifications = study.write_study_plans(
        runtime, ("no_introgression", "introgression")
    )

    for scenario, expected_count in (("no_introgression", 9), ("introgression", 27)):
        root = study._study_dir(tmp_path, scenario)  # noqa: SLF001
        design = json.loads((root / "study_design.json").read_text(encoding="utf-8"))
        table = pd.read_csv(root / "specifications.tsv", sep="\t")
        assert design["schema"] == study.SCHEMA_VERSION
        assert design["specification_count"] == expected_count
        assert design["rng"]["base_seed"] == study.BASE_SEED
        assert design["tmrca_thresholds_years"] == [
            1_000,
            4_500,
            10_000,
            20_000,
            30_000,
            40_000,
            50_000,
        ]
        assert len(table) == len(specifications[scenario]) == expected_count

    contracts = json.loads(
        (tmp_path / study.NO_INTROGRESSION_DIR / "origin_age_contracts.json").read_text(
            encoding="utf-8"
        )
    )
    assert len(contracts) == 9
    for contract in contracts.values():
        source = contract["demography"]["source"]
        assert source["source_sha256_expected"] == digest
        assert source["source_sha256_actual"] == digest
        assert source["source_sha256_verified"] is True
        assert source["source_path"] == "resources/EAS.npz"

    intro_design = json.loads(
        (tmp_path / study.INTROGRESSION_DIR / "study_design.json").read_text(
            encoding="utf-8"
        )
    )
    scaled_times = intro_design["scaled_introgression_event_times_generations_ago"]
    assert scaled_times["requested_catalog_han_split"] == 2_016
    assert scaled_times["realized_demographic_han_split"] == 2_010
    assert scaled_times["selection_population_transition"] == 2_010
    assert intro_design["scaled_time_rules"]["han_split_alignment_validated"] is True


def test_accepted_tree_validation_checks_sample_af_and_all_genotype_classes(tmp_path):
    runtime = _runtime(tmp_path)
    specification = _accepted_specification()

    pair_table, validation = study._validate_accepted_tree(  # noqa: SLF001
        _accepted_tiny_tree_sequence(), specification, runtime
    )

    assert pair_table["genotype_class"].tolist() == [
        "hom_ref",
        "hom_ref",
        "heterozygous",
        "hom_alt",
        "hom_alt",
    ]
    assert validation == {
        "accepted": True,
        "achieved_allele_frequency": 0.5,
        "achieved_allele_frequency_scope": "sample_estimate",
        "population_frequency_condition_lower": 0.45,
        "population_frequency_condition_upper": 0.55,
        "n_sample_diploids": 5,
        "n_hom_ref": 2,
        "n_heterozygous": 1,
        "n_hom_alt": 2,
        "focal_site_id": 0,
        "focal_mutation_count": 1,
        "focal_mutation_node_population": "EAS",
        "focal_mutation_declared_origin_population": "EAS",
        "focal_mutation_node_population_is_origin_evidence": False,
        "focal_mutation_origin_evidence": "serialized_forward_event_contract",
        "n_segregating_nucleotide_snps": 1,
        "vcf_allele_encoding": "single_nucleotide_ACGT",
    }

    too_narrow = {**specification, "lower_allele_frequency": 0.51}
    _, rejected = study._validate_accepted_tree(  # noqa: SLF001
        _accepted_tiny_tree_sequence(), too_narrow, runtime
    )
    assert rejected["accepted"] is False


def test_simplified_mutation_node_population_is_not_used_as_origin_evidence(tmp_path):
    runtime = _runtime(tmp_path)
    specification = {
        **_accepted_specification(),
        "study_type": "introgression",
        "lower_allele_frequency": 0.10,
        "upper_allele_frequency": 0.95,
    }

    _, validation = study._validate_accepted_tree(  # noqa: SLF001
        _accepted_tiny_tree_sequence(), specification, runtime
    )

    assert validation["accepted"] is True
    assert validation["focal_mutation_node_population"] == "EAS"
    assert validation["focal_mutation_declared_origin_population"] == "Neanderthal"
    assert validation["focal_mutation_node_population_is_origin_evidence"] is False


def _install_fake_simulator(monkeypatch, tree_sequence):
    calls = []

    class FakeEngine:
        def simulate(self, *args, **kwargs):
            calls.append((args, kwargs))
            return tree_sequence

    sweep = SimpleNamespace(contig=object(), extended_events=())
    monkeypatch.setattr(
        study,
        "_software_record",
        lambda **kwargs: {
            "python": "test",
            "stdpopsim": "0.3.0",
            "tskit": tskit.__version__,
            "slim_path": str(kwargs.get("slim_path")),
            "slim_sha256": "test-double",
            "slim_version_output": "SLiM 4.2.2 test double",
        },
    )
    monkeypatch.setattr(
        study,
        "_build_simulation_objects",
        lambda specification, runtime: (
            object(),
            sweep,
            {"EAS": runtime.sample_diploids},
            {"schema": "tiny-selected-site/v1", "population": "EAS"},
        ),
    )
    monkeypatch.setattr(study.stdpopsim, "get_engine", lambda name: FakeEngine())
    monkeypatch.setattr(
        study,
        "tree_truth_profiles",
        lambda *args, **kwargs: pd.DataFrame(
            {
                "position_0based": [5_000_000],
                "genotype_class": ["overall"],
                "threshold_years": [1_000],
                "p_tmrca_lt_threshold": [0.5],
            }
        ),
    )
    monkeypatch.setattr(
        study,
        "internal_genotype_contrasts",
        lambda frame: pd.DataFrame({"position_0based": [5_000_000], "contrast": [0.0]}),
    )
    return calls


def test_completed_simulation_is_reused_without_another_slim_draw(
    tmp_path, monkeypatch
):
    runtime = _runtime(tmp_path)
    specification = _accepted_specification()
    calls = _install_fake_simulator(monkeypatch, _accepted_tiny_tree_sequence())
    task = {
        "specification": specification,
        "runtime": asdict(runtime),
        "slim_path": sys.executable,
        "specification_dir": str(tmp_path / "work" / specification["specification_id"]),
    }

    first = study._simulate_specification_task(task)  # noqa: SLF001
    second = study._simulate_specification_task(task)  # noqa: SLF001

    assert first["status"] == "complete"
    assert second["status"] == "cached"
    assert len(calls) == 1


def test_accepted_tree_checkpoint_resumes_materialization_without_new_draw(
    tmp_path, monkeypatch
):
    runtime = _runtime(tmp_path)
    specification = _accepted_specification()
    calls = _install_fake_simulator(monkeypatch, _accepted_tiny_tree_sequence())
    specification_dir = tmp_path / "work" / specification["specification_id"]
    task = {
        "specification": specification,
        "runtime": asdict(runtime),
        "slim_path": sys.executable,
        "specification_dir": str(specification_dir),
    }
    atomic_frame = study._atomic_frame  # noqa: SLF001
    interrupted = False

    def interrupt_after_checkpoint(path, frame):
        nonlocal interrupted
        if (
            not interrupted
            and str(path).endswith("external_draws.tsv")
            and (specification_dir / "accepted_checkpoint.trees").is_file()
            and not frame.empty
        ):
            interrupted = True
            raise RuntimeError("materialization interruption")
        return atomic_frame(path, frame)

    monkeypatch.setattr(study, "_atomic_frame", interrupt_after_checkpoint)
    with pytest.raises(RuntimeError, match="materialization interruption"):
        study._simulate_specification_task(task)  # noqa: SLF001
    monkeypatch.setattr(study, "_atomic_frame", atomic_frame)

    result = study._simulate_specification_task(task)  # noqa: SLF001

    assert result["status"] == "complete"
    assert len(calls) == 1
    assert not (specification_dir / "accepted_checkpoint.trees").exists()
    attempts = pd.read_csv(specification_dir / "external_draws.tsv", sep="\t")
    assert attempts["external_draw_zero_based"].tolist() == [0]


def test_accepted_checkpoint_reconciles_crash_after_attempt_ledger(
    tmp_path, monkeypatch
):
    runtime = _runtime(tmp_path)
    specification = _accepted_specification()
    calls = _install_fake_simulator(monkeypatch, _accepted_tiny_tree_sequence())
    specification_dir = tmp_path / "work" / specification["specification_id"]
    task = {
        "specification": specification,
        "runtime": asdict(runtime),
        "slim_path": sys.executable,
        "specification_dir": str(specification_dir),
    }
    atomic_json = study._atomic_json  # noqa: SLF001
    interrupted = False

    def interrupt_state_after_attempt_ledger(path, payload):
        nonlocal interrupted
        if (
            not interrupted
            and str(path).endswith("external_draw_state.json")
            and (specification_dir / "external_draws.tsv").is_file()
            and payload.get("completed_external_sample_panel_draws") == 1
            and payload.get("active_internal_conditioned_trajectory") is None
        ):
            interrupted = True
            raise RuntimeError("crash after accepted attempt ledger")
        return atomic_json(path, payload)

    monkeypatch.setattr(study, "_atomic_json", interrupt_state_after_attempt_ledger)
    with pytest.raises(RuntimeError, match="crash after accepted attempt ledger"):
        study._simulate_specification_task(task)  # noqa: SLF001
    monkeypatch.setattr(study, "_atomic_json", atomic_json)

    result = study._simulate_specification_task(task)  # noqa: SLF001

    assert result["status"] == "complete"
    assert len(calls) == 1
    state = json.loads(
        (specification_dir / "external_draw_state.json").read_text(encoding="utf-8")
    )
    assert state["active_internal_conditioned_trajectory"] is None
    assert state["completed_external_sample_panel_draws"] == 1
    assert state["interrupted_internal_conditioned_trajectories"] == []


def test_launched_seed_cap_counts_interrupted_internal_trajectories(
    tmp_path, monkeypatch
):
    runtime = _runtime(
        tmp_path,
        max_external_draws=20,
        max_launched_draws=2,
    )
    specification = _accepted_specification()
    _install_fake_simulator(monkeypatch, _accepted_tiny_tree_sequence())
    seeds = []

    class InterruptedEngine:
        def simulate(self, *args, **kwargs):
            seeds.append(kwargs["seed"])
            raise RuntimeError("interrupted")

    monkeypatch.setattr(study.stdpopsim, "get_engine", lambda name: InterruptedEngine())
    specification_dir = tmp_path / "work" / specification["specification_id"]
    task = {
        "specification": specification,
        "runtime": asdict(runtime),
        "slim_path": sys.executable,
        "specification_dir": str(specification_dir),
    }

    for _ in range(2):
        with pytest.raises(RuntimeError, match="interrupted"):
            study._simulate_specification_task(task)  # noqa: SLF001
    with pytest.raises(study.SimulationBudgetExhausted, match="2 launched seeds"):
        study._simulate_specification_task(task)  # noqa: SLF001

    assert seeds == [
        study._stable_seed(  # noqa: SLF001
            runtime.base_seed, specification["specification_id"], draw
        )
        for draw in range(2)
    ]


def test_draw_budget_is_reported_as_nonretryable_exhaustion(tmp_path, monkeypatch):
    runtime = _runtime(tmp_path, max_external_draws=1)
    specification = _accepted_specification()
    _install_fake_simulator(monkeypatch, _accepted_tiny_tree_sequence((4, 6, 8, 9)))

    result = study.simulate_studies(
        runtime,
        {"no_introgression": [specification]},
        slim_path=sys.executable,
    )[specification["specification_id"]]

    assert result["status"] == "exhausted"
    assert result["failure"]["error_type"] == "SimulationBudgetExhausted"
    assert result["failure"]["retryable"] is False


def test_completion_and_history_record_effective_timeout(tmp_path, monkeypatch):
    runtime = _runtime(tmp_path, max_specification_seconds=300.0)
    specification = _accepted_specification()
    _install_fake_simulator(monkeypatch, _accepted_tiny_tree_sequence())
    specification_dir = tmp_path / "work" / specification["specification_id"]
    task = {
        "specification": specification,
        "runtime": asdict(runtime),
        "slim_path": sys.executable,
        "specification_dir": str(specification_dir),
        "effective_timeout_seconds": 75.0,
    }

    result = study._simulate_specification_task(task)  # noqa: SLF001
    history = pd.read_csv(specification_dir / "simulation_invocations.tsv", sep="\t")

    assert result["completion"]["declared_timeout_seconds"] == 300.0
    assert result["completion"]["effective_timeout_seconds"] == 75.0
    assert history["declared_timeout_seconds"].tolist() == [300.0]
    assert history["effective_timeout_seconds"].tolist() == [75.0]
    assert history["simulation_contract_sha256"].str.fullmatch(r"[0-9a-f]{64}").all()
    assert history["core_source_sha256"].str.fullmatch(r"[0-9a-f]{64}").all()
    assert history["orchestrator_source_sha256"].str.fullmatch(r"[0-9a-f]{64}").all()
    assert (
        result["completion"]["implementation"]["sources"]["core"]["sha256"]
        == history.loc[0, "core_source_sha256"]
    )


def test_simulation_cache_requires_all_recorded_outputs(tmp_path, monkeypatch):
    runtime = _runtime(tmp_path)
    specification = _accepted_specification()
    _install_fake_simulator(monkeypatch, _accepted_tiny_tree_sequence())
    specification_dir = tmp_path / "work" / specification["specification_id"]
    task = {
        "specification": specification,
        "runtime": asdict(runtime),
        "slim_path": sys.executable,
        "specification_dir": str(specification_dir),
    }
    study._simulate_specification_task(task)  # noqa: SLF001
    contract = json.loads(
        (specification_dir / "simulation_contract.json").read_text(encoding="utf-8")
    )

    assert (
        study._valid_simulation_cache(  # noqa: SLF001
            specification_dir, contract, specification, runtime
        )
        is not None
    )
    (specification_dir / "truth_profiles.tsv.gz").unlink()
    assert (
        study._valid_simulation_cache(  # noqa: SLF001
            specification_dir, contract, specification, runtime
        )
        is None
    )


def test_simulation_cache_accepts_legacy_windows_pair_separators(tmp_path, monkeypatch):
    runtime = _runtime(tmp_path)
    specification = _accepted_specification()
    _install_fake_simulator(monkeypatch, _accepted_tiny_tree_sequence())
    specification_dir = tmp_path / "work" / specification["specification_id"]
    task = {
        "specification": specification,
        "runtime": asdict(runtime),
        "slim_path": sys.executable,
        "specification_dir": str(specification_dir),
    }
    study._simulate_specification_task(task)  # noqa: SLF001
    contract = json.loads(
        (specification_dir / "simulation_contract.json").read_text(encoding="utf-8")
    )
    completion_path = specification_dir / "simulation_complete.json"
    completion = json.loads(completion_path.read_text(encoding="utf-8"))
    for pair_record in completion["pair_files"].values():
        pair_record["path"] = pair_record["path"].replace("/", "\\")
    completion_path.write_text(json.dumps(completion), encoding="utf-8")

    assert (
        study._valid_simulation_cache(  # noqa: SLF001
            specification_dir, contract, specification, runtime
        )
        is not None
    )


def test_outer_seed_progress_distinguishes_interrupted_internal_trajectory(
    tmp_path, monkeypatch
):
    runtime = _runtime(tmp_path, max_external_draws=2)
    specification = _accepted_specification()
    _install_fake_simulator(monkeypatch, _accepted_tiny_tree_sequence())
    outcomes = [
        _accepted_tiny_tree_sequence((4, 6, 8, 9)),
        RuntimeError("interrupted internal conditioned trajectory"),
        _accepted_tiny_tree_sequence(),
    ]
    seeds = []

    class SequencedEngine:
        def simulate(self, *args, **kwargs):
            seeds.append(kwargs["seed"])
            outcome = outcomes.pop(0)
            if isinstance(outcome, Exception):
                raise outcome
            return outcome

    engine = SequencedEngine()
    monkeypatch.setattr(study.stdpopsim, "get_engine", lambda name: engine)
    specification_dir = tmp_path / "work" / specification["specification_id"]
    task = {
        "specification": specification,
        "runtime": asdict(runtime),
        "slim_path": sys.executable,
        "specification_dir": str(specification_dir),
    }

    with pytest.raises(RuntimeError, match="interrupted internal"):
        study._simulate_specification_task(task)  # noqa: SLF001
    result = study._simulate_specification_task(task)  # noqa: SLF001

    expected_seeds = [
        study._stable_seed(  # noqa: SLF001
            runtime.base_seed, specification["specification_id"], draw
        )
        for draw in range(3)
    ]
    assert seeds == expected_seeds
    attempts = pd.read_csv(specification_dir / "external_draws.tsv", sep="\t")
    assert attempts["external_draw_zero_based"].tolist() == [0, 2]
    assert result["completion"]["accepted_external_draw_zero_based"] == 2
    assert result["completion"]["attempts_to_accept"] == 2
    state = json.loads(
        (specification_dir / "external_draw_state.json").read_text(encoding="utf-8")
    )
    assert state["completed_external_sample_panel_draws"] == 2
    assert len(state["interrupted_internal_conditioned_trajectories"]) == 1
    interrupted = state["interrupted_internal_conditioned_trajectories"][0]
    assert interrupted["external_draw_zero_based"] == 1
    assert interrupted["seed"] == expected_seeds[1]
    assert interrupted["outcome"] == (
        "interrupted_stdpopsim_internal_conditioned_trajectory"
    )
    assert interrupted["external_sample_panel_completed"] is False


def test_legacy_draw_logs_advance_seed_only_under_unchanged_contract(
    tmp_path, monkeypatch
):
    runtime = _runtime(tmp_path, max_external_draws=2)
    specification = _accepted_specification()
    _install_fake_simulator(monkeypatch, _accepted_tiny_tree_sequence())
    outcomes = [
        _accepted_tiny_tree_sequence((4, 6, 8, 9)),
        RuntimeError("legacy timeout"),
        _accepted_tiny_tree_sequence(),
    ]
    seeds = []

    class SequencedEngine:
        def simulate(self, *args, **kwargs):
            seeds.append(kwargs["seed"])
            outcome = outcomes.pop(0)
            if isinstance(outcome, Exception):
                raise outcome
            return outcome

    monkeypatch.setattr(study.stdpopsim, "get_engine", lambda name: SequencedEngine())
    specification_dir = tmp_path / "work" / specification["specification_id"]
    task = {
        "specification": specification,
        "runtime": asdict(runtime),
        "slim_path": sys.executable,
        "specification_dir": str(specification_dir),
    }
    with pytest.raises(RuntimeError, match="legacy timeout"):
        study._simulate_specification_task(task)  # noqa: SLF001
    (specification_dir / "external_draw_state.json").unlink()
    (specification_dir / "slim_draw_000.csv").write_text("legacy\n", encoding="utf-8")
    (specification_dir / "slim_draw_001.csv").write_text("legacy\n", encoding="utf-8")

    result = study._simulate_specification_task(task)  # noqa: SLF001

    assert result["completion"]["accepted_external_draw_zero_based"] == 2
    assert seeds == [
        study._stable_seed(  # noqa: SLF001
            runtime.base_seed, specification["specification_id"], draw
        )
        for draw in range(3)
    ]
    state = json.loads(
        (specification_dir / "external_draw_state.json").read_text(encoding="utf-8")
    )
    assert state["interrupted_internal_conditioned_trajectories"] == [
        {
            "external_draw_zero_based": 1,
            "seed": seeds[1],
            "stage": "legacy_uncheckpointed_stdpopsim_conditioning",
            "outcome": (
                "legacy_uncheckpointed_stdpopsim_internal_conditioned_trajectory"
            ),
            "external_sample_panel_completed": False,
        }
    ]


def test_invocation_history_filters_contracts_and_fails_closed_on_corruption(
    tmp_path,
):
    specification_dir = tmp_path / "work" / "spec"
    specification_dir.mkdir(parents=True)
    current_sha = "a" * 64
    stale_sha = "b" * 64
    history_path = specification_dir / "simulation_invocations.tsv"
    frame = pd.DataFrame(
        [
            {
                "invocation_index": 1,
                "recorded_at_utc": "2026-08-12T00:00:00+00:00",
                "simulation_contract_sha256": stale_sha,
                "status": "timed_out",
                "elapsed_seconds": 300.0,
                "declared_timeout_seconds": 300.0,
                "effective_timeout_seconds": 300.0,
                "completed_external_sample_panel_draws": 1,
                "interrupted_internal_conditioned_trajectories": 1,
                "error_type": "SpecificationTimeout",
            },
            {
                "invocation_index": 2,
                "recorded_at_utc": "2026-08-12T00:05:00+00:00",
                "simulation_contract_sha256": current_sha,
                "status": "complete",
                "elapsed_seconds": 10.0,
                "declared_timeout_seconds": 300.0,
                "effective_timeout_seconds": 25.0,
                "completed_external_sample_panel_draws": 2,
                "interrupted_internal_conditioned_trajectories": 1,
                "error_type": None,
            },
        ]
    )
    frame.to_csv(history_path, sep="\t", index=False)

    current = study._simulation_invocation_history(  # noqa: SLF001
        specification_dir, expected_contract_sha256=current_sha
    )

    assert current["invocation_index"].tolist() == [2]
    assert current["elapsed_seconds"].tolist() == [10.0]
    frame["elapsed_seconds"] = frame["elapsed_seconds"].astype(object)
    frame.loc[1, "elapsed_seconds"] = "corrupt"
    frame.to_csv(history_path, sep="\t", index=False)
    with pytest.raises(ValueError, match="invalid elapsed_seconds"):
        study._simulation_invocation_history(  # noqa: SLF001
            specification_dir, expected_contract_sha256=current_sha
        )


def test_valid_legacy_invocation_history_is_migrated_to_current_contract(tmp_path):
    specification_dir = tmp_path / "work" / "spec"
    specification_dir.mkdir(parents=True)
    history_path = specification_dir / "simulation_invocations.tsv"
    legacy = pd.DataFrame(
        [
            {
                "invocation_index": 1,
                "recorded_at_utc": "2026-08-12T00:00:00+00:00",
                "status": "timed_out",
                "elapsed_seconds": 12.5,
                "declared_timeout_seconds": 20.0,
                "completed_external_sample_panel_draws": 1,
                "interrupted_internal_conditioned_trajectories": 1,
                "error_type": "SpecificationTimeout",
            }
        ]
    )
    legacy.to_csv(history_path, sep="\t", index=False)
    contract_sha256 = "c" * 64

    with pytest.raises(ValueError, match="incompatible schema"):
        study._simulation_invocation_history(  # noqa: SLF001
            specification_dir, expected_contract_sha256=contract_sha256
        )
    migrated = study._simulation_invocation_history(  # noqa: SLF001
        specification_dir,
        expected_contract_sha256=contract_sha256,
        allow_legacy_migration=True,
    )

    assert migrated["simulation_contract_sha256"].tolist() == [contract_sha256]
    assert migrated["effective_timeout_seconds"].tolist() == [20.0]
    persisted = pd.read_csv(history_path, sep="\t")
    assert list(persisted.columns) == list(study.SIMULATION_INVOCATION_COLUMNS)


def test_specification_lock_rejects_live_owner_and_recovers_stale_owner(
    tmp_path, monkeypatch
):
    specification_dir = tmp_path / "work" / "spec"
    with study._specification_execution_lock(specification_dir):  # noqa: SLF001
        with pytest.raises(study.SimulationSpecificationLocked, match="already locked"):
            with study._specification_execution_lock(  # noqa: SLF001
                specification_dir
            ):
                pytest.fail("a second live owner acquired the specification lock")
    assert not (specification_dir / ".simulation.lock").exists()

    stale_lock = specification_dir / ".simulation.lock"
    stale_lock.write_text(
        json.dumps(
            {
                "pid": 999_999_999,
                "hostname": socket.gethostname(),
                "created_at_utc": "2026-08-12T00:00:00+00:00",
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr(study, "_pid_is_alive", lambda pid: False)
    with study._specification_execution_lock(specification_dir):  # noqa: SLF001
        assert stale_lock.is_file()
    assert not stale_lock.exists()
    assert list(specification_dir.glob(".simulation.lock.stale.*"))


def test_atomic_replace_retries_access_denied_and_remains_fail_closed(
    tmp_path, monkeypatch
):
    destination = tmp_path / "state.json"
    real_replace = study.os.replace
    calls = []
    delays = []

    def transient_replace(source, target):
        calls.append((source, target))
        if len(calls) < 3:
            raise PermissionError(13, "transient OneDrive access denied")
        return real_replace(source, target)

    monkeypatch.setattr(study.os, "replace", transient_replace)
    monkeypatch.setattr(study, "sleep", delays.append)
    study._atomic_json(destination, {"status": "complete"})  # noqa: SLF001

    assert json.loads(destination.read_text(encoding="utf-8"))["status"] == "complete"
    assert len(calls) == 3
    assert delays == [0.05, 0.1]

    monkeypatch.setattr(
        study.os,
        "replace",
        lambda *args: (_ for _ in ()).throw(PermissionError(13, "still denied")),
    )
    with pytest.raises(PermissionError, match="still denied"):
        study._replace_with_access_retry(  # noqa: SLF001
            tmp_path / "source", tmp_path / "target", max_attempts=2
        )


def test_status_recovers_valid_partial_draw_progress(tmp_path, monkeypatch):
    runtime = _runtime(tmp_path, max_external_draws=2)
    specification = _accepted_specification()
    _install_fake_simulator(monkeypatch, _accepted_tiny_tree_sequence())

    class InterruptedEngine:
        def simulate(self, *args, **kwargs):
            raise RuntimeError("interrupted conditioned trajectory")

    monkeypatch.setattr(study.stdpopsim, "get_engine", lambda name: InterruptedEngine())
    specification_dir = (
        tmp_path
        / study.NO_INTROGRESSION_DIR
        / "work"
        / specification["specification_id"]
    )
    task = {
        "specification": specification,
        "runtime": asdict(runtime),
        "slim_path": sys.executable,
        "specification_dir": str(specification_dir),
    }
    with pytest.raises(RuntimeError, match="interrupted conditioned"):
        study._simulate_specification_task(task)  # noqa: SLF001
    contract = json.loads(
        (specification_dir / "simulation_contract.json").read_text(encoding="utf-8")
    )
    monkeypatch.setattr(
        study,
        "_recorded_simulation_contract_if_current",
        lambda *args, **kwargs: contract,
    )

    study._write_simulation_status_tables(  # noqa: SLF001
        runtime, {"no_introgression": [specification]}, {}
    )
    status = pd.read_csv(
        tmp_path / study.NO_INTROGRESSION_DIR / "results" / "simulation_status.tsv",
        sep="\t",
    ).iloc[0]

    assert status["status"] == "partial"
    assert status["completed_external_draws"] == 0
    assert status["launched_external_draws"] == 1
    assert "available for resume" in status["detail"]


def test_watchdog_prefers_completion_that_races_with_timeout(tmp_path, monkeypatch):
    runtime = _runtime(tmp_path, max_specification_seconds=1.0)
    specification = _accepted_specification()
    task = {
        "specification": specification,
        "runtime": asdict(runtime),
        "slim_path": sys.executable,
        "specification_dir": str(tmp_path / "work" / specification["specification_id"]),
    }

    class FakeQueue:
        def get_nowait(self):
            raise queue.Empty

        def close(self):
            return None

        def join_thread(self):
            return None

    class FakeProcess:
        pid = 12345
        exitcode = None

        def __init__(self, *args, **kwargs):
            self.alive = True

        def start(self):
            return None

        def join(self, timeout=None):
            return None

        def is_alive(self):
            return self.alive

    process = FakeProcess()

    class FakeContext:
        def Queue(self):
            return FakeQueue()

        def Process(self, *args, **kwargs):
            return process

    clock = iter([0.0, 2.0])
    monkeypatch.setattr(
        study.multiprocessing, "get_context", lambda name: FakeContext()
    )
    monkeypatch.setattr(study, "perf_counter", lambda: next(clock))
    terminated = []

    def terminate(fake_process):
        terminated.append(fake_process.pid)
        fake_process.alive = False

    monkeypatch.setattr(study, "_terminate_process_tree", terminate)
    completion_checks = iter([None, {"status": "complete", "tree_sha256": "valid"}])
    monkeypatch.setattr(
        study, "_valid_task_completion", lambda current_task: next(completion_checks)
    )
    monkeypatch.setattr(
        study,
        "_failed_simulation_result",
        lambda *args, **kwargs: pytest.fail("valid completion was overwritten"),
    )

    result = study._simulate_tasks_with_timeouts([task], runtime)  # noqa: SLF001

    assert terminated == [12345]
    assert result[specification["specification_id"]]["status"] == "complete"
    assert (
        result[specification["specification_id"]]["completion"]["tree_sha256"]
        == "valid"
    )


def test_status_and_decode_ignore_stale_simulation_contract(tmp_path, monkeypatch):
    runtime = _runtime(tmp_path)
    specification = _accepted_specification()
    _install_fake_simulator(monkeypatch, _accepted_tiny_tree_sequence())
    specification_dir = (
        tmp_path
        / study.NO_INTROGRESSION_DIR
        / "work"
        / specification["specification_id"]
    )
    task = {
        "specification": specification,
        "runtime": asdict(runtime),
        "slim_path": sys.executable,
        "specification_dir": str(specification_dir),
    }
    study._simulate_specification_task(task)  # noqa: SLF001
    changed_runtime = _runtime(tmp_path, output_stride_bp=20_000)

    study._write_simulation_status_tables(  # noqa: SLF001
        changed_runtime, {"no_introgression": [specification]}, {}
    )
    status_path = (
        tmp_path / study.NO_INTROGRESSION_DIR / "results" / "simulation_status.tsv"
    )
    status = pd.read_csv(status_path, sep="\t").iloc[0]
    assert not any(line.endswith("\t") for line in status_path.read_text().splitlines())
    assert status["status"] == "not_requested"
    assert "stale or incompatible" in status["detail"]
    assert status["completed_external_draws"] == 0
    monkeypatch.setattr(
        study,
        "run_within_decoder",
        lambda *args, **kwargs: pytest.fail("stale simulation reached the decoder"),
    )
    with pytest.raises(ValueError, match="incompatible with the current"):
        study.decode_specification(
            specification, changed_runtime, decoder_path=sys.executable
        )


def test_plot_decoder_validation_is_explicit_or_uses_recorded_contract(
    tmp_path, monkeypatch
):
    runtime = _runtime(tmp_path)
    specification = _accepted_specification()
    repo_decoder = tmp_path / "bin" / "gamma_smc"
    repo_decoder.parent.mkdir(parents=True)
    repo_decoder.write_bytes(b"repo decoder")
    custom_decoder = tmp_path / "custom_gamma_smc"
    custom_decoder.write_bytes(b"custom decoder")
    seen_decoder_paths = []

    monkeypatch.setattr(
        study,
        "_valid_recorded_simulation_cache",
        lambda *args, **kwargs: {"status": "complete"},
    )

    def fake_decode_cache(*args, decoder_path=None, **kwargs):
        seen_decoder_paths.append(decoder_path)
        if decoder_path is None or decoder_path == custom_decoder.resolve():
            return {"status": "complete"}
        return None

    monkeypatch.setattr(study, "_valid_recorded_decode_cache", fake_decode_cache)

    _, portable_decodes, portable_status_path = study._validated_plot_inputs(  # noqa: SLF001
        runtime,
        "no_introgression",
        [specification],
        decoder_path=None,
    )
    portable_status = pd.read_csv(portable_status_path, sep="\t").iloc[0]
    assert not any(
        line.endswith("\t") for line in portable_status_path.read_text().splitlines()
    )
    assert seen_decoder_paths[-1] is None
    assert specification["specification_id"] in portable_decodes
    assert portable_status["status"] == "complete"
    assert portable_status["decoder_validation_mode"] == "recorded_decode_contract"

    _, explicit_decodes, explicit_status_path = study._validated_plot_inputs(  # noqa: SLF001
        runtime,
        "no_introgression",
        [specification],
        decoder_path=custom_decoder,
    )
    explicit_status = pd.read_csv(explicit_status_path, sep="\t").iloc[0]
    assert seen_decoder_paths[-1] == custom_decoder.resolve()
    assert specification["specification_id"] in explicit_decodes
    assert explicit_status["status"] == "complete"
    assert explicit_status["decoder_validation_mode"] == "explicit_decoder_sha256"
    assert (
        explicit_status["expected_decoder_sha256"]
        == hashlib.sha256(custom_decoder.read_bytes()).hexdigest()
    )

    _, mismatched_decodes, mismatched_status_path = study._validated_plot_inputs(  # noqa: SLF001
        runtime,
        "no_introgression",
        [specification],
        decoder_path=repo_decoder,
    )
    mismatched_status = pd.read_csv(mismatched_status_path, sep="\t").iloc[0]
    assert not mismatched_decodes
    assert mismatched_status["status"] == "not_requested"


def test_main_passes_used_decoder_to_all_but_plot_without_override_is_portable(
    tmp_path, monkeypatch
):
    custom_decoder = tmp_path / "custom_gamma_smc"
    custom_decoder.write_bytes(b"custom decoder")
    monkeypatch.setattr(study, "write_study_plans", lambda *args, **kwargs: {})
    monkeypatch.setattr(study, "simulate_studies", lambda *args, **kwargs: {})
    monkeypatch.setattr(study, "decode_studies", lambda *args, **kwargs: {})
    plot_decoder_paths = []

    def fake_plot(*args, decoder_path=None, **kwargs):
        plot_decoder_paths.append(decoder_path)
        return {}

    monkeypatch.setattr(study, "plot_and_aggregate_studies", fake_plot)

    assert (
        study.main(
            [
                "all",
                "--repo-root",
                str(tmp_path),
                "--decoder-bin",
                str(custom_decoder),
            ]
        )
        == 0
    )
    assert plot_decoder_paths[-1] == custom_decoder

    assert study.main(["plot", "--repo-root", str(tmp_path)]) == 0
    assert plot_decoder_paths[-1] is None


def _write_decode_cache_fixture(
    tmp_path,
    *,
    legacy=False,
    threads=4,
):
    runtime = _runtime(tmp_path)
    specification = _accepted_specification()
    specification_dir = tmp_path / "work" / specification["specification_id"]
    decoded_dir = specification_dir / "decoded"
    pairs_dir = specification_dir / "pairs"
    decoded_dir.mkdir(parents=True)
    pairs_dir.mkdir(parents=True)
    tree_path = specification_dir / "selected.trees"
    pairs_path = pairs_dir / "overall.pairs.tsv"
    output = decoded_dir / "overall.summary.tsv"
    decoder_path = tmp_path / "bin" / "gamma_smc"
    decoder_path.parent.mkdir(parents=True)
    decoder_path.write_bytes(b"decoder fixture")
    tree_path.write_bytes(b"tree fixture")
    pairs_path.write_text("left\tright\n0\t1\n", encoding="utf-8")
    pd.DataFrame(
        {"position_0based": [5_000_000], "mean_p_tmrca_lt_threshold": [0.5]}
    ).to_csv(output, sep="\t", index=False)
    builder = (
        study._legacy_decode_contract_with_sha256  # noqa: SLF001
        if legacy
        else study._decode_contract_with_sha256  # noqa: SLF001
    )
    contract = builder(
        specification,
        runtime,
        tree_path=tree_path,
        pairs_path=pairs_path,
        decoder_sha256=hashlib.sha256(decoder_path.read_bytes()).hexdigest(),
        genotype_class="overall",
    )
    contract_path = output.with_suffix(output.suffix + ".contract.json")
    contract_path.write_text(json.dumps(contract), encoding="utf-8")
    thresholds = ",".join(str(value) for value in contract["thresholds_years"])
    producer_script = (
        "import sys; "
        "from gamma_smc_aou.tree_sequence import stream_tree_sequence_vcf; "
        "stream_tree_sequence_vcf(sys.argv[1], sys.stdout, input_format=sys.argv[2])"
    )
    command = [
        str(decoder_path),
        "--input",
        "/dev/stdin",
        "--input_format",
        "vcf",
        "--scaled_mutation_rate",
        str(contract["scaled_mutation_rate"]),
        "--recombination_to_mutation_ratio",
        str(contract["recombination_to_mutation_ratio"]),
        "--unscaled_mutation_rate",
        str(contract["unscaled_mutation_rate"]),
        "--recent_threshold_years",
        thresholds,
        "--generation_time",
        str(contract["generation_time_years"]),
        "--recent_summary",
        str(output),
        "--recent_call",
        "median",
        "--recent_call_probability",
        "0.5",
        "--output_at_hets=false",
        "--output_at_stride",
        str(contract["output_stride_bp"]),
        "--cache_size",
        str(contract["cache_size_bp"]),
        "--threads",
        str(threads),
        "--pair_block",
        "256",
        "--backward_alignment",
        "fixed",
        "--exp10",
        "accurate",
        "--pairs_file",
        str(pairs_path),
    ]
    run = {
        "command": command,
        "input_path": str(tree_path),
        "input_format": "trees",
        "tree_sequence_vcf_producer_command": [
            sys.executable,
            "-c",
            producer_script,
            str(tree_path),
            "trees",
        ],
        "stdout": "",
        "stderr": "",
        "decode_seconds": 1.0,
        "stride_bp": contract["output_stride_bp"],
        "cache_size_bp": contract["cache_size_bp"],
        "n_output_positions": 1,
        "pairs_manifest": None,
        "n_pairs_recorded": None,
    }
    run_path = output.with_suffix(output.suffix + ".run.json")
    run_path.write_text(json.dumps(run), encoding="utf-8")
    completion = {
        "contract_sha256": study._contract_sha256(contract),  # noqa: SLF001
        "summary_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
        "run_sha256": hashlib.sha256(run_path.read_bytes()).hexdigest(),
    }
    if not legacy:
        completion["execution"] = {
            "threads": threads,
            "threads_cache_role": "provenance_only",
            "run_settings_validated": True,
        }
    completion_path = output.with_suffix(output.suffix + ".complete.json")
    completion_path.write_text(json.dumps(completion), encoding="utf-8")
    return {
        "runtime": runtime,
        "specification": specification,
        "tree_path": tree_path,
        "pairs_path": pairs_path,
        "output": output,
        "decoder_path": decoder_path,
        "contract_path": contract_path,
        "completion_path": completion_path,
        "run_path": run_path,
    }


def test_portable_path_matches_wsl_drvfs_and_windows_spelling():
    windows = (
        r"C:\Users\Lenovo\OneDrive\Documents\cxt_demo\gamma_smc_ts"
        r"\work\intro_s0p001_af10__rep001\decoded\overall.summary.tsv"
    )
    wsl = (
        "/mnt/c/Users/Lenovo/OneDrive/Documents/cxt_demo/gamma_smc_ts/"
        "work/intro_s0p001_af10__rep001/decoded/overall.summary.tsv"
    )

    assert study._portable_path_matches(wsl, windows)  # noqa: SLF001
    assert study._portable_path_matches(  # noqa: SLF001
        wsl, windows.replace("C:", "c:", 1)
    )
    assert not study._portable_path_matches(  # noqa: SLF001
        wsl.replace("/mnt/c/", "/mnt/d/"), windows
    )
    assert not study._portable_path_matches(  # noqa: SLF001
        wsl.replace("rep001", "rep002"), windows
    )
    assert not study._portable_path_matches(  # noqa: SLF001
        wsl.replace("/mnt/c/", "/mnt/C/"), windows
    )
    assert not study._portable_path_matches(  # noqa: SLF001
        wsl.replace("/Users/", "/users/"), windows
    )
    assert not study._portable_path_matches(  # noqa: SLF001
        "work/intro_s0p001_af10__rep001/decoded/overall.summary.tsv", windows
    )
    assert not study._portable_path_matches(  # noqa: SLF001
        "/unrelated/work/intro_s0p001_af10__rep001/decoded/overall.summary.tsv",
        windows,
    )
    assert not study._portable_path_matches(  # noqa: SLF001
        wsl.replace("/work/", "/work/../work/"), windows
    )


@pytest.mark.skipif(os.name != "nt", reason="requires native Windows path spellings")
def test_decode_cache_accepts_equivalent_wsl_run_paths(tmp_path):
    fixture = _write_decode_cache_fixture(tmp_path)
    run = json.loads(fixture["run_path"].read_text(encoding="utf-8"))

    def as_wsl(path):
        resolved = Path(path).resolve()
        return f"/mnt/{resolved.drive[0].lower()}{resolved.as_posix()[2:]}"

    run["command"][0] = as_wsl(fixture["decoder_path"])
    run["command"][run["command"].index("--recent_summary") + 1] = as_wsl(
        fixture["output"]
    )
    run["command"][run["command"].index("--pairs_file") + 1] = as_wsl(
        fixture["pairs_path"]
    )
    run["tree_sequence_vcf_producer_command"][3] = as_wsl(fixture["tree_path"])
    run["input_path"] = as_wsl(fixture["tree_path"])
    fixture["run_path"].write_text(json.dumps(run), encoding="utf-8")
    completion = json.loads(fixture["completion_path"].read_text(encoding="utf-8"))
    completion["run_sha256"] = hashlib.sha256(
        fixture["run_path"].read_bytes()
    ).hexdigest()
    fixture["completion_path"].write_text(json.dumps(completion), encoding="utf-8")

    contract = json.loads(fixture["contract_path"].read_text(encoding="utf-8"))
    assert study._valid_decode_cache(  # noqa: SLF001
        fixture["output"],
        contract,
        tree_path=fixture["tree_path"],
        pairs_path=fixture["pairs_path"],
        decoder_path=fixture["decoder_path"],
    )


def test_decode_lock_is_separate_and_recovers_only_a_dead_local_owner(
    tmp_path, monkeypatch
):
    specification_dir = tmp_path / "work" / "spec"
    with study._decode_execution_lock(specification_dir):  # noqa: SLF001
        assert (specification_dir / ".decode.lock").is_file()
        assert not (specification_dir / ".simulation.lock").exists()
        with pytest.raises(study.DecodeSpecificationLocked, match="already locked"):
            with study._decode_execution_lock(specification_dir):  # noqa: SLF001
                pytest.fail("concurrent decode acquired a live lock")
    assert not (specification_dir / ".decode.lock").exists()

    stale_lock = specification_dir / ".decode.lock"
    stale_lock.write_text(
        json.dumps(
            {
                "pid": 999_999_999,
                "hostname": socket.gethostname(),
                "created_at_utc": "2026-08-12T00:00:00+00:00",
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr(study, "_pid_is_alive", lambda pid: False)
    with study._decode_execution_lock(specification_dir):  # noqa: SLF001
        assert stale_lock.is_file()
    assert list(specification_dir.glob(".decode.lock.stale.*"))


def test_decode_contract_binds_fixed_settings_and_threads_are_provenance_only(
    tmp_path,
):
    fixture = _write_decode_cache_fixture(tmp_path, threads=4)
    contract = json.loads(fixture["contract_path"].read_text(encoding="utf-8"))
    assert contract["schema"] == study.DECODE_CONTRACT_SCHEMA
    assert contract["decode_settings"] == {
        **study.DECODE_FIXED_SETTINGS,
        "threads": {
            "cache_role": "provenance_only",
            "required_minimum": 1,
            "determinism_contract": "output_invariant_across_thread_count",
        },
    }
    assert contract["python_vcf_producer"]["metadata_status"] == "recorded"
    assert study._valid_decode_cache(  # noqa: SLF001
        fixture["output"],
        contract,
        tree_path=fixture["tree_path"],
        pairs_path=fixture["pairs_path"],
        decoder_path=fixture["decoder_path"],
    )

    run = json.loads(fixture["run_path"].read_text(encoding="utf-8"))
    thread_index = run["command"].index("--threads") + 1
    run["command"][thread_index] = "8"
    fixture["run_path"].write_text(json.dumps(run), encoding="utf-8")
    completion = json.loads(fixture["completion_path"].read_text(encoding="utf-8"))
    completion["run_sha256"] = hashlib.sha256(
        fixture["run_path"].read_bytes()
    ).hexdigest()
    completion["execution"] = {
        "threads": 8,
        "threads_cache_role": "provenance_only",
        "run_settings_validated": True,
    }
    fixture["completion_path"].write_text(json.dumps(completion), encoding="utf-8")
    assert study._valid_decode_cache(  # noqa: SLF001
        fixture["output"],
        contract,
        tree_path=fixture["tree_path"],
        pairs_path=fixture["pairs_path"],
        decoder_path=fixture["decoder_path"],
    )

    run["command"][run["command"].index("--exp10") + 1] = "fast"
    fixture["run_path"].write_text(json.dumps(run), encoding="utf-8")
    completion["run_sha256"] = hashlib.sha256(
        fixture["run_path"].read_bytes()
    ).hexdigest()
    fixture["completion_path"].write_text(json.dumps(completion), encoding="utf-8")
    assert not study._valid_decode_cache(  # noqa: SLF001
        fixture["output"],
        contract,
        tree_path=fixture["tree_path"],
        pairs_path=fixture["pairs_path"],
        decoder_path=fixture["decoder_path"],
    )


def test_legacy_decode_cache_migrates_only_after_run_semantics_validate(tmp_path):
    fixture = _write_decode_cache_fixture(tmp_path, legacy=True)
    migrated = study._validated_or_migrated_class_decode_cache(  # noqa: SLF001
        fixture["specification"],
        fixture["runtime"],
        tree_path=fixture["tree_path"],
        pairs_path=fixture["pairs_path"],
        output=fixture["output"],
        decoder_sha256=hashlib.sha256(fixture["decoder_path"].read_bytes()).hexdigest(),
        genotype_class="overall",
        decoder_path=fixture["decoder_path"],
        migrate_legacy=True,
    )
    assert migrated is not None
    contract, execution, was_migrated = migrated
    assert was_migrated is True
    assert execution["threads"] == 4
    assert contract["schema"] == study.DECODE_CONTRACT_SCHEMA
    assert contract["python_vcf_producer"]["metadata_status"] == "legacy_not_recorded"
    assert contract["migration"]["validation"] == (
        "legacy_contract_hash_and_run_semantics"
    )
    assert (
        json.loads(fixture["completion_path"].read_text(encoding="utf-8"))["execution"]
        == execution
    )

    invalid = _write_decode_cache_fixture(tmp_path / "invalid", legacy=True)
    run = json.loads(invalid["run_path"].read_text(encoding="utf-8"))
    run["command"][run["command"].index("--pair_block") + 1] = "128"
    invalid["run_path"].write_text(json.dumps(run), encoding="utf-8")
    invalid_completion = json.loads(
        invalid["completion_path"].read_text(encoding="utf-8")
    )
    invalid_completion["run_sha256"] = hashlib.sha256(
        invalid["run_path"].read_bytes()
    ).hexdigest()
    invalid["completion_path"].write_text(
        json.dumps(invalid_completion), encoding="utf-8"
    )
    assert (
        study._validated_or_migrated_class_decode_cache(  # noqa: SLF001
            invalid["specification"],
            invalid["runtime"],
            tree_path=invalid["tree_path"],
            pairs_path=invalid["pairs_path"],
            output=invalid["output"],
            decoder_sha256=hashlib.sha256(
                invalid["decoder_path"].read_bytes()
            ).hexdigest(),
            genotype_class="overall",
            decoder_path=invalid["decoder_path"],
            migrate_legacy=True,
        )
        is None
    )
    assert (
        json.loads(invalid["contract_path"].read_text(encoding="utf-8"))["schema"]
        == study.SCHEMA_VERSION
    )


def test_simulation_status_prefers_valid_completion_over_invocation_failure(
    tmp_path, monkeypatch
):
    runtime = _runtime(tmp_path)
    specification = _accepted_specification()
    _install_fake_simulator(monkeypatch, _accepted_tiny_tree_sequence())
    specification_dir = (
        tmp_path
        / study.NO_INTROGRESSION_DIR
        / "work"
        / specification["specification_id"]
    )
    task = {
        "specification": specification,
        "runtime": asdict(runtime),
        "slim_path": sys.executable,
        "specification_dir": str(specification_dir),
    }
    study._simulate_specification_task(task)  # noqa: SLF001
    recorded_contract = json.loads(
        (specification_dir / "simulation_contract.json").read_text(encoding="utf-8")
    )
    monkeypatch.setattr(
        study,
        "_recorded_simulation_contract_if_current",
        lambda *args, **kwargs: recorded_contract,
    )
    study._write_simulation_status_tables(  # noqa: SLF001
        runtime,
        {"no_introgression": [specification]},
        {
            specification["specification_id"]: {
                "status": "failed",
                "failure": {"error": "concurrent invocation saw a simulation lock"},
            }
        },
    )
    status = pd.read_csv(
        tmp_path / study.NO_INTROGRESSION_DIR / "results" / "simulation_status.tsv",
        sep="\t",
    ).iloc[0]
    assert status["status"] == "complete"
    assert status["latest_invocation_status"] == "failed"
    assert status["latest_invocation_error"] == (
        "concurrent invocation saw a simulation lock"
    )


def test_source_epoch_mismatch_fails_before_work_directory_mutation(tmp_path):
    runtime = _runtime(tmp_path)
    specification = _accepted_specification()
    specification_dir = tmp_path / "work" / specification["specification_id"]
    expected = study._implementation_sha256s(  # noqa: SLF001
        study._implementation_provenance()  # noqa: SLF001
    )
    expected["core"] = "0" * 64
    task = {
        "specification": specification,
        "runtime": asdict(runtime),
        "slim_path": sys.executable,
        "specification_dir": str(specification_dir),
        "expected_implementation_sha256s": expected,
    }
    with pytest.raises(
        study.ImplementationProvenanceMismatch, match="phase-start source snapshot"
    ):
        study._simulate_specification_task(task)  # noqa: SLF001
    assert not specification_dir.exists()


def test_external_interruption_recovery_uses_lock_interval_and_refuses_live_owner(
    tmp_path, monkeypatch
):
    specification_dir = tmp_path / "work" / "spec"
    specification_dir.mkdir(parents=True)
    lock_path = specification_dir / ".simulation.lock"
    lock_path.write_text(
        json.dumps(
            {
                "pid": 12345,
                "hostname": socket.gethostname(),
                "created_at_utc": "2026-08-12T12:00:00+00:00",
            }
        ),
        encoding="utf-8",
    )
    task = {
        "specification": {"specification_id": "spec"},
        "specification_dir": str(specification_dir),
    }
    monkeypatch.setattr(study, "_valid_task_completion", lambda task: None)
    monkeypatch.setattr(study, "_pid_is_alive", lambda pid: False)
    recorded = {}

    def fake_failure(task, **kwargs):
        recorded.update(kwargs)
        return {"status": kwargs["status"]}

    monkeypatch.setattr(study, "_failed_simulation_result", fake_failure)
    result = study.record_external_simulation_interruption(
        task, stopped_at_utc="2026-08-12T12:03:30+00:00"
    )
    assert result["status"] == "failed"
    assert recorded["elapsed_seconds"] == 210.0
    assert recorded["error_type"] == "ExternalProcessInterruption"

    monkeypatch.setattr(study, "_pid_is_alive", lambda pid: True)
    with pytest.raises(study.SimulationSpecificationLocked, match="live"):
        study.record_external_simulation_interruption(
            task, stopped_at_utc="2026-08-12T12:03:30+00:00"
        )


def test_parser_defaults_match_execution_contract(tmp_path):
    args = study.build_parser().parse_args(["plan", "--repo-root", str(tmp_path)])

    assert args.phase == "plan"
    assert args.scenario == "both"
    assert args.sample_diploids == study.DEFAULT_SAMPLE_DIPLOIDS == 100
    assert args.workers == study.DEFAULT_SIMULATION_WORKERS == 4
    assert args.threads == study.DEFAULT_THREADS == 4
    assert args.slim_scaling_factor == study.DEFAULT_SLIM_SCALING_FACTOR == 10.0
    assert args.slim_burn_in == study.DEFAULT_SLIM_BURN_IN
    assert args.max_external_draws == study.DEFAULT_MAX_EXTERNAL_DRAWS == 20
    assert args.max_launched_draws == 0
    assert args.spec_timeout_minutes == 0.0
    assert args.cumulative_spec_timeout_minutes == 0.0
    assert args.spec == []
    assert study._resolve_scenarios(args.scenario) == (  # noqa: SLF001
        "no_introgression",
        "introgression",
    )


def test_runtime_rejects_invalid_execution_values(tmp_path):
    with pytest.raises(ValueError, match="threads must be positive"):
        _runtime(tmp_path, threads=0).validate()
    with pytest.raises(ValueError, match="exactly 10 Mb"):
        _runtime(tmp_path, sequence_length_bp=100).validate()
    with pytest.raises(ValueError, match="max_specification_seconds"):
        _runtime(tmp_path, max_specification_seconds=-1).validate()
    with pytest.raises(ValueError, match="max_launched_draws"):
        _runtime(tmp_path, max_launched_draws=-1).validate()
    with pytest.raises(ValueError, match="max_cumulative_specification_seconds"):
        _runtime(tmp_path, max_cumulative_specification_seconds=-1).validate()
    with pytest.raises(ValueError, match="campaign_subdirectory"):
        _runtime(tmp_path, campaign_subdirectory="../escape").validate()
