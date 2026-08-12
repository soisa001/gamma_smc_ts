from __future__ import annotations

import hashlib
import json
import sys
from dataclasses import asdict
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
    assert args.spec_timeout_minutes == 0.0
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
