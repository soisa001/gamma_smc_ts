from __future__ import annotations

import copy
import json
import math
from contextlib import nullcontext
from pathlib import Path
from types import SimpleNamespace

import msprime
import numpy as np
import pyslim
import pytest
import stdpopsim
import tskit
from gamma_smc_aou import focused_selection_campaign as campaign
from gamma_smc_aou import focused_selection_simulation as simulation


def _slim_focal_test_tree(
    *,
    mutation_type: int | None,
    subpopulation: int | None = None,
    sample_diploids: int = 4,
):
    ancestry = msprime.sim_ancestry(
        samples=sample_diploids,
        ploidy=2,
        sequence_length=10_000_000,
        recombination_rate=0,
        population_size=10_000,
        random_seed=23,
    )
    annotated = pyslim.annotate(
        ancestry,
        model_type="WF",
        tick=5_000,
        stage="late",
    )
    tables = annotated.dump_tables()
    top_metadata = copy.deepcopy(tables.metadata)
    top_metadata["SLiM"]["user_metadata"] = {"Q": [10.0]}
    tables.metadata = top_metadata
    if mutation_type is None:
        return tables.tree_sequence()
    site = tables.sites.add_row(position=5_000_000, ancestral_state="A")
    tables.mutations.add_row(
        site=site,
        node=int(annotated.samples()[0]),
        time=10,
        derived_state="G",
        metadata={
            "mutation_list": [
                {
                    "mutation_type": mutation_type,
                    "selection_coeff": 0.0,
                    "subpopulation": (
                        int(subpopulation)
                        if subpopulation is not None
                        else (0 if mutation_type == 1 else -1)
                    ),
                    "slim_time": 4_999,
                    "nucleotide": 2,
                }
            ]
        },
    )
    tables.sort()
    return tables.tree_sequence()


def _test_selected_focal_expectation() -> dict:
    return {
        "schema": simulation.SELECTED_FOCAL_IDENTITY_SCHEMA,
        "single_site_id": "eas_selected_site",
        "focal_position_0based": 5_000_000,
        "expected_mutation_type": 1,
        "expected_selection_coeff": 0.0,
        "source_population": "EAS",
        "source_subpopulation": 0,
        "requested_origin_time_generations": 10.0,
        "realized_origin_time_generations": 10.0,
        "slim_scaling_factor": 10.0,
        "slim_time_rule": "cycle - realized_origin_time / Q",
    }


def _neutral_semantic_test_tree():
    tables = tskit.TableCollection(sequence_length=10_000_000)
    population = tables.populations.add_row()
    samples = []
    for _ in range(6):
        individual = tables.individuals.add_row()
        samples.extend(
            tables.nodes.add_row(
                flags=tskit.NODE_IS_SAMPLE,
                time=0,
                population=population,
                individual=individual,
            )
            for _ in range(2)
        )
    selected_ancestor = tables.nodes.add_row(time=1, population=population)
    root = tables.nodes.add_row(time=2, population=population)
    for sample in samples[:5]:
        tables.edges.add_row(0, tables.sequence_length, selected_ancestor, sample)
    tables.edges.add_row(0, tables.sequence_length, root, selected_ancestor)
    for sample in samples[5:]:
        tables.edges.add_row(0, tables.sequence_length, root, sample)
    site = tables.sites.add_row(position=5_000_000, ancestral_state="A")
    tables.mutations.add_row(
        site=site,
        node=selected_ancestor,
        time=1.5,
        derived_state="G",
    )
    tables.sort()
    return tables.tree_sequence()


def _introgressed_neutral_test_tree(
    *,
    mutation_time: float = 2_400,
    mutation_population: str = "Neanderthal",
    include_pulse: bool = True,
):
    tables = tskit.TableCollection(sequence_length=10_000_000)
    tables.populations.metadata_schema = tskit.MetadataSchema.permissive_json()
    populations = {
        name: tables.populations.add_row(metadata={"name": name, "id": name})
        for name in ("Han", "Loschbour", "Neanderthal")
    }
    samples = []
    for _ in range(8):
        individual = tables.individuals.add_row()
        samples.extend(
            tables.nodes.add_row(
                flags=tskit.NODE_IS_SAMPLE,
                time=0,
                population=populations["Han"],
                individual=individual,
            )
            for _ in range(2)
        )
    pulse_lineage = tables.nodes.add_row(
        time=2_016, population=populations["Loschbour"]
    )
    selected_ancestor = tables.nodes.add_row(
        time=2_300, population=populations[mutation_population]
    )
    root = tables.nodes.add_row(time=30_000, population=populations["Loschbour"])
    selected_samples = [samples[index] for index in (0, 1, 2, 3, 4, 6)]
    for sample in selected_samples:
        tables.edges.add_row(0, tables.sequence_length, pulse_lineage, sample)
    tables.edges.add_row(0, tables.sequence_length, selected_ancestor, pulse_lineage)
    tables.edges.add_row(0, tables.sequence_length, root, selected_ancestor)
    for sample in samples:
        if sample in selected_samples:
            continue
        tables.edges.add_row(0, tables.sequence_length, root, sample)
    if include_pulse:
        tables.migrations.add_row(
            0,
            tables.sequence_length,
            pulse_lineage,
            populations["Loschbour"],
            populations["Neanderthal"],
            2_272,
        )
    site = tables.sites.add_row(position=5_000_000, ancestral_state="A")
    tables.mutations.add_row(
        site=site,
        node=selected_ancestor,
        time=mutation_time,
        derived_state="G",
    )
    tables.sort()
    return tables.tree_sequence()


def _introgressed_candidate_pool_tree():
    """One 500-diploid Han tree with a panel-feasible 18% pulse branch."""

    tables = tskit.TableCollection(sequence_length=20_000_000)
    tables.populations.metadata_schema = tskit.MetadataSchema.permissive_json()
    populations = {
        name: tables.populations.add_row(metadata={"name": name, "id": name})
        for name in ("Han", "Loschbour", "Neanderthal")
    }
    samples = []
    for _ in range(500):
        individual = tables.individuals.add_row()
        samples.append(
            (
                tables.nodes.add_row(
                    flags=tskit.NODE_IS_SAMPLE,
                    time=0,
                    population=populations["Han"],
                    individual=individual,
                ),
                tables.nodes.add_row(
                    flags=tskit.NODE_IS_SAMPLE,
                    time=0,
                    population=populations["Han"],
                    individual=individual,
                ),
            )
        )
    pulse_lineage = tables.nodes.add_row(
        time=2_016, population=populations["Loschbour"]
    )
    origin_branch = tables.nodes.add_row(
        time=2_300, population=populations["Neanderthal"]
    )
    root = tables.nodes.add_row(time=30_000, population=populations["Loschbour"])
    selected_samples = {node for pair in samples[:80] for node in pair} | {
        pair[0] for pair in samples[80:100]
    }
    for pair in samples:
        for sample in pair:
            parent = pulse_lineage if sample in selected_samples else root
            tables.edges.add_row(0, tables.sequence_length, parent, sample)
    tables.edges.add_row(0, tables.sequence_length, origin_branch, pulse_lineage)
    tables.edges.add_row(0, tables.sequence_length, root, origin_branch)
    tables.migrations.add_row(
        0,
        tables.sequence_length,
        pulse_lineage,
        populations["Loschbour"],
        populations["Neanderthal"],
        2_272,
    )
    tables.sort()
    return tables.tree_sequence()


def _unit(tmp_path: Path, demography_id: str, af: float) -> dict:
    plan = campaign.FocusedCampaignPlan(
        repo_root=tmp_path,
        campaign_dir=tmp_path / campaign.DEFAULT_CAMPAIGN_DIR,
        selected_replicates=1,
        neutral_replicates=1,
    )
    frame = campaign.build_execution_units(plan)
    return (
        frame[
            (frame["demography_id"] == demography_id)
            & (frame["target_allele_frequency"] == af)
            & (frame["simulation_class"] == "neutral")
        ]
        .iloc[0]
        .to_dict()
    )


def _write_frozen_han_calibration(tmp_path: Path) -> tuple[Path, Path]:
    campaign_dir = tmp_path / campaign.DEFAULT_CAMPAIGN_DIR
    campaign_dir.mkdir(parents=True, exist_ok=True)
    plan = campaign.FocusedCampaignPlan(
        repo_root=tmp_path,
        campaign_dir=campaign_dir,
        selected_replicates=1,
        neutral_replicates=1,
    )
    units = campaign.build_execution_units(plan)
    units_path = campaign_dir / "execution_units.tsv"
    units.to_csv(units_path, sep="\t", index=False, lineterminator="\n")
    calibration_dir = campaign_dir / "calibration"
    calibration_dir.mkdir(parents=True)
    selected_rows = []
    mapping = {}
    for coefficient in (0.005, 0.01):
        for target, endpoint in zip((0.1, 0.2, 0.3), (1_800, 1_600, 1_400)):
            key = simulation._han_endpoint_key(coefficient, target)
            entry = {
                "selection_coefficient": coefficient,
                "target_allele_frequency": target,
                "duration_multiplier": 0.5,
                "requested_selection_end_generations_ago": endpoint + 1.2,
                "realized_selection_end_generations_ago": endpoint,
                "band_hit_rate": 0.25,
                "n_band_hits": 5,
                "n_evaluable": 20,
            }
            mapping[key] = entry
            selected_rows.append(entry)
    selected = simulation.pd.DataFrame(selected_rows)
    artifacts = {
        "han_selection_end_plan.tsv": simulation.pd.DataFrame(
            {"calibration_id": ["test"]}
        ),
        "han_selection_end_terminal_draws.tsv": simulation.pd.DataFrame(
            {"calibration_id": ["test"], "evaluable": [True]}
        ),
        "han_selection_end_summary.tsv": simulation.pd.DataFrame(
            {"selection_coefficient": [0.01], "band_hit_rate": [0.25]}
        ),
        "han_selection_end_selected.tsv": selected,
    }
    hashes = {}
    for name, frame in artifacts.items():
        path = calibration_dir / name
        frame.to_csv(path, sep="\t", index=False, lineterminator="\n")
        hashes[name] = simulation.sha256_file(path)
    frozen = calibration_dir / "han_selection_end_frozen.json"
    payload = {
        "schema": simulation.HAN_SELECTION_CALIBRATION_SCHEMA,
        "status": "frozen",
        "minimum_band_hit_rate": 0.20,
        "production_seed_digest": simulation._production_seed_digest(units),
        "production_contract": {
            "schema": simulation.HAN_SELECTION_PRODUCTION_CONTRACT_SCHEMA,
            "slim_scaling_factor": 5.0,
            "candidate_pool_diploids": 500,
            "af_half_width": 0.025,
            "selection_coefficients": [0.01, 0.005],
            "target_allele_frequencies": [0.1, 0.2, 0.3],
            "selection_endpoint_mode": (
                "fixed_endpoint_from_disjoint_terminal_af_calibration"
            ),
            "post_pulse_recipient_nonloss_conditioning": True,
        },
        "mapping": mapping,
        "artifacts": hashes,
    }
    frozen.write_text(json.dumps(payload), encoding="utf-8")
    return frozen, units_path


def test_frozen_han_calibration_binds_q_seed_grid_and_endpoint_table(tmp_path):
    frozen, units_path = _write_frozen_han_calibration(tmp_path)
    binding = simulation.load_frozen_han_selection_calibration(
        frozen,
        repo_root=tmp_path,
        execution_units_path=units_path,
        slim_scaling_factor=5.0,
        candidate_pool_diploids=500,
    )
    endpoint = simulation.resolve_han_selection_endpoint(
        binding, selection_coefficient=0.01, target_frequency=0.2
    )
    assert binding["sha256"] == simulation.sha256_file(frozen)
    assert endpoint["selection_end_generations_ago"] == 1_600
    assert (
        endpoint["selected_endpoint_table_sha256"]
        == binding["artifact_sha256"]["han_selection_end_selected.tsv"]
    )

    with pytest.raises(ValueError, match="Q differs"):
        simulation.load_frozen_han_selection_calibration(
            frozen,
            repo_root=tmp_path,
            execution_units_path=units_path,
            slim_scaling_factor=10.0,
            candidate_pool_diploids=500,
        )


def test_frozen_han_calibration_fails_closed_on_endpoint_table_tampering(tmp_path):
    frozen, units_path = _write_frozen_han_calibration(tmp_path)
    selected = frozen.parent / "han_selection_end_selected.tsv"
    selected.write_text(selected.read_text(encoding="utf-8") + "#tampered\n")
    with pytest.raises(ValueError, match="checksum failed"):
        simulation.load_frozen_han_selection_calibration(
            frozen,
            repo_root=tmp_path,
            execution_units_path=units_path,
            slim_scaling_factor=5.0,
            candidate_pool_diploids=500,
        )


def test_exact_af_and_shared_neutral_grid_contract(tmp_path):
    plan = campaign.FocusedCampaignPlan(
        repo_root=tmp_path,
        campaign_dir=tmp_path / campaign.DEFAULT_CAMPAIGN_DIR,
    )
    units = campaign.build_execution_units(plan)
    neutral = units[units["simulation_class"] == "neutral"]
    selected = units[units["simulation_class"] == "selected"]

    assert set(neutral["selection_coefficient"]) == {0.0}
    assert set(selected["selection_coefficient"]) == {0.005, 0.01}
    assert set(neutral["exact_sample_alt_count"]) == {20, 40, 60}
    assert len(neutral) == 600
    assert len(selected) == 120


def test_shared_candidate_pool_gate_has_inclusive_integer_bounds(tmp_path):
    unit = _unit(tmp_path, "eas_phlash_median", 0.1)

    def counts_for_alt_count(alt_count: int) -> np.ndarray:
        n_hom_alt = 2
        n_het = alt_count - 2 * n_hom_alt
        n_hom_ref = 500 - n_hom_alt - n_het
        return np.asarray(
            [0] * n_hom_ref + [1] * n_het + [2] * n_hom_alt,
            dtype=int,
        )

    lower = simulation.candidate_pool_frequency_gate(unit, counts_for_alt_count(75))
    upper = simulation.candidate_pool_frequency_gate(unit, counts_for_alt_count(125))
    below = simulation.candidate_pool_frequency_gate(unit, counts_for_alt_count(74))
    above = simulation.candidate_pool_frequency_gate(unit, counts_for_alt_count(126))

    assert lower["minimum_alt_count_inclusive"] == 75
    assert lower["maximum_alt_count_inclusive"] == 125
    assert lower["passed"] is True
    assert upper["passed"] is True
    assert below["passed"] is False
    assert above["passed"] is False
    triples = simulation.feasible_exact_panel_genotype_triples(
        counts_for_alt_count(75), sample_diploids=100, exact_alt_count=20
    )
    assert [
        (record["hom_ref"], record["heterozygous"], record["hom_alt"])
        for record in triples
    ] == [(82, 16, 2)]
    assert math.isfinite(triples[0]["log_number_of_subsets"])
    assert (
        simulation.feasible_exact_panel_genotype_triples(
            [0] * 450 + [1] * 50,
            sample_diploids=100,
            exact_alt_count=20,
        )
        == []
    )


def test_neutral_candidates_use_500_pool_band_and_panel_feasibility(tmp_path):
    unit = _unit(tmp_path, "eas_phlash_median", 0.1)
    demography = msprime.Demography()
    demography.add_population(name="EAS", initial_size=10_000)
    ancestry = msprime.sim_ancestry(
        samples={"EAS": 500},
        ploidy=2,
        demography=demography,
        sequence_length=10_000_000,
        recombination_rate=0,
        random_seed=43,
        model="smc_prime",
    )

    candidates = simulation.find_neutral_candidates(ancestry, unit, demography)

    assert candidates
    assert all(
        sum(candidate.candidate_pool_genotype_counts) == 500 for candidate in candidates
    )
    assert all(
        75 <= candidate.candidate_pool_alt_count <= 125 for candidate in candidates
    )
    assert all(
        candidate.eligible_genotype_count_triples > 0 for candidate in candidates
    )
    assert all(candidate.left == 5_000_000 for candidate in candidates)
    assert all(candidate.right == 5_000_001 for candidate in candidates)
    candidate = candidates[0]
    tables = ancestry.dump_tables()
    site = tables.sites.add_row(position=5_000_000, ancestral_state="A")
    tables.mutations.add_row(
        site=site,
        node=candidate.node,
        derived_state="G",
        time=(candidate.mutation_time_lower + candidate.mutation_time_upper) / 2,
    )
    tables.sort()
    pool_ts = tables.tree_sequence()
    pool_manifest = simulation.build_diploid_pair_table(pool_ts, 5_000_000)
    pool_gate = simulation.candidate_pool_frequency_gate(
        unit, pool_manifest["focal_selected_allele_count"]
    )
    panel_ts, panel_record = simulation.select_exact_sample_panel(
        pool_ts,
        sample_diploids=100,
        exact_alt_count=20,
        seed=19,
    )
    panel_manifest = simulation.build_diploid_pair_table(panel_ts, 5_000_000)

    assert pool_gate["passed"] is True
    assert panel_record["candidate_pool_diploids"] == 500
    assert len(panel_manifest) == 100
    assert int(panel_manifest["focal_selected_allele_count"].sum()) == 20
    assert set(panel_manifest["genotype_class"]) == {
        "hom_ref",
        "heterozygous",
        "hom_alt",
    }


def test_eas_neutral_constructor_produces_real_exact_focal_variant(repo_root, tmp_path):
    unit = _unit(repo_root, "eas_phlash_median", 0.1)
    ts, pairs, provenance, attempts = simulation.simulate_neutral_unit(
        unit,
        repo_root,
        tmp_path / "neutral-unit",
        max_attempts=10,
    )

    assert ts.sequence_length == 10_000_000
    assert int(pairs["focal_selected_allele_count"].sum()) == 20
    assert set(pairs["genotype_class"]) == {"hom_ref", "heterozygous", "hom_alt"}
    assert provenance["mutation_origin_population"] == "EAS"
    assert provenance["archaic_specific_no_ils"] is False
    assert provenance["candidate_pool_diploids"] == 500
    assert provenance["candidate_pool_frequency_gate"]["passed"] is True
    assert provenance["sample_panel"]["candidate_pool_diploids"] == 500
    assert provenance["model"] == "DTWF200ThenSmcPrimeApproxCoalescent"
    assert provenance["model_keywords"] == ["dtwf", "smc_prime"]
    assert provenance["recent_dtwf_duration_generations"] == 200
    assert provenance["ascertainment_length_bp"] == 10_000_000
    assert provenance["candidate_index"] == "marginal_tree_at_fixed_5Mb_site"
    assert any(record["accepted"] for record in attempts)


def test_introgression_candidates_require_pulse_migration(repo_root):
    unit = _unit(repo_root, "ancient_eurasia_han_introgression", 0.2)
    demography = msprime.Demography()
    for name in ("Han", "Loschbour", "Neanderthal"):
        demography.add_population(name=name, initial_size=10_000)
    ancestry = _introgressed_candidate_pool_tree()
    candidates = simulation.find_neutral_candidates(ancestry, unit, demography)

    assert candidates
    assert all(item.mutation_origin_population == "Neanderthal" for item in candidates)
    assert all(item.mutation_time_lower == 2_400 for item in candidates)
    assert all(item.mutation_time_upper == 2_400 for item in candidates)
    assert all(175 <= item.candidate_pool_alt_count <= 225 for item in candidates)
    assert all(sum(item.candidate_pool_genotype_counts) == 500 for item in candidates)
    assert all(item.eligible_genotype_count_triples > 0 for item in candidates)
    assert sum(item.weight for item in candidates) == 10_000_000
    assert all(item.n_introgressed_descendant_lineages == 1 for item in candidates)
    tables = ancestry.dump_tables()
    tables.migrations.clear()
    assert (
        simulation.find_neutral_candidates(tables.tree_sequence(), unit, demography)
        == []
    )


def test_introgressed_neutral_identity_requires_age_source_and_pulse_route():
    record = {"demography_kind": "introgression"}
    expectation = simulation._neutral_focal_expectation(record)

    valid = simulation._validate_neutral_focal_identity(
        _introgressed_neutral_test_tree(), expectation
    )
    assert valid["origin_time_generations"] == 2_400
    assert valid["observed_source_population"] == "Neanderthal"
    assert valid["pulse_route_evidence"]["n_descendant_lineages"] == 1
    assert valid["archaic_specific_no_ils_by_construction"] is True

    with pytest.raises(ValueError, match="age differs"):
        simulation._validate_neutral_focal_identity(
            _introgressed_neutral_test_tree(mutation_time=2_500), expectation
        )
    with pytest.raises(ValueError, match="not resident"):
        simulation._validate_neutral_focal_identity(
            _introgressed_neutral_test_tree(mutation_population="Loschbour"),
            expectation,
        )
    with pytest.raises(ValueError, match="pulse-route evidence"):
        simulation._validate_neutral_focal_identity(
            _introgressed_neutral_test_tree(include_pulse=False), expectation
        )


def test_exact_panel_preserves_migrations_but_only_chosen_sample_flags():
    pool = _introgressed_neutral_test_tree()
    original_migrations = [
        (
            migration.left,
            migration.right,
            migration.node,
            migration.source,
            migration.dest,
            migration.time,
        )
        for migration in pool.migrations()
    ]

    panel, panel_record = simulation.select_exact_sample_panel(
        pool,
        sample_diploids=6,
        exact_alt_count=6,
        seed=19,
    )
    manifest = simulation.build_diploid_pair_table(panel, 5_000_000)

    assert panel.num_samples == 12
    assert len(manifest) == 6
    assert int(manifest["focal_selected_allele_count"].sum()) == 6
    assert panel.num_nodes == pool.num_nodes
    assert panel.num_individuals == pool.num_individuals
    assert panel.num_populations == pool.num_populations
    assert [
        (
            migration.left,
            migration.right,
            migration.node,
            migration.source,
            migration.dest,
            migration.time,
        )
        for migration in panel.migrations()
    ] == original_migrations
    assert panel_record["candidate_pool_diploids"] == 8


def test_simulation_rejects_onedrive_output(repo_root):
    unit = _unit(repo_root, "eas_phlash_median", 0.1)
    with pytest.raises(ValueError, match="OneDrive"):
        simulation.simulate_unit(
            unit,
            repo_root,
            repo_root / "OneDrive - Institution" / "campaign",
        )


def test_selected_defaults_make_calibrated_acceptance_exhaustion_negligible():
    assert simulation.DEFAULT_SLIM_SCALING_FACTOR == 5.0
    assert simulation.DEFAULT_SELECTED_MAX_DRAWS == 50
    assert simulation.DEFAULT_SELECTED_CUMULATIVE_TIMEOUT_SECONDS == 45 * 60
    assert (1 - 0.20) ** simulation.DEFAULT_SELECTED_MAX_DRAWS < 1.5e-5


def test_introgression_selected_objects_use_explicit_calibration_endpoint(repo_root):
    plan = campaign.FocusedCampaignPlan(
        repo_root=repo_root,
        campaign_dir=repo_root / campaign.DEFAULT_CAMPAIGN_DIR,
        selected_replicates=1,
        neutral_replicates=1,
    )
    units = campaign.build_execution_units(plan)
    unit = units[
        (units["demography_id"] == "ancient_eurasia_han_introgression")
        & (units["target_allele_frequency"] == 0.2)
        & (units["simulation_class"] == "selected")
        & (units["selection_coefficient"] == 0.01)
    ].iloc[0]

    _, sweep, samples, model_record = simulation._selected_objects(
        unit.to_dict(),
        repo_root,
        5.0,
        500,
        han_selection_end_generations_ago=1_500.0,
    )
    fitness = [
        event
        for event in sweep.extended_events
        if isinstance(event, stdpopsim.ChangeMutationFitness)
    ]
    draws = [
        event
        for event in sweep.extended_events
        if isinstance(event, stdpopsim.DrawMutation)
    ]
    conditions = [
        event
        for event in sweep.extended_events
        if isinstance(event, stdpopsim.ConditionOnAlleleFrequency)
    ]
    by_population = {str(event.population): event for event in fitness}
    assert samples == {"Han": 500}
    assert set(by_population) == {"Loschbour", "Han"}
    assert float(by_population["Loschbour"].start_time) == 2_272
    assert float(by_population["Loschbour"].end_time) == 2_015
    assert float(by_population["Han"].start_time) == 2_015
    assert float(by_population["Han"].end_time) == 1_500
    assert len(draws) == 1
    assert str(draws[0].population) == "Neanderthal"
    assert float(draws[0].time) == 2_400
    assert len(conditions) == 4
    assert {str(event.population) for event in conditions} == {
        "Neanderthal",
        "Loschbour",
        "Han",
    }
    assert {str(event.op) for event in conditions} == {">", ">="}
    assert not any(
        str(event.population) == "Han"
        and float(event.start_time) == 0
        and float(event.end_time) == 0
        for event in conditions
    )
    assert model_record["sweep"]["mutation_origin"] == {
        "population": "Neanderthal",
        "age_generations": 2_400.0,
        "age_years": 60_000.0,
        "after_human_neanderthal_split": True,
        "before_introgression_pulse": True,
        "archaic_specific_no_ils_by_construction": True,
    }
    assert model_record["sweep"]["selection"]["mode"].endswith("finite_episode")
    assert model_record["selection_endpoint"]["source"] == (
        "explicit_calibration_candidate"
    )
    assert (
        model_record["selection_endpoint"]["binding"]["selection_end_generations_ago"]
        == 1_500
    )
    nonloss = model_record["sweep"]["post_pulse_recipient_nonloss"]
    assert nonloss["conditioned_on_introgressed_copy_reaching_recipient"] is True
    endpoint = model_record["sweep"]["focused_sample_endpoint_ascertainment"]
    assert endpoint["population_frequency_conditioning_in_slim"] is False
    assert endpoint["target_sample_allele_frequency"] == 0.2
    assert endpoint["sample_diploids"] == 100
    assert endpoint["exact_sample_alt_count"] == 40
    assert endpoint["candidate_pool_diploids"] == 500
    endpoint_gate = endpoint["candidate_pool_frequency_gate"]
    assert endpoint_gate["schema"] == simulation.CANDIDATE_POOL_FREQUENCY_GATE_SCHEMA
    assert endpoint_gate["lower_inclusive"] == pytest.approx(0.175)
    assert endpoint_gate["upper_inclusive"] == pytest.approx(0.225)
    assert endpoint_gate["applied_after_completed_slim_trajectory"] is True
    assert endpoint_gate["applied_before_exact_sample_panel"] is True
    assert endpoint["minimum_hom_alt_diploids"] == 2
    present = model_record["sweep"]["present_han_frequency"]
    assert present["population_conditioning_in_slim"] is False
    assert present["removed_terminal_condition_count"] == 2
    assert len(model_record["sweep"]["extended_events"]) == len(sweep.extended_events)


def test_eas_selected_objects_remove_terminal_frequency_conditions(repo_root):
    plan = campaign.FocusedCampaignPlan(
        repo_root=repo_root,
        campaign_dir=repo_root / campaign.DEFAULT_CAMPAIGN_DIR,
        selected_replicates=1,
        neutral_replicates=1,
    )
    units = campaign.build_execution_units(plan)
    unit = units[
        (units["demography_id"] == "eas_phlash_median")
        & (units["target_allele_frequency"] == 0.2)
        & (units["simulation_class"] == "selected")
        & (units["selection_coefficient"] == 0.01)
    ].iloc[0]
    _, sweep, samples, model_record = simulation._selected_objects(
        unit.to_dict(), repo_root, 5.0, 500
    )
    conditions = [
        event
        for event in sweep.extended_events
        if isinstance(event, stdpopsim.ConditionOnAlleleFrequency)
    ]
    assert samples == {"EAS": 500}
    assert len(conditions) == 1
    assert str(conditions[0].population) == "EAS"
    assert str(conditions[0].op) == ">"
    assert not any(
        math.isclose(float(event.start_time), 0.0, abs_tol=1e-12)
        and math.isclose(float(event.end_time), 0.0, abs_tol=1e-12)
        for event in conditions
    )
    present = model_record["sweep"]["present_frequency_interval"]
    assert present["population_conditioning_in_slim"] is False
    assert present["removed_terminal_condition_count"] == 2


def test_unit_lock_blocks_a_second_owner_and_recovers_without_deletion(tmp_path):
    lock_path = tmp_path / "unit" / ".simulation.lock"
    with simulation.exclusive_unit_lock(
        lock_path, unit_id="test-unit", phase="simulation"
    ):
        assert simulation.unit_lock_is_held(lock_path)
        with (
            pytest.raises(simulation.UnitLockHeld, match="already held"),
            simulation.exclusive_unit_lock(
                lock_path, unit_id="test-unit", phase="simulation"
            ),
        ):
            pass

    assert lock_path.is_file()
    assert not simulation.unit_lock_is_held(lock_path)
    with simulation.exclusive_unit_lock(
        lock_path, unit_id="test-unit", phase="simulation"
    ):
        assert lock_path.is_file()


@pytest.mark.parametrize(
    ("mutation_type", "expected_rejection"),
    [
        (None, "selected_focal_mutation_absent"),
        (2, "selected_focal_mutation_absent_background_only"),
    ],
)
def test_selected_draw_loss_is_recorded_and_retried(
    monkeypatch, tmp_path, mutation_type, expected_rejection
):
    plan = campaign.FocusedCampaignPlan(
        repo_root=tmp_path,
        campaign_dir=tmp_path / campaign.DEFAULT_CAMPAIGN_DIR,
        selected_replicates=1,
        neutral_replicates=1,
    )
    unit = (
        campaign.build_execution_units(plan)
        .query(
            "demography_id == 'ancient_eurasia_han_introgression' "
            "and simulation_class == 'selected' "
            "and selection_coefficient == 0.01 "
            "and target_allele_frequency == 0.1"
        )
        .iloc[0]
        .to_dict()
    )

    ancestry = _slim_focal_test_tree(mutation_type=mutation_type)

    class Engine:
        def simulate(self, *args, **kwargs):
            return ancestry

    monkeypatch.setattr(simulation.stdpopsim, "get_engine", lambda _: Engine())
    monkeypatch.setattr(simulation, "_wall_clock_timeout", lambda _: nullcontext())
    monkeypatch.setattr(
        simulation,
        "_selected_focal_expectation",
        lambda *args, **kwargs: _test_selected_focal_expectation(),
    )
    monkeypatch.setattr(
        simulation,
        "_selected_objects",
        lambda *args, **kwargs: (
            object(),
            SimpleNamespace(contig=object(), extended_events=()),
            {},
            {},
        ),
    )

    with pytest.raises(simulation.SimulationExhausted, match="no exact sampled AF"):
        simulation.simulate_selected_unit(
            unit,
            tmp_path,
            tmp_path / "unit",
            slim_path=tmp_path / "slim",
            max_draws=2,
            draw_timeout_seconds=10,
            cumulative_timeout_seconds=30,
        )

    ledger = simulation.pd.read_csv(tmp_path / "unit" / "attempts.tsv", sep="\t")
    assert len(ledger) == 2
    assert set(ledger["status"]) == {"complete"}
    assert set(ledger["candidate_pool_alt_count"]) == {0}
    assert set(ledger["candidate_pool_af"]) == {0.0}
    assert set(ledger["accepted"]) == {False}
    assert set(ledger["rejection_reason"]) == {expected_rejection}
    assert ledger["seed"].tolist() == [
        simulation._stable_seed(int(unit["seed"]), f"selected:{attempt}")
        for attempt in range(2)
    ]


def test_selected_candidate_pool_must_pass_band_before_panel_sampling(
    monkeypatch, tmp_path
):
    plan = campaign.FocusedCampaignPlan(
        repo_root=tmp_path,
        campaign_dir=tmp_path / campaign.DEFAULT_CAMPAIGN_DIR,
        selected_replicates=1,
        neutral_replicates=1,
    )
    unit = (
        campaign.build_execution_units(plan)
        .query(
            "demography_id == 'eas_phlash_median' "
            "and simulation_class == 'selected' "
            "and selection_coefficient == 0.01 "
            "and target_allele_frequency == 0.1"
        )
        .iloc[0]
        .to_dict()
    )
    ancestry = _slim_focal_test_tree(mutation_type=1, sample_diploids=10)

    class Engine:
        def simulate(self, *args, **kwargs):
            return ancestry

    monkeypatch.setattr(simulation.stdpopsim, "get_engine", lambda _: Engine())
    monkeypatch.setattr(simulation, "_wall_clock_timeout", lambda _: nullcontext())
    monkeypatch.setattr(
        simulation,
        "_selected_focal_expectation",
        lambda *args, **kwargs: _test_selected_focal_expectation(),
    )
    monkeypatch.setattr(
        simulation,
        "_selected_objects",
        lambda *args, **kwargs: (
            object(),
            SimpleNamespace(contig=object(), extended_events=()),
            {},
            {},
        ),
    )

    def panel_must_not_run(*args, **kwargs):
        pytest.fail("exact-panel sampler ran before the candidate-pool AF gate")

    monkeypatch.setattr(simulation, "select_exact_sample_panel", panel_must_not_run)
    with pytest.raises(simulation.SimulationExhausted, match="no exact sampled AF"):
        simulation.simulate_selected_unit(
            unit,
            tmp_path,
            tmp_path / "selected-outside-pool-band",
            slim_path=tmp_path / "slim",
            selected_pool_diploids=10,
            max_draws=1,
            draw_timeout_seconds=10,
            cumulative_timeout_seconds=30,
        )

    ledger = simulation.pd.read_csv(
        tmp_path / "selected-outside-pool-band" / "attempts.tsv", sep="\t"
    )
    assert ledger["candidate_pool_alt_count"].tolist() == [1]
    assert ledger["candidate_pool_af"].tolist() == [0.05]
    assert ledger["candidate_pool_frequency_gate_passed"].tolist() == [False]
    assert ledger["accepted"].tolist() == [False]
    assert ledger["rejection_reason"].tolist() == [
        "candidate_pool_af_outside_prespecified_band"
    ]


def test_selected_focal_identity_requires_type_source_time_and_metadata():
    expectation = _test_selected_focal_expectation()
    valid = simulation.validate_selected_focal_identity(
        _slim_focal_test_tree(mutation_type=1), expectation
    )
    assert valid["status"] == "valid"
    assert valid["observed_mutation_type"] == 1
    assert valid["observed_source_subpopulation"] == 0
    assert valid["observed_origin_time_generations"] == 10

    background = simulation.inspect_selected_focal_identity(
        _slim_focal_test_tree(mutation_type=2), expectation
    )
    assert background["status"] == "selected_focal_mutation_absent_background_only"
    assert background["observed_mutation_types"] == [2]
    with pytest.raises(simulation.SelectedFocalIdentityError, match="background-only"):
        simulation.validate_selected_focal_identity(
            _slim_focal_test_tree(mutation_type=2), expectation
        )

    with pytest.raises(simulation.SelectedFocalIdentityError, match="subpopulation"):
        simulation.validate_selected_focal_identity(
            _slim_focal_test_tree(mutation_type=1, subpopulation=7), expectation
        )


def test_semantic_validation_requires_exact_tree_manifest_and_af():
    ts = _neutral_semantic_test_tree()
    pairs = simulation.build_diploid_pair_table(ts, 5_000_000)
    record = {
        "simulation_class": "neutral",
        "sample_diploids": 6,
        "exact_sample_alt_count": 5,
    }
    contract = {
        "schema": simulation.FOCAL_SEMANTIC_VALIDATION_SCHEMA,
        "simulation_class": "neutral",
        "sample_diploids": 6,
        "exact_sample_alt_count": 5,
        "minimum_hom_ref_diploids": 2,
        "minimum_heterozygous_diploids": 1,
        "minimum_hom_alt_diploids": 2,
        "identity": {
            "schema": "gamma-smc.neutral-focal-identity/v1",
            "focal_position_0based": 5_000_000,
            "ancestral_state": "A",
            "derived_state": "G",
            "single_nonrecurrent_mutation": True,
        },
        "tree_manifest_exact_match": True,
    }

    result = simulation._validate_unit_semantics(ts, pairs, record, contract)
    assert result["sample_alt_count"] == 5
    assert result["genotype_counts"] == {
        "hom_alt": 2,
        "hom_ref": 3,
        "heterozygous": 1,
    }

    wrong_manifest = pairs.copy()
    wrong_manifest.loc[0, "sample_node_0"] = 999
    with pytest.raises(ValueError, match="does not exactly match"):
        simulation._validate_unit_semantics(ts, wrong_manifest, record, contract)

    with pytest.raises(ValueError, match="exact AF"):
        simulation._validate_unit_semantics(
            ts,
            pairs,
            {**record, "exact_sample_alt_count": 6},
            {**contract, "exact_sample_alt_count": 6},
        )


def test_completion_cache_revalidates_checksum_consistent_manifest(tmp_path):
    ts = _neutral_semantic_test_tree()
    pair_table = simulation.build_diploid_pair_table(ts, 5_000_000)
    record = {
        "unit_id": "synthetic-neutral",
        "demography_id": "eas_phlash_median",
        "demography_kind": "no_introgression",
        "population": "EAS",
        "simulation_class": "neutral",
        "selection_coefficient": 0.0,
        "target_allele_frequency": 5 / 12,
        "population_af_lower": 0.39,
        "population_af_upper": 0.44,
        "sample_diploids": 6,
        "exact_sample_alt_count": 5,
        "replicate_index": 0,
        "seed": 1,
        "sequence_length_bp": 10_000_000,
        "focal_position_bp": 5_000_000,
    }
    validation_contract = {
        "schema": simulation.FOCAL_SEMANTIC_VALIDATION_SCHEMA,
        "simulation_class": "neutral",
        "sample_diploids": 6,
        "exact_sample_alt_count": 5,
        "minimum_hom_ref_diploids": 2,
        "minimum_heterozygous_diploids": 1,
        "minimum_hom_alt_diploids": 2,
        "identity": {
            "schema": "gamma-smc.neutral-focal-identity/v1",
            "focal_position_0based": 5_000_000,
            "ancestral_state": "A",
            "derived_state": "G",
            "single_nonrecurrent_mutation": True,
        },
        "tree_manifest_exact_match": True,
    }
    contract = {
        "schema": simulation.SCHEMA_VERSION,
        "unit": record,
        "parameters": {
            "candidate_pool_diploids": 500,
            "focal_semantic_validation": validation_contract,
        },
    }
    artifacts = simulation._completion_artifacts(tmp_path)
    artifacts.tree_path.parent.mkdir(parents=True, exist_ok=True)
    ts.dump(artifacts.tree_path)
    pair_table.to_csv(artifacts.pair_table_path, sep="\t", index=False)
    for path in (
        artifacts.overall_pairs_path,
        artifacts.truth_profiles_path,
        artifacts.truth_class_summaries_path,
    ):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"test")
    semantic = simulation._validate_unit_semantics(
        ts, pair_table, record, validation_contract
    )
    outputs = {
        label: {
            "path": path.relative_to(tmp_path).as_posix(),
            "sha256": simulation.sha256_file(path),
        }
        for label, path in {
            "tree": artifacts.tree_path,
            "pair_table": artifacts.pair_table_path,
            "overall_pairs": artifacts.overall_pairs_path,
            "truth_profiles": artifacts.truth_profiles_path,
            "truth_class_summaries": artifacts.truth_class_summaries_path,
        }.items()
    }
    pool_counts = np.asarray([0] * 85 + [1] * 413 + [2] * 2, dtype=int)
    completion = {
        "schema": simulation.SCHEMA_VERSION,
        "status": "complete",
        "contract": contract,
        "contract_sha256": simulation._canonical_sha256(contract),
        "provenance": {
            "candidate_pool_frequency_gate": (
                simulation.candidate_pool_frequency_gate(record, pool_counts)
            ),
            "sample_panel": {
                "candidate_pool_diploids": 500,
                "candidate_pool_genotype_counts": {
                    "hom_ref": 85,
                    "heterozygous": 413,
                    "hom_alt": 2,
                },
                "eligible_genotype_count_triples": 1,
                "chosen_sample_genotype_counts": {
                    "hom_ref": 3,
                    "heterozygous": 1,
                    "hom_alt": 2,
                },
                "selection_method": (
                    "uniform_over_all_size_N_exact_alt_count_subsets_passing_"
                    "genotype_qc"
                ),
                "panel_seed": 1,
            },
            "focal_semantic_validation": semantic,
            "neutral_focal_identity": semantic["identity"],
        },
        "outputs": outputs,
    }
    (tmp_path / "simulation_contract.json").write_text(
        json.dumps(contract), encoding="utf-8"
    )
    artifacts.completion_path.write_text(json.dumps(completion), encoding="utf-8")
    assert simulation._validate_completion(artifacts, contract).cache_hit is True

    invalid_gate_completion = copy.deepcopy(completion)
    invalid_gate_completion["provenance"]["candidate_pool_frequency_gate"][
        "observed_alt_count"
    ] = 389
    invalid_gate_completion["provenance"]["candidate_pool_frequency_gate"][
        "observed_af"
    ] = 0.389
    artifacts.completion_path.write_text(
        json.dumps(invalid_gate_completion), encoding="utf-8"
    )
    with pytest.raises(ValueError, match="observed count is outside bounds"):
        simulation._validate_completion(artifacts, contract)

    pair_table.loc[0, "sample_node_0"] = 999
    pair_table.to_csv(artifacts.pair_table_path, sep="\t", index=False)
    completion["outputs"]["pair_table"]["sha256"] = simulation.sha256_file(
        artifacts.pair_table_path
    )
    artifacts.completion_path.write_text(json.dumps(completion), encoding="utf-8")
    with pytest.raises(ValueError, match="does not exactly match"):
        simulation._validate_completion(artifacts, contract)


def test_completion_cache_rejects_old_contract_before_artifact_reuse(tmp_path):
    artifacts = simulation._completion_artifacts(tmp_path)
    old_contract = {"schema": simulation.SCHEMA_VERSION, "unit": {}}
    current_contract = {
        "schema": simulation.SCHEMA_VERSION,
        "unit": {},
        "parameters": {"focal_semantic_validation": {}},
    }
    artifacts.completion_path.parent.mkdir(parents=True, exist_ok=True)
    artifacts.completion_path.write_text(
        json.dumps(
            {
                "schema": simulation.SCHEMA_VERSION,
                "status": "complete",
                "contract": old_contract,
                "contract_sha256": simulation._canonical_sha256(old_contract),
            }
        ),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="contract is incompatible"):
        simulation._validate_completion(artifacts, current_contract)


@pytest.fixture
def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]
