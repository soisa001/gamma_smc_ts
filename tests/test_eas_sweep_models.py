from __future__ import annotations

import hashlib
import json
import math

import msprime
import numpy as np
import pytest
import stdpopsim

from gamma_smc_aou.eas_sweep_models import (
    ARCHAIC_POPULATION,
    DE_NOVO_DEFAULT_LOWER_FREQUENCY,
    DE_NOVO_DEFAULT_UPPER_FREQUENCY,
    DEFAULT_ARCHAIC_SOURCE_MINIMUM_FREQUENCY,
    FOCAL_POSITION_BP,
    FOCAL_SITE_ID,
    GENERATION_TIME_YEARS,
    HAN_SPLIT_GENERATIONS,
    HUMAN_NEANDERTHAL_SPLIT_GENERATIONS,
    INTROGRESSION_PULSE_GENERATIONS,
    INTROGRESSION_RECIPIENT_POPULATION,
    INTROGRESSION_TARGET_POPULATION,
    MUTATION_RATE,
    PHLASH_MVN_SCHEMA,
    PRE_PULSE_SOURCE_CHECK_GENERATIONS,
    RECOMBINATION_RATE,
    SELECTION_COEFFICIENTS,
    SEQUENCE_LENGTH_BP,
    TMRCA_THRESHOLDS_GENERATIONS,
    TMRCA_THRESHOLDS_YEARS,
    build_de_novo_origin_age_grid,
    build_eas_demography_models,
    build_introgression_sweep_spec,
    build_no_introgression_sweep_spec,
    estimate_de_novo_origin_age,
    load_ancient_eurasia_model,
    load_phlash_eas_npz,
)


def _write_eas_npz(tmp_path, **overrides):
    time = np.array([100.0, 1_000.0, 10_000.0, 40_000.0])
    bootstrap_ne = np.array(
        [
            [8_000.0, 9_000.0, 10_000.0, 11_000.0],
            [10_000.0, 11_000.0, 12_000.0, 13_000.0],
            [12_000.0, 13_000.0, 14_000.0, 15_000.0],
            [14_000.0, 15_000.0, 16_000.0, 17_000.0],
        ]
    )
    payload = {
        "schema": np.asarray(PHLASH_MVN_SCHEMA),
        "population": np.asarray("EAS"),
        "time": time,
        "mean_log_ne": np.log(np.median(bootstrap_ne, axis=0)),
        "covariance_factor": np.full((2, time.size), 0.01),
        "bootstrap_ne": bootstrap_ne,
        "jitter": np.asarray(0.001),
    }
    payload.update(overrides)
    path = tmp_path / "EAS.npz"
    np.savez_compressed(path, **payload)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return path, digest, bootstrap_ne


def _events_of_type(spec, event_type):
    return [event for event in spec.extended_events if isinstance(event, event_type)]


def test_constants_match_study_contract():
    assert TMRCA_THRESHOLDS_YEARS == (
        1_000,
        4_500,
        10_000,
        20_000,
        30_000,
        40_000,
        50_000,
    )
    assert TMRCA_THRESHOLDS_GENERATIONS == (40, 180, 400, 800, 1_200, 1_600, 2_000)
    assert SELECTION_COEFFICIENTS == (0.01, 0.005, 0.001)
    assert GENERATION_TIME_YEARS == 25
    assert SEQUENCE_LENGTH_BP == 10_000_000
    assert MUTATION_RATE == 1.25e-8
    assert RECOMBINATION_RATE == 1e-8


def test_load_eas_artifact_verifies_identity_arrays_and_pointwise_quantiles(tmp_path):
    path, digest, bootstrap = _write_eas_npz(tmp_path)

    artifact = load_phlash_eas_npz(path, expected_sha256=digest.upper())

    expected = np.quantile(bootstrap, [0.025, 0.5, 0.975], axis=0)
    np.testing.assert_allclose(artifact.q025_ne, expected[0])
    np.testing.assert_allclose(artifact.median_ne, expected[1])
    np.testing.assert_allclose(artifact.q975_ne, expected[2])
    assert artifact.schema == PHLASH_MVN_SCHEMA
    assert artifact.population == "EAS"
    assert artifact.actual_sha256 == digest
    assert artifact.expected_sha256 == digest
    record = artifact.provenance_record()
    assert record["source_sha256_verified"] is True
    assert record["number_of_bootstrap_fits"] == 4
    json.dumps(record)


def test_load_eas_artifact_rejects_hash_schema_population_and_bad_arrays(tmp_path):
    path, digest, _ = _write_eas_npz(tmp_path)
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        load_phlash_eas_npz(path, expected_sha256="0" * 64)

    path, digest, _ = _write_eas_npz(tmp_path, schema=np.asarray("wrong/v1"))
    with pytest.raises(ValueError, match="artifact schema"):
        load_phlash_eas_npz(path, expected_sha256=digest)

    path, digest, _ = _write_eas_npz(tmp_path, population=np.asarray("EUR"))
    with pytest.raises(ValueError, match="artifact population"):
        load_phlash_eas_npz(path, expected_sha256=digest)

    path, digest, _ = _write_eas_npz(
        tmp_path,
        bootstrap_ne=np.array([[8_000.0, np.nan, 10_000.0, 11_000.0]] * 2),
    )
    with pytest.raises(ValueError, match="bootstrap_ne must be positive and finite"):
        load_phlash_eas_npz(path, expected_sha256=digest)


def test_demographic_models_use_quantiles_and_explicit_zero_to_100_extrapolation(
    tmp_path,
):
    path, digest, bootstrap = _write_eas_npz(tmp_path)
    artifact = load_phlash_eas_npz(path, expected_sha256=digest)

    models = build_eas_demography_models(artifact)

    assert list(models) == ["q025", "median", "q975"]
    expected = np.quantile(bootstrap, [0.025, 0.5, 0.975], axis=0)
    for index, label in enumerate(models):
        bundle = models[label]
        assert bundle.time_generations.tolist() == [0, 100, 1_000, 10_000, 40_000]
        assert bundle.ne[0] == pytest.approx(bundle.ne[1])
        np.testing.assert_allclose(bundle.ne[1:], expected[index])
        assert bundle.msprime_demography.populations[0].name == "EAS"
        assert [event.time for event in bundle.msprime_demography.events] == [
            100,
            1_000,
            10_000,
            40_000,
        ]
        assert bundle.stdpopsim_model.generation_time == 25
        assert bundle.stdpopsim_model.mutation_rate == MUTATION_RATE
        assert bundle.stdpopsim_model.recombination_rate == RECOMBINATION_RATE
        extrapolation = bundle.record["presentward_extrapolation"]
        assert extrapolation["from_generations_ago"] == 100
        assert extrapolation["to_generations_ago"] == 0
        json.dumps(bundle.provenance_record())


def test_origin_age_uses_h_times_s_logit_approximation_and_builds_full_grid(tmp_path):
    times = np.array([0.0, 100.0, 40_000.0])
    sizes = np.array([10_000.0, 10_000.0, 10_000.0])

    estimate = estimate_de_novo_origin_age(
        times,
        sizes,
        selection_coefficient=0.01,
        trajectory_label="median",
    )

    expected_age = (
        math.log(0.5 / (1 - 0.5)) - math.log((1 / 20_000) / (1 - 1 / 20_000))
    ) / (0.5 * 0.01)
    assert estimate.age_generations == pytest.approx(expected_age, abs=0.51)
    assert estimate.age_years == estimate.age_generations * 25
    assert estimate.initial_frequency == pytest.approx(1 / 20_000)
    assert estimate.slim_scaling_factor == 1
    assert estimate.to_record()["is_exact"] is False
    assert "h*s" in estimate.to_record()["equation"]

    scaled_estimate = estimate_de_novo_origin_age(
        times,
        sizes,
        selection_coefficient=0.01,
        trajectory_label="median",
        slim_scaling_factor=10,
    )
    scaled_expected_age = (
        math.log(0.5 / (1 - 0.5)) - math.log((10 / 20_000) / (1 - 10 / 20_000))
    ) / (0.5 * 0.01)
    assert scaled_estimate.age_generations == pytest.approx(
        scaled_expected_age, abs=0.51
    )
    assert scaled_estimate.age_generations < estimate.age_generations
    assert scaled_estimate.initial_frequency == pytest.approx(10 / 20_000)
    scaled_record = scaled_estimate.to_record()
    assert scaled_record["schema"] == "gamma-smc.de-novo-origin-age/v2"
    assert scaled_record["slim_scaling_factor"] == 10
    assert "Q/(2*Ne" in scaled_record["equation"]

    path, digest, _ = _write_eas_npz(tmp_path)
    models = build_eas_demography_models(
        load_phlash_eas_npz(path, expected_sha256=digest)
    )
    grid = build_de_novo_origin_age_grid(models, slim_scaling_factor=10)
    assert len(grid) == 3 * 3
    assert [(x.selection_coefficient, x.trajectory_label) for x in grid] == [
        (selection, label)
        for selection in SELECTION_COEFFICIENTS
        for label in ("q025", "median", "q975")
    ]
    assert all(item.to_record()["target_frequency_is_placeholder"] for item in grid)
    assert all(item.slim_scaling_factor == 10 for item in grid)


def test_no_introgression_spec_draws_once_selects_from_origin_and_bounds_present_af():
    spec = build_no_introgression_sweep_spec(
        origin_age_generations=2_000,
        selection_coefficient=0.01,
    )

    draws = _events_of_type(spec, stdpopsim.DrawMutation)
    fitness = _events_of_type(spec, stdpopsim.ChangeMutationFitness)
    conditions = _events_of_type(spec, stdpopsim.ConditionOnAlleleFrequency)
    assert len(draws) == 1
    assert draws[0].time == 2_000
    assert draws[0].population == "EAS"
    assert len(fitness) == 1
    assert (fitness[0].start_time, fitness[0].end_time) == (2_000, 0)
    assert fitness[0].population == "EAS"
    present = [condition for condition in conditions if condition.start_time == 0]
    assert [(event.op, event.allele_frequency) for event in present] == [
        (">=", DE_NOVO_DEFAULT_LOWER_FREQUENCY),
        ("<=", DE_NOVO_DEFAULT_UPPER_FREQUENCY),
    ]

    assert spec.contig.length == SEQUENCE_LENGTH_BP
    assert spec.contig.mutation_rate == MUTATION_RATE
    assert spec.contig.recombination_map.mean_rate == RECOMBINATION_RATE
    focal_index = [
        i for i, dfe in enumerate(spec.contig.dfe_list) if dfe.id == FOCAL_SITE_ID
    ]
    assert len(focal_index) == 1
    np.testing.assert_array_equal(
        spec.contig.interval_list[focal_index[0]],
        [[FOCAL_POSITION_BP, FOCAL_POSITION_BP + 1]],
    )
    assert spec.record["contig"]["recurrence_at_focal_base"] is False
    assert spec.record["target_frequency_is_placeholder"] is True
    json.dumps(spec.provenance_record())


def test_ancient_eurasia_contract_and_introgression_event_timing():
    model = load_ancient_eurasia_model()
    assert model.id == "AncientEurasia_9K19"
    assert model.generation_time == 25
    assert {population.name for population in model.populations}.issuperset(
        {
            ARCHAIC_POPULATION,
            INTROGRESSION_RECIPIENT_POPULATION,
            INTROGRESSION_TARGET_POPULATION,
        }
    )
    mass_migration_times = {
        event.time
        for event in model.model.events
        if isinstance(event, msprime.MassMigration)
    }
    assert {
        HAN_SPLIT_GENERATIONS,
        INTROGRESSION_PULSE_GENERATIONS,
        HUMAN_NEANDERTHAL_SPLIT_GENERATIONS,
    }.issubset(mass_migration_times)

    spec = build_introgression_sweep_spec(
        mutation_age_generations=10_000,
        selection_coefficient=0.005,
        source_minimum_frequency=1.0,
        minimum_han_frequency=0.2,
        maximum_han_frequency=0.3,
    )
    draws = _events_of_type(spec, stdpopsim.DrawMutation)
    fitness = _events_of_type(spec, stdpopsim.ChangeMutationFitness)
    conditions = _events_of_type(spec, stdpopsim.ConditionOnAlleleFrequency)
    assert len(draws) == 1
    assert draws[0].population == ARCHAIC_POPULATION
    assert draws[0].time == 10_000
    assert all(event.population != ARCHAIC_POPULATION for event in fitness)
    assert [
        (event.population, event.start_time, event.end_time) for event in fitness
    ] == [
        (
            INTROGRESSION_RECIPIENT_POPULATION,
            INTROGRESSION_PULSE_GENERATIONS,
            HAN_SPLIT_GENERATIONS,
        ),
        (INTROGRESSION_TARGET_POPULATION, HAN_SPLIT_GENERATIONS, 0),
    ]
    assert spec.record["selection"]["requested_han_split_generations_ago"] == 2_016
    assert spec.record["selection"]["realized_han_split_generations_ago"] == 2_016

    source_check = [
        event
        for event in conditions
        if event.population == ARCHAIC_POPULATION
        and event.start_time == PRE_PULSE_SOURCE_CHECK_GENERATIONS
        and event.end_time == PRE_PULSE_SOURCE_CHECK_GENERATIONS
    ]
    assert len(source_check) == 1
    assert (source_check[0].op, source_check[0].allele_frequency) == (">=", 1.0)
    # Conditions spanning the Loschbour/Han restore and split are unsafe in
    # stdpopsim 0.3.0's generated SLiM: those subpopulations are not defined
    # for the entire callback interval. Only fitness callbacks span the split.
    assert not any(
        event.population == INTROGRESSION_RECIPIENT_POPULATION for event in conditions
    )
    assert not any(
        event.population == INTROGRESSION_TARGET_POPULATION and event.start_time != 0
        for event in conditions
    )
    present_han = [
        event
        for event in conditions
        if event.population == INTROGRESSION_TARGET_POPULATION and event.start_time == 0
    ]
    assert [(event.op, event.allele_frequency) for event in present_han] == [
        (">=", 0.2),
        ("<=", 0.3),
    ]
    assert spec.record["mutation_origin"]["archaic_specific_no_ils_by_construction"]
    assert spec.record["selection"]["mode"] == (
        "standing_variation_at_introgression_onset"
    )
    json.dumps(spec.provenance_record())


def test_q10_generated_script_aligns_han_split_and_fitness_ticks(capsys):
    realized_split = math.floor(HAN_SPLIT_GENERATIONS / 10) * 10
    spec = build_introgression_sweep_spec(
        mutation_age_generations=2_400,
        selection_coefficient=0.005,
        source_minimum_frequency=0.1,
        minimum_han_frequency=0.5,
        maximum_han_frequency=0.95,
        realized_han_split_generations=realized_split,
    )

    with pytest.warns(stdpopsim.SLiMScalingFactorWarning):
        stdpopsim.get_engine("slim").simulate(
            spec.demographic_model,
            spec.contig,
            {INTROGRESSION_TARGET_POPULATION: 2},
            seed=123,
            extended_events=list(spec.extended_events),
            slim_script=True,
            slim_scaling_factor=10,
            slim_burn_in=0.1,
        )
    script = capsys.readouterr().out

    # The catalog split remains 2,016 generations (50,400 years) in _T and
    # reaches tick 201 through generated Eidos asInteger() truncation. Both
    # fitness callbacks use 2,010 generations (50,250 years), also tick 201.
    assert "51500, 50400, 50000" in script
    assert "c(_T[5], 5, _N[6,5], 3)" in script
    assert "c(56750, 50250, 1, 3, 0.005, 0.5)" in script
    assert "c(50250, 0, 1, 5, 0.005, 0.5)" in script
    demographic_tick = int((HAN_SPLIT_GENERATIONS * GENERATION_TIME_YEARS) / 25 / 10)
    fitness_tick = int((realized_split * GENERATION_TIME_YEARS) / 25 / 10)
    assert demographic_tick == fitness_tick == 201
    assert spec.record["selection"]["requested_han_split_generations_ago"] == 2_016
    assert spec.record["selection"]["realized_han_split_generations_ago"] == 2_010


def test_introgression_optional_upper_bound_and_origin_must_be_on_archaic_branch():
    spec = build_introgression_sweep_spec(
        selection_coefficient=0.001,
        minimum_han_frequency=0.9,
    )
    assert spec.record["source_frequency_condition"]["minimum_inclusive"] == (
        DEFAULT_ARCHAIC_SOURCE_MINIMUM_FREQUENCY
    )
    present_han = [
        event
        for event in _events_of_type(spec, stdpopsim.ConditionOnAlleleFrequency)
        if event.population == INTROGRESSION_TARGET_POPULATION and event.start_time == 0
    ]
    assert [(event.op, event.allele_frequency) for event in present_han] == [
        (">=", 0.9)
    ]

    with pytest.raises(ValueError, match="after the human-Neanderthal split"):
        build_introgression_sweep_spec(
            selection_coefficient=0.001,
            minimum_han_frequency=0.1,
            mutation_age_generations=30_000,
        )
