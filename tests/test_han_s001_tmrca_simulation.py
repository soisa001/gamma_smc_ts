import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import stdpopsim
import tskit

from gamma_smc_aou import han_s001_tmrca_simulation as study


def test_exact_tmrca_and_natural_class_contract() -> None:
    assert study.TMRCA_THRESHOLDS_YEARS == (
        1_000,
        4_500,
        10_000,
        20_000,
        30_000,
        40_000,
        50_000,
    )
    assert study.DEFAULT_SELECTION_COEFFICIENTS == {
        "neutral": 0.0,
        "selected": 0.001,
    }
    assert study.DEFAULT_POOL_DIPLOIDS == 500
    assert study.DEFAULT_PANEL_DIPLOIDS == 100
    assert study.WORK_RELATIVE_DIR == Path(
        "focused_selection_EAS_sim/work/han_s001_tmrca"
    )


def test_plan_is_idempotent_source_bound_and_has_fixed_units(
    monkeypatch, tmp_path
) -> None:
    monkeypatch.setattr(study, "_source_hashes", lambda _root: {"source.py": "a" * 64})
    monkeypatch.setattr(
        study,
        "_slim_patch_contract",
        lambda: {"stdpopsim_version": "0.3.0", "contract": "fixed"},
    )
    monkeypatch.setattr(study, "_runtime_versions", lambda: {"python": "test"})
    path = study.write_plan(tmp_path, n_neutral=2, n_selected=3)
    assert path == study.plan_path(tmp_path)
    payload = study.load_plan(tmp_path)
    assert (payload["n_neutral"], payload["n_selected"]) == (2, 3)
    assert [unit["unit_id"] for unit in payload["units"]] == [
        "neutral_0000",
        "neutral_0001",
        "selected_0000",
        "selected_0001",
        "selected_0002",
    ]
    assert len({unit["seed"] for unit in payload["units"]}) == 5
    assert study.write_plan(tmp_path, n_neutral=2, n_selected=3) == path
    with pytest.raises(ValueError, match="existing study plan differs"):
        study.write_plan(tmp_path, n_neutral=2, n_selected=4)

    tampered = json.loads(path.read_text(encoding="utf-8"))
    tampered["design"]["terminal_frequency_conditioning"] = True
    path.write_text(json.dumps(tampered), encoding="utf-8")
    with pytest.raises(ValueError, match="contract checksum failed"):
        study.load_plan(tmp_path)


def test_donor_fixed_patch_is_singleton_scoped_and_restores_exact_object() -> None:
    original = study.slim_engine._slim_functions
    with study._scoped_donor_fixed_patch(5):
        patched = study.slim_engine._slim_functions
        assert patched is not original
        assert patched.count(study._DONOR_FIXED_ADD_MUT) == 1
        assert "mut = targets.addNewDrawnMutation" in patched
        assert "remaining = setDifference(pop.genomes, targets)" in patched
        assert "remaining.addMutations(mut)" in patched
        assert "mut_type = m1" in patched
        assert "sim.setValue" not in study._DONOR_FIXED_ADD_MUT
    assert study.slim_engine._slim_functions is original

    with pytest.raises(RuntimeError, match="sentinel"):
        with study._scoped_donor_fixed_patch(5):
            raise RuntimeError("sentinel")
    assert study.slim_engine._slim_functions is original


def test_extended_events_condition_only_on_survival() -> None:
    neutral = study._build_extended_events(0.0)
    selected = study._build_extended_events(0.001)
    assert sum(isinstance(event, stdpopsim.DrawMutation) for event in neutral) == 1
    neutral_conditions = [
        event
        for event in neutral
        if isinstance(event, stdpopsim.ConditionOnAlleleFrequency)
    ]
    assert len(neutral_conditions) == 2
    assert {
        (event.population, event.op, event.allele_frequency)
        for event in neutral_conditions
    } == {
        ("Loschbour", ">", 0.0),
        ("Han", ">", 0.0),
    }
    assert not any(
        isinstance(event, stdpopsim.ChangeMutationFitness) for event in neutral
    )
    fitness = [
        event
        for event in selected
        if isinstance(event, stdpopsim.ChangeMutationFitness)
    ]
    assert len(fitness) == 2
    assert {event.population for event in fitness} == {"Loschbour", "Han"}
    assert {event.selection_coeff for event in fitness} == {0.001}
    assert {event.dominance_coeff for event in fitness} == {0.5}


def _site_free_pool() -> tskit.TreeSequence:
    tables = tskit.TableCollection(sequence_length=study.SEQUENCE_LENGTH_BP)
    population = tables.populations.add_row()
    root = tables.nodes.add_row(time=1, population=population)
    for _ in range(study.DEFAULT_POOL_DIPLOIDS):
        individual = tables.individuals.add_row()
        for _haplotype in range(2):
            node = tables.nodes.add_row(
                flags=tskit.NODE_IS_SAMPLE,
                time=0,
                population=population,
                individual=individual,
            )
            tables.edges.add_row(0, study.SEQUENCE_LENGTH_BP, root, node)
    tables.sort()
    return tables.tree_sequence()


def test_pool_and_panel_undetected_are_accepted() -> None:
    pool = _site_free_pool()
    counts = study._focal_counts_or_zero(pool)
    assert counts.shape == (study.DEFAULT_POOL_DIPLOIDS,)
    assert counts.sum() == 0
    panel, panel_counts, chosen = study._uniform_panel(pool, counts, seed=19)
    assert len(chosen) == study.DEFAULT_PANEL_DIPLOIDS
    assert panel_counts.sum() == 0
    assert panel.num_samples == 2 * study.DEFAULT_PANEL_DIPLOIDS
    manifest = study.build_diploid_pair_table(
        panel,
        study.FOCAL_POSITION_BP,
        focal_genotype_counts=panel_counts,
    )
    assert set(manifest["genotype_class"]) == {"hom_ref"}


def _strict_metadata() -> dict:
    values = {
        "donor_fixed_source_af_at_draw": [1.0],
        "donor_fixed_final_population_alt_count": [2],
        "donor_fixed_final_population_total_count": [100],
        "donor_fixed_final_population_af": [0.02],
        "donor_fixed_final_population_survived": [True],
        "donor_fixed_final_population_fixed": [False],
    }
    return {"SLiM": {"user_metadata": values}}


def test_slim_metadata_uses_one_exact_path_and_singletons() -> None:
    metadata = _strict_metadata()
    assert study._metadata_scalar(metadata, "donor_fixed_final_population_af") == 0.02
    assert (
        study._metadata_integer(metadata, "donor_fixed_final_population_alt_count") == 2
    )

    missing = _strict_metadata()
    del missing["SLiM"]["user_metadata"]["donor_fixed_final_population_fixed"]
    with pytest.raises(ValueError, match="lacks donor/final keys"):
        study._metadata_scalar(missing, "donor_fixed_final_population_af")

    duplicate_path = _strict_metadata()
    duplicate_path["donor_fixed_final_population_af"] = 0.02
    with pytest.raises(ValueError, match="occurs outside its one exact path"):
        study._metadata_scalar(duplicate_path, "donor_fixed_final_population_af")

    nonscalar = _strict_metadata()
    nonscalar["SLiM"]["user_metadata"]["donor_fixed_final_population_af"] = [0.02, 0.03]
    with pytest.raises(ValueError, match="not scalar"):
        study._metadata_scalar(nonscalar, "donor_fixed_final_population_af")


def _truth_profile() -> pd.DataFrame:
    rows = []
    for position in range(0, study.SEQUENCE_LENGTH_BP, study.DEFAULT_OUTPUT_STRIDE_BP):
        for threshold_index, threshold in enumerate(study.TMRCA_THRESHOLDS_YEARS):
            rows.append(
                {
                    "position_0based": position,
                    "position_1based": position + 1,
                    "threshold_years": threshold,
                    "threshold_generations": (threshold / study.GENERATION_TIME_YEARS),
                    "generation_time_years": study.GENERATION_TIME_YEARS,
                    "n_pairs": study.DEFAULT_PANEL_DIPLOIDS,
                    "p_tmrca_lt_threshold": 0.05 + 0.1 * threshold_index,
                    "genotype_class": "overall",
                }
            )
    return pd.DataFrame(rows)


def _gamma_profile() -> pd.DataFrame:
    frame = pd.DataFrame(
        {
            "position_0based": np.arange(
                0,
                study.SEQUENCE_LENGTH_BP,
                study.DEFAULT_OUTPUT_STRIDE_BP,
            ),
            "position_1based": np.arange(
                0,
                study.SEQUENCE_LENGTH_BP,
                study.DEFAULT_OUTPUT_STRIDE_BP,
            )
            + 1,
            "n_pairs": study.DEFAULT_PANEL_DIPLOIDS,
            "mean_tmrca_generations": 100.0,
        }
    )
    for threshold_index, threshold in enumerate(study.TMRCA_THRESHOLDS_YEARS):
        frame[f"mean_p_lt_{threshold}"] = 0.05 + 0.1 * threshold_index
    return frame


@pytest.mark.parametrize("mode", ("missing", "duplicate", "off_grid"))
@pytest.mark.parametrize("source", ("truth", "gamma"))
def test_score_reducers_reject_grid_tampering(source, mode) -> None:
    frame = _truth_profile() if source == "truth" else _gamma_profile()
    if mode == "missing":
        frame = frame.iloc[1:].copy()
    elif mode == "duplicate":
        frame = pd.concat([frame, frame.iloc[[0]]], ignore_index=True)
    else:
        frame = frame.copy()
        frame.loc[0, "position_0based"] = 1
    reducer = study._truth_scores if source == "truth" else study._gamma_scores
    with pytest.raises(ValueError, match="exactly one row"):
        reducer(frame)


def test_score_reducers_use_one_focal_and_21_local_positions() -> None:
    truth = study._truth_scores(_truth_profile())
    gamma = study._gamma_scores(_gamma_profile())
    assert len(truth) == len(gamma) == 2 * len(study.TMRCA_THRESHOLDS_YEARS)
    assert (
        set(truth["window"])
        == set(gamma["window"])
        == {
            "focal",
            "local_100kb",
        }
    )


@pytest.mark.parametrize("source", ("truth", "gamma"))
def test_score_reducers_require_all_100_within_diploid_pairs(source) -> None:
    frame = _truth_profile() if source == "truth" else _gamma_profile()
    frame.loc[0, "n_pairs"] = 99
    reducer = study._truth_scores if source == "truth" else study._gamma_scores
    with pytest.raises(ValueError, match="exactly 100 pairs"):
        reducer(frame)


@pytest.mark.parametrize("source", ("truth", "gamma"))
def test_score_reducers_reject_inconsistent_one_based_positions(source) -> None:
    frame = _truth_profile() if source == "truth" else _gamma_profile()
    frame.loc[0, "position_1based"] += 1
    reducer = study._truth_scores if source == "truth" else study._gamma_scores
    with pytest.raises(ValueError, match="one-based positions"):
        reducer(frame)


@pytest.mark.parametrize("kind", ("status", "unit_id", "extra"))
@pytest.mark.parametrize("phase", ("simulation", "decode"))
def test_completion_validator_rejects_top_level_tampering(
    tmp_path, phase, kind
) -> None:
    expected_contract = {"unit": {"unit_id": "unit_1"}}
    if phase == "simulation":
        destination = tmp_path
        filename = "simulation_complete.json"
        completion = {
            "schema": study.SIMULATION_SCHEMA,
            "status": "complete",
            "unit_id": "unit_1",
            "contract": expected_contract,
            "contract_sha256": study._canonical_sha256(expected_contract),
            "metadata": {},
            "elapsed_seconds": 1.0,
            "outputs": {},
        }
        validator = study._validate_simulation_completion
    else:
        destination = tmp_path / "decoded"
        destination.mkdir()
        filename = "decode_complete.json"
        completion = {
            "schema": study.DECODE_SCHEMA,
            "status": "complete",
            "unit_id": "unit_1",
            "contract": expected_contract,
            "contract_sha256": study._canonical_sha256(expected_contract),
            "elapsed_seconds": 1.0,
            "outputs": {},
        }
        validator = study._validate_decode_completion
    if kind == "status":
        completion["status"] = "failed"
        message = "status is not complete"
    elif kind == "unit_id":
        completion["unit_id"] = "wrong"
        message = "unit_id changed"
    else:
        completion["unexpected"] = True
        message = "top-level fields changed"
    (destination / filename).write_text(json.dumps(completion), encoding="utf-8")
    with pytest.raises(ValueError, match=message):
        validator(tmp_path, expected_contract)


def test_output_records_reject_extra_fields(tmp_path) -> None:
    path = tmp_path / "artifact.tsv"
    path.write_text("x\n1\n", encoding="utf-8")
    record = study._output_record(path, rows=1)
    record["unexpected"] = True
    with pytest.raises(ValueError, match="record fields changed"):
        study._validate_output_record(record, path, expected_rows=1)


def test_binary_contract_identity_must_be_repo_relative(tmp_path) -> None:
    inside = tmp_path / "bin/tool"
    inside.parent.mkdir()
    inside.write_bytes(b"binary")
    assert study._portable_binary_identity(inside, tmp_path) == "bin/tool"
    assert study._resolve_binary_identity("bin/tool", tmp_path, "test") == inside
    with pytest.raises(ValueError, match="inside repo_root"):
        study._portable_binary_identity(Path("/").resolve(), tmp_path)
    with pytest.raises(ValueError, match="repo relative"):
        study._resolve_binary_identity("../tool", tmp_path, "test")
