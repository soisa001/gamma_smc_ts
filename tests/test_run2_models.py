"""Configuration and model-construction tests for the run2 study.

None of these need the SLiM binary: stdpopsim can generate the script it would
have run, which is enough to assert every scheduling invariant.
"""

from __future__ import annotations

import json

import msprime
import numpy as np
import pytest

stdpopsim = pytest.importorskip("stdpopsim")

from gamma_smc_aou import run2_config as config  # noqa: E402
from gamma_smc_aou import run2_models as models  # noqa: E402
from gamma_smc_aou import run2_workflow as workflow  # noqa: E402

pytestmark = pytest.mark.skipif(
    stdpopsim.__version__ != config.REQUIRED_STDPOPSIM_VERSION,
    reason=f"run2 requires stdpopsim {config.REQUIRED_STDPOPSIM_VERSION}",
)


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------


def test_seeds_are_unique_across_arms_and_modes():
    seeds = [
        arm.seed(mode, index)
        for arm in config.ARMS
        for mode in config.MODES
        for index in range(config.N_REPLICATES)
    ]
    assert len(set(seeds)) == len(seeds) == 400


def test_seed_is_a_pure_function_of_arm_mode_index():
    assert config.CHB_ARM.seed("selected", 7) == config.CHB_ARM.seed("selected", 7)
    assert config.CHB_ARM.seed("selected", 7) != config.CHB_ARM.seed("neutral", 7)
    assert config.CHB_ARM.seed("selected", 7) != config.EAS_ARM.seed("selected", 7)


def test_seed_rejects_unknown_mode():
    with pytest.raises(ValueError):
        config.CHB_ARM.seed("conditioned", 0)


def test_arm_config_records_declare_no_conditioning():
    for arm in config.ARMS:
        record = config.arm_config_record(arm)
        shared = record["shared"]
        assert shared["present_day_af_conditioning"] is None
        assert shared["allele_survival_conditioning"] is None
        assert shared["null_conditioning"] is None
        assert shared["standing_frequency"] == 0.0296
        json.dumps(record)


def test_config_files_are_written(tmp_path):
    written = config.write_arm_configs(tmp_path)
    assert set(written) == {arm.arm_id for arm in config.ARMS}
    for path in written.values():
        assert json.loads(path.read_text(encoding="utf-8"))["schema"]


# ---------------------------------------------------------------------------
# Tick arithmetic
# ---------------------------------------------------------------------------


def test_snap_generations_matches_stdpopsim_rounding():
    assert models.snap_generations(2272.0) == 2270.0
    assert models.snap_generations(2016.0) == 2015.0
    assert models.snap_generations(stdpopsim.GenerationAfter(2272.0)) == 2265.0


def test_standing_variation_lands_one_tick_after_the_pulse():
    for arm in config.ARMS:
        schedule = models.tick_schedule(arm)
        offsets = schedule["tick_offsets_before_final_tick"]
        assert offsets["standing_variation"] == offsets["pulse"] - 1
        assert schedule["realized_generations_ago"]["pulse"] == 2270.0
        assert schedule["realized_generations_ago"]["standing_variation"] == 2265.0


def test_han_split_tick_is_shared_by_both_fitness_callbacks():
    schedule = models.tick_schedule(config.CHB_ARM)
    intervals = schedule["selection_intervals_generations_ago"]
    assert intervals[0]["population"] == "Loschbour"
    assert intervals[1]["population"] == "Han"
    assert intervals[0]["end"] == intervals[1]["start"] == 2015.0
    assert intervals[1]["end"] == 0.0


def test_selection_covers_the_whole_post_introgression_interval():
    for arm in config.ARMS:
        intervals = models.tick_schedule(arm)["selection_intervals_generations_ago"]
        assert intervals[0]["start"] == 2265.0
        assert intervals[-1]["end"] == 0.0
        for earlier, later in zip(intervals, intervals[1:]):
            assert earlier["end"] == later["start"], "selection must not have gaps"


# ---------------------------------------------------------------------------
# Events
# ---------------------------------------------------------------------------


def test_neutral_mode_has_no_extended_events():
    for arm in config.ARMS:
        assert models.build_extended_events(arm, "neutral") == ()


def test_selected_mode_has_no_conditioning_events():
    for arm in config.ARMS:
        events = models.build_extended_events(arm, "selected")
        assert not any(
            isinstance(event, stdpopsim.ConditionOnAlleleFrequency) for event in events
        )
        draws = [e for e in events if isinstance(e, stdpopsim.DrawMutation)]
        fitness = [e for e in events if isinstance(e, stdpopsim.ChangeMutationFitness)]
        assert len(draws) == 1
        assert len(fitness) == (2 if arm.introgression else 1)
        assert all(e.selection_coeff == 0.01 for e in fitness)
        assert all(e.dominance_coeff == 0.5 for e in fitness)


def test_focal_draw_targets_the_pulse_recipient():
    events = models.build_extended_events(config.CHB_ARM, "selected")
    draw = next(e for e in events if isinstance(e, stdpopsim.DrawMutation))
    assert draw.population == "Loschbour"
    assert float(draw.time) == 2272.0
    assert isinstance(draw.time, stdpopsim.GenerationAfter)


def test_build_extended_events_rejects_unknown_mode():
    with pytest.raises(ValueError):
        models.build_extended_events(config.CHB_ARM, "conditioned")


# ---------------------------------------------------------------------------
# Focal-base guard
# ---------------------------------------------------------------------------


def test_focal_mask_removes_exactly_one_base_of_mutational_mass():
    rate = 1.25e-8
    uniform = msprime.RateMap(
        position=np.array([0.0, 1e7]), rate=np.array([rate])
    )
    masked, changed = models.mask_focal_rate_map(uniform, 5_000_000)
    assert changed
    total = masked.get_cumulative_mass(1e7)
    assert total == pytest.approx(rate * (1e7 - 1), rel=0, abs=1e-18)
    focal = masked.get_cumulative_mass(5_000_001.0) - masked.get_cumulative_mass(
        5_000_000.0
    )
    assert focal == 0.0


def test_focal_mask_is_idempotent_and_skips_zero_rate_maps():
    uniform = msprime.RateMap(position=np.array([0.0, 1e7]), rate=np.array([1e-8]))
    masked, _ = models.mask_focal_rate_map(uniform, 5_000_000)
    _, changed = models.mask_focal_rate_map(masked, 5_000_000)
    assert not changed
    zeroed = msprime.RateMap(position=np.array([0.0, 1e7]), rate=np.array([0.0]))
    _, changed = models.mask_focal_rate_map(zeroed, 5_000_000)
    assert not changed


def test_focal_overlay_patch_restores_msprime_on_exit():
    original = msprime.sim_mutations
    with models.scoped_focal_overlay_patch() as receipt:
        assert msprime.sim_mutations is not original
        assert receipt["intercepted_calls"] == 0
    assert msprime.sim_mutations is original
    assert receipt["restored"] is True


def test_focal_overlay_patch_restores_on_exception():
    original = msprime.sim_mutations
    with pytest.raises(RuntimeError):
        with models.scoped_focal_overlay_patch():
            raise RuntimeError("boom")
    assert msprime.sim_mutations is original


# ---------------------------------------------------------------------------
# SLiM patch scoping
# ---------------------------------------------------------------------------


def test_slim_patch_restores_stdpopsim_helpers():
    from stdpopsim import slim_engine

    original = slim_engine._slim_functions
    with models.scoped_slim_patch(config.CHB_ARM, "selected", 5) as record:
        assert slim_engine._slim_functions is not original
        assert record["add_mut_replaced"] and record["end_replaced"]
    assert slim_engine._slim_functions is original


def test_slim_patch_is_a_no_op_for_the_neutral_mode():
    from stdpopsim import slim_engine

    original = slim_engine._slim_functions
    with models.scoped_slim_patch(config.CHB_ARM, "neutral", 5) as record:
        assert slim_engine._slim_functions is original
        assert not record["add_mut_replaced"]
        assert not record["end_replaced"]


# ---------------------------------------------------------------------------
# Generated script
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def chb_selected_script(repo_root):
    return models.generate_slim_script(config.CHB_ARM, "selected", repo_root)


def test_generated_script_never_conditions(repo_root):
    for arm in config.ARMS:
        for mode in config.MODES:
            script = models.generate_slim_script(arm, mode, repo_root)
            checks = workflow.check_slim_script(arm, mode, script)
            assert checks["rejection_sampling"] is False
            assert checks["passed"]


def test_introgression_arm_selects_carriers_by_migrant_status(chb_selected_script):
    assert "carriers = inds[inds.migrant];" in chb_selected_script
    assert "run2_standing_frequency" in chb_selected_script


def test_non_introgression_arm_places_the_standing_frequency_directly(repo_root):
    script = models.generate_slim_script(config.EAS_ARM, "selected", repo_root)
    assert "0.0296 * size(all_genomes)" in script
    assert "inds.migrant" not in script


def test_census_end_handles_loss_and_fixation(chb_selected_script):
    assert '"lost"' in chb_selected_script
    assert '"fixed_as_substitution"' in chb_selected_script
    assert "run2_final_census_af" in chb_selected_script


def test_model_record_is_serializable(repo_root):
    for arm in config.ARMS:
        for mode in config.MODES:
            record = models.model_record(arm, mode, repo_root)
            assert record["conditioning_events"] == 0
            json.dumps(record)
