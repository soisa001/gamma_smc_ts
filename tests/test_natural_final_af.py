from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import gamma_smc_aou.natural_final_af as natural_af
from gamma_smc_aou.natural_final_af import (
    ANALYSIS_OUTPUT_NAMES,
    ATTEMPT_COLUMNS,
    DOMINANCE_COEFFICIENT,
    EAS_ORIGIN_GENERATIONS,
    HAN_LOSCHBOUR_N,
    HAN_POST_SPLIT_GENERATIONS,
    HAN_PRE_SPLIT_GENERATIONS,
    HAN_PRESENT_N,
    PANEL_DIPLOIDS,
    TrajectorySchedule,
    build_eas_schedule,
    build_han_schedule,
    selected_gamete_frequency,
    simulate_trajectory_batch,
    summarize_attempts,
)


REPO_ROOT = Path(__file__).resolve().parents[1]


def _cell(coefficient: float) -> dict:
    return {
        "cell_id": f"test_s{coefficient}",
        "demography_id": "test",
        "population": "EAS",
        "simulation_class": "neutral" if coefficient == 0 else "selected",
        "selection_coefficient": coefficient,
        "dominance_coefficient": DOMINANCE_COEFFICIENT,
        "attempts": 8,
        "batch_size": 8,
        "n_batches": 1,
    }


def _short_schedule(initial_mode: str = "de_novo_one_copy") -> TrajectorySchedule:
    return TrajectorySchedule(
        demography_id="test",
        ages=np.asarray([3, 2, 1, 0], dtype=np.int64),
        diploid_sizes=np.asarray([20, 20, 30, 30], dtype=np.int64),
        populations=("EAS", "EAS", "EAS", "EAS"),
        initial_mode=initial_mode,
        initial_frequency=0.0296 if initial_mode != "de_novo_one_copy" else None,
        contract={"integer_schedule_sha256": "test"},
    )


def test_additive_selection_transition_is_valid_and_increases_intermediate_p():
    p = np.asarray([0.0, 0.01, 0.25, 0.5, 1.0])
    neutral = selected_gamete_frequency(p, 0.0)
    selected = selected_gamete_frequency(p, 0.01, 0.5)
    np.testing.assert_allclose(neutral, p)
    assert np.all(selected[1:-1] > p[1:-1])
    assert selected[0] == 0
    assert selected[-1] == 1


def test_han_schedule_is_exact_q1_recipient_lineage():
    schedule = build_han_schedule(validate_catalog=False)
    assert len(schedule.ages) == 2_273
    assert len(schedule.diploid_sizes) == 2_273
    assert np.all(
        schedule.diploid_sizes[: HAN_PRE_SPLIT_GENERATIONS + 1] == HAN_LOSCHBOUR_N
    )
    assert np.all(
        schedule.diploid_sizes[HAN_PRE_SPLIT_GENERATIONS + 1 :] == HAN_PRESENT_N
    )
    assert schedule.ages[0] == HAN_PRE_SPLIT_GENERATIONS + HAN_POST_SPLIT_GENERATIONS
    assert schedule.ages[-1] == 0
    assert schedule.contract["initial_count_distribution"] == "Binomial(2*2340, 0.0296)"
    assert schedule.contract["pulse_state_convention"] == (
        "immediately post-pulse at generation 2272"
    )
    assert schedule.contract["stdpopsim_version"] == "0.3.0"


def test_eas_schedule_matches_pointwise_median_step_contract():
    schedule = build_eas_schedule(REPO_ROOT)
    assert len(schedule.ages) == EAS_ORIGIN_GENERATIONS + 1
    assert schedule.diploid_sizes[0] == 3_632
    assert schedule.diploid_sizes[-1] == 37_637
    assert schedule.contract["integer_schedule_rows"] == 2_001
    assert len(schedule.contract["integer_schedule_sha256"]) == 64


def test_batch_is_deterministic_retains_losses_and_draws_one_panel():
    schedule = _short_schedule()
    first = simulate_trajectory_batch(
        _cell(0.0),
        schedule,
        attempt_start=0,
        attempt_stop=8,
        batch_index=0,
        batch_seed=12345,
    )
    second = simulate_trajectory_batch(
        _cell(0.0),
        schedule,
        attempt_start=0,
        attempt_stop=8,
        batch_index=0,
        batch_seed=12345,
    )
    assert first.equals(second)
    assert len(first) == 8
    assert np.array_equal(first["population_survived"], first["final_alt_count"] > 0)
    assert np.array_equal(first["sample_detected"], first["sample_alt_count"] > 0)
    assert set(first["panel_diploids"]) == {PANEL_DIPLOIDS}
    assert first["sample_alt_count"].between(0, 2 * PANEL_DIPLOIDS).all()


def test_han_initial_counts_are_drawn_in_loschbour_not_han():
    schedule = TrajectorySchedule(
        demography_id="test",
        ages=np.asarray([1, 0], dtype=np.int64),
        diploid_sizes=np.asarray([HAN_LOSCHBOUR_N, HAN_PRESENT_N], dtype=np.int64),
        populations=("Loschbour", "Han"),
        initial_mode="donor_fixed_introgression_pulse",
        initial_frequency=0.0296,
        contract={"integer_schedule_sha256": "test"},
    )
    result = simulate_trajectory_batch(
        {
            **_cell(0.0),
            "demography_id": "ancient_eurasia_han_introgression",
            "population": "Han",
        },
        schedule,
        attempt_start=0,
        attempt_stop=2_000,
        batch_index=0,
        batch_seed=54321,
    )
    assert set(result["initial_diploid_n"]) == {HAN_LOSCHBOUR_N}
    assert result["initial_alt_count"].mean() == pytest.approx(
        2 * HAN_LOSCHBOUR_N * 0.0296, rel=0.03
    )


def test_summary_filters_only_on_population_survival():
    neutral = simulate_trajectory_batch(
        _cell(0.0),
        _short_schedule(),
        attempt_start=0,
        attempt_stop=8,
        batch_index=0,
        batch_seed=101,
    )
    selected = simulate_trajectory_batch(
        _cell(0.01),
        _short_schedule(),
        attempt_start=0,
        attempt_stop=8,
        batch_index=0,
        batch_seed=202,
    )
    attempts = pd.concat([neutral, selected], ignore_index=True)
    summary, comparisons = summarize_attempts(attempts)
    assert set(summary["n_attempted"]) == {8}
    expected = attempts.groupby("simulation_class")["population_survived"].sum()
    observed = summary.set_index("simulation_class")["n_population_survived"]
    assert observed.to_dict() == expected.astype(int).to_dict()
    assert len(comparisons) == 1
    assert np.array_equal(
        summary["n_lost"] + summary["n_segregating"] + summary["n_fixed"],
        summary["n_attempted"],
    )
    assert summary["population_survival_wilson_ci_low"].between(0, 1).all()
    assert summary["population_survival_wilson_ci_high"].between(0, 1).all()
    by_class = summary.set_index("simulation_class")
    assert np.isfinite(by_class.loc["neutral", "neutral_martingale_z"])
    assert np.isnan(by_class.loc["selected", "neutral_martingale_z"])


def _write_completed_batch(
    tmp_path: Path,
    frame: pd.DataFrame,
    contract: dict,
) -> tuple[Path, Path]:
    table = (
        tmp_path / "batches" / str(contract["cell"]["cell_id"]) / "batch_00000.tsv.gz"
    )
    completion = table.with_name("batch_00000.completion.json")
    natural_af._atomic_frame(table, frame)
    payload = {
        "schema": natural_af.BATCH_SCHEMA,
        "status": "complete",
        "contract": contract,
        "contract_sha256": natural_af._canonical_sha256(contract),
        "output": {
            "path": contract["output_relative_path"],
            "sha256": natural_af.sha256_file(table),
            "size_bytes": table.stat().st_size,
            "rows": len(frame),
            "columns": list(frame.columns),
        },
    }
    natural_af._atomic_json(completion, payload)
    return table, completion


def test_deterministic_gzip_has_no_temporary_filename_header(tmp_path):
    frame = pd.DataFrame({"x": [1, 2], "value": [0.1, 0.2]})
    first = tmp_path / "first.tsv.gz"
    second = tmp_path / "second.tsv.gz"
    natural_af._atomic_frame(first, frame)
    natural_af._atomic_frame(second, frame)
    assert first.read_bytes() == second.read_bytes()
    assert first.read_bytes()[3] & 0x08 == 0
    assert natural_af._tsv_payload_sha256(first) == natural_af._frame_tsv_sha256(frame)


def test_canonical_source_hash_is_crlf_portable_and_rejects_bad_text(tmp_path):
    lf = tmp_path / "lf.py"
    crlf = tmp_path / "crlf.py"
    lf.write_bytes(b"x = 1\ny = 2\n")
    crlf.write_bytes(b"x = 1\r\ny = 2\r\n")
    assert natural_af.canonical_lf_sha256(lf) == natural_af.canonical_lf_sha256(crlf)
    crlf.write_bytes(b"x = 1\ry = 2\n")
    with pytest.raises(ValueError, match="bare carriage return"):
        natural_af.canonical_lf_sha256(crlf)
    crlf.write_bytes(b"\xef\xbb\xbfx = 1\n")
    with pytest.raises(ValueError, match="UTF-8 BOM"):
        natural_af.canonical_lf_sha256(crlf)


def test_batch_validation_rejects_semantic_tampering_with_fresh_checksum(tmp_path):
    cell = _cell(0.0)
    schedule = _short_schedule()
    plan_payload = {"schema": natural_af.SCHEMA_VERSION, "rng": {"base_seed": 7}}
    contract = natural_af._batch_contract(
        plan_payload,
        cell,
        schedule,
        batch_index=0,
        attempt_start=0,
        attempt_stop=8,
        batch_seed=12345,
    )
    frame = simulate_trajectory_batch(
        cell,
        schedule,
        attempt_start=0,
        attempt_stop=8,
        batch_index=0,
        batch_seed=12345,
    )
    table, completion = _write_completed_batch(tmp_path, frame, contract)
    natural_af._validate_batch(table, completion, contract)
    tampered = frame.copy()
    tampered.loc[0, "final_population_af"] = 0.987654321
    table, completion = _write_completed_batch(tmp_path, tampered, contract)
    with pytest.raises(ValueError, match="terminal count/AF"):
        natural_af._validate_batch(table, completion, contract)


def test_table_only_batch_is_restartably_regenerated(tmp_path, monkeypatch):
    cell = _cell(0.0)
    schedule = _short_schedule()
    payload = {"schema": natural_af.SCHEMA_VERSION, "rng": {"base_seed": 13}}
    work = tmp_path / "work"
    table, completion = natural_af._batch_paths(work, str(cell["cell_id"]), 0)
    natural_af._atomic_frame(table, pd.DataFrame(columns=ATTEMPT_COLUMNS))
    monkeypatch.setattr(natural_af, "_schedule_for_cell", lambda *_: schedule)
    result = natural_af._simulate_batch_task(
        {
            "repo_root": str(tmp_path),
            "work_dir": str(work),
            "plan_payload": payload,
            "cell": cell,
            "batch_index": 0,
        }
    )
    assert result["status"] == "complete"
    assert completion.is_file()
    assert len(pd.read_csv(table, sep="\t")) == 8


def test_unconditional_fate_figure_is_letter_png_and_pdf(tmp_path):
    summary = pd.DataFrame(
        [
            {
                "demography_id": demography,
                "simulation_class": simulation_class,
                "n_attempted": 10,
                "lost_fraction": 0.6,
                "segregating_fraction": 0.3,
                "fixed_fraction": 0.1,
                "population_survival_fraction": 0.4,
                "fixation_fraction_among_survivors": 0.25,
            }
            for demography in (
                "eas_phlash_median",
                "ancient_eurasia_han_introgression",
            )
            for simulation_class in ("neutral", "selected")
        ]
    )
    paths = natural_af._plot_unconditional_fates(summary, tmp_path)
    assert {path.name for path in paths} == {
        ANALYSIS_OUTPUT_NAMES["fate_png"],
        ANALYSIS_OUTPUT_NAMES["fate_pdf"],
    }
    assert natural_af._figure_geometry(
        tmp_path / ANALYSIS_OUTPUT_NAMES["fate_png"]
    ) == {
        "width_px": 3_300,
        "height_px": 2_550,
    }
    assert natural_af._figure_geometry(
        tmp_path / ANALYSIS_OUTPUT_NAMES["fate_pdf"]
    ) == {
        "width_points": 792.0,
        "height_points": 612.0,
        "page_count": 1,
    }
