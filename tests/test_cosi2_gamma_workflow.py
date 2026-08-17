from pathlib import Path

import pandas as pd
import pytest
from gamma_smc_aou import cosi2_gamma_analysis as analysis
from gamma_smc_aou import cosi2_gamma_workflow as workflow


def _unit_inventory() -> pd.DataFrame:
    rows = []
    for demography, settings in workflow.DEMOGRAPHY_SETTINGS.items():
        seed_start = (
            workflow.neutral.EAS_SEED_START
            if demography == "EAS"
            else workflow.neutral.HAN_SEED_START
        )
        for index in range(workflow.EXPECTED_NEUTRAL_PER_DEMOGRAPHY):
            count = index % 21
            rows.append(
                {
                    "unit_id": f"{settings['neutral_model_id']}_r{index:03d}",
                    "demography": demography,
                    "demography_id": settings["demography_id"],
                    "model_id": settings["neutral_model_id"],
                    "simulation_class": "neutral",
                    "selection_coefficient": 0.0,
                    "seed": seed_start + index,
                    "generation_time_years": settings["generation_time_years"],
                    "present_ne": settings["present_ne"],
                    "final_population_af": 0.001 + index / 1000,
                    "sample_alt_count": count,
                    "sample_af": count / analysis.PANEL_HAPLOTYPES,
                    "ms_path": f"/read-only/{demography}/neutral_{index}.ms",
                    "ms_sha256": f"{index + 1000 * (demography == 'Han'):064x}"[-64:],
                    "source_completion_path": (
                        f"/read-only/{demography}/neutral_{index}.json"
                    ),
                    "source_completion_sha256": f"{index + 100:064x}"[-64:],
                }
            )
        selected_count = 40 if demography == "EAS" else 32
        rows.append(
            {
                "unit_id": settings["selected_unit_id"],
                "demography": demography,
                "demography_id": settings["demography_id"],
                "model_id": settings["selected_cell_id"],
                "simulation_class": "selected",
                "selection_coefficient": workflow.SELECTION_COEFFICIENT,
                "seed": settings["selected_seed"],
                "generation_time_years": settings["generation_time_years"],
                "present_ne": settings["present_ne"],
                "final_population_af": 0.2,
                "sample_alt_count": selected_count,
                "sample_af": selected_count / analysis.PANEL_HAPLOTYPES,
                "ms_path": f"/read-only/{demography}/selected.ms",
                "ms_sha256": "a" * 64,
                "source_completion_path": f"/read-only/{demography}/selected.json",
                "source_completion_sha256": "b" * 64,
            }
        )
    return workflow.validate_unit_inventory(pd.DataFrame(rows))


def _scores(units: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for unit in units.to_dict(orient="records"):
        for statistic_index, statistic in enumerate(analysis.SUPPORTED_STATISTICS):
            for window_index, window in enumerate(analysis.SUPPORTED_WINDOWS):
                for threshold_index, threshold in enumerate(
                    analysis.TMRCA_THRESHOLDS_YEARS
                ):
                    rows.append(
                        {
                            **{
                                column: unit[column]
                                for column in analysis.UNIT_METADATA_COLUMNS
                            },
                            "statistic": statistic,
                            "window": window,
                            "threshold_years": threshold,
                            "score": min(
                                0.99,
                                0.03
                                + threshold_index * 0.08
                                + statistic_index * 0.01
                                + window_index * 0.005
                                + 0.12 * (unit["simulation_class"] == "selected"),
                            ),
                        }
                    )
    return analysis.validate_scores(pd.DataFrame(rows))


def _decode_records(units: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "unit_id": units["unit_id"],
            "status": "complete",
            "decode_completion_sha256": [
                f"{index + 10_000:064x}" for index in range(len(units))
            ],
            "decode_full_sample_alt_count": units["sample_alt_count"].astype(int),
            "decode_pair_panel_alt_count": (
                units["sample_alt_count"].astype(int) // 2
            ).clip(upper=100),
            "decode_focal_count_status": "unique_exact_focal_site",
        }
    )


def _fake_plot(scores, output_dir, *, stem):
    assert len(scores) == workflow.EXPECTED_TOTAL_UNITS * 28
    destination = Path(output_dir)
    png = destination / f"{stem}.png"
    pdf = destination / f"{stem}.pdf"
    png.write_bytes(b"\x89PNG\r\n\x1a\nsynthetic")
    pdf.write_bytes(b"%PDF-1.4\nsynthetic\n%%EOF\n")
    return {"png": png, "pdf": pdf}


def test_exact_unit_inventory_and_sample_count_convention():
    units = _unit_inventory()
    counts = units.groupby(["demography", "simulation_class"]).size().to_dict()
    assert counts == {
        ("EAS", "neutral"): 100,
        ("EAS", "selected"): 1,
        ("Han", "neutral"): 100,
        ("Han", "selected"): 1,
    }
    assert (
        workflow._neutral_sample_alt_count(
            {
                "present_population_frequency": 0.01,
                "reported_focal_position_count": 0,
                "focal_alt_count_if_unambiguous": None,
            }
        )
        == 0
    )
    assert (
        workflow._neutral_sample_alt_count(
            {
                "present_population_frequency": 1.0,
                "reported_focal_position_count": 0,
                "focal_alt_count_if_unambiguous": None,
            }
        )
        == 200
    )
    with pytest.raises(ValueError, match="exactly 202"):
        workflow.validate_unit_inventory(units.iloc[:-1])


def test_atomic_result_bundle_cache_and_tamper_detection(monkeypatch, tmp_path):
    units = _unit_inventory()
    scores = _scores(units)
    decodes = _decode_records(units)
    cosi2 = tmp_path / "coalescent"
    decoder = tmp_path / "gamma_smc"
    cosi2.write_bytes(b"synthetic cosi2")
    decoder.write_bytes(b"synthetic gamma")
    monkeypatch.setattr(analysis, "plot_score_profiles", _fake_plot)
    destination = tmp_path / "results"

    completion = workflow.write_results_bundle(
        tmp_path,
        results_dir=destination,
        units=units,
        scores=scores,
        decode_records=decodes,
        decoder_bin=decoder,
        cosi2_binary=cosi2,
    )
    assert completion["status"] == "complete"
    assert completion["cache_hit"] is False
    assert {path.name for path in destination.iterdir()} == {
        workflow.RESULTS_COMPLETION_FILENAME,
        *workflow.RESULT_OUTPUTS.values(),
    }
    assert len(list(destination.glob("*.png"))) == 1
    assert len(list(destination.glob("*.pdf"))) == 1
    persisted_scores = pd.read_csv(
        destination / workflow.RESULT_OUTPUTS["scores"], sep="\t"
    )
    assert len(persisted_scores) == workflow.EXPECTED_TOTAL_UNITS * 28
    pointwise = pd.read_csv(
        destination / workflow.RESULT_OUTPUTS["pointwise_pvalues"], sep="\t"
    )
    assert set(pointwise["n_neutral"]) == {100}
    assert pointwise["minimum_attainable_p"].unique() == pytest.approx([1 / 101])
    unit_metadata = pd.read_csv(
        destination / workflow.RESULT_OUTPUTS["unit_metadata"], sep="\t"
    )
    selected = unit_metadata[unit_metadata["simulation_class"] == "selected"]
    assert set(selected["decode_full_sample_alt_count"]) == {32, 40}
    assert set(selected["decode_pair_panel_alt_count"]) == {16, 20}

    cached = workflow.write_results_bundle(
        tmp_path,
        results_dir=destination,
        units=units,
        scores=scores,
        decode_records=decodes,
        decoder_bin=decoder,
        cosi2_binary=cosi2,
    )
    assert cached["cache_hit"] is True

    target = destination / workflow.RESULT_OUTPUTS["scores"]
    target.write_bytes(target.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="checksum failed"):
        workflow.verify_results_bundle(tmp_path, results_dir=destination)


def test_result_bundle_rejects_nonempty_unowned_destination(monkeypatch, tmp_path):
    units = _unit_inventory()
    scores = _scores(units)
    cosi2 = tmp_path / "coalescent"
    decoder = tmp_path / "gamma_smc"
    cosi2.write_bytes(b"cosi")
    decoder.write_bytes(b"gamma")
    destination = tmp_path / "results"
    destination.mkdir()
    (destination / "orphan.txt").write_text("partial\n", encoding="utf-8")
    monkeypatch.setattr(analysis, "plot_score_profiles", _fake_plot)
    with pytest.raises((ValueError, TypeError)):
        workflow.write_results_bundle(
            tmp_path,
            results_dir=destination,
            units=units,
            scores=scores,
            decode_records=_decode_records(units),
            decoder_bin=decoder,
            cosi2_binary=cosi2,
        )


def test_cli_exposes_all_restartable_phases():
    action = next(
        action for action in workflow.build_parser()._actions if action.dest == "phase"
    )
    assert set(action.choices) == {
        "plan",
        "simulate",
        "status",
        "decode",
        "analyze",
        "verify",
        "all",
    }
