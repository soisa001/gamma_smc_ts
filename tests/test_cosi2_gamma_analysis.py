from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from gamma_smc_aou import cosi2_gamma_analysis as analysis
from scipy.stats import beta


def _synthetic_scores(*, n_neutral: int = 4) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for demography_index, demography in enumerate(("EAS", "HAN")):
        units = [
            (
                f"{demography}_neutral_{neutral_index:03d}",
                "neutral",
                0.0,
                10_000 * (demography_index + 1) + neutral_index,
                0.01 + 0.001 * neutral_index,
                2,
            )
            for neutral_index in range(n_neutral)
        ]
        units.append(
            (
                f"{demography}_selected_000",
                "selected",
                0.01,
                10_000 * (demography_index + 1) + n_neutral,
                0.20,
                40 if demography == "EAS" else 32,
            )
        )
        for (
            unit_id,
            simulation_class,
            coefficient,
            seed,
            population_af,
            sample_alt_count,
        ) in units:
            neutral_index = int(unit_id.rsplit("_", 1)[-1])
            for statistic_index, statistic in enumerate(analysis.SUPPORTED_STATISTICS):
                for window_index, window in enumerate(analysis.SUPPORTED_WINDOWS):
                    offset = 0.01 * statistic_index + 0.02 * window_index
                    for threshold_index, threshold in enumerate(
                        analysis.TMRCA_THRESHOLDS_YEARS
                    ):
                        if simulation_class == "neutral":
                            score = (
                                0.10
                                + 0.10 * threshold_index
                                + 0.01 * neutral_index
                                + offset
                            )
                        else:
                            score = 0.135 + 0.10 * threshold_index + offset
                        rows.append(
                            {
                                "unit_id": unit_id,
                                "demography": demography,
                                "simulation_class": simulation_class,
                                "selection_coefficient": coefficient,
                                "seed": seed,
                                "generation_time_years": (
                                    25.0 if demography == "EAS" else 29.0
                                ),
                                "statistic": statistic,
                                "window": window,
                                "threshold_years": threshold,
                                "score": score,
                                "final_population_af": population_af,
                                "sample_alt_count": sample_alt_count,
                                "sample_af": sample_alt_count
                                / analysis.PANEL_HAPLOTYPES,
                            }
                        )
    return pd.DataFrame(rows)


def test_validate_scores_accepts_only_complete_monotone_unit_curves() -> None:
    source = _synthetic_scores()
    validated = analysis.validate_scores(source.sample(frac=1.0, random_state=17))
    assert tuple(validated.columns) == analysis.INPUT_COLUMNS
    assert len(validated) == len(source)
    assert set(validated["statistic"]) == set(analysis.SUPPORTED_STATISTICS)
    assert set(validated["window"]) == set(analysis.SUPPORTED_WINDOWS)
    assert set(validated["threshold_years"].astype(int)) == set(
        analysis.TMRCA_THRESHOLDS_YEARS
    )
    rows_per_unit = (
        len(analysis.SUPPORTED_STATISTICS)
        * len(analysis.SUPPORTED_WINDOWS)
        * len(analysis.TMRCA_THRESHOLDS_YEARS)
    )
    assert set(validated.groupby("unit_id").size()) == {rows_per_unit}


@pytest.mark.parametrize(
    ("mutation", "message"),
    (
        ("missing_column", "missing columns"),
        ("unknown_statistic", "two supported score definitions"),
        ("missing_curve_row", "exact complete score grid"),
        ("nonmonotone", "nonmonotone curve"),
        ("changing_metadata", "unit metadata changes"),
        ("sample_af_mismatch", "sample_af disagrees"),
        ("lost_allele", "must be positive"),
        ("demography_generation_time_change", "changes within demography"),
    ),
)
def test_validate_scores_rejects_malformed_inputs(mutation: str, message: str) -> None:
    frame = _synthetic_scores()
    if mutation == "missing_column":
        frame = frame.drop(columns="seed")
    elif mutation == "unknown_statistic":
        frame.loc[
            frame["statistic"] == analysis.SUPPORTED_STATISTICS[0], "statistic"
        ] = "unsupported"
    elif mutation == "missing_curve_row":
        frame = frame.drop(index=frame.index[0])
    elif mutation == "nonmonotone":
        unit_id = frame.iloc[0]["unit_id"]
        mask = (
            (frame["unit_id"] == unit_id)
            & (frame["statistic"] == analysis.SUPPORTED_STATISTICS[0])
            & (frame["window"] == analysis.SUPPORTED_WINDOWS[0])
            & (frame["threshold_years"] == 4_500)
        )
        frame.loc[mask, "score"] = 0.0
    elif mutation == "changing_metadata":
        frame.loc[frame.index[0], "final_population_af"] = 0.8
    elif mutation == "sample_af_mismatch":
        frame.loc[frame.index[0], "sample_af"] = 0.5
    elif mutation == "lost_allele":
        frame.loc[frame["unit_id"] == "EAS_neutral_000", "final_population_af"] = 0.0
    elif mutation == "demography_generation_time_change":
        frame.loc[frame["unit_id"] == "HAN_neutral_000", "generation_time_years"] = 28.0
    else:  # pragma: no cover - protects the test fixture itself
        raise AssertionError(mutation)
    with pytest.raises(ValueError, match=message):
        analysis.validate_scores(frame)


def test_clopper_pearson_interval_has_exact_boundaries() -> None:
    low, high = analysis.clopper_pearson_interval(0, 4)
    assert low == 0.0
    assert high == pytest.approx(beta.ppf(0.975, 1, 4))
    low, high = analysis.clopper_pearson_interval(4, 4)
    assert low == pytest.approx(beta.ppf(0.025, 4, 1))
    assert high == 1.0
    with pytest.raises(ValueError, match="0 <= successes"):
        analysis.clopper_pearson_interval(5, 4)


def test_pointwise_pvalues_use_independent_units_and_conservative_ties() -> None:
    frame = _synthetic_scores(n_neutral=4)
    result = analysis.pointwise_pvalues(frame)
    expected_rows = (
        2
        * len(analysis.SUPPORTED_STATISTICS)
        * len(analysis.SUPPORTED_WINDOWS)
        * len(analysis.TMRCA_THRESHOLDS_YEARS)
    )
    assert len(result) == expected_rows
    first = result[
        (result["demography"] == "EAS")
        & (result["statistic"] == analysis.SUPPORTED_STATISTICS[0])
        & (result["window"] == analysis.PRIMARY_WINDOW)
        & (result["threshold_years"] == 1_000)
    ].iloc[0]
    assert first["n_neutral"] == 4
    assert first["neutral_upper_exceedance_count"] == 0
    assert first["neutral_lower_exceedance_count"] == 4
    assert first["mc_p_upper"] == pytest.approx(1 / 5)
    assert first["mc_p_lower"] == 1.0
    assert first["mc_p_two_sided"] == pytest.approx(2 / 5)
    assert first["minimum_attainable_p"] == pytest.approx(1 / 5)
    assert first["generation_time_years"] == 25.0
    assert first["threshold_generations"] == 40.0
    assert first["neutral_midrank_descriptive_auc"] == 1.0
    assert first["upper_tail_cp95_low"] == 0.0
    assert first["upper_tail_cp95_high"] == pytest.approx(beta.ppf(0.975, 1, 4))
    assert bool(first["primary_window"])

    han = result[
        (result["demography"] == "HAN")
        & (result["statistic"] == analysis.SUPPORTED_STATISTICS[0])
        & (result["window"] == analysis.PRIMARY_WINDOW)
        & (result["threshold_years"] == 1_000)
    ].iloc[0]
    assert han["generation_time_years"] == 29.0
    assert han["threshold_generations"] == pytest.approx(1_000 / 29)

    selected_id = "EAS_selected_000"
    tie_mask = (
        (frame["unit_id"] == selected_id)
        & (frame["statistic"] == analysis.SUPPORTED_STATISTICS[0])
        & (frame["window"] == analysis.PRIMARY_WINDOW)
        & (frame["threshold_years"] == 1_000)
    )
    frame.loc[tie_mask, "score"] = 0.12
    tied = analysis.pointwise_pvalues(frame)
    tied = tied[
        (tied["unit_id"] == selected_id)
        & (tied["statistic"] == analysis.SUPPORTED_STATISTICS[0])
        & (tied["window"] == analysis.PRIMARY_WINDOW)
        & (tied["threshold_years"] == 1_000)
    ].iloc[0]
    assert tied["neutral_upper_exceedance_count"] == 2
    assert tied["neutral_lower_exceedance_count"] == 3
    assert tied["neutral_tie_count"] == 1
    assert tied["mc_p_upper"] == pytest.approx(3 / 5)
    assert tied["mc_p_lower"] == pytest.approx(4 / 5)
    assert tied["mc_p_two_sided"] == 1.0
    assert tied["neutral_midrank_descriptive_auc"] == pytest.approx(0.625)


def test_minp_omnibus_retains_complete_seven_threshold_unit_vectors() -> None:
    result = analysis.minp_omnibus(_synthetic_scores(n_neutral=4))
    assert len(result) == 2 * len(analysis.SUPPORTED_STATISTICS) * len(
        analysis.SUPPORTED_WINDOWS
    )
    assert set(result["n_neutral"]) == {4}
    assert set(result["n_thresholds"]) == {7}
    assert set(result["thresholds_years"]) == {
        "1000,4500,10000,20000,30000,40000,50000"
    }
    assert np.allclose(result["omnibus_mc_p_upper"], 1 / 5)
    assert np.allclose(result["omnibus_mc_p_lower"], 1.0)
    assert np.allclose(result["omnibus_mc_p_two_sided"], 2 / 5)
    assert set(result["threshold_years_at_min_pointwise_p_upper"]) == {1_000}
    assert set(result["neutral_minp_as_or_more_extreme_count_upper"]) == {0}
    assert set(result["omnibus_method"]) == {
        "exchangeable_whole_unit_minP_over_seven_thresholds"
    }


def test_analyze_scores_returns_only_unit_level_inference_tables() -> None:
    scores = _synthetic_scores()
    result = analysis.analyze_scores(scores)
    assert set(result) == {"scores", "pointwise_pvalues", "minp_omnibus"}
    assert result["scores"]["unit_id"].nunique() == 10
    assert set(result["pointwise_pvalues"]["n_neutral"]) == {4}
    assert set(result["minp_omnibus"]["n_neutral"]) == {4}


def test_plot_score_profiles_builds_letter_landscape_without_rendering(
    tmp_path, monkeypatch
) -> None:
    frame = _synthetic_scores()
    captured: dict[str, object] = {}

    def capture(figure, output_dir, stem):
        captured["size"] = tuple(figure.get_size_inches())
        captured["axis_count"] = len(figure.axes)
        captured["legend_labels"] = tuple(
            text.get_text() for text in figure.legends[0].get_texts()
        )
        return {
            "png": output_dir / f"{stem}.png",
            "pdf": output_dir / f"{stem}.pdf",
        }

    monkeypatch.setattr(analysis, "_atomic_save_figure_pair", capture)
    outputs = analysis.plot_score_profiles(frame, tmp_path)
    assert set(outputs) == {"png", "pdf"}
    assert captured["size"] == pytest.approx((11.0, 8.5))
    assert captured["axis_count"] == 4
    assert captured["legend_labels"] == (
        "+/-100 kb (primary)",
        "Focal",
        "Selected median",
        "Neutral median",
        "Neutral 95% range",
    )
    assert not outputs["png"].exists()
    assert not outputs["pdf"].exists()
