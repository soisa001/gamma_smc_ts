from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from gamma_smc_aou import han_s001_tmrca_analysis as analysis


def _synthetic_scores() -> pd.DataFrame:
    rows = []
    specifications = (
        ("neutral", 0.0, (0, 10, 20, 30)),
        ("selected", 0.001, (10, 20, 40)),
    )
    for simulation_class, coefficient, sample_counts in specifications:
        for replicate_index, sample_count in enumerate(sample_counts, start=1):
            unit_id = f"{simulation_class}__rep{replicate_index:03d}"
            final_af = sample_count / 200 + 0.01
            for source_index, source in enumerate(analysis.ALLOWED_SOURCES):
                for window_index, window in enumerate(analysis.ALLOWED_WINDOWS):
                    panel_offset = 0.002 * source_index + 0.001 * window_index
                    for threshold_index, threshold in enumerate(
                        analysis.THRESHOLDS_YEARS
                    ):
                        if simulation_class == "neutral":
                            score = (
                                0.10 * threshold_index
                                + 0.01 * replicate_index
                                + panel_offset
                            )
                        else:
                            score = (
                                0.11 * threshold_index
                                + 0.002 * replicate_index
                                + panel_offset
                            )
                        rows.append(
                            {
                                "unit_id": unit_id,
                                "simulation_class": simulation_class,
                                "selection_coefficient": coefficient,
                                "seed": (
                                    replicate_index
                                    + (100 if simulation_class == "selected" else 0)
                                ),
                                "source": source,
                                "window": window,
                                "threshold_years": threshold,
                                "score": score,
                                "final_population_af": final_af,
                                "sample_alt_count": sample_count,
                                "sample_af": sample_count / 200,
                                "sample_detected": sample_count > 0,
                            }
                        )
    return pd.DataFrame(rows, columns=analysis.SCORE_COLUMNS)


def _complete_empirical() -> pd.DataFrame:
    result = analysis.empirical_input_template()
    result["empirical_id"] = "han_locus_1"
    result["variant_id"] = "rs-test"
    result["score"] = [
        0.05 + 0.10 * analysis.THRESHOLDS_YEARS.index(int(threshold))
        for threshold in result["threshold_years"]
    ]
    result["sample_alt_count"] = 20
    result["sample_af"] = 0.1
    result["sample_detected"] = True
    for index, column in enumerate(
        (
            "source_data_sha256",
            "sample_manifest_sha256",
            "decoder_binary_sha256",
            "decoder_parameters_sha256",
            "ne_history_sha256",
        ),
        start=1,
    ):
        result[column] = str(index) * 64
    result["qc_status"] = "pass"
    result["data_status"] = "empirical_observed"
    return result


def test_primary_auc_is_unflipped_shared_and_deterministic():
    scores = analysis.validate_replicate_scores(_synthetic_scores())
    assert len(scores) == 7 * 4 * 7

    auc, omnibus = analysis.auc_by_threshold(
        scores, bootstrap_draws=39, permutation_draws=79, base_seed=91
    )
    cell = auc[
        (auc["source"] == "gamma_smc") & (auc["window"] == "local_100kb")
    ].sort_values("threshold_years")
    assert tuple(cell["threshold_years"].astype(int)) == analysis.THRESHOLDS_YEARS
    assert cell.iloc[0]["roc_auc"] == pytest.approx(0.0)
    assert cell.iloc[-1]["roc_auc"] == pytest.approx(1.0)
    assert cell["permutation_seed"].nunique() == 1
    assert set(cell["permutation_family"]) == {
        "shared_whole_unit_labels_across_seven_thresholds"
    }
    assert set(cell["score_direction"]) == {"larger_is_more_selected_no_flipping"}
    assert len(omnibus) == 4
    assert omnibus["max_auc_permutation_p_upper"].between(0, 1).all()

    repeated_auc, repeated_omnibus = analysis.auc_by_threshold(
        scores, bootstrap_draws=39, permutation_draws=79, base_seed=91
    )
    pd.testing.assert_frame_equal(auc, repeated_auc)
    pd.testing.assert_frame_equal(omnibus, repeated_omnibus)


def test_selected_tail_af_benchmarks_and_secondary_stratification():
    scores = _synthetic_scores()
    pvalues = analysis.selected_unit_neutral_tail_pvalues(scores)
    first = pvalues[
        (pvalues["source"] == "gamma_smc")
        & (pvalues["window"] == "local_100kb")
        & (pvalues["threshold_years"] == 1_000)
    ]
    last = pvalues[
        (pvalues["source"] == "gamma_smc")
        & (pvalues["window"] == "local_100kb")
        & (pvalues["threshold_years"] == 50_000)
    ]
    assert set(first["mc_p_upper"]) == {1.0}
    assert set(last["mc_p_upper"]) == {0.2}
    assert set(last["minimum_attainable_p"]) == {0.2}

    benchmarks = analysis.af_only_benchmarks(
        scores, bootstrap_draws=29, permutation_draws=59, base_seed=12
    )
    assert set(benchmarks["metric"]) == {"final_population_af", "sample_af"}
    assert benchmarks["roc_auc"].between(0, 1).all()
    assert benchmarks["interpretation"].str.contains("benchmark_only").all()

    stratified, omnibus = analysis.sample_af_stratified_auc(
        scores,
        af_stratum_width=0.05,
        bootstrap_draws=29,
        permutation_draws=59,
        base_seed=12,
    )
    assert len(stratified) == 4 * 7
    assert len(omnibus) == 4
    assert set(stratified["estimand"]) == {
        "secondary_sample_af_stratified_common_support_pair_weighted"
    }
    assert stratified["interpretation"].str.startswith("secondary").all()
    assert stratified["roc_auc"].between(0, 1).all()
    assert (
        stratified.groupby(["source", "window"])["permutation_seed"].nunique() == 1
    ).all()


@pytest.mark.parametrize(
    ("mutation", "message"),
    (
        ("duplicate", "duplicate key"),
        ("missing_threshold", "exact threshold curve"),
        ("wrong_neutral_s", "neutral rows"),
        ("wrong_selected_s", "selected rows"),
        ("bad_sample_af", "sample_af disagrees"),
        ("bad_detection", "sample_detected disagrees"),
        ("nonmonotone", "nonmonotone"),
        ("duplicate_seed", "unique seeds"),
        ("missing_unit_id", "unit_id contains a missing value"),
        ("unexpected_column", "unexpected"),
    ),
)
def test_validation_fails_closed(mutation: str, message: str):
    frame = _synthetic_scores()
    if mutation == "duplicate":
        frame = pd.concat([frame, frame.iloc[[0]]], ignore_index=True)
    elif mutation == "missing_threshold":
        row = frame.iloc[0]
        frame = frame[
            ~(
                (frame["unit_id"] == row["unit_id"])
                & (frame["source"] == row["source"])
                & (frame["window"] == row["window"])
                & (frame["threshold_years"] == row["threshold_years"])
            )
        ]
    elif mutation == "wrong_neutral_s":
        frame.loc[frame["simulation_class"] == "neutral", "selection_coefficient"] = (
            0.001
        )
    elif mutation == "wrong_selected_s":
        frame.loc[frame["simulation_class"] == "selected", "selection_coefficient"] = (
            0.005
        )
    elif mutation == "bad_sample_af":
        frame.loc[frame.index[0], "sample_af"] = 0.9
    elif mutation == "bad_detection":
        frame.loc[frame["sample_alt_count"] == 0, "sample_detected"] = True
    elif mutation == "nonmonotone":
        mask = (
            (frame["unit_id"] == "neutral__rep001")
            & (frame["source"] == "gamma_smc")
            & (frame["window"] == "local_100kb")
            & (frame["threshold_years"] == 4_500)
        )
        frame.loc[mask, "score"] = 0.0
    elif mutation == "duplicate_seed":
        frame.loc[frame["unit_id"] == "selected__rep001", "seed"] = 1
    elif mutation == "missing_unit_id":
        frame.loc[frame.index[0], "unit_id"] = pd.NA
    elif mutation == "unexpected_column":
        frame["not_in_contract"] = 1
    with pytest.raises(ValueError, match=message):
        analysis.validate_replicate_scores(frame)


def test_empirical_input_is_unavailable_until_explicit_and_validated():
    scores = _synthetic_scores()
    unavailable = analysis.empirical_neutral_pvalues(None, scores)
    assert set(unavailable["status"]) == {"unavailable_no_explicit_empirical_table"}
    assert set(unavailable["estimand"]) == {
        analysis.PRIMARY_EMPIRICAL_ESTIMAND,
        analysis.SECONDARY_EMPIRICAL_ESTIMAND,
    }
    assert unavailable["mc_p_upper"].isna().all()

    unavailable_omnibus = analysis.empirical_neutral_omnibus_pvalues(None, scores)
    assert set(unavailable_omnibus["status"]) == {
        "unavailable_no_explicit_empirical_table"
    }
    assert unavailable_omnibus["omnibus_mc_p_upper"].isna().all()

    template = analysis.empirical_input_template()
    assert len(template) == 2 * 7
    with pytest.raises(ValueError, match="unset|finite"):
        analysis.validate_empirical_input(template)

    empirical = _complete_empirical()
    validated = analysis.validate_empirical_input(empirical)
    assert len(validated) == 14
    results = analysis.empirical_neutral_pvalues(
        validated, scores, af_match_tolerance=0.025
    )
    assert set(results["status"]) == {"available_explicit_empirical_input"}
    assert len(results) == 2 * len(validated)
    assert not results.duplicated(
        [
            "empirical_id",
            "source",
            "window",
            "threshold_years",
            "estimand",
            "null_model",
        ]
    ).any()
    primary = results[results["estimand"] == analysis.PRIMARY_EMPIRICAL_ESTIMAND]
    secondary = results[results["estimand"] == analysis.SECONDARY_EMPIRICAL_ESTIMAND]
    assert set(primary["null_model"]) == {analysis.PRIMARY_EMPIRICAL_NULL_MODEL}
    assert set(secondary["null_model"]) == {analysis.SECONDARY_EMPIRICAL_NULL_MODEL}
    assert set(primary["n_matched_neutral"]) == {4}
    assert set(secondary["n_matched_neutral"]) == {1}
    assert primary["af_match_tolerance"].isna().all()
    assert not primary["sample_detection_matched"].any()
    assert secondary["sample_detection_matched"].all()
    assert (results["n_matched_neutral"] > 0).all()
    assert results["mc_p_upper"].between(0, 1).all()
    assert set(results["source"]) == {"gamma_smc"}

    omnibus = analysis.empirical_neutral_omnibus_pvalues(
        validated, scores, af_match_tolerance=0.025
    )
    assert len(omnibus) == 2 * len(analysis.ALLOWED_WINDOWS)
    assert set(omnibus["estimand"]) == {
        analysis.PRIMARY_EMPIRICAL_ESTIMAND,
        analysis.SECONDARY_EMPIRICAL_ESTIMAND,
    }
    assert set(
        omnibus.loc[
            omnibus["estimand"] == analysis.PRIMARY_EMPIRICAL_ESTIMAND,
            "n_shared_neutral_profiles",
        ]
    ) == {4}
    assert set(
        omnibus.loc[
            omnibus["estimand"] == analysis.SECONDARY_EMPIRICAL_ESTIMAND,
            "n_shared_neutral_profiles",
        ]
    ) == {1}
    assert omnibus["omnibus_mc_p_upper"].between(0, 1).all()

    invalid = empirical.copy()
    invalid["source"] = "tree_truth"
    with pytest.raises(ValueError, match="source must be gamma_smc"):
        analysis.validate_empirical_input(invalid)


def test_empirical_validation_requires_both_windows_for_every_id():
    first = _complete_empirical()
    second = _complete_empirical()
    second["empirical_id"] = "han_locus_2"
    second["variant_id"] = "rs-test-2"
    incomplete = pd.concat(
        [first[first["window"] == "local_100kb"], second], ignore_index=True
    )
    with pytest.raises(ValueError, match="empirical_id han_locus_1 lacks both windows"):
        analysis.validate_empirical_input(incomplete)


@pytest.mark.parametrize("column", ("empirical_id", "variant_id"))
def test_empirical_validation_rejects_missing_identifiers(column: str):
    empirical = _complete_empirical()
    empirical.loc[empirical.index[0], column] = pd.NA
    with pytest.raises(
        ValueError, match=f"empirical provenance field is missing: {column}"
    ):
        analysis.validate_empirical_input(empirical)


def test_af_strata_separate_undetected_zero_from_detected_low_af():
    strata = analysis._af_strata(
        np.asarray([0.0, 0.005]),
        0.05,
        np.asarray([False, True]),
    )
    assert tuple(strata) == (-1, 0)
    assert analysis._format_af_strata(strata) == (
        "sample_undetected_af_zero,sample_detected_af_bin_0"
    )


@pytest.mark.parametrize(
    ("argument", "value"),
    (
        ("bootstrap_draws", True),
        ("bootstrap_draws", 0),
        ("bootstrap_draws", 1.5),
        ("permutation_draws", False),
        ("permutation_draws", 2.5),
    ),
)
def test_draw_counts_are_strict_positive_integers(argument: str, value):
    arguments = {"bootstrap_draws": 3, "permutation_draws": 3}
    arguments[argument] = value
    with pytest.raises(ValueError, match=f"{argument} must be a positive integer"):
        analysis.auc_by_threshold(_synthetic_scores(), **arguments)


@pytest.mark.parametrize("value", (True, 1.5, -1))
def test_public_base_seed_is_a_strict_nonnegative_integer(value, tmp_path):
    scores = _synthetic_scores()
    calls = (
        lambda: analysis.auc_by_threshold(
            scores, bootstrap_draws=3, permutation_draws=3, base_seed=value
        ),
        lambda: analysis.af_only_benchmarks(
            scores, bootstrap_draws=3, permutation_draws=3, base_seed=value
        ),
        lambda: analysis.sample_af_stratified_auc(
            scores, bootstrap_draws=3, permutation_draws=3, base_seed=value
        ),
        lambda: analysis.analyze_replicate_scores(
            scores, bootstrap_draws=3, permutation_draws=3, base_seed=value
        ),
        lambda: analysis.plot_selected_neutral_curves(
            scores, tmp_path, bootstrap_draws=3, base_seed=value
        ),
    )
    for call in calls:
        with pytest.raises(ValueError, match="base_seed must be a nonnegative integer"):
            call()


def test_analysis_and_atomic_table_interface_without_rendering(tmp_path):
    scores = _synthetic_scores()
    results = analysis.analyze_replicate_scores(
        scores,
        bootstrap_draws=19,
        permutation_draws=29,
        base_seed=7,
        include_af_stratified=True,
    )
    assert set(results) == set(analysis.TABLE_FILENAMES)
    assert results["empirical_pvalues"].iloc[0]["status"].startswith("unavailable")
    assert results["empirical_omnibus"].iloc[0]["status"].startswith("unavailable")
    assert set(results["empirical_input"]["data_status"]) == {"template_unfilled"}
    paths = analysis.write_analysis_tables(results, tmp_path)
    assert set(paths) == set(analysis.TABLE_FILENAMES)
    assert all(path.is_file() and path.stat().st_size > 0 for path in paths.values())
    restored = pd.read_csv(paths["selected_unit_pvalues"], sep="\t")
    assert len(restored) == len(results["selected_unit_pvalues"])
    assert not list(tmp_path.glob("*.tmp.*"))

    input_path = tmp_path / "replicate_scores.tsv"
    scores.to_csv(input_path, sep="\t", index=False)
    loaded = analysis.load_replicate_scores(input_path)
    assert len(loaded) == len(scores)

    second_gzip = tmp_path / "same_content_different_name.tsv.gz"
    analysis._atomic_tsv(results["selected_unit_pvalues"], second_gzip)
    assert paths["selected_unit_pvalues"].read_bytes() == second_gzip.read_bytes()


def test_analysis_preserves_validated_explicit_empirical_provenance():
    empirical = _complete_empirical()
    expected = analysis.validate_empirical_input(empirical)
    results = analysis.analyze_replicate_scores(
        _synthetic_scores(),
        empirical=empirical,
        bootstrap_draws=3,
        permutation_draws=3,
    )
    pd.testing.assert_frame_equal(results["empirical_input"], expected)
    assert set(results["empirical_input"]["data_status"]) == {"empirical_observed"}
    assert set(results["empirical_pvalues"]["status"]) == {
        "available_explicit_empirical_input"
    }


def test_plot_functions_build_letter_layout_without_rendering(monkeypatch, tmp_path):
    captured = []

    def capture(figure, output_dir, stem):
        captured.append((figure, output_dir, stem))
        return tmp_path / f"{stem}.png", tmp_path / f"{stem}.pdf"

    monkeypatch.setattr(analysis, "_save_figure_pair", capture)
    scores = _synthetic_scores()
    analysis.plot_selected_neutral_curves(
        scores, tmp_path, bootstrap_draws=5, base_seed=3
    )
    results = analysis.analyze_replicate_scores(
        scores,
        bootstrap_draws=5,
        permutation_draws=9,
        base_seed=3,
    )
    analysis.plot_auc_profiles(
        results["auc_by_threshold"],
        results["af_only_benchmarks"],
        tmp_path,
        stratified_auc=results["sample_af_stratified_auc"],
    )
    assert len(captured) == 2
    for figure, output_dir, stem in captured:
        assert tuple(figure.get_size_inches()) == pytest.approx((11, 8.5))
        assert output_dir == tmp_path
        assert stem.startswith("han_s001_")
        assert figure._suptitle.get_fontsize() >= 20
        assert max(axis.get_position().y1 for axis in figure.axes) <= 0.83
        analysis.plt.close(figure)


def test_atomic_figure_disables_tight_cropping_without_rendering(tmp_path):
    calls = []

    class FakeFigure:
        def savefig(self, path, **kwargs):
            calls.append(kwargs)
            path.write_bytes(b"synthetic image bytes")

    output = tmp_path / "letter.png"
    analysis._atomic_figure(FakeFigure(), output)
    assert output.read_bytes() == b"synthetic image bytes"
    pdf_output = tmp_path / "letter.pdf"
    analysis._atomic_figure(FakeFigure(), pdf_output)
    assert pdf_output.read_bytes() == b"synthetic image bytes"
    assert calls == [
        {"format": "png", "dpi": 300, "bbox_inches": None},
        {
            "format": "pdf",
            "dpi": None,
            "bbox_inches": None,
            "metadata": {"CreationDate": None, "ModDate": None},
        },
    ]
