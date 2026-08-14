from __future__ import annotations

import hashlib
import json
import struct

import numpy as np
import pandas as pd
import pytest
from gamma_smc_aou import focused_selection_analysis as analysis

THRESHOLDS = (1_000, 4_500, 10_000, 20_000, 30_000, 40_000, 50_000)
CLASSES = ("overall", "hom_ref", "heterozygous", "hom_alt")


def _synthetic_classes() -> pd.DataFrame:
    rows = []
    neutral_contrasts = (-0.01, 0.0, 0.005, 0.01)
    for simulation_class, contrasts in (
        ("selected", (0.20, 0.25)),
        ("neutral", neutral_contrasts),
    ):
        for replicate_index, contrast in enumerate(contrasts, start=1):
            unit_id = f"eas__af10__{simulation_class}__rep{replicate_index:03d}"
            for threshold_index, threshold in enumerate(THRESHOLDS):
                reference = 0.08 + 0.10 * threshold_index + 0.001 * replicate_index
                values = {
                    "hom_ref": reference,
                    "hom_alt": min(reference + contrast, 0.99),
                }
                values["heterozygous"] = 0.5 * (values["hom_ref"] + values["hom_alt"])
                values["overall"] = values["heterozygous"]
                for genotype_class in CLASSES:
                    selected = simulation_class == "selected"
                    tmrca_ref = 1_000.0
                    tmrca_alt = 500.0 if selected else 1_000.0 - 100 * contrast
                    tmrca = {
                        "hom_ref": tmrca_ref,
                        "hom_alt": tmrca_alt,
                        "heterozygous": 0.5 * (tmrca_ref + tmrca_alt),
                        "overall": 0.5 * (tmrca_ref + tmrca_alt),
                    }[genotype_class]
                    rows.append(
                        {
                            "unit_id": unit_id,
                            "demography_id": "eas_phlash_median",
                            "simulation_class": simulation_class,
                            "selection_coefficient": (
                                0.01 if simulation_class == "selected" else 0.0
                            ),
                            "target_allele_frequency": 0.1,
                            "replicate_index": replicate_index,
                            "source": "gamma_smc",
                            "genotype_class": genotype_class,
                            "n_pairs": {
                                "overall": 100,
                                "hom_ref": 82,
                                "heterozygous": 16,
                                "hom_alt": 2,
                            }[genotype_class],
                            "threshold_years": threshold,
                            "region_mean_p_tmrca_lt_threshold": values[genotype_class],
                            "focal_mean_p_tmrca_lt_threshold": values[genotype_class],
                            "region_mean_tmrca_generations": tmrca,
                            "focal_mean_tmrca_generations": tmrca,
                        }
                    )
    return pd.DataFrame(rows)


def _complete_empirical_rows(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    hashes = {
        "source_data_sha256": "1" * 64,
        "sample_manifest_sha256": "2" * 64,
        "pair_table_sha256": "3" * 64,
        "reference_decoder_binary_sha256": "4" * 64,
        "empirical_decoder_binary_sha256": "4" * 64,
        "reference_decoder_parameters_sha256": "5" * 64,
        "empirical_decoder_parameters_sha256": "5" * 64,
        "reference_ne_history_sha256": "6" * 64,
        "empirical_ne_history_sha256": "6" * 64,
        "mask_manifest_sha256": "7" * 64,
    }
    for column, value in hashes.items():
        result[column] = value
    result["empirical_id"] = "empirical_locus_1"
    result["variant_id"] = "rs-test-1"
    result["source_label"] = "empirical_test_dataset"
    result["demography_assignment"] = "EAS ancestry, PHLASH-median null stratum"
    result["population"] = "EAS"
    result["observed_allele_frequency"] = 0.1
    result["n_diploid_samples"] = 100
    result["n_overall_pairs"] = 100
    result["n_hom_ref_pairs"] = 82
    result["n_heterozygous_pairs"] = 16
    result["n_hom_alt_pairs"] = 2
    result["alt_allele_count"] = 20
    result["decoder_version"] = "test-version"
    result["genome_build"] = "GRCh38"
    result["contig"] = "chr1"
    result["region_start_0based"] = 10_000_000
    result["region_end_0based"] = 20_000_000
    result["focal_position_0based"] = 15_000_000
    result["position_1based"] = 15_000_001
    result["ref_allele"] = "A"
    result["alt_allele"] = "G"
    result["observed_p_hom_ref"] = 0.10
    result["observed_p_heterozygous"] = 0.25
    result["observed_p_hom_alt"] = 0.40
    result["observed_p_overall"] = 0.20
    observed = {
        analysis.PRIMARY_STATISTIC: 0.30,
        "p_hom_alt": 0.40,
        "p_hom_ref": 0.10,
        "p_heterozygous": 0.25,
        "p_overall": 0.20,
    }
    result["observed_value"] = result["statistic"].map(observed)
    result["qc_status"] = "pass"
    return result


def _synthetic_pairs() -> pd.DataFrame:
    classes = _synthetic_classes()
    rows = []
    selected = classes[classes["simulation_class"] == "selected"]
    for unit_id, unit in selected.groupby("unit_id", sort=True):
        metadata = unit.iloc[0][
            [
                "unit_id",
                "demography_id",
                "simulation_class",
                "selection_coefficient",
                "target_allele_frequency",
                "replicate_index",
            ]
        ].to_dict()
        for threshold, threshold_group in unit.groupby("threshold_years"):
            indexed = threshold_group.set_index("genotype_class")
            for pair_index, genotype_class, offset in (
                (0, "hom_ref", -0.001),
                (1, "hom_ref", 0.001),
                (2, "hom_alt", -0.001),
                (3, "hom_alt", 0.001),
                (4, "heterozygous", 0.0),
            ):
                record = indexed.loc[genotype_class]
                rows.append(
                    {
                        **metadata,
                        "pair_index": pair_index,
                        "genotype_class": genotype_class,
                        "threshold_years": threshold,
                        "region_mean_p_tmrca_lt_threshold": record[
                            "region_mean_p_tmrca_lt_threshold"
                        ]
                        + offset,
                        "focal_p_tmrca_lt_threshold": record[
                            "focal_mean_p_tmrca_lt_threshold"
                        ]
                        + offset,
                        "region_mean_tmrca_generations": record[
                            "region_mean_tmrca_generations"
                        ],
                        "focal_tmrca_generations": record[
                            "focal_mean_tmrca_generations"
                        ],
                    }
                )
    return pd.DataFrame(rows)


def _synthetic_spatial(source: str = "gamma_smc") -> pd.DataFrame:
    classes = _synthetic_classes().copy()
    classes["source"] = source
    unit_table = classes[
        [
            "unit_id",
            "demography_id",
            "simulation_class",
            "selection_coefficient",
            "target_allele_frequency",
        ]
    ].drop_duplicates()
    counts = (
        unit_table.groupby(
            [
                "demography_id",
                "simulation_class",
                "selection_coefficient",
                "target_allele_frequency",
            ]
        )
        .size()
        .to_dict()
    )
    rows = []
    for cell, n_units in counts.items():
        demography, simulation_class, coefficient, af = cell
        for genotype_index, genotype_class in enumerate(CLASSES):
            for threshold_index, threshold in enumerate(THRESHOLDS):
                for position in (0, 5_000_000, 9_990_000):
                    probability = min(
                        0.99,
                        0.03
                        + 0.08 * threshold_index
                        + 0.01 * genotype_index
                        + (0.05 if simulation_class == "selected" else 0),
                    )
                    rows.append(
                        {
                            "source": source,
                            "demography_id": demography,
                            "simulation_class": simulation_class,
                            "selection_coefficient": coefficient,
                            "target_allele_frequency": af,
                            "genotype_class": genotype_class,
                            "threshold_years": threshold,
                            "position_0based": position,
                            "n_units": n_units,
                            "mean_p_tmrca_lt_threshold": probability,
                            "sd_p_tmrca_lt_threshold": 0.01,
                            "sem_p_tmrca_lt_threshold": 0.01 / np.sqrt(n_units),
                            "mean_tmrca_generations": 900.0 - 10 * genotype_index,
                            "sd_tmrca_generations": 20.0,
                            "sem_tmrca_generations": 20.0 / np.sqrt(n_units),
                        }
                    )
    return pd.DataFrame(rows)


def test_matched_null_pvalues_loo_auc_and_within_selected_effects():
    classes = analysis.validate_class_summaries(
        _synthetic_classes(),
        minimum_selected_per_cell=2,
        minimum_neutral_per_cell=4,
    )
    scores = analysis.build_replicate_statistics(classes)
    primary = scores[
        (scores["statistic_scope"] == "focal_nearest")
        & (scores["metric"] == "p_tmrca_lt_threshold")
        & (scores["statistic"] == analysis.PRIMARY_STATISTIC)
        & (scores["threshold_years"] == 1_000)
    ]
    np.testing.assert_allclose(
        primary.loc[primary["simulation_class"] == "selected", "selection_score"],
        [0.20, 0.25],
    )

    pvalues = analysis.selected_vs_neutral_pvalues(scores)
    selected = pvalues[
        (pvalues["statistic_scope"] == "focal_nearest")
        & (pvalues["statistic"] == analysis.PRIMARY_STATISTIC)
        & (pvalues["threshold_years"] == 1_000)
    ]
    assert set(selected["mc_p_upper"]) == {0.2}
    assert set(selected["n_null"]) == {4}
    assert (selected["null_percentile_midrank"] == 1).all()
    power = analysis.conditional_power_summary(pvalues, alpha=0.2)
    primary_power = power[
        (power["statistic_scope"] == "focal_nearest")
        & (power["statistic"] == analysis.PRIMARY_STATISTIC)
        & (power["threshold_years"] == 1_000)
    ].iloc[0]
    assert primary_power["n_selected"] == 2
    assert primary_power["n_detected"] == 2
    assert primary_power["conditional_power"] == 1
    assert primary_power["conditional_power_raw"] == 1
    assert primary_power["conditional_power_bh_within_7"] == 1
    assert primary_power["bh_family_size"] == len(THRESHOLDS)
    assert 0 < primary_power["wilson_ci95_low"] < 1
    adjusted_fixture = pvalues.copy()
    adjusted_mask = (
        (adjusted_fixture["statistic_scope"] == "focal_nearest")
        & (adjusted_fixture["statistic"] == analysis.PRIMARY_STATISTIC)
        & (adjusted_fixture["threshold_years"] == 1_000)
    )
    adjusted_fixture.loc[adjusted_mask, "mc_p_upper"] = 0.01
    adjusted_fixture.loc[adjusted_mask, "bh_q_upper"] = 0.07
    adjusted_power = analysis.conditional_power_summary(adjusted_fixture, alpha=0.05)
    adjusted_primary = adjusted_power[
        (adjusted_power["statistic_scope"] == "focal_nearest")
        & (adjusted_power["statistic"] == analysis.PRIMARY_STATISTIC)
        & (adjusted_power["threshold_years"] == 1_000)
    ].iloc[0]
    assert adjusted_primary["conditional_power_raw"] == 1
    assert adjusted_primary["conditional_power_bh_within_7"] == 0

    loo, calibration = analysis.neutral_leave_one_out(scores)
    loo_primary = loo[
        (loo["statistic_scope"] == "focal_nearest")
        & (loo["statistic"] == analysis.PRIMARY_STATISTIC)
        & (loo["threshold_years"] == 1_000)
    ]
    assert len(loo_primary) == 4
    assert set(loo_primary["n_loo_null"]) == {3}
    assert set(loo_primary["diagnostic_role"]) == {"descriptive_only_not_out_of_sample"}
    assert not calibration.empty

    crossfit, crossfit_calibration = analysis.neutral_crossfit_calibration(
        scores, n_folds=2, base_seed=17
    )
    crossfit_again, _ = analysis.neutral_crossfit_calibration(
        scores, n_folds=2, base_seed=17
    )
    pd.testing.assert_frame_equal(crossfit, crossfit_again)
    assert set(crossfit["n_training_null"]) == {2}
    assert (
        crossfit.groupby("unit_id")["crossfit_fold"].nunique().to_numpy() == 1
    ).all()
    assert not crossfit_calibration.empty

    omnibus = analysis.primary_threshold_omnibus(scores)
    assert len(omnibus) == 2 * 2
    assert set(omnibus["n_thresholds"]) == {len(THRESHOLDS)}
    assert set(omnibus["omnibus_p_upper"]) == {0.2}
    assert set(omnibus["omnibus_method"]) == {
        "pooled_leave_one_out_minP_exact_unit_exchangeability"
    }

    comparisons = analysis.cell_roc_comparisons(
        scores, bootstrap_draws=50, permutation_draws=29, base_seed=7
    )
    primary_auc = comparisons[
        (comparisons["statistic_scope"] == "focal_nearest")
        & (comparisons["statistic"] == analysis.PRIMARY_STATISTIC)
        & (comparisons["threshold_years"] == 1_000)
    ].iloc[0]
    assert primary_auc["roc_auc"] == 1
    assert primary_auc["roc_auc_ci95_low"] == 1
    assert primary_auc["mean_selected_minus_neutral"] > 0
    assert 0 < primary_auc["mean_unit_label_permutation_p_upper"] <= 1
    assert 0 < primary_auc["auc_unit_label_permutation_p_upper"] <= 1
    assert primary_auc["unit_label_permutation_draws"] == 29
    repeated_comparisons = analysis.cell_roc_comparisons(
        scores, bootstrap_draws=50, permutation_draws=29, base_seed=7
    )
    pd.testing.assert_series_equal(
        comparisons["mean_unit_label_permutation_p_upper"],
        repeated_comparisons["mean_unit_label_permutation_p_upper"],
    )
    pd.testing.assert_series_equal(
        comparisons["auc_unit_label_permutation_p_upper"],
        repeated_comparisons["auc_unit_label_permutation_p_upper"],
    )

    within = analysis.within_selected_summary(scores, bootstrap_draws=50, base_seed=7)
    within_primary = within[
        (within["statistic_scope"] == "focal_nearest")
        & (within["statistic"] == analysis.PRIMARY_STATISTIC)
        & (within["threshold_years"] == 1_000)
    ].iloc[0]
    assert within_primary["n_positive"] == 2
    assert within_primary["sign_test_p_upper"] == 0.25


def test_two_selection_coefficients_share_null_without_pooling_selected(tmp_path):
    frame = _synthetic_classes()
    second = frame[frame["simulation_class"] == "selected"].copy()
    second["selection_coefficient"] = 0.005
    second["unit_id"] = second["unit_id"].str.replace(
        "__selected__", "__s005__selected__", regex=False
    )
    combined = analysis.validate_class_summaries(
        pd.concat([frame, second], ignore_index=True),
        minimum_selected_per_cell=2,
        minimum_neutral_per_cell=4,
    )
    scores = analysis.build_replicate_statistics(combined)
    comparisons = analysis.cell_roc_comparisons(
        scores, bootstrap_draws=20, permutation_draws=19, base_seed=11
    )
    primary = comparisons[
        (comparisons["statistic_scope"] == "focal_nearest")
        & (comparisons["statistic"] == analysis.PRIMARY_STATISTIC)
        & (comparisons["threshold_years"] == 1_000)
    ]
    assert set(primary["selection_coefficient"]) == {0.005, 0.01}
    assert set(primary["n_selected"]) == {2}
    assert set(primary["n_neutral"]) == {4}
    assert set(primary["neutral_bank_scope"]) == {
        "shared_by_demography_and_af_across_selection_coefficients"
    }
    figures = analysis.plot_primary_auc(comparisons, tmp_path)
    figures += analysis.plot_primary_pvalues(comparisons, tmp_path)
    figures += analysis.plot_probability_curves(combined, tmp_path)
    figures += analysis.plot_primary_contrast_curves(scores, tmp_path)
    assert figures
    assert all(path.is_file() and path.stat().st_size > 0 for path in figures)


def test_spatial_profiles_validate_counts_and_write_png_pdf_with_source_label(tmp_path):
    classes = analysis.validate_class_summaries(
        _synthetic_classes(),
        minimum_selected_per_cell=2,
        minimum_neutral_per_cell=4,
    )
    spatial = analysis.validate_spatial_summaries(_synthetic_spatial(), classes)
    assert analysis._source_display_name(spatial) == "Gamma-SMC"

    paths = analysis.plot_spatial_probability_profiles(spatial, tmp_path)
    assert paths
    assert {path.suffix for path in paths} == {".png", ".pdf"}
    assert all(path.is_file() and path.stat().st_size > 0 for path in paths)
    png = next(path for path in paths if path.suffix == ".png")
    with png.open("rb") as stream:
        stream.seek(16)
        width, height = struct.unpack(">II", stream.read(8))
    assert (width, height) == (3_300, 2_550)

    bad = spatial.copy()
    bad.loc[bad.index[0], "n_units"] += 1
    with pytest.raises(ValueError, match="unit count"):
        analysis.validate_spatial_summaries(bad, classes)


def test_spatial_profiles_use_one_shared_x_label_for_three_af_panels(
    tmp_path, monkeypatch
):
    spatial = pd.concat(
        [
            _synthetic_spatial().assign(target_allele_frequency=af)
            for af in (0.1, 0.2, 0.3)
        ],
        ignore_index=True,
    )
    captured = []

    def capture_layout(figure, _path, **save_options):
        if not captured:
            captured.append(
                {
                    "figure_labels": [text.get_text() for text in figure.texts],
                    "axis_labels": [axis.get_xlabel() for axis in figure.axes],
                    "axes_count": len(figure.axes),
                    "subplot_bottom": figure.subplotpars.bottom,
                    "save_options": save_options,
                }
            )

    monkeypatch.setattr(analysis, "_atomic_figure", capture_layout)
    paths = analysis.plot_spatial_probability_profiles(spatial, tmp_path)

    assert paths
    assert captured[0]["axes_count"] == 3
    assert captured[0]["figure_labels"].count("Position in 10-Mb region (Mb)") == 1
    assert captured[0]["axis_labels"] == ["", "", ""]
    assert captured[0]["subplot_bottom"] == pytest.approx(0.28)
    assert captured[0]["save_options"] == {"bbox_inches": None}


def test_analysis_source_hash_is_portable_across_crlf_checkouts(tmp_path):
    lf = tmp_path / "lf.py"
    crlf = tmp_path / "crlf.py"
    lf.write_bytes(b"print('one')\nprint('two')\n")
    crlf.write_bytes(b"print('one')\r\nprint('two')\r\n")

    expected = hashlib.sha256(lf.read_bytes()).hexdigest()
    assert analysis._canonical_text_sha256(lf) == expected
    assert analysis._canonical_text_sha256(crlf) == expected


@pytest.mark.parametrize(
    "payload",
    (b"\xef\xbb\xbfprint('bom')\n", b"print('bare')\rnext\n", b"\xff\n"),
)
def test_analysis_source_hash_rejects_nonportable_text(tmp_path, payload):
    source = tmp_path / "source.py"
    source.write_bytes(payload)

    with pytest.raises(ValueError):
        analysis._canonical_text_sha256(source)


def test_analysis_implementation_paths_are_repo_relative():
    contract = analysis._analysis_implementation_contract()

    assert contract["path_identity"] == "repo-relative POSIX for in-repository paths"
    assert all(
        record["path"].startswith("python/gamma_smc_aou/")
        for record in contract["sources"].values()
    )


def test_empirical_schema_matches_exact_cell_and_uses_plus_one_null():
    classes = analysis.validate_class_summaries(
        _synthetic_classes(),
        minimum_selected_per_cell=2,
        minimum_neutral_per_cell=4,
    )
    scores = analysis.build_replicate_statistics(classes)
    template = analysis.empirical_input_template(classes)
    assert len(template) == 2 * 5 * len(THRESHOLDS)
    empirical = _complete_empirical_rows(template.iloc[[0]])
    empirical["observed_value"] = 0.30
    result = analysis.empirical_null_pvalues(empirical, scores)
    assert result.iloc[0]["n_null"] == 4
    assert result.iloc[0]["mc_p_upper"] == 0.2
    assert result.iloc[0]["bh_q_upper"] == 0.2
    assert result.iloc[0]["provenance_qc_status"] == "pass"
    assert result.iloc[0]["genotype_derived_alt_count"] == 20

    incomplete = _complete_empirical_rows(template.iloc[[0]])
    incomplete["observed_value"] = 0.30
    incomplete["source_data_sha256"] = "replace_with_sha256"
    with pytest.raises(ValueError, match="provenance field is unset"):
        analysis.empirical_null_pvalues(incomplete, scores)

    bad = empirical.copy()
    bad["demography_id"] = "unmatched"
    with pytest.raises(ValueError, match="no neutral reference settings"):
        analysis.empirical_null_pvalues(bad, scores)

    bad_counts = empirical.copy()
    bad_counts["alt_allele_count"] = 19
    with pytest.raises(ValueError, match="alt count disagrees"):
        analysis.empirical_null_pvalues(bad_counts, scores)

    bad_statistic = empirical.copy()
    bad_statistic["observed_value"] = 0.29
    with pytest.raises(ValueError, match="disagrees with genotype classes"):
        analysis.empirical_null_pvalues(bad_statistic, scores)

    bad_decoder = empirical.copy()
    bad_decoder["decoder_name"] = "different_decoder"
    with pytest.raises(ValueError, match="decoder name does not match"):
        analysis.empirical_null_pvalues(bad_decoder, scores)

    mismatched_settings = empirical.copy()
    mismatched_settings["empirical_decoder_parameters_sha256"] = "8" * 64
    with pytest.raises(ValueError, match="does not match its neutral reference"):
        analysis.empirical_null_pvalues(mismatched_settings, scores)

    inconsistent_rows = _complete_empirical_rows(template.iloc[:2])
    inconsistent_rows["observed_value"] = 0.30
    inconsistent_rows.loc[inconsistent_rows.index[1], "pair_table_sha256"] = "9" * 64
    with pytest.raises(ValueError, match="provenance changes across statistic rows"):
        analysis.empirical_null_pvalues(inconsistent_rows, scores)

    nonmonotone = _complete_empirical_rows(template.iloc[:2])
    nonmonotone.loc[nonmonotone.index[1], "observed_p_hom_alt"] = 0.30
    nonmonotone.loc[nonmonotone.index[1], "observed_value"] = 0.20
    with pytest.raises(ValueError, match="not nondecreasing"):
        analysis.empirical_null_pvalues(nonmonotone, scores)


def test_within_selected_unit_pair_auc_and_exact_label_permutation():
    result = analysis.within_selected_unit_pair_tests(
        _synthetic_pairs(), permutation_draws=19, base_seed=3
    )
    assert len(result) == 2 * 2 * (len(THRESHOLDS) + 2)
    primary = result[
        (result["statistic_scope"] == "focal_nearest")
        & (result["metric"] == "p_tmrca_lt_threshold")
        & (result["threshold_years"] == 1_000)
    ]
    assert set(primary["within_pair_roc_auc"]) == {1.0}
    assert set(primary["permutation_method"]) == {"exact_all_label_assignments"}
    assert set(primary["permutation_draws_or_assignments"]) == {6}
    np.testing.assert_allclose(primary["genotype_label_permutation_p_upper"], 1 / 6)


def test_analysis_writes_letter_plots_both_formats_and_validates_cache(
    tmp_path, monkeypatch
):
    input_path = tmp_path / "combined_class_summaries.tsv.gz"
    _synthetic_classes().to_csv(input_path, sep="\t", index=False, compression="gzip")
    pair_path = tmp_path / "combined_pair_summaries.tsv.gz"
    _synthetic_pairs().to_csv(pair_path, sep="\t", index=False, compression="gzip")
    spatial_path = tmp_path / "aggregated_spatial_profiles.tsv.gz"
    _synthetic_spatial().to_csv(spatial_path, sep="\t", index=False, compression="gzip")
    output = tmp_path / "analysis"
    paths = analysis.run_focused_analysis(
        input_path,
        output,
        pair_summaries_path=pair_path,
        spatial_summaries_path=spatial_path,
        minimum_selected_per_cell=2,
        minimum_neutral_per_cell=4,
        bootstrap_draws=20,
        permutation_draws=19,
        crossfit_folds=2,
        base_seed=9,
        make_plots=True,
    )
    assert paths["completion"].is_file()
    completion = json.loads(paths["completion"].read_text(encoding="utf-8"))
    implementation = completion["contract"]["implementation"]
    assert set(implementation["sources"]) == {
        path.name for path in analysis.ANALYSIS_IMPLEMENTATION_SOURCES
    }
    assert set(implementation["software"]) == {
        "python",
        "matplotlib",
        "numpy",
        "pandas",
        "scipy",
    }
    assert all(
        len(record["sha256"]) == 64 for record in implementation["sources"].values()
    )
    png = sorted((output / "figures").glob("*.png"))
    pdf = sorted((output / "figures").glob("*.pdf"))
    assert png and len(png) == len(pdf)
    table_mtime = paths["cell_comparisons"].stat().st_mtime_ns
    cached = analysis.run_focused_analysis(
        input_path,
        output,
        pair_summaries_path=pair_path,
        spatial_summaries_path=spatial_path,
        minimum_selected_per_cell=2,
        minimum_neutral_per_cell=4,
        bootstrap_draws=20,
        permutation_draws=19,
        crossfit_folds=2,
        base_seed=9,
        make_plots=True,
    )
    assert cached["cell_comparisons"].stat().st_mtime_ns == table_mtime

    implementation_contract = analysis._analysis_implementation_contract
    changed = implementation_contract()
    changed["software"] = {**changed["software"], "python": "changed-for-test"}
    monkeypatch.setattr(analysis, "_analysis_implementation_contract", lambda: changed)
    with pytest.raises(ValueError, match="contract is incompatible"):
        analysis.run_focused_analysis(
            input_path,
            output,
            pair_summaries_path=pair_path,
            spatial_summaries_path=spatial_path,
            minimum_selected_per_cell=2,
            minimum_neutral_per_cell=4,
            bootstrap_draws=20,
            permutation_draws=19,
            crossfit_folds=2,
            base_seed=9,
            make_plots=True,
        )
    monkeypatch.setattr(
        analysis, "_analysis_implementation_contract", implementation_contract
    )

    with pytest.raises(ValueError, match="contract is incompatible"):
        analysis.run_focused_analysis(
            input_path,
            output,
            pair_summaries_path=pair_path,
            spatial_summaries_path=spatial_path,
            minimum_selected_per_cell=2,
            minimum_neutral_per_cell=4,
            bootstrap_draws=20,
            permutation_draws=23,
            crossfit_folds=2,
            base_seed=9,
            make_plots=True,
        )

    paths["cell_comparisons"].write_text("corrupt\n", encoding="utf-8")
    with pytest.raises(ValueError, match="checksum"):
        analysis.run_focused_analysis(
            input_path,
            output,
            pair_summaries_path=pair_path,
            spatial_summaries_path=spatial_path,
            minimum_selected_per_cell=2,
            minimum_neutral_per_cell=4,
            bootstrap_draws=20,
            permutation_draws=19,
            crossfit_folds=2,
            base_seed=9,
            make_plots=True,
        )


def test_tree_truth_analysis_without_pair_table_is_cache_idempotent(tmp_path):
    input_path = tmp_path / "combined_truth_class_summaries.tsv.gz"
    classes = _synthetic_classes().copy()
    classes["source"] = "tree_truth"
    classes.to_csv(input_path, sep="\t", index=False, compression="gzip")
    spatial_path = tmp_path / "aggregated_truth_spatial_profiles.tsv.gz"
    _synthetic_spatial("tree_truth").to_csv(
        spatial_path, sep="\t", index=False, compression="gzip"
    )
    output = tmp_path / "analysis_tree_truth"
    first = analysis.run_focused_analysis(
        input_path,
        output,
        spatial_summaries_path=spatial_path,
        minimum_selected_per_cell=2,
        minimum_neutral_per_cell=4,
        bootstrap_draws=20,
        permutation_draws=19,
        crossfit_folds=2,
        base_seed=9,
        make_plots=False,
    )
    table_mtime = first["unit_pair_tests"].stat().st_mtime_ns
    second = analysis.run_focused_analysis(
        input_path,
        output,
        spatial_summaries_path=spatial_path,
        minimum_selected_per_cell=2,
        minimum_neutral_per_cell=4,
        bootstrap_draws=20,
        permutation_draws=19,
        crossfit_folds=2,
        base_seed=9,
        make_plots=False,
    )

    assert second["unit_pair_tests"].stat().st_mtime_ns == table_mtime


def test_analysis_rejects_incomplete_cell_counts_and_duplicate_rows():
    frame = _synthetic_classes()
    with pytest.raises(ValueError, match="at least 100"):
        analysis.validate_class_summaries(frame, minimum_selected_per_cell=2)
    duplicate = pd.concat([frame, frame.iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError, match="duplicate"):
        analysis.validate_class_summaries(
            duplicate,
            minimum_selected_per_cell=2,
            minimum_neutral_per_cell=4,
        )


def test_truth_gamma_comparison_is_checksummed_and_fail_closed(tmp_path, monkeypatch):
    truth_path = tmp_path / "truth.tsv.gz"
    gamma_path = tmp_path / "gamma.tsv.gz"
    truth = _synthetic_classes()
    gamma = truth.copy()
    gamma["focal_mean_p_tmrca_lt_threshold"] = np.minimum(
        1.0, gamma["focal_mean_p_tmrca_lt_threshold"] + 0.01
    )
    truth.to_csv(truth_path, sep="\t", index=False, compression="gzip")
    gamma.to_csv(gamma_path, sep="\t", index=False, compression="gzip")
    output = tmp_path / "truth_gamma"
    paths = analysis.run_truth_gamma_comparison(
        truth_path,
        gamma_path,
        output,
        minimum_selected_per_cell=2,
        minimum_neutral_per_cell=4,
    )
    assert paths["summary"].is_file()
    original_mtime = paths["summary"].stat().st_mtime_ns
    cached = analysis.run_truth_gamma_comparison(
        truth_path,
        gamma_path,
        output,
        minimum_selected_per_cell=2,
        minimum_neutral_per_cell=4,
    )
    assert cached["summary"].stat().st_mtime_ns == original_mtime
    completion = json.loads(paths["completion"].read_text(encoding="utf-8"))
    assert completion["contract"]["implementation"]["sources"]
    assert completion["contract"]["implementation"]["software"]

    implementation_contract = analysis._analysis_implementation_contract
    changed = implementation_contract()
    changed["software"] = {**changed["software"], "scipy": "changed-for-test"}
    monkeypatch.setattr(analysis, "_analysis_implementation_contract", lambda: changed)
    with pytest.raises(ValueError, match="contract is incompatible"):
        analysis.run_truth_gamma_comparison(
            truth_path,
            gamma_path,
            output,
            minimum_selected_per_cell=2,
            minimum_neutral_per_cell=4,
        )
    monkeypatch.setattr(
        analysis, "_analysis_implementation_contract", implementation_contract
    )

    paths["summary"].write_text("corrupt\n", encoding="utf-8")
    with pytest.raises(ValueError, match="checksum"):
        analysis.run_truth_gamma_comparison(
            truth_path,
            gamma_path,
            output,
            minimum_selected_per_cell=2,
            minimum_neutral_per_cell=4,
        )
