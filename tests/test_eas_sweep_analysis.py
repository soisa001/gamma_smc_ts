import numpy as np
import pandas as pd
import pytest
import tskit

from gamma_smc_aou.eas_sweep_analysis import (
    build_diploid_pair_table,
    focal_genotype_counts,
    gamma_pair_tables,
    gamma_summaries_to_long,
    gamma_summary_to_long,
    internal_genotype_contrasts,
    normalized_cdf_auc,
    plot_acceptance_effort,
    plot_ancient_eurasia_introgression_timeline,
    plot_center_internal_contrast_heatmaps,
    plot_eas_ne_origin_design,
    plot_simulation_status_grid,
    tree_truth_profiles,
    write_cross_spec_study_plots,
    write_gamma_pair_files,
    write_publication_plots,
)


def tiny_diploid_tree_sequence():
    """Three diploids with known within-pair TMRCA in two spatial trees."""
    tables = tskit.TableCollection(sequence_length=100)
    individuals = [tables.individuals.add_row() for _ in range(3)]
    samples = [
        tables.nodes.add_row(
            flags=tskit.NODE_IS_SAMPLE,
            time=0,
            individual=individuals[index // 2],
        )
        for index in range(6)
    ]

    first_parents = [tables.nodes.add_row(time=time) for time in (100, 20, 10)]
    first_root = tables.nodes.add_row(time=200)
    second_parents = [tables.nodes.add_row(time=time) for time in (80, 25, 5)]
    second_root = tables.nodes.add_row(time=200)
    for interval, parents, root in (
        ((0, 50), first_parents, first_root),
        ((50, 100), second_parents, second_root),
    ):
        left, right = interval
        for individual_index, parent in enumerate(parents):
            for child in samples[2 * individual_index : 2 * individual_index + 2]:
                tables.edges.add_row(left, right, parent, child)
            tables.edges.add_row(left, right, root, parent)

    site = tables.sites.add_row(position=25, ancestral_state="A")
    for node in (samples[2], samples[4], samples[5]):
        tables.mutations.add_row(site=site, node=node, derived_state="G")
    tables.sort()
    return tables.tree_sequence()


def test_focal_genotypes_map_deterministically_to_gamma_pair_files(tmp_path):
    ts = tiny_diploid_tree_sequence()
    np.testing.assert_array_equal(focal_genotype_counts(ts, 25), [0, 1, 2])

    table = build_diploid_pair_table(ts, 25)
    assert table["tree_sequence_individual_id"].tolist() == [0, 1, 2]
    assert table["gamma_smc_haplotype_0"].tolist() == [0, 2, 4]
    assert table["gamma_smc_haplotype_1"].tolist() == [1, 3, 5]
    assert table["genotype_class"].tolist() == [
        "hom_ref",
        "heterozygous",
        "hom_alt",
    ]

    pair_tables = gamma_pair_tables(table)
    assert pair_tables["overall"].to_numpy().tolist() == [[0, 1], [2, 3], [4, 5]]
    assert pair_tables["hom_alt"].to_numpy().tolist() == [[4, 5]]

    paths = write_gamma_pair_files(table, tmp_path, prefix="tiny_")
    assert set(paths) == {"overall", "hom_ref", "heterozygous", "hom_alt"}
    written = pd.read_csv(
        paths["overall"], sep="\t", comment="#", header=None, dtype=int
    )
    assert written.to_numpy().tolist() == [[0, 1], [2, 3], [4, 5]]
    assert "# genotype_class\toverall" in paths["overall"].read_text(encoding="utf-8")


def test_tree_truth_profiles_use_declared_year_to_generation_conversion():
    ts = tiny_diploid_tree_sequence()
    pair_table = build_diploid_pair_table(ts, 25, [0, 1, 2])
    profile = tree_truth_profiles(
        ts,
        [25, 75],
        pair_table,
        thresholds_years=[1_000, 4_500],
        generation_time_years=25,
    )
    assert len(profile) == 2 * 4 * 2

    center = profile[profile["position_0based"] == 25]
    overall_1k = center[
        (center["genotype_class"] == "overall") & (center["threshold_years"] == 1_000)
    ].iloc[0]
    assert overall_1k["threshold_generations"] == 40
    assert overall_1k["p_tmrca_lt_threshold"] == pytest.approx(2 / 3)

    ref = center[center["genotype_class"] == "hom_ref"]
    alt = center[center["genotype_class"] == "hom_alt"]
    assert ref["mean_tmrca_generations"].unique().tolist() == [100]
    assert alt["mean_tmrca_generations"].unique().tolist() == [10]
    assert alt.sort_values("threshold_years")["p_tmrca_lt_threshold"].tolist() == [1, 1]
    assert ref.sort_values("threshold_years")["p_tmrca_lt_threshold"].tolist() == [0, 1]


def test_normalized_cdf_auc_includes_origin_and_contrast_signs_are_explicit():
    ts = tiny_diploid_tree_sequence()
    profile = tree_truth_profiles(
        ts,
        [25],
        build_diploid_pair_table(ts, 25),
        thresholds_years=[1_000, 4_500],
    )
    auc = normalized_cdf_auc(profile)
    ref_auc = auc.loc[auc["genotype_class"] == "hom_ref"].iloc[0]
    alt_auc = auc.loc[auc["genotype_class"] == "hom_alt"].iloc[0]
    assert bool(ref_auc["cdf_auc_includes_origin_0_0"])
    assert ref_auc["normalized_cdf_auc"] == pytest.approx(1_750 / 4_500)
    assert alt_auc["normalized_cdf_auc"] == pytest.approx(4_000 / 4_500)

    contrast = internal_genotype_contrasts(profile).iloc[0]
    assert contrast["contrast_scope"] == "internal_genotype_contrast_no_neutral_null"
    assert contrast["alt_minus_ref_normalized_cdf_auc"] == pytest.approx(0.5)
    assert contrast["alt_minus_ref_mean_tmrca_generations"] == -90
    assert contrast["alt_minus_ref_mean_tmrca_years"] == -2_250
    assert contrast["expected_sweep_sign_normalized_cdf_auc"] == "positive"
    assert contrast["expected_sweep_sign_mean_tmrca"] == "negative"


def _gamma_frame(*, mean_tmrca, p_1k, p_4k5, recent_1k, recent_4k5):
    return pd.DataFrame(
        {
            "position_0based": [25],
            "position_1based": [26],
            "n_pairs": [10],
            "mean_tmrca_generations": [mean_tmrca],
            "mean_p_lt_1000": [p_1k],
            "mean_p_lt_4500": [p_4k5],
            "frac_recent_1000": [recent_1k],
            "frac_recent_4500": [recent_4k5],
        }
    )


def test_gamma_multi_threshold_summary_is_long_without_conflating_call_fraction():
    frame = _gamma_frame(
        mean_tmrca=100,
        p_1k=0.1,
        p_4k5=0.5,
        recent_1k=0.0,
        recent_4k5=0.4,
    )
    long = gamma_summary_to_long(frame, genotype_class="hom_ref")
    assert long["threshold_years"].tolist() == [1_000, 4_500]
    assert long["threshold_generations"].tolist() == [40, 180]
    assert long["p_tmrca_lt_threshold"].tolist() == [0.1, 0.5]
    assert long["frac_recent"].tolist() == [0.0, 0.4]
    assert long["mean_tmrca_years"].tolist() == [2_500, 2_500]

    combined = gamma_summaries_to_long(
        {
            "hom_ref": frame,
            "hom_alt": _gamma_frame(
                mean_tmrca=20,
                p_1k=0.5,
                p_4k5=0.9,
                recent_1k=0.4,
                recent_4k5=0.8,
            ),
        }
    )
    contrast = internal_genotype_contrasts(combined).iloc[0]
    assert contrast["alt_minus_ref_normalized_cdf_auc"] > 0
    assert contrast["alt_minus_ref_mean_tmrca_generations"] == -80


def test_nonmonotone_cdf_is_rejected():
    frame = _gamma_frame(
        mean_tmrca=100,
        p_1k=0.8,
        p_4k5=0.2,
        recent_1k=0.8,
        recent_4k5=0.2,
    )
    long = gamma_summary_to_long(frame, genotype_class="hom_ref")
    with pytest.raises(ValueError, match="nondecreasing"):
        normalized_cdf_auc(long)


def test_publication_plot_bundle_writes_png_and_pdf(tmp_path):
    ts = tiny_diploid_tree_sequence()
    profile = tree_truth_profiles(
        ts,
        [25, 75],
        build_diploid_pair_table(ts, 25),
        thresholds_years=[1_000, 4_500],
    )
    outputs = write_publication_plots(
        profile,
        focal_coordinate=25,
        output_dir=tmp_path,
        prefix="tiny",
    )
    assert set(outputs) == {
        "center_cdf",
        "tmrca_profiles",
        "spatial_heatmaps_tree_truth",
    }
    for formats in outputs.values():
        assert set(formats) == {"png", "pdf"}
        for suffix, path in formats.items():
            assert path.suffix == f".{suffix}"
            assert path.is_file()
            assert path.stat().st_size > 1_000


def _compact_cross_spec_tables():
    acceptance = pd.DataFrame(
        [
            {
                "study_type": "no_introgression",
                "specification_id": "no_lower_s01",
                "selection_coefficient": 0.01,
                "demography_curve": "lower_ci",
                "target_allele_frequency": 0.50,
                "achieved_allele_frequency": 0.48,
                "attempts_to_accept": 3,
            },
            {
                "study_type": "no_introgression",
                "specification_id": "no_median_s005",
                "selection_coefficient": 0.005,
                "demography_curve": "median",
                "target_allele_frequency": 0.50,
                "achieved_allele_frequency": 0.52,
                "attempts_to_accept": 17,
            },
            {
                "study_type": "introgression",
                "specification_id": "intro_af10",
                "selection_coefficient": 0.01,
                "demography_curve": np.nan,
                "target_allele_frequency": 0.10,
                "achieved_allele_frequency": 0.12,
                "attempts_to_accept": 2,
            },
            {
                "study_type": "introgression",
                "specification_id": "intro_af30",
                "selection_coefficient": 0.005,
                "demography_curve": np.nan,
                "target_allele_frequency": 0.30,
                "achieved_allele_frequency": 0.31,
                "attempts_to_accept": 43,
            },
        ]
    )
    center_rows = []
    for row in acceptance.to_dict(orient="records"):
        for source, scale in (("tree_truth", 1.0), ("gamma_smc", 0.8)):
            center_rows.append(
                {
                    **{
                        key: row[key]
                        for key in (
                            "study_type",
                            "specification_id",
                            "selection_coefficient",
                            "demography_curve",
                            "target_allele_frequency",
                        )
                    },
                    "source": source,
                    "alt_minus_ref_normalized_cdf_auc": 0.2 * scale,
                    "alt_minus_ref_mean_tmrca_generations": -80 * scale,
                }
            )
    return pd.DataFrame(center_rows), acceptance


def test_cross_spec_plot_bundle_supports_study_type_fallback(tmp_path):
    center, acceptance = _compact_cross_spec_tables()
    outputs = write_cross_spec_study_plots(
        center,
        acceptance,
        tmp_path,
        prefix="compact",
    )
    assert set(outputs) == {"introgression", "no_introgression"}
    for scenario_outputs in outputs.values():
        assert set(scenario_outputs) == {
            "acceptance_effort",
            "center_internal_contrasts",
        }
        for formats in scenario_outputs.values():
            assert set(formats) == {"png", "pdf"}
            for path in formats.values():
                assert path.is_file()
                assert path.stat().st_size > 1_000


def test_cross_spec_plot_primitives_reject_ambiguous_rows(tmp_path):
    center, acceptance = _compact_cross_spec_tables()
    no_intro_acceptance = acceptance[acceptance["study_type"] == "no_introgression"]
    acceptance_paths = plot_acceptance_effort(
        no_intro_acceptance,
        tmp_path / "acceptance",
    )
    assert all(path.is_file() for path in acceptance_paths.values())

    no_intro_center = center[center["study_type"] == "no_introgression"]
    contrast_paths = plot_center_internal_contrast_heatmaps(
        no_intro_center,
        tmp_path / "contrasts",
    )
    assert all(path.is_file() for path in contrast_paths.values())

    duplicate = pd.concat([no_intro_center, no_intro_center.iloc[[0]]])
    with pytest.raises(ValueError, match="one row per specification/source"):
        plot_center_internal_contrast_heatmaps(
            duplicate,
            tmp_path / "duplicate",
        )


def test_cross_spec_bundle_rejects_mismatched_specification_sets(tmp_path):
    center, acceptance = _compact_cross_spec_tables()
    acceptance = acceptance[acceptance["specification_id"] != "intro_af30"]
    with pytest.raises(
        ValueError, match="different center and acceptance specifications"
    ):
        write_cross_spec_study_plots(center, acceptance, tmp_path)


def _tidy_eas_design_tables():
    times = [0, 100, 1_000, 10_000, 40_000]
    sizes = {
        "q025": [30_000, 30_000, 8_000, 12_000, 18_000],
        "median": [36_000, 36_000, 10_000, 15_000, 22_000],
        "q975": [44_000, 44_000, 13_000, 19_000, 28_000],
    }
    trajectories = pd.DataFrame(
        [
            {
                "demography_curve": curve,
                "time_generations": time,
                "effective_population_size": ne,
            }
            for curve, curve_sizes in sizes.items()
            for time, ne in zip(times, curve_sizes, strict=True)
        ]
    )
    ages = {
        0.01: [1_748, 1_764, 1_782],
        0.005: [3_697, 3_768, 3_838],
        0.001: [20_190, 20_391, 20_677],
    }
    curves = ["q025", "median", "q975"]
    origins = pd.DataFrame(
        [
            {
                "demography_curve": curve,
                "selection_coefficient": coefficient,
                "origin_age_generations": age,
            }
            for coefficient, coefficient_ages in ages.items()
            for curve, age in zip(curves, coefficient_ages, strict=True)
        ]
    )
    return trajectories, origins


def test_eas_ne_origin_design_writes_letter_png_and_pdf(tmp_path):
    trajectories, origins = _tidy_eas_design_tables()
    paths = plot_eas_ne_origin_design(
        trajectories,
        origins,
        tmp_path / "eas_ne_origins",
    )
    assert set(paths) == {"png", "pdf"}
    for path in paths.values():
        assert path.is_file()
        assert path.stat().st_size > 1_000

    with pytest.raises(ValueError, match="nine curve-by-selection rows"):
        plot_eas_ne_origin_design(
            trajectories,
            origins.iloc[:-1],
            tmp_path / "incomplete_origins",
        )


def test_ancient_eurasia_timeline_writes_png_pdf_and_checks_event_order(tmp_path):
    paths = plot_ancient_eurasia_introgression_timeline(
        tmp_path / "ancient_eurasia_timeline",
        realized_introgression_pulse_generations=2_270,
        realized_han_split_generations=2_010,
    )
    assert set(paths) == {"png", "pdf"}
    for path in paths.values():
        assert path.is_file()
        assert path.stat().st_size > 1_000

    with pytest.raises(ValueError, match="split > mutation origin > introgression"):
        plot_ancient_eurasia_introgression_timeline(
            tmp_path / "invalid_timeline",
            mutation_origin_generations=2_200,
        )

    with pytest.raises(ValueError, match="realized Han split"):
        plot_ancient_eurasia_introgression_timeline(
            tmp_path / "invalid_realized_split",
            realized_han_split_generations=2_300,
        )

    with pytest.raises(ValueError, match="realized introgression pulse"):
        plot_ancient_eurasia_introgression_timeline(
            tmp_path / "invalid_realized_pulse",
            realized_introgression_pulse_generations=2_500,
        )


def _simulation_status_frame(mode):
    selections = [0.01, 0.005, 0.001]
    statuses = ["complete", "cached", "timed_out", "failed", "not_requested"]
    rows = []
    if mode == "no_introgression":
        columns = [(curve, 0.5) for curve in ("q025", "median", "q975")]
    else:
        columns = [(np.nan, target / 10) for target in range(1, 10)]
    for selection_index, selection in enumerate(selections):
        for column_index, (curve, target) in enumerate(columns):
            index = selection_index * len(columns) + column_index
            status = statuses[index % len(statuses)]
            draws = index + 1 if status != "not_requested" else 0
            rows.append(
                {
                    "specification_id": f"{mode}_{selection_index}_{column_index}",
                    "selection_coefficient": selection,
                    "demography_curve": curve,
                    "target_allele_frequency": target,
                    "status": status,
                    "elapsed_seconds": (
                        float((index + 1) * 75) if status != "not_requested" else 0.0
                    ),
                    "completed_external_draws": draws,
                    "last_completed_draw_sample_af": (
                        min(0.99, target + 0.01) if draws else np.nan
                    ),
                }
            )
    return pd.DataFrame(rows)


@pytest.mark.parametrize("mode", ["no_introgression", "introgression"])
def test_simulation_status_grid_writes_all_design_cells(mode, tmp_path):
    frame = _simulation_status_frame(mode)
    paths = plot_simulation_status_grid(
        frame,
        tmp_path / f"{mode}_status_grid",
        scenario_label=mode.replace("_", " ").title(),
    )
    assert set(paths) == {"png", "pdf"}
    for path in paths.values():
        assert path.is_file()
        assert path.stat().st_size > 1_000


def test_simulation_status_grid_rejects_missing_design_cell(tmp_path):
    incomplete = _simulation_status_frame("no_introgression").iloc[:-1]
    with pytest.raises(ValueError, match="design grid is incomplete"):
        plot_simulation_status_grid(incomplete, tmp_path / "incomplete_status")
