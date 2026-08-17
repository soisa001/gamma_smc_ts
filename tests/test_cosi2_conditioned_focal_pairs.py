import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import zstandard
from gamma_smc_aou import cosi2_conditioned_focal_pairs as focal
from matplotlib.backends.backend_agg import FigureCanvasAgg
from scipy.special import gammainc


def _study_frames() -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, object]] = []
    eligibility: list[dict[str, object]] = []
    thresholds = np.asarray(focal.THRESHOLDS_YEARS)
    base = np.arange(1, len(thresholds) + 1, dtype=float) * 0.05

    units = [
        *((f"eas_n{index}", "EAS", "neutral", True, index) for index in range(4)),
        ("eas_selected", "EAS", "selected", True, 10),
        ("han_n0", "Han", "neutral", True, 0),
        ("han_n1", "Han", "neutral", True, 1),
        ("han_ref_only0", "Han", "neutral", False, 2),
        ("han_ref_only1", "Han", "neutral", False, 3),
        ("han_selected", "Han", "selected", True, 10),
    ]
    for unit_id, demography, simulation_class, has_alt, offset in units:
        is_selected = simulation_class == "selected"
        ref = base if is_selected else base + 0.01 * offset
        delta_selected = 0.15 if demography == "EAS" else 0.10
        delta = (
            np.full(len(base), delta_selected)
            if is_selected
            else np.full(len(base), 0.02 + 0.01 * offset)
        )
        curves = {"p_hom_ref": ref}
        if has_alt:
            curves["p_hom_alt"] = ref + delta
            curves["p_hom_alt_minus_hom_ref"] = delta
        for statistic, curve in curves.items():
            if statistic == "p_hom_ref":
                count_fields = {
                    "n_pairs": 100,
                    "n_valid_pairs": 100,
                    "n_ref_pairs": 100,
                    "n_alt_pairs": 0,
                    "n_valid_ref_pairs": 100,
                    "n_valid_alt_pairs": 0,
                }
            elif statistic == "p_hom_alt":
                count_fields = {
                    "n_pairs": 10,
                    "n_valid_pairs": 10,
                    "n_ref_pairs": 0,
                    "n_alt_pairs": 10,
                    "n_valid_ref_pairs": 0,
                    "n_valid_alt_pairs": 10,
                }
            else:
                count_fields = {
                    "n_pairs": -1,
                    "n_valid_pairs": -1,
                    "n_ref_pairs": 100,
                    "n_alt_pairs": 10,
                    "n_valid_ref_pairs": 100,
                    "n_valid_alt_pairs": 10,
                }
            for threshold, score in zip(thresholds, curve, strict=True):
                rows.append(
                    {
                        "schema": focal.SCHEMA_VERSION,
                        "unit_id": unit_id,
                        "demography": demography,
                        "simulation_class": simulation_class,
                        "selection_coefficient": 0.01 if is_selected else 0.0,
                        "seed": offset,
                        "generation_time_years": 25.0,
                        "source": "synthetic_float32_posterior",
                        "pairing_unit": "unordered_haplotype_pair_pseudodiploid",
                        "statistic": statistic,
                        "threshold_years": int(threshold),
                        "score": float(score),
                        **count_fields,
                    }
                )
        eligibility.append(
            {
                "unit_id": unit_id,
                "demography": demography,
                "simulation_class": simulation_class,
                "selection_coefficient": 0.01 if is_selected else 0.0,
                "seed": offset,
                "generation_time_years": 25.0,
                "focal_status": (
                    "unique_exact_focal_site" if has_alt else "absent_from_sample"
                ),
                "focal_alt_count_full_sample": 20 if has_alt else 0,
                "focal_alt_count_pair_panel": 10 if has_alt else 0,
                "n_hom_ref_pairs": 100,
                "n_heterozygous_pairs": 20 if has_alt else 0,
                "n_hom_alt_pairs": 10 if has_alt else 0,
                "hom_ref_eligible": True,
                "heterozygous_eligible": has_alt,
                "hom_alt_eligible": has_alt,
                "primary_shared_homozygote_eligible": has_alt,
                "exclusion_reason": (
                    ""
                    if has_alt
                    else "fewer_than_two_alt_haplotypes_in_fixed_pair_panel"
                ),
            }
        )
    eligibility.append(
        {
            "unit_id": "eas_ambiguous",
            "demography": "EAS",
            "simulation_class": "neutral",
            "selection_coefficient": 0.0,
            "seed": 99,
            "generation_time_years": 25.0,
            "focal_status": "ambiguous_exact_coordinate_collision",
            "focal_alt_count_full_sample": -1,
            "focal_alt_count_pair_panel": -1,
            "n_hom_ref_pairs": 0,
            "n_heterozygous_pairs": 0,
            "n_hom_alt_pairs": 0,
            "hom_ref_eligible": False,
            "heterozygous_eligible": False,
            "hom_alt_eligible": False,
            "primary_shared_homozygote_eligible": False,
            "exclusion_reason": "ambiguous_exact_coordinate_collision",
        }
    )
    return pd.DataFrame(rows), pd.DataFrame(eligibility)


def _reconciliation(eligibility: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    eligible = eligibility[
        eligibility["focal_status"] != "ambiguous_exact_coordinate_collision"
    ]
    for unit in eligible.to_dict(orient="records"):
        for index, threshold in enumerate(focal.THRESHOLDS_YEARS):
            decoder = 0.1 + index * 0.01
            signed = 1e-6
            rows.append(
                {
                    "schema": focal.SCHEMA_VERSION,
                    "unit_id": unit["unit_id"],
                    "demography": unit["demography"],
                    "simulation_class": unit["simulation_class"],
                    "selection_coefficient": unit["selection_coefficient"],
                    "seed": unit["seed"],
                    "generation_time_years": unit["generation_time_years"],
                    "threshold_years": threshold,
                    "posthoc_overall_score": decoder + signed,
                    "decoder_cpp_overall_score": decoder,
                    "signed_difference": signed,
                    "absolute_difference": abs(signed),
                }
            )
    return pd.DataFrame(rows)


def _write_focal_vcf(path: Path, *, duplicate_target: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    all_ref = ["0|0"] * 100
    target = ["1|1", "1|1", *(["0|0"] * 98)]
    records = [
        ["1", "5000001", "cosi_s00000006", "A", "G", ".", "PASS", ".", "GT", *all_ref],
        ["1", "5000001", "cosi_s00000007", "A", "G", ".", "PASS", ".", "GT", *target],
    ]
    if duplicate_target:
        records.append(records[-1])
    text = (
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n"
    )
    text += "\n".join("\t".join(record) for record in records) + "\n"
    path.write_text(text, encoding="ascii")


def test_exact_focal_vcf_id_wins_over_same_position_decoy(tmp_path):
    unit = tmp_path / "unit"
    vcf = unit / "converted" / "input.vcf"
    _write_focal_vcf(vcf)
    metadata = {
        "exact_focal_site_count": 1,
        "exact_focal_site_index": 7,
        "exact_focal_sample_alt_count": 4,
    }

    alleles, status = focal._focal_alleles(
        unit, metadata, recorded_status="unique_exact_focal_site"
    )

    assert status == "unique_exact_focal_site"
    assert alleles is not None
    assert int(alleles.sum()) == 4
    _write_focal_vcf(vcf, duplicate_target=True)
    with pytest.raises(ValueError, match="record ID is not unique"):
        focal._focal_alleles(unit, metadata, recorded_status="unique_exact_focal_site")


def test_pair_order_requires_exact_panel_combinations():
    panel = np.arange(100, dtype=np.int64)
    pairs = np.asarray(list(itertools.combinations(panel.tolist(), 2)), dtype=np.int64)
    focal._validate_fixed_pair_order(pairs, panel, unit_id="synthetic")

    reordered = pairs.copy()
    reordered[[0, 1]] = reordered[[1, 0]]
    duplicated = pairs.copy()
    duplicated[1] = duplicated[0]
    self_pair = pairs.copy()
    self_pair[0] = [0, 0]
    for corrupt in (reordered, duplicated, self_pair):
        with pytest.raises(ValueError, match="all-pairs panel order"):
            focal._validate_fixed_pair_order(corrupt, panel, unit_id="synthetic")


def test_focal_alt_counts_bind_vcf_conversion_and_study_metadata():
    alleles = np.zeros(200, dtype=np.int8)
    alleles[[0, 2, 101, 103]] = 1
    panel = np.arange(100, dtype=np.int64)
    conversion = {
        "exact_focal_sample_alt_count": 4,
        "exact_focal_pair_panel_alt_count": 2,
    }
    record = {
        "decode_full_sample_alt_count": 4,
        "decode_pair_panel_alt_count": 2,
    }
    assert focal._validate_focal_alt_counts(
        alleles,
        panel,
        status="unique_exact_focal_site",
        conversion_metadata=conversion,
        unit_record=record,
        unit_id="unique",
    ) == (4, 2)

    with pytest.raises(ValueError, match="recorded full-sample"):
        focal._validate_focal_alt_counts(
            np.zeros(200, dtype=np.int8),
            panel,
            status="absent_from_sample",
            conversion_metadata={
                "exact_focal_sample_alt_count": None,
                "exact_focal_pair_panel_alt_count": None,
            },
            unit_record={
                "decode_full_sample_alt_count": 1,
                "decode_pair_panel_alt_count": 0,
            },
            unit_id="absent",
        )


def test_streaming_focal_reducer_reads_only_focal_row(tmp_path):
    alpha0 = np.asarray([[99, 99], [2, 3]], dtype=np.float32)
    beta0 = np.asarray([[99, 99], [1, 1]], dtype=np.float32)
    alpha1 = np.asarray([[99, 99], [4, 0]], dtype=np.float32)
    beta1 = np.asarray([[99, 99], [2, 0]], dtype=np.float32)
    raw = b"".join(array.tobytes(order="C") for array in (alpha0, beta0, alpha1, beta1))
    posterior = tmp_path / "posterior.zst"
    posterior.write_bytes(zstandard.ZstdCompressor().compress(raw))
    metadata = {
        "num_pairs": 3,
        "chunk_size": 2,
        "sequence_length": 2,
        "output_positions": [0, focal.FOCAL_POSITION_0BASED],
        "scaled_mutation_rate": 2 * focal.MUTATION_RATE * 10_000,
    }
    pair_classes = np.asarray([0, 2, 0], dtype=np.int8)

    means, counts, overall = focal._stream_focal_scores(
        posterior,
        metadata,
        pair_classes,
        generation_time_years=1.0,
    )

    scaled = np.asarray(focal.THRESHOLDS_YEARS, dtype=float) / 10_000
    ref_expected = (gammainc(2.0, scaled) + gammainc(4.0, 2.0 * scaled)) / 2
    alt_expected = gammainc(3.0, scaled)
    np.testing.assert_allclose(means["p_hom_ref"], ref_expected, rtol=0, atol=1e-14)
    np.testing.assert_allclose(means["p_hom_alt"], alt_expected, rtol=0, atol=1e-14)
    np.testing.assert_allclose(
        overall, (2 * ref_expected + alt_expected) / 3, rtol=0, atol=1e-14
    )
    assert counts == {"p_hom_ref": 2, "p_heterozygous": 0, "p_hom_alt": 1}


def test_simulation_pvalues_use_mixed_maximal_and_shared_banks():
    scores, eligibility = _study_frames()

    checked_scores = focal.validate_scores(scores)
    contrast = checked_scores[checked_scores["statistic"] == "p_hom_alt_minus_hom_ref"]
    assert set(contrast["n_pairs"]) == {-1}
    assert set(contrast["n_valid_pairs"]) == {-1}
    assert set(contrast["n_ref_pairs"]) == {100}
    assert set(contrast["n_alt_pairs"]) == {10}

    pointwise = focal.pointwise_simulation_pvalues(scores, eligibility)
    omnibus = focal.minp_omnibus(scores, eligibility)

    expected_ns = {
        ("EAS", "p_hom_ref"): 4,
        ("EAS", "p_hom_alt"): 4,
        ("EAS", "p_hom_alt_minus_hom_ref"): 4,
        ("Han", "p_hom_ref"): 4,
        ("Han", "p_hom_alt"): 2,
        ("Han", "p_hom_alt_minus_hom_ref"): 2,
    }
    for key, expected_n in expected_ns.items():
        cell = pointwise[
            (pointwise["demography"] == key[0]) & (pointwise["statistic"] == key[1])
        ]
        assert len(cell) == len(focal.THRESHOLDS_YEARS)
        assert set(cell["n_neutral"]) == {expected_n}
    han_ref = pointwise[
        (pointwise["demography"] == "Han") & (pointwise["statistic"] == "p_hom_ref")
    ]
    han_delta = pointwise[
        (pointwise["demography"] == "Han")
        & (pointwise["statistic"] == "p_hom_alt_minus_hom_ref")
    ]
    assert set(han_ref["null_bank"]) == {"maximal_eligible_hom_ref_bank"}
    assert set(han_delta["null_bank"]) == {"shared_complete_case_ref_ref_and_alt_alt"}
    np.testing.assert_allclose(han_delta["mc_p_upper"], 1 / 3, rtol=0, atol=1e-15)
    assert set(omnibus["n_neutral"]) == {2, 4}
    assert set(omnibus["n_thresholds"]) == {7}


def test_reconciliation_requires_complete_finite_exact_cartesian_product():
    _, eligibility = _study_frames()
    reconciliation = _reconciliation(eligibility)

    checked = focal.validate_reconciliation(reconciliation, eligibility)

    assert len(checked) == 10 * len(focal.THRESHOLDS_YEARS)
    assert checked[["unit_id", "threshold_years"]].duplicated().sum() == 0

    one_row = reconciliation.iloc[:1].copy()
    with pytest.raises(ValueError, match="coverage"):
        focal.validate_reconciliation(one_row, eligibility)
    duplicate = pd.concat([reconciliation, reconciliation.iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError, match="duplicate"):
        focal.validate_reconciliation(duplicate, eligibility)
    nonfinite = reconciliation.copy()
    nonfinite.loc[0, "absolute_difference"] = np.nan
    with pytest.raises(ValueError, match="nonfinite"):
        focal.validate_reconciliation(nonfinite, eligibility)
    inconsistent = reconciliation.copy()
    inconsistent.loc[0, "signed_difference"] = 0.01
    with pytest.raises(ValueError, match="signed differences"):
        focal.validate_reconciliation(inconsistent, eligibility)
    wrong_schema = reconciliation.copy()
    wrong_schema.loc[0, "schema"] = "wrong"
    with pytest.raises(ValueError, match="schema"):
        focal.validate_reconciliation(wrong_schema, eligibility)


def test_figure_validator_rejects_fake_pdf_and_wrong_png_dimensions(tmp_path):
    fake_pdf = tmp_path / "fake.pdf"
    fake_pdf.write_bytes(b"%PDF-synthetic-not-a-real-document" + b"x" * 300)
    with pytest.raises(ValueError, match="PDF"):
        focal._validate_figure(fake_pdf, "pdf")

    wrong_png = tmp_path / "wrong.png"
    wrong_png.write_bytes(
        b"\x89PNG\r\n\x1a\n"
        + (13).to_bytes(4, "big")
        + b"IHDR"
        + (100).to_bytes(4, "big")
        + (100).to_bytes(4, "big")
        + b"\x08\x06\x00\x00\x00"
        + b"\x00\x00\x00\x00"
        + b"\x00\x00\x00\x00IEND\x00\x00\x00\x00"
    )
    with pytest.raises(ValueError, match="3300x2550"):
        focal._validate_figure(wrong_png, "png")


def test_plot_has_letter_layout_curves_denominators_and_minp(tmp_path, monkeypatch):
    scores, eligibility = _study_frames()
    pointwise = focal.pointwise_simulation_pvalues(scores, eligibility)
    omnibus = focal.minp_omnibus(scores, eligibility)
    captured: dict[str, object] = {}

    def fake_save(figure, output_dir, stem):
        canvas = FigureCanvasAgg(figure)
        canvas.draw()
        renderer = canvas.get_renderer()
        top_left, top_right, bottom_left, bottom_right = figure.axes
        captured["layout_clear"] = (
            bottom_left.title.get_window_extent(renderer).y1 < top_left.bbox.y0
            and bottom_right.title.get_window_extent(renderer).y1 < top_right.bbox.y0
            and figure.legends[0].get_window_extent(renderer).y0
            > max(top_left.bbox.y1, top_right.bbox.y1)
        )
        captured["annotations_inside"] = all(
            axis.bbox.contains(*corner)
            for axis in (top_left, top_right)
            for text in axis.texts
            for corner in (
                (
                    text.get_bbox_patch().get_window_extent(renderer).x0,
                    text.get_bbox_patch().get_window_extent(renderer).y0,
                ),
                (
                    text.get_bbox_patch().get_window_extent(renderer).x1,
                    text.get_bbox_patch().get_window_extent(renderer).y1,
                ),
            )
        )
        captured["size"] = tuple(figure.get_size_inches())
        captured["axes"] = [
            {
                "title": axis.get_title(),
                "texts": [text.get_text() for text in axis.texts],
                "text_has_box": [
                    text.get_bbox_patch() is not None for text in axis.texts
                ],
                "line_count": len(axis.lines),
            }
            for axis in figure.axes
        ]
        return {
            "png": output_dir / f"{stem}.png",
            "pdf": output_dir / f"{stem}.pdf",
        }

    monkeypatch.setattr(focal, "_atomic_save_figure_pair", fake_save)
    outputs = focal.plot_focal_pair_results(
        scores, eligibility, pointwise, omnibus, tmp_path
    )

    assert captured["size"] == pytest.approx((11, 8.5))
    assert bool(captured["layout_clear"])
    assert bool(captured["annotations_inside"])
    axes = captured["axes"]
    assert isinstance(axes, list)
    assert len(axes) == 4
    combined_titles = " ".join(str(axis["title"]) for axis in axes)
    assert "neutral n: REF/REF=4, ALT/ALT=2, delta=2" in combined_titles
    assert "delta 7-threshold minP" in combined_titles
    lower_titles = [
        str(axis["title"])
        for axis in axes
        if str(axis["title"]).startswith("Pointwise")
    ]
    assert len(lower_titles) == 2
    assert all(title.count("\n") == 2 for title in lower_titles)
    top_annotations = [
        text
        for axis in axes
        for text in axis["texts"]
        if str(text).startswith("neutral n:")
    ]
    assert len(top_annotations) == 2
    assert all("\nselected pairs:" in str(text) for text in top_annotations)
    assert all(all(axis["text_has_box"]) for axis in axes if axis["texts"])
    assert all(axis["line_count"] >= 2 for axis in axes)
    assert outputs["pdf"].suffix == ".pdf"
    assert not outputs["pdf"].exists()


def test_atomic_bundle_custom_base_cache_inventory_and_tamper(tmp_path, monkeypatch):
    repo = tmp_path / "repo"
    repo.mkdir()
    fake_module = repo / "python" / "gamma_smc_aou" / "focal.py"
    fake_module.parent.mkdir(parents=True)
    fake_module.write_text("# synthetic module binding\n", encoding="utf-8")
    monkeypatch.setattr(focal, "MODULE_PATH", fake_module)
    for name, module in (
        ("conditioned_workflow.py", focal.conditioned_workflow),
        ("gamma_analysis.py", focal.gamma_analysis),
        ("gamma_decode.py", focal.gamma_decode),
    ):
        module_path = fake_module.parent / name
        module_path.write_text(f"# synthetic {name} binding\n", encoding="utf-8")
        monkeypatch.setattr(module, "__file__", str(module_path))

    custom_base = repo / "custom" / "conditioned_results"
    custom_base.mkdir(parents=True)
    custom_work = repo / "custom" / "decode_work"
    custom_work.mkdir()
    base_contract_sha = "a" * 64
    base_completion = custom_base / "RESULTS_COMPLETION.json"
    base_completion.write_text(
        json.dumps({"contract_sha256": base_contract_sha}) + "\n",
        encoding="utf-8",
    )
    source = custom_work / "unit" / "source.bin"
    source.parent.mkdir()
    source.write_bytes(b"immutable synthetic input")
    input_manifest = pd.DataFrame(
        [
            {
                "unit_id": "synthetic",
                "input_label": "posterior",
                "path": source.relative_to(repo).as_posix(),
                "sha256": focal.sha256_file(source),
                "size_bytes": source.stat().st_size,
            }
        ]
    )
    scores, eligibility = _study_frames()
    reconciliation = _reconciliation(eligibility)

    def fake_save(figure, output_dir, stem):
        del figure
        output_dir.mkdir(parents=True, exist_ok=True)
        png = output_dir / f"{stem}.png"
        pdf = output_dir / f"{stem}.pdf"
        png.write_bytes(b"\x89PNG\r\n\x1a\nsynthetic")
        pdf.write_bytes(b"%PDF-synthetic-not-a-real-document")
        return {"png": png, "pdf": pdf}

    monkeypatch.setattr(focal, "_atomic_save_figure_pair", fake_save)
    validation_calls: list[tuple[str, str]] = []

    def fake_validate(path, extension):
        assert path.is_file()
        validation_calls.append((path.name, extension))

    monkeypatch.setattr(focal, "_validate_figure", fake_validate)
    results = repo / "results" / "focal"
    arguments = {
        "results_dir": results,
        "work_dir": custom_work,
        "base_results_dir": custom_base,
        "base_results_completion_sha256": focal.sha256_file(base_completion),
        "base_results_contract_sha256": base_contract_sha,
        "scores": scores,
        "eligibility": eligibility,
        "input_manifest": input_manifest,
        "reconciliation": reconciliation,
    }

    first = focal.write_results_bundle(repo, **arguments)
    assert first["cache_hit"] is False
    assert validation_calls
    assert {path.name for path in results.iterdir()} == {
        focal.RESULTS_COMPLETION_FILENAME,
        *focal.RESULT_OUTPUTS.values(),
    }
    assert (results / focal.RESULT_OUTPUTS["figure_pdf"]).read_bytes() == (
        b"%PDF-synthetic-not-a-real-document"
    )
    second = focal.write_results_bundle(repo, **arguments)
    assert second["cache_hit"] is True
    assert second["contract"]["source"]["base_results_path"] == (
        custom_base.relative_to(repo).as_posix()
    )
    assert second["contract"]["source"]["work_path"] == (
        custom_work.relative_to(repo).as_posix()
    )
    assert second["contract"]["source"]["reconciliation_expected_unit_count"] == 10
    assert second["contract"]["source"]["reconciliation_expected_row_count"] == 70
    assert set(second["contract"]["implementation"]) == {
        "focal_pairs",
        "conditioned_workflow",
        "gamma_analysis",
        "gamma_decode",
    }
    assert second["contract"]["scientific_contract"][
        "mutation_rate_per_bp_per_generation"
    ] == pytest.approx(focal.MUTATION_RATE)

    wrong_work = repo / "custom" / "wrong_work"
    wrong_work.mkdir()
    with pytest.raises(ValueError, match="work path differs|work root"):
        focal.verify_results_bundle(
            repo,
            results_dir=results,
            expected_work_dir=wrong_work,
        )

    other_base = repo / "custom" / "wrong_results"
    other_base.mkdir()
    (other_base / "RESULTS_COMPLETION.json").write_text(
        base_completion.read_text(encoding="utf-8"), encoding="utf-8"
    )
    with pytest.raises(ValueError, match="path differs"):
        focal.verify_results_bundle(
            repo,
            results_dir=results,
            expected_base_results_dir=other_base,
        )

    score_path = results / focal.RESULT_OUTPUTS["scores"]
    score_path.write_text(
        score_path.read_text(encoding="utf-8") + "\n", encoding="utf-8"
    )
    with pytest.raises(ValueError, match="checksum"):
        focal.verify_results_bundle(repo, results_dir=results)


def test_verify_cli_honors_requested_work_root(tmp_path, monkeypatch, capsys):
    repo = tmp_path / "repo"
    repo.mkdir()
    work = repo / "custom_work"
    base = repo / "custom_base"
    results = repo / "custom_results"
    captured: dict[str, object] = {}

    monkeypatch.setattr(
        focal.conditioned_workflow,
        "verify_results_bundle",
        lambda *_args, **_kwargs: {"status": "complete"},
    )

    def fake_verify(*_args, **kwargs):
        captured.update(kwargs)
        return {
            "status": "complete",
            "cache_hit": True,
            "contract_sha256": "a" * 64,
        }

    monkeypatch.setattr(focal, "verify_results_bundle", fake_verify)
    assert (
        focal.main(
            [
                "verify",
                "--repo-root",
                str(repo),
                "--work-dir",
                str(work),
                "--base-results-dir",
                str(base),
                "--results-dir",
                str(results),
            ]
        )
        == 0
    )

    assert captured["expected_work_dir"] == work.resolve()
    assert captured["expected_base_results_dir"] == base.resolve()
    assert json.loads(capsys.readouterr().out)["cache_hit"] is True
