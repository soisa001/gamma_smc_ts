from __future__ import annotations

import hashlib
import json
import struct
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from gamma_smc_aou import focused_selection_supplemental_plots as plots


def _synthetic_auc() -> pd.DataFrame:
    rows = []
    for demography_index, demography in enumerate(plots.EXPECTED_DEMOGRAPHIES):
        for af in plots.EXPECTED_AFS:
            for coefficient in plots.EXPECTED_SELECTION_COEFFICIENTS:
                for threshold in plots.TMRCA_THRESHOLDS_YEARS:
                    auc = min(
                        0.96,
                        0.42
                        + 0.06 * demography_index
                        + 0.75 * af
                        + 8 * coefficient
                        + 0.03 * np.log10(threshold),
                    )
                    rows.append(
                        {
                            "demography_id": demography,
                            "target_allele_frequency": af,
                            "selection_coefficient": coefficient,
                            "threshold_years": threshold,
                            "statistic_scope": "focal_nearest",
                            "metric": "p_tmrca_lt_threshold",
                            "statistic": plots.PRIMARY_STATISTIC,
                            "n_selected": 10,
                            "n_neutral": 100,
                            "roc_auc": auc,
                            "roc_auc_ci95_low": max(0, auc - 0.08),
                            "roc_auc_ci95_high": min(1, auc + 0.08),
                        }
                    )
    return pd.DataFrame(rows)


def _synthetic_spatial(source: str) -> pd.DataFrame:
    positions = np.arange(4_900_000, 5_100_001, 10_000)
    rows = []
    for demography_index, demography in enumerate(plots.EXPECTED_DEMOGRAPHIES):
        for simulation_class, coefficients, n_units in (
            ("neutral", (0.0,), 100),
            ("selected", plots.EXPECTED_SELECTION_COEFFICIENTS, 10),
        ):
            for coefficient in coefficients:
                for af in plots.EXPECTED_AFS:
                    for genotype_index, genotype_class in enumerate(
                        plots.EXPECTED_GENOTYPE_CLASSES
                    ):
                        for position in positions:
                            offset = (position - plots.FOCAL_POSITION_BP) / 50_000
                            peak = np.exp(-(offset**2))
                            mean = (
                                0.08
                                + 0.01 * demography_index
                                + 0.015 * genotype_index
                                + 0.02 * af
                                + (
                                    (0.06 + 5 * coefficient) * peak
                                    if simulation_class == "selected"
                                    else 0
                                )
                            )
                            rows.append(
                                {
                                    "source": source,
                                    "demography_id": demography,
                                    "simulation_class": simulation_class,
                                    "selection_coefficient": coefficient,
                                    "target_allele_frequency": af,
                                    "genotype_class": genotype_class,
                                    "threshold_years": plots.DEFAULT_THRESHOLD_YEARS,
                                    "position_0based": position,
                                    "n_units": n_units,
                                    "mean_p_tmrca_lt_threshold": mean,
                                    "sem_p_tmrca_lt_threshold": 0.01,
                                }
                            )
    return pd.DataFrame(rows)


def _bundle(spec: plots.SourceSpec) -> plots.SourceBundle:
    return plots.SourceBundle(
        spec=spec,
        auc=_synthetic_auc(),
        spatial=_synthetic_spatial(spec.spatial_source),
        input_records={
            "analysis_completion": {
                "path": spec.analysis_completion,
                "sha256": "1" * 64,
                "size_bytes": 100,
                "contract_sha256": "2" * 64,
            },
            "auc_table": {
                "path": spec.auc_table,
                "sha256": "3" * 64,
                "size_bytes": 200,
            },
            "spatial_table": {
                "path": spec.spatial_table,
                "sha256": "4" * 64,
                "size_bytes": 300,
            },
        },
        upstream_contract_sha256="2" * 64,
    )


def test_auc_and_localized_layouts_use_letter_geometry_and_external_legends(
    tmp_path, monkeypatch
):
    captured = []

    def capture_pair(figure, directory, stem):
        captured.append((figure, stem))
        return [directory / f"{stem}.png", directory / f"{stem}.pdf"]

    monkeypatch.setattr(plots, "_save_pair", capture_pair)
    bundle = _bundle(plots.SOURCE_SPECS[0])

    auc_paths = plots.plot_auc_profiles(bundle, tmp_path)
    local_paths = plots.plot_localized_profiles(
        bundle,
        tmp_path,
        focal_half_window_bp=100_000,
        threshold_years=10_000,
    )

    assert [path.name for path in auc_paths] == [
        "primary_auc_profiles.png",
        "primary_auc_profiles.pdf",
    ]
    assert len(local_paths) == 4
    assert len(captured) == 3

    auc_figure, auc_stem = captured[0]
    assert auc_stem == "primary_auc_profiles"
    assert tuple(auc_figure.get_size_inches()) == pytest.approx((11, 8.5))
    assert len(auc_figure.axes) == 6
    assert [text.get_text() for text in auc_figure.legends[0].texts] == [
        "s=0.005",
        "s=0.01",
        "AUC=0.5",
    ]
    assert auc_figure.legends[0]._ncols == 3
    for axis in auc_figure.axes:
        assert len(axis.lines) == 3
        assert len(axis.collections) == 2
        assert axis.get_ylim() == pytest.approx((0, 1.03))
    assert auc_figure.subplotpars.top == pytest.approx(0.76)

    local_figure, local_stem = captured[1]
    assert local_stem.endswith("_t10000_window100kb")
    assert tuple(local_figure.get_size_inches()) == pytest.approx((11, 8.5))
    assert len(local_figure.axes) == 6
    assert local_figure.legends[0]._ncols == 4
    assert len(local_figure.legends[0].texts) == 8
    assert "+/-100 kb" in local_figure._suptitle.get_text()
    assert local_figure.subplotpars.top == pytest.approx(0.70)
    for index, axis in enumerate(local_figure.axes):
        assert axis.get_xlim() == pytest.approx((-100, 100))
        if index >= 3:
            assert [tick.get_text() for tick in axis.get_xticklabels()] == [
                "-100",
                "-50",
                "0",
                "50",
                "100",
            ]
        assert any(
            len(line.get_xdata()) == 2 and np.allclose(line.get_xdata(), [0, 0])
            for line in axis.lines
        )
        lower, upper = axis.get_ylim()
        assert lower == pytest.approx(0)
        for collection in axis.collections:
            for path in collection.get_paths():
                y = path.vertices[:, 1]
                assert np.nanmin(y) >= lower - 1e-12
                assert np.nanmax(y) <= upper + 1e-12


def test_generate_verify_cache_and_corruption_fail_closed(tmp_path, monkeypatch):
    bundles = {spec.key: _bundle(spec) for spec in plots.SOURCE_SPECS}
    monkeypatch.setattr(
        plots, "_load_source_bundle", lambda _repo, spec: bundles[spec.key]
    )
    output = tmp_path / "supplemental"
    repo_root = Path(__file__).resolve().parents[1]

    generated = plots.generate_supplemental_plots(repo_root, output_dir=output)

    assert generated["status"] == "generated"
    assert generated["output_count"] == 13
    assert generated["png_count"] == 6
    assert generated["pdf_count"] == 6
    completion_path = output / "supplemental_plots_completion.json"
    completion = json.loads(completion_path.read_text(encoding="utf-8"))
    assert completion["status"] == "complete"
    assert completion["contract"]["parameters"]["focal_half_window_bp"] == 100_000
    assert completion["contract"]["parameters"]["localized_threshold_years"] == 10_000
    assert completion["contract"]["parameters"]["randomness"].startswith("none")
    pngs = sorted(output.rglob("*.png"))
    pdfs = sorted(output.rglob("*.pdf"))
    assert len(pngs) == len(pdfs) == 6
    for png in pngs:
        with png.open("rb") as stream:
            stream.seek(16)
            assert struct.unpack(">II", stream.read(8)) == (3_300, 2_550)
    for pdf in pdfs:
        assert plots._pdf_geometry(pdf) == pytest.approx((1, 792, 612))
        record = completion["outputs"][pdf.relative_to(output).as_posix()]
        assert record["page_count"] == 1
        assert record["width_points"] == pytest.approx(792)
        assert record["height_points"] == pytest.approx(612)
    completion_mtime = completion_path.stat().st_mtime_ns

    cached = plots.generate_supplemental_plots(repo_root, output_dir=output)
    verified = plots.generate_supplemental_plots(
        repo_root, output_dir=output, verify_only=True
    )

    assert cached["status"] == "verified_cached"
    assert verified["status"] == "verified_cached"
    assert completion_path.stat().st_mtime_ns == completion_mtime

    original_completion = completion_path.read_bytes()
    tampered_completion = json.loads(original_completion)
    first_output = min(tampered_completion["outputs"])
    tampered_completion["outputs"][first_output]["path"] = "wrong/path"
    completion_path.write_text(
        json.dumps(tampered_completion, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="record"):
        plots.generate_supplemental_plots(
            repo_root, output_dir=output, verify_only=True
        )
    completion_path.write_bytes(original_completion)

    tampered_completion = json.loads(original_completion)
    tampered_completion["outputs"].pop(first_output)
    completion_path.write_text(
        json.dumps(tampered_completion, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="output set"):
        plots.generate_supplemental_plots(
            repo_root, output_dir=output, verify_only=True
        )
    completion_path.write_bytes(original_completion)

    outside_completion = tmp_path / "outside-completion.json"
    outside_completion.write_bytes(original_completion)
    completion_path.unlink()
    completion_path.symlink_to(outside_completion)
    with pytest.raises(ValueError, match="completion path"):
        plots.generate_supplemental_plots(
            repo_root, output_dir=output, verify_only=True
        )
    completion_path.unlink()
    completion_path.write_bytes(original_completion)

    readme = output / "README.md"
    original_readme = readme.read_bytes()
    outside = tmp_path / "outside-readme.md"
    outside.write_bytes(original_readme)
    readme.unlink()
    readme.symlink_to(outside)
    with pytest.raises(ValueError, match="symlink"):
        plots.generate_supplemental_plots(
            repo_root, output_dir=output, verify_only=True
        )
    readme.unlink()
    readme.write_bytes(original_readme)

    unexpected = output / "unexpected"
    unexpected.mkdir()
    with pytest.raises(ValueError, match="directory"):
        plots.generate_supplemental_plots(
            repo_root, output_dir=output, verify_only=True
        )
    unexpected.rmdir()

    readme.write_text(readme.read_text(encoding="utf-8") + "tamper\n", encoding="utf-8")
    with pytest.raises(ValueError, match="size|hash"):
        plots.generate_supplemental_plots(
            repo_root, output_dir=output, verify_only=True
        )


@pytest.mark.parametrize("window", [0, -10_000, 15_000, True, 100_000.0])
def test_invalid_windows_fail_before_input_loading(tmp_path, monkeypatch, window):
    monkeypatch.setattr(
        plots,
        "_load_source_bundle",
        lambda *_args, **_kwargs: pytest.fail("inputs must not load"),
    )
    with pytest.raises((TypeError, ValueError), match="half-window"):
        plots.generate_supplemental_plots(
            Path(__file__).resolve().parents[1],
            output_dir=tmp_path / "output",
            focal_half_window_bp=window,
        )


@pytest.mark.parametrize("threshold", [0, 5_000, True, 10_000.0])
def test_invalid_thresholds_fail_before_input_loading(tmp_path, monkeypatch, threshold):
    monkeypatch.setattr(
        plots,
        "_load_source_bundle",
        lambda *_args, **_kwargs: pytest.fail("inputs must not load"),
    )
    with pytest.raises((TypeError, ValueError), match="threshold"):
        plots.generate_supplemental_plots(
            Path(__file__).resolve().parents[1],
            output_dir=tmp_path / "output",
            threshold_years=threshold,
        )


def test_auc_validation_rejects_direction_or_interval_rewriting():
    auc = _synthetic_auc()
    valid = plots._validate_auc_table(auc)
    assert len(valid) == 84
    assert valid["roc_auc"].min() < 1

    broken = auc.copy()
    broken.loc[broken.index[0], "roc_auc_ci95_low"] = 0.99
    with pytest.raises(ValueError, match="outside"):
        plots._validate_auc_table(broken)


@pytest.mark.parametrize(
    ("column", "value", "message"),
    [
        ("threshold_years", 1_000.5, "threshold"),
        ("n_selected", 10.9, "selected"),
        ("n_neutral", 100.9, "neutral"),
    ],
)
def test_auc_validation_rejects_fractional_integer_contract_fields(
    column, value, message
):
    broken = _synthetic_auc()
    broken[column] = broken[column].astype(float)
    broken.loc[broken.index[0], column] = value

    with pytest.raises(ValueError, match=message):
        plots._validate_auc_table(broken)


@pytest.mark.parametrize(
    ("values", "label"),
    [
        (pd.Series([10_000, 20_000.5]), "spatial positions"),
        (pd.Series([10, 10.9]), "spatial selected counts"),
        (pd.Series([10_000, np.inf]), "spatial thresholds"),
    ],
)
def test_exact_integer_contract_fields_reject_fractional_or_nonfinite(values, label):
    with pytest.raises(ValueError, match=label):
        plots._exact_integer_array(values, label)


def test_source_hash_is_portable_across_crlf_checkouts(tmp_path):
    lf = tmp_path / "lf.py"
    crlf = tmp_path / "crlf.py"
    lf.write_bytes(b"print('one')\nprint('two')\n")
    crlf.write_bytes(b"print('one')\r\nprint('two')\r\n")

    expected = hashlib.sha256(lf.read_bytes()).hexdigest()
    assert plots._canonical_text_sha256(lf) == expected
    assert plots._canonical_text_sha256(crlf) == expected


def test_required_input_rejects_in_repository_symlink(tmp_path):
    root = tmp_path / "repo"
    root.mkdir()
    actual = root / "actual.tsv"
    actual.write_text("value\n", encoding="utf-8")
    linked = root / "linked.tsv"
    linked.symlink_to(actual)

    with pytest.raises(ValueError, match="symlink"):
        plots._require_regular_file(linked, root)


def test_pdf_geometry_rejects_non_pdf_and_wrong_media_box(tmp_path):
    non_pdf = tmp_path / "not.pdf"
    non_pdf.write_bytes(b"not a PDF")
    with pytest.raises(ValueError, match="invalid"):
        plots._pdf_geometry(non_pdf)

    wrong_size = tmp_path / "wrong.pdf"
    wrong_size.write_bytes(b"%PDF-1.4\n<< /Type /Page /MediaBox [ 0 0 700 612 ] >>\n")
    with pytest.raises(ValueError, match="landscape letter"):
        plots._pdf_geometry(wrong_size)


def test_failed_generation_removes_temporary_tree(tmp_path, monkeypatch):
    bundles = {spec.key: _bundle(spec) for spec in plots.SOURCE_SPECS}
    monkeypatch.setattr(
        plots, "_load_source_bundle", lambda _repo, spec: bundles[spec.key]
    )
    monkeypatch.setattr(
        plots,
        "plot_auc_profiles",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("injected")),
    )
    output = tmp_path / "supplemental"

    with pytest.raises(RuntimeError, match="injected"):
        plots.generate_supplemental_plots(
            Path(__file__).resolve().parents[1], output_dir=output
        )

    assert not output.exists()
    assert not output.with_name(output.name + f".tmp.{plots.os.getpid()}").exists()


def test_cli_exposes_only_generate_and_verify():
    parser = plots._parser()
    action = next(item for item in parser._actions if item.dest == "action")
    assert tuple(action.choices) == ("generate", "verify")
    help_text = parser.format_help().lower()
    assert "simulate" not in help_text
    assert "decode" not in help_text
