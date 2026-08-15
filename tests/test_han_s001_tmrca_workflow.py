import json
from pathlib import Path

import pandas as pd
import pytest

from gamma_smc_aou import han_s001_tmrca_analysis as analysis
from gamma_smc_aou import han_s001_tmrca_workflow as workflow


def _synthetic_scores() -> pd.DataFrame:
    rows = []
    metadata = (
        ("neutral_0000", "neutral", 0.0, 100, 0.10, 10),
        ("neutral_0001", "neutral", 0.0, 101, 0.20, 20),
        ("selected_0000", "selected", 0.001, 200, 0.12, 10),
        ("selected_0001", "selected", 0.001, 201, 0.22, 20),
    )
    for (
        unit_id,
        simulation_class,
        coefficient,
        seed,
        final_af,
        sample_count,
    ) in metadata:
        for source_index, source in enumerate(analysis.ALLOWED_SOURCES):
            for window_index, window in enumerate(analysis.ALLOWED_WINDOWS):
                for threshold_index, threshold in enumerate(analysis.THRESHOLDS_YEARS):
                    rows.append(
                        {
                            "unit_id": unit_id,
                            "simulation_class": simulation_class,
                            "selection_coefficient": coefficient,
                            "seed": seed,
                            "source": source,
                            "window": window,
                            "threshold_years": threshold,
                            "score": (
                                0.04
                                + 0.07 * threshold_index
                                + 0.10 * (simulation_class == "selected")
                                + 0.01 * source_index
                                + 0.005 * window_index
                            ),
                            "final_population_af": final_af,
                            "sample_alt_count": sample_count,
                            "sample_af": sample_count / 200,
                            "sample_detected": True,
                        }
                    )
    return pd.DataFrame(rows, columns=analysis.SCORE_COLUMNS)


def _fake_plot_pair(*args, stem, **_kwargs):
    output_dir = args[-1]
    destination = Path(output_dir).resolve()
    png = destination / f"{stem}.png"
    pdf = destination / f"{stem}.pdf"
    png.write_bytes(b"\x89PNG\r\n\x1a\nsynthetic non-rendered IEND placeholder")
    pdf.write_bytes(b"%PDF-1.4\nsynthetic non-rendered placeholder\n%%EOF\n")
    return png, pdf


def _install_fake_study(monkeypatch, repo_root: Path) -> pd.DataFrame:
    scores = _synthetic_scores()
    plan_path = repo_root / "study" / "study_plan.json"
    plan_path.parent.mkdir(parents=True)
    plan_payload = {
        "schema": "synthetic-plan/v1",
        "n_neutral": 2,
        "n_selected": 2,
    }
    plan_path.write_text(json.dumps(plan_payload, sort_keys=True) + "\n")
    score_path = repo_root / "study" / "replicate_scores.tsv"
    scores.to_csv(score_path, sep="\t", index=False)
    scores = analysis.load_replicate_scores(score_path)

    def write_plan(_root, *, n_neutral, n_selected, force=False):
        assert Path(_root).resolve() == repo_root.resolve()
        assert (n_neutral, n_selected, force) == (2, 2, False)
        return plan_path

    monkeypatch.setattr(workflow.simulation, "write_plan", write_plan)
    monkeypatch.setattr(workflow.simulation, "plan_path", lambda _root: plan_path)
    monkeypatch.setattr(
        workflow.simulation,
        "SCORES_RELATIVE_PATH",
        Path("study/replicate_scores.tsv"),
    )
    monkeypatch.setattr(workflow.simulation, "load_plan", lambda _root: plan_payload)
    monkeypatch.setattr(
        workflow.simulation,
        "verify_study",
        lambda _root: {
            "status": "complete",
            "simulation_completion_sha256": "a" * 64,
            "decode_completion_sha256": "b" * 64,
        },
    )
    monkeypatch.setattr(
        workflow.simulation,
        "collect_replicate_scores",
        lambda _root: scores.copy(),
    )
    monkeypatch.setattr(analysis, "plot_selected_neutral_curves", _fake_plot_pair)
    monkeypatch.setattr(analysis, "plot_auc_profiles", _fake_plot_pair)
    return scores


def _analyze(monkeypatch, tmp_path: Path, *, output_name: str = "analysis", **kwargs):
    _install_fake_study(monkeypatch, tmp_path)
    return workflow.analyze_study(
        tmp_path,
        output_dir=tmp_path / output_name,
        n_neutral=2,
        n_selected=2,
        bootstrap_draws=5,
        permutation_draws=7,
        base_seed=11,
        **kwargs,
    )


def test_atomic_analysis_bundle_is_exact_and_cache_validates(monkeypatch, tmp_path):
    completion = _analyze(monkeypatch, tmp_path)
    destination = tmp_path / "analysis"
    assert completion["status"] == "complete"
    assert completion["cache_hit"] is False
    expected = {
        workflow.COMPLETION_FILENAME,
        *workflow._output_specs().values(),
    }
    assert {path.name for path in destination.iterdir()} == expected
    manifest = json.loads(
        (destination / workflow.COMPLETION_FILENAME).read_text(encoding="utf-8")
    )
    assert set(manifest["outputs"]) == set(workflow._output_specs())
    assert (
        manifest["contract"]["inputs"]["empirical"]["claim_status"]
        == "no_real_empirical_p_value_claimed"
    )
    assert set(manifest["contract"]["implementation"]["software_versions"]) == {
        "python",
        "numpy",
        "pandas",
        "scipy",
        "matplotlib",
    }
    empirical = pd.read_csv(
        destination / analysis.TABLE_FILENAMES["empirical_pvalues"], sep="\t"
    )
    assert empirical["status"].str.startswith("unavailable").all()
    assert not list(tmp_path.glob(".han_s001_tmrca_analysis.staging.*"))

    cached = workflow.analyze_study(
        tmp_path,
        output_dir=destination,
        n_neutral=2,
        n_selected=2,
        bootstrap_draws=5,
        permutation_draws=7,
        base_seed=11,
    )
    assert cached["cache_hit"] is True
    assert cached["contract_sha256"] == completion["contract_sha256"]


@pytest.mark.parametrize("failure_mode", ("tampered", "extra"))
def test_cached_verify_rejects_tampered_or_extra_outputs(
    monkeypatch, tmp_path, failure_mode
):
    _analyze(monkeypatch, tmp_path)
    monkeypatch.setattr(
        workflow.simulation,
        "collect_replicate_scores",
        lambda _root: pytest.fail("verify must not rewrite collected scores"),
    )
    destination = tmp_path / "analysis"
    if failure_mode == "tampered":
        target = destination / analysis.TABLE_FILENAMES["auc_by_threshold"]
        target.write_bytes(target.read_bytes() + b"\n")
        message = "checksum failed"
    else:
        (destination / "undeclared.txt").write_text("not in manifest\n")
        message = "inventory mismatch"
    with pytest.raises(ValueError, match=message):
        workflow.verify_analysis(
            tmp_path,
            output_dir=destination,
            n_neutral=2,
            n_selected=2,
            bootstrap_draws=5,
            permutation_draws=7,
            base_seed=11,
        )


def test_nonempty_partial_analysis_is_not_reused(monkeypatch, tmp_path):
    _install_fake_study(monkeypatch, tmp_path)
    destination = tmp_path / "partial"
    destination.mkdir()
    (destination / "orphan.tsv").write_text("partial\n")
    with pytest.raises(ValueError, match="nonempty without a valid completion"):
        workflow.analyze_study(
            tmp_path,
            output_dir=destination,
            n_neutral=2,
            n_selected=2,
            bootstrap_draws=3,
            permutation_draws=3,
        )


def test_analysis_fails_closed_on_incomplete_verification_or_wrong_counts(
    monkeypatch, tmp_path
):
    scores = _install_fake_study(monkeypatch, tmp_path)
    monkeypatch.setattr(
        workflow.simulation,
        "verify_study",
        lambda _root: {"status": "partial"},
    )
    with pytest.raises(ValueError, match="verification is not complete"):
        workflow.analyze_study(
            tmp_path,
            output_dir=tmp_path / "incomplete",
            n_neutral=2,
            n_selected=2,
            bootstrap_draws=3,
            permutation_draws=3,
        )

    monkeypatch.setattr(
        workflow.simulation,
        "verify_study",
        lambda _root: {"status": "complete"},
    )
    incomplete_scores = scores[scores["unit_id"] != "neutral_0001"].copy()
    monkeypatch.setattr(
        workflow.simulation,
        "collect_replicate_scores",
        lambda _root: incomplete_scores,
    )
    with pytest.raises(ValueError, match="replicate counts do not match"):
        workflow.analyze_study(
            tmp_path,
            output_dir=tmp_path / "wrong_counts",
            n_neutral=2,
            n_selected=2,
            bootstrap_draws=3,
            permutation_draws=3,
        )


def test_fixed_all_table_contract_rejects_disabled_af_stratification(tmp_path):
    with pytest.raises(ValueError, match="fixed all-table contract"):
        workflow.analyze_study(tmp_path, include_af_stratified=False)


def test_explicit_empirical_input_is_provenance_bound(monkeypatch, tmp_path):
    _install_fake_study(monkeypatch, tmp_path)
    empirical = analysis.empirical_input_template()
    empirical["empirical_id"] = "han_observed_1"
    empirical["variant_id"] = "chr1:5000001:A:G"
    empirical["score"] = empirical["threshold_years"] / 100_000
    empirical["sample_alt_count"] = 20
    empirical["sample_af"] = 0.1
    empirical["sample_detected"] = True
    for column in (
        "source_data_sha256",
        "sample_manifest_sha256",
        "decoder_binary_sha256",
        "decoder_parameters_sha256",
        "ne_history_sha256",
    ):
        empirical[column] = "c" * 64
    empirical["qc_status"] = "pass"
    empirical["data_status"] = "empirical_observed"
    empirical_path = tmp_path / "explicit_empirical.tsv"
    empirical.to_csv(empirical_path, sep="\t", index=False)

    completion = workflow.analyze_study(
        tmp_path,
        output_dir=tmp_path / "explicit_analysis",
        empirical_path=empirical_path,
        n_neutral=2,
        n_selected=2,
        bootstrap_draws=3,
        permutation_draws=3,
    )
    record = completion["contract"]["inputs"]["empirical"]
    assert record["mode"] == "explicit_observation_with_self_declared_provenance_fields"
    assert record["sha256"] == workflow._sha256_file(empirical_path)
    result = pd.read_csv(
        tmp_path / "explicit_analysis" / analysis.TABLE_FILENAMES["empirical_pvalues"],
        sep="\t",
    )
    assert set(result["status"]) == {"available_explicit_empirical_input"}


def test_cli_exposes_all_phases_and_all_dispatches_in_order(monkeypatch, tmp_path):
    phases = ("plan", "simulate", "status", "decode", "analyze", "verify", "all")
    for phase in phases:
        assert workflow.build_parser().parse_args([phase]).phase == phase

    calls = []
    plan_path = tmp_path / "plan.json"
    plan_path.write_text("{}\n")
    slim_bin = tmp_path / "slim"
    decoder_bin = tmp_path / "gamma_smc"
    slim_bin.write_bytes(b"synthetic")
    decoder_bin.write_bytes(b"synthetic")

    def write_plan(_root, *, n_neutral, n_selected, force):
        calls.append(("plan", n_neutral, n_selected, force))
        return plan_path

    monkeypatch.setattr(workflow.simulation, "write_plan", write_plan)
    monkeypatch.setattr(workflow.simulation, "plan_path", lambda _root: plan_path)
    monkeypatch.setattr(
        workflow.simulation,
        "load_plan",
        lambda _root: {"n_neutral": 100, "n_selected": 100},
    )
    monkeypatch.setattr(
        workflow.simulation,
        "simulate_study",
        lambda *_args, **_kwargs: (
            calls.append(("simulate",)) or pd.DataFrame({"status": ["complete"]})
        ),
    )
    monkeypatch.setattr(
        workflow.simulation,
        "decode_study",
        lambda *_args, **_kwargs: (
            calls.append(("decode",)) or pd.DataFrame({"status": ["complete"]})
        ),
    )
    monkeypatch.setattr(
        workflow,
        "analyze_study",
        lambda *_args, **_kwargs: (
            calls.append(("analyze",)) or {"status": "complete", "cache_hit": False}
        ),
    )
    monkeypatch.setattr(
        workflow.simulation,
        "verify_study",
        lambda *_args, **_kwargs: (
            calls.append(("simulation_verify",)) or {"status": "complete"}
        ),
    )
    monkeypatch.setattr(
        workflow,
        "verify_analysis",
        lambda *_args, **_kwargs: (
            calls.append(("analysis_verify",))
            or {"status": "complete", "contract_sha256": "d" * 64}
        ),
    )
    assert (
        workflow.main(
            [
                "all",
                "--repo-root",
                str(tmp_path),
                "--slim-bin",
                str(slim_bin),
                "--decoder-bin",
                str(decoder_bin),
            ]
        )
        == 0
    )
    assert [call[0] for call in calls] == [
        "plan",
        "simulate",
        "decode",
        "analyze",
        "simulation_verify",
        "analysis_verify",
    ]
    assert calls[0][1:3] == (100, 100)


def test_cli_status_is_read_only_and_all_stops_after_failed_simulation(
    monkeypatch, tmp_path
):
    plan_path = tmp_path / "plan.json"
    plan_path.write_text("{}\n")
    slim_bin = tmp_path / "slim"
    decoder_bin = tmp_path / "gamma_smc"
    slim_bin.write_bytes(b"synthetic")
    decoder_bin.write_bytes(b"synthetic")
    monkeypatch.setattr(workflow.simulation, "plan_path", lambda _root: plan_path)
    monkeypatch.setattr(
        workflow.simulation,
        "load_plan",
        lambda _root: {"n_neutral": 2, "n_selected": 2},
    )
    monkeypatch.setattr(
        workflow.simulation,
        "write_plan",
        lambda *_args, **_kwargs: pytest.fail("status must not write a plan"),
    )
    monkeypatch.setattr(
        workflow.simulation,
        "simulation_status",
        lambda _root: pd.DataFrame({"status": ["complete"]}),
    )
    monkeypatch.setattr(
        workflow.simulation,
        "decode_status",
        lambda _root: pd.DataFrame({"status": ["complete"]}),
    )
    assert workflow.main(["status", "--repo-root", str(tmp_path)]) == 0

    monkeypatch.setattr(
        workflow.simulation,
        "write_plan",
        lambda *_args, **_kwargs: plan_path,
    )
    monkeypatch.setattr(
        workflow.simulation,
        "simulate_study",
        lambda *_args, **_kwargs: pd.DataFrame({"status": ["failed"]}),
    )
    monkeypatch.setattr(
        workflow.simulation,
        "decode_study",
        lambda *_args, **_kwargs: pytest.fail(
            "decode must not run after failed simulation"
        ),
    )
    with pytest.raises(ValueError, match="simulation did not complete"):
        workflow.main(
            [
                "all",
                "--repo-root",
                str(tmp_path),
                "--slim-bin",
                str(slim_bin),
                "--decoder-bin",
                str(decoder_bin),
            ]
        )


def test_cli_plan_uses_real_simulation_plan_api(tmp_path):
    source_repo = Path(__file__).resolve().parents[1]
    for relative in workflow.simulation.SOURCE_PATHS:
        source = source_repo / relative
        destination = tmp_path / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(source.read_bytes())

    assert (
        workflow.main(
            [
                "plan",
                "--repo-root",
                str(tmp_path),
                "--neutral-replicates",
                "2",
                "--selected-replicates",
                "3",
            ]
        )
        == 0
    )
    plan = workflow.simulation.load_plan(tmp_path)
    assert plan["n_neutral"] == 2
    assert plan["n_selected"] == 3
    assert Path(workflow.simulation.plan_path(tmp_path)).is_file()
