from __future__ import annotations

import hashlib
import json
import os
import socket
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from matplotlib import image as mpl_image

import gamma_smc_aou.eas_replicate_analysis as replicate
from gamma_smc_aou.eas_sweep_study import StudyRuntime


THRESHOLDS = (1_000, 4_500, 10_000, 20_000, 30_000, 40_000, 50_000)


def _runtime(tmp_path: Path) -> StudyRuntime:
    return StudyRuntime(repo_root=str(tmp_path), sample_diploids=10)


def _base(identifier: str, selection: float, frequency: float) -> dict:
    return {
        "scenario": "neanderthal_introgression_standing_variation",
        "study_type": "introgression",
        "specification_id": identifier,
        "selection_coefficient": selection,
        "target_allele_frequency": frequency,
        "lower_allele_frequency": frequency,
        "upper_allele_frequency": 0.95,
        "origin_age_generations": 2_400.0,
        "source_minimum_frequency": 0.9,
    }


def _expanded(base: dict, index: int) -> dict:
    identifier = f"{base['specification_id']}_rep{index:03d}"
    return {
        **base,
        "base_specification_id": base["specification_id"],
        "replicate_index": index,
        "replicate_id": identifier,
        "specification_id": identifier,
    }


def _profile(source: str, *, shift: float = 0.0) -> pd.DataFrame:
    reference = np.array([0.00, 0.08, 0.18, 0.30, 0.42, 0.52, 0.60])
    alternate = np.array([0.03, 0.18, 0.34, 0.52, 0.67, 0.79, 0.88]) + shift
    rows = []
    for genotype_class, probabilities, mean_tmrca in (
        ("overall", (reference + alternate) / 2, 800.0),
        ("hom_ref", reference, 1_100.0),
        ("heterozygous", (reference + alternate) / 2, 800.0),
        ("hom_alt", alternate, 500.0),
    ):
        for threshold, probability in zip(THRESHOLDS, probabilities, strict=True):
            rows.append(
                {
                    "source": source,
                    "position_0based": 5_000_000.0,
                    "position_1based": 5_000_001.0,
                    "genotype_class": genotype_class,
                    "n_pairs": 10,
                    "generation_time_years": 25.0,
                    "threshold_years": float(threshold),
                    "threshold_generations": threshold / 25,
                    "p_tmrca_lt_threshold": float(np.clip(probability, 0, 1)),
                    "mean_p_tmrca_lt_threshold": float(np.clip(probability, 0, 1)),
                    "frac_recent": float(np.clip(probability, 0, 1)),
                    "mean_tmrca_generations": mean_tmrca,
                    "median_tmrca_generations": mean_tmrca,
                    "mean_tmrca_years": mean_tmrca * 25,
                    "median_tmrca_years": mean_tmrca * 25,
                }
            )
    return pd.DataFrame(rows)


def _write_profiles(campaign: Path, specification_id: str) -> None:
    specification_dir = campaign / "work" / specification_id
    specification_dir.mkdir(parents=True)
    _profile("tree_truth").to_csv(
        specification_dir / "truth_profiles.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    _profile("gamma_smc", shift=-0.01).to_csv(
        specification_dir / "decoded_profiles.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )


def _completion(af: float, seed: int) -> dict:
    return {
        "accepted_seed": seed,
        "attempts_to_accept": 2,
        "completed_external_sample_panel_draws": 2,
        "interrupted_internal_conditioned_trajectories": 0,
        "validation": {
            "achieved_allele_frequency": af,
            "achieved_allele_frequency_scope": "sample_estimate",
            "n_hom_ref": 3,
            "n_heterozygous": 4,
            "n_hom_alt": 3,
        },
    }


def test_replicate_specifications_reject_duplicate_indices_and_changed_biology():
    base = _base("intro_s0p010_af10", 0.01, 0.1)
    first = _expanded(base, 0)
    duplicate = {
        **_expanded(base, 0),
        "specification_id": "another_id",
        "replicate_id": "another_id",
    }

    with pytest.raises(ValueError, match="duplicate replicate_index"):
        replicate._normalise_specifications([base], [first, duplicate])  # noqa: SLF001

    changed = {**_expanded(base, 1), "selection_coefficient": 0.005}
    with pytest.raises(ValueError, match="changes biological field"):
        replicate._normalise_specifications([base], [changed])  # noqa: SLF001


def test_wilson_and_bootstrap_are_bounded_and_deterministic():
    low, high = replicate._wilson_interval(7, 10)  # noqa: SLF001
    assert 0 < low < 0.7 < high < 1
    assert all(
        np.isnan(value)
        for value in replicate._wilson_interval(0, 0)  # noqa: SLF001
    )

    values = np.array([1.0, 2.0, 4.0])
    first = replicate._bootstrap_mean_interval(values, seed=123, draws=500)  # noqa: SLF001
    second = replicate._bootstrap_mean_interval(values, seed=123, draws=500)  # noqa: SLF001
    assert first == second
    assert first[0] <= values.mean() <= first[1]


def test_profile_aggregation_preserves_exact_biological_cell_and_handles_nan(
    tmp_path,
):
    base_a = _base("intro_s0p010_af10", 0.01, 0.1)
    base_b = _base("intro_s0p010_af20", 0.01, 0.2)
    frames = []
    for base, value in ((base_a, 0.1), (base_b, 0.9)):
        frame = _profile("tree_truth")
        frame["p_tmrca_lt_threshold"] = value
        frame["base_specification_id"] = base["specification_id"]
        frame["replicate_id"] = f"{base['specification_id']}_rep000"
        frame["replicate_index"] = 0
        frame["specification_id"] = f"{base['specification_id']}_rep000"
        frame["selection_coefficient"] = base["selection_coefficient"]
        frame["target_allele_frequency"] = base["target_allele_frequency"]
        frames.append(frame)

    aggregated = replicate._aggregate_profiles(  # noqa: SLF001
        pd.concat(frames, ignore_index=True), _runtime(tmp_path)
    )

    assert set(aggregated["base_specification_id"]) == {
        base_a["specification_id"],
        base_b["specification_id"],
    }
    assert set(aggregated["target_allele_frequency"]) == {0.1, 0.2}
    overall = aggregated[
        (aggregated["genotype_class"] == "overall")
        & (aggregated["threshold_years"] == THRESHOLDS[0])
    ].set_index("target_allele_frequency")
    assert overall.loc[0.1, "mean_p_tmrca_lt_threshold"] == pytest.approx(0.1)
    assert overall.loc[0.2, "mean_p_tmrca_lt_threshold"] == pytest.approx(0.9)

    nan_frame = frames[0].copy()
    nan_frame["p_tmrca_lt_threshold"] = np.nan
    nan_aggregate = replicate._aggregate_profiles(  # noqa: SLF001
        nan_frame, _runtime(tmp_path)
    )
    assert (nan_aggregate["n_observed_replicates"] == 0).all()
    assert nan_aggregate["mean_p_tmrca_lt_threshold"].isna().all()
    assert nan_aggregate["bootstrap_mean_ci95_low"].isna().all()


def test_paired_accuracy_does_not_count_jointly_wrong_signs_as_detection():
    base = _base("intro_s0p010_af10", 0.01, 0.1)
    paired = pd.DataFrame(
        {
            "base_specification_id": [base["specification_id"]] * 3,
            "selection_coefficient": [0.01] * 3,
            "target_allele_frequency": [0.1] * 3,
            "metric": ["alt_minus_ref_normalized_cdf_auc"] * 3,
            # right/right, wrong/wrong, and right/wrong
            "tree_truth": [1.0, -1.0, 1.0],
            "gamma_smc": [2.0, -2.0, -1.0],
        }
    )

    summary = replicate._paired_accuracy([base], paired)  # noqa: SLF001
    row = summary[
        (summary["summary_scope"] == "overall_micro_trajectory_weighted")
        & (summary["metric"] == "alt_minus_ref_normalized_cdf_auc")
    ].iloc[0]
    assert row["truth_expected_sign_rate"] == pytest.approx(2 / 3)
    assert row["gamma_expected_sign_rate"] == pytest.approx(1 / 3)
    assert row["raw_truth_gamma_sign_agreement_rate"] == pytest.approx(2 / 3)
    assert row["truth_expected_sign_denominator"] == 2
    assert row["gamma_detection_given_truth_expected_sign_count"] == 1
    assert row["gamma_detection_given_truth_expected_sign_rate"] == pytest.approx(0.5)


def test_headline_accuracy_is_cell_equal_not_trajectory_weighted():
    base_a = _base("intro_s0p010_af10", 0.01, 0.1)
    base_b = _base("intro_s0p001_af90", 0.001, 0.9)
    paired = pd.DataFrame(
        {
            "base_specification_id": [base_a["specification_id"]] * 10
            + [base_b["specification_id"]] * 2,
            "metric": ["alt_minus_ref_normalized_cdf_auc"] * 12,
            "tree_truth": [1.0] * 12,
            "gamma_smc": [1.0] * 10 + [-1.0] * 2,
        }
    )

    summary = replicate._paired_accuracy([base_a, base_b], paired)  # noqa: SLF001
    metric = summary[summary["metric"] == "alt_minus_ref_normalized_cdf_auc"].set_index(
        "summary_scope"
    )
    assert metric.loc[
        "overall_micro_trajectory_weighted", "gamma_expected_sign_rate"
    ] == pytest.approx(10 / 12)
    assert metric.loc[
        "overall_macro_cell_equal", "gamma_expected_sign_rate"
    ] == pytest.approx(0.5)
    assert bool(metric.loc["overall_macro_cell_equal", "is_headline"])
    assert (
        metric.loc["overall_macro_cell_equal", "n_biological_cells_meeting_minimum"]
        == 2
    )


def test_decoder_sha_validation_rejects_mixed_and_explicit_mismatch(tmp_path):
    bases = [
        _base("intro_s0p010_af10", 0.01, 0.1),
        _base("intro_s0p005_af20", 0.005, 0.2),
    ]
    records = [_expanded(base, 0) for base in bases]
    for record, decoder_sha in zip(records, ("a" * 64, "b" * 64), strict=True):
        directory = tmp_path / "work" / record["specification_id"]
        (directory / "decoded").mkdir(parents=True)
        (directory / "decode_complete.json").write_text("{}\n", encoding="utf-8")
        (directory / "decoded" / "overall.contract.json").write_text(
            json.dumps({"decoder_sha256": decoder_sha}) + "\n", encoding="utf-8"
        )

    with pytest.raises(RuntimeError, match="mixed Gamma-SMC"):
        replicate._observed_decoder_sha256s(  # noqa: SLF001
            tmp_path, records, None
        )

    second = tmp_path / "work" / records[1]["specification_id"]
    (second / "decoded" / "overall.contract.json").write_text(
        json.dumps({"decoder_sha256": "a" * 64}) + "\n", encoding="utf-8"
    )
    with pytest.raises(RuntimeError, match="explicit Gamma-SMC decoder SHA"):
        replicate._observed_decoder_sha256s(  # noqa: SLF001
            tmp_path, records, "c" * 64
        )
    assert replicate._observed_decoder_sha256s(  # noqa: SLF001
        tmp_path, records, "a" * 64
    ) == ["a" * 64]


def test_quiescence_rejects_stale_task_lock_and_unterminated_phase(tmp_path):
    specification_dir = tmp_path / "work" / "rep001"
    specification_dir.mkdir(parents=True)
    (specification_dir / ".simulation.lock").write_text(
        json.dumps(
            {
                "pid": 2_147_483_000,
                "hostname": socket.gethostname(),
            }
        )
        + "\n",
        encoding="utf-8",
    )
    with pytest.raises(RuntimeError, match="stale_unreconciled"):
        replicate._assert_campaign_quiescent(tmp_path)  # noqa: SLF001
    (specification_dir / ".simulation.lock").unlink()

    event_dir = tmp_path / "phase_invocations"
    event_dir.mkdir()
    (event_dir / "started.json").write_text(
        json.dumps(
            {
                "invocation_id": "unfinished",
                "event": "started",
                "process_id": 2_147_483_000,
            }
        )
        + "\n",
        encoding="utf-8",
    )
    with pytest.raises(RuntimeError, match="unterminated campaign phase"):
        replicate._assert_campaign_quiescent(tmp_path)  # noqa: SLF001


def test_aggregation_lock_recovers_dead_local_owner(tmp_path):
    lock_path = tmp_path / ".replicate_aggregation.lock"
    tmp_path.mkdir(exist_ok=True)
    lock_path.write_text(
        json.dumps(
            {
                "pid": 2_147_483_000,
                "hostname": socket.gethostname(),
                "created_at_utc": "2000-01-01T00:00:00+00:00",
                "token": "old",
            }
        )
        + "\n",
        encoding="utf-8",
    )
    with replicate._campaign_aggregation_lock(tmp_path):  # noqa: SLF001
        current = json.loads(lock_path.read_text(encoding="utf-8"))
        assert current["pid"] == os.getpid()
    assert not lock_path.exists()
    assert list(tmp_path.glob(".replicate_aggregation.lock.stale.*"))


def test_snapshot_copies_compact_coverage_inputs_and_hashes_tree(tmp_path):
    base = _base("intro_s0p010_af10", 0.01, 0.1)
    record = _expanded(base, 0)
    campaign = tmp_path / "campaign"
    specification_dir = campaign / "work" / record["specification_id"]
    specification_dir.mkdir(parents=True)
    expected = {
        "external_draw_state.json": '{"active_internal_conditioned_trajectory": null}\n',
        "external_draws.tsv": "seed\tstatus\n1\trejected\n",
        "simulation_invocations.tsv": "status\telapsed_seconds\ntimed_out\t3\n",
        "simulation_failed.json": '{"status": "timed_out"}\n',
    }
    for filename, contents in expected.items():
        (specification_dir / filename).write_text(contents, encoding="utf-8")
    (specification_dir / "selected.trees").write_bytes(b"tree bytes")
    event_dir = campaign / "phase_invocations"
    event_dir.mkdir(parents=True)
    (event_dir / "phase.json").write_text(
        '{"invocation_id":"x","event":"completed"}\n', encoding="utf-8"
    )

    snapshot = replicate._build_input_snapshot(campaign, [record])  # noqa: SLF001
    records = {row["original_path"]: row for row in snapshot["manifest"]["files"]}
    for filename in expected:
        relative = f"work/{record['specification_id']}/{filename}"
        assert records[relative]["snapshot_mode"] == "copied"
        copied = Path(snapshot["root"]) / relative
        assert copied.read_text(encoding="utf-8") == expected[filename]
    tree_relative = f"work/{record['specification_id']}/selected.trees"
    assert records[tree_relative]["snapshot_mode"] == "hash_only"
    assert records[tree_relative]["sha256"] == hashlib.sha256(b"tree bytes").hexdigest()
    assert records["phase_invocations/phase.json"]["snapshot_mode"] == "copied"


def test_decode_status_evidence_distinguishes_failure_partial_and_stale(tmp_path):
    assert replicate._decode_status_evidence(  # noqa: SLF001
        tmp_path, True, None, {"status": "failed", "error": "boom"}
    ) == ("failed", None, "boom")
    (tmp_path / "decoded").mkdir()
    (tmp_path / "decoded" / "partial.tsv").write_text("x\n", encoding="utf-8")
    assert (
        replicate._decode_status_evidence(  # noqa: SLF001
            tmp_path, True, None, None
        )[0]
        == "incomplete_or_corrupt"
    )
    (tmp_path / "decode_complete.json").write_text("{}\n", encoding="utf-8")
    assert (
        replicate._decode_status_evidence(  # noqa: SLF001
            tmp_path, True, None, None
        )[0]
        == "stale_or_incompatible"
    )


def test_active_draw_is_launched_but_not_interrupted_until_terminal():
    launched, interrupted = replicate._external_draw_counts(  # noqa: SLF001
        {
            "next_external_draw_zero_based": 4,
            "active_internal_conditioned_trajectory": {"seed": 123},
            "interrupted_internal_conditioned_trajectories": [{"seed": 122}],
        },
        [{}, {}, {}, {}],
    )
    assert launched == 5
    assert interrupted == 1


def test_empty_campaign_handles_all_nan_summaries_and_plots(tmp_path, monkeypatch):
    base = _base("intro_s0p001_af90", 0.001, 0.9)
    record = _expanded(base, 0)
    campaign = tmp_path / "introgression_EAS_sim" / "empty_campaign"
    monkeypatch.setattr(
        replicate,
        "_valid_recorded_simulation_cache",
        lambda specification_dir, specification, runtime: None,
    )
    monkeypatch.setattr(
        replicate,
        "_recorded_simulation_contract_if_current",
        lambda specification_dir, specification, runtime: None,
    )

    artifacts = replicate.aggregate_introgression_replicates(
        _runtime(tmp_path), [base], [record], campaign
    )

    effects = pd.read_csv(artifacts["cell_effect_summary_tsv"], sep="\t")
    assert effects["mean"].isna().all()
    assert (effects["n_observed_replicates"] == 0).all()
    accuracy = pd.read_csv(artifacts["truth_gamma_accuracy_tsv"], sep="\t")
    assert (accuracy["n_paired_replicates"] == 0).all()
    assert accuracy["gamma_detection_given_truth_expected_sign_rate"].isna().all()
    assert Path(artifacts["figure_effect_sign_png"]).is_file()
    assert "figure_tmrca_cdf_s0p001_af90_90_pdf" in artifacts


def test_aggregate_introgression_replicates_writes_valid_outputs(tmp_path, monkeypatch):
    base_a = _base("intro_s0p010_af10", 0.01, 0.1)
    base_b = _base("intro_s0p005_af20", 0.005, 0.2)
    records = [
        _expanded(base_a, 0),
        _expanded(base_a, 1),
        _expanded(base_b, 0),
    ]
    campaign = tmp_path / "introgression_EAS_sim" / "replicate_study"
    campaign.mkdir(parents=True)
    (campaign / "study_design.json").write_text(
        json.dumps({"schema": "test-plan"}) + "\n", encoding="utf-8"
    )
    pd.DataFrame(records).to_csv(campaign / "specifications.tsv", sep="\t", index=False)
    for record in records[:2]:
        _write_profiles(campaign, record["specification_id"])
        specification_dir = campaign / "work" / record["specification_id"]
        (specification_dir / "simulation_contract.json").write_text(
            json.dumps({"specification_id": record["specification_id"]}) + "\n",
            encoding="utf-8",
        )
        (specification_dir / "simulation_complete.json").write_text(
            json.dumps({"status": "complete"}) + "\n", encoding="utf-8"
        )
    decoded_dir = campaign / "work" / records[0]["specification_id"] / "decoded"
    decoded_dir.mkdir()
    (
        campaign / "work" / records[0]["specification_id"] / "decode_contract.json"
    ).write_text(json.dumps({"phase": "decode_bundle"}) + "\n", encoding="utf-8")
    (
        campaign / "work" / records[0]["specification_id"] / "decode_complete.json"
    ).write_text(json.dumps({"status": "complete"}) + "\n", encoding="utf-8")
    (decoded_dir / "overall.summary.tsv.contract.json").write_text(
        json.dumps({"decoder_sha256": "a" * 64}) + "\n", encoding="utf-8"
    )

    completions = {
        records[0]["specification_id"]: _completion(0.45, 1001),
        records[1]["specification_id"]: _completion(0.60, 1002),
    }
    decodes = {records[0]["specification_id"]: {"status": "complete"}}

    monkeypatch.setattr(
        replicate,
        "_valid_recorded_simulation_cache",
        lambda specification_dir, specification, runtime: completions.get(
            specification["specification_id"]
        ),
    )
    monkeypatch.setattr(
        replicate,
        "_valid_recorded_decode_cache",
        lambda specification_dir, specification, runtime, decoder_path=None: (
            decodes.get(specification["specification_id"])
        ),
    )
    current_contract = {"schema": "current-test-contract"}
    expected_contract_sha256 = replicate._contract_sha256(current_contract)  # noqa: SLF001
    monkeypatch.setattr(
        replicate,
        "_recorded_simulation_contract_if_current",
        lambda specification_dir, specification, runtime: current_contract,
    )
    invocation_calls = []
    monkeypatch.setattr(
        replicate,
        "_simulation_invocation_history",
        lambda specification_dir, *, expected_contract_sha256, allow_legacy_migration: (
            invocation_calls.append((expected_contract_sha256, allow_legacy_migration))
            or pd.DataFrame()
        ),
    )

    artifacts = replicate.aggregate_introgression_replicates(
        _runtime(tmp_path), [base_a, base_b], records, campaign
    )

    status = pd.read_csv(artifacts["replicate_simulation_status_tsv"], sep="\t")
    assert status.set_index("specification_id")["simulation_status"].to_dict() == {
        records[0]["specification_id"]: "complete",
        records[1]["specification_id"]: "complete",
        records[2]["specification_id"]: "not_requested",
    }
    coverage = pd.read_csv(artifacts["cell_coverage_tsv"], sep="\t")
    cell_a = coverage.set_index("base_specification_id").loc[base_a["specification_id"]]
    assert (
        cell_a["target_replicates"],
        cell_a["accepted_replicates"],
        cell_a["decoded_replicates"],
    ) == (2, 2, 1)

    contrasts = pd.read_csv(artifacts["replicate_center_contrasts_tsv"], sep="\t")
    assert len(contrasts) == 3
    assert contrasts.groupby("source").size().to_dict() == {
        "gamma_smc": 1,
        "tree_truth": 2,
    }
    effects = pd.read_csv(artifacts["cell_effect_summary_tsv"], sep="\t")
    truth_cdf = effects[
        (effects["base_specification_id"] == base_a["specification_id"])
        & (effects["source"] == "tree_truth")
        & (effects["metric"] == "alt_minus_ref_normalized_cdf_auc")
    ].iloc[0]
    assert truth_cdf["n_observed_replicates"] == 2
    assert truth_cdf["expected_sign_rate"] == pytest.approx(1.0)
    assert 0 < truth_cdf["expected_sign_wilson95_low"] < 1

    accuracy = pd.read_csv(artifacts["truth_gamma_accuracy_tsv"], sep="\t")
    overall_cdf = accuracy[
        (accuracy["summary_scope"] == "overall_micro_trajectory_weighted")
        & (accuracy["metric"] == "alt_minus_ref_normalized_cdf_auc")
    ].iloc[0]
    assert overall_cdf["n_paired_replicates"] == 1
    assert overall_cdf["truth_expected_sign_rate"] == pytest.approx(1.0)
    assert overall_cdf["gamma_expected_sign_rate"] == pytest.approx(1.0)
    assert overall_cdf["raw_truth_gamma_sign_agreement_rate"] == pytest.approx(1.0)
    assert overall_cdf[
        "gamma_detection_given_truth_expected_sign_rate"
    ] == pytest.approx(1.0)
    headline_cdf = accuracy[
        (accuracy["summary_scope"] == "overall_macro_cell_equal")
        & (accuracy["metric"] == "alt_minus_ref_normalized_cdf_auc")
    ].iloc[0]
    assert bool(headline_cdf["is_headline"])
    assert headline_cdf["minimum_paired_replicates_per_cell"] == 2
    assert headline_cdf["n_biological_cells_meeting_minimum"] == 0
    assert np.isnan(headline_cdf["gamma_expected_sign_rate"])
    assert invocation_calls == [(expected_contract_sha256, True)] * len(records)

    aggregate_profiles = pd.read_csv(artifacts["aggregate_cdf_profiles_tsv"], sep="\t")
    assert {
        "base_specification_id",
        "target_allele_frequency",
        "bootstrap_mean_ci95_low",
        "bootstrap_mean_ci95_high",
    }.issubset(aggregate_profiles.columns)
    paired = pd.read_csv(artifacts["truth_gamma_paired_replicates_tsv"], sep="\t")
    assert {
        "tree_truth_expected_sign",
        "gamma_smc_expected_sign",
        "raw_truth_gamma_sign_agreement",
        "gamma_detected_given_truth_expected_sign",
    }.issubset(paired.columns)

    summary = json.loads(
        Path(artifacts["replicate_results_summary_json"]).read_text(encoding="utf-8")
    )
    assert summary["requested_trajectories"] == 3
    assert summary["accepted_trajectories"] == 2
    assert summary["decoded_trajectories"] == 1
    assert any("No neutral-null" in line for line in summary["interpretation"])

    manifest = json.loads(
        Path(artifacts["replicate_results_manifest_json"]).read_text(encoding="utf-8")
    )
    assert manifest["replicate_unit"].endswith("slim_trajectory")
    assert manifest["conditioning_missingness"] == "nonrandom"
    assert (
        manifest["inputs"]["study_design_json"]["sha256"]
        == hashlib.sha256((campaign / "study_design.json").read_bytes()).hexdigest()
    )
    assert (
        manifest["inputs"]["specifications_tsv"]["sha256"]
        == hashlib.sha256((campaign / "specifications.tsv").read_bytes()).hexdigest()
    )
    provenance = manifest["validated_contracts_and_completions"]
    assert len(provenance["validated_simulation_completions"]) == 2
    assert len(provenance["validated_decode_completions"]) == 1
    assert provenance["decoder_sha256s"] == ["a" * 64]
    assert manifest["decoder_sha256s"] == ["a" * 64]
    for record in manifest["artifacts"].values():
        path = campaign / record["path"]
        assert path.is_file()
        assert hashlib.sha256(path.read_bytes()).hexdigest() == record["sha256"]
    for key, path in artifacts.items():
        assert Path(path).is_file(), key
    assert Path(artifacts["figure_coverage_png"]).stat().st_size > 10_000
    assert Path(artifacts["figure_coverage_pdf"]).stat().st_size > 1_000
    coverage_png = mpl_image.imread(artifacts["figure_coverage_png"])
    assert coverage_png.shape[:2] == (1_870, 2_420)
    assert "figure_tmrca_cdf_s0p010_af10_10_png" in artifacts
    assert "figure_tmrca_cdf_s0p005_af20_20_pdf" in artifacts
