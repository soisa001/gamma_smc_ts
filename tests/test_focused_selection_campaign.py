from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
import zstandard
from gamma_smc_aou import focused_selection_campaign as focused
from gamma_smc_aou import focused_selection_decode as decode
from gamma_smc_aou import focused_selection_simulation as simulation


def _plan(tmp_path: Path, **overrides) -> focused.FocusedCampaignPlan:
    values = {
        "repo_root": tmp_path,
        "campaign_dir": tmp_path / focused.DEFAULT_CAMPAIGN_DIR,
    }
    values.update(overrides)
    return focused.FocusedCampaignPlan(**values)


def test_default_grid_has_twelve_selected_cells_and_720_unique_units(tmp_path):
    units = focused.build_execution_units(_plan(tmp_path))

    assert len(units) == 720
    assert units["unit_id"].nunique() == 720
    assert units["seed"].nunique() == 720
    counts = units.groupby(
        ["demography_id", "target_allele_frequency", "simulation_class"]
    ).size()
    assert set(counts.xs("selected", level="simulation_class")) == {20}
    assert set(counts.xs("neutral", level="simulation_class")) == {100}
    assert set(units["target_allele_frequency"]) == {0.1, 0.2, 0.3}
    assert set(units["selection_coefficient"]) == {0.0, 0.005, 0.01}
    selected = units[units["simulation_class"] == "selected"]
    assert set(
        selected.groupby(
            ["demography_id", "target_allele_frequency", "selection_coefficient"]
        ).size()
    ) == {10}
    assert set(units["exact_sample_alt_count"]) == {20, 40, 60}


def test_seed_identity_is_stable_when_replicate_targets_extend(tmp_path):
    small = focused.build_execution_units(
        _plan(tmp_path, selected_replicates=2, neutral_replicates=3)
    ).set_index("unit_id")
    large = focused.build_execution_units(
        _plan(tmp_path, selected_replicates=4, neutral_replicates=5)
    ).set_index("unit_id")

    assert small.index.isin(large.index).all()
    pd.testing.assert_series_equal(
        small["seed"], large.loc[small.index, "seed"], check_names=False
    )


def test_plan_is_idempotent_and_fails_closed_on_corruption(tmp_path):
    plan = _plan(tmp_path)
    table, design = focused.write_plan(plan)
    table_mtime = table.stat().st_mtime_ns
    design_mtime = design.stat().st_mtime_ns
    focused.write_plan(plan)
    assert table.stat().st_mtime_ns == table_mtime
    assert design.stat().st_mtime_ns == design_mtime

    frame = pd.read_csv(table, sep="\t")
    frame.loc[0, "seed"] += 1
    frame.to_csv(table, sep="\t", index=False)
    with pytest.raises(ValueError, match="incompatible|checksum"):
        focused.write_plan(plan)


def test_design_records_thresholds_raw_profile_key_and_one_decode(tmp_path):
    plan = _plan(tmp_path)
    units = focused.build_execution_units(plan)
    design = focused.build_study_design(plan, units)

    assert design["scientific_grid"]["total_units"] == 720
    assert design["scientific_grid"]["selected_units"] == 120
    assert design["scientific_grid"]["neutral_units"] == 600
    assert design["time"]["thresholds_years"] == [
        1_000,
        4_500,
        10_000,
        20_000,
        30_000,
        40_000,
        50_000,
    ]
    assert design["decode_contract"]["gamma_raw_invocations_per_unit"] == 1
    assert design["decode_contract"]["flow_field_cache_bp"] == 1_000
    assert design["decode_contract"]["raw_posterior_retained_by_default"] is True
    assert "unit_id" in design["analysis_contract"]["raw_profile_key"]
    assert design["analysis_contract"]["p_value"].startswith(
        "one-sided upper-tail Monte Carlo"
    )
    simulation_contract = design["simulation_contract"]
    assert simulation_contract["candidate_pool_diploids"] == 500
    assert simulation_contract["selected_candidate_pool_diploids"] == 500
    assert simulation_contract["neutral_candidate_pool_diploids"] == 500
    assert simulation_contract["candidate_pool_frequency_gate_scope"] == (
        "selected_and_neutral"
    )
    assert simulation_contract["candidate_pool_frequency_gate_order"] == (
        "before_exact_sample_panel"
    )
    assert simulation_contract["neutral_engine"] == (
        "msprime DTWF for 200 generations then SmcPrimeApproxCoalescent "
        "under matched demography"
    )
    assert simulation_contract["neutral_recent_dtwf_duration_generations"] == 200
    assert simulation_contract["neutral_coalescent_approximation"] == (
        "SMC-prime after 200 generations"
    )
    assert "fixed 5-Mb marginal tree" in simulation_contract["neutral_ascertainment"]
    assert "pulse migrations" in simulation_contract["neutral_ascertainment"]
    json.dumps(design, allow_nan=False)


def test_focused_decode_uses_validated_one_kb_cache_and_records_every_setting():
    settings = focused._decode_settings(6_300, threads=1)

    assert settings["cache_size"] == 1_000
    assert settings["output_at_stride"] == 10_000
    assert settings["output_at_hets"] is False
    assert settings["pair_selector"] == "explicit_within_individual_pairs"
    assert settings["recent_call"] == "median"
    assert settings["threads"] == 1


def test_spatial_accumulator_preserves_cell_mean_sd_and_position_grid():
    unit = {
        "unit_id": "unit-a",
        "demography_id": "eas_phlash_median",
        "simulation_class": "selected",
        "selection_coefficient": 0.005,
        "target_allele_frequency": 0.1,
    }
    rows = []
    for genotype_class in ("overall", "hom_ref", "heterozygous", "hom_alt"):
        for threshold in focused.TMRCA_THRESHOLDS_YEARS:
            for position, probability in ((0, 0.1), (5_000_000, 0.2)):
                rows.append(
                    {
                        **unit,
                        "source": "gamma_smc",
                        "genotype_class": genotype_class,
                        "threshold_years": threshold,
                        "position_0based": position,
                        "mean_p_tmrca_lt_threshold": probability,
                        "mean_tmrca_generations": 500.0,
                    }
                )
    accumulator = {}
    focused._update_spatial_accumulator(accumulator, pd.DataFrame(rows), unit)
    second = pd.DataFrame(rows)
    second["mean_p_tmrca_lt_threshold"] += 0.2
    focused._update_spatial_accumulator(accumulator, second, unit)

    result = focused._spatial_accumulator_frame(accumulator)
    assert set(result["n_units"]) == {2}
    assert set(result["position_0based"]) == {0, 5_000_000}
    assert set(result["mean_p_tmrca_lt_threshold"].round(6)) == {0.2, 0.3}
    assert (result["sd_p_tmrca_lt_threshold"] > 0).all()


def test_campaign_must_be_in_repo_and_not_onedrive(tmp_path):
    with pytest.raises(ValueError, match="inside repo_root"):
        _plan(tmp_path, campaign_dir=tmp_path.parent / "outside").validate()
    one_drive = tmp_path / "OneDrive" / focused.DEFAULT_CAMPAIGN_DIR
    with pytest.raises(ValueError, match="OneDrive"):
        _plan(tmp_path, campaign_dir=one_drive).validate()
    institution = tmp_path / "OneDrive - Institution" / focused.DEFAULT_CAMPAIGN_DIR
    with pytest.raises(ValueError, match="OneDrive"):
        _plan(tmp_path, campaign_dir=institution).validate()


def test_status_phase_validates_plan_and_reports_not_started(tmp_path):
    assert (
        focused.main(
            [
                "status",
                "--repo-root",
                str(tmp_path),
                "--selected-replicates",
                "1",
                "--neutral-replicates",
                "1",
            ]
        )
        == 0
    )
    assert (tmp_path / focused.DEFAULT_CAMPAIGN_DIR / "study_design.json").is_file()
    status = pd.read_csv(
        tmp_path / focused.DEFAULT_CAMPAIGN_DIR / "results" / "simulation_status.tsv",
        sep="\t",
    )
    assert set(status["status"]) == {"not_started"}


def test_status_validates_direct_completion_outputs_and_reports_corruption(tmp_path):
    plan = _plan(tmp_path, selected_replicates=1, neutral_replicates=1)
    units = focused.build_execution_units(plan).iloc[[0]].reset_index(drop=True)
    unit = units.iloc[0].to_dict()
    directory = plan.campaign_dir / "work" / unit["unit_id"]
    paths = {
        "tree": directory / "simulation.trees",
        "pair_table": directory / "sample_manifest.tsv",
        "overall_pairs": directory / "pairs" / "overall.pairs.tsv",
        "truth_profiles": directory / "truth_profiles.tsv.gz",
        "truth_class_summaries": directory / "truth_class_summaries.tsv",
    }
    for label, path in paths.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(label.encode("utf-8"))
    contract = {
        "schema": simulation.SCHEMA_VERSION,
        "unit": {key: unit[key] for key in focused.EXECUTION_COLUMNS},
    }
    completion = {
        "schema": simulation.SCHEMA_VERSION,
        "status": "complete",
        "contract": contract,
        "contract_sha256": focused._canonical_sha256(contract),
        "outputs": {
            label: {
                "path": path.relative_to(directory).as_posix(),
                "sha256": focused._sha256_file(path),
            }
            for label, path in paths.items()
        },
    }
    completion_path = directory / "simulation_complete.json"
    completion_path.write_text(json.dumps(completion), encoding="utf-8")

    status = focused._simulation_status(plan.campaign_dir, units)
    assert status.iloc[0]["status"] == "complete"
    assert status.iloc[0]["status_error"] == ""

    paths["tree"].write_bytes(b"corrupt")
    status = focused._simulation_status(plan.campaign_dir, units)
    assert status.iloc[0]["status"] == "invalid_completion"
    assert "tree" in status.iloc[0]["status_error"]


def test_status_reports_active_lock_before_stale_failure(tmp_path):
    plan = _plan(tmp_path, selected_replicates=1, neutral_replicates=1)
    units = focused.build_execution_units(plan).iloc[[0]].reset_index(drop=True)
    unit_id = str(units.iloc[0]["unit_id"])
    directory = plan.campaign_dir / "work" / unit_id
    directory.mkdir(parents=True)
    (directory / "simulation_failed.json").write_text("{}", encoding="utf-8")

    with simulation.exclusive_unit_lock(
        directory / ".simulation.lock", unit_id=unit_id, phase="simulation"
    ):
        status = focused._simulation_status(plan.campaign_dir, units)
        assert status.iloc[0]["status"] == "running"

    status = focused._simulation_status(plan.campaign_dir, units)
    assert status.iloc[0]["status"] == "failed_or_exhausted"


def test_decode_task_rejects_a_concurrent_owner(tmp_path, monkeypatch):
    unit = {"unit_id": "unit-a"}
    payload = {"unit": unit, "campaign_dir": str(tmp_path)}
    lock_path = tmp_path / "work" / "unit-a" / "decoded" / ".decode.lock"
    monkeypatch.setattr(
        focused,
        "_decode_task_unlocked",
        lambda ignored: {"unit_id": "unit-a", "status": "complete"},
    )

    with (
        simulation.exclusive_unit_lock(lock_path, unit_id="unit-a", phase="decode"),
        pytest.raises(simulation.UnitLockHeld, match="already held"),
    ):
        focused._decode_task(payload)

    assert focused._decode_task(payload)["status"] == "complete"


def test_filters_keep_both_s_but_one_shared_neutral_bank(tmp_path):
    parser = focused.build_parser()
    args = parser.parse_args(["simulate", "--af", "0.1"])
    units = focused.build_execution_units(_plan(tmp_path))
    selected = focused._filtered_units(units, args)

    assert len(selected[selected["simulation_class"] == "neutral"]) == 200
    assert len(selected[selected["simulation_class"] == "selected"]) == 40
    assert set(selected["selection_coefficient"]) == {0.0, 0.005, 0.01}


def test_aggregate_revalidates_decode_contract_and_builds_spatial_outputs(
    tmp_path, monkeypatch
):
    plan = _plan(tmp_path, selected_replicates=1, neutral_replicates=1)
    units = focused.build_execution_units(plan).iloc[[0]].reset_index(drop=True)
    unit = units.iloc[0].to_dict()
    unit_root = plan.campaign_dir / "work" / str(unit["unit_id"])
    decoded = unit_root / "decoded"
    (unit_root / "pairs").mkdir(parents=True)
    decoded.mkdir(parents=True)

    tree_path = unit_root / "simulation.trees"
    pair_table_path = unit_root / "sample_manifest.tsv"
    overall_pairs_path = unit_root / "pairs" / "overall.pairs.tsv"
    truth_profiles_path = unit_root / "truth_profiles.tsv.gz"
    truth_classes_path = unit_root / "truth_class_summaries.tsv"
    tree_path.write_bytes(b"synthetic tree fixture")
    pair_table = pd.DataFrame(
        {
            "gamma_smc_haplotype_0": [0, 2, 4],
            "gamma_smc_haplotype_1": [1, 3, 5],
            "genotype_class": ["hom_ref", "heterozygous", "hom_alt"],
        }
    )
    pair_table.to_csv(pair_table_path, sep="\t", index=False)
    overall_pairs_path.write_text("0\t1\n2\t3\n4\t5\n", encoding="utf-8")

    positions = (0, focused.FOCAL_POSITION_BP, 9_990_000)
    truth_rows = []
    truth_class_rows = []
    for genotype_class in ("overall", "hom_ref", "heterozygous", "hom_alt"):
        for threshold in focused.TMRCA_THRESHOLDS_YEARS:
            prefix = {
                "unit_id": unit["unit_id"],
                "demography_id": unit["demography_id"],
                "simulation_class": unit["simulation_class"],
                "selection_coefficient": unit["selection_coefficient"],
                "target_allele_frequency": unit["target_allele_frequency"],
                "source": "tree_truth",
                "genotype_class": genotype_class,
                "threshold_years": threshold,
            }
            truth_rows.extend(
                {
                    **prefix,
                    "position_0based": position,
                    "mean_p_tmrca_lt_threshold": 0.25,
                    "mean_tmrca_generations": 1_000.0,
                }
                for position in positions
            )
            truth_class_rows.append(
                {
                    **prefix,
                    "region_mean_p_tmrca_lt_threshold": 0.25,
                    "focal_mean_p_tmrca_lt_threshold": 0.25,
                    "region_mean_tmrca_generations": 1_000.0,
                    "focal_mean_tmrca_generations": 1_000.0,
                    "n_pairs": 3 if genotype_class == "overall" else 1,
                    "focal_output_position_0based": focused.FOCAL_POSITION_BP,
                    "focal_offset_bp": 0,
                }
            )
    pd.DataFrame(truth_rows).to_csv(
        truth_profiles_path, sep="\t", index=False, compression="gzip"
    )
    pd.DataFrame(truth_class_rows).to_csv(truth_classes_path, sep="\t", index=False)

    simulation_contract = {
        "schema": simulation.SCHEMA_VERSION,
        "unit": {key: unit[key] for key in focused.EXECUTION_COLUMNS},
        "implementation": {"slim": {"path": "bin/slim", "sha256": "b" * 64}},
    }
    simulation_outputs = {
        "tree": tree_path,
        "pair_table": pair_table_path,
        "overall_pairs": overall_pairs_path,
        "truth_profiles": truth_profiles_path,
        "truth_class_summaries": truth_classes_path,
    }
    simulation_completion = {
        "schema": simulation.SCHEMA_VERSION,
        "status": "complete",
        "contract": simulation_contract,
        "contract_sha256": focused._canonical_sha256(simulation_contract),
        "outputs": {
            label: {
                "path": path.relative_to(unit_root).as_posix(),
                "sha256": focused._sha256_file(path),
                "size_bytes": path.stat().st_size,
            }
            for label, path in simulation_outputs.items()
        },
    }
    simulation_completion_path = unit_root / "simulation_complete.json"
    simulation_completion_path.write_text(
        json.dumps(simulation_completion), encoding="utf-8"
    )

    decoder_path = tmp_path / "bin" / "gamma_smc"
    decoder_path.parent.mkdir()
    decoder_path.write_bytes(b"synthetic decoder")
    implementation_sources = {"synthetic-source": "a" * 64}
    present_ne = 6_300.0
    monkeypatch.setattr(
        focused, "_decode_source_hashes", lambda ignored: implementation_sources
    )
    monkeypatch.setattr(
        simulation,
        "load_cell_demography",
        lambda ignored_unit, ignored_root: (None, "EAS", present_ne),
    )
    monkeypatch.setattr(
        simulation,
        "_simulate_unit_unlocked",
        lambda ignored_unit, ignored_root, ignored_campaign, **ignored_kwargs: (
            SimpleNamespace(cache_hit=True)
        ),
    )
    provenance = {
        "decoder_sha256": focused._sha256_file(decoder_path),
        "decoder_path": focused._portable_decoder_path(decoder_path, tmp_path),
        "simulation_completion_sha256": focused._sha256_file(
            simulation_completion_path
        ),
        "tree_sha256": focused._sha256_file(tree_path),
        "pair_table_sha256": focused._sha256_file(pair_table_path),
        "overall_pairs_sha256": focused._sha256_file(overall_pairs_path),
        "implementation_sources": implementation_sources,
        "settings": focused._decode_settings(present_ne, threads=1),
    }

    raw_path = decoded / "posterior.zst"
    metadata = {
        "scaled_mutation_rate": provenance["settings"]["scaled_mutation_rate"],
        "sequence_length": len(positions),
        "chunk_size": 4,
        "num_pairs": 3,
        "output_positions": list(positions),
        "pairs": [[0, 1], [2, 3], [4, 5]],
    }
    raw_path.with_name(raw_path.name + ".meta").write_text(
        json.dumps(metadata), encoding="utf-8"
    )
    alpha = np.ones((len(positions), 4), dtype=np.float32)
    beta = np.ones((len(positions), 4), dtype=np.float32)
    raw_path.write_bytes(
        zstandard.ZstdCompressor().compress(alpha.tobytes() + beta.tobytes())
    )
    auxiliary_paths = {
        label: decoded / filename
        for label, filename in focused.DECODE_AUXILIARY_OUTPUTS.items()
    }
    auxiliary_paths["overall_summary"].write_text("position\tvalue\n", encoding="utf-8")
    auxiliary_paths["run_manifest"].write_text("{}\n", encoding="utf-8")
    auxiliary_paths["decoder_stdout"].write_bytes(
        zstandard.ZstdCompressor().compress(b"stdout")
    )
    auxiliary_paths["decoder_stderr"].write_bytes(
        zstandard.ZstdCompressor().compress(b"")
    )
    decode.postprocess_raw_posteriors(
        raw_path,
        pair_table_path,
        decoded,
        unit_record={
            **unit,
            "present_ne": present_ne,
            "decode_provenance": provenance,
        },
        focal_position_0based=focused.FOCAL_POSITION_BP,
        auxiliary_output_paths=auxiliary_paths,
    )

    outputs = focused._aggregate_compact_decode(plan.campaign_dir, units)
    assert len(outputs) == 5
    assert all(path.is_file() for path in outputs)
    gamma_spatial = pd.read_csv(outputs[3], sep="\t")
    truth_spatial = pd.read_csv(outputs[4], sep="\t")
    assert set(gamma_spatial["source"]) == {"gamma_smc"}
    assert set(truth_spatial["source"]) == {"tree_truth"}
    assert set(gamma_spatial["n_units"]) == {1}

    decoder_path.write_bytes(b"changed decoder")
    with pytest.raises(ValueError, match="current decoder binary"):
        focused._aggregate_compact_decode(plan.campaign_dir, units)
