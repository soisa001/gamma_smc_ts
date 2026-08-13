from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from gamma_smc_aou import focused_selection_analysis as analysis
from gamma_smc_aou import focused_selection_report as report
from gamma_smc_aou.focused_selection_campaign import (
    FocusedCampaignPlan,
    write_plan,
)
from gamma_smc_aou.focused_selection_campaign import (
    main as campaign_main,
)
from gamma_smc_aou.focused_selection_decode import (
    SCHEMA_VERSION as DECODE_SCHEMA_VERSION,
)
from gamma_smc_aou.focused_selection_report import generate_run_report
from gamma_smc_aou.focused_selection_simulation import (
    SCHEMA_VERSION as SIMULATION_SCHEMA_VERSION,
)

THRESHOLDS = (1_000.0, 4_500.0, 10_000.0, 20_000.0, 30_000.0, 40_000.0, 50_000.0)
CLASSES = ("overall", "hom_ref", "heterozygous", "hom_alt")


def _canonical_sha256(value: object) -> str:
    payload = json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _campaign(
    tmp_path: Path, *, neutral_replicates: int = 1
) -> tuple[Path, pd.DataFrame]:
    root = tmp_path / "repo"
    root.mkdir()
    campaign = root / "focused_selection_EAS_sim"
    plan = FocusedCampaignPlan(
        repo_root=root,
        campaign_dir=campaign,
        selected_replicates=1,
        neutral_replicates=neutral_replicates,
        base_seed=17,
    )
    table, _ = write_plan(plan)
    return campaign, pd.read_csv(table, sep="\t")


def _json_record(row: pd.Series) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in row.to_dict().items():
        if isinstance(value, np.generic):
            value = value.item()
        result[str(key)] = value
    return result


def _write_simulation_completion(campaign: Path, row: pd.Series) -> Path:
    unit = _json_record(row)
    target = int(unit["exact_sample_alt_count"])
    hom_alt = 2
    heterozygous = target - 2 * hom_alt
    hom_ref = int(unit["sample_diploids"]) - hom_alt - heterozygous
    unit_dir = campaign / "work" / str(unit["unit_id"])
    outputs: dict[str, dict[str, object]] = {}
    for label, relative in report.SIMULATION_OUTPUTS.items():
        output_path = unit_dir / relative
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_bytes(f"synthetic {label}\n".encode())
        outputs[label] = {
            "path": Path(relative).as_posix(),
            "sha256": hashlib.sha256(output_path.read_bytes()).hexdigest(),
        }
    contract = {"schema": SIMULATION_SCHEMA_VERSION, "unit": unit}
    completion = {
        "schema": SIMULATION_SCHEMA_VERSION,
        "status": "complete",
        "contract": contract,
        "contract_sha256": _canonical_sha256(contract),
        "attempts_completed": 3,
        "elapsed_seconds": 12.5,
        "genotype_counts": {
            "hom_ref": hom_ref,
            "heterozygous": heterozygous,
            "hom_alt": hom_alt,
        },
        "sample_alt_count": target,
        "sample_af": target / (2 * int(unit["sample_diploids"])),
        "outputs": outputs,
    }
    path = unit_dir / "simulation_complete.json"
    (unit_dir / "simulation_contract.json").write_text(
        json.dumps(contract), encoding="utf-8"
    )
    path.write_text(json.dumps(completion), encoding="utf-8")
    return path


def _write_decode_completion(campaign: Path, row: pd.Series) -> Path:
    unit = _json_record(row)
    decoded = campaign / "work" / str(unit["unit_id"]) / "decoded"
    outputs: dict[str, dict[str, object]] = {}
    for label, filename in report.DECODE_OUTPUTS.items():
        output_path = decoded / filename
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_bytes(f"synthetic {label}\n".encode())
        outputs[label] = {
            "path": filename,
            "sha256": hashlib.sha256(output_path.read_bytes()).hexdigest(),
            "size_bytes": output_path.stat().st_size,
        }
    contract = {"static": {"unit_record": unit}}
    completion = {
        "schema": DECODE_SCHEMA_VERSION,
        "status": "complete",
        "contract": contract,
        "contract_sha256": _canonical_sha256(contract),
        "outputs": outputs,
    }
    path = decoded / "completion.json"
    path.write_text(json.dumps(completion), encoding="utf-8")
    return path


def _class_summaries(units: pd.DataFrame, source: str) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    metadata = [
        "unit_id",
        "demography_id",
        "simulation_class",
        "selection_coefficient",
        "target_allele_frequency",
        "replicate_index",
    ]
    for unit in units.to_dict(orient="records"):
        selected = str(unit["simulation_class"]) == "selected"
        coefficient = float(unit["selection_coefficient"])
        for threshold_index, threshold in enumerate(THRESHOLDS):
            baseline = 0.03 + 0.07 * threshold_index
            class_values = {
                "overall": baseline + (0.02 if selected else 0),
                "hom_ref": baseline,
                "heterozygous": baseline + (0.03 if selected else 0.005),
                "hom_alt": baseline + ((0.12 + 5 * coefficient) if selected else 0.01),
            }
            counts = {
                "overall": 100,
                "hom_ref": 82,
                "heterozygous": 16,
                "hom_alt": 2,
            }
            for genotype_class in CLASSES:
                value = class_values[genotype_class]
                rows.append(
                    {
                        **{key: unit[key] for key in metadata},
                        "source": source,
                        "genotype_class": genotype_class,
                        "threshold_years": threshold,
                        "region_mean_p_tmrca_lt_threshold": value,
                        "focal_mean_p_tmrca_lt_threshold": min(1.0, value + 0.01),
                        "region_mean_tmrca_generations": 900
                        - 20 * CLASSES.index(genotype_class),
                        "focal_mean_tmrca_generations": 850
                        - 20 * CLASSES.index(genotype_class),
                        "n_pairs": counts[genotype_class],
                        "focal_output_position_0based": 5_000_000,
                        "focal_offset_bp": 0,
                    }
                )
    return pd.DataFrame(rows)


def _pair_summaries(units: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    selected = units[units["simulation_class"] == "selected"]
    metadata = [
        "unit_id",
        "demography_id",
        "simulation_class",
        "selection_coefficient",
        "target_allele_frequency",
        "replicate_index",
    ]
    genotype_classes = (
        "hom_ref",
        "hom_ref",
        "heterozygous",
        "heterozygous",
        "hom_alt",
        "hom_alt",
    )
    for unit in selected.to_dict(orient="records"):
        for pair_index, genotype_class in enumerate(genotype_classes):
            offset = {"hom_ref": 0.0, "heterozygous": 0.04, "hom_alt": 0.15}[
                genotype_class
            ]
            for threshold_index, threshold in enumerate(THRESHOLDS):
                probability = 0.02 + 0.07 * threshold_index + offset
                rows.append(
                    {
                        **{key: unit[key] for key in metadata},
                        "pair_index": pair_index,
                        "genotype_class": genotype_class,
                        "threshold_years": threshold,
                        "region_mean_p_tmrca_lt_threshold": probability,
                        "focal_p_tmrca_lt_threshold": min(1.0, probability + 0.01),
                        "region_mean_tmrca_generations": 800 - 100 * offset,
                        "focal_tmrca_generations": 780 - 100 * offset,
                    }
                )
    return pd.DataFrame(rows)


def _finalize_analyses(campaign: Path, units: pd.DataFrame) -> None:
    results = campaign / "results"
    results.mkdir(exist_ok=True)
    gamma_path = results / "gamma.tsv.gz"
    truth_path = results / "truth.tsv.gz"
    pair_path = results / "pairs.tsv.gz"
    gamma = _class_summaries(units, "gamma_smc")
    truth = _class_summaries(units, "tree_truth")
    gamma.to_csv(gamma_path, sep="\t", index=False, compression="gzip")
    truth.to_csv(truth_path, sep="\t", index=False, compression="gzip")
    _pair_summaries(units).to_csv(pair_path, sep="\t", index=False, compression="gzip")
    analysis.run_focused_analysis(
        gamma_path,
        results / "analysis",
        pair_summaries_path=pair_path,
        minimum_selected_per_cell=1,
        minimum_neutral_per_cell=2,
        bootstrap_draws=9,
        permutation_draws=9,
        crossfit_folds=2,
        base_seed=17,
        make_plots=False,
    )
    analysis.run_focused_analysis(
        truth_path,
        results / "analysis_tree_truth",
        minimum_selected_per_cell=1,
        minimum_neutral_per_cell=2,
        bootstrap_draws=9,
        permutation_draws=9,
        crossfit_folds=2,
        base_seed=17,
        make_plots=False,
    )
    analysis.run_truth_gamma_comparison(
        truth_path,
        gamma_path,
        results / "truth_gamma_comparison",
        minimum_selected_per_cell=1,
        minimum_neutral_per_cell=2,
    )


def test_partial_report_is_deterministic_and_includes_validated_qc(tmp_path: Path):
    campaign, units = _campaign(tmp_path)
    _write_simulation_completion(campaign, units.iloc[0])

    path = generate_run_report(campaign)
    first = path.read_bytes()
    assert "Report status: **partial**" in first.decode()
    assert "1/18 simulations" in first.decode()
    assert "Exact AF and genotype QC" in first.decode()
    assert "Gamma-SMC analysis is not finalized yet" in first.decode()

    assert generate_run_report(campaign).read_bytes() == first


def test_campaign_report_phase_writes_partial_report(tmp_path: Path):
    assert (
        campaign_main(
            [
                "report",
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
    report = tmp_path / "focused_selection_EAS_sim" / "RUN_RESULTS.md"
    assert report.is_file()
    assert "Report status: **partial**" in report.read_text(encoding="utf-8")


def test_report_rejects_completed_unit_with_inexact_af(tmp_path: Path):
    campaign, units = _campaign(tmp_path)
    completion_path = _write_simulation_completion(campaign, units.iloc[0])
    completion = json.loads(completion_path.read_text(encoding="utf-8"))
    completion["sample_alt_count"] += 1
    completion_path.write_text(json.dumps(completion), encoding="utf-8")

    with pytest.raises(ValueError, match="exact target"):
        generate_run_report(campaign)


def test_report_rejects_stale_simulation_schema_and_corrupt_output(tmp_path: Path):
    campaign, units = _campaign(tmp_path)
    completion_path = _write_simulation_completion(campaign, units.iloc[0])
    completion = json.loads(completion_path.read_text(encoding="utf-8"))
    completion["schema"] = "gamma-smc.focused-selection-simulation/v4"
    completion_path.write_text(json.dumps(completion), encoding="utf-8")

    with pytest.raises(ValueError, match="schema/status"):
        generate_run_report(campaign)

    completion["schema"] = SIMULATION_SCHEMA_VERSION
    completion_path.write_text(json.dumps(completion), encoding="utf-8")
    (completion_path.parent / "simulation.trees").write_bytes(b"corrupt\n")
    with pytest.raises(ValueError, match="checksum"):
        generate_run_report(campaign)


def test_report_rejects_stale_decode_schema_and_corrupt_auxiliary(tmp_path: Path):
    campaign, units = _campaign(tmp_path)
    _write_simulation_completion(campaign, units.iloc[0])
    completion_path = _write_decode_completion(campaign, units.iloc[0])
    completion = json.loads(completion_path.read_text(encoding="utf-8"))
    completion["schema"] = "gamma-smc.focused-selection-decode/v0"
    completion_path.write_text(json.dumps(completion), encoding="utf-8")

    with pytest.raises(ValueError, match="schema/status"):
        generate_run_report(campaign)

    completion["schema"] = DECODE_SCHEMA_VERSION
    completion_path.write_text(json.dumps(completion), encoding="utf-8")
    (completion_path.parent / report.DECODE_OUTPUTS["decoder_stderr"]).write_bytes(
        b"corrupt\n"
    )
    with pytest.raises(ValueError, match="checksum"):
        generate_run_report(campaign)


def test_complete_analysis_report_summarizes_primary_contracts(tmp_path: Path):
    campaign, units = _campaign(tmp_path, neutral_replicates=2)
    _finalize_analyses(campaign, units)

    text = generate_run_report(campaign).read_text(encoding="utf-8")
    assert "Primary matched-neutral results" in text
    assert "Gamma-SMC | ancient_eurasia_han_introgression" in text
    assert "Tree truth | eas_phlash_median" in text
    assert "Neutral cross-fit calibration" in text
    assert "Gamma-SMC versus tree truth" in text
    assert "Empirical status: **not provided**" in text
    assert "No empirical selection p-value can be reported" in text
    assert "Q=5" in text
    assert "constant-s selection" in text
    assert "fixed s-by-AF endpoint from the disjoint terminal-AF calibration" in text
    assert "Q=10" not in text
    assert "logistic selection episode" not in text


def test_report_fails_closed_on_corrupt_finalized_table(tmp_path: Path):
    campaign, units = _campaign(tmp_path, neutral_replicates=2)
    _finalize_analyses(campaign, units)
    table = campaign / "results" / "analysis" / "cell_roc_auc.tsv"
    table.write_text("corrupt\n", encoding="utf-8")

    with pytest.raises(ValueError, match="checksum"):
        generate_run_report(campaign)


def test_report_rejects_manifested_but_incomplete_primary_grid(tmp_path: Path):
    campaign, units = _campaign(tmp_path, neutral_replicates=2)
    _finalize_analyses(campaign, units)
    analysis_dir = campaign / "results" / "analysis"
    table = analysis_dir / "cell_roc_auc.tsv"
    frame = pd.read_csv(table, sep="\t")
    primary_index = frame.index[
        (frame["metric"] == "p_tmrca_lt_threshold")
        & (frame["statistic"] == analysis.PRIMARY_STATISTIC)
    ][0]
    frame = frame.drop(index=primary_index).copy()
    frame.to_csv(table, sep="\t", index=False)
    completion_path = analysis_dir / "analysis_completion.json"
    completion = json.loads(completion_path.read_text(encoding="utf-8"))
    record = completion["outputs"]["cell_comparisons"]
    record["sha256"] = hashlib.sha256(table.read_bytes()).hexdigest()
    record["size_bytes"] = table.stat().st_size
    record["rows"] = len(frame)
    record["columns"] = list(frame.columns)
    completion_path.write_text(json.dumps(completion), encoding="utf-8")

    with pytest.raises(ValueError, match="primary grid is incomplete"):
        generate_run_report(campaign)
