"""Deterministic Markdown reporting for the focused selection campaign.

The report is intentionally usable while simulations are still running.  Unit
coverage and simulation QC are read directly from the immutable execution plan
and completed-unit records.  Inferential results are included only when their
analysis directory has a valid completion manifest.  Once an analysis claims
completion, every manifested output and the primary result grid are validated
before ``RUN_RESULTS.md`` is replaced.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .focused_selection_analysis import PRIMARY_STATISTIC, SCHEMA_VERSION
from .focused_selection_decode import SCHEMA_VERSION as DECODE_SCHEMA_VERSION
from .focused_selection_decode import sha256_file
from .focused_selection_simulation import (
    SCHEMA_VERSION as SIMULATION_SCHEMA_VERSION,
)

REPORT_SCHEMA_VERSION = "gamma-smc.focused-selection-report/v1"
ANALYSIS_TABLES = {
    "cell_comparisons": "cell_roc_auc.tsv",
    "selected_pvalues": "selected_vs_neutral_pvalues.tsv.gz",
    "conditional_power": "conditional_power_summary.tsv",
    "within_selected": "within_selected_summary.tsv",
    "unit_pair_tests": "within_selected_unit_pair_tests.tsv.gz",
    "neutral_crossfit_calibration": "neutral_crossfit_calibration.tsv",
    "primary_omnibus": "primary_7_threshold_omnibus.tsv.gz",
    "empirical_pvalues": "empirical_pvalues.tsv",
}
SCOPES = ("region_mean", "focal_nearest")
SIMULATION_OUTPUTS = {
    "tree": "simulation.trees",
    "pair_table": "sample_manifest.tsv",
    "overall_pairs": "pairs/overall.pairs.tsv",
    "truth_profiles": "truth_profiles.tsv.gz",
    "truth_class_summaries": "truth_class_summaries.tsv",
}
DECODE_OUTPUTS = {
    "spatial_profiles": "spatial_profiles.tsv.gz",
    "pair_summaries": "pair_summaries.tsv.gz",
    "class_summaries": "class_summaries.tsv",
    "overall_summary": "overall.summary.tsv",
    "run_manifest": "overall.summary.tsv.run.json",
    "decoder_stdout": "decoder.stdout.zst",
    "decoder_stderr": "decoder.stderr.zst",
}


def _canonical_sha256(value: Any) -> str:
    payload = json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _read_json(path: Path, label: str) -> Mapping[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"{label} is unreadable: {path}") from error
    if not isinstance(value, Mapping):
        raise TypeError(f"{label} is not a JSON object: {path}")
    return value


def _atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        temporary.write_text(text, encoding="utf-8", newline="\n")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _isclose(left: Any, right: Any) -> bool:
    try:
        return bool(np.isclose(float(left), float(right), rtol=0, atol=1e-12))
    except (TypeError, ValueError):
        return str(left) == str(right)


def _format_number(value: Any, digits: int = 3) -> str:
    if value is None:
        return "--"
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not math.isfinite(number):
        return "--"
    if number == 0:
        return "0"
    if abs(number) < 0.001:
        return f"{number:.2e}"
    return f"{number:.{digits}f}".rstrip("0").rstrip(".")


def _format_af(value: Any) -> str:
    return f"{100 * float(value):g}%"


def _format_s(value: Any) -> str:
    number = float(value)
    return "--" if np.isclose(number, 0) else f"{number:g}"


def _format_threshold(value: Any) -> str:
    years = float(value)
    return f"{years / 1_000:g}k"


def _han_model_contract_text(design: Mapping[str, Any]) -> str:
    simulation = design.get("simulation_contract")
    if not isinstance(simulation, Mapping):
        raise TypeError("focused study design lacks a simulation contract")
    episode = simulation.get("introgression_selection_episode")
    if not isinstance(episode, Mapping):
        raise TypeError("focused study design lacks the Han selection episode")
    try:
        scaling = float(simulation["slim_scaling_factor"])
        burn_in = float(simulation["slim_burn_in"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(
            "focused study design has invalid SLiM scaling/burn-in"
        ) from error
    if (
        not np.isfinite(scaling)
        or scaling <= 0
        or not np.isfinite(burn_in)
        or burn_in < 0
    ):
        raise ValueError("focused study design has invalid SLiM scaling/burn-in")
    cessation_rule = str(episode.get("cessation_rule", "")).strip()
    if not cessation_rule:
        raise ValueError("focused study design lacks the Han cessation rule")
    return (
        f"The Han selected model uses Q={scaling:g} and SLiM burn-in {burn_in:g}. "
        "It applies constant-s selection from introgression until a prespecified "
        "s- and AF-specific fixed endpoint, followed by neutrality. The recorded "
        f"cessation rule is: {cessation_rule}."
    )


def _markdown_table(headers: Sequence[str], rows: Sequence[Sequence[Any]]) -> str:
    def escape(value: Any) -> str:
        return str(value).replace("|", "\\|").replace("\n", " ")

    lines = [
        "| " + " | ".join(escape(value) for value in headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    lines.extend(
        "| " + " | ".join(escape(value) for value in row) + " |" for row in rows
    )
    return "\n".join(lines)


def _read_table(path: Path, label: str) -> pd.DataFrame:
    try:
        return pd.read_csv(path, sep="\t")
    except pd.errors.EmptyDataError:
        return pd.DataFrame()
    except Exception as error:
        raise ValueError(f"finalized table is unreadable: {label}") from error


def _safe_manifest_path(root: Path, relative: Any, label: str) -> Path:
    path = (root / str(relative)).resolve()
    try:
        path.relative_to(root.resolve())
    except ValueError as error:
        raise ValueError(
            f"finalized manifest path escapes its directory: {label}"
        ) from error
    return path


def _validate_output_record(root: Path, label: str, record: Mapping[str, Any]) -> Path:
    path = _safe_manifest_path(root, record.get("path", ""), label)
    if not path.is_file() or sha256_file(path) != record.get("sha256"):
        raise ValueError(f"finalized output checksum failed: {label}")
    if "size_bytes" in record and path.stat().st_size != int(record["size_bytes"]):
        raise ValueError(f"finalized output size failed: {label}")
    if "rows" in record:
        frame = _read_table(path, label)
        expected_rows = int(record["rows"])
        expected_columns = list(record.get("columns", []))
        if len(frame) != expected_rows or list(frame.columns) != expected_columns:
            raise ValueError(f"finalized output table schema failed: {label}")
    return path


def _validate_finalized_analysis(directory: Path) -> dict[str, Path] | None:
    """Return manifested output paths, or ``None`` for an unfinished analysis."""

    completion_path = directory / "analysis_completion.json"
    if not completion_path.is_file():
        return None
    completion = _read_json(completion_path, "analysis completion")
    if completion.get("schema") != SCHEMA_VERSION:
        raise ValueError(f"finalized analysis schema is incompatible: {directory}")
    if completion.get("status") != "complete":
        raise ValueError(f"finalized analysis status is not complete: {directory}")
    contract = completion.get("contract")
    if not isinstance(contract, Mapping) or completion.get(
        "contract_sha256"
    ) != _canonical_sha256(contract):
        raise ValueError(
            f"finalized analysis contract checksum is invalid: {directory}"
        )
    outputs = completion.get("outputs")
    if not isinstance(outputs, Mapping) or not outputs:
        raise ValueError(
            f"finalized analysis output manifest is malformed: {directory}"
        )
    paths: dict[str, Path] = {}
    for label, record in outputs.items():
        if not isinstance(record, Mapping):
            raise TypeError(f"finalized analysis output record is malformed: {label}")
        paths[str(label)] = _validate_output_record(directory, str(label), record)
    for label, filename in ANALYSIS_TABLES.items():
        if label not in paths or paths[label].name != filename:
            raise ValueError(f"finalized analysis lacks required table: {label}")
    paths["completion"] = completion_path
    return paths


def _validate_truth_gamma(directory: Path) -> dict[str, Path] | None:
    completion_path = directory / "truth_gamma_completion.json"
    if not completion_path.is_file():
        return None
    completion = _read_json(completion_path, "truth-Gamma completion")
    if completion.get("schema") != SCHEMA_VERSION:
        raise ValueError("finalized truth-Gamma schema is incompatible")
    if completion.get("status") != "complete":
        raise ValueError("finalized truth-Gamma status is not complete")
    contract = completion.get("contract")
    if not isinstance(contract, Mapping) or completion.get(
        "contract_sha256"
    ) != _canonical_sha256(contract):
        raise ValueError("finalized truth-Gamma contract checksum is invalid")
    outputs = completion.get("outputs")
    if not isinstance(outputs, Mapping):
        raise TypeError("finalized truth-Gamma output manifest is malformed")
    paths: dict[str, Path] = {}
    for label in ("paired", "summary"):
        record = outputs.get(label)
        if not isinstance(record, Mapping):
            raise TypeError(f"finalized truth-Gamma output is missing: {label}")
        paths[label] = _validate_output_record(directory, label, record)
    figures = outputs.get("figures", [])
    if not isinstance(figures, list):
        raise TypeError("finalized truth-Gamma figure manifest is malformed")
    for index, record in enumerate(figures):
        if not isinstance(record, Mapping):
            raise TypeError("finalized truth-Gamma figure record is malformed")
        _validate_output_record(directory, f"figure_{index}", record)
    paths["completion"] = completion_path
    return paths


def _validate_plan(campaign_dir: Path) -> tuple[pd.DataFrame, Mapping[str, Any]]:
    if any("onedrive" in part.casefold() for part in campaign_dir.resolve().parts):
        raise ValueError("focused report must not read or write a OneDrive campaign")
    plan_path = campaign_dir / "execution_units.tsv"
    design_path = campaign_dir / "study_design.json"
    if not plan_path.is_file() or not design_path.is_file():
        raise ValueError("focused execution plan and study design are required")
    try:
        units = pd.read_csv(plan_path, sep="\t")
    except Exception as error:
        raise ValueError("focused execution plan is unreadable") from error
    required = {
        "unit_id",
        "demography_id",
        "simulation_class",
        "selection_coefficient",
        "target_allele_frequency",
        "exact_sample_alt_count",
        "sample_diploids",
        "replicate_index",
    }
    missing = sorted(required.difference(units.columns))
    if missing or units.empty or units["unit_id"].duplicated().any():
        detail = ", ".join(missing) if missing else "empty or duplicate unit IDs"
        raise ValueError(f"focused execution plan is invalid: {detail}")
    design = _read_json(design_path, "focused study design")
    expected_plan_sha = design.get("hashes", {}).get("execution_units_tsv_sha256")
    actual_plan_sha = hashlib.sha256(
        plan_path.read_bytes().replace(b"\r\n", b"\n")
    ).hexdigest()
    if actual_plan_sha != expected_plan_sha:
        raise ValueError("focused execution plan checksum is invalid")
    scientific_grid = design.get("scientific_grid", {})
    if int(scientific_grid.get("total_units", -1)) != len(units):
        raise ValueError("focused study design total does not match execution plan")
    return units.sort_values("unit_id").reset_index(drop=True), design


def _unit_matches_plan(recorded: Mapping[str, Any], planned: Mapping[str, Any]) -> bool:
    for key, expected in planned.items():
        if key not in recorded or not _isclose(recorded[key], expected):
            return False
    return True


def _validate_simulation_completion(
    path: Path, planned: Mapping[str, Any]
) -> dict[str, Any]:
    completion = _read_json(path, "simulation completion")
    contract = completion.get("contract")
    if (
        completion.get("schema") != SIMULATION_SCHEMA_VERSION
        or completion.get("status") != "complete"
        or not isinstance(contract, Mapping)
    ):
        raise ValueError(
            f"simulation completion schema/status is incompatible: {planned['unit_id']}"
        )
    if completion.get("contract_sha256") != _canonical_sha256(contract):
        raise ValueError(
            f"simulation completion contract checksum is invalid: {planned['unit_id']}"
        )
    recorded_unit = contract.get("unit")
    if not isinstance(recorded_unit, Mapping) or not _unit_matches_plan(
        recorded_unit, planned
    ):
        raise ValueError(
            f"simulation completion unit disagrees with plan: {planned['unit_id']}"
        )
    contract_path = path.parent / "simulation_contract.json"
    disk_contract = _read_json(contract_path, "simulation contract")
    if _canonical_sha256(disk_contract) != _canonical_sha256(contract):
        raise ValueError(
            f"simulation contract file disagrees with completion: {planned['unit_id']}"
        )
    outputs = completion.get("outputs")
    if not isinstance(outputs, Mapping) or set(outputs) != set(SIMULATION_OUTPUTS):
        raise ValueError(
            f"simulation output manifest is malformed: {planned['unit_id']}"
        )
    for label, relative in SIMULATION_OUTPUTS.items():
        record = outputs.get(label)
        if not isinstance(record, Mapping):
            raise TypeError(
                f"simulation output is missing for {planned['unit_id']}: {label}"
            )
        output_path = _validate_output_record(
            path.parent, f"{planned['unit_id']}:{label}", record
        )
        if output_path != (path.parent / relative).resolve():
            raise ValueError(
                f"simulation output path is invalid for {planned['unit_id']}: {label}"
            )
    sample_diploids = int(planned["sample_diploids"])
    expected_alt = int(planned["exact_sample_alt_count"])
    alt_count = int(completion.get("sample_alt_count", -1))
    sample_af = float(completion.get("sample_af", np.nan))
    if alt_count != expected_alt or not np.isclose(
        sample_af, expected_alt / (2 * sample_diploids), rtol=0, atol=1e-12
    ):
        raise ValueError(
            f"completed sample AF fails exact target: {planned['unit_id']}"
        )
    counts = completion.get("genotype_counts")
    if not isinstance(counts, Mapping) or set(counts) != {
        "hom_ref",
        "heterozygous",
        "hom_alt",
    }:
        raise ValueError(
            f"completed genotype counts are malformed: {planned['unit_id']}"
        )
    parsed_counts = {key: int(counts[key]) for key in counts}
    if (
        sum(parsed_counts.values()) != sample_diploids
        or parsed_counts["heterozygous"] + 2 * parsed_counts["hom_alt"] != alt_count
        or parsed_counts["hom_ref"] < 2
        or parsed_counts["heterozygous"] < 1
        or parsed_counts["hom_alt"] < 2
    ):
        raise ValueError(f"completed genotype counts fail QC: {planned['unit_id']}")
    attempts = int(completion.get("attempts_completed", -1))
    elapsed = float(completion.get("elapsed_seconds", np.nan))
    if attempts < 1 or not np.isfinite(elapsed) or elapsed < 0:
        raise ValueError(
            f"completed runtime/attempt record is invalid: {planned['unit_id']}"
        )
    return {
        "sample_alt_count": alt_count,
        "sample_af": sample_af,
        "hom_ref": parsed_counts["hom_ref"],
        "heterozygous": parsed_counts["heterozygous"],
        "hom_alt": parsed_counts["hom_alt"],
        "attempts_completed": attempts,
        "elapsed_seconds": elapsed,
    }


def _validate_decode_completion(path: Path, planned: Mapping[str, Any]) -> None:
    unit_id = str(planned["unit_id"])
    completion = _read_json(path, "decode completion")
    contract = completion.get("contract")
    if (
        completion.get("schema") != DECODE_SCHEMA_VERSION
        or completion.get("status") != "complete"
        or not isinstance(contract, Mapping)
    ):
        raise ValueError(f"decode completion schema/status is incompatible: {unit_id}")
    if completion.get("contract_sha256") != _canonical_sha256(contract):
        raise ValueError(f"decode completion contract checksum is invalid: {unit_id}")
    recorded_unit = contract.get("static", {}).get("unit_record", {})
    if not isinstance(recorded_unit, Mapping) or not _unit_matches_plan(
        recorded_unit, planned
    ):
        raise ValueError(f"decode completion unit disagrees with plan: {unit_id}")
    outputs = completion.get("outputs")
    if not isinstance(outputs, Mapping) or set(outputs) != set(DECODE_OUTPUTS):
        raise ValueError(f"decode output manifest is malformed: {unit_id}")
    for label, filename in DECODE_OUTPUTS.items():
        record = outputs.get(label)
        if not isinstance(record, Mapping):
            raise TypeError(f"decode output is missing for {unit_id}: {label}")
        output_path = _validate_output_record(path.parent, f"{unit_id}:{label}", record)
        if output_path != (path.parent / filename).resolve():
            raise ValueError(f"decode output path is invalid for {unit_id}: {label}")


def _collect_unit_status(campaign_dir: Path, units: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for planned in units.to_dict(orient="records"):
        unit_id = str(planned["unit_id"])
        unit_dir = campaign_dir / "work" / unit_id
        simulation_completion = unit_dir / "simulation_complete.json"
        failure = unit_dir / "simulation_failed.json"
        decoded_completion = unit_dir / "decoded" / "completion.json"
        detail: dict[str, Any] = {}
        if simulation_completion.is_file():
            detail = _validate_simulation_completion(simulation_completion, planned)
            simulation_status = "complete"
        elif failure.is_file():
            simulation_status = "failed_or_exhausted"
        elif unit_dir.exists():
            simulation_status = "partial"
        else:
            simulation_status = "not_started"
        if decoded_completion.is_file():
            if simulation_status != "complete":
                raise ValueError(f"decoded unit lacks completed simulation: {unit_id}")
            _validate_decode_completion(decoded_completion, planned)
            decode_status = "complete"
        elif (unit_dir / "decoded").exists():
            decode_status = "partial"
        else:
            decode_status = "not_started"
        rows.append(
            {
                **planned,
                "simulation_status": simulation_status,
                "decode_status": decode_status,
                **detail,
            }
        )
    return pd.DataFrame(rows)


def _expected_primary_grid(
    units: pd.DataFrame, thresholds: Sequence[float]
) -> set[tuple[Any, ...]]:
    selected = units[units["simulation_class"].astype(str) == "selected"]
    cells = selected[
        ["demography_id", "target_allele_frequency", "selection_coefficient"]
    ].drop_duplicates()
    return {
        (
            str(row.demography_id),
            float(row.target_allele_frequency),
            float(row.selection_coefficient),
            scope,
            float(threshold),
        )
        for row in cells.itertuples(index=False)
        for scope in SCOPES
        for threshold in thresholds
    }


def _primary(frame: pd.DataFrame) -> pd.DataFrame:
    required = {"metric", "statistic", "threshold_years"}
    if not required.issubset(frame.columns):
        raise ValueError("finalized analysis table lacks primary-statistic columns")
    return frame[
        (frame["metric"].astype(str) == "p_tmrca_lt_threshold")
        & (frame["statistic"].astype(str) == PRIMARY_STATISTIC)
    ].copy()


def _validate_primary_analysis(
    paths: Mapping[str, Path],
    units: pd.DataFrame,
    thresholds: Sequence[float],
    *,
    require_pair_tests: bool,
) -> dict[str, pd.DataFrame]:
    frames = {label: _read_table(paths[label], label) for label in ANALYSIS_TABLES}
    expected = _expected_primary_grid(units, thresholds)
    selected_units = units[units["simulation_class"].astype(str) == "selected"]
    neutral_units = units[units["simulation_class"].astype(str) == "neutral"]
    completion = _read_json(paths["completion"], "analysis completion")
    contract = completion.get("contract", {})
    parameters = contract.get("parameters", {})
    inputs = contract.get("inputs", {})
    selected_cell_sizes = set(
        selected_units.groupby(
            ["demography_id", "target_allele_frequency", "selection_coefficient"]
        ).size()
    )
    neutral_cell_sizes = set(
        neutral_units.groupby(["demography_id", "target_allele_frequency"]).size()
    )
    if len(selected_cell_sizes) != 1 or len(neutral_cell_sizes) != 1:
        raise ValueError("focused plan has unequal per-cell replicate counts")
    if (
        tuple(float(value) for value in parameters.get("thresholds_years", []))
        != tuple(float(value) for value in thresholds)
        or int(parameters.get("minimum_selected_per_cell", -1))
        != int(next(iter(selected_cell_sizes)))
        or int(parameters.get("minimum_neutral_per_cell", -1))
        != int(next(iter(neutral_cell_sizes)))
        or str(parameters.get("primary_statistic", "")) != PRIMARY_STATISTIC
    ):
        raise ValueError("finalized analysis contract disagrees with focused plan")
    pair_input_present = bool(inputs.get("pair_summaries_path")) and bool(
        inputs.get("pair_summaries_sha256")
    )
    if pair_input_present != bool(require_pair_tests):
        raise ValueError("finalized analysis pair-input contract is inconsistent")

    for label in ("cell_comparisons", "conditional_power", "within_selected"):
        primary = _primary(frames[label])
        key_columns = [
            "demography_id",
            "target_allele_frequency",
            "selection_coefficient",
            "statistic_scope",
            "threshold_years",
        ]
        if (
            not set(key_columns).issubset(primary.columns)
            or primary.duplicated(key_columns).any()
        ):
            raise ValueError(
                f"finalized primary grid is duplicated or malformed: {label}"
            )
        observed = {
            (str(row[0]), float(row[1]), float(row[2]), str(row[3]), float(row[4]))
            for row in primary[key_columns].itertuples(index=False, name=None)
        }
        if observed != expected:
            raise ValueError(f"finalized primary grid is incomplete: {label}")

    comparisons = _primary(frames["cell_comparisons"])
    for column in (
        "roc_auc",
        "roc_auc_ci95_low",
        "roc_auc_ci95_high",
        "auc_unit_label_permutation_p_upper",
    ):
        values = pd.to_numeric(comparisons[column], errors="coerce")
        if not np.isfinite(values).all() or ((values < 0) | (values > 1)).any():
            raise ValueError(f"finalized primary comparison is invalid: {column}")
    expected_selected_counts = (
        selected_units.groupby(
            ["demography_id", "target_allele_frequency", "selection_coefficient"]
        )
        .size()
        .to_dict()
    )
    expected_neutral_counts = (
        neutral_units.groupby(["demography_id", "target_allele_frequency"])
        .size()
        .to_dict()
    )
    for row in comparisons.itertuples(index=False):
        selected_key = (
            str(row.demography_id),
            float(row.target_allele_frequency),
            float(row.selection_coefficient),
        )
        neutral_key = (str(row.demography_id), float(row.target_allele_frequency))
        if int(row.n_selected) != int(expected_selected_counts[selected_key]) or int(
            row.n_neutral
        ) != int(expected_neutral_counts[neutral_key]):
            raise ValueError(
                "finalized primary comparison unit counts disagree with plan"
            )

    power = _primary(frames["conditional_power"])
    for column in ("conditional_power_raw", "conditional_power_bh_within_7"):
        values = pd.to_numeric(power[column], errors="coerce")
        if not np.isfinite(values).all() or ((values < 0) | (values > 1)).any():
            raise ValueError(f"finalized conditional power is invalid: {column}")

    selected_pvalues = _primary(frames["selected_pvalues"])
    selected_key = ["unit_id", "statistic_scope", "threshold_years"]
    expected_selected_pvalues = {
        (str(row.unit_id), scope, float(threshold))
        for row in selected_units.itertuples(index=False)
        for scope in SCOPES
        for threshold in thresholds
    }
    observed_selected_pvalues = {
        (str(row[0]), str(row[1]), float(row[2]))
        for row in selected_pvalues[selected_key].itertuples(index=False, name=None)
    }
    if (
        selected_pvalues.duplicated(selected_key).any()
        or observed_selected_pvalues != expected_selected_pvalues
    ):
        raise ValueError("finalized selected-unit p-value grid is incomplete")
    for column in ("mc_p_upper", "bh_q_upper"):
        values = pd.to_numeric(selected_pvalues[column], errors="coerce")
        if not np.isfinite(values).all() or ((values < 0) | (values > 1)).any():
            raise ValueError(f"finalized selected-unit p-values are invalid: {column}")

    omnibus = frames["primary_omnibus"].copy()
    required_omnibus = {
        "unit_id",
        "demography_id",
        "target_allele_frequency",
        "selection_coefficient",
        "statistic_scope",
        "omnibus_p_upper",
        "n_thresholds",
    }
    if (
        not required_omnibus.issubset(omnibus.columns)
        or omnibus.duplicated(["unit_id", "statistic_scope"]).any()
    ):
        raise ValueError("finalized primary omnibus table is malformed")
    expected_omnibus = {
        (str(row.unit_id), scope)
        for row in selected_units.itertuples(index=False)
        for scope in SCOPES
    }
    observed_omnibus = set(
        omnibus[["unit_id", "statistic_scope"]]
        .astype(str)
        .itertuples(index=False, name=None)
    )
    if observed_omnibus != expected_omnibus or set(
        omnibus["n_thresholds"].astype(int)
    ) != {len(thresholds)}:
        raise ValueError("finalized primary omnibus grid is incomplete")
    omnibus_pvalues = pd.to_numeric(omnibus["omnibus_p_upper"], errors="coerce")
    if (
        not np.isfinite(omnibus_pvalues).all()
        or ((omnibus_pvalues < 0) | (omnibus_pvalues > 1)).any()
    ):
        raise ValueError("finalized primary omnibus p-values are invalid")

    crossfit = _primary(frames["neutral_crossfit_calibration"])
    crossfit_key = [
        "demography_id",
        "target_allele_frequency",
        "statistic_scope",
        "threshold_years",
    ]
    expected_crossfit = {
        (
            str(row.demography_id),
            float(row.target_allele_frequency),
            scope,
            float(threshold),
        )
        for row in neutral_units[["demography_id", "target_allele_frequency"]]
        .drop_duplicates()
        .itertuples(index=False)
        for scope in SCOPES
        for threshold in thresholds
    }
    observed_crossfit = {
        (str(row[0]), float(row[1]), str(row[2]), float(row[3]))
        for row in crossfit[crossfit_key].itertuples(index=False, name=None)
    }
    if (
        crossfit.duplicated(crossfit_key).any()
        or observed_crossfit != expected_crossfit
    ):
        raise ValueError("finalized neutral cross-fit grid is incomplete")
    for column in ("mean_crossfit_p", "fraction_crossfit_p_le_0_05", "ks_uniform_p"):
        values = pd.to_numeric(crossfit[column], errors="coerce")
        if not np.isfinite(values).all() or ((values < 0) | (values > 1)).any():
            raise ValueError(
                f"finalized neutral cross-fit diagnostic is invalid: {column}"
            )
    for row in crossfit.itertuples(index=False):
        neutral_key = (str(row.demography_id), float(row.target_allele_frequency))
        if int(row.n_neutral) != int(expected_neutral_counts[neutral_key]):
            raise ValueError("finalized neutral cross-fit count disagrees with plan")

    raw_pair_tests = frames["unit_pair_tests"]
    pair_tests = (
        _primary(raw_pair_tests) if not raw_pair_tests.empty else raw_pair_tests.copy()
    )
    if require_pair_tests:
        pair_key = ["unit_id", "statistic_scope", "threshold_years"]
        observed_pair_tests = {
            (str(row[0]), str(row[1]), float(row[2]))
            for row in pair_tests[pair_key].itertuples(index=False, name=None)
        }
        if (
            pair_tests.duplicated(pair_key).any()
            or observed_pair_tests != expected_selected_pvalues
        ):
            raise ValueError("finalized within-unit pair-test grid is incomplete")
    elif not pair_tests.empty:
        raise ValueError("analysis without pair input unexpectedly contains pair tests")
    return frames


def _series_by_threshold(group: pd.DataFrame, column: str) -> str:
    ordered = group.sort_values("threshold_years")
    return "; ".join(
        f"{_format_threshold(threshold)}:{_format_number(value)}"
        for threshold, value in ordered[["threshold_years", column]].itertuples(
            index=False, name=None
        )
    )


def _coverage_section(status: pd.DataFrame) -> str:
    group_columns = [
        "demography_id",
        "target_allele_frequency",
        "simulation_class",
        "selection_coefficient",
    ]
    rows: list[list[str]] = []
    for labels, group in status.groupby(group_columns, sort=True, dropna=False):
        counts = group["simulation_status"].value_counts()
        decode_counts = group["decode_status"].value_counts()
        rows.append(
            [
                str(labels[0]),
                _format_af(labels[1]),
                str(labels[2]),
                _format_s(labels[3]),
                str(len(group)),
                str(int(counts.get("complete", 0))),
                str(int(counts.get("failed_or_exhausted", 0))),
                str(int(counts.get("partial", 0))),
                str(int(counts.get("not_started", 0))),
                str(int(decode_counts.get("complete", 0))),
            ]
        )
    return _markdown_table(
        [
            "Demography",
            "AF",
            "Class",
            "s",
            "Planned",
            "Sim complete",
            "Failed",
            "Partial",
            "Not started",
            "Decode complete",
        ],
        rows,
    )


def _qc_runtime_section(status: pd.DataFrame) -> tuple[str, str]:
    complete = status[status["simulation_status"] == "complete"].copy()
    if complete.empty:
        return (
            "No completed simulations are available for exact-AF/genotype QC.",
            "No completed simulations are available for runtime summaries.",
        )
    group_columns = [
        "demography_id",
        "target_allele_frequency",
        "simulation_class",
        "selection_coefficient",
    ]
    qc_rows: list[list[str]] = []
    runtime_rows: list[list[str]] = []
    for labels, group in complete.groupby(group_columns, sort=True, dropna=False):
        expected_alt = int(group["exact_sample_alt_count"].iloc[0])
        qc_rows.append(
            [
                str(labels[0]),
                _format_af(labels[1]),
                str(labels[2]),
                _format_s(labels[3]),
                str(len(group)),
                str(expected_alt),
                _format_number(group["sample_af"].min()),
                _format_number(group["sample_af"].max()),
                f"{int(group['hom_ref'].min())}--{int(group['hom_ref'].max())}",
                f"{int(group['heterozygous'].min())}--{int(group['heterozygous'].max())}",
                f"{int(group['hom_alt'].min())}--{int(group['hom_alt'].max())}",
                "pass",
            ]
        )
        runtime_rows.append(
            [
                str(labels[0]),
                _format_af(labels[1]),
                str(labels[2]),
                _format_s(labels[3]),
                str(len(group)),
                _format_number(group["attempts_completed"].median(), 1),
                f"{int(group['attempts_completed'].min())}--{int(group['attempts_completed'].max())}",
                _format_number(group["elapsed_seconds"].median(), 1),
                _format_number(group["elapsed_seconds"].sum() / 60, 1),
            ]
        )
    qc = _markdown_table(
        [
            "Demography",
            "AF",
            "Class",
            "s",
            "n",
            "Exact ALT copies",
            "AF min",
            "AF max",
            "ref/ref range",
            "ref/alt range",
            "alt/alt range",
            "QC",
        ],
        qc_rows,
    )
    runtime = _markdown_table(
        [
            "Demography",
            "AF",
            "Class",
            "s",
            "n",
            "Median attempts",
            "Attempt range",
            "Median seconds",
            "Sum completed unit-minutes",
        ],
        runtime_rows,
    )
    return qc, runtime


def _primary_results_section(source: str, frames: Mapping[str, pd.DataFrame]) -> str:
    comparisons = _primary(frames["cell_comparisons"])
    power = _primary(frames["conditional_power"])
    merged = comparisons.merge(
        power[
            [
                "demography_id",
                "target_allele_frequency",
                "selection_coefficient",
                "statistic_scope",
                "threshold_years",
                "conditional_power_raw",
                "conditional_power_bh_within_7",
            ]
        ],
        on=[
            "demography_id",
            "target_allele_frequency",
            "selection_coefficient",
            "statistic_scope",
            "threshold_years",
        ],
        how="inner",
        validate="one_to_one",
    )
    selected_pvalues = _primary(frames["selected_pvalues"])
    pvalue_summary = (
        selected_pvalues.groupby(
            [
                "demography_id",
                "target_allele_frequency",
                "selection_coefficient",
                "statistic_scope",
                "threshold_years",
            ],
            sort=True,
        )["mc_p_upper"]
        .agg(["median", "min"])
        .reset_index()
        .rename(columns={"median": "median_unit_p", "min": "minimum_unit_p"})
    )
    merged = merged.merge(
        pvalue_summary,
        on=[
            "demography_id",
            "target_allele_frequency",
            "selection_coefficient",
            "statistic_scope",
            "threshold_years",
        ],
        how="inner",
        validate="one_to_one",
    )
    omnibus = frames["primary_omnibus"]
    rows: list[list[str]] = []
    group_columns = [
        "demography_id",
        "target_allele_frequency",
        "selection_coefficient",
        "statistic_scope",
    ]
    for labels, group in merged.groupby(group_columns, sort=True):
        unit_omnibus = omnibus[
            (omnibus["demography_id"].astype(str) == str(labels[0]))
            & np.isclose(omnibus["target_allele_frequency"], float(labels[1]))
            & np.isclose(omnibus["selection_coefficient"], float(labels[2]))
            & (omnibus["statistic_scope"].astype(str) == str(labels[3]))
        ]
        rows.append(
            [
                source,
                str(labels[0]),
                _format_af(labels[1]),
                _format_s(labels[2]),
                str(labels[3]),
                _series_by_threshold(group, "roc_auc"),
                _series_by_threshold(group, "auc_unit_label_permutation_p_upper"),
                "; ".join(
                    f"{_format_threshold(threshold)}:{_format_number(median)} "
                    f"({_format_number(minimum)})"
                    for threshold, median, minimum in group.sort_values(
                        "threshold_years"
                    )[
                        ["threshold_years", "median_unit_p", "minimum_unit_p"]
                    ].itertuples(index=False, name=None)
                ),
                _series_by_threshold(group, "conditional_power_raw"),
                _series_by_threshold(group, "conditional_power_bh_within_7"),
                (
                    f"median={_format_number(unit_omnibus['omnibus_p_upper'].median())}; "
                    f"p<=.05={int((unit_omnibus['omnibus_p_upper'] <= 0.05).sum())}/"
                    f"{len(unit_omnibus)}"
                ),
            ]
        )
    return _markdown_table(
        [
            "Source",
            "Demography",
            "AF",
            "s",
            "Scope",
            "AUC by x",
            "AUC permutation p by x",
            "Median unit p (min) by x",
            "Raw power by x",
            "BH power by x",
            "7-x omnibus",
        ],
        rows,
    )


def _within_section(frames: Mapping[str, pd.DataFrame]) -> str:
    within = _primary(frames["within_selected"])
    pair_tests = _primary(frames["unit_pair_tests"])
    rows: list[list[str]] = []
    group_columns = [
        "demography_id",
        "target_allele_frequency",
        "selection_coefficient",
        "statistic_scope",
    ]
    for labels, group in within.groupby(group_columns, sort=True):
        tests = pair_tests[
            (pair_tests["demography_id"].astype(str) == str(labels[0]))
            & np.isclose(pair_tests["target_allele_frequency"], float(labels[1]))
            & np.isclose(pair_tests["selection_coefficient"], float(labels[2]))
            & (pair_tests["statistic_scope"].astype(str) == str(labels[3]))
        ]
        auc_text = "--"
        p_text = "--"
        if not tests.empty:
            auc_text = _format_number(tests["within_pair_roc_auc"].median())
            p_text = (
                f"{int((tests['genotype_label_permutation_p_upper'] <= 0.05).sum())}/"
                f"{len(tests)}"
            )
        rows.append(
            [
                str(labels[0]),
                _format_af(labels[1]),
                _format_s(labels[2]),
                str(labels[3]),
                _series_by_threshold(group, "mean_within_selection_score"),
                _series_by_threshold(group, "sign_test_p_upper"),
                auc_text,
                p_text,
            ]
        )
    return _markdown_table(
        [
            "Demography",
            "AF",
            "s",
            "Scope",
            "Mean D by x",
            "Sign-test p by x",
            "Median within-pair AUC",
            "Pair permutation p<=.05",
        ],
        rows,
    )


def _crossfit_section(frames: Mapping[str, pd.DataFrame]) -> str:
    crossfit = _primary(frames["neutral_crossfit_calibration"])
    rows: list[list[str]] = []
    group_columns = ["demography_id", "target_allele_frequency", "statistic_scope"]
    for labels, group in crossfit.groupby(group_columns, sort=True):
        rows.append(
            [
                str(labels[0]),
                _format_af(labels[1]),
                str(labels[2]),
                _series_by_threshold(group, "mean_crossfit_p"),
                _series_by_threshold(group, "fraction_crossfit_p_le_0_05"),
                _series_by_threshold(group, "ks_uniform_p"),
            ]
        )
    return _markdown_table(
        [
            "Demography",
            "AF",
            "Scope",
            "Mean held-out p by x",
            "Type-I fraction p<=.05 by x",
            "KS-uniform p by x",
        ],
        rows,
    )


def _truth_gamma_section(
    path: Path, units: pd.DataFrame, thresholds: Sequence[float]
) -> str:
    frame = _read_table(path, "truth-Gamma accuracy")
    required = {
        "simulation_class",
        "demography_id",
        "target_allele_frequency",
        "selection_coefficient",
        "statistic_scope",
        "metric",
        "statistic",
        "threshold_years",
        "n_units",
        "bias_gamma_minus_truth",
        "mae",
        "pearson_r",
        "spearman_rho",
        "sign_agreement_fraction",
    }
    if not required.issubset(frame.columns):
        raise ValueError("finalized truth-Gamma accuracy table is malformed")
    primary = frame[
        (frame["simulation_class"].astype(str) == "selected")
        & (frame["metric"].astype(str) == "p_tmrca_lt_threshold")
        & (frame["statistic"].astype(str) == PRIMARY_STATISTIC)
    ].copy()
    key_columns = [
        "demography_id",
        "target_allele_frequency",
        "selection_coefficient",
        "statistic_scope",
        "threshold_years",
    ]
    observed = {
        (str(row[0]), float(row[1]), float(row[2]), str(row[3]), float(row[4]))
        for row in primary[key_columns].itertuples(index=False, name=None)
    }
    if primary.duplicated(key_columns).any() or observed != _expected_primary_grid(
        units, thresholds
    ):
        raise ValueError("finalized truth-Gamma primary accuracy grid is incomplete")
    selected_counts = (
        units[units["simulation_class"].astype(str) == "selected"]
        .groupby(["demography_id", "target_allele_frequency", "selection_coefficient"])
        .size()
        .to_dict()
    )
    for row in primary.itertuples(index=False):
        key = (
            str(row.demography_id),
            float(row.target_allele_frequency),
            float(row.selection_coefficient),
        )
        if int(row.n_units) != int(selected_counts[key]):
            raise ValueError("finalized truth-Gamma unit count disagrees with plan")
    mae = pd.to_numeric(primary["mae"], errors="coerce")
    sign = pd.to_numeric(primary["sign_agreement_fraction"], errors="coerce")
    if (
        not np.isfinite(mae).all()
        or (mae < 0).any()
        or not np.isfinite(sign).all()
        or ((sign < 0) | (sign > 1)).any()
    ):
        raise ValueError("finalized truth-Gamma accuracy values are invalid")
    rows: list[list[str]] = []
    group_columns = [
        "demography_id",
        "target_allele_frequency",
        "selection_coefficient",
        "statistic_scope",
    ]
    for labels, group in primary.groupby(group_columns, sort=True):
        rows.append(
            [
                str(labels[0]),
                _format_af(labels[1]),
                _format_s(labels[2]),
                str(labels[3]),
                _series_by_threshold(group, "bias_gamma_minus_truth"),
                _series_by_threshold(group, "mae"),
                _series_by_threshold(group, "spearman_rho"),
                _series_by_threshold(group, "sign_agreement_fraction"),
            ]
        )
    return _markdown_table(
        [
            "Demography",
            "AF",
            "s",
            "Scope",
            "Gamma-truth bias by x",
            "MAE by x",
            "Spearman rho by x",
            "Sign agreement by x",
        ],
        rows,
    )


def _relative_link(report_path: Path, target: Path, label: str) -> str:
    relative = os.path.relpath(target, report_path.parent).replace("\\", "/")
    return f"[{label}]({relative})"


def generate_run_report(
    campaign_dir: str | Path,
    *,
    output_path: str | Path | None = None,
) -> Path:
    """Validate available campaign state and atomically write ``RUN_RESULTS.md``."""

    campaign = Path(campaign_dir).resolve()
    report_path = (
        Path(output_path).resolve()
        if output_path is not None
        else campaign / "RUN_RESULTS.md"
    )
    try:
        report_path.relative_to(campaign)
    except ValueError as error:
        raise ValueError(
            "focused report output must remain inside the campaign"
        ) from error
    units, design = _validate_plan(campaign)
    status = _collect_unit_status(campaign, units)
    thresholds = tuple(
        float(value) for value in design.get("time", {}).get("thresholds_years", [])
    )
    if len(thresholds) != 7 or sorted(thresholds) != list(thresholds):
        raise ValueError("focused study design must contain seven ordered thresholds")

    analysis_root = campaign / "results" / "analysis"
    truth_root = campaign / "results" / "analysis_tree_truth"
    comparison_root = campaign / "results" / "truth_gamma_comparison"
    gamma_paths = _validate_finalized_analysis(analysis_root)
    truth_paths = _validate_finalized_analysis(truth_root)
    comparison_paths = _validate_truth_gamma(comparison_root)
    gamma_frames = (
        _validate_primary_analysis(
            gamma_paths, units, thresholds, require_pair_tests=True
        )
        if gamma_paths is not None
        else None
    )
    truth_frames = (
        _validate_primary_analysis(
            truth_paths, units, thresholds, require_pair_tests=False
        )
        if truth_paths is not None
        else None
    )

    completed_simulations = int((status["simulation_status"] == "complete").sum())
    completed_decodes = int((status["decode_status"] == "complete").sum())
    analysis_complete = all(
        value is not None for value in (gamma_paths, truth_paths, comparison_paths)
    )
    campaign_complete = (
        completed_simulations == len(units)
        and completed_decodes == len(units)
        and analysis_complete
    )
    qc, runtime = _qc_runtime_section(status)
    empirical_status = "unavailable"
    empirical_rows = 0
    if gamma_frames is not None:
        empirical_rows = len(gamma_frames["empirical_pvalues"])
        empirical_status = "provided" if empirical_rows else "not provided"

    lines = [
        "# Focused EAS Gamma-SMC simulation results",
        "",
        f"Report status: **{'complete' if campaign_complete else 'partial'}**.",
        "",
        (
            f"Validated coverage is {completed_simulations}/{len(units)} simulations and "
            f"{completed_decodes}/{len(units)} Gamma-SMC decodes. "
            "This file is deterministic for the validated campaign state and contains "
            "no wall-clock generation timestamp."
        ),
        "",
        "## Fixed study contract",
        "",
        (
            "The plan compares EAS PHLASH-median no-introgression and Han "
            "AncientEurasia_9K19 introgression demographies at exact sampled AF "
            "10%, 20%, and 30%. Selected cells use s=0.01 and s=0.005 with 10 "
            "accepted replicates; each demography-by-AF cell has 100 independently "
            "seeded neutral replicates. Sequences are 10 Mb with a focal site at 5 "
            "Mb, and x is 1, 4.5, 10, 20, 30, 40, or 50 kya at 25 years/generation."
        ),
        "",
        (
            "The primary selection-directed statistic is D_x = P(TMRCA<x | alt/alt) "
            "- P(TMRCA<x | ref/ref). Larger D_x is younger and more sweep-like. The "
            "one-sided simulation p-value is the upper-tail matched-neutral value with "
            "the +1 correction; the smallest attainable value with 100 nulls is 1/101."
        ),
        "",
        "## Coverage",
        "",
        _coverage_section(status),
        "",
        "## Exact AF and genotype QC",
        "",
        qc,
        "",
        (
            "Every row marked pass has exactly the planned 20, 40, or 60 alternate "
            "copies among 200 sampled haplotypes and at least two ref/ref, one ref/alt, "
            "and two alt/alt diploids. Ranges summarize completed units only."
        ),
        "",
        "## Simulation attempts and runtime",
        "",
        runtime,
        "",
        (
            "Elapsed time is the completed unit's recorded simulation time. The total "
            "does not include failed/partial attempts from units without a completion "
            "record and is therefore not total campaign wall time."
        ),
        "",
        "## Primary matched-neutral results",
        "",
    ]
    if gamma_frames is None:
        lines.extend(
            [
                (
                    "Gamma-SMC analysis is not finalized yet. Any unmanifested partial "
                    "tables are intentionally ignored."
                ),
                "",
            ]
        )
    else:
        lines.extend(
            [
                _primary_results_section("Gamma-SMC", gamma_frames),
                "",
                (
                    "AUC is Pr(D_selected > D_neutral), with ties receiving half "
                    "weight. The AUC p-value is the seeded unit-label permutation "
                    "test. Unit p-values are one-sided matched-neutral values with "
                    "the +1 correction; their median and minimum are shown. Raw "
                    "power is the fraction of the 10 accepted selected replicates "
                    "with p<=0.05; BH power uses the within-seven-threshold adjustment."
                ),
                "",
            ]
        )
    lines.extend(["## Tree-truth primary results", ""])
    if truth_frames is None:
        lines.extend(
            [
                (
                    "Tree-truth analysis is not finalized yet. Any unmanifested partial "
                    "tables are intentionally ignored."
                ),
                "",
            ]
        )
    else:
        lines.extend([_primary_results_section("Tree truth", truth_frames), ""])
    lines.extend(["## Within-selected-simulation association", ""])
    if gamma_frames is None:
        lines.extend(["Within-simulation Gamma-SMC tests are not finalized yet.", ""])
    else:
        lines.extend(
            [
                _within_section(gamma_frames),
                "",
                (
                    "Within-pair label permutations test genotype association inside an "
                    "accepted realization. They are not demographic selection p-values "
                    "and should not replace the matched-neutral result."
                ),
                "",
            ]
        )
    lines.extend(["## Neutral cross-fit calibration", ""])
    if gamma_frames is None:
        lines.extend(["Held-out neutral calibration is not finalized yet.", ""])
    else:
        lines.extend(
            [
                _crossfit_section(gamma_frames),
                "",
                (
                    "These are deterministic held-out-fold diagnostics. The separate "
                    "leave-one-out file is descriptive only because it reuses the neutral "
                    "bank."
                ),
                "",
            ]
        )
    lines.extend(["## Gamma-SMC versus tree truth", ""])
    if comparison_paths is None:
        lines.extend(
            ["The finalized truth-versus-Gamma comparison is not available yet.", ""]
        )
    else:
        lines.extend(
            [
                _truth_gamma_section(comparison_paths["summary"], units, thresholds),
                "",
            ]
        )
    lines.extend(
        [
            "## Empirical p-values",
            "",
            (
                f"Empirical status: **{empirical_status}** ({empirical_rows} result "
                "rows). No empirical selection p-value can be reported without the "
                "user's locus-level Gamma-SMC values, sample/genotype counts, and "
                "decoder/Ne provenance. The finalized analysis writes a fail-closed "
                "input template for those data."
            ),
            "",
            "## Interpretation and limitations",
            "",
            (
                "This report deliberately does not turn the grid into a single binary "
                "'works/does not work' label. Evidence for the method is strongest when "
                "Gamma-SMC D_x is positive, selected-versus-neutral AUC is above 0.5, "
                "power is reproducible across accepted replicates, held-out neutral "
                "type-I error is near 0.05, and Gamma-SMC agrees with tree truth. The "
                "tables above retain threshold-, AF-, s-, demography-, and scope-specific "
                "heterogeneity needed for that judgment."
            ),
            "",
            "The selected results are conditional on trajectories that passed exact "
            "final-AF and genotype QC; acceptance cost must be interpreted alongside "
            "power. "
            + _han_model_contract_text(design)
            + " It is exploratory pending scaling and burn-in sensitivity analyses. "
            "The EAS "
            "model uses the PHLASH median only. There are 10 selected replicates per "
            "cell, confidence intervals are correspondingly wide, the seven x values "
            "are correlated, and no empirical data were supplied.",
            "",
            "## Machine-readable outputs",
            "",
        ]
    )
    links: list[str] = []
    if gamma_paths is not None:
        links.extend(
            [
                _relative_link(
                    report_path, gamma_paths["cell_comparisons"], "Gamma AUC table"
                ),
                _relative_link(
                    report_path,
                    gamma_paths["selected_pvalues"],
                    "selected-unit p-value table",
                ),
                _relative_link(
                    report_path, gamma_paths["conditional_power"], "Gamma power table"
                ),
                _relative_link(
                    report_path, gamma_paths["primary_omnibus"], "Gamma omnibus table"
                ),
                _relative_link(
                    report_path,
                    gamma_paths["empirical_template"],
                    "empirical input template",
                ),
            ]
        )
    if truth_paths is not None:
        links.append(
            _relative_link(
                report_path, truth_paths["cell_comparisons"], "tree-truth AUC table"
            )
        )
    if comparison_paths is not None:
        links.append(
            _relative_link(
                report_path, comparison_paths["summary"], "truth-Gamma accuracy table"
            )
        )
    if links:
        lines.append("- " + "\n- ".join(links))
    else:
        lines.append("Finalized inferential tables are not available yet.")
    lines.extend(
        [
            "",
            f"Report schema: `{REPORT_SCHEMA_VERSION}`.",
            "",
        ]
    )
    _atomic_text(report_path, "\n".join(lines))
    return report_path


__all__ = ["REPORT_SCHEMA_VERSION", "generate_run_report"]
