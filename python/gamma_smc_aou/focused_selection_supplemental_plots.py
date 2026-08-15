"""Additive AUC and focal-window figures for the focused selection study.

This module deliberately consumes the frozen, published analysis tables instead
of changing the canonical analysis implementation.  Existing simulations,
decodes, aggregates, statistical tables, figures, completion manifests, and the
campaign report therefore remain byte-for-byte reusable.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import re
import shutil
import struct
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from .eas_sweep_models import FOCAL_POSITION_BP, TMRCA_THRESHOLDS_YEARS
from .focused_selection_analysis import PRIMARY_STATISTIC
from .focused_selection_decode import sha256_file

SCHEMA_VERSION = "gamma-smc.focused-selection-supplemental-plots/v1"
ANALYSIS_SCHEMA_VERSION = "gamma-smc.focused-selection-analysis/v3"
DEFAULT_FOCAL_HALF_WINDOW_BP = 100_000
DEFAULT_THRESHOLD_YEARS = 10_000
SEQUENCE_LENGTH_BP = 10_000_000
SPATIAL_GRID_STEP_BP = 10_000
EXPECTED_DEMOGRAPHIES = (
    "eas_phlash_median",
    "ancient_eurasia_han_introgression",
)
EXPECTED_AFS = (0.1, 0.2, 0.3)
EXPECTED_SELECTION_COEFFICIENTS = (0.005, 0.01)
EXPECTED_GENOTYPE_CLASSES = ("overall", "hom_ref", "heterozygous", "hom_alt")


@dataclass(frozen=True)
class SourceSpec:
    key: str
    label: str
    spatial_source: str
    analysis_completion: str
    auc_table: str
    spatial_table: str


SOURCE_SPECS = (
    SourceSpec(
        key="gamma_smc",
        label="Gamma-SMC",
        spatial_source="gamma_smc",
        analysis_completion=(
            "focused_selection_EAS_sim/results/analysis/analysis_completion.json"
        ),
        auc_table="focused_selection_EAS_sim/results/analysis/cell_roc_auc.tsv",
        spatial_table=(
            "focused_selection_EAS_sim/results/aggregated_spatial_profiles.tsv.gz"
        ),
    ),
    SourceSpec(
        key="tree_truth",
        label="Tree truth",
        spatial_source="tree_truth",
        analysis_completion=(
            "focused_selection_EAS_sim/results/analysis_tree_truth/"
            "analysis_completion.json"
        ),
        auc_table=(
            "focused_selection_EAS_sim/results/analysis_tree_truth/cell_roc_auc.tsv"
        ),
        spatial_table=(
            "focused_selection_EAS_sim/results/aggregated_truth_spatial_profiles.tsv.gz"
        ),
    ),
)

AUC_REQUIRED_COLUMNS = {
    "demography_id",
    "target_allele_frequency",
    "statistic_scope",
    "metric",
    "statistic",
    "threshold_years",
    "selection_coefficient",
    "n_selected",
    "n_neutral",
    "roc_auc",
    "roc_auc_ci95_low",
    "roc_auc_ci95_high",
}
SPATIAL_REQUIRED_COLUMNS = {
    "source",
    "demography_id",
    "simulation_class",
    "selection_coefficient",
    "target_allele_frequency",
    "genotype_class",
    "threshold_years",
    "position_0based",
    "n_units",
    "mean_p_tmrca_lt_threshold",
    "sem_p_tmrca_lt_threshold",
}


@dataclass
class SourceBundle:
    spec: SourceSpec
    auc: pd.DataFrame
    spatial: pd.DataFrame
    input_records: dict[str, Any]
    upstream_contract_sha256: str


def _canonical_sha256(value: Any) -> str:
    payload = json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _canonical_text_sha256(path: str | Path) -> str:
    source_path = Path(path)
    try:
        raw = source_path.read_bytes()
    except OSError as error:
        raise ValueError(f"source is unreadable: {source_path}") from error
    if raw.startswith(b"\xef\xbb\xbf"):
        raise ValueError(f"source has a UTF-8 BOM: {source_path}")
    try:
        text = raw.decode("utf-8", errors="strict")
    except UnicodeDecodeError as error:
        raise ValueError(f"source is not strict UTF-8: {source_path}") from error
    if "\r" in text.replace("\r\n", ""):
        raise ValueError(f"source contains a bare carriage return: {source_path}")
    return hashlib.sha256(text.replace("\r\n", "\n").encode("utf-8")).hexdigest()


def _path_identity(path: Path, repo_root: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(repo_root.resolve()).as_posix()
    except ValueError:
        return str(resolved)


def _require_regular_file(path: Path, repo_root: Path) -> None:
    resolved_root = repo_root.resolve()
    lexical_path = Path(os.path.abspath(path))
    try:
        relative = lexical_path.relative_to(resolved_root)
    except ValueError as error:
        raise ValueError(f"required input is outside the repository: {path}") from error
    lexical = resolved_root
    for component in relative.parts:
        lexical = lexical / component
        if lexical.is_symlink():
            raise ValueError(f"required input traverses a symlink: {path}")
    resolved = lexical_path.resolve()
    try:
        resolved.relative_to(resolved_root)
    except ValueError as error:
        raise ValueError(
            f"required input resolves outside the repository: {path}"
        ) from error
    if not resolved.is_file():
        raise ValueError(f"required input is absent: {path}")


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"JSON input is unreadable: {path}") from error
    if not isinstance(value, dict):
        raise TypeError(f"JSON input is not an object: {path}")
    return value


def _atomic_json(path: Path, value: Any) -> None:
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        temporary.write_text(
            json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_figure(figure: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        figure.savefig(
            temporary,
            format=path.suffix.removeprefix("."),
            dpi=300 if path.suffix == ".png" else None,
            bbox_inches=None,
        )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _strict_int(value: Any, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{label} must be an integer")
    return value


def _validate_window(window_bp: int) -> int:
    value = _strict_int(window_bp, "focal half-window")
    if value <= 0:
        raise ValueError("focal half-window must be positive")
    if value % SPATIAL_GRID_STEP_BP:
        raise ValueError("focal half-window must align to the 10-kb spatial grid")
    lower = FOCAL_POSITION_BP - value
    upper = FOCAL_POSITION_BP + value
    if lower < 0 or upper >= SEQUENCE_LENGTH_BP:
        raise ValueError("focal half-window extends outside the decoded region")
    return value


def _validate_threshold(threshold_years: int) -> int:
    value = _strict_int(threshold_years, "localized threshold")
    if value not in TMRCA_THRESHOLDS_YEARS:
        raise ValueError(
            "localized threshold must be one of the seven decoded thresholds"
        )
    return value


def _float_set(values: pd.Series) -> tuple[float, ...]:
    try:
        numeric = pd.to_numeric(values, errors="raise").to_numpy(dtype=float)
    except (TypeError, ValueError) as error:
        raise ValueError("numeric contract field is invalid") from error
    if not np.isfinite(numeric).all():
        raise ValueError("numeric contract field is nonfinite")
    return tuple(sorted(float(value) for value in np.unique(numeric)))


def _exact_integer_array(values: pd.Series, label: str) -> np.ndarray:
    try:
        numeric = pd.to_numeric(values, errors="raise").to_numpy(dtype=float)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} is not numeric") from error
    if not np.isfinite(numeric).all() or np.any(numeric != np.floor(numeric)):
        raise ValueError(f"{label} must contain exact finite integers")
    return numeric.astype(np.int64)


def _validate_auc_table(frame: pd.DataFrame) -> pd.DataFrame:
    missing = AUC_REQUIRED_COLUMNS - set(frame.columns)
    if missing:
        raise ValueError(f"AUC table lacks columns: {sorted(missing)}")
    primary = frame[
        (frame["statistic_scope"] == "focal_nearest")
        & (frame["metric"] == "p_tmrca_lt_threshold")
        & (frame["statistic"] == PRIMARY_STATISTIC)
    ].copy()
    keys = [
        "demography_id",
        "target_allele_frequency",
        "selection_coefficient",
        "threshold_years",
    ]
    if primary.duplicated(keys).any():
        raise ValueError("AUC primary rows contain duplicate cells")
    if len(primary) != 84:
        raise ValueError(f"AUC primary row count is {len(primary)}, expected 84")
    if tuple(sorted(primary["demography_id"].unique())) != tuple(
        sorted(EXPECTED_DEMOGRAPHIES)
    ):
        raise ValueError("AUC demography coverage is incompatible")
    if _float_set(primary["target_allele_frequency"]) != EXPECTED_AFS:
        raise ValueError("AUC allele-frequency coverage is incompatible")
    if _float_set(primary["selection_coefficient"]) != EXPECTED_SELECTION_COEFFICIENTS:
        raise ValueError("AUC selection-coefficient coverage is incompatible")
    thresholds = _exact_integer_array(primary["threshold_years"], "AUC thresholds")
    if tuple(sorted(np.unique(thresholds))) != tuple(TMRCA_THRESHOLDS_YEARS):
        raise ValueError("AUC threshold coverage is incompatible")
    numeric = primary[["roc_auc", "roc_auc_ci95_low", "roc_auc_ci95_high"]].to_numpy(
        dtype=float
    )
    if not np.isfinite(numeric).all() or ((numeric < 0) | (numeric > 1)).any():
        raise ValueError("AUC estimates or confidence intervals are invalid")
    auc = primary["roc_auc"].to_numpy(dtype=float)
    low = primary["roc_auc_ci95_low"].to_numpy(dtype=float)
    high = primary["roc_auc_ci95_high"].to_numpy(dtype=float)
    if np.any(low > auc) or np.any(auc > high):
        raise ValueError("AUC estimate lies outside its bootstrap interval")
    if set(_exact_integer_array(primary["n_selected"], "AUC selected counts")) != {10}:
        raise ValueError("AUC selected sample count is incompatible")
    if set(_exact_integer_array(primary["n_neutral"], "AUC neutral counts")) != {100}:
        raise ValueError("AUC neutral sample count is incompatible")
    return primary.sort_values(keys).reset_index(drop=True)


def _validate_spatial_table(frame: pd.DataFrame, expected_source: str) -> pd.DataFrame:
    missing = SPATIAL_REQUIRED_COLUMNS - set(frame.columns)
    if missing:
        raise ValueError(f"spatial table lacks columns: {sorted(missing)}")
    if len(frame) != 504_000:
        raise ValueError(f"spatial row count is {len(frame)}, expected 504000")
    if set(frame["source"].astype(str)) != {expected_source}:
        raise ValueError("spatial source label is incompatible")
    if tuple(sorted(frame["demography_id"].unique())) != tuple(
        sorted(EXPECTED_DEMOGRAPHIES)
    ):
        raise ValueError("spatial demography coverage is incompatible")
    if _float_set(frame["target_allele_frequency"]) != EXPECTED_AFS:
        raise ValueError("spatial allele-frequency coverage is incompatible")
    if tuple(sorted(frame["genotype_class"].unique())) != tuple(
        sorted(EXPECTED_GENOTYPE_CLASSES)
    ):
        raise ValueError("spatial genotype coverage is incompatible")
    thresholds = _exact_integer_array(frame["threshold_years"], "spatial thresholds")
    if tuple(sorted(np.unique(thresholds))) != tuple(TMRCA_THRESHOLDS_YEARS):
        raise ValueError("spatial threshold coverage is incompatible")
    expected_positions = np.arange(
        0, SEQUENCE_LENGTH_BP, SPATIAL_GRID_STEP_BP, dtype=np.int64
    )
    positions = _exact_integer_array(frame["position_0based"], "spatial positions")
    observed_positions = np.sort(np.unique(positions))
    if not np.array_equal(observed_positions, expected_positions):
        raise ValueError("spatial position grid is incompatible")
    classes = set(frame["simulation_class"].astype(str))
    if classes != {"neutral", "selected"}:
        raise ValueError("spatial simulation-class coverage is incompatible")
    selected = frame[frame["simulation_class"] == "selected"]
    neutral = frame[frame["simulation_class"] == "neutral"]
    if _float_set(selected["selection_coefficient"]) != EXPECTED_SELECTION_COEFFICIENTS:
        raise ValueError("spatial selected coefficients are incompatible")
    if _float_set(neutral["selection_coefficient"]) != (0.0,):
        raise ValueError("spatial neutral coefficient is incompatible")
    if set(_exact_integer_array(selected["n_units"], "spatial selected counts")) != {
        10
    }:
        raise ValueError("spatial selected unit count is incompatible")
    if set(_exact_integer_array(neutral["n_units"], "spatial neutral counts")) != {100}:
        raise ValueError("spatial neutral unit count is incompatible")
    key_columns = [
        "source",
        "demography_id",
        "simulation_class",
        "selection_coefficient",
        "target_allele_frequency",
        "genotype_class",
        "threshold_years",
        "position_0based",
    ]
    if frame.duplicated(key_columns).any():
        raise ValueError("spatial table contains duplicate grid rows")
    probability = frame["mean_p_tmrca_lt_threshold"].to_numpy(dtype=float)
    sem = frame["sem_p_tmrca_lt_threshold"].to_numpy(dtype=float)
    if not np.isfinite(probability).all() or np.any(
        (probability < 0) | (probability > 1)
    ):
        raise ValueError("spatial probabilities are invalid")
    if not np.isfinite(sem).all() or np.any(sem < 0):
        raise ValueError("spatial standard errors are invalid")
    return frame


def _validate_upstream_sources(manifest: dict[str, Any], repo_root: Path) -> None:
    implementation = manifest.get("contract", {}).get("implementation", {})
    sources = implementation.get("sources")
    if not isinstance(sources, dict) or not sources:
        raise ValueError("analysis manifest lacks implementation sources")
    for record in sources.values():
        if not isinstance(record, dict):
            raise TypeError("analysis source record is malformed")
        relative = record.get("path")
        expected_sha = record.get("sha256")
        if not isinstance(relative, str) or not isinstance(expected_sha, str):
            raise TypeError("analysis source record is malformed")
        source_path = repo_root / relative
        _require_regular_file(source_path, repo_root)
        if _canonical_text_sha256(source_path) != expected_sha:
            raise ValueError(f"analysis source hash is incompatible: {relative}")


def _load_source_bundle(repo_root: Path, spec: SourceSpec) -> SourceBundle:
    completion_path = repo_root / spec.analysis_completion
    auc_path = repo_root / spec.auc_table
    spatial_path = repo_root / spec.spatial_table
    for path in (completion_path, auc_path, spatial_path):
        _require_regular_file(path, repo_root)
    manifest = _read_json(completion_path)
    if (
        manifest.get("schema") != ANALYSIS_SCHEMA_VERSION
        or manifest.get("status") != "complete"
    ):
        raise ValueError(f"{spec.key} analysis completion is incompatible")
    contract = manifest.get("contract")
    if not isinstance(contract, dict):
        raise TypeError(f"{spec.key} analysis contract is absent")
    contract_sha = _canonical_sha256(contract)
    if manifest.get("contract_sha256") != contract_sha:
        raise ValueError(f"{spec.key} analysis contract hash is incompatible")
    _validate_upstream_sources(manifest, repo_root)
    auc_record = manifest.get("outputs", {}).get("cell_comparisons")
    if not isinstance(auc_record, dict) or auc_record.get("path") != auc_path.name:
        raise ValueError(f"{spec.key} AUC output record is incompatible")
    auc_sha = sha256_file(auc_path)
    if (
        auc_record.get("sha256") != auc_sha
        or auc_record.get("size_bytes") != auc_path.stat().st_size
    ):
        raise ValueError(f"{spec.key} AUC table hash or size is incompatible")
    inputs = contract.get("inputs")
    if not isinstance(inputs, dict):
        raise TypeError(f"{spec.key} analysis input record is absent")
    spatial_identity = _path_identity(spatial_path, repo_root)
    spatial_sha = sha256_file(spatial_path)
    if (
        inputs.get("spatial_summaries_path") != spatial_identity
        or inputs.get("spatial_summaries_sha256") != spatial_sha
    ):
        raise ValueError(f"{spec.key} spatial input binding is incompatible")
    try:
        auc = pd.read_csv(auc_path, sep="\t")
        spatial = pd.read_csv(spatial_path, sep="\t")
    except Exception as error:
        raise ValueError(
            f"{spec.key} supplemental plotting input is unreadable"
        ) from error
    primary_auc = _validate_auc_table(auc)
    valid_spatial = _validate_spatial_table(spatial, spec.spatial_source)
    return SourceBundle(
        spec=spec,
        auc=primary_auc,
        spatial=valid_spatial,
        upstream_contract_sha256=contract_sha,
        input_records={
            "analysis_completion": {
                "path": spec.analysis_completion,
                "sha256": sha256_file(completion_path),
                "size_bytes": completion_path.stat().st_size,
                "contract_sha256": contract_sha,
            },
            "auc_table": {
                "path": spec.auc_table,
                "sha256": auc_sha,
                "size_bytes": auc_path.stat().st_size,
            },
            "spatial_table": {
                "path": spec.spatial_table,
                "sha256": spatial_sha,
                "size_bytes": spatial_path.stat().st_size,
            },
        },
    )


def _implementation_contract(repo_root: Path) -> dict[str, Any]:
    source_paths = (
        Path(__file__).resolve(),
        repo_root / "scripts/run_focused_selection_supplemental_plots.py",
    )
    sources: dict[str, dict[str, str]] = {}
    for path in source_paths:
        _require_regular_file(path, repo_root)
        sources[path.name] = {
            "path": _path_identity(path, repo_root),
            "sha256": _canonical_text_sha256(path),
        }
    return {
        "sources": sources,
        "source_hashing": {
            "encoding": "strict UTF-8",
            "newline_canonicalization": "CRLF to LF",
            "unicode_normalization": "none",
            "reject_bom": True,
            "reject_bare_carriage_return": True,
        },
        "path_identity": "repo-relative POSIX for in-repository paths",
        "software": {
            "python": platform.python_version(),
            "matplotlib": matplotlib.__version__,
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
    }


def _build_contract(
    repo_root: Path,
    bundles: list[SourceBundle],
    *,
    focal_half_window_bp: int,
    threshold_years: int,
) -> dict[str, Any]:
    return {
        "schema": SCHEMA_VERSION,
        "implementation": _implementation_contract(repo_root),
        "inputs": {bundle.spec.key: bundle.input_records for bundle in bundles},
        "parameters": {
            "auc_scope": "focal_nearest",
            "auc_metric": "p_tmrca_lt_threshold",
            "auc_statistic": PRIMARY_STATISTIC,
            "auc_intervals": "existing pointwise bootstrap 95% intervals",
            "focal_position_0based": FOCAL_POSITION_BP,
            "focal_half_window_bp": focal_half_window_bp,
            "localized_threshold_years": threshold_years,
            "spatial_grid_step_bp": SPATIAL_GRID_STEP_BP,
            "spatial_interval": "mean +/- 1.96 * SEM, clipped to [0, 1]",
            "randomness": "none; published estimates and intervals are reused",
            "figure_size_inches": [11, 8.5],
            "png_dpi": 300,
            "pdf_bbox_inches": None,
        },
    }


def _save_pair(figure: plt.Figure, directory: Path, stem: str) -> list[Path]:
    outputs: list[Path] = []
    for extension in ("png", "pdf"):
        path = directory / f"{stem}.{extension}"
        _atomic_figure(figure, path)
        outputs.append(path)
    return outputs


def plot_auc_profiles(bundle: SourceBundle, output_dir: Path) -> list[Path]:
    """Plot all seven AUC thresholds with pointwise bootstrap intervals."""

    figure, axes = plt.subplots(
        2, 3, figsize=(11, 8.5), sharex=True, sharey=True, squeeze=False
    )
    colors = {0.005: "#E69F00", 0.01: "#0072B2"}
    markers = {0.005: "o", 0.01: "s"}
    thresholds = np.asarray(TMRCA_THRESHOLDS_YEARS, dtype=float) / 1_000
    for row, demography in enumerate(EXPECTED_DEMOGRAPHIES):
        for column, af in enumerate(EXPECTED_AFS):
            axis = axes[row, column]
            cell = bundle.auc[
                (bundle.auc["demography_id"] == demography)
                & np.isclose(bundle.auc["target_allele_frequency"], af)
            ]
            for coefficient in EXPECTED_SELECTION_COEFFICIENTS:
                curve = cell[
                    np.isclose(cell["selection_coefficient"], coefficient)
                ].sort_values("threshold_years")
                if len(curve) != len(TMRCA_THRESHOLDS_YEARS):
                    raise ValueError("AUC plot lacks a complete threshold curve")
                y = curve["roc_auc"].to_numpy(dtype=float)
                low = curve["roc_auc_ci95_low"].to_numpy(dtype=float)
                high = curve["roc_auc_ci95_high"].to_numpy(dtype=float)
                axis.plot(
                    thresholds,
                    y,
                    color=colors[coefficient],
                    marker=markers[coefficient],
                    markersize=5.5,
                    linewidth=2.5,
                )
                axis.fill_between(
                    thresholds,
                    low,
                    high,
                    color=colors[coefficient],
                    alpha=0.17,
                )
            axis.axhline(0.5, color="#666666", linestyle="--", linewidth=1.6)
            axis.set_ylim(0, 1.03)
            axis.set_xticks(
                thresholds,
                [f"{value:g}" for value in thresholds],
                rotation=45,
                ha="right",
                fontsize=10.5,
            )
            axis.tick_params(axis="y", labelsize=12)
            axis.grid(alpha=0.22)
            axis.set_title(f"Final AF {100 * af:.0f}%", fontsize=15)
        demography_label = (
            "EAS median" if demography == "eas_phlash_median" else "Han introgression"
        )
        axes[row, 0].set_ylabel(f"{demography_label}\nROC AUC", fontsize=14)
    figure.suptitle(
        f"{bundle.spec.label} selected-vs-neutral ROC AUC\n"
        "Focal alt/alt-minus-ref/ref score",
        fontsize=20,
        y=0.98,
    )
    legend_handles = [
        Line2D(
            [0],
            [0],
            color=colors[coefficient],
            marker=markers[coefficient],
            linewidth=2.5,
            label=f"s={coefficient:g}",
        )
        for coefficient in EXPECTED_SELECTION_COEFFICIENTS
    ]
    legend_handles.append(
        Line2D(
            [0],
            [0],
            color="#666666",
            linestyle="--",
            linewidth=1.6,
            label="AUC=0.5",
        )
    )
    figure.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.885),
        ncol=3,
        fontsize=14,
        frameon=True,
    )
    figure.supxlabel("TMRCA threshold (kya)", fontsize=16, y=0.055)
    figure.subplots_adjust(
        left=0.11, right=0.98, bottom=0.16, top=0.76, wspace=0.16, hspace=0.32
    )
    paths = _save_pair(figure, output_dir, "primary_auc_profiles")
    plt.close(figure)
    return paths


def _window_slug(window_bp: int) -> str:
    if window_bp % 1_000 == 0:
        return f"{window_bp // 1_000}kb"
    return f"{window_bp}bp"


def _demography_title(demography: str) -> str:
    return (
        "EAS median demography"
        if demography == "eas_phlash_median"
        else "Han AncientEurasia introgression"
    )


def plot_localized_profiles(
    bundle: SourceBundle,
    output_dir: Path,
    *,
    focal_half_window_bp: int,
    threshold_years: int,
) -> list[Path]:
    """Plot selected and matched-neutral profiles around the focal site."""

    window_bp = _validate_window(focal_half_window_bp)
    threshold = _validate_threshold(threshold_years)
    lower = FOCAL_POSITION_BP - window_bp
    upper = FOCAL_POSITION_BP + window_bp
    local = bundle.spatial[
        (bundle.spatial["position_0based"] >= lower)
        & (bundle.spatial["position_0based"] <= upper)
        & (bundle.spatial["threshold_years"] == threshold)
    ]
    expected_positions = np.arange(
        lower, upper + SPATIAL_GRID_STEP_BP, SPATIAL_GRID_STEP_BP
    )
    observed_positions = np.sort(
        np.unique(_exact_integer_array(local["position_0based"], "localized positions"))
    )
    if not np.array_equal(observed_positions, expected_positions):
        raise ValueError("localized spatial window has incomplete grid coverage")
    colors = {
        "overall": "#000000",
        "hom_ref": "#595959",
        "heterozygous": "#E69F00",
        "hom_alt": "#0072B2",
    }
    outputs: list[Path] = []
    window_kb = window_bp / 1_000
    ticks = np.linspace(-window_kb, window_kb, 5)
    for demography in EXPECTED_DEMOGRAPHIES:
        demographic = local[local["demography_id"] == demography]
        upper_band = np.clip(
            demographic["mean_p_tmrca_lt_threshold"].to_numpy(dtype=float)
            + 1.96 * demographic["sem_p_tmrca_lt_threshold"].to_numpy(dtype=float),
            0,
            1,
        )
        maximum_upper_band = float(np.max(upper_band))
        y_maximum = min(1.0, maximum_upper_band + max(0.03, 0.06 * maximum_upper_band))
        figure, axes = plt.subplots(
            2, 3, figsize=(11, 8.5), sharex=True, sharey=True, squeeze=False
        )
        for row, coefficient in enumerate(EXPECTED_SELECTION_COEFFICIENTS):
            comparison = demographic[
                (demographic["simulation_class"] == "neutral")
                | (
                    (demographic["simulation_class"] == "selected")
                    & np.isclose(demographic["selection_coefficient"], coefficient)
                )
            ]
            for column, af in enumerate(EXPECTED_AFS):
                axis = axes[row, column]
                cell = comparison[np.isclose(comparison["target_allele_frequency"], af)]
                for simulation_class, linestyle in (
                    ("selected", "-"),
                    ("neutral", "--"),
                ):
                    for genotype_class in EXPECTED_GENOTYPE_CLASSES:
                        curve = cell[
                            (cell["simulation_class"] == simulation_class)
                            & (cell["genotype_class"] == genotype_class)
                        ].sort_values("position_0based")
                        if len(curve) != len(expected_positions):
                            raise ValueError(
                                "localized plot lacks a selected/neutral genotype curve"
                            )
                        x = (
                            curve["position_0based"].to_numpy(dtype=float)
                            - FOCAL_POSITION_BP
                        ) / 1_000
                        y = curve["mean_p_tmrca_lt_threshold"].to_numpy(dtype=float)
                        sem = curve["sem_p_tmrca_lt_threshold"].to_numpy(dtype=float)
                        color = colors[genotype_class]
                        axis.plot(
                            x,
                            y,
                            color=color,
                            linestyle=linestyle,
                            linewidth=2.1,
                            label=(
                                f"{simulation_class.capitalize()} "
                                f"{genotype_class.replace('_', '/')}"
                            ),
                        )
                        axis.fill_between(
                            x,
                            np.clip(y - 1.96 * sem, 0, 1),
                            np.clip(y + 1.96 * sem, 0, 1),
                            color=color,
                            alpha=0.07 if simulation_class == "selected" else 0.035,
                        )
                axis.axvline(0, color="#CC3311", linewidth=1.6, alpha=0.9)
                axis.set_xlim(-window_kb, window_kb)
                axis.set_xticks(
                    ticks,
                    [f"{value:g}" for value in ticks],
                    fontsize=11,
                )
                axis.tick_params(axis="y", labelsize=11)
                axis.grid(alpha=0.20)
                if row == 0:
                    axis.set_title(f"Final AF {100 * af:.0f}%", fontsize=15)
            axes[row, 0].set_ylabel(
                f"s={coefficient:g}\nMean P(TMRCA < {threshold / 1_000:g} kya)",
                fontsize=13,
            )
        axes[0, 0].set_ylim(0, y_maximum)
        handles, labels = axes[-1, -1].get_legend_handles_labels()
        figure.suptitle(
            f"{bundle.spec.label} focal-window genomic profiles\n"
            f"{_demography_title(demography)}; +/-{window_kb:g} kb",
            fontsize=19,
            y=0.985,
        )
        figure.legend(
            handles,
            labels,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.875),
            ncol=4,
            fontsize=11.5,
            frameon=True,
        )
        figure.supxlabel("Distance from focal site (kb)", fontsize=16, y=0.055)
        figure.subplots_adjust(
            left=0.105,
            right=0.97,
            bottom=0.13,
            top=0.70,
            wspace=0.16,
            hspace=0.28,
        )
        stem = (
            f"localized_spatial_{demography}_t{threshold}_"
            f"window{_window_slug(window_bp)}"
        )
        outputs.extend(_save_pair(figure, output_dir, stem))
        plt.close(figure)
    return outputs


def _png_dimensions(path: Path) -> tuple[int, int]:
    with path.open("rb") as stream:
        signature = stream.read(24)
    if len(signature) != 24 or signature[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError(f"PNG output is invalid: {path}")
    return struct.unpack(">II", signature[16:24])


def _pdf_geometry(path: Path) -> tuple[int, float, float]:
    try:
        payload = path.read_bytes()
    except OSError as error:
        raise ValueError(f"PDF output is unreadable: {path}") from error
    if not payload.startswith(b"%PDF-"):
        raise ValueError(f"PDF output is invalid: {path}")
    page_count = len(re.findall(rb"/Type\s*/Page(?!s)\b", payload))
    media_boxes = re.findall(
        rb"/MediaBox\s*\[\s*([-+0-9.]+)\s+([-+0-9.]+)\s+"
        rb"([-+0-9.]+)\s+([-+0-9.]+)\s*\]",
        payload,
    )
    if page_count != 1 or len(media_boxes) != 1:
        raise ValueError(f"PDF output page contract is incompatible: {path}")
    try:
        x0, y0, x1, y1 = (float(value) for value in media_boxes[0])
    except ValueError as error:
        raise ValueError(f"PDF output media box is invalid: {path}") from error
    width = x1 - x0
    height = y1 - y0
    if not np.isclose(width, 792) or not np.isclose(height, 612):
        raise ValueError(f"PDF output is not exact landscape letter size: {path}")
    return page_count, width, height


def _output_record(path: Path, root: Path) -> dict[str, Any]:
    record: dict[str, Any] = {
        "path": path.relative_to(root).as_posix(),
        "sha256": sha256_file(path),
        "size_bytes": path.stat().st_size,
    }
    if path.suffix == ".png":
        width, height = _png_dimensions(path)
        record["width_px"] = width
        record["height_px"] = height
    elif path.suffix == ".pdf":
        page_count, width, height = _pdf_geometry(path)
        record["page_count"] = page_count
        record["width_points"] = width
        record["height_points"] = height
    return record


def _readme_text(contract: dict[str, Any]) -> str:
    parameters = contract["parameters"]
    return (
        "# Focused selection supplemental AUC and focal-window plots\n\n"
        "These figures are additive views of the frozen published analysis. No "
        "simulation, decoding, aggregation, statistical testing, or canonical "
        "analysis output was rerun.\n\n"
        "- AUC figures show the prespecified focal alt/alt-minus-ref/ref score "
        "for selected versus matched-neutral replicates across all seven TMRCA "
        "thresholds. Bands are the already-published pointwise bootstrap 95% "
        "intervals. AUC below 0.5 is retained as a direction reversal.\n"
        f"- Localized figures show +/-{parameters['focal_half_window_bp'] // 1000} "
        "kb around the focal site on the existing 10-kb grid at "
        f"{parameters['localized_threshold_years'] / 1000:g} kya. Bands are "
        "mean +/- 1.96 SEM and do not represent a covariance-aware contrast "
        "interval.\n"
        "- Selected cells contain 10 replicates and matched neutral cells contain "
        "100 replicates. Thresholds and selection-coefficient comparisons reuse "
        "correlated inputs.\n"
    )


def _validate_cached_output(
    output_dir: Path, completion: dict[str, Any], contract: dict[str, Any]
) -> dict[str, Any]:
    if (
        completion.get("schema") != SCHEMA_VERSION
        or completion.get("status") != "complete"
    ):
        raise ValueError("supplemental completion is incompatible")
    if completion.get("contract") != contract:
        raise ValueError("supplemental plot contract is incompatible")
    contract_sha = _canonical_sha256(contract)
    if completion.get("contract_sha256") != contract_sha:
        raise ValueError("supplemental contract hash is incompatible")
    outputs = completion.get("outputs")
    if not isinstance(outputs, dict) or not outputs:
        raise ValueError("supplemental completion lacks outputs")
    parameters = contract["parameters"]
    window_slug = _window_slug(parameters["focal_half_window_bp"])
    threshold = parameters["localized_threshold_years"]
    expected_paths = {"README.md"}
    for spec in SOURCE_SPECS:
        for extension in ("png", "pdf"):
            expected_paths.add(f"{spec.key}/primary_auc_profiles.{extension}")
            for demography in EXPECTED_DEMOGRAPHIES:
                expected_paths.add(
                    f"{spec.key}/localized_spatial_{demography}_t{threshold}_"
                    f"window{window_slug}.{extension}"
                )
    if set(outputs) != expected_paths:
        raise ValueError("supplemental completion output set is incompatible")
    completion_path = output_dir / "supplemental_plots_completion.json"
    if completion_path.is_symlink() or not completion_path.is_file():
        raise ValueError("supplemental completion path is incompatible")
    allowed_directories = {spec.key for spec in SOURCE_SPECS}
    for path in output_dir.rglob("*"):
        relative = path.relative_to(output_dir).as_posix()
        if path.is_symlink():
            raise ValueError(f"supplemental output traverses a symlink: {relative}")
        if path.is_dir() and relative not in allowed_directories:
            raise ValueError(f"supplemental output directory is unexpected: {relative}")
    actual_paths = {
        path.relative_to(output_dir).as_posix()
        for path in output_dir.rglob("*")
        if path.is_file() and path.name != "supplemental_plots_completion.json"
    }
    if actual_paths != expected_paths:
        raise ValueError("supplemental output disk set is incompatible")
    for relative, record in outputs.items():
        if not isinstance(record, dict) or record.get("path") != relative:
            raise ValueError(f"supplemental output record is incompatible: {relative}")
        path = output_dir / relative
        if path.is_symlink() or not path.is_file():
            raise ValueError(f"supplemental output is absent: {relative}")
        if record.get("size_bytes") != path.stat().st_size:
            raise ValueError(f"supplemental output size is incompatible: {relative}")
        if record.get("sha256") != sha256_file(path):
            raise ValueError(f"supplemental output hash is incompatible: {relative}")
        if path.suffix == ".png":
            width, height = _png_dimensions(path)
            if (width, height) != (3_300, 2_550):
                raise ValueError(
                    f"supplemental PNG geometry is incompatible: {relative}"
                )
            if record.get("width_px") != width or record.get("height_px") != height:
                raise ValueError(f"supplemental PNG record is incompatible: {relative}")
        elif path.suffix == ".pdf":
            page_count, width, height = _pdf_geometry(path)
            if (
                record.get("page_count") != page_count
                or record.get("width_points") != width
                or record.get("height_points") != height
            ):
                raise ValueError(f"supplemental PDF record is incompatible: {relative}")
    return {
        "status": "verified_cached",
        "output_dir": str(output_dir),
        "contract_sha256": contract_sha,
        "output_count": len(outputs),
        "png_count": sum(path.endswith(".png") for path in outputs),
        "pdf_count": sum(path.endswith(".pdf") for path in outputs),
    }


def generate_supplemental_plots(
    repo_root: str | Path,
    *,
    output_dir: str | Path | None = None,
    focal_half_window_bp: int = DEFAULT_FOCAL_HALF_WINDOW_BP,
    threshold_years: int = DEFAULT_THRESHOLD_YEARS,
    verify_only: bool = False,
) -> dict[str, Any]:
    """Generate or verify the additive supplemental figure bundle."""

    root = Path(repo_root).resolve()
    if not root.is_dir():
        raise ValueError(f"repository root is absent: {root}")
    window_bp = _validate_window(focal_half_window_bp)
    threshold = _validate_threshold(threshold_years)
    bundles = [_load_source_bundle(root, spec) for spec in SOURCE_SPECS]
    contract = _build_contract(
        root,
        bundles,
        focal_half_window_bp=window_bp,
        threshold_years=threshold,
    )
    destination = (
        Path(output_dir)
        if output_dir is not None
        else root / "focused_selection_EAS_sim/results/supplemental_auc_localized"
    )
    if not destination.is_absolute():
        destination = root / destination
    destination = destination.resolve()
    completion_path = destination / "supplemental_plots_completion.json"
    if completion_path.is_file():
        return _validate_cached_output(
            destination, _read_json(completion_path), contract
        )
    if verify_only:
        raise ValueError("supplemental completion is absent")
    if destination.exists() and any(destination.iterdir()):
        raise ValueError(
            "supplemental output directory is nonempty without completion; "
            "preserve and inspect the partial directory before retrying"
        )
    temporary = destination.with_name(destination.name + f".tmp.{os.getpid()}")
    if temporary.exists():
        raise ValueError(
            f"supplemental temporary directory already exists: {temporary}"
        )
    temporary.parent.mkdir(parents=True, exist_ok=True)
    temporary.mkdir()
    try:
        paths: list[Path] = []
        for bundle in bundles:
            source_dir = temporary / bundle.spec.key
            paths.extend(plot_auc_profiles(bundle, source_dir))
            paths.extend(
                plot_localized_profiles(
                    bundle,
                    source_dir,
                    focal_half_window_bp=window_bp,
                    threshold_years=threshold,
                )
            )
        readme = temporary / "README.md"
        readme.write_text(_readme_text(contract), encoding="utf-8")
        paths.append(readme)
        if len([path for path in paths if path.suffix == ".png"]) != 6:
            raise ValueError("supplemental PNG output count is incompatible")
        if len([path for path in paths if path.suffix == ".pdf"]) != 6:
            raise ValueError("supplemental PDF output count is incompatible")
        output_records = {
            path.relative_to(temporary).as_posix(): _output_record(path, temporary)
            for path in sorted(paths)
        }
        completion = {
            "schema": SCHEMA_VERSION,
            "status": "complete",
            "contract": contract,
            "contract_sha256": _canonical_sha256(contract),
            "outputs": output_records,
            "interpretation": {
                "auc": (
                    "selected-versus-matched-neutral discrimination of the focal "
                    "alt/alt-minus-ref/ref score; not empirical prediction"
                ),
                "localized_profiles": (
                    "published aggregate means and marginal SEMs on the existing "
                    "10-kb grid; not position-specific AUC"
                ),
            },
        }
        _atomic_json(temporary / "supplemental_plots_completion.json", completion)
        if destination.exists():
            destination.rmdir()
        os.replace(temporary, destination)
    except Exception:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    result = _validate_cached_output(destination, _read_json(completion_path), contract)
    result["status"] = "generated"
    return result


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate or verify focused-study AUC and focal-window plots."
    )
    parser.add_argument("action", choices=("generate", "verify"))
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument(
        "--focal-half-window-bp",
        type=int,
        default=DEFAULT_FOCAL_HALF_WINDOW_BP,
    )
    parser.add_argument("--threshold-years", type=int, default=DEFAULT_THRESHOLD_YEARS)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    result = generate_supplemental_plots(
        args.repo_root,
        output_dir=args.output_dir,
        focal_half_window_bp=args.focal_half_window_bp,
        threshold_years=args.threshold_years,
        verify_only=args.action == "verify",
    )
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


__all__ = [
    "DEFAULT_FOCAL_HALF_WINDOW_BP",
    "DEFAULT_THRESHOLD_YEARS",
    "SCHEMA_VERSION",
    "generate_supplemental_plots",
    "main",
    "plot_auc_profiles",
    "plot_localized_profiles",
]
