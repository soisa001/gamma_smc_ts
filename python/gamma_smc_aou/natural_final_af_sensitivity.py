"""Additive final-AF sensitivity study for weaker selection coefficients.

This module deliberately leaves :mod:`gamma_smc_aou.natural_final_af` and its
published ``s=0``/``s=0.01`` result bundle unchanged.  It simulates only the
four new demography-by-selection cells for ``s=0.001`` and ``s=0.005``.  The
analysis then combines those trajectories with the checksum-validated base
summary and population-survivor table.  Population survival
(``final_alt_count > 0``) remains the only analysis filter.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import tempfile
from collections.abc import Callable, Mapping, Sequence
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from . import natural_final_af as base
from .eas_sweep_models import GENERATION_TIME_YEARS
from .eas_sweep_study import BASE_SEED


SCHEMA_VERSION = "gamma-smc.natural-final-af-selection-sensitivity/v1"
ANALYSIS_SCHEMA = "gamma-smc.natural-final-af-selection-sensitivity-analysis/v1"
MODULE_PATH = "python/gamma_smc_aou/natural_final_af_sensitivity.py"
WRAPPER_PATH = "scripts/run_natural_final_af_sensitivity.py"
DEFAULT_CAMPAIGN_DIR = base.DEFAULT_CAMPAIGN_DIR
DEFAULT_WORK_SUBDIR = "work/natural_final_af_selection_sensitivity"
DEFAULT_RESULTS_SUBDIR = "results/natural_final_af_selection_sensitivity"
BASE_RESULTS_SUBDIR = base.DEFAULT_RESULTS_SUBDIR

NEW_SELECTION_COEFFICIENTS = (0.001, 0.005)
COMBINED_SELECTION_COEFFICIENTS = (0.0, 0.001, 0.005, 0.01)

# Fixed attempted trajectories; the run never stops on a survivor target.
DEFAULT_EAS_S001_ATTEMPTS = 1_500_000
DEFAULT_EAS_S005_ATTEMPTS = 500_000
DEFAULT_HAN_S001_ATTEMPTS = 20_000
DEFAULT_HAN_S005_ATTEMPTS = 20_000
DEFAULT_EAS_BATCH_SIZE = 50_000
DEFAULT_HAN_BATCH_SIZE = 1_000

ATTEMPT_COLUMNS = base.ATTEMPT_COLUMNS
SUMMARY_COLUMNS = base.SUMMARY_COLUMNS
COMPARISON_COLUMNS = (
    "demography_id",
    "neutral_selection_coefficient",
    "selected_selection_coefficient",
    "statistic",
    "population_af_roc_auc_survivors",
    "n_neutral_survivors",
    "n_selected_survivors",
)

ANALYSIS_OUTPUT_NAMES = {
    "attempts": "new_selection_attempts.tsv.gz",
    "survivors": "final_af_population_survivors.tsv.gz",
    "summary": "final_af_summary.tsv",
    "comparisons": "final_af_comparisons.tsv",
    "readme": "README.md",
    "fate_png": "unconditional_fate_fractions.png",
    "fate_pdf": "unconditional_fate_fractions.pdf",
    "population_ecdf_png": "final_population_af_survivor_ecdf.png",
    "population_ecdf_pdf": "final_population_af_survivor_ecdf.pdf",
    "population_sample_ecdf_png": "population_vs_sample_af_survivor_ecdf.png",
    "population_sample_ecdf_pdf": "population_vs_sample_af_survivor_ecdf.pdf",
}

TABLE_COLUMNS = {
    "attempts": ATTEMPT_COLUMNS,
    "survivors": ATTEMPT_COLUMNS,
    "summary": SUMMARY_COLUMNS,
    "comparisons": COMPARISON_COLUMNS,
}

DEMOGRAPHY_LABELS = {
    "eas_phlash_median": "EAS age-matched de novo (50 kya)",
    "ancient_eurasia_han_introgression": "Han donor-fixed pulse (~2.96%)",
}
COEFFICIENT_COLORS = {
    0.0: "#4C78A8",
    0.001: "#59A14F",
    0.005: "#F28E2B",
    0.01: "#D1495B",
}


@dataclass(frozen=True)
class SensitivityPlan:
    """Immutable run configuration for the additive sensitivity extension."""

    repo_root: Path
    campaign_dir: Path
    eas_s001_attempts: int = DEFAULT_EAS_S001_ATTEMPTS
    eas_s005_attempts: int = DEFAULT_EAS_S005_ATTEMPTS
    han_s001_attempts: int = DEFAULT_HAN_S001_ATTEMPTS
    han_s005_attempts: int = DEFAULT_HAN_S005_ATTEMPTS
    eas_batch_size: int = DEFAULT_EAS_BATCH_SIZE
    han_batch_size: int = DEFAULT_HAN_BATCH_SIZE
    base_seed: int = BASE_SEED

    def validate(self) -> None:
        root = self.repo_root.resolve()
        campaign = self.campaign_dir.resolve()
        try:
            campaign.relative_to(root)
        except ValueError as error:
            raise ValueError("campaign_dir must be inside repo_root") from error
        if any("onedrive" in part.casefold() for part in campaign.parts):
            raise ValueError("sensitivity work must not be inside OneDrive")
        for field in (
            "eas_s001_attempts",
            "eas_s005_attempts",
            "han_s001_attempts",
            "han_s005_attempts",
            "eas_batch_size",
            "han_batch_size",
            "base_seed",
        ):
            value = getattr(self, field)
            if isinstance(value, bool) or int(value) != value or int(value) < 1:
                raise ValueError(f"{field} must be a positive integer")


def _coefficient_label(value: float) -> str:
    if math.isclose(float(value), 0.0, rel_tol=0, abs_tol=1e-15):
        return "s=0"
    return f"s={float(value):.3f}".rstrip("0")


def _coefficient_key(value: float) -> int:
    """Return an exact integer key at the study's 1e-3 coefficient resolution."""

    scaled = float(value) * 1_000
    rounded = int(round(scaled))
    if not math.isclose(scaled, rounded, rel_tol=0, abs_tol=1e-12):
        raise ValueError(f"unsupported selection coefficient: {value}")
    return rounded


def _cell_id(demography: str, coefficient: float) -> str:
    return f"{demography}__selected_s{coefficient:.3f}".replace(".", "p")


def build_cells(plan: SensitivityPlan) -> list[dict[str, Any]]:
    """Build exactly the four new positive-selection cells."""

    plan.validate()
    specifications = (
        (
            "eas_phlash_median",
            "EAS",
            0.001,
            plan.eas_s001_attempts,
            plan.eas_batch_size,
        ),
        (
            "eas_phlash_median",
            "EAS",
            0.005,
            plan.eas_s005_attempts,
            plan.eas_batch_size,
        ),
        (
            "ancient_eurasia_han_introgression",
            "Han",
            0.001,
            plan.han_s001_attempts,
            plan.han_batch_size,
        ),
        (
            "ancient_eurasia_han_introgression",
            "Han",
            0.005,
            plan.han_s005_attempts,
            plan.han_batch_size,
        ),
    )
    cells = []
    for demography, population, coefficient, attempts, batch_size in specifications:
        cells.append(
            {
                "cell_id": _cell_id(demography, coefficient),
                "demography_id": demography,
                "population": population,
                "simulation_class": "selected",
                "selection_coefficient": coefficient,
                "dominance_coefficient": base.DOMINANCE_COEFFICIENT,
                "attempts": int(attempts),
                "batch_size": int(batch_size),
                "n_batches": int(math.ceil(attempts / batch_size)),
            }
        )
    if {_coefficient_key(cell["selection_coefficient"]) for cell in cells} != {1, 5}:
        raise AssertionError("sensitivity plan must contain only s=0.001 and s=0.005")
    return cells


def _base_result_dir(plan: SensitivityPlan) -> Path:
    return plan.campaign_dir.resolve() / BASE_RESULTS_SUBDIR


def _base_cells_expected(plan: SensitivityPlan) -> list[dict[str, Any]]:
    base_plan = base.NaturalAfPlan(
        repo_root=plan.repo_root.resolve(),
        campaign_dir=plan.campaign_dir.resolve(),
        base_seed=plan.base_seed,
    )
    return base.build_cells(base_plan)


def _base_bundle_contract(plan: SensitivityPlan) -> dict[str, Any]:
    """Validate and describe the immutable published s=0/s=0.01 input bundle."""

    plan.validate()
    root = plan.repo_root.resolve()
    result_dir = _base_result_dir(plan)
    completion_path = result_dir / "analysis_completion.json"
    if (
        not result_dir.is_dir()
        or result_dir.is_symlink()
        or not completion_path.is_file()
        or completion_path.is_symlink()
    ):
        raise ValueError("published natural-final-AF base bundle is unavailable")
    completion_raw = json.loads(completion_path.read_text(encoding="utf-8"))
    recorded_contract = completion_raw.get("contract")
    if not isinstance(recorded_contract, Mapping):
        raise ValueError("base analysis completion has no valid contract")
    completion = base._validate_analysis_bundle(result_dir, recorded_contract)
    if completion.get("schema") != base.ANALYSIS_SCHEMA:
        raise ValueError("base analysis schema is incompatible")

    current_implementation = base._implementation_contract(root)
    if base._canonical_sha256(recorded_contract.get("implementation")) != (
        base._canonical_sha256(current_implementation)
    ):
        raise ValueError("base implementation no longer matches its published bundle")
    expected_cells = _base_cells_expected(plan)
    if base._canonical_sha256(recorded_contract.get("cells")) != (
        base._canonical_sha256(expected_cells)
    ):
        raise ValueError("base bundle cells do not match the frozen production design")

    summary_path = result_dir / base.ANALYSIS_OUTPUT_NAMES["summary"]
    survivors_path = result_dir / base.ANALYSIS_OUTPUT_NAMES["survivors"]
    summary = pd.read_csv(summary_path, sep="\t")
    survivors = pd.read_csv(survivors_path, sep="\t")
    if tuple(summary.columns) != SUMMARY_COLUMNS or tuple(survivors.columns) != (
        ATTEMPT_COLUMNS
    ):
        raise ValueError("base summary/survivor columns are incompatible")
    expected_keys = {
        (
            str(cell["demography_id"]),
            _coefficient_key(float(cell["selection_coefficient"])),
        )
        for cell in expected_cells
    }
    summary_keys = {
        (str(row.demography_id), _coefficient_key(float(row.selection_coefficient)))
        for row in summary.itertuples(index=False)
    }
    if len(summary) != 4 or summary_keys != expected_keys:
        raise ValueError("base summary does not contain the exact four frozen cells")
    survived = base._strict_boolean(
        survivors["population_survived"], label="base population_survived"
    )
    if not np.all(survived) or not np.all(
        survivors["final_alt_count"].to_numpy(dtype=int) > 0
    ):
        raise ValueError("base survivor table violates the population-survival filter")
    survivor_keys = {
        (str(row.demography_id), _coefficient_key(float(row.selection_coefficient)))
        for row in survivors.itertuples(index=False)
    }
    if not survivor_keys.issubset(expected_keys):
        raise ValueError("base survivor table contains an unexpected cell")
    if len(survivors) != int(summary["n_population_survived"].sum()):
        raise ValueError("base survivor and summary row counts disagree")

    outputs = completion["outputs"]
    relative = result_dir.relative_to(root).as_posix()
    return {
        "schema": base.ANALYSIS_SCHEMA,
        "relative_path": relative,
        "completion_sha256": base.sha256_file(completion_path),
        "contract_sha256": str(completion["contract_sha256"]),
        "plan_sha256": str(recorded_contract["plan_sha256"]),
        "cells": expected_cells,
        "summary": dict(outputs["summary"]),
        "survivors": dict(outputs["survivors"]),
        "implementation": current_implementation,
    }


def _load_bound_base_tables(
    plan: SensitivityPlan, expected: Mapping[str, Any]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    current = _base_bundle_contract(plan)
    if base._canonical_sha256(current) != base._canonical_sha256(expected):
        raise ValueError("published base bundle changed after sensitivity planning")
    result_dir = _base_result_dir(plan)
    summary = pd.read_csv(result_dir / base.ANALYSIS_OUTPUT_NAMES["summary"], sep="\t")
    survivors = pd.read_csv(
        result_dir / base.ANALYSIS_OUTPUT_NAMES["survivors"], sep="\t"
    )
    return summary, survivors


def _implementation_contract(repo_root: Path) -> dict[str, Any]:
    result: dict[str, Any] = {
        "hash_method": (
            "sha256(strict UTF-8 without BOM or bare CR; CRLF canonicalized to LF)"
        ),
        "base": base._implementation_contract(repo_root),
    }
    for relative in (MODULE_PATH, WRAPPER_PATH):
        path = repo_root / relative
        if not path.is_file() or path.is_symlink():
            raise FileNotFoundError(path)
        result[relative] = base.canonical_lf_sha256(path)
    return result


def _build_plan_payload(plan: SensitivityPlan) -> dict[str, Any]:
    plan.validate()
    root = plan.repo_root.resolve()
    eas = base.build_eas_schedule(root)
    han = base.build_han_schedule(validate_catalog=True)
    return {
        "schema": SCHEMA_VERSION,
        "status": "planned",
        "scientific_contract": {
            "generation_time_years": GENERATION_TIME_YEARS,
            "new_selection_coefficients": list(NEW_SELECTION_COEFFICIENTS),
            "combined_selection_coefficients": list(COMBINED_SELECTION_COEFFICIENTS),
            "dominance_coefficient": base.DOMINANCE_COEFFICIENT,
            "wright_fisher_scaling_factor": 1,
            "terminal_population_af_conditioning": False,
            "attempt_denominator": "all fixed attempted trajectories per new cell",
            "analysis_filter": "final population alternate count > 0",
            "fixation_counts_as_survival": True,
            "present_panel": {
                "diploids": base.PANEL_DIPLOIDS,
                "draw": "one Binomial(200, final_population_af), without resampling",
                "sample_detection_is_not_an_acceptance_condition": True,
            },
            "eas": eas.contract,
            "han": han.contract,
            "selection_timing": "first reproduction after origin or pulse through present",
        },
        "rng": {
            "base_seed": int(plan.base_seed),
            "batch_seed": "sha256(base_seed:cell_id:batch:batch_index)",
            "panel_seed": "sha256(batch_seed:present-panel)",
        },
        "cells": build_cells(plan),
        "base_bundle": _base_bundle_contract(plan),
        "paths": {
            "campaign_dir": plan.campaign_dir.resolve().relative_to(root).as_posix(),
            "work_subdir": DEFAULT_WORK_SUBDIR,
            "results_subdir": DEFAULT_RESULTS_SUBDIR,
            "base_results_subdir": BASE_RESULTS_SUBDIR,
        },
        "implementation": _implementation_contract(root),
        "software": base._software_versions(),
    }


def write_plan(plan: SensitivityPlan) -> Path:
    payload = _build_plan_payload(plan)
    path = plan.campaign_dir.resolve() / DEFAULT_WORK_SUBDIR / "plan.json"
    if path.is_file():
        if path.is_symlink():
            raise ValueError("sensitivity plan must not be a symlink")
        existing = json.loads(path.read_text(encoding="utf-8"))
        if base._canonical_sha256(existing) != base._canonical_sha256(payload):
            raise ValueError("existing sensitivity plan is incompatible")
    else:
        base._atomic_json(path, payload)
    return path


def _load_plan(path: str | Path) -> dict[str, Any]:
    plan_path = Path(path)
    if not plan_path.is_file() or plan_path.is_symlink():
        raise ValueError("sensitivity plan is unavailable")
    payload = json.loads(plan_path.read_text(encoding="utf-8"))
    if payload.get("schema") != SCHEMA_VERSION or payload.get("status") != "planned":
        raise ValueError("sensitivity plan schema/status is incompatible")
    cells = payload.get("cells")
    if not isinstance(cells, list) or len(cells) != 4:
        raise ValueError("sensitivity plan must contain exactly four new cells")
    return payload


def simulate_study(plan: SensitivityPlan, *, workers: int = 4) -> pd.DataFrame:
    if not 1 <= int(workers) <= 24:
        raise ValueError("workers must be between 1 and 24")
    plan_path = write_plan(plan)
    payload = _load_plan(plan_path)
    work_dir = plan.campaign_dir.resolve() / DEFAULT_WORK_SUBDIR
    tasks = [
        {
            "repo_root": str(plan.repo_root.resolve()),
            "work_dir": str(work_dir),
            "plan_payload": payload,
            "cell": cell,
            "batch_index": batch_index,
        }
        for cell in payload["cells"]
        for batch_index in range(int(cell["n_batches"]))
    ]
    if int(workers) == 1:
        rows = [base._simulate_batch_task(task) for task in tasks]
    else:
        rows = []
        with ProcessPoolExecutor(max_workers=int(workers)) as pool:
            futures = {
                pool.submit(base._simulate_batch_task, task): task for task in tasks
            }
            for future in as_completed(futures):
                rows.append(future.result())
    status = pd.DataFrame(rows).sort_values(["cell_id", "batch_index"])
    base._atomic_frame(work_dir / "last_simulation_status.tsv", status)
    return status.reset_index(drop=True)


def _expected_new_batch_inputs(
    plan: SensitivityPlan, payload: Mapping[str, Any]
) -> tuple[list[pd.DataFrame], list[dict[str, Any]]]:
    work_dir = plan.campaign_dir.resolve() / DEFAULT_WORK_SUBDIR
    frames: list[pd.DataFrame] = []
    records: list[dict[str, Any]] = []
    for cell in payload["cells"]:
        schedule = base._schedule_for_cell(plan.repo_root.resolve(), cell)
        for batch_index in range(int(cell["n_batches"])):
            start = batch_index * int(cell["batch_size"])
            stop = min(int(cell["attempts"]), start + int(cell["batch_size"]))
            seed = base._stable_seed(
                int(payload["rng"]["base_seed"]),
                f"{cell['cell_id']}:batch:{batch_index}",
            )
            contract = base._batch_contract(
                payload,
                cell,
                schedule,
                batch_index=batch_index,
                attempt_start=start,
                attempt_stop=stop,
                batch_seed=seed,
            )
            table, completion = base._batch_paths(
                work_dir, str(cell["cell_id"]), batch_index
            )
            frames.append(base._validate_batch(table, completion, contract))
            records.append(
                {
                    "cell_id": cell["cell_id"],
                    "batch_index": batch_index,
                    "table_sha256": base.sha256_file(table),
                    "completion_sha256": base.sha256_file(completion),
                }
            )
    return frames, records


def _validate_new_attempts(
    attempts: pd.DataFrame, cells: Sequence[Mapping[str, Any]]
) -> None:
    if tuple(attempts.columns) != ATTEMPT_COLUMNS:
        raise ValueError("new-attempt columns are incompatible")
    if attempts["trajectory_id"].duplicated().any():
        raise ValueError("new trajectory IDs are not unique")
    if set(attempts["simulation_class"].astype(str)) != {"selected"}:
        raise ValueError("sensitivity attempts must all be selected")
    keys = {
        (str(row.demography_id), _coefficient_key(float(row.selection_coefficient)))
        for row in attempts.itertuples(index=False)
    }
    expected = {
        (
            str(cell["demography_id"]),
            _coefficient_key(float(cell["selection_coefficient"])),
        )
        for cell in cells
    }
    if keys != expected or {key[1] for key in keys} != {1, 5}:
        raise ValueError("new attempts contain an unexpected coefficient or demography")
    expected_rows = {str(cell["cell_id"]): int(cell["attempts"]) for cell in cells}
    observed_rows = attempts.groupby("cell_id", sort=False).size().to_dict()
    if observed_rows != expected_rows:
        raise ValueError("new attempt counts do not match the sensitivity plan")


def build_comparisons(survivors: pd.DataFrame) -> pd.DataFrame:
    """Compute coefficient-specific survivor-AF AUCs against neutral."""

    survived = base._strict_boolean(
        survivors["population_survived"], label="population_survived"
    )
    if not np.all(survived):
        raise ValueError("comparison input must contain only population survivors")
    rows = []
    for demography in DEMOGRAPHY_LABELS:
        group = survivors[survivors["demography_id"] == demography]
        neutral = np.sort(
            group.loc[
                np.isclose(
                    group["selection_coefficient"].to_numpy(float),
                    0.0,
                    rtol=0,
                    atol=1e-15,
                ),
                "final_population_af",
            ].to_numpy(float)
        )
        for coefficient in NEW_SELECTION_COEFFICIENTS + (0.01,):
            selected = group.loc[
                np.isclose(
                    group["selection_coefficient"].to_numpy(float),
                    coefficient,
                    rtol=0,
                    atol=1e-15,
                ),
                "final_population_af",
            ].to_numpy(float)
            if len(neutral) and len(selected):
                less = np.searchsorted(neutral, selected, side="left")
                equal = np.searchsorted(neutral, selected, side="right") - less
                auc = float(np.mean((less + 0.5 * equal) / len(neutral)))
            else:
                auc = np.nan
            rows.append(
                {
                    "demography_id": demography,
                    "neutral_selection_coefficient": 0.0,
                    "selected_selection_coefficient": coefficient,
                    "statistic": (
                        "P(AF_selected > AF_s0 | population survival) + half ties"
                    ),
                    "population_af_roc_auc_survivors": auc,
                    "n_neutral_survivors": len(neutral),
                    "n_selected_survivors": len(selected),
                }
            )
    return pd.DataFrame(rows, columns=COMPARISON_COLUMNS)


def combine_analysis_tables(
    new_attempts: pd.DataFrame,
    base_summary: pd.DataFrame,
    base_survivors: pd.DataFrame,
    cells: Sequence[Mapping[str, Any]],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return combined summary, combined survivors, and non-pooled AUC rows."""

    _validate_new_attempts(new_attempts, cells)
    if (
        tuple(base_summary.columns) != SUMMARY_COLUMNS
        or tuple(base_survivors.columns) != ATTEMPT_COLUMNS
    ):
        raise ValueError("base tables have incompatible columns")
    new_summary, _ = base.summarize_attempts(new_attempts)
    new_survived = base._strict_boolean(
        new_attempts["population_survived"], label="new population_survived"
    )
    new_survivors = new_attempts[new_survived].copy()
    combined_summary = (
        pd.concat([base_summary, new_summary], ignore_index=True)
        .sort_values(["demography_id", "selection_coefficient"])
        .reset_index(drop=True)
    )
    combined_survivors = (
        pd.concat([base_survivors, new_survivors], ignore_index=True)
        .sort_values(["cell_id", "attempt_index"])
        .reset_index(drop=True)
    )
    if combined_survivors["trajectory_id"].duplicated().any():
        raise ValueError("combined survivor trajectory IDs are not unique")
    expected_summary_keys = {
        (demography, _coefficient_key(coefficient))
        for demography in DEMOGRAPHY_LABELS
        for coefficient in COMBINED_SELECTION_COEFFICIENTS
    }
    observed_summary_keys = {
        (str(row.demography_id), _coefficient_key(float(row.selection_coefficient)))
        for row in combined_summary.itertuples(index=False)
    }
    if len(combined_summary) != 8 or observed_summary_keys != expected_summary_keys:
        raise ValueError("combined summary must contain the exact eight study cells")
    survived = base._strict_boolean(
        combined_survivors["population_survived"],
        label="combined population_survived",
    )
    if not np.all(survived) or not np.all(
        combined_survivors["final_alt_count"].to_numpy(dtype=int) > 0
    ):
        raise ValueError("combined survivor table has a non-surviving trajectory")
    if len(combined_survivors) != int(combined_summary["n_population_survived"].sum()):
        raise ValueError("combined survivor and summary counts disagree")
    comparisons = build_comparisons(combined_survivors)
    return combined_summary, combined_survivors, comparisons


SavePair = Callable[[plt.Figure, Path, str], list[Path]]


def _save_figure_pair(figure: plt.Figure, result_dir: Path, stem: str) -> list[Path]:
    paths = []
    for extension in ("png", "pdf"):
        path = result_dir / f"{stem}.{extension}"
        temporary = path.with_name(f".{path.stem}.tmp.{os.getpid()}{path.suffix}")
        figure.savefig(
            temporary,
            dpi=base.PNG_DPI if extension == "png" else None,
            metadata={"Creator": "gamma_smc_aou.natural_final_af_sensitivity"},
        )
        base._validate_letter_geometry(temporary)
        os.replace(temporary, path)
        paths.append(path)
    return paths


def _plot_unconditional_fates(
    summary: pd.DataFrame,
    result_dir: Path,
    *,
    save_pair: SavePair = _save_figure_pair,
) -> list[Path]:
    fate_columns = (
        ("lost_fraction", "Lost", "#B9B9B9"),
        ("segregating_fraction", "Segregating", "#4C78A8"),
        ("fixed_fraction", "Fixed", "#D1495B"),
    )
    figure, axes = plt.subplots(1, 2, figsize=base.FIGURE_SIZE_INCHES, sharey=True)
    x = np.arange(len(COMBINED_SELECTION_COEFFICIENTS))
    for axis, demography in zip(axes, DEMOGRAPHY_LABELS, strict=True):
        cell = summary[summary["demography_id"] == demography].set_index(
            "selection_coefficient"
        )
        bottom = np.zeros(len(x))
        for column, _, color in fate_columns:
            values = np.asarray(
                [
                    float(cell.loc[coefficient, column])
                    for coefficient in COMBINED_SELECTION_COEFFICIENTS
                ]
            )
            axis.bar(x, values, bottom=bottom, width=0.70, color=color)
            bottom += values
        axis.set_xticks(
            x,
            tuple(
                _coefficient_label(value) for value in COMBINED_SELECTION_COEFFICIENTS
            ),
            fontsize=13,
        )
        axis.set_ylim(0, 1.0)
        axis.set_title(DEMOGRAPHY_LABELS[demography], fontsize=18)
        axis.tick_params(axis="y", labelsize=13)
        axis.grid(axis="y", alpha=0.2)
        for x_value, coefficient in enumerate(COMBINED_SELECTION_COEFFICIENTS):
            survival = float(cell.loc[coefficient, "population_survival_fraction"])
            fixed = float(cell.loc[coefficient, "fixation_fraction_among_survivors"])
            axis.text(
                x_value,
                0.025,
                f"survive\n{100 * survival:.3g}%\nfixed|survive\n{100 * fixed:.3g}%",
                ha="center",
                va="bottom",
                fontsize=9.5,
            )
        attempt_pairs = (
            COMBINED_SELECTION_COEFFICIENTS[:2],
            COMBINED_SELECTION_COEFFICIENTS[2:],
        )
        attempts = "\n".join(
            "   ".join(
                f"{_coefficient_label(coefficient)}: "
                f"n={int(cell.loc[coefficient, 'n_attempted']):,}"
                for coefficient in pair
            )
            for pair in attempt_pairs
        )
        axis.text(
            0.5,
            -0.16,
            attempts,
            transform=axis.transAxes,
            ha="center",
            fontsize=10.5,
        )
    axes[0].set_ylabel("Fraction of all attempted trajectories", fontsize=16)
    figure.legend(
        handles=[
            Patch(facecolor=color, label=label) for _, label, color in fate_columns
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, 0.91),
        ncol=3,
        fontsize=14,
        frameon=False,
    )
    figure.suptitle("Unconditional focal-allele fates", fontsize=20, y=0.98)
    figure.subplots_adjust(top=0.84, bottom=0.25, left=0.10, right=0.98, wspace=0.16)
    paths = save_pair(figure, result_dir, "unconditional_fate_fractions")
    plt.close(figure)
    return paths


def _ecdf(
    axis: plt.Axes, values: np.ndarray, *, color: str, linestyle: str = "-"
) -> None:
    ordered = np.sort(np.asarray(values, dtype=float))
    if len(ordered):
        x = np.append(ordered, 1.0)
        y = np.append(np.arange(1, len(ordered) + 1) / len(ordered), 1.0)
        axis.step(x, y, where="post", lw=2.7, ls=linestyle, color=color)


def _plot_ecdf(
    survivors: pd.DataFrame,
    result_dir: Path,
    *,
    save_pair: SavePair = _save_figure_pair,
) -> list[Path]:
    figure, axes = plt.subplots(1, 2, figsize=base.FIGURE_SIZE_INCHES, sharey=True)
    for axis, demography in zip(axes, DEMOGRAPHY_LABELS, strict=True):
        group = survivors[survivors["demography_id"] == demography]
        annotation = []
        for coefficient in COMBINED_SELECTION_COEFFICIENTS:
            values = group.loc[
                np.isclose(
                    group["selection_coefficient"].to_numpy(float),
                    coefficient,
                    rtol=0,
                    atol=1e-15,
                ),
                "final_population_af",
            ].to_numpy(float)
            _ecdf(axis, values, color=COEFFICIENT_COLORS[coefficient])
            fixed = float(np.mean(values == 1.0)) if len(values) else np.nan
            annotation.append(
                f"{_coefficient_label(coefficient)}: n={len(values):,}; "
                f"fixed={100 * fixed:.2f}%"
            )
        axis.set_xscale("symlog", linthresh=1e-4, base=10)
        axis.set_xlim(0, 1)
        axis.set_ylim(0, 1.01)
        axis.set_title(DEMOGRAPHY_LABELS[demography], fontsize=18)
        axis.tick_params(labelsize=13)
        axis.grid(alpha=0.2)
        axis.text(
            0.03,
            0.96,
            "\n".join(annotation),
            transform=axis.transAxes,
            va="top",
            fontsize=10.5,
        )
    axes[0].set_ylabel("Empirical CDF among population survivors", fontsize=16)
    figure.legend(
        handles=[
            Line2D(
                [0],
                [0],
                color=COEFFICIENT_COLORS[coefficient],
                lw=3,
                label=_coefficient_label(coefficient),
            )
            for coefficient in COMBINED_SELECTION_COEFFICIENTS
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, 0.87),
        ncol=4,
        fontsize=13,
        frameon=False,
    )
    figure.suptitle(
        "Natural final-AF selection sensitivity\n(no terminal-frequency conditioning)",
        fontsize=20,
        y=0.98,
    )
    figure.supxlabel("Present population allele frequency", fontsize=16, y=0.045)
    figure.subplots_adjust(top=0.79, bottom=0.15, left=0.10, right=0.98, wspace=0.16)
    paths = save_pair(figure, result_dir, "final_population_af_survivor_ecdf")
    plt.close(figure)
    return paths


def _plot_population_vs_sample(
    survivors: pd.DataFrame,
    result_dir: Path,
    *,
    save_pair: SavePair = _save_figure_pair,
) -> list[Path]:
    figure, axes = plt.subplots(1, 2, figsize=base.FIGURE_SIZE_INCHES, sharey=True)
    for axis, demography in zip(axes, DEMOGRAPHY_LABELS, strict=True):
        group = survivors[survivors["demography_id"] == demography]
        annotation = []
        for coefficient in COMBINED_SELECTION_COEFFICIENTS:
            cell = group[
                np.isclose(
                    group["selection_coefficient"].to_numpy(float),
                    coefficient,
                    rtol=0,
                    atol=1e-15,
                )
            ]
            annotation.append(f"{_coefficient_label(coefficient)}: n={len(cell):,}")
            _ecdf(
                axis,
                cell["final_population_af"].to_numpy(float),
                color=COEFFICIENT_COLORS[coefficient],
            )
            _ecdf(
                axis,
                cell["sample_af"].to_numpy(float),
                color=COEFFICIENT_COLORS[coefficient],
                linestyle="--",
            )
        axis.set_xscale("symlog", linthresh=1e-4, base=10)
        axis.set_xlim(0, 1)
        axis.set_ylim(0, 1.01)
        axis.set_title(DEMOGRAPHY_LABELS[demography], fontsize=18)
        axis.tick_params(labelsize=13)
        axis.grid(alpha=0.2)
        axis.text(
            0.03,
            0.96,
            "\n".join(annotation),
            transform=axis.transAxes,
            va="top",
            fontsize=10.5,
        )
    axes[0].set_ylabel("Empirical CDF among population survivors", fontsize=16)
    handles = [
        Line2D(
            [0],
            [0],
            color=COEFFICIENT_COLORS[coefficient],
            lw=2.7,
            label=_coefficient_label(coefficient),
        )
        for coefficient in COMBINED_SELECTION_COEFFICIENTS
    ]
    handles.extend(
        (
            Line2D([0], [0], color="black", lw=2.7, ls="-", label="population"),
            Line2D(
                [0],
                [0],
                color="black",
                lw=2.7,
                ls="--",
                label="100-diploid panel",
            ),
        )
    )
    figure.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.91),
        ncol=3,
        fontsize=11.5,
        frameon=False,
    )
    figure.suptitle(
        "Population AF versus one unconditioned 100-diploid panel", fontsize=20, y=0.98
    )
    figure.supxlabel("Present allele frequency", fontsize=16, y=0.045)
    figure.subplots_adjust(top=0.80, bottom=0.15, left=0.10, right=0.98, wspace=0.16)
    paths = save_pair(figure, result_dir, "population_vs_sample_af_survivor_ecdf")
    plt.close(figure)
    return paths


def _frame_tsv_sha256(frame: pd.DataFrame) -> str:
    return base._frame_tsv_sha256(frame)


def _describe_output(
    label: str, path: Path, *, table: pd.DataFrame | None = None
) -> dict[str, Any]:
    record: dict[str, Any] = {
        "path": path.name,
        "sha256": base.sha256_file(path),
        "size_bytes": path.stat().st_size,
    }
    if label in TABLE_COLUMNS:
        if table is None:
            rows, columns = base._tsv_shape(path)
        else:
            rows, columns = len(table), list(table.columns)
        record.update({"rows": int(rows), "columns": list(columns)})
    elif path.suffix in {".png", ".pdf"}:
        record["geometry"] = base._validate_letter_geometry(path)
    return record


def _validate_analysis_bundle(
    plan: SensitivityPlan,
    result_dir: Path,
    contract: Mapping[str, Any],
) -> dict[str, Any]:
    completion_path = result_dir / "analysis_completion.json"
    if (
        not result_dir.is_dir()
        or result_dir.is_symlink()
        or not completion_path.is_file()
        or completion_path.is_symlink()
    ):
        raise ValueError("sensitivity result bundle is incomplete")
    completion = json.loads(completion_path.read_text(encoding="utf-8"))
    outputs = completion.get("outputs")
    expected_files = set(ANALYSIS_OUTPUT_NAMES.values()) | {completion_path.name}
    if (
        completion.get("schema") != ANALYSIS_SCHEMA
        or completion.get("status") != "complete"
        or not isinstance(completion.get("contract"), Mapping)
        or completion.get("contract_sha256")
        != base._canonical_sha256(completion["contract"])
        or base._canonical_sha256(completion["contract"])
        != base._canonical_sha256(contract)
        or not isinstance(outputs, Mapping)
        or set(outputs) != set(ANALYSIS_OUTPUT_NAMES)
        or {entry.name for entry in result_dir.iterdir()} != expected_files
    ):
        raise ValueError("sensitivity analysis completion contract is invalid")
    if base._canonical_sha256(_base_bundle_contract(plan)) != base._canonical_sha256(
        contract["base_bundle"]
    ):
        raise ValueError("sensitivity base-bundle dependency changed")
    if base._canonical_sha256(_implementation_contract(plan.repo_root.resolve())) != (
        base._canonical_sha256(contract["implementation"])
    ):
        raise ValueError("sensitivity implementation changed")

    table_shapes: dict[str, tuple[int, list[str]]] = {}
    for label, expected_name in ANALYSIS_OUTPUT_NAMES.items():
        output = outputs[label]
        path = result_dir / expected_name
        if (
            not isinstance(output, Mapping)
            or output.get("path") != expected_name
            or not path.is_file()
            or path.is_symlink()
            or output.get("sha256") != base.sha256_file(path)
            or int(output.get("size_bytes", -1)) != path.stat().st_size
        ):
            raise ValueError(f"cached sensitivity output failed: {label}")
        if label in TABLE_COLUMNS:
            rows, columns = base._tsv_shape(path)
            table_shapes[label] = (rows, columns)
            if (
                int(output.get("rows", -1)) != rows
                or output.get("columns") != columns
                or columns != list(TABLE_COLUMNS[label])
                or base._tsv_payload_sha256(path)
                != contract["table_payload_sha256"][label]
            ):
                raise ValueError(f"cached sensitivity table failed: {label}")
        elif path.suffix in {".png", ".pdf"}:
            geometry = base._validate_letter_geometry(path)
            if output.get("geometry") != geometry:
                raise ValueError(f"cached sensitivity figure failed: {label}")

    attempts = pd.read_csv(result_dir / ANALYSIS_OUTPUT_NAMES["attempts"], sep="\t")
    survivors = pd.read_csv(result_dir / ANALYSIS_OUTPUT_NAMES["survivors"], sep="\t")
    summary = pd.read_csv(result_dir / ANALYSIS_OUTPUT_NAMES["summary"], sep="\t")
    comparisons = pd.read_csv(
        result_dir / ANALYSIS_OUTPUT_NAMES["comparisons"], sep="\t"
    )
    _validate_new_attempts(attempts, contract["cells"])
    if (
        table_shapes["attempts"][0]
        != sum(int(cell["attempts"]) for cell in contract["cells"])
        or table_shapes["summary"][0] != 8
        or table_shapes["comparisons"][0] != 6
        or table_shapes["survivors"][0] != int(summary["n_population_survived"].sum())
    ):
        raise ValueError("sensitivity table row contract failed")
    if not np.array_equal(
        summary["n_lost"].to_numpy(dtype=int)
        + summary["n_segregating"].to_numpy(dtype=int)
        + summary["n_fixed"].to_numpy(dtype=int),
        summary["n_attempted"].to_numpy(dtype=int),
    ) or not np.array_equal(
        summary["n_segregating"].to_numpy(dtype=int)
        + summary["n_fixed"].to_numpy(dtype=int),
        summary["n_population_survived"].to_numpy(dtype=int),
    ):
        raise ValueError("sensitivity summary fate totals disagree")
    expected_cells = list(contract["base_bundle"]["cells"]) + list(contract["cells"])
    expected_summary = {
        (
            str(cell["demography_id"]),
            _coefficient_key(float(cell["selection_coefficient"])),
        ): (int(cell["attempts"]), str(cell["simulation_class"]))
        for cell in expected_cells
    }
    observed_summary = {
        (
            str(row.demography_id),
            _coefficient_key(float(row.selection_coefficient)),
        ): (int(row.n_attempted), str(row.simulation_class))
        for row in summary.itertuples(index=False)
    }
    if observed_summary != expected_summary:
        raise ValueError("sensitivity summary cells do not match their bound inputs")
    survived = base._strict_boolean(
        survivors["population_survived"], label="population_survived"
    )
    if not np.all(survived) or not np.all(
        survivors["final_alt_count"].to_numpy(dtype=int) > 0
    ):
        raise ValueError("sensitivity survivor filter is invalid")
    if survivors["trajectory_id"].duplicated().any():
        raise ValueError("sensitivity survivor trajectory IDs are not unique")
    expected_survivor_counts = {
        (
            str(row.demography_id),
            _coefficient_key(float(row.selection_coefficient)),
        ): int(row.n_population_survived)
        for row in summary.itertuples(index=False)
    }
    observed_survivor_counts = {
        (str(demography), _coefficient_key(float(coefficient))): int(count)
        for (demography, coefficient), count in survivors.groupby(
            ["demography_id", "selection_coefficient"], sort=False
        )
        .size()
        .items()
    }
    if any(
        observed_survivor_counts.get(key, 0) != count
        for key, count in expected_survivor_counts.items()
    ) or not set(observed_survivor_counts).issubset(expected_survivor_counts):
        raise ValueError("sensitivity survivor counts disagree by cell")
    expected_pairs = {
        (demography, _coefficient_key(coefficient))
        for demography in DEMOGRAPHY_LABELS
        for coefficient in NEW_SELECTION_COEFFICIENTS + (0.01,)
    }
    observed_pairs = {
        (
            str(row.demography_id),
            _coefficient_key(float(row.selected_selection_coefficient)),
        )
        for row in comparisons.itertuples(index=False)
    }
    if observed_pairs != expected_pairs or not np.allclose(
        comparisons["neutral_selection_coefficient"].to_numpy(float), 0.0
    ):
        raise ValueError("sensitivity comparison rows are pooled or incomplete")
    rebuilt_comparisons = build_comparisons(survivors).sort_values(
        ["demography_id", "selected_selection_coefficient"]
    )
    cached_comparisons = comparisons.sort_values(
        ["demography_id", "selected_selection_coefficient"]
    )
    if not np.array_equal(
        rebuilt_comparisons[["n_neutral_survivors", "n_selected_survivors"]].to_numpy(
            dtype=int
        ),
        cached_comparisons[["n_neutral_survivors", "n_selected_survivors"]].to_numpy(
            dtype=int
        ),
    ) or not np.allclose(
        rebuilt_comparisons["population_af_roc_auc_survivors"].to_numpy(float),
        cached_comparisons["population_af_roc_auc_survivors"].to_numpy(float),
        rtol=0,
        atol=5e-12,
        equal_nan=True,
    ):
        raise ValueError("sensitivity comparison values do not match survivor AFs")
    return dict(completion)


def _remove_owned_staging(stage: Path, parent: Path) -> None:
    resolved_stage = stage.resolve()
    resolved_parent = parent.resolve()
    if resolved_stage.parent != resolved_parent or not resolved_stage.name.startswith(
        ".natural_final_af_selection_sensitivity.staging."
    ):
        raise RuntimeError("refusing to remove an unowned sensitivity staging path")
    if resolved_stage.exists():
        shutil.rmtree(resolved_stage)


def analyze_study(plan: SensitivityPlan) -> dict[str, Any]:
    plan_path = write_plan(plan)
    payload = _load_plan(plan_path)
    frames, batch_inputs = _expected_new_batch_inputs(plan, payload)
    new_attempts = pd.concat(frames, ignore_index=True).sort_values(
        ["cell_id", "attempt_index"]
    )
    base_summary, base_survivors = _load_bound_base_tables(plan, payload["base_bundle"])
    summary, survivors, comparisons = combine_analysis_tables(
        new_attempts, base_summary, base_survivors, payload["cells"]
    )
    tables = {
        "attempts": new_attempts,
        "survivors": survivors,
        "summary": summary,
        "comparisons": comparisons,
    }
    result_dir = plan.campaign_dir.resolve() / DEFAULT_RESULTS_SUBDIR
    contract = {
        "schema": ANALYSIS_SCHEMA,
        "plan_sha256": base.sha256_file(plan_path),
        "base_bundle": payload["base_bundle"],
        "batch_inputs": batch_inputs,
        "cells": payload["cells"],
        "combined_selection_coefficients": list(COMBINED_SELECTION_COEFFICIENTS),
        "survival_filter": "final_alt_count > 0",
        "sample_panel_is_not_a_filter": True,
        "implementation": _implementation_contract(plan.repo_root.resolve()),
        "software": base._software_versions(),
        "table_payload_sha256": {
            label: _frame_tsv_sha256(frame) for label, frame in tables.items()
        },
    }
    if result_dir.exists():
        if not (result_dir / "analysis_completion.json").is_file():
            raise ValueError("incompatible partial sensitivity result bundle exists")
        completion = _validate_analysis_bundle(plan, result_dir, contract)
        return {**completion, "cache_hit": True}
    parent = result_dir.parent
    parent.mkdir(parents=True, exist_ok=True)
    stage = Path(
        tempfile.mkdtemp(
            prefix=".natural_final_af_selection_sensitivity.staging.", dir=parent
        )
    )
    try:
        for label, frame in tables.items():
            base._atomic_frame(stage / ANALYSIS_OUTPUT_NAMES[label], frame)
        _plot_unconditional_fates(summary, stage)
        _plot_ecdf(survivors, stage)
        _plot_population_vs_sample(survivors, stage)
        base._atomic_text(
            stage / ANALYSIS_OUTPUT_NAMES["readme"],
            "# Natural final-AF selection sensitivity\n\n"
            "This additive extension simulates only `s=0.001` and `s=0.005`. "
            "The published `natural_final_af` bundle supplies checksum-validated "
            "`s=0` and `s=0.01` summary and survivor rows; its source and artifacts "
            "are not modified or resimulated. The all-attempt table in this directory "
            "therefore contains only the four new cells, while the survivor, summary, "
            "comparison, and figure artifacts combine all four coefficients.\n\n"
            "All fixed attempted Wright--Fisher trajectories are retained. Population "
            "survival (`final_alt_count > 0`) is the only analysis filter; there is no "
            "terminal-AF gate, fixation counts as survival, and the one Binomial(200, "
            "population AF) panel draw is never an acceptance condition. EAS begins "
            "from one de novo copy 2,000 generations ago. Han uses the donor-fixed "
            "pulse approximation at 2.96%. Selection is additive (`h=0.5`) from "
            "origin/pulse through the present. AUC rows compare each positive "
            "coefficient separately with neutral among population survivors.\n",
        )
        outputs = {
            label: _describe_output(label, stage / name, table=tables.get(label))
            for label, name in ANALYSIS_OUTPUT_NAMES.items()
        }
        completion = {
            "schema": ANALYSIS_SCHEMA,
            "status": "complete",
            "contract": contract,
            "contract_sha256": base._canonical_sha256(contract),
            "outputs": outputs,
        }
        base._atomic_json(stage / "analysis_completion.json", completion)
        _validate_analysis_bundle(plan, stage, contract)
        if result_dir.exists():
            raise ValueError("sensitivity result bundle appeared during publication")
        os.replace(stage, result_dir)
        published = _validate_analysis_bundle(plan, result_dir, contract)
        return {**published, "cache_hit": False}
    finally:
        if stage.exists():
            _remove_owned_staging(stage, parent)


def status_table(plan: SensitivityPlan) -> pd.DataFrame:
    plan_path = write_plan(plan)
    payload = _load_plan(plan_path)
    work_dir = plan.campaign_dir.resolve() / DEFAULT_WORK_SUBDIR
    rows = []
    for cell in payload["cells"]:
        schedule = base._schedule_for_cell(plan.repo_root.resolve(), cell)
        for batch_index in range(int(cell["n_batches"])):
            table, completion = base._batch_paths(
                work_dir, str(cell["cell_id"]), batch_index
            )
            start = batch_index * int(cell["batch_size"])
            stop = min(int(cell["attempts"]), start + int(cell["batch_size"]))
            seed = base._stable_seed(
                int(payload["rng"]["base_seed"]),
                f"{cell['cell_id']}:batch:{batch_index}",
            )
            contract = base._batch_contract(
                payload,
                cell,
                schedule,
                batch_index=batch_index,
                attempt_start=start,
                attempt_stop=stop,
                batch_seed=seed,
            )
            error = ""
            if table.is_file() and not completion.is_file():
                status = "restartable_table_only"
            elif completion.is_file() and not table.is_file():
                status = "invalid"
                error = "completion exists without table"
            elif not table.is_file() and not completion.is_file():
                status = "missing"
            else:
                try:
                    base._validate_batch(table, completion, contract)
                    status = "complete"
                except (ValueError, OSError, json.JSONDecodeError) as exception:
                    status = "invalid"
                    error = str(exception)
            rows.append(
                {
                    "cell_id": cell["cell_id"],
                    "batch_index": batch_index,
                    "status": status,
                    "error": error,
                }
            )
    return pd.DataFrame(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Additive natural final-AF selection sensitivity"
    )
    parser.add_argument(
        "phase", choices=("plan", "simulate", "analyze", "status", "all")
    )
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--campaign-dir", type=Path)
    parser.add_argument(
        "--eas-s001-attempts", type=int, default=DEFAULT_EAS_S001_ATTEMPTS
    )
    parser.add_argument(
        "--eas-s005-attempts", type=int, default=DEFAULT_EAS_S005_ATTEMPTS
    )
    parser.add_argument(
        "--han-s001-attempts", type=int, default=DEFAULT_HAN_S001_ATTEMPTS
    )
    parser.add_argument(
        "--han-s005-attempts", type=int, default=DEFAULT_HAN_S005_ATTEMPTS
    )
    parser.add_argument("--eas-batch-size", type=int, default=DEFAULT_EAS_BATCH_SIZE)
    parser.add_argument("--han-batch-size", type=int, default=DEFAULT_HAN_BATCH_SIZE)
    parser.add_argument("--base-seed", type=int, default=BASE_SEED)
    parser.add_argument("--workers", type=int, default=4)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    root = args.repo_root.resolve()
    campaign = (args.campaign_dir or root / DEFAULT_CAMPAIGN_DIR).resolve()
    plan = SensitivityPlan(
        repo_root=root,
        campaign_dir=campaign,
        eas_s001_attempts=args.eas_s001_attempts,
        eas_s005_attempts=args.eas_s005_attempts,
        han_s001_attempts=args.han_s001_attempts,
        han_s005_attempts=args.han_s005_attempts,
        eas_batch_size=args.eas_batch_size,
        han_batch_size=args.han_batch_size,
        base_seed=args.base_seed,
    )
    write_plan(plan)
    if args.phase == "plan":
        return 0
    if args.phase in {"simulate", "all"}:
        status = simulate_study(plan, workers=args.workers)
        print(status.groupby("status").size().to_string())
    if args.phase in {"analyze", "all"}:
        completion = analyze_study(plan)
        print(
            json.dumps(
                {"status": completion["status"], "cache_hit": completion["cache_hit"]}
            )
        )
    if args.phase == "status":
        status = status_table(plan)
        print(status.groupby("status").size().to_string())
    return 0


__all__ = [
    "SensitivityPlan",
    "analyze_study",
    "build_cells",
    "build_comparisons",
    "build_parser",
    "combine_analysis_tables",
    "main",
    "simulate_study",
    "status_table",
    "write_plan",
]
