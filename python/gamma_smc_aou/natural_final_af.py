"""Unconditioned single-locus Wright--Fisher final-AF sensitivity study.

The four cells cross EAS versus Han demography with ``s=0`` versus ``s=0.01``.
Every attempted trajectory is retained.  Population survival (terminal allele
count greater than zero) is an analysis filter, never a simulation acceptance
gate.  A single independent 100-diploid panel is sampled at the terminal
population frequency without rejection or resampling.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import math
import os
import platform
import re
import shutil
import struct
import tempfile
from collections.abc import Mapping, Sequence
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from importlib.metadata import version as distribution_version
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from .eas_sweep_models import (
    GENERATION_TIME_YEARS,
    load_ancient_eurasia_model,
    load_phlash_eas_npz,
)
from .eas_sweep_study import BASE_SEED, EAS_RESOURCE_SHA256


SCHEMA_VERSION = "gamma-smc.natural-final-af/v1"
BATCH_SCHEMA = "gamma-smc.natural-final-af-batch/v1"
ANALYSIS_SCHEMA = "gamma-smc.natural-final-af-analysis/v1"
MODULE_PATH = "python/gamma_smc_aou/natural_final_af.py"
WRAPPER_PATH = "scripts/run_natural_final_af.py"
EAS_RESOURCE_PATH = "no_introgression_EAS_sim/resources/EAS.npz"
DEFAULT_CAMPAIGN_DIR = "focused_selection_EAS_sim"
DEFAULT_WORK_SUBDIR = "work/natural_final_af"
DEFAULT_RESULTS_SUBDIR = "results/natural_final_af"

EAS_ORIGIN_GENERATIONS = 2_000
HAN_PULSE_GENERATIONS = 2_272
HAN_PRE_SPLIT_GENERATIONS = 256
HAN_POST_SPLIT_GENERATIONS = 2_016
HAN_LOSCHBOUR_N = 2_340
HAN_PRESENT_N = 6_300
HAN_PULSE_PROPORTION = 0.0296
DOMINANCE_COEFFICIENT = 0.5
PANEL_DIPLOIDS = 100
SELECTION_COEFFICIENTS = (0.0, 0.01)
FIGURE_SIZE_INCHES = (11.0, 8.5)
PNG_DPI = 300

# These are fixed attempted trajectories, not requested survivor counts.
DEFAULT_EAS_NEUTRAL_ATTEMPTS = 2_000_000
DEFAULT_EAS_SELECTED_ATTEMPTS = 200_000
DEFAULT_HAN_NEUTRAL_ATTEMPTS = 20_000
DEFAULT_HAN_SELECTED_ATTEMPTS = 20_000
DEFAULT_EAS_NEUTRAL_BATCH_SIZE = 50_000
DEFAULT_EAS_SELECTED_BATCH_SIZE = 10_000
DEFAULT_HAN_BATCH_SIZE = 1_000

ATTEMPT_COLUMNS = (
    "trajectory_id",
    "cell_id",
    "demography_id",
    "population",
    "simulation_class",
    "selection_coefficient",
    "dominance_coefficient",
    "attempt_index",
    "batch_index",
    "batch_seed",
    "panel_seed",
    "origin_generations_ago",
    "initial_population",
    "initial_diploid_n",
    "initial_alt_count",
    "initial_population_af",
    "present_diploid_n",
    "final_alt_count",
    "final_population_af",
    "focal_allele_outcome",
    "population_survived",
    "absorption_generations_ago",
    "panel_diploids",
    "sample_alt_count",
    "sample_af",
    "sample_detected",
)

SUMMARY_COLUMNS = (
    "demography_id",
    "population",
    "simulation_class",
    "selection_coefficient",
    "n_attempted",
    "n_lost",
    "n_segregating",
    "n_fixed",
    "lost_fraction",
    "segregating_fraction",
    "fixed_fraction",
    "n_population_survived",
    "population_survival_fraction",
    "population_survival_wilson_ci_low",
    "population_survival_wilson_ci_high",
    "fixation_fraction_among_survivors",
    "initial_population_af_mean_all_attempts",
    "final_population_af_mean_all_attempts",
    "neutral_martingale_difference",
    "neutral_martingale_mc_se",
    "neutral_martingale_z",
    "n_sample_detected_among_survivors",
    "sample_detection_fraction_among_survivors",
    "population_af_mean_survivors",
    "population_af_median_survivors",
    "population_af_q025_survivors",
    "population_af_q25_survivors",
    "population_af_q75_survivors",
    "population_af_q975_survivors",
    "sample_af_mean_survivors",
    "sample_af_median_survivors",
)

COMPARISON_COLUMNS = (
    "demography_id",
    "statistic",
    "population_af_roc_auc_survivors",
    "n_neutral_survivors",
    "n_selected_survivors",
)

SCHEDULE_COLUMNS = (
    "demography_id",
    "generations_ago",
    "population",
    "diploid_n",
)

ANALYSIS_OUTPUT_NAMES = {
    "attempts": "final_af_attempts.tsv.gz",
    "survivors": "final_af_population_survivors.tsv.gz",
    "summary": "final_af_summary.tsv",
    "comparisons": "final_af_comparisons.tsv",
    "schedules": "integer_demography_schedules.tsv.gz",
    "readme": "README.md",
    "fate_png": "unconditional_fate_fractions.png",
    "fate_pdf": "unconditional_fate_fractions.pdf",
    "population_ecdf_png": "final_population_af_survivor_ecdf.png",
    "population_ecdf_pdf": "final_population_af_survivor_ecdf.pdf",
    "population_sample_ecdf_png": "population_vs_sample_af_survivor_ecdf.png",
    "population_sample_ecdf_pdf": "population_vs_sample_af_survivor_ecdf.pdf",
}


@dataclass(frozen=True)
class TrajectorySchedule:
    demography_id: str
    ages: np.ndarray
    diploid_sizes: np.ndarray
    populations: tuple[str, ...]
    initial_mode: str
    initial_frequency: float | None
    contract: dict[str, Any]

    @property
    def origin_generations_ago(self) -> int:
        return int(self.ages[0])


@dataclass(frozen=True)
class NaturalAfPlan:
    repo_root: Path
    campaign_dir: Path
    eas_neutral_attempts: int = DEFAULT_EAS_NEUTRAL_ATTEMPTS
    eas_selected_attempts: int = DEFAULT_EAS_SELECTED_ATTEMPTS
    han_neutral_attempts: int = DEFAULT_HAN_NEUTRAL_ATTEMPTS
    han_selected_attempts: int = DEFAULT_HAN_SELECTED_ATTEMPTS
    eas_neutral_batch_size: int = DEFAULT_EAS_NEUTRAL_BATCH_SIZE
    eas_selected_batch_size: int = DEFAULT_EAS_SELECTED_BATCH_SIZE
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
            raise ValueError("natural-final-AF work must not be inside OneDrive")
        integer_fields = (
            "eas_neutral_attempts",
            "eas_selected_attempts",
            "han_neutral_attempts",
            "han_selected_attempts",
            "eas_neutral_batch_size",
            "eas_selected_batch_size",
            "han_batch_size",
            "base_seed",
        )
        for field in integer_fields:
            value = getattr(self, field)
            if isinstance(value, bool) or int(value) != value or int(value) < 1:
                raise ValueError(f"{field} must be a positive integer")


def _canonical_sha256(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(
            value, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode()
    ).hexdigest()


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_lf_sha256(path: str | Path) -> str:
    """Hash strict UTF-8 source after canonicalizing CRLF to LF.

    BOMs and bare carriage returns are rejected rather than normalized so the
    portable source identity cannot conceal an encoding or newline defect.
    """

    source = Path(path)
    raw = source.read_bytes()
    if raw.startswith(b"\xef\xbb\xbf"):
        raise ValueError(f"source has a UTF-8 BOM: {source}")
    try:
        value = raw.decode("utf-8", errors="strict")
    except UnicodeDecodeError as error:
        raise ValueError(f"source is not strict UTF-8: {source}") from error
    if "\r" in value.replace("\r\n", ""):
        raise ValueError(f"source contains a bare carriage return: {source}")
    canonical = value.replace("\r\n", "\n")
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def _software_versions() -> dict[str, str]:
    return {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "matplotlib": matplotlib.__version__,
        "msprime": distribution_version("msprime"),
        "stdpopsim": distribution_version("stdpopsim"),
    }


def _stable_seed(base_seed: int, label: str) -> int:
    digest = hashlib.sha256(f"{int(base_seed)}:{label}".encode()).digest()
    return int(1 + int.from_bytes(digest[:8], "big") % (2**31 - 2))


def _atomic_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        temporary.write_text(
            json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_frame(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        if path.suffix == ".gz":
            # Empty gzip member name plus mtime=0 makes bytes independent of the
            # PID-bearing temporary path and therefore reproducible across reruns.
            with temporary.open("wb") as raw:
                with gzip.GzipFile(
                    filename="",
                    mode="wb",
                    compresslevel=6,
                    fileobj=raw,
                    mtime=0,
                ) as compressed:
                    with io.TextIOWrapper(
                        compressed, encoding="utf-8", newline=""
                    ) as text_handle:
                        frame.to_csv(
                            text_handle,
                            sep="\t",
                            index=False,
                            float_format="%.12g",
                            na_rep="NA",
                            lineterminator="\n",
                        )
        else:
            frame.to_csv(
                temporary,
                sep="\t",
                index=False,
                float_format="%.12g",
                na_rep="NA",
                lineterminator="\n",
            )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        temporary.write_text(value, encoding="utf-8", newline="\n")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _integer_sizes(values: np.ndarray) -> np.ndarray:
    sizes = np.floor(np.asarray(values, dtype=float) + 0.5).astype(np.int64)
    if sizes.ndim != 1 or np.any(sizes < 2):
        raise ValueError("diploid Ne schedule must contain positive integer sizes")
    return sizes


def _strict_boolean(series: pd.Series, *, label: str) -> np.ndarray:
    normalized = series.astype(str).str.casefold()
    if not normalized.isin(("true", "false")).all():
        raise ValueError(f"{label} must contain only strict true/false values")
    return normalized.eq("true").to_numpy(dtype=bool)


def build_eas_schedule(repo_root: str | Path) -> TrajectorySchedule:
    """Build the exact integer, pointwise-median piecewise-constant EAS schedule."""

    root = Path(repo_root).resolve()
    resource = root / EAS_RESOURCE_PATH
    artifact = load_phlash_eas_npz(resource, expected_sha256=EAS_RESOURCE_SHA256)
    source = artifact.quantile_trajectories()["median"]
    # Match build_eas_demography_models(): copy the first supported Ne to time zero,
    # then use searchsorted(side="right") for piecewise-constant backward epochs.
    times = np.concatenate(([0.0], artifact.time_generations.astype(float)))
    ne = np.concatenate(([float(source[0])], source.astype(float)))
    ages = np.arange(EAS_ORIGIN_GENERATIONS, -1, -1, dtype=np.int64)
    epoch = np.searchsorted(times, ages.astype(float), side="right") - 1
    sizes = _integer_sizes(ne[np.clip(epoch, 0, len(ne) - 1)])
    schedule_records = [[int(age), int(size)] for age, size in zip(ages, sizes)]
    contract = {
        "source": "PHLASH EAS pointwise median",
        "source_path": EAS_RESOURCE_PATH,
        "source_sha256": sha256_file(resource),
        "trajectory_label": "median",
        "piecewise_constant_lookup": "searchsorted(time, integer_age, side=right)-1",
        "presentward_extrapolation": "copy first supported Ne to generation zero",
        "integer_rounding": "floor(Ne + 0.5)",
        "origin_generations_ago": EAS_ORIGIN_GENERATIONS,
        "origin_diploid_n": int(sizes[0]),
        "present_diploid_n": int(sizes[-1]),
        "integer_schedule_sha256": _canonical_sha256(schedule_records),
        "integer_schedule_rows": len(schedule_records),
    }
    return TrajectorySchedule(
        demography_id="eas_phlash_median",
        ages=ages,
        diploid_sizes=sizes,
        populations=tuple("EAS" for _ in ages),
        initial_mode="de_novo_one_copy",
        initial_frequency=None,
        contract=contract,
    )


def _ancient_eurasia_catalog_contract() -> dict[str, Any]:
    return {
        "catalog_model": "AncientEurasia_9K19",
        "stdpopsim_version": distribution_version("stdpopsim"),
        "validated_pulse_generations_ago": HAN_PULSE_GENERATIONS,
        "validated_pulse_proportion": HAN_PULSE_PROPORTION,
        "validated_loschbour_diploid_n": HAN_LOSCHBOUR_N,
        "validated_han_diploid_n": HAN_PRESENT_N,
        "validated_han_split_generations_ago": HAN_POST_SPLIT_GENERATIONS,
    }


def _validate_ancient_eurasia_contract() -> dict[str, Any]:
    expected = _ancient_eurasia_catalog_contract()
    if expected["stdpopsim_version"] != "0.3.0":
        raise RuntimeError("natural-final-AF requires stdpopsim 0.3.0")
    model = load_ancient_eurasia_model().model
    populations = {population.name: population for population in model.populations}
    population_names = [population.name for population in model.populations]

    def population_name(value: int | str) -> str:
        return value if isinstance(value, str) else population_names[int(value)]

    if not math.isclose(float(populations["Han"].initial_size), HAN_PRESENT_N):
        raise RuntimeError("AncientEurasia Han present size changed")
    pulse = [
        event
        for event in model.events
        if type(event).__name__ == "MassMigration"
        and math.isclose(float(event.time), HAN_PULSE_GENERATIONS)
        and population_name(event.source) == "Loschbour"
        and population_name(event.dest) == "Neanderthal"
    ]
    split = [
        event
        for event in model.events
        if type(event).__name__ == "MassMigration"
        and math.isclose(float(event.time), HAN_POST_SPLIT_GENERATIONS)
        and population_name(event.source) == "Han"
        and population_name(event.dest) == "Loschbour"
    ]
    loschbour_change = [
        event
        for event in model.events
        if type(event).__name__ == "PopulationParametersChange"
        and math.isclose(float(event.time), HAN_POST_SPLIT_GENERATIONS)
        and population_name(event.population) == "Loschbour"
    ]
    if len(pulse) != 1 or not math.isclose(
        float(pulse[0].proportion), HAN_PULSE_PROPORTION
    ):
        raise RuntimeError("AncientEurasia introgression pulse changed")
    if len(split) != 1 or not math.isclose(float(split[0].proportion), 1.0):
        raise RuntimeError("AncientEurasia Han split changed")
    if len(loschbour_change) != 1 or not math.isclose(
        float(loschbour_change[0].initial_size), HAN_LOSCHBOUR_N
    ):
        raise RuntimeError("AncientEurasia ancestral Loschbour size changed")
    return expected


def build_han_schedule(*, validate_catalog: bool = True) -> TrajectorySchedule:
    """Build 256 Loschbour then 2,016 Han Wright--Fisher transitions."""

    catalog = (
        _validate_ancient_eurasia_contract()
        if validate_catalog
        else _ancient_eurasia_catalog_contract()
    )
    ages = np.arange(HAN_PULSE_GENERATIONS, -1, -1, dtype=np.int64)
    # Initial pulse state plus 256 generations at N=2,340; the following 2,016
    # states are Han at N=6,300.  There are exactly 2,272 transitions.
    sizes = np.concatenate(
        (
            np.full(HAN_PRE_SPLIT_GENERATIONS + 1, HAN_LOSCHBOUR_N, dtype=np.int64),
            np.full(HAN_POST_SPLIT_GENERATIONS, HAN_PRESENT_N, dtype=np.int64),
        )
    )
    populations = tuple(
        ["Loschbour"] * (HAN_PRE_SPLIT_GENERATIONS + 1)
        + ["Han"] * HAN_POST_SPLIT_GENERATIONS
    )
    if len(ages) != len(sizes):
        raise RuntimeError("Han trajectory schedule length is inconsistent")
    schedule_records = [
        [int(age), str(population), int(size)]
        for age, population, size in zip(ages, populations, sizes)
    ]
    contract = {
        **catalog,
        "origin": "donor-fixed pulse into recipient Loschbour",
        "pulse_state_convention": "immediately post-pulse at generation 2272",
        "initial_count_distribution": "Binomial(2*2340, 0.0296)",
        "pulse_generations_ago": HAN_PULSE_GENERATIONS,
        "pulse_proportion": HAN_PULSE_PROPORTION,
        "loschbour_transitions": HAN_PRE_SPLIT_GENERATIONS,
        "loschbour_diploid_n": HAN_LOSCHBOUR_N,
        "han_transitions": HAN_POST_SPLIT_GENERATIONS,
        "han_diploid_n": HAN_PRESENT_N,
        "integer_schedule_sha256": _canonical_sha256(schedule_records),
        "integer_schedule_rows": len(schedule_records),
    }
    return TrajectorySchedule(
        demography_id="ancient_eurasia_han_introgression",
        ages=ages,
        diploid_sizes=sizes,
        populations=populations,
        initial_mode="donor_fixed_introgression_pulse",
        initial_frequency=HAN_PULSE_PROPORTION,
        contract=contract,
    )


def selected_gamete_frequency(
    allele_frequency: np.ndarray,
    selection_coefficient: float,
    dominance_coefficient: float = DOMINANCE_COEFFICIENT,
) -> np.ndarray:
    """Return post-viability-selection A frequency for diploid fitness 1,1+hs,1+s."""

    p = np.asarray(allele_frequency, dtype=float)
    s = float(selection_coefficient)
    h = float(dominance_coefficient)
    if np.any((p < 0) | (p > 1)) or s < 0 or not 0 <= h <= 1:
        raise ValueError("invalid allele frequency, s, or h")
    heterozygote = 1.0 + h * s
    homozygote = 1.0 + s
    numerator = p * p * homozygote + p * (1.0 - p) * heterozygote
    denominator = (
        (1.0 - p) ** 2 + 2.0 * p * (1.0 - p) * heterozygote + p * p * homozygote
    )
    return np.divide(
        numerator, denominator, out=np.zeros_like(p), where=denominator > 0
    )


def build_cells(plan: NaturalAfPlan) -> list[dict[str, Any]]:
    plan.validate()
    specifications = (
        (
            "eas_phlash_median",
            "EAS",
            0.0,
            plan.eas_neutral_attempts,
            plan.eas_neutral_batch_size,
        ),
        (
            "eas_phlash_median",
            "EAS",
            0.01,
            plan.eas_selected_attempts,
            plan.eas_selected_batch_size,
        ),
        (
            "ancient_eurasia_han_introgression",
            "Han",
            0.0,
            plan.han_neutral_attempts,
            plan.han_batch_size,
        ),
        (
            "ancient_eurasia_han_introgression",
            "Han",
            0.01,
            plan.han_selected_attempts,
            plan.han_batch_size,
        ),
    )
    cells = []
    for demography, population, coefficient, attempts, batch_size in specifications:
        simulation_class = "neutral" if coefficient == 0 else "selected"
        cell_id = f"{demography}__{simulation_class}_s{coefficient:.3f}".replace(
            ".", "p"
        )
        cells.append(
            {
                "cell_id": cell_id,
                "demography_id": demography,
                "population": population,
                "simulation_class": simulation_class,
                "selection_coefficient": coefficient,
                "dominance_coefficient": DOMINANCE_COEFFICIENT,
                "attempts": int(attempts),
                "batch_size": int(batch_size),
                "n_batches": int(math.ceil(attempts / batch_size)),
            }
        )
    return cells


def _implementation_contract(repo_root: Path) -> dict[str, str]:
    result = {
        "hash_method": (
            "sha256(strict UTF-8 without BOM or bare CR; CRLF canonicalized to LF)"
        )
    }
    for relative in (
        MODULE_PATH,
        WRAPPER_PATH,
        "python/gamma_smc_aou/eas_sweep_models.py",
    ):
        path = repo_root / relative
        if not path.is_file() or path.is_symlink():
            raise FileNotFoundError(path)
        result[relative] = canonical_lf_sha256(path)
    return result


def _build_plan_payload(plan: NaturalAfPlan) -> dict[str, Any]:
    root = plan.repo_root.resolve()
    eas = build_eas_schedule(root)
    han = build_han_schedule(validate_catalog=True)
    return {
        "schema": SCHEMA_VERSION,
        "status": "planned",
        "scientific_contract": {
            "generation_time_years": GENERATION_TIME_YEARS,
            "selection_coefficients": list(SELECTION_COEFFICIENTS),
            "dominance_coefficient": DOMINANCE_COEFFICIENT,
            "wright_fisher_scaling_factor": 1,
            "terminal_population_af_conditioning": False,
            "attempt_denominator": "all fixed attempted trajectories",
            "analysis_filter": "final population alternate count > 0",
            "fixation_counts_as_survival": True,
            "present_panel": {
                "diploids": PANEL_DIPLOIDS,
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
        "paths": {
            "campaign_dir": plan.campaign_dir.resolve().relative_to(root).as_posix(),
            "work_subdir": DEFAULT_WORK_SUBDIR,
            "results_subdir": DEFAULT_RESULTS_SUBDIR,
        },
        "implementation": _implementation_contract(root),
        "software": _software_versions(),
    }


def write_plan(plan: NaturalAfPlan) -> Path:
    payload = _build_plan_payload(plan)
    path = plan.campaign_dir.resolve() / DEFAULT_WORK_SUBDIR / "plan.json"
    if path.is_file():
        existing = json.loads(path.read_text(encoding="utf-8"))
        if _canonical_sha256(existing) != _canonical_sha256(payload):
            raise ValueError("existing natural-final-AF plan is incompatible")
    else:
        _atomic_json(path, payload)
    return path


def _load_plan(plan_path: str | Path) -> dict[str, Any]:
    path = Path(plan_path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema") != SCHEMA_VERSION or payload.get("status") != "planned":
        raise ValueError("natural-final-AF plan schema/status is incompatible")
    return payload


def simulate_trajectory_batch(
    cell: Mapping[str, Any],
    schedule: TrajectorySchedule,
    *,
    attempt_start: int,
    attempt_stop: int,
    batch_index: int,
    batch_seed: int,
) -> pd.DataFrame:
    """Simulate one deterministic batch and retain losses, survivors, and fixations."""

    n = int(attempt_stop) - int(attempt_start)
    if n < 1 or len(schedule.ages) != len(schedule.diploid_sizes):
        raise ValueError("invalid batch bounds or trajectory schedule")
    rng = np.random.default_rng(int(batch_seed))
    initial_n = int(schedule.diploid_sizes[0])
    if schedule.initial_mode == "de_novo_one_copy":
        counts = np.ones(n, dtype=np.int64)
    elif schedule.initial_mode == "donor_fixed_introgression_pulse":
        counts = rng.binomial(
            2 * initial_n, float(schedule.initial_frequency), size=n
        ).astype(np.int64)
    else:
        raise ValueError("unknown initial allele-count mode")
    initial_counts = counts.copy()
    absorption_step = np.full(n, -1, dtype=np.int64)
    absorption_step[counts == 0] = 0
    absorption_step[counts == 2 * initial_n] = 0
    current_n = initial_n
    coefficient = float(cell["selection_coefficient"])
    dominance = float(cell["dominance_coefficient"])
    for step, next_n_raw in enumerate(schedule.diploid_sizes[1:], start=1):
        next_n = int(next_n_raw)
        lost = counts == 0
        fixed = counts == 2 * current_n
        active = ~(lost | fixed)
        next_counts = np.zeros(n, dtype=np.int64)
        next_counts[fixed] = 2 * next_n
        if np.any(active):
            p = counts[active].astype(float) / (2.0 * current_n)
            post_selection = selected_gamete_frequency(p, coefficient, dominance)
            next_counts[active] = rng.binomial(2 * next_n, post_selection)
            newly_absorbed = active & ((next_counts == 0) | (next_counts == 2 * next_n))
            absorption_step[(absorption_step < 0) & newly_absorbed] = step
        counts = next_counts
        current_n = next_n
    final_af = counts.astype(float) / (2.0 * current_n)
    population_survived = counts > 0
    outcome = np.full(n, "segregating", dtype=object)
    outcome[counts == 0] = "lost"
    outcome[counts == 2 * current_n] = "fixed"
    panel_seed = _stable_seed(batch_seed, "present-panel")
    panel_rng = np.random.default_rng(panel_seed)
    sample_alt_count = panel_rng.binomial(2 * PANEL_DIPLOIDS, final_af).astype(np.int64)
    absorption_age = np.full(n, np.nan, dtype=float)
    absorbed = absorption_step >= 0
    absorption_age[absorbed] = (
        schedule.origin_generations_ago - absorption_step[absorbed]
    )
    attempt_indices = np.arange(attempt_start, attempt_stop, dtype=np.int64)
    cell_id = str(cell["cell_id"])
    records = pd.DataFrame(
        {
            "trajectory_id": [
                f"{cell_id}__attempt{value:08d}" for value in attempt_indices
            ],
            "cell_id": cell_id,
            "demography_id": str(cell["demography_id"]),
            "population": str(cell["population"]),
            "simulation_class": str(cell["simulation_class"]),
            "selection_coefficient": coefficient,
            "dominance_coefficient": dominance,
            "attempt_index": attempt_indices,
            "batch_index": int(batch_index),
            "batch_seed": int(batch_seed),
            "panel_seed": int(panel_seed),
            "origin_generations_ago": schedule.origin_generations_ago,
            "initial_population": schedule.populations[0],
            "initial_diploid_n": initial_n,
            "initial_alt_count": initial_counts,
            "initial_population_af": initial_counts.astype(float) / (2.0 * initial_n),
            "present_diploid_n": current_n,
            "final_alt_count": counts,
            "final_population_af": final_af,
            "focal_allele_outcome": outcome,
            "population_survived": population_survived,
            "absorption_generations_ago": absorption_age,
            "panel_diploids": PANEL_DIPLOIDS,
            "sample_alt_count": sample_alt_count,
            "sample_af": sample_alt_count.astype(float) / (2.0 * PANEL_DIPLOIDS),
            "sample_detected": sample_alt_count > 0,
        },
        columns=ATTEMPT_COLUMNS,
    )
    return records


def _batch_paths(work_dir: Path, cell_id: str, batch_index: int) -> tuple[Path, Path]:
    directory = work_dir / "batches" / cell_id
    stem = f"batch_{int(batch_index):05d}"
    return directory / f"{stem}.tsv.gz", directory / f"{stem}.completion.json"


def _batch_contract(
    plan_payload: Mapping[str, Any],
    cell: Mapping[str, Any],
    schedule: TrajectorySchedule,
    *,
    batch_index: int,
    attempt_start: int,
    attempt_stop: int,
    batch_seed: int,
) -> dict[str, Any]:
    cell_id = str(cell["cell_id"])
    output_relative_path = f"batches/{cell_id}/batch_{int(batch_index):05d}.tsv.gz"
    return {
        "schema": BATCH_SCHEMA,
        "plan_sha256": _canonical_sha256(plan_payload),
        "cell": dict(cell),
        "schedule": schedule.contract,
        "trajectory_state": {
            "origin_generations_ago": schedule.origin_generations_ago,
            "initial_population": schedule.populations[0],
            "initial_diploid_n": int(schedule.diploid_sizes[0]),
            "initial_mode": schedule.initial_mode,
            "initial_frequency": schedule.initial_frequency,
            "present_population": schedule.populations[-1],
            "present_diploid_n": int(schedule.diploid_sizes[-1]),
        },
        "batch_index": int(batch_index),
        "attempt_start_inclusive": int(attempt_start),
        "attempt_stop_exclusive": int(attempt_stop),
        "batch_seed": int(batch_seed),
        "panel_seed": _stable_seed(batch_seed, "present-panel"),
        "output_relative_path": output_relative_path,
        "software": _software_versions(),
    }


def _validate_batch(
    table_path: Path, completion_path: Path, expected_contract: Mapping[str, Any]
) -> pd.DataFrame:
    if (
        not table_path.is_file()
        or not completion_path.is_file()
        or table_path.is_symlink()
        or completion_path.is_symlink()
    ):
        raise ValueError("natural-final-AF batch is partial")
    completion = json.loads(completion_path.read_text(encoding="utf-8"))
    contract = completion.get("contract")
    output = completion.get("output", {})
    if (
        completion.get("schema") != BATCH_SCHEMA
        or completion.get("status") != "complete"
        or not isinstance(contract, Mapping)
        or completion.get("contract_sha256") != _canonical_sha256(contract)
        or _canonical_sha256(contract) != _canonical_sha256(expected_contract)
        or output.get("path") != expected_contract["output_relative_path"]
        or output.get("sha256") != sha256_file(table_path)
        or int(output.get("size_bytes", -1)) != table_path.stat().st_size
        or output.get("columns") != list(ATTEMPT_COLUMNS)
    ):
        raise ValueError("natural-final-AF batch completion/checksum is invalid")
    frame = pd.read_csv(table_path, sep="\t")
    if tuple(frame.columns) != ATTEMPT_COLUMNS:
        raise ValueError("natural-final-AF batch columns are invalid")
    start = int(contract["attempt_start_inclusive"])
    stop = int(contract["attempt_stop_exclusive"])
    if len(frame) != stop - start or int(output.get("rows", -1)) != len(frame):
        raise ValueError("natural-final-AF batch row count is invalid")
    if not np.array_equal(
        frame["attempt_index"].to_numpy(dtype=int), np.arange(start, stop)
    ):
        raise ValueError("natural-final-AF batch attempts are not contiguous")

    cell = contract["cell"]
    state = contract["trajectory_state"]
    static_values = {
        "cell_id": str(cell["cell_id"]),
        "demography_id": str(cell["demography_id"]),
        "population": str(cell["population"]),
        "simulation_class": str(cell["simulation_class"]),
        "selection_coefficient": float(cell["selection_coefficient"]),
        "dominance_coefficient": float(cell["dominance_coefficient"]),
        "batch_index": int(contract["batch_index"]),
        "batch_seed": int(contract["batch_seed"]),
        "panel_seed": int(contract["panel_seed"]),
        "origin_generations_ago": int(state["origin_generations_ago"]),
        "initial_population": str(state["initial_population"]),
        "initial_diploid_n": int(state["initial_diploid_n"]),
        "present_diploid_n": int(state["present_diploid_n"]),
        "panel_diploids": PANEL_DIPLOIDS,
    }
    for column, expected in static_values.items():
        values = frame[column].to_numpy()
        if isinstance(expected, float):
            valid = np.isfinite(values.astype(float)) & np.isclose(
                values.astype(float), expected, rtol=0, atol=1e-15
            )
        else:
            valid = values.astype(str) == str(expected)
        if not np.all(valid):
            raise ValueError(f"natural-final-AF batch {column} contract is invalid")

    expected_ids = np.asarray(
        [f"{cell['cell_id']}__attempt{index:08d}" for index in range(start, stop)]
    )
    if not np.array_equal(frame["trajectory_id"].astype(str).to_numpy(), expected_ids):
        raise ValueError("natural-final-AF trajectory IDs are invalid")

    integer_columns = (
        "initial_alt_count",
        "final_alt_count",
        "sample_alt_count",
    )
    integers: dict[str, np.ndarray] = {}
    for column in integer_columns:
        values = pd.to_numeric(frame[column], errors="coerce").to_numpy(dtype=float)
        if not np.all(np.isfinite(values) & (values == np.floor(values))):
            raise ValueError(f"natural-final-AF {column} is not finite integer data")
        integers[column] = values.astype(np.int64)

    initial_count = integers["initial_alt_count"]
    initial_n = int(state["initial_diploid_n"])
    if np.any((initial_count < 0) | (initial_count > 2 * initial_n)):
        raise ValueError("natural-final-AF initial counts are out of range")
    if state["initial_mode"] == "de_novo_one_copy" and not np.all(initial_count == 1):
        raise ValueError("natural-final-AF de novo initial count must be one")
    initial_af = pd.to_numeric(
        frame["initial_population_af"], errors="coerce"
    ).to_numpy(dtype=float)
    if not np.all(
        np.isfinite(initial_af)
        & np.isclose(initial_af, initial_count / (2 * initial_n), rtol=0, atol=5e-12)
    ):
        raise ValueError("natural-final-AF initial AF semantics are invalid")

    final_count = integers["final_alt_count"]
    present_n = int(state["present_diploid_n"])
    denominator = 2 * present_n
    final_af = pd.to_numeric(frame["final_population_af"], errors="coerce").to_numpy(
        dtype=float
    )
    if np.any((final_count < 0) | (final_count > denominator)) or not np.all(
        np.isfinite(final_af)
        & (final_af >= 0)
        & (final_af <= 1)
        & np.isclose(final_af, final_count / denominator, rtol=0, atol=5e-12)
    ):
        raise ValueError("natural-final-AF terminal count/AF semantics are invalid")
    expected_fate = np.full(len(frame), "segregating", dtype=object)
    expected_fate[final_count == 0] = "lost"
    expected_fate[final_count == denominator] = "fixed"
    if not np.array_equal(
        frame["focal_allele_outcome"].astype(str).to_numpy(), expected_fate
    ):
        raise ValueError("natural-final-AF focal allele fate is invalid")
    survived = _strict_boolean(
        frame["population_survived"], label="population_survived"
    )
    if not np.array_equal(survived, final_count > 0):
        raise ValueError("natural-final-AF population survival semantics are invalid")

    sample_count = integers["sample_alt_count"]
    sample_af = pd.to_numeric(frame["sample_af"], errors="coerce").to_numpy(dtype=float)
    if np.any((sample_count < 0) | (sample_count > 2 * PANEL_DIPLOIDS)) or not np.all(
        np.isfinite(sample_af)
        & (sample_af >= 0)
        & (sample_af <= 1)
        & np.isclose(
            sample_af,
            sample_count / (2 * PANEL_DIPLOIDS),
            rtol=0,
            atol=5e-12,
        )
    ):
        raise ValueError("natural-final-AF panel count/AF semantics are invalid")
    detected = _strict_boolean(frame["sample_detected"], label="sample_detected")
    if not np.array_equal(detected, sample_count > 0):
        raise ValueError("natural-final-AF sample detection semantics are invalid")

    absorption_age = pd.to_numeric(
        frame["absorption_generations_ago"], errors="coerce"
    ).to_numpy(dtype=float)
    absorbed = expected_fate != "segregating"
    if (
        np.any(~np.isnan(absorption_age[~absorbed]))
        or np.any(~np.isfinite(absorption_age[absorbed]))
        or np.any(absorption_age[absorbed] < 0)
        or np.any(absorption_age[absorbed] > int(state["origin_generations_ago"]))
    ):
        raise ValueError("natural-final-AF absorption timing is invalid")
    return frame


def _schedule_for_cell(repo_root: Path, cell: Mapping[str, Any]) -> TrajectorySchedule:
    if str(cell["demography_id"]) == "eas_phlash_median":
        return build_eas_schedule(repo_root)
    if str(cell["demography_id"]) == "ancient_eurasia_han_introgression":
        return build_han_schedule(validate_catalog=False)
    raise ValueError(f"unknown demography: {cell['demography_id']}")


def _simulate_batch_task(task: Mapping[str, Any]) -> dict[str, Any]:
    root = Path(str(task["repo_root"]))
    work_dir = Path(str(task["work_dir"]))
    plan_payload = dict(task["plan_payload"])
    cell = dict(task["cell"])
    batch_index = int(task["batch_index"])
    start = batch_index * int(cell["batch_size"])
    stop = min(int(cell["attempts"]), start + int(cell["batch_size"]))
    seed = _stable_seed(
        int(plan_payload["rng"]["base_seed"]), f"{cell['cell_id']}:batch:{batch_index}"
    )
    schedule = _schedule_for_cell(root, cell)
    contract = _batch_contract(
        plan_payload,
        cell,
        schedule,
        batch_index=batch_index,
        attempt_start=start,
        attempt_stop=stop,
        batch_seed=seed,
    )
    table_path, completion_path = _batch_paths(
        work_dir, str(cell["cell_id"]), batch_index
    )
    if completion_path.exists() and not table_path.exists():
        raise ValueError("natural-final-AF completion exists without its batch table")
    if completion_path.exists():
        _validate_batch(table_path, completion_path, contract)
        return {
            "cell_id": cell["cell_id"],
            "batch_index": batch_index,
            "status": "cached",
        }
    # A table without a completion marker is an interrupted, unpublished write.
    # Regenerate deterministically and atomically replace only that owned table.
    frame = simulate_trajectory_batch(
        cell,
        schedule,
        attempt_start=start,
        attempt_stop=stop,
        batch_index=batch_index,
        batch_seed=seed,
    )
    _atomic_frame(table_path, frame)
    completion = {
        "schema": BATCH_SCHEMA,
        "status": "complete",
        "contract": contract,
        "contract_sha256": _canonical_sha256(contract),
        "output": {
            "path": contract["output_relative_path"],
            "sha256": sha256_file(table_path),
            "size_bytes": table_path.stat().st_size,
            "rows": len(frame),
            "columns": list(frame.columns),
        },
    }
    _atomic_json(completion_path, completion)
    _validate_batch(table_path, completion_path, contract)
    return {
        "cell_id": cell["cell_id"],
        "batch_index": batch_index,
        "status": "complete",
    }


def simulate_study(plan: NaturalAfPlan, *, workers: int = 4) -> pd.DataFrame:
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
        rows = [_simulate_batch_task(task) for task in tasks]
    else:
        rows = []
        with ProcessPoolExecutor(max_workers=int(workers)) as pool:
            futures = {pool.submit(_simulate_batch_task, task): task for task in tasks}
            for future in as_completed(futures):
                rows.append(future.result())
    status = pd.DataFrame(rows).sort_values(["cell_id", "batch_index"])
    _atomic_frame(work_dir / "last_simulation_status.tsv", status)
    return status.reset_index(drop=True)


def _expected_batch_inputs(
    plan: NaturalAfPlan, payload: Mapping[str, Any]
) -> tuple[list[pd.DataFrame], list[dict[str, Any]]]:
    work_dir = plan.campaign_dir.resolve() / DEFAULT_WORK_SUBDIR
    frames: list[pd.DataFrame] = []
    records: list[dict[str, Any]] = []
    for cell in payload["cells"]:
        schedule = _schedule_for_cell(plan.repo_root.resolve(), cell)
        for batch_index in range(int(cell["n_batches"])):
            start = batch_index * int(cell["batch_size"])
            stop = min(int(cell["attempts"]), start + int(cell["batch_size"]))
            seed = _stable_seed(
                int(payload["rng"]["base_seed"]),
                f"{cell['cell_id']}:batch:{batch_index}",
            )
            contract = _batch_contract(
                payload,
                cell,
                schedule,
                batch_index=batch_index,
                attempt_start=start,
                attempt_stop=stop,
                batch_seed=seed,
            )
            table, completion = _batch_paths(
                work_dir, str(cell["cell_id"]), batch_index
            )
            frames.append(_validate_batch(table, completion, contract))
            records.append(
                {
                    "cell_id": cell["cell_id"],
                    "batch_index": batch_index,
                    "table_sha256": sha256_file(table),
                    "completion_sha256": sha256_file(completion),
                }
            )
    return frames, records


def _wilson_interval(successes: int, total: int) -> tuple[float, float]:
    if total < 1:
        return np.nan, np.nan
    z = 1.959963984540054
    proportion = successes / total
    denominator = 1.0 + z * z / total
    center = (proportion + z * z / (2.0 * total)) / denominator
    half_width = (
        z
        * math.sqrt(
            proportion * (1.0 - proportion) / total + z * z / (4.0 * total * total)
        )
        / denominator
    )
    return center - half_width, center + half_width


def summarize_attempts(attempts: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    population_survived = _strict_boolean(
        attempts["population_survived"], label="population_survived"
    )
    sample_detected = _strict_boolean(
        attempts["sample_detected"], label="sample_detected"
    )
    working = attempts.copy()
    working["_population_survived"] = population_survived
    working["_sample_detected"] = sample_detected
    summaries = []
    for keys, group in working.groupby(
        ["demography_id", "population", "simulation_class", "selection_coefficient"],
        sort=True,
    ):
        survived = group[group["_population_survived"]]
        detected = survived[survived["_sample_detected"]]
        af = survived["final_population_af"].to_numpy(dtype=float)
        sample_af = survived["sample_af"].to_numpy(dtype=float)
        initial_af_all = group["initial_population_af"].to_numpy(dtype=float)
        final_af_all = group["final_population_af"].to_numpy(dtype=float)
        fate_counts = group["focal_allele_outcome"].value_counts()
        n_lost = int(fate_counts.get("lost", 0))
        n_segregating = int(fate_counts.get("segregating", 0))
        n_fixed = int(fate_counts.get("fixed", 0))
        if n_lost + n_segregating + n_fixed != len(group):
            raise ValueError("unexpected focal allele fate in attempts")
        ci_low, ci_high = _wilson_interval(len(survived), len(group))
        initial_mean = float(np.mean(initial_af_all))
        final_mean = float(np.mean(final_af_all))
        paired_difference = final_af_all - initial_af_all
        final_mc_se = float(np.std(paired_difference, ddof=1) / math.sqrt(len(group)))
        if str(keys[2]) == "neutral":
            martingale_difference = final_mean - initial_mean
            martingale_z = (
                martingale_difference / final_mc_se if final_mc_se > 0 else np.nan
            )
        else:
            martingale_difference = np.nan
            final_mc_se = np.nan
            martingale_z = np.nan
        summaries.append(
            {
                "demography_id": keys[0],
                "population": keys[1],
                "simulation_class": keys[2],
                "selection_coefficient": keys[3],
                "n_attempted": len(group),
                "n_lost": n_lost,
                "n_segregating": n_segregating,
                "n_fixed": n_fixed,
                "lost_fraction": n_lost / len(group),
                "segregating_fraction": n_segregating / len(group),
                "fixed_fraction": n_fixed / len(group),
                "n_population_survived": len(survived),
                "population_survival_fraction": len(survived) / len(group),
                "population_survival_wilson_ci_low": ci_low,
                "population_survival_wilson_ci_high": ci_high,
                "fixation_fraction_among_survivors": (
                    n_fixed / len(survived) if len(survived) else np.nan
                ),
                "initial_population_af_mean_all_attempts": initial_mean,
                "final_population_af_mean_all_attempts": final_mean,
                "neutral_martingale_difference": martingale_difference,
                "neutral_martingale_mc_se": final_mc_se,
                "neutral_martingale_z": martingale_z,
                "n_sample_detected_among_survivors": len(detected),
                "sample_detection_fraction_among_survivors": (
                    len(detected) / len(survived) if len(survived) else np.nan
                ),
                "population_af_mean_survivors": np.mean(af) if len(af) else np.nan,
                "population_af_median_survivors": np.median(af) if len(af) else np.nan,
                "population_af_q025_survivors": np.quantile(af, 0.025)
                if len(af)
                else np.nan,
                "population_af_q25_survivors": np.quantile(af, 0.25)
                if len(af)
                else np.nan,
                "population_af_q75_survivors": np.quantile(af, 0.75)
                if len(af)
                else np.nan,
                "population_af_q975_survivors": np.quantile(af, 0.975)
                if len(af)
                else np.nan,
                "sample_af_mean_survivors": np.mean(sample_af)
                if len(sample_af)
                else np.nan,
                "sample_af_median_survivors": np.median(sample_af)
                if len(sample_af)
                else np.nan,
            }
        )
    summary = pd.DataFrame(summaries, columns=SUMMARY_COLUMNS)
    comparisons = []
    for demography, group in working[working["_population_survived"]].groupby(
        "demography_id"
    ):
        neutral = np.sort(
            group.loc[
                group["simulation_class"] == "neutral", "final_population_af"
            ].to_numpy(float)
        )
        selected = group.loc[
            group["simulation_class"] == "selected", "final_population_af"
        ].to_numpy(float)
        if len(neutral) and len(selected):
            less = np.searchsorted(neutral, selected, side="left")
            equal = np.searchsorted(neutral, selected, side="right") - less
            auc = float(np.mean((less + 0.5 * equal) / len(neutral)))
        else:
            auc = np.nan
        comparisons.append(
            {
                "demography_id": demography,
                "statistic": "P(AF_s0.01 > AF_s0 | population survival) + half ties",
                "population_af_roc_auc_survivors": auc,
                "n_neutral_survivors": len(neutral),
                "n_selected_survivors": len(selected),
            }
        )
    return summary, pd.DataFrame(comparisons, columns=COMPARISON_COLUMNS)


def _figure_geometry(path: Path) -> dict[str, float | int]:
    if path.suffix == ".png":
        header = path.read_bytes()[:24]
        if len(header) != 24 or header[:8] != b"\x89PNG\r\n\x1a\n":
            raise ValueError(f"invalid PNG output: {path}")
        width, height = struct.unpack(">II", header[16:24])
        return {"width_px": int(width), "height_px": int(height)}
    if path.suffix == ".pdf":
        payload = path.read_bytes()
        match = re.search(
            rb"/MediaBox\s*\[\s*0(?:\.0+)?\s+0(?:\.0+)?\s+"
            rb"([0-9.]+)\s+([0-9.]+)\s*\]",
            payload,
        )
        if not payload.startswith(b"%PDF-") or match is None:
            raise ValueError(f"invalid PDF output: {path}")
        return {
            "width_points": float(match.group(1)),
            "height_points": float(match.group(2)),
            "page_count": len(re.findall(rb"/Type\s*/Page\b", payload)),
        }
    raise ValueError(f"unsupported figure type: {path}")


def _validate_letter_geometry(path: Path) -> dict[str, float | int]:
    geometry = _figure_geometry(path)
    if path.suffix == ".png":
        expected = {
            "width_px": int(FIGURE_SIZE_INCHES[0] * PNG_DPI),
            "height_px": int(FIGURE_SIZE_INCHES[1] * PNG_DPI),
        }
        if geometry != expected:
            raise ValueError(f"PNG is not 11 x 8.5 inches at {PNG_DPI} dpi: {path}")
    else:
        if not (
            math.isclose(float(geometry["width_points"]), 792.0, abs_tol=0.1)
            and math.isclose(float(geometry["height_points"]), 612.0, abs_tol=0.1)
            and int(geometry["page_count"]) == 1
        ):
            raise ValueError(f"PDF is not one-page landscape US letter: {path}")
    return geometry


def _save_figure_pair(figure: plt.Figure, result_dir: Path, stem: str) -> list[Path]:
    paths = []
    for extension in ("png", "pdf"):
        path = result_dir / f"{stem}.{extension}"
        temporary = path.with_name(f".{path.stem}.tmp.{os.getpid()}{path.suffix}")
        figure.savefig(
            temporary,
            dpi=PNG_DPI if extension == "png" else None,
            metadata={"Creator": "gamma_smc_aou.natural_final_af"},
        )
        _validate_letter_geometry(temporary)
        os.replace(temporary, path)
        paths.append(path)
    return paths


def _plot_unconditional_fates(summary: pd.DataFrame, result_dir: Path) -> list[Path]:
    labels = {
        "eas_phlash_median": "EAS age-matched de novo (50 kya)",
        "ancient_eurasia_han_introgression": "Han donor-fixed pulse (~2.96%)",
    }
    fate_columns = (
        ("lost_fraction", "Lost", "#B9B9B9"),
        ("segregating_fraction", "Segregating", "#4C78A8"),
        ("fixed_fraction", "Fixed", "#D1495B"),
    )
    figure, axes = plt.subplots(1, 2, figsize=FIGURE_SIZE_INCHES, sharey=True)
    for axis, demography in zip(axes, labels, strict=True):
        cell = summary[summary["demography_id"] == demography].set_index(
            "simulation_class"
        )
        x = np.arange(2)
        bottom = np.zeros(2)
        for column, _, color in fate_columns:
            values = np.asarray(
                [
                    cell.loc[simulation_class, column]
                    for simulation_class in ("neutral", "selected")
                ]
            )
            axis.bar(x, values, bottom=bottom, width=0.64, color=color)
            bottom += values
        axis.set_xticks(x, ("s=0", "s=0.01"), fontsize=15)
        axis.set_ylim(0, 1.0)
        axis.set_title(labels[demography], fontsize=18)
        axis.tick_params(axis="y", labelsize=13)
        axis.grid(axis="y", alpha=0.2)
        n_text = "   ".join(
            f"{simulation_class}: n={int(cell.loc[simulation_class, 'n_attempted']):,}"
            for simulation_class in ("neutral", "selected")
        )
        axis.text(
            0.5, -0.13, n_text, transform=axis.transAxes, ha="center", fontsize=12
        )
        for x_value, simulation_class in enumerate(("neutral", "selected")):
            survival = float(cell.loc[simulation_class, "population_survival_fraction"])
            fixed_given_survival = float(
                cell.loc[simulation_class, "fixation_fraction_among_survivors"]
            )
            axis.text(
                x_value,
                0.035,
                f"survive {100 * survival:.3g}%\n"
                f"fixed | survive {100 * fixed_given_survival:.3g}%",
                ha="center",
                va="bottom",
                fontsize=11,
            )
    axes[0].set_ylabel("Fraction of all attempted trajectories", fontsize=16)
    handles = [Patch(facecolor=color, label=label) for _, label, color in fate_columns]
    figure.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.91),
        ncol=3,
        fontsize=14,
        frameon=False,
    )
    figure.suptitle("Unconditional focal-allele fates", fontsize=20, y=0.98)
    figure.subplots_adjust(top=0.84, bottom=0.19, left=0.10, right=0.98, wspace=0.16)
    paths = _save_figure_pair(figure, result_dir, "unconditional_fate_fractions")
    plt.close(figure)
    return paths


def _plot_ecdf(attempts: pd.DataFrame, result_dir: Path) -> list[Path]:
    colors = {"neutral": "#4C78A8", "selected": "#D1495B"}
    labels = {
        "eas_phlash_median": "EAS age-matched de novo (50 kya)",
        "ancient_eurasia_han_introgression": "Han donor-fixed pulse (~2.96%)",
    }
    figure, axes = plt.subplots(1, 2, figsize=FIGURE_SIZE_INCHES, sharey=True)
    survived = _strict_boolean(
        attempts["population_survived"], label="population_survived"
    )
    survivors = attempts[survived]
    for axis, demography in zip(axes, labels, strict=True):
        panel_counts = {}
        panel_fixation = {}
        for simulation_class in ("neutral", "selected"):
            values = np.sort(
                survivors.loc[
                    (survivors["demography_id"] == demography)
                    & (survivors["simulation_class"] == simulation_class),
                    "final_population_af",
                ].to_numpy(float)
            )
            panel_counts[simulation_class] = len(values)
            panel_fixation[simulation_class] = (
                float(np.mean(values == 1.0)) if len(values) else np.nan
            )
            if len(values):
                ecdf_x = np.append(values, 1.0)
                ecdf_y = np.append(np.arange(1, len(values) + 1) / len(values), 1.0)
                axis.step(
                    ecdf_x,
                    ecdf_y,
                    where="post",
                    lw=3,
                    color=colors[simulation_class],
                )
        axis.set_xscale("symlog", linthresh=1e-4, base=10)
        axis.set_xlim(0, 1)
        axis.set_ylim(0, 1.01)
        axis.set_title(labels[demography], fontsize=18)
        axis.tick_params(labelsize=13)
        axis.grid(alpha=0.2)
        axis.text(
            0.03,
            0.96,
            f"survivors: s=0 n={panel_counts['neutral']:,}\n"
            f"s=0.01 n={panel_counts['selected']:,}\n"
            f"fixed | survive: {100 * panel_fixation['neutral']:.2f}% / "
            f"{100 * panel_fixation['selected']:.2f}%",
            transform=axis.transAxes,
            va="top",
            fontsize=12,
        )
    axes[0].set_ylabel("Empirical CDF among population survivors", fontsize=16)
    handles = [
        Line2D([0], [0], color=colors["neutral"], lw=3, label="s=0"),
        Line2D([0], [0], color=colors["selected"], lw=3, label="s=0.01"),
    ]
    figure.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.86),
        ncol=2,
        fontsize=14,
        frameon=False,
    )
    figure.suptitle(
        "Natural final-AF distributions\n(no terminal-frequency conditioning)",
        fontsize=20,
        y=0.98,
    )
    figure.supxlabel("Present population allele frequency", fontsize=16, y=0.045)
    figure.subplots_adjust(top=0.78, bottom=0.15, left=0.10, right=0.98, wspace=0.16)
    paths = _save_figure_pair(figure, result_dir, "final_population_af_survivor_ecdf")
    plt.close(figure)
    return paths


def _plot_population_vs_sample(attempts: pd.DataFrame, result_dir: Path) -> list[Path]:
    colors = {"neutral": "#4C78A8", "selected": "#D1495B"}
    labels = {
        "eas_phlash_median": "EAS age-matched de novo (50 kya)",
        "ancient_eurasia_han_introgression": "Han donor-fixed pulse (~2.96%)",
    }
    figure, axes = plt.subplots(1, 2, figsize=FIGURE_SIZE_INCHES, sharey=True)
    survived = _strict_boolean(
        attempts["population_survived"], label="population_survived"
    )
    survivors = attempts[survived]
    for axis, demography in zip(axes, labels, strict=True):
        panel_counts = {}
        for simulation_class in ("neutral", "selected"):
            cell = survivors[
                (survivors["demography_id"] == demography)
                & (survivors["simulation_class"] == simulation_class)
            ]
            panel_counts[simulation_class] = len(cell)
            for column, linestyle in (
                ("final_population_af", "-"),
                ("sample_af", "--"),
            ):
                values = np.sort(cell[column].to_numpy(float))
                if len(values):
                    ecdf_x = np.append(values, 1.0)
                    ecdf_y = np.append(np.arange(1, len(values) + 1) / len(values), 1.0)
                    axis.step(
                        ecdf_x,
                        ecdf_y,
                        where="post",
                        lw=2.7,
                        ls=linestyle,
                        color=colors[simulation_class],
                    )
        axis.set_xscale("symlog", linthresh=1e-4, base=10)
        axis.set_xlim(0, 1)
        axis.set_ylim(0, 1.01)
        axis.set_title(labels[demography], fontsize=18)
        axis.tick_params(labelsize=13)
        axis.grid(alpha=0.2)
        axis.text(
            0.03,
            0.96,
            f"survivors: s=0 n={panel_counts['neutral']:,}\n"
            f"s=0.01 n={panel_counts['selected']:,}",
            transform=axis.transAxes,
            va="top",
            fontsize=12,
        )
    axes[0].set_ylabel("Empirical CDF among population survivors", fontsize=16)
    handles = [
        Line2D([0], [0], color=colors["neutral"], lw=2.7, label="s=0"),
        Line2D([0], [0], color=colors["selected"], lw=2.7, label="s=0.01"),
        Line2D([0], [0], color="black", lw=2.7, ls="-", label="population"),
        Line2D(
            [0],
            [0],
            color="black",
            lw=2.7,
            ls="--",
            label="100-diploid panel",
        ),
    ]
    figure.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.91),
        ncol=4,
        fontsize=12,
        frameon=False,
    )
    figure.suptitle(
        "Population AF versus one unconditioned 100-diploid panel", fontsize=20, y=0.98
    )
    figure.supxlabel("Present allele frequency", fontsize=16, y=0.045)
    figure.subplots_adjust(top=0.82, bottom=0.15, left=0.10, right=0.98, wspace=0.16)
    paths = _save_figure_pair(
        figure, result_dir, "population_vs_sample_af_survivor_ecdf"
    )
    plt.close(figure)
    return paths


def _schedule_table(repo_root: Path) -> pd.DataFrame:
    rows = []
    for schedule in (
        build_eas_schedule(repo_root),
        build_han_schedule(validate_catalog=False),
    ):
        rows.extend(
            {
                "demography_id": schedule.demography_id,
                "generations_ago": int(age),
                "population": population,
                "diploid_n": int(size),
            }
            for age, population, size in zip(
                schedule.ages, schedule.populations, schedule.diploid_sizes
            )
        )
    return pd.DataFrame(rows, columns=SCHEDULE_COLUMNS)


def _tsv_shape(path: Path) -> tuple[int, list[str]]:
    if path.suffix == ".gz":
        handle = gzip.open(path, mode="rt", encoding="utf-8", newline="")
    else:
        handle = path.open(mode="rt", encoding="utf-8", newline="")
    with handle:
        header = handle.readline().rstrip("\r\n")
        rows = sum(1 for _ in handle)
    return rows, header.split("\t") if header else []


def _analysis_table_columns() -> dict[str, tuple[str, ...]]:
    return {
        "attempts": ATTEMPT_COLUMNS,
        "survivors": ATTEMPT_COLUMNS,
        "summary": SUMMARY_COLUMNS,
        "comparisons": COMPARISON_COLUMNS,
        "schedules": SCHEDULE_COLUMNS,
    }


class _HashingTextSink:
    """Minimal text sink for streaming a canonical TSV into SHA-256."""

    def __init__(self) -> None:
        self.digest = hashlib.sha256()

    def write(self, value: str) -> int:
        self.digest.update(value.encode("utf-8"))
        return len(value)

    def flush(self) -> None:
        return None


def _frame_tsv_sha256(frame: pd.DataFrame) -> str:
    sink = _HashingTextSink()
    frame.to_csv(
        sink,
        sep="\t",
        index=False,
        float_format="%.12g",
        na_rep="NA",
        lineterminator="\n",
        chunksize=50_000,
    )
    return sink.digest.hexdigest()


def _tsv_payload_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    if path.suffix == ".gz":
        handle = gzip.open(path, mode="rb")
    else:
        handle = path.open(mode="rb")
    with handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _describe_analysis_output(
    label: str,
    path: Path,
    *,
    table: pd.DataFrame | None = None,
) -> dict[str, Any]:
    record: dict[str, Any] = {
        "path": path.name,
        "sha256": sha256_file(path),
        "size_bytes": path.stat().st_size,
    }
    if label in _analysis_table_columns():
        if table is None:
            rows, columns = _tsv_shape(path)
        else:
            rows, columns = len(table), list(table.columns)
        record.update({"rows": int(rows), "columns": list(columns)})
    elif path.suffix in {".png", ".pdf"}:
        record["geometry"] = _validate_letter_geometry(path)
    return record


def _validate_analysis_bundle(
    result_dir: Path,
    contract: Mapping[str, Any],
) -> dict[str, Any]:
    completion_path = result_dir / "analysis_completion.json"
    if not result_dir.is_dir() or not completion_path.is_file():
        raise ValueError("natural-final-AF result bundle is incomplete")
    completion = json.loads(completion_path.read_text(encoding="utf-8"))
    outputs = completion.get("outputs")
    expected_files = set(ANALYSIS_OUTPUT_NAMES.values()) | {completion_path.name}
    actual_entries = {entry.name for entry in result_dir.iterdir()}
    if (
        completion.get("schema") != ANALYSIS_SCHEMA
        or completion.get("status") != "complete"
        or not isinstance(completion.get("contract"), Mapping)
        or completion.get("contract_sha256")
        != _canonical_sha256(completion["contract"])
        or _canonical_sha256(completion["contract"]) != _canonical_sha256(contract)
        or not isinstance(outputs, Mapping)
        or set(outputs) != set(ANALYSIS_OUTPUT_NAMES)
        or actual_entries != expected_files
    ):
        raise ValueError("natural-final-AF analysis completion contract is invalid")

    table_shapes: dict[str, tuple[int, list[str]]] = {}
    for label, expected_name in ANALYSIS_OUTPUT_NAMES.items():
        output = outputs[label]
        path = result_dir / expected_name
        if (
            not isinstance(output, Mapping)
            or output.get("path") != expected_name
            or not path.is_file()
            or path.is_symlink()
            or output.get("sha256") != sha256_file(path)
            or int(output.get("size_bytes", -1)) != path.stat().st_size
        ):
            raise ValueError(f"cached natural-final-AF output failed: {label}")
        if label in _analysis_table_columns():
            rows, columns = _tsv_shape(path)
            table_shapes[label] = (rows, columns)
            if (
                int(output.get("rows", -1)) != rows
                or output.get("columns") != columns
                or columns != list(_analysis_table_columns()[label])
                or _tsv_payload_sha256(path) != contract["table_payload_sha256"][label]
            ):
                raise ValueError(f"cached natural-final-AF table shape failed: {label}")
        elif path.suffix in {".png", ".pdf"}:
            geometry = _validate_letter_geometry(path)
            if output.get("geometry") != geometry:
                raise ValueError(
                    f"cached natural-final-AF figure geometry failed: {label}"
                )

    expected_attempts = sum(int(cell["attempts"]) for cell in contract["cells"])
    expected_schedules = EAS_ORIGIN_GENERATIONS + 1 + HAN_PULSE_GENERATIONS + 1
    if (
        table_shapes["attempts"][0] != expected_attempts
        or table_shapes["summary"][0] != len(contract["cells"])
        or table_shapes["comparisons"][0] != 2
        or table_shapes["schedules"][0] != expected_schedules
    ):
        raise ValueError("natural-final-AF result table row contract failed")
    summary = pd.read_csv(result_dir / ANALYSIS_OUTPUT_NAMES["summary"], sep="\t")
    if (
        int(summary["n_attempted"].sum()) != table_shapes["attempts"][0]
        or int(summary["n_population_survived"].sum()) != table_shapes["survivors"][0]
        or not np.array_equal(
            summary["n_lost"].to_numpy(dtype=int)
            + summary["n_segregating"].to_numpy(dtype=int)
            + summary["n_fixed"].to_numpy(dtype=int),
            summary["n_attempted"].to_numpy(dtype=int),
        )
        or not np.array_equal(
            summary["n_segregating"].to_numpy(dtype=int)
            + summary["n_fixed"].to_numpy(dtype=int),
            summary["n_population_survived"].to_numpy(dtype=int),
        )
    ):
        raise ValueError("natural-final-AF summary/output row totals disagree")
    return dict(completion)


def _remove_owned_staging(stage: Path, parent: Path) -> None:
    resolved_stage = stage.resolve()
    resolved_parent = parent.resolve()
    if resolved_stage.parent != resolved_parent or not resolved_stage.name.startswith(
        ".natural_final_af.staging."
    ):
        raise RuntimeError("refusing to remove an unowned analysis staging path")
    if resolved_stage.exists():
        shutil.rmtree(resolved_stage)


def analyze_study(plan: NaturalAfPlan) -> dict[str, Any]:
    plan_path = write_plan(plan)
    payload = _load_plan(plan_path)
    frames, batch_inputs = _expected_batch_inputs(plan, payload)
    attempts = pd.concat(frames, ignore_index=True).sort_values(
        ["cell_id", "attempt_index"]
    )
    survived = _strict_boolean(
        attempts["population_survived"], label="population_survived"
    )
    survivors = attempts[survived].copy()
    summary, comparisons = summarize_attempts(attempts)
    schedules = _schedule_table(plan.repo_root.resolve())
    tables = {
        "attempts": attempts,
        "survivors": survivors,
        "summary": summary,
        "comparisons": comparisons,
        "schedules": schedules,
    }
    result_dir = plan.campaign_dir.resolve() / DEFAULT_RESULTS_SUBDIR
    contract = {
        "schema": ANALYSIS_SCHEMA,
        "plan_sha256": sha256_file(plan_path),
        "batch_inputs": batch_inputs,
        "cells": payload["cells"],
        "survival_filter": "final_alt_count > 0",
        "sample_panel_is_not_a_filter": True,
        "implementation": _implementation_contract(plan.repo_root.resolve()),
        "software": _software_versions(),
        "table_payload_sha256": {
            label: _frame_tsv_sha256(frame) for label, frame in tables.items()
        },
    }
    if result_dir.exists():
        if not (result_dir / "analysis_completion.json").is_file():
            raise ValueError(
                "incompatible partial natural-final-AF result bundle exists"
            )
        completion = _validate_analysis_bundle(result_dir, contract)
        return {**completion, "cache_hit": True}
    result_parent = result_dir.parent
    result_parent.mkdir(parents=True, exist_ok=True)
    stage = Path(
        tempfile.mkdtemp(prefix=".natural_final_af.staging.", dir=result_parent)
    )
    try:
        for label, frame in tables.items():
            _atomic_frame(stage / ANALYSIS_OUTPUT_NAMES[label], frame)
        _plot_unconditional_fates(summary, stage)
        _plot_ecdf(attempts, stage)
        _plot_population_vs_sample(attempts, stage)
        _atomic_text(
            stage / ANALYSIS_OUTPUT_NAMES["readme"],
            "# Natural final-AF sensitivity\n\n"
            "All fixed attempted Wright--Fisher trajectories are retained. Population "
            "survival (`final_alt_count > 0`) is the only analysis filter; there is no "
            "terminal-AF gate, and fixation counts as survival. Each trajectory also has "
            "one independent Binomial(200, population AF) panel draw with no resampling.\n\n"
            "The EAS variant begins as one copy 2,000 generations (50 kya) ago under the "
            "PHLASH pointwise-median piecewise-constant Ne history; this is an "
            "age-matched de novo neutral comparison, not a mixture over neutral mutation "
            "ages. The Han donor-fixed pulse approximation sets K0 to the immediately-"
            "post-pulse Loschbour state Binomial(4,680, 0.0296) at generation "
            "2,272, followed by 256 Loschbour and 2,016 Han transitions. Selection is "
            "additive (`h=0.5`) and runs from origin/pulse through the present for the "
            "`s=0.01` class. Survivor ECDFs estimate P(final AF | K_present > 0); "
            "unconditional fate fractions preserve the loss mass and fixation atom. The "
            "all-attempt table is the audit/replot source artifact.\n",
        )
        outputs = {
            label: _describe_analysis_output(
                label,
                stage / name,
                table=tables.get(label),
            )
            for label, name in ANALYSIS_OUTPUT_NAMES.items()
        }
        completion = {
            "schema": ANALYSIS_SCHEMA,
            "status": "complete",
            "contract": contract,
            "contract_sha256": _canonical_sha256(contract),
            "outputs": outputs,
        }
        _atomic_json(stage / "analysis_completion.json", completion)
        _validate_analysis_bundle(stage, contract)
        if result_dir.exists():
            raise ValueError(
                "natural-final-AF final bundle appeared during publication"
            )
        os.replace(stage, result_dir)
        published = _validate_analysis_bundle(result_dir, contract)
        return {**published, "cache_hit": False}
    finally:
        if stage.exists():
            _remove_owned_staging(stage, result_parent)


def status_table(plan: NaturalAfPlan) -> pd.DataFrame:
    plan_path = write_plan(plan)
    payload = _load_plan(plan_path)
    work_dir = plan.campaign_dir.resolve() / DEFAULT_WORK_SUBDIR
    rows = []
    for cell in payload["cells"]:
        schedule = _schedule_for_cell(plan.repo_root.resolve(), cell)
        for batch_index in range(int(cell["n_batches"])):
            table, completion = _batch_paths(
                work_dir, str(cell["cell_id"]), batch_index
            )
            start = batch_index * int(cell["batch_size"])
            stop = min(int(cell["attempts"]), start + int(cell["batch_size"]))
            seed = _stable_seed(
                int(payload["rng"]["base_seed"]),
                f"{cell['cell_id']}:batch:{batch_index}",
            )
            contract = _batch_contract(
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
                    _validate_batch(table, completion, contract)
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
        description="Natural final-AF Wright--Fisher study"
    )
    parser.add_argument(
        "phase", choices=("plan", "simulate", "analyze", "status", "all")
    )
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--campaign-dir", type=Path)
    parser.add_argument(
        "--eas-neutral-attempts", type=int, default=DEFAULT_EAS_NEUTRAL_ATTEMPTS
    )
    parser.add_argument(
        "--eas-selected-attempts", type=int, default=DEFAULT_EAS_SELECTED_ATTEMPTS
    )
    parser.add_argument(
        "--han-neutral-attempts", type=int, default=DEFAULT_HAN_NEUTRAL_ATTEMPTS
    )
    parser.add_argument(
        "--han-selected-attempts", type=int, default=DEFAULT_HAN_SELECTED_ATTEMPTS
    )
    parser.add_argument(
        "--eas-neutral-batch-size", type=int, default=DEFAULT_EAS_NEUTRAL_BATCH_SIZE
    )
    parser.add_argument(
        "--eas-selected-batch-size", type=int, default=DEFAULT_EAS_SELECTED_BATCH_SIZE
    )
    parser.add_argument("--han-batch-size", type=int, default=DEFAULT_HAN_BATCH_SIZE)
    parser.add_argument("--base-seed", type=int, default=BASE_SEED)
    parser.add_argument("--workers", type=int, default=4)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    root = args.repo_root.resolve()
    campaign = (args.campaign_dir or root / DEFAULT_CAMPAIGN_DIR).resolve()
    plan = NaturalAfPlan(
        repo_root=root,
        campaign_dir=campaign,
        eas_neutral_attempts=args.eas_neutral_attempts,
        eas_selected_attempts=args.eas_selected_attempts,
        han_neutral_attempts=args.han_neutral_attempts,
        han_selected_attempts=args.han_selected_attempts,
        eas_neutral_batch_size=args.eas_neutral_batch_size,
        eas_selected_batch_size=args.eas_selected_batch_size,
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
    "NaturalAfPlan",
    "TrajectorySchedule",
    "analyze_study",
    "build_cells",
    "build_eas_schedule",
    "build_han_schedule",
    "build_parser",
    "main",
    "selected_gamete_frequency",
    "simulate_study",
    "simulate_trajectory_batch",
    "summarize_attempts",
    "write_plan",
]
