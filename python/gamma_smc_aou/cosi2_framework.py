"""Render and validate CoSi2 EAS/Han sweep-model framework inputs.

This module deliberately stops before production simulation.  It converts the
checked-in EAS PHLASH pointwise-median trajectory and the stdpopsim 0.3.0
``OutOfAfricaArchaicAdmixture_5R19`` catalog model into CoSi2 v2.4.0 parameter
files, emits an explicit selected-allele timing contract, and makes a
publication-sized demographic-model figure.

CoSi2 conditions a stochastic one-locus trajectory on *population* frequency
at the present.  It does not infer mutation age from ``s`` and endpoint AF, nor
does it natively condition an intermediate introgressed-allele frequency.  The
EAS birth ages are deterministic initialization centers; accepted CoSi2
trajectories remain stochastic.  For Han/CHB, selection starts at the integer
generation when 5R19's continuous Neanderthal migration ends.  Saved
trajectories must additionally pass the user-requested operational CHB
allele-frequency screen at that boundary.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import math
import os
import re
import shutil
import sys
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from .eas_sweep_models import (
    GENERATION_TIME_YEARS,
    MUTATION_RATE,
    RECOMBINATION_RATE,
    SEQUENCE_LENGTH_BP,
    load_phlash_eas_npz,
)
from .eas_sweep_study import BASE_SEED, EAS_RESOURCE_SHA256

matplotlib.use("Agg", force=True)


SCHEMA = "gamma-smc.cosi2-framework/v2"
COMPLETION_SCHEMA = "gamma-smc.cosi2-framework-completion/v2"
MODULE_PATH = "python/gamma_smc_aou/cosi2_framework.py"
WRAPPER_PATH = "scripts/run_cosi2_framework.py"
EAS_RESOURCE_PATH = "no_introgression_EAS_sim/resources/EAS.npz"
DEFAULT_OUTPUT_RELATIVE = "focused_selection_EAS_sim/cosi2_framework"
RUN_ROOT_RELATIVE = "../cosi2_framework_runs"

COSI2_VERSION = "2.4.0"
COSI2_TAG = "v2.4.0"
COSI2_COMMIT = "5da717cb8bf7324fb961e772a0c9042385de381a"
COSI2_REPOSITORY = "https://github.com/broadinstitute/cosi2.git"

FOCAL_RELATIVE_POSITION = 0.5
DEFAULT_SELECTION_COEFFICIENTS = (0.01,)
DOMINANCE_COEFFICIENT = 0.5  # hard-coded by CoSi2 sweep_mult_standing
DEFAULT_PRESENT_AF = 0.20
DEFAULT_PRESENT_AF_TOLERANCE = 0.005
DEFAULT_SAMPLE_HAPLOIDS = 200
DEFAULT_HAN_MUTATION_BIRTH_GENERATIONS = 2_400
HAN_MODEL_ID = "OutOfAfricaArchaicAdmixture_5R19"
HAN_DEMOGRAPHY_ID = "out_of_africa_archaic_admixture_5r19_chb"
HAN_GENERATION_TIME_YEARS = 29.0
HAN_ARCHAIC_ANCESTRY_ESTIMATE = 0.012
HAN_ARCHAIC_ANCESTRY_UNCERTAINTY = 0.006
DEFAULT_HAN_SELECTION_START_AF_TARGET = 0.012
DEFAULT_HAN_SELECTION_START_AF_TOLERANCE = 0.001
HAN_NEANDERTHAL_MIGRATION_RATE = 0.825e-5
HAN_YRI_ARCHAIC_AFR_MIGRATION_RATE = 1.98e-5
HAN_YRI_CEU_MIGRATION_RATE = 2.48e-5
HAN_CEU_CHB_MIGRATION_RATE = 11.3e-5
HAN_YRI_OOA_MIGRATION_RATE = 52.2e-5
HAN_MIGRATION_END_CATALOG_GENERATIONS = 18_700 / HAN_GENERATION_TIME_YEARS
HAN_CHB_SPLIT_CATALOG_GENERATIONS = 36_000 / HAN_GENERATION_TIME_YEARS
HAN_OOA_SPLIT_CATALOG_GENERATIONS = 60_700 / HAN_GENERATION_TIME_YEARS
HAN_ARCHAIC_AFR_MIGRATION_START_CATALOG_GENERATIONS = (
    125_000 / HAN_GENERATION_TIME_YEARS
)
HAN_ANCESTRAL_SIZE_CATALOG_GENERATIONS = 300_000 / HAN_GENERATION_TIME_YEARS
HAN_ARCHAIC_AFR_SPLIT_CATALOG_GENERATIONS = 499_000 / HAN_GENERATION_TIME_YEARS
HAN_NEANDERTHAL_SPLIT_CATALOG_GENERATIONS = 559_000 / HAN_GENERATION_TIME_YEARS
HAN_NEANDERTHAL_CONTACT_DURATION_GENERATIONS = (
    HAN_OOA_SPLIT_CATALOG_GENERATIONS - HAN_MIGRATION_END_CATALOG_GENERATIONS
)
HAN_SYMMETRIC_CONTACT_ANCESTRY_EXPECTATION = (
    1.0
    - math.exp(
        -2.0
        * HAN_NEANDERTHAL_MIGRATION_RATE
        * HAN_NEANDERTHAL_CONTACT_DURATION_GENERATIONS
    )
) / 2.0
HAN_MIGRATION_END_GENERATIONS = 645
HAN_CHB_SPLIT_GENERATIONS = 1_241
HAN_OOA_SPLIT_GENERATIONS = 2_093
HAN_ARCHAIC_AFR_MIGRATION_START_GENERATIONS = 4_310
HAN_ANCESTRAL_SIZE_GENERATIONS = 10_345
HAN_ARCHAIC_AFR_SPLIT_GENERATIONS = 17_207
HAN_NEANDERTHAL_SPLIT_GENERATIONS = 19_276
FIGURE_SIZE_INCHES = (11.0, 8.5)
PNG_DPI = 300
SELECTION_PALETTE = ("#e68613", "#2171b5", "#6a51a3", "#238b45", "#cb181d")
SELECTION_MARKERS = ("o", "s", "^", "D", "P", "X")


@dataclass(frozen=True)
class FrameworkPlan:
    """Resolved inputs for the additive CoSi2 framework."""

    selection_coefficients: tuple[float, ...] = DEFAULT_SELECTION_COEFFICIENTS
    present_af: float = DEFAULT_PRESENT_AF
    present_af_tolerance: float = DEFAULT_PRESENT_AF_TOLERANCE
    sample_haploids: int = DEFAULT_SAMPLE_HAPLOIDS
    han_mutation_birth_generations: int = DEFAULT_HAN_MUTATION_BIRTH_GENERATIONS
    han_selection_start_af_target: float = DEFAULT_HAN_SELECTION_START_AF_TARGET
    han_selection_start_af_tolerance: float = DEFAULT_HAN_SELECTION_START_AF_TOLERANCE
    eas_generation_time_years: float = GENERATION_TIME_YEARS
    han_generation_time_years: float = HAN_GENERATION_TIME_YEARS
    sequence_length_bp: int = SEQUENCE_LENGTH_BP
    mutation_rate: float = MUTATION_RATE
    recombination_rate: float = RECOMBINATION_RATE
    focal_relative_position: float = FOCAL_RELATIVE_POSITION
    base_seed: int = BASE_SEED

    def validate(self) -> None:
        coefficients = tuple(float(value) for value in self.selection_coefficients)
        if not coefficients or len(set(coefficients)) != len(coefficients):
            raise ValueError("selection coefficients must be nonempty and unique")
        if any(not math.isfinite(value) or value <= 0 for value in coefficients):
            raise ValueError("selection coefficients must be finite and positive")
        if not math.isfinite(self.present_af) or not 0 < self.present_af < 1:
            raise ValueError("present_af must be in (0, 1)")
        if (
            not math.isfinite(self.present_af_tolerance)
            or self.present_af_tolerance <= 0
            or self.present_af - self.present_af_tolerance <= 0
            or self.present_af + self.present_af_tolerance >= 1
        ):
            raise ValueError("present_af_tolerance must define an interval in (0, 1)")
        if isinstance(self.sample_haploids, bool) or self.sample_haploids <= 0:
            raise ValueError("sample_haploids must be a positive integer")
        if int(self.sample_haploids) != self.sample_haploids:
            raise ValueError("sample_haploids must be an integer")
        if (
            isinstance(self.han_mutation_birth_generations, bool)
            or not math.isfinite(float(self.han_mutation_birth_generations))
            or int(self.han_mutation_birth_generations)
            != self.han_mutation_birth_generations
        ):
            raise ValueError("Han mutation birth must be an integer generation")
        if not (
            HAN_OOA_SPLIT_GENERATIONS
            < self.han_mutation_birth_generations
            < HAN_NEANDERTHAL_SPLIT_GENERATIONS
        ):
            raise ValueError(
                "Han mutation birth must predate Neanderthal-Eurasian contact and "
                "postdate the human-Neanderthal split"
            )
        if not math.isfinite(self.han_selection_start_af_target) or not (
            0 < self.han_selection_start_af_target < 1
        ):
            raise ValueError("Han selection-start allele-AF target must be in (0, 1)")
        if (
            not math.isfinite(self.han_selection_start_af_tolerance)
            or self.han_selection_start_af_tolerance <= 0
            or self.han_selection_start_af_target
            - self.han_selection_start_af_tolerance
            <= 0
            or self.han_selection_start_af_target
            + self.han_selection_start_af_tolerance
            >= 1
        ):
            raise ValueError("Han selection-start allele-AF gate must lie in (0, 1)")
        if (
            not math.isfinite(self.eas_generation_time_years)
            or self.eas_generation_time_years != GENERATION_TIME_YEARS
        ):
            raise ValueError(
                "the checked-in EAS model fixes generation time at 25 years"
            )
        if (
            not math.isfinite(self.han_generation_time_years)
            or self.han_generation_time_years != HAN_GENERATION_TIME_YEARS
        ):
            raise ValueError("5R19 fixes Han generation time at 29 years")
        if self.sequence_length_bp <= 0:
            raise ValueError("sequence length must be positive")
        if self.mutation_rate <= 0 or self.recombination_rate <= 0:
            raise ValueError("mutation and recombination rates must be positive")
        if not 0 < self.focal_relative_position < 1:
            raise ValueError("focal relative position must lie in (0, 1)")
        if isinstance(self.base_seed, bool) or int(self.base_seed) != self.base_seed:
            raise ValueError("base_seed must be an integer")

    def to_record(self) -> dict[str, Any]:
        self.validate()
        record = asdict(self)
        record["selection_coefficients"] = [
            float(value) for value in self.selection_coefficients
        ]
        return record


def _canonical_json(value: Any) -> str:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    )


def _canonical_sha256(value: Any) -> str:
    return hashlib.sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_source_sha256(path: str | Path) -> str:
    text = Path(path).read_text(encoding="utf-8")
    return hashlib.sha256(text.replace("\r\n", "\n").encode("utf-8")).hexdigest()


def _half_up(value: float) -> int:
    if not math.isfinite(value) or value <= 0:
        raise ValueError("population size must be finite and positive")
    return math.floor(value + 0.5)


def _fmt(value: float) -> str:
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    return format(float(value), ".12g")


def _selection_coefficient_index(selection_coefficient: float) -> int:
    """Return a stable coefficient index while preserving legacy encodings."""

    token = _fmt(selection_coefficient)
    legacy_offsets = {"0.005": 0, "0.01": 1}
    if token in legacy_offsets:
        return legacy_offsets[token]
    digest = hashlib.sha256(token.encode("ascii")).digest()
    return 1_000 + int.from_bytes(digest[:4], "big") % 1_000_000


def _selection_plot_style(selection_coefficient: float) -> tuple[str, str]:
    index = _selection_coefficient_index(selection_coefficient)
    return (
        SELECTION_PALETTE[index % len(SELECTION_PALETTE)],
        SELECTION_MARKERS[index % len(SELECTION_MARKERS)],
    )


def _han_bar_offsets(cell_count: int) -> np.ndarray:
    if cell_count <= 0:
        raise ValueError("Han plot requires at least one model cell")
    if cell_count == 1:
        return np.zeros(1, dtype=float)
    return np.linspace(-0.11, 0.11, cell_count)


def _atomic_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        temporary.write_text(value, encoding="utf-8", newline="\n")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_frame(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        frame.to_csv(temporary, sep="\t", index=False, lineterminator="\n")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _package_version(name: str) -> str:
    return importlib.metadata.version(name)


def load_eas_history(repo_root: str | Path) -> pd.DataFrame:
    """Return exact PHLASH float and half-up integer EAS histories."""

    root = Path(repo_root).resolve()
    resource = root / EAS_RESOURCE_PATH
    artifact = load_phlash_eas_npz(resource, expected_sha256=EAS_RESOURCE_SHA256)
    source = artifact.quantile_trajectories()["median"].astype(float)
    times = np.concatenate(([0.0], artifact.time_generations.astype(float)))
    values = np.concatenate(([float(source[0])], source))
    integer = np.floor(values + 0.5).astype(np.int64)
    frame = pd.DataFrame(
        {
            "time_generations": times,
            "time_years": times * GENERATION_TIME_YEARS,
            "raw_median_ne": values,
            "cosi2_integer_ne": integer,
        }
    )
    if len(frame) != 10_001 or not frame["time_generations"].is_monotonic_increasing:
        raise RuntimeError("unexpected EAS PHLASH time grid")
    if int(frame.iloc[0]["cosi2_integer_ne"]) != 37_637:
        raise RuntimeError("EAS present integer Ne changed")
    return frame


def rle_eas_history(history: pd.DataFrame) -> pd.DataFrame:
    """Collapse only consecutive identical integer Ne values."""

    required = {"time_generations", "raw_median_ne", "cosi2_integer_ne"}
    if set(history.columns).issuperset(required) is False:
        raise ValueError("EAS history schema is incomplete")
    keep = history["cosi2_integer_ne"].ne(history["cosi2_integer_ne"].shift())
    rle = history.loc[keep].reset_index(drop=True).copy()
    rle.insert(0, "state_index", np.arange(len(rle), dtype=np.int64))
    if len(rle) != 4_667:
        raise RuntimeError(f"unexpected EAS integer RLE state count: {len(rle)}")
    if int(rle.iloc[-1]["cosi2_integer_ne"]) != 13_138:
        raise RuntimeError("EAS ancient terminal integer Ne changed")
    return rle


def _logit(value: float) -> float:
    if not 0 < value < 1:
        raise ValueError("logit input must lie in (0, 1)")
    return math.log(value / (1.0 - value))


def _expit(value: np.ndarray | float) -> np.ndarray | float:
    array = np.asarray(value, dtype=float)
    result = 1.0 / (1.0 + np.exp(-array))
    return float(result) if result.ndim == 0 else result


def deterministic_final_af(
    initial_af: np.ndarray | float,
    selection_coefficient: float,
    selected_generations: np.ndarray | float,
) -> np.ndarray | float:
    """Additive logistic initialization approximation used only for age centers."""

    initial = np.asarray(initial_af, dtype=float)
    generations = np.asarray(selected_generations, dtype=float)
    if np.any((initial <= 0) | (initial >= 1)) or np.any(generations < 0):
        raise ValueError("initial AF and selected generations are invalid")
    s = float(selection_coefficient)
    if not math.isfinite(s) or s <= 0:
        raise ValueError("selection coefficient must be positive")
    return _expit(
        np.log(initial / (1.0 - initial)) + DOMINANCE_COEFFICIENT * s * generations
    )


def derive_eas_birth_age(
    history: pd.DataFrame,
    *,
    selection_coefficient: float,
    target_af: float,
) -> dict[str, Any]:
    """Solve the one-copy PHLASH fixed-point age on integer generations."""

    max_age = math.floor(float(history["time_generations"].max()))
    ages = np.arange(1, max_age + 1, dtype=np.int64)
    times = history["time_generations"].to_numpy(dtype=float)
    sizes = history["cosi2_integer_ne"].to_numpy(dtype=np.int64)
    indices = np.searchsorted(times, ages.astype(float), side="right") - 1
    ne = sizes[np.clip(indices, 0, len(sizes) - 1)]
    initial = 1.0 / (2.0 * ne.astype(float))
    estimated_age = (_logit(float(target_af)) - np.log(initial / (1.0 - initial))) / (
        DOMINANCE_COEFFICIENT * float(selection_coefficient)
    )
    error = np.abs(ages.astype(float) - estimated_age)
    best = int(np.argmin(error))
    age = int(ages[best])
    record = {
        "age_generations": age,
        "age_years": age * GENERATION_TIME_YEARS,
        "age_kya": age * GENERATION_TIME_YEARS / 1_000.0,
        "diploid_ne_at_birth": int(ne[best]),
        "initial_af_one_copy": float(initial[best]),
        "fixed_point_residual_generations": float(error[best]),
        "deterministic_final_af": float(
            deterministic_final_af(initial[best], selection_coefficient, age)
        ),
    }
    return record


def derive_han_selection_onset(
    *,
    selection_coefficient: float,
    initial_af: float,
    target_af: float,
) -> dict[str, Any]:
    """Return the duration required for a deterministic Han AF bridge.

    This is a diagnostic only.  The 5R19 framework fixes selection onset at the
    end of archaic migration (generation 645) for every coefficient.
    """

    exact = (_logit(target_af) - _logit(initial_af)) / (
        DOMINANCE_COEFFICIENT * float(selection_coefficient)
    )
    age = math.floor(exact + 0.5)
    if age <= 0:
        raise ValueError("derived Han selection duration must be positive")
    return {
        "age_generations": age,
        "age_generations_unrounded": float(exact),
        "age_years": age * HAN_GENERATION_TIME_YEARS,
        "age_kya": age * HAN_GENERATION_TIME_YEARS / 1_000.0,
        "initial_af_assumption": float(initial_af),
        "deterministic_final_af": float(
            deterministic_final_af(initial_af, selection_coefficient, age)
        ),
    }


def _present_af_token(plan: FrameworkPlan) -> str:
    lower = plan.present_af - plan.present_af_tolerance
    upper = plan.present_af + plan.present_af_tolerance
    return f"{_fmt(lower)}-{_fmt(upper)}"


def _recombination_map_text(plan: FrameworkPlan) -> str:
    """Render a constant map accepted by CoSi2 v2.4's strict-positive parser."""

    return f"1\t{_fmt(plan.recombination_rate)}\n"


def build_model_cells(history: pd.DataFrame, plan: FrameworkPlan) -> pd.DataFrame:
    plan.validate()
    rows: list[dict[str, Any]] = []
    coefficients = sorted(plan.selection_coefficients)
    seed_offsets = [_selection_coefficient_index(value) for value in coefficients]
    if len(seed_offsets) != len(set(seed_offsets)):
        raise RuntimeError("selection coefficients produced a seed-offset collision")
    for coefficient, seed_offset in zip(coefficients, seed_offsets, strict=True):
        eas = derive_eas_birth_age(
            history,
            selection_coefficient=coefficient,
            target_af=plan.present_af,
        )
        rows.append(
            {
                "cell_id": f"eas_s{str(coefficient).replace('.', 'p')}",
                "demography_id": "eas_phlash_median",
                "sample_population": "EAS",
                "birth_population": "EAS",
                "selection_population": "EAS",
                "selection_coefficient": coefficient,
                "dominance_coefficient": DOMINANCE_COEFFICIENT,
                "mutation_birth_generations_ago": eas["age_generations"],
                "selection_onset_generations_ago": eas["age_generations"],
                "selection_onset_kya": eas["age_kya"],
                "source_generation_time_years": plan.eas_generation_time_years,
                "initial_af_timing_assumption": eas["initial_af_one_copy"],
                "deterministic_present_af_from_gate_center": plan.present_af,
                "target_bridge_duration_generations": eas["age_generations"],
                "present_af_target": plan.present_af,
                "present_af_lower": plan.present_af - plan.present_af_tolerance,
                "present_af_upper": plan.present_af + plan.present_af_tolerance,
                "archaic_migration_end_generations_ago": np.nan,
                "archaic_migration_start_generations_ago": np.nan,
                "archaic_migration_rate_per_generation": np.nan,
                "introgression_ancestry_fraction": np.nan,
                "introgressed_allele_af_target": np.nan,
                "introgressed_allele_af_lower": np.nan,
                "introgressed_allele_af_upper": np.nan,
                "introgression_enforcement": "not_applicable",
                "seed": int(plan.base_seed + seed_offset),
                "native_cosi2_status": "native_endpoint_conditioned",
            }
        )
        bridge = derive_han_selection_onset(
            selection_coefficient=coefficient,
            initial_af=plan.han_selection_start_af_target,
            target_af=plan.present_af,
        )
        fixed_endpoint = deterministic_final_af(
            plan.han_selection_start_af_target,
            coefficient,
            HAN_MIGRATION_END_GENERATIONS,
        )
        rows.append(
            {
                "cell_id": f"han_introgressed_s{str(coefficient).replace('.', 'p')}",
                "demography_id": HAN_DEMOGRAPHY_ID,
                "sample_population": "Han",
                "birth_population": "Neanderthal",
                "selection_population": "CHB",
                "selection_coefficient": coefficient,
                "dominance_coefficient": DOMINANCE_COEFFICIENT,
                "mutation_birth_generations_ago": (plan.han_mutation_birth_generations),
                "selection_onset_generations_ago": HAN_MIGRATION_END_GENERATIONS,
                "selection_onset_kya": (
                    HAN_MIGRATION_END_GENERATIONS
                    * plan.han_generation_time_years
                    / 1_000.0
                ),
                "source_generation_time_years": plan.han_generation_time_years,
                "initial_af_timing_assumption": plan.han_selection_start_af_target,
                "deterministic_present_af_from_gate_center": float(fixed_endpoint),
                "target_bridge_duration_generations": bridge["age_generations"],
                "present_af_target": plan.present_af,
                "present_af_lower": plan.present_af - plan.present_af_tolerance,
                "present_af_upper": plan.present_af + plan.present_af_tolerance,
                "archaic_migration_end_generations_ago": (
                    HAN_MIGRATION_END_GENERATIONS
                ),
                "archaic_migration_start_generations_ago": HAN_OOA_SPLIT_GENERATIONS,
                "archaic_migration_rate_per_generation": (
                    HAN_NEANDERTHAL_MIGRATION_RATE
                ),
                "introgression_ancestry_fraction": HAN_ARCHAIC_ANCESTRY_ESTIMATE,
                "introgressed_allele_af_target": (plan.han_selection_start_af_target),
                "introgressed_allele_af_lower": (
                    plan.han_selection_start_af_target
                    - plan.han_selection_start_af_tolerance
                ),
                "introgressed_allele_af_upper": (
                    plan.han_selection_start_af_target
                    + plan.han_selection_start_af_tolerance
                ),
                "introgression_enforcement": (
                    "posthoc_COSI_SAVE_TRAJ_CHB_AF_gate_at_migration_end"
                ),
                "seed": int(plan.base_seed + 100 + seed_offset),
                "native_cosi2_status": (
                    "native_donor_origin_continuous_migration_and_endpoint; "
                    "migration_end_AF_requires_gate"
                ),
            }
        )
    frame = pd.DataFrame(rows).sort_values(
        ["demography_id", "selection_coefficient"], kind="mergesort"
    )
    if len(frame) != 2 * len(plan.selection_coefficients):
        raise RuntimeError("model cell cardinality is inconsistent")
    return frame.reset_index(drop=True)


def build_age_response(history: pd.DataFrame, plan: FrameworkPlan) -> pd.DataFrame:
    """Build the deterministic age/onset response plotted by the framework."""

    rows: list[dict[str, Any]] = []
    times = history["time_generations"].to_numpy(dtype=float)
    sizes = history["cosi2_integer_ne"].to_numpy(dtype=np.int64)
    for coefficient in sorted(plan.selection_coefficients):
        target_age = int(
            derive_eas_birth_age(
                history,
                selection_coefficient=coefficient,
                target_af=plan.present_af,
            )["age_generations"]
        )
        maximum_age = min(
            math.floor(float(times[-1])),
            max(3_500, math.ceil(1.08 * target_age)),
        )
        eas_ages = np.unique(
            np.append(
                np.arange(1, maximum_age + 1, 10, dtype=np.int64),
                [target_age, maximum_age],
            )
        )
        indices = np.searchsorted(times, eas_ages.astype(float), side="right") - 1
        ne = sizes[np.clip(indices, 0, len(sizes) - 1)]
        initial = 1.0 / (2.0 * ne.astype(float))
        endpoints = deterministic_final_af(initial, coefficient, eas_ages)
        rows.extend(
            {
                "demography_id": "eas_phlash_median",
                "selection_coefficient": coefficient,
                "time_variable": "mutation_birth_and_selection_onset_generations_ago",
                "time_generations": int(age),
                "time_kya": float(age * GENERATION_TIME_YEARS / 1_000.0),
                "initial_af_assumption": float(p0),
                "deterministic_final_af": float(endpoint),
            }
            for age, p0, endpoint in zip(eas_ages, initial, endpoints)
        )
        han_ages = np.arange(0, HAN_MIGRATION_END_GENERATIONS + 1, 2, dtype=np.int64)
        if han_ages[-1] != HAN_MIGRATION_END_GENERATIONS:
            han_ages = np.append(han_ages, HAN_MIGRATION_END_GENERATIONS)
        elapsed = HAN_MIGRATION_END_GENERATIONS - han_ages
        endpoints = deterministic_final_af(
            plan.han_selection_start_af_target, coefficient, elapsed
        )
        rows.extend(
            {
                "demography_id": HAN_DEMOGRAPHY_ID,
                "selection_coefficient": coefficient,
                "time_variable": "fixed_selection_path_generations_before_present",
                "time_generations": int(age),
                "time_kya": float(age * plan.han_generation_time_years / 1_000.0),
                "initial_af_assumption": plan.han_selection_start_af_target,
                "deterministic_final_af": float(endpoint),
            }
            for age, endpoint in zip(han_ages, endpoints)
        )
    return pd.DataFrame(rows)


_HAN_POPULATIONS = (
    (1, "YRI", 13_900.0, 13_900),
    (2, "CEU_OOA", 10_855.080951853866, 10_855),
    (3, "Han_CHB", 65_834.77001122756, 65_835),
    (4, "Neanderthal", 3_600.0, 3_600),
    (5, "ArchaicAFR", 3_600.0, 3_600),
)


def build_han_event_ledger(plan: FrameworkPlan) -> pd.DataFrame:
    """Return the integer-realized stdpopsim 5R19 event ledger."""

    plan.validate()

    def event(
        event_type: str,
        label: str,
        catalog_time: float,
        cosi2_time: float,
        population: int,
        other_population: int | None = None,
        *,
        migration_rate: float | None = None,
        integer_size: int | None = None,
        catalog_start_time: float | None = None,
        cosi2_start_time: float | None = None,
        integer_start_size: int | None = None,
    ) -> dict[str, Any]:
        nominal = math.floor(cosi2_time + 0.5)
        return {
            "event_type": event_type,
            "label": label,
            "catalog_time_years": catalog_time * HAN_GENERATION_TIME_YEARS,
            "catalog_time_generations": catalog_time,
            "cosi2_nominal_integer_generation": nominal,
            "cosi2_time_generations": cosi2_time,
            "rounding_delta_generations": nominal - catalog_time,
            "rounding_delta_years": (nominal - catalog_time)
            * HAN_GENERATION_TIME_YEARS,
            "population_id": population,
            "other_population_id": other_population,
            "migration_rate": migration_rate,
            "proportion": None,
            "integer_size": integer_size,
            "catalog_start_time_generations": catalog_start_time,
            "cosi2_start_time_generations": cosi2_start_time,
            "integer_start_size": integer_start_size,
            "rounding_policy": (
                "nearest_integer_half_up; same-generation net-equivalent order retained"
            ),
        }

    events = [
        event(
            "exp_change_size",
            "CEU_growth",
            0.0,
            0.0,
            2,
            integer_size=10_855,
            catalog_start_time=HAN_CHB_SPLIT_CATALOG_GENERATIONS,
            cosi2_start_time=HAN_CHB_SPLIT_GENERATIONS,
            integer_start_size=2_300,
        ),
        event(
            "exp_change_size",
            "CHB_growth",
            0.0,
            0.0,
            3,
            integer_size=65_835,
            catalog_start_time=HAN_CHB_SPLIT_CATALOG_GENERATIONS,
            cosi2_start_time=HAN_CHB_SPLIT_GENERATIONS,
            integer_start_size=650,
        ),
        event(
            "migration_rate",
            "YRI_CEU_recent_a",
            0.0,
            0.0,
            1,
            2,
            migration_rate=HAN_YRI_CEU_MIGRATION_RATE,
        ),
        event(
            "migration_rate",
            "YRI_CEU_recent_b",
            0.0,
            0.0,
            2,
            1,
            migration_rate=HAN_YRI_CEU_MIGRATION_RATE,
        ),
        event(
            "migration_rate",
            "CEU_CHB_recent_a",
            0.0,
            0.0,
            2,
            3,
            migration_rate=HAN_CEU_CHB_MIGRATION_RATE,
        ),
        event(
            "migration_rate",
            "CEU_CHB_recent_b",
            0.0,
            0.0,
            3,
            2,
            migration_rate=HAN_CEU_CHB_MIGRATION_RATE,
        ),
        event(
            "migration_rate",
            "YRI_ArchaicAFR_contact_a",
            HAN_MIGRATION_END_CATALOG_GENERATIONS,
            HAN_MIGRATION_END_GENERATIONS,
            1,
            5,
            migration_rate=HAN_YRI_ARCHAIC_AFR_MIGRATION_RATE,
        ),
        event(
            "migration_rate",
            "YRI_ArchaicAFR_contact_b",
            HAN_MIGRATION_END_CATALOG_GENERATIONS,
            HAN_MIGRATION_END_GENERATIONS,
            5,
            1,
            migration_rate=HAN_YRI_ARCHAIC_AFR_MIGRATION_RATE,
        ),
        event(
            "migration_rate",
            "CEU_Neanderthal_contact_a",
            HAN_MIGRATION_END_CATALOG_GENERATIONS,
            HAN_MIGRATION_END_GENERATIONS,
            2,
            4,
            migration_rate=HAN_NEANDERTHAL_MIGRATION_RATE,
        ),
        event(
            "migration_rate",
            "CEU_Neanderthal_contact_b",
            HAN_MIGRATION_END_CATALOG_GENERATIONS,
            HAN_MIGRATION_END_GENERATIONS,
            4,
            2,
            migration_rate=HAN_NEANDERTHAL_MIGRATION_RATE,
        ),
        event(
            "migration_rate",
            "CHB_Neanderthal_contact_a",
            HAN_MIGRATION_END_CATALOG_GENERATIONS,
            HAN_MIGRATION_END_GENERATIONS,
            3,
            4,
            migration_rate=HAN_NEANDERTHAL_MIGRATION_RATE,
        ),
        event(
            "migration_rate",
            "CHB_Neanderthal_contact_b",
            HAN_MIGRATION_END_CATALOG_GENERATIONS,
            HAN_MIGRATION_END_GENERATIONS,
            4,
            3,
            migration_rate=HAN_NEANDERTHAL_MIGRATION_RATE,
        ),
        event(
            "split",
            "CHB_from_CEU_OOA",
            HAN_CHB_SPLIT_CATALOG_GENERATIONS,
            HAN_CHB_SPLIT_GENERATIONS,
            2,
            3,
        ),
        event(
            "change_size",
            "OOA_bottleneck",
            HAN_CHB_SPLIT_CATALOG_GENERATIONS,
            HAN_CHB_SPLIT_GENERATIONS,
            2,
            integer_size=880,
        ),
        event(
            "migration_rate",
            "YRI_OOA_a",
            HAN_CHB_SPLIT_CATALOG_GENERATIONS,
            HAN_CHB_SPLIT_GENERATIONS,
            1,
            2,
            migration_rate=HAN_YRI_OOA_MIGRATION_RATE,
        ),
        event(
            "migration_rate",
            "YRI_OOA_b",
            HAN_CHB_SPLIT_CATALOG_GENERATIONS,
            HAN_CHB_SPLIT_GENERATIONS,
            2,
            1,
            migration_rate=HAN_YRI_OOA_MIGRATION_RATE,
        ),
        event(
            "split",
            "OOA_from_YRI",
            HAN_OOA_SPLIT_CATALOG_GENERATIONS,
            HAN_OOA_SPLIT_GENERATIONS,
            1,
            2,
        ),
        event(
            "migration_rate",
            "YRI_ArchaicAFR_off_a",
            HAN_ARCHAIC_AFR_MIGRATION_START_CATALOG_GENERATIONS,
            HAN_ARCHAIC_AFR_MIGRATION_START_GENERATIONS,
            1,
            5,
            migration_rate=0.0,
        ),
        event(
            "migration_rate",
            "YRI_ArchaicAFR_off_b",
            HAN_ARCHAIC_AFR_MIGRATION_START_CATALOG_GENERATIONS,
            HAN_ARCHAIC_AFR_MIGRATION_START_GENERATIONS,
            5,
            1,
            migration_rate=0.0,
        ),
        event(
            "change_size",
            "ancestral_modern_N",
            HAN_ANCESTRAL_SIZE_CATALOG_GENERATIONS,
            HAN_ANCESTRAL_SIZE_GENERATIONS,
            1,
            integer_size=3_600,
        ),
        event(
            "split",
            "ArchaicAFR_from_modern",
            HAN_ARCHAIC_AFR_SPLIT_CATALOG_GENERATIONS,
            HAN_ARCHAIC_AFR_SPLIT_GENERATIONS,
            1,
            5,
        ),
        event(
            "split",
            "Neanderthal_from_modern",
            HAN_NEANDERTHAL_SPLIT_CATALOG_GENERATIONS,
            HAN_NEANDERTHAL_SPLIT_GENERATIONS,
            1,
            4,
        ),
    ]
    frame = pd.DataFrame(events).sort_values(
        ["cosi2_time_generations", "label"], kind="mergesort"
    )
    frame.insert(0, "demography_id", HAN_DEMOGRAPHY_ID)
    if not frame["cosi2_time_generations"].is_monotonic_increasing:
        raise RuntimeError("Han CoSi2 event times are not increasing")
    return frame.reset_index(drop=True)


def _event_directive(row: Mapping[str, Any]) -> str | None:
    event = str(row["event_type"])
    label = str(row["label"])
    time = _fmt(float(row["cosi2_time_generations"]))
    population = int(row["population_id"])
    other = row.get("other_population_id")
    if event == "migration_rate":
        return (
            f'pop_event migration_rate "{label}" {population} {int(other)} '
            f"{time} {_fmt(float(row['migration_rate']))}"
        )
    if event == "split":
        return f'pop_event split "{label}" {population} {int(other)} {time}'
    if event == "change_size":
        return (
            f'pop_event change_size "{label}" {population} {time} '
            f"{int(row['integer_size'])}"
        )
    if event == "exp_change_size":
        return (
            f'pop_event exp_change_size "{label}" {population} {time} '
            f"{_fmt(float(row['cosi2_start_time_generations']))} "
            f"{int(row['integer_size'])} {int(row['integer_start_size'])}"
        )
    raise ValueError(f"unsupported event type {event}")


def render_eas_parameter_file(
    rle: pd.DataFrame, cell: Mapping[str, Any], plan: FrameworkPlan
) -> str:
    lines = [
        "# Generated CoSi2 v2.4.0 EAS de novo model; do not hand-edit.",
        f"# CoSi2 commit: {COSI2_COMMIT}",
        "length " + str(plan.sequence_length_bp),
        "mutation_rate " + _fmt(plan.mutation_rate),
        "recomb_file recombination_map.tsv",
        "gene_conversion_relative_rate 0",
        "infinite_sites yes",
        "pop_define 1 EAS",
        f"pop_size 1 {int(rle.iloc[0]['cosi2_integer_ne'])}",
        f"sample_size 1 {plan.sample_haploids}",
    ]
    for row in rle.iloc[1:].itertuples(index=False):
        lines.append(
            f'pop_event change_size "eas_ne_{int(row.state_index):06d}" 1 '
            f"{_fmt(row.time_generations)} {int(row.cosi2_integer_ne)}"
        )
    lines.extend(
        [
            (
                f'pop_event sweep_mult_standing "{cell["cell_id"]}" 1 '
                f"{int(cell['mutation_birth_generations_ago'])} "
                f"{_fmt(cell['selection_coefficient'])} "
                f"{_fmt(plan.focal_relative_position)} {_present_af_token(plan)} 1 "
                f"{int(cell['selection_onset_generations_ago'])}"
            ),
            f"random_seed {int(cell['seed'])}",
        ]
    )
    return "\n".join(lines) + "\n"


def render_han_parameter_file(
    event_ledger: pd.DataFrame, cell: Mapping[str, Any], plan: FrameworkPlan
) -> str:
    gate_lower = (
        plan.han_selection_start_af_target - plan.han_selection_start_af_tolerance
    )
    gate_upper = (
        plan.han_selection_start_af_target + plan.han_selection_start_af_tolerance
    )
    lines = [
        "# Generated CoSi2 v2.4.0 5R19 Han/CHB model; do not hand-edit.",
        f"# CoSi2 commit: {COSI2_COMMIT}",
        "# stdpopsim model: OutOfAfricaArchaicAdmixture_5R19 (v0.3.0).",
        "# 5R19 uses continuous symmetric migration, not an admixture pulse.",
        (
            "# The operational CHB allele-AF gate at generation "
            f"{HAN_MIGRATION_END_GENERATIONS} is "
            f"{_fmt(gate_lower)}-{_fmt(gate_upper)} and is applied post hoc."
        ),
        "# Mutation/recombination values below are study overrides; 5R19 mutation_rate=None.",
        f"length {plan.sequence_length_bp}",
        f"mutation_rate {_fmt(plan.mutation_rate)}",
        "recomb_file recombination_map.tsv",
        "gene_conversion_relative_rate 0",
        "infinite_sites yes",
    ]
    for identifier, name, _raw_size, integer_size in _HAN_POPULATIONS:
        lines.append(f"pop_define {identifier} {name}")
        lines.append(f"pop_size {identifier} {integer_size}")
        sample_size = plan.sample_haploids if name == "Han_CHB" else 0
        lines.append(f"sample_size {identifier} {sample_size}")
    for _, row in event_ledger.iterrows():
        directive = _event_directive(row)
        if directive is not None:
            lines.append(directive)
    lines.extend(
        [
            (
                f'pop_event sweep_mult_standing "{cell["cell_id"]}" 4 '
                f"{int(cell['mutation_birth_generations_ago'])} "
                f"{_fmt(cell['selection_coefficient'])} "
                f"{_fmt(plan.focal_relative_position)} {_present_af_token(plan)} 3 "
                f"{int(cell['selection_onset_generations_ago'])}"
            ),
            f"random_seed {int(cell['seed'])}",
        ]
    )
    return "\n".join(lines) + "\n"


def lint_parameter_text(
    text: str, *, expected_cell: Mapping[str, Any]
) -> dict[str, Any]:
    """Strictly lint the subset of CoSi2 syntax emitted by this module."""

    lines = [line.strip() for line in text.splitlines() if line.strip()]
    directives = [line for line in lines if not line.startswith("#")]
    sweep = [
        line for line in directives if line.startswith("pop_event sweep_mult_standing ")
    ]
    if len(sweep) != 1:
        raise ValueError(
            "parameter file must contain exactly one modern sweep directive"
        )
    tokens = next(csv.reader([sweep[0]], delimiter=" ", quotechar='"'))
    tokens = [token for token in tokens if token]
    if len(tokens) != 10 or tokens[:3] != [
        "pop_event",
        "sweep_mult_standing",
        expected_cell["cell_id"],
    ]:
        raise ValueError("modern sweep directive has the wrong field inventory")
    birth_age = float(tokens[4])
    birth_population = int(tokens[3])
    coefficient = float(tokens[5])
    endpoint = tokens[7]
    selection_population = int(tokens[8])
    onset_age = float(tokens[9])
    expected_population_ids = {
        "EAS": (1, 1),
        "Han": (4, 3),
    }
    population_key = str(expected_cell["sample_population"])
    if population_key not in expected_population_ids:
        raise ValueError(f"unsupported sample population {population_key}")
    expected_birth, expected_selection = expected_population_ids[population_key]
    if (birth_population, selection_population) != (
        expected_birth,
        expected_selection,
    ):
        raise ValueError("sweep birth/selection population IDs are incorrect")
    if not math.isclose(
        birth_age, float(expected_cell["mutation_birth_generations_ago"])
    ):
        raise ValueError("sweep birth age differs from the model cell")
    if not math.isclose(coefficient, float(expected_cell["selection_coefficient"])):
        raise ValueError("sweep coefficient differs from the model cell")
    if not math.isclose(
        onset_age, float(expected_cell["selection_onset_generations_ago"])
    ):
        raise ValueError("selection onset differs from the model cell")
    if "-" not in endpoint:
        raise ValueError("framework requires an explicit reachable present-AF interval")
    endpoint_parts = endpoint.split("-")
    if len(endpoint_parts) != 2:
        raise ValueError("present-AF interval has the wrong syntax")
    try:
        endpoint_lower, endpoint_upper = map(float, endpoint_parts)
    except ValueError as exc:
        raise ValueError("present-AF interval is not numeric") from exc
    if not (
        math.isclose(endpoint_lower, float(expected_cell["present_af_lower"]))
        and math.isclose(endpoint_upper, float(expected_cell["present_af_upper"]))
    ):
        raise ValueError("present-AF interval differs from the model cell")
    labels = [line for line in directives if line.startswith("pop_define ")]
    samples = [line for line in directives if line.startswith("sample_size ")]
    label_ids: set[int] = set()
    for line in labels:
        fields = line.split()
        if len(fields) != 3:
            raise ValueError("population definition has the wrong syntax")
        try:
            population_id = int(fields[1])
        except ValueError as exc:
            raise ValueError("population definition ID is not an integer") from exc
        if population_id in label_ids:
            raise ValueError("population definition ID is duplicated")
        label_ids.add(population_id)
    sample_sizes: dict[int, int] = {}
    for line in samples:
        fields = line.split()
        if len(fields) != 3:
            raise ValueError("sample-size request has the wrong syntax")
        try:
            population_id, sample_size = map(int, fields[1:])
        except ValueError as exc:
            raise ValueError("sample-size request is not integer-valued") from exc
        if population_id in sample_sizes or sample_size < 0:
            raise ValueError("sample-size request is duplicated or negative")
        sample_sizes[population_id] = sample_size
    if population_key == "EAS":
        if label_ids != {1} or set(sample_sizes) != {1} or sample_sizes[1] <= 0:
            raise ValueError("EAS population or sample definition is invalid")
    else:
        expected_han_ids = {1, 2, 3, 4, 5}
        if label_ids != expected_han_ids or set(sample_sizes) != expected_han_ids:
            raise ValueError(
                "5R19 requires one sample_size request for every defined population"
            )
        if sample_sizes[3] <= 0 or any(
            sample_sizes[identifier] != 0 for identifier in expected_han_ids - {3}
        ):
            raise ValueError(
                "5R19 must sample CHB only and explicitly request zero samples elsewhere"
            )
    if "recomb_file recombination_map.tsv" not in directives:
        raise ValueError("parameter file must use the framework recombination map")
    migration_directives = [
        line for line in directives if line.startswith("pop_event migration_rate ")
    ]
    if population_key == "Han":
        if any(line.startswith("pop_event admix ") for line in directives):
            raise ValueError("5R19 must not be rendered as a pulse admixture model")
        if len(migration_directives) != 14:
            raise ValueError("5R19 migration schedule is incomplete")
        expected_migrations = {
            f'pop_event migration_rate "YRI_CEU_recent_a" 1 2 0 {_fmt(HAN_YRI_CEU_MIGRATION_RATE)}',
            f'pop_event migration_rate "YRI_CEU_recent_b" 2 1 0 {_fmt(HAN_YRI_CEU_MIGRATION_RATE)}',
            f'pop_event migration_rate "CEU_CHB_recent_a" 2 3 0 {_fmt(HAN_CEU_CHB_MIGRATION_RATE)}',
            f'pop_event migration_rate "CEU_CHB_recent_b" 3 2 0 {_fmt(HAN_CEU_CHB_MIGRATION_RATE)}',
            f'pop_event migration_rate "YRI_ArchaicAFR_contact_a" 1 5 645 {_fmt(HAN_YRI_ARCHAIC_AFR_MIGRATION_RATE)}',
            f'pop_event migration_rate "YRI_ArchaicAFR_contact_b" 5 1 645 {_fmt(HAN_YRI_ARCHAIC_AFR_MIGRATION_RATE)}',
            f'pop_event migration_rate "CEU_Neanderthal_contact_a" 2 4 645 {_fmt(HAN_NEANDERTHAL_MIGRATION_RATE)}',
            f'pop_event migration_rate "CEU_Neanderthal_contact_b" 4 2 645 {_fmt(HAN_NEANDERTHAL_MIGRATION_RATE)}',
            f'pop_event migration_rate "CHB_Neanderthal_contact_a" 3 4 645 {_fmt(HAN_NEANDERTHAL_MIGRATION_RATE)}',
            f'pop_event migration_rate "CHB_Neanderthal_contact_b" 4 3 645 {_fmt(HAN_NEANDERTHAL_MIGRATION_RATE)}',
            f'pop_event migration_rate "YRI_OOA_a" 1 2 1241 {_fmt(HAN_YRI_OOA_MIGRATION_RATE)}',
            f'pop_event migration_rate "YRI_OOA_b" 2 1 1241 {_fmt(HAN_YRI_OOA_MIGRATION_RATE)}',
            'pop_event migration_rate "YRI_ArchaicAFR_off_a" 1 5 4310 0',
            'pop_event migration_rate "YRI_ArchaicAFR_off_b" 5 1 4310 0',
        }
        if set(migration_directives) != expected_migrations:
            raise ValueError("5R19 migration schedule differs from the catalog mapping")
    return {
        "n_directives": len(directives),
        "n_populations": len(labels),
        "n_sample_requests": len(sample_sizes),
        "n_positive_sample_requests": sum(value > 0 for value in sample_sizes.values()),
        "n_change_size": sum(" change_size " in line for line in directives),
        "n_migration_rate": len(migration_directives),
        "endpoint_token": endpoint,
    }


def validate_saved_trajectory(
    trajectory_path: str | Path,
    *,
    cell: Mapping[str, Any],
) -> dict[str, Any]:
    """Validate a modern trajectory and the Han migration-end AF screen."""

    path = Path(trajectory_path)
    if not path.is_file() or path.is_symlink():
        raise ValueError("trajectory must be a regular file")
    frame = pd.read_csv(path, sep="\t")
    required = {"sim", "gen"}
    if not required.issubset(frame.columns):
        raise ValueError("trajectory lacks sim/gen columns")
    frequency_columns = [column for column in frame if column.startswith("selfreq_")]
    if not frequency_columns:
        raise ValueError("trajectory lacks population-frequency columns")
    numeric = frame[frequency_columns].apply(pd.to_numeric, errors="raise")
    if not np.isfinite(numeric.to_numpy(dtype=float)).all():
        raise ValueError("trajectory contains nonfinite frequencies")
    if ((numeric < 0) | (numeric > 1)).any(axis=None):
        raise ValueError("trajectory frequency lies outside [0,1]")
    if frame["sim"].isna().any():
        raise ValueError("trajectory contains a missing simulation identifier")
    numeric_generations = pd.to_numeric(frame["gen"], errors="raise")
    if not np.isfinite(numeric_generations.to_numpy(dtype=float)).all():
        raise ValueError("trajectory contains a nonfinite generation")
    is_present = numeric_generations == 0
    is_boundary = np.isclose(
        numeric_generations.to_numpy(dtype=float),
        float(cell.get("archaic_migration_end_generations_ago", np.nan)),
        rtol=0.0,
        atol=1e-9,
    )
    for _sim, group in frame.groupby("sim", sort=False):
        generations = pd.to_numeric(group["gen"], errors="raise").to_numpy(dtype=float)
        if len(generations) < 2 or not np.all(np.diff(generations) < 0):
            raise ValueError("trajectory generations must be strictly descending")
        group_indices = group.index.to_numpy()
        if int(is_present.loc[group_indices].sum()) != 1:
            raise ValueError("each simulation must contain exactly one present-day row")
        if (
            str(cell["sample_population"]) == "Han"
            and int(np.count_nonzero(is_boundary[group_indices])) != 1
        ):
            raise ValueError(
                "each Han simulation must contain exactly one migration-end row"
            )
    sample_population = str(cell["sample_population"])
    sample_population_ids = {"EAS": "1", "Han": "3"}
    if sample_population not in sample_population_ids:
        raise ValueError(f"unsupported sample population: {sample_population}")
    endpoint_candidates = [
        column
        for column in frequency_columns
        if column.casefold()
        in {
            f"selfreq_{sample_population}".casefold(),
            f"selfreq_{sample_population_ids[sample_population]}",
        }
    ]
    if len(endpoint_candidates) != 1:
        raise ValueError("cannot identify the sample-population trajectory column")
    present = frame.loc[is_present]
    endpoint = pd.to_numeric(present[endpoint_candidates[0]], errors="raise")
    lower = float(cell["present_af_lower"])
    upper = float(cell["present_af_upper"])
    if not endpoint.between(lower, upper, inclusive="both").all():
        raise ValueError("trajectory violates the present population-AF interval")
    record: dict[str, Any] = {
        "status": "valid",
        "n_rows": len(frame),
        "n_simulations": int(frame["sim"].nunique()),
        "present_af_min": float(endpoint.min()),
        "present_af_max": float(endpoint.max()),
        "selection_start_gate_required": str(cell["sample_population"]) == "Han",
    }
    if str(cell["sample_population"]) == "Han":
        chb_candidates = [
            column
            for column in frequency_columns
            if column.casefold() in {"selfreq_chb", "selfreq_han_chb", "selfreq_3"}
        ]
        if len(chb_candidates) != 1:
            raise ValueError("cannot identify the CHB trajectory column")
        boundary_rows = frame.loc[is_boundary]
        boundary_af = pd.to_numeric(boundary_rows[chb_candidates[0]], errors="raise")
        boundary_lower = float(cell["introgressed_allele_af_lower"])
        boundary_upper = float(cell["introgressed_allele_af_upper"])
        if not boundary_af.between(
            boundary_lower, boundary_upper, inclusive="both"
        ).all():
            raise ValueError("trajectory violates the migration-end allele-AF gate")
        record.update(
            {
                "selection_start_af_min": float(boundary_af.min()),
                "selection_start_af_max": float(boundary_af.max()),
                "selection_start_gate_passed": True,
            }
        )
    return record


def _plot_framework(
    output_png: Path,
    output_pdf: Path,
    *,
    history: pd.DataFrame,
    age_response: pd.DataFrame,
    cells: pd.DataFrame,
    plan: FrameworkPlan,
) -> None:
    plt.rcParams.update(
        {
            "font.size": 12,
            "axes.titlesize": 16,
            "axes.labelsize": 14,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "legend.fontsize": 12,
            "figure.titlesize": 20,
            "pdf.fonttype": 42,
        }
    )
    figure, axes = plt.subplots(2, 2, figsize=FIGURE_SIZE_INCHES)
    coefficients = sorted(float(value) for value in plan.selection_coefficients)
    colors = {
        coefficient: _selection_plot_style(coefficient)[0]
        for coefficient in coefficients
    }
    markers = {
        coefficient: _selection_plot_style(coefficient)[1]
        for coefficient in coefficients
    }
    figure.suptitle(
        "CoSi2 sweep specifications: EAS PHLASH + 5R19 Han/CHB",
        y=0.975,
        fontsize=18.5,
    )

    ax = axes[0, 0]
    eas_cells = cells[cells["sample_population"] == "EAS"]
    eas_x_max = max(
        3.5,
        1.08 * float(eas_cells["mutation_birth_generations_ago"].max()) / 1_000.0,
    )
    recent = history[history["time_generations"] <= eas_x_max * 1_000.0]
    ax.step(
        recent["time_generations"] / 1_000.0,
        recent["cosi2_integer_ne"],
        where="post",
        color="#343a40",
        linewidth=1.8,
        label="PHLASH median (integer CoSi2 schedule)",
    )
    for index, row in enumerate(eas_cells.itertuples(index=False)):
        coefficient = float(row.selection_coefficient)
        x = float(row.mutation_birth_generations_ago) / 1_000.0
        ax.axvline(x, color=colors[coefficient], linewidth=2.0, alpha=0.9)
        offset = (-6, 8) if index % 2 == 0 else (6, 8)
        ax.annotate(
            f"s={coefficient:g}\n{int(row.mutation_birth_generations_ago):,} g",
            xy=(x, 3_350),
            xytext=offset,
            textcoords="offset points",
            color=colors[coefficient],
            ha="right" if index % 2 == 0 else "left",
            va="bottom",
            fontsize=11,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.88, "pad": 1.5},
        )
    ax.set_yscale("log")
    ax.set_xlim(0, eas_x_max)
    ax.set_ylim(2_800, 50_000)
    ax.set_title("A. EAS de novo - PHLASH median $N_e(t)$")
    ax.set_xlabel("Thousands of generations before present (present = 0)")
    ax.set_ylabel("Diploid effective population size, $N_e$")
    ax.grid(axis="y", alpha=0.25)

    ax = axes[0, 1]
    eas_response = age_response[age_response["demography_id"] == "eas_phlash_median"]
    for coefficient, group in eas_response.groupby("selection_coefficient", sort=True):
        coefficient = float(coefficient)
        ax.plot(
            group["time_generations"] / 1_000.0,
            group["deterministic_final_af"],
            color=colors[coefficient],
            linewidth=2.2,
        )
        row = eas_cells[
            np.isclose(eas_cells["selection_coefficient"], coefficient)
        ].iloc[0]
        ax.scatter(
            row["mutation_birth_generations_ago"] / 1_000.0,
            plan.present_af,
            s=65,
            marker=markers[coefficient],
            color=colors[coefficient],
            edgecolor="white",
            linewidth=0.8,
            zorder=4,
        )
    ax.axhspan(
        plan.present_af - plan.present_af_tolerance,
        plan.present_af + plan.present_af_tolerance,
        color="#666666",
        alpha=0.12,
    )
    ax.axhline(plan.present_af, color="#555555", linestyle="--", linewidth=1.3)
    ax.text(
        0.03 * eas_x_max,
        plan.present_af + 0.028,
        f"present population AF = {100 * plan.present_af:g}%",
        fontsize=11,
    )
    ax.set_xlim(0, eas_x_max)
    ax.set_ylim(0, 1.02)
    ax.set_title("B. EAS age initialization curve")
    ax.set_xlabel("Allele age (10$^3$ generations BP)")
    ax.set_ylabel("Deterministic final AF")
    ax.grid(alpha=0.25)

    ax = axes[1, 0]
    ax.set_xlim(0, 2.55)
    ax.set_ylim(-0.72, 3.62)
    lanes = {"Neanderthal": 3.0, "Eurasian B": 2.0, "CEU": 1.0, "Han / CHB": 0.0}
    migration_end_x = HAN_MIGRATION_END_GENERATIONS / 1_000.0
    chb_split_x = HAN_CHB_SPLIT_GENERATIONS / 1_000.0
    ooa_split_x = HAN_OOA_SPLIT_GENERATIONS / 1_000.0
    ax.hlines(lanes["Neanderthal"], 0, 2.55, color="#666666", linewidth=2.2)
    ax.hlines(
        lanes["Eurasian B"], chb_split_x, ooa_split_x, color="#666666", linewidth=2.2
    )
    ax.hlines(lanes["CEU"], 0, chb_split_x, color="#666666", linewidth=2.2)
    ax.hlines(lanes["Han / CHB"], 0, chb_split_x, color="#666666", linewidth=2.2)
    ax.vlines(
        chb_split_x,
        lanes["Han / CHB"],
        lanes["Eurasian B"],
        color="#777777",
        linewidth=1.3,
        linestyle=":",
    )
    ax.axvspan(migration_end_x, ooa_split_x, color="#7a5195", alpha=0.09)
    ax.annotate(
        "",
        xy=(0.87, 2.92),
        xytext=(0.87, 1.08),
        arrowprops={"arrowstyle": "<->", "linewidth": 1.6, "color": "#7a5195"},
    )
    ax.annotate(
        "",
        xy=(1.02, 2.92),
        xytext=(1.02, 0.15),
        arrowprops={
            "arrowstyle": "<->",
            "linewidth": 1.6,
            "color": "#7a5195",
            "connectionstyle": "arc3,rad=-0.10",
        },
    )
    ax.annotate(
        "",
        xy=(1.69, 2.92),
        xytext=(1.69, 2.08),
        arrowprops={"arrowstyle": "<->", "linewidth": 1.6, "color": "#7a5195"},
    )
    ax.text(
        1.48,
        3.38,
        "continuous symmetric Neanderthal-Eurasian migration\n"
        "$m=8.25\\times10^{-6}$/gen; recipients: CEU/CHB, then Eurasian B",
        color="#6a3d7c",
        fontsize=10.2,
        ha="center",
        va="center",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.86, "pad": 1.5},
    )
    ax.axvline(migration_end_x, color="#7a5195", linestyle="--", linewidth=1.8)
    ax.annotate(
        (
            "migration ends; selection begins\n"
            "operational CHB allele-AF screen\n"
            f"{100 * plan.han_selection_start_af_target:.1f}% "
            f"({100 * (plan.han_selection_start_af_target - plan.han_selection_start_af_tolerance):.1f}%–"
            f"{100 * (plan.han_selection_start_af_target + plan.han_selection_start_af_tolerance):.1f}%)"
        ),
        xy=(migration_end_x, 0.02),
        xytext=(0.75, -0.52),
        arrowprops={"arrowstyle": "->", "color": "#7a5195", "linewidth": 1.2},
        color="#6a3d7c",
        fontsize=10.2,
        ha="center",
        va="center",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.93, "pad": 1.7},
    )
    birth_x = plan.han_mutation_birth_generations / 1_000.0
    ax.scatter(
        [birth_x], [lanes["Neanderthal"]], marker="*", s=130, color="#7a5195", zorder=4
    )
    ax.annotate(
        f"donor mutation birth\n{plan.han_mutation_birth_generations:,} g",
        xy=(birth_x, lanes["Neanderthal"]),
        xytext=(2.48, 2.72),
        textcoords="data",
        fontsize=10.2,
        ha="right",
        va="center",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.9, "pad": 1.4},
    )
    ax.text(
        chb_split_x + 0.035, 1.48, "CEU-CHB split\n1,241 g", fontsize=10.2, ha="left"
    )
    ax.vlines(
        ooa_split_x,
        lanes["Eurasian B"] - 0.36,
        lanes["Eurasian B"] + 0.36,
        color="#777777",
        linestyle=":",
        linewidth=1.3,
    )
    ax.text(
        ooa_split_x - 0.035,
        2.30,
        "OOA split\n2,093 g",
        fontsize=10.2,
        ha="right",
        va="center",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.9, "pad": 1.2},
    )
    han_cells = cells[cells["sample_population"] == "Han"]
    y_offsets = _han_bar_offsets(len(han_cells))
    for index, row in enumerate(han_cells.itertuples(index=False)):
        coefficient = float(row.selection_coefficient)
        onset = float(row.selection_onset_generations_ago) / 1_000.0
        y = lanes["Han / CHB"] + float(y_offsets[index])
        ax.hlines(y, onset, 0, color=colors[coefficient], linewidth=4.0)
        ax.scatter(
            [onset],
            [y],
            marker=markers[coefficient],
            s=60,
            color=colors[coefficient],
            edgecolor="white",
            linewidth=0.8,
            zorder=4,
        )
    ax.text(
        2.50,
        0.55,
        "YRI, ArchaicAFR, and older splits\nretained in full event ledger",
        fontsize=9.4,
        ha="right",
        va="center",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.88, "pad": 1.2},
    )
    ax.set_yticks(list(lanes.values()), labels=list(lanes.keys()))
    ax.set_title("C. 5R19 Han/CHB demography (recent)")
    ax.set_xlabel("Thousands of generations before present (present = 0)")
    ax.grid(axis="x", alpha=0.2)

    ax = axes[1, 1]
    han_response = age_response[age_response["demography_id"] == HAN_DEMOGRAPHY_ID]
    for coefficient, group in han_response.groupby("selection_coefficient", sort=True):
        coefficient = float(coefficient)
        ax.plot(
            group["time_generations"] / 1_000.0,
            group["deterministic_final_af"],
            color=colors[coefficient],
            linewidth=2.2,
        )
        ax.scatter(
            0.0,
            float(
                group.loc[group["time_generations"].idxmin(), "deterministic_final_af"]
            ),
            s=65,
            marker=markers[coefficient],
            color=colors[coefficient],
            edgecolor="white",
            linewidth=0.8,
            zorder=4,
        )
        endpoint = float(
            group.loc[group["time_generations"].idxmin(), "deterministic_final_af"]
        )
        ax.annotate(
            f"{endpoint:.3f}",
            xy=(0.0, endpoint),
            xytext=(8, 0),
            textcoords="offset points",
            color=colors[coefficient],
            fontsize=10.5,
            va="center",
        )
    ax.axhspan(
        plan.present_af - plan.present_af_tolerance,
        plan.present_af + plan.present_af_tolerance,
        color="#666666",
        alpha=0.12,
    )
    ax.axhline(plan.present_af, color="#555555", linestyle="--", linewidth=1.3)
    ax.axvline(migration_end_x, color="#7a5195", linestyle="--", linewidth=1.5)
    ax.text(
        migration_end_x - 0.015,
        0.268,
        "migration ends / selection starts: 645 g",
        color="#6a3d7c",
        fontsize=10.2,
        ha="right",
        va="top",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.9, "pad": 1.3},
    )
    ax.set_xlim(0, 0.72)
    han_y_max = max(
        0.28,
        float(han_response["deterministic_final_af"].max()) * 1.12,
        (plan.present_af + plan.present_af_tolerance) * 1.18,
    )
    ax.set_ylim(0, min(1.02, han_y_max))
    ax.set_title("D. Han/CHB allele-frequency paths")
    ax.set_xlabel("Time before present (10$^3$ generations)")
    ax.set_ylabel("Deterministic allele frequency")
    ax.grid(alpha=0.25)

    handles = [
        Line2D(
            [0],
            [0],
            color=colors[coefficient],
            marker=markers[coefficient],
            linewidth=2.2,
            label=f"$s={coefficient:g}$",
        )
        for coefficient in coefficients
    ] + [
        Line2D(
            [0], [0], color="#555555", linestyle="--", linewidth=1.5, label="AF target"
        ),
        Line2D(
            [0],
            [0],
            color="#7a5195",
            linestyle="--",
            linewidth=1.8,
            marker="*",
            label="5R19 contact / donor origin",
        ),
    ]
    figure.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.937),
        ncol=len(handles),
        frameon=False,
    )
    figure.subplots_adjust(
        left=0.12, right=0.97, bottom=0.11, top=0.85, hspace=0.46, wspace=0.42
    )
    output_png.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(
        output_png,
        dpi=PNG_DPI,
        bbox_inches=None,
        metadata={"Software": "gamma_smc_aou.cosi2_framework"},
    )
    figure.savefig(
        output_pdf,
        bbox_inches=None,
        metadata={"Creator": "gamma_smc_aou.cosi2_framework"},
    )
    plt.close(figure)


def _readme(plan: FrameworkPlan, cells: pd.DataFrame) -> str:
    rows = []
    for row in cells.itertuples(index=False):
        rows.append(
            f"| `{row.cell_id}` | {row.selection_coefficient:g} | "
            f"{int(row.mutation_birth_generations_ago):,} | "
            f"{int(row.selection_onset_generations_ago):,} | "
            f"{row.deterministic_present_af_from_gate_center:.3f} | "
            f"{row.present_af_lower:.3f}-{row.present_af_upper:.3f} |"
        )
    table = "\n".join(rows)
    introgression_af_lower = (
        plan.han_selection_start_af_target - plan.han_selection_start_af_tolerance
    )
    introgression_af_upper = (
        plan.han_selection_start_af_target + plan.han_selection_start_af_tolerance
    )
    han_endpoint_summary = "; ".join(
        f"s={row.selection_coefficient:g}: p={row.deterministic_present_af_from_gate_center:.3f}"
        for row in cells[cells["sample_population"] == "Han"].itertuples(index=False)
    )
    return f"""# CoSi2 EAS and Han framework

This additive bundle renders parameter files for CoSi2 {COSI2_VERSION}, pinned
to commit `{COSI2_COMMIT}`. It does **not** contain production simulations.

## Canonical cells

| Cell | s | mutation birth (generations ago) | selection onset (generations ago) | deterministic endpoint from gate center | present population AF gate |
|---|---:|---:|---:|---:|---:|
{table}

The EAS mutation-birth ages are one-copy deterministic fixed-point centers on
the checked-in PHLASH median, not inferred ages. EAS uses 25 years/generation.
The Han/CHB model uses 29 years/generation, and the mutation is born as one copy
in Neanderthal at the explicitly configurable donor age. Selection starts in
CHB at generation {HAN_MIGRATION_END_GENERATIONS:,}, the integer realization of
18.7 kya, when 5R19's continuous archaic migration stops.

## Critical Han distinction

`{HAN_MODEL_ID}` does **not** contain a 1.2% pulse. It models symmetric,
continuous Neanderthal-Eurasian migration at {HAN_NEANDERTHAL_MIGRATION_RATE:g}
per generation from the 60.7-kya OOA split until 18.7 kya. The
Ragsdale-Gravel estimate of {100 * HAN_ARCHAIC_ANCESTRY_ESTIMATE:.1f}% +/-
{100 * HAN_ARCHAIC_ANCESTRY_UNCERTAINTY:.1f}% is a genome-wide ancestry result,
not a per-locus migration parameter. The simple symmetric two-state contact
expectation for the catalog interval is
{100 * HAN_SYMMETRIC_CONTACT_ANCESTRY_EXPECTATION:.3f}%, consistent with that
estimate.

Following the requested ascertainment, this framework separately applies an
operational CHB allele-AF screen at generation {HAN_MIGRATION_END_GENERATIONS:,} around
{100 * plan.han_selection_start_af_target:.1f}%
({100 * introgression_af_lower:.1f}%-{100 * introgression_af_upper:.1f}%). Native
`sweep_mult_standing` enforces only the present CHB population-AF range. Set
`COSI_SAVE_TRAJ`, then retain a Han output only when `selfreq_3` passes this
migration-end screen. This is a user-specified locus ascertainment rule, not a
claim that genome-wide ancestry and every introgressed allele have equal AF.
Use the CLI `validate-trajectory` phase before replaying an accepted trajectory.

The endpoint gate is a **population** AF interval. CoSi2 subsequently draws the
sampled carrier count binomially, so a {plan.sample_haploids}-haplotype sample is not constrained
to exactly {100 * plan.present_af:g}%.

CoSi2's `MSweep::makeSweepModel` constructs selected/unselected siblings for
every population in the base model and looks up a sampling request for each one.
The Han parameter file therefore emits explicit `sample_size 0` requests for
YRI, CEU/OOA, Neanderthal, and ArchaicAFR, while sampling
{plan.sample_haploids} haplotypes from CHB. The zero requests are required CoSi2
bookkeeping; they do not add sampled chromosomes.

The default 2,400-generation donor birth is a transparent age-grid starting
point, not a calibrated age. It predates the full Neanderthal-contact interval.
Scan donor birth ages or supply a separately validated trajectory after the
smoke phase. With the fixed 645-generation selection duration and configured
gate center p={plan.han_selection_start_af_target:g}, the deterministic endpoint
summary is {han_endpoint_summary}. Joint intermediate-plus-endpoint acceptance
can be extremely low. For the narrow 1.1%-1.3% default Han boundary gate, the
current feasibility estimate is roughly one joint acceptance per 28 million
trajectory draws; treat this as a planning estimate, not a calibrated rate.
`run_commands.sh` is a syntax smoke, does not implement an outer rejection loop,
and does not promise that its one Han draw will pass the boundary validator.

stdpopsim records `mutation_rate=None` for 5R19 because the model was calibrated
using recombination. The mutation rate and recombination map in these parameter
files are explicitly labeled study simulation overrides. The one-line constant
map starts at bp 1 because CoSi2 v2.4 rejects a zero start; its first rate covers
the region start and remains in force through the sequence end.

## Intended phases

1. `plan`: render and strictly validate this bundle; run no CoSi2 simulation.
2. Smoke: run one fixed-seed trajectory/coalescent replicate per cell from
   `run_commands.sh` and inspect stderr plus the saved trajectory.
3. Han gate: validate CHB AF at the migration-end/selection-start boundary.
4. Replay: strip `popsize_*` columns, validate every population-frequency
   column, and use `COSI_LOAD_TRAJ` for a later full-region replicate.
5. Scale only after smoke acceptance and runtime profiling. Keep raw ARG to
   tskit conversion as a separately tested downstream phase.

## Build CoSi2 and render this framework

```bash
git clone git@github.com:soisa001/gamma_smc_ts.git
cd gamma_smc_ts
git pull --ff-only origin AOU_run_opt
uv sync --frozen --extra selection
uv run python scripts/run_cosi2_framework.py plan

git clone --branch {COSI2_TAG} --depth 1 {COSI2_REPOSITORY} external/cosi2
test "$(git -C external/cosi2 rev-parse HEAD)" = "{COSI2_COMMIT}"
cd external/cosi2 && ./configure && make -j4 && cd ../..
```

Run commands are in `run_commands.sh`. They intentionally request one smoke
replicate per cell, use fixed seeds, and run from this directory because CoSi2
resolves `recomb_file` relative to the current working directory. All mutable
outputs are written outside this exact, immutable bundle under
`{RUN_ROOT_RELATIVE}/smoke/`; this lets `validate-trajectory` revalidate the
bundle inventory after a run. A smoke output is not accepted until the
trajectory and endpoint validators pass. Increase replicate counts only after
measuring the joint gate acceptance rate.

CoSi2 emits ms-format haplotypes and optional raw ARG edges, not native tskit
tree sequences. Any ARG-to-tskit bridge should be a separate validated phase.
"""


def _run_commands(cells: pd.DataFrame) -> str:
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "COALESCENT=${COSI2_COALESCENT:-../../external/cosi2/coalescent}",
        "MAXATTEMPTS=${COSI_MAXATTEMPTS:-10000}",
        'test -x "$COALESCENT"',
        f"RUN_ROOT={RUN_ROOT_RELATIVE}",
        'mkdir -p "$RUN_ROOT/smoke"',
        "# Syntax smoke only: the narrow Han boundary gate may reject the one draw.",
        "# This script does not implement the required outer rejection loop.",
    ]
    for row in cells.itertuples(index=False):
        lines.extend(
            [
                f'mkdir -p "$RUN_ROOT/smoke/{row.cell_id}"',
                (
                    f'COSI_MAXATTEMPTS="$MAXATTEMPTS" '
                    f'COSI_SAVE_TRAJ="$RUN_ROOT/smoke/{row.cell_id}/trajectory.tsv" '
                    f'"$COALESCENT" -p configs/{row.cell_id}.par -n 1 '
                    f"-r {int(row.seed)} -u 1 -m -M -e "
                    f'> "$RUN_ROOT/smoke/{row.cell_id}/simulation.ms" '
                    f'2> "$RUN_ROOT/smoke/{row.cell_id}/simulation.stderr.log"'
                ),
            ]
        )
        lines.append(
            "uv run python ../../scripts/run_cosi2_framework.py "
            f"--repo-root ../.. validate-trajectory --cell-id {row.cell_id} "
            f'--trajectory "$RUN_ROOT/smoke/{row.cell_id}/trajectory.tsv"'
        )
    return "\n".join(lines) + "\n"


def _output_records(root: Path, *, exclude: Iterable[str] = ()) -> dict[str, Any]:
    excluded = set(exclude)
    records: dict[str, Any] = {}
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        relative = path.relative_to(root).as_posix()
        if relative in excluded:
            continue
        records[relative] = {
            "path": relative,
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
    return records


def _framework_contract(repo_root: Path, plan: FrameworkPlan) -> dict[str, Any]:
    resource = repo_root / EAS_RESOURCE_PATH
    module = repo_root / MODULE_PATH
    wrapper = repo_root / WRAPPER_PATH
    return {
        "schema": SCHEMA,
        "plan": plan.to_record(),
        "cosi2": {
            "version_label": COSI2_VERSION,
            "tag": COSI2_TAG,
            "commit": COSI2_COMMIT,
            "repository": COSI2_REPOSITORY,
            "version_command_is_stale_in_tag": True,
        },
        "source_hashing": "strict_UTF8_canonical_LF",
        "sources": {
            MODULE_PATH: canonical_source_sha256(module),
            WRAPPER_PATH: canonical_source_sha256(wrapper),
        },
        "inputs": {
            EAS_RESOURCE_PATH: {
                "sha256": sha256_file(resource),
                "expected_sha256": EAS_RESOURCE_SHA256,
            },
            "stdpopsim_catalog_model": HAN_MODEL_ID,
            "stdpopsim_version": _package_version("stdpopsim"),
        },
        "runtime": {
            "python": sys.version.split()[0],
            "numpy": _package_version("numpy"),
            "pandas": _package_version("pandas"),
            "matplotlib": _package_version("matplotlib"),
        },
        "semantics": {
            "population_sizes": "explicit_half_up_integer_rounding",
            "present_af": "population_frequency_interval_not_sample_AF",
            "dominance": "CoSi2_hardcoded_additive_h_0.5",
            "han_introgression": (
                "5R19_continuous_symmetric_Neanderthal_migration; paper_1.2pct_"
                "is_genomewide_ancestry_estimate; operational_CHB_AF_screen_"
                "validated_posthoc_at_integer_generation_645"
            ),
            "han_operational_allele_af_screen": {
                "generation": HAN_MIGRATION_END_GENERATIONS,
                "target": plan.han_selection_start_af_target,
                "tolerance": plan.han_selection_start_af_tolerance,
                "lower": (
                    plan.han_selection_start_af_target
                    - plan.han_selection_start_af_tolerance
                ),
                "upper": (
                    plan.han_selection_start_af_target
                    + plan.han_selection_start_af_tolerance
                ),
                "is_catalog_parameter": False,
            },
            "han_catalog_generation_time_years": HAN_GENERATION_TIME_YEARS,
            "han_catalog_mutation_rate": None,
            "han_symmetric_contact_ancestry_expectation": (
                HAN_SYMMETRIC_CONTACT_ANCESTRY_EXPECTATION
            ),
            "study_mutation_rate_override": plan.mutation_rate,
            "study_recombination_map_override": plan.recombination_rate,
            "recombination_map_first_position_bp": 1,
            "recombination_map_coverage": (
                "first_rate_applies_from_region_start_through_sequence_end"
            ),
            "catalog_to_cosi2_time_policy": (
                "nearest_integer_half_up_with_catalog_and_rounding_delta_recorded"
            ),
            "selection_duration": "one_population_from_onset_to_present",
            "han_population_sampling_requests": (
                "explicit_zero_for_unsampled_5R19_populations_required_by_CoSi2_MSweep"
            ),
            "time_axis_orientation": "present_generation_0_on_left",
            "bundle_mutability": "immutable_exact_inventory",
            "run_artifact_root_relative_to_bundle": RUN_ROOT_RELATIVE,
        },
    }


def generate_framework(
    repo_root: str | Path,
    *,
    output_dir: str | Path | None = None,
    plan: FrameworkPlan | None = None,
) -> dict[str, Any]:
    """Atomically generate or verify the framework bundle."""

    root = Path(repo_root).resolve()
    resolved_plan = plan or FrameworkPlan()
    resolved_plan.validate()
    destination = (
        Path(output_dir).resolve()
        if output_dir is not None
        else root / DEFAULT_OUTPUT_RELATIVE
    )
    if destination.exists():
        try:
            result = validate_framework(root, output_dir=destination)
        except Exception as error:
            entries = (
                list(destination.iterdir()) if destination.is_dir() else [destination]
            )
            if entries:
                raise RuntimeError(
                    "refusing to overwrite a nonempty incompatible CoSi2 framework "
                    f"bundle: {destination}"
                ) from error
        else:
            expected = _framework_contract(root, resolved_plan)
            if result["contract_sha256"] != _canonical_sha256(expected):
                raise RuntimeError(
                    "cached framework plan differs from requested inputs"
                )
            return {**result, "status": "verified_cached"}
    staging = destination.with_name(destination.name + f".staging.{os.getpid()}")
    if staging.exists():
        raise RuntimeError(f"staging directory already exists: {staging}")
    staging.mkdir(parents=True)
    try:
        history = load_eas_history(root)
        rle = rle_eas_history(history)
        cells = build_model_cells(history, resolved_plan)
        age_response = build_age_response(history, resolved_plan)
        han_events = build_han_event_ledger(resolved_plan)
        eas_events = rle.iloc[1:][
            ["state_index", "time_generations", "raw_median_ne", "cosi2_integer_ne"]
        ].copy()
        eas_events.insert(0, "demography_id", "eas_phlash_median")
        eas_events.insert(1, "event_type", "change_size")
        eas_events.insert(
            2,
            "label",
            eas_events["state_index"].map(lambda value: f"eas_ne_{int(value):06d}"),
        )
        eas_events = eas_events.rename(
            columns={
                "time_generations": "catalog_time_generations",
                "cosi2_integer_ne": "integer_size",
            }
        )
        eas_events["cosi2_time_generations"] = eas_events["catalog_time_generations"]
        eas_events["catalog_time_years"] = (
            eas_events["catalog_time_generations"]
            * resolved_plan.eas_generation_time_years
        )
        eas_events["cosi2_nominal_integer_generation"] = np.floor(
            eas_events["cosi2_time_generations"] + 0.5
        ).astype(np.int64)
        eas_events["rounding_delta_generations"] = (
            eas_events["cosi2_time_generations"]
            - eas_events["catalog_time_generations"]
        )
        eas_events["rounding_delta_years"] = 0.0
        eas_events["rounding_policy"] = "exact_PHLASH_grid_time"
        eas_events["population_id"] = 1
        eas_events["other_population_id"] = np.nan
        eas_events["migration_rate"] = np.nan
        eas_events["proportion"] = np.nan
        eas_events["catalog_start_time_generations"] = np.nan
        eas_events["cosi2_start_time_generations"] = np.nan
        eas_events["integer_start_size"] = np.nan
        eas_events["source_model"] = "PHLASH_EAS_pointwise_median"
        han_events = han_events.copy()
        han_events["raw_median_ne"] = np.nan
        han_events["source_model"] = "stdpopsim_OutOfAfricaArchaicAdmixture_5R19"
        ledger_columns = [
            "demography_id",
            "source_model",
            "event_type",
            "label",
            "catalog_time_years",
            "catalog_time_generations",
            "cosi2_nominal_integer_generation",
            "cosi2_time_generations",
            "rounding_delta_generations",
            "rounding_delta_years",
            "rounding_policy",
            "population_id",
            "other_population_id",
            "migration_rate",
            "proportion",
            "raw_median_ne",
            "integer_size",
            "catalog_start_time_generations",
            "cosi2_start_time_generations",
            "integer_start_size",
        ]
        demography_events = pd.concat(
            [
                eas_events[ledger_columns],
                han_events[ledger_columns],
            ],
            ignore_index=True,
        )
        _atomic_frame(staging / "model_cells.tsv", cells)
        _atomic_frame(staging / "allele_age_response.tsv", age_response)
        _atomic_frame(staging / "demography_events.tsv", demography_events)
        _atomic_frame(staging / "eas_phlash_integer_history.tsv", history)
        _atomic_text(
            staging / "recombination_map.tsv",
            _recombination_map_text(resolved_plan),
        )
        for row in cells.to_dict(orient="records"):
            if row["sample_population"] == "EAS":
                text = render_eas_parameter_file(rle, row, resolved_plan)
            else:
                text = render_han_parameter_file(han_events, row, resolved_plan)
            lint_parameter_text(text, expected_cell=row)
            _atomic_text(staging / "configs" / f"{row['cell_id']}.par", text)
        _atomic_text(staging / "README.md", _readme(resolved_plan, cells))
        _atomic_text(staging / "run_commands.sh", _run_commands(cells))
        contract = _framework_contract(root, resolved_plan)
        plan_record = {
            "contract": contract,
            "contract_sha256": _canonical_sha256(contract),
            "cells": (
                cells.astype(object)
                .where(pd.notna(cells), None)
                .to_dict(orient="records")
            ),
            "eas_history": {
                "raw_rows": len(history),
                "integer_rle_states": len(rle),
                "change_size_directives": len(rle) - 1,
                "present_integer_ne": int(rle.iloc[0]["cosi2_integer_ne"]),
                "ancient_terminal_integer_ne": int(rle.iloc[-1]["cosi2_integer_ne"]),
            },
            "han_5r19": {
                "model_id": HAN_MODEL_ID,
                "generation_time_years": HAN_GENERATION_TIME_YEARS,
                "neanderthal_migration_rate_per_generation": (
                    HAN_NEANDERTHAL_MIGRATION_RATE
                ),
                "migration_end_catalog_generations": (
                    HAN_MIGRATION_END_CATALOG_GENERATIONS
                ),
                "migration_end_cosi2_generation": HAN_MIGRATION_END_GENERATIONS,
                "paper_genomewide_ancestry_estimate": (HAN_ARCHAIC_ANCESTRY_ESTIMATE),
                "paper_genomewide_ancestry_uncertainty": (
                    HAN_ARCHAIC_ANCESTRY_UNCERTAINTY
                ),
                "symmetric_contact_ancestry_expectation": (
                    HAN_SYMMETRIC_CONTACT_ANCESTRY_EXPECTATION
                ),
                "operational_allele_gate_is_user_ascertainment": True,
            },
        }
        _atomic_text(
            staging / "framework_plan.json",
            json.dumps(plan_record, sort_keys=True, indent=2, allow_nan=False) + "\n",
        )
        _plot_framework(
            staging / "figures" / "cosi2_demography_specification.png",
            staging / "figures" / "cosi2_demography_specification.pdf",
            history=history,
            age_response=age_response,
            cells=cells,
            plan=resolved_plan,
        )
        outputs = _output_records(staging)
        completion = {
            "schema": COMPLETION_SCHEMA,
            "status": "complete",
            "contract": contract,
            "contract_sha256": _canonical_sha256(contract),
            "outputs": outputs,
            "output_count": len(outputs),
        }
        _atomic_text(
            staging / "framework_completion.json",
            json.dumps(completion, sort_keys=True, indent=2, allow_nan=False) + "\n",
        )
        if destination.exists():
            if destination.is_dir() and not any(destination.iterdir()):
                destination.rmdir()
            else:
                raise RuntimeError("destination became nonempty during generation")
        destination.parent.mkdir(parents=True, exist_ok=True)
        os.replace(staging, destination)
    finally:
        if staging.exists():
            shutil.rmtree(staging)
    return validate_framework(root, output_dir=destination)


def _png_dimensions(path: Path) -> tuple[int, int]:
    payload = path.read_bytes()
    if payload[:8] != b"\x89PNG\r\n\x1a\n" or payload[12:16] != b"IHDR":
        raise ValueError("invalid PNG signature")
    return int.from_bytes(payload[16:20], "big"), int.from_bytes(payload[20:24], "big")


def _pdf_letter_geometry(path: Path) -> tuple[int, float, float]:
    payload = path.read_bytes()
    if not payload.startswith(b"%PDF-"):
        raise ValueError("invalid PDF signature")
    page_count = len(re.findall(rb"/Type\s*/Page(?!s)\b", payload))
    boxes = re.findall(
        rb"/MediaBox\s*\[\s*0(?:\.0+)?\s+0(?:\.0+)?\s+"
        rb"([0-9.]+)\s+([0-9.]+)\s*\]",
        payload,
    )
    if page_count != 1 or not boxes:
        raise ValueError("framework PDF must contain exactly one parseable page")
    width = float(boxes[0][0])
    height = float(boxes[0][1])
    if not math.isclose(width, 792.0, abs_tol=0.01) or not math.isclose(
        height, 612.0, abs_tol=0.01
    ):
        raise ValueError(f"framework PDF is not landscape Letter: {width} x {height}")
    return page_count, width, height


def validate_framework(
    repo_root: str | Path,
    *,
    output_dir: str | Path | None = None,
) -> dict[str, Any]:
    root = Path(repo_root).resolve()
    bundle = (
        Path(output_dir).resolve()
        if output_dir is not None
        else root / DEFAULT_OUTPUT_RELATIVE
    )
    completion_path = bundle / "framework_completion.json"
    if not completion_path.is_file() or completion_path.is_symlink():
        raise ValueError("framework completion is absent or not a regular file")
    completion = json.loads(completion_path.read_text(encoding="utf-8"))
    if set(completion) != {
        "schema",
        "status",
        "contract",
        "contract_sha256",
        "outputs",
        "output_count",
    }:
        raise ValueError("framework completion field inventory changed")
    if completion["schema"] != COMPLETION_SCHEMA or completion["status"] != "complete":
        raise ValueError("framework completion schema/status is invalid")
    resolved_plan = FrameworkPlan(**completion["contract"]["plan"])
    current_contract = _framework_contract(root, resolved_plan)
    current_sha = _canonical_sha256(current_contract)
    if (
        completion["contract"] != current_contract
        or completion["contract_sha256"] != current_sha
    ):
        raise ValueError("framework completion contract is stale")
    outputs = completion["outputs"]
    if not isinstance(outputs, dict) or completion["output_count"] != len(outputs):
        raise ValueError("framework output manifest is invalid")
    coefficients = {float(value) for value in resolved_plan.selection_coefficients}
    expected_cell_ids = {
        f"{prefix}_s{str(coefficient).replace('.', 'p')}"
        for coefficient in coefficients
        for prefix in ("eas", "han_introgressed")
    }
    expected_output_labels = {
        "README.md",
        "allele_age_response.tsv",
        "demography_events.tsv",
        "eas_phlash_integer_history.tsv",
        "figures/cosi2_demography_specification.pdf",
        "figures/cosi2_demography_specification.png",
        "framework_plan.json",
        "model_cells.tsv",
        "recombination_map.tsv",
        "run_commands.sh",
        *(f"configs/{cell_id}.par" for cell_id in expected_cell_ids),
    }
    if set(outputs) != expected_output_labels:
        raise ValueError("framework output manifest label set is invalid")
    allowed_directories = {"configs", "figures"}
    for path in bundle.rglob("*"):
        relative = path.relative_to(bundle).as_posix()
        if path.is_symlink():
            raise ValueError(f"framework bundle contains a symlink: {relative}")
        if path.is_dir() and relative not in allowed_directories:
            raise ValueError(
                f"framework bundle contains an unexpected directory: {relative}"
            )
    actual = {
        path.relative_to(bundle).as_posix()
        for path in bundle.rglob("*")
        if path.is_file()
    }
    expected = set(outputs) | {"framework_completion.json"}
    if actual != expected:
        raise ValueError(
            f"framework disk inventory mismatch: missing={sorted(expected - actual)}, "
            f"extra={sorted(actual - expected)}"
        )
    for label, record in outputs.items():
        if set(record) != {"path", "bytes", "sha256"} or record["path"] != label:
            raise ValueError(f"invalid output record {label}")
        path = bundle / label
        if path.is_symlink() or not path.is_file():
            raise ValueError(f"output is not a regular file: {label}")
        if (
            path.stat().st_size != record["bytes"]
            or sha256_file(path) != record["sha256"]
        ):
            raise ValueError(f"output hash/size mismatch: {label}")
    recombination_map = (bundle / "recombination_map.tsv").read_text(encoding="utf-8")
    if recombination_map != _recombination_map_text(resolved_plan):
        raise ValueError("framework recombination map content is invalid")
    cells = pd.read_csv(bundle / "model_cells.tsv", sep="\t")
    expected_cells = 2 * len(coefficients)
    observed_pairs = set(
        zip(
            cells["sample_population"],
            cells["selection_coefficient"].astype(float),
            strict=True,
        )
    )
    expected_pairs = {
        (population, coefficient)
        for population in ("EAS", "Han")
        for coefficient in coefficients
    }
    if (
        len(cells) != expected_cells
        or observed_pairs != expected_pairs
        or set(cells["cell_id"]) != expected_cell_ids
    ):
        raise ValueError("model cell table is incomplete")
    for row in cells.to_dict(orient="records"):
        text = (bundle / "configs" / f"{row['cell_id']}.par").read_text(
            encoding="utf-8"
        )
        lint_parameter_text(text, expected_cell=row)
    png = bundle / "figures" / "cosi2_demography_specification.png"
    pdf = bundle / "figures" / "cosi2_demography_specification.pdf"
    if _png_dimensions(png) != (3_300, 2_550):
        raise ValueError("framework PNG must be exactly 3300 x 2550")
    pages, width, height = _pdf_letter_geometry(pdf)
    return {
        "status": "complete",
        "bundle": str(bundle),
        "contract_sha256": current_sha,
        "output_count": len(outputs),
        "cell_count": len(cells),
        "png_dimensions": [3_300, 2_550],
        "pdf_pages": pages,
        "pdf_points": [width, height],
    }


def _plan_from_args(args: argparse.Namespace) -> FrameworkPlan:
    coefficients = tuple(float(value) for value in args.selection_coefficients)
    return FrameworkPlan(
        selection_coefficients=coefficients,
        present_af=float(args.present_af),
        present_af_tolerance=float(args.present_af_tolerance),
        sample_haploids=int(args.sample_haploids),
        han_mutation_birth_generations=int(args.han_mutation_birth_generations),
        han_selection_start_af_target=float(args.han_selection_start_af_target),
        han_selection_start_af_tolerance=float(args.han_selection_start_af_tolerance),
        base_seed=int(args.base_seed),
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--output-dir", type=Path)
    subparsers = parser.add_subparsers(dest="phase", required=True)
    plan = subparsers.add_parser("plan", help="atomically render configs and figures")
    plan.add_argument(
        "--selection-coefficients",
        nargs="+",
        type=float,
        default=list(DEFAULT_SELECTION_COEFFICIENTS),
    )
    plan.add_argument("--present-af", type=float, default=DEFAULT_PRESENT_AF)
    plan.add_argument(
        "--present-af-tolerance", type=float, default=DEFAULT_PRESENT_AF_TOLERANCE
    )
    plan.add_argument("--sample-haploids", type=int, default=DEFAULT_SAMPLE_HAPLOIDS)
    plan.add_argument(
        "--han-mutation-birth-generations",
        type=int,
        default=DEFAULT_HAN_MUTATION_BIRTH_GENERATIONS,
    )
    plan.add_argument(
        "--han-selection-start-af-target",
        type=float,
        default=DEFAULT_HAN_SELECTION_START_AF_TARGET,
    )
    plan.add_argument(
        "--han-selection-start-af-tolerance",
        type=float,
        default=DEFAULT_HAN_SELECTION_START_AF_TOLERANCE,
    )
    plan.add_argument("--base-seed", type=int, default=BASE_SEED)
    subparsers.add_parser("validate", help="read-only strict bundle validation")
    trajectory = subparsers.add_parser(
        "validate-trajectory", help="validate one saved modern CoSi2 trajectory"
    )
    trajectory.add_argument("--trajectory", type=Path, required=True)
    trajectory.add_argument("--cell-id", required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    repo_root = args.repo_root.resolve()
    output_dir = args.output_dir
    if args.phase == "plan":
        result = generate_framework(
            repo_root,
            output_dir=output_dir,
            plan=_plan_from_args(args),
        )
    elif args.phase == "validate":
        result = validate_framework(repo_root, output_dir=output_dir)
    else:
        bundle = (
            output_dir.resolve()
            if output_dir is not None
            else repo_root / DEFAULT_OUTPUT_RELATIVE
        )
        validate_framework(repo_root, output_dir=bundle)
        cells = pd.read_csv(bundle / "model_cells.tsv", sep="\t")
        matches = cells[cells["cell_id"].astype(str) == str(args.cell_id)]
        if len(matches) != 1:
            raise ValueError(f"unknown or duplicate cell id: {args.cell_id}")
        result = validate_saved_trajectory(
            args.trajectory,
            cell=matches.iloc[0].to_dict(),
        )
    print(json.dumps(result, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
