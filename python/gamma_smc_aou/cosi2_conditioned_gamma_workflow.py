"""Restartable conditioned-neutral CoSi2 to Gamma-SMC study workflow.

This additive workflow compares the two accepted ``s=0.01`` selected loci
with two deliberately different neutral reference banks:

* EAS: 100 neutral loci conditioned on the exact selected present population
  count (14,704 / 75,274), with mutation ages sampled with replacement from
  the fixed global importance-weighted reverse-time proposal bank implemented by
  :mod:`cosi2_conditioned_null`.
* Han: 100 neutral Neanderthal-origin loci in the fixed 5R19 demography,
  gated to CHB AF 0.011--0.013 at generation 645 and required only to remain
  segregating in CHB at the present.

Simulation, decoding, and reporting artifacts live in roots separate from the
completed natural-neutral study.  Every durable directory is checksum-bound
and promoted from an adjacent staging directory only after validation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time
import platform
from collections.abc import Mapping, Sequence
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np
import pandas as pd
import scipy
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from . import cosi2_conditioned_null as conditioned
from . import cosi2_gamma_analysis as analysis
from . import cosi2_gamma_decode as gamma_decode
from . import cosi2_gamma_workflow as natural_workflow
from . import cosi2_neutral_simulation as natural_neutral

matplotlib.use("Agg")


MODULE_PATH = Path(__file__).resolve()

SCHEMA_VERSION = "gamma-smc.cosi2-conditioned-gamma-workflow/v1"
EAS_BANK_WORKFLOW_SCHEMA = "gamma-smc.cosi2-conditioned-gamma-eas-bank/v2"
PLAN_SCHEMA = "gamma-smc.cosi2-conditioned-gamma-plan/v7"
LEGACY_V6_PLAN_SCHEMA = "gamma-smc.cosi2-conditioned-gamma-plan/v6"
LEGACY_V5_PLAN_SCHEMA = "gamma-smc.cosi2-conditioned-gamma-plan/v5"
LEGACY_V4_PLAN_SCHEMA = "gamma-smc.cosi2-conditioned-gamma-plan/v4"
LEGACY_V3_PLAN_SCHEMA = "gamma-smc.cosi2-conditioned-gamma-plan/v3"
LEGACY_V2_PLAN_SCHEMA = "gamma-smc.cosi2-conditioned-gamma-plan/v2"
EAS_PARALLEL_REPLAY_PLAN_SCHEMA = "gamma-smc.cosi2-eas-parallel-replay-plan/v5"
LEGACY_V4_PARALLEL_REPLAY_PLAN_SCHEMA = "gamma-smc.cosi2-eas-parallel-replay-plan/v4"
LEGACY_V3_PARALLEL_REPLAY_PLAN_SCHEMA = "gamma-smc.cosi2-eas-parallel-replay-plan/v3"
EAS_PARALLEL_REPLAY_BLOCK_SCHEMA = "gamma-smc.cosi2-eas-parallel-replay-block/v5"
LEGACY_V4_PARALLEL_REPLAY_BLOCK_SCHEMA = "gamma-smc.cosi2-eas-parallel-replay-block/v4"
LEGACY_V3_PARALLEL_REPLAY_BLOCK_SCHEMA = "gamma-smc.cosi2-eas-parallel-replay-block/v3"
EAS_PARALLEL_REPLAY_COMPLETION_SCHEMA = (
    "gamma-smc.cosi2-eas-parallel-replay-completion/v5"
)
LEGACY_V4_PARALLEL_REPLAY_COMPLETION_SCHEMA = (
    "gamma-smc.cosi2-eas-parallel-replay-completion/v4"
)
EAS_SELECTION_SCHEMA = "gamma-smc.cosi2-eas-global-selection/v5"
RESULTS_SCHEMA = "gamma-smc.cosi2-conditioned-gamma-results/v1"

DEFAULT_WORK_RELATIVE = Path(
    "focused_selection_EAS_sim/cosi2_conditioned_gamma_study_work"
)
DEFAULT_RESULTS_RELATIVE = Path(
    "focused_selection_EAS_sim/results/cosi2_gamma_conditioned_s0p01"
)
DEFAULT_SELECTED_RELATIVE = natural_workflow.DEFAULT_SELECTED_RELATIVE
DEFAULT_NATURAL_WORK_RELATIVE = natural_workflow.DEFAULT_WORK_RELATIVE
DEFAULT_COSI2_BINARY_RELATIVE = natural_workflow.DEFAULT_COSI2_BINARY_RELATIVE
DEFAULT_GAMMA_BINARY_RELATIVE = natural_workflow.DEFAULT_GAMMA_BINARY_RELATIVE

# Preserve the completed million-proposal pilot plan byte-for-byte.  The v2
# plan is a sibling bundle whose fixed 10-million ledger and base-bank
# attestation must exist before any extension block is generated.
LEGACY_PLAN_DIRNAME = "conditioned_plan"
LEGACY_V2_PLAN_DIRNAME = "conditioned_plan_v2"
LEGACY_V3_PLAN_DIRNAME = "conditioned_plan_v3"
LEGACY_V4_PLAN_DIRNAME = "conditioned_plan_v4"
LEGACY_V5_PLAN_DIRNAME = "conditioned_plan_v5"
LEGACY_V6_PLAN_DIRNAME = "conditioned_plan_v6"
PLAN_DIRNAME = "conditioned_plan_v7"
EAS_PROPOSAL_DIRNAME = "eas_proposal_blocks"
LEGACY_V3_PARALLEL_REPLAY_DIRNAME = "eas_parallel_replays_v3"
LEGACY_V4_PARALLEL_REPLAY_DIRNAME = "eas_parallel_replays_v4"
LEGACY_V5_PARALLEL_REPLAY_DIRNAME = "eas_parallel_replays_v5"
EAS_PARALLEL_REPLAY_DIRNAME = LEGACY_V5_PARALLEL_REPLAY_DIRNAME
LEGACY_V4_SELECTION_DIRNAME = "eas_global_selection"
LEGACY_V5_SELECTION_DIRNAME = "eas_global_selection_v5"
EAS_SELECTION_DIRNAME = LEGACY_V5_SELECTION_DIRNAME
HAN_SCREEN_DIRNAME = "han_screen"
LEGACY_V6_SIMULATION_DIRNAME = "conditioned_simulations"
SIMULATION_DIRNAME = "conditioned_simulations_v7"
SIMULATION_FAILURE_DIRNAME = "conditioned_simulation_failures_v7"
DECODE_DIRNAME = "decode_units"

PLAN_COMPLETION_FILENAME = "PLAN_COMPLETION.json"
EAS_SELECTION_COMPLETION = "selection_completion.json"
RESULTS_COMPLETION_FILENAME = "RESULTS_COMPLETION.json"

MAX_WORKERS = 20
LEGACY_V2_WORKFLOW_SHA256 = (
    "a905e9aa7a830e81f1af96556ad10f065879cbafac3153efb7aa10fcf5b06be1"
)
LEGACY_V2_PLAN_COMPLETION_SHA256 = (
    "c328a7bed54704e55dd668e8a8713e1b3e2573a46ed0465dd3e7e69bef024d63"
)
LEGACY_V3_WORKFLOW_SHA256 = (
    "827e06c636ab0aea9cc33645fa2570c7ae3e7dcf4abf30fe604150887c837b0d"
)
LEGACY_V3_PLAN_COMPLETION_SHA256 = (
    "2cb35c370528349b4da47f74ac1da03d1dcdd06096e8d4c6b0d6e56c84f76d68"
)
LEGACY_V3_REPLAY_PLAN_FILE_SHA256 = (
    "aadd2d1b0681a045b47796dc9cf8a23eb250f5919b2b39d24af9e4b618a399a2"
)
LEGACY_V3_REPLAY_PLAN_PAYLOAD_SHA256 = (
    "15c5c8103bb8293cc3cd8b3e0bc83bfc209609645fbde2b8853ea86cde421f3e"
)
LEGACY_V3_REPLAY_BLOCK_INDEXES = (5, 9, 11, 18, 19, 21, 24, 25, 34, 39, 41, 46, 53)
LEGACY_V3_REPLAY_MANIFEST_SHA256 = (
    "bb334be92318d9b87fb380184602e2eac5ee9a29941f6720e3b436334dfbcbf3"
)
LEGACY_V4_WORKFLOW_SHA256 = (
    "5e08e4990951e8b0c7c9fda17ddf2edb86875f637c2710e5c921cc0a88f5baf4"
)
LEGACY_V4_PLAN_COMPLETION_SHA256 = (
    "8093cca6336790dbd6879c384f0f312e515fc3cd208d980e04706f8738ca2ce3"
)
LEGACY_V4_REPLAY_PLAN_FILE_SHA256 = (
    "22e7fd16e38a80828054b6a0104ac282d5efe89d751f29fc75374c925bcfc88e"
)
LEGACY_V4_REPLAY_PLAN_PAYLOAD_SHA256 = (
    "7e744620422dbfaf4503edb219e8db4caded946d620b57c82a0d96b3afda2098"
)
LEGACY_V4_REPLAY_COMPLETION_SHA256 = (
    "2c6c692c7dd33618eb03055e0253372cc3ad1f2453b8c91f96e8ddb5235e4c85"
)
LEGACY_V4_REPLAY_ROOT_MANIFEST_SHA256 = (
    "b01007a3a6175fafd327380dbde5b9e1af5c29247994e18cead4848476acb4b1"
)
LEGACY_V5_WORKFLOW_SHA256 = (
    "8adfbff74ead615ee555693bb716e2374191cfad62f9b06b9304fd8e7e9048f1"
)
LEGACY_V5_PLAN_COMPLETION_SHA256 = (
    "9c9dc43e972f9c93b0b489a237e8f5aa27aaaadb6c6ec65d144f12009bfc4c1d"
)
LEGACY_V5_PLAN_CONTRACT_SHA256 = (
    "08826f51ea29f72a5d89627412afe903c09d43c031c6d8c4736cf50dc15aece7"
)
LEGACY_V5_REPLAY_PLAN_FILE_SHA256 = (
    "55de0f08320fa2a2cf6dd4d76a917b87c24e3bad8ce677e2023ede0b16fcac79"
)
LEGACY_V5_REPLAY_PLAN_PAYLOAD_SHA256 = (
    "8c12502b2244babd204f1f8f0bc77133d97230026025ea13e56a7cbe99e6ed48"
)
LEGACY_V5_SELECTION_COMPLETION_SHA256 = (
    "d0d2c9a394f2fb8e9d6c53d9ea4d79f66e90d93fc76fab57ddbfeea21ea39836"
)
LEGACY_V5_SELECTION_CONTRACT_SHA256 = (
    "c2ce251078ea3721ebd6d93c4e1e55becc313c47a5b4d0e765092651ceac939a"
)
LEGACY_V5_SELECTION_ROOT_FILE_COUNT = 503
LEGACY_V5_SELECTION_ROOT_MANIFEST_SHA256 = (
    "d9f8bd496c6eb89f2d662efbb47e7dd5d6d68940e9421c0034a544f8045963d1"
)
LEGACY_V6_WORKFLOW_SHA256 = (
    "382f7e8f88d20c6b44ed98633ecb0092ca14817d8bbc79e100d04a6d374815dd"
)
LEGACY_V6_PLAN_COMPLETION_SHA256 = (
    "43f11a3cdcfb13f64f3fc3da09bde751444f695ad8f999d5172c274aab88ff9b"
)
LEGACY_V6_PLAN_CONTRACT_SHA256 = (
    "62d73deaaefd2797039bfa7f05af4b78ed881cae5b9613cc67b305c595229552"
)
LEGACY_V6_HAN_ACCEPTED_SHA256 = (
    "1d726454bc4578bbf2c4701267d26b823399ea418b71f4f504a3d4a9c00f1dd3"
)
LEGACY_V6_HAN_SELECTION_COMPLETION_SHA256 = (
    "d8885305ebd723d2401ad366c83d646d669b06abc551a6e8cbea8b547faa953f"
)
LEGACY_V6_HAN_SELECTION_CONTRACT_SHA256 = (
    "f289249471011c060be46483ec203dfbcd2028db4ca3988c8e34e7b2673c2cdc"
)
LEGACY_V6_EAS_SIMULATION_FILE_COUNT = 500
LEGACY_V6_EAS_SIMULATION_ROOT_MANIFEST_SHA256 = (
    "32965ffd5cd3be6cf2bf24c5a2b3f95b0f583c040b736ca0ee0c0177874c6fce"
)
HAN_LOADER_TRAJECTORY_SCHEMA = "gamma-smc.cosi2-han-loader-trajectory/v1"
SIMULATION_FAILURE_SCHEMA = "gamma-smc.cosi2-conditioned-simulation-failure/v1"
SIMULATION_UNIT_SCHEMA = "gamma-smc.cosi2-conditioned-null-unit/v7"
SIMULATION_ROOT_SCHEMA = "gamma-smc.cosi2-conditioned-simulation-root/v7"
SIMULATION_ROOT_COMPLETION_FILENAME = "SIMULATION_COMPLETION.json"
RECOMBINATION_MAP_RELATIVE = Path(
    "focused_selection_EAS_sim/cosi2_framework/recombination_map.tsv"
)
FRAMEWORK_COMPLETION_RELATIVE = Path(
    "focused_selection_EAS_sim/cosi2_framework/framework_completion.json"
)
FRAMEWORK_COMPLETION_SHA256 = (
    "e2108ec5ba4ad415d84181e97c3eda6e6ad774fbf873005e0265434700a4f9d5"
)
RECOMBINATION_MAP_SIZE_BYTES = 8
RECOMBINATION_MAP_SHA256 = (
    "a1411e5610b5d40ce39fa2b6cd2cc63d1199dab0b6cd85a7c83ac859171953eb"
)
LEGACY_V4_PARALLEL_REPLAY_METHOD = (
    "one_process_task_per_selected_unique_block_core_exact_replay_validation"
)
EAS_PARALLEL_REPLAY_METHOD = "atomic_exact_v4_replay_root_migration_no_replay"
EAS_PARALLEL_REPLAY_COMPLETION_FILENAME = "REPLAY_COMPLETION.json"
EXPECTED_NEUTRAL_PER_DEMOGRAPHY = 100
EXPECTED_SELECTED_PER_DEMOGRAPHY = 1
EXPECTED_UNITS_PER_DEMOGRAPHY = 101
EXPECTED_TOTAL_UNITS = 202
SELECTION_COEFFICIENT = 0.01

DEMOGRAPHY_SETTINGS: dict[str, dict[str, Any]] = {
    "EAS": {
        **natural_workflow.DEMOGRAPHY_SETTINGS["EAS"],
        "neutral_model_id": conditioned.EAS_MODEL_ID,
    },
    "Han": {
        **natural_workflow.DEMOGRAPHY_SETTINGS["Han"],
        "neutral_model_id": conditioned.HAN_MODEL_ID,
    },
}

UNIT_COLUMNS = (
    *natural_workflow.UNIT_COLUMNS,
    "conditioning",
    "birth_generation",
    "allele_age_years",
    "eas_selected_log_importance_weight",
    "eas_global_effective_sample_size",
    "eas_global_weight_sum_over_max",
    "eas_global_max_normalized_weight",
    "eas_time_reversal_is_approximation",
    "eas_source_particle_id",
    "eas_draw_occurrence",
    "eas_particle_multiplicity",
    "han_candidate_id",
    "han_candidate_order",
    "han_candidate_source_kind",
    "han_generation_645_chb_af",
)

PROFILE_FIGURE_STEM = "cosi2_gamma_conditioned_selected_vs_neutral_profiles"
CONDITIONING_FIGURE_STEM = "cosi2_conditioning_diagnostics"
RESULT_OUTPUTS = {
    "scores": "cosi2_gamma_conditioned_scores.tsv",
    "conditioned_metadata": "conditioned_neutral_metadata.tsv",
    "unit_metadata": "unit_metadata.tsv",
    "eas_diagnostics": "eas_age_importance_diagnostics.tsv",
    "han_distribution": "han_gate_final_af_distribution.tsv",
    "conditioning_contrasts": "conditioning_simulation_pvalues.tsv",
    "pointwise_pvalues": "pointwise_simulation_pvalues.tsv",
    "minp_omnibus": "minp_omnibus.tsv",
    "run_results": "RUN_RESULTS.md",
    "profile_png": f"{PROFILE_FIGURE_STEM}.png",
    "profile_pdf": f"{PROFILE_FIGURE_STEM}.pdf",
    "conditioning_png": f"{CONDITIONING_FIGURE_STEM}.png",
    "conditioning_pdf": f"{CONDITIONING_FIGURE_STEM}.pdf",
}

_SHA256_RE = natural_workflow._SHA256_RE


def sha256_file(path: str | Path) -> str:
    """Return the streaming SHA-256 digest for a regular file."""

    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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


def _atomic_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.tmp.", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(value)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_json(path: Path, value: Mapping[str, Any]) -> None:
    _atomic_text(
        path,
        json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n",
    )


def _atomic_frame(path: Path, frame: pd.DataFrame) -> None:
    _atomic_text(
        path,
        frame.to_csv(
            sep="\t",
            index=False,
            lineterminator="\n",
            float_format="%.17g",
        ),
    )


def _read_json(path: Path, label: str) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"{label} is absent: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"{label} is unreadable: {path}") from error
    if not isinstance(value, dict):
        raise TypeError(f"{label} must contain a JSON object")
    return value


def _portable_frame_sha256(frame: pd.DataFrame) -> str:
    text = frame.to_csv(
        sep="\t",
        index=False,
        lineterminator="\n",
        float_format="%.17g",
    )
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _portable_record(record: Mapping[str, Any]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in record.items():
        if isinstance(value, np.generic):
            value = value.item()
        if value is None or (isinstance(value, float) and math.isnan(value)):
            result[str(key)] = None
        elif isinstance(value, Path):
            result[str(key)] = str(value)
        else:
            result[str(key)] = value
    return result


def _strict_boolean(series: pd.Series, label: str) -> pd.Series:
    def convert(value: Any) -> bool:
        if isinstance(value, (bool, np.bool_)):
            return bool(value)
        if isinstance(value, str):
            token = value.strip().casefold()
            if token == "true":
                return True
            if token == "false":
                return False
        raise ValueError(f"{label} must contain exact booleans")

    return series.map(convert)


def _output_record(
    path: Path, root: Path, *, rows: int | None = None
) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"output is absent or is a symlink: {path}")
    record: dict[str, Any] = {
        "path": path.relative_to(root).as_posix(),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }
    if rows is not None:
        record["rows"] = int(rows)
    return record


def _positive_workers(value: Any) -> int:
    if isinstance(value, bool):
        raise ValueError(f"workers must be an integer in [1, {MAX_WORKERS}]")
    workers = int(value)
    if workers != value or not 1 <= workers <= MAX_WORKERS:
        raise ValueError(f"workers must be an integer in [1, {MAX_WORKERS}]")
    return workers


def _relative(path: Path, root: Path) -> str:
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return str(path.resolve())


def _resolve_paths(
    repo_root: Path,
    work_dir: str | Path | None,
    results_dir: str | Path | None,
) -> dict[str, Path]:
    work = (
        Path(work_dir).resolve()
        if work_dir is not None
        else repo_root / DEFAULT_WORK_RELATIVE
    )
    results = (
        Path(results_dir).resolve()
        if results_dir is not None
        else repo_root / DEFAULT_RESULTS_RELATIVE
    )
    return {
        "work": work,
        "plan": work / PLAN_DIRNAME,
        "legacy_v2_plan": work / LEGACY_V2_PLAN_DIRNAME,
        "legacy_v3_plan": work / LEGACY_V3_PLAN_DIRNAME,
        "legacy_v4_plan": work / LEGACY_V4_PLAN_DIRNAME,
        "legacy_v5_plan": work / LEGACY_V5_PLAN_DIRNAME,
        "legacy_v6_plan": work / LEGACY_V6_PLAN_DIRNAME,
        "eas_proposals": work / EAS_PROPOSAL_DIRNAME,
        "legacy_v3_eas_parallel_replays": (work / LEGACY_V3_PARALLEL_REPLAY_DIRNAME),
        "legacy_v4_eas_parallel_replays": (work / LEGACY_V4_PARALLEL_REPLAY_DIRNAME),
        "eas_parallel_replays": work / EAS_PARALLEL_REPLAY_DIRNAME,
        "eas_selection": work / EAS_SELECTION_DIRNAME,
        "han_screen": work / HAN_SCREEN_DIRNAME,
        "legacy_v6_simulations": work / LEGACY_V6_SIMULATION_DIRNAME,
        "simulations": work / SIMULATION_DIRNAME,
        "simulation_failures": work / SIMULATION_FAILURE_DIRNAME,
        "decode": work / DECODE_DIRNAME,
        "results": results,
    }


def _selected_with_conditioning_columns(selected: pd.DataFrame) -> pd.DataFrame:
    result = selected.copy()
    result["conditioning"] = "selected_endpoint_population_af_0p195_to_0p205"
    result["birth_generation"] = result["demography"].map(
        {"EAS": 1_475, "Han": conditioned.HAN_BIRTH_GENERATION}
    )
    result["allele_age_years"] = (
        result["birth_generation"] * result["generation_time_years"]
    )
    result["eas_selected_log_importance_weight"] = np.nan
    result["eas_global_effective_sample_size"] = np.nan
    result["eas_global_weight_sum_over_max"] = np.nan
    result["eas_global_max_normalized_weight"] = np.nan
    result["eas_time_reversal_is_approximation"] = False
    result["eas_source_particle_id"] = ""
    result["eas_draw_occurrence"] = -1
    result["eas_particle_multiplicity"] = -1
    result["han_candidate_id"] = ""
    result["han_candidate_order"] = -1
    result["han_candidate_source_kind"] = ""
    result["han_generation_645_chb_af"] = np.nan
    return result.loc[:, UNIT_COLUMNS]


def validate_unit_inventory(units: pd.DataFrame) -> pd.DataFrame:
    """Fail closed unless the exact conditioned 100+1 per-model bank exists."""

    if not isinstance(units, pd.DataFrame):
        raise TypeError("units must be a pandas DataFrame")
    missing = sorted(set(UNIT_COLUMNS).difference(units.columns))
    if missing:
        raise ValueError("unit inventory lacks columns: " + ", ".join(missing))
    result = units.loc[:, UNIT_COLUMNS].copy()
    if len(result) != EXPECTED_TOTAL_UNITS or result["unit_id"].duplicated().any():
        raise ValueError("unit inventory must contain exactly 202 unique units")
    if set(result["demography"].astype(str)) != set(DEMOGRAPHY_SETTINGS):
        raise ValueError("unit inventory must contain exactly EAS and Han")
    if result[["demography", "seed"]].duplicated().any():
        raise ValueError("seeds must be unique within demography")

    for demography, settings in DEMOGRAPHY_SETTINGS.items():
        group = result[result["demography"].astype(str) == demography]
        classes = group["simulation_class"].astype(str).value_counts().to_dict()
        if len(group) != EXPECTED_UNITS_PER_DEMOGRAPHY or classes != {
            "neutral": EXPECTED_NEUTRAL_PER_DEMOGRAPHY,
            "selected": EXPECTED_SELECTED_PER_DEMOGRAPHY,
        }:
            raise ValueError(
                f"{demography} must contain exactly 100 neutral + 1 selected"
            )
        if set(group["demography_id"].astype(str)) != {settings["demography_id"]}:
            raise ValueError(f"{demography} demography mapping changed")
        if set(pd.to_numeric(group["present_ne"], errors="raise")) != {
            settings["present_ne"]
        }:
            raise ValueError(f"{demography} present Ne changed")
        if set(pd.to_numeric(group["generation_time_years"], errors="raise")) != {
            settings["generation_time_years"]
        }:
            raise ValueError(f"{demography} generation time changed")

        neutral = group[group["simulation_class"].astype(str) == "neutral"]
        selected = group[group["simulation_class"].astype(str) == "selected"].iloc[0]
        if set(neutral["model_id"].astype(str)) != {settings["neutral_model_id"]}:
            raise ValueError(f"{demography} conditioned-neutral model changed")
        if (
            str(selected["unit_id"]) != settings["selected_unit_id"]
            or str(selected["model_id"]) != settings["selected_cell_id"]
            or int(selected["seed"]) != settings["selected_seed"]
        ):
            raise ValueError(f"{demography} selected unit mapping changed")
        neutral_ids = set(neutral["unit_id"].astype(str))
        expected_ids = {
            f"{settings['neutral_model_id']}_r{index:03d}"
            for index in range(EXPECTED_NEUTRAL_PER_DEMOGRAPHY)
        }
        if neutral_ids != expected_ids:
            raise ValueError(f"{demography} conditioned-neutral unit IDs changed")

    coefficients = pd.to_numeric(result["selection_coefficient"], errors="raise")
    classes = result["simulation_class"].astype(str)
    if not np.allclose(coefficients[classes == "neutral"], 0.0, rtol=0, atol=0):
        raise ValueError("neutral units must have s=0")
    if not np.allclose(
        coefficients[classes == "selected"], SELECTION_COEFFICIENT, rtol=0, atol=0
    ):
        raise ValueError("selected units must have s=0.01")

    numeric = (
        "seed",
        "generation_time_years",
        "present_ne",
        "final_population_af",
        "sample_alt_count",
        "sample_af",
        "birth_generation",
        "allele_age_years",
    )
    for column in numeric:
        values = pd.to_numeric(result[column], errors="coerce").to_numpy(dtype=float)
        if not np.all(np.isfinite(values)):
            raise ValueError(f"unit {column} must be finite")
    seeds = pd.to_numeric(result["seed"]).to_numpy(dtype=float)
    births = pd.to_numeric(result["birth_generation"]).to_numpy(dtype=float)
    counts = pd.to_numeric(result["sample_alt_count"]).to_numpy(dtype=float)
    if (
        not np.array_equal(seeds, np.floor(seeds))
        or np.any(seeds < 0)
        or not np.array_equal(births, np.floor(births))
        or np.any(births <= 0)
        or not np.array_equal(counts, np.floor(counts))
        or np.any((counts < 0) | (counts > analysis.PANEL_HAPLOTYPES))
    ):
        raise ValueError("seed/birth/sample-count integer contract changed")
    population_af = pd.to_numeric(result["final_population_af"]).to_numpy(float)
    sample_af = pd.to_numeric(result["sample_af"]).to_numpy(float)
    if np.any((population_af <= 0) | (population_af >= 1)):
        raise ValueError("all units must be strictly segregating in the population")
    if not np.allclose(
        sample_af,
        counts / analysis.PANEL_HAPLOTYPES,
        rtol=0,
        atol=1e-12,
    ):
        raise ValueError("sample AF is inconsistent with sample count")
    expected_ages = births * pd.to_numeric(result["generation_time_years"]).to_numpy(
        float
    )
    if not np.allclose(
        pd.to_numeric(result["allele_age_years"]).to_numpy(float),
        expected_ages,
        rtol=0,
        atol=1e-9,
    ):
        raise ValueError("allele age is inconsistent with birth generation")

    for column in ("ms_sha256", "source_completion_sha256"):
        if (
            not result[column]
            .astype(str)
            .map(lambda value: bool(_SHA256_RE.fullmatch(value)))
            .all()
        ):
            raise ValueError(f"{column} must contain lowercase SHA-256 digests")
    for column in ("ms_path", "source_completion_path", "conditioning"):
        if (
            result[column].isna().any()
            or (result[column].astype(str).str.strip() == "").any()
        ):
            raise ValueError(f"{column} must contain nonempty values")

    eas = result[
        (result["demography"].astype(str) == "EAS")
        & (result["simulation_class"].astype(str) == "neutral")
    ]
    if set(pd.to_numeric(eas["seed"], errors="raise")) != set(
        range(
            conditioned.EAS_SEED_START,
            conditioned.EAS_SEED_START + EXPECTED_NEUTRAL_PER_DEMOGRAPHY,
        )
    ):
        raise ValueError("EAS conditioned-neutral seed stream changed")
    if not np.allclose(
        pd.to_numeric(eas["final_population_af"]),
        conditioned.EAS_SELECTED_PRESENT_AF,
        rtol=0,
        atol=1e-15,
    ):
        raise ValueError("EAS neutral endpoint no longer matches the selected count")
    eas_importance_columns = (
        "eas_selected_log_importance_weight",
        "eas_global_effective_sample_size",
        "eas_global_weight_sum_over_max",
        "eas_global_max_normalized_weight",
    )
    eas_importance = eas.loc[:, eas_importance_columns].apply(
        pd.to_numeric, errors="coerce"
    )
    if not np.isfinite(eas_importance.to_numpy(dtype=float)).all():
        raise ValueError("EAS importance diagnostics must be finite")
    eas_approximation = _strict_boolean(
        eas["eas_time_reversal_is_approximation"],
        "eas_time_reversal_is_approximation",
    )
    non_eas_approximation = _strict_boolean(
        result.loc[
            result.index.difference(eas.index), "eas_time_reversal_is_approximation"
        ],
        "eas_time_reversal_is_approximation",
    )
    if (
        not eas_approximation.all()
        or non_eas_approximation.any()
        or (
            eas["eas_global_effective_sample_size"] < conditioned.EAS_MIN_GLOBAL_ESS
        ).any()
        or (
            eas["eas_global_weight_sum_over_max"]
            < conditioned.EAS_MIN_WEIGHT_SUM_OVER_MAX
        ).any()
        or (
            eas["eas_global_max_normalized_weight"]
            > conditioned.EAS_MAX_NORMALIZED_WEIGHT
        ).any()
    ):
        raise ValueError("EAS global-bank approximation diagnostics fail their gate")
    if (
        eas["eas_source_particle_id"].astype(str).str.strip().eq("").any()
        or set(pd.to_numeric(eas["eas_draw_occurrence"], errors="coerce"))
        != set(range(EXPECTED_NEUTRAL_PER_DEMOGRAPHY))
        or (pd.to_numeric(eas["eas_particle_multiplicity"], errors="coerce") < 1).any()
    ):
        raise ValueError("EAS particle/draw occurrence provenance changed")
    observed_multiplicity = (
        eas["eas_source_particle_id"]
        .astype(str)
        .map(eas["eas_source_particle_id"].astype(str).value_counts())
    )
    if not np.array_equal(
        observed_multiplicity.to_numpy(dtype=int),
        pd.to_numeric(eas["eas_particle_multiplicity"]).to_numpy(dtype=int),
    ):
        raise ValueError("EAS particle multiplicities disagree with occurrence rows")

    han = result[
        (result["demography"].astype(str) == "Han")
        & (result["simulation_class"].astype(str) == "neutral")
    ]
    han_gate = pd.to_numeric(han["han_generation_645_chb_af"], errors="coerce")
    han_order = pd.to_numeric(han["han_candidate_order"], errors="coerce")
    if (
        not np.isfinite(han_gate.to_numpy(float)).all()
        or not han_gate.between(
            conditioned.HAN_GATE_LOWER, conditioned.HAN_GATE_UPPER, inclusive="both"
        ).all()
        or han["han_candidate_id"].astype(str).str.strip().eq("").any()
        or han["han_candidate_id"].duplicated().any()
        or not np.array_equal(han_order, np.floor(han_order))
        or han_order.duplicated().any()
        or set(pd.to_numeric(han["birth_generation"]))
        != {conditioned.HAN_BIRTH_GENERATION}
    ):
        raise ValueError("Han gate/candidate provenance contract changed")
    return result.sort_values(
        ["demography", "simulation_class", "unit_id"], kind="mergesort"
    ).reset_index(drop=True)


LEGACY_V5_PLAN_OUTPUTS = {
    "eas_units": "eas_units.tsv",
    "eas_proposal_blocks": "eas_proposal_blocks.tsv",
    "eas_base_bank_manifest": "eas_base_bank_manifest.json",
    "eas_extension_plan": "eas_extension_plan.json",
    "eas_parallel_replay_plan": "eas_parallel_replay_plan.json",
    "han_candidates": "han_candidates.tsv",
    "han_screen_config": "han_screen.params",
    "han_replay_config": "han_replay.params",
}
LEGACY_V2_PLAN_OUTPUTS = {
    key: value
    for key, value in LEGACY_V5_PLAN_OUTPUTS.items()
    if key != "eas_parallel_replay_plan"
}
LEGACY_V3_PLAN_OUTPUTS = dict(LEGACY_V5_PLAN_OUTPUTS)
LEGACY_V4_PLAN_OUTPUTS = dict(LEGACY_V5_PLAN_OUTPUTS)
LEGACY_V6_PLAN_OUTPUTS = {
    **LEGACY_V5_PLAN_OUTPUTS,
    "recombination_map": "recombination_map.tsv",
}
PLAN_OUTPUTS = dict(LEGACY_V6_PLAN_OUTPUTS)


def _plan_payloads(
    repo_root: Path,
    *,
    eas_base_bank_manifest: Mapping[str, Any],
    eas_extension_plan: Mapping[str, Any],
    eas_parallel_replay_plan: Mapping[str, Any] | None = None,
    natural_run_dir: str | Path | None = None,
) -> dict[str, pd.DataFrame | str | dict[str, Any]]:
    natural_root = (
        Path(natural_run_dir).resolve()
        if natural_run_dir is not None
        else repo_root
        / DEFAULT_NATURAL_WORK_RELATIVE
        / natural_workflow.NEUTRAL_RUN_DIRNAME
    )
    tables = conditioned.build_plan_tables(existing_han_run_dir=natural_root)
    blocks = pd.DataFrame(
        {
            "block_index": np.arange(conditioned.EAS_PROPOSAL_BLOCK_COUNT),
            "block_seed": np.arange(
                conditioned.EAS_PROPOSAL_BLOCK_SEED_START,
                conditioned.EAS_PROPOSAL_BLOCK_SEED_START
                + conditioned.EAS_PROPOSAL_BLOCK_COUNT,
            ),
            "proposals": conditioned.EAS_PROPOSALS_PER_BLOCK,
        }
    )
    blocks["proposal_start"] = (
        blocks["block_index"] * conditioned.EAS_PROPOSALS_PER_BLOCK
    )
    blocks["proposal_stop_exclusive"] = (
        blocks["proposal_start"] + conditioned.EAS_PROPOSALS_PER_BLOCK
    )
    blocks["source_contract"] = np.where(
        blocks["block_index"] < conditioned.EAS_BASE_PROPOSAL_BLOCK_COUNT,
        "frozen_v1_pilot_hash_attested",
        "post_pilot_v2_extension_current_source",
    )
    blocks["proposal_bank_contract_id"] = conditioned.EAS_PROPOSAL_BANK_CONTRACT_ID
    payloads: dict[str, pd.DataFrame | str | dict[str, Any]] = {
        "eas_units": tables["eas_units"],
        "eas_proposal_blocks": blocks,
        "eas_base_bank_manifest": dict(eas_base_bank_manifest),
        "eas_extension_plan": dict(eas_extension_plan),
        "han_candidates": tables["han_candidates"],
        "han_screen_config": conditioned.render_han_parameter_text(
            repo_root, screen=True
        ),
        "han_replay_config": conditioned.render_han_parameter_text(
            repo_root, screen=False
        ),
    }
    if eas_parallel_replay_plan is not None:
        payloads["eas_parallel_replay_plan"] = dict(eas_parallel_replay_plan)
    return payloads


def _plan_source_binding(
    repo_root: Path,
    *,
    workflow_sha256: str | None = None,
) -> dict[str, Any]:
    return {
        "conditioned_workflow_sha256": (
            str(workflow_sha256)
            if workflow_sha256 is not None
            else sha256_file(MODULE_PATH)
        ),
        "conditioned_null_sha256": sha256_file(conditioned.MODULE_PATH),
        "framework_sha256": sha256_file(conditioned.framework.MODULE_PATH),
        "eas_resource": conditioned.framework.EAS_RESOURCE_PATH,
        "eas_resource_sha256": sha256_file(
            repo_root / conditioned.framework.EAS_RESOURCE_PATH
        ),
    }


def _eas_bank_plan_contract(
    base_manifest: Mapping[str, Any],
    extension_plan: Mapping[str, Any],
    *,
    source_binding: Mapping[str, Any],
) -> dict[str, Any]:
    bank_contract = conditioned.EasProposalBankContract().to_record()
    return {
        "schema": EAS_BANK_WORKFLOW_SCHEMA,
        "accepted_per_demography": EXPECTED_NEUTRAL_PER_DEMOGRAPHY,
        "eas_proposal_bank_contract": bank_contract,
        "eas_proposal_bank_contract_sha256": _canonical_sha256(bank_contract),
        "eas_total_proposals": conditioned.EAS_TOTAL_PROPOSALS,
        "eas_global_resample_seed": conditioned.EAS_CATEGORICAL_RESAMPLE_SEED,
        "eas_resampling": "fixed_seed_categorical_multinomial_with_replacement",
        "eas_fixed_bank_rule": (
            "exact_blocks_0_through_199_evaluated_once_no_adaptive_stopping"
        ),
        "eas_base_manifest_payload_sha256": base_manifest["manifest_payload_sha256"],
        "eas_extension_plan_payload_sha256": extension_plan["plan_payload_sha256"],
        "eas_post_pilot_revision": {
            "base_total_proposals": conditioned.EAS_BASE_TOTAL_PROPOSALS,
            "pilot_failed_quality_gates": base_manifest["base_failure_diagnostics"][
                "pilot_failed_quality_gates"
            ],
            "v2_failed_quality_gates_at_pilot_size": base_manifest[
                "base_failure_diagnostics"
            ]["v2_failed_quality_gates"],
            "base_proposal_bank_binding_sha256": base_manifest[
                "base_failure_diagnostics"
            ]["proposal_bank_binding_sha256"],
            "extension_block_count": conditioned.EAS_V2_EXTENSION_PROPOSAL_BLOCK_COUNT,
            "target_total_proposals": conditioned.EAS_V2_TOTAL_PROPOSALS,
        },
        "han_candidate_count": conditioned.HAN_CANDIDATE_COUNT,
        "han_selection_rule": (
            "first_100_candidate_order_joint_g645_gate_and_present_segregation"
        ),
        "source_binding": dict(source_binding),
    }


def _legacy_v2_source_binding(repo_root: Path) -> dict[str, Any]:
    binding = _plan_source_binding(repo_root)
    binding["conditioned_workflow_sha256"] = LEGACY_V2_WORKFLOW_SHA256
    return binding


def _verify_legacy_v2_plan(
    repo_root: Path,
    *,
    legacy_plan_dir: Path,
    proposal_dir: Path,
    natural_run_dir: str | Path | None = None,
    current_replay_workflow_sha256: str | None = None,
) -> dict[str, Any]:
    """Attest the exact frozen A905 v2 plan without rewriting its bytes."""

    bundle = legacy_plan_dir.resolve()
    expected_names = {PLAN_COMPLETION_FILENAME, *LEGACY_V2_PLAN_OUTPUTS.values()}
    children = list(bundle.iterdir()) if bundle.is_dir() else []
    if (
        bundle.is_symlink()
        or {path.name for path in children} != expected_names
        or any(path.is_symlink() or path.is_dir() for path in children)
    ):
        raise ValueError("legacy v2 conditioned plan inventory changed")
    completion_path = bundle / PLAN_COMPLETION_FILENAME
    if sha256_file(completion_path) != LEGACY_V2_PLAN_COMPLETION_SHA256:
        raise ValueError("legacy v2 conditioned plan completion hash changed")
    completion = _read_json(completion_path, "legacy v2 plan completion")
    outputs = completion.get("outputs")
    contract = completion.get("contract")
    if (
        completion.get("schema") != LEGACY_V2_PLAN_SCHEMA
        or completion.get("status") != "complete"
        or not isinstance(outputs, dict)
        or set(outputs) != set(LEGACY_V2_PLAN_OUTPUTS)
        or completion.get("output_count") != len(LEGACY_V2_PLAN_OUTPUTS)
        or not isinstance(contract, dict)
        or completion.get("contract_sha256") != _canonical_sha256(contract)
    ):
        raise ValueError("legacy v2 conditioned plan completion changed")
    for label, filename in LEGACY_V2_PLAN_OUTPUTS.items():
        path = bundle / filename
        record = outputs[label]
        if (
            record.get("path") != filename
            or path.stat().st_size != int(record.get("size_bytes", -1))
            or sha256_file(path) != record.get("sha256")
        ):
            raise ValueError(f"legacy v2 conditioned plan checksum failed: {label}")
    base_manifest = _read_json(
        bundle / LEGACY_V2_PLAN_OUTPUTS["eas_base_bank_manifest"],
        "legacy v2 EAS base-bank manifest",
    )
    extension_plan = _read_json(
        bundle / LEGACY_V2_PLAN_OUTPUTS["eas_extension_plan"],
        "legacy v2 EAS extension plan",
    )
    conditioned.validate_eas_base_bank_manifest(
        repo_root,
        base_manifest,
        proposal_dir=proposal_dir,
    )
    conditioned.validate_eas_extension_plan(
        repo_root,
        extension_plan,
        base_manifest=base_manifest,
    )
    expected_contract = _eas_bank_plan_contract(
        base_manifest,
        extension_plan,
        source_binding=_legacy_v2_source_binding(repo_root),
    )
    if contract != expected_contract:
        raise ValueError("legacy v2 conditioned plan contract changed")
    expected_payloads = _plan_payloads(
        repo_root,
        eas_base_bank_manifest=base_manifest,
        eas_extension_plan=extension_plan,
        natural_run_dir=natural_run_dir,
    )
    for label, filename in LEGACY_V2_PLAN_OUTPUTS.items():
        expected = expected_payloads[label]
        path = bundle / filename
        if isinstance(expected, pd.DataFrame):
            observed = pd.read_csv(path, sep="\t")
            pd.testing.assert_frame_equal(
                observed,
                expected,
                check_dtype=False,
                check_exact=False,
                rtol=0,
                atol=1e-15,
            )
            if int(outputs[label].get("rows", -1)) != len(expected):
                raise ValueError(f"legacy v2 plan row count changed: {label}")
        elif isinstance(expected, dict):
            if _read_json(path, f"legacy v2 plan {label}") != expected:
                raise ValueError(f"legacy v2 plan JSON changed: {label}")
        elif path.read_text(encoding="utf-8") != expected:
            raise ValueError(f"legacy v2 plan config changed: {label}")
    current_replay_sha256 = (
        str(current_replay_workflow_sha256)
        if current_replay_workflow_sha256 is not None
        else sha256_file(MODULE_PATH)
    )
    attestation: dict[str, Any] = {
        "schema": "gamma-smc.cosi2-conditioned-gamma-v2-plan-attestation/v3",
        "status": "exact_frozen_predecessor",
        "legacy_plan_dirname": LEGACY_V2_PLAN_DIRNAME,
        "legacy_plan_completion_sha256": sha256_file(completion_path),
        "legacy_plan_contract_sha256": completion["contract_sha256"],
        "legacy_plan_producer_workflow_sha256": LEGACY_V2_WORKFLOW_SHA256,
        "current_replay_workflow_sha256": current_replay_sha256,
        "base_manifest_payload_sha256": base_manifest["manifest_payload_sha256"],
        "extension_plan_payload_sha256": extension_plan["plan_payload_sha256"],
        "legacy_output_sha256": {
            label: str(record["sha256"]) for label, record in outputs.items()
        },
        "legacy_output_bindings": {
            label: {
                "size_bytes": int(record["size_bytes"]),
                "sha256": str(record["sha256"]),
            }
            for label, record in outputs.items()
        },
    }
    attestation["attestation_payload_sha256"] = _canonical_sha256(attestation)
    return {
        "attestation": attestation,
        "base_manifest": base_manifest,
        "extension_plan": extension_plan,
    }


def _verify_legacy_v3_replay_root(replay_dir: Path) -> dict[str, Any]:
    """Attest the exact 13 unpublished v3 replay bundles without changing them."""

    root = replay_dir.resolve()
    expected_names = {f"block_{index:02d}" for index in LEGACY_V3_REPLAY_BLOCK_INDEXES}
    children = list(root.iterdir()) if root.is_dir() else []
    if (
        root.is_symlink()
        or {path.name for path in children} != expected_names
        or any(path.is_symlink() or not path.is_dir() for path in children)
    ):
        raise ValueError("legacy v3 replay root inventory changed")
    entries: list[dict[str, Any]] = []
    for block_index in LEGACY_V3_REPLAY_BLOCK_INDEXES:
        bundle = root / f"block_{block_index:02d}"
        files = list(bundle.iterdir())
        if {path.name for path in files} != {
            "completion.json",
            "selected_rows.tsv",
            "paths.npz",
        } or any(path.is_symlink() or not path.is_file() for path in files):
            raise ValueError(f"legacy v3 replay block inventory changed: {block_index}")
        entries.append(
            {
                "block_dir": bundle.name,
                "files": {
                    path.name: {
                        "size_bytes": path.stat().st_size,
                        "sha256": sha256_file(path),
                    }
                    for path in sorted(files, key=lambda value: value.name)
                },
            }
        )
    manifest_sha256 = _canonical_sha256(entries)
    if manifest_sha256 != LEGACY_V3_REPLAY_MANIFEST_SHA256:
        raise ValueError("legacy v3 replay bundle hashes changed")
    return {
        "schema": "gamma-smc.cosi2-eas-legacy-v3-replay-manifest/v4",
        "status": "exact_unpublished_predecessor",
        "legacy_replay_dirname": LEGACY_V3_PARALLEL_REPLAY_DIRNAME,
        "legacy_block_indexes": list(LEGACY_V3_REPLAY_BLOCK_INDEXES),
        "legacy_block_count": len(entries),
        "legacy_block_bindings": entries,
        "legacy_replay_manifest_sha256": manifest_sha256,
        "legacy_root_completion_present": False,
        "legacy_selection_published": False,
    }


def _verify_legacy_v3_plan(
    repo_root: Path,
    *,
    legacy_plan_dir: Path,
    legacy_replay_dir: Path,
    proposal_dir: Path,
    natural_run_dir: str | Path | None = None,
    current_replay_workflow_sha256: str | None = None,
) -> dict[str, Any]:
    """Attest the exact 827e v3 plan and its exact unpublished replay cache."""

    bundle = legacy_plan_dir.resolve()
    expected_names = {PLAN_COMPLETION_FILENAME, *LEGACY_V3_PLAN_OUTPUTS.values()}
    children = list(bundle.iterdir()) if bundle.is_dir() else []
    if (
        bundle.is_symlink()
        or {path.name for path in children} != expected_names
        or any(path.is_symlink() or not path.is_file() for path in children)
    ):
        raise ValueError("legacy v3 conditioned plan inventory changed")
    completion_path = bundle / PLAN_COMPLETION_FILENAME
    if sha256_file(completion_path) != LEGACY_V3_PLAN_COMPLETION_SHA256:
        raise ValueError("legacy v3 conditioned plan completion hash changed")
    completion = _read_json(completion_path, "legacy v3 plan completion")
    outputs = completion.get("outputs")
    contract = completion.get("contract")
    if (
        completion.get("schema") != LEGACY_V3_PLAN_SCHEMA
        or completion.get("status") != "complete"
        or not isinstance(outputs, dict)
        or set(outputs) != set(LEGACY_V3_PLAN_OUTPUTS)
        or completion.get("output_count") != len(LEGACY_V3_PLAN_OUTPUTS)
        or not isinstance(contract, dict)
        or completion.get("contract_sha256") != _canonical_sha256(contract)
        or contract.get("source_binding")
        != _plan_source_binding(
            repo_root,
            workflow_sha256=LEGACY_V3_WORKFLOW_SHA256,
        )
    ):
        raise ValueError("legacy v3 conditioned plan completion changed")
    for label, filename in LEGACY_V3_PLAN_OUTPUTS.items():
        path = bundle / filename
        record = outputs[label]
        if (
            record.get("path") != filename
            or path.stat().st_size != int(record.get("size_bytes", -1))
            or sha256_file(path) != record.get("sha256")
        ):
            raise ValueError(f"legacy v3 conditioned plan checksum failed: {label}")
    v2_migration = _verify_legacy_v2_plan(
        repo_root,
        legacy_plan_dir=bundle.parent / LEGACY_V2_PLAN_DIRNAME,
        proposal_dir=proposal_dir,
        natural_run_dir=natural_run_dir,
        current_replay_workflow_sha256=LEGACY_V3_WORKFLOW_SHA256,
    )
    base_manifest = _read_json(
        bundle / LEGACY_V3_PLAN_OUTPUTS["eas_base_bank_manifest"],
        "legacy v3 EAS base-bank manifest",
    )
    extension_plan = _read_json(
        bundle / LEGACY_V3_PLAN_OUTPUTS["eas_extension_plan"],
        "legacy v3 EAS extension plan",
    )
    if (
        base_manifest != v2_migration["base_manifest"]
        or extension_plan != v2_migration["extension_plan"]
    ):
        raise ValueError("legacy v3 inherited bank payload changed")
    replay_plan_path = bundle / LEGACY_V3_PLAN_OUTPUTS["eas_parallel_replay_plan"]
    if sha256_file(replay_plan_path) != LEGACY_V3_REPLAY_PLAN_FILE_SHA256:
        raise ValueError("legacy v3 replay-plan file hash changed")
    replay_plan = _read_json(replay_plan_path, "legacy v3 parallel replay plan")
    replay_payload = {
        key: value for key, value in replay_plan.items() if key != "plan_payload_sha256"
    }
    if (
        replay_plan.get("schema") != LEGACY_V3_PARALLEL_REPLAY_PLAN_SCHEMA
        or replay_plan.get("status") != "complete"
        or replay_plan.get("plan_payload_sha256")
        != LEGACY_V3_REPLAY_PLAN_PAYLOAD_SHA256
        or replay_plan.get("plan_payload_sha256") != _canonical_sha256(replay_payload)
        or replay_plan.get("legacy_v2_plan_attestation") != v2_migration["attestation"]
        or replay_plan.get("source_binding", {}).get("workflow_sha256")
        != LEGACY_V3_WORKFLOW_SHA256
        or replay_plan.get("source_binding", {}).get("core_sha256")
        != sha256_file(conditioned.MODULE_PATH)
        or contract.get("legacy_v2_plan_attestation") != v2_migration["attestation"]
        or contract.get("eas_parallel_replay_plan_payload_sha256")
        != LEGACY_V3_REPLAY_PLAN_PAYLOAD_SHA256
    ):
        raise ValueError("legacy v3 replay-plan provenance changed")
    for label in LEGACY_V2_PLAN_OUTPUTS:
        if {
            "size_bytes": int(outputs[label]["size_bytes"]),
            "sha256": str(outputs[label]["sha256"]),
        } != v2_migration["attestation"]["legacy_output_bindings"][label]:
            raise ValueError(f"legacy v3 inherited v2 bytes changed: {label}")
    replay_manifest = _verify_legacy_v3_replay_root(legacy_replay_dir)
    attestation: dict[str, Any] = {
        "schema": "gamma-smc.cosi2-conditioned-gamma-v3-plan-attestation/v4",
        "status": "exact_frozen_predecessor_with_unpublished_replay_cache",
        "legacy_plan_dirname": LEGACY_V3_PLAN_DIRNAME,
        "legacy_plan_completion_sha256": sha256_file(completion_path),
        "legacy_plan_contract_sha256": str(completion["contract_sha256"]),
        "legacy_plan_producer_workflow_sha256": LEGACY_V3_WORKFLOW_SHA256,
        "current_replay_workflow_sha256": (
            str(current_replay_workflow_sha256)
            if current_replay_workflow_sha256 is not None
            else sha256_file(MODULE_PATH)
        ),
        "base_manifest_payload_sha256": base_manifest["manifest_payload_sha256"],
        "extension_plan_payload_sha256": extension_plan["plan_payload_sha256"],
        "legacy_replay_plan_file_sha256": sha256_file(replay_plan_path),
        "legacy_replay_plan_payload_sha256": replay_plan["plan_payload_sha256"],
        "legacy_output_bindings": {
            label: {
                "size_bytes": int(record["size_bytes"]),
                "sha256": str(record["sha256"]),
            }
            for label, record in outputs.items()
        },
        "legacy_v2_plan_attestation": v2_migration["attestation"],
        "legacy_replay_manifest": replay_manifest,
    }
    attestation["attestation_payload_sha256"] = _canonical_sha256(attestation)
    return {
        "attestation": attestation,
        "base_manifest": base_manifest,
        "extension_plan": extension_plan,
        "legacy_replay_plan": replay_plan,
        "legacy_replay_manifest": replay_manifest,
    }


def _legacy_v4_replay_root_manifest(replay_dir: Path) -> dict[str, Any]:
    """Hash-attest the exact completed v4 replay root without reinterpreting it."""

    root = replay_dir.resolve()
    completion_path = root / EAS_PARALLEL_REPLAY_COMPLETION_FILENAME
    if root.is_symlink() or not root.is_dir():
        raise ValueError("legacy v4 replay root is absent")
    if sha256_file(completion_path) != LEGACY_V4_REPLAY_COMPLETION_SHA256:
        raise ValueError("legacy v4 replay completion hash changed")
    completion = _read_json(completion_path, "legacy v4 replay completion")
    contract = completion.get("contract")
    bindings = (
        contract.get("selected_block_bindings") if isinstance(contract, dict) else None
    )
    if (
        completion.get("schema") != LEGACY_V4_PARALLEL_REPLAY_COMPLETION_SCHEMA
        or completion.get("status") != "complete"
        or not isinstance(contract, dict)
        or completion.get("contract_sha256") != _canonical_sha256(contract)
        or contract.get("workflow_sha256") != LEGACY_V4_WORKFLOW_SHA256
        or contract.get("plan_completion_sha256") != LEGACY_V4_PLAN_COMPLETION_SHA256
        or contract.get("parallel_replay_plan_payload_sha256")
        != LEGACY_V4_REPLAY_PLAN_PAYLOAD_SHA256
        or contract.get("selected_count") != EXPECTED_NEUTRAL_PER_DEMOGRAPHY
        or contract.get("unique_selected_block_count") != 81
        or not isinstance(bindings, list)
        or len(bindings) != 81
    ):
        raise ValueError("legacy v4 replay completion contract changed")
    block_indexes = [int(binding["block_index"]) for binding in bindings]
    if block_indexes != sorted(set(block_indexes)):
        raise ValueError("legacy v4 replay block ledger changed")
    expected_names = {
        EAS_PARALLEL_REPLAY_COMPLETION_FILENAME,
        *(f"block_{index:02d}" for index in block_indexes),
    }
    children = list(root.iterdir())
    if {path.name for path in children} != expected_names or any(
        path.is_symlink()
        or (path.name == EAS_PARALLEL_REPLAY_COMPLETION_FILENAME and not path.is_file())
        or (path.name != EAS_PARALLEL_REPLAY_COMPLETION_FILENAME and not path.is_dir())
        for path in children
    ):
        raise ValueError("legacy v4 replay root inventory changed")
    for block_index in block_indexes:
        bundle = root / f"block_{block_index:02d}"
        files = list(bundle.iterdir())
        if {path.name for path in files} != {
            "completion.json",
            "selected_rows.tsv",
            "paths.npz",
        } or any(path.is_symlink() or not path.is_file() for path in files):
            raise ValueError(f"legacy v4 replay block inventory changed: {block_index}")
    file_bindings = [
        {
            "path": path.relative_to(root).as_posix(),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        for path in sorted(
            root.rglob("*"), key=lambda value: value.relative_to(root).as_posix()
        )
        if path.is_file()
    ]
    manifest_sha256 = _canonical_sha256(file_bindings)
    if (
        len(file_bindings) != 244
        or manifest_sha256 != LEGACY_V4_REPLAY_ROOT_MANIFEST_SHA256
    ):
        raise ValueError("legacy v4 replay root file hashes changed")
    return {
        "schema": "gamma-smc.cosi2-eas-legacy-v4-replay-root-manifest/v5",
        "status": "exact_complete_unpublished_predecessor",
        "legacy_replay_dirname": LEGACY_V4_PARALLEL_REPLAY_DIRNAME,
        "legacy_replay_completion_sha256": LEGACY_V4_REPLAY_COMPLETION_SHA256,
        "legacy_block_indexes": block_indexes,
        "legacy_block_count": len(block_indexes),
        "legacy_file_count": len(file_bindings),
        "legacy_file_bindings": file_bindings,
        "legacy_replay_root_manifest_sha256": manifest_sha256,
    }


def _verify_legacy_v4_unpublished_state(work_root: Path) -> dict[str, Any]:
    """Require old selection/staging absence while permitting v5 publication."""

    root = work_root.resolve()
    forbidden = [
        path
        for path in root.iterdir()
        if path.name == LEGACY_V4_SELECTION_DIRNAME
        or path.name.startswith(f".{LEGACY_V4_SELECTION_DIRNAME}.staging.")
        or path.name.startswith(f".{LEGACY_V4_PARALLEL_REPLAY_DIRNAME}.staging.")
    ]
    if forbidden:
        raise ValueError("legacy v4 selection or replay staging was published")
    return {
        "legacy_selection_absent": True,
        "legacy_selection_staging_absent": True,
        "legacy_replay_staging_absent": True,
        "current_v5_selection_allowed": True,
    }


def _verify_legacy_v4_plan(
    repo_root: Path,
    *,
    legacy_plan_dir: Path,
    legacy_replay_dir: Path,
    proposal_dir: Path,
    natural_run_dir: str | Path | None = None,
    current_assembly_workflow_sha256: str | None = None,
) -> dict[str, Any]:
    """Attest the exact 5e08 v4 plan and completed unpublished replay root."""

    bundle = legacy_plan_dir.resolve()
    work_root = bundle.parent
    expected_names = {PLAN_COMPLETION_FILENAME, *LEGACY_V4_PLAN_OUTPUTS.values()}
    children = list(bundle.iterdir()) if bundle.is_dir() else []
    if (
        bundle.is_symlink()
        or {path.name for path in children} != expected_names
        or any(path.is_symlink() or not path.is_file() for path in children)
    ):
        raise ValueError("legacy v4 conditioned plan inventory changed")
    unpublished_state = _verify_legacy_v4_unpublished_state(work_root)
    completion_path = bundle / PLAN_COMPLETION_FILENAME
    if sha256_file(completion_path) != LEGACY_V4_PLAN_COMPLETION_SHA256:
        raise ValueError("legacy v4 conditioned plan completion hash changed")
    completion = _read_json(completion_path, "legacy v4 plan completion")
    outputs = completion.get("outputs")
    contract = completion.get("contract")
    if (
        completion.get("schema") != LEGACY_V4_PLAN_SCHEMA
        or completion.get("status") != "complete"
        or not isinstance(outputs, dict)
        or set(outputs) != set(LEGACY_V4_PLAN_OUTPUTS)
        or completion.get("output_count") != len(LEGACY_V4_PLAN_OUTPUTS)
        or not isinstance(contract, dict)
        or completion.get("contract_sha256") != _canonical_sha256(contract)
        or contract.get("source_binding")
        != _plan_source_binding(
            repo_root,
            workflow_sha256=LEGACY_V4_WORKFLOW_SHA256,
        )
    ):
        raise ValueError("legacy v4 conditioned plan completion changed")
    for label, filename in LEGACY_V4_PLAN_OUTPUTS.items():
        path = bundle / filename
        record = outputs[label]
        if (
            record.get("path") != filename
            or path.stat().st_size != int(record.get("size_bytes", -1))
            or sha256_file(path) != record.get("sha256")
        ):
            raise ValueError(f"legacy v4 conditioned plan checksum failed: {label}")
    v3_migration = _verify_legacy_v3_plan(
        repo_root,
        legacy_plan_dir=bundle.parent / LEGACY_V3_PLAN_DIRNAME,
        legacy_replay_dir=bundle.parent / LEGACY_V3_PARALLEL_REPLAY_DIRNAME,
        proposal_dir=proposal_dir,
        natural_run_dir=natural_run_dir,
        current_replay_workflow_sha256=LEGACY_V4_WORKFLOW_SHA256,
    )
    base_manifest = _read_json(
        bundle / LEGACY_V4_PLAN_OUTPUTS["eas_base_bank_manifest"],
        "legacy v4 EAS base-bank manifest",
    )
    extension_plan = _read_json(
        bundle / LEGACY_V4_PLAN_OUTPUTS["eas_extension_plan"],
        "legacy v4 EAS extension plan",
    )
    if (
        base_manifest != v3_migration["base_manifest"]
        or extension_plan != v3_migration["extension_plan"]
    ):
        raise ValueError("legacy v4 inherited bank payload changed")
    replay_plan_path = bundle / LEGACY_V4_PLAN_OUTPUTS["eas_parallel_replay_plan"]
    if sha256_file(replay_plan_path) != LEGACY_V4_REPLAY_PLAN_FILE_SHA256:
        raise ValueError("legacy v4 replay-plan file hash changed")
    replay_plan = _read_json(replay_plan_path, "legacy v4 parallel replay plan")
    replay_payload = {
        key: value for key, value in replay_plan.items() if key != "plan_payload_sha256"
    }
    if (
        replay_plan.get("schema") != LEGACY_V4_PARALLEL_REPLAY_PLAN_SCHEMA
        or replay_plan.get("status") != "complete"
        or replay_plan.get("plan_payload_sha256")
        != LEGACY_V4_REPLAY_PLAN_PAYLOAD_SHA256
        or replay_plan.get("plan_payload_sha256") != _canonical_sha256(replay_payload)
        or replay_plan.get("legacy_v3_plan_attestation") != v3_migration["attestation"]
        or replay_plan.get("source_binding", {}).get("workflow_sha256")
        != LEGACY_V4_WORKFLOW_SHA256
        or replay_plan.get("source_binding", {}).get("core_sha256")
        != sha256_file(conditioned.MODULE_PATH)
        or contract.get("legacy_v3_plan_attestation") != v3_migration["attestation"]
        or contract.get("eas_parallel_replay_plan_payload_sha256")
        != LEGACY_V4_REPLAY_PLAN_PAYLOAD_SHA256
    ):
        raise ValueError("legacy v4 replay-plan provenance changed")
    for label in LEGACY_V2_PLAN_OUTPUTS:
        if {
            "size_bytes": int(outputs[label]["size_bytes"]),
            "sha256": str(outputs[label]["sha256"]),
        } != v3_migration["attestation"]["legacy_output_bindings"][label]:
            raise ValueError(f"legacy v4 inherited v3 bytes changed: {label}")
    replay_manifest = _legacy_v4_replay_root_manifest(legacy_replay_dir)
    attestation: dict[str, Any] = {
        "schema": "gamma-smc.cosi2-conditioned-gamma-v4-plan-attestation/v5",
        "status": "exact_frozen_predecessor_with_complete_unpublished_replay",
        "legacy_plan_dirname": LEGACY_V4_PLAN_DIRNAME,
        "legacy_plan_completion_sha256": LEGACY_V4_PLAN_COMPLETION_SHA256,
        "legacy_plan_contract_sha256": str(completion["contract_sha256"]),
        "legacy_plan_producer_workflow_sha256": LEGACY_V4_WORKFLOW_SHA256,
        "current_assembly_workflow_sha256": (
            str(current_assembly_workflow_sha256)
            if current_assembly_workflow_sha256 is not None
            else sha256_file(MODULE_PATH)
        ),
        "base_manifest_payload_sha256": base_manifest["manifest_payload_sha256"],
        "extension_plan_payload_sha256": extension_plan["plan_payload_sha256"],
        "legacy_replay_plan_file_sha256": LEGACY_V4_REPLAY_PLAN_FILE_SHA256,
        "legacy_replay_plan_payload_sha256": LEGACY_V4_REPLAY_PLAN_PAYLOAD_SHA256,
        "legacy_output_bindings": {
            label: {
                "size_bytes": int(record["size_bytes"]),
                "sha256": str(record["sha256"]),
            }
            for label, record in outputs.items()
        },
        "legacy_v3_plan_attestation": v3_migration["attestation"],
        "legacy_replay_root_manifest": replay_manifest,
        "legacy_unpublished_state": unpublished_state,
    }
    attestation["attestation_payload_sha256"] = _canonical_sha256(attestation)
    return {
        "attestation": attestation,
        "base_manifest": base_manifest,
        "extension_plan": extension_plan,
        "legacy_replay_plan": replay_plan,
        "legacy_replay_manifest": replay_manifest,
        "legacy_replay_completion": _read_json(
            legacy_replay_dir / EAS_PARALLEL_REPLAY_COMPLETION_FILENAME,
            "legacy v4 replay completion",
        ),
    }


def _build_eas_parallel_replay_plan(
    legacy_v4_attestation: Mapping[str, Any],
    *,
    workflow_sha256: str | None = None,
) -> dict[str, Any]:
    producer_sha256 = (
        str(workflow_sha256)
        if workflow_sha256 is not None
        else sha256_file(MODULE_PATH)
    )
    plan: dict[str, Any] = {
        "schema": EAS_PARALLEL_REPLAY_PLAN_SCHEMA,
        "status": "complete",
        "execution_method": EAS_PARALLEL_REPLAY_METHOD,
        "migration_contract": {
            "source": "exact_completed_eas_parallel_replays_v4",
            "source_file_count": 244,
            "source_selected_block_count": 81,
            "copy_scope": "entire_root_including_legacy_v4_root_completion",
            "promotion": "adjacent_atomic_directory_replace_after_full_verification",
            "destination_bytes": "byte_identical_to_attested_v4_source_root",
            "proposal_replay_or_resampling": False,
            "partial_nonempty_destination_allowed": False,
            "source_selection_and_staging_required_absent": True,
            "legacy_completion_interpretation": (
                "verify_only_under_hard_bound_v4_producer_plan_and_root_contract"
            ),
        },
        "selection_contract": {
            "method": "fixed_seed_categorical_multinomial_with_replacement",
            "seed": conditioned.EAS_CATEGORICAL_RESAMPLE_SEED,
            "selected_count": EXPECTED_NEUTRAL_PER_DEMOGRAPHY,
            "selection_recomputed_from_exact_fixed_200_block_bank": True,
            "selected_row_adapter": (
                "unit_id_preserved_with_set_index_drop_false_then_real_core_writer"
            ),
        },
        "revision_scope": {
            "proposal_target_inherited_from_pre_extension_v2_plan": True,
            "proposal_generation_complete_before_v5": True,
            "v4_parallel_replay_complete_before_v5": True,
            "v5_changes_selection_row_adapter_and_migration_provenance_only": True,
            "legacy_v2_v3_v4_plan_or_replay_bytes_rewritten": False,
        },
        "source_binding": {
            "workflow_path": _relative(MODULE_PATH, MODULE_PATH.parents[2]),
            "workflow_sha256": producer_sha256,
            "core_path": _relative(
                Path(conditioned.MODULE_PATH), Path(conditioned.MODULE_PATH).parents[2]
            ),
            "core_sha256": sha256_file(conditioned.MODULE_PATH),
            "migration_helper": "_migrate_exact_legacy_v4_replay_root",
            "runner_helper": "run_eas_parallel_selected_replays",
            "selection_row_helper": "_selected_particle_row",
        },
        "legacy_v4_plan_attestation": dict(legacy_v4_attestation),
        "proposal_bytes_mutable": False,
        "v2_v3_or_v4_plan_and_replay_bytes_mutable": False,
        "reuse_exact_attested_complete_unpublished_v4_replay_root": True,
    }
    plan["plan_payload_sha256"] = _canonical_sha256(plan)
    return plan


def _validate_eas_parallel_replay_plan(
    plan: Mapping[str, Any],
    *,
    legacy_v4_attestation: Mapping[str, Any],
    workflow_sha256: str | None = None,
) -> dict[str, Any]:
    expected = _build_eas_parallel_replay_plan(
        legacy_v4_attestation,
        workflow_sha256=workflow_sha256,
    )
    if dict(plan) != expected:
        raise ValueError("EAS replay migration plan differs from its v5 contract")
    return {
        "status": "valid",
        "plan_payload_sha256": str(plan["plan_payload_sha256"]),
        "execution_method": EAS_PARALLEL_REPLAY_METHOD,
        "source_selected_block_count": 81,
    }


def _file_tree_manifest(directory: Path) -> dict[str, Any]:
    """Return a canonical exact-file manifest while rejecting every symlink."""

    root = directory.resolve()
    if directory.is_symlink() or not root.is_dir():
        raise ValueError(f"artifact root is absent: {directory}")
    descendants = sorted(
        root.rglob("*"), key=lambda value: value.relative_to(root).as_posix()
    )
    if any(path.is_symlink() for path in descendants):
        raise ValueError(f"artifact root contains a symlink: {directory}")
    bindings = [
        {
            "path": path.relative_to(root).as_posix(),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        for path in descendants
        if path.is_file()
    ]
    return {
        "file_count": len(bindings),
        "file_bindings": bindings,
        "manifest_sha256": _canonical_sha256(bindings),
    }


def _legacy_v5_plan_contract(
    base_manifest: Mapping[str, Any],
    extension_plan: Mapping[str, Any],
    *,
    legacy_v4_attestation: Mapping[str, Any],
    replay_plan: Mapping[str, Any],
    repo_root: Path,
) -> dict[str, Any]:
    contract = _eas_bank_plan_contract(
        base_manifest,
        extension_plan,
        source_binding=_plan_source_binding(
            repo_root,
            workflow_sha256=LEGACY_V5_WORKFLOW_SHA256,
        ),
    )
    contract.update(
        {
            "eas_parallel_replay_plan_payload_sha256": replay_plan[
                "plan_payload_sha256"
            ],
            "legacy_v4_plan_attestation": dict(legacy_v4_attestation),
            "replay_only_post_extension_revision": {
                "proposal_target_inherited_from_pre_extension_v2": True,
                "proposal_blocks_rewritten": False,
                "legacy_v2_plan_rewritten": False,
                "legacy_v3_or_v4_plan_or_replay_rewritten": False,
                "selected_tsv_float_parser": "pandas_float_precision_round_trip",
                "canonical_float_round_trip_is_bit_exact": True,
                "exact_v4_replay_root_migrated_without_replay": True,
                "selection_unit_id_adapter_fixed": True,
            },
        }
    )
    return contract


def _verify_legacy_v5_eas_state(
    repo_root: Path,
    *,
    replay_dir: Path,
    selection_dir: Path,
) -> dict[str, Any]:
    """Hash-attest the immutable completed v5 replay and selection roots."""

    work_root = selection_dir.resolve().parent
    stale = [
        path.name
        for path in work_root.iterdir()
        if path.name.startswith(f".{LEGACY_V5_PARALLEL_REPLAY_DIRNAME}.staging.")
        or path.name.startswith(f".{LEGACY_V5_SELECTION_DIRNAME}.staging.")
    ]
    if stale:
        raise ValueError("legacy v5 EAS state has stale staging: " + ", ".join(stale))
    replay_manifest = _legacy_v4_replay_root_manifest(replay_dir)
    if (
        replay_manifest["legacy_file_count"] != 244
        or replay_manifest["legacy_replay_root_manifest_sha256"]
        != LEGACY_V4_REPLAY_ROOT_MANIFEST_SHA256
    ):
        raise ValueError("legacy v5 replay root bytes changed")

    bundle = selection_dir.resolve()
    top = list(bundle.iterdir()) if bundle.is_dir() else []
    if (
        bundle.is_symlink()
        or {path.name for path in top}
        != {
            EAS_SELECTION_COMPLETION,
            "selected_occurrences.tsv",
            "global_diagnostics.json",
            "units",
        }
        or any(path.is_symlink() for path in top)
    ):
        raise ValueError("legacy v5 EAS selection top-level inventory changed")
    tree = _file_tree_manifest(bundle)
    if (
        tree["file_count"] != LEGACY_V5_SELECTION_ROOT_FILE_COUNT
        or tree["manifest_sha256"] != LEGACY_V5_SELECTION_ROOT_MANIFEST_SHA256
    ):
        raise ValueError("legacy v5 EAS selection file hashes changed")
    completion_path = bundle / EAS_SELECTION_COMPLETION
    if sha256_file(completion_path) != LEGACY_V5_SELECTION_COMPLETION_SHA256:
        raise ValueError("legacy v5 EAS selection completion hash changed")
    completion = _read_json(completion_path, "legacy v5 EAS selection")
    contract = completion.get("contract")
    if (
        completion.get("schema") != EAS_SELECTION_SCHEMA
        or completion.get("status") != "complete"
        or not isinstance(contract, dict)
        or completion.get("contract_sha256") != _canonical_sha256(contract)
        or completion.get("contract_sha256") != LEGACY_V5_SELECTION_CONTRACT_SHA256
        or contract.get("workflow_sha256") != LEGACY_V5_WORKFLOW_SHA256
        or contract.get("plan_completion_sha256") != LEGACY_V5_PLAN_COMPLETION_SHA256
        or contract.get("parallel_replay_plan_payload_sha256")
        != LEGACY_V5_REPLAY_PLAN_PAYLOAD_SHA256
        or contract.get("parallel_replay_completion_sha256")
        != LEGACY_V4_REPLAY_COMPLETION_SHA256
        or contract.get("conditioned_null_sha256")
        != sha256_file(conditioned.MODULE_PATH)
        or contract.get("selected_count") != EXPECTED_NEUTRAL_PER_DEMOGRAPHY
    ):
        raise ValueError("legacy v5 EAS selection contract changed")

    selected = _read_eas_selected_frame(bundle / "selected_occurrences.tsv")
    expected_ids = [
        f"{conditioned.EAS_MODEL_ID}_r{index:03d}"
        for index in range(EXPECTED_NEUTRAL_PER_DEMOGRAPHY)
    ]
    if (
        len(selected) != EXPECTED_NEUTRAL_PER_DEMOGRAPHY
        or selected["unit_id"].astype(str).tolist() != expected_ids
        or pd.to_numeric(selected["accepted_index"], errors="raise")
        .astype(int)
        .tolist()
        != list(range(EXPECTED_NEUTRAL_PER_DEMOGRAPHY))
        or pd.to_numeric(selected["cosi2_seed"], errors="raise").astype(int).tolist()
        != [
            conditioned.EAS_SEED_START + index
            for index in range(EXPECTED_NEUTRAL_PER_DEMOGRAPHY)
        ]
    ):
        raise ValueError("legacy v5 EAS selected occurrence ledger changed")
    for record in selected.to_dict(orient="records"):
        conditioned._eas_selected_particle_record(_portable_record(record))
    diagnostics = _read_json(
        bundle / "global_diagnostics.json", "legacy v5 EAS diagnostics"
    )
    if (
        diagnostics.get("selected_count") != EXPECTED_NEUTRAL_PER_DEMOGRAPHY
        or float(diagnostics.get("global_effective_sample_size", -1))
        < conditioned.EAS_MIN_GLOBAL_ESS
        or float(diagnostics.get("global_weight_sum_over_max", -1))
        < conditioned.EAS_MIN_WEIGHT_SUM_OVER_MAX
        or float(diagnostics.get("global_max_normalized_weight", math.inf))
        > conditioned.EAS_MAX_NORMALIZED_WEIGHT
        or int(diagnostics.get("cap_hit_count", -1)) != 0
    ):
        raise ValueError("legacy v5 EAS selection diagnostics fail their fixed gate")
    outputs = completion.get("outputs", {})
    unit_outputs = completion.get("unit_outputs", {})
    if set(outputs) != {"selected_occurrences", "global_diagnostics"} or set(
        unit_outputs
    ) != set(expected_ids):
        raise ValueError("legacy v5 EAS selection output ledger changed")
    for label, filename in {
        "selected_occurrences": "selected_occurrences.tsv",
        "global_diagnostics": "global_diagnostics.json",
    }.items():
        path = bundle / filename
        record = outputs[label]
        if (
            record.get("path") != filename
            or path.stat().st_size != int(record.get("size_bytes", -1))
            or sha256_file(path) != record.get("sha256")
        ):
            raise ValueError(f"legacy v5 EAS selection output changed: {label}")
    unit_root = bundle / "units"
    unit_dirs = list(unit_root.iterdir()) if unit_root.is_dir() else []
    if {path.name for path in unit_dirs} != set(expected_ids) or any(
        path.is_symlink() or not path.is_dir() for path in unit_dirs
    ):
        raise ValueError("legacy v5 EAS selection unit inventory changed")
    expected_unit_files = {
        "bridge_counts.tsv",
        "trajectory.tsv",
        "parameters.par",
        "selection.json",
        "completion.json",
    }
    for unit_id in expected_ids:
        unit_dir = unit_root / unit_id
        files = list(unit_dir.iterdir())
        if {path.name for path in files} != expected_unit_files or any(
            path.is_symlink() or not path.is_file() for path in files
        ):
            raise ValueError(f"legacy v5 EAS unit inventory changed: {unit_id}")
        if (
            sha256_file(unit_dir / "completion.json")
            != unit_outputs[unit_id]["completion_sha256"]
            or sha256_file(unit_dir / "trajectory.tsv")
            != unit_outputs[unit_id]["trajectory_sha256"]
        ):
            raise ValueError(f"legacy v5 EAS unit binding changed: {unit_id}")
    return {
        "replay_manifest": replay_manifest,
        "selection_tree": {
            "file_count": tree["file_count"],
            "manifest_sha256": tree["manifest_sha256"],
        },
        "selection_completion": completion,
        "selected": selected,
        "diagnostics": diagnostics,
    }


def _verify_legacy_v5_plan(
    repo_root: Path,
    *,
    legacy_plan_dir: Path,
    proposal_dir: Path,
    legacy_replay_dir: Path,
    legacy_selection_dir: Path,
    natural_run_dir: str | Path | None = None,
    current_v6_workflow_sha256: str | None = None,
) -> dict[str, Any]:
    """Verify the exact 8adf v5 plan and completed EAS artifacts."""

    bundle = legacy_plan_dir.resolve()
    children = list(bundle.iterdir()) if bundle.is_dir() else []
    expected_names = {PLAN_COMPLETION_FILENAME, *LEGACY_V5_PLAN_OUTPUTS.values()}
    if (
        bundle.is_symlink()
        or {path.name for path in children} != expected_names
        or any(path.is_symlink() or not path.is_file() for path in children)
    ):
        raise ValueError("legacy v5 conditioned plan inventory changed")
    completion_path = bundle / PLAN_COMPLETION_FILENAME
    if sha256_file(completion_path) != LEGACY_V5_PLAN_COMPLETION_SHA256:
        raise ValueError("legacy v5 conditioned plan completion hash changed")
    completion = _read_json(completion_path, "legacy v5 plan completion")
    contract = completion.get("contract")
    outputs = completion.get("outputs")
    if (
        completion.get("schema") != LEGACY_V5_PLAN_SCHEMA
        or completion.get("status") != "complete"
        or not isinstance(contract, dict)
        or completion.get("contract_sha256") != _canonical_sha256(contract)
        or completion.get("contract_sha256") != LEGACY_V5_PLAN_CONTRACT_SHA256
        or not isinstance(outputs, dict)
        or set(outputs) != set(LEGACY_V5_PLAN_OUTPUTS)
        or completion.get("output_count") != len(LEGACY_V5_PLAN_OUTPUTS)
    ):
        raise ValueError("legacy v5 conditioned plan completion changed")
    for label, filename in LEGACY_V5_PLAN_OUTPUTS.items():
        path = bundle / filename
        record = outputs[label]
        if (
            record.get("path") != filename
            or path.stat().st_size != int(record.get("size_bytes", -1))
            or sha256_file(path) != record.get("sha256")
        ):
            raise ValueError(f"legacy v5 conditioned plan checksum failed: {label}")
    v4 = _verify_legacy_v4_plan(
        repo_root,
        legacy_plan_dir=bundle.parent / LEGACY_V4_PLAN_DIRNAME,
        legacy_replay_dir=bundle.parent / LEGACY_V4_PARALLEL_REPLAY_DIRNAME,
        proposal_dir=proposal_dir,
        natural_run_dir=natural_run_dir,
        current_assembly_workflow_sha256=LEGACY_V5_WORKFLOW_SHA256,
    )
    base_manifest = _read_json(
        bundle / LEGACY_V5_PLAN_OUTPUTS["eas_base_bank_manifest"],
        "legacy v5 EAS base-bank manifest",
    )
    extension_plan = _read_json(
        bundle / LEGACY_V5_PLAN_OUTPUTS["eas_extension_plan"],
        "legacy v5 EAS extension plan",
    )
    replay_plan_path = bundle / LEGACY_V5_PLAN_OUTPUTS["eas_parallel_replay_plan"]
    replay_plan = _read_json(replay_plan_path, "legacy v5 EAS replay plan")
    expected_replay_plan = _build_eas_parallel_replay_plan(
        v4["attestation"],
        workflow_sha256=LEGACY_V5_WORKFLOW_SHA256,
    )
    if (
        base_manifest != v4["base_manifest"]
        or extension_plan != v4["extension_plan"]
        or sha256_file(replay_plan_path) != LEGACY_V5_REPLAY_PLAN_FILE_SHA256
        or replay_plan.get("plan_payload_sha256")
        != LEGACY_V5_REPLAY_PLAN_PAYLOAD_SHA256
        or replay_plan != expected_replay_plan
        or contract
        != _legacy_v5_plan_contract(
            base_manifest,
            extension_plan,
            legacy_v4_attestation=v4["attestation"],
            replay_plan=replay_plan,
            repo_root=repo_root,
        )
    ):
        raise ValueError("legacy v5 plan or transitive v4 provenance changed")
    for label in LEGACY_V2_PLAN_OUTPUTS:
        record = outputs[label]
        if {
            "size_bytes": int(record["size_bytes"]),
            "sha256": str(record["sha256"]),
        } != v4["attestation"]["legacy_output_bindings"][label]:
            raise ValueError(f"legacy v5 inherited v4 bytes changed: {label}")
    eas_state = _verify_legacy_v5_eas_state(
        repo_root,
        replay_dir=legacy_replay_dir,
        selection_dir=legacy_selection_dir,
    )
    attestation: dict[str, Any] = {
        "schema": "gamma-smc.cosi2-conditioned-gamma-v5-state-attestation/v6",
        "status": "exact_frozen_predecessor_with_complete_eas_state",
        "legacy_plan_dirname": LEGACY_V5_PLAN_DIRNAME,
        "legacy_plan_completion_sha256": LEGACY_V5_PLAN_COMPLETION_SHA256,
        "legacy_plan_contract_sha256": LEGACY_V5_PLAN_CONTRACT_SHA256,
        "legacy_plan_producer_workflow_sha256": LEGACY_V5_WORKFLOW_SHA256,
        "current_v6_workflow_sha256": (
            str(current_v6_workflow_sha256)
            if current_v6_workflow_sha256 is not None
            else sha256_file(MODULE_PATH)
        ),
        "legacy_output_bindings": {
            label: {
                "size_bytes": int(record["size_bytes"]),
                "sha256": str(record["sha256"]),
                **({"rows": int(record["rows"])} if "rows" in record else {}),
            }
            for label, record in outputs.items()
        },
        "legacy_replay_dirname": LEGACY_V5_PARALLEL_REPLAY_DIRNAME,
        "legacy_replay_completion_sha256": LEGACY_V4_REPLAY_COMPLETION_SHA256,
        "legacy_replay_root_manifest_sha256": (LEGACY_V4_REPLAY_ROOT_MANIFEST_SHA256),
        "legacy_selection_dirname": LEGACY_V5_SELECTION_DIRNAME,
        "legacy_selection_completion_sha256": (LEGACY_V5_SELECTION_COMPLETION_SHA256),
        "legacy_selection_contract_sha256": LEGACY_V5_SELECTION_CONTRACT_SHA256,
        "legacy_selection_file_count": LEGACY_V5_SELECTION_ROOT_FILE_COUNT,
        "legacy_selection_root_manifest_sha256": (
            LEGACY_V5_SELECTION_ROOT_MANIFEST_SHA256
        ),
        "legacy_v4_attestation_payload_sha256": v4["attestation"][
            "attestation_payload_sha256"
        ],
        "proposal_replay_or_resampling_in_v6": False,
        "eas_artifacts_rewritten_in_v6": False,
    }
    attestation["attestation_payload_sha256"] = _canonical_sha256(attestation)
    return {
        "attestation": attestation,
        "base_manifest": base_manifest,
        "extension_plan": extension_plan,
        "legacy_replay_plan": replay_plan,
        "selected": eas_state["selected"],
        "diagnostics": eas_state["diagnostics"],
        "selection_completion": eas_state["selection_completion"],
    }


def _recombination_map_binding(repo_root: Path) -> dict[str, Any]:
    source = repo_root / RECOMBINATION_MAP_RELATIVE
    framework_completion_path = repo_root / FRAMEWORK_COMPLETION_RELATIVE
    if (
        source.is_symlink()
        or not source.is_file()
        or source.stat().st_size != RECOMBINATION_MAP_SIZE_BYTES
        or sha256_file(source) != RECOMBINATION_MAP_SHA256
        or framework_completion_path.is_symlink()
        or not framework_completion_path.is_file()
        or sha256_file(framework_completion_path) != FRAMEWORK_COMPLETION_SHA256
    ):
        raise ValueError("canonical CoSi2 recombination map is missing or changed")
    framework_completion = _read_json(
        framework_completion_path, "CoSi2 framework completion"
    )
    framework_record = framework_completion.get("outputs", {}).get(
        "recombination_map.tsv"
    )
    if framework_record != {
        "bytes": RECOMBINATION_MAP_SIZE_BYTES,
        "path": "recombination_map.tsv",
        "sha256": RECOMBINATION_MAP_SHA256,
    }:
        raise ValueError("CoSi2 framework completion no longer binds the map")
    return {
        "source_path": RECOMBINATION_MAP_RELATIVE.as_posix(),
        "size_bytes": RECOMBINATION_MAP_SIZE_BYTES,
        "sha256": RECOMBINATION_MAP_SHA256,
        "framework_completion_path": FRAMEWORK_COMPLETION_RELATIVE.as_posix(),
        "framework_completion_sha256": FRAMEWORK_COMPLETION_SHA256,
        "framework_output_record": framework_record,
        "plan_relative_path": PLAN_OUTPUTS["recombination_map"],
        "cosi_resolution": "recomb_file relative to exact plan working directory",
    }


def _v6_plan_contract(
    base_manifest: Mapping[str, Any],
    extension_plan: Mapping[str, Any],
    *,
    legacy_v5_attestation: Mapping[str, Any],
    recombination_map_binding: Mapping[str, Any],
    repo_root: Path,
    workflow_sha256: str | None = None,
) -> dict[str, Any]:
    contract = _eas_bank_plan_contract(
        base_manifest,
        extension_plan,
        source_binding=_plan_source_binding(repo_root, workflow_sha256=workflow_sha256),
    )
    contract.update(
        {
            "legacy_v5_state_attestation": dict(legacy_v5_attestation),
            "recombination_map_binding": dict(recombination_map_binding),
            "v6_revision_scope": {
                "reason": "CoSi2 resolves recomb_file relative to plan cwd",
                "add_exact_plan_root_recombination_map": True,
                "han_configs_ledger_and_seed_stream_inherited_byte_identically": True,
                "eas_v5_replay_and_selection_reused_byte_identically": True,
                "proposal_replay_or_resampling": False,
                "legacy_plan_replay_selection_bytes_rewritten": False,
            },
        }
    )
    return contract


def _verify_legacy_v6_plan(
    repo_root: Path,
    *,
    legacy_plan_dir: Path,
    proposal_dir: Path,
    natural_run_dir: str | Path | None = None,
) -> dict[str, Any]:
    """Verify the exact map-bearing v6 plan without current-source substitution."""

    bundle = legacy_plan_dir.resolve()
    children = list(bundle.iterdir()) if bundle.is_dir() else []
    expected_names = {PLAN_COMPLETION_FILENAME, *LEGACY_V6_PLAN_OUTPUTS.values()}
    if (
        bundle.is_symlink()
        or {path.name for path in children} != expected_names
        or any(path.is_symlink() or not path.is_file() for path in children)
    ):
        raise ValueError("legacy v6 conditioned plan inventory changed")
    completion_path = bundle / PLAN_COMPLETION_FILENAME
    if sha256_file(completion_path) != LEGACY_V6_PLAN_COMPLETION_SHA256:
        raise ValueError("legacy v6 conditioned plan completion hash changed")
    completion = _read_json(completion_path, "legacy v6 plan completion")
    contract = completion.get("contract")
    outputs = completion.get("outputs")
    if (
        completion.get("schema") != LEGACY_V6_PLAN_SCHEMA
        or completion.get("status") != "complete"
        or not isinstance(contract, dict)
        or completion.get("contract_sha256") != _canonical_sha256(contract)
        or completion.get("contract_sha256") != LEGACY_V6_PLAN_CONTRACT_SHA256
        or not isinstance(outputs, dict)
        or set(outputs) != set(LEGACY_V6_PLAN_OUTPUTS)
        or completion.get("output_count") != len(LEGACY_V6_PLAN_OUTPUTS)
    ):
        raise ValueError("legacy v6 conditioned plan completion changed")
    for label, filename in LEGACY_V6_PLAN_OUTPUTS.items():
        path = bundle / filename
        record = outputs[label]
        if (
            record.get("path") != filename
            or path.stat().st_size != int(record.get("size_bytes", -1))
            or sha256_file(path) != record.get("sha256")
        ):
            raise ValueError(f"legacy v6 conditioned plan checksum failed: {label}")
    v5 = _verify_legacy_v5_plan(
        repo_root,
        legacy_plan_dir=bundle.parent / LEGACY_V5_PLAN_DIRNAME,
        proposal_dir=proposal_dir,
        legacy_replay_dir=bundle.parent / LEGACY_V5_PARALLEL_REPLAY_DIRNAME,
        legacy_selection_dir=bundle.parent / LEGACY_V5_SELECTION_DIRNAME,
        natural_run_dir=natural_run_dir,
        current_v6_workflow_sha256=LEGACY_V6_WORKFLOW_SHA256,
    )
    map_binding = _recombination_map_binding(repo_root)
    expected_contract = _v6_plan_contract(
        v5["base_manifest"],
        v5["extension_plan"],
        legacy_v5_attestation=v5["attestation"],
        recombination_map_binding=map_binding,
        repo_root=repo_root,
        workflow_sha256=LEGACY_V6_WORKFLOW_SHA256,
    )
    if contract != expected_contract:
        raise ValueError("legacy v6 plan transitive provenance changed")
    inherited = v5["attestation"]["legacy_output_bindings"]
    for label, record in inherited.items():
        observed = outputs[label]
        if {
            "size_bytes": int(observed["size_bytes"]),
            "sha256": str(observed["sha256"]),
        } != {
            "size_bytes": int(record["size_bytes"]),
            "sha256": str(record["sha256"]),
        }:
            raise ValueError(f"legacy v6 inherited v5 bytes changed: {label}")
    attestation: dict[str, Any] = {
        "schema": "gamma-smc.cosi2-conditioned-gamma-v6-plan-attestation/v7",
        "status": "exact_frozen_v6_plan",
        "legacy_plan_dirname": LEGACY_V6_PLAN_DIRNAME,
        "legacy_plan_completion_sha256": LEGACY_V6_PLAN_COMPLETION_SHA256,
        "legacy_plan_contract_sha256": LEGACY_V6_PLAN_CONTRACT_SHA256,
        "legacy_plan_producer_workflow_sha256": LEGACY_V6_WORKFLOW_SHA256,
        "current_v7_workflow_sha256": sha256_file(MODULE_PATH),
        "legacy_output_bindings": {
            label: {
                "size_bytes": int(record["size_bytes"]),
                "sha256": str(record["sha256"]),
                **({"rows": int(record["rows"])} if "rows" in record else {}),
            }
            for label, record in outputs.items()
        },
        "legacy_v5_attestation_payload_sha256": v5["attestation"][
            "attestation_payload_sha256"
        ],
        "v6_plan_or_eas_artifacts_rewritten_in_v7": False,
    }
    attestation["attestation_payload_sha256"] = _canonical_sha256(attestation)
    return {
        "completion": completion,
        "attestation": attestation,
        "base_manifest": v5["base_manifest"],
        "extension_plan": v5["extension_plan"],
        "selected": v5["selected"],
        "diagnostics": v5["diagnostics"],
    }


def _v7_plan_contract(
    predecessor: Mapping[str, Any],
    execution_state: Mapping[str, Any],
    *,
    repo_root: Path,
) -> dict[str, Any]:
    contract = _eas_bank_plan_contract(
        predecessor["base_manifest"],
        predecessor["extension_plan"],
        source_binding=_plan_source_binding(repo_root),
    )
    contract.update(
        {
            "legacy_v6_plan_attestation": dict(predecessor["attestation"]),
            "legacy_v6_execution_state": dict(execution_state),
            "v7_revision_scope": {
                "reason": "CoSi2 trajectory loader rejects Han popsize columns",
                "loader_projection_schema": HAN_LOADER_TRAJECTORY_SCHEMA,
                "loader_projection": (
                    "lexical_tab_projection_sim_gen_selfreq_1_through_5_canonical_lf"
                ),
                "full_source_trajectory_retained_and_checksum_validated": True,
                "loader_trajectory_retained_and_checksum_bound": True,
                "failed_attempt_evidence_published_outside_completed_unit_root": True,
                "legacy_v6_plan_han_selection_and_eas_units_rewritten": False,
                "native_v7_simulation_root_absent_at_plan_migration": True,
            },
        }
    )
    return contract


def _generate_plan_v5_legacy_implementation(
    repo_root: str | Path,
    *,
    output_dir: str | Path,
    proposal_dir: str | Path | None = None,
    natural_run_dir: str | Path | None = None,
) -> dict[str, Any]:
    """Create or reuse the assembly-only v5 plan from the exact frozen v4 state."""

    root = Path(repo_root).resolve()
    destination = Path(output_dir).resolve()
    proposal_root = (
        Path(proposal_dir).resolve()
        if proposal_dir is not None
        else destination.parent / EAS_PROPOSAL_DIRNAME
    )
    if destination.exists() and any(destination.iterdir()):
        return verify_plan(
            root,
            plan_dir=destination,
            proposal_dir=proposal_root,
            natural_run_dir=natural_run_dir,
        )
    legacy_plan_root = destination.parent / LEGACY_V4_PLAN_DIRNAME
    migration = _verify_legacy_v4_plan(
        root,
        legacy_plan_dir=legacy_plan_root,
        legacy_replay_dir=(destination.parent / LEGACY_V4_PARALLEL_REPLAY_DIRNAME),
        proposal_dir=proposal_root,
        natural_run_dir=natural_run_dir,
    )
    base_manifest = migration["base_manifest"]
    extension_plan = migration["extension_plan"]
    replay_plan = _build_eas_parallel_replay_plan(migration["attestation"])
    payloads = _plan_payloads(
        root,
        eas_base_bank_manifest=base_manifest,
        eas_extension_plan=extension_plan,
        eas_parallel_replay_plan=replay_plan,
        natural_run_dir=natural_run_dir,
    )
    if destination.exists() and not destination.is_dir():
        raise ValueError("plan destination exists but is not a directory")
    if destination.exists():
        destination.rmdir()
    destination.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(
        tempfile.mkdtemp(prefix=f".{destination.name}.staging.", dir=destination.parent)
    )
    try:
        row_counts: dict[str, int] = {}
        for label, payload in payloads.items():
            path = stage / PLAN_OUTPUTS[label]
            if label in LEGACY_V2_PLAN_OUTPUTS:
                shutil.copyfile(legacy_plan_root / LEGACY_V4_PLAN_OUTPUTS[label], path)
                if isinstance(payload, pd.DataFrame):
                    row_counts[label] = len(payload)
            elif isinstance(payload, dict):
                _atomic_json(path, payload)
            else:
                raise RuntimeError(f"unexpected new v5 plan payload type: {label}")
        outputs = {
            label: _output_record(
                stage / filename,
                stage,
                rows=row_counts.get(label),
            )
            for label, filename in PLAN_OUTPUTS.items()
        }
        contract = _eas_bank_plan_contract(
            base_manifest,
            extension_plan,
            source_binding=_plan_source_binding(root),
        )
        contract.update(
            {
                "eas_parallel_replay_plan_payload_sha256": replay_plan[
                    "plan_payload_sha256"
                ],
                "legacy_v4_plan_attestation": migration["attestation"],
                "replay_only_post_extension_revision": {
                    "proposal_target_inherited_from_pre_extension_v2": True,
                    "proposal_blocks_rewritten": False,
                    "legacy_v2_plan_rewritten": False,
                    "legacy_v3_or_v4_plan_or_replay_rewritten": False,
                    "selected_tsv_float_parser": "pandas_float_precision_round_trip",
                    "canonical_float_round_trip_is_bit_exact": True,
                    "exact_v4_replay_root_migrated_without_replay": True,
                    "selection_unit_id_adapter_fixed": True,
                },
            }
        )
        completion = {
            "schema": PLAN_SCHEMA,
            "status": "complete",
            "contract": contract,
            "contract_sha256": _canonical_sha256(contract),
            "outputs": outputs,
            "output_count": len(outputs),
        }
        _atomic_json(stage / PLAN_COMPLETION_FILENAME, completion)
        expected = {PLAN_COMPLETION_FILENAME, *PLAN_OUTPUTS.values()}
        if {path.name for path in stage.iterdir()} != expected:
            raise ValueError("staged plan has an unexpected inventory")
        os.replace(stage, destination)
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    verified = verify_plan(
        root,
        plan_dir=destination,
        proposal_dir=proposal_root,
        natural_run_dir=natural_run_dir,
    )
    return {**verified, "cache_hit": False}


def _verify_plan_v5_legacy_implementation(
    repo_root: str | Path,
    *,
    plan_dir: str | Path,
    proposal_dir: str | Path | None = None,
    natural_run_dir: str | Path | None = None,
) -> dict[str, Any]:
    """Verify the assembly-only v5 plan and exact frozen v4 predecessor."""

    root = Path(repo_root).resolve()
    bundle = Path(plan_dir).resolve()
    if bundle.is_symlink() or not bundle.is_dir():
        raise ValueError(f"conditioned plan is absent: {bundle}")
    expected_names = {PLAN_COMPLETION_FILENAME, *PLAN_OUTPUTS.values()}
    children = list(bundle.iterdir())
    if (
        any(path.is_symlink() or path.is_dir() for path in children)
        or {path.name for path in children} != expected_names
    ):
        raise ValueError("conditioned plan inventory mismatch")
    completion = _read_json(bundle / PLAN_COMPLETION_FILENAME, "plan completion")
    contract = completion.get("contract")
    proposal_root = (
        Path(proposal_dir).resolve()
        if proposal_dir is not None
        else bundle.parent / EAS_PROPOSAL_DIRNAME
    )
    migration = _verify_legacy_v4_plan(
        root,
        legacy_plan_dir=bundle.parent / LEGACY_V4_PLAN_DIRNAME,
        legacy_replay_dir=bundle.parent / LEGACY_V4_PARALLEL_REPLAY_DIRNAME,
        proposal_dir=proposal_root,
        natural_run_dir=natural_run_dir,
    )
    if (
        completion.get("schema") != PLAN_SCHEMA
        or completion.get("status") != "complete"
        or not isinstance(contract, dict)
        or completion.get("contract_sha256") != _canonical_sha256(contract)
        or contract.get("source_binding") != _plan_source_binding(root)
    ):
        raise ValueError("conditioned plan completion is invalid or stale")
    outputs = completion.get("outputs")
    if (
        not isinstance(outputs, dict)
        or set(outputs) != set(PLAN_OUTPUTS)
        or completion.get("output_count") != len(PLAN_OUTPUTS)
    ):
        raise ValueError("conditioned plan output manifest changed")
    base_manifest = _read_json(
        bundle / PLAN_OUTPUTS["eas_base_bank_manifest"],
        "conditioned plan EAS base-bank manifest",
    )
    extension_plan = _read_json(
        bundle / PLAN_OUTPUTS["eas_extension_plan"],
        "conditioned plan EAS extension plan",
    )
    replay_plan = _read_json(
        bundle / PLAN_OUTPUTS["eas_parallel_replay_plan"],
        "conditioned plan EAS parallel replay plan",
    )
    conditioned.validate_eas_base_bank_manifest(
        root,
        base_manifest,
        proposal_dir=proposal_root,
    )
    conditioned.validate_eas_extension_plan(
        root,
        extension_plan,
        base_manifest=base_manifest,
    )
    _validate_eas_parallel_replay_plan(
        replay_plan,
        legacy_v4_attestation=migration["attestation"],
    )
    expected_contract = _eas_bank_plan_contract(
        base_manifest,
        extension_plan,
        source_binding=_plan_source_binding(root),
    )
    expected_contract.update(
        {
            "eas_parallel_replay_plan_payload_sha256": replay_plan[
                "plan_payload_sha256"
            ],
            "legacy_v4_plan_attestation": migration["attestation"],
            "replay_only_post_extension_revision": {
                "proposal_target_inherited_from_pre_extension_v2": True,
                "proposal_blocks_rewritten": False,
                "legacy_v2_plan_rewritten": False,
                "legacy_v3_or_v4_plan_or_replay_rewritten": False,
                "selected_tsv_float_parser": "pandas_float_precision_round_trip",
                "canonical_float_round_trip_is_bit_exact": True,
                "exact_v4_replay_root_migrated_without_replay": True,
                "selection_unit_id_adapter_fixed": True,
            },
        }
    )
    if contract != expected_contract:
        raise ValueError("conditioned v5 plan replay/bank provenance is stale")
    payloads = _plan_payloads(
        root,
        eas_base_bank_manifest=base_manifest,
        eas_extension_plan=extension_plan,
        eas_parallel_replay_plan=replay_plan,
        natural_run_dir=natural_run_dir,
    )
    for label, filename in PLAN_OUTPUTS.items():
        path = bundle / filename
        record = outputs[label]
        if (
            record.get("path") != filename
            or path.stat().st_size != int(record.get("size_bytes", -1))
            or sha256_file(path) != record.get("sha256")
        ):
            raise ValueError(f"conditioned plan checksum failed: {label}")
        if (
            label in LEGACY_V2_PLAN_OUTPUTS
            and {
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            != migration["attestation"]["legacy_output_bindings"][label]
        ):
            raise ValueError(f"conditioned v5 inherited bytes changed: {label}")
        expected = payloads[label]
        if isinstance(expected, pd.DataFrame):
            observed = pd.read_csv(path, sep="\t")
            pd.testing.assert_frame_equal(
                observed,
                expected,
                check_dtype=False,
                check_exact=False,
                rtol=0,
                atol=1e-15,
            )
            if int(record.get("rows", -1)) != len(expected):
                raise ValueError(f"conditioned plan row count changed: {label}")
        elif isinstance(expected, dict):
            if _read_json(path, f"conditioned plan {label}") != expected:
                raise ValueError(f"conditioned plan JSON changed: {label}")
        elif path.read_text(encoding="utf-8") != expected:
            raise ValueError(f"conditioned plan config changed: {label}")
    return {**completion, "cache_hit": True}


def generate_plan(
    repo_root: str | Path,
    *,
    output_dir: str | Path,
    proposal_dir: str | Path | None = None,
    natural_run_dir: str | Path | None = None,
) -> dict[str, Any]:
    """Create or reuse the Han loader-adapter v7 plan from exact frozen v6 state."""

    root = Path(repo_root).resolve()
    destination = Path(output_dir).resolve()
    proposal_root = (
        Path(proposal_dir).resolve()
        if proposal_dir is not None
        else destination.parent / EAS_PROPOSAL_DIRNAME
    )
    if destination.exists() and any(destination.iterdir()):
        return verify_plan(
            root,
            plan_dir=destination,
            proposal_dir=proposal_root,
            natural_run_dir=natural_run_dir,
        )
    predecessor = _verify_legacy_v6_plan(
        root,
        legacy_plan_dir=destination.parent / LEGACY_V6_PLAN_DIRNAME,
        proposal_dir=proposal_root,
        natural_run_dir=natural_run_dir,
    )
    native_simulation_root = destination.parent / SIMULATION_DIRNAME
    if native_simulation_root.exists() and (
        not native_simulation_root.is_dir() or any(native_simulation_root.iterdir())
    ):
        raise ValueError(
            "v7 simulation root must be absent at immutable plan migration"
        )
    execution_state = _verify_legacy_v6_execution_state(
        root,
        predecessor=predecessor,
        legacy_plan_dir=destination.parent / LEGACY_V6_PLAN_DIRNAME,
        screen_dir=destination.parent / HAN_SCREEN_DIRNAME,
        simulation_dir=destination.parent / LEGACY_V6_SIMULATION_DIRNAME,
        cosi2_binary=root / DEFAULT_COSI2_BINARY_RELATIVE,
    )
    if destination.exists() and not destination.is_dir():
        raise ValueError("v7 plan destination exists but is not a directory")
    if destination.exists():
        destination.rmdir()
    destination.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(
        tempfile.mkdtemp(prefix=f".{destination.name}.staging.", dir=destination.parent)
    )
    try:
        legacy_root = destination.parent / LEGACY_V6_PLAN_DIRNAME
        for label, filename in LEGACY_V6_PLAN_OUTPUTS.items():
            shutil.copyfile(legacy_root / filename, stage / PLAN_OUTPUTS[label])
        outputs: dict[str, dict[str, Any]] = {}
        legacy_outputs = predecessor["attestation"]["legacy_output_bindings"]
        for label, filename in PLAN_OUTPUTS.items():
            record = _output_record(stage / filename, stage)
            if label in legacy_outputs:
                source_record = legacy_outputs[label]
                if "rows" in source_record:
                    record["rows"] = int(source_record["rows"])
            outputs[label] = record
        contract = _v7_plan_contract(
            predecessor,
            execution_state,
            repo_root=root,
        )
        completion = {
            "schema": PLAN_SCHEMA,
            "status": "complete",
            "contract": contract,
            "contract_sha256": _canonical_sha256(contract),
            "outputs": outputs,
            "output_count": len(outputs),
        }
        _atomic_json(stage / PLAN_COMPLETION_FILENAME, completion)
        expected = {PLAN_COMPLETION_FILENAME, *PLAN_OUTPUTS.values()}
        if {path.name for path in stage.iterdir()} != expected:
            raise ValueError("staged v7 plan has an unexpected inventory")
        os.replace(stage, destination)
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    verified = verify_plan(
        root,
        plan_dir=destination,
        proposal_dir=proposal_root,
        natural_run_dir=natural_run_dir,
    )
    return {**verified, "cache_hit": False}


def verify_plan(
    repo_root: str | Path,
    *,
    plan_dir: str | Path,
    proposal_dir: str | Path | None = None,
    natural_run_dir: str | Path | None = None,
) -> dict[str, Any]:
    """Verify v7 loader semantics and the exact frozen v6 execution boundary."""

    root = Path(repo_root).resolve()
    bundle = Path(plan_dir).resolve()
    if bundle.is_symlink() or not bundle.is_dir():
        raise ValueError(f"conditioned v7 plan is absent: {bundle}")
    children = list(bundle.iterdir())
    expected_names = {PLAN_COMPLETION_FILENAME, *PLAN_OUTPUTS.values()}
    if {path.name for path in children} != expected_names or any(
        path.is_symlink() or not path.is_file() for path in children
    ):
        raise ValueError("conditioned v7 plan inventory mismatch")
    proposal_root = (
        Path(proposal_dir).resolve()
        if proposal_dir is not None
        else bundle.parent / EAS_PROPOSAL_DIRNAME
    )
    predecessor = _verify_legacy_v6_plan(
        root,
        legacy_plan_dir=bundle.parent / LEGACY_V6_PLAN_DIRNAME,
        proposal_dir=proposal_root,
        natural_run_dir=natural_run_dir,
    )
    execution_state = _verify_legacy_v6_execution_state(
        root,
        predecessor=predecessor,
        legacy_plan_dir=bundle.parent / LEGACY_V6_PLAN_DIRNAME,
        screen_dir=bundle.parent / HAN_SCREEN_DIRNAME,
        simulation_dir=bundle.parent / LEGACY_V6_SIMULATION_DIRNAME,
        cosi2_binary=root / DEFAULT_COSI2_BINARY_RELATIVE,
    )
    completion = _read_json(bundle / PLAN_COMPLETION_FILENAME, "v7 plan completion")
    contract = completion.get("contract")
    outputs = completion.get("outputs")
    expected_contract = _v7_plan_contract(
        predecessor,
        execution_state,
        repo_root=root,
    )
    map_binding = _recombination_map_binding(root)
    if (
        completion.get("schema") != PLAN_SCHEMA
        or completion.get("status") != "complete"
        or not isinstance(contract, dict)
        or contract != expected_contract
        or completion.get("contract_sha256") != _canonical_sha256(contract)
        or not isinstance(outputs, dict)
        or set(outputs) != set(PLAN_OUTPUTS)
        or completion.get("output_count") != len(PLAN_OUTPUTS)
    ):
        raise ValueError("conditioned v7 plan completion is invalid or stale")
    inherited = predecessor["attestation"]["legacy_output_bindings"]
    for label, filename in PLAN_OUTPUTS.items():
        path = bundle / filename
        record = outputs[label]
        if (
            record.get("path") != filename
            or path.stat().st_size != int(record.get("size_bytes", -1))
            or sha256_file(path) != record.get("sha256")
        ):
            raise ValueError(f"conditioned v7 plan checksum failed: {label}")
        if label in inherited and {
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        } != {
            "size_bytes": int(inherited[label]["size_bytes"]),
            "sha256": str(inherited[label]["sha256"]),
        }:
            raise ValueError(f"conditioned v7 inherited v6 bytes changed: {label}")
    map_path = bundle / PLAN_OUTPUTS["recombination_map"]
    if (
        map_path.stat().st_size != map_binding["size_bytes"]
        or sha256_file(map_path) != map_binding["sha256"]
        or map_path.read_bytes() != (root / RECOMBINATION_MAP_RELATIVE).read_bytes()
    ):
        raise ValueError("conditioned v7 plan recombination map changed")
    return {**completion, "cache_hit": True}


def _eas_block_dir(proposal_root: Path, block_index: int) -> Path:
    return proposal_root / f"block_{block_index:02d}"


def load_eas_proposal_blocks(
    repo_root: str | Path,
    *,
    proposal_dir: str | Path,
    base_manifest: Mapping[str, Any],
    extension_plan: Mapping[str, Any],
) -> list[dict[str, Any]]:
    root = Path(repo_root).resolve()
    proposal_root = Path(proposal_dir).resolve()
    expected_names = {
        f"block_{index:02d}" for index in range(conditioned.EAS_PROPOSAL_BLOCK_COUNT)
    }
    children = list(proposal_root.iterdir()) if proposal_root.is_dir() else []
    if (
        proposal_root.is_symlink()
        or {path.name for path in children} != expected_names
        or any(path.is_symlink() or not path.is_dir() for path in children)
    ):
        raise ValueError(
            "EAS proposal root does not contain the exact fixed-200 block inventory"
        )
    return conditioned.load_eas_extended_proposal_bank(
        root,
        proposal_root,
        base_manifest=base_manifest,
        extension_plan=extension_plan,
    )


def _run_eas_extension_block_process(
    repo_root: str,
    proposal_dir: str,
    block_index: int,
    base_manifest: Mapping[str, Any],
    extension_plan: Mapping[str, Any],
) -> dict[str, Any]:
    """Generate or verify one extension block inside a worker process."""

    started = time.time()
    destination = Path(proposal_dir).resolve()
    unit_dir = _eas_block_dir(destination, int(block_index))
    cache_hit = (unit_dir / "completion.json").is_file()
    block = conditioned.write_eas_proposal_block(
        Path(repo_root).resolve(),
        unit_dir,
        block_index=int(block_index),
        base_manifest=base_manifest,
        extension_plan=extension_plan,
    )
    return {
        "block_index": int(block_index),
        "status": "complete",
        "cache_hit": bool(cache_hit),
        "worker_process_id": os.getpid(),
        "elapsed_seconds": time.time() - started,
        "valid_proposals": int(
            np.count_nonzero(np.asarray(block["outcomes"]) == "valid")
        ),
        "proposal_block_sha256": block["proposal_block_sha256"],
        "completion_sha256": block["completion_sha256"],
    }


def run_eas_proposal_blocks(
    repo_root: str | Path,
    *,
    proposal_dir: str | Path,
    base_manifest: Mapping[str, Any],
    extension_plan: Mapping[str, Any],
    max_workers: int = MAX_WORKERS,
    block_indexes: Sequence[int] | None = None,
) -> pd.DataFrame:
    """Write fixed extension blocks in restartable worker processes."""

    root = Path(repo_root).resolve()
    destination = Path(proposal_dir).resolve()
    workers = _positive_workers(max_workers)
    requested = (
        list(dict.fromkeys(int(value) for value in block_indexes))
        if block_indexes is not None
        else list(
            range(
                conditioned.EAS_BASE_PROPOSAL_BLOCK_COUNT,
                conditioned.EAS_PROPOSAL_BLOCK_COUNT,
            )
        )
    )
    if any(
        value < conditioned.EAS_BASE_PROPOSAL_BLOCK_COUNT
        or value >= conditioned.EAS_PROPOSAL_BLOCK_COUNT
        for value in requested
    ):
        raise ValueError(
            "EAS extension request must lie within immutable block indexes 20..199"
        )
    if not requested:
        return pd.DataFrame(
            columns=(
                "block_index",
                "status",
                "cache_hit",
                "worker_process_id",
                "elapsed_seconds",
                "valid_proposals",
                "proposal_block_sha256",
                "completion_sha256",
            )
        )
    destination.mkdir(parents=True, exist_ok=True)
    conditioned.validate_eas_base_bank_manifest(
        root, base_manifest, proposal_dir=destination
    )
    conditioned.validate_eas_extension_plan(
        root,
        extension_plan,
        base_manifest=base_manifest,
    )
    allowed_names = {
        f"block_{index:02d}" for index in range(conditioned.EAS_PROPOSAL_BLOCK_COUNT)
    }
    children = list(destination.iterdir())
    if any(
        path.name not in allowed_names or path.is_symlink() or not path.is_dir()
        for path in children
    ):
        raise ValueError("EAS proposal root contains an unexpected partial inventory")
    run_started = time.time()
    print(
        json.dumps(
            {
                "event": "eas_extension_process_pool_start",
                "proposal_dir": str(destination),
                "requested_block_count": len(requested),
                "first_requested_block": min(requested),
                "last_requested_block": max(requested),
                "max_worker_processes": min(workers, len(requested)),
                "extension_plan_payload_sha256": extension_plan["plan_payload_sha256"],
            },
            sort_keys=True,
        ),
        file=sys.stderr,
        flush=True,
    )

    rows: list[dict[str, Any]] = []
    failures: list[tuple[int, str]] = []
    with ProcessPoolExecutor(max_workers=min(workers, len(requested))) as pool:
        futures = {
            pool.submit(
                _run_eas_extension_block_process,
                str(root),
                str(destination),
                index,
                dict(base_manifest),
                dict(extension_plan),
            ): index
            for index in requested
        }
        for future in as_completed(futures):
            try:
                row = future.result()
                rows.append(row)
                print(
                    json.dumps(
                        {"event": "eas_extension_block_complete", **row},
                        sort_keys=True,
                    ),
                    file=sys.stderr,
                    flush=True,
                )
            except Exception as error:  # noqa: BLE001 - preserve block failures
                failures.append((futures[future], str(error)))
                print(
                    json.dumps(
                        {
                            "event": "eas_extension_block_failed",
                            "block_index": futures[future],
                            "error": str(error),
                        },
                        sort_keys=True,
                    ),
                    file=sys.stderr,
                    flush=True,
                )
    print(
        json.dumps(
            {
                "event": "eas_extension_process_pool_finished",
                "requested_block_count": len(requested),
                "completed_block_count": len(rows),
                "failed_block_count": len(failures),
                "elapsed_seconds": time.time() - run_started,
            },
            sort_keys=True,
        ),
        file=sys.stderr,
        flush=True,
    )
    if failures:
        detail = "; ".join(f"block {index}: {error}" for index, error in failures)
        raise RuntimeError(f"{len(failures)} EAS proposal blocks failed: {detail}")
    return pd.DataFrame(rows).sort_values("block_index").reset_index(drop=True)


def _eas_parallel_replay_block_dir(replay_root: Path, block_index: int) -> Path:
    return replay_root / f"block_{int(block_index):02d}"


def _selected_block_frame(
    selected: pd.DataFrame,
    block_index: int,
) -> pd.DataFrame:
    group = selected[
        pd.to_numeric(selected["block_index"], errors="raise").astype(int)
        == int(block_index)
    ].copy()
    if group.empty or set(pd.to_numeric(group["block_index"]).astype(int)) != {
        int(block_index)
    }:
        raise ValueError(f"selected EAS replay block group is empty: {block_index}")
    return group.reset_index(drop=True)


def _read_eas_selected_frame(path: Path) -> pd.DataFrame:
    return pd.read_csv(
        path,
        sep="\t",
        float_precision="round_trip",
        dtype={
            "unit_id": "string",
            "selection_method": "string",
            "source_proposal_block_sha256": "string",
            "source_proposal_completion_sha256": "string",
        },
    )


def _validate_selected_group_against_block(
    selected_group: pd.DataFrame,
    block: Mapping[str, Any],
) -> None:
    required = {
        "accepted_index",
        "unit_id",
        "cosi2_seed",
        "block_index",
        "block_seed",
        "proposal_index",
        "global_proposal_index",
        "birth_generation",
        "log_forward",
        "log_reverse",
        "log_mutation",
        "log_weight",
        "source_proposal_block_sha256",
        "source_proposal_completion_sha256",
    }
    if (
        not required.issubset(selected_group.columns)
        or selected_group["unit_id"].duplicated().any()
    ):
        raise ValueError("parallel EAS selected block group lacks provenance")
    for row in selected_group.to_dict(orient="records"):
        conditioned._eas_selected_particle_record(row)
    block_index = int(block["block_index"])
    proposal_count = int(block["proposal_count"])
    proposal_indexes = pd.to_numeric(
        selected_group["proposal_index"], errors="raise"
    ).to_numpy(dtype=np.int64)
    global_indexes = pd.to_numeric(
        selected_group["global_proposal_index"], errors="raise"
    ).to_numpy(dtype=np.int64)
    if (
        set(pd.to_numeric(selected_group["block_index"]).astype(int)) != {block_index}
        or set(pd.to_numeric(selected_group["block_seed"]).astype(int))
        != {int(block["block_seed"])}
        or np.any((proposal_indexes < 0) | (proposal_indexes >= proposal_count))
        or not np.array_equal(
            global_indexes,
            block_index * proposal_count + proposal_indexes,
        )
        or set(selected_group["source_proposal_block_sha256"].astype(str))
        != {str(block["proposal_block_sha256"])}
        or set(selected_group["source_proposal_completion_sha256"].astype(str))
        != {str(block["completion_sha256"])}
    ):
        raise ValueError("parallel EAS selected group/source block provenance changed")
    outcomes = np.asarray(block["outcomes"]).astype(str)[proposal_indexes]
    if set(outcomes) != {"valid"}:
        raise ValueError("parallel EAS selected group includes an invalid proposal")
    comparisons = {
        "birth_generation": np.asarray(block["birth_generations"])[proposal_indexes],
        "log_forward": np.asarray(block["log_forward"])[proposal_indexes],
        "log_reverse": np.asarray(block["log_reverse"])[proposal_indexes],
        "log_mutation": np.asarray(block["log_mutation"])[proposal_indexes],
        "log_weight": np.asarray(block["log_weights"])[proposal_indexes],
    }
    for label, observed in comparisons.items():
        expected = pd.to_numeric(selected_group[label], errors="raise").to_numpy()
        if label == "birth_generation":
            equal = np.array_equal(observed.astype(np.int64), expected.astype(np.int64))
        else:
            equal = np.array_equal(observed.astype(float), expected.astype(float))
        if not equal:
            raise ValueError(f"parallel EAS selected group changed {label}")


def _eas_parallel_replay_block_contract(
    *,
    selected_group: pd.DataFrame,
    block: Mapping[str, Any],
    plan_completion_sha256: str,
    parallel_replay_plan: Mapping[str, Any],
    base_manifest: Mapping[str, Any],
    extension_plan: Mapping[str, Any],
) -> dict[str, Any]:
    block_index = int(block["block_index"])
    _validate_selected_group_against_block(selected_group, block)
    return {
        "schema": EAS_PARALLEL_REPLAY_BLOCK_SCHEMA,
        "execution_method": EAS_PARALLEL_REPLAY_METHOD,
        "worker_helper": "_write_eas_parallel_replay_block",
        "workflow_sha256": sha256_file(MODULE_PATH),
        "core_sha256": sha256_file(conditioned.MODULE_PATH),
        "plan_completion_sha256": str(plan_completion_sha256),
        "legacy_v2_plan_completion_sha256": LEGACY_V2_PLAN_COMPLETION_SHA256,
        "legacy_v3_plan_completion_sha256": LEGACY_V3_PLAN_COMPLETION_SHA256,
        "parallel_replay_plan_payload_sha256": parallel_replay_plan[
            "plan_payload_sha256"
        ],
        "base_manifest_payload_sha256": base_manifest["manifest_payload_sha256"],
        "extension_plan_payload_sha256": extension_plan["plan_payload_sha256"],
        "block_index": block_index,
        "block_seed": int(block["block_seed"]),
        "proposal_count": int(block["proposal_count"]),
        "source_proposal_block_sha256": str(block["proposal_block_sha256"]),
        "source_proposal_completion_sha256": str(block["completion_sha256"]),
        "selected_group_sha256": _portable_frame_sha256(selected_group),
        "selected_occurrence_count": len(selected_group),
        "selected_unique_proposal_count": int(
            selected_group["proposal_index"].nunique()
        ),
        "selected_unit_ids_in_original_order": selected_group["unit_id"]
        .astype(str)
        .tolist(),
        "artifact_provenance_kind": "current_v4",
        "selected_tsv_float_parser": "pandas_float_precision_round_trip",
        "canonical_selected_tsv_sha256": _portable_frame_sha256(selected_group),
        "core_replay_validation": (
            "full_outcomes_and_numeric_arrays_equal_with_nan_then_selected_rows_"
            "exact_then_full_path_log_components_abs_tolerance_1e-9"
        ),
    }


def _legacy_v3_parallel_replay_block_contract(
    *,
    selected_group: pd.DataFrame,
    block: Mapping[str, Any],
    legacy_parallel_replay_plan: Mapping[str, Any],
    base_manifest: Mapping[str, Any],
    extension_plan: Mapping[str, Any],
) -> dict[str, Any]:
    block_index = int(block["block_index"])
    _validate_selected_group_against_block(selected_group, block)
    if (
        legacy_parallel_replay_plan.get("plan_payload_sha256")
        != LEGACY_V3_REPLAY_PLAN_PAYLOAD_SHA256
    ):
        raise ValueError("legacy v3 replay-plan payload changed")
    return {
        "schema": LEGACY_V3_PARALLEL_REPLAY_BLOCK_SCHEMA,
        "execution_method": LEGACY_V4_PARALLEL_REPLAY_METHOD,
        "worker_helper": "_write_eas_parallel_replay_block",
        "workflow_sha256": LEGACY_V3_WORKFLOW_SHA256,
        "core_sha256": sha256_file(conditioned.MODULE_PATH),
        "plan_completion_sha256": LEGACY_V3_PLAN_COMPLETION_SHA256,
        "legacy_v2_plan_completion_sha256": LEGACY_V2_PLAN_COMPLETION_SHA256,
        "parallel_replay_plan_payload_sha256": (LEGACY_V3_REPLAY_PLAN_PAYLOAD_SHA256),
        "base_manifest_payload_sha256": base_manifest["manifest_payload_sha256"],
        "extension_plan_payload_sha256": extension_plan["plan_payload_sha256"],
        "block_index": block_index,
        "block_seed": int(block["block_seed"]),
        "proposal_count": int(block["proposal_count"]),
        "source_proposal_block_sha256": str(block["proposal_block_sha256"]),
        "source_proposal_completion_sha256": str(block["completion_sha256"]),
        "selected_group_sha256": _portable_frame_sha256(selected_group),
        "selected_occurrence_count": len(selected_group),
        "selected_unique_proposal_count": int(
            selected_group["proposal_index"].nunique()
        ),
        "selected_unit_ids_in_original_order": selected_group["unit_id"]
        .astype(str)
        .tolist(),
        "core_replay_validation": (
            "full_outcomes_and_numeric_arrays_equal_with_nan_then_selected_rows_"
            "exact_then_full_path_log_components_abs_tolerance_1e-9"
        ),
    }


def _legacy_v4_parallel_replay_block_contract(
    *,
    selected_group: pd.DataFrame,
    block: Mapping[str, Any],
    legacy_parallel_replay_plan: Mapping[str, Any],
    base_manifest: Mapping[str, Any],
    extension_plan: Mapping[str, Any],
) -> dict[str, Any]:
    """Reconstruct the frozen 5e08 v4 block contract exactly."""

    block_index = int(block["block_index"])
    _validate_selected_group_against_block(selected_group, block)
    if (
        legacy_parallel_replay_plan.get("plan_payload_sha256")
        != LEGACY_V4_REPLAY_PLAN_PAYLOAD_SHA256
    ):
        raise ValueError("legacy v4 replay-plan payload changed")
    return {
        "schema": LEGACY_V4_PARALLEL_REPLAY_BLOCK_SCHEMA,
        "execution_method": LEGACY_V4_PARALLEL_REPLAY_METHOD,
        "worker_helper": "_write_eas_parallel_replay_block",
        "workflow_sha256": LEGACY_V4_WORKFLOW_SHA256,
        "core_sha256": sha256_file(conditioned.MODULE_PATH),
        "plan_completion_sha256": LEGACY_V4_PLAN_COMPLETION_SHA256,
        "legacy_v2_plan_completion_sha256": LEGACY_V2_PLAN_COMPLETION_SHA256,
        "legacy_v3_plan_completion_sha256": LEGACY_V3_PLAN_COMPLETION_SHA256,
        "parallel_replay_plan_payload_sha256": (LEGACY_V4_REPLAY_PLAN_PAYLOAD_SHA256),
        "base_manifest_payload_sha256": base_manifest["manifest_payload_sha256"],
        "extension_plan_payload_sha256": extension_plan["plan_payload_sha256"],
        "block_index": block_index,
        "block_seed": int(block["block_seed"]),
        "proposal_count": int(block["proposal_count"]),
        "source_proposal_block_sha256": str(block["proposal_block_sha256"]),
        "source_proposal_completion_sha256": str(block["completion_sha256"]),
        "selected_group_sha256": _portable_frame_sha256(selected_group),
        "selected_occurrence_count": len(selected_group),
        "selected_unique_proposal_count": int(
            selected_group["proposal_index"].nunique()
        ),
        "selected_unit_ids_in_original_order": selected_group["unit_id"]
        .astype(str)
        .tolist(),
        "artifact_provenance_kind": "current_v4",
        "selected_tsv_float_parser": "pandas_float_precision_round_trip",
        "canonical_selected_tsv_sha256": _portable_frame_sha256(selected_group),
        "core_replay_validation": (
            "full_outcomes_and_numeric_arrays_equal_with_nan_then_selected_rows_"
            "exact_then_full_path_log_components_abs_tolerance_1e-9"
        ),
    }


def _pack_selected_paths(
    selected_group: pd.DataFrame,
    paths: Mapping[str, Sequence[int]],
) -> dict[str, np.ndarray]:
    unit_ids = selected_group["unit_id"].astype(str).tolist()
    if set(paths) != set(unit_ids):
        raise RuntimeError("parallel EAS replay path inventory is incomplete")
    offsets = [0]
    flattened: list[int] = []
    for unit_id in unit_ids:
        path = tuple(int(value) for value in paths[unit_id])
        flattened.extend(path)
        offsets.append(len(flattened))
    return {
        "unit_ids": np.asarray(unit_ids, dtype=str),
        "offsets": np.asarray(offsets, dtype=np.int64),
        "derived_counts": np.asarray(flattened, dtype=np.int64),
    }


def _unpack_selected_paths(
    selected_group: pd.DataFrame,
    arrays: Mapping[str, np.ndarray],
) -> dict[str, tuple[int, ...]]:
    if set(arrays) != {"unit_ids", "offsets", "derived_counts"}:
        raise ValueError("parallel EAS replay path-array inventory changed")
    expected_ids = selected_group["unit_id"].astype(str).tolist()
    observed_ids = np.asarray(arrays["unit_ids"]).astype(str).tolist()
    offsets = np.asarray(arrays["offsets"], dtype=np.int64)
    counts = np.asarray(arrays["derived_counts"], dtype=np.int64)
    if (
        observed_ids != expected_ids
        or offsets.shape != (len(expected_ids) + 1,)
        or int(offsets[0]) != 0
        or int(offsets[-1]) != len(counts)
        or np.any(np.diff(offsets) <= 0)
    ):
        raise ValueError("parallel EAS replay packed paths changed")
    return {
        unit_id: tuple(counts[offsets[index] : offsets[index + 1]].tolist())
        for index, unit_id in enumerate(expected_ids)
    }


def _verify_eas_parallel_replay_block_contents(
    repo_root: str | Path,
    *,
    bundle_dir: str | Path,
    selected_group: pd.DataFrame,
    expected_schema: str,
    expected_contract: Mapping[str, Any],
    provenance_kind: str,
) -> dict[str, Any]:
    """Verify exact replay-bundle bytes and paths under a supplied contract."""

    root = Path(repo_root).resolve()
    bundle = Path(bundle_dir).resolve()
    block_indexes = set(
        pd.to_numeric(selected_group["block_index"], errors="raise").astype(int)
    )
    if len(block_indexes) != 1:
        raise ValueError("parallel EAS replay bundle must contain one block")
    block_index = next(iter(block_indexes))
    expected_names = {"completion.json", "selected_rows.tsv", "paths.npz"}
    children = list(bundle.iterdir()) if bundle.is_dir() else []
    if (
        bundle.is_symlink()
        or {path.name for path in children} != expected_names
        or any(path.is_symlink() or not path.is_file() for path in children)
    ):
        raise ValueError(f"parallel EAS replay block inventory changed: {block_index}")
    completion = _read_json(bundle / "completion.json", "parallel replay block")
    if (
        completion.get("schema") != expected_schema
        or completion.get("status") != "complete"
        or completion.get("contract") != expected_contract
        or completion.get("contract_sha256") != _canonical_sha256(expected_contract)
    ):
        raise ValueError(f"parallel EAS replay block contract changed: {block_index}")
    outputs = completion.get("outputs")
    if not isinstance(outputs, dict) or set(outputs) != {
        "selected_rows",
        "paths",
    }:
        raise ValueError("parallel EAS replay block output manifest changed")
    for label, filename in {
        "selected_rows": "selected_rows.tsv",
        "paths": "paths.npz",
    }.items():
        path = bundle / filename
        record = outputs[label]
        if (
            record.get("path") != filename
            or path.stat().st_size != int(record.get("size_bytes", -1))
            or sha256_file(path) != record.get("sha256")
        ):
            raise ValueError(f"parallel EAS replay checksum failed: {block_index}")
    observed_selected = _read_eas_selected_frame(bundle / "selected_rows.tsv")
    if sha256_file(bundle / "selected_rows.tsv") != expected_contract.get(
        "selected_group_sha256"
    ):
        raise ValueError("parallel EAS canonical selected-row bytes changed")
    pd.testing.assert_frame_equal(
        observed_selected,
        selected_group,
        check_dtype=False,
        check_exact=True,
    )
    if int(outputs["selected_rows"].get("rows", -1)) != len(selected_group):
        raise ValueError("parallel EAS replay selected-row count changed")
    with np.load(bundle / "paths.npz", allow_pickle=False) as archive:
        arrays = {name: np.asarray(archive[name]) for name in archive.files}
    paths = _unpack_selected_paths(selected_group, arrays)
    integer_ne = conditioned.eas_integer_ne_by_generation(root)
    for row in selected_group.to_dict(orient="records"):
        path = paths[str(row["unit_id"])]
        if len(path) != int(row["birth_generation"]) + 1 or path[-1] != 1:
            raise ValueError("parallel EAS replay selected path changed")
        conditioned.eas_selected_path_frames(integer_ne, path)
        recomputed = conditioned._reverse_path_log_components(integer_ne, path)
        for label, observed in recomputed.items():
            if not math.isclose(
                float(observed),
                float(row[label]),
                rel_tol=0.0,
                abs_tol=1e-9,
            ):
                raise ValueError(
                    f"parallel EAS cached replay path changed {label}: {row['unit_id']}"
                )
    return {
        "status": "complete",
        "cache_hit": True,
        "block_index": block_index,
        "selected_occurrence_count": len(selected_group),
        "selected_unique_proposal_count": int(
            selected_group["proposal_index"].nunique()
        ),
        "completion_sha256": sha256_file(bundle / "completion.json"),
        "artifact_provenance_kind": provenance_kind,
        "source_replay_plan_completion_sha256": str(
            expected_contract["plan_completion_sha256"]
        ),
        "original_legacy_completion_sha256": (
            sha256_file(bundle / "completion.json")
            if provenance_kind == "exact_legacy_v3_migration"
            else ""
        ),
        "paths": paths,
    }


def verify_eas_parallel_replay_block(
    repo_root: str | Path,
    *,
    proposal_dir: str | Path,
    bundle_dir: str | Path,
    selected_group: pd.DataFrame,
    base_manifest: Mapping[str, Any],
    extension_plan: Mapping[str, Any],
    parallel_replay_plan: Mapping[str, Any],
    plan_completion_sha256: str,
) -> dict[str, Any]:
    """Verify one current-schema restartable selected-block replay bundle."""

    block_indexes = set(
        pd.to_numeric(selected_group["block_index"], errors="raise").astype(int)
    )
    if len(block_indexes) != 1:
        raise ValueError("parallel EAS replay bundle must contain one block")
    block_index = next(iter(block_indexes))
    block = conditioned.verify_eas_proposal_block(
        Path(repo_root).resolve(),
        _eas_block_dir(Path(proposal_dir).resolve(), block_index),
        expected_block_index=block_index,
        base_manifest=base_manifest,
        extension_plan=extension_plan,
    )
    expected_contract = _eas_parallel_replay_block_contract(
        selected_group=selected_group,
        block=block,
        plan_completion_sha256=plan_completion_sha256,
        parallel_replay_plan=parallel_replay_plan,
        base_manifest=base_manifest,
        extension_plan=extension_plan,
    )
    return _verify_eas_parallel_replay_block_contents(
        repo_root,
        bundle_dir=bundle_dir,
        selected_group=selected_group,
        expected_schema=EAS_PARALLEL_REPLAY_BLOCK_SCHEMA,
        expected_contract=expected_contract,
        provenance_kind="current_v4",
    )


def verify_legacy_v3_parallel_replay_block(
    repo_root: str | Path,
    *,
    proposal_dir: str | Path,
    bundle_dir: str | Path,
    selected_group: pd.DataFrame,
    base_manifest: Mapping[str, Any],
    extension_plan: Mapping[str, Any],
    legacy_parallel_replay_plan: Mapping[str, Any],
    expected_bundle_binding: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Verify one exact unpublished v3 replay bundle under its frozen contract."""

    block_indexes = set(
        pd.to_numeric(selected_group["block_index"], errors="raise").astype(int)
    )
    if len(block_indexes) != 1:
        raise ValueError("legacy v3 replay bundle must contain one block")
    block_index = next(iter(block_indexes))
    if block_index not in LEGACY_V3_REPLAY_BLOCK_INDEXES:
        raise ValueError("legacy v3 replay block is not in the frozen cache manifest")
    block = conditioned.verify_eas_proposal_block(
        Path(repo_root).resolve(),
        _eas_block_dir(Path(proposal_dir).resolve(), block_index),
        expected_block_index=block_index,
        base_manifest=base_manifest,
        extension_plan=extension_plan,
    )
    expected_contract = _legacy_v3_parallel_replay_block_contract(
        selected_group=selected_group,
        block=block,
        legacy_parallel_replay_plan=legacy_parallel_replay_plan,
        base_manifest=base_manifest,
        extension_plan=extension_plan,
    )
    if expected_bundle_binding is not None:
        bundle = Path(bundle_dir).resolve()
        observed_binding = {
            "block_dir": f"block_{block_index:02d}",
            "files": {
                path.name: {
                    "size_bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
                for path in sorted(bundle.iterdir(), key=lambda value: value.name)
            },
        }
        if observed_binding != dict(expected_bundle_binding):
            raise ValueError("legacy v3 replay bundle differs from frozen manifest")
    return _verify_eas_parallel_replay_block_contents(
        repo_root,
        bundle_dir=bundle_dir,
        selected_group=selected_group,
        expected_schema=LEGACY_V3_PARALLEL_REPLAY_BLOCK_SCHEMA,
        expected_contract=expected_contract,
        provenance_kind="exact_legacy_v3_migration",
    )


def verify_legacy_v4_parallel_replay_block(
    repo_root: str | Path,
    *,
    proposal_dir: str | Path,
    bundle_dir: str | Path,
    selected_group: pd.DataFrame,
    base_manifest: Mapping[str, Any],
    extension_plan: Mapping[str, Any],
    legacy_parallel_replay_plan: Mapping[str, Any],
) -> dict[str, Any]:
    """Verify one frozen v4-native replay block without current reinterpretation."""

    block_indexes = set(
        pd.to_numeric(selected_group["block_index"], errors="raise").astype(int)
    )
    if len(block_indexes) != 1:
        raise ValueError("legacy v4 replay bundle must contain one block")
    block_index = next(iter(block_indexes))
    block = conditioned.verify_eas_proposal_block(
        Path(repo_root).resolve(),
        _eas_block_dir(Path(proposal_dir).resolve(), block_index),
        expected_block_index=block_index,
        base_manifest=base_manifest,
        extension_plan=extension_plan,
    )
    expected_contract = _legacy_v4_parallel_replay_block_contract(
        selected_group=selected_group,
        block=block,
        legacy_parallel_replay_plan=legacy_parallel_replay_plan,
        base_manifest=base_manifest,
        extension_plan=extension_plan,
    )
    return _verify_eas_parallel_replay_block_contents(
        repo_root,
        bundle_dir=bundle_dir,
        selected_group=selected_group,
        expected_schema=LEGACY_V4_PARALLEL_REPLAY_BLOCK_SCHEMA,
        expected_contract=expected_contract,
        provenance_kind="exact_legacy_v4_root_migration",
    )


def _verify_exact_legacy_v4_replay_root_semantics(
    repo_root: str | Path,
    *,
    proposal_dir: str | Path,
    replay_dir: str | Path,
    selected: pd.DataFrame,
    base_manifest: Mapping[str, Any],
    extension_plan: Mapping[str, Any],
    legacy_v4_parallel_replay_plan: Mapping[str, Any],
    legacy_v3_parallel_replay_plan: Mapping[str, Any],
) -> dict[str, Any]:
    """Deep-verify the hard-bound v4 root, its blocks, and selected paths."""

    root = Path(repo_root).resolve()
    replay_root = Path(replay_dir).resolve()
    manifest = _legacy_v4_replay_root_manifest(replay_root)
    block_indexes = sorted(
        set(pd.to_numeric(selected["block_index"], errors="raise").astype(int))
    )
    if block_indexes != manifest["legacy_block_indexes"]:
        raise ValueError("legacy v4 replay block indexes differ from selection")
    paths_by_unit: dict[str, tuple[int, ...]] = {}
    semantic_bindings: list[dict[str, Any]] = []
    for block_index in block_indexes:
        bundle = _eas_parallel_replay_block_dir(replay_root, block_index)
        selected_group = _selected_block_frame(selected, block_index)
        block_completion = _read_json(
            bundle / "completion.json", "legacy v4 replay block"
        )
        if block_completion.get("schema") == LEGACY_V3_PARALLEL_REPLAY_BLOCK_SCHEMA:
            verified = verify_legacy_v3_parallel_replay_block(
                root,
                proposal_dir=proposal_dir,
                bundle_dir=bundle,
                selected_group=selected_group,
                base_manifest=base_manifest,
                extension_plan=extension_plan,
                legacy_parallel_replay_plan=legacy_v3_parallel_replay_plan,
                expected_bundle_binding=_legacy_v3_replay_bundle_binding(
                    legacy_v4_parallel_replay_plan, block_index
                ),
            )
        elif block_completion.get("schema") == LEGACY_V4_PARALLEL_REPLAY_BLOCK_SCHEMA:
            verified = verify_legacy_v4_parallel_replay_block(
                root,
                proposal_dir=proposal_dir,
                bundle_dir=bundle,
                selected_group=selected_group,
                base_manifest=base_manifest,
                extension_plan=extension_plan,
                legacy_parallel_replay_plan=legacy_v4_parallel_replay_plan,
            )
        else:
            raise ValueError(f"legacy v4 replay block schema changed: {block_index}")
        overlap = set(paths_by_unit).intersection(verified["paths"])
        if overlap:
            raise ValueError("legacy v4 replay unit occurs in multiple blocks")
        paths_by_unit.update(verified["paths"])
        semantic_bindings.append(
            {
                "block_index": block_index,
                "selected_occurrence_count": verified["selected_occurrence_count"],
                "selected_unique_proposal_count": verified[
                    "selected_unique_proposal_count"
                ],
                "replay_completion_sha256": verified["completion_sha256"],
                "artifact_provenance_kind": (
                    "exact_legacy_v3_migration"
                    if block_completion.get("schema")
                    == LEGACY_V3_PARALLEL_REPLAY_BLOCK_SCHEMA
                    else "current_v4"
                ),
                "source_replay_plan_completion_sha256": (
                    LEGACY_V3_PLAN_COMPLETION_SHA256
                    if block_completion.get("schema")
                    == LEGACY_V3_PARALLEL_REPLAY_BLOCK_SCHEMA
                    else LEGACY_V4_PLAN_COMPLETION_SHA256
                ),
                "original_legacy_completion_sha256": (
                    verified["completion_sha256"]
                    if block_completion.get("schema")
                    == LEGACY_V3_PARALLEL_REPLAY_BLOCK_SCHEMA
                    else ""
                ),
            }
        )
    completion = _read_json(
        replay_root / EAS_PARALLEL_REPLAY_COMPLETION_FILENAME,
        "legacy v4 replay completion",
    )
    contract = completion["contract"]
    expected_ids = selected["unit_id"].astype(str).tolist()
    if (
        contract.get("selected_block_bindings") != semantic_bindings
        or contract.get("selected_occurrences_sha256")
        != _portable_frame_sha256(selected)
        or contract.get("selected_unit_ids_in_original_order") != expected_ids
        or set(paths_by_unit) != set(expected_ids)
        or contract.get("core_sha256") != sha256_file(conditioned.MODULE_PATH)
        or contract.get("legacy_v3_replay_manifest_sha256")
        != LEGACY_V3_REPLAY_MANIFEST_SHA256
    ):
        raise ValueError("legacy v4 replay root semantic contract changed")
    return {
        "status": "complete",
        "cache_hit": True,
        "completion_sha256": LEGACY_V4_REPLAY_COMPLETION_SHA256,
        "unique_selected_block_count": len(block_indexes),
        "block_bindings": semantic_bindings,
        "paths": {unit_id: paths_by_unit[unit_id] for unit_id in expected_ids},
        "legacy_v4_replay_root_manifest": manifest,
        "completion_interpretation": "hard_bound_legacy_v4_only",
    }


def _legacy_v3_replay_bundle_binding(
    parallel_replay_plan: Mapping[str, Any],
    block_index: int,
) -> dict[str, Any]:
    attestation = parallel_replay_plan.get("legacy_v3_plan_attestation", {})
    manifest = attestation.get("legacy_replay_manifest", {})
    entries = manifest.get("legacy_block_bindings", [])
    expected_name = f"block_{int(block_index):02d}"
    matches = [entry for entry in entries if entry.get("block_dir") == expected_name]
    if len(matches) != 1:
        raise ValueError("legacy v3 replay binding is absent from the v4 plan")
    return dict(matches[0])


def _write_eas_parallel_replay_block(
    repo_root: str,
    proposal_dir: str,
    bundle_dir: str,
    selected_records: Sequence[Mapping[str, Any]],
    base_manifest: Mapping[str, Any],
    extension_plan: Mapping[str, Any],
    parallel_replay_plan: Mapping[str, Any],
    plan_completion_sha256: str,
    legacy_v3_bundle_dir: str,
    legacy_v3_parallel_replay_plan: Mapping[str, Any],
) -> dict[str, Any]:
    """Replay one selected source block, exactly validate it, and persist paths."""

    started = time.time()
    root = Path(repo_root).resolve()
    destination = Path(bundle_dir).resolve()
    selected_group = pd.DataFrame.from_records(selected_records)
    block_indexes = set(
        pd.to_numeric(selected_group["block_index"], errors="raise").astype(int)
    )
    if len(block_indexes) != 1:
        raise ValueError("parallel EAS replay worker received multiple blocks")
    block_index = next(iter(block_indexes))
    if (destination / "completion.json").is_file():
        observed_completion = _read_json(
            destination / "completion.json", "parallel replay block"
        )
        if observed_completion.get("schema") == LEGACY_V3_PARALLEL_REPLAY_BLOCK_SCHEMA:
            result = verify_legacy_v3_parallel_replay_block(
                root,
                proposal_dir=proposal_dir,
                bundle_dir=destination,
                selected_group=selected_group,
                base_manifest=base_manifest,
                extension_plan=extension_plan,
                legacy_parallel_replay_plan=legacy_v3_parallel_replay_plan,
                expected_bundle_binding=_legacy_v3_replay_bundle_binding(
                    parallel_replay_plan, block_index
                ),
            )
        else:
            result = verify_eas_parallel_replay_block(
                root,
                proposal_dir=proposal_dir,
                bundle_dir=destination,
                selected_group=selected_group,
                base_manifest=base_manifest,
                extension_plan=extension_plan,
                parallel_replay_plan=parallel_replay_plan,
                plan_completion_sha256=plan_completion_sha256,
            )
        return {key: value for key, value in result.items() if key != "paths"} | {
            "worker_process_id": os.getpid(),
            "elapsed_seconds": time.time() - started,
        }
    if destination.exists() and any(destination.iterdir()):
        raise ValueError("parallel EAS replay destination is incomplete and nonempty")
    legacy_source = Path(legacy_v3_bundle_dir).resolve()
    if (legacy_source / "completion.json").is_file():
        source = verify_legacy_v3_parallel_replay_block(
            root,
            proposal_dir=proposal_dir,
            bundle_dir=legacy_source,
            selected_group=selected_group,
            base_manifest=base_manifest,
            extension_plan=extension_plan,
            legacy_parallel_replay_plan=legacy_v3_parallel_replay_plan,
            expected_bundle_binding=_legacy_v3_replay_bundle_binding(
                parallel_replay_plan, block_index
            ),
        )
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists():
            destination.rmdir()
        stage = Path(
            tempfile.mkdtemp(
                prefix=f".{destination.name}.legacy-v3-staging.",
                dir=destination.parent,
            )
        )
        try:
            for filename in ("completion.json", "selected_rows.tsv", "paths.npz"):
                shutil.copyfile(legacy_source / filename, stage / filename)
                if sha256_file(stage / filename) != sha256_file(
                    legacy_source / filename
                ):
                    raise RuntimeError("legacy v3 replay atomic copy changed bytes")
            os.replace(stage, destination)
        finally:
            if stage.exists():
                shutil.rmtree(stage, ignore_errors=True)
        copied = verify_legacy_v3_parallel_replay_block(
            root,
            proposal_dir=proposal_dir,
            bundle_dir=destination,
            selected_group=selected_group,
            base_manifest=base_manifest,
            extension_plan=extension_plan,
            legacy_parallel_replay_plan=legacy_v3_parallel_replay_plan,
            expected_bundle_binding=_legacy_v3_replay_bundle_binding(
                parallel_replay_plan, block_index
            ),
        )
        if copied["completion_sha256"] != source["completion_sha256"]:
            raise RuntimeError("legacy v3 replay migration completion changed")
        return {key: value for key, value in copied.items() if key != "paths"} | {
            "cache_hit": True,
            "legacy_v3_migration_cache_hit": True,
            "worker_process_id": os.getpid(),
            "elapsed_seconds": time.time() - started,
        }
    block = conditioned.verify_eas_proposal_block(
        root,
        _eas_block_dir(Path(proposal_dir).resolve(), block_index),
        expected_block_index=block_index,
        base_manifest=base_manifest,
        extension_plan=extension_plan,
    )
    integer_ne = conditioned.eas_integer_ne_by_generation(root)
    paths = conditioned.replay_eas_global_selection(
        integer_ne,
        [block],
        selected_group,
    )
    packed = _pack_selected_paths(selected_group, paths)
    contract = _eas_parallel_replay_block_contract(
        selected_group=selected_group,
        block=block,
        plan_completion_sha256=plan_completion_sha256,
        parallel_replay_plan=parallel_replay_plan,
        base_manifest=base_manifest,
        extension_plan=extension_plan,
    )
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        destination.rmdir()
    stage = Path(
        tempfile.mkdtemp(prefix=f".{destination.name}.staging.", dir=destination.parent)
    )
    try:
        _atomic_frame(stage / "selected_rows.tsv", selected_group)
        np.savez_compressed(stage / "paths.npz", **packed)
        outputs = {
            "selected_rows": _output_record(
                stage / "selected_rows.tsv", stage, rows=len(selected_group)
            ),
            "paths": _output_record(stage / "paths.npz", stage),
        }
        completion = {
            "schema": EAS_PARALLEL_REPLAY_BLOCK_SCHEMA,
            "status": "complete",
            "contract": contract,
            "contract_sha256": _canonical_sha256(contract),
            "outputs": outputs,
        }
        _atomic_json(stage / "completion.json", completion)
        os.replace(stage, destination)
    finally:
        if stage.exists():
            shutil.rmtree(stage, ignore_errors=True)
    result = verify_eas_parallel_replay_block(
        root,
        proposal_dir=proposal_dir,
        bundle_dir=destination,
        selected_group=selected_group,
        base_manifest=base_manifest,
        extension_plan=extension_plan,
        parallel_replay_plan=parallel_replay_plan,
        plan_completion_sha256=plan_completion_sha256,
    )
    return {key: value for key, value in result.items() if key != "paths"} | {
        "cache_hit": False,
        "worker_process_id": os.getpid(),
        "elapsed_seconds": time.time() - started,
    }


def _eas_parallel_replay_root_contract(
    *,
    selected: pd.DataFrame,
    block_bindings: Sequence[Mapping[str, Any]],
    plan_completion_sha256: str,
    parallel_replay_plan: Mapping[str, Any],
) -> dict[str, Any]:
    return {
        "schema": EAS_PARALLEL_REPLAY_COMPLETION_SCHEMA,
        "execution_method": EAS_PARALLEL_REPLAY_METHOD,
        "workflow_sha256": sha256_file(MODULE_PATH),
        "core_sha256": sha256_file(conditioned.MODULE_PATH),
        "plan_completion_sha256": str(plan_completion_sha256),
        "legacy_v2_plan_completion_sha256": LEGACY_V2_PLAN_COMPLETION_SHA256,
        "legacy_v3_plan_completion_sha256": LEGACY_V3_PLAN_COMPLETION_SHA256,
        "legacy_v3_replay_manifest_sha256": LEGACY_V3_REPLAY_MANIFEST_SHA256,
        "parallel_replay_plan_payload_sha256": parallel_replay_plan[
            "plan_payload_sha256"
        ],
        "selected_occurrences_sha256": _portable_frame_sha256(selected),
        "selected_count": len(selected),
        "selected_unit_ids_in_original_order": selected["unit_id"].astype(str).tolist(),
        "unique_selected_block_count": len(block_bindings),
        "selected_block_bindings": [dict(value) for value in block_bindings],
        "legacy_v3_migrated_block_count": sum(
            binding.get("artifact_provenance_kind") == "exact_legacy_v3_migration"
            for binding in block_bindings
        ),
        "current_v4_block_count": sum(
            binding.get("artifact_provenance_kind") == "current_v4"
            for binding in block_bindings
        ),
        "task_contract": "exactly_one_process_task_per_selected_unique_block",
        "adaptive_selection_or_replay": False,
    }


def _verify_any_eas_parallel_replay_block(
    repo_root: str | Path,
    *,
    proposal_dir: str | Path,
    bundle_dir: str | Path,
    selected_group: pd.DataFrame,
    base_manifest: Mapping[str, Any],
    extension_plan: Mapping[str, Any],
    parallel_replay_plan: Mapping[str, Any],
    legacy_v3_parallel_replay_plan: Mapping[str, Any],
    plan_completion_sha256: str,
) -> dict[str, Any]:
    bundle = Path(bundle_dir).resolve()
    completion = _read_json(bundle / "completion.json", "parallel replay block")
    block_index = int(selected_group["block_index"].iloc[0])
    if completion.get("schema") == LEGACY_V3_PARALLEL_REPLAY_BLOCK_SCHEMA:
        return verify_legacy_v3_parallel_replay_block(
            repo_root,
            proposal_dir=proposal_dir,
            bundle_dir=bundle,
            selected_group=selected_group,
            base_manifest=base_manifest,
            extension_plan=extension_plan,
            legacy_parallel_replay_plan=legacy_v3_parallel_replay_plan,
            expected_bundle_binding=_legacy_v3_replay_bundle_binding(
                parallel_replay_plan, block_index
            ),
        )
    return verify_eas_parallel_replay_block(
        repo_root,
        proposal_dir=proposal_dir,
        bundle_dir=bundle,
        selected_group=selected_group,
        base_manifest=base_manifest,
        extension_plan=extension_plan,
        parallel_replay_plan=parallel_replay_plan,
        plan_completion_sha256=plan_completion_sha256,
    )


def verify_eas_parallel_replay_completion(
    repo_root: str | Path,
    *,
    proposal_dir: str | Path,
    replay_dir: str | Path,
    selected: pd.DataFrame,
    base_manifest: Mapping[str, Any],
    extension_plan: Mapping[str, Any],
    parallel_replay_plan: Mapping[str, Any],
    legacy_v3_parallel_replay_plan: Mapping[str, Any],
    plan_completion_sha256: str,
    legacy_v4_parallel_replay_plan: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Verify every replay block and assemble paths in selected-row order."""

    root = Path(repo_root).resolve()
    replay_root = Path(replay_dir).resolve()
    if parallel_replay_plan.get("schema") == EAS_PARALLEL_REPLAY_PLAN_SCHEMA:
        if legacy_v4_parallel_replay_plan is None:
            raise ValueError("v5 replay verification requires the exact v4 replay plan")
        attestation = parallel_replay_plan.get("legacy_v4_plan_attestation")
        if (
            not isinstance(attestation, dict)
            or legacy_v4_parallel_replay_plan.get("plan_payload_sha256")
            != LEGACY_V4_REPLAY_PLAN_PAYLOAD_SHA256
            or attestation.get("legacy_plan_completion_sha256")
            != LEGACY_V4_PLAN_COMPLETION_SHA256
            or attestation.get("legacy_plan_producer_workflow_sha256")
            != LEGACY_V4_WORKFLOW_SHA256
            or attestation.get("legacy_replay_root_manifest")
            != _legacy_v4_replay_root_manifest(replay_root)
        ):
            raise ValueError("v5 replay root is not the exact attested v4 predecessor")
        _validate_eas_parallel_replay_plan(
            parallel_replay_plan,
            legacy_v4_attestation=attestation,
        )
        return _verify_exact_legacy_v4_replay_root_semantics(
            root,
            proposal_dir=proposal_dir,
            replay_dir=replay_root,
            selected=selected,
            base_manifest=base_manifest,
            extension_plan=extension_plan,
            legacy_v4_parallel_replay_plan=legacy_v4_parallel_replay_plan,
            legacy_v3_parallel_replay_plan=legacy_v3_parallel_replay_plan,
        )
    block_indexes = sorted(
        set(pd.to_numeric(selected["block_index"], errors="raise").astype(int))
    )
    expected_names = {
        EAS_PARALLEL_REPLAY_COMPLETION_FILENAME,
        *(f"block_{index:02d}" for index in block_indexes),
    }
    children = list(replay_root.iterdir()) if replay_root.is_dir() else []
    if (
        replay_root.is_symlink()
        or {path.name for path in children} != expected_names
        or any(
            path.is_symlink()
            or (
                path.name == EAS_PARALLEL_REPLAY_COMPLETION_FILENAME
                and not path.is_file()
            )
            or (
                path.name != EAS_PARALLEL_REPLAY_COMPLETION_FILENAME
                and not path.is_dir()
            )
            for path in children
        )
    ):
        raise ValueError("parallel EAS replay root inventory changed")
    paths_by_unit: dict[str, tuple[int, ...]] = {}
    block_bindings: list[dict[str, Any]] = []
    for block_index in block_indexes:
        selected_group = _selected_block_frame(selected, block_index)
        verified = _verify_any_eas_parallel_replay_block(
            root,
            proposal_dir=proposal_dir,
            bundle_dir=_eas_parallel_replay_block_dir(replay_root, block_index),
            selected_group=selected_group,
            base_manifest=base_manifest,
            extension_plan=extension_plan,
            parallel_replay_plan=parallel_replay_plan,
            legacy_v3_parallel_replay_plan=legacy_v3_parallel_replay_plan,
            plan_completion_sha256=plan_completion_sha256,
        )
        overlap = set(paths_by_unit).intersection(verified["paths"])
        if overlap:
            raise RuntimeError("parallel EAS replay unit appears in multiple blocks")
        paths_by_unit.update(verified["paths"])
        block_bindings.append(
            {
                "block_index": block_index,
                "selected_occurrence_count": verified["selected_occurrence_count"],
                "selected_unique_proposal_count": verified[
                    "selected_unique_proposal_count"
                ],
                "replay_completion_sha256": verified["completion_sha256"],
                "artifact_provenance_kind": verified["artifact_provenance_kind"],
                "source_replay_plan_completion_sha256": verified[
                    "source_replay_plan_completion_sha256"
                ],
                "original_legacy_completion_sha256": verified[
                    "original_legacy_completion_sha256"
                ],
            }
        )
    expected_contract = _eas_parallel_replay_root_contract(
        selected=selected,
        block_bindings=block_bindings,
        plan_completion_sha256=plan_completion_sha256,
        parallel_replay_plan=parallel_replay_plan,
    )
    completion = _read_json(
        replay_root / EAS_PARALLEL_REPLAY_COMPLETION_FILENAME,
        "parallel EAS replay completion",
    )
    if (
        completion.get("schema") != EAS_PARALLEL_REPLAY_COMPLETION_SCHEMA
        or completion.get("status") != "complete"
        or completion.get("contract") != expected_contract
        or completion.get("contract_sha256") != _canonical_sha256(expected_contract)
    ):
        raise ValueError("parallel EAS replay completion contract changed")
    expected_ids = selected["unit_id"].astype(str).tolist()
    if set(paths_by_unit) != set(expected_ids):
        raise RuntimeError("parallel EAS replay final path inventory is incomplete")
    ordered_paths = {unit_id: paths_by_unit[unit_id] for unit_id in expected_ids}
    return {
        "status": "complete",
        "cache_hit": True,
        "completion_sha256": sha256_file(
            replay_root / EAS_PARALLEL_REPLAY_COMPLETION_FILENAME
        ),
        "unique_selected_block_count": len(block_indexes),
        "block_bindings": block_bindings,
        "paths": ordered_paths,
    }


def _write_eas_parallel_replay_completion(
    repo_root: str | Path,
    *,
    proposal_dir: str | Path,
    replay_dir: str | Path,
    selected: pd.DataFrame,
    base_manifest: Mapping[str, Any],
    extension_plan: Mapping[str, Any],
    parallel_replay_plan: Mapping[str, Any],
    legacy_v3_parallel_replay_plan: Mapping[str, Any],
    plan_completion_sha256: str,
) -> dict[str, Any]:
    replay_root = Path(replay_dir).resolve()
    completion_path = replay_root / EAS_PARALLEL_REPLAY_COMPLETION_FILENAME
    if completion_path.is_file():
        return verify_eas_parallel_replay_completion(
            repo_root,
            proposal_dir=proposal_dir,
            replay_dir=replay_root,
            selected=selected,
            base_manifest=base_manifest,
            extension_plan=extension_plan,
            parallel_replay_plan=parallel_replay_plan,
            legacy_v3_parallel_replay_plan=legacy_v3_parallel_replay_plan,
            plan_completion_sha256=plan_completion_sha256,
        )
    block_indexes = sorted(
        set(pd.to_numeric(selected["block_index"], errors="raise").astype(int))
    )
    block_bindings: list[dict[str, Any]] = []
    for block_index in block_indexes:
        verified = _verify_any_eas_parallel_replay_block(
            repo_root,
            proposal_dir=proposal_dir,
            bundle_dir=_eas_parallel_replay_block_dir(replay_root, block_index),
            selected_group=_selected_block_frame(selected, block_index),
            base_manifest=base_manifest,
            extension_plan=extension_plan,
            parallel_replay_plan=parallel_replay_plan,
            legacy_v3_parallel_replay_plan=legacy_v3_parallel_replay_plan,
            plan_completion_sha256=plan_completion_sha256,
        )
        block_bindings.append(
            {
                "block_index": block_index,
                "selected_occurrence_count": verified["selected_occurrence_count"],
                "selected_unique_proposal_count": verified[
                    "selected_unique_proposal_count"
                ],
                "replay_completion_sha256": verified["completion_sha256"],
                "artifact_provenance_kind": verified["artifact_provenance_kind"],
                "source_replay_plan_completion_sha256": verified[
                    "source_replay_plan_completion_sha256"
                ],
                "original_legacy_completion_sha256": verified[
                    "original_legacy_completion_sha256"
                ],
            }
        )
    contract = _eas_parallel_replay_root_contract(
        selected=selected,
        block_bindings=block_bindings,
        plan_completion_sha256=plan_completion_sha256,
        parallel_replay_plan=parallel_replay_plan,
    )
    _atomic_json(
        completion_path,
        {
            "schema": EAS_PARALLEL_REPLAY_COMPLETION_SCHEMA,
            "status": "complete",
            "contract": contract,
            "contract_sha256": _canonical_sha256(contract),
        },
    )
    return verify_eas_parallel_replay_completion(
        repo_root,
        proposal_dir=proposal_dir,
        replay_dir=replay_root,
        selected=selected,
        base_manifest=base_manifest,
        extension_plan=extension_plan,
        parallel_replay_plan=parallel_replay_plan,
        legacy_v3_parallel_replay_plan=legacy_v3_parallel_replay_plan,
        plan_completion_sha256=plan_completion_sha256,
    )


def _run_eas_parallel_selected_replays_v4_legacy_implementation(
    repo_root: str | Path,
    *,
    plan_dir: str | Path,
    proposal_dir: str | Path,
    replay_dir: str | Path,
    selected: pd.DataFrame,
    blocks: Sequence[Mapping[str, Any]],
    base_manifest: Mapping[str, Any],
    extension_plan: Mapping[str, Any],
    max_workers: int = MAX_WORKERS,
) -> dict[str, Any]:
    """Replay selected source blocks in independent restartable processes."""

    root = Path(repo_root).resolve()
    plan_root = Path(plan_dir).resolve()
    replay_root = Path(replay_dir).resolve()
    workers = _positive_workers(max_workers)
    verify_plan(root, plan_dir=plan_root, proposal_dir=proposal_dir)
    parallel_replay_plan = _read_json(
        plan_root / PLAN_OUTPUTS["eas_parallel_replay_plan"],
        "EAS parallel replay plan",
    )
    legacy_v3_parallel_replay_plan = _read_json(
        plan_root.parent
        / LEGACY_V3_PLAN_DIRNAME
        / LEGACY_V3_PLAN_OUTPUTS["eas_parallel_replay_plan"],
        "legacy v3 EAS parallel replay plan",
    )
    legacy_v3_replay_root = (
        plan_root.parent / LEGACY_V3_PARALLEL_REPLAY_DIRNAME
    ).resolve()
    plan_completion_sha256 = sha256_file(plan_root / PLAN_COMPLETION_FILENAME)
    if (
        len(selected) != EXPECTED_NEUTRAL_PER_DEMOGRAPHY
        or selected["unit_id"].duplicated().any()
        or selected["accepted_index"].astype(int).tolist()
        != list(range(EXPECTED_NEUTRAL_PER_DEMOGRAPHY))
        or selected["unit_id"].astype(str).tolist()
        != [
            f"{conditioned.EAS_MODEL_ID}_r{index:03d}"
            for index in range(EXPECTED_NEUTRAL_PER_DEMOGRAPHY)
        ]
        or pd.to_numeric(selected["cosi2_seed"], errors="raise").astype(int).tolist()
        != [
            conditioned.EAS_SEED_START + index
            for index in range(EXPECTED_NEUTRAL_PER_DEMOGRAPHY)
        ]
    ):
        raise ValueError("parallel EAS replay selected occurrence ledger changed")
    block_indexes = [int(block["block_index"]) for block in blocks]
    if block_indexes != list(range(conditioned.EAS_PROPOSAL_BLOCK_COUNT)):
        raise ValueError("parallel EAS replay requires the exact fixed-200 block bank")
    grouped = [
        (int(block_index), group.reset_index(drop=True))
        for block_index, group in selected.groupby("block_index", sort=True)
    ]
    indexed_blocks = {int(block["block_index"]): block for block in blocks}
    for block_index, group in grouped:
        if block_index not in indexed_blocks:
            raise ValueError(
                f"selected EAS replay source block is absent: {block_index}"
            )
        _validate_selected_group_against_block(group, indexed_blocks[block_index])
    expected_block_names = {f"block_{index:02d}" for index, _group in grouped}
    replay_root.mkdir(parents=True, exist_ok=True)
    children = list(replay_root.iterdir())
    allowed_names = {
        *expected_block_names,
        EAS_PARALLEL_REPLAY_COMPLETION_FILENAME,
    }
    if any(
        path.name not in allowed_names
        or path.is_symlink()
        or (path.name == EAS_PARALLEL_REPLAY_COMPLETION_FILENAME and not path.is_file())
        or (path.name != EAS_PARALLEL_REPLAY_COMPLETION_FILENAME and not path.is_dir())
        for path in children
    ):
        raise ValueError("parallel EAS replay root contains unexpected state")
    if (replay_root / EAS_PARALLEL_REPLAY_COMPLETION_FILENAME).is_file():
        cached = verify_eas_parallel_replay_completion(
            root,
            proposal_dir=proposal_dir,
            replay_dir=replay_root,
            selected=selected,
            base_manifest=base_manifest,
            extension_plan=extension_plan,
            parallel_replay_plan=parallel_replay_plan,
            legacy_v3_parallel_replay_plan=legacy_v3_parallel_replay_plan,
            plan_completion_sha256=plan_completion_sha256,
        )
        return {**cached, "worker_results": pd.DataFrame()}
    started = time.time()
    print(
        json.dumps(
            {
                "event": "eas_parallel_replay_process_pool_start",
                "replay_dir": str(replay_root),
                "selected_occurrence_count": len(selected),
                "unique_selected_block_count": len(grouped),
                "max_worker_processes": min(workers, len(grouped)),
                "plan_completion_sha256": plan_completion_sha256,
                "execution_method": EAS_PARALLEL_REPLAY_METHOD,
                "legacy_v3_cache_block_count": len(LEGACY_V3_REPLAY_BLOCK_INDEXES),
                "legacy_v3_replay_manifest_sha256": (LEGACY_V3_REPLAY_MANIFEST_SHA256),
            },
            sort_keys=True,
        ),
        file=sys.stderr,
        flush=True,
    )
    rows: list[dict[str, Any]] = []
    failures: list[tuple[int, str]] = []
    with ProcessPoolExecutor(max_workers=min(workers, len(grouped))) as pool:
        futures = {
            pool.submit(
                _write_eas_parallel_replay_block,
                str(root),
                str(Path(proposal_dir).resolve()),
                str(_eas_parallel_replay_block_dir(replay_root, block_index)),
                [
                    _portable_record(record)
                    for record in group.to_dict(orient="records")
                ],
                dict(base_manifest),
                dict(extension_plan),
                dict(parallel_replay_plan),
                plan_completion_sha256,
                str(_eas_parallel_replay_block_dir(legacy_v3_replay_root, block_index)),
                dict(legacy_v3_parallel_replay_plan),
            ): block_index
            for block_index, group in grouped
        }
        for future in as_completed(futures):
            block_index = futures[future]
            try:
                row = future.result()
                rows.append(row)
                print(
                    json.dumps(
                        {"event": "eas_parallel_replay_block_complete", **row},
                        sort_keys=True,
                    ),
                    file=sys.stderr,
                    flush=True,
                )
            except Exception as error:  # noqa: BLE001
                failures.append((block_index, str(error)))
                print(
                    json.dumps(
                        {
                            "event": "eas_parallel_replay_block_failed",
                            "block_index": block_index,
                            "error": str(error),
                        },
                        sort_keys=True,
                    ),
                    file=sys.stderr,
                    flush=True,
                )
    print(
        json.dumps(
            {
                "event": "eas_parallel_replay_process_pool_finished",
                "unique_selected_block_count": len(grouped),
                "completed_block_count": len(rows),
                "failed_block_count": len(failures),
                "elapsed_seconds": time.time() - started,
            },
            sort_keys=True,
        ),
        file=sys.stderr,
        flush=True,
    )
    if failures:
        detail = "; ".join(f"block {key}: {value}" for key, value in failures)
        raise RuntimeError(f"{len(failures)} parallel EAS replays failed: {detail}")
    result = _write_eas_parallel_replay_completion(
        root,
        proposal_dir=proposal_dir,
        replay_dir=replay_root,
        selected=selected,
        base_manifest=base_manifest,
        extension_plan=extension_plan,
        parallel_replay_plan=parallel_replay_plan,
        legacy_v3_parallel_replay_plan=legacy_v3_parallel_replay_plan,
        plan_completion_sha256=plan_completion_sha256,
    )
    return {
        **result,
        "cache_hit": bool(rows) and all(bool(row["cache_hit"]) for row in rows),
        "worker_results": pd.DataFrame(rows)
        .sort_values("block_index")
        .reset_index(drop=True),
    }


def _migrate_exact_legacy_v4_replay_root(
    repo_root: str | Path,
    *,
    proposal_dir: str | Path,
    source_replay_dir: str | Path,
    destination_replay_dir: str | Path,
    selected: pd.DataFrame,
    base_manifest: Mapping[str, Any],
    extension_plan: Mapping[str, Any],
    current_replay_plan: Mapping[str, Any],
    legacy_v4_parallel_replay_plan: Mapping[str, Any],
    legacy_v3_parallel_replay_plan: Mapping[str, Any],
) -> dict[str, Any]:
    """Atomically copy and deep-verify the entire exact v4 replay root."""

    source = Path(source_replay_dir).resolve()
    destination = Path(destination_replay_dir).resolve()
    stale_stages = list(
        destination.parent.glob(f".{destination.name}.legacy-v4-staging.*")
    )
    if stale_stages:
        raise ValueError("v5 replay migration has stale staging state")
    source_manifest = _legacy_v4_replay_root_manifest(source)
    attested_manifest = current_replay_plan.get("legacy_v4_plan_attestation", {}).get(
        "legacy_replay_root_manifest"
    )
    if source_manifest != attested_manifest:
        raise ValueError("legacy v4 replay source differs from the v5 plan attestation")
    if destination.exists() and any(destination.iterdir()):
        destination_manifest = _legacy_v4_replay_root_manifest(destination)
        if destination_manifest != source_manifest:
            raise ValueError("v5 replay migration destination bytes changed")
        verified = _verify_exact_legacy_v4_replay_root_semantics(
            repo_root,
            proposal_dir=proposal_dir,
            replay_dir=destination,
            selected=selected,
            base_manifest=base_manifest,
            extension_plan=extension_plan,
            legacy_v4_parallel_replay_plan=legacy_v4_parallel_replay_plan,
            legacy_v3_parallel_replay_plan=legacy_v3_parallel_replay_plan,
        )
        return {**verified, "cache_hit": True, "migration_file_count": 244}
    if destination.exists() and not destination.is_dir():
        raise ValueError("v5 replay migration destination is not a directory")
    if destination.exists():
        destination.rmdir()
    _verify_exact_legacy_v4_replay_root_semantics(
        repo_root,
        proposal_dir=proposal_dir,
        replay_dir=source,
        selected=selected,
        base_manifest=base_manifest,
        extension_plan=extension_plan,
        legacy_v4_parallel_replay_plan=legacy_v4_parallel_replay_plan,
        legacy_v3_parallel_replay_plan=legacy_v3_parallel_replay_plan,
    )
    destination.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(
        tempfile.mkdtemp(
            prefix=f".{destination.name}.legacy-v4-staging.",
            dir=destination.parent,
        )
    )
    try:
        for binding in source_manifest["legacy_file_bindings"]:
            relative = Path(str(binding["path"]))
            source_file = source / relative
            target_file = stage / relative
            target_file.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source_file, target_file)
            if (
                target_file.stat().st_size != int(binding["size_bytes"])
                or sha256_file(target_file) != binding["sha256"]
            ):
                raise RuntimeError("v5 replay migration copy changed source bytes")
        staged_manifest = _legacy_v4_replay_root_manifest(stage)
        if staged_manifest != source_manifest:
            raise RuntimeError("v5 staged replay root differs from its v4 source")
        _verify_exact_legacy_v4_replay_root_semantics(
            repo_root,
            proposal_dir=proposal_dir,
            replay_dir=stage,
            selected=selected,
            base_manifest=base_manifest,
            extension_plan=extension_plan,
            legacy_v4_parallel_replay_plan=legacy_v4_parallel_replay_plan,
            legacy_v3_parallel_replay_plan=legacy_v3_parallel_replay_plan,
        )
        if _legacy_v4_replay_root_manifest(source) != source_manifest:
            raise RuntimeError("legacy v4 replay source changed during migration")
        os.replace(stage, destination)
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    verified = _verify_exact_legacy_v4_replay_root_semantics(
        repo_root,
        proposal_dir=proposal_dir,
        replay_dir=destination,
        selected=selected,
        base_manifest=base_manifest,
        extension_plan=extension_plan,
        legacy_v4_parallel_replay_plan=legacy_v4_parallel_replay_plan,
        legacy_v3_parallel_replay_plan=legacy_v3_parallel_replay_plan,
    )
    return {**verified, "cache_hit": False, "migration_file_count": 244}


def run_eas_parallel_selected_replays(
    repo_root: str | Path,
    *,
    plan_dir: str | Path,
    proposal_dir: str | Path,
    replay_dir: str | Path,
    selected: pd.DataFrame,
    blocks: Sequence[Mapping[str, Any]],
    base_manifest: Mapping[str, Any],
    extension_plan: Mapping[str, Any],
    max_workers: int = MAX_WORKERS,
) -> dict[str, Any]:
    """Reuse the exact completed v4 replay root; never replay a proposal block."""

    root = Path(repo_root).resolve()
    plan_root = Path(plan_dir).resolve()
    _positive_workers(max_workers)
    verify_plan(root, plan_dir=plan_root, proposal_dir=proposal_dir)
    if plan_root.name == PLAN_DIRNAME:
        expected_replay = plan_root.parent / LEGACY_V5_PARALLEL_REPLAY_DIRNAME
        if Path(replay_dir).resolve() != expected_replay.resolve():
            raise ValueError("v6 requires the exact immutable v5 EAS replay root")
        _verify_legacy_v5_eas_state(
            root,
            replay_dir=expected_replay,
            selection_dir=plan_root.parent / LEGACY_V5_SELECTION_DIRNAME,
        )
        raise ValueError(
            "v6 never reruns EAS replay; verify the immutable v5 selection instead"
        )
    replay_plan = _read_json(
        plan_root / PLAN_OUTPUTS["eas_parallel_replay_plan"],
        "v5 EAS replay migration plan",
    )
    legacy_v4_plan_root = plan_root.parent / LEGACY_V4_PLAN_DIRNAME
    legacy_v4_replay_plan = _read_json(
        legacy_v4_plan_root / LEGACY_V4_PLAN_OUTPUTS["eas_parallel_replay_plan"],
        "legacy v4 EAS parallel replay plan",
    )
    legacy_v3_replay_plan = _read_json(
        plan_root.parent
        / LEGACY_V3_PLAN_DIRNAME
        / LEGACY_V3_PLAN_OUTPUTS["eas_parallel_replay_plan"],
        "legacy v3 EAS parallel replay plan",
    )
    if (
        len(selected) != EXPECTED_NEUTRAL_PER_DEMOGRAPHY
        or selected["unit_id"].duplicated().any()
        or selected["accepted_index"].astype(int).tolist()
        != list(range(EXPECTED_NEUTRAL_PER_DEMOGRAPHY))
        or selected["unit_id"].astype(str).tolist()
        != [
            f"{conditioned.EAS_MODEL_ID}_r{index:03d}"
            for index in range(EXPECTED_NEUTRAL_PER_DEMOGRAPHY)
        ]
        or pd.to_numeric(selected["cosi2_seed"], errors="raise").astype(int).tolist()
        != [
            conditioned.EAS_SEED_START + index
            for index in range(EXPECTED_NEUTRAL_PER_DEMOGRAPHY)
        ]
    ):
        raise ValueError("v5 EAS replay selected occurrence ledger changed")
    block_indexes = [int(block["block_index"]) for block in blocks]
    if block_indexes != list(range(conditioned.EAS_PROPOSAL_BLOCK_COUNT)):
        raise ValueError("v5 EAS replay migration requires the fixed-200 block bank")
    indexed_blocks = {int(block["block_index"]): block for block in blocks}
    for block_index, group in selected.groupby("block_index", sort=True):
        _validate_selected_group_against_block(
            group.reset_index(drop=True), indexed_blocks[int(block_index)]
        )
    started = time.time()
    print(
        json.dumps(
            {
                "event": "eas_v5_exact_v4_replay_root_migration_start",
                "source_replay_dir": str(
                    plan_root.parent / LEGACY_V4_PARALLEL_REPLAY_DIRNAME
                ),
                "destination_replay_dir": str(Path(replay_dir).resolve()),
                "source_file_count": 244,
                "source_selected_block_count": 81,
                "legacy_v4_plan_completion_sha256": (LEGACY_V4_PLAN_COMPLETION_SHA256),
                "legacy_v4_replay_completion_sha256": (
                    LEGACY_V4_REPLAY_COMPLETION_SHA256
                ),
                "proposal_replay": False,
            },
            sort_keys=True,
        ),
        file=sys.stderr,
        flush=True,
    )
    result = _migrate_exact_legacy_v4_replay_root(
        root,
        proposal_dir=proposal_dir,
        source_replay_dir=(plan_root.parent / LEGACY_V4_PARALLEL_REPLAY_DIRNAME),
        destination_replay_dir=replay_dir,
        selected=selected,
        base_manifest=base_manifest,
        extension_plan=extension_plan,
        current_replay_plan=replay_plan,
        legacy_v4_parallel_replay_plan=legacy_v4_replay_plan,
        legacy_v3_parallel_replay_plan=legacy_v3_replay_plan,
    )
    print(
        json.dumps(
            {
                "event": "eas_v5_exact_v4_replay_root_migration_complete",
                "cache_hit": bool(result["cache_hit"]),
                "source_file_count": result["migration_file_count"],
                "selected_occurrence_count": len(selected),
                "unique_selected_block_count": result["unique_selected_block_count"],
                "elapsed_seconds": time.time() - started,
            },
            sort_keys=True,
        ),
        file=sys.stderr,
        flush=True,
    )
    return {**result, "worker_results": pd.DataFrame()}


def _selected_particle_row(
    selected_index: pd.DataFrame,
    unit_id: str,
) -> dict[str, Any]:
    """Return one complete core-selected row while preserving its indexed ID."""

    if "unit_id" not in selected_index.columns:
        raise ValueError("selected-particle index was created with drop=True")
    row = selected_index.loc[str(unit_id)]
    if not isinstance(row, pd.Series):
        raise ValueError(f"selected-particle unit ID is not unique: {unit_id}")
    record = _portable_record(row.to_dict())
    if record.get("unit_id") != str(unit_id):
        raise ValueError("selected-particle indexed unit ID changed")
    conditioned._eas_selected_particle_record(record)
    return record


def write_eas_global_selection(
    repo_root: str | Path,
    *,
    plan_dir: str | Path,
    proposal_dir: str | Path,
    selection_dir: str | Path,
    replay_dir: str | Path | None = None,
    max_workers: int = MAX_WORKERS,
) -> dict[str, Any]:
    """Select 100 occurrence units and persist only their replayed paths."""

    root = Path(repo_root).resolve()
    plan_root = Path(plan_dir).resolve()
    destination = Path(selection_dir).resolve()
    replay_root = (
        Path(replay_dir).resolve()
        if replay_dir is not None
        else destination.parent / EAS_PARALLEL_REPLAY_DIRNAME
    )
    if list(destination.parent.glob(f".{destination.name}.staging.*")):
        raise ValueError("EAS selection has stale staging state")
    verify_plan(root, plan_dir=plan_root, proposal_dir=proposal_dir)
    if plan_root.name == PLAN_DIRNAME:
        expected_selection = plan_root.parent / LEGACY_V5_SELECTION_DIRNAME
        expected_replay = plan_root.parent / LEGACY_V5_PARALLEL_REPLAY_DIRNAME
        if (
            destination != expected_selection.resolve()
            or replay_root != expected_replay.resolve()
        ):
            raise ValueError("v6 requires the exact immutable v5 EAS roots")
        if not (destination / EAS_SELECTION_COMPLETION).is_file():
            raise ValueError("v6 immutable v5 EAS selection is absent")
        return verify_eas_global_selection(
            root,
            plan_dir=plan_root,
            proposal_dir=proposal_dir,
            selection_dir=destination,
            replay_dir=replay_root,
        )
    base_manifest = _read_json(
        plan_root / PLAN_OUTPUTS["eas_base_bank_manifest"],
        "EAS base-bank manifest",
    )
    extension_plan = _read_json(
        plan_root / PLAN_OUTPUTS["eas_extension_plan"],
        "EAS extension plan",
    )
    parallel_replay_plan = _read_json(
        plan_root / PLAN_OUTPUTS["eas_parallel_replay_plan"],
        "EAS parallel replay plan",
    )
    blocks = load_eas_proposal_blocks(
        root,
        proposal_dir=proposal_dir,
        base_manifest=base_manifest,
        extension_plan=extension_plan,
    )
    selected, diagnostics = conditioned.select_eas_global_proposals(
        blocks,
        repo_root=root,
        base_manifest=base_manifest,
        extension_plan=extension_plan,
    )
    if destination.exists() and any(destination.iterdir()):
        return verify_eas_global_selection(
            root,
            plan_dir=plan_root,
            proposal_dir=proposal_dir,
            selection_dir=destination,
            replay_dir=replay_root,
        )
    if destination.exists() and not destination.is_dir():
        raise ValueError("EAS selection destination is not a directory")
    if destination.exists():
        destination.rmdir()
    replay = run_eas_parallel_selected_replays(
        root,
        plan_dir=plan_root,
        proposal_dir=proposal_dir,
        replay_dir=replay_root,
        selected=selected,
        blocks=blocks,
        base_manifest=base_manifest,
        extension_plan=extension_plan,
        max_workers=max_workers,
    )
    destination.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(
        tempfile.mkdtemp(prefix=f".{destination.name}.staging.", dir=destination.parent)
    )
    try:
        paths = replay["paths"]
        _atomic_frame(stage / "selected_occurrences.tsv", selected)
        _atomic_json(stage / "global_diagnostics.json", diagnostics)
        unit_records: dict[str, dict[str, Any]] = {}
        selected_index = selected.set_index("unit_id", drop=False)
        for unit_id, counts in paths.items():
            unit_dir = stage / "units" / unit_id
            row = _selected_particle_row(selected_index, unit_id)
            conditioned.write_eas_selected_path_bundle(
                root,
                unit_dir,
                selected_row=row,
                derived_counts_ascending=counts,
                global_diagnostics=diagnostics,
            )
            unit_records[unit_id] = {
                "completion_sha256": sha256_file(unit_dir / "completion.json"),
                "trajectory_sha256": sha256_file(unit_dir / "trajectory.tsv"),
            }
        contract = {
            "schema": EAS_BANK_WORKFLOW_SCHEMA,
            "proposal_bank_contract": (
                conditioned.EasProposalBankContract().to_record()
            ),
            "plan_completion_sha256": sha256_file(plan_root / PLAN_COMPLETION_FILENAME),
            "legacy_v2_plan_completion_sha256": (LEGACY_V2_PLAN_COMPLETION_SHA256),
            "legacy_v3_plan_completion_sha256": (LEGACY_V3_PLAN_COMPLETION_SHA256),
            "legacy_v4_plan_completion_sha256": (LEGACY_V4_PLAN_COMPLETION_SHA256),
            "legacy_v4_replay_completion_sha256": (LEGACY_V4_REPLAY_COMPLETION_SHA256),
            "legacy_v4_replay_root_manifest_sha256": (
                LEGACY_V4_REPLAY_ROOT_MANIFEST_SHA256
            ),
            "base_manifest_payload_sha256": base_manifest["manifest_payload_sha256"],
            "extension_plan_payload_sha256": extension_plan["plan_payload_sha256"],
            "parallel_replay_plan_payload_sha256": parallel_replay_plan[
                "plan_payload_sha256"
            ],
            "parallel_replay_execution_method": EAS_PARALLEL_REPLAY_METHOD,
            "parallel_replay_completion_interpretation": (
                "hard_bound_legacy_v4_bytes_not_current_v5_completion"
            ),
            "parallel_replay_completion_sha256": replay["completion_sha256"],
            "parallel_replay_unique_selected_block_count": replay[
                "unique_selected_block_count"
            ],
            "workflow_sha256": sha256_file(MODULE_PATH),
            "proposal_bank_binding_sha256": diagnostics["proposal_bank_binding_sha256"],
            "proposal_completion_sha256": {
                str(block["block_index"]): str(block["completion_sha256"])
                for block in blocks
            },
            "conditioned_null_sha256": sha256_file(conditioned.MODULE_PATH),
            "selection_method": ("fixed_seed_categorical_multinomial_with_replacement"),
            "resample_seed": conditioned.EAS_CATEGORICAL_RESAMPLE_SEED,
            "selected_count": EXPECTED_NEUTRAL_PER_DEMOGRAPHY,
            "selected_row_adapter": "set_index_drop_false_preserve_unit_id",
        }
        outputs = {
            "selected_occurrences": _output_record(
                stage / "selected_occurrences.tsv", stage, rows=len(selected)
            ),
            "global_diagnostics": _output_record(
                stage / "global_diagnostics.json", stage
            ),
        }
        completion = {
            "schema": EAS_SELECTION_SCHEMA,
            "status": "complete",
            "contract": contract,
            "contract_sha256": _canonical_sha256(contract),
            "outputs": outputs,
            "unit_outputs": unit_records,
        }
        _atomic_json(stage / EAS_SELECTION_COMPLETION, completion)
        os.replace(stage, destination)
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    return verify_eas_global_selection(
        root,
        plan_dir=plan_root,
        proposal_dir=proposal_dir,
        selection_dir=destination,
        replay_dir=replay_root,
    )


def verify_eas_global_selection(
    repo_root: str | Path,
    *,
    plan_dir: str | Path,
    proposal_dir: str | Path,
    selection_dir: str | Path,
    replay_dir: str | Path | None = None,
) -> dict[str, Any]:
    """Verify selected particles, global diagnostics, and all path checksums."""

    root = Path(repo_root).resolve()
    plan_root = Path(plan_dir).resolve()
    verify_plan(root, plan_dir=plan_root, proposal_dir=proposal_dir)
    bundle = Path(selection_dir).resolve()
    replay_root = (
        Path(replay_dir).resolve()
        if replay_dir is not None
        else bundle.parent / EAS_PARALLEL_REPLAY_DIRNAME
    )
    if plan_root.name == PLAN_DIRNAME:
        expected_selection = plan_root.parent / LEGACY_V5_SELECTION_DIRNAME
        expected_replay = plan_root.parent / LEGACY_V5_PARALLEL_REPLAY_DIRNAME
        if (
            bundle != expected_selection.resolve()
            or replay_root != expected_replay.resolve()
        ):
            raise ValueError("v6 requires the exact immutable v5 EAS roots")
        state = _verify_legacy_v5_eas_state(
            root,
            replay_dir=expected_replay,
            selection_dir=expected_selection,
        )
        return {
            **state["selection_completion"],
            "cache_hit": True,
            "selected": state["selected"],
            "diagnostics": state["diagnostics"],
            "legacy_v5_reuse": True,
        }
    base_manifest = _read_json(
        plan_root / PLAN_OUTPUTS["eas_base_bank_manifest"],
        "EAS base-bank manifest",
    )
    extension_plan = _read_json(
        plan_root / PLAN_OUTPUTS["eas_extension_plan"],
        "EAS extension plan",
    )
    parallel_replay_plan = _read_json(
        plan_root / PLAN_OUTPUTS["eas_parallel_replay_plan"],
        "EAS parallel replay plan",
    )
    legacy_v4_parallel_replay_plan = _read_json(
        plan_root.parent
        / LEGACY_V4_PLAN_DIRNAME
        / LEGACY_V4_PLAN_OUTPUTS["eas_parallel_replay_plan"],
        "legacy v4 EAS parallel replay plan",
    )
    legacy_v3_parallel_replay_plan = _read_json(
        plan_root.parent
        / LEGACY_V3_PLAN_DIRNAME
        / LEGACY_V3_PLAN_OUTPUTS["eas_parallel_replay_plan"],
        "legacy v3 EAS parallel replay plan",
    )
    if bundle.is_symlink() or not bundle.is_dir():
        raise ValueError(f"EAS global selection is absent: {bundle}")
    completion = _read_json(bundle / EAS_SELECTION_COMPLETION, "EAS selection")
    if {path.name for path in bundle.iterdir()} != {
        EAS_SELECTION_COMPLETION,
        "selected_occurrences.tsv",
        "global_diagnostics.json",
        "units",
    } or any(path.is_symlink() for path in bundle.iterdir()):
        raise ValueError("EAS global selection top-level inventory changed")
    contract = completion.get("contract")
    if (
        completion.get("schema") != EAS_SELECTION_SCHEMA
        or completion.get("status") != "complete"
        or not isinstance(contract, dict)
        or completion.get("contract_sha256") != _canonical_sha256(contract)
        or contract.get("conditioned_null_sha256")
        != sha256_file(conditioned.MODULE_PATH)
    ):
        raise ValueError("EAS global selection contract is invalid")
    blocks = load_eas_proposal_blocks(
        root,
        proposal_dir=proposal_dir,
        base_manifest=base_manifest,
        extension_plan=extension_plan,
    )
    selected, diagnostics = conditioned.select_eas_global_proposals(
        blocks,
        repo_root=root,
        base_manifest=base_manifest,
        extension_plan=extension_plan,
    )
    replay = verify_eas_parallel_replay_completion(
        root,
        proposal_dir=proposal_dir,
        replay_dir=replay_root,
        selected=selected,
        base_manifest=base_manifest,
        extension_plan=extension_plan,
        parallel_replay_plan=parallel_replay_plan,
        legacy_v3_parallel_replay_plan=legacy_v3_parallel_replay_plan,
        legacy_v4_parallel_replay_plan=legacy_v4_parallel_replay_plan,
        plan_completion_sha256=sha256_file(plan_root / PLAN_COMPLETION_FILENAME),
    )
    expected_contract = {
        "schema": EAS_BANK_WORKFLOW_SCHEMA,
        "proposal_bank_contract": conditioned.EasProposalBankContract().to_record(),
        "plan_completion_sha256": sha256_file(plan_root / PLAN_COMPLETION_FILENAME),
        "legacy_v2_plan_completion_sha256": LEGACY_V2_PLAN_COMPLETION_SHA256,
        "legacy_v3_plan_completion_sha256": LEGACY_V3_PLAN_COMPLETION_SHA256,
        "legacy_v4_plan_completion_sha256": LEGACY_V4_PLAN_COMPLETION_SHA256,
        "legacy_v4_replay_completion_sha256": LEGACY_V4_REPLAY_COMPLETION_SHA256,
        "legacy_v4_replay_root_manifest_sha256": (
            LEGACY_V4_REPLAY_ROOT_MANIFEST_SHA256
        ),
        "base_manifest_payload_sha256": base_manifest["manifest_payload_sha256"],
        "extension_plan_payload_sha256": extension_plan["plan_payload_sha256"],
        "parallel_replay_plan_payload_sha256": parallel_replay_plan[
            "plan_payload_sha256"
        ],
        "parallel_replay_execution_method": EAS_PARALLEL_REPLAY_METHOD,
        "parallel_replay_completion_interpretation": (
            "hard_bound_legacy_v4_bytes_not_current_v5_completion"
        ),
        "parallel_replay_completion_sha256": replay["completion_sha256"],
        "parallel_replay_unique_selected_block_count": replay[
            "unique_selected_block_count"
        ],
        "workflow_sha256": sha256_file(MODULE_PATH),
        "proposal_bank_binding_sha256": diagnostics["proposal_bank_binding_sha256"],
        "proposal_completion_sha256": {
            str(block["block_index"]): str(block["completion_sha256"])
            for block in blocks
        },
        "conditioned_null_sha256": sha256_file(conditioned.MODULE_PATH),
        "selection_method": "fixed_seed_categorical_multinomial_with_replacement",
        "resample_seed": conditioned.EAS_CATEGORICAL_RESAMPLE_SEED,
        "selected_count": EXPECTED_NEUTRAL_PER_DEMOGRAPHY,
        "selected_row_adapter": "set_index_drop_false_preserve_unit_id",
    }
    if contract != expected_contract:
        raise ValueError("EAS global selection no longer matches the fixed bank")
    observed = _read_eas_selected_frame(bundle / "selected_occurrences.tsv")
    pd.testing.assert_frame_equal(
        observed,
        selected,
        check_dtype=False,
        check_exact=True,
    )
    if _read_json(bundle / "global_diagnostics.json", "EAS diagnostics") != diagnostics:
        raise ValueError("EAS global diagnostics changed")
    outputs = completion.get("outputs", {})
    unit_outputs = completion.get("unit_outputs", {})
    if set(outputs) != {"selected_occurrences", "global_diagnostics"} or set(
        unit_outputs
    ) != set(selected["unit_id"].astype(str)):
        raise ValueError("EAS global selection output inventory changed")
    unit_root = bundle / "units"
    if (
        not unit_root.is_dir()
        or {path.name for path in unit_root.iterdir()} != set(unit_outputs)
        or any(path.is_symlink() or not path.is_dir() for path in unit_root.iterdir())
    ):
        raise ValueError("EAS selected path directory inventory changed")
    for label, filename in {
        "selected_occurrences": "selected_occurrences.tsv",
        "global_diagnostics": "global_diagnostics.json",
    }.items():
        record = outputs[label]
        path = bundle / filename
        if sha256_file(path) != record.get("sha256"):
            raise ValueError(f"EAS selection checksum failed: {label}")
    for row in selected.to_dict(orient="records"):
        unit_id = str(row["unit_id"])
        unit_dir = bundle / "units" / unit_id
        verified_path = conditioned.verify_eas_selected_path_bundle(
            root,
            unit_dir,
            expected_selected_row=row,
            expected_counts=replay["paths"][unit_id],
            expected_global_diagnostics=diagnostics,
        )
        if (
            sha256_file(unit_dir / "completion.json")
            != unit_outputs[unit_id]["completion_sha256"]
            or sha256_file(unit_dir / "trajectory.tsv")
            != unit_outputs[unit_id]["trajectory_sha256"]
        ):
            raise ValueError(f"EAS path completion changed: {unit_id}")
        if int(verified_path["birth_generation"]) != int(row["birth_generation"]):
            raise ValueError(f"EAS path age changed: {unit_id}")
    return {
        **completion,
        "cache_hit": True,
        "selected": selected,
        "diagnostics": diagnostics,
    }


HAN_ACCEPTED_FILENAME = "accepted_units.tsv"
HAN_SELECTION_COMPLETION = "selection_completion.json"


def _file_inventory(directory: Path, *, exclude: Sequence[str] = ()) -> dict[str, Any]:
    excluded = set(exclude)
    return {
        path.name: {
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        for path in sorted(directory.iterdir(), key=lambda value: value.name)
        if path.name not in excluded and path.is_file() and not path.is_symlink()
    }


def _han_candidate_dir(screen_root: Path, candidate_id: str) -> Path:
    return screen_root / "candidates" / candidate_id


def _expected_han_candidate_source_binding(
    candidate: Mapping[str, Any],
    *,
    plan_dir: Path,
    expected_plan_schema: str = PLAN_SCHEMA,
) -> dict[str, Any]:
    source_kind = str(candidate["source_kind"])
    if source_kind == "existing_complete_natural_neutral_candidate":
        source_dir = Path(str(candidate["source_unit_dir"])).resolve()
        completion_path = source_dir / natural_neutral.UNIT_COMPLETION_FILENAME
        source_completion = _read_json(
            completion_path, "existing natural-neutral completion"
        )
        trajectory = source_dir / "trajectory.tsv"
        source_outputs = source_completion.get("outputs", {})
        if (
            source_completion.get("status") != "complete"
            or source_completion.get("unit_id") != candidate["source_unit_id"]
            or not trajectory.is_file()
            or sha256_file(trajectory)
            != source_outputs.get("trajectory.tsv", {}).get("sha256")
        ):
            raise ValueError("existing Han neutral candidate failed source binding")
        return {
            "source_completion_path": str(completion_path.resolve()),
            "source_completion_sha256": sha256_file(completion_path),
            "source_trajectory_sha256": sha256_file(trajectory),
        }
    if source_kind != "new_screen_then_replay":
        raise ValueError(f"unknown Han candidate source kind: {source_kind}")
    plan_completion_path = plan_dir / PLAN_COMPLETION_FILENAME
    plan_completion = _read_json(plan_completion_path, "Han screen plan completion")
    outputs = plan_completion.get("outputs", {})
    config = plan_dir / PLAN_OUTPUTS["han_screen_config"]
    recombination_map = plan_dir / PLAN_OUTPUTS["recombination_map"]
    config_record = outputs.get("han_screen_config", {})
    map_record = outputs.get("recombination_map", {})
    if (
        plan_completion.get("schema") != expected_plan_schema
        or config_record.get("path") != PLAN_OUTPUTS["han_screen_config"]
        or config_record.get("sha256") != sha256_file(config)
        or map_record.get("path") != PLAN_OUTPUTS["recombination_map"]
        or map_record.get("sha256") != RECOMBINATION_MAP_SHA256
        or map_record.get("sha256") != sha256_file(recombination_map)
        or recombination_map.stat().st_size != RECOMBINATION_MAP_SIZE_BYTES
    ):
        raise ValueError("Han plan/config/recombination-map binding changed")
    return {
        "plan_completion_path": str(plan_completion_path.resolve()),
        "plan_completion_sha256": sha256_file(plan_completion_path),
        "config_path": PLAN_OUTPUTS["han_screen_config"],
        "config_sha256": sha256_file(config),
        "recombination_map_path": PLAN_OUTPUTS["recombination_map"],
        "recombination_map_size_bytes": RECOMBINATION_MAP_SIZE_BYTES,
        "recombination_map_sha256": RECOMBINATION_MAP_SHA256,
        "working_directory": str(plan_dir.resolve()),
    }


def verify_han_screen_candidate(
    candidate: Mapping[str, Any],
    *,
    plan_dir: str | Path,
    screen_dir: str | Path,
    cosi2_binary: str | Path,
    expected_plan_schema: str = PLAN_SCHEMA,
) -> dict[str, Any]:
    candidate = _portable_record(candidate)
    candidate_id = str(candidate["candidate_id"])
    plan_root = Path(plan_dir).resolve()
    unit_dir = _han_candidate_dir(Path(screen_dir).resolve(), candidate_id)
    completion = _read_json(unit_dir / "completion.json", "Han screen completion")
    source_kind = str(candidate["source_kind"])
    expected_source_binding = _expected_han_candidate_source_binding(
        candidate,
        plan_dir=plan_root,
        expected_plan_schema=expected_plan_schema,
    )
    if (
        completion.get("schema") != conditioned.HAN_SCREEN_SCHEMA
        or completion.get("status") != "complete"
        or completion.get("candidate_id") != candidate_id
        or int(completion.get("seed", -1)) != int(candidate["seed"])
        or completion.get("candidate_record_sha256") != _canonical_sha256(candidate)
        or completion.get("source_kind") != source_kind
        or completion.get("source_binding") != expected_source_binding
        or completion.get("cosi2_binary_sha256")
        != sha256_file(Path(cosi2_binary).resolve())
    ):
        raise ValueError(f"Han screen completion changed: {candidate_id}")
    actual = _file_inventory(unit_dir, exclude=("completion.json",))
    children = list(unit_dir.iterdir())
    if any(path.is_symlink() or not path.is_file() for path in children) or {
        path.name for path in children
    } != {"completion.json", *actual}:
        raise ValueError(f"Han screen file inventory changed: {candidate_id}")
    if completion.get("outputs") != actual:
        raise ValueError(f"Han screen output inventory changed: {candidate_id}")
    command = completion.get("command")
    if source_kind == "existing_complete_natural_neutral_candidate":
        if command is not None or set(actual) != {"trajectory.tsv"}:
            raise ValueError(f"inherited Han command/output changed: {candidate_id}")
    else:
        parameter_path = (
            Path(str(command[2]))
            if isinstance(command, list) and len(command) > 2
            else Path()
        )
        expected_tail = [
            "-n",
            "1",
            "-r",
            str(int(candidate["seed"])),
            "-u",
            "1",
            "-m",
        ]
        if (
            not isinstance(command, list)
            or len(command) != 10
            or command[0] != str(Path(cosi2_binary).resolve())
            or command[1] != "-p"
            or command[3:] != expected_tail
            or parameter_path.name != "model.params"
            or not parameter_path.parent.name.startswith(f".{candidate_id}.stage.")
            or parameter_path.parent.parent.resolve() != unit_dir.parent.resolve()
            or sha256_file(unit_dir / "model.params")
            != expected_source_binding["config_sha256"]
            or set(actual)
            != {"model.params", "screen.ms", "screen.stderr.log", "trajectory.tsv"}
        ):
            raise ValueError(f"generated Han command/output changed: {candidate_id}")
    validation = conditioned.validate_han_trajectory(
        unit_dir / "trajectory.tsv", require_gate=False
    )
    if completion.get("trajectory_validation") != validation:
        raise ValueError(f"Han screen validation changed: {candidate_id}")
    return {
        "candidate_id": candidate_id,
        "status": "complete",
        "trajectory_path": str((unit_dir / "trajectory.tsv").resolve()),
        "trajectory_sha256": actual["trajectory.tsv"]["sha256"],
        "present_chb_af": validation["present_chb_af"],
        "han_generation_645_chb_af": validation["han_generation_645_chb_af"],
        "gate_passed": validation["han_generation_645_gate_passed"],
        "strictly_segregating": validation["present_strictly_segregating"],
        "screen_completion_path": str((unit_dir / "completion.json").resolve()),
        "screen_completion_sha256": sha256_file(unit_dir / "completion.json"),
    }


def _run_han_screen_candidate(
    repo_root: Path,
    candidate: Mapping[str, Any],
    *,
    plan_dir: Path,
    screen_root: Path,
    binary: Path,
    max_attempts: int,
    timeout_seconds: float | None,
) -> dict[str, Any]:
    candidate = _portable_record(candidate)
    candidate_id = str(candidate["candidate_id"])
    final_dir = _han_candidate_dir(screen_root, candidate_id)
    if final_dir.exists():
        return {
            **verify_han_screen_candidate(
                candidate,
                plan_dir=plan_dir,
                screen_dir=screen_root,
                cosi2_binary=binary,
            ),
            "cache_hit": True,
        }
    final_dir.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(
        tempfile.mkdtemp(prefix=f".{candidate_id}.stage.", dir=final_dir.parent)
    )
    started = time.monotonic()
    try:
        source_kind = str(candidate["source_kind"])
        source_binding: dict[str, Any]
        command: list[str] | None
        if source_kind == "existing_complete_natural_neutral_candidate":
            source_dir = Path(str(candidate["source_unit_dir"])).resolve()
            source_completion = _read_json(
                source_dir / natural_neutral.UNIT_COMPLETION_FILENAME,
                "existing natural-neutral completion",
            )
            source_outputs = source_completion.get("outputs", {})
            trajectory_source = source_dir / "trajectory.tsv"
            if (
                source_completion.get("status") != "complete"
                or source_completion.get("unit_id") != candidate["source_unit_id"]
                or sha256_file(trajectory_source)
                != source_outputs.get("trajectory.tsv", {}).get("sha256")
            ):
                raise ValueError("existing Han neutral candidate failed source binding")
            shutil.copyfile(trajectory_source, stage / "trajectory.tsv")
            command = None
            source_binding = _expected_han_candidate_source_binding(
                candidate, plan_dir=plan_dir
            )
        else:
            config = plan_dir / PLAN_OUTPUTS["han_screen_config"]
            shutil.copyfile(config, stage / "model.params")
            command = [
                str(binary),
                "-p",
                str((stage / "model.params").resolve()),
                "-n",
                "1",
                "-r",
                str(int(candidate["seed"])),
                "-u",
                "1",
                "-m",
            ]
            environment = {
                key: value
                for key, value in os.environ.items()
                if not key.startswith("COSI_")
            }
            environment.update(
                {
                    "COSI_MAXATTEMPTS": str(max_attempts),
                    "COSI_SAVE_TRAJ": str((stage / "trajectory.tsv").resolve()),
                }
            )
            with (
                (stage / "screen.ms").open("wb") as stdout_handle,
                (stage / "screen.stderr.log").open("wb") as stderr_handle,
            ):
                process = subprocess.run(
                    command,
                    cwd=plan_dir,
                    env=environment,
                    stdout=stdout_handle,
                    stderr=stderr_handle,
                    check=False,
                    timeout=timeout_seconds,
                )
            if process.returncode != 0:
                raise RuntimeError(f"CoSi2 Han screen exited {process.returncode}")
            source_binding = _expected_han_candidate_source_binding(
                candidate, plan_dir=plan_dir
            )
        validation = conditioned.validate_han_trajectory(
            stage / "trajectory.tsv", require_gate=False
        )
        outputs = _file_inventory(stage)
        completion = {
            "schema": conditioned.HAN_SCREEN_SCHEMA,
            "status": "complete",
            "candidate_id": candidate_id,
            "candidate_order": int(candidate["candidate_order"]),
            "seed": int(candidate["seed"]),
            "candidate_record_sha256": _canonical_sha256(candidate),
            "source_kind": source_kind,
            "source_binding": source_binding,
            "cosi2_binary_path": str(binary),
            "cosi2_binary_sha256": sha256_file(binary),
            "command": command,
            "max_attempts": max_attempts,
            "elapsed_seconds": time.monotonic() - started,
            "trajectory_validation": validation,
            "outputs": outputs,
        }
        _atomic_json(stage / "completion.json", completion)
        os.replace(stage, final_dir)
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    return {
        **verify_han_screen_candidate(
            candidate,
            plan_dir=plan_dir,
            screen_dir=screen_root,
            cosi2_binary=binary,
        ),
        "cache_hit": False,
    }


def run_han_screens(
    repo_root: str | Path,
    *,
    plan_dir: str | Path,
    screen_dir: str | Path,
    cosi2_binary: str | Path,
    max_workers: int = MAX_WORKERS,
    max_attempts: int = natural_neutral.DEFAULT_COSI_MAXATTEMPTS,
    timeout_seconds: float | None = None,
    candidate_ids: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Run or reuse all requested candidates from the fixed Han ledger."""

    root = Path(repo_root).resolve()
    plan_root = Path(plan_dir).resolve()
    verify_plan(root, plan_dir=plan_root)
    ledger = pd.read_csv(plan_root / PLAN_OUTPUTS["han_candidates"], sep="\t")
    requested = (
        list(dict.fromkeys(str(value) for value in candidate_ids))
        if candidate_ids is not None
        else list(ledger["candidate_id"].astype(str))
    )
    unknown = sorted(set(requested).difference(ledger["candidate_id"].astype(str)))
    if unknown:
        raise ValueError("unknown Han candidate IDs: " + ", ".join(unknown))
    order = {candidate_id: index for index, candidate_id in enumerate(requested)}
    selected = ledger[ledger["candidate_id"].astype(str).isin(requested)].copy()
    selected["_order"] = selected["candidate_id"].map(order)
    selected = selected.sort_values("_order", kind="mergesort")
    binary = Path(cosi2_binary).resolve()
    screen_root = Path(screen_dir).resolve()
    workers = _positive_workers(max_workers)
    rows: list[dict[str, Any]] = []
    failures: list[tuple[str, str]] = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(
                _run_han_screen_candidate,
                root,
                record,
                plan_dir=plan_root,
                screen_root=screen_root,
                binary=binary,
                max_attempts=max_attempts,
                timeout_seconds=timeout_seconds,
            ): str(record["candidate_id"])
            for record in selected.drop(columns="_order").to_dict(orient="records")
        }
        for future in as_completed(futures):
            try:
                rows.append(future.result())
            except Exception as error:  # noqa: BLE001
                failures.append((futures[future], str(error)))
    if failures:
        detail = "; ".join(f"{key}: {value}" for key, value in failures[:10])
        raise RuntimeError(f"{len(failures)} Han screens failed: {detail}")
    result = pd.DataFrame(rows)
    result["_order"] = result["candidate_id"].map(order)
    return result.sort_values("_order").drop(columns="_order").reset_index(drop=True)


def freeze_han_acceptances(
    repo_root: str | Path,
    *,
    plan_dir: str | Path,
    screen_dir: str | Path,
    cosi2_binary: str | Path,
) -> pd.DataFrame:
    """Verify all screens and freeze the first 100 joint passes."""

    del repo_root
    plan_root = Path(plan_dir).resolve()
    screen_root = Path(screen_dir).resolve()
    ledger = pd.read_csv(plan_root / PLAN_OUTPUTS["han_candidates"], sep="\t")
    rows: list[dict[str, Any]] = []
    first_missing: int | None = None
    ledger_records = ledger.to_dict(orient="records")
    for index, row in enumerate(ledger_records):
        completion_path = (
            _han_candidate_dir(screen_root, str(row["candidate_id"]))
            / "completion.json"
        )
        if not completion_path.is_file():
            first_missing = index
            break
        rows.append(
            verify_han_screen_candidate(
                row,
                plan_dir=plan_root,
                screen_dir=screen_root,
                cosi2_binary=cosi2_binary,
            )
        )
    if first_missing is not None:
        later = ledger.iloc[first_missing + 1 :]["candidate_id"].map(
            lambda value: (
                _han_candidate_dir(screen_root, str(value)) / "completion.json"
            ).is_file()
        )
        if later.any():
            raise ValueError("Han completed screens do not form a contiguous prefix")
    results = pd.DataFrame(rows)
    candidate_root = screen_root / "candidates"
    expected_candidate_dirs = set(ledger.iloc[: len(rows)]["candidate_id"].astype(str))
    actual_candidate_children = (
        list(candidate_root.iterdir()) if candidate_root.is_dir() else []
    )
    if {
        path.name for path in actual_candidate_children
    } != expected_candidate_dirs or any(
        path.is_symlink() or not path.is_dir() for path in actual_candidate_children
    ):
        raise ValueError("Han candidate directory inventory is not the exact prefix")
    accepted = conditioned.select_han_acceptances(ledger, results)
    destination = screen_root / HAN_ACCEPTED_FILENAME
    completion_path = screen_root / HAN_SELECTION_COMPLETION
    contract = {
        "schema": conditioned.HAN_SELECTION_SCHEMA,
        "plan_completion_sha256": sha256_file(plan_root / PLAN_COMPLETION_FILENAME),
        "screen_completion_sha256": {
            str(row["candidate_id"]): str(row["screen_completion_sha256"])
            for row in results.to_dict(orient="records")
        },
        "selection_rule": (
            "first_100_candidate_order_joint_g645_gate_and_present_segregation"
        ),
    }
    if completion_path.is_file():
        completion = _read_json(completion_path, "Han selection completion")
        observed = pd.read_csv(destination, sep="\t")
        pd.testing.assert_frame_equal(
            observed, accepted, check_dtype=False, check_exact=False, rtol=0, atol=1e-15
        )
        if completion.get("contract_sha256") != _canonical_sha256(
            contract
        ) or completion.get("accepted_sha256") != sha256_file(destination):
            raise ValueError("Han accepted selection is stale")
        if {path.name for path in screen_root.iterdir()} != {
            "candidates",
            HAN_ACCEPTED_FILENAME,
            HAN_SELECTION_COMPLETION,
        } or any(path.is_symlink() for path in screen_root.iterdir()):
            raise ValueError("Han screen root inventory changed")
        return accepted
    _atomic_frame(destination, accepted)
    _atomic_json(
        completion_path,
        {
            "schema": conditioned.HAN_SELECTION_SCHEMA,
            "status": "complete",
            "contract": contract,
            "contract_sha256": _canonical_sha256(contract),
            "accepted_sha256": sha256_file(destination),
            "accepted_rows": len(accepted),
        },
    )
    if {path.name for path in screen_root.iterdir()} != {
        "candidates",
        HAN_ACCEPTED_FILENAME,
        HAN_SELECTION_COMPLETION,
    }:
        raise ValueError("Han screen root inventory changed after selection")
    return accepted


def _verify_legacy_v6_han_state(
    *,
    legacy_plan_dir: Path,
    screen_dir: Path,
    cosi2_binary: Path,
) -> dict[str, Any]:
    """Deep-verify the exact completed Han screen and v6-bound selection."""

    plan_root = legacy_plan_dir.resolve()
    screen_root = screen_dir.resolve()
    completion_path = screen_root / HAN_SELECTION_COMPLETION
    accepted_path = screen_root / HAN_ACCEPTED_FILENAME
    if (
        sha256_file(completion_path) != LEGACY_V6_HAN_SELECTION_COMPLETION_SHA256
        or sha256_file(accepted_path) != LEGACY_V6_HAN_ACCEPTED_SHA256
    ):
        raise ValueError("legacy v6 Han selection bytes changed")
    completion = _read_json(completion_path, "legacy v6 Han selection completion")
    contract = completion.get("contract")
    screen_bindings = (
        contract.get("screen_completion_sha256", {})
        if isinstance(contract, dict)
        else {}
    )
    if (
        completion.get("schema") != conditioned.HAN_SELECTION_SCHEMA
        or completion.get("status") != "complete"
        or completion.get("accepted_rows") != EXPECTED_NEUTRAL_PER_DEMOGRAPHY
        or completion.get("accepted_sha256") != LEGACY_V6_HAN_ACCEPTED_SHA256
        or completion.get("contract_sha256") != LEGACY_V6_HAN_SELECTION_CONTRACT_SHA256
        or completion.get("contract_sha256") != _canonical_sha256(contract)
        or contract.get("plan_completion_sha256") != LEGACY_V6_PLAN_COMPLETION_SHA256
        or not isinstance(screen_bindings, dict)
    ):
        raise ValueError("legacy v6 Han selection contract changed")
    ledger = pd.read_csv(plan_root / LEGACY_V6_PLAN_OUTPUTS["han_candidates"], sep="\t")
    screened_count = len(screen_bindings)
    if screened_count < EXPECTED_NEUTRAL_PER_DEMOGRAPHY or screened_count > len(ledger):
        raise ValueError("legacy v6 Han screened prefix size changed")
    prefix = ledger.iloc[:screened_count].copy()
    expected_ids = list(prefix["candidate_id"].astype(str))
    if list(screen_bindings) != expected_ids:
        raise ValueError("legacy v6 Han completion ledger is not the exact prefix")
    candidate_root = screen_root / "candidates"
    children = list(candidate_root.iterdir()) if candidate_root.is_dir() else []
    if (
        {path.name for path in children} != set(expected_ids)
        or any(path.is_symlink() or not path.is_dir() for path in children)
        or {path.name for path in screen_root.iterdir()}
        != {"candidates", HAN_ACCEPTED_FILENAME, HAN_SELECTION_COMPLETION}
        or any(path.is_symlink() for path in screen_root.iterdir())
    ):
        raise ValueError("legacy v6 Han screen inventory changed")
    observed_accepted = pd.read_csv(accepted_path, sep="\t")
    accepted_ids = set(observed_accepted["candidate_id"].astype(str))
    if len(accepted_ids) != EXPECTED_NEUTRAL_PER_DEMOGRAPHY:
        raise ValueError("legacy v6 Han accepted candidate IDs changed")
    if not accepted_ids.issubset(set(expected_ids)):
        raise ValueError(
            "legacy v6 Han accepted candidates are outside screened prefix"
        )
    if (
        list(pd.to_numeric(observed_accepted["accepted_index"]).astype(int))
        != list(range(EXPECTED_NEUTRAL_PER_DEMOGRAPHY))
        or not pd.to_numeric(
            observed_accepted["candidate_order"]
        ).is_monotonic_increasing
    ):
        raise ValueError("legacy v6 Han accepted order changed")
    completion_manifest: dict[str, str] = {
        str(key): str(value) for key, value in screen_bindings.items()
    }
    trajectory_manifest: dict[str, str] = {}
    ledger_by_id = ledger.set_index("candidate_id", drop=False)
    for accepted_row in observed_accepted.to_dict(orient="records"):
        candidate_id = str(accepted_row["candidate_id"])
        candidate_row = ledger_by_id.loc[candidate_id].to_dict()
        verified = verify_han_screen_candidate(
            candidate_row,
            plan_dir=plan_root,
            screen_dir=screen_root,
            cosi2_binary=cosi2_binary,
            expected_plan_schema=LEGACY_V6_PLAN_SCHEMA,
        )
        if (
            verified["screen_completion_sha256"] != screen_bindings[candidate_id]
            or verified["trajectory_sha256"] != accepted_row["trajectory_sha256"]
            or Path(verified["trajectory_path"]).resolve()
            != Path(str(accepted_row["trajectory_path"])).resolve()
            or not np.isclose(
                float(verified["present_chb_af"]),
                float(accepted_row["present_chb_af"]),
                rtol=0,
                atol=1e-15,
            )
            or not np.isclose(
                float(verified["han_generation_645_chb_af"]),
                float(accepted_row["han_generation_645_chb_af"]),
                rtol=0,
                atol=1e-15,
            )
            or not bool(verified["gate_passed"])
            or not bool(verified["strictly_segregating"])
        ):
            raise ValueError(f"legacy v6 accepted Han path changed: {candidate_id}")
        trajectory_manifest[candidate_id] = str(verified["trajectory_sha256"])
    return {
        "schema": "gamma-smc.cosi2-conditioned-gamma-v6-han-state/v7",
        "status": "exact_complete_v6_han_selection",
        "screened_candidate_count": screened_count,
        "accepted_rows": len(observed_accepted),
        "selection_completion_sha256": LEGACY_V6_HAN_SELECTION_COMPLETION_SHA256,
        "selection_contract_sha256": LEGACY_V6_HAN_SELECTION_CONTRACT_SHA256,
        "accepted_units_sha256": LEGACY_V6_HAN_ACCEPTED_SHA256,
        "candidate_completion_manifest_sha256": _canonical_sha256(completion_manifest),
        "accepted_candidate_trajectory_manifest_sha256": _canonical_sha256(
            trajectory_manifest
        ),
        "accepted": observed_accepted,
    }


def _legacy_v6_eas_ledger(
    predecessor: Mapping[str, Any],
    *,
    selection_dir: Path,
    han_accepted: pd.DataFrame,
) -> pd.DataFrame:
    eas = predecessor["selected"].copy()
    eas["demography"] = "EAS"
    eas["model_id"] = conditioned.EAS_MODEL_ID
    eas["seed"] = eas["cosi2_seed"].astype(int)
    eas["trajectory_path"] = eas["unit_id"].map(
        lambda value: str(
            (selection_dir / "units" / str(value) / "trajectory.tsv").resolve()
        )
    )
    eas["trajectory_sha256"] = eas["trajectory_path"].map(sha256_file)
    han = han_accepted.copy()
    han["demography"] = "Han"
    han["model_id"] = conditioned.HAN_MODEL_ID
    combined = pd.concat([eas, han], ignore_index=True, sort=False)
    combined = combined.sort_values(["demography", "unit_id"]).reset_index(drop=True)
    return combined[combined["demography"].astype(str) == "EAS"].copy()


def _verify_legacy_v6_eas_simulations(
    repo_root: Path,
    *,
    predecessor: Mapping[str, Any],
    legacy_plan_dir: Path,
    simulation_dir: Path,
    cosi2_binary: Path,
    han_accepted: pd.DataFrame,
) -> dict[str, Any]:
    """Verify and bind the original EAS-only v6 simulation root in place."""

    run_root = simulation_dir.resolve()
    if (
        run_root.is_symlink()
        or not run_root.is_dir()
        or {path.name for path in run_root.iterdir()} != {"units"}
        or any(path.is_symlink() for path in run_root.iterdir())
    ):
        raise ValueError("legacy v6 simulation root is not exact EAS-only inventory")
    selection_root = legacy_plan_dir.parent / LEGACY_V5_SELECTION_DIRNAME
    ledger = _legacy_v6_eas_ledger(
        predecessor,
        selection_dir=selection_root,
        han_accepted=han_accepted,
    )
    unit_root = run_root / "units"
    children = list(unit_root.iterdir()) if unit_root.is_dir() else []
    expected_ids = set(ledger["unit_id"].astype(str))
    if (
        len(expected_ids) != EXPECTED_NEUTRAL_PER_DEMOGRAPHY
        or {path.name for path in children} != expected_ids
        or any(path.is_symlink() or not path.is_dir() for path in children)
    ):
        raise ValueError("legacy v6 EAS unit directory inventory changed")
    for row in ledger.to_dict(orient="records"):
        _verify_legacy_v6_conditioned_unit(
            row,
            repo_root=repo_root,
            plan_dir=legacy_plan_dir,
            run_dir=run_root,
            cosi2_binary=cosi2_binary,
        )
    tree = _file_tree_manifest(unit_root)
    if (
        tree["file_count"] != LEGACY_V6_EAS_SIMULATION_FILE_COUNT
        or tree["manifest_sha256"] != LEGACY_V6_EAS_SIMULATION_ROOT_MANIFEST_SHA256
    ):
        raise ValueError("legacy v6 EAS simulation bytes changed")
    return {
        "schema": "gamma-smc.cosi2-conditioned-gamma-v6-eas-simulations/v7",
        "status": "exact_complete_v6_eas_only_root",
        "legacy_simulation_dirname": LEGACY_V6_SIMULATION_DIRNAME,
        "legacy_plan_completion_sha256": LEGACY_V6_PLAN_COMPLETION_SHA256,
        "unit_count": len(ledger),
        "file_count": tree["file_count"],
        "root_manifest_sha256": tree["manifest_sha256"],
    }


def _verify_legacy_v6_execution_state(
    repo_root: Path,
    *,
    predecessor: Mapping[str, Any],
    legacy_plan_dir: Path,
    screen_dir: Path,
    simulation_dir: Path,
    cosi2_binary: Path,
) -> dict[str, Any]:
    han = _verify_legacy_v6_han_state(
        legacy_plan_dir=legacy_plan_dir,
        screen_dir=screen_dir,
        cosi2_binary=cosi2_binary,
    )
    eas = _verify_legacy_v6_eas_simulations(
        repo_root,
        predecessor=predecessor,
        legacy_plan_dir=legacy_plan_dir,
        simulation_dir=simulation_dir,
        cosi2_binary=cosi2_binary,
        han_accepted=han["accepted"],
    )
    return {
        "schema": "gamma-smc.cosi2-conditioned-gamma-v6-execution-state/v7",
        "status": "exact_pre_v7_boundary",
        "han": {key: value for key, value in han.items() if key != "accepted"},
        "eas": eas,
        "v6_han_output_or_staging_present": False,
        "v6_eas_units_rewritten_or_copied_in_v7": False,
    }


def build_conditioned_simulation_ledger(
    repo_root: str | Path,
    *,
    plan_dir: str | Path,
    proposal_dir: str | Path,
    selection_dir: str | Path,
    screen_dir: str | Path,
    cosi2_binary: str | Path,
    verified_eas_bundle: Mapping[str, Any] | None = None,
) -> pd.DataFrame:
    """Bind the exact 100 EAS occurrences and 100 accepted Han candidates."""

    root = Path(repo_root).resolve()
    plan_root = Path(plan_dir).resolve()
    eas_bundle = (
        dict(verified_eas_bundle)
        if verified_eas_bundle is not None
        else verify_eas_global_selection(
            root,
            plan_dir=plan_root,
            proposal_dir=proposal_dir,
            selection_dir=selection_dir,
        )
    )
    if "selected" not in eas_bundle or "diagnostics" not in eas_bundle:
        raise ValueError("verified EAS selection bundle is incomplete")
    eas = eas_bundle["selected"].copy()
    eas["demography"] = "EAS"
    eas["model_id"] = conditioned.EAS_MODEL_ID
    eas["seed"] = eas["cosi2_seed"].astype(int)
    eas["trajectory_path"] = eas["unit_id"].map(
        lambda value: str(
            (Path(selection_dir).resolve() / "units" / str(value) / "trajectory.tsv")
        )
    )
    eas["trajectory_sha256"] = eas["trajectory_path"].map(sha256_file)

    if plan_root.name == PLAN_DIRNAME:
        han = _verify_legacy_v6_han_state(
            legacy_plan_dir=plan_root.parent / LEGACY_V6_PLAN_DIRNAME,
            screen_dir=Path(screen_dir).resolve(),
            cosi2_binary=Path(cosi2_binary).resolve(),
        )["accepted"].copy()
    else:
        han = freeze_han_acceptances(
            root,
            plan_dir=plan_root,
            screen_dir=screen_dir,
            cosi2_binary=cosi2_binary,
        ).copy()
    han["demography"] = "Han"
    han["model_id"] = conditioned.HAN_MODEL_ID
    required = {
        "unit_id",
        "demography",
        "model_id",
        "seed",
        "trajectory_path",
        "trajectory_sha256",
    }
    if not required.issubset(eas.columns) or not required.issubset(han.columns):
        raise ValueError("conditioned simulation ledger lacks trajectory provenance")
    ledger = pd.concat([eas, han], ignore_index=True, sort=False)
    if len(ledger) != 200 or ledger["unit_id"].duplicated().any():
        raise ValueError("conditioned simulation ledger must contain 200 unique units")
    return ledger.sort_values(["demography", "unit_id"]).reset_index(drop=True)


HAN_SOURCE_TRAJECTORY_COLUMNS = (
    "sim",
    "gen",
    "selfreq_1",
    "selfreq_2",
    "selfreq_3",
    "selfreq_4",
    "selfreq_5",
    "popsize_1",
    "popsize_2",
    "popsize_3",
    "popsize_4",
    "popsize_5",
)
HAN_LOADER_TRAJECTORY_COLUMNS = HAN_SOURCE_TRAJECTORY_COLUMNS[:7]


def _han_loader_trajectory_projection(
    source_trajectory: str | Path,
) -> tuple[bytes, dict[str, Any]]:
    """Lexically project a validated 12-column Han path to CoSi's 7-column ABI."""

    source = Path(source_trajectory).resolve()
    raw = source.read_bytes()
    if not raw.endswith(b"\n") or b"\r" in raw or b"\x00" in raw:
        raise ValueError("Han source trajectory must be canonical UTF-8/LF text")
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValueError("Han source trajectory is not UTF-8") from error
    lines = text[:-1].split("\n")
    if not lines or tuple(lines[0].split("\t")) != HAN_SOURCE_TRAJECTORY_COLUMNS:
        raise ValueError("Han source trajectory header is not the exact 12-column ABI")
    expected_rows = conditioned.HAN_BIRTH_GENERATION + 1
    if len(lines) != expected_rows + 1:
        raise ValueError("Han source trajectory row count changed")
    projected = ["\t".join(HAN_LOADER_TRAJECTORY_COLUMNS)]
    for row_index, line in enumerate(lines[1:]):
        fields = line.split("\t")
        if len(fields) != len(HAN_SOURCE_TRAJECTORY_COLUMNS):
            raise ValueError(f"Han source trajectory row {row_index + 2} is malformed")
        if (
            fields[0] != "1"
            or int(fields[1]) != conditioned.HAN_BIRTH_GENERATION - row_index
        ):
            raise ValueError("Han source trajectory sim/generation sequence changed")
        numeric = np.asarray([float(token) for token in fields[2:]], dtype=float)
        if (
            not np.isfinite(numeric).all()
            or (numeric[:5] < 0).any()
            or (numeric[:5] > 1).any()
            or (numeric[5:] <= 0).any()
        ):
            raise ValueError(
                "Han source trajectory contains an invalid frequency/popsize"
            )
        projected.append("\t".join(fields[:7]))
    loader_bytes = ("\n".join(projected) + "\n").encode("utf-8")
    metadata = {
        "schema": HAN_LOADER_TRAJECTORY_SCHEMA,
        "projection_method": (
            "lexical_tab_projection_fields_0_through_6_no_float_reformat_canonical_lf"
        ),
        "source_columns": list(HAN_SOURCE_TRAJECTORY_COLUMNS),
        "loader_columns": list(HAN_LOADER_TRAJECTORY_COLUMNS),
        "data_rows": expected_rows,
        "source_size_bytes": len(raw),
        "source_sha256": hashlib.sha256(raw).hexdigest(),
        "loader_size_bytes": len(loader_bytes),
        "loader_sha256": hashlib.sha256(loader_bytes).hexdigest(),
    }
    return loader_bytes, metadata


def _write_han_loader_trajectory(
    source_trajectory: str | Path, destination: str | Path
) -> dict[str, Any]:
    loader_bytes, metadata = _han_loader_trajectory_projection(source_trajectory)
    target = Path(destination)
    _atomic_text(target, loader_bytes.decode("utf-8"))
    if (
        target.read_bytes() != loader_bytes
        or sha256_file(target) != metadata["loader_sha256"]
    ):
        raise RuntimeError("Han loader trajectory changed during atomic publication")
    return metadata


def _validate_han_loader_trajectory(
    source_trajectory: str | Path, loader_trajectory: str | Path
) -> dict[str, Any]:
    expected, metadata = _han_loader_trajectory_projection(source_trajectory)
    loader = Path(loader_trajectory)
    if (
        loader.read_bytes() != expected
        or sha256_file(loader) != metadata["loader_sha256"]
    ):
        raise ValueError(
            "Han loader trajectory is not the exact lexical source projection"
        )
    return metadata


def _verify_legacy_v6_conditioned_unit(
    unit: Mapping[str, Any],
    *,
    repo_root: str | Path,
    plan_dir: str | Path,
    run_dir: str | Path,
    cosi2_binary: str | Path,
) -> dict[str, Any]:
    root = Path(repo_root).resolve()
    plan_root = Path(plan_dir).resolve()
    unit_id = str(unit["unit_id"])
    unit_dir = Path(run_dir).resolve() / "units" / unit_id
    completion = _read_json(unit_dir / "completion.json", "conditioned unit completion")
    plan_map = plan_root / PLAN_OUTPUTS["recombination_map"]
    if (
        completion.get("schema") != conditioned.UNIT_SCHEMA
        or completion.get("status") != "complete"
        or completion.get("unit_id") != unit_id
        or completion.get("ledger_row_sha256")
        != _canonical_sha256(_portable_record(unit))
        or completion.get("cosi2_binary_sha256")
        != sha256_file(Path(cosi2_binary).resolve())
        or completion.get("plan_completion_sha256")
        != sha256_file(plan_root / PLAN_COMPLETION_FILENAME)
        or completion.get("source_trajectory_path")
        != str(Path(str(unit["trajectory_path"])).resolve())
        or completion.get("source_trajectory_sha256") != str(unit["trajectory_sha256"])
        or completion.get("working_directory") != str(plan_root)
        or completion.get("plan_recombination_map_sha256") != RECOMBINATION_MAP_SHA256
        or completion.get("plan_recombination_map_sha256") != sha256_file(plan_map)
    ):
        raise ValueError(f"conditioned unit completion changed: {unit_id}")
    outputs = _file_inventory(unit_dir, exclude=("completion.json",))
    children = list(unit_dir.iterdir())
    if any(path.is_symlink() or not path.is_file() for path in children) or {
        path.name for path in children
    } != {"completion.json", *outputs}:
        raise ValueError(f"conditioned unit contains an unexpected file: {unit_id}")
    if completion.get("outputs") != outputs:
        raise ValueError(f"conditioned unit output inventory changed: {unit_id}")
    expected_parameters = (
        conditioned.render_eas_replay_parameter_text(
            root, birth_generation=int(unit["birth_generation"])
        )
        if str(unit["demography"]) == "EAS"
        else (plan_root / PLAN_OUTPUTS["han_replay_config"]).read_text(encoding="utf-8")
    )
    if (unit_dir / "model.params").read_text(encoding="utf-8") != expected_parameters:
        raise ValueError(f"conditioned unit parameters changed: {unit_id}")
    command = completion.get("command")
    parameter_path = (
        Path(str(command[2]))
        if isinstance(command, list) and len(command) > 2
        else Path()
    )
    expected_tail = [
        "-n",
        "1",
        "-r",
        str(int(unit["seed"])),
        "-u",
        "1",
        "-m",
    ]
    if (
        "recomb_file recombination_map.tsv" not in expected_parameters
        or not isinstance(command, list)
        or len(command) != 10
        or command[0] != str(Path(cosi2_binary).resolve())
        or command[1] != "-p"
        or command[3:] != expected_tail
        or parameter_path.name != "model.params"
        or not parameter_path.parent.name.startswith(f".{unit_id}.stage.")
        or parameter_path.parent.parent.resolve() != unit_dir.parent.resolve()
    ):
        raise ValueError(f"conditioned unit command/cwd/map changed: {unit_id}")
    environment = completion.get("environment_contract")
    if not isinstance(environment, dict):
        raise ValueError(f"conditioned unit environment is absent: {unit_id}")
    max_attempts = environment.get("COSI_MAXATTEMPTS")
    expected_environment = {
        "inherited_cosi_variables_removed": True,
        "COSI_MAXATTEMPTS": max_attempts,
        "COSI_LOAD_TRAJ": "checksum_bound_source_trajectory",
        "trajectory_copy": "checksum_bound_source_to_unit_staging",
    }
    if (
        isinstance(max_attempts, bool)
        or not isinstance(max_attempts, int)
        or max_attempts <= 0
        or environment != expected_environment
    ):
        raise ValueError(f"conditioned unit environment contract changed: {unit_id}")
    if str(unit["demography"]) == "EAS":
        source = pd.read_csv(str(unit["trajectory_path"]), sep="\t")
        observed = pd.read_csv(unit_dir / "trajectory.tsv", sep="\t")
        pd.testing.assert_frame_equal(
            observed, source, check_dtype=False, check_exact=False, rtol=0, atol=5e-16
        )
        trajectory_validation = {
            "status": "valid",
            "birth_generation": int(unit["birth_generation"]),
            "present_population_af": conditioned.EAS_SELECTED_PRESENT_AF,
        }
    else:
        trajectory_validation = conditioned.validate_han_trajectory(
            unit_dir / "trajectory.tsv", require_gate=True
        )
    if completion.get("trajectory_validation") != trajectory_validation:
        raise ValueError(f"conditioned trajectory validation changed: {unit_id}")
    ms_validation = natural_neutral.validate_cosi_ms(
        unit_dir / "simulation.ms",
        unit={"sample_haploids": analysis.PANEL_HAPLOTYPES, "seed": int(unit["seed"])},
    )
    if completion.get("ms_validation") != ms_validation:
        raise ValueError(f"conditioned ms validation changed: {unit_id}")
    return {
        **completion,
        "completion_path": str((unit_dir / "completion.json").resolve()),
        "completion_sha256": sha256_file(unit_dir / "completion.json"),
        "simulation_ms_path": str((unit_dir / "simulation.ms").resolve()),
        "simulation_ms_sha256": outputs["simulation.ms"]["sha256"],
    }


def verify_conditioned_unit(
    unit: Mapping[str, Any],
    *,
    repo_root: str | Path,
    plan_dir: str | Path,
    run_dir: str | Path,
    cosi2_binary: str | Path,
) -> dict[str, Any]:
    """Verify a legacy-v6 EAS unit or a native-v7 Han unit, never as each other."""

    root = Path(repo_root).resolve()
    plan_root = Path(plan_dir).resolve()
    unit = _portable_record(unit)
    if str(unit["demography"]) == "EAS":
        if plan_root.name != PLAN_DIRNAME:
            raise ValueError("v7 EAS verification requires the exact v7 plan boundary")
        return _verify_legacy_v6_conditioned_unit(
            unit,
            repo_root=root,
            plan_dir=plan_root.parent / LEGACY_V6_PLAN_DIRNAME,
            run_dir=plan_root.parent / LEGACY_V6_SIMULATION_DIRNAME,
            cosi2_binary=cosi2_binary,
        )

    unit_id = str(unit["unit_id"])
    unit_dir = Path(run_dir).resolve() / "units" / unit_id
    completion = _read_json(unit_dir / "completion.json", "v7 Han unit completion")
    plan_map = plan_root / PLAN_OUTPUTS["recombination_map"]
    if (
        str(unit["demography"]) != "Han"
        or completion.get("schema") != SIMULATION_UNIT_SCHEMA
        or completion.get("status") != "complete"
        or completion.get("unit_id") != unit_id
        or completion.get("ledger_row_sha256")
        != _canonical_sha256(_portable_record(unit))
        or completion.get("cosi2_binary_sha256")
        != sha256_file(Path(cosi2_binary).resolve())
        or completion.get("plan_completion_sha256")
        != sha256_file(plan_root / PLAN_COMPLETION_FILENAME)
        or completion.get("source_trajectory_path")
        != str(Path(str(unit["trajectory_path"])).resolve())
        or completion.get("source_trajectory_sha256") != str(unit["trajectory_sha256"])
        or completion.get("working_directory") != str(plan_root)
        or completion.get("plan_recombination_map_sha256") != RECOMBINATION_MAP_SHA256
        or completion.get("plan_recombination_map_sha256") != sha256_file(plan_map)
    ):
        raise ValueError(f"v7 Han unit completion changed: {unit_id}")
    children = list(unit_dir.iterdir())
    expected_output_names = {
        "model.params",
        "simulation.ms",
        "simulation.stderr.log",
        "trajectory.tsv",
        "loader_trajectory.tsv",
    }
    if any(path.is_symlink() or not path.is_file() for path in children) or {
        path.name for path in children
    } != {"completion.json", *expected_output_names}:
        raise ValueError(f"v7 Han unit inventory changed: {unit_id}")
    outputs = _file_inventory(unit_dir, exclude=("completion.json",))
    if completion.get("outputs") != outputs or set(outputs) != expected_output_names:
        raise ValueError(f"v7 Han unit output binding changed: {unit_id}")
    expected_parameters = (plan_root / PLAN_OUTPUTS["han_replay_config"]).read_text(
        encoding="utf-8"
    )
    if (unit_dir / "model.params").read_text(encoding="utf-8") != expected_parameters:
        raise ValueError(f"v7 Han unit parameters changed: {unit_id}")
    command = completion.get("command")
    parameter_path = (
        Path(str(command[2]))
        if isinstance(command, list) and len(command) > 2
        else Path()
    )
    expected_tail = [
        "-n",
        "1",
        "-r",
        str(int(unit["seed"])),
        "-u",
        "1",
        "-m",
    ]
    if (
        "recomb_file recombination_map.tsv" not in expected_parameters
        or not isinstance(command, list)
        or len(command) != 10
        or command[0] != str(Path(cosi2_binary).resolve())
        or command[1] != "-p"
        or command[3:] != expected_tail
        or parameter_path.name != "model.params"
        or not parameter_path.parent.name.startswith(f".{unit_id}.stage.")
        or parameter_path.parent.parent.resolve() != unit_dir.parent.resolve()
    ):
        raise ValueError(f"v7 Han unit command/cwd/map changed: {unit_id}")
    environment = completion.get("environment_contract")
    max_attempts = (
        environment.get("COSI_MAXATTEMPTS") if isinstance(environment, dict) else None
    )
    expected_environment = {
        "inherited_cosi_variables_removed": True,
        "COSI_MAXATTEMPTS": max_attempts,
        "COSI_LOAD_TRAJ": "checksum_bound_lexical_loader_trajectory",
        "source_trajectory_copy": "checksum_bound_full_source_to_unit_staging",
        "loader_trajectory_retained": True,
    }
    if (
        isinstance(max_attempts, bool)
        or not isinstance(max_attempts, int)
        or max_attempts <= 0
        or environment != expected_environment
    ):
        raise ValueError(f"v7 Han unit environment contract changed: {unit_id}")
    source_trajectory = Path(str(unit["trajectory_path"])).resolve()
    if (
        sha256_file(source_trajectory) != str(unit["trajectory_sha256"])
        or (unit_dir / "trajectory.tsv").read_bytes() != source_trajectory.read_bytes()
    ):
        raise ValueError(f"v7 Han full source trajectory changed: {unit_id}")
    trajectory_validation = conditioned.validate_han_trajectory(
        unit_dir / "trajectory.tsv", require_gate=True
    )
    loader_validation = _validate_han_loader_trajectory(
        unit_dir / "trajectory.tsv", unit_dir / "loader_trajectory.tsv"
    )
    if (
        completion.get("trajectory_validation") != trajectory_validation
        or completion.get("loader_trajectory_validation") != loader_validation
    ):
        raise ValueError(f"v7 Han trajectory validation changed: {unit_id}")
    ms_validation = natural_neutral.validate_cosi_ms(
        unit_dir / "simulation.ms",
        unit={"sample_haploids": analysis.PANEL_HAPLOTYPES, "seed": int(unit["seed"])},
    )
    if completion.get("ms_validation") != ms_validation:
        raise ValueError(f"v7 Han ms validation changed: {unit_id}")
    return {
        **completion,
        "completion_path": str((unit_dir / "completion.json").resolve()),
        "completion_sha256": sha256_file(unit_dir / "completion.json"),
        "simulation_ms_path": str((unit_dir / "simulation.ms").resolve()),
        "simulation_ms_sha256": outputs["simulation.ms"]["sha256"],
    }


def _publish_conditioned_simulation_failure(
    stage: Path,
    *,
    failure_root: Path,
    unit: Mapping[str, Any],
    plan_dir: Path,
    binary: Path,
    error: Exception,
    command: Sequence[str] | None,
    environment_contract: Mapping[str, Any] | None,
) -> Path:
    """Atomically retain one failed attempt outside the completed-unit inventory."""

    unit_id = str(unit["unit_id"])
    unit_failure_root = failure_root / unit_id
    unit_failure_root.mkdir(parents=True, exist_ok=True)
    children = sorted(unit_failure_root.iterdir(), key=lambda value: value.name)
    expected = [f"attempt_{index:04d}" for index in range(1, len(children) + 1)]
    if [path.name for path in children] != expected or any(
        path.is_symlink() or not path.is_dir() for path in children
    ):
        raise ValueError(f"conditioned failure attempts are not contiguous: {unit_id}")
    attempt = unit_failure_root / f"attempt_{len(children) + 1:04d}"
    stage_outputs = _file_inventory(stage)
    failure = {
        "schema": SIMULATION_FAILURE_SCHEMA,
        "status": "failed",
        "unit_id": unit_id,
        "demography": str(unit["demography"]),
        "seed": int(unit["seed"]),
        "ledger_row_sha256": _canonical_sha256(_portable_record(unit)),
        "plan_completion_sha256": sha256_file(plan_dir / PLAN_COMPLETION_FILENAME),
        "working_directory": str(plan_dir.resolve()),
        "cosi2_binary_path": str(binary.resolve()),
        "cosi2_binary_sha256": sha256_file(binary),
        "source_trajectory_path": str(Path(str(unit["trajectory_path"])).resolve()),
        "source_trajectory_sha256": str(unit["trajectory_sha256"]),
        "command": list(command) if command is not None else None,
        "environment_contract": (
            dict(environment_contract) if environment_contract is not None else None
        ),
        "error_type": type(error).__name__,
        "error_message": str(error),
        "outputs": stage_outputs,
    }
    _atomic_json(stage / "failure.json", failure)
    os.replace(stage, attempt)
    return attempt


def _run_conditioned_unit(
    repo_root: Path,
    unit: Mapping[str, Any],
    *,
    plan_dir: Path,
    run_root: Path,
    binary: Path,
    failure_root: Path,
    max_attempts: int,
    timeout_seconds: float | None,
) -> dict[str, Any]:
    unit = _portable_record(unit)
    unit_id = str(unit["unit_id"])
    if str(unit["demography"]) == "EAS":
        return {
            **verify_conditioned_unit(
                unit,
                repo_root=repo_root,
                plan_dir=plan_dir,
                run_dir=run_root,
                cosi2_binary=binary,
            ),
            "cache_hit": True,
            "artifact_generation": "exact_legacy_v6_eas_reference",
        }
    final_dir = run_root / "units" / unit_id
    if final_dir.exists():
        return {
            **verify_conditioned_unit(
                unit,
                repo_root=repo_root,
                plan_dir=plan_dir,
                run_dir=run_root,
                cosi2_binary=binary,
            ),
            "cache_hit": True,
        }
    final_dir.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{unit_id}.stage.", dir=final_dir.parent))
    started = time.monotonic()
    command: list[str] | None = None
    environment_contract: dict[str, Any] | None = None
    try:
        if str(unit["demography"]) != "Han":
            raise ValueError(f"v7 native simulation is Han-only: {unit_id}")
        parameter_text = (plan_dir / PLAN_OUTPUTS["han_replay_config"]).read_text(
            encoding="utf-8"
        )
        _atomic_text(stage / "model.params", parameter_text)
        source_trajectory = Path(str(unit["trajectory_path"])).resolve()
        if sha256_file(source_trajectory) != str(unit["trajectory_sha256"]):
            raise ValueError(f"conditioned source trajectory changed: {unit_id}")
        source_validation = conditioned.validate_han_trajectory(
            source_trajectory, require_gate=True
        )
        shutil.copyfile(source_trajectory, stage / "trajectory.tsv")
        loader_validation = _write_han_loader_trajectory(
            stage / "trajectory.tsv", stage / "loader_trajectory.tsv"
        )
        command = [
            str(binary),
            "-p",
            str((stage / "model.params").resolve()),
            "-n",
            "1",
            "-r",
            str(int(unit["seed"])),
            "-u",
            "1",
            "-m",
        ]
        environment = {
            key: value
            for key, value in os.environ.items()
            if not key.startswith("COSI_")
        }
        environment_contract = {
            "inherited_cosi_variables_removed": True,
            "COSI_MAXATTEMPTS": max_attempts,
            "COSI_LOAD_TRAJ": "checksum_bound_lexical_loader_trajectory",
            "source_trajectory_copy": "checksum_bound_full_source_to_unit_staging",
            "loader_trajectory_retained": True,
        }
        environment.update(
            {
                "COSI_MAXATTEMPTS": str(max_attempts),
                "COSI_LOAD_TRAJ": str((stage / "loader_trajectory.tsv").resolve()),
            }
        )
        with (
            (stage / "simulation.ms").open("wb") as stdout_handle,
            (stage / "simulation.stderr.log").open("wb") as stderr_handle,
        ):
            process = subprocess.run(
                command,
                cwd=plan_dir,
                env=environment,
                stdout=stdout_handle,
                stderr=stderr_handle,
                check=False,
                timeout=timeout_seconds,
            )
        if process.returncode != 0:
            raise RuntimeError(f"CoSi2 conditioned replay exited {process.returncode}")
        if (stage / "trajectory.tsv").read_bytes() != source_trajectory.read_bytes():
            raise RuntimeError(f"conditioned source copy changed: {unit_id}")
        trajectory_validation = conditioned.validate_han_trajectory(
            stage / "trajectory.tsv", require_gate=True
        )
        if trajectory_validation != source_validation:
            raise RuntimeError(f"conditioned source validation changed: {unit_id}")
        loader_validation = _validate_han_loader_trajectory(
            stage / "trajectory.tsv", stage / "loader_trajectory.tsv"
        )
        ms_validation = natural_neutral.validate_cosi_ms(
            stage / "simulation.ms",
            unit={
                "sample_haploids": analysis.PANEL_HAPLOTYPES,
                "seed": int(unit["seed"]),
            },
        )
        outputs = _file_inventory(stage)
        completion = {
            "schema": SIMULATION_UNIT_SCHEMA,
            "status": "complete",
            "unit_id": unit_id,
            "demography": str(unit["demography"]),
            "model_id": str(unit["model_id"]),
            "seed": int(unit["seed"]),
            "ledger_row_sha256": _canonical_sha256(unit),
            "source_trajectory_path": str(source_trajectory),
            "source_trajectory_sha256": str(unit["trajectory_sha256"]),
            "plan_completion_sha256": sha256_file(plan_dir / PLAN_COMPLETION_FILENAME),
            "working_directory": str(plan_dir.resolve()),
            "plan_recombination_map_sha256": sha256_file(
                plan_dir / PLAN_OUTPUTS["recombination_map"]
            ),
            "cosi2_binary_path": str(binary),
            "cosi2_binary_sha256": sha256_file(binary),
            "command": command,
            "environment_contract": environment_contract,
            "elapsed_seconds": time.monotonic() - started,
            "trajectory_validation": trajectory_validation,
            "loader_trajectory_validation": loader_validation,
            "ms_validation": ms_validation,
            "outputs": outputs,
        }
        _atomic_json(stage / "completion.json", completion)
        os.replace(stage, final_dir)
    except Exception as error:
        if stage.exists():
            _publish_conditioned_simulation_failure(
                stage,
                failure_root=failure_root,
                unit=unit,
                plan_dir=plan_dir,
                binary=binary,
                error=error,
                command=command,
                environment_contract=environment_contract,
            )
        raise
    return {
        **verify_conditioned_unit(
            unit,
            repo_root=repo_root,
            plan_dir=plan_dir,
            run_dir=run_root,
            cosi2_binary=binary,
        ),
        "cache_hit": False,
    }


def _validate_v7_simulation_root_inventory(
    run_root: Path, *, expected_han_ids: set[str], require_complete: bool
) -> None:
    if not run_root.exists():
        if require_complete:
            raise ValueError("v7 Han simulation root is absent")
        return
    children = list(run_root.iterdir())
    allowed_root_names = {"units", SIMULATION_ROOT_COMPLETION_FILENAME}
    if (
        run_root.is_symlink()
        or not run_root.is_dir()
        or not {path.name for path in children}.issubset(allowed_root_names)
        or any(path.is_symlink() for path in children)
        or any(
            path.name == "units"
            and not path.is_dir()
            or path.name == SIMULATION_ROOT_COMPLETION_FILENAME
            and not path.is_file()
            for path in children
        )
    ):
        raise ValueError("v7 Han simulation root inventory changed")
    unit_root = run_root / "units"
    unit_children = list(unit_root.iterdir()) if unit_root.is_dir() else []
    actual_ids = {path.name for path in unit_children}
    if (
        not actual_ids.issubset(expected_han_ids)
        or any(path.is_symlink() or not path.is_dir() for path in unit_children)
        or any(path.name.startswith(".") for path in unit_children)
    ):
        raise ValueError("v7 Han simulation unit inventory changed")
    completion_present = (run_root / SIMULATION_ROOT_COMPLETION_FILENAME).is_file()
    if completion_present and actual_ids != expected_han_ids:
        raise ValueError("v7 Han simulation completion exists before all units")
    if require_complete and (actual_ids != expected_han_ids or not completion_present):
        raise ValueError("v7 Han simulation root is incomplete")


def _simulation_root_contract(
    *, plan_dir: Path, ledger: pd.DataFrame, verified_rows: Sequence[Mapping[str, Any]]
) -> dict[str, Any]:
    han_rows = ledger[ledger["demography"].astype(str) == "Han"].copy()
    bindings = {
        str(row["unit_id"]): str(row["completion_sha256"])
        for row in verified_rows
        if str(row.get("demography")) == "Han"
    }
    expected_ids = set(han_rows["unit_id"].astype(str))
    if set(bindings) != expected_ids:
        raise ValueError("v7 Han completion bindings are incomplete")
    return {
        "schema": SIMULATION_ROOT_SCHEMA,
        "plan_completion_sha256": sha256_file(plan_dir / PLAN_COMPLETION_FILENAME),
        "legacy_v6_eas_root_manifest_sha256": (
            LEGACY_V6_EAS_SIMULATION_ROOT_MANIFEST_SHA256
        ),
        "legacy_v6_eas_unit_count": EXPECTED_NEUTRAL_PER_DEMOGRAPHY,
        "native_v7_han_unit_count": EXPECTED_NEUTRAL_PER_DEMOGRAPHY,
        "native_v7_han_completion_sha256": bindings,
        "union_rule": "reference_exact_legacy_v6_eas_root_plus_native_v7_han_root",
    }


def run_conditioned_simulations(
    repo_root: str | Path,
    *,
    plan_dir: str | Path,
    proposal_dir: str | Path,
    selection_dir: str | Path,
    screen_dir: str | Path,
    run_dir: str | Path,
    cosi2_binary: str | Path,
    max_workers: int = MAX_WORKERS,
    max_attempts: int = natural_neutral.DEFAULT_COSI_MAXATTEMPTS,
    timeout_seconds: float | None = None,
    unit_ids: Sequence[str] | None = None,
) -> pd.DataFrame:
    root = Path(repo_root).resolve()
    eas_bundle = verify_eas_global_selection(
        root,
        plan_dir=plan_dir,
        proposal_dir=proposal_dir,
        selection_dir=selection_dir,
    )
    ledger = build_conditioned_simulation_ledger(
        root,
        plan_dir=plan_dir,
        proposal_dir=proposal_dir,
        selection_dir=selection_dir,
        screen_dir=screen_dir,
        cosi2_binary=cosi2_binary,
        verified_eas_bundle=eas_bundle,
    )
    requested = (
        list(dict.fromkeys(str(value) for value in unit_ids))
        if unit_ids is not None
        else list(ledger["unit_id"].astype(str))
    )
    unknown = sorted(set(requested).difference(ledger["unit_id"].astype(str)))
    if unknown:
        raise ValueError("unknown conditioned unit IDs: " + ", ".join(unknown))
    order = {unit_id: index for index, unit_id in enumerate(requested)}
    selected = ledger[ledger["unit_id"].astype(str).isin(requested)].copy()
    selected["_order"] = selected["unit_id"].map(order)
    selected = selected.sort_values("_order")
    binary = Path(cosi2_binary).resolve()
    run_root = Path(run_dir).resolve()
    failure_root = run_root.parent / SIMULATION_FAILURE_DIRNAME
    expected_han_ids = set(
        ledger.loc[ledger["demography"].astype(str) == "Han", "unit_id"].astype(str)
    )
    _validate_v7_simulation_root_inventory(
        run_root, expected_han_ids=expected_han_ids, require_complete=False
    )
    (run_root / "units").mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    failures: list[tuple[str, str]] = []
    with ThreadPoolExecutor(max_workers=_positive_workers(max_workers)) as pool:
        futures = {
            pool.submit(
                _run_conditioned_unit,
                root,
                record,
                plan_dir=Path(plan_dir).resolve(),
                run_root=run_root,
                binary=binary,
                failure_root=failure_root,
                max_attempts=max_attempts,
                timeout_seconds=timeout_seconds,
            ): str(record["unit_id"])
            for record in selected.drop(columns="_order").to_dict(orient="records")
        }
        for future in as_completed(futures):
            try:
                rows.append(future.result())
            except Exception as error:  # noqa: BLE001
                failures.append((futures[future], str(error)))
    if failures:
        detail = "; ".join(f"{key}: {value}" for key, value in failures[:10])
        raise RuntimeError(f"{len(failures)} conditioned units failed: {detail}")
    result = pd.DataFrame(rows)
    result["_order"] = result["unit_id"].map(order)
    result = result.sort_values("_order").drop(columns="_order").reset_index(drop=True)
    completed_han_ids = {
        path.name
        for path in (run_root / "units").iterdir()
        if path.is_dir() and not path.is_symlink()
    }
    if completed_han_ids == expected_han_ids:
        verified_all: list[dict[str, Any]] = []
        for row in ledger.to_dict(orient="records"):
            verified_all.append(
                verify_conditioned_unit(
                    row,
                    repo_root=root,
                    plan_dir=plan_dir,
                    run_dir=run_root,
                    cosi2_binary=binary,
                )
            )
        contract = _simulation_root_contract(
            plan_dir=Path(plan_dir).resolve(), ledger=ledger, verified_rows=verified_all
        )
        completion_path = run_root / SIMULATION_ROOT_COMPLETION_FILENAME
        expected_completion = {
            "schema": SIMULATION_ROOT_SCHEMA,
            "status": "complete",
            "contract": contract,
            "contract_sha256": _canonical_sha256(contract),
        }
        if completion_path.is_file():
            if (
                _read_json(completion_path, "v7 simulation completion")
                != expected_completion
            ):
                raise ValueError("v7 simulation completion is stale")
        else:
            _atomic_json(completion_path, expected_completion)
        _validate_v7_simulation_root_inventory(
            run_root, expected_han_ids=expected_han_ids, require_complete=True
        )
    return result


def collect_conditioned_units(
    repo_root: str | Path,
    *,
    plan_dir: str | Path,
    proposal_dir: str | Path,
    selection_dir: str | Path,
    screen_dir: str | Path,
    run_dir: str | Path,
    cosi2_binary: str | Path,
) -> pd.DataFrame:
    root = Path(repo_root).resolve()
    eas_bundle = verify_eas_global_selection(
        root,
        plan_dir=plan_dir,
        proposal_dir=proposal_dir,
        selection_dir=selection_dir,
    )
    ledger = build_conditioned_simulation_ledger(
        root,
        plan_dir=plan_dir,
        proposal_dir=proposal_dir,
        selection_dir=selection_dir,
        screen_dir=screen_dir,
        cosi2_binary=cosi2_binary,
        verified_eas_bundle=eas_bundle,
    )
    eas_diagnostics = eas_bundle["diagnostics"]
    expected_han_ids = set(
        ledger.loc[ledger["demography"].astype(str) == "Han", "unit_id"].astype(str)
    )
    _validate_v7_simulation_root_inventory(
        Path(run_dir).resolve(),
        expected_han_ids=expected_han_ids,
        require_complete=True,
    )
    rows: list[dict[str, Any]] = []
    verified_units: list[dict[str, Any]] = []
    for unit in ledger.to_dict(orient="records"):
        verified = verify_conditioned_unit(
            unit,
            repo_root=root,
            plan_dir=plan_dir,
            run_dir=run_dir,
            cosi2_binary=cosi2_binary,
        )
        verified_units.append(verified)
        ms = verified["ms_validation"]
        sample_count = (
            int(ms["focal_alt_count_if_unambiguous"])
            if int(ms["reported_focal_position_count"]) == 1
            else 0
        )
        if unit["demography"] == "EAS":
            final_af = conditioned.EAS_SELECTED_PRESENT_AF
            birth = int(unit["birth_generation"])
            gate_af = np.nan
            candidate_id = ""
            candidate_order = -1
            candidate_source = ""
            source_particle = str(int(unit["global_proposal_index"]))
            occurrence = int(unit["accepted_index"])
            multiplicity = int(unit["source_particle_draw_multiplicity"])
            log_weight = float(unit["log_weight"])
            approximation = True
            conditioning_label = "exact_present_population_count_global_importance_bank"
        else:
            validation = verified["trajectory_validation"]
            final_af = float(validation["present_chb_af"])
            birth = conditioned.HAN_BIRTH_GENERATION
            gate_af = float(validation["han_generation_645_chb_af"])
            candidate_id = str(unit["candidate_id"])
            candidate_order = int(unit["candidate_order"])
            candidate_source = str(unit["source_kind"])
            source_particle = ""
            occurrence = -1
            multiplicity = -1
            log_weight = np.nan
            approximation = False
            conditioning_label = "g645_chb_af_gate_and_present_strict_segregation"
        settings = DEMOGRAPHY_SETTINGS[str(unit["demography"])]
        rows.append(
            {
                "unit_id": str(unit["unit_id"]),
                "demography": str(unit["demography"]),
                "demography_id": settings["demography_id"],
                "model_id": str(unit["model_id"]),
                "simulation_class": "neutral",
                "selection_coefficient": 0.0,
                "seed": int(unit["seed"]),
                "generation_time_years": settings["generation_time_years"],
                "present_ne": settings["present_ne"],
                "final_population_af": final_af,
                "sample_alt_count": sample_count,
                "sample_af": sample_count / analysis.PANEL_HAPLOTYPES,
                "ms_path": verified["simulation_ms_path"],
                "ms_sha256": verified["simulation_ms_sha256"],
                "source_completion_path": verified["completion_path"],
                "source_completion_sha256": verified["completion_sha256"],
                "conditioning": conditioning_label,
                "birth_generation": birth,
                "allele_age_years": birth * settings["generation_time_years"],
                "eas_selected_log_importance_weight": log_weight,
                "eas_global_effective_sample_size": (
                    eas_diagnostics["global_effective_sample_size"]
                    if unit["demography"] == "EAS"
                    else np.nan
                ),
                "eas_global_weight_sum_over_max": (
                    eas_diagnostics["global_weight_sum_over_max"]
                    if unit["demography"] == "EAS"
                    else np.nan
                ),
                "eas_global_max_normalized_weight": (
                    eas_diagnostics["global_max_normalized_weight"]
                    if unit["demography"] == "EAS"
                    else np.nan
                ),
                "eas_time_reversal_is_approximation": approximation,
                "eas_source_particle_id": source_particle,
                "eas_draw_occurrence": occurrence,
                "eas_particle_multiplicity": multiplicity,
                "han_candidate_id": candidate_id,
                "han_candidate_order": candidate_order,
                "han_candidate_source_kind": candidate_source,
                "han_generation_645_chb_af": gate_af,
            }
        )
    simulation_contract = _simulation_root_contract(
        plan_dir=Path(plan_dir).resolve(),
        ledger=ledger,
        verified_rows=verified_units,
    )
    expected_completion = {
        "schema": SIMULATION_ROOT_SCHEMA,
        "status": "complete",
        "contract": simulation_contract,
        "contract_sha256": _canonical_sha256(simulation_contract),
    }
    observed_completion = _read_json(
        Path(run_dir).resolve() / SIMULATION_ROOT_COMPLETION_FILENAME,
        "v7 simulation completion",
    )
    if observed_completion != expected_completion:
        raise ValueError("v7 simulation root completion changed")
    return pd.DataFrame(rows, columns=UNIT_COLUMNS)


def build_unit_inventory(
    repo_root: str | Path,
    *,
    plan_dir: str | Path,
    proposal_dir: str | Path,
    selection_dir: str | Path,
    screen_dir: str | Path,
    run_dir: str | Path,
    cosi2_binary: str | Path,
    selected_dir: str | Path | None = None,
) -> pd.DataFrame:
    selected = _selected_with_conditioning_columns(
        natural_workflow.bind_selected_units(repo_root, selected_dir=selected_dir)
    )
    neutral = collect_conditioned_units(
        repo_root,
        plan_dir=plan_dir,
        proposal_dir=proposal_dir,
        selection_dir=selection_dir,
        screen_dir=screen_dir,
        run_dir=run_dir,
        cosi2_binary=cosi2_binary,
    )
    return validate_unit_inventory(pd.concat([neutral, selected], ignore_index=True))


def plan_study(
    repo_root: str | Path,
    *,
    work_dir: str | Path | None = None,
    natural_run_dir: str | Path | None = None,
) -> dict[str, Any]:
    root = Path(repo_root).resolve()
    paths = _resolve_paths(root, work_dir, None)
    natural_workflow.bind_selected_units(root)
    completion = generate_plan(
        root,
        output_dir=paths["plan"],
        proposal_dir=paths["eas_proposals"],
        natural_run_dir=natural_run_dir,
    )
    return {
        "status": "planned",
        "plan_dir": str(paths["plan"]),
        "cache_hit": completion["cache_hit"],
        "eas_proposal_blocks": conditioned.EAS_PROPOSAL_BLOCK_COUNT,
        "han_candidates": conditioned.HAN_CANDIDATE_COUNT,
    }


def eas_bridge_study(
    repo_root: str | Path,
    *,
    work_dir: str | Path | None = None,
    max_workers: int = MAX_WORKERS,
    block_indexes: Sequence[int] | None = None,
) -> dict[str, Any]:
    root = Path(repo_root).resolve()
    paths = _resolve_paths(root, work_dir, None)
    verify_plan(
        root,
        plan_dir=paths["plan"],
        proposal_dir=paths["eas_proposals"],
    )
    base_manifest = _read_json(
        paths["plan"] / PLAN_OUTPUTS["eas_base_bank_manifest"],
        "EAS base-bank manifest",
    )
    extension_plan = _read_json(
        paths["plan"] / PLAN_OUTPUTS["eas_extension_plan"],
        "EAS extension plan",
    )
    status = run_eas_proposal_blocks(
        root,
        proposal_dir=paths["eas_proposals"],
        base_manifest=base_manifest,
        extension_plan=extension_plan,
        max_workers=max_workers,
        block_indexes=block_indexes,
    )
    if block_indexes is not None:
        return {
            "status": "partial",
            "blocks_completed_this_call": len(status),
            "selection": "not_attempted_for_partial_block_request",
        }
    completion = write_eas_global_selection(
        root,
        plan_dir=paths["plan"],
        proposal_dir=paths["eas_proposals"],
        selection_dir=paths["eas_selection"],
        replay_dir=paths["eas_parallel_replays"],
        max_workers=max_workers,
    )
    return {
        "status": "complete",
        "proposal_blocks": conditioned.EAS_PROPOSAL_BLOCK_COUNT,
        "extension_blocks_processed_this_call": len(status),
        "selected_occurrences": len(completion["selected"]),
        "global_diagnostics": completion["diagnostics"],
    }


def han_screen_study(
    repo_root: str | Path,
    *,
    work_dir: str | Path | None = None,
    cosi2_binary: str | Path | None = None,
    max_workers: int = MAX_WORKERS,
    max_attempts: int = natural_neutral.DEFAULT_COSI_MAXATTEMPTS,
    timeout_seconds: float | None = None,
    candidate_ids: Sequence[str] | None = None,
) -> dict[str, Any]:
    root = Path(repo_root).resolve()
    paths = _resolve_paths(root, work_dir, None)
    binary = (
        Path(cosi2_binary).resolve()
        if cosi2_binary is not None
        else root / DEFAULT_COSI2_BINARY_RELATIVE
    )
    verify_plan(
        root,
        plan_dir=paths["plan"],
        proposal_dir=paths["eas_proposals"],
    )
    plan_completion_path = paths["plan"] / PLAN_COMPLETION_FILENAME
    current_plan_schema = (
        _read_json(plan_completion_path, "Han-screen plan completion").get("schema")
        if plan_completion_path.is_file()
        else None
    )
    if current_plan_schema == PLAN_SCHEMA:
        if candidate_ids is not None:
            raise ValueError(
                "v7 freezes the completed v6 Han screen; no new candidates"
            )
        state = _verify_legacy_v6_han_state(
            legacy_plan_dir=paths["legacy_v6_plan"],
            screen_dir=paths["han_screen"],
            cosi2_binary=binary,
        )
        return {
            "status": "complete_legacy_v6_han_screen_reused",
            "screened": int(state["screened_candidate_count"]),
            "accepted": int(state["accepted_rows"]),
            "selection_completion_sha256": state["selection_completion_sha256"],
        }
    if candidate_ids is None:
        ledger = pd.read_csv(paths["plan"] / PLAN_OUTPUTS["han_candidates"], sep="\t")
        candidate_root = paths["han_screen"] / "candidates"
        completed_names = (
            {
                path.name
                for path in candidate_root.iterdir()
                if path.is_dir() and (path / "completion.json").is_file()
            }
            if candidate_root.is_dir()
            else set()
        )
        order_by_id = dict(
            zip(
                ledger["candidate_id"].astype(str),
                ledger["candidate_order"].astype(int),
                strict=True,
            )
        )
        unknown_completed = sorted(completed_names.difference(order_by_id))
        if unknown_completed:
            raise ValueError(
                "Han screen root contains unknown completed candidates: "
                + ", ".join(unknown_completed[:10])
            )
        greatest_completed_order = max(
            (order_by_id[value] for value in completed_names), default=-1
        )
        completed: list[pd.DataFrame] = []
        for start in range(0, len(ledger), 100):
            batch_ids = list(
                ledger.iloc[start : start + 100]["candidate_id"].astype(str)
            )
            batch = run_han_screens(
                root,
                plan_dir=paths["plan"],
                screen_dir=paths["han_screen"],
                cosi2_binary=binary,
                max_workers=max_workers,
                max_attempts=max_attempts,
                timeout_seconds=timeout_seconds,
                candidate_ids=batch_ids,
            )
            completed.append(batch)
            prefix = pd.concat(completed, ignore_index=True)
            if (
                int(prefix["gate_passed"].sum()) >= EXPECTED_NEUTRAL_PER_DEMOGRAPHY
                and start + len(batch) - 1 >= greatest_completed_order
            ):
                break
        status = pd.concat(completed, ignore_index=True)
    else:
        status = run_han_screens(
            root,
            plan_dir=paths["plan"],
            screen_dir=paths["han_screen"],
            cosi2_binary=binary,
            max_workers=max_workers,
            max_attempts=max_attempts,
            timeout_seconds=timeout_seconds,
            candidate_ids=candidate_ids,
        )
    result = {
        "status": "complete" if candidate_ids is None else "partial",
        "candidates_completed_this_call": len(status),
        "gate_passed_this_call": int(status["gate_passed"].sum()),
    }
    if candidate_ids is None:
        accepted = freeze_han_acceptances(
            root,
            plan_dir=paths["plan"],
            screen_dir=paths["han_screen"],
            cosi2_binary=binary,
        )
        result["accepted"] = len(accepted)
        result["screening_rule"] = (
            "fixed_ledger_contiguous_100_candidate_batches_until_100th_pass"
        )
    return result


def simulate_study(
    repo_root: str | Path,
    *,
    work_dir: str | Path | None = None,
    cosi2_binary: str | Path | None = None,
    max_workers: int = MAX_WORKERS,
    max_attempts: int = natural_neutral.DEFAULT_COSI_MAXATTEMPTS,
    timeout_seconds: float | None = None,
    unit_ids: Sequence[str] | None = None,
) -> pd.DataFrame:
    root = Path(repo_root).resolve()
    paths = _resolve_paths(root, work_dir, None)
    binary = (
        Path(cosi2_binary).resolve()
        if cosi2_binary is not None
        else root / DEFAULT_COSI2_BINARY_RELATIVE
    )
    return run_conditioned_simulations(
        root,
        plan_dir=paths["plan"],
        proposal_dir=paths["eas_proposals"],
        selection_dir=paths["eas_selection"],
        screen_dir=paths["han_screen"],
        run_dir=paths["simulations"],
        cosi2_binary=binary,
        max_workers=max_workers,
        max_attempts=max_attempts,
        timeout_seconds=timeout_seconds,
        unit_ids=unit_ids,
    )


def decode_study(
    repo_root: str | Path,
    *,
    work_dir: str | Path | None = None,
    cosi2_binary: str | Path | None = None,
    decoder_bin: str | Path | None = None,
    max_workers: int = MAX_WORKERS,
    unit_ids: Sequence[str] | None = None,
) -> pd.DataFrame:
    root = Path(repo_root).resolve()
    paths = _resolve_paths(root, work_dir, None)
    cosi2 = (
        Path(cosi2_binary).resolve()
        if cosi2_binary is not None
        else root / DEFAULT_COSI2_BINARY_RELATIVE
    )
    decoder = (
        Path(decoder_bin).resolve()
        if decoder_bin is not None
        else root / DEFAULT_GAMMA_BINARY_RELATIVE
    )
    units = build_unit_inventory(
        root,
        plan_dir=paths["plan"],
        proposal_dir=paths["eas_proposals"],
        selection_dir=paths["eas_selection"],
        screen_dir=paths["han_screen"],
        run_dir=paths["simulations"],
        cosi2_binary=cosi2,
    )
    return run_decodes(
        units,
        decode_root=paths["decode"],
        decoder_bin=decoder,
        max_workers=max_workers,
        unit_ids=unit_ids,
    )


def analyze_study(
    repo_root: str | Path,
    *,
    work_dir: str | Path | None = None,
    results_dir: str | Path | None = None,
    cosi2_binary: str | Path | None = None,
    decoder_bin: str | Path | None = None,
) -> dict[str, Any]:
    root = Path(repo_root).resolve()
    paths = _resolve_paths(root, work_dir, results_dir)
    cosi2 = (
        Path(cosi2_binary).resolve()
        if cosi2_binary is not None
        else root / DEFAULT_COSI2_BINARY_RELATIVE
    )
    decoder = (
        Path(decoder_bin).resolve()
        if decoder_bin is not None
        else root / DEFAULT_GAMMA_BINARY_RELATIVE
    )
    units = build_unit_inventory(
        root,
        plan_dir=paths["plan"],
        proposal_dir=paths["eas_proposals"],
        selection_dir=paths["eas_selection"],
        screen_dir=paths["han_screen"],
        run_dir=paths["simulations"],
        cosi2_binary=cosi2,
    )
    scores, decode_records = collect_decoded_scores(
        units, decode_root=paths["decode"], decoder_bin=decoder
    )
    return write_results_bundle(
        root,
        results_dir=paths["results"],
        units=units,
        scores=scores,
        decode_records=decode_records,
        decoder_bin=decoder,
        cosi2_binary=cosi2,
    )


def workflow_status(
    repo_root: str | Path,
    *,
    work_dir: str | Path | None = None,
    results_dir: str | Path | None = None,
    cosi2_binary: str | Path | None = None,
    decoder_bin: str | Path | None = None,
) -> dict[str, Any]:
    root = Path(repo_root).resolve()
    paths = _resolve_paths(root, work_dir, results_dir)
    binary = (
        Path(cosi2_binary).resolve()
        if cosi2_binary is not None
        else root / DEFAULT_COSI2_BINARY_RELATIVE
    )
    decoder = (
        Path(decoder_bin).resolve()
        if decoder_bin is not None
        else root / DEFAULT_GAMMA_BINARY_RELATIVE
    )
    result: dict[str, Any] = {
        "plan": "absent",
        "legacy_pilot_plan": (
            "present_immutable"
            if (paths["work"] / LEGACY_PLAN_DIRNAME).is_dir()
            else "absent"
        ),
        "legacy_v2_plan": (
            "present_immutable" if paths["legacy_v2_plan"].is_dir() else "absent"
        ),
        "eas_proposal_blocks": {
            "complete": 0,
            "planned": conditioned.EAS_PROPOSAL_BLOCK_COUNT,
            "base_complete": 0,
            "base_planned": conditioned.EAS_BASE_PROPOSAL_BLOCK_COUNT,
            "extension_complete": 0,
            "extension_planned": conditioned.EAS_V2_EXTENSION_PROPOSAL_BLOCK_COUNT,
        },
        "eas_parallel_replay": {
            "complete_selected_blocks": 0,
            "root_completion": "absent",
        },
        "eas_selection": "absent",
        "han_screens": {"complete": 0, "planned": conditioned.HAN_CANDIDATE_COUNT},
        "conditioned_simulations": {
            "complete": 0,
            "planned": 200,
            "legacy_v6_eas_complete": 0,
            "native_v7_han_complete": 0,
            "native_v7_root_completion": "absent",
        },
        "decode": {"complete": 0, "planned": EXPECTED_TOTAL_UNITS},
        "results": "absent",
    }
    if paths["plan"].is_dir():
        try:
            verify_plan(root, plan_dir=paths["plan"])
        except Exception as error:  # noqa: BLE001
            result["plan"] = f"invalid: {error}"
            return result
        result["plan"] = "complete"
    base_complete = sum(
        (_eas_block_dir(paths["eas_proposals"], index) / "completion.json").is_file()
        for index in range(conditioned.EAS_BASE_PROPOSAL_BLOCK_COUNT)
    )
    extension_complete = sum(
        (_eas_block_dir(paths["eas_proposals"], index) / "completion.json").is_file()
        for index in range(
            conditioned.EAS_BASE_PROPOSAL_BLOCK_COUNT,
            conditioned.EAS_PROPOSAL_BLOCK_COUNT,
        )
    )
    result["eas_proposal_blocks"].update(
        {
            "complete": base_complete + extension_complete,
            "base_complete": base_complete,
            "extension_complete": extension_complete,
        }
    )
    replay_root = paths["eas_parallel_replays"]
    result["eas_parallel_replay"]["complete_selected_blocks"] = (
        len(list(replay_root.glob("block_*/completion.json")))
        if replay_root.is_dir()
        else 0
    )
    if (replay_root / EAS_PARALLEL_REPLAY_COMPLETION_FILENAME).is_file():
        result["eas_parallel_replay"]["root_completion"] = (
            "complete_present_not_deep_verified"
        )
    if (paths["eas_selection"] / EAS_SELECTION_COMPLETION).is_file():
        result["eas_selection"] = "complete_present_not_deep_verified"
    candidate_root = paths["han_screen"] / "candidates"
    result["han_screens"]["complete"] = (
        len(list(candidate_root.glob("*/completion.json")))
        if candidate_root.is_dir()
        else 0
    )
    legacy_run_root = paths["legacy_v6_simulations"] / "units"
    run_root = paths["simulations"] / "units"
    legacy_eas_complete = (
        len(list(legacy_run_root.glob("eas_*/completion.json")))
        if legacy_run_root.is_dir()
        else 0
    )
    native_han_complete = (
        len(list(run_root.glob("han_*/completion.json"))) if run_root.is_dir() else 0
    )
    result["conditioned_simulations"].update(
        {
            "complete": legacy_eas_complete + native_han_complete,
            "legacy_v6_eas_complete": legacy_eas_complete,
            "native_v7_han_complete": native_han_complete,
            "native_v7_root_completion": (
                "complete_present_not_deep_verified"
                if (
                    paths["simulations"] / SIMULATION_ROOT_COMPLETION_FILENAME
                ).is_file()
                else "absent"
            ),
        }
    )
    if result["conditioned_simulations"]["complete"] == 200:
        try:
            units = build_unit_inventory(
                root,
                plan_dir=paths["plan"],
                proposal_dir=paths["eas_proposals"],
                selection_dir=paths["eas_selection"],
                screen_dir=paths["han_screen"],
                run_dir=paths["simulations"],
                cosi2_binary=binary,
            )
            decodes = decode_status(
                units, decode_root=paths["decode"], decoder_bin=decoder
            )
            result["decode"] = decodes["status"].value_counts().to_dict()
        except Exception as error:  # noqa: BLE001
            result["decode"] = {"invalid_upstream": str(error)}
    if paths["results"].exists():
        try:
            verify_results_bundle(root, results_dir=paths["results"])
        except Exception as error:  # noqa: BLE001
            result["results"] = f"invalid: {error}"
        else:
            result["results"] = "complete"
    return result


def verify_study(
    repo_root: str | Path,
    *,
    work_dir: str | Path | None = None,
    results_dir: str | Path | None = None,
    cosi2_binary: str | Path | None = None,
    decoder_bin: str | Path | None = None,
) -> dict[str, Any]:
    root = Path(repo_root).resolve()
    paths = _resolve_paths(root, work_dir, results_dir)
    cosi2 = (
        Path(cosi2_binary).resolve()
        if cosi2_binary is not None
        else root / DEFAULT_COSI2_BINARY_RELATIVE
    )
    decoder = (
        Path(decoder_bin).resolve()
        if decoder_bin is not None
        else root / DEFAULT_GAMMA_BINARY_RELATIVE
    )
    verify_plan(
        root,
        plan_dir=paths["plan"],
        proposal_dir=paths["eas_proposals"],
    )
    base_manifest = _read_json(
        paths["plan"] / PLAN_OUTPUTS["eas_base_bank_manifest"],
        "EAS base-bank manifest",
    )
    extension_plan = _read_json(
        paths["plan"] / PLAN_OUTPUTS["eas_extension_plan"],
        "EAS extension plan",
    )
    blocks = load_eas_proposal_blocks(
        root,
        proposal_dir=paths["eas_proposals"],
        base_manifest=base_manifest,
        extension_plan=extension_plan,
    )
    verify_eas_global_selection(
        root,
        plan_dir=paths["plan"],
        proposal_dir=paths["eas_proposals"],
        selection_dir=paths["eas_selection"],
    )
    expected_han_ids = {
        f"{conditioned.HAN_MODEL_ID}_r{index:03d}"
        for index in range(EXPECTED_NEUTRAL_PER_DEMOGRAPHY)
    }
    _validate_v7_simulation_root_inventory(
        paths["simulations"],
        expected_han_ids=expected_han_ids,
        require_complete=True,
    )
    units = build_unit_inventory(
        root,
        plan_dir=paths["plan"],
        proposal_dir=paths["eas_proposals"],
        selection_dir=paths["eas_selection"],
        screen_dir=paths["han_screen"],
        run_dir=paths["simulations"],
        cosi2_binary=cosi2,
    )
    expected_decode_ids = set(units["unit_id"].astype(str))
    decode_children = (
        list(paths["decode"].iterdir()) if paths["decode"].is_dir() else []
    )
    if {path.name for path in decode_children} != expected_decode_ids or any(
        path.is_symlink() or not path.is_dir() for path in decode_children
    ):
        raise ValueError("decode root has an unexpected exact-202 inventory")
    scores, decode_records = collect_decoded_scores(
        units, decode_root=paths["decode"], decoder_bin=decoder
    )
    contract = _result_contract(
        root,
        units,
        scores,
        decode_records,
        decoder_bin=decoder,
        cosi2_binary=cosi2,
    )
    results = verify_results_bundle(
        root, results_dir=paths["results"], expected_contract=contract
    )
    return {
        "status": "complete",
        "eas_proposal_blocks": len(blocks),
        "neutral_units": 200,
        "selected_units": 2,
        "decoded_units": len(decode_records),
        "score_rows": len(scores),
        "results_contract_sha256": results["contract_sha256"],
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "phase",
        choices=(
            "plan",
            "eas-bridge",
            "han-screen",
            "simulate",
            "decode",
            "analyze",
            "status",
            "verify",
            "all",
        ),
    )
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--work-dir", type=Path)
    parser.add_argument("--results-dir", type=Path)
    parser.add_argument("--cosi2-binary", type=Path)
    parser.add_argument("--decoder-bin", type=Path)
    parser.add_argument("--workers", type=int, default=MAX_WORKERS)
    parser.add_argument(
        "--max-attempts", type=int, default=natural_neutral.DEFAULT_COSI_MAXATTEMPTS
    )
    parser.add_argument("--timeout-seconds", type=float)
    parser.add_argument("--block-index", type=int, action="append")
    parser.add_argument("--candidate-id", action="append")
    parser.add_argument("--unit-id", action="append")
    return parser


def _print_status(frame: pd.DataFrame, label: str) -> None:
    counts = frame.get("status", pd.Series(dtype=str)).value_counts().to_dict()
    print(json.dumps({label: counts}, sort_keys=True))


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    root = args.repo_root.resolve()
    if args.phase == "plan":
        output = plan_study(root, work_dir=args.work_dir)
    elif args.phase == "eas-bridge":
        output = eas_bridge_study(
            root,
            work_dir=args.work_dir,
            max_workers=args.workers,
            block_indexes=args.block_index,
        )
    elif args.phase == "han-screen":
        output = han_screen_study(
            root,
            work_dir=args.work_dir,
            cosi2_binary=args.cosi2_binary,
            max_workers=args.workers,
            max_attempts=args.max_attempts,
            timeout_seconds=args.timeout_seconds,
            candidate_ids=args.candidate_id,
        )
    elif args.phase == "simulate":
        frame = simulate_study(
            root,
            work_dir=args.work_dir,
            cosi2_binary=args.cosi2_binary,
            max_workers=args.workers,
            max_attempts=args.max_attempts,
            timeout_seconds=args.timeout_seconds,
            unit_ids=args.unit_id,
        )
        _print_status(frame, "simulate")
        return 0
    elif args.phase == "decode":
        frame = decode_study(
            root,
            work_dir=args.work_dir,
            cosi2_binary=args.cosi2_binary,
            decoder_bin=args.decoder_bin,
            max_workers=args.workers,
            unit_ids=args.unit_id,
        )
        _print_status(frame, "decode")
        return 0
    elif args.phase == "analyze":
        completion = analyze_study(
            root,
            work_dir=args.work_dir,
            results_dir=args.results_dir,
            cosi2_binary=args.cosi2_binary,
            decoder_bin=args.decoder_bin,
        )
        output = {"status": completion["status"], "cache_hit": completion["cache_hit"]}
    elif args.phase == "status":
        output = workflow_status(
            root,
            work_dir=args.work_dir,
            results_dir=args.results_dir,
            cosi2_binary=args.cosi2_binary,
            decoder_bin=args.decoder_bin,
        )
    elif args.phase == "verify":
        output = verify_study(
            root,
            work_dir=args.work_dir,
            results_dir=args.results_dir,
            cosi2_binary=args.cosi2_binary,
            decoder_bin=args.decoder_bin,
        )
    else:
        plan_study(root, work_dir=args.work_dir)
        eas_bridge_study(root, work_dir=args.work_dir, max_workers=args.workers)
        han_screen_study(
            root,
            work_dir=args.work_dir,
            cosi2_binary=args.cosi2_binary,
            max_workers=args.workers,
            max_attempts=args.max_attempts,
            timeout_seconds=args.timeout_seconds,
        )
        simulate_study(
            root,
            work_dir=args.work_dir,
            cosi2_binary=args.cosi2_binary,
            max_workers=args.workers,
            max_attempts=args.max_attempts,
            timeout_seconds=args.timeout_seconds,
        )
        decode_study(
            root,
            work_dir=args.work_dir,
            cosi2_binary=args.cosi2_binary,
            decoder_bin=args.decoder_bin,
            max_workers=args.workers,
        )
        analyze_study(
            root,
            work_dir=args.work_dir,
            results_dir=args.results_dir,
            cosi2_binary=args.cosi2_binary,
            decoder_bin=args.decoder_bin,
        )
        output = verify_study(
            root,
            work_dir=args.work_dir,
            results_dir=args.results_dir,
            cosi2_binary=args.cosi2_binary,
            decoder_bin=args.decoder_bin,
        )
    print(json.dumps(output, indent=2, sort_keys=True))
    return 0


__all__ = [
    "DEFAULT_RESULTS_RELATIVE",
    "DEFAULT_WORK_RELATIVE",
    "RESULT_OUTPUTS",
    "RESULTS_COMPLETION_FILENAME",
    "analyze_study",
    "build_parser",
    "build_unit_inventory",
    "conditioning_contrasts",
    "decode_study",
    "eas_bridge_study",
    "han_screen_study",
    "main",
    "plan_study",
    "run_eas_parallel_selected_replays",
    "simulate_study",
    "validate_unit_inventory",
    "verify_results_bundle",
    "verify_study",
    "workflow_status",
    "write_results_bundle",
]


def run_decodes(
    units: pd.DataFrame,
    *,
    decode_root: str | Path,
    decoder_bin: str | Path,
    max_workers: int = MAX_WORKERS,
    unit_ids: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Decode the exact inventory with one-thread, all-pairs Gamma processes."""

    inventory = validate_unit_inventory(units)
    workers = _positive_workers(max_workers)
    requested = (
        list(dict.fromkeys(str(value) for value in unit_ids))
        if unit_ids is not None
        else list(inventory["unit_id"].astype(str))
    )
    unknown = sorted(set(requested).difference(inventory["unit_id"].astype(str)))
    if unknown:
        raise ValueError("unknown decode unit IDs: " + ", ".join(unknown))
    selected = inventory[inventory["unit_id"].astype(str).isin(requested)].copy()
    order = {unit_id: index for index, unit_id in enumerate(requested)}
    selected["_order"] = selected["unit_id"].map(order)
    selected = selected.sort_values("_order", kind="mergesort")
    root = Path(decode_root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    binary = Path(decoder_bin).resolve()

    def run_one(record: Mapping[str, Any]) -> dict[str, Any]:
        artifacts = gamma_decode.decode_cosi_unit(
            record["ms_path"],
            root / str(record["unit_id"]),
            binary,
            unit_id=str(record["unit_id"]),
            demography_id=str(record["demography_id"]),
            present_ne=float(record["present_ne"]),
            generation_time_years=float(record["generation_time_years"]),
            threads=1,
            expected_ms_sha256=str(record["ms_sha256"]),
        )
        verified = natural_workflow.verify_decode_unit(
            record,
            decode_root=root,
            decoder_bin=binary,
        )
        return {**verified, "cache_hit": artifacts.cache_hit}

    rows: list[dict[str, Any]] = []
    failures: list[tuple[str, str]] = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(run_one, record): str(record["unit_id"])
            for record in selected.drop(columns="_order").to_dict(orient="records")
        }
        for future in as_completed(futures):
            try:
                rows.append(future.result())
            except Exception as error:  # noqa: BLE001 - preserve all unit failures
                failures.append((futures[future], str(error)))
    if failures:
        detail = "; ".join(f"{unit_id}: {error}" for unit_id, error in failures[:10])
        raise RuntimeError(f"{len(failures)} Gamma decode units failed: {detail}")
    result = pd.DataFrame(rows)
    result["_order"] = result["unit_id"].map(order)
    return (
        result.sort_values("_order", kind="mergesort")
        .drop(columns="_order")
        .reset_index(drop=True)
    )


def decode_status(
    units: pd.DataFrame,
    *,
    decode_root: str | Path,
    decoder_bin: str | Path,
) -> pd.DataFrame:
    """Return read-only status for every expected decode unit."""

    inventory = validate_unit_inventory(units)
    root = Path(decode_root).resolve()
    rows: list[dict[str, Any]] = []
    for record in inventory.to_dict(orient="records"):
        unit_id = str(record["unit_id"])
        try:
            verified = natural_workflow.verify_decode_unit(
                record,
                decode_root=root,
                decoder_bin=decoder_bin,
            )
        except Exception as error:  # noqa: BLE001 - status reports invalid caches
            rows.append(
                {
                    "unit_id": unit_id,
                    "status": "invalid" if (root / unit_id).exists() else "pending",
                    "detail": str(error),
                }
            )
        else:
            rows.append({**verified, "detail": "validated immutable decode"})
    return pd.DataFrame(rows)


def collect_decoded_scores(
    units: pd.DataFrame,
    *,
    decode_root: str | Path,
    decoder_bin: str | Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Verify all 202 decodes and assemble both exact Gamma score grids."""

    inventory = validate_unit_inventory(units)
    records = decode_status(
        inventory,
        decode_root=decode_root,
        decoder_bin=decoder_bin,
    )
    if len(records) != EXPECTED_TOTAL_UNITS or set(records["status"]) != {"complete"}:
        raise RuntimeError(
            "Gamma decode bank is incomplete or invalid: "
            + str(records["status"].value_counts().to_dict())
        )
    metadata = inventory.set_index("unit_id")
    rows: list[dict[str, Any]] = []
    statistic_mapping = {
        "hard_mean_call_fraction": "paper_frac_posterior_mean_below",
        "soft_mean_cdf_score": "mean_p_tmrca_lt",
    }
    for record in records.to_dict(orient="records"):
        unit_id = str(record["unit_id"])
        unit = metadata.loc[unit_id]
        reduced = pd.read_csv(record["scores_path"], sep="\t")
        if len(reduced) != 2 * len(gamma_decode.TMRCA_THRESHOLDS_YEARS):
            raise ValueError(f"decode score geometry changed: {unit_id}")
        for source_column, statistic in statistic_mapping.items():
            if source_column not in reduced:
                raise ValueError(f"decode score column is absent: {source_column}")
            for score in reduced.to_dict(orient="records"):
                rows.append(
                    {
                        "unit_id": unit_id,
                        "demography": unit["demography"],
                        "simulation_class": unit["simulation_class"],
                        "selection_coefficient": unit["selection_coefficient"],
                        "seed": unit["seed"],
                        "generation_time_years": unit["generation_time_years"],
                        "statistic": statistic,
                        "window": score["window"],
                        "threshold_years": score["threshold_years"],
                        "score": score[source_column],
                        "final_population_af": unit["final_population_af"],
                        "sample_alt_count": unit["sample_alt_count"],
                        "sample_af": unit["sample_af"],
                    }
                )
    scores = analysis.validate_scores(
        pd.DataFrame(rows, columns=analysis.INPUT_COLUMNS)
    )
    expected_rows = EXPECTED_TOTAL_UNITS * 28
    if len(scores) != expected_rows:
        raise ValueError("assembled score table has the wrong cardinality")
    return scores, records.sort_values("unit_id", kind="mergesort").reset_index(
        drop=True
    )


def _atomic_save_figure_pair(
    figure: Figure, output_dir: Path, stem: str
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "png": output_dir / f"{stem}.png",
        "pdf": output_dir / f"{stem}.pdf",
    }
    for extension, path in outputs.items():
        temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
        try:
            figure.savefig(
                temporary,
                format=extension,
                dpi=300 if extension == "png" else None,
                facecolor="white",
                metadata=(
                    {"Software": "gamma_smc_aou.cosi2_conditioned_gamma_workflow"}
                    if extension == "png"
                    else {
                        "Creator": "gamma_smc_aou.cosi2_conditioned_gamma_workflow",
                        "CreationDate": None,
                        "ModDate": None,
                    }
                ),
            )
            os.replace(temporary, path)
        finally:
            temporary.unlink(missing_ok=True)
    return outputs


def plot_conditioned_score_profiles(
    scores: pd.DataFrame,
    output_dir: str | Path,
    *,
    stem: str = PROFILE_FIGURE_STEM,
) -> dict[str, Path]:
    """Plot selected curves over demography-specific conditioned nulls."""

    frame = analysis.validate_scores(scores)
    destination = Path(output_dir).resolve()
    demographies = ("EAS", "Han")
    statistic_labels = {
        "paper_frac_posterior_mean_below": (
            "Fraction with posterior mean\nTMRCA below threshold"
        ),
        "mean_p_tmrca_lt": "Mean posterior P(TMRCA < threshold)",
    }
    colors = {"local_100kb": "#1764ab", "focal": "#d95f02"}
    x = np.asarray(analysis.TMRCA_THRESHOLDS_YEARS, dtype=float) / 1_000.0
    figure = Figure(figsize=(11.0, 8.5))
    axes = figure.subplots(2, 2, sharex=True, sharey=True, squeeze=False)
    try:
        for row_index, statistic in enumerate(analysis.SUPPORTED_STATISTICS):
            for column_index, demography in enumerate(demographies):
                axis = axes[row_index, column_index]
                cell = frame[
                    (frame["demography"] == demography)
                    & (frame["statistic"] == statistic)
                ]
                for window in analysis.SUPPORTED_WINDOWS:
                    window_frame = cell[cell["window"] == window]
                    profiles = window_frame.pivot(
                        index="unit_id", columns="threshold_years", values="score"
                    ).reindex(columns=analysis.TMRCA_THRESHOLDS_YEARS)
                    classes = (
                        window_frame[["unit_id", "simulation_class"]]
                        .drop_duplicates()
                        .set_index("unit_id")["simulation_class"]
                    )
                    neutral_values = profiles.loc[
                        classes.index[classes == "neutral"]
                    ].to_numpy(float)
                    selected_values = profiles.loc[
                        classes.index[classes == "selected"]
                    ].to_numpy(float)
                    low, median, high = np.quantile(
                        neutral_values, [0.025, 0.5, 0.975], axis=0
                    )
                    color = colors[window]
                    axis.fill_between(x, low, high, color=color, alpha=0.14)
                    axis.plot(x, median, color=color, linestyle="--", linewidth=2.3)
                    axis.plot(
                        x,
                        selected_values[0],
                        color=color,
                        linestyle="-",
                        linewidth=2.8,
                        marker="o" if window == analysis.PRIMARY_WINDOW else "s",
                        markersize=5.5,
                    )
                axis.set_title(
                    f"{demography} - {statistic_labels[statistic]}",
                    fontsize=15,
                    pad=10,
                )
                axis.set_ylim(0, 1)
                axis.set_xticks(x)
                axis.set_xticklabels(["1", "4.5", "10", "20", "30", "40", "50"])
                axis.tick_params(labelsize=12)
                axis.grid(True, color="#d9d9d9", linewidth=0.7)
                axis.text(
                    0.03,
                    0.95,
                    "neutral n=100; selected n=1",
                    transform=axis.transAxes,
                    va="top",
                    fontsize=11.5,
                )
                if column_index == 0:
                    axis.set_ylabel("Unit-level score", fontsize=14)
                if row_index == 1:
                    axis.set_xlabel("TMRCA threshold (kya)", fontsize=14)
        handles: list[Any] = [
            Line2D([0], [0], color=colors["local_100kb"], lw=3, label="+/-100 kb"),
            Line2D([0], [0], color=colors["focal"], lw=3, label="Focal"),
            Line2D([0], [0], color="#333333", lw=2.8, label="Selected"),
            Line2D([0], [0], color="#333333", lw=2.3, ls="--", label="Neutral median"),
            Patch(color="#999999", alpha=0.2, label="Neutral 95% range"),
        ]
        figure.legend(
            handles=handles,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.985),
            ncol=5,
            frameon=False,
            fontsize=12.5,
        )
        figure.suptitle(
            "Gamma-SMC: selected versus conditioned-neutral loci",
            fontsize=20,
            y=0.92,
        )
        figure.subplots_adjust(
            left=0.08,
            right=0.985,
            bottom=0.09,
            top=0.80,
            wspace=0.17,
            hspace=0.34,
        )
        return _atomic_save_figure_pair(figure, destination, stem)
    finally:
        figure.clear()


def plot_conditioning_diagnostics(
    unit_metadata: pd.DataFrame,
    output_dir: str | Path,
    *,
    stem: str = CONDITIONING_FIGURE_STEM,
) -> dict[str, Path]:
    """Plot age, importance, Han gate, and Han endpoint diagnostics."""

    units = validate_unit_inventory(unit_metadata)
    eas = units[
        (units["demography"] == "EAS") & (units["simulation_class"] == "neutral")
    ]
    han = units[
        (units["demography"] == "Han") & (units["simulation_class"] == "neutral")
    ]
    selected = units[units["simulation_class"] == "selected"].set_index("demography")
    figure = Figure(figsize=(11.0, 8.5))
    axes = figure.subplots(2, 2, squeeze=False)
    try:
        ages_kya = eas["allele_age_years"].to_numpy(float) / 1_000.0
        axes[0, 0].hist(ages_kya, bins=18, color="#3182bd", alpha=0.85)
        axes[0, 0].axvline(
            float(selected.loc["EAS", "allele_age_years"]) / 1_000.0,
            color="#cb181d",
            linewidth=3,
            label="Selected EAS age",
        )
        axes[0, 0].set_title("EAS neutral mutation-age draws", fontsize=15)
        axes[0, 0].set_xlabel("Allele age (kya)", fontsize=13)
        axes[0, 0].set_ylabel("Loci", fontsize=13)

        axes[0, 1].scatter(
            eas["birth_generation"],
            eas["eas_selected_log_importance_weight"],
            s=34,
            color="#6a51a3",
            alpha=0.8,
        )
        axes[0, 1].set_title("EAS selected-particle importance weights", fontsize=15)
        axes[0, 1].set_xlabel("Birth generation", fontsize=13)
        axes[0, 1].set_ylabel("Log importance weight", fontsize=13)

        axes[1, 0].hist(
            100 * han["han_generation_645_chb_af"].to_numpy(float),
            bins=14,
            color="#31a354",
            alpha=0.85,
        )
        axes[1, 0].axvspan(
            100 * conditioned.HAN_GATE_LOWER,
            100 * conditioned.HAN_GATE_UPPER,
            color="#006d2c",
            alpha=0.12,
            label="Required gate",
        )
        axes[1, 0].set_title("Han CHB AF at introgression cessation", fontsize=15)
        axes[1, 0].set_xlabel("CHB AF at generation 645 (%)", fontsize=13)
        axes[1, 0].set_ylabel("Loci", fontsize=13)

        axes[1, 1].hist(
            100 * han["final_population_af"].to_numpy(float),
            bins=18,
            color="#fd8d3c",
            alpha=0.85,
        )
        axes[1, 1].axvline(
            100 * float(selected.loc["Han", "final_population_af"]),
            color="#cb181d",
            linewidth=3,
            label="Selected Han endpoint",
        )
        axes[1, 1].set_title("Han neutral present-day CHB AF", fontsize=15)
        axes[1, 1].set_xlabel("Present CHB AF (%)", fontsize=13)
        axes[1, 1].set_ylabel("Loci", fontsize=13)

        for axis in axes.flat:
            axis.tick_params(labelsize=11.5)
            axis.grid(True, color="#dddddd", linewidth=0.7, alpha=0.8)
        handles = [
            Line2D([0], [0], color="#cb181d", lw=3, label="Selected reference"),
            Patch(color="#006d2c", alpha=0.2, label="Han 1.1%-1.3% gate"),
        ]
        figure.legend(
            handles=handles,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.985),
            ncol=2,
            frameon=False,
            fontsize=13,
        )
        figure.suptitle(
            "Conditioned-neutral ascertainment diagnostics", fontsize=20, y=0.92
        )
        figure.subplots_adjust(
            left=0.085,
            right=0.985,
            bottom=0.09,
            top=0.80,
            wspace=0.25,
            hspace=0.38,
        )
        return _atomic_save_figure_pair(figure, Path(output_dir).resolve(), stem)
    finally:
        figure.clear()


def conditioning_contrasts(unit_metadata: pd.DataFrame) -> pd.DataFrame:
    """Return prespecified descriptive age/AF simulation-reference contrasts."""

    units = validate_unit_inventory(unit_metadata)
    rows: list[dict[str, Any]] = []
    for demography in ("EAS", "Han"):
        group = units[units["demography"] == demography]
        neutral = group[group["simulation_class"] == "neutral"]
        selected = group[group["simulation_class"] == "selected"].iloc[0]
        if demography == "EAS":
            observed = float(selected["birth_generation"])
            values = neutral["birth_generation"].to_numpy(float)
            extreme = int(np.count_nonzero(values <= observed))
            statistic = "birth_generation"
            tail = "lower"
            interpretation = "selected_younger_than_endpoint_matched_neutral"
        else:
            observed = float(selected["final_population_af"])
            values = neutral["final_population_af"].to_numpy(float)
            extreme = int(np.count_nonzero(values >= observed))
            statistic = "present_chb_population_af"
            tail = "upper"
            interpretation = "selected_endpoint_higher_than_surviving_neutral"
        rows.append(
            {
                "schema": SCHEMA_VERSION,
                "demography": demography,
                "statistic": statistic,
                "tail": tail,
                "observed_selected": observed,
                "n_neutral": len(values),
                "neutral_extreme_count": extreme,
                "plus_one_simulation_p": (1 + extreme) / (1 + len(values)),
                "neutral_mean": float(np.mean(values)),
                "neutral_median": float(np.median(values)),
                "interpretation": interpretation,
                "inference_scope": "descriptive_simulation_reference_diagnostic",
            }
        )
    return pd.DataFrame(rows)


def _derive_inference(scores: pd.DataFrame) -> dict[str, pd.DataFrame]:
    derived = analysis.analyze_scores(scores)
    pointwise = derived["pointwise_pvalues"].copy()
    pointwise.loc[pointwise["demography"] == "EAS", "neutral_bank"] = (
        "same_demography_exact_present_population_count_"
        "importance_resampled_neutral_mutation_age"
    )
    pointwise.loc[pointwise["demography"] == "Han", "neutral_bank"] = (
        "same_demography_introgressed_g645_af_gate_present_segregating_neutral"
    )
    return {
        "scores": derived["scores"],
        "pointwise_pvalues": pointwise,
        "minp_omnibus": derived["minp_omnibus"],
    }


def _eas_diagnostic_table(units: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "unit_id",
        "seed",
        "birth_generation",
        "allele_age_years",
        "final_population_af",
        "eas_source_particle_id",
        "eas_draw_occurrence",
        "eas_particle_multiplicity",
        "eas_selected_log_importance_weight",
        "eas_global_effective_sample_size",
        "eas_global_weight_sum_over_max",
        "eas_global_max_normalized_weight",
        "eas_time_reversal_is_approximation",
        "source_completion_sha256",
    ]
    return (
        units[(units["demography"] == "EAS") & (units["simulation_class"] == "neutral")]
        .loc[:, columns]
        .reset_index(drop=True)
    )


def _han_distribution_table(units: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "unit_id",
        "seed",
        "han_candidate_id",
        "han_candidate_order",
        "han_candidate_source_kind",
        "han_generation_645_chb_af",
        "final_population_af",
        "sample_alt_count",
        "sample_af",
        "source_completion_sha256",
    ]
    return (
        units[(units["demography"] == "Han") & (units["simulation_class"] == "neutral")]
        .loc[:, columns]
        .reset_index(drop=True)
    )


def _validate_figure(path: Path, extension: str) -> None:
    if path.is_symlink() or not path.is_file() or path.stat().st_size < 8:
        raise ValueError(f"result figure is missing or empty: {path}")
    header = path.read_bytes()[:8]
    if extension == "png" and header != b"\x89PNG\r\n\x1a\n":
        raise ValueError("result PNG signature is invalid")
    if extension == "pdf" and not header.startswith(b"%PDF-"):
        raise ValueError("result PDF signature is invalid")


def _result_contract(
    repo_root: Path,
    units: pd.DataFrame,
    scores: pd.DataFrame,
    decode_records: pd.DataFrame,
    *,
    decoder_bin: Path,
    cosi2_binary: Path,
) -> dict[str, Any]:
    neutral = units[units["simulation_class"] == "neutral"]
    selected = units[units["simulation_class"] == "selected"]
    return {
        "schema": SCHEMA_VERSION,
        "scientific_contract": {
            "units": EXPECTED_TOTAL_UNITS,
            "neutral_per_demography": EXPECTED_NEUTRAL_PER_DEMOGRAPHY,
            "selected_per_demography": EXPECTED_SELECTED_PER_DEMOGRAPHY,
            "selection_coefficient": SELECTION_COEFFICIENT,
            "eas_conditioning": (
                "exact_population_count_14704_of_75274_global_importance_"
                "categorical_resampling_with_replacement"
            ),
            "eas_reverse_time_is_approximation": True,
            "han_conditioning": (
                "neanderthal_birth_g2400_chb_af_0p011_to_0p013_g645_"
                "strict_present_segregation"
            ),
            "han_final_af_matched": False,
            "sample_af_conditioned": False,
            "pair_panel": {
                "scheme": ("fixed_sha256_ranked_100_haplotypes_all_unordered_pairs"),
                "seed": gamma_decode.DEFAULT_PAIR_PANEL_SEED,
                "haplotypes": gamma_decode.DEFAULT_PAIR_PANEL_HAPLOTYPES,
                "unordered_pairs": gamma_decode.DEFAULT_PAIR_COUNT,
            },
            "thresholds_years": list(gamma_decode.TMRCA_THRESHOLDS_YEARS),
            "primary_window": "local_100kb",
        },
        "input_payloads": {
            "unit_metadata_sha256": _portable_frame_sha256(units),
            "conditioned_metadata_sha256": _portable_frame_sha256(neutral),
            "scores_sha256": _portable_frame_sha256(scores),
            "selected_run_completion_sha256": {
                str(row["demography"]): str(row["source_completion_sha256"])
                for row in selected.to_dict(orient="records")
            },
            "neutral_completion_sha256": {
                str(row["unit_id"]): str(row["source_completion_sha256"])
                for row in neutral.to_dict(orient="records")
            },
            "decode_completion_sha256": {
                str(row["unit_id"]): str(row["decode_completion_sha256"])
                for row in decode_records.to_dict(orient="records")
            },
        },
        "binaries": {
            "cosi2": {
                "path": _relative(cosi2_binary, repo_root),
                "sha256": sha256_file(cosi2_binary),
            },
            "gamma_smc": {
                "path": _relative(decoder_bin, repo_root),
                "sha256": sha256_file(decoder_bin),
            },
        },
        "implementation": {
            label: {
                "path": _relative(path, repo_root),
                "sha256": sha256_file(path),
            }
            for label, path in {
                "workflow": MODULE_PATH,
                "conditioned_null": Path(conditioned.__file__).resolve(),
                "gamma_decode": Path(gamma_decode.__file__).resolve(),
                "gamma_analysis": Path(analysis.__file__).resolve(),
            }.items()
        },
        "runtime_versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
        },
    }


def _run_results_text(
    pointwise: pd.DataFrame,
    omnibus: pd.DataFrame,
    contrasts: pd.DataFrame,
) -> str:
    primary = pointwise[
        pointwise["primary_window"]
        & (pointwise["statistic"] == "paper_frac_posterior_mean_below")
    ].sort_values(["demography", "threshold_years"])
    primary_omnibus = omnibus[
        omnibus["primary_window"]
        & (omnibus["statistic"] == "paper_frac_posterior_mean_below")
    ].sort_values("demography")
    lines = [
        "# CoSi2 Gamma-SMC conditioned-neutral results",
        "",
        "EAS neutral loci match the selected present population count exactly, ",
        "but their mutation ages come from a fixed global importance-resampled ",
        "reverse-time approximation. Han neutral loci match the 1.1%-1.3% CHB ",
        "gate at generation 645 and survive to the present, without matching the ",
        "selected present-day AF. Neither null conditions on sampled AF/carriers.",
        "",
        "## Primary Gamma-SMC pointwise results",
        "",
        "| Demography | Threshold years | Selected score | Upper p | Two-sided p | Descriptive AUC |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in primary.to_dict(orient="records"):
        lines.append(
            "| {demography} | {threshold_years} | {observed_score:.6g} | "
            "{mc_p_upper:.6g} | {mc_p_two_sided:.6g} | "
            "{neutral_midrank_descriptive_auc:.6g} |".format(**row)
        )
    lines.extend(
        [
            "",
            "## Correlation-aware seven-threshold minP",
            "",
            "| Demography | Upper omnibus p | Two-sided omnibus p |",
            "|---|---:|---:|",
        ]
    )
    for row in primary_omnibus.to_dict(orient="records"):
        lines.append(
            "| {demography} | {omnibus_mc_p_upper:.6g} | "
            "{omnibus_mc_p_two_sided:.6g} |".format(**row)
        )
    lines.extend(
        [
            "",
            "## Conditioning diagnostics",
            "",
            "| Demography | Statistic | Tail | Selected | Extreme / 100 | Plus-one p |",
            "|---|---|---|---:|---:|---:|",
        ]
    )
    for row in contrasts.to_dict(orient="records"):
        lines.append(
            "| {demography} | {statistic} | {tail} | {observed_selected:.6g} | "
            "{neutral_extreme_count} | {plus_one_simulation_p:.6g} |".format(**row)
        )
    lines.extend(
        [
            "",
            "With one selected locus per demography, AUC is its neutral-bank ",
            "midrank, not a stable ROC power estimate. Gamma TMRCA is a genealogy ",
            "quantity and is not the same as the EAS mutation-age diagnostic.",
            "",
        ]
    )
    return "\n".join(lines)


def write_results_bundle(
    repo_root: str | Path,
    *,
    results_dir: str | Path,
    units: pd.DataFrame,
    scores: pd.DataFrame,
    decode_records: pd.DataFrame,
    decoder_bin: str | Path,
    cosi2_binary: str | Path,
) -> dict[str, Any]:
    """Create or strictly reuse one atomic results directory."""

    root = Path(repo_root).resolve()
    destination = Path(results_dir).resolve()
    inventory = validate_unit_inventory(units)
    score_frame = analysis.validate_scores(scores)
    decoder = Path(decoder_bin).resolve()
    cosi2 = Path(cosi2_binary).resolve()
    decode_counts = decode_records[
        [
            "unit_id",
            "decode_full_sample_alt_count",
            "decode_pair_panel_alt_count",
            "decode_focal_count_status",
        ]
    ]
    if (
        len(decode_counts) != EXPECTED_TOTAL_UNITS
        or decode_counts["unit_id"].duplicated().any()
    ):
        raise ValueError("decode panel-count inventory is incomplete")
    enriched = inventory.merge(
        decode_counts, on="unit_id", how="left", validate="one_to_one"
    )
    contract = _result_contract(
        root,
        inventory,
        score_frame,
        decode_records,
        decoder_bin=decoder,
        cosi2_binary=cosi2,
    )
    if destination.exists() and any(destination.iterdir()):
        return verify_results_bundle(
            root, results_dir=destination, expected_contract=contract
        )
    if destination.exists() and not destination.is_dir():
        raise ValueError("results destination exists but is not a directory")
    if destination.exists():
        destination.rmdir()
    destination.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(
        tempfile.mkdtemp(prefix=f".{destination.name}.staging.", dir=destination.parent)
    )
    try:
        derived = _derive_inference(score_frame)
        pointwise = derived["pointwise_pvalues"]
        omnibus = derived["minp_omnibus"]
        conditioned_metadata = enriched[
            enriched["simulation_class"] == "neutral"
        ].reset_index(drop=True)
        eas_diagnostics = _eas_diagnostic_table(inventory)
        han_distribution = _han_distribution_table(inventory)
        contrasts = conditioning_contrasts(inventory)
        frames = {
            "scores": score_frame,
            "conditioned_metadata": conditioned_metadata,
            "unit_metadata": enriched,
            "eas_diagnostics": eas_diagnostics,
            "han_distribution": han_distribution,
            "conditioning_contrasts": contrasts,
            "pointwise_pvalues": pointwise,
            "minp_omnibus": omnibus,
        }
        for label, frame in frames.items():
            _atomic_frame(stage / RESULT_OUTPUTS[label], frame)
        _atomic_text(
            stage / RESULT_OUTPUTS["run_results"],
            _run_results_text(pointwise, omnibus, contrasts),
        )
        profile = plot_conditioned_score_profiles(score_frame, stage)
        conditioning_plot = plot_conditioning_diagnostics(inventory, stage)
        expected_figures = {
            "profile_png": profile["png"],
            "profile_pdf": profile["pdf"],
            "conditioning_png": conditioning_plot["png"],
            "conditioning_pdf": conditioning_plot["pdf"],
        }
        for label, path in expected_figures.items():
            expected = stage / RESULT_OUTPUTS[label]
            if Path(path) != expected:
                raise ValueError("analysis figure helper returned unexpected path")
            _validate_figure(expected, "png" if label.endswith("png") else "pdf")
        outputs = {
            label: _output_record(
                stage / filename,
                stage,
                rows=len(frames[label]) if label in frames else None,
            )
            for label, filename in RESULT_OUTPUTS.items()
        }
        completion = {
            "schema": RESULTS_SCHEMA,
            "status": "complete",
            "contract": contract,
            "contract_sha256": _canonical_sha256(contract),
            "outputs": outputs,
            "output_count": len(outputs),
        }
        _atomic_json(stage / RESULTS_COMPLETION_FILENAME, completion)
        if {path.name for path in stage.iterdir()} != {
            RESULTS_COMPLETION_FILENAME,
            *RESULT_OUTPUTS.values(),
        }:
            raise ValueError("staged results have an unexpected inventory")
        os.replace(stage, destination)
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    return {
        **_read_json(destination / RESULTS_COMPLETION_FILENAME, "results completion"),
        "cache_hit": False,
    }


def verify_results_bundle(
    repo_root: str | Path,
    *,
    results_dir: str | Path,
    expected_contract: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Strictly verify the result inventory, checksums, and recomputable tables."""

    root = Path(repo_root).resolve()
    destination = Path(results_dir).resolve()
    if destination.is_symlink() or not destination.is_dir():
        raise ValueError(f"results directory is absent: {destination}")
    expected_names = {RESULTS_COMPLETION_FILENAME, *RESULT_OUTPUTS.values()}
    children = list(destination.iterdir())
    if (
        any(path.is_symlink() or path.is_dir() for path in children)
        or {path.name for path in children} != expected_names
    ):
        raise ValueError("results inventory mismatch")
    completion = _read_json(
        destination / RESULTS_COMPLETION_FILENAME, "results completion"
    )
    contract = completion.get("contract")
    if (
        completion.get("schema") != RESULTS_SCHEMA
        or completion.get("status") != "complete"
        or not isinstance(contract, dict)
        or completion.get("contract_sha256") != _canonical_sha256(contract)
    ):
        raise ValueError("results completion contract is invalid")
    if expected_contract is not None and _canonical_sha256(
        contract
    ) != _canonical_sha256(expected_contract):
        raise ValueError("results completion is stale for current inputs")
    expected_implementation = {
        "workflow": MODULE_PATH,
        "conditioned_null": Path(conditioned.__file__).resolve(),
        "gamma_decode": Path(gamma_decode.__file__).resolve(),
        "gamma_analysis": Path(analysis.__file__).resolve(),
    }
    for label, path in expected_implementation.items():
        if contract.get("implementation", {}).get(label, {}).get(
            "sha256"
        ) != sha256_file(path):
            raise ValueError(f"results implementation binding is stale: {label}")
    expected_runtime = {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "matplotlib": matplotlib.__version__,
    }
    if contract.get("runtime_versions") != expected_runtime:
        raise ValueError("results runtime-version binding is stale")
    for label, record in contract.get("binaries", {}).items():
        recorded = Path(str(record.get("path", "")))
        binary = recorded if recorded.is_absolute() else root / recorded
        if not binary.is_file() or sha256_file(binary) != record.get("sha256"):
            raise ValueError(f"results binary binding is stale: {label}")
    outputs = completion.get("outputs")
    if (
        not isinstance(outputs, dict)
        or set(outputs) != set(RESULT_OUTPUTS)
        or completion.get("output_count") != len(RESULT_OUTPUTS)
    ):
        raise ValueError("results output manifest changed")
    for label, filename in RESULT_OUTPUTS.items():
        path = destination / filename
        record = outputs[label]
        if (
            record.get("path") != filename
            or path.stat().st_size != int(record.get("size_bytes", -1))
            or sha256_file(path) != record.get("sha256")
        ):
            raise ValueError(f"results output checksum failed: {label}")
    scores = analysis.validate_scores(
        pd.read_csv(destination / RESULT_OUTPUTS["scores"], sep="\t")
    )
    persisted_units = pd.read_csv(
        destination / RESULT_OUTPUTS["unit_metadata"], sep="\t"
    )
    units = validate_unit_inventory(persisted_units)
    derived = _derive_inference(scores)
    expected_frames = {
        "pointwise_pvalues": derived["pointwise_pvalues"],
        "minp_omnibus": derived["minp_omnibus"],
        "eas_diagnostics": _eas_diagnostic_table(units),
        "han_distribution": _han_distribution_table(units),
        "conditioning_contrasts": conditioning_contrasts(units),
    }
    for label, expected in expected_frames.items():
        observed = pd.read_csv(destination / RESULT_OUTPUTS[label], sep="\t")
        pd.testing.assert_frame_equal(
            observed,
            expected,
            check_dtype=False,
            check_exact=False,
            rtol=0,
            atol=1e-12,
        )
        if int(outputs[label].get("rows", -1)) != len(expected):
            raise ValueError(f"results row count changed: {label}")
    for label in ("profile_png", "conditioning_png"):
        _validate_figure(destination / RESULT_OUTPUTS[label], "png")
    for label in ("profile_pdf", "conditioning_pdf"):
        _validate_figure(destination / RESULT_OUTPUTS[label], "pdf")
    if (
        not (destination / RESULT_OUTPUTS["run_results"])
        .read_text(encoding="utf-8")
        .strip()
    ):
        raise ValueError("RUN_RESULTS.md is empty")
    return {**completion, "cache_hit": True}
