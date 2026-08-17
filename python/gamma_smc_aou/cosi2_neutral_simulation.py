"""Restartable natural-neutral CoSi2 simulations for the EAS and Han models.

The neutral null is deliberately *not* frequency matched to the selected
replicate.  It keeps the selected model's allele origin and birth age, sets
``s=0``, and conditions only on a positive population frequency at generation
zero.  The Han allele is born in Neanderthal, so survival in CHB also requires
introgression through the native 5R19 continuous-migration history; there is no
generation-645 CHB-frequency gate.

This module binds read-only to :mod:`gamma_smc_aou.cosi2_framework`.  It never
modifies the selected framework or selected simulation artifacts.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import platform
import re
import shutil
import subprocess
import tempfile
import time
import uuid
from collections.abc import Iterable, Mapping, Sequence
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from . import cosi2_framework as framework

PLAN_SCHEMA = "gamma-smc.cosi2-natural-neutral-plan/v1"
PLAN_COMPLETION_SCHEMA = "gamma-smc.cosi2-natural-neutral-plan-completion/v1"
UNIT_COMPLETION_SCHEMA = "gamma-smc.cosi2-natural-neutral-unit-completion/v1"
UNIT_FAILURE_SCHEMA = "gamma-smc.cosi2-natural-neutral-unit-failure/v1"

MODULE_PATH = "python/gamma_smc_aou/cosi2_neutral_simulation.py"
DEFAULT_PLAN_RELATIVE = "focused_selection_EAS_sim/cosi2_natural_neutral_plan"
DEFAULT_RUN_RELATIVE = "focused_selection_EAS_sim/cosi2_natural_neutral_runs"

EAS_MODEL_ID = "eas_neutral_survival"
HAN_MODEL_ID = "han_introgressed_neutral_survival"
MODEL_IDS = (EAS_MODEL_ID, HAN_MODEL_ID)

REPLICATES_PER_MODEL = 100
TOTAL_UNITS = 2 * REPLICATES_PER_MODEL
EAS_SEED_START = 2_026_081_600
HAN_SEED_START = 2_026_082_600
SEED_FORMULAS = {
    EAS_MODEL_ID: f"seed={EAS_SEED_START}+replicate_index",
    HAN_MODEL_ID: f"seed={HAN_SEED_START}+replicate_index",
}
SELECTED_SEEDS_EXCLUDED = (20_240_524, 20_240_658)

NEUTRAL_SELECTION_COEFFICIENT = 0.0
SURVIVAL_ENDPOINT_LOWER = 1e-12
SURVIVAL_ENDPOINT_UPPER = 1.0
# CoSi2 v2.4 ValRange splits on every hyphen, so the mathematically equivalent
# literal ``1e-12-1`` is not parseable.  Decimal notation is required here.
SURVIVAL_ENDPOINT_TOKEN = "0.000000000001-1"
SOURCE_SELECTED_COEFFICIENT = 0.01

DEFAULT_MAX_WORKERS = 20
MAX_WORKERS = 20
DEFAULT_COSI_MAXATTEMPTS = 1_000_000
SAMPLE_HAPLOIDS = framework.DEFAULT_SAMPLE_HAPLOIDS

PLAN_FILENAME = "neutral_units.tsv"
MANIFEST_FILENAME = "neutral_plan.json"
PLAN_COMPLETION_FILENAME = "neutral_plan_completion.json"
RECOMBINATION_MAP_FILENAME = "recombination_map.tsv"
UNIT_COMPLETION_FILENAME = "completion.json"


@dataclass(frozen=True)
class NeutralSimulationPlan:
    """The immutable scientific and seed contract for the natural-neutral bank."""

    replicates_per_model: int = REPLICATES_PER_MODEL
    eas_seed_start: int = EAS_SEED_START
    han_seed_start: int = HAN_SEED_START
    endpoint_lower: float = SURVIVAL_ENDPOINT_LOWER
    endpoint_upper: float = SURVIVAL_ENDPOINT_UPPER
    selection_coefficient: float = NEUTRAL_SELECTION_COEFFICIENT
    source_selected_coefficient: float = SOURCE_SELECTED_COEFFICIENT
    sample_haploids: int = SAMPLE_HAPLOIDS

    def validate(self) -> None:
        """Reject any drift from the predeclared 100 + 100 null contract."""

        expected = NeutralSimulationPlan()
        if self != expected:
            raise ValueError(
                "the natural-neutral plan is immutable: use the canonical 100 + 100 "
                "unit, seed, endpoint, and age-matching contract"
            )
        seeds = [
            *(
                self.eas_seed_start + index
                for index in range(self.replicates_per_model)
            ),
            *(
                self.han_seed_start + index
                for index in range(self.replicates_per_model)
            ),
        ]
        if len(seeds) != TOTAL_UNITS or len(set(seeds)) != TOTAL_UNITS:
            raise RuntimeError("canonical neutral seed streams overlap")
        if set(seeds).intersection(SELECTED_SEEDS_EXCLUDED):
            raise RuntimeError("canonical neutral seeds overlap a selected replicate")
        if any(seed <= 0 or seed >= 2**32 for seed in seeds):
            raise RuntimeError("canonical neutral seeds must be unique uint32 values")

    def to_record(self) -> dict[str, Any]:
        self.validate()
        return asdict(self)


def _utc_now() -> str:
    return datetime.now(UTC).isoformat().replace("+00:00", "Z")


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


def _sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _atomic_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}.{uuid.uuid4().hex}")
    try:
        temporary.write_text(value, encoding="utf-8", newline="\n")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_json(path: Path, value: Mapping[str, Any]) -> None:
    _atomic_text(
        path,
        json.dumps(value, sort_keys=True, indent=2, ensure_ascii=False, allow_nan=False)
        + "\n",
    )


def _frame_text(frame: pd.DataFrame) -> str:
    return frame.to_csv(sep="\t", index=False, lineterminator="\n")


def _atomic_frame(path: Path, frame: pd.DataFrame) -> None:
    _atomic_text(path, _frame_text(frame))


def _regular_file(path: Path, label: str) -> Path:
    if path.is_symlink():
        raise ValueError(f"{label} must not be a symlink: {path}")
    resolved = path.resolve()
    if not resolved.is_file():
        raise ValueError(f"{label} must be a regular, non-symlink file: {path}")
    return resolved


def _regular_directory(path: Path, label: str) -> Path:
    if path.is_symlink():
        raise ValueError(f"{label} must not be a symlink: {path}")
    resolved = path.resolve()
    if not resolved.is_dir():
        raise ValueError(f"{label} must be a regular directory: {path}")
    return resolved


def _json_file(path: Path, *, label: str) -> dict[str, Any]:
    _regular_file(path, label)
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"{label} is not valid UTF-8 JSON: {path}") from error
    if not isinstance(value, dict):
        raise TypeError(f"{label} must contain a JSON object")
    return value


def _portable_record(record: Mapping[str, Any]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in record.items():
        if isinstance(value, np.generic):
            value = value.item()
        if pd.isna(value):
            value = None
        result[str(key)] = value
    return result


def _selected_model_cells(
    repo_root: Path,
) -> tuple[framework.FrameworkPlan, pd.DataFrame]:
    source_plan = framework.FrameworkPlan(
        selection_coefficients=(SOURCE_SELECTED_COEFFICIENT,)
    )
    source_plan.validate()
    history = framework.load_eas_history(repo_root)
    cells = framework.build_model_cells(history, source_plan)
    if set(cells["cell_id"].astype(str)) != {
        "eas_s0p01",
        "han_introgressed_s0p01",
    }:
        raise RuntimeError("the bound s=0.01 framework cells changed unexpectedly")
    return source_plan, cells


def build_neutral_model_cells(repo_root: str | Path) -> pd.DataFrame:
    """Return the two age/origin-matched model definitions used by all units."""

    root = Path(repo_root).resolve()
    _source_plan, selected = _selected_model_cells(root)
    rows: list[dict[str, Any]] = []
    for selected_row in selected.to_dict(orient="records"):
        sample_population = str(selected_row["sample_population"])
        model_id = EAS_MODEL_ID if sample_population == "EAS" else HAN_MODEL_ID
        row = dict(selected_row)
        row.update(
            {
                "cell_id": model_id,
                "model_id": model_id,
                "source_selected_cell_id": selected_row["cell_id"],
                "source_selected_coefficient": SOURCE_SELECTED_COEFFICIENT,
                "selection_coefficient": NEUTRAL_SELECTION_COEFFICIENT,
                "present_af_target": None,
                "present_af_lower": SURVIVAL_ENDPOINT_LOWER,
                "present_af_upper": SURVIVAL_ENDPOINT_UPPER,
                "endpoint_conditioning": "population_survival_only",
                "final_population_af_conditioned": False,
                "sample_af_conditioned": False,
                "sample_carrier_conditioned": False,
                "fixation_allowed": True,
                "adaptive_stopping": False,
                "seed": 0,
            }
        )
        if sample_population == "Han":
            row.update(
                {
                    "introgressed_allele_af_target": None,
                    "introgressed_allele_af_lower": None,
                    "introgressed_allele_af_upper": None,
                    "introgression_enforcement": (
                        "Neanderthal_birth_plus_native_5R19_migration_plus_CHB_survival"
                    ),
                    "han_intermediate_af_conditioned": False,
                    "han_migration_end_af_role": "descriptive_only",
                }
            )
        else:
            row.update(
                {
                    "han_intermediate_af_conditioned": False,
                    "han_migration_end_af_role": "not_applicable",
                }
            )
        rows.append(row)
    result = pd.DataFrame(rows).sort_values("model_id", kind="mergesort")
    if len(result) != 2:
        raise RuntimeError("natural-neutral model cardinality changed")
    return result.reset_index(drop=True)


def _replace_parameter_contract(
    base_text: str, *, cell: Mapping[str, Any], sweep_line: str
) -> str:
    lines: list[str] = []
    saw_sweep = False
    saw_seed = False
    for line in base_text.splitlines():
        if line.startswith("# Generated CoSi2"):
            lines.extend(
                [
                    line.replace(" model;", " natural-neutral model;"),
                    "# Null: same allele origin and birth age as the s=0.01 selected case.",
                    "# Ascertainment: present population survival only; fixation is allowed.",
                    "# No final target AF, sampled AF, or sampled-carrier conditioning.",
                ]
            )
        elif line.startswith("# The operational CHB allele-AF gate"):
            lines.append(
                "# No generation-645 CHB AF gate; that frequency is descriptive only."
            )
        elif line.startswith("pop_event sweep_mult_standing "):
            if saw_sweep:
                raise RuntimeError(
                    "base parameter text contains duplicate sweep directives"
                )
            lines.append(sweep_line)
            saw_sweep = True
        elif line.startswith("random_seed "):
            lines.extend(
                [
                    "# Per-unit fixed seed is mandatory via the runner's -r option.",
                    "random_seed 0",
                ]
            )
            saw_seed = True
        else:
            lines.append(line)
    if not saw_sweep or not saw_seed:
        raise RuntimeError("base parameter text lacks the sweep or seed directive")
    text = "\n".join(lines) + "\n"
    lint_neutral_parameter_text(text, expected_cell=cell)
    return text


def render_neutral_parameter_files(repo_root: str | Path) -> dict[str, str]:
    """Render the two parser-safe ``s=0`` survival-only CoSi2 configs."""

    root = Path(repo_root).resolve()
    source_plan, _selected = _selected_model_cells(root)
    history = framework.load_eas_history(root)
    rle = framework.rle_eas_history(history)
    han_events = framework.build_han_event_ledger(source_plan)
    texts: dict[str, str] = {}
    for cell in build_neutral_model_cells(root).to_dict(orient="records"):
        model_id = str(cell["model_id"])
        if model_id == EAS_MODEL_ID:
            base = framework.render_eas_parameter_file(rle, cell, source_plan)
            sweep = (
                f'pop_event sweep_mult_standing "{model_id}" 1 '
                f"{int(cell['mutation_birth_generations_ago'])} 0 "
                f"{framework.FOCAL_RELATIVE_POSITION:g} {SURVIVAL_ENDPOINT_TOKEN} 1 "
                f"{int(cell['selection_onset_generations_ago'])}"
            )
        elif model_id == HAN_MODEL_ID:
            base = framework.render_han_parameter_file(han_events, cell, source_plan)
            sweep = (
                f'pop_event sweep_mult_standing "{model_id}" 4 '
                f"{int(cell['mutation_birth_generations_ago'])} 0 "
                f"{framework.FOCAL_RELATIVE_POSITION:g} {SURVIVAL_ENDPOINT_TOKEN} 3 "
                f"{int(cell['selection_onset_generations_ago'])}"
            )
        else:  # pragma: no cover - protected by build_neutral_model_cells
            raise RuntimeError(f"unsupported neutral model: {model_id}")
        texts[model_id] = _replace_parameter_contract(base, cell=cell, sweep_line=sweep)
    return texts


def lint_neutral_parameter_text(
    text: str, *, expected_cell: Mapping[str, Any]
) -> dict[str, Any]:
    """Lint native syntax and fail closed on any frequency-matching condition."""

    audit = framework.lint_parameter_text(text, expected_cell=expected_cell)
    if audit["endpoint_token"] != SURVIVAL_ENDPOINT_TOKEN:
        raise ValueError("neutral config does not use the survival-only endpoint")
    sweep_lines = [
        line.strip()
        for line in text.splitlines()
        if line.strip().startswith("pop_event sweep_mult_standing ")
    ]
    tokens = [
        token
        for token in next(csv.reader([sweep_lines[0]], delimiter=" ", quotechar='"'))
        if token
    ]
    if float(tokens[5]) != 0.0:
        raise ValueError("natural-neutral config must set s=0")
    forbidden = ("0.195-0.205", "0.011-0.013", "applied post hoc")
    if any(token in text for token in forbidden):
        raise ValueError("neutral config retains a selected or intermediate AF gate")
    if "random_seed 0" not in text:
        raise ValueError("neutral template must require a per-unit CLI seed override")
    audit.update(
        {
            "selection_coefficient": 0.0,
            "survival_only": True,
            "fixation_allowed": True,
            "sample_af_conditioned": False,
            "han_intermediate_af_conditioned": False,
        }
    )
    return audit


def _config_path(model_id: str) -> str:
    return f"configs/{model_id}.par"


def build_neutral_simulation_plan(repo_root: str | Path) -> pd.DataFrame:
    """Build all 200 predeclared units; outcomes never alter this ledger."""

    root = Path(repo_root).resolve()
    plan = NeutralSimulationPlan()
    plan.validate()
    cells = build_neutral_model_cells(root).set_index("model_id")
    configs = render_neutral_parameter_files(root)
    config_hashes = {key: _sha256_text(value) for key, value in configs.items()}
    rows: list[dict[str, Any]] = []
    for model_order, model_id in enumerate(MODEL_IDS):
        cell = cells.loc[model_id]
        seed_start = EAS_SEED_START if model_id == EAS_MODEL_ID else HAN_SEED_START
        for replicate_index in range(REPLICATES_PER_MODEL):
            seed = seed_start + replicate_index
            rows.append(
                {
                    "plan_order": model_order * REPLICATES_PER_MODEL + replicate_index,
                    "unit_id": f"{model_id}_r{replicate_index:03d}",
                    "model_id": model_id,
                    "replicate_index": replicate_index,
                    "seed": seed,
                    "seed_formula": SEED_FORMULAS[model_id],
                    "config_path": _config_path(model_id),
                    "config_sha256": config_hashes[model_id],
                    "demography_id": cell["demography_id"],
                    "sample_population": cell["sample_population"],
                    "sample_haploids": SAMPLE_HAPLOIDS,
                    "birth_population": cell["birth_population"],
                    "selection_population": cell["selection_population"],
                    "mutation_birth_generations_ago": int(
                        cell["mutation_birth_generations_ago"]
                    ),
                    "selection_onset_generations_ago": int(
                        cell["selection_onset_generations_ago"]
                    ),
                    "selection_coefficient": 0.0,
                    "endpoint_lower": SURVIVAL_ENDPOINT_LOWER,
                    "endpoint_upper": SURVIVAL_ENDPOINT_UPPER,
                    "endpoint_conditioning": "population_survival_only",
                    "final_population_af_conditioned": False,
                    "sample_af_conditioned": False,
                    "sample_carrier_conditioned": False,
                    "han_intermediate_af_conditioned": False,
                    "fixation_allowed": True,
                    "adaptive_stopping": False,
                    "source_selected_cell_id": cell["source_selected_cell_id"],
                    "source_selected_coefficient": SOURCE_SELECTED_COEFFICIENT,
                }
            )
    frame = pd.DataFrame(rows).sort_values("plan_order", kind="mergesort")
    if len(frame) != TOTAL_UNITS or frame["unit_id"].duplicated().any():
        raise RuntimeError("natural-neutral plan has the wrong unit inventory")
    if frame["seed"].duplicated().any() or set(frame["seed"]).intersection(
        SELECTED_SEEDS_EXCLUDED
    ):
        raise RuntimeError("natural-neutral plan has a seed collision")
    return frame.reset_index(drop=True)


def _source_binding(repo_root: Path) -> dict[str, Any]:
    paths = (
        framework.MODULE_PATH,
        MODULE_PATH,
        framework.EAS_RESOURCE_PATH,
    )
    return {
        "cosi2_version": framework.COSI2_VERSION,
        "cosi2_commit": framework.COSI2_COMMIT,
        "framework_plan": NeutralSimulationPlan().to_record(),
        "source_selected_coefficient": SOURCE_SELECTED_COEFFICIENT,
        "files": {relative: sha256_file(repo_root / relative) for relative in paths},
    }


def _inventory(root: Path, *, exclude: Iterable[str] = ()) -> dict[str, Any]:
    excluded = set(exclude)
    records: dict[str, Any] = {}
    items = sorted(root.rglob("*"))
    symlinks = [
        item.relative_to(root).as_posix() for item in items if item.is_symlink()
    ]
    if symlinks:
        raise ValueError(
            "artifact inventory contains symlinks: " + ", ".join(symlinks[:10])
        )
    for path in (item for item in items if item.is_file()):
        relative = path.relative_to(root).as_posix()
        if relative in excluded:
            continue
        records[relative] = {
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
    return records


def _plan_manifest(
    repo_root: Path,
    *,
    units: pd.DataFrame,
    configs: Mapping[str, str],
    plan_tsv_sha256: str,
) -> dict[str, Any]:
    return {
        "schema": PLAN_SCHEMA,
        "status": "planned",
        "unit_count": TOTAL_UNITS,
        "units_per_model": REPLICATES_PER_MODEL,
        "model_ids": list(MODEL_IDS),
        "plan_tsv": PLAN_FILENAME,
        "plan_tsv_sha256": plan_tsv_sha256,
        "plan_order_is_immutable": True,
        "all_predeclared_units_are_retained": True,
        "adaptive_stopping": False,
        "seed_contract": {
            "formulas": SEED_FORMULAS,
            "eas_range_inclusive": [EAS_SEED_START, EAS_SEED_START + 99],
            "han_range_inclusive": [HAN_SEED_START, HAN_SEED_START + 99],
            "selected_seeds_excluded": list(SELECTED_SEEDS_EXCLUDED),
            "cosi2_effective_uint32_uniqueness_checked": True,
        },
        "neutral_semantics": {
            "selection_coefficient": 0.0,
            "source_selected_coefficient": SOURCE_SELECTED_COEFFICIENT,
            "source_selected_origin_and_birth_age_retained": True,
            "endpoint_population_af_interval": [
                SURVIVAL_ENDPOINT_LOWER,
                SURVIVAL_ENDPOINT_UPPER,
            ],
            "endpoint_parser_token": SURVIVAL_ENDPOINT_TOKEN,
            "endpoint_conditioning": "survival_only",
            "final_population_af_conditioned": False,
            "sample_af_conditioned": False,
            "sample_carrier_conditioned": False,
            "han_intermediate_af_conditioned": False,
            "fixation_allowed": True,
        },
        "configs": {
            model_id: {
                "path": _config_path(model_id),
                "sha256": _sha256_text(configs[model_id]),
            }
            for model_id in MODEL_IDS
        },
        "source_binding": _source_binding(repo_root),
        "unit_record_sha256": {
            str(row["unit_id"]): _canonical_sha256(_portable_record(row))
            for row in units.to_dict(orient="records")
        },
    }


def generate_neutral_plan(
    repo_root: str | Path, *, output_dir: str | Path | None = None
) -> dict[str, Any]:
    """Atomically create, or strictly verify, the immutable neutral plan bundle."""

    root = Path(repo_root).resolve()
    destination = (
        Path(output_dir).resolve()
        if output_dir is not None
        else root / DEFAULT_PLAN_RELATIVE
    )
    if destination.exists():
        try:
            result = verify_neutral_plan(root, plan_dir=destination)
        except Exception as error:
            if destination.is_dir() and not any(destination.iterdir()):
                destination.rmdir()
            else:
                raise RuntimeError(
                    f"refusing to overwrite incompatible neutral plan: {destination}"
                ) from error
        else:
            return {**result, "status": "verified_cached"}

    destination.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(prefix=f".{destination.name}.staging.", dir=destination.parent)
    )
    try:
        units = build_neutral_simulation_plan(root)
        configs = render_neutral_parameter_files(root)
        for model_id, text in configs.items():
            _atomic_text(staging / _config_path(model_id), text)
        source_plan = framework.FrameworkPlan(
            selection_coefficients=(SOURCE_SELECTED_COEFFICIENT,)
        )
        _atomic_text(
            staging / RECOMBINATION_MAP_FILENAME,
            framework._recombination_map_text(source_plan),
        )
        _atomic_frame(staging / PLAN_FILENAME, units)
        plan_tsv_sha256 = sha256_file(staging / PLAN_FILENAME)
        manifest = _plan_manifest(
            root,
            units=units,
            configs=configs,
            plan_tsv_sha256=plan_tsv_sha256,
        )
        _atomic_json(staging / MANIFEST_FILENAME, manifest)
        outputs = _inventory(staging)
        completion = {
            "schema": PLAN_COMPLETION_SCHEMA,
            "status": "complete",
            "plan_manifest_sha256": sha256_file(staging / MANIFEST_FILENAME),
            "plan_tsv_sha256": plan_tsv_sha256,
            "outputs": outputs,
            "output_count": len(outputs),
        }
        _atomic_json(staging / PLAN_COMPLETION_FILENAME, completion)
        if destination.exists():
            raise RuntimeError("neutral plan destination appeared during generation")
        os.replace(staging, destination)
    finally:
        if staging.exists() and staging.parent == destination.parent:
            shutil.rmtree(staging)
    return verify_neutral_plan(root, plan_dir=destination)


def _validate_plan_frame(frame: pd.DataFrame, manifest: Mapping[str, Any]) -> None:
    required = {
        "plan_order",
        "unit_id",
        "model_id",
        "replicate_index",
        "seed",
        "seed_formula",
        "config_path",
        "config_sha256",
        "selection_coefficient",
        "endpoint_lower",
        "endpoint_upper",
        "sample_af_conditioned",
        "han_intermediate_af_conditioned",
        "adaptive_stopping",
    }
    if not required.issubset(frame.columns) or len(frame) != TOTAL_UNITS:
        raise ValueError("neutral unit plan has the wrong schema or cardinality")
    if frame["unit_id"].duplicated().any() or frame["seed"].duplicated().any():
        raise ValueError("neutral unit IDs and seeds must be unique")
    expected_order = np.arange(TOTAL_UNITS, dtype=np.int64)
    np.testing.assert_array_equal(
        pd.to_numeric(frame["plan_order"], errors="raise").to_numpy(), expected_order
    )
    for model_id, seed_start in (
        (EAS_MODEL_ID, EAS_SEED_START),
        (HAN_MODEL_ID, HAN_SEED_START),
    ):
        group = frame[frame["model_id"].astype(str) == model_id]
        if len(group) != REPLICATES_PER_MODEL:
            raise ValueError(f"neutral plan does not contain 100 {model_id} units")
        indexes = pd.to_numeric(group["replicate_index"], errors="raise").astype(int)
        seeds = pd.to_numeric(group["seed"], errors="raise").astype(int)
        np.testing.assert_array_equal(indexes.to_numpy(), np.arange(100))
        np.testing.assert_array_equal(seeds.to_numpy(), seed_start + np.arange(100))
        if set(group["seed_formula"].astype(str)) != {SEED_FORMULAS[model_id]}:
            raise ValueError(f"{model_id} seed formula changed")
    if set(pd.to_numeric(frame["seed"], errors="raise")).intersection(
        SELECTED_SEEDS_EXCLUDED
    ):
        raise ValueError("neutral plan overlaps a selected seed")
    numeric_s = pd.to_numeric(frame["selection_coefficient"], errors="raise")
    if not np.allclose(numeric_s, 0.0, rtol=0.0, atol=0.0):
        raise ValueError("neutral plan contains nonzero selection")
    lower = pd.to_numeric(frame["endpoint_lower"], errors="raise")
    upper = pd.to_numeric(frame["endpoint_upper"], errors="raise")
    if not np.allclose(lower, SURVIVAL_ENDPOINT_LOWER, rtol=0.0, atol=0.0):
        raise ValueError("neutral plan lower endpoint changed")
    if not np.allclose(upper, 1.0, rtol=0.0, atol=0.0):
        raise ValueError("neutral plan no longer allows fixation")
    for column in (
        "sample_af_conditioned",
        "han_intermediate_af_conditioned",
        "adaptive_stopping",
    ):
        values = frame[column].astype(str).str.casefold()
        if not values.isin({"false", "0"}).all():
            raise ValueError(f"neutral plan unexpectedly enables {column}")
    hashes = manifest.get("unit_record_sha256")
    if not isinstance(hashes, dict) or set(hashes) != set(frame["unit_id"].astype(str)):
        raise ValueError("neutral manifest unit hash inventory is incomplete")
    for row in frame.to_dict(orient="records"):
        if hashes[str(row["unit_id"])] != _canonical_sha256(_portable_record(row)):
            raise ValueError(f"neutral unit record changed: {row['unit_id']}")


def verify_neutral_plan(
    repo_root: str | Path, *, plan_dir: str | Path | None = None
) -> dict[str, Any]:
    """Strictly verify plan inventory, source binding, configs, seeds, and units."""

    root = Path(repo_root).resolve()
    bundle_input = (
        Path(plan_dir) if plan_dir is not None else root / DEFAULT_PLAN_RELATIVE
    )
    bundle = _regular_directory(bundle_input, "neutral plan")
    completion = _json_file(
        bundle / PLAN_COMPLETION_FILENAME, label="neutral plan completion"
    )
    manifest = _json_file(bundle / MANIFEST_FILENAME, label="neutral plan manifest")
    if (
        completion.get("schema") != PLAN_COMPLETION_SCHEMA
        or completion.get("status") != "complete"
    ):
        raise ValueError("neutral plan completion schema/status is invalid")
    if manifest.get("schema") != PLAN_SCHEMA or manifest.get("status") != "planned":
        raise ValueError("neutral plan manifest schema/status is invalid")
    if completion.get("plan_manifest_sha256") != sha256_file(
        bundle / MANIFEST_FILENAME
    ):
        raise ValueError("neutral plan manifest checksum differs from completion")
    actual_outputs = _inventory(bundle, exclude=(PLAN_COMPLETION_FILENAME,))
    if completion.get("outputs") != actual_outputs or completion.get(
        "output_count"
    ) != len(actual_outputs):
        raise ValueError("neutral plan output inventory is incomplete or corrupt")
    source_binding = manifest.get("source_binding")
    if (
        not isinstance(source_binding, dict)
        or source_binding.get("cosi2_commit") != framework.COSI2_COMMIT
    ):
        raise ValueError("neutral plan CoSi2 binding changed")
    expected_sources = _source_binding(root)
    if source_binding != expected_sources:
        raise ValueError("neutral plan source binding is stale")
    plan_path = bundle / PLAN_FILENAME
    if manifest.get("plan_tsv_sha256") != sha256_file(plan_path) or completion.get(
        "plan_tsv_sha256"
    ) != sha256_file(plan_path):
        raise ValueError("neutral unit plan checksum mismatch")
    units = pd.read_csv(plan_path, sep="\t")
    _validate_plan_frame(units, manifest)
    cells = build_neutral_model_cells(root).set_index("model_id")
    for model_id in MODEL_IDS:
        config_path = bundle / _config_path(model_id)
        text = _regular_file(config_path, f"{model_id} config").read_text(
            encoding="utf-8"
        )
        lint_neutral_parameter_text(text, expected_cell=cells.loc[model_id].to_dict())
        expected_hash = manifest["configs"][model_id]["sha256"]
        if sha256_file(config_path) != expected_hash:
            raise ValueError(f"{model_id} config checksum mismatch")
    source_plan = framework.FrameworkPlan(
        selection_coefficients=(SOURCE_SELECTED_COEFFICIENT,)
    )
    if (bundle / RECOMBINATION_MAP_FILENAME).read_text(
        encoding="utf-8"
    ) != framework._recombination_map_text(source_plan):
        raise ValueError("neutral recombination map changed")
    return {
        "schema": PLAN_COMPLETION_SCHEMA,
        "status": "complete",
        "plan_dir": str(bundle),
        "unit_count": len(units),
        "model_counts": {
            str(key): int(value)
            for key, value in units["model_id"].value_counts().sort_index().items()
        },
        "plan_tsv_sha256": sha256_file(plan_path),
        "plan_manifest_sha256": sha256_file(bundle / MANIFEST_FILENAME),
        "plan_completion_sha256": sha256_file(bundle / PLAN_COMPLETION_FILENAME),
    }


def _read_trajectory(path: Path) -> pd.DataFrame:
    _regular_file(path, "saved trajectory")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as error:
            raise ValueError("saved trajectory is empty") from error
    if not header or len(header) != len(set(header)):
        raise ValueError("saved trajectory header is empty or duplicated")
    try:
        return pd.read_csv(path, sep="\t")
    except (OSError, UnicodeDecodeError, pd.errors.ParserError) as error:
        raise ValueError("saved trajectory is not a valid TSV") from error


def validate_neutral_trajectory(
    trajectory_path: str | Path, *, unit: Mapping[str, Any]
) -> dict[str, Any]:
    """Validate a complete trajectory, conditioning only on population survival."""

    path = Path(trajectory_path)
    frame = _read_trajectory(path)
    model_id = str(unit["model_id"])
    if model_id not in MODEL_IDS:
        raise ValueError(f"unsupported neutral model: {model_id}")
    population_ids = (1,) if model_id == EAS_MODEL_ID else (1, 2, 3, 4, 5)
    expected_frequency_columns = {
        f"selfreq_{identifier}" for identifier in population_ids
    }
    expected_size_columns = {f"popsize_{identifier}" for identifier in population_ids}
    if not {"sim", "gen", *expected_frequency_columns, *expected_size_columns}.issubset(
        frame.columns
    ):
        raise ValueError("saved trajectory lacks required population columns")
    if set(frame["sim"].astype(str)) not in ({"1"}, {"1.0"}):
        raise ValueError("each unit trajectory must contain exactly simulation 1")
    generations = pd.to_numeric(frame["gen"], errors="raise").to_numpy(dtype=float)
    birth_generation = int(unit["mutation_birth_generations_ago"])
    expected_generations = np.arange(birth_generation, -1, -1, dtype=float)
    if len(generations) != len(expected_generations) or not np.array_equal(
        generations, expected_generations
    ):
        raise ValueError(
            "saved trajectory must contain every integer generation from birth to present"
        )
    frequencies = frame[sorted(expected_frequency_columns)].apply(
        pd.to_numeric, errors="raise"
    )
    sizes = frame[sorted(expected_size_columns)].apply(pd.to_numeric, errors="raise")
    if not np.isfinite(frequencies.to_numpy(dtype=float)).all() or (
        (frequencies < 0) | (frequencies > 1)
    ).any(axis=None):
        raise ValueError("saved trajectory contains an invalid frequency")
    if not np.isfinite(sizes.to_numpy(dtype=float)).all() or (sizes <= 0).any(
        axis=None
    ):
        raise ValueError("saved trajectory contains an invalid population size")
    sample_id = 1 if model_id == EAS_MODEL_ID else 3
    birth_id = 1 if model_id == EAS_MODEL_ID else 4
    endpoint = float(frequencies.iloc[-1][f"selfreq_{sample_id}"])
    if not SURVIVAL_ENDPOINT_LOWER <= endpoint <= SURVIVAL_ENDPOINT_UPPER:
        raise ValueError("neutral trajectory does not survive in the sample population")
    birth_frequency = float(frequencies.iloc[0][f"selfreq_{birth_id}"])
    birth_size = float(sizes.iloc[0][f"popsize_{birth_id}"])
    expected_birth_frequency = 1.0 / (2.0 * birth_size)
    if not math.isclose(
        birth_frequency, expected_birth_frequency, rel_tol=1e-9, abs_tol=1e-15
    ):
        raise ValueError("neutral trajectory does not begin from one donor copy")
    other_birth = [
        float(frequencies.iloc[0][f"selfreq_{identifier}"])
        for identifier in population_ids
        if identifier != birth_id
    ]
    if any(value != 0.0 for value in other_birth):
        raise ValueError(
            "neutral allele is present outside its birth population at birth"
        )
    record: dict[str, Any] = {
        "status": "valid",
        "model_id": model_id,
        "rows": len(frame),
        "birth_generation": birth_generation,
        "birth_population_frequency": birth_frequency,
        "present_population_frequency": endpoint,
        "survival_condition_passed": True,
        "final_population_af_conditioned": False,
        "sample_af_conditioned": False,
        "sample_carrier_conditioned": False,
        "han_intermediate_af_conditioned": False,
    }
    if model_id == HAN_MODEL_ID:
        boundary = framework.HAN_MIGRATION_END_GENERATIONS
        boundary_row = frame.loc[generations == float(boundary)]
        if len(boundary_row) != 1:
            raise ValueError("Han trajectory lacks one generation-645 row")
        record["han_generation_645_chb_af_descriptive"] = float(
            boundary_row.iloc[0]["selfreq_3"]
        )
    return record


def validate_cosi_ms(ms_path: str | Path, *, unit: Mapping[str, Any]) -> dict[str, Any]:
    """Stream-validate one compact ``-m`` CoSi2 replicate without loading haplotypes."""

    path = _regular_file(Path(ms_path), "CoSi2 ms output")
    sample_haploids = int(unit["sample_haploids"])
    expected_seed = int(unit["seed"])
    segsites: int | None = None
    positions: list[float] | None = None
    haplotype_count = 0
    derived_counts: np.ndarray | None = None
    saw_ms_header = False
    saw_rng_header = False
    saw_rep_header = False
    phase = "header"
    with path.open("r", encoding="utf-8", newline="") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith(("E ", "muttimes:")):
                raise ValueError("CoSi2 output contains forbidden -e or -M content")
            if phase == "header":
                if stripped.startswith("ms "):
                    if saw_ms_header:
                        raise ValueError("CoSi2 ms header is duplicated")
                    fields = stripped.split()
                    if fields != ["ms", str(sample_haploids), "1"]:
                        raise ValueError(
                            "CoSi2 ms sample/replicate header is incorrect"
                        )
                    saw_ms_header = True
                elif stripped.startswith("cosi_rand "):
                    if not saw_ms_header or int(stripped.split()[1]) != expected_seed:
                        raise ValueError("CoSi2 RNG header seed differs from the plan")
                    saw_rng_header = True
                elif stripped.startswith("// seed="):
                    if (
                        not saw_rng_header
                        or int(stripped.split("=", 1)[1]) != expected_seed
                    ):
                        raise ValueError("CoSi2 replicate seed differs from the plan")
                    saw_rep_header = True
                elif stripped.startswith("segsites:"):
                    if not (saw_ms_header and saw_rng_header and saw_rep_header):
                        raise ValueError("CoSi2 ms headers are incomplete")
                    fields = stripped.split()
                    if len(fields) != 2:
                        raise ValueError("CoSi2 segsites line is malformed")
                    segsites = int(fields[1])
                    if segsites <= 0:
                        raise ValueError("10-Mb CoSi2 output has no segregating sites")
                    derived_counts = np.zeros(segsites, dtype=np.int64)
                    phase = "positions"
                else:
                    raise ValueError(
                        f"unexpected compact CoSi2 header line: {stripped[:40]}"
                    )
            elif phase == "positions":
                if not stripped.startswith("positions:"):
                    raise ValueError("CoSi2 positions line is missing")
                assert segsites is not None
                fields = stripped.split()[1:]
                if len(fields) != segsites:
                    raise ValueError("CoSi2 position count differs from segsites")
                try:
                    positions = [float(value) for value in fields]
                except ValueError as error:
                    raise ValueError(
                        "CoSi2 positions contain a nonnumeric value"
                    ) from error
                values = np.asarray(positions, dtype=float)
                if (
                    not np.isfinite(values).all()
                    or (values <= 0).any()
                    or (values > 1).any()
                    or (np.diff(values) < 0).any()
                ):
                    raise ValueError(
                        "CoSi2 positions are nonfinite, out of range, or unsorted"
                    )
                phase = "haplotypes"
            else:
                assert segsites is not None and derived_counts is not None
                if len(stripped) != segsites or set(stripped).difference({"0", "1"}):
                    raise ValueError("CoSi2 haplotype row is malformed")
                if haplotype_count >= sample_haploids:
                    raise ValueError("CoSi2 output has too many haplotype rows")
                derived_counts += np.fromiter(
                    (character == "1" for character in stripped),
                    dtype=np.int8,
                    count=segsites,
                )
                haplotype_count += 1
    if phase != "haplotypes" or positions is None or derived_counts is None:
        raise ValueError("CoSi2 compact ms output is truncated")
    if haplotype_count != sample_haploids:
        raise ValueError("CoSi2 haplotype count differs from the plan")
    focal_indexes = np.flatnonzero(
        np.asarray(positions, dtype=float) == framework.FOCAL_RELATIVE_POSITION
    )
    focal_alt_count: int | None = None
    if len(focal_indexes) == 1:
        focal_alt_count = int(derived_counts[int(focal_indexes[0])])
    return {
        "status": "valid",
        "sample_haploids": sample_haploids,
        "segsites": len(positions),
        "haplotype_rows": haplotype_count,
        "reported_focal_position_count": len(focal_indexes),
        "focal_alt_count_if_unambiguous": focal_alt_count,
        "sample_carrier_conditioned": False,
        "compact_m_only": True,
    }


def _unit_record_sha256(unit: Mapping[str, Any]) -> str:
    return _canonical_sha256(_portable_record(unit))


def _unit_directory(run_root: Path, unit_id: str) -> Path:
    return run_root / "units" / unit_id


def _completion_expected(
    completion: Mapping[str, Any],
    *,
    unit: Mapping[str, Any],
    plan_tsv_sha256: str,
    binary_sha256: str | None,
) -> None:
    if (
        completion.get("schema") != UNIT_COMPLETION_SCHEMA
        or completion.get("status") != "complete"
    ):
        raise ValueError("neutral unit completion schema/status is invalid")
    expected = {
        "unit_id": str(unit["unit_id"]),
        "unit_record_sha256": _unit_record_sha256(unit),
        "plan_tsv_sha256": plan_tsv_sha256,
        "config_sha256": str(unit["config_sha256"]),
        "seed": int(unit["seed"]),
        "seed_formula": str(unit["seed_formula"]),
        "replicate_index": int(unit["replicate_index"]),
    }
    for key, value in expected.items():
        if completion.get(key) != value:
            raise ValueError(f"neutral unit completion has stale {key}")
    if (
        binary_sha256 is not None
        and completion.get("cosi2_binary_sha256") != binary_sha256
    ):
        raise ValueError("neutral unit was produced by a different CoSi2 binary")


def verify_neutral_unit(
    plan_dir: str | Path,
    run_dir: str | Path,
    *,
    unit: Mapping[str, Any],
    cosi2_binary_sha256: str | None = None,
) -> dict[str, Any]:
    """Verify one completion and every checksum/semantic output it references."""

    bundle = Path(plan_dir).resolve()
    run_root = Path(run_dir).resolve()
    unit_dir = _unit_directory(run_root, str(unit["unit_id"]))
    if not unit_dir.is_dir() or unit_dir.is_symlink():
        raise ValueError(f"neutral unit directory is absent: {unit_dir}")
    completion_path = unit_dir / UNIT_COMPLETION_FILENAME
    completion = _json_file(completion_path, label="neutral unit completion")
    plan_sha = sha256_file(bundle / PLAN_FILENAME)
    _completion_expected(
        completion,
        unit=unit,
        plan_tsv_sha256=plan_sha,
        binary_sha256=cosi2_binary_sha256,
    )
    expected_outputs = completion.get("outputs")
    if not isinstance(expected_outputs, dict):
        raise TypeError("neutral unit completion lacks an output inventory")
    actual_outputs = _inventory(unit_dir, exclude=(UNIT_COMPLETION_FILENAME,))
    if expected_outputs != actual_outputs:
        raise ValueError("neutral unit output inventory is incomplete or corrupt")
    trajectory_audit = validate_neutral_trajectory(
        unit_dir / "trajectory.tsv", unit=unit
    )
    ms_audit = validate_cosi_ms(unit_dir / "simulation.ms", unit=unit)
    if completion.get("trajectory_validation") != trajectory_audit:
        raise ValueError("neutral trajectory validation record is stale")
    if completion.get("ms_validation") != ms_audit:
        raise ValueError("neutral ms validation record is stale")
    binary_path = completion.get("cosi2_binary_path")
    if not isinstance(binary_path, str) or not binary_path:
        raise ValueError("neutral completion lacks the CoSi2 binary path")
    expected_command = [
        binary_path,
        "-p",
        str(unit["config_path"]),
        "-n",
        "1",
        "-r",
        str(int(unit["seed"])),
        "-u",
        "1",
        "-m",
    ]
    if completion.get("command") != expected_command:
        raise ValueError(
            "neutral completion command differs from the exact run contract"
        )
    max_attempts = completion.get("max_attempts")
    if (
        isinstance(max_attempts, bool)
        or not isinstance(max_attempts, int)
        or max_attempts <= 0
    ):
        raise ValueError("neutral completion max_attempts is invalid")
    expected_environment = {
        "inherited_cosi_variables_removed": True,
        "COSI_MAXATTEMPTS": max_attempts,
        "COSI_SAVE_TRAJ": "unit_staging/trajectory.tsv",
    }
    if completion.get("environment_contract") != expected_environment:
        raise ValueError(
            "neutral completion environment contract is stale or incomplete"
        )
    expected_semantics = {
        "survival_only": True,
        "final_population_af_conditioned": False,
        "sample_af_conditioned": False,
        "sample_carrier_conditioned": False,
        "han_intermediate_af_conditioned": False,
        "adaptive_stopping": False,
    }
    if completion.get("neutral_semantics") != expected_semantics:
        raise ValueError("neutral completion semantics differ from the natural null")
    return {
        **completion,
        "completion_path": str(completion_path),
        "completion_sha256": sha256_file(completion_path),
        "unit_dir": str(unit_dir),
    }


def _failure_destination(run_root: Path, unit_id: str) -> Path:
    timestamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    return (
        run_root
        / "failures"
        / unit_id
        / f"attempt_{timestamp}_{os.getpid()}_{uuid.uuid4().hex[:10]}"
    )


def _quarantine_existing(run_root: Path, unit_dir: Path) -> Path:
    quarantine = (
        run_root
        / "quarantine"
        / f"{unit_dir.name}.invalid.{datetime.now(UTC).strftime('%Y%m%dT%H%M%SZ')}."
        f"{uuid.uuid4().hex[:10]}"
    )
    quarantine.parent.mkdir(parents=True, exist_ok=True)
    os.replace(unit_dir, quarantine)
    return quarantine


def _runtime_identity() -> dict[str, str]:
    boot_id = ""
    boot_id_path = Path("/proc/sys/kernel/random/boot_id")
    if boot_id_path.is_file() and not boot_id_path.is_symlink():
        try:
            boot_id = boot_id_path.read_text(encoding="ascii").strip()
        except OSError:
            boot_id = ""
    return {
        "system": platform.system(),
        "node": platform.node(),
        "boot_id": boot_id,
    }


def _pid_is_live(pid: int) -> bool | None:
    """Return True/False for a probed PID, or None when liveness is unknowable."""

    if pid <= 0:
        return None
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    except OSError as error:
        # ESRCH means the process is absent; EPERM means it exists but is owned by
        # another principal.  Anything else is not safe evidence for lock theft.
        if getattr(error, "errno", None) == 3:
            return False
        if getattr(error, "errno", None) == 1:
            return True
        return None
    return True


def _read_lock_record(lock: Path) -> tuple[dict[str, Any] | None, int | None]:
    try:
        text = lock.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return None, None
    try:
        value = json.loads(text)
    except json.JSONDecodeError:
        # A killed writer can leave partial JSON.  Recover a PID only when the
        # literal is unambiguous, then still require an explicit dead-PID probe.
        matches = re.findall(r'"pid"\s*:\s*([0-9]+)', text)
        return None, int(matches[0]) if len(matches) == 1 else None
    if not isinstance(value, dict):
        return None, None
    pid = value.get("pid")
    if isinstance(pid, bool) or not isinstance(pid, int):
        return value, None
    return value, pid


def _archive_reclaimed_lock(
    run_root: Path,
    lock: Path,
    *,
    unit_id: str,
    record: Mapping[str, Any] | None,
    pid: int | None,
    reason: str,
) -> None:
    provenance_dir = run_root / "lock_provenance" / unit_id
    provenance_dir.mkdir(parents=True, exist_ok=True)
    stem = (
        f"reclaimed_{datetime.now(UTC).strftime('%Y%m%dT%H%M%SZ')}_"
        f"{uuid.uuid4().hex[:10]}"
    )
    archived_lock = provenance_dir / f"{stem}.lock"
    os.replace(lock, archived_lock)
    audit = {
        "status": "reclaimed_verified_dead_lock",
        "unit_id": unit_id,
        "reason": reason,
        "pid": pid,
        "record": dict(record) if record is not None else None,
        "current_runtime": _runtime_identity(),
        "archived_lock": archived_lock.name,
        "archived_lock_sha256": sha256_file(archived_lock),
        "reclaimed_utc": _utc_now(),
    }
    _atomic_json(provenance_dir / f"{stem}.json", audit)


def _reclaim_dead_lock(run_root: Path, lock: Path, *, unit_id: str) -> bool:
    """Archive a demonstrably dead lock; never remove a live/unknown lock."""

    record, pid = _read_lock_record(lock)
    if pid is None:
        return False
    current_runtime = _runtime_identity()
    recorded_runtime = record.get("runtime") if isinstance(record, dict) else None
    if isinstance(recorded_runtime, dict):
        same_node = recorded_runtime.get("node") == current_runtime["node"]
        same_system = recorded_runtime.get("system") == current_runtime["system"]
        recorded_boot = str(recorded_runtime.get("boot_id", ""))
        current_boot = current_runtime["boot_id"]
        if not same_node or not same_system:
            return False
        if recorded_boot and current_boot and recorded_boot != current_boot:
            _archive_reclaimed_lock(
                run_root,
                lock,
                unit_id=unit_id,
                record=record,
                pid=pid,
                reason="same_host_prior_boot_id",
            )
            return True
    live = _pid_is_live(pid)
    if live is not False:
        return False
    _archive_reclaimed_lock(
        run_root,
        lock,
        unit_id=unit_id,
        record=record,
        pid=pid,
        reason="same_runtime_pid_absent",
    )
    return True


def _acquire_lock(run_root: Path, unit_id: str) -> tuple[int, Path]:
    lock = run_root / "locks" / f"{unit_id}.lock"
    lock.parent.mkdir(parents=True, exist_ok=True)
    for _attempt in range(3):
        try:
            descriptor = os.open(lock, os.O_WRONLY | os.O_CREAT | os.O_EXCL)
        except FileExistsError as error:
            try:
                reclaimed = _reclaim_dead_lock(run_root, lock, unit_id=unit_id)
            except FileNotFoundError:
                continue
            if reclaimed:
                continue
            raise RuntimeError(
                f"neutral unit has a live or unverifiable lock: {unit_id}"
            ) from error
        record = {
            "unit_id": unit_id,
            "pid": os.getpid(),
            "runtime": _runtime_identity(),
            "started_utc": _utc_now(),
        }
        os.write(descriptor, _canonical_json(record).encode("utf-8"))
        os.fsync(descriptor)
        return descriptor, lock
    raise RuntimeError(f"could not acquire neutral unit lock after retries: {unit_id}")


def _run_neutral_unit(
    *,
    bundle: Path,
    run_root: Path,
    unit: Mapping[str, Any],
    binary: Path,
    binary_sha256: str,
    max_attempts: int,
    timeout_seconds: float | None,
) -> dict[str, Any]:
    unit = _portable_record(unit)
    unit_id = str(unit["unit_id"])
    final_dir = _unit_directory(run_root, unit_id)
    if final_dir.exists():
        try:
            cached = verify_neutral_unit(
                bundle,
                run_root,
                unit=unit,
                cosi2_binary_sha256=binary_sha256,
            )
        except Exception:  # noqa: BLE001 -- corrupt cache is quarantined fail-closed
            _quarantine_existing(run_root, final_dir)
        else:
            return {
                "unit_id": unit_id,
                "model_id": unit["model_id"],
                "replicate_index": int(unit["replicate_index"]),
                "seed": int(unit["seed"]),
                "status": "cached",
                "completion_sha256": cached["completion_sha256"],
                "detail": "validated immutable completion reused",
            }
    descriptor: int | None = None
    lock: Path | None = None
    staging: Path | None = None
    started = _utc_now()
    start_clock = time.monotonic()
    try:
        descriptor, lock = _acquire_lock(run_root, unit_id)
        staging_root = run_root / ".staging"
        staging_root.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(prefix=f"{unit_id}.", dir=staging_root))
        trajectory_path = (staging / "trajectory.tsv").resolve()
        ms_path = staging / "simulation.ms"
        stderr_path = staging / "simulation.stderr.log"
        command = [
            str(binary),
            "-p",
            str(unit["config_path"]),
            "-n",
            "1",
            "-r",
            str(int(unit["seed"])),
            "-u",
            "1",
            "-m",
        ]
        if "-e" in command or "-M" in command or command.count("-m") != 1:
            raise RuntimeError("internal error: neutral command is not compact -m only")
        environment = {
            key: value
            for key, value in os.environ.items()
            if not key.startswith("COSI_")
        }
        environment.update(
            {
                "COSI_MAXATTEMPTS": str(max_attempts),
                "COSI_SAVE_TRAJ": str(trajectory_path),
            }
        )
        with (
            ms_path.open("wb") as stdout_handle,
            stderr_path.open("wb") as stderr_handle,
        ):
            process = subprocess.run(
                command,
                cwd=bundle,
                env=environment,
                stdout=stdout_handle,
                stderr=stderr_handle,
                check=False,
                timeout=timeout_seconds,
            )
        if process.returncode != 0:
            raise RuntimeError(f"CoSi2 exited with status {process.returncode}")
        trajectory_audit = validate_neutral_trajectory(trajectory_path, unit=unit)
        ms_audit = validate_cosi_ms(ms_path, unit=unit)
        outputs = _inventory(staging)
        completion = {
            "schema": UNIT_COMPLETION_SCHEMA,
            "status": "complete",
            "unit_id": unit_id,
            "model_id": unit["model_id"],
            "replicate_index": int(unit["replicate_index"]),
            "seed": int(unit["seed"]),
            "seed_formula": unit["seed_formula"],
            "unit_record_sha256": _unit_record_sha256(unit),
            "plan_tsv_sha256": sha256_file(bundle / PLAN_FILENAME),
            "config_sha256": unit["config_sha256"],
            "cosi2_binary_path": str(binary),
            "cosi2_binary_sha256": binary_sha256,
            "cosi2_commit_expected": framework.COSI2_COMMIT,
            "command": command,
            "environment_contract": {
                "inherited_cosi_variables_removed": True,
                "COSI_MAXATTEMPTS": max_attempts,
                "COSI_SAVE_TRAJ": "unit_staging/trajectory.tsv",
            },
            "max_attempts": max_attempts,
            "started_utc": started,
            "completed_utc": _utc_now(),
            "elapsed_seconds": time.monotonic() - start_clock,
            "trajectory_validation": trajectory_audit,
            "ms_validation": ms_audit,
            "neutral_semantics": {
                "survival_only": True,
                "final_population_af_conditioned": False,
                "sample_af_conditioned": False,
                "sample_carrier_conditioned": False,
                "han_intermediate_af_conditioned": False,
                "adaptive_stopping": False,
            },
            "outputs": outputs,
        }
        _atomic_json(staging / UNIT_COMPLETION_FILENAME, completion)
        final_dir.parent.mkdir(parents=True, exist_ok=True)
        if final_dir.exists():
            raise RuntimeError("neutral unit destination appeared during execution")
        os.replace(staging, final_dir)
        staging = None
        verified = verify_neutral_unit(
            bundle,
            run_root,
            unit=unit,
            cosi2_binary_sha256=binary_sha256,
        )
        return {
            "unit_id": unit_id,
            "model_id": unit["model_id"],
            "replicate_index": int(unit["replicate_index"]),
            "seed": int(unit["seed"]),
            "status": "complete",
            "completion_sha256": verified["completion_sha256"],
            "detail": "CoSi2 output completed and validated",
        }
    # A worker boundary must preserve provenance for any Python/subprocess failure
    # while allowing every other predeclared unit to finish.
    except Exception as error:  # noqa: BLE001
        if staging is not None and staging.exists():
            failure = {
                "schema": UNIT_FAILURE_SCHEMA,
                "status": "failed",
                "unit_id": unit_id,
                "model_id": unit["model_id"],
                "replicate_index": int(unit["replicate_index"]),
                "seed": int(unit["seed"]),
                "seed_formula": unit["seed_formula"],
                "unit_record_sha256": _unit_record_sha256(unit),
                "plan_tsv_sha256": sha256_file(bundle / PLAN_FILENAME),
                "config_sha256": unit["config_sha256"],
                "cosi2_binary_sha256": binary_sha256,
                "max_attempts": max_attempts,
                "started_utc": started,
                "failed_utc": _utc_now(),
                "elapsed_seconds": time.monotonic() - start_clock,
                "error_type": type(error).__name__,
                "error": str(error),
            }
            _atomic_json(staging / "failure.json", failure)
            destination = _failure_destination(run_root, unit_id)
            destination.parent.mkdir(parents=True, exist_ok=True)
            os.replace(staging, destination)
            staging = None
            failure_path = str(destination)
        else:
            failure_path = ""
        return {
            "unit_id": unit_id,
            "model_id": unit["model_id"],
            "replicate_index": int(unit["replicate_index"]),
            "seed": int(unit["seed"]),
            "status": "failed",
            "completion_sha256": "",
            "detail": str(error),
            "failure_path": failure_path,
        }
    finally:
        if (
            staging is not None
            and staging.exists()
            and staging.parent.name == ".staging"
        ):
            shutil.rmtree(staging)
        if descriptor is not None:
            os.close(descriptor)
        if lock is not None:
            lock.unlink(missing_ok=True)


def _load_plan_units(bundle: Path) -> pd.DataFrame:
    manifest = _json_file(bundle / MANIFEST_FILENAME, label="neutral plan manifest")
    units = pd.read_csv(bundle / PLAN_FILENAME, sep="\t")
    _validate_plan_frame(units, manifest)
    return units


def run_neutral_simulations(
    repo_root: str | Path,
    *,
    cosi2_binary: str | Path,
    plan_dir: str | Path | None = None,
    run_dir: str | Path | None = None,
    max_workers: int = DEFAULT_MAX_WORKERS,
    max_attempts: int = DEFAULT_COSI_MAXATTEMPTS,
    timeout_seconds: float | None = None,
    unit_ids: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Run every requested predeclared unit; never stop based on observed outcomes.

    ``unit_ids`` supports fixed batches and smoke/resume operation.  It selects
    from, but never rewrites, the immutable 200-unit ledger.
    """

    if isinstance(max_workers, bool) or not 1 <= int(max_workers) <= MAX_WORKERS:
        raise ValueError(f"max_workers must be an integer in [1, {MAX_WORKERS}]")
    if isinstance(max_attempts, bool) or int(max_attempts) <= 0:
        raise ValueError("max_attempts must be a positive integer")
    if timeout_seconds is not None and (
        not math.isfinite(timeout_seconds) or timeout_seconds <= 0
    ):
        raise ValueError("timeout_seconds must be positive when provided")
    root = Path(repo_root).resolve()
    bundle = (
        Path(plan_dir).resolve()
        if plan_dir is not None
        else root / DEFAULT_PLAN_RELATIVE
    )
    verify_neutral_plan(root, plan_dir=bundle)
    run_root = (
        Path(run_dir).resolve() if run_dir is not None else root / DEFAULT_RUN_RELATIVE
    )
    run_root.mkdir(parents=True, exist_ok=True)
    binary = _regular_file(Path(cosi2_binary), "CoSi2 binary")
    binary_hash = sha256_file(binary)
    units = _load_plan_units(bundle)
    if unit_ids is not None:
        requested = list(dict.fromkeys(str(value) for value in unit_ids))
        unknown = sorted(set(requested).difference(units["unit_id"].astype(str)))
        if unknown:
            raise ValueError(f"unknown neutral unit IDs: {', '.join(unknown)}")
        order = {value: index for index, value in enumerate(requested)}
        units = units[units["unit_id"].astype(str).isin(requested)].copy()
        units["_requested_order"] = units["unit_id"].map(order)
        units = units.sort_values("_requested_order", kind="mergesort").drop(
            columns="_requested_order"
        )
    records = units.to_dict(orient="records")
    results: list[dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=int(max_workers)) as pool:
        futures = {
            pool.submit(
                _run_neutral_unit,
                bundle=bundle,
                run_root=run_root,
                unit=unit,
                binary=binary,
                binary_sha256=binary_hash,
                max_attempts=int(max_attempts),
                timeout_seconds=timeout_seconds,
            ): int(unit["plan_order"])
            for unit in records
        }
        for future in as_completed(futures):
            result = future.result()
            result["plan_order"] = futures[future]
            results.append(result)
    if not results:
        return pd.DataFrame(
            columns=[
                "plan_order",
                "unit_id",
                "model_id",
                "replicate_index",
                "seed",
                "status",
                "completion_sha256",
                "detail",
            ]
        )
    return (
        pd.DataFrame(results)
        .sort_values("plan_order", kind="mergesort")
        .reset_index(drop=True)
    )


def neutral_simulation_status(
    repo_root: str | Path,
    *,
    plan_dir: str | Path | None = None,
    run_dir: str | Path | None = None,
    cosi2_binary: str | Path | None = None,
) -> pd.DataFrame:
    """Return one fail-closed status row for every predeclared neutral unit."""

    root = Path(repo_root).resolve()
    bundle = (
        Path(plan_dir).resolve()
        if plan_dir is not None
        else root / DEFAULT_PLAN_RELATIVE
    )
    verify_neutral_plan(root, plan_dir=bundle)
    run_root = (
        Path(run_dir).resolve() if run_dir is not None else root / DEFAULT_RUN_RELATIVE
    )
    binary_hash = (
        sha256_file(_regular_file(Path(cosi2_binary), "CoSi2 binary"))
        if cosi2_binary is not None
        else None
    )
    rows: list[dict[str, Any]] = []
    for unit in _load_plan_units(bundle).to_dict(orient="records"):
        unit_id = str(unit["unit_id"])
        unit_dir = _unit_directory(run_root, unit_id)
        lock = run_root / "locks" / f"{unit_id}.lock"
        status = "pending"
        detail = "no completion found"
        completion_sha = ""
        if unit_dir.exists():
            try:
                verified = verify_neutral_unit(
                    bundle,
                    run_root,
                    unit=unit,
                    cosi2_binary_sha256=binary_hash,
                )
            except Exception as error:  # noqa: BLE001 -- status reports corruption
                status = "invalid_completion"
                detail = str(error)
            else:
                status = "complete"
                detail = "validated immutable completion"
                completion_sha = verified["completion_sha256"]
        elif lock.is_file():
            status = "running"
            detail = "unit lock exists; inspect PID before treating it as stale"
        else:
            failures = sorted((run_root / "failures" / unit_id).glob("attempt_*"))
            if failures:
                status = "failed"
                detail = f"latest preserved failure: {failures[-1]}"
        rows.append(
            {
                "plan_order": int(unit["plan_order"]),
                "unit_id": unit_id,
                "model_id": unit["model_id"],
                "replicate_index": int(unit["replicate_index"]),
                "seed": int(unit["seed"]),
                "status": status,
                "detail": detail,
                "completion_sha256": completion_sha,
            }
        )
    return (
        pd.DataFrame(rows)
        .sort_values("plan_order", kind="mergesort")
        .reset_index(drop=True)
    )


def collect_neutral_metadata(
    repo_root: str | Path,
    *,
    plan_dir: str | Path | None = None,
    run_dir: str | Path | None = None,
    cosi2_binary: str | Path | None = None,
    require_all: bool = False,
) -> pd.DataFrame:
    """Collect strictly verified completion metadata without reading haplotypes again."""

    root = Path(repo_root).resolve()
    bundle = (
        Path(plan_dir).resolve()
        if plan_dir is not None
        else root / DEFAULT_PLAN_RELATIVE
    )
    verify_neutral_plan(root, plan_dir=bundle)
    run_root = (
        Path(run_dir).resolve() if run_dir is not None else root / DEFAULT_RUN_RELATIVE
    )
    binary_hash = (
        sha256_file(_regular_file(Path(cosi2_binary), "CoSi2 binary"))
        if cosi2_binary is not None
        else None
    )
    rows: list[dict[str, Any]] = []
    incomplete: list[str] = []
    for unit in _load_plan_units(bundle).to_dict(orient="records"):
        try:
            completion = verify_neutral_unit(
                bundle,
                run_root,
                unit=unit,
                cosi2_binary_sha256=binary_hash,
            )
        except Exception:  # noqa: BLE001 -- collection skips every invalid completion
            incomplete.append(str(unit["unit_id"]))
            continue
        trajectory = completion["trajectory_validation"]
        ms_audit = completion["ms_validation"]
        rows.append(
            {
                "plan_order": int(unit["plan_order"]),
                "unit_id": unit["unit_id"],
                "model_id": unit["model_id"],
                "replicate_index": int(unit["replicate_index"]),
                "seed": int(unit["seed"]),
                "seed_formula": unit["seed_formula"],
                "present_population_frequency": trajectory[
                    "present_population_frequency"
                ],
                "han_generation_645_chb_af_descriptive": trajectory.get(
                    "han_generation_645_chb_af_descriptive"
                ),
                "segsites": ms_audit["segsites"],
                "reported_focal_position_count": ms_audit[
                    "reported_focal_position_count"
                ],
                "focal_alt_count_if_unambiguous": ms_audit[
                    "focal_alt_count_if_unambiguous"
                ],
                "sample_carrier_conditioned": False,
                "completion_sha256": completion["completion_sha256"],
                "simulation_ms_sha256": completion["outputs"]["simulation.ms"][
                    "sha256"
                ],
                "trajectory_sha256": completion["outputs"]["trajectory.tsv"]["sha256"],
                "cosi2_binary_sha256": completion["cosi2_binary_sha256"],
                "max_attempts": completion["max_attempts"],
                "elapsed_seconds": completion["elapsed_seconds"],
            }
        )
    if require_all and incomplete:
        raise RuntimeError(
            f"neutral bank is incomplete: {len(incomplete)} of {TOTAL_UNITS} units missing/invalid"
        )
    if not rows:
        return pd.DataFrame(
            columns=[
                "plan_order",
                "unit_id",
                "model_id",
                "replicate_index",
                "seed",
                "present_population_frequency",
            ]
        )
    return (
        pd.DataFrame(rows)
        .sort_values("plan_order", kind="mergesort")
        .reset_index(drop=True)
    )


def verify_neutral_simulations(
    repo_root: str | Path,
    *,
    plan_dir: str | Path | None = None,
    run_dir: str | Path | None = None,
    cosi2_binary: str | Path | None = None,
    require_all: bool = True,
) -> dict[str, Any]:
    """Verify the whole bank and summarize complete/pending/failed units."""

    status = neutral_simulation_status(
        repo_root,
        plan_dir=plan_dir,
        run_dir=run_dir,
        cosi2_binary=cosi2_binary,
    )
    counts = {
        str(key): int(value)
        for key, value in status["status"].value_counts().sort_index().items()
    }
    complete = int(counts.get("complete", 0))
    if require_all and complete != TOTAL_UNITS:
        raise RuntimeError(
            f"neutral bank is not complete: {complete}/{TOTAL_UNITS} validated units; "
            f"status counts={counts}"
        )
    return {
        "status": "complete" if complete == TOTAL_UNITS else "partial",
        "planned_units": TOTAL_UNITS,
        "validated_complete_units": complete,
        "status_counts": counts,
        "all_predeclared_units_retained": len(status) == TOTAL_UNITS,
        "adaptive_stopping": False,
    }


__all__ = [
    "DEFAULT_COSI_MAXATTEMPTS",
    "DEFAULT_MAX_WORKERS",
    "EAS_MODEL_ID",
    "EAS_SEED_START",
    "HAN_MODEL_ID",
    "HAN_SEED_START",
    "MAX_WORKERS",
    "REPLICATES_PER_MODEL",
    "SEED_FORMULAS",
    "SELECTED_SEEDS_EXCLUDED",
    "SURVIVAL_ENDPOINT_LOWER",
    "SURVIVAL_ENDPOINT_TOKEN",
    "SURVIVAL_ENDPOINT_UPPER",
    "TOTAL_UNITS",
    "NeutralSimulationPlan",
    "build_neutral_model_cells",
    "build_neutral_simulation_plan",
    "collect_neutral_metadata",
    "generate_neutral_plan",
    "lint_neutral_parameter_text",
    "neutral_simulation_status",
    "render_neutral_parameter_files",
    "run_neutral_simulations",
    "sha256_file",
    "validate_cosi_ms",
    "validate_neutral_trajectory",
    "verify_neutral_plan",
    "verify_neutral_simulations",
    "verify_neutral_unit",
]
