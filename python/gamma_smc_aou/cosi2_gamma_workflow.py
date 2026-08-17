"""Restartable CoSi2 natural-neutral to Gamma-SMC study workflow.

The workflow deliberately keeps the two accepted selected simulations and all
neutral simulation outputs read-only.  Mutable conversion/decode artifacts live
under ``cosi2_gamma_study_work`` and the final statistical report is promoted as
one checksum-bound directory only after every table and figure is complete.

The inferential inventory is fixed at 202 independent CoSi2 loci: one selected
and 100 natural-neutral survivors for each of EAS and Han.  Decoder processes
are parallelised across loci; every Gamma-SMC process itself uses one thread.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import shutil
import tempfile
from collections.abc import Mapping, Sequence
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from . import cosi2_gamma_analysis as analysis
from . import cosi2_gamma_decode as gamma_decode
from . import cosi2_neutral_simulation as neutral

SCHEMA_VERSION = "gamma-smc.cosi2-gamma-workflow/v1"
RESULTS_SCHEMA = "gamma-smc.cosi2-gamma-results/v1"

DEFAULT_WORK_RELATIVE = Path("focused_selection_EAS_sim/cosi2_gamma_study_work")
DEFAULT_RESULTS_RELATIVE = Path("focused_selection_EAS_sim/results/cosi2_gamma_s0p01")
DEFAULT_SELECTED_RELATIVE = Path("focused_selection_EAS_sim/cosi2_framework_runs")
DEFAULT_COSI2_BINARY_RELATIVE = Path("external/cosi2/coalescent")
DEFAULT_GAMMA_BINARY_RELATIVE = Path("bin/gamma_smc")

NEUTRAL_PLAN_DIRNAME = "neutral_plan"
NEUTRAL_RUN_DIRNAME = "neutral_simulations"
DECODE_DIRNAME = "decode_units"

MAX_WORKERS = 20
EXPECTED_NEUTRAL_PER_DEMOGRAPHY = 100
EXPECTED_SELECTED_PER_DEMOGRAPHY = 1
EXPECTED_UNITS_PER_DEMOGRAPHY = 101
EXPECTED_TOTAL_UNITS = 202
SELECTION_COEFFICIENT = 0.01

DEMOGRAPHY_SETTINGS: dict[str, dict[str, Any]] = {
    "EAS": {
        "demography_id": "EAS_PHLASH_median",
        "neutral_model_id": neutral.EAS_MODEL_ID,
        "selected_cell_id": "eas_s0p01",
        "selected_unit_id": "eas_selected_s0p01",
        "selected_seed": 20_240_524,
        "present_ne": 37_637.0,
        "generation_time_years": 25.0,
        "selected_relative_path": "smoke/eas_s0p01",
        "manifest_key": "eas_manifest",
    },
    "Han": {
        "demography_id": "OutOfAfricaArchaicAdmixture_5R19_CHB",
        "neutral_model_id": neutral.HAN_MODEL_ID,
        "selected_cell_id": "han_introgressed_s0p01",
        "selected_unit_id": "han_selected_s0p01",
        "selected_seed": 20_240_658,
        "present_ne": 65_835.0,
        "generation_time_years": 29.0,
        "selected_relative_path": (
            "pilot_native/han_introgressed_s0p01/final_canonical_seed_20240658"
        ),
        "manifest_key": "han_manifest",
    },
}

UNIT_COLUMNS = (
    "unit_id",
    "demography",
    "demography_id",
    "model_id",
    "simulation_class",
    "selection_coefficient",
    "seed",
    "generation_time_years",
    "present_ne",
    "final_population_af",
    "sample_alt_count",
    "sample_af",
    "ms_path",
    "ms_sha256",
    "source_completion_path",
    "source_completion_sha256",
)

PROFILE_FIGURE_STEM = "cosi2_gamma_selected_vs_neutral_profiles"
RESULT_OUTPUTS = {
    "scores": "cosi2_gamma_scores.tsv",
    "neutral_metadata": "neutral_metadata.tsv",
    "unit_metadata": "unit_metadata.tsv",
    "pointwise_pvalues": "pointwise_simulation_pvalues.tsv",
    "minp_omnibus": "minp_omnibus.tsv",
    "run_results": "RUN_RESULTS.md",
    "profile_png": f"{PROFILE_FIGURE_STEM}.png",
    "profile_pdf": f"{PROFILE_FIGURE_STEM}.pdf",
}
RESULTS_COMPLETION_FILENAME = "RESULTS_COMPLETION.json"
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


def sha256_file(path: str | Path) -> str:
    """Return a streaming SHA-256 digest."""

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
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    try:
        temporary.write_text(value, encoding="utf-8", newline="\n")
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
        raise ValueError(f"{label} must be a regular file: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"{label} is not valid UTF-8 JSON: {path}") from error
    if not isinstance(value, dict):
        raise TypeError(f"{label} must contain a JSON object")
    return value


def _relative(path: Path, root: Path) -> str:
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return str(path.resolve())


def _portable_frame_sha256(frame: pd.DataFrame) -> str:
    """Hash table values stably across TSV numeric dtype inference."""

    rows: list[list[Any]] = []
    for record in frame.itertuples(index=False, name=None):
        portable: list[Any] = []
        for value in record:
            if isinstance(value, np.generic):
                value = value.item()
            if value is None or pd.isna(value):
                portable.append(None)
            elif isinstance(value, (bool, np.bool_)):
                portable.append("true" if value else "false")
            elif isinstance(value, int):
                portable.append(str(value))
            elif isinstance(value, float):
                if not math.isfinite(value):
                    raise ValueError("cannot hash a table containing nonfinite values")
                portable.append(format(value, ".12g"))
            else:
                portable.append(str(value))
        rows.append(portable)
    return _canonical_sha256({"columns": list(frame.columns), "rows": rows})


def _positive_workers(value: Any) -> int:
    if isinstance(value, bool):
        raise TypeError("workers must be an integer")
    workers = int(value)
    if workers != value or not 1 <= workers <= MAX_WORKERS:
        raise ValueError(f"workers must be an integer in [1, {MAX_WORKERS}]")
    return workers


def _validate_manifest_file(directory: Path, manifest_path: Path) -> dict[str, str]:
    if manifest_path.is_symlink() or not manifest_path.is_file():
        raise ValueError(f"selected SHA256SUMS is absent: {manifest_path}")
    records: dict[str, str] = {}
    for line_number, raw in enumerate(
        manifest_path.read_text(encoding="utf-8").splitlines(), start=1
    ):
        parts = raw.split(maxsplit=1)
        if len(parts) != 2 or len(parts[0]) != 64:
            raise ValueError(f"malformed SHA256SUMS line {line_number}")
        digest, name = parts[0], parts[1].lstrip(" *")
        if name in records or Path(name).name != name:
            raise ValueError("selected SHA256SUMS names must be unique basenames")
        path = directory / name
        if path.is_symlink() or not path.is_file() or sha256_file(path) != digest:
            raise ValueError(f"selected checksum failed: {path}")
        records[name] = digest
    if "simulation.ms" not in records:
        raise ValueError("selected SHA256SUMS omits simulation.ms")
    return records


def bind_selected_units(
    repo_root: str | Path,
    *,
    selected_dir: str | Path | None = None,
) -> pd.DataFrame:
    """Validate and bind the two accepted selected CoSi2 outputs read-only."""

    root = Path(repo_root).resolve()
    selected_root = (
        Path(selected_dir).resolve()
        if selected_dir is not None
        else root / DEFAULT_SELECTED_RELATIVE
    )
    if selected_root.is_symlink() or not selected_root.is_dir():
        raise ValueError(f"selected run root is absent: {selected_root}")
    completion_path = selected_root / "RUN_COMPLETION.json"
    summary_path = selected_root / "simulation_summary.tsv"
    completion = _read_json(completion_path, "selected RUN_COMPLETION")
    if completion.get("status") != "complete":
        raise ValueError("selected RUN_COMPLETION is not complete")
    required_contract = {
        "selection_coefficient": SELECTION_COEFFICIENT,
        "cosi2_commit": "5da717cb8bf7324fb961e772a0c9042385de381a",
        "present_population_af_lower": 0.195,
        "present_population_af_upper": 0.205,
    }
    for key, expected in required_contract.items():
        if completion.get(key) != expected:
            raise ValueError(f"selected RUN_COMPLETION has incompatible {key}")
    outputs = completion.get("outputs")
    if not isinstance(outputs, dict):
        raise TypeError("selected RUN_COMPLETION lacks outputs")
    expected_output_keys = {
        "run_results",
        "simulation_summary",
        "eas_manifest",
        "han_manifest",
    }
    if set(outputs) != expected_output_keys:
        raise ValueError("selected RUN_COMPLETION output inventory changed")
    for label, record in outputs.items():
        if not isinstance(record, dict):
            raise TypeError(f"selected output record is malformed: {label}")
        output = selected_root / str(record.get("path", ""))
        if output.is_symlink() or not output.is_file():
            raise ValueError(f"selected output is absent: {output}")
        if sha256_file(output) != record.get("sha256"):
            raise ValueError(f"selected output checksum failed: {label}")
    if Path(outputs["simulation_summary"]["path"]) != Path("simulation_summary.tsv"):
        raise ValueError("selected summary path changed")
    summary = pd.read_csv(summary_path, sep="\t")
    if len(summary) != 2 or summary["cell_id"].duplicated().any():
        raise ValueError("selected summary must contain exactly two unique cells")
    if int(outputs["simulation_summary"].get("rows", -1)) != len(summary):
        raise ValueError("selected summary row count disagrees with completion")

    rows: list[dict[str, Any]] = []
    indexed = summary.set_index("cell_id", drop=False)
    if set(indexed.index.astype(str)) != {
        settings["selected_cell_id"] for settings in DEMOGRAPHY_SETTINGS.values()
    }:
        raise ValueError("selected summary contains the wrong model mapping")
    for demography, settings in DEMOGRAPHY_SETTINGS.items():
        record = indexed.loc[settings["selected_cell_id"]]
        source_relative = Path(settings["selected_relative_path"])
        source_dir = selected_root / source_relative
        if source_dir.is_symlink() or not source_dir.is_dir():
            raise ValueError(f"accepted selected directory is absent: {source_dir}")
        sentinel = source_dir / "COMPLETE"
        if sentinel.is_symlink() or not sentinel.is_file():
            raise ValueError(f"selected COMPLETE sentinel is absent: {source_dir}")
        recorded_relative = Path(str(record["output_relative_path"]))
        expected_relative = Path("cosi2_framework_runs") / source_relative
        if recorded_relative != expected_relative:
            raise ValueError(f"selected {demography} summary path changed")
        numeric_checks = {
            "selection_coefficient": SELECTION_COEFFICIENT,
            "trajectory_seed": int(settings["selected_seed"]),
            "sample_haplotypes": analysis.PANEL_HAPLOTYPES,
        }
        for column, expected in numeric_checks.items():
            if float(record[column]) != float(expected):
                raise ValueError(f"selected {demography} has incompatible {column}")
        if str(record["demography"]) != settings["demography_id"]:
            raise ValueError(f"selected {demography} demography mapping changed")
        if str(record["status"]) != "complete":
            raise ValueError(f"selected {demography} summary is not complete")
        final_af = float(record["present_population_af"])
        sample_count = int(record["present_sample_alt_count"])
        sample_af = float(record["present_sample_af"])
        if not 0.195 <= final_af <= 0.205:
            raise ValueError(f"selected {demography} endpoint AF is outside its gate")
        if not math.isclose(
            sample_af,
            sample_count / analysis.PANEL_HAPLOTYPES,
            rel_tol=0.0,
            abs_tol=1e-12,
        ):
            raise ValueError(f"selected {demography} sampled AF is inconsistent")
        manifest_record = outputs[settings["manifest_key"]]
        manifest_path = selected_root / str(manifest_record["path"])
        if manifest_path != source_dir / "SHA256SUMS":
            raise ValueError(f"selected {demography} manifest path changed")
        files = _validate_manifest_file(source_dir, manifest_path)
        ms_path = source_dir / "simulation.ms"
        if ms_path.stat().st_size != int(record["ms_bytes"]):
            raise ValueError(f"selected {demography} ms byte count changed")
        rows.append(
            {
                "unit_id": settings["selected_unit_id"],
                "demography": demography,
                "demography_id": settings["demography_id"],
                "model_id": settings["selected_cell_id"],
                "simulation_class": "selected",
                "selection_coefficient": SELECTION_COEFFICIENT,
                "seed": int(settings["selected_seed"]),
                "generation_time_years": settings["generation_time_years"],
                "present_ne": settings["present_ne"],
                "final_population_af": final_af,
                "sample_alt_count": sample_count,
                "sample_af": sample_af,
                "ms_path": str(ms_path.resolve()),
                "ms_sha256": files["simulation.ms"],
                "source_completion_path": str(completion_path.resolve()),
                "source_completion_sha256": sha256_file(completion_path),
            }
        )
    return (
        pd.DataFrame(rows, columns=UNIT_COLUMNS)
        .sort_values("demography", kind="mergesort")
        .reset_index(drop=True)
    )


def _neutral_sample_alt_count(record: Mapping[str, Any]) -> int:
    """Apply the prespecified absent-focal-site sample-count convention."""

    population_af = float(record["present_population_frequency"])
    focal_count = record.get("focal_alt_count_if_unambiguous")
    focal_sites = int(record.get("reported_focal_position_count", 0))
    if focal_sites == 1 and focal_count is not None and not pd.isna(focal_count):
        count = int(focal_count)
    elif math.isclose(population_af, 1.0, rel_tol=0.0, abs_tol=1e-12):
        count = analysis.PANEL_HAPLOTYPES
    else:
        count = 0
    if not 0 <= count <= analysis.PANEL_HAPLOTYPES:
        raise ValueError("neutral focal sample count lies outside [0, 200]")
    return count


def collect_neutral_units(
    repo_root: str | Path,
    *,
    plan_dir: str | Path,
    run_dir: str | Path,
    cosi2_binary: str | Path,
) -> pd.DataFrame:
    """Verify all 200 neutral units and return their decode metadata."""

    root = Path(repo_root).resolve()
    plan_root = Path(plan_dir).resolve()
    run_root = Path(run_dir).resolve()
    neutral.verify_neutral_simulations(
        root,
        plan_dir=plan_root,
        run_dir=run_root,
        cosi2_binary=cosi2_binary,
        require_all=True,
    )
    metadata = neutral.collect_neutral_metadata(
        root,
        plan_dir=plan_root,
        run_dir=run_root,
        cosi2_binary=cosi2_binary,
        require_all=True,
    )
    if len(metadata) != neutral.TOTAL_UNITS:
        raise ValueError("neutral metadata does not contain exactly 200 units")
    rows: list[dict[str, Any]] = []
    model_to_demography = {
        settings["neutral_model_id"]: name
        for name, settings in DEMOGRAPHY_SETTINGS.items()
    }
    for record in metadata.to_dict(orient="records"):
        model_id = str(record["model_id"])
        if model_id not in model_to_demography:
            raise ValueError(f"unknown neutral model mapping: {model_id}")
        demography = model_to_demography[model_id]
        settings = DEMOGRAPHY_SETTINGS[demography]
        sample_count = _neutral_sample_alt_count(record)
        unit_id = str(record["unit_id"])
        unit_dir = run_root / "units" / unit_id
        ms_path = unit_dir / "simulation.ms"
        completion_path = unit_dir / neutral.UNIT_COMPLETION_FILENAME
        if sha256_file(ms_path) != record["simulation_ms_sha256"]:
            raise ValueError(f"neutral ms checksum changed: {unit_id}")
        rows.append(
            {
                "unit_id": unit_id,
                "demography": demography,
                "demography_id": settings["demography_id"],
                "model_id": model_id,
                "simulation_class": "neutral",
                "selection_coefficient": 0.0,
                "seed": int(record["seed"]),
                "generation_time_years": settings["generation_time_years"],
                "present_ne": settings["present_ne"],
                "final_population_af": float(record["present_population_frequency"]),
                "sample_alt_count": sample_count,
                "sample_af": sample_count / analysis.PANEL_HAPLOTYPES,
                "ms_path": str(ms_path.resolve()),
                "ms_sha256": str(record["simulation_ms_sha256"]),
                "source_completion_path": str(completion_path.resolve()),
                "source_completion_sha256": str(record["completion_sha256"]),
            }
        )
    return (
        pd.DataFrame(rows, columns=UNIT_COLUMNS)
        .sort_values(["demography", "unit_id"], kind="mergesort")
        .reset_index(drop=True)
    )


def validate_unit_inventory(units: pd.DataFrame) -> pd.DataFrame:
    """Fail closed unless the exact two-demography 100+1 inventory is present."""

    if not isinstance(units, pd.DataFrame):
        raise TypeError("units must be a pandas DataFrame")
    missing = sorted(set(UNIT_COLUMNS).difference(units.columns))
    if missing:
        raise ValueError("unit inventory lacks columns: " + ", ".join(missing))
    result = units.loc[:, UNIT_COLUMNS].copy()
    if len(result) != EXPECTED_TOTAL_UNITS or result["unit_id"].duplicated().any():
        raise ValueError("unit inventory must contain exactly 202 unique units")
    if (
        result["unit_id"].isna().any()
        or (result["unit_id"].astype(str).str.strip() == "").any()
    ):
        raise ValueError("unit IDs must be nonempty")
    if set(result["demography"].astype(str)) != set(DEMOGRAPHY_SETTINGS):
        raise ValueError("unit inventory must contain exactly EAS and Han")
    if result[["demography", "seed"]].duplicated().any():
        raise ValueError("seeds must be unique within demography")
    for demography, settings in DEMOGRAPHY_SETTINGS.items():
        group = result[result["demography"].astype(str) == demography]
        classes = group["simulation_class"].astype(str).value_counts().to_dict()
        expected_classes = {
            "neutral": EXPECTED_NEUTRAL_PER_DEMOGRAPHY,
            "selected": EXPECTED_SELECTED_PER_DEMOGRAPHY,
        }
        if len(group) != EXPECTED_UNITS_PER_DEMOGRAPHY or classes != expected_classes:
            raise ValueError(
                f"{demography} must contain exactly 100 neutral + 1 selected"
            )
        if set(group["demography_id"].astype(str)) != {settings["demography_id"]}:
            raise ValueError(f"{demography} demography_id mapping changed")
        if set(pd.to_numeric(group["present_ne"], errors="raise")) != {
            settings["present_ne"]
        }:
            raise ValueError(f"{demography} present Ne changed")
        if set(pd.to_numeric(group["generation_time_years"], errors="raise")) != {
            settings["generation_time_years"]
        }:
            raise ValueError(f"{demography} generation time changed")
        neutral_group = group[group["simulation_class"].astype(str) == "neutral"]
        selected_group = group[group["simulation_class"].astype(str) == "selected"]
        if set(neutral_group["model_id"].astype(str)) != {settings["neutral_model_id"]}:
            raise ValueError(f"{demography} neutral model mapping changed")
        selected = selected_group.iloc[0]
        if (
            str(selected["unit_id"]) != settings["selected_unit_id"]
            or str(selected["model_id"]) != settings["selected_cell_id"]
            or int(selected["seed"]) != settings["selected_seed"]
        ):
            raise ValueError(f"{demography} selected unit mapping changed")
        expected_seed_start = (
            neutral.EAS_SEED_START if demography == "EAS" else neutral.HAN_SEED_START
        )
        expected_neutral_seeds = set(
            range(
                expected_seed_start,
                expected_seed_start + EXPECTED_NEUTRAL_PER_DEMOGRAPHY,
            )
        )
        if set(pd.to_numeric(neutral_group["seed"], errors="raise")) != (
            expected_neutral_seeds
        ):
            raise ValueError(f"{demography} neutral seed stream changed")
        selected_af = float(selected["final_population_af"])
        if not 0.195 <= selected_af <= 0.205:
            raise ValueError(f"{demography} selected endpoint AF changed")
    coefficients = pd.to_numeric(result["selection_coefficient"], errors="raise")
    classes = result["simulation_class"].astype(str)
    if not np.allclose(coefficients[classes == "neutral"], 0.0, rtol=0, atol=0):
        raise ValueError("neutral units must have s=0")
    if not np.allclose(
        coefficients[classes == "selected"], SELECTION_COEFFICIENT, rtol=0, atol=0
    ):
        raise ValueError("selected units must have s=0.01")
    seeds = pd.to_numeric(result["seed"], errors="coerce").to_numpy(dtype=float)
    if (
        not np.all(np.isfinite(seeds))
        or not np.array_equal(seeds, np.floor(seeds))
        or np.any(seeds < 0)
    ):
        raise ValueError("unit seeds must be nonnegative integers")
    numeric_columns = (
        "final_population_af",
        "sample_alt_count",
        "sample_af",
    )
    for column in numeric_columns:
        values = pd.to_numeric(result[column], errors="coerce").to_numpy(dtype=float)
        if not np.all(np.isfinite(values)):
            raise ValueError(f"unit {column} must be finite")
    population_af = pd.to_numeric(result["final_population_af"])
    if ((population_af <= 0) | (population_af > 1)).any():
        raise ValueError(
            "all units must have a surviving population frequency in (0, 1]"
        )
    counts = pd.to_numeric(result["sample_alt_count"]).to_numpy(dtype=float)
    sample_af = pd.to_numeric(result["sample_af"]).to_numpy(dtype=float)
    if not np.array_equal(counts, np.floor(counts)) or np.any(
        (counts < 0) | (counts > analysis.PANEL_HAPLOTYPES)
    ):
        raise ValueError("sample_alt_count must be an integer in [0, 200]")
    if not np.allclose(
        sample_af,
        counts / analysis.PANEL_HAPLOTYPES,
        rtol=0,
        atol=1e-12,
    ):
        raise ValueError("sample AF is inconsistent with sample count")
    for column in ("ms_sha256", "source_completion_sha256"):
        if (
            not result[column]
            .astype(str)
            .map(lambda value: bool(_SHA256_RE.fullmatch(value)))
            .all()
        ):
            raise ValueError(f"{column} must contain lowercase SHA-256 digests")
    for column in ("ms_path", "source_completion_path"):
        if (
            result[column].isna().any()
            or (result[column].astype(str).str.strip() == "").any()
        ):
            raise ValueError(f"{column} must contain nonempty provenance paths")
    return result.sort_values(
        ["demography", "simulation_class", "unit_id"], kind="mergesort"
    ).reset_index(drop=True)


def build_unit_inventory(
    repo_root: str | Path,
    *,
    plan_dir: str | Path,
    run_dir: str | Path,
    cosi2_binary: str | Path,
    selected_dir: str | Path | None = None,
) -> pd.DataFrame:
    selected = bind_selected_units(repo_root, selected_dir=selected_dir)
    natural = collect_neutral_units(
        repo_root,
        plan_dir=plan_dir,
        run_dir=run_dir,
        cosi2_binary=cosi2_binary,
    )
    return validate_unit_inventory(pd.concat([natural, selected], ignore_index=True))


def _verify_output_manifest(
    directory: Path,
    completion_path: Path,
    *,
    expected_schema: str,
    expected_labels: set[str],
) -> dict[str, Any]:
    completion = _read_json(completion_path, "decode completion")
    if (
        completion.get("schema") != expected_schema
        or completion.get("status") != "complete"
    ):
        raise ValueError(f"decode completion schema/status failed: {completion_path}")
    contract = completion.get("contract")
    if not isinstance(contract, dict) or completion.get(
        "contract_sha256"
    ) != _canonical_sha256(contract):
        raise ValueError(
            f"decode completion contract checksum failed: {completion_path}"
        )
    outputs = completion.get("outputs")
    if not isinstance(outputs, dict) or set(outputs) != expected_labels:
        raise ValueError(f"decode completion output labels changed: {completion_path}")
    expected_files = {completion_path.name}
    for label, record in outputs.items():
        if not isinstance(record, dict):
            raise TypeError(f"decode output record is malformed: {label}")
        relative = Path(str(record.get("path", "")))
        if relative.is_absolute() or len(relative.parts) != 1:
            raise ValueError("decode output path must be a basename")
        output = directory / relative
        expected_files.add(relative.name)
        if (
            output.is_symlink()
            or not output.is_file()
            or output.stat().st_size != int(record.get("size_bytes", -1))
            or sha256_file(output) != record.get("sha256")
        ):
            raise ValueError(f"decode output checksum failed: {output}")
    actual_files = {path.name for path in directory.iterdir() if path.is_file()}
    directories = [path for path in directory.iterdir() if path.is_dir()]
    symlinks = [path for path in directory.iterdir() if path.is_symlink()]
    if actual_files != expected_files or directories or symlinks:
        raise ValueError(f"decode directory has an unexpected inventory: {directory}")
    return completion


def verify_decode_unit(
    unit: Mapping[str, Any],
    *,
    decode_root: str | Path,
    decoder_bin: str | Path,
) -> dict[str, Any]:
    """Strictly and read-only verify one conversion plus Gamma decode."""

    unit_id = str(unit["unit_id"])
    root = Path(decode_root).resolve() / unit_id
    binary = Path(decoder_bin).resolve()
    if root.is_symlink() or not root.is_dir():
        raise ValueError(f"decode unit directory is absent: {root}")
    children = {path.name for path in root.iterdir()}
    if children != {"converted", "decoded"}:
        raise ValueError(f"decode unit has an unexpected inventory: {unit_id}")
    ms_path = Path(str(unit["ms_path"])).resolve()
    if not ms_path.is_file() or sha256_file(ms_path) != str(unit["ms_sha256"]):
        raise ValueError(f"source ms checksum failed: {unit_id}")
    if not binary.is_file():
        raise ValueError(f"Gamma-SMC decoder is absent: {binary}")

    converted = root / "converted"
    conversion = _verify_output_manifest(
        converted,
        converted / "conversion_complete.json",
        expected_schema=gamma_decode.CONVERSION_SCHEMA,
        expected_labels={"vcf", "pairs", "metadata"},
    )
    conversion_contract = conversion["contract"]
    conversion_input = conversion_contract.get("input", {})
    conversion_parameters = conversion_contract.get("parameters", {})
    expected_panel_indices = gamma_decode._pair_panel_indices()
    if (
        Path(str(conversion_input.get("path", ""))).resolve() != ms_path
        or conversion_input.get("sha256") != unit["ms_sha256"]
        or int(conversion_parameters.get("expected_haplotypes", -1))
        != analysis.PANEL_HAPLOTYPES
        or int(conversion_parameters.get("pair_panel_haplotypes", -1))
        != gamma_decode.DEFAULT_PAIR_PANEL_HAPLOTYPES
        or int(conversion_parameters.get("pair_panel_seed", -1))
        != gamma_decode.DEFAULT_PAIR_PANEL_SEED
        or int(conversion_parameters.get("pair_count", -1))
        != gamma_decode.DEFAULT_PAIR_COUNT
        or tuple(conversion_parameters.get("pair_panel_indices", ()))
        != expected_panel_indices
        or conversion_parameters.get("pairing_scheme")
        != "fixed_sha256_ranked_100_haplotypes_all_unordered_pairs"
        or int(conversion_parameters.get("sequence_length_bp", -1))
        != gamma_decode.DEFAULT_SEQUENCE_LENGTH_BP
        or conversion_contract.get("implementation", {}).get("module_sha256")
        != sha256_file(Path(gamma_decode.__file__).resolve())
    ):
        raise ValueError(f"conversion contract differs for {unit_id}")
    if int(conversion["outputs"]["pairs"].get("rows", -1)) != (
        gamma_decode.DEFAULT_PAIR_COUNT
    ):
        raise ValueError(f"conversion pair-panel row count differs for {unit_id}")
    expected_pairs = tuple(
        (left, right)
        for left_index, left in enumerate(expected_panel_indices)
        for right in expected_panel_indices[left_index + 1 :]
    )
    observed_pairs = gamma_decode._data_pair_rows(converted / "overall.pairs.tsv")
    if observed_pairs != expected_pairs:
        raise ValueError(f"conversion pair-panel content differs for {unit_id}")
    conversion_metadata = _read_json(
        converted / "conversion_metadata.json", "conversion metadata"
    )
    focal_site_count = int(conversion_metadata.get("exact_focal_site_count", -1))
    raw_full_count = conversion_metadata.get("exact_focal_sample_alt_count")
    raw_panel_count = conversion_metadata.get("exact_focal_pair_panel_alt_count")
    if focal_site_count == 0:
        full_sample_alt_count = 0
        pair_panel_alt_count = 0
        focal_count_status = "absent_from_sample"
    elif focal_site_count == 1:
        full_sample_alt_count = int(raw_full_count)
        pair_panel_alt_count = int(raw_panel_count)
        focal_count_status = "unique_exact_focal_site"
    elif focal_site_count > 1:
        full_sample_alt_count = -1
        pair_panel_alt_count = -1
        focal_count_status = "ambiguous_exact_coordinate_collision"
    else:
        raise ValueError(f"conversion focal-site count is invalid for {unit_id}")
    if (
        conversion_metadata.get("schema") != gamma_decode.CONVERSION_SCHEMA
        or conversion_metadata.get("contract_sha256")
        != conversion.get("contract_sha256")
        or (
            focal_count_status != "ambiguous_exact_coordinate_collision"
            and full_sample_alt_count != int(unit["sample_alt_count"])
        )
        or (
            focal_count_status == "ambiguous_exact_coordinate_collision"
            and int(unit["sample_alt_count"]) != 0
        )
        or not -1 <= pair_panel_alt_count <= gamma_decode.DEFAULT_PAIR_PANEL_HAPLOTYPES
        or (pair_panel_alt_count == -1)
        != (focal_count_status == "ambiguous_exact_coordinate_collision")
        or int(conversion_metadata.get("pair_count", -1))
        != gamma_decode.DEFAULT_PAIR_COUNT
        or conversion_metadata.get("pair_panel_indices_sha256")
        != gamma_decode._pair_panel_sha256(expected_panel_indices)
    ):
        raise ValueError(f"conversion focal/pair-panel metadata differs for {unit_id}")

    decoded = root / "decoded"
    completion_path = decoded / "decode_complete.json"
    decode_completion = _verify_output_manifest(
        decoded,
        completion_path,
        expected_schema=gamma_decode.DECODE_SCHEMA,
        expected_labels={
            "summary",
            "scores",
            "posterior",
            "posterior_metadata",
            "pair_manifest",
            "run_manifest",
        },
    )
    contract = decode_completion["contract"]
    settings = contract.get("settings", {})
    if (
        contract.get("unit_id") != unit_id
        or contract.get("demography_id") != unit["demography_id"]
        or contract.get("conversion", {}).get("completion_sha256")
        != sha256_file(converted / "conversion_complete.json")
        or contract.get("decoder_binary", {}).get("sha256") != sha256_file(binary)
        or float(settings.get("present_ne", -1)) != float(unit["present_ne"])
        or float(settings.get("generation_time_years", -1))
        != float(unit["generation_time_years"])
        or int(settings.get("threads", -1)) != 1
        or settings.get("recent_call") != "mean"
        or settings.get("pair_selector")
        != "fixed_sha256_ranked_100_haplotypes_all_unordered_pairs"
        or int(settings.get("pair_panel_haplotypes", -1))
        != gamma_decode.DEFAULT_PAIR_PANEL_HAPLOTYPES
        or int(settings.get("pair_panel_seed", -1))
        != gamma_decode.DEFAULT_PAIR_PANEL_SEED
        or int(settings.get("pair_count", -1)) != gamma_decode.DEFAULT_PAIR_COUNT
        or int(settings.get("n_pairs", -1)) != gamma_decode.DEFAULT_PAIR_COUNT
        or tuple(settings.get("thresholds_years", ()))
        != gamma_decode.TMRCA_THRESHOLDS_YEARS
        or contract.get("implementation", {}).get("module_sha256")
        != sha256_file(Path(gamma_decode.__file__).resolve())
    ):
        raise ValueError(f"Gamma decode contract differs for {unit_id}")
    if int(decode_completion["outputs"]["pair_manifest"].get("rows", -1)) != (
        gamma_decode.DEFAULT_PAIR_COUNT
    ):
        raise ValueError(f"decoded pair-manifest row count differs for {unit_id}")
    if gamma_decode._data_pair_rows(decoded / "decoded_pairs.tsv") != expected_pairs:
        raise ValueError(f"decoded pair-panel content differs for {unit_id}")
    reduced = gamma_decode.reduce_gamma_summary(decoded / "overall.summary.tsv")
    durable = pd.read_csv(decoded / "unit_scores.tsv", sep="\t")
    expected = reduced.copy()
    expected.insert(0, "demography_id", str(unit["demography_id"]))
    expected.insert(0, "unit_id", unit_id)
    expected.insert(2, "generation_time_years", float(unit["generation_time_years"]))
    pd.testing.assert_frame_equal(
        durable,
        expected,
        check_dtype=False,
        check_exact=False,
        rtol=0.0,
        atol=1e-12,
    )
    return {
        "unit_id": unit_id,
        "status": "complete",
        "conversion_completion_sha256": sha256_file(
            converted / "conversion_complete.json"
        ),
        "decode_completion_sha256": sha256_file(completion_path),
        "scores_path": str((decoded / "unit_scores.tsv").resolve()),
        "scores_sha256": sha256_file(decoded / "unit_scores.tsv"),
        "decode_full_sample_alt_count": full_sample_alt_count,
        "decode_pair_panel_alt_count": pair_panel_alt_count,
        "decode_focal_count_status": focal_count_status,
    }


def run_decodes(
    units: pd.DataFrame,
    *,
    decode_root: str | Path,
    decoder_bin: str | Path,
    max_workers: int = MAX_WORKERS,
    unit_ids: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Launch up to 20 concurrent one-thread Gamma-SMC subprocesses."""

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
        verified = verify_decode_unit(
            record,
            decode_root=root,
            decoder_bin=binary,
        )
        return {**verified, "cache_hit": artifacts.cache_hit}

    results: list[dict[str, Any]] = []
    failures: list[tuple[str, str]] = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(run_one, record): str(record["unit_id"])
            for record in selected.drop(columns="_order").to_dict(orient="records")
        }
        for future in as_completed(futures):
            try:
                results.append(future.result())
            except Exception as error:  # noqa: BLE001 - collect all worker failures
                failures.append((futures[future], str(error)))
    if failures:
        detail = "; ".join(f"{unit_id}: {error}" for unit_id, error in failures[:10])
        raise RuntimeError(f"{len(failures)} Gamma decode units failed: {detail}")
    result = pd.DataFrame(results)
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
    """Return read-only decode status for an exact unit inventory."""

    inventory = validate_unit_inventory(units)
    rows: list[dict[str, Any]] = []
    root = Path(decode_root).resolve()
    for record in inventory.to_dict(orient="records"):
        unit_id = str(record["unit_id"])
        try:
            verified = verify_decode_unit(
                record,
                decode_root=root,
                decoder_bin=decoder_bin,
            )
        except Exception as error:  # noqa: BLE001 - status must report corrupt units
            status = "invalid" if (root / unit_id).exists() else "pending"
            rows.append({"unit_id": unit_id, "status": status, "detail": str(error)})
        else:
            rows.append({**verified, "detail": "validated immutable decode"})
    return pd.DataFrame(rows)


def verify_decode_bank(
    units: pd.DataFrame,
    *,
    decode_root: str | Path,
    decoder_bin: str | Path,
) -> pd.DataFrame:
    status = decode_status(
        units,
        decode_root=decode_root,
        decoder_bin=decoder_bin,
    )
    if len(status) != EXPECTED_TOTAL_UNITS or set(status["status"]) != {"complete"}:
        counts = status["status"].value_counts().to_dict()
        raise RuntimeError(f"Gamma decode bank is incomplete or invalid: {counts}")
    return status.sort_values("unit_id", kind="mergesort").reset_index(drop=True)


def collect_decoded_scores(
    units: pd.DataFrame,
    *,
    decode_root: str | Path,
    decoder_bin: str | Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Verify all decode completions and map both decoder scores exactly."""

    inventory = validate_unit_inventory(units)
    decode_records = verify_decode_bank(
        inventory,
        decode_root=decode_root,
        decoder_bin=decoder_bin,
    )
    metadata = inventory.set_index("unit_id")
    rows: list[dict[str, Any]] = []
    statistic_mapping = {
        "hard_mean_call_fraction": "paper_frac_posterior_mean_below",
        "soft_mean_cdf_score": "mean_p_tmrca_lt",
    }
    for record in decode_records.to_dict(orient="records"):
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
    expected_rows = (
        EXPECTED_TOTAL_UNITS
        * len(analysis.SUPPORTED_STATISTICS)
        * len(analysis.SUPPORTED_WINDOWS)
        * len(analysis.TMRCA_THRESHOLDS_YEARS)
    )
    if len(scores) != expected_rows:
        raise ValueError("assembled score table has the wrong cardinality")
    return scores, decode_records


def _output_record(
    path: Path, root: Path, *, rows: int | None = None
) -> dict[str, Any]:
    record: dict[str, Any] = {
        "path": path.relative_to(root).as_posix(),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }
    if rows is not None:
        record["rows"] = int(rows)
    return record


def _result_contract(
    repo_root: Path,
    units: pd.DataFrame,
    scores: pd.DataFrame,
    decode_records: pd.DataFrame,
    *,
    decoder_bin: Path,
    cosi2_binary: Path,
) -> dict[str, Any]:
    selected = units[units["simulation_class"] == "selected"]
    neutral_units = units[units["simulation_class"] == "neutral"]
    implementation_paths = {
        "workflow": Path(__file__).resolve(),
        "neutral_simulation": Path(neutral.__file__).resolve(),
        "gamma_decode": Path(gamma_decode.__file__).resolve(),
        "gamma_analysis": Path(analysis.__file__).resolve(),
    }
    return {
        "schema": SCHEMA_VERSION,
        "scientific_contract": {
            "units": EXPECTED_TOTAL_UNITS,
            "demographies": list(DEMOGRAPHY_SETTINGS),
            "neutral_per_demography": EXPECTED_NEUTRAL_PER_DEMOGRAPHY,
            "selected_per_demography": EXPECTED_SELECTED_PER_DEMOGRAPHY,
            "neutral_conditioning": "present_population_survival_only_no_AF_gate",
            "selected_coefficient": SELECTION_COEFFICIENT,
            "thresholds_years": list(analysis.TMRCA_THRESHOLDS_YEARS),
            "statistics": list(analysis.SUPPORTED_STATISTICS),
            "windows": list(analysis.SUPPORTED_WINDOWS),
            "primary_window": analysis.PRIMARY_WINDOW,
            "paper_statistic": "paper_frac_posterior_mean_below",
            "decoder_settings": {
                "sequence_length_bp": gamma_decode.DEFAULT_SEQUENCE_LENGTH_BP,
                "focal_position_0based": gamma_decode.DEFAULT_FOCAL_POSITION_0BASED,
                "output_stride_bp": gamma_decode.DEFAULT_OUTPUT_STRIDE_BP,
                "local_half_width_bp": gamma_decode.DEFAULT_LOCAL_HALF_WIDTH_BP,
                "mutation_rate": gamma_decode.DEFAULT_MUTATION_RATE,
                "recombination_rate": gamma_decode.DEFAULT_RECOMBINATION_RATE,
                "recent_call": "mean",
            },
            "demography_settings": {
                name: {
                    "present_ne": settings["present_ne"],
                    "generation_time_years": settings["generation_time_years"],
                }
                for name, settings in DEMOGRAPHY_SETTINGS.items()
            },
            "pair_panel": {
                "selection": "fixed_SHA256_ranked_allele_blind",
                "seed": gamma_decode.DEFAULT_PAIR_PANEL_SEED,
                "haplotypes": gamma_decode.DEFAULT_PAIR_PANEL_HAPLOTYPES,
                "unordered_pairs": gamma_decode.DEFAULT_PAIR_COUNT,
            },
            "decoder_threads_per_unit": 1,
        },
        "input_payloads": {
            "unit_metadata_sha256": _portable_frame_sha256(units),
            "neutral_metadata_sha256": _portable_frame_sha256(neutral_units),
            "scores_sha256": _portable_frame_sha256(scores),
            "selected_run_completion_sha256": {
                str(row["demography"]): str(row["source_completion_sha256"])
                for row in selected.to_dict(orient="records")
            },
            "neutral_completion_sha256": {
                str(row["unit_id"]): str(row["source_completion_sha256"])
                for row in neutral_units.to_dict(orient="records")
            },
            "decode_completion_sha256": {
                str(row["unit_id"]): str(row["decode_completion_sha256"])
                for row in decode_records.to_dict(orient="records")
            },
            "decode_pair_panel_alt_count": {
                str(row["unit_id"]): int(row["decode_pair_panel_alt_count"])
                for row in decode_records.to_dict(orient="records")
            },
            "decode_focal_count_status": {
                str(row["unit_id"]): str(row["decode_focal_count_status"])
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
            for label, path in implementation_paths.items()
        },
    }


def _run_results_text(
    pointwise: pd.DataFrame,
    omnibus: pd.DataFrame,
    unit_metadata: pd.DataFrame,
) -> str:
    primary_pointwise = pointwise[
        (pointwise["primary_window"])
        & (pointwise["statistic"] == "paper_frac_posterior_mean_below")
    ].sort_values(["demography", "threshold_years"])
    primary_omnibus = omnibus[
        (omnibus["primary_window"])
        & (omnibus["statistic"] == "paper_frac_posterior_mean_below")
    ].sort_values("demography")
    lines = [
        "# CoSi2 Gamma-SMC selected versus natural-neutral results",
        "",
        "This bundle compares one accepted s=0.01 selected locus with 100 ",
        "natural-neutral surviving loci in each demography. Neutral loci were ",
        "not conditioned on present-day or sampled allele frequency.",
        "",
        "The primary statistic is the fraction of pairwise posterior-mean TMRCA ",
        "calls below each threshold across all 4,950 unordered pairs among a ",
        "fixed allele-blind SHA-256-ranked panel of 100 haplotypes, averaged over ",
        "the 21 positions spanning +/-100 kb around the focal site. The soft mean ",
        "posterior CDF and the single focal position are secondary.",
        "",
        "## Fixed pair-panel focal counts",
        "",
        "| Demography | Selected full sample | Selected fixed panel |",
        "|---|---:|---:|",
    ]
    for row in (
        unit_metadata[unit_metadata["simulation_class"].astype(str) == "selected"]
        .sort_values("demography")
        .to_dict(orient="records")
    ):
        lines.append(
            f"| {row['demography']} | {int(row['decode_full_sample_alt_count'])}/200 | "
            f"{int(row['decode_pair_panel_alt_count'])}/100 |"
        )
    lines.extend(
        [
            "",
            "## Primary pointwise simulation-reference p-values",
            "",
            "| Demography | Threshold (years) | Selected score | Upper p | Two-sided p | Descriptive AUC |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in primary_pointwise.to_dict(orient="records"):
        lines.append(
            "| {demography} | {threshold_years} | {observed_score:.6g} | "
            "{mc_p_upper:.6g} | {mc_p_two_sided:.6g} | "
            "{neutral_midrank_descriptive_auc:.6g} |".format(**row)
        )
    lines.extend(
        [
            "",
            "## Correlation-aware seven-threshold minP omnibus",
            "",
            "| Demography | Upper omnibus p | Two-sided omnibus p | Minimum attainable p |",
            "|---|---:|---:|---:|",
        ]
    )
    for row in primary_omnibus.to_dict(orient="records"):
        lines.append(
            "| {demography} | {omnibus_mc_p_upper:.6g} | "
            "{omnibus_mc_p_two_sided:.6g} | {minimum_attainable_p:.6g} |".format(**row)
        )
    lines.extend(
        [
            "",
            "These are simulation-reference model-contrast p-values, not p-values ",
            "conditional on a 20% endpoint AF. With one selected simulation per ",
            "demography, the reported midrank/AUC is a descriptive percentile and ",
            "must not be interpreted as a stable ROC power estimate.",
            "",
        ]
    )
    return "\n".join(lines)


def _validate_figure(path: Path, extension: str) -> None:
    if path.is_symlink() or not path.is_file() or path.stat().st_size < 8:
        raise ValueError(f"result figure is missing or empty: {path}")
    header = path.read_bytes()[:8]
    if extension == "png" and header != b"\x89PNG\r\n\x1a\n":
        raise ValueError("result PNG signature is invalid")
    if extension == "pdf" and not header.startswith(b"%PDF-"):
        raise ValueError("result PDF signature is invalid")


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
    """Create or strictly validate the atomic final analysis bundle."""

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
    if decode_counts["unit_id"].duplicated().any() or len(decode_counts) != len(
        inventory
    ):
        raise ValueError("decode panel-count inventory is incomplete")
    enriched_inventory = inventory.merge(
        decode_counts, on="unit_id", how="left", validate="one_to_one"
    )
    if (
        enriched_inventory[
            [
                "decode_full_sample_alt_count",
                "decode_pair_panel_alt_count",
                "decode_focal_count_status",
            ]
        ]
        .isna()
        .any(axis=None)
    ):
        raise ValueError("decode panel counts are absent after unit join")
    contract = _result_contract(
        root,
        enriched_inventory,
        score_frame,
        decode_records,
        decoder_bin=decoder,
        cosi2_binary=cosi2,
    )
    if destination.exists() and any(destination.iterdir()):
        verified = verify_results_bundle(
            root,
            results_dir=destination,
            expected_contract=contract,
        )
        return {**verified, "cache_hit": True}
    if destination.exists() and not destination.is_dir():
        raise ValueError("results destination exists but is not a directory")
    if destination.exists():
        destination.rmdir()
    destination.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(
        tempfile.mkdtemp(prefix=f".{destination.name}.staging.", dir=destination.parent)
    )
    try:
        derived = analysis.analyze_scores(score_frame)
        pointwise = derived["pointwise_pvalues"]
        omnibus = derived["minp_omnibus"]
        neutral_metadata = enriched_inventory[
            enriched_inventory["simulation_class"] == "neutral"
        ].reset_index(drop=True)
        _atomic_frame(stage / RESULT_OUTPUTS["scores"], score_frame)
        _atomic_frame(stage / RESULT_OUTPUTS["neutral_metadata"], neutral_metadata)
        _atomic_frame(stage / RESULT_OUTPUTS["unit_metadata"], enriched_inventory)
        _atomic_frame(stage / RESULT_OUTPUTS["pointwise_pvalues"], pointwise)
        _atomic_frame(stage / RESULT_OUTPUTS["minp_omnibus"], omnibus)
        _atomic_text(
            stage / RESULT_OUTPUTS["run_results"],
            _run_results_text(pointwise, omnibus, enriched_inventory),
        )
        figures = analysis.plot_score_profiles(
            score_frame,
            stage,
            stem=PROFILE_FIGURE_STEM,
        )
        expected_figures = {
            "png": stage / RESULT_OUTPUTS["profile_png"],
            "pdf": stage / RESULT_OUTPUTS["profile_pdf"],
        }
        if {key: Path(value) for key, value in figures.items()} != expected_figures:
            raise ValueError("analysis figure helper returned an unexpected inventory")
        _validate_figure(expected_figures["png"], "png")
        _validate_figure(expected_figures["pdf"], "pdf")
        row_counts = {
            "scores": len(score_frame),
            "neutral_metadata": len(neutral_metadata),
            "unit_metadata": len(inventory),
            "pointwise_pvalues": len(pointwise),
            "minp_omnibus": len(omnibus),
        }
        outputs = {
            label: _output_record(
                stage / filename,
                stage,
                rows=row_counts.get(label),
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
        expected_names = {RESULTS_COMPLETION_FILENAME, *RESULT_OUTPUTS.values()}
        if {path.name for path in stage.iterdir()} != expected_names:
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
    """Strict, read-only validation of the final exact-inventory bundle."""

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
        raise ValueError("results completion is stale for the current inputs")
    scientific = contract.get("scientific_contract", {})
    if (
        contract.get("schema") != SCHEMA_VERSION
        or scientific.get("units") != EXPECTED_TOTAL_UNITS
        or scientific.get("neutral_per_demography") != EXPECTED_NEUTRAL_PER_DEMOGRAPHY
        or scientific.get("selected_per_demography") != EXPECTED_SELECTED_PER_DEMOGRAPHY
        or scientific.get("pair_panel", {}).get("seed")
        != gamma_decode.DEFAULT_PAIR_PANEL_SEED
        or scientific.get("pair_panel", {}).get("haplotypes")
        != gamma_decode.DEFAULT_PAIR_PANEL_HAPLOTYPES
        or scientific.get("pair_panel", {}).get("unordered_pairs")
        != gamma_decode.DEFAULT_PAIR_COUNT
    ):
        raise ValueError("results scientific contract is incompatible")
    outputs = completion.get("outputs")
    if (
        not isinstance(outputs, dict)
        or set(outputs) != set(RESULT_OUTPUTS)
        or completion.get("output_count") != len(RESULT_OUTPUTS)
    ):
        raise ValueError("results output manifest changed")
    for label, filename in RESULT_OUTPUTS.items():
        record = outputs[label]
        path = destination / filename
        if (
            not isinstance(record, dict)
            or record.get("path") != filename
            or path.stat().st_size != int(record.get("size_bytes", -1))
            or sha256_file(path) != record.get("sha256")
        ):
            raise ValueError(f"results output checksum failed: {label}")
    scores = analysis.validate_scores(
        pd.read_csv(destination / RESULT_OUTPUTS["scores"], sep="\t")
    )
    persisted_unit_metadata = pd.read_csv(
        destination / RESULT_OUTPUTS["unit_metadata"], sep="\t"
    )
    units = validate_unit_inventory(persisted_unit_metadata)
    panel_columns = {
        "decode_full_sample_alt_count",
        "decode_pair_panel_alt_count",
        "decode_focal_count_status",
    }
    if not panel_columns.issubset(persisted_unit_metadata.columns):
        raise ValueError("persisted unit metadata lacks decode panel counts")
    full_counts = pd.to_numeric(
        persisted_unit_metadata["decode_full_sample_alt_count"], errors="raise"
    )
    panel_counts = pd.to_numeric(
        persisted_unit_metadata["decode_pair_panel_alt_count"], errors="raise"
    )
    if np.any((panel_counts < -1) | (panel_counts > 100)) or np.any(
        panel_counts != np.floor(panel_counts)
    ):
        raise ValueError("persisted decode panel counts are invalid")
    statuses = persisted_unit_metadata["decode_focal_count_status"].astype(str)
    allowed_statuses = {
        "absent_from_sample",
        "unique_exact_focal_site",
        "ambiguous_exact_coordinate_collision",
    }
    unambiguous = statuses != "ambiguous_exact_coordinate_collision"
    if (
        not set(statuses).issubset(allowed_statuses)
        or np.any(
            (panel_counts == -1) != statuses.eq("ambiguous_exact_coordinate_collision")
        )
        or not np.array_equal(
            full_counts[unambiguous].to_numpy(),
            units.loc[unambiguous, "sample_alt_count"].to_numpy(),
        )
    ):
        raise ValueError("persisted focal-count statuses are invalid")
    neutral_metadata = pd.read_csv(
        destination / RESULT_OUTPUTS["neutral_metadata"], sep="\t"
    )
    if len(neutral_metadata) != 2 * EXPECTED_NEUTRAL_PER_DEMOGRAPHY or set(
        neutral_metadata["simulation_class"].astype(str)
    ) != {"neutral"}:
        raise ValueError("persisted neutral metadata has the wrong inventory")
    expected_metadata = persisted_unit_metadata[
        persisted_unit_metadata["simulation_class"] == "neutral"
    ].reset_index(drop=True)
    pd.testing.assert_frame_equal(
        neutral_metadata,
        expected_metadata,
        check_dtype=False,
        check_exact=False,
        rtol=0,
        atol=1e-12,
    )
    input_payloads = contract.get("input_payloads", {})
    expected_payloads = {
        "scores_sha256": _portable_frame_sha256(scores),
        "unit_metadata_sha256": _portable_frame_sha256(persisted_unit_metadata),
        "neutral_metadata_sha256": _portable_frame_sha256(expected_metadata),
    }
    for label, expected in expected_payloads.items():
        if input_payloads.get(label) != expected:
            raise ValueError(f"results contract has stale {label}")
    expected_binding_counts = {
        "selected_run_completion_sha256": 2,
        "neutral_completion_sha256": 2 * EXPECTED_NEUTRAL_PER_DEMOGRAPHY,
        "decode_completion_sha256": EXPECTED_TOTAL_UNITS,
    }
    for label, expected_count in expected_binding_counts.items():
        binding = input_payloads.get(label)
        if not isinstance(binding, dict) or len(binding) != expected_count:
            raise ValueError(f"results contract has incomplete {label}")
    expected_row_counts = {
        "scores": len(scores),
        "unit_metadata": len(units),
        "neutral_metadata": len(neutral_metadata),
    }
    for label, expected in expected_row_counts.items():
        if int(outputs[label].get("rows", -1)) != expected:
            raise ValueError(f"results output has stale row count: {label}")
    implementation = contract.get("implementation", {})
    expected_modules = {
        "workflow": Path(__file__).resolve(),
        "neutral_simulation": Path(neutral.__file__).resolve(),
        "gamma_decode": Path(gamma_decode.__file__).resolve(),
        "gamma_analysis": Path(analysis.__file__).resolve(),
    }
    for label, module_path in expected_modules.items():
        record = implementation.get(label, {})
        if record.get("sha256") != sha256_file(module_path):
            raise ValueError(f"results implementation binding is stale: {label}")
    binaries = contract.get("binaries", {})
    if not isinstance(binaries, dict) or set(binaries) != {"cosi2", "gamma_smc"}:
        raise ValueError("results binary binding inventory changed")
    for label, record in binaries.items():
        recorded_path = Path(str(record.get("path", "")))
        binary_path = (
            recorded_path if recorded_path.is_absolute() else root / recorded_path
        )
        if not binary_path.is_file() or sha256_file(binary_path) != record.get(
            "sha256"
        ):
            raise ValueError(f"results binary binding is stale: {label}")
    derived = analysis.analyze_scores(scores)
    persisted_pointwise = pd.read_csv(
        destination / RESULT_OUTPUTS["pointwise_pvalues"], sep="\t"
    )
    persisted_omnibus = pd.read_csv(
        destination / RESULT_OUTPUTS["minp_omnibus"], sep="\t"
    )
    pd.testing.assert_frame_equal(
        persisted_pointwise,
        derived["pointwise_pvalues"],
        check_dtype=False,
        check_exact=False,
        rtol=0,
        atol=1e-12,
    )
    if int(outputs["pointwise_pvalues"].get("rows", -1)) != len(
        persisted_pointwise
    ) or int(outputs["minp_omnibus"].get("rows", -1)) != len(persisted_omnibus):
        raise ValueError("results inference table row count is stale")
    pd.testing.assert_frame_equal(
        persisted_omnibus,
        derived["minp_omnibus"],
        check_dtype=False,
        check_exact=False,
        rtol=0,
        atol=1e-12,
    )
    _validate_figure(destination / RESULT_OUTPUTS["profile_png"], "png")
    _validate_figure(destination / RESULT_OUTPUTS["profile_pdf"], "pdf")
    if (
        not (destination / RESULT_OUTPUTS["run_results"])
        .read_text(encoding="utf-8")
        .strip()
    ):
        raise ValueError("RUN_RESULTS.md is empty")
    return {**completion, "cache_hit": True}


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
        "plan": work / NEUTRAL_PLAN_DIRNAME,
        "neutral_runs": work / NEUTRAL_RUN_DIRNAME,
        "decode": work / DECODE_DIRNAME,
        "results": results,
    }


def plan_study(
    repo_root: str | Path,
    *,
    work_dir: str | Path | None = None,
) -> dict[str, Any]:
    root = Path(repo_root).resolve()
    paths = _resolve_paths(root, work_dir, None)
    selected = bind_selected_units(root)
    plan = neutral.generate_neutral_plan(root, output_dir=paths["plan"])
    return {
        "status": "planned",
        "selected_units": len(selected),
        "neutral_units": neutral.TOTAL_UNITS,
        "neutral_plan": str(paths["plan"]),
        "neutral_plan_status": plan.get("status"),
    }


def simulate_study(
    repo_root: str | Path,
    *,
    work_dir: str | Path | None = None,
    cosi2_binary: str | Path | None = None,
    max_workers: int = MAX_WORKERS,
    max_attempts: int = neutral.DEFAULT_COSI_MAXATTEMPTS,
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
    bind_selected_units(root)
    neutral.generate_neutral_plan(root, output_dir=paths["plan"])
    return neutral.run_neutral_simulations(
        root,
        cosi2_binary=binary,
        plan_dir=paths["plan"],
        run_dir=paths["neutral_runs"],
        max_workers=_positive_workers(max_workers),
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
        run_dir=paths["neutral_runs"],
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
        run_dir=paths["neutral_runs"],
        cosi2_binary=cosi2,
    )
    scores, decode_records = collect_decoded_scores(
        units,
        decode_root=paths["decode"],
        decoder_bin=decoder,
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
    """Return compact phase status without creating or modifying artifacts."""

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
    selected = bind_selected_units(root)
    result: dict[str, Any] = {
        "selected": {"complete": len(selected), "planned": 2},
        "neutral_plan": "absent",
        "neutral_simulations": {"pending": neutral.TOTAL_UNITS},
        "decode": {"pending": EXPECTED_TOTAL_UNITS},
        "results": "absent",
    }
    if paths["plan"].is_dir():
        neutral.verify_neutral_plan(root, plan_dir=paths["plan"])
        result["neutral_plan"] = "complete"
        status = neutral.neutral_simulation_status(
            root,
            plan_dir=paths["plan"],
            run_dir=paths["neutral_runs"],
            cosi2_binary=cosi2,
        )
        result["neutral_simulations"] = {
            str(key): int(value)
            for key, value in status["status"].value_counts().sort_index().items()
        }
        if (status["status"] == "complete").all():
            units = build_unit_inventory(
                root,
                plan_dir=paths["plan"],
                run_dir=paths["neutral_runs"],
                cosi2_binary=cosi2,
            )
            decodes = decode_status(
                units,
                decode_root=paths["decode"],
                decoder_bin=decoder,
            )
            result["decode"] = {
                str(key): int(value)
                for key, value in decodes["status"].value_counts().sort_index().items()
            }
    if paths["results"].exists():
        try:
            verify_results_bundle(root, results_dir=paths["results"])
        except Exception as error:  # noqa: BLE001 - status reports invalid artifacts
            result["results"] = f"invalid: {error}"
        else:
            result["results"] = "complete_internal_inventory"
    return result


def verify_study(
    repo_root: str | Path,
    *,
    work_dir: str | Path | None = None,
    results_dir: str | Path | None = None,
    cosi2_binary: str | Path | None = None,
    decoder_bin: str | Path | None = None,
) -> dict[str, Any]:
    """Read-only verification of selected, neutral, decode, and result phases."""

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
        run_dir=paths["neutral_runs"],
        cosi2_binary=cosi2,
    )
    scores, decode_records = collect_decoded_scores(
        units,
        decode_root=paths["decode"],
        decoder_bin=decoder,
    )
    enriched_units = units.merge(
        decode_records[
            [
                "unit_id",
                "decode_full_sample_alt_count",
                "decode_pair_panel_alt_count",
                "decode_focal_count_status",
            ]
        ],
        on="unit_id",
        how="left",
        validate="one_to_one",
    )
    contract = _result_contract(
        root,
        enriched_units,
        scores,
        decode_records,
        decoder_bin=decoder,
        cosi2_binary=cosi2,
    )
    results = verify_results_bundle(
        root,
        results_dir=paths["results"],
        expected_contract=contract,
    )
    return {
        "status": "complete",
        "selected_units": 2,
        "neutral_units": neutral.TOTAL_UNITS,
        "decoded_units": len(decode_records),
        "score_rows": len(scores),
        "results_contract_sha256": results["contract_sha256"],
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "phase",
        choices=("plan", "simulate", "status", "decode", "analyze", "verify", "all"),
    )
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--work-dir", type=Path)
    parser.add_argument("--results-dir", type=Path)
    parser.add_argument("--cosi2-binary", type=Path)
    parser.add_argument("--decoder-bin", type=Path)
    parser.add_argument("--workers", type=int, default=MAX_WORKERS)
    parser.add_argument(
        "--max-attempts", type=int, default=neutral.DEFAULT_COSI_MAXATTEMPTS
    )
    parser.add_argument("--timeout-seconds", type=float)
    parser.add_argument(
        "--unit-id",
        action="append",
        help="Run only a fixed predeclared neutral/decode unit; repeat as needed.",
    )
    return parser


def _print_frame_summary(label: str, frame: pd.DataFrame) -> None:
    counts = frame["status"].value_counts().sort_index().to_dict()
    print(json.dumps({label: counts}, sort_keys=True))


def _require_unit_phase_success(frame: pd.DataFrame, *, label: str) -> None:
    if frame.empty or not frame["status"].isin({"complete", "cached"}).all():
        counts = frame["status"].value_counts().sort_index().to_dict()
        raise RuntimeError(f"{label} phase did not complete: {counts}")


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    root = args.repo_root.resolve()
    common = {
        "work_dir": args.work_dir,
        "cosi2_binary": args.cosi2_binary,
    }
    if args.unit_id and args.phase in {"analyze", "verify", "all", "status", "plan"}:
        raise ValueError("--unit-id is supported only for simulate or decode")
    if args.phase == "plan":
        print(json.dumps(plan_study(root, work_dir=args.work_dir), indent=2))
    elif args.phase == "simulate":
        result = simulate_study(
            root,
            **common,
            max_workers=args.workers,
            max_attempts=args.max_attempts,
            timeout_seconds=args.timeout_seconds,
            unit_ids=args.unit_id,
        )
        _print_frame_summary("simulate", result)
        _require_unit_phase_success(result, label="simulate")
    elif args.phase == "status":
        print(
            json.dumps(
                workflow_status(
                    root,
                    **common,
                    results_dir=args.results_dir,
                    decoder_bin=args.decoder_bin,
                ),
                indent=2,
                sort_keys=True,
            )
        )
    elif args.phase == "decode":
        result = decode_study(
            root,
            **common,
            decoder_bin=args.decoder_bin,
            max_workers=args.workers,
            unit_ids=args.unit_id,
        )
        _print_frame_summary("decode", result)
    elif args.phase == "analyze":
        result = analyze_study(
            root,
            **common,
            results_dir=args.results_dir,
            decoder_bin=args.decoder_bin,
        )
        print(
            json.dumps({"status": result["status"], "cache_hit": result["cache_hit"]})
        )
    elif args.phase == "verify":
        print(
            json.dumps(
                verify_study(
                    root,
                    **common,
                    results_dir=args.results_dir,
                    decoder_bin=args.decoder_bin,
                ),
                indent=2,
                sort_keys=True,
            )
        )
    else:
        plan_study(root, work_dir=args.work_dir)
        simulations = simulate_study(
            root,
            **common,
            max_workers=args.workers,
            max_attempts=args.max_attempts,
            timeout_seconds=args.timeout_seconds,
        )
        _require_unit_phase_success(simulations, label="simulate")
        decode_study(
            root,
            **common,
            decoder_bin=args.decoder_bin,
            max_workers=args.workers,
        )
        analyze_study(
            root,
            **common,
            results_dir=args.results_dir,
            decoder_bin=args.decoder_bin,
        )
        print(
            json.dumps(
                verify_study(
                    root,
                    **common,
                    results_dir=args.results_dir,
                    decoder_bin=args.decoder_bin,
                ),
                indent=2,
                sort_keys=True,
            )
        )
    return 0


__all__ = [
    "DEFAULT_RESULTS_RELATIVE",
    "DEFAULT_WORK_RELATIVE",
    "DEMOGRAPHY_SETTINGS",
    "EXPECTED_TOTAL_UNITS",
    "MAX_WORKERS",
    "RESULTS_COMPLETION_FILENAME",
    "RESULT_OUTPUTS",
    "analyze_study",
    "bind_selected_units",
    "build_parser",
    "build_unit_inventory",
    "collect_decoded_scores",
    "decode_status",
    "decode_study",
    "main",
    "plan_study",
    "run_decodes",
    "simulate_study",
    "validate_unit_inventory",
    "verify_decode_bank",
    "verify_decode_unit",
    "verify_results_bundle",
    "verify_study",
    "workflow_status",
    "write_results_bundle",
]
