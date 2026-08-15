"""Restartable workflow for the natural-AF Han ``s=0.001`` TMRCA study.

Simulation, decoding, and score collection are owned by
``han_s001_tmrca_simulation``.  This module is the additive integration layer:
it dispatches those phases and atomically publishes the complete statistical
analysis bundle.  A cached bundle is accepted only when its plan, verified
simulation/decode inputs, source files, parameters, exact file inventory, and
all output checksums still match.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import shutil
import tempfile
from collections.abc import Mapping, Sequence
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

import pandas as pd

from . import eas_sweep_models, eas_sweep_study
from . import han_s001_tmrca_analysis as analysis
from . import han_s001_tmrca_simulation as simulation


BASE_SEED = eas_sweep_study.BASE_SEED


SCHEMA_VERSION = "gamma-smc.han-s001-tmrca-workflow/v1"
DEFAULT_NEUTRAL_REPLICATES = 100
DEFAULT_SELECTED_REPLICATES = 100
DEFAULT_WORKERS = 20
DEFAULT_THREADS = 1
DEFAULT_ANALYSIS_DIR = Path("focused_selection_EAS_sim/results/han_s001_tmrca/analysis")
COMPLETION_FILENAME = "analysis_completion.json"
REPLICATE_SCORES_FILENAME = "han_s001_replicate_scores.tsv"
SELECTED_NEUTRAL_FIGURE_STEM = "han_s001_selected_vs_neutral_tmrca_curves"
AUC_FIGURE_STEM = "han_s001_tmrca_roc_auc_profiles"


def _sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _jsonable(value: Any, *, repo_root: Path | None = None) -> Any:
    """Convert provenance objects to deterministic strict-JSON values."""

    if isinstance(value, Path):
        resolved = value.resolve()
        if repo_root is not None:
            try:
                return resolved.relative_to(repo_root).as_posix()
            except ValueError:
                pass
        return resolved.as_posix()
    if isinstance(value, pd.DataFrame):
        return {
            "columns": [str(column) for column in value.columns],
            "records": [
                {
                    str(key): _jsonable(item, repo_root=repo_root)
                    for key, item in record.items()
                }
                for record in value.to_dict(orient="records")
            ],
        }
    if isinstance(value, Mapping):
        return {
            str(key): _jsonable(item, repo_root=repo_root)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, (list, tuple)):
        return [_jsonable(item, repo_root=repo_root) for item in value]
    if isinstance(value, (set, frozenset)):
        converted = [_jsonable(item, repo_root=repo_root) for item in value]
        return sorted(converted, key=_canonical_json)
    if hasattr(value, "item") and callable(value.item):
        return _jsonable(value.item(), repo_root=repo_root)
    if isinstance(value, float):
        if math.isnan(value):
            return {"__nonfinite_float__": "nan"}
        if math.isinf(value):
            return {"__nonfinite_float__": "inf" if value > 0 else "-inf"}
        return value
    if value is None or isinstance(value, (str, int, bool)):
        return value
    raise TypeError(f"unsupported provenance value: {type(value).__name__}")


def _canonical_json(value: Any) -> str:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    )


def _canonical_sha256(value: Any) -> str:
    return hashlib.sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def _atomic_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        temporary.write_text(
            json.dumps(
                payload,
                indent=2,
                sort_keys=True,
                ensure_ascii=True,
                allow_nan=False,
            )
            + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        frame.to_csv(temporary, sep="\t", index=False, lineterminator="\n")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _positive_integer(value: Any, *, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{label} must be a positive integer")
    return value


def _nonnegative_integer(value: Any, *, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{label} must be a nonnegative integer")
    return value


def _resolve_analysis_dir(repo_root: Path, output_dir: str | Path | None) -> Path:
    if output_dir is not None:
        candidate = Path(output_dir)
        return (
            (repo_root / candidate).resolve()
            if not candidate.is_absolute()
            else candidate.resolve()
        )
    return (repo_root / DEFAULT_ANALYSIS_DIR).resolve()


def _portable_path(path: Path, repo_root: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(repo_root).as_posix()
    except ValueError:
        return resolved.as_posix()


def _source_record(module: Any, repo_root: Path) -> dict[str, Any]:
    source = Path(module.__file__).resolve()
    if source.suffix == ".pyc" and source.with_suffix(".py").is_file():
        source = source.with_suffix(".py")
    if not source.is_file():
        raise ValueError(f"source file is absent: {source}")
    try:
        normalized = (
            source.read_text(encoding="utf-8").replace("\r\n", "\n").encode("utf-8")
        )
    except (OSError, UnicodeDecodeError) as error:
        raise ValueError(
            f"source file is not canonical UTF-8 text: {source}"
        ) from error
    return {
        "path": _portable_path(source, repo_root),
        "sha256_lf_utf8": hashlib.sha256(normalized).hexdigest(),
        "normalized_size_bytes": len(normalized),
        "hash_semantics": "UTF-8 text with CRLF normalized to LF",
    }


def _implementation_contract(repo_root: Path) -> dict[str, Any]:
    import sys

    package_versions: dict[str, str] = {}
    for package in ("numpy", "pandas", "scipy", "matplotlib"):
        try:
            package_versions[package] = version(package)
        except PackageNotFoundError as error:
            raise ValueError(
                f"required analysis package is absent: {package}"
            ) from error
    return {
        "workflow": _source_record(sys.modules[__name__], repo_root),
        "simulation": _source_record(simulation, repo_root),
        "analysis": _source_record(analysis, repo_root),
        "eas_sweep_models": _source_record(eas_sweep_models, repo_root),
        "eas_sweep_study": _source_record(eas_sweep_study, repo_root),
        "software_versions": {
            "python": platform.python_version(),
            **package_versions,
        },
    }


def _output_specs() -> dict[str, str]:
    specs = {"replicate_scores": REPLICATE_SCORES_FILENAME}
    specs.update(
        {
            f"table_{label}": filename
            for label, filename in analysis.TABLE_FILENAMES.items()
        }
    )
    specs.update(
        {
            "figure_selected_neutral_png": f"{SELECTED_NEUTRAL_FIGURE_STEM}.png",
            "figure_selected_neutral_pdf": f"{SELECTED_NEUTRAL_FIGURE_STEM}.pdf",
            "figure_auc_png": f"{AUC_FIGURE_STEM}.png",
            "figure_auc_pdf": f"{AUC_FIGURE_STEM}.pdf",
        }
    )
    if len(set(specs.values())) != len(specs):
        raise RuntimeError("analysis output filenames are not unique")
    for relative in specs.values():
        path = Path(relative)
        if path.is_absolute() or ".." in path.parts or len(path.parts) != 1:
            raise RuntimeError(f"unsafe analysis output path: {relative}")
    return specs


def _frame_payload_sha256(frame: pd.DataFrame) -> str:
    payload = frame.to_csv(None, sep="\t", index=False, lineterminator="\n").encode(
        "utf-8"
    )
    return hashlib.sha256(payload).hexdigest()


def _plan_semantics(loaded: Any, repo_root: Path) -> Any:
    if not isinstance(loaded, (Mapping, pd.DataFrame)):
        raise TypeError("simulation load_plan must return a mapping or DataFrame")
    return _jsonable(loaded, repo_root=repo_root)


def _plan_counts(loaded: Any) -> dict[str, int]:
    if isinstance(loaded, Mapping):
        try:
            return {
                "neutral": _positive_integer(loaded["n_neutral"], label="n_neutral"),
                "selected": _positive_integer(loaded["n_selected"], label="n_selected"),
            }
        except KeyError as error:
            raise ValueError("study plan lacks replicate counts") from error
    if isinstance(loaded, pd.DataFrame) and {
        "unit_id",
        "simulation_class",
    }.issubset(loaded.columns):
        units = loaded.loc[:, ["unit_id", "simulation_class"]].drop_duplicates()
        if units["unit_id"].duplicated(keep=False).any():
            raise ValueError("study plan changes a unit simulation class")
        counts = units.groupby("simulation_class").size().to_dict()
        if set(counts) != {"neutral", "selected"}:
            raise ValueError("study plan lacks a simulation class")
        return {
            key: _positive_integer(int(value), label=f"n_{key}")
            for key, value in counts.items()
        }
    raise TypeError("simulation load_plan must expose planned replicate counts")


def _load_durable_replicate_scores(repo_root: Path) -> pd.DataFrame:
    """Load the canonical score payload published by the simulation workflow."""

    relative_scores = getattr(simulation, "SCORES_RELATIVE_PATH", None)
    if relative_scores is None:
        raise ValueError("simulation module does not expose its verified score path")
    return analysis.load_replicate_scores(repo_root / Path(relative_scores))


def _empirical_input(
    empirical_path: str | Path | None, repo_root: Path
) -> tuple[pd.DataFrame | None, dict[str, Any]]:
    if empirical_path is None:
        return None, {
            "mode": "template_only_no_empirical_observation",
            "path": None,
            "sha256": None,
            "claim_status": "no_real_empirical_p_value_claimed",
        }
    candidate = Path(empirical_path)
    path = (
        (repo_root / candidate).resolve()
        if not candidate.is_absolute()
        else candidate.resolve()
    )
    if not path.is_file():
        raise ValueError("requested explicit empirical input is absent")
    try:
        empirical = pd.read_csv(path, sep="\t")
    except Exception as error:
        raise ValueError("requested explicit empirical input is unreadable") from error
    validated = analysis.validate_empirical_input(empirical)
    return validated, {
        "mode": "explicit_observation_with_self_declared_provenance_fields",
        "path": _portable_path(path, repo_root),
        "sha256": _sha256_file(path),
        "rows": len(validated),
        "columns": list(validated.columns),
        "claim_status": (
            "empirical_p_values_available_for_explicit_input; provenance_hash_"
            "fields_are_format_validated_but_not_independently_resolved"
        ),
    }


def _analysis_inputs(
    repo_root: Path,
    *,
    n_neutral: int | None,
    n_selected: int | None,
    empirical_path: str | Path | None,
    collect_scores: bool,
) -> tuple[pd.DataFrame, pd.DataFrame | None, dict[str, Any]]:
    loaded_plan = simulation.load_plan(repo_root)
    plan_path = Path(simulation.plan_path(repo_root)).resolve()
    if not plan_path.is_file():
        raise ValueError("validated simulation plan path is absent")
    planned_counts = _plan_counts(loaded_plan)
    requested_counts = {"neutral": n_neutral, "selected": n_selected}
    for simulation_class, requested in requested_counts.items():
        if requested is not None and requested != planned_counts[simulation_class]:
            raise ValueError(
                f"requested {simulation_class} replicate count differs from the "
                "existing immutable study plan"
            )
    verification = simulation.verify_study(repo_root)
    if not isinstance(verification, Mapping):
        raise TypeError("simulation verify_study must return a mapping")
    if verification.get("status") != "complete":
        raise ValueError("simulation/decode verification is not complete")
    if collect_scores:
        simulation.collect_replicate_scores(repo_root)
    scores = _load_durable_replicate_scores(repo_root)
    unit_classes = scores.loc[:, ["unit_id", "simulation_class"]].drop_duplicates()
    if unit_classes["unit_id"].duplicated(keep=False).any():
        raise ValueError("collected unit changes simulation class")
    observed_counts = unit_classes.groupby("simulation_class").size().to_dict()
    if observed_counts != planned_counts:
        raise ValueError(
            "collected replicate counts do not match the immutable plan: "
            f"expected={planned_counts}, observed={observed_counts}"
        )
    empirical, empirical_record = _empirical_input(empirical_path, repo_root)
    normalized_verification = _jsonable(verification, repo_root=repo_root)
    inputs = {
        "plan": {
            "path": _portable_path(plan_path, repo_root),
            "sha256": _sha256_file(plan_path),
            "size_bytes": plan_path.stat().st_size,
            "semantic_sha256": _canonical_sha256(
                _plan_semantics(loaded_plan, repo_root)
            ),
            "replicate_counts": planned_counts,
        },
        "verified_simulation_decode_completions": normalized_verification,
        "verified_simulation_decode_completions_sha256": _canonical_sha256(
            normalized_verification
        ),
        "replicate_score_payload_sha256": _frame_payload_sha256(scores),
        "replicate_score_rows": len(scores),
        "replicate_score_columns": list(scores.columns),
        "empirical": empirical_record,
    }
    return scores, empirical, inputs


def _analysis_contract(
    repo_root: Path,
    inputs: Mapping[str, Any],
    *,
    bootstrap_draws: int,
    permutation_draws: int,
    base_seed: int,
    include_af_stratified: bool,
    af_stratum_width: float,
    empirical_af_match_tolerance: float,
) -> dict[str, Any]:
    return {
        "schema": SCHEMA_VERSION,
        "study": "Han donor-fixed introgression; population survival only",
        "comparison": "selected s=0.001 versus natural neutral introgressed loci",
        "terminal_frequency_conditioning": False,
        "score": "overall-pair mean P(TMRCA < x)",
        "thresholds_years": list(analysis.THRESHOLDS_YEARS),
        "inputs": dict(inputs),
        "implementation": _implementation_contract(repo_root),
        "parameters": {
            "bootstrap_draws": bootstrap_draws,
            "permutation_draws": permutation_draws,
            "base_seed": base_seed,
            "include_af_stratified": bool(include_af_stratified),
            "af_stratum_width": float(af_stratum_width),
            "empirical_af_match_tolerance": float(empirical_af_match_tolerance),
            "generation_time_years": float(eas_sweep_models.GENERATION_TIME_YEARS),
            "selected_selection_coefficient": float(
                analysis.EXPECTED_SELECTION_COEFFICIENT
            ),
            "auc_direction": "selected_score_greater_than_neutral_score_no_flipping",
            "monte_carlo_pvalue": "one_sided_upper_tail_plus_one",
        },
        "expected_outputs": _output_specs(),
        "empirical_claim_guard": (
            "No real empirical p value is claimed unless an explicit empirical TSV "
            "is supplied. Its schema and self-declared provenance hash fields are "
            "validated, but referenced artifacts are not independently resolved."
        ),
    }


def _read_table(path: Path) -> pd.DataFrame:
    try:
        return pd.read_csv(path, sep="\t")
    except Exception as error:
        raise ValueError(f"analysis table is unreadable: {path.name}") from error


def _output_record(
    path: Path, root: Path, frame: pd.DataFrame | None
) -> dict[str, Any]:
    if not path.is_file() or path.stat().st_size < 1:
        raise ValueError(f"analysis output is absent or empty: {path.name}")
    record: dict[str, Any] = {
        "path": path.relative_to(root).as_posix(),
        "sha256": _sha256_file(path),
        "size_bytes": path.stat().st_size,
    }
    if frame is not None:
        record["rows"] = len(frame)
        record["columns"] = list(frame.columns)
    return record


def _validate_figure_signature(path: Path) -> None:
    with path.open("rb") as handle:
        header = handle.read(8)
        handle.seek(max(0, path.stat().st_size - 2_048))
        trailer = handle.read()
    if path.suffix == ".png":
        if header != b"\x89PNG\r\n\x1a\n" or b"IEND" not in trailer:
            raise ValueError(f"analysis PNG is structurally invalid: {path.name}")
    elif path.suffix == ".pdf":
        if not header.startswith(b"%PDF-") or b"%%EOF" not in trailer:
            raise ValueError(f"analysis PDF is structurally invalid: {path.name}")
    else:
        raise ValueError(f"unexpected figure extension: {path.name}")


def _validate_analysis_bundle(
    output_dir: Path,
    *,
    expected_contract: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    completion_path = output_dir / COMPLETION_FILENAME
    if not completion_path.is_file():
        raise ValueError("analysis completion manifest is absent")
    try:
        completion = json.loads(completion_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError("analysis completion manifest is unreadable") from error
    if completion.get("schema") != SCHEMA_VERSION:
        raise ValueError("analysis completion schema is incompatible")
    if completion.get("status") != "complete":
        raise ValueError("analysis completion status is not complete")
    if set(completion) != {
        "schema",
        "status",
        "contract",
        "contract_sha256",
        "outputs",
    }:
        raise ValueError("analysis completion manifest has missing or extra fields")
    contract = completion.get("contract")
    if not isinstance(contract, Mapping):
        raise ValueError("analysis completion contract is malformed")
    if completion.get("contract_sha256") != _canonical_sha256(contract):
        raise ValueError("analysis completion contract checksum is invalid")
    if expected_contract is not None and _canonical_sha256(
        contract
    ) != _canonical_sha256(expected_contract):
        raise ValueError("cached analysis contract is incompatible with current inputs")

    expected_specs = _output_specs()
    if contract.get("expected_outputs") != expected_specs:
        raise ValueError("analysis contract output inventory is incompatible")
    outputs = completion.get("outputs")
    if not isinstance(outputs, Mapping) or set(outputs) != set(expected_specs):
        raise ValueError("analysis output manifest has missing or extra labels")

    expected_files = {COMPLETION_FILENAME, *expected_specs.values()}
    actual_entries: set[str] = set()
    for entry in output_dir.rglob("*"):
        relative = entry.relative_to(output_dir).as_posix()
        if entry.is_symlink():
            raise ValueError(f"analysis bundle contains a symlink: {relative}")
        if entry.is_dir():
            raise ValueError(
                f"analysis bundle contains an unexpected directory: {relative}"
            )
        actual_entries.add(relative)
    if actual_entries != expected_files:
        missing = sorted(expected_files.difference(actual_entries))
        extra = sorted(actual_entries.difference(expected_files))
        raise ValueError(
            "analysis file inventory mismatch"
            + (f"; missing={missing}" if missing else "")
            + (f"; extra={extra}" if extra else "")
        )

    for label, expected_relative in expected_specs.items():
        record = outputs[label]
        if not isinstance(record, Mapping):
            raise ValueError(f"analysis output record is malformed: {label}")
        expected_record_fields = {"path", "sha256", "size_bytes"}
        if label == "replicate_scores" or label.startswith("table_"):
            expected_record_fields.update({"rows", "columns"})
        if set(record) != expected_record_fields:
            raise ValueError(f"analysis output record fields are invalid: {label}")
        if record.get("path") != expected_relative:
            raise ValueError(f"analysis output path is incompatible: {label}")
        path = output_dir / expected_relative
        if _sha256_file(path) != record.get("sha256"):
            raise ValueError(f"analysis output checksum failed: {label}")
        if path.stat().st_size != record.get("size_bytes"):
            raise ValueError(f"analysis output size failed: {label}")
        if label == "replicate_scores" or label.startswith("table_"):
            frame = _read_table(path)
            if len(frame) != record.get("rows") or list(frame.columns) != record.get(
                "columns"
            ):
                raise ValueError(f"analysis output table contract failed: {label}")
        elif label.startswith("figure_"):
            _validate_figure_signature(path)
    return {
        **completion,
        "completion": completion_path,
        "cache_hit": expected_contract is not None,
    }


def _remove_owned_staging(stage: Path, parent: Path) -> None:
    resolved_stage = stage.resolve()
    resolved_parent = parent.resolve()
    if resolved_stage.parent != resolved_parent or not resolved_stage.name.startswith(
        ".han_s001_tmrca_analysis.staging."
    ):
        raise RuntimeError("refusing to remove an unowned analysis staging path")
    if resolved_stage.exists():
        shutil.rmtree(resolved_stage)


def analyze_study(
    repo_root: str | Path,
    *,
    output_dir: str | Path | None = None,
    empirical_path: str | Path | None = None,
    n_neutral: int | None = None,
    n_selected: int | None = None,
    bootstrap_draws: int = analysis.DEFAULT_BOOTSTRAP_DRAWS,
    permutation_draws: int = analysis.DEFAULT_PERMUTATION_DRAWS,
    base_seed: int = BASE_SEED,
    include_af_stratified: bool = True,
    af_stratum_width: float = analysis.DEFAULT_AF_STRATUM_WIDTH,
    empirical_af_match_tolerance: float = 0.025,
) -> dict[str, Any]:
    """Analyze collected scores and transactionally publish the exact bundle."""

    root = Path(repo_root).resolve()
    if n_neutral is not None:
        n_neutral = _positive_integer(n_neutral, label="n_neutral")
    if n_selected is not None:
        n_selected = _positive_integer(n_selected, label="n_selected")
    bootstrap_draws = _positive_integer(bootstrap_draws, label="bootstrap_draws")
    permutation_draws = _positive_integer(permutation_draws, label="permutation_draws")
    base_seed = _nonnegative_integer(base_seed, label="base_seed")
    if include_af_stratified is not True:
        raise ValueError(
            "include_af_stratified must remain true for the fixed all-table contract"
        )
    destination = _resolve_analysis_dir(root, output_dir)
    scores, empirical, inputs = _analysis_inputs(
        root,
        n_neutral=n_neutral,
        n_selected=n_selected,
        empirical_path=empirical_path,
        collect_scores=True,
    )
    contract = _analysis_contract(
        root,
        inputs,
        bootstrap_draws=bootstrap_draws,
        permutation_draws=permutation_draws,
        base_seed=base_seed,
        include_af_stratified=include_af_stratified,
        af_stratum_width=af_stratum_width,
        empirical_af_match_tolerance=empirical_af_match_tolerance,
    )
    if (destination / COMPLETION_FILENAME).is_file():
        return _validate_analysis_bundle(destination, expected_contract=contract)
    if destination.exists() and any(destination.iterdir()):
        raise ValueError(
            "analysis directory is nonempty without a valid completion manifest; "
            "quarantine the partial directory before retrying"
        )

    results = analysis.analyze_replicate_scores(
        scores,
        empirical=empirical,
        bootstrap_draws=bootstrap_draws,
        permutation_draws=permutation_draws,
        base_seed=base_seed,
        include_af_stratified=include_af_stratified,
        af_stratum_width=af_stratum_width,
        empirical_af_match_tolerance=empirical_af_match_tolerance,
    )
    parent = destination.parent
    parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        destination.rmdir()
    stage = Path(
        tempfile.mkdtemp(prefix=".han_s001_tmrca_analysis.staging.", dir=parent)
    )
    try:
        _atomic_tsv(scores, stage / REPLICATE_SCORES_FILENAME)
        table_paths = analysis.write_analysis_tables(results, stage)
        selected_neutral_paths = analysis.plot_selected_neutral_curves(
            scores,
            stage,
            bootstrap_draws=bootstrap_draws,
            base_seed=base_seed,
            stem=SELECTED_NEUTRAL_FIGURE_STEM,
        )
        auc_paths = analysis.plot_auc_profiles(
            results["auc_by_threshold"],
            results["af_only_benchmarks"],
            stage,
            stratified_auc=results["sample_af_stratified_auc"],
            stem=AUC_FIGURE_STEM,
        )
        expected_pairs = (
            (
                stage / f"{SELECTED_NEUTRAL_FIGURE_STEM}.png",
                stage / f"{SELECTED_NEUTRAL_FIGURE_STEM}.pdf",
            ),
            (
                stage / f"{AUC_FIGURE_STEM}.png",
                stage / f"{AUC_FIGURE_STEM}.pdf",
            ),
        )
        if tuple(map(Path, selected_neutral_paths)) != expected_pairs[0]:
            raise ValueError("selected-versus-neutral plot returned unexpected paths")
        if tuple(map(Path, auc_paths)) != expected_pairs[1]:
            raise ValueError("AUC plot returned unexpected paths")
        if {
            label: path.name for label, path in table_paths.items()
        } != analysis.TABLE_FILENAMES:
            raise ValueError("analysis table writer returned an incompatible inventory")

        frames = {"replicate_scores": scores}
        frames.update({f"table_{key}": value for key, value in results.items()})
        outputs = {
            label: _output_record(
                stage / relative,
                stage,
                frames.get(label),
            )
            for label, relative in _output_specs().items()
        }
        completion = {
            "schema": SCHEMA_VERSION,
            "status": "complete",
            "contract": contract,
            "contract_sha256": _canonical_sha256(contract),
            "outputs": outputs,
        }
        _atomic_json(stage / COMPLETION_FILENAME, completion)
        _validate_analysis_bundle(stage, expected_contract=contract)
        if destination.exists():
            raise ValueError("analysis destination appeared during atomic publication")
        os.replace(stage, destination)
        published = _validate_analysis_bundle(destination, expected_contract=contract)
        published["cache_hit"] = False
        return published
    finally:
        if stage.exists():
            _remove_owned_staging(stage, parent)


def verify_analysis(
    repo_root: str | Path,
    *,
    output_dir: str | Path | None = None,
    empirical_path: str | Path | None = None,
    n_neutral: int | None = None,
    n_selected: int | None = None,
    bootstrap_draws: int = analysis.DEFAULT_BOOTSTRAP_DRAWS,
    permutation_draws: int = analysis.DEFAULT_PERMUTATION_DRAWS,
    base_seed: int = BASE_SEED,
    include_af_stratified: bool = True,
    af_stratum_width: float = analysis.DEFAULT_AF_STRATUM_WIDTH,
    empirical_af_match_tolerance: float = 0.025,
) -> dict[str, Any]:
    """Rebuild the current contract and strictly verify a cached analysis."""

    root = Path(repo_root).resolve()
    if include_af_stratified is not True:
        raise ValueError(
            "include_af_stratified must remain true for the fixed all-table contract"
        )
    scores, _empirical, inputs = _analysis_inputs(
        root,
        n_neutral=(
            _positive_integer(n_neutral, label="n_neutral")
            if n_neutral is not None
            else None
        ),
        n_selected=(
            _positive_integer(n_selected, label="n_selected")
            if n_selected is not None
            else None
        ),
        empirical_path=empirical_path,
        collect_scores=False,
    )
    del scores
    contract = _analysis_contract(
        root,
        inputs,
        bootstrap_draws=_positive_integer(bootstrap_draws, label="bootstrap_draws"),
        permutation_draws=_positive_integer(
            permutation_draws, label="permutation_draws"
        ),
        base_seed=_nonnegative_integer(base_seed, label="base_seed"),
        include_af_stratified=include_af_stratified,
        af_stratum_width=af_stratum_width,
        empirical_af_match_tolerance=empirical_af_match_tolerance,
    )
    return _validate_analysis_bundle(
        _resolve_analysis_dir(root, output_dir), expected_contract=contract
    )


def status_tables(repo_root: str | Path) -> dict[str, pd.DataFrame]:
    root = Path(repo_root).resolve()
    return {
        "simulation": simulation.simulation_status(root),
        "decode": simulation.decode_status(root),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Natural-AF Han s=0.001 full-locus TMRCA workflow"
    )
    parser.add_argument(
        "phase",
        choices=("plan", "simulate", "status", "decode", "analyze", "verify", "all"),
    )
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--analysis-dir", type=Path)
    parser.add_argument("--neutral-replicates", type=int)
    parser.add_argument("--selected-replicates", type=int)
    parser.add_argument("--base-seed", type=int, default=BASE_SEED)
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument("--threads", type=int, default=DEFAULT_THREADS)
    parser.add_argument("--slim-bin", type=Path)
    parser.add_argument("--decoder-bin", type=Path)
    parser.add_argument("--empirical-tsv", type=Path)
    parser.add_argument(
        "--bootstrap-draws", type=int, default=analysis.DEFAULT_BOOTSTRAP_DRAWS
    )
    parser.add_argument(
        "--permutation-draws", type=int, default=analysis.DEFAULT_PERMUTATION_DRAWS
    )
    parser.add_argument(
        "--af-stratum-width", type=float, default=analysis.DEFAULT_AF_STRATUM_WIDTH
    )
    parser.add_argument("--empirical-af-match-tolerance", type=float, default=0.025)
    parser.add_argument("--force-plan", action="store_true")
    return parser


def _print_status(label: str, frame: pd.DataFrame) -> None:
    print(f"[{label}]")
    if "status" in frame.columns:
        print(frame.groupby("status", dropna=False).size().to_string())
    else:
        print(frame.to_string(index=False))


def _require_complete_phase(frame: pd.DataFrame, *, label: str) -> None:
    if not isinstance(frame, pd.DataFrame) or frame.empty:
        raise ValueError(f"{label} returned no status rows")
    if "status" not in frame.columns:
        raise ValueError(f"{label} status table lacks a status column")
    statuses = set(frame["status"].astype(str))
    incomplete = sorted(statuses.difference({"complete", "cached"}))
    if incomplete:
        raise ValueError(f"{label} did not complete successfully: {incomplete}")


def _preflight_analysis_arguments(args: argparse.Namespace, root: Path) -> None:
    _positive_integer(args.bootstrap_draws, label="bootstrap_draws")
    _positive_integer(args.permutation_draws, label="permutation_draws")
    _nonnegative_integer(args.base_seed, label="base_seed")
    width = float(args.af_stratum_width)
    if not math.isfinite(width) or width <= 0 or width > 1:
        raise ValueError("af_stratum_width must evenly divide [0, 1]")
    reciprocal = 1.0 / width
    if not math.isclose(reciprocal, round(reciprocal), rel_tol=0.0, abs_tol=1e-10):
        raise ValueError("af_stratum_width must evenly divide [0, 1]")
    tolerance = float(args.empirical_af_match_tolerance)
    if not math.isfinite(tolerance) or tolerance < 0 or tolerance > 0.05:
        raise ValueError("empirical_af_match_tolerance must lie in [0, 0.05]")
    if args.empirical_tsv is not None:
        _empirical_input(args.empirical_tsv, root)


def _validate_requested_counts_against_plan(
    loaded_plan: Any, *, n_neutral: int | None, n_selected: int | None
) -> None:
    planned = _plan_counts(loaded_plan)
    requested = {"neutral": n_neutral, "selected": n_selected}
    for simulation_class, value in requested.items():
        if value is not None and value != planned[simulation_class]:
            raise ValueError(
                f"requested {simulation_class} replicate count differs from the "
                "existing immutable study plan"
            )


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    root = args.repo_root.resolve()
    if args.force_plan and args.phase != "plan":
        raise ValueError("--force-plan is allowed only with the plan phase")
    if args.neutral_replicates is not None:
        _positive_integer(args.neutral_replicates, label="neutral_replicates")
    if args.selected_replicates is not None:
        _positive_integer(args.selected_replicates, label="selected_replicates")
    _positive_integer(args.workers, label="workers")
    _positive_integer(args.threads, label="threads")
    if args.phase in {"analyze", "verify", "all"}:
        _preflight_analysis_arguments(args, root)
    slim_bin = (args.slim_bin or root / ".native-stdpopsim/bin/slim").resolve()
    decoder_bin = (args.decoder_bin or root / "bin/gamma_smc").resolve()
    if args.phase in {"simulate", "all"} and not slim_bin.is_file():
        raise ValueError(f"SLiM binary is absent: {slim_bin}")
    if args.phase in {"decode", "all"} and not decoder_bin.is_file():
        raise ValueError(f"Gamma-SMC binary is absent: {decoder_bin}")

    existing_plan_path = Path(simulation.plan_path(root)).resolve()
    should_create_plan = args.phase in {"plan", "all"} or (
        args.phase == "simulate" and not existing_plan_path.is_file()
    )
    if should_create_plan:
        if existing_plan_path.is_file() and not args.force_plan:
            existing_plan = simulation.load_plan(root)
            existing_counts = _plan_counts(existing_plan)
        else:
            existing_counts = {
                "neutral": DEFAULT_NEUTRAL_REPLICATES,
                "selected": DEFAULT_SELECTED_REPLICATES,
            }
        plan_path = simulation.write_plan(
            root,
            n_neutral=(
                args.neutral_replicates
                if args.neutral_replicates is not None
                else existing_counts["neutral"]
            ),
            n_selected=(
                args.selected_replicates
                if args.selected_replicates is not None
                else existing_counts["selected"]
            ),
            force=bool(args.force_plan),
        )
    else:
        loaded_plan = simulation.load_plan(root)
        _validate_requested_counts_against_plan(
            loaded_plan,
            n_neutral=args.neutral_replicates,
            n_selected=args.selected_replicates,
        )
        plan_path = existing_plan_path
    if args.phase == "plan":
        print(Path(plan_path).resolve())
        return 0

    if args.phase in {"simulate", "all"}:
        simulation_frame = simulation.simulate_study(
            root, slim_bin, workers=args.workers
        )
        _print_status(
            "simulation",
            simulation_frame,
        )
        _require_complete_phase(simulation_frame, label="simulation")
    if args.phase in {"decode", "all"}:
        decode_frame = simulation.decode_study(
            root,
            decoder_bin,
            workers=args.workers,
            threads=args.threads,
        )
        _print_status(
            "decode",
            decode_frame,
        )
        _require_complete_phase(decode_frame, label="decode")
    if args.phase == "status":
        for label, frame in status_tables(root).items():
            _print_status(label, frame)
        return 0

    analysis_options = {
        "output_dir": args.analysis_dir,
        "empirical_path": args.empirical_tsv,
        "n_neutral": args.neutral_replicates,
        "n_selected": args.selected_replicates,
        "bootstrap_draws": args.bootstrap_draws,
        "permutation_draws": args.permutation_draws,
        "base_seed": args.base_seed,
        "include_af_stratified": True,
        "af_stratum_width": args.af_stratum_width,
        "empirical_af_match_tolerance": args.empirical_af_match_tolerance,
    }
    if args.phase in {"analyze", "all"}:
        completion = analyze_study(root, **analysis_options)
        print(
            json.dumps(
                {
                    "analysis_status": completion["status"],
                    "analysis_cache_hit": completion["cache_hit"],
                },
                sort_keys=True,
            )
        )
    if args.phase in {"verify", "all"}:
        simulation_verification = simulation.verify_study(root)
        analysis_verification = verify_analysis(root, **analysis_options)
        print(
            json.dumps(
                {
                    "simulation_decode": _jsonable(
                        simulation_verification, repo_root=root
                    ),
                    "analysis": {
                        "status": analysis_verification["status"],
                        "contract_sha256": analysis_verification["contract_sha256"],
                    },
                },
                sort_keys=True,
            )
        )
    return 0


__all__ = [
    "DEFAULT_NEUTRAL_REPLICATES",
    "DEFAULT_ANALYSIS_DIR",
    "DEFAULT_SELECTED_REPLICATES",
    "SCHEMA_VERSION",
    "analyze_study",
    "build_parser",
    "main",
    "status_tables",
    "verify_analysis",
]
