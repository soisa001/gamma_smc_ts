"""Fail-closed EAS selected-production adapter and aggregate readiness audit.

EAS and Han selected simulations have deliberately separate production
executors.  This module authorizes and runs only the 60 EAS units from the
frozen EAS calibration.  Han is produced by ``focused_han_direct_v3``.  A
read-only readiness audit binds the two independently validated campaigns only
after Han's pilot and all-unit audits exist.

Fresh-clone validation never needs the ignored calibration work tree.  The
signed adapter sidecars bind tracked checksum tables containing the 840 raw
completion checksums; this module authenticates those compact mappings without
opening any raw calibration completion, trajectory, or tree file.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import threading
from collections.abc import Iterator, Mapping, Sequence
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from . import focused_eas_selected_audit_adapter as eas_audit
from . import focused_eas_selected_calibration as eas_calibration
from . import focused_selection_campaign as campaign_plan
from . import focused_selection_simulation as simulation
from .eas_sweep_models import EAS_POPULATION

SCHEMA_VERSION = "gamma-smc.focused-selected-production-adapter/v2"
INTEGRATION_MANIFEST_SCHEMA = "gamma-smc.focused-eas-selected-production/v2"
AGGREGATE_READINESS_SCHEMA = "gamma-smc.focused-selected-aggregate-readiness/v1"
QUARANTINE_PLAN_SCHEMA = "gamma-smc.focused-selected-quarantine-plan/v1"
EAS_ATTEMPT_LEDGER_SCHEMA = "gamma-smc.eas-selected-production-attempt-ledger/v1"

MODULE_SOURCE_PATH = "python/gamma_smc_aou/focused_selected_production.py"
WRAPPER_SOURCE_PATH = "scripts/run_focused_selected_production.py"
DEFAULT_CAMPAIGN_RELATIVE_PATH = "focused_selection_EAS_sim"
DEFAULT_INTEGRATION_RELATIVE_PATH = (
    "focused_selection_EAS_sim/calibration/focused_selected_production/"
    "integration_manifest.json"
)
DEFAULT_EAS_ADAPTER_AUTHORIZATION_RELATIVE_PATH = (
    "focused_selection_EAS_sim/calibration/eas_selected/"
    "eas_selected_calibration_frozen.adapter.json"
)
DEFAULT_EAS_FROZEN_RELATIVE_PATH = (
    "focused_selection_EAS_sim/calibration/eas_selected/"
    "eas_selected_calibration_frozen.json"
)
DEFAULT_HAN_V3_OUTPUT_RELATIVE_PATH = (
    "focused_selection_EAS_sim/calibration/han_direct_production_v3"
)
DEFAULT_HAN_V3_MANIFEST_RELATIVE_PATH = (
    f"{DEFAULT_HAN_V3_OUTPUT_RELATIVE_PATH}/han_direct_v3_manifest.json"
)
DEFAULT_HAN_V3_PILOT_AUDIT_RELATIVE_PATH = (
    f"{DEFAULT_HAN_V3_OUTPUT_RELATIVE_PATH}/pilot_audit.json"
)
DEFAULT_HAN_V3_FINAL_AUDIT_RELATIVE_PATH = (
    f"{DEFAULT_HAN_V3_OUTPUT_RELATIVE_PATH}/all_units_audit.json"
)
DEFAULT_HAN_V3_AUTHORIZATION_RELATIVE_PATH = (
    f"{DEFAULT_HAN_V3_OUTPUT_RELATIVE_PATH}/han_direct_v3_authorization.json"
)
DEFAULT_HAN_V3_RECOVERY_AUTHORIZATION_RELATIVE_PATH = (
    f"{DEFAULT_HAN_V3_OUTPUT_RELATIVE_PATH}/"
    "han_direct_v3_dtype_recovery_authorization.json"
)
DEFAULT_HAN_V3_RECOVERY_AUDIT_RELATIVE_PATH = (
    f"{DEFAULT_HAN_V3_OUTPUT_RELATIVE_PATH}/han_direct_v3_dtype_recovery_audit.json"
)
DEFAULT_AGGREGATE_READINESS_RELATIVE_PATH = (
    "focused_selection_EAS_sim/calibration/focused_selected_production/"
    "aggregate_readiness.json"
)

DEFAULT_EAS_MAX_DRAWS = 68
DEFAULT_DRAW_TIMEOUT_SECONDS = 30.0 * 60.0
DEFAULT_EXPECTED_STALE_SELECTED_DIRS = 24

INTEGRATION_IMPLEMENTATION_PATHS = (
    "python/gamma_smc_aou/.gitattributes",
    "scripts/.gitattributes",
    "pyproject.toml",
    "uv.lock",
    simulation.SIMULATION_SOURCE_PATH,
    "python/gamma_smc_aou/eas_sweep_models.py",
    "python/gamma_smc_aou/focused_selection_campaign.py",
    "python/gamma_smc_aou/focused_han_direct_v3.py",
    eas_calibration.MODULE_SOURCE_PATH,
    eas_audit.MODULE_SOURCE_PATH,
    MODULE_SOURCE_PATH,
    WRAPPER_SOURCE_PATH,
    DEFAULT_INTEGRATION_RELATIVE_PATH,
)

_PATCH_LOCK = threading.RLock()
_WORKER_BUNDLE: SelectedProductionBundle | None = None
_WORKER_PATCH_CONTEXT: Any = None


@dataclass(frozen=True)
class SelectedProductionBundle:
    """Fully validated inputs needed by one EAS selected worker process."""

    repo_root: Path
    campaign_dir: Path
    integration_manifest_path: Path
    slim_path: Path
    execution_units_path: Path
    manifest: Mapping[str, Any]
    eas_authorization: Mapping[str, Any]
    eas_frozen: Mapping[str, Any]
    implementation_paths: tuple[str, ...]


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_sha256(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(
            value, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode("utf-8")
    ).hexdigest()


def _inside_repo(path: str | Path, repo_root: Path, *, label: str) -> Path:
    candidate = Path(path)
    if not candidate.is_absolute():
        candidate = repo_root / candidate
    resolved = candidate.resolve()
    try:
        resolved.relative_to(repo_root)
    except ValueError as error:
        raise ValueError(f"{label} must remain inside repo_root") from error
    if any("onedrive" in part.casefold() for part in resolved.parts):
        raise ValueError(f"{label} must not be inside OneDrive")
    return resolved


def _file_record(path: Path, repo_root: Path) -> dict[str, Any]:
    resolved = _inside_repo(path, repo_root, label="bound file")
    if not resolved.is_file():
        raise ValueError(f"bound file is absent: {resolved}")
    return {
        "path": resolved.relative_to(repo_root).as_posix(),
        "sha256": sha256_file(resolved),
        "size_bytes": resolved.stat().st_size,
    }


def _validate_hashed_payload(
    payload: Mapping[str, Any],
    *,
    schema: str,
    status_key: str,
    expected_status: str,
    label: str,
) -> None:
    unhashed = dict(payload)
    digest = unhashed.pop("payload_sha256", None)
    if (
        payload.get("schema") != schema
        or payload.get(status_key) != expected_status
        or digest != _canonical_sha256(unhashed)
    ):
        raise ValueError(f"{label} is incompatible or corrupt")


def _read_json(path: Path, *, label: str) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"{label} is unreadable") from error
    if not isinstance(payload, dict):
        raise TypeError(f"{label} must be a JSON object")
    return payload


def _exact_file_record(
    record: Any, expected_path: Path, repo_root: Path, *, label: str
) -> None:
    expected = _file_record(expected_path, repo_root)
    if not isinstance(record, Mapping) or any(
        record.get(key) != value for key, value in expected.items()
    ):
        raise ValueError(f"{label} file binding differs")


def _load_eas_authorization_bundle(
    *,
    repo_root: Path,
    campaign_dir: Path,
    slim_path: Path,
    execution_units_path: Path,
    adapter_authorization_path: Path,
    frozen_path: Path,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Validate the compact EAS authorization chain without raw work files."""

    frozen = eas_calibration.load_frozen_eas_selected_calibration(
        frozen_path,
        repo_root=repo_root,
        execution_units_path=execution_units_path,
        slim_path=slim_path,
    )
    frozen_slim = frozen.get("slim", {})
    production_contract = frozen.get("production_contract", {})
    if (
        float(frozen_slim.get("slim_scaling_factor", math.nan))
        != simulation.DEFAULT_SLIM_SCALING_FACTOR
        or float(frozen_slim.get("slim_burn_in", math.nan))
        != simulation.DEFAULT_SLIM_BURN_IN
        or production_contract.get("selection_origin")
        != "de_novo_at_prespecified_Q5_origin_age"
        or production_contract.get("selection_episode")
        != "continuous_from_origin_to_present"
        or production_contract.get("terminal_conditions")
        != "lower_and_upper_inclusive_target_plus_minus_0.025"
        or int(production_contract.get("candidate_pool_diploids", -1))
        != simulation.DEFAULT_SELECTED_POOL_DIPLOIDS
        or int(production_contract.get("sample_diploids", -1)) != 100
    ):
        raise ValueError("EAS frozen production semantics differ")
    authorization = _read_json(
        adapter_authorization_path, label="EAS adapter authorization"
    )
    _validate_hashed_payload(
        authorization,
        schema=eas_audit.AUTHORIZATION_SCHEMA_VERSION,
        status_key="status",
        expected_status="authorized",
        label="EAS adapter authorization",
    )
    adapter_record = {
        "path": eas_audit.MODULE_SOURCE_PATH,
        "sha256": sha256_file(repo_root / eas_audit.MODULE_SOURCE_PATH),
    }
    if authorization.get("adapter") != adapter_record:
        raise ValueError("EAS adapter authorization source binding differs")
    _exact_file_record(
        authorization.get("frozen_authorization"),
        frozen_path,
        repo_root,
        label="EAS frozen authorization",
    )
    if authorization["frozen_authorization"].get("payload_sha256") != frozen.get(
        "payload_sha256"
    ):
        raise ValueError("EAS adapter and frozen payload bindings differ")

    if authorization.get("legacy_source_sha256") != frozen.get(
        "source_sha256"
    ) or authorization.get("plan_manifest_contract_sha256") != frozen.get(
        "plan_manifest_contract_sha256"
    ):
        raise ValueError("EAS adapter plan/source authorization differs")
    output = _inside_repo(
        campaign_dir / "calibration" / "eas_selected",
        repo_root,
        label="EAS compact calibration output",
    )
    expected_phases = set(eas_calibration.PHASE_ORDER)
    phase_records = authorization.get("phase_audits")
    if not isinstance(phase_records, Mapping) or set(phase_records) != expected_phases:
        raise ValueError("EAS adapter authorization lacks exact phase coverage")
    adapter_sha256 = adapter_record["sha256"]
    for phase in eas_calibration.PHASE_ORDER:
        sidecar_path = output / f"{phase}_adapter_audit.json"
        sidecar = _load_compact_phase_sidecar(
            phase,
            output=output,
            repo_root=repo_root,
            adapter_sha256=adapter_sha256,
            expected_source_sha256=frozen["source_sha256"],
            expected_plan_contract_sha256=str(frozen["plan_manifest_contract_sha256"]),
            slim_path=slim_path,
        )
        _exact_file_record(
            phase_records[phase], sidecar_path, repo_root, label=f"EAS {phase} audit"
        )
        if phase_records[phase].get("payload_sha256") != sidecar.get("payload_sha256"):
            raise ValueError(f"EAS {phase} audit payload binding differs")
    return authorization, frozen


def _compact_completion_mapping(
    audit: pd.DataFrame, *, phase: str, repo_root: Path
) -> dict[str, Any]:
    """Canonicalize tracked checksum rows without dereferencing raw paths."""

    columns = (
        "calibration_id",
        "completion_path",
        "completion_sha256",
        "trajectory_sha256",
        "panel_tree_sha256",
        "panel_manifest_sha256",
    )
    if set(columns).difference(audit.columns):
        raise ValueError(f"{phase} EAS checksum table columns differ")
    canonical = audit.loc[:, columns].copy()
    for column in ("panel_tree_sha256", "panel_manifest_sha256"):
        canonical[column] = canonical[column].fillna("")
    if canonical["calibration_id"].astype(str).duplicated().any():
        raise ValueError(f"{phase} EAS checksum calibration IDs are duplicated")
    mapping = canonical.sort_values("calibration_id").to_dict(orient="records")
    expected_prefix = (
        f"{DEFAULT_CAMPAIGN_RELATIVE_PATH}/work/eas_selected_calibration/{phase}/"
    )
    hexadecimal = set("0123456789abcdef")
    for row in mapping:
        completion_path = str(row["completion_path"])
        path = Path(completion_path)
        if (
            path.is_absolute()
            or ".." in path.parts
            or not completion_path.startswith(expected_prefix)
            or not completion_path.endswith("/completion.json")
        ):
            raise ValueError(f"{phase} EAS completion provenance path differs")
        # Resolve only lexically to reject escapes.  Do not stat/open the raw file.
        _inside_repo(repo_root / completion_path, repo_root, label="EAS raw provenance")
        for column in (
            "completion_sha256",
            "trajectory_sha256",
            "panel_tree_sha256",
            "panel_manifest_sha256",
        ):
            value = str(row[column])
            if value and (len(value) != 64 or set(value).difference(hexadecimal)):
                raise ValueError(f"{phase} EAS {column} is not a SHA-256 digest")
        if not str(row["completion_sha256"]) or not str(row["trajectory_sha256"]):
            raise ValueError(f"{phase} EAS required completion digests are absent")
    return {
        "rows": len(mapping),
        "canonical_sha256": _canonical_sha256(mapping),
        "entries": mapping,
    }


def _load_compact_phase_sidecar(
    phase: str,
    *,
    repo_root: Path,
    output: Path,
    adapter_sha256: str,
    expected_source_sha256: Mapping[str, Any],
    expected_plan_contract_sha256: str,
    slim_path: Path,
) -> dict[str, Any]:
    """Authenticate one signed sidecar and its tracked checksum mapping."""

    sidecar_path = output / f"{phase}_adapter_audit.json"
    payload = _read_json(sidecar_path, label=f"{phase} EAS adapter sidecar")
    _validate_hashed_payload(
        payload,
        schema=eas_audit.SCHEMA_VERSION,
        status_key="status",
        expected_status="validated",
        label=f"{phase} EAS adapter sidecar",
    )
    if payload.get("phase") != phase or payload.get("adapter") != {
        "path": eas_audit.MODULE_SOURCE_PATH,
        "sha256": adapter_sha256,
    }:
        raise ValueError(f"{phase} EAS adapter/phase binding differs")
    if payload.get("legacy_source_sha256") != expected_source_sha256:
        raise ValueError(f"{phase} EAS source binding differs")

    plan_manifest_path = output / "eas_selected_calibration_plan.json"
    plan_manifest_record = dict(payload.get("plan_manifest", {}))
    if (
        plan_manifest_record.pop("contract_sha256", None)
        != expected_plan_contract_sha256
    ):
        raise ValueError(f"{phase} EAS plan contract binding differs")
    _exact_file_record(
        plan_manifest_record,
        plan_manifest_path,
        repo_root,
        label=f"{phase} EAS plan manifest",
    )
    _exact_file_record(
        payload.get("plan"),
        output / f"{phase}_plan.tsv",
        repo_root,
        label=f"{phase} EAS plan",
    )
    expected_outputs = {
        "rebuilt_ledger": output / f"{phase}_ledger.tsv",
        "summary": output / f"{phase}_summary.tsv",
        "completion_checksums": output / f"{phase}_completion_checksums.tsv",
    }
    if set(payload.get("outputs", {})) != set(expected_outputs):
        raise ValueError(f"{phase} EAS compact output set differs")
    for key, path in expected_outputs.items():
        _exact_file_record(
            payload["outputs"].get(key), path, repo_root, label=f"{phase} EAS {key}"
        )
    _exact_file_record(
        payload.get("normalization", {}).get("input_ledger"),
        expected_outputs["rebuilt_ledger"],
        repo_root,
        label=f"{phase} EAS normalized ledger",
    )
    _exact_file_record(
        payload.get("slim"), slim_path, repo_root, label=f"{phase} EAS SLiM"
    )
    _exact_file_record(
        payload.get("eas_resource"),
        repo_root / eas_calibration.EAS_RESOURCE_PATH,
        repo_root,
        label=f"{phase} EAS resource",
    )

    expected_rows = 6 * int(eas_calibration.PHASE_DRAWS[phase])
    counts = payload.get("counts", {})
    if counts != {
        "planned": expected_rows,
        "validated_completions": expected_rows,
        "completion_checksum_rows": expected_rows,
        "cells": 6,
        "passing_cells": 6,
    }:
        raise ValueError(f"{phase} EAS sidecar counts differ")
    checksum_table = pd.read_csv(
        expected_outputs["completion_checksums"],
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )
    compact_mapping = _compact_completion_mapping(
        checksum_table, phase=phase, repo_root=repo_root
    )
    if compact_mapping != payload.get("completion_mapping"):
        raise ValueError(f"{phase} EAS compact completion mapping differs")
    return payload


def _validate_execution_units(units: pd.DataFrame) -> dict[str, Any]:
    required = set(campaign_plan.EXECUTION_COLUMNS)
    if required.difference(units.columns) or len(units) != 720:
        raise ValueError("focused execution units must contain the exact 720-unit plan")
    if units["unit_id"].duplicated().any() or units["seed"].duplicated().any():
        raise ValueError("focused execution-unit identifiers or seeds are duplicated")
    classes = units.groupby("simulation_class").size().to_dict()
    if classes != {"neutral": 600, "selected": 120}:
        raise ValueError("focused execution-unit class counts differ")
    selected = units[units["simulation_class"].astype(str).eq("selected")].copy()
    neutral = units[units["simulation_class"].astype(str).eq("neutral")].copy()
    demographies = {
        "ancient_eurasia_han_introgression",
        "eas_phlash_median",
    }
    if set(selected["demography_id"].astype(str)) != demographies:
        raise ValueError("selected execution-unit demographies differ")
    selected_cells = selected.groupby(
        ["demography_id", "selection_coefficient", "target_allele_frequency"]
    ).size()
    neutral_cells = neutral.groupby(["demography_id", "target_allele_frequency"]).size()
    if (
        len(selected_cells) != 12
        or set(selected_cells.astype(int)) != {10}
        or len(neutral_cells) != 6
        or set(neutral_cells.astype(int)) != {100}
    ):
        raise ValueError("focused execution-unit cell replication differs")
    if set(selected["selection_coefficient"].astype(float)) != {0.005, 0.01}:
        raise ValueError("selected coefficient grid differs")
    if set(selected["target_allele_frequency"].astype(float)) != {0.1, 0.2, 0.3}:
        raise ValueError("selected target-frequency grid differs")
    if not np.allclose(
        selected["population_af_upper"].astype(float)
        - selected["target_allele_frequency"].astype(float),
        0.025,
        rtol=0.0,
        atol=1e-12,
    ) or not np.allclose(
        selected["target_allele_frequency"].astype(float)
        - selected["population_af_lower"].astype(float),
        0.025,
        rtol=0.0,
        atol=1e-12,
    ):
        raise ValueError("selected candidate-pool AF bands differ")
    expected_alt = np.rint(
        2
        * selected["sample_diploids"].astype(int)
        * selected["target_allele_frequency"].astype(float)
    ).astype(int)
    if set(selected["sample_diploids"].astype(int)) != {100} or not np.array_equal(
        selected["exact_sample_alt_count"].astype(int).to_numpy(),
        expected_alt.to_numpy(),
    ):
        raise ValueError("selected exact-k panel contract differs")
    return {
        "total_units": 720,
        "selected_units": 120,
        "neutral_units": 600,
        "selected_cells": 12,
        "selected_replicates_per_cell": 10,
        "neutral_cells": 6,
        "neutral_replicates_per_cell": 100,
    }


def _implementation_records(repo_root: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    for relative in INTEGRATION_IMPLEMENTATION_PATHS[:-1]:
        path = _inside_repo(repo_root / relative, repo_root, label="implementation")
        if not path.is_file():
            raise ValueError(f"selected implementation source is absent: {relative}")
        records[relative] = sha256_file(path)
    return records


def _binding_summary(
    *,
    eas_authorization: Mapping[str, Any],
    eas_frozen: Mapping[str, Any],
) -> dict[str, Any]:
    return {
        "adapter_payload_sha256": eas_authorization["payload_sha256"],
        "frozen_payload_sha256": eas_frozen["payload_sha256"],
        "cell_authorization_count": len(eas_frozen["cell_authorizations"]),
        "production_units": 60,
        "terminal_population_af_conditioning_in_slim": True,
        "continuous_selection_to_present": True,
        "raw_calibration_completion_files_required": False,
    }


def build_integration_manifest(
    *,
    repo_root: str | Path,
    campaign_dir: str | Path,
    slim_path: str | Path,
    eas_adapter_authorization_path: str | Path,
    eas_frozen_path: str | Path,
    integration_manifest_path: str | Path = DEFAULT_INTEGRATION_RELATIVE_PATH,
) -> dict[str, Any]:
    """Build, but do not publish, the independent EAS authorization."""

    root = Path(repo_root).resolve()
    if any("onedrive" in part.casefold() for part in root.parts):
        raise ValueError("selected production repo_root must not be inside OneDrive")
    campaign = _inside_repo(campaign_dir, root, label="campaign directory")
    units_path = _inside_repo(
        campaign / "execution_units.tsv", root, label="execution units"
    )
    slim = _inside_repo(slim_path, root, label="SLiM binary")
    eas_adapter_path = _inside_repo(
        eas_adapter_authorization_path, root, label="EAS adapter authorization"
    )
    eas_frozen = _inside_repo(eas_frozen_path, root, label="EAS frozen authorization")
    integration = _inside_repo(
        integration_manifest_path, root, label="integration manifest"
    )
    if integration.relative_to(root).as_posix() != DEFAULT_INTEGRATION_RELATIVE_PATH:
        raise ValueError("selected integration manifest path is frozen")
    units = pd.read_csv(units_path, sep="\t")
    counts = _validate_execution_units(units)
    eas_authorization, eas_payload = _load_eas_authorization_bundle(
        repo_root=root,
        campaign_dir=campaign,
        slim_path=slim,
        execution_units_path=units_path,
        adapter_authorization_path=eas_adapter_path,
        frozen_path=eas_frozen,
    )
    payload: dict[str, Any] = {
        "schema": INTEGRATION_MANIFEST_SCHEMA,
        "status": "eas_ready",
        "adapter_schema": SCHEMA_VERSION,
        "execution_units": {
            **_file_record(units_path, root),
            **counts,
        },
        "runtime": {
            "slim": _file_record(slim, root),
            "slim_scaling_factor": simulation.DEFAULT_SLIM_SCALING_FACTOR,
            "slim_burn_in": simulation.DEFAULT_SLIM_BURN_IN,
            "candidate_pool_diploids": simulation.DEFAULT_SELECTED_POOL_DIPLOIDS,
            "sample_diploids": 100,
            "eas_max_draws": DEFAULT_EAS_MAX_DRAWS,
            "per_draw_timeout_seconds": DEFAULT_DRAW_TIMEOUT_SECONDS,
            "cumulative_timeout_seconds": (
                DEFAULT_EAS_MAX_DRAWS * DEFAULT_DRAW_TIMEOUT_SECONDS
            ),
        },
        "authorization_artifacts": {
            "eas_adapter": _file_record(eas_adapter_path, root),
            "eas_frozen": _file_record(eas_frozen, root),
        },
        "authorization_bindings": _binding_summary(
            eas_authorization=eas_authorization,
            eas_frozen=eas_payload,
        ),
        "implementation": {
            "mode": "eas_selected_worker_process_local_patch",
            "neutral_executor_modified": False,
            "campaign_plan_modified": False,
            "selected_objects_patch_only": True,
            "han_executor_or_artifacts_modified": False,
            "contract_paths": list(INTEGRATION_IMPLEMENTATION_PATHS),
            "source_sha256": _implementation_records(root),
        },
        "scientific_contract": {
            "eas": (
                "Q5 de-novo origin; continuous selection to present; terminal "
                "population AF bounds retained in SLiM; external 500-diploid AF "
                "gate and exact-k N100 panel retained; exactly 60 production units"
            ),
            "han": "separate direct-production-v3 executor and authorization",
            "gamma_smc_statistics_used_for_calibration": False,
        },
        "quarantine": {
            "planning_only": True,
            "moves_or_deletions_authorized": False,
            "expected_preintegration_stale_eas_selected_directories": (
                DEFAULT_EXPECTED_STALE_SELECTED_DIRS
            ),
            "han_directories_in_scope": False,
            "neutral_directories_in_scope": False,
        },
    }
    payload["payload_sha256"] = _canonical_sha256(payload)
    return payload


def write_integration_manifest(path: str | Path, payload: Mapping[str, Any]) -> Path:
    """Publish an immutable integration manifest idempotently."""

    destination = Path(path).resolve()
    if any("onedrive" in part.casefold() for part in destination.parts):
        raise ValueError("selected integration manifest must not be inside OneDrive")
    _validate_hashed_payload(
        payload,
        schema=INTEGRATION_MANIFEST_SCHEMA,
        status_key="status",
        expected_status="eas_ready",
        label="selected integration manifest",
    )
    encoded = (
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")
    if destination.exists():
        if destination.read_bytes() != encoded:
            raise FileExistsError(
                f"integration manifest already exists with different bytes: {destination}"
            )
        return destination
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.{os.getpid()}.tmp")
    temporary.write_bytes(encoded)
    os.replace(temporary, destination)
    return destination


def _records_match(
    left: Mapping[str, Any], right: Mapping[str, Any], columns: Sequence[str]
) -> bool:
    for column in columns:
        a = left[column]
        b = right[column]
        if isinstance(a, (float, np.floating)) or isinstance(b, (float, np.floating)):
            if not math.isclose(float(a), float(b), rel_tol=0.0, abs_tol=1e-12):
                return False
        elif str(a) != str(b):
            return False
    return True


def load_selected_production_bundle(
    path: str | Path,
    *,
    repo_root: str | Path,
    slim_path: str | Path,
) -> SelectedProductionBundle:
    """Strictly reload the independent EAS production authorization."""

    root = Path(repo_root).resolve()
    manifest_path = _inside_repo(path, root, label="integration manifest")
    if manifest_path.relative_to(root).as_posix() != DEFAULT_INTEGRATION_RELATIVE_PATH:
        raise ValueError("selected integration manifest path differs")
    payload = _read_json(manifest_path, label="selected integration manifest")
    _validate_hashed_payload(
        payload,
        schema=INTEGRATION_MANIFEST_SCHEMA,
        status_key="status",
        expected_status="eas_ready",
        label="EAS selected integration manifest",
    )
    slim = _inside_repo(slim_path, root, label="SLiM binary")
    runtime = payload.get("runtime", {})
    _exact_file_record(runtime.get("slim"), slim, root, label="SLiM")
    if (
        float(runtime.get("slim_scaling_factor", math.nan))
        != simulation.DEFAULT_SLIM_SCALING_FACTOR
        or float(runtime.get("slim_burn_in", math.nan))
        != simulation.DEFAULT_SLIM_BURN_IN
        or int(runtime.get("candidate_pool_diploids", -1))
        != simulation.DEFAULT_SELECTED_POOL_DIPLOIDS
        or int(runtime.get("sample_diploids", -1)) != 100
        or int(runtime.get("eas_max_draws", -1)) != DEFAULT_EAS_MAX_DRAWS
        or float(runtime.get("per_draw_timeout_seconds", math.nan))
        != DEFAULT_DRAW_TIMEOUT_SECONDS
        or float(runtime.get("cumulative_timeout_seconds", math.nan))
        != DEFAULT_EAS_MAX_DRAWS * DEFAULT_DRAW_TIMEOUT_SECONDS
    ):
        raise ValueError("selected integration runtime contract differs")
    units_path = _inside_repo(
        root / str(payload.get("execution_units", {}).get("path", "")),
        root,
        label="execution units",
    )
    unit_record = dict(payload["execution_units"])
    for key in (
        "total_units",
        "selected_units",
        "neutral_units",
        "selected_cells",
        "selected_replicates_per_cell",
        "neutral_cells",
        "neutral_replicates_per_cell",
    ):
        unit_record.pop(key, None)
    _exact_file_record(unit_record, units_path, root, label="execution units")
    counts = _validate_execution_units(pd.read_csv(units_path, sep="\t"))
    if any(
        payload["execution_units"].get(key) != value for key, value in counts.items()
    ):
        raise ValueError("selected integration execution-unit counts differ")
    campaign = units_path.parent

    artifacts = payload.get("authorization_artifacts")
    if not isinstance(artifacts, Mapping) or set(artifacts) != {
        "eas_adapter",
        "eas_frozen",
    }:
        raise ValueError("selected integration authorization artifact set differs")
    artifact_paths: dict[str, Path] = {}
    for key, record in artifacts.items():
        if not isinstance(record, Mapping):
            raise TypeError(f"selected integration artifact record is invalid: {key}")
        artifact_path = _inside_repo(
            root / str(record.get("path", "")), root, label=f"{key} artifact"
        )
        _exact_file_record(record, artifact_path, root, label=key)
        artifact_paths[key] = artifact_path

    eas_authorization, eas_frozen = _load_eas_authorization_bundle(
        repo_root=root,
        campaign_dir=campaign,
        slim_path=slim,
        execution_units_path=units_path,
        adapter_authorization_path=artifact_paths["eas_adapter"],
        frozen_path=artifact_paths["eas_frozen"],
    )
    if payload.get("authorization_bindings") != _binding_summary(
        eas_authorization=eas_authorization,
        eas_frozen=eas_frozen,
    ):
        raise ValueError("selected integration authorization bindings differ")

    implementation = payload.get("implementation")
    if (
        not isinstance(implementation, Mapping)
        or implementation.get("mode") != "eas_selected_worker_process_local_patch"
        or implementation.get("neutral_executor_modified") is not False
        or implementation.get("campaign_plan_modified") is not False
        or implementation.get("selected_objects_patch_only") is not True
        or implementation.get("han_executor_or_artifacts_modified") is not False
        or tuple(implementation.get("contract_paths", ()))
        != INTEGRATION_IMPLEMENTATION_PATHS
        or implementation.get("source_sha256") != _implementation_records(root)
    ):
        raise ValueError("selected integration implementation binding differs")
    if (
        INTEGRATION_IMPLEMENTATION_PATHS[-1]
        != manifest_path.relative_to(root).as_posix()
    ):
        raise ValueError("integration manifest is absent from worker source contract")
    return SelectedProductionBundle(
        repo_root=root,
        campaign_dir=campaign,
        integration_manifest_path=manifest_path,
        slim_path=slim,
        execution_units_path=units_path,
        manifest=payload,
        eas_authorization=eas_authorization,
        eas_frozen=eas_frozen,
        implementation_paths=INTEGRATION_IMPLEMENTATION_PATHS,
    )


def _eas_cell_key(coefficient: float, target: float) -> str:
    return f"s={coefficient:.6f}|af={target:.6f}"


def _authorized_eas_selected_objects(
    unit: Mapping[str, Any],
    *,
    bundle: SelectedProductionBundle,
    slim_scaling_factor: float,
    selected_pool_diploids: int,
) -> tuple[Any, Any, dict[str, int], dict[str, Any]]:
    record = simulation._validate_unit(unit)
    coefficient = float(record["selection_coefficient"])
    target = float(record["target_allele_frequency"])
    lower = float(record["population_af_lower"])
    upper = float(record["population_af_upper"])
    if (
        not math.isclose(
            float(slim_scaling_factor),
            simulation.DEFAULT_SLIM_SCALING_FACTOR,
            rel_tol=0.0,
            abs_tol=1e-12,
        )
        or int(selected_pool_diploids) != simulation.DEFAULT_SELECTED_POOL_DIPLOIDS
        or not math.isclose(target - lower, 0.025, rel_tol=0.0, abs_tol=1e-12)
        or not math.isclose(upper - target, 0.025, rel_tol=0.0, abs_tol=1e-12)
    ):
        raise ValueError("EAS selected production Q, pool, or AF band differs")
    model, sweep, _, origin, event = eas_calibration._load_eas_cell(
        bundle.repo_root,
        selection_coefficient=coefficient,
        target_frequency=target,
        af_half_width=0.025,
        slim_scaling_factor=slim_scaling_factor,
    )
    key = _eas_cell_key(coefficient, target)
    authorization = bundle.eas_frozen.get("cell_authorizations", {}).get(key)
    if not isinstance(authorization, Mapping):
        raise TypeError(f"EAS selected authorization lacks cell {key}")
    if (
        str(authorization.get("origin_contract_sha256"))
        != eas_calibration._canonical_sha256(origin)
        or str(authorization.get("event_contract_sha256"))
        != str(event["event_contract_sha256"])
        or not math.isclose(
            float(authorization.get("origin_age_generations", math.nan)),
            float(origin["age_generations"]),
            rel_tol=0.0,
            abs_tol=1e-12,
        )
    ):
        raise ValueError(f"EAS selected scientific authorization differs: {key}")
    serialized_events = eas_calibration.serialize_extended_events(sweep.extended_events)
    terminal = [
        item
        for item in serialized_events
        if item.get("event_type") == "ConditionOnAlleleFrequency"
        and math.isclose(
            float(item.get("start_time", {}).get("generations_ago", math.nan)),
            0.0,
            abs_tol=1e-12,
        )
        and math.isclose(
            float(item.get("end_time", {}).get("generations_ago", math.nan)),
            0.0,
            abs_tol=1e-12,
        )
    ]
    if len(terminal) != 2 or {str(item.get("operator")) for item in terminal} != {
        ">=",
        "<=",
    }:
        raise ValueError("EAS selected terminal population-AF conditions were lost")
    return (
        model,
        sweep,
        {EAS_POPULATION: int(selected_pool_diploids)},
        {
            "sweep": sweep.provenance_record(),
            "origin_age": origin,
            "selected_production_authorization": {
                "schema": SCHEMA_VERSION,
                "demography_id": "eas_phlash_median",
                "cell_key": key,
                "integration_manifest_path": bundle.integration_manifest_path.relative_to(
                    bundle.repo_root
                ).as_posix(),
                "integration_payload_sha256": bundle.manifest["payload_sha256"],
                "eas_adapter_payload_sha256": bundle.eas_authorization[
                    "payload_sha256"
                ],
                "eas_frozen_payload_sha256": bundle.eas_frozen["payload_sha256"],
                "terminal_population_af_conditioning_in_slim": True,
                "external_candidate_pool_and_exact_panel_retained": True,
                "process_group_supervision": {
                    "implementation": "shared_hardened_han_v3_supervisor",
                    "implementation_path": (
                        "python/gamma_smc_aou/focused_han_direct_v3.py"
                    ),
                    "draw_timeout_seconds": DEFAULT_DRAW_TIMEOUT_SECONDS,
                    "term_then_kill_full_process_group": True,
                    "direct_child_and_adopted_descendants_reaped": True,
                    "attempt_lifecycle_journaled": True,
                    "accepted_precompletion_crash_reconciled": True,
                },
            },
        },
    )


@contextmanager
def patched_selected_executor(bundle: SelectedProductionBundle) -> Iterator[None]:
    """Patch selected construction and source bindings in this process only."""

    with _PATCH_LOCK:
        original_objects = simulation._selected_objects
        original_paths = simulation.IMPLEMENTATION_PATHS

        def selected_objects(
            unit: Mapping[str, Any],
            repo_root: str | Path,
            slim_scaling_factor: float,
            selected_pool_diploids: int,
            *,
            han_selection_calibration: Mapping[str, Any] | None = None,
            han_selection_end_generations_ago: float | None = None,
        ) -> tuple[Any, Any, dict[str, int], dict[str, Any]]:
            record = simulation._validate_unit(unit)
            if Path(repo_root).resolve() != bundle.repo_root:
                raise ValueError("selected worker repo_root differs from integration")
            if record["simulation_class"] != "selected":
                raise ValueError("selected production adapter refuses neutral units")
            if record["demography_id"] == "eas_phlash_median":
                if (
                    han_selection_end_generations_ago is not None
                    or han_selection_calibration is not None
                ):
                    raise ValueError("EAS selected adapter received a Han endpoint")
                return _authorized_eas_selected_objects(
                    record,
                    bundle=bundle,
                    slim_scaling_factor=slim_scaling_factor,
                    selected_pool_diploids=selected_pool_diploids,
                )
            if record["demography_id"] == "ancient_eurasia_han_introgression":
                raise ValueError(
                    "Han selected production is owned by focused_han_direct_v3"
                )
            raise ValueError(
                f"selected production adapter rejects demography {record['demography_id']!r}"
            )

        simulation._selected_objects = selected_objects
        simulation.IMPLEMENTATION_PATHS = bundle.implementation_paths
        try:
            yield
        finally:
            simulation._selected_objects = original_objects
            simulation.IMPLEMENTATION_PATHS = original_paths


def _validate_requested_units(
    requested: pd.DataFrame, canonical: pd.DataFrame
) -> pd.DataFrame:
    expected = _eas_production_units(canonical)
    if len(requested) != 60:
        raise ValueError("EAS selected production requires exactly 60 units")
    if set(requested["simulation_class"].astype(str)) != {"selected"} or set(
        requested["demography_id"].astype(str)
    ) != {"eas_phlash_median"}:
        raise ValueError("EAS selected production refuses neutral or Han units")
    if requested["unit_id"].astype(str).duplicated().any():
        raise ValueError("EAS selected production contains duplicate unit IDs")
    if set(requested["unit_id"].astype(str)) != set(expected["unit_id"].astype(str)):
        raise ValueError("EAS selected production unit set differs from the exact 60")
    canonical_by_id = {
        str(row["unit_id"]): row for row in canonical.to_dict(orient="records")
    }
    for row in requested.to_dict(orient="records"):
        expected = canonical_by_id.get(str(row.get("unit_id", "")))
        if expected is None or not _records_match(
            row, expected, campaign_plan.EXECUTION_COLUMNS
        ):
            raise ValueError(
                f"selected requested unit differs from frozen plan: {row.get('unit_id')}"
            )
    return requested.reset_index(drop=True)


def _eas_production_units(canonical: pd.DataFrame) -> pd.DataFrame:
    """Return the immutable 2 x 3 x 10 EAS selected-production grid."""

    _validate_execution_units(canonical)
    eas = canonical[
        canonical["simulation_class"].astype(str).eq("selected")
        & canonical["demography_id"].astype(str).eq("eas_phlash_median")
    ].copy()
    cells = eas.groupby(["selection_coefficient", "target_allele_frequency"]).size()
    if (
        len(eas) != 60
        or len(cells) != 6
        or set(cells.astype(int)) != {10}
        or eas["unit_id"].astype(str).duplicated().any()
        or eas["seed"].astype(int).duplicated().any()
    ):
        raise ValueError("EAS selected production grid differs from 60 immutable units")
    return eas.sort_values("unit_id").reset_index(drop=True)


def _worker_configuration(
    bundle: SelectedProductionBundle,
) -> dict[str, str]:
    return {
        "repo_root": str(bundle.repo_root),
        "integration_manifest_path": str(bundle.integration_manifest_path),
        "slim_path": str(bundle.slim_path),
    }


def _eas_supervisor_schedule(unit: Mapping[str, Any]) -> pd.DataFrame:
    """Freeze the exact 68-attempt EAS simulation and panel seed schedule."""

    record = simulation._validate_unit(unit)
    if (
        record["simulation_class"] != "selected"
        or record["demography_id"] != "eas_phlash_median"
    ):
        raise ValueError("EAS supervisor schedule refuses non-EAS units")
    base_seed = int(record["seed"])
    return pd.DataFrame(
        [
            {
                "attempt_zero_based": attempt,
                "simulation_seed": simulation._stable_seed(
                    base_seed, f"selected:{attempt}"
                ),
                "panel_seed": simulation._stable_seed(
                    base_seed, f"selected-panel:{attempt}"
                ),
            }
            for attempt in range(DEFAULT_EAS_MAX_DRAWS)
        ]
    )


@contextmanager
def _patched_eas_process_supervisor(
    bundle: SelectedProductionBundle, unit: Mapping[str, Any]
) -> Iterator[Any]:
    """Apply the hardened process-group and crash journal to one EAS unit."""

    from . import focused_han_direct_v3 as hardened_supervisor

    record = simulation._validate_unit(unit)
    unit_dir = bundle.campaign_dir / "work" / str(record["unit_id"])
    unit_dir.mkdir(parents=True, exist_ok=True)
    schedule = _eas_supervisor_schedule(record)
    by_attempt = schedule.set_index("attempt_zero_based")
    with _PATCH_LOCK:
        original_atomic_frame = simulation._atomic_frame
        original_wall_clock_timeout = simulation._wall_clock_timeout
        original_get_engine = simulation.stdpopsim.get_engine
        original_cumulative = hardened_supervisor.CUMULATIVE_TIMEOUT_SECONDS
        hardened_supervisor.CUMULATIVE_TIMEOUT_SECONDS = (
            DEFAULT_EAS_MAX_DRAWS * DEFAULT_DRAW_TIMEOUT_SECONDS
        )
        supervisor = hardened_supervisor._SupervisedSlimEngine(
            original_get_engine("slim"), unit_dir=unit_dir, schedule=schedule
        )

        def atomic_frame(path: Path, frame: pd.DataFrame) -> None:
            if Path(path).name != "attempts.tsv" or "attempt_zero_based" not in frame:
                original_atomic_frame(path, frame)
                return
            enriched = frame.copy()
            for raw_attempt in enriched["attempt_zero_based"]:
                attempt = int(raw_attempt)
                if attempt not in by_attempt.index or attempt not in supervisor.records:
                    raise ValueError(
                        "EAS attempt lacks its immutable seed or supervisor record"
                    )
                if int(
                    enriched.loc[
                        enriched["attempt_zero_based"].astype(int).eq(attempt), "seed"
                    ].iloc[0]
                ) != int(by_attempt.loc[attempt, "simulation_seed"]):
                    raise ValueError("EAS attempt simulation seed differs")
            columns = (
                ("eas_supervisor_schema", "schema"),
                ("eas_supervisor_execution_id", "execution_id"),
                ("eas_supervisor_status", "supervisor_status"),
                ("eas_supervisor_pid", "supervisor_pid"),
                (
                    "eas_supervisor_process_identity",
                    "supervisor_process_identity",
                ),
                ("eas_supervisor_started_unix_ns", "supervisor_started_unix_ns"),
                ("eas_child_pid", "child_pid"),
                ("eas_child_pgid", "child_pgid"),
                ("eas_child_sid", "child_sid"),
                ("eas_child_process_identity", "child_process_identity"),
                ("eas_child_identity_observed", "child_identity_observed"),
                ("eas_engine_wrapper_pid", "engine_wrapper_pid"),
                (
                    "eas_engine_wrapper_process_identity",
                    "engine_wrapper_process_identity",
                ),
                ("eas_parent_death_watchdog_armed", "parent_death_watchdog_armed"),
                ("eas_supervisor_timeout_seconds", "timeout_seconds"),
                ("eas_supervisor_term_grace_seconds", "term_grace_seconds"),
                ("eas_sigterm_sent", "sigterm_sent"),
                ("eas_sigkill_sent", "sigkill_sent"),
                ("eas_direct_child_reaped", "direct_child_reaped"),
                ("eas_adopted_children_reaped", "adopted_children_reaped"),
                ("eas_process_group_gone", "process_group_gone"),
                ("eas_attempt_temp_directory", "temp_directory"),
                ("eas_attempt_temp_tree_created", "temp_tree_created"),
                ("eas_attempt_temp_tree_loaded", "temp_tree_loaded"),
                ("eas_attempt_temp_tree_removed", "temp_tree_removed"),
                ("eas_supervisor_elapsed_seconds", "supervisor_elapsed_seconds"),
            )
            for output_column, source_column in columns:
                enriched[output_column] = [
                    supervisor.records[int(attempt)].get(source_column)
                    for attempt in enriched["attempt_zero_based"]
                ]
            enriched["eas_planned_simulation_seed"] = [
                int(by_attempt.loc[int(attempt), "simulation_seed"])
                for attempt in enriched["attempt_zero_based"]
            ]
            enriched["eas_planned_panel_seed"] = [
                int(by_attempt.loc[int(attempt), "panel_seed"])
                for attempt in enriched["attempt_zero_based"]
            ]
            enriched["eas_panel_seed_used"] = (
                enriched.get("accepted", pd.Series(False, index=enriched.index))
                .fillna(False)
                .astype(bool)
            )
            if "rejection_reason" not in enriched:
                enriched["rejection_reason"] = ""
            enriched["rejection_reason"] = enriched["rejection_reason"].astype("object")
            timed_out = (
                enriched.get("status", pd.Series("", index=enriched.index))
                .astype(str)
                .eq("timed_out")
            )
            blank_reason = enriched["rejection_reason"].fillna("").astype(str).eq("")
            enriched.loc[timed_out & blank_reason, "rejection_reason"] = "draw_timeout"
            enriched["eas_attempt_ledger_schema"] = EAS_ATTEMPT_LEDGER_SCHEMA
            enriched["eas_supervisor_source_sha256"] = sha256_file(
                bundle.repo_root / hardened_supervisor.MODULE_SOURCE_PATH
            )
            original_atomic_frame(path, enriched)
            supervisor.mark_base_attempt_ledger_committed(enriched)

        def get_engine(name: str, *args: Any, **kwargs: Any) -> Any:
            if name == "slim":
                if args or kwargs:
                    raise ValueError("EAS SLiM engine lookup arguments differ")
                return supervisor
            return original_get_engine(name, *args, **kwargs)

        simulation._atomic_frame = atomic_frame
        simulation._wall_clock_timeout = hardened_supervisor._supervised_draw_timeout
        simulation.stdpopsim.get_engine = get_engine
        completed = False
        try:
            yield supervisor
            completed = True
        finally:
            try:
                if completed:
                    supervisor.finalize_unit_completion()
            finally:
                simulation._atomic_frame = original_atomic_frame
                simulation._wall_clock_timeout = original_wall_clock_timeout
                simulation.stdpopsim.get_engine = original_get_engine
                hardened_supervisor.CUMULATIVE_TIMEOUT_SECONDS = original_cumulative


def _initialize_selected_worker(configuration: Mapping[str, str]) -> None:
    global _WORKER_BUNDLE, _WORKER_PATCH_CONTEXT
    _WORKER_BUNDLE = load_selected_production_bundle(
        configuration["integration_manifest_path"],
        repo_root=configuration["repo_root"],
        slim_path=configuration["slim_path"],
    )
    _WORKER_PATCH_CONTEXT = patched_selected_executor(_WORKER_BUNDLE)
    _WORKER_PATCH_CONTEXT.__enter__()


def _simulate_selected_task(task: Mapping[str, Any]) -> dict[str, Any]:
    try:
        if _WORKER_BUNDLE is None:
            raise RuntimeError("selected production worker was not initialized")
        unit = task["unit"]
        unit_dir = _WORKER_BUNDLE.campaign_dir / "work" / str(unit["unit_id"])
        with (
            simulation.exclusive_unit_lock(
                unit_dir / ".simulation.lock",
                unit_id=str(unit["unit_id"]),
                phase="eas_selected_production_supervised",
            ),
            _patched_eas_process_supervisor(_WORKER_BUNDLE, unit),
        ):
            artifacts = simulation._simulate_unit_unlocked(
                unit,
                _WORKER_BUNDLE.repo_root,
                _WORKER_BUNDLE.campaign_dir,
                slim_path=_WORKER_BUNDLE.slim_path,
                selected_max_draws=int(task["selected_max_draws"]),
                slim_scaling_factor=simulation.DEFAULT_SLIM_SCALING_FACTOR,
                slim_burn_in=simulation.DEFAULT_SLIM_BURN_IN,
                selected_draw_timeout_seconds=float(task["draw_timeout_seconds"]),
                selected_cumulative_timeout_seconds=float(
                    task["cumulative_timeout_seconds"]
                ),
            )
        return {
            "unit_id": unit["unit_id"],
            "status": "cached" if artifacts.cache_hit else "complete",
            "completion_path": str(artifacts.completion_path),
            "error_type": "",
            "error": "",
        }
    except Exception as error:  # noqa: BLE001 - worker must report all failures
        return {
            "unit_id": task["unit"]["unit_id"],
            "status": "failed",
            "completion_path": "",
            "error_type": type(error).__name__,
            "error": str(error),
        }


def simulate_selected_units(
    units: pd.DataFrame,
    *,
    repo_root: str | Path,
    integration_manifest_path: str | Path,
    slim_path: str | Path,
    workers: int = 4,
    eas_max_draws: int = DEFAULT_EAS_MAX_DRAWS,
    draw_timeout_seconds: float = DEFAULT_DRAW_TIMEOUT_SECONDS,
) -> pd.DataFrame:
    """Run exactly 60 EAS units through isolated, resumable worker processes."""

    if int(workers) < 1 or int(workers) > 24:
        raise ValueError("selected production workers must be between 1 and 24")
    if int(eas_max_draws) != DEFAULT_EAS_MAX_DRAWS:
        raise ValueError("EAS selected retry budget is frozen at 68 draws")
    if not math.isclose(
        float(draw_timeout_seconds),
        DEFAULT_DRAW_TIMEOUT_SECONDS,
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        raise ValueError("EAS selected per-draw timeout is frozen at 30 minutes")
    bundle = load_selected_production_bundle(
        integration_manifest_path,
        repo_root=repo_root,
        slim_path=slim_path,
    )
    if int(eas_max_draws) != int(bundle.manifest["runtime"]["eas_max_draws"]):
        raise ValueError("EAS selected retry budget differs from authorization")
    canonical = pd.read_csv(bundle.execution_units_path, sep="\t")
    selected = _validate_requested_units(units.copy(), canonical)
    tasks = []
    for row in selected.to_dict(orient="records"):
        tasks.append(
            {
                "unit": row,
                "selected_max_draws": int(eas_max_draws),
                "draw_timeout_seconds": float(draw_timeout_seconds),
                "cumulative_timeout_seconds": (
                    float(draw_timeout_seconds) * int(eas_max_draws)
                ),
            }
        )
    configuration = _worker_configuration(bundle)
    rows: list[dict[str, Any]] = []
    if int(workers) == 1:
        _initialize_selected_worker(configuration)
        try:
            rows = [_simulate_selected_task(task) for task in tasks]
        finally:
            global _WORKER_BUNDLE, _WORKER_PATCH_CONTEXT
            if _WORKER_PATCH_CONTEXT is not None:
                _WORKER_PATCH_CONTEXT.__exit__(None, None, None)
            _WORKER_PATCH_CONTEXT = None
            _WORKER_BUNDLE = None
    else:
        with ProcessPoolExecutor(
            max_workers=int(workers),
            initializer=_initialize_selected_worker,
            initargs=(configuration,),
        ) as pool:
            futures = {
                pool.submit(_simulate_selected_task, task): task for task in tasks
            }
            for future in as_completed(futures):
                rows.append(future.result())
    return pd.DataFrame(rows).sort_values("unit_id").reset_index(drop=True)


def _validate_eas_attempt_provenance(
    bundle: SelectedProductionBundle,
    unit: Mapping[str, Any],
    completion: Mapping[str, Any],
) -> dict[str, Any]:
    """Validate durable retry, seed, supervisor, and no-orphan provenance."""

    from . import focused_han_direct_v3 as hardened_supervisor

    record = simulation._validate_unit(unit)
    unit_id = str(record["unit_id"])
    unit_dir = bundle.campaign_dir / "work" / unit_id
    attempts_path = unit_dir / "attempts.tsv"
    supervision_path = unit_dir / "supervision.tsv"
    if not attempts_path.is_file() or not supervision_path.is_file():
        raise ValueError("EAS selected supervisor attempt provenance is absent")
    attempts = pd.read_csv(attempts_path, sep="\t")
    supervision = pd.read_csv(supervision_path, sep="\t")
    copied_supervisor_fields = {
        "schema": "eas_supervisor_schema",
        "execution_id": "eas_supervisor_execution_id",
        "supervisor_status": "eas_supervisor_status",
        "supervisor_pid": "eas_supervisor_pid",
        "supervisor_process_identity": "eas_supervisor_process_identity",
        "supervisor_started_unix_ns": "eas_supervisor_started_unix_ns",
        "child_pid": "eas_child_pid",
        "child_pgid": "eas_child_pgid",
        "child_sid": "eas_child_sid",
        "child_process_identity": "eas_child_process_identity",
        "child_identity_observed": "eas_child_identity_observed",
        "engine_wrapper_pid": "eas_engine_wrapper_pid",
        "engine_wrapper_process_identity": "eas_engine_wrapper_process_identity",
        "parent_death_watchdog_armed": "eas_parent_death_watchdog_armed",
        "timeout_seconds": "eas_supervisor_timeout_seconds",
        "term_grace_seconds": "eas_supervisor_term_grace_seconds",
        "sigterm_sent": "eas_sigterm_sent",
        "sigkill_sent": "eas_sigkill_sent",
        "direct_child_reaped": "eas_direct_child_reaped",
        "adopted_children_reaped": "eas_adopted_children_reaped",
        "process_group_gone": "eas_process_group_gone",
        "temp_directory": "eas_attempt_temp_directory",
        "temp_tree_created": "eas_attempt_temp_tree_created",
        "temp_tree_loaded": "eas_attempt_temp_tree_loaded",
        "temp_tree_removed": "eas_attempt_temp_tree_removed",
        "supervisor_elapsed_seconds": "eas_supervisor_elapsed_seconds",
    }
    required_attempt_columns = {
        "attempt_zero_based",
        "seed",
        "status",
        "timeout_seconds",
        "elapsed_seconds",
        "accepted",
        "rejection_reason",
        "eas_planned_simulation_seed",
        "eas_planned_panel_seed",
        "eas_panel_seed_used",
        "eas_attempt_ledger_schema",
        "eas_supervisor_source_sha256",
        *copied_supervisor_fields.values(),
    }
    if required_attempt_columns.difference(attempts.columns):
        raise ValueError("EAS selected supervisor attempt columns differ")
    if copied_supervisor_fields.keys() - set(supervision.columns):
        raise ValueError("EAS selected process-supervision columns differ")

    def exact_integers(
        frame: pd.DataFrame, column: str, *, label: str, minimum: int
    ) -> pd.Series:
        values = pd.to_numeric(frame[column], errors="coerce")
        if (
            values.isna().any()
            or not np.isfinite(values.astype(float)).all()
            or not np.equal(values.astype(float), np.rint(values.astype(float))).all()
            or (values < minimum).any()
        ):
            raise ValueError(f"EAS selected {label} integer {column} differs")
        return values.astype(np.int64)

    def strict_booleans(frame: pd.DataFrame, column: str, *, label: str) -> pd.Series:
        values = frame[column]
        if values.dtype == bool:
            return values
        normalized = values.fillna("").astype(str).str.strip().str.casefold()
        if not normalized.isin({"true", "false"}).all():
            raise ValueError(f"EAS selected {label} boolean {column} differs")
        return normalized.eq("true")

    def normalized_values(values: pd.Series) -> list[str]:
        return [
            "" if pd.isna(value) else str(value).strip().casefold() for value in values
        ]

    for column, minimum in {
        "attempt_zero_based": 0,
        "seed": 1,
        "eas_planned_simulation_seed": 1,
        "eas_planned_panel_seed": 1,
        "eas_supervisor_pid": 1,
        "eas_supervisor_started_unix_ns": 1,
        "eas_adopted_children_reaped": 0,
    }.items():
        attempts[column] = exact_integers(
            attempts, column, label="attempt ledger", minimum=minimum
        )
    observed = list(attempts["attempt_zero_based"])
    if (
        not observed
        or observed != list(range(len(attempts)))
        or len(attempts) > DEFAULT_EAS_MAX_DRAWS
        or int(completion.get("attempts_completed", -1)) != len(attempts)
    ):
        raise ValueError("EAS selected attempt count/contiguity differs")
    schedule = _eas_supervisor_schedule(record).set_index("attempt_zero_based")
    source_sha256 = sha256_file(
        bundle.repo_root / hardened_supervisor.MODULE_SOURCE_PATH
    )
    expected_simulation_seeds = schedule.loc[observed, "simulation_seed"].astype(int)
    expected_panel_seeds = schedule.loc[observed, "panel_seed"].astype(int)
    if (
        not np.array_equal(attempts["seed"].to_numpy(), expected_simulation_seeds)
        or not np.array_equal(
            attempts["eas_planned_simulation_seed"].to_numpy(),
            expected_simulation_seeds,
        )
        or not np.array_equal(
            attempts["eas_planned_panel_seed"].to_numpy(), expected_panel_seeds
        )
        or set(attempts["eas_attempt_ledger_schema"].astype(str))
        != {EAS_ATTEMPT_LEDGER_SCHEMA}
        or set(attempts["eas_supervisor_source_sha256"].astype(str)) != {source_sha256}
    ):
        raise ValueError("EAS selected attempt seed/source binding differs")

    accepted = strict_booleans(attempts, "accepted", label="attempt ledger")
    panel_seed_used = strict_booleans(
        attempts, "eas_panel_seed_used", label="attempt ledger"
    )
    if int(accepted.sum()) != 1 or not bool(accepted.iloc[-1]):
        raise ValueError("EAS selected attempt ledger lacks one terminal acceptance")
    if not panel_seed_used.equals(accepted):
        raise ValueError("EAS selected panel-seed-used flags differ")

    observed_children = pd.to_numeric(
        attempts["eas_child_pid"], errors="coerce"
    ).notna()
    child_identity_observed = strict_booleans(
        attempts, "eas_child_identity_observed", label="attempt ledger"
    )
    direct_children_reaped = strict_booleans(
        attempts, "eas_direct_child_reaped", label="attempt ledger"
    )
    process_groups_gone = strict_booleans(
        attempts, "eas_process_group_gone", label="attempt ledger"
    )
    watchdog_armed = strict_booleans(
        attempts, "eas_parent_death_watchdog_armed", label="attempt ledger"
    )
    sigterm_sent = strict_booleans(attempts, "eas_sigterm_sent", label="attempt ledger")
    sigkill_sent = strict_booleans(attempts, "eas_sigkill_sent", label="attempt ledger")
    temp_created = strict_booleans(
        attempts, "eas_attempt_temp_tree_created", label="attempt ledger"
    )
    temp_loaded = strict_booleans(
        attempts, "eas_attempt_temp_tree_loaded", label="attempt ledger"
    )
    temp_removed = strict_booleans(
        attempts, "eas_attempt_temp_tree_removed", label="attempt ledger"
    )
    if (
        not observed_children.equals(child_identity_observed)
        or not direct_children_reaped[observed_children].all()
        or direct_children_reaped[~observed_children].any()
        or not process_groups_gone.all()
        or not watchdog_armed[observed_children].all()
        or watchdog_armed[~observed_children].any()
        or sigterm_sent[~observed_children].any()
        or sigkill_sent[~observed_children].any()
        or not temp_removed.all()
    ):
        raise ValueError("EAS selected child/process-group lifecycle differs")

    child_pid = pd.to_numeric(attempts["eas_child_pid"], errors="coerce")
    child_pgid = pd.to_numeric(attempts["eas_child_pgid"], errors="coerce")
    child_sid = pd.to_numeric(attempts["eas_child_sid"], errors="coerce")
    observed_child_values = pd.concat([child_pid, child_pgid, child_sid], axis=1)[
        observed_children
    ]
    if (
        observed_child_values.isna().any().any()
        or not np.equal(
            observed_child_values.astype(float),
            np.rint(observed_child_values.astype(float)),
        )
        .all()
        .all()
        or (observed_child_values < 1).any().any()
        or pd.concat([child_pid, child_pgid, child_sid], axis=1)[~observed_children]
        .notna()
        .any()
        .any()
        or not np.array_equal(
            child_pid[observed_children].astype(int),
            child_pgid[observed_children].astype(int),
        )
        or not np.array_equal(
            child_pid[observed_children].astype(int),
            child_sid[observed_children].astype(int),
        )
    ):
        raise ValueError("EAS selected child PID provenance differs")

    wrapper_pid = pd.to_numeric(attempts["eas_engine_wrapper_pid"], errors="coerce")
    wrapper_observed = wrapper_pid.notna()
    if (wrapper_pid[wrapper_observed] < 1).any() or not np.equal(
        wrapper_pid[wrapper_observed].astype(float),
        np.rint(wrapper_pid[wrapper_observed].astype(float)),
    ).all():
        raise ValueError("EAS selected engine-wrapper PID differs")

    status = attempts["status"].astype(str)
    timed_out = status.eq("timed_out")
    completed = status.eq("complete")
    rejection_reason = attempts["rejection_reason"].fillna("").astype(str)
    scientific_rejections = {
        "selected_focal_mutation_absent",
        "selected_focal_mutation_absent_background_only",
        "candidate_pool_af_outside_prespecified_band",
        "candidate_pool_has_no_exact_af_genotype_qc_panel",
    }
    interrupted_rejections = {
        "worker_interrupted_before_base_attempt_commit",
        "worker_interrupted_after_base_attempt_commit",
        "worker_interrupted_after_acceptance_before_unit_completion",
    }
    if (
        not status.isin({"complete", "timed_out"}).all()
        or not (completed[accepted] & rejection_reason[accepted].eq("")).all()
        or not rejection_reason[completed & ~accepted].isin(scientific_rejections).all()
        or not rejection_reason[timed_out]
        .isin({"draw_timeout", *interrupted_rejections})
        .all()
        or accepted[timed_out].any()
    ):
        raise ValueError("EAS selected attempt status/rejection matrix differs")

    supervisor_status = attempts["eas_supervisor_status"].astype(str)
    completed_supervisor = supervisor_status.eq("complete")
    allowed_timeout_supervisor = supervisor_status.isin(
        {
            "timed_out",
            "timed_out_after_terminal_message",
            "interrupted_by_worker_death",
            "startup_interrupted_before_child_observed",
            "engine_complete_without_base_attempt_commit",
        }
    )
    if (
        not (
            completed_supervisor.eq(~timed_out)
            & (completed_supervisor | allowed_timeout_supervisor)
        ).all()
        or not temp_created[completed_supervisor].all()
        or not temp_loaded[completed_supervisor].all()
    ):
        raise ValueError("EAS selected supervisor status matrix differs")

    base_elapsed = pd.to_numeric(attempts["elapsed_seconds"], errors="coerce")
    supervisor_elapsed = pd.to_numeric(
        attempts["eas_supervisor_elapsed_seconds"], errors="coerce"
    )
    supervisor_timeouts = pd.to_numeric(
        attempts["eas_supervisor_timeout_seconds"], errors="coerce"
    )
    base_timeouts = pd.to_numeric(attempts["timeout_seconds"], errors="coerce")
    if (
        base_elapsed.isna().any()
        or supervisor_elapsed.isna().any()
        or supervisor_timeouts.isna().any()
        or base_timeouts.isna().any()
        or not np.isfinite(base_elapsed.astype(float)).all()
        or not np.isfinite(supervisor_elapsed.astype(float)).all()
        or (base_elapsed < 0).any()
        or (supervisor_elapsed <= 0).any()
        or (base_elapsed + 1e-9 < supervisor_elapsed).any()
        or not np.allclose(supervisor_timeouts, base_timeouts, rtol=0.0, atol=1e-9)
        or (base_timeouts <= 0).any()
        or (base_timeouts > DEFAULT_DRAW_TIMEOUT_SECONDS + 1e-9).any()
        or set(attempts["eas_supervisor_term_grace_seconds"].astype(float))
        != {hardened_supervisor.SUPERVISOR_TERM_GRACE_SECONDS}
    ):
        raise ValueError("EAS selected elapsed/timeout accounting differs")

    supervision["attempt_zero_based"] = exact_integers(
        supervision,
        "attempt_zero_based",
        label="supervision ledger",
        minimum=0,
    )
    supervision["simulation_seed"] = exact_integers(
        supervision, "simulation_seed", label="supervision ledger", minimum=1
    )
    execution_ids = supervision["execution_id"].fillna("").astype(str)
    if (
        len(supervision) != len(attempts)
        or list(supervision["attempt_zero_based"].astype(int)) != observed
        or set(supervision["schema"].astype(str))
        != {hardened_supervisor.SUPERVISOR_LEDGER_SCHEMA}
        or not np.array_equal(
            supervision["simulation_seed"].to_numpy(), expected_simulation_seeds
        )
        or execution_ids.eq("").any()
        or execution_ids.duplicated().any()
        or not execution_ids.str.fullmatch(r"[0-9a-f]{64}").all()
    ):
        raise ValueError("EAS selected process-supervision ledger differs")

    for source_column, attempt_column in copied_supervisor_fields.items():
        if normalized_values(supervision[source_column]) != normalized_values(
            attempts[attempt_column]
        ):
            raise ValueError(f"EAS selected supervisor field {source_column} differs")

    reconciliation_path = unit_dir / "supervision_reconciliation.tsv"
    reconciliation_sha256: str | None = None
    interrupted_charge = 0.0
    reconciliation_rows = 0
    if reconciliation_path.is_file():
        reconciliation = pd.read_csv(reconciliation_path, sep="\t")
        required_reconciliation = {
            "schema",
            "reason",
            "execution_id",
            "attempt_zero_based",
            "simulation_seed",
            "prior_supervisor_pid",
            "prior_child_pid",
            "prior_child_pgid",
            "process_group_gone",
            "temp_directory_removed",
            "charged_elapsed_seconds",
            "transferred_to_attempt_elapsed_seconds",
            "atomic_temp_files_removed",
            "atomic_temp_file_count_removed",
            "reconciled_unix_ns",
        }
        if required_reconciliation.difference(reconciliation.columns):
            raise ValueError("EAS selected reconciliation fields differ")
        reconciliation["attempt_zero_based"] = exact_integers(
            reconciliation,
            "attempt_zero_based",
            label="reconciliation ledger",
            minimum=0,
        )
        reconciliation["simulation_seed"] = exact_integers(
            reconciliation,
            "simulation_seed",
            label="reconciliation ledger",
            minimum=1,
        )
        reconciliation["prior_supervisor_pid"] = exact_integers(
            reconciliation,
            "prior_supervisor_pid",
            label="reconciliation ledger",
            minimum=1,
        )
        reconciliation["atomic_temp_file_count_removed"] = exact_integers(
            reconciliation,
            "atomic_temp_file_count_removed",
            label="reconciliation ledger",
            minimum=0,
        )
        reconciliation["reconciled_unix_ns"] = exact_integers(
            reconciliation,
            "reconciled_unix_ns",
            label="reconciliation ledger",
            minimum=1,
        )
        charges = pd.to_numeric(
            reconciliation["charged_elapsed_seconds"], errors="coerce"
        )
        transferred = pd.to_numeric(
            reconciliation["transferred_to_attempt_elapsed_seconds"],
            errors="coerce",
        )
        reconciliation_attempts = reconciliation["attempt_zero_based"]
        reconciliation_execution_ids = (
            reconciliation["execution_id"].fillna("").astype(str)
        )
        reasons = set(reconciliation["reason"].astype(str))
        if (
            set(reconciliation["schema"].astype(str))
            != {hardened_supervisor.SUPERVISOR_RECONCILIATION_SCHEMA}
            or charges.isna().any()
            or transferred.isna().any()
            or not np.isfinite(charges.astype(float)).all()
            or not np.isfinite(transferred.astype(float)).all()
            or (charges < 0).any()
            or (transferred < 0).any()
            or not reasons.issubset(
                {
                    "worker_interrupted_active_engine_call",
                    "engine_complete_without_base_attempt_commit",
                    "attempt_interrupted_after_base_commit",
                    "accepted_attempt_interrupted_before_unit_completion",
                }
            )
            or (reconciliation_attempts >= len(attempts)).any()
            or reconciliation_attempts.duplicated().any()
            or reconciliation["reconciled_unix_ns"].duplicated().any()
            or reconciliation_execution_ids.eq("").any()
            or reconciliation_execution_ids.duplicated().any()
            or not reconciliation_execution_ids.str.fullmatch(r"[0-9a-f]{64}").all()
            or list(reconciliation["simulation_seed"])
            != [
                int(schedule.loc[attempt, "simulation_seed"])
                for attempt in reconciliation_attempts
            ]
            or not strict_booleans(
                reconciliation, "process_group_gone", label="reconciliation ledger"
            ).all()
            or not strict_booleans(
                reconciliation,
                "temp_directory_removed",
                label="reconciliation ledger",
            ).all()
        ):
            raise ValueError("EAS selected reconciliation ledger differs")

        special_dispositions: set[int] = set()
        supervisor_execution_by_attempt = {
            int(row["attempt_zero_based"]): str(row["execution_id"])
            for row in supervision.to_dict(orient="records")
        }
        for row_index, row in reconciliation.iterrows():
            attempt = int(row["attempt_zero_based"])
            reason = str(row["reason"])
            attempt_row = attempts.iloc[attempt]
            supervisor_row = supervision.iloc[attempt]
            prior_child_pid = pd.to_numeric(
                pd.Series([row.get("prior_child_pid")]), errors="coerce"
            ).iloc[0]
            prior_child_pgid = pd.to_numeric(
                pd.Series([row.get("prior_child_pgid")]), errors="coerce"
            ).iloc[0]
            expected_child_pid = pd.to_numeric(
                pd.Series([attempt_row.get("eas_child_pid")]), errors="coerce"
            ).iloc[0]
            expected_child_pgid = pd.to_numeric(
                pd.Series([attempt_row.get("eas_child_pgid")]), errors="coerce"
            ).iloc[0]
            if (
                str(row["execution_id"]) != supervisor_execution_by_attempt.get(attempt)
                or str(attempt_row["eas_supervisor_execution_id"])
                != str(row["execution_id"])
                or int(row["prior_supervisor_pid"])
                != int(attempt_row["eas_supervisor_pid"])
                or not (
                    (pd.isna(prior_child_pid) and pd.isna(expected_child_pid))
                    or int(prior_child_pid) == int(expected_child_pid)
                )
                or not (
                    (pd.isna(prior_child_pgid) and pd.isna(expected_child_pgid))
                    or int(prior_child_pgid) == int(expected_child_pgid)
                )
                or str(attempt_row["status"]) != "timed_out"
                or bool(accepted.iloc[attempt])
                or bool(panel_seed_used.iloc[attempt])
            ):
                raise ValueError("EAS selected reconciliation linkage differs")
            charge = float(charges.iloc[row_index])
            transfer = float(transferred.iloc[row_index])
            supervisor_lifecycle = str(supervisor_row["supervisor_status"])
            rejection = str(attempt_row["rejection_reason"])
            try:
                removed_atomic_temps = json.loads(str(row["atomic_temp_files_removed"]))
            except json.JSONDecodeError as error:
                raise ValueError(
                    "EAS selected reconciliation atomic-temp record differs"
                ) from error
            if (
                not isinstance(removed_atomic_temps, list)
                or any(not isinstance(name, str) for name in removed_atomic_temps)
                or len(removed_atomic_temps) != len(set(removed_atomic_temps))
                or len(removed_atomic_temps)
                != int(row["atomic_temp_file_count_removed"])
                or any(
                    not hardened_supervisor._is_unit_atomic_temp_name(name)
                    for name in removed_atomic_temps
                )
            ):
                raise ValueError(
                    "EAS selected reconciliation atomic-temp record differs"
                )
            if reason == "worker_interrupted_active_engine_call":
                valid = (
                    rejection == "worker_interrupted_before_base_attempt_commit"
                    and supervisor_lifecycle
                    in {
                        "interrupted_by_worker_death",
                        "startup_interrupted_before_child_observed",
                    }
                    and transfer > 0
                    and charge == 0
                    and math.isclose(
                        float(base_elapsed.iloc[attempt]),
                        transfer,
                        rel_tol=0.0,
                        abs_tol=1e-9,
                    )
                )
            elif reason == "engine_complete_without_base_attempt_commit":
                valid = (
                    rejection == "worker_interrupted_before_base_attempt_commit"
                    and supervisor_lifecycle
                    == "engine_complete_without_base_attempt_commit"
                    and transfer > 0
                    and charge == 0
                    and math.isclose(
                        float(base_elapsed.iloc[attempt]),
                        transfer,
                        rel_tol=0.0,
                        abs_tol=1e-9,
                    )
                )
            elif reason == "attempt_interrupted_after_base_commit":
                valid = (
                    rejection == "worker_interrupted_after_base_attempt_commit"
                    and supervisor_lifecycle == "interrupted_by_worker_death"
                    and transfer == 0
                    and charge > 0
                )
            else:
                valid = (
                    rejection
                    == "worker_interrupted_after_acceptance_before_unit_completion"
                    and supervisor_lifecycle == "interrupted_by_worker_death"
                    and transfer == 0
                    and charge > 0
                )
            if not valid:
                raise ValueError("EAS selected reconciliation reason matrix differs")
            special_dispositions.add(attempt)
        expected_special_dispositions = {
            index
            for index, value in enumerate(rejection_reason)
            if value in interrupted_rejections
        }
        if special_dispositions != expected_special_dispositions:
            raise ValueError("EAS selected reconciliation dispositions differ")
        interrupted_charge = float(charges.sum())
        reconciliation_rows = len(reconciliation)
        reconciliation_sha256 = sha256_file(reconciliation_path)
    elif rejection_reason.isin(interrupted_rejections).any():
        raise ValueError("EAS selected interrupted attempt lacks reconciliation")

    active = list(unit_dir.glob(".han_v3_active_*.json"))
    scratch = list(unit_dir.glob(".han_v3_attempt_*"))
    atomic_temps = hardened_supervisor._unit_atomic_temp_paths(unit_dir)
    if active or scratch or atomic_temps:
        raise ValueError("EAS selected completed unit retains active supervisor state")
    if simulation.unit_lock_is_held(unit_dir / ".simulation.lock"):
        raise ValueError("EAS selected completed unit lock is still held")

    for index, value in enumerate(attempts["eas_supervisor_process_identity"]):
        hardened_supervisor._validate_linux_process_identity(
            value, expected_pid=int(attempts.iloc[index]["eas_supervisor_pid"])
        )
    for index, value in enumerate(attempts["eas_child_process_identity"]):
        if bool(observed_children.iloc[index]):
            identity = hardened_supervisor._validate_linux_process_identity(
                value, expected_pid=int(child_pid.iloc[index])
            )
            if hardened_supervisor._same_linux_process(identity):
                raise ValueError("EAS selected supervised child identity is still live")
        elif hardened_supervisor._parse_linux_process_identity(value) is not None or (
            isinstance(value, str) and value.strip()
        ):
            raise ValueError("EAS selected unobserved child identity differs")
    for index, value in enumerate(attempts["eas_engine_wrapper_process_identity"]):
        if bool(wrapper_observed.iloc[index]):
            identity = hardened_supervisor._validate_linux_process_identity(
                value, expected_pid=int(wrapper_pid.iloc[index])
            )
            if hardened_supervisor._same_linux_process(identity):
                raise ValueError("EAS selected engine-wrapper identity is still live")
        elif hardened_supervisor._parse_linux_process_identity(value) is not None or (
            isinstance(value, str) and value.strip()
        ):
            raise ValueError("EAS selected unobserved engine-wrapper identity differs")

    committed_elapsed = float(base_elapsed.sum())
    cumulative_elapsed = committed_elapsed + interrupted_charge
    cumulative_budget = DEFAULT_EAS_MAX_DRAWS * DEFAULT_DRAW_TIMEOUT_SECONDS
    if (
        interrupted_charge >= cumulative_budget
        or cumulative_elapsed > cumulative_budget + DEFAULT_DRAW_TIMEOUT_SECONDS + 1e-6
    ):
        raise ValueError("EAS selected cumulative elapsed accounting differs")
    return {
        "attempts_completed": len(attempts),
        "accepted_attempt_zero_based": observed[-1],
        "attempts_sha256": sha256_file(attempts_path),
        "supervision_sha256": sha256_file(supervision_path),
        "all_process_groups_gone_and_children_reaped": True,
        "fixed_seed_schedule_validated": True,
        "supervisor_fields_equal_attempt_ledger": True,
        "one_terminal_disposition_per_attempted_seed": True,
        "crash_reconciliation_rows": reconciliation_rows,
        "reconciliation_sha256": reconciliation_sha256,
        "interrupted_elapsed_seconds_charged": interrupted_charge,
        "committed_elapsed_seconds": committed_elapsed,
        "accounted_cumulative_elapsed_seconds": cumulative_elapsed,
    }


def validate_eas_cached_units_read_only(
    bundle: SelectedProductionBundle,
) -> dict[str, Any]:
    """Fully validate all 60 EAS caches without creating or changing a file."""

    canonical = pd.read_csv(bundle.execution_units_path, sep="\t")
    units = _eas_production_units(canonical)
    rows: list[dict[str, Any]] = []
    with patched_selected_executor(bundle):
        for unit in units.to_dict(orient="records"):
            unit_dir = bundle.campaign_dir / "work" / str(unit["unit_id"])
            completion_path = unit_dir / "simulation_complete.json"
            if not completion_path.is_file():
                raise ValueError(
                    f"EAS selected completion is absent: {unit['unit_id']}"
                )
            before = {
                path.relative_to(unit_dir).as_posix(): (
                    path.stat().st_size,
                    path.stat().st_mtime_ns,
                )
                for path in unit_dir.rglob("*")
                if path.is_file()
            }
            artifacts = simulation._simulate_unit_unlocked(
                unit,
                bundle.repo_root,
                bundle.campaign_dir,
                slim_path=bundle.slim_path,
                selected_max_draws=DEFAULT_EAS_MAX_DRAWS,
                slim_scaling_factor=simulation.DEFAULT_SLIM_SCALING_FACTOR,
                slim_burn_in=simulation.DEFAULT_SLIM_BURN_IN,
                selected_draw_timeout_seconds=DEFAULT_DRAW_TIMEOUT_SECONDS,
                selected_cumulative_timeout_seconds=(
                    DEFAULT_EAS_MAX_DRAWS * DEFAULT_DRAW_TIMEOUT_SECONDS
                ),
            )
            if artifacts.cache_hit is not True:
                raise RuntimeError("EAS read-only validation unexpectedly simulated")
            after = {
                path.relative_to(unit_dir).as_posix(): (
                    path.stat().st_size,
                    path.stat().st_mtime_ns,
                )
                for path in unit_dir.rglob("*")
                if path.is_file()
            }
            if before != after:
                raise RuntimeError("EAS read-only validation changed cached artifacts")
            completion = _read_json(
                completion_path, label=f"EAS completion {unit['unit_id']}"
            )
            attempt_audit = _validate_eas_attempt_provenance(bundle, unit, completion)
            rows.append(
                {
                    "unit_id": str(unit["unit_id"]),
                    "completion_path": completion_path.relative_to(
                        bundle.repo_root
                    ).as_posix(),
                    "completion_sha256": sha256_file(completion_path),
                    "contract_sha256": str(completion["contract_sha256"]),
                    **attempt_audit,
                }
            )
    return {
        "units": 60,
        "cells": 6,
        "replicates_per_cell": 10,
        "cache_validated_read_only": True,
        "rows_canonical_sha256": _canonical_sha256(rows),
        "rows": rows,
    }


def _load_han_v3_stage_audit(
    path: Path,
    *,
    repo_root: Path,
    stage: str,
    expected_units: int,
    expected_replicates: Sequence[int],
    manifest_path: Path,
    seed_plan_path: Path,
    expected_unit_ids: set[str],
    expected_seed_plan_audit: Mapping[str, Any],
) -> dict[str, Any]:
    """Validate one signed Han-v3 stage audit without touching Han caches."""

    from . import focused_han_direct_v3 as han_v3

    payload = _read_json(path, label=f"Han v3 {stage} audit")
    _validate_hashed_payload(
        payload,
        schema=han_v3.STAGE_AUDIT_SCHEMA,
        status_key="status",
        expected_status="passed",
        label=f"Han v3 {stage} audit",
    )
    required_true = (
        "all_current_completion_budgets_match_authorization",
        "all_current_completion_counts_match_attempt_ledgers",
        "all_current_completion_sources_match_pinned_runtime",
        "all_current_completions_cache_validated_read_only",
        "all_rejections_persisted",
        "all_unit_authorizations_match_manifest_and_seed_plan",
        "all_units_strict_type1_identity",
        "all_units_use_exact_final_contract",
        "pilot_units_count_toward_final",
    )
    if (
        payload.get("stage") != stage
        or int(payload.get("expected_units", -1)) != expected_units
        or int(payload.get("accepted_units", -1)) != expected_units
        or list(payload.get("per_cell_replicates", ())) != list(expected_replicates)
        or payload.get("manifest_sha256") != sha256_file(manifest_path)
        or payload.get("seed_plan_sha256") != sha256_file(seed_plan_path)
        or payload.get("seed_plan_audit") != expected_seed_plan_audit
        or payload.get("failed_v2_inferential_use") is not False
        or payload.get("gamma_smc_statistics_used") is not False
        or any(payload.get(key) is not True for key in required_true)
    ):
        raise ValueError(f"Han v3 {stage} audit contract differs")
    rows = payload.get("rows")
    if not isinstance(rows, list) or len(rows) != expected_units:
        raise ValueError(f"Han v3 {stage} audit row coverage differs")
    unit_ids = [str(row.get("unit_id", "")) for row in rows if isinstance(row, Mapping)]
    if (
        len(unit_ids) != expected_units
        or len(set(unit_ids)) != expected_units
        or set(unit_ids) != expected_unit_ids
    ):
        raise ValueError(f"Han v3 {stage} audit unit IDs differ")
    return payload


def _load_han_v3_final_authorization(
    path: Path,
    *,
    manifest_path: Path,
    seed_plan_path: Path,
    manifest: Mapping[str, Any],
    final_audit: Mapping[str, Any],
) -> dict[str, Any]:
    """Require the final Han-v3 operational/recovery authorization."""

    from . import focused_han_direct_v3 as han_v3

    payload = _read_json(path, label="Han v3 final authorization")
    expected: dict[str, Any] = {
        "schema": han_v3.AUTHORIZATION_SCHEMA,
        "status": "authorized",
        "manifest_sha256": sha256_file(manifest_path),
        "seed_plan_sha256": sha256_file(seed_plan_path),
        "all_units_audit_payload_sha256": final_audit["payload_sha256"],
        "accepted_units": 60,
        "cells": 6,
        "replicates_per_cell": 10,
        "pilot_units_included": 6,
        "validated_completion_count": 60,
        "validated_attempt_count": sum(
            int(row["attempts_completed"]) for row in final_audit["rows"]
        ),
        "attempts_per_unit_budget": han_v3.ATTEMPTS_PER_UNIT,
        "draw_timeout_seconds": han_v3.DRAW_TIMEOUT_SECONDS,
        "cumulative_timeout_seconds": han_v3.CUMULATIVE_TIMEOUT_SECONDS,
        "all_completion_counts_budgets_and_authorizations_validated": True,
        "all_process_groups_gone_and_children_reaped": True,
        "pinned_runtime": manifest["runtime"],
        "direct_operational_estimand": True,
        "conditional_inference": True,
        "gamma_smc_statistics_used": False,
        "failed_v2_inferential_use": False,
        "scientific_caveats": manifest["scientific_caveats"],
    }
    expected["payload_sha256"] = _canonical_sha256(expected)
    if payload != expected:
        raise ValueError("Han v3 final operational authorization differs")
    return payload


def _han_v3_recovery_evidence_present(output_dir: Path) -> bool:
    names = (
        "recovery_history",
        "recovery_authorizations",
        "recovery_run_snapshots",
        "recovery_completion_audits",
        "han_direct_v3_dtype_recovery_authorization.json",
        "han_direct_v3_dtype_recovery_audit.json",
    )
    return any((output_dir / name).exists() for name in names)


def _load_han_v3_recovery_binding(
    *,
    repo_root: Path,
    slim_path: Path,
    han_manifest_path: Path,
    final_audit: Mapping[str, Any],
    authorization_path: Path,
    audit_path: Path,
) -> dict[str, Any]:
    """Strictly bind the additive dtype-recovery execution seam."""

    from . import focused_han_direct_v3_recovery as recovery

    bundle = recovery.load_recovery_bundle(
        han_manifest_path, repo_root=repo_root, slim_path=slim_path
    )
    authorization = recovery._load_authorization(bundle, authorization_path)
    audit = _read_json(audit_path, label="Han v3 recovery completion audit")
    _validate_hashed_payload(
        audit,
        schema=recovery.RECOVERY_COMPLETION_AUDIT_SCHEMA,
        status_key="status",
        expected_status="passed",
        label="Han v3 recovery completion audit",
    )
    raw_authorized_ids = [
        str(value) for value in authorization.get("recoverable_unit_ids", ())
    ]
    authorized_ids = sorted(raw_authorized_ids)
    recovered_ids = sorted(
        str(value) for value in audit.get("post_recovery_recovered_unit_ids", ())
    )
    pre_recovery_ids = sorted(
        str(value) for value in audit.get("pre_recovery_incomplete_unit_ids", ())
    )
    rows = audit.get("rows")
    row_ids = (
        sorted(str(row.get("unit_id", "")) for row in rows if isinstance(row, Mapping))
        if isinstance(rows, list)
        else []
    )
    process_local_patch = {
        "input": "attempts.tsv frame passed to frozen v3 callback",
        "column": "rejection_reason",
        "change": "cast existing or new column to object dtype",
        "planned_seed_change": False,
        "retry_change": False,
        "timeout_change": False,
        "budget_change": False,
        "reconciliation_change": False,
        "scientific_model_change": False,
    }
    history_path = _inside_repo(
        repo_root / str(authorization.get("history_manifest_path", "")),
        repo_root,
        label="Han v3 recovery history manifest",
    )
    history = _read_json(history_path, label="Han v3 recovery history manifest")
    _validate_hashed_payload(
        history,
        schema=recovery.RECOVERY_HISTORY_SCHEMA,
        status_key="status",
        expected_status="immutable_stopped_state_snapshot",
        label="Han v3 recovery history manifest",
    )
    history_files = history.get("files")
    if not isinstance(history_files, list):
        raise TypeError("Han v3 recovery history file inventory differs")
    history_snapshot_root = history_path.parent.resolve()
    observed_history_bytes = 0
    observed_history_paths: set[str] = set()
    for file_record in history_files:
        if not isinstance(file_record, Mapping):
            raise TypeError("Han v3 recovery history file record differs")
        relative_snapshot = str(file_record.get("snapshot_path", ""))
        snapshot = (history_snapshot_root / relative_snapshot).resolve()
        try:
            snapshot.relative_to(history_snapshot_root)
        except ValueError as error:
            raise ValueError(
                "Han v3 recovery snapshot escapes its history root"
            ) from error
        if relative_snapshot in observed_history_paths:
            raise ValueError("Han v3 recovery history snapshot path is duplicated")
        observed_history_paths.add(relative_snapshot)
        if (
            not snapshot.is_file()
            or snapshot.is_symlink()
            or snapshot.stat().st_size != int(file_record.get("size_bytes", -1))
            or sha256_file(snapshot) != str(file_record.get("sha256", ""))
        ):
            raise ValueError("Han v3 recovery history snapshot checksum differs")
        observed_history_bytes += snapshot.stat().st_size
    authorized_stage_ids = set(
        recovery.v3._stage_units(bundle.direct, recovery.v3.STAGE_SIMULATE)[
            "unit_id"
        ].astype(str)
    )
    canonical_authorization_relative = authorization_path.relative_to(
        repo_root
    ).as_posix()
    expected_true = (
        audit.get(
            "original_attempt_seed_timeout_budget_and_reconciliation_semantics_preserved"
        )
        is True
    )
    if (
        len(raw_authorized_ids) != recovery.EXPECTED_STOPPED_RECOVERABLE_UNITS
        or len(set(raw_authorized_ids)) != len(raw_authorized_ids)
        or raw_authorized_ids != sorted(raw_authorized_ids)
        or any(unit_id not in authorized_stage_ids for unit_id in authorized_ids)
        or recovered_ids != authorized_ids
        or pre_recovery_ids != authorized_ids
        or row_ids != authorized_ids
        or int(audit.get("initial_authorized_recovery_units", -1))
        != len(authorized_ids)
        or int(audit.get("validated_recovered_completions", -1)) != len(authorized_ids)
        or int(audit.get("remaining_incomplete_units", -1)) != 0
        or audit.get("post_recovery_incomplete_unit_ids") != []
        or audit.get("manifest_sha256") != sha256_file(han_manifest_path)
        or audit.get("seed_plan_sha256") != sha256_file(bundle.direct.seed_plan_path)
        or audit.get("original_implementation_sources")
        != bundle.direct.manifest["implementation_sources"]
        or audit.get("recovery_sources") != dict(bundle.recovery_sources)
        or authorization.get("process_local_patch") != process_local_patch
        or audit.get("dry_plan_payload_sha256")
        != authorization.get("dry_plan_payload_sha256")
        or audit.get("history_manifest_path")
        != authorization.get("history_manifest_path")
        or audit.get("history_manifest_sha256")
        != authorization.get("history_manifest_sha256")
        or audit.get("history_payload_sha256")
        != authorization.get("history_payload_sha256")
        or sha256_file(history_path) != authorization.get("history_manifest_sha256")
        or history.get("payload_sha256") != authorization.get("history_payload_sha256")
        or history.get("dry_plan_payload_sha256")
        != authorization.get("dry_plan_payload_sha256")
        or int(history.get("recoverable_units", -1)) != len(authorized_ids)
        or int(history.get("file_count", -1)) != len(history_files)
        or int(history.get("total_bytes", -1)) != observed_history_bytes
        or history.get("manifest_sha256") != audit.get("manifest_sha256")
        or history.get("seed_plan_sha256") != audit.get("seed_plan_sha256")
        or history.get("original_implementation_sources")
        != audit.get("original_implementation_sources")
        or history.get("recovery_sources") != audit.get("recovery_sources")
        or audit.get("canonical_authorization_path") != canonical_authorization_relative
        or audit.get("canonical_authorization_sha256")
        != sha256_file(authorization_path)
        or audit.get("post_run_base_audit_status") != "passed"
        or int(audit.get("post_run_base_audit_accepted_units", -1)) != 60
        or audit.get("post_run_base_audit_payload_sha256")
        != final_audit["payload_sha256"]
        or audit.get("post_run_base_audit") != final_audit
        or not expected_true
    ):
        raise ValueError("Han v3 dtype-recovery audit binding differs")
    for row in rows:
        unit_id = str(row["unit_id"])
        completion = (
            bundle.direct.campaign_dir / "work" / unit_id / "simulation_complete.json"
        )
        if not completion.is_file() or str(
            row.get("completion_sha256", "")
        ) != sha256_file(completion):
            raise ValueError(f"Han v3 recovered completion checksum differs: {unit_id}")
    return {
        "authorization": {
            **_file_record(authorization_path, repo_root),
            "payload_sha256": authorization["payload_sha256"],
        },
        "audit": {
            **_file_record(audit_path, repo_root),
            "payload_sha256": audit["payload_sha256"],
        },
        "recovered_units": len(authorized_ids),
        "recovered_unit_ids_sha256": _canonical_sha256(authorized_ids),
        "history": {
            **_file_record(history_path, repo_root),
            "payload_sha256": history["payload_sha256"],
        },
        "process_local_patch_only": True,
        "base_final_audit_recomputed_after_recovery": True,
    }


def build_aggregate_readiness(
    *,
    repo_root: str | Path,
    eas_integration_manifest_path: str | Path,
    slim_path: str | Path,
    han_v3_manifest_path: str | Path = DEFAULT_HAN_V3_MANIFEST_RELATIVE_PATH,
    han_v3_pilot_audit_path: str | Path = DEFAULT_HAN_V3_PILOT_AUDIT_RELATIVE_PATH,
    han_v3_final_audit_path: str | Path = DEFAULT_HAN_V3_FINAL_AUDIT_RELATIVE_PATH,
    han_v3_authorization_path: str | Path = DEFAULT_HAN_V3_AUTHORIZATION_RELATIVE_PATH,
    han_v3_recovery_authorization_path: str | Path = (
        DEFAULT_HAN_V3_RECOVERY_AUTHORIZATION_RELATIVE_PATH
    ),
    han_v3_recovery_audit_path: str
    | Path = DEFAULT_HAN_V3_RECOVERY_AUDIT_RELATIVE_PATH,
    require_final: bool = True,
    validate_eas_cache: bool = True,
) -> dict[str, Any]:
    """Bind current EAS and Han-v3 evidence for standard aggregation."""

    from . import focused_han_direct_v3 as han_v3

    root = Path(repo_root).resolve()
    if any("onedrive" in part.casefold() for part in root.parts):
        raise ValueError("selected aggregate repo_root must not be inside OneDrive")
    slim = _inside_repo(slim_path, root, label="SLiM binary")
    eas_manifest = _inside_repo(
        eas_integration_manifest_path, root, label="EAS integration manifest"
    )
    eas_bundle = load_selected_production_bundle(
        eas_manifest, repo_root=root, slim_path=slim
    )
    eas_cache = (
        validate_eas_cached_units_read_only(eas_bundle)
        if validate_eas_cache
        else {
            "units": 60,
            "cache_validated_read_only": False,
            "validation_deferred_for_status_only": True,
        }
    )

    han_manifest = _inside_repo(han_v3_manifest_path, root, label="Han v3 manifest")
    if (
        han_manifest.relative_to(root).as_posix()
        != DEFAULT_HAN_V3_MANIFEST_RELATIVE_PATH
    ):
        raise ValueError("Han v3 canonical manifest path differs")
    han_bundle = han_v3.load_bundle(han_manifest, repo_root=root, slim_path=slim)
    pilot_path = _inside_repo(han_v3_pilot_audit_path, root, label="Han v3 pilot audit")
    if (
        pilot_path.relative_to(root).as_posix()
        != DEFAULT_HAN_V3_PILOT_AUDIT_RELATIVE_PATH
    ):
        raise ValueError("Han v3 canonical pilot-audit path differs")
    pilot = _load_han_v3_stage_audit(
        pilot_path,
        repo_root=root,
        stage=han_v3.STAGE_PILOT,
        expected_units=6,
        expected_replicates=[han_v3.PILOT_REPLICATE_INDEX],
        manifest_path=han_manifest,
        seed_plan_path=han_bundle.seed_plan_path,
        expected_unit_ids=set(
            han_bundle.units[
                han_bundle.units["replicate_index"]
                .astype(int)
                .eq(han_v3.PILOT_REPLICATE_INDEX)
            ]["unit_id"].astype(str)
        ),
        expected_seed_plan_audit=han_v3._seed_plan_audit(han_bundle),
    )
    if han_v3.audit_stage(han_bundle, han_v3.STAGE_PILOT, write=False) != pilot:
        raise ValueError("Han v3 pilot audit differs from current validated caches")
    final_path = _inside_repo(
        han_v3_final_audit_path, root, label="Han v3 all-unit audit"
    )
    if (
        final_path.relative_to(root).as_posix()
        != DEFAULT_HAN_V3_FINAL_AUDIT_RELATIVE_PATH
    ):
        raise ValueError("Han v3 canonical all-unit-audit path differs")
    final: dict[str, Any] | None = None
    final_authorization: dict[str, Any] | None = None
    authorization_path = _inside_repo(
        han_v3_authorization_path, root, label="Han v3 final authorization"
    )
    if (
        authorization_path.relative_to(root).as_posix()
        != DEFAULT_HAN_V3_AUTHORIZATION_RELATIVE_PATH
    ):
        raise ValueError("Han v3 canonical final-authorization path differs")
    if final_path.is_file():
        final = _load_han_v3_stage_audit(
            final_path,
            repo_root=root,
            stage=han_v3.STAGE_ALL,
            expected_units=60,
            expected_replicates=list(range(1, han_v3.REPLICATES_PER_CELL + 1)),
            manifest_path=han_manifest,
            seed_plan_path=han_bundle.seed_plan_path,
            expected_unit_ids=set(han_bundle.units["unit_id"].astype(str)),
            expected_seed_plan_audit=han_v3._seed_plan_audit(han_bundle),
        )
        if han_v3.audit_stage(han_bundle, han_v3.STAGE_SIMULATE, write=False) != final:
            raise ValueError(
                "Han v3 all-unit audit differs from current validated caches"
            )
        if authorization_path.is_file():
            final_authorization = _load_han_v3_final_authorization(
                authorization_path,
                manifest_path=han_manifest,
                seed_plan_path=han_bundle.seed_plan_path,
                manifest=han_bundle.manifest,
                final_audit=final,
            )
        elif require_final:
            raise ValueError(
                "Han v3 final operational authorization is required for aggregation"
            )
    elif require_final:
        raise ValueError("Han v3 all-unit audit is required for aggregation")

    recovery_evidence_present = _han_v3_recovery_evidence_present(han_bundle.output_dir)
    recovery_binding: dict[str, Any] | None = None
    if recovery_evidence_present:
        recovery_authorization_path = _inside_repo(
            han_v3_recovery_authorization_path,
            root,
            label="Han v3 recovery authorization",
        )
        recovery_audit_path = _inside_repo(
            han_v3_recovery_audit_path, root, label="Han v3 recovery audit"
        )
        if (
            recovery_authorization_path.relative_to(root).as_posix()
            != DEFAULT_HAN_V3_RECOVERY_AUTHORIZATION_RELATIVE_PATH
            or recovery_audit_path.relative_to(root).as_posix()
            != DEFAULT_HAN_V3_RECOVERY_AUDIT_RELATIVE_PATH
        ):
            raise ValueError("Han v3 canonical dtype-recovery paths differ")
        if (
            not recovery_authorization_path.is_file()
            or not recovery_audit_path.is_file()
        ):
            raise ValueError(
                "Han v3 recovery evidence requires canonical authorization and audit"
            )
        if final is None:
            raise ValueError("Han v3 recovery evidence requires the final base audit")
        recovery_binding = _load_han_v3_recovery_binding(
            repo_root=root,
            slim_path=slim,
            han_manifest_path=han_manifest,
            final_audit=final,
            authorization_path=recovery_authorization_path,
            audit_path=recovery_audit_path,
        )

    ready = (
        final is not None
        and final_authorization is not None
        and (not recovery_evidence_present or recovery_binding is not None)
        and eas_cache.get("cache_validated_read_only") is True
    )
    status = (
        "ready_for_aggregate"
        if ready
        else (
            "final_audit_passed_authorization_pending"
            if final is not None
            else "pilot_passed_final_pending"
        )
    )
    payload: dict[str, Any] = {
        "schema": AGGREGATE_READINESS_SCHEMA,
        "status": status,
        "eas": {
            "integration_manifest": _file_record(eas_manifest, root),
            "integration_payload_sha256": eas_bundle.manifest["payload_sha256"],
            "cache_audit": eas_cache,
        },
        "han_v3": {
            "manifest": _file_record(han_manifest, root),
            "manifest_payload_sha256": han_bundle.manifest["payload_sha256"],
            "pilot_audit": {
                **_file_record(pilot_path, root),
                "payload_sha256": pilot["payload_sha256"],
            },
            "all_units_audit": (
                {
                    **_file_record(final_path, root),
                    "payload_sha256": final["payload_sha256"],
                }
                if final is not None
                else None
            ),
            "final_operational_authorization": (
                {
                    **_file_record(authorization_path, root),
                    "payload_sha256": final_authorization["payload_sha256"],
                }
                if final_authorization is not None
                else None
            ),
            "dtype_recovery": {
                "evidence_present": recovery_evidence_present,
                "binding": recovery_binding,
            },
        },
        "counts": {
            "eas_selected": 60,
            "han_selected": 60 if final is not None else 6,
            "selected_ready": 120 if ready else 60,
            "neutral_preserved": 600,
            "campaign_total_when_ready": 720,
        },
        "aggregate_contract": {
            "standard_han_v3_cache_dispatch_required": True,
            "eas_selected_patch_nested_only_for_eas": True,
            "neutral_cache_dispatch_unchanged": True,
            "gamma_smc_statistics_used_for_authorization": False,
        },
    }
    payload["payload_sha256"] = _canonical_sha256(payload)
    if require_final and not ready:
        raise ValueError("selected campaign is not ready for aggregation")
    return payload


def write_aggregate_readiness(path: str | Path, payload: Mapping[str, Any]) -> Path:
    """Publish only a complete, immutable selected aggregate authorization."""

    _validate_hashed_payload(
        payload,
        schema=AGGREGATE_READINESS_SCHEMA,
        status_key="status",
        expected_status="ready_for_aggregate",
        label="selected aggregate readiness",
    )
    destination = Path(path).resolve()
    if any("onedrive" in part.casefold() for part in destination.parts):
        raise ValueError("selected aggregate readiness must not be inside OneDrive")
    content = (
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        if destination.read_bytes() != content:
            raise FileExistsError(f"aggregate readiness already differs: {destination}")
        return destination
    temporary = destination.with_name(f".{destination.name}.{os.getpid()}.tmp")
    temporary.write_bytes(content)
    os.replace(temporary, destination)
    return destination


def load_aggregate_readiness(
    path: str | Path,
    *,
    repo_root: str | Path,
    slim_path: str | Path,
) -> dict[str, Any]:
    """Rebuild and compare the complete readiness contract before aggregation."""

    root = Path(repo_root).resolve()
    readiness_path = _inside_repo(path, root, label="aggregate readiness")
    if (
        readiness_path.relative_to(root).as_posix()
        != DEFAULT_AGGREGATE_READINESS_RELATIVE_PATH
    ):
        raise ValueError("selected aggregate-readiness path differs")
    payload = _read_json(readiness_path, label="selected aggregate readiness")
    _validate_hashed_payload(
        payload,
        schema=AGGREGATE_READINESS_SCHEMA,
        status_key="status",
        expected_status="ready_for_aggregate",
        label="selected aggregate readiness",
    )
    eas_manifest = _inside_repo(
        root
        / str(payload.get("eas", {}).get("integration_manifest", {}).get("path", "")),
        root,
        label="readiness EAS manifest",
    )
    han_manifest = _inside_repo(
        root / str(payload.get("han_v3", {}).get("manifest", {}).get("path", "")),
        root,
        label="readiness Han manifest",
    )
    pilot = _inside_repo(
        root / str(payload.get("han_v3", {}).get("pilot_audit", {}).get("path", "")),
        root,
        label="readiness Han pilot audit",
    )
    final_record = payload.get("han_v3", {}).get("all_units_audit")
    if not isinstance(final_record, Mapping):
        raise TypeError("selected aggregate readiness lacks Han all-unit audit")
    final = _inside_repo(
        root / str(final_record.get("path", "")),
        root,
        label="readiness Han all-unit audit",
    )
    authorization_record = payload.get("han_v3", {}).get(
        "final_operational_authorization"
    )
    if not isinstance(authorization_record, Mapping):
        raise TypeError("selected aggregate readiness lacks Han final authorization")
    authorization = _inside_repo(
        root / str(authorization_record.get("path", "")),
        root,
        label="readiness Han final authorization",
    )
    recovery = payload.get("han_v3", {}).get("dtype_recovery", {})
    recovery_authorization: Path | None = None
    recovery_audit: Path | None = None
    if recovery.get("evidence_present") is True:
        binding = recovery.get("binding")
        if not isinstance(binding, Mapping):
            raise TypeError("selected readiness lacks Han recovery binding")
        recovery_authorization = _inside_repo(
            root / str(binding.get("authorization", {}).get("path", "")),
            root,
            label="readiness Han recovery authorization",
        )
        recovery_audit = _inside_repo(
            root / str(binding.get("audit", {}).get("path", "")),
            root,
            label="readiness Han recovery audit",
        )
    current = build_aggregate_readiness(
        repo_root=root,
        eas_integration_manifest_path=eas_manifest,
        slim_path=slim_path,
        han_v3_manifest_path=han_manifest,
        han_v3_pilot_audit_path=pilot,
        han_v3_final_audit_path=final,
        han_v3_authorization_path=authorization,
        han_v3_recovery_authorization_path=(
            recovery_authorization
            or root / DEFAULT_HAN_V3_RECOVERY_AUTHORIZATION_RELATIVE_PATH
        ),
        han_v3_recovery_audit_path=(
            recovery_audit or root / DEFAULT_HAN_V3_RECOVERY_AUDIT_RELATIVE_PATH
        ),
        require_final=True,
        validate_eas_cache=True,
    )
    if current != payload:
        raise ValueError("selected aggregate readiness differs from current evidence")
    return payload


def aggregate_with_readiness(
    *,
    repo_root: str | Path,
    slim_path: str | Path,
    readiness_path: str | Path = DEFAULT_AGGREGATE_READINESS_RELATIVE_PATH,
) -> tuple[Path, Path, Path, Path, Path]:
    """Run standard aggregation only through the validated EAS plus Han-v3 dispatch."""

    from . import focused_han_direct_v3 as han_v3

    root = Path(repo_root).resolve()
    readiness = load_aggregate_readiness(
        readiness_path, repo_root=root, slim_path=slim_path
    )
    eas_manifest = root / readiness["eas"]["integration_manifest"]["path"]
    han_manifest = root / readiness["han_v3"]["manifest"]["path"]
    han_bundle = han_v3.load_bundle(han_manifest, repo_root=root, slim_path=slim_path)
    units = pd.read_csv(han_bundle.campaign_dir / "execution_units.tsv", sep="\t")
    _validate_execution_units(units)
    return han_v3.aggregate_decoded_outputs_with_v3_dispatch(
        han_bundle,
        units,
        selected_integration_manifest=eas_manifest,
    )


def analyze_with_readiness(
    *,
    repo_root: str | Path,
    slim_path: str | Path,
    readiness_path: str | Path = DEFAULT_AGGREGATE_READINESS_RELATIVE_PATH,
    empirical_path: str | Path | None = None,
) -> dict[str, Any]:
    """Aggregate, analyze, compare, and report without leaving v3 dispatch."""

    from .focused_selection_analysis import (
        run_focused_analysis,
        run_truth_gamma_comparison,
    )
    from .focused_selection_report import generate_run_report

    root = Path(repo_root).resolve()
    campaign = root / DEFAULT_CAMPAIGN_RELATIVE_PATH
    classes, pairs, truth, spatial, truth_spatial = aggregate_with_readiness(
        repo_root=root,
        slim_path=slim_path,
        readiness_path=readiness_path,
    )
    empirical = (
        None
        if empirical_path is None
        else _inside_repo(empirical_path, root, label="empirical comparison table")
    )
    analysis = run_focused_analysis(
        classes,
        campaign / "results" / "analysis",
        pair_summaries_path=pairs,
        spatial_summaries_path=spatial,
        empirical_path=empirical,
    )
    truth_analysis = run_focused_analysis(
        truth,
        campaign / "results" / "analysis_tree_truth",
        spatial_summaries_path=truth_spatial,
        minimum_selected_per_cell=campaign_plan.DEFAULT_SELECTED_REPLICATES,
        minimum_neutral_per_cell=campaign_plan.DEFAULT_NEUTRAL_REPLICATES,
    )
    comparison = run_truth_gamma_comparison(
        truth,
        classes,
        campaign / "results" / "truth_gamma_comparison",
    )
    report = generate_run_report(campaign)

    def relative_paths(paths: Mapping[str, Path]) -> dict[str, str]:
        return {
            str(label): Path(path).resolve().relative_to(root).as_posix()
            for label, path in paths.items()
        }

    return {
        "status": "analyzed_and_reported",
        "aggregate_outputs": [
            path.relative_to(root).as_posix()
            for path in (classes, pairs, truth, spatial, truth_spatial)
        ],
        "analysis": relative_paths(analysis),
        "truth_analysis": relative_paths(truth_analysis),
        "truth_gamma_comparison": relative_paths(comparison),
        "report": Path(report).resolve().relative_to(root).as_posix(),
        "all_simulation_cache_validation_dispatched_through_readiness": True,
    }


def report_with_readiness(
    *,
    repo_root: str | Path,
    slim_path: str | Path,
    readiness_path: str | Path = DEFAULT_AGGREGATE_READINESS_RELATIVE_PATH,
) -> dict[str, Any]:
    """Generate the report only after revalidating current aggregate readiness."""

    from .focused_selection_report import generate_run_report

    root = Path(repo_root).resolve()
    readiness = load_aggregate_readiness(
        readiness_path, repo_root=root, slim_path=slim_path
    )
    report = Path(generate_run_report(root / DEFAULT_CAMPAIGN_RELATIVE_PATH)).resolve()
    return {
        "status": "reported",
        "readiness_payload_sha256": readiness["payload_sha256"],
        "report": report.relative_to(root).as_posix(),
        "report_sha256": sha256_file(report),
        "simulation_cache_validation_dispatched_through_readiness": True,
        "unsafe_standard_reaggregation_used": False,
    }


def _inventory_directory(directory: Path, repo_root: Path) -> dict[str, Any]:
    files: list[dict[str, Any]] = []
    empty_directories: list[str] = []
    for path in sorted(directory.rglob("*"), key=lambda value: value.as_posix()):
        if path.is_symlink():
            raise ValueError(f"quarantine source contains a symlink: {path}")
        relative = path.relative_to(directory).as_posix()
        if path.is_dir():
            if not any(path.iterdir()):
                empty_directories.append(relative)
            continue
        if not path.is_file():
            raise ValueError(f"quarantine source contains an unsupported entry: {path}")
        files.append(
            {
                "path": relative,
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    record = {
        "source": directory.relative_to(repo_root).as_posix(),
        "files": files,
        "empty_directories": empty_directories,
    }
    return {
        "file_count": len(files),
        "total_size_bytes": sum(int(item["size_bytes"]) for item in files),
        "inventory_sha256": _canonical_sha256(record),
        "inventory": record,
    }


def _contract_implementation_sources(directory: Path) -> dict[str, str]:
    candidates = (
        directory / "simulation_contract.json",
        directory / "simulation_complete.json",
    )
    for path in candidates:
        if not path.is_file():
            continue
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
            contract = payload.get("contract", payload)
            sources = contract["implementation"]["sources"]
            if not isinstance(sources, Mapping):
                return {}
            return {str(key): str(value) for key, value in sources.items()}
        except (OSError, json.JSONDecodeError, KeyError, TypeError):
            return {}
    return {}


def build_quarantine_plan(
    *,
    repo_root: str | Path,
    campaign_dir: str | Path,
    integration_manifest_path: str | Path = DEFAULT_INTEGRATION_RELATIVE_PATH,
    expected_stale_count: int = DEFAULT_EXPECTED_STALE_SELECTED_DIRS,
) -> dict[str, Any]:
    """Inventory stale EAS selected directories without moving anything."""

    root = Path(repo_root).resolve()
    campaign = _inside_repo(campaign_dir, root, label="campaign directory")
    units_path = _inside_repo(
        campaign / "execution_units.tsv", root, label="execution units"
    )
    units = pd.read_csv(units_path, sep="\t")
    _validate_execution_units(units)
    integration_relative = (
        _inside_repo(integration_manifest_path, root, label="integration manifest")
        .relative_to(root)
        .as_posix()
    )
    required_source_paths = {MODULE_SOURCE_PATH, integration_relative}
    required_source_sha256 = {
        relative: sha256_file(root / relative)
        for relative in required_source_paths
        if (root / relative).is_file()
    }
    selected = units[
        units["simulation_class"].astype(str).eq("selected")
        & units["demography_id"].astype(str).eq("eas_phlash_median")
    ]
    work = _inside_repo(campaign / "work", root, label="campaign work")
    destination_root = _inside_repo(
        campaign / "quarantine" / "pre_focused_eas_selected_production",
        root,
        label="quarantine destination",
    )
    rows: list[dict[str, Any]] = []
    inventories: dict[str, Any] = {}
    for unit in selected.to_dict(orient="records"):
        unit_id = str(unit["unit_id"])
        source = _inside_repo(work / unit_id, root, label="selected unit directory")
        if not source.is_dir() or not any(source.iterdir()):
            continue
        old_sources = _contract_implementation_sources(source)
        if set(required_source_sha256) == required_source_paths and all(
            old_sources.get(relative) == expected
            for relative, expected in required_source_sha256.items()
        ):
            continue
        if simulation.unit_lock_is_held(source / ".simulation.lock"):
            raise ValueError(f"selected unit lock is currently held: {unit_id}")
        destination = _inside_repo(
            destination_root / unit_id, root, label="quarantine unit destination"
        )
        if destination.exists():
            raise ValueError(f"quarantine destination already exists: {destination}")
        if source.anchor.casefold() != destination.anchor.casefold():
            raise ValueError("quarantine source and destination are not on one volume")
        inventory = _inventory_directory(source, root)
        inventories[unit_id] = inventory["inventory"]
        rows.append(
            {
                "unit_id": unit_id,
                "simulation_class": "selected",
                "demography_id": str(unit["demography_id"]),
                "source": source.relative_to(root).as_posix(),
                "destination": destination.relative_to(root).as_posix(),
                "operation": "same_volume_rename_not_executed",
                "stale_reason": (
                    "eas_selected_directory_predates_frozen_production_adapter_binding"
                ),
                "completion_present": (source / "simulation_complete.json").is_file(),
                "partial_contract_present": (
                    source / "simulation_contract.json"
                ).is_file(),
                "file_count": inventory["file_count"],
                "total_size_bytes": inventory["total_size_bytes"],
                "inventory_sha256": inventory["inventory_sha256"],
                "lock_held": False,
            }
        )
    rows.sort(key=lambda row: row["unit_id"])
    if len(rows) != int(expected_stale_count):
        raise ValueError(
            "stale planned EAS selected-directory count differs: "
            f"observed {len(rows)}, expected {int(expected_stale_count)}"
        )
    if any(row["simulation_class"] != "selected" for row in rows):
        raise RuntimeError("neutral unit entered selected quarantine plan")
    if any(row["demography_id"] != "eas_phlash_median" for row in rows):
        raise RuntimeError("Han unit entered EAS selected quarantine plan")
    payload: dict[str, Any] = {
        "schema": QUARANTINE_PLAN_SCHEMA,
        "status": "planned_not_executed",
        "execution_supported_by_this_module": False,
        "scope": "stale_eas_selected_directories_only",
        "moves_performed": 0,
        "deletions_performed": 0,
        "expected_stale_selected_directories": int(expected_stale_count),
        "observed_stale_selected_directories": len(rows),
        "execution_units": _file_record(units_path, root),
        "required_new_contract_sources": {
            relative: required_source_sha256.get(relative, "not_yet_published")
            for relative in sorted(required_source_paths)
        },
        "rows": rows,
        "inventories": inventories,
        "neutral_directories_in_scope": False,
        "han_directories_in_scope": False,
    }
    payload["payload_sha256"] = _canonical_sha256(payload)
    return payload


def write_quarantine_plan(
    output_dir: str | Path, payload: Mapping[str, Any]
) -> tuple[Path, Path]:
    """Write planning artifacts only; this function cannot move unit data."""

    _validate_hashed_payload(
        payload,
        schema=QUARANTINE_PLAN_SCHEMA,
        status_key="status",
        expected_status="planned_not_executed",
        label="selected quarantine plan",
    )
    directory = Path(output_dir).resolve()
    if any("onedrive" in part.casefold() for part in directory.parts):
        raise ValueError("selected quarantine plan must not be inside OneDrive")
    directory.mkdir(parents=True, exist_ok=True)
    json_path = directory / "selected_quarantine_plan.json"
    table_path = directory / "selected_quarantine_plan.tsv"
    json_bytes = (
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")
    table_bytes = (
        pd.DataFrame(payload["rows"])
        .to_csv(sep="\t", index=False, lineterminator="\n")
        .encode("utf-8")
    )
    for path, content in ((json_path, json_bytes), (table_path, table_bytes)):
        if path.exists():
            if path.read_bytes() != content:
                raise FileExistsError(f"quarantine plan differs: {path}")
            continue
        temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
        temporary.write_bytes(content)
        os.replace(temporary, path)
    return json_path, table_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "action",
        choices=(
            "build-manifest",
            "validate",
            "quarantine-plan",
            "simulate-eas",
            "aggregate-status",
            "aggregate-validate",
            "aggregate",
            "analyze",
            "report",
        ),
    )
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--campaign-dir", type=Path)
    parser.add_argument("--slim-bin", type=Path)
    parser.add_argument("--integration-manifest", type=Path)
    parser.add_argument("--eas-adapter-authorization", type=Path)
    parser.add_argument("--eas-frozen", type=Path)
    parser.add_argument("--han-v3-manifest", type=Path)
    parser.add_argument("--han-v3-pilot-audit", type=Path)
    parser.add_argument("--han-v3-final-audit", type=Path)
    parser.add_argument("--han-v3-authorization", type=Path)
    parser.add_argument("--han-v3-recovery-authorization", type=Path)
    parser.add_argument("--han-v3-recovery-audit", type=Path)
    parser.add_argument("--aggregate-readiness", type=Path)
    parser.add_argument("--empirical-tsv", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument(
        "--expected-stale-count",
        type=int,
        default=DEFAULT_EXPECTED_STALE_SELECTED_DIRS,
    )
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument(
        "--status-skip-eas-cache",
        action="store_true",
        help="For aggregate-status only; report Han stage without validating EAS caches.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    root = args.repo_root.resolve()
    campaign = (args.campaign_dir or root / DEFAULT_CAMPAIGN_RELATIVE_PATH).resolve()
    integration = (
        args.integration_manifest or root / DEFAULT_INTEGRATION_RELATIVE_PATH
    ).resolve()
    if args.action == "quarantine-plan":
        payload = build_quarantine_plan(
            repo_root=root,
            campaign_dir=campaign,
            integration_manifest_path=integration,
            expected_stale_count=args.expected_stale_count,
        )
        output = args.output_dir or integration.parent
        write_quarantine_plan(output, payload)
        return 0
    if args.slim_bin is None:
        raise ValueError(f"{args.action} requires --slim-bin")
    if args.action == "validate":
        load_selected_production_bundle(
            integration, repo_root=root, slim_path=args.slim_bin
        )
        return 0
    if args.action == "simulate-eas":
        bundle = load_selected_production_bundle(
            integration, repo_root=root, slim_path=args.slim_bin
        )
        units = _eas_production_units(
            pd.read_csv(bundle.execution_units_path, sep="\t")
        )
        status = simulate_selected_units(
            units,
            repo_root=root,
            integration_manifest_path=integration,
            slim_path=args.slim_bin,
            workers=args.workers,
        )
        destination = campaign / "results/last_selected_production_invocation.tsv"
        destination.parent.mkdir(parents=True, exist_ok=True)
        temporary = destination.with_name(f".{destination.name}.{os.getpid()}.tmp")
        status.to_csv(temporary, sep="\t", index=False, lineterminator="\n")
        os.replace(temporary, destination)
        print(
            json.dumps(
                {
                    "status": (
                        "failed" if status["status"].eq("failed").any() else "complete"
                    ),
                    "eas_selected_units": len(status),
                    "fixed_max_draws": DEFAULT_EAS_MAX_DRAWS,
                    "fixed_draw_timeout_seconds": DEFAULT_DRAW_TIMEOUT_SECONDS,
                    "workers": int(args.workers),
                    "invocation_table": destination.relative_to(root).as_posix(),
                    "invocation_table_sha256": sha256_file(destination),
                },
                indent=2,
                sort_keys=True,
            )
        )
        return int(status["status"].eq("failed").any())
    if args.action in {"aggregate-status", "aggregate-validate"}:
        require_final = args.action == "aggregate-validate"
        payload = build_aggregate_readiness(
            repo_root=root,
            eas_integration_manifest_path=integration,
            slim_path=args.slim_bin,
            han_v3_manifest_path=(
                args.han_v3_manifest or root / DEFAULT_HAN_V3_MANIFEST_RELATIVE_PATH
            ),
            han_v3_pilot_audit_path=(
                args.han_v3_pilot_audit
                or root / DEFAULT_HAN_V3_PILOT_AUDIT_RELATIVE_PATH
            ),
            han_v3_final_audit_path=(
                args.han_v3_final_audit
                or root / DEFAULT_HAN_V3_FINAL_AUDIT_RELATIVE_PATH
            ),
            han_v3_authorization_path=(
                args.han_v3_authorization
                or root / DEFAULT_HAN_V3_AUTHORIZATION_RELATIVE_PATH
            ),
            han_v3_recovery_authorization_path=(
                args.han_v3_recovery_authorization
                or root / DEFAULT_HAN_V3_RECOVERY_AUTHORIZATION_RELATIVE_PATH
            ),
            han_v3_recovery_audit_path=(
                args.han_v3_recovery_audit
                or root / DEFAULT_HAN_V3_RECOVERY_AUDIT_RELATIVE_PATH
            ),
            require_final=require_final,
            validate_eas_cache=(require_final or not bool(args.status_skip_eas_cache)),
        )
        if require_final:
            readiness = (
                args.aggregate_readiness
                or root / DEFAULT_AGGREGATE_READINESS_RELATIVE_PATH
            )
            write_aggregate_readiness(readiness, payload)
        print(json.dumps(payload, indent=2, sort_keys=True))
        return 0
    if args.action in {"aggregate", "analyze", "report"}:
        readiness = (
            args.aggregate_readiness or root / DEFAULT_AGGREGATE_READINESS_RELATIVE_PATH
        )
        if args.action == "analyze":
            payload = analyze_with_readiness(
                repo_root=root,
                slim_path=args.slim_bin,
                readiness_path=readiness,
                empirical_path=args.empirical_tsv,
            )
            print(json.dumps(payload, indent=2, sort_keys=True))
            return 0
        if args.action == "report":
            payload = report_with_readiness(
                repo_root=root,
                slim_path=args.slim_bin,
                readiness_path=readiness,
            )
            print(json.dumps(payload, indent=2, sort_keys=True))
            return 0
        outputs = aggregate_with_readiness(
            repo_root=root,
            slim_path=args.slim_bin,
            readiness_path=readiness,
        )
        print(
            json.dumps(
                {
                    "status": "aggregated",
                    "readiness": Path(readiness).resolve().relative_to(root).as_posix(),
                    "outputs": [path.relative_to(root).as_posix() for path in outputs],
                    "eas_selected_validated": 60,
                    "han_v3_selected_validated": 60,
                    "neutral_cache_dispatch_unchanged": True,
                },
                indent=2,
                sort_keys=True,
            )
        )
        return 0
    payload = build_integration_manifest(
        repo_root=root,
        campaign_dir=campaign,
        slim_path=args.slim_bin,
        eas_adapter_authorization_path=(
            args.eas_adapter_authorization
            or root / DEFAULT_EAS_ADAPTER_AUTHORIZATION_RELATIVE_PATH
        ),
        eas_frozen_path=args.eas_frozen or root / DEFAULT_EAS_FROZEN_RELATIVE_PATH,
        integration_manifest_path=integration,
    )
    write_integration_manifest(integration, payload)
    return 0


__all__ = [
    "AGGREGATE_READINESS_SCHEMA",
    "DEFAULT_AGGREGATE_READINESS_RELATIVE_PATH",
    "DEFAULT_DRAW_TIMEOUT_SECONDS",
    "DEFAULT_EAS_MAX_DRAWS",
    "DEFAULT_EXPECTED_STALE_SELECTED_DIRS",
    "DEFAULT_HAN_V3_AUTHORIZATION_RELATIVE_PATH",
    "DEFAULT_HAN_V3_FINAL_AUDIT_RELATIVE_PATH",
    "DEFAULT_HAN_V3_MANIFEST_RELATIVE_PATH",
    "DEFAULT_HAN_V3_PILOT_AUDIT_RELATIVE_PATH",
    "DEFAULT_HAN_V3_RECOVERY_AUDIT_RELATIVE_PATH",
    "DEFAULT_HAN_V3_RECOVERY_AUTHORIZATION_RELATIVE_PATH",
    "DEFAULT_INTEGRATION_RELATIVE_PATH",
    "INTEGRATION_IMPLEMENTATION_PATHS",
    "INTEGRATION_MANIFEST_SCHEMA",
    "MODULE_SOURCE_PATH",
    "QUARANTINE_PLAN_SCHEMA",
    "SCHEMA_VERSION",
    "SelectedProductionBundle",
    "aggregate_with_readiness",
    "analyze_with_readiness",
    "build_aggregate_readiness",
    "build_integration_manifest",
    "build_quarantine_plan",
    "load_aggregate_readiness",
    "load_selected_production_bundle",
    "main",
    "patched_selected_executor",
    "report_with_readiness",
    "simulate_selected_units",
    "validate_eas_cached_units_read_only",
    "write_aggregate_readiness",
    "write_integration_manifest",
    "write_quarantine_plan",
]
