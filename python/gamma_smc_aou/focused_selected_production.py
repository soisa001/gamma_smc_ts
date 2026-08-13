"""Fail-closed production adapter for the focused selected simulations.

This module is deliberately additive.  It does not change the frozen focused
campaign, neutral executor, or either calibration implementation.  Instead it
validates the two independent selected-model authorizations, writes one
immutable integration manifest, and patches the selected-object constructor
only inside a selected worker process.

The EAS adapter retains the calibrated Q=5 de-novo origin, continuous
selection, and both terminal population-frequency conditions inside SLiM.  The
existing 500-diploid external AF gate and exact-k 100-diploid panel are still
applied by the unchanged executor.  Han uses the fixed-v2 authorization and its
loader-compatible frozen endpoint artifact.
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
from . import focused_han_fixed_v2_calibration as han_fixed_v2
from . import focused_selection_campaign as campaign_plan
from . import focused_selection_simulation as simulation
from .eas_sweep_models import EAS_POPULATION

SCHEMA_VERSION = "gamma-smc.focused-selected-production-adapter/v1"
INTEGRATION_MANIFEST_SCHEMA = "gamma-smc.focused-selected-production-integration/v1"
QUARANTINE_PLAN_SCHEMA = "gamma-smc.focused-selected-quarantine-plan/v1"

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
DEFAULT_HAN_AUXILIARY_RELATIVE_PATH = (
    "focused_selection_EAS_sim/calibration/han_fixed_cessation_v2/"
    "han_fixed_v2_authorization.json"
)
DEFAULT_HAN_FROZEN_RELATIVE_PATH = han_fixed_v2.DEFAULT_FROZEN_RELATIVE_PATH
DEFAULT_HAN_MANIFEST_RELATIVE_PATH = (
    "focused_selection_EAS_sim/calibration/han_fixed_cessation_v2/"
    "han_fixed_v2_manifest.json"
)

DEFAULT_EAS_MAX_DRAWS = 68
DEFAULT_HAN_MAX_DRAWS = 100
DEFAULT_DRAW_TIMEOUT_SECONDS = 30.0 * 60.0
DEFAULT_EXPECTED_STALE_SELECTED_DIRS = 24

INTEGRATION_IMPLEMENTATION_PATHS = (
    ".gitattributes",
    "pyproject.toml",
    "uv.lock",
    simulation.SIMULATION_SOURCE_PATH,
    "python/gamma_smc_aou/eas_sweep_models.py",
    "python/gamma_smc_aou/eas_sweep_analysis.py",
    "python/gamma_smc_aou/eas_sweep_study.py",
    "python/gamma_smc_aou/focused_selection_campaign.py",
    eas_calibration.MODULE_SOURCE_PATH,
    eas_audit.MODULE_SOURCE_PATH,
    "python/gamma_smc_aou/focused_han_fixed_v2_calibration.py",
    MODULE_SOURCE_PATH,
    WRAPPER_SOURCE_PATH,
    DEFAULT_INTEGRATION_RELATIVE_PATH,
)

_PATCH_LOCK = threading.RLock()
_WORKER_BUNDLE: SelectedProductionBundle | None = None
_WORKER_PATCH_CONTEXT: Any = None


@dataclass(frozen=True)
class SelectedProductionBundle:
    """Fully validated inputs needed by one selected worker process."""

    repo_root: Path
    campaign_dir: Path
    integration_manifest_path: Path
    slim_path: Path
    execution_units_path: Path
    manifest: Mapping[str, Any]
    eas_authorization: Mapping[str, Any]
    eas_frozen: Mapping[str, Any]
    han_auxiliary: Mapping[str, Any]
    han_binding: Mapping[str, Any]
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
    """Validate both layers of the EAS adapter authorization."""

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

    output, _, _, plan_manifest = eas_audit._load_bound_inputs(
        repo_root=repo_root,
        campaign_dir=campaign_dir,
        slim_path=slim_path,
    )
    if authorization.get("legacy_source_sha256") != plan_manifest.get(
        "source_sha256"
    ) or authorization.get("plan_manifest_contract_sha256") != plan_manifest.get(
        "contract_sha256"
    ):
        raise ValueError("EAS adapter plan/source authorization differs")
    expected_phases = set(eas_calibration.PHASE_ORDER)
    phase_records = authorization.get("phase_audits")
    if not isinstance(phase_records, Mapping) or set(phase_records) != expected_phases:
        raise ValueError("EAS adapter authorization lacks exact phase coverage")
    adapter_sha256 = adapter_record["sha256"]
    for phase in eas_calibration.PHASE_ORDER:
        sidecar_path = output / f"{phase}_adapter_audit.json"
        sidecar = eas_audit._load_phase_sidecar(
            phase,
            output=output,
            repo_root=repo_root,
            adapter_sha256=adapter_sha256,
            manifest=plan_manifest,
        )
        _exact_file_record(
            phase_records[phase], sidecar_path, repo_root, label=f"EAS {phase} audit"
        )
        if phase_records[phase].get("payload_sha256") != sidecar.get("payload_sha256"):
            raise ValueError(f"EAS {phase} audit payload binding differs")
    return authorization, frozen


def _load_han_authorization_bundle(
    *,
    repo_root: Path,
    execution_units_path: Path,
    auxiliary_path: Path,
    frozen_path: Path,
    manifest_path: Path,
    slim_scaling_factor: float,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Validate fixed-v2 auxiliary and its loader-compatible frozen artifact."""

    auxiliary = han_fixed_v2.load_frozen_fixed_v2(
        auxiliary_path,
        repo_root=repo_root,
        manifest_path=manifest_path,
    )
    _exact_file_record(
        {
            "path": str(auxiliary.get("loader_compatible_frozen_path", "")),
            "sha256": str(auxiliary.get("loader_compatible_frozen_sha256", "")),
            "size_bytes": frozen_path.stat().st_size if frozen_path.is_file() else -1,
        },
        frozen_path,
        repo_root,
        label="Han loader-compatible frozen authorization",
    )
    binding = simulation.load_frozen_han_selection_calibration(
        frozen_path,
        repo_root=repo_root,
        execution_units_path=execution_units_path,
        slim_scaling_factor=slim_scaling_factor,
        candidate_pool_diploids=simulation.DEFAULT_SELECTED_POOL_DIPLOIDS,
    )
    if _canonical_sha256(binding) != auxiliary.get("production_loader_binding_sha256"):
        raise ValueError("Han auxiliary and production loader binding differ")
    return auxiliary, binding


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
    han_auxiliary: Mapping[str, Any],
    han_binding: Mapping[str, Any],
) -> dict[str, Any]:
    return {
        "eas": {
            "adapter_payload_sha256": eas_authorization["payload_sha256"],
            "frozen_payload_sha256": eas_frozen["payload_sha256"],
            "cell_authorization_count": len(eas_frozen["cell_authorizations"]),
            "terminal_population_af_conditioning_in_slim": True,
            "continuous_selection_to_present": True,
        },
        "han": {
            "fixed_v2_manifest_sha256": han_auxiliary["manifest_sha256"],
            "loader_binding_sha256": _canonical_sha256(han_binding),
            "cell_authorization_count": len(han_auxiliary["cells"]),
            "post_pulse_recipient_nonloss_conditioning": True,
            "fixed_selection_cessation": True,
        },
    }


def build_integration_manifest(
    *,
    repo_root: str | Path,
    campaign_dir: str | Path,
    slim_path: str | Path,
    eas_adapter_authorization_path: str | Path,
    eas_frozen_path: str | Path,
    han_auxiliary_path: str | Path,
    han_frozen_path: str | Path,
    han_manifest_path: str | Path,
    integration_manifest_path: str | Path = DEFAULT_INTEGRATION_RELATIVE_PATH,
) -> dict[str, Any]:
    """Build, but do not publish, the selected integration authorization."""

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
    han_aux = _inside_repo(han_auxiliary_path, root, label="Han fixed-v2 authorization")
    han_frozen = _inside_repo(han_frozen_path, root, label="Han frozen authorization")
    han_manifest = _inside_repo(han_manifest_path, root, label="Han fixed-v2 manifest")
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
    han_auxiliary, han_binding = _load_han_authorization_bundle(
        repo_root=root,
        execution_units_path=units_path,
        auxiliary_path=han_aux,
        frozen_path=han_frozen,
        manifest_path=han_manifest,
        slim_scaling_factor=simulation.DEFAULT_SLIM_SCALING_FACTOR,
    )
    payload: dict[str, Any] = {
        "schema": INTEGRATION_MANIFEST_SCHEMA,
        "status": "ready",
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
            "han_max_draws": DEFAULT_HAN_MAX_DRAWS,
            "minimum_per_draw_timeout_seconds": DEFAULT_DRAW_TIMEOUT_SECONDS,
        },
        "authorization_artifacts": {
            "eas_adapter": _file_record(eas_adapter_path, root),
            "eas_frozen": _file_record(eas_frozen, root),
            "han_fixed_v2": _file_record(han_aux, root),
            "han_loader_compatible_frozen": _file_record(han_frozen, root),
            "han_fixed_v2_manifest": _file_record(han_manifest, root),
        },
        "authorization_bindings": _binding_summary(
            eas_authorization=eas_authorization,
            eas_frozen=eas_payload,
            han_auxiliary=han_auxiliary,
            han_binding=han_binding,
        ),
        "implementation": {
            "mode": "selected_worker_process_local_patch",
            "neutral_executor_modified": False,
            "campaign_plan_modified": False,
            "selected_objects_patch_only": True,
            "contract_paths": list(INTEGRATION_IMPLEMENTATION_PATHS),
            "source_sha256": _implementation_records(root),
        },
        "scientific_contract": {
            "eas": (
                "Q5 de-novo origin; continuous selection to present; terminal "
                "population AF bounds retained in SLiM; external 500-diploid AF "
                "gate and exact-k N100 panel retained"
            ),
            "han": (
                "Q5 archaic-specific origin; post-pulse recipient nonloss; "
                "fixed-v2 calibrated cessation; external 500-diploid AF gate "
                "and exact-k N100 panel retained"
            ),
            "gamma_smc_statistics_used_for_calibration": False,
        },
        "quarantine": {
            "planning_only": True,
            "moves_or_deletions_authorized": False,
            "expected_preintegration_stale_selected_directories": (
                DEFAULT_EXPECTED_STALE_SELECTED_DIRS
            ),
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
        expected_status="ready",
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
    """Strictly reload the integration manifest and all four authorizations."""

    root = Path(repo_root).resolve()
    manifest_path = _inside_repo(path, root, label="integration manifest")
    if manifest_path.relative_to(root).as_posix() != DEFAULT_INTEGRATION_RELATIVE_PATH:
        raise ValueError("selected integration manifest path differs")
    payload = _read_json(manifest_path, label="selected integration manifest")
    _validate_hashed_payload(
        payload,
        schema=INTEGRATION_MANIFEST_SCHEMA,
        status_key="status",
        expected_status="ready",
        label="selected integration manifest",
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
        or int(runtime.get("han_max_draws", -1)) != DEFAULT_HAN_MAX_DRAWS
        or float(runtime.get("minimum_per_draw_timeout_seconds", math.nan))
        != DEFAULT_DRAW_TIMEOUT_SECONDS
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
        "han_fixed_v2",
        "han_loader_compatible_frozen",
        "han_fixed_v2_manifest",
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
    han_auxiliary, han_binding = _load_han_authorization_bundle(
        repo_root=root,
        execution_units_path=units_path,
        auxiliary_path=artifact_paths["han_fixed_v2"],
        frozen_path=artifact_paths["han_loader_compatible_frozen"],
        manifest_path=artifact_paths["han_fixed_v2_manifest"],
        slim_scaling_factor=simulation.DEFAULT_SLIM_SCALING_FACTOR,
    )
    if payload.get("authorization_bindings") != _binding_summary(
        eas_authorization=eas_authorization,
        eas_frozen=eas_frozen,
        han_auxiliary=han_auxiliary,
        han_binding=han_binding,
    ):
        raise ValueError("selected integration authorization bindings differ")

    implementation = payload.get("implementation")
    if (
        not isinstance(implementation, Mapping)
        or implementation.get("mode") != "selected_worker_process_local_patch"
        or implementation.get("neutral_executor_modified") is not False
        or implementation.get("campaign_plan_modified") is not False
        or implementation.get("selected_objects_patch_only") is not True
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
        han_auxiliary=han_auxiliary,
        han_binding=han_binding,
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
                if han_selection_end_generations_ago is not None:
                    raise ValueError("EAS selected adapter received a Han endpoint")
                return _authorized_eas_selected_objects(
                    record,
                    bundle=bundle,
                    slim_scaling_factor=slim_scaling_factor,
                    selected_pool_diploids=selected_pool_diploids,
                )
            if record["demography_id"] == "ancient_eurasia_han_introgression":
                if (
                    han_selection_end_generations_ago is not None
                    or han_selection_calibration is None
                    or _canonical_sha256(han_selection_calibration)
                    != _canonical_sha256(bundle.han_binding)
                ):
                    raise ValueError(
                        "Han selected worker binding differs from fixed-v2"
                    )
                model, sweep, samples, provenance = original_objects(
                    record,
                    repo_root,
                    slim_scaling_factor,
                    selected_pool_diploids,
                    han_selection_calibration=bundle.han_binding,
                )
                provenance = dict(provenance)
                provenance["selected_production_authorization"] = {
                    "schema": SCHEMA_VERSION,
                    "demography_id": "ancient_eurasia_han_introgression",
                    "integration_manifest_path": (
                        bundle.integration_manifest_path.relative_to(
                            bundle.repo_root
                        ).as_posix()
                    ),
                    "integration_payload_sha256": bundle.manifest["payload_sha256"],
                    "han_fixed_v2_manifest_sha256": bundle.han_auxiliary[
                        "manifest_sha256"
                    ],
                    "han_loader_binding_sha256": _canonical_sha256(bundle.han_binding),
                    "post_pulse_recipient_nonloss_conditioning": True,
                    "external_candidate_pool_and_exact_panel_retained": True,
                }
                return model, sweep, samples, provenance
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
    if requested.empty:
        raise ValueError("selected production request contains no units")
    if set(requested["simulation_class"].astype(str)) != {"selected"}:
        raise ValueError("selected production adapter refuses neutral units")
    if requested["unit_id"].astype(str).duplicated().any():
        raise ValueError("selected production request contains duplicate unit IDs")
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


def _worker_configuration(
    bundle: SelectedProductionBundle,
) -> dict[str, str]:
    return {
        "repo_root": str(bundle.repo_root),
        "integration_manifest_path": str(bundle.integration_manifest_path),
        "slim_path": str(bundle.slim_path),
    }


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
        artifacts = simulation.simulate_unit(
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
            han_selection_calibration=_WORKER_BUNDLE.han_binding,
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
    han_max_draws: int = DEFAULT_HAN_MAX_DRAWS,
    draw_timeout_seconds: float = DEFAULT_DRAW_TIMEOUT_SECONDS,
) -> pd.DataFrame:
    """Run only selected units through isolated, process-local adapter workers."""

    if int(workers) < 1 or int(workers) > 24:
        raise ValueError("selected production workers must be between 1 and 24")
    if int(eas_max_draws) < 1 or int(han_max_draws) < 1:
        raise ValueError("selected production retry budgets must be positive")
    if float(draw_timeout_seconds) < DEFAULT_DRAW_TIMEOUT_SECONDS:
        raise ValueError("selected production per-draw timeout must be at least 30 min")
    bundle = load_selected_production_bundle(
        integration_manifest_path,
        repo_root=repo_root,
        slim_path=slim_path,
    )
    if int(eas_max_draws) != int(bundle.manifest["runtime"]["eas_max_draws"]) or int(
        han_max_draws
    ) != int(bundle.manifest["runtime"]["han_max_draws"]):
        raise ValueError("selected production retry budgets differ from integration")
    canonical = pd.read_csv(bundle.execution_units_path, sep="\t")
    selected = _validate_requested_units(units.copy(), canonical)
    tasks = []
    for row in selected.to_dict(orient="records"):
        is_eas = str(row["demography_id"]) == "eas_phlash_median"
        max_draws = int(eas_max_draws if is_eas else han_max_draws)
        tasks.append(
            {
                "unit": row,
                "selected_max_draws": max_draws,
                "draw_timeout_seconds": float(draw_timeout_seconds),
                "cumulative_timeout_seconds": float(draw_timeout_seconds) * max_draws,
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
    """Inventory stale planned selected directories without moving anything."""

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
    selected = units[units["simulation_class"].astype(str).eq("selected")]
    work = _inside_repo(campaign / "work", root, label="campaign work")
    destination_root = _inside_repo(
        campaign / "quarantine" / "pre_focused_selected_production",
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
                    "selected_directory_predates_frozen_production_adapter_binding"
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
            "stale planned selected-directory count differs: "
            f"observed {len(rows)}, expected {int(expected_stale_count)}"
        )
    if any(row["simulation_class"] != "selected" for row in rows):
        raise RuntimeError("neutral unit entered selected quarantine plan")
    payload: dict[str, Any] = {
        "schema": QUARANTINE_PLAN_SCHEMA,
        "status": "planned_not_executed",
        "execution_supported_by_this_module": False,
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
        choices=("build-manifest", "validate", "quarantine-plan", "simulate"),
    )
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--campaign-dir", type=Path)
    parser.add_argument("--slim-bin", type=Path)
    parser.add_argument("--integration-manifest", type=Path)
    parser.add_argument("--eas-adapter-authorization", type=Path)
    parser.add_argument("--eas-frozen", type=Path)
    parser.add_argument("--han-fixed-v2-authorization", type=Path)
    parser.add_argument("--han-frozen", type=Path)
    parser.add_argument("--han-fixed-v2-manifest", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument(
        "--expected-stale-count",
        type=int,
        default=DEFAULT_EXPECTED_STALE_SELECTED_DIRS,
    )
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument(
        "--demography",
        action="append",
        choices=("eas_phlash_median", "ancient_eurasia_han_introgression"),
    )
    parser.add_argument("--af", action="append", type=float)
    parser.add_argument("--selection-coefficient", action="append", type=float)
    parser.add_argument("--unit", action="append")
    parser.add_argument("--eas-max-draws", type=int, default=DEFAULT_EAS_MAX_DRAWS)
    parser.add_argument("--han-max-draws", type=int, default=DEFAULT_HAN_MAX_DRAWS)
    parser.add_argument(
        "--draw-timeout-minutes",
        type=float,
        default=DEFAULT_DRAW_TIMEOUT_SECONDS / 60.0,
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
    if args.action == "simulate":
        bundle = load_selected_production_bundle(
            integration, repo_root=root, slim_path=args.slim_bin
        )
        units = pd.read_csv(bundle.execution_units_path, sep="\t")
        units = units[units["simulation_class"].astype(str).eq("selected")].copy()
        if args.demography:
            units = units[units["demography_id"].astype(str).isin(args.demography)]
        if args.af:
            units = units[
                np.logical_or.reduce(
                    [
                        np.isclose(
                            units["target_allele_frequency"].astype(float), value
                        )
                        for value in args.af
                    ]
                )
            ]
        if args.selection_coefficient:
            units = units[
                np.logical_or.reduce(
                    [
                        np.isclose(units["selection_coefficient"].astype(float), value)
                        for value in args.selection_coefficient
                    ]
                )
            ]
        if args.unit:
            units = units[
                units["unit_id"]
                .astype(str)
                .map(lambda value: any(token in value for token in args.unit))
            ]
        status = simulate_selected_units(
            units,
            repo_root=root,
            integration_manifest_path=integration,
            slim_path=args.slim_bin,
            workers=args.workers,
            eas_max_draws=args.eas_max_draws,
            han_max_draws=args.han_max_draws,
            draw_timeout_seconds=args.draw_timeout_minutes * 60.0,
        )
        destination = campaign / "results/last_selected_production_invocation.tsv"
        destination.parent.mkdir(parents=True, exist_ok=True)
        temporary = destination.with_name(f".{destination.name}.{os.getpid()}.tmp")
        status.to_csv(temporary, sep="\t", index=False, lineterminator="\n")
        os.replace(temporary, destination)
        return int(status["status"].eq("failed").any())
    payload = build_integration_manifest(
        repo_root=root,
        campaign_dir=campaign,
        slim_path=args.slim_bin,
        eas_adapter_authorization_path=(
            args.eas_adapter_authorization
            or root / DEFAULT_EAS_ADAPTER_AUTHORIZATION_RELATIVE_PATH
        ),
        eas_frozen_path=args.eas_frozen or root / DEFAULT_EAS_FROZEN_RELATIVE_PATH,
        han_auxiliary_path=(
            args.han_fixed_v2_authorization
            or root / DEFAULT_HAN_AUXILIARY_RELATIVE_PATH
        ),
        han_frozen_path=args.han_frozen or root / DEFAULT_HAN_FROZEN_RELATIVE_PATH,
        han_manifest_path=(
            args.han_fixed_v2_manifest or root / DEFAULT_HAN_MANIFEST_RELATIVE_PATH
        ),
        integration_manifest_path=integration,
    )
    write_integration_manifest(integration, payload)
    return 0


__all__ = [
    "DEFAULT_DRAW_TIMEOUT_SECONDS",
    "DEFAULT_EAS_MAX_DRAWS",
    "DEFAULT_EXPECTED_STALE_SELECTED_DIRS",
    "DEFAULT_HAN_MAX_DRAWS",
    "DEFAULT_INTEGRATION_RELATIVE_PATH",
    "INTEGRATION_IMPLEMENTATION_PATHS",
    "INTEGRATION_MANIFEST_SCHEMA",
    "MODULE_SOURCE_PATH",
    "QUARANTINE_PLAN_SCHEMA",
    "SCHEMA_VERSION",
    "SelectedProductionBundle",
    "build_integration_manifest",
    "build_quarantine_plan",
    "load_selected_production_bundle",
    "main",
    "patched_selected_executor",
    "simulate_selected_units",
    "write_integration_manifest",
    "write_quarantine_plan",
]
