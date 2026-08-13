"""Provenance-preserving audit/freeze adapter for the bound EAS calibration.

The original calibration plan and every draw completion bind the byte identity
of :mod:`focused_eas_selected_calibration`.  This adapter therefore leaves that
module untouched.  It repairs only the lossy pandas TSV read of three nullable
text fields before delegating every scientific and artifact check to the bound
implementation.
"""

from __future__ import annotations

import argparse
import json
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import pandas as pd

from . import focused_eas_selected_calibration as legacy

SCHEMA_VERSION = "gamma-smc.eas-selected-audit-adapter/v1"
AUTHORIZATION_SCHEMA_VERSION = "gamma-smc.eas-selected-audit-adapter-authorization/v1"
MODULE_SOURCE_PATH = "python/gamma_smc_aou/focused_eas_selected_audit_adapter.py"
TEXT_COLUMNS_WITH_STRUCTURAL_EMPTY_STRING = (
    "panel_tree_sha256",
    "panel_manifest_sha256",
    "error",
)


def _adapter_sha256(repo_root: Path) -> str:
    path = repo_root / MODULE_SOURCE_PATH
    if not path.is_file():
        raise ValueError(f"EAS audit adapter source is absent: {path}")
    return legacy.sha256_file(path)


def _normalized_ledger(path: Path) -> pd.DataFrame:
    """Read a legacy ledger while restoring only its structural text blanks."""

    ledger = pd.read_csv(path, sep="\t")
    missing = [
        column
        for column in TEXT_COLUMNS_WITH_STRUCTURAL_EMPTY_STRING
        if column not in ledger
    ]
    if missing:
        raise ValueError(
            "EAS calibration ledger lacks adapter text columns: " + ", ".join(missing)
        )
    for column in TEXT_COLUMNS_WITH_STRUCTURAL_EMPTY_STRING:
        ledger[column] = ledger[column].fillna("")
    return ledger


def _canonical_file_record(path: Path, repo_root: Path) -> dict[str, Any]:
    resolved = path.resolve()
    return {
        "path": resolved.relative_to(repo_root).as_posix(),
        "sha256": legacy.sha256_file(resolved),
        "size_bytes": resolved.stat().st_size,
    }


def _completion_mapping_record(audit: pd.DataFrame, repo_root: Path) -> dict[str, Any]:
    columns = (
        "calibration_id",
        "completion_path",
        "completion_sha256",
        "trajectory_sha256",
        "panel_tree_sha256",
        "panel_manifest_sha256",
    )
    canonical = audit.loc[:, columns].copy()
    for column in ("panel_tree_sha256", "panel_manifest_sha256"):
        canonical[column] = canonical[column].fillna("")
    mapping = canonical.sort_values("calibration_id").to_dict(orient="records")
    for row in mapping:
        completion = legacy._inside_repo(
            repo_root / str(row["completion_path"]),
            repo_root,
            label="EAS completion",
        )
        if not completion.is_file() or legacy.sha256_file(completion) != str(
            row["completion_sha256"]
        ):
            raise ValueError("EAS completion mapping changed during adapter audit")
    return {
        "rows": len(mapping),
        "canonical_sha256": legacy._canonical_sha256(mapping),
        "entries": mapping,
    }


def _load_bound_inputs(
    *, repo_root: Path, campaign_dir: Path, slim_path: Path
) -> tuple[Path, Path, dict[str, pd.DataFrame], dict[str, Any]]:
    output = legacy._inside_repo(
        campaign_dir / "calibration" / "eas_selected",
        repo_root,
        label="EAS calibration output",
    )
    work = legacy._inside_repo(
        campaign_dir / "work" / "eas_selected_calibration",
        repo_root,
        label="EAS calibration work root",
    )
    plans, manifest = legacy.read_plan_bundle(output)
    legacy._manifest_current_binding(manifest, repo_root, slim_path)
    if manifest.get("source_sha256") != legacy._source_hashes(repo_root):
        raise ValueError("EAS calibration source binding differs from signed plan")
    for phase in legacy.PHASE_ORDER:
        legacy._validate_plan_manifest_binding(plans[phase], phase, manifest)
    return output, work, plans, manifest


def validate_phase_with_adapter(
    phase: str,
    *,
    repo_root: str | Path,
    campaign_dir: str | Path,
    slim_path: str | Path,
) -> dict[str, Any]:
    """Audit one phase and atomically write canonical products plus a sidecar."""

    if phase not in legacy.PHASE_ORDER:
        raise ValueError(f"unknown EAS calibration phase: {phase}")
    root = Path(repo_root).resolve()
    campaign = legacy._inside_repo(
        campaign_dir, root, label="focused campaign directory"
    )
    slim = Path(slim_path).resolve()
    output, work, plans, manifest = _load_bound_inputs(
        repo_root=root, campaign_dir=campaign, slim_path=slim
    )
    ledger_path = output / f"{phase}_ledger.tsv"
    input_ledger = _canonical_file_record(ledger_path, root)
    normalized = _normalized_ledger(ledger_path)
    rebuilt, completion_audit = legacy.audit_phase_completions(
        plans[phase],
        normalized,
        repo_root=root,
        work_root=work,
        slim_path=slim,
        manifest=manifest,
    )
    summary = legacy.summarize_phase(plans[phase], rebuilt, phase)
    if not summary["passed"].map(legacy._coerce_boolean).all():
        raise ValueError(f"{phase} does not pass every EAS calibration cell")

    summary_path = output / f"{phase}_summary.tsv"
    checksums_path = output / f"{phase}_completion_checksums.tsv"
    legacy._atomic_frame(ledger_path, rebuilt)
    legacy._atomic_frame(summary_path, summary)
    legacy._atomic_frame(checksums_path, completion_audit)

    mapping = _completion_mapping_record(completion_audit, root)
    sidecar_path = output / f"{phase}_adapter_audit.json"
    payload: dict[str, Any] = {
        "schema": SCHEMA_VERSION,
        "status": "validated",
        "phase": phase,
        "adapter": {
            "path": MODULE_SOURCE_PATH,
            "sha256": _adapter_sha256(root),
        },
        "legacy_source_sha256": dict(manifest["source_sha256"]),
        "plan_manifest": {
            **_canonical_file_record(
                output / "eas_selected_calibration_plan.json", root
            ),
            "contract_sha256": str(manifest["contract_sha256"]),
        },
        "plan": _canonical_file_record(output / f"{phase}_plan.tsv", root),
        "slim": _canonical_file_record(slim, root),
        "eas_resource": _canonical_file_record(root / legacy.EAS_RESOURCE_PATH, root),
        "normalization": {
            "policy": "fill_default_na_only_for_structural_text_fields",
            "columns": list(TEXT_COLUMNS_WITH_STRUCTURAL_EMPTY_STRING),
            "input_ledger": input_ledger,
        },
        "outputs": {
            "rebuilt_ledger": _canonical_file_record(ledger_path, root),
            "summary": _canonical_file_record(summary_path, root),
            "completion_checksums": _canonical_file_record(checksums_path, root),
        },
        "completion_mapping": mapping,
        "counts": {
            "planned": len(plans[phase]),
            "validated_completions": len(rebuilt),
            "completion_checksum_rows": len(completion_audit),
            "cells": len(summary),
            "passing_cells": int(summary["passed"].map(legacy._coerce_boolean).sum()),
        },
    }
    payload["payload_sha256"] = legacy._canonical_sha256(payload)
    legacy._atomic_json(sidecar_path, payload)
    return payload


def _load_phase_sidecar(
    phase: str,
    *,
    output: Path,
    repo_root: Path,
    adapter_sha256: str,
    manifest: Mapping[str, Any],
) -> dict[str, Any]:
    path = output / f"{phase}_adapter_audit.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    unhashed = dict(payload)
    digest = unhashed.pop("payload_sha256", None)
    if (
        payload.get("schema") != SCHEMA_VERSION
        or payload.get("status") != "validated"
        or payload.get("phase") != phase
        or digest != legacy._canonical_sha256(unhashed)
    ):
        raise ValueError(f"{phase} EAS adapter audit sidecar is incompatible")
    if payload.get("adapter") != {
        "path": MODULE_SOURCE_PATH,
        "sha256": adapter_sha256,
    }:
        raise ValueError(f"{phase} EAS adapter source binding differs")
    if payload.get("legacy_source_sha256") != manifest.get("source_sha256"):
        raise ValueError(f"{phase} EAS legacy source binding differs")
    expected_outputs = {
        "rebuilt_ledger": output / f"{phase}_ledger.tsv",
        "summary": output / f"{phase}_summary.tsv",
        "completion_checksums": output / f"{phase}_completion_checksums.tsv",
    }
    for key, path_value in expected_outputs.items():
        if payload.get("outputs", {}).get(key) != _canonical_file_record(
            path_value, repo_root
        ):
            raise ValueError(f"{phase} EAS adapter output binding differs: {key}")
    audit = pd.read_csv(expected_outputs["completion_checksums"], sep="\t")
    if _completion_mapping_record(audit, repo_root) != payload.get(
        "completion_mapping"
    ):
        raise ValueError(f"{phase} EAS completion mapping binding differs")
    return payload


def freeze_with_adapter(
    *,
    repo_root: str | Path,
    campaign_dir: str | Path,
    slim_path: str | Path,
) -> dict[str, Any]:
    """Require all adapter audits, then freeze and strictly reload authorization."""

    root = Path(repo_root).resolve()
    campaign = legacy._inside_repo(
        campaign_dir, root, label="focused campaign directory"
    )
    slim = Path(slim_path).resolve()
    output, work, plans, manifest = _load_bound_inputs(
        repo_root=root, campaign_dir=campaign, slim_path=slim
    )
    adapter_sha256 = _adapter_sha256(root)
    sidecars = {
        phase: _load_phase_sidecar(
            phase,
            output=output,
            repo_root=root,
            adapter_sha256=adapter_sha256,
            manifest=manifest,
        )
        for phase in legacy.PHASE_ORDER
    }
    ledgers = {
        phase: _normalized_ledger(output / f"{phase}_ledger.tsv")
        for phase in legacy.PHASE_ORDER
    }
    reservations = legacy.collect_seed_reservations(campaign)
    units_path = campaign / "execution_units.tsv"
    frozen = legacy.freeze_calibration(
        plans,
        ledgers,
        manifest,
        repo_root=root,
        output_dir=output,
        work_root=work,
        slim_path=slim,
        execution_units_path=units_path,
        seed_reservations=reservations,
    )
    if frozen.get("status") != "frozen":
        raise ValueError("legacy EAS calibration did not freeze")
    frozen_path = output / "eas_selected_calibration_frozen.json"
    loaded = legacy.load_frozen_eas_selected_calibration(
        frozen_path,
        repo_root=root,
        execution_units_path=units_path,
        slim_path=slim,
    )
    if loaded.get("payload_sha256") != frozen.get("payload_sha256"):
        raise ValueError("strict EAS frozen-loader result differs")

    # The legacy freezer rewrites the canonical phase products.  Require them
    # to remain byte-identical to the already authorized adapter outputs.
    sidecars = {
        phase: _load_phase_sidecar(
            phase,
            output=output,
            repo_root=root,
            adapter_sha256=adapter_sha256,
            manifest=manifest,
        )
        for phase in legacy.PHASE_ORDER
    }

    authorization_path = output / "eas_selected_calibration_frozen.adapter.json"
    payload: dict[str, Any] = {
        "schema": AUTHORIZATION_SCHEMA_VERSION,
        "status": "authorized",
        "adapter": {"path": MODULE_SOURCE_PATH, "sha256": adapter_sha256},
        "legacy_source_sha256": dict(manifest["source_sha256"]),
        "plan_manifest_contract_sha256": str(manifest["contract_sha256"]),
        "phase_audits": {
            phase: {
                **_canonical_file_record(output / f"{phase}_adapter_audit.json", root),
                "payload_sha256": str(sidecars[phase]["payload_sha256"]),
            }
            for phase in legacy.PHASE_ORDER
        },
        "frozen_authorization": {
            **_canonical_file_record(frozen_path, root),
            "payload_sha256": str(frozen["payload_sha256"]),
        },
    }
    payload["payload_sha256"] = legacy._canonical_sha256(payload)
    legacy._atomic_json(authorization_path, payload)
    return payload


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("validate", "freeze"))
    parser.add_argument("--phase", choices=legacy.PHASE_ORDER)
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--campaign-dir", type=Path)
    parser.add_argument("--slim-bin", type=Path, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    root = args.repo_root.resolve()
    campaign = args.campaign_dir or root / "focused_selection_EAS_sim"
    if args.action == "validate":
        if args.phase is None:
            raise ValueError("validate requires --phase")
        validate_phase_with_adapter(
            args.phase,
            repo_root=root,
            campaign_dir=campaign,
            slim_path=args.slim_bin,
        )
        return 0
    if args.phase is not None:
        raise ValueError("freeze does not accept --phase")
    freeze_with_adapter(
        repo_root=root,
        campaign_dir=campaign,
        slim_path=args.slim_bin,
    )
    return 0


__all__ = [
    "AUTHORIZATION_SCHEMA_VERSION",
    "MODULE_SOURCE_PATH",
    "SCHEMA_VERSION",
    "TEXT_COLUMNS_WITH_STRUCTURAL_EMPTY_STRING",
    "freeze_with_adapter",
    "main",
    "validate_phase_with_adapter",
]
