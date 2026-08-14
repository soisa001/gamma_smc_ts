"""Additive final-audit adapter for completed Han direct-v3 recovery outputs.

The frozen direct-v3 implementation correctly represents units without an
interruption journal with ``None`` values.  When units with and without such a
journal are combined in one pandas frame, pandas promotes those two structural
absences to floating-point ``NaN``.  The frozen canonical JSON encoder then
rejects the all-unit audit, as intended.

This adapter preserves the frozen simulator and dtype-recovery sources.  It
normalizes only the two interruption-journal reference fields, and only for a
row whose already-validated interruption count is zero.  Every other
non-finite value fails closed.  Finalization also refuses to run unless every
unit in the existing recovery authorization is already complete, so this
module cannot launch or retry a simulation.
"""

from __future__ import annotations

import argparse
import json
import math
import numbers
import os
import re
import shutil
from collections.abc import Iterator, Mapping, Sequence
from contextlib import contextmanager
from pathlib import Path
from time import time_ns
from typing import Any

import pandas as pd

from . import focused_han_direct_v3 as v3
from . import focused_han_direct_v3_recovery as recovery

AUDIT_NORMALIZATION_SCHEMA = (
    "gamma-smc.han-direct-production-v3-recovery-audit-normalization/v1"
)
MODULE_SOURCE_PATH = "python/gamma_smc_aou/focused_han_direct_v3_recovery_audit.py"
WRAPPER_SOURCE_PATH = "scripts/run_focused_han_direct_v3_recovery_audit.py"
ADAPTER_SOURCE_PATHS = (MODULE_SOURCE_PATH, WRAPPER_SOURCE_PATH)
CANONICAL_NORMALIZATION_AUDIT_NAME = "han_direct_v3_recovery_audit_normalization.json"
STRUCTURAL_ABSENCE_FIELDS = (
    "reconciliation_path",
    "reconciliation_sha256",
)
_SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")


def _canonical_sha256(value: Any) -> str:
    return v3._canonical_sha256(value)


def _atomic_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f"{path.name}.tmp.{os.getpid()}")
    try:
        temporary.write_text(
            json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _immutable_json(path: Path, payload: Any) -> None:
    if path.is_file():
        try:
            current = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            raise ValueError(
                f"Han v3 audit-adapter JSON is unreadable: {path}"
            ) from error
        if current != payload:
            raise ValueError(f"immutable Han v3 audit-adapter JSON differs: {path}")
        return
    _atomic_json(path, payload)


def _is_nan(value: Any) -> bool:
    return (
        isinstance(value, numbers.Real)
        and not isinstance(value, (bool, numbers.Integral))
        and math.isnan(float(value))
    )


def _is_nonfinite(value: Any) -> bool:
    return (
        isinstance(value, numbers.Real)
        and not isinstance(value, (bool, numbers.Integral))
        and not math.isfinite(float(value))
    )


def _normalize_reconciliation_absence(frame: pd.DataFrame) -> pd.DataFrame:
    """Normalize only validated, structural interruption-journal absences."""

    required = {
        "unit_id",
        "interrupted_nonterminal_dispositions",
        *STRUCTURAL_ABSENCE_FIELDS,
    }
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(
            "Han v3 recovery-audit frame lacks required fields: "
            + ", ".join(sorted(missing))
        )
    normalized = frame.copy()
    for field in STRUCTURAL_ABSENCE_FIELDS:
        normalized[field] = normalized[field].astype("object")

    for index, row in frame.iterrows():
        unit_id = str(row["unit_id"])
        raw_count = row["interrupted_nonterminal_dispositions"]
        if isinstance(raw_count, bool) or not isinstance(raw_count, numbers.Integral):
            raise ValueError(
                f"Han v3 recovery-audit interruption count differs: {unit_id}"
            )
        count = int(raw_count)
        if count < 0:
            raise ValueError(
                f"Han v3 recovery-audit interruption count differs: {unit_id}"
            )

        for column, value in row.items():
            if _is_nonfinite(value) and not (
                column in STRUCTURAL_ABSENCE_FIELDS and count == 0 and _is_nan(value)
            ):
                raise ValueError(
                    "Han v3 recovery-audit unexpected nonfinite value: "
                    f"{unit_id}.{column}"
                )

        path_value = row["reconciliation_path"]
        sha_value = row["reconciliation_sha256"]
        if count == 0:
            if not (
                (path_value is None or _is_nan(path_value))
                and (sha_value is None or _is_nan(sha_value))
            ):
                raise ValueError(
                    f"Han v3 recovery-audit structural absence differs: {unit_id}"
                )
            normalized.at[index, "reconciliation_path"] = None
            normalized.at[index, "reconciliation_sha256"] = None
        elif (
            not isinstance(path_value, str)
            or not path_value.endswith("/supervision_reconciliation.tsv")
            or not isinstance(sha_value, str)
            or _SHA256_PATTERN.fullmatch(sha_value) is None
        ):
            raise ValueError(
                f"Han v3 recovery-audit reconciliation reference differs: {unit_id}"
            )

    for record in normalized.to_dict(orient="records"):
        for column, value in record.items():
            if _is_nonfinite(value):
                raise ValueError(
                    "Han v3 recovery-audit normalization retained a nonfinite value: "
                    f"{record['unit_id']}.{column}"
                )
    return normalized


@contextmanager
def patched_recovery_audit_aggregation() -> Iterator[None]:
    """Patch only the in-memory aggregate returned to the frozen audit code."""

    frozen_aggregate = v3.aggregate_v3_cached_completions

    def normalized_aggregate(
        bundle: v3.DirectV3Bundle,
        stage: str = v3.STAGE_ALL,
    ) -> pd.DataFrame:
        return _normalize_reconciliation_absence(frozen_aggregate(bundle, stage))

    v3.aggregate_v3_cached_completions = normalized_aggregate
    try:
        yield
    finally:
        v3.aggregate_v3_cached_completions = frozen_aggregate


def _read_hashed_json(
    path: Path,
    *,
    schema: str,
    status: str,
    label: str,
) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"Han v3 {label} is unreadable") from error
    digest = payload.get("payload_sha256")
    unhashed = dict(payload)
    unhashed.pop("payload_sha256", None)
    if (
        payload.get("schema") != schema
        or payload.get("status") != status
        or not isinstance(digest, str)
        or digest != _canonical_sha256(unhashed)
    ):
        raise ValueError(f"Han v3 {label} is corrupt or incompatible")
    return payload


def _adapter_source_records(repo_root: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    for relative in ADAPTER_SOURCE_PATHS:
        path = repo_root / relative
        if not path.is_file():
            raise ValueError(f"Han v3 recovery-audit source is absent: {relative}")
        records[relative] = v3.sha256_file(path)
    return records


def _audit_adapter_binding(bundle: recovery.RecoveryBundle) -> dict[str, Any]:
    return {
        "sources": _adapter_source_records(bundle.repo_root),
        "audit_only_no_simulation_route": True,
        "normalized_fields": list(STRUCTURAL_ABSENCE_FIELDS),
        "normalization_rule": (
            "pandas NaN to JSON null only when "
            "interrupted_nonterminal_dispositions equals zero"
        ),
        "all_other_nonfinite_values_rejected": True,
        "scientific_rows_changed": False,
        "seed_retry_timeout_budget_or_reconciliation_semantics_changed": False,
    }


def _validate_recovery_completion_audit(
    bundle: recovery.RecoveryBundle,
    authorization: Mapping[str, Any],
    base_audit: Mapping[str, Any],
    *,
    require_adapter_binding: bool = True,
) -> dict[str, Any]:
    audit_path = bundle.output_dir / recovery.CANONICAL_AUDIT_NAME
    audit = _read_hashed_json(
        audit_path,
        schema=recovery.RECOVERY_COMPLETION_AUDIT_SCHEMA,
        status="passed",
        label="dtype-recovery completion audit",
    )
    authorized_ids = sorted(
        str(value) for value in authorization["recoverable_unit_ids"]
    )
    rows = audit.get("rows")
    row_ids = (
        sorted(str(row.get("unit_id", "")) for row in rows if isinstance(row, Mapping))
        if isinstance(rows, list)
        else []
    )
    if (
        row_ids != authorized_ids
        or sorted(
            str(value) for value in audit.get("post_recovery_recovered_unit_ids", ())
        )
        != authorized_ids
        or audit.get("post_recovery_incomplete_unit_ids") != []
        or int(audit.get("remaining_incomplete_units", -1)) != 0
        or audit.get("post_run_base_audit") != dict(base_audit)
        or audit.get("post_run_base_audit_payload_sha256")
        != base_audit["payload_sha256"]
        or audit.get("recovery_sources") != dict(bundle.recovery_sources)
        or audit.get("canonical_authorization_sha256")
        != v3.sha256_file(bundle.output_dir / recovery.CANONICAL_AUTHORIZATION_NAME)
        or (
            require_adapter_binding
            and audit.get("audit_adapter") != _audit_adapter_binding(bundle)
        )
    ):
        raise ValueError("Han v3 dtype-recovery completion audit binding differs")
    return audit


def _publish_completed_recovery_audit(
    bundle: recovery.RecoveryBundle,
    authorization: Mapping[str, Any],
    completed: list[dict[str, Any]],
    base_audit: Mapping[str, Any],
    *,
    superseded_unbound_audit: Mapping[str, Any] | None = None,
    replace_canonical: bool = False,
) -> dict[str, Any]:
    """Publish the frozen recovery schema without entering a simulation route."""

    authorization_path = bundle.output_dir / recovery.CANONICAL_AUTHORIZATION_NAME
    history_path = bundle.repo_root / str(authorization["history_manifest_path"])
    history = _read_hashed_json(
        history_path,
        schema=recovery.RECOVERY_HISTORY_SCHEMA,
        status="immutable_stopped_state_snapshot",
        label="dtype-recovery history manifest",
    )
    if (
        v3.sha256_file(history_path) != authorization["history_manifest_sha256"]
        or history["payload_sha256"] != authorization["history_payload_sha256"]
    ):
        raise ValueError("Han v3 recovery-audit history binding differs")

    authorized_ids = sorted(
        str(value) for value in authorization["recoverable_unit_ids"]
    )
    completed_ids = sorted(str(row["unit_id"]) for row in completed)
    if completed_ids != authorized_ids:
        raise ValueError("Han v3 recovery-audit completed unit IDs differ")
    invocation_started_ns = time_ns()
    invocation_prestate = {
        "authorization_payload_sha256": authorization["payload_sha256"],
        "already_completed": completed,
        "still_incomplete": [],
        "invocation_started_unix_ns": invocation_started_ns,
        "worker_count": 0,
        "audit_only_no_simulation_route": True,
    }
    invocation_id = _canonical_sha256(invocation_prestate)
    result_columns = (
        "unit_id",
        "status",
        "completion_path",
        "elapsed_seconds",
        "error_type",
        "error",
        "recovery_authorization_sha256",
    )
    snapshot = pd.DataFrame(columns=result_columns)
    snapshot_path = (
        bundle.output_dir / "recovery_run_snapshots" / f"audit_only_{invocation_id}.tsv"
    )
    v3._atomic_frame(snapshot_path, snapshot)
    payload: dict[str, Any] = {
        "schema": recovery.RECOVERY_COMPLETION_AUDIT_SCHEMA,
        "status": "passed",
        "invocation_id": invocation_id,
        "invocation_started_unix_ns": invocation_started_ns,
        "dry_plan_payload_sha256": authorization["dry_plan_payload_sha256"],
        "authorization_path": authorization_path.relative_to(
            bundle.repo_root
        ).as_posix(),
        "authorization_sha256": v3.sha256_file(authorization_path),
        "canonical_authorization_path": authorization_path.relative_to(
            bundle.repo_root
        ).as_posix(),
        "canonical_authorization_sha256": v3.sha256_file(authorization_path),
        "history_manifest_path": history_path.relative_to(bundle.repo_root).as_posix(),
        "history_manifest_sha256": v3.sha256_file(history_path),
        "history_payload_sha256": history["payload_sha256"],
        "history_file_count": int(history["file_count"]),
        "history_total_bytes": int(history["total_bytes"]),
        "manifest_sha256": authorization["manifest_sha256"],
        "seed_plan_sha256": authorization["seed_plan_sha256"],
        "original_implementation_sources": authorization[
            "original_implementation_sources"
        ],
        "recovery_sources": authorization["recovery_sources"],
        "audit_adapter": _audit_adapter_binding(bundle),
        "run_snapshot_path": snapshot_path.relative_to(bundle.repo_root).as_posix(),
        "run_snapshot_sha256": v3.sha256_file(snapshot_path),
        "initial_authorized_recovery_units": len(authorized_ids),
        "invocation_scheduled_units": 0,
        "pre_invocation_already_recovered_unit_ids": completed_ids,
        "pre_invocation_incomplete_unit_ids": [],
        "pre_recovery_incomplete_unit_ids": authorized_ids,
        "post_recovery_recovered_unit_ids": completed_ids,
        "post_recovery_incomplete_unit_ids": [],
        "invocation_failed_unit_ids": [],
        "validated_recovered_completions": len(completed),
        "remaining_incomplete_units": 0,
        "invocation_failed_units": 0,
        "rows": completed,
        "post_run_base_audit_status": base_audit["status"],
        "post_run_base_audit_accepted_units": int(base_audit["accepted_units"]),
        "post_run_base_audit_payload_sha256": base_audit["payload_sha256"],
        "post_run_base_audit": dict(base_audit),
        "original_attempt_seed_timeout_budget_and_reconciliation_semantics_preserved": (
            True
        ),
    }
    if superseded_unbound_audit is not None:
        payload["superseded_unbound_canonical_audit"] = dict(superseded_unbound_audit)
    payload["payload_sha256"] = _canonical_sha256(payload)
    invocation_audit_path = (
        bundle.output_dir
        / "recovery_completion_audits"
        / f"audit_only_{invocation_id}.json"
    )
    _immutable_json(invocation_audit_path, payload)
    canonical_path = bundle.output_dir / recovery.CANONICAL_AUDIT_NAME
    if replace_canonical:
        _atomic_json(canonical_path, payload)
    else:
        _immutable_json(canonical_path, payload)
    return payload


def _archive_unbound_canonical_audit(
    bundle: recovery.RecoveryBundle,
    authorization: Mapping[str, Any],
    base_audit: Mapping[str, Any],
) -> dict[str, Any]:
    """Preserve and bind the one adapter-generated audit lacking source fields."""

    canonical_path = bundle.output_dir / recovery.CANONICAL_AUDIT_NAME
    prior = _validate_recovery_completion_audit(
        bundle,
        authorization,
        base_audit,
        require_adapter_binding=False,
    )
    if (
        "audit_adapter" in prior
        or int(prior.get("invocation_scheduled_units", -1)) != 0
        or prior.get("pre_invocation_incomplete_unit_ids") != []
        or not Path(str(prior.get("run_snapshot_path", ""))).name.startswith(
            "audit_only_"
        )
    ):
        raise ValueError(
            "Han v3 recovery-audit refuses to replace an unexpected canonical audit"
        )
    prior_file_sha256 = v3.sha256_file(canonical_path)
    archive_path = (
        bundle.output_dir
        / "recovery_completion_audits"
        / f"superseded_unbound_{prior_file_sha256}.json"
    )
    if archive_path.is_file():
        if v3.sha256_file(archive_path) != prior_file_sha256:
            raise ValueError("Han v3 recovery-audit superseded archive differs")
    else:
        archive_path.parent.mkdir(parents=True, exist_ok=True)
        temporary = archive_path.with_name(f"{archive_path.name}.tmp.{os.getpid()}")
        try:
            shutil.copyfile(canonical_path, temporary)
            if v3.sha256_file(temporary) != prior_file_sha256:
                raise ValueError(
                    "Han v3 recovery-audit superseded copy checksum differs"
                )
            os.replace(temporary, archive_path)
        finally:
            temporary.unlink(missing_ok=True)
    return {
        "reason": "canonical audit predated direct audit-adapter source binding",
        "path": archive_path.relative_to(bundle.repo_root).as_posix(),
        "sha256": prior_file_sha256,
        "payload_sha256": prior["payload_sha256"],
        "preserved_byte_for_byte": True,
    }


def finalize_recovery_audit(
    bundle: recovery.RecoveryBundle,
    *,
    expected_plan_sha256: str,
) -> dict[str, Any]:
    """Finalize audits for an already-complete authorized recovery, idempotently."""

    authorization_path = bundle.output_dir / recovery.CANONICAL_AUTHORIZATION_NAME
    authorization = recovery._load_authorization(bundle, authorization_path)
    if authorization.get("dry_plan_payload_sha256") != str(expected_plan_sha256):
        raise ValueError("Han v3 recovery-audit requested dry-plan SHA differs")
    completed, incomplete = recovery._authorized_current_state(bundle, authorization)
    authorized_ids = sorted(
        str(value) for value in authorization["recoverable_unit_ids"]
    )
    completed_ids = sorted(str(row["unit_id"]) for row in completed)
    if incomplete or completed_ids != authorized_ids:
        raise ValueError(
            "Han v3 recovery-audit finalizer requires every authorized unit complete; "
            "it refuses to launch simulations"
        )

    canonical_recovery_audit = bundle.output_dir / recovery.CANONICAL_AUDIT_NAME
    with patched_recovery_audit_aggregation():
        base_audit = v3.audit_stage(bundle.direct, v3.STAGE_SIMULATE, write=True)

    if not canonical_recovery_audit.is_file():
        _publish_completed_recovery_audit(bundle, authorization, completed, base_audit)
    else:
        current = _read_hashed_json(
            canonical_recovery_audit,
            schema=recovery.RECOVERY_COMPLETION_AUDIT_SCHEMA,
            status="passed",
            label="dtype-recovery completion audit",
        )
        if current.get("audit_adapter") != _audit_adapter_binding(bundle):
            superseded = _archive_unbound_canonical_audit(
                bundle, authorization, base_audit
            )
            _publish_completed_recovery_audit(
                bundle,
                authorization,
                completed,
                base_audit,
                superseded_unbound_audit=superseded,
                replace_canonical=True,
            )

    recovery_audit = _validate_recovery_completion_audit(
        bundle, authorization, base_audit
    )
    base_audit_path = bundle.output_dir / "all_units_audit.json"
    base_authorization_path = bundle.output_dir / "han_direct_v3_authorization.json"
    absent_rows = sum(
        row["reconciliation_path"] is None and row["reconciliation_sha256"] is None
        for row in base_audit["rows"]
    )
    present_rows = len(base_audit["rows"]) - absent_rows
    normalization: dict[str, Any] = {
        "schema": AUDIT_NORMALIZATION_SCHEMA,
        "status": "passed",
        "manifest_sha256": v3.sha256_file(bundle.direct.manifest_path),
        "seed_plan_sha256": v3.sha256_file(bundle.direct.seed_plan_path),
        "adapter_sources": _adapter_source_records(bundle.repo_root),
        "frozen_recovery_sources": dict(bundle.recovery_sources),
        "recovery_authorization_path": authorization_path.relative_to(
            bundle.repo_root
        ).as_posix(),
        "recovery_authorization_sha256": v3.sha256_file(authorization_path),
        "recovery_authorization_payload_sha256": authorization["payload_sha256"],
        "recovery_completion_audit_path": canonical_recovery_audit.relative_to(
            bundle.repo_root
        ).as_posix(),
        "recovery_completion_audit_sha256": v3.sha256_file(canonical_recovery_audit),
        "recovery_completion_audit_payload_sha256": recovery_audit["payload_sha256"],
        "base_all_units_audit_path": base_audit_path.relative_to(
            bundle.repo_root
        ).as_posix(),
        "base_all_units_audit_sha256": v3.sha256_file(base_audit_path),
        "base_all_units_audit_payload_sha256": base_audit["payload_sha256"],
        "base_authorization_path": base_authorization_path.relative_to(
            bundle.repo_root
        ).as_posix(),
        "base_authorization_sha256": v3.sha256_file(base_authorization_path),
        "normalized_fields": list(STRUCTURAL_ABSENCE_FIELDS),
        "normalization_rule": (
            "pandas NaN to JSON null only when "
            "interrupted_nonterminal_dispositions equals zero"
        ),
        "structurally_absent_rows": int(absent_rows),
        "reconciliation_present_rows": int(present_rows),
        "validated_rows": len(base_audit["rows"]),
        "all_other_nonfinite_values_rejected": True,
        "simulation_launched": False,
        "scientific_rows_changed": False,
        "seed_retry_timeout_budget_or_reconciliation_semantics_changed": False,
    }
    normalization["payload_sha256"] = _canonical_sha256(normalization)
    normalization_path = bundle.output_dir / CANONICAL_NORMALIZATION_AUDIT_NAME
    _immutable_json(normalization_path, normalization)
    return normalization


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Finalize Han v3 recovery audits after every authorized unit is complete."
        )
    )
    parser.add_argument("command", choices=("finalize",))
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument(
        "--output-dir", type=Path, default=Path(v3.DEFAULT_OUTPUT_RELATIVE_PATH)
    )
    parser.add_argument(
        "--slim-bin", type=Path, default=Path(".native-stdpopsim/bin/slim")
    )
    parser.add_argument("--expected-plan-sha256", required=True)
    return parser


def _default_paths(args: argparse.Namespace) -> tuple[Path, Path, Path]:
    root = Path(args.repo_root).resolve()
    output = Path(args.output_dir)
    if not output.is_absolute():
        output = root / output
    slim = Path(args.slim_bin)
    if not slim.is_absolute():
        slim = root / slim
    return root, output.resolve(), slim.resolve()


def main(argv: Sequence[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    root, output, slim = _default_paths(args)
    bundle = recovery.load_recovery_bundle(
        output / v3.DEFAULT_PLAN_MANIFEST_NAME,
        repo_root=root,
        slim_path=slim,
    )
    payload = finalize_recovery_audit(
        bundle,
        expected_plan_sha256=args.expected_plan_sha256,
    )
    print(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
