"""Hash-bound, simulation-free selected-campaign postproduction adapter.

The frozen selected-production readiness builder re-audits the Han direct-v3
cache on every postproduction action.  Mixed interruption-journal presence can
make pandas represent the two structural absences as ``NaN`` before the frozen
JSON encoder sees them.  This additive adapter authenticates the completed Han
normalization evidence and enters the narrow recovery-audit normalization
context for every readiness build or rebuild.

Only authorization, validation, aggregation, analysis, and reporting are
exposed.  There is deliberately no simulation route.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from collections.abc import Iterator, Mapping, Sequence
from contextlib import contextmanager
from pathlib import Path
from typing import Any

from . import focused_han_direct_v3 as han_v3
from . import focused_han_direct_v3_recovery as recovery
from . import focused_han_direct_v3_recovery_audit as recovery_audit
from . import focused_selected_production as selected

POSTPRODUCTION_AUTHORIZATION_SCHEMA = (
    "gamma-smc.focused-selected-postproduction-authorization/v1"
)
MODULE_SOURCE_PATH = "python/gamma_smc_aou/focused_selected_postproduction.py"
WRAPPER_SOURCE_PATH = "scripts/run_focused_selected_postproduction.py"
POSTPRODUCTION_SOURCE_PATHS = (MODULE_SOURCE_PATH, WRAPPER_SOURCE_PATH)
DELEGATED_SOURCE_PATHS = (
    selected.MODULE_SOURCE_PATH,
    selected.WRAPPER_SOURCE_PATH,
)
DEFAULT_AUTHORIZATION_RELATIVE_PATH = (
    "focused_selection_EAS_sim/calibration/focused_selected_production/"
    "postproduction_authorization.json"
)
DEFAULT_NORMALIZATION_RELATIVE_PATH = (
    f"{selected.DEFAULT_HAN_V3_OUTPUT_RELATIVE_PATH}/"
    f"{recovery_audit.CANONICAL_NORMALIZATION_AUDIT_NAME}"
)
POSTPRODUCTION_ACTIONS = (
    "authorize",
    "validate",
    "aggregate-validate",
    "aggregate",
    "analyze",
    "report",
)
NORMALIZATION_RULE = (
    "pandas NaN to JSON null only when interrupted_nonterminal_dispositions equals zero"
)
EXPECTED_STRUCTURALLY_ABSENT_ROWS = 35
EXPECTED_RECONCILIATION_PRESENT_ROWS = 25


def _canonical_sha256(value: Any) -> str:
    return han_v3._canonical_sha256(value)


def canonical_text_sha256(path: str | Path) -> str:
    """Hash strict UTF-8 text with CRLF canonicalized to LF.

    UTF-8 BOMs, invalid UTF-8, and bare carriage returns are rejected.  No
    Unicode normalization is performed.
    """

    source_path = Path(path)
    try:
        raw = source_path.read_bytes()
    except OSError as error:
        raise ValueError(
            f"postproduction source is unreadable: {source_path}"
        ) from error
    if raw.startswith(b"\xef\xbb\xbf"):
        raise ValueError(f"postproduction source has a UTF-8 BOM: {source_path}")
    try:
        text = raw.decode("utf-8", errors="strict")
    except UnicodeDecodeError as error:
        raise ValueError(
            f"postproduction source is not strict UTF-8: {source_path}"
        ) from error
    without_crlf = text.replace("\r\n", "")
    if "\r" in without_crlf:
        raise ValueError(
            f"postproduction source contains a bare carriage return: {source_path}"
        )
    canonical = text.replace("\r\n", "\n").encode("utf-8")
    return hashlib.sha256(canonical).hexdigest()


def _inside_repo(path: str | Path, repo_root: Path, *, label: str) -> Path:
    candidate = Path(path)
    if not candidate.is_absolute():
        candidate = repo_root / candidate
    lexical = Path(os.path.abspath(candidate))
    try:
        relative = lexical.relative_to(repo_root)
    except ValueError as error:
        raise ValueError(f"{label} must be inside repo_root") from error
    current = repo_root
    for part in relative.parts:
        current = current / part
        if current.is_symlink():
            raise ValueError(f"{label} path contains a symlink")
    resolved = lexical.resolve()
    try:
        resolved.relative_to(repo_root)
    except ValueError as error:
        raise ValueError(f"{label} must resolve inside repo_root") from error
    return resolved


def _canonical_source_records(
    repo_root: Path, relative_paths: Sequence[str]
) -> dict[str, str]:
    records: dict[str, str] = {}
    for relative in relative_paths:
        path = _inside_repo(repo_root / relative, repo_root, label="source")
        if not path.is_file() or path.is_symlink():
            raise ValueError(f"postproduction source is absent or unsafe: {relative}")
        records[relative] = canonical_text_sha256(path)
    return records


def _validate_canonical_source_records(
    records: Any,
    *,
    repo_root: Path,
    relative_paths: Sequence[str],
    label: str,
) -> dict[str, str]:
    expected = _canonical_source_records(repo_root, relative_paths)
    if not isinstance(records, Mapping) or dict(records) != expected:
        raise ValueError(f"{label} canonical source binding differs")
    return expected


def _read_json(path: Path, *, label: str) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"{label} is unreadable") from error
    if not isinstance(payload, dict):
        raise TypeError(f"{label} is not a JSON object")
    return payload


def _validate_hashed_payload(
    payload: Mapping[str, Any],
    *,
    schema: str,
    status: str,
    label: str,
) -> None:
    digest = payload.get("payload_sha256")
    unhashed = dict(payload)
    unhashed.pop("payload_sha256", None)
    if (
        payload.get("schema") != schema
        or payload.get("status") != status
        or not isinstance(digest, str)
        or digest != _canonical_sha256(unhashed)
    ):
        raise ValueError(f"{label} is corrupt or incompatible")


def _file_record(path: Path, repo_root: Path) -> dict[str, Any]:
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"postproduction artifact is absent or unsafe: {path}")
    return {
        "path": path.relative_to(repo_root).as_posix(),
        "size_bytes": path.stat().st_size,
        "sha256": han_v3.sha256_file(path),
    }


def _canonical_artifact(
    normalization: Mapping[str, Any],
    *,
    path_key: str,
    sha_key: str,
    repo_root: Path,
    expected_relative: str,
    label: str,
) -> Path:
    if normalization.get(path_key) != expected_relative:
        raise ValueError(f"Han normalization {label} path differs")
    path = _inside_repo(repo_root / expected_relative, repo_root, label=label)
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"Han normalization {label} is absent or unsafe")
    if normalization.get(sha_key) != han_v3.sha256_file(path):
        raise ValueError(f"Han normalization {label} checksum differs")
    return path


def _validate_normalization_partition(normalization: Mapping[str, Any]) -> None:
    expected = {
        "structurally_absent_rows": EXPECTED_STRUCTURALLY_ABSENT_ROWS,
        "reconciliation_present_rows": EXPECTED_RECONCILIATION_PRESENT_ROWS,
    }
    for key, expected_value in expected.items():
        value = normalization.get(key)
        if type(value) is not int or value != expected_value:
            raise ValueError("Han recovery-audit normalization row partition differs")


def _base_reconciliation_partition(
    base_audit_payload: Mapping[str, Any],
) -> tuple[int, int]:
    rows = base_audit_payload.get("rows")
    if not isinstance(rows, list) or len(rows) != 60:
        raise ValueError("Han base all-unit reconciliation rows differ")
    absent = 0
    present = 0
    for row in rows:
        if not isinstance(row, Mapping):
            raise TypeError("Han base all-unit reconciliation row differs")
        count = row.get("interrupted_nonterminal_dispositions")
        path = row.get("reconciliation_path")
        digest = row.get("reconciliation_sha256")
        if type(count) is not int or count < 0:
            raise ValueError("Han base all-unit interruption count differs")
        if path is None and digest is None and count == 0:
            absent += 1
            continue
        if (
            count > 0
            and isinstance(path, str)
            and path.endswith("/supervision_reconciliation.tsv")
            and isinstance(digest, str)
            and len(digest) == 64
            and all(character in "0123456789abcdef" for character in digest)
        ):
            present += 1
            continue
        raise ValueError("Han base all-unit reconciliation reference differs")
    if (
        absent != EXPECTED_STRUCTURALLY_ABSENT_ROWS
        or present != EXPECTED_RECONCILIATION_PRESENT_ROWS
    ):
        raise ValueError("Han base all-unit reconciliation partition differs")
    return absent, present


def _validate_normalization_evidence_for_bundle(
    *,
    repo_root: Path,
    normalization_path: Path,
    slim_path: Path,
) -> dict[str, Any]:
    """Validate normalization and every current canonical Han artifact."""

    normalization = _read_json(
        normalization_path, label="Han recovery-audit normalization"
    )
    _validate_hashed_payload(
        normalization,
        schema=recovery_audit.AUDIT_NORMALIZATION_SCHEMA,
        status="passed",
        label="Han recovery-audit normalization",
    )
    if normalization_path.relative_to(repo_root).as_posix() != (
        DEFAULT_NORMALIZATION_RELATIVE_PATH
    ):
        raise ValueError("Han canonical normalization path differs")
    adapter_sources = _validate_canonical_source_records(
        normalization.get("adapter_sources"),
        repo_root=repo_root,
        relative_paths=recovery_audit.ADAPTER_SOURCE_PATHS,
        label="Han recovery-audit adapter",
    )

    manifest_path = _inside_repo(
        repo_root / selected.DEFAULT_HAN_V3_MANIFEST_RELATIVE_PATH,
        repo_root,
        label="Han manifest",
    )
    direct = han_v3.load_bundle(manifest_path, repo_root=repo_root, slim_path=slim_path)
    recovery_bundle = recovery.load_recovery_bundle(
        manifest_path, repo_root=repo_root, slim_path=slim_path
    )
    _validate_normalization_partition(normalization)
    if (
        normalization.get("manifest_sha256") != han_v3.sha256_file(manifest_path)
        or normalization.get("seed_plan_sha256")
        != han_v3.sha256_file(direct.seed_plan_path)
        or normalization.get("frozen_recovery_sources")
        != dict(recovery_bundle.recovery_sources)
        or normalization.get("normalized_fields")
        != list(recovery_audit.STRUCTURAL_ABSENCE_FIELDS)
        or normalization.get("normalization_rule") != NORMALIZATION_RULE
        or normalization.get("validated_rows") != 60
        or normalization.get("all_other_nonfinite_values_rejected") is not True
        or normalization.get("simulation_launched") is not False
        or normalization.get("scientific_rows_changed") is not False
        or normalization.get(
            "seed_retry_timeout_budget_or_reconciliation_semantics_changed"
        )
        is not False
    ):
        raise ValueError("Han recovery-audit normalization contract differs")

    recovery_authorization = _canonical_artifact(
        normalization,
        path_key="recovery_authorization_path",
        sha_key="recovery_authorization_sha256",
        repo_root=repo_root,
        expected_relative=selected.DEFAULT_HAN_V3_RECOVERY_AUTHORIZATION_RELATIVE_PATH,
        label="recovery authorization",
    )
    recovery_completion = _canonical_artifact(
        normalization,
        path_key="recovery_completion_audit_path",
        sha_key="recovery_completion_audit_sha256",
        repo_root=repo_root,
        expected_relative=selected.DEFAULT_HAN_V3_RECOVERY_AUDIT_RELATIVE_PATH,
        label="recovery completion audit",
    )
    base_audit = _canonical_artifact(
        normalization,
        path_key="base_all_units_audit_path",
        sha_key="base_all_units_audit_sha256",
        repo_root=repo_root,
        expected_relative=selected.DEFAULT_HAN_V3_FINAL_AUDIT_RELATIVE_PATH,
        label="base all-unit audit",
    )
    base_authorization = _canonical_artifact(
        normalization,
        path_key="base_authorization_path",
        sha_key="base_authorization_sha256",
        repo_root=repo_root,
        expected_relative=selected.DEFAULT_HAN_V3_AUTHORIZATION_RELATIVE_PATH,
        label="base authorization",
    )

    recovery_authorization_payload = _read_json(
        recovery_authorization, label="Han recovery authorization"
    )
    recovery_completion_payload = _read_json(
        recovery_completion, label="Han recovery completion audit"
    )
    base_audit_payload = _read_json(base_audit, label="Han base all-unit audit")
    base_authorization_payload = _read_json(
        base_authorization, label="Han base authorization"
    )
    observed_absent, observed_present = _base_reconciliation_partition(
        base_audit_payload
    )
    if (
        normalization["structurally_absent_rows"] != observed_absent
        or normalization["reconciliation_present_rows"] != observed_present
    ):
        raise ValueError("Han normalization and base-audit row partitions differ")
    for payload, schema, status, label in (
        (
            recovery_authorization_payload,
            recovery.RECOVERY_AUTHORIZATION_SCHEMA,
            "authorized_after_immutable_history_snapshot",
            "Han recovery authorization",
        ),
        (
            recovery_completion_payload,
            recovery.RECOVERY_COMPLETION_AUDIT_SCHEMA,
            "passed",
            "Han recovery completion audit",
        ),
        (
            base_audit_payload,
            han_v3.STAGE_AUDIT_SCHEMA,
            "passed",
            "Han base all-unit audit",
        ),
        (
            base_authorization_payload,
            han_v3.AUTHORIZATION_SCHEMA,
            "authorized",
            "Han base authorization",
        ),
    ):
        _validate_hashed_payload(payload, schema=schema, status=status, label=label)

    recorded_adapter = recovery_completion_payload.get("audit_adapter")
    if (
        not isinstance(recorded_adapter, Mapping)
        or recorded_adapter.get("sources") != adapter_sources
        or recorded_adapter.get("normalized_fields")
        != list(recovery_audit.STRUCTURAL_ABSENCE_FIELDS)
        or recorded_adapter.get("normalization_rule") != NORMALIZATION_RULE
        or recorded_adapter.get("audit_only_no_simulation_route") is not True
        or recorded_adapter.get("all_other_nonfinite_values_rejected") is not True
        or recorded_adapter.get("scientific_rows_changed") is not False
        or recorded_adapter.get(
            "seed_retry_timeout_budget_or_reconciliation_semantics_changed"
        )
        is not False
        or normalization.get("recovery_authorization_payload_sha256")
        != recovery_authorization_payload["payload_sha256"]
        or normalization.get("recovery_completion_audit_payload_sha256")
        != recovery_completion_payload["payload_sha256"]
        or normalization.get("base_all_units_audit_payload_sha256")
        != base_audit_payload["payload_sha256"]
        or recovery_completion_payload.get("post_run_base_audit") != base_audit_payload
        or recovery_completion_payload.get("post_run_base_audit_payload_sha256")
        != base_audit_payload["payload_sha256"]
        or base_authorization_payload.get("all_units_audit_payload_sha256")
        != base_audit_payload["payload_sha256"]
    ):
        raise ValueError("Han normalized canonical artifact binding differs")

    return {
        "normalization": {
            **_file_record(normalization_path, repo_root),
            "payload_sha256": normalization["payload_sha256"],
        },
        "recovery_audit_adapter_sources": adapter_sources,
        "canonical_han_artifacts": {
            "manifest": _file_record(manifest_path, repo_root),
            "seed_plan": _file_record(direct.seed_plan_path, repo_root),
            "recovery_authorization": {
                **_file_record(recovery_authorization, repo_root),
                "payload_sha256": recovery_authorization_payload["payload_sha256"],
            },
            "recovery_completion_audit": {
                **_file_record(recovery_completion, repo_root),
                "payload_sha256": recovery_completion_payload["payload_sha256"],
            },
            "base_all_units_audit": {
                **_file_record(base_audit, repo_root),
                "payload_sha256": base_audit_payload["payload_sha256"],
            },
            "base_authorization": {
                **_file_record(base_authorization, repo_root),
                "payload_sha256": base_authorization_payload["payload_sha256"],
            },
        },
    }


def _validate_preflight_binding(
    preflight: Mapping[str, Any], evidence: Mapping[str, Any]
) -> None:
    artifacts = evidence.get("canonical_han_artifacts", {})
    dtype_binding = preflight.get("han_v3", {}).get("dtype_recovery", {}).get("binding")
    if (
        not isinstance(artifacts, Mapping)
        or not isinstance(dtype_binding, Mapping)
        or dtype_binding.get("authorization", {}).get("sha256")
        != artifacts.get("recovery_authorization", {}).get("sha256")
        or dtype_binding.get("audit", {}).get("sha256")
        != artifacts.get("recovery_completion_audit", {}).get("sha256")
        or preflight.get("han_v3", {}).get("all_units_audit", {}).get("sha256")
        != artifacts.get("base_all_units_audit", {}).get("sha256")
        or preflight.get("han_v3", {})
        .get("final_operational_authorization", {})
        .get("sha256")
        != artifacts.get("base_authorization", {}).get("sha256")
    ):
        raise ValueError("Han selected-readiness preflight binding differs")


def build_postproduction_authorization(
    *,
    repo_root: str | Path,
    slim_path: str | Path,
    integration_manifest_path: str | Path = selected.DEFAULT_INTEGRATION_RELATIVE_PATH,
    normalization_path: str | Path = DEFAULT_NORMALIZATION_RELATIVE_PATH,
) -> dict[str, Any]:
    """Build a deterministic authorization for simulation-free postproduction."""

    root = Path(repo_root).resolve()
    if any("onedrive" in part.casefold() for part in root.parts):
        raise ValueError(
            "selected postproduction repo_root must not be inside OneDrive"
        )
    slim = _inside_repo(slim_path, root, label="SLiM binary")
    integration = _inside_repo(
        integration_manifest_path, root, label="selected integration manifest"
    )
    normalization = _inside_repo(
        normalization_path, root, label="Han normalization artifact"
    )

    evidence = _validate_normalization_evidence_for_bundle(
        repo_root=root,
        normalization_path=normalization,
        slim_path=slim,
    )
    with recovery_audit.patched_recovery_audit_aggregation():
        preflight = selected.build_aggregate_readiness(
            repo_root=root,
            eas_integration_manifest_path=integration,
            slim_path=slim,
            require_final=False,
            validate_eas_cache=False,
        )
    if (
        preflight.get("counts", {}).get("han_selected") != 60
        or preflight.get("han_v3", {}).get("all_units_audit") is None
        or preflight.get("han_v3", {}).get("final_operational_authorization") is None
        or preflight.get("han_v3", {}).get("dtype_recovery", {}).get("evidence_present")
        is not True
        or preflight.get("han_v3", {}).get("dtype_recovery", {}).get("binding") is None
    ):
        raise ValueError("selected postproduction Han preflight is incomplete")

    _validate_preflight_binding(preflight, evidence)
    integration_payload = _read_json(integration, label="selected integration manifest")
    payload: dict[str, Any] = {
        "schema": POSTPRODUCTION_AUTHORIZATION_SCHEMA,
        "status": "authorized",
        "postproduction_sources": _canonical_source_records(
            root, POSTPRODUCTION_SOURCE_PATHS
        ),
        "delegated_selected_sources": _canonical_source_records(
            root, DELEGATED_SOURCE_PATHS
        ),
        "selected_integration": {
            **_file_record(integration, root),
            "payload_sha256": integration_payload["payload_sha256"],
            "strict_loaded": True,
        },
        "han_normalization_evidence": evidence,
        "preflight_payload_sha256": preflight["payload_sha256"],
        "preflight_eas_cache_validation_deferred": True,
        "normalization_context": (
            "focused_han_direct_v3_recovery_audit.patched_recovery_audit_aggregation"
        ),
        "actions": [
            "validate",
            "aggregate-validate",
            "aggregate",
            "analyze",
            "report",
        ],
        "simulation_route_present": False,
        "standard_selected_simulate_eas_untouched": True,
        "every_readiness_revalidation_inside_normalization_context": True,
        "source_hash_semantics": {
            "encoding": "strict UTF-8",
            "newline_canonicalization": "CRLF to LF",
            "utf8_bom_rejected": True,
            "bare_cr_rejected": True,
            "unicode_normalization": False,
        },
    }
    payload["payload_sha256"] = _canonical_sha256(payload)
    return payload


def write_postproduction_authorization(
    path: str | Path, payload: Mapping[str, Any]
) -> Path:
    _validate_hashed_payload(
        payload,
        schema=POSTPRODUCTION_AUTHORIZATION_SCHEMA,
        status="authorized",
        label="selected postproduction authorization",
    )
    destination = Path(path).resolve()
    if any("onedrive" in part.casefold() for part in destination.parts):
        raise ValueError("postproduction authorization must not be inside OneDrive")
    content = (
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        if destination.read_bytes() != content:
            raise FileExistsError(
                f"immutable postproduction authorization differs: {destination}"
            )
        return destination
    temporary = destination.with_name(f".{destination.name}.{os.getpid()}.tmp")
    try:
        temporary.write_bytes(content)
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)
    return destination


def load_postproduction_authorization(
    path: str | Path,
    *,
    repo_root: str | Path,
    slim_path: str | Path,
) -> dict[str, Any]:
    root = Path(repo_root).resolve()
    authorization = _inside_repo(path, root, label="postproduction authorization")
    if (
        authorization.relative_to(root).as_posix()
        != DEFAULT_AUTHORIZATION_RELATIVE_PATH
    ):
        raise ValueError("canonical postproduction authorization path differs")
    payload = _read_json(authorization, label="selected postproduction authorization")
    _validate_hashed_payload(
        payload,
        schema=POSTPRODUCTION_AUTHORIZATION_SCHEMA,
        status="authorized",
        label="selected postproduction authorization",
    )
    current = build_postproduction_authorization(
        repo_root=root,
        slim_path=slim_path,
        integration_manifest_path=(
            root / str(payload.get("selected_integration", {}).get("path", ""))
        ),
        normalization_path=(
            root
            / str(
                payload.get("han_normalization_evidence", {})
                .get("normalization", {})
                .get("path", "")
            )
        ),
    )
    if current != payload:
        raise ValueError(
            "selected postproduction authorization differs from current evidence"
        )
    return payload


@contextmanager
def authorized_postproduction_context(
    *,
    repo_root: str | Path,
    slim_path: str | Path,
    authorization_path: str | Path = DEFAULT_AUTHORIZATION_RELATIVE_PATH,
) -> Iterator[dict[str, Any]]:
    authorization = load_postproduction_authorization(
        authorization_path, repo_root=repo_root, slim_path=slim_path
    )
    with recovery_audit.patched_recovery_audit_aggregation():
        yield authorization


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=POSTPRODUCTION_ACTIONS)
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument(
        "--slim-bin", type=Path, default=Path(".native-stdpopsim/bin/slim")
    )
    parser.add_argument("--authorization", type=Path)
    parser.add_argument("--integration-manifest", type=Path)
    parser.add_argument("--normalization", type=Path)
    parser.add_argument("--aggregate-readiness", type=Path)
    parser.add_argument("--empirical-tsv", type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    root = args.repo_root.resolve()
    slim = args.slim_bin
    if not slim.is_absolute():
        slim = root / slim
    authorization = args.authorization or root / DEFAULT_AUTHORIZATION_RELATIVE_PATH
    integration = (
        args.integration_manifest or root / selected.DEFAULT_INTEGRATION_RELATIVE_PATH
    )
    normalization = args.normalization or root / DEFAULT_NORMALIZATION_RELATIVE_PATH
    readiness = (
        args.aggregate_readiness
        or root / selected.DEFAULT_AGGREGATE_READINESS_RELATIVE_PATH
    )

    if args.action == "authorize":
        payload = build_postproduction_authorization(
            repo_root=root,
            slim_path=slim,
            integration_manifest_path=integration,
            normalization_path=normalization,
        )
        write_postproduction_authorization(authorization, payload)
        print(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False))
        return 0

    with authorized_postproduction_context(
        repo_root=root, slim_path=slim, authorization_path=authorization
    ) as postproduction_authorization:
        if args.action == "validate":
            print(
                json.dumps(
                    {
                        "status": "validated",
                        "postproduction_authorization_payload_sha256": (
                            postproduction_authorization["payload_sha256"]
                        ),
                        "simulation_route_present": False,
                    },
                    indent=2,
                    sort_keys=True,
                )
            )
            return 0
        if args.action == "aggregate-validate":
            payload = selected.build_aggregate_readiness(
                repo_root=root,
                eas_integration_manifest_path=integration,
                slim_path=slim,
                require_final=True,
                validate_eas_cache=True,
            )
            selected.write_aggregate_readiness(readiness, payload)
            print(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False))
            return 0
        if args.action == "aggregate":
            outputs = selected.aggregate_with_readiness(
                repo_root=root, slim_path=slim, readiness_path=readiness
            )
            print(
                json.dumps(
                    {
                        "status": "aggregated",
                        "outputs": [
                            path.resolve().relative_to(root).as_posix()
                            for path in outputs
                        ],
                        "postproduction_authorization_payload_sha256": (
                            postproduction_authorization["payload_sha256"]
                        ),
                    },
                    indent=2,
                    sort_keys=True,
                )
            )
            return 0
        if args.action == "analyze":
            payload = selected.analyze_with_readiness(
                repo_root=root,
                slim_path=slim,
                readiness_path=readiness,
                empirical_path=args.empirical_tsv,
            )
            print(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False))
            return 0
        if args.action == "report":
            payload = selected.report_with_readiness(
                repo_root=root, slim_path=slim, readiness_path=readiness
            )
            print(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False))
            return 0
    raise AssertionError(f"unhandled postproduction action: {args.action}")


__all__ = [
    "DEFAULT_AUTHORIZATION_RELATIVE_PATH",
    "DEFAULT_NORMALIZATION_RELATIVE_PATH",
    "MODULE_SOURCE_PATH",
    "POSTPRODUCTION_ACTIONS",
    "POSTPRODUCTION_AUTHORIZATION_SCHEMA",
    "WRAPPER_SOURCE_PATH",
    "authorized_postproduction_context",
    "build_postproduction_authorization",
    "canonical_text_sha256",
    "load_postproduction_authorization",
    "main",
    "write_postproduction_authorization",
]


if __name__ == "__main__":
    raise SystemExit(main())
