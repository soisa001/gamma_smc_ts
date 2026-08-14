"""Hash-bound, simulation-free selected-campaign postproduction adapter.

The frozen selected-production readiness builder re-audits the Han direct-v3
cache on every postproduction action.  Mixed interruption-journal presence can
make pandas represent the two structural absences as ``NaN`` before the frozen
JSON encoder sees them.  This additive adapter authenticates the completed Han
normalization evidence and enters the narrow recovery-audit normalization
context for every readiness build or rebuild.

Only provenance attestation, authorization, validation, aggregation, decoding,
analysis, and reporting are exposed.  There is deliberately no simulation
route.
"""

from __future__ import annotations

import argparse
import base64
import copy
import hashlib
import json
import os
import threading
from collections.abc import Iterator, Mapping, Sequence
from contextlib import contextmanager
from pathlib import Path
from typing import Any

import pandas as pd

from . import focused_han_direct_v3 as han_v3
from . import focused_han_direct_v3_recovery as recovery
from . import focused_han_direct_v3_recovery_audit as recovery_audit
from . import focused_selection_campaign as campaign
from . import focused_selection_simulation as simulation
from . import focused_selected_production as selected

POSTPRODUCTION_AUTHORIZATION_SCHEMA = (
    "gamma-smc.focused-selected-postproduction-authorization/v2"
)
EAS_CONTRACT_ATTESTATION_SCHEMA = (
    "gamma-smc.focused-eas-selected-contract-attestation/v1"
)
LEGACY_POSTPRODUCTION_AUTHORIZATION_SCHEMA = (
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
DEFAULT_EAS_ATTESTATION_RELATIVE_PATH = (
    "focused_selection_EAS_sim/calibration/focused_selected_production/"
    "eas_contract_attestation.json"
)
DEFAULT_EAS_EVIDENCE_RELATIVE_DIR = (
    "focused_selection_EAS_sim/calibration/focused_selected_production/"
    "eas_contract_attestation_sources"
)
DEFAULT_SUPERSEDED_AUTHORIZATION_RELATIVE_DIR = (
    "focused_selection_EAS_sim/calibration/focused_selected_production/"
    "superseded_postproduction_authorizations"
)
DEFAULT_DECODER_RELATIVE_PATH = "bin/gamma_smc"
DEFAULT_DECODE_WORKERS = 20
DEFAULT_DECODE_THREADS = 1
POSTPRODUCTION_ACTIONS = (
    "attest-eas",
    "supersede-legacy",
    "authorize",
    "validate",
    "aggregate-validate",
    "aggregate",
    "decode",
    "analyze",
    "report",
)
NORMALIZATION_RULE = (
    "pandas NaN to JSON null only when interrupted_nonterminal_dispositions equals zero"
)
EXPECTED_STRUCTURALLY_ABSENT_ROWS = 35
EXPECTED_RECONCILIATION_PRESENT_ROWS = 25

EAS_ATTRIBUTE_SOURCE_PATHS = (
    "python/gamma_smc_aou/.gitattributes",
    "scripts/.gitattributes",
)
EAS_ATTRIBUTE_CONTRACT_PATHS = tuple(
    ("implementation", "sources", relative)
    for relative in EAS_ATTRIBUTE_SOURCE_PATHS
)
EAS_CANONICAL_ATTRIBUTE_HASH_PAIR = (
    "3b9653aa38de61bc0517a81d48d20b8c008b6b3c936ac7de949bf1e0a08160b9",
    "43e4e5442a9a3117ecc24c01c5fa896e5e045e352ed2d2b4435d62887204ff0b",
)
EAS_TRANSIENT_ATTRIBUTE_HASH_PAIR = (
    "05590d060f3095ed805868a185463d12846c18f0608f822a70a2762a0f50f185",
    "6fd1d6724508c3574e07045726ce91f8a70ee51977256bc9135c0755018b2bd2",
)
EAS_EXPECTED_VARIANT_COUNTS = {"canonical": 20, "transient_audit_text": 40}
EAS_TRANSIENT_ATTRIBUTE_EVIDENCE = {
    "python/gamma_smc_aou/.gitattributes": {
        "archive_name": "python_gamma_smc_aou.gitattributes.transient",
        "content_base64": (
            "Zm9jdXNlZF9zZWxlY3RlZF9wcm9kdWN0aW9uLnB5IC10ZXh0CmZvY3VzZWRfaGFu"
            "X2RpcmVjdF92M19yZWNvdmVyeS5weSAtdGV4dApmb2N1c2VkX2hhbl9kaXJlY3Rf"
            "djNfcmVjb3ZlcnlfYXVkaXQucHkgLXRleHQK"
        ),
        "sha256": EAS_TRANSIENT_ATTRIBUTE_HASH_PAIR[0],
        "size_bytes": 123,
    },
    "scripts/.gitattributes": {
        "archive_name": "scripts.gitattributes.transient",
        "content_base64": (
            "cnVuX2ZvY3VzZWRfc2VsZWN0ZWRfcHJvZHVjdGlvbi5weSAtdGV4dApydW5fZm9j"
            "dXNlZF9oYW5fZGlyZWN0X3YzX3JlY292ZXJ5LnB5IC10ZXh0CnJ1bl9mb2N1c2Vk"
            "X2hhbl9kaXJlY3RfdjNfcmVjb3ZlcnlfYXVkaXQucHkgLXRleHQK"
        ),
        "sha256": EAS_TRANSIENT_ATTRIBUTE_HASH_PAIR[1],
        "size_bytes": 135,
    },
}
_PATCH_LOCK = threading.RLock()


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


def _load_superseded_legacy_authorization(repo_root: Path) -> dict[str, Any]:
    directory = _inside_repo(
        repo_root / DEFAULT_SUPERSEDED_AUTHORIZATION_RELATIVE_DIR,
        repo_root,
        label="superseded authorization directory",
    )
    candidates = (
        sorted(directory.glob("postproduction_authorization_v1__*.json"))
        if directory.is_dir()
        else []
    )
    if len(candidates) != 1:
        raise ValueError(
            "exactly one superseded v1 postproduction authorization is required"
        )
    path = candidates[0]
    if path.is_symlink() or not path.is_file():
        raise ValueError("superseded v1 postproduction authorization is unsafe")
    payload = _read_json(path, label="superseded v1 postproduction authorization")
    _validate_hashed_payload(
        payload,
        schema=LEGACY_POSTPRODUCTION_AUTHORIZATION_SCHEMA,
        status="authorized",
        label="superseded v1 postproduction authorization",
    )
    expected_name = (
        "postproduction_authorization_v1__"
        f"{payload['payload_sha256']}.json"
    )
    if path.name != expected_name:
        raise ValueError("superseded v1 authorization filename binding differs")
    return {
        **_file_record(path, repo_root),
        "payload_sha256": payload["payload_sha256"],
        "schema": LEGACY_POSTPRODUCTION_AUTHORIZATION_SCHEMA,
        "status": "superseded_by_v2_without_byte_changes",
        "original_path": DEFAULT_AUTHORIZATION_RELATIVE_PATH,
        "preservation_operation": "same_volume_atomic_rename",
        "source_bytes_modified": False,
        "bytes_deleted": False,
        "original_path_vacated_by_rename": True,
    }


def supersede_legacy_postproduction_authorization(
    *, repo_root: str | Path
) -> dict[str, Any]:
    """Recoverably archive the exact v1 bytes before publishing v2."""

    root = Path(repo_root).resolve()
    canonical = _inside_repo(
        root / DEFAULT_AUTHORIZATION_RELATIVE_PATH,
        root,
        label="canonical postproduction authorization",
    )
    if canonical.is_file():
        if canonical.is_symlink():
            raise ValueError("canonical postproduction authorization is unsafe")
        payload = _read_json(canonical, label="canonical postproduction authorization")
        schema = payload.get("schema")
        if schema == POSTPRODUCTION_AUTHORIZATION_SCHEMA:
            _validate_hashed_payload(
                payload,
                schema=POSTPRODUCTION_AUTHORIZATION_SCHEMA,
                status="authorized",
                label="canonical v2 postproduction authorization",
            )
            return _load_superseded_legacy_authorization(root)
        if schema != LEGACY_POSTPRODUCTION_AUTHORIZATION_SCHEMA:
            raise ValueError("canonical postproduction authorization schema is unknown")
        _validate_hashed_payload(
            payload,
            schema=LEGACY_POSTPRODUCTION_AUTHORIZATION_SCHEMA,
            status="authorized",
            label="canonical v1 postproduction authorization",
        )
        archive = _inside_repo(
            root
            / DEFAULT_SUPERSEDED_AUTHORIZATION_RELATIVE_DIR
            / (
                "postproduction_authorization_v1__"
                f"{payload['payload_sha256']}.json"
            ),
            root,
            label="superseded v1 authorization archive",
        )
        if archive.exists():
            raise FileExistsError(
                "superseded v1 archive already exists; refusing to overwrite or "
                "remove either copy"
            )
        if canonical.anchor.casefold() != archive.anchor.casefold():
            raise ValueError("legacy authorization archive must remain on one volume")
        archive.parent.mkdir(parents=True, exist_ok=True)
        before = canonical.read_bytes()
        os.replace(canonical, archive)
        if canonical.exists() or archive.read_bytes() != before:
            raise RuntimeError("legacy authorization atomic preservation failed")
    return _load_superseded_legacy_authorization(root)


def _file_record(path: Path, repo_root: Path) -> dict[str, Any]:
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"postproduction artifact is absent or unsafe: {path}")
    return {
        "path": path.relative_to(repo_root).as_posix(),
        "size_bytes": path.stat().st_size,
        "sha256": han_v3.sha256_file(path),
    }


def _immutable_bytes(path: Path, content: bytes, *, label: str) -> Path:
    """Publish exact evidence bytes once, or prove an identical prior write."""

    if path.exists():
        if not path.is_file() or path.is_symlink() or path.read_bytes() != content:
            raise FileExistsError(f"immutable {label} differs: {path}")
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        temporary.write_bytes(content)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)
    return path


def _transient_attribute_evidence_records(repo_root: Path) -> dict[str, Any]:
    """Return the fixed byte archives that explain the transient hash pair."""

    records: dict[str, Any] = {}
    evidence_root = repo_root / DEFAULT_EAS_EVIDENCE_RELATIVE_DIR
    for source_path in EAS_ATTRIBUTE_SOURCE_PATHS:
        fixed = EAS_TRANSIENT_ATTRIBUTE_EVIDENCE[source_path]
        try:
            content = base64.b64decode(
                str(fixed["content_base64"]), validate=True
            )
        except (ValueError, TypeError) as error:
            raise RuntimeError("embedded transient attribute evidence is invalid") from error
        digest = hashlib.sha256(content).hexdigest()
        if digest != fixed["sha256"] or len(content) != fixed["size_bytes"]:
            raise RuntimeError("embedded transient attribute evidence checksum differs")
        archive_path = evidence_root / str(fixed["archive_name"])
        records[source_path] = {
            "path": archive_path.relative_to(repo_root).as_posix(),
            "size_bytes": len(content),
            "sha256": digest,
            "encoding": "UTF-8",
            "line_endings": "LF",
            "content_base64_sha256": hashlib.sha256(
                str(fixed["content_base64"]).encode("ascii")
            ).hexdigest(),
        }
    return records


def _validate_current_canonical_attributes(repo_root: Path) -> dict[str, Any]:
    records: dict[str, Any] = {}
    observed: list[str] = []
    for relative in EAS_ATTRIBUTE_SOURCE_PATHS:
        path = _inside_repo(repo_root / relative, repo_root, label="attribute source")
        if not path.is_file() or path.is_symlink():
            raise ValueError(f"canonical attribute source is absent or unsafe: {relative}")
        observed.append(han_v3.sha256_file(path))
        records[relative] = _file_record(path, repo_root)
    if tuple(observed) != EAS_CANONICAL_ATTRIBUTE_HASH_PAIR:
        raise ValueError("current canonical EAS attribute hash pair differs")
    return records


def _contract_attribute_pair(contract: Mapping[str, Any]) -> tuple[str, str]:
    implementation = contract.get("implementation")
    if not isinstance(implementation, Mapping):
        raise TypeError("EAS stored contract implementation is malformed")
    sources = implementation.get("sources")
    if not isinstance(sources, Mapping):
        raise TypeError("EAS stored contract source manifest is malformed")
    pair: list[str] = []
    for relative in EAS_ATTRIBUTE_SOURCE_PATHS:
        value = sources.get(relative)
        if not isinstance(value, str):
            raise ValueError(f"EAS stored contract lacks attribute hash: {relative}")
        pair.append(value)
    return tuple(pair)  # type: ignore[return-value]


def _deep_difference_paths(
    left: Any, right: Any, prefix: tuple[str, ...] = ()
) -> list[tuple[str, ...]]:
    """Enumerate exact JSON-structure differences without normalizing values."""

    if isinstance(left, Mapping) and isinstance(right, Mapping):
        paths: list[tuple[str, ...]] = []
        keys = sorted(set(left) | set(right), key=str)
        for key in keys:
            child = prefix + (str(key),)
            if key not in left or key not in right:
                paths.append(child)
            else:
                paths.extend(_deep_difference_paths(left[key], right[key], child))
        return paths
    if isinstance(left, list) and isinstance(right, list):
        paths = []
        for index in range(max(len(left), len(right))):
            child = prefix + (str(index),)
            if index >= len(left) or index >= len(right):
                paths.append(child)
            else:
                paths.extend(_deep_difference_paths(left[index], right[index], child))
        return paths
    return [] if type(left) is type(right) and left == right else [prefix]


def _normalized_eas_contract(
    contract: Mapping[str, Any], expected: Mapping[str, Any]
) -> dict[str, Any]:
    normalized = copy.deepcopy(dict(contract))
    expected_pair = _contract_attribute_pair(expected)
    if expected_pair != EAS_CANONICAL_ATTRIBUTE_HASH_PAIR:
        raise ValueError("rebuilt EAS contract does not use the canonical attribute pair")
    sources = normalized.get("implementation", {}).get("sources")
    if not isinstance(sources, dict):
        raise TypeError("EAS contract source manifest cannot be normalized")
    for relative, digest in zip(EAS_ATTRIBUTE_SOURCE_PATHS, expected_pair, strict=True):
        sources[relative] = digest
    return normalized


def _classify_eas_contract_variant(
    pair: tuple[str, str], differences: Sequence[tuple[str, ...]]
) -> str:
    difference_set = set(differences)
    allowed_set = set(EAS_ATTRIBUTE_CONTRACT_PATHS)
    if pair == EAS_CANONICAL_ATTRIBUTE_HASH_PAIR:
        if difference_set:
            raise ValueError("canonical EAS contract differs from current contract")
        return "canonical"
    if pair == EAS_TRANSIENT_ATTRIBUTE_HASH_PAIR:
        if difference_set != allowed_set:
            raise ValueError("transient EAS contract differs outside two attributes")
        return "transient_audit_text"
    raise ValueError("unknown or mixed EAS attribute hash pair")


def _unit_file_state(unit_dir: Path) -> dict[str, tuple[int, int]]:
    return {
        path.relative_to(unit_dir).as_posix(): (
            path.stat().st_size,
            path.stat().st_mtime_ns,
        )
        for path in unit_dir.rglob("*")
        if path.is_file()
    }


def _rebuild_current_expected_contract(
    bundle: selected.SelectedProductionBundle,
    unit: Mapping[str, Any],
) -> dict[str, Any]:
    """Capture the contract independently rebuilt by the current base executor.

    The completion validator alone is replaced, after a separate base validation
    has already succeeded.  Since the terminal completion exists, no simulation
    or artifact-writing branch can be entered.
    """

    captured: dict[str, Any] = {}
    original = simulation._validate_completion

    def capture(
        artifacts: simulation.SimulationArtifacts, contract: Mapping[str, Any]
    ) -> simulation.SimulationArtifacts:
        if captured:
            raise RuntimeError("current EAS contract was rebuilt more than once")
        captured["contract"] = copy.deepcopy(dict(contract))
        return simulation.SimulationArtifacts(
            **{**artifacts.__dict__, "cache_hit": True}
        )

    with _PATCH_LOCK:
        simulation._validate_completion = capture
        try:
            artifacts = simulation._simulate_unit_unlocked(
                unit,
                bundle.repo_root,
                bundle.campaign_dir,
                slim_path=bundle.slim_path,
                selected_max_draws=selected.DEFAULT_EAS_MAX_DRAWS,
                slim_scaling_factor=simulation.DEFAULT_SLIM_SCALING_FACTOR,
                slim_burn_in=simulation.DEFAULT_SLIM_BURN_IN,
                selected_draw_timeout_seconds=selected.DEFAULT_DRAW_TIMEOUT_SECONDS,
                selected_cumulative_timeout_seconds=(
                    selected.DEFAULT_EAS_MAX_DRAWS
                    * selected.DEFAULT_DRAW_TIMEOUT_SECONDS
                ),
            )
        finally:
            simulation._validate_completion = original
    if artifacts.cache_hit is not True or set(captured) != {"contract"}:
        raise RuntimeError("current EAS contract rebuild left the cache-only route")
    return captured["contract"]


def _attested_output_records(
    *,
    unit_dir: Path,
    completion: Mapping[str, Any],
    repo_root: Path,
) -> dict[str, Any]:
    artifacts = simulation._completion_artifacts(unit_dir)
    paths = {
        "tree": artifacts.tree_path,
        "pair_table": artifacts.pair_table_path,
        "overall_pairs": artifacts.overall_pairs_path,
        "truth_profiles": artifacts.truth_profiles_path,
        "truth_class_summaries": artifacts.truth_class_summaries_path,
    }
    outputs = completion.get("outputs")
    if not isinstance(outputs, Mapping) or set(outputs) != set(paths):
        raise ValueError("EAS completion output manifest differs")
    records: dict[str, Any] = {}
    for label, path in paths.items():
        record = outputs.get(label)
        expected_relative = path.relative_to(unit_dir).as_posix()
        if (
            not isinstance(record, Mapping)
            or str(record.get("path", "")).replace("\\", "/")
            != expected_relative
            or not path.is_file()
            or path.is_symlink()
        ):
            raise ValueError(f"EAS completion output record differs: {label}")
        actual = _file_record(path, repo_root)
        if actual["sha256"] != record.get("sha256"):
            raise ValueError(f"EAS completion output checksum differs: {label}")
        records[label] = actual
    return records


def audit_eas_contract_attestation(
    bundle: selected.SelectedProductionBundle,
) -> dict[str, Any]:
    """Strictly audit all 60 EAS selected completions without writing anything."""

    _validate_current_canonical_attributes(bundle.repo_root)
    _transient_attribute_evidence_records(bundle.repo_root)
    canonical = pd.read_csv(bundle.execution_units_path, sep="\t")
    units = selected._eas_production_units(canonical)
    rows: list[dict[str, Any]] = []
    variant_counts = {key: 0 for key in EAS_EXPECTED_VARIANT_COUNTS}
    with selected.patched_selected_executor(bundle):
        for unit in units.to_dict(orient="records"):
            unit_id = str(unit["unit_id"])
            unit_dir = bundle.campaign_dir / "work" / unit_id
            completion_path = unit_dir / "simulation_complete.json"
            contract_path = unit_dir / "simulation_contract.json"
            if not completion_path.is_file() or not contract_path.is_file():
                raise ValueError(f"EAS attestation terminal artifacts are absent: {unit_id}")
            before = _unit_file_state(unit_dir)
            completion = _read_json(completion_path, label=f"EAS completion {unit_id}")
            stored_contract = completion.get("contract")
            if not isinstance(stored_contract, Mapping):
                raise TypeError(f"EAS completion contract is malformed: {unit_id}")
            disk_contract = _read_json(contract_path, label=f"EAS contract {unit_id}")
            stored_sha256 = _canonical_sha256(stored_contract)
            if (
                completion.get("contract_sha256") != stored_sha256
                or _canonical_sha256(disk_contract) != stored_sha256
            ):
                raise ValueError(f"EAS completion/contract self-binding differs: {unit_id}")

            artifacts = simulation._completion_artifacts(unit_dir)
            validated = simulation._validate_completion(artifacts, stored_contract)
            if validated.cache_hit is not True:
                raise RuntimeError("base EAS stored-contract validation was not cache-only")
            expected = _rebuild_current_expected_contract(bundle, unit)
            pair = _contract_attribute_pair(stored_contract)
            differences = _deep_difference_paths(stored_contract, expected)
            try:
                variant = _classify_eas_contract_variant(pair, differences)
            except ValueError as error:
                raise ValueError(f"{error}: {unit_id}") from error
            normalized = _normalized_eas_contract(stored_contract, expected)
            if _canonical_sha256(normalized) != _canonical_sha256(expected):
                raise ValueError(f"normalized EAS contract differs scientifically: {unit_id}")

            outputs = _attested_output_records(
                unit_dir=unit_dir,
                completion=completion,
                repo_root=bundle.repo_root,
            )
            attempt_audit = selected._validate_eas_attempt_provenance(
                bundle, unit, completion
            )
            after = _unit_file_state(unit_dir)
            if before != after:
                raise RuntimeError(f"EAS contract attestation changed a cache: {unit_id}")
            variant_counts[variant] += 1
            rows.append(
                {
                    "unit_id": unit_id,
                    "contract_variant": variant,
                    "completion": _file_record(completion_path, bundle.repo_root),
                    "contract_file": _file_record(contract_path, bundle.repo_root),
                    "stored_contract_sha256": stored_sha256,
                    "normalized_contract_sha256": _canonical_sha256(normalized),
                    "current_expected_contract_sha256": _canonical_sha256(expected),
                    "attribute_source_hashes": dict(
                        zip(EAS_ATTRIBUTE_SOURCE_PATHS, pair, strict=True)
                    ),
                    "contract_difference_paths": [
                        list(path) for path in sorted(differences)
                    ],
                    "outputs": outputs,
                    "outputs_canonical_sha256": _canonical_sha256(outputs),
                    "base_validation_against_stored_contract": True,
                    "current_contract_rebuilt_without_simulation": True,
                    "normalized_equality_except_allowed_attribute_hashes": True,
                    "attempt_provenance": attempt_audit,
                }
            )
    if len(rows) != 60 or variant_counts != EAS_EXPECTED_VARIANT_COUNTS:
        raise ValueError(
            "EAS contract variant counts differ: "
            f"observed {variant_counts}, expected {EAS_EXPECTED_VARIANT_COUNTS}"
        )
    return {
        "units": 60,
        "cells": 6,
        "replicates_per_cell": 10,
        "cache_validated_read_only": True,
        "base_validated_against_each_stored_contract": True,
        "current_contract_independently_rebuilt_for_each_unit": True,
        "normalized_differences_limited_to_two_attribute_hashes": True,
        "simulation_launched": False,
        "artifacts_changed": False,
        "variant_counts": variant_counts,
        "rows_canonical_sha256": _canonical_sha256(rows),
        "rows": rows,
    }


def build_eas_contract_attestation(
    *,
    repo_root: str | Path,
    slim_path: str | Path,
    integration_manifest_path: str | Path = selected.DEFAULT_INTEGRATION_RELATIVE_PATH,
) -> dict[str, Any]:
    root = Path(repo_root).resolve()
    if any("onedrive" in part.casefold() for part in root.parts):
        raise ValueError("EAS contract attestation repo_root must not be inside OneDrive")
    slim = _inside_repo(slim_path, root, label="SLiM binary")
    integration = _inside_repo(
        integration_manifest_path, root, label="selected integration manifest"
    )
    bundle = selected.load_selected_production_bundle(
        integration, repo_root=root, slim_path=slim
    )
    payload: dict[str, Any] = {
        "schema": EAS_CONTRACT_ATTESTATION_SCHEMA,
        "status": "passed",
        "selected_integration": {
            **_file_record(integration, root),
            "payload_sha256": bundle.manifest["payload_sha256"],
        },
        "execution_units": _file_record(bundle.execution_units_path, root),
        "attestation_sources": _canonical_source_records(
            root, POSTPRODUCTION_SOURCE_PATHS
        ),
        "canonical_attribute_sources": _validate_current_canonical_attributes(root),
        "allowed_attribute_hash_pairs": {
            "canonical": dict(
                zip(
                    EAS_ATTRIBUTE_SOURCE_PATHS,
                    EAS_CANONICAL_ATTRIBUTE_HASH_PAIR,
                    strict=True,
                )
            ),
            "transient_audit_text": dict(
                zip(
                    EAS_ATTRIBUTE_SOURCE_PATHS,
                    EAS_TRANSIENT_ATTRIBUTE_HASH_PAIR,
                    strict=True,
                )
            ),
        },
        "transient_attribute_evidence": _transient_attribute_evidence_records(root),
        "cache_audit": audit_eas_contract_attestation(bundle),
        "science_source_hash_changes_allowed": False,
        "only_normalized_fields": [list(path) for path in EAS_ATTRIBUTE_CONTRACT_PATHS],
        "mixed_attribute_pairs_allowed": False,
        "unknown_attribute_hashes_allowed": False,
        "completion_contracts_or_outputs_modified": False,
        "simulation_route_present": False,
    }
    payload["payload_sha256"] = _canonical_sha256(payload)
    return payload


def write_eas_contract_attestation(
    path: str | Path,
    payload: Mapping[str, Any],
    *,
    repo_root: str | Path,
) -> Path:
    _validate_hashed_payload(
        payload,
        schema=EAS_CONTRACT_ATTESTATION_SCHEMA,
        status="passed",
        label="EAS contract attestation",
    )
    root = Path(repo_root).resolve()
    destination = _inside_repo(path, root, label="EAS contract attestation")
    if destination.relative_to(root).as_posix() != DEFAULT_EAS_ATTESTATION_RELATIVE_PATH:
        raise ValueError("canonical EAS contract-attestation path differs")
    if any("onedrive" in part.casefold() for part in destination.parts):
        raise ValueError("EAS contract attestation must not be inside OneDrive")
    expected_records = _transient_attribute_evidence_records(root)
    if payload.get("transient_attribute_evidence") != expected_records:
        raise ValueError("EAS transient attribute evidence binding differs")
    for source_path, record in expected_records.items():
        fixed = EAS_TRANSIENT_ATTRIBUTE_EVIDENCE[source_path]
        content = base64.b64decode(str(fixed["content_base64"]), validate=True)
        archive = _inside_repo(root / record["path"], root, label="attribute archive")
        _immutable_bytes(archive, content, label="transient attribute evidence")
    content = (
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")
    return _immutable_bytes(destination, content, label="EAS contract attestation")


def load_eas_contract_attestation(
    path: str | Path,
    *,
    repo_root: str | Path,
    slim_path: str | Path,
    audit_current: bool = True,
) -> dict[str, Any]:
    root = Path(repo_root).resolve()
    attestation_path = _inside_repo(path, root, label="EAS contract attestation")
    if (
        attestation_path.relative_to(root).as_posix()
        != DEFAULT_EAS_ATTESTATION_RELATIVE_PATH
    ):
        raise ValueError("canonical EAS contract-attestation path differs")
    payload = _read_json(attestation_path, label="EAS contract attestation")
    _validate_hashed_payload(
        payload,
        schema=EAS_CONTRACT_ATTESTATION_SCHEMA,
        status="passed",
        label="EAS contract attestation",
    )
    expected_records = _transient_attribute_evidence_records(root)
    if payload.get("transient_attribute_evidence") != expected_records:
        raise ValueError("EAS transient attribute evidence binding differs")
    for source_path, record in expected_records.items():
        fixed = EAS_TRANSIENT_ATTRIBUTE_EVIDENCE[source_path]
        expected_content = base64.b64decode(
            str(fixed["content_base64"]), validate=True
        )
        archive = _inside_repo(root / record["path"], root, label="attribute archive")
        if (
            not archive.is_file()
            or archive.is_symlink()
            or archive.read_bytes() != expected_content
            or _file_record(archive, root) != {
                "path": record["path"],
                "size_bytes": record["size_bytes"],
                "sha256": record["sha256"],
            }
        ):
            raise ValueError("archived transient attribute evidence differs")
    if audit_current:
        integration = root / str(payload.get("selected_integration", {}).get("path", ""))
        current = build_eas_contract_attestation(
            repo_root=root,
            slim_path=slim_path,
            integration_manifest_path=integration,
        )
        if current != payload:
            raise ValueError("EAS contract attestation differs from current caches")
    return payload


def _attested_eas_cache_audit(
    payload: Mapping[str, Any], *, path: Path, repo_root: Path
) -> dict[str, Any]:
    audit = payload.get("cache_audit")
    if not isinstance(audit, Mapping):
        raise TypeError("EAS contract attestation cache audit is malformed")
    return {
        **copy.deepcopy(dict(audit)),
        "contract_attestation": {
            **_file_record(path, repo_root),
            "payload_sha256": payload["payload_sha256"],
        },
    }


@contextmanager
def patched_eas_contract_attestation_validation(
    payload: Mapping[str, Any],
    *,
    path: str | Path,
    repo_root: str | Path,
) -> Iterator[None]:
    """Narrowly replace only EAS cache validation for postproduction routes."""

    root = Path(repo_root).resolve()
    attestation_path = _inside_repo(path, root, label="EAS contract attestation")
    expected_integration = payload.get("selected_integration", {})
    expected_units = payload.get("execution_units", {})
    original = selected.validate_eas_cached_units_read_only

    def validate(bundle: selected.SelectedProductionBundle) -> dict[str, Any]:
        if (
            not isinstance(expected_integration, Mapping)
            or not isinstance(expected_units, Mapping)
            or _file_record(bundle.integration_manifest_path, root)
            != {
                key: expected_integration[key]
                for key in ("path", "size_bytes", "sha256")
            }
            or _file_record(bundle.execution_units_path, root) != dict(expected_units)
        ):
            raise ValueError("EAS attestation bundle binding differs")
        return _attested_eas_cache_audit(
            payload, path=attestation_path, repo_root=root
        )

    with _PATCH_LOCK:
        selected.validate_eas_cached_units_read_only = validate
        try:
            yield
        finally:
            selected.validate_eas_cached_units_read_only = original


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


def _validate_eas_preflight_binding(
    preflight: Mapping[str, Any],
    *,
    attestation: Mapping[str, Any],
    attestation_path: Path,
    repo_root: Path,
) -> None:
    cache_audit = preflight.get("eas", {}).get("cache_audit")
    expected = _attested_eas_cache_audit(
        attestation, path=attestation_path, repo_root=repo_root
    )
    if cache_audit != expected:
        raise ValueError("EAS selected-readiness attestation binding differs")


def build_postproduction_authorization(
    *,
    repo_root: str | Path,
    slim_path: str | Path,
    integration_manifest_path: str | Path = selected.DEFAULT_INTEGRATION_RELATIVE_PATH,
    normalization_path: str | Path = DEFAULT_NORMALIZATION_RELATIVE_PATH,
    eas_attestation_path: str | Path = DEFAULT_EAS_ATTESTATION_RELATIVE_PATH,
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
    eas_attestation_file = _inside_repo(
        eas_attestation_path, root, label="EAS contract attestation"
    )
    superseded_legacy = _load_superseded_legacy_authorization(root)

    evidence = _validate_normalization_evidence_for_bundle(
        repo_root=root,
        normalization_path=normalization,
        slim_path=slim,
    )
    eas_attestation = load_eas_contract_attestation(
        eas_attestation_file,
        repo_root=root,
        slim_path=slim,
        audit_current=True,
    )
    with (
        recovery_audit.patched_recovery_audit_aggregation(),
        patched_eas_contract_attestation_validation(
            eas_attestation,
            path=eas_attestation_file,
            repo_root=root,
        ),
    ):
        preflight = selected.build_aggregate_readiness(
            repo_root=root,
            eas_integration_manifest_path=integration,
            slim_path=slim,
            require_final=False,
            validate_eas_cache=True,
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
    _validate_eas_preflight_binding(
        preflight,
        attestation=eas_attestation,
        attestation_path=eas_attestation_file,
        repo_root=root,
    )
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
        "eas_contract_attestation": {
            **_file_record(eas_attestation_file, root),
            "payload_sha256": eas_attestation["payload_sha256"],
            "validated_units": 60,
            "variant_counts": dict(
                eas_attestation["cache_audit"]["variant_counts"]
            ),
        },
        "superseded_legacy_authorization": superseded_legacy,
        "han_normalization_evidence": evidence,
        "preflight_payload_sha256": preflight["payload_sha256"],
        "preflight_eas_cache_validation_deferred": False,
        "contexts": [
            "focused_han_direct_v3_recovery_audit.patched_recovery_audit_aggregation",
            "focused_selected_postproduction."
            "patched_eas_contract_attestation_validation",
        ],
        "actions": [
            "validate",
            "aggregate-validate",
            "aggregate",
            "decode",
            "analyze",
            "report",
        ],
        "simulation_route_present": False,
        "standard_selected_simulate_eas_untouched": True,
        "every_postproduction_route_inside_both_contexts": True,
        "eas_cache_validator_patched_only_inside_postproduction_context": True,
        "decode_contract": {
            "delegate": "focused_selection_campaign.main decode",
            "all_units": 720,
            "default_workers": DEFAULT_DECODE_WORKERS,
            "default_threads_per_decoder": DEFAULT_DECODE_THREADS,
            "default_decoder": DEFAULT_DECODER_RELATIVE_PATH,
            "raw_posterior_retained": True,
        },
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
        eas_attestation_path=(
            root
            / str(payload.get("eas_contract_attestation", {}).get("path", ""))
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
    root = Path(repo_root).resolve()
    attestation_path = root / str(
        authorization.get("eas_contract_attestation", {}).get("path", "")
    )
    attestation = load_eas_contract_attestation(
        attestation_path,
        repo_root=root,
        slim_path=slim_path,
        audit_current=False,
    )
    with (
        recovery_audit.patched_recovery_audit_aggregation(),
        patched_eas_contract_attestation_validation(
            attestation,
            path=attestation_path,
            repo_root=root,
        ),
    ):
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
    parser.add_argument("--eas-attestation", type=Path)
    parser.add_argument("--aggregate-readiness", type=Path)
    parser.add_argument("--empirical-tsv", type=Path)
    parser.add_argument(
        "--decoder-bin", type=Path, default=Path(DEFAULT_DECODER_RELATIVE_PATH)
    )
    parser.add_argument("--workers", type=int, default=DEFAULT_DECODE_WORKERS)
    parser.add_argument("--threads", type=int, default=DEFAULT_DECODE_THREADS)
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
    eas_attestation = (
        args.eas_attestation or root / DEFAULT_EAS_ATTESTATION_RELATIVE_PATH
    )
    readiness = (
        args.aggregate_readiness
        or root / selected.DEFAULT_AGGREGATE_READINESS_RELATIVE_PATH
    )

    if args.action == "attest-eas":
        payload = build_eas_contract_attestation(
            repo_root=root,
            slim_path=slim,
            integration_manifest_path=integration,
        )
        write_eas_contract_attestation(
            eas_attestation, payload, repo_root=root
        )
        print(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False))
        return 0

    if args.action == "supersede-legacy":
        payload = supersede_legacy_postproduction_authorization(repo_root=root)
        print(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False))
        return 0

    if args.action == "authorize":
        supersede_legacy_postproduction_authorization(repo_root=root)
        payload = build_postproduction_authorization(
            repo_root=root,
            slim_path=slim,
            integration_manifest_path=integration,
            normalization_path=normalization,
            eas_attestation_path=eas_attestation,
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
                        "eas_contract_attestation_payload_sha256": (
                            postproduction_authorization[
                                "eas_contract_attestation"
                            ]["payload_sha256"]
                        ),
                        "simulation_route_present": False,
                    },
                    indent=2,
                    sort_keys=True,
                )
            )
            return 0
        if args.action == "decode":
            if not 1 <= int(args.workers) <= 24:
                raise ValueError("postproduction decode workers must be between 1 and 24")
            if int(args.threads) < 1:
                raise ValueError("postproduction decode threads must be positive")
            decoder = args.decoder_bin
            if not decoder.is_absolute():
                decoder = root / decoder
            decoder = _inside_repo(decoder, root, label="Gamma-SMC decoder")
            if not decoder.is_file() or decoder.is_symlink():
                raise ValueError("Gamma-SMC decoder is absent or unsafe")
            result = campaign.main(
                [
                    "decode",
                    "--repo-root",
                    str(root),
                    "--campaign-dir",
                    str(root / selected.DEFAULT_CAMPAIGN_RELATIVE_PATH),
                    "--decoder-bin",
                    str(decoder),
                    "--workers",
                    str(int(args.workers)),
                    "--threads",
                    str(int(args.threads)),
                ]
            )
            if result != 0:
                return int(result)
            status_path = (
                root
                / selected.DEFAULT_CAMPAIGN_RELATIVE_PATH
                / "results"
                / "decode_status.tsv"
            )
            status = pd.read_csv(status_path, sep="\t")
            print(
                json.dumps(
                    {
                        "status": "decoded",
                        "units": len(status),
                        "status_counts": {
                            str(key): int(value)
                            for key, value in status.groupby("status").size().items()
                        },
                        "workers": int(args.workers),
                        "threads_per_decoder": int(args.threads),
                        "raw_posterior_retained": True,
                        "decoder": _file_record(decoder, root),
                        "postproduction_authorization_payload_sha256": (
                            postproduction_authorization["payload_sha256"]
                        ),
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
    "DEFAULT_EAS_ATTESTATION_RELATIVE_PATH",
    "DEFAULT_NORMALIZATION_RELATIVE_PATH",
    "EAS_CONTRACT_ATTESTATION_SCHEMA",
    "LEGACY_POSTPRODUCTION_AUTHORIZATION_SCHEMA",
    "MODULE_SOURCE_PATH",
    "POSTPRODUCTION_ACTIONS",
    "POSTPRODUCTION_AUTHORIZATION_SCHEMA",
    "WRAPPER_SOURCE_PATH",
    "authorized_postproduction_context",
    "audit_eas_contract_attestation",
    "build_eas_contract_attestation",
    "build_postproduction_authorization",
    "canonical_text_sha256",
    "load_eas_contract_attestation",
    "load_postproduction_authorization",
    "main",
    "patched_eas_contract_attestation_validation",
    "supersede_legacy_postproduction_authorization",
    "write_eas_contract_attestation",
    "write_postproduction_authorization",
]


if __name__ == "__main__":
    raise SystemExit(main())
