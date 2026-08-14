"""Attested, simulation-free aggregation routes for the focused campaign.

The completed EAS selected caches contain two already-attested source-contract
variants.  The standard decoded-output aggregator rebuilds the current
simulation contract and therefore cannot directly validate the 40 completions
whose only differences are the two archived ``.gitattributes`` hashes.  This
additive adapter keeps the immutable postproduction authorization untouched,
enters its authenticated EAS and Han contexts, and narrowly validates each EAS
completion against its attestation row during nested aggregation.

Only aggregate, analyze, and report are exposed.  There is no simulation or
decode route.
"""

from __future__ import annotations

import argparse
import copy
import json
import re
import threading
from collections.abc import Iterator, Mapping, Sequence
from contextlib import contextmanager
from pathlib import Path
from typing import Any

from . import focused_selected_postproduction as postproduction
from . import focused_selected_production as selected
from . import focused_selection_simulation as simulation

MODULE_SOURCE_PATH = "python/gamma_smc_aou/focused_selected_postdecode.py"
WRAPPER_SOURCE_PATH = "scripts/run_focused_selected_postdecode.py"
ADAPTER_SOURCE_PATHS = (MODULE_SOURCE_PATH, WRAPPER_SOURCE_PATH)
POSTDECODE_AUTHORIZATION_SCHEMA = (
    "gamma-smc.focused-selected-postdecode-authorization/v1"
)
DEFAULT_POSTDECODE_AUTHORIZATION_RELATIVE_PATH = (
    "focused_selection_EAS_sim/calibration/focused_selected_production/"
    "postdecode_authorization.json"
)
POSTDECODE_ROUTES = ("aggregate", "analyze", "report")
POSTDECODE_ACTIONS = ("authorize", *POSTDECODE_ROUTES)

_SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")
_PATCH_LOCK = threading.RLock()


def _strict_positive_integer(value: Any, *, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{label} must be a strict positive integer")
    return value


def _strict_sha256(value: Any, *, label: str) -> str:
    if not isinstance(value, str) or _SHA256_PATTERN.fullmatch(value) is None:
        raise ValueError(f"{label} SHA-256 is malformed")
    return value


def _validate_file_record(
    value: Any, *, expected_path: str, label: str
) -> dict[str, Any]:
    if not isinstance(value, Mapping) or set(value) != {
        "path",
        "sha256",
        "size_bytes",
    }:
        raise ValueError(f"{label} file record is malformed")
    if value.get("path") != expected_path:
        raise ValueError(f"{label} path differs")
    _strict_sha256(value.get("sha256"), label=label)
    _strict_positive_integer(value.get("size_bytes"), label=f"{label} size")
    return dict(value)


def _attested_eas_rows(payload: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    """Validate and index the exact 60 authenticated EAS cache rows."""

    audit = payload.get("cache_audit")
    if not isinstance(audit, Mapping):
        raise TypeError("EAS aggregate attestation cache audit is malformed")
    rows = audit.get("rows")
    if not isinstance(rows, list) or len(rows) != 60:
        raise ValueError("EAS aggregate attestation must contain exactly 60 rows")
    if (
        isinstance(audit.get("units"), bool)
        or audit.get("units") != 60
        or audit.get("cells") != 6
        or audit.get("replicates_per_cell") != 10
        or audit.get("cache_validated_read_only") is not True
        or audit.get("simulation_launched") is not False
        or audit.get("artifacts_changed") is not False
        or audit.get("base_validated_against_each_stored_contract") is not True
        or audit.get("current_contract_independently_rebuilt_for_each_unit") is not True
        or audit.get("normalized_differences_limited_to_two_attribute_hashes")
        is not True
        or audit.get("variant_counts") != postproduction.EAS_EXPECTED_VARIANT_COUNTS
        or audit.get("rows_canonical_sha256") != postproduction._canonical_sha256(rows)
    ):
        raise ValueError("EAS aggregate attestation cache-audit contract differs")

    expected_output_paths = {
        "tree": "simulation.trees",
        "pair_table": "sample_manifest.tsv",
        "overall_pairs": "pairs/overall.pairs.tsv",
        "truth_profiles": "truth_profiles.tsv.gz",
        "truth_class_summaries": "truth_class_summaries.tsv",
    }
    indexed: dict[str, dict[str, Any]] = {}
    observed_variants = {key: 0 for key in postproduction.EAS_EXPECTED_VARIANT_COUNTS}
    for raw_row in rows:
        if not isinstance(raw_row, Mapping):
            raise TypeError("EAS aggregate attestation row is malformed")
        row = dict(raw_row)
        unit_id = row.get("unit_id")
        if (
            not isinstance(unit_id, str)
            or not unit_id.startswith("eas_phlash_median__")
            or unit_id in indexed
        ):
            raise ValueError("EAS aggregate attestation unit IDs differ")
        unit_root = f"{selected.DEFAULT_CAMPAIGN_RELATIVE_PATH}/work/{unit_id}"
        _validate_file_record(
            row.get("completion"),
            expected_path=f"{unit_root}/simulation_complete.json",
            label=f"EAS {unit_id} completion",
        )
        _validate_file_record(
            row.get("contract_file"),
            expected_path=f"{unit_root}/simulation_contract.json",
            label=f"EAS {unit_id} contract",
        )
        outputs = row.get("outputs")
        if not isinstance(outputs, Mapping) or set(outputs) != set(
            expected_output_paths
        ):
            raise ValueError(f"EAS {unit_id} output records differ")
        for label, relative in expected_output_paths.items():
            _validate_file_record(
                outputs.get(label),
                expected_path=f"{unit_root}/{relative}",
                label=f"EAS {unit_id} {label}",
            )
        if row.get("outputs_canonical_sha256") != postproduction._canonical_sha256(
            outputs
        ):
            raise ValueError(f"EAS {unit_id} output-record digest differs")
        for key in (
            "stored_contract_sha256",
            "normalized_contract_sha256",
            "current_expected_contract_sha256",
        ):
            _strict_sha256(row.get(key), label=f"EAS {unit_id} {key}")
        variant = row.get("contract_variant")
        if variant not in observed_variants:
            raise ValueError(f"EAS {unit_id} contract variant is unknown")
        expected_pair = (
            postproduction.EAS_CANONICAL_ATTRIBUTE_HASH_PAIR
            if variant == "canonical"
            else postproduction.EAS_TRANSIENT_ATTRIBUTE_HASH_PAIR
        )
        expected_sources = dict(
            zip(postproduction.EAS_ATTRIBUTE_SOURCE_PATHS, expected_pair, strict=True)
        )
        expected_differences = (
            []
            if variant == "canonical"
            else [
                list(path)
                for path in sorted(postproduction.EAS_ATTRIBUTE_CONTRACT_PATHS)
            ]
        )
        if (
            row.get("attribute_source_hashes") != expected_sources
            or row.get("contract_difference_paths") != expected_differences
            or row.get("base_validation_against_stored_contract") is not True
            or row.get("current_contract_rebuilt_without_simulation") is not True
            or row.get("normalized_equality_except_allowed_attribute_hashes")
            is not True
        ):
            raise ValueError(f"EAS {unit_id} attested contract semantics differ")
        observed_variants[str(variant)] += 1
        indexed[unit_id] = row
    if observed_variants != postproduction.EAS_EXPECTED_VARIANT_COUNTS:
        raise ValueError("EAS aggregate attestation row variants differ")
    return indexed


def _is_eas_selected_contract(contract: Mapping[str, Any]) -> bool:
    unit = contract.get("unit")
    return isinstance(unit, Mapping) and (
        str(unit.get("simulation_class")) == "selected"
        and str(unit.get("demography_id")) == "eas_phlash_median"
    )


def _validate_attested_eas_completion(
    artifacts: simulation.SimulationArtifacts,
    current_contract: Mapping[str, Any],
    *,
    row: Mapping[str, Any],
    repo_root: Path,
    original_validator: Any,
) -> simulation.SimulationArtifacts:
    unit = current_contract.get("unit")
    if not isinstance(unit, Mapping):
        raise TypeError("current EAS aggregate unit contract is malformed")
    record = simulation._validate_unit(unit)
    unit_id = str(record["unit_id"])
    expected_unit_dir = (
        repo_root / selected.DEFAULT_CAMPAIGN_RELATIVE_PATH / "work" / unit_id
    )
    safe_unit_dir = postproduction._inside_repo(
        artifacts.unit_dir, repo_root, label=f"EAS aggregate unit {unit_id}"
    )
    safe_expected_unit_dir = postproduction._inside_repo(
        expected_unit_dir, repo_root, label=f"expected EAS aggregate unit {unit_id}"
    )
    if safe_unit_dir != safe_expected_unit_dir:
        raise ValueError(f"EAS aggregate unit directory differs: {unit_id}")
    if postproduction._canonical_sha256(current_contract) != row.get(
        "current_expected_contract_sha256"
    ):
        raise ValueError(f"current EAS aggregate contract differs: {unit_id}")

    completion = postproduction._read_json(
        artifacts.completion_path, label=f"EAS aggregate completion {unit_id}"
    )
    stored_contract = completion.get("contract")
    if not isinstance(stored_contract, Mapping):
        raise TypeError(f"stored EAS aggregate contract is malformed: {unit_id}")
    stored_sha256 = postproduction._canonical_sha256(stored_contract)
    if (
        completion.get("contract_sha256") != stored_sha256
        or stored_sha256 != row.get("stored_contract_sha256")
        or postproduction._file_record(artifacts.completion_path, repo_root)
        != row.get("completion")
    ):
        raise ValueError(f"stored EAS aggregate completion differs: {unit_id}")
    contract_path = artifacts.unit_dir / "simulation_contract.json"
    disk_contract = postproduction._read_json(
        contract_path, label=f"EAS aggregate contract file {unit_id}"
    )
    if postproduction._canonical_sha256(
        disk_contract
    ) != stored_sha256 or postproduction._file_record(
        contract_path, repo_root
    ) != row.get("contract_file"):
        raise ValueError(f"disk EAS aggregate contract differs: {unit_id}")

    pair = postproduction._contract_attribute_pair(stored_contract)
    differences = postproduction._deep_difference_paths(
        stored_contract, current_contract
    )
    variant = postproduction._classify_eas_contract_variant(pair, differences)
    normalized = postproduction._normalized_eas_contract(
        stored_contract, current_contract
    )
    if (
        variant != row.get("contract_variant")
        or dict(zip(postproduction.EAS_ATTRIBUTE_SOURCE_PATHS, pair, strict=True))
        != row.get("attribute_source_hashes")
        or [list(path) for path in sorted(differences)]
        != row.get("contract_difference_paths")
        or postproduction._canonical_sha256(normalized)
        != row.get("normalized_contract_sha256")
        or postproduction._canonical_sha256(normalized)
        != postproduction._canonical_sha256(current_contract)
    ):
        raise ValueError(f"normalized EAS aggregate contract differs: {unit_id}")
    outputs = postproduction._attested_output_records(
        unit_dir=artifacts.unit_dir,
        completion=completion,
        repo_root=repo_root,
    )
    if outputs != row.get("outputs") or postproduction._canonical_sha256(
        outputs
    ) != row.get("outputs_canonical_sha256"):
        raise ValueError(f"attested EAS aggregate outputs differ: {unit_id}")

    validated = original_validator(artifacts, stored_contract)
    if validated.cache_hit is not True:
        raise RuntimeError(f"EAS aggregate validation was not cache-only: {unit_id}")
    return validated


@contextmanager
def patched_attested_eas_aggregate_validation(
    payload: Mapping[str, Any],
    *,
    path: str | Path,
    repo_root: str | Path,
    require_complete_coverage: bool,
) -> Iterator[dict[str, Any]]:
    """Patch only EAS per-unit completion validation for aggregation."""

    root = Path(repo_root).resolve()
    attestation_path = postproduction._inside_repo(
        path, root, label="EAS aggregate attestation"
    )
    if (
        attestation_path.relative_to(root).as_posix()
        != postproduction.DEFAULT_EAS_ATTESTATION_RELATIVE_PATH
    ):
        raise ValueError("canonical EAS aggregate-attestation path differs")
    postproduction._validate_hashed_payload(
        payload,
        schema=postproduction.EAS_CONTRACT_ATTESTATION_SCHEMA,
        status="passed",
        label="EAS aggregate attestation",
    )
    on_disk = postproduction._read_json(
        attestation_path, label="EAS aggregate attestation"
    )
    if on_disk != dict(payload):
        raise ValueError("EAS aggregate attestation file differs from loaded payload")
    rows = _attested_eas_rows(payload)
    coverage: dict[str, Any] = {
        "required": bool(require_complete_coverage),
        "validated_unit_ids": [],
    }

    with _PATCH_LOCK:
        original = simulation._validate_completion

        def validate(
            artifacts: simulation.SimulationArtifacts,
            contract: Mapping[str, Any],
        ) -> simulation.SimulationArtifacts:
            if not _is_eas_selected_contract(contract):
                return original(artifacts, contract)
            unit = contract.get("unit")
            if not isinstance(unit, Mapping):  # pragma: no cover - guarded above
                raise TypeError("current EAS aggregate unit contract is malformed")
            unit_id = str(unit.get("unit_id", ""))
            row = rows.get(unit_id)
            if row is None:
                raise ValueError(f"EAS aggregate unit lacks attestation: {unit_id}")
            validated = _validate_attested_eas_completion(
                artifacts,
                contract,
                row=row,
                repo_root=root,
                original_validator=original,
            )
            coverage["validated_unit_ids"].append(unit_id)
            return validated

        simulation._validate_completion = validate
        succeeded = False
        try:
            yield coverage
            succeeded = True
        finally:
            simulation._validate_completion = original
    if succeeded and require_complete_coverage:
        observed = list(coverage["validated_unit_ids"])
        if len(observed) != 60 or set(observed) != set(rows):
            raise ValueError(
                "postdecode aggregation did not validate every EAS attested unit "
                "exactly once"
            )


def _route_binding(
    authorization: Mapping[str, Any],
    readiness: Mapping[str, Any],
    attestation: Mapping[str, Any],
    *,
    attestation_path: Path,
    repo_root: Path,
) -> dict[str, Any]:
    """Require the immutable authorization, readiness, and attestation to agree."""

    attestation_record = {
        **postproduction._file_record(attestation_path, repo_root),
        "payload_sha256": attestation["payload_sha256"],
    }
    authorization_attestation_record = {
        **attestation_record,
        "validated_units": 60,
        "variant_counts": dict(postproduction.EAS_EXPECTED_VARIANT_COUNTS),
    }
    readiness_eas = readiness.get("eas")
    if not isinstance(readiness_eas, Mapping):
        raise TypeError("aggregate readiness EAS binding is malformed")
    readiness_audit = readiness_eas.get("cache_audit")
    if not isinstance(readiness_audit, Mapping):
        raise TypeError("aggregate readiness EAS cache audit is malformed")
    readiness_audit_without_record = copy.deepcopy(dict(readiness_audit))
    readiness_attestation = readiness_audit_without_record.pop(
        "contract_attestation", None
    )
    authorization_integration = authorization.get("selected_integration")
    readiness_integration = readiness_eas.get("integration_manifest")
    if not isinstance(authorization_integration, Mapping) or not isinstance(
        readiness_integration, Mapping
    ):
        raise TypeError("selected integration route binding is malformed")
    expected_readiness_integration = {
        **dict(readiness_integration),
        "payload_sha256": readiness_eas.get("integration_payload_sha256"),
        "strict_loaded": True,
    }
    if (
        authorization.get("preflight_payload_sha256") != readiness.get("payload_sha256")
        or authorization.get("eas_contract_attestation")
        != authorization_attestation_record
        or readiness_attestation != attestation_record
        or readiness_audit_without_record != attestation.get("cache_audit")
        or dict(authorization_integration) != expected_readiness_integration
        or authorization.get("every_postproduction_route_inside_both_contexts")
        is not True
        or authorization.get(
            "eas_cache_validator_patched_only_inside_postproduction_context"
        )
        is not True
        or not set(POSTDECODE_ROUTES).issubset(set(authorization.get("actions", ())))
    ):
        raise ValueError("postdecode authorization/readiness route binding differs")
    return {
        "postproduction_authorization_payload_sha256": authorization["payload_sha256"],
        "aggregate_readiness_payload_sha256": readiness["payload_sha256"],
        "eas_contract_attestation_payload_sha256": attestation["payload_sha256"],
        "eas_attested_units": len(_attested_eas_rows(attestation)),
        "han_nan_normalization_context": True,
        "simulation_route_present": False,
        "decode_route_present": False,
    }


def build_postdecode_authorization(
    *,
    repo_root: str | Path,
    upstream_authorization: Mapping[str, Any],
    readiness: Mapping[str, Any],
    attestation: Mapping[str, Any],
    upstream_authorization_path: str | Path,
    readiness_path: str | Path,
    attestation_path: str | Path,
) -> dict[str, Any]:
    """Build a source- and upstream-bound authorization without writing it."""

    root = Path(repo_root).resolve()
    upstream_path = postproduction._inside_repo(
        upstream_authorization_path,
        root,
        label="upstream postproduction authorization",
    )
    ready_path = postproduction._inside_repo(
        readiness_path, root, label="aggregate readiness"
    )
    eas_path = postproduction._inside_repo(
        attestation_path, root, label="EAS contract attestation"
    )
    route = _route_binding(
        upstream_authorization,
        readiness,
        attestation,
        attestation_path=eas_path,
        repo_root=root,
    )
    payload: dict[str, Any] = {
        "schema": POSTDECODE_AUTHORIZATION_SCHEMA,
        "status": "authorized",
        "adapter_sources": postproduction._canonical_source_records(
            root, ADAPTER_SOURCE_PATHS
        ),
        "upstream_artifacts": {
            "postproduction_authorization": {
                **postproduction._file_record(upstream_path, root),
                "payload_sha256": upstream_authorization["payload_sha256"],
            },
            "aggregate_readiness": {
                **postproduction._file_record(ready_path, root),
                "payload_sha256": readiness["payload_sha256"],
            },
            "eas_contract_attestation": {
                **postproduction._file_record(eas_path, root),
                "payload_sha256": attestation["payload_sha256"],
            },
        },
        "route_binding": route,
        "authorization_actions": ["authorize"],
        "postdecode_routes": list(POSTDECODE_ROUTES),
        "eas_per_unit_validation": {
            "patch_target": "focused_selection_simulation._validate_completion",
            "scope": "eas_phlash_median selected units only",
            "current_contract_bound_to_exact_attestation_row": True,
            "stored_contract_and_outputs_bound_to_exact_attestation_row": True,
            "base_validator_called_against_stored_contract": True,
            "non_eas_delegated_unchanged": True,
            "aggregate_and_analyze_require_exactly_60_unique_units": True,
            "report_requires_per_unit_coverage": False,
        },
        "outer_contexts": [
            "focused_han_direct_v3_recovery_audit.patched_recovery_audit_aggregation",
            "focused_selected_postproduction.patched_eas_contract_attestation_validation",
        ],
        "global_patches_restored_on_all_exits": True,
        "simulation_route_present": False,
        "decode_route_present": False,
        "recovery_route_present": False,
        "only_preauthorization_write": DEFAULT_POSTDECODE_AUTHORIZATION_RELATIVE_PATH,
    }
    payload["payload_sha256"] = postproduction._canonical_sha256(payload)
    return payload


def write_postdecode_authorization(
    path: str | Path,
    payload: Mapping[str, Any],
    *,
    repo_root: str | Path,
) -> Path:
    """Publish exact postdecode authorization bytes once."""

    postproduction._validate_hashed_payload(
        payload,
        schema=POSTDECODE_AUTHORIZATION_SCHEMA,
        status="authorized",
        label="postdecode authorization",
    )
    root = Path(repo_root).resolve()
    destination = postproduction._inside_repo(
        path, root, label="postdecode authorization"
    )
    if (
        destination.relative_to(root).as_posix()
        != DEFAULT_POSTDECODE_AUTHORIZATION_RELATIVE_PATH
    ):
        raise ValueError("canonical postdecode-authorization path differs")
    if any("onedrive" in part.casefold() for part in destination.parts):
        raise ValueError("postdecode authorization must not be inside OneDrive")
    content = (
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")
    return postproduction._immutable_bytes(
        destination, content, label="postdecode authorization"
    )


def load_postdecode_authorization(
    path: str | Path,
    *,
    repo_root: str | Path,
    expected: Mapping[str, Any],
) -> dict[str, Any]:
    """Load the canonical authorization and require an exact current rebuild."""

    root = Path(repo_root).resolve()
    authorization_path = postproduction._inside_repo(
        path, root, label="postdecode authorization"
    )
    if (
        authorization_path.relative_to(root).as_posix()
        != DEFAULT_POSTDECODE_AUTHORIZATION_RELATIVE_PATH
    ):
        raise ValueError("canonical postdecode-authorization path differs")
    payload = postproduction._read_json(
        authorization_path, label="postdecode authorization"
    )
    postproduction._validate_hashed_payload(
        payload,
        schema=POSTDECODE_AUTHORIZATION_SCHEMA,
        status="authorized",
        label="postdecode authorization",
    )
    if payload != dict(expected):
        raise ValueError(
            "postdecode authorization differs from current sources/evidence"
        )
    return payload


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=POSTDECODE_ACTIONS)
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument(
        "--slim-bin", type=Path, default=Path(".native-stdpopsim/bin/slim")
    )
    parser.add_argument("--empirical-tsv", type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    root = args.repo_root.resolve()
    slim = args.slim_bin
    if not slim.is_absolute():
        slim = root / slim
    authorization_path = root / postproduction.DEFAULT_AUTHORIZATION_RELATIVE_PATH
    readiness_path = root / selected.DEFAULT_AGGREGATE_READINESS_RELATIVE_PATH
    postdecode_authorization_path = (
        root / DEFAULT_POSTDECODE_AUTHORIZATION_RELATIVE_PATH
    )

    with postproduction.authorized_postproduction_context(
        repo_root=root,
        slim_path=slim,
        authorization_path=authorization_path,
    ) as authorization:
        attestation_path = root / str(
            authorization.get("eas_contract_attestation", {}).get("path", "")
        )
        attestation = postproduction.load_eas_contract_attestation(
            attestation_path,
            repo_root=root,
            slim_path=slim,
            audit_current=False,
        )
        readiness = selected.load_aggregate_readiness(
            readiness_path, repo_root=root, slim_path=slim
        )
        route = _route_binding(
            authorization,
            readiness,
            attestation,
            attestation_path=attestation_path,
            repo_root=root,
        )
        expected_postdecode_authorization = build_postdecode_authorization(
            repo_root=root,
            upstream_authorization=authorization,
            readiness=readiness,
            attestation=attestation,
            upstream_authorization_path=authorization_path,
            readiness_path=readiness_path,
            attestation_path=attestation_path,
        )
        if args.action == "authorize":
            written = write_postdecode_authorization(
                postdecode_authorization_path,
                expected_postdecode_authorization,
                repo_root=root,
            )
            print(
                json.dumps(
                    {
                        "status": "authorized",
                        "path": written.relative_to(root).as_posix(),
                        "payload_sha256": expected_postdecode_authorization[
                            "payload_sha256"
                        ],
                        "simulation_route_present": False,
                        "decode_route_present": False,
                    },
                    indent=2,
                    sort_keys=True,
                    allow_nan=False,
                )
            )
            return 0
        postdecode_authorization = load_postdecode_authorization(
            postdecode_authorization_path,
            repo_root=root,
            expected=expected_postdecode_authorization,
        )
        with patched_attested_eas_aggregate_validation(
            attestation,
            path=attestation_path,
            repo_root=root,
            require_complete_coverage=args.action in {"aggregate", "analyze"},
        ) as coverage:
            if args.action == "aggregate":
                outputs = selected.aggregate_with_readiness(
                    repo_root=root,
                    slim_path=slim,
                    readiness_path=readiness_path,
                )
                payload: dict[str, Any] = {
                    "status": "aggregated",
                    "outputs": [
                        path.resolve().relative_to(root).as_posix() for path in outputs
                    ],
                }
            elif args.action == "analyze":
                payload = selected.analyze_with_readiness(
                    repo_root=root,
                    slim_path=slim,
                    readiness_path=readiness_path,
                    empirical_path=args.empirical_tsv,
                )
            elif args.action == "report":
                payload = selected.report_with_readiness(
                    repo_root=root,
                    slim_path=slim,
                    readiness_path=readiness_path,
                )
            else:  # pragma: no cover - argparse constrains the action
                raise AssertionError(f"unhandled postdecode action: {args.action}")
    result = {
        **payload,
        "postdecode_route_binding": route,
        "postdecode_authorization": {
            **postproduction._file_record(postdecode_authorization_path, root),
            "payload_sha256": postdecode_authorization["payload_sha256"],
        },
        "eas_aggregate_validation": {
            "coverage_required": coverage["required"],
            "validated_units": len(coverage["validated_unit_ids"]),
            "validated_unit_ids_sha256": postproduction._canonical_sha256(
                coverage["validated_unit_ids"]
            ),
        },
    }
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


__all__ = [
    "ADAPTER_SOURCE_PATHS",
    "DEFAULT_POSTDECODE_AUTHORIZATION_RELATIVE_PATH",
    "MODULE_SOURCE_PATH",
    "POSTDECODE_ACTIONS",
    "POSTDECODE_AUTHORIZATION_SCHEMA",
    "POSTDECODE_ROUTES",
    "WRAPPER_SOURCE_PATH",
    "build_postdecode_authorization",
    "load_postdecode_authorization",
    "main",
    "patched_attested_eas_aggregate_validation",
    "write_postdecode_authorization",
]
