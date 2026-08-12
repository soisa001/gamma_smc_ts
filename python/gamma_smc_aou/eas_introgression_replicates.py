"""Replicate campaign orchestration for the EAS introgressed-sweep study.

The single-trajectory study remains the source of the biological model and
the simulation/decode implementation.  This module only expands each of its
27 introgression cells into stable, independently seeded execution units and
places them below an isolated campaign directory.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import io
import json
import os
import re
import subprocess
import sys
import time
import uuid
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

import pandas as pd

from . import eas_sweep_study as study
from .eas_sweep_models import DOMINANCE_COEFFICIENT

REPLICATE_SCHEMA = "gamma-smc.eas-introgression-replicates/v1"
DEFAULT_REPLICATES = 10
DEFAULT_CAMPAIGN = "replicate_study_n10"
REQUIRED_SLIM_VERSION = "4.2.2"
BIOLOGICAL_SPECIFICATION_COUNT = 27
PHASE_INVOCATION_SCHEMA = "gamma-smc.eas-introgression-phase-invocation/v1"
PHASE_INVOCATION_DIRECTORY = "phase_invocations"
_PHASE_EVENT_RESERVED_FIELDS = frozenset(
    {
        "schema",
        "invocation_id",
        "event",
        "recorded_at_utc",
        "process_id",
        "elapsed_seconds",
        "tool_provenance",
        "error",
    }
)
_PLAN_CREATION_LOCK = ".study_plan_creation.lock"
_PLAN_RUNTIME_FIELDS = (
    "sample_diploids",
    "sequence_length_bp",
    "focal_position_bp",
    "mutation_rate",
    "recombination_rate",
    "generation_time_years",
    "output_stride_bp",
    "cache_size_bp",
    "slim_scaling_factor",
    "slim_burn_in",
    "minimum_homozygous_diploids",
    "base_seed",
    "campaign_subdirectory",
)


class ExternalProcessInterruption(RuntimeError):
    """A prior phase ended outside the orchestrator's terminal-event path."""


def _canonical_sha256(value: Any) -> str:
    payload = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _validate_replicate_count(replicates: int) -> int:
    if isinstance(replicates, bool) or int(replicates) != replicates:
        raise ValueError("replicates must be a positive integer")
    replicate_count = int(replicates)
    if replicate_count < 1:
        raise ValueError("replicates must be a positive integer")
    return replicate_count


def expand_introgression_specifications(
    base_specifications: Sequence[Mapping[str, Any]],
    replicates: int = DEFAULT_REPLICATES,
) -> list[dict[str, Any]]:
    """Expand biological cells in replicate-major order with stable unit IDs.

    ``replicate_id`` is globally unique and deliberately equals the expanded
    ``specification_id``.  The core seed derivation therefore implements
    ``sha256(base_seed:replicate_id:external_draw)`` without a second seed
    namespace.  No total-replicate-count field is stored on an execution unit,
    so increasing a campaign from *n* to *m* leaves the first *n* contracts
    byte-for-byte compatible.
    """
    replicate_count = _validate_replicate_count(replicates)
    base_records = [copy.deepcopy(dict(record)) for record in base_specifications]
    if len(base_records) != BIOLOGICAL_SPECIFICATION_COUNT:
        raise ValueError(
            "the introgression replicate campaign requires exactly "
            f"{BIOLOGICAL_SPECIFICATION_COUNT} biological specifications; "
            f"found {len(base_records)}"
        )
    base_ids = [str(record["specification_id"]) for record in base_records]
    if len(set(base_ids)) != len(base_ids):
        raise ValueError("base introgression specification IDs must be unique")
    if any(record.get("study_type") != "introgression" for record in base_records):
        raise ValueError("all base specifications must be introgression models")

    expanded: list[dict[str, Any]] = []
    for replicate_index in range(1, replicate_count + 1):
        suffix = f"rep{replicate_index:03d}"
        for base_record in base_records:
            base_id = str(base_record["specification_id"])
            replicate_id = f"{base_id}__{suffix}"
            record = copy.deepcopy(base_record)
            record.update(
                {
                    "base_specification_id": base_id,
                    "replicate_index": replicate_index,
                    "replicate_id": replicate_id,
                    "specification_id": replicate_id,
                }
            )
            expanded.append(record)
    return expanded


def build_introgression_replicates(
    runtime: study.StudyRuntime,
    replicates: int = DEFAULT_REPLICATES,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Build the fixed biological grid and its independent replicate slots."""
    runtime.validate()
    base = study.build_specifications(runtime, ("introgression",))["introgression"]
    expanded = expand_introgression_specifications(base, replicates)
    return base, expanded


def _compact_specification_frame(
    specifications: Sequence[Mapping[str, Any]],
) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    leading = (
        "base_specification_id",
        "replicate_index",
        "replicate_id",
        "specification_id",
        "study_type",
        "scenario",
        "selection_coefficient",
        "target_allele_frequency",
        "lower_allele_frequency",
        "upper_allele_frequency",
    )
    for specification in specifications:
        scalar = {
            key: value
            for key, value in specification.items()
            if not isinstance(value, (dict, list, tuple, set))
        }
        records.append(scalar)
    frame = pd.DataFrame(records)
    ordered = [column for column in leading if column in frame.columns]
    ordered.extend(column for column in frame.columns if column not in ordered)
    return frame.loc[:, ordered]


def _campaign_contract(
    runtime: study.StudyRuntime,
    base_specifications: Sequence[Mapping[str, Any]],
    replicates: int,
) -> dict[str, Any]:
    targets = sorted(
        {float(record["target_allele_frequency"]) for record in base_specifications}
    )
    selections = sorted(
        {float(record["selection_coefficient"]) for record in base_specifications},
        reverse=True,
    )
    return {
        "scenario": "introgression",
        "campaign_subdirectory": runtime.campaign_subdirectory,
        "biological_specification_count": BIOLOGICAL_SPECIFICATION_COUNT,
        "replicates_per_biological_specification": replicates,
        "target_execution_slots": BIOLOGICAL_SPECIFICATION_COUNT * replicates,
        "selection_coefficients": selections,
        "han_final_frequency_floors": targets,
        "sequence_length_bp": runtime.sequence_length_bp,
        "focal_position_bp": runtime.focal_position_bp,
        "sample_diploids": runtime.sample_diploids,
        "acceptance": {
            "han_population_frequency": {
                "minimum_inclusive": "cell target_allele_frequency",
                "maximum_inclusive": float(
                    max(
                        record["upper_allele_frequency"]
                        for record in base_specifications
                    )
                ),
                "conditioned_during_slim": True,
            },
            "sample_frequency": {
                "minimum_inclusive": "cell target_allele_frequency",
                "maximum_inclusive": float(
                    max(
                        record["upper_allele_frequency"]
                        for record in base_specifications
                    )
                ),
                "validated_after_sampling": True,
            },
            "genotype_classes": {
                "hom_ref_minimum_diploids": runtime.minimum_homozygous_diploids,
                "heterozygous_minimum_diploids": 1,
                "hom_alt_minimum_diploids": runtime.minimum_homozygous_diploids,
            },
            "archaic_specific_no_ils": {
                "mutation_declared_origin_population": "Neanderthal",
                "origin_after_human_neanderthal_split": True,
                "origin_before_introgression_pulse": True,
                "source_frequency_minimum_inclusive": float(
                    min(
                        record["source_minimum_frequency"]
                        for record in base_specifications
                    )
                ),
                "han_entry_route": "introgression_pulse_only",
                "tree_node_population_is_origin_evidence": False,
                "origin_evidence": "serialized_forward_event_contract",
            },
        },
        "selection": {
            "mode": "standing_variation_at_introgression_onset",
            "dominance_coefficient": DOMINANCE_COEFFICIENT,
        },
    }


def _portable_base_specifications(
    base: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    return [
        {
            key: value
            for key, value in record.items()
            if key not in {"model_record", "origin_record"}
        }
        for record in base
    ]


def _specification_table_sha256(frame: pd.DataFrame) -> str:
    """Hash the exact UTF-8 TSV representation used by ``_atomic_frame``."""
    payload = frame.to_csv(sep="\t", index=False).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _campaign_design(
    runtime: study.StudyRuntime,
    base: Sequence[Mapping[str, Any]],
    replicates: int,
    *,
    specification_table_sha256: str,
) -> dict[str, Any]:
    base_portable = _portable_base_specifications(base)
    contract = _campaign_contract(runtime, base, replicates)
    return {
        "schema": REPLICATE_SCHEMA,
        "phase": "plan",
        "scenario": "introgression",
        "campaign_contract": contract,
        "campaign_contract_sha256": _canonical_sha256(contract),
        "runtime": study._portable_runtime_record(runtime),  # noqa: SLF001
        "rng": {
            "base_seed": runtime.base_seed,
            "replicate_id_semantics": (
                "globally unique expanded specification ID: <base>__repNNN"
            ),
            "external_draw_index_origin": 0,
            "derivation": "sha256(base_seed:replicate_id:external_draw)",
            "integer_mapping": "1 + first_8_digest_bytes mod (2^31 - 2)",
        },
        "cache_semantics": {
            "unit_identity": "replicate_id",
            "replicate_count_stored_on_unit": False,
            "extension_rule": (
                "increasing target slots preserves existing repNNN IDs and work caches"
            ),
        },
        "hashes": {
            "base_specifications_sha256": _canonical_sha256(base_portable),
            "specifications_tsv_sha256": specification_table_sha256,
            "orchestrator_source_sha256": study._sha256(__file__),  # noqa: SLF001
            "core_orchestrator_source_sha256": study._sha256(  # noqa: SLF001
                Path(study.__file__)
            ),
        },
        "software": study._software_record(),  # noqa: SLF001
        "fixed_model_contract": {
            "demography": "stdpopsim AncientEurasia_9K19",
            "mutation_mode": "single archaic-specific de novo origin",
            "selection_onset": "introgression onset; standing variation thereafter",
            "neutral_null_simulations": False,
        },
        "slim_compatibility": {
            "required_version": REQUIRED_SLIM_VERSION,
            "slim_5_2": (
                "rejected: incompatible with the stdpopsim 0.3 generated model"
            ),
        },
    }


def _require_equal_plan_field(
    recorded: Mapping[str, Any], expected: Mapping[str, Any], field: str
) -> None:
    if _canonical_sha256(recorded.get(field)) != _canonical_sha256(expected[field]):
        raise ValueError(
            f"existing replicate campaign plan is incompatible in {field!r}; "
            "use its original parameters or choose a new --campaign directory"
        )


def _validate_existing_campaign_plan(
    runtime: study.StudyRuntime,
    base: Sequence[Mapping[str, Any]],
    expanded: Sequence[Mapping[str, Any]],
    replicates: int,
    campaign_dir: Path,
) -> None:
    """Validate an immutable plan without comparing drift-prone source hashes."""
    table_path = campaign_dir / "specifications.tsv"
    design_path = campaign_dir / "study_design.json"
    if not table_path.is_file() or not design_path.is_file():
        present = [path.name for path in (table_path, design_path) if path.exists()]
        raise ValueError(
            "replicate campaign plan is incomplete; expected both "
            f"specifications.tsv and study_design.json, found {present or 'neither'}"
        )
    try:
        recorded = json.loads(design_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(
            f"replicate campaign design is unreadable: {design_path}"
        ) from error

    expected_frame = _compact_specification_frame(expanded)
    expected_table_sha256 = _specification_table_sha256(expected_frame)
    expected = _campaign_design(
        runtime,
        base,
        replicates,
        specification_table_sha256=expected_table_sha256,
    )
    for field in (
        "schema",
        "phase",
        "scenario",
        "campaign_contract",
        "campaign_contract_sha256",
        "rng",
        "cache_semantics",
        "fixed_model_contract",
        "slim_compatibility",
    ):
        _require_equal_plan_field(recorded, expected, field)

    recorded_runtime = recorded.get("runtime")
    if not isinstance(recorded_runtime, dict):
        raise ValueError("replicate campaign design has no valid runtime record")
    expected_runtime = expected["runtime"]
    for field in _PLAN_RUNTIME_FIELDS:
        recorded_value = recorded_runtime.get(field)
        expected_value = expected_runtime[field]
        if (
            isinstance(recorded_value, (int, float))
            and not isinstance(recorded_value, bool)
            and isinstance(expected_value, (int, float))
            and not isinstance(expected_value, bool)
        ):
            compatible = float(recorded_value) == float(expected_value)
        else:
            compatible = recorded_value == expected_value
        if not compatible:
            raise ValueError(
                "existing replicate campaign plan is incompatible in runtime field "
                f"{field!r}; choose a new --campaign directory"
            )

    hashes = recorded.get("hashes")
    if not isinstance(hashes, dict):
        raise ValueError("replicate campaign design has no valid hashes record")
    recorded_table_sha256 = hashes.get("specifications_tsv_sha256")
    actual_table_sha256 = study._sha256(table_path)  # noqa: SLF001
    if recorded_table_sha256 != actual_table_sha256:
        raise ValueError(
            "replicate campaign specifications.tsv checksum does not match "
            "the immutable study design"
        )
    try:
        actual_frame = pd.read_csv(table_path, sep="\t")
        expected_roundtrip = pd.read_csv(
            io.StringIO(expected_frame.to_csv(sep="\t", index=False)), sep="\t"
        )
        pd.testing.assert_frame_equal(
            actual_frame,
            expected_roundtrip,
            check_dtype=False,
            check_exact=False,
            rtol=1e-12,
            atol=1e-12,
        )
    except (OSError, ValueError, pd.errors.ParserError, AssertionError) as error:
        raise ValueError(
            "existing replicate campaign specifications are incompatible with "
            "the requested biological grid or replicate count"
        ) from error
    expected_base_sha256 = expected["hashes"]["base_specifications_sha256"]
    if hashes.get("base_specifications_sha256") != expected_base_sha256:
        raise ValueError(
            "existing replicate campaign base-specification checksum is incompatible"
        )


def write_introgression_replicate_plan(
    runtime: study.StudyRuntime,
    replicates: int = DEFAULT_REPLICATES,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Create a plan once, then checksum-validate it without rewriting it."""
    replicate_count = _validate_replicate_count(replicates)
    if runtime.campaign_subdirectory is None:
        raise ValueError("replicate campaigns require campaign_subdirectory")
    base, expanded = build_introgression_replicates(runtime, replicate_count)
    campaign_dir = study._runtime_study_dir(runtime, "introgression")  # noqa: SLF001
    campaign_dir.mkdir(parents=True, exist_ok=True)
    table_path = campaign_dir / "specifications.tsv"
    design_path = campaign_dir / "study_design.json"

    if table_path.exists() or design_path.exists():
        _validate_existing_campaign_plan(
            runtime, base, expanded, replicate_count, campaign_dir
        )
        return base, expanded

    lock_path = campaign_dir / _PLAN_CREATION_LOCK
    try:
        descriptor = os.open(
            lock_path,
            os.O_CREAT | os.O_EXCL | os.O_WRONLY,
            0o644,
        )
    except FileExistsError as error:
        raise RuntimeError(
            "replicate campaign plan creation is already in progress or its lock "
            f"is stale: {lock_path}"
        ) from error
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(
                {
                    "pid": os.getpid(),
                    "created_at_utc": datetime.now(timezone.utc).isoformat(),
                },
                handle,
            )
            handle.write("\n")
        frame = _compact_specification_frame(expanded)
        expected_table_sha256 = _specification_table_sha256(frame)
        written_table = study._atomic_frame(table_path, frame)  # noqa: SLF001
        if study._sha256(written_table) != expected_table_sha256:  # noqa: SLF001
            raise RuntimeError("new replicate specification table failed checksum")
        design = _campaign_design(
            runtime,
            base,
            replicate_count,
            specification_table_sha256=expected_table_sha256,
        )
        study._atomic_json(design_path, design)  # noqa: SLF001
        _validate_existing_campaign_plan(
            runtime, base, expanded, replicate_count, campaign_dir
        )
    finally:
        lock_path.unlink(missing_ok=True)
    return base, expanded


def verify_slim_version(
    slim_path: str | Path,
    *,
    required_version: str = REQUIRED_SLIM_VERSION,
) -> dict[str, str]:
    """Require the stdpopsim-compatible SLiM 4.2.2 executable."""
    executable = Path(slim_path).resolve()
    if not executable.is_file():
        raise FileNotFoundError(f"SLiM executable not found: {executable}")
    completed = subprocess.run(
        [str(executable), "-v"],
        text=True,
        capture_output=True,
        check=False,
        timeout=30,
    )
    output = (completed.stdout + completed.stderr).strip()
    if completed.returncode != 0:
        raise RuntimeError(
            f"could not verify SLiM version (exit {completed.returncode}): {output}"
        )
    match = re.search(r"\bSLiM(?:\s+version)?\s+(\d+\.\d+\.\d+)\b", output)
    observed = match.group(1) if match else None
    if observed != required_version:
        raise RuntimeError(
            "this campaign requires SLiM "
            f"{required_version} for stdpopsim 0.3 compatibility; observed "
            f"{observed or 'an unparseable version'} from {executable}. "
            "SLiM 5.2 is intentionally rejected because the current generated "
            "model is not compatible with it."
        )
    return {
        "path": str(executable),
        "sha256": study._sha256(executable),  # noqa: SLF001
        "version": observed,
        "version_output": output,
    }


def _phase_invocation_context(
    args: argparse.Namespace,
    runtime: study.StudyRuntime,
    campaign_dir: Path,
    *,
    slim_path: Path,
    decoder_path: Path,
) -> dict[str, Any]:
    design_path = campaign_dir / "study_design.json"
    table_path = campaign_dir / "specifications.tsv"
    wrapper_path = (
        Path(runtime.repo_root) / "scripts" / "run_eas_introgression_replicates.py"
    ).resolve()
    design = json.loads(design_path.read_text(encoding="utf-8"))
    context = {
        "phase": args.phase,
        "phase_actions": (
            ["simulate", "decode", "plot"] if args.phase == "all" else [args.phase]
        ),
        "campaign": args.campaign,
        "replicates": int(args.replicates),
        "specification_filters": list(args.spec),
        "repo_root": str(Path(runtime.repo_root).resolve()),
        "campaign_directory": str(campaign_dir.resolve()),
        "runtime": study._portable_runtime_record(runtime),  # noqa: SLF001
        "plan": {
            "study_design_path": str(design_path.resolve()),
            "study_design_sha256": study._sha256(design_path),  # noqa: SLF001
            "specifications_path": str(table_path.resolve()),
            "specifications_sha256": study._sha256(table_path),  # noqa: SLF001
            "campaign_contract_sha256": design["campaign_contract_sha256"],
        },
        "requested_tools": {
            "slim_path": str(slim_path.resolve()),
            "decoder_path": str(decoder_path.resolve()),
        },
        "sources": {
            "orchestrator_path": str(Path(__file__).resolve()),
            "orchestrator_sha256": study._sha256(__file__),  # noqa: SLF001
            "core_orchestrator_path": str(Path(study.__file__).resolve()),
            "core_orchestrator_sha256": study._sha256(  # noqa: SLF001
                Path(study.__file__)
            ),
            "wrapper_path": str(wrapper_path),
            "wrapper_sha256": (
                study._sha256(wrapper_path) if wrapper_path.is_file() else None  # noqa: SLF001
            ),
        },
        "python": {
            "executable": str(Path(sys.executable).resolve()),
            "version": sys.version.split()[0],
        },
    }
    if args.phase == "recover":
        context["recovery_request"] = {
            "interruption_stop_utc": args.interruption_stop_utc,
            "interruption_note": args.interruption_note,
        }
    return context


def _append_phase_invocation_event(
    campaign_dir: Path,
    *,
    invocation_id: str,
    event: str,
    context: Mapping[str, Any],
    elapsed_seconds: float | None = None,
    tool_provenance: Mapping[str, Any] | None = None,
    error: BaseException | None = None,
) -> Path:
    """Write one immutable event file; existing provenance is never rewritten."""
    if event not in {"started", "completed", "failed"}:
        raise ValueError(f"unsupported phase invocation event: {event!r}")
    reserved = _PHASE_EVENT_RESERVED_FIELDS.intersection(context)
    if reserved:
        raise ValueError(
            "phase invocation context contains reserved event fields: "
            + ", ".join(sorted(reserved))
        )
    recorded_at = datetime.now(timezone.utc)
    payload: dict[str, Any] = {
        "schema": PHASE_INVOCATION_SCHEMA,
        "invocation_id": invocation_id,
        "event": event,
        "recorded_at_utc": recorded_at.isoformat(),
        "process_id": os.getpid(),
        **dict(context),
    }
    if elapsed_seconds is not None:
        payload["elapsed_seconds"] = float(elapsed_seconds)
    if tool_provenance:
        payload["tool_provenance"] = dict(tool_provenance)
    if error is not None:
        payload["error"] = {
            "type": type(error).__name__,
            "message": str(error),
        }
    invocation_dir = campaign_dir / PHASE_INVOCATION_DIRECTORY
    invocation_dir.mkdir(parents=True, exist_ok=True)
    stamp = recorded_at.strftime("%Y%m%dT%H%M%S.%fZ")
    destination = invocation_dir / f"{stamp}_{invocation_id}_{event}.json"
    if destination.exists():
        raise RuntimeError(f"phase invocation provenance collision: {destination}")
    return study._atomic_json(destination, payload)  # noqa: SLF001


def _parse_aware_utc(value: str, *, label: str) -> datetime:
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as error:
        raise ValueError(f"{label} is not a valid ISO-8601 timestamp") from error
    if parsed.tzinfo is None:
        raise ValueError(f"{label} must include a UTC offset")
    return parsed.astimezone(timezone.utc)


def _phase_event_context(record: Mapping[str, Any]) -> dict[str, Any]:
    """Recover only the original context, never event-owned reserved fields."""
    return {
        key: value
        for key, value in record.items()
        if key not in _PHASE_EVENT_RESERVED_FIELDS
    }


def _close_unmatched_phase_invocations(
    campaign_dir: Path,
    *,
    current_invocation_id: str,
    interruption_stop: datetime,
    interruption_note: str | None,
) -> dict[str, Any]:
    invocation_dir = campaign_dir / PHASE_INVOCATION_DIRECTORY
    if not invocation_dir.is_dir():
        return {
            "unmatched_started_count": 0,
            "closed_invocation_ids": [],
            "started_after_stop_invocation_ids": [],
        }
    grouped: dict[str, list[dict[str, Any]]] = {}
    for path in sorted(invocation_dir.glob("*.json")):
        try:
            record = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            raise ValueError(f"phase invocation event is unreadable: {path}") from error
        if (
            not isinstance(record, dict)
            or record.get("schema") != PHASE_INVOCATION_SCHEMA
            or record.get("event") not in {"started", "completed", "failed"}
            or not isinstance(record.get("invocation_id"), str)
        ):
            raise ValueError(f"phase invocation event has an invalid schema: {path}")
        grouped.setdefault(record["invocation_id"], []).append(record)

    closed: list[str] = []
    after_stop: list[str] = []
    unmatched_count = 0
    for invocation_id, records in sorted(grouped.items()):
        if invocation_id == current_invocation_id:
            continue
        starts = [record for record in records if record["event"] == "started"]
        terminals = [
            record for record in records if record["event"] in {"completed", "failed"}
        ]
        if terminals:
            continue
        if len(starts) != 1:
            raise ValueError(
                "unmatched phase invocation must have exactly one started event: "
                f"{invocation_id}"
            )
        unmatched_count += 1
        started_at = _parse_aware_utc(
            str(starts[0].get("recorded_at_utc")),
            label=f"started event {invocation_id}",
        )
        if started_at > interruption_stop:
            after_stop.append(invocation_id)
            continue
        elapsed_seconds = (interruption_stop - started_at).total_seconds()
        note_suffix = f" Note: {interruption_note}" if interruption_note else ""
        interruption = ExternalProcessInterruption(
            "phase process ended without a terminal event at "
            f"{interruption_stop.isoformat()}.{note_suffix}"
        )
        _append_phase_invocation_event(
            campaign_dir,
            invocation_id=invocation_id,
            event="failed",
            context=_phase_event_context(starts[0]),
            elapsed_seconds=elapsed_seconds,
            tool_provenance={
                "interruption_recovery": {
                    "recovered_by_invocation_id": current_invocation_id,
                    "interruption_stop_utc": interruption_stop.isoformat(),
                    "interruption_note": interruption_note,
                }
            },
            error=interruption,
        )
        closed.append(invocation_id)
    return {
        "unmatched_started_count": unmatched_count,
        "closed_invocation_ids": closed,
        "started_after_stop_invocation_ids": after_stop,
    }


def _recorded_slim_path(specification_dir: Path, fallback: Path) -> Path:
    contract_path = specification_dir / "simulation_contract.json"
    try:
        contract = json.loads(contract_path.read_text(encoding="utf-8"))
        recorded = contract["software"]["slim_path"]
        if isinstance(recorded, str) and recorded:
            return Path(recorded)
    except (KeyError, OSError, TypeError, json.JSONDecodeError):
        pass
    return fallback


def _recover_interrupted_campaign(
    runtime: study.StudyRuntime,
    expanded: Sequence[Mapping[str, Any]],
    campaign_dir: Path,
    *,
    slim_path: Path,
    current_invocation_id: str,
    interruption_stop: datetime,
    interruption_note: str | None,
) -> dict[str, Any]:
    """Recover stale simulation locks and append missing phase terminals."""
    recovered_locks: list[dict[str, Any]] = []
    for specification in expanded:
        specification_id = str(specification["specification_id"])
        specification_dir = campaign_dir / "work" / specification_id
        lock_path = specification_dir / ".simulation.lock"
        if not lock_path.is_file():
            continue
        task = {
            "specification": dict(specification),
            "runtime": asdict(runtime),
            "slim_path": str(
                _recorded_slim_path(specification_dir, slim_path).resolve()
            ),
            "specification_dir": str(specification_dir),
        }
        result = study.record_external_simulation_interruption(
            task,
            stopped_at_utc=interruption_stop.isoformat(),
        )
        # A valid completion is authoritative in the core helper and can make
        # it return before consuming the stale lock.  Entering the same core
        # lock once archives a dead local owner and still refuses live/remote
        # owners, leaving this recovery idempotent.
        if lock_path.exists():
            with study._specification_execution_lock(specification_dir):  # noqa: SLF001
                pass
        payload = result.get("failure", result.get("completion", {}))
        recovered_locks.append(
            {
                "specification_id": specification_id,
                "result_status": result.get("status"),
                "elapsed_seconds": payload.get("elapsed_seconds"),
                "error_type": payload.get("error_type"),
            }
        )

    phase_recovery = _close_unmatched_phase_invocations(
        campaign_dir,
        current_invocation_id=current_invocation_id,
        interruption_stop=interruption_stop,
        interruption_note=interruption_note,
    )
    return {
        "interruption_stop_utc": interruption_stop.isoformat(),
        "interruption_note": interruption_note,
        "recovered_simulation_lock_count": len(recovered_locks),
        "recovered_simulation_locks": recovered_locks,
        "phase_invocations": phase_recovery,
        "no_op": not recovered_locks and not phase_recovery["closed_invocation_ids"],
    }


def _decoder_provenance(decoder_path: Path) -> dict[str, Any]:
    resolved = decoder_path.resolve()
    record: dict[str, Any] = {
        "path": str(resolved),
        "exists": resolved.is_file(),
    }
    if resolved.is_file():
        record["sha256"] = study._sha256(resolved)  # noqa: SLF001
    return record


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the replicated EAS introgressed-sweep campaign"
    )
    parser.add_argument(
        "phase", choices=("plan", "recover", "simulate", "decode", "plot", "all")
    )
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--campaign", default=DEFAULT_CAMPAIGN)
    parser.add_argument("--replicates", type=int, default=DEFAULT_REPLICATES)
    parser.add_argument("--slim-bin", type=Path)
    parser.add_argument("--decoder-bin", type=Path)
    parser.add_argument(
        "--sample-diploids", type=int, default=study.DEFAULT_SAMPLE_DIPLOIDS
    )
    parser.add_argument("--workers", type=int, default=study.DEFAULT_SIMULATION_WORKERS)
    parser.add_argument("--threads", type=int, default=study.DEFAULT_THREADS)
    parser.add_argument(
        "--slim-scaling-factor", type=float, default=study.DEFAULT_SLIM_SCALING_FACTOR
    )
    parser.add_argument(
        "--slim-burn-in", type=float, default=study.DEFAULT_SLIM_BURN_IN
    )
    parser.add_argument(
        "--max-external-draws", type=int, default=study.DEFAULT_MAX_EXTERNAL_DRAWS
    )
    parser.add_argument(
        "--max-launched-draws",
        type=int,
        default=40,
        help=(
            "cap launched seed attempts, including interrupted trajectories, "
            "for each replicate slot"
        ),
    )
    parser.add_argument(
        "--spec-timeout-minutes",
        type=float,
        default=5.0,
        help="terminate and provenance-mark each replicate slot after this wall time",
    )
    parser.add_argument(
        "--cumulative-spec-timeout-minutes",
        type=float,
        default=30.0,
        help=(
            "cap cumulative wall time across resumptions of each replicate slot; "
            "0 is unlimited"
        ),
    )
    parser.add_argument(
        "--spec",
        action="append",
        default=[],
        help="run only expanded specification IDs containing this text (repeatable)",
    )
    parser.add_argument(
        "--interruption-stop-utc",
        help="required for recover; aware ISO-8601 time when the external run stopped",
    )
    parser.add_argument(
        "--interruption-note",
        help="optional operator note recorded with recovered terminal events",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.spec and args.phase in {"recover", "plot", "all"}:
        parser.error(
            "--spec cannot be used with recover, plot, or all because those are "
            "campaign-wide phases"
        )
    if args.phase == "recover" and not args.interruption_stop_utc:
        parser.error("recover requires --interruption-stop-utc")
    if args.phase != "recover" and (
        args.interruption_stop_utc is not None or args.interruption_note is not None
    ):
        parser.error(
            "--interruption-stop-utc and --interruption-note are only valid with recover"
        )
    interruption_stop = (
        _parse_aware_utc(
            args.interruption_stop_utc,
            label="--interruption-stop-utc",
        )
        if args.phase == "recover"
        else None
    )
    repo_root = args.repo_root.resolve()
    runtime = study.StudyRuntime(
        repo_root=str(repo_root),
        sample_diploids=args.sample_diploids,
        threads=args.threads,
        simulation_workers=args.workers,
        slim_scaling_factor=args.slim_scaling_factor,
        slim_burn_in=args.slim_burn_in,
        max_external_draws=args.max_external_draws,
        max_launched_draws=args.max_launched_draws,
        max_specification_seconds=args.spec_timeout_minutes * 60.0,
        max_cumulative_specification_seconds=(
            args.cumulative_spec_timeout_minutes * 60.0
        ),
        campaign_subdirectory=args.campaign,
    )
    base, expanded = write_introgression_replicate_plan(runtime, args.replicates)
    specifications = {"introgression": expanded}
    slim_path = Path(
        args.slim_bin or study._default_slim(repo_root)  # noqa: SLF001
    )
    decoder_path = Path(args.decoder_bin or (repo_root / "bin" / "gamma_smc"))
    campaign_dir = study._runtime_study_dir(  # noqa: SLF001
        runtime, "introgression"
    )
    invocation_id = uuid.uuid4().hex
    invocation_context = _phase_invocation_context(
        args,
        runtime,
        campaign_dir,
        slim_path=slim_path,
        decoder_path=decoder_path,
    )
    _append_phase_invocation_event(
        campaign_dir,
        invocation_id=invocation_id,
        event="started",
        context=invocation_context,
    )
    started = time.perf_counter()
    tool_provenance: dict[str, Any] = {}
    try:
        if args.phase == "recover":
            if interruption_stop is None:  # pragma: no cover - parser invariant
                raise RuntimeError("recover phase has no interruption stop time")
            tool_provenance["recovery"] = _recover_interrupted_campaign(
                runtime,
                expanded,
                campaign_dir,
                slim_path=slim_path,
                current_invocation_id=invocation_id,
                interruption_stop=interruption_stop,
                interruption_note=args.interruption_note,
            )
        if args.phase in {"simulate", "all"}:
            tool_provenance["slim"] = verify_slim_version(slim_path)
            study.simulate_studies(
                runtime,
                specifications,
                slim_path=slim_path,
                specification_filters=args.spec,
            )
        if args.phase in {"decode", "all"}:
            tool_provenance["decoder"] = _decoder_provenance(decoder_path)
            study.decode_studies(
                runtime,
                specifications,
                decoder_path=decoder_path,
                specification_filters=args.spec,
            )
        if args.phase in {"plot", "all"}:
            from .eas_replicate_analysis import aggregate_introgression_replicates

            aggregate_introgression_replicates(
                runtime,
                base,
                expanded,
                campaign_dir,
                decoder_path=(
                    decoder_path if args.decoder_bin or args.phase == "all" else None
                ),
            )
    except BaseException as error:
        _append_phase_invocation_event(
            campaign_dir,
            invocation_id=invocation_id,
            event="failed",
            context=invocation_context,
            elapsed_seconds=time.perf_counter() - started,
            tool_provenance=tool_provenance,
            error=error,
        )
        raise
    _append_phase_invocation_event(
        campaign_dir,
        invocation_id=invocation_id,
        event="completed",
        context=invocation_context,
        elapsed_seconds=time.perf_counter() - started,
        tool_provenance=tool_provenance,
    )
    return 0


__all__ = [
    "BIOLOGICAL_SPECIFICATION_COUNT",
    "DEFAULT_CAMPAIGN",
    "DEFAULT_REPLICATES",
    "ExternalProcessInterruption",
    "PHASE_INVOCATION_DIRECTORY",
    "PHASE_INVOCATION_SCHEMA",
    "REPLICATE_SCHEMA",
    "REQUIRED_SLIM_VERSION",
    "build_introgression_replicates",
    "build_parser",
    "expand_introgression_specifications",
    "main",
    "verify_slim_version",
    "write_introgression_replicate_plan",
]
