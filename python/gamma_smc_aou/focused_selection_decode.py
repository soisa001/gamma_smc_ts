"""One-pass Gamma-SMC posterior reduction for the focused EAS campaign.

Gamma-SMC writes chunk-major Gamma shape/rate parameters for every decoded
pair.  This module reads that stream once and derives all focal-genotype
classes from the same posterior.  The compact outputs are checksum contracted
so that a completed unit can be validated after its large raw stream is
removed.
"""

from __future__ import annotations

import hashlib
import json
import os
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scipy
import zstandard
from scipy.special import gammainc

from .eas_sweep_models import (
    GENERATION_TIME_YEARS,
    MUTATION_RATE,
    TMRCA_THRESHOLDS_YEARS,
)

SCHEMA_VERSION = "gamma-smc.focused-selection-decode/v1"
GENOTYPE_CLASSES = ("overall", "hom_ref", "heterozygous", "hom_alt")
PAIR_CLASSES = GENOTYPE_CLASSES[1:]
STANDARD_UNIT_COLUMNS = (
    "unit_id",
    "demography_id",
    "demography_kind",
    "population",
    "simulation_class",
    "selection_coefficient",
    "target_allele_frequency",
    "replicate_index",
    "seed",
)
OUTPUT_FILENAMES = {
    "spatial_profiles": "spatial_profiles.tsv.gz",
    "pair_summaries": "pair_summaries.tsv.gz",
    "class_summaries": "class_summaries.tsv",
}
POSTPROCESSOR_SOURCE = Path(__file__).resolve()


@dataclass(frozen=True)
class DecodeArtifacts:
    """Validated compact outputs for one decoded simulation unit."""

    output_dir: Path
    spatial_profiles: Path
    pair_summaries: Path
    class_summaries: Path
    completion: Path
    cache_hit: bool


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_sha256(value: Any) -> str:
    payload = json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _atomic_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        temporary.write_text(
            json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_frame(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    compression: str | Mapping[str, Any] | None = None
    if path.suffix == ".gz":
        compression = {"method": "gzip", "compresslevel": 6, "mtime": 0}
    try:
        frame.to_csv(
            temporary,
            sep="\t",
            index=False,
            compression=compression,
            float_format="%.12g",
        )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _read_exact(reader, size: int) -> bytes:
    chunks: list[bytes] = []
    remaining = size
    while remaining:
        block = reader.read(remaining)
        if not block:
            raise EOFError(f"posterior stream ended {remaining} bytes early")
        chunks.append(block)
        remaining -= len(block)
    return b"".join(chunks)


def _coerce_positive_float(value: Any, label: str) -> float:
    result = float(value)
    if not np.isfinite(result) or result <= 0:
        raise ValueError(f"{label} must be finite and positive")
    return result


def _coerce_positive_int(value: Any, label: str) -> int:
    numeric = float(value)
    if not np.isfinite(numeric) or numeric < 1 or numeric != np.floor(numeric):
        raise ValueError(f"{label} must be a positive integer")
    return int(numeric)


def _coerce_thresholds(values: Sequence[float]) -> tuple[float, ...]:
    thresholds = tuple(float(value) for value in values)
    if not thresholds:
        raise ValueError("at least one TMRCA threshold is required")
    array = np.asarray(thresholds, dtype=float)
    if (
        not np.all(np.isfinite(array))
        or np.any(array <= 0)
        or np.any(np.diff(array) <= 0)
    ):
        raise ValueError("TMRCA thresholds must be finite, positive, and increasing")
    return thresholds


def _read_pair_table(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise ValueError(f"pair table is absent: {path}")
    frame = pd.read_csv(path, sep="\t")
    required = {
        "gamma_smc_haplotype_0",
        "gamma_smc_haplotype_1",
        "genotype_class",
    }
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"pair table is missing columns: {', '.join(missing)}")
    if frame.empty:
        raise ValueError("pair table contains no diploids")
    pairs = frame[["gamma_smc_haplotype_0", "gamma_smc_haplotype_1"]].to_numpy(
        dtype=np.int64
    )
    numeric = frame[["gamma_smc_haplotype_0", "gamma_smc_haplotype_1"]].to_numpy(
        dtype=float
    )
    if not np.all(np.isfinite(numeric)) or not np.array_equal(numeric, pairs):
        raise ValueError("pair-table haplotypes must be finite integers")
    if np.any(pairs < 0) or np.any(pairs[:, 0] >= pairs[:, 1]):
        raise ValueError("pair-table haplotypes must be nonnegative ordered pairs")
    if len(np.unique(pairs, axis=0)) != len(pairs):
        raise ValueError("pair table contains duplicate haplotype pairs")
    classes = frame["genotype_class"].astype(str)
    unexpected = sorted(set(classes).difference(PAIR_CLASSES))
    if unexpected:
        raise ValueError(f"unexpected genotype classes: {', '.join(unexpected)}")
    absent = sorted(set(PAIR_CLASSES).difference(classes))
    if absent:
        raise ValueError(
            "focused comparison requires every focal genotype class; absent: "
            + ", ".join(absent)
        )
    result = frame.copy()
    result["genotype_class"] = classes
    return result


def _read_raw_metadata(raw_path: Path, pair_table: pd.DataFrame) -> dict[str, Any]:
    metadata_path = raw_path.with_name(raw_path.name + ".meta")
    if not raw_path.is_file() or not metadata_path.is_file():
        raise ValueError("raw posterior and its .meta file are required")
    try:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError("raw posterior metadata is unreadable") from error
    for key in (
        "num_pairs",
        "chunk_size",
        "sequence_length",
        "output_positions",
        "pairs",
        "scaled_mutation_rate",
    ):
        if key not in metadata:
            raise ValueError(f"raw posterior metadata is missing {key}")
    n_pairs = _coerce_positive_int(metadata["num_pairs"], "metadata num_pairs")
    chunk_size = _coerce_positive_int(metadata["chunk_size"], "metadata chunk_size")
    n_positions = _coerce_positive_int(
        metadata["sequence_length"], "metadata sequence_length"
    )
    if n_pairs != len(pair_table) or chunk_size < 1 or n_positions < 1:
        raise ValueError("raw posterior metadata has inconsistent dimensions")
    positions = np.asarray(metadata["output_positions"], dtype=np.int64)
    numeric_positions = np.asarray(metadata["output_positions"], dtype=float)
    if (
        positions.ndim != 1
        or len(positions) != n_positions
        or not np.array_equal(positions, numeric_positions)
        or np.any(positions < 0)
        or np.any(np.diff(positions) <= 0)
    ):
        raise ValueError(
            "raw posterior output positions must be nonnegative increasing integers"
        )
    numeric_raw_pairs = np.asarray(metadata["pairs"], dtype=float)
    raw_pairs = np.asarray(metadata["pairs"], dtype=np.int64)
    expected_pairs = pair_table[
        ["gamma_smc_haplotype_0", "gamma_smc_haplotype_1"]
    ].to_numpy(dtype=np.int64)
    if (
        raw_pairs.shape != expected_pairs.shape
        or not np.all(np.isfinite(numeric_raw_pairs))
        or not np.array_equal(numeric_raw_pairs, raw_pairs)
        or not np.array_equal(raw_pairs, expected_pairs)
    ):
        raise ValueError(
            "raw posterior pair order does not exactly match the diploid pair table"
        )
    _coerce_positive_float(metadata["scaled_mutation_rate"], "scaled mutation rate")
    return metadata


def _unit_prefix(unit_record: Mapping[str, Any]) -> dict[str, Any]:
    prefix: dict[str, Any] = {}
    for column in STANDARD_UNIT_COLUMNS:
        if column in unit_record:
            value = unit_record[column]
            if isinstance(value, np.generic):
                value = value.item()
            prefix[column] = value
    if not str(prefix.get("unit_id", "")).strip():
        raise ValueError("unit_record must contain a nonempty unit_id")
    return prefix


def _static_contract(
    *,
    unit_record: Mapping[str, Any],
    pair_table_path: Path,
    pair_table_sha256: str,
    mutation_rate: float,
    generation_time_years: float,
    thresholds_years: tuple[float, ...],
    focal_position_0based: float,
    auxiliary_output_paths: Mapping[str, Path],
) -> dict[str, Any]:
    unit = {
        str(key): (value.item() if isinstance(value, np.generic) else value)
        for key, value in unit_record.items()
    }
    return {
        "schema": SCHEMA_VERSION,
        "unit_record": unit,
        "pair_table": {
            "path": str(pair_table_path.resolve()),
            "sha256": pair_table_sha256,
        },
        "parameters": {
            "unscaled_mutation_rate_per_bp_per_generation": mutation_rate,
            "generation_time_years": generation_time_years,
            "thresholds_years": list(thresholds_years),
            "focal_position_0based": focal_position_0based,
        },
        "implementation": {
            "postprocessor_path": str(POSTPROCESSOR_SOURCE),
            "postprocessor_sha256": sha256_file(POSTPROCESSOR_SOURCE),
        },
        "auxiliary_output_paths": {
            str(label): str(path.resolve())
            for label, path in sorted(auxiliary_output_paths.items())
        },
    }


def _output_record(path: Path, frame: pd.DataFrame | None = None) -> dict[str, Any]:
    record: dict[str, Any] = {
        "path": path.name,
        "sha256": sha256_file(path),
        "size_bytes": path.stat().st_size,
    }
    if frame is not None:
        record["rows"] = len(frame)
        record["columns"] = list(frame.columns)
    return record


def _validate_cached_outputs(
    *,
    output_dir: Path,
    completion: Mapping[str, Any],
    static_contract: Mapping[str, Any],
    raw_path: Path,
) -> DecodeArtifacts:
    if completion.get("schema") != SCHEMA_VERSION:
        raise ValueError("decode completion schema is incompatible")
    if completion.get("status") != "complete":
        raise ValueError("decode completion status is not complete")
    recorded_contract = completion.get("contract")
    if not isinstance(recorded_contract, Mapping) or completion.get(
        "contract_sha256"
    ) != _canonical_sha256(recorded_contract):
        raise ValueError("decode completion contract checksum is invalid")
    recorded_static = recorded_contract.get("static")
    if _canonical_sha256(recorded_static) != _canonical_sha256(static_contract):
        raise ValueError("existing decode contract is incompatible")
    raw_metadata_path = raw_path.with_name(raw_path.name + ".meta")
    if raw_path.is_file() != raw_metadata_path.is_file():
        raise ValueError("only one of the raw posterior and metadata remains")
    recorded_inputs = recorded_contract.get("raw_inputs", {})
    if not isinstance(recorded_inputs, Mapping):
        raise TypeError("decode completion input manifest is malformed")
    if raw_path.is_file():
        if sha256_file(raw_path) != recorded_inputs.get("posterior_sha256"):
            raise ValueError("raw posterior checksum differs from completed decode")
        if sha256_file(raw_metadata_path) != recorded_inputs.get("metadata_sha256"):
            raise ValueError("raw posterior metadata checksum differs from completion")
    auxiliary_paths = {
        str(label): Path(path)
        for label, path in static_contract.get("auxiliary_output_paths", {}).items()
    }
    expected_outputs = {*OUTPUT_FILENAMES, *auxiliary_paths}
    outputs = completion.get("outputs")
    if not isinstance(outputs, Mapping) or set(outputs) != expected_outputs:
        raise ValueError("decode completion output manifest is malformed")
    paths: dict[str, Path] = {}
    for label, filename in OUTPUT_FILENAMES.items():
        record = outputs[label]
        if not isinstance(record, Mapping):
            raise TypeError(f"decode completion output record is malformed: {label}")
        path = output_dir / filename
        paths[label] = path
        if record.get("path") != filename or not path.is_file():
            raise ValueError(f"completed decode output is absent: {filename}")
        if sha256_file(path) != record.get("sha256"):
            raise ValueError(f"completed decode output checksum failed: {filename}")
        if path.stat().st_size != int(record.get("size_bytes", -1)):
            raise ValueError(f"completed decode output size failed: {filename}")
        try:
            frame = pd.read_csv(path, sep="\t")
        except Exception as error:
            raise ValueError(
                f"completed decode output is unreadable: {filename}"
            ) from error
        if len(frame) != int(record.get("rows", -1)) or list(frame.columns) != list(
            record.get("columns", [])
        ):
            raise ValueError(f"completed decode output schema failed: {filename}")
    for label, path in auxiliary_paths.items():
        record = outputs[label]
        if not isinstance(record, Mapping):
            raise TypeError(f"decode auxiliary output record is malformed: {label}")
        if record.get("path") != path.name or not path.is_file():
            raise ValueError(f"completed decode auxiliary output is absent: {label}")
        if sha256_file(path) != record.get("sha256") or path.stat().st_size != int(
            record.get("size_bytes", -1)
        ):
            raise ValueError(f"decode auxiliary output checksum failed: {label}")
    return DecodeArtifacts(
        output_dir=output_dir,
        spatial_profiles=paths["spatial_profiles"],
        pair_summaries=paths["pair_summaries"],
        class_summaries=paths["class_summaries"],
        completion=output_dir / "completion.json",
        cache_hit=True,
    )


def _posterior_arrays(
    raw_path: Path,
    metadata: Mapping[str, Any],
    *,
    mutation_rate: float,
    generation_time_years: float,
    thresholds_years: tuple[float, ...],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    n_pairs = int(metadata["num_pairs"])
    chunk_size = int(metadata["chunk_size"])
    n_positions = int(metadata["sequence_length"])
    scaled_rate = float(metadata["scaled_mutation_rate"])
    two_ne_generations = scaled_rate / (2.0 * mutation_rate)
    threshold_scaled = (
        np.asarray(thresholds_years, dtype=np.float64)
        / generation_time_years
        / two_ne_generations
    )
    probabilities = np.full(
        (len(thresholds_years), n_positions, n_pairs), np.nan, dtype=np.float64
    )
    tmrca = np.full((n_positions, n_pairs), np.nan, dtype=np.float64)
    valid = np.zeros((n_positions, n_pairs), dtype=bool)
    block_bytes = n_positions * chunk_size * np.dtype(np.float32).itemsize
    with (
        raw_path.open("rb") as compressed,
        zstandard.ZstdDecompressor().stream_reader(compressed) as reader,
    ):
        for pair_start in range(0, n_pairs, chunk_size):
            real_pairs = min(chunk_size, n_pairs - pair_start)
            alpha = np.frombuffer(
                _read_exact(reader, block_bytes), dtype=np.float32
            ).reshape(n_positions, chunk_size)[:, :real_pairs]
            beta = np.frombuffer(
                _read_exact(reader, block_bytes), dtype=np.float32
            ).reshape(n_positions, chunk_size)[:, :real_pairs]
            chunk_valid = (
                np.isfinite(alpha) & np.isfinite(beta) & (alpha > 0) & (beta > 0)
            )
            pair_end = pair_start + real_pairs
            valid[:, pair_start:pair_end] = chunk_valid
            chunk_tmrca = np.full(alpha.shape, np.nan, dtype=np.float64)
            chunk_tmrca[chunk_valid] = (
                alpha[chunk_valid].astype(np.float64)
                / beta[chunk_valid].astype(np.float64)
                * two_ne_generations
            )
            tmrca[:, pair_start:pair_end] = chunk_tmrca
            for threshold_index, scaled_threshold in enumerate(threshold_scaled):
                chunk_probability = np.full(alpha.shape, np.nan, dtype=np.float64)
                chunk_probability[chunk_valid] = gammainc(
                    alpha[chunk_valid].astype(np.float64),
                    beta[chunk_valid].astype(np.float64) * scaled_threshold,
                )
                probabilities[threshold_index, :, pair_start:pair_end] = (
                    chunk_probability
                )
        if reader.read(1):
            raise ValueError("raw posterior stream has trailing decompressed bytes")
    if np.any(valid.sum(axis=0) == 0):
        bad = np.flatnonzero(valid.sum(axis=0) == 0).tolist()
        raise ValueError(f"raw posterior contains all-invalid pairs: {bad[:10]}")
    if np.any(valid.sum(axis=1) == 0):
        bad = np.flatnonzero(valid.sum(axis=1) == 0).tolist()
        raise ValueError(f"raw posterior contains all-invalid positions: {bad[:10]}")
    return probabilities, tmrca, valid, float(two_ne_generations)


def _build_frames(
    *,
    unit_record: Mapping[str, Any],
    pair_table: pd.DataFrame,
    metadata: Mapping[str, Any],
    probabilities: np.ndarray,
    tmrca: np.ndarray,
    valid: np.ndarray,
    thresholds_years: tuple[float, ...],
    focal_position_0based: float,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    prefix = _unit_prefix(unit_record)
    positions = np.asarray(metadata["output_positions"], dtype=np.int64)
    focal_matches = np.flatnonzero(
        positions.astype(np.float64) == float(focal_position_0based)
    )
    if len(focal_matches) != 1:
        raise ValueError(
            "raw posterior output grid must contain the exact focal position once"
        )
    focal_index = int(focal_matches[0])
    focal_output_position = int(positions[focal_index])
    classes = pair_table["genotype_class"].astype(str).to_numpy()
    haplotypes = pair_table[
        ["gamma_smc_haplotype_0", "gamma_smc_haplotype_1"]
    ].to_numpy(dtype=np.int64)

    pair_rows: list[dict[str, Any]] = []
    for pair_index, genotype_class in enumerate(classes):
        pair_valid = valid[:, pair_index]
        for threshold_index, threshold in enumerate(thresholds_years):
            values = probabilities[threshold_index, :, pair_index]
            pair_rows.append(
                {
                    **prefix,
                    "source": "gamma_smc",
                    "pair_index": pair_index,
                    "gamma_smc_haplotype_0": int(haplotypes[pair_index, 0]),
                    "gamma_smc_haplotype_1": int(haplotypes[pair_index, 1]),
                    "genotype_class": genotype_class,
                    "threshold_years": threshold,
                    "region_mean_p_tmrca_lt_threshold": float(np.nanmean(values)),
                    "focal_p_tmrca_lt_threshold": float(values[focal_index]),
                    "region_mean_tmrca_generations": float(
                        np.nanmean(tmrca[:, pair_index])
                    ),
                    "focal_tmrca_generations": float(tmrca[focal_index, pair_index]),
                    "n_valid_positions": int(pair_valid.sum()),
                    "focal_output_position_0based": focal_output_position,
                    "focal_offset_bp": float(
                        focal_output_position - focal_position_0based
                    ),
                }
            )
    pair_summaries = pd.DataFrame(pair_rows)

    class_rows: list[dict[str, Any]] = []
    spatial_rows: list[dict[str, Any]] = []
    for genotype_class in GENOTYPE_CLASSES:
        indices = (
            np.arange(len(classes), dtype=np.int64)
            if genotype_class == "overall"
            else np.flatnonzero(classes == genotype_class)
        )
        if not len(indices):
            raise ValueError(f"genotype class has no decoded pairs: {genotype_class}")
        class_valid = valid[:, indices]
        if not class_valid[focal_index, :].any():
            raise ValueError(
                f"genotype class has no valid posterior at the focal output: "
                f"{genotype_class}"
            )
        for threshold_index, threshold in enumerate(thresholds_years):
            # Index the threshold first so NumPy does not move the advanced
            # pair index ahead of the position dimension.
            class_probability = probabilities[threshold_index][:, indices]
            class_rows.append(
                {
                    **prefix,
                    "source": "gamma_smc",
                    "genotype_class": genotype_class,
                    "threshold_years": threshold,
                    "region_mean_p_tmrca_lt_threshold": float(
                        np.nanmean(class_probability)
                    ),
                    "focal_mean_p_tmrca_lt_threshold": float(
                        np.nanmean(class_probability[focal_index, :])
                    ),
                    "region_mean_tmrca_generations": float(
                        np.nanmean(tmrca[:, indices])
                    ),
                    "focal_mean_tmrca_generations": float(
                        np.nanmean(tmrca[focal_index, indices])
                    ),
                    "n_pairs": len(indices),
                    "n_valid_posterior_cells_region": int(class_valid.sum()),
                    "n_valid_pairs_focal": int(class_valid[focal_index, :].sum()),
                    "focal_output_position_0based": focal_output_position,
                    "focal_offset_bp": float(
                        focal_output_position - focal_position_0based
                    ),
                }
            )
            for position_index, position in enumerate(positions):
                position_probability = class_probability[position_index, :]
                position_tmrca = tmrca[position_index, indices]
                spatial_rows.append(
                    {
                        **prefix,
                        "source": "gamma_smc",
                        "genotype_class": genotype_class,
                        "threshold_years": threshold,
                        "position_0based": int(position),
                        "mean_p_tmrca_lt_threshold": float(
                            np.nanmean(position_probability)
                        ),
                        "mean_tmrca_generations": float(np.nanmean(position_tmrca)),
                        "n_valid_pairs": int(np.isfinite(position_probability).sum()),
                    }
                )
    class_summaries = pd.DataFrame(class_rows)
    spatial_profiles = pd.DataFrame(spatial_rows)
    center = {
        "requested_focal_position_0based": focal_position_0based,
        "output_position_0based": focal_output_position,
        "offset_bp": float(focal_output_position - focal_position_0based),
        "output_position_index": focal_index,
    }
    return spatial_profiles, pair_summaries, class_summaries, center


def postprocess_raw_posteriors(
    raw_path: str | Path,
    pair_table_path: str | Path,
    output_dir: str | Path,
    *,
    unit_record: Mapping[str, Any],
    mutation_rate: float = MUTATION_RATE,
    generation_time_years: float = GENERATION_TIME_YEARS,
    thresholds_years: Sequence[float] = TMRCA_THRESHOLDS_YEARS,
    focal_position_0based: float,
    remove_raw_after_success: bool = False,
    auxiliary_output_paths: Mapping[str, str | Path] | None = None,
) -> DecodeArtifacts:
    """Reduce one raw Gamma stream to all genotype classes in one pass.

    A pre-existing completion is treated as an immutable cache entry.  Every
    compact file is re-read and checksum validated.  The raw stream may be
    absent only after such a completion exists; this supports intentional
    compaction without weakening restart validation.
    """

    raw_path = Path(raw_path).resolve()
    pair_table_path = Path(pair_table_path).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    completion_path = output_dir / "completion.json"
    mutation_rate = _coerce_positive_float(mutation_rate, "mutation rate")
    generation_time_years = _coerce_positive_float(
        generation_time_years, "generation time"
    )
    focal_position_0based = float(focal_position_0based)
    if not np.isfinite(focal_position_0based) or focal_position_0based < 0:
        raise ValueError("focal position must be finite and nonnegative")
    thresholds = _coerce_thresholds(thresholds_years)
    _unit_prefix(unit_record)
    pair_table = _read_pair_table(pair_table_path)
    auxiliary_paths = {
        str(label): Path(path).resolve()
        for label, path in (auxiliary_output_paths or {}).items()
    }
    if set(auxiliary_paths).intersection(OUTPUT_FILENAMES):
        raise ValueError("auxiliary decode output labels overlap compact outputs")
    for label, path in auxiliary_paths.items():
        if not label.strip() or not path.is_file():
            raise ValueError(f"decode auxiliary output is absent: {label}")
    static_contract = _static_contract(
        unit_record=unit_record,
        pair_table_path=pair_table_path,
        pair_table_sha256=sha256_file(pair_table_path),
        mutation_rate=mutation_rate,
        generation_time_years=generation_time_years,
        thresholds_years=thresholds,
        focal_position_0based=focal_position_0based,
        auxiliary_output_paths=auxiliary_paths,
    )

    if completion_path.is_file():
        try:
            completion = json.loads(completion_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            raise ValueError("decode completion is unreadable") from error
        artifacts = _validate_cached_outputs(
            output_dir=output_dir,
            completion=completion,
            static_contract=static_contract,
            raw_path=raw_path,
        )
        if remove_raw_after_success:
            raw_path.unlink(missing_ok=True)
            raw_path.with_name(raw_path.name + ".meta").unlink(missing_ok=True)
        return artifacts

    metadata = _read_raw_metadata(raw_path, pair_table)
    raw_metadata_path = raw_path.with_name(raw_path.name + ".meta")
    sequence_length_bp = unit_record.get("sequence_length_bp")
    if sequence_length_bp is not None:
        length = _coerce_positive_float(sequence_length_bp, "sequence length")
        positions = np.asarray(metadata["output_positions"], dtype=float)
        if np.any(positions >= length) or focal_position_0based >= length:
            raise ValueError("raw output or focal position lies outside the sequence")
    probabilities, tmrca, valid, two_ne_generations = _posterior_arrays(
        raw_path,
        metadata,
        mutation_rate=mutation_rate,
        generation_time_years=generation_time_years,
        thresholds_years=thresholds,
    )
    spatial, pairs, classes, center = _build_frames(
        unit_record=unit_record,
        pair_table=pair_table,
        metadata=metadata,
        probabilities=probabilities,
        tmrca=tmrca,
        valid=valid,
        thresholds_years=thresholds,
        focal_position_0based=focal_position_0based,
    )
    frames = {
        "spatial_profiles": spatial,
        "pair_summaries": pairs,
        "class_summaries": classes,
    }
    paths = {
        label: output_dir / filename for label, filename in OUTPUT_FILENAMES.items()
    }
    for label, frame in frames.items():
        _atomic_frame(paths[label], frame)
    contract = {
        "static": static_contract,
        "raw_inputs": {
            "posterior_path": str(raw_path),
            "posterior_sha256": sha256_file(raw_path),
            "metadata_path": str(raw_metadata_path),
            "metadata_sha256": sha256_file(raw_metadata_path),
        },
    }
    completion = {
        "schema": SCHEMA_VERSION,
        "status": "complete",
        "contract": contract,
        "contract_sha256": _canonical_sha256(contract),
        "posterior_scale": {
            "scaled_mutation_rate": float(metadata["scaled_mutation_rate"]),
            "unscaled_mutation_rate": mutation_rate,
            "two_ne_generations": two_ne_generations,
        },
        "posterior_dimensions": {
            "num_pairs": int(metadata["num_pairs"]),
            "num_output_positions": int(metadata["sequence_length"]),
            "chunk_size": int(metadata["chunk_size"]),
            "valid_cells": int(valid.sum()),
            "total_cells": int(valid.size),
        },
        "focal_center": center,
        "software": {
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "zstandard": zstandard.__version__,
        },
        "outputs": {
            **{label: _output_record(paths[label], frames[label]) for label in frames},
            **{
                label: _output_record(path)
                for label, path in sorted(auxiliary_paths.items())
            },
        },
    }
    _atomic_json(completion_path, completion)
    artifacts = DecodeArtifacts(
        output_dir=output_dir,
        spatial_profiles=paths["spatial_profiles"],
        pair_summaries=paths["pair_summaries"],
        class_summaries=paths["class_summaries"],
        completion=completion_path,
        cache_hit=False,
    )
    if remove_raw_after_success:
        raw_path.unlink()
        raw_metadata_path.unlink()
    return artifacts


__all__ = [
    "GENOTYPE_CLASSES",
    "OUTPUT_FILENAMES",
    "PAIR_CLASSES",
    "SCHEMA_VERSION",
    "DecodeArtifacts",
    "postprocess_raw_posteriors",
    "sha256_file",
]
