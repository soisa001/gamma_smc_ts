"""Post-hoc focal pair-class analysis for the conditioned CoSi2 study.

This module is deliberately additive.  It consumes the immutable Gamma-SMC
posterior streams and conditioned-study result receipt without changing the
decode or result modules whose hashes are already frozen.  The fixed
100-haplotype all-pairs panel is split by the allele at the focal mutation;
``hom_ref`` and ``hom_alt`` therefore describe unordered *haplotype-pair*
classes, not the original biological diploid sample pairing.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import os
import platform
import re
import shutil
import struct
import tempfile
from collections.abc import Mapping, Sequence
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

# Workbench can export the inline notebook backend even when matplotlib-inline
# is absent.  Set the headless backend before importing any Matplotlib module.
os.environ["MPLBACKEND"] = "Agg"

import matplotlib
import numpy as np
import pandas as pd
import scipy
import zstandard
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.special import gammainc

from . import cosi2_conditioned_gamma_workflow as conditioned_workflow
from . import cosi2_gamma_analysis as gamma_analysis
from . import cosi2_gamma_decode as gamma_decode

matplotlib.use("Agg", force=True)

SCHEMA_VERSION = "gamma-smc.cosi2-conditioned-focal-pairs/v2"
RESULTS_SCHEMA = "gamma-smc.cosi2-conditioned-focal-pair-results/v2"
MODULE_PATH = Path(__file__).resolve()

DEFAULT_BASE_RESULTS_RELATIVE = conditioned_workflow.DEFAULT_RESULTS_RELATIVE
DEFAULT_WORK_RELATIVE = conditioned_workflow.DEFAULT_WORK_RELATIVE
DEFAULT_RESULTS_RELATIVE = Path(
    "focused_selection_EAS_sim/results/cosi2_gamma_conditioned_s0p01_focal_pairs"
)
DEFAULT_MAX_WORKERS = 4
MAX_WORKERS = 24
FOCAL_POSITION_0BASED = gamma_decode.DEFAULT_FOCAL_POSITION_0BASED
MUTATION_RATE = gamma_decode.DEFAULT_MUTATION_RATE
THRESHOLDS_YEARS = tuple(int(value) for value in gamma_decode.TMRCA_THRESHOLDS_YEARS)
PAIR_CLASS_STATISTICS = ("p_hom_ref", "p_heterozygous", "p_hom_alt")
PRIMARY_STATISTICS = ("p_hom_alt_minus_hom_ref",)
SECONDARY_STATISTICS = ("p_hom_ref", "p_hom_alt")
PRIMARY_ROLE = "primary_shared_homozygote_bank"
SECONDARY_ROLE = "secondary_class_specific_bank"
EPSILON = 1e-14
RECONCILIATION_TOLERANCE = 5e-5
PNG_WIDTH_PIXELS = 3300
PNG_HEIGHT_PIXELS = 2550
PDF_WIDTH_POINTS = 792.0
PDF_HEIGHT_POINTS = 612.0
FIGURE_STEM = "cosi2_gamma_conditioned_focal_pair_profiles_and_pvalues"
RESULTS_COMPLETION_FILENAME = "RESULTS_COMPLETION.json"
RESULT_OUTPUTS = {
    "scores": "focal_pair_unit_scores.tsv",
    "eligibility": "focal_pair_eligibility.tsv",
    "input_manifest": "focal_pair_input_manifest.tsv",
    "reconciliation": "focal_pair_overall_reconciliation.tsv",
    "pointwise": "focal_pair_simulation_pvalues.tsv",
    "omnibus": "focal_pair_minp_omnibus.tsv",
    "run_results": "RUN_RESULTS.md",
    "figure_png": f"{FIGURE_STEM}.png",
    "figure_pdf": f"{FIGURE_STEM}.pdf",
}
SCORE_COLUMNS = (
    "schema",
    "unit_id",
    "demography",
    "simulation_class",
    "selection_coefficient",
    "seed",
    "generation_time_years",
    "source",
    "pairing_unit",
    "statistic",
    "threshold_years",
    "score",
    "n_pairs",
    "n_valid_pairs",
    "n_ref_pairs",
    "n_alt_pairs",
    "n_valid_ref_pairs",
    "n_valid_alt_pairs",
)
RECONCILIATION_COLUMNS = (
    "schema",
    "unit_id",
    "demography",
    "simulation_class",
    "selection_coefficient",
    "seed",
    "generation_time_years",
    "threshold_years",
    "posthoc_overall_score",
    "decoder_cpp_overall_score",
    "signed_difference",
    "absolute_difference",
)


def sha256_file(path: str | Path) -> str:
    """Return a streaming SHA-256 digest."""

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


def _portable_frame_sha256(frame: pd.DataFrame) -> str:
    payload = frame.to_csv(
        sep="\t", index=False, lineterminator="\n", float_format="%.17g"
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _read_json(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"{label} is absent or unreadable: {path}") from error
    if not isinstance(value, dict):
        raise TypeError(f"{label} must be a JSON object: {path}")
    return value


def _atomic_json(path: Path, value: Any) -> None:
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    try:
        temporary.write_text(
            json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_text(path: Path, value: str) -> None:
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    try:
        temporary.write_text(value, encoding="utf-8", newline="\n")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_frame(path: Path, frame: pd.DataFrame) -> None:
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    try:
        frame.to_csv(
            temporary,
            sep="\t",
            index=False,
            lineterminator="\n",
            float_format="%.17g",
        )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _relative(path: Path, root: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(root.resolve()).as_posix()
    except ValueError as error:
        raise ValueError(f"input lies outside the repository: {resolved}") from error


def _output_record(path: Path, *, rows: int | None = None) -> dict[str, Any]:
    record: dict[str, Any] = {
        "path": path.name,
        "sha256": sha256_file(path),
        "size_bytes": path.stat().st_size,
    }
    if rows is not None:
        record["rows"] = int(rows)
    return record


def _validate_record(path: Path, record: Mapping[str, Any], label: str) -> None:
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"{label} is absent or not a regular file: {path}")
    if path.stat().st_size != int(record.get("size_bytes", -1)) or sha256_file(
        path
    ) != record.get("sha256"):
        raise ValueError(f"{label} checksum or size differs: {path}")


def _data_pair_rows(path: Path) -> np.ndarray:
    rows: list[tuple[int, int]] = []
    with path.open("r", encoding="ascii") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split("\t")
            if len(fields) not in (2, 4):
                raise ValueError(f"invalid pair row {line_number}: {path}")
            try:
                left, right = int(fields[0]), int(fields[1])
            except ValueError as error:
                raise ValueError(
                    f"noninteger pair row {line_number}: {path}"
                ) from error
            rows.append((left, right))
    result = np.asarray(rows, dtype=np.int64)
    if result.shape != (gamma_decode.DEFAULT_PAIR_COUNT, 2):
        raise ValueError(f"pair manifest has unexpected geometry: {path}")
    return result


def _validated_unit_inputs(
    repo_root: Path,
    unit_dir: Path,
    *,
    expected_decode_sha256: str,
) -> tuple[dict[str, Any], dict[str, Any], list[dict[str, Any]]]:
    decoded = unit_dir / "decoded"
    converted = unit_dir / "converted"
    decode_completion_path = decoded / "decode_complete.json"
    if sha256_file(decode_completion_path) != expected_decode_sha256:
        raise ValueError(f"decode completion differs for {unit_dir.name}")
    decode_completion = _read_json(decode_completion_path, "decode completion")
    if (
        decode_completion.get("schema") != gamma_decode.DECODE_SCHEMA
        or decode_completion.get("status") != "complete"
    ):
        raise ValueError(f"decode completion contract is invalid: {unit_dir.name}")

    manifest: list[dict[str, Any]] = []

    def register(label: str, path: Path, record: Mapping[str, Any]) -> None:
        _validate_record(path, record, f"{unit_dir.name} {label}")
        manifest.append(
            {
                "unit_id": unit_dir.name,
                "input_label": label,
                "path": _relative(path, repo_root),
                "sha256": record["sha256"],
                "size_bytes": int(record["size_bytes"]),
            }
        )

    manifest.append(
        {
            "unit_id": unit_dir.name,
            "input_label": "decode_completion",
            "path": _relative(decode_completion_path, repo_root),
            "sha256": expected_decode_sha256,
            "size_bytes": decode_completion_path.stat().st_size,
        }
    )
    decode_outputs = decode_completion.get("outputs", {})
    decoded_paths = {
        "posterior": decoded / "posterior.zst",
        "posterior_metadata": decoded / "posterior.zst.meta",
        "pair_manifest": decoded / "decoded_pairs.tsv",
    }
    for label, path in decoded_paths.items():
        record = decode_outputs.get(label)
        if not isinstance(record, dict):
            raise TypeError(f"decode output record is absent: {unit_dir.name} {label}")
        register(label, path, record)

    conversion_completion_path = converted / "conversion_complete.json"
    expected_conversion_sha = (
        decode_completion.get("contract", {})
        .get("conversion", {})
        .get("completion_sha256")
    )
    if (
        not isinstance(expected_conversion_sha, str)
        or sha256_file(conversion_completion_path) != expected_conversion_sha
    ):
        raise ValueError(f"conversion completion differs for {unit_dir.name}")
    conversion_completion = _read_json(
        conversion_completion_path, "conversion completion"
    )
    manifest.append(
        {
            "unit_id": unit_dir.name,
            "input_label": "conversion_completion",
            "path": _relative(conversion_completion_path, repo_root),
            "sha256": expected_conversion_sha,
            "size_bytes": conversion_completion_path.stat().st_size,
        }
    )
    conversion_outputs = conversion_completion.get("outputs", {})
    converted_paths = {
        "conversion_metadata": converted / "conversion_metadata.json",
        "vcf": converted / "input.vcf",
        "requested_pairs": converted / "overall.pairs.tsv",
    }
    conversion_labels = {
        "conversion_metadata": "metadata",
        "vcf": "vcf",
        "requested_pairs": "pairs",
    }
    for label, path in converted_paths.items():
        record = conversion_outputs.get(conversion_labels[label])
        if not isinstance(record, dict):
            raise TypeError(
                f"conversion output record is absent: {unit_dir.name} {label}"
            )
        register(label, path, record)
    return decode_completion, conversion_completion, manifest


def _focal_alleles(
    unit_dir: Path,
    conversion_metadata: Mapping[str, Any],
    *,
    recorded_status: str,
) -> tuple[np.ndarray | None, str]:
    site_count = int(conversion_metadata.get("exact_focal_site_count", -1))
    if site_count > 1:
        if recorded_status != "ambiguous_exact_coordinate_collision":
            raise ValueError(f"focal collision status differs: {unit_dir.name}")
        return None, recorded_status
    if site_count == 0:
        if recorded_status != "absent_from_sample":
            raise ValueError(f"absent focal status differs: {unit_dir.name}")
        return np.zeros(
            gamma_decode.DEFAULT_SAMPLE_HAPLOTYPES, dtype=np.int8
        ), recorded_status
    if site_count != 1 or recorded_status != "unique_exact_focal_site":
        raise ValueError(f"unique focal status differs: {unit_dir.name}")

    site_index = conversion_metadata.get("exact_focal_site_index")
    if site_index is None or int(site_index) < 0:
        raise ValueError(f"exact focal site index is absent: {unit_dir.name}")
    expected_id = f"cosi_s{int(site_index):08d}"
    matches: list[list[str]] = []
    expected_position = FOCAL_POSITION_0BASED + 1
    with (unit_dir / "converted" / "input.vcf").open("r", encoding="ascii") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3 and fields[2] == expected_id:
                matches.append(fields)
    if len(matches) != 1:
        raise ValueError(f"exact focal VCF record ID is not unique: {unit_dir.name}")
    if matches[0][0] != "1" or int(matches[0][1]) != expected_position:
        raise ValueError(f"exact focal VCF record position differs: {unit_dir.name}")
    genotypes = matches[0][9:]
    alleles: list[int] = []
    for genotype in genotypes:
        fields = genotype.split(":", maxsplit=1)[0].split("|")
        if len(fields) != 2 or any(value not in {"0", "1"} for value in fields):
            raise ValueError(f"focal genotype is not phased biallelic: {unit_dir.name}")
        alleles.extend(int(value) for value in fields)
    result = np.asarray(alleles, dtype=np.int8)
    if result.shape != (gamma_decode.DEFAULT_SAMPLE_HAPLOTYPES,):
        raise ValueError(f"focal sample geometry differs: {unit_dir.name}")
    expected = conversion_metadata.get("exact_focal_sample_alt_count")
    if expected is None or int(result.sum()) != int(expected):
        raise ValueError(f"focal sample ALT count differs: {unit_dir.name}")
    return result, recorded_status


def _read_exact(reader: Any, size: int) -> bytes:
    chunks: list[bytes] = []
    remaining = size
    while remaining:
        block = reader.read(remaining)
        if not block:
            raise EOFError(f"posterior stream ended {remaining} bytes early")
        chunks.append(block)
        remaining -= len(block)
    return b"".join(chunks)


def _stream_focal_scores(
    posterior_path: Path,
    metadata: Mapping[str, Any],
    pair_classes: np.ndarray,
    *,
    generation_time_years: float,
) -> tuple[dict[str, np.ndarray], dict[str, int], np.ndarray]:
    n_pairs = int(metadata["num_pairs"])
    chunk_size = int(metadata["chunk_size"])
    n_positions = int(metadata["sequence_length"])
    positions = np.asarray(metadata["output_positions"], dtype=np.int64)
    focal_matches = np.flatnonzero(positions == FOCAL_POSITION_0BASED)
    if len(focal_matches) != 1 or n_pairs != len(pair_classes):
        raise ValueError("posterior focal position or pair geometry differs")
    focal_index = int(focal_matches[0])
    two_ne_generations = float(metadata["scaled_mutation_rate"]) / (2.0 * MUTATION_RATE)
    scaled_thresholds = (
        np.asarray(THRESHOLDS_YEARS, dtype=np.float64)
        / float(generation_time_years)
        / two_ne_generations
    )
    labels = {0: "p_hom_ref", 1: "p_heterozygous", 2: "p_hom_alt"}
    sums = {label: np.zeros(len(THRESHOLDS_YEARS)) for label in labels.values()}
    counts = {label: 0 for label in labels.values()}
    overall_sum = np.zeros(len(THRESHOLDS_YEARS))
    overall_count = 0
    block_bytes = n_positions * chunk_size * np.dtype(np.float32).itemsize
    with (
        posterior_path.open("rb") as compressed,
        zstandard.ZstdDecompressor().stream_reader(compressed) as reader,
    ):
        for pair_start in range(0, n_pairs, chunk_size):
            real_pairs = min(chunk_size, n_pairs - pair_start)
            alpha = np.frombuffer(
                _read_exact(reader, block_bytes), dtype=np.float32
            ).reshape(n_positions, chunk_size)[focal_index, :real_pairs]
            beta = np.frombuffer(
                _read_exact(reader, block_bytes), dtype=np.float32
            ).reshape(n_positions, chunk_size)[focal_index, :real_pairs]
            valid = np.isfinite(alpha) & np.isfinite(beta) & (alpha > 0) & (beta > 0)
            if not np.any(valid):
                continue
            alpha64 = alpha[valid].astype(np.float64)
            beta64 = beta[valid].astype(np.float64)
            probabilities = np.stack(
                [
                    gammainc(alpha64, beta64 * threshold)
                    for threshold in scaled_thresholds
                ]
            )
            overall_sum += probabilities.sum(axis=1)
            overall_count += int(valid.sum())
            chunk_classes = pair_classes[pair_start : pair_start + real_pairs][valid]
            for code, label in labels.items():
                mask = chunk_classes == code
                if np.any(mask):
                    sums[label] += probabilities[:, mask].sum(axis=1)
                    counts[label] += int(mask.sum())
        if reader.read(1):
            raise ValueError("posterior stream contains trailing decompressed bytes")
    if overall_count != n_pairs:
        raise ValueError("focal posterior contains invalid pair cells")
    means = {
        label: values / counts[label]
        for label, values in sums.items()
        if counts[label] > 0
    }
    return means, counts, overall_sum / overall_count


def _validate_fixed_pair_order(
    pairs: np.ndarray, panel: np.ndarray, *, unit_id: str
) -> None:
    """Require the decoder's exact C(100, 2) panel-order enumeration."""

    if panel.shape != (gamma_decode.DEFAULT_PAIR_PANEL_HAPLOTYPES,) or len(
        np.unique(panel)
    ) != len(panel):
        raise ValueError(f"fixed pair panel differs: {unit_id}")
    expected = np.asarray(
        list(itertools.combinations(panel.tolist(), 2)), dtype=np.int64
    )
    if pairs.shape != expected.shape or not np.array_equal(pairs, expected):
        raise ValueError(f"fixed all-pairs panel order differs: {unit_id}")


def _validate_focal_alt_counts(
    alleles: np.ndarray,
    panel: np.ndarray,
    *,
    status: str,
    conversion_metadata: Mapping[str, Any],
    unit_record: Mapping[str, Any],
    unit_id: str,
) -> tuple[int, int]:
    """Cross-check focal ALT counts across VCF, conversion, and study metadata."""

    sample_alt = int(alleles.sum())
    panel_alt = int(alleles[panel].sum())
    recorded_sample_alt = int(unit_record["decode_full_sample_alt_count"])
    recorded_panel_alt = int(unit_record["decode_pair_panel_alt_count"])
    if status == "unique_exact_focal_site":
        expected_sample = conversion_metadata.get("exact_focal_sample_alt_count")
        expected_panel = conversion_metadata.get("exact_focal_pair_panel_alt_count")
        if expected_sample is None or sample_alt != int(expected_sample):
            raise ValueError(f"full-sample ALT count differs: {unit_id}")
        if expected_panel is None or panel_alt != int(expected_panel):
            raise ValueError(f"pair-panel ALT count differs: {unit_id}")
    elif status == "absent_from_sample":
        if sample_alt != 0 or panel_alt != 0:
            raise ValueError(f"absent focal allele counts are nonzero: {unit_id}")
        for key in (
            "exact_focal_sample_alt_count",
            "exact_focal_pair_panel_alt_count",
        ):
            value = conversion_metadata.get(key)
            if value is not None and int(value) != 0:
                raise ValueError(f"absent focal conversion count differs: {unit_id}")
    else:
        raise ValueError(f"unsupported count-bearing focal status: {unit_id}")
    if recorded_sample_alt != sample_alt:
        raise ValueError(f"recorded full-sample ALT count differs: {unit_id}")
    if recorded_panel_alt != panel_alt:
        raise ValueError(f"recorded pair-panel ALT count differs: {unit_id}")
    return sample_alt, panel_alt


def reduce_focal_unit(
    repo_root: str | Path,
    unit_record: Mapping[str, Any],
    *,
    unit_dir: str | Path,
    expected_decode_sha256: str,
) -> tuple[pd.DataFrame, dict[str, Any], pd.DataFrame, pd.DataFrame]:
    """Reduce one retained posterior to focal pair-class probabilities."""

    root = Path(repo_root).resolve()
    directory = Path(unit_dir).resolve()
    unit_id = str(unit_record["unit_id"])
    if directory.name != unit_id:
        raise ValueError("decode unit directory and unit ID disagree")
    decode_completion, _, manifest_rows = _validated_unit_inputs(
        root, directory, expected_decode_sha256=expected_decode_sha256
    )
    conversion_metadata = _read_json(
        directory / "converted" / "conversion_metadata.json", "conversion metadata"
    )
    recorded_status = str(unit_record["decode_focal_count_status"])
    alleles, status = _focal_alleles(
        directory, conversion_metadata, recorded_status=recorded_status
    )
    metadata = _read_json(
        directory / "decoded" / "posterior.zst.meta", "posterior metadata"
    )
    pairs = np.asarray(metadata.get("pairs", []), dtype=np.int64)
    decoded_pairs = _data_pair_rows(directory / "decoded" / "decoded_pairs.tsv")
    requested_pairs = _data_pair_rows(directory / "converted" / "overall.pairs.tsv")
    if (
        pairs.shape != decoded_pairs.shape
        or not np.array_equal(pairs, decoded_pairs)
        or not np.array_equal(pairs, requested_pairs)
    ):
        raise ValueError(f"posterior pair order differs: {unit_id}")
    panel = np.asarray(
        conversion_metadata.get("pair_panel_indices", []), dtype=np.int64
    )
    _validate_fixed_pair_order(pairs, panel, unit_id=unit_id)

    metadata_columns = {
        "unit_id": unit_id,
        "demography": str(unit_record["demography"]),
        "simulation_class": str(unit_record["simulation_class"]),
        "selection_coefficient": float(unit_record["selection_coefficient"]),
        "seed": int(unit_record["seed"]),
        "generation_time_years": float(unit_record["generation_time_years"]),
    }
    eligibility: dict[str, Any] = {
        **metadata_columns,
        "focal_status": status,
        "focal_alt_count_full_sample": -1,
        "focal_alt_count_pair_panel": -1,
        "n_hom_ref_pairs": 0,
        "n_heterozygous_pairs": 0,
        "n_hom_alt_pairs": 0,
        "hom_ref_eligible": False,
        "heterozygous_eligible": False,
        "hom_alt_eligible": False,
        "primary_shared_homozygote_eligible": False,
        "exclusion_reason": "ambiguous_exact_coordinate_collision",
    }
    if alleles is None:
        return (
            pd.DataFrame(),
            eligibility,
            pd.DataFrame(),
            pd.DataFrame(manifest_rows),
        )

    pair_alleles = alleles[pairs]
    allele_sums = pair_alleles.sum(axis=1)
    if np.any((allele_sums < 0) | (allele_sums > 2)):
        raise ValueError(f"pair allele class is invalid: {unit_id}")
    pair_classes = allele_sums.astype(np.int8)
    sample_alt, panel_alt = _validate_focal_alt_counts(
        alleles,
        panel,
        status=status,
        conversion_metadata=conversion_metadata,
        unit_record=unit_record,
        unit_id=unit_id,
    )
    expected_pair_counts = {
        "p_hom_ref": math.comb(len(panel) - panel_alt, 2),
        "p_heterozygous": panel_alt * (len(panel) - panel_alt),
        "p_hom_alt": math.comb(panel_alt, 2),
    }
    observed_pair_counts = {
        "p_hom_ref": int(np.count_nonzero(pair_classes == 0)),
        "p_heterozygous": int(np.count_nonzero(pair_classes == 1)),
        "p_hom_alt": int(np.count_nonzero(pair_classes == 2)),
    }
    if observed_pair_counts != expected_pair_counts:
        raise ValueError(f"pair-class combinatorics differ: {unit_id}")

    means, valid_counts, overall = _stream_focal_scores(
        directory / "decoded" / "posterior.zst",
        metadata,
        pair_classes,
        generation_time_years=float(unit_record["generation_time_years"]),
    )
    rows: list[dict[str, Any]] = []
    for statistic, values in means.items():
        is_ref = statistic == "p_hom_ref"
        is_alt = statistic == "p_hom_alt"
        for threshold, score in zip(THRESHOLDS_YEARS, values, strict=True):
            rows.append(
                {
                    "schema": SCHEMA_VERSION,
                    **metadata_columns,
                    "source": "retained_float32_gamma_posterior_posthoc",
                    "pairing_unit": "unordered_haplotype_pair_pseudodiploid",
                    "statistic": statistic,
                    "threshold_years": threshold,
                    "score": float(score),
                    "n_pairs": observed_pair_counts[statistic],
                    "n_valid_pairs": valid_counts[statistic],
                    "n_ref_pairs": (observed_pair_counts["p_hom_ref"] if is_ref else 0),
                    "n_alt_pairs": (observed_pair_counts["p_hom_alt"] if is_alt else 0),
                    "n_valid_ref_pairs": (valid_counts["p_hom_ref"] if is_ref else 0),
                    "n_valid_alt_pairs": (valid_counts["p_hom_alt"] if is_alt else 0),
                }
            )
    if "p_hom_ref" in means and "p_hom_alt" in means:
        contrast = means["p_hom_alt"] - means["p_hom_ref"]
        for threshold, score in zip(THRESHOLDS_YEARS, contrast, strict=True):
            rows.append(
                {
                    "schema": SCHEMA_VERSION,
                    **metadata_columns,
                    "source": "retained_float32_gamma_posterior_posthoc",
                    "pairing_unit": "unordered_haplotype_pair_pseudodiploid",
                    "statistic": "p_hom_alt_minus_hom_ref",
                    "threshold_years": threshold,
                    "score": float(score),
                    "n_pairs": -1,
                    "n_valid_pairs": -1,
                    "n_ref_pairs": observed_pair_counts["p_hom_ref"],
                    "n_alt_pairs": observed_pair_counts["p_hom_alt"],
                    "n_valid_ref_pairs": valid_counts["p_hom_ref"],
                    "n_valid_alt_pairs": valid_counts["p_hom_alt"],
                }
            )
    eligibility.update(
        {
            "focal_alt_count_full_sample": sample_alt,
            "focal_alt_count_pair_panel": panel_alt,
            "n_hom_ref_pairs": observed_pair_counts["p_hom_ref"],
            "n_heterozygous_pairs": observed_pair_counts["p_heterozygous"],
            "n_hom_alt_pairs": observed_pair_counts["p_hom_alt"],
            "hom_ref_eligible": observed_pair_counts["p_hom_ref"] > 0,
            "heterozygous_eligible": observed_pair_counts["p_heterozygous"] > 0,
            "hom_alt_eligible": observed_pair_counts["p_hom_alt"] > 0,
            "primary_shared_homozygote_eligible": (
                observed_pair_counts["p_hom_ref"] > 0
                and observed_pair_counts["p_hom_alt"] > 0
            ),
            "exclusion_reason": ""
            if observed_pair_counts["p_hom_alt"] > 0
            else ("fewer_than_two_alt_haplotypes_in_fixed_pair_panel"),
        }
    )
    reconciliation = pd.DataFrame(
        {
            "schema": SCHEMA_VERSION,
            **{
                key: [value] * len(THRESHOLDS_YEARS)
                for key, value in metadata_columns.items()
            },
            "threshold_years": THRESHOLDS_YEARS,
            "posthoc_overall_score": overall,
        }
    )
    # Bind the reducer to the same decode contract even though the result table
    # is compared to the persisted conditioned aggregate in the parent process.
    if decode_completion.get("contract", {}).get("unit_id") != unit_id:
        raise ValueError(f"decode contract unit ID differs: {unit_id}")
    return pd.DataFrame(rows), eligibility, reconciliation, pd.DataFrame(manifest_rows)


def validate_scores(scores: pd.DataFrame) -> pd.DataFrame:
    required = set(SCORE_COLUMNS)
    missing = sorted(required.difference(scores.columns))
    if missing:
        raise ValueError("focal score columns are absent: " + ", ".join(missing))
    result = scores.loc[:, list(SCORE_COLUMNS)].copy()
    if result.empty or set(result["schema"]) != {SCHEMA_VERSION}:
        raise ValueError("focal score schema differs")
    allowed = {*PAIR_CLASS_STATISTICS, "p_hom_alt_minus_hom_ref"}
    if set(result["statistic"]) - allowed:
        raise ValueError("unknown focal pair statistic")
    result["threshold_years"] = pd.to_numeric(
        result["threshold_years"], errors="raise"
    ).astype(int)
    if set(result["threshold_years"]) != set(THRESHOLDS_YEARS):
        raise ValueError("focal score thresholds differ")
    result["score"] = pd.to_numeric(result["score"], errors="raise")
    if not np.all(np.isfinite(result["score"])):
        raise ValueError("focal scores must be finite")
    direct = result["statistic"].isin(PAIR_CLASS_STATISTICS)
    if ((result.loc[direct, "score"] < 0) | (result.loc[direct, "score"] > 1)).any():
        raise ValueError("focal probabilities lie outside [0,1]")
    if ((result.loc[~direct, "score"] < -1) | (result.loc[~direct, "score"] > 1)).any():
        raise ValueError("focal contrasts lie outside [-1,1]")
    count_columns = (
        "n_pairs",
        "n_valid_pairs",
        "n_ref_pairs",
        "n_alt_pairs",
        "n_valid_ref_pairs",
        "n_valid_alt_pairs",
    )
    for column in count_columns:
        values = pd.to_numeric(result[column], errors="raise")
        if not np.all(np.isfinite(values)) or not np.all(values == np.floor(values)):
            raise ValueError(f"focal pair count is not a finite integer: {column}")
        result[column] = values.astype(int)
    if (
        (result.loc[direct, "n_pairs"] <= 0).any()
        or (result.loc[direct, "n_valid_pairs"] != result.loc[direct, "n_pairs"]).any()
        or (result.loc[direct, "n_ref_pairs"] < 0).any()
        or (result.loc[direct, "n_alt_pairs"] < 0).any()
    ):
        raise ValueError("direct focal pair denominators are invalid")
    contrast = ~direct
    if (
        (result.loc[contrast, "n_pairs"] != -1).any()
        or (result.loc[contrast, "n_valid_pairs"] != -1).any()
        or (result.loc[contrast, "n_ref_pairs"] <= 0).any()
        or (result.loc[contrast, "n_alt_pairs"] <= 0).any()
        or (
            result.loc[contrast, "n_valid_ref_pairs"]
            != result.loc[contrast, "n_ref_pairs"]
        ).any()
        or (
            result.loc[contrast, "n_valid_alt_pairs"]
            != result.loc[contrast, "n_alt_pairs"]
        ).any()
    ):
        raise ValueError("contrast focal pair denominators are invalid")
    direct_denominators = {
        "p_hom_ref": ("n_ref_pairs", "n_valid_ref_pairs"),
        "p_hom_alt": ("n_alt_pairs", "n_valid_alt_pairs"),
    }
    for statistic, (pair_column, valid_column) in direct_denominators.items():
        mask = result["statistic"] == statistic
        other_pair = "n_alt_pairs" if statistic == "p_hom_ref" else "n_ref_pairs"
        other_valid = (
            "n_valid_alt_pairs" if statistic == "p_hom_ref" else "n_valid_ref_pairs"
        )
        if (
            (result.loc[mask, pair_column] != result.loc[mask, "n_pairs"]).any()
            or (
                result.loc[mask, valid_column] != result.loc[mask, "n_valid_pairs"]
            ).any()
            or (result.loc[mask, other_pair] != 0).any()
            or (result.loc[mask, other_valid] != 0).any()
        ):
            raise ValueError(f"{statistic} explicit denominators differ")
    heterozygous = result["statistic"] == "p_heterozygous"
    if (
        result.loc[
            heterozygous,
            [
                "n_ref_pairs",
                "n_alt_pairs",
                "n_valid_ref_pairs",
                "n_valid_alt_pairs",
            ],
        ]
        != 0
    ).any(axis=None):
        raise ValueError("heterozygous rows must not claim homozygous denominators")
    key = ["unit_id", "statistic", "threshold_years"]
    if result.duplicated(key).any():
        raise ValueError("duplicate focal score rows")
    for _, group in result.groupby(["unit_id", "statistic"], sort=False):
        if tuple(sorted(group["threshold_years"])) != THRESHOLDS_YEARS:
            raise ValueError("focal score curve is incomplete")
        if group.iloc[0]["statistic"] in PAIR_CLASS_STATISTICS:
            curve = group.sort_values("threshold_years")["score"].to_numpy(float)
            if np.any(np.diff(curve) < -1e-12):
                raise ValueError("focal probability curve is not monotone")
    return result.sort_values(key, kind="mergesort").reset_index(drop=True)


def _expected_reconciliation_unit_ids(eligibility: pd.DataFrame) -> tuple[str, ...]:
    required = {
        "unit_id",
        "demography",
        "simulation_class",
        "selection_coefficient",
        "seed",
        "generation_time_years",
        "focal_status",
    }
    missing = sorted(required.difference(eligibility.columns))
    if missing:
        raise ValueError("focal eligibility columns are absent: " + ", ".join(missing))
    if eligibility.empty or eligibility["unit_id"].astype(str).duplicated().any():
        raise ValueError("focal eligibility unit inventory is empty or duplicated")
    statuses = set(eligibility["focal_status"].astype(str))
    allowed_statuses = {
        "unique_exact_focal_site",
        "absent_from_sample",
        "ambiguous_exact_coordinate_collision",
    }
    if not statuses or statuses - allowed_statuses:
        raise ValueError("focal eligibility contains an unknown status")
    unit_ids = tuple(
        sorted(
            eligibility.loc[
                eligibility["focal_status"].astype(str)
                != "ambiguous_exact_coordinate_collision",
                "unit_id",
            ].astype(str)
        )
    )
    if not unit_ids:
        raise ValueError("no nonambiguous focal units are eligible for reconciliation")
    return unit_ids


def validate_reconciliation(
    reconciliation: pd.DataFrame, eligibility: pd.DataFrame
) -> pd.DataFrame:
    """Require complete, finite, algebraically consistent overall-score checks."""

    if set(reconciliation.columns) != set(RECONCILIATION_COLUMNS):
        raise ValueError("focal reconciliation columns differ")
    result = reconciliation.loc[:, list(RECONCILIATION_COLUMNS)].copy()
    if result.empty or set(result["schema"].astype(str)) != {SCHEMA_VERSION}:
        raise ValueError("focal reconciliation schema differs")
    expected_ids = _expected_reconciliation_unit_ids(eligibility)
    numeric_columns = (
        "selection_coefficient",
        "seed",
        "generation_time_years",
        "threshold_years",
        "posthoc_overall_score",
        "decoder_cpp_overall_score",
        "signed_difference",
        "absolute_difference",
    )
    for column in numeric_columns:
        result[column] = pd.to_numeric(result[column], errors="raise")
        if not np.all(np.isfinite(result[column])):
            raise ValueError(f"focal reconciliation contains nonfinite {column}")
    for column in ("seed", "threshold_years"):
        if not np.all(result[column] == np.floor(result[column])):
            raise ValueError(f"focal reconciliation contains noninteger {column}")
        result[column] = result[column].astype(int)
    key = ["unit_id", "threshold_years"]
    if result.duplicated(key).any():
        raise ValueError("duplicate focal reconciliation rows")
    result["unit_id"] = result["unit_id"].astype(str)
    expected_index = pd.MultiIndex.from_product(
        [expected_ids, THRESHOLDS_YEARS], names=key
    )
    observed_index = pd.MultiIndex.from_frame(result.loc[:, key])
    if len(result) != len(expected_index) or set(observed_index) != set(expected_index):
        raise ValueError("focal reconciliation unit-threshold coverage differs")

    eligibility_metadata = eligibility.loc[
        eligibility["unit_id"].astype(str).isin(expected_ids),
        [
            "unit_id",
            "demography",
            "simulation_class",
            "selection_coefficient",
            "seed",
            "generation_time_years",
        ],
    ].copy()
    eligibility_metadata["unit_id"] = eligibility_metadata["unit_id"].astype(str)
    comparison = result.merge(
        eligibility_metadata,
        on="unit_id",
        how="left",
        suffixes=("", "_expected"),
        validate="many_to_one",
    )
    for column in ("demography", "simulation_class"):
        if not (
            comparison[column].astype(str)
            == comparison[f"{column}_expected"].astype(str)
        ).all():
            raise ValueError(f"focal reconciliation {column} differs from eligibility")
    for column in (
        "selection_coefficient",
        "seed",
        "generation_time_years",
    ):
        expected = pd.to_numeric(comparison[f"{column}_expected"], errors="raise")
        if not np.array_equal(comparison[column].to_numpy(), expected.to_numpy()):
            raise ValueError(f"focal reconciliation {column} differs from eligibility")

    computed_signed = (
        result["posthoc_overall_score"] - result["decoder_cpp_overall_score"]
    ).to_numpy(float)
    signed = result["signed_difference"].to_numpy(float)
    absolute = result["absolute_difference"].to_numpy(float)
    if not np.allclose(signed, computed_signed, rtol=0, atol=1e-15):
        raise ValueError("focal reconciliation signed differences are inconsistent")
    if not np.allclose(absolute, np.abs(signed), rtol=0, atol=1e-15):
        raise ValueError("focal reconciliation absolute differences are inconsistent")
    if float(absolute.max()) > RECONCILIATION_TOLERANCE:
        raise ValueError("posthoc/C++ overall reconciliation exceeds tolerance")
    return result.sort_values(key, kind="mergesort").reset_index(drop=True)


def _tail_statistics(observed: float, null: np.ndarray) -> dict[str, Any]:
    if not np.isfinite(observed) or not len(null) or not np.all(np.isfinite(null)):
        raise ValueError("simulation-reference scores must be finite and nonempty")
    upper_count = int(np.count_nonzero(null >= observed - EPSILON))
    lower_count = int(np.count_nonzero(null <= observed + EPSILON))
    upper = (1 + upper_count) / (1 + len(null))
    lower = (1 + lower_count) / (1 + len(null))
    return {
        "neutral_upper_extreme_count": upper_count,
        "neutral_lower_extreme_count": lower_count,
        "mc_p_upper": float(upper),
        "mc_p_lower": float(lower),
        "mc_p_two_sided": float(min(1.0, 2.0 * min(upper, lower))),
        "neutral_midrank_descriptive_auc": float(
            (
                np.count_nonzero(null < observed - EPSILON)
                + 0.5 * np.count_nonzero(np.abs(null - observed) <= EPSILON)
            )
            / len(null)
        ),
    }


def pointwise_simulation_pvalues(
    scores: pd.DataFrame, eligibility: pd.DataFrame
) -> pd.DataFrame:
    """Return primary shared-bank and secondary class-specific p-values."""

    frame = validate_scores(scores)
    eligible = eligibility.copy()
    rows: list[dict[str, Any]] = []
    configurations = [
        (PRIMARY_ROLE, PRIMARY_STATISTICS, "primary_shared_homozygote_eligible"),
        (SECONDARY_ROLE, SECONDARY_STATISTICS, None),
    ]
    class_columns = {
        "p_hom_ref": "hom_ref_eligible",
        "p_heterozygous": "heterozygous_eligible",
        "p_hom_alt": "hom_alt_eligible",
    }
    for demography in ("EAS", "Han"):
        demography_eligibility = eligible[eligible["demography"] == demography]
        for role, statistics, shared_column in configurations:
            for statistic in statistics:
                eligibility_column = shared_column or class_columns[statistic]
                neutral_ids = set(
                    demography_eligibility.loc[
                        (demography_eligibility["simulation_class"] == "neutral")
                        & demography_eligibility[eligibility_column].astype(bool),
                        "unit_id",
                    ].astype(str)
                )
                selected_ids = set(
                    demography_eligibility.loc[
                        (demography_eligibility["simulation_class"] == "selected")
                        & demography_eligibility[eligibility_column].astype(bool),
                        "unit_id",
                    ].astype(str)
                )
                if len(selected_ids) != 1 or not neutral_ids:
                    raise ValueError(
                        f"incomplete {role} bank for {demography} {statistic}"
                    )
                cell = frame[
                    (frame["demography"] == demography)
                    & (frame["statistic"] == statistic)
                ]
                for threshold in THRESHOLDS_YEARS:
                    threshold_frame = cell[cell["threshold_years"] == threshold]
                    selected = threshold_frame[
                        threshold_frame["unit_id"].astype(str).isin(selected_ids)
                    ]
                    neutral = threshold_frame[
                        threshold_frame["unit_id"].astype(str).isin(neutral_ids)
                    ]
                    if len(selected) != 1 or len(neutral) != len(neutral_ids):
                        raise ValueError(
                            f"incomplete score bank for {demography} {statistic}"
                        )
                    observed = float(selected.iloc[0]["score"])
                    null = neutral["score"].to_numpy(float)
                    rows.append(
                        {
                            "schema": SCHEMA_VERSION,
                            "analysis_role": role,
                            "demography": demography,
                            "statistic": statistic,
                            "threshold_years": threshold,
                            "observed_selected": observed,
                            "n_neutral": len(null),
                            "minimum_attainable_p": 1 / (1 + len(null)),
                            "neutral_mean": float(np.mean(null)),
                            "neutral_median": float(np.median(null)),
                            "neutral_q025": float(np.quantile(null, 0.025)),
                            "neutral_q975": float(np.quantile(null, 0.975)),
                            **_tail_statistics(observed, null),
                            "tail_contract": "plus_one_conservative_ties",
                            "null_bank": (
                                "shared_complete_case_ref_ref_and_alt_alt"
                                if role == PRIMARY_ROLE
                                else f"maximal_eligible_{statistic.removeprefix('p_')}_bank"
                            ),
                            "pairing_unit": ("unordered_haplotype_pair_pseudodiploid"),
                        }
                    )
    return (
        pd.DataFrame(rows)
        .sort_values(
            ["analysis_role", "demography", "statistic", "threshold_years"],
            kind="mergesort",
        )
        .reset_index(drop=True)
    )


def minp_omnibus(scores: pd.DataFrame, eligibility: pd.DataFrame) -> pd.DataFrame:
    """Calibrate seven-threshold primary minP using whole neutral vectors."""

    frame = validate_scores(scores)
    rows: list[dict[str, Any]] = []
    for demography in ("EAS", "Han"):
        eligible = eligibility[
            (eligibility["demography"] == demography)
            & eligibility["primary_shared_homozygote_eligible"].astype(bool)
        ]
        neutral_ids = set(
            eligible.loc[eligible["simulation_class"] == "neutral", "unit_id"].astype(
                str
            )
        )
        selected_ids = set(
            eligible.loc[eligible["simulation_class"] == "selected", "unit_id"].astype(
                str
            )
        )
        if len(selected_ids) != 1 or not neutral_ids:
            raise ValueError(f"incomplete primary minP bank for {demography}")
        for statistic in PRIMARY_STATISTICS:
            cell = frame[
                (frame["demography"] == demography) & (frame["statistic"] == statistic)
            ]
            profiles = cell.pivot(
                index="unit_id", columns="threshold_years", values="score"
            ).reindex(columns=THRESHOLDS_YEARS)
            neutral = profiles.loc[sorted(neutral_ids)].to_numpy(float)
            observed = profiles.loc[next(iter(selected_ids))].to_numpy(float)
            inference = gamma_analysis._exchangeable_minp(neutral, observed)
            output: dict[str, Any] = {
                "schema": SCHEMA_VERSION,
                "analysis_role": PRIMARY_ROLE,
                "demography": demography,
                "statistic": statistic,
                "n_neutral": len(neutral),
                "n_thresholds": len(THRESHOLDS_YEARS),
                "thresholds_years": ",".join(map(str, THRESHOLDS_YEARS)),
                "minimum_attainable_p": 1 / (1 + len(neutral)),
                "omnibus_method": (
                    "exchangeable_whole_unit_minP_over_seven_thresholds"
                ),
            }
            for tail in ("upper", "lower", "two_sided"):
                index = int(inference[f"threshold_index_at_min_pointwise_p_{tail}"])
                output[f"observed_min_pointwise_p_{tail}"] = float(
                    inference[f"observed_min_pointwise_p_{tail}"]
                )
                output[f"threshold_years_at_min_pointwise_p_{tail}"] = THRESHOLDS_YEARS[
                    index
                ]
                output[f"observed_score_at_min_pointwise_p_{tail}"] = float(
                    observed[index]
                )
                output[f"neutral_minp_as_or_more_extreme_count_{tail}"] = int(
                    inference[f"neutral_minp_as_or_more_extreme_count_{tail}"]
                )
                output[f"omnibus_mc_p_{tail}"] = float(
                    inference[f"omnibus_mc_p_{tail}"]
                )
            rows.append(output)
    return (
        pd.DataFrame(rows)
        .sort_values(["demography", "statistic"], kind="mergesort")
        .reset_index(drop=True)
    )


def _atomic_save_figure_pair(
    figure: Figure, output_dir: Path, stem: str
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = {"png": output_dir / f"{stem}.png", "pdf": output_dir / f"{stem}.pdf"}
    for extension, path in outputs.items():
        temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
        try:
            figure.savefig(
                temporary,
                format=extension,
                dpi=300 if extension == "png" else None,
                facecolor="white",
                metadata=(
                    {"Software": "gamma_smc_aou.cosi2_conditioned_focal_pairs"}
                    if extension == "png"
                    else {
                        "Creator": "gamma_smc_aou.cosi2_conditioned_focal_pairs",
                        "CreationDate": None,
                        "ModDate": None,
                    }
                ),
            )
            os.replace(temporary, path)
        finally:
            temporary.unlink(missing_ok=True)
    return outputs


def plot_focal_pair_results(
    scores: pd.DataFrame,
    eligibility: pd.DataFrame,
    pointwise: pd.DataFrame,
    omnibus: pd.DataFrame,
    output_dir: str | Path,
    *,
    stem: str = FIGURE_STEM,
) -> dict[str, Path]:
    """Create the letter-landscape focal profile and p-value figure pair."""

    frame = validate_scores(scores)
    plot_p = pointwise[
        (
            (pointwise["analysis_role"] == SECONDARY_ROLE)
            & pointwise["statistic"].isin(SECONDARY_STATISTICS)
        )
        | (
            (pointwise["analysis_role"] == PRIMARY_ROLE)
            & (pointwise["statistic"] == "p_hom_alt_minus_hom_ref")
        )
    ]
    x = np.asarray(THRESHOLDS_YEARS, dtype=float) / 1_000
    colors = {
        "p_hom_ref": "#595959",
        "p_hom_alt": "#0072B2",
        "p_hom_alt_minus_hom_ref": "#D55E00",
    }
    labels = {
        "p_hom_ref": "REF/REF pair",
        "p_hom_alt": "ALT/ALT pair",
        "p_hom_alt_minus_hom_ref": "ALT/ALT - REF/REF",
    }
    figure = Figure(figsize=(11, 8.5))
    axes = figure.subplots(2, 2, sharex="col", squeeze=False)
    try:
        for column, demography in enumerate(("EAS", "Han")):
            top = axes[0, column]
            bottom = axes[1, column]
            class_ns: dict[str, int] = {}
            for statistic in ("p_hom_ref", "p_hom_alt"):
                eligibility_column = (
                    "hom_ref_eligible"
                    if statistic == "p_hom_ref"
                    else "hom_alt_eligible"
                )
                eligible_ids = set(
                    eligibility.loc[
                        (eligibility["demography"] == demography)
                        & (eligibility["simulation_class"] == "neutral")
                        & eligibility[eligibility_column].astype(bool),
                        "unit_id",
                    ].astype(str)
                )
                class_ns[statistic] = len(eligible_ids)
                cell = frame[
                    (frame["demography"] == demography)
                    & (frame["statistic"] == statistic)
                ]
                neutral = (
                    cell[
                        (cell["simulation_class"] == "neutral")
                        & cell["unit_id"].astype(str).isin(eligible_ids)
                    ]
                    .pivot(index="unit_id", columns="threshold_years", values="score")
                    .reindex(columns=THRESHOLDS_YEARS)
                    .to_numpy(float)
                )
                selected = (
                    cell[cell["simulation_class"] == "selected"]
                    .set_index("threshold_years")
                    .reindex(THRESHOLDS_YEARS)["score"]
                    .to_numpy(float)
                )
                low, median, high = np.quantile(neutral, [0.025, 0.5, 0.975], axis=0)
                color = colors[statistic]
                top.fill_between(x, low, high, color=color, alpha=0.13)
                top.plot(x, median, color=color, linestyle="--", linewidth=2.4)
                top.plot(
                    x,
                    selected,
                    color=color,
                    linewidth=3.0,
                    marker="o" if statistic == "p_hom_alt" else "s",
                    markersize=6,
                )
            selected_eligibility = eligibility[
                (eligibility["demography"] == demography)
                & (eligibility["simulation_class"] == "selected")
            ].iloc[0]
            top.set_title(
                f"{demography}: focal pair-class posterior",
                fontsize=17,
                pad=9,
            )
            top.text(
                0.03,
                0.97,
                (
                    f"neutral n: REF/REF={class_ns['p_hom_ref']}, "
                    f"ALT/ALT={class_ns['p_hom_alt']}\nselected pairs: "
                    f"REF/REF={int(selected_eligibility['n_hom_ref_pairs'])}, "
                    f"ALT/ALT={int(selected_eligibility['n_hom_alt_pairs'])}"
                ),
                transform=top.transAxes,
                va="top",
                fontsize=10.5,
                linespacing=1.25,
                clip_on=True,
                bbox={
                    "boxstyle": "round,pad=0.3",
                    "facecolor": "white",
                    "edgecolor": "#bdbdbd",
                    "alpha": 0.82,
                },
            )
            top.set_ylim(0, 1)
            top.grid(True, color="#d9d9d9", linewidth=0.7)
            top.tick_params(labelsize=12)
            if column == 0:
                top.set_ylabel("Mean posterior P(TMRCA < x)", fontsize=14)

            p_ns: dict[str, int] = {}
            for statistic in (*SECONDARY_STATISTICS, *PRIMARY_STATISTICS):
                curve = plot_p[
                    (plot_p["demography"] == demography)
                    & (plot_p["statistic"] == statistic)
                ].sort_values("threshold_years")
                if len(curve) != len(THRESHOLDS_YEARS):
                    raise ValueError(
                        f"pointwise plot curve is incomplete: {demography} {statistic}"
                    )
                p_ns[statistic] = int(curve.iloc[0]["n_neutral"])
                bottom.plot(
                    x,
                    curve["mc_p_upper"].to_numpy(float),
                    color=colors[statistic],
                    linewidth=2.6,
                    marker="o",
                    markersize=5.5,
                )
            bottom.axhline(0.05, color="#9E2A2B", linestyle=":", linewidth=2)
            bottom.set_ylim(0, 1)
            bottom.set_xticks(x)
            bottom.set_xticklabels(["1", "4.5", "10", "20", "30", "40", "50"])
            bottom.set_xlabel("TMRCA threshold x (kya)", fontsize=14)
            delta_omnibus = omnibus[
                (omnibus["demography"] == demography)
                & (omnibus["statistic"] == "p_hom_alt_minus_hom_ref")
            ]
            if len(delta_omnibus) != 1:
                raise ValueError(f"delta omnibus result is incomplete: {demography}")
            bottom.set_title(
                (
                    "Pointwise upper-tail simulation p-values\n"
                    f"neutral n: REF/REF={p_ns['p_hom_ref']}, "
                    f"ALT/ALT={p_ns['p_hom_alt']}, "
                    f"delta={p_ns['p_hom_alt_minus_hom_ref']}\n"
                    "delta 7-threshold minP: "
                    f"p={float(delta_omnibus.iloc[0]['omnibus_mc_p_upper']):.3g}"
                ),
                fontsize=11.5,
                pad=8,
            )
            bottom.grid(True, color="#d9d9d9", linewidth=0.7)
            bottom.tick_params(labelsize=12)
            if column == 0:
                bottom.set_ylabel("Upper-tail plus-one p-value", fontsize=14)
        handles: list[Any] = [
            Line2D(
                [0], [0], color=colors["p_hom_ref"], lw=3, label=labels["p_hom_ref"]
            ),
            Line2D(
                [0], [0], color=colors["p_hom_alt"], lw=3, label=labels["p_hom_alt"]
            ),
            Line2D(
                [0],
                [0],
                color=colors["p_hom_alt_minus_hom_ref"],
                lw=3,
                label=labels["p_hom_alt_minus_hom_ref"],
            ),
            Line2D([0], [0], color="#333333", lw=3, label="Selected"),
            Line2D([0], [0], color="#333333", lw=2.4, ls="--", label="Neutral median"),
            Patch(color="#999999", alpha=0.2, label="Neutral 95% range"),
        ]
        figure.legend(
            handles=handles,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.93),
            ncol=3,
            frameon=False,
            fontsize=12.5,
        )
        figure.suptitle(
            "Gamma-SMC focal pair classes and simulation-reference p-values",
            fontsize=20,
            y=0.985,
        )
        figure.subplots_adjust(
            left=0.09,
            right=0.985,
            bottom=0.085,
            top=0.79,
            wspace=0.17,
            hspace=0.52,
        )
        return _atomic_save_figure_pair(figure, Path(output_dir).resolve(), stem)
    finally:
        figure.clear()


def _validate_figure(path: Path, extension: str) -> None:
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"result figure is missing or empty: {path}")
    payload = path.read_bytes()
    if extension == "png":
        if (
            len(payload) < 45
            or payload[:8] != b"\x89PNG\r\n\x1a\n"
            or struct.unpack(">I", payload[8:12])[0] != 13
            or payload[12:16] != b"IHDR"
            or payload[-12:-8] != b"\x00\x00\x00\x00"
            or payload[-8:-4] != b"IEND"
        ):
            raise ValueError("result PNG structure is invalid")
        width, height = struct.unpack(">II", payload[16:24])
        if (width, height) != (PNG_WIDTH_PIXELS, PNG_HEIGHT_PIXELS):
            raise ValueError("result PNG dimensions are not 3300x2550")
        return
    if extension != "pdf":
        raise ValueError(f"unknown result figure extension: {extension}")
    if len(payload) < 200 or re.match(rb"%PDF-\d\.\d\s", payload) is None:
        raise ValueError("result PDF signature or size is invalid")
    end_match = re.search(rb"startxref\s+(\d+)\s+%%EOF\s*$", payload)
    if end_match is None:
        raise ValueError("result PDF trailer is invalid")
    xref_offset = int(end_match.group(1))
    if xref_offset < 0 or payload[xref_offset : xref_offset + 4] != b"xref":
        raise ValueError("result PDF cross-reference pointer is invalid")
    if (
        re.search(rb"/Type\s*/Catalog\b", payload) is None
        or re.search(rb"trailer\s*<<[\s\S]*?/Root\s+\d+\s+\d+\s+R", payload) is None
    ):
        raise ValueError("result PDF catalog is invalid")
    page_objects = re.findall(rb"/Type\s*/Page(?!s)\b", payload)
    if len(page_objects) != 1 or re.search(rb"/Count\s+1\b", payload) is None:
        raise ValueError("result PDF must contain exactly one page")
    number = rb"[-+]?(?:\d+(?:\.\d*)?|\.\d+)"
    media_boxes = re.findall(
        rb"/MediaBox\s*\[\s*("
        + number
        + rb")\s+("
        + number
        + rb")\s+("
        + number
        + rb")\s+("
        + number
        + rb")\s*\]",
        payload,
    )
    if len(media_boxes) != 1:
        raise ValueError("result PDF must declare one page MediaBox")
    media_box = tuple(float(value) for value in media_boxes[0])
    expected_media_box = (0.0, 0.0, PDF_WIDTH_POINTS, PDF_HEIGHT_POINTS)
    if not np.allclose(media_box, expected_media_box, rtol=0, atol=0.01):
        raise ValueError("result PDF is not letter-landscape")


def _run_results_text(pointwise: pd.DataFrame, omnibus: pd.DataFrame) -> str:
    displayed = pointwise[
        (
            (pointwise["analysis_role"] == SECONDARY_ROLE)
            & pointwise["statistic"].isin(SECONDARY_STATISTICS)
        )
        | (pointwise["analysis_role"] == PRIMARY_ROLE)
    ]
    lines = [
        "# CoSi2 Gamma-SMC conditioned focal pair results",
        "",
        "REF/REF and ALT/ALT denote unordered pairs from the allele-blind fixed ",
        "100-haplotype all-pairs panel; they are pseudo-diploid pair classes, ",
        "not the original biological diploid pairing. Probabilities were reduced ",
        "post hoc from retained float32 Gamma shape/rate streams with SciPy.",
        "",
        "The primary ALT/ALT-minus-REF/REF contrast uses the shared neutral bank ",
        "in which both homozygous pair classes exist. Direct REF/REF and ALT/ALT ",
        "tests use their maximal eligible banks and are secondary comparisons.",
        "",
        "## Pointwise results plotted",
        "",
        "| Role | Demography | Statistic | Best upper-tail threshold | Selected score | Upper p | Two-sided p | Neutral n |",
        "|---|---|---|---:|---:|---:|---:|---:|",
    ]
    for (role, demography, statistic), group in displayed.groupby(
        ["analysis_role", "demography", "statistic"], sort=True
    ):
        row = group.sort_values(["mc_p_upper", "threshold_years"]).iloc[0]
        lines.append(
            f"| {role} | {demography} | {statistic} | {int(row['threshold_years'])} | "
            f"{float(row['observed_selected']):.6g} | {float(row['mc_p_upper']):.6g} | "
            f"{float(row['mc_p_two_sided']):.6g} | {int(row['n_neutral'])} |"
        )
    lines.extend(
        [
            "",
            "## Whole-vector seven-threshold minP",
            "",
            "| Demography | Statistic | Upper omnibus p | Two-sided omnibus p |",
            "|---|---|---:|---:|",
        ]
    )
    for row in omnibus.to_dict(orient="records"):
        lines.append(
            f"| {row['demography']} | {row['statistic']} | "
            f"{float(row['omnibus_mc_p_upper']):.6g} | "
            f"{float(row['omnibus_mc_p_two_sided']):.6g} |"
        )
    lines.extend(
        [
            "",
            "Pointwise p-values use conservative ties and the plus-one correction. ",
            "The minP test retains each neutral unit's complete seven-threshold vector.",
            "",
        ]
    )
    return "\n".join(lines)


def _implementation_bindings(repo_root: Path) -> dict[str, dict[str, Any]]:
    paths = {
        "focal_pairs": MODULE_PATH,
        "conditioned_workflow": Path(conditioned_workflow.__file__).resolve(),
        "gamma_analysis": Path(gamma_analysis.__file__).resolve(),
        "gamma_decode": Path(gamma_decode.__file__).resolve(),
    }
    return {
        label: {"path": _relative(path, repo_root), "sha256": sha256_file(path)}
        for label, path in paths.items()
    }


def _scientific_contract() -> dict[str, Any]:
    return {
        "focal_position_0based": FOCAL_POSITION_0BASED,
        "thresholds_years": list(THRESHOLDS_YEARS),
        "mutation_rate_per_bp_per_generation": MUTATION_RATE,
        "sample_haplotypes": gamma_decode.DEFAULT_SAMPLE_HAPLOTYPES,
        "pair_panel_haplotypes": gamma_decode.DEFAULT_PAIR_PANEL_HAPLOTYPES,
        "pair_panel_seed": gamma_decode.DEFAULT_PAIR_PANEL_SEED,
        "pair_count": gamma_decode.DEFAULT_PAIR_COUNT,
        "pairing_scheme": "fixed_sha256_ranked_100_haplotypes_all_unordered_pairs",
        "pairing_unit": "unordered_haplotype_pair_pseudodiploid",
        "posterior_reduction": "retained_float32_alpha_beta_scipy_gammainc",
        "primary_bank": "shared_units_with_ref_ref_and_alt_alt_pairs",
        "pointwise_pvalue": "plus_one_conservative_ties",
        "omnibus": "exchangeable_whole_unit_minP_over_seven_thresholds",
        "reconciliation_tolerance": RECONCILIATION_TOLERANCE,
    }


def _result_contract(
    repo_root: Path,
    *,
    work_dir: Path,
    base_results_dir: Path,
    base_results_completion_sha256: str,
    base_results_contract_sha256: str,
    scores: pd.DataFrame,
    eligibility: pd.DataFrame,
    input_manifest: pd.DataFrame,
    reconciliation: pd.DataFrame,
) -> dict[str, Any]:
    reconciliation_ids = _expected_reconciliation_unit_ids(eligibility)
    return {
        "schema": SCHEMA_VERSION,
        "scientific_contract": _scientific_contract(),
        "source": {
            "work_path": _relative(work_dir, repo_root),
            "base_results_path": _relative(base_results_dir, repo_root),
            "base_results_completion_sha256": base_results_completion_sha256,
            "base_results_contract_sha256": base_results_contract_sha256,
            "input_manifest_sha256": _portable_frame_sha256(input_manifest),
            "scores_sha256": _portable_frame_sha256(scores),
            "eligibility_sha256": _portable_frame_sha256(eligibility),
            "reconciliation_sha256": _portable_frame_sha256(reconciliation),
            "reconciliation_expected_unit_count": len(reconciliation_ids),
            "reconciliation_expected_row_count": len(reconciliation_ids)
            * len(THRESHOLDS_YEARS),
            "reconciliation_expected_unit_ids_sha256": _canonical_sha256(
                list(reconciliation_ids)
            ),
        },
        "implementation": _implementation_bindings(repo_root),
        "runtime_versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "zstandard": zstandard.__version__,
            "matplotlib": matplotlib.__version__,
        },
    }


def write_results_bundle(
    repo_root: str | Path,
    *,
    results_dir: str | Path,
    work_dir: str | Path,
    base_results_dir: str | Path,
    base_results_completion_sha256: str,
    base_results_contract_sha256: str,
    scores: pd.DataFrame,
    eligibility: pd.DataFrame,
    input_manifest: pd.DataFrame,
    reconciliation: pd.DataFrame,
) -> dict[str, Any]:
    """Atomically create or strictly reuse the focal-pair sibling bundle."""

    root = Path(repo_root).resolve()
    destination = Path(results_dir).resolve()
    frame = validate_scores(scores)
    checked_reconciliation = validate_reconciliation(reconciliation, eligibility)
    pointwise = pointwise_simulation_pvalues(frame, eligibility)
    omnibus = minp_omnibus(frame, eligibility)
    contract = _result_contract(
        root,
        work_dir=Path(work_dir).resolve(),
        base_results_dir=Path(base_results_dir).resolve(),
        base_results_completion_sha256=base_results_completion_sha256,
        base_results_contract_sha256=base_results_contract_sha256,
        scores=frame,
        eligibility=eligibility,
        input_manifest=input_manifest,
        reconciliation=checked_reconciliation,
    )
    if destination.exists() and not destination.is_dir():
        raise ValueError("focal result destination exists but is not a directory")
    if destination.exists() and any(destination.iterdir()):
        return verify_results_bundle(
            root,
            results_dir=destination,
            expected_contract=contract,
            expected_base_results_dir=base_results_dir,
            expected_work_dir=work_dir,
            verify_inputs=True,
        )
    if destination.exists():
        destination.rmdir()
    destination.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(
        tempfile.mkdtemp(prefix=f".{destination.name}.staging.", dir=destination.parent)
    )
    try:
        frames = {
            "scores": frame,
            "eligibility": eligibility,
            "input_manifest": input_manifest,
            "reconciliation": checked_reconciliation,
            "pointwise": pointwise,
            "omnibus": omnibus,
        }
        for label, value in frames.items():
            _atomic_frame(stage / RESULT_OUTPUTS[label], value)
        _atomic_text(
            stage / RESULT_OUTPUTS["run_results"],
            _run_results_text(pointwise, omnibus),
        )
        figures = plot_focal_pair_results(frame, eligibility, pointwise, omnibus, stage)
        if figures != {
            "png": stage / RESULT_OUTPUTS["figure_png"],
            "pdf": stage / RESULT_OUTPUTS["figure_pdf"],
        }:
            raise ValueError("focal figure helper returned unexpected paths")
        _validate_figure(figures["png"], "png")
        _validate_figure(figures["pdf"], "pdf")
        outputs = {
            label: _output_record(
                stage / filename, rows=len(frames[label]) if label in frames else None
            )
            for label, filename in RESULT_OUTPUTS.items()
        }
        completion = {
            "schema": RESULTS_SCHEMA,
            "status": "complete",
            "contract": contract,
            "contract_sha256": _canonical_sha256(contract),
            "outputs": outputs,
            "output_count": len(outputs),
        }
        _atomic_json(stage / RESULTS_COMPLETION_FILENAME, completion)
        if {path.name for path in stage.iterdir()} != {
            RESULTS_COMPLETION_FILENAME,
            *RESULT_OUTPUTS.values(),
        }:
            raise ValueError("staged focal result inventory differs")
        os.replace(stage, destination)
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    return {
        **_read_json(destination / RESULTS_COMPLETION_FILENAME, "focal completion"),
        "cache_hit": False,
    }


def verify_results_bundle(
    repo_root: str | Path,
    *,
    results_dir: str | Path,
    expected_contract: Mapping[str, Any] | None = None,
    expected_base_results_dir: str | Path | None = None,
    expected_work_dir: str | Path | None = None,
    verify_inputs: bool = True,
) -> dict[str, Any]:
    """Verify exact inventory, checksums, inputs, and recomputable inference."""

    root = Path(repo_root).resolve()
    destination = Path(results_dir).resolve()
    expected_names = {RESULTS_COMPLETION_FILENAME, *RESULT_OUTPUTS.values()}
    children = list(destination.iterdir()) if destination.is_dir() else []
    if (
        destination.is_symlink()
        or {path.name for path in children} != expected_names
        or any(path.is_symlink() or path.is_dir() for path in children)
    ):
        raise ValueError("focal results inventory mismatch")
    completion = _read_json(
        destination / RESULTS_COMPLETION_FILENAME, "focal results completion"
    )
    contract = completion.get("contract")
    if (
        completion.get("schema") != RESULTS_SCHEMA
        or completion.get("status") != "complete"
        or not isinstance(contract, dict)
        or completion.get("contract_sha256") != _canonical_sha256(contract)
    ):
        raise ValueError("focal results completion contract is invalid")
    if expected_contract is not None and _canonical_sha256(
        contract
    ) != _canonical_sha256(expected_contract):
        raise ValueError("focal results are stale for current inputs")
    if contract.get("scientific_contract") != _scientific_contract():
        raise ValueError("focal result scientific binding is stale")
    if contract.get("implementation") != _implementation_bindings(root):
        raise ValueError("focal result implementation binding is stale")
    expected_runtime = {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "zstandard": zstandard.__version__,
        "matplotlib": matplotlib.__version__,
    }
    if contract.get("runtime_versions") != expected_runtime:
        raise ValueError("focal result runtime binding is stale")
    recorded_relative = str(contract.get("source", {}).get("base_results_path", ""))
    if not recorded_relative:
        raise ValueError("focal result base-results path is absent")
    recorded_base = (root / Path(recorded_relative)).resolve()
    if _relative(recorded_base, root) != Path(recorded_relative).as_posix():
        raise ValueError("focal result base-results path is not normalized")
    if (
        expected_base_results_dir is not None
        and recorded_base.resolve() != Path(expected_base_results_dir).resolve()
    ):
        raise ValueError(
            "focal result base-results path differs from the requested base"
        )
    recorded_work_relative = str(contract.get("source", {}).get("work_path", ""))
    if not recorded_work_relative:
        raise ValueError("focal result decode-work path is absent")
    recorded_work = (root / Path(recorded_work_relative)).resolve()
    if _relative(recorded_work, root) != Path(recorded_work_relative).as_posix():
        raise ValueError("focal result decode-work path is not normalized")
    if (
        expected_work_dir is not None
        and recorded_work != Path(expected_work_dir).resolve()
    ):
        raise ValueError(
            "focal result decode-work path differs from the requested work root"
        )
    base_completion = recorded_base / "RESULTS_COMPLETION.json"
    if sha256_file(base_completion) != contract.get("source", {}).get(
        "base_results_completion_sha256"
    ):
        raise ValueError("conditioned base result receipt changed")
    base_completion_payload = _read_json(
        base_completion, "conditioned base result receipt"
    )
    if base_completion_payload.get("contract_sha256") != contract.get("source", {}).get(
        "base_results_contract_sha256"
    ):
        raise ValueError("conditioned base result contract binding differs")
    outputs = completion.get("outputs")
    if (
        not isinstance(outputs, dict)
        or set(outputs) != set(RESULT_OUTPUTS)
        or completion.get("output_count") != len(RESULT_OUTPUTS)
    ):
        raise ValueError("focal output manifest differs")
    for label, filename in RESULT_OUTPUTS.items():
        _validate_record(
            destination / filename, outputs[label], f"focal output {label}"
        )

    scores = validate_scores(
        pd.read_csv(
            destination / RESULT_OUTPUTS["scores"],
            sep="\t",
            float_precision="round_trip",
        )
    )
    eligibility = pd.read_csv(
        destination / RESULT_OUTPUTS["eligibility"],
        sep="\t",
        keep_default_na=False,
        float_precision="round_trip",
    )
    input_manifest = pd.read_csv(
        destination / RESULT_OUTPUTS["input_manifest"],
        sep="\t",
        keep_default_na=False,
        float_precision="round_trip",
    )
    if _portable_frame_sha256(scores) != contract["source"]["scores_sha256"]:
        raise ValueError("focal score payload binding differs")
    if _portable_frame_sha256(eligibility) != contract["source"]["eligibility_sha256"]:
        raise ValueError("focal eligibility payload binding differs")
    if (
        _portable_frame_sha256(input_manifest)
        != contract["source"]["input_manifest_sha256"]
    ):
        raise ValueError("focal input manifest binding differs")
    if verify_inputs:
        for row in input_manifest.to_dict(orient="records"):
            path = root / Path(str(row["path"]))
            if (
                path.is_symlink()
                or not path.is_file()
                or path.stat().st_size != int(row["size_bytes"])
                or sha256_file(path) != str(row["sha256"])
            ):
                raise ValueError(
                    f"focal source input differs: {row['unit_id']} {row['input_label']}"
                )
    expected_frames = {
        "pointwise": pointwise_simulation_pvalues(scores, eligibility),
        "omnibus": minp_omnibus(scores, eligibility),
    }
    for label, expected in expected_frames.items():
        observed = pd.read_csv(
            destination / RESULT_OUTPUTS[label],
            sep="\t",
            float_precision="round_trip",
        )
        pd.testing.assert_frame_equal(
            observed,
            expected,
            check_dtype=False,
            check_exact=False,
            rtol=0,
            atol=1e-12,
        )
    reconciliation = pd.read_csv(
        destination / RESULT_OUTPUTS["reconciliation"],
        sep="\t",
        float_precision="round_trip",
    )
    reconciliation = validate_reconciliation(reconciliation, eligibility)
    reconciliation_ids = _expected_reconciliation_unit_ids(eligibility)
    source_contract = contract.get("source", {})
    if (
        _portable_frame_sha256(reconciliation)
        != source_contract.get("reconciliation_sha256")
        or len(reconciliation_ids)
        != source_contract.get("reconciliation_expected_unit_count")
        or len(reconciliation)
        != source_contract.get("reconciliation_expected_row_count")
        or _canonical_sha256(list(reconciliation_ids))
        != source_contract.get("reconciliation_expected_unit_ids_sha256")
    ):
        raise ValueError("focal reconciliation payload binding differs")
    _validate_figure(destination / RESULT_OUTPUTS["figure_png"], "png")
    _validate_figure(destination / RESULT_OUTPUTS["figure_pdf"], "pdf")
    if (
        not (destination / RESULT_OUTPUTS["run_results"])
        .read_text(encoding="utf-8")
        .strip()
    ):
        raise ValueError("focal RUN_RESULTS.md is empty")
    return {**completion, "cache_hit": True}


def _reduce_worker(arguments: tuple[Any, ...]) -> tuple[Any, ...]:
    return reduce_focal_unit(
        arguments[0],
        arguments[1],
        unit_dir=arguments[2],
        expected_decode_sha256=arguments[3],
    )


def analyze_study(
    repo_root: str | Path,
    *,
    work_dir: str | Path | None = None,
    base_results_dir: str | Path | None = None,
    results_dir: str | Path | None = None,
    max_workers: int = DEFAULT_MAX_WORKERS,
) -> dict[str, Any]:
    """Consume the completed study and atomically publish focal pair results."""

    root = Path(repo_root).resolve()
    work = (
        Path(work_dir).resolve()
        if work_dir is not None
        else root / DEFAULT_WORK_RELATIVE
    )
    base_results = (
        Path(base_results_dir).resolve()
        if base_results_dir is not None
        else root / DEFAULT_BASE_RESULTS_RELATIVE
    )
    destination = (
        Path(results_dir).resolve()
        if results_dir is not None
        else root / DEFAULT_RESULTS_RELATIVE
    )
    if destination == base_results or base_results in destination.parents:
        raise ValueError(
            "focal outputs must be a sibling of the immutable base results"
        )
    workers = int(max_workers)
    if workers < 1 or workers > MAX_WORKERS:
        raise ValueError(f"max_workers must lie in [1,{MAX_WORKERS}]")
    conditioned_workflow.verify_results_bundle(root, results_dir=base_results)
    base_completion_path = base_results / "RESULTS_COMPLETION.json"
    base_completion = _read_json(base_completion_path, "conditioned base completion")
    if destination.exists() and not destination.is_dir():
        raise ValueError("focal result destination exists but is not a directory")
    if destination.exists() and any(destination.iterdir()):
        return verify_results_bundle(
            root,
            results_dir=destination,
            expected_base_results_dir=base_results,
            expected_work_dir=work,
            verify_inputs=True,
        )

    units = pd.read_csv(base_results / "unit_metadata.tsv", sep="\t")
    if len(units) != conditioned_workflow.EXPECTED_TOTAL_UNITS:
        raise ValueError("conditioned unit inventory cardinality differs")
    decode_hashes = base_completion["contract"]["input_payloads"][
        "decode_completion_sha256"
    ]
    arguments = [
        (
            root,
            record,
            work / conditioned_workflow.DECODE_DIRNAME / str(record["unit_id"]),
            str(decode_hashes[str(record["unit_id"])]),
        )
        for record in units.to_dict(orient="records")
    ]
    reduced: list[tuple[Any, ...] | None] = [None] * len(arguments)
    if workers == 1:
        for index, argument in enumerate(arguments):
            reduced[index] = _reduce_worker(argument)
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(_reduce_worker, argument): index
                for index, argument in enumerate(arguments)
            }
            for future in as_completed(futures):
                reduced[futures[future]] = future.result()
    if any(value is None for value in reduced):
        raise RuntimeError("focal reducer did not return every unit")
    complete = [value for value in reduced if value is not None]
    scores = validate_scores(
        pd.concat(
            [value[0] for value in complete if not value[0].empty], ignore_index=True
        )
    )
    eligibility = (
        pd.DataFrame([value[1] for value in complete])
        .sort_values("unit_id", kind="mergesort")
        .reset_index(drop=True)
    )
    reconciliation = pd.concat(
        [value[2] for value in complete if not value[2].empty], ignore_index=True
    )
    input_manifest = (
        pd.concat([value[3] for value in complete], ignore_index=True)
        .sort_values(["unit_id", "input_label"], kind="mergesort")
        .reset_index(drop=True)
    )

    base_scores = pd.read_csv(
        base_results / "cosi2_gamma_conditioned_scores.tsv", sep="\t"
    )
    official = base_scores[
        (base_scores["statistic"] == "mean_p_tmrca_lt")
        & (base_scores["window"] == "focal")
    ][["unit_id", "threshold_years", "score"]].rename(
        columns={"score": "decoder_cpp_overall_score"}
    )
    reconciliation = reconciliation.merge(
        official,
        on=["unit_id", "threshold_years"],
        how="left",
        validate="one_to_one",
    )
    if reconciliation["decoder_cpp_overall_score"].isna().any():
        raise ValueError("overall reconciliation lacks official focal scores")
    reconciliation["signed_difference"] = (
        reconciliation["posthoc_overall_score"]
        - reconciliation["decoder_cpp_overall_score"]
    )
    reconciliation["absolute_difference"] = reconciliation["signed_difference"].abs()
    if reconciliation["absolute_difference"].max() > RECONCILIATION_TOLERANCE:
        raise ValueError("posthoc/C++ overall reconciliation exceeds tolerance")
    return write_results_bundle(
        root,
        results_dir=destination,
        work_dir=work,
        base_results_dir=base_results,
        base_results_completion_sha256=sha256_file(base_completion_path),
        base_results_contract_sha256=str(base_completion["contract_sha256"]),
        scores=scores,
        eligibility=eligibility,
        input_manifest=input_manifest,
        reconciliation=reconciliation,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=("analyze", "verify"))
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--work-dir", type=Path)
    parser.add_argument("--base-results-dir", type=Path)
    parser.add_argument("--results-dir", type=Path)
    parser.add_argument("--workers", type=int, default=DEFAULT_MAX_WORKERS)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    root = args.repo_root.resolve()
    results = (
        args.results_dir.resolve()
        if args.results_dir is not None
        else root / DEFAULT_RESULTS_RELATIVE
    )
    if args.phase == "analyze":
        completion = analyze_study(
            root,
            work_dir=args.work_dir,
            base_results_dir=args.base_results_dir,
            results_dir=results,
            max_workers=args.workers,
        )
    else:
        work = (
            args.work_dir.resolve()
            if args.work_dir is not None
            else root / DEFAULT_WORK_RELATIVE
        )
        base_results = (
            args.base_results_dir.resolve()
            if args.base_results_dir is not None
            else root / DEFAULT_BASE_RESULTS_RELATIVE
        )
        conditioned_workflow.verify_results_bundle(root, results_dir=base_results)
        completion = verify_results_bundle(
            root,
            results_dir=results,
            expected_base_results_dir=base_results,
            expected_work_dir=work,
            verify_inputs=True,
        )
    print(
        json.dumps(
            {
                "status": completion["status"],
                "cache_hit": completion["cache_hit"],
                "results_dir": str(results),
                "contract_sha256": completion["contract_sha256"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


__all__ = [
    "DEFAULT_RESULTS_RELATIVE",
    "FIGURE_STEM",
    "PRIMARY_ROLE",
    "RESULTS_COMPLETION_FILENAME",
    "RESULT_OUTPUTS",
    "SCHEMA_VERSION",
    "THRESHOLDS_YEARS",
    "analyze_study",
    "build_parser",
    "main",
    "minp_omnibus",
    "plot_focal_pair_results",
    "pointwise_simulation_pvalues",
    "reduce_focal_unit",
    "validate_reconciliation",
    "validate_scores",
    "verify_results_bundle",
    "write_results_bundle",
]
