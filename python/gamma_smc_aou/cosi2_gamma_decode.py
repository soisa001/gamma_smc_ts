"""Checksum-bound CoSi2 ms-to-VCF and Gamma-SMC decode helpers.

CoSi2 writes ms-format haplotypes, optionally preceded by raw ARG edge records.
Gamma-SMC does not consume ms format.  This module keeps that format boundary
explicit: it validates one complete CoSi2 replicate, emits a deterministic
phased VCF and an explicit allele-blind pair panel, and then runs the
repository Gamma-SMC decoder through its VCF interface.

The primary pair panel is deliberately independent of both CoSi2 output order
and focal-allele sample frequency.  A fixed SHA-256 ranking selects 100 of the
200 haplotypes, and all C(100, 2) = 4,950 unordered pairs are decoded in every
selected and natural-neutral simulation.  This matches the Gamma-SMC paper's
100-haplotype all-pairs statistic while remaining defined when a
population-surviving neutral allele is absent from the sampled haplotypes.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import shutil
import tempfile
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from itertools import combinations, pairwise
from pathlib import Path
from typing import Any, TextIO

import numpy as np
import pandas as pd

SCHEMA_VERSION = "gamma-smc.cosi2-gamma-decode/v2"
CONVERSION_SCHEMA = "gamma-smc.cosi2-vcf-conversion/v3"
DECODE_SCHEMA = "gamma-smc.cosi2-unit-decode/v3"

DEFAULT_SEQUENCE_LENGTH_BP = 10_000_000
DEFAULT_FOCAL_POSITION_0BASED = 5_000_000
DEFAULT_SAMPLE_HAPLOTYPES = 200
DEFAULT_SAMPLE_DIPLOIDS = 100
DEFAULT_PAIR_PANEL_HAPLOTYPES = 100
DEFAULT_PAIR_PANEL_SEED = 1_729
DEFAULT_PAIR_COUNT = math.comb(DEFAULT_PAIR_PANEL_HAPLOTYPES, 2)
DEFAULT_OUTPUT_STRIDE_BP = 10_000
DEFAULT_LOCAL_HALF_WIDTH_BP = 100_000
DEFAULT_MUTATION_RATE = 1.25e-8
DEFAULT_RECOMBINATION_RATE = 1.0e-8
DEFAULT_CACHE_SIZE_BP = 1_000
DEFAULT_PAIR_BLOCK = 256
TMRCA_THRESHOLDS_YEARS = (1_000, 4_500, 10_000, 20_000, 30_000, 40_000, 50_000)

MODULE_PATH = Path(__file__).resolve()
_HAPLOTYPE_RE = re.compile(r"^[01]+$")


@dataclass(frozen=True)
class ParsedCosiMs:
    """One validated CoSi2 ms replicate held in a compact in-memory form."""

    sample_haplotypes: int
    replicate_count: int
    header_seed: int
    replicate_seed: int | None
    arg_edge_count: int
    positions: tuple[float, ...]
    mutation_times_generations: tuple[float, ...] | None
    haplotypes: tuple[str, ...]

    @property
    def segregating_sites(self) -> int:
        return len(self.positions)


@dataclass(frozen=True)
class ConversionArtifacts:
    output_dir: Path
    vcf_path: Path
    pairs_path: Path
    metadata_path: Path
    completion_path: Path
    cache_hit: bool


@dataclass(frozen=True)
class DecodeArtifacts:
    unit_dir: Path
    conversion: ConversionArtifacts
    decoded_dir: Path
    summary_path: Path
    scores_path: Path
    posterior_path: Path
    pair_manifest_path: Path
    run_manifest_path: Path
    completion_path: Path
    cache_hit: bool


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_sha256(value: Any) -> str:
    payload = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _atomic_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        temporary.write_text(
            json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
            newline="\n",
        )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_frame(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        frame.to_csv(
            temporary,
            sep="\t",
            index=False,
            lineterminator="\n",
            float_format="%.12g",
        )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _positive_integer(value: Any, label: str) -> int:
    if isinstance(value, bool):
        raise TypeError(f"{label} must be a positive integer")
    try:
        numeric = int(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be a positive integer") from error
    if numeric < 1 or numeric != value:
        raise ValueError(f"{label} must be a positive integer")
    return numeric


def _nonnegative_integer(value: Any, label: str) -> int:
    if isinstance(value, bool):
        raise TypeError(f"{label} must be a nonnegative integer")
    try:
        numeric = int(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be a nonnegative integer") from error
    if numeric < 0 or numeric != value:
        raise ValueError(f"{label} must be a nonnegative integer")
    return numeric


def _positive_float(value: Any, label: str) -> float:
    try:
        numeric = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be finite and positive") from error
    if not math.isfinite(numeric) or numeric <= 0:
        raise ValueError(f"{label} must be finite and positive")
    return numeric


def _parse_int(token: str, label: str, *, minimum: int = 0) -> int:
    try:
        value = int(token)
    except ValueError as error:
        raise ValueError(f"invalid {label}: {token!r}") from error
    if value < minimum:
        raise ValueError(f"{label} must be >= {minimum}")
    return value


def _parse_floats(tokens: Sequence[str], label: str) -> tuple[float, ...]:
    try:
        values = tuple(float(token) for token in tokens)
    except ValueError as error:
        raise ValueError(f"{label} contains a nonnumeric token") from error
    if not all(math.isfinite(value) for value in values):
        raise ValueError(f"{label} contains a nonfinite value")
    return values


def _validate_arg_edge(line: str, line_number: int) -> None:
    fields = line.split()
    if len(fields) < 8 or fields[0] != "E" or fields[1] not in {"R", "G", "C"}:
        raise ValueError(f"malformed ARG edge at line {line_number}")
    if (len(fields) - 6) % 2:
        raise ValueError(
            f"ARG edge has an unpaired segment bound at line {line_number}"
        )
    _parse_int(fields[2], "ARG node_1")
    _parse_int(fields[3], "ARG node_2")
    generations = _parse_floats(fields[4:6], "ARG node generations")
    if generations[0] < 0 or generations[1] < generations[0]:
        raise ValueError(f"ARG edge generations are invalid at line {line_number}")
    bounds = _parse_floats(fields[6:], "ARG segment bounds")
    for begin, end in zip(bounds[::2], bounds[1::2], strict=True):
        # CoSi2 can print a zero-width interval when two distinct internal
        # breakpoints round to the same text value.  It is harmless here
        # because ARG edges are syntax-checked and intentionally discarded.
        if not 0 <= begin <= end <= 1:
            raise ValueError(f"ARG segment bounds are invalid at line {line_number}")


def _next_line(handle: TextIO, line_number: int) -> tuple[str, int]:
    line = handle.readline()
    if line == "":
        raise ValueError("unexpected end of CoSi2 ms file")
    return line.rstrip("\r\n"), line_number + 1


def parse_cosi_ms(
    path: str | Path,
    *,
    expected_haplotypes: int = DEFAULT_SAMPLE_HAPLOTYPES,
) -> ParsedCosiMs:
    """Strictly stream and validate exactly one CoSi2 ms replicate.

    Optional raw ``E`` ARG records are validated and discarded.  Optional
    ``muttimes:`` values are retained.  Only the site positions and 200 binary
    haplotype strings are held in memory.
    """

    source = Path(path).resolve()
    expected = _positive_integer(expected_haplotypes, "expected_haplotypes")
    if not source.is_file():
        raise ValueError(f"CoSi2 ms input is absent: {source}")

    with source.open("r", encoding="ascii", newline="") as handle:
        line_number = 0
        header, line_number = _next_line(handle, line_number)
        fields = header.split()
        if len(fields) != 3 or fields[0] != "ms":
            raise ValueError(
                "CoSi2 ms header must be exactly 'ms <samples> <replicates>'"
            )
        samples = _parse_int(fields[1], "ms sample count", minimum=1)
        replicates = _parse_int(fields[2], "ms replicate count", minimum=1)
        if samples != expected:
            raise ValueError(
                f"CoSi2 ms sample count is {samples}; expected {expected} haplotypes"
            )
        if replicates != 1:
            raise ValueError("CoSi2 ms input must contain exactly one replicate")

        random_line, line_number = _next_line(handle, line_number)
        random_fields = random_line.split()
        if len(random_fields) != 2 or random_fields[0] != "cosi_rand":
            raise ValueError("CoSi2 ms input lacks the exact cosi_rand seed header")
        header_seed = _parse_int(random_fields[1], "cosi_rand seed")

        delimiter = ""
        while not delimiter:
            delimiter, line_number = _next_line(handle, line_number)
        if not delimiter.startswith("//"):
            raise ValueError("CoSi2 ms replicate lacks the // delimiter")
        replicate_seed: int | None = None
        match = re.search(r"(?:^|\s)seed=(\d+)(?:\s|$)", delimiter[2:].strip())
        if match is not None:
            replicate_seed = _parse_int(match.group(1), "replicate seed")
            if replicate_seed != header_seed:
                raise ValueError("cosi_rand and replicate seeds disagree")

        edge_count = 0
        while True:
            line, line_number = _next_line(handle, line_number)
            if not line:
                continue
            if line.startswith("E "):
                _validate_arg_edge(line, line_number)
                edge_count += 1
                continue
            if not line.startswith("segsites:"):
                raise ValueError(
                    f"unexpected record before segsites at line {line_number}: {line[:40]!r}"
                )
            segsite_fields = line.split()
            if len(segsite_fields) != 2 or segsite_fields[0] != "segsites:":
                raise ValueError("malformed segsites record")
            segregating_sites = _parse_int(
                segsite_fields[1], "segregating-site count", minimum=1
            )
            break

        position_line, line_number = _next_line(handle, line_number)
        position_fields = position_line.split()
        if not position_fields or position_fields[0] != "positions:":
            raise ValueError(
                "CoSi2 ms input lacks positions immediately after segsites"
            )
        positions = _parse_floats(position_fields[1:], "positions")
        if len(positions) != segregating_sites:
            raise ValueError("positions count differs from segsites")
        if any(value < 0 or value > 1 for value in positions):
            raise ValueError("normalized positions must lie in [0, 1]")
        if any(right < left for left, right in pairwise(positions)):
            raise ValueError("normalized positions are not nondecreasing")

        first_data, line_number = _next_line(handle, line_number)
        mutation_times: tuple[float, ...] | None = None
        if first_data.startswith("muttimes:"):
            mutation_fields = first_data.split()
            mutation_times = _parse_floats(mutation_fields[1:], "mutation times")
            if len(mutation_times) != segregating_sites:
                raise ValueError("mutation-time count differs from segsites")
            if any(value < 0 for value in mutation_times):
                raise ValueError("mutation times must be nonnegative")
            first_data, line_number = _next_line(handle, line_number)

        haplotypes = [first_data]
        for _ in range(1, samples):
            haplotype, line_number = _next_line(handle, line_number)
            haplotypes.append(haplotype)
        for index, haplotype in enumerate(haplotypes):
            if (
                len(haplotype) != segregating_sites
                or _HAPLOTYPE_RE.fullmatch(haplotype) is None
            ):
                raise ValueError(
                    f"haplotype {index} is not a {segregating_sites}-site binary string"
                )

        for trailing in handle:
            line_number += 1
            if trailing.strip():
                raise ValueError(
                    f"unexpected trailing record at line {line_number}; input may contain "
                    "more than one replicate"
                )

    return ParsedCosiMs(
        sample_haplotypes=samples,
        replicate_count=replicates,
        header_seed=header_seed,
        replicate_seed=replicate_seed,
        arg_edge_count=edge_count,
        positions=positions,
        mutation_times_generations=mutation_times,
        haplotypes=tuple(haplotypes),
    )


def _discrete_positions(
    normalized: Sequence[float], sequence_length_bp: int
) -> np.ndarray:
    length = _positive_integer(sequence_length_bp, "sequence_length_bp")
    normalized_values = np.asarray(normalized, dtype=float)
    values = np.floor(normalized_values * length).astype(np.int64)
    # CoSi2 can emit the closed right endpoint exactly as normalized position
    # 1.0. VCF coordinates are closed 1..L, so map that endpoint to 0-based
    # L-1 rather than rejecting a valid terminal-locus mutation.
    values[normalized_values == 1.0] = length - 1
    if np.any(values < 0) or np.any(values >= length):
        raise ValueError("discretized CoSi2 positions fall outside the contig")
    if np.any(np.diff(values) < 0):
        raise ValueError("discretized CoSi2 positions are not nondecreasing")
    return values


def _write_vcf(
    path: Path,
    parsed: ParsedCosiMs,
    *,
    sequence_length_bp: int,
    contig: str,
) -> np.ndarray:
    discrete = _discrete_positions(parsed.positions, sequence_length_bp)
    sample_names = [f"D{index:03d}" for index in range(parsed.sample_haplotypes // 2)]
    with path.open("w", encoding="ascii", newline="\n") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("##source=gamma_smc_aou.cosi2_gamma_decode\n")
        handle.write(f"##contig=<ID={contig},length={sequence_length_bp}>\n")
        handle.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        handle.write(
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased genotype">\n'
        )
        handle.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
            + "\t".join(sample_names)
            + "\n"
        )
        haplotypes = parsed.haplotypes
        for site_index, position_0based in enumerate(discrete):
            genotypes = [
                f"{haplotypes[2 * sample][site_index]}|"
                f"{haplotypes[2 * sample + 1][site_index]}"
                for sample in range(len(sample_names))
            ]
            handle.write(
                f"{contig}\t{int(position_0based) + 1}\t"
                f"cosi_s{site_index:08d}\tA\tG\t.\tPASS\t.\tGT\t"
                + "\t".join(genotypes)
                + "\n"
            )
    return discrete


def _pair_panel_indices(
    sample_haplotypes: int = DEFAULT_SAMPLE_HAPLOTYPES,
    panel_haplotypes: int = DEFAULT_PAIR_PANEL_HAPLOTYPES,
    seed: int = DEFAULT_PAIR_PANEL_SEED,
) -> tuple[int, ...]:
    """Select a cross-platform deterministic allele-blind haplotype panel."""

    total = _positive_integer(sample_haplotypes, "sample_haplotypes")
    panel = _positive_integer(panel_haplotypes, "panel_haplotypes")
    panel_seed = _nonnegative_integer(seed, "seed")
    if panel > total:
        raise ValueError("panel_haplotypes cannot exceed sample_haplotypes")
    ranked = sorted(
        range(total),
        key=lambda index: (
            hashlib.sha256(f"{panel_seed}:{index}".encode("ascii")).digest(),
            index,
        ),
    )
    selected = tuple(sorted(ranked[:panel]))
    if len(selected) != panel or len(set(selected)) != panel:
        raise AssertionError("pair panel selection is not unique")
    return selected


def _pair_panel_sha256(indices: Sequence[int]) -> str:
    return _canonical_sha256([int(index) for index in indices])


def _write_pair_table(path: Path, indices: Sequence[int], *, seed: int) -> None:
    selected = tuple(int(index) for index in indices)
    if (
        len(selected) != DEFAULT_PAIR_PANEL_HAPLOTYPES
        or len(set(selected)) != len(selected)
        or min(selected) < 0
        or max(selected) >= DEFAULT_SAMPLE_HAPLOTYPES
    ):
        raise ValueError("pair-panel indices violate the fixed 100-of-200 contract")
    pairs = tuple(combinations(selected, 2))
    if len(pairs) != DEFAULT_PAIR_COUNT:
        raise AssertionError("pair-panel cardinality is not C(100, 2)")
    with path.open("w", encoding="ascii", newline="\n") as handle:
        handle.write("# Gamma-SMC explicit haplotype pairs\n")
        handle.write(
            "# pairing_scheme\tfixed_sha256_ranked_100_haplotypes_all_unordered_pairs\n"
        )
        handle.write(f"# panel_seed\t{seed}\n")
        handle.write(f"# panel_indices_sha256\t{_pair_panel_sha256(selected)}\n")
        handle.write(f"# n_panel_haplotypes\t{len(selected)}\n")
        handle.write(f"# n_pairs\t{len(pairs)}\n")
        for left, right in pairs:
            handle.write(f"{left}\t{right}\n")


def _data_pair_rows(path: Path) -> tuple[tuple[int, int], ...]:
    rows: list[tuple[int, int]] = []
    try:
        with path.open("r", encoding="ascii") as handle:
            for line_number, raw in enumerate(handle, start=1):
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                fields = line.split("\t")
                if len(fields) not in {2, 4}:
                    raise ValueError(
                        f"pair manifest row {line_number} must have two or four columns"
                    )
                left = _parse_int(fields[0], "left haplotype")
                right = _parse_int(fields[1], "right haplotype")
                if left >= right or right >= DEFAULT_SAMPLE_HAPLOTYPES:
                    raise ValueError(
                        f"pair manifest row {line_number} has invalid indices"
                    )
                if len(fields) == 4:
                    expected_labels = (
                        f"D{left // 2:03d}.{left % 2}",
                        f"D{right // 2:03d}.{right % 2}",
                    )
                    if tuple(fields[2:]) != expected_labels:
                        raise ValueError(
                            f"pair manifest row {line_number} has invalid labels"
                        )
                rows.append((left, right))
    except UnicodeDecodeError as error:
        raise ValueError(f"pair manifest is not ASCII: {path}") from error
    return tuple(rows)


def _validate_decoded_pair_manifest(decoded: Path, requested: Path) -> None:
    expected = _data_pair_rows(requested)
    observed = _data_pair_rows(decoded)
    if len(expected) != DEFAULT_PAIR_COUNT or observed != expected:
        raise ValueError(
            "decoded pair manifest does not equal the fixed all-pairs panel"
        )


def _output_record(
    path: Path, root: Path, *, rows: int | None = None
) -> dict[str, Any]:
    record: dict[str, Any] = {
        "path": path.relative_to(root).as_posix(),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }
    if rows is not None:
        record["rows"] = int(rows)
    return record


def _validate_completion(
    root: Path,
    completion_path: Path,
    *,
    schema: str,
    contract: Mapping[str, Any],
    expected_outputs: set[str],
) -> dict[str, Any]:
    try:
        completion = json.loads(completion_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(
            f"completion record is unreadable: {completion_path}"
        ) from error
    recorded_contract = completion.get("contract")
    if (
        completion.get("schema") != schema
        or completion.get("status") != "complete"
        or not isinstance(recorded_contract, Mapping)
        or completion.get("contract_sha256") != _canonical_sha256(recorded_contract)
        or _canonical_sha256(recorded_contract) != _canonical_sha256(contract)
    ):
        raise ValueError(f"completion contract is incompatible: {completion_path}")
    outputs = completion.get("outputs")
    if not isinstance(outputs, Mapping) or set(outputs) != expected_outputs:
        raise ValueError(
            f"completion output manifest is incompatible: {completion_path}"
        )
    for label, record in outputs.items():
        if not isinstance(record, Mapping):
            raise TypeError(f"completion output record is malformed: {label}")
        output = root / str(record.get("path", ""))
        if (
            not output.is_file()
            or output.stat().st_size != int(record.get("size_bytes", -1))
            or sha256_file(output) != record.get("sha256")
        ):
            raise ValueError(f"completion output checksum failed: {label}")
    return completion


def prepare_cosi_gamma_input(
    ms_path: str | Path,
    output_dir: str | Path,
    *,
    sequence_length_bp: int = DEFAULT_SEQUENCE_LENGTH_BP,
    focal_position_0based: int = DEFAULT_FOCAL_POSITION_0BASED,
    expected_haplotypes: int = DEFAULT_SAMPLE_HAPLOTYPES,
    pair_panel_haplotypes: int = DEFAULT_PAIR_PANEL_HAPLOTYPES,
    pair_panel_seed: int = DEFAULT_PAIR_PANEL_SEED,
    contig: str = "1",
    expected_ms_sha256: str | None = None,
) -> ConversionArtifacts:
    """Convert one checksum-bound CoSi2 ms replicate to a deterministic VCF."""

    source = Path(ms_path).resolve()
    destination = Path(output_dir).resolve()
    length = _positive_integer(sequence_length_bp, "sequence_length_bp")
    focal = _positive_integer(focal_position_0based, "focal_position_0based")
    samples = _positive_integer(expected_haplotypes, "expected_haplotypes")
    panel_size = _positive_integer(pair_panel_haplotypes, "pair_panel_haplotypes")
    panel_seed = _nonnegative_integer(pair_panel_seed, "pair_panel_seed")
    if samples % 2:
        raise ValueError("expected_haplotypes must be even")
    if samples != DEFAULT_SAMPLE_HAPLOTYPES:
        raise ValueError("the paper-compatible panel requires exactly 200 haplotypes")
    if panel_size != DEFAULT_PAIR_PANEL_HAPLOTYPES:
        raise ValueError("the paper-compatible panel requires exactly 100 haplotypes")
    panel_indices = _pair_panel_indices(samples, panel_size, panel_seed)
    panel_sha256 = _pair_panel_sha256(panel_indices)
    if focal >= length:
        raise ValueError("focal_position_0based lies outside the contig")
    if not contig or any(character.isspace() for character in contig):
        raise ValueError("contig must be a nonempty token")
    if not source.is_file():
        raise ValueError(f"CoSi2 ms input is absent: {source}")
    observed_ms_sha256 = sha256_file(source)
    if expected_ms_sha256 is not None and observed_ms_sha256 != expected_ms_sha256:
        raise ValueError("CoSi2 ms input checksum differs from expected_ms_sha256")
    contract = {
        "schema": CONVERSION_SCHEMA,
        "input": {
            "path": str(source),
            "size_bytes": source.stat().st_size,
            "sha256": observed_ms_sha256,
        },
        "parameters": {
            "sequence_length_bp": length,
            "focal_position_0based": focal,
            "expected_haplotypes": samples,
            "sample_diploids": samples // 2,
            "pair_panel_haplotypes": panel_size,
            "pair_panel_seed": panel_seed,
            "pair_panel_indices": list(panel_indices),
            "pair_panel_indices_sha256": panel_sha256,
            "pair_count": DEFAULT_PAIR_COUNT,
            "contig": contig,
            "coordinate_mapping": (
                "vcf_pos=min(floor(normalized_position*length),length-1)+1"
            ),
            "duplicate_position_policy": "retain_all_records_with_unique_ids",
            "pairing_scheme": (
                "fixed_sha256_ranked_100_haplotypes_all_unordered_pairs"
            ),
        },
        "implementation": {
            "module_path": str(MODULE_PATH),
            "module_sha256": sha256_file(MODULE_PATH),
        },
    }
    completion_path = destination / "conversion_complete.json"
    expected_outputs = {"vcf", "pairs", "metadata"}
    if completion_path.is_file():
        _validate_completion(
            destination,
            completion_path,
            schema=CONVERSION_SCHEMA,
            contract=contract,
            expected_outputs=expected_outputs,
        )
        return ConversionArtifacts(
            output_dir=destination,
            vcf_path=destination / "input.vcf",
            pairs_path=destination / "overall.pairs.tsv",
            metadata_path=destination / "conversion_metadata.json",
            completion_path=completion_path,
            cache_hit=True,
        )
    if destination.exists() and any(destination.iterdir()):
        raise ValueError("conversion directory is nonempty without valid completion")
    destination.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(
        tempfile.mkdtemp(prefix=f".{destination.name}.stage.", dir=destination.parent)
    )
    try:
        parsed = parse_cosi_ms(source, expected_haplotypes=samples)
        vcf = stage / "input.vcf"
        pairs = stage / "overall.pairs.tsv"
        metadata_path = stage / "conversion_metadata.json"
        discrete = _write_vcf(
            vcf,
            parsed,
            sequence_length_bp=length,
            contig=contig,
        )
        _write_pair_table(pairs, panel_indices, seed=panel_seed)
        exact_focal_indices = [
            index
            for index, value in enumerate(parsed.positions)
            if value * length == focal
        ]
        focal_alt_count = (
            sum(
                int(haplotype[exact_focal_indices[0]])
                for haplotype in parsed.haplotypes
            )
            if len(exact_focal_indices) == 1
            else None
        )
        metadata = {
            "schema": CONVERSION_SCHEMA,
            "ms_header_seed": parsed.header_seed,
            "ms_replicate_seed": parsed.replicate_seed,
            "sample_haplotypes": parsed.sample_haplotypes,
            "sample_diploids": parsed.sample_haplotypes // 2,
            "pair_panel_haplotypes": panel_size,
            "pair_panel_seed": panel_seed,
            "pair_panel_indices": list(panel_indices),
            "pair_panel_indices_sha256": panel_sha256,
            "pair_count": DEFAULT_PAIR_COUNT,
            "segregating_sites": parsed.segregating_sites,
            "arg_edge_count_skipped": parsed.arg_edge_count,
            "mutation_times_present": parsed.mutation_times_generations is not None,
            "discrete_unique_positions": len(np.unique(discrete)),
            "discrete_position_collisions": int(
                len(discrete) - len(np.unique(discrete))
            ),
            "normalized_terminal_position_count": int(
                np.sum(np.asarray(parsed.positions, dtype=float) == 1.0)
            ),
            "exact_focal_site_present": bool(exact_focal_indices),
            "exact_focal_site_count": len(exact_focal_indices),
            "exact_focal_site_index": (
                int(exact_focal_indices[0]) if exact_focal_indices else None
            ),
            "exact_focal_sample_alt_count": focal_alt_count,
            "exact_focal_pair_panel_alt_count": (
                sum(
                    int(parsed.haplotypes[index][exact_focal_indices[0]])
                    for index in panel_indices
                )
                if len(exact_focal_indices) == 1
                else None
            ),
            "contract_sha256": _canonical_sha256(contract),
        }
        _atomic_json(metadata_path, metadata)
        outputs = {
            "vcf": _output_record(vcf, stage, rows=parsed.segregating_sites),
            "pairs": _output_record(pairs, stage, rows=DEFAULT_PAIR_COUNT),
            "metadata": _output_record(metadata_path, stage),
        }
        _atomic_json(
            stage / "conversion_complete.json",
            {
                "schema": CONVERSION_SCHEMA,
                "status": "complete",
                "contract": contract,
                "contract_sha256": _canonical_sha256(contract),
                "outputs": outputs,
            },
        )
        if destination.exists():
            destination.rmdir()
        os.replace(stage, destination)
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    return ConversionArtifacts(
        output_dir=destination,
        vcf_path=destination / "input.vcf",
        pairs_path=destination / "overall.pairs.tsv",
        metadata_path=destination / "conversion_metadata.json",
        completion_path=destination / "conversion_complete.json",
        cache_hit=False,
    )


def _thresholds(values: Sequence[int | float]) -> tuple[int, ...]:
    result: list[int] = []
    for value in values:
        numeric = _positive_float(value, "TMRCA threshold")
        integer = int(numeric)
        if integer != numeric:
            raise ValueError("TMRCA thresholds must be integer years")
        result.append(integer)
    if not result or len(set(result)) != len(result) or result != sorted(result):
        raise ValueError("TMRCA thresholds must be nonempty, unique, and increasing")
    return tuple(result)


def reduce_gamma_summary(
    summary: str | Path | pd.DataFrame,
    *,
    thresholds_years: Sequence[int | float] = TMRCA_THRESHOLDS_YEARS,
    sequence_length_bp: int = DEFAULT_SEQUENCE_LENGTH_BP,
    focal_position_0based: int = DEFAULT_FOCAL_POSITION_0BASED,
    output_stride_bp: int = DEFAULT_OUTPUT_STRIDE_BP,
    local_half_width_bp: int = DEFAULT_LOCAL_HALF_WIDTH_BP,
    expected_pairs: int = DEFAULT_PAIR_COUNT,
) -> pd.DataFrame:
    """Reduce a strict Gamma summary grid to focal and +/-100 kb scores.

    ``hard_mean_call_fraction`` is the paper-style fraction of pairwise calls
    with posterior *mean* TMRCA below the threshold (the decoder must be run
    with ``recent_call='mean'``).  ``soft_mean_cdf_score`` averages the exact
    posterior CDF values over the same pairs and output positions.
    """

    thresholds = _thresholds(thresholds_years)
    length = _positive_integer(sequence_length_bp, "sequence_length_bp")
    focal = _positive_integer(focal_position_0based, "focal_position_0based")
    stride = _positive_integer(output_stride_bp, "output_stride_bp")
    half_width = _positive_integer(local_half_width_bp, "local_half_width_bp")
    n_pairs = _positive_integer(expected_pairs, "expected_pairs")
    if focal >= length or focal % stride:
        raise ValueError("focal position must lie on the requested output stride")
    if half_width % stride:
        raise ValueError("local half-width must be an integer number of strides")
    frame = (
        summary.copy()
        if isinstance(summary, pd.DataFrame)
        else pd.read_csv(Path(summary), sep="\t")
    )
    required = {"position_0based", "position_1based", "n_pairs"}
    for threshold in thresholds:
        required.update(
            {
                f"n_recent_{threshold}",
                f"frac_recent_{threshold}",
                f"mean_p_lt_{threshold}",
            }
        )
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError("Gamma summary lacks columns: " + ", ".join(missing))
    positions = pd.to_numeric(frame["position_0based"], errors="coerce").to_numpy(
        dtype=float
    )
    positions_1 = pd.to_numeric(frame["position_1based"], errors="coerce").to_numpy(
        dtype=float
    )
    expected_positions = np.arange(0, length, stride, dtype=np.int64)
    if (
        len(frame) != len(expected_positions)
        or not np.all(np.isfinite(positions))
        or not np.array_equal(positions, expected_positions)
        or not np.array_equal(positions_1, expected_positions + 1)
    ):
        raise ValueError("Gamma summary does not contain the exact ordered stride grid")
    observed_pairs = pd.to_numeric(frame["n_pairs"], errors="coerce").to_numpy(
        dtype=float
    )
    if not np.all(observed_pairs == n_pairs):
        raise ValueError(f"Gamma summary must use exactly {n_pairs} pairs")

    hard_columns: list[np.ndarray] = []
    soft_columns: list[np.ndarray] = []
    for threshold in thresholds:
        counts = pd.to_numeric(
            frame[f"n_recent_{threshold}"], errors="coerce"
        ).to_numpy(dtype=float)
        hard = pd.to_numeric(
            frame[f"frac_recent_{threshold}"], errors="coerce"
        ).to_numpy(dtype=float)
        soft = pd.to_numeric(frame[f"mean_p_lt_{threshold}"], errors="coerce").to_numpy(
            dtype=float
        )
        if (
            not np.all(np.isfinite(counts))
            or not np.all(counts == np.floor(counts))
            or np.any((counts < 0) | (counts > n_pairs))
            or not np.all(np.isfinite(hard))
            or not np.all(np.isfinite(soft))
            or np.any((hard < 0) | (hard > 1))
            or np.any((soft < 0) | (soft > 1))
            # Gamma-SMC emits summary floats with ten significant digits.
            # This tolerance is below one recent-pair increment (1/4,950)
            # while accepting the largest possible decimal-rounding error.
            or not np.allclose(hard, counts / n_pairs, rtol=0.0, atol=1e-10)
        ):
            raise ValueError(f"Gamma summary threshold {threshold} values are invalid")
        hard_columns.append(hard)
        soft_columns.append(soft)
    hard_matrix = np.column_stack(hard_columns)
    soft_matrix = np.column_stack(soft_columns)
    if np.any(np.diff(hard_matrix, axis=1) < -1e-12) or np.any(
        np.diff(soft_matrix, axis=1) < -1e-12
    ):
        raise ValueError(
            "Gamma recent probabilities are not monotone across thresholds"
        )
    if "mean_p_tmrca_lt_threshold" in frame.columns:
        legacy = pd.to_numeric(
            frame["mean_p_tmrca_lt_threshold"], errors="coerce"
        ).to_numpy(dtype=float)
        if not np.allclose(legacy, soft_matrix[:, 0], rtol=0.0, atol=1e-12):
            raise ValueError("Gamma first-threshold compatibility column disagrees")

    focal_mask = expected_positions == focal
    local_mask = (expected_positions >= focal - half_width) & (
        expected_positions <= focal + half_width
    )
    expected_local_positions = 2 * (half_width // stride) + 1
    if focal_mask.sum() != 1 or local_mask.sum() != expected_local_positions:
        raise ValueError("focal/local output-position geometry is incomplete")
    rows: list[dict[str, Any]] = []
    for window, mask in (("focal", focal_mask), ("local_100kb", local_mask)):
        for threshold_index, threshold in enumerate(thresholds):
            rows.append(
                {
                    "window": window,
                    "threshold_years": threshold,
                    "hard_mean_call_fraction": float(
                        np.mean(hard_matrix[mask, threshold_index])
                    ),
                    "soft_mean_cdf_score": float(
                        np.mean(soft_matrix[mask, threshold_index])
                    ),
                    "n_pairs": n_pairs,
                    "n_output_positions": int(mask.sum()),
                    "focal_position_0based": focal,
                    "local_half_width_bp": half_width,
                }
            )
    return pd.DataFrame(rows)


def _decode_contract(
    *,
    unit_id: str,
    demography_id: str,
    conversion: ConversionArtifacts,
    decoder_bin: Path,
    present_ne: float,
    generation_time_years: float,
    thresholds_years: tuple[int, ...],
    mutation_rate: float,
    recombination_rate: float,
    sequence_length_bp: int,
    focal_position_0based: int,
    output_stride_bp: int,
    local_half_width_bp: int,
    threads: int,
) -> dict[str, Any]:
    return {
        "schema": DECODE_SCHEMA,
        "unit_id": unit_id,
        "demography_id": demography_id,
        "conversion": {
            "completion_path": str(conversion.completion_path),
            "completion_sha256": sha256_file(conversion.completion_path),
            "vcf_sha256": sha256_file(conversion.vcf_path),
            "pairs_sha256": sha256_file(conversion.pairs_path),
        },
        "decoder_binary": {
            "path": str(decoder_bin),
            "size_bytes": decoder_bin.stat().st_size,
            "sha256": sha256_file(decoder_bin),
        },
        "settings": {
            "input_format": "vcf",
            "pair_selector": ("fixed_sha256_ranked_100_haplotypes_all_unordered_pairs"),
            "pair_panel_haplotypes": DEFAULT_PAIR_PANEL_HAPLOTYPES,
            "pair_panel_seed": DEFAULT_PAIR_PANEL_SEED,
            "pair_count": DEFAULT_PAIR_COUNT,
            "n_pairs": DEFAULT_PAIR_COUNT,
            "present_ne": present_ne,
            "scaled_mutation_rate": 4 * present_ne * mutation_rate,
            "unscaled_mutation_rate": mutation_rate,
            "recombination_rate": recombination_rate,
            "recombination_to_mutation_ratio": recombination_rate / mutation_rate,
            "generation_time_years": generation_time_years,
            "thresholds_years": list(thresholds_years),
            "output_at_stride": output_stride_bp,
            "output_at_hets": False,
            "recent_call": "mean",
            "recent_call_probability": 0.5,
            "cache_size": DEFAULT_CACHE_SIZE_BP,
            "pair_block": DEFAULT_PAIR_BLOCK,
            "exp10": "accurate",
            "backward_alignment": "fixed",
            "threads": threads,
            "mask": None,
        },
        "score_geometry": {
            "sequence_length_bp": sequence_length_bp,
            "focal_position_0based": focal_position_0based,
            "local_half_width_bp": local_half_width_bp,
        },
        "implementation": {
            "module_path": str(MODULE_PATH),
            "module_sha256": sha256_file(MODULE_PATH),
        },
    }


def decode_cosi_unit(
    ms_path: str | Path,
    unit_dir: str | Path,
    decoder_bin: str | Path,
    *,
    unit_id: str,
    demography_id: str,
    present_ne: float,
    generation_time_years: float,
    threads: int = 4,
    thresholds_years: Sequence[int | float] = TMRCA_THRESHOLDS_YEARS,
    mutation_rate: float = DEFAULT_MUTATION_RATE,
    recombination_rate: float = DEFAULT_RECOMBINATION_RATE,
    sequence_length_bp: int = DEFAULT_SEQUENCE_LENGTH_BP,
    focal_position_0based: int = DEFAULT_FOCAL_POSITION_0BASED,
    output_stride_bp: int = DEFAULT_OUTPUT_STRIDE_BP,
    local_half_width_bp: int = DEFAULT_LOCAL_HALF_WIDTH_BP,
    expected_ms_sha256: str | None = None,
    runner: Callable[..., Mapping[str, Any]] | None = None,
) -> DecodeArtifacts:
    """Prepare and decode one immutable CoSi2 unit, or validate its cache."""

    root = Path(unit_dir).resolve()
    if not unit_id.strip() or not demography_id.strip():
        raise ValueError("unit_id and demography_id must be nonempty")
    binary = Path(decoder_bin).resolve()
    if not binary.is_file():
        raise ValueError(f"Gamma-SMC decoder is absent: {binary}")
    ne = _positive_float(present_ne, "present_ne")
    generation_time = _positive_float(generation_time_years, "generation_time_years")
    thread_count = _positive_integer(threads, "threads")
    mutation = _positive_float(mutation_rate, "mutation_rate")
    recombination = _positive_float(recombination_rate, "recombination_rate")
    length = _positive_integer(sequence_length_bp, "sequence_length_bp")
    focal = _positive_integer(focal_position_0based, "focal_position_0based")
    stride = _positive_integer(output_stride_bp, "output_stride_bp")
    half_width = _positive_integer(local_half_width_bp, "local_half_width_bp")
    thresholds = _thresholds(thresholds_years)
    root.mkdir(parents=True, exist_ok=True)
    conversion = prepare_cosi_gamma_input(
        ms_path,
        root / "converted",
        sequence_length_bp=length,
        focal_position_0based=focal,
        expected_ms_sha256=expected_ms_sha256,
    )
    contract = _decode_contract(
        unit_id=unit_id,
        demography_id=demography_id,
        conversion=conversion,
        decoder_bin=binary,
        present_ne=ne,
        generation_time_years=generation_time,
        thresholds_years=thresholds,
        mutation_rate=mutation,
        recombination_rate=recombination,
        sequence_length_bp=length,
        focal_position_0based=focal,
        output_stride_bp=stride,
        local_half_width_bp=half_width,
        threads=thread_count,
    )
    decoded = root / "decoded"
    completion_path = decoded / "decode_complete.json"
    expected_outputs = {
        "summary",
        "scores",
        "posterior",
        "posterior_metadata",
        "pair_manifest",
        "run_manifest",
    }
    if completion_path.is_file():
        _validate_completion(
            decoded,
            completion_path,
            schema=DECODE_SCHEMA,
            contract=contract,
            expected_outputs=expected_outputs,
        )
        _validate_decoded_pair_manifest(
            decoded / "decoded_pairs.tsv", conversion.pairs_path
        )
        reduce_gamma_summary(
            decoded / "overall.summary.tsv",
            thresholds_years=thresholds,
            sequence_length_bp=length,
            focal_position_0based=focal,
            output_stride_bp=stride,
            local_half_width_bp=half_width,
        )
        return DecodeArtifacts(
            unit_dir=root,
            conversion=conversion,
            decoded_dir=decoded,
            summary_path=decoded / "overall.summary.tsv",
            scores_path=decoded / "unit_scores.tsv",
            posterior_path=decoded / "posterior.zst",
            pair_manifest_path=decoded / "decoded_pairs.tsv",
            run_manifest_path=decoded / "decoder_run.json",
            completion_path=completion_path,
            cache_hit=True,
        )
    if decoded.exists() and any(decoded.iterdir()):
        raise ValueError("decode directory is nonempty without valid completion")

    stage = Path(tempfile.mkdtemp(prefix=".decoded.stage.", dir=root))
    try:
        summary_path = stage / "overall.summary.tsv"
        posterior_path = stage / "posterior.zst"
        pair_manifest_path = stage / "decoded_pairs.tsv"
        decode_runner = runner
        if decode_runner is None:
            from .decoder import run_within_decoder

            decode_runner = run_within_decoder
        run = dict(
            decode_runner(
                binary,
                conversion.vcf_path,
                summary_path,
                scaled_mutation_rate=4 * ne * mutation,
                recombination_to_mutation_ratio=recombination / mutation,
                mutation_rate=mutation,
                threshold_years=thresholds,
                generation_time=generation_time,
                input_format="vcf",
                raw_output=posterior_path,
                output_at_stride=stride,
                output_at_hets=False,
                only_within=False,
                pairs_file=conversion.pairs_path,
                pairs_manifest=pair_manifest_path,
                recent_call="mean",
                recent_call_probability=0.5,
                threads=thread_count,
                cache_size=DEFAULT_CACHE_SIZE_BP,
                pair_block=DEFAULT_PAIR_BLOCK,
                exp10="accurate",
                backward_alignment="fixed",
            )
        )
        required = {
            "summary": summary_path,
            "posterior": posterior_path,
            "posterior_metadata": Path(str(posterior_path) + ".meta"),
            "pair_manifest": pair_manifest_path,
        }
        missing = [label for label, path in required.items() if not path.is_file()]
        if missing:
            raise ValueError("Gamma-SMC decode lacks outputs: " + ", ".join(missing))
        _validate_decoded_pair_manifest(pair_manifest_path, conversion.pairs_path)
        scores = reduce_gamma_summary(
            summary_path,
            thresholds_years=thresholds,
            sequence_length_bp=length,
            focal_position_0based=focal,
            output_stride_bp=stride,
            local_half_width_bp=half_width,
        )
        scores.insert(0, "demography_id", demography_id)
        scores.insert(0, "unit_id", unit_id)
        scores.insert(2, "generation_time_years", generation_time)
        scores_path = stage / "unit_scores.tsv"
        _atomic_frame(scores_path, scores)
        automatic_run_path = summary_path.with_suffix(summary_path.suffix + ".run.json")
        run_manifest_path = stage / "decoder_run.json"
        if automatic_run_path.is_file():
            automatic_run_path.unlink()
        _atomic_json(
            run_manifest_path,
            {
                **run,
                "contract_sha256": _canonical_sha256(contract),
                "execution_paths_use_ephemeral_atomic_stage": True,
                "atomic_stage_promoted_to_decoded_directory": True,
            },
        )
        outputs = {
            "summary": _output_record(
                summary_path,
                stage,
                rows=len(np.arange(0, length, stride)),
            ),
            "scores": _output_record(scores_path, stage, rows=len(scores)),
            "posterior": _output_record(posterior_path, stage),
            "posterior_metadata": _output_record(
                Path(str(posterior_path) + ".meta"), stage
            ),
            "pair_manifest": _output_record(
                pair_manifest_path, stage, rows=DEFAULT_PAIR_COUNT
            ),
            "run_manifest": _output_record(run_manifest_path, stage),
        }
        _atomic_json(
            stage / "decode_complete.json",
            {
                "schema": DECODE_SCHEMA,
                "status": "complete",
                "contract": contract,
                "contract_sha256": _canonical_sha256(contract),
                "outputs": outputs,
            },
        )
        if decoded.exists():
            decoded.rmdir()
        os.replace(stage, decoded)
        (root / "decode_failed.json").unlink(missing_ok=True)
    except Exception as error:
        shutil.rmtree(stage, ignore_errors=True)
        _atomic_json(
            root / "decode_failed.json",
            {
                "schema": DECODE_SCHEMA,
                "status": "failed",
                "unit_id": unit_id,
                "error_type": type(error).__name__,
                "error": str(error),
                "contract_sha256": _canonical_sha256(contract),
            },
        )
        raise
    return DecodeArtifacts(
        unit_dir=root,
        conversion=conversion,
        decoded_dir=decoded,
        summary_path=decoded / "overall.summary.tsv",
        scores_path=decoded / "unit_scores.tsv",
        posterior_path=decoded / "posterior.zst",
        pair_manifest_path=decoded / "decoded_pairs.tsv",
        run_manifest_path=decoded / "decoder_run.json",
        completion_path=decoded / "decode_complete.json",
        cache_hit=False,
    )


__all__ = [
    "CONVERSION_SCHEMA",
    "DECODE_SCHEMA",
    "DEFAULT_FOCAL_POSITION_0BASED",
    "DEFAULT_LOCAL_HALF_WIDTH_BP",
    "DEFAULT_OUTPUT_STRIDE_BP",
    "DEFAULT_SAMPLE_DIPLOIDS",
    "DEFAULT_SAMPLE_HAPLOTYPES",
    "DEFAULT_SEQUENCE_LENGTH_BP",
    "SCHEMA_VERSION",
    "TMRCA_THRESHOLDS_YEARS",
    "ConversionArtifacts",
    "DecodeArtifacts",
    "ParsedCosiMs",
    "decode_cosi_unit",
    "parse_cosi_ms",
    "prepare_cosi_gamma_input",
    "reduce_gamma_summary",
    "sha256_file",
]
