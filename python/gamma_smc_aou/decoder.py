from __future__ import annotations

import json
import subprocess
import sys
from collections.abc import Sequence
from pathlib import Path
from time import perf_counter

from .defaults import (
    DEFAULT_CACHE_SIZE,
    DEFAULT_GENERATION_TIME,
    DEFAULT_MUTATION_RATE,
    DEFAULT_OUTPUT_STRIDE,
    DEFAULT_RECOMBINATION_TO_MUTATION_RATIO,
    DEFAULT_SCALED_MUTATION_RATE,
)


TREE_SEQUENCE_SUFFIXES = {".trees", ".ts", ".tsz"}


def _streams_tree_sequence(input_path: str | Path, input_format: str) -> bool:
    return input_format in {"trees", "tsz"} or (
        input_format == "auto" and Path(input_path).suffix.lower() in TREE_SEQUENCE_SUFFIXES
    )


def _vcf_producer_command(input_path: str | Path, input_format: str) -> list[str]:
    script = (
        "import sys; "
        "from gamma_smc_aou.tree_sequence import stream_tree_sequence_vcf; "
        "stream_tree_sequence_vcf(sys.argv[1], sys.stdout, input_format=sys.argv[2])"
    )
    return [sys.executable, "-c", script, str(input_path), input_format]


def _decode_failure(
    command: list[str],
    completed: subprocess.CompletedProcess[str],
    *,
    producer_command: list[str] | None,
    producer_returncode: int,
    producer_stderr: str,
) -> RuntimeError:
    parts = [
        "Gamma-SMC decoder pipeline failed",
        f"decoder return code: {completed.returncode}",
        f"decoder command: {command!r}",
    ]
    if completed.stdout:
        parts.extend(["decoder stdout:", completed.stdout.rstrip()])
    if completed.stderr:
        parts.extend(["decoder stderr:", completed.stderr.rstrip()])
    if producer_command is not None:
        parts.extend(
            [
                f"tree-sequence VCF producer return code: {producer_returncode}",
                f"tree-sequence VCF producer command: {producer_command!r}",
            ]
        )
        if producer_stderr:
            parts.extend(["tree-sequence VCF producer stderr:", producer_stderr.rstrip()])
    return RuntimeError("\n".join(parts))


def run_within_decoder(
    executable: str | Path,
    input_path: str | Path,
    output_summary: str | Path,
    *,
    scaled_mutation_rate: float = DEFAULT_SCALED_MUTATION_RATE,
    recombination_to_mutation_ratio: float = DEFAULT_RECOMBINATION_TO_MUTATION_RATIO,
    mutation_rate: float = DEFAULT_MUTATION_RATE,
    threshold_years: float | Sequence[float] = 4500,
    generation_time: float = DEFAULT_GENERATION_TIME,
    input_format: str = "auto",
    raw_output: str | Path | None = None,
    bitmatrix_output: str | Path | None = None,
    mask: str | Path | None = None,
    masks_per_sample: str | Path | None = None,
    samples: str | Path | None = None,
    output_positions_file: str | Path | None = None,
    output_at_stride: int = DEFAULT_OUTPUT_STRIDE,
    output_at_hets: bool = False,
    only_within: bool = True,
    n_random_pairs: int = 0,
    pairs_seed: int = 1729,
    pairs_file: str | Path | None = None,
    pairs_manifest: str | Path | None = None,
    exclude_within: bool = False,
    recent_call: str = "median",
    recent_call_probability: float = 0.5,
    threads: int = 0,
    cache_size: int = DEFAULT_CACHE_SIZE,
    pair_block: int = 256,
    exp10: str = "accurate",
    backward_alignment: str = "fixed",
    extra_args: Sequence[str] | None = None,
) -> dict:
    """Run the decoder and write a per-position recent-coalescence summary.

    ``threshold_years`` accepts several thresholds; the summary gains an
    ``n_recent_<years>``/``frac_recent_<years>``/``mean_p_lt_<years>`` block for
    each, and ``mean_p_tmrca_lt_threshold`` stays as the first threshold's mean
    probability so existing consumers keep working.

    Pair selection is mutually exclusive: ``only_within`` (one pair per
    diploid), ``n_random_pairs`` (a uniform sample over all haplotype pairs,
    within-individual pairs included unless ``exclude_within``), or
    ``pairs_file``.
    """
    output_summary = Path(output_summary)
    output_summary.parent.mkdir(parents=True, exist_ok=True)

    thresholds = (
        [threshold_years]
        if isinstance(threshold_years, (int, float))
        else list(threshold_years)
    )
    if not thresholds:
        raise ValueError("at least one threshold is required")

    selectors = sum([bool(only_within), n_random_pairs > 0, pairs_file is not None])
    if selectors > 1:
        raise ValueError("only_within, n_random_pairs and pairs_file are mutually exclusive")
    if selectors == 0:
        raise ValueError(
            "choose a pair set: only_within, n_random_pairs or pairs_file "
            "(the exhaustive default is O(n^2) haplotype pairs)"
        )

    streams_tree_sequence = _streams_tree_sequence(input_path, input_format)
    decoder_input = "/dev/stdin" if streams_tree_sequence else str(input_path)
    decoder_input_format = "vcf" if streams_tree_sequence else input_format
    command = [
        str(executable), "--input", decoder_input, "--input_format", decoder_input_format,
        "--scaled_mutation_rate", str(scaled_mutation_rate),
        "--recombination_to_mutation_ratio", str(recombination_to_mutation_ratio),
        "--unscaled_mutation_rate", str(mutation_rate),
        "--recent_threshold_years", ",".join(str(value) for value in thresholds),
        "--generation_time", str(generation_time),
        "--recent_summary", str(output_summary),
        "--recent_call", recent_call,
        "--recent_call_probability", str(recent_call_probability),
        # cxxopts boolean options use an implicit ``true`` value when the flag
        # is a standalone token.  Keep the value attached so ``false`` is not
        # misread as an unused positional argument.
        f"--output_at_hets={str(output_at_hets).lower()}",
        "--output_at_stride", str(output_at_stride),
        "--cache_size", str(cache_size),
        "--threads", str(threads),
        "--pair_block", str(pair_block),
        "--backward_alignment", backward_alignment,
        "--exp10", exp10,
    ]
    if only_within:
        command.append("--only_within")
    if n_random_pairs > 0:
        command.extend(["--n_random_pairs", str(n_random_pairs), "--pairs_seed", str(pairs_seed)])
        if exclude_within:
            command.append("--exclude_within")
    if pairs_file is not None:
        command.extend(["--pairs_file", str(pairs_file)])
    if pairs_manifest is not None:
        Path(pairs_manifest).parent.mkdir(parents=True, exist_ok=True)
        command.extend(["--pairs_manifest", str(pairs_manifest)])
    if raw_output is not None:
        command.extend(["--output", str(raw_output)])
    if bitmatrix_output is not None:
        Path(bitmatrix_output).parent.mkdir(parents=True, exist_ok=True)
        command.extend(["--recent_bitmatrix", str(bitmatrix_output)])
    if mask is not None:
        command.extend(["--mask", str(mask)])
    if masks_per_sample is not None:
        command.extend(["--masks_per_sample", str(masks_per_sample)])
    if samples is not None:
        command.extend(["--samples", str(samples)])
    if output_positions_file is not None:
        command.extend(["--output_positions_file", str(output_positions_file)])
    if extra_args:
        command.extend(str(value) for value in extra_args)

    started = perf_counter()
    producer_command = (
        _vcf_producer_command(input_path, input_format) if streams_tree_sequence else None
    )
    producer = None
    producer_returncode = 0
    producer_stderr = ""
    if producer_command is not None:
        producer = subprocess.Popen(
            producer_command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if producer.stdout is None or producer.stderr is None:
            raise RuntimeError("could not open tree-sequence VCF streaming pipes")
    try:
        completed = subprocess.run(
            command,
            check=False,
            text=True,
            capture_output=True,
            stdin=producer.stdout if producer is not None else None,
        )
    finally:
        if producer is not None:
            producer.stdout.close()
            producer_stderr = producer.stderr.read().decode("utf-8", errors="replace")
            producer_returncode = producer.wait()
    decode_seconds = perf_counter() - started
    if completed.returncode != 0 or producer_returncode != 0:
        raise _decode_failure(
            command,
            completed,
            producer_command=producer_command,
            producer_returncode=producer_returncode,
            producer_stderr=producer_stderr,
        )
    with output_summary.open(encoding="utf-8") as handle:
        n_output_positions = max(0, sum(1 for _ in handle) - 1)

    # Record where the draw was written. Without an explicit --pairs_manifest
    # the binary derives one next to the output whenever it sampled at random,
    # so resolve the same way rather than reporting nothing.
    manifest = Path(pairs_manifest) if pairs_manifest is not None else None
    if manifest is None and n_random_pairs > 0:
        anchor = output_summary or bitmatrix_output or raw_output
        if anchor is not None:
            manifest = Path(str(anchor) + ".pairs.tsv")

    run = {
        "command": command,
        "input_path": str(input_path),
        "input_format": input_format,
        "tree_sequence_vcf_producer_command": producer_command,
        "stdout": completed.stdout,
        "stderr": completed.stderr,
        "decode_seconds": float(decode_seconds),
        "stride_bp": int(output_at_stride),
        "cache_size_bp": int(cache_size),
        "n_output_positions": int(n_output_positions),
        "pairs_manifest": str(manifest) if manifest is not None else None,
        "n_pairs_recorded": (
            sum(
                1 for line in manifest.read_text().splitlines()
                if line.strip() and not line.startswith("#")
            )
            if manifest is not None and manifest.exists() else None
        ),
    }
    with output_summary.with_suffix(output_summary.suffix + ".run.json").open("w", encoding="utf-8") as handle:
        json.dump(run, handle, indent=2)
    return run
