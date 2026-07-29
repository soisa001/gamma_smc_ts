from __future__ import annotations

import json
import subprocess
from collections.abc import Sequence
from pathlib import Path


def run_within_decoder(
    executable: str | Path,
    input_path: str | Path,
    output_summary: str | Path,
    *,
    scaled_mutation_rate: float = 0.00075,
    recombination_to_mutation_ratio: float = 0.8,
    mutation_rate: float = 1.29e-8,
    threshold_years: float | Sequence[float] = 4500,
    generation_time: float = 25,
    input_format: str = "auto",
    raw_output: str | Path | None = None,
    bitmatrix_output: str | Path | None = None,
    mask: str | Path | None = None,
    masks_per_sample: str | Path | None = None,
    output_at_stride: int = -1,
    output_at_hets: bool = True,
    only_within: bool = True,
    n_random_pairs: int = 0,
    pairs_seed: int = 1729,
    pairs_file: str | Path | None = None,
    pairs_manifest: str | Path | None = None,
    exclude_within: bool = False,
    recent_call: str = "median",
    recent_call_probability: float = 0.5,
    threads: int = 0,
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

    command = [
        str(executable), "--input", str(input_path), "--input_format", input_format,
        "--scaled_mutation_rate", str(scaled_mutation_rate),
        "--recombination_to_mutation_ratio", str(recombination_to_mutation_ratio),
        "--unscaled_mutation_rate", str(mutation_rate),
        "--recent_threshold_years", ",".join(str(value) for value in thresholds),
        "--generation_time", str(generation_time),
        "--recent_summary", str(output_summary),
        "--recent_call", recent_call,
        "--recent_call_probability", str(recent_call_probability),
        "--output_at_hets", str(output_at_hets).lower(),
        "--output_at_stride", str(output_at_stride),
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
    if extra_args:
        command.extend(str(value) for value in extra_args)

    completed = subprocess.run(command, check=True, text=True, capture_output=True)

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
        "stdout": completed.stdout,
        "stderr": completed.stderr,
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
