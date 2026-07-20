from __future__ import annotations

import json
import subprocess
from pathlib import Path


def run_within_decoder(
    executable: str | Path,
    input_path: str | Path,
    output_summary: str | Path,
    *,
    scaled_mutation_rate: float,
    recombination_to_mutation_ratio: float,
    mutation_rate: float,
    threshold_years: float = 4500,
    generation_time: float = 30,
    input_format: str = "auto",
    raw_output: str | Path | None = None,
    mask: str | Path | None = None,
    masks_per_sample: str | Path | None = None,
    output_at_stride: int = -1,
    output_at_hets: bool = True,
) -> dict:
    output_summary = Path(output_summary)
    output_summary.parent.mkdir(parents=True, exist_ok=True)
    command = [
        str(executable), "--input", str(input_path), "--input_format", input_format,
        "--only_within", "--scaled_mutation_rate", str(scaled_mutation_rate),
        "--recombination_to_mutation_ratio", str(recombination_to_mutation_ratio),
        "--unscaled_mutation_rate", str(mutation_rate),
        "--recent_threshold_years", str(threshold_years),
        "--generation_time", str(generation_time),
        "--recent_summary", str(output_summary),
        "--output_at_hets", str(output_at_hets).lower(),
        "--output_at_stride", str(output_at_stride),
    ]
    if raw_output is not None:
        command.extend(["--output", str(raw_output)])
    if mask is not None:
        command.extend(["--mask", str(mask)])
    if masks_per_sample is not None:
        command.extend(["--masks_per_sample", str(masks_per_sample)])
    completed = subprocess.run(command, check=True, text=True, capture_output=True)
    run = {"command": command, "stdout": completed.stdout, "stderr": completed.stderr}
    with output_summary.with_suffix(output_summary.suffix + ".run.json").open("w", encoding="utf-8") as handle:
        json.dump(run, handle, indent=2)
    return run
