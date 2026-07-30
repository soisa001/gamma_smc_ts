from __future__ import annotations

import json
import re
import shutil
import subprocess
from pathlib import Path
from time import perf_counter

import numpy as np
import pandas as pd
import zstandard
from scipy.special import gammainc

from .defaults import (
    DEFAULT_GENERATION_TIME,
    DEFAULT_MUTATION_RATE,
    DEFAULT_OUTPUT_STRIDE,
    DEFAULT_RECOMBINATION_TO_MUTATION_RATIO,
    DEFAULT_SCALED_MUTATION_RATE,
)


DEFAULT_IMAGE = "docker.io/regevsch/gamma_smc:v0.2"
ANSI_ESCAPE = re.compile(r"\x1b\[[0-?]*[ -/]*[@-~]")


def find_container_runtime(requested: str = "auto") -> str:
    """Find a supported Gamma-SMC container runtime."""
    choices = ("apptainer", "singularity", "docker")
    if requested != "auto":
        if requested not in choices:
            raise ValueError(f"unsupported container runtime: {requested}")
        if shutil.which(requested) is None:
            raise FileNotFoundError(f"{requested} is not on PATH")
        return requested
    for candidate in choices:
        if shutil.which(candidate) is not None:
            return candidate
    raise FileNotFoundError(
        "Gamma-SMC v0.2 needs Apptainer/Singularity or Docker; none was found"
    )


def _container_image(runtime: str, image: str) -> str:
    if runtime in ("apptainer", "singularity") and not image.startswith("docker://"):
        return f"docker://{image}"
    return image


def build_container_command(
    *,
    runtime: str,
    image: str,
    work_dir: str | Path,
    input_name: str,
    output_name: str,
    scaled_mutation_rate: float,
    recombination_to_mutation_ratio: float,
    stride: int,
) -> list[str]:
    """Build the official v0.2 invocation using one mounted work directory."""
    if stride < 1:
        raise ValueError("stride must be positive")
    work_dir = str(Path(work_dir).resolve())
    arguments = [
        "--input", f"/work/{input_name}",
        "--output", f"/work/{output_name}",
        "--only_within",
        "--scaled_mutation_rate", str(scaled_mutation_rate),
        "--recombination_to_mutation_ratio", str(recombination_to_mutation_ratio),
        "--output_at_hets=false",
        "--output_at_stride", str(stride),
    ]
    if runtime == "docker":
        return [
            "docker", "run", "--rm", "-v", f"{work_dir}:/work",
            _container_image(runtime, image), *arguments,
        ]
    return [
        runtime, "run", "--bind", f"{work_dir}:/work",
        _container_image(runtime, image), *arguments,
    ]


def _read_exact(reader, n_bytes: int) -> bytes:
    chunks: list[bytes] = []
    remaining = n_bytes
    while remaining:
        chunk = reader.read(remaining)
        if not chunk:
            raise EOFError(f"posterior stream ended {remaining} bytes early")
        chunks.append(chunk)
        remaining -= len(chunk)
    return b"".join(chunks)


def _concise_log(value: str) -> str:
    """Keep informational lines while dropping animated progress-bar frames."""
    clean = ANSI_ESCAPE.sub("", value).replace("\r", "\n")
    lines = [
        line.rstrip()
        for line in clean.splitlines()
        if line.strip() and not line.lstrip().startswith("[")
    ]
    return "\n".join(lines[-200:])


def summarize_posteriors(
    posterior_path: str | Path,
    *,
    mutation_rate: float,
    threshold_years: float,
    generation_time: float,
) -> pd.DataFrame:
    """Stream a v0.2 posterior and average P(TMRCA<threshold) over pairs."""
    posterior_path = Path(posterior_path)
    with posterior_path.with_name(posterior_path.name + ".meta").open(
        encoding="utf-8"
    ) as handle:
        meta = json.load(handle)
    n_pairs = int(meta["num_pairs"])
    chunk_size = int(meta["chunk_size"])
    n_positions = int(meta["sequence_length"])
    positions = np.asarray(meta["output_positions"], dtype=np.int64)
    if len(positions) != n_positions:
        raise ValueError("metadata output_positions length is inconsistent")
    theta = float(meta["scaled_mutation_rate"])
    if theta <= 0 or mutation_rate <= 0 or generation_time <= 0:
        raise ValueError("rates and generation time must be positive")
    two_ne_generations = theta / (2 * mutation_rate)
    threshold_scaled = (threshold_years / generation_time) / two_ne_generations
    sum_probability = np.zeros(n_positions, dtype=np.float64)
    sum_tmrca_generations = np.zeros(n_positions, dtype=np.float64)
    n_valid = np.zeros(n_positions, dtype=np.int64)
    block_bytes = n_positions * chunk_size * np.dtype(np.float32).itemsize

    with posterior_path.open("rb") as compressed:
        with zstandard.ZstdDecompressor().stream_reader(compressed) as reader:
            for starting_pair in range(0, n_pairs, chunk_size):
                real_pairs = min(chunk_size, n_pairs - starting_pair)
                alpha = np.frombuffer(
                    _read_exact(reader, block_bytes), dtype=np.float32
                ).reshape(n_positions, chunk_size)[:, :real_pairs]
                beta = np.frombuffer(
                    _read_exact(reader, block_bytes), dtype=np.float32
                ).reshape(n_positions, chunk_size)[:, :real_pairs]
                valid = (alpha > 0) & (beta > 0) & np.isfinite(alpha) & np.isfinite(beta)
                probability = np.zeros_like(alpha, dtype=np.float64)
                tmrca = np.zeros_like(alpha, dtype=np.float64)
                probability[valid] = gammainc(
                    alpha[valid].astype(np.float64),
                    beta[valid].astype(np.float64) * threshold_scaled,
                )
                tmrca[valid] = (
                    alpha[valid].astype(np.float64)
                    / beta[valid].astype(np.float64)
                    * two_ne_generations
                )
                sum_probability += probability.sum(axis=1)
                sum_tmrca_generations += tmrca.sum(axis=1)
                n_valid += valid.sum(axis=1)

    denom = np.where(n_valid > 0, n_valid, 1)
    return pd.DataFrame({
        "position_0based": positions,
        "position_1based": positions + 1,
        "n_pairs": n_valid,
        "mean_p_tmrca_lt_threshold": sum_probability / denom,
        "mean_tmrca_generations": sum_tmrca_generations / denom,
    })


def run_container_decoder(
    input_path: str | Path,
    output_summary: str | Path,
    *,
    scaled_mutation_rate: float = DEFAULT_SCALED_MUTATION_RATE,
    recombination_to_mutation_ratio: float = DEFAULT_RECOMBINATION_TO_MUTATION_RATIO,
    mutation_rate: float = DEFAULT_MUTATION_RATE,
    threshold_years: float = 4500,
    generation_time: float = DEFAULT_GENERATION_TIME,
    stride: int = DEFAULT_OUTPUT_STRIDE,
    runtime: str = "auto",
    image: str = DEFAULT_IMAGE,
    keep_raw: bool = False,
) -> dict:
    """Run official Gamma-SMC v0.2 and reduce its output to a small TSV."""
    input_path = Path(input_path).resolve()
    output_summary = Path(output_summary).resolve()
    output_summary.parent.mkdir(parents=True, exist_ok=True)
    if input_path.parent != output_summary.parent:
        raise ValueError("input and output summary must share one container work directory")
    selected_runtime = find_container_runtime(runtime)
    posterior_path = output_summary.with_suffix(".posteriors.zst")
    command = build_container_command(
        runtime=selected_runtime,
        image=image,
        work_dir=output_summary.parent,
        input_name=input_path.name,
        output_name=posterior_path.name,
        scaled_mutation_rate=scaled_mutation_rate,
        recombination_to_mutation_ratio=recombination_to_mutation_ratio,
        stride=stride,
    )
    started = perf_counter()
    completed = subprocess.run(command, check=True, text=True, capture_output=True)
    decode_seconds = perf_counter() - started
    started_summary = perf_counter()
    profile = summarize_posteriors(
        posterior_path,
        mutation_rate=mutation_rate,
        threshold_years=threshold_years,
        generation_time=generation_time,
    )
    summary_seconds = perf_counter() - started_summary
    profile.to_csv(output_summary, sep="\t", index=False)
    raw_bytes = posterior_path.stat().st_size
    metadata_path = posterior_path.with_name(posterior_path.name + ".meta")
    run = {
        "command": command,
        "runtime": selected_runtime,
        "image": image,
        "stride_bp": int(stride),
        "n_output_positions": int(len(profile)),
        "n_pairs": int(profile["n_pairs"].max()),
        "decode_seconds": float(decode_seconds),
        "summary_seconds": float(summary_seconds),
        "posterior_bytes": int(raw_bytes),
        "stdout_summary": _concise_log(completed.stdout),
        "stderr_summary": _concise_log(completed.stderr),
    }
    with output_summary.with_name(output_summary.name + ".run.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(run, handle, indent=2)
    if not keep_raw:
        posterior_path.unlink()
        metadata_path.unlink(missing_ok=True)
    return run
