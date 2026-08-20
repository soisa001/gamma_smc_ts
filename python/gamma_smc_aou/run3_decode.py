"""Gamma-SMC decoding for the run3 study.

Runs the decoder over the same tree sequences the truth analysis reads, with the
same within-individual pairs, the same 10 kb output stride and the same
threshold grid, so the decoded and true profiles are directly comparable
position by position and pair set by pair set.

The decoder emits two per-position quantities per threshold:

* ``frac_recent_<years>`` -- the fraction of pairs whose *called* TMRCA is below
  the threshold. This is the direct analogue of the tree-truth statistic.
* ``mean_p_lt_<years>`` -- the mean posterior probability that the TMRCA is
  below the threshold. Softer, and it uses the whole posterior rather than a
  point call.

Both are carried through, because which one the decoder does better on is part
of what run3 is measuring.
"""

from __future__ import annotations

import json
import os
import traceback
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter
from typing import Any, Sequence

import pandas as pd

from .decoder import run_within_decoder
from .run3_config import (
    DECODER_SETTINGS,
    MODES,
    Run3Arm,
    TMRCA_THRESHOLDS_YEARS,
    get_arm,
    scaled_mutation_rate,
)
from .run3_simulate import REPLICATE_OUTPUTS

DECODE_OUTPUTS = {
    "summary": "gamma_smc_summary.tsv",
    "receipt": "gamma_smc_receipt.json",
}


@dataclass(frozen=True)
class DecodeTask:
    arm_id: str
    mode: str
    replicate_index: int
    directory: str
    decoder_path: str
    threads: int


def build_decode_tasks(
    arm: Run3Arm,
    mode: str,
    *,
    study_root: str | Path,
    decoder_path: str | Path,
    replicates: int,
    threads: int = 1,
    indices: Sequence[int] | None = None,
) -> list[DecodeTask]:
    if mode not in MODES:
        raise ValueError(f"mode must be one of {MODES}, got {mode!r}")
    chosen = list(range(replicates)) if indices is None else list(indices)
    base = Path(study_root) / arm.arm_id / mode / "replicates"
    return [
        DecodeTask(
            arm_id=arm.arm_id,
            mode=mode,
            replicate_index=index,
            directory=str(base / f"rep{index:03d}"),
            decoder_path=str(decoder_path),
            threads=threads,
        )
        for index in chosen
    ]


def decode_replicate(task: DecodeTask) -> dict[str, Any]:
    """Decode one replicate's tree sequence. Never raises."""

    directory = Path(task.directory)
    started = perf_counter()
    receipt: dict[str, Any]
    try:
        arm = get_arm(task.arm_id)
        tree_path = directory / REPLICATE_OUTPUTS["tree"]
        if not tree_path.is_file():
            raise FileNotFoundError(f"no tree sequence at {tree_path}")
        summary_path = directory / DECODE_OUTPUTS["summary"]
        settings = dict(DECODER_SETTINGS)
        settings["threads"] = task.threads
        result = run_within_decoder(
            task.decoder_path,
            tree_path,
            summary_path,
            scaled_mutation_rate=scaled_mutation_rate(arm),
            threshold_years=list(TMRCA_THRESHOLDS_YEARS),
            generation_time=arm.generation_time,
            **settings,
        )
        frame = pd.read_csv(summary_path, sep="\t")
        receipt = {
            "schema": "gamma-smc.run3-decode/v1",
            "arm_id": task.arm_id,
            "mode": task.mode,
            "replicate_index": task.replicate_index,
            "status": "completed",
            "elapsed_seconds": perf_counter() - started,
            "scaled_mutation_rate": scaled_mutation_rate(arm),
            "generation_time_years": arm.generation_time,
            "thresholds_years": list(TMRCA_THRESHOLDS_YEARS),
            "settings": settings,
            "summary_rows": int(len(frame)),
            "summary_columns": list(frame.columns),
            "decoder_result": result,
            "summary": str(summary_path),
        }
    except Exception as error:  # noqa: BLE001
        receipt = {
            "schema": "gamma-smc.run3-decode/v1",
            "arm_id": task.arm_id,
            "mode": task.mode,
            "replicate_index": task.replicate_index,
            "status": "failed",
            "elapsed_seconds": perf_counter() - started,
            "error": f"{type(error).__name__}: {error}",
            "traceback": traceback.format_exc(),
        }

    directory.mkdir(parents=True, exist_ok=True)
    path = directory / DECODE_OUTPUTS["receipt"]
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    temporary.write_text(
        json.dumps(receipt, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)
    return {
        "arm_id": receipt["arm_id"],
        "mode": receipt["mode"],
        "replicate_index": receipt["replicate_index"],
        "status": receipt["status"],
        "elapsed_seconds": receipt["elapsed_seconds"],
        "error": receipt.get("error", ""),
    }


def run_decode_tasks(tasks: Sequence[DecodeTask], *, workers: int = 1) -> pd.DataFrame:
    columns = ["arm_id", "mode", "replicate_index", "status", "elapsed_seconds", "error"]
    if not tasks:
        return pd.DataFrame(columns=columns)
    if workers == 1:
        rows = [decode_replicate(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            rows = list(executor.map(decode_replicate, tasks))
    return pd.DataFrame(rows)


__all__ = [
    "DECODE_OUTPUTS",
    "DecodeTask",
    "build_decode_tasks",
    "decode_replicate",
    "run_decode_tasks",
]
