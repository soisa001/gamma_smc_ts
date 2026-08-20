"""Replicate execution for the run4 study.

The important difference from run2 is that **raw per-pair TMRCA values are
stored**.  run2 summarised to a fixed threshold grid at simulation time, which
meant that discovering the grid was wrong cost a full 23-core-hour re-run.  Here
each replicate writes the complete ``(n_positions, n_pairs)`` TMRCA matrix as
float32 plus the exact focal-base column as float64, so any threshold, window,
genotype stratification or summary is a post-hoc calculation.

Tree sequences are also kept, because the Gamma-SMC decode phase consumes them.
"""

from __future__ import annotations

import json
import os
import platform
import traceback
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd
import stdpopsim
import tskit

from .run4_config import (
    FOCAL_POSITION_BP,
    MODES,
    PROFILE_STRIDE_BP,
    Run4Arm,
    SAMPLE_DIPLOIDS,
    SEQUENCE_LENGTH_BP,
    SLIM_BURN_IN,
    SLIM_SCALING_FACTOR,
    get_arm,
)
from .run4_models import (
    build_contig,
    build_demographic_model,
    build_extended_events,
    model_record,
    population_id,
    scoped_focal_overlay_patch,
    scoped_slim_patch,
)

GENOTYPE_CLASSES = ("overall", "hom_carrier", "het", "hom_noncarrier", "any_carrier")

REPLICATE_OUTPUTS = {
    "endpoint": "endpoint.json",
    "tmrca_matrix": "tmrca_matrix.npz",
    "focal_tmrca": "focal_tmrca.npy",
    "focal_genotypes": "focal_genotypes.npy",
    "tree": "simulation.trees",
}


# ---------------------------------------------------------------------------
# Tree-sequence helpers
# ---------------------------------------------------------------------------


def sampled_diploid_pairs(ts: tskit.TreeSequence) -> np.ndarray:
    """Return the ``(n, 2)`` array of within-individual sample node pairs.

    The order matches ``TreeSequence.write_vcf`` individual order, which is the
    order Gamma-SMC assigns to consecutive haplotype pairs, so truth and decoded
    profiles are directly comparable pair by pair.
    """
    sample_nodes = set(int(node) for node in ts.samples())
    pairs: list[tuple[int, int]] = []
    for individual in ts.individuals():
        nodes = [int(node) for node in individual.nodes if int(node) in sample_nodes]
        if not nodes:
            continue
        if len(nodes) != 2:
            raise ValueError(
                f"individual {individual.id} has {len(nodes)} sample nodes; expected 2"
            )
        pairs.append((nodes[0], nodes[1]))
    if not pairs:
        raise ValueError("tree sequence has no sampled diploid individuals")
    flat = [node for pair in pairs for node in pair]
    if len(set(flat)) != len(flat) or set(flat) != sample_nodes:
        raise ValueError("sample nodes do not partition cleanly into diploids")
    return np.asarray(pairs, dtype=np.int64)


def focal_allele_counts(
    ts: tskit.TreeSequence,
    pairs: np.ndarray,
    focal_position: int = FOCAL_POSITION_BP,
) -> np.ndarray:
    """Return per-diploid focal allele counts; all zero when the site is absent.

    In run4 the allele exists in both modes, but archaic ancestry is frequently
    absent at any single base, so an empty focal site is a common and legitimate
    outcome rather than a failure.
    """
    matches = [
        site
        for site in ts.sites()
        if np.isclose(float(site.position), float(focal_position), rtol=0.0, atol=1e-9)
    ]
    if not matches:
        return np.zeros(len(pairs), dtype=np.int8)
    if len(matches) > 1:
        raise ValueError(f"found {len(matches)} sites at the reserved focal base")
    ordered = [int(node) for pair in pairs for node in pair]
    variant = tskit.Variant(ts, samples=ordered)
    variant.decode(matches[0].id)
    genotypes = np.asarray(variant.genotypes, dtype=np.int64)
    if np.any(genotypes < 0):
        raise ValueError("missing focal genotypes cannot be classified")
    return (genotypes.reshape(len(pairs), 2) != 0).sum(axis=1).astype(np.int8)


def profile_positions(
    sequence_length: int = SEQUENCE_LENGTH_BP,
    stride: int = PROFILE_STRIDE_BP,
    focal_position: int = FOCAL_POSITION_BP,
) -> np.ndarray:
    grid = np.arange(0, sequence_length, stride, dtype=float)
    return np.unique(np.concatenate([grid, [float(focal_position)]]))


def pairwise_tmrca_matrix(
    ts: tskit.TreeSequence, pairs: np.ndarray, positions: np.ndarray
) -> np.ndarray:
    """Return the ``(n_positions, n_pairs)`` within-diploid TMRCA matrix."""
    values = np.full((len(positions), len(pairs)), np.nan, dtype=float)
    index = 0
    for tree in ts.trees():
        _, right = tree.interval
        while index < len(positions) and positions[index] < right:
            values[index] = [
                tree.tmrca(int(a), int(b)) for a, b in pairs
            ]
            index += 1
        if index >= len(positions):
            break
    if index != len(positions) or not np.all(np.isfinite(values)):
        raise ValueError("failed to evaluate a finite TMRCA at every profile position")
    return values


# ---------------------------------------------------------------------------
# Replicate execution
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ReplicateTask:
    arm_id: str
    mode: str
    replicate_index: int
    seed: int
    unit_id: str
    directory: str
    repo_root: str
    slim_path: str


def build_tasks(
    arm: Run4Arm,
    mode: str,
    *,
    study_root: str | Path,
    repo_root: str | Path,
    slim_path: str | Path,
    replicates: int,
    indices: Sequence[int] | None = None,
) -> list[ReplicateTask]:
    if mode not in MODES:
        raise ValueError(f"mode must be one of {MODES}, got {mode!r}")
    chosen = list(range(replicates)) if indices is None else list(indices)
    base = Path(study_root) / arm.arm_id / mode / "replicates"
    return [
        ReplicateTask(
            arm_id=arm.arm_id,
            mode=mode,
            replicate_index=index,
            seed=arm.seed(mode, index),
            unit_id=arm.unit_id(mode, index),
            directory=str(base / f"rep{index:03d}"),
            repo_root=str(repo_root),
            slim_path=str(slim_path),
        )
        for index in chosen
    ]


def _runtime_versions() -> dict[str, str]:
    return {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "stdpopsim": stdpopsim.__version__,
        "tskit": tskit.__version__,
        "numpy": np.__version__,
        "pandas": pd.__version__,
    }


def _metadata_block(ts: tskit.TreeSequence) -> Mapping[str, Any]:
    metadata = ts.metadata
    slim = metadata.get("SLiM") if isinstance(metadata, Mapping) else None
    user = slim.get("user_metadata") if isinstance(slim, Mapping) else None
    return user if isinstance(user, Mapping) else {}


def _scalar(block: Mapping[str, Any], key: str) -> Any:
    if key not in block:
        return None
    value = block[key]
    while isinstance(value, (list, tuple)) and len(value) == 1:
        value = value[0]
    return value


def run_replicate(task: ReplicateTask) -> dict[str, Any]:
    """Run one replicate and write its raw outputs. Never raises."""

    directory = Path(task.directory)
    directory.mkdir(parents=True, exist_ok=True)
    started = perf_counter()
    try:
        arm = get_arm(task.arm_id)
        model = build_demographic_model(arm, task.repo_root)
        contig = build_contig()
        events = build_extended_events(arm, task.mode)
        target_id = population_id(model, arm.target_population)
        engine = stdpopsim.get_engine("slim")

        with (
            scoped_slim_patch(arm, task.mode, target_id) as slim_patch,
            scoped_focal_overlay_patch() as overlay,
        ):
            ts = engine.simulate(
                model,
                contig,
                {arm.target_population: SAMPLE_DIPLOIDS},
                seed=task.seed,
                extended_events=list(events),
                slim_path=task.slim_path,
                slim_scaling_factor=SLIM_SCALING_FACTOR,
                slim_burn_in=SLIM_BURN_IN,
                keep_mutation_ids_as_alleles=False,
            )
        if overlay["masked_calls"] < 1:
            raise RuntimeError("the reserved-base guard never fired")

        pairs = sampled_diploid_pairs(ts)
        if len(pairs) != SAMPLE_DIPLOIDS:
            raise RuntimeError(f"expected {SAMPLE_DIPLOIDS} diploids, got {len(pairs)}")
        counts = focal_allele_counts(ts, pairs)
        positions = profile_positions()
        tmrca = pairwise_tmrca_matrix(ts, pairs, positions)

        focal_index = int(np.searchsorted(positions, float(FOCAL_POSITION_BP)))
        if positions[focal_index] != float(FOCAL_POSITION_BP):
            raise RuntimeError("the focal base is missing from the profile grid")

        # Raw outputs: everything downstream is recomputable from these.
        ts.dump(directory / REPLICATE_OUTPUTS["tree"])
        np.savez_compressed(
            directory / REPLICATE_OUTPUTS["tmrca_matrix"],
            positions=positions.astype(np.float64),
            tmrca_generations=tmrca.astype(np.float32),
            pairs=pairs.astype(np.int32),
            focal_index=np.int64(focal_index),
        )
        np.save(directory / REPLICATE_OUTPUTS["focal_tmrca"], tmrca[focal_index])
        np.save(directory / REPLICATE_OUTPUTS["focal_genotypes"], counts)

        block = _metadata_block(ts)
        sample_alt = int(counts.sum())
        endpoint: dict[str, Any] = {
            "schema": "gamma-smc.run4-endpoint/v1",
            "unit_id": task.unit_id,
            "arm_id": task.arm_id,
            "mode": task.mode,
            "replicate_index": task.replicate_index,
            "seed": task.seed,
            "status": "completed",
            "elapsed_seconds": perf_counter() - started,
            "generation_time_years": arm.generation_time,
            "model": model_record(arm, task.mode, task.repo_root),
            "slim_patch": slim_patch,
            "focal_overlay_guard": dict(overlay),
            "runtime_versions": _runtime_versions(),
            "tree_sequence": {
                "num_samples": int(ts.num_samples),
                "num_sites": int(ts.num_sites),
                "num_trees": int(ts.num_trees),
                "sequence_length": float(ts.sequence_length),
            },
            "placement": {
                "tick": _scalar(block, "run4_placement_tick"),
                "population_id": _scalar(block, "run4_placement_population_id"),
                "mode": _scalar(block, "run4_placement_mode"),
                "carrier_genomes": _scalar(block, "run4_placement_carrier_genomes"),
                "population_genomes": _scalar(block, "run4_placement_population_genomes"),
                "frequency": _scalar(block, "run4_placement_frequency"),
            },
            "final_allele_frequency": {
                "census_alt_count": _scalar(block, "run4_final_census_alt_count"),
                "census_total_count": _scalar(block, "run4_final_census_total_count"),
                "census_af": _scalar(block, "run4_final_census_af"),
                "census_state": _scalar(block, "run4_final_census_state"),
                "sample_alt_count": sample_alt,
                "sample_total_count": 2 * len(pairs),
                "sample_af": sample_alt / (2 * len(pairs)),
                "sample_carrier_diploids": int((counts > 0).sum()),
                "sample_hom_carrier_diploids": int((counts == 2).sum()),
                "archaic_present_at_focal_base": bool(sample_alt > 0),
            },
            "focal_tmrca_summary": {
                "n_pairs": int(tmrca.shape[1]),
                "mean_generations": float(tmrca[focal_index].mean()),
                "median_generations": float(np.median(tmrca[focal_index])),
                "min_generations": float(tmrca[focal_index].min()),
                "max_generations": float(tmrca[focal_index].max()),
            },
            "outputs": {
                key: str(directory / name) for key, name in REPLICATE_OUTPUTS.items()
            },
        }
    except Exception as error:  # noqa: BLE001
        endpoint = {
            "schema": "gamma-smc.run4-endpoint/v1",
            "unit_id": task.unit_id,
            "arm_id": task.arm_id,
            "mode": task.mode,
            "replicate_index": task.replicate_index,
            "seed": task.seed,
            "status": "failed",
            "elapsed_seconds": perf_counter() - started,
            "error": f"{type(error).__name__}: {error}",
            "traceback": traceback.format_exc(),
        }

    path = directory / REPLICATE_OUTPUTS["endpoint"]
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    temporary.write_text(
        json.dumps(endpoint, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)
    return {
        "unit_id": endpoint["unit_id"],
        "arm_id": endpoint["arm_id"],
        "mode": endpoint["mode"],
        "replicate_index": endpoint["replicate_index"],
        "seed": endpoint["seed"],
        "status": endpoint["status"],
        "elapsed_seconds": endpoint["elapsed_seconds"],
        "error": endpoint.get("error", ""),
    }


def run_tasks(tasks: Sequence[ReplicateTask], *, workers: int = 1) -> pd.DataFrame:
    if workers < 1:
        raise ValueError("workers must be positive")
    columns = [
        "unit_id", "arm_id", "mode", "replicate_index",
        "seed", "status", "elapsed_seconds", "error",
    ]
    if not tasks:
        return pd.DataFrame(columns=columns)
    if workers == 1:
        rows = [run_replicate(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            rows = list(executor.map(run_replicate, tasks))
    return pd.DataFrame(rows)


def check_environment(slim_path: str | Path, decoder_path: str | Path | None = None):
    """Return the environment record, raising if a hard requirement is unmet."""
    import subprocess

    from .run4_config import REQUIRED_SLIM_VERSION, REQUIRED_STDPOPSIM_VERSION

    if stdpopsim.__version__ != REQUIRED_STDPOPSIM_VERSION:
        raise RuntimeError(
            f"run4 requires stdpopsim {REQUIRED_STDPOPSIM_VERSION}; "
            f"found {stdpopsim.__version__}"
        )
    slim = Path(slim_path)
    if not slim.exists():
        raise RuntimeError(f"SLiM is unavailable at {slim}")
    reported = subprocess.check_output([str(slim), "-v"], text=True).strip()
    if REQUIRED_SLIM_VERSION not in reported:
        raise RuntimeError(
            f"run4 requires SLiM {REQUIRED_SLIM_VERSION}; {slim} reported: {reported}"
        )
    record = {
        "slim_path": str(slim),
        "slim_version_reported": reported,
        "cpu_count": os.cpu_count(),
        **_runtime_versions(),
    }
    if decoder_path is not None:
        decoder = Path(decoder_path)
        if not decoder.exists():
            raise RuntimeError(f"Gamma-SMC is unavailable at {decoder}")
        record["gamma_smc_path"] = str(decoder)
    return record


__all__ = [
    "GENOTYPE_CLASSES",
    "REPLICATE_OUTPUTS",
    "ReplicateTask",
    "build_tasks",
    "check_environment",
    "focal_allele_counts",
    "pairwise_tmrca_matrix",
    "profile_positions",
    "run_replicate",
    "run_tasks",
    "sampled_diploid_pairs",
]
