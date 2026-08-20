"""Replicate execution for the run2 study.

One replicate is one SLiM run plus the tree-truth statistics derived from it.
There is no ledger, no adapter and no recovery layer: a replicate directory is
self-describing, and a replicate is re-runnable in isolation because its seed is
a pure function of ``(arm, mode, index)``.

Nothing here conditions on anything.  A selected replicate that loses the focal
allele is recorded as lost and kept.
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

from .run2_config import (
    FOCAL_POSITION_BP,
    GENERATION_TIME_YEARS,
    MODES,
    PROFILE_STRIDE_BP,
    Run2Arm,
    SAMPLE_DIPLOIDS,
    SEQUENCE_LENGTH_BP,
    SLIM_BURN_IN,
    SLIM_SCALING_FACTOR,
    TMRCA_THRESHOLDS_YEARS,
    get_arm,
)
from .run2_models import (
    build_contig,
    build_demographic_model,
    build_extended_events,
    model_record,
    population_id,
    scoped_focal_overlay_patch,
    scoped_slim_patch,
)

GENOTYPE_CLASSES = ("overall", "hom_carrier", "het", "hom_noncarrier")
_CLASS_BY_COUNT = {0: "hom_noncarrier", 1: "het", 2: "hom_carrier"}

REPLICATE_OUTPUTS = {
    "endpoint": "endpoint.json",
    "focal_scores": "focal_scores.tsv",
    "spatial_profile": "spatial_profile.tsv.gz",
    "tree": "simulation.trees",
}


# ---------------------------------------------------------------------------
# Tree-sequence helpers
# ---------------------------------------------------------------------------


def sampled_diploid_pairs(ts: tskit.TreeSequence) -> np.ndarray:
    """Return the ``(n, 2)`` array of within-individual sample node pairs."""

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
) -> np.ndarray | None:
    """Return per-diploid focal allele counts, or ``None`` when the site is absent.

    The focal base is reserved from every neutral overlay, so a site here can
    only be the selected mutation.  Its absence means the allele was lost, which
    is a legitimate outcome that run2 records rather than rejects.
    """
    matches = [
        site
        for site in ts.sites()
        if np.isclose(float(site.position), float(focal_position), rtol=0.0, atol=1e-9)
    ]
    if not matches:
        return None
    if len(matches) > 1:
        raise ValueError(f"found {len(matches)} sites at the reserved focal base")
    ordered = [int(node) for pair in pairs for node in pair]
    variant = tskit.Variant(ts, samples=ordered)
    variant.decode(matches[0].id)
    genotypes = np.asarray(variant.genotypes, dtype=np.int64)
    if np.any(genotypes < 0):
        raise ValueError("missing focal genotypes cannot be classified")
    # Allele 0 is ancestral; anything else is the selected allele.
    return (genotypes.reshape(len(pairs), 2) != 0).sum(axis=1).astype(np.int8)


def profile_positions(
    sequence_length: int = SEQUENCE_LENGTH_BP,
    stride: int = PROFILE_STRIDE_BP,
    focal_position: int = FOCAL_POSITION_BP,
) -> np.ndarray:
    """Return the sorted stride grid, guaranteed to include the focal base."""
    grid = np.arange(0, sequence_length, stride, dtype=float)
    return np.unique(np.concatenate([grid, [float(focal_position)]]))


def pairwise_tmrca_matrix(
    ts: tskit.TreeSequence, pairs: np.ndarray, positions: np.ndarray
) -> np.ndarray:
    """Return the ``(n_positions, n_pairs)`` within-diploid TMRCA matrix.

    Trees are walked once in order rather than seeking per position, which is
    the difference between seconds and minutes at a 10 kb stride.
    """
    values = np.full((len(positions), len(pairs)), np.nan, dtype=float)
    index = 0
    for tree in ts.trees():
        _, right = tree.interval
        while index < len(positions) and positions[index] < right:
            values[index] = [
                tree.tmrca(int(left_node), int(right_node))
                for left_node, right_node in pairs
            ]
            index += 1
        if index >= len(positions):
            break
    if index != len(positions) or not np.all(np.isfinite(values)):
        raise ValueError("failed to evaluate a finite TMRCA at every profile position")
    return values


def _class_masks(counts: np.ndarray | None, n_pairs: int) -> dict[str, np.ndarray]:
    masks = {"overall": np.ones(n_pairs, dtype=bool)}
    for value, name in _CLASS_BY_COUNT.items():
        if counts is None:
            masks[name] = (
                np.ones(n_pairs, dtype=bool)
                if name == "hom_noncarrier"
                else np.zeros(n_pairs, dtype=bool)
            )
        else:
            masks[name] = counts == value
    return masks


def focal_scores(
    tmrca: np.ndarray,
    positions: np.ndarray,
    counts: np.ndarray | None,
    *,
    focal_position: int = FOCAL_POSITION_BP,
    thresholds_years: Sequence[float] = TMRCA_THRESHOLDS_YEARS,
    generation_time: float = GENERATION_TIME_YEARS,
) -> pd.DataFrame:
    """Return P(TMRCA < x) at the focal base, per threshold and genotype class.

    The ``overall`` class is the statistic the test uses; the carrier classes are
    descriptive only, and are undefined (NaN) whenever the class is empty --
    which is always the case for the neutral mode and for selected replicates
    that lost the allele.
    """
    focal_index = int(np.searchsorted(positions, float(focal_position)))
    if focal_index >= len(positions) or positions[focal_index] != float(focal_position):
        raise ValueError("the focal position is missing from the profile grid")
    times = tmrca[focal_index]
    masks = _class_masks(counts, len(times))
    rows = []
    for genotype_class in GENOTYPE_CLASSES:
        mask = masks[genotype_class]
        selected = times[mask]
        for threshold_years in thresholds_years:
            threshold_generations = float(threshold_years) / generation_time
            rows.append(
                {
                    "genotype_class": genotype_class,
                    "n_pairs": int(mask.sum()),
                    "threshold_years": float(threshold_years),
                    "threshold_generations": threshold_generations,
                    "p_tmrca_lt_threshold": (
                        float(np.mean(selected < threshold_generations))
                        if selected.size
                        else float("nan")
                    ),
                    "mean_tmrca_generations": (
                        float(selected.mean()) if selected.size else float("nan")
                    ),
                    "median_tmrca_generations": (
                        float(np.median(selected)) if selected.size else float("nan")
                    ),
                }
            )
    return pd.DataFrame.from_records(rows)


def spatial_profile(
    tmrca: np.ndarray,
    positions: np.ndarray,
    *,
    thresholds_years: Sequence[float] = TMRCA_THRESHOLDS_YEARS,
    generation_time: float = GENERATION_TIME_YEARS,
) -> pd.DataFrame:
    """Return the overall-class P(TMRCA < x) along the contig, wide by threshold."""

    frame = pd.DataFrame({"position_0based": positions})
    for threshold_years in thresholds_years:
        threshold_generations = float(threshold_years) / generation_time
        frame[f"p_lt_{int(threshold_years)}y"] = np.mean(
            tmrca < threshold_generations, axis=1
        )
    frame["mean_tmrca_generations"] = tmrca.mean(axis=1)
    frame["median_tmrca_generations"] = np.median(tmrca, axis=1)
    return frame


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
    save_tree: bool


def build_tasks(
    arm: Run2Arm,
    mode: str,
    *,
    study_root: str | Path,
    repo_root: str | Path,
    slim_path: str | Path,
    replicates: int,
    save_trees: bool = False,
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
            save_tree=bool(save_trees),
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
    """Read a SLiM user-metadata scalar, which arrives wrapped in a list."""
    if key not in block:
        return None
    value = block[key]
    while isinstance(value, (list, tuple)) and len(value) == 1:
        value = value[0]
    return value


def run_replicate(task: ReplicateTask) -> dict[str, Any]:
    """Run one replicate and write its outputs. Never raises; reports status."""

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
                extended_events=list(events) or None,
                slim_path=task.slim_path,
                slim_scaling_factor=SLIM_SCALING_FACTOR,
                slim_burn_in=SLIM_BURN_IN,
                keep_mutation_ids_as_alleles=False,
            )
        if overlay["masked_calls"] < 1:
            raise RuntimeError(
                "the focal base was never masked; the reserved-base guard did "
                "not fire in any neutral overlay"
            )

        pairs = sampled_diploid_pairs(ts)
        if len(pairs) != SAMPLE_DIPLOIDS:
            raise RuntimeError(
                f"expected {SAMPLE_DIPLOIDS} sampled diploids, got {len(pairs)}"
            )
        counts = (
            focal_allele_counts(ts, pairs) if task.mode == "selected" else None
        )
        positions = profile_positions()
        tmrca = pairwise_tmrca_matrix(ts, pairs, positions)

        scores = focal_scores(tmrca, positions, counts)
        profile = spatial_profile(tmrca, positions)
        scores.insert(0, "unit_id", task.unit_id)
        scores.insert(1, "replicate_index", task.replicate_index)
        scores.to_csv(directory / REPLICATE_OUTPUTS["focal_scores"], sep="\t", index=False)
        profile.to_csv(
            directory / REPLICATE_OUTPUTS["spatial_profile"], sep="\t", index=False
        )
        if task.save_tree:
            ts.dump(directory / REPLICATE_OUTPUTS["tree"])

        block = _metadata_block(ts)
        census_alt = _scalar(block, "run2_final_census_alt_count")
        census_total = _scalar(block, "run2_final_census_total_count")
        sample_alt = int(counts.sum()) if counts is not None else 0
        sample_total = 2 * len(pairs)

        endpoint: dict[str, Any] = {
            "schema": "gamma-smc.run2-endpoint/v1",
            "unit_id": task.unit_id,
            "arm_id": task.arm_id,
            "mode": task.mode,
            "replicate_index": task.replicate_index,
            "seed": task.seed,
            "status": "completed",
            "elapsed_seconds": perf_counter() - started,
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
            "standing_variation": {
                "tick": _scalar(block, "run2_standing_tick"),
                "population_id": _scalar(block, "run2_standing_population_id"),
                "placement": _scalar(block, "run2_standing_placement"),
                "carrier_genomes": _scalar(block, "run2_standing_carrier_genomes"),
                "population_genomes": _scalar(block, "run2_standing_population_genomes"),
                "population_individuals": _scalar(
                    block, "run2_standing_population_individuals"
                ),
                "realized_frequency": _scalar(block, "run2_standing_frequency"),
            },
            "final_allele_frequency": {
                "census_alt_count": census_alt,
                "census_total_count": census_total,
                "census_af": _scalar(block, "run2_final_census_af"),
                "census_state": _scalar(block, "run2_final_census_state"),
                "focal_site_present_in_sample": counts is not None,
                "sample_alt_count": sample_alt,
                "sample_total_count": sample_total,
                "sample_af": sample_alt / sample_total,
                "sample_carrier_diploids": (
                    int((counts > 0).sum()) if counts is not None else 0
                ),
            },
            "outputs": {
                key: str(directory / name)
                for key, name in REPLICATE_OUTPUTS.items()
                if key != "tree" or task.save_tree
            },
        }
    except Exception as error:  # noqa: BLE001 - a failed replicate must not stop the arm
        endpoint = {
            "schema": "gamma-smc.run2-endpoint/v1",
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
        "endpoint": str(path),
    }


def run_tasks(tasks: Sequence[ReplicateTask], *, workers: int = 1) -> pd.DataFrame:
    """Run replicates, optionally in parallel, and return the status table."""

    if workers < 1:
        raise ValueError("workers must be positive")
    if not tasks:
        return pd.DataFrame(
            columns=[
                "unit_id",
                "arm_id",
                "mode",
                "replicate_index",
                "seed",
                "status",
                "elapsed_seconds",
                "error",
                "endpoint",
            ]
        )
    if workers == 1:
        rows = [run_replicate(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            rows = list(executor.map(run_replicate, tasks))
    return pd.DataFrame(rows)


def check_environment(slim_path: str | Path) -> dict[str, Any]:
    """Return the G0 environment record, raising if a hard requirement is unmet."""

    import subprocess

    from .run2_config import REQUIRED_SLIM_VERSION, REQUIRED_STDPOPSIM_VERSION

    if stdpopsim.__version__ != REQUIRED_STDPOPSIM_VERSION:
        raise RuntimeError(
            f"run2 requires stdpopsim {REQUIRED_STDPOPSIM_VERSION}; "
            f"found {stdpopsim.__version__}"
        )
    path = Path(slim_path)
    if not path.exists():
        raise RuntimeError(f"SLiM is unavailable at {path}")
    reported = subprocess.check_output([str(path), "-v"], text=True).strip()
    if REQUIRED_SLIM_VERSION not in reported:
        raise RuntimeError(
            f"run2 requires SLiM {REQUIRED_SLIM_VERSION}; {path} reported: {reported}"
        )
    return {
        "slim_path": str(path),
        "slim_version_reported": reported,
        "cpu_count": os.cpu_count(),
        **_runtime_versions(),
    }


__all__ = [
    "GENOTYPE_CLASSES",
    "ReplicateTask",
    "build_tasks",
    "check_environment",
    "focal_allele_counts",
    "focal_scores",
    "pairwise_tmrca_matrix",
    "profile_positions",
    "run_replicate",
    "run_tasks",
    "sampled_diploid_pairs",
    "spatial_profile",
]
