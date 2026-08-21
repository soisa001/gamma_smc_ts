#!/usr/bin/env python
"""Simulate run7: selection on an introgressed allele in the grafted EAS model.

    python scripts/run_run7.py --slim-bin .native-stdpopsim/bin/slim --workers 8
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "python"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import stdpopsim  # noqa: E402

from gamma_smc_aou.run2_models import scoped_focal_overlay_patch  # noqa: E402
from gamma_smc_aou.run5_simulate import (  # noqa: E402
    focal_allele_counts,
    pairwise_tmrca_matrix,
    profile_positions,
    sampled_diploid_pairs,
)
from gamma_smc_aou.run7_config import (  # noqa: E402
    EAS_INTROGRESSED,
    SELECTION_COEFFICIENTS,
    arm_id_for,
    EAS_POPULATION,
    FOCAL_POSITION_BP,
    MODES,
    N_REPLICATES,
    PROFILE_STRIDE_BP,
    SAMPLE_DIPLOIDS,
    SEQUENCE_LENGTH_BP,
    SLIM_BURN_IN,
    SLIM_SCALING_FACTOR,
    write_arm_configs,
)
from gamma_smc_aou.run7_models import (  # noqa: E402
    build_contig,
    build_extended_events,
    build_stdpopsim_model,
    scoped_slim_patch,
)

STUDY_ROOT = REPO / "sim_results_run7"


@dataclass(frozen=True)
class Task:
    mode: str
    replicate_index: int
    seed: int
    directory: str
    repo_root: str
    slim_path: str
    selection_coefficient: float
    arm_id: str


def run_one(task: Task) -> dict:
    directory = Path(task.directory)
    directory.mkdir(parents=True, exist_ok=True)
    started = perf_counter()
    arm = EAS_INTROGRESSED
    try:
        model = build_stdpopsim_model(task.repo_root)
        contig = build_contig()
        events = build_extended_events(task.mode, task.selection_coefficient)
        engine = stdpopsim.get_engine("slim")

        target_id = [p.name for p in model.model.populations].index(EAS_POPULATION)
        with (
            scoped_slim_patch(task.mode, target_id) as patch,
            scoped_focal_overlay_patch() as overlay,
        ):
            ts = engine.simulate(
                model,
                contig,
                {EAS_POPULATION: SAMPLE_DIPLOIDS},
                seed=task.seed,
                extended_events=list(events) or None,
                slim_path=task.slim_path,
                slim_scaling_factor=SLIM_SCALING_FACTOR,
                slim_burn_in=SLIM_BURN_IN,
                keep_mutation_ids_as_alleles=False,
            )
        if overlay["masked_calls"] < 1:
            raise RuntimeError("the reserved-base guard never fired")

        pairs = sampled_diploid_pairs(ts)
        counts = focal_allele_counts(ts, pairs)
        positions = profile_positions(
            SEQUENCE_LENGTH_BP, PROFILE_STRIDE_BP, FOCAL_POSITION_BP
        )
        tmrca = pairwise_tmrca_matrix(ts, pairs, positions)
        focal_index = int(np.searchsorted(positions, float(FOCAL_POSITION_BP)))

        np.savez_compressed(
            directory / "tmrca_matrix.npz",
            positions=positions.astype(np.float64),
            tmrca_generations=tmrca.astype(np.float32),
            pairs=pairs.astype(np.int32),
            focal_index=np.int64(focal_index),
        )
        np.save(directory / "focal_tmrca.npy", tmrca[focal_index])
        np.save(directory / "focal_genotypes.npy", counts)

        metadata = ts.metadata.get("SLiM", {}).get("user_metadata", {})

        def scalar(key):
            value = metadata.get(key)
            while isinstance(value, (list, tuple)) and len(value) == 1:
                value = value[0]
            return value

        endpoint = {
            "schema": "gamma-smc.run7-endpoint/v1",
            "unit_id": f"{task.arm_id}__{task.mode}__rep{task.replicate_index:03d}",
            "arm_id": task.arm_id,
            "selection_coefficient": task.selection_coefficient,
            "mode": task.mode,
            "replicate_index": task.replicate_index,
            "seed": task.seed,
            "status": "completed",
            "elapsed_seconds": perf_counter() - started,
            "focal_overlay_guard": dict(overlay),
            "slim_patch": patch,
            "final_allele_frequency": {
                "census_af": scalar("run7_final_census_af"),
                "census_state": scalar("run7_final_census_state"),
                "sample_af": float(counts.sum()) / (2 * len(pairs)),
                "sample_carrier_diploids": int((counts > 0).sum()),
                "sample_hom_carrier_diploids": int((counts == 2).sum()),
            },
            "placement": {
                "mode": scalar("run7_placement_mode"),
                "frequency": scalar("run7_placement_frequency"),
                "carrier_genomes": scalar("run7_placement_carrier_genomes"),
            },
        }
    except Exception as error:  # noqa: BLE001
        endpoint = {
            "schema": "gamma-smc.run7-endpoint/v1",
            "arm_id": task.arm_id,
            "selection_coefficient": task.selection_coefficient,
            "mode": task.mode,
            "replicate_index": task.replicate_index,
            "seed": task.seed,
            "status": "failed",
            "elapsed_seconds": perf_counter() - started,
            "error": f"{type(error).__name__}: {error}",
            "traceback": traceback.format_exc(),
        }

    path = directory / "endpoint.json"
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    temporary.write_text(
        json.dumps(endpoint, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)
    return {
        "arm_id": endpoint["arm_id"],
        "selection_coefficient": endpoint["selection_coefficient"],
        "mode": endpoint["mode"],
        "replicate_index": endpoint["replicate_index"],
        "status": endpoint["status"],
        "elapsed_seconds": endpoint["elapsed_seconds"],
        "error": endpoint.get("error", ""),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--slim-bin", required=True)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--replicates", type=int, default=N_REPLICATES)
    parser.add_argument("--repo-root", default=str(REPO))
    args = parser.parse_args(argv)

    write_arm_configs(STUDY_ROOT / "config")
    arm = EAS_INTROGRESSED
    tasks: list[Task] = []
    # The neutral mode does not depend on s, so it is simulated once and shared
    # by every selected arm rather than repeated four times over.
    for index in range(args.replicates):
        tasks.append(
            Task(
                mode="neutral",
                replicate_index=index,
                seed=arm.seed("neutral", index),
                directory=str(
                    STUDY_ROOT / "neutral" / "replicates" / f"rep{index:03d}"
                ),
                repo_root=args.repo_root,
                slim_path=args.slim_bin,
                selection_coefficient=0.0,
                arm_id="neutral",
            )
        )
    for s_index, coefficient in enumerate(SELECTION_COEFFICIENTS):
        arm_id = arm_id_for(coefficient)
        for index in range(args.replicates):
            tasks.append(
                Task(
                    mode="selected",
                    replicate_index=index,
                    seed=arm.seed("selected", index, s_index),
                    directory=str(
                        STUDY_ROOT / arm_id / "selected" / "replicates" / f"rep{index:03d}"
                    ),
                    repo_root=args.repo_root,
                    slim_path=args.slim_bin,
                    selection_coefficient=float(coefficient),
                    arm_id=arm_id,
                )
            )
    if args.workers == 1:
        rows = [run_one(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            rows = list(executor.map(run_one, tasks))
    status = pd.DataFrame(rows)
    (STUDY_ROOT / "logs").mkdir(parents=True, exist_ok=True)
    status.to_csv(STUDY_ROOT / "logs" / "simulation_status.tsv", sep="\t", index=False)
    failed = int((status["status"] != "completed").sum())
    print(json.dumps({"attempted": len(status), "failed": failed}, indent=2))
    if failed:
        print(status[status["status"] != "completed"].head(3).to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
