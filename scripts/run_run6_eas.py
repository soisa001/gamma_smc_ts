#!/usr/bin/env python
"""Simulate run6's EAS arm: SLiM, full PHLASH history, marked selected allele.

Only the EAS arm is simulated here. The CHB side of run6 reuses run5's
``chb_from_introgression`` replicates, whose raw TMRCA matrices and focal
genotypes were retained and are exactly what the scan consumes.

    python scripts/run_run6_eas.py --slim-bin .native-stdpopsim/bin/slim --workers 8
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

from gamma_smc_aou.eas_sweep_models import (  # noqa: E402
    build_eas_demography_models,
    load_phlash_eas_npz,
)
from gamma_smc_aou.run2_models import (  # noqa: E402
    mask_focal_rate_map,
    scoped_focal_overlay_patch,
)
from gamma_smc_aou.run6_config import (  # noqa: E402
    DOMINANCE_COEFFICIENT,
    EAS_PLACEMENT_FREQUENCY,
    EAS_STANDING,
    FOCAL_POSITION_BP,
    FOCAL_SITE_ID,
    MODES,
    MUTATION_RATE,
    N_REPLICATES,
    PHLASH_EAS_RELATIVE_PATH,
    PHLASH_EAS_SHA256,
    PROFILE_STRIDE_BP,
    RECOMBINATION_RATE,
    SAMPLE_DIPLOIDS,
    SELECTION_COEFFICIENT,
    SEQUENCE_LENGTH_BP,
    SLIM_BURN_IN,
    SLIM_SCALING_FACTOR,
    write_arm_configs,
)
from gamma_smc_aou.run6_models import (  # noqa: E402
    build_contig,
    build_eas_model,
    build_extended_events,
    scoped_slim_patch,
)
# Pure tree-sequence helpers, reused rather than duplicated. They depend only on
# the tree sequence and explicit arguments, not on run5's configuration.
from gamma_smc_aou.run5_simulate import (  # noqa: E402
    focal_allele_counts,
    pairwise_tmrca_matrix,
    profile_positions,
    sampled_diploid_pairs,
)

STUDY_ROOT = REPO / "sim_results_run6"


@dataclass(frozen=True)
class Task:
    mode: str
    replicate_index: int
    seed: int
    directory: str
    repo_root: str
    slim_path: str


def run_one(task: Task) -> dict:
    directory = Path(task.directory)
    directory.mkdir(parents=True, exist_ok=True)
    started = perf_counter()
    try:
        arm = EAS_STANDING
        model = build_eas_model(task.repo_root)
        contig = build_contig()
        events = build_extended_events(task.mode)
        engine = stdpopsim.get_engine("slim")
        with (
            scoped_slim_patch(task.mode) as patch,
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
            "schema": "gamma-smc.run6-endpoint/v1",
            "unit_id": arm.unit_id(task.mode, task.replicate_index),
            "arm_id": arm.arm_id,
            "mode": task.mode,
            "replicate_index": task.replicate_index,
            "seed": task.seed,
            "status": "completed",
            "elapsed_seconds": perf_counter() - started,
            "generation_time_years": arm.generation_time,
            "focal_overlay_guard": dict(overlay),
            "slim_patch": patch,
            "final_allele_frequency": {
                "census_af": scalar("run6_final_census_af"),
                "census_state": scalar("run6_final_census_state"),
                "sample_af": float(counts.sum()) / (2 * len(pairs)),
                "sample_carrier_diploids": int((counts > 0).sum()),
                "sample_hom_carrier_diploids": int((counts == 2).sum()),
            },
            "placement": {
                "mode": scalar("run6_placement_mode"),
                "frequency": scalar("run6_placement_frequency"),
            },
        }
    except Exception as error:  # noqa: BLE001
        endpoint = {
            "schema": "gamma-smc.run6-endpoint/v1",
            "arm_id": "eas_standing",
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
    arm = EAS_STANDING
    tasks = [
        Task(
            mode=mode,
            replicate_index=index,
            seed=arm.seed(mode, index),
            directory=str(
                STUDY_ROOT / arm.arm_id / mode / "replicates" / f"rep{index:03d}"
            ),
            repo_root=args.repo_root,
            slim_path=args.slim_bin,
        )
        for mode in MODES
        for index in range(args.replicates)
    ]
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
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
