#!/usr/bin/env python
"""Simulate run9 and decode each replicate with Gamma-SMC on the same pairs.

Order matters here. The tree sequence is written first and decoded *before* truth
is computed, because the decoder's random pair draw is the authority on which
pairs exist: its manifest is read back and tree truth is evaluated on exactly
those haplotype pairs. Truth and decode are then comparable pair-for-pair rather
than merely in aggregate, which is what makes "decode versus truth" a statement
about the decoder instead of about two different pair sets.

The tree is deleted once decoded -- nothing reads it twice, and 450 replicates of
a 10 Mb contig would otherwise be several gigabytes of dead intermediate.

    python scripts/run_run9.py --slim-bin .native-stdpopsim/bin/slim \\
        --decoder bin/gamma_smc --workers 20
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
import tskit  # noqa: E402

from gamma_smc_aou.decoder import run_within_decoder  # noqa: E402
from gamma_smc_aou.run2_models import scoped_focal_overlay_patch  # noqa: E402
from gamma_smc_aou.run5_simulate import (  # noqa: E402
    pairwise_tmrca_matrix,
    profile_positions,
)
from gamma_smc_aou.run7_config import (  # noqa: E402
    EAS_POPULATION,
    FOCAL_POSITION_BP,
    GENERATION_TIME_YEARS,
    PROFILE_STRIDE_BP,
    SEQUENCE_LENGTH_BP,
    SLIM_BURN_IN,
    SLIM_SCALING_FACTOR,
)
from gamma_smc_aou.run9_config import (  # noqa: E402
    DECODER_SETTINGS,
    N_NEUTRAL,
    SAMPLE_DIPLOIDS,
    SCENARIOS,
    THRESHOLDS_YEARS,
    get_scenario,
    write_config,
)
from gamma_smc_aou.run9_models import (  # noqa: E402
    build_contig,
    build_events,
    build_stdpopsim_model,
    scoped_slim_patch,
)

STUDY_ROOT = REPO / "sim_results_run9"


def haplotype_nodes(ts: tskit.TreeSequence) -> np.ndarray:
    """Sample nodes in VCF haplotype order: individual ``i`` occupies ``2i, 2i+1``.

    Gamma-SMC indexes haplotypes the same way, so this array is the dictionary
    that turns a manifest's integer pair into a pair of tree-sequence nodes.
    """
    sample_nodes = set(int(node) for node in ts.samples())
    ordered: list[int] = []
    for individual in ts.individuals():
        nodes = [int(n) for n in individual.nodes if int(n) in sample_nodes]
        if len(nodes) != 2:
            raise ValueError(
                f"individual {individual.id} has {len(nodes)} sample nodes; expected 2"
            )
        ordered.extend(nodes)
    if len(set(ordered)) != len(ordered) or set(ordered) != sample_nodes:
        raise ValueError("sample nodes do not partition cleanly into diploids")
    return np.asarray(ordered, dtype=np.int64)


def read_pairs_manifest(path: Path) -> np.ndarray:
    """Parse the decoder's manifest into an ``(n, 2)`` array of haplotype indices."""
    pairs: list[tuple[int, int]] = []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.replace(",", " ").split()
        if len(fields) < 2:
            continue
        try:
            left, right = int(fields[0]), int(fields[1])
        except ValueError:
            # A textual header row; the manifest is reusable as --pairs_file.
            continue
        pairs.append((left, right))
    if not pairs:
        raise ValueError(f"no haplotype pairs parsed from {path}")
    return np.asarray(pairs, dtype=np.int64)


def focal_haplotype_carriers(
    ts: tskit.TreeSequence, nodes: np.ndarray, focal_position: int = FOCAL_POSITION_BP
) -> np.ndarray:
    """Per-haplotype carrier status at the focal base, in VCF haplotype order.

    Returns all-False when the site is absent, which is the normal outcome for a
    neutral replicate and for any replicate where the allele was lost.
    """
    matches = [
        site
        for site in ts.sites()
        if np.isclose(float(site.position), float(focal_position), rtol=0.0, atol=1e-9)
    ]
    if not matches:
        return np.zeros(nodes.size, dtype=bool)
    if len(matches) > 1:
        raise ValueError(f"found {len(matches)} sites at the reserved focal base")
    variant = tskit.Variant(ts, samples=[int(n) for n in nodes])
    variant.decode(matches[0].id)
    genotypes = np.asarray(variant.genotypes, dtype=np.int64)
    if np.any(genotypes < 0):
        raise ValueError("missing focal genotypes cannot be classified")
    return genotypes != 0


@dataclass(frozen=True)
class Task:
    scenario_id: str
    mode: str
    replicate_index: int
    seed: int
    directory: str
    repo_root: str
    slim_path: str
    decoder_path: str


def run_one(task: Task) -> dict:
    directory = Path(task.directory)
    directory.mkdir(parents=True, exist_ok=True)
    started = perf_counter()
    tree_path = directory / "simulation.trees"
    scenario = get_scenario(task.scenario_id)
    try:
        model = build_stdpopsim_model(task.repo_root)
        contig = build_contig()
        events = build_events(scenario, task.mode)
        engine = stdpopsim.get_engine("slim")
        target_id = [p.name for p in model.model.populations].index(EAS_POPULATION)

        with (
            scoped_slim_patch(
                task.mode,
                target_id,
                patch_placement=scenario.places_allele_in_archaic,
            ) as patch,
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

        # ---- decode first: its manifest defines the pair set ----
        ts.dump(tree_path)
        summary_path = directory / "gamma_smc_summary.tsv"
        manifest_path = directory / "gamma_smc_pairs.tsv"
        decode_result = run_within_decoder(
            task.decoder_path,
            tree_path,
            summary_path,
            threshold_years=list(THRESHOLDS_YEARS),
            generation_time=GENERATION_TIME_YEARS,
            pairs_manifest=manifest_path,
            **DECODER_SETTINGS,
        )
        decoded = pd.read_csv(summary_path, sep="\t")
        pair_indices = read_pairs_manifest(manifest_path)

        # ---- truth on exactly those pairs ----
        nodes = haplotype_nodes(ts)
        if pair_indices.max() >= nodes.size:
            raise ValueError("manifest references a haplotype outside the panel")
        node_pairs = nodes[pair_indices]
        positions = profile_positions(
            SEQUENCE_LENGTH_BP, PROFILE_STRIDE_BP, FOCAL_POSITION_BP
        )
        tmrca = pairwise_tmrca_matrix(ts, node_pairs, positions)
        focal_index = int(np.searchsorted(positions, float(FOCAL_POSITION_BP)))
        carriers = focal_haplotype_carriers(ts, nodes)

        np.savez_compressed(
            directory / "tmrca_matrix.npz",
            positions=positions.astype(np.float64),
            tmrca_generations=tmrca.astype(np.float32),
            pair_haplotypes=pair_indices.astype(np.int32),
            focal_index=np.int64(focal_index),
        )
        np.save(directory / "focal_haplotype_carriers.npy", carriers)

        metadata = ts.metadata.get("SLiM", {}).get("user_metadata", {})

        def scalar(key):
            value = metadata.get(key)
            while isinstance(value, (list, tuple)) and len(value) == 1:
                value = value[0]
            return value

        # Carrier status per pair, which is the stratification on haplotype pairs.
        left = carriers[pair_indices[:, 0]]
        right = carriers[pair_indices[:, 1]]
        both = int(np.sum(left & right))
        one = int(np.sum(left ^ right))
        neither = int(np.sum(~left & ~right))

        endpoint = {
            "schema": "gamma-smc.run9-endpoint/v1",
            "scenario_id": task.scenario_id,
            "selection_coefficient": scenario.selection_coefficient,
            "origin": scenario.origin,
            "onset_generations": scenario.onset_generations,
            "mode": task.mode,
            "replicate_index": task.replicate_index,
            "seed": task.seed,
            "status": "completed",
            "elapsed_seconds": perf_counter() - started,
            "slim_patch": patch,
            "final_allele_frequency": {
                "census_af": scalar("run7_final_census_af"),
                "census_state": scalar("run7_final_census_state"),
                "sample_af": float(carriers.mean()),
                "carrier_haplotypes": int(carriers.sum()),
                "n_haplotypes": int(carriers.size),
            },
            "pair_classes": {
                "n_pairs": int(pair_indices.shape[0]),
                "carrier_carrier": both,
                "carrier_noncarrier": one,
                "noncarrier_noncarrier": neither,
            },
            "decode": {
                "rows": int(len(decoded)),
                "columns": list(decoded.columns),
                "theta": DECODER_SETTINGS["scaled_mutation_rate"],
                "rho": DECODER_SETTINGS["scaled_mutation_rate"]
                * DECODER_SETTINGS["recombination_to_mutation_ratio"],
                "recent_call": DECODER_SETTINGS["recent_call"],
                "result": decode_result,
            },
        }
    except Exception as error:  # noqa: BLE001
        endpoint = {
            "schema": "gamma-smc.run9-endpoint/v1",
            "scenario_id": task.scenario_id,
            "mode": task.mode,
            "replicate_index": task.replicate_index,
            "seed": task.seed,
            "status": "failed",
            "elapsed_seconds": perf_counter() - started,
            "error": f"{type(error).__name__}: {error}",
            "traceback": traceback.format_exc(),
        }
    finally:
        if tree_path.is_file():
            tree_path.unlink()

    path = directory / "endpoint.json"
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    temporary.write_text(
        json.dumps(endpoint, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)
    return {
        "scenario_id": task.scenario_id,
        "mode": task.mode,
        "replicate_index": task.replicate_index,
        "status": endpoint["status"],
        "elapsed_seconds": endpoint["elapsed_seconds"],
        "error": endpoint.get("error", ""),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--slim-bin", required=True)
    parser.add_argument("--decoder", required=True)
    parser.add_argument("--workers", type=int, default=20)
    parser.add_argument("--replicates", type=int, default=None)
    parser.add_argument("--neutral", type=int, default=None)
    parser.add_argument("--repo-root", default=str(REPO))
    parser.add_argument("--redo", action="store_true")
    args = parser.parse_args(argv)

    write_config(STUDY_ROOT / "config")
    tasks: list[Task] = []

    n_neutral = args.neutral if args.neutral is not None else N_NEUTRAL
    for index in range(n_neutral):
        tasks.append(
            Task(
                scenario_id=SCENARIOS[0].scenario_id,
                mode="neutral",
                replicate_index=index,
                seed=SCENARIOS[0].seed("neutral", index, 0),
                directory=str(STUDY_ROOT / "neutral" / "replicates" / f"rep{index:03d}"),
                repo_root=args.repo_root,
                slim_path=args.slim_bin,
                decoder_path=args.decoder,
            )
        )
    for scenario_index, scenario in enumerate(SCENARIOS):
        count = args.replicates if args.replicates is not None else scenario.n_replicates
        for index in range(count):
            tasks.append(
                Task(
                    scenario_id=scenario.scenario_id,
                    mode="selected",
                    replicate_index=index,
                    seed=scenario.seed("selected", index, scenario_index),
                    directory=str(
                        STUDY_ROOT
                        / scenario.scenario_id
                        / "selected"
                        / "replicates"
                        / f"rep{index:03d}"
                    ),
                    repo_root=args.repo_root,
                    slim_path=args.slim_bin,
                    decoder_path=args.decoder,
                )
            )

    if not args.redo:
        keep = []
        for task in tasks:
            endpoint = Path(task.directory) / "endpoint.json"
            if endpoint.is_file():
                try:
                    if json.loads(endpoint.read_text())["status"] == "completed":
                        continue
                except (ValueError, KeyError):
                    pass
            keep.append(task)
        print(json.dumps({"queued": len(keep), "skipped": len(tasks) - len(keep)}))
        tasks = keep
    if not tasks:
        print("nothing to do")
        return 0

    if args.workers == 1:
        rows = [run_one(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            rows = list(executor.map(run_one, tasks))
    status = pd.DataFrame(rows)
    (STUDY_ROOT / "logs").mkdir(parents=True, exist_ok=True)
    path = STUDY_ROOT / "logs" / "simulation_status.tsv"
    if path.is_file():
        status = pd.concat([pd.read_csv(path, sep="\t"), status], ignore_index=True)
    status.to_csv(path, sep="\t", index=False)
    failed = int((status["status"] != "completed").sum())
    print(json.dumps({"attempted": len(status), "failed": failed}, indent=2))
    if failed:
        print(status[status["status"] != "completed"].head(3).to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
