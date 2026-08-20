"""Bounded Han endpoint pilot; never used by the production campaign.

This pilot removes linked-sequence recording because endpoint calibration uses
only the terminal selected-locus frequency.  The single-locus Wright-Fisher
process, demography, pulse, selection, sample size, Q, and burn-in are unchanged.
"""

from __future__ import annotations

import argparse
import json
import math
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import replace
from pathlib import Path
from time import perf_counter

import pandas as pd
import stdpopsim
import numpy as np

import gamma_smc_aou.focused_selection_simulation as focused_simulation

from gamma_smc_aou.eas_sweep_analysis import build_diploid_pair_table
from gamma_smc_aou.eas_sweep_models import (
    FOCAL_SITE_ID,
    INTROGRESSION_PULSE_GENERATIONS,
    build_focal_contig,
)
from gamma_smc_aou.focused_selection_calibration import apply_han_selection_end
from gamma_smc_aou.focused_selection_simulation import (
    DEFAULT_SLIM_BURN_IN,
    _selected_focal_expectation,
    _selected_objects,
    _stable_seed,
    _wall_clock_timeout,
    candidate_pool_frequency_gate,
    feasible_exact_panel_genotype_triples,
    inspect_selected_focal_identity,
)


def _unit(s: float, target: float, seed: int, label: str) -> dict[str, object]:
    return {
        "unit_id": label,
        "demography_id": "ancient_eurasia_han_introgression",
        "demography_kind": "introgression",
        "population": "Han",
        "source_model": "stdpopsim AncientEurasia_9K19",
        "selection_origin": "archaic_specific_introgressed_standing_variation",
        "simulation_class": "selected",
        "selection_coefficient": s,
        "target_allele_frequency": target,
        "population_af_lower": target - 0.025,
        "population_af_upper": target + 0.025,
        "exact_sample_alt_count": round(200 * target),
        "sample_diploids": 100,
        "replicate_index": 1,
        "seed": seed,
        "sequence_length_bp": 10_000_000,
        "focal_position_bp": 5_000_000,
    }


def _add_successful_transfer_condition(sweep):
    post_pulse = stdpopsim.GenerationAfter(INTROGRESSION_PULSE_GENERATIONS)
    events = list(sweep.extended_events)
    events.append(
        stdpopsim.ConditionOnAlleleFrequency(
            start_time=post_pulse,
            end_time=post_pulse,
            single_site_id=FOCAL_SITE_ID,
            population="Loschbour",
            op=">",
            allele_frequency=0.0,
        )
    )
    return replace(sweep, extended_events=tuple(events))


def _add_recipient_nonloss_conditions(sweep, q: float):
    """Reject extinction along the unique post-pulse recipient lineage.

    With a nonrecurrent focal allele, terminal AF > 0 implies that every
    accepted path already satisfies these conditions.  They therefore alter
    rejection-sampling efficiency, not the target set of accepted paths.
    """

    realized_split = math.floor(2016.0 / q) * q
    events = list(sweep.extended_events)
    events.extend(
        (
            stdpopsim.ConditionOnAlleleFrequency(
                start_time=stdpopsim.GenerationAfter(
                    INTROGRESSION_PULSE_GENERATIONS
                ),
                end_time=realized_split,
                single_site_id=FOCAL_SITE_ID,
                population="Loschbour",
                op=">",
                allele_frequency=0.0,
            ),
            stdpopsim.ConditionOnAlleleFrequency(
                start_time=stdpopsim.GenerationAfter(realized_split),
                end_time=0.0,
                single_site_id=FOCAL_SITE_ID,
                population="Han",
                op=">",
                allele_frequency=0.0,
            ),
        )
    )
    return replace(sweep, extended_events=tuple(events))


def _draw(task: dict[str, object]) -> dict[str, object]:
    root = Path(str(task["repo_root"]))
    q = float(task["q"])
    s = float(task["s"])
    target = float(task["target"])
    multiplier = float(task["multiplier"])
    draw_index = int(task["draw_index"])
    mode = str(task["mode"])
    label = (
        f"pilot__{mode}__s{s:.3f}__af{target:.2f}__"
        f"m{multiplier:.3f}__draw{draw_index:03d}"
    )
    seed = int(task.get("override_seed", _stable_seed(20260813, label)))
    unit = _unit(s, target, seed, label)
    # A forked worker can be reused after a previous task changed this pilot-
    # local module global; restore the production coordinate before validation.
    focused_simulation.FOCAL_POSITION_BP = 5_000_000
    model, sweep, samples, _ = _selected_objects(unit, root, q, 500)
    analytic_end = 2272.0 - multiplier * (
        2272.0
        - (
            2272.0
            - (
                math.log(target / (1 - target))
                - math.log((0.0296 / 16) / (1 - 0.0296 / 16))
            )
            / (0.5 * s)
        )
    )
    endpoint = math.floor(analytic_end / q + 0.5) * q
    sweep, _ = apply_han_selection_end(sweep, endpoint)
    if mode == "successful_transfer":
        sweep = _add_successful_transfer_condition(sweep)
    elif mode == "recipient_nonloss":
        sweep = _add_recipient_nonloss_conditions(sweep, q)
    elif mode != "source_presence_only":
        raise ValueError(mode)

    # Terminal-AF calibration is exactly a one-locus task.  Keeping the same
    # 10 Mb map would only record linked neutral genealogies that are never read.
    sweep = replace(
        sweep,
        contig=build_focal_contig(sequence_length=2, focal_position=1),
    )
    # The production inspector intentionally binds to the production focal
    # coordinate.  Point its module-level coordinate at this pilot-only locus
    # inside the isolated worker process.
    focused_simulation.FOCAL_POSITION_BP = 1
    expectation = _selected_focal_expectation(model, sweep, q)
    started = perf_counter()
    try:
        with _wall_clock_timeout(float(task["timeout_seconds"])):
            ts = stdpopsim.get_engine("slim").simulate(
                model,
                sweep.contig,
                samples,
                seed=seed,
                extended_events=list(sweep.extended_events),
                slim_path=str(task["slim_path"]),
                slim_scaling_factor=q,
                slim_burn_in=DEFAULT_SLIM_BURN_IN,
                keep_mutation_ids_as_alleles=False,
            )
    except Exception as error:  # noqa: BLE001 - pilot records every failure
        return {
            **task,
            "label": label,
            "seed": seed,
            "endpoint": endpoint,
            "status": "failed",
            "terminal_class": "failed",
            "pool_af": None,
            "band_hit": False,
            "panel_feasible": False,
            "elapsed_seconds": perf_counter() - started,
            "error": f"{type(error).__name__}: {error}",
        }
    identity = inspect_selected_focal_identity(ts, expectation)
    if identity["status"] != "valid":
        terminal = "loss"
        af = 0.0
        hit = False
        feasible = False
    else:
        counts = build_diploid_pair_table(ts, 1)[
            "focal_selected_allele_count"
        ].to_numpy(dtype=int)
        alt_count = int(counts.sum())
        af = alt_count / 1000
        terminal = "fixed" if alt_count == 1000 else "segregating"
        hit = bool(
            target - 0.025 <= af <= target + 0.025
        )
        feasible = bool(
            feasible_exact_panel_genotype_triples(
                counts,
                sample_diploids=100,
                exact_alt_count=round(200 * target),
            )
        )
    return {
        **task,
        "label": label,
        "seed": seed,
        "endpoint": endpoint,
        "status": "complete",
        "terminal_class": terminal,
        "pool_af": af,
        "band_hit": hit and feasible,
        "panel_feasible": feasible,
        "identity": json.dumps(identity, sort_keys=True),
        "elapsed_seconds": perf_counter() - started,
        "error": "",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--slim-path", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--mode",
        choices=(
            "source_presence_only",
            "successful_transfer",
            "recipient_nonloss",
        ),
        default="successful_transfer",
    )
    parser.add_argument("--s", type=float, default=0.01)
    parser.add_argument("--target", type=float, default=0.10)
    parser.add_argument("--multipliers", type=float, nargs="+", default=(0.25, 0.4, 0.55, 0.7))
    parser.add_argument("--draws", type=int, default=8)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--q", type=float, default=5.0)
    parser.add_argument("--timeout-seconds", type=float, default=600.0)
    parser.add_argument("--seed-mode", choices=("by_cell", "common"), default="by_cell")
    args = parser.parse_args()
    tasks = [
        {
            "repo_root": str(args.repo_root.resolve()),
            "slim_path": str(args.slim_path.resolve()),
            "q": args.q,
            "mode": args.mode,
            "s": args.s,
            "target": args.target,
            "multiplier": multiplier,
            "draw_index": index,
            "timeout_seconds": args.timeout_seconds,
        }
        for multiplier in args.multipliers
        for index in range(args.draws)
    ]
    if args.seed_mode == "common":
        # Reuse the same stochastic seed at all endpoints within this one cell.
        # This common-random-number design is useful for estimating the endpoint
        # response with fewer pilot trajectories; it does not touch production.
        common_seeds = [
            _stable_seed(
                20260813,
                f"pilot__{args.mode}__s{args.s:.3f}__af{args.target:.2f}__"
                f"common__draw{index:03d}",
            )
            for index in range(args.draws)
        ]
        for task in tasks:
            task["override_seed"] = common_seeds[int(task["draw_index"])]
    rows = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(_draw, task) for task in tasks]
        for future in as_completed(futures):
            row = future.result()
            rows.append(row)
            print(json.dumps(row, sort_keys=True), flush=True)
    frame = pd.DataFrame(rows).sort_values(["multiplier", "draw_index"])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.output, sep="\t", index=False)
    print(
        frame.groupby("multiplier", as_index=False).agg(
            n=("status", "size"),
            n_failed=("status", lambda x: int((x != "complete").sum())),
            n_loss=("terminal_class", lambda x: int((x == "loss").sum())),
            n_fixed=("terminal_class", lambda x: int((x == "fixed").sum())),
            n_band_hits=("band_hit", "sum"),
            median_af=("pool_af", "median"),
            mean_seconds=("elapsed_seconds", "mean"),
        ).to_csv(sep="\t", index=False),
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
