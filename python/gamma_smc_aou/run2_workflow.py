"""Phase orchestration for the run2 study.

Phases are ``config``, ``validate``, ``simulate``, ``analyze``, ``plot`` and
``all``.

``validate`` is the interesting one: it needs stdpopsim but **not** the SLiM
binary, so the entire event schedule can be generated and asserted before a
single simulation is launched.  It is what turns the framework's gates G1, G2
and G4 from prose into checks.
"""

from __future__ import annotations

import json
import os
import re
from pathlib import Path
from time import perf_counter
from typing import Any, Iterable, Sequence

import msprime
import numpy as np
import pandas as pd

from . import run2_analysis, run2_plots
from .run2_config import (
    ARMS,
    FOCAL_POSITION_BP,
    MODES,
    N_REPLICATES,
    Run2Arm,
    SEQUENCE_LENGTH_BP,
    STANDING_FREQUENCY,
    get_arm,
    shared_parameter_record,
    write_arm_configs,
)
from .run2_models import (
    generate_slim_script,
    mask_focal_rate_map,
    model_record,
    tick_schedule,
)
from .run2_simulate import build_tasks, check_environment, run_tasks

DEFAULT_STUDY_DIRNAME = "sim_results_run2"

PHASES = ("config", "validate", "simulate", "analyze", "plot", "all")


def study_root(repo_root: str | Path) -> Path:
    return Path(repo_root) / DEFAULT_STUDY_DIRNAME


def _atomic_json(path: Path, payload: Any) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)
    return path


# ---------------------------------------------------------------------------
# config
# ---------------------------------------------------------------------------


def phase_config(repo_root: str | Path) -> dict[str, Any]:
    root = study_root(repo_root)
    written = write_arm_configs(root / "config")
    return {"phase": "config", "written": {k: str(v) for k, v in written.items()}}


# ---------------------------------------------------------------------------
# validate
# ---------------------------------------------------------------------------

_CONSTANT_PATTERN = 'defineConstant("{name}"'


def _extract_constant(script: str, name: str) -> str:
    marker = _CONSTANT_PATTERN.format(name=name)
    start = script.find(marker)
    if start < 0:
        raise AssertionError(f"generated SLiM script has no constant {name!r}")
    depth = 0
    for index in range(start, len(script)):
        character = script[index]
        if character == "(":
            depth += 1
        elif character == ")":
            depth -= 1
            if depth == 0:
                return script[start : index + 1]
    raise AssertionError(f"unterminated constant {name!r} in the generated script")


def _matrix_row_count(block: str) -> int:
    """Return the declared row count of a stdpopsim ``array(c(...), c(r, n))``."""
    match = re.search(r"c\(\s*(\d+)\s*,\s*(\d+)\s*\)\s*\)+\s*;?\s*$", block.strip())
    if not match:
        raise AssertionError(f"could not read the row count from: {block.strip()[:200]}")
    return int(match.group(2))


def check_slim_script(arm: Run2Arm, mode: str, script: str) -> dict[str, Any]:
    """Assert the run2 invariants against a generated SLiM script."""

    schedule = tick_schedule(arm)
    checks: dict[str, Any] = {"arm_id": arm.arm_id, "mode": mode}

    conditioning = _extract_constant(script, "condition_on_allele_frequency")
    checks["conditioning_block"] = conditioning.strip()
    if "c()" not in conditioning:
        raise AssertionError(
            "run2 must not condition on allele frequency, but the generated "
            f"script for {arm.arm_id}/{mode} contains: {conditioning}"
        )
    checks["rejection_sampling"] = False

    drawn = _extract_constant(script, "drawn_mutations")
    fitness = _extract_constant(script, "fitness_callbacks")
    drawn_rows = 0 if "c()" in drawn else _matrix_row_count(drawn)
    fitness_rows = 0 if "c()" in fitness else _matrix_row_count(fitness)
    checks["drawn_mutation_rows"] = drawn_rows
    checks["fitness_callback_rows"] = fitness_rows

    if mode == "selected":
        expected_fitness = 2 if arm.introgression else 1
        if drawn_rows != 1:
            raise AssertionError(f"expected one drawn mutation, found {drawn_rows}")
        if fitness_rows != expected_fitness:
            raise AssertionError(
                f"expected {expected_fitness} fitness callbacks, found {fitness_rows}"
            )
        if "run2_standing_frequency" not in script:
            raise AssertionError("the standing-variation patch is not in the script")
        if "run2_final_census_af" not in script:
            raise AssertionError("the census end patch is not in the script")
        if arm.introgression and "inds.migrant" not in script:
            raise AssertionError(
                "the introgression arm must select carriers by migrant status"
            )
        if arm.introgression:
            # The allele is placed in the Han ancestor, so the frequency the Han
            # sweep starts from must be measured at the split rather than assumed.
            split_years = (
                schedule["realized_generations_ago"]["han_split"]
                * schedule["generation_time_years"]
            )
            expected_call = f"time_to_tick({split_years:.1f})"
            if "run2_record_entry" not in script:
                raise AssertionError("the entry recorder is not in the script")
            if script.count(expected_call) != 2:
                raise AssertionError(
                    f"the entry recorder is not registered at {expected_call}"
                )
            checks["entry_recorder_years_ago"] = split_years
        elif "run2_record_entry" in script:
            raise AssertionError(
                "the non-introgression arm has no founding event to record"
            )
        if not arm.introgression and f"{STANDING_FREQUENCY:.10g}" not in script:
            raise AssertionError(
                "the non-introgression arm must place the standing frequency directly"
            )
        # The focal draw must be scheduled one tick after the pulse.
        offsets = schedule["tick_offsets_before_final_tick"]
        draw_years = float(
            re.search(r"c\(([0-9.]+), \d+, \d+, \d+, \d+\)", drawn).group(1)
        )
        expected_years = (
            schedule["realized_generations_ago"]["standing_variation"]
            * schedule["generation_time_years"]
        )
        if not np.isclose(draw_years, expected_years):
            raise AssertionError(
                f"the focal draw is scheduled at {draw_years} years, expected "
                f"{expected_years}"
            )
        checks["draw_years_ago"] = draw_years
        checks["standing_variation_tick_offset"] = offsets["standing_variation"]
        checks["pulse_tick_offset"] = offsets["pulse"]
    else:
        if drawn_rows or fitness_rows:
            raise AssertionError(
                "the neutral mode must have no drawn mutations and no fitness events"
            )
        if "run2_standing_frequency" in script or "run2_final_census_af" in script:
            raise AssertionError("the neutral mode must not be patched")
        checks["standing_variation_tick_offset"] = None

    checks["passed"] = True
    return checks


def check_focal_mask(focal_position: int = FOCAL_POSITION_BP) -> dict[str, Any]:
    """Assert the reserved-base guard zeroes exactly one base and nothing else."""

    rate = 1.25e-8
    uniform = msprime.RateMap(
        position=np.array([0.0, float(SEQUENCE_LENGTH_BP)]), rate=np.array([rate])
    )
    masked, changed = mask_focal_rate_map(uniform, focal_position)
    if not changed:
        raise AssertionError("the focal mask did not fire on a uniform rate map")
    span = masked.get_cumulative_mass(float(SEQUENCE_LENGTH_BP))
    expected = rate * (SEQUENCE_LENGTH_BP - 1)
    if not np.isclose(span, expected, rtol=0, atol=1e-18):
        raise AssertionError(
            f"masked integrated rate {span} differs from the expected {expected}"
        )
    focal_mass = masked.get_cumulative_mass(
        float(focal_position) + 1.0
    ) - masked.get_cumulative_mass(float(focal_position))
    if focal_mass != 0.0:
        raise AssertionError(f"the focal base retains mutational mass {focal_mass}")
    # A second application must be a no-op.
    _, again = mask_focal_rate_map(masked, focal_position)
    if again:
        raise AssertionError("the focal mask is not idempotent")
    return {
        "integrated_rate_after_mask": float(span),
        "expected_integrated_rate": float(expected),
        "focal_base_mass": float(focal_mass),
        "idempotent": True,
        "passed": True,
    }


def check_seeds() -> dict[str, Any]:
    seeds = [
        arm.seed(mode, index)
        for arm in ARMS
        for mode in MODES
        for index in range(N_REPLICATES)
    ]
    if len(set(seeds)) != len(seeds):
        raise AssertionError("replicate seeds are not unique across arms and modes")
    return {"n_seeds": len(seeds), "unique": True, "passed": True}


def phase_validate(
    repo_root: str | Path,
    *,
    arms: Iterable[Run2Arm] = ARMS,
    slim_path: str | Path | None = None,
    write_scripts: bool = True,
) -> dict[str, Any]:
    """Run every pre-flight check that does not need the SLiM binary."""

    root = study_root(repo_root)
    started = perf_counter()
    report: dict[str, Any] = {
        "schema": "gamma-smc.run2-validation/v1",
        "shared_parameters": shared_parameter_record(),
        "seeds": check_seeds(),
        "focal_mask": check_focal_mask(),
        "arms": {},
    }

    for arm in arms:
        arm_report: dict[str, Any] = {
            "label": arm.label,
            "tick_schedule": tick_schedule(arm),
            "modes": {},
        }
        directory = root / arm.arm_id / "validation"
        for mode in MODES:
            script = generate_slim_script(arm, mode, repo_root)
            arm_report["modes"][mode] = {
                "model": model_record(arm, mode, repo_root),
                "script_checks": check_slim_script(arm, mode, script),
                "script_characters": len(script),
            }
            if write_scripts:
                directory.mkdir(parents=True, exist_ok=True)
                path = directory / f"{mode}.slim"
                path.write_text(script, encoding="utf-8")
                arm_report["modes"][mode]["script_path"] = str(path)
        report["arms"][arm.arm_id] = arm_report
        _atomic_json(directory / "validation_report.json", arm_report)

    if slim_path is not None:
        report["environment"] = check_environment(slim_path)
    else:
        report["environment"] = {
            "slim_checked": False,
            "note": "SLiM was not checked; pass --slim-bin to run gate G0",
        }

    report["elapsed_seconds"] = perf_counter() - started
    report["passed"] = True
    _atomic_json(root / "validation_report.json", report)
    return report


# ---------------------------------------------------------------------------
# simulate
# ---------------------------------------------------------------------------


def phase_simulate(
    repo_root: str | Path,
    *,
    slim_path: str | Path,
    arms: Iterable[Run2Arm] = ARMS,
    modes: Sequence[str] = MODES,
    replicates: int = N_REPLICATES,
    workers: int = 1,
    save_trees: bool = False,
    indices: Sequence[int] | None = None,
) -> dict[str, Any]:
    root = study_root(repo_root)
    environment = check_environment(slim_path)
    started = perf_counter()
    frames = []
    for arm in arms:
        for mode in modes:
            tasks = build_tasks(
                arm,
                mode,
                study_root=root,
                repo_root=repo_root,
                slim_path=slim_path,
                replicates=replicates,
                save_trees=save_trees,
                indices=indices,
            )
            frames.append(run_tasks(tasks, workers=workers))
    status = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    (root / "logs").mkdir(parents=True, exist_ok=True)
    status_path = root / "logs" / "simulation_status.tsv"
    status.to_csv(status_path, sep="\t", index=False)

    failures = (
        int((status["status"] != "completed").sum()) if not status.empty else 0
    )
    summary = {
        "phase": "simulate",
        "environment": environment,
        "replicates_attempted": int(len(status)),
        "replicates_failed": failures,
        "elapsed_seconds": perf_counter() - started,
        "status_table": str(status_path),
    }
    _atomic_json(root / "logs" / "simulation_summary.json", summary)
    return summary


# ---------------------------------------------------------------------------
# analyze and plot
# ---------------------------------------------------------------------------


def phase_analyze(
    repo_root: str | Path, *, arms: Iterable[Run2Arm] = ARMS
) -> dict[str, Any]:
    root = study_root(repo_root)
    written: dict[str, Any] = {}
    for arm in arms:
        written[arm.arm_id] = {
            key: str(path) for key, path in run2_analysis.analyse_arm(root, arm).items()
        }
    cross = run2_analysis.cross_arm_comparison(root, ARMS)
    cross_dir = root / "cross_arm" / "results"
    cross_dir.mkdir(parents=True, exist_ok=True)
    cross_path = cross_dir / "power_by_arm.tsv"
    cross.to_csv(cross_path, sep="\t", index=False)
    written["cross_arm"] = {"power_by_arm": str(cross_path)}
    return {"phase": "analyze", "written": written}


def phase_plot(repo_root: str | Path, *, arms: Iterable[Run2Arm] = ARMS) -> dict[str, Any]:
    root = study_root(repo_root)
    written: dict[str, Any] = {}
    for arm in arms:
        written[arm.arm_id] = run2_plots.plot_arm(root, arm)
    try:
        written["cross_arm"] = run2_plots.plot_cross_arm(root)
    except RuntimeError as error:
        written["cross_arm"] = {"skipped": str(error)}
    return {"phase": "plot", "written": written}


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def run_phase(
    phase: str,
    repo_root: str | Path,
    *,
    slim_path: str | Path | None = None,
    arm_ids: Sequence[str] | None = None,
    modes: Sequence[str] = MODES,
    replicates: int = N_REPLICATES,
    workers: int = 1,
    save_trees: bool = False,
    indices: Sequence[int] | None = None,
) -> dict[str, Any]:
    if phase not in PHASES:
        raise ValueError(f"phase must be one of {PHASES}, got {phase!r}")
    arms = ARMS if arm_ids is None else tuple(get_arm(arm_id) for arm_id in arm_ids)

    if phase == "config":
        return phase_config(repo_root)
    if phase == "validate":
        return phase_validate(repo_root, arms=arms, slim_path=slim_path)
    if phase == "simulate":
        if slim_path is None:
            raise ValueError("the simulate phase requires --slim-bin")
        return phase_simulate(
            repo_root,
            slim_path=slim_path,
            arms=arms,
            modes=modes,
            replicates=replicates,
            workers=workers,
            save_trees=save_trees,
            indices=indices,
        )
    if phase == "analyze":
        return phase_analyze(repo_root, arms=arms)
    if phase == "plot":
        return phase_plot(repo_root, arms=arms)

    results = [phase_config(repo_root), phase_validate(repo_root, arms=arms, slim_path=slim_path)]
    if slim_path is None:
        raise ValueError("the all phase requires --slim-bin")
    results.append(
        phase_simulate(
            repo_root,
            slim_path=slim_path,
            arms=arms,
            modes=modes,
            replicates=replicates,
            workers=workers,
            save_trees=save_trees,
            indices=indices,
        )
    )
    results.append(phase_analyze(repo_root, arms=arms))
    results.append(phase_plot(repo_root, arms=arms))
    return {"phase": "all", "steps": results}


__all__ = [
    "DEFAULT_STUDY_DIRNAME",
    "PHASES",
    "check_focal_mask",
    "check_seeds",
    "check_slim_script",
    "phase_analyze",
    "phase_config",
    "phase_plot",
    "phase_simulate",
    "phase_validate",
    "run_phase",
    "study_root",
]
