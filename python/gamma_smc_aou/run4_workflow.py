"""Phase orchestration for the run4 study.

Phases: ``config``, ``validate``, ``simulate``, ``decode``, ``analyze``,
``plot``, ``all``.  ``validate`` needs stdpopsim but not SLiM; ``simulate``
needs SLiM; ``decode`` needs the Gamma-SMC binary and the tree sequences that
``simulate`` wrote.
"""

from __future__ import annotations

import json
import os
import re
from pathlib import Path
from time import perf_counter
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd

from . import run4_analysis
from .run4_config import (
    ARMS,
    MODES,
    N_REPLICATES,
    Run4Arm,
    get_arm,
    shared_parameter_record,
    write_arm_configs,
)
from .run4_decode import build_decode_tasks, run_decode_tasks
from .run4_models import (
    generate_slim_script,
    model_record,
    tick_schedule,
)
from .run4_simulate import build_tasks, check_environment, run_tasks

DEFAULT_STUDY_DIRNAME = "sim_results_run4"
PHASES = ("config", "validate", "simulate", "decode", "analyze", "plot", "all")


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
# validate
# ---------------------------------------------------------------------------


def _extract_constant(script: str, name: str) -> str:
    marker = f'defineConstant("{name}"'
    start = script.find(marker)
    if start < 0:
        raise AssertionError(f"generated SLiM script has no constant {name!r}")
    depth = 0
    for index in range(start, len(script)):
        if script[index] == "(":
            depth += 1
        elif script[index] == ")":
            depth -= 1
            if depth == 0:
                return script[start : index + 1]
    raise AssertionError(f"unterminated constant {name!r}")


def _matrix_row_count(block: str) -> int:
    match = re.search(r"c\(\s*(\d+)\s*,\s*(\d+)\s*\)\s*\)+\s*;?\s*$", block.strip())
    if not match:
        raise AssertionError(f"could not read the row count from {block.strip()[:160]}")
    return int(match.group(2))


def check_slim_script(arm: Run4Arm, mode: str, script: str) -> dict[str, Any]:
    """Assert the run4 invariants against a generated SLiM script."""
    schedule = tick_schedule(arm)
    checks: dict[str, Any] = {"arm_id": arm.arm_id, "mode": mode}

    conditioning = _extract_constant(script, "condition_on_allele_frequency")
    conditioned = "c()" not in conditioning
    # run4 conditions the SELECTED mode on the sweep surviving, because a single
    # beneficial copy is lost with probability ~1-2hs. The neutral mode is
    # unconditioned and carries no allele at all.
    if conditioned != (mode == "selected"):
        raise AssertionError(
            f"{mode} conditioning is wrong: conditioned={conditioned}"
        )
    checks["conditioned_on_survival"] = conditioned

    drawn = _extract_constant(script, "drawn_mutations")
    fitness = _extract_constant(script, "fitness_callbacks")
    drawn_rows = 0 if "c()" in drawn else _matrix_row_count(drawn)
    fitness_rows = 0 if "c()" in fitness else _matrix_row_count(fitness)

    # run4's neutral mode has NO allele at all.
    expected_drawn = 1 if mode == "selected" else 0
    expected_fitness = 1 if mode == "selected" else 0
    if drawn_rows != expected_drawn:
        raise AssertionError(
            f"{mode} expects {expected_drawn} drawn mutations, found {drawn_rows}"
        )
    if fitness_rows != expected_fitness:
        raise AssertionError(
            f"{mode} expects {expected_fitness} fitness callbacks, found {fitness_rows}"
        )
    offsets = schedule["tick_offsets_before_final_tick"]
    if mode == "neutral":
        if "run4_placement_frequency" in script or "run4_final_census_af" in script:
            raise AssertionError("the neutral mode must not be patched")
        checks["passed"] = True
        checks["drawn_mutation_rows"] = drawn_rows
        checks["fitness_callback_rows"] = fitness_rows
        return checks

    if "run4_placement_frequency" not in script:
        raise AssertionError("the placement patch is not in the script")
    if "run4_final_census_af" not in script:
        raise AssertionError("the census patch is not in the script")
    if arm.archaic:
        if "carrier_genomes = pop.genomes;" not in script:
            raise AssertionError("the archaic arm must fix the allele in the source")
        if "run4_single_founder" not in script:
            raise AssertionError("the archaic arm must reduce to a single founder")
        if offsets["archaic_fixation"] != offsets["archaic_split"] - 1:
            raise AssertionError("fixation must be one tick after the archaic split")
        if offsets["selection_onset"] >= offsets["archaic_migration_end"]:
            raise AssertionError("selection must start after archaic migration ends")
        checks["archaic_split_offset"] = offsets["archaic_split"]
    else:
        if "sample(pop.genomes, 1)" not in script:
            raise AssertionError("a de novo arm must place a single copy")
    checks["selection_onset_offset"] = offsets["selection_onset"]

    checks["drawn_mutation_rows"] = drawn_rows
    checks["fitness_callback_rows"] = fitness_rows
    checks["passed"] = True
    return checks


def check_seeds() -> dict[str, Any]:
    seeds = [
        arm.seed(mode, i)
        for arm in ARMS
        for mode in MODES
        for i in range(N_REPLICATES)
    ]
    if len(set(seeds)) != len(seeds):
        raise AssertionError("replicate seeds are not unique")
    return {"n_seeds": len(seeds), "unique": True, "passed": True}


def phase_config(repo_root: str | Path) -> dict[str, Any]:
    written = write_arm_configs(study_root(repo_root) / "config")
    return {"phase": "config", "written": {k: str(v) for k, v in written.items()}}


def phase_validate(
    repo_root: str | Path,
    *,
    arms: Iterable[Run4Arm] = ARMS,
    slim_path: str | Path | None = None,
    decoder_path: str | Path | None = None,
) -> dict[str, Any]:
    root = study_root(repo_root)
    started = perf_counter()
    report: dict[str, Any] = {
        "schema": "gamma-smc.run4-validation/v1",
        "shared_parameters": shared_parameter_record(),
        "seeds": check_seeds(),
        "arms": {},
    }
    for arm in arms:
        arm_report: dict[str, Any] = {
            "label": arm.label,
            "generation_time_years": arm.generation_time,
            "tick_schedule": tick_schedule(arm),
            "modes": {},
        }
        directory = root / arm.arm_id / "validation"
        directory.mkdir(parents=True, exist_ok=True)
        for mode in MODES:
            script = generate_slim_script(arm, mode, repo_root)
            (directory / f"{mode}.slim").write_text(script, encoding="utf-8")
            arm_report["modes"][mode] = {
                "model": model_record(arm, mode, repo_root),
                "script_checks": check_slim_script(arm, mode, script),
                "script_path": str(directory / f"{mode}.slim"),
            }
        report["arms"][arm.arm_id] = arm_report
        _atomic_json(directory / "validation_report.json", arm_report)

    if slim_path is not None:
        report["environment"] = check_environment(slim_path, decoder_path)
    else:
        report["environment"] = {"slim_checked": False}
    report["elapsed_seconds"] = perf_counter() - started
    report["passed"] = True
    _atomic_json(root / "validation_report.json", report)
    return report


# ---------------------------------------------------------------------------
# simulate / decode / analyze
# ---------------------------------------------------------------------------


def phase_simulate(
    repo_root: str | Path,
    *,
    slim_path: str | Path,
    arms: Iterable[Run4Arm] = ARMS,
    modes: Sequence[str] = MODES,
    replicates: int = N_REPLICATES,
    workers: int = 1,
    indices: Sequence[int] | None = None,
) -> dict[str, Any]:
    root = study_root(repo_root)
    environment = check_environment(slim_path)
    started = perf_counter()
    frames = []
    for arm in arms:
        for mode in modes:
            frames.append(
                run_tasks(
                    build_tasks(
                        arm, mode,
                        study_root=root, repo_root=repo_root, slim_path=slim_path,
                        replicates=replicates, indices=indices,
                    ),
                    workers=workers,
                )
            )
    status = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    (root / "logs").mkdir(parents=True, exist_ok=True)
    status.to_csv(root / "logs" / "simulation_status.tsv", sep="\t", index=False)
    summary = {
        "phase": "simulate",
        "environment": environment,
        "attempted": int(len(status)),
        "failed": int((status["status"] != "completed").sum()) if len(status) else 0,
        "elapsed_seconds": perf_counter() - started,
    }
    _atomic_json(root / "logs" / "simulation_summary.json", summary)
    return summary


def phase_decode(
    repo_root: str | Path,
    *,
    decoder_path: str | Path,
    arms: Iterable[Run4Arm] = ARMS,
    modes: Sequence[str] = MODES,
    replicates: int = N_REPLICATES,
    workers: int = 1,
    threads: int = 1,
    indices: Sequence[int] | None = None,
) -> dict[str, Any]:
    root = study_root(repo_root)
    if not Path(decoder_path).exists():
        raise RuntimeError(f"Gamma-SMC is unavailable at {decoder_path}")
    started = perf_counter()
    frames = []
    for arm in arms:
        for mode in modes:
            frames.append(
                run_decode_tasks(
                    build_decode_tasks(
                        arm, mode,
                        study_root=root, decoder_path=decoder_path,
                        replicates=replicates, threads=threads, indices=indices,
                    ),
                    workers=workers,
                )
            )
    status = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    (root / "logs").mkdir(parents=True, exist_ok=True)
    status.to_csv(root / "logs" / "decode_status.tsv", sep="\t", index=False)
    summary = {
        "phase": "decode",
        "decoder_path": str(decoder_path),
        "attempted": int(len(status)),
        "failed": int((status["status"] != "completed").sum()) if len(status) else 0,
        "elapsed_seconds": perf_counter() - started,
    }
    _atomic_json(root / "logs" / "decode_summary.json", summary)
    return summary


def phase_analyze(repo_root: str | Path, *, arms: Iterable[Run4Arm] = ARMS):
    root = study_root(repo_root)
    written: dict[str, Any] = {}
    for arm in arms:
        written[arm.arm_id] = {
            k: str(v) for k, v in run4_analysis.analyse_arm(root, arm).items()
        }
    cross = run4_analysis.cross_arm_power(root, ARMS)
    directory = root / "cross_arm" / "results"
    directory.mkdir(parents=True, exist_ok=True)
    cross.to_csv(directory / "power_by_arm.tsv", sep="\t", index=False)
    written["cross_arm"] = str(directory / "power_by_arm.tsv")
    return {"phase": "analyze", "written": written}


def phase_plot(repo_root: str | Path, *, arms: Iterable[Run4Arm] = ARMS):
    from . import run4_plots

    root = study_root(repo_root)
    written = {arm.arm_id: run4_plots.plot_arm(root, arm) for arm in arms}
    try:
        written["cross_arm"] = run4_plots.plot_cross_arm(root)
    except RuntimeError as error:
        written["cross_arm"] = {"skipped": str(error)}
    return {"phase": "plot", "written": written}


def run_phase(
    phase: str,
    repo_root: str | Path,
    *,
    slim_path: str | Path | None = None,
    decoder_path: str | Path | None = None,
    arm_ids: Sequence[str] | None = None,
    modes: Sequence[str] = MODES,
    replicates: int = N_REPLICATES,
    workers: int = 1,
    threads: int = 1,
    indices: Sequence[int] | None = None,
) -> dict[str, Any]:
    if phase not in PHASES:
        raise ValueError(f"phase must be one of {PHASES}, got {phase!r}")
    arms = ARMS if arm_ids is None else tuple(get_arm(a) for a in arm_ids)

    if phase == "config":
        return phase_config(repo_root)
    if phase == "validate":
        return phase_validate(
            repo_root, arms=arms, slim_path=slim_path, decoder_path=decoder_path
        )
    if phase == "simulate":
        if slim_path is None:
            raise ValueError("simulate requires --slim-bin")
        return phase_simulate(
            repo_root, slim_path=slim_path, arms=arms, modes=modes,
            replicates=replicates, workers=workers, indices=indices,
        )
    if phase == "decode":
        if decoder_path is None:
            raise ValueError("decode requires --decoder-bin")
        return phase_decode(
            repo_root, decoder_path=decoder_path, arms=arms, modes=modes,
            replicates=replicates, workers=workers, threads=threads, indices=indices,
        )
    if phase == "analyze":
        return phase_analyze(repo_root, arms=arms)
    if phase == "plot":
        return phase_plot(repo_root, arms=arms)

    if slim_path is None or decoder_path is None:
        raise ValueError("the all phase requires --slim-bin and --decoder-bin")
    steps = [
        phase_config(repo_root),
        phase_validate(repo_root, arms=arms, slim_path=slim_path, decoder_path=decoder_path),
        phase_simulate(
            repo_root, slim_path=slim_path, arms=arms, modes=modes,
            replicates=replicates, workers=workers, indices=indices,
        ),
        phase_decode(
            repo_root, decoder_path=decoder_path, arms=arms, modes=modes,
            replicates=replicates, workers=workers, threads=threads, indices=indices,
        ),
        phase_analyze(repo_root, arms=arms),
        phase_plot(repo_root, arms=arms),
    ]
    return {"phase": "all", "steps": steps}


__all__ = [
    "DEFAULT_STUDY_DIRNAME",
    "PHASES",
    "check_seeds",
    "check_slim_script",
    "phase_analyze",
    "phase_config",
    "phase_decode",
    "phase_plot",
    "phase_simulate",
    "phase_validate",
    "run_phase",
    "study_root",
]
