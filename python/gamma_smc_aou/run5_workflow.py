"""Phase orchestration for the run5 study.

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

from . import run5_analysis
from .run5_config import (
    ARMS,
    MODES,
    N_REPLICATES,
    Run5Arm,
    get_arm,
    shared_parameter_record,
    write_arm_configs,
)
from .run5_decode import build_decode_tasks, run_decode_tasks
from .run5_models import (
    generate_slim_script,
    model_record,
    tick_schedule,
)
from .run5_simulate import build_tasks, check_environment, run_tasks

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


def check_slim_script(arm: Run5Arm, mode: str, script: str) -> dict[str, Any]:
    """Assert the run5 invariants against a generated SLiM script."""
    schedule = tick_schedule(arm)
    offsets = schedule["tick_offsets_before_final_tick"]
    checks: dict[str, Any] = {"arm_id": arm.arm_id, "mode": mode}

    conditioning = _extract_constant(script, "condition_on_allele_frequency")
    conditioned = "c()" not in conditioning
    # The selected mode is ascertained on the allele actually having
    # introgressed to an appreciable frequency by the end of the migration
    # window. The neutral mode is unconditioned and carries no allele at all.
    if conditioned != (mode == "selected"):
        raise AssertionError(f"{mode} conditioning is wrong: conditioned={conditioned}")
    checks["ascertained_on_introgressed_frequency"] = conditioned

    drawn = _extract_constant(script, "drawn_mutations")
    fitness = _extract_constant(script, "fitness_callbacks")
    drawn_rows = 0 if "c()" in drawn else _matrix_row_count(drawn)
    fitness_rows = 0 if "c()" in fitness else _matrix_row_count(fitness)
    checks["drawn_mutation_rows"] = drawn_rows
    checks["fitness_callback_rows"] = fitness_rows

    expected_drawn = 1 if mode == "selected" else 0
    if drawn_rows != expected_drawn:
        raise AssertionError(
            f"{mode} expects {expected_drawn} drawn mutations, found {drawn_rows}"
        )

    if mode == "neutral":
        if "run5_placement_frequency" in script or "run5_final_census_af" in script:
            raise AssertionError("the neutral mode must not be patched")
        if fitness_rows:
            raise AssertionError("the neutral mode must have no fitness callbacks")
        checks["passed"] = True
        return checks

    if "run5_placement_frequency" not in script:
        raise AssertionError("the placement patch is not in the script")
    if "run5_final_census_af" not in script:
        raise AssertionError("the census patch is not in the script")

    if arm.nea_selection:
        if "sample(all_genomes, n_target)" not in script:
            raise AssertionError(
                "the Neanderthal-selection arm must start from archaic standing variation"
            )
    elif "carrier_genomes = pop.genomes;" not in script:
        raise AssertionError("the arm must fix the allele in the archaic source")

    expected_fitness = 2 if arm.nea_selection else 1
    if fitness_rows != expected_fitness:
        raise AssertionError(
            f"expected {expected_fitness} fitness callbacks, found {fitness_rows}"
        )
    if offsets["archaic_placement"] != offsets["archaic_split"] - 1:
        raise AssertionError("placement must be one tick after the archaic split")
    if offsets["chb_founding_and_selection_onset"] <= offsets["archaic_migration_end"]:
        raise AssertionError("selection must start before archaic migration ends")

    checks["selection_onset_offset"] = offsets["chb_founding_and_selection_onset"]
    checks["archaic_migration_end_offset"] = offsets["archaic_migration_end"]
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
    arms: Iterable[Run5Arm] = ARMS,
    slim_path: str | Path | None = None,
    decoder_path: str | Path | None = None,
) -> dict[str, Any]:
    root = study_root(repo_root)
    started = perf_counter()
    report: dict[str, Any] = {
        "schema": "gamma-smc.run5-validation/v1",
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
            entry: dict[str, Any] = {"model": model_record(arm, mode, repo_root)}
            if arm.engine == "slim":
                script = generate_slim_script(arm, mode, repo_root)
                (directory / f"{mode}.slim").write_text(script, encoding="utf-8")
                entry["script_checks"] = check_slim_script(arm, mode, script)
                entry["script_path"] = str(directory / f"{mode}.slim")
            else:
                entry["engine"] = "msprime"
                entry["script_checks"] = {"engine": "msprime", "passed": True}
            arm_report["modes"][mode] = entry
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
    arms: Iterable[Run5Arm] = ARMS,
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
    arms: Iterable[Run5Arm] = ARMS,
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


def phase_analyze(repo_root: str | Path, *, arms: Iterable[Run5Arm] = ARMS):
    root = study_root(repo_root)
    written: dict[str, Any] = {}
    for arm in arms:
        written[arm.arm_id] = {
            k: str(v) for k, v in run5_analysis.analyse_arm(root, arm).items()
        }
    cross = run5_analysis.cross_arm_power(root, ARMS)
    directory = root / "cross_arm" / "results"
    directory.mkdir(parents=True, exist_ok=True)
    cross.to_csv(directory / "power_by_arm.tsv", sep="\t", index=False)
    written["cross_arm"] = str(directory / "power_by_arm.tsv")
    return {"phase": "analyze", "written": written}


def phase_plot(repo_root: str | Path, *, arms: Iterable[Run5Arm] = ARMS):
    from . import run5_plots

    root = study_root(repo_root)
    written = {arm.arm_id: run5_plots.plot_arm(root, arm) for arm in arms}
    try:
        written["cross_arm"] = run5_plots.plot_cross_arm(root)
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
