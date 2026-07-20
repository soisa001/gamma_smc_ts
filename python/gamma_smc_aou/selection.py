from __future__ import annotations

import json
import os
import shutil
import subprocess
import warnings
from importlib import resources
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import msprime
import numpy as np
import pandas as pd
import tskit

from .calibration import monte_carlo_pvalue


def slim_executable(explicit: str | Path | None = None) -> str | None:
    if explicit is not None:
        return str(explicit)
    return os.environ.get("SLIM_BIN") or shutil.which("slim")


def within_individual_tmrca_grid(
    ts, positions: np.ndarray, threshold_generations: float
) -> pd.DataFrame:
    sample_nodes = set(ts.samples())
    pairs = []
    for individual in ts.individuals():
        nodes = [node for node in individual.nodes if node in sample_nodes]
        if len(nodes) == 2:
            pairs.append(tuple(nodes))
    if not pairs:
        raise ValueError("tree sequence has no sampled diploid individuals")
    rows = []
    for position in positions:
        tree = ts.at(float(position))
        times = np.fromiter((tree.tmrca(a, b) for a, b in pairs), dtype=float)
        rows.append({
            "position_0based": float(position),
            "n_pairs": len(times),
            "mean_p_tmrca_lt_threshold": float(np.mean(times < threshold_generations)),
            "mean_tmrca_generations": float(np.mean(times)),
        })
    return pd.DataFrame(rows)


def run_slim_hard_sweep(
    output_path: str | Path,
    *,
    executable: str | Path | None = None,
    population_size: int = 200,
    sequence_length: int = 100_000,
    sweep_position: int | None = None,
    selection_coefficient: float = 0.5,
    recombination_rate: float = 1e-7,
    burnin: int | None = None,
    max_tick: int = 100_000,
    seed: int = 24681357,
):
    executable = slim_executable(executable)
    if executable is None:
        raise RuntimeError("SLiM executable not found; set SLIM_BIN or install SLiM")
    output_path = Path(output_path).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    sweep_position = sequence_length // 2 if sweep_position is None else sweep_position
    burnin = 5 * population_size if burnin is None else burnin
    script = resources.files("gamma_smc_aou").joinpath("slim/hard_sweep.slim")
    definitions = {
        "OUT": f'"{output_path.as_posix()}"',
        "POPULATION_SIZE": population_size,
        "LENGTH": sequence_length,
        "SWEEP_POSITION": sweep_position,
        "SELECTION_COEFFICIENT": selection_coefficient,
        "RECOMBINATION_RATE": recombination_rate,
        "BURNIN": burnin,
        "MAX_TICK": max_tick,
    }
    command = [str(executable), "-s", str(seed)]
    for key, value in definitions.items():
        command.extend(["-d", f"{key}={value}"])
    command.append(str(script))
    completed = subprocess.run(command, check=True, capture_output=True, text=True)
    if "SWEEP_FIXED" not in completed.stdout or not output_path.exists():
        raise RuntimeError(f"SLiM did not produce a fixed sweep:\n{completed.stdout}\n{completed.stderr}")
    ts = tskit.load(output_path)
    if any(tree.num_roots > 1 for tree in ts.trees()):
        try:
            import pyslim
        except ImportError as error:
            raise RuntimeError("Recapitating a SLiM tree sequence requires pyslim") from error
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=msprime.TimeUnitsMismatchWarning)
            ts = pyslim.recapitate(
                ts,
                ancestral_Ne=population_size,
                recombination_rate=recombination_rate,
                random_seed=seed + 1,
            )
        ts.dump(output_path)
    return ts, completed.stdout


def validate_slim_hard_sweep(
    output_dir: str | Path,
    *,
    executable: str | Path | None = None,
    population_size: int = 200,
    sequence_length: int = 100_000,
    selection_coefficient: float = 0.5,
    recombination_rate: float = 1e-7,
    threshold_generations: float = 150,
    neutral_replicates: int = 39,
    seed: int = 24681357,
) -> dict[str, float]:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    sweep_position = sequence_length // 2
    ts, slim_stdout = run_slim_hard_sweep(
        output_dir / "hard_sweep.trees",
        executable=executable,
        population_size=population_size,
        sequence_length=sequence_length,
        sweep_position=sweep_position,
        selection_coefficient=selection_coefficient,
        recombination_rate=recombination_rate,
        seed=seed,
    )
    positions = np.unique(np.append(np.linspace(0, sequence_length - 1, 201), sweep_position))
    selected = within_individual_tmrca_grid(ts, positions, threshold_generations)
    selected.to_csv(output_dir / "hard_sweep_truth.tsv", sep="\t", index=False)
    center_index = int(np.argmin(np.abs(selected["position_0based"] - sweep_position)))
    center = selected.iloc[center_index]
    flank = selected[
        (selected["position_0based"] <= sequence_length * 0.2)
        | (selected["position_0based"] >= sequence_length * 0.8)
    ]

    neutral_stats = []
    for replicate in range(neutral_replicates):
        neutral = msprime.sim_ancestry(
            samples=[msprime.SampleSet(population_size, ploidy=2)],
            population_size=population_size,
            sequence_length=sequence_length,
            recombination_rate=recombination_rate,
            model=msprime.StandardCoalescent(),
            random_seed=seed + 10 + replicate,
        )
        neutral_row = within_individual_tmrca_grid(
            neutral, np.asarray([sweep_position]), threshold_generations
        ).iloc[0]
        neutral_stats.append(float(neutral_row["mean_p_tmrca_lt_threshold"]))
    pvalue = monte_carlo_pvalue(
        float(center["mean_p_tmrca_lt_threshold"]), neutral_stats
    )
    metrics = {
        "population_size": population_size,
        "selection_coefficient": selection_coefficient,
        "sweep_position": sweep_position,
        "center_mean_tmrca_generations": float(center["mean_tmrca_generations"]),
        "flank_mean_tmrca_generations": float(flank["mean_tmrca_generations"].mean()),
        "center_to_flank_tmrca_ratio": float(
            center["mean_tmrca_generations"] / flank["mean_tmrca_generations"].mean()
        ),
        "center_fraction_recent": float(center["mean_p_tmrca_lt_threshold"]),
        "flank_fraction_recent": float(flank["mean_p_tmrca_lt_threshold"].mean()),
        "neutral_mean_fraction_recent": float(np.mean(neutral_stats)),
        "neutral_replicates": neutral_replicates,
        "mc_p_upper": float(pvalue),
    }
    with (output_dir / "hard_sweep_metrics.json").open("w", encoding="utf-8") as handle:
        json.dump({**metrics, "slim_stdout": slim_stdout}, handle, indent=2)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2), constrained_layout=True)
    axes[0].plot(selected["position_0based"], selected["mean_tmrca_generations"], lw=1.0)
    axes[0].axvline(sweep_position, color="firebrick", ls="--")
    axes[0].set(xlabel="Position (bp)", ylabel="Mean within-pair TMRCA", title="SLiM hard-sweep genealogy")
    axes[1].plot(selected["position_0based"], selected["mean_p_tmrca_lt_threshold"], lw=1.0)
    axes[1].axvline(sweep_position, color="firebrick", ls="--")
    axes[1].set(xlabel="Position (bp)", ylabel="Fraction recent", title="Recent-coalescence signal")
    axes[2].hist(neutral_stats, bins=15, alpha=0.8)
    axes[2].axvline(center["mean_p_tmrca_lt_threshold"], color="firebrick", label="selected center")
    axes[2].set(xlabel="Neutral fraction recent", ylabel="Replicates", title=f"Monte Carlo p={pvalue:.4g}")
    axes[2].legend()
    fig.savefig(output_dir / "hard_sweep_validation.png", dpi=180)
    plt.close(fig)
    return metrics
