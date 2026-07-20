from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import warnings
from concurrent.futures import ThreadPoolExecutor
from importlib import resources
from pathlib import Path
from time import perf_counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import msprime
import numpy as np
import pandas as pd
import tskit

from .calibration import calibration_metrics, monte_carlo_pvalue


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


def _define_path(path: Path) -> str:
    return f'"{path.resolve().as_posix()}"'


def _sample_diploids(ts, n_diploids: int, seed: int):
    import pyslim

    alive = np.asarray(pyslim.individuals_alive_at(ts, 0), dtype=int)
    if n_diploids > len(alive):
        raise ValueError(
            f"requested {n_diploids} diploids but only {len(alive)} are alive"
        )
    chosen = np.random.default_rng(seed).choice(alive, size=n_diploids, replace=False)
    nodes = np.concatenate([ts.individual(int(i)).nodes for i in chosen])
    return ts.simplify(nodes, filter_individuals=True)


def run_slim_recent_sweep(
    output_path: str | Path,
    *,
    executable: str | Path | None = None,
    population_size: int = 10_000,
    sample_diploids: int = 2_000,
    sequence_length: int = 10_000_000,
    sweep_position: int | None = None,
    selection_coefficient: float = 0.01,
    age_generations: int = 180,
    mutation_rate: float = 1.25e-8,
    recombination_rate: float = 1e-8,
    seed: int = 24681357,
):
    """Simulate one unconditional recent single-origin selected trajectory.

    The neutral burn-in uses msprime's StandardCoalescent. SLiM then runs exactly
    ``age_generations`` Wright-Fisher generations. The focal mutation may be lost,
    segregating, or fixed at the present. The coalescence scan is evaluated at its
    coordinate in every case. Neutral mutations are overlaid after sampling.
    """
    if selection_coefficient < 0:
        raise ValueError("selection_coefficient must be nonnegative")
    if age_generations < 1 or population_size < 2 or sample_diploids < 1:
        raise ValueError("age and population/sample sizes must be positive")
    executable = slim_executable(executable)
    if executable is None:
        raise RuntimeError("SLiM executable not found; set SLIM_BIN or install SLiM")
    try:
        import pyslim
    except ImportError as error:
        raise RuntimeError("recent sweep simulation requires pyslim") from error

    output_path = Path(output_path).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    sweep_position = sequence_length // 2 if sweep_position is None else sweep_position
    initial_path = output_path.with_suffix(".initial.trees")
    raw_path = output_path.with_suffix(".raw.trees")

    started = perf_counter()
    initial = msprime.sim_ancestry(
        samples=[msprime.SampleSet(population_size, ploidy=2)],
        population_size=population_size,
        sequence_length=sequence_length,
        recombination_rate=recombination_rate,
        model=msprime.StandardCoalescent(),
        random_seed=seed,
    )
    pyslim.annotate(initial, model_type="WF", tick=1, stage="early").dump(initial_path)
    ancestry_seconds = perf_counter() - started

    script = resources.files("gamma_smc_aou").joinpath("slim/recent_sweep.slim")
    definitions = {
        "IN": _define_path(initial_path),
        "OUT": _define_path(raw_path),
        "LENGTH": sequence_length,
        "SWEEP_POSITION": sweep_position,
        "SELECTION_COEFFICIENT": selection_coefficient,
        "RECOMBINATION_RATE": recombination_rate,
        "AGE_GENERATIONS": age_generations,
    }
    command = [str(executable), "-s", str(seed + 1)]
    for key, value in definitions.items():
        command.extend(["-d", f"{key}={value}"])
    command.append(str(script))
    slim_started = perf_counter()
    completed = subprocess.run(command, check=True, capture_output=True, text=True)
    slim_seconds = perf_counter() - slim_started
    match = re.search(
        r"SWEEP_COMPLETE frequency=([0-9.eE+-]+) outcome=([a-z]+)",
        completed.stdout,
    )
    if match is None or not raw_path.exists():
        raise RuntimeError(
            f"SLiM did not produce a recent-sweep trajectory:\n"
            f"{completed.stdout}\n{completed.stderr}"
        )
    allele_frequency = float(match.group(1))
    focal_allele_outcome = match.group(2)
    sampled = _sample_diploids(tskit.load(raw_path), sample_diploids, seed + 2)
    tables = sampled.dump_tables()
    tables.sites.clear()
    tables.mutations.clear()
    sampled = tables.tree_sequence()
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=msprime.TimeUnitsMismatchWarning)
        sampled = msprime.sim_mutations(
            sampled,
            rate=mutation_rate,
            model=msprime.BinaryMutationModel(),
            random_seed=seed + 3,
        )
    sampled.dump(output_path)
    for path in (initial_path, raw_path):
        path.unlink(missing_ok=True)
    return sampled, {
        "selection_coefficient": selection_coefficient,
        "age_generations": age_generations,
        "population_size": population_size,
        "sample_diploids": sample_diploids,
        "realized_population_allele_frequency": allele_frequency,
        "focal_allele_outcome": focal_allele_outcome,
        "ancestry_seconds": ancestry_seconds,
        "slim_seconds": slim_seconds,
        "n_sites": sampled.num_sites,
        "n_trees": sampled.num_trees,
        "slim_stdout": completed.stdout,
    }


def _sweep_grid(sequence_length: int, sweep_position: int) -> np.ndarray:
    coarse = np.arange(0, sequence_length, 50_000, dtype=float)
    fine = np.arange(sweep_position - 50_000, sweep_position + 50_001, 1_000, dtype=float)
    return np.unique(np.clip(np.append(coarse, fine), 0, sequence_length - 1))


def _profile_statistics(profile: pd.DataFrame, center: int) -> dict[str, float]:
    distance = np.abs(profile["position_0based"].to_numpy(dtype=float) - center)
    exact = profile.iloc[int(np.argmin(distance))]
    result = {
        "center_fraction_recent": float(exact["mean_p_tmrca_lt_threshold"]),
        "center_mean_tmrca_generations": float(exact["mean_tmrca_generations"]),
    }
    for width in (10_000, 100_000):
        subset = profile[distance <= width / 2]
        result[f"window_{width}_fraction_recent"] = float(
            subset["mean_p_tmrca_lt_threshold"].mean()
        )
    flank = profile[distance >= 2_000_000]
    result["neutral_plateau_fraction_recent"] = float(
        flank["mean_p_tmrca_lt_threshold"].mean()
    )
    result["neutral_plateau_mean_tmrca_generations"] = float(
        flank["mean_tmrca_generations"].mean()
    )
    return result


def summarize_recent_sweep_results(
    output_dir: str | Path,
    *,
    population_size: int,
    age_generations: int,
) -> dict:
    """Write null-calibration, error, and selected-power summaries."""
    from scipy.stats import mannwhitneyu, spearmanr

    output_dir = Path(output_dir)
    stats = pd.read_csv(output_dir / "replicate_statistics.tsv", sep="\t")
    null = stats[stats["selection_coefficient"] == 0].reset_index(drop=True)
    if len(null) < 2:
        raise ValueError("at least two neutral replicates are required for calibration")
    stat_columns = (
        "center_fraction_recent",
        "window_10000_fraction_recent",
        "window_100000_fraction_recent",
    )
    loo = {}
    for column in stat_columns:
        values = null[column].to_numpy(dtype=float)
        pvalues = np.asarray([
            monte_carlo_pvalue(value, np.delete(values, index))
            for index, value in enumerate(values)
        ])
        loo[column] = {
            **calibration_metrics(pvalues),
            "fraction_p_le_0_05": float(np.mean(pvalues <= 0.05)),
        }

    expected_recent = float(1 - np.exp(-age_generations / (2 * population_size)))
    recent_errors = null["center_fraction_recent"].to_numpy(dtype=float) - expected_recent
    tmrca_errors = (
        null["center_mean_tmrca_generations"].to_numpy(dtype=float)
        - 2 * population_size
    )
    null_error = {
        "theoretical_center_fraction_recent": expected_recent,
        "observed_center_fraction_recent_mean": float(
            null["center_fraction_recent"].mean()
        ),
        "center_fraction_recent_bias": float(np.mean(recent_errors)),
        "center_fraction_recent_rmse": float(np.sqrt(np.mean(recent_errors**2))),
        "theoretical_mean_tmrca_generations": float(2 * population_size),
        "observed_center_mean_tmrca_generations": float(
            null["center_mean_tmrca_generations"].mean()
        ),
        "center_mean_tmrca_bias_generations": float(np.mean(tmrca_errors)),
        "center_mean_tmrca_rmse_generations": float(
            np.sqrt(np.mean(tmrca_errors**2))
        ),
    }

    power_rows = []
    selected = stats[stats["selection_coefficient"] > 0]
    for coefficient, group in selected.groupby("selection_coefficient"):
        frequency = group["realized_population_allele_frequency"].to_numpy(dtype=float)
        recent = group["center_fraction_recent"].to_numpy(dtype=float)
        if np.unique(frequency).size > 1 and np.unique(recent).size > 1:
            rho, rho_p = spearmanr(frequency, recent)
        else:
            rho, rho_p = np.nan, np.nan
        outcomes = group["focal_allele_outcome"].value_counts(normalize=True)
        base = {
            "selection_coefficient": float(coefficient),
            "n_replicates": int(len(group)),
            "median_realized_allele_frequency": float(np.median(frequency)),
            "range_realized_allele_frequency_min": float(np.min(frequency)),
            "range_realized_allele_frequency_max": float(np.max(frequency)),
            "fraction_lost": float(outcomes.get("lost", 0.0)),
            "fraction_segregating": float(outcomes.get("segregating", 0.0)),
            "fraction_fixed": float(outcomes.get("fixed", 0.0)),
            "spearman_frequency_vs_center_recent": float(rho),
            "spearman_p": float(rho_p),
        }
        for column in stat_columns:
            p_column = f"{column}_mc_p_upper"
            u = mannwhitneyu(
                group[column], null[column], alternative="greater", method="auto"
            ).statistic
            base[f"{column}_mean"] = float(group[column].mean())
            base[f"{column}_median_p"] = float(group[p_column].median())
            base[f"{column}_power_p_le_0_05"] = float(
                np.mean(group[p_column] <= 0.05)
            )
            base[f"{column}_auc"] = float(u / (len(group) * len(null)))
        power_rows.append(base)
    power = pd.DataFrame(power_rows)
    power.to_csv(output_dir / "power_summary.tsv", sep="\t", index=False)
    result = {
        "null_leave_one_out_calibration": loo,
        "null_truth_error": null_error,
        "selected_power": power.to_dict("records"),
    }
    with (output_dir / "statistical_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2)
    return result


def validate_recent_sweep_grid(
    output_dir: str | Path,
    *,
    executable: str | Path | None = None,
    population_size: int = 10_000,
    sample_diploids: int = 2_000,
    sequence_length: int = 10_000_000,
    selection_coefficients: tuple[float, ...] = (0.0, 0.001, 0.01),
    age_generations: int = 180,
    mutation_rate: float = 1.25e-8,
    recombination_rate: float = 1e-8,
    neutral_replicates: int = 100,
    selected_replicates: int = 1,
    workers: int = 1,
    save_trees: bool = True,
    seed: int = 271828,
) -> dict:
    """Validate the central recent-coalescence test on a full 10-Mb region."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    tree_dir = output_dir / "trees"
    tree_dir.mkdir(exist_ok=True)
    center = sequence_length // 2
    positions = _sweep_grid(sequence_length, center)
    threshold_generations = float(age_generations)
    profiles = []
    rows = []
    started = perf_counter()
    if workers < 1:
        raise ValueError("workers must be positive")

    def _trajectory_task(task):
        coefficient, replicate = task
        run_seed = seed + 1_000_000 + int(coefficient * 1e8) + replicate * 10
        destination = tree_dir / f"s_{coefficient:g}_rep_{replicate:03d}.trees"
        ts, run = run_slim_recent_sweep(
            destination,
            executable=executable,
            population_size=population_size,
            sample_diploids=sample_diploids,
            sequence_length=sequence_length,
            sweep_position=center,
            selection_coefficient=coefficient,
            age_generations=age_generations,
            mutation_rate=mutation_rate,
            recombination_rate=recombination_rate,
            seed=run_seed,
        )
        profile = within_individual_tmrca_grid(ts, positions, threshold_generations)
        profile["selection_coefficient"] = coefficient
        profile["replicate"] = replicate
        row = {
            "selection_coefficient": coefficient,
            "replicate": replicate,
            "realized_population_allele_frequency": run[
                "realized_population_allele_frequency"
            ],
            "focal_allele_outcome": run["focal_allele_outcome"],
            "ancestry_seconds": run["ancestry_seconds"],
            "slim_seconds": run["slim_seconds"],
            **_profile_statistics(profile, center),
        }
        if not save_trees:
            destination.unlink(missing_ok=True)
        return profile, row

    trajectory_tasks = [
        (float(coefficient), replicate)
        for coefficient in selection_coefficients
        for replicate in range(
            neutral_replicates if coefficient == 0 else selected_replicates
        )
    ]
    if workers == 1:
        trajectory_results = map(_trajectory_task, trajectory_tasks)
        for profile, row in trajectory_results:
            profiles.append(profile)
            rows.append(row)
    else:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            for profile, row in executor.map(_trajectory_task, trajectory_tasks):
                profiles.append(profile)
                rows.append(row)

    profile_table = pd.concat(profiles, ignore_index=True)
    stats = pd.DataFrame(rows)
    null = stats[stats["selection_coefficient"] == 0]
    for column in (
        "center_fraction_recent",
        "window_10000_fraction_recent",
        "window_100000_fraction_recent",
    ):
        stats[f"{column}_mc_p_upper"] = [
            np.nan if coefficient == 0 else monte_carlo_pvalue(value, null[column])
            for coefficient, value in zip(stats["selection_coefficient"], stats[column])
        ]
    stats.to_csv(
        output_dir / "replicate_statistics.tsv", sep="\t", index=False, na_rep="NA"
    )
    profile_table.to_csv(output_dir / "truth_profiles.tsv", sep="\t", index=False)

    summary = (
        profile_table.groupby(["selection_coefficient", "position_0based"], as_index=False)
        .agg(
            mean_tmrca_generations=("mean_tmrca_generations", "mean"),
            q25_tmrca=("mean_tmrca_generations", lambda x: x.quantile(0.25)),
            q75_tmrca=("mean_tmrca_generations", lambda x: x.quantile(0.75)),
            mean_fraction_recent=("mean_p_tmrca_lt_threshold", "mean"),
            q25_fraction_recent=("mean_p_tmrca_lt_threshold", lambda x: x.quantile(0.25)),
            q75_fraction_recent=("mean_p_tmrca_lt_threshold", lambda x: x.quantile(0.75)),
        )
    )
    summary.to_csv(output_dir / "profile_summary.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(2, 2, figsize=(14, 8), constrained_layout=True)
    colors = {0.0: "0.35", 0.001: "#3b82f6", 0.01: "#dc2626"}
    for coefficient, group in summary.groupby("selection_coefficient"):
        label = f"s={coefficient:g}"
        color = colors.get(float(coefficient))
        x = group["position_0based"].to_numpy() / 1e6
        axes[0, 0].plot(x, group["mean_tmrca_generations"], label=label, color=color)
        axes[0, 0].fill_between(
            x, group["q25_tmrca"], group["q75_tmrca"], color=color, alpha=0.15
        )
        axes[0, 1].plot(x, group["mean_fraction_recent"], label=label, color=color)
        axes[0, 1].fill_between(
            x, group["q25_fraction_recent"], group["q75_fraction_recent"],
            color=color, alpha=0.15,
        )
    for axis in axes[0]:
        axis.axvline(center / 1e6, color="black", ls="--", lw=0.8)
        axis.set_xlim(0, sequence_length / 1e6)
        axis.legend()
    axes[0, 0].set(
        xlabel="Position (Mbp)", ylabel="Mean within-diploid TMRCA (generations)",
        title="Full-region genealogy: selected center and neutral plateaus",
    )
    axes[0, 1].set(
        xlabel="Position (Mbp)", ylabel=f"P(TMRCA < {age_generations})",
        title="Recent-coalescence statistic",
    )
    axes[1, 0].hist(null["center_fraction_recent"], bins=15, color="0.55", alpha=0.8)
    for coefficient, group in stats[stats["selection_coefficient"] > 0].groupby(
        "selection_coefficient"
    ):
        axes[1, 0].axvline(
            group["center_fraction_recent"].mean(), color=colors.get(float(coefficient)),
            label=f"s={coefficient:g}",
        )
    axes[1, 0].set(
        xlabel=f"Neutral center P(TMRCA < {age_generations})", ylabel="Replicates",
        title="Matched neutral Monte Carlo null",
    )
    axes[1, 0].legend()
    selected_stats = stats[stats["selection_coefficient"] > 0]
    axes[1, 1].scatter(
        selected_stats["selection_coefficient"],
        selected_stats["realized_population_allele_frequency"],
        c=[colors.get(float(x), "C0") for x in selected_stats["selection_coefficient"]],
        alpha=0.55,
    )
    max_frequency = max(
        float(selected_stats["realized_population_allele_frequency"].max()), 0.001
    )
    for coefficient, group in selected_stats.groupby("selection_coefficient"):
        lost = int(np.count_nonzero(group["focal_allele_outcome"] == "lost"))
        axes[1, 1].text(
            coefficient, max_frequency * 0.92,
            f"lost {lost}/{len(group)}", ha="center", va="top",
            color=colors.get(float(coefficient)),
        )
    axes[1, 1].set(
        xlabel="Selection coefficient", ylabel="Realized population allele frequency",
        title=f"Unconditional single-origin outcome after {age_generations} generations",
    )
    fig.savefig(output_dir / "recent_sweep_10mb_validation.png", dpi=180)
    plt.close(fig)

    metrics = {
        "population_size": population_size,
        "sample_diploids": sample_diploids,
        "sequence_length": sequence_length,
        "sweep_position": center,
        "selection_coefficients": list(selection_coefficients),
        "age_generations": age_generations,
        "generation_time_years": 25,
        "threshold_years": age_generations * 25,
        "mutation_rate": mutation_rate,
        "recombination_rate": recombination_rate,
        "neutral_replicates": neutral_replicates,
        "selected_replicates_per_s": selected_replicates,
        "trajectory_workers": workers,
        "conditioning": "matched unconditional single-origin SLiM trajectories for s=0 and s>0; lost alleles remain frequency 0 and are still scanned at 5 Mb",
        "elapsed_seconds": perf_counter() - started,
        "results_by_s": stats.groupby("selection_coefficient").mean(numeric_only=True).reset_index().to_dict("records"),
    }
    with (output_dir / "metrics.json").open("w", encoding="utf-8") as handle:
        json.dump(metrics, handle, indent=2)
    summarize_recent_sweep_results(
        output_dir,
        population_size=population_size,
        age_generations=age_generations,
    )
    return metrics


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
