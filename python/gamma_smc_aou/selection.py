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
from tempfile import TemporaryDirectory
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


def within_individual_tmrca_details(
    ts,
    position: float,
    threshold_generations: float,
    focal_carrier_counts: np.ndarray | None = None,
) -> pd.DataFrame:
    """Return one exact-position TMRCA row per sampled diploid individual."""
    sample_nodes = set(ts.samples())
    tree = ts.at(float(position))
    rows = []
    for individual in ts.individuals():
        nodes = [node for node in individual.nodes if node in sample_nodes]
        if len(nodes) != 2:
            continue
        tmrca = float(tree.tmrca(*nodes))
        row = {
            "individual_id": int(individual.id),
            "tmrca_generations": tmrca,
            "tmrca_lt_threshold": bool(tmrca < threshold_generations),
        }
        if focal_carrier_counts is not None:
            row["focal_carrier_copies"] = int(focal_carrier_counts[individual.id])
        rows.append(row)
    if not rows:
        raise ValueError("tree sequence has no sampled diploid individuals")
    return pd.DataFrame(rows)


def _focal_carrier_counts(ts, sweep_position: int) -> np.ndarray:
    """Count focal derived copies in each sampled diploid before site removal."""
    focal_variants = [
        variant
        for variant in ts.variants()
        if np.isclose(variant.site.position, sweep_position)
    ]
    if len(focal_variants) != 1:
        raise ValueError(
            f"expected one retained focal variant at {sweep_position}, "
            f"found {len(focal_variants)}"
        )
    variant = focal_variants[0]
    node_genotype = {
        int(node): int(genotype)
        for node, genotype in zip(ts.samples(), variant.genotypes)
    }
    sample_nodes = set(node_genotype)
    counts = np.full(ts.num_individuals, -1, dtype=np.int8)
    for individual in ts.individuals():
        nodes = [node for node in individual.nodes if node in sample_nodes]
        if len(nodes) == 2:
            counts[individual.id] = sum(node_genotype[node] > 0 for node in nodes)
    if np.any(counts < 0):
        raise ValueError("not every retained individual has two focal genotypes")
    return counts


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
    capture_focal_genotypes: bool = False,
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
    focal_carrier_counts = (
        _focal_carrier_counts(sampled, sweep_position)
        if capture_focal_genotypes
        else None
    )
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
    run = {
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
    if focal_carrier_counts is not None:
        run["focal_carrier_counts"] = focal_carrier_counts
        run["sample_focal_allele_frequency"] = float(
            focal_carrier_counts.sum() / (2 * len(focal_carrier_counts))
        )
    return sampled, run


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
        lost_group = group[group["focal_allele_outcome"] == "lost"]
        present_group = group[
            group["focal_allele_outcome"].isin(("segregating", "fixed"))
        ]
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
            "n_focal_allele_lost": int(len(lost_group)),
            "n_focal_allele_present": int(len(present_group)),
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
            if len(present_group):
                present_u = mannwhitneyu(
                    present_group[column],
                    null[column],
                    alternative="greater",
                    method="auto",
                ).statistic
                base[f"{column}_mean_given_focal_allele_present"] = float(
                    present_group[column].mean()
                )
                base[f"{column}_median_p_given_focal_allele_present"] = float(
                    present_group[p_column].median()
                )
                base[
                    f"{column}_power_p_le_0_05_given_focal_allele_present"
                ] = float(np.mean(present_group[p_column] <= 0.05))
                base[f"{column}_auc_given_focal_allele_present"] = float(
                    present_u / (len(present_group) * len(null))
                )
            else:
                base[f"{column}_mean_given_focal_allele_present"] = np.nan
                base[f"{column}_median_p_given_focal_allele_present"] = np.nan
                base[
                    f"{column}_power_p_le_0_05_given_focal_allele_present"
                ] = np.nan
                base[f"{column}_auc_given_focal_allele_present"] = np.nan
            base[f"{column}_false_positive_p_le_0_05_given_focal_allele_lost"] = (
                float(np.mean(lost_group[p_column] <= 0.05))
                if len(lost_group)
                else np.nan
            )
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
    reuse_null_from: str | Path | None = None,
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
    if reuse_null_from is None and not any(x == 0 for x in selection_coefficients):
        raise ValueError("selection_coefficients must include 0 unless reusing a null")

    reused_null_path = None
    if reuse_null_from is not None:
        reused_null_path = Path(reuse_null_from).resolve()
        if reused_null_path == output_dir.resolve():
            raise ValueError("reuse_null_from must differ from output_dir")
        with (reused_null_path / "metrics.json").open(encoding="utf-8") as handle:
            source_metrics = json.load(handle)
        expected = {
            "population_size": population_size,
            "sample_diploids": sample_diploids,
            "sequence_length": sequence_length,
            "age_generations": age_generations,
        }
        for key, value in expected.items():
            if source_metrics.get(key) != value:
                raise ValueError(
                    f"reused null {key}={source_metrics.get(key)!r}, expected {value!r}"
                )
        for key, value in {
            "mutation_rate": mutation_rate,
            "recombination_rate": recombination_rate,
        }.items():
            if not np.isclose(float(source_metrics.get(key, np.nan)), value):
                raise ValueError(
                    f"reused null {key}={source_metrics.get(key)!r}, expected {value!r}"
                )
        source_stats = pd.read_csv(
            reused_null_path / "replicate_statistics.tsv", sep="\t"
        )
        source_profiles = pd.read_csv(reused_null_path / "truth_profiles.tsv", sep="\t")
        source_stats = source_stats[source_stats["selection_coefficient"] == 0].copy()
        source_profiles = source_profiles[
            source_profiles["selection_coefficient"] == 0
        ].copy()
        if source_stats.empty or source_profiles.empty:
            raise ValueError("reuse_null_from contains no s=0 replicates")
        neutral_replicates = int(source_stats["replicate"].nunique())
        profiles.append(source_profiles)
        rows.extend(source_stats.to_dict("records"))

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
        if not (reuse_null_from is not None and coefficient == 0)
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
    colors = {
        0.0: "0.35", 0.001: "#3b82f6", 0.01: "#dc2626", 0.1: "#7c3aed"
    }
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
        "selection_coefficients": sorted(
            float(x) for x in stats["selection_coefficient"].unique()
        ),
        "age_generations": age_generations,
        "generation_time_years": 25,
        "threshold_years": age_generations * 25,
        "mutation_rate": mutation_rate,
        "recombination_rate": recombination_rate,
        "neutral_replicates": neutral_replicates,
        "selected_replicates_per_s": selected_replicates,
        "trajectory_workers": workers,
        "seed": int(seed),
        "neutral_reused_from": (
            str(reused_null_path) if reused_null_path is not None else None
        ),
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


def retained_sweep_calibration_table(
    stats: pd.DataFrame,
    *,
    selection_coefficient: float,
    retained_replicates: int,
    sample_diploids: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Select retained trajectories and compare each with the full neutral null."""
    neutral = stats[np.isclose(stats["selection_coefficient"], 0)].copy()
    eligible = stats[
        np.isclose(stats["selection_coefficient"], selection_coefficient)
        & stats["focal_allele_outcome"].isin(("segregating", "fixed"))
    ].sort_values("replicate")
    if len(neutral) < 2:
        raise ValueError("at least two neutral replicates are required")
    if len(eligible) < retained_replicates:
        raise ValueError(
            f"requested {retained_replicates} retained trajectories but only "
            f"{len(eligible)} are available"
        )
    selected = eligible.head(retained_replicates).copy()
    neutral["n_pairs_tmrca_lt_threshold"] = np.rint(
        neutral["center_fraction_recent"] * sample_diploids
    ).astype(int)
    rows = []
    null_values = neutral["center_fraction_recent"].to_numpy(dtype=float)
    for row in selected.itertuples(index=False):
        value = float(row.center_fraction_recent)
        exceedances = int(np.count_nonzero(null_values >= value))
        rows.append({
            "replicate": int(row.replicate),
            "realized_population_allele_frequency": float(
                row.realized_population_allele_frequency
            ),
            "center_fraction_recent": value,
            "n_pairs_tmrca_lt_threshold": int(round(value * sample_diploids)),
            "neutral_exceedances": exceedances,
            "neutral_replicates": int(len(neutral)),
            "mc_p_upper": float((1 + exceedances) / (1 + len(neutral))),
        })
    return neutral, pd.DataFrame(rows)


def _plot_retained_calibration(
    neutral: pd.DataFrame,
    selected: pd.DataFrame,
    *,
    threshold_generations: int,
    output_path: Path,
) -> None:
    null_counts = neutral["n_pairs_tmrca_lt_threshold"].to_numpy(dtype=int)
    selected_counts = selected["n_pairs_tmrca_lt_threshold"].to_numpy(dtype=int)
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), constrained_layout=True)

    bins = np.arange(null_counts.min() - 0.5, null_counts.max() + 1.5)
    axes[0].hist(null_counts, bins=bins, color="0.55", edgecolor="white")
    axes[0].axvline(
        np.quantile(null_counts, 0.95), color="#dc2626", ls="--", lw=1.2,
        label="neutral 95th percentile",
    )
    axes[0].set(
        xlabel=f"Number of pairs with TMRCA < {threshold_generations}",
        ylabel="Neutral simulations",
        title=f"Full neutral null (n={len(neutral)})",
    )
    axes[0].legend(fontsize=8)

    thresholds = np.unique(np.r_[1, null_counts, selected_counts])
    exceedance_curve = np.asarray([
        np.count_nonzero(null_counts >= threshold) for threshold in thresholds
    ])
    axes[1].step(thresholds, exceedance_curve, where="post", color="0.25")
    axes[1].scatter(
        selected_counts,
        selected["neutral_exceedances"],
        color="#7c3aed",
        s=34,
        zorder=3,
        label="retained s=0.1",
    )
    axes[1].set_xscale("log")
    axes[1].set(
        xlabel=f"Observed number of pairs with TMRCA < {threshold_generations}",
        ylabel=f"Number of {len(neutral)} neutral simulations >= observed",
        title="Empirical upper-tail count",
    )
    axes[1].legend(fontsize=8)

    x = np.arange(len(selected))
    axes[2].axhspan(
        null_counts.min(), null_counts.max(), color="0.75", alpha=0.45,
        label="neutral range",
    )
    axes[2].scatter(x, selected_counts, color="#7c3aed", s=42, zorder=3)
    for index, row in enumerate(selected.itertuples(index=False)):
        axes[2].annotate(
            f"{row.neutral_exceedances}/{row.neutral_replicates}\n"
            f"p={row.mc_p_upper:.4f}",
            (index, row.n_pairs_tmrca_lt_threshold),
            xytext=(0, 6),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=7,
        )
    axes[2].set_yscale("log")
    axes[2].set_ylim(
        max(1, min(null_counts) * 0.7),
        max(selected_counts) * 1.8,
    )
    axes[2].set_xticks(x, [str(value) for value in selected["replicate"]])
    axes[2].set(
        xlabel="Retained selected replicate",
        ylabel=f"Pairs with TMRCA < {threshold_generations}",
        title="Selected statistics and Monte Carlo p-values",
    )
    axes[2].legend(fontsize=8)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def _plot_carrier_tmrca(
    details: pd.DataFrame,
    selected: pd.DataFrame,
    *,
    threshold_generations: int,
    output_path: Path,
) -> None:
    fig, axes = plt.subplots(5, 2, figsize=(13, 17), sharex=True, sharey=True)
    axes = axes.ravel()
    for axis, selected_row in zip(axes, selected.itertuples(index=False)):
        group = details[details["replicate"] == selected_row.replicate]
        for copies, label, color in (
            (0, "noncarrier / noncarrier", "0.4"),
            (1, "carrier / noncarrier", "#60a5fa"),
            (2, "carrier / carrier", "#7c3aed"),
        ):
            values = np.sort(
                group.loc[
                    group["focal_carrier_copies"] == copies,
                    "tmrca_generations",
                ].to_numpy(dtype=float)
            )
            if not len(values):
                continue
            ecdf = np.arange(1, len(values) + 1) / len(values)
            axis.step(values, ecdf, where="post", color=color, lw=1.4, label=label)
        recent_by_copy = group.groupby("focal_carrier_copies")[
            "tmrca_lt_threshold"
        ].mean()
        n_by_copy = group["focal_carrier_copies"].value_counts()
        axis.axvline(threshold_generations, color="#dc2626", ls="--", lw=0.9)
        axis.text(
            0.98,
            0.04,
            "P(recent), copies 0/1/2: "
            + "/".join(f"{recent_by_copy.get(i, np.nan):.3f}" for i in range(3))
            + "\nn, copies 0/1/2: "
            + "/".join(str(int(n_by_copy.get(i, 0))) for i in range(3)),
            transform=axis.transAxes,
            ha="right",
            va="bottom",
            fontsize=8,
        )
        axis.set_title(
            f"replicate {selected_row.replicate}; population AF="
            f"{selected_row.realized_population_allele_frequency:.3f}"
        )
        axis.grid(alpha=0.15)
    for axis in axes:
        axis.set_xscale("log")
        axis.set_ylim(0, 1)
    for axis in axes[::2]:
        axis.set_ylabel("Cumulative fraction of pairs")
    for axis in axes[-2:]:
        axis.set_xlabel("Within-diploid TMRCA (generations, log scale)")
    axes[0].legend(loc="upper left", fontsize=8)
    fig.suptitle(
        "Within each selected simulation: carrier-carrier versus noncarrier pairs",
        fontsize=14,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def summarize_carrier_tmrca(
    details: pd.DataFrame,
    selected: pd.DataFrame,
    *,
    threshold_generations: int,
    output_dir: str | Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Write carrier-copy summaries and the direct 2-copy versus 0-copy test."""
    from scipy.stats import mannwhitneyu

    output_dir = Path(output_dir)
    genotype_summary = (
        details.groupby(["replicate", "focal_carrier_copies"], as_index=False)
        .agg(
            n_pairs=("tmrca_generations", "size"),
            mean_tmrca_generations=("tmrca_generations", "mean"),
            median_tmrca_generations=("tmrca_generations", "median"),
            fraction_tmrca_lt_threshold=("tmrca_lt_threshold", "mean"),
        )
    )
    genotype_summary["pair_class"] = genotype_summary[
        "focal_carrier_copies"
    ].map({
        0: "noncarrier_noncarrier",
        1: "carrier_noncarrier",
        2: "carrier_carrier",
    })
    genotype_summary.to_csv(
        output_dir / "carrier_copy_summary.tsv", sep="\t", index=False
    )
    direct_summary = genotype_summary[
        genotype_summary["focal_carrier_copies"].isin((0, 2))
    ].copy()
    direct_summary.to_csv(
        output_dir / "carrier_vs_noncarrier_summary.tsv", sep="\t", index=False
    )
    effect_rows = []
    for replicate, group in details.groupby("replicate"):
        noncarrier = group.loc[
            group["focal_carrier_copies"] == 0, "tmrca_generations"
        ]
        mixed = group.loc[
            group["focal_carrier_copies"] == 1, "tmrca_generations"
        ]
        carrier = group.loc[
            group["focal_carrier_copies"] == 2, "tmrca_generations"
        ]
        u_result = mannwhitneyu(
            carrier, noncarrier, alternative="less", method="auto"
        )
        effect_rows.append({
            "replicate": int(replicate),
            "n_noncarrier_noncarrier_pairs": int(len(noncarrier)),
            "n_mixed_pairs": int(len(mixed)),
            "n_carrier_carrier_pairs": int(len(carrier)),
            "noncarrier_fraction_tmrca_lt_threshold": float(
                np.mean(noncarrier < threshold_generations)
            ),
            "carrier_fraction_tmrca_lt_threshold": float(
                np.mean(carrier < threshold_generations)
            ),
            "mann_whitney_u": float(u_result.statistic),
            "mann_whitney_p_carrier_lower": float(u_result.pvalue),
        })
    effects = pd.DataFrame(effect_rows)
    effects.to_csv(output_dir / "carrier_tmrca_effect_tests.tsv", sep="\t", index=False)
    _plot_carrier_tmrca(
        details,
        selected,
        threshold_generations=threshold_generations,
        output_path=output_dir / "carrier_vs_noncarrier_tmrca_ecdf.png",
    )
    return genotype_summary, direct_summary, effects


def analyze_retained_recent_sweeps(
    source_dir: str | Path,
    output_dir: str | Path,
    *,
    executable: str | Path | None = None,
    selection_coefficient: float = 0.1,
    retained_replicates: int = 10,
    workers: int = 20,
    seed: int = 424242,
) -> dict:
    """Calibrate retained sweeps and reconstruct carrier-resolved pair TMRCAs."""
    if workers < 1 or retained_replicates < 1:
        raise ValueError("workers and retained_replicates must be positive")
    source_dir = Path(source_dir).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    with (source_dir / "metrics.json").open(encoding="utf-8") as handle:
        source_metrics = json.load(handle)
    if source_metrics.get("seed") is not None and int(source_metrics["seed"]) != seed:
        raise ValueError(
            f"source seed is {source_metrics['seed']}, but reconstruction seed is {seed}"
        )
    stats = pd.read_csv(source_dir / "replicate_statistics.tsv", sep="\t")
    sample_diploids = int(source_metrics["sample_diploids"])
    threshold_generations = int(source_metrics["age_generations"])
    neutral, selected = retained_sweep_calibration_table(
        stats,
        selection_coefficient=selection_coefficient,
        retained_replicates=retained_replicates,
        sample_diploids=sample_diploids,
    )
    neutral.to_csv(output_dir / "neutral_center_statistics.tsv", sep="\t", index=False)
    selected.to_csv(
        output_dir / "retained_selected_calibration.tsv", sep="\t", index=False
    )
    _plot_retained_calibration(
        neutral,
        selected,
        threshold_generations=threshold_generations,
        output_path=output_dir / "neutral_vs_retained_selected_calibration.png",
    )

    started = perf_counter()
    selected_by_replicate = selected.set_index("replicate")
    center = int(source_metrics["sweep_position"])

    with TemporaryDirectory(prefix="gamma_smc_retained_") as temporary_directory:
        temporary_directory = Path(temporary_directory)

        def reconstruct(replicate: int) -> pd.DataFrame:
            expected = selected_by_replicate.loc[replicate]
            run_seed = (
                seed
                + 1_000_000
                + int(selection_coefficient * 1e8)
                + int(replicate) * 10
            )
            ts, run = run_slim_recent_sweep(
                temporary_directory / f"replicate_{replicate}.trees",
                executable=executable,
                population_size=int(source_metrics["population_size"]),
                sample_diploids=sample_diploids,
                sequence_length=int(source_metrics["sequence_length"]),
                sweep_position=center,
                selection_coefficient=selection_coefficient,
                age_generations=threshold_generations,
                mutation_rate=float(source_metrics["mutation_rate"]),
                recombination_rate=float(source_metrics["recombination_rate"]),
                seed=run_seed,
                capture_focal_genotypes=True,
            )
            if run["focal_allele_outcome"] not in ("segregating", "fixed"):
                raise RuntimeError(
                    f"replicate {replicate} did not reproduce a retained mutation"
                )
            if not np.isclose(
                run["realized_population_allele_frequency"],
                expected["realized_population_allele_frequency"],
                atol=1e-9,
                rtol=0,
            ):
                raise RuntimeError(f"replicate {replicate} allele frequency changed")
            details = within_individual_tmrca_details(
                ts,
                center,
                threshold_generations,
                focal_carrier_counts=run["focal_carrier_counts"],
            )
            observed_recent = float(details["tmrca_lt_threshold"].mean())
            if not np.isclose(
                observed_recent, expected["center_fraction_recent"], atol=1e-12, rtol=0
            ):
                raise RuntimeError(f"replicate {replicate} center statistic changed")
            details["replicate"] = int(replicate)
            details["population_focal_allele_frequency"] = float(
                run["realized_population_allele_frequency"]
            )
            details["sample_focal_allele_frequency"] = float(
                run["sample_focal_allele_frequency"]
            )
            details["pair_has_focal_allele"] = details["focal_carrier_copies"] > 0
            return details

        replicate_ids = selected["replicate"].astype(int).tolist()
        if workers == 1:
            detail_tables = [reconstruct(replicate) for replicate in replicate_ids]
        else:
            with ThreadPoolExecutor(max_workers=workers) as executor:
                detail_tables = list(executor.map(reconstruct, replicate_ids))

    details = pd.concat(detail_tables, ignore_index=True)
    details.to_csv(output_dir / "carrier_pair_tmrca.tsv", sep="\t", index=False)
    _, _, effects = summarize_carrier_tmrca(
        details,
        selected,
        threshold_generations=threshold_generations,
        output_dir=output_dir,
    )

    result = {
        "source_dir": str(source_dir),
        "selection_coefficient": float(selection_coefficient),
        "retained_replicates": replicate_ids,
        "neutral_replicates": int(len(neutral)),
        "sample_diploids": sample_diploids,
        "threshold_generations": threshold_generations,
        "threshold_years": int(source_metrics["threshold_years"]),
        "workers_requested": int(workers),
        "workers_used": int(min(workers, retained_replicates)),
        "base_seed": int(seed),
        "deterministic_reconstruction_verified": True,
        "elapsed_seconds": float(perf_counter() - started),
        "calibration": selected.to_dict("records"),
        "carrier_effect_tests": effects.to_dict("records"),
    }
    with (output_dir / "metrics.json").open("w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2)
    return result


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
