from __future__ import annotations

import json
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass
from pathlib import Path
from time import perf_counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import msprime
import numpy as np
import pandas as pd


@dataclass(frozen=True)
class SimulationConfig:
    n_replicates: int = 1000
    n_diploids: int = 2000
    sequence_length: int = 1_000_000
    effective_size: float = 10_000
    mutation_rate: float = 1.25e-8
    recombination_rate: float = 1.0e-8
    threshold_years: float = 4500
    generation_time: float = 30
    seed: int = 1729
    save_trees: bool = False

    @property
    def threshold_generations(self) -> float:
        return self.threshold_years / self.generation_time


def demography_from_history(times: list[float], sizes: list[float]) -> msprime.Demography:
    if len(times) != len(sizes) or not times or times[0] != 0:
        raise ValueError("history requires equal time/size arrays beginning at time 0")
    if np.any(np.diff(times) <= 0) or np.any(np.asarray(sizes) <= 0):
        raise ValueError("history times must increase and sizes must be positive")
    demography = msprime.Demography()
    demography.add_population(name="pop", initial_size=float(sizes[0]))
    for time, size in zip(times[1:], sizes[1:]):
        demography.add_population_parameters_change(
            time=float(time), initial_size=float(size), population="pop"
        )
    return demography


def draw_log_ne_histories(
    mean_log_ne: np.ndarray,
    covariance: np.ndarray,
    *,
    n: int,
    seed: int,
    minimum_ne: float = 100.0,
    maximum_ne: float = 1e7,
) -> np.ndarray:
    """Draw valid piecewise Ne histories from a PHLASH log-Ne MVN."""
    rng = np.random.default_rng(seed)
    draws = rng.multivariate_normal(mean_log_ne, covariance, size=n)
    return np.clip(np.exp(draws), minimum_ne, maximum_ne)


def _rate_map(path: str | Path | None, default: float, sequence_length: int):
    if path is None:
        return default
    table = pd.read_csv(path, sep="\t")
    if not {"position", "rate"}.issubset(table.columns):
        raise ValueError("rate map TSV must contain position and rate columns")
    positions = table["position"].to_numpy(dtype=float)
    rates = table["rate"].to_numpy(dtype=float)
    if positions[0] != 0 or positions[-1] != sequence_length:
        raise ValueError("rate map positions must begin at 0 and end at sequence_length")
    if len(rates) == len(positions):
        rates = rates[:-1]
    if len(rates) != len(positions) - 1:
        raise ValueError("rate map needs one rate per interval")
    return msprime.RateMap(position=positions, rate=rates)


def within_individual_truth(ts, threshold_generations: float) -> pd.DataFrame:
    pairs = []
    for individual in ts.individuals():
        if len(individual.nodes) == 2:
            pairs.append(tuple(individual.nodes))
    if not pairs:
        raise ValueError("tree sequence has no diploid individuals with exactly two sample nodes")
    rows = []
    for site in ts.sites():
        tree = ts.at(site.position)
        times = np.fromiter((tree.tmrca(a, b) for a, b in pairs), dtype=float)
        rows.append(
            {
                "position_0based": int(site.position),
                "position_1based": int(site.position) + 1,
                "n_pairs": len(times),
                "mean_p_tmrca_lt_threshold": float(np.mean(times < threshold_generations)),
                "mean_tmrca_generations": float(np.mean(times)),
            }
        )
    return pd.DataFrame.from_records(rows)


def _qc_tables(ts) -> tuple[np.ndarray, pd.DataFrame]:
    genotypes = ts.genotype_matrix()
    if genotypes.size:
        sfs = np.bincount(np.sum(genotypes > 0, axis=1), minlength=ts.num_samples + 1)
    else:
        sfs = np.zeros(ts.num_samples + 1, dtype=int)
    positions = np.fromiter((site.position for site in ts.sites()), dtype=float)
    ld_rows = []
    for index in range(max(0, len(positions) - 1)):
        x = genotypes[index].astype(float)
        y = genotypes[index + 1].astype(float)
        if x.std() > 0 and y.std() > 0:
            r = np.corrcoef(x, y)[0, 1]
            ld_rows.append({"distance_bp": positions[index + 1] - positions[index], "r2": r * r})
    return sfs, pd.DataFrame(ld_rows)


def _plot_qc(ts, summary: pd.DataFrame, demography, destination: Path) -> None:
    sfs, ld = _qc_tables(ts)
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
    dbg = demography.debug()
    epochs = list(dbg.epochs)
    times = [epoch.start_time for epoch in epochs]
    sizes = [epoch.populations[0].start_size for epoch in epochs]
    plot_end = max(float(np.max(ts.tables.nodes.time)), 1.0)
    times.append(plot_end)
    sizes.append(sizes[-1])
    axes[0, 0].step(times, sizes, where="post")
    axes[0, 0].set(xlabel="Generations ago", ylabel="Ne", title="Simulated demographic history", yscale="log")
    nonzero = np.flatnonzero(sfs)
    axes[0, 1].scatter(nonzero, sfs[nonzero], s=10)
    axes[0, 1].set(xlabel="Derived allele count", ylabel="Sites", title="Site-frequency spectrum", yscale="log")
    if not ld.empty:
        axes[1, 0].scatter(ld["distance_bp"], ld["r2"], s=8, alpha=0.5)
    axes[1, 0].set(xlabel="Adjacent-site distance (bp)", ylabel="$r^2$", title="Adjacent-site LD")
    if not summary.empty:
        axes[1, 1].plot(summary["position_1based"], summary["mean_tmrca_generations"], lw=0.7)
    axes[1, 1].set(xlabel="Position (bp)", ylabel="Mean within-pair TMRCA", title="True first-coalescence curve")
    fig.savefig(destination, dpi=160)
    plt.close(fig)


def _simulate_one(task):
    (
        replicate, config, output_dir, histories, mutation_map,
        recombination_map, ancestry_seed, mutation_seed,
    ) = task
    output_dir = Path(output_dir)
    if histories:
        times, sizes = histories[replicate % len(histories)]
        demography = demography_from_history(times, sizes)
        history_id = replicate % len(histories)
    else:
        demography = msprime.Demography.isolated_model([config.effective_size])
        demography.populations[0].name = "pop"
        history_id = -1
    ts = msprime.sim_ancestry(
        samples=[msprime.SampleSet(config.n_diploids, population="pop", ploidy=2)],
        demography=demography,
        sequence_length=config.sequence_length,
        recombination_rate=_rate_map(
            recombination_map, config.recombination_rate, config.sequence_length
        ),
        model=msprime.StandardCoalescent(),
        random_seed=ancestry_seed,
    )
    mts = msprime.sim_mutations(
        ts,
        rate=_rate_map(mutation_map, config.mutation_rate, config.sequence_length),
        model=msprime.BinaryMutationModel(),
        random_seed=mutation_seed,
    )
    summary = within_individual_truth(mts, config.threshold_generations)
    summary_path = output_dir / "truth_summaries" / f"replicate_{replicate:05d}.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    if config.save_trees:
        mts.dump(output_dir / "trees" / f"replicate_{replicate:05d}.trees")
    if replicate == 0:
        _plot_qc(mts, summary, demography, output_dir / "plots" / "simulation_qc.png")
    return {
        "replicate": replicate,
        "model": "StandardCoalescent",
        "history_id": history_id,
        "ancestry_seed": ancestry_seed,
        "mutation_seed": mutation_seed,
        "n_diploids": config.n_diploids,
        "n_sites": mts.num_sites,
        "n_trees": mts.num_trees,
        "summary": str(summary_path),
    }


def simulate_replicates(
    config: SimulationConfig,
    output_dir: str | Path,
    *,
    histories: list[tuple[list[float], list[float]]] | None = None,
    mutation_map: str | Path | None = None,
    recombination_map: str | Path | None = None,
    workers: int = 1,
) -> pd.DataFrame:
    """Simulate neutral replicates using only ``msprime.StandardCoalescent``."""
    output_dir = Path(output_dir)
    summary_dir = output_dir / "truth_summaries"
    plot_dir = output_dir / "plots"
    tree_dir = output_dir / "trees"
    summary_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)
    if config.save_trees:
        tree_dir.mkdir(parents=True, exist_ok=True)

    if workers < 1:
        raise ValueError("workers must be positive")
    seed_sequence = np.random.SeedSequence(config.seed)
    seeds = seed_sequence.generate_state(config.n_replicates * 2, dtype=np.uint32).reshape(-1, 2)
    started = perf_counter()
    tasks = [
        (
            replicate, config, str(output_dir), histories,
            str(mutation_map) if mutation_map is not None else None,
            str(recombination_map) if recombination_map is not None else None,
            int(seeds[replicate, 0]), int(seeds[replicate, 1]),
        )
        for replicate in range(config.n_replicates)
    ]
    if workers == 1:
        manifest_rows = [_simulate_one(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            manifest_rows = list(executor.map(_simulate_one, tasks))
    manifest = pd.DataFrame(manifest_rows)
    manifest["elapsed_seconds_total"] = perf_counter() - started
    manifest.to_csv(output_dir / "manifest.tsv", sep="\t", index=False)
    with (output_dir / "config.json").open("w", encoding="utf-8") as handle:
        json.dump(
            {
                **asdict(config),
                "ancestry_model": "StandardCoalescent",
                "mutation_model": "BinaryMutationModel",
                "workers": workers,
            },
            handle,
            indent=2,
        )
    return manifest
