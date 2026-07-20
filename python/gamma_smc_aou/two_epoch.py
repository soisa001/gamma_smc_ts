from __future__ import annotations

import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from tempfile import TemporaryDirectory
from time import perf_counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .calibration import (
    calibration_metrics,
    monte_carlo_pvalue,
    randomized_rank_pvalue,
)
from .carrier_profiles import (
    _position_grids,
    pair_tmrca_profile_by_focal_copy,
    plot_carrier_profile_figure,
)
from .selection import run_slim_recent_sweep, within_individual_tmrca_details


LARGE_FONTS = {
    "title": 24,
    "axis": 20,
    "tick": 17,
    "legend": 16,
    "annotation": 18,
}


def two_epoch_recent_probability(
    *,
    ancestral_population_size: int,
    present_population_size: int,
    size_change_generations_ago: int,
    threshold_generations: int,
) -> float:
    """Pairwise coalescence probability below a threshold under two epochs."""
    recent_duration = min(threshold_generations, size_change_generations_ago)
    ancestral_duration = max(0, threshold_generations - recent_duration)
    cumulative_hazard = (
        recent_duration / (2 * present_population_size)
        + ancestral_duration / (2 * ancestral_population_size)
    )
    return float(1 - np.exp(-cumulative_hazard))


def _plot_demography_and_variant(
    output_path: Path,
    *,
    ancestral_population_size: int,
    present_population_size: int,
    size_change_generations_ago: int,
    variant_age_generations: int,
    generation_time_years: float,
    selection_coefficient: float,
) -> None:
    oldest = max(variant_age_generations + 70, size_change_generations_ago + 100)
    x = np.asarray([0, size_change_generations_ago, size_change_generations_ago, oldest])
    y = np.asarray([
        present_population_size,
        present_population_size,
        ancestral_population_size,
        ancestral_population_size,
    ])
    fig, axis = plt.subplots(figsize=(14, 7.5), constrained_layout=True)
    axis.plot(x, y, color="#0072b2", lw=2.1, label="population size")
    axis.fill_between(x, y, step="pre", alpha=0.12, color="#56b4e9")
    axis.axvline(
        size_change_generations_ago,
        color="#009e73",
        ls="--",
        lw=1.3,
        label=f"growth: {size_change_generations_ago} generations ago",
    )
    axis.axvline(
        variant_age_generations,
        color="#e69f00",
        ls="--",
        lw=1.3,
        label=(
            f"selected variant added: {variant_age_generations} generations "
            f"({variant_age_generations * generation_time_years:g} years) ago"
        ),
    )
    axis.scatter(
        [variant_age_generations],
        [ancestral_population_size],
        color="#e69f00",
        zorder=3,
    )
    axis.text(
        size_change_generations_ago / 2,
        present_population_size * 1.025,
        f"present epoch: N={present_population_size:,}",
        ha="center",
        va="bottom",
        fontsize=LARGE_FONTS["annotation"],
    )
    axis.text(
        (size_change_generations_ago + oldest) / 2,
        ancestral_population_size * 0.95,
        f"ancestral epoch: N={ancestral_population_size:,}",
        ha="center",
        va="top",
        fontsize=LARGE_FONTS["annotation"],
    )
    axis.set(
        xlim=(0, oldest),
        ylim=(0, present_population_size * 1.22),
        xlabel="Generations before present",
        ylabel="Population size",
        title=f"Two-epoch recent growth and selected-variant timing (s={selection_coefficient:g})",
    )
    axis.set_title(
        f"Two-epoch recent growth and selected-variant timing (s={selection_coefficient:g})",
        fontsize=LARGE_FONTS["title"],
    )
    axis.set_xlabel("Generations before present", fontsize=LARGE_FONTS["axis"])
    axis.set_ylabel("Population size", fontsize=LARGE_FONTS["axis"])
    axis.tick_params(axis="both", labelsize=LARGE_FONTS["tick"])
    axis.legend(loc="lower right", fontsize=LARGE_FONTS["legend"])
    axis.grid(axis="y", alpha=0.18)
    fig.savefig(output_path, dpi=190)
    plt.close(fig)


def _plot_null_calibration(
    neutral: pd.DataFrame,
    selected_fraction_recent: float,
    *,
    theoretical_fraction_recent: float,
    threshold_generations: int,
    generation_time_years: float,
    selection_coefficient: float,
    output_path: Path,
) -> tuple[int, float]:
    values = neutral["center_fraction_recent"].to_numpy(dtype=float)
    exceedances = int(np.count_nonzero(values >= selected_fraction_recent))
    pvalue = float((1 + exceedances) / (1 + len(values)))
    fig, axes = plt.subplots(1, 2, figsize=(18, 7), constrained_layout=True)
    axes[0].hist(values, bins=16, color="0.55", edgecolor="white", alpha=0.9)
    axes[0].axvline(
        theoretical_fraction_recent,
        color="#009e73",
        ls="--",
        lw=1.3,
        label=f"two-epoch expectation={theoretical_fraction_recent:.4f}",
    )
    axes[0].axvline(
        selected_fraction_recent,
        color="#7c3aed",
        lw=2,
        label=f"retained s={selection_coefficient:g} observation={selected_fraction_recent:.4f}",
    )
    threshold_years = threshold_generations * generation_time_years
    axes[0].set(
        xlabel=(
            f"P(TMRCA < {threshold_years:g} years; "
            f"{threshold_generations} generations) across 2,000 pairs"
        ),
        ylabel="Neutral simulations",
        title=f"Two-epoch neutral null (n={len(neutral)})",
    )
    axes[0].legend(fontsize=LARGE_FONTS["legend"])

    thresholds = np.unique(np.r_[values, selected_fraction_recent])
    counts = np.asarray([np.count_nonzero(values >= value) for value in thresholds])
    axes[1].step(thresholds, counts, where="post", color="0.25", lw=1.4)
    axes[1].scatter(
        [selected_fraction_recent],
        [exceedances],
        color="#7c3aed",
        s=55,
        zorder=3,
    )
    axes[1].annotate(
        f"{exceedances}/{len(neutral)} neutral >= observed\np={pvalue:.4f}",
        (selected_fraction_recent, exceedances),
        xytext=(-8, 10),
        textcoords="offset points",
        ha="right",
        va="bottom",
        fontsize=LARGE_FONTS["annotation"],
    )
    axes[1].set(
        xlabel=f"Observed P(TMRCA < {threshold_years:g} years)",
        ylabel="Neutral simulations >= observed",
        title="Empirical upper-tail Monte Carlo test",
    )
    axes[1].grid(alpha=0.15)
    for axis in axes:
        axis.title.set_fontsize(LARGE_FONTS["title"])
        axis.xaxis.label.set_fontsize(LARGE_FONTS["axis"])
        axis.yaxis.label.set_fontsize(LARGE_FONTS["axis"])
        axis.tick_params(axis="both", labelsize=LARGE_FONTS["tick"])
    fig.savefig(output_path, dpi=190)
    plt.close(fig)
    return exceedances, pvalue


def _plot_allele_frequency_trajectory(
    trajectory: pd.DataFrame,
    *,
    population_allele_frequency: float,
    sample_allele_frequency: float,
    variant_age_generations: int,
    size_change_generations_ago: int,
    generation_time_years: float,
    selection_coefficient: float,
    output_path: Path,
) -> None:
    """Plot the realized selected-allele frequency from origin to sampling."""
    ordered = trajectory.sort_values("generations_ago", ascending=False)
    fig, axis = plt.subplots(figsize=(14, 7.5), constrained_layout=True)
    axis.plot(
        ordered["generations_ago"],
        ordered["population_allele_frequency"],
        color="#7c3aed",
        lw=2.8,
        label="population allele frequency",
    )
    axis.axvline(
        size_change_generations_ago,
        color="#009e73",
        ls="--",
        lw=1.8,
        label=f"growth: {size_change_generations_ago} generations ago",
    )
    axis.scatter(
        [variant_age_generations, 0],
        [
            ordered.iloc[0]["population_allele_frequency"],
            population_allele_frequency,
        ],
        color=["#e69f00", "#7c3aed"],
        s=90,
        zorder=3,
    )
    axis.annotate(
        "single-copy origin\n"
        f"{variant_age_generations * generation_time_years:g} years ago",
        (
            variant_age_generations,
            ordered.iloc[0]["population_allele_frequency"],
        ),
        xytext=(-18, 20),
        textcoords="offset points",
        ha="right",
        fontsize=LARGE_FONTS["annotation"],
    )
    axis.annotate(
        f"present population AF={population_allele_frequency:.4f}\n"
        f"sample AF={sample_allele_frequency:.4f}",
        (0, population_allele_frequency),
        xytext=(18, -6),
        textcoords="offset points",
        ha="left",
        va="center",
        fontsize=LARGE_FONTS["annotation"],
    )
    axis.set_xlim(-15, variant_age_generations + 5)
    axis.set_ylim(0, max(0.12, population_allele_frequency * 1.35))
    axis.set_title(
        f"Realized selected-allele trajectory (s={selection_coefficient:g})",
        fontsize=LARGE_FONTS["title"],
    )
    axis.set_xlabel(
        "Generations before present",
        fontsize=LARGE_FONTS["axis"],
    )
    axis.set_ylabel("Selected-allele frequency", fontsize=LARGE_FONTS["axis"])
    axis.tick_params(axis="both", labelsize=LARGE_FONTS["tick"])
    axis.legend(loc="upper left", fontsize=LARGE_FONTS["legend"])
    axis.grid(alpha=0.18)
    fig.savefig(output_path, dpi=190)
    plt.close(fig)


def validate_two_epoch_growth(
    output_dir: str | Path,
    *,
    executable: str | Path | None = None,
    ancestral_population_size: int = 10_000,
    present_population_size: int = 20_000,
    size_change_generations_ago: int = 100,
    sample_diploids: int = 2_000,
    sequence_length: int = 10_000_000,
    variant_age_generations: int = 180,
    generation_time_years: float = 25,
    selection_coefficient: float = 0.05,
    mutation_rate: float = 1.25e-8,
    recombination_rate: float = 1e-8,
    neutral_replicates: int = 100,
    workers: int = 20,
    max_selected_attempts: int = 1_000,
    minimum_hom_alt_pairs: int = 2,
    full_step: int = 50_000,
    zoom_half_width: int = 500_000,
    zoom_step: int = 5_000,
    seed: int = 515151,
) -> dict:
    """Run a two-epoch neutral null and rejection-sampled selected validation."""
    if workers < 1 or neutral_replicates < 2 or max_selected_attempts < 1:
        raise ValueError("workers, replicates, and attempts must be positive")
    if not 0 < size_change_generations_ago < variant_age_generations:
        raise ValueError("population growth must occur after variant addition")
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    center = sequence_length // 2
    _, _, profile_positions = _position_grids(
        sequence_length,
        center,
        full_step=full_step,
        zoom_half_width=zoom_half_width,
        zoom_step=zoom_step,
    )
    started = perf_counter()

    with TemporaryDirectory(prefix="gamma_smc_two_epoch_") as temporary_directory:
        temporary_directory = Path(temporary_directory)

        def simulate_neutral(replicate: int) -> dict:
            run_seed = seed + 1_000_000 + replicate * 10
            ts, run = run_slim_recent_sweep(
                temporary_directory / f"neutral_{replicate}.trees",
                executable=executable,
                population_size=ancestral_population_size,
                present_population_size=present_population_size,
                size_change_generations_ago=size_change_generations_ago,
                sample_diploids=sample_diploids,
                sequence_length=sequence_length,
                sweep_position=center,
                selection_coefficient=0.0,
                age_generations=variant_age_generations,
                mutation_rate=mutation_rate,
                recombination_rate=recombination_rate,
                seed=run_seed,
            )
            details = within_individual_tmrca_details(
                ts, center, variant_age_generations
            )
            return {
                "replicate": int(replicate),
                "center_fraction_recent": float(
                    details["tmrca_lt_threshold"].mean()
                ),
                "center_n_pairs_recent": int(
                    details["tmrca_lt_threshold"].sum()
                ),
                "center_mean_tmrca_generations": float(
                    details["tmrca_generations"].mean()
                ),
                "focal_allele_outcome": run["focal_allele_outcome"],
                "focal_population_allele_frequency": float(
                    run["realized_population_allele_frequency"]
                ),
                "ancestry_seconds": float(run["ancestry_seconds"]),
                "slim_seconds": float(run["slim_seconds"]),
            }

        neutral_ids = list(range(neutral_replicates))
        if workers == 1:
            neutral_rows = [simulate_neutral(value) for value in neutral_ids]
        else:
            with ThreadPoolExecutor(max_workers=workers) as executor:
                neutral_rows = list(executor.map(simulate_neutral, neutral_ids))
        neutral = pd.DataFrame(neutral_rows).sort_values("replicate")
        neutral.to_csv(output_dir / "neutral_statistics.tsv", sep="\t", index=False)

        def simulate_selected(
            attempt: int,
        ) -> tuple[dict, pd.DataFrame | None, pd.DataFrame]:
            run_seed = seed + 10_000_000 + attempt * 10
            ts, run = run_slim_recent_sweep(
                temporary_directory / f"selected_attempt_{attempt}.trees",
                executable=executable,
                population_size=ancestral_population_size,
                present_population_size=present_population_size,
                size_change_generations_ago=size_change_generations_ago,
                sample_diploids=sample_diploids,
                sequence_length=sequence_length,
                sweep_position=center,
                selection_coefficient=selection_coefficient,
                age_generations=variant_age_generations,
                mutation_rate=mutation_rate,
                recombination_rate=recombination_rate,
                seed=run_seed,
                capture_focal_genotypes=True,
            )
            carrier_counts = run["focal_carrier_counts"]
            n_hom_ref = int(np.count_nonzero(carrier_counts == 0))
            n_heterozygous = int(np.count_nonzero(carrier_counts == 1))
            n_hom_alt = int(np.count_nonzero(carrier_counts == 2))
            details = within_individual_tmrca_details(
                ts,
                center,
                variant_age_generations,
                focal_carrier_counts=carrier_counts,
            )
            accepted = bool(
                run["focal_allele_outcome"] in ("segregating", "fixed")
                and n_hom_alt >= minimum_hom_alt_pairs
                and n_hom_ref >= 2
            )
            row = {
                "attempt": int(attempt),
                "accepted": accepted,
                "focal_allele_outcome": run["focal_allele_outcome"],
                "population_allele_frequency": float(
                    run["realized_population_allele_frequency"]
                ),
                "sample_allele_frequency": float(
                    run["sample_focal_allele_frequency"]
                ),
                "n_hom_ref_pairs": n_hom_ref,
                "n_heterozygous_pairs": n_heterozygous,
                "n_hom_alt_pairs": n_hom_alt,
                "center_fraction_recent": float(
                    details["tmrca_lt_threshold"].mean()
                ),
                "center_n_pairs_recent": int(
                    details["tmrca_lt_threshold"].sum()
                ),
                "center_mean_tmrca_generations": float(
                    details["tmrca_generations"].mean()
                ),
                "ancestry_seconds": float(run["ancestry_seconds"]),
                "slim_seconds": float(run["slim_seconds"]),
            }
            profile = None
            if accepted:
                profile = pair_tmrca_profile_by_focal_copy(
                    ts, profile_positions, carrier_counts
                )
                profile["selected_attempt"] = int(attempt)
                profile["population_focal_allele_frequency"] = float(
                    run["realized_population_allele_frequency"]
                )
                profile["sample_focal_allele_frequency"] = float(
                    run["sample_focal_allele_frequency"]
                )
            trajectory = pd.DataFrame(run["allele_frequency_trajectory"])
            trajectory["selected_attempt"] = int(attempt)
            return row, profile, trajectory

        attempts = []
        accepted_row = None
        accepted_profile = None
        accepted_trajectory = None
        next_attempt = 0
        while accepted_row is None and next_attempt < max_selected_attempts:
            batch = list(
                range(
                    next_attempt,
                    min(next_attempt + workers, max_selected_attempts),
                )
            )
            if workers == 1:
                results = [simulate_selected(value) for value in batch]
            else:
                with ThreadPoolExecutor(max_workers=workers) as executor:
                    results = list(executor.map(simulate_selected, batch))
            for row, profile, trajectory in results:
                attempts.append(row)
                if accepted_row is None and row["accepted"]:
                    accepted_row = row
                    accepted_profile = profile
                    accepted_trajectory = trajectory
            next_attempt += len(batch)
        if (
            accepted_row is None
            or accepted_profile is None
            or accepted_trajectory is None
        ):
            raise RuntimeError(
                f"no analyzable retained s={selection_coefficient:g} allele in "
                f"{max_selected_attempts} attempts"
            )

    attempts_table = pd.DataFrame(attempts).sort_values("attempt")
    attempts_table["used_for_rejection_decision"] = (
        attempts_table["attempt"] <= int(accepted_row["attempt"])
    )
    attempts_table.to_csv(
        output_dir / "selected_rejection_attempts.tsv", sep="\t", index=False
    )
    pd.DataFrame([accepted_row]).to_csv(
        output_dir / "selected_observation.tsv", sep="\t", index=False
    )
    accepted_profile.to_csv(
        output_dir / "selected_hom_alt_vs_hom_ref_tmrca_profile.tsv",
        sep="\t",
        index=False,
    )
    accepted_trajectory.to_csv(
        output_dir / "selected_allele_frequency_trajectory.tsv",
        sep="\t",
        index=False,
    )

    theoretical_recent = two_epoch_recent_probability(
        ancestral_population_size=ancestral_population_size,
        present_population_size=present_population_size,
        size_change_generations_ago=size_change_generations_ago,
        threshold_generations=variant_age_generations,
    )
    _plot_demography_and_variant(
        output_dir / "demography_and_variant_timing.png",
        ancestral_population_size=ancestral_population_size,
        present_population_size=present_population_size,
        size_change_generations_ago=size_change_generations_ago,
        variant_age_generations=variant_age_generations,
        generation_time_years=generation_time_years,
        selection_coefficient=selection_coefficient,
    )
    exceedances, selected_pvalue = _plot_null_calibration(
        neutral,
        float(accepted_row["center_fraction_recent"]),
        theoretical_fraction_recent=theoretical_recent,
        threshold_generations=variant_age_generations,
        generation_time_years=generation_time_years,
        selection_coefficient=selection_coefficient,
        output_path=output_dir / "neutral_null_and_selected_pvalue.png",
    )
    plot_carrier_profile_figure(
        accepted_profile,
        replicate=int(accepted_row["attempt"]),
        population_allele_frequency=float(
            accepted_row["population_allele_frequency"]
        ),
        sequence_length=sequence_length,
        center=center,
        zoom_half_width=zoom_half_width,
        output_path=output_dir / "selected_hom_alt_vs_hom_ref_tmrca.png",
    )
    _plot_allele_frequency_trajectory(
        accepted_trajectory,
        population_allele_frequency=float(
            accepted_row["population_allele_frequency"]
        ),
        sample_allele_frequency=float(accepted_row["sample_allele_frequency"]),
        variant_age_generations=variant_age_generations,
        size_change_generations_ago=size_change_generations_ago,
        generation_time_years=generation_time_years,
        selection_coefficient=selection_coefficient,
        output_path=output_dir / "selected_allele_frequency_trajectory.png",
    )

    null_values = neutral["center_fraction_recent"].to_numpy(dtype=float)
    loo_pvalues = np.asarray([
        monte_carlo_pvalue(value, np.delete(null_values, index))
        for index, value in enumerate(null_values)
    ])
    tie_rng = np.random.default_rng(seed + 99)
    randomized_loo_pvalues = np.asarray([
        randomized_rank_pvalue(
            value, np.delete(null_values, index), tie_rng.random()
        )
        for index, value in enumerate(null_values)
    ])
    null_error = null_values - theoretical_recent
    result = {
        "demography": {
            "ancestral_population_size": int(ancestral_population_size),
            "present_population_size": int(present_population_size),
            "size_change_generations_ago": int(size_change_generations_ago),
            "interpretation": "forward growth from 10000 to 20000; backward-time size is 20000 from 0-100 generations and 10000 earlier",
        },
        "variant_age_generations": int(variant_age_generations),
        "generation_time_years": float(generation_time_years),
        "variant_age_years": float(
            variant_age_generations * generation_time_years
        ),
        "selection_coefficient": float(selection_coefficient),
        "sample_diploids": int(sample_diploids),
        "sequence_length": int(sequence_length),
        "mutation_rate": float(mutation_rate),
        "recombination_rate": float(recombination_rate),
        "neutral_replicates": int(neutral_replicates),
        "neutral_theoretical_fraction_recent": theoretical_recent,
        "neutral_observed_mean_fraction_recent": float(null_values.mean()),
        "neutral_fraction_recent_bias": float(null_error.mean()),
        "neutral_fraction_recent_rmse": float(np.sqrt(np.mean(null_error**2))),
        "neutral_leave_one_out_calibration": {
            **calibration_metrics(loo_pvalues),
            "fraction_p_le_0_05": float(np.mean(loo_pvalues <= 0.05)),
        },
        "neutral_randomized_tie_calibration_diagnostic": {
            **calibration_metrics(randomized_loo_pvalues),
            "fraction_p_le_0_05": float(
                np.mean(randomized_loo_pvalues <= 0.05)
            ),
        },
        "selected_rejection": {
            "acceptance_rule": f"population allele retained and at least {minimum_hom_alt_pairs} sampled hom-alt plus 2 hom-ref diploid pairs",
            "accepted_attempt_zero_based": int(accepted_row["attempt"]),
            "attempts_to_accept_in_seed_order": int(accepted_row["attempt"] + 1),
            "trajectories_computed_in_parallel_batches": int(len(attempts_table)),
            "population_allele_frequency": float(
                accepted_row["population_allele_frequency"]
            ),
            "sample_allele_frequency": float(
                accepted_row["sample_allele_frequency"]
            ),
            "trajectory_file": "selected_allele_frequency_trajectory.tsv",
            "trajectory_generations_logged": int(len(accepted_trajectory)),
            "initial_population_allele_frequency": float(
                accepted_trajectory.iloc[0]["population_allele_frequency"]
            ),
            "n_hom_ref_pairs": int(accepted_row["n_hom_ref_pairs"]),
            "n_heterozygous_pairs": int(
                accepted_row["n_heterozygous_pairs"]
            ),
            "n_hom_alt_pairs": int(accepted_row["n_hom_alt_pairs"]),
            "observed_fraction_recent": float(
                accepted_row["center_fraction_recent"]
            ),
            "neutral_exceedances": int(exceedances),
            "monte_carlo_p_upper": float(selected_pvalue),
        },
        "confidence_interval": "two-sided 95% log-Wald interval for the positive mean across diploid pairs, using a Student-t critical value",
        "workers_requested": int(workers),
        "workers_used": int(min(workers, neutral_replicates)),
        "base_seed": int(seed),
        "elapsed_seconds": float(perf_counter() - started),
    }
    center_profile = accepted_profile[
        np.isclose(accepted_profile["position_0based"], center)
    ]
    result["selected_center_carrier_tmrca"] = center_profile[
        [
            "pair_class",
            "n_pairs",
            "mean_tmrca_generations",
            "ci95_lower_generations",
            "ci95_upper_generations",
        ]
    ].to_dict("records")
    with (output_dir / "metrics.json").open("w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2)
    return result
