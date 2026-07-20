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

from .selection import (
    retained_sweep_calibration_table,
    run_slim_recent_sweep,
    within_individual_tmrca_details,
)


PAIR_STYLES = {
    0: ("hom ref", "#595959", "#bdbdbd"),
    2: ("hom alt", "#0072b2", "#56b4e9"),
}

PLOT_FONTS = {
    "suptitle": 22,
    "title": 20,
    "axis": 18,
    "tick": 15,
    "legend": 15,
}


def pair_tmrca_profile_by_focal_copy(
    ts,
    positions: np.ndarray,
    focal_carrier_counts: np.ndarray,
    *,
    confidence_level: float = 0.95,
) -> pd.DataFrame:
    """Mean and positive log-Wald CI for hom-ref and hom-alt diploid pairs."""
    from scipy.stats import t

    positions = np.asarray(positions, dtype=float)
    focal_carrier_counts = np.asarray(focal_carrier_counts, dtype=np.int8)
    if positions.ndim != 1 or not len(positions):
        raise ValueError("positions must be a nonempty one-dimensional array")
    if len(focal_carrier_counts) != ts.num_individuals:
        raise ValueError("focal carrier counts must have one value per individual")
    if not 0 < confidence_level < 1:
        raise ValueError("confidence_level must be between zero and one")

    sample_nodes = set(ts.samples())
    pairs_by_copy: dict[int, list[tuple[int, int]]] = {0: [], 2: []}
    for individual in ts.individuals():
        nodes = tuple(node for node in individual.nodes if node in sample_nodes)
        copies = int(focal_carrier_counts[individual.id])
        if len(nodes) == 2 and copies in pairs_by_copy:
            pairs_by_copy[copies].append(nodes)
    for copies, pairs in pairs_by_copy.items():
        if not pairs:
            raise ValueError(f"no sampled diploid pairs have {copies} focal copies")

    alpha = 1 - confidence_level
    rows = []
    for position in np.unique(positions):
        tree = ts.at(float(position))
        for copies, pairs in pairs_by_copy.items():
            times = np.fromiter(
                (tree.tmrca(node_a, node_b) for node_a, node_b in pairs),
                dtype=float,
                count=len(pairs),
            )
            mean = float(times.mean())
            if len(times) > 1:
                standard_error = float(times.std(ddof=1) / np.sqrt(len(times)))
                log_half_width = float(
                    t.ppf(1 - alpha / 2, df=len(times) - 1)
                    * standard_error
                    / mean
                )
            else:
                standard_error = np.nan
                log_half_width = np.nan
            rows.append({
                "position_0based": float(position),
                "focal_carrier_copies": copies,
                "pair_class": "hom_ref" if copies == 0 else "hom_alt",
                "n_pairs": int(len(times)),
                "mean_tmrca_generations": mean,
                "standard_error_generations": standard_error,
                "ci95_lower_generations": mean * np.exp(-log_half_width),
                "ci95_upper_generations": mean * np.exp(log_half_width),
            })
    return pd.DataFrame(rows)


def _position_grids(
    sequence_length: int,
    center: int,
    *,
    full_step: int,
    zoom_half_width: int,
    zoom_step: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if min(full_step, zoom_half_width, zoom_step) < 1:
        raise ValueError("profile steps and zoom width must be positive")
    full = np.arange(0, sequence_length, full_step, dtype=float)
    full = np.unique(np.r_[full, center, sequence_length - 1])
    zoom_start = max(0, center - zoom_half_width)
    zoom_end = min(sequence_length - 1, center + zoom_half_width)
    zoom = np.arange(zoom_start, zoom_end + 1, zoom_step, dtype=float)
    zoom = np.unique(np.r_[zoom, center, zoom_end])
    return full, zoom, np.unique(np.r_[full, zoom])


def _draw_profile_panel(
    axis,
    profile: pd.DataFrame,
    *,
    center: int,
    xlim: tuple[float, float],
    title: str,
    show_ylabel: bool,
    font_scale: float = 1.0,
) -> None:
    for copies in (0, 2):
        label, line_color, fill_color = PAIR_STYLES[copies]
        group = profile[profile["focal_carrier_copies"] == copies].sort_values(
            "position_0based"
        )
        x = group["position_0based"].to_numpy(dtype=float) / 1e6
        mean = group["mean_tmrca_generations"].to_numpy(dtype=float)
        lower = group["ci95_lower_generations"].to_numpy(dtype=float)
        upper = group["ci95_upper_generations"].to_numpy(dtype=float)
        axis.fill_between(x, lower, upper, color=fill_color, alpha=0.28, linewidth=0)
        axis.plot(x, mean, color=line_color, lw=1.35, label=f"{label} mean")
    axis.axvline(center / 1e6, color="#e69f00", lw=1.0, ls="--", label="selected site")
    axis.set_yscale("log")
    axis.set_xlim(xlim[0] / 1e6, xlim[1] / 1e6)
    axis.set_xlabel("Position (Mb)", fontsize=PLOT_FONTS["axis"] * font_scale)
    if show_ylabel:
        axis.set_ylabel(
            "Mean pairwise TMRCA (generations)",
            fontsize=PLOT_FONTS["axis"] * font_scale,
        )
    axis.set_title(title, fontsize=PLOT_FONTS["title"] * font_scale)
    axis.tick_params(axis="both", labelsize=PLOT_FONTS["tick"] * font_scale)
    axis.grid(alpha=0.16, linewidth=0.6)


def plot_carrier_profile_figure(
    profile: pd.DataFrame,
    *,
    replicate: int,
    population_allele_frequency: float,
    sequence_length: int,
    center: int,
    zoom_half_width: int,
    output_path: str | Path,
) -> None:
    """Draw a two-panel full-region and selected-site zoom profile."""
    output_path = Path(output_path)
    fig, axes = plt.subplots(1, 2, figsize=(19, 7.2), sharey=True, constrained_layout=True)
    _draw_profile_panel(
        axes[0],
        profile,
        center=center,
        xlim=(0, sequence_length),
        title="Full 10 Mb region",
        show_ylabel=True,
    )
    zoom_start = max(0, center - zoom_half_width)
    zoom_end = min(sequence_length, center + zoom_half_width)
    zoom = profile[
        profile["position_0based"].between(zoom_start, zoom_end, inclusive="both")
    ]
    _draw_profile_panel(
        axes[1],
        zoom,
        center=center,
        xlim=(zoom_start, zoom_end),
        title=f"Selected site +/-{zoom_half_width / 1e3:g} kb",
        show_ylabel=False,
    )
    n_ref = int(
        profile.loc[profile["focal_carrier_copies"] == 0, "n_pairs"].iloc[0]
    )
    n_alt = int(
        profile.loc[profile["focal_carrier_copies"] == 2, "n_pairs"].iloc[0]
    )
    axes[0].legend(loc="best", fontsize=PLOT_FONTS["legend"])
    fig.suptitle(
        f"Selected replicate {replicate}: hom-alt versus hom-ref pairwise TMRCA "
        f"(population AF={population_allele_frequency:.3f}; "
        f"n alt={n_alt}, n ref={n_ref}; mean +/- 95% CI)",
        fontsize=PLOT_FONTS["suptitle"],
    )
    fig.savefig(output_path, dpi=190)
    plt.close(fig)


def _plot_profile_overview(
    profiles: pd.DataFrame,
    selected: pd.DataFrame,
    *,
    sequence_length: int,
    center: int,
    zoom_half_width: int,
    output_path: Path,
) -> None:
    fig, axes = plt.subplots(5, 4, figsize=(26, 30), sharey=True, constrained_layout=True)
    for index, selected_row in enumerate(selected.itertuples(index=False)):
        row = index // 2
        column = (index % 2) * 2
        profile = profiles[profiles["replicate"] == selected_row.replicate]
        _draw_profile_panel(
            axes[row, column],
            profile,
            center=center,
            xlim=(0, sequence_length),
            title=f"rep {selected_row.replicate}: full 10 Mb",
            show_ylabel=column == 0,
            font_scale=0.72,
        )
        zoom_start = max(0, center - zoom_half_width)
        zoom_end = min(sequence_length, center + zoom_half_width)
        zoom = profile[
            profile["position_0based"].between(zoom_start, zoom_end, inclusive="both")
        ]
        _draw_profile_panel(
            axes[row, column + 1],
            zoom,
            center=center,
            xlim=(zoom_start, zoom_end),
            title=f"rep {selected_row.replicate}: +/-{zoom_half_width / 1e3:g} kb",
            show_ylabel=False,
            font_scale=0.72,
        )
    axes[0, 0].legend(loc="best", fontsize=PLOT_FONTS["legend"] * 0.72)
    fig.suptitle(
        "Hom-alt versus hom-ref within-diploid TMRCA: mean and 95% CI",
        fontsize=PLOT_FONTS["suptitle"],
    )
    fig.savefig(output_path, dpi=170)
    plt.close(fig)


def plot_retained_carrier_tmrca_profiles(
    source_dir: str | Path,
    output_dir: str | Path,
    *,
    executable: str | Path | None = None,
    selection_coefficient: float = 0.1,
    retained_replicates: int = 10,
    workers: int = 20,
    seed: int = 424242,
    full_step: int = 50_000,
    zoom_half_width: int = 500_000,
    zoom_step: int = 5_000,
) -> dict:
    """Reconstruct retained sweeps and plot hom-alt/hom-ref TMRCA profiles."""
    if workers < 1 or retained_replicates < 1:
        raise ValueError("workers and retained_replicates must be positive")
    source_dir = Path(source_dir).resolve()
    output_dir = Path(output_dir).resolve()
    figures_dir = output_dir / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    with (source_dir / "metrics.json").open(encoding="utf-8") as handle:
        source_metrics = json.load(handle)
    if source_metrics.get("seed") is not None and int(source_metrics["seed"]) != seed:
        raise ValueError(
            f"source seed is {source_metrics['seed']}, but reconstruction seed is {seed}"
        )
    stats = pd.read_csv(source_dir / "replicate_statistics.tsv", sep="\t")
    sample_diploids = int(source_metrics["sample_diploids"])
    _, selected = retained_sweep_calibration_table(
        stats,
        selection_coefficient=selection_coefficient,
        retained_replicates=retained_replicates,
        sample_diploids=sample_diploids,
    )
    selected_by_replicate = selected.set_index("replicate")
    sequence_length = int(source_metrics["sequence_length"])
    center = int(source_metrics["sweep_position"])
    full_positions, zoom_positions, positions = _position_grids(
        sequence_length,
        center,
        full_step=full_step,
        zoom_half_width=zoom_half_width,
        zoom_step=zoom_step,
    )
    started = perf_counter()

    with TemporaryDirectory(prefix="gamma_smc_carrier_profile_") as temporary_directory:
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
                sequence_length=sequence_length,
                sweep_position=center,
                selection_coefficient=selection_coefficient,
                age_generations=int(source_metrics["age_generations"]),
                mutation_rate=float(source_metrics["mutation_rate"]),
                recombination_rate=float(source_metrics["recombination_rate"]),
                seed=run_seed,
                capture_focal_genotypes=True,
            )
            if not np.isclose(
                run["realized_population_allele_frequency"],
                expected["realized_population_allele_frequency"],
                atol=1e-9,
                rtol=0,
            ):
                raise RuntimeError(f"replicate {replicate} allele frequency changed")
            center_details = within_individual_tmrca_details(
                ts,
                center,
                int(source_metrics["age_generations"]),
                focal_carrier_counts=run["focal_carrier_counts"],
            )
            if not np.isclose(
                center_details["tmrca_lt_threshold"].mean(),
                expected["center_fraction_recent"],
                atol=1e-12,
                rtol=0,
            ):
                raise RuntimeError(f"replicate {replicate} center statistic changed")
            profile = pair_tmrca_profile_by_focal_copy(
                ts, positions, run["focal_carrier_counts"]
            )
            profile["replicate"] = int(replicate)
            profile["population_focal_allele_frequency"] = float(
                run["realized_population_allele_frequency"]
            )
            profile["sample_focal_allele_frequency"] = float(
                run["sample_focal_allele_frequency"]
            )
            return profile

        replicate_ids = selected["replicate"].astype(int).tolist()
        if workers == 1:
            profile_tables = [reconstruct(replicate) for replicate in replicate_ids]
        else:
            with ThreadPoolExecutor(max_workers=workers) as executor:
                profile_tables = list(executor.map(reconstruct, replicate_ids))

    profiles = pd.concat(profile_tables, ignore_index=True)
    profiles.to_csv(output_dir / "hom_alt_vs_hom_ref_tmrca_profiles.tsv", sep="\t", index=False)
    for selected_row in selected.itertuples(index=False):
        profile = profiles[profiles["replicate"] == selected_row.replicate]
        plot_carrier_profile_figure(
            profile,
            replicate=int(selected_row.replicate),
            population_allele_frequency=float(
                selected_row.realized_population_allele_frequency
            ),
            sequence_length=sequence_length,
            center=center,
            zoom_half_width=zoom_half_width,
            output_path=figures_dir
            / f"replicate_{int(selected_row.replicate):03d}_hom_alt_vs_hom_ref_tmrca.png",
        )
    _plot_profile_overview(
        profiles,
        selected,
        sequence_length=sequence_length,
        center=center,
        zoom_half_width=zoom_half_width,
        output_path=output_dir / "hom_alt_vs_hom_ref_tmrca_overview.png",
    )
    result = {
        "source_dir": str(source_dir),
        "selection_coefficient": float(selection_coefficient),
        "retained_replicates": replicate_ids,
        "sample_diploids": sample_diploids,
        "sequence_length": sequence_length,
        "sweep_position": center,
        "full_step_bp": int(full_step),
        "zoom_half_width_bp": int(zoom_half_width),
        "zoom_step_bp": int(zoom_step),
        "full_positions": int(len(full_positions)),
        "zoom_positions": int(len(zoom_positions)),
        "unique_positions_evaluated": int(len(positions)),
        "confidence_interval": "two-sided 95% log-Wald interval for the positive mean across diploid pairs, using a Student-t critical value",
        "workers_requested": int(workers),
        "workers_used": int(min(workers, retained_replicates)),
        "base_seed": int(seed),
        "deterministic_reconstruction_verified": True,
        "elapsed_seconds": float(perf_counter() - started),
    }
    with (output_dir / "metrics.json").open("w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2)
    return result
