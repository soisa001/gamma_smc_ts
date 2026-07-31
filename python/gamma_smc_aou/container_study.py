from __future__ import annotations

import json
import os
import shutil
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path
from time import perf_counter

import numpy as np
import pandas as pd
import tskit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .container_decoder import DEFAULT_IMAGE, run_container_decoder
from .decoder import run_within_decoder
from .defaults import DEFAULT_CACHE_SIZE, DEFAULT_OUTPUT_STRIDE
from .selection import run_slim_recent_sweep, within_individual_tmrca_grid
from .spatial_scan import (
    _plot_null_spatial_calibration,
    _plot_observed_profile,
    calibrate_spatial_windows,
    significant_regions,
)
from .tree_sequence import tree_sequence_to_vcf


def scan_density_calibration(
    scan: pd.DataFrame,
    neutral: pd.DataFrame,
    *,
    center: int,
    local_half_width: int,
    output_dir: Path,
    statistic_column: str = "mean_p_tmrca_lt_threshold",
    statistic_label: str = "recent-coalescence statistic",
) -> dict:
    """Calibrate the number of pointwise-positive windows against null scans."""
    positions = scan["position_0based"].to_numpy(dtype=float)
    pivot = neutral[neutral["position_0based"].isin(positions)].pivot(
        index="replicate",
        columns="position_0based",
        values=statistic_column,
    ).reindex(columns=positions)
    if pivot.isna().any().any():
        raise ValueError("scan-density calibration requires a complete null grid")
    n_null = len(pivot)
    # Leave-one-out upper-tail p-values. With minimum ranks, ties are conservative.
    ranks = pivot.rank(axis=0, method="min").to_numpy(dtype=float)
    null_pvalues = (1 + n_null - ranks) / n_null
    local = np.abs(positions - center) <= local_half_width
    null_global_counts = np.count_nonzero(null_pvalues < 0.05, axis=1)
    null_local_counts = np.count_nonzero(null_pvalues[:, local] < 0.05, axis=1)
    observed_global = int(np.count_nonzero(scan["p_upper"].to_numpy() < 0.05))
    observed_local = int(
        np.count_nonzero(scan.loc[local, "p_upper"].to_numpy() < 0.05)
    )
    global_p = float(
        (1 + np.count_nonzero(null_global_counts >= observed_global))
        / (1 + n_null)
    )
    local_p = float(
        (1 + np.count_nonzero(null_local_counts >= observed_local))
        / (1 + n_null)
    )
    table = pd.DataFrame({
        "replicate": pivot.index.to_numpy(dtype=int),
        "global_n_pointwise_p_lt_0_05": null_global_counts,
        "local_n_pointwise_p_lt_0_05": null_local_counts,
    })
    table.to_csv(
        output_dir / "neutral_leave_one_out_scan_density.tsv", sep="\t", index=False
    )
    fig, axes = plt.subplots(1, 2, figsize=(18, 7), constrained_layout=True)
    panels = [
        (axes[0], null_global_counts, observed_global, global_p, "Full 10 Mb region"),
        (
            axes[1],
            null_local_counts,
            observed_local,
            local_p,
            f"Selected-site +/-{local_half_width / 1e3:g} kb",
        ),
    ]
    for axis, null_counts, observed_count, pvalue, title in panels:
        axis.hist(null_counts, bins=18, color="0.55", edgecolor="white")
        axis.axvline(
            observed_count,
            color="#7c3aed",
            lw=2.6,
            label=f"selected pseudo-data: {observed_count:,} windows",
        )
        axis.set_title(title, fontsize=23)
        axis.set_xlabel("Windows with pointwise p<0.05", fontsize=20)
        axis.set_ylabel("Neutral leave-one-out scans", fontsize=20)
        axis.tick_params(axis="both", labelsize=17)
        axis.annotate(
            f"selected: {observed_count:,} windows",
            xy=(observed_count, axis.get_ylim()[1]),
            xytext=(8, -8),
            textcoords="offset points",
            ha="left",
            va="top",
            fontsize=16,
            color="#7c3aed",
        )
        axis.text(
            0.97,
            0.78,
            f"regional-density p={pvalue:.4f}\n"
            f"null median={np.median(null_counts):,.0f}",
            transform=axis.transAxes,
            ha="right",
            va="top",
            fontsize=18,
        )
        axis.grid(axis="y", alpha=0.16)
    fig.suptitle(
        f"Simulation-calibrated density of {statistic_label} windows",
        fontsize=25,
    )
    fig.savefig(
        output_dir / "decoded_scan_density_calibration.png",
        dpi=190,
        bbox_inches="tight",
    )
    plt.close(fig)
    return {
        "global": {
            "n_windows": int(len(positions)),
            "selected_positive_windows": observed_global,
            "null_mean_positive_windows": float(null_global_counts.mean()),
            "null_median_positive_windows": float(np.median(null_global_counts)),
            "monte_carlo_p_upper": global_p,
        },
        "local": {
            "half_width_bp": int(local_half_width),
            "n_windows": int(np.count_nonzero(local)),
            "selected_positive_windows": observed_local,
            "null_mean_positive_windows": float(null_local_counts.mean()),
            "null_median_positive_windows": float(np.median(null_local_counts)),
            "monte_carlo_p_upper": local_p,
        },
    }


def _plot_decoded_center_calibration(
    neutral: pd.DataFrame,
    *,
    center_position: float,
    observed_value: float,
    pvalue: float,
    output_path: Path,
    statistic_column: str = "mean_p_tmrca_lt_threshold",
    statistic_label: str = "Mean inferred P(TMRCA < 4,500 years) across 2,000 pairs",
) -> None:
    center_rows = neutral[np.isclose(neutral["position_0based"], center_position)]
    values = center_rows[statistic_column].to_numpy(dtype=float)
    exceedances = int(np.count_nonzero(values >= observed_value))
    fig, axis = plt.subplots(figsize=(14, 7.5), constrained_layout=True)
    axis.hist(values, bins=16, color="0.55", edgecolor="white", alpha=0.9)
    axis.axvline(
        observed_value,
        color="#7c3aed",
        lw=2.8,
        label=f"selected Gamma-SMC estimate={observed_value:.4f}",
    )
    axis.set_title("Gamma-SMC center statistic versus decoded neutral null", fontsize=24)
    axis.set_xlabel(statistic_label, fontsize=20)
    axis.set_ylabel("Decoded neutral simulations", fontsize=20)
    axis.tick_params(axis="both", labelsize=17)
    if np.ptp(values) == 0:
        value = float(values[0])
        pad = max(0.005, abs(value) * 0.5)
        axis.set_xlim(value - pad, value + pad)
        axis.text(
            0.03,
            0.93,
            f"all {len(values)} neutral simulations = {value:g}",
            transform=axis.transAxes,
            ha="left",
            va="top",
            fontsize=18,
        )
    axis.text(
        0.97,
        0.93,
        f"{exceedances}/{len(values)} neutral >= selected\nMonte Carlo p={pvalue:.4f}",
        transform=axis.transAxes,
        ha="right",
        va="top",
        fontsize=18,
    )
    handles, labels = axis.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper left",
        bbox_to_anchor=(1.005, 0.88),
        fontsize=16,
    )
    axis.grid(axis="y", alpha=0.16)
    fig.savefig(output_path, dpi=190, bbox_inches="tight")
    plt.close(fig)


def _threshold_suffix(years: float) -> str:
    return str(int(years)) if float(years).is_integer() else str(years)


def _calibration_column(calibration_statistic: str, threshold_years: float) -> str:
    if calibration_statistic == "mean-posterior-probability":
        return "mean_p_tmrca_lt_threshold"
    if calibration_statistic == "called-fraction":
        return f"frac_recent_{_threshold_suffix(threshold_years)}"
    raise ValueError(
        "calibration_statistic must be mean-posterior-probability or called-fraction"
    )


def _statistic_label(
    calibration_statistic: str,
    recent_call: str,
) -> str:
    if calibration_statistic == "mean-posterior-probability":
        return "Mean posterior P"
    return f"Fraction called recent by posterior {recent_call}"


def _assert_matched_soft_profiles(
    first: pd.DataFrame,
    second: pd.DataFrame,
) -> None:
    """Ensure two call-rule decodes used identical data and posteriors."""
    columns = [
        "position_0based",
        "position_1based",
        "n_pairs",
        "mean_p_tmrca_lt_threshold",
        "mean_tmrca_generations",
    ]
    probability_columns = sorted(
        column for column in first if column.startswith("mean_p_lt_")
    )
    columns.extend(probability_columns)
    if list(first.columns) != list(second.columns) or len(first) != len(second):
        raise RuntimeError("recent-call comparison profiles have different schemas")
    for column in columns:
        left = first[column].to_numpy()
        right = second[column].to_numpy()
        if not np.allclose(left, right, rtol=0, atol=1e-12, equal_nan=True):
            raise RuntimeError(
                f"recent-call comparison changed posterior output column {column}"
            )


def _write_recent_call_comparison(
    observed_by_call: dict[str, pd.DataFrame],
    neutral_by_call: dict[str, pd.DataFrame],
    truth: pd.DataFrame,
    *,
    threshold_years: float,
    center: int,
    sequence_length: int,
    stride: int,
    output_dir: Path,
) -> dict:
    """Compare soft probability, posterior-mean calls, and median calls."""
    required_calls = {"mean", "median"}
    if required_calls.difference(observed_by_call):
        raise ValueError("comparison requires both mean and median selected profiles")
    if required_calls.difference(neutral_by_call):
        raise ValueError("comparison requires both mean and median neutral profiles")

    called_column = f"frac_recent_{_threshold_suffix(threshold_years)}"
    definitions = [
        (
            "soft_probability",
            "mean",
            "mean_p_tmrca_lt_threshold",
            "Mean posterior P(TMRCA < 4,500 years)",
        ),
        (
            "posterior_mean_call",
            "mean",
            called_column,
            "Fraction: posterior mean TMRCA < 4,500 years",
        ),
        (
            "posterior_median_call",
            "median",
            called_column,
            "Fraction: posterior median TMRCA < 4,500 years",
        ),
    ]
    truth_profile = truth.set_index("position_0based")["truth_fraction_recent"]
    metrics: dict[str, dict] = {}
    center_rows = []
    plot_rows = []
    for slug, call, statistic_column, label in definitions:
        observed = observed_by_call[call]
        neutral = neutral_by_call[call]
        scan = calibrate_spatial_windows(
            observed,
            neutral,
            statistic_column=statistic_column,
        )
        regions = significant_regions(
            scan,
            sequence_length=sequence_length,
            window_size=stride,
        )
        scan.to_csv(
            output_dir / f"comparison_{slug}_spatial_calibration.tsv",
            sep="\t",
            index=False,
        )
        regions.to_csv(
            output_dir / f"comparison_{slug}_significant_regions.tsv",
            sep="\t",
            index=False,
        )
        center_row = scan.iloc[
            np.argmin(
                np.abs(scan["position_0based"].to_numpy(dtype=float) - center)
            )
        ]
        aligned_truth = truth_profile.reindex(
            observed["position_0based"].to_numpy(dtype=float)
        ).to_numpy(dtype=float)
        estimate = observed[statistic_column].to_numpy(dtype=float)
        error = estimate - aligned_truth
        correlation = (
            float(np.corrcoef(estimate, aligned_truth)[0, 1])
            if np.std(estimate) > 0 and np.std(aligned_truth) > 0
            else None
        )
        item = {
            "recent_call": call if slug != "soft_probability" else None,
            "statistic_column": statistic_column,
            "statistic_label": label,
            "center_position_0based": int(center_row["position_0based"]),
            "center_selected": float(center_row["observed_fraction_recent"]),
            "center_null_mean": float(center_row["neutral_mean_fraction_recent"]),
            "center_null_median": float(
                center_row["neutral_median_fraction_recent"]
            ),
            "center_null_ci95_lower": float(center_row["neutral_ci95_lower"]),
            "center_null_ci95_upper": float(center_row["neutral_ci95_upper"]),
            "center_neutral_exceedances": int(
                center_row["neutral_exceedances"]
            ),
            "center_monte_carlo_p_upper": float(center_row["p_upper"]),
            "profile_truth_bias": float(np.mean(error)),
            "profile_truth_mae": float(np.mean(np.abs(error))),
            "profile_truth_rmse": float(np.sqrt(np.mean(error**2))),
            "profile_truth_correlation": correlation,
            "n_pointwise_p_lt_0_05_windows": int(
                np.count_nonzero(scan["p_upper"] < 0.05)
            ),
            "n_bh_q_lt_0_05_windows": int(
                np.count_nonzero(scan["q_bh"] < 0.05)
            ),
            "n_significant_regions": int(len(regions)),
        }
        metrics[slug] = item
        center_rows.append({"statistic": slug, **item})
        plot_rows.append((slug, label, observed, neutral, scan))

    truth_center_position = float(
        truth.iloc[
            np.argmin(
                np.abs(truth["position_0based"].to_numpy(dtype=float) - center)
            )
        ]["position_0based"]
    )
    truth_center = float(
        truth.loc[
            np.isclose(truth["position_0based"], truth_center_position),
            "truth_fraction_recent",
        ].iloc[0]
    )
    metrics["truth"] = {
        "center_position_0based": int(truth_center_position),
        "center_fraction_recent": truth_center,
    }
    pd.DataFrame(center_rows).to_csv(
        output_dir / "posterior_summary_rule_center_comparison.tsv",
        sep="\t",
        index=False,
    )

    fig, axes = plt.subplots(
        len(plot_rows),
        2,
        figsize=(22, 16),
        gridspec_kw={"width_ratios": [1.6, 1]},
        constrained_layout=True,
    )
    for row_index, (slug, label, _observed, neutral, scan) in enumerate(plot_rows):
        position = scan["position_0based"].to_numpy(dtype=float)
        left = axes[row_index, 0]
        left.fill_between(
            position / 1e6,
            scan["neutral_ci95_lower"].to_numpy(dtype=float),
            scan["neutral_ci95_upper"].to_numpy(dtype=float),
            color="0.78",
            alpha=0.55,
            linewidth=0,
            label="neutral 95% interval",
        )
        left.plot(
            position / 1e6,
            scan["neutral_median_fraction_recent"],
            color="0.35",
            lw=1.2,
            label="neutral median",
        )
        left.plot(
            position / 1e6,
            scan["observed_fraction_recent"],
            color="#7c3aed",
            lw=1.6,
            label="selected pseudo-data",
        )
        if slug != "soft_probability":
            left.plot(
                truth["position_0based"].to_numpy(dtype=float) / 1e6,
                truth["truth_fraction_recent"].to_numpy(dtype=float),
                color="#009e73",
                lw=1.8,
                label="tree-sequence truth",
            )
            left.set_ylim(
                0,
                max(
                    0.01,
                    1.08
                    * float(
                        np.nanmax(
                            truth["truth_fraction_recent"].to_numpy(dtype=float)
                        )
                    ),
                ),
            )
        else:
            left.set_ylim(bottom=0)
        left.axvline(center / 1e6, color="#e69f00", ls="--", lw=1.4)
        left.set_xlim(0, sequence_length / 1e6)
        left.set_xlabel("Position (Mb)", fontsize=15)
        left.set_ylabel(label, fontsize=15)
        left.tick_params(axis="both", labelsize=12)
        left.grid(alpha=0.16)

        center_row = scan.iloc[
            np.argmin(
                np.abs(scan["position_0based"].to_numpy(dtype=float) - center)
            )
        ]
        center_position = float(center_row["position_0based"])
        # Pull the null values through the calibrated scan's source definition.
        definition = next(value for value in definitions if value[0] == slug)
        null_values = neutral.loc[
            np.isclose(neutral["position_0based"], center_position),
            definition[2],
        ].to_numpy(dtype=float)
        right = axes[row_index, 1]
        if np.ptp(null_values) == 0:
            value = float(null_values[0])
            width = max(0.001, abs(value) * 0.2)
            right.bar(
                [value],
                [len(null_values)],
                width=width,
                color="0.58",
                edgecolor="white",
            )
            pad = max(0.005, abs(value) * 0.5)
            right.set_xlim(value - pad, value + pad)
            right.text(
                0.03,
                0.92,
                f"all {len(null_values)} nulls = {value:g}",
                transform=right.transAxes,
                ha="left",
                va="top",
                fontsize=14,
            )
        else:
            right.hist(null_values, bins=16, color="0.58", edgecolor="white")
        right.axvline(
            float(center_row["observed_fraction_recent"]),
            color="#7c3aed",
            lw=2.5,
            label="selected pseudo-data",
        )
        right.set_xlabel(f"Center: {label}", fontsize=15)
        right.set_ylabel("Neutral simulations", fontsize=15)
        right.tick_params(axis="both", labelsize=12)
        right.grid(axis="y", alpha=0.16)
        right.text(
            0.97,
            0.92,
            f"p={float(center_row['p_upper']):.4f}\n"
            f"truth fraction={truth_center:.4f}",
            transform=right.transAxes,
            ha="right",
            va="top",
            fontsize=14,
        )
    handles, labels = [], []
    for axis in axes[:, 0]:
        for handle, label in zip(*axis.get_legend_handles_labels()):
            if label not in labels:
                handles.append(handle)
                labels.append(label)
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.995),
        ncol=len(labels),
        fontsize=13,
    )
    fig.suptitle(
        "Matched posterior-summary comparison on the same selected and null simulations",
        fontsize=23,
        y=1.02,
    )
    fig.savefig(
        output_dir / "posterior_mean_median_rule_comparison.png",
        dpi=190,
        bbox_inches="tight",
    )
    plt.close(fig)
    with (output_dir / "posterior_summary_rule_comparison_metrics.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(metrics, handle, indent=2)
    return metrics


def _load_design(source_dir: Path) -> dict:
    with (source_dir / "metrics.json").open(encoding="utf-8") as handle:
        metrics = json.load(handle)
    selected = pd.read_csv(source_dir / "selected_observation.tsv", sep="\t").iloc[0]
    selected_tree_file = str(
        metrics.get("selected_tree_file", "selected_s0p05_af30.trees")
    )
    if Path(selected_tree_file).name != selected_tree_file:
        raise ValueError("selected_tree_file must be a file name, not a path")
    return {
        "ancestral_size": int(metrics["demography"]["ancestral_population_size"]),
        "present_size": int(metrics["demography"]["present_population_size"]),
        "change_age": int(metrics["demography"]["size_change_generations_ago"]),
        "sample_diploids": int(metrics["sample_diploids"]),
        "sequence_length": int(metrics["sequence_length"]),
        "variant_age": int(metrics["variant_age_generations"]),
        "generation_time": float(metrics["generation_time_years"]),
        "selection_coefficient": float(metrics["selection_coefficient"]),
        "mutation_rate": float(metrics["mutation_rate"]),
        "recombination_rate": float(metrics["recombination_rate"]),
        "seed": int(metrics["base_seed"]),
        "selected_attempt": int(selected["attempt"]),
        "selected_population_af": float(selected["population_allele_frequency"]),
        "selected_tree_file": selected_tree_file,
    }


def _simulate(
    path: Path,
    design: dict,
    *,
    selection_coefficient: float,
    seed: int,
    capture_focal_genotypes: bool = False,
):
    return run_slim_recent_sweep(
        path,
        population_size=design["ancestral_size"],
        present_population_size=design["present_size"],
        size_change_generations_ago=design["change_age"],
        sample_diploids=design["sample_diploids"],
        sequence_length=design["sequence_length"],
        sweep_position=design["sequence_length"] // 2,
        selection_coefficient=selection_coefficient,
        age_generations=design["variant_age"],
        mutation_rate=design["mutation_rate"],
        recombination_rate=design["recombination_rate"],
        seed=seed,
        capture_focal_genotypes=capture_focal_genotypes,
    )


def _prepare_neutral_input(task: dict) -> dict:
    """Simulate and write one neutral VCF in an independent process."""
    replicate = int(task["replicate"])
    tree_path = Path(task["tree_path"])
    vcf_path = Path(task["vcf_path"])
    neutral_ts, _ = _simulate(
        tree_path,
        task["design"],
        selection_coefficient=0.0,
        seed=int(task["seed"]),
    )
    del neutral_ts
    tree_sequence_to_vcf(tree_path, vcf_path)
    return {
        "replicate": replicate,
        "tree_path": str(tree_path),
        "vcf_path": str(vcf_path),
    }


def _run_stride_study(
    source_dir: str | Path,
    output_dir: str | Path,
    *,
    neutral_replicates: int = 100,
    stride: int = DEFAULT_OUTPUT_STRIDE,
    runtime: str = "auto",
    image: str = DEFAULT_IMAGE,
    executable: str | Path | None = None,
    cache_size: int = DEFAULT_CACHE_SIZE,
    threads: int = 1,
    keep_vcfs: bool = False,
    workers: int = 1,
    recent_call: str = "median",
    comparison_recent_call: str | None = None,
    calibration_statistic: str = "mean-posterior-probability",
) -> dict:
    """Decode a retained sweep and matched nulls with one decoder backend."""
    if (
        neutral_replicates < 1
        or stride < 1
        or workers < 1
        or cache_size < 1
        or threads < 1
    ):
        raise ValueError(
            "replicate count, stride, workers, cache size, and threads must be positive"
        )
    if recent_call not in {"mean", "median"}:
        raise ValueError("study recent_call must be mean or median")
    if comparison_recent_call not in {None, "mean", "median"}:
        raise ValueError("comparison_recent_call must be mean, median, or omitted")
    if comparison_recent_call == recent_call:
        raise ValueError("comparison_recent_call must differ from recent_call")
    if calibration_statistic not in {
        "mean-posterior-probability",
        "called-fraction",
    }:
        raise ValueError(
            "calibration_statistic must be mean-posterior-probability or called-fraction"
        )
    if executable is None and (
        recent_call != "median"
        or comparison_recent_call is not None
        or calibration_statistic != "mean-posterior-probability"
    ):
        raise ValueError(
            "posterior mean/median call comparisons require the optimized native decoder"
        )
    source_dir = Path(source_dir).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    work_dir = output_dir / "work"
    work_dir.mkdir(exist_ok=True)
    design = _load_design(source_dir)
    center = design["sequence_length"] // 2
    theta = 4 * design["ancestral_size"] * design["mutation_rate"]
    rho_over_theta = design["recombination_rate"] / design["mutation_rate"]
    threshold_years = design["variant_age"] * design["generation_time"]
    statistic_column = _calibration_column(
        calibration_statistic,
        threshold_years,
    )
    call_rules = [recent_call]
    if comparison_recent_call is not None:
        call_rules.append(comparison_recent_call)
    started = perf_counter()
    backend = "container_v0.2" if executable is None else "native_optimized"

    def decode(vcf_path: Path, summary_path: Path, call_rule: str) -> dict:
        if executable is None:
            return run_container_decoder(
                vcf_path,
                summary_path,
                scaled_mutation_rate=theta,
                recombination_to_mutation_ratio=rho_over_theta,
                mutation_rate=design["mutation_rate"],
                threshold_years=threshold_years,
                generation_time=design["generation_time"],
                stride=stride,
                runtime=runtime,
                image=image,
            )
        return run_within_decoder(
            executable,
            vcf_path,
            summary_path,
            scaled_mutation_rate=theta,
            recombination_to_mutation_ratio=rho_over_theta,
            mutation_rate=design["mutation_rate"],
            threshold_years=threshold_years,
            generation_time=design["generation_time"],
            input_format="vcf",
            output_at_stride=stride,
            output_at_hets=False,
            only_within=True,
            recent_call=call_rule,
            cache_size=cache_size,
            threads=threads,
        )

    selected_tree_path = work_dir / "selected.trees"
    prepared_selected_tree = source_dir / design["selected_tree_file"]
    if prepared_selected_tree.exists():
        shutil.copy2(prepared_selected_tree, selected_tree_path)
        selected_ts = tskit.load(selected_tree_path)
        if not np.isclose(
            selected_ts.sequence_length,
            design["sequence_length"],
            atol=0,
            rtol=0,
        ):
            raise RuntimeError("prepared selected tree has the wrong sequence length")
        n_diploids = sum(
            len(individual.nodes) == 2 for individual in selected_ts.individuals()
        )
        if n_diploids != design["sample_diploids"]:
            raise RuntimeError(
                "prepared selected tree has the wrong number of diploid individuals"
            )
        selected_origin = "prepared_selected_tree"
    else:
        selected_ts, selected_run = _simulate(
            selected_tree_path,
            design,
            selection_coefficient=design["selection_coefficient"],
            seed=design["seed"] + 10_000_000 + design["selected_attempt"] * 10,
            capture_focal_genotypes=True,
        )
        if not np.isclose(
            selected_run["realized_population_allele_frequency"],
            design["selected_population_af"],
            atol=1e-12,
            rtol=0,
        ):
            raise RuntimeError("selected replicate did not reproduce the retained allele")
        selected_origin = "reproduced_from_seed"
    selected_vcf = work_dir / "selected.vcf.gz"
    tree_sequence_to_vcf(selected_tree_path, selected_vcf)
    selected_profiles: dict[str, pd.DataFrame] = {}
    selected_decodes: dict[str, dict] = {}
    for call_rule in call_rules:
        selected_summary = (
            work_dir / "selected.tsv"
            if call_rule == recent_call
            else work_dir / f"selected_{call_rule}_call.tsv"
        )
        selected_decodes[call_rule] = decode(
            selected_vcf,
            selected_summary,
            call_rule,
        )
        selected_profiles[call_rule] = pd.read_csv(selected_summary, sep="\t")
    observed = selected_profiles[recent_call]
    selected_decode = selected_decodes[recent_call]
    if comparison_recent_call is not None:
        _assert_matched_soft_profiles(
            observed,
            selected_profiles[comparison_recent_call],
        )
        for call_rule, profile in selected_profiles.items():
            profile.to_csv(
                output_dir / f"selected_decoded_{call_rule}_call_profile.tsv",
                sep="\t",
                index=False,
            )
    observed.to_csv(
        output_dir / "selected_decoded_recent_probability_profile.tsv",
        sep="\t",
        index=False,
    )
    truth = within_individual_tmrca_grid(
        selected_ts,
        observed["position_0based"].to_numpy(dtype=float),
        design["variant_age"],
    ).rename(columns={"mean_p_tmrca_lt_threshold": "truth_fraction_recent"})
    comparison = observed.merge(
        truth[["position_0based", "truth_fraction_recent", "mean_tmrca_generations"]],
        on="position_0based",
        how="left",
        suffixes=("_decoded", "_truth"),
        validate="one_to_one",
    )
    error = (
        comparison["mean_p_tmrca_lt_threshold"]
        - comparison["truth_fraction_recent"]
    )
    comparison.to_csv(
        output_dir / "selected_decoded_vs_truth.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    def neutral_task(
        prepared: dict,
    ) -> tuple[dict[str, pd.DataFrame], dict[str, dict]]:
        replicate = int(prepared["replicate"])
        tree_path = Path(prepared["tree_path"])
        vcf_path = Path(prepared["vcf_path"])
        profiles: dict[str, pd.DataFrame] = {}
        runs: dict[str, dict] = {}
        summary_paths: list[Path] = []
        for call_rule in call_rules:
            summary_path = (
                work_dir / f"neutral_{replicate:04d}.tsv"
                if call_rule == recent_call
                else work_dir / f"neutral_{replicate:04d}_{call_rule}_call.tsv"
            )
            summary_paths.append(summary_path)
            runs[call_rule] = {
                "replicate": replicate,
                **decode(vcf_path, summary_path, call_rule),
            }
            profile = pd.read_csv(summary_path, sep="\t")
            profile["replicate"] = replicate
            profiles[call_rule] = profile
        if comparison_recent_call is not None:
            _assert_matched_soft_profiles(
                profiles[recent_call],
                profiles[comparison_recent_call],
            )
        if not keep_vcfs:
            vcf_path.unlink(missing_ok=True)
            tree_path.unlink(missing_ok=True)
            for summary_path in summary_paths:
                summary_path.unlink(missing_ok=True)
                summary_path.with_name(summary_path.name + ".run.json").unlink(
                    missing_ok=True
                )
        return profiles, runs

    neutral_tasks = [
        {
            "replicate": replicate,
            "tree_path": str(work_dir / f"neutral_{replicate:04d}.trees"),
            "vcf_path": str(work_dir / f"neutral_{replicate:04d}.vcf.gz"),
            "design": design,
            "seed": design["seed"] + 1_000_000 + replicate * 10,
        }
        for replicate in range(neutral_replicates)
    ]
    if workers == 1:
        neutral_results = [
            neutral_task(_prepare_neutral_input(task)) for task in neutral_tasks
        ]
    else:
        worker_count = min(workers, neutral_replicates)
        neutral_results = []
        with (
            ProcessPoolExecutor(max_workers=worker_count) as simulation_executor,
            ThreadPoolExecutor(max_workers=worker_count) as decoder_executor,
        ):
            for batch_start in range(0, neutral_replicates, worker_count):
                batch = neutral_tasks[batch_start : batch_start + worker_count]
                prepared = list(
                    simulation_executor.map(_prepare_neutral_input, batch)
                )
                neutral_results.extend(
                    decoder_executor.map(neutral_task, prepared)
                )
    neutral_by_call = {
        call_rule: pd.concat(
            [result[0][call_rule] for result in neutral_results],
            ignore_index=True,
        )
        for call_rule in call_rules
    }
    neutral_runs_by_call = {
        call_rule: [result[1][call_rule] for result in neutral_results]
        for call_rule in call_rules
    }
    neutral = neutral_by_call[recent_call]
    neutral_runs = neutral_runs_by_call[recent_call]
    if comparison_recent_call is not None:
        _assert_matched_soft_profiles(
            neutral,
            neutral_by_call[comparison_recent_call],
        )
        for call_rule, profile in neutral_by_call.items():
            profile.to_csv(
                output_dir / f"neutral_decoded_{call_rule}_call_profiles.tsv.gz",
                sep="\t",
                index=False,
                compression="gzip",
            )
    neutral.to_csv(
        output_dir / "neutral_decoded_recent_probability_profiles.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    scan = calibrate_spatial_windows(
        observed,
        neutral,
        statistic_column=statistic_column,
    )
    regions = significant_regions(
        scan,
        sequence_length=design["sequence_length"],
        window_size=stride,
    )
    scan.to_csv(
        output_dir / "decoded_spatial_pointwise_calibration.tsv",
        sep="\t",
        index=False,
    )
    regions.to_csv(
        output_dir / "decoded_spatial_significant_regions.tsv",
        sep="\t",
        index=False,
    )
    _plot_observed_profile(
        scan,
        center=center,
        sequence_length=design["sequence_length"],
        zoom_half_width=500_000,
        threshold_years=threshold_years,
        window_size=stride,
        output_path=output_dir / "selected_decoded_recent_probability_spatial.png",
        pair_count=design["sample_diploids"],
        series_label=f"selected pseudo-data ({recent_call} call)",
        statistic_label=_statistic_label(calibration_statistic, recent_call),
        profile_title=(
            "Recent-coalescence call profile"
            if calibration_statistic == "called-fraction"
            else "Recent-coalescence probability profile"
        ),
    )
    _plot_null_spatial_calibration(
        scan,
        regions,
        center=center,
        sequence_length=design["sequence_length"],
        zoom_half_width=500_000,
        threshold_years=threshold_years,
        window_size=stride,
        output_path=output_dir / "selected_decoded_vs_neutral_pvalues.png",
        pair_count=design["sample_diploids"],
        series_label=f"selected pseudo-data ({recent_call} call)",
        statistic_label=_statistic_label(calibration_statistic, recent_call),
    )
    density_calibration = scan_density_calibration(
        scan,
        neutral,
        center=center,
        local_half_width=500_000,
        output_dir=output_dir,
        statistic_column=statistic_column,
        statistic_label=_statistic_label(calibration_statistic, recent_call),
    )
    center_row = scan.iloc[
        np.argmin(np.abs(scan["position_0based"].to_numpy(dtype=float) - center))
    ]
    _plot_decoded_center_calibration(
        neutral,
        center_position=float(center_row["position_0based"]),
        observed_value=float(center_row["observed_fraction_recent"]),
        pvalue=float(center_row["p_upper"]),
        output_path=output_dir / "gamma_smc_center_null_and_selected_pvalue.png",
        statistic_column=statistic_column,
        statistic_label=(
            f"{_statistic_label(calibration_statistic, recent_call)} across "
            f"{design['sample_diploids']:,} pairs"
        ),
    )
    comparison_metrics = None
    if comparison_recent_call is not None:
        comparison_metrics = _write_recent_call_comparison(
            selected_profiles,
            neutral_by_call,
            truth,
            threshold_years=threshold_years,
            center=center,
            sequence_length=design["sequence_length"],
            stride=stride,
            output_dir=output_dir,
        )
    decode_times = np.asarray(
        [
            run["decode_seconds"]
            for run in selected_decodes.values()
        ]
        + [
            run["decode_seconds"]
            for runs in neutral_runs_by_call.values()
            for run in runs
        ]
    )
    result = {
        "data_interpretation": "retained selected simulation treated as pseudo-empirical data",
        "selected_input_origin": selected_origin,
        "decoder_backend": backend,
        "container_image": image if executable is None else None,
        "container_runtime": selected_decode.get("runtime"),
        "native_executable": str(Path(executable).resolve()) if executable is not None else None,
        "cache_size_bp": int(cache_size) if executable is not None else None,
        "decoder_threads_per_replicate": int(threads) if executable is not None else None,
        "recent_call": recent_call,
        "comparison_recent_call": comparison_recent_call,
        "calibration_statistic": calibration_statistic,
        "calibration_statistic_column": statistic_column,
        "stride_bp": int(stride),
        "sequence_length": design["sequence_length"],
        "sample_diploids": design["sample_diploids"],
        "within_individual_pairs": design["sample_diploids"],
        "neutral_replicates": int(neutral_replicates),
        "neutral_workers": int(min(workers, neutral_replicates)),
        "neutral_preparation_executor": (
            "serial" if workers == 1 else "process_pool_batched"
        ),
        "neutral_decode_executor": (
            "serial" if workers == 1 else "thread_pool_batched"
        ),
        "neutral_batch_size": int(min(workers, neutral_replicates)),
        "observed_output_positions": int(len(observed)),
        "calibrated_complete_positions": int(len(scan)),
        "dropped_incomplete_null_positions": int(len(observed) - len(scan)),
        "scaled_mutation_rate": float(theta),
        "recombination_to_mutation_ratio": float(rho_over_theta),
        "threshold_years": float(threshold_years),
        "selected_decode_seconds": float(selected_decode["decode_seconds"]),
        "neutral_decode_seconds_mean": float(
            np.mean([run["decode_seconds"] for run in neutral_runs])
        ),
        "all_decode_seconds_sum": float(decode_times.sum()),
        "selected_profile_truth_bias": float(np.mean(error)),
        "selected_profile_truth_mae": float(np.mean(np.abs(error))),
        "selected_profile_truth_rmse": float(np.sqrt(np.mean(error**2))),
        "selected_profile_truth_correlation": float(
            np.corrcoef(
                comparison["mean_p_tmrca_lt_threshold"],
                comparison["truth_fraction_recent"],
            )[0, 1]
        ),
        "center_position_0based": int(center_row["position_0based"]),
        "center_observed_statistic": float(center_row["observed_fraction_recent"]),
        "center_null_mean_statistic": float(
            center_row["neutral_mean_fraction_recent"]
        ),
        "center_null_median_statistic": float(
            center_row["neutral_median_fraction_recent"]
        ),
        "center_null_ci95_lower": float(center_row["neutral_ci95_lower"]),
        "center_null_ci95_upper": float(center_row["neutral_ci95_upper"]),
        "center_neutral_exceedances": int(center_row["neutral_exceedances"]),
        "center_monte_carlo_p_upper": float(center_row["p_upper"]),
        "n_pointwise_p_lt_0_05_windows": int(np.count_nonzero(scan["p_upper"] < 0.05)),
        "n_bh_q_lt_0_05_windows": int(np.count_nonzero(scan["q_bh"] < 0.05)),
        "n_significant_regions": int(len(regions)),
        "significant_regions": regions.to_dict("records"),
        "scan_density_calibration": density_calibration,
        "posterior_summary_rule_comparison": comparison_metrics,
        "multiple_testing_note": "regions use raw pointwise p<0.05; BH q-values are in the calibration TSV",
        "elapsed_seconds": float(perf_counter() - started),
    }
    if calibration_statistic == "mean-posterior-probability":
        result.update({
            "center_observed_mean_p_recent": float(
                center_row["observed_fraction_recent"]
            ),
            "center_null_mean_p_recent": float(
                center_row["neutral_mean_fraction_recent"]
            ),
        })
    else:
        result.update({
            "center_observed_called_fraction_recent": float(
                center_row["observed_fraction_recent"]
            ),
            "center_null_mean_called_fraction_recent": float(
                center_row["neutral_mean_fraction_recent"]
            ),
        })
    if comparison_recent_call is not None:
        comparison_selected_decode = selected_decodes[comparison_recent_call]
        comparison_neutral_runs = neutral_runs_by_call[comparison_recent_call]
        result.update({
            "comparison_selected_decode_seconds": float(
                comparison_selected_decode["decode_seconds"]
            ),
            "comparison_neutral_decode_seconds_mean": float(
                np.mean(
                    [run["decode_seconds"] for run in comparison_neutral_runs]
                )
            ),
        })
    with (output_dir / "decoded_study_metrics.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(result, handle, indent=2)
    if not keep_vcfs:
        selected_vcf.unlink(missing_ok=True)
        selected_tree_path.unlink(missing_ok=True)
    return result


def run_container_stride_study(
    source_dir: str | Path,
    output_dir: str | Path,
    *,
    neutral_replicates: int = 100,
    stride: int = DEFAULT_OUTPUT_STRIDE,
    runtime: str = "auto",
    image: str = DEFAULT_IMAGE,
    keep_vcfs: bool = False,
    workers: int = 1,
) -> dict:
    """Decode a retained sweep and matched nulls with official Gamma-SMC v0.2."""
    return _run_stride_study(
        source_dir,
        output_dir,
        neutral_replicates=neutral_replicates,
        stride=stride,
        runtime=runtime,
        image=image,
        keep_vcfs=keep_vcfs,
        workers=workers,
    )


def run_native_stride_study(
    source_dir: str | Path,
    output_dir: str | Path,
    *,
    executable: str | Path,
    neutral_replicates: int = 100,
    stride: int = DEFAULT_OUTPUT_STRIDE,
    cache_size: int = DEFAULT_CACHE_SIZE,
    threads: int = 1,
    keep_vcfs: bool = False,
    workers: int = 12,
    recent_call: str = "median",
    comparison_recent_call: str | None = None,
    calibration_statistic: str = "mean-posterior-probability",
) -> dict:
    """Decode a retained sweep and matched nulls with the optimized binary."""
    return _run_stride_study(
        source_dir,
        output_dir,
        neutral_replicates=neutral_replicates,
        stride=stride,
        executable=executable,
        cache_size=cache_size,
        threads=threads,
        keep_vcfs=keep_vcfs,
        workers=workers,
        recent_call=recent_call,
        comparison_recent_call=comparison_recent_call,
        calibration_statistic=calibration_statistic,
    )


def finalize_container_stride_study(
    source_dir: str | Path,
    output_dir: str | Path,
    *,
    stride: int = DEFAULT_OUTPUT_STRIDE,
    workflow_elapsed_seconds: float | None = None,
    neutral_profiles_path: str | Path | None = None,
) -> dict:
    """Finalize plots and p-values from already decoded selected/null profiles."""
    source_dir = Path(source_dir).resolve()
    output_dir = Path(output_dir).resolve()
    design = _load_design(source_dir)
    center = design["sequence_length"] // 2
    threshold_years = design["variant_age"] * design["generation_time"]
    observed = pd.read_csv(
        output_dir / "selected_decoded_recent_probability_profile.tsv", sep="\t"
    )
    neutral_profiles_path = (
        output_dir / "neutral_decoded_recent_probability_profiles.tsv.gz"
        if neutral_profiles_path is None
        else Path(neutral_profiles_path).resolve()
    )
    neutral = pd.read_csv(neutral_profiles_path, sep="\t")
    comparison = pd.read_csv(
        output_dir / "selected_decoded_vs_truth.tsv.gz", sep="\t"
    )
    scan = calibrate_spatial_windows(observed, neutral)
    regions = significant_regions(
        scan,
        sequence_length=design["sequence_length"],
        window_size=stride,
    )
    scan.to_csv(
        output_dir / "decoded_spatial_pointwise_calibration.tsv",
        sep="\t",
        index=False,
    )
    regions.to_csv(
        output_dir / "decoded_spatial_significant_regions.tsv",
        sep="\t",
        index=False,
    )
    _plot_observed_profile(
        scan,
        center=center,
        sequence_length=design["sequence_length"],
        zoom_half_width=500_000,
        threshold_years=threshold_years,
        window_size=stride,
        output_path=output_dir / "selected_decoded_recent_probability_spatial.png",
        pair_count=design["sample_diploids"],
    )
    _plot_null_spatial_calibration(
        scan,
        regions,
        center=center,
        sequence_length=design["sequence_length"],
        zoom_half_width=500_000,
        threshold_years=threshold_years,
        window_size=stride,
        output_path=output_dir / "selected_decoded_vs_neutral_pvalues.png",
        pair_count=design["sample_diploids"],
    )
    density_calibration = scan_density_calibration(
        scan,
        neutral,
        center=center,
        local_half_width=500_000,
        output_dir=output_dir,
    )
    center_row = scan.iloc[
        np.argmin(np.abs(scan["position_0based"].to_numpy(dtype=float) - center))
    ]
    _plot_decoded_center_calibration(
        neutral,
        center_position=float(center_row["position_0based"]),
        observed_value=float(center_row["observed_fraction_recent"]),
        pvalue=float(center_row["p_upper"]),
        output_path=output_dir / "gamma_smc_center_null_and_selected_pvalue.png",
    )
    error = (
        comparison["mean_p_tmrca_lt_threshold"]
        - comparison["truth_fraction_recent"]
    )
    selected_run = {}
    selected_run_paths = [
        output_dir / "work" / "selected.tsv.run.json",
        output_dir / "selected_decoded_recent_probability_profile.tsv.run.json",
    ]
    for selected_run_path in selected_run_paths:
        if selected_run_path.exists():
            with selected_run_path.open(encoding="utf-8") as handle:
                selected_run = json.load(handle)
            break
    result = {
        "data_interpretation": "retained selected simulation treated as pseudo-empirical data",
        "recovered_from_complete_decoded_profiles": True,
        "container_image": selected_run.get("image", DEFAULT_IMAGE),
        "container_runtime": selected_run.get("runtime"),
        "stride_bp": int(stride),
        "sequence_length": design["sequence_length"],
        "sample_diploids": design["sample_diploids"],
        "within_individual_pairs": design["sample_diploids"],
        "neutral_replicates": int(neutral["replicate"].nunique()),
        "neutral_decoded_profiles_source": os.path.relpath(
            neutral_profiles_path, output_dir
        ),
        "observed_output_positions": int(len(observed)),
        "calibrated_complete_positions": int(len(scan)),
        "dropped_incomplete_null_positions": int(len(observed) - len(scan)),
        "scaled_mutation_rate": float(
            4 * design["ancestral_size"] * design["mutation_rate"]
        ),
        "recombination_to_mutation_ratio": float(
            design["recombination_rate"] / design["mutation_rate"]
        ),
        "threshold_years": float(threshold_years),
        "selected_decode_seconds": selected_run.get("decode_seconds"),
        "selected_profile_truth_bias": float(np.mean(error)),
        "selected_profile_truth_mae": float(np.mean(np.abs(error))),
        "selected_profile_truth_rmse": float(np.sqrt(np.mean(error**2))),
        "selected_profile_truth_correlation": float(
            np.corrcoef(
                comparison["mean_p_tmrca_lt_threshold"],
                comparison["truth_fraction_recent"],
            )[0, 1]
        ),
        "center_position_0based": int(center_row["position_0based"]),
        "center_observed_mean_p_recent": float(center_row["observed_fraction_recent"]),
        "center_null_mean_p_recent": float(center_row["neutral_mean_fraction_recent"]),
        "center_null_ci95_lower": float(center_row["neutral_ci95_lower"]),
        "center_null_ci95_upper": float(center_row["neutral_ci95_upper"]),
        "center_neutral_exceedances": int(center_row["neutral_exceedances"]),
        "center_monte_carlo_p_upper": float(center_row["p_upper"]),
        "n_pointwise_p_lt_0_05_windows": int(np.count_nonzero(scan["p_upper"] < 0.05)),
        "n_bh_q_lt_0_05_windows": int(np.count_nonzero(scan["q_bh"] < 0.05)),
        "n_significant_regions": int(len(regions)),
        "significant_regions": regions.to_dict("records"),
        "scan_density_calibration": density_calibration,
        "multiple_testing_note": "regions use raw pointwise p<0.05; BH q-values are in the calibration TSV",
        "workflow_elapsed_seconds_to_postprocessing_error": workflow_elapsed_seconds,
    }
    with (output_dir / "decoded_study_metrics.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(result, handle, indent=2)
    return result
