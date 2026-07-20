from __future__ import annotations

import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from time import perf_counter

import numpy as np
import pandas as pd

from .container_decoder import DEFAULT_IMAGE, run_container_decoder
from .selection import run_slim_recent_sweep, within_individual_tmrca_grid
from .spatial_scan import (
    _plot_null_spatial_calibration,
    _plot_observed_profile,
    calibrate_spatial_windows,
    significant_regions,
)
from .tree_sequence import tree_sequence_to_vcf


def _load_design(source_dir: Path) -> dict:
    with (source_dir / "metrics.json").open(encoding="utf-8") as handle:
        metrics = json.load(handle)
    selected = pd.read_csv(source_dir / "selected_observation.tsv", sep="\t").iloc[0]
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


def run_container_stride_study(
    source_dir: str | Path,
    output_dir: str | Path,
    *,
    neutral_replicates: int = 100,
    stride: int = 1_000,
    runtime: str = "auto",
    image: str = DEFAULT_IMAGE,
    keep_vcfs: bool = False,
    workers: int = 1,
) -> dict:
    """Decode a retained sweep and matched nulls with official Gamma-SMC v0.2."""
    if neutral_replicates < 1 or stride < 1 or workers < 1:
        raise ValueError("replicate count, stride, and workers must be positive")
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
    started = perf_counter()

    selected_tree_path = work_dir / "selected.trees"
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
    selected_vcf = work_dir / "selected.vcf.gz"
    tree_sequence_to_vcf(selected_tree_path, selected_vcf)
    selected_summary = work_dir / "selected.tsv"
    selected_decode = run_container_decoder(
        selected_vcf,
        selected_summary,
        scaled_mutation_rate=theta,
        recombination_to_mutation_ratio=rho_over_theta,
        mutation_rate=design["mutation_rate"],
        threshold_years=threshold_years,
        generation_time=design["generation_time"],
        stride=stride,
        runtime=runtime,
        image=image,
    )
    observed = pd.read_csv(selected_summary, sep="\t")
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

    def neutral_task(replicate: int) -> tuple[pd.DataFrame, dict]:
        tree_path = work_dir / f"neutral_{replicate:04d}.trees"
        neutral_ts, _ = _simulate(
            tree_path,
            design,
            selection_coefficient=0.0,
            seed=design["seed"] + 1_000_000 + replicate * 10,
        )
        del neutral_ts
        vcf_path = work_dir / f"neutral_{replicate:04d}.vcf.gz"
        tree_sequence_to_vcf(tree_path, vcf_path)
        summary_path = work_dir / f"neutral_{replicate:04d}.tsv"
        run = run_container_decoder(
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
        profile = pd.read_csv(summary_path, sep="\t")
        profile["replicate"] = replicate
        if not keep_vcfs:
            vcf_path.unlink(missing_ok=True)
            tree_path.unlink(missing_ok=True)
            summary_path.unlink(missing_ok=True)
            summary_path.with_name(summary_path.name + ".run.json").unlink(
                missing_ok=True
            )
        return profile, {"replicate": replicate, **run}

    if workers == 1:
        neutral_results = [neutral_task(value) for value in range(neutral_replicates)]
    else:
        with ThreadPoolExecutor(max_workers=min(workers, neutral_replicates)) as executor:
            neutral_results = list(executor.map(neutral_task, range(neutral_replicates)))
    neutral_profiles = [result[0] for result in neutral_results]
    neutral_runs = [result[1] for result in neutral_results]

    neutral = pd.concat(neutral_profiles, ignore_index=True)
    neutral.to_csv(
        output_dir / "neutral_decoded_recent_probability_profiles.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
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
    center_row = scan.iloc[
        np.argmin(np.abs(scan["position_0based"].to_numpy(dtype=float) - center))
    ]
    decode_times = np.asarray(
        [selected_decode["decode_seconds"]]
        + [run["decode_seconds"] for run in neutral_runs]
    )
    result = {
        "data_interpretation": "retained selected simulation treated as pseudo-empirical data",
        "container_image": image,
        "container_runtime": selected_decode["runtime"],
        "stride_bp": int(stride),
        "sequence_length": design["sequence_length"],
        "sample_diploids": design["sample_diploids"],
        "within_individual_pairs": design["sample_diploids"],
        "neutral_replicates": int(neutral_replicates),
        "neutral_workers": int(min(workers, neutral_replicates)),
        "scaled_mutation_rate": float(theta),
        "recombination_to_mutation_ratio": float(rho_over_theta),
        "threshold_years": float(threshold_years),
        "selected_decode_seconds": float(selected_decode["decode_seconds"]),
        "neutral_decode_seconds_mean": float(
            np.mean([run["decode_seconds"] for run in neutral_runs])
        ),
        "all_decode_seconds_sum": float(decode_times.sum()),
        "selected_profile_truth_mae": float(np.mean(np.abs(error))),
        "selected_profile_truth_rmse": float(np.sqrt(np.mean(error**2))),
        "center_position_0based": int(center_row["position_0based"]),
        "center_observed_mean_p_recent": float(center_row["observed_fraction_recent"]),
        "center_null_mean_p_recent": float(center_row["neutral_mean_fraction_recent"]),
        "center_neutral_exceedances": int(center_row["neutral_exceedances"]),
        "center_monte_carlo_p_upper": float(center_row["p_upper"]),
        "n_pointwise_p_lt_0_05_windows": int(np.count_nonzero(scan["p_upper"] < 0.05)),
        "n_significant_regions": int(len(regions)),
        "significant_regions": regions.to_dict("records"),
        "multiple_testing_note": "regions use raw pointwise p<0.05; BH q-values are in the calibration TSV",
        "elapsed_seconds": float(perf_counter() - started),
    }
    with (output_dir / "decoded_study_metrics.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(result, handle, indent=2)
    if not keep_vcfs:
        selected_vcf.unlink(missing_ok=True)
        selected_tree_path.unlink(missing_ok=True)
    return result
