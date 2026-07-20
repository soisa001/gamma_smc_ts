from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .calibration import calibrate_sites, calibration_metrics, monte_carlo_pvalue, randomized_rank_pvalue
from .carrier_profiles import plot_retained_carrier_tmrca_profiles
from .container_decoder import DEFAULT_IMAGE, run_container_decoder
from .container_study import run_container_stride_study
from .decoder import run_within_decoder
from .evaluation import evaluate_pairs
from .plotting import plot_scan, plot_truth_tmrca_relationship
from .simulation import SimulationConfig, simulate_replicates
from .selection import (
    analyze_retained_recent_sweeps,
    validate_recent_sweep_grid,
    validate_slim_hard_sweep,
)
from .spatial_scan import analyze_two_epoch_spatial_truth
from .tree_sequence import tree_sequence_to_vcf
from .two_epoch import validate_two_epoch_growth


def _histories(path: str | None):
    if path is None:
        return None
    table = pd.read_csv(path, sep="\t")
    required = {"draw", "time_generations", "ne"}
    if not required.issubset(table.columns):
        raise ValueError(f"history TSV needs columns: {sorted(required)}")
    return [
        (group["time_generations"].tolist(), group["ne"].tolist())
        for _, group in table.sort_values(["draw", "time_generations"]).groupby("draw", sort=True)
    ]


def command_simulate(args):
    config = SimulationConfig(
        n_replicates=args.replicates,
        n_diploids=args.diploids,
        sequence_length=args.length,
        effective_size=args.ne,
        mutation_rate=args.mutation_rate,
        recombination_rate=args.recombination_rate,
        threshold_years=args.threshold_years,
        generation_time=args.generation_time,
        seed=args.seed,
        save_trees=args.save_trees,
    )
    simulate_replicates(
        config, args.output_dir, histories=_histories(args.histories),
        mutation_map=args.mutation_map, recombination_map=args.recombination_map,
        workers=args.workers,
    )


def _read_summaries(pattern: str) -> list[pd.DataFrame]:
    files = sorted(Path().glob(pattern)) if not Path(pattern).is_absolute() else sorted(Path(pattern).parent.glob(Path(pattern).name))
    if not files:
        raise FileNotFoundError(f"no summary TSVs matched {pattern}")
    return [pd.read_csv(path, sep="\t") for path in files]


def command_calibrate(args):
    observed = pd.read_csv(args.observed, sep="\t")
    simulations = _read_summaries(args.sim_glob)
    result, null = calibrate_sites(
        observed, simulations, observed_length=args.observed_length,
        simulation_lengths=[args.simulation_length] * len(simulations), match=args.match,
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output, sep="\t", index=False)
    np.save(output.with_suffix(".null.npy"), null)
    plot_scan(result, null, output.parent / "plots")


def command_validate_null(args):
    simulations = _read_summaries(args.sim_glob)
    values = []
    for frame in simulations:
        if frame.empty:
            values.append(np.nan)
            continue
        positions = frame["position_0based"].to_numpy(dtype=float) / args.sequence_length
        values.append(float(frame.iloc[np.argmin(np.abs(positions - args.relative_position))]["mean_p_tmrca_lt_threshold"]))
    values = np.asarray(values)
    pvalues = np.asarray([
        monte_carlo_pvalue(values[i], np.delete(values, i)) for i in range(len(values))
    ])
    rng = np.random.default_rng(args.tie_seed)
    randomized = np.asarray([
        randomized_rank_pvalue(values[i], np.delete(values, i), rng.random()) for i in range(len(values))
    ])
    metrics = {
        "conservative_upper_tail": calibration_metrics(pvalues),
        "randomized_tie_diagnostic": calibration_metrics(randomized),
        "n_unique_statistics": int(np.unique(values[np.isfinite(values)]).size),
    }
    destination = Path(args.output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({
        "replicate": np.arange(len(values)), "statistic": values,
        "leave_one_out_p": pvalues, "randomized_tie_p": randomized,
    }).to_csv(
        destination / "null_rank_calibration.tsv", sep="\t", index=False
    )
    with (destination / "null_rank_metrics.json").open("w", encoding="utf-8") as handle:
        json.dump(metrics, handle, indent=2)
    fig, ax = plt.subplots(figsize=(5, 5), constrained_layout=True)
    ordered = np.sort(pvalues[np.isfinite(pvalues)])
    ordered_random = np.sort(randomized[np.isfinite(randomized)])
    expected = (np.arange(len(ordered)) + 0.5) / len(ordered)
    ax.scatter(expected, ordered, s=12)
    ax.scatter(expected, ordered_random, s=12, alpha=0.7, label="randomized tie diagnostic")
    ax.plot([0, 1], [0, 1], color="grey")
    ax.set(xlabel="Expected uniform quantile", ylabel="Leave-one-out p", title="Neutral rank calibration")
    ax.legend()
    fig.savefig(destination / "null_rank_qq.png", dpi=160)
    plt.close(fig)


def command_convert(args):
    tree_sequence_to_vcf(args.input, args.output)


def command_decode(args):
    run_within_decoder(
        args.executable, args.input, args.output,
        scaled_mutation_rate=args.theta,
        recombination_to_mutation_ratio=args.rho_over_theta,
        mutation_rate=args.mutation_rate,
        threshold_years=args.threshold_years,
        generation_time=args.generation_time,
        input_format=args.input_format,
        raw_output=args.raw_output,
        mask=args.mask,
        masks_per_sample=args.masks_per_sample,
        output_at_stride=args.output_at_stride,
        output_at_hets=not args.no_output_at_hets,
    )


def command_decode_container(args):
    run_container_decoder(
        args.input,
        args.output,
        scaled_mutation_rate=args.theta,
        recombination_to_mutation_ratio=args.rho_over_theta,
        mutation_rate=args.mutation_rate,
        threshold_years=args.threshold_years,
        generation_time=args.generation_time,
        stride=args.output_at_stride,
        runtime=args.runtime,
        image=args.image,
        keep_raw=args.keep_raw,
    )


def command_container_study(args):
    run_container_stride_study(
        args.source_dir,
        args.output_dir,
        neutral_replicates=args.neutral_replicates,
        stride=args.output_at_stride,
        runtime=args.runtime,
        image=args.image,
        keep_vcfs=args.keep_vcfs,
        workers=args.workers,
    )


def command_evaluate(args):
    truth = sorted(Path(args.truth_dir).glob("*.tsv"))
    decoded = sorted(Path(args.decoded_dir).glob("*.tsv"))
    truth_by_name = {path.name: path for path in truth}
    decoded_by_name = {path.name: path for path in decoded}
    names = sorted(truth_by_name.keys() & decoded_by_name.keys())
    if not names:
        raise FileNotFoundError("no identically named truth and decoded TSVs")
    evaluate_pairs([(truth_by_name[name], decoded_by_name[name]) for name in names], args.output_dir)


def command_plot_truth(args):
    summaries = sorted(Path(args.summary_dir).glob("*.tsv"))
    if not summaries:
        raise FileNotFoundError("no truth summary TSVs")
    plot_truth_tmrca_relationship(
        summaries,
        sequence_length=args.sequence_length,
        effective_size=args.ne,
        threshold_generations=args.threshold_years / args.generation_time,
        relative_position=args.relative_position,
        output_dir=args.output_dir,
    )


def command_validate_sweep(args):
    validate_slim_hard_sweep(
        args.output_dir,
        executable=args.slim,
        population_size=args.population_size,
        sequence_length=args.sequence_length,
        selection_coefficient=args.selection_coefficient,
        recombination_rate=args.recombination_rate,
        threshold_generations=args.threshold_years / args.generation_time,
        neutral_replicates=args.neutral_replicates,
        seed=args.seed,
    )


def command_validate_recent_sweep(args):
    validate_recent_sweep_grid(
        args.output_dir,
        executable=args.slim,
        population_size=args.population_size,
        sample_diploids=args.sample_diploids,
        sequence_length=args.sequence_length,
        selection_coefficients=tuple(args.selection_coefficients),
        age_generations=args.age_generations,
        mutation_rate=args.mutation_rate,
        recombination_rate=args.recombination_rate,
        neutral_replicates=args.neutral_replicates,
        selected_replicates=args.selected_replicates,
        workers=args.workers,
        reuse_null_from=args.reuse_null_from,
        save_trees=not args.no_save_trees,
        seed=args.seed,
    )


def command_analyze_retained_sweeps(args):
    analyze_retained_recent_sweeps(
        args.source_dir,
        args.output_dir,
        executable=args.slim,
        selection_coefficient=args.selection_coefficient,
        retained_replicates=args.retained_replicates,
        workers=args.workers,
        seed=args.seed,
    )


def command_plot_retained_carrier_profiles(args):
    plot_retained_carrier_tmrca_profiles(
        args.source_dir,
        args.output_dir,
        executable=args.slim,
        selection_coefficient=args.selection_coefficient,
        retained_replicates=args.retained_replicates,
        workers=args.workers,
        seed=args.seed,
        full_step=args.full_step,
        zoom_half_width=args.zoom_half_width,
        zoom_step=args.zoom_step,
    )


def command_validate_two_epoch_growth(args):
    validate_two_epoch_growth(
        args.output_dir,
        executable=args.slim,
        ancestral_population_size=args.ancestral_population_size,
        present_population_size=args.present_population_size,
        size_change_generations_ago=args.size_change_generations_ago,
        sample_diploids=args.sample_diploids,
        sequence_length=args.sequence_length,
        variant_age_generations=args.variant_age_generations,
        generation_time_years=args.generation_time_years,
        selection_coefficient=args.selection_coefficient,
        mutation_rate=args.mutation_rate,
        recombination_rate=args.recombination_rate,
        neutral_replicates=args.neutral_replicates,
        workers=args.workers,
        max_selected_attempts=args.max_selected_attempts,
        minimum_hom_alt_pairs=args.minimum_hom_alt_pairs,
        full_step=args.full_step,
        zoom_half_width=args.zoom_half_width,
        zoom_step=args.zoom_step,
        seed=args.seed,
    )


def command_two_epoch_spatial_truth(args):
    analyze_two_epoch_spatial_truth(
        args.source_dir,
        args.output_dir,
        executable=args.slim,
        workers=args.workers,
        window_size=args.window_size,
        zoom_half_width=args.zoom_half_width,
    )


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(prog="gamma-smc-aou")
    commands = root.add_subparsers(required=True)
    sim = commands.add_parser("simulate", help="simulate neutral regions with StandardCoalescent")
    sim.add_argument("--output-dir", required=True)
    sim.add_argument("--replicates", type=int, default=1000)
    sim.add_argument("--diploids", type=int, default=2000)
    sim.add_argument("--length", type=int, default=1_000_000)
    sim.add_argument("--ne", type=float, default=10_000)
    sim.add_argument("--mutation-rate", type=float, default=1.25e-8)
    sim.add_argument("--recombination-rate", type=float, default=1e-8)
    sim.add_argument("--mutation-map")
    sim.add_argument("--recombination-map")
    sim.add_argument("--histories", help="TSV with draw,time_generations,ne (e.g. PHLASH MVN draws)")
    sim.add_argument("--threshold-years", type=float, default=4500)
    sim.add_argument("--generation-time", type=float, default=30)
    sim.add_argument("--seed", type=int, default=1729)
    sim.add_argument("--save-trees", action="store_true")
    sim.add_argument("--workers", type=int, default=1, help="independent simulation processes")
    sim.set_defaults(func=command_simulate)

    cal = commands.add_parser("calibrate", help="calculate site-level simulation p-values")
    cal.add_argument("--observed", required=True)
    cal.add_argument("--sim-glob", required=True)
    cal.add_argument("--observed-length", type=float, required=True)
    cal.add_argument("--simulation-length", type=float, required=True)
    cal.add_argument("--match", choices=["nearest", "region_max"], default="nearest")
    cal.add_argument("--output", required=True)
    cal.set_defaults(func=command_calibrate)

    validate = commands.add_parser("validate-null", help="leave-one-replicate-out p-value calibration")
    validate.add_argument("--sim-glob", required=True)
    validate.add_argument("--sequence-length", type=float, required=True)
    validate.add_argument("--relative-position", type=float, default=0.5)
    validate.add_argument("--tie-seed", type=int, default=8675309)
    validate.add_argument("--output-dir", required=True)
    validate.set_defaults(func=command_validate_null)

    convert = commands.add_parser("convert-ts", help="convert .trees/.tsz to diploid VCF")
    convert.add_argument("--input", required=True)
    convert.add_argument("--output", required=True)
    convert.set_defaults(func=command_convert)

    decode = commands.add_parser("decode", help="run streaming within-individual Gamma-SMC summary")
    decode.add_argument(
        "--executable",
        default=os.environ.get("GAMMA_SMC_BIN", "gamma_smc"),
        help="Gamma-SMC binary (default: GAMMA_SMC_BIN or gamma_smc on PATH)",
    )
    decode.add_argument("--input", required=True)
    decode.add_argument("--input-format", choices=["auto", "vcf", "trees", "tsz"], default="auto")
    decode.add_argument("--output", required=True)
    decode.add_argument("--raw-output")
    decode.add_argument("--mask", help="global BED mask; use the same callable-region policy in data and simulations")
    decode.add_argument("--masks-per-sample", help="sample-to-BED TSV")
    decode.add_argument("--theta", type=float, required=True)
    decode.add_argument("--rho-over-theta", type=float, required=True)
    decode.add_argument("--mutation-rate", type=float, required=True)
    decode.add_argument("--threshold-years", type=float, default=4500)
    decode.add_argument("--generation-time", type=float, default=30)
    decode.add_argument("--output-at-stride", type=int, default=-1)
    decode.add_argument("--no-output-at-hets", action="store_true")
    decode.set_defaults(func=command_decode)

    container = commands.add_parser(
        "decode-container",
        help="decode a VCF with the official Gamma-SMC v0.2 container",
    )
    container.add_argument("--input", required=True)
    container.add_argument("--output", required=True)
    container.add_argument("--theta", type=float, required=True)
    container.add_argument("--rho-over-theta", type=float, required=True)
    container.add_argument("--mutation-rate", type=float, required=True)
    container.add_argument("--threshold-years", type=float, default=4500)
    container.add_argument("--generation-time", type=float, default=25)
    container.add_argument("--output-at-stride", type=int, default=1000)
    container.add_argument(
        "--runtime", choices=["auto", "apptainer", "singularity", "docker"], default="auto"
    )
    container.add_argument("--image", default=DEFAULT_IMAGE)
    container.add_argument("--keep-raw", action="store_true")
    container.set_defaults(func=command_decode_container)

    study = commands.add_parser(
        "run-container-study",
        help="decode the retained selected simulation and matched nulls at fixed stride",
    )
    study.add_argument("--source-dir", required=True)
    study.add_argument("--output-dir", required=True)
    study.add_argument("--neutral-replicates", type=int, default=100)
    study.add_argument("--output-at-stride", type=int, default=1000)
    study.add_argument(
        "--runtime", choices=["auto", "apptainer", "singularity", "docker"], default="auto"
    )
    study.add_argument("--image", default=DEFAULT_IMAGE)
    study.add_argument("--keep-vcfs", action="store_true")
    study.add_argument("--workers", type=int, default=1)
    study.set_defaults(func=command_container_study)

    evaluate = commands.add_parser("evaluate-decoder", help="compare decoded simulations with tree-sequence truth")
    evaluate.add_argument("--truth-dir", required=True)
    evaluate.add_argument("--decoded-dir", required=True)
    evaluate.add_argument("--output-dir", required=True)
    evaluate.set_defaults(func=command_evaluate)

    truth_plot = commands.add_parser("plot-truth", help="plot true recent fraction against simulated mean TMRCA")
    truth_plot.add_argument("--summary-dir", required=True)
    truth_plot.add_argument("--sequence-length", type=float, required=True)
    truth_plot.add_argument("--ne", type=float, required=True)
    truth_plot.add_argument("--threshold-years", type=float, default=4500)
    truth_plot.add_argument("--generation-time", type=float, default=30)
    truth_plot.add_argument("--relative-position", type=float, default=0.5)
    truth_plot.add_argument("--output-dir", required=True)
    truth_plot.set_defaults(func=command_plot_truth)

    sweep = commands.add_parser("validate-sweep", help="validate recent-coalescence power with a SLiM hard sweep")
    sweep.add_argument("--output-dir", required=True)
    sweep.add_argument("--slim", help="SLiM executable; otherwise use SLIM_BIN/PATH")
    sweep.add_argument("--population-size", type=int, default=200)
    sweep.add_argument("--sequence-length", type=int, default=100_000)
    sweep.add_argument("--selection-coefficient", type=float, default=0.5)
    sweep.add_argument("--recombination-rate", type=float, default=1e-7)
    sweep.add_argument("--threshold-years", type=float, default=4500)
    sweep.add_argument("--generation-time", type=float, default=30)
    sweep.add_argument("--neutral-replicates", type=int, default=39)
    sweep.add_argument("--seed", type=int, default=24681357)
    sweep.set_defaults(func=command_validate_sweep)

    recent = commands.add_parser(
        "validate-recent-sweep",
        help="validate a 10-Mb, age-controlled recent SLiM sweep against a neutral null",
    )
    recent.add_argument("--output-dir", required=True)
    recent.add_argument("--slim", help="SLiM executable; otherwise use SLIM_BIN/PATH")
    recent.add_argument("--population-size", type=int, default=10_000)
    recent.add_argument("--sample-diploids", type=int, default=2_000)
    recent.add_argument("--sequence-length", type=int, default=10_000_000)
    recent.add_argument(
        "--selection-coefficients", type=float, nargs="+", default=[0.0, 0.001, 0.01]
    )
    recent.add_argument("--age-generations", type=int, default=180)
    recent.add_argument("--mutation-rate", type=float, default=1.25e-8)
    recent.add_argument("--recombination-rate", type=float, default=1e-8)
    recent.add_argument("--neutral-replicates", type=int, default=100)
    recent.add_argument("--selected-replicates", type=int, default=1)
    recent.add_argument("--workers", type=int, default=1, help="parallel selected SLiM trajectories")
    recent.add_argument(
        "--reuse-null-from",
        help="reuse compatible s=0 replicate statistics/profiles from this result directory",
    )
    recent.add_argument("--no-save-trees", action="store_true")
    recent.add_argument("--seed", type=int, default=271828)
    recent.set_defaults(func=command_validate_recent_sweep)

    retained = commands.add_parser(
        "analyze-retained-sweeps",
        help="compare retained sweeps to a saved null and reconstruct carrier TMRCAs",
    )
    retained.add_argument("--source-dir", required=True)
    retained.add_argument("--output-dir", required=True)
    retained.add_argument("--slim", help="SLiM executable; otherwise use SLIM_BIN/PATH")
    retained.add_argument("--selection-coefficient", type=float, default=0.1)
    retained.add_argument("--retained-replicates", type=int, default=10)
    retained.add_argument("--workers", type=int, default=20)
    retained.add_argument(
        "--seed",
        type=int,
        default=424242,
        help="base seed used for the source validate-recent-sweep run",
    )
    retained.set_defaults(func=command_analyze_retained_sweeps)

    carrier_profiles = commands.add_parser(
        "plot-retained-carrier-profiles",
        help="plot hom-alt and hom-ref TMRCA profiles for retained sweeps",
    )
    carrier_profiles.add_argument("--source-dir", required=True)
    carrier_profiles.add_argument("--output-dir", required=True)
    carrier_profiles.add_argument(
        "--slim", help="SLiM executable; otherwise use SLIM_BIN/PATH"
    )
    carrier_profiles.add_argument("--selection-coefficient", type=float, default=0.1)
    carrier_profiles.add_argument("--retained-replicates", type=int, default=10)
    carrier_profiles.add_argument("--workers", type=int, default=20)
    carrier_profiles.add_argument("--seed", type=int, default=424242)
    carrier_profiles.add_argument("--full-step", type=int, default=50_000)
    carrier_profiles.add_argument("--zoom-half-width", type=int, default=500_000)
    carrier_profiles.add_argument("--zoom-step", type=int, default=5_000)
    carrier_profiles.set_defaults(func=command_plot_retained_carrier_profiles)

    two_epoch = commands.add_parser(
        "validate-two-epoch-growth",
        help="run a two-epoch neutral null and rejection-sampled selected validation",
    )
    two_epoch.add_argument("--output-dir", required=True)
    two_epoch.add_argument(
        "--slim", help="SLiM executable; otherwise use SLIM_BIN/PATH"
    )
    two_epoch.add_argument("--ancestral-population-size", type=int, default=10_000)
    two_epoch.add_argument("--present-population-size", type=int, default=20_000)
    two_epoch.add_argument("--size-change-generations-ago", type=int, default=100)
    two_epoch.add_argument("--sample-diploids", type=int, default=2_000)
    two_epoch.add_argument("--sequence-length", type=int, default=10_000_000)
    two_epoch.add_argument("--variant-age-generations", type=int, default=180)
    two_epoch.add_argument("--generation-time-years", type=float, default=25)
    two_epoch.add_argument("--selection-coefficient", type=float, default=0.05)
    two_epoch.add_argument("--mutation-rate", type=float, default=1.25e-8)
    two_epoch.add_argument("--recombination-rate", type=float, default=1e-8)
    two_epoch.add_argument("--neutral-replicates", type=int, default=100)
    two_epoch.add_argument("--workers", type=int, default=20)
    two_epoch.add_argument("--max-selected-attempts", type=int, default=1_000)
    two_epoch.add_argument("--minimum-hom-alt-pairs", type=int, default=2)
    two_epoch.add_argument("--full-step", type=int, default=50_000)
    two_epoch.add_argument("--zoom-half-width", type=int, default=500_000)
    two_epoch.add_argument("--zoom-step", type=int, default=5_000)
    two_epoch.add_argument("--seed", type=int, default=515151)
    two_epoch.set_defaults(func=command_validate_two_epoch_growth)

    spatial = commands.add_parser(
        "analyze-two-epoch-spatial-truth",
        help="reconstruct two-epoch selected/null truth profiles and spatial p-values",
    )
    spatial.add_argument("--source-dir", required=True)
    spatial.add_argument("--output-dir")
    spatial.add_argument(
        "--slim", help="SLiM executable; otherwise use SLIM_BIN/PATH"
    )
    spatial.add_argument("--workers", type=int, default=20)
    spatial.add_argument("--window-size", type=int, default=20_000)
    spatial.add_argument("--zoom-half-width", type=int, default=500_000)
    spatial.set_defaults(func=command_two_epoch_spatial_truth)
    return root


def main(argv=None):
    args = parser().parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
