from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .calibration import calibrate_sites, calibration_metrics, monte_carlo_pvalue, randomized_rank_pvalue
from .decoder import run_within_decoder
from .evaluation import evaluate_pairs
from .plotting import plot_scan, plot_truth_tmrca_relationship
from .simulation import SimulationConfig, simulate_replicates
from .selection import validate_slim_hard_sweep
from .tree_sequence import tree_sequence_to_vcf


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
    decode.add_argument("--executable", default="gamma_smc")
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
    decode.set_defaults(func=command_decode)

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
    return root


def main(argv=None):
    args = parser().parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
