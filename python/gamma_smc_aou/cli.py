from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .bitmatrix import to_frame as bitmatrix_to_frame
from .calibration import (
    calibrate_sites,
    calibration_metrics,
    monte_carlo_pvalue,
    randomized_rank_pvalue,
)
from .carrier_profiles import plot_retained_carrier_tmrca_profiles
from .container_decoder import DEFAULT_IMAGE, run_container_decoder
from .container_study import (
    finalize_container_stride_study,
    run_container_stride_study,
    run_native_stride_study,
)
from .defaults import (
    DEFAULT_CACHE_SIZE,
    DEFAULT_GENERATION_TIME,
    DEFAULT_MUTATION_RATE,
    DEFAULT_OUTPUT_STRIDE,
    DEFAULT_RECOMBINATION_TO_MUTATION_RATIO,
    DEFAULT_SCALED_MUTATION_RATE,
)
from .high_af_study import finalize_high_af_selected, prepare_high_af_selected
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
from .workbench import (
    AUTOSOMES,
    POPULATIONS,
    build_workbench_callable_mask,
    build_workbench_sample_list,
    build_workbench_contract,
    plot_workbench_population,
    summarize_workbench_run,
    validate_workbench_completion,
    write_workbench_completion,
)
from .workbench_candidates import (
    analyze_candidate_regions,
    build_candidate_regions,
    write_empty_candidate_analysis,
)


def _histories(path: str | None):
    if path is None:
        return None
    table = pd.read_csv(path, sep="\t")
    required = {"draw", "time_generations", "ne"}
    if not required.issubset(table.columns):
        raise ValueError(f"history TSV needs columns: {sorted(required)}")
    return [
        (group["time_generations"].tolist(), group["ne"].tolist())
        for _, group in table.sort_values(["draw", "time_generations"]).groupby(
            "draw", sort=True
        )
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
        config,
        args.output_dir,
        histories=_histories(args.histories),
        mutation_map=args.mutation_map,
        recombination_map=args.recombination_map,
        workers=args.workers,
    )


def _read_summaries(pattern: str) -> list[pd.DataFrame]:
    files = (
        sorted(Path().glob(pattern))
        if not Path(pattern).is_absolute()
        else sorted(Path(pattern).parent.glob(Path(pattern).name))
    )
    if not files:
        raise FileNotFoundError(f"no summary TSVs matched {pattern}")
    return [pd.read_csv(path, sep="\t") for path in files]


def command_calibrate(args):
    observed = pd.read_csv(args.observed, sep="\t")
    simulations = _read_summaries(args.sim_glob)
    result, null = calibrate_sites(
        observed,
        simulations,
        observed_length=args.observed_length,
        simulation_lengths=[args.simulation_length] * len(simulations),
        match=args.match,
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
        positions = (
            frame["position_0based"].to_numpy(dtype=float) / args.sequence_length
        )
        values.append(
            float(
                frame.iloc[np.argmin(np.abs(positions - args.relative_position))][
                    "mean_p_tmrca_lt_threshold"
                ]
            )
        )
    values = np.asarray(values)
    pvalues = np.asarray(
        [
            monte_carlo_pvalue(values[i], np.delete(values, i))
            for i in range(len(values))
        ]
    )
    rng = np.random.default_rng(args.tie_seed)
    randomized = np.asarray(
        [
            randomized_rank_pvalue(values[i], np.delete(values, i), rng.random())
            for i in range(len(values))
        ]
    )
    metrics = {
        "conservative_upper_tail": calibration_metrics(pvalues),
        "randomized_tie_diagnostic": calibration_metrics(randomized),
        "n_unique_statistics": int(np.unique(values[np.isfinite(values)]).size),
    }
    destination = Path(args.output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        {
            "replicate": np.arange(len(values)),
            "statistic": values,
            "leave_one_out_p": pvalues,
            "randomized_tie_p": randomized,
        }
    ).to_csv(destination / "null_rank_calibration.tsv", sep="\t", index=False)
    with (destination / "null_rank_metrics.json").open("w", encoding="utf-8") as handle:
        json.dump(metrics, handle, indent=2)
    fig, ax = plt.subplots(figsize=(5, 5), constrained_layout=True)
    ordered = np.sort(pvalues[np.isfinite(pvalues)])
    ordered_random = np.sort(randomized[np.isfinite(randomized)])
    expected = (np.arange(len(ordered)) + 0.5) / len(ordered)
    ax.scatter(expected, ordered, s=12)
    ax.scatter(
        expected, ordered_random, s=12, alpha=0.7, label="randomized tie diagnostic"
    )
    ax.plot([0, 1], [0, 1], color="grey")
    ax.set(
        xlabel="Expected uniform quantile",
        ylabel="Leave-one-out p",
        title="Neutral rank calibration",
    )
    ax.legend()
    fig.savefig(destination / "null_rank_qq.png", dpi=160)
    plt.close(fig)


def command_convert(args):
    tree_sequence_to_vcf(args.input, args.output)


def command_decode(args):
    run_within_decoder(
        args.executable,
        args.input,
        args.output,
        scaled_mutation_rate=args.theta,
        recombination_to_mutation_ratio=args.rho_over_theta,
        mutation_rate=args.mutation_rate,
        threshold_years=args.threshold_years,
        generation_time=args.generation_time,
        input_format=args.input_format,
        raw_output=args.raw_output,
        bitmatrix_output=args.bitmatrix,
        mask=args.mask,
        masks_per_sample=args.masks_per_sample,
        samples=args.samples,
        output_positions_file=args.output_positions_file,
        output_at_stride=args.output_at_stride,
        output_at_hets=not args.no_output_at_hets,
        only_within=(args.n_random_pairs <= 0 and args.pairs_file is None),
        n_random_pairs=args.n_random_pairs,
        pairs_seed=args.pairs_seed,
        pairs_file=args.pairs_file,
        exclude_within=args.exclude_within,
        recent_call=args.recent_call,
        recent_call_probability=args.recent_call_probability,
        threads=args.threads,
        cache_size=args.cache_size,
        pair_block=args.pair_block,
        pairs_manifest=args.pairs_manifest,
        exp10=args.exp10,
        backward_alignment=args.backward_alignment,
    )


def command_bitmatrix_summary(args):
    frame = bitmatrix_to_frame(args.input, pair_indices=_pair_indices(args.pairs))
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.output, sep="\t", index=False)
    print(f"wrote {len(frame)} positions to {args.output}")


def _pair_indices(path: str | None):
    if path is None:
        return None
    return np.loadtxt(path, dtype=np.int64, ndmin=1)


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


def _workbench_contract_from_args(args) -> dict:
    return build_workbench_contract(
        population=args.population,
        chromosome=args.chromosome,
        input_uri=args.input_uri,
        input_fingerprint=args.input_fingerprint,
        local_input=args.local_input,
        index_uri=args.index_uri,
        index_fingerprint=args.index_fingerprint,
        local_index=args.local_index,
        output_summary=args.summary,
        pairs_manifest=args.pairs_manifest,
        sample_list=args.sample_list,
        sample_audit=args.sample_audit,
        ancestry_uri=args.ancestry_uri,
        ancestry_fingerprint=args.ancestry_fingerprint,
        local_ancestry=args.local_ancestry,
        qc_exclusions_uri=args.qc_exclusions_uri,
        qc_exclusions_fingerprint=args.qc_exclusions_fingerprint,
        local_qc_exclusions=args.local_qc_exclusions,
        relatedness_exclusions_uri=args.relatedness_exclusions_uri,
        relatedness_exclusions_fingerprint=args.relatedness_exclusions_fingerprint,
        local_relatedness_exclusions=args.local_relatedness_exclusions,
        mask_uri=args.mask_uri,
        mask_fingerprint=args.mask_fingerprint,
        local_mask_source=args.local_mask_source,
        local_mask=args.local_mask,
        mask_audit=args.mask_audit,
        theta=args.theta,
        rho_over_theta=args.rho_over_theta,
        mutation_rate=args.mutation_rate,
        generation_time=args.generation_time,
        threshold_years=args.threshold_years,
        recent_call=args.recent_call,
        stride=args.output_at_stride,
        cache_size=args.cache_size,
        threads=args.threads,
        pair_block=args.pair_block,
        exp10=args.exp10,
        backward_alignment=args.backward_alignment,
        code_commit=args.code_commit,
        bitmatrix=args.bitmatrix,
        candidate_regions=args.candidate_regions,
        candidate_manifest=args.candidate_manifest,
        n_random_pairs=args.n_random_pairs,
        pairs_seed=args.pairs_seed,
        exclude_within=args.exclude_within,
        signal_fraction=args.signal_fraction,
        merge_gap=args.merge_gap,
        profile_half_width=args.profile_half_width,
        variant_half_width=args.variant_half_width,
        minimum_genotype_pairs=args.minimum_genotype_pairs,
    )


def command_workbench_samples(args):
    audit = build_workbench_sample_list(
        ancestry_path=args.ancestry,
        qc_exclusions_path=args.qc_exclusions,
        relatedness_exclusions_path=args.relatedness_exclusions,
        bcf_samples_path=args.bcf_samples,
        population=args.population,
        output_path=args.output,
        audit_path=args.audit_output,
    )
    print(
        f"selected {audit['counts']['selected_samples']} {audit['population']} samples "
        f"from {audit['counts']['bcf_samples']} BCF samples"
    )


def command_workbench_mask(args):
    audit = build_workbench_callable_mask(
        hardmask_path=args.hardmask,
        contig=args.contig,
        sequence_length=args.sequence_length,
        output_path=args.output,
        audit_path=args.audit_output,
    )
    print(
        f"prepared {audit['contig']} callable mask: "
        f"{audit['counts']['callable_bases']} / {audit['sequence_length']} bp callable"
    )


def command_workbench_validate(args):
    contract = _workbench_contract_from_args(args)
    try:
        if args.check_only:
            completion = validate_workbench_completion(
                summary_path=args.summary,
                run_json_path=args.run_json,
                completion_path=args.completion,
                contract=contract,
            )
        else:
            completion = write_workbench_completion(
                summary_path=args.summary,
                run_json_path=args.run_json,
                completion_path=args.completion,
                contract=contract,
            )
    except (KeyError, OSError, TypeError, ValueError, json.JSONDecodeError) as error:
        if args.check_only:
            print(f"incomplete: {error}", file=sys.stderr)
            raise SystemExit(1) from error
        raise
    print(
        f"validated {contract['population']} chr{contract['chromosome']}: "
        f"{completion['summary']['n_output_positions']} output positions"
    )


def command_workbench_plot(args):
    summary_dir = Path(args.summary_dir)
    summary_paths = {
        chromosome: summary_dir / f"chr{chromosome}.gamma_smc.tsv"
        for chromosome in args.chromosomes
    }
    result = plot_workbench_population(
        summary_paths,
        population=args.population,
        output_dir=args.output_dir,
        threshold_years=args.threshold_years,
        whole_genome=args.whole_genome,
        top_n=args.top_n,
        signal_fraction=args.signal_fraction,
    )
    scope = (
        "whole genome" if result["whole_genome_complete"] else "requested chromosomes"
    )
    action = "reused" if result.get("reused") else "plotted"
    print(
        f"{action} {result['population']} {scope} in {Path(args.output_dir).resolve()}"
    )


def command_workbench_report(args):
    result = summarize_workbench_run(
        args.results_root,
        populations=args.populations,
        chromosomes=args.chromosomes,
        output_dir=args.output_dir,
        threshold_years=args.threshold_years,
        signal_fraction=args.signal_fraction,
        whole_genome=args.whole_genome,
    )
    action = "reused" if result.get("reused") else "wrote"
    print(f"{action} Workbench report: {Path(args.output_dir).resolve()}")
    for population in result["populations"]:
        print(
            f"  {population}: {result['population_region_counts'][population]} "
            "candidate region(s)"
        )
    print(f"  total: {result['total_regions']} candidate region(s)")


def command_workbench_regions(args):
    regions = build_candidate_regions(
        args.summary,
        population=args.population,
        chromosome=args.chromosome,
        sequence_length=args.sequence_length,
        output_path=args.output,
        positions_path=args.positions_output,
        threshold_years=args.threshold_years,
        minimum_fraction=args.minimum_fraction,
        merge_gap=args.merge_gap,
        stride=args.output_at_stride,
        profile_half_width=args.profile_half_width,
    )
    print(
        f"found {len(regions)} {args.population.upper()} chr{args.chromosome} "
        f"regions above {args.minimum_fraction:.3g}"
    )


def command_workbench_candidates(args):
    result = analyze_candidate_regions(
        summary_path=args.summary,
        regions_path=args.regions,
        raw_path=args.raw_posteriors,
        pairs_manifest_path=args.pairs_manifest,
        bcf_path=args.input,
        sample_list_path=args.sample_list,
        bcftools=args.bcftools,
        contig=args.contig,
        population=args.population,
        chromosome=args.chromosome,
        output_dir=args.output_dir,
        mutation_rate=args.mutation_rate,
        threshold_years=args.threshold_years,
        profile_half_width=args.profile_half_width,
        variant_half_width=args.variant_half_width,
        minimum_genotype_pairs=args.minimum_genotype_pairs,
        minimum_fraction=args.minimum_fraction,
        merge_gap=args.merge_gap,
    )
    print(
        f"analyzed {result['n_regions']} candidate regions; "
        f"{result['n_regions_with_representative_variant']} have a representative variant"
    )


def command_workbench_candidates_empty(args):
    write_empty_candidate_analysis(
        regions_path=args.regions,
        population=args.population,
        chromosome=args.chromosome,
        pairs_manifest_path=args.pairs_manifest,
        output_dir=args.output_dir,
    )
    print(f"no {args.population.upper()} chr{args.chromosome} candidate regions")


def command_native_study(args):
    run_native_stride_study(
        args.source_dir,
        args.output_dir,
        executable=args.executable,
        neutral_replicates=args.neutral_replicates,
        stride=args.output_at_stride,
        cache_size=args.cache_size,
        threads=args.threads,
        keep_vcfs=args.keep_vcfs,
        workers=args.workers,
        recent_call=args.recent_call,
        comparison_recent_call=args.compare_recent_call,
        calibration_statistic=args.calibration_statistic,
    )


def command_finalize_container_study(args):
    finalize_container_stride_study(
        args.source_dir,
        args.output_dir,
        stride=args.output_at_stride,
        workflow_elapsed_seconds=args.workflow_elapsed_seconds,
        neutral_profiles_path=args.neutral_profiles,
    )


def command_prepare_high_af_selected(args):
    prepare_high_af_selected(
        args.null_truth_dir,
        args.output_dir,
        executable=args.slim,
        minimum_population_af=args.minimum_population_af,
        selection_coefficient=args.selection_coefficient,
        workers=args.workers,
        max_attempts=args.max_attempts,
        seed=args.seed,
        selected_attempt=args.selected_attempt,
        stride=args.output_at_stride,
    )


def command_finalize_high_af_selected(args):
    finalize_high_af_selected(
        args.source_dir,
        args.neutral_decoded_profiles,
        stride=args.output_at_stride,
    )


def command_evaluate(args):
    truth = sorted(Path(args.truth_dir).glob("*.tsv"))
    decoded = sorted(Path(args.decoded_dir).glob("*.tsv"))
    truth_by_name = {path.name: path for path in truth}
    decoded_by_name = {path.name: path for path in decoded}
    names = sorted(truth_by_name.keys() & decoded_by_name.keys())
    if not names:
        raise FileNotFoundError("no identically named truth and decoded TSVs")
    evaluate_pairs(
        [(truth_by_name[name], decoded_by_name[name]) for name in names],
        args.output_dir,
    )


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
    sim = commands.add_parser(
        "simulate", help="simulate neutral regions with StandardCoalescent"
    )
    sim.add_argument("--output-dir", required=True)
    sim.add_argument("--replicates", type=int, default=1000)
    sim.add_argument("--diploids", type=int, default=2000)
    sim.add_argument("--length", type=int, default=1_000_000)
    sim.add_argument("--ne", type=float, default=10_000)
    sim.add_argument("--mutation-rate", type=float, default=DEFAULT_MUTATION_RATE)
    sim.add_argument("--recombination-rate", type=float, default=1e-8)
    sim.add_argument("--mutation-map")
    sim.add_argument("--recombination-map")
    sim.add_argument(
        "--histories", help="TSV with draw,time_generations,ne (e.g. PHLASH MVN draws)"
    )
    sim.add_argument("--threshold-years", type=float, default=4500)
    sim.add_argument("--generation-time", type=float, default=DEFAULT_GENERATION_TIME)
    sim.add_argument("--seed", type=int, default=1729)
    sim.add_argument("--save-trees", action="store_true")
    sim.add_argument(
        "--workers", type=int, default=1, help="independent simulation processes"
    )
    sim.set_defaults(func=command_simulate)

    cal = commands.add_parser(
        "calibrate", help="calculate site-level simulation p-values"
    )
    cal.add_argument("--observed", required=True)
    cal.add_argument("--sim-glob", required=True)
    cal.add_argument("--observed-length", type=float, required=True)
    cal.add_argument("--simulation-length", type=float, required=True)
    cal.add_argument("--match", choices=["nearest", "region_max"], default="nearest")
    cal.add_argument("--output", required=True)
    cal.set_defaults(func=command_calibrate)

    validate = commands.add_parser(
        "validate-null", help="leave-one-replicate-out p-value calibration"
    )
    validate.add_argument("--sim-glob", required=True)
    validate.add_argument("--sequence-length", type=float, required=True)
    validate.add_argument("--relative-position", type=float, default=0.5)
    validate.add_argument("--tie-seed", type=int, default=8675309)
    validate.add_argument("--output-dir", required=True)
    validate.set_defaults(func=command_validate_null)

    convert = commands.add_parser(
        "convert-ts", help="convert .trees/.tsz to diploid VCF"
    )
    convert.add_argument("--input", required=True)
    convert.add_argument("--output", required=True)
    convert.set_defaults(func=command_convert)

    decode = commands.add_parser(
        "decode", help="run a streaming Gamma-SMC recent-coalescence summary"
    )
    decode.add_argument(
        "--executable",
        default=os.environ.get("GAMMA_SMC_BIN", "gamma_smc"),
        help="Gamma-SMC binary (default: GAMMA_SMC_BIN or gamma_smc on PATH)",
    )
    decode.add_argument("--input", required=True)
    decode.add_argument(
        "--input-format", choices=["auto", "vcf", "trees", "tsz"], default="auto"
    )
    decode.add_argument("--output", required=True)
    decode.add_argument(
        "--raw-output", help="raw alpha/beta posteriors; 8 bytes per pair per position"
    )
    decode.add_argument(
        "--bitmatrix",
        help="packed per-pair recent-coalescence calls, one bit per pair/position/threshold",
    )
    decode.add_argument(
        "--mask",
        help="global BED mask; use the same callable-region policy in data and simulations",
    )
    decode.add_argument("--masks-per-sample", help="sample-to-BED TSV")
    decode.add_argument(
        "--samples",
        help="one sample ID per line; IDs must be present in the VCF/BCF header",
    )
    decode.add_argument(
        "--output-positions-file",
        help="one exact 0-based posterior output position per line",
    )
    decode.add_argument(
        "--theta",
        type=float,
        default=DEFAULT_SCALED_MUTATION_RATE,
        help="scaled mutation rate; fixed by default so coalescent-time units, "
        "and therefore P(T<t), are comparable across datasets",
    )
    decode.add_argument(
        "--rho-over-theta",
        type=float,
        default=DEFAULT_RECOMBINATION_TO_MUTATION_RATIO,
        help="0.8 with the default theta gives the reference rho = 0.0006",
    )
    decode.add_argument("--mutation-rate", type=float, default=DEFAULT_MUTATION_RATE)
    decode.add_argument(
        "--threshold-years",
        type=float,
        nargs="+",
        default=[4500],
        help="one or more thresholds, e.g. --threshold-years 4500 10000",
    )
    decode.add_argument(
        "--generation-time", type=float, default=DEFAULT_GENERATION_TIME
    )
    decode.add_argument(
        "--output-at-stride",
        type=int,
        default=DEFAULT_OUTPUT_STRIDE,
        help="10 kb balances scan resolution, output size, and Workbench runtime",
    )
    decode.add_argument("--no-output-at-hets", action="store_true")
    decode.add_argument(
        "--n-random-pairs",
        type=int,
        default=0,
        help="sample this many haplotype pairs uniformly instead of one per diploid",
    )
    decode.add_argument("--pairs-seed", type=int, default=1729)
    decode.add_argument(
        "--pairs-file",
        help="explicit pair list: two 0-based haplotype indices per line",
    )
    decode.add_argument(
        "--exclude-within",
        action="store_true",
        help="drop within-individual pairs when sampling at random",
    )
    decode.add_argument(
        "--recent-call",
        choices=["median", "mean", "prob"],
        default="median",
        help="per-pair call rule: posterior median below the threshold (default), posterior mean, or P>=p",
    )
    decode.add_argument("--recent-call-probability", type=float, default=0.5)
    decode.add_argument(
        "--threads", type=int, default=0, help="0 uses every available core"
    )
    decode.add_argument(
        "--cache-size",
        type=int,
        default=DEFAULT_CACHE_SIZE,
        help="maximum transition-cache segment in bp; 1 kb is the validated default",
    )
    decode.add_argument("--pair-block", type=int, default=256)
    decode.add_argument(
        "--pairs-manifest",
        help="write the decoded pair list here (reusable as --pairs-file); "
        "derived next to the output automatically when sampling at random",
    )
    decode.add_argument(
        "--exp10",
        choices=["accurate", "fast"],
        default="accurate",
        help="fast reproduces upstream's 10^x approximation, which carries -3.9%%..+2.0%% relative error",
    )
    decode.add_argument(
        "--backward-alignment",
        choices=["fixed", "legacy"],
        default="fixed",
        help="legacy reproduces upstream's one-output-position shift of the backward message",
    )
    decode.set_defaults(func=command_decode)

    workbench_samples = commands.add_parser(
        "workbench-samples",
        help="build one population's BCF-ordered, QC-filtered Workbench sample list",
    )
    workbench_samples.add_argument("--population", type=str.upper, required=True)
    workbench_samples.add_argument("--ancestry", required=True)
    workbench_samples.add_argument("--qc-exclusions", required=True)
    workbench_samples.add_argument("--relatedness-exclusions", required=True)
    workbench_samples.add_argument("--bcf-samples", required=True)
    workbench_samples.add_argument("--output", required=True)
    workbench_samples.add_argument("--audit-output", required=True)
    workbench_samples.set_defaults(func=command_workbench_samples)

    workbench_mask = commands.add_parser(
        "workbench-mask",
        help="complement an exclusion hard mask into Gamma-SMC callable intervals",
    )
    workbench_mask.add_argument("--hardmask", required=True)
    workbench_mask.add_argument("--contig", required=True)
    workbench_mask.add_argument("--sequence-length", type=int, required=True)
    workbench_mask.add_argument("--output", required=True)
    workbench_mask.add_argument("--audit-output", required=True)
    workbench_mask.set_defaults(func=command_workbench_mask)

    workbench_validate = commands.add_parser(
        "workbench-validate",
        help="write or verify a restart-safe Workbench chromosome completion record",
    )
    workbench_validate.add_argument("--population", type=str.upper, required=True)
    workbench_validate.add_argument(
        "--chromosome", type=int, choices=AUTOSOMES, required=True
    )
    workbench_validate.add_argument("--input-uri", required=True)
    workbench_validate.add_argument("--input-fingerprint", required=True)
    workbench_validate.add_argument("--local-input", required=True)
    workbench_validate.add_argument("--index-uri", required=True)
    workbench_validate.add_argument("--index-fingerprint", required=True)
    workbench_validate.add_argument("--local-index", required=True)
    workbench_validate.add_argument("--summary", required=True)
    workbench_validate.add_argument("--run-json", required=True)
    workbench_validate.add_argument("--pairs-manifest", required=True)
    workbench_validate.add_argument("--bitmatrix")
    workbench_validate.add_argument("--candidate-regions")
    workbench_validate.add_argument("--candidate-manifest")
    workbench_validate.add_argument("--sample-list", required=True)
    workbench_validate.add_argument("--sample-audit", required=True)
    workbench_validate.add_argument("--ancestry-uri", required=True)
    workbench_validate.add_argument("--ancestry-fingerprint", required=True)
    workbench_validate.add_argument("--local-ancestry", required=True)
    workbench_validate.add_argument("--qc-exclusions-uri", required=True)
    workbench_validate.add_argument("--qc-exclusions-fingerprint", required=True)
    workbench_validate.add_argument("--local-qc-exclusions", required=True)
    workbench_validate.add_argument("--relatedness-exclusions-uri", required=True)
    workbench_validate.add_argument(
        "--relatedness-exclusions-fingerprint", required=True
    )
    workbench_validate.add_argument("--local-relatedness-exclusions", required=True)
    workbench_validate.add_argument("--completion", required=True)
    workbench_validate.add_argument("--mask-uri")
    workbench_validate.add_argument("--mask-fingerprint")
    workbench_validate.add_argument("--local-mask-source")
    workbench_validate.add_argument("--local-mask")
    workbench_validate.add_argument("--mask-audit")
    workbench_validate.add_argument(
        "--theta", type=float, default=DEFAULT_SCALED_MUTATION_RATE
    )
    workbench_validate.add_argument(
        "--rho-over-theta",
        type=float,
        default=DEFAULT_RECOMBINATION_TO_MUTATION_RATIO,
    )
    workbench_validate.add_argument(
        "--mutation-rate", type=float, default=DEFAULT_MUTATION_RATE
    )
    workbench_validate.add_argument(
        "--generation-time", type=float, default=DEFAULT_GENERATION_TIME
    )
    workbench_validate.add_argument("--threshold-years", type=float, default=4500)
    workbench_validate.add_argument(
        "--recent-call", choices=["mean", "median"], default="mean"
    )
    workbench_validate.add_argument(
        "--output-at-stride", type=int, default=DEFAULT_OUTPUT_STRIDE
    )
    workbench_validate.add_argument(
        "--cache-size", type=int, default=DEFAULT_CACHE_SIZE
    )
    workbench_validate.add_argument("--threads", type=int, default=12)
    workbench_validate.add_argument("--pair-block", type=int, default=256)
    workbench_validate.add_argument("--n-random-pairs", type=int, default=0)
    workbench_validate.add_argument("--pairs-seed", type=int, default=1729)
    workbench_validate.add_argument("--exclude-within", action="store_true")
    workbench_validate.add_argument("--signal-fraction", type=float, default=0.05)
    workbench_validate.add_argument("--merge-gap", type=int, default=20_000)
    workbench_validate.add_argument("--profile-half-width", type=int, default=500_000)
    workbench_validate.add_argument("--variant-half-width", type=int, default=100_000)
    workbench_validate.add_argument("--minimum-genotype-pairs", type=int, default=20)
    workbench_validate.add_argument(
        "--exp10", choices=["accurate", "fast"], default="accurate"
    )
    workbench_validate.add_argument(
        "--backward-alignment", choices=["fixed", "legacy"], default="fixed"
    )
    workbench_validate.add_argument("--code-commit", required=True)
    workbench_validate.add_argument("--check-only", action="store_true")
    workbench_validate.set_defaults(func=command_workbench_validate)

    workbench_plot = commands.add_parser(
        "workbench-plot",
        help="make separate chromosome and whole-genome plots for one population",
    )
    workbench_plot.add_argument("--population", type=str.upper, required=True)
    workbench_plot.add_argument("--summary-dir", required=True)
    workbench_plot.add_argument(
        "--chromosomes", type=int, choices=AUTOSOMES, nargs="+", required=True
    )
    workbench_plot.add_argument("--output-dir", required=True)
    workbench_plot.add_argument("--threshold-years", type=float, default=4500)
    workbench_plot.add_argument("--whole-genome", action="store_true")
    workbench_plot.add_argument("--top-n", type=int, default=100)
    workbench_plot.add_argument("--signal-fraction", type=float, default=0.05)
    workbench_plot.set_defaults(func=command_workbench_plot)

    workbench_report = commands.add_parser(
        "workbench-report",
        help="summarize population regions and make a combined genome-wide plot",
    )
    workbench_report.add_argument("--results-root", required=True)
    workbench_report.add_argument(
        "--populations", type=str.upper, choices=POPULATIONS, nargs="+", required=True
    )
    workbench_report.add_argument(
        "--chromosomes", type=int, choices=AUTOSOMES, nargs="+", required=True
    )
    workbench_report.add_argument("--output-dir", required=True)
    workbench_report.add_argument("--threshold-years", type=float, default=4500)
    workbench_report.add_argument("--signal-fraction", type=float, default=0.05)
    workbench_report.add_argument("--whole-genome", action="store_true")
    workbench_report.set_defaults(func=command_workbench_report)

    workbench_regions = commands.add_parser(
        "workbench-regions",
        help="merge >5%% recent-coalescence windows and select +/-500 kb positions",
    )
    workbench_regions.add_argument("--population", type=str.upper, required=True)
    workbench_regions.add_argument(
        "--chromosome", type=int, choices=AUTOSOMES, required=True
    )
    workbench_regions.add_argument("--summary", required=True)
    workbench_regions.add_argument("--sequence-length", type=int, required=True)
    workbench_regions.add_argument("--output", required=True)
    workbench_regions.add_argument("--positions-output", required=True)
    workbench_regions.add_argument("--threshold-years", type=float, default=4500)
    workbench_regions.add_argument("--minimum-fraction", type=float, default=0.05)
    workbench_regions.add_argument("--merge-gap", type=int, default=20_000)
    workbench_regions.add_argument(
        "--output-at-stride", type=int, default=DEFAULT_OUTPUT_STRIDE
    )
    workbench_regions.add_argument("--profile-half-width", type=int, default=500_000)
    workbench_regions.set_defaults(func=command_workbench_regions)

    workbench_candidates = commands.add_parser(
        "workbench-candidates",
        help="rank nearby variants and plot manifest-ordered ref/ref and alt/alt pairs",
    )
    workbench_candidates.add_argument("--population", type=str.upper, required=True)
    workbench_candidates.add_argument(
        "--chromosome", type=int, choices=AUTOSOMES, required=True
    )
    workbench_candidates.add_argument("--summary", required=True)
    workbench_candidates.add_argument("--regions", required=True)
    workbench_candidates.add_argument("--raw-posteriors", required=True)
    workbench_candidates.add_argument("--pairs-manifest", required=True)
    workbench_candidates.add_argument("--input", required=True)
    workbench_candidates.add_argument("--sample-list", required=True)
    workbench_candidates.add_argument("--bcftools", required=True)
    workbench_candidates.add_argument("--contig", required=True)
    workbench_candidates.add_argument("--output-dir", required=True)
    workbench_candidates.add_argument(
        "--mutation-rate", type=float, default=DEFAULT_MUTATION_RATE
    )
    workbench_candidates.add_argument("--threshold-years", type=float, default=4500)
    workbench_candidates.add_argument("--profile-half-width", type=int, default=500_000)
    workbench_candidates.add_argument("--variant-half-width", type=int, default=100_000)
    workbench_candidates.add_argument("--minimum-genotype-pairs", type=int, default=20)
    workbench_candidates.add_argument("--minimum-fraction", type=float, default=0.05)
    workbench_candidates.add_argument("--merge-gap", type=int, default=20_000)
    workbench_candidates.set_defaults(func=command_workbench_candidates)

    workbench_empty = commands.add_parser(
        "workbench-candidates-empty",
        help="write a complete empty candidate-analysis manifest",
    )
    workbench_empty.add_argument("--population", type=str.upper, required=True)
    workbench_empty.add_argument(
        "--chromosome", type=int, choices=AUTOSOMES, required=True
    )
    workbench_empty.add_argument("--regions", required=True)
    workbench_empty.add_argument("--pairs-manifest", required=True)
    workbench_empty.add_argument("--output-dir", required=True)
    workbench_empty.set_defaults(func=command_workbench_candidates_empty)

    bits = commands.add_parser(
        "bitmatrix-summary",
        help="per-position recent-coalescence counts from a packed bit matrix",
    )
    bits.add_argument(
        "--input",
        required=True,
        help="path to the .bits file (its .meta must sit alongside)",
    )
    bits.add_argument("--output", required=True)
    bits.add_argument(
        "--pairs",
        help="optional file of 0-based pair indices to restrict the counts to",
    )
    bits.set_defaults(func=command_bitmatrix_summary)

    container = commands.add_parser(
        "decode-container",
        help="decode a VCF with the official Gamma-SMC v0.2 container",
    )
    container.add_argument("--input", required=True)
    container.add_argument("--output", required=True)
    container.add_argument("--theta", type=float, default=DEFAULT_SCALED_MUTATION_RATE)
    container.add_argument(
        "--rho-over-theta",
        type=float,
        default=DEFAULT_RECOMBINATION_TO_MUTATION_RATIO,
    )
    container.add_argument("--mutation-rate", type=float, default=DEFAULT_MUTATION_RATE)
    container.add_argument("--threshold-years", type=float, default=4500)
    container.add_argument(
        "--generation-time", type=float, default=DEFAULT_GENERATION_TIME
    )
    container.add_argument(
        "--output-at-stride", type=int, default=DEFAULT_OUTPUT_STRIDE
    )
    container.add_argument(
        "--runtime",
        choices=["auto", "apptainer", "singularity", "docker"],
        default="auto",
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
    study.add_argument("--output-at-stride", type=int, default=DEFAULT_OUTPUT_STRIDE)
    study.add_argument(
        "--runtime",
        choices=["auto", "apptainer", "singularity", "docker"],
        default="auto",
    )
    study.add_argument("--image", default=DEFAULT_IMAGE)
    study.add_argument("--keep-vcfs", action="store_true")
    study.add_argument("--workers", type=int, default=1)
    study.set_defaults(func=command_container_study)

    native_study = commands.add_parser(
        "run-native-study",
        help="decode the retained selected simulation and matched nulls with the optimized binary",
    )
    native_study.add_argument("--source-dir", required=True)
    native_study.add_argument("--output-dir", required=True)
    native_study.add_argument(
        "--executable",
        default=os.environ.get("GAMMA_SMC_BIN", "gamma_smc"),
    )
    native_study.add_argument("--neutral-replicates", type=int, default=100)
    native_study.add_argument(
        "--output-at-stride", type=int, default=DEFAULT_OUTPUT_STRIDE
    )
    native_study.add_argument("--cache-size", type=int, default=DEFAULT_CACHE_SIZE)
    native_study.add_argument("--threads", type=int, default=1)
    native_study.add_argument("--keep-vcfs", action="store_true")
    native_study.add_argument("--workers", type=int, default=12)
    native_study.add_argument(
        "--recent-call",
        choices=["mean", "median"],
        default="median",
        help="per-pair posterior summary used for n_recent/frac_recent",
    )
    native_study.add_argument(
        "--compare-recent-call",
        choices=["mean", "median"],
        help="also decode the same selected and null simulations with this call rule",
    )
    native_study.add_argument(
        "--calibration-statistic",
        choices=["mean-posterior-probability", "called-fraction"],
        default="mean-posterior-probability",
        help="statistic calibrated against the neutral simulations",
    )
    native_study.set_defaults(func=command_native_study)

    finalize = commands.add_parser(
        "finalize-container-study",
        help="calculate plots and p-values from existing selected/null decoded profiles",
    )
    finalize.add_argument("--source-dir", required=True)
    finalize.add_argument("--output-dir", required=True)
    finalize.add_argument("--output-at-stride", type=int, default=DEFAULT_OUTPUT_STRIDE)
    finalize.add_argument("--workflow-elapsed-seconds", type=float)
    finalize.add_argument("--neutral-profiles")
    finalize.set_defaults(func=command_finalize_container_study)

    high_af = commands.add_parser(
        "prepare-high-af-selected",
        help="rejection-sample one high-frequency selected replicate and prepare its VCF",
    )
    high_af.add_argument("--null-truth-dir", required=True)
    high_af.add_argument("--output-dir", required=True)
    high_af.add_argument("--slim", help="SLiM executable; otherwise use SLIM_BIN/PATH")
    high_af.add_argument("--minimum-population-af", type=float, default=0.30)
    high_af.add_argument("--selection-coefficient", type=float, default=0.05)
    high_af.add_argument("--workers", type=int, default=12)
    high_af.add_argument("--max-attempts", type=int, default=2_000)
    high_af.add_argument("--seed", type=int, default=910_241)
    high_af.add_argument(
        "--selected-attempt",
        type=int,
        help="materialize and validate one known deterministic attempt instead of re-screening",
    )
    high_af.add_argument("--output-at-stride", type=int, default=DEFAULT_OUTPUT_STRIDE)
    high_af.set_defaults(func=command_prepare_high_af_selected)

    finish_high_af = commands.add_parser(
        "finalize-high-af-selected",
        help="calibrate one selected container decode against an existing decoded null",
    )
    finish_high_af.add_argument("--source-dir", required=True)
    finish_high_af.add_argument("--neutral-decoded-profiles", required=True)
    finish_high_af.add_argument(
        "--output-at-stride", type=int, default=DEFAULT_OUTPUT_STRIDE
    )
    finish_high_af.set_defaults(func=command_finalize_high_af_selected)

    evaluate = commands.add_parser(
        "evaluate-decoder", help="compare decoded simulations with tree-sequence truth"
    )
    evaluate.add_argument("--truth-dir", required=True)
    evaluate.add_argument("--decoded-dir", required=True)
    evaluate.add_argument("--output-dir", required=True)
    evaluate.set_defaults(func=command_evaluate)

    truth_plot = commands.add_parser(
        "plot-truth", help="plot true recent fraction against simulated mean TMRCA"
    )
    truth_plot.add_argument("--summary-dir", required=True)
    truth_plot.add_argument("--sequence-length", type=float, required=True)
    truth_plot.add_argument("--ne", type=float, required=True)
    truth_plot.add_argument("--threshold-years", type=float, default=4500)
    truth_plot.add_argument(
        "--generation-time", type=float, default=DEFAULT_GENERATION_TIME
    )
    truth_plot.add_argument("--relative-position", type=float, default=0.5)
    truth_plot.add_argument("--output-dir", required=True)
    truth_plot.set_defaults(func=command_plot_truth)

    sweep = commands.add_parser(
        "validate-sweep",
        help="validate recent-coalescence power with a SLiM hard sweep",
    )
    sweep.add_argument("--output-dir", required=True)
    sweep.add_argument("--slim", help="SLiM executable; otherwise use SLIM_BIN/PATH")
    sweep.add_argument("--population-size", type=int, default=200)
    sweep.add_argument("--sequence-length", type=int, default=100_000)
    sweep.add_argument("--selection-coefficient", type=float, default=0.5)
    sweep.add_argument("--recombination-rate", type=float, default=1e-7)
    sweep.add_argument("--threshold-years", type=float, default=4500)
    sweep.add_argument("--generation-time", type=float, default=DEFAULT_GENERATION_TIME)
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
    recent.add_argument("--mutation-rate", type=float, default=DEFAULT_MUTATION_RATE)
    recent.add_argument("--recombination-rate", type=float, default=1e-8)
    recent.add_argument("--neutral-replicates", type=int, default=100)
    recent.add_argument("--selected-replicates", type=int, default=1)
    recent.add_argument(
        "--workers", type=int, default=1, help="parallel selected SLiM trajectories"
    )
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
    two_epoch.add_argument(
        "--generation-time-years", type=float, default=DEFAULT_GENERATION_TIME
    )
    two_epoch.add_argument("--selection-coefficient", type=float, default=0.05)
    two_epoch.add_argument("--mutation-rate", type=float, default=DEFAULT_MUTATION_RATE)
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
    spatial.add_argument("--slim", help="SLiM executable; otherwise use SLIM_BIN/PATH")
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
