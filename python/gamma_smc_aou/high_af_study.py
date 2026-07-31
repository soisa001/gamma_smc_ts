from __future__ import annotations

import json
import os
import shutil
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from tempfile import TemporaryDirectory
from time import perf_counter, sleep

import matplotlib
matplotlib.use("Agg")
import msprime
import numpy as np
import pandas as pd

from .defaults import DEFAULT_OUTPUT_STRIDE
import pyslim
import tskit

from .carrier_profiles import (
    _position_grids,
    pair_tmrca_profile_by_focal_copy,
    plot_carrier_profile_figure,
)
from .container_study import finalize_container_stride_study
from .selection import (
    run_slim_recent_sweep,
    within_individual_tmrca_details,
    within_individual_tmrca_grid,
)
from .spatial_scan import _plot_observed_profile
from .tree_sequence import tree_sequence_to_vcf
from .two_epoch import (
    _plot_allele_frequency_trajectory,
    _plot_demography_and_variant,
    _plot_null_calibration,
    two_epoch_recent_probability,
)


def focal_carrier_pair_table(ts, focal_carrier_counts: np.ndarray) -> pd.DataFrame:
    """Map saved focal genotypes to the VCF haplotype-pair ordering.

    ``tree_sequence_to_vcf`` writes diploid individuals in tree-sequence
    individual order. Gamma-SMC numbers the resulting haplotypes consecutively,
    so VCF diploid row ``k`` is pair ``(2*k, 2*k + 1)``.
    """
    counts = np.asarray(focal_carrier_counts, dtype=np.int8)
    sample_nodes = set(ts.samples())
    rows = []
    for individual in ts.individuals():
        nodes = [node for node in individual.nodes if node in sample_nodes]
        if not nodes:
            continue
        vcf_diploid_index = len(rows)
        if len(nodes) != 2:
            raise ValueError(
                f"individual {individual.id} has {len(nodes)} sample nodes; expected 2"
            )
        if individual.id >= len(counts):
            raise ValueError("focal carrier counts do not cover every VCF individual")
        copies = int(counts[individual.id])
        if copies not in {0, 1, 2}:
            raise ValueError("focal carrier counts must contain only 0, 1, or 2")
        rows.append({
            "vcf_diploid_index": vcf_diploid_index,
            "tree_sequence_individual_id": individual.id,
            "sample_node_0": nodes[0],
            "sample_node_1": nodes[1],
            "gamma_smc_haplotype_0": 2 * vcf_diploid_index,
            "gamma_smc_haplotype_1": 2 * vcf_diploid_index + 1,
            "focal_carrier_copies": copies,
            "focal_genotype_class": (
                "hom_ref"
                if copies == 0
                else "heterozygous"
                if copies == 1
                else "hom_alt"
            ),
        })
    table = pd.DataFrame(rows)
    if len(table) != len(counts):
        raise ValueError(
            f"VCF contains {len(table)} diploids but there are {len(counts)} carrier counts"
        )
    return table


def _simulate_high_af_attempt(task: dict) -> dict:
    """Run one independent trajectory in a worker process."""
    attempt = int(task["attempt"])
    run_seed = int(task["seed"] + attempt * 10)
    tree_path = Path(task["temporary"]) / f"selected_attempt_{attempt:04d}.trees"
    ts, run = run_slim_recent_sweep(
        tree_path,
        executable=task["executable"],
        population_size=task["ancestral_population_size"],
        present_population_size=task["present_population_size"],
        size_change_generations_ago=task["size_change_generations_ago"],
        sample_diploids=task["sample_diploids"],
        sequence_length=task["sequence_length"],
        sweep_position=task["center"],
        selection_coefficient=task["selection_coefficient"],
        age_generations=task["variant_age_generations"],
        mutation_rate=task["mutation_rate"],
        recombination_rate=task["recombination_rate"],
        seed=run_seed,
        capture_focal_genotypes=not task["screen_only"],
        initial_annotated_path=task["initial_annotated_path"],
        screen_only=task["screen_only"],
    )
    population_af = float(run["realized_population_allele_frequency"])
    if task["screen_only"]:
        return {
            "row": {
                "attempt": attempt,
                "seed": run_seed,
                "accepted": bool(
                    population_af >= task["minimum_population_af"]
                    and run["focal_allele_outcome"] == "segregating"
                ),
                "focal_allele_outcome": run["focal_allele_outcome"],
                "population_allele_frequency": population_af,
                "ancestry_seconds": float(run["ancestry_seconds"]),
                "slim_seconds": float(run["slim_seconds"]),
            },
            "trajectory": run["allele_frequency_trajectory"],
        }
    counts = np.asarray(run["focal_carrier_counts"], dtype=np.int8)
    details = within_individual_tmrca_details(
        ts,
        task["center"],
        task["variant_age_generations"],
        focal_carrier_counts=counts,
    )
    row = {
        "attempt": attempt,
        "seed": run_seed,
        "accepted": bool(
            population_af >= task["minimum_population_af"]
            and np.count_nonzero(counts == 0) >= 2
            and np.count_nonzero(counts == 2) >= 2
        ),
        "focal_allele_outcome": run["focal_allele_outcome"],
        "population_allele_frequency": population_af,
        "sample_allele_frequency": float(run["sample_focal_allele_frequency"]),
        "n_hom_ref_pairs": int(np.count_nonzero(counts == 0)),
        "n_heterozygous_pairs": int(np.count_nonzero(counts == 1)),
        "n_hom_alt_pairs": int(np.count_nonzero(counts == 2)),
        "center_fraction_recent": float(details["tmrca_lt_threshold"].mean()),
        "center_n_pairs_recent": int(details["tmrca_lt_threshold"].sum()),
        "center_mean_tmrca_generations": float(details["tmrca_generations"].mean()),
        "ancestry_seconds": float(run["ancestry_seconds"]),
        "slim_seconds": float(run["slim_seconds"]),
    }
    return {
        "row": row,
        "tree_path": str(tree_path),
        "carrier_counts": counts,
        "trajectory": run["allele_frequency_trajectory"],
    }


def _check_null_compatibility(null_metrics: dict, design: dict) -> None:
    checks = {
        "sample_diploids": design["sample_diploids"],
        "sequence_length": design["sequence_length"],
        "variant_age_generations": design["variant_age_generations"],
        "generation_time_years": design["generation_time_years"],
        "mutation_rate": design["mutation_rate"],
        "recombination_rate": design["recombination_rate"],
    }
    for key, expected in checks.items():
        if not np.isclose(
            float(null_metrics[key]), float(expected), rtol=1e-12, atol=0.0
        ):
            raise ValueError(
                f"null design mismatch for {key}: {null_metrics[key]} != {expected}"
            )
    demography = null_metrics["demography"]
    for key, expected in {
        "ancestral_population_size": design["ancestral_population_size"],
        "present_population_size": design["present_population_size"],
        "size_change_generations_ago": design["size_change_generations_ago"],
    }.items():
        if int(demography[key]) != int(expected):
            raise ValueError(
                f"null demography mismatch for {key}: {demography[key]} != {expected}"
            )


def _resolve_neutral_statistics(
    null_truth_dir: Path,
    null_metrics: dict,
) -> Path:
    """Find the 100-replicate truth table, following recorded provenance."""
    local = null_truth_dir / "neutral_statistics.tsv"
    if local.is_file():
        return local
    for key in ("neutral_statistics_source", "neutral_truth_source"):
        recorded = null_metrics.get(key)
        if not recorded:
            continue
        source = (null_truth_dir / str(recorded)).resolve()
        candidate = (
            source
            if source.name == "neutral_statistics.tsv"
            else source / "neutral_statistics.tsv"
        )
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(
        f"{local} is absent and metrics.json does not resolve to a "
        "neutral_statistics.tsv provenance source"
    )


def _replace_with_retry(
    source: Path,
    destination: Path,
    *,
    attempts: int = 20,
) -> None:
    """Atomically replace a file despite brief Windows/OneDrive read locks."""
    for attempt in range(attempts):
        try:
            os.replace(source, destination)
            return
        except PermissionError:
            if attempt + 1 == attempts:
                raise
            sleep(min(0.05 * (attempt + 1), 0.5))


def prepare_high_af_selected(
    null_truth_dir: str | Path,
    output_dir: str | Path,
    *,
    executable: str | Path | None = None,
    minimum_population_af: float = 0.30,
    selection_coefficient: float = 0.05,
    workers: int = 12,
    max_attempts: int = 2_000,
    seed: int = 910_241,
    selected_attempt: int | None = None,
    stride: int = DEFAULT_OUTPUT_STRIDE,
    full_step: int = 50_000,
    zoom_half_width: int = 500_000,
    zoom_step: int = 5_000,
) -> dict:
    """Rejection-sample one high-frequency sweep while reusing a fixed null."""
    if not 0 < minimum_population_af <= 1:
        raise ValueError("minimum_population_af must be in (0, 1]")
    if workers < 1 or max_attempts < 1 or stride < 1:
        raise ValueError("workers, max_attempts, and stride must be positive")
    if selected_attempt is not None and selected_attempt < 0:
        raise ValueError("selected_attempt must be non-negative")
    null_truth_dir = Path(null_truth_dir).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    with (null_truth_dir / "metrics.json").open(encoding="utf-8") as handle:
        null_metrics = json.load(handle)
    design = {
        "ancestral_population_size": int(
            null_metrics["demography"]["ancestral_population_size"]
        ),
        "present_population_size": int(
            null_metrics["demography"]["present_population_size"]
        ),
        "size_change_generations_ago": int(
            null_metrics["demography"]["size_change_generations_ago"]
        ),
        "sample_diploids": int(null_metrics["sample_diploids"]),
        "sequence_length": int(null_metrics["sequence_length"]),
        "variant_age_generations": int(null_metrics["variant_age_generations"]),
        "generation_time_years": float(null_metrics["generation_time_years"]),
        "mutation_rate": float(null_metrics["mutation_rate"]),
        "recombination_rate": float(null_metrics["recombination_rate"]),
    }
    _check_null_compatibility(null_metrics, design)
    neutral_statistics_path = _resolve_neutral_statistics(
        null_truth_dir,
        null_metrics,
    )
    neutral = pd.read_csv(neutral_statistics_path, sep="\t")
    if len(neutral) != 100:
        raise ValueError(f"expected the saved 100-replicate null; found {len(neutral)}")
    center = design["sequence_length"] // 2
    _, _, carrier_positions = _position_grids(
        design["sequence_length"],
        center,
        full_step=full_step,
        zoom_half_width=zoom_half_width,
        zoom_step=zoom_step,
    )
    started = perf_counter()
    accepted = None
    checkpoint_path = output_dir / "selected_rejection_screen_checkpoint.tsv"
    checkpoint_temporary = checkpoint_path.with_suffix(".tsv.tmp")
    contract_path = output_dir / "selected_rejection_screen_contract.json"
    screen_contract = {
        "schema_version": 1,
        "ancestral_population_size": design["ancestral_population_size"],
        "present_population_size": design["present_population_size"],
        "size_change_generations_ago": design["size_change_generations_ago"],
        "sample_diploids": design["sample_diploids"],
        "sequence_length": design["sequence_length"],
        "variant_age_generations": design["variant_age_generations"],
        "mutation_rate": design["mutation_rate"],
        "recombination_rate": design["recombination_rate"],
        "minimum_population_af": float(minimum_population_af),
        "selection_coefficient": float(selection_coefficient),
        "seed": int(seed),
        "stride": int(stride),
        "null_truth_dir": os.path.relpath(null_truth_dir, output_dir),
        "neutral_statistics_source": os.path.relpath(
            neutral_statistics_path,
            output_dir,
        ),
    }
    attempt_rows: list[dict] = []
    next_attempt = 0
    resume_accepted_attempt: int | None = None
    checkpoint_required = {
        "attempt",
        "seed",
        "accepted",
        "focal_allele_outcome",
        "population_allele_frequency",
    }
    checkpoint_candidates: list[pd.DataFrame] = []
    for candidate_path in (checkpoint_path, checkpoint_temporary):
        if not candidate_path.exists():
            continue
        try:
            candidate = pd.read_csv(candidate_path, sep="\t")
            if checkpoint_required.difference(candidate.columns):
                continue
            candidate.sort_values("attempt", inplace=True)
        except (KeyError, OSError, ValueError, pd.errors.ParserError):
            continue
        candidate_attempts = candidate["attempt"].to_numpy(dtype=int)
        if np.array_equal(
            candidate_attempts,
            np.arange(len(candidate), dtype=int),
        ):
            checkpoint_candidates.append(candidate)
    if checkpoint_candidates:
        if not contract_path.exists():
            raise RuntimeError(
                f"{checkpoint_path} exists without {contract_path}; refusing "
                "to mix an unverified rejection screen"
            )
        with contract_path.open(encoding="utf-8") as handle:
            saved_contract = json.load(handle)
        if saved_contract != screen_contract:
            raise RuntimeError(
                "rejection-screen checkpoint contract does not match this run"
            )
        checkpoint = max(checkpoint_candidates, key=len)
        observed_attempts = checkpoint["attempt"].to_numpy(dtype=int)
        if not np.array_equal(
            observed_attempts,
            np.arange(len(checkpoint), dtype=int),
        ):
            raise RuntimeError(
                "rejection-screen checkpoint attempts are not contiguous from zero"
            )
        expected_seeds = seed + observed_attempts * 10
        if not np.array_equal(
            checkpoint["seed"].to_numpy(dtype=int),
            expected_seeds,
        ):
            raise RuntimeError("rejection-screen checkpoint seeds do not match")
        expected_accepted = (
            checkpoint["population_allele_frequency"].to_numpy(dtype=float)
            >= minimum_population_af
        ) & checkpoint["focal_allele_outcome"].eq("segregating").to_numpy()
        observed_accepted = (
            checkpoint["accepted"]
            .astype(str)
            .str.lower()
            .eq("true")
            .to_numpy()
        )
        if not np.array_equal(observed_accepted, expected_accepted):
            raise RuntimeError(
                "rejection-screen checkpoint acceptance decisions do not match"
            )
        attempt_rows = checkpoint.to_dict("records")
        next_attempt = len(checkpoint)
        if np.any(expected_accepted):
            resume_accepted_attempt = int(observed_attempts[expected_accepted][0])
    else:
        contract_temporary = contract_path.with_suffix(".json.tmp")
        with contract_temporary.open("w", encoding="utf-8") as handle:
            json.dump(screen_contract, handle, indent=2)
            handle.write("\n")
        _replace_with_retry(contract_temporary, contract_path)

    with TemporaryDirectory(prefix="gamma_smc_high_af_") as temporary:
        temporary = Path(temporary)
        initial_annotated_path = temporary / "shared_neutral_ancestry.trees"
        ancestry_started = perf_counter()
        initial = msprime.sim_ancestry(
            samples=[
                msprime.SampleSet(
                    design["ancestral_population_size"], ploidy=2
                )
            ],
            population_size=design["ancestral_population_size"],
            sequence_length=design["sequence_length"],
            recombination_rate=design["recombination_rate"],
            model=msprime.StandardCoalescent(),
            random_seed=seed - 1,
        )
        pyslim.annotate(initial, model_type="WF", tick=1, stage="early").dump(
            initial_annotated_path
        )
        shared_ancestry_seconds = perf_counter() - ancestry_started

        def attempt_task(attempt: int, *, screen_only: bool) -> dict:
            return {
                "attempt": attempt,
                "seed": seed,
                "temporary": str(temporary),
                "executable": None if executable is None else str(executable),
                "ancestral_population_size": design["ancestral_population_size"],
                "present_population_size": design["present_population_size"],
                "size_change_generations_ago": design[
                    "size_change_generations_ago"
                ],
                "sample_diploids": design["sample_diploids"],
                "sequence_length": design["sequence_length"],
                "center": center,
                "selection_coefficient": selection_coefficient,
                "variant_age_generations": design["variant_age_generations"],
                "mutation_rate": design["mutation_rate"],
                "recombination_rate": design["recombination_rate"],
                "minimum_population_af": minimum_population_af,
                "initial_annotated_path": str(initial_annotated_path),
                "screen_only": screen_only,
            }

        if selected_attempt is not None:
            accepted = _simulate_high_af_attempt(
                attempt_task(selected_attempt, screen_only=False)
            )
            attempt_rows.append(accepted["row"])
            if not accepted["row"]["accepted"]:
                raise RuntimeError(
                    f"explicit attempt {selected_attempt} did not reach population "
                    f"AF >= {minimum_population_af:g}"
                )
        else:
            if resume_accepted_attempt is not None:
                accepted = {"row": {"attempt": resume_accepted_attempt}}
            executor = (
                None
                if workers == 1
                else ProcessPoolExecutor(max_workers=workers)
            )
            queued_attempts = workers if executor is None else workers * 10
            try:
                while accepted is None and next_attempt < max_attempts:
                    batch = list(
                        range(
                            next_attempt,
                            min(next_attempt + queued_attempts, max_attempts),
                        )
                    )
                    tasks = [
                        attempt_task(attempt, screen_only=True)
                        for attempt in batch
                    ]
                    results = (
                        [_simulate_high_af_attempt(task) for task in tasks]
                        if executor is None
                        else list(executor.map(_simulate_high_af_attempt, tasks))
                    )
                    for result in results:
                        attempt_rows.append(result["row"])
                        if accepted is None and result["row"]["accepted"]:
                            accepted = result
                    pd.DataFrame(attempt_rows).sort_values("attempt").to_csv(
                        checkpoint_temporary,
                        sep="\t",
                        index=False,
                    )
                    _replace_with_retry(
                        checkpoint_temporary,
                        checkpoint_path,
                    )
                    next_attempt += len(batch)
            finally:
                if executor is not None:
                    executor.shutdown()
            if accepted is None:
                failed = pd.DataFrame(attempt_rows).sort_values("attempt")
                failed["used_for_rejection_decision"] = True
                failed.to_csv(
                    output_dir / "selected_rejection_attempts.tsv",
                    sep="\t",
                    index=False,
                )
                raise RuntimeError(
                    f"no s={selection_coefficient:g} trajectory reached population AF "
                    f">={minimum_population_af:g} in {max_attempts} attempts"
                )
            accepted = _simulate_high_af_attempt(
                attempt_task(int(accepted["row"]["attempt"]), screen_only=False)
            )
            if not accepted["row"]["accepted"]:
                raise RuntimeError(
                    "accepted trajectory was not reproducible on materialization"
                )
        selected_tree = output_dir / "selected.trees"
        shutil.copy2(Path(accepted["tree_path"]), selected_tree)

    accepted_row = accepted["row"]
    attempts = pd.DataFrame(attempt_rows).sort_values("attempt")
    attempts["used_for_rejection_decision"] = (
        attempts["attempt"] <= int(accepted_row["attempt"])
    )
    attempts.to_csv(output_dir / "selected_rejection_attempts.tsv", sep="\t", index=False)
    pd.DataFrame([accepted_row]).to_csv(
        output_dir / "selected_observation.tsv", sep="\t", index=False
    )
    trajectory = pd.DataFrame(accepted["trajectory"])
    trajectory["selected_attempt"] = int(accepted_row["attempt"])
    trajectory.to_csv(
        output_dir / "selected_allele_frequency_trajectory.tsv", sep="\t", index=False
    )
    ts = tskit.load(selected_tree)
    carrier_pairs = focal_carrier_pair_table(ts, accepted["carrier_counts"])
    carrier_pairs.to_csv(
        output_dir / "selected_focal_carrier_pairs.tsv",
        sep="\t",
        index=False,
    )
    carrier_profile = pair_tmrca_profile_by_focal_copy(
        ts, carrier_positions, accepted["carrier_counts"]
    )
    carrier_profile.to_csv(
        output_dir / "selected_hom_alt_vs_hom_ref_tmrca_profile.tsv",
        sep="\t",
        index=False,
    )
    truth_positions = np.arange(0, design["sequence_length"], stride, dtype=float)
    truth = within_individual_tmrca_grid(
        ts, truth_positions, design["variant_age_generations"]
    )
    truth.to_csv(
        output_dir / "selected_truth_recent_probability_profile.tsv.gz",
        sep="\t",
        index=False,
    )
    truth_scan = truth.rename(
        columns={"mean_p_tmrca_lt_threshold": "observed_fraction_recent"}
    )
    _plot_observed_profile(
        truth_scan,
        center=center,
        sequence_length=design["sequence_length"],
        zoom_half_width=zoom_half_width,
        threshold_years=(
            design["variant_age_generations"] * design["generation_time_years"]
        ),
        window_size=stride,
        output_path=output_dir / "selected_truth_recent_probability_spatial.png",
        pair_count=design["sample_diploids"],
        series_label="tree-sequence truth",
        statistic_label="True P",
    )
    _plot_demography_and_variant(
        output_dir / "demography_and_variant_timing.png",
        ancestral_population_size=design["ancestral_population_size"],
        present_population_size=design["present_population_size"],
        size_change_generations_ago=design["size_change_generations_ago"],
        variant_age_generations=design["variant_age_generations"],
        generation_time_years=design["generation_time_years"],
        selection_coefficient=selection_coefficient,
    )
    _plot_allele_frequency_trajectory(
        trajectory,
        population_allele_frequency=accepted_row["population_allele_frequency"],
        sample_allele_frequency=accepted_row["sample_allele_frequency"],
        variant_age_generations=design["variant_age_generations"],
        size_change_generations_ago=design["size_change_generations_ago"],
        generation_time_years=design["generation_time_years"],
        selection_coefficient=selection_coefficient,
        output_path=output_dir / "selected_allele_frequency_trajectory.png",
    )
    plot_carrier_profile_figure(
        carrier_profile,
        replicate=int(accepted_row["attempt"]),
        population_allele_frequency=accepted_row["population_allele_frequency"],
        sequence_length=design["sequence_length"],
        center=center,
        zoom_half_width=zoom_half_width,
        output_path=output_dir / "selected_hom_alt_vs_hom_ref_tmrca.png",
    )
    theoretical = two_epoch_recent_probability(
        ancestral_population_size=design["ancestral_population_size"],
        present_population_size=design["present_population_size"],
        size_change_generations_ago=design["size_change_generations_ago"],
        threshold_generations=design["variant_age_generations"],
    )
    exceedances, truth_p = _plot_null_calibration(
        neutral,
        accepted_row["center_fraction_recent"],
        theoretical_fraction_recent=theoretical,
        threshold_generations=design["variant_age_generations"],
        generation_time_years=design["generation_time_years"],
        selection_coefficient=selection_coefficient,
        output_path=output_dir / "truth_null_and_selected_pvalue.png",
    )
    selected_vcf = output_dir / "selected.vcf.gz"
    tree_sequence_to_vcf(selected_tree, selected_vcf)
    metrics = {
        "demography": {
            "ancestral_population_size": design["ancestral_population_size"],
            "present_population_size": design["present_population_size"],
            "size_change_generations_ago": design["size_change_generations_ago"],
        },
        **{key: value for key, value in design.items() if key not in {
            "ancestral_population_size", "present_population_size",
            "size_change_generations_ago",
        }},
        "selection_coefficient": float(selection_coefficient),
        "minimum_population_allele_frequency": float(minimum_population_af),
        "neutral_replicates": int(len(neutral)),
        "neutral_truth_source": os.path.relpath(null_truth_dir, output_dir),
        "neutral_statistics_source": os.path.relpath(
            neutral_statistics_path,
            output_dir,
        ),
        "neutral_theoretical_fraction_recent": float(theoretical),
        "selected_tree_file": selected_tree.name,
        "selected_vcf_file": selected_vcf.name,
        "selected_rejection": {
            "search_mode": (
                "explicit_deterministic_attempt"
                if selected_attempt is not None
                else "seed_order_rejection_search"
            ),
            "accepted_attempt_zero_based": int(accepted_row["attempt"]),
            "attempts_to_accept_in_seed_order": int(accepted_row["attempt"] + 1),
            "trajectories_computed_in_parallel_batches": int(len(attempts)),
            "population_allele_frequency": float(
                accepted_row["population_allele_frequency"]
            ),
            "sample_allele_frequency": float(accepted_row["sample_allele_frequency"]),
            "n_hom_ref_pairs": int(accepted_row["n_hom_ref_pairs"]),
            "n_heterozygous_pairs": int(accepted_row["n_heterozygous_pairs"]),
            "n_hom_alt_pairs": int(accepted_row["n_hom_alt_pairs"]),
            "truth_center_fraction_recent": float(
                accepted_row["center_fraction_recent"]
            ),
            "truth_center_neutral_exceedances": int(exceedances),
            "truth_center_monte_carlo_p_upper": float(truth_p),
        },
        "workers_requested": int(workers),
        "workers_used": int(1 if selected_attempt is not None else workers),
        "screen_attempts_per_batch": int(
            1
            if selected_attempt is not None
            else workers
            if workers == 1
            else workers * 10
        ),
        "shared_neutral_ancestry_seconds": float(shared_ancestry_seconds),
        "base_seed": int(seed),
        "stride_bp": int(stride),
        "elapsed_seconds": float(perf_counter() - started),
    }
    with (output_dir / "metrics.json").open("w", encoding="utf-8") as handle:
        json.dump(metrics, handle, indent=2)
    return metrics


def finalize_high_af_selected(
    source_dir: str | Path,
    neutral_decoded_profiles: str | Path,
    *,
    stride: int = DEFAULT_OUTPUT_STRIDE,
) -> dict:
    """Join the selected official decode to truth and reuse decoded nulls."""
    source_dir = Path(source_dir).resolve()
    observed_path = source_dir / "selected_decoded_recent_probability_profile.tsv"
    observed = pd.read_csv(observed_path, sep="\t")
    truth = pd.read_csv(
        source_dir / "selected_truth_recent_probability_profile.tsv.gz", sep="\t"
    )
    comparison = observed.merge(
        truth.rename(columns={
            "mean_p_tmrca_lt_threshold": "truth_fraction_recent",
            "mean_tmrca_generations": "mean_tmrca_generations_truth",
        }),
        on=["position_0based", "n_pairs"],
        how="inner",
        validate="one_to_one",
    ).rename(columns={
        "mean_tmrca_generations": "mean_tmrca_generations_decoded"
    })
    comparison.to_csv(
        source_dir / "selected_decoded_vs_truth.tsv.gz", sep="\t", index=False
    )
    result = finalize_container_stride_study(
        source_dir,
        source_dir,
        stride=stride,
        neutral_profiles_path=neutral_decoded_profiles,
    )
    with (source_dir / "metrics.json").open(encoding="utf-8") as handle:
        preparation = json.load(handle)
    result["truth_center_fraction_recent"] = preparation["selected_rejection"][
        "truth_center_fraction_recent"
    ]
    result["truth_center_monte_carlo_p_upper"] = preparation[
        "selected_rejection"
    ]["truth_center_monte_carlo_p_upper"]
    result["population_allele_frequency"] = preparation["selected_rejection"][
        "population_allele_frequency"
    ]
    result["sample_allele_frequency"] = preparation["selected_rejection"][
        "sample_allele_frequency"
    ]
    with (source_dir / "decoded_study_metrics.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(result, handle, indent=2)
    return result
