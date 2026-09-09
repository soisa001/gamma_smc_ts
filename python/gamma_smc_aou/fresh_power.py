"""Fresh, resumable decoded-power experiment; no old simulation inputs.

SLiM selected replicates condition on population survival, then reject a draw
if the selected allele is absent from the sampled panel. Fixation is retained.
Both engines simulate 11 Mb and are decoded after cropping to 10 Mb. The neutral
crop places the nearest segregating site at the focal coordinate without an AF
floor. Exact pair and position manifests are shared across both engines.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, FIRST_COMPLETED, wait
from contextlib import redirect_stdout, redirect_stderr
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import random
import shutil
import sys
import time
import traceback

os.environ["MPLBACKEND"] = "Agg"
for _name in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[_name] = "1"

import msprime
import numpy as np
import pandas as pd
import stdpopsim
import tskit

from .decoder import run_within_decoder
from .eas_sweep_models import build_eas_demography_models, load_phlash_eas_npz
from .run2_models import scoped_focal_overlay_patch
from .run7_models import scoped_slim_patch

REPO = Path(__file__).resolve().parents[2]
SITE_ID = "fresh_power_selected_site"
ARTIFACTS = (
    "simulation.trees",
    "decoded_input.trees",
    "frac_recent.tsv",
    "truth_frac_recent.tsv",
    "focal_carriers.npy",
)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def canonical_hash(value):
    return hashlib.sha256(
        json.dumps(
            value, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode()
    ).hexdigest()


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def seed_for(cfg, task, attempt=0, component="simulation"):
    # Content-addressed arm identity: adding/reordering coefficients never
    # changes a prior replicate's seed. The fixed study seed is never advanced.
    key = [
        cfg["seed_base"],
        cfg["population"],
        cfg["pulse_years"],
        cfg["selection_onset_years"],
        task["mode"],
        task["s"],
        task["replicate"],
        attempt,
        component,
    ]
    return int(canonical_hash(key)[:16], 16) % (2**31 - 2) + 1


def task_id(task):
    arm = (
        "neutral"
        if task["mode"] == "neutral"
        else f"selected_s{task['s']:.3f}".replace(".", "p")
    )
    return f"{arm}/rep{task['replicate']:04d}"


def plan(cfg):
    # Interleave neutral and selected work; no arm can monopolise the queue.
    tasks = []
    for i in range(
        max(cfg["neutral_replicates"], cfg["selected_replicates_per_coefficient"])
    ):
        if i < cfg["neutral_replicates"]:
            tasks.append(dict(mode="neutral", s=0.0, replicate=i))
        if i < cfg["selected_replicates_per_coefficient"]:
            tasks.extend(
                dict(mode="selected", s=float(s), replicate=i)
                for s in cfg["selection_coefficients"]
            )
    return tasks


def validate_config(cfg):
    if (
        cfg["population"] != "EAS"
        or cfg["statistic"] != "frac_recent_T"
        or cfg["recent_call"] != "mean"
    ):
        raise ValueError(
            "This version implements EAS and frac_recent_T from posterior means only"
        )
    if not 0 < cfg["selection_onset_years"] <= cfg["pulse_years"]:
        raise ValueError(
            "Selection must start at or after introgression, before present"
        )
    if cfg["scored_length_bp"] != 2 * cfg["focal_position_bp"]:
        raise ValueError("Scored region must be centered on the focal coordinate")
    if cfg["simulated_length_bp"] <= cfg["scored_length_bp"]:
        raise ValueError("Simulation requires flanking space for site alignment")
    possible = 2 * cfg["sample_diploids"] * (2 * cfg["sample_diploids"] - 1) // 2
    if not 0 < cfg["haplotype_pairs"] <= possible:
        raise ValueError("Requested pair count is not possible for this panel")
    if cfg["tmrca_cutoffs_years"] != sorted(set(cfg["tmrca_cutoffs_years"])):
        raise ValueError("Cutoffs must be increasing and unique")


def model(cfg):
    artifact = load_phlash_eas_npz(
        REPO / cfg["phlash_resource"], expected_sha256=cfg["phlash_sha256"]
    )
    history = build_eas_demography_models(artifact)["median"]
    dem = msprime.Demography()
    dem.add_population(name="EAS", initial_size=float(history.ne[0]))
    dem.add_population(name="Neanderthal", initial_size=cfg["archaic_effective_size"])
    for t, ne in zip(history.time_generations[1:], history.ne[1:]):
        dem.add_population_parameters_change(
            time=float(t), initial_size=float(ne), population="EAS"
        )
    pulse = cfg["pulse_years"] / cfg["generation_time_years"]
    dem.add_mass_migration(
        time=pulse,
        source="EAS",
        dest="Neanderthal",
        proportion=cfg["introgression_proportion"],
    )
    # Preserve the earlier study's archaic branch on BOTH sides. Decoding can
    # depend on older genealogy even for a TMRCA cutoff at/below the pulse.
    dem.add_population_parameters_change(
        time=pulse,
        initial_size=cfg["archaic_bottleneck_size"],
        population="Neanderthal",
    )
    dem.add_population_parameters_change(
        time=pulse + cfg["archaic_bottleneck_generations"],
        initial_size=cfg["archaic_effective_size"],
        population="Neanderthal",
    )
    dem.add_mass_migration(
        time=cfg["archaic_split_years"] / cfg["generation_time_years"],
        source="Neanderthal",
        dest="EAS",
        proportion=1,
    )
    dem.sort_events()
    dem.validate()
    return stdpopsim.DemographicModel(
        id="FreshEASQ02",
        description="PHLASH EAS with matched archaic pulse",
        long_description="Fresh survival-conditioned power simulation",
        model=dem,
        generation_time=cfg["generation_time_years"],
        mutation_rate=cfg["mutation_rate"],
        recombination_rate=cfg["recombination_rate"],
    )


def selected_events(cfg, s):
    pulse = cfg["pulse_years"] / cfg["generation_time_years"]
    onset = cfg["selection_onset_years"] / cfg["generation_time_years"]
    split = cfg["archaic_split_years"] / cfg["generation_time_years"]
    return [
        stdpopsim.DrawMutation(
            time=stdpopsim.GenerationAfter(split),
            single_site_id=SITE_ID,
            population="Neanderthal",
        ),
        stdpopsim.ChangeMutationFitness(
            start_time=onset,
            end_time=0,
            single_site_id=SITE_ID,
            population="EAS",
            selection_coeff=s,
            dominance_coeff=cfg["dominance"],
        ),
        stdpopsim.ConditionOnAlleleFrequency(
            start_time=pulse - cfg["slim_scaling_factor"],
            end_time=0,
            single_site_id=SITE_ID,
            population="EAS",
            op=">",
            allele_frequency=0,
        ),
    ]


def ordered_nodes(ts):
    samples = set(map(int, ts.samples()))
    groups = [
        [int(n) for n in ind.nodes if int(n) in samples] for ind in ts.individuals()
    ]
    groups = [g for g in groups if g]
    if any(len(g) != 2 for g in groups):
        raise ValueError("Non-diploid sampled individual")
    nodes = np.array([n for g in groups for n in g], dtype=int)
    if len(nodes) != len(samples) or set(nodes) != samples:
        raise ValueError(
            "Sample ordering does not cover the sampled panel exactly once"
        )
    return nodes


def focal_variant(ts, nodes, position=None):
    center = ts.sequence_length / 2
    best = None
    for variant in ts.variants(samples=nodes):
        pos = float(variant.site.position)
        if position is not None and pos != position:
            continue
        genotypes = np.asarray(variant.genotypes)
        if np.any(genotypes < 0) or genotypes.max() > 1:
            continue
        af = float(np.mean(genotypes != 0))
        if not (af > 0 and (position is not None or af < 1)):
            continue
        distance = abs(pos - center)
        if best is None or distance < best["distance_bp"]:
            best = dict(
                position=pos, af=af, distance_bp=distance, carriers=genotypes != 0
            )
        if position is not None or (
            best is not None and pos - center > best["distance_bp"]
        ):
            break
    return best


def crop_at_site(ts, position, cfg):
    offset = position - cfg["focal_position_bp"]
    end = offset + cfg["scored_length_bp"]
    if offset < 0 or end > ts.sequence_length:
        raise ValueError("Focal-site shift exceeds the simulated flanks")
    # SLiM attaches an ancestral reference string; tskit cannot trim it
    # implicitly. Crop that string explicitly while preserving allele states.
    tables = ts.dump_tables()
    reference = tables.reference_sequence.data
    tables.reference_sequence.clear()
    cropped = (
        tables.tree_sequence().keep_intervals([[offset, end]], simplify=True).trim()
    )
    if reference:
        tables = cropped.dump_tables()
        tables.reference_sequence.data = reference[int(offset) : int(end)]
        cropped = tables.tree_sequence()
    if cropped.sequence_length != cfg["scored_length_bp"]:
        raise ValueError("Cropped length mismatch")
    return cropped, offset


def read_profile(path, cfg):
    frame = pd.read_csv(path, sep="\t")
    grid = np.arange(0, cfg["scored_length_bp"], cfg["stride_bp"])
    if not np.array_equal(frame["position_0based"].to_numpy(), grid):
        raise ValueError(f"Profile has incorrect positions: {path}")
    columns = [f"frac_recent_{t}" for t in cfg["tmrca_cutoffs_years"]]
    values = frame[columns].to_numpy(dtype=float)
    if (
        not np.isfinite(values).all()
        or np.any((values < 0) | (values > 1))
        or np.any(np.diff(values, axis=1) < -1e-12)
    ):
        raise ValueError(f"Invalid probability or cutoff ordering in {path}")
    if not np.allclose(
        values * cfg["haplotype_pairs"],
        np.round(values * cfg["haplotype_pairs"]),
        atol=1e-6,
        rtol=0,
    ):
        raise ValueError("frac_recent is not a fraction of the declared pair count")
    return frame[["position_0based"] + columns]


def truth_profile(ts, pairs, cfg):
    nodes = ordered_nodes(ts)
    node_pairs = nodes[pairs]
    grid = np.arange(0, cfg["scored_length_bp"], cfg["stride_bp"])
    cutoffs = np.array(cfg["tmrca_cutoffs_years"]) / cfg["generation_time_years"]
    values = np.empty((len(grid), len(cutoffs)))
    row = 0
    # Compute once per marginal tree, without retaining a 1000 x 10000 matrix.
    for tree in ts.trees():
        if row == len(grid):
            break
        if grid[row] >= tree.interval.right:
            continue
        times = np.fromiter(
            (tree.tmrca(int(a), int(b)) for a, b in node_pairs), dtype=float
        )
        stats = (times[:, None] < cutoffs).mean(axis=0)
        while row < len(grid) and grid[row] < tree.interval.right:
            values[row] = stats
            row += 1
    if row != len(grid):
        raise ValueError("Truth profile did not cover the entire region")
    frame = pd.DataFrame(
        values, columns=[f"frac_recent_{t}" for t in cfg["tmrca_cutoffs_years"]]
    )
    frame.insert(0, "position_0based", grid)
    return frame


def completed_record(directory, fingerprint, cfg):
    path = directory / "complete.json"
    if not path.exists():
        return None
    record = json.loads(path.read_text())
    if record["fingerprint"] != fingerprint:
        raise ValueError(
            f"Different scientific inputs in existing replicate: {directory}"
        )
    for name in ARTIFACTS:
        if (
            not (directory / name).is_file()
            or digest(directory / name) != record["artifacts"][name]["sha256"]
        ):
            raise ValueError(
                f"Completed artifact is missing or corrupt: {directory / name}"
            )
    read_profile(directory / "frac_recent.tsv", cfg)
    return record


def run_one(payload):
    cfg, task, out, fingerprint, slim, decoder = payload
    out = Path(out)
    directory = out / task_id(task)
    directory.mkdir(parents=True, exist_ok=True)
    if (directory / "failed.json").exists():
        history = directory / "failure_history"
        history.mkdir(exist_ok=True)
        os.replace(
            directory / "failed.json", history / f"failure_{time.time_ns()}.json"
        )
    started = time.perf_counter()
    try:
        with (
            (directory / "execution.log").open("a", buffering=1) as log,
            redirect_stdout(log),
            redirect_stderr(log),
        ):
            receipt = directory / "simulated.json"
            if receipt.exists():
                sim = json.loads(receipt.read_text())
                if sim["fingerprint"] != fingerprint:
                    raise ValueError(
                        "Simulation inputs differ from the existing phase receipt"
                    )
                for name, expected in sim["simulation_hashes"].items():
                    if digest(directory / name) != expected:
                        raise ValueError(f"Simulation file corrupted: {name}")
                cropped = tskit.load(directory / "decoded_input.trees")
            else:
                built = model(cfg)
                attempts = []
                for attempt in range(cfg["maximum_sample_attempts"]):
                    raw_receipt = directory / "raw_simulated.json"
                    if raw_receipt.exists():
                        raw = json.loads(raw_receipt.read_text())
                        if (
                            raw["fingerprint"] != fingerprint
                            or digest(directory / "simulation.trees") != raw["sha256"]
                        ):
                            raise ValueError(
                                "Raw simulation receipt or hash does not match"
                            )
                        ts = tskit.load(directory / "simulation.trees")
                        attempts = raw["attempts"]
                        seed = raw["seed"]
                        position = (
                            None
                            if task["mode"] == "neutral"
                            else cfg["simulated_length_bp"] // 2
                        )
                        site = focal_variant(ts, ordered_nodes(ts), position)
                        if site is None:
                            raise ValueError(
                                "Saved accepted simulation lacks its focal allele"
                            )
                        break
                    seed = seed_for(cfg, task, attempt)
                    random.seed(seed)
                    np.random.seed(seed)
                    tick = time.perf_counter()
                    if task["mode"] == "neutral":
                        ts = msprime.sim_ancestry(
                            samples={"EAS": cfg["sample_diploids"]},
                            demography=built.model,
                            sequence_length=cfg["simulated_length_bp"],
                            recombination_rate=cfg["recombination_rate"],
                            ploidy=2,
                            random_seed=seed,
                        )
                        mutation_seed = seed_for(cfg, task, attempt, "mutations")
                        ts = msprime.sim_mutations(
                            ts, rate=cfg["mutation_rate"], random_seed=mutation_seed
                        )
                        site = focal_variant(ts, ordered_nodes(ts))
                    else:
                        center = cfg["simulated_length_bp"] // 2
                        contig = stdpopsim.get_species("HomSap").get_contig(
                            length=cfg["simulated_length_bp"],
                            mutation_rate=cfg["mutation_rate"],
                            recombination_rate=cfg["recombination_rate"],
                        )
                        contig.add_single_site(id=SITE_ID, coordinate=center)
                        with (
                            scoped_slim_patch("selected", 0, patch_placement=True),
                            scoped_focal_overlay_patch(center) as overlay,
                        ):
                            ts = stdpopsim.get_engine("slim").simulate(
                                built,
                                contig,
                                {"EAS": cfg["sample_diploids"]},
                                seed=seed,
                                extended_events=selected_events(cfg, task["s"]),
                                slim_path=slim,
                                slim_scaling_factor=cfg["slim_scaling_factor"],
                                slim_burn_in=cfg["slim_burn_in"],
                                keep_mutation_ids_as_alleles=False,
                            )
                        if overlay["masked_calls"] < 1:
                            raise ValueError(
                                "The reserved selected-site mutation guard did not run"
                            )
                        site = focal_variant(ts, ordered_nodes(ts), center)
                    attempts.append(
                        dict(
                            attempt=attempt,
                            seed=seed,
                            accepted=site is not None,
                            seconds=time.perf_counter() - tick,
                        )
                    )
                    atomic_json(directory / "attempts.json", attempts)
                    if site is not None:
                        break
                else:
                    raise RuntimeError(
                        "No sampled surviving selected allele within the attempt limit"
                    )
                raw_temp = directory / "simulation.trees.tmp"
                ts.dump(raw_temp)
                os.replace(raw_temp, directory / "simulation.trees")
                atomic_json(
                    directory / "raw_simulated.json",
                    dict(
                        fingerprint=fingerprint,
                        seed=seed,
                        attempts=attempts,
                        sha256=digest(directory / "simulation.trees"),
                    ),
                )
                cropped, offset = crop_at_site(ts, site["position"], cfg)
                aligned = focal_variant(
                    cropped, ordered_nodes(cropped), cfg["focal_position_bp"]
                )
                if aligned is None or aligned["af"] != site["af"]:
                    raise ValueError(
                        "Cropping changed the focal allele or its sample frequency"
                    )
                for name, tree in (("decoded_input.trees", cropped),):
                    temporary = directory / (name + ".tmp")
                    tree.dump(temporary)
                    if tskit.load(temporary).num_samples != 2 * cfg["sample_diploids"]:
                        raise ValueError("Saved tree has the wrong sample size")
                    os.replace(temporary, directory / name)
                np.save(
                    directory / "focal_carriers.npy",
                    aligned["carriers"],
                    allow_pickle=False,
                )
                sim = dict(
                    task=task,
                    fingerprint=fingerprint,
                    seed=seed,
                    attempts=attempts,
                    sample_af=site["af"],
                    fixed=site["af"] == 1,
                    original_focal_position=site["position"],
                    crop_offset=offset,
                    simulation_seconds=time.perf_counter() - started,
                    simulation_hashes={
                        n: digest(directory / n)
                        for n in (
                            "simulation.trees",
                            "decoded_input.trees",
                            "focal_carriers.npy",
                        )
                    },
                )
                atomic_json(receipt, sim)
            result = run_within_decoder(
                decoder,
                directory / "decoded_input.trees",
                directory / "gamma_smc_summary.tsv",
                threshold_years=cfg["tmrca_cutoffs_years"],
                generation_time=cfg["generation_time_years"],
                mutation_rate=cfg["mutation_rate"],
                scaled_mutation_rate=cfg["decoder_scaled_mutation_rate"],
                recombination_to_mutation_ratio=cfg[
                    "decoder_recombination_to_mutation_ratio"
                ],
                output_at_stride=-1,
                output_at_hets=False,
                output_positions_file=out / "positions.txt",
                vcf_position_transform="one_based",
                only_within=False,
                pairs_file=out / "pairs.tsv",
                recent_call="mean",
                threads=1,
                cache_size=cfg["decoder_cache_size"],
                pair_block=cfg["decoder_pair_block"],
                exp10=cfg["decoder_exp10"],
                backward_alignment=cfg["decoder_backward_alignment"],
            )
            for channel in ("stdout", "stderr"):
                (directory / f"decoder.{channel}.log").write_text(
                    result.pop(channel), encoding="utf-8"
                )
            atomic_json(directory / "gamma_smc_summary.tsv.run.json", result)
            frame = read_profile(directory / "gamma_smc_summary.tsv", cfg)
            frame.to_csv(directory / "frac_recent.tsv", sep="\t", index=False)
            pairs = np.loadtxt(out / "pairs.tsv", dtype=int)
            truth_profile(cropped, pairs, cfg).to_csv(
                directory / "truth_frac_recent.tsv", sep="\t", index=False
            )
            read_profile(directory / "truth_frac_recent.tsv", cfg)
            record = dict(
                sim,
                decode=result,
                elapsed_seconds=time.perf_counter() - started,
                completed_at=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                artifacts={
                    n: dict(
                        sha256=digest(directory / n),
                        bytes=(directory / n).stat().st_size,
                    )
                    for n in ARTIFACTS
                },
            )
            atomic_json(directory / "complete.json", record)
        return dict(
            task_id=task_id(task),
            status="completed",
            seconds=time.perf_counter() - started,
            af=record["sample_af"],
        )
    except Exception as error:
        failed = dict(
            task_id=task_id(task),
            status="failed",
            error=str(error),
            traceback=traceback.format_exc(),
        )
        atomic_json(directory / "failed.json", failed)
        return dict(
            task_id=task_id(task),
            status="failed",
            error=str(error).splitlines()[0],
            details=str(directory / "failed.json"),
        )


def empirical_p(selected, neutral):
    neutral = np.asarray(neutral, dtype=float)
    selected = np.asarray(selected, dtype=float)
    if (
        neutral.size == 0
        or not np.isfinite(neutral).all()
        or not np.isfinite(selected).all()
    ):
        raise ValueError("Finite nonempty reference required")
    return (
        1 + neutral.size - np.searchsorted(np.sort(neutral), selected, side="left")
    ) / (neutral.size + 1)


def analyse(out, cfg, fingerprint):
    rows = []
    for task in plan(cfg):
        directory = out / task_id(task)
        record = completed_record(directory, fingerprint, cfg)
        if record is None:
            raise ValueError(f"Cannot publish final power: incomplete {task_id(task)}")
        frame = read_profile(directory / "frac_recent.tsv", cfg)
        focal = frame[frame.position_0based == cfg["focal_position_bp"]].iloc[0]
        rows.append(
            dict(
                task_id=task_id(task),
                **task,
                sample_af=record["sample_af"],
                fixed=record["fixed"],
                **{
                    f"frac_recent_{t}": float(focal[f"frac_recent_{t}"])
                    for t in cfg["tmrca_cutoffs_years"]
                },
            )
        )
    data = pd.DataFrame(rows)
    null = data[data["mode"] == "neutral"]
    selected = data[data["mode"] == "selected"].copy()
    summary = []
    for cutoff in cfg["tmrca_cutoffs_years"]:
        column = f"frac_recent_{cutoff}"
        selected[f"p_{cutoff}"] = empirical_p(selected[column], null[column])
        for s, group in selected.groupby("s"):
            hits = int((group[f"p_{cutoff}"] < cfg["alpha"]).sum())
            n = len(group)
            power = hits / n
            # Wilson binomial interval: conditional on the estimated null.
            z = 1.959963984540054
            denominator = 1 + z * z / n
            middle = (power + z * z / (2 * n)) / denominator
            radius = (
                z * np.sqrt(power * (1 - power) / n + z * z / (4 * n * n)) / denominator
            )
            summary.append(
                dict(
                    s=s,
                    cutoff_years=cutoff,
                    n=n,
                    detected=hits,
                    power=power,
                    ci_low=max(0, middle - radius),
                    ci_high=min(1, middle + radius),
                    alpha=cfg["alpha"],
                    sample_af_mean=group.sample_af.mean(),
                    n_fixed=int(group.fixed.sum()),
                )
            )
    analysis = out / "analysis"
    analysis.mkdir(exist_ok=True)
    data.to_csv(analysis / "focal_statistics.csv", index=False)
    selected.to_csv(analysis / "selected_pvalues.csv", index=False)
    power = pd.DataFrame(summary)
    power.to_csv(analysis / "power.csv", index=False)
    import matplotlib

    matplotlib.use("Agg")
    matplotlib.rcParams.update(
        {
            "font.size": 18,
            "axes.labelsize": 20,
            "axes.titlesize": 20,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(11, 8.5))
    for cutoff, group in power.groupby("cutoff_years"):
        group = group.sort_values("s")
        ax.plot(group.s, group.power, "o-", label=f"{cutoff / 1000:g} kya")
    ax.set(
        xlabel="Selection coefficient",
        ylabel="Power at the focal stride",
        ylim=(-0.02, 1.02),
    )
    ax.legend(
        ncol=3,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.02),
        frameon=False,
        fontsize=15,
    )
    ax.spines[["top", "right"]].set_visible(False)
    fig.subplots_adjust(top=0.77, bottom=0.13, left=0.13, right=0.97)
    fig.savefig(analysis / "power.png", dpi=200)
    fig.savefig(analysis / "power.pdf")
    plt.close(fig)
    atomic_json(
        analysis / "provenance.json",
        dict(
            fingerprint=fingerprint,
            config=cfg,
            neutral_replicates=len(null),
            selected_replicates=len(selected),
            interval="95% Wilson intervals conditional on the estimated null; null Monte Carlo uncertainty is not included",
            scope="Prespecified focal stride in isolated 10 Mb regions; no genome-wide error-control claim",
        ),
    )


def initialise(out, cfg, slim, decoder):
    validate_config(cfg)
    out.mkdir(parents=True, exist_ok=True)
    np.random.seed(cfg["seed_base"])
    random.seed(cfg["seed_base"])
    versions = {
        p: importlib.metadata.version(p)
        for p in ("msprime", "tskit", "stdpopsim", "numpy", "pyslim")
    }
    if versions["stdpopsim"] != "0.3.0":
        raise ValueError("Selected-event patches require stdpopsim 0.3.0")
    identity = {
        k: v
        for k, v in cfg.items()
        if k
        not in (
            "workers",
            "memory_limit_gb",
            "storage_limit_bytes",
            "selected_replicates_per_coefficient",
            "neutral_replicates",
            "alpha",
            "selection_coefficients",
        )
    }
    identity.update(
        versions=versions, slim_sha256=digest(slim), decoder_sha256=digest(decoder)
    )
    fingerprint = canonical_hash(identity)
    if (out / "manifest.json").exists():
        previous = json.loads((out / "manifest.json").read_text())
        if previous["fingerprint"] != fingerprint:
            raise ValueError(
                "Existing study has different scientific parameters; use a new output directory"
            )
    possible = np.column_stack(np.triu_indices(2 * cfg["sample_diploids"], k=1))
    indices = np.random.default_rng(cfg["pairs_seed"]).choice(
        len(possible), cfg["haplotype_pairs"], replace=False
    )
    pairs = possible[np.sort(indices)]
    for name, values in (
        ("pairs.tsv", pairs),
        ("positions.txt", np.arange(0, cfg["scored_length_bp"], cfg["stride_bp"])),
    ):
        path = out / name
        if path.exists():
            if not np.array_equal(np.loadtxt(path, dtype=int), values):
                raise ValueError(f"Existing sample/position manifest differs: {path}")
        else:
            np.savetxt(path, values, fmt="%d", delimiter="\t")
    atomic_json(
        out / "manifest.json",
        dict(
            config=cfg,
            fingerprint=fingerprint,
            versions=versions,
            python=sys.version,
            platform=platform.platform(),
            output_root=str(out),
            repo_root=str(REPO),
            slim=str(slim),
            decoder=str(decoder),
            scientific_identity=identity,
            pairs_sha256=digest(out / "pairs.tsv"),
            positions_sha256=digest(out / "positions.txt"),
            source_hashes={
                str(p.relative_to(REPO)): digest(p)
                for p in (
                    Path(__file__),
                    REPO / "python/gamma_smc_aou/decoder.py",
                    REPO / "python/gamma_smc_aou/run7_models.py",
                    REPO / "python/gamma_smc_aou/run2_models.py",
                )
            },
            notes=[
                "No reused simulation outputs",
                "Both engines simulate 11 Mb and decode a site-aligned 10 Mb crop",
                "2% demographic pulse, without conditioning on realised ancestry or a pulse AF band",
                "Selected sample AF > 0; AF=1 retained; neutral nearest segregating site has no MAF floor",
                "Archaic bottleneck Ne=10 for 100 generations retained from the prior design and matched in the null",
            ],
        ),
    )
    atomic_json(
        out / "sample_manifest.json",
        [dict(task_id=task_id(t), **t, seed=seed_for(cfg, t)) for t in plan(cfg)],
    )
    return fingerprint


def directory_bytes(path):
    return sum(p.stat().st_size for p in Path(path).rglob("*") if p.is_file())


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path, default=REPO / "configs/eas_q02_50k.json"
    )
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--slim", type=Path, required=True)
    parser.add_argument("--decoder", type=Path, default=REPO / "bin/gamma_smc")
    parser.add_argument("--workers", type=int)
    parser.add_argument(
        "--phase",
        choices=("smoke", "simulate", "run", "analyse", "status"),
        default="run",
    )
    args = parser.parse_args(argv)
    cfg = json.loads(args.config.read_text(encoding="utf-8-sig"))
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    workers = args.workers or cfg["workers"]
    if not 1 <= workers <= 24:
        raise ValueError("Local worker count must be between 1 and 24")
    import fcntl

    with (out / "runner.lock").open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            print(
                "A runner already holds this study lock; no duplicate work was started."
            )
            if (out / "run_status.json").exists():
                print((out / "run_status.json").read_text())
            return 0
        fingerprint = initialise(out, cfg, args.slim.resolve(), args.decoder.resolve())
        tasks = plan(cfg)
        if args.phase == "smoke":
            tasks = [
                t
                for t in tasks
                if (t["mode"] == "neutral" and t["replicate"] < 3)
                or (
                    t["mode"] == "selected"
                    and t["replicate"] == 0
                    and t["s"] in (0.001, 0.005, 0.010)
                )
            ]
        done = []
        todo = []
        for task in tasks:
            if completed_record(out / task_id(task), fingerprint, cfg):
                done.append(task)
            else:
                todo.append(task)

        def status(state, failed=None):
            compact_failures = [
                dict(task_id=f.get("task_id"), error=f["error"].splitlines()[0])
                for f in (failed or [])
            ]
            value = dict(
                state=state,
                phase=args.phase,
                completed=len(done),
                target=len(tasks),
                workers=workers,
                failed=compact_failures,
                updated_at=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            )
            atomic_json(out / "run_status.json", value)
            print(json.dumps(value), flush=True)

        status("validated")
        if args.phase == "status":
            return 0
        if args.phase == "analyse":
            analyse(out, cfg, fingerprint)
            status("complete")
            return 0
        # All temporary SLiM checkpoints and logs are on the data disk too.
        temporary = out / "tmp"
        temporary.mkdir(exist_ok=True)
        os.environ["TMPDIR"] = str(temporary)
        import tempfile

        tempfile.tempdir = str(temporary)
        failures = []
        pending = {}
        cursor = 0
        with ProcessPoolExecutor(max_workers=workers) as pool:
            while cursor < len(todo) or pending:
                while cursor < len(todo) and len(pending) < workers and not failures:
                    if cursor % 20 == 0:
                        # Whole-volume usage is a conservative upper bound on
                        # simulation storage. Avoid repeatedly stat-ing every
                        # old artifact over the Windows/WSL filesystem bridge.
                        volume = shutil.disk_usage(out)
                        reserve = 20_000_000_000
                        if (
                            volume.used + reserve > cfg["storage_limit_bytes"]
                            or volume.free < reserve
                        ):
                            failures.append(
                                dict(
                                    error="Storage reserve reached; stopping new submissions"
                                )
                            )
                            break
                    task = todo[cursor]
                    cursor += 1
                    future = pool.submit(
                        run_one,
                        (
                            cfg,
                            task,
                            str(out),
                            fingerprint,
                            str(args.slim.resolve()),
                            str(args.decoder.resolve()),
                        ),
                    )
                    pending[future] = task
                if not pending:
                    break
                finished, _ = wait(pending, timeout=30, return_when=FIRST_COMPLETED)
                if not finished:
                    status("running", failures)
                    continue
                for future in finished:
                    task = pending.pop(future)
                    try:
                        result = future.result()
                    except Exception as error:
                        result = dict(
                            status="failed", error=str(error), task_id=task_id(task)
                        )
                    print(json.dumps(result), flush=True)
                    if result["status"] == "completed":
                        done.append(task)
                    else:
                        failures.append(result)
                status("running", failures)
                if failures and not pending:
                    break
        if failures:
            status("failed", failures)
            return 1
        if len(done) != len(tasks):
            status("incomplete")
            return 1
        if args.phase == "run":
            analyse(out, cfg, fingerprint)
        if args.phase == "smoke":
            atomic_json(
                out / "smoke_validation.json",
                dict(
                    status="passed",
                    tasks=[task_id(t) for t in tasks],
                    fingerprint=fingerprint,
                ),
            )
        status("complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
