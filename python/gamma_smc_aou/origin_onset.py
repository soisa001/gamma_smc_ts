"""Matched focal-allele origin/onset pilot, with independent null test targets.

Simulation, audit, decode, and analysis are explicit, restartable phases.
All three cohort roles use the same focal-allele simulation engine/events.
"""
from __future__ import annotations

import os
os.environ["MPLBACKEND"] = "Agg"
for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_key] = "1"

import argparse
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from contextlib import contextmanager, redirect_stdout, redirect_stderr
import importlib.metadata
import io
import json
from pathlib import Path
import random
import shutil
import sys
import time
import traceback

import numpy as np
import pandas as pd
import stdpopsim
from stdpopsim import slim_engine
import tskit

from . import fresh_power as legacy
from .decoder import run_within_decoder
from .pair_class_profiles import pair_classes, posterior_arrays, summarize_decoded, summarize_truth, build_frame
from .run2_models import scoped_focal_overlay_patch
from .run7_models import scoped_slim_patch

REPO = legacy.REPO
SIM_CONTRACT = "origin-onset-focal/v1"
FAMILIES = {"I50": ("introgressed", 50000), "D50": ("de_novo_EAS", 50000), "D10": ("de_novo_EAS", 10000)}
SIM_KEYS = (
    "population", "introgression_proportion", "pulse_years", "sample_diploids",
    "mutation_rate", "recombination_rate", "generation_time_years", "simulated_length_bp",
    "scored_length_bp", "focal_position_bp", "dominance", "archaic_split_years",
    "archaic_effective_size", "archaic_bottleneck_size", "archaic_bottleneck_generations",
    "slim_burn_in", "seed_base", "phlash_resource", "phlash_sha256",
)


def tasks(cfg):
    result = []
    for rep in range(max(cfg["null_replicates"], cfg["target_replicates"])):
        for family, (origin, introduction) in FAMILIES.items():
            common = dict(family=family, origin=origin, introduction_years=introduction,
                          replicate=rep, scaling_factor=cfg["scaling_by_origin"][origin])
            if rep < cfg["null_replicates"]:
                result.append(dict(common, id=f"{family}/null/rep{rep:04d}", role="null", s=0.0, onset_years=0))
            if rep < cfg["target_replicates"]:
                result.append(dict(common, id=f"{family}/neutral_target/rep{rep:04d}", role="neutral_target", s=0.0, onset_years=0))
                for onset in ([50000, 10000] if origin == "introgressed" else [introduction]):
                    for s in cfg["selection_coefficients"]:
                        arm = f"onset{onset}_s{s:.3f}".replace(".", "p")
                        result.append(dict(common, id=f"{family}/{arm}/rep{rep:04d}", role="selected", s=s, onset_years=onset))
    return result


def seed_for(cfg, task, attempt=0):
    identity = [SIM_CONTRACT, cfg["seed_base"], task["family"], task["role"],
                task["onset_years"], task["s"], task["replicate"], attempt]
    return int(legacy.canonical_hash(identity)[:16], 16) % (2**31 - 2) + 1


def events(cfg, task):
    gen = cfg["generation_time_years"]
    if task["origin"] == "introgressed":
        draw = stdpopsim.GenerationAfter(cfg["archaic_split_years"] / gen)
        population = "Neanderthal"
        check = cfg["pulse_years"] / gen - task["scaling_factor"]
    else:
        draw = task["introduction_years"] / gen
        population = "EAS"
        check = stdpopsim.GenerationAfter(draw)
    result = [stdpopsim.DrawMutation(time=draw, single_site_id=legacy.SITE_ID, population=population)]
    if task["s"] > 0:
        result.append(stdpopsim.ChangeMutationFitness(
            start_time=task["onset_years"] / gen, end_time=0,
            single_site_id=legacy.SITE_ID, population="EAS",
            selection_coeff=task["s"], dominance_coeff=cfg["dominance"]))
    result.append(stdpopsim.ConditionOnAlleleFrequency(
        start_time=check, end_time=0, single_site_id=legacy.SITE_ID,
        population="EAS", op=">", allele_frequency=0))
    return result


@contextmanager
def simulation_patch(cfg, task, directory=None):
    """Record placement and the accepted EAS trajectory, including de novo census.

The draw event resets trajectory vectors on every internal survival restart.
This is a conditional trajectory, never an estimate of establishment probability.
"""
    original_makescript = slim_engine.slim_makescript
    with scoped_slim_patch("selected", 0, patch_placement=task["origin"] == "introgressed") as receipt:
        functions = slim_engine._slim_functions
        insertion = """
    metadata.setValue("run7_focal_mutation_type_id", mut_type.id);
    metadata.setValue("origin_placement_tick", community.tick);
    metadata.setValue("origin_placement_population", pop.id);
    metadata.setValue("origin_placement_total_copies", size(pop.genomes));
    metadata.setValue("origin_placement_alt_copies", sum(pop.genomes.containsMutations(sim.mutationsOfType(mut_type)[0])));
    metadata.setValue("origin_ticks", integer(0));
    metadata.setValue("origin_af", float(0));
"""
        if task["origin"] == "introgressed":
            anchor = '    metadata.setValue("run7_placement_mode", "fixed_in_archaic_population");'
        else:
            anchor = '   targets.addNewDrawnMutation(mut_type, pos);'
        if functions.count(anchor) != 1:
            raise ValueError("SLiM placement recording contract changed")
        functions = functions.replace(anchor, anchor + insertion)
        recorder = f"""
function (void)origin_record(void) {{
    if (community.tick < time_to_tick({task['introduction_years']})) return;
    mt_id = metadata.getValue("run7_focal_mutation_type_id");
    if (isNULL(mt_id)) return;
    pop = sim.subpopulations[sim.subpopulations.id == 0];
    if (size(pop) != 1) return;
    mt = sim.mutationTypes[sim.mutationTypes.id == mt_id];
    if (size(mt) != 1) return;
    muts = sim.mutationsOfType(mt);
    value = 0.0;
    if (size(muts) == 1) value = sim.mutationFrequencies(pop, muts[0]);
    else if (size(sim.substitutions[sim.substitutions.mutationType == mt]) == 1) value = 1.0;
    ticks = metadata.getValue("origin_ticks");
    values = metadata.getValue("origin_af");
    if (isNULL(ticks)) {{ ticks = integer(0); values = float(0); }}
    keep = ticks < community.tick;
    metadata.setValue("origin_ticks", c(ticks[keep], community.tick));
    metadata.setValue("origin_af", c(values[keep], value));
}}
"""
        functions += recorder
        anchor = '    metadata.setValue("run7_final_census_state", state);'
        functions = functions.replace(anchor, anchor + '\n    origin_record();')
        slim_engine._slim_functions = functions

        def makescript(stream, *args, **kwargs):
            buffer = io.StringIO()
            result = original_makescript(buffer, *args, **kwargs)
            script = buffer.getvalue() + '\n1: late() { origin_record(); }\n'
            stream.write(script)
            if directory is not None:
                (directory / "model.slim").write_text(script)
            return result
        slim_engine.slim_makescript = makescript
        try:
            yield receipt
        finally:
            slim_engine.slim_makescript = original_makescript


def scalar(value):
    while isinstance(value, (list, tuple)) and len(value) == 1:
        value = value[0]
    return value


def simulation_identity(cfg, task, runtime):
    return dict(contract=SIM_CONTRACT, parameters={k: cfg[k] for k in SIM_KEYS},
                task=task, runtime=runtime)


def artifact_specs(directory, names):
    return {name: dict(bytes=(directory / name).stat().st_size, sha256=legacy.digest(directory / name)) for name in names}


def read_receipt(directory, name, identity):
    path = directory / name
    if not path.exists():
        return None
    record = json.loads(path.read_text())
    if record["identity"] != identity:
        raise ValueError(f"Incompatible cache identity: {path}")
    for file, spec in record["artifacts"].items():
        source = directory / file
        if source.stat().st_size != spec["bytes"] or legacy.digest(source) != spec["sha256"]:
            raise ValueError(f"Corrupt or missing artifact: {source}")
    return record


def validate_trees(directory, cfg, record):
    carriers = np.load(directory / "focal_carriers.npy", allow_pickle=False)
    assert carriers.dtype == np.bool_ and carriers.shape == (2 * cfg["sample_diploids"],)
    assert 0 < carriers.mean() <= 1 and float(carriers.mean()) == record["sample_af"]
    for name, length, position in (
        ("simulation.trees", cfg["simulated_length_bp"], cfg["simulated_length_bp"] / 2),
        ("decoded_input.trees", cfg["scored_length_bp"], cfg["focal_position_bp"]),
    ):
        ts = tskit.load(directory / name)
        assert ts.num_samples == len(carriers) and ts.sequence_length == length
        focal = legacy.exact_focal_variant(ts, legacy.ordered_nodes(ts), position)
        assert focal is not None
        np.testing.assert_array_equal(focal["carriers"], carriers)
    return True


def maybe_reuse(cfg, task, directory, identity, legacy_root):
    if not legacy_root or not (task["role"] == "selected" and task["family"] == "I50"
                              and task["onset_years"] == 50000 and task["s"] <= .010):
        return None
    root = Path(legacy_root)
    old_task = dict(mode="selected", s=task["s"], replicate=task["replicate"])
    source = root / legacy.task_id(old_task)
    if not (source / "simulated.json").exists():
        return None
    manifest = json.loads((root / "manifest.json").read_text())
    old_cfg = manifest["config"]
    for key in SIM_KEYS:
        if old_cfg[key] != cfg[key]:
            raise ValueError(f"Legacy reuse mismatch: {key}")
    assert old_cfg["slim_scaling_factor"] == task["scaling_factor"]
    assert old_cfg["selection_onset_years"] == task["onset_years"]
    assert manifest["versions"] == identity["runtime"]["versions"]
    assert manifest["scientific_identity"]["slim_sha256"] == identity["runtime"]["slim_sha256"]
    old = legacy.simulated_record(source, manifest["fingerprint"], old_cfg, old_task)
    for name in legacy.SIMULATION_ARTIFACTS:
        # Copies preserve the archived originals when new derived files are written.
        shutil.copy2(source / name, directory / name)
    record = dict(identity=identity, task=task, seed=old["seed"], sample_af=old["sample_af"],
                  reused_from=str(source), original_receipt_sha256=legacy.digest(source / "simulated.json"),
                  trajectory_available=False, artifacts=artifact_specs(directory, legacy.SIMULATION_ARTIFACTS))
    validate_trees(directory, cfg, record)
    legacy.atomic_json(directory / "simulation.json", record)
    return record


def simulate(cfg, task, directory, identity, slim, legacy_root):
    found = read_receipt(directory, "simulation.json", identity)
    if found:
        validate_trees(directory, cfg, found)
        return "cached"
    if maybe_reuse(cfg, task, directory, identity, legacy_root):
        return "reused"
    model = legacy.model(cfg)
    center = cfg["simulated_length_bp"] // 2
    contig = stdpopsim.get_species("HomSap").get_contig(
        length=cfg["simulated_length_bp"], mutation_rate=cfg["mutation_rate"], recombination_rate=cfg["recombination_rate"])
    contig.add_single_site(id=legacy.SITE_ID, coordinate=center)
    attempts = []
    for attempt in range(cfg["maximum_sample_attempts"]):
        seed = seed_for(cfg, task, attempt)
        random.seed(seed)
        np.random.seed(seed)
        began = time.monotonic()
        with simulation_patch(cfg, task, directory), scoped_focal_overlay_patch(center) as overlay:
            ts = stdpopsim.get_engine("slim").simulate(
                model, contig, {"EAS": cfg["sample_diploids"]}, seed=seed,
                extended_events=events(cfg, task), slim_path=slim,
                slim_scaling_factor=task["scaling_factor"], slim_burn_in=cfg["slim_burn_in"],
                keep_mutation_ids_as_alleles=False)
        if overlay["masked_calls"] < 1:
            raise ValueError("Focal mutation overlay guard did not run")
        focal = legacy.exact_focal_variant(ts, legacy.ordered_nodes(ts), center)
        attempts.append(dict(attempt=attempt, seed=seed, accepted=focal is not None, seconds=time.monotonic() - began))
        legacy.atomic_json(directory / "attempts.json", attempts)
        if focal is not None:
            break
    else:
        raise RuntimeError("No observed focal allele within sample-attempt limit")
    metadata = ts.metadata["SLiM"]["user_metadata"]
    copies = int(scalar(metadata["origin_placement_alt_copies"]))
    total = int(scalar(metadata["origin_placement_total_copies"]))
    assert copies == (total if task["origin"] == "introgressed" else 1)
    ticks = np.array(metadata["origin_ticks"], dtype=float).ravel()
    frequencies = np.array(metadata["origin_af"], dtype=float).ravel()
    assert ticks.shape == frequencies.shape and len(ticks) > 0 and np.all(np.diff(ticks) > 0)
    assert np.isfinite(frequencies).all() and ((frequencies >= 0) & (frequencies <= 1)).all()
    years = (ticks[-1] - ticks) * task["scaling_factor"] * cfg["generation_time_years"]
    pd.DataFrame(dict(years_ago=years, population_af=frequencies)).to_csv(directory / "trajectory.csv", index=False)
    cropped, offset = legacy.crop_at_site(ts, center, cfg)
    for name, tree in (("simulation.trees", ts), ("decoded_input.trees", cropped)):
        temporary = directory / (name + ".tmp")
        tree.dump(temporary)
        temporary.replace(directory / name)
    np.save(directory / "focal_carriers.npy", focal["carriers"])
    record = dict(identity=identity, task=task, seed=seed, sample_af=focal["af"], attempts=attempts,
                  crop_offset=offset, trajectory_available=True, placement_alt_copies=copies,
                  placement_total_copies=total, final_census_af=float(scalar(metadata["run7_final_census_af"])),
                  conditional_on="population survival and focal observation in the sample",
                  artifacts=artifact_specs(directory, (*legacy.SIMULATION_ARTIFACTS, "trajectory.csv", "model.slim")))
    validate_trees(directory, cfg, record)
    legacy.atomic_json(directory / "simulation.json", record)
    return "simulated"


def decode(cfg, task, directory, identity, out, decoder):
    sim = read_receipt(directory, "simulation.json", identity)
    if sim is None:
        raise ValueError(f"Missing simulation: {directory}")
    decode_identity = dict(contract="origin-focal-decode/v1",
        tree=sim["artifacts"]["decoded_input.trees"], carriers=sim["artifacts"]["focal_carriers.npy"],
        pairs=legacy.digest(out / "pairs.tsv"), positions=legacy.digest(out / "positions.txt"),
        decoder=legacy.digest(decoder), parameters={k: v for k, v in cfg.items() if k.startswith("decoder_") or k in (
            "tmrca_cutoffs_years", "generation_time_years", "mutation_rate", "recent_call")})
    if read_receipt(directory, "decoded.json", decode_identity):
        return "cached"
    pairs = np.loadtxt(out / "pairs.tsv", dtype=int, ndmin=2)
    positions = np.loadtxt(out / "positions.txt", dtype=int, ndmin=1)
    raw = directory / "posteriors.zst"
    result = run_within_decoder(decoder, directory / "decoded_input.trees", directory / "native_mean.tsv",
        raw_output=raw, threshold_years=cfg["tmrca_cutoffs_years"], generation_time=cfg["generation_time_years"],
        mutation_rate=cfg["mutation_rate"], scaled_mutation_rate=cfg["decoder_scaled_mutation_rate"],
        recombination_to_mutation_ratio=cfg["decoder_recombination_to_mutation_ratio"],
        output_at_stride=-1, output_at_hets=False, only_within=False,
        output_positions_file=out / "positions.txt", pairs_file=out / "pairs.tsv",
        vcf_position_transform="one_based", recent_call="mean", threads=1,
        cache_size=cfg["decoder_cache_size"], pair_block=cfg["decoder_pair_block"],
        exp10=cfg["decoder_exp10"], backward_alignment=cfg["decoder_backward_alignment"],
        extra_args=["--no_recent_probability"])
    for channel in ("stdout", "stderr"):
        (directory / f"decoder.{channel}.log").write_text(result.pop(channel))
    meta = json.loads(Path(str(raw) + ".meta").read_text())
    np.testing.assert_array_equal(meta["pairs"], pairs)
    np.testing.assert_array_equal(meta["output_positions"], positions)
    alpha, beta = posterior_arrays(raw, len(pairs), len(positions))
    classes = pair_classes(np.load(directory / "focal_carriers.npy", allow_pickle=False), pairs)
    ts = tskit.load(directory / "decoded_input.trees")
    frames = []
    for source, values in (("decoded", summarize_decoded(alpha, beta, classes, cfg)),
                           ("truth", summarize_truth(ts, pairs, classes, positions, cfg))):
        counts, sums = values
        frames.append(build_frame(counts, sums, classes, positions, cfg, source))
    profile = pd.concat(frames, ignore_index=True)
    native = pd.read_csv(directory / "native_mean.tsv", sep="\t")
    cols = [f"frac_recent_{t}" for t in cfg["tmrca_cutoffs_years"]]
    observed = profile[(profile.source == "decoded") & (profile.pair_class == "all pairs")]
    np.testing.assert_allclose(observed[cols], native[cols], atol=1e-7, rtol=0)
    values = profile[cols].to_numpy()
    assert ((values[np.isfinite(values)] >= 0) & (values[np.isfinite(values)] <= 1)).all()
    assert np.all(np.nan_to_num(np.diff(values, axis=1)) >= 0)
    profile.to_csv(directory / "profiles.csv", index=False)
    legacy.atomic_json(directory / "decoded.json", dict(identity=decode_identity, decoder=result,
        artifacts=artifact_specs(directory, ("profiles.csv", "native_mean.tsv", "posteriors.zst", "posteriors.zst.meta"))))
    return "decoded"


def worker(payload):
    cfg, task, runtime, out_path, slim, decoder, phase, legacy_root = payload
    directory = Path(out_path) / "regions" / task["id"]
    directory.mkdir(parents=True, exist_ok=True)
    os.environ["TMPDIR"] = str(directory)
    import tempfile
    tempfile.tempdir = str(directory)
    import resource
    limit = int(cfg["memory_limit_gb"] * 1e9 / cfg["workers"])
    resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
    identity = simulation_identity(cfg, task, runtime)
    started = time.time()
    legacy.atomic_json(directory / "status.json", dict(task=task, phase=phase, state="running", pid=os.getpid(), started=started))
    try:
        usage = shutil.disk_usage(directory)
        if usage.free < 20_000_000_000 or usage.used >= cfg["storage_limit_bytes"]:
            raise RuntimeError("Storage budget/reserve reached; completed artifacts remain reusable")
        with (directory / f"{phase}.stdout.log").open("a") as stdout, (directory / f"{phase}.stderr.log").open("a") as stderr:
            with redirect_stdout(stdout), redirect_stderr(stderr):
                if phase == "simulate":
                    state = simulate(cfg, task, directory, identity, slim, legacy_root)
                elif phase == "audit":
                    record = read_receipt(directory, "simulation.json", identity)
                    if record is None:
                        raise ValueError(f"Missing simulation: {directory}")
                    validate_trees(directory, cfg, record)
                    state = "audited"
                else:
                    state = decode(cfg, task, directory, identity, Path(out_path), decoder)
        result = dict(task_id=task["id"], phase=phase, state=state, seconds=time.time() - started)
    except Exception:
        result = dict(task_id=task["id"], phase=phase, state="failed", seconds=time.time() - started, error=traceback.format_exc())
    legacy.atomic_json(directory / "status.json", result)
    return result


def initialise(cfg, out, slim, decoder):
    assert cfg["population"] == "EAS" and cfg["recent_call"] == "mean"
    assert cfg["target_replicates"] > 0 and cfg["null_replicates"] > 0
    assert 1 <= cfg["workers"] <= 24 and cfg["memory_limit_gb"] > 0
    assert cfg["tmrca_cutoffs_years"] == sorted(set(cfg["tmrca_cutoffs_years"]))
    assert cfg["scaling_by_origin"]["de_novo_EAS"] == 1
    runtime = dict(versions={p: importlib.metadata.version(p) for p in ("msprime", "numpy", "pyslim", "stdpopsim", "tskit")},
                   slim_sha256=legacy.digest(slim))
    assert runtime["versions"]["stdpopsim"] == "0.3.0"
    plan = tasks(cfg)
    assert len({t["id"] for t in plan}) == len(plan)
    assert len({seed_for(cfg, t) for t in plan}) == len(plan)
    possible = np.column_stack(np.triu_indices(2 * cfg["sample_diploids"], k=1))
    choices = np.random.default_rng(cfg["pairs_seed"]).choice(len(possible), cfg["haplotype_pairs"], replace=False)
    pairs = possible[np.sort(choices)]
    positions = np.arange(0, cfg["scored_length_bp"], cfg["stride_bp"])
    assert cfg["focal_position_bp"] in positions
    for name, values in (("pairs.tsv", pairs), ("positions.txt", positions)):
        path = out / name
        if path.exists():
            np.testing.assert_array_equal(np.loadtxt(path, dtype=int), values)
        else:
            np.savetxt(path, values, fmt="%d", delimiter="\t")
    legacy.versioned_json(out / "manifest.json", dict(config=cfg, runtime=runtime, tasks=plan,
        python=sys.version, repo=str(REPO), output=str(out), slim=slim, decoder=decoder,
        pairs_sha256=legacy.digest(out / "pairs.tsv"), positions_sha256=legacy.digest(out / "positions.txt"),
        source_sha256=legacy.digest(Path(__file__)), decoder_sha256=legacy.digest(decoder)))
    return runtime, plan


def run_phase(cfg, out, slim, decoder, phase, runtime, plan, legacy_root):
    def status(results, active, state):
        legacy.atomic_json(out / f"{phase}_status.json", dict(phase=phase, state=state, target=len(plan),
            completed=len(results), failed=[r for r in results if r["state"] == "failed"],
            active=list(active), updated=time.time(), workers=cfg["workers"], results=results))
    remaining = iter(plan)
    results = []
    futures = {}
    with ProcessPoolExecutor(max_workers=cfg["workers"]) as pool:
        def submit_next():
            task = next(remaining, None)
            if task is not None:
                payload = (cfg, task, runtime, str(out), slim, decoder, phase, legacy_root)
                futures[pool.submit(worker, payload)] = task["id"]
        for _ in range(cfg["workers"]):
            submit_next()
        status(results, futures.values(), "running")
        failed = False
        while futures:
            done, _ = wait(futures, timeout=30, return_when=FIRST_COMPLETED)
            for future in done:
                task_id = futures.pop(future)
                try:
                    result = future.result()
                except Exception:
                    result = dict(task_id=task_id, state="failed", error=traceback.format_exc())
                results.append(result)
                print(json.dumps(result), flush=True)
                failed |= result["state"] == "failed"
            if not failed:
                while len(futures) < cfg["workers"]:
                    count = len(futures)
                    submit_next()
                    if len(futures) == count:
                        break
            status(results, futures.values(), "stopping_after_failure" if failed else "running")
    status(results, [], "failed" if failed else "complete")
    if failed:
        raise RuntimeError(f"{phase} failed; see {out / (phase + '_status.json')}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--slim", required=True)
    parser.add_argument("--decoder", required=True)
    parser.add_argument("--legacy-root")
    parser.add_argument("--phase", choices=("plan", "simulate", "audit", "decode", "analyze", "run"), required=True)
    parser.add_argument("--task", action="append", help="Exact task ID; repeat for a bounded smoke run")
    args = parser.parse_args()
    cfg = json.loads(args.config.read_text())
    args.out.mkdir(parents=True, exist_ok=True)
    import fcntl
    with (args.out / "study.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        runtime, plan = initialise(cfg, args.out, args.slim, args.decoder)
        if args.task:
            keep = set(args.task)
            plan = [task for task in plan if task["id"] in keep]
            if len(plan) != len(keep):
                raise ValueError("Unknown task ID")
        if args.phase == "plan":
            print(json.dumps(dict(planned=len(plan), output=str(args.out)), indent=2))
            return
        phases = ("simulate", "audit", "decode", "analyze") if args.phase == "run" else (args.phase,)
        began=time.time()
        state=dict(command_phase=args.phase, target=len(plan), pid=os.getpid(), started=began)
        try:
            for phase in phases:
                legacy.atomic_json(args.out/"run_status.json",dict(state,phase=phase,state="running",updated=time.time()))
                if phase == "analyze":
                    if args.task:
                        raise ValueError("Analysis requires the entire planned cohort")
                    from .origin_onset_analysis import analyze
                    analyze(args.out, cfg, plan)
                else:
                    run_phase(cfg, args.out, args.slim, args.decoder, phase, runtime, plan, args.legacy_root)
        except Exception:
            legacy.atomic_json(args.out/"run_status.json",dict(state,phase=phase,state="failed",updated=time.time(),error=traceback.format_exc()))
            raise
        legacy.atomic_json(args.out/"run_status.json",dict(state,phase=args.phase,state="complete",updated=time.time()))


if __name__ == "__main__":
    main()
