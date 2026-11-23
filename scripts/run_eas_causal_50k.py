"""Simulate the EAS strength grid at the causal coordinate, retaining fixation."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import os
from pathlib import Path
import shutil
import tempfile
import time

import run_eas_strength_series as series
import numpy as np

arrays = series.arrays
base = arrays.base
RAW = Path.home() / "sim/eas_causal_50k_20261123"
OUT = series.OUT.parent / "eas_causal_50k_20261123"


def prepare(raw):
    old = json.loads((series.OUT / "simulation_manifest.json").read_text())
    configs = {}
    for window, previous in old["configs"].items():
        cfg = dict(previous, allele_age_years=50000, retain_fixation=True,
                   seed_namespace=f"eas-causal-50k-v1:{window}:",
                   selected_condition="population focal copies > 0 at present, including fixation",
                   null_condition="age-matched 50 kya neutral bank", stride=None,
                   source_sha256=base.digest(arrays.__file__))
        cfg["slim"] = os.environ.get("EAS_SLIM", cfg["slim"] if Path(cfg["slim"]).exists()
                                     else str(base.REPO/".native/bin/slim"))
        cfg["decoder"] = os.environ.get("EAS_DECODER", str(base.REPO/"bin/gamma_smc"))
        for field in ("slim", "decoder"):
            assert base.digest(cfg[field]) == cfg[field + "_sha256"]
        assert base.digest(base.REPO / "no_introgression_EAS_sim/resources/EAS.npz") == cfg["phlash_sha256"]
        configs[window] = cfg
    manifest = dict(configs=configs, target=100, strengths=series.STRENGTHS, focal_bp=500000,
                    origins=["de_novo", "introgression"], source_bank=base.digest(series.OUT/"mass.npz"),
                    pairing="Common seed per attempt across strengths within each origin, window and independent replicate slot",
                    estimator="1000 fixed haplotype pairs, posterior-weighted cumulative TMRCA at 1..50 kya",
                    founder_rule="One independent de novo ancestry per slot; fresh introgression founder per proposal")
    raw.mkdir(parents=True, exist_ok=True)
    path = raw/"manifest.json"
    if path.exists():
        assert json.loads(path.read_text()) == json.loads(json.dumps(manifest))
    else:
        base.save(path, manifest)
    pairs = raw/"pairs.tsv"
    if not pairs.exists():
        possible = np.column_stack(np.triu_indices(200, 1))
        subset = np.sort(np.random.default_rng(1729).choice(len(possible), 1000, replace=False))
        np.savetxt(pairs, possible[subset], fmt="%d", delimiter="\t")
    assert base.digest(pairs) == old["input_sha256"]["pairs.tsv"]
    positions = raw/"positions.txt"
    if not positions.exists():
        positions.write_text("500000\n")
    np.testing.assert_array_equal(np.loadtxt(positions, ndmin=1), [500000])
    return configs


def group(raw, origin, window, slot, cfg):
    work = raw/"groups"/origin/window/f"{slot:03d}"
    work.mkdir(parents=True, exist_ok=True)
    prefix = f"{origin}/{window}/slot_{slot:03d}"
    status = work/"progress.json"
    attempt = json.loads(status.read_text())["next_attempt"] if status.exists() else 0
    pending = []
    for s in series.STRENGTHS:
        destination = work/f"s{s:g}"
        if (destination/"complete.json").exists():
            rec = json.loads((destination/"complete.json").read_text())
            for file, field in [("trajectory.tsv","trajectory_sha256"),("scores.npz","scores_sha256")]:
                assert base.digest(destination/file) == rec[field]
        else:
            pending.append(s)
    if not pending:
        return
    cached = series.grid.cached
    old_onset = cached.ORIGINAL_ONSET
    try:
        if origin == "de_novo":
            founder = work/"founder"
            founder.mkdir(exist_ok=True)
            label = cfg["seed_namespace"] + prefix + ":founder"
            if not (founder/"founder.json").exists():
                ids, sizes = old_onset(dict(cfg, ancestry_start_generations=2000), origin, label, founder)
                assert len(ids) == 1
                base.save(founder/"founder.json", dict(carrier_ids=ids, sizes=sizes, age_generations=2000,
                          label=label, onset_sha256=base.digest(founder/"onset.trees")))
            meta = json.loads((founder/"founder.json").read_text())
            assert meta["onset_sha256"] == base.digest(founder/"onset.trees")
        while pending:
            unit = f"{prefix}/attempt_{attempt:07d}"
            with tempfile.TemporaryDirectory(dir=work) as temp:
                proposal = Path(temp)
                if origin == "introgression":
                    founder = proposal/"founder"
                    founder.mkdir()
                    ids, sizes = old_onset(dict(cfg, ancestry_start_generations=2000), origin,
                                           cfg["seed_namespace"] + unit, founder)
                else:
                    ids, sizes = meta["carrier_ids"], meta["sizes"]
                def install(config, pop, seed, directory):
                    for name in ("onset.trees", "onset_provenance.json"):
                        (directory/name).symlink_to(founder/name)
                    return ids, sizes
                arrays.ancestor.onset = install
                for s in pending.copy():
                    directory = proposal/f"s{s:g}"
                    directory.mkdir()
                    task = dict(unit=unit, origin=origin, cohort="selected", selection=s)
                    record = arrays.simulate(cfg, task, directory)
                    if not record["accepted"]:
                        assert record["rejection"] == "loss"
                        # Identical random state and s=0 in every arm through the recent onset.
                        if window == "10k" and record["stop_years_ago"] >= 10000:
                            break
                        continue
                    assert record["trajectory_rows"] == 2001 and record["stop_years_ago"] == 0
                    record.update(base.score(cfg, directory, raw))
                    record.update(task, slot=slot, attempt=attempt, window=window, position=500000)
                    destination = work/directory.name
                    destination.mkdir(exist_ok=True)
                    for name in ("scores.npz","trajectory.tsv","observed.trees","model.slim","onset_provenance.json"):
                        shutil.copyfile(directory/name, destination/name)
                    base.save(destination/"complete.json", record)
                    pending.remove(s)
                attempt += 1
                base.save(status, dict(next_attempt=attempt, remaining=pending, updated=time.time()))
    finally:
        arrays.ancestor.onset = old_onset


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--raw", type=Path, default=RAW)
    p.add_argument("--workers", type=int, default=4, choices=range(1,21))
    a = p.parse_args()
    raw = a.raw.resolve()
    cfgs = prepare(raw)
    for src in (Path(__file__), Path(arrays.__file__), Path(base.__file__), Path(arrays.ancestor.__file__)):
        dest = raw/"executed_sources"/(src.stem+"_"+base.digest(src)[:12]+".py")
        dest.parent.mkdir(exist_ok=True)
        if not dest.exists():
            shutil.copyfile(src, dest)
    os.environ["LD_LIBRARY_PATH"] = str(Path(cfgs["10k"]["slim"]).parent.parent/"lib")+":"+os.environ.get("LD_LIBRARY_PATH","")
    print(f"Started {time.ctime()} | {raw} | {a.workers} workers | 7200 selected histories", flush=True)
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        jobs = [pool.submit(group, raw, origin, win, slot, cfg)
                for slot in range(100) for win,cfg in cfgs.items() for origin in ("introgression","de_novo")]
        for i,j in enumerate(as_completed(jobs), 1):
            j.result()
            print(f"{time.ctime()} | {i}/600 groups finished", flush=True)


if __name__ == "__main__":
    main()
