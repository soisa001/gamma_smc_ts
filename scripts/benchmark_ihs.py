"""Resumable iHS benchmark on saved haplotypes; no simulation or decoding."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import fcntl
import hashlib
import importlib.metadata
import inspect
import json
import os
from pathlib import Path
import time

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[name] = "1"
os.environ["MPLBACKEND"] = "Agg"

import allel
import numpy as np
import pandas as pd
import tskit

REPO = Path(__file__).resolve().parents[1]
METHODS = ("nearest_core_abs_ihs", "focal_100kb_extreme_fraction")


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def json_write(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w") as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def ordered_nodes(ts):
    samples = set(map(int, ts.samples()))
    groups = [[int(n) for n in individual.nodes if int(n) in samples] for individual in ts.individuals()]
    groups = [g for g in groups if g]
    if any(len(g) != 2 for g in groups):
        raise ValueError("Non-diploid sample ordering")
    nodes = np.array([n for g in groups for n in g], dtype=int)
    if len(nodes) != len(samples) or set(nodes) != samples:
        raise ValueError("Sample ordering is incomplete")
    return nodes


def haplotypes(ts):
    nodes = ordered_nodes(ts)
    h = np.empty((ts.num_sites, len(nodes)), dtype=np.int8)
    positions = np.empty(ts.num_sites, dtype=np.int64)
    counts = np.empty(ts.num_sites, dtype=np.int32)
    kept = 0
    for variant in ts.variants(samples=nodes):
        if len(variant.alleles) != 2:
            continue
        g = variant.genotypes
        if np.any(g < 0) or np.any(g > 1):
            raise ValueError("Missing or nonbinary genotypes at biallelic site")
        count = int(g.sum())
        if count == 0 or count == len(nodes):
            continue
        if int(variant.site.position) != variant.site.position:
            raise ValueError("Noninteger variant coordinate")
        h[kept] = g
        positions[kept] = int(variant.site.position)
        counts[kept] = count
        kept += 1
    if kept < 2 or not np.all(np.diff(positions[:kept]) > 0):
        raise ValueError("Insufficient or unordered segregating sites")
    return h[:kept], positions[:kept], counts[:kept]


def raw_ihs(h, positions, cfg):
    with np.errstate(divide="ignore", invalid="ignore"):
        return allel.ihs(h, positions, map_pos=positions*cfg["recombination_rate"],
            min_ehh=cfg["min_ehh"], min_maf=cfg["min_maf"],
            include_edges=cfg["include_edges"], gap_scale=cfg["gap_scale"],
            max_gap=cfg["max_gap"], use_threads=False)


def raw_identity(item, cfg):
    settings = {k: cfg[k] for k in ("recombination_rate", "min_ehh", "min_maf", "include_edges", "gap_scale", "max_gap")}
    implementation = "\n".join(inspect.getsource(f) for f in (ordered_nodes, haplotypes, raw_ihs))
    return dict(schema="ihs-raw/v1", simulation_sha256=item["simulation_sha256"],
        settings=settings, versions={m: importlib.metadata.version(m) for m in ("scikit-allel", "numpy", "tskit")},
        implementation_sha256=hashlib.sha256(implementation.encode()).hexdigest())


def raw_one(payload):
    item, cfg, out = payload
    directory = out / "raw" / item["task_id"]
    directory.mkdir(parents=True, exist_ok=True)
    receipt_path = directory / "complete.json"
    identity = raw_identity(item, cfg)
    if receipt_path.exists():
        receipt = json.loads(receipt_path.read_text())
        if receipt["identity"] != identity:
            raise ValueError(f"Raw input/config changed: {directory}")
        if digest(directory / "ihs.npz") != receipt["output_sha256"]:
            raise ValueError(f"Corrupt cached iHS: {directory}")
        return dict(task_id=item["task_id"], status="cached", seconds=receipt["seconds"])
    started = time.monotonic()
    if digest(item["simulation_path"]) != item["simulation_sha256"]:
        raise ValueError(f"Corrupt simulation: {item['task_id']}")
    ts = tskit.load(item["simulation_path"])
    if ts.num_samples != 400 or ts.sequence_length != 11_000_000:
        raise ValueError("Expected 400 haplotypes across 11 Mb")
    h, positions, counts = haplotypes(ts)
    raw = raw_ihs(h, positions, cfg)
    if np.isinf(raw).any():
        raise ValueError("Infinite raw iHS")
    eligible = np.minimum(counts, h.shape[1]-counts) / h.shape[1] >= cfg["min_maf"]
    if np.any(np.isfinite(raw[~eligible])):
        raise ValueError("iHS unexpectedly scored a core below minimum MAF")
    temporary = directory / "ihs.npz.tmp"
    with temporary.open("wb") as stream:
        np.savez_compressed(stream, positions=positions[eligible], counts=counts[eligible], raw=raw[eligible])
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, directory / "ihs.npz")
    receipt = dict(identity=identity, task_id=item["task_id"], all_segregating_biallelic_sites=len(positions),
        eligible_core_sites=int(eligible.sum()), finite_core_sites=int(np.isfinite(raw).sum()),
        excluded_edge_cores=int((eligible & ~np.isfinite(raw)).sum()), seconds=time.monotonic()-started,
        output_sha256=digest(directory / "ihs.npz"), source_tree=item["simulation_path"])
    json_write(receipt_path, receipt)
    return dict(task_id=item["task_id"], status="computed", seconds=receipt["seconds"],
                cores=receipt["finite_core_sites"])


def bin_ids(counts, bins):
    return np.minimum(counts*bins//400, bins-1)


def moments(counts, values, bins):
    ids = bin_ids(counts, bins)
    return np.array([np.bincount(ids, minlength=bins),
                     np.bincount(ids, weights=values, minlength=bins),
                     np.bincount(ids, weights=values**2, minlength=bins)])


def load_raw(payload):
    item, cfg, out = payload
    directory = out / "raw" / item["task_id"]
    receipt = json.loads((directory / "complete.json").read_text())
    if receipt["identity"] != raw_identity(item, cfg) or digest(directory / "ihs.npz") != receipt["output_sha256"]:
        raise ValueError(f"Bad raw receipt: {directory}")
    simulation = json.loads((Path(item["simulation_path"]).parent / "simulated.json").read_text())
    with np.load(directory / "ihs.npz", allow_pickle=False) as saved:
        pos = saved["positions"] - simulation["crop_offset"]
        keep = (pos >= 0) & (pos < cfg["scored_length_bp"])
        pos, counts, raw = pos[keep], saved["counts"][keep], saved["raw"][keep]
    finite = np.isfinite(raw)
    focal = (pos == cfg["focal_position_bp"]) & finite
    return dict(pos=pos[finite], counts=counts[finite], raw=raw[finite],
        moments=moments(counts[finite], raw[finite], cfg["frequency_bins"]),
        cores=int(keep.sum()), finite=int(finite.sum()), focal_valid=bool(focal.any()),
        focal_af=float(simulation["sample_af"]), raw_receipt_sha256=digest(directory / "complete.json"))


def score_region(pos, z, cfg):
    absolute = np.abs(z)
    maximum = float(np.max(absolute, initial=0))
    starts = np.arange(0, cfg["scored_length_bp"]-cfg["window_bp"]+1, cfg["window_stride_bp"])
    left = np.searchsorted(pos, starts)
    right = np.searchsorted(pos, starts+cfg["window_bp"])
    number = right-left
    extreme = np.r_[0, np.cumsum(absolute > cfg["extreme_abs_ihs"])]
    counts = extreme[right]-extreme[left]
    fractions = np.divide(counts, number, out=np.zeros(len(starts)), where=number >= cfg["minimum_window_sites"])
    peak = float(fractions.max(initial=0))
    nearest = int(np.argmin(abs(pos-cfg["focal_position_bp"]))) if len(pos) else None
    focal_max = (float(absolute[nearest]) if nearest is not None and
                 abs(pos[nearest]-cfg["focal_position_bp"]) <= cfg["core_radius_bp"] else 0.)
    center = cfg["focal_position_bp"]-cfg["window_bp"]//2
    focal_window = float(fractions[starts == center][0])
    return np.array([maximum, peak]), np.array([focal_max, focal_window])


def wilson(k, n):
    z = 1.959963984540054
    c = (k/n+z*z/(2*n))/(1+z*z/n)
    r = z*np.sqrt(k/n*(1-k/n)/n+z*z/(4*n*n))/(1+z*z/n)
    return float(c-r), float(c+r)


def analyze(items, cfg, out, fold_path, workers):
    reference_manifest = json.loads((fold_path.parent / "artifact_manifest.json").read_text())
    if digest(fold_path) != reference_manifest[fold_path.name]["sha256"]:
        raise ValueError("Corrupt reference folds")
    old_folds = pd.read_csv(fold_path).set_index("key").fold.to_dict()
    folds = np.array([old_folds[i["task_id"] if i["mode"] == "neutral" else f'onset50000/rep{i["replicate"]:04d}'] for i in items])
    neutral = np.array([i["mode"] == "neutral" for i in items])
    with ThreadPoolExecutor(max_workers=workers) as pool:
        loaded = list(pool.map(load_raw, [(i, cfg, out) for i in items]))
    bin_moments = np.array([data["moments"] for data in loaded])
    scores = np.zeros((5, len(items), len(METHODS)))
    focal_scores = np.zeros_like(scores)
    normalizations = []
    for fold in range(5):
        fit = neutral & np.isin(folds, [(fold+1)%5, (fold+2)%5])
        cal = neutral & np.isin(folds, [(fold+3)%5, (fold+4)%5])
        test = neutral & (folds == fold)
        assert (int(fit.sum()), int(cal.sum()), int(test.sum())) == (400, 400, 200)
        assert not np.any(fit & cal) and not np.any(fit & test) and not np.any(cal & test)
        n, sums, squares = bin_moments[fit].sum(axis=0)
        mean = np.divide(sums, n, out=np.zeros_like(n), where=n>0)
        variance = np.divide(squares, n, out=np.ones_like(n), where=n>0)-mean**2
        std = np.sqrt(np.maximum(variance, 0))
        for j, data in enumerate(loaded):
            ids = bin_ids(data["counts"], cfg["frequency_bins"])
            if np.any(n[ids] < cfg["minimum_normalization_sites"]) or np.any(std[ids] <= 0):
                raise ValueError("Insufficient neutral frequency-bin normalization")
            z = (data["raw"]-mean[ids])/std[ids]
            scores[fold, j], focal_scores[fold, j] = score_region(data["pos"], z, cfg)
        normalizations.append(dict(fold=fold, fit_ids=[items[k]["task_id"] for k in np.flatnonzero(fit)],
            calibration_ids=[items[k]["task_id"] for k in np.flatnonzero(cal)],
            test_neutral_ids=[items[k]["task_id"] for k in np.flatnonzero(test)],
            counts=n.astype(int).tolist(), means=mean.tolist(), stds=std.tolist()))
    records = []
    for fold in range(5):
        cal = neutral & np.isin(folds, [(fold+3)%5, (fold+4)%5])
        test_ids = np.flatnonzero(folds == fold)
        for k, method in enumerate(METHODS):
            reference = focal_scores[fold, cal, k]
            endpoint = ("positional", "focal_100kb_window")[k]
            for j in test_ids:
                value = focal_scores[fold, j, k]
                p = float((1+np.count_nonzero(reference >= value))/401)
                for alpha in cfg["alphas"]:
                    records.append(dict(task_id=items[j]["task_id"], s=items[j]["s"],
                        mode=items[j]["mode"], fold=fold, endpoint=endpoint, method=method,
                        alpha=alpha, score=value, p=p, called=bool(p <= alpha)))
    predictions = pd.DataFrame(records)
    predictions.to_csv(out / "predictions.csv.gz", index=False, compression=dict(method="gzip", mtime=0))
    metrics = []
    for (endpoint, method, alpha, s), group in predictions[predictions["mode"] == "selected"].groupby(["endpoint", "method", "alpha", "s"]):
        k, n = int(group.called.sum()), len(group)
        background = predictions[(predictions["mode"] == "neutral") & (predictions.method == method) & (predictions.alpha == alpha)]
        assert len(background) == 1000
        lo, hi = wilson(k, n)
        metrics.append(dict(endpoint=endpoint, method=method, alpha=alpha, s=s, selected_regions=n,
            selected_called=k, power=k/n, wilson_low=lo, wilson_high=hi,
            neutral_regions=len(background), neutral_called=int(background.called.sum()),
            local_fpr=float(background.called.mean())))
    metrics = pd.DataFrame(metrics)
    metrics.to_csv(out / "metrics.csv", index=False)
    diagnostics = pd.DataFrame([dict(task_id=item["task_id"], mode=item["mode"], s=item["s"],
        replicate=item["replicate"], fold=int(folds[j]), eligible_cores=d["cores"], finite_cores=d["finite"],
        focal_valid=d["focal_valid"], focal_af=d["focal_af"]) for j, (item, d) in enumerate(zip(items, loaded))])
    diagnostics.to_csv(out / "regions.csv", index=False)
    np.savez_compressed(out / "scores.npz", scores=scores, focal_scores=focal_scores, folds=folds)
    json_write(out / "normalization.json", normalizations)
    # Independent rank reconstruction from persisted scores and predictions.
    with np.load(out / "scores.npz", allow_pickle=False) as saved:
        stored_focal = saved["focal_scores"]
    check = pd.read_csv(out / "predictions.csv.gz")
    lookup = {item["task_id"]: j for j, item in enumerate(items)}
    for (fold, endpoint, method), group in check.groupby(["fold", "endpoint", "method"]):
        k = METHODS.index(method)
        ids = [lookup[t] for t in group.task_id]
        cal = neutral & np.isin(folds, [(fold+3)%5, (fold+4)%5])
        v = stored_focal[fold, ids, k]
        expected = (1+(stored_focal[fold, cal, k, None] >= v).sum(axis=0))/401
        np.testing.assert_allclose(expected, group.p, atol=1e-14, rtol=0)
        np.testing.assert_array_equal(expected <= group.alpha, group.called)
    assert np.all(focal_scores <= scores+1e-14)
    json_write(out / "audit.json", dict(status="passed", regions=len(items), neutral_regions=int(neutral.sum()),
        selected_regions=int((~neutral).sum()), folds_disjoint=True, selected_excluded_from_normalization=True,
        prediction_ranks_recomputed=True, focal_scores_bounded_by_region_scores=True,
        zero_or_undefined_local_scores_retained=True, minimum_p=1/401,
        null_is_same_local_position_or_window=True, regional_maximum_used_for_pvalues=False))
    json_write(out / "analysis_provenance.json", dict(config=cfg, input_manifest_sha256=digest(out / "sample_manifest.json"),
        source_sha256=digest(Path(__file__)), fold_reference_sha256=digest(fold_path),
        raw_receipt_hashes={item["task_id"]: d["raw_receipt_sha256"] for item, d in zip(items, loaded)},
        versions={m: importlib.metadata.version(m) for m in ("scikit-allel", "numpy", "pandas", "tskit")},
        simulations_started=False, gamma_smc_decoding_started=False,
        selected_cohort_incomplete=any(sum(i["mode"] == "selected" and i["s"] == s for i in items) != 100
            for s in sorted({i["s"] for i in items if i["mode"] == "selected"}))))
    names = ("predictions.csv.gz", "metrics.csv", "regions.csv", "scores.npz", "normalization.json", "audit.json", "analysis_provenance.json")
    json_write(out / "artifact_manifest.json", {name: dict(sha256=digest(out/name), bytes=(out/name).stat().st_size) for name in names})
    print(metrics[metrics.alpha == .05].to_string(index=False), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, default=Path("/mnt/d/phase2simselection/sim/eas_q02_h400"))
    parser.add_argument("--out", type=Path, default=Path("/mnt/d/phase2simselection/sim/eas_ihs_h400"))
    parser.add_argument("--config", type=Path, default=REPO / "configs/ihs_eas_h400.json")
    parser.add_argument("--folds", type=Path, default=Path("/mnt/d/phase2simselection/sim/eas_allele_class_ablation/regions.csv"))
    parser.add_argument("--phase", choices=("smoke", "raw", "analyse", "run"), default="run")
    args = parser.parse_args()
    cfg = json.loads(args.config.read_text())
    args.out.mkdir(parents=True, exist_ok=True)
    inventory = pd.read_csv(args.study / "simulation_inventory.csv")
    status = json.loads((args.study / "simulation_status.json").read_text())
    if status["completed"] != len(inventory) or status["failed"]:
        raise ValueError("Simulation inventory differs from validated status")
    items = inventory.to_dict("records")
    workers = cfg["workers"]
    if not 1 <= workers <= 24:
        raise ValueError("Local workers must be between 1 and 24")
    with (args.out / "runner.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        json_write(args.out / "sample_manifest.json", items)
        if args.phase != "analyse":
            tasks = items
            if args.phase == "smoke":
                tasks = [i for i in items if (i["mode"] == "neutral" and i["replicate"] < 2)
                    or (i["mode"] == "selected" and i["replicate"] == 0 and i["s"] in (.001, .005))
                    or (i["mode"] == "selected" and i["replicate"] == 1 and i["s"] == .010)]
            records = []
            started = time.monotonic()
            with ProcessPoolExecutor(max_workers=workers) as pool:
                futures = [pool.submit(raw_one, (i, cfg, args.out)) for i in tasks]
                for future in as_completed(futures):
                    record = future.result()
                    records.append(record)
                    print(json.dumps(record), flush=True)
                    json_write(args.out / ("smoke_status.json" if args.phase == "smoke" else "raw_status.json"),
                        dict(completed=len(records), target=len(tasks), workers=workers,
                             state="complete" if len(records) == len(tasks) else "running", seconds=time.monotonic()-started))
        if args.phase in ("analyse", "run"):
            analyze(items, cfg, args.out, args.folds, workers)


if __name__ == "__main__":
    main()
