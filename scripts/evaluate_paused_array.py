"""Evaluate saved selected-site truth against existing neutral region maxima.

This is a lower bound on full regional power, not a replacement for a blind
regional scan. It runs neither simulations nor decoding. All outputs go to the
specified directory; the simulation archive and old analyses are read-only.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import json
import os
from pathlib import Path
import time

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "1"
os.environ["MPLBACKEND"] = "Agg"

import numpy as np
import pandas as pd
from scipy.stats import binomtest, spearmanr
import tskit

from gamma_smc_aou.fresh_power import atomic_json, digest, ordered_nodes


def ancestral_groups(tree, nodes, cutoff_generations):
    """Same group iff pair TMRCA is strictly below the cutoff (not <=)."""
    labels = []
    for sample in nodes:
        node = int(sample)
        parent = tree.parent(node)
        while parent != tskit.NULL and tree.time(parent) < cutoff_generations:
            node = parent
            parent = tree.parent(node)
        labels.append(node)
    return np.asarray(labels)


def verify_files(directory, names):
    manifest = json.loads((directory / "artifact_manifest.json").read_text())
    for name in names:
        if digest(directory / name) != manifest[name]["sha256"]:
            raise ValueError(f"Corrupt reference artifact: {directory / name}")
    return {name: manifest[name]["sha256"] for name in names}


def focal_one(payload):
    row, cfg, pairs, out = payload
    cache = out / "cache" / (row["task_id"].replace("/", "_") + ".json")
    identity = {k: row[k] for k in ("task_id", "decoder_input_sha256", "carriers_sha256")}
    identity.update(pairs_sha256=digest(Path(row["decoder_input_path"]).parents[2] / "pairs.tsv"),
                    source_sha256=digest(Path(__file__)), cutoffs=cfg["tmrca_cutoffs_years"],
                    generation_time_years=cfg["generation_time_years"],
                    focal_position_bp=cfg["focal_position_bp"])
    if cache.exists():
        saved = json.loads(cache.read_text())
        if saved["identity"] != identity:
            raise ValueError(f"Changed focal cache inputs: {cache}")
        if digest(cache) != json.loads(cache.with_suffix(".sha256.json").read_text())["sha256"]:
            raise ValueError(f"Corrupt focal cache: {cache}")
        return saved["result"]
    ts = tskit.load(row["decoder_input_path"])
    nodes = ordered_nodes(ts)
    carriers = np.load(row["carriers_path"], allow_pickle=False)
    if ts.num_samples != 400 or len(pairs) != cfg["haplotype_pairs"]:
        raise ValueError("Unexpected sample or pair panel")
    tree = ts.at(cfg["focal_position_bp"])
    ages = np.fromiter((tree.tmrca(int(nodes[a]), int(nodes[b])) for a, b in pairs), float)
    ages *= cfg["generation_time_years"]
    aa = carriers[pairs].all(axis=1)
    record = dict(task_id=row["task_id"], s=row["s"], replicate=int(row["replicate"]),
                  sample_af=float(carriers.mean()), n_alt=int(carriers.sum()),
                  n_alt_alt=int(aa.sum()), n_pairs=len(pairs),
                  accepted_attempt=int(row["accepted_attempt"]), fixed=bool(carriers.all()))
    if record["sample_af"] != row["sample_af"]:
        raise ValueError("Carrier mask differs from verified inventory")
    for cutoff in cfg["tmrca_cutoffs_years"]:
        recent = ages < cutoff
        record[f"all_{cutoff}"] = float(recent.mean())
        record[f"mass_{cutoff}"] = float((recent & aa).mean())
        record[f"aa_{cutoff}"] = float(recent[aa].mean()) if aa.any() else None
    groups = ancestral_groups(tree, nodes, cfg["pulse_years"] / cfg["generation_time_years"])
    same = groups[pairs[:, 0]] == groups[pairs[:, 1]]
    np.testing.assert_array_equal(same, ages < cfg["pulse_years"])
    _, sizes = np.unique(groups[carriers], return_counts=True)
    weights = sizes / sizes.sum()
    record.update(carrier_lineages_at_50k=int(len(sizes)),
                  effective_carrier_lineages_at_50k=float(1 / np.sum(weights**2)),
                  largest_carrier_lineage_fraction=float(weights.max()),
                  exact_all_pairs_carrier_fraction_50k=(float(np.sum(sizes * (sizes - 1)) /
                    (carriers.sum() * (carriers.sum() - 1))) if carriers.sum() > 1 else None))
    cache.parent.mkdir(parents=True, exist_ok=True)
    atomic_json(cache, dict(identity=identity, result=record))
    atomic_json(cache.with_suffix(".sha256.json"), dict(sha256=digest(cache)))
    return record


def wilson(k, n):
    z = 1.959963984540054
    center = (k / n + z*z / (2*n)) / (1 + z*z/n)
    radius = z * np.sqrt(k/n * (1-k/n)/n + z*z/(4*n*n)) / (1 + z*z/n)
    return float(center-radius), float(center+radius)


def evaluate(args):
    start = time.monotonic()
    args.out.mkdir(parents=True, exist_ok=True)
    cfg = json.loads((args.study / "manifest.json").read_text())["config"]
    status = json.loads((args.study / "simulation_status.json").read_text())
    if status["phase"] not in ("simulation-status", "audit-simulations") or status["state"] not in ("incomplete", "complete"):
        raise ValueError("First finish a validation-only archive pass")
    inventory = pd.read_csv(args.study / "simulation_inventory.csv")
    if len(inventory) != status["completed"] or status["failed"]:
        raise ValueError("Archive status/inventory mismatch")
    ref_names = ["scores.npz", "methods.json", "regions.csv", "provenance.json", "metrics.csv"]
    ref_hashes = verify_files(args.reference, ref_names)
    old_focal_hashes = verify_files(args.focal_reference, ["scores.npz", "methods.json", "regions.csv"])
    profile = json.loads(args.profiles.read_text())
    provenance = json.loads((args.reference / "provenance.json").read_text())
    if provenance["identity"]["profile_manifest"] != digest(args.profiles):
        raise ValueError("Reference profile manifest changed")
    for key, value in cfg.items():
        if key not in ("selection_coefficients", "selected_replicates_per_coefficient", "workers"):
            if profile["config"][key] != value:
                raise ValueError(f"Model differs from neutral reference: {key}")
    by_task = inventory.set_index("task_id")
    for item in profile["items"]:
        if item["onset_years"] not in (0, 50000):
            continue
        task = item["key"].replace("onset50000/", "selected_s0p005/")
        for name, column in (("simulation.trees", "simulation_sha256"),
                             ("decoded_input.trees", "decoder_input_sha256"),
                             ("focal_carriers.npy", "carriers_sha256")):
            if by_task.loc[task, column] != item["artifacts"][name]["sha256"]:
                raise ValueError(f"Reference uses a different simulation: {task}")
    selected = inventory[inventory["mode"] == "selected"].sort_values(["s", "replicate"])
    pairs = np.loadtxt(args.study / "pairs.tsv", dtype=int)
    payloads = [(row, cfg, pairs, args.out) for row in selected.to_dict("records")]
    results = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for result in pool.map(focal_one, payloads):
            results.append(result)
            if len(results) % 100 == 0:
                print(json.dumps(dict(phase="focal-truth", completed=len(results), target=len(payloads))), flush=True)
    focal = pd.DataFrame(results)
    focal.to_csv(args.out / "focal_truth.csv", index=False)
    old_focal_regions = pd.read_csv(args.focal_reference / "regions.csv")
    old_focal_methods = json.loads((args.focal_reference / "methods.json").read_text())
    with np.load(args.focal_reference / "scores.npz", allow_pickle=False) as saved:
        old_focal_scores = saved["scores"][0]
    parity = 0
    for row in focal[focal.s == .005].to_dict("records"):
        idx = old_focal_regions.index[old_focal_regions.key == f'onset50000/rep{row["replicate"]:04d}'][0]
        for j, method in enumerate(old_focal_methods):
            family = method["family"]
            if family not in ("af", "all", "mass", "aa"):
                continue
            value = row["sample_af"] if family == "af" else row[f'{family}_{method["years"]}']
            np.testing.assert_allclose(value, old_focal_scores[idx, j], rtol=0, atol=1e-14)
            parity += 1
    print(f"Verified {parity} focal scores against the earlier s=0.005 analysis", flush=True)
    regions = pd.read_csv(args.reference / "regions.csv")
    ref_methods = json.loads((args.reference / "methods.json").read_text())
    method_index = {m["name"]: j for j, m in enumerate(ref_methods)}
    with np.load(args.reference / "scores.npz", allow_pickle=False) as saved:
        # [fold, grid/site, truth/decoded, region, method]. Ungated scores have no fold dependence.
        scores = saved["scores"][0, 0, 0]
    neutral = regions.onset_years.to_numpy() == 0
    folds = regions.fold.to_numpy()
    replicate_folds = {int(r.key.split("rep")[-1]): int(r.fold)
                       for r in regions[~neutral].itertuples()}
    focal["fold"] = focal.replicate.map(replicate_folds)
    chosen = [dict(name="af_r1", family="af", years=0)]
    chosen += [dict(name=f"{family}_{'g0_' if family == 'mass' else ''}r1_T{cutoff}",
                    family=family, years=cutoff)
               for cutoff in cfg["tmrca_cutoffs_years"] for family in ("all", "mass")]
    predictions = []
    neutral_metrics = []
    old_metrics = pd.read_csv(args.reference / "metrics.csv")
    for method in chosen:
        ref = scores[:, method_index[method["name"]]]
        values = focal.sample_af.to_numpy() if method["family"] == "af" else focal[f'{method["family"]}_{method["years"]}'].to_numpy()
        p = np.zeros(len(focal))
        neutral_p = np.zeros(neutral.sum())
        for fold in range(5):
            cal = neutral & np.isin(folds, [(fold+3) % 5, (fold+4) % 5])
            assert cal.sum() == 400
            mask = focal.fold.to_numpy() == fold
            p[mask] = (1 + (ref[cal, None] >= values[mask]).sum(axis=0)) / 401
            heldout = folds[neutral] == fold
            neutral_p[heldout] = (1 + (ref[cal, None] >= ref[neutral][heldout]).sum(axis=0)) / 401
        for alpha in (.01, .05):
            false_calls = int((neutral_p <= alpha).sum())
            baseline = old_metrics[(old_metrics.scheme == "grid10kb") & (old_metrics.source == "truth")
                & (old_metrics.method == method["name"]) & (old_metrics.alpha == alpha)].iloc[0]
            assert false_calls == baseline.neutral_called
            neutral_metrics.append(dict(method=method["name"], alpha=alpha, neutral_regions=int(neutral.sum()),
                                        neutral_called=false_calls, neutral_region_call_fraction=false_calls/neutral.sum()))
            for k, row in focal.iterrows():
                predictions.append(dict(task_id=row.task_id, s=row.s, replicate=int(row.replicate),
                    fold=int(row.fold), method=method["name"], years=method["years"], alpha=alpha,
                    score=float(values[k]), p=float(p[k]), called=bool(p[k] <= alpha)))
    predictions = pd.DataFrame(predictions)
    predictions.to_csv(args.out / "selected_site_predictions.csv.gz", index=False,
                       compression=dict(method="gzip", mtime=0))
    pd.DataFrame(neutral_metrics).to_csv(args.out / "neutral_region_metrics.csv", index=False)
    summaries = []
    common = set.intersection(*(set(g.replicate) for _, g in focal.groupby("s")))
    for subset, data in (("all_saved", predictions), ("common_replicates", predictions[predictions.replicate.isin(common)])):
        for (s, method, years, alpha), group in data.groupby(["s", "method", "years", "alpha"]):
            k, n = int(group.called.sum()), len(group)
            lo, hi = wilson(k, n)
            summaries.append(dict(subset=subset, s=s, method=method, years=years, alpha=alpha,
                selected_regions=n, selected_site_called=k, regional_power_lower_bound=k/n,
                wilson_low=lo, wilson_high=hi))
    pd.DataFrame(summaries).to_csv(args.out / "selected_site_metrics.csv", index=False)
    af_summary = []
    for s, group in focal.groupby("s"):
        af = group.sample_af
        af_summary.append(dict(s=s, saved=len(group), target=100, missing=100-len(group),
            mean_af=af.mean(), median_af=af.median(), q10_af=af.quantile(.1), q90_af=af.quantile(.9),
            af_below_50pct=int((af < .5).sum()), fixed=int(group.fixed.sum()),
            median_aa_frac_recent_50000=group.aa_50000.median(),
            median_carrier_lineages_at_50k=group.carrier_lineages_at_50k.median(),
            single_carrier_lineage_at_50k=int((group.carrier_lineages_at_50k == 1).sum()),
            median_effective_carrier_lineages_at_50k=group.effective_carrier_lineages_at_50k.median(),
            spearman_af_mass_50000=float(spearmanr(af, group.mass_50000).statistic)))
    pd.DataFrame(af_summary).to_csv(args.out / "arm_summary.csv", index=False)
    comparisons = []
    primary = predictions[(predictions.alpha == .05) & predictions.method.isin(("af_r1", "mass_g0_r1_T50000"))]
    for s, group in primary.groupby("s"):
        wide = group.pivot(index="task_id", columns="method", values="called")
        af, mass = wide.af_r1, wide.mass_g0_r1_T50000
        ao, mo = int((af & ~mass).sum()), int((mass & ~af).sum())
        comparisons.append(dict(s=s, n=len(wide), both=int((af & mass).sum()), af_only=ao,
            mass_only=mo, neither=int((~af & ~mass).sum()),
            paired_exact_p=binomtest(mo, ao+mo).pvalue if ao+mo else 1.0))
    pd.DataFrame(comparisons).to_csv(args.out / "af_mass_paired_comparison.csv", index=False)
    atomic_json(args.out / "provenance.json", dict(study=str(args.study), config=cfg,
        inventory_sha256=digest(args.study / "simulation_inventory.csv"),
        sample_manifest_sha256=digest(args.study / "sample_manifest.json"),
        profile_manifest_sha256=digest(args.profiles), reference_hashes=ref_hashes,
        focal_reference_hashes=old_focal_hashes, source_sha256=digest(Path(__file__)),
        workers=args.workers, seconds=time.monotonic()-start, archive_status=status,
        endpoint="Selected causal-site truth versus neutral whole-region grid maxima; lower bound on full-scan power",
        lineage_definition="Distinct sampled-carrier ancestral branches just before 50-kya cutoff; not an identified count of introgressing founders",
        common_replicate_ids=sorted(common), neutral_calibration_per_fold=400,
        minimum_rank_p=1/401, comparisons_exploratory=True, new_simulations=False, decoding=False))
    atomic_json(args.out / "audit.json", dict(status="passed", selected_regions=len(focal),
        neutral_regions=int(neutral.sum()), saved_s005_focal_values_matched=parity,
        lineage_partition_vs_pairwise_tmrca_verified=True, neutral_call_rates_match_previous_analysis=True,
        original_reference_tree_hashes_matched=True, common_replicates_per_arm=len(common)))
    names = [p.name for p in args.out.iterdir() if p.is_file() and p.name != "artifact_manifest.json"]
    atomic_json(args.out / "artifact_manifest.json", {name: dict(bytes=(args.out/name).stat().st_size,
        sha256=digest(args.out/name)) for name in sorted(names)})
    print(pd.DataFrame(af_summary).to_string(index=False), flush=True)
    table = pd.DataFrame(summaries)
    print(table[(table.subset == "all_saved") & (table.alpha == .05) &
                table.method.isin(("af_r1", "mass_g0_r1_T50000", "all_r1_T50000"))].to_string(index=False), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("study", "reference", "focal-reference", "profiles", "out"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--workers", type=int, default=20)
    args = parser.parse_args()
    if not 1 <= args.workers <= 24:
        parser.error("Use 1 to 24 workers")
    evaluate(args)
