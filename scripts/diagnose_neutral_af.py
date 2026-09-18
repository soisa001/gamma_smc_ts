"""Separate pulse ancestry, surviving-marker frequency and regional maxima.

Reads previously saved ancestry replays only; never simulates or decodes.
"""
from concurrent.futures import ProcessPoolExecutor
import json
import os
from pathlib import Path

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "1"
import numpy as np
import pandas as pd
import tskit

from gamma_smc_aou.fresh_power import atomic_json, digest
from gamma_smc_aou.joint_scan_labels import ancestry_mask

ROOT = Path("/mnt/d/phase2simselection/sim")
OUT = ROOT / "eas_neutral_af_diagnostic_20260918"


def summary(x):
    x = np.asarray(x, dtype=float)
    return dict(n=len(x), mean=float(x.mean()), median=float(np.median(x)),
                q05=float(np.quantile(x, .05)), q95=float(np.quantile(x, .95)),
                maximum=float(x.max()), fraction_ge_018=float(np.mean(x >= .18)))


def fixed_coordinate(row):
    path = ROOT / "eas_joint_scan/labels" / row["task_id"]
    receipt = json.loads((path / "complete.json").read_text())
    assert receipt["seed"] == int(row["seed"])
    for name in ("ancestry_mapping.npz", "ancestry_replay.trees"):
        if digest(path / name) != receipt["outputs"][name]["sha256"]:
            raise ValueError(f"Corrupt saved ancestry: {path / name}")
    ts = tskit.load(path / "ancestry_replay.trees")
    with np.load(path / "ancestry_mapping.npz", allow_pickle=False) as saved:
        samples, ancestors = saved["samples"], saved["ancestors"]
    assert len(samples) == 400 and ts.sequence_length == 11_000_000
    # Prespecified coordinate in the original 11 Mb sequence, before any
    # movement to an observed variant. Zero-ancestry regions remain included.
    position = 5_500_000
    tree = ts.at(position)
    lookup = np.full(ts.num_nodes, -1, dtype=int)
    lookup[samples] = np.arange(len(samples))
    mask = ancestry_mask(tree, ancestors, lookup, position)
    groups = []
    for sample in samples:
        node = int(sample)
        while tree.parent(node) != tskit.NULL and tree.time(tree.parent(node)) < 2000:
            node = tree.parent(node)
        groups.append(node)
    groups = np.asarray(groups)
    sizes = np.unique(groups[mask], return_counts=True)[1]
    return dict(task_id=row["task_id"], ancestry_af=float(mask.mean()),
        ancestral_lineages_at_pulse=len(np.unique(groups)),
        carrier_lineages_at_pulse=len(sizes), carriers=int(mask.sum()),
        carrier_frac_recent_50000=(float(np.sum(sizes*(sizes-1))/(mask.sum()*(mask.sum()-1)))
                                  if mask.sum() > 1 else None),
        replay_sha256=receipt["outputs"]["ancestry_replay.trees"]["sha256"],
        mapping_sha256=receipt["outputs"]["ancestry_mapping.npz"]["sha256"])


def verify_reference(directory, names):
    manifest = json.loads((directory / "artifact_manifest.json").read_text())
    for name in names:
        if digest(directory / name) != manifest[name]["sha256"]:
            raise ValueError(f"Corrupt reference: {directory / name}")


def main():
    OUT.mkdir(exist_ok=True)
    ref = ROOT / "eas_allele_class_ablation"
    verify_reference(ref, ("regions.csv", "scores.npz", "methods.json"))
    regions = pd.read_csv(ref / "regions.csv")
    methods = json.loads((ref / "methods.json").read_text())
    neutral = regions.onset_years.to_numpy() == 0
    with np.load(ref / "scores.npz", allow_pickle=False) as saved:
        scores = saved["scores"][0, 0, 0]
    index = {m["name"]: j for j, m in enumerate(methods)}
    result = dict(regional_neutral_maxima={}, fold_thresholds={}, nominal_pulse=.02)
    for method in ("af_r1", "mass_g0_r1_T50000"):
        values = scores[:, index[method]]
        result["regional_neutral_maxima"][method] = summary(values[neutral])
        boundaries = []
        for fold in range(5):
            calibration = neutral & np.isin(regions.fold, [(fold+3) % 5, (fold+4) % 5])
            assert calibration.sum() == 400
            boundaries.append(float(np.sort(values[calibration])[-20]))
        result["fold_thresholds"][method] = dict(alpha=.05, strict_greater_than=boundaries)
    focal_dir = ROOT / "eas_allele_focal_comparison"
    verify_reference(focal_dir, ("regions.csv",))
    focal = pd.read_csv(focal_dir / "regions.csv")
    focal = focal[focal.onset_years == 0]
    result["nearest_observed_archaic_marker"] = summary(focal.loc[focal.marker_position.notna(), "af"])
    result["neutral_regions_without_markers"] = int(focal.marker_position.isna().sum())
    inventory_path = ROOT / "eas_q02_h400/simulation_inventory.csv"
    inventory = pd.read_csv(inventory_path)
    rows = inventory[inventory["mode"] == "neutral"].to_dict("records")
    assert len(rows) == 1000
    with ProcessPoolExecutor(max_workers=20) as pool:
        data = []
        for row in pool.map(fixed_coordinate, rows):
            data.append(row)
            if len(data) % 100 == 0:
                print(f"Read {len(data)}/1000 saved neutral ancestry replays", flush=True)
    data = pd.DataFrame(data)
    data.to_csv(OUT / "neutral_fixed_coordinate.csv", index=False)
    result["fixed_coordinate_all_regions"] = summary(data.ancestry_af)
    present = data.carriers > 0
    result["fixed_coordinate_conditional_on_ancestry_present"] = summary(data.loc[present, "ancestry_af"])
    result["probability_ancestry_present_at_fixed_coordinate"] = float(present.mean())
    result["conditional_single_carrier_lineage_fraction"] = float((data.loc[present, "carrier_lineages_at_pulse"] == 1).mean())
    result["mean_lineages_at_pulse"] = float(data.ancestral_lineages_at_pulse.mean())
    # The conditional mean identity is exact for the empirical distribution.
    np.testing.assert_allclose(data.ancestry_af.mean(), present.mean()*data.loc[present,"ancestry_af"].mean(), rtol=1e-14)
    selected_path = ROOT / "eas_h400_pause_evaluation_20260917/focal_truth.csv"
    selected = pd.read_csv(selected_path)
    result["selected_s001_survival_conditioned_af"] = summary(selected.loc[selected.s == .001, "sample_af"])
    result["provenance"] = dict(source_sha256=digest(Path(__file__)),
        inventory_sha256=digest(inventory_path), selected_focal_sha256=digest(selected_path),
        reference_artifact_manifest_sha256=digest(ref / "artifact_manifest.json"),
        focal_reference_artifact_manifest_sha256=digest(focal_dir / "artifact_manifest.json"),
        workers=20, coordinate_raw_bp=5_500_000, haplotypes=400, pulse_years=50000,
        generation_time_years=25, simulations_started=False, decoding_started=False,
        fixed_coordinate_endpoint="Local archaic ancestry fraction, including zero-ancestry regions; not a marker-ascertained AF")
    atomic_json(OUT / "summary.json", result)
    atomic_json(OUT / "artifact_manifest.json", {name:dict(sha256=digest(OUT/name),bytes=(OUT/name).stat().st_size)
                for name in ("summary.json", "neutral_fixed_coordinate.csv")})
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
