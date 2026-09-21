"""Pointwise power/FPR along saved s=0.005 profiles; no spatial maximum."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "1"
import numpy as np
import pandas as pd

ROOT = Path("/mnt/d/phase2simselection/sim")
OUT = ROOT / "eas_positional_distance_h400"
TIMES = (5000, 10000, 20000, 30000, 40000, 50000)
METHODS = ["af"]+[f"{family}_{t}" for t in TIMES for family in ("all", "mass")]


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def read(item):
    directory = Path(item.get("profile_root", ROOT / "eas_joint_scan/profiles")) / item["key"]
    receipt = json.loads((directory / "complete.json").read_text())
    sha = digest(directory / "features.npz")
    if sha != receipt["outputs"]["features.npz"]["sha256"]:
        raise ValueError("Corrupt saved spatial profile")
    with np.load(directory / "features.npz", allow_pickle=False) as features:
        if "archaic_labels_scope" in features and str(features["archaic_labels_scope"]) == "focal_only":
            if not item.get("all_pairs_only", False):
                raise ValueError("Spatial carrier inference requires labels across the region; only focal labels are available")
        np.testing.assert_array_equal(features["grid"], np.arange(1000)*10000)
        np.testing.assert_array_equal(features["cutoffs"], TIMES)
        markers = features["grid_markers"]
        valid = markers >= 0
        counts, n = features["grid_counts"], features["grid_n"][:, 3]
        assert np.all(n == 10000)
        values = np.zeros((2, 1000, len(METHODS)))
        values[:, valid, 0] = features["site_af"][markers[valid]]
        for k, _ in enumerate(TIMES):
            values[:, :, 1+2*k] = counts[:, :, 3, k]/n
            values[:, :, 2+2*k] = counts[:, :, 2, k]/n
        if item.get("all_pairs_only", False):
            values = values[:, :, 1::2]
        if not item.get("all_pairs_only", False):
            assert not values[:, ~valid, 2::2].any()
    return values, sha


def main(out=OUT, selected_s=.005, selected_profiles=None, focal_dir=ROOT / "eas_positional_h400", all_pairs_only=False):
    OUT = out
    methods = [f"all_{t}" for t in TIMES] if all_pairs_only else METHODS
    OUT.mkdir(parents=True, exist_ok=True)
    fold_path = ROOT / "eas_allele_class_ablation/regions.csv"
    manifest = json.loads((fold_path.parent / "artifact_manifest.json").read_text())
    assert digest(fold_path) == manifest[fold_path.name]["sha256"]
    regions = pd.read_csv(fold_path)
    regions["all_pairs_only"] = all_pairs_only
    if selected_profiles is not None:
        regions["profile_root"] = str(ROOT / "eas_joint_scan/profiles")
        regions.loc[regions.key.str.startswith("onset50000/"), "profile_root"] = str(selected_profiles)
    with ThreadPoolExecutor(max_workers=2) as pool:
        loaded = list(pool.map(read, regions.to_dict("records")))
    values = np.array([entry[0] for entry in loaded])
    neutral = regions.key.str.startswith("neutral/").to_numpy()
    folds = regions.fold.to_numpy()
    assert neutral.sum() == 1000 and (~neutral).sum() == 100
    pvalues = np.zeros_like(values)
    for fold in range(5):
        calibration = neutral & np.isin(folds, [(fold+3)%5, (fold+4)%5])
        test = np.flatnonzero(folds == fold)
        assert calibration.sum() == 400
        reference = np.sort(values[calibration], axis=0)
        for source in range(2):
            for method in range(len(methods)):
                for position in range(1000):
                    scores = values[test, source, position, method]
                    below = np.searchsorted(reference[:, source, position, method], scores, side="left")
                    pvalues[test, source, position, method] = (401-below)/401
    metrics = []
    for source, source_name in enumerate(("truth", "decoded")):
        for method, method_name in enumerate(methods):
            for alpha in (.01, .05):
                calls = pvalues[:, source, :, method] <= alpha
                positives = calls[~neutral].sum(axis=0)
                negatives = calls[neutral].sum(axis=0)
                for j in range(1000):
                    metrics.append(dict(source=source_name, method=method_name, alpha=alpha,
                        position_bp=j*10000, distance_bp=j*10000-5_000_000,
                        selected_positions=100, selected_called=int(positives[j]), power=positives[j]/100,
                        neutral_positions=1000, neutral_called=int(negatives[j]), positional_fpr=negatives[j]/1000))
    metrics = pd.DataFrame(metrics)
    metrics.to_csv(OUT / "metrics.csv", index=False)
    chosen = (0, 10000, 20000, 50000, 100000, 200000, 500000, 1000000, 2000000)
    metrics[metrics.distance_bp.abs().isin(chosen)].to_csv(OUT / "selected_distances.csv", index=False)
    # The center must independently reproduce the completed focal evaluation.
    focal = pd.read_csv(focal_dir / "metrics.csv")
    focal = focal[focal.s == selected_s]
    merged = metrics[metrics.distance_bp == 0].merge(focal, on=["source", "method", "alpha"], suffixes=("_spatial", "_focal"))
    assert len(merged) == 4*len(methods)
    for column in ("selected_called", "neutral_called"):
        np.testing.assert_array_equal(merged[column+"_spatial"], merged[column+"_focal"])
    audit = dict(status="passed", grid_positions=1000, neutral_regions=1000, selected_regions=100,
        selection_coefficient=selected_s, matched_coordinate_in_neutral_calibration=True,
        all_pairs_only=all_pairs_only,
        no_spatial_or_time_maximum=True, focal_counts_match_independent_evaluation=True,
        interpretation="Pointwise rejection probability at each distance, not gene-level call probability or FDR.")
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2)+"\n")
    provenance = dict(source_sha256=digest(Path(__file__)), fold_reference_sha256=digest(fold_path),
        feature_hashes={r.key: item[1] for r, item in zip(regions.itertuples(), loaded)},
        stride_bp=10000, marker_radius_bp=5000, cutoffs_years=TIMES, calibration_per_fold=400,
        simulations_started=False, decoding_started=False)
    (OUT / "provenance.json").write_text(json.dumps(provenance, indent=2)+"\n")
    names = ("metrics.csv", "selected_distances.csv", "audit.json", "provenance.json")
    (OUT / "artifact_manifest.json").write_text(json.dumps(
        {name: dict(sha256=digest(OUT/name), bytes=(OUT/name).stat().st_size) for name in names}, indent=2)+"\n")
    print(json.dumps(audit, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--selected-s", type=float, default=.005)
    parser.add_argument("--selected-profiles", type=Path)
    parser.add_argument("--focal-dir", type=Path, default=ROOT / "eas_positional_h400")
    parser.add_argument("--all-pairs-only", action="store_true")
    args = parser.parse_args()
    main(args.out, args.selected_s, args.selected_profiles, args.focal_dir, args.all_pairs_only)
