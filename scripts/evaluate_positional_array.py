"""Matched-position AF and carrier-TMRCA calibration, without regional maxima."""
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "1"
import numpy as np
import pandas as pd

from gamma_smc_aou.fresh_power import atomic_json, digest

ROOT = Path("/mnt/d/phase2simselection/sim")
OUT = ROOT / "eas_positional_h400"
CUTOFFS = (5000, 10000, 20000, 30000, 40000, 50000)
METHODS = ["af"] + [f"{family}_{t}" for t in CUTOFFS for family in ("all", "mass")]


def read_point(item):
    directory = ROOT / "eas_joint_scan/profiles" / item["key"]
    receipt = json.loads((directory / "complete.json").read_text())
    if digest(directory / "features.npz") != receipt["outputs"]["features.npz"]["sha256"]:
        raise ValueError(f"Corrupt saved profile: {directory}")
    with np.load(directory / "features.npz", allow_pickle=False) as features:
        index = 500  # 5 Mb on the existing 10 kb grid.
        marker = int(features["grid_markers"][index])
        af = float(features["site_af"][marker]) if marker >= 0 else 0.
        counts, denominators = features["grid_counts"][:, index], features["grid_n"][index]
        assert denominators[3] == 10000
        values = np.zeros((2, len(METHODS)))
        values[:, 0] = af
        for j, t in enumerate(CUTOFFS):
            values[:, 1+2*j] = counts[:, 3, j]/denominators[3]
            values[:, 2+2*j] = counts[:, 2, j]/denominators[3]
        if marker < 0:
            assert af == 0 and not counts[:, 2].any()
        coordinate = float(features["sites"][marker]) if marker >= 0 else None
        if coordinate is not None:
            assert abs(coordinate-5_000_000) <= 5000
    return dict(key=item["key"], fold=item["fold"], marker_present=marker >= 0,
                marker_position=coordinate, values=values, profile_sha256=digest(directory / "features.npz"))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    old = ROOT / "eas_allele_class_ablation"
    artifact = json.loads((old / "artifact_manifest.json").read_text())
    if digest(old / "regions.csv") != artifact["regions.csv"]["sha256"]:
        raise ValueError("Corrupt fold assignments")
    regions = pd.read_csv(old / "regions.csv")
    # Read-only feature extraction; simulation generation remains paused.
    with ThreadPoolExecutor(max_workers=4) as pool:
        points = list(pool.map(read_point, regions.to_dict("records")))
    by_key = {p["key"]: p for p in points}
    neutral_items = [p for p in points if p["key"].startswith("neutral/")]
    neutral_values = np.array([p["values"] for p in neutral_items])
    neutral_folds = np.array([p["fold"] for p in neutral_items])
    focal_dir = ROOT / "eas_h400_pause_evaluation_20260917"
    focal_artifacts = json.loads((focal_dir / "artifact_manifest.json").read_text())
    if digest(focal_dir / "focal_truth.csv") != focal_artifacts["focal_truth.csv"]["sha256"]:
        raise ValueError("Corrupt selected focal truth")
    selected = pd.read_csv(focal_dir / "focal_truth.csv")
    values = np.zeros((len(selected), len(METHODS)))
    values[:, 0] = selected.sample_af
    for k, method in enumerate(METHODS[1:], 1):
        values[:, k] = selected[method]
    selected_folds = np.array([by_key[f"onset50000/rep{int(r):04d}"]["fold"] for r in selected.replicate])
    for j, row in selected.iterrows():
        if row.s == .005:
            np.testing.assert_allclose(values[j], by_key[f"onset50000/rep{int(row.replicate):04d}"]["values"][0], atol=1e-14, rtol=0)
    predictions = []
    for src, source in enumerate(("truth", "decoded")):
        source_selected = np.flatnonzero(np.ones(len(selected), dtype=bool) if source == "truth" else selected.s == .005)
        for fold in range(5):
            cal = np.isin(neutral_folds, [(fold+3)%5, (fold+4)%5])
            test = neutral_folds == fold
            assert cal.sum() == 400 and test.sum() == 200
            for k, method in enumerate(METHODS):
                reference = neutral_values[cal, src, k]
                for j in np.flatnonzero(test):
                    score = neutral_values[j, src, k]
                    p = float((1+(reference >= score).sum())/401)
                    for alpha in (.01, .05):
                        predictions.append(dict(task_id=neutral_items[j]["key"], mode="neutral", s=0., source=source,
                            method=method, fold=fold, alpha=alpha, score=score, p=p, called=bool(p <= alpha)))
                for j in source_selected[selected_folds[source_selected] == fold]:
                    row = selected.iloc[j]
                    score = values[j, k] if source == "truth" else by_key[f"onset50000/rep{int(row.replicate):04d}"]["values"][1, k]
                    p = float((1+(reference >= score).sum())/401)
                    for alpha in (.01, .05):
                        predictions.append(dict(task_id=row.task_id, mode="selected", s=row.s, source=source,
                            method=method, fold=fold, alpha=alpha, score=score, p=p, called=bool(p <= alpha)))
    predictions = pd.DataFrame(predictions)
    predictions.to_csv(OUT / "predictions.csv.gz", index=False, compression=dict(method="gzip", mtime=0))
    metrics = []
    for (source, method, alpha, s), group in predictions[predictions["mode"] == "selected"].groupby(["source", "method", "alpha", "s"]):
        background = predictions[(predictions["mode"] == "neutral") & (predictions.source == source)
            & (predictions.method == method) & (predictions.alpha == alpha)]
        assert len(background) == 1000
        metrics.append(dict(source=source, method=method, alpha=alpha, s=s, selected_regions=len(group),
            selected_called=int(group.called.sum()), power=float(group.called.mean()), neutral_positions=1000,
            neutral_called=int(background.called.sum()), positional_fpr=float(background.called.mean())))
    metrics = pd.DataFrame(metrics)
    metrics.to_csv(OUT / "metrics.csv", index=False)
    pd.DataFrame([dict(key=p["key"], fold=p["fold"], marker_present=p["marker_present"],
        marker_position=p["marker_position"], af=p["values"][0, 0], profile_sha256=p["profile_sha256"])
        for p in neutral_items]).to_csv(OUT / "neutral_positions.csv", index=False)
    np.savez_compressed(OUT / "scores.npz", neutral=neutral_values, neutral_folds=neutral_folds,
                        selected_truth=values, selected_folds=selected_folds)
    # Check saved prediction ranks against an independent matrix calculation.
    for (source, method, fold), group in predictions.groupby(["source", "method", "fold"]):
        reference = neutral_values[np.isin(neutral_folds, [(fold+3)%5, (fold+4)%5]), ("truth", "decoded").index(source), METHODS.index(method)]
        expected = (1+(reference[:, None] >= group.score.to_numpy()).sum(axis=0))/401
        np.testing.assert_allclose(expected, group.p, atol=1e-14, rtol=0)
        np.testing.assert_array_equal(expected <= group.alpha, group.called)
    atomic_json(OUT / "audit.json", dict(status="passed", neutral_positions=1000,
        neutral_positions_with_archaic_marker=sum(p["marker_present"] for p in neutral_items),
        truth_selected_positions=len(selected), decoded_selected_positions=100,
        s005_truth_matches_existing_profiles=True, pvalues_recomputed=True, regional_maximum_used=False))
    atomic_json(OUT / "provenance.json", dict(source_sha256=digest(Path(__file__)), focal_truth_sha256=digest(focal_dir / "focal_truth.csv"),
        fold_reference_sha256=digest(old / "regions.csv"), marker_radius_bp=5000, tmrca_position_bp=5_000_000,
        cutoffs_years=CUTOFFS, methods=METHODS, calibration_regions_per_fold=400, holdout_regions_per_fold=200,
        statistic="frac_recent_T; carrier mass = n_recent_ALT_ALT/10000", null="same fixed coordinate in every neutral region; absent nearby archaic markers score zero",
        input_profile_hashes={p["key"]:p["profile_sha256"] for p in points},
        new_simulations=False, new_decoding=False, scope="pointwise at a pre-specified coordinate, not family-wise over a genomic scan"))
    names = ("predictions.csv.gz", "metrics.csv", "neutral_positions.csv", "scores.npz", "audit.json", "provenance.json")
    atomic_json(OUT / "artifact_manifest.json", {name: dict(bytes=(OUT/name).stat().st_size, sha256=digest(OUT/name)) for name in names})
    print(metrics[(metrics.alpha == .05) & metrics.method.isin(("af", "mass_50000", "all_50000"))].to_string(index=False))


if __name__ == "__main__":
    main()
