"""Checked local-null thresholds, paired AF comparison, and pooled p=0.001 ranks."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest


def digest(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def main(directory):
    manifest = json.loads((directory / "artifact_manifest.json").read_text())
    for name in ("predictions.csv.gz", "scores.npz", "neutral_positions.csv", "provenance.json"):
        if digest(directory / name) != manifest[name]["sha256"]:
            raise ValueError(f"Corrupt input: {name}")
    predictions = pd.read_csv(directory / "predictions.csv.gz", float_precision="round_trip")
    positions = pd.read_csv(directory / "neutral_positions.csv", float_precision="round_trip")
    methods = json.loads((directory / "provenance.json").read_text())["methods"]
    with np.load(directory / "scores.npz", allow_pickle=False) as saved:
        neutral = saved["neutral"]
        folds = saved["neutral_folds"]
    thresholds, comparisons, pooled = [], [], []
    for (source, method), group in predictions[predictions.alpha == .05].groupby(["source", "method"]):
        reference = neutral[:, ("truth", "decoded").index(source), methods.index(method)]
        for fold in range(5):
            values = np.sort(reference[np.isin(folds, [(fold+3)%5, (fold+4)%5])])
            for alpha in (.01, .05):
                # With 400 calibration values, p <= alpha iff score is strictly
                # above this order statistic. Ties are deliberately conservative.
                max_exceedances = int(np.floor(alpha*(len(values)+1)-1))
                cutoff = values[-max_exceedances-1]
                tested = predictions[(predictions.source == source) & (predictions.method == method)
                    & (predictions.fold == fold) & (predictions.alpha == alpha)]
                np.testing.assert_array_equal(tested.score > cutoff, tested.called)
                thresholds.append(dict(source=source, method=method, fold=fold, alpha=alpha,
                    calibration_positions=len(values), strict_score_threshold=cutoff))
        for row in group.itertuples(index=False):
            is_neutral = row.mode == "neutral"
            # Each neutral is compared with the other 999 neutrals. Selected
            # positions use all 1,000. No fitted normalization is needed here.
            exceedances = int(np.count_nonzero(reference >= row.score))-int(is_neutral)
            n = len(reference)-int(is_neutral)
            p = (1+exceedances)/(n+1)
            pooled.append(dict(task_id=row.task_id, mode=row.mode, source=source, method=method,
                s=row.s, score=row.score, calibration_positions=n, p=p, alpha=.001,
                called=bool(p <= .001)))
    selected = predictions[predictions["mode"] == "selected"]
    for (source, alpha, s), group in selected.groupby(["source", "alpha", "s"]):
        af = group[group.method == "af"].set_index("task_id").called
        for method, mass in group[group.method.str.startswith("mass_")].groupby("method"):
            calls = mass.set_index("task_id").called.reindex(af.index)
            assert not calls.isna().any()
            gain = int((calls & ~af).sum())
            loss = int((~calls & af).sum())
            comparisons.append(dict(source=source, alpha=alpha, s=s, method=method, n=len(af),
                af_called=int(af.sum()), mass_called=int(calls.sum()), mass_only=gain, af_only=loss,
                power_difference=float(calls.mean()-af.mean()),
                paired_exact_p=binomtest(gain, gain+loss, .5).pvalue if gain+loss else 1.,
                interpretation="exploratory; unadjusted across time cutoffs and selection arms"))
    pd.DataFrame(thresholds).to_csv(directory / "calibration_thresholds.csv", index=False)
    pd.DataFrame(comparisons).to_csv(directory / "paired_af_comparison.csv", index=False)
    pooled = pd.DataFrame(pooled)
    pooled.to_csv(directory / "pooled_rank_predictions.csv.gz", index=False,
                  compression=dict(method="gzip", mtime=0))
    metrics = []
    for (source, method, s), group in pooled[pooled["mode"] == "selected"].groupby(["source", "method", "s"]):
        background = pooled[(pooled["mode"] == "neutral") & (pooled.source == source) & (pooled.method == method)]
        assert len(background) == 1000
        metrics.append(dict(source=source, method=method, alpha=.001, s=s,
            selected_positions=len(group), selected_called=int(group.called.sum()), power=float(group.called.mean()),
            neutral_positions=len(background), neutral_called=int(background.called.sum()),
            positional_fpr=float(background.called.mean()), selected_calibration_positions=1000,
            neutral_calibration_positions=999, calibration="pooled neutral ranks; leave-one-out for neutral evaluation"))
    pd.DataFrame(metrics).to_csv(directory / "pooled_rank_metrics.csv", index=False)
    af = positions.af.to_numpy()
    summary = dict(status="passed", source_sha256=digest(Path(__file__)),
        source_manifest_sha256=digest(directory / "artifact_manifest.json"),
        neutral_positions=1000, neutral_positions_with_marker=int(positions.marker_present.sum()),
        neutral_positions_af_at_least_018=int(np.count_nonzero(af >= .18)),
        mean_positional_marker_af=float(af.mean()),
        pooled_p001_caveat="Only one extreme neutral rank is available; this is not precise validation of a 0.1% tail.",
        comparisons="Paired tests are exploratory and not adjusted for multiple comparisons; T was not tuned per replicate.")
    (directory / "supplement_audit.json").write_text(json.dumps(summary, indent=2)+"\n")
    names = ("calibration_thresholds.csv", "paired_af_comparison.csv", "pooled_rank_predictions.csv.gz",
             "pooled_rank_metrics.csv", "supplement_audit.json")
    (directory / "supplement_manifest.json").write_text(json.dumps(
        {n: dict(sha256=digest(directory/n), bytes=(directory/n).stat().st_size) for n in names}, indent=2)+"\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=Path("/mnt/d/phase2simselection/sim/eas_positional_h400"))
    main(parser.parse_args().directory)
