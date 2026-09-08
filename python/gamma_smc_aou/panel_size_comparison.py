"""Compare fresh 200/400-haplotype panels at a fixed number of pairs.

Every region is one evaluation unit. Folds split whole replicates; neutral
normalization and null calibration use disjoint training subsets. All selected
replicates, including sample fixation, remain in training and evaluation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path

os.environ["MPLBACKEND"] = "Agg"
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "1"

import numpy as np
import pandas as pd

from .fresh_region_metrics import load_profiles
from .partial_sweep_cv import make_folds


def empirical_upper_p(reference, values):
    ordered = np.sort(reference)
    return (1 + len(ordered) - np.searchsorted(ordered, values, side="left")) / (
        len(ordered) + 1
    )


def target_threshold(training_selected, target):
    if not 0 < target <= 1 or not len(training_selected):
        raise ValueError("A sensitivity target and nonempty training set are required")
    return float(
        np.sort(training_selected)[-math.ceil(target * len(training_selected))]
    )


def fold_roles(meta, folds, fold):
    test = np.flatnonzero(folds == fold)
    neutral_train = np.flatnonzero((folds != fold) & (meta["mode"] == "neutral"))
    selected_train = np.flatnonzero((folds != fold) & (meta["mode"] == "selected"))
    # Each reference subset excludes every test region. Alternation is fixed
    # before reading statistic values and used identically for both panels.
    fit, calibrate = neutral_train[::2], neutral_train[1::2]
    return test, fit, calibrate, selected_train


def evaluate(values, meta, cutoffs, seed, haplotypes, source):
    folds = make_folds(meta, seed, 5)
    rows, choices = [], []
    for fold in range(5):
        test, fit, calibrate, selected_train = fold_roles(meta, folds, fold)
        mean = values[fit].mean(axis=0)
        sd = np.maximum(values[fit].std(axis=0), 1e-4)
        for method in ("maximum_frac_recent", "maximum_standardized_frac_recent"):
            scores = (
                values if method == "maximum_frac_recent" else (values - mean) / sd
            ).max(axis=1)
            for t, cutoff in enumerate(cutoffs):
                pvalues = empirical_upper_p(scores[calibrate, t], scores[test, t])
                threshold = target_threshold(scores[selected_train, t], 0.7)
                choices.append(
                    dict(
                        haplotypes=haplotypes,
                        source=source,
                        method=method,
                        cutoff_years=cutoff,
                        fold=fold,
                        selected_target=0.7,
                        threshold=threshold,
                        fit_neutral_n=len(fit),
                        calibration_neutral_n=len(calibrate),
                        training_selected_n=len(selected_train),
                    )
                )
                for j, i in enumerate(test):
                    m = meta.iloc[i]
                    base = dict(
                        haplotypes=haplotypes,
                        source=source,
                        method=method,
                        cutoff_years=cutoff,
                        task_id=m.task_id,
                        fold=fold,
                        mode=m["mode"],
                        sample_af=float(m.sample_af),
                        fixed=bool(m.fixed),
                    )
                    for alpha in (0.01, 0.05, 0.1):
                        rows.append(
                            dict(
                                base,
                                policy="neutral_target",
                                target=alpha,
                                called=int(pvalues[j] <= alpha),
                                p=float(pvalues[j]),
                            )
                        )
                    rows.append(
                        dict(
                            base,
                            policy="selected_target",
                            target=0.7,
                            called=int(scores[i, t] >= threshold),
                            p=float(pvalues[j]),
                        )
                    )
    return rows, choices


def figures(summary, diagnostics, out):
    import matplotlib
    import matplotlib.pyplot as plt

    matplotlib.rcParams.update(
        {"font.size": 17, "axes.labelsize": 19, "pdf.fonttype": 42, "ps.fonttype": 42}
    )
    colors = {200: "#316aab", 400: "#ca672c"}
    fig, axes = plt.subplots(1, 2, figsize=(11, 8.5))
    fig.subplots_adjust(left=0.1, right=0.98, bottom=0.14, top=0.78, wspace=0.30)
    for haplotypes in sorted(diagnostics.haplotypes.unique()):
        for source, linestyle in [("decoded", "-"), ("truth", "--")]:
            block = diagnostics[
                (diagnostics.haplotypes == haplotypes) & (diagnostics.source == source)
            ]
            neutral = block[block["mode"] == "neutral"].set_index("cutoff_years")
            selected = block[block["mode"] == "selected"].set_index("cutoff_years")
            x = neutral.index.to_numpy() / 1000
            label = f"{haplotypes} haplotypes, {source}"
            axes[0].plot(
                x,
                neutral.focal_sd,
                color=colors[haplotypes],
                linestyle=linestyle,
                linewidth=2.5,
                label=label,
            )
            axes[1].plot(
                x,
                selected.focal_mean - neutral.focal_mean,
                color=colors[haplotypes],
                linestyle=linestyle,
                linewidth=2.5,
            )
    axes[0].set_ylabel("Neutral SD of frac_recent_T")
    axes[1].set_ylabel("Selected minus neutral mean")
    for ax in axes:
        ax.set_xlabel("TMRCA cutoff (kya)")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_title("Variability at the focal position")
    axes[1].set_title("Signal at the focal position")
    fig.legend(
        *axes[0].get_legend_handles_labels(),
        loc="upper center",
        ncol=2,
        bbox_to_anchor=(0.5, 0.94),
        frameon=False,
        fontsize=15,
    )
    fig.suptitle("Panel size at s = 0.005; 10,000 pairs", y=0.99, fontsize=21)
    fig.savefig(out / "panel_noise.png", dpi=200)
    fig.savefig(out / "panel_noise.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 8.5), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.11, right=0.98, bottom=0.15, top=0.80, wspace=0.10)
    for ax, source in zip(axes, ("decoded", "truth")):
        for haplotypes in sorted(summary.haplotypes.unique()):
            block = summary[
                (summary.haplotypes == haplotypes) & (summary.source == source)
            ]
            key = ["method", "cutoff_years", "policy", "target"]
            selected = block[block["mode"] == "selected"].set_index(key)
            neutral = block[block["mode"] == "neutral"].set_index(key)
            ax.scatter(
                selected.region_call_rate,
                neutral.loc[selected.index].region_call_rate,
                color=colors[haplotypes],
                s=44,
                alpha=0.55,
                label=f"{haplotypes} haplotypes",
            )
        ax.axvline(0.7, color="gray", linestyle="--", linewidth=1.3)
        ax.axhline(0.05, color="gray", linestyle=":", linewidth=1.3)
        ax.set(xlim=(-0.01, 1.01), ylim=(-0.01, 1.01), title=source.capitalize())
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Fraction of neutral regions called")
    fig.supxlabel("Fraction of selected regions called", y=0.04, fontsize=20)
    fig.legend(
        *axes[0].get_legend_handles_labels(),
        loc="upper center",
        ncol=2,
        bbox_to_anchor=(0.5, 0.925),
        frameon=False,
    )
    fig.suptitle("Held-out calls anywhere in a 10 Mb region", y=0.99, fontsize=21)
    fig.savefig(out / "panel_region_calls.png", dpi=200)
    fig.savefig(out / "panel_region_calls.pdf")
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--larger", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=20)
    args = parser.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    outputs, diagnostics, errors, choices, provenance = [], [], [], [], []
    previous_cfg = None
    for root in (args.baseline.resolve(), args.larger.resolve()):
        paired = {}
        for source, filename in (
            ("decoded", "frac_recent.tsv"),
            ("truth", "truth_frac_recent.tsv"),
        ):
            cfg, fingerprint, positions, meta, counts, inputs, manifest_hash = (
                load_profiles(root, args.workers, filename)
            )
            if previous_cfg is not None:
                allowed = {
                    "sample_diploids",
                    "selection_coefficients",
                    "workers",
                    "memory_limit_gb",
                    "storage_limit_bytes",
                }
                assert {k: v for k, v in cfg.items() if k not in allowed} == {
                    k: v for k, v in previous_cfg.items() if k not in allowed
                }, "Panel comparison changes parameters beyond sample size"
            previous_cfg = cfg
            use = (meta["mode"] == "neutral") | (
                (meta["mode"] == "selected") & (meta.s == 0.005)
            )
            meta = meta[use].reset_index(drop=True)
            values = counts[use].astype(np.float64) / cfg["haplotype_pairs"]
            paired[source] = values
            haplotypes = 2 * cfg["sample_diploids"]
            focal = int(np.flatnonzero(positions == cfg["focal_position_bp"])[0])
            for mode in ("neutral", "selected"):
                v = values[meta["mode"] == mode]
                for t, cutoff in enumerate(cfg["tmrca_cutoffs_years"]):
                    diagnostics.append(
                        dict(
                            haplotypes=haplotypes,
                            source=source,
                            mode=mode,
                            cutoff_years=cutoff,
                            n=len(v),
                            focal_mean=float(v[:, focal, t].mean()),
                            focal_sd=float(v[:, focal, t].std(ddof=1)),
                            region_mean=float(v[:, :, t].mean()),
                            region_mean_sd=float(v[:, :, t].mean(axis=1).std(ddof=1)),
                            average_position_sd=float(
                                v[:, :, t].std(axis=0, ddof=1).mean()
                            ),
                            profile_roughness=float(
                                np.abs(np.diff(v[:, :, t], axis=1)).mean()
                            ),
                        )
                    )
            result, fitted = evaluate(
                values,
                meta,
                cfg["tmrca_cutoffs_years"],
                cfg["seed_base"],
                haplotypes,
                source,
            )
            outputs.extend(result)
            choices.extend(fitted)
            provenance.append(
                dict(
                    root=str(root),
                    haplotypes=haplotypes,
                    source=source,
                    fingerprint=fingerprint,
                    config=cfg,
                    sample_manifest_sha256=manifest_hash,
                    profile_hashes=inputs,
                )
            )
            print(f"Scored {haplotypes} haplotypes, {source}", flush=True)
        difference = paired["decoded"] - paired["truth"]
        for mode in ("neutral", "selected"):
            v = difference[meta["mode"] == mode]
            for t, cutoff in enumerate(cfg["tmrca_cutoffs_years"]):
                errors.append(
                    dict(
                        haplotypes=haplotypes,
                        mode=mode,
                        cutoff_years=cutoff,
                        n=len(v),
                        mean_region_rmse=float(
                            np.sqrt((v[:, :, t] ** 2).mean(axis=1)).mean()
                        ),
                        mean_bias=float(v[:, :, t].mean()),
                        error_roughness=float(
                            np.abs(np.diff(v[:, :, t], axis=1)).mean()
                        ),
                    )
                )
    raw = pd.DataFrame(outputs)
    keys = [
        "haplotypes",
        "source",
        "method",
        "cutoff_years",
        "policy",
        "target",
        "mode",
    ]
    summary = (
        raw.groupby(keys, sort=True).called.agg(["count", "sum", "mean"]).reset_index()
    )
    summary = summary.rename(
        columns={"count": "n", "sum": "called", "mean": "region_call_rate"}
    )
    summary.to_csv(args.out / "region_calls.csv", index=False)
    raw.to_csv(args.out / "held_out_regions.csv.gz", index=False)
    pd.DataFrame(choices).to_csv(args.out / "training_thresholds.csv", index=False)
    diagnostics = pd.DataFrame(diagnostics)
    diagnostics.to_csv(args.out / "profile_diagnostics.csv", index=False)
    pd.DataFrame(errors).to_csv(args.out / "decoded_error.csv", index=False)
    figures(summary, diagnostics, args.out)
    record = dict(
        primary_endpoint="At least one call anywhere in each isolated 10 Mb region",
        sample_fixation="Retained in both training and evaluation",
        comparison="Regenerated panels with unchanged seed convention; not nested samples",
        null_calibration="Five folds: 400 neutral fit, 400 calibration, 200 held out per fold",
        selected_calibration="80 training selected regions; 20 held out per fold; no method selection",
        interpretation="Each cutoff and method evaluated separately; choosing among results is exploratory",
        seed=previous_cfg["seed_base"],
        inputs=provenance,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        helper_source_sha256={
            name: hashlib.sha256(
                Path(__file__).with_name(name).read_bytes()
            ).hexdigest()
            for name in ("fresh_region_metrics.py", "partial_sweep_cv.py")
        },
        outputs={
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in args.out.iterdir()
            if p.is_file() and p.name != "provenance.json"
        },
    )
    (args.out / "provenance.json").write_text(json.dumps(record, indent=2) + "\n")
    print(f"Complete: {args.out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
