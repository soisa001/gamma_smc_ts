"""Localization FDR/F1 from the saved fresh-study decoded profiles.

No simulation or decoding is performed. The target interval is an explicitly
chosen localization tolerance, not a claim about the causal sweep footprint.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import io
import json
import os
from pathlib import Path
import platform

os.environ["MPLBACKEND"] = "Agg"
for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np
import pandas as pd


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def aligned_pvalues(neutral, selected):
    """Upper-tail ranks at matched positions; remove self for neutral ranks.

    Arrays have shape (replicate, position, cutoff). Work with integer pair
    counts so ties are identical to the underlying discrete statistic.
    """
    n = len(neutral)
    selected_p = np.empty(selected.shape, dtype=float)
    neutral_p = np.empty(neutral.shape, dtype=float)
    for position in range(neutral.shape[1]):
        for cutoff in range(neutral.shape[2]):
            reference = np.sort(neutral[:, position, cutoff])
            selected_p[:, position, cutoff] = (
                1
                + n
                - np.searchsorted(reference, selected[:, position, cutoff], side="left")
            ) / (n + 1)
            # n - rank includes this neutral replicate; equivalently,
            # 1 + the number of OTHER replicates at least as large.
            neutral_p[:, position, cutoff] = (
                n
                - np.searchsorted(reference, neutral[:, position, cutoff], side="left")
            ) / n
    return neutral_p, selected_p


def runs(mask, minimum=1):
    padded = np.pad(np.asarray(mask, dtype=np.int8), (1, 1))
    changes = np.diff(padded)
    starts, ends = np.flatnonzero(changes == 1), np.flatnonzero(changes == -1)
    return [(int(a), int(b)) for a, b in zip(starts, ends) if b - a >= minimum]


def score_strides(mask, target):
    tp = int(np.count_nonzero(mask & target))
    fp = int(np.count_nonzero(mask & ~target))
    fn = int(np.count_nonzero(~mask & target))
    return score_counts(tp, fp, fn)


def score_peaks(found, target):
    # One selected locus is one event. A maximum of one called peak can match;
    # additional calls, including duplicate overlapping calls, are false calls.
    tp = int(any(target[a:b].any() for a, b in found))
    return score_counts(tp, len(found) - tp, 1 - tp)


def score_counts(tp, fp, fn):
    return dict(
        tp=tp,
        fp=fp,
        fn=fn,
        calls=tp + fp,
        fdp=fp / (tp + fp) if tp + fp else 0.0,
        f1=2 * tp / (2 * tp + fp + fn) if 2 * tp + fp + fn else 0.0,
        recall=tp / (tp + fn) if tp + fn else 0.0,
        detected=int(tp > 0),
    )


def load_profiles(root, workers):
    manifest_bytes = (root / "sample_manifest.json").read_bytes()
    manifest = json.loads(manifest_bytes)
    study = json.loads((root / "manifest.json").read_text())
    cfg, fingerprint = study["config"], study["fingerprint"]
    positions = np.arange(0, cfg["scored_length_bp"], cfg["stride_bp"])
    columns = [f"frac_recent_{t}" for t in cfg["tmrca_cutoffs_years"]]
    if len({item["task_id"] for item in manifest}) != len(manifest):
        raise ValueError("Duplicate task IDs in sample manifest")

    def read(item):
        directory = root / item["task_id"]
        record = json.loads((directory / "complete.json").read_text())
        if record["fingerprint"] != fingerprint:
            raise ValueError(f"Scientific fingerprint mismatch: {directory}")
        for key in ("mode", "s", "replicate"):
            if item[key] != record["task"][key]:
                raise ValueError(f"Manifest/receipt task mismatch: {directory}")
        data = (directory / "frac_recent.tsv").read_bytes()
        digest = sha256(data)
        artifact = record["artifacts"]["frac_recent.tsv"]
        if digest != artifact["sha256"] or len(data) != artifact["bytes"]:
            raise ValueError(f"Corrupt decoded profile: {directory}")
        frame = pd.read_csv(io.BytesIO(data), sep="\t")
        values = frame[columns].to_numpy(float)
        if not np.array_equal(frame.position_0based, positions):
            raise ValueError(f"Position mismatch: {directory}")
        if (
            not np.isfinite(values).all()
            or np.any(values < 0)
            or np.any(values > 1)
            or np.any(np.diff(values, axis=1) < -1e-12)
        ):
            raise ValueError(f"Invalid frac_recent values: {directory}")
        counts = values * cfg["haplotype_pairs"]
        if not np.allclose(counts, np.rint(counts), atol=1e-6, rtol=0):
            raise ValueError(f"Not an integer pair count: {directory}")
        metadata = dict(item, sample_af=record["sample_af"], fixed=record["fixed"])
        return metadata, np.rint(counts).astype(np.int32), digest

    print(
        f"Validating {len(manifest)} decoded profiles with {workers} readers",
        flush=True,
    )
    with ThreadPoolExecutor(max_workers=workers) as pool:
        loaded = list(pool.map(read, manifest))
    metadata = pd.DataFrame([row[0] for row in loaded])
    counts = np.stack([row[1] for row in loaded])
    if metadata[metadata["mode"] == "neutral"].shape[0] != cfg[
        "neutral_replicates"
    ] or set(metadata.loc[metadata["mode"] == "selected", "s"]) != set(
        cfg["selection_coefficients"]
    ):
        raise ValueError("Study replicate inventory does not match configuration")
    if not (
        metadata[metadata["mode"] == "selected"].groupby("s").size()
        == cfg["selected_replicates_per_coefficient"]
    ).all():
        raise ValueError("Incomplete selected arm")
    inputs = [{"task_id": row[0]["task_id"], "sha256": row[2]} for row in loaded]
    return cfg, fingerprint, positions, metadata, counts, inputs, sha256(manifest_bytes)


def summarize(frame, seed):
    rows = []
    rng = np.random.default_rng(seed)
    keys = ["metric", "min_strides", "alpha", "cutoff_years", "half_width_bp", "s"]
    for key, group in frame.groupby(keys, sort=True):
        tp, fp, fn = group[["tp", "fp", "fn"]].sum()
        row = dict(zip(keys, key))
        row.update(
            n=len(group),
            detected=int(group.detected.sum()),
            detection=group.detected.mean(),
            fdr_mean=group.fdp.mean(),
            f1_mean=group.f1.mean(),
            recall_mean=group.recall.mean(),
            fp_mean=group.fp.mean(),
            calls_mean=group.calls.mean(),
            no_calls=int((group.calls == 0).sum()),
            tp_total=int(tp),
            fp_total=int(fp),
            fn_total=int(fn),
            fdp_pooled=fp / (tp + fp) if tp + fp else np.nan,
            f1_pooled=2 * tp / (2 * tp + fp + fn),
        )
        # Resample entire selected replicates, never linked positions.
        indices = rng.integers(0, len(group), size=(2000, len(group)))
        for metric in ("fdp", "f1"):
            means = group[metric].to_numpy()[indices].mean(axis=1)
            row[f"{metric}_ci_low"], row[f"{metric}_ci_high"] = np.quantile(
                means, [0.025, 0.975]
            )
        rows.append(row)
    return pd.DataFrame(rows)


def figures(summary, out):
    import matplotlib
    import matplotlib.pyplot as plt

    matplotlib.rcParams.update(
        {
            "font.size": 17,
            "axes.labelsize": 19,
            "axes.titlesize": 19,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    for metric, minimum, title in (
        ("stride", 1, "All significant strides"),
        ("peak", 5, "Peaks of at least 5 consecutive strides"),
    ):
        block = summary[
            (summary.metric == metric)
            & (summary.min_strides == minimum)
            & (summary.alpha == 0.05)
            & (summary.half_width_bp == 100000)
        ]
        fig, axes = plt.subplots(1, 2, figsize=(11, 8.5), layout="constrained")
        for ax, column, label in zip(axes, ("fdr_mean", "f1_mean"), ("FDR", "F1")):
            table = block.pivot(index="s", columns="cutoff_years", values=column)
            values = table.to_numpy()
            ax.imshow(values, vmin=0, vmax=1, cmap="viridis", aspect="auto")
            for y, x in np.ndindex(values.shape):
                ax.text(
                    x,
                    y,
                    f"{values[y, x]:.2f}",
                    ha="center",
                    va="center",
                    fontsize=13,
                    color="white" if values[y, x] < 0.52 else "black",
                )
            ax.set_xticks(
                range(len(table.columns)), [str(t // 1000) for t in table.columns]
            )
            ax.set_yticks(range(len(table.index)), [f"{s:.3f}" for s in table.index])
            ax.set(xlabel="TMRCA cutoff (kya)", title=label)
        axes[0].set_ylabel("Selection coefficient")
        fig.suptitle(
            f"{title}\np < 0.05; target within ±100 kb of selected allele", fontsize=18
        )
        fig.savefig(out / f"{metric}_fdr_f1.png", dpi=200)
        fig.savefig(out / f"{metric}_fdr_f1.pdf")
        plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=20)
    parser.add_argument("--alphas", type=float, nargs="+", default=[0.05, 0.01])
    parser.add_argument(
        "--half-widths", type=int, nargs="+", default=[50000, 100000, 250000, 500000]
    )
    args = parser.parse_args(argv)
    if not 1 <= args.workers <= 20:
        parser.error("Use 1 to 20 workers")
    root = args.study.resolve()
    out = root / "analysis" / "region"
    out.mkdir(parents=True, exist_ok=True)
    cfg, fingerprint, positions, metadata, counts, inputs, manifest_hash = (
        load_profiles(root, args.workers)
    )
    is_neutral = metadata["mode"].to_numpy() == "neutral"
    neutral, selected = counts[is_neutral], counts[~is_neutral]
    if min(args.alphas) <= 1 / len(neutral):
        raise ValueError(
            "Neutral leave-one-out p-values cannot resolve this strict alpha"
        )
    neutral_p, selected_p = aligned_pvalues(neutral, selected)
    selected_meta = metadata[~is_neutral].reset_index(drop=True)
    neutral_meta = metadata[is_neutral].reset_index(drop=True)

    # Match the published focal analysis exactly before extending to the region.
    focal = int(np.flatnonzero(positions == cfg["focal_position_bp"])[0])
    published = pd.read_csv(root / "analysis" / "selected_pvalues.csv").set_index(
        "task_id"
    )
    for i, cutoff in enumerate(cfg["tmrca_cutoffs_years"]):
        if not np.allclose(
            selected_p[:, focal, i],
            published.loc[selected_meta.task_id, f"p_{cutoff}"],
            atol=1e-14,
            rtol=0,
        ):
            raise ValueError(
                "Regional calibration fails to reproduce published focal p-values"
            )
    print("All 6,000 focal p-values reproduced; scoring full regions", flush=True)

    selected_rows, neutral_rows = [], []
    targets = {
        width: np.abs(positions - cfg["focal_position_bp"]) <= width
        for width in args.half_widths
    }
    for i, cutoff in enumerate(cfg["tmrca_cutoffs_years"]):
        for alpha in args.alphas:
            for j, meta in selected_meta.iterrows():
                mask = selected_p[j, :, i] < alpha
                found = runs(mask)
                long = [pair for pair in found if pair[1] - pair[0] >= 5]
                base = dict(
                    task_id=meta.task_id,
                    s=meta.s,
                    sample_af=meta.sample_af,
                    fixed=meta.fixed,
                    cutoff_years=cutoff,
                    alpha=alpha,
                )
                for width, target in targets.items():
                    for metric, minimum, score in (
                        ("stride", 1, score_strides(mask, target)),
                        ("peak", 1, score_peaks(found, target)),
                        ("peak", 5, score_peaks(long, target)),
                    ):
                        selected_rows.append(
                            dict(
                                base,
                                metric=metric,
                                min_strides=minimum,
                                half_width_bp=width,
                                **score,
                            )
                        )
            for j, meta in neutral_meta.iterrows():
                mask = neutral_p[j, :, i] < alpha
                found = runs(mask)
                for metric, minimum, calls in (
                    ("stride", 1, int(mask.sum())),
                    ("peak", 1, len(found)),
                    ("peak", 5, sum(b - a >= 5 for a, b in found)),
                ):
                    neutral_rows.append(
                        dict(
                            task_id=meta.task_id,
                            cutoff_years=cutoff,
                            alpha=alpha,
                            metric=metric,
                            min_strides=minimum,
                            false_calls=calls,
                            any_call=int(calls > 0),
                        )
                    )
        print(f"Scored cutoff {cutoff} years", flush=True)
    raw = pd.DataFrame(selected_rows)
    summary = summarize(raw, cfg["seed_base"])
    null_raw = pd.DataFrame(neutral_rows)
    null_summary = (
        null_raw.groupby(["metric", "min_strides", "alpha", "cutoff_years"])
        .agg(
            n=("any_call", "size"),
            any_call_rate=("any_call", "mean"),
            false_calls_mean=("false_calls", "mean"),
            false_calls_median=("false_calls", "median"),
        )
        .reset_index()
    )
    for name, frame in (
        ("selected_per_replicate.csv.gz", raw),
        ("selected_summary.csv", summary),
        ("neutral_per_replicate.csv.gz", null_raw),
        ("neutral_summary.csv", null_summary),
    ):
        temporary = out / ("temporary_" + name)
        frame.to_csv(temporary, index=False)
        temporary.replace(out / name)
    figures(summary, out)
    provenance = dict(
        created_utc=datetime.now(timezone.utc).isoformat(),
        study_path=str(root),
        fingerprint=fingerprint,
        sample_manifest_sha256=manifest_hash,
        source_sha256=sha256(Path(__file__).read_bytes()),
        python=platform.python_version(),
        numpy=np.__version__,
        pandas=pd.__version__,
        workers=args.workers,
        statistic="frac_recent_T from decoded posterior mean TMRCA",
        cutoffs=cfg["tmrca_cutoffs_years"],
        alphas=args.alphas,
        half_widths_bp=args.half_widths,
        selected_n=len(selected),
        neutral_n=len(neutral),
        calibration="Position-matched independent neutral replicates; ties included; strict p < alpha",
        selected_p="(1 + count of neutral values >= selected value) / (N + 1)",
        neutral_p="Leave one replicate out: (1 + count of other neutral values >= value) / N",
        target="Profile positions within inclusive ±half_width of 5 Mb; localization tolerance, not causal truth",
        fdr="Mean replicate FP / max(TP + FP, 1); no-call replicates contribute zero",
        f1="Mean replicate 2 TP / (2 TP + FP + FN); no-call selected replicates contribute zero",
        peak_matching="Consecutive significant strides only; at most one peak matches the single target event",
        bootstrap="2,000 selected-replicate resamples; seed from study; intervals conditional on estimated null",
        scope="Isolated 10 Mb regions; no genome-wide error control or prevalence-mixed FDR claim",
        notes="Wholly neutral FDR equals probability of any call. Neutral F1 is not a useful endpoint.",
        profile_hashes=inputs,
    )
    provenance["outputs"] = {
        p.name: sha256(p.read_bytes())
        for p in sorted(out.iterdir())
        if p.name != "provenance.json" and p.is_file()
    }
    temporary = out / "temporary_provenance.json"
    temporary.write_text(json.dumps(provenance, indent=2) + "\n")
    temporary.replace(out / "provenance.json")
    print(f"Complete: {out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
