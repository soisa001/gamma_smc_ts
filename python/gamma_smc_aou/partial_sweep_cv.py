"""Reduce regional false peaks at a specified partial-sweep sensitivity.

Five outer folds hold out entire simulations. Training uses only nonfixed
s=0.005 selected replicates and neutral replicates outside the test fold.
Peak scoring is translation invariant; the selected coordinate is used only
for training a temporal template and for localization labels, never to favor
that coordinate during scanning. All features derive from decoded frac_recent_T.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path

os.environ["MPLBACKEND"] = "Agg"
for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np
import pandas as pd
from scipy.ndimage import uniform_filter1d
from scipy.signal import find_peaks

from .fresh_region_metrics import load_profiles, runs, score_counts


def write_json(path, value):
    temporary = path.with_name("temporary_" + path.name)
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def make_folds(metadata, seed, n_folds):
    rng = np.random.default_rng(seed)
    folds = np.empty(len(metadata), dtype=int)
    for _, group in metadata.groupby(["mode", "s"], sort=True):
        ordered = group.sort_values("replicate").index.to_numpy()
        folds[rng.permutation(ordered)] = np.arange(len(ordered)) % n_folds
    return folds


def normalize(features, neutral_train, floor=1e-4):
    reference = features[neutral_train]
    mean = reference.mean(axis=0, dtype=np.float64)
    sd = np.maximum(reference.std(axis=0, dtype=np.float64), floor)
    return ((features - mean) / sd).astype(np.float32)


def fit_linear(features, selected_train, neutral_train, focal, ridge):
    """Shrinkage linear discriminant with six or twelve profile features."""
    normalized = normalize(features, neutral_train)
    # Sampling positions estimates covariance; evaluation units remain regions.
    reference = normalized[neutral_train, ::10].reshape(-1, features.shape[-1])
    covariance = np.cov(reference, rowvar=False)
    contrast = normalized[selected_train, focal].mean(axis=0, dtype=np.float64)
    weights = np.linalg.solve(covariance + ridge * np.eye(len(contrast)), contrast)
    weights /= max(float(np.sqrt(weights @ covariance @ weights)), 1e-12)
    scores = np.einsum("ijk,k->ij", normalized, weights, optimize=False)
    return scores, weights.tolist()


def summit_candidates(scores, separation):
    """Unthresholded maxima, including boundaries, with fixed separation.

    A summit is one discovery point. No wider true-target allowance is gained
    by increasing smoothing width or by lowering the calling threshold.
    Columns are strength, summit, first position, final position (inclusive).
    """
    rows = []
    for profile in scores:
        summits = (
            find_peaks(np.r_[-np.inf, profile, -np.inf], distance=separation)[0] - 1
        )
        summits = summits[np.argsort(-profile[summits], kind="stable")]
        rows.append(np.column_stack((profile[summits], summits, summits, summits)))
    return rows


def cluster_candidates(scores, reference, quantile, minimum, separation):
    threshold = np.quantile(scores[reference], quantile, axis=0)
    rows = []
    for profile in scores:
        excess = profile - threshold
        candidates = []
        for a, b in runs(excess > 0, minimum):
            summit = a + int(np.argmax(excess[a:b]))
            candidates.append((float(excess[a:b].sum()), summit, a, b - 1))
        retained = []
        for peak in sorted(candidates, key=lambda item: (-item[0], item[1])):
            if all(abs(peak[1] - other[1]) >= separation for other in retained):
                retained.append(peak)
        rows.append(np.asarray(retained, dtype=float).reshape(-1, 4))
    return rows


def baseline_candidates(values, reference, cutoff, alpha, minimum):
    counts = np.rint(values[:, :, cutoff] * 10000).astype(np.int32)
    n = len(reference)
    pvalues = np.empty(counts.shape)
    for position in range(counts.shape[1]):
        null = np.sort(counts[reference, position])
        pvalues[:, position] = (
            1 + n - np.searchsorted(null, counts[:, position], side="left")
        ) / (n + 1)
    peaks = []
    for profile, p in zip(values[:, :, cutoff], pvalues):
        found = []
        for a, b in runs(p < alpha, minimum):
            summit = a + int(np.argmax(profile[a:b]))
            found.append((1.0, summit, a, b - 1))
        peaks.append(np.asarray(found, dtype=float).reshape(-1, 4))
    return peaks


def threshold_scores(peaks, threshold, limit, low, high, neutral=False):
    called = peaks[peaks[:, 0] >= threshold]
    if limit:
        called = called[:limit]
    hit = int(not neutral and np.any((called[:, 1] >= low) & (called[:, 1] <= high)))
    overlap = int(
        not neutral and np.any((called[:, 3] >= low) & (called[:, 2] <= high))
    )
    if neutral:
        result = dict(
            tp=0,
            fp=len(called),
            fn=0,
            calls=len(called),
            fdp=float(len(called) > 0),
            f1=0.0,
            recall=0.0,
            detected=0,
        )
    else:
        result = score_counts(hit, len(called) - hit, 1 - hit)
    result.update(
        overlap_detected=overlap,
        overlap_fp=len(called) - overlap,
        overlap_fdp=(len(called) - overlap) / len(called) if len(called) else 0.0,
        overlap_f1=2 * overlap / (len(called) + 1) if not neutral else 0.0,
    )
    return result


def select_threshold(
    peaks,
    selected_train,
    neutral_train,
    target,
    limits,
    low,
    high,
    endpoint="localization",
):
    """Select using training labels only; minimize neutral-region call rate."""
    candidates = []
    for limit in limits:
        best_true = []
        for index in selected_train:
            p = peaks[index][:limit] if limit else peaks[index]
            good = (
                p[:, 0]
                if endpoint == "region_any_call"
                else p[(p[:, 1] >= low) & (p[:, 1] <= high), 0]
            )
            best_true.append(float(np.max(good)) if len(good) else -np.inf)
        ordered = np.sort(best_true)[::-1]
        for wanted in sorted(
            set(
                [target] + [p for p in (0.75, 0.8, 0.85, 0.9, 0.95, 1.0) if p >= target]
            )
        ):
            threshold = ordered[int(np.ceil(wanted * len(ordered))) - 1]
            if not np.isfinite(threshold):
                continue
            scored = [
                threshold_scores(peaks[i], threshold, limit, low, high)
                for i in selected_train
            ]
            tp = sum(row["tp"] for row in scored)
            fp = sum(row["fp"] for row in scored)
            detected = (
                sum(row["calls"] > 0 for row in scored)
                if endpoint == "region_any_call"
                else tp
            )
            null_calls = [
                int(np.count_nonzero(peaks[i][:, 0] >= threshold))
                for i in neutral_train
            ]
            if limit:
                null_calls = [min(n, limit) for n in null_calls]
            candidates.append(
                dict(
                    threshold=float(threshold),
                    call_limit=limit,
                    training_power=detected / len(selected_train),
                    training_pooled_fdr=fp / (tp + fp),
                    training_mean_fdr=float(np.mean([r["fdp"] for r in scored])),
                    training_neutral_any=float(np.mean(np.asarray(null_calls) > 0)),
                    training_neutral_calls=float(np.mean(null_calls)),
                )
            )
    if not candidates:
        return None
    return min(
        candidates,
        key=lambda c: (
            c["training_neutral_any"],
            c["training_neutral_calls"],
            -c["training_power"],
            c["training_pooled_fdr"],
        ),
    )


def summarize(frame):
    rows = []
    group_keys = ["policy", "training_target", "method", "mode", "s", "subset"]
    for key, group in frame.groupby(group_keys, sort=True):
        tp, fp, fn = group[["tp", "fp", "fn"]].sum()
        row = dict(zip(group_keys, key))
        row.update(
            n=len(group),
            detected=int(tp),
            power=float(group.detected.mean()),
            fdr_mean=float(group.fdp.mean()),
            fdr_pooled=float(fp / (tp + fp)) if tp + fp else None,
            f1_mean=float(group.f1.mean()),
            fp_mean=float(group.fp.mean()),
            any_call=float((group.calls > 0).mean()),
            overlap_power=float(group.overlap_detected.mean()),
            overlap_fdr_mean=float(group.overlap_fdp.mean()),
            overlap_f1_mean=float(group.overlap_f1.mean()),
            tp_total=int(tp),
            fp_total=int(fp),
            fn_total=int(fn),
        )
        rows.append(row)
    return pd.DataFrame(rows)


def method_specs(spec, cutoffs):
    methods = []
    for width in spec["smoothing_strides"]:
        for t in range(len(cutoffs)):
            methods.append(
                dict(
                    name=f"smooth_w{width}_T{cutoffs[t]}",
                    family="Smoothed cutoff",
                    width=width,
                    t=t,
                )
            )
    for width in spec["contrast_strides"]:
        for t in range(1, len(cutoffs)):
            methods.append(
                dict(
                    name=f"contrast_w{width}_T{cutoffs[t]}",
                    family="Local contrast",
                    width=width,
                    t=t,
                )
            )
    for width in spec["lda_smoothing_strides"]:
        methods.append(
            dict(name=f"joint_time_w{width}", family="Joint time", width=width)
        )
        methods.append(
            dict(
                name=f"joint_shape_w{width}", family="Joint time and shape", width=width
            )
        )
    for t in range(1, len(cutoffs)):
        methods.append(
            dict(
                name=f"cluster_mass_T{cutoffs[t]}", family="Cluster mass", width=1, t=t
            )
        )
    return methods


def make_figure(summary, methods, out, target_s):
    import matplotlib
    import matplotlib.pyplot as plt

    matplotlib.rcParams.update(
        {"font.size": 18, "axes.labelsize": 20, "pdf.fonttype": 42, "ps.fonttype": 42}
    )
    families = {m["name"]: m["family"] for m in methods}
    colors = dict(zip(sorted(set(families.values())), plt.cm.tab10.colors))
    b = summary[
        (summary["mode"] == "selected")
        & (summary.s == target_s)
        & (summary.subset == "nonfixed")
    ]
    null = summary[(summary["mode"] == "neutral") & (summary.subset == "all")]
    b = b.merge(
        null[["policy", "training_target", "method", "any_call"]],
        on=["policy", "training_target", "method"],
        suffixes=("", "_neutral"),
    )
    fig, ax = plt.subplots(figsize=(11, 8.5), layout="constrained")
    for family, color in colors.items():
        group = b[(b.policy == "method") & b.method.map(families).eq(family)]
        # Require all nonfixed replicates; infeasible fold/methods are not silently dropped.
        group = group[group.n == b[b.policy == "adaptive"].n.max()]
        ax.scatter(
            group.any_call,
            group.any_call_neutral,
            color=color,
            s=55,
            alpha=0.7,
            label=family,
        )
    for policy, marker, color in [
        ("adaptive", "*", "black"),
        ("baseline", "X", "#cf2b27"),
    ]:
        group = b[b.policy == policy]
        for _, row in group.iterrows():
            ax.scatter(
                row.any_call,
                row.any_call_neutral,
                marker=marker,
                color=color,
                s=230,
                zorder=4,
            )
            label = (
                "Baseline"
                if policy == "baseline"
                else f"Target {row.training_target:.0%}"
            )
            ax.annotate(
                label,
                (row.any_call, row.any_call_neutral),
                xytext=(
                    (-110, -28)
                    if row.training_target == 0.7
                    else (75, 30)
                    if row.training_target == 0.8
                    else (40, -30)
                ),
                textcoords="offset points",
                fontsize=15,
                arrowprops=dict(arrowstyle="-", color=color, lw=1),
            )
    ax.axvline(0.7, color="gray", linestyle="--", linewidth=1.5)
    ax.set(
        xlim=(0, 1.04),
        ylim=(0, 1.04),
        xlabel="Fraction of nonfixed selected regions called",
        ylabel="Neutral 10 Mb regions with at least one call",
        title="Power–neutral false-call trade-off",
    )
    ax.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 1.08),
        ncol=2,
        fontsize=13,
        frameon=False,
    )
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(out / "tradeoff.png", dpi=200)
    fig.savefig(out / "tradeoff.pdf")
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument("--spec", type=Path, required=True)
    parser.add_argument("--plots-only", action="store_true")
    args = parser.parse_args(argv)
    spec_bytes = args.spec.read_bytes()
    spec = json.loads(spec_bytes)
    root = args.study.resolve()
    out = (
        root
        / "analysis"
        / (
            "partial_sweep_region_cv"
            if spec.get("power_endpoint") == "region_any_call"
            else "partial_sweep_cv"
        )
    )
    out.mkdir(parents=True, exist_ok=True)
    if args.plots_only:
        provenance = json.loads((out / "provenance.json").read_text())
        if provenance["analysis_spec"] != spec:
            raise ValueError("Plot specification differs from completed analysis")
        for name, digest in provenance["outputs"].items():
            if hashlib.sha256((out / name).read_bytes()).hexdigest() != digest:
                raise ValueError(f"Corrupt completed output: {name}")
        summary = pd.read_csv(out / "summary.csv")
        summary["region_call_rate"] = summary["any_call"]
        summary["localization_rate"] = summary["power"]
        summary.to_csv(out / "summary.csv", index=False)
        methods = provenance["method_grid"]
        make_figure(summary, methods, out, spec["target_selection_coefficient"])
        provenance["primary_metric"] = (
            "region_call_rate / any_call, same anywhere-in-region endpoint in neutral and selected sets"
        )
        provenance["localization_metric"] = (
            "localization_rate / power is a separate diagnostic, not the primary sensitivity"
        )
        provenance["plot_source_sha256"] = hashlib.sha256(
            Path(__file__).read_bytes()
        ).hexdigest()
        provenance["outputs"] = {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in out.iterdir()
            if p.is_file() and p.name != "provenance.json"
        }
        write_json(out / "provenance.json", provenance)
        return 0
    cfg, fingerprint, positions, meta, counts, inputs, manifest_hash = load_profiles(
        root, spec["workers"]
    )
    values = counts.astype(np.float32) / cfg["haplotype_pairs"]
    del counts
    folds = make_folds(meta, spec["seed"], spec["outer_folds"])
    meta["fold"] = folds
    meta.to_csv(out / "fold_manifest.csv", index=False)
    widths = sorted(
        set(
            spec["smoothing_strides"]
            + spec["contrast_strides"]
            + spec["lda_smoothing_strides"]
            + [spec["background_strides"]]
        )
    )
    print("Preparing spatial features", flush=True)
    with ThreadPoolExecutor(max_workers=spec["workers"]) as pool:
        smoothed = dict(
            pool.map(
                lambda w: (w, uniform_filter1d(values, w, axis=1, mode="nearest")),
                widths,
            )
        )
    focal = int(np.flatnonzero(positions == cfg["focal_position_bp"])[0])
    half = spec["target_half_width_bp"] // cfg["stride_bp"]
    low, high = focal - half, focal + half
    methods = method_specs(spec, cfg["tmrca_cutoffs_years"])
    run_hash = hashlib.sha256(
        spec_bytes
        + Path(__file__).read_bytes()
        + Path(__file__).with_name("fresh_region_metrics.py").read_bytes()
        + json.dumps(inputs, sort_keys=True).encode()
    ).hexdigest()
    all_rows, choice_rows, coefficients = [], [], []

    for fold in range(spec["outer_folds"]):
        neutral_train = np.flatnonzero((meta["mode"] == "neutral") & (folds != fold))
        selected_train = np.flatnonzero(
            (meta["mode"] == "selected")
            & (meta.s == spec["target_selection_coefficient"])
            & (folds != fold)
            & (~meta.fixed if spec["exclude_sample_fixed_from_training"] else True)
        )
        test = np.flatnonzero(folds == fold)
        assert not set(test) & (set(selected_train) | set(neutral_train))
        fold_path = out / f"fold{fold}.json"
        if fold_path.exists():
            cached = json.loads(fold_path.read_text())
            if cached["run_hash"] == run_hash:
                payload = {
                    key: cached[key] for key in ("rows", "choices", "coefficients")
                }
                if (
                    hashlib.sha256(
                        json.dumps(payload, sort_keys=True).encode()
                    ).hexdigest()
                    != cached["payload_sha256"]
                ):
                    raise ValueError(f"Corrupt fold cache: {fold_path}")
                all_rows.extend(cached["rows"])
                choice_rows.extend(cached["choices"])
                coefficients.extend(cached["coefficients"])
                print(
                    f"Validated cached fold {fold + 1}/{spec['outer_folds']}",
                    flush=True,
                )
                continue
        print(
            f"Fold {fold + 1}/{spec['outer_folds']}: {len(selected_train)} partial training sweeps; {len(test)} held-out regions",
            flush=True,
        )
        baseline = baseline_candidates(
            values,
            neutral_train,
            cfg["tmrca_cutoffs_years"].index(spec["baseline_cutoff_years"]),
            spec["baseline_alpha"],
            spec["baseline_min_strides"],
        )

        def evaluate_method(method):
            width = method["width"]
            family = method["family"]
            weights = None
            if family in ("Smoothed cutoff", "Cluster mass"):
                scores = normalize(smoothed[width][:, :, method["t"]], neutral_train)
            elif family == "Local contrast":
                scores = normalize(
                    smoothed[width][:, :, method["t"]]
                    - smoothed[spec["background_strides"]][:, :, method["t"]],
                    neutral_train,
                )
            else:
                features = smoothed[width]
                if family == "Joint time and shape":
                    features = np.concatenate(
                        (features, features - smoothed[spec["background_strides"]]),
                        axis=2,
                    )
                scores, weights = fit_linear(
                    features,
                    selected_train,
                    neutral_train,
                    focal,
                    spec["lda_covariance_ridge"],
                )
            if family == "Cluster mass":
                peaks = cluster_candidates(
                    scores,
                    neutral_train,
                    spec["cluster_forming_quantile"],
                    spec["cluster_min_strides"],
                    spec["peak_separation_strides"],
                )
            else:
                peaks = summit_candidates(scores, spec["peak_separation_strides"])
            results = []
            for target in spec["power_targets_in_training"]:
                choice = select_threshold(
                    peaks,
                    selected_train,
                    neutral_train,
                    target,
                    spec["call_limits"],
                    low,
                    high,
                    endpoint=spec.get("power_endpoint", "localization"),
                )
                if choice is not None:
                    results.append((target, choice))
            return method, peaks, results, weights

        with ThreadPoolExecutor(max_workers=spec["workers"]) as pool:
            evaluated = list(pool.map(evaluate_method, methods))
        fold_rows, fold_choices, fold_coefficients = [], [], []

        def append_scores(policy, name, target, choice, peaks):
            for index in test:
                item = meta.iloc[index]
                scores = threshold_scores(
                    peaks[index],
                    choice["threshold"],
                    choice["call_limit"],
                    low,
                    high,
                    item["mode"] == "neutral",
                )
                fold_rows.append(
                    dict(
                        task_id=item.task_id,
                        fold=fold,
                        mode=item["mode"],
                        s=float(item.s),
                        sample_af=float(item.sample_af),
                        fixed=bool(item.fixed),
                        policy=policy,
                        method=name,
                        training_target=target,
                        **scores,
                    )
                )

        append_scores(
            "baseline",
            f"T{spec['baseline_cutoff_years']}_p{spec['baseline_alpha']}_min{spec['baseline_min_strides']}",
            0.0,
            dict(threshold=0.0, call_limit=0),
            baseline,
        )
        for method, peaks, results, weights in evaluated:
            if weights is not None:
                fold_coefficients.append(
                    dict(fold=fold, method=method["name"], weights=weights)
                )
            for target, choice in results:
                fold_choices.append(
                    dict(
                        fold=fold,
                        method=method["name"],
                        training_target=target,
                        **choice,
                    )
                )
                append_scores("method", method["name"], target, choice, peaks)
        for target in spec["power_targets_in_training"]:
            eligible = [
                (method, peaks, choice)
                for method, peaks, results, _ in evaluated
                for wanted, choice in results
                if wanted == target
            ]
            if not eligible:
                raise ValueError(
                    f"No training rule reaches {target} power in fold {fold}"
                )
            method, peaks, choice = min(
                eligible,
                key=lambda x: (
                    x[2]["training_neutral_any"],
                    x[2]["training_neutral_calls"],
                    -x[2]["training_power"],
                    x[2]["training_pooled_fdr"],
                    x[0]["name"],
                ),
            )
            fold_choices.append(
                dict(
                    fold=fold,
                    method="ADAPTIVE",
                    selected_method=method["name"],
                    training_target=target,
                    **choice,
                )
            )
            append_scores("adaptive", "ADAPTIVE", target, choice, peaks)
            print(
                f"  target {target:.0%}: {method['name']}, train power {choice['training_power']:.1%}, neutral any-call {choice['training_neutral_any']:.1%}",
                flush=True,
            )
        cached = dict(
            run_hash=run_hash,
            rows=fold_rows,
            choices=fold_choices,
            coefficients=fold_coefficients,
        )
        cached["payload_sha256"] = hashlib.sha256(
            json.dumps(
                {key: cached[key] for key in ("rows", "choices", "coefficients")},
                sort_keys=True,
            ).encode()
        ).hexdigest()
        write_json(fold_path, cached)
        all_rows.extend(fold_rows)
        choice_rows.extend(fold_choices)
        coefficients.extend(fold_coefficients)

    frame = pd.DataFrame(all_rows)
    frame["subset"] = "all"
    nonfixed = frame[(frame["mode"] == "selected") & ~frame.fixed].copy()
    nonfixed["subset"] = "nonfixed"
    strata = frame[
        (frame["mode"] == "selected")
        & (frame.s == spec["target_selection_coefficient"])
    ].copy()
    strata["subset"] = pd.cut(
        strata.sample_af,
        [-0.001, 0.5, 0.8, 0.999999, 1.0],
        labels=["AF<=0.5", "0.5<AF<=0.8", "0.8<AF<1", "AF=1"],
    ).astype(str)
    combined = pd.concat([frame, nonfixed, strata], ignore_index=True)
    summary = summarize(combined)
    summary["region_call_rate"] = summary["any_call"]
    summary["localization_rate"] = summary["power"]
    frame.to_csv(out / "out_of_fold_per_replicate.csv.gz", index=False)
    summary.to_csv(out / "summary.csv", index=False)
    pd.DataFrame(choice_rows).to_csv(out / "training_choices.csv", index=False)
    write_json(out / "linear_coefficients.json", coefficients)
    make_figure(summary, methods, out, spec["target_selection_coefficient"])
    write_json(
        out / "provenance.json",
        dict(
            created_utc=datetime.now(timezone.utc).isoformat(),
            study_path=str(root),
            study_fingerprint=fingerprint,
            sample_manifest_sha256=manifest_hash,
            run_hash=run_hash,
            analysis_spec=spec,
            method_grid=methods,
            source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            profile_hashes=inputs,
            validation="Five outer folds of whole replicates; test fold excluded from transformations, weights, thresholds, and method choice",
            limitations="Internal exploratory cross-validation of previously summarized simulations; no formal FDR-control guarantee; no independent confirmatory sample",
            outputs={
                p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                for p in out.iterdir()
                if p.is_file() and p.name != "provenance.json"
            },
        ),
    )
    print(f"Complete: {out}", flush=True)
    print(
        summary[
            (summary.policy.isin(["baseline", "adaptive"]))
            & (summary.s == 0.005)
            & (summary.subset == "nonfixed")
        ][
            [
                "policy",
                "training_target",
                "n",
                "region_call_rate",
                "localization_rate",
                "fdr_mean",
                "fdr_pooled",
                "f1_mean",
                "overlap_power",
            ]
        ].to_string(index=False),
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
