"""Matched true/decoded peak rules, including the 0.001 rank boundary.

The inclusive p <= alpha convention makes the neutral leave-one-out minimum
1/1000 eligible at alpha=0.001. Selected minimum is 1/1001. No zero p-values or
independence assumption across linked strides is introduced.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
import multiprocessing
import os
from pathlib import Path

os.environ["MPLBACKEND"] = "Agg"
for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np
import pandas as pd

from .fresh_region_metrics import aligned_pvalues, load_profiles, score_counts

DATA = {}
MIN_RUNS = [1, 3, 5, 10, 20]
STEPS = [1, 2, 5]
GAPS = [0, 1]
ALPHAS = [0.05, 0.01, 0.001]
RADII = [50000, 100000, 250000, 500000]


def peak_groups(mask, values, max_gap):
    hits = np.flatnonzero(mask)
    groups = (
        np.split(hits, np.flatnonzero(np.diff(hits) > max_gap + 1) + 1)
        if len(hits)
        else []
    )
    return [
        (len(g), int(g[np.argmax(values[g])]), int(g[0]), int(g[-1])) for g in groups
    ]


def grid_job(job):
    source, t, alpha = job
    item = DATA[source]
    values, pvalues, metadata, cutoff = (
        item["values"],
        item["pvalues"],
        item["meta"],
        item["cutoffs"][t],
    )
    rows, distances = [], []
    for step in STEPS:
        positions = item["positions"][::step]
        v, p = values[:, ::step, t], pvalues[:, ::step, t]
        for gap in GAPS:
            for i, meta in metadata.iterrows():
                found = peak_groups(p[i] <= alpha, v[i], gap)
                neutral = meta["mode"] == "neutral"
                for minimum in MIN_RUNS:
                    retained = [peak for peak in found if peak[0] >= minimum]
                    summit_distance = np.array(
                        [abs(positions[peak[1]] - item["focal"]) for peak in retained]
                    )
                    base = dict(
                        source=source,
                        cutoff_years=cutoff,
                        alpha=alpha,
                        stride_bp=step * 10000,
                        max_gap_strides=gap,
                        min_significant_strides=minimum,
                        task_id=meta.task_id,
                        mode=meta["mode"],
                        sample_af=float(meta.sample_af),
                        fixed=bool(meta.fixed),
                        focal_significant=int(
                            p[i, np.flatnonzero(positions == item["focal"])[0]] <= alpha
                        ),
                    )
                    tp = int(not neutral and np.any(summit_distance <= 100000))
                    score = score_counts(tp, len(retained) - tp, 1 - tp)
                    if neutral:
                        score.update(fn=0, f1=0.0, recall=0.0)
                    extra = {}
                    for radius in RADII:
                        hit = int(not neutral and np.any(summit_distance <= radius))
                        overlap = int(
                            not neutral
                            and any(
                                positions[b] >= item["focal"] - radius
                                and positions[a] <= item["focal"] + radius
                                for _, _, a, b in retained
                            )
                        )
                        extra[f"hit_{radius}"] = hit
                        extra[f"overlap_{radius}"] = overlap
                    rows.append(dict(base, **score, **extra))
                    # Counts describe spatial decay without calling the tail biologically false.
                    bins = np.histogram(
                        summit_distance,
                        [0, 50000, 100000, 250000, 500000, 1000000, 2000000, 5000001],
                    )[0]
                    distances.append(
                        dict(base, **{f"bin_{j}": int(n) for j, n in enumerate(bins)})
                    )
    return pd.DataFrame(rows), pd.DataFrame(distances)


def aggregate(raw):
    keys = [
        "source",
        "cutoff_years",
        "alpha",
        "stride_bp",
        "max_gap_strides",
        "min_significant_strides",
        "mode",
        "subset",
    ]
    rows = []
    for key, group in raw.groupby(keys, sort=True):
        tp, fp = group[["tp", "fp"]].sum()
        row = dict(zip(keys, key))
        row.update(
            n=len(group),
            power=float(group.detected.mean()),
            fdr_mean=float(group.fdp.mean()),
            fdr_pooled=float(fp / (tp + fp)) if tp + fp else np.nan,
            f1_mean=float(group.f1.mean()),
            fp_mean=float(group.fp.mean()),
            any_call=float((group.calls > 0).mean()),
            focal_significant=float(group.focal_significant.mean()),
            tp_total=int(tp),
            fp_total=int(fp),
        )
        for radius in RADII:
            hits = group[f"hit_{radius}"].to_numpy()
            overlap = group[f"overlap_{radius}"].to_numpy()
            calls = group.calls.to_numpy()
            row[f"power_{radius}"] = float(hits.mean())
            row[f"overlap_power_{radius}"] = float(overlap.mean())
            row[f"overlap_fdr_{radius}"] = float(
                np.mean((calls - overlap) / np.maximum(calls, 1))
            )
        rows.append(row)
    return pd.DataFrame(rows)


def figures(summary, out):
    import matplotlib
    import matplotlib.pyplot as plt

    matplotlib.rcParams.update(
        {"font.size": 18, "axes.labelsize": 20, "pdf.fonttype": 42, "ps.fonttype": 42}
    )
    b = summary[(summary["mode"] == "selected") & (summary.subset == "nonfixed")]
    keys = [
        "source",
        "cutoff_years",
        "alpha",
        "stride_bp",
        "max_gap_strides",
        "min_significant_strides",
    ]
    neutral = summary[(summary["mode"] == "neutral") & (summary.subset == "all")]
    b = b.merge(neutral[keys + ["any_call"]], on=keys, suffixes=("", "_neutral"))
    fig, axes = plt.subplots(1, 2, figsize=(11, 8.5), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.12, right=0.98, bottom=0.13, top=0.80, wspace=0.13)
    colors = {a: c for a, c in zip(ALPHAS, ["#4277ac", "#dc8039", "#42804d"])}
    for ax, source in zip(axes, ["decoded", "truth"]):
        for alpha in ALPHAS:
            bb = b[(b.source == source) & (b.alpha == alpha)]
            ax.scatter(
                bb.any_call,
                bb.any_call_neutral,
                color=colors[alpha],
                s=30,
                alpha=0.5,
                label=f"p ≤ {alpha:g}",
            )
        ax.axvline(0.7, color="gray", linestyle="--")
        ax.set(
            title=source.capitalize(),
            xlim=(0, 1.02),
            ylim=(0, 1.02),
        )
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Fraction of neutral regions called")
    fig.supxlabel("Fraction of nonfixed selected regions called", fontsize=20, y=0.025)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.925),
        ncol=3,
        frameon=False,
        fontsize=15,
    )
    fig.suptitle(
        "Region calls at s = 0.005",
        fontsize=20,
        y=0.99,
    )
    fig.savefig(out / "truth_decoded_tradeoff.png", dpi=200)
    fig.savefig(out / "truth_decoded_tradeoff.pdf")
    plt.close(fig)


def distance_statistics(out):
    """Average each region within distance bins before estimating uncertainty."""
    import matplotlib.pyplot as plt

    records = []
    for source, item in DATA.items():
        distance = np.abs(item["positions"] - item["focal"])
        bins = np.where(distance == 0, 0, (distance - 1) // 50000 + 1)
        for b in np.unique(bins):
            take = bins == b
            per_region = item["values"][:, take].mean(axis=1) / 10000
            for mode in ["neutral", "selected"]:
                use = item["meta"]["mode"].to_numpy() == mode
                if mode == "selected":
                    use &= ~item["meta"].fixed.to_numpy()
                for t, cutoff in enumerate(item["cutoffs"]):
                    v = per_region[use, t]
                    records.append(
                        dict(
                            source=source,
                            distance_bp=float(distance[take].mean()),
                            mode=mode,
                            cutoff_years=cutoff,
                            n=len(v),
                            mean=float(v.mean()),
                            sem=float(v.std(ddof=1) / np.sqrt(len(v))),
                        )
                    )
    if records:
        frame = pd.DataFrame(records)
        frame.to_csv(out / "distance_statistics.csv", index=False)
    else:
        frame = pd.read_csv(out / "distance_statistics.csv")
    fig, axes = plt.subplots(1, 2, figsize=(11, 8.5), sharey=True)
    fig.subplots_adjust(left=0.12, right=0.98, bottom=0.13, top=0.80, wspace=0.13)
    for ax, source in zip(axes, ["decoded", "truth"]):
        for cutoff, color in [(20000, "#4277ac"), (50000, "#b65832")]:
            block = frame[(frame.source == source) & (frame.cutoff_years == cutoff)]
            selected = block[block["mode"] == "selected"].set_index("distance_bp")
            neutral = block[block["mode"] == "neutral"].set_index("distance_bp")
            x = selected.index.to_numpy() / 1e6
            delta = selected["mean"].to_numpy() - neutral["mean"].to_numpy()
            error = 1.96 * np.sqrt(
                selected["sem"].to_numpy() ** 2 + neutral["sem"].to_numpy() ** 2
            )
            ax.plot(
                x, delta, color=color, linewidth=2.5, label=f"T = {cutoff // 1000} kya"
            )
            ax.fill_between(x, delta - error, delta + error, color=color, alpha=0.15)
        ax.axhline(0, color="gray", linewidth=1)
        ax.axvline(0.1, color="gray", linestyle="--", linewidth=1.3)
        ax.set(xlim=(0, 1), title=source.capitalize())
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Mean frac_recent_T excess over neutral")
    fig.supxlabel("Distance from selected allele (Mb)", fontsize=20, y=0.025)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.925),
        ncol=2,
        frameon=False,
    )
    fig.suptitle(
        "Nonfixed s = 0.005: decay of the regional signal", fontsize=18, y=0.99
    )
    fig.savefig(out / "signal_decay.png", dpi=200)
    fig.savefig(out / "signal_decay.pdf")
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=20)
    parser.add_argument("--plots-only", action="store_true")
    args = parser.parse_args(argv)
    root = args.study.resolve()
    out = root / "analysis" / "paired_peak_grid"
    out.mkdir(parents=True, exist_ok=True)
    if args.plots_only:
        provenance = json.loads((out / "provenance.json").read_text())
        for name, expected in provenance["outputs"].items():
            if hashlib.sha256((out / name).read_bytes()).hexdigest() != expected:
                raise ValueError(f"Output checksum mismatch: {out / name}")
        figures(pd.read_csv(out / "summary.csv"), out)
        distance_statistics(out)
        provenance["plot_source_sha256"] = hashlib.sha256(
            Path(__file__).read_bytes()
        ).hexdigest()
        provenance["primary_plot_endpoint"] = (
            "Same any-call endpoint in selected and neutral 10 Mb regions"
        )
        provenance["primary_metric"] = (
            "Fraction of 10 Mb regions with at least one qualifying call anywhere"
        )
        provenance["truth_definition"] = (
            "Secondary localization only: summit within +/-100 kb and interval "
            "overlap at four radii; neither is a causal footprint label"
        )
        provenance["outputs"] = {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in out.iterdir()
            if p.is_file() and p.name != "provenance.json"
        }
        (out / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
        return 0
    inputs = {}
    for source, name in [
        ("decoded", "frac_recent.tsv"),
        ("truth", "truth_frac_recent.tsv"),
    ]:
        cfg, fingerprint, positions, meta, counts, hashes, manifest_hash = (
            load_profiles(root, args.workers, name)
        )
        use = (meta["mode"] == "neutral") | (
            (meta["mode"] == "selected") & (meta.s == 0.005)
        )
        counts, meta = counts[use], meta[use].reset_index(drop=True)
        neutral = meta["mode"].to_numpy() == "neutral"
        null_p, selected_p = aligned_pvalues(counts[neutral], counts[~neutral])
        pvalues = np.empty(counts.shape)
        pvalues[neutral], pvalues[~neutral] = null_p, selected_p
        DATA[source] = dict(
            values=counts,
            pvalues=pvalues,
            meta=meta,
            cutoffs=cfg["tmrca_cutoffs_years"],
            focal=cfg["focal_position_bp"],
            positions=positions,
        )
        inputs[source] = hashes
        print(f"Loaded and calibrated {source}", flush=True)
    jobs = [(source, t, alpha) for source in DATA for t in range(6) for alpha in ALPHAS]
    all_rows, all_distances = [], []
    with ProcessPoolExecutor(
        max_workers=args.workers, mp_context=multiprocessing.get_context("fork")
    ) as pool:
        for job, (rows, distances) in zip(jobs, pool.map(grid_job, jobs)):
            all_rows.append(rows)
            all_distances.append(distances)
            print(
                f"Scored {job[0]} T={cfg['tmrca_cutoffs_years'][job[1]]} alpha={job[2]}",
                flush=True,
            )
    raw = pd.concat(all_rows, ignore_index=True)
    raw.to_csv(out / "per_replicate.csv.gz", index=False)
    raw["subset"] = "all"
    nonfixed = raw[(raw["mode"] == "selected") & ~raw.fixed].copy()
    nonfixed["subset"] = "nonfixed"
    summary = aggregate(pd.concat([raw, nonfixed], ignore_index=True))
    summary.to_csv(out / "summary.csv", index=False)
    distance = pd.concat(all_distances, ignore_index=True)
    keys = [
        "source",
        "cutoff_years",
        "alpha",
        "stride_bp",
        "max_gap_strides",
        "min_significant_strides",
        "mode",
    ]
    distance.groupby(keys)[[f"bin_{j}" for j in range(7)]].mean().reset_index().to_csv(
        out / "distance_profile.csv", index=False
    )
    # Quantify decoder disagreement with the true-statistic calls; not causal FDR.
    agreement = []
    for t, cutoff in enumerate(cfg["tmrca_cutoffs_years"]):
        for alpha in ALPHAS:
            d = DATA["decoded"]["pvalues"][:, :, t] <= alpha
            v = DATA["truth"]["pvalues"][:, :, t] <= alpha
            for mode in ["neutral", "selected"]:
                subset = DATA["decoded"]["meta"]["mode"].to_numpy() == mode
                both = int(np.count_nonzero(d[subset] & v[subset]))
                only_d = int(np.count_nonzero(d[subset] & ~v[subset]))
                only_v = int(np.count_nonzero(~d[subset] & v[subset]))
                agreement.append(
                    dict(
                        cutoff_years=cutoff,
                        alpha=alpha,
                        mode=mode,
                        both=both,
                        decoded_only=only_d,
                        truth_only=only_v,
                        jaccard=both / (both + only_d + only_v)
                        if both + only_d + only_v
                        else np.nan,
                    )
                )
    pd.DataFrame(agreement).to_csv(
        out / "truth_decoded_stride_agreement.csv", index=False
    )
    figures(summary, out)
    distance_statistics(out)
    provenance = dict(
        created_utc=datetime.now(timezone.utc).isoformat(),
        study_fingerprint=fingerprint,
        sample_manifest_sha256=manifest_hash,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        alpha_rule="p <= alpha on both sides; neutral self removed; ties retained",
        alphas=ALPHAS,
        reporting_strides_bp=[step * 10000 for step in STEPS],
        run_requirements=MIN_RUNS,
        max_gap_strides=GAPS,
        statistic_sources=[
            "decoded posterior-mean frac_recent_T",
            "true TMRCA frac_recent_T",
        ],
        scope="Exploratory grid; 100 selected s=.005 and 1000 neutral; no formal FDR-control claim",
        primary_metric="Fraction of 10 Mb regions with at least one qualifying call anywhere",
        truth_definition="Secondary localization only: summit within +/-100 kb and interval overlap at four radii; neither is a causal footprint label",
        downsampling="Reporting stride changes use saved 10 kb profiles; the decoder was not rerun",
        profile_hashes=inputs,
        outputs={
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in out.iterdir()
            if p.is_file() and p.name != "provenance.json"
        },
    )
    (out / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"Complete: {out}", flush=True)
    rule_keys = [
        "source",
        "cutoff_years",
        "alpha",
        "stride_bp",
        "max_gap_strides",
        "min_significant_strides",
    ]
    neutral_rates = summary[(summary["mode"] == "neutral") & (summary.subset == "all")][
        rule_keys + ["any_call"]
    ].rename(columns={"any_call": "neutral_any_call"})
    for source in DATA:
        b = summary[
            (summary.source == source)
            & (summary["mode"] == "selected")
            & (summary.subset == "all")
            & (summary.any_call >= 0.7)
        ].merge(neutral_rates, on=rule_keys, validate="one_to_one")
        print(
            source, "exploratory rules with >=70% selected regions called", flush=True
        )
        print(
            b.sort_values("neutral_any_call")
            .head(5)[
                [
                    "cutoff_years",
                    "alpha",
                    "stride_bp",
                    "max_gap_strides",
                    "min_significant_strides",
                    "any_call",
                    "neutral_any_call",
                ]
            ]
            .to_string(index=False),
            flush=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
