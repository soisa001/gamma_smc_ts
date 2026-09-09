"""Matched hard-call decoding and 10/50 kya selection-onset comparison.

Simulation inputs are immutable. Every hard-call rule has its own matched null.
No posterior probabilities are averaged into the test statistic.
"""

from __future__ import annotations

import argparse
from concurrent.futures import (
    ProcessPoolExecutor,
    ThreadPoolExecutor,
    wait,
    FIRST_COMPLETED,
)
import json
import os
from pathlib import Path
import shutil
import time

os.environ["MPLBACKEND"] = "Agg"
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "1"

import numpy as np
import pandas as pd

from .fresh_power import atomic_json, digest, read_profile
from .fresh_region_metrics import aligned_pvalues
from .panel_size_comparison import evaluate

from .posterior_replay import (
    RULES,
    cached_decode,
    checked_source,
    decode_identity,
    decode_bundle,
)

SOURCES = ["mean", "median", "prob80", "prob90", "truth"]


def inventory(baseline, later, decoder, allow_incomplete=False):
    items, configs, studies = [], [], []
    replay_helper_sha = digest(decoder.with_name("summarize_recent_rules"))
    for root in (baseline, later):
        study = json.loads((root / "manifest.json").read_text())
        cfg = study["config"]
        configs.append(cfg)
        if study["scientific_identity"]["decoder_sha256"] != digest(decoder):
            raise ValueError("Decoder binary differs from the original mean decoder")
        for name, field in (
            ("pairs.tsv", "pairs_sha256"),
            ("positions.txt", "positions_sha256"),
        ):
            if digest(root / name) != study[field]:
                raise ValueError(f"Corrupt {name} in {root}")
        studies.append(
            dict(
                root=str(root),
                manifest_sha256=digest(root / "manifest.json"),
                sample_manifest_sha256=digest(root / "sample_manifest.json"),
                fingerprint=study["fingerprint"],
            )
        )
        manifest = json.loads((root / "sample_manifest.json").read_text())
        for entry in manifest:
            if entry["mode"] == "selected" and entry["s"] != 0.005:
                continue
            directory = root / entry["task_id"]
            path = directory / "complete.json"
            if allow_incomplete and not path.exists():
                continue
            record = json.loads(path.read_text())
            if record["fingerprint"] != study["fingerprint"]:
                raise ValueError(f"Mismatched simulation receipt: {directory}")
            if any(
                record["task"][key] != entry[key] for key in ("mode", "s", "replicate")
            ):
                raise ValueError(f"Mismatched task identity: {directory}")
            onset = 0 if entry["mode"] == "neutral" else cfg["selection_onset_years"]
            key = (
                f"neutral/rep{entry['replicate']:04d}"
                if onset == 0
                else f"onset{onset}/rep{entry['replicate']:04d}"
            )
            items.append(
                dict(
                    key=key,
                    directory=str(directory),
                    root=str(root),
                    cfg=cfg,
                    task_id=key,
                    mode=entry["mode"],
                    s=entry["s"],
                    replicate=entry["replicate"],
                    onset_years=onset,
                    sample_af=record["sample_af"],
                    fixed=record["fixed"],
                    artifacts=record["artifacts"],
                    receipt_sha256=digest(path),
                    replay_helper_sha256=replay_helper_sha,
                )
            )
    allowed = {"selection_onset_years", "neutral_replicates"}
    if {k: v for k, v in configs[0].items() if k not in allowed} != {
        k: v for k, v in configs[1].items() if k not in allowed
    }:
        raise ValueError("Scientific settings differ beyond selection onset")
    if (
        configs[0]["selection_onset_years"] != 50000
        or configs[1]["selection_onset_years"] != 10000
    ):
        raise ValueError("Expected selection onsets 50 and 10 kya")
    if (
        configs[1]["neutral_replicates"] != 0
        or configs[0]["neutral_replicates"] != 1000
    ):
        raise ValueError("Expected one shared 1000-neutral panel")
    if len({i["key"] for i in items}) != len(items):
        raise ValueError("Duplicate simulation identities")
    if not allow_incomplete and len(items) != 1200:
        raise ValueError("Expected 1000 neutral and 100 selected per onset")
    if digest(baseline / "pairs.tsv") != digest(later / "pairs.tsv"):
        raise ValueError("Pair panels differ")
    if max(configs[0]["tmrca_cutoffs_years"]) > 50000:
        raise ValueError("Cutoffs must not exceed 50 kya")
    return items, configs[0], studies


def decode(items, cfg, out, decoder, workers, smoke=False):
    if smoke:
        items = [i for i in items if i["replicate"] == 0]
    binary_sha = digest(decoder)
    tasks = [(i, str(out), str(decoder), binary_sha) for i in items]
    done, pending, failures, cursor = 0, {}, [], 0

    def status(state):
        result = dict(
            state=state,
            completed=done,
            target=len(tasks),
            workers=workers,
            failed=failures,
            updated_at=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        )
        atomic_json(out / "decode_status.json", result)
        print(json.dumps(result), flush=True)

    with ProcessPoolExecutor(max_workers=min(workers, len(tasks))) as pool:
        while cursor < len(tasks) or pending:
            while cursor < len(tasks) and len(pending) < workers and not failures:
                usage = shutil.disk_usage(out)
                if usage.used + 20_000_000_000 > cfg["storage_limit_bytes"]:
                    failures.append(dict(error="Storage reserve reached"))
                    break
                pending[pool.submit(decode_bundle, tasks[cursor])] = cursor
                cursor += 1
            if not pending:
                break
            finished, _ = wait(pending, timeout=30, return_when=FIRST_COMPLETED)
            for future in finished:
                pending.pop(future)
                result = future.result()
                if result["status"] == "failed":
                    failures.append(result)
                else:
                    done += 1
                print(json.dumps(result), flush=True)
            status("running")
    status("failed" if failures else "complete")
    if failures or done != len(tasks):
        raise RuntimeError("Decoding incomplete; inspect decode_status.json")
    if smoke:
        checks = []
        for item in items:
            profiles = {
                rule: read_profile(
                    out / "decoded" / rule / item["key"] / "frac_recent.tsv", cfg
                )
                .iloc[:, 1:]
                .to_numpy()
                for rule in RULES
            }
            mean = (
                read_profile(checked_source(item, "frac_recent.tsv"), cfg)
                .iloc[:, 1:]
                .to_numpy()
            )
            for lower, upper in (
                (profiles["prob90"], profiles["prob80"]),
                (profiles["prob80"], profiles["median"]),
                (mean, profiles["median"]),
            ):
                if np.any(lower > upper + 1e-12):
                    raise ValueError("Smoke hard-call nesting failed")
            checks.append(item["key"])
        atomic_json(
            out / "smoke_validation.json",
            dict(
                status="passed",
                samples=checks,
                checks=[
                    "mean <= median",
                    "prob90 <= prob80 <= median",
                    "all pairs retained",
                ],
            ),
        )


def longest_run(mask, gap):
    """Number of significant strides in the largest run with <=gap misses."""
    count = np.zeros(mask.shape[0], dtype=int)
    best = count.copy()
    missing = count.copy()
    for column in mask.T:
        missing = np.where(column, 0, missing + 1)
        count = np.where(missing > gap, 0, count + column)
        best = np.maximum(best, count)
    return best


def grid_summary(values, meta, cfg, source, onset):
    null = meta["mode"].to_numpy() == "neutral"
    pn, ps = aligned_pvalues(
        np.rint(values[null] * cfg["haplotype_pairs"]).astype(int),
        np.rint(values[~null] * cfg["haplotype_pairs"]).astype(int),
    )
    p = np.empty_like(values)
    p[null], p[~null] = pn, ps
    rows = []
    for t, cutoff in enumerate(cfg["tmrca_cutoffs_years"]):
        for alpha in (0.05, 0.01, 0.001):
            for step in (1, 2, 5):
                for gap in (0, 1):
                    lengths = longest_run(p[:, ::step, t] <= alpha, gap)
                    for minimum in (1, 3, 5, 10, 20):
                        for label, subset in (("neutral", null), ("selected", ~null)):
                            called = lengths[subset] >= minimum
                            rows.append(
                                dict(
                                    source=source,
                                    onset_years=onset,
                                    cutoff_years=cutoff,
                                    alpha=alpha,
                                    stride_bp=cfg["stride_bp"] * step,
                                    max_gap_strides=gap,
                                    min_significant_strides=minimum,
                                    mode=label,
                                    n=int(subset.sum()),
                                    called=int(called.sum()),
                                    region_call_rate=float(called.mean()),
                                )
                            )
    return rows


def load_rule(items, source, out, decoder_sha, workers):
    def read(item):
        if source in ("mean", "truth"):
            path = checked_source(
                item, "frac_recent.tsv" if source == "mean" else "truth_frac_recent.tsv"
            )
        else:
            directory = out / "decoded" / source / item["key"]
            if not cached_decode(
                directory, decode_identity(item, source, decoder_sha), item["cfg"]
            ):
                raise ValueError(f"Missing decoding: {directory}")
            path = directory / "frac_recent.tsv"
        values = read_profile(path, item["cfg"]).iloc[:, 1:].to_numpy()
        # Recover the exact integer count before dividing, as in the baseline
        # analyses, so CSV parser rounding cannot alter tied score ranks.
        values = (
            np.rint(values * item["cfg"]["haplotype_pairs"])
            / item["cfg"]["haplotype_pairs"]
        )
        return values, dict(key=item["key"], source=source, sha256=digest(path))

    with ThreadPoolExecutor(max_workers=workers) as pool:
        loaded = list(pool.map(read, items))
    return np.stack([row[0] for row in loaded]), [row[1] for row in loaded]


def analyse(items, cfg, out, decoder, workers):
    dest = out / "analysis"
    dest.mkdir(exist_ok=True)
    meta = pd.DataFrame(
        [
            {
                k: i[k]
                for k in (
                    "task_id",
                    "mode",
                    "s",
                    "replicate",
                    "onset_years",
                    "sample_af",
                    "fixed",
                )
            }
            for i in items
        ]
    )
    meta.to_csv(dest / "replicates.csv", index=False)
    truth, inputs = load_rule(items, "truth", out, digest(decoder), workers)
    all_rows, choices, grid, diagnostics, prior = [], [], [], [], None
    for source in SOURCES:
        values, provenance = load_rule(items, source, out, digest(decoder), workers)
        inputs.extend(provenance)
        if source in ("prob80", "prob90") and np.any(values > prior + 1e-12):
            raise ValueError("Posterior mass calls not nested by required probability")
        if source in ("median", "prob80"):
            prior = values.copy()
        for onset in (50000, 10000):
            use = (meta["mode"] == "neutral") | (meta.onset_years == onset)
            m, v = meta[use].reset_index(drop=True), values[use]
            rows, trained = evaluate(
                v, m, cfg["tmrca_cutoffs_years"], cfg["seed_base"], 400, source
            )
            for block in (rows, trained):
                for row in block:
                    row["onset_years"] = onset
            all_rows.extend(rows)
            choices.extend(trained)
            grid.extend(grid_summary(v, m, cfg, source, onset))
        for onset in (0, 10000, 50000):
            use = meta.onset_years.to_numpy() == onset
            for j, cutoff in enumerate(cfg["tmrca_cutoffs_years"]):
                value = values[use, :, j]
                diagnostics.append(
                    dict(
                        source=source,
                        onset_years=onset,
                        cutoff_years=cutoff,
                        n=int(use.sum()),
                        region_mean=float(value.mean()),
                        average_position_sd=float(value.std(axis=0, ddof=1).mean()),
                        focal_mean=float(value[:, 500].mean()),
                        rmse_to_truth=float(
                            np.sqrt(
                                ((value - truth[use, :, j]) ** 2).mean(axis=1)
                            ).mean()
                        ),
                    )
                )
        print(f"Analysed {source}", flush=True)
    raw = pd.DataFrame(all_rows)
    raw.to_csv(
        dest / "held_out_regions.csv.gz",
        index=False,
        compression=dict(method="gzip", mtime=0),
    )
    keys = [
        "onset_years",
        "source",
        "method",
        "cutoff_years",
        "policy",
        "target",
        "mode",
    ]
    summary = (
        raw.groupby(keys)
        .agg(
            n=("called", "size"),
            called=("called", "sum"),
            region_call_rate=("called", "mean"),
        )
        .reset_index()
    )
    summary.to_csv(dest / "region_calls.csv", index=False)
    pd.DataFrame(choices).to_csv(dest / "training_thresholds.csv", index=False)
    pd.DataFrame(grid).to_csv(dest / "paired_grid.csv", index=False)
    pd.DataFrame(diagnostics).to_csv(dest / "diagnostics.csv", index=False)
    figures(summary, dest)
    atomic_json(
        dest / "provenance.json",
        dict(
            inputs=inputs,
            seed=cfg["seed_base"],
            rules=RULES,
            statistic="frac_recent_T",
            config=cfg,
            neutral_panel_shared=True,
            scope="Internal evaluation; per-rule results, no winner independently confirmed",
            primary_endpoint="Any call anywhere in each 10 Mb region; all selected replicates retained",
            evaluation="Five whole-replicate folds; 400 fit neutrals, 400 calibration neutrals, 200 test neutrals; 80 train/20 test selected",
            source_hashes={
                str(p.name): digest(p)
                for p in (
                    Path(__file__),
                    Path(__file__).with_name("panel_size_comparison.py"),
                    Path(__file__).with_name("partial_sweep_cv.py"),
                    Path(__file__).with_name("posterior_replay.py"),
                    Path(__file__).resolve().parents[2]
                    / "src/summarize_recent_rules.cpp",
                    Path(__file__).resolve().parents[2] / "src/recent_stats.h",
                )
            },
            decoder_sha256=digest(decoder),
            replay_helper_sha256=digest(decoder.with_name("summarize_recent_rules")),
            outputs={
                p.name: dict(bytes=p.stat().st_size, sha256=digest(p))
                for p in dest.iterdir()
                if p.is_file() and p.name != "provenance.json"
            },
        ),
    )


def figures(summary, out):
    import matplotlib
    import matplotlib.pyplot as plt

    matplotlib.rcParams.update(
        {"font.size": 17, "axes.labelsize": 19, "pdf.fonttype": 42, "ps.fonttype": 42}
    )
    colors = dict(zip(SOURCES, ("#377eb8", "#e58224", "#39945b", "#964da5", "#333333")))
    layout_checks = []
    for alpha in (0.05, 0.1):
        fig, axes = plt.subplots(1, 2, figsize=(11, 8.5), sharey=True)
        fig.subplots_adjust(left=0.11, right=0.98, bottom=0.15, top=0.76, wspace=0.13)
        for ax, onset in zip(axes, (50000, 10000)):
            for source in SOURCES:
                d = summary[
                    (summary.source == source)
                    & (summary.onset_years == onset)
                    & (summary.method == "maximum_frac_recent")
                    & (summary.policy == "neutral_target")
                    & (summary.target == alpha)
                    & (summary["mode"] == "selected")
                ].sort_values("cutoff_years")
                ax.plot(
                    d.cutoff_years / 1000,
                    d.region_call_rate,
                    "o-",
                    color=colors[source],
                    label=source,
                    linewidth=2,
                )
            ax.set(
                title=f"Selection onset {onset // 1000} kya",
                xlabel="TMRCA cutoff (kya)",
                ylim=(-0.02, 1.02),
                xticks=[5, 10, 20, 30, 40, 50],
                yticks=np.linspace(0, 1, 6),
            )
            ax.axhline(0.7, color="gray", linestyle=":")
            ax.spines[["top", "right"]].set_visible(False)
        axes[0].set_ylabel("Fraction of selected regions called")
        legend = fig.legend(
            *axes[0].get_legend_handles_labels(),
            loc="upper center",
            bbox_to_anchor=(0.5, 0.90),
            ncol=5,
            frameon=False,
            fontsize=14,
        )
        fig.suptitle(
            f"Regional maximum; nominal neutral error {alpha:.0%}", y=0.99, fontsize=20
        )
        # Inspect layout bounds on the headless canvas without emitting images.
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        legend_box = legend.get_window_extent(renderer)
        assert all(
            not legend_box.overlaps(ax.get_window_extent(renderer)) for ax in axes
        )
        texts = [fig._suptitle, *legend.get_texts()]
        for ax in axes:
            texts.extend(
                [
                    ax.title,
                    ax.xaxis.label,
                    ax.yaxis.label,
                    *ax.get_xticklabels(),
                    *ax.get_yticklabels(),
                ]
            )
        for text in texts:
            if not text.get_visible() or not text.get_text():
                continue
            box = text.get_window_extent(renderer)
            assert (
                box.x0 >= 0
                and box.y0 >= 0
                and box.x1 <= fig.bbox.width
                and box.y1 <= fig.bbox.height
            ), text.get_text()
        fig.savefig(out / f"power_neutral_{alpha:g}.png", dpi=200)
        pdf = out / f"power_neutral_{alpha:g}.pdf"
        fig.savefig(pdf)
        assert b"/Subtype /Image" not in pdf.read_bytes(), (
            "Unexpected rasterized PDF content"
        )
        layout_checks.append(
            dict(
                alpha=alpha,
                all_text_in_canvas=True,
                legend_outside_panels=True,
                pdf_contains_no_images=True,
                pdf_fonttype=42,
                displayed_inline=False,
            )
        )
        plt.close(fig)
    atomic_json(out / "figure_layout_checks.json", layout_checks)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--later", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--decoder", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=20)
    parser.add_argument(
        "--phase", choices=("smoke", "decode", "analyse", "run"), default="run"
    )
    args = parser.parse_args(argv)
    if not 1 <= args.workers <= 20:
        raise ValueError("Use 1-20 workers")
    args.out.mkdir(parents=True, exist_ok=True)
    import fcntl

    with (args.out / "runner.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        items, cfg, studies = inventory(
            args.baseline, args.later, args.decoder, args.phase == "smoke"
        )
        atomic_json(
            args.out / "manifest.json",
            dict(
                studies=studies,
                rules=RULES,
                config=cfg,
                workers=args.workers,
                selected_onsets=[50000, 10000],
                items=items,
            ),
        )
        if args.phase in ("smoke", "decode", "run"):
            decode(
                items, cfg, args.out, args.decoder, args.workers, args.phase == "smoke"
            )
        if args.phase in ("analyse", "run"):
            analyse(items, cfg, args.out, args.decoder, args.workers)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
