"""Focal archaic-allele pair classes on the completed 400-haplotype study.

Classes describe arbitrary haplotype pairs, not diploid homozygous individuals.
Membership is fixed at the focal selected archaic allele while position varies.
This is a selected-cohort diagnostic, not a calibrated neutral scan.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import importlib.metadata
import inspect
import json
import os
from pathlib import Path
import shutil
import time
import traceback

os.environ["MPLBACKEND"] = "Agg"
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "1"

import numpy as np
import pandas as pd
import tskit
import zstandard

from .fresh_power import (
    atomic_json,
    canonical_hash,
    digest,
    ordered_nodes,
    read_profile,
)
from .posterior_replay import checked_source
from .recent_call_comparison import inventory

CLASSES = ("ref/ref", "alt/ref", "alt/alt", "all pairs")
COLORS = {
    "alt/alt": "#228833",
    "ref/ref": "#4477AA",
    "alt/ref": "#CC8833",
    "all pairs": "#555555",
}


def pair_classes(carriers, pairs):
    carriers, pairs = np.asarray(carriers), np.asarray(pairs)
    if carriers.ndim != 1 or not np.isin(carriers, (0, 1)).all():
        raise ValueError("Expected one binary archaic-allele state per haplotype")
    if (
        pairs.ndim != 2
        or pairs.shape[1] != 2
        or not np.issubdtype(pairs.dtype, np.integer)
    ):
        raise ValueError("Expected a two-column integer haplotype-pair manifest")
    if (
        np.any(pairs < 0)
        or np.any(pairs >= len(carriers))
        or np.any(pairs[:, 0] == pairs[:, 1])
    ):
        raise ValueError("Invalid haplotype indices")
    return carriers.astype(np.int8)[pairs].sum(axis=1)


def posterior_arrays(path, n_pairs, n_positions, chunk_size=8):
    expected = (
        ((n_pairs + chunk_size - 1) // chunk_size) * 2 * n_positions * chunk_size * 4
    )
    with (
        path.open("rb") as stream,
        zstandard.ZstdDecompressor().stream_reader(stream) as reader,
    ):
        raw = reader.read(expected + 1)
    if len(raw) != expected:
        raise ValueError("Truncated or oversized posterior stream")
    data = np.frombuffer(raw, dtype="<f4").reshape(-1, 2, n_positions, chunk_size)
    arrays = tuple(
        data[:, k].transpose(1, 0, 2).reshape(n_positions, -1)[:, :n_pairs]
        for k in (0, 1)
    )
    if any(not np.isfinite(a).all() or np.any(a <= 0) for a in arrays):
        raise ValueError("Invalid real-pair posterior parameters")
    return arrays


def summarize_decoded(alpha, beta, classes, cfg):
    two_ne = float(np.float32(cfg["decoder_scaled_mutation_rate"])) / (
        2 * cfg["mutation_rate"]
    )
    units = two_ne * cfg["generation_time_years"]
    ages = alpha.astype(np.float64) / beta.astype(np.float64) * units
    # The native mean-rule LUT is linear in alpha. Preserve its float32
    # arithmetic, endpoint clamping and >= tie convention exactly.
    cutoff_alpha = np.clip(alpha, np.float32(2**-16), np.float32(2**16))
    cutoffs = cfg["tmrca_cutoffs_years"]
    counts = np.zeros((4, alpha.shape[0], len(cutoffs)), dtype=np.int64)
    sums = np.zeros((4, alpha.shape[0]))
    for c in range(3):
        sums[c] = ages[:, classes == c].sum(axis=1)
    for j, years in enumerate(cutoffs):
        threshold = np.float32((years / cfg["generation_time_years"]) / two_ne)
        recent = beta * threshold >= cutoff_alpha
        for c in range(3):
            counts[c, :, j] = recent[:, classes == c].sum(axis=1)
    counts[3], sums[3] = counts[:3].sum(axis=0), sums[:3].sum(axis=0)
    return counts, sums


def summarize_truth(ts, pairs, classes, positions, cfg):
    nodes = ordered_nodes(ts)[pairs]
    cutoffs = np.asarray(cfg["tmrca_cutoffs_years"]) / cfg["generation_time_years"]
    counts = np.zeros((4, len(positions), len(cutoffs)), dtype=np.int64)
    sums = np.zeros((4, len(positions)))
    row = 0
    for tree in ts.trees():
        if row == len(positions):
            break
        if positions[row] >= tree.interval.right:
            continue
        times = np.fromiter((tree.tmrca(int(a), int(b)) for a, b in nodes), dtype=float)
        if not np.isfinite(times).all() or np.any(times < 0):
            raise ValueError("Nonfinite or negative true TMRCA")
        recent = times[:, None] < cutoffs
        local_counts = np.array([recent[classes == c].sum(axis=0) for c in range(3)])
        local_sums = (
            np.array([times[classes == c].sum() for c in range(3)])
            * cfg["generation_time_years"]
        )
        while row < len(positions) and positions[row] < tree.interval.right:
            counts[:3, row], sums[:3, row] = local_counts, local_sums
            row += 1
    if row != len(positions):
        raise ValueError("Tree sequence did not cover the output positions")
    counts[3], sums[3] = counts[:3].sum(axis=0), sums[:3].sum(axis=0)
    return counts, sums


def build_frame(counts, sums, classes, positions, cfg, source):
    denominators = [int(np.sum(classes == c)) for c in range(3)] + [len(classes)]
    frames = []
    for c, label in enumerate(CLASSES):
        n = denominators[c]
        block = pd.DataFrame(
            {
                "source": source,
                "pair_class": label,
                "position_0based": positions,
                "n_pairs": n,
                "mean_tmrca_years": sums[c] / n if n else np.nan,
            }
        )
        for j, cutoff in enumerate(cfg["tmrca_cutoffs_years"]):
            block[f"count_recent_{cutoff}"] = counts[c, :, j]
            block[f"frac_recent_{cutoff}"] = counts[c, :, j] / n if n else np.nan
        frames.append(block)
    return pd.concat(frames, ignore_index=True)


def profile_one(payload):
    item, comparison, out, computation_hash = payload
    cfg = item["cfg"]
    dest = Path(out) / "profiles" / item["key"]
    dest.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    try:
        raw_dir = Path(comparison) / "posteriors" / item["key"]
        posterior_receipt = json.loads(
            (raw_dir / "posterior_complete.json").read_text()
        )
        identity = dict(
            schema="archaic-pair-classes/v1",
            simulation=item["receipt_sha256"],
            posterior=digest(raw_dir / "posterior_complete.json"),
            pairs=digest(Path(item["root"]) / "pairs.tsv"),
            positions=digest(Path(item["root"]) / "positions.txt"),
            computation=computation_hash,
        )
        fingerprint = canonical_hash(identity)
        receipt = dest / "complete.json"
        profile = dest / "profile.csv.gz"
        if receipt.exists():
            old = json.loads(receipt.read_text())
            if old["fingerprint"] != fingerprint:
                raise ValueError("Existing profile has different inputs or computation")
            if (
                profile.stat().st_size != old["profile_bytes"]
                or digest(profile) != old["profile_sha256"]
            ):
                raise ValueError("Corrupt saved class profile")
            return dict(key=item["key"], status="cached")
        tree = checked_source(item, "decoded_input.trees")
        carriers = np.load(
            checked_source(item, "focal_carriers.npy"), allow_pickle=False
        )
        pairs = np.loadtxt(Path(item["root"]) / "pairs.tsv", dtype=int)
        positions = np.loadtxt(Path(item["root"]) / "positions.txt", dtype=int)
        ts = tskit.load(tree)
        focal = next(
            ts.variants(
                samples=ordered_nodes(ts),
                left=cfg["focal_position_bp"],
                right=cfg["focal_position_bp"] + 1,
            )
        )
        if focal.site.position != cfg["focal_position_bp"] or not np.array_equal(
            focal.genotypes != 0, carriers
        ):
            raise ValueError("Focal allele and saved haplotype mask disagree")
        if not np.isclose(carriers.mean(), item["sample_af"], rtol=0, atol=1e-12):
            raise ValueError("Carrier frequency differs from simulation receipt")
        classes = pair_classes(carriers, pairs)
        for name, spec in posterior_receipt["outputs"].items():
            p = raw_dir / name
            if p.stat().st_size != spec["bytes"] or digest(p) != spec["sha256"]:
                raise ValueError(f"Corrupt posterior input {name}")
        meta = json.loads((raw_dir / "posteriors.zst.meta").read_text())
        if (
            meta["num_pairs"] != len(pairs)
            or meta["chunk_size"] != 8
            or not np.array_equal(meta["pairs"], pairs)
            or not np.array_equal(meta["output_positions"], positions)
        ):
            raise ValueError(
                "Posterior ordering does not match pair/position manifests"
            )
        alpha, beta = posterior_arrays(
            raw_dir / "posteriors.zst", len(pairs), len(positions)
        )
        decoded = summarize_decoded(alpha, beta, classes, cfg)
        del alpha, beta
        truth = summarize_truth(ts, pairs, classes, positions, cfg)
        frames = []
        for source, values, old_name in (
            ("decoded", decoded, "frac_recent.tsv"),
            ("truth", truth, "truth_frac_recent.tsv"),
        ):
            counts, sums = values
            old = (
                read_profile(checked_source(item, old_name), cfg).iloc[:, 1:].to_numpy()
            )
            if not np.array_equal(
                counts[3], np.rint(old * len(pairs)).astype(np.int64)
            ):
                raise ValueError(
                    f"All-pair {source} counts differ from established baseline"
                )
            if not np.array_equal(counts[:3].sum(axis=0), counts[3]):
                raise ValueError("Pair classes do not partition the all-pair calls")
            frames.append(build_frame(counts, sums, classes, positions, cfg, source))
        frame = pd.concat(frames, ignore_index=True)
        frame["replicate"], frame["onset_years"], frame["sample_af"] = (
            item["replicate"],
            item["onset_years"],
            item["sample_af"],
        )
        temporary = profile.with_name("profile.tmp.gz")
        frame.to_csv(temporary, index=False, compression=dict(method="gzip", mtime=0))
        reread = pd.read_csv(temporary)
        if (
            len(reread) != 8 * len(positions)
            or reread.duplicated(["source", "pair_class", "position_0based"]).any()
        ):
            raise ValueError("Invalid written profile rows")
        temporary.replace(profile)
        atomic_json(
            receipt,
            dict(
                fingerprint=fingerprint,
                identity=identity,
                key=item["key"],
                onset_years=item["onset_years"],
                replicate=item["replicate"],
                sample_af=item["sample_af"],
                n_alt_haplotypes=int(carriers.sum()),
                n_ref_haplotypes=int((~carriers.astype(bool)).sum()),
                pair_counts={CLASSES[c]: int(np.sum(classes == c)) for c in range(3)},
                profile_sha256=digest(profile),
                profile_bytes=profile.stat().st_size,
                exact_all_pair_counts_reproduced=True,
                seconds=time.perf_counter() - started,
            ),
        )
        return dict(
            key=item["key"], status="complete", seconds=time.perf_counter() - started
        )
    except Exception:
        failure = dict(
            key=item["key"], status="failed", traceback=traceback.format_exc()
        )
        atomic_json(dest / "failed.json", failure)
        return failure


def analyse(items, out, cfg, workers):
    dest = out / "analysis"
    dest.mkdir(exist_ok=True)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        tables = list(
            pool.map(
                lambda i: pd.read_csv(out / "profiles" / i["key"] / "profile.csv.gz"),
                items,
            )
        )
    profiles = pd.concat(tables, ignore_index=True)
    del tables
    metrics = [
        "mean_tmrca_years",
        *(f"frac_recent_{t}" for t in cfg["tmrca_cutoffs_years"]),
    ]
    keys = ["onset_years", "source", "pair_class", "position_0based"]
    groups = profiles.groupby(keys, sort=True)
    mean, sd, n = groups[metrics].mean(), groups[metrics].std(), groups[metrics].count()
    from scipy.stats import t as student_t

    summary = mean.reset_index()[keys]
    for metric in metrics:
        error = (
            student_t.ppf(0.975, n[metric].to_numpy() - 1)
            * sd[metric].to_numpy()
            / np.sqrt(n[metric].to_numpy())
        )
        summary[metric] = mean[metric].to_numpy()
        summary[metric + "_n_regions"] = n[metric].to_numpy()
        summary[metric + "_lower95"] = np.maximum(0, mean[metric].to_numpy() - error)
        summary[metric + "_upper95"] = np.minimum(
            1 if metric.startswith("frac") else np.inf, mean[metric].to_numpy() + error
        )
    summary.to_csv(dest / "mean_profiles.csv", index=False, float_format="%.9g")
    focal = profiles[profiles.position_0based == cfg["focal_position_bp"]].copy()
    focal.to_csv(dest / "focal_by_replicate.csv", index=False)
    summary[summary.position_0based == cfg["focal_position_bp"]].to_csv(
        dest / "focal_summary.csv", index=False, float_format="%.9g"
    )
    receipts = [
        json.loads((out / "profiles" / i["key"] / "complete.json").read_text())
        for i in items
    ]
    pd.DataFrame(
        [
            dict(
                onset_years=r["onset_years"],
                replicate=r["replicate"],
                sample_af=r["sample_af"],
                n_alt_haplotypes=r["n_alt_haplotypes"],
                n_ref_haplotypes=r["n_ref_haplotypes"],
                **r["pair_counts"],
            )
            for r in receipts
        ]
    ).to_csv(dest / "pair_counts.csv", index=False)
    figures(summary, dest, cfg)
    atomic_json(
        dest / "provenance.json",
        dict(
            statistic="frac_recent_T",
            hard_call="Native posterior-mean rule; exact all-pair baseline agreement",
            pair_definition="Both haplotypes' alleles at the known focal selected archaic SNP; membership held fixed across coordinates",
            populations="EAS selected s=0.005, onset 50k and 10k; no neutral inference",
            region_weighting="Each nonempty class contributes one equally weighted region; empty classes retained as missing",
            uncertainty="Pointwise Student-t 95% intervals across independent simulated regions; pairs are not independent observations",
            input_profiles=receipts,
            source_sha256=digest(Path(__file__)),
            versions={
                x: importlib.metadata.version(x)
                for x in (
                    "numpy",
                    "pandas",
                    "tskit",
                    "zstandard",
                    "scipy",
                    "matplotlib",
                )
            },
            outputs={
                p.name: dict(bytes=p.stat().st_size, sha256=digest(p))
                for p in dest.iterdir()
                if p.is_file() and p.name != "provenance.json"
            },
        ),
    )


def figures(summary, out, cfg):
    import matplotlib
    import matplotlib.pyplot as plt

    matplotlib.rcParams.update(
        {
            "font.size": 17,
            "axes.titlesize": 17,
            "axes.labelsize": 19,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    checks = []
    for metric in (
        "mean_tmrca_years",
        *(f"frac_recent_{t}" for t in cfg["tmrca_cutoffs_years"]),
    ):
        age = metric == "mean_tmrca_years"
        for scope, extent in (("full", 5), ("zoom", 0.5)):
            fig, axes = plt.subplots(2, 2, figsize=(11, 8.5), sharex=True, sharey=True)
            fig.subplots_adjust(
                left=0.12, right=0.98, bottom=0.14, top=0.79, wspace=0.14, hspace=0.28
            )
            for row, source in enumerate(("truth", "decoded")):
                for column, onset in enumerate((50000, 10000)):
                    ax = axes[row, column]
                    for label in ("alt/alt", "ref/ref", "alt/ref", "all pairs"):
                        d = summary[
                            (summary.source == source)
                            & (summary.onset_years == onset)
                            & (summary.pair_class == label)
                        ].sort_values("position_0based")
                        x = (
                            d.position_0based.to_numpy() - cfg["focal_position_bp"]
                        ) / 1e6
                        use = abs(x) <= extent
                        scale = 1000 if age else 1
                        ax.plot(
                            x[use],
                            d[metric].to_numpy()[use] / scale,
                            color=COLORS[label],
                            linestyle="--" if label == "all pairs" else "-",
                            linewidth=2.3 if label in ("alt/alt", "ref/ref") else 1.6,
                            label=label,
                        )
                        ax.fill_between(
                            x[use],
                            d[metric + "_lower95"].to_numpy()[use] / scale,
                            d[metric + "_upper95"].to_numpy()[use] / scale,
                            color=COLORS[label],
                            alpha=0.10,
                            linewidth=0,
                        )
                    ax.axvline(0, color="gray", linestyle=":", linewidth=1)
                    ax.set_title(
                        f"{'True' if source == 'truth' else 'Decoded'} TMRCA · onset {onset // 1000} kya"
                    )
                    ax.set_xlim(-extent, extent)
                    ax.set_xticks(
                        [-extent, 0, extent] if extent < 1 else [-5, -2.5, 0, 2.5, 5]
                    )
                    ax.spines[["top", "right"]].set_visible(False)
                    if age:
                        ax.set_yscale("log")
                    else:
                        ax.set_ylim(0, 1)
                        ax.set_yticks([0, 0.5, 1])
            if age:
                positive = summary[metric + "_lower95"].to_numpy() / 1000
                lo = 10 ** np.floor(np.log10(positive[positive > 0].min()))
                hi = 10 ** np.ceil(np.log10(summary[metric + "_upper95"].max() / 1000))
                ticks = 10.0 ** np.arange(int(np.log10(lo)), int(np.log10(hi)) + 1)
                axes[0, 0].set_ylim(lo, hi)
                axes[0, 0].set_yticks(ticks, [f"{v:g}" for v in ticks])
                title, ylabel = (
                    "Pairwise ages by focal archaic allele",
                    "Mean pairwise TMRCA (kya)",
                )
            else:
                title, ylabel = (
                    f"Recent pairs: T = {int(metric.split('_')[-1]) // 1000} kya",
                    "Fraction of pairs recent within class",
                )
            suptitle = fig.suptitle(title, y=0.99, fontsize=21)
            xlabel = fig.supxlabel(
                "Distance from the archaic SNP (Mb)", y=0.035, fontsize=19
            )
            ylabel_obj = fig.supylabel(ylabel, x=0.02, fontsize=19)
            legend = fig.legend(
                *axes[0, 0].get_legend_handles_labels(),
                loc="upper center",
                bbox_to_anchor=(0.55, 0.94),
                ncol=4,
                frameon=False,
                fontsize=16,
            )
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
            assert all(
                not legend.get_window_extent(renderer).overlaps(
                    ax.get_window_extent(renderer)
                )
                for ax in axes.flat
            )
            text_items = [suptitle, xlabel, ylabel_obj, *legend.get_texts()]
            for ax in axes.flat:
                text_items.extend(
                    [ax.title, *ax.get_xticklabels(), *ax.get_yticklabels()]
                )
            for text in text_items:
                if not text.get_visible() or not text.get_text():
                    continue
                box = text.get_window_extent(renderer)
                assert (
                    box.x0 >= 0
                    and box.y0 >= 0
                    and box.x1 <= fig.bbox.width
                    and box.y1 <= fig.bbox.height
                ), text.get_text()
            stem = f"{metric}_{scope}"
            fig.savefig(out / f"{stem}.png", dpi=200)
            fig.savefig(
                out / f"{stem}.pdf", metadata={"CreationDate": None, "ModDate": None}
            )
            assert b"/Subtype /Image" not in (out / f"{stem}.pdf").read_bytes()
            checks.append(
                dict(
                    figure=stem,
                    text_within_canvas=True,
                    legend_outside_panels=True,
                    vector_pdf=True,
                    displayed_inline=False,
                )
            )
            plt.close(fig)
    atomic_json(out / "figure_checks.json", checks)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--later", type=Path, required=True)
    parser.add_argument("--comparison", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--decoder", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=20)
    parser.add_argument("--limit-per-onset", type=int)
    parser.add_argument(
        "--phase", choices=("profiles", "analyse", "run"), default="run"
    )
    args = parser.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    import fcntl

    with (args.out / "runner.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        items, cfg, studies = inventory(args.baseline, args.later, args.decoder)
        items = [
            i
            for i in items
            if i["mode"] == "selected"
            and (args.limit_per_onset is None or i["replicate"] < args.limit_per_onset)
        ]
        if not items or args.workers < 1:
            raise ValueError("No selected regions or invalid worker count")
        code_hash = canonical_hash(
            [
                inspect.getsource(f)
                for f in (
                    pair_classes,
                    posterior_arrays,
                    summarize_decoded,
                    summarize_truth,
                    build_frame,
                )
            ]
        )
        atomic_json(
            args.out / "manifest.json",
            dict(
                studies=studies,
                config=cfg,
                selected_regions=len(items),
                sample_ids=[i["key"] for i in items],
                computation_hash=code_hash,
                source_sha256=digest(Path(__file__)),
                workers=args.workers,
                new_rng="None; original simulation/pair seeds retained",
            ),
        )
        if args.phase != "analyse":
            if (
                shutil.disk_usage(args.out).used + 20_000_000_000
                > cfg["storage_limit_bytes"]
            ):
                raise ValueError("Storage reserve exceeded")
            completed, failures = 0, []
            with ProcessPoolExecutor(max_workers=args.workers) as pool:
                futures = [
                    pool.submit(
                        profile_one, (i, str(args.comparison), str(args.out), code_hash)
                    )
                    for i in items
                ]
                for future in as_completed(futures):
                    result = future.result()
                    if result["status"] == "failed":
                        failures.append(result)
                    else:
                        completed += 1
                    print(json.dumps(result), flush=True)
                    atomic_json(
                        args.out / "status.json",
                        dict(
                            completed=completed,
                            target=len(items),
                            failed=failures,
                            state="running",
                        ),
                    )
            atomic_json(
                args.out / "status.json",
                dict(
                    completed=completed,
                    target=len(items),
                    failed=failures,
                    state="failed" if failures else "complete",
                ),
            )
            if failures:
                raise RuntimeError("Class profiling failed; inspect status.json")
        if args.phase != "profiles":
            for item in items:
                receipt = json.loads(
                    (args.out / "profiles" / item["key"] / "complete.json").read_text()
                )
                if (
                    receipt["identity"]["computation"] != code_hash
                    or receipt["identity"]["simulation"] != item["receipt_sha256"]
                    or receipt["identity"]["posterior"]
                    != digest(
                        args.comparison
                        / "posteriors"
                        / item["key"]
                        / "posterior_complete.json"
                    )
                    or digest(args.out / "profiles" / item["key"] / "profile.csv.gz")
                    != receipt["profile_sha256"]
                ):
                    raise ValueError("Class profile hash or computation differs")
            analyse(items, args.out, cfg, args.workers)
        print("Pair-class phase complete", flush=True)


if __name__ == "__main__":
    main()
