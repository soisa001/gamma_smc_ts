"""Package synthetic benchmark inputs, summaries, provenance, and static plots."""

import argparse
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

os.environ["MPLBACKEND"] = "Agg"
import matplotlib

matplotlib.rcParams.update(
    {
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.size": 15,
        "axes.titlesize": 17,
        "axes.labelsize": 16,
        "xtick.labelsize": 13,
        "ytick.labelsize": 13,
        "legend.fontsize": 14,
    }
)
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def write_json(p, obj):
    p.write_text(json.dumps(obj, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source", type=Path, default=ROOT.parent / "sim/sv_phylogeny_20261028"
    )
    parser.add_argument(
        "--out", type=Path, default=ROOT / "docs/results/sv_phylogeny_20261028"
    )
    args = parser.parse_args()
    source, out = args.source, args.out
    out.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(source / "estimates.tsv", sep="\t")
    summary = pd.read_csv(source / "summary.tsv", sep="\t")
    assert len(summary) == 8 and (summary.events == 12).all()
    assert json.loads((source / "audit.json").read_text())["passed"]
    for name in (
        "estimates.tsv",
        "summary.tsv",
        "beast_diagnostics.tsv",
        "windows.tsv",
        "by_history.tsv",
        "panel.json",
        "audit.json",
    ):
        if (source / name).exists():
            shutil.copyfile(source / name, out / name)
            assert digest(source / name) == digest(out / name)
    colors = {
        "post_introgression": "#0072B2",
        "imported_old": "#D55E00",
        "ancestral_retained": "#009E73",
    }
    labels = {
        "post_introgression": "Post-introgression",
        "imported_old": "Imported",
        "ancestral_retained": "Ancestral retained",
    }
    methods = [
        ("IQTREE_LSD2", "w10000", "IQ-TREE + LSD2"),
        ("BEAST2", "w10000", "BEAST 2"),
        ("Relate", "500kb", "Relate"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(11, 8.5), sharex=True, sharey=True)
    for ax, (method, window, label) in zip(axes, methods):
        subset = frame[
            (frame.method == method) & (frame.window == window) & frame.usable
        ]
        for history, group in subset.groupby("history"):
            x = group.true_age_years.to_numpy() / 1000
            y = group.estimate_years.to_numpy() / 1000
            lo = np.maximum(group.lower_years.to_numpy() / 1000, 0.1)
            hi = group.upper_years.to_numpy() / 1000
            ax.errorbar(
                x,
                y,
                yerr=np.maximum(np.vstack([y - lo, hi - y]), 0),
                fmt="o",
                markersize=8,
                capsize=3,
                color=colors[history],
                alpha=0.85,
            )
        ax.plot([0.1, 10000], [0.1, 10000], color="0.5", ls=":", lw=1.5)
        ax.axvline(50, color="0.65", ls="--", lw=1)
        ax.axhline(50, color="0.65", ls="--", lw=1)
        ax.set(
            xscale="log",
            yscale="log",
            xlim=(1, 6000),
            ylim=(0.5, 6000),
            xlabel="True age (kyr)",
            title=f"{label}\n{len(subset)}/12 placements",
        )
        ax.set_xticks([1, 10, 100, 1000], ["1", "10", "100", "1,000"])
        ax.grid(alpha=0.15)
    axes[0].set_ylabel("Estimated age (kyr)")
    handles = [
        Line2D(
            [0], [0], marker="o", color=c, linestyle="", label=labels[h], markersize=8
        )
        for h, c in colors.items()
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        ncol=3,
        bbox_to_anchor=(0.5, 0.995),
        frameon=False,
    )
    fig.subplots_adjust(left=0.10, right=0.98, bottom=0.14, top=0.84, wspace=0.12)
    for ext in ("png", "pdf"):
        fig.savefig(out / f"age_comparison.{ext}", dpi=180)
    plt.close(fig)
    # Show every event, including failures, so conditional age accuracy cannot
    # hide loci for which a method has no supported placement.
    panel = json.loads((source / "panel.json").read_text())
    ordered = sorted(panel, key=lambda x: (list(colors).index(x["history"]), x["id"]))
    columns = [
        ("IQTREE_LSD2", "oracle", "IQ\noracle"),
        ("IQTREE_LSD2", "w1000", "IQ\n1 kb"),
        ("IQTREE_LSD2", "w10000", "IQ\n10 kb"),
        ("IQTREE_LSD2", "w25000", "IQ\n25 kb"),
        ("BEAST2", "oracle", "BEAST\noracle"),
        ("BEAST2", "w10000", "BEAST\n10 kb"),
        ("Relate", "500kb", "Relate\n500 kb"),
    ]
    matrix = np.zeros((len(ordered), len(columns)))
    for i, event in enumerate(ordered):
        for j, (method, window, _) in enumerate(columns):
            row = frame[
                (frame.id == event["id"])
                & (frame.method == method)
                & (frame.window == window)
            ].iloc[0]
            matrix[i, j] = 0 if not row.usable else (2 if row.age_contained else 1)
    from matplotlib.colors import ListedColormap

    fig, ax = plt.subplots(figsize=(11, 8.5))
    ax.pcolormesh(
        np.arange(len(columns) + 1) - 0.5,
        np.arange(len(ordered) + 1) - 0.5,
        matrix,
        cmap=ListedColormap(["#eeeeee", "#E69F00", "#56B4E9"]),
        vmin=-0.5,
        vmax=2.5,
        shading="flat",
        rasterized=False,
    )
    ax.invert_yaxis()
    ax.set_xticks(range(len(columns)), [c[2] for c in columns])
    ax.set_yticks(
        range(len(ordered)),
        [
            f"{e['id']}  { {'post_introgression': 'New', 'imported_old': 'Imported', 'ancestral_retained': 'Ancestral'}[e['history']] }"
            for e in ordered
        ],
    )
    for i in range(len(ordered)):
        for j in range(len(columns)):
            ax.text(
                j,
                i,
                ["—", "Outside", "Inside"][int(matrix[i, j])],
                ha="center",
                va="center",
                fontsize=12,
            )
    ax.set_xticks(np.arange(-0.5, len(columns), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(ordered), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=2)
    ax.tick_params(which="minor", bottom=False, left=False)
    handles = [
        Line2D([0], [0], marker="s", markersize=12, linestyle="", color=c, label=l)
        for c, l in zip(
            ["#eeeeee", "#E69F00", "#56B4E9"],
            ["No usable placement", "Age outside bounds", "Age inside bounds"],
        )
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.5, 1.0),
    )
    fig.subplots_adjust(left=0.29, right=0.98, bottom=0.11, top=0.91)
    for ext in ("png", "pdf"):
        fig.savefig(out / f"all_locus_outcomes.{ext}", dpi=180)
    plt.close(fig)
    files = {"panel.json": source / "panel.json"}
    for event in panel:
        d = source / event["id"]
        for name in [
            "input.json",
            "evaluation.json",
            "human.haps",
            "human.sample",
            "genetic.map",
            "rate.txt",
        ] + [p.name for p in d.glob("*.fasta")]:
            files[event["id"] + "/" + name] = d / name
    archive = out / "benchmark_inputs.zip"
    with zipfile.ZipFile(
        archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9
    ) as z:
        for name, p in sorted(files.items()):
            info = zipfile.ZipInfo(name, date_time=(2026, 10, 28, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            z.writestr(info, p.read_bytes())
    write_json(
        out / "inputs_manifest.json",
        {
            "synthetic_only": True,
            "zip_sha256": digest(archive),
            "files": {name: digest(p) for name, p in files.items()},
        },
    )
    tools = ROOT / ".tools/phylo"
    provenance = {
        "source_directory": str(source),
        "all_inputs_synthetic": True,
        "config": json.loads((ROOT / "configs/sv_phylogeny.json").read_text()),
        "raw_receipts": {},
        "raw_commands": {},
    }
    provenance["runtime"] = {
        "python": sys.version,
        "platform": platform.platform(),
        "java": subprocess.run(
            ["java", "-version"], capture_output=True, text=True
        ).stderr,
        "python_packages": subprocess.check_output(
            ["uv", "pip", "freeze", "--python", sys.executable], text=True
        ),
        "iqtree_version": "2.4.0",
        "beast_base_and_app_version": "2.7.8",
        "beast_launcher_version": "2.7.7",
    }
    for event in panel:
        for p in (source / event["id"]).glob("*.receipt.json"):
            provenance["raw_receipts"][event["id"] + "/" + p.name] = json.loads(
                p.read_text()
            )
        for p in (source / event["id"]).glob("*.command.json"):
            provenance["raw_commands"][event["id"] + "/" + p.name] = json.loads(
                p.read_text()
            )
    provenance["relate_windows_build"] = (
        json.loads((tools / "relate-build/build.json").read_text())
        if (tools / "relate-build/build.json").exists()
        else None
    )
    provenance["source_hashes"] = {
        str(p.relative_to(ROOT)): digest(p)
        for p in [
            ROOT / "scripts/sv_phylogeny_benchmark.py",
            ROOT / "scripts/package_sv_phylogeny.py",
            ROOT / "configs/sv_phylogeny.json",
        ]
    }
    write_json(out / "provenance.json", provenance)
    for p in out.glob("*.pdf"):
        data = p.read_bytes()
        assert data.startswith(b"%PDF-") and b"%%EOF" in data[-1024:]
        assert b"/FontFile2" in data and b"/Subtype /Image" not in data
    write_json(
        out / "manifest.json",
        {
            "files": {
                p.name: {"sha256": digest(p), "bytes": p.stat().st_size}
                for p in sorted(out.iterdir())
                if p.name != "manifest.json" and p.is_file()
            }
        },
    )
    print(f"Packaged {len(files)} synthetic inputs and results in {out}")


if __name__ == "__main__":
    main()
