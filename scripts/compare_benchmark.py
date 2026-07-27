"""Compare a benchmark run of this branch against the AOU_run baseline.

Reads the stdout logs and `/usr/bin/time -v` reports written by
.github/workflows/decoder-benchmark.yml and prints an ablation, so the speedup
is attributed rather than just asserted:

    base   AOU_run, one thread
    exact  this branch, one thread, upstream numerics, --exact_recent_stats
           (same boost::math::gamma_p per element as upstream, so the difference
           from `base` is the CPU/memory work alone: hoisted emission fill,
           bit-packed genotypes, thread-local flow-field scratch)
    legacy this branch, one thread, upstream numerics
           (the difference from `exact` is the lookup tables)
    optN   this branch, every core, upstream numerics
           (the difference from `legacy` is threading)
    fixed  this branch, every core, shipping defaults
           (corrected exp10 and backward alignment; the difference from `optN`
           is how far those two corrections move the statistic)

Total wall time is a poor measure at benchmark scale: it is dominated by fixed
startup -- converting the tree sequence through a Python subprocess, reading it,
and building the 490 MB flow-field cache -- which is identical for every run and
hides the decode. Both binaries print per-phase timers which sum to the decode
when single-threaded, so that is what is compared.
"""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd


PHASES = (
    "Emissions preparation time",
    "Forward pass time",
    "Backward pass time",
    "Output time",
)

# Columns both binaries write, and which must agree.
SHARED_COLUMNS = ("mean_p_tmrca_lt_threshold", "mean_tmrca_generations")


def total_time(path: Path) -> tuple[float, float]:
    """Wall seconds and peak RSS in GB from a `/usr/bin/time -v` report."""
    text = path.read_text()
    match = re.search(r"Elapsed \(wall clock\) time.*?:\s*([\d:.]+)", text)
    seconds = 0.0
    for part in match.group(1).split(":"):
        seconds = seconds * 60 + float(part)
    rss = float(re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", text).group(1))
    return seconds, rss / 1048576.0


def decode_cpu(path: Path) -> float:
    text = path.read_text()
    total = 0.0
    for label in PHASES:
        match = re.search(re.escape(label) + r":\s*([\d.]+)\s*secs", text)
        if match:
            total += float(match.group(1))
    return total


def decode_wall(path: Path) -> float | None:
    match = re.search(r"Decoding wall time:\s*([\d.]+)\s*secs", path.read_text())
    return float(match.group(1)) if match else None


def phase_breakdown(path: Path) -> dict[str, float]:
    text = path.read_text()
    out = {}
    for label in PHASES:
        match = re.search(re.escape(label) + r":\s*([\d.]+)\s*secs", text)
        out[label] = float(match.group(1)) if match else 0.0
    return out


def compare_outputs(left: Path, right: Path, label: str) -> None:
    a = pd.read_csv(left, sep="\t")
    b = pd.read_csv(right, sep="\t")
    assert len(a) == len(b), (label, len(a), len(b))
    print(f"  {label}")
    for column in SHARED_COLUMNS:
        if column not in a.columns or column not in b.columns:
            continue
        difference = np.abs(a[column].to_numpy() - b[column].to_numpy())
        scale = np.maximum(np.abs(a[column].to_numpy()), 1e-12)
        print(f"    {column:28} max |delta| = {np.nanmax(difference):.3e}"
              f"   max rel = {np.nanmax(difference / scale):.3e}")
    # Counts only exist on the new binary; compare them when both sides have them.
    for column in a.columns:
        if column.startswith("n_recent_") and column in b.columns:
            differing = int((a[column].to_numpy() != b[column].to_numpy()).sum())
            total = int(a[column].to_numpy().sum())
            print(f"    {column:28} positions differing: {differing}/{len(a)}"
                  f"   (total calls {total})")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--directory", default=".")
    parser.add_argument("--length", type=int, required=True)
    parser.add_argument("--random-pairs", type=int, required=True)
    parser.add_argument("--pairs", type=int, required=True, help="pairs in the ablation runs")
    args = parser.parse_args()

    root = Path(args.directory)
    names = ("base", "exact", "legacy", "optN", "fixed", "scale")
    totals = {name: total_time(root / f"{name}.time") for name in names}
    cpu = {name: decode_cpu(root / f"{name}.log") for name in names}
    wall = {name: decode_wall(root / f"{name}.log") for name in names}
    cores = os.cpu_count() or 1

    print(f"{'run':7} {'total s':>9} {'decode CPU s':>13} {'decode wall s':>14} {'peak GB':>9}")
    for name in names:
        shown = f"{wall[name]:.2f}" if wall[name] is not None else "n/a"
        print(f"{name:7} {totals[name][0]:9.1f} {cpu[name]:13.2f} {shown:>14} {totals[name][1]:9.2f}")
    print()

    print("Decode-time ablation, single thread unless stated:")
    print(f"  AOU_run baseline                        {cpu['base']:8.2f} s")
    print(f"  + CPU/memory work, exact gamma_p        {cpu['exact']:8.2f} s"
          f"   ({cpu['base'] / cpu['exact']:5.2f}x)")
    print(f"  + lookup tables                         {cpu['legacy']:8.2f} s"
          f"   ({cpu['exact'] / cpu['legacy']:5.2f}x)")
    if wall["optN"]:
        print(f"  + {cores} threads                            {wall['optN']:8.2f} s"
              f"   ({wall['legacy'] / wall['optN']:5.2f}x)")
        print(f"  overall                                          "
              f"   ({cpu['base'] / wall['optN']:5.2f}x)")
    print()

    print("Where the single-threaded time goes (seconds):")
    print(f"  {'phase':32} {'base':>8} {'exact':>8} {'legacy':>8}")
    base_phases = phase_breakdown(root / "base.log")
    exact_phases = phase_breakdown(root / "exact.log")
    opt_phases = phase_breakdown(root / "legacy.log")
    for label in PHASES:
        print(f"  {label:32} {base_phases[label]:8.2f} {exact_phases[label]:8.2f} "
              f"{opt_phases[label]:8.2f}")
    print()

    print("Output agreement:")
    compare_outputs(root / "base.tsv", root / "legacy.tsv",
                    "AOU_run vs this branch at upstream numerics (must be ~0: no accidental change)")
    compare_outputs(root / "exact.tsv", root / "legacy.tsv",
                    "exact gamma_p vs lookup tables, same binary (lookup-table error alone)")
    compare_outputs(root / "legacy.tsv", root / "optN.tsv", f"1 thread vs {cores} threads")
    compare_outputs(root / "optN.tsv", root / "fixed.tsv",
                    "upstream numerics vs shipping defaults (what the two corrections change)")
    print()

    per_gbp_pair = wall["scale"] / (args.random_pairs * args.length / 1e9)
    hours = per_gbp_pair * 3.1 * 1e5 / 3600
    print(f"Throughput: {per_gbp_pair:.4f} s per Gbp per pair, wall, on {cores} vCPU")
    print(f"  scale run used {cpu['scale']:.1f} s CPU in {wall['scale']:.2f} s wall "
          f"({cpu['scale'] / wall['scale']:.2f}x on {cores} cores)")
    print(f"  3.1 Gbp x 100,000 pairs at this rate:  {hours:.1f} h")
    print(f"  at the same per-core efficiency on 32 cores: {hours * cores / 32:.1f} h")

    bits = (root / "scale.bits").stat().st_size
    positions = args.length // 1000
    raw = args.random_pairs * 2 * positions / 8
    print()
    print(f"Bit matrix: {bits / 1048576:.1f} MB compressed, {raw / 1048576:.1f} MB raw "
          f"({raw / max(bits, 1):.0f}x) for {args.random_pairs} pairs x 2 thresholds "
          f"x {positions} positions")


if __name__ == "__main__":
    main()
