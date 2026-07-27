"""Compare a benchmark run of this branch against the AOU_run baseline.

Reads the stdout logs and /usr/bin/time -v reports written by
.github/workflows/decoder-benchmark.yml and prints the comparison.

Total wall time is a poor measure at benchmark scale: it is dominated by fixed
startup -- converting the tree sequence through a Python subprocess, reading it,
and building the 490 MB flow-field cache -- which is identical for both
binaries and hides the decode entirely. Both binaries print per-phase timers,
and single-threaded those four sum to the decode, so that is what is compared.
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--directory", default=".")
    parser.add_argument("--length", type=int, required=True)
    parser.add_argument("--random-pairs", type=int, required=True)
    args = parser.parse_args()

    root = Path(args.directory)
    names = ("base", "opt1", "optN", "scale")
    totals = {name: total_time(root / f"{name}.time") for name in names}
    cpu = {name: decode_cpu(root / f"{name}.log") for name in names}
    wall = {name: decode_wall(root / f"{name}.log") for name in names}
    cores = os.cpu_count() or 1

    print(f"{'run':6} {'total s':>9} {'decode CPU s':>13} {'decode wall s':>14} {'peak GB':>9}")
    for name in names:
        shown = f"{wall[name]:.2f}" if wall[name] is not None else "n/a"
        print(f"{name:6} {totals[name][0]:9.1f} {cpu[name]:13.2f} {shown:>14} {totals[name][1]:9.2f}")
    print()

    print(f"Decode only, single thread vs AOU_run: {cpu['base'] / cpu['opt1']:.2f}x")
    if wall["optN"]:
        print(f"Thread scaling on {cores} vCPU:             {wall['opt1'] / wall['optN']:.2f}x")
        print(f"Decode only, {cores} vCPU vs AOU_run:       {cpu['base'] / wall['optN']:.2f}x")
    print()

    # The two binaries must agree on the statistic itself.
    base = pd.read_csv(root / "base.tsv", sep="\t")
    opt = pd.read_csv(root / "opt1.tsv", sep="\t")
    assert len(base) == len(opt), (len(base), len(opt))
    delta = np.abs(base["mean_p_tmrca_lt_threshold"] - opt["mean_p_tmrca_lt_threshold"])
    print(f"mean_p_tmrca_lt_threshold vs baseline: max |delta| = {np.nanmax(delta):.3e}")
    difference = np.abs(base["mean_tmrca_generations"] - opt["mean_tmrca_generations"])
    relative = difference / np.maximum(np.abs(base["mean_tmrca_generations"]), 1e-12)
    print(f"mean_tmrca_generations:                max rel     = {np.nanmax(relative):.3e}")
    print()

    per_gbp_pair = wall["scale"] / (args.random_pairs * args.length / 1e9)
    hours = per_gbp_pair * 3.1 * 1e5 / 3600
    print(f"Throughput: {per_gbp_pair:.4f} s per Gbp per pair, wall, on {cores} vCPU")
    print(f"3.1 Gbp x 100,000 pairs at this rate:  {hours:.1f} h")
    print(f"...at the same per-core efficiency on 32 cores: {hours * cores / 32:.1f} h")

    bits = (root / "scale.bits").stat().st_size
    positions = args.length // 1000
    print()
    print(f"Bit matrix: {bits / 1048576:.1f} MB compressed for "
          f"{args.random_pairs} pairs x 2 thresholds x {positions} positions "
          f"({args.random_pairs * 2 * positions / 8 / 1048576:.1f} MB raw)")


if __name__ == "__main__":
    main()
