# /// script
# requires-python = ">=3.11"
# dependencies = ["matplotlib==3.10.3", "numpy==2.2.6"]
# ///
"""Plot discovery ALT counts and lengths; never infer reference genotypes.

Input columns (no header): CHROM POS ALT SVTYPE SVLEN END AC.
Only aggregate tables/figures are written. No notebook display is used.
"""

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path

os.environ["MPLBACKEND"] = "Agg"
for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(key, "4")

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter
import numpy as np

GROUPS = ("DEL", "INS", "Complex / other")
COLORS = ("#2563A6", "#D46A27", "#7756A3")


def digest(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def parse_row(line, nsamples):
    chrom, pos, alt, svtype, svlen, end, ac = line.rstrip("\n").split("\t")
    if "," in alt:
        raise ValueError("Multiallelic record: split alleles before plotting")
    count = int(ac)
    if not 0 <= count <= 2 * nsamples:
        raise ValueError("AC outside [0, 2 * sample count]")
    if svtype in ("", "."):
        svtype = alt[1:-1] if alt.startswith("<") and alt.endswith(">") else "UNKNOWN"
    base = svtype.split(":")[0].upper()
    group = base if base in ("DEL", "INS") else "Complex / other"
    length, method = None, "unavailable"
    # Breakend spans are not event lengths. Multivalued complex SVLEN is ambiguous.
    if base not in ("BND", "TRA", "CTX"):
        if svlen not in ("", ".") and "," not in svlen:
            value = int(svlen)
            if value:
                length, method = abs(value), "abs(SVLEN)"
        if length is None and base in ("DEL", "DUP", "INV") and end not in ("", "."):
            span = int(end) - int(pos)
            if span > 0:
                length, method = span, "END-POS"
    return group, svtype, count, length, method


def render(values, edges, xlabel, title, stem, out):
    fig, axes = plt.subplots(3, 1, figsize=(8.5, 11), layout="constrained")
    fig.suptitle(title, fontsize=21)
    rows = ["group\tbin_left\tbin_right\trecords"]
    for group, color, ax in zip(GROUPS, COLORS, axes):
        data = np.asarray(values[group])
        counts, _ = np.histogram(data, bins=edges)
        if int(counts.sum()) != len(data):
            raise ValueError("Histogram lost observations")
        ax.stairs(counts, edges, fill=True, color=color, alpha=0.85, linewidth=1)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(edges[0], edges[-1])
        ax.set_ylim(0.8, max(2, int(counts.max()) * 2.5))
        ax.set_title(f"{group}   |   n = {len(data):,}", loc="left", fontsize=18)
        ax.set_ylabel("SV records", fontsize=17)
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:,.3g}"))
        ax.yaxis.set_major_locator(LogLocator(base=10, numticks=4))
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.tick_params(labelsize=15)
        ax.grid(axis="y", which="major", alpha=0.2)
        ax.spines[["right", "top"]].set_visible(False)
        if not len(data):
            ax.text(0.5, 0.5, "No eligible records", ha="center", transform=ax.transAxes)
        for left, right, count in zip(edges[:-1], edges[1:], counts):
            rows.append(f"{group}\t{left:.12g}\t{right:.12g}\t{int(count)}")
    axes[-1].set_xlabel(xlabel, fontsize=18)
    for ext in ("png", "pdf"):
        fig.savefig(out / f"{stem}.{ext}", dpi=180)
    plt.close(fig)
    (out / f"{stem}.bins.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sites-tsv", type=Path, required=True)
    parser.add_argument("--samples-file", type=Path, required=True)
    parser.add_argument("--source-bcf", required=True, help="Provenance label; not opened")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--length-bins", type=int, default=30)
    args = parser.parse_args()
    if args.length_bins < 1:
        parser.error("--length-bins must be positive")
    samples = args.samples_file.read_text(encoding="utf-8").splitlines()
    if not samples or any(not s for s in samples) or len(set(samples)) != len(samples):
        raise ValueError("Sample manifest is empty or has blank/duplicate sample names")
    nsamples = len(samples)
    acs, lengths = ({g: [] for g in GROUPS} for _ in range(2))
    types, missing, zero_ac, methods = Counter(), Counter(), Counter(), Counter()
    nrecords = 0
    with args.sites_tsv.open(encoding="utf-8") as handle:
        for lineno, line in enumerate(handle, 1):
            try:
                group, svtype, ac, length, method = parse_row(line, nsamples)
            except ValueError as exc:
                raise ValueError(f"Invalid sites TSV at line {lineno}: {exc}") from exc
            nrecords += 1
            types[svtype] += 1
            if ac:
                acs[group].append(ac)
            else:
                zero_ac[group] += 1
            if length is None:
                missing[group] += 1
            else:
                lengths[group].append(length / 1000)
            methods[method] += 1
    if not nrecords:
        raise ValueError("No input records")
    args.out.mkdir(parents=True, exist_ok=True)
    matplotlib.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "DejaVu Sans"})
    # Exact singleton/doubleton bins, then intervals 3-4, 5-8, 9-16, ...
    top = max([1] + [max(v) for v in acs.values() if v])
    upper = [1]
    while upper[-1] < top:
        upper.append(upper[-1] * 2)
    ac_edges = np.array([0.5] + [x + 0.5 for x in upper])
    render(acs, ac_edges, "Observed ALT allele count (log scale)",
           "Discovery ALT-count spectrum", "discovery_afs", args.out)
    all_lengths = [x for group in GROUPS for x in lengths[group]]
    if all_lengths:
        low, high = min(all_lengths), max(all_lengths)
        if low == high:
            low, high = low / 1.1, high * 1.1
        length_edges = np.geomspace(low * 0.999, high * 1.001, args.length_bins + 1)
    else:
        length_edges = np.array([1.0, 10.0])
    render(lengths, length_edges, "SV length (kb; log scale)",
           "SV length distribution", "length_distribution", args.out)
    exact = ["group\tobserved_alt_count\trecords"]
    for group in GROUPS:
        for count, records in sorted(Counter(acs[group]).items()):
            exact.append(f"{group}\t{count}\t{records}")
    (args.out / "discovery_afs.exact.tsv").write_text("\n".join(exact) + "\n", encoding="utf-8")
    summary = {
        "utc": datetime.now(timezone.utc).isoformat(),
        "parameters": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
        "script_sha256": digest(Path(__file__)), "sites_sha256": digest(args.sites_tsv),
        "samples_sha256": digest(args.samples_file), "samples": nsamples, "records": nrecords,
        "svtype_records": dict(types), "length_methods": dict(methods),
        "grouping": "DEL and INS including subtypes; all remaining types are Complex / other",
        "afs": "INFO/AC as observed discovery dosage, not AC/AN or validated cohort frequency",
        "length": "abs(SVLEN); END-POS fallback for DEL/DUP/INV; exclude BND/TRA/CTX lengths",
        "bins": "AC: 1, 2, 3-4, 5-8, ...; length: shared geometric bins; counts not densities",
        "length_scope": "All input records with defined positive length; no additional length filter",
        "versions": {"numpy": np.__version__, "matplotlib": matplotlib.__version__},
        "threads": {k: os.environ[k] for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")},
        "randomness": "None; deterministic summaries and plots",
        "groups": {g: {"afs_records": len(acs[g]), "length_records": len(lengths[g]),
                       "length_unavailable": missing[g], "zero_AC_excluded_from_afs": zero_ac[g],
                       "length_le_10kb": sum(x <= 10 for x in lengths[g])} for g in GROUPS},
    }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"samples": nsamples, "records": nrecords, "groups": summary["groups"]}, indent=2))
    print(f"Saved PNG, PDF, aggregate TSVs, and summary.json to {args.out.resolve()}")


if __name__ == "__main__":
    main()
