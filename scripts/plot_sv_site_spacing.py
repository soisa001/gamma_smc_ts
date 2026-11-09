# /// script
# requires-python = ">=3.11"
# dependencies = ["matplotlib==3.10.3", "numpy==2.2.6"]
# ///
"""Count ALT alleles per record and plot consecutive SV start-position distances.

Input is the seven-column sites.tsv exported by run_sv_discovery_spectra.sh.
Distances count records, including same-position records; never cross contigs.
The class-specific panels measure spacing within each class, not to any SV.
"""

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import json
from pathlib import Path

from plot_sv_discovery_spectra import COLORS, GROUPS, digest, matplotlib, np, render

SPACING_GROUPS = ("All SVs",) + GROUPS


def classify(svtype, alt):
    if svtype not in ("", "."):
        bases = {part.split(":")[0].upper() for part in svtype.split(",")}
    else:
        bases = {part[1:-1].split(":")[0].upper()
                 if part.startswith("<") and part.endswith(">") else "UNKNOWN"
                 for part in alt.split(",")}
    if bases == {"DEL"}:
        return "DEL"
    if bases == {"INS"}:
        return "INS"
    return "Complex / other"


def audit(path):
    positions = {g: defaultdict(list) for g in SPACING_GROUPS}
    alt_hist = Counter()
    multi_by_group = Counter()
    nrecords = 0
    with path.open(encoding="utf-8") as handle:
        for lineno, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 7:
                raise ValueError(f"Expected seven fields at line {lineno}")
            chrom, pos, alt, svtype, *_ = fields
            try:
                pos = int(pos)
            except ValueError as exc:
                raise ValueError(f"Invalid POS at line {lineno}") from exc
            if not chrom or chrom == "." or pos < 1:
                raise ValueError(f"Invalid chromosome/POS at line {lineno}")
            if not alt or any(not a for a in alt.split(",")):
                raise ValueError(f"Empty ALT at line {lineno}")
            nalt = 0 if alt == "." else len(alt.split(","))
            alt_hist[nalt] += 1
            nrecords += 1
            # Reference-only rows are counted by the audit but are not SV starts.
            if not nalt:
                continue
            group = classify(svtype, alt)
            positions["All SVs"][chrom].append(pos)
            positions[group][chrom].append(pos)
            if nalt > 1:
                multi_by_group[group] += 1
    if not nrecords:
        raise ValueError("No input records")
    distances, details = {}, {}
    for group, contigs in positions.items():
        gaps = []
        shared_positions = records_at_shared = n = 0
        for starts in contigs.values():
            ordered = sorted(starts)
            gaps.extend(b - a for a, b in zip(ordered, ordered[1:]))
            multiplicities = Counter(starts)
            shared_positions += sum(k > 1 for k in multiplicities.values())
            records_at_shared += sum(k for k in multiplicities.values() if k > 1)
            n += len(starts)
        assert len(gaps) == n - len(contigs)
        zero = sum(g == 0 for g in gaps)
        positive = [g for g in gaps if g > 0]
        distances[group] = positive
        details[group] = {
            "records": n, "contigs": len(contigs), "adjacent_pairs": len(gaps),
            "zero_distance_pairs": zero, "positive_distance_pairs": len(positive),
            "positions_with_multiple_records": shared_positions,
            "records_at_shared_positions": records_at_shared,
            "minimum_positive_distance_bp": min(positive) if positive else None,
            "median_positive_distance_bp": float(np.median(positive)) if positive else None,
            "maximum_distance_bp": max(gaps) if gaps else None,
        }
    multi = sum(n for k, n in alt_hist.items() if k > 1)
    summary = {
        "records": nrecords, "multiallelic_records": multi,
        "multiallelic_record_percent": 100 * multi / nrecords,
        "biallelic_records": alt_hist[1], "records_without_ALT": alt_hist[0],
        "ALT_count_histogram": dict(sorted(alt_hist.items())),
        "multiallelic_records_by_group": {g: multi_by_group[g] for g in GROUPS},
        "spacing": details,
    }
    return distances, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sites-tsv", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--distance-bins", type=int, default=30)
    args = parser.parse_args()
    if args.distance_bins < 1:
        parser.error("--distance-bins must be positive")
    distances, summary = audit(args.sites_tsv)
    args.out.mkdir(parents=True, exist_ok=True)
    positive = [g for values in distances.values() for g in values]
    low, high = (min(positive), max(positive)) if positive else (1, 10)
    if low == high:
        low, high = low / 1.1, high * 1.1
    edges = np.geomspace(low * 0.999, high * 1.001, args.distance_bins + 1)
    matplotlib.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "DejaVu Sans"})
    render(distances, edges, "Consecutive SV start distance (bp; log scale)",
           "Inter-SV start distances", "inter_sv_distance", args.out,
           groups=SPACING_GROUPS, colors=("#455A64",) + COLORS,
           ylabel="Adjacent pairs", zero_counts={g: summary["spacing"][g]["zero_distance_pairs"]
                                                for g in SPACING_GROUPS})
    summary.update({
        "utc": datetime.now(timezone.utc).isoformat(),
        "sites_tsv": str(args.sites_tsv.resolve()), "sites_sha256": digest(args.sites_tsv),
        "script_sha256": digest(Path(__file__)),
        "plot_helper_sha256": digest(Path(__file__).with_name("plot_sv_discovery_spectra.py")),
        "distance_bins": args.distance_bins,
        "multiallelic_definition": "More than one comma-separated ALT entry on a single VCF record",
        "distance_definition": "Sorted consecutive POS differences within a chromosome; class panels use same-class neighbors",
        "distance_unit": "bp",
        "zero_handling": "Kept in summary and panel labels; excluded from log-scale histogram",
        "scope": "All ALT-bearing input records, one position per row including multiallelic rows; no deduplication or length filtering",
        "shared_positions": "Repeated CHROM/POS are counted separately; not proof of equivalent or multiallelic biological events",
        "versions": {"numpy": np.__version__, "matplotlib": matplotlib.__version__},
        "randomness": "None; deterministic sorting and summaries",
    })
    (args.out / "site_spacing_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    rows = ["ALT_entries\trecords"] + [f"{k}\t{v}" for k, v in sorted(summary["ALT_count_histogram"].items())]
    (args.out / "multiallelic_counts.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")
    print(f"Records: {summary['records']:,}")
    print(f"Multiallelic records: {summary['multiallelic_records']:,} ({summary['multiallelic_record_percent']:.4f}%)")
    for group, detail in summary["spacing"].items():
        print(f"{group}: {detail['adjacent_pairs']:,} adjacent pairs; "
              f"{detail['zero_distance_pairs']:,} zero distances; "
              f"{detail['positions_with_multiple_records']:,} positions with multiple records")
    print(f"Saved PNG/PDF, aggregate counts, and summary to {args.out.resolve()}")


if __name__ == "__main__":
    main()
