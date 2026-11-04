#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = ["matplotlib==3.10.8"]
# ///
"""Compare completed INS/DEL outputs without reading cohort genotypes/new.vcf.gz.

Uses the prepare allele manifest and the four kept/removed VCFs. No cloud I/O,
merging, inline plots, or changes to existing workflow artifacts.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import json
import logging
import os
from pathlib import Path
import sys
import time

import run_chr1_sv_merge as pilot

VERSION = 2
LITERAL_STATUSES = {"eligible", "above_sizemax", "ambiguous_sequence", "complex_replacement"}
SCOPES = ("merge_eligible", "all_ins_del")


def read_source(path):
    rows = {}
    with pilot.open_text(path) as handle:
        for record in csv.DictReader(handle, delimiter="\t"):
            if record["type"] not in {"INS", "DEL"}:
                continue
            if record["id"] in rows:
                raise ValueError("Duplicate source ID in source_alleles.tsv.gz")
            if record["status"] not in LITERAL_STATUSES | {"symbolic"}:
                raise ValueError("Unrecognized source allele status; cannot infer raw length safely")
            for key in ("pos", "start0", "end0"):
                record[key] = int(record[key])
            # In prepare schema 1, literal size was computed from abs(len(ALT)-len(REF)).
            # Symbolic SVLEN/END is deliberately NOT used as raw allele length here.
            record["raw_bp"] = int(record["size"]) if record["status"] in LITERAL_STATUSES else None
            record["merge_submitted"] = bool(int(record.get("merge_submitted", record["status"] == "eligible")))
            if record["raw_bp"] is not None and record["raw_bp"] <= 0:
                raise ValueError("Invalid raw allele length")
            rows[record["id"]] = record
    return rows


def allele_key(row):
    return row["chrom"], row["pos"], row["ref_sha256"], row["alt_sha256"]


def read_partition(path, kind, source, expected_samples_sha256):
    """Read site fields only; consume gzip to EOF to check CRC/truncation."""
    found, sample_header = {}, False
    previous = None
    with pilot.open_text(path) as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                samples = line.rstrip("\r\n").split("\t")[9:]
                actual = hashlib.sha256("".join(s + "\n" for s in samples).encode()).hexdigest()
                if actual != expected_samples_sha256:
                    raise ValueError(f"Sample manifest mismatch: {path.name}")
                sample_header = True
            if line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t", 8)[:8]
            if len(fields) != 8 or "," in fields[4]:
                raise ValueError(f"Invalid biallelic VCF row: {path.name}")
            chrom, position, ident, ref, alt = fields[:5]
            row = source.get(ident)
            if row is None or row["type"] != kind or not row["merge_submitted"] or ident in found:
                raise ValueError(f"Unknown/duplicate/wrong-type source allele: {path.name}")
            actual_key = chrom, int(position), hashlib.sha256(ref.encode()).hexdigest(), hashlib.sha256(alt.encode()).hexdigest()
            delta = len(alt) - len(ref)
            if actual_key != allele_key(row) or (row["raw_bp"] is not None and
                    delta != row["raw_bp"] * (1 if kind == "INS" else -1)):
                raise ValueError(f"REF/ALT, position, or raw length changed from source: {path.name}")
            coordinate = (chrom, int(position))
            if previous is not None and coordinate < previous:
                raise ValueError(f"Unsorted VCF: {path.name}")
            previous = coordinate
            found[ident] = pilot.info_parse(fields[7])
    if not sample_header:
        raise ValueError(f"Missing VCF header: {path.name}")
    return found


def validate_membership(source, partitions):
    retained, removed = set(), set()
    for kind in ("INS", "DEL"):
        kept, gone = partitions[(kind, "kept")], partitions[(kind, "removed")]
        expected = {ident for ident, row in source.items() if row["type"] == kind and row["merge_submitted"]}
        if set(kept) & set(gone) or set(kept) | set(gone) != expected:
            raise ValueError(f"{kind}: original eligible IDs != disjoint kept + removed IDs")
        clusters = {}
        expected_sizes = {}
        for ident, info in kept.items():
            n = int(info.get("NumCollapsed", 0))
            if n:
                cid = info.get("CollapseId")
                if not cid or cid in clusters:
                    raise ValueError(f"{kind}: missing/duplicate cluster ID")
                clusters[cid] = 0
                expected_sizes[cid] = n
        for info in gone.values():
            cid = info.get("MatchId")
            if cid not in clusters:
                raise ValueError(f"{kind}: removed allele has no retained representative")
            clusters[cid] += 1
        if clusters != expected_sizes:
            raise ValueError(f"{kind}: membership differs from NumCollapsed")
        retained.update(kept)
        removed.update(gone)
    return retained, removed


def union_intervals(intervals):
    """Merge 0-based half-open spans, including touching intervals, per contig."""
    result = []
    for chrom, start, end in sorted(intervals):
        if start < 0 or end < start:
            raise ValueError("Invalid deletion interval")
        if start == end:
            continue
        if result and result[-1][0] == chrom and start <= result[-1][2]:
            old_chrom, old_start, old_end = result[-1]
            result[-1] = old_chrom, old_start, max(old_end, end)
        else:
            result.append((chrom, start, end))
    return result


def summarize(rows, kind):
    unique = {allele_key(row): row for row in rows}
    literal = [row for row in rows if row["raw_bp"] is not None]
    unique_literal = [row for row in unique.values() if row["raw_bp"] is not None]
    # A simple deletion has changed-reference span == raw net deletion length.
    # Complex replacements have no uniquely defined deleted-base set from net size alone.
    deletions = [r for r in unique_literal if kind == "DEL" and r["end0"] - r["start0"] == r["raw_bp"]]
    intervals = union_intervals((r["chrom"], r["start0"], r["end0"]) for r in deletions)
    result = {
        "sv_records": len(rows), "unique_alleles": len(unique),
        "exact_duplicate_records": len(rows) - len(unique),
        "literal_records": len(literal), "symbolic_records_excluded_from_raw_bp": len(rows) - len(literal),
        "raw_bp_sum_records": sum(r["raw_bp"] for r in literal),
        "raw_bp_sum_unique_alleles": sum(r["raw_bp"] for r in unique_literal),
        "unique_inserted_bp": sum(r["raw_bp"] for r in unique_literal) if kind == "INS" else None,
        "unique_deleted_reference_bp": sum(end - start for _, start, end in intervals) if kind == "DEL" else None,
        "unique_deletion_intervals_used": len(deletions) if kind == "DEL" else None,
        "complex_del_alleles_excluded_from_union": len(unique_literal) - len(deletions) if kind == "DEL" else None,
    }
    return result, intervals


def spacing(rows, unique=False):
    if unique:
        rows = list({allele_key(r): r for r in rows}.values())
    ordered = sorted((r["chrom"], r["pos"]) for r in rows)
    # No gaps between contigs; ties between distinct same-position alleles remain zero.
    return [b[1] - a[1] for a, b in zip(ordered, ordered[1:]) if a[0] == b[0]]


def figures(output, summary, gaps):
    os.environ["MPLBACKEND"] = "Agg"
    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    matplotlib.rcParams.update({"font.size": 18, "axes.labelsize": 18, "axes.titlesize": 21,
                               "xtick.labelsize": 17, "ytick.labelsize": 17,
                               "pdf.fonttype": 42, "ps.fonttype": 42})
    colors = ["#0072B2", "#D55E00"]
    lookup = {(r["scope"], r["type"], r["callset"]): r for r in summary}
    paths = []

    def save(fig, stem):
        for ext in ("png", "pdf"):
            path = output / f"{stem}.{ext}"
            fig.savefig(path, dpi=180)
            paths.append(path)
        plt.close(fig)

    for scope in SCOPES:
        title = "Alleles submitted to Truvari" if scope == "merge_eligible" else "All INS/DEL alleles"
        for metric in ("sv_records", "unique_affected_bp"):
            fig, axes = plt.subplots(1, 2, figsize=(11, 8.5), layout="constrained")
            fig.suptitle(title, fontsize=22)
            for ax, kind in zip(axes, ("DEL", "INS")):
                field = metric if metric == "sv_records" else ("unique_deleted_reference_bp" if kind == "DEL" else "unique_inserted_bp")
                values = [lookup[(scope, kind, name)][field] for name in ("before", "after")]
                ax.bar([0, 1], values, color=colors)
                ax.set(xticks=[0, 1], xticklabels=["Before", "After"], title=kind)
                ax.set_ylabel("SV allele records" if metric == "sv_records" else
                              ("Deleted reference bp (union)" if kind == "DEL" else "Inserted bp (unique alleles)"))
                ax.yaxis.set_major_locator(MaxNLocator(integer=True))
                ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 6), useMathText=True)
            save(fig, f"{scope}.{metric}")
        for kind in ("DEL", "INS"):
            fig, ax = plt.subplots(figsize=(11, 8.5), layout="constrained")
            largest = max([10] + [max(gaps[(scope, kind, name, "records")], default=0)
                                  for name in ("before", "after")])
            tail = largest * 1.1
            for name, color in zip(("before", "after"), colors):
                values = sorted(gaps[(scope, kind, name, "records")])
                if values:
                    ax.step([0] + values + [tail], [0] + [(i + 1) / len(values) for i in range(len(values))] + [1],
                            where="post", color=color, linewidth=2.5, label=name.capitalize())
                else:
                    ax.plot([], [], color=color, label=f"{name.capitalize()} (no adjacent pairs)")
            ax.set_xscale("symlog", linthresh=1)
            ax.set(xlim=(0, tail), ylim=(0, 1.02), xlabel="Adjacent SV POS distance (bp; zero retained)",
                   ylabel="Cumulative fraction", title=f"{kind}: {title.lower()}")
            ax.legend(loc="lower center", bbox_to_anchor=(.5, 1.02), ncols=2, frameon=False)
            save(fig, f"{scope}.{kind}.inter_sv_distance")
    return paths


def run(args):
    root, output = args.root, args.output
    prepared_path = root / "prepare" / "done.json"
    prepared = json.loads(prepared_path.read_text())
    if prepared["schema"] != 1:
        raise ValueError("Expected prepare schema 1")
    source_path = root / "prepare" / "source_alleles.tsv.gz"
    sample_path = root / "prepare" / "samples.txt"
    paths = {"source_alleles": source_path, "samples": sample_path}
    for kind in ("INS", "DEL"):
        for role in ("kept", "removed"):
            path = root / "merge" / f"{kind}.{role}.vcf.gz"
            if not Path(str(path) + ".tbi").is_file():
                raise ValueError(f"Indexed output not yet present: {path.name}. Wait for indexing to finish.")
            paths[f"{kind}.{role}"] = path
    before_stats = {name: pilot.common.identity(path) for name, path in paths.items()}
    logging.info("Hashing only allele manifest, sample list, and separate INS/DEL outputs")
    with ThreadPoolExecutor(max_workers=args.threads) as pool:
        digests = dict(zip(paths, pool.map(pilot.common.digest, paths.values())))
    for name, path in (("source_alleles", source_path), ("samples", sample_path)):
        if digests[name] != prepared["artifacts"][path.name]:
            raise ValueError(f"Preparation artifact checksum mismatch: {path.name}")
    key = {"comparison_version": VERSION, "inputs_sha256": digests,
           "source_parameters": prepared["key"], "seed": pilot.SEED}
    if pilot.valid_manifest(output, key):
        if any(pilot.common.identity(path) != before_stats[name] for name, path in paths.items()):
            raise ValueError("An input changed during cache verification; wait for completion and rerun")
        logging.info("Reusing checksum-verified comparison")
        return
    (output / "done.json").unlink(missing_ok=True)
    logging.info("Loading original allele metadata (no cohort genotype parsing)")
    source = read_source(source_path)
    sample_header_sha256 = hashlib.sha256(sample_path.read_text(encoding="utf-8").encode()).hexdigest()
    tasks = [(kind, role) for kind in ("INS", "DEL") for role in ("kept", "removed")]
    def load(task):
        kind, role = task
        return read_partition(paths[f"{kind}.{role}"], kind, source, sample_header_sha256)
    with ThreadPoolExecutor(max_workers=args.threads) as pool:
        partitions = dict(zip(tasks, pool.map(load, tasks)))
    kept, removed = validate_membership(source, partitions)
    summary, distance_rows, gaps, artifacts = [], [], {}, []
    for scope in SCOPES:
        for kind in ("DEL", "INS"):
            before = [r for r in source.values() if r["type"] == kind and
                      (scope == "all_ins_del" or r["merge_submitted"])]
            after = [r for r in before if r["id"] not in removed]
            for name, rows in (("before", before), ("after", after)):
                stats, intervals = summarize(rows, kind)
                summary.append({"scope": scope, "type": kind, "callset": name, **stats})
                if kind == "DEL":
                    bed = output / f"{scope}.{name}.deleted_union.bed"
                    with open(bed, "w", encoding="utf-8") as handle:
                        for chrom, start, end in intervals:
                            handle.write(f"{chrom}\t{start}\t{end}\n")
                    artifacts.append(bed)
                for basis in ("records", "unique_alleles"):
                    values = spacing(rows, unique=basis == "unique_alleles")
                    gaps[(scope, kind, name, basis)] = values
                    distance_rows.append({"scope": scope, "type": kind, "callset": name, "basis": basis,
                                          "n_pairs": len(values), "n_zero": values.count(0),
                                          "n_le_500_bp": sum(v <= 500 for v in values),
                                          "n_le_1000_bp": sum(v <= 1000 for v in values), **pilot.quantiles(values)})
    delta_rows = []
    for before, after in zip(summary[::2], summary[1::2]):
        for metric in ("sv_records", "unique_alleles", "raw_bp_sum_records", "raw_bp_sum_unique_alleles",
                       "unique_inserted_bp", "unique_deleted_reference_bp"):
            a, b = before[metric], after[metric]
            if a is not None:
                delta_rows.append({"scope": before["scope"], "type": before["type"], "metric": metric,
                                   "before": a, "after": b, "reduction": a - b,
                                   "percent_reduction": 100 * (a - b) / a if a else None})
    for name, rows in (("before_after.tsv", summary), ("changes.tsv", delta_rows), ("inter_sv_distance_summary.tsv", distance_rows)):
        path = output / name
        pilot.tsv(path, rows, list(rows[0]))
        artifacts.append(path)
    artifacts.extend(figures(output, summary, gaps))
    for name, path in paths.items():
        if pilot.common.identity(path) != before_stats[name]:
            raise ValueError(f"Input changed during comparison: {path.name}; wait for it to finish and rerun")
    report = {
        "comparison_version": VERSION, "conservation_passed": True,
        "eligible_kept": len(kept), "eligible_removed": len(removed), "summary": summary,
        "input_paths": {name: str(path) for name, path in paths.items()}, "cache_key": key,
        "definitions": {
            "merge_eligible": "All alleles submitted to Truvari (merge_submitted); legacy manifests use status=eligible. Native tool filters may retain records unchanged.",
            "unique_allele": "Exact CHROM, POS, REF, ALT; sequence identity stored as SHA-256. No additional merging.",
            "raw_bp": "INS: len(ALT)-len(REF); DEL: len(REF)-len(ALT). No AC/AF weighting or INFO/SVLEN substitution.",
            "unique_inserted_bp": "Sum of raw lengths once per distinct literal insertion allele at its locus.",
            "unique_deleted_reference_bp": "Union of simple literal deletion intervals, anchor/padding excluded, each reference bp counted once.",
            "exclusions": "Symbolic lengths are unavailable from literal strings; complex DELs count in raw sums but cannot define a deletion union from net length alone.",
            "all_ins_del": "Original manifest alleles minus removed IDs; includes unchanged pass-through alleles without scanning their genotype VCF.",
            "spacing": "Adjacent POS differences separately within type and contig; zero gaps retained; plots use record basis.",
        },
    }
    pilot.common.json_write(output / "summary.json", report)
    artifacts.append(output / "summary.json")
    pilot.manifest(output, key, artifacts)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args(argv)
    if args.threads < 1:
        parser.error("--threads must be positive")
    args.root = args.root.resolve()
    args.output = args.output.resolve() if args.output else args.root / "compare_ins_del"
    if args.output == args.root or args.output in [args.root / x for x in ("prepare", "merge", "inputs", "compare")]:
        parser.error("Use a separate comparison output directory")
    args.output.mkdir(parents=True, exist_ok=True)
    lock = args.output / ".running.lock"
    descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    os.write(descriptor, f"{os.getpid()}\n".encode())
    os.close(descriptor)
    logging.basicConfig(filename=args.output / "comparison.log", level=logging.INFO, force=True,
                        format="%(asctime)s %(levelname)s %(message)s")
    try:
        logging.info("PARAMETERS %s", json.dumps(vars(args), default=str))
        logging.info("PROVENANCE python=%s script_sha256=%s seed=%s", sys.version,
                     pilot.common.digest(__file__), pilot.SEED)
        start = time.monotonic()
        print("Comparing completed INS/DEL outputs; tables and plots stay on VM disk.", flush=True)
        run(args)
        logging.info("Completed elapsed_seconds=%.3f", time.monotonic() - start)
        print(f"Comparison ready: {args.output}", flush=True)
    except Exception:
        logging.exception("Comparison failed")
        raise SystemExit(f"Comparison failed; inspect {args.output / 'comparison.log'} on the VM") from None
    finally:
        lock.unlink(missing_ok=True)
        logging.shutdown()


if __name__ == "__main__":
    main()
