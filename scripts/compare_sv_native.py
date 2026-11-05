#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = []
# ///
"""Summarize completed native SV merges without emitting sample/allele data."""
from __future__ import annotations

import argparse
import csv
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import logging
import os
from pathlib import Path
import sys

import compare_chr1_sv_alleles as metrics
import run_chr1_sv_merge as pilot

VERSION = 2


def read_sites(args, path, label):
    rows, membership = [], Counter()
    cmd = [args.bcftools, "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n", path]
    with pilot.stream(args, cmd, label) as handle:
        for line in handle:
            chrom, pos, ref, alt = line.rstrip("\r\n").split("\t")
            pos = int(pos)
            if chrom != args.chrom or "," in alt:
                raise ValueError("Expected one chromosome and biallelic input")
            ref_hash = hashlib.sha256(ref.encode()).hexdigest()
            alt_hash = hashlib.sha256(alt.encode()).hexdigest()
            membership[(chrom, pos, ref_hash, alt_hash)] += 1
            parsed = pilot.allele(ref, alt, pos, {}, minimum=1, maximum=-1)
            if parsed is None or parsed[0] not in {"INS", "DEL"}:
                continue
            kind, size, start, end, status = parsed
            rows.append({"chrom": chrom, "pos": pos, "type": kind,
                "ref_sha256": ref_hash, "alt_sha256": alt_hash,
                "start0": start, "end0": end,
                "raw_bp": None if status == "symbolic" else size})
    return rows, membership


def compare(root, chrom, bcftools, threads=4, similarity=.7):
    if similarity not in (.7, .9):
        raise ValueError("Supported similarities are 0.7 and 0.9")
    tag = f"p{round(similarity * 100)}"
    base = root / chrom
    merge = base / tag
    paths = {"before": base / "prepare/input.vcf.gz", "after": merge / f"{chrom}.vcf.gz",
             "removed": merge / f"{chrom}.collapsed.vcf.gz",
             "prepare_manifest": base / "prepare/done.json", "merge_manifest": merge / "done.json"}
    for path in paths.values():
        if not path.is_file():
            raise ValueError(f"Not completed: {path}")
    prepared = json.loads(paths["prepare_manifest"].read_text())
    merged = json.loads(paths["merge_manifest"].read_text())
    if merged["key"].get("genotype_consolidation") != "native" or merged["key"]["similarity"] != similarity:
        raise ValueError(f"Expected completed native {similarity} merge")
    if merged["key"]["prepare"] != prepared:
        raise ValueError("Merge does not match current preparation manifest")
    final_minimum = merged["key"].get("final_min_raw_bp")
    for label in ("before", "after", "removed"):
        paths[label + "_index"] = Path(str(paths[label]) + ".tbi")
    snapshots = {k: pilot.common.identity(p) for k, p in paths.items()}
    output = base / f"compare_native_{tag}"
    (output / "logs").mkdir(parents=True, exist_ok=True)
    args = argparse.Namespace(root=output, chrom=chrom, bcftools=bcftools)
    logging.basicConfig(filename=output / "comparison.log", level=logging.INFO, force=True,
                        format="%(asctime)s %(levelname)s %(message)s")
    version = pilot.command(args, [bcftools, "--version"], "bcftools-version", True).splitlines()[0]
    key = {"version": VERSION, "input_snapshots": snapshots, "prepare": prepared, "merge": merged,
           "bcftools": version, "seed": pilot.SEED}
    logging.info("PARAMETERS chrom=%s threads=%s paths=%s", chrom, threads, paths)
    logging.info("PROVENANCE python=%s script_sha256=%s seed=%s", sys.version,
                 pilot.common.digest(__file__), pilot.SEED)
    # Input identities + completion manifests avoid hashing huge cohort files again.
    # A fresh scan checks decoder success and exact before = after + removed membership.
    if pilot.valid_manifest(output, key):
        return output
    lock = output / ".running.lock"
    fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    os.write(fd, f"{os.getpid()}\n".encode())
    os.close(fd)
    try:
        (output / "done.json").unlink(missing_ok=True)
        labels = ("before", "after", "removed")
        with ThreadPoolExecutor(max_workers=min(threads, 3)) as pool:
            scans = dict(zip(labels, pool.map(lambda role: read_sites(args, paths[role], role), labels)))
        if scans["before"][1] != scans["after"][1] + scans["removed"][1]:
            raise ValueError("Original site multiset differs from retained + collapsed sites")
        if sum(scans["before"][1].values()) != prepared["key"]["records"]:
            raise ValueError("Original count differs from preparation manifest")
        stats, changes, distances = [], [], []
        for kind in ("INS", "DEL"):
            pair = []
            for role in ("before", "after"):
                rows = [r for r in scans[role][0] if r["type"] == kind and
                        (final_minimum is None or (r["raw_bp"] is not None and r["raw_bp"] >= final_minimum))]
                summary, _ = metrics.summarize(rows, kind)
                summary.update({"raw_length_bp_" + k: v for k, v in
                    pilot.quantiles([r["raw_bp"] for r in rows]).items()})
                stats.append({"chrom": chrom, "type": kind, "callset": role, **summary})
                pair.append(summary)
                for basis in ("records", "unique_alleles"):
                    gaps = metrics.spacing(rows, unique=basis == "unique_alleles")
                    distances.append({"chrom": chrom, "type": kind, "callset": role, "basis": basis,
                        "n_pairs": len(gaps), "n_zero": gaps.count(0),
                        "n_le_500_bp": sum(g <= 500 for g in gaps),
                        "n_le_1000_bp": sum(g <= 1000 for g in gaps), **pilot.quantiles(gaps)})
            for field in pair[0]:
                a, b = pair[0][field], pair[1][field]
                if a is not None and b is not None:
                    changes.append({"chrom": chrom, "type": kind, "metric": field,
                        "before": a, "after": b, "reduction": a - b,
                        "percent_reduction": 100 * (a - b) / a if a else None})
        artifacts = []
        for filename, rows in (("changes.tsv", changes), ("before_after.tsv", stats),
                               ("inter_sv_distance_summary.tsv", distances)):
            path = output / filename
            pilot.tsv(path, rows, list(rows[0]))
            artifacts.append(path)
        report = {"chrom": chrom, "similarity": similarity, "conservation_passed": True,
            "posthoc_min_raw_bp": final_minimum,
            "all_record_counts": {role: sum(scans[role][1].values()) for role in labels},
            "ins_del_summary": stats,
            "definitions": {
                "baseline": "Prepared input after chromosome selection and multiallelic splitting. If posthoc_min_raw_bp is set, apply that raw-length cutoff to BOTH baseline and native retained representatives for INS/DEL statistics; conservation still covers ALL records.",
                "unique_allele": "Exact CHROM/POS/REF/ALT, sequences hashed with SHA-256; no AC weighting.",
                "raw_bp": "Absolute len(ALT)-len(REF) for literal alleles; symbolic lengths excluded.",
                "unique_inserted_bp": "Sum raw lengths once per exact insertion allele at its locus, not unique sequence k-mers.",
                "unique_deleted_reference_bp": "Union of simple literal deletions after shared prefix/suffix trimming. Complex replacements and symbolic DELs excluded.",
                "spacing": "Adjacent POS distances within type, including zero-distance ties.",
                "validation": "bcftools decode to EOF; full site-multiset conservation; recorded input sizes/mtimes unchanged. No new whole-VCF SHA verification."}}
        pilot.common.json_write(output / "summary.json", report)
        artifacts.append(output / "summary.json")
        if snapshots != {k: pilot.common.identity(p) for k, p in paths.items()}:
            raise ValueError("An input changed while computing comparison")
        pilot.manifest(output, key, artifacts)
        return output
    finally:
        lock.unlink(missing_ok=True)


def side_by_side(root, chrom, outputs):
    """Join verified per-threshold summaries without another cohort VCF scan."""
    keys, tables, reports = [], {}, []
    for similarity in (.7, .9):
        output = outputs[similarity]
        completed = pilot.require_manifest(output)
        keys.append(completed["key"]["prepare"])
        reports.append(completed)
        with open(output / "changes.tsv") as handle:
            tables[similarity] = {(r["type"], r["metric"]): r for r in csv.DictReader(handle, delimiter="\t")}
    if keys[0] != keys[1]:
        raise ValueError("Threshold comparisons must use the same original input")
    rows = []
    for kind, metric in sorted(set(tables[.7]) | set(tables[.9])):
        a, b = tables[.7].get((kind, metric)), tables[.9].get((kind, metric))
        if a and b and float(a["before"]) != float(b["before"]):
            raise ValueError("Threshold baselines disagree")
        rows.append({"chrom": chrom, "type": kind, "metric": metric, "before": (a or b)["before"],
            "after_p70": a["after"] if a else None, "after_p90": b["after"] if b else None,
            "percent_reduction_p70": a["percent_reduction"] if a else None,
            "percent_reduction_p90": b["percent_reduction"] if b else None,
            "p90_minus_p70": float(b["after"]) - float(a["after"]) if a and b else None})
    output = root / chrom / "compare_p70_p90"
    output.mkdir(exist_ok=True)
    key = {"version": 1, "comparison_manifests": reports}
    if not pilot.valid_manifest(output, key):
        (output / "done.json").unlink(missing_ok=True)
        path = output / "threshold_comparison.tsv"
        pilot.tsv(path, rows, list(rows[0]))
        pilot.manifest(output, key, [path])
    return output / "threshold_comparison.tsv"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--chromosomes", nargs="+", default=[f"chr{i}" for i in range(1, 22)])
    parser.add_argument("--bcftools", default=os.environ.get("BCFTOOLS_BIN", "bcftools"))
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--similarities", nargs="+", type=float, choices=(.7, .9), default=[.7])
    args = parser.parse_args()
    if args.threads < 1 or any(c not in {f"chr{i}" for i in range(1, 22)} for c in args.chromosomes):
        parser.error("Use positive threads and chr1–chr21")
    root = args.root.resolve()
    for chrom in args.chromosomes:
        outputs = {}
        for similarity in args.similarities:
            tag = f"p{round(similarity * 100)}"
            if not (root / chrom / tag / "done.json").is_file():
                print(f"SKIP {chrom}: no completed {tag} result", flush=True)
                continue
            print(f"Comparing {chrom} {tag}; aggregate tables only", flush=True)
            output = compare(root, chrom, args.bcftools, args.threads, similarity)
            outputs[similarity] = output
            print(f"Results: {output}", flush=True)
            if len(args.similarities) == 1:
                print((output / "changes.tsv").read_text(), end="", flush=True)
        if .7 in outputs and .9 in outputs:
            path = side_by_side(root, chrom, outputs)
            print(f"Side-by-side comparison: {path}", flush=True)
            print(path.read_text(), end="", flush=True)


if __name__ == "__main__":
    main()
