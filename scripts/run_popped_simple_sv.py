#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = ["truvari==5.5.0", "google-crc32c==1.8.0"]
# ///
"""Merge only trimmed pure indels from the popped panel; all outputs stay local."""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from contextlib import ExitStack
import copy
import gzip
import importlib.metadata
import json
import logging
import os
from pathlib import Path
import sys
import time

import pysam

import inspect_sv_representation as representation
import run_chr1_sv_merge as pilot
import run_sv_phased_batch as batch

VERSION = 1
KINDS = ("INS", "DEL")
# No length threshold on merge candidates. Equal-length literal alleles cannot
# be pure indels or pass the final net-length cutoff, so avoid decoding their GTs.
EXTRACT = 'ILEN!=0 || ALT~"<" || TYPE~"bnd"'


def classify(pos, ref, alt):
    """Return class, net size and minimally anchored alleles; never decompose replacements."""
    if "," in alt:
        raise ValueError("Expected the confirmed biallelic popped panel; multiallelic input is unsupported")
    desc, excluded = representation.describe(ref, alt)
    if excluded:
        return "OTHER", 0, pos, ref, alt, False
    size, kind = desc["net_bp"], desc["type"]
    simple = kind in KINDS and not desc["complex"] and not desc["ambiguous"]
    if not simple:
        return kind, size, pos, ref, alt, False
    prefix, suffix = desc["shared_prefix_bp"], desc["shared_suffix_bp"]
    # Keep one shared base as the VCF anchor. A right anchor already in the
    # source is retained when there is no shared prefix; no reference is fetched.
    start = prefix - 1 if prefix else 0
    end_ref = len(ref) - suffix + (0 if prefix else 1)
    end_alt = len(alt) - suffix + (0 if prefix else 1)
    return kind, size, pos + start, ref[start:end_ref], alt[start:end_alt], True


def final_allele(ref, alt):
    desc, excluded = representation.describe(ref, alt)
    if desc:
        return desc["type"] in KINDS and desc["net_bp"] >= 50
    # Symbolic SVs/BNDs have no literal net length: preserve them unchanged.
    return excluded == "symbolic_or_breakend"


def bump(counts, kind, size):
    if kind in KINDS and size >= 50:
        counts[kind] += 1
        counts[kind + "_bp"] += size


def finish(args, partial, final, label, sort=False):
    count = batch.check_vcf(args, partial, label, allow_sort=sort)
    partial.replace(final)
    Path(str(partial) + ".tbi").replace(Path(str(final) + ".tbi"))
    return count


def artifacts(paths):
    return [p for path in paths for p in (path, Path(str(path) + ".tbi"))]


def prepare(args):
    directory = args.root / "prepare"
    directory.mkdir(exist_ok=True)
    source = args.input_bcf or args.root / "inputs" / args.source_uri.rsplit("/", 1)[1]
    if not args.input_bcf:
        pilot.stage(args)
    source = source.resolve()
    source_snapshot = pilot.common.identity(source)
    settings = {"version": VERSION, "source": {"path": str(source), "sha256": pilot.common.digest(source)},
                "chrom": args.chrom, "extract": EXTRACT, "final_min_net_bp": 50,
                "merge": "trimmed_pure_indels_separate_types", "no_input_size_filter": True}
    cached = pilot.valid_manifest(directory)
    if cached and cached["key"]["settings"] == settings:
        logging.info("Reuse verified simple-indel preparation")
        return cached
    (directory / "done.json").unlink(missing_ok=True)
    samples = pilot.command(args, [args.bcftools, "query", "-l", source], "source-samples", True)
    if not samples.strip():
        raise ValueError("Original sample genotypes are required for --keep common --gt all")
    (directory / "samples.txt").write_text(samples)
    counts = {k: Counter() for k in ("before", "simple", "passthrough")}
    records, trimmed, written = 0, 0, Counter()
    parts = (*KINDS, "passthrough.sv50")
    with ExitStack() as stack:
        handles = {k: stack.enter_context(pysam.BGZFile(str(directory / f"{k}.partial.vcf.gz"), "w")) for k in parts}
        mapping = stack.enter_context(gzip.open(directory / "trimmed_alleles.tsv.gz", "wt"))
        mapping.write("chrom\tid\toriginal_pos\toriginal_ref\toriginal_alt\ttrimmed_pos\ttrimmed_ref\ttrimmed_alt\n")
        lines = stack.enter_context(pilot.stream(args, [args.bcftools, "view", "--no-version", "--threads",
            args.threads, "-r", args.chrom, "-i", EXTRACT, "-Ov", source], "partition"))
        for line in lines:
            if line.startswith("#"):
                for handle in handles.values():
                    handle.write(line.encode())
                continue
            records += 1
            # Leave the entire FORMAT/sample tail opaque and unchanged.
            fields = line.rstrip("\n").split("\t", 8)
            pos, ref, alt = int(fields[1]), fields[3], fields[4]
            kind, size, newpos, newref, newalt, simple = classify(pos, ref, alt)
            bump(counts["before"], kind, size)
            if simple:
                info = pilot.info_parse(fields[7])
                if ("SVTYPE" in info and info["SVTYPE"] != kind) or ("SVLEN" in info and abs(int(info["SVLEN"])) != size):
                    raise ValueError("A simple allele has inconsistent INFO/SVTYPE or SVLEN; inspect input metadata on VM")
                bump(counts["simple"], kind, size)
                if (pos, ref, alt) != (newpos, newref, newalt):
                    trimmed += 1
                    mapping.write("\t".join(map(str, (fields[0], fields[2], pos, ref, alt, newpos, newref, newalt))) + "\n")
                    fields[1], fields[3], fields[4] = str(newpos), newref, newalt
                    if "END" in info:
                        info["END"] = str(newpos + len(newref) - 1)
                        fields[7] = pilot.info_string(info)
                handles[kind].write(("\t".join(fields) + "\n").encode())
                written[kind] += 1
            elif final_allele(ref, alt):
                bump(counts["passthrough"], kind, size)
                handles["passthrough.sv50"].write(line.encode())
                written["passthrough.sv50"] += 1
            if records % 100000 == 0:
                logging.info("PARTITION_PROGRESS indel_or_symbolic_records=%d written=%s", records, dict(written))
    if records == 0:
        raise ValueError("No indel candidates found on chr1; verify the input chromosome and index")
    paths = []
    for kind in parts:
        target = directory / f"{kind}.vcf.gz"
        n = finish(args, directory / f"{kind}.partial.vcf.gz", target, "prepare-" + kind, sort=True)
        if n != written[kind]:
            raise ValueError("Partition count changed during compression/sorting")
        paths.append(target)
    report = {"counts_sv50": counts, "partition_records_all_sizes": dict(written),
              "scanned_indel_or_symbolic_records": records, "trimmed_records": trimmed,
              "input": str(source), "samples": len(samples.splitlines())}
    if pilot.common.identity(source) != source_snapshot:
        raise ValueError("Source changed during preparation; rerun from a stable indexed input")
    pilot.common.json_write(directory / "summary.json", report)
    result = pilot.manifest(directory, {"settings": settings, **report}, artifacts(paths) +
        [directory / "samples.txt", directory / "summary.json", directory / "trimmed_alleles.tsv.gz"])
    print("Prepared simple indels; aggregate counts: " + json.dumps(report["counts_sv50"]), flush=True)
    return result


def site_counts(args, path, label):
    counts = Counter()
    with pilot.stream(args, [args.bcftools, "query", "-f", "%REF\t%ALT\n", path], label) as lines:
        for line in lines:
            ref, alt = line.rstrip("\n").split("\t")
            desc, _ = representation.describe(ref, alt)
            if desc:
                bump(counts, desc["type"], desc["net_bp"])
    return counts


def collapse(args, prepared, similarity, kind):
    directory = args.root / f"p{round(similarity * 100)}" / kind
    directory.mkdir(parents=True, exist_ok=True)
    key = {"version": VERSION, "prepare": prepared, "similarity": similarity, "type": kind,
           "truvari": "5.5.0", "keep": "common", "gt": "all", "refdist": 500, "pctovl": 0,
           "sizemin": 0, "sizemax": -1, "hap": False, "chain": False, "genotypes": "native"}
    if pilot.valid_manifest(directory, key):
        return json.loads((directory / "summary.json").read_text())
    (directory / "done.json").unlink(missing_ok=True)
    kept, removed = [directory / f"{name}.partial.vcf.gz" for name in ("retained", "collapsed")]
    pilot.command(args, ["truvari", "collapse", "-i", args.root / "prepare" / f"{kind}.vcf.gz",
        "-o", kept, "-c", removed, "-p", similarity, "-P", similarity, "-r", 500, "-O", 0,
        "-s", 0, "-S", -1, "--keep", "common", "--gt", "all"], f"p{round(similarity*100)}-{kind}-collapse")
    paths = [directory / f"{name}.vcf.gz" for name in ("retained", "collapsed")]
    ns = [finish(args, temp, path, directory.parent.name + "-" + kind + "-" + path.stem, sort=True)
          for temp, path in zip((kept, removed), paths)]
    if sum(ns) != prepared["key"]["partition_records_all_sizes"].get(kind, 0):
        raise ValueError("Native retained + collapsed counts differ from prepared partition")
    samples = (args.root / "prepare/samples.txt").read_text()
    for path in paths:
        if pilot.command(args, [args.bcftools, "query", "-l", path], directory.parent.name + "-" + kind + "-" + path.stem + "-samples", True) != samples:
            raise ValueError("Output sample order changed")
    report = {"all_sizes": {"retained": ns[0], "collapsed": ns[1]},
              "retained_sv50": site_counts(args, paths[0], directory.parent.name + "-" + kind + "-retained-counts"),
              "collapsed_sv50": site_counts(args, paths[1], directory.parent.name + "-" + kind + "-collapsed-counts")}
    expected = prepared["key"]["counts_sv50"]["simple"]
    for metric in (kind, kind + "_bp"):
        if report["retained_sv50"][metric] + report["collapsed_sv50"][metric] != expected.get(metric, 0):
            raise ValueError("Native retained + collapsed SV50 count/bp differs from prepared simple alleles")
    pilot.common.json_write(directory / "summary.json", report)
    pilot.manifest(directory, key, artifacts(paths) + [directory / "summary.json"])
    return report


def merge(args, prepared):
    rows = []
    for similarity in args.similarities:
        tag = f"p{round(similarity * 100)}"
        print(f"Starting {tag}: native collapse separately for simple INS and DEL.", flush=True)
        jobs = min(2, args.threads)
        part_args = copy.copy(args)
        part_args.memory_gb = max(4, args.memory_gb // jobs)
        part_args.threads = max(1, args.threads // jobs)
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            reports = dict(zip(KINDS, pool.map(lambda kind: collapse(part_args, prepared, similarity, kind), KINDS)))
        directory = args.root / tag
        key = {"version": VERSION, "prepare": prepared, "parts": reports, "similarity": similarity,
               "final_min_net_bp": 50, "symbolic": "unchanged"}
        before = prepared["key"]["counts_sv50"]
        after = Counter(before["passthrough"])
        for report in reports.values():
            after.update(report["retained_sv50"])
        for scope, b, a in (("all_INS_DEL", before["before"], after),
                            ("simple_only", before["simple"], sum((Counter(r["retained_sv50"]) for r in reports.values()), Counter())),
                            ("unchanged_passthrough", before["passthrough"], before["passthrough"])):
            for kind in (*KINDS, "TOTAL"):
                kinds = KINDS if kind == "TOTAL" else (kind,)
                n0, n1 = (sum(c.get(k, 0) for k in kinds) for c in (b, a))
                rows.append({"chrom": args.chrom, "similarity": similarity, "scope": scope, "type": kind,
                    "before_sv50": n0, "after_sv50": n1, "merged_out": n0 - n1,
                    "percent_merged_out": 100 * (n0 - n1) / n0 if n0 else 0,
                    "before_raw_bp": sum(b.get(k + "_bp", 0) for k in kinds),
                    "after_raw_bp": sum(a.get(k + "_bp", 0) for k in kinds)})
        if not pilot.valid_manifest(directory, key):
            (directory / "done.json").unlink(missing_ok=True)
            # -a interleaves sorted partitions; no giant uncompressed intermediate.
            partial = directory / f"{args.chrom}.sv50.partial.vcf.gz"
            batch.pipeline(args, [
                [args.bcftools, "concat", "--no-version", "-a", "-Ou",
                 directory / "INS/retained.vcf.gz", directory / "DEL/retained.vcf.gz",
                 args.root / "prepare/passthrough.sv50.vcf.gz"],
                [args.bcftools, "view", "--no-version", "-i", batch.FINAL_FILTER + ' || ALT~"<" || TYPE~"bnd"',
                 "--threads", args.threads, "-Oz", "-o", partial, "-"],
            ], tag + "-combine")
            target = directory / f"{args.chrom}.sv50.vcf.gz"
            total = finish(args, partial, target, tag + "-final", sort=True)
            expected = (sum(sum(r["retained_sv50"].get(k, 0) for k in KINDS) for r in reports.values())
                        + prepared["key"]["partition_records_all_sizes"].get("passthrough.sv50", 0))
            if total != expected or site_counts(args, target, tag + "-final-counts") != after:
                raise ValueError("Combined final counts differ from partition totals")
            if pilot.command(args, [args.bcftools, "query", "-l", target], tag + "-final-samples", True) != (args.root / "prepare/samples.txt").read_text():
                raise ValueError("Combined sample order changed")
            pilot.manifest(directory, key, artifacts([target]))
        print(f"Completed {tag}.", flush=True)
    pilot.tsv(args.root / "comparison.tsv", rows, list(rows[0]))
    print((args.root / "comparison.tsv").read_text(), end="")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--input-bcf", type=Path, help="Existing indexed popped panel; never removed")
    parser.add_argument("--source-uri", default=representation.SOURCE)
    parser.add_argument("--chrom", choices=("chr1",), default="chr1")
    parser.add_argument("--similarities", nargs="+", type=float, choices=(.7, .9), default=[.7, .9])
    parser.add_argument("--phase", choices=("prepare", "merge", "all"), default="all")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--memory-gb", type=int, default=26)
    parser.add_argument("--bcftools", default=os.environ.get("BCFTOOLS_BIN", "bcftools"))
    parser.add_argument("--billing-project", default=None)
    opts = parser.parse_args(argv)
    if len(set(opts.similarities)) != len(opts.similarities):
        parser.error("Similarity thresholds must be unique")
    args = pilot.parse_args(["--root", str(opts.root), "--source-uri", opts.source_uri,
        "--confirm-source-uri", opts.source_uri, "--threads", str(opts.threads),
        "--memory-gb", str(opts.memory_gb), "--bcftools", opts.bcftools])
    args.input_bcf = opts.input_bcf.resolve() if opts.input_bcf else None
    args.similarities, args.phase, args.billing_project = opts.similarities, opts.phase, opts.billing_project
    (args.root / "logs").mkdir(parents=True, exist_ok=True)
    logging.basicConfig(filename=args.root / "logs/workflow.log", level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s", force=True)
    import fcntl
    with open(args.root / ".simple_merge.lock", "a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise SystemExit("Another simple merge is running against this root") from None
        provenance = {"parameters": vars(args), "seed": pilot.SEED, "script_sha256": pilot.common.digest(__file__),
                      "python": sys.version, "truvari": importlib.metadata.version("truvari"),
                      "bcftools": pilot.command(args, [args.bcftools, "--version"], "bcftools-version", True)}
        if provenance["truvari"] != "5.5.0":
            raise ValueError("Use the pinned uv environment with Truvari 5.5.0")
        logging.info("PROVENANCE %s", json.dumps(provenance, default=str))
        (args.root / "logs" / f"run-{time.time_ns()}.json").write_text(json.dumps(provenance, default=str, indent=2))
        try:
            prepared = prepare(args) if args.phase in {"prepare", "all"} else pilot.require_manifest(args.root / "prepare")
            if prepared["key"]["settings"]["version"] != VERSION:
                raise ValueError("Preparation version changed; run prepare again")
            if args.phase in {"merge", "all"}:
                merge(args, prepared)
        except Exception:
            logging.exception("Simple SV workflow failed")
            raise SystemExit(f"Failed; inspect {args.root / 'logs/workflow.log'} on the VM") from None


if __name__ == "__main__":
    main()
