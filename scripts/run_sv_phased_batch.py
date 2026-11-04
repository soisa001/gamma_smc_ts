#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = ["google-crc32c==1.8.0"]
# ///
"""chr1–chr21 native Truvari outputs at 0.7/0.9, up to four chromosome jobs.

Unmodified --keep common --gt all clustering and genotype consolidation.
Run via uv run --with truvari==5.5.0. No uploads.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import gzip
import json
import logging
import os
from pathlib import Path
import subprocess
import sys
import time

import run_chr1_sv_merge as pilot

VERSION = 1
TEMPLATE = pilot.SOURCE.replace(".chr1.", ".{chrom}.")
THRESHOLDS = (.7, .9)
# Only distinguish structural alleles from the panel's small variants.
# No pure-allele, sequence-alphabet, FILTER, SV-type or upper-size restriction.
SV_FILTER = 'ILEN>=50 || ILEN<=-50 || TYPE~"bnd" || (ALT~"<" && ALT!="<NON_REF>" && ALT!="<NONREF>")'


def pipeline(args, commands, label):
    """Run a checked streaming pipeline; record data-bearing output only on disk."""
    processes = []
    with open(args.root / "logs" / f"{label}.log", "a") as log:
        try:
            for i, command in enumerate(commands):
                logging.info("PIPELINE %s", json.dumps(list(map(str, command))))
                previous = processes[-1].stdout if processes else None
                processes.append(subprocess.Popen(list(map(str, command)), stdin=previous,
                    stdout=subprocess.PIPE if i < len(commands) - 1 else log, stderr=log))
                if previous is not None:
                    previous.close()
            codes = [proc.wait() for proc in processes]
            if any(codes):
                raise RuntimeError(f"{label} failed {codes}; inspect {args.root / 'logs' / (label + '.log')}")
        finally:
            for proc in processes:
                if proc.poll() is None:
                    proc.terminate()
            for proc in processes:
                proc.wait()


def local_args(args, root, chrom):
    return pilot.parse_args(["--root", str(root), "--chrom", chrom, "--threads", "1",
        "--memory-gb", str(max(4, args.memory_gb // args.jobs)), "--bcftools", args.bcftools,
        "--source-uri", args.source_template.format(chrom=chrom), "--confirm-source-uri",
        args.source_template.format(chrom=chrom), "--gcloud", args.gcloud])


def check_vcf(args, path, label, allow_sort=False):
    try:
        pilot.command(args, [args.bcftools, "index", "-f", "-t", path], label + "-index")
    except RuntimeError:
        if not allow_sort:
            raise
        sorted_path = path.with_name(path.name + ".sorted.gz")
        pilot.command(args, [args.bcftools, "sort", "-m", f"{max(1, args.memory_gb // 2)}G",
                      "-T", args.root / (label + "-sort"), "-Oz", "-o", sorted_path, path], label + "-sort")
        sorted_path.replace(path)
        pilot.command(args, [args.bcftools, "index", "-f", "-t", path], label + "-index-sorted")
    with gzip.open(path, "rb") as handle:
        for _ in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            pass
    return int(pilot.command(args, [args.bcftools, "index", "-n", path], label + "-count", True).strip())


def prepare(args, chrom, root):
    pa = local_args(args, root, chrom)
    directory = root / "prepare"
    directory.mkdir(exist_ok=True)
    uri = args.source_template.format(chrom=chrom)
    source = Path(args.input_template.format(chrom=chrom)).resolve() if args.input_template else None
    origin = {"path": str(source), "sha256": pilot.common.digest(source)} if source else {"uri": uri}
    settings = {"version": VERSION, "origin": origin, "chrom": chrom,
                "sv_filter": SV_FILTER, "split_multiallelic": True, "multi_overlaps": "0", "genotypes": True}
    cached = pilot.valid_manifest(directory)
    if cached and cached["key"]["settings"] == settings:
        logging.info("Reuse verified phased SV input; source snapshot recorded in manifest")
        return cached
    (directory / "done.json").unlink(missing_ok=True)
    if source is None:
        pilot.stage(pa)
        provenance = pilot.require_manifest(root / "inputs")
        source = root / "inputs" / uri.rsplit("/", 1)[1]
    else:
        provenance = origin
    partial = directory / "input.partial.vcf.gz"
    # Streaming conversion/decomposition only; no padded-allele rewriting or genotype parsing in Python.
    pipeline(pa, [
        [args.bcftools, "view", "--no-version", "-r", chrom, "-i", SV_FILTER, "-Ou", source],
        [args.bcftools, "norm", "--no-version", "-m", "-any", "--multi-overlaps", "0", "-Ou", "-"],
        [args.bcftools, "view", "--no-version", "-i", SV_FILTER, "-Oz", "-o", partial, "-"],
    ], "prepare-input")
    samples = pilot.command(pa, [args.bcftools, "query", "-l", partial], "sample-manifest", True)
    if not samples.strip():
        raise ValueError("--keep common / --gt all require the original cohort genotypes, not sites-only input")
    (directory / "samples.txt").write_text(samples)
    count = check_vcf(pa, partial, "input")
    if count == 0:
        raise ValueError(f"No SVs found on {chrom}")
    target = directory / "input.vcf.gz"
    partial.replace(target)
    Path(str(partial) + ".tbi").replace(Path(str(target) + ".tbi"))
    version = pilot.command(pa, [args.bcftools, "--version"], "bcftools-version", True).splitlines()[0]
    result = pilot.manifest(directory, {"settings": settings, "source": provenance,
        "bcftools": version, "records": count}, [target, Path(str(target) + ".tbi"), directory / "samples.txt"])
    if not args.input_template and not args.keep_staged_bcf:
        staged = (root / "inputs").resolve()
        for name in (uri.rsplit("/", 1)[1], uri.rsplit("/", 1)[1] + ".csi"):
            path = staged / name
            if path.resolve().parent != staged:
                raise ValueError("Unexpected staging path; refusing cleanup")
            path.unlink(missing_ok=True)
        (staged / "done.json").unlink(missing_ok=True)
    return result


def merge(args, chrom, root, prepared):
    rows = []
    for similarity in THRESHOLDS:
        directory = root / f"p{round(similarity * 100)}"
        (directory / "logs").mkdir(parents=True, exist_ok=True)
        pa = local_args(args, directory, chrom)
        key = {"version": VERSION, "prepare": prepared, "similarity": similarity,
               "refdist": 500, "keep": "common", "gt": "all", "sizemin": 50, "sizemax": -1,
               "chain": False, "typeignore": False, "hap": False, "truvari": "5.5.0",
               "genotype_consolidation": "native"}
        if pilot.valid_manifest(directory, key):
            logging.info("Reuse %s %.1f native merged VCF", chrom, similarity)
            rows.append(json.loads((directory / "summary.json").read_text())["counts"])
            continue
        (directory / "done.json").unlink(missing_ok=True)
        kept = directory / f"{chrom}.partial.vcf.gz"
        removed = directory / f"{chrom}.collapsed.partial.vcf.gz"
        cmd = ["truvari", "collapse", "-i", root / "prepare/input.vcf.gz", "-o", kept, "-c", removed,
               "-p", similarity, "-P", similarity, "-r", "500", "-s", "50", "-S", "-1",
               "--keep", "common", "--gt", "all"]
        pilot.command(pa, cmd, "collapse-native")
        kept_count = check_vcf(pa, kept, "kept", allow_sort=True)
        removed_count = check_vcf(pa, removed, "collapsed", allow_sort=True)
        if kept_count + removed_count != prepared["key"]["records"]:
            raise ValueError("Input count differs from retained + collapsed counts")
        sample_text = (root / "prepare/samples.txt").read_text()
        for name, path in (("kept", kept), ("collapsed", removed)):
            if pilot.command(pa, [args.bcftools, "query", "-l", path], name + "-samples", True) != sample_text:
                raise ValueError("Output sample manifest changed")
        artifacts = []
        for partial, final in ((kept, directory / f"{chrom}.vcf.gz"),
                               (removed, directory / f"{chrom}.collapsed.vcf.gz")):
            partial.replace(final)
            Path(str(partial) + ".tbi").replace(Path(str(final) + ".tbi"))
            artifacts.extend([final, Path(str(final) + ".tbi")])
        counts = {"chrom": chrom, "similarity": similarity, "before": prepared["key"]["records"],
                  "after": kept_count, "collapsed": removed_count}
        pilot.common.json_write(directory / "summary.json", {"counts": counts, "parameters": key,
            "phase_policy": "Unmodified Truvari 5.5.0 genotype consolidation; no custom GT or phase transfer."})
        artifacts.append(directory / "summary.json")
        pilot.manifest(directory, key, artifacts)
        rows.append(counts)
    pilot.tsv(root / "counts.tsv", rows, list(rows[0]))


def worker(args, chrom):
    root = args.root / chrom
    (root / "logs").mkdir(parents=True, exist_ok=True)
    logging.basicConfig(filename=root / "logs/batch.log", level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s", force=True)
    lock = root / ".phased_batch.lock"
    fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    os.write(fd, f"{os.getpid()}\n".encode())
    os.close(fd)
    try:
        logging.info("PARAMETERS %s", json.dumps(vars(args), default=str, sort_keys=True))
        logging.info("PROVENANCE seed=%s script_sha256=%s python=%s", pilot.SEED, pilot.common.digest(__file__), sys.version)
        prepared = prepare(args, chrom, root) if args.phase in {"prepare", "all"} else pilot.require_manifest(root / "prepare")
        if args.phase in {"merge", "all"}:
            merge(args, chrom, root, prepared)
    except Exception:
        logging.exception("Chromosome workflow failed")
        raise
    finally:
        lock.unlink(missing_ok=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--phase", choices=("plan", "prepare", "merge", "all"), default="plan")
    parser.add_argument("--chromosomes", nargs="+", default=[f"chr{i}" for i in range(1, 22)])
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--memory-gb", type=int, default=26, help="Total budget divided between jobs; not an OS cap")
    parser.add_argument("--source-template", default=TEMPLATE)
    parser.add_argument("--input-template", help="Optional local indexed BCF/VCF template containing {chrom}")
    parser.add_argument("--keep-staged-bcf", action="store_true")
    parser.add_argument("--bcftools", default=os.environ.get("BCFTOOLS_BIN", "bcftools"))
    parser.add_argument("--gcloud", default="gcloud")
    parser.add_argument("--worker", help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    args.root = args.root.resolve()
    if not 1 <= args.jobs <= 4 or args.memory_gb < 4 * args.jobs:
        parser.error("Use 1–4 jobs and at least 4 GB per job")
    if not set(args.chromosomes) <= {f"chr{i}" for i in range(1, 22)} or len(args.chromosomes) != len(set(args.chromosomes)):
        parser.error("Chromosomes must be a unique subset of chr1–chr21")
    if "{chrom}" not in args.source_template or (args.input_template and "{chrom}" not in args.input_template):
        parser.error("Input templates must contain {chrom}")
    if args.worker:
        if args.worker not in args.chromosomes:
            parser.error("Worker chromosome not requested")
        return worker(args, args.worker)
    if args.phase == "plan":
        print(json.dumps({"parameters": vars(args), "thresholds": THRESHOLDS, "keep": "common", "gt": "all",
            "refdist": 500, "hap": False, "sample_genotypes": True, "genotype_consolidation": "native", "uploads": False,
            "sources": [args.source_template.format(chrom=c) for c in args.chromosomes]}, default=str, indent=2))
        return
    (args.root / "logs").mkdir(parents=True, exist_ok=True)
    pilot.common.json_write(args.root / f"run-{time.time_ns()}.json", {
        "parameters": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
        "thresholds": THRESHOLDS, "seed": pilot.SEED, "script_sha256": pilot.common.digest(__file__)})
    forwarded = list(argv if argv is not None else sys.argv[1:])
    def run_one(chrom):
        print(f"Starting {chrom}; native outputs at similarity 0.7 and 0.9.", flush=True)
        with open(args.root / "logs" / f"{chrom}.worker.log", "a") as log:
            result = subprocess.run([sys.executable, __file__, *forwarded, "--worker", chrom], stdout=log, stderr=log)
        if result.returncode:
            raise RuntimeError(f"{chrom} failed; inspect logs/{chrom}.worker.log and {chrom}/logs/batch.log")
        print(f"Completed {chrom}.", flush=True)
    failures = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(run_one, chrom) for chrom in args.chromosomes]
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as error:
                failures.append(str(error))
    if failures:
        raise SystemExit("\n".join(failures))
    if args.phase in {"merge", "all"}:
        rows = []
        for chrom in args.chromosomes:
            with open(args.root / chrom / "counts.tsv") as handle:
                rows.extend(pilot.csv.DictReader(handle, delimiter="\t"))
        pilot.tsv(args.root / "counts.tsv", rows, list(rows[0]))
    print(f"Completed {args.phase}; outputs remain in {args.root}")


if __name__ == "__main__":
    main()
