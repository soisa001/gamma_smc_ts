"""Parallel interval preparation with sparse candidates and one panel merge.

Imported by phase_singletons. All positions and interval endpoints are 1-based,
inclusive. Halos supply evidence only; disjoint cores own output records.
"""
from __future__ import annotations

import argparse
from bisect import bisect_right
from concurrent.futures import ProcessPoolExecutor
import gzip
import hashlib
import io
import json
import logging
import os
from pathlib import Path
import re
import subprocess
import sys
import time

import phase_singletons as ph

VERSION = 1


def read_header(path, contig, bcftools):
    if Path(path).suffix == ".bcf":
        data = subprocess.check_output([bcftools, "view", "--no-version", "-h", str(path)], text=True)
        return ph.VCF(io.StringIO(data), contig)
    with ph.text_input(path, bcftools) as handle:
        return ph.VCF(handle, contig)


def partition(job, args, samples):
    lengths = []
    for name in ("scaffold", "hiphase"):
        header = read_header(job[name], job["contig"], args.bcftools)
        if header.samples != samples:
            raise ValueError("Sample identity/order differs from manifest")
        definitions = [line for line in header.meta if ph.header_id(line) == ("contig", job["contig"])]
        match = re.search(r"(?:<|,)length=([0-9]+)(?:,|>)", definitions[0]) if len(definitions) == 1 else None
        if not match or int(match[1]) < 1:
            raise ValueError("Chunking requires a positive contig length in both input headers")
        lengths.append(int(match[1]))
    if lengths[0] != lengths[1]:
        raise ValueError("Contig lengths differ between inputs")
    length = lengths[0]
    count = min(args.chunks_per_chrom, length)
    result = []
    for i in range(count):
        start, end = i * length // count + 1, (i + 1) * length // count
        result.append({**job, "chunk": i, "core": [start, end],
                       "region": [max(1, start - args.overlap_bp), min(length, end + args.overlap_bp)]})
    return result


def cache(marker, signature, directory):
    try:
        previous = json.loads(marker.read_text(encoding="utf-8"))
        if previous["signature"] == signature and all((directory / name).is_file() and ph.digest(directory / name) == sha for name, sha in previous["artifacts"].items()):
            return previous
    except (OSError, ValueError, KeyError):
        pass
    return None


def extract_worker(payload):
    job, root, params, samples, signature = payload
    directory = Path(root) / f"chunk_{job['chunk']:03d}"
    directory.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(filename=directory / "workflow.log", level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", force=True)
    signature = {**signature, "chunk_job": job}
    marker = directory / "extract.done.json"
    previous = cache(marker, signature, directory)
    if previous:
        return {"directory": str(directory), "cached": True, **previous}
    marker.unlink(missing_ok=True)
    logging.info("EXTRACT %s", json.dumps({"job": job, "parameters": params, "signature": signature}))
    stats = ph.prepare(job, directory, argparse.Namespace(**params), samples, sparse=True, defer_phase=True)
    stats["process_peak_memory_mib"] = ph.peak_memory_mib()
    files = ("raw_candidates.jsonl.gz", "anchors.sqlite", "duplicate_sites.tsv.gz", "merged.header.vcf", "header_merge.json")
    result = {"signature": signature, "stats": stats, "artifacts": {name: ph.digest(directory / name) for name in files}}
    ph.json_write(marker, result)
    logging.info("EXTRACT COMPLETE %s", json.dumps(stats))
    return {"directory": str(directory), "cached": False, **result}


def phase_worker(payload):
    job, directory, params, samples, signatures, indexes = payload
    directory = Path(directory)
    logging.basicConfig(filename=directory / "workflow.log", level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", force=True)
    signature = {"evidence": signatures, "chunk": job["chunk"], "min_anchors": params["min_anchors"]}
    marker = directory / "phase.done.json"
    previous = cache(marker, signature, directory)
    if previous:
        return {"cached": True, **previous}
    marker.unlink(missing_ok=True)
    started = time.monotonic()
    stats = dict(ph.resolve_candidates(job, directory, argparse.Namespace(**params), samples, sparse=True, other_indexes=indexes))
    stats["phase_seconds"] = time.monotonic() - started
    stats["process_peak_memory_mib"] = ph.peak_memory_mib()
    files = ("added_singletons.jsonl.gz", "excluded_singletons.tsv.gz")
    result = {"signature": signature, "stats": stats, "artifacts": {name: ph.digest(directory / name) for name in files}}
    ph.json_write(marker, result)
    logging.info("PHASE COMPLETE %s", json.dumps(stats))
    return {"cached": False, **result}


def sparse_records(directory, chunks, samples):
    """Expand one sparse record at a time, without a cohort-wide addition spool."""
    for job in chunks:
        with gzip.open(directory / f"chunk_{job['chunk']:03d}" / "added_singletons.jsonl.gz", "rt", encoding="utf-8") as handle:
            last = 0
            for line in handle:
                value = json.loads(line)
                fields, carrier = value["fields"], value["carrier"]
                pos = int(fields[1])
                if fields[0] != job["contig"] or not job["core"][0] <= pos <= job["core"][1] or pos < last:
                    raise ValueError("Sparse addition outside its core or out of order")
                if not 0 <= carrier < len(samples) or value["gt"] not in ("0|1", "1|0") or value["background"] not in ("0|0", "1|1"):
                    raise ValueError("Invalid sparse singleton genotype")
                last = pos
                cells = [value["background"]] * len(samples)
                cells[carrier] = value["gt"]
                row = fields + ["GT"] + cells
                yield row, "\t".join(row)


def merge(job, directory, args, samples, chunks, extracted, phased):
    temp = directory / "augmented.partial.vcf.gz"
    output = directory / "augmented.vcf.gz"
    hashes = [hashlib.sha256() for _ in chunks]
    counts = [0] * len(chunks)
    starts = [c["core"][0] for c in chunks]
    header = (Path(extracted[0]["directory"]) / "merged.header.vcf").read_text(encoding="utf-8")
    if any((Path(e["directory"]) / "merged.header.vcf").read_text(encoding="utf-8") != header for e in extracted):
        raise ValueError("Chunk headers differ")
    (directory / "merged.header.vcf").write_text(header, encoding="utf-8")
    additions = iter(sparse_records(directory, chunks, samples))
    with ph.text_input(job["scaffold"], args.bcftools, job["contig"]) as handle, gzip.open(temp, "wt", encoding="utf-8", newline="\n") as out:
        vcf = ph.VCF(handle, job["contig"])
        if vcf.samples != samples:
            raise ValueError("Scaffold sample order changed")
        out.write(header)
        addition = next(additions, None)
        for row, raw in vcf.records():
            pos = int(row[1])
            if pos > chunks[-1]["core"][1]:
                raise ValueError("Scaffold POS exceeds declared contig length")
            while addition is not None and int(addition[0][1]) < pos:
                out.write(addition[1] + "\n")
                addition = next(additions, None)
            out.write(raw + "\n")
            owner = bisect_right(starts, pos) - 1
            hashes[owner].update((raw + "\n").encode())
            counts[owner] += 1
        while addition is not None:
            out.write(addition[1] + "\n")
            addition = next(additions, None)
    for i, evidence in enumerate(extracted):
        if counts[i] != evidence["stats"].get("scaffold_records", 0) or hashes[i].hexdigest() != evidence["stats"]["scaffold_record_sha256"]:
            raise ValueError("Scaffold changed since chunk extraction")
    if not sum(counts):
        raise ValueError("Scaffold has no records for requested chromosome")
    audit = ph.audit_preservation(job["scaffold"], temp, job["contig"], args.bcftools)
    totals = {k: sum(p["stats"].get(k, 0) for p in phased) for k in ("physical", "excluded")}
    if audit["added_records"] != totals["physical"] or sum(totals.values()) != sum(e["stats"].get("singleton_candidates", 0) for e in extracted):
        raise ValueError("Chunk candidate/output count mismatch")
    # Compact audit files can be concatenated in core order; no full VCF spool.
    for name in ("added_singletons.jsonl.gz", "excluded_singletons.tsv.gz", "duplicate_sites.tsv.gz"):
        with gzip.open(directory / name, "wt", encoding="utf-8", newline="\n") as dest:
            for i, evidence in enumerate(extracted):
                with gzip.open(Path(evidence["directory"]) / name, "rt", encoding="utf-8") as source:
                    if i and name.endswith("tsv.gz"):
                        next(source)
                    for line in source:
                        dest.write(line)
    os.replace(temp, output)
    return {"final_counts": totals, "preservation": audit, "output": str(output)}


def job_run(payload):
    job, root, params, samples, manifest_hash = payload
    args = argparse.Namespace(**params)
    directory = Path(root) / job["contig"]
    directory.mkdir(parents=True, exist_ok=True)
    lock = directory / ".lock"
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        raise RuntimeError(f"Job locked: {lock}. Remove only after verifying no process uses this job.") from None
    try:
        with os.fdopen(fd, "w") as handle:
            handle.write(str(os.getpid()))
        logging.basicConfig(filename=directory / "workflow.log", level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", force=True)
        inputs = {name: ph.identity(job[name]) for name in ("scaffold", "hiphase")}
        signature = {"chunk_version": VERSION, "phase_version": ph.VERSION, "manifest_sha256": manifest_hash,
                     "inputs": inputs, "parameters": {k: v for k, v in params.items() if k not in ("workers", "memory_gb")}}
        marker = directory / "done.json"
        previous = cache(marker, signature, directory)
        if previous:
            return {"contig": job["contig"], "cached": True, **previous["result"]}
        marker.unlink(missing_ok=True)
        started = time.monotonic()
        logging.info("START %s", json.dumps({"signature": signature, "execution_parameters": params, "python": sys.version}))
        chunks = partition(job, args, samples)
        ph.json_write(directory / "chunks.json", chunks)
        # One pool across the two phases; total workers never multiplies by chromosomes.
        with ProcessPoolExecutor(max_workers=min(args.workers, len(chunks))) as pool:
            extracted = list(pool.map(extract_worker, [(c, str(directory), params, samples, signature) for c in chunks]))
            if not sum(e["stats"].get("hiphase_records", 0) for e in extracted):
                raise ValueError("HiPhase has no records for requested chromosome")
            evidence = [{"signature": e["signature"], "artifacts": e["artifacts"]} for e in extracted]
            phase_payloads = [(c, e["directory"], params, samples, evidence,
                               [str(Path(other["directory"]) / "anchors.sqlite") for j, other in enumerate(extracted) if j != i])
                              for i, (c, e) in enumerate(zip(chunks, extracted))]
            phased = list(pool.map(phase_worker, phase_payloads))
        result = {**merge(job, directory, args, samples, chunks, extracted, phased),
                  "chunks": [{"core": c["core"], "region": c["region"], "extract": e["stats"], "phase": p["stats"],
                              "extract_cached": e["cached"], "phase_cached": p["cached"]} for c, e, p in zip(chunks, extracted, phased)],
                  "elapsed_seconds": time.monotonic() - started, "merge_process_peak_memory_mib": ph.peak_memory_mib()}
        if inputs != {name: ph.identity(job[name]) for name in inputs}:
            raise ValueError("Input files changed during processing")
        ph.json_write(directory / "audit.json", {"signature": signature, "execution_parameters": params, **result})
        names = ("augmented.vcf.gz", "added_singletons.jsonl.gz", "excluded_singletons.tsv.gz", "duplicate_sites.tsv.gz", "merged.header.vcf", "chunks.json", "audit.json")
        artifacts = {name: ph.digest(directory / name) for name in names}
        for i, (e, p) in enumerate(zip(extracted, phased)):
            artifacts.update({f"chunk_{i:03d}/{name}": sha for name, sha in {**e["artifacts"], **p["artifacts"]}.items()})
        ph.json_write(marker, {"signature": signature, "result": result, "artifacts": artifacts})
        logging.info("COMPLETE %s", json.dumps(result))
        return {"contig": job["contig"], "cached": False, **result}
    except BaseException:
        logging.exception("Chunk workflow failed; no new valid completion marker")
        raise
    finally:
        lock.unlink(missing_ok=True)
