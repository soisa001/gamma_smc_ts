#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = ["google-crc32c==1.8.0", "matplotlib==3.10.8"]
# ///
"""Phased, local-output Truvari pilot. See docs/chr1_sv_merge_workbench.md.

Truvari 5.5.0 is installed separately (Linux); add --with truvari==5.5.0 to
uv run for the merge/all phase. No upload operation exists in this program.
"""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager, ExitStack
import csv
import gzip
import hashlib
import importlib.metadata
import json
import logging
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

import phase_singletons as common
import run_singleton_chr as storage

SCHEMA = 1  # Bump for changes to analysis semantics, not documentation/commits.
SEED = 1729  # Existing repository seed; this workflow does not sample.
SOURCE = ("gs://rw-long-reads-transfer-2026-06-17/v9/lrWGS/panel/panel/"
          "panel_bubble_split_vcf/aou_lr_phase2_v1.chr1.bubble.split.bcf")
OLD_TRUVARI = set("NumCollapsed NumConsolidated CollapseId CollapseStart CollapseEnd CollapseSize "
                  "TruScore PctSeqSimilarity PctSizeSimilarity PctRecOverlap StartDistance "
                  "EndDistance SizeDiff GTMatch MatchId Multi".split())
REPLACED = OLD_TRUVARI | {"AC", "AN", "AF", "SVLEN", "SVTYPE"}
HEADERS = [
    '##INFO=<ID=PM_ROW,Number=1,Type=Integer,Description="Source candidate row before ALT splitting">',
    '##INFO=<ID=PM_STATUS,Number=1,Type=String,Description="Allele representation diagnostic; not merge eligibility under native policy">',
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Pilot allele SV type">',
    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Signed net allele size; unknown omitted">',
    '##INFO=<ID=AC,Number=A,Type=Integer,Description="ALT count recomputed from GT">',
    '##INFO=<ID=AN,Number=1,Type=Integer,Description="Called allele count including partial GT calls">',
    '##INFO=<ID=AF,Number=A,Type=Float,Description="AC divided by AN; missing when AN is zero">',
]


def info_parse(value):
    return dict(token.split("=", 1) if "=" in token else (token, True)
                for token in value.split(";") if token and token != ".")


def info_string(info):
    return ";".join(k if v is True else f"{k}={v}" for k, v in info.items()) or "."


def open_text(path, mode="rt"):
    return (gzip.open(path, mode, encoding="utf-8", compresslevel=1) if str(path).endswith(".gz")
            else open(path, mode, encoding="utf-8"))


def command(args, cmd, label, capture=False, env=None):
    """Keep tool stderr and data-bearing stdout on disk, not notebook output."""
    cmd = [str(x) for x in cmd]
    logging.info("COMMAND %s", json.dumps(cmd))
    with open(args.root / "logs" / f"{label}.log", "a", encoding="utf-8") as log:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE if capture else log,
                              stderr=log, text=True, env={**os.environ, **(env or {})})
        if capture and proc.stdout:
            log.write(proc.stdout)
    if proc.returncode:
        raise RuntimeError(f"{label} failed; inspect logs/{label}.log on the VM")
    return proc.stdout if capture else None


@contextmanager
def stream(args, cmd, label):
    logging.info("STREAM %s", json.dumps([str(x) for x in cmd]))
    with open(args.root / "logs" / f"{label}.log", "a", encoding="utf-8") as log:
        proc = subprocess.Popen([str(x) for x in cmd], stdout=subprocess.PIPE, stderr=log,
                                text=True, encoding="utf-8")
        try:
            yield proc.stdout
            proc.stdout.close()
            if proc.wait():
                raise RuntimeError(f"{label} failed; inspect logs/{label}.log on the VM")
        finally:
            if proc.poll() is None:
                proc.terminate()
                proc.wait()


def manifest(directory, key, paths):
    data = {"schema": SCHEMA, "key": key,
            "artifacts": {str(p.relative_to(directory)): common.digest(p) for p in paths}}
    common.json_write(directory / "done.json", data)
    return data


def valid_manifest(directory, key=None):
    try:
        data = json.loads((directory / "done.json").read_text())
        if data["schema"] != SCHEMA or (key is not None and data["key"] != key):
            return None
        if not data["artifacts"]:
            return None
        for name, digest in data["artifacts"].items():
            if common.digest(directory / name) != digest:
                return None
        return data
    except (OSError, ValueError, KeyError):
        return None


def require_manifest(directory):
    data = valid_manifest(directory)
    if data is None:
        raise ValueError(f"Missing/corrupt/stale phase artifacts; rerun phase: {directory.name}")
    return data


@contextmanager
def rsync_lock():
    """Serialize this workflow's rsync calls sharing gcloud's listing cache."""
    import fcntl  # Staging runs on the Linux Workbench/HPC host.
    cache = Path(os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache")
    path = cache / "gamma_smc_ts" / "gcloud-rsync.lock"
    path.parent.mkdir(parents=True, exist_ok=True)
    # Keep the inode: unlinking a lock file can let waiters lock different files.
    with path.open("a") as handle:
        started = time.monotonic()
        logging.info("WAIT_RSYNC_LOCK %s", path)
        fcntl.flock(handle, fcntl.LOCK_EX)
        try:
            logging.info("ACQUIRED_RSYNC_LOCK %s wait_seconds=%.3f", path, time.monotonic() - started)
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)
            logging.info("RELEASED_RSYNC_LOCK %s", path)


def stage(args):
    if args.confirm_source_uri != args.source_uri:
        raise ValueError("Staging requires --confirm-source-uri equal to --source-uri")
    directory = args.root / "inputs"
    directory.mkdir(exist_ok=True)
    remote = storage.metadata_for(args, args.source_uri)
    key = {"source": remote}
    if valid_manifest(directory, key):
        logging.info("Reuse checksum-verified staged input")
        return
    (directory / "done.json").unlink(missing_ok=True)
    need = sum(x["size"] for x in remote.values() if x) + args.scratch_reserve_gb * 2**30
    if shutil.disk_usage(directory).free < need:
        raise ValueError("Insufficient free disk for input replacement plus scratch reserve")
    with rsync_lock():
        if shutil.disk_usage(directory).free < need:
            raise ValueError("Insufficient free disk after waiting for staging lock")
        command(args, storage.sync_command(args.gcloud, args.source_uri, directory, args.billing_project),
                "rsync", env={"CLOUDSDK_STORAGE_PROCESS_COUNT": "1",
                              "CLOUDSDK_STORAGE_THREAD_COUNT": str(args.threads)})
    path = directory / args.source_uri.rsplit("/", 1)[1]
    storage.verify_object(path, remote["bcf"])
    index = Path(str(path) + ".csi")
    if remote["csi"]:
        storage.verify_object(index, remote["csi"])
    else:
        command(args, [args.bcftools, "index", "-f", "-c", "--threads", args.threads, path], "stage-index")
    if storage.metadata_for(args, args.source_uri) != remote:
        raise ValueError("Source generation changed during staging; rerun stage")
    count = command(args, [args.bcftools, "index", "--nrecords", path], "stage-check", True).strip()
    if not count.isdigit():
        raise ValueError("Index did not return a record count")
    manifest(directory, key, [path, index])


def allele(ref, alt, pos, info, minimum=50, maximum=100000):
    """Return type, size, changed-reference half-open interval, merge status.

    Literal alleles use net length (not padded REF length); BND has unknown size.
    Complex replacements and unresolved calls remain in both catalogues.
    """
    if alt in {".", "*", "<NON_REF>", "<NONREF>"}:
        return None
    start = pos - 1
    if "[" in alt or "]" in alt or (len(alt) > 1 and (alt.startswith(".") or alt.endswith("."))):
        return "BND", None, start, start, "breakend"
    if alt.startswith("<"):
        kind = alt.strip("<>").split(":")[0]
        raw_size = info.get("SVLEN", ".").split(",")[0]
        size = abs(int(raw_size)) if raw_size not in {".", ""} else None
        end = int(info.get("END", pos))
        if size is None and kind in {"DEL", "DUP", "INV", "CNV"} and end > pos:
            size = end - pos
        if size is not None and size < minimum:
            return None
        return kind, size, start, max(start, end) if kind != "INS" else start, "symbolic"
    delta = len(alt) - len(ref)
    if abs(delta) < minimum:
        return None
    kind, size = ("INS" if delta > 0 else "DEL"), abs(delta)
    prefix = 0
    while prefix < min(len(ref), len(alt)) and ref[prefix] == alt[prefix]:
        prefix += 1
    r, a = ref[prefix:], alt[prefix:]
    suffix = 0
    while suffix < min(len(r), len(a)) and r[-1 - suffix] == a[-1 - suffix]:
        suffix += 1
    rlen, alen = len(r) - suffix, len(a) - suffix
    status = "eligible"
    if set((ref + alt).upper()) - set("ACGT"):
        status = "ambiguous_sequence"
    elif rlen and alen:
        status = "complex_replacement"
    elif maximum != -1 and size > maximum:
        status = "above_sizemax"
    return kind, size, start + prefix, start + prefix + rlen, status


def genotypes(fields):
    n = max(0, len(fields) - 9)
    fmt = fields[8].split(":") if n else []
    if n and "GT" not in fmt:
        raise ValueError("Cohort records require FORMAT/GT")
    ix = fmt.index("GT") if n else 0
    ac = an = missing = carriers = 0
    bits = bytearray((n + 7) // 8)
    gt_digest = hashlib.sha256()
    for i, sample in enumerate(fields[9:]):
        parts = sample.split(":")
        gt = parts[ix] if ix < len(parts) else "."
        gt_digest.update((gt + "\n").encode())
        alleles = re.split(r"[/|]", gt)
        if any(a not in {"0", "1", "."} for a in alleles):
            raise ValueError("Non-biallelic GT after decomposition")
        ac += alleles.count("1")
        an += sum(a != "." for a in alleles)
        missing += int("." in alleles)
        if "1" in alleles:
            carriers += 1
            bits[i // 8] |= 1 << (i % 8)
    return {"ac": ac, "an": an, "af": ac / an if an else None,
            "samples": n, "missing_samples": missing, "carriers": carriers,
            "carrier_bits": int.from_bytes(bits, "little"), "gt_sha256": gt_digest.hexdigest()}


def merge_submitted(kind, status, policy):
    """Representation diagnostics must not veto Truvari's own matching rules."""
    return kind in {"INS", "DEL"} and (policy == "native" or status == "eligible")


def reuse_prepare(args):
    """Repartition the existing decomposed catalogue without repeating extraction."""
    source = args.reuse_prepared.resolve() / "prepare"
    directory = args.root / "prepare"
    if source == directory.resolve():
        raise ValueError("--reuse-prepared requires a different output root")
    parent = json.loads((source / "done.json").read_text())
    if parent["schema"] != SCHEMA:
        raise ValueError("Unsupported prepared catalogue schema")
    for name in ("chrom", "min_size"):
        if parent["key"][name] != getattr(args, name):
            raise ValueError(f"Cannot reuse catalogue with different {name}")
    if args.merge_policy == "strict" and parent["key"]["max_size"] != args.max_size:
        raise ValueError("Strict reuse requires the original max_size; use native policy for broader submission")
    names = ("original.vcf.gz", "original.vcf.gz.tbi", "source_alleles.tsv.gz", "samples.txt")
    snapshots = {name: common.identity(source / name) for name in names}
    def verify(name):
        if common.digest(source / name) != parent["artifacts"][name]:
            raise ValueError(f"Corrupt prepared input: {name}")
    with ThreadPoolExecutor(max_workers=args.threads) as pool:
        list(pool.map(verify, names))
    version = command(args, [args.bcftools, "--version"], "bcftools-version", True).splitlines()[0]
    key = {"preparation_version": 2, "reused_artifacts": {n: parent["artifacts"][n] for n in names},
           "chrom": args.chrom, "min_size": args.min_size, "max_size": args.max_size,
           "merge_policy": args.merge_policy, "bcftools": version}
    directory.mkdir(exist_ok=True)
    if valid_manifest(directory, key):
        logging.info("Reuse verified repartitioned catalogue")
        return
    (directory / "done.json").unlink(missing_ok=True)
    metadata = []
    with open_text(source / "source_alleles.tsv.gz") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            row["merge_submitted"] = int(merge_submitted(row["type"], row["status"], args.merge_policy))
            metadata.append(row)
    expected = {row["id"]: row for row in metadata}
    if len(expected) != len(metadata):
        raise ValueError("Duplicate source allele IDs")
    counts = Counter()
    with ExitStack() as stack:
        handles = {name: stack.enter_context(open_text(directory / f"{name}.staging.vcf.gz", "wt"))
                   for name in ("INS", "DEL", "passthrough")}
        lines = stack.enter_context(open_text(source / "original.vcf.gz"))
        for line in lines:
            if line.startswith("#"):
                for handle in handles.values():
                    handle.write(line)
                continue
            f = line.split("\t", 8)
            row = expected.pop(f[2], None)
            if row is None:
                raise ValueError("Unexpected/duplicate source VCF ID")
            if (f[0] != row["chrom"] or int(f[1]) != int(row["pos"])
                    or hashlib.sha256(f[3].encode()).hexdigest() != row["ref_sha256"]
                    or hashlib.sha256(f[4].encode()).hexdigest() != row["alt_sha256"]):
                raise ValueError("Source VCF and allele manifest disagree")
            route = row["type"] if row["merge_submitted"] else "passthrough"
            handles[route].write(line)
            counts[route] += 1
            if sum(counts.values()) % 10000 == 0:
                logging.info("REPARTITION records=%d partitions=%s", sum(counts.values()), dict(counts))
    if expected:
        raise ValueError("Source VCF is missing manifest alleles")
    artifacts = []
    for name in handles:
        target = directory / f"{name}.vcf.gz"
        staging = directory / f"{name}.staging.vcf.gz"
        pack(args, staging, target, "repartition-" + name)
        staging.unlink()
        artifacts.extend([target, Path(str(target) + ".tbi")])
    for name in ("original.vcf.gz", "original.vcf.gz.tbi", "samples.txt"):
        target = directory / name
        target.unlink(missing_ok=True)
        if name == "samples.txt":
            shutil.copyfile(source / name, target)
        else:
            try:
                os.link(source / name, target)
            except OSError:
                shutil.copyfile(source / name, target)
        if common.digest(target) != parent["artifacts"][name]:
            raise ValueError(f"Reused file checksum mismatch: {name}")
        artifacts.append(target)
    if any(common.identity(source / name) != snapshots[name] for name in names):
        raise ValueError("Source catalogue changed during repartition")
    tsv(directory / "source_alleles.tsv.gz", metadata, list(metadata[0]))
    common.json_write(directory / "summary.json", {
        "original_sv_alleles": len(metadata), "partition_counts": dict(counts),
        "status_counts": dict(Counter(row["status"] for row in metadata)),
        "reused_prepared_root": str(args.reuse_prepared.resolve()), "parameters": key,
        "sample_manifest_sha256": common.digest(directory / "samples.txt")})
    artifacts.extend(directory / name for name in ("source_alleles.tsv.gz", "summary.json"))
    manifest(directory, key, artifacts)


def update_counts(fields, info):
    gt = genotypes(fields)
    info.update(AC=str(gt["ac"]), AN=str(gt["an"]),
                AF=f'{gt["af"]:.9g}' if gt["af"] is not None else ".")
    return gt


def pack(args, source, target, label, threads=None, sort=False):
    threads = args.threads if threads is None else threads
    cmd = [args.bcftools, "sort", "-m", f"{max(1, args.memory_gb // 4)}G",
           "-T", args.root / f"sort-{label}"] if sort else [args.bcftools, "view", "--threads", threads]
    partial = Path(str(target) + ".partial")
    command(args, cmd + ["-Oz", "-o", partial, source], label)
    # Check complete BGZF stream, rather than trusting index creation alone.
    with gzip.open(partial, "rb") as handle:
        for _ in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            pass
    # Atomic replacement also protects hardlinked catalogues reused by other roots.
    partial.replace(target)
    index = Path(str(target) + ".tbi")
    index_partial = Path(str(index) + ".partial")
    command(args, [args.bcftools, "index", "-f", "-t", "--threads", threads,
                   "-o", index_partial, target], label + "-index")
    index_partial.replace(index)


def tsv(path, rows, columns):
    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def prepare(args):
    if args.reuse_prepared:
        return reuse_prepare(args)
    directory = args.root / "prepare"
    directory.mkdir(exist_ok=True)
    if args.input_bcf:
        source = args.input_bcf
    else:
        staged = require_manifest(args.root / "inputs")
        if staged["key"]["source"]["bcf"]["uri"] != args.source_uri:
            raise ValueError("Staged source differs from requested URI")
        source = args.root / "inputs" / args.source_uri.rsplit("/", 1)[1]
    version = command(args, [args.bcftools, "--version"], "bcftools-version", True).splitlines()[0]
    key = {"input_sha256": common.digest(source), "chrom": args.chrom, "min_size": args.min_size,
           "max_size": args.max_size, "bcftools": version, "multi_overlaps": "0", "filter": "all",
           "preparation_version": 2, "merge_policy": args.merge_policy}
    if valid_manifest(directory, key):
        logging.info("Reuse verified prepare artifacts")
        return
    (directory / "done.json").unlink(missing_ok=True)
    candidates = directory / "candidates.vcf.gz"
    expression = f'ILEN>={args.min_size} || ILEN<=-{args.min_size} || TYPE~"bnd" || ALT~"<"'
    count = 0
    samples = None
    with stream(args, [args.bcftools, "view", "--no-version", "--threads", args.threads,
                       "-r", args.chrom, "-i", expression, "-Ov", source], "extract") as lines, open_text(candidates, "wt") as dest:
        for line in lines:
            if line.startswith("##INFO=<ID=PM_"):
                raise ValueError("Input already has reserved PM_ tags; use original panel")
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                if len(samples) != len(set(samples)):
                    raise ValueError("Duplicate sample names")
                dest.write(HEADERS[0] + "\n")
            if not line.startswith("#"):
                count += 1
                fields = line.rstrip("\n").split("\t")
                info = info_parse(fields[7])
                info["PM_ROW"] = str(count)
                fields[7] = info_string(info)
                line = "\t".join(fields) + "\n"
            dest.write(line)
    if samples is None:
        raise ValueError("Missing VCF sample header")
    (directory / "samples.txt").write_text("".join(s + "\n" for s in samples), encoding="utf-8")
    metadata, source_rows, statuses = [], set(), Counter()
    with ExitStack() as stack:
        handles = {name: stack.enter_context(open_text(directory / f"{name}.staging.vcf.gz", "wt"))
                   for name in ("original", "INS", "DEL", "passthrough")}
        lines = stack.enter_context(stream(args, [args.bcftools, "norm", "--no-version", "-m", "-any",
                                                 "--multi-overlaps", "0", "-Ov", candidates], "split"))
        for line in lines:
            if line.startswith("#"):
                tag = re.match(r"##INFO=<ID=([^,]+)", line)
                if tag and tag[1] in REPLACED:
                    continue
                if line.startswith("#CHROM"):
                    for h in handles.values():
                        h.write("\n".join(HEADERS[1:]) + "\n")
                for h in handles.values():
                    h.write(line)
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] != args.chrom or "," in f[4] or len(f) != (9 + len(samples) if samples else 8):
                raise ValueError("Unexpected chromosome/ALT/sample count after splitting")
            info = info_parse(f[7])
            result = allele(f[3], f[4], int(f[1]), info, args.min_size, args.max_size)
            if result is None:
                continue
            kind, size, start, end, status = result
            source_row = info["PM_ROW"]
            source_rows.add(source_row)
            statuses[status] += 1
            original_id = f[2]
            f[2] = f"PM_{args.chrom}_{len(metadata) + 1:08d}"
            info = {k: v for k, v in info.items() if k not in REPLACED}
            info.update(SVTYPE=kind, PM_STATUS=status)
            if size is not None:
                info["SVLEN"] = str(-size if kind == "DEL" else size)
            update_counts(f, info)
            f[7] = info_string(info)
            line = "\t".join(f) + "\n"
            handles["original"].write(line)
            submitted = merge_submitted(kind, status, args.merge_policy)
            handles[kind if submitted else "passthrough"].write(line)
            metadata.append({"id": f[2], "source_row": source_row, "source_id": original_id,
                             "chrom": f[0], "pos": int(f[1]), "type": kind, "size": size,
                             "start0": start, "end0": end, "status": status, "merge_submitted": int(submitted),
                             "ref_sha256": hashlib.sha256(f[3].encode()).hexdigest(),
                             "alt_sha256": hashlib.sha256(f[4].encode()).hexdigest()})
    if not metadata:
        raise ValueError("No SV alleles found; check contig name and input")
    tsv(directory / "source_alleles.tsv.gz", metadata, list(metadata[0]))
    common.json_write(directory / "summary.json", {"candidate_source_rows": count,
                      "source_rows_with_sv": len(source_rows), "original_sv_alleles": len(metadata),
                      "samples": len(samples), "status_counts": dict(statuses),
                      "sample_manifest_sha256": common.digest(directory / "samples.txt"),
                      "input": str(source.resolve()), "parameters": key})
    artifacts = [directory / x for x in ("source_alleles.tsv.gz", "summary.json", "samples.txt")]
    # Compression is parallel internally; full genomic order is retained.
    for name in handles:
        target = directory / f"{name}.vcf.gz"
        pack(args, directory / f"{name}.staging.vcf.gz", target, "pack-" + name)
        artifacts += [target, Path(str(target) + ".tbi")]
        (directory / f"{name}.staging.vcf.gz").unlink()
    candidates.unlink()
    manifest(directory, key, artifacts)


def collapse_command(args, kind, output_dir):
    return [args.truvari, "collapse", "-i", args.root / "prepare" / f"{kind}.vcf.gz",
            "-o", output_dir / f"{kind}.kept.vcf", "-c", output_dir / f"{kind}.removed.vcf",
            "-p", str(args.threshold), "-P", str(args.threshold), "-r", str(args.refdist),
            "-O", "0", "-s", str(args.min_size), "-S", str(args.max_size),
            "-k", "first", "--gt", "off", "--median-info"]


def records(path):
    with open_text(path) as handle:
        for line in handle:
            if not line.startswith("#"):
                yield line.rstrip("\n").split("\t")


def vcf_headers(path):
    with open_text(path) as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            yield line


def merge(args):
    prepared = require_manifest(args.root / "prepare")
    for name, value in (("chrom", args.chrom), ("min_size", args.min_size), ("max_size", args.max_size)):
        if prepared["key"][name] != value:
            raise ValueError(f"Prepare parameter {name} differs; rerun prepare")
    if prepared["key"].get("merge_policy", "strict") != args.merge_policy:
        raise ValueError("Prepare merge policy differs; use a new root and --reuse-prepared")
    directory = args.root / "merge"
    directory.mkdir(exist_ok=True)
    version = command(args, [args.truvari, "version"], "truvari-version", True).strip()
    if not re.search(r"(?<![\d.])5\.5\.0(?![\d.])", version):
        raise ValueError("Expected Truvari 5.5.0; inspect logs/truvari-version.log")
    key = {"prepare": prepared, "threshold": args.threshold, "refdist": args.refdist,
           "truvari": version, "keep": "first", "gt": "off", "chain": False,
           "consolidate": True, "typeignore": False, "pctovl": 0, "median_info": True,
           "separate_only": args.separate_only}
    if valid_manifest(directory, key):
        logging.info("Reuse verified merge artifacts")
        return
    (directory / "done.json").unlink(missing_ok=True)
    workers = min(2, args.threads)

    def run_type(kind):
        source = args.root / "prepare" / f"{kind}.vcf.gz"
        if next(records(source), None) is None:
            for role in ("kept", "removed"):
                with open(directory / f"{kind}.{role}.vcf", "w", encoding="utf-8") as out:
                    out.writelines(vcf_headers(source))
        else:
            command(args, collapse_command(args, kind, directory), f"collapse-{kind}")
        for role in ("kept", "removed"):
            path = directory / f"{kind}.{role}.vcf"
            pack(args, path, Path(str(path) + ".gz"), f"pack-{kind}-{role}", args.threads // workers, sort=True)
            path.unlink()

    with ThreadPoolExecutor(max_workers=workers) as pool:
        list(pool.map(run_type, ("INS", "DEL")))
    if args.separate_only:
        manifest(directory, key, [directory / f"{kind}.{role}.vcf.gz{ext}"
                 for kind in ("INS", "DEL") for role in ("kept", "removed") for ext in ("", ".tbi")])
        return
    parts = [directory / f"{kind}.kept.vcf.gz" for kind in ("INS", "DEL")]
    parts.append(args.root / "prepare" / "passthrough.vcf.gz")
    headers, sample_header = {}, None
    for part in parts:
        for line in vcf_headers(part):
            if line.startswith("#CHROM"):
                if sample_header is not None and sample_header != line:
                    raise ValueError("Sample order changed in merge")
                sample_header = line
            else:
                match = re.match(r"##(INFO|FORMAT|FILTER|contig)=<ID=([^,>]+)", line)
                ident = tuple(match.groups()) if match else line
                if ident in headers and headers[ident] != line:
                    raise ValueError("Conflicting VCF header definitions")
                headers[ident] = line
    with open(directory / "new.unsorted.vcf", "w", encoding="utf-8") as out:
        out.writelines(headers.values())
        out.write(sample_header)
        for part in parts:
            for f in records(part):
                info = info_parse(f[7])
                update_counts(f, info)
                f[7] = info_string(info)
                out.write("\t".join(f) + "\n")
    target = directory / "new.vcf.gz"
    pack(args, directory / "new.unsorted.vcf", target, "pack-new", sort=True)
    (directory / "new.unsorted.vcf").unlink()
    artifacts = [target, Path(str(target) + ".tbi")]
    artifacts += [directory / f"{kind}.{role}.vcf.gz{ext}" for kind in ("INS", "DEL")
                  for role in ("kept", "removed") for ext in ("", ".tbi")]
    manifest(directory, key, artifacts)


def load_catalogue(path):
    result, previous = {}, -1
    for f in records(path):
        info = info_parse(f[7])
        classified = allele(f[3], f[4], int(f[1]), info, 1, -1)
        if classified is None or f[2] in result:
            raise ValueError("Invalid or duplicate pilot ID in catalogue")
        kind, size, start, end, _ = classified
        if int(f[1]) < previous:
            raise ValueError("Unsorted catalogue")
        previous = int(f[1])
        result[f[2]] = {"id": f[2], "pos": int(f[1]), "type": kind, "size": size,
                        "start0": start, "end0": end, "status": info["PM_STATUS"],
                        "collapse_id": info.get("CollapseId"),
                        "num_collapsed": int(info.get("NumCollapsed", 0)), **genotypes(f)}
    return result


def quantiles(values):
    values = sorted(v for v in values if v is not None)
    if not values:
        return {k: None for k in ("min", "p25", "median", "p75", "p95", "max", "mean")}
    def q(p):
        ix = (len(values) - 1) * p
        low, high = math.floor(ix), math.ceil(ix)
        return values[low] + (values[high] - values[low]) * (ix - low)
    return dict(zip(("min", "p25", "median", "p75", "p95", "max", "mean"),
                    (values[0], q(.25), q(.5), q(.75), q(.95), values[-1], sum(values) / len(values))))


def distances(rows):
    """Adjacent POS difference and nearest POS neighbor; interval gap to union
    of preceding changed-reference intervals (handles nested deletions).
    """
    rows = sorted(rows, key=lambda x: (x["pos"], x["id"]))
    positions = [r["pos"] for r in rows]
    gaps = [b - a for a, b in zip(positions, positions[1:])]
    nearest = [min(([gaps[i - 1]] if i else []) + ([gaps[i]] if i < len(gaps) else []))
               for i in range(len(rows))] if len(rows) > 1 else []
    ordered = sorted(rows, key=lambda x: (x["start0"], x["end0"], x["id"]))
    interval_gaps, overlaps, right = [], 0, None
    for r in ordered:
        if right is not None:
            interval_gaps.append(max(0, r["start0"] - right))
            overlaps += int(r["start0"] < right)
        right = max(right if right is not None else r["end0"], r["end0"])
    return gaps, nearest, interval_gaps, overlaps


def audit_clusters(original, new, removed_paths):
    kept_ids = set(new)
    if not kept_ids <= set(original):
        raise ValueError("New catalogue contains unknown source IDs")
    lookup = {(r["type"], r["collapse_id"]): r["id"] for r in new.values() if r["collapse_id"]}
    members = {ident: [ident] for ident in new}
    links = [{"source_id": ident, "representative_id": ident, "role": "kept",
              "sequence_similarity": None, "size_similarity": None} for ident in new]
    seen = set(kept_ids)
    for kind, path in removed_paths.items():
        for f in records(path):
            info = info_parse(f[7])
            ident = f[2]
            if ident in seen or ident not in original:
                raise ValueError("Removed ID duplicated, retained, or unknown")
            representative = lookup.get((kind, info.get("MatchId")))
            if representative is None or original[ident]["type"] != new[representative]["type"]:
                raise ValueError("Unresolved or cross-type cluster membership")
            members[representative].append(ident)
            seen.add(ident)
            links.append({"source_id": ident, "representative_id": representative, "role": "removed",
                          "sequence_similarity": info.get("PctSeqSimilarity"),
                          "size_similarity": info.get("PctSizeSimilarity")})
    if seen != set(original):
        raise ValueError("Conservation failed: original != new + removed")
    clusters = []
    for ident, ids in members.items():
        group, representative = [original[i] for i in ids], new[ident]
        if len(ids) - 1 != representative["num_collapsed"]:
            raise ValueError("NumCollapsed disagrees with membership")
        union = repeated = 0
        for row in group:
            repeated |= union & row["carrier_bits"]
            union |= row["carrier_bits"]
        sizes = [r["size"] for r in group if r["size"] is not None]
        clusters.append({"representative_id": ident, "type": representative["type"],
                         "members": len(ids), "position_span_bp": max(r["pos"] for r in group) - min(r["pos"] for r in group),
                         "min_size": min(sizes) if sizes else None, "max_size": max(sizes) if sizes else None,
                         "samples_carrying_multiple_members": repeated.bit_count(),
                         "union_carriers": union.bit_count(), "new_carriers": representative["carriers"],
                         "lost_carriers": (union & ~representative["carrier_bits"]).bit_count(),
                         "gained_carriers": (representative["carrier_bits"] & ~union).bit_count(),
                         "representative_ac_before": original[ident]["ac"], "ac_after": representative["ac"],
                         "representative_gt_changed": int(original[ident]["gt_sha256"] != representative["gt_sha256"])})
    return links, clusters


def make_plots(directory, catalogues, counts):
    os.environ["MPLBACKEND"] = "Agg"
    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    matplotlib.rcParams.update({"font.size": 18, "axes.titlesize": 20, "axes.labelsize": 18,
                               "xtick.labelsize": 16, "ytick.labelsize": 16,
                               "pdf.fonttype": 42, "ps.fonttype": 42})
    colors = {"original": "#0072B2", "new": "#D55E00"}
    paths = []
    fig, ax = plt.subplots(figsize=(11, 8.5), layout="constrained")
    labels = [c["type"] for c in counts]
    for j, name in enumerate(catalogues):
        ax.bar([i + (j - .5) * .36 for i in range(len(labels))],
               [c[name] for c in counts], width=.36, label=name.capitalize(), color=colors[name])
    ax.set(xticks=range(len(labels)), xticklabels=labels, ylabel="SV allele records", title="chr1 SV counts")
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.legend(loc="lower center", bbox_to_anchor=(.5, 1.01), ncols=2, frameon=False)
    for ext in ("png", "pdf"):
        path = directory / f"counts.{ext}"
        fig.savefig(path, dpi=180)
        paths.append(path)
    plt.close(fig)
    for kind in ("ALL", "INS", "DEL"):
        fig, ax = plt.subplots(figsize=(11, 8.5), layout="constrained")
        xmax = 10
        for name, rows in catalogues.items():
            selected = [r for r in rows.values() if kind == "ALL" or r["type"] == kind]
            values = sorted(distances(selected)[0])
            if values:
                tail = max(10, values[-1] * 1.1)
                xmax = max(xmax, tail)
                ax.step([0] + values + [tail], [0] + [(i + 1) / len(values) for i in range(len(values))] + [1],
                        where="post", label=name.capitalize(), color=colors[name], linewidth=2.5)
        ax.set_xscale("symlog", linthresh=1)
        ax.set(xlabel="Adjacent VCF POS distance (bp; zero retained)", ylabel="Cumulative fraction",
               title=f"chr1 {kind}: inter-SV spacing", ylim=(0, 1.02), xlim=(0, xmax))
        if ax.lines:
            ax.legend(loc="lower center", bbox_to_anchor=(.5, 1.01), ncols=2, frameon=False)
        else:
            ax.text(.5, .5, "Fewer than two records", ha="center", transform=ax.transAxes)
        for ext in ("png", "pdf"):
            path = directory / f"spacing_{kind}.{ext}"
            fig.savefig(path, dpi=180)
            paths.append(path)
        plt.close(fig)
    return paths


def compare(args):
    prepared = require_manifest(args.root / "prepare")
    merged = require_manifest(args.root / "merge")
    if merged["key"].get("separate_only"):
        raise ValueError("Use compare_chr1_sv_alleles.py for --separate-only merge outputs")
    if merged["key"]["prepare"] != prepared:
        raise ValueError("Merge is stale relative to prepare; rerun merge")
    directory = args.root / "compare"
    directory.mkdir(exist_ok=True)
    key = {"prepare": prepared, "merge": merged, "density_bin_bp": args.density_bin_bp}
    if valid_manifest(directory, key):
        logging.info("Reuse verified comparison artifacts")
        return
    (directory / "done.json").unlink(missing_ok=True)
    with ThreadPoolExecutor(max_workers=min(2, args.threads)) as pool:
        catalogues = dict(zip(("original", "new"), pool.map(load_catalogue,
            (args.root / "prepare" / "original.vcf.gz", args.root / "merge" / "new.vcf.gz"))))
    original, new = catalogues["original"], catalogues["new"]
    links, clusters = audit_clusters(original, new, {k: args.root / "merge" / f"{k}.removed.vcf.gz" for k in ("INS", "DEL")})
    kinds = ["ALL", "INS", "DEL", "OTHER"]
    def selected(rows, kind):
        return [r for r in rows.values() if kind == "ALL" or r["type"] == kind or (kind == "OTHER" and r["type"] not in {"INS", "DEL"})]
    counts, statistics, density, size_counts = [], [], [], []
    for kind in kinds:
        before, after = len(selected(original, kind)), len(selected(new, kind))
        counts.append({"type": kind, "original": before, "new": after, "removed": before - after,
                       "percent_reduction": 100 * (before - after) / before if before else None,
                       "eligible_original": sum(merge_submitted(r["type"], r["status"],
                           prepared["key"].get("merge_policy", "strict")) for r in selected(original, kind))})
    for name, catalogue in catalogues.items():
        exported = []
        for r in catalogue.values():
            exported.append({**r, "missing_fraction": r["missing_samples"] / r["samples"] if r["samples"] else None})
        cols = [c for c in exported[0] if c not in {"carrier_bits", "gt_sha256"}] if exported else ["id"]
        tsv(directory / f"{name}.records.tsv.gz", exported, cols)
        for kind in kinds:
            rows = selected(catalogue, kind)
            gaps, nearest, interval, overlaps = distances(rows)
            for metric, values in (("adjacent_pos_bp", gaps), ("nearest_pos_bp", nearest),
                                   ("preceding_interval_union_gap_bp", interval), ("sv_size_bp", [r["size"] for r in rows]),
                                   ("allele_frequency", [r["af"] for r in rows]),
                                   ("missing_sample_fraction", [r["missing_samples"] / r["samples"] for r in rows if r["samples"]])):
                usable = [v for v in values if v is not None]
                statistics.append({"callset": name, "type": kind, "metric": metric, "n": len(usable),
                                   "n_zero": sum(v == 0 for v in usable),
                                   "n_le_50": sum(v <= 50 for v in usable) if metric.endswith("_bp") else None,
                                   "n_le_500": sum(v <= 500 for v in usable) if metric.endswith("_bp") else None,
                                   "n_le_1000": sum(v <= 1000 for v in usable) if metric.endswith("_bp") else None,
                                   "overlapping_preceding_intervals": overlaps if metric == "preceding_interval_union_gap_bp" else None,
                                   **quantiles(usable)})
            bins = Counter((r["pos"] - 1) // args.density_bin_bp for r in rows)
            # Sparse occupied bins; absence is zero, not missing data.
            density += [{"callset": name, "type": kind, "start0": b * args.density_bin_bp,
                         "end0": (b + 1) * args.density_bin_bp, "count": n} for b, n in sorted(bins.items())]
            for label, low, high in (("50-99", 50, 100), ("100-499", 100, 500), ("500-999", 500, 1000),
                                     ("1k-9.9k", 1000, 10000), ("10k-99.9k", 10000, 100000),
                                     (">=100k", 100000, math.inf), ("unknown", None, None)):
                n = sum((r["size"] is None if low is None else r["size"] is not None and low <= r["size"] < high) for r in rows)
                size_counts.append({"callset": name, "type": kind, "size_bin": label, "count": n})
    tables = {"counts.tsv": counts, "statistics.tsv": statistics, "size_counts.tsv": size_counts,
              "density.tsv": density, "membership.tsv.gz": links, "clusters.tsv.gz": clusters}
    for filename, rows in tables.items():
        tsv(directory / filename, rows, list(rows[0]) if rows else ["empty"])
    summary = {"schema": SCHEMA, "counts": counts, "conservation_passed": True,
               "merged_clusters": sum(c["members"] > 1 for c in clusters),
               "clusters_with_same_sample_cocarriage": sum(c["samples_carrying_multiple_members"] > 0 for c in clusters),
               "representatives_with_gt_changes": sum(c["representative_gt_changed"] for c in clusters),
               "clusters_with_lost_carriers": sum(c["lost_carriers"] > 0 for c in clusters),
               "clusters_with_gained_carriers": sum(c["gained_carriers"] > 0 for c in clusters),
               "cluster_size": quantiles([c["members"] for c in clusters]),
               "notes": ["Counts are decomposed SV alleles, not all panel rows.",
                         "Adjacent/nearest distances use VCF POS, within each displayed type.",
                         "Zero distance retains distinct same-position alleles.",
                         "Consolidated GT is exploratory; it is not validated grouped allele dosage.",
                         "AF is recomputed; other inherited annotations may be representative-specific."]}
    common.json_write(directory / "summary.json", summary)
    plots = make_plots(directory, catalogues, counts)
    artifacts = [directory / name for name in tables]
    artifacts += [directory / f"{name}.records.tsv.gz" for name in catalogues]
    artifacts += [directory / "summary.json"] + plots
    manifest(directory, key, artifacts)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--phase", choices=("plan", "stage", "prepare", "merge", "compare", "all"), default="plan")
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--source-uri", default=SOURCE)
    p.add_argument("--confirm-source-uri", help="Exact source URI approved for staging")
    p.add_argument("--input-bcf", type=Path, help="Already-local indexed BCF/VCF; bypass cloud staging")
    p.add_argument("--reuse-prepared", type=Path, help="Existing run root; repartition its prepared catalogue without staging/extraction")
    p.add_argument("--merge-policy", choices=("native", "strict"), default="native",
                   help="native submits all INS/DEL representations; strict reproduces the original representation filters")
    p.add_argument("--separate-only", action="store_true", help="Merge outputs per type only; run compare_chr1_sv_alleles.py afterward")
    p.add_argument("--chrom", default="chr1")
    p.add_argument("--threshold", type=float, default=.7)
    p.add_argument("--refdist", type=int, default=500)
    p.add_argument("--min-size", type=int, default=50)
    p.add_argument("--max-size", type=int, default=-1, help="Truvari size ceiling; -1 disables the ceiling")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--memory-gb", type=int, default=26, help="Resource budget; bcftools sort uses one quarter per process; not an OS cap")
    p.add_argument("--scratch-reserve-gb", type=int, default=64, help="Free disk floor beyond source size, not a full-run prediction")
    p.add_argument("--density-bin-bp", type=int, default=1000000)
    p.add_argument("--billing-project", default=None,
                   help="Optional explicit override; otherwise use existing gcloud configuration")
    p.add_argument("--bcftools", default=os.environ.get("BCFTOOLS_BIN", "bcftools"))
    p.add_argument("--truvari", default="truvari")
    p.add_argument("--gcloud", default="gcloud")
    args = p.parse_args(argv)
    if not (0 < args.threshold <= 1 and args.threads >= 1 and args.memory_gb >= 4 and args.refdist >= 0
            and args.min_size >= 50 and (args.max_size == -1 or args.max_size >= args.min_size)
            and args.density_bin_bp >= 1 and args.scratch_reserve_gb >= 0):
        p.error("Invalid thresholds/resource settings; minimum SV size must be at least 50")
    if not args.source_uri.startswith("gs://") or not args.source_uri.endswith(".bcf"):
        p.error("Source must be an explicit gs:// BCF object")
    if args.input_bcf and args.phase == "stage":
        p.error("--input-bcf cannot be combined with cloud stage phase")
    if args.reuse_prepared and (args.input_bcf or args.phase == "stage"):
        p.error("--reuse-prepared cannot be combined with --input-bcf or cloud staging")
    if args.separate_only and args.phase in {"all", "compare"}:
        p.error("--separate-only requires separate prepare/merge phases and compare_chr1_sv_alleles.py")
    args.root = args.root.resolve()
    return args


def main(argv=None):
    args = parse_args(argv)
    if args.phase == "plan":
        print(json.dumps({"parameters": vars(args), "source_index": args.source_uri + ".csi",
                          "phases": ["stage", "prepare", "merge", "compare"],
                          "collapse_commands": [collapse_command(args, k, args.root / "merge") for k in ("INS", "DEL")],
                          "uploads": False, "seed": SEED}, default=str, indent=2))
        return
    (args.root / "logs").mkdir(parents=True, exist_ok=True)
    logging.basicConfig(filename=args.root / "logs" / "workflow.log", level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s", force=True)
    lock = args.root / ".running.lock"
    try:
        descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        raise SystemExit("Run lock exists. Check the VM process; remove a stale .running.lock only after it has stopped.")
    os.write(descriptor, f"{os.getpid()}\n".encode())
    os.close(descriptor)
    try:
        logging.info("PARAMETERS %s", json.dumps(vars(args), default=str, sort_keys=True))
        logging.info("PROVENANCE %s", json.dumps({"python": sys.version, "script_sha256": common.digest(__file__),
                     "seed": SEED, "schema": SCHEMA, "platform": sys.platform}))
        common.json_write(args.root / "logs" / f"run-{time.time_ns()}.json", {
            "parameters": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
            "python": sys.version, "script_sha256": common.digest(__file__), "seed": SEED,
            "schema": SCHEMA, "packages": dict(sorted((d.metadata["Name"], d.version)
                                                     for d in importlib.metadata.distributions()))})
        phases = ["stage", "prepare", "merge", "compare"] if args.phase == "all" else [args.phase]
        if (args.input_bcf or args.reuse_prepared) and "stage" in phases:
            phases.remove("stage")
        for phase in phases:
            start = time.monotonic()
            print(f"Starting {phase}; details stay in VM logs.", flush=True)
            globals()[phase](args)
            logging.info("PHASE_COMPLETE %s elapsed_seconds=%.3f", phase, time.monotonic() - start)
            print(f"Completed {phase}.", flush=True)
    except Exception:
        logging.exception("Workflow failed")
        raise SystemExit(f"Workflow failed. Inspect {args.root / 'logs' / 'workflow.log'} on the VM.") from None
    finally:
        lock.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
