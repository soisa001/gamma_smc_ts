#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = ["google-crc32c==1.8.0"]
# ///
"""Stage one AoU chromosome, transfer singleton phase, retain results locally.

Defaults to chr1. No upload operation exists. Run in the authorized data VM.
"""
from __future__ import annotations

import argparse
import base64
import json
import logging
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

import phase_singletons as ph

PREFIX = "gs://rw-long-reads-transfer-2026-06-17/v9/lrWGS"
STAGING_VERSION = 1


def source_paths(chrom):
    return {
        "hiphase": f"{PREFIX}/phasing/hiphase/aou_lr_phase2_v1.hiphase.short_sv.{chrom}.bcf",
        "scaffold": f"{PREFIX}/panel/panel/panel_bubble_split_vcf/aou_lr_phase2_v1.{chrom}.bubble.split.bcf",
    }


def sync_command(gcloud, uri, dest, billing_project):
    parent, name = uri.rsplit("/", 1)
    # rsync takes a DIRECTORY source. Exclude everything except this exact BCF
    # and its optional CSI, so chr1 never downloads chr10 or another chromosome.
    exclude = rf"^(?!{re.escape(name)}(?:\.csi)?$).*"
    return [gcloud, "storage", "rsync", parent + "/", str(dest),
            "--exclude", exclude, "--checksums-only", "--do-not-decompress",
            "--billing-project", billing_project]


def command(cmd, capture=False, env=None):
    logging.info("COMMAND %s", json.dumps([str(x) for x in cmd]))
    result = subprocess.run([str(x) for x in cmd], text=True, capture_output=capture, env={**os.environ, **env} if env else None)
    if result.returncode:
        raise subprocess.CalledProcessError(result.returncode, cmd, result.stdout, result.stderr)
    return result.stdout if capture else None


def describe(gcloud, uri, billing_project, optional=False):
    cmd = [gcloud, "storage", "objects", "describe", uri, "--raw", "--format=json", "--billing-project", billing_project]
    try:
        data = json.loads(command(cmd, capture=True))
    except subprocess.CalledProcessError as error:
        if optional and re.search(r"(?:\b404\b|\bNOT_FOUND\b)", error.stderr or ""):
            return None
        raise
    if not all(data.get(k) is not None for k in ("size", "generation", "crc32c")):
        raise ValueError(f"Object metadata lacks size/generation/CRC32C: {uri}")
    if int(data["size"]) <= 0:
        raise ValueError(f"Empty cloud object: {uri}")
    return {"uri": uri, "size": int(data["size"]), "generation": str(data["generation"]), "crc32c": data["crc32c"]}


def crc32c(path):
    import google_crc32c
    checksum = google_crc32c.Checksum()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            checksum.update(block)
    return base64.b64encode(checksum.digest()).decode("ascii")


def verify_object(path, metadata):
    if not path.is_file() or path.stat().st_size != metadata["size"]:
        raise ValueError(f"Staged file size mismatch: {path}")
    actual = crc32c(path)
    if actual != metadata["crc32c"]:
        raise ValueError(f"Staged file CRC32C mismatch: {path}")
    return {"identity": ph.identity(path), "crc32c": actual}


def metadata_for(args, uri):
    return {"bcf": describe(args.gcloud, uri, args.billing_project),
            "csi": describe(args.gcloud, uri + ".csi", args.billing_project, optional=True)}


def cached_stage(marker, metadata, bcf, csi):
    if not marker.exists():
        return False
    try:
        old = json.loads(marker.read_text())
        return (old["staging_version"] == STAGING_VERSION and old["remote"] == metadata
                and all(p.is_file() and old["files"][name]["identity"] == ph.identity(p)
                        for name, p in (("bcf", bcf), ("csi", csi))))
    except (ValueError, KeyError, OSError):
        return False


def stage_one(args, name, metadata):
    uri = source_paths(args.chrom)[name]
    dest = args.root / "inputs" / name
    dest.mkdir(parents=True, exist_ok=True)
    bcf = dest / uri.rsplit("/", 1)[1]
    csi = Path(str(bcf) + ".csi")
    marker = dest / "stage.done.json"
    if cached_stage(marker, metadata, bcf, csi):
        logging.info("Validated staging manifest; reuse %s", bcf)
        return bcf
    marker.unlink(missing_ok=True)
    print(f"Staging {args.chrom} {name} ({metadata['bcf']['size'] / 2**30:.1f} GiB)", flush=True)
    command(sync_command(args.gcloud, uri, dest, args.billing_project), env={
        "CLOUDSDK_STORAGE_PROCESS_COUNT": "1", "CLOUDSDK_STORAGE_THREAD_COUNT": str(args.io_threads)})
    files = {"bcf": verify_object(bcf, metadata["bcf"])}
    if metadata["csi"] is None:
        command([args.bcftools, "index", "--force", "--csi", "--threads", str(args.io_threads), str(bcf)])
        files["csi"] = {"identity": ph.identity(csi), "sha256": ph.digest(csi)}
    else:
        files["csi"] = verify_object(csi, metadata["csi"])
    # Fail if the source changed while rsync was running rather than mixing versions.
    if metadata_for(args, uri) != metadata:
        raise ValueError(f"Cloud BCF/index changed during staging: {uri}; rerun to stage a consistent pair")
    count = command([args.bcftools, "index", "--nrecords", str(bcf)], capture=True).strip()
    if not count.isdigit():
        raise ValueError(f"Could not read BCF index count: {bcf}")
    ph.json_write(marker, {"staging_version": STAGING_VERSION, "remote": metadata, "files": files, "indexed_records": int(count)})
    return bcf


def stage(args):
    remotes = {name: metadata_for(args, uri) for name, uri in source_paths(args.chrom).items()}
    required = 0
    for name, meta in remotes.items():
        dest = args.root / "inputs" / name
        bcf = dest / source_paths(args.chrom)[name].rsplit("/", 1)[1]
        # A full object may be needed for an atomic replacement. Do not count
        # an unvalidated partial download as reusable free capacity.
        if not cached_stage(dest / "stage.done.json", meta, bcf, Path(str(bcf) + ".csi")):
            required += sum(metadata["size"] for metadata in meta.values() if metadata)
    free = shutil.disk_usage(args.root).free
    floor = required + args.scratch_reserve_gb * 2**30
    if free < floor:
        raise RuntimeError(f"Need at least {floor / 2**30:.1f} GiB free for staging + scratch reserve; have {free / 2**30:.1f} GiB. Reserve is a floor, not full-run sizing.")
    paths = {name: stage_one(args, name, meta) for name, meta in remotes.items()}
    ph.json_write(args.root / "staging.json", {"chrom": args.chrom, "remote": remotes, "local": {k: str(v) for k, v in paths.items()}})
    return paths


def staged_paths(args):
    paths = {}
    for name, uri in source_paths(args.chrom).items():
        dest = args.root / "inputs" / name
        bcf = dest / uri.rsplit("/", 1)[1]
        marker = dest / "stage.done.json"
        if not marker.exists():
            raise ValueError(f"Stage phase has not completed: {marker}")
        metadata = json.loads(marker.read_text())["remote"]
        if not cached_stage(marker, metadata, bcf, Path(str(bcf) + ".csi")):
            raise ValueError(f"Staged files changed; rerun --phase stage: {bcf}")
        paths[name] = bcf
    return paths


def plan(args):
    return {
        "chromosome": args.chrom, "phase": args.phase, "genome_build": args.genome_build,
        "root": str(args.root), "source_bcf": source_paths(args.chrom),
        "optional_source_indexes": {k: v + ".csi" for k, v in source_paths(args.chrom).items()},
        "staging_commands": [sync_command(args.gcloud, uri, args.root / "inputs" / name, args.billing_project or "$GOOGLE_PROJECT") for name, uri in source_paths(args.chrom).items()],
        "output": str(args.root / "results" / args.chrom / "augmented.vcf.gz"),
        "workers": args.workers, "active_chromosome_workers": 1,
        "chunks_per_chrom": args.chunks_per_chrom, "overlap_bp": args.overlap_bp,
        "active_chunk_workers": min(args.workers, args.chunks_per_chrom),
        "io_threads": args.io_threads, "memory_gb": args.memory_gb,
        "gcloud_processes": 1, "gcloud_threads_per_process": args.io_threads,
        "min_anchors": args.min_anchors, "seed": ph.SEED, "uploads": False,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chrom", choices=[f"chr{i}" for i in range(1, 23)], default="chr1")
    parser.add_argument("--root", type=Path, default=Path.home() / "singleton_transfer_chr1")
    parser.add_argument("--phase", choices=("all", "stage", "prepare", "run"), default="all")
    parser.add_argument("--genome-build", default="GRCh38")
    parser.add_argument("--billing-project", default=os.environ.get("GOOGLE_PROJECT", ""))
    parser.add_argument("--gcloud", default="gcloud")
    repo = Path(__file__).resolve().parents[1]
    local_bcftools = next((p for p in (repo / ".native/bin/bcftools", repo / ".tools/singleton-native/bin/bcftools") if p.is_file()), None)
    parser.add_argument("--bcftools", default=os.environ.get("BCFTOOLS_BIN") or (str(local_bcftools) if local_bcftools else "bcftools"))
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--chunks-per-chrom", type=int, default=4)
    parser.add_argument("--overlap-bp", type=int, default=1000000)
    parser.add_argument("--io-threads", type=int, default=4, help="gcloud download threads and index-creation threads; does not parallelize phasing")
    parser.add_argument("--memory-gb", type=float, default=26)
    parser.add_argument("--min-anchors", type=int, choices=(1, 2), default=2)
    parser.add_argument("--scratch-reserve-gb", type=float, default=64, help="Minimum free space in addition to objects needing staging; not full-run sizing")
    parser.add_argument("--dry-run", action="store_true", help="Print resolved commands/paths with no tool/cloud calls or filesystem writes")
    args = parser.parse_args(argv)
    if not 1 <= args.workers <= 32 or not 1 <= args.io_threads <= 32 or args.memory_gb <= 0 or args.scratch_reserve_gb < 0 or not 1 <= args.chunks_per_chrom <= 32 or args.overlap_bp < 0:
        parser.error("Invalid resource setting")
    args.root = args.root.resolve()
    if args.dry_run:
        print(json.dumps(plan(args), indent=2))
        return
    for executable in ([args.gcloud] if args.phase in ("all", "stage") else []) + [args.bcftools]:
        if shutil.which(executable) is None:
            parser.error(f"Required executable unavailable: {executable}. Run in the authorized VM with gcloud and bcftools installed.")
    if args.phase in ("all", "stage") and not args.billing_project:
        parser.error("Set GOOGLE_PROJECT or --billing-project for requester-pays staging")
    args.root.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(filename=args.root / "runner.log", level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", force=True)
    lock = args.root / ".runner.lock"
    try:
        fd = os.open(lock, os.O_WRONLY | os.O_CREAT | os.O_EXCL)
    except FileExistsError:
        raise RuntimeError(f"Runner is already active or left a stale lock: {lock}") from None
    try:
        with os.fdopen(fd, "w") as handle:
            handle.write(str(os.getpid()))
        logging.info("RUN %s", json.dumps(plan(args)))
        logging.info("BCFTOOLS VERSION %s", command([args.bcftools, "--version"], capture=True))
        if args.phase in ("all", "stage"):
            logging.info("GCLOUD VERSION %s", command([args.gcloud, "--version"], capture=True))
        paths = stage(args) if args.phase in ("all", "stage") else staged_paths(args)
        if args.phase == "stage":
            print(f"Staging complete: {args.root}")
            return
        manifest = args.root / "manifest.json"
        scripts = Path(__file__).resolve().parent
        # Recreate only the small manifest; same sample/input identity preserves
        # phasing caches, regardless of current repository commit.
        command([sys.executable, scripts / "make_singleton_manifest.py", "--scaffold", paths["scaffold"], "--hiphase", paths["hiphase"],
                 "--contigs", args.chrom, "--genome-build", args.genome_build, "--assert-chromosome-phase", "--bcftools", args.bcftools, "--output", manifest])
        if args.phase == "prepare":
            print(f"Manifest ready: {manifest}")
            return
        command([sys.executable, scripts / "phase_singletons.py", "--manifest", manifest, "--out", args.root / "results",
                 "--workers", str(args.workers), "--chunks-per-chrom", str(args.chunks_per_chrom), "--overlap-bp", str(args.overlap_bp),
                 "--memory-gb", str(args.memory_gb), "--min-anchors", str(args.min_anchors), "--bcftools", args.bcftools])
        print(f"Local output: {args.root / 'results' / args.chrom / 'augmented.vcf.gz'}")
    except BaseException:
        logging.exception("Runner failed; completed stages can be reused on restart")
        raise
    finally:
        lock.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
