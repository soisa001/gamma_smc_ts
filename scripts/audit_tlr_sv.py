#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = ["truvari==5.5.0", "google-crc32c==1.8.0"]
# ///
"""Local TLR-locus pair scores and native cluster audit; no individual output on stdout."""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from contextlib import ExitStack, redirect_stdout
import copy
import gzip
import importlib.metadata
from itertools import combinations
import json
import logging
import os
from pathlib import Path
import re
import sys

import edlib
import pysam
import truvari
from truvari.collapse import CollapsedCalls

import run_chr1_sv_merge as pilot
import run_popped_simple_sv as simple

SOURCE = simple.representation.SOURCE.replace(".chr1.", ".chr4.")
VERSION = 1


def prepare(args):
    directory = args.root / "prepare"
    directory.mkdir(exist_ok=True)
    source = args.input_bcf
    if source is None:
        pilot.stage(args)
        source = args.root / "inputs" / args.source_uri.rsplit("/", 1)[1]
    source = source.resolve()
    indexes = [p for suffix in (".csi", ".tbi") if (p := Path(str(source) + suffix)).exists()]
    if not indexes:
        raise ValueError("An indexed original popped BCF/VCF is required")
    origin = [pilot.common.identity(p) for p in [source, *indexes]]
    bounds = max(1, args.start - args.flank), args.end + args.flank
    key = {"version": VERSION, "origin": origin, "chrom": args.chrom, "start": args.start, "end": args.end,
           "flank": args.flank, "merge": "simple_trim_only", "extract": simple.EXTRACT,
           "label_policy": "unique audit ID in working copies only; original ID retained in allele map"}
    cached = pilot.valid_manifest(directory)
    if cached and cached["key"]["settings"] == key:
        return cached
    (directory / "done.json").unlink(missing_ok=True)
    samples = pilot.command(args, [args.bcftools, "query", "-l", source], "samples", True)
    if not samples.strip():
        raise ValueError("Cohort genotypes are required to evaluate --gt all")
    (directory / "samples.txt").write_text(samples)
    counts = {k: Counter() for k in ("before", "simple", "passthrough")}
    written = Counter()
    parts = ("INS", "DEL", "passthrough.sv50")
    with ExitStack() as stack:
        handles = {k: stack.enter_context(pysam.BGZFile(str(directory / f"{k}.partial.vcf.gz"), "w")) for k in parts}
        meta = stack.enter_context(gzip.open(directory / "alleles.jsonl.gz", "wt"))
        lines = stack.enter_context(pilot.stream(args, [args.bcftools, "view", "--no-version", "--threads", args.threads,
            "-r", f"{args.chrom}:{bounds[0]}-{bounds[1]}", "-i", simple.EXTRACT, "-Ov", source], "extract-locus"))
        number = 0
        for line in lines:
            if line.startswith("#"):
                for handle in handles.values():
                    handle.write(line.encode())
                continue
            number += 1
            fields = line.rstrip("\n").split("\t", 8)
            pos, ref, alt, original_id = int(fields[1]), fields[3], fields[4], fields[2]
            kind, size, newpos, newref, newalt, eligible = simple.classify(pos, ref, alt)
            simple.bump(counts["before"], kind, size)
            if not eligible and not simple.final_allele(ref, alt):
                continue
            label = f"TLR{number:07d}"
            fields[2] = label
            if eligible:
                info = pilot.info_parse(fields[7])
                if ("SVTYPE" in info and info["SVTYPE"] != kind) or ("SVLEN" in info and abs(int(info["SVLEN"])) != size):
                    raise ValueError("Simple allele has inconsistent SVTYPE/SVLEN; inspect source on VM")
                fields[1], fields[3], fields[4] = str(newpos), newref, newalt
                if "END" in info:
                    info["END"] = str(newpos + len(newref) - 1)
                    fields[7] = pilot.info_string(info)
                part = kind
                simple.bump(counts["simple"], kind, size)
            else:
                part = "passthrough.sv50"
                simple.bump(counts["passthrough"], kind, size)
            handles[part].write(("\t".join(fields) + "\n").encode())
            written[part] += 1
            meta.write(json.dumps({"label": label, "original_id": original_id, "chrom": fields[0],
                "original_pos": pos, "original_ref": ref, "original_alt": alt, "pos": newpos,
                "ref": newref, "alt": newalt, "kind": kind, "net_bp": size, "simple": eligible,
                "in_focus": pos <= args.end and pos + len(ref) - 1 >= args.start}) + "\n")
    if not sum(written.values()):
        raise ValueError("No SV/indel records found; check callset, region and assembly")
    paths = []
    for kind in parts:
        final = directory / f"{kind}.vcf.gz"
        count = simple.finish(args, directory / f"{kind}.partial.vcf.gz", final, "prepare-" + kind, sort=True)
        if count != written[kind]:
            raise ValueError("Partition counts changed")
        paths.append(final)
    if [pilot.common.identity(p) for p in [source, *indexes]] != origin:
        raise ValueError("Input changed during regional extraction")
    report = {"settings": key, "counts_sv50": counts, "partition_records_all_sizes": dict(written), "samples": len(samples.splitlines())}
    return pilot.manifest(directory, report, simple.artifacts(paths) + [directory / "samples.txt", directory / "alleles.jsonl.gz"])


def memberships(root, kind):
    representatives, lookup = {}, {}
    with pysam.VariantFile(root / kind / "retained.vcf.gz", drop_samples=True) as handle:
        for record in handle:
            lookup[record.id] = record.id
            if "CollapseId" in record.info:
                representatives[str(record.info["CollapseId"])] = record.id
    with pysam.VariantFile(root / kind / "collapsed.vcf.gz", drop_samples=True) as handle:
        for record in handle:
            ids = record.info["MatchId"]
            ids = ids if isinstance(ids, tuple) else ids.split(",")
            matches = [representatives[i] for i in ids if i in representatives]
            if len(matches) != 1:
                raise ValueError("Removed record lacks a unique representative link")
            lookup[record.id] = matches[0]
    return lookup


def pair_scores(task):
    a, b, max_sequence_bp = task
    ma, va, maska = a
    mb, vb, maskb = b
    same_type = ma["kind"] == mb["kind"]
    simple_pair = ma["simple"] and mb["simple"] and same_type
    shared = int((maska & maskb).sum())
    size_similarity = va.sizesim(vb)[0]
    sequence_similarity = None
    payload_similarity = None
    within_limit = max(map(len, (va.ref, va.alts[0], vb.ref, vb.alts[0]))) <= max_sequence_bp
    if same_type and within_limit:
        sequence_similarity = va.seqsim(vb)
        if simple_pair:
            payloads = []
            for record in (va, vb):
                desc, _ = simple.representation.describe(record.ref, record.alts[0])
                seq = record.ref if ma["kind"] == "DEL" else record.alts[0]
                start, suffix = desc["shared_prefix_bp"], desc["shared_suffix_bp"]
                payloads.append(seq[start:len(seq)-suffix if suffix else None])
            distance = edlib.align(*payloads, mode="NW", task="distance")["editDistance"]
            payload_similarity = 1 - distance / max(map(len, payloads))
    ba, ea = va.boundaries()
    bb, eb = vb.boundaries()
    distance_pass = truvari.overlaps(ba-500, ea+500, bb, eb)
    rows = []
    for threshold in (.7, .9):
        for policy in ("current_r500", "same_start"):
            failures = []
            if not ma["simple"] or not mb["simple"]:
                failures.append("non_simple_passthrough")
            if not same_type:
                failures.append("different_type")
            if not distance_pass:
                failures.append("outside_native_refdist")
            if policy == "same_start" and ma["pos"] != mb["pos"]:
                failures.append("different_trimmed_pos")
            if size_similarity < threshold:
                failures.append("size_similarity")
            if same_type and not within_limit:
                failures.append("sequence_not_assessed_length_limit")
            elif same_type and sequence_similarity < threshold:
                failures.append("sequence_similarity")
            if shared:
                failures.append("shared_carriers_gt_all")
            rows.append({"allele_a": ma["label"], "allele_b": mb["label"], "type_a": ma["kind"], "type_b": mb["kind"],
                "original_pos_a": ma["original_pos"], "original_pos_b": mb["original_pos"],
                "trimmed_pos_a": ma["pos"], "trimmed_pos_b": mb["pos"], "net_bp_a": ma["net_bp"], "net_bp_b": mb["net_bp"],
                "native_size_a": va.var_size(), "native_size_b": vb.var_size(), "size_similarity": size_similarity,
                "truvari_sequence_similarity": sequence_similarity, "payload_edit_similarity_max_length": payload_similarity,
                "reciprocal_overlap": va.recovl(vb) if same_type else None,
                "carriers_a": int(maska.sum()), "carriers_b": int(maskb.sum()), "shared_carriers": shared,
                "policy": policy, "threshold": threshold, "pair_passes": not failures,
                "failed_criteria": ";".join(failures) or "none"})
    return rows


def audit(args, prepared):
    directory = args.root / "audit"
    directory.mkdir(exist_ok=True)
    key = {"version": VERSION, "prepare": prepared, "thresholds": [.7, .9], "max_sequence_bp": args.max_sequence_bp,
           "max_focus_alleles": args.max_focus_alleles, "left_aligned": False, "truvari": "5.5.0"}
    cached = pilot.valid_manifest(directory, key)
    if cached:
        print((directory / "summary.json").read_text())
        return
    (directory / "done.json").unlink(missing_ok=True)
    with gzip.open(args.root / "prepare/alleles.jsonl.gz", "rt") as handle:
        metadata = {row["label"]: row for line in handle if (row := json.loads(line))}
    focus = {k for k, v in metadata.items() if v["in_focus"] and v["kind"] in ("INS", "DEL") and v["net_bp"] >= 50}
    if len(focus) > args.max_focus_alleles:
        raise ValueError("Too many focus alleles for an all-pairs quick audit; narrow --region or raise --max-focus-alleles")
    if not focus:
        raise ValueError("No >=50bp literal INS/DEL in the requested focus window")
    cluster_maps = {}
    for policy in ("current_r500", "same_start"):
        local = copy.copy(args)
        local.root = args.root / policy
        (local.root / "logs").mkdir(parents=True, exist_ok=True)
        link = local.root / "prepare"
        if not link.exists():
            link.symlink_to(args.root / "prepare", target_is_directory=True)
        elif link.resolve() != (args.root / "prepare").resolve():
            raise ValueError("Unexpected working preparation path")
        local.same_start = policy == "same_start"
        local.similarities = [.7, .9]
        print(f"Testing regional {policy} at 0.7 and 0.9; native GT constraints enabled.", flush=True)
        with open(args.root / "logs" / f"{policy}-progress.log", "a") as log, redirect_stdout(log):
            simple.merge(local, prepared)
        for threshold in (.7, .9):
            mapping = {k: k for k in metadata if not metadata[k]["simple"]}
            for kind in ("INS", "DEL"):
                mapping.update(memberships(local.root / f"p{round(threshold*100)}", kind))
            cluster_maps[policy, threshold] = mapping
    params = truvari.VariantParams(pctseq=.7, pctsize=.7, refdist=500, pctovl=0, sizemin=0, sizemax=-1,
                                  skip_gt=True, short_circuit=False)
    records = []
    for kind in ("INS", "DEL", "passthrough.sv50"):
        with truvari.VariantFile(args.root / "prepare" / f"{kind}.vcf.gz", params=params) as handle:
            for record in handle:
                if record.id in focus:
                    records.append((metadata[record.id], record, CollapsedCalls.make_genotype_mask(record, "all")))
    pairs = list(combinations(records, 2))
    results = []
    with ThreadPoolExecutor(max_workers=args.threads) as pool:
        for rows in pool.map(pair_scores, ((a, b, args.max_sequence_bp) for a, b in pairs)):
            for row in rows:
                mapping = cluster_maps[row["policy"], row["threshold"]]
                ra, rb = mapping[row["allele_a"]], mapping[row["allele_b"]]
                row.update({"representative_a": ra, "representative_b": rb, "same_local_cluster": ra == rb,
                    "representative_a_passes_final50": metadata[ra]["net_bp"] >= 50,
                    "representative_b_passes_final50": metadata[rb]["net_bp"] >= 50})
                results.append(row)
    if results:
        pilot.tsv(directory / "pair_scores.tsv", results, list(results[0]))
    else:
        (directory / "pair_scores.tsv").write_text("allele_a\tallele_b\tpair_passes\n")
    columns = ("label", "original_id", "chrom", "original_pos", "pos", "kind", "net_bp", "simple")
    allele_rows = [{k: m[k] for k in columns} for m, _, _ in records]
    pilot.tsv(directory / "alleles.tsv", allele_rows, columns)
    assignment_rows = [{"policy": policy, "threshold": threshold, "allele": label, "representative": mapping[label],
        "representative_passes_final50": metadata[mapping[label]]["net_bp"] >= 50}
        for (policy, threshold), mapping in cluster_maps.items() for label in sorted(focus)]
    pilot.tsv(directory / "local_cluster_membership.tsv", assignment_rows, list(assignment_rows[0]))
    summary = {"focus_region": args.region, "context_flank_bp": args.flank, "focus_sv_alleles": len(focus),
        "samples": prepared["key"]["samples"], "focus_types": dict(Counter(metadata[k]["kind"] for k in focus)),
        "focus_simple": sum(metadata[k]["simple"] for k in focus), "distinct_pairs": len(pairs),
        "comparisons": [{"policy": policy, "threshold": threshold,
            "pairs_pass": sum(r["pair_passes"] for r in results if r["policy"] == policy and r["threshold"] == threshold),
            "pairs_same_local_cluster": sum(r["same_local_cluster"] for r in results if r["policy"] == policy and r["threshold"] == threshold)}
            for policy in ("current_r500", "same_start") for threshold in (.7, .9)],
        "notes": ["Working copies only; original callset is unchanged. Unique audit IDs map to original alleles in prepare/alleles.jsonl.gz.",
                  "Current trimming-only policies; reference-based left-alignment is not performed.",
                  "Pair pass is direct eligibility. Regional cluster outcome also depends on representatives and other alleles.",
                  "Local collapse includes all simple indel sizes in the window and flanks. Whole-chromosome outcomes can differ through additional outside records.",
                  "Genotype constraints use all samples in this input, not just plotted introgressed haplotypes.",
                  "The figure does not identify alleles unambiguously; match local allele positions/IDs to the plot's source records.",
                  "Payload edit similarity uses 1-distance/max_length; Truvari sequence similarity uses native scoring. They are different metrics."]}
    pilot.common.json_write(directory / "summary.json", summary)
    pilot.manifest(directory, key, [directory / name for name in ("summary.json", "pair_scores.tsv", "alleles.tsv", "local_cluster_membership.tsv")])
    print(json.dumps(summary, indent=2))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--input-bcf", type=Path)
    parser.add_argument("--source-uri", default=SOURCE)
    parser.add_argument("--region", default="chr4:38660000-38960000")
    parser.add_argument("--flank", type=int, default=2000)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--memory-gb", type=int, default=26)
    parser.add_argument("--max-focus-alleles", type=int, default=1000)
    parser.add_argument("--max-sequence-bp", type=int, default=100000)
    parser.add_argument("--phase", choices=("prepare", "audit", "all"), default="all")
    parser.add_argument("--bcftools", default=os.environ.get("BCFTOOLS_BIN", "bcftools"))
    options = parser.parse_args(argv)
    match = re.fullmatch(r"(chr\d+):(\d+)-(\d+)", options.region)
    if not match or int(match[2]) < 1 or int(match[3]) < int(match[2]) or options.flank < 0:
        parser.error("Expected chrN:start-end and a nonnegative flank")
    if not 1 <= options.threads <= 4 or options.max_focus_alleles < 1 or options.max_sequence_bp < 1:
        parser.error("Use 1–4 threads and positive audit limits")
    args = pilot.parse_args(["--root", str(options.root), "--chrom", match[1], "--source-uri", options.source_uri,
        "--confirm-source-uri", options.source_uri, "--threads", str(options.threads), "--memory-gb", str(options.memory_gb),
        "--bcftools", options.bcftools])
    for name in ("region", "flank", "max_focus_alleles", "max_sequence_bp", "phase"):
        setattr(args, name, getattr(options, name))
    args.start, args.end = int(match[2]), int(match[3])
    args.input_bcf = options.input_bcf.resolve() if options.input_bcf else None
    (args.root / "logs").mkdir(parents=True, exist_ok=True)
    logging.basicConfig(filename=args.root / "logs/audit.log", level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s", force=True)
    import fcntl
    with open(args.root / ".audit.lock", "a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if importlib.metadata.version("truvari") != "5.5.0":
            raise ValueError("Truvari 5.5.0 is required")
        logging.info("PROVENANCE parameters=%s seed=%s python=%s script_sha256=%s truvari=5.5.0", vars(args), pilot.SEED,
                     sys.version, pilot.common.digest(__file__))
        pilot.command(args, [args.bcftools, "--version"], "bcftools-version")
        try:
            prepared = prepare(args) if args.phase in {"prepare", "all"} else pilot.require_manifest(args.root / "prepare")
            settings = prepared["key"]["settings"]
            if (settings["chrom"], settings["start"], settings["end"], settings["flank"]) != (args.chrom, args.start, args.end, args.flank):
                raise ValueError("Region differs from prepared extraction; run phase all")
            if args.phase in {"audit", "all"}:
                audit(args, prepared)
        except Exception:
            logging.exception("TLR audit failed")
            raise SystemExit(f"Audit failed; inspect {args.root / 'logs/audit.log'} on the VM") from None


if __name__ == "__main__":
    main()
