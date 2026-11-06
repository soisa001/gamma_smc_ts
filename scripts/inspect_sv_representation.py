#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = ["google-crc32c==1.8.0"]
# ///
"""Stage and inspect popped chr1 alleles; no merging, uploads, or sample output."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import logging
import math
import os
from pathlib import Path
import re
import sys

import run_chr1_sv_merge as pilot

VERSION = 1
SOURCE = ("gs://rw-long-reads-transfer-2026-06-17/v9/lrWGS/panel/panel/"
          "panel_popped_vcf/aou_lr_phase2_v1.chr1.popped.bcf")


def describe(ref, alt):
    if alt in {".", "*", "<NON_REF>", "<NONREF>"}:
        return None, "missing_or_spanning"
    if alt.startswith("<") or "[" in alt or "]" in alt or alt.startswith(".") or alt.endswith("."):
        return None, "symbolic_or_breakend"
    ref, alt = ref.upper(), alt.upper()
    delta = len(alt) - len(ref)
    kind = "INS" if delta > 0 else "DEL" if delta < 0 else "SNP" if len(ref) == 1 else "REPL"
    prefix = 0
    limit = min(len(ref), len(alt))
    while prefix < limit and ref[prefix] == alt[prefix]:
        prefix += 1
    suffix = 0
    while suffix < limit - prefix and ref[-suffix - 1] == alt[-suffix - 1]:
        suffix += 1
    rleft, aleft = len(ref) - prefix - suffix, len(alt) - prefix - suffix
    if not rleft and not aleft:
        kind = "REF_EQUAL"
    return {"type": kind, "net_bp": abs(delta), "ref_bp": len(ref), "alt_bp": len(alt),
            "shared_prefix_bp": prefix, "shared_suffix_bp": suffix,
            "shared_flank_bp": prefix + suffix, "residual_ref_bp": rleft, "residual_alt_bp": aleft,
            "shared_fraction_per_mille": round(1000 * (prefix + suffix) / limit) if limit else 0,
            "complex": bool(rleft and aleft),
            "ambiguous": bool(set(ref + alt) - set("ACGT"))}, None


def quantiles(hist):
    n = sum(hist.values())
    if not n:
        return {key: None for key in ("min", "median", "p95", "max", "mean")}
    ordered = sorted(hist.items())
    def at(index):
        count = 0
        for value, frequency in ordered:
            count += frequency
            if count > index:
                return value
        raise AssertionError("Invalid histogram rank")
    def q(fraction):
        rank = (n - 1) * fraction
        lo, hi = math.floor(rank), math.ceil(rank)
        return at(lo) + (at(hi) - at(lo)) * (rank - lo)
    return {"min": ordered[0][0], "median": q(.5), "p95": q(.95), "max": ordered[-1][0],
            "mean": sum(value * count for value, count in ordered) / n}


def summarize(lines, chrom):
    sites, alt_counts, groups, histograms = Counter(), Counter(), defaultdict(Counter), {}
    # Track exact duplicates and shared REF groups at one POS at a time.
    locus = defaultdict(lambda: defaultdict(Counter))
    previous = None

    def flush():
        for group, refs in locus.items():
            counts = groups[group]
            counts["distinct_positions"] += 1
            counts["distinct_pos_ref_groups"] += len(refs)
            for alts in refs.values():
                n = sum(alts.values())
                counts["unique_alleles"] += len(alts)
                counts["exact_duplicate_alleles"] += n - len(alts)
                if len(alts) > 1:
                    counts["pos_ref_groups_with_multiple_unique_alts"] += 1
                    counts["alleles_in_multiallelic_pos_ref_groups"] += n
        locus.clear()

    for line in lines:
        contig, pos, ref, alt_string = line.rstrip("\r\n").split("\t")
        pos = int(pos)
        if contig != chrom or (previous is not None and pos < previous):
            raise ValueError("Expected sorted records from the requested chromosome")
        if previous != pos:
            flush()
        previous = pos
        sites["records"] += 1
        alts = alt_string.split(",")
        alt_counts[len(alts)] += 1
        sites["multiallelic_records"] += len(alts) > 1
        ref_hash = hashlib.sha256(ref.encode()).hexdigest()
        for alt in alts:
            sites["alt_entries"] += 1
            row, excluded = describe(ref, alt)
            if excluded:
                sites[excluded] += 1
                continue
            sites["literal_alt_entries"] += 1
            scopes = ["all_literal"]
            if row["net_bp"] >= 50:
                scopes.append("net_ge50")
            elif row["complex"] and max(row["residual_ref_bp"], row["residual_alt_bp"]) >= 50:
                sites["complex_net_lt50_but_residual_span_ge50"] += 1
            for scope in scopes:
                group = (scope, row["type"])
                counts = groups[group]
                counts["alleles"] += 1
                counts["complex_alleles"] += row["complex"]
                counts["simple_indel_alleles"] += row["type"] in {"INS", "DEL"} and not row["complex"]
                counts["ambiguous_alleles"] += row["ambiguous"]
                counts["both_ref_alt_gt1"] += row["ref_bp"] > 1 and row["alt_bp"] > 1
                counts["shared_flank_gt1"] += row["shared_flank_bp"] > 1
                for field, value in row.items():
                    if field not in {"type", "complex", "ambiguous"}:
                        histograms.setdefault((*group, field), Counter())[value] += 1
                locus[group][ref_hash][hashlib.sha256(alt.encode()).hexdigest()] += 1
        if sites["records"] % 100000 == 0:
            logging.info("PROGRESS records=%d literal_alts=%d", sites["records"], sites["literal_alt_entries"])
    flush()
    rows, distributions = [], []
    for (scope, kind), counts in sorted(groups.items()):
        n = counts["alleles"]
        row = {"scope": scope, "type": kind, **{field: counts[field] for field in (
            "alleles", "unique_alleles", "exact_duplicate_alleles", "simple_indel_alleles", "complex_alleles",
            "ambiguous_alleles", "both_ref_alt_gt1", "shared_flank_gt1", "distinct_positions",
            "distinct_pos_ref_groups", "pos_ref_groups_with_multiple_unique_alts",
            "alleles_in_multiallelic_pos_ref_groups")},
            "complex_percent": 100 * counts["complex_alleles"] / n,
            "shared_flank_gt1_percent": 100 * counts["shared_flank_gt1"] / n}
        for (s, t, field), hist in histograms.items():
            if (s, t) != (scope, kind):
                continue
            stats = quantiles(hist)
            distributions.append({"scope": scope, "type": kind, "metric": field, **stats})
            if field in {"ref_bp", "alt_bp", "shared_flank_bp"}:
                row[field + "_median"] = stats["median"]
        rows.append(row)
    return {"site_counts": dict(sites), "alts_per_record": quantiles(alt_counts),
            "allele_summary": rows, "distributions": distributions}


def inspect(args, label, path):
    directory = args.root / label
    (directory / "logs").mkdir(parents=True, exist_ok=True)
    local = argparse.Namespace(root=directory, bcftools=args.bcftools)
    snapshot = pilot.common.identity(path)
    indexes = [Path(str(path) + suffix) for suffix in (".csi", ".tbi") if Path(str(path) + suffix).is_file()]
    if not indexes:
        raise ValueError(f"An existing index is required (input is not modified): {path}")
    index_snapshots = [pilot.common.identity(index) for index in indexes]
    key = {"schema": VERSION, "input": snapshot, "chrom": "chr1",
           "indexes": index_snapshots,
           "bcftools": pilot.command(local, [args.bcftools, "--version"], "version", True).splitlines()[0]}
    if pilot.valid_manifest(directory, key):
        logging.info("Reuse completed representation report: %s", label)
        return json.loads((directory / "summary.json").read_text())
    (directory / "done.json").unlink(missing_ok=True)
    header = pilot.command(local, [args.bcftools, "view", "-h", path], "header", True)
    (directory / "header.txt").write_text(header)
    cmd = [args.bcftools, "query", "-r", "chr1", "-f", "%CHROM\t%POS\t%REF\t%ALT\n", path]
    with pilot.stream(local, cmd, "scan") as lines:
        result = summarize(lines, "chr1")
    if not result["site_counts"].get("records"):
        raise ValueError(f"No chr1 records: {path}")
    if snapshot != pilot.common.identity(path) or index_snapshots != [pilot.common.identity(index) for index in indexes]:
        raise ValueError("Input changed while scanning")
    result.update({"label": label, "input": str(path), "seed": pilot.SEED,
                   "header_info_ids": re.findall(r"^##INFO=<ID=([^,>]+)", header, re.M),
                   "definitions": {
                       "scope": "all_literal covers every literal ALT; net_ge50 uses abs(len(ALT)-len(REF)) >= 50 only for reporting, not input selection",
                       "complex": "Both REF and ALT have residual bases after maximal identical prefix then suffix trimming; not proof of multiple biological events",
                       "padding": "Shared flank bp is removable common prefix plus suffix, including any anchor; shared_fraction_per_mille divides this by min(REF length, ALT length). A simple anchor can give 1000 per mille, so use shared_flank_gt1 and absolute bp to assess substantial padding. Does not measure shared internal sequence or ALT-to-ALT similarity",
                       "groups": "CHROM/POS/REF groups, not proof of original graph-bubble identity",
                       "validation": "Indexed chr1 site-field scan to EOF; unchanged input identity; no sample or genotype analysis",
                       "comparison": "An already split/size-filtered baseline has a different record structure; compare net_ge50 allele metrics, not raw multiallelic-record percentages"}})
    pilot.common.json_write(directory / "summary.json", result)
    artifacts = [directory / "summary.json", directory / "header.txt"]
    for name, rows in (("allele_summary.tsv", result["allele_summary"]), ("length_distributions.tsv", result["distributions"])):
        if rows:
            pilot.tsv(directory / name, rows, list(rows[0]))
            artifacts.append(directory / name)
    pilot.manifest(directory, key, artifacts)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--phase", choices=("stage", "inspect", "all"), default="all")
    parser.add_argument("--input-bcf", type=Path, help="Local indexed popped BCF/VCF instead of staging")
    parser.add_argument("--compare-input", type=Path, help="Optional local indexed bubble-split baseline; never changed")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--memory-gb", type=int, default=26)
    parser.add_argument("--bcftools", default="bcftools")
    parser.add_argument("--gcloud", default="gcloud")
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("threads must be positive")
    args.root = args.root.resolve()
    (args.root / "logs").mkdir(parents=True, exist_ok=True)
    logging.basicConfig(filename=args.root / "inspection.log", level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s", force=True)
    logging.info("PARAMETERS %s", json.dumps(vars(args), default=str, sort_keys=True))
    logging.info("PROVENANCE python=%s seed=%s script_sha256=%s", sys.version, pilot.SEED, pilot.common.digest(__file__))
    import fcntl
    with (args.root / ".inspection.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        path = args.input_bcf.resolve() if args.input_bcf else args.root / "inputs" / SOURCE.rsplit("/", 1)[1]
        if not args.input_bcf:
            if args.phase in {"stage", "all"}:
                stage_args = pilot.parse_args(["--root", str(args.root), "--source-uri", SOURCE,
                    "--confirm-source-uri", SOURCE, "--threads", str(args.threads),
                    "--memory-gb", str(args.memory_gb), "--scratch-reserve-gb", "4",
                    "--bcftools", args.bcftools, "--gcloud", args.gcloud])
                print("Staging checksum-verified popped chr1 BCF/index; transfer log stays on VM.", flush=True)
                pilot.stage(stage_args)
            else:
                pilot.require_manifest(args.root / "inputs")
        if args.phase == "stage":
            return
        inputs = [("popped", path)]
        if args.compare_input:
            inputs.append(("bubble_baseline", args.compare_input.resolve()))
        print("Scanning allele strings; no genotypes or individual variants will be printed.", flush=True)
        with ThreadPoolExecutor(max_workers=min(args.threads, len(inputs))) as pool:
            reports = list(pool.map(lambda item: inspect(args, *item), inputs))
        combined = [{"dataset": report["label"], **row} for report in reports
                    for row in report["allele_summary"] if row["scope"] == "net_ge50"]
        if combined:
            pilot.tsv(args.root / "net_ge50_comparison.tsv", combined, list(combined[0]))
            print((args.root / "net_ge50_comparison.tsv").read_text(), end="")
        else:
            (args.root / "net_ge50_comparison.tsv").unlink(missing_ok=True)
            print("No literal INS/DEL alleles with net length >=50 bp in the inspected inputs.")
        for report in reports:
            print(json.dumps({"dataset": report["label"], "site_counts": report["site_counts"],
                              "alts_per_record": report["alts_per_record"]}))
        print(f"Reports and headers saved on VM: {args.root}")


if __name__ == "__main__":
    main()
