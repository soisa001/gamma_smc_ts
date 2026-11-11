#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = ["google-crc32c==1.8.0"]
# ///
"""Compare current vs same-start runs; site fields only, aggregate output only."""
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from contextlib import ExitStack
import argparse
import hashlib
import json
import logging
import os
from pathlib import Path
import time

import inspect_sv_representation as representation
import run_chr1_sv_merge as pilot


def union(intervals):
    result = []
    for start, end in sorted(intervals):
        if result and start <= result[-1][1]:
            result[-1] = (result[-1][0], max(end, result[-1][1]))
        else:
            result.append((start, end))
    return result


def coverage(intervals):
    return sum(end - start for start, end in intervals)


def intersection_bp(left, right):
    i = j = total = 0
    while i < len(left) and j < len(right):
        total += max(0, min(left[i][1], right[j][1]) - max(left[i][0], right[j][0]))
        if left[i][1] < right[j][1]:
            i += 1
        else:
            j += 1
    return total


def scan(args, paths, label):
    stats = {kind: Counter({metric: 0 for metric in ("sv_records", "raw_bp_sum_records", "unique_alleles",
                          "raw_bp_sum_unique_alleles", "exact_duplicate_records")}) for kind in ("INS", "DEL")}
    stats["DEL"]["alleles_excluded_from_deleted_union"] = 0
    keys, lengths, intervals = Counter(), {}, []
    for index, path in enumerate(paths):
        snapshot = pilot.common.identity(path)
        with pilot.stream(args, [args.bcftools, "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n", path], f"{label}-{index}") as lines:
            for line in lines:
                chrom, pos, ref, alt = line.rstrip("\n").split("\t")
                if chrom != "chr1" or "," in alt:
                    raise ValueError("Comparison requires biallelic chr1 records")
                desc, _ = representation.describe(ref, alt)
                if not desc or desc["type"] not in stats or desc["net_bp"] < 50:
                    continue
                kind, size = desc["type"], desc["net_bp"]
                key = (kind, int(pos), hashlib.sha256(ref.encode()).hexdigest(), hashlib.sha256(alt.encode()).hexdigest())
                keys[key] += 1
                lengths[key] = size
                stats[kind]["sv_records"] += 1
                stats[kind]["raw_bp_sum_records"] += size
                if kind == "DEL":
                    if desc["complex"] or desc["ambiguous"]:
                        stats[kind]["alleles_excluded_from_deleted_union"] += 1
                    else:
                        start = int(pos) - 1 + desc["shared_prefix_bp"]
                        intervals.append((start, start + desc["residual_ref_bp"]))
        if pilot.common.identity(path) != snapshot:
            raise ValueError("A compared file changed while scanning")
    for key, size in lengths.items():
        stats[key[0]]["unique_alleles"] += 1
        stats[key[0]]["raw_bp_sum_unique_alleles"] += size
    for kind in stats:
        stats[kind]["exact_duplicate_records"] = stats[kind]["sv_records"] - stats[kind]["unique_alleles"]
    intervals = union(intervals)
    stats["DEL"]["unique_deleted_reference_bp"] = coverage(intervals)
    stats["INS"]["unique_inserted_allele_bp"] = stats["INS"]["raw_bp_sum_unique_alleles"]
    return {"stats": stats, "keys": keys, "union": intervals}


def percentage(removed, before):
    return 100 * removed / before if before else 0.0


def compare(args):
    manifests = {}
    prepared = pilot.require_manifest(args.current_root / "prepare")
    if pilot.require_manifest(args.same_start_root / "prepare") != prepared:
        raise ValueError("Runs do not share exactly the same verified preparation")
    paths = {"before": [args.current_root / "prepare" / f"{kind}.vcf.gz"
                         for kind in ("INS", "DEL", "passthrough.sv50")]}
    for policy, root in (("current_r500", args.current_root), ("same_start", args.same_start_root)):
        for similarity in (.7, .9):
            tag = f"p{round(similarity*100)}"
            manifest = pilot.require_manifest(root / tag)
            if manifest["key"]["prepare"] != prepared:
                raise ValueError("Final output and preparation manifests disagree")
            expected_policy = "same_trimmed_CHROM_POS_v1" if policy == "same_start" else None
            if manifest["key"].get("candidate_grouping") != expected_policy:
                raise ValueError("Unexpected merge policy; check --current-root and --same-start-root")
            if manifest["key"]["similarity"] != similarity:
                raise ValueError("Unexpected threshold in merge manifest")
            name = policy + "_" + tag
            manifests[name] = manifest
            paths[name] = [root / tag / "chr1.sv50.vcf.gz"]
    key = {"version": 1, "prepare": prepared, "merges": manifests, "min_net_bp": 50,
           "paths": {name: list(map(str, values)) for name, values in paths.items()}}
    if pilot.valid_manifest(args.root, key):
        print((args.root / "comparison.tsv").read_text(), end="")
        return
    (args.root / "done.json").unlink(missing_ok=True)
    print("Scanning site fields for all four results and their shared baseline; no sample output.", flush=True)
    with ThreadPoolExecutor(max_workers=args.threads) as pool:
        results = dict(zip(paths, pool.map(lambda item: scan(args, item[1], item[0]), paths.items())))
    original = results["before"]
    expected = prepared["key"]["counts_sv50"]["before"]
    for kind in ("INS", "DEL"):
        if original["stats"][kind]["sv_records"] != expected.get(kind, 0) or original["stats"][kind]["raw_bp_sum_records"] != expected.get(kind + "_bp", 0):
            raise ValueError("Baseline counts or bp differ from preparation manifest")
    metrics, overview = [], []
    for name, result in results.items():
        if name == "before":
            continue
        if result["keys"] - original["keys"]:
            raise ValueError("Final alleles are not a submultiset of the prepared baseline")
        policy, tag = name.rsplit("_", 1)
        similarity = int(tag[1:]) / 100
        shared = intersection_bp(original["union"], result["union"])
        if coverage(result["union"]) != shared:
            raise ValueError("Final deletion union gained reference bases absent from baseline")
        for kind in ("INS", "DEL", "TOTAL"):
            if kind == "TOTAL":
                before = sum((original["stats"][k] for k in ("INS", "DEL")), Counter())
                after = sum((result["stats"][k] for k in ("INS", "DEL")), Counter())
                # No meaningful combined INS-sequence/DEL-reference union metric.
                for counts in (before, after):
                    counts.pop("unique_deleted_reference_bp", None)
                    counts.pop("unique_inserted_allele_bp", None)
            else:
                before, after = original["stats"][kind], result["stats"][kind]
            for metric in sorted(set(before) | set(after)):
                b, a = before[metric], after[metric]
                metrics.append({"policy": policy, "similarity": similarity, "type": kind,
                    "metric": metric, "before": b, "after": a, "removed": b-a, "percent_removed": percentage(b-a, b)})
            ub = "unique_deleted_reference_bp" if kind == "DEL" else "unique_inserted_allele_bp"
            b, a = before["sv_records"], after["sv_records"]
            overview.append({"policy": policy, "similarity": similarity, "type": kind,
                "before_alleles": b, "after_alleles": a, "alleles_removed": b-a,
                "percent_alleles_removed": percentage(b-a, b),
                "unique_bp_definition": ub if kind != "TOTAL" else "not_combined",
                "before_unique_bp": before[ub] if kind != "TOTAL" else "",
                "after_unique_bp": after[ub] if kind != "TOTAL" else "",
                "unique_bp_removed": before[ub]-after[ub] if kind != "TOTAL" else "",
                "percent_unique_bp_removed": percentage(before[ub]-after[ub], before[ub]) if kind != "TOTAL" else ""})
    pilot.tsv(args.root / "comparison.tsv", overview, list(overview[0]))
    pilot.tsv(args.root / "metrics.tsv", metrics, list(metrics[0]))
    paired = []
    lookup = {(row["policy"], row["similarity"], row["type"], row["metric"]): row for row in metrics}
    for row in metrics:
        if row["policy"] != "same_start":
            continue
        current = lookup["current_r500", row["similarity"], row["type"], row["metric"]]
        paired.append({"similarity": row["similarity"], "type": row["type"], "metric": row["metric"],
            "before": row["before"], "after_current_r500": current["after"], "after_same_start": row["after"],
            "same_start_minus_current": row["after"]-current["after"]})
    pilot.tsv(args.root / "same_start_vs_current.tsv", paired, list(paired[0]))
    pilot.common.json_write(args.root / "definitions.json", {
        "seed": pilot.SEED, "parameters": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
        "scope": "Literal net >=50bp INS/DEL, including unchanged replacements; prepared representation used on both sides",
        "unique_alleles": "Exact CHROM/POS/REF/ALT after preparation, each counted once",
        "unique_deleted_reference_bp": "Union of pure unambiguous deletion intervals after removing shared flanks; anchor excluded",
        "unique_bp_removed_DEL": "Before union minus after union; verified no newly covered reference bp",
        "unique_inserted_allele_bp": "Sum abs(len(ALT)-len(REF)) once per exact insertion allele; not nonredundant sequence",
        "complex_ambiguous_DEL": "Excluded from reference union and explicitly counted in metrics.tsv",
        "input_snapshots": key})
    pilot.manifest(args.root, key, [args.root / name for name in
        ("comparison.tsv", "metrics.tsv", "same_start_vs_current.tsv", "definitions.json")])
    print((args.root / "comparison.tsv").read_text(), end="")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--current-root", type=Path, required=True)
    parser.add_argument("--same-start-root", type=Path, required=True)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--bcftools", default=os.environ.get("BCFTOOLS_BIN", "bcftools"))
    args = parser.parse_args(argv)
    if not 1 <= args.threads <= 4:
        parser.error("Use 1–4 scan threads")
    for name in ("root", "current_root", "same_start_root"):
        setattr(args, name, getattr(args, name).resolve())
    if args.current_root == args.same_start_root or args.root in {args.current_root, args.same_start_root}:
        parser.error("Use distinct run and report directories")
    (args.root / "logs").mkdir(parents=True, exist_ok=True)
    logging.basicConfig(filename=args.root / "logs/comparison.log", level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s", force=True)
    logging.info("PROVENANCE parameters=%s seed=%s script_sha256=%s time=%s", vars(args), pilot.SEED,
                 pilot.common.digest(__file__), time.time())
    import fcntl
    with ExitStack() as stack:
        for root, mode in ((args.current_root, fcntl.LOCK_SH), (args.same_start_root, fcntl.LOCK_SH), (args.root, fcntl.LOCK_EX)):
            lock = stack.enter_context(open(root / ".simple_merge.lock", "a"))
            fcntl.flock(lock, mode | fcntl.LOCK_NB)
        compare(args)


if __name__ == "__main__":
    main()
