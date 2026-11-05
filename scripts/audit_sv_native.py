#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = []
# ///
"""Read-only site/cluster audit of completed native Truvari output; no GT export."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
import logging
import os
from pathlib import Path

import run_chr1_sv_merge as pilot

VERSION = 1
TAGS = ("SVTYPE", "SVLEN", "CollapseId", "MatchId", "NumCollapsed", "PctSeqSimilarity", "PctSizeSimilarity")


def number(value):
    return None if value in (".", "") else float(value.split(",")[0])


def records(args, path, label):
    fmt = "%CHROM\t%POS\t%REF\t%ALT\t" + "\t".join("%INFO/" + tag for tag in TAGS) + "\n"
    with pilot.stream(args, [args.bcftools, "query", "-u", "-f", fmt, path], label) as handle:
        for line in handle:
            chrom, pos, ref, alt, *values = line.rstrip("\r\n").split("\t")
            if chrom != args.chrom or "," in alt:
                raise ValueError("Unexpected chromosome or multiallelic record")
            info = dict(zip(TAGS, values))
            literal = alt not in {".", "*"} and not any(c in alt for c in "<>[].")
            delta = len(alt) - len(ref) if literal else 0
            kind = ("INS" if delta > 0 else "DEL" if delta < 0 else "OTHER") if literal else "OTHER"
            parsed = pilot.allele(ref, alt, int(pos), {}, minimum=1, maximum=-1)
            raw = abs(delta) if literal else None
            complex_allele = bool(kind in {"INS", "DEL"} and parsed and parsed[3] - parsed[2] != (raw if kind == "DEL" else 0))
            rh, ah = hashlib.sha256(ref.encode()).hexdigest(), hashlib.sha256(alt.encode()).hexdigest()
            yield {"key": (chrom, int(pos), rh, ah), "pos": int(pos), "ref": rh, "alt": ah,
                "kind": kind, "complex": complex_allele, "raw": raw,
                "ref_len": len(ref), "alt_len": len(alt), "info": info}


def run(args):
    similarity = getattr(args, "similarity", .7)
    tag = f"p{round(similarity * 100)}"
    base = args.batch_root / args.chrom
    paths = {"before": base / "prepare/input.vcf.gz", "after": base / f"{tag}/{args.chrom}.vcf.gz",
             "removed": base / f"{tag}/{args.chrom}.collapsed.vcf.gz"}
    manifests = {"prepare": base / "prepare/done.json", "merge": base / tag / "done.json"}
    for path in (*paths.values(), *manifests.values()):
        if not path.is_file():
            raise ValueError(f"Completed input/output required: {path}")
    metadata = {name: json.loads(path.read_text()) for name, path in manifests.items()}
    if metadata["merge"]["key"]["prepare"] != metadata["prepare"]:
        raise ValueError("Preparation and merge manifests disagree")
    if metadata["merge"]["key"].get("chain"):
        raise ValueError("This representative-pair audit expects --chain off")
    snapshots = {name: pilot.common.identity(path) for name, path in {**paths, **manifests}.items()}
    cache = {"version": VERSION, "inputs": snapshots, "manifests": metadata, "seed": pilot.SEED}
    (args.root / "logs").mkdir(parents=True, exist_ok=True)
    logging.basicConfig(filename=args.root / "audit.log", level=logging.INFO, force=True,
                        format="%(asctime)s %(levelname)s %(message)s")
    logging.info("PARAMETERS %s", vars(args))
    logging.info("PROVENANCE sha256=%s seed=%s", pilot.common.digest(__file__), pilot.SEED)
    if pilot.valid_manifest(args.root, cache):
        return json.loads((args.root / "audit.json").read_text())
    (args.root / "done.json").unlink(missing_ok=True)
    original, accounted = Counter(), Counter()
    categories = []
    representatives, groups = {}, {}
    pair_counts = defaultdict(Counter)
    pair_values = defaultdict(lambda: defaultdict(list))
    totals = {}
    for role in ("before", "after", "removed"):
        print(f"Scanning {args.chrom} {role} site fields...", flush=True)
        counts, positions, pos_refs = defaultdict(Counter), defaultdict(set), defaultdict(set)
        total = 0
        for row in records(args, paths[role], role):
            total += 1
            (original if role == "before" else accounted)[row["key"]] += 1
            kind, info = row["kind"], row["info"]
            counts[kind]["records"] += 1
            counts[kind]["complex_records"] += row["complex"]
            counts[kind]["ref_length_gt_1"] += row["ref_len"] > 1
            counts[kind]["alt_length_gt_1"] += row["alt_len"] > 1
            positions[kind].add(row["pos"])
            pos_refs[kind].add((row["pos"], row["ref"]))
            if kind in {"INS", "DEL"}:
                if info["SVTYPE"] != ".":
                    counts[kind]["svtype_present"] += 1
                    counts[kind]["svtype_disagrees_with_net_length"] += info["SVTYPE"].split(",")[0] != kind
                if info["SVLEN"] != ".":
                    counts[kind]["svlen_present"] += 1
                    counts[kind]["svlen_disagrees_with_raw_length"] += abs(number(info["SVLEN"])) != row["raw"]
            if role == "after" and info["CollapseId"] != ".":
                cid = info["CollapseId"]
                if cid in representatives:
                    raise ValueError("Duplicate representative CollapseId")
                representatives[cid] = row
                groups[cid] = {"collapse_id": cid, "type": kind, "representative_pos": row["pos"],
                    "expected_removed": int(number(info["NumCollapsed"]) or 0), "observed_removed": 0,
                    "min_pos": row["pos"], "max_pos": row["pos"],
                    "same_pos_ref_removed": 0, "complex_removed": 0}
            if role != "removed":
                continue
            matches = [cid for cid in info["MatchId"].split(",") if cid in representatives]
            if len(matches) != 1:
                raise ValueError("Removed record has no unique representative MatchId")
            cid = matches[0]
            rep, group = representatives[cid], groups[cid]
            group["observed_removed"] += 1
            group["min_pos"] = min(group["min_pos"], row["pos"])
            group["max_pos"] = max(group["max_pos"], row["pos"])
            same = row["pos"] == rep["pos"] and row["ref"] == rep["ref"]
            group["same_pos_ref_removed"] += same
            group["complex_removed"] += row["complex"]
            tally = pair_counts[kind]
            tally["pairs"] += 1
            tally["same_pos_and_ref"] += same
            tally["same_pos_ref_but_different_alt"] += same and row["alt"] != rep["alt"]
            tally["complex_same_pos_ref_different_alt"] += same and row["complex"] and row["alt"] != rep["alt"]
            tally["representative_net_type_differs"] += kind != rep["kind"]
            pair_values[kind]["abs_pos_distance_bp"].append(abs(row["pos"] - rep["pos"]))
            if row["raw"] and rep["raw"]:
                ratio = min(row["raw"], rep["raw"]) / max(row["raw"], rep["raw"])
                pair_values[kind]["raw_length_ratio"].append(ratio)
                tally["raw_length_ratio_below_0.7"] += ratio < .7 - 1e-8
            for tag in ("PctSeqSimilarity", "PctSizeSimilarity"):
                val = number(info[tag])
                if val is None:
                    tally[tag + "_missing"] += 1
                    continue
                pair_values[kind][tag].append(val)
                # Scores in the VCF are rounded to four decimals.
                tally[tag + "_below_0.7"] += val < .7 - .00005
                tally[tag + "_below_0.9"] += val < .9 - .00005
                tally[tag + "_equal_1"] += val == 1
        totals[role] = total
        for kind, values in sorted(counts.items()):
            categories.append({"callset": role, "type": kind, **dict(values),
                "unique_positions": len(positions[kind]), "unique_pos_ref_groups": len(pos_refs[kind])})
    cluster_rows = []
    for group in groups.values():
        group["members"] = group["observed_removed"] + 1
        group["position_span_bp"] = group["max_pos"] - group["min_pos"]
        cluster_rows.append(group)
    cluster_rows.sort(key=lambda r: (-r["members"], r["representative_pos"]))
    summary = {"chrom": args.chrom, "record_counts": totals,
        "exact_site_multiset_conserved": original == accounted,
        "before_count_matches_manifest": totals["before"] == metadata["prepare"]["key"]["records"],
        "numcollapsed_matches_membership": all(g["expected_removed"] == g["observed_removed"] for g in cluster_rows),
        "representation_summary": categories,
        "removed_vs_representative": {kind: {"counts": dict(pair_counts[kind]),
            "distributions": {name: pilot.quantiles(values) for name, values in pair_values[kind].items()}}
            for kind in sorted(pair_counts)},
        "merged_cluster_sizes": {kind: pilot.quantiles([g["members"] for g in cluster_rows
            if g["type"] == kind and g["observed_removed"]]) for kind in ("INS", "DEL")},
        "merged_cluster_position_span_bp": {kind: pilot.quantiles([g["position_span_bp"] for g in cluster_rows
            if g["type"] == kind and g["observed_removed"]]) for kind in ("INS", "DEL")},
        "limitations": "Site audit only, not genotype transfer or biological equivalence validation. Stored scores are native rounded annotations. POS/REF groups are representation groups, not proven biological loci."}
    paths_out = [args.root / "audit.json"]
    pilot.common.json_write(paths_out[0], summary)
    if cluster_rows:
        path = args.root / "clusters.tsv"
        pilot.tsv(path, cluster_rows, list(cluster_rows[0]))
        paths_out.append(path)
    if snapshots != {name: pilot.common.identity(path) for name, path in {**paths, **manifests}.items()}:
        raise ValueError("Input changed during audit")
    pilot.manifest(args.root, cache, paths_out)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch-root", type=Path, required=True)
    parser.add_argument("--chrom", default="chr8", choices=[f"chr{i}" for i in range(1, 22)])
    parser.add_argument("--bcftools", default=os.environ.get("BCFTOOLS_BIN", "bcftools"))
    parser.add_argument("--similarity", type=float, choices=(.7, .9), default=.7)
    args = parser.parse_args()
    args.batch_root = args.batch_root.resolve()
    args.root = args.batch_root / args.chrom / f"audit_native_p{round(args.similarity * 100)}"
    args.root.mkdir(parents=True, exist_ok=True)
    import fcntl
    with (args.root / ".audit.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        summary = run(args)
    print(json.dumps(summary, indent=2))
    print(f"Detailed cluster table stays on VM: {args.root / 'clusters.tsv'}")


if __name__ == "__main__":
    main()
