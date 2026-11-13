"""Read-only audit of prepared >10 kb DEL/INS pairs using Truvari 5.5.0.

Run with the existing uv-built Truvari environment. All pair details stay on disk;
only aggregate counts are printed. No merging, staging, or uploads.
"""

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime, timezone
import gzip
from itertools import combinations, groupby
import json
from pathlib import Path
import sys

import compare_stringent_recollapse as workflow
import pysam
import truvari
import truvari.collapse as native

VERSION = 1
STATUSES = ("eligible_pair", "similarity_only", "shared_carriers_only", "both",
            "native_filter", "other_matching_rule")
FAILURES = ("fails_size", "fails_sequence", "fails_overlap")
PAIR_FIELDS = ("threshold", "svtype", "mode", "chrom", "pos", "id_a", "id_b", "size_a", "size_b",
               "size_similarity", "sequence_similarity", "reciprocal_overlap", "carriers_a", "carriers_b",
               "shared_carriers", *FAILURES, "native_filtered_a", "native_filtered_b", "status")


def params_from_receipt(receipt):
    command = receipt["command"]
    args = native.parse_args(command[command.index("-i"):])
    params = truvari.VariantParams(args=args, short_circuit=False, skip_gt=True,
                                  sizefilt=args.sizemin, chunksize=args.refdist)
    if args.gt != "all" or args.dynthresh is not None:
        raise ValueError("Audit requires gt=all and fixed similarity thresholds")
    return params


def audit_locus(task):
    records, configs, kind, mode = task
    params = configs[0.7]
    entries = [truvari.VariantRecord(record, params=params) for record in records]
    if kind == "INS" and any(not e.is_resolved() for e in entries):
        raise ValueError("Sequence-resolved INS partition contains unresolved records")
    masks = [native.CollapsedCalls.make_genotype_mask(e, "all") for e in entries]
    filtered = [e.filter_call(base=True) or e.filter_size(base=True) for e in entries]
    rows = []
    for i, j in combinations(range(len(entries)), 2):
        a, b = entries[i], entries[j]
        # short_circuit=False computes every similarity even if size already
        # fails, so reasons are independent rather than first-failure counts.
        match = a.var_match(b)
        shared = int((masks[i] & masks[j]).sum())
        sa, ea = a.boundaries()
        sb, eb = b.boundaries()
        for threshold, config in sorted(configs.items(), reverse=True):
            other_pass = ((a.chrom, a.pos) == (b.chrom, b.pos) and a.same_type(b)
                          and truvari.overlaps(sa-config.refdist, ea+config.refdist, sb, eb))
            fails_size = match.sizesim < config.pctsize
            fails_sequence = config.pctseq > 0 and match.seqsim < config.pctseq
            fails_overlap = match.ovlpct < config.pctovl
            sim_pass = not (fails_size or fails_sequence or fails_overlap)
            if threshold == 0.7 and bool(match.state) != bool(sim_pass and other_pass):
                raise ValueError("Audit disagrees with native Truvari pair matching")
            if filtered[i] or filtered[j]:
                status = "native_filter"
            elif not other_pass:
                status = "other_matching_rule"
            elif not sim_pass and shared:
                status = "both"
            elif not sim_pass:
                status = "similarity_only"
            elif shared:
                status = "shared_carriers_only"
            else:
                status = "eligible_pair"
            rows.append(dict(threshold=threshold, svtype=kind, mode=mode, chrom=a.chrom, pos=a.pos,
                id_a=a.id, id_b=b.id, size_a=a.var_size(), size_b=b.var_size(),
                size_similarity=match.sizesim, sequence_similarity=match.seqsim if config.pctseq else "NA",
                reciprocal_overlap=match.ovlpct, carriers_a=int(masks[i].sum()),
                carriers_b=int(masks[j].sum()), shared_carriers=shared,
                fails_size=int(fails_size), fails_sequence=int(fails_sequence), fails_overlap=int(fails_overlap),
                native_filtered_a=int(filtered[i]), native_filtered_b=int(filtered[j]), status=status))
    return rows


def show(rows):
    fields = ("threshold", "svtype", "mode", "records", "repeated_positions", "same_position_pairs", *STATUSES)
    print("\t".join(fields))
    for row in rows:
        print("\t".join(str(row[k]) for k in fields))
    if any(row["eligible_pair"] for row in rows):
        print("Eligible pairs exist: similarity/shared-carrier checks alone do not explain zero merges. "
              "Inspect native grouping and cluster behavior next; pair counts are not predicted merge counts.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True, help="Completed lenient same-start run directory")
    parser.add_argument("--out", type=Path, help="Default: RUN/pair_audit")
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    if truvari.__version__ != "5.5.0" or args.workers < 1:
        parser.error("Requires Truvari 5.5.0 and positive workers")
    root = args.run.resolve()
    out = args.out.resolve() if args.out else root / "pair_audit"
    if out == root or out == root / "partitions":
        parser.error("Audit output must be a separate directory")
    out.mkdir(parents=True, exist_ok=True)
    manifest_path = root / "partitions/manifest.json"
    prepared = json.loads(manifest_path.read_text())
    sample_hash = workflow.totals_hash(prepared["samples"])
    jobs, excluded, inputs = [], [], {}
    print("Phase 1: validate prepared inputs and completed run parameters.", flush=True)
    for slug, meta in sorted(prepared["partitions"].items()):
        kind, mode = meta["svtype"].split(":")[0], meta["mode"]
        eligible = (kind == "DEL" and mode in ("sequence_resolved", "symbolic_size_position")) or (
                    kind == "INS" and mode == "sequence_resolved")
        if not eligible:
            excluded.append(dict(partition=slug, records=meta["records"], reason="passthrough_policy"))
            continue
        output_mode = ("symbolic_size_overlap" if mode == "symbolic_size_position" else "literal_size_overlap") if kind == "DEL" else mode
        output_slug = slug.removesuffix(mode) + output_mode
        source = root / "partitions" / meta["filename"]
        source_hash = workflow.digest(source)
        if source_hash != prepared["hashes"][source.name]:
            raise ValueError(f"Corrupt or changed partition: {source}")
        with pysam.VariantFile(source) as vcf:
            if workflow.totals_hash(list(vcf.header.samples)) != sample_hash:
                raise ValueError("Sample manifest/order mismatch")
        receipts = {}
        for threshold in (0.9, 0.7):
            folder = root / f"pct{round(threshold*100)}" / output_slug
            receipt = json.loads((folder / "complete.json").read_text())
            recipe = receipt["recipe"]
            if (not recipe.get("same_start") or recipe["truvari"] != "5.5.0"
                    or recipe["input_sha256"] != source_hash or recipe["samples_sha256"] != sample_hash):
                raise ValueError("Completed run does not match prepared same-position input")
            config = params_from_receipt(receipt)
            if (config.pctsize != threshold or config.pctseq != (threshold if kind == "INS" else 0)
                    or config.pctovl != (threshold if kind == "DEL" else 0)
                    or config.sizemin != 10001 or config.sizemax != -1 or config.typeignore):
                raise ValueError("Unexpected run criteria; audit must be updated to match them")
            receipts[threshold] = receipt
        inputs[slug] = dict(sha256=source_hash, recipes={str(t): r["recipe"] for t, r in receipts.items()},
                           native_parameters={str(t): vars(params_from_receipt(r)) for t, r in receipts.items()})
        jobs.append((source, meta, output_mode, receipts))
    key = dict(version=VERSION, script_sha256=workflow.digest(__file__), truvari=truvari.__version__,
               samples_sha256=sample_hash, inputs=inputs, manifest_sha256=workflow.digest(manifest_path))
    receipt_path = out / "complete.json"
    if receipt_path.exists():
        cached = json.loads(receipt_path.read_text())
        if cached["key"] == key and workflow.files_valid(out, cached["hashes"]):
            print("Reusing verified pair audit", flush=True)
            show(cached["summary"])
            return
    summary = []
    detail = out / "pairs.tsv.gz"
    with gzip.open(detail, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=PAIR_FIELDS, delimiter="\t")
        writer.writeheader()
        for source, meta, mode, receipts in jobs:
            print(f"Phase 2: count shared positions in {meta['svtype']}/{mode}.", flush=True)
            with pysam.VariantFile(source, drop_samples=True) as vcf:
                counts = Counter((r.chrom, r.pos) for r in vcf)
            if sum(counts.values()) != meta["records"]:
                raise ValueError("Partition record count disagrees with manifest")
            repeated = {k for k, n in counts.items() if n > 1}
            pairs = sum(n*(n-1)//2 for n in counts.values())
            print(f"{meta['records']:,} records; {len(repeated):,} repeated positions; {pairs:,} pairs.", flush=True)
            configs = {t: params_from_receipt(r) for t, r in receipts.items()}
            tallies = {t: Counter() for t in configs}
            def tasks(vcf):
                for locus, records in groupby(vcf, key=lambda r: (r.chrom, r.pos)):
                    if locus in repeated:
                        group = list(records)
                        if len(group) != counts[locus]:
                            raise ValueError("Same-position records are not contiguous/sorted")
                        yield group, configs, meta["svtype"].split(":")[0], mode
            print(f"Phase 3: score all pairs and carrier intersections with {args.workers} workers.", flush=True)
            processed = 0
            with pysam.VariantFile(source) as vcf, ThreadPoolExecutor(max_workers=args.workers) as pool:
                for rows in pool.map(audit_locus, tasks(vcf)):
                    for row in rows:
                        writer.writerow(row)
                        tally = tallies[row["threshold"]]
                        tally[row["status"]] += 1
                        for field in FAILURES:
                            tally[field] += row[field]
                    processed += len(rows)//2
                    print(f"Compared {processed:,}/{pairs:,} pairs in {meta['svtype']}/{mode}", flush=True)
            for t, tally in sorted(tallies.items(), reverse=True):
                if sum(tally[s] for s in STATUSES) != pairs:
                    raise ValueError("Pair accounting failed")
                summary.append(dict(threshold=t, svtype=meta["svtype"], mode=mode, records=meta["records"],
                    unique_positions=len(counts), repeated_positions=len(repeated),
                    records_at_repeated_positions=sum(counts[k] for k in repeated), same_position_pairs=pairs,
                    **{s: tally[s] for s in (*STATUSES, *FAILURES)},
                    observed_collapsed_records=receipts[t]["summary"]["collapsed_records"]))
            if workflow.digest(source) != inputs[source.stem.removesuffix(".vcf")]["sha256"]:
                raise ValueError("Input changed during audit")
    if not summary:
        raise ValueError("No eligible DEL/INS partitions")
    table = out / "summary.tsv"
    with table.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary)
    provenance = out / "provenance.json"
    workflow.write_json(provenance, dict(utc=datetime.now(timezone.utc).isoformat(), run=str(root),
        workers=args.workers, command=sys.argv, key=key, excluded=excluded, random_seed=None,
        randomness="No random operations", shared_carrier="Native gt=all mask: ALT allele 1 in both GTs",
        interpretation="Pairwise eligibility is not a predicted cluster/merge count; failure columns overlap",
        privacy="No sample names or allele sequences written to pair table; locus/ID details remain local"))
    workflow.write_json(receipt_path, dict(key=key, summary=summary,
        hashes={p.name: workflow.digest(p) for p in (detail, table, provenance)}))
    show(summary)
    print(f"Saved aggregate summary and local pair details to {out}", flush=True)


if __name__ == "__main__":
    main()
