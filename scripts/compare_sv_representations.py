#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = ["truvari==5.5.0", "google-crc32c==1.8.0"]
# ///
"""Compare chr1 bubble, popped and imputed allele representations on the VM."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from collections import Counter
import hashlib
import itertools
import json
import logging
from pathlib import Path
import re
import sqlite3
import sys

import edlib
import pysam
import truvari
import inspect_sv_representation as representation
import run_chr1_sv_merge as pilot

VERSION = 1
SOURCES = {
    "bubble": pilot.SOURCE,
    "popped": representation.SOURCE,
    "imputed": "gs://long-reads-phase-2-imputation-vwb/batch-001/aou2_50k_init.batch-001.chr1.prod-05fea60.glimpse2.popped.bcf",
}


def digest(*parts):
    return hashlib.sha256(json.dumps(parts, separators=(",", ":")).encode()).hexdigest()


def trim(pos, ref, alt):
    ref, alt = ref.upper(), alt.upper()
    info, _ = representation.describe(ref, alt)
    prefix, suffix = info["shared_prefix_bp"], info["shared_suffix_bp"]
    return pos + prefix, ref[prefix:len(ref) - suffix if suffix else None], alt[prefix:len(alt) - suffix if suffix else None], info


def scan_rows(local, path, regions=None):
    cmd = [local.bcftools, "query", "-u", "-r", regions or "chr1", "-f",
           "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/SVTYPE\t%INFO/SVLEN\n", path]
    with pilot.stream(local, cmd, "sites" if regions is None else "examples") as lines:
        for line in lines:
            chrom, pos, ident, ref, alts, svtype, svlen = line.rstrip("\n").split("\t")
            if chrom != "chr1":
                raise ValueError("Unexpected chromosome")
            yield int(pos), ident, ref, alts.split(","), svtype, svlen


def index_input(args, label, path):
    root = args.root / label
    (root / "logs").mkdir(parents=True, exist_ok=True)
    local = argparse.Namespace(root=root, bcftools=args.bcftools)
    indexes = [Path(str(path) + ext) for ext in (".csi", ".tbi") if Path(str(path) + ext).is_file()]
    if not indexes:
        raise ValueError(f"Missing index for {path}")
    identities = [pilot.common.identity(p) for p in [path, *indexes]]
    key = {"schema": VERSION, "inputs": identities, "min_net_bp": args.min_net_bp,
           "bcftools": pilot.command(local, [args.bcftools, "--version"], "version", True).splitlines()[0]}
    if pilot.valid_manifest(root, key):
        print(f"Reuse verified {label} allele index.", flush=True)
        return root / "alleles.sqlite"
    (root / "done.json").unlink(missing_ok=True)
    partial = root / "alleles.partial.sqlite"
    partial.unlink(missing_ok=True)
    db = sqlite3.connect(partial)
    db.execute("PRAGMA journal_mode=DELETE")
    db.execute("PRAGMA cache_size=-65536")
    db.executescript("""
        CREATE TABLE exact(k TEXT PRIMARY KEY, pos INTEGER, kind TEXT, complex INTEGER,
                           refhash TEXT, canon TEXT, n INTEGER);
        CREATE TABLE canonical(k TEXT PRIMARY KEY, kind TEXT, n INTEGER);
    """)
    counts = Counter()
    previous = 0
    try:
        for pos, ident, ref, alts, svtype, svlen in scan_rows(local, path):
            if pos < previous:
                raise ValueError("Input must be coordinate sorted")
            previous = pos
            counts["input_records"] += 1
            counts["multiallelic_records"] += len(alts) > 1
            for alt in alts:
                counts["input_alt_entries"] += 1
                info, excluded = representation.describe(ref, alt)
                if excluded or info["type"] not in {"INS", "DEL"} or info["net_bp"] < args.min_net_bp:
                    continue
                newpos, r, a, _ = trim(pos, ref, alt)
                k, c = digest("chr1", pos, ref, alt), digest("chr1", newpos, r, a)
                db.execute("INSERT INTO exact VALUES(?,?,?,?,?,?,1) ON CONFLICT(k) DO UPDATE SET n=n+1",
                           (k, pos, info["type"], int(info["complex"]), digest(ref), c))
                db.execute("INSERT INTO canonical VALUES(?,?,1) ON CONFLICT(k) DO UPDATE SET n=n+1", (c, info["type"]))
                counts["selected_alleles"] += 1
            if counts["input_records"] % 100000 == 0:
                db.commit()
                logging.info("PROGRESS %s %s", label, dict(counts))
                print(f"{label}: scanned {counts['input_records']:,} records.", flush=True)
        db.execute("CREATE INDEX locus ON exact(pos,refhash,kind)")
        db.commit()
        if db.execute("PRAGMA integrity_check").fetchone()[0] != "ok":
            raise ValueError("SQLite integrity check failed")
    finally:
        db.close()
    if identities != [pilot.common.identity(p) for p in [path, *indexes]]:
        raise ValueError("Input changed during scan")
    header = pilot.command(local, [args.bcftools, "view", "-h", path], "header", True)
    columns = next((line.split("\t") for line in header.splitlines() if line.startswith("#CHROM")), [])
    counts["sample_columns"] = max(0, len(columns) - 9)
    (root / "header.txt").write_text(header)
    pilot.common.json_write(root / "counts.json", {"counts": dict(counts), "input": str(path)})
    target = root / "alleles.sqlite"
    partial.replace(target)
    pilot.manifest(root, key, [target, root / "counts.json", root / "header.txt"])
    return target


def compare_databases(paths, examples):
    rows, requests, pair_requests = [], {label: {} for label in paths}, []
    labels = list(paths)
    db = sqlite3.connect(":memory:")
    for i, label in enumerate(labels):
        db.execute(f"ATTACH DATABASE ? AS d{i}", (str(paths[label]),))
    for i, j in itertools.combinations(range(len(labels)), 2):
        left, right = labels[i], labels[j]
        for kind in ("INS", "DEL"):
            for table in ("exact", "canonical"):
                a, b = f"d{i}.{table}", f"d{j}.{table}"
                na, ra = db.execute(f"SELECT count(*),coalesce(sum(n),0) FROM {a} WHERE kind=?", (kind,)).fetchone()
                nb, rb = db.execute(f"SELECT count(*),coalesce(sum(n),0) FROM {b} WHERE kind=?", (kind,)).fetchone()
                shared, matched = db.execute(f"SELECT count(*),coalesce(sum(min(a.n,b.n)),0) FROM {a} a JOIN {b} b USING(k) WHERE a.kind=?", (kind,)).fetchone()
                rows.append({"left": left, "right": right, "type": kind, "comparison": table,
                    "left_unique": na, "right_unique": nb, "shared_unique": shared,
                    "left_only_unique": na - shared, "right_only_unique": nb - shared,
                    "left_records": ra, "right_records": rb, "shared_record_multiplicity": matched,
                    "left_shared_percent": 100 * shared / na if na else None,
                    "right_shared_percent": 100 * shared / nb if nb else None,
                    "sets_equal": na == nb == shared, "multisets_equal": ra == rb == matched})
            a, b = f"d{i}.exact", f"d{j}.exact"
            for category, join, condition in (("shared", "JOIN", ""), ("left_only", "LEFT JOIN", "AND b.k IS NULL")):
                for complex_ in (0, 1):
                    selected = db.execute(f"SELECT a.k,a.pos FROM {a} a {join} {b} b USING(k) WHERE a.kind=? AND a.complex=? {condition} ORDER BY a.k LIMIT ?", (kind, complex_, examples))
                    for k, pos in selected:
                        requests[left].setdefault(k, {"pos": pos, "categories": []})["categories"].append(f"{left}_vs_{right}:{category}:{kind}:{complex_}")
                        if category == "shared":
                            requests[right].setdefault(k, {"pos": pos, "categories": []})["categories"].append(f"shared_with_{left}")
            selected = db.execute(f"SELECT b.k,b.pos FROM {b} b LEFT JOIN {a} a USING(k) WHERE b.kind=? AND a.k IS NULL ORDER BY b.k LIMIT ?", (kind, examples))
            for k, pos in selected:
                requests[right].setdefault(k, {"pos": pos, "categories": []})["categories"].append(f"only_vs_{left}:{kind}")
    # Diagnostic pairs within a POS/REF group: two hash-ordered distinct ALTs.
    for i, label in enumerate(labels):
        for kind in ("INS", "DEL"):
            groups = db.execute(f"SELECT pos,refhash FROM d{i}.exact WHERE kind=? AND complex=1 GROUP BY pos,refhash HAVING count(*)>1 ORDER BY refhash,pos LIMIT ?", (kind, examples)).fetchall()
            for pos, refhash in groups:
                keys = [r[0] for r in db.execute(f"SELECT k FROM d{i}.exact WHERE kind=? AND pos=? AND refhash=? ORDER BY k LIMIT 2", (kind, pos, refhash))]
                pair_requests.append((label, *keys))
                for k in keys:
                    requests[label].setdefault(k, {"pos": pos, "categories": []})["categories"].append(f"within_site_pair:{kind}")
    db.close()
    return rows, requests, pair_requests


def alignment(ref, alt, maximum):
    if max(len(ref), len(alt)) > maximum:
        return {"skipped": f"sequence exceeds {maximum} bp alignment limit"}
    result = edlib.align(alt.upper(), ref.upper(), mode="NW", task="path")
    if not ref or not alt:
        result["cigar"] = f"{len(ref)}D" if ref else f"{len(alt)}I" if alt else ""
    operations = Counter()
    for count, op in re.findall(r"(\d+)([=XID])", result["cigar"] or ""):
        operations[op] += int(count)
    return {"edit_distance": result["editDistance"], "cigar_ALT_query_vs_REF_target": result["cigar"],
            "operation_bp": dict(operations), "caution": "One optimal alignment; repeats can permit alternative edit placements"}


def fetch_examples(args, label, path, requested):
    result = {}
    if not requested:
        return result
    local = argparse.Namespace(root=args.root / label, bcftools=args.bcftools)
    regions = ",".join(f"chr1:{p}-{p}" for p in sorted({r["pos"] for r in requested.values()}))
    for pos, ident, ref, alts, svtype, svlen in scan_rows(local, path, regions):
        for index, alt in enumerate(alts):
            k = digest("chr1", pos, ref, alt)
            if k not in requested or k in result:
                continue
            newpos, r, a, info = trim(pos, ref, alt)
            def annotation(value):
                values = value.split(",")
                return values[index] if len(values) == len(alts) else values[0] if len(values) == 1 else "."
            result[k] = {"key": k, "dataset": label, "chrom": "chr1", "pos": pos, "id": ident,
                "ref": ref, "alt": alt, "original_alt_index": index + 1, "original_alt_count": len(alts),
                "svtype": annotation(svtype), "svlen": annotation(svlen),
                "categories": requested[k]["categories"], "metrics": info,
                "trimmed": {"pos": newpos, "ref": r, "alt": a},
                "alignment": alignment(r, a, args.max_align_bp)}
    if set(result) != set(requested):
        raise ValueError(f"Could not recover selected examples in {label}; input/index may have changed")
    return result


def native_record(row):
    header = pysam.VariantHeader()
    header.contigs.add("chr1", length=1000000000)
    header.info.add("SVTYPE", 1, "String", "Source annotation")
    header.info.add("SVLEN", 1, "Integer", "Source annotation")
    record = header.new_record(contig="chr1", start=row["pos"] - 1, alleles=(row["ref"], row["alt"]))
    if row["svtype"] != ".":
        record.info["SVTYPE"] = row["svtype"]
    if row["svlen"] != ".":
        record.info["SVLEN"] = int(row["svlen"])
    return truvari.VariantRecord(record)


def pair_scores(a, b, maximum):
    if max(len(r[s]) for r in (a, b) for s in ("ref", "alt")) > maximum:
        return {"skipped": f"full allele exceeds {maximum} bp score limit"}
    first, second = native_record(a), native_record(b)
    seq, size = first.seqsim(second), first.sizesim(second)[0]
    result = {"native_type_a": first.var_type().name, "native_type_b": second.var_type().name,
              "native_seqsim": seq, "native_sizesim": size,
              "passes_p70_sequence_size_only": seq >= .7 and size >= .7,
              "passes_p90_sequence_size_only": seq >= .9 and size >= .9}
    for field in ("ref", "alt"):
        x, y = a["trimmed"][field], b["trimmed"][field]
        distance = edlib.align(x, y)["editDistance"]
        result[f"trimmed_{field}_direct_score"] = 1 - distance / (len(x) + len(y)) if x or y else 1.0
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--phase", choices=("stage", "compare", "all"), default="all")
    for label in SOURCES:
        parser.add_argument("--" + label, type=Path, help=f"Existing indexed {label} file; otherwise stage the original BCF")
    parser.add_argument("--bcftools", default="bcftools")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--memory-gb", type=int, default=26)
    parser.add_argument("--min-net-bp", type=int, default=50)
    parser.add_argument("--examples", type=int, default=3, help="Examples per category, deterministic SHA256 order")
    parser.add_argument("--max-align-bp", type=int, default=10000)
    args = parser.parse_args()
    if min(args.threads, args.examples, args.max_align_bp, args.min_net_bp) < 1:
        parser.error("Numeric settings must be positive")
    args.root = args.root.resolve()
    args.root.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(filename=args.root / "comparison.log", level=logging.INFO, force=True,
                        format="%(asctime)s %(levelname)s %(message)s")
    logging.info("PARAMETERS %s", json.dumps(vars(args), default=str, sort_keys=True))
    logging.info("PROVENANCE python=%s truvari=%s seed=%s script=%s", sys.version, truvari.__version__, pilot.SEED, pilot.common.digest(__file__))
    import fcntl
    with (args.root / ".comparison.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        inputs = {}
        for label, uri in SOURCES.items():
            supplied = getattr(args, label)
            if supplied:
                inputs[label] = supplied.resolve()
                continue
            staging = args.root / "staging" / label
            (staging / "logs").mkdir(parents=True, exist_ok=True)
            if args.phase in {"all", "stage"}:
                pa = pilot.parse_args(["--root", str(staging), "--source-uri", uri, "--confirm-source-uri", uri,
                    "--threads", str(args.threads), "--memory-gb", str(args.memory_gb), "--scratch-reserve-gb", "4", "--bcftools", args.bcftools])
                print(f"Stage/verify {label}; logs stay on VM.", flush=True)
                pilot.stage(pa)
            else:
                pilot.require_manifest(staging / "inputs")
            inputs[label] = staging / "inputs" / uri.rsplit("/", 1)[1]
        if args.phase == "stage":
            return
        with ThreadPoolExecutor(max_workers=min(args.threads, 3)) as pool:
            databases = dict(zip(inputs, pool.map(lambda item: index_input(args, *item), inputs.items())))
        comparison_key = {"schema": VERSION, "datasets": {label: json.loads((args.root / label / "done.json").read_text())["key"] for label in inputs},
                          "examples": args.examples, "max_align_bp": args.max_align_bp, "truvari": truvari.__version__}
        if pilot.valid_manifest(args.root, comparison_key):
            print("Reuse verified three-file comparison.", flush=True)
            print((args.root / "record_structure.tsv").read_text(), end="")
            print((args.root / "pairwise_summary.tsv").read_text(), end="")
            print(f"Individual examples stay on VM: {args.root}")
            return
        (args.root / "done.json").unlink(missing_ok=True)
        rows, requested, pairs = compare_databases(databases, args.examples)
        examples = {label: fetch_examples(args, label, path, requested[label]) for label, path in inputs.items()}
        with (args.root / "allele_examples.jsonl").open("w") as out:
            for values in examples.values():
                for row in values.values():
                    out.write(json.dumps(row) + "\n")
        def preview(seq):
            return seq if len(seq) <= 160 else seq[:80] + f" ... [{len(seq) - 160} bases omitted] ... " + seq[-80:]
        with (args.root / "allele_examples.txt").open("w") as out:
            out.write("Individual alleles: keep this file on the VM. Full strings are in allele_examples.jsonl.\n"
                      "Trimming removes common flanks only. Alignment is one possible minimum-edit explanation.\n\n")
            for values in examples.values():
                for row in values.values():
                    out.write(f"{row['dataset']} {row['chrom']}:{row['pos']} key={row['key']}\n")
                    out.write(f"Original ID={row['id']} ALT index={row['original_alt_index']}/{row['original_alt_count']}\n")
                    out.write("Categories: " + ", ".join(row["categories"]) + "\n")
                    out.write("Metrics: " + json.dumps(row["metrics"]) + "\n")
                    for field in ("ref", "alt"):
                        out.write(f"Original {field.upper()}: {preview(row[field])}\n")
                    out.write(f"After trimming, start={row['trimmed']['pos']} (empty REF/ALT permitted in this diagnostic):\n")
                    for field in ("ref", "alt"):
                        out.write(f"  {field.upper()}: {preview(row['trimmed'][field]) or '(empty)'}\n")
                    out.write("Alignment: " + json.dumps(row["alignment"]) + "\n\n")
        with (args.root / "within_site_pairs.jsonl").open("w") as out:
            for label, ka, kb in pairs:
                scores = pair_scores(examples[label][ka], examples[label][kb], args.max_align_bp)
                out.write(json.dumps({"dataset": label, "key_a": ka, "key_b": kb, **scores,
                    "caution": "Scores only: no genotypes, GT conflict test, or actual collapse; trimmed direct scores are diagnostic, not native Truvari after normalization"}) + "\n")
        pilot.tsv(args.root / "pairwise_summary.tsv", rows, list(rows[0]))
        structures = []
        for label in inputs:
            counts = json.loads((args.root / label / "counts.json").read_text())["counts"]
            structures.append({"dataset": label, **{field: counts.get(field, 0) for field in
                ("input_records", "multiallelic_records", "input_alt_entries", "selected_alleles", "sample_columns")}})
        pilot.tsv(args.root / "record_structure.tsv", structures, list(structures[0]))
        pilot.common.json_write(args.root / "comparison_metadata.json", {
            "inputs": {k: str(v) for k, v in inputs.items()}, "min_net_bp": args.min_net_bp,
            "exact": "SHA256 of CHROM/POS/REF/ALT; ALT entries compared independently of row packing, ID, INFO, GT, AC and AF",
            "canonical": "Uppercase sequences; remove common prefix/suffix and shift POS, allowing empty strings; no repeat left alignment, atomization or general haplotype equivalence",
            "examples": "Deterministic SHA256 selection per category, illustrative rather than prevalence estimates",
            "seed": pilot.SEED, "truvari": truvari.__version__})
        for label in inputs:
            for expected in comparison_key["datasets"][label]["inputs"]:
                if expected != pilot.common.identity(expected["path"]):
                    raise ValueError("An input or index changed during example extraction")
        pilot.manifest(args.root, comparison_key, [args.root / name for name in
            ("pairwise_summary.tsv", "record_structure.tsv", "allele_examples.jsonl", "allele_examples.txt",
             "within_site_pairs.jsonl", "comparison_metadata.json")])
        print((args.root / "record_structure.tsv").read_text(), end="")
        print((args.root / "pairwise_summary.tsv").read_text(), end="")
        print(f"Individual sequences, edit alignments and pair scores saved ONLY on VM: {args.root}")


if __name__ == "__main__":
    main()
