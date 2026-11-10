# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
"""Compare four v3 callsets locally, skipping BND alleles; no VCF modification.

Streams annotations through bcftools, stores exact allele keys in SQLite, and
keeps sequence-resolved matches separate from symbolic metadata matches.
No sample genotypes or sequences are printed or uploaded.
"""

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime, timezone
import hashlib
from itertools import combinations
import json
from pathlib import Path
import re
import sqlite3
import subprocess

LABELS = ("main", "ultralong_all", "ultralong_lenient", "ultralong_stringent")
CACHE_SCHEMA = 1
BINS = ("lt_50bp", "50_999bp", "1kb_to_10kb", "gt_10kb_to_100kb", "gt_100kb", "unknown")


def digest_bytes(value):
    return hashlib.sha256(value).hexdigest()


def stat_signature(path):
    stat = path.stat()
    return {"path": str(path.resolve()), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def write_tsv(path, rows, columns):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def query(command, errors):
    with errors.open("w", encoding="utf-8") as stderr:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=stderr,
                                   text=True, encoding="utf-8")
        try:
            yield from process.stdout
            if process.wait():
                raise RuntimeError(f"bcftools failed; see {errors}")
        finally:
            process.stdout.close()
            if process.poll() is None:
                process.terminate()
            process.wait()


def allele_value(value, index, count, number):
    if value in (".", ""):
        return None
    values = value.split(",")
    if number == "A":
        return values[index] if len(values) == count and values[index] != "." else None
    if number == "R":
        return values[index + 1] if len(values) == count + 1 and values[index + 1] != "." else None
    if len(values) == 1:
        return values[0]
    return values[index] if len(values) == count and values[index] != "." else None


def length_bin(length):
    if length is None:
        return "unknown"
    if length < 50:
        return "lt_50bp"
    if length < 1000:
        return "50_999bp"
    if length <= 10000:
        return "1kb_to_10kb"
    if length <= 100000:
        return "gt_10kb_to_100kb"
    return "gt_100kb"


def describe_allele(chrom, pos, ref, alt, svtype, svlen, end):
    ref, alt = ref.upper(), alt.upper()
    if alt in (".", "*", "<NON_REF>", "<*>"):
        return {"skip": "placeholder"}
    symbolic = alt.startswith("<") and alt.endswith(">")
    kind = (svtype or (alt[1:-1] if symbolic else "")).split(":")[0].upper()
    symbolic_kind = alt[1:-1].split(":")[0] if symbolic else ""
    if kind in ("BND", "TRA", "CTX") or symbolic_kind in ("BND", "TRA", "CTX") or "[" in alt or "]" in alt or alt.startswith(".") or alt.endswith("."):
        return {"skip": "BND"}
    literal = bool(re.fullmatch("[ACGTN]+", ref) and re.fullmatch("[ACGTN]+", alt))
    if not kind:
        if literal:
            kind = "INS" if len(alt) > len(ref) else "DEL" if len(ref) > len(alt) else "SNV" if len(ref) == 1 else "MNV"
        else:
            kind = "UNKNOWN"
    length = abs(int(svlen)) if svlen is not None else None
    stop = int(end) if end is not None else None
    if not length:
        length = None
        if literal and len(ref) != len(alt):
            length = abs(len(alt) - len(ref))
        elif literal and kind in ("SNV", "MNV"):
            length = len(ref)
        elif kind in ("DEL", "DUP", "INV", "CPX") and stop is not None and stop > pos:
            length = stop - pos
    if literal:
        # Same literal alleles match despite differences in INFO annotation.
        evidence = "sequence_resolved" if "N" not in ref + alt else "literal_with_N"
        identity = [chrom, pos, ref, alt, kind]
    elif symbolic and length is not None:
        evidence = "symbolic_metadata"
        # END is uninformative for INS. Derive a missing span end for DEL/DUP/INV.
        effective_end = None if kind == "INS" else stop
        if effective_end is None and kind in ("DEL", "DUP", "INV"):
            effective_end = pos + length
        identity = [chrom, pos, ref, alt, kind, effective_end, length]
    else:
        evidence = "unmatchable"
        identity = None
    key = hashlib.sha256(json.dumps(identity, separators=(",", ":")).encode()).digest() if identity else None
    return dict(skip=None, kind=kind, length=length, length_bin=length_bin(length), evidence=evidence, key=key)


def parse_line(line, numbers):
    fields = line.rstrip("\n").split("\t")
    if len(fields) != 8:
        raise ValueError("Expected eight bcftools annotation fields")
    chrom, pos, ref, alts, types, lengths, ends, filters = fields
    pos = int(pos)
    if not chrom or pos < 1:
        raise ValueError("Invalid chromosome or position")
    alt_list = alts.split(",")
    alleles = []
    for index, alt in enumerate(alt_list):
        values = [allele_value(raw, index, len(alt_list), numbers.get(tag))
                  for raw, tag in ((types, "SVTYPE"), (lengths, "SVLEN"), (ends, "END"))]
        if len(alt_list) > 1 and "," not in lengths:
            # A scalar SVLEN on a multiallelic row cannot describe each allele reliably.
            values[1] = None
        alleles.append(describe_allele(chrom, pos, ref, alt, *values))
    return alleles, len(alt_list) > 1, filters


def valid_cache(db, metadata, fingerprint):
    if not db.exists() or not metadata.exists():
        return None
    try:
        saved = json.loads(metadata.read_text(encoding="utf-8"))
        if saved["fingerprint"] != fingerprint:
            return None
        with sqlite3.connect(f"{db.resolve().as_uri()}?mode=ro", uri=True) as conn:
            if conn.execute("PRAGMA quick_check").fetchone()[0] != "ok":
                return None
            if conn.execute("SELECT count(*) FROM alleles").fetchone()[0] != saved["summary"]["unique_matchable_keys"]:
                return None
        return saved
    except (OSError, ValueError, KeyError, sqlite3.Error):
        return None


def scan_callset(task):
    label, source, out, version = task
    cache = out / "cache" / label
    cache.mkdir(parents=True, exist_ok=True)
    before = stat_signature(source)
    header_cmd = ["bcftools", "view", "--no-version", "-h", str(source)]
    header = "".join(query(header_cmd, cache / "header.stderr.log"))
    columns = [line for line in header.splitlines() if line.startswith("#CHROM\t")]
    if len(columns) != 1:
        raise ValueError(f"{label}: missing VCF column header")
    samples = columns[0].split("\t")[9:]
    if len(set(samples)) != len(samples):
        raise ValueError(f"{label}: duplicate sample names")
    sample_text = "".join(s + "\n" for s in samples)
    fingerprint = dict(input=before, schema=CACHE_SCHEMA, bcftools_version=version,
                       header_sha256=digest_bytes(header.encode()), sample_manifest_sha256=digest_bytes(sample_text.encode()))
    db, metadata = cache / "alleles.sqlite", cache / "scan.json"
    saved = valid_cache(db, metadata, fingerprint)
    if saved:
        print(f"{label}: reused validated annotation cache", flush=True)
        return label, db, saved, samples
    print(f"{label}: streaming annotations (all FILTER values; BND skipped)", flush=True)
    (cache / "header.vcf").write_text(header, encoding="utf-8")
    (cache / "samples.txt").write_text(sample_text, encoding="utf-8")
    numbers = dict(re.findall(r"^##INFO=<ID=([^,]+),Number=([^,]+),", header, re.MULTILINE))
    tmp = cache / "alleles.build.sqlite"
    if tmp.exists():
        tmp.unlink()
    stats, types, bins, filters, evidence = Counter(), Counter(), Counter(), Counter(), Counter()
    conn = sqlite3.connect(tmp)
    conn.execute("PRAGMA cache_size=-65536")
    conn.execute("CREATE TABLE alleles (key BLOB PRIMARY KEY, svtype TEXT NOT NULL, evidence TEXT NOT NULL) WITHOUT ROWID")
    command = ["bcftools", "query", "-u", "-f", "%CHROM\t%POS\t%REF\t%ALT\t%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/END\t%FILTER\n", str(source)]
    batch = []
    try:
        for line in query(command, cache / "query.stderr.log"):
            alleles, multi, filt = parse_line(line, numbers)
            stats["input_records"] += 1
            stats["multiallelic_input_records"] += int(multi)
            kept = [a for a in alleles if not a["skip"]]
            stats["records_with_BND"] += int(any(a["skip"] == "BND" for a in alleles))
            stats["BND_alleles_skipped"] += sum(a["skip"] == "BND" for a in alleles)
            stats["placeholder_alleles_skipped"] += sum(a["skip"] == "placeholder" for a in alleles)
            if kept:
                stats["retained_records"] += 1
                stats["retained_multiallelic_records"] += int(len(kept) > 1)
                filters[filt] += 1
            for allele in kept:
                stats["retained_ALT_alleles"] += 1
                types[allele["kind"]] += 1
                bins[allele["kind"], allele["length_bin"]] += 1
                evidence[allele["evidence"]] += 1
                if allele["length"] is not None and allele["length"] > 10000:
                    stats["ALT_alleles_gt10kb"] += 1
                if allele["key"] is not None:
                    batch.append((allele["key"], allele["kind"], allele["evidence"]))
            if len(batch) >= 10000:
                conn.executemany("INSERT OR IGNORE INTO alleles VALUES (?,?,?)", batch)
                conn.commit()
                batch.clear()
            if stats["input_records"] % 100000 == 0:
                print(f"{label}: decoded {stats['input_records']:,} records", flush=True)
        conn.executemany("INSERT OR IGNORE INTO alleles VALUES (?,?,?)", batch)
        conn.commit()
        stats["unique_matchable_keys"] = conn.execute("SELECT count(*) FROM alleles").fetchone()[0]
        stats["unmatchable_ALT_alleles"] = evidence["unmatchable"]
        stats["repeated_matchable_keys"] = stats["retained_ALT_alleles"] - evidence["unmatchable"] - stats["unique_matchable_keys"]
        if conn.execute("PRAGMA quick_check").fetchone()[0] != "ok":
            raise ValueError(f"{label}: SQLite integrity check failed")
    finally:
        conn.close()
    if before != stat_signature(source):
        raise ValueError(f"{label}: input changed during scan")
    tmp.replace(db)
    summary = dict(callset=label, samples=len(samples), **stats)
    for key in ("input_records", "multiallelic_input_records", "records_with_BND", "BND_alleles_skipped",
                "placeholder_alleles_skipped", "ALT_alleles_gt10kb", "retained_records", "retained_ALT_alleles", "retained_multiallelic_records"):
        summary.setdefault(key, 0)
    saved = dict(fingerprint=fingerprint, summary=summary, types=dict(types),
                 length_bins=[dict(svtype=t, length_bin=b, ALT_alleles=n) for (t, b), n in sorted(bins.items())],
                 filters=dict(filters), evidence=dict(evidence), commands=[header_cmd, command])
    metadata.write_text(json.dumps(saved, indent=2) + "\n", encoding="utf-8")
    print(f"{label}: complete; {stats['retained_ALT_alleles']:,} retained ALT alleles", flush=True)
    return label, db, saved, samples


def overlap_rows(memberships, summaries):
    rows = []
    for i, j in combinations(range(len(LABELS)), 2):
        a, b = LABELS[i], LABELS[j]
        shared = sum(n for mask, _, _, n in memberships if mask & (1 << i) and mask & (1 << j))
        na, nb = summaries[a]["unique_matchable_keys"], summaries[b]["unique_matchable_keys"]
        seq = sum(n for mask, _, evidence, n in memberships
                  if mask & (1 << i) and mask & (1 << j) and evidence == "sequence_resolved")
        symbolic = sum(n for mask, _, evidence, n in memberships
                       if mask & (1 << i) and mask & (1 << j) and evidence == "symbolic_metadata")
        rows.append(dict(callset_a=a, callset_b=b, unique_a=na, unique_b=nb, shared=shared,
                         only_a=na-shared, only_b=nb-shared, shared_sequence_resolved=seq,
                         shared_symbolic_metadata=symbolic, shared_literal_with_N=shared-seq-symbolic,
                         percent_a_shared=round(100*shared/na, 4) if na else None,
                         percent_b_shared=round(100*shared/nb, 4) if nb else None,
                         a_subset_of_b="yes" if shared == na and not summaries[a]["unmatchable_ALT_alleles"] else "no_or_unresolved",
                         b_subset_of_a="yes" if shared == nb and not summaries[b]["unmatchable_ALT_alleles"] else "no_or_unresolved"))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")
    args.input_dir, args.out = args.input_dir.resolve(), args.out.resolve()
    sources = {label: args.input_dir / f"v3_{label}.bcf" for label in LABELS}
    for source in sources.values():
        if not source.is_file():
            raise FileNotFoundError(source)
    args.out.mkdir(parents=True, exist_ok=True)
    version = subprocess.check_output(["bcftools", "--version"], text=True).splitlines()[0]
    tasks = [(label, sources[label], args.out, version) for label in LABELS]
    with ThreadPoolExecutor(max_workers=min(args.workers, len(LABELS))) as pool:
        scans = list(pool.map(scan_callset, tasks))
    summaries = {label: saved["summary"] for label, _, saved, _ in scans}
    membership_db = args.out / "membership.build.sqlite"
    if membership_db.exists():
        membership_db.unlink()
    conn = sqlite3.connect(membership_db)
    try:
        conn.execute("PRAGMA cache_size=-65536")
        conn.execute("CREATE TABLE membership (key BLOB PRIMARY KEY, svtype TEXT, evidence TEXT, mask INTEGER) WITHOUT ROWID")
        for index, (_, db, _, _) in enumerate(scans):
            conn.execute("ATTACH DATABASE ? AS incoming", (str(db),))
            conn.execute("INSERT INTO membership SELECT key,svtype,evidence,? FROM incoming.alleles WHERE 1 "
                         "ON CONFLICT(key) DO UPDATE SET mask=membership.mask | excluded.mask", (1 << index,))
            conn.commit()
            conn.execute("DETACH DATABASE incoming")
        memberships = list(conn.execute("SELECT mask,svtype,evidence,count(*) FROM membership GROUP BY mask,svtype,evidence ORDER BY mask,svtype,evidence"))
        if conn.execute("PRAGMA quick_check").fetchone()[0] != "ok":
            raise ValueError("Membership SQLite integrity check failed")
    finally:
        conn.close()
    membership_db.replace(args.out / "membership.sqlite")
    columns = ["callset", "samples", "input_records", "records_with_BND", "BND_alleles_skipped",
               "retained_records", "retained_ALT_alleles", "multiallelic_input_records", "retained_multiallelic_records",
               "ALT_alleles_gt10kb", "unique_matchable_keys", "repeated_matchable_keys", "unmatchable_ALT_alleles"]
    write_tsv(args.out / "callset_summary.tsv", summaries.values(), columns)
    type_rows, bin_rows, filter_rows, evidence_rows = [], [], [], []
    for label, _, saved, _ in scans:
        type_rows.extend(dict(callset=label, svtype=t, ALT_alleles=n) for t, n in sorted(saved["types"].items()))
        bin_rows.extend(dict(callset=label, **row) for row in saved["length_bins"])
        filter_rows.extend(dict(callset=label, FILTER=f, records=n) for f, n in sorted(saved["filters"].items()))
        evidence_rows.extend(dict(callset=label, evidence=e, ALT_alleles=n) for e, n in sorted(saved["evidence"].items()))
    write_tsv(args.out / "svtype_counts.tsv", type_rows, ["callset", "svtype", "ALT_alleles"])
    write_tsv(args.out / "length_bins.tsv", bin_rows, ["callset", "svtype", "length_bin", "ALT_alleles"])
    write_tsv(args.out / "filter_counts.tsv", filter_rows, ["callset", "FILTER", "records"])
    write_tsv(args.out / "allele_evidence.tsv", evidence_rows, ["callset", "evidence", "ALT_alleles"])
    overlaps = overlap_rows(memberships, summaries)
    write_tsv(args.out / "pairwise_overlap.tsv", overlaps, list(overlaps[0]))
    masks = [dict(present_in=",".join(label for i, label in enumerate(LABELS) if mask & (1 << i)),
                  svtype=t, evidence=e, unique_keys=n) for mask, t, e, n in memberships]
    write_tsv(args.out / "membership_by_type.tsv", masks, ["present_in", "svtype", "evidence", "unique_keys"])
    type_overlaps = []
    for i, j in combinations(range(4), 2):
        for kind in sorted({t for _, t, _, _ in memberships}):
            na = sum(n for mask, t, _, n in memberships if t == kind and mask & (1 << i))
            nb = sum(n for mask, t, _, n in memberships if t == kind and mask & (1 << j))
            shared = sum(n for mask, t, _, n in memberships if t == kind and mask & (1 << i) and mask & (1 << j))
            type_overlaps.append(dict(callset_a=LABELS[i], callset_b=LABELS[j], svtype=kind,
                                      unique_a=na, unique_b=nb, shared=shared, only_a=na-shared, only_b=nb-shared))
    write_tsv(args.out / "overlap_by_type.tsv", type_overlaps, ["callset_a", "callset_b", "svtype", "unique_a", "unique_b", "shared", "only_a", "only_b"])
    sample_rows = []
    for i, j in combinations(range(4), 2):
        sa, sb = scans[i][3], scans[j][3]
        sample_rows.append(dict(callset_a=LABELS[i], callset_b=LABELS[j], samples_a=len(sa), samples_b=len(sb),
                                shared_samples=len(set(sa) & set(sb)), same_order=sa == sb))
    write_tsv(args.out / "sample_overlap.tsv", sample_rows, list(sample_rows[0]))
    manifest = {
        "utc": datetime.now(timezone.utc).isoformat(), "schema": CACHE_SCHEMA,
        "parameters": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
        "script_sha256": digest_bytes(Path(__file__).read_bytes()), "bcftools_version": version,
        "input_fingerprints": {label: saved["fingerprint"] for label, _, saved, _ in scans},
        "comparison": "Per-ALT exact keys. Literal: CHROM POS REF ALT SVTYPE; symbolic additionally includes effective END and abs(SVLEN). No left normalization, sequence trimming, or fuzzy clustering.",
        "symbolic_warning": "Symbolic metadata equality does not establish sequence identity. Unmatchable alleles are excluded from overlap denominators and reported separately.",
        "selection": "All FILTER values and all lengths retained. BND/TRA/CTX, bracket/single breakends, and placeholder alleles excluded; no carrier or genotype filtering.",
        "counting": "Summary/type/length tables count ALT occurrences; overlap tables count distinct matchable keys. Genotype concordance is not assessed.",
        "length_bins": "<50bp; 50-999bp; 1000-10000bp; 10001-100000bp; >100000bp; unknown",
        "cache": "Input path/size/mtime, header hash, sample-manifest hash, algorithm schema, bcftools version, SQLite quick_check and row count; no Git commit invalidation",
        "randomness": "None",
    }
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print("\nCALLSET COUNTS (BND excluded)")
    display = ["callset", "samples", "retained_records", "retained_ALT_alleles", "ALT_alleles_gt10kb", "unmatchable_ALT_alleles"]
    print("\t".join(display))
    for row in summaries.values():
        print("\t".join(str(row[c]) for c in display))
    for filename in ("svtype_counts.tsv", "pairwise_overlap.tsv", "sample_overlap.tsv"):
        print(f"\n{filename}")
        print((args.out / filename).read_text(encoding="utf-8"), end="")
    print(f"\nSaved comparison tables and provenance to {args.out}")


if __name__ == "__main__":
    main()
