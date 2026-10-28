#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
"""Add cohort singleton SNVs without rewriting any existing scaffold record.

Use agreeing nearest informative SNVs on either side in the carrier's HiPhase
block. Unsupported candidates are excluded and recorded in a separate audit.
Standard-library implementation; BCF input and indexed VCF input require bcftools.
See docs/singleton_phasing.md for the input contract and confidence limitations.
"""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from contextlib import contextmanager
from functools import lru_cache
import gzip
import hashlib
import itertools
import json
import logging
import os
from pathlib import Path
import re
import sqlite3
import subprocess
import sys
import time

VERSION = 4  # Preserve repeated panel rows; resolve duplicate source evidence conservatively.
SEED = 1729  # Existing repository seed.
TAGS = {
    "SP_METHOD": '##INFO=<ID=SP_METHOD,Number=1,Type=String,Description="Singleton phase: PHYSICAL">',
    "SP_REASON": '##INFO=<ID=SP_REASON,Number=1,Type=String,Description="Singleton phase evidence: agreeing flanks or explicit one-sided mode">',
    "SP_ANCHORS": '##INFO=<ID=SP_ANCHORS,Number=2,Type=Integer,Description="Selected nearest SNV anchor observations supporting original and flipped orientation">',
    "SP_LEFT_POS": '##INFO=<ID=SP_LEFT_POS,Number=1,Type=Integer,Description="Nearest informative shared SNV position strictly left in carrier HiPhase block">',
    "SP_RIGHT_POS": '##INFO=<ID=SP_RIGHT_POS,Number=1,Type=Integer,Description="Nearest informative shared SNV position strictly right in carrier HiPhase block">',
    "SP_SOURCE_PS": '##INFO=<ID=SP_SOURCE_PS,Number=1,Type=String,Description="Singleton carrier original HiPhase PS">',
    "SP_SOURCE_GT": '##INFO=<ID=SP_SOURCE_GT,Number=1,Type=String,Description="Singleton carrier original HiPhase GT">',
}


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for data in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(data)
    return h.hexdigest()


def json_write(path, value):
    path = Path(path)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temp, path)


def identity(path):
    path = Path(path).resolve()
    stat = path.stat()
    return {"path": str(path), "bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def peak_memory_mib():
    if os.name == "nt":
        import ctypes
        from ctypes import wintypes
        class Counters(ctypes.Structure):
            _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD)] + [(name, ctypes.c_size_t) for name in
                ("PeakWorkingSetSize", "WorkingSetSize", "QuotaPeakPagedPoolUsage", "QuotaPagedPoolUsage", "QuotaPeakNonPagedPoolUsage", "QuotaNonPagedPoolUsage", "PagefileUsage", "PeakPagefileUsage")]
        counters = Counters()
        counters.cb = ctypes.sizeof(counters)
        current = ctypes.windll.kernel32.GetCurrentProcess
        current.restype = wintypes.HANDLE
        query = ctypes.windll.psapi.GetProcessMemoryInfo
        query.argtypes = [wintypes.HANDLE, ctypes.POINTER(Counters), wintypes.DWORD]
        if query(current(), ctypes.byref(counters), counters.cb):
            return counters.PeakWorkingSetSize / 2**20
        return None
    import resource
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value / (2**20 if sys.platform == "darwin" else 1024)


@contextmanager
def text_input(path, bcftools="bcftools", contig=None, region=None):
    """Indexed BCF/VCF via bcftools; ordinary VCF/gzip works on Windows too."""
    path = Path(path)
    indexed = any(Path(str(path) + ext).exists() for ext in (".csi", ".tbi"))
    if path.suffix == ".bcf" or (contig and indexed):
        cmd = [bcftools, "view", "--no-version", "-Ov"]
        if contig:
            cmd += ["-r", f"{contig}:{region[0]}-{region[1]}" if region else contig, "--regions-overlap", "0"]
        cmd.append(str(path))
        logging.info("VCF decode command: %s", json.dumps(cmd))
        # Never capture unbounded stderr in memory or risk a full pipe deadlock.
        with subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, encoding="utf-8") as proc:
            try:
                yield proc.stdout
            except BaseException:
                proc.terminate()
                proc.wait()
                raise
            finally:
                proc.stdout.close()
            if proc.wait() != 0:
                raise RuntimeError(f"VCF decoding failed: {cmd}")
    else:
        opener = gzip.open if path.suffix in (".gz", ".bgz") else open
        with opener(path, "rt", encoding="utf-8", newline="") as handle:
            yield handle


class VCF:
    def __init__(self, handle, contig, region=None):
        self.handle, self.contig, self.region = handle, contig, region
        self.meta = []
        for line in handle:
            line = line.rstrip("\r\n")
            if line.startswith("##"):
                self.meta.append(line)
            elif line.startswith("#CHROM\t"):
                self.columns = line
                self.samples = line.split("\t")[9:]
                if not self.samples or len(self.samples) != len(set(self.samples)):
                    raise ValueError("VCF needs unique cohort sample names")
                break
            else:
                raise ValueError("Invalid VCF header")
        else:
            raise ValueError("Missing #CHROM header")

    def records(self):
        last = 0
        for line in self.handle:
            raw = line.rstrip("\r\n")
            row = raw.split("\t")
            if len(row) != 9 + len(self.samples):
                raise ValueError("Truncated VCF record / sample column count mismatch")
            if row[0] != self.contig:
                raise ValueError(f"Input must be a {self.contig} shard or indexed: found {row[0]}")
            pos = int(row[1])
            if pos < last or pos < 1:
                raise ValueError("VCF is not sorted by POS")
            last = pos
            if self.region and not self.region[0] <= pos <= self.region[1]:
                continue
            yield row, raw


def key(row):
    return row[0], int(row[1]), row[3], row[4]


def groups(vcf):
    for pos, rows in itertools.groupby(vcf.records(), lambda x: int(x[0][1])):
        yield pos, list(rows)


def allele_groups(records):
    result = {}
    for row, raw in records:
        result.setdefault(key(row), []).append((row, raw))
    return result


def snv(row):
    return len(row[3]) == len(row[4]) == 1 and row[3] in "ACGT" and row[4] in "ACGT" and row[3] != row[4]


def anchor_compatible(a, b):
    return key(a) == key(b) and snv(a) and snv(b)


def calls(row):
    fmt = row[8].split(":")
    if "GT" not in fmt or len(fmt) != len(set(fmt)):
        raise ValueError("Missing GT or duplicate FORMAT key")
    gi = fmt.index("GT")
    pi = fmt.index("PS") if "PS" in fmt else None
    max_allele = row[4].count(",") + 1
    result = []
    for cell in row[9:]:
        fields = cell.split(":")
        if len(fields) > len(fmt):
            raise ValueError("Too many FORMAT values")
        fields += ["."] * (len(fmt) - len(fields))
        gt = fields[gi]
        alleles = parse_gt(gt)
        if any(a is not None and (a < 0 or a > max_allele) for a in alleles):
            raise ValueError("GT allele outside REF/ALT range")
        result.append((alleles, "|" in gt, fields[pi] if pi is not None else "."))
    return result


@lru_cache(maxsize=1024)
def parse_gt(gt):
    return tuple(None if a == "." else int(a) for a in re.split(r"[/|]", gt))


def hetero(call):
    return len(call[0]) == 2 and None not in call[0] and call[0][0] != call[0][1]


def singleton(cs):
    # Conservative cohort MAC: missing, haploid, or multiallelic genotypes excluded.
    if any(len(c[0]) != 2 or None in c[0] for c in cs):
        return None
    ac = sum(sum(c[0]) for c in cs)
    an = 2 * len(cs)
    if min(ac, an - ac) != 1:
        return None
    carrier = next(i for i, c in enumerate(cs) if hetero(c))
    return carrier


def call_signature(call, physical=True):
    # Homozygous phase/PS does not change the haplotypes. Panel PS is irrelevant
    # because the panel is asserted to be oriented across the chromosome.
    alleles, phased, ps = call
    return alleles, phased if hetero(call) else None, ps if physical and hetero(call) else None


def anchor_evidence(source_calls, panel_calls):
    """One vote per allele/sample/PS, or a conflict marker for ambiguous copies."""
    if not source_calls or not panel_calls:
        return
    if len(source_calls) == len(panel_calls) == 1:
        for sample, (a, b) in enumerate(zip(source_calls[0], panel_calls[0])):
            if hetero(a) and hetero(b) and a[1] and b[1] and a[2] != "." and sorted(a[0]) == sorted(b[0]):
                flipped = int(a[0] != b[0])
                yield sample, a[2], (1 - flipped, flipped), False
        return
    for sample, (sources, targets) in enumerate(zip(zip(*source_calls), zip(*panel_calls))):
        ambiguous = len({call_signature(c) for c in sources}) > 1 or len({call_signature(c, physical=False) for c in targets}) > 1
        votes = {}
        for a in sources:
            if not (hetero(a) and a[1] and a[2] != "."):
                continue
            for b in targets:
                if hetero(b) and b[1] and sorted(a[0]) == sorted(b[0]):
                    votes.setdefault(a[2], set()).add(int(a[0] != b[0]))
        for ps, orientations in votes.items():
            # Keep ambiguity at this position: never search past it for agreement.
            counts = (1, 1) if ambiguous else (int(0 in orientations), int(1 in orientations))
            yield sample, ps, counts, ambiguous


def header_id(line):
    match = re.match(r"##(INFO|FORMAT|FILTER|contig)=<ID=([^,>]+)", line)
    return match.groups() if match else None


def merged_header(scaffold, hiphase):
    return header_merge_plan(scaffold, hiphase)[0]


def header_merge_plan(scaffold, hiphase):
    """Keep every panel definition; rename incompatible source annotations."""
    result = list(scaffold.meta)
    definitions = {header_id(x): x for x in result if header_id(x)}
    occupied = {header_id(x) for x in scaffold.meta + hiphase.meta if header_id(x)}
    renames = {"INFO": {}, "FORMAT": {}}
    for name in TAGS:
        if ("INFO", name) in occupied:
            raise ValueError(f"Reserved INFO field already exists: {name}")
    for line in hiphase.meta:
        tag = header_id(line)
        if tag is None:
            continue  # Keep scaffold provenance; full input headers hashed in audit.
        if tag in definitions and definitions[tag] != line:
            # Description-only differences do not alter the data schema.
            attrs = ("Number", "Type") if tag[0] in ("INFO", "FORMAT") else ("length",)
            for attr in attrs:
                pattern = rf"(?:<|,){attr}=([^,>]+)"
                a, b = re.search(pattern, definitions[tag]), re.search(pattern, line)
                if (a.group(1) if a else None) != (b.group(1) if b else None):
                    if tag[0] not in renames or tag == ("FORMAT", "GT"):
                        raise ValueError(f"Incompatible header definition: {tag}; panel={definitions[tag]}; HiPhase={line}")
                    alias = f"SP_HIPHASE_{tag[1]}"
                    suffix = 2
                    while (tag[0], alias) in occupied:
                        alias = f"SP_HIPHASE_{tag[1]}_{suffix}"
                        suffix += 1
                    occupied.add((tag[0], alias))
                    renames[tag[0]][tag[1]] = alias
                    result.append(line.replace(f"<ID={tag[1]}", f"<ID={alias}", 1))
                    logging.warning("Header schema conflict %s -> %s; panel=%s; HiPhase=%s", tag, alias, definitions[tag], line)
                    break
        elif tag not in definitions:
            result.append(line)
            definitions[tag] = line
    for name, line in TAGS.items():
        if ("INFO", name) in definitions:
            raise ValueError(f"Reserved INFO field already exists: {name}")
        result.append(line)
    result.append(scaffold.columns)
    return "\n".join(result) + "\n", renames


def rename_source_annotations(row, renames):
    """Apply only to additions, after extracting GT/PS and removing physical PS/PQ."""
    row = list(row)
    if row[7] != ".":
        values = []
        for item in row[7].split(";"):
            name, separator, value = item.partition("=")
            values.append(renames["INFO"].get(name, name) + separator + value)
        row[7] = ";".join(values)
    row[8] = ":".join(renames["FORMAT"].get(name, name) for name in row[8].split(":"))
    return row


def write_header(handle, vcf):
    handle.write("\n".join(vcf.meta + [vcf.columns]) + "\n")


def annotate(row, carrier, gt, method, reason, anchors, left, right):
    row = list(row)
    source_ps = info_value(row, "SP_SOURCE_PS")
    source_gt = info_value(row, "SP_SOURCE_GT")
    if source_gt is None:
        original = calls(row[:9] + [row[9 + carrier]])[0]
        source_gt = ("|" if original[1] else "/").join(map(str, original[0]))
        source_ps = original[2]
    old = row[8].split(":")
    # HiPhase PS/PQ describes a different coordinate system. Never export it as
    # the final scaffold PS/PQ. Original scaffold FORMAT is never touched.
    keep = [i for i, field in enumerate(old) if field not in ("PS", "PQ")]
    new = [old[i] for i in keep]
    gi = new.index("GT")
    for sample, cell in enumerate(row[9:]):
        values = cell.split(":") + ["."] * (len(old) - len(cell.split(":")))
        values = [values[i] for i in keep]
        if sample == carrier:
            values[gi] = gt
        row[9 + sample] = ":".join(values)
    row[8] = ":".join(new)
    info = [] if row[7] == "." else row[7].split(";")
    info = [x for x in info if x.split("=")[0] not in TAGS]
    info += [f"SP_METHOD={method}", f"SP_REASON={reason}", f"SP_ANCHORS={anchors[0]},{anchors[1]}"]
    info += [f"SP_SOURCE_PS={source_ps}", f"SP_SOURCE_GT={source_gt}"]
    info += [f"SP_LEFT_POS={left[0] if left else '.'}", f"SP_RIGHT_POS={right[0] if right else '.'}"]
    row[7] = ";".join(info)
    return row


def info_value(row, name):
    for item in row[7].split(";"):
        if item.startswith(name + "="):
            return item.split("=", 1)[1]
    return None


def nearest_anchors(db, sample, ps, pos):
    # PK(sample,ps,pos) supports two bounded index seeks, not a block-wide scan.
    left = db.execute("SELECT pos,same,flipped FROM anchors WHERE sample=? AND ps=? AND pos<? ORDER BY pos DESC LIMIT 1", (sample, ps, pos)).fetchone()
    right = db.execute("SELECT pos,same,flipped FROM anchors WHERE sample=? AND ps=? AND pos>? ORDER BY pos ASC LIMIT 1", (sample, ps, pos)).fetchone()
    return left, right


def anchor_decision(left, right, min_anchors):
    selected = [a for a in (left, right) if a is not None]
    counts = tuple(sum(a[i] for a in selected) for i in (1, 2))
    # Do not search past a conflicting informative marker for an agreement.
    if min(counts) > 0:
        return "CONFLICTING_NEAREST_ANCHORS", counts
    if not selected:
        return "NO_INFORMATIVE_ANCHOR", counts
    if len(selected) < min_anchors:
        return "MISSING_LEFT_ANCHOR" if left is None else "MISSING_RIGHT_ANCHOR", counts
    return ("AGREEING_FLANKS" if len(selected) == 2 else "ONE_SIDED_ANCHOR"), counts


def prepare(job, directory, args, expected_samples, sparse=False, defer_phase=False):
    """One HiPhase pass; disk-backed block evidence, then candidate-only pass."""
    stats = Counter()
    started = time.monotonic()
    db_path = directory / "anchors.sqlite"
    db_path.unlink(missing_ok=True)
    db = sqlite3.connect(db_path)
    db.execute("PRAGMA journal_mode=WAL")
    db.execute("PRAGMA cache_size=-65536")
    db.execute("CREATE TABLE anchors(sample INTEGER, ps TEXT, pos INTEGER, same INTEGER, flipped INTEGER, PRIMARY KEY(sample,ps,pos)) WITHOUT ROWID")
    db.execute("CREATE TABLE candidate_conflicts(pos INTEGER, ref TEXT, alt TEXT, PRIMARY KEY(pos,ref,alt)) WITHOUT ROWID")
    pending = {}

    def flush():
        db.executemany("INSERT INTO anchors VALUES(?,?,?,?,?) ON CONFLICT(sample,ps,pos) DO UPDATE SET same=same+excluded.same, flipped=flipped+excluded.flipped",
                       [(sample, ps, pos, counts[0], counts[1]) for (sample, ps, pos), counts in pending.items()])
        db.commit()
        pending.clear()

    region, core = job.get("region"), job.get("core")
    with text_input(job["scaffold"], args.bcftools, job["contig"], region) as sf, text_input(job["hiphase"], args.bcftools, job["contig"], region) as hf:
        scaffold, hiphase = VCF(sf, job["contig"], region), VCF(hf, job["contig"], region)
        if scaffold.samples != hiphase.samples or scaffold.samples != expected_samples:
            raise ValueError("Sample identity/order differs from the manifest or between inputs")
        header, renames = header_merge_plan(scaffold, hiphase)
        json_write(directory / "header_merge.json", renames)
        (directory / "merged.header.vcf").write_text(header, encoding="utf-8")
        stats["samples"] = len(scaffold.samples)
        stats["renamed_info_fields"] = len(renames["INFO"])
        stats["renamed_format_fields"] = len(renames["FORMAT"])
        stats["scaffold_header_sha256"] = hashlib.sha256(("\n".join(scaffold.meta + [scaffold.columns]) + "\n").encode()).hexdigest()
        stats["hiphase_header_sha256"] = hashlib.sha256(("\n".join(hiphase.meta + [hiphase.columns]) + "\n").encode()).hexdigest()
        scaffold_hash, hiphase_hash = hashlib.sha256(), hashlib.sha256()
        sg = iter(groups(scaffold))
        current = next(sg, None)
        candidate_name = "raw_candidates.jsonl.gz" if sparse else "raw_candidates.vcf.gz"
        with gzip.open(directory / candidate_name, "wt", encoding="utf-8", newline="\n") as candidates, gzip.open(directory / "duplicate_sites.tsv.gz", "wt", encoding="utf-8", newline="\n") as duplicates:
            if not sparse:
                write_header(candidates, hiphase)
            duplicates.write("input\tchrom\tpos\tref\talt\trecord_count\taction\n")

            def record_duplicate(source, rows, action):
                row = rows[0][0]
                if len(rows) < 2 or (core and not core[0] <= int(row[1]) <= core[1]):
                    return
                duplicates.write("\t".join(map(str, [source, *key(row), len(rows), action])) + "\n")
                stats[f"{source}_duplicate_keys"] += 1
                stats[f"{source}_extra_duplicate_rows"] += len(rows) - 1

            def save_scaffold(group):
                for rows in allele_groups(group[1]).values():
                    record_duplicate("scaffold", rows, "PANEL_ROWS_PRESERVED")
                for row, raw in group[1]:
                    if not core or core[0] <= int(row[1]) <= core[1]:
                        scaffold_hash.update((raw + "\n").encode())
                        stats["scaffold_records"] += 1
                    cs = calls(row)
                    if any(hetero(c) and not c[1] for c in cs):
                        raise ValueError(f"Unphased scaffold heterozygote at {key(row)}")

            next_progress = 100000
            for pos, hr in groups(hiphase):
                while current is not None and current[0] < pos:
                    save_scaffold(current)
                    current = next(sg, None)
                shared = allele_groups(current[1]) if current and current[0] == pos else {}
                for row, raw in hr:
                    stats["hiphase_records"] += 1
                    hiphase_hash.update((raw + "\n").encode())
                for allele_key, source_rows in allele_groups(hr).items():
                    row, raw = source_rows[0]
                    if not snv(row):
                        record_duplicate("hiphase", source_rows, "IGNORED_NON_SNV")
                        stats["excluded_incompatible_shared_anchor" if allele_key in shared else "excluded_non_biallelic_snv"] += len(source_rows)
                        continue
                    eligible = [(r, raw) for r, raw in source_rows if r[6] in ("PASS", ".")]
                    if allele_key in shared:
                        record_duplicate("hiphase", source_rows, "ANCHOR_EVIDENCE_COMBINED")
                        stats["shared_anchor_sites"] += 1
                        stats["shared_snvs"] += 1
                        sources = [calls(r) for r, _ in eligible]
                        panels = [calls(r) for r, _ in shared[allele_key] if r[6] in ("PASS", ".")]
                        for sample, ps, counts, ambiguous in anchor_evidence(sources, panels):
                            block = (sample, ps, pos)
                            if block not in pending:
                                pending[block] = [0, 0]
                            pending[block][0] += counts[0]
                            pending[block][1] += counts[1]
                            stats["anchor_observations"] += sum(counts)
                            stats["ambiguous_duplicate_anchor_blocks"] += int(ambiguous)
                        if len(pending) >= 50000:
                            flush()
                        continue
                    if core and not core[0] <= pos <= core[1]:
                        continue  # Halo records supply anchors, never candidates.
                    if not eligible:
                        record_duplicate("hiphase", source_rows, "FILTERED_SNV")
                        stats["excluded_filtered_snv"] += 1
                        continue
                    source_calls = [calls(r) for r, _ in eligible]
                    carriers = [singleton(cs) for cs in source_calls]
                    first = next((i for i, carrier in enumerate(carriers) if carrier is not None), None)
                    if first is not None:
                        row, raw = eligible[first]
                        cs, carrier = source_calls[first], carriers[first]
                        conflict = any(any(call_signature(a) != call_signature(b) for a, b in zip(cs, other)) for other in source_calls if other is not cs)
                        if conflict:
                            db.execute("INSERT INTO candidate_conflicts VALUES(?,?,?)", (pos, row[3], row[4]))
                        record_duplicate("hiphase", source_rows, "CONFLICTING_DUPLICATE_SNV" if conflict else "SINGLETON_COPIES_COLLAPSED")
                        if sparse:
                            background = "0|0" if sum(sum(c[0]) for c in cs) == 1 else "1|1"
                            candidates.write(json.dumps({"row": row[:9] + [row[9 + carrier]], "carrier": carrier, "background": background}, separators=(",", ":")) + "\n")
                        else:
                            candidates.write(raw + "\n")
                        stats["singleton_candidates"] += 1
                    else:
                        record_duplicate("hiphase", source_rows, "NON_SINGLETON_OR_INCOMPLETE")
                        stats["excluded_non_singleton_or_incomplete_diploid"] += 1
                if current and current[0] == pos:
                    save_scaffold(current)
                    current = next(sg, None)
                if stats["hiphase_records"] >= next_progress:
                    logging.info("%s scanned %s HiPhase rows", job["contig"], stats["hiphase_records"])
                    next_progress = (stats["hiphase_records"] // 100000 + 1) * 100000
            while current:
                save_scaffold(current)
                current = next(sg, None)
        flush()
        if not core and (not stats["scaffold_records"] or not stats["hiphase_records"]):
            raise ValueError("An input has no records for the requested contig")
        if not stats["shared_anchor_sites"]:
            logging.warning("No local exact shared SNV anchors; other chunk indexes will be queried." if defer_phase else
                            "No exact shared SNV anchors; candidates will be excluded. Check build/allele representations.")
    stats["anchor_positions"] = db.execute("SELECT count(*) FROM anchors").fetchone()[0]
    db.close()
    stats["prepare_seconds"] = time.monotonic() - started
    stats["scaffold_record_sha256"] = scaffold_hash.hexdigest()
    stats["hiphase_record_sha256"] = hiphase_hash.hexdigest()
    if not defer_phase:
        stats.update(resolve_candidates(job, directory, args, expected_samples, sparse))
    return dict(stats)


def resolve_candidates(job, directory, args, expected_samples, sparse=False, other_indexes=()):
    """Resolve compact candidates; indexes in other chunks close halo boundaries."""
    stats = Counter()
    renames = json.loads((directory / "header_merge.json").read_text(encoding="utf-8"))
    paths = [directory / "anchors.sqlite", *other_indexes]
    databases = [sqlite3.connect(Path(p).resolve().as_uri() + "?mode=ro", uri=True) for p in paths]
    try:
        for db in databases:
            db.execute("PRAGMA cache_size=-8192")
        candidate_name = "raw_candidates.jsonl.gz" if sparse else "raw_candidates.vcf.gz"
        addition_name = "added_singletons.jsonl.gz" if sparse else "added_singletons.vcf.gz"
        with text_input(directory / candidate_name) as raw_handle, gzip.open(directory / addition_name, "wt", encoding="utf-8", newline="\n") as output, gzip.open(directory / "excluded_singletons.tsv.gz", "wt", encoding="utf-8", newline="\n") as excluded:
            if sparse:
                candidates = (json.loads(line) for line in raw_handle)
            else:
                candidates = ({"row": row, "carrier": singleton(calls(row))} for row, _ in VCF(raw_handle, job["contig"]).records())
                output.write((directory / "merged.header.vcf").read_text(encoding="utf-8"))
            excluded.write("chrom\tpos\tref\talt\tsample\tsource_ps\treason\tleft_pos\tright_pos\tsame\tflipped\n")
            for candidate in candidates:
                row, carrier = candidate["row"], candidate["carrier"]
                row_carrier = 0 if sparse else carrier
                alleles, phased, ps = calls(row)[row_carrier]
                duplicate_conflict = databases[0].execute("SELECT 1 FROM candidate_conflicts WHERE pos=? AND ref=? AND alt=?", (int(row[1]), row[3], row[4])).fetchone() is not None
                can_anchor = phased and ps != "." and not duplicate_conflict
                left, right = nearest_anchors(databases[0], carrier, ps, int(row[1])) if can_anchor else (None, None)
                if can_anchor and (left is None or right is None) and other_indexes:
                    stats["cross_chunk_lookups"] += 1
                    # Every index contains a complete contiguous interval. A flank
                    # found locally is already nearest. Missing flanks need all
                    # other indexes; overlapping copies are selected, not summed.
                    pairs = [nearest_anchors(db, carrier, ps, int(row[1])) for db in databases[1:]]
                    if left is None:
                        left = max((a for a, _ in pairs if a is not None), default=None, key=lambda a: a[0])
                    if right is None:
                        right = min((b for _, b in pairs if b is not None), default=None, key=lambda b: b[0])
                if duplicate_conflict:
                    reason, counts = "CONFLICTING_DUPLICATE_SNV", (0, 0)
                elif not phased or ps == ".":
                    reason, counts = "NO_PHYSICAL_PS", (0, 0)
                else:
                    reason, counts = anchor_decision(left, right, args.min_anchors)
                if reason in ("AGREEING_FLANKS", "ONE_SIDED_ANCHOR"):
                    gt = "|".join(map(str, alleles[::-1] if counts[1] else alleles))
                    annotated = annotate(row, row_carrier, gt, "PHYSICAL", reason, counts, left, right)
                    annotated = rename_source_annotations(annotated, renames)
                    if sparse:
                        output.write(json.dumps({"fields": annotated[:8], "carrier": carrier, "gt": gt, "background": candidate["background"]}, separators=(",", ":")) + "\n")
                    else:
                        output.write("\t".join(annotated) + "\n")
                    stats["physical"] += 1
                else:
                    excluded.write("\t".join(map(str, [row[0], row[1], row[3], row[4], "." if duplicate_conflict else expected_samples[carrier], "." if duplicate_conflict else ps, reason, left[0] if left else ".", right[0] if right else ".", *counts])) + "\n")
                    stats["excluded"] += 1
                stats[reason.lower()] += 1
    finally:
        for db in databases:
            db.close()
    return stats


def finalize(job, directory, args, expected_samples, stats):
    added = directory / "added_singletons.vcf.gz"
    counts = {method: stats.get(method, 0) for method in ("physical", "excluded")}
    output = directory / "augmented.vcf.gz"
    temp = directory / "augmented.partial.vcf.gz"
    h = hashlib.sha256()
    n = 0
    with text_input(job["scaffold"], args.bcftools, job["contig"]) as sh, text_input(added) as ah, gzip.open(temp, "wt", encoding="utf-8", newline="\n") as out:
        sv, av = VCF(sh, job["contig"]), VCF(ah, job["contig"])
        if sv.samples != expected_samples:
            raise ValueError("Scaffold sample manifest changed")
        out.write((directory / "merged.header.vcf").read_text(encoding="utf-8"))
        ai = iter(av.records())
        a = next(ai, None)
        for row, raw in sv.records():
            while a is not None and int(a[0][1]) < int(row[1]):
                out.write(a[1] + "\n")
                a = next(ai, None)
            out.write(raw + "\n")
            h.update((raw + "\n").encode())
            n += 1
        while a:
            out.write(a[1] + "\n")
            a = next(ai, None)
    if h.hexdigest() != stats["scaffold_record_sha256"] or n != stats["scaffold_records"]:
        raise ValueError("Scaffold changed since preparation")
    audit = audit_preservation(job["scaffold"], temp, job["contig"], args.bcftools)
    if audit["added_records"] != counts["physical"] or sum(counts.values()) != stats.get("singleton_candidates", 0):
        raise ValueError("Unexpected added-record count")
    os.replace(temp, output)
    return {"final_counts": dict(counts), "preservation": audit, "output": str(output)}


def audit_preservation(scaffold_path, output_path, contig, bcftools="bcftools"):
    """Independent streaming check: every original header/record survives exactly."""
    with text_input(scaffold_path, bcftools, contig) as sh, text_input(output_path, bcftools, contig) as oh:
        sv, ov = VCF(sh, contig), VCF(oh, contig)
        if sv.samples != ov.samples or ov.meta[:len(sv.meta)] != sv.meta:
            raise ValueError("Scaffold header/sample order changed")
        originals = iter(groups(sv))
        current = next(originals, None)
        preserved = added = 0
        h = hashlib.sha256()
        for pos, rows in groups(ov):
            if current and current[0] < pos:
                raise ValueError("Scaffold record missing from final output")
            expected = current[1] if current and current[0] == pos else []
            original_keys = {key(r) for r, _ in expected}
            added_keys = set()
            seen = []
            for row, raw in rows:
                if key(row) in original_keys:
                    if len(seen) >= len(expected) or raw != expected[len(seen)][1]:
                        raise ValueError(f"Original scaffold record changed: {key(row)}")
                    seen.append(raw)
                    h.update((raw + "\n").encode())
                    preserved += 1
                else:
                    if key(row) in added_keys:
                        raise ValueError(f"Duplicate singleton addition: {key(row)}")
                    added_keys.add(key(row))
                    cs = calls(row) if snv(row) else []
                    carrier = singleton(cs) if cs else None
                    if carrier is None or not cs[carrier][1] or info_value(row, "SP_METHOD") != "PHYSICAL":
                        raise ValueError("Non-singleton or unaudited record added")
                    added += 1
            if seen != [raw for _, raw in expected]:
                raise ValueError("Scaffold records missing/reordered")
            if expected:
                current = next(originals, None)
        if current:
            raise ValueError("Trailing scaffold records missing")
    return {"preserved_records": preserved, "added_records": added, "scaffold_record_sha256": h.hexdigest(), "exact_record_preservation": True}


def job_run(payload):
    job, root, params, samples, manifest_hash = payload
    args = argparse.Namespace(**params)
    directory = Path(root) / job["contig"]
    directory.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(filename=directory / "workflow.log", level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", force=True)
    signature = {"algorithm_version": VERSION, "manifest_sha256": manifest_hash,
                 "sample_manifest_sha256": hashlib.sha256(("\n".join(samples) + "\n").encode()).hexdigest(),
                 "inputs": {k: identity(job[k]) for k in ("scaffold", "hiphase")},
                 "parameters": {k: v for k, v in params.items() if k not in ("workers", "memory_gb")}, "job": job}
    lock = directory / ".lock"
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        raise RuntimeError(f"Job locked: {lock}. Remove only after verifying no process is using this job.") from None
    try:
        with os.fdopen(fd, "w") as handle:
            handle.write(str(os.getpid()))
        done = directory / "done.json"
        if done.exists():
            previous = json.loads(done.read_text())
            if previous["signature"] == signature and all((directory / name).exists() and digest(directory / name) == sha for name, sha in previous["artifacts"].items()):
                logging.info("Validated cached outputs; skipped")
                return {"contig": job["contig"], "cached": True, **previous["result"]}
        done.unlink(missing_ok=True)
        started = time.monotonic()
        logging.info("START %s", json.dumps({"signature": signature, "execution_parameters": params, "python": sys.version}, sort_keys=True))
        prep_done = directory / "prepare.done.json"
        prep_names = ("added_singletons.vcf.gz", "excluded_singletons.tsv.gz", "duplicate_sites.tsv.gz", "merged.header.vcf", "header_merge.json", "anchors.sqlite")
        prep_signature = {"algorithm_version": VERSION, "inputs": {k: v for k, v in signature["inputs"].items() if k != "map"},
                          "contig": job["contig"], "sample_manifest_sha256": signature["sample_manifest_sha256"], "min_anchors": args.min_anchors}
        if prep_done.exists() and (previous := json.loads(prep_done.read_text()))["signature"] == prep_signature and all((directory / name).exists() and digest(directory / name) == sha for name, sha in previous["artifacts"].items()):
            stats = previous["stats"]
        else:
            stats = prepare(job, directory, args, samples)
            json_write(prep_done, {"signature": prep_signature, "stats": stats, "artifacts": {name: digest(directory / name) for name in prep_names}})
            (directory / "raw_candidates.vcf.gz").unlink(missing_ok=True)
        result = {"prepare": stats, **finalize(job, directory, args, samples, stats)}
        result["elapsed_seconds"] = time.monotonic() - started
        result["process_peak_memory_mib"] = peak_memory_mib()
        # Check source identity again to catch mutations during the run.
        if signature["inputs"] != {k: identity(job[k]) for k in signature["inputs"]}:
            raise ValueError("Input files changed during processing")
        json_write(directory / "audit.json", {"signature": signature, "execution_parameters": params, **result})
        artifacts = {name: digest(directory / name) for name in ("augmented.vcf.gz", "added_singletons.vcf.gz", "excluded_singletons.tsv.gz", "duplicate_sites.tsv.gz", "anchors.sqlite", "header_merge.json", "audit.json")}
        json_write(done, {"signature": signature, "result": result, "artifacts": artifacts})
        logging.info("COMPLETE %s", json.dumps(result))
        return {"contig": job["contig"], "cached": False, **result}
    except BaseException:
        logging.exception("Job failed; no valid new completion marker written")
        raise
    finally:
        lock.unlink(missing_ok=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=4, help="Maximum worker processes: chunks when enabled, otherwise chromosomes")
    parser.add_argument("--chunks-per-chrom", type=int, default=1, help="Disjoint coordinate cores per chromosome; chr1 runner defaults to 4")
    parser.add_argument("--overlap-bp", type=int, default=1000000, help="Anchor halo on each side; missing flanks query other chunk indexes")
    parser.add_argument("--memory-gb", type=float, default=26, help="Declared TOTAL scheduler memory budget, not an OS-enforced cap")
    parser.add_argument("--min-anchors", type=int, choices=(1, 2), default=2, help="2 requires agreeing left/right SNVs; 1 explicitly permits a single-sided anchor")
    parser.add_argument("--bcftools", default="bcftools")
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args(argv)
    if args.workers < 1 or args.workers > 32 or args.min_anchors < 1 or args.memory_gb <= 0 or not 1 <= args.chunks_per_chrom <= 32 or args.overlap_bp < 0:
        parser.error("Invalid resource or anchor-threshold parameter")
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    if not manifest.get("genome_build") or not manifest.get("jobs") or not manifest.get("samples"):
        parser.error("Manifest needs genome_build, samples (ordered sample ID list), and jobs")
    if manifest.get("scaffold_phase_scope") != "chromosome":
        parser.error("Set scaffold_phase_scope=chromosome only for globally oriented statistical haplotypes; independent PS blocks are not a valid global scaffold")
    samples = manifest["samples"]
    if not isinstance(samples, list) or len(samples) < 2 or any(not isinstance(x, str) or not x for x in samples) or len(samples) != len(set(samples)):
        parser.error("Manifest samples must be a nonempty unique ordered list")
    jobs = manifest["jobs"]
    contigs = [job["contig"] for job in jobs]
    if len(contigs) != len(set(contigs)) or any(not re.fullmatch(r"[A-Za-z0-9_.-]+", c) or c in (".", "..") for c in contigs):
        parser.error("Jobs need distinct safe contig names")
    for job in jobs:
        for name in ("scaffold", "hiphase"):
            path = Path(job[name])
            job[name] = str((args.manifest.resolve().parent / path).resolve())
            if not Path(job[name]).is_file():
                parser.error(f"Missing {name}: {job[name]}")
            if Path(job[name]).is_relative_to((args.out / job["contig"]).resolve()):
                parser.error("Input files must be outside the generated job directory")
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "summary.json").unlink(missing_ok=True)
    params = vars(args).copy()
    del params["manifest"], params["out"]
    manifest_hash = hashlib.sha256(json.dumps(manifest, sort_keys=True).encode()).hexdigest()
    payloads = [(job, str(args.out.resolve()), params, samples, manifest_hash) for job in jobs]
    if args.chunks_per_chrom > 1:
        from singleton_chunks import job_run as chunk_job_run
        results = [chunk_job_run(p) for p in payloads]
    elif args.workers == 1 or len(jobs) == 1:
        results = [job_run(p) for p in payloads]
    else:
        with ProcessPoolExecutor(max_workers=min(args.workers, len(jobs))) as pool:
            results = list(pool.map(job_run, payloads))
    json_write(args.out / "summary.json", {"genome_build": manifest["genome_build"], "samples": len(samples), "results": results})
    print(f"Finished {len(results)} contig(s); audit and provenance: {args.out.resolve() / 'summary.json'}")


if __name__ == "__main__":
    main()
