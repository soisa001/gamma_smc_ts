# /// script
# requires-python = ">=3.12,<3.13"
# dependencies = ["truvari==5.5.0", "pywfa==0.5.1"]
# ///
"""Independent 90% and 70% Truvari re-collapse experiments on stringent SVs.

Separate each type and representation. DEL uses size plus reciprocal overlap;
only sequence-resolved INS can collapse. No uploads.
"""

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
import copy
import csv
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

# Four independent jobs, each with one numerical-library thread.
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "1"

import pysam
import truvari

SCHEMA = 2
OLD_COLLAPSE = ("CollapseId", "NumCollapsed", "NumConsolidated", "CollapseStart", "CollapseEnd", "CollapseSize")


def digest(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def write_json(path, obj):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(obj, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def run(command, log):
    with log.open("a", encoding="utf-8") as handle:
        handle.write("COMMAND " + json.dumps([str(c) for c in command]) + "\n")
        handle.flush()
        subprocess.run([str(c) for c in command], stdout=handle, stderr=handle, check=True)


def files_valid(root, hashes):
    return bool(hashes) and all((root / f).is_file() and digest(root / f) == h for f, h in hashes.items())


def classify(record):
    if not record.alts or len(record.alts) != 1:
        raise ValueError("Expected biallelic ALT records; do not split genotypes automatically")
    alt = record.alts[0]
    symbolic = alt.startswith("<") and alt.endswith(">")
    kind = str(record.info.get("SVTYPE", alt[1:-1] if symbolic else "UNKNOWN")).upper()
    base = kind.split(":")[0]
    if base in ("BND", "TRA", "CTX") or any(x in alt for x in "[]") or alt.startswith(".") or alt.endswith("."):
        return "BND", "excluded", 0
    raw = record.info.get("SVLEN")
    if isinstance(raw, tuple):
        raw = raw[0] if len(raw) == 1 else None
    length = abs(int(raw)) if raw else None
    if length is None and not symbolic:
        length = abs(len(record.ref) - len(alt))
    if not length and base in ("DEL", "DUP", "INV"):
        length = record.stop - record.pos
    if not length or length <= 10000:
        raise ValueError(f"Expected a known >10 kb length at {record.chrom}:{record.pos}")
    if base not in ("DEL", "INS", "DUP", "INV"):
        return kind, "unresolved_passthrough", length
    # Full symbolic label keeps e.g. DUP:TANDEM separate from generic DUP.
    if symbolic:
        label = alt[1:-1].upper()
        # pysam >=0.24 can infer stop from SVLEN even when END is invalid.
        # Inspect the coordinate span without that inference before accepting it.
        span_record = record.copy()
        if "SVLEN" in span_record.info:
            del span_record.info["SVLEN"]
        bad_span = base in ("DEL", "DUP", "INV") and (
            span_record.stop <= record.pos or span_record.stop != record.stop)
        if label.split(":")[0] != base or bad_span:
            return kind, "unresolved_passthrough", length
        return label, "symbolic_size_position", length
    # Sequence comparison is conservative: ambiguous bases and replacements go
    # to passthrough, never to the symbolic size-only branch.
    if base not in ("DEL", "INS") or not re.fullmatch("[ACGTacgt]+", record.ref + alt):
        return kind, "unresolved_passthrough", length
    ref, allele = record.ref.upper(), alt.upper()
    start = 0
    while start < min(len(ref), len(allele)) and ref[start] == allele[start]:
        start += 1
    ref, allele = ref[start:], allele[start:]
    stop = 0
    while stop < min(len(ref), len(allele)) and ref[-stop-1] == allele[-stop-1]:
        stop += 1
    if stop:
        ref, allele = ref[:-stop], allele[:-stop]
    payload = allele if base == "INS" and not ref else ref if base == "DEL" and not allele else None
    mode = "sequence_resolved" if payload and len(payload) == length else "unresolved_passthrough"
    return kind, mode, length


def genotype_totals(record, dosage, carriers):
    ac = an = ns = nc = 0
    for index, sample in enumerate(record.samples.values()):
        gt = sample.get("GT") or ()
        if any(a not in (None, 0, 1) for a in gt):
            raise ValueError("Unexpected allele index")
        called = sum(a is not None for a in gt)
        copies = gt.count(1)
        dosage[index] += copies
        carriers[index] += int(copies > 0)
        ac += copies
        an += called
        ns += int(called > 0)
        nc += int(copies > 0)
    return ac, an, ns, nc


def totals_hash(values):
    return hashlib.sha256(json.dumps(values, separators=(",", ":")).encode()).hexdigest()


def prepare(source, out):
    folder = out / "partitions"
    folder.mkdir(parents=True, exist_ok=True)
    receipt = folder / "manifest.json"
    with pysam.VariantFile(source) as vcf:
        sample_names = list(vcf.header.samples)
        if not sample_names or len(set(sample_names)) != len(sample_names):
            raise ValueError("Sample manifest empty or contains duplicate names")
        signature = dict(path=str(source), size=source.stat().st_size, mtime_ns=source.stat().st_mtime_ns,
                         header_sha256=hashlib.sha256(str(vcf.header).encode()).hexdigest(),
                         samples_sha256=totals_hash(sample_names), schema=SCHEMA)
        if receipt.exists():
            saved = json.loads(receipt.read_text())
            if saved["input"] == signature and files_valid(folder, saved["hashes"]):
                print("Reusing validated partitions and genotype totals", flush=True)
                return saved
        header = vcf.header.copy()
        for key, number, typ, description in [
            ("RC_ORIGINAL_ID", "1", "String", "Original record ID before this recollapse experiment"),
            ("RC_SOURCE_DISCOVERY", "1", "Integer", "Representative source discovery count before recollapse"),
            ("RC_N_CARRIERS", "1", "Integer", "Current subset samples with at least one ALT allele"),
            ("AC", "A", "Integer", "ALT allele count recalculated from GT"),
            ("AN", "1", "Integer", "Number of called alleles"),
            ("AF", "A", "Float", "AC/AN among called alleles; not cohort frequency for discovery-only data"),
            ("NS", "1", "Integer", "Samples with any called allele"),
        ]:
            if key not in header.info:
                header.add_line(f'##INFO=<ID={key},Number={number},Type={typ},Description="{description}">')
        handles, metadata, ends, original_rows = {}, {}, defaultdict(dict), []
        skipped = Counter()
        try:
            for row, record in enumerate(vcf, 1):
                kind, mode, length = classify(record)
                if mode == "excluded":
                    skipped[kind] += 1
                    continue
                slug = re.sub(r"[^A-Za-z0-9_.-]", "_", kind) + "." + mode
                if slug not in handles:
                    filename = slug + ".vcf.gz"
                    handles[slug] = pysam.VariantFile(folder / filename, "wz", header=header)
                    metadata[slug] = dict(svtype=kind, mode=mode, filename=filename, records=0,
                                          gt50kb=0, max_length=0, dosage=[0]*len(sample_names), carriers=[0]*len(sample_names))
                m = metadata[slug]
                record.translate(header)
                original_id = record.id or "."
                record.id = f"RC{row:08d}"
                record.info["RC_ORIGINAL_ID"] = original_id
                for key in OLD_COLLAPSE:
                    if key in record.info:
                        del record.info[key]
                if "N_DISCOVERY_SAMPLES" in record.info:
                    record.info["RC_SOURCE_DISCOVERY"] = record.info["N_DISCOVERY_SAMPLES"]
                    del record.info["N_DISCOVERY_SAMPLES"]
                ac, an, ns, nc = genotype_totals(record, m["dosage"], m["carriers"])
                record.info["AC"], record.info["AN"] = (ac,), an
                record.info["AF"], record.info["NS"] = ((ac/an if an else 0.0),), ns
                record.info["RC_N_CARRIERS"] = nc
                # Confirm that our >10 kb eligibility agrees with Truvari's size.
                actual_size = truvari.VariantRecord(record).var_size()
                if actual_size != length:
                    raise ValueError(f"Truvari size disagrees with audit length for {record.id}")
                handles[slug].write(record)
                m["records"] += 1
                m["gt50kb"] += int(length > 50000)
                m["max_length"] = max(m["max_length"], length)
                ends[slug][record.chrom] = max(ends[slug].get(record.chrom, 0), record.stop, record.pos+length+1)
                original_rows.append((record.id, original_id, kind, mode, length))
                if row % 1000 == 0:
                    print(f"Partitioned {row:,} records", flush=True)
        finally:
            for handle in handles.values():
                handle.close()
    if source.stat().st_size != signature["size"] or source.stat().st_mtime_ns != signature["mtime_ns"]:
        raise ValueError("Input changed during partitioning")
    hashes = {}
    for slug, m in metadata.items():
        path = folder / m["filename"]
        # Sorting also provides a check that every record can be decoded.
        sorted_path = folder / (slug + ".sorted.vcf.gz")
        run(["bcftools", "sort", "-Oz", "-o", sorted_path, path], folder / "prepare.log")
        sorted_path.replace(path)
        pysam.tabix_index(str(path), preset="vcf", force=True)
        bed = folder / (slug + ".bed")
        bed.write_text("".join(f"{chrom}\t0\t{end+1000}\n" for chrom, end in ends[slug].items()))
        for f in (path, Path(str(path)+".tbi"), bed):
            hashes[f.name] = digest(f)
        m["bed"] = bed.name
    with (folder / "source_records.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["audit_ID", "original_ID", "SVTYPE", "mode", "length_bp"])
        writer.writerows(original_rows)
    hashes["source_records.tsv"] = digest(folder / "source_records.tsv")
    saved = dict(input=signature, samples=sample_names, partitions=metadata, hashes=hashes, skipped=dict(skipped))
    write_json(receipt, saved)
    return saved


def validate_and_finalize(raw, final, partition, samples):
    dosage, carriers = [0]*len(samples), [0]*len(samples)
    count = gt50 = 0
    ids = set()
    with pysam.VariantFile(raw) as vcf:
        if list(vcf.header.samples) != samples:
            raise ValueError("Sample order changed")
        with pysam.VariantFile(final, "wz", header=vcf.header) as writer:
            for record in vcf:
                if record.id in ids:
                    raise ValueError("Duplicate representative ID")
                ids.add(record.id)
                ac, an, ns, nc = genotype_totals(record, dosage, carriers)
                record.info["AC"], record.info["AN"] = (ac,), an
                record.info["AF"], record.info["NS"] = ((ac/an if an else 0.0),), ns
                record.info["RC_N_CARRIERS"] = nc
                writer.write(record)
                count += 1
                gt50 += int(truvari.VariantRecord(record).var_size() > 50000)
    if dosage != partition["dosage"] or carriers != partition["carriers"]:
        raise ValueError("Per-sample ALT dosage or carrier-record counts changed unexpectedly")
    pysam.tabix_index(str(final), preset="vcf", force=True)
    return count, gt50, ids, totals_hash(dosage), totals_hash(carriers)


def trim_literal(record):
    """Minimal VCF anchoring, matching the TRUVARI chat; no left alignment."""
    ref, alt = record.ref, record.alts[0]
    prefix = 0
    while prefix < min(len(ref), len(alt)) and ref[prefix].upper() == alt[prefix].upper():
        prefix += 1
    suffix = 0
    while suffix < min(len(ref), len(alt)) - prefix and ref[-suffix-1].upper() == alt[-suffix-1].upper():
        suffix += 1
    start = prefix - 1 if prefix else 0
    end_ref = len(ref) - suffix + (0 if prefix else 1)
    end_alt = len(alt) - suffix + (0 if prefix else 1)
    newref, newalt = ref[start:end_ref], alt[start:end_alt]
    if not newref or not newalt or min(len(newref), len(newalt)) != 1:
        raise ValueError("Expected a pure, anchored literal indel")
    pos = record.pos + start
    record.alleles = (newref, newalt)
    record.pos = pos
    record.stop = pos - 1 + len(newref)


def prepare_same_start(prepared, source_folder, out):
    """Reuse symbolic partitions and derive only the trimmed literal partitions."""
    folder = out / "partitions"
    folder.mkdir(parents=True, exist_ok=True)
    if folder.resolve() == source_folder.resolve():
        raise ValueError("Same-start results must use a different output directory from the baseline")
    key = dict(source=prepared["input"], hashes=prepared["hashes"], trim_version=1,
               samples_sha256=totals_hash(prepared["samples"]))
    receipt = folder / "manifest.json"
    if receipt.exists():
        cached = json.loads(receipt.read_text())
        if cached["input"] == key and files_valid(folder, cached["hashes"]):
            print("Reusing validated same-start partitions", flush=True)
            return cached
    metadata = copy.deepcopy(prepared["partitions"])
    hashes, audit = {}, []
    for slug, m in metadata.items():
        source = source_folder / m["filename"]
        final = folder / m["filename"]
        if m["mode"] == "sequence_resolved":
            with pysam.VariantFile(source) as vcf:
                header = vcf.header.copy()
                if "RC_PRETRIM_POS" not in header.info:
                    header.add_line('##INFO=<ID=RC_PRETRIM_POS,Number=1,Type=Integer,Description="POS before minimal allele trimming">')
                with pysam.VariantFile(final, "wz", header=header) as writer:
                    for record in vcf:
                        record.translate(header)
                        oldpos = record.pos
                        oldlengths = (len(record.ref), len(record.alts[0]))
                        record.info["RC_PRETRIM_POS"] = oldpos
                        trim_literal(record)
                        audit.append((record.id, record.chrom, oldpos, record.pos, *oldlengths,
                                      len(record.ref), len(record.alts[0])))
                        writer.write(record)
            sorted_path = folder / (slug + ".sorted.vcf.gz")
            run(["bcftools", "sort", "-Oz", "-o", sorted_path, final], folder / "prepare.log")
            sorted_path.replace(final)
            pysam.tabix_index(str(final), preset="vcf", force=True)
            # This verifies record decoding, order, and all sample dosage/carrier
            # totals after trimming before the adapter ever sees these alleles.
            check = folder / (slug + ".checked.vcf.gz")
            count, _, _, _, _ = validate_and_finalize(final, check, m, prepared["samples"])
            if count != m["records"]:
                raise ValueError("Trimming changed the number of records")
            check.replace(final)
            Path(str(check) + ".tbi").replace(Path(str(final) + ".tbi"))
        else:
            shutil.copy2(source, final)
            shutil.copy2(Path(str(source) + ".tbi"), Path(str(final) + ".tbi"))
        shutil.copy2(source_folder / m["bed"], folder / m["bed"])
        for path in (final, Path(str(final) + ".tbi"), folder / m["bed"]):
            hashes[path.name] = digest(path)
    shutil.copy2(source_folder / "source_records.tsv", folder / "source_records.tsv")
    with (folder / "trimmed_positions.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["audit_ID", "CHROM", "original_POS", "trimmed_POS", "original_REF_length",
                         "original_ALT_length", "trimmed_REF_length", "trimmed_ALT_length"])
        writer.writerows(audit)
    for name in ("source_records.tsv", "trimmed_positions.tsv"):
        hashes[name] = digest(folder / name)
    if not files_valid(source_folder, prepared["hashes"]):
        raise ValueError("Baseline partitions changed during preparation")
    saved = dict(input=key, samples=prepared["samples"], partitions=metadata, hashes=hashes,
                 skipped=prepared["skipped"])
    write_json(receipt, saved)
    print(f"Prepared exact-position trial; trimmed {sum(r[2] != r[3] or r[4:6] != r[6:8] for r in audit):,} literal records", flush=True)
    return saved


def collapse_job(task):
    root, slug, m, samples, threshold, refdist, same_start = task
    kind = m["svtype"].split(":")[0]
    mode = m["mode"]
    passthrough = mode == "unresolved_passthrough" or (kind == "INS" and mode != "sequence_resolved")
    pctseq = threshold if kind == "INS" and not passthrough else 0.0
    pctovl = threshold if kind == "DEL" and not passthrough else 0.0
    if kind == "INS" and mode == "symbolic_size_position":
        mode = "symbolic_passthrough"
    elif kind == "DEL" and not passthrough:
        mode = "symbolic_size_overlap" if mode == "symbolic_size_position" else "literal_size_overlap"
    # Keep previous size-only DEL/symbolic-INS results as a baseline. Input
    # partitions are identical and remain reusable without reading all GTs again.
    output_slug = slug.removesuffix(m["mode"]) + mode
    job = root / f"pct{int(threshold*100)}" / output_slug
    job.mkdir(parents=True, exist_ok=True)
    source = root / "partitions" / m["filename"]
    bed = root / "partitions" / m["bed"]
    # Truvari 5.5 accepts -1 to disable the default 50 kb upper bound.
    sizemax = -1
    recipe = dict(schema=SCHEMA, input_sha256=digest(source), samples_sha256=totals_hash(samples),
                  pctsize=threshold, pctseq=pctseq, refdist=refdist, sizemin=10001,
                  sizemax=sizemax, gt="all", keep="common", chain=False,
                  typeignore=False, pctovl=pctovl, truvari=truvari.__version__, mode=mode,
                  bed_sha256=digest(bed))
    adapter = Path(__file__).with_name("collapse_same_start.py")
    if same_start:
        recipe.update(same_start=True, adapter_sha256=digest(adapter),
                      position_policy="identical CHROM/POS after minimal literal trimming; symbolic POS unchanged")
    receipt = job / "complete.json"
    if receipt.exists():
        saved = json.loads(receipt.read_text())
        if saved["recipe"] == recipe and files_valid(job, saved["hashes"]):
            print(f"Reused {job.parent.name}/{output_slug}", flush=True)
            return dict(saved["summary"], pctovl=pctovl, action="passthrough" if passthrough else "collapse")
    sequence_rule = f">={pctseq:g}" if pctseq else "off"
    overlap_rule = f">={pctovl:g}" if pctovl else "off"
    criteria = "PASSTHROUGH (no collapse)" if passthrough else (
        f"size>={threshold:g}; sequence={sequence_rule}; reciprocal_overlap={overlap_rule}; "
        f"refdist={refdist}; gt=all; keep=common; size max={sizemax}")
    if same_start:
        criteria += "; exact same CHROM/POS required"
    print(f"Running {job.parent.name}/{output_slug}: {m['records']:,} records; {criteria}", flush=True)
    raw, removed = job / "retained.raw.vcf", job / "removed.vcf"
    command = [sys.executable, "-m", "truvari", "collapse", "-i", str(source),
               "-o", str(raw), "-c", str(removed), "--bed", str(bed),
               "--gt", "all", "--keep", "common", "--pctsize", str(threshold),
               "--pctseq", str(pctseq), "--refdist", str(refdist),
               "--sizemin", "10001", "--sizemax", str(sizemax), "--pctovl", str(pctovl)]
    if same_start:
        command = [sys.executable, str(adapter), *command[4:]]
    if passthrough:
        raw = source
        with pysam.VariantFile(source) as vcf:
            removed.write_text(str(vcf.header))
    else:
        run(command, job / "truvari.log")
    final = job / "retained.vcf.gz"
    count, gt50, ids, dosage_hash, carrier_hash = validate_and_finalize(raw, final, m, samples)
    removed_ids = set()
    removed_count = removed_gt50 = 0
    with pysam.VariantFile(removed) as vcf:
        for record in vcf:
            removed_count += 1
            removed_gt50 += int(truvari.VariantRecord(record).var_size() > 50000)
            removed_ids.add(record.id)
    with pysam.VariantFile(source) as vcf:
        source_ids = {record.id for record in vcf}
    if count + removed_count != m["records"] or ids & removed_ids or ids | removed_ids != source_ids:
        raise ValueError("Retained/removed records do not partition input IDs")
    removed_gz = job / "removed.vcf.gz"
    pysam.tabix_compress(str(removed), str(removed_gz), force=True)
    # Removed order follows clusters, so no tabix index is promised for this file.
    summary = dict(threshold=threshold, svtype=m["svtype"], mode=mode, input_records=m["records"],
                   retained_records=count, collapsed_records=removed_count,
                   input_gt50kb=m["gt50kb"], retained_gt50kb=gt50, collapsed_gt50kb=removed_gt50,
                   pctsize=threshold, pctseq=pctseq, pctovl=pctovl, refdist=refdist, sizemax=sizemax,
                   action="passthrough" if passthrough else "collapse",
                   per_sample_ALT_dosage_preserved=True, per_sample_carrier_counts_preserved=True,
                   dosage_sha256=dosage_hash, carriers_sha256=carrier_hash)
    hashes = {f.name: digest(f) for f in (final, Path(str(final)+".tbi"), removed_gz)}
    write_json(receipt, dict(recipe=recipe, command=command if not passthrough else None,
                            summary=summary, hashes=hashes))
    print(f"Completed {job.parent.name}/{output_slug}: {m['records']:,} -> {count:,}", flush=True)
    return summary


def main():
    if truvari.__version__ != "5.5.0":
        raise RuntimeError("This workflow requires Truvari 5.5.0; run build_truvari55.py first")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bcf", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--refdist", type=int, default=500)
    parser.add_argument("--same-start", action="store_true", help="Only collapse identical CHROM/POS after minimal literal trimming")
    parser.add_argument("--reuse-prepared", type=Path, help="Baseline output directory whose verified partitions can be reused for --same-start")
    args = parser.parse_args()
    if args.workers < 1 or args.refdist < 0:
        parser.error("Invalid worker count or reference distance")
    if args.reuse_prepared and not args.same_start:
        parser.error("--reuse-prepared requires --same-start")
    args.bcf, args.out = args.bcf.resolve(), args.out.resolve()
    if args.reuse_prepared and args.reuse_prepared.resolve() == args.out:
        parser.error("Use a separate --out for the same-start trial")
    args.out.mkdir(parents=True, exist_ok=True)
    manifest = dict(utc=datetime.now(timezone.utc).isoformat(), parameters={k: str(v) for k, v in vars(args).items()},
                    script_sha256=digest(Path(__file__)), versions={p: importlib.metadata.version(p) for p in ("truvari", "pysam", "numpy", "pandas")},
                    bcftools=subprocess.check_output(["bcftools", "--version"], text=True).splitlines()[0],
                    thresholds=[0.9, 0.7], random_seed=None, randomness="No random operations",
                    position_policy="same CHROM/POS after minimal trimming; no left alignment" if args.same_start else "native proximity matching",
                    deletion="Size AND reciprocal overlap at each threshold; pctseq=0; representations stay separate",
                    insertion="Only literal unambiguous pure INS collapse: size AND sequence at each threshold; symbolic/ambiguous INS pass through",
                    other_symbolic="DUP/INV retain separate size/position matching; pctseq=0, pctovl=0",
                    genotype="gt all, keep common; AC/AN/AF/NS and RC_N_CARRIERS recalculated; missing alleles not imputed",
                    provenance="RC IDs map to original input; old collapse INFO cleared, source discovery counts retained as RC_SOURCE_DISCOVERY",
                    scope="Each threshold starts from identical partitions; no cross-type or cross-representation merges; symbolic INS and unresolved calls pass through",
                    note="Already-collapsed representatives cannot reconstruct every original call or cluster")
    write_json(args.out / "run_manifest.json", manifest)
    if args.same_start:
        base = args.reuse_prepared.resolve() if args.reuse_prepared else args.out / "untrimmed"
        prepared = prepare(args.bcf, base)
        prepared = prepare_same_start(prepared, base / "partitions", args.out)
    else:
        prepared = prepare(args.bcf, args.out)
    tasks = [(args.out, slug, m, prepared["samples"], threshold, args.refdist, args.same_start)
             for threshold in (0.9, 0.7) for slug, m in sorted(prepared["partitions"].items())]
    results = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        results = list(pool.map(collapse_job, tasks))
    if not results:
        raise ValueError("No eligible non-BND records")
    summary_path = args.out / "comparison.tsv"
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(results)
    if args.reuse_prepared and (args.reuse_prepared / "comparison.tsv").exists():
        def row_key(row):
            return float(row["threshold"]), row["svtype"], row["mode"]
        with (args.reuse_prepared / "comparison.tsv").open() as handle:
            previous = {row_key(r): r for r in csv.DictReader(handle, delimiter="\t")}
        compared = []
        for row in results:
            old = previous.get(row_key(row))
            if old is None or int(old["input_records"]) != row["input_records"]:
                continue
            compared.append(dict(threshold=row["threshold"], svtype=row["svtype"], mode=row["mode"],
                input_records=row["input_records"], baseline_retained=int(old["retained_records"]),
                same_start_retained=row["retained_records"],
                additional_retained=row["retained_records"] - int(old["retained_records"])))
        if compared:
            with (args.out / "comparison_vs_baseline.tsv").open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(compared[0]), delimiter="\t")
                writer.writeheader()
                writer.writerows(compared)
            print("Saved comparison_vs_baseline.tsv for matching baseline partitions", flush=True)
    print("\nthreshold\tSVTYPE\tmode\tinput\tretained\tcollapsed\tcollapsed_gt50kb")
    for r in results:
        print("\t".join(str(r[k]) for k in ("threshold", "svtype", "mode", "input_records", "retained_records", "collapsed_records", "collapsed_gt50kb")))
    print(f"\nSaved independent per-type experiments and comparison.tsv to {args.out}")


if __name__ == "__main__":
    main()
