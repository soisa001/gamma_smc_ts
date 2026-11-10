# /// script
# requires-python = ">=3.11"
# dependencies = ["rapidfuzz==3.14.1"]
# ///
"""Read-only same-position DEL/INS merge screen. Never rewrites a VCF.

Length similarity is min/max. Sequence similarity is 1-global_edit_distance/max
on the inserted/deleted sequence after removing common REF/ALT prefix/suffix.
This explicit payload metric is not Truvari's full haplotype comparison.
"""

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
import csv
from datetime import datetime, timezone
import hashlib
from itertools import combinations
import json
import math
from pathlib import Path
import re
import subprocess

import rapidfuzz
from rapidfuzz.distance import Levenshtein


def sha256(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def query_lines(command, stderr_path):
    with stderr_path.open("a", encoding="utf-8") as errors:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=errors,
                                   text=True, encoding="utf-8")
        try:
            yield from process.stdout
            if process.wait():
                raise RuntimeError(f"bcftools failed; see {stderr_path}")
        finally:
            process.stdout.close()
            if process.poll() is None:
                process.terminate()
            process.wait()


def payload(ref, alt, kind):
    """Return sequence, trimmed start offset, and why sequence is unavailable."""
    ref, alt = ref.upper(), alt.upper()
    if not re.fullmatch("[ACGTN]+", ref) or not re.fullmatch("[ACGTN]+", alt):
        return None, None, "symbolic_or_nonliteral_allele"
    prefix = 0
    while prefix < min(len(ref), len(alt)) and ref[prefix] == alt[prefix]:
        prefix += 1
    ref, alt = ref[prefix:], alt[prefix:]
    suffix = 0
    while suffix < min(len(ref), len(alt)) and ref[-suffix-1] == alt[-suffix-1]:
        suffix += 1
    if suffix:
        ref, alt = ref[:-suffix], alt[:-suffix]
    seq = alt if kind == "INS" and not ref else ref if kind == "DEL" and not alt else None
    if not seq:
        return None, prefix, "not_a_pure_insertion_or_deletion"
    if "N" in seq:
        return None, prefix, "ambiguous_sequence"
    return seq, prefix, "available"


def parse_record(line, row, nsamples):
    fields = line.rstrip("\n").split("\t")
    if len(fields) != 8 + nsamples:
        raise ValueError(f"Unexpected query field count at selected row {row}")
    chrom, pos, ident, ref, alt, kind, svlen, end = fields[:8]
    if "," in alt:
        raise ValueError("This audit expects the previously verified biallelic records")
    if kind == "." and alt.startswith("<") and alt.endswith(">"):
        kind = alt[1:-1]
    kind = kind.split(":")[0].upper()
    seq, offset, reason = payload(ref, alt, kind)
    length = abs(int(svlen)) if svlen not in (".", "") and "," not in svlen else None
    if not length:
        length = len(seq) if seq else None
    if not length and kind == "DEL" and end not in (".", ""):
        span = int(end) - int(pos)
        length = span if span > 0 else None
    calls, carriers, partial = {}, set(), 0
    for index, gt in enumerate(fields[8:]):
        alleles = re.split(r"[/|]", gt)
        if any(a not in ("0", "1", ".") for a in alleles):
            raise ValueError(f"Invalid biallelic GT at selected row {row}")
        known = [int(a) for a in alleles if a != "."]
        if 1 in known:
            carriers.add(index)
        if len(known) == len(alleles):
            calls[index] = (len(known), sum(known))
        elif known:
            partial += 1
    return dict(row=row, chrom=chrom, pos=int(pos), id=ident, kind=kind,
                length=length, seq=seq, offset=offset, sequence_status=reason,
                length_conflict=bool(seq and length != len(seq)),
                calls=calls, carriers=carriers, partial=partial)


def compare(task):
    a, b, size_threshold, seq_threshold, max_bp = task
    overlap = a["carriers"] & b["carriers"]
    common = a["calls"].keys() & b["calls"].keys()
    conflicts = sum(a["calls"][s] != b["calls"][s] for s in common)
    result = dict(chrom=a["chrom"], pos=a["pos"], row_a=a["row"], row_b=b["row"],
                  id_a=a["id"], id_b=b["id"], type_a=a["kind"], type_b=b["kind"],
                  length_a=a["length"], length_b=b["length"], length_similarity=None,
                  sequence_status_a=a["sequence_status"], sequence_status_b=b["sequence_status"],
                  sequence_similarity=None, sequence_below_threshold=False,
                  carriers_a=len(a["carriers"]), carriers_b=len(b["carriers"]),
                  shared_carriers=len(overlap), conflicting_called_genotypes=conflicts,
                  partial_genotypes=a["partial"]+b["partial"],
                  length_annotation_conflict=a["length_conflict"] or b["length_conflict"],
                  similarity_pass=False, status="")
    if a["kind"] != b["kind"]:
        status = "different_type"
    elif a["kind"] not in ("DEL", "INS"):
        status = "unsupported_type"
    elif not a["length"] or not b["length"]:
        status = "length_unavailable"
    else:
        ratio = min(a["length"], b["length"]) / max(a["length"], b["length"])
        result["length_similarity"] = ratio
        if ratio + 1e-12 < size_threshold:
            status = "length_below_threshold"
        elif not a["seq"] or not b["seq"]:
            status = "sequence_unavailable"
        elif result["length_annotation_conflict"] or a["offset"] != b["offset"]:
            status = "representation_review"
        elif max(len(a["seq"]), len(b["seq"])) > max_bp:
            status = "sequence_length_limit"
        else:
            longest = max(len(a["seq"]), len(b["seq"]))
            cutoff = math.floor((1-seq_threshold) * longest + 1e-9)
            distance = Levenshtein.distance(a["seq"], b["seq"], score_cutoff=cutoff)
            if distance > cutoff:
                result["sequence_below_threshold"] = True
                status = "sequence_below_threshold"
            else:
                result["sequence_similarity"] = 1-distance/longest
                result["similarity_pass"] = True
                review = overlap or conflicts or result["partial_genotypes"] or not a["carriers"] or not b["carriers"]
                status = "candidate_genotype_review" if review else "candidate_disjoint_carriers"
    result["status"] = status
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bcf", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--length-similarity", type=float, default=0.7)
    parser.add_argument("--sequence-similarity", type=float, default=0.7)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--max-sequence-bp", type=int, default=1000000)
    args = parser.parse_args()
    if not (0 <= args.length_similarity <= 1 and 0 <= args.sequence_similarity <= 1):
        parser.error("Similarity thresholds must be between 0 and 1")
    if args.workers < 1 or args.max_sequence_bp < 1:
        parser.error("Workers and maximum sequence length must be positive")
    args.bcf = args.bcf.resolve()
    args.out.mkdir(parents=True, exist_ok=True)
    args.out = args.out.resolve()
    stat = args.bcf.stat()
    errors = args.out / "bcftools.stderr.log"
    errors.write_text("", encoding="utf-8")
    commands = []

    def run(command):
        commands.append(command)
        return query_lines(command, errors)

    print("Phase 1: scan positions and identify repeated starts.", flush=True)
    samples = list(run(["bcftools", "query", "-l", str(args.bcf)]))
    if not samples or len(set(samples)) != len(samples):
        raise ValueError("Empty or duplicate sample manifest")
    counts = Counter()
    for line in run(["bcftools", "query", "-f", "%CHROM\t%POS\n", str(args.bcf)]):
        chrom, pos = line.rstrip("\n").split("\t")
        counts[chrom, int(pos)] += 1
    repeated = {key: n for key, n in counts.items() if n > 1}
    targets = args.out / "shared_positions.tsv"
    targets.write_text("".join(f"{c}\t{p}\n" for c, p in sorted(repeated)), encoding="utf-8")
    print(f"{sum(counts.values()):,} records; {len(repeated):,} repeated positions; "
          f"{sum(repeated.values()):,} records to inspect.", flush=True)
    records = defaultdict(list)
    print("Phase 2: read alleles and genotypes at repeated starts.", flush=True)
    if repeated:
        fmt = "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/END[\t%GT]\n"
        command = ["bcftools", "query", "-u", "-T", str(targets), "--targets-overlap", "0", "-f", fmt, str(args.bcf)]
        for row, line in enumerate(run(command), 1):
            rec = parse_record(line, row, len(samples))
            records[rec["chrom"], rec["pos"]].append(rec)
    if {key: len(value) for key, value in records.items()} != repeated:
        raise ValueError("Selected records do not match the position scan")
    current_stat = args.bcf.stat()
    if (stat.st_size, stat.st_mtime_ns) != (current_stat.st_size, current_stat.st_mtime_ns):
        raise ValueError("BCF changed during inspection")
    tasks = [(a, b, args.length_similarity, args.sequence_similarity, args.max_sequence_bp)
             for key in sorted(records) for a, b in combinations(records[key], 2)]
    print(f"Phase 3: compare {len(tasks):,} pairs with {args.workers} workers.", flush=True)
    results = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for i, result in enumerate(pool.map(compare, tasks, chunksize=1), 1):
            results.append(result)
            if i % 100 == 0:
                print(f"Compared {i:,}/{len(tasks):,} pairs.", flush=True)
    pair_path = args.out / "pair_audit.tsv"
    with pair_path.open("w", encoding="utf-8", newline="") as handle:
        if results:
            writer = csv.DictWriter(handle, fieldnames=list(results[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(results)
        else:
            handle.write("chrom\tpos\tstatus\n")
    per_locus = defaultdict(list)
    for result in results:
        per_locus[result["chrom"], result["pos"]].append(result)
    all_pass = sum(all(r["similarity_pass"] for r in group) for group in per_locus.values())
    disjoint_groups = [key for key, group in per_locus.items()
                       if all(r["status"] == "candidate_disjoint_carriers" for r in group)]
    summary = {
        "utc": datetime.now(timezone.utc).isoformat(),
        "parameters": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
        "input_size_bytes": stat.st_size, "input_mtime_ns": stat.st_mtime_ns,
        "script_sha256": sha256(Path(__file__)), "targets_sha256": sha256(targets),
        "samples_sha256": hashlib.sha256("".join(samples).encode()).hexdigest(),
        "samples": len(samples), "records": sum(counts.values()),
        "shared_positions": len(repeated), "records_at_shared_positions": sum(repeated.values()),
        "sequence_availability_records": dict(Counter(r["sequence_status"] for group in records.values() for r in group)),
        "pairs": len(results), "pair_status_counts": dict(Counter(r["status"] for r in results)),
        "pairs_passing_both_similarities": sum(r["similarity_pass"] for r in results),
        "loci_with_any_similarity_pass": sum(any(r["similarity_pass"] for r in group) for group in per_locus.values()),
        "loci_with_all_pairs_passing_similarity": all_pass,
        "loci_with_all_pairs_passing_and_disjoint_carriers": len(disjoint_groups),
        "possible_reduction_in_those_disjoint_groups": sum(len(records[key])-1 for key in disjoint_groups),
        "length_metric": "min(abs(SVLEN))/max(abs(SVLEN)); sequence-length or DEL END-POS fallback when absent",
        "sequence_metric": "1 - global Levenshtein distance/max(payload lengths); same reference-forward orientation; trimmed pure DEL/INS payload only",
        "below_threshold_scores": "Bounded alignment: exact score only for passing pairs; failed score is blank, not zero",
        "genotypes": "Full calls compared by ploidy and ALT dosage; partial calls and shared carriers trigger review; phase not resolved",
        "scope": "Identical CHROM/POS only; symbolic or ambiguous sequence is unavailable, not a sequence match",
        "interpretation": "Candidate screen only; no transitive clustering, no genotype reconciliation, no VCF modification",
        "commands": commands, "rapidfuzz_version": rapidfuzz.__version__,
        "bcftools_version": subprocess.check_output(["bcftools", "--version"], text=True).splitlines()[0],
        "randomness": "None; deterministic ordering",
    }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: summary[k] for k in ("shared_positions", "records_at_shared_positions", "sequence_availability_records", "pairs", "pair_status_counts", "pairs_passing_both_similarities", "loci_with_all_pairs_passing_similarity", "loci_with_all_pairs_passing_and_disjoint_carriers", "possible_reduction_in_those_disjoint_groups")}, indent=2), flush=True)
    print(f"Read-only audit complete: {args.out}", flush=True)


if __name__ == "__main__":
    main()
