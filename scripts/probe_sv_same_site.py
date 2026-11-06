#!/usr/bin/env python3
"""Synthetic feasibility probe, NOT a production panel merging pipeline.

Run with uv run --with truvari==5.5.0 python scripts/probe_sv_same_site.py.
Uses direct ALT edit similarity 1-distance/max(lengths), strictly > threshold.
All examples have one common POS/REF; no controlled data are read.
"""
from collections import Counter
import json
from pathlib import Path
import tempfile

import edlib
import pysam


def cluster(ref, alts, genotypes, threshold, keep_disjoint=True):
    """Proof of within-record clustering and allele-index mapping only."""
    # Trim common context from ALL alleles and REF together, not independently.
    sequences = [ref, *alts]
    while all(sequences) and len({s[0] for s in sequences}) == 1:
        sequences = [s[1:] for s in sequences]
    while all(sequences) and len({s[-1] for s in sequences}) == 1:
        sequences = [s[:-1] for s in sequences]
    seq = dict(enumerate(sequences[1:], 1))
    ac = Counter(a for gt in genotypes for a in gt if a is not None and a > 0)
    carriers = {a: {i for i, gt in enumerate(genotypes) if a in gt} for a in seq}
    reps, members, mapping = [], [], {0: 0, None: None}
    for a in sorted(seq, key=lambda a: (-ac[a], a)):
        for i, rep in enumerate(reps):
            if keep_disjoint and carriers[a] & members[i]:
                continue
            denominator = max(len(seq[a]), len(seq[rep]))
            similarity = 1 if not denominator else 1 - edlib.align(seq[a], seq[rep])["editDistance"] / denominator
            # Retain net DEL/INS classes and net-size compatibility.
            size_a, size_r = len(alts[a - 1]) - len(ref), len(alts[rep - 1]) - len(ref)
            if size_a * size_r <= 0:
                continue
            size_similarity = min(abs(size_a), abs(size_r)) / max(abs(size_a), abs(size_r))
            if similarity > threshold and size_similarity >= threshold:
                mapping[a] = i + 1
                members[i].update(carriers[a])
                break
        else:
            reps.append(a)
            members.append(set(carriers[a]))
            mapping[a] = len(reps)
    return reps, mapping


def main():
    ref = "A" + "G" * 200 + "G"
    alts = ("A" + "C" * 100 + "G", "A" + "C" * 75 + "T" * 25 + "G", "A" + "T" * 100 + "G")
    gts = [(1, 0), (0, 2), (1, 3), (None, 0)]
    reps, mapping = cluster(ref, alts, gts, .7)
    assert reps == [1, 3] and mapping[1] == mapping[2] == 1
    assert len(cluster(ref, alts, gts, .9)[0]) == 3
    assert len(cluster(ref, alts, gts + [(1, 2)], .7)[0]) == 3
    # Strict >0.7: exactly 0.7 must not merge.
    boundary = (alts[0], "A" + "C" * 70 + "T" * 30 + "G")
    assert len(cluster(ref, boundary, [(1, 0), (0, 2)], .7)[0]) == 2
    header = pysam.VariantHeader()
    header.contigs.add("chr1", length=1000000)
    header.formats.add("GT", 1, "String", "Genotype")
    header.formats.add("PS", 1, "Integer", "Phase set")
    for i in range(len(gts)):
        header.add_sample(f"synthetic_{i}")
    with tempfile.TemporaryDirectory() as temp:
        path = Path(temp) / "merged.vcf.gz"
        with pysam.VariantFile(path, "wz", header=header) as out:
            record = out.new_record(contig="chr1", start=999, alleles=(ref, *(alts[a-1] for a in reps)))
            for sample, gt in zip(record.samples.values(), gts):
                sample["GT"] = tuple(mapping[a] for a in gt)
                sample.phased = True
                sample["PS"] = 500
            out.write(record)
        pysam.tabix_index(str(path), preset="vcf", force=True)
        with pysam.VariantFile(path) as inp:
            result = next(inp)
            actual = [sample["GT"] for sample in result.samples.values()]
            assert actual == [(1, 0), (0, 1), (1, 2), (None, 0)]
            assert all(s.phased and s["PS"] == 500 for s in result.samples.values())
    report = {"synthetic_only": True, "seed": 1729, "randomness_used": False,
              "score": "1 - ALT edit distance / max ALT lengths after group-common flank trimming",
              "threshold": ">0.7", "alleles_before": 3, "alleles_after": 2,
              "p90_alleles_after": 3, "phase_and_missingness_roundtrip": True,
              "same_sample_merge_prevented": True, "exact_threshold_rejected": True,
              "production_pipeline_implemented": False}
    output = Path(__file__).resolve().parents[1] / "docs/sv_same_site_synthetic_results.json"
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
