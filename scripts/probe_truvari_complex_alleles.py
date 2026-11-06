#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = ["truvari==5.5.0"]
# ///
"""Synthetic, deterministic probes of native complex-allele matching (no panel data)."""
import argparse
import hashlib
import inspect
import json
from pathlib import Path
import subprocess
import tempfile

import pysam
import truvari
import truvari.comparisons


def run(output):
    assert truvari.__version__ == "5.5.0"
    cases = [
        ("simple_DEL_different_deleted_sequence", "A" + "C" * 60, "A", "A" + "G" * 60, "A"),
        ("complex_DEL_same_REF_different_ALT", "A" * 80, "A" + "C" * 19, "A" * 80, "A" + "G" * 19),
        ("complex_DEL_different_net_size", "A" * 80, "A" + "C" * 19, "A" * 80, "A" + "G" * 9),
        ("simple_INS_different_inserted_sequence", "A", "A" + "C" * 100, "A", "A" + "G" * 100),
        ("padded_complex_INS", "A" * 1000, "A" * 999 + "C" * 101,
         "A" * 1000, "A" * 999 + "G" * 101),
        ("same_complex_INS_shared_prefix_trimmed", "A", "C" * 101, "A", "G" * 101),
        ("equal_length_replacement", "A" * 80, "C" * 80, "A" * 80, "G" * 80),
    ]
    results = []
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        for name, ref1, alt1, ref2, alt2 in cases:
            header = pysam.VariantHeader()
            header.contigs.add("chr1", length=100000)
            header.formats.add("GT", 1, "String", "Genotype")
            for sample in ("s1", "s2"):
                header.add_sample(sample)
            source = root / f"{name}.vcf.gz"
            # Shift both trimmed records to preserve their genomic haplotypes.
            starts = [1999, 1999] if name == "same_complex_INS_shared_prefix_trimmed" else [1000, 1000]
            # Distinct deleted reference sequences occupy distinct loci in this control.
            if name == "simple_DEL_different_deleted_sequence":
                starts = [1000, 1200]
            with pysam.VariantFile(source, "wz", header=header) as handle:
                for i, (ref, alt) in enumerate(((ref1, alt1), (ref2, alt2))):
                    record = handle.new_record(contig="chr1", start=starts[i], alleles=(ref, alt), id=f"v{i}")
                    for j, sample in enumerate(("s1", "s2")):
                        record.samples[sample]["GT"] = (1, 0) if i == j else (0, 0)
                        record.samples[sample].phased = True
                    handle.write(record)
            pysam.tabix_index(str(source), preset="vcf", force=True)
            with truvari.VariantFile(str(source)) as handle:
                a, b = list(handle)
                row = {"case": name, "type": a.var_type().name,
                       "start0": starts,
                       "ref_lengths": [len(ref1), len(ref2)], "alt_lengths": [len(alt1), len(alt2)],
                       "native_sizes": [a.var_size(), b.var_size()],
                       "sequence_similarity": a.seqsim(b), "size_similarity": a.sizesim(b)[0],
                       "reciprocal_overlap": a.recovl(b)}
            for threshold in (.7, .9):
                kept = root / f"{name}.{threshold}.kept.vcf.gz"
                removed = root / f"{name}.{threshold}.removed.vcf.gz"
                cmd = ["truvari", "collapse", "-i", str(source), "-o", str(kept), "-c", str(removed),
                       "-p", str(threshold), "-P", str(threshold), "-r", "500", "-O", "0",
                       "-s", "0", "-S", "-1", "--keep", "common", "--gt", "all"]
                proc = subprocess.run(cmd, text=True, capture_output=True)
                if proc.returncode:
                    raise RuntimeError(proc.stderr)
                # Index to silence reader warnings and validate ordering.
                pysam.tabix_index(str(kept), preset="vcf", force=True)
                pysam.tabix_index(str(removed), preset="vcf", force=True)
                with pysam.VariantFile(kept) as handle:
                    count = sum(1 for _ in handle)
                with pysam.VariantFile(removed) as handle:
                    removed_count = sum(1 for _ in handle)
                assert count + removed_count == 2
                row[f"retained_p{round(threshold * 100)}"] = count
            results.append(row)
    by_name = {r["case"]: r for r in results}
    assert by_name["complex_DEL_same_REF_different_ALT"]["sequence_similarity"] == 1
    assert by_name["complex_DEL_same_REF_different_ALT"]["retained_p90"] == 1
    assert by_name["complex_DEL_different_net_size"]["retained_p70"] == 1
    assert by_name["complex_DEL_different_net_size"]["retained_p90"] == 2
    assert by_name["padded_complex_INS"]["retained_p90"] == 1
    assert by_name["same_complex_INS_shared_prefix_trimmed"]["retained_p70"] == 2
    provenance = {"truvari": truvari.__version__, "seed": 1729, "data": "synthetic only; no RNG",
                  "source_sha256": {Path(inspect.getfile(obj)).name:
                      hashlib.sha256(Path(inspect.getfile(obj)).read_bytes()).hexdigest()
                      for obj in (truvari.VariantRecord, truvari.comparisons)}}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps({"provenance": provenance, "results": results}, indent=2) + "\n")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args().output)
