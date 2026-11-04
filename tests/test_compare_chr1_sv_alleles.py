from __future__ import annotations

import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import compare_chr1_sv_alleles as compare


def sha(value):
    return hashlib.sha256(value.encode()).hexdigest()


def record(ident, pos, ref, alt, status="eligible", kind=None):
    info = {"SVLEN": "500"} if status == "symbolic" else {}
    result = compare.pilot.allele(ref, alt, pos, info, 1, -1)
    derived, size, start, end, _ = result
    return {"id": ident, "source_row": ident, "source_id": ident, "chrom": "chr1", "pos": pos,
            "type": kind or derived, "size": size, "start0": start, "end0": end, "status": status,
            "ref_sha256": sha(ref), "alt_sha256": sha(alt)}, (ref, alt)


def fixture(root, native=False):
    (root / "prepare").mkdir()
    (root / "merge").mkdir()
    definitions = [
        record("d1", 101, "A" + "T" * 100, "A"),
        record("d2", 151, "A" + "T" * 100, "A"),
        record("d3", 101, "A" + "T" * 100, "A"),
        record("d4", 301, "A" + "T" * 101, "A" + "C" * 50, "complex_replacement"),
        record("ds", 401, "N", "<DEL>", "symbolic"),
        record("i1", 501, "A", "A" + "C" * 100),
        record("i2", 501, "A", "A" + "C" * 100),
        record("i3", 501, "A", "A" + "G" * 100),
        record("i4", 601, "A", "A" + "C" * 100),
        record("i5", 701, "AT", "A" + "C" * 51, "complex_replacement"),
    ]
    rows = {r["id"]: r for r, _ in definitions}
    if native:
        for row in rows.values():
            row["merge_submitted"] = 1
    seq = {r["id"]: s for r, s in definitions}
    source = root / "prepare" / "source_alleles.tsv.gz"
    samples = root / "prepare" / "samples.txt"
    samples.write_text("s1\ns2\n")
    compare.pilot.tsv(source, list(rows.values()), list(next(iter(rows.values()))))
    compare.pilot.manifest(root / "prepare", {"chrom": "chr1", "min_size": 50, "max_size": 100000}, [source, samples])
    parts = {("DEL", "kept"): [("d1", "CollapseId=1.0;NumCollapsed=2")],
             ("DEL", "removed"): [("d3", "MatchId=1.0"), ("d2", "MatchId=1.0")],
             ("INS", "kept"): [("i1", "CollapseId=1.0;NumCollapsed=1"), ("i3", "."), ("i4", ".")],
             ("INS", "removed"): [("i2", "MatchId=1.0")]}
    if native:
        parts[("DEL", "kept")].extend([("d4", "."), ("ds", ".")])
        parts[("INS", "kept")].append(("i5", "."))
    for (kind, role), ids in parts.items():
        path = root / "merge" / f"{kind}.{role}.vcf.gz"
        with gzip.open(path, "wt") as handle:
            handle.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n")
            for ident, info in ids:
                ref, alt = seq[ident]
                # Deliberately conflicting SVLEN proves summary uses sequence lengths.
                info = (info + ";" if info != "." else "") + "SVLEN=99999"
                handle.write(f"chr1\t{rows[ident]['pos']}\t{ident}\t{ref}\t{alt}\t.\tPASS\t{info}\tGT\t0/1\t0/0\n")
        Path(str(path) + ".tbi").write_bytes(b"synthetic index presence marker")
    return rows


class MetricsTests(unittest.TestCase):
    def test_native_comparison_includes_complex_and_symbolic_members(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture(root, native=True)
            output = root / "compare_ins_del"
            output.mkdir()
            compare.run(SimpleNamespace(root=root, output=output, threads=4))
            report = json.loads((output / "summary.json").read_text())
            rows = {(r["scope"], r["type"], r["callset"]): r for r in report["summary"]}
            for kind, after_count in (("DEL", 3), ("INS", 4)):
                for callset in ("before", "after"):
                    submitted = rows[("merge_eligible", kind, callset)]
                    all_rows = rows[("all_ins_del", kind, callset)]
                    self.assertEqual(submitted["sv_records"], 5 if callset == "before" else after_count)
                    self.assertEqual(submitted["raw_bp_sum_records"], all_rows["raw_bp_sum_records"])
            self.assertEqual(rows[("merge_eligible", "DEL", "after")]["symbolic_records_excluded_from_raw_bp"], 1)

    def test_deletion_union_handles_nested_touching_and_contigs(self):
        intervals = [("chr1", 100, 200), ("chr1", 150, 250), ("chr1", 170, 180),
                     ("chr1", 250, 300), ("chr2", 100, 200)]
        self.assertEqual(compare.union_intervals(intervals), [("chr1", 100, 300), ("chr2", 100, 200)])

    def test_unique_insertions_keep_distinct_sequence_and_position(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture(root)
            rows = compare.read_source(root / "prepare" / "source_alleles.tsv.gz")
            subset = [r for r in rows.values() if r["type"] == "INS" and r["status"] == "eligible"]
            summary, _ = compare.summarize(subset, "INS")
            self.assertEqual(summary["sv_records"], 4)
            self.assertEqual(summary["unique_alleles"], 3)
            self.assertEqual(summary["raw_bp_sum_records"], 400)
            self.assertEqual(summary["unique_inserted_bp"], 300)
            self.assertEqual(compare.spacing(subset), [0, 0, 100])
            self.assertEqual(compare.spacing(subset, unique=True), [0, 100])

    def test_full_comparison_without_combined_vcf_and_cache(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture(root)
            output = root / "compare_ins_del"
            output.mkdir()
            args = SimpleNamespace(root=root, output=output, threads=4)
            compare.run(args)
            report = json.loads((output / "summary.json").read_text())
            rows = {(r["scope"], r["type"], r["callset"]): r for r in report["summary"]}
            before = rows[("merge_eligible", "DEL", "before")]
            after = rows[("merge_eligible", "DEL", "after")]
            self.assertEqual((before["sv_records"], after["sv_records"]), (3, 1))
            self.assertEqual((before["unique_deleted_reference_bp"], after["unique_deleted_reference_bp"]), (150, 100))
            self.assertEqual(before["raw_bp_sum_records"], 300)
            self.assertEqual(before["raw_bp_sum_unique_alleles"], 200)
            full = rows[("all_ins_del", "DEL", "before")]
            self.assertEqual(full["sv_records"], 5)
            self.assertEqual(full["raw_bp_sum_records"], 351)
            self.assertEqual(full["symbolic_records_excluded_from_raw_bp"], 1)
            self.assertEqual(full["complex_del_alleles_excluded_from_union"], 1)
            self.assertEqual(full["unique_deleted_reference_bp"], 150)
            self.assertEqual(rows[("all_ins_del", "INS", "before")]["unique_inserted_bp"], 350)
            self.assertEqual(report["eligible_removed"], 3)
            self.assertFalse((root / "merge" / "new.vcf.gz").exists())
            mtime = (output / "before_after.tsv").stat().st_mtime_ns
            compare.run(args)
            self.assertEqual((output / "before_after.tsv").stat().st_mtime_ns, mtime)
            if os.environ.get("SV_ALLELE_QA_DIR"):
                shutil.copytree(output, os.environ["SV_ALLELE_QA_DIR"], dirs_exist_ok=True)

    def test_missing_membership_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture(root)
            source = compare.read_source(root / "prepare" / "source_alleles.tsv.gz")
            partitions = {(k, r): {} for k in ("INS", "DEL") for r in ("kept", "removed")}
            with self.assertRaisesRegex(ValueError, "disjoint"):
                compare.validate_membership(source, partitions)

    def test_truncated_gzip_and_changed_allele_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture(root)
            source = compare.read_source(root / "prepare" / "source_alleles.tsv.gz")
            sample_sha = sha("s1\ns2\n")
            path = root / "merge" / "INS.kept.vcf.gz"
            source["i1"]["raw_bp"] = 101
            with self.assertRaisesRegex(ValueError, "raw length changed"):
                compare.read_partition(path, "INS", source, sample_sha)
            source["i1"]["raw_bp"] = 100
            path.write_bytes(path.read_bytes()[:-8])
            with self.assertRaises(EOFError):
                compare.read_partition(path, "INS", source, sample_sha)

    def test_empty_type_and_no_cross_contig_distance(self):
        summary, intervals = compare.summarize([], "DEL")
        self.assertEqual(summary["unique_deleted_reference_bp"], 0)
        self.assertEqual(intervals, [])
        self.assertEqual(compare.spacing([]), [])
        self.assertEqual(compare.spacing([{"chrom": "chr1", "pos": 100}, {"chrom": "chr2", "pos": 200}]), [])


if __name__ == "__main__":
    unittest.main()
