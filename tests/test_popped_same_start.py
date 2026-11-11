"""Synthetic exact-position clustering and deletion-union regression tests."""
import csv
import gzip
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import run_popped_simple_sv as simple
import compare_popped_merge_modes as compare


class UnionTests(unittest.TestCase):
    def test_overlap_and_union_loss(self):
        before = compare.union([(100, 200), (100, 180), (105, 205), (300, 350)])
        after = compare.union([(100, 200), (300, 350)])
        self.assertEqual(compare.coverage(before), 155)
        self.assertEqual(compare.coverage(after), 150)
        self.assertEqual(compare.coverage(before) - compare.intersection_bp(before, after), 5)


@unittest.skipUnless(os.environ.get("SV_PILOT_INTEGRATION") == "1", "Requires Linux native tools")
class SameStartTests(unittest.TestCase):
    def test_run_compare_and_resume(self):
        bcftools = os.environ.get("BCFTOOLS_BIN", "bcftools")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raw = root / "input.vcf"
            header = ('##fileformat=VCFv4.2\n##contig=<ID=chr1>\n'
                      '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
                      '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="PS">\n'
                      '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\tc\n')
            gt = [("1|0:1", "0|0:1", "0|0:1"), ("0|0:1", "0|1:1", "0|0:1"), ("0|0:1", "0|0:1", "1|0:1")]
            definitions = [
                (1000, "d100", "A" + "C"*100, "A", 0),
                (1000, "d80", "A" + "C"*80, "A", 1),
                (1005, "d_shift", "A" + "C"*100, "A", 2),
                (2000, "i100", "A", "A" + "T"*100, 0),
                (2000, "i80", "A", "A" + "T"*80, 1),
                (2005, "i_shift", "A", "A" + "T"*100, 2),
                (3000, "padded", "GGA", "GGA" + "T"*60, 0),
                (3002, "already_trimmed", "A", "A" + "T"*60, 1),
                (4000, "replacement", "A", "T"*101, 0),
                (5000, "small49", "A", "A" + "G"*49, 0),
                (5000, "large50", "A", "A" + "G"*50, 1),
                (6000, "right_del", "T"*60 + "A", "A", 0),
            ]
            raw.write_text(header + "".join("\t".join(map(str, ("chr1", pos, name, ref, alt, ".", "PASS", ".", "GT:PS", *gt[sample]))) + "\n"
                                           for pos, name, ref, alt, sample in definitions))
            compressed = root / "input.vcf.gz"
            subprocess.run([bcftools, "view", "-Oz", "-o", str(compressed), str(raw)], check=True, capture_output=True)
            subprocess.run([bcftools, "index", "-t", str(compressed)], check=True, capture_output=True)
            current, same = root / "current", root / "same"
            simple.main(["--root", str(current), "--input-bcf", str(compressed), "--bcftools", bcftools])
            baseline_times = {p: p.stat().st_mtime_ns for p in current.rglob("*") if p.is_file()}
            args = ["--root", str(same), "--reuse-prepared", str(current), "--same-start", "--phase", "merge", "--bcftools", bcftools]
            try:
                simple.main(args)
            except BaseException:
                for path in (same / "logs").glob("*.log"):
                    print(path.name, path.read_text()[-5000:])
                raise
            self.assertTrue((same / "prepare").is_symlink())
            self.assertTrue(all(p.stat().st_mtime_ns == stamp for p, stamp in baseline_times.items()))
            import pysam
            for tag in ("p70", "p90"):
                for kind in ("INS", "DEL"):
                    base = same / tag / kind
                    with pysam.VariantFile(base / "retained.vcf.gz") as handle:
                        kept = list(handle)
                    with pysam.VariantFile(base / "collapsed.vcf.gz") as handle:
                        removed = list(handle)
                    self.assertTrue({"i_shift" if kind == "INS" else "d_shift"} <= {r.id for r in kept})
                    # Validate native links independently of the adapter assertion.
                    reps = {str(r.info["CollapseId"]): r for r in kept if "CollapseId" in r.info}
                    for record in removed:
                        mid = record.info["MatchId"]
                        cluster = mid[0] if isinstance(mid, tuple) else mid.split(",")[0]
                        self.assertEqual(record.pos, reps[cluster].pos)
                    if kind == "DEL":
                        self.assertEqual(len(kept), 3 if tag == "p70" else 4)
                    if kind == "INS":
                        merged = {r.id for r in removed}
                        self.assertTrue(bool({"padded", "already_trimmed"} & merged))
                        self.assertTrue(bool({"small49", "large50"} & merged))
                with pysam.VariantFile(same / tag / "chr1.sv50.vcf.gz") as handle:
                    replacements = [str(r) for r in handle if r.id == "replacement"]
                self.assertEqual(len(replacements), 1)
                self.assertIn("1|0:1\t0|0:1\t0|0:1", replacements[0])
            report = root / "report"
            cmd = ["--current-root", str(current), "--same-start-root", str(same), "--root", str(report), "--bcftools", bcftools]
            compare.main(cmd)
            with open(report / "comparison.tsv") as handle:
                rows = {(r["policy"], r["similarity"], r["type"]): r for r in csv.DictReader(handle, delimiter="\t")}
            # Enumerate deleted bases independently, including right-anchored DEL.
            before = set(range(1000, 1105)) | set(range(5999, 6059))
            for tag, threshold in (("p70", "0.7"), ("p90", "0.9")):
                for policy, directory in (("current_r500", current), ("same_start", same)):
                    covered = set()
                    with pysam.VariantFile(directory / tag / "chr1.sv50.vcf.gz") as handle:
                        for r in handle:
                            if len(r.ref) > len(r.alts[0]):
                                start = r.pos-1 if r.id == "right_del" else r.pos
                                covered.update(range(start, start + len(r.ref)-len(r.alts[0])))
                    row = rows[policy, threshold, "DEL"]
                    self.assertEqual(int(row["before_unique_bp"]), len(before))
                    self.assertEqual(int(row["after_unique_bp"]), len(covered))
                    self.assertEqual(int(row["unique_bp_removed"]), len(before-covered))
                self.assertGreater(int(rows["same_start", threshold, "DEL"]["after_alleles"]),
                                   int(rows["current_r500", threshold, "DEL"]["after_alleles"]))
            stamps = {p: p.stat().st_mtime_ns for p in same.glob("p*/*/*.vcf.gz")}
            simple.main(args)
            self.assertTrue(all(p.stat().st_mtime_ns == stamp for p, stamp in stamps.items()))
            stamp = (report / "comparison.tsv").stat().st_mtime_ns
            compare.main(cmd)
            self.assertEqual((report / "comparison.tsv").stat().st_mtime_ns, stamp)


if __name__ == "__main__":
    unittest.main()
