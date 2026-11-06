"""Representation metrics and real bcftools decoding on synthetic alleles."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
import inspect_sv_representation as inspect


def records():
    return [
        (100, "ACCCC", "A,AT,ACCCC"),
        (200, "A", "ATTT,<INS>,*"),
        (300, "A" * 100, "A" + "C" * 50),
        (400, "A", "A" + "G" * 60),
        (500, "A" * 100, "A"),
        (600, "A" * 100, "A" + "C" * 20),
        (600, "A" * 100, "A" + "C" * 20),
        (600, "A" * 100, "A" + "G" * 20),
    ]


class MetricsTests(unittest.TestCase):
    def test_padding_is_not_complexity(self):
        pure, _ = inspect.describe("ACCCGT", "AGT")
        self.assertFalse(pure["complex"])
        self.assertEqual((pure["shared_prefix_bp"], pure["shared_suffix_bp"]), (1, 2))
        self.assertEqual((pure["residual_ref_bp"], pure["residual_alt_bp"]), (3, 0))
        complex_, _ = inspect.describe("A" * 1000, "A" * 999 + "C" * 101)
        self.assertTrue(complex_["complex"])
        self.assertEqual(complex_["net_bp"], 100)
        self.assertEqual(complex_["shared_flank_bp"], 999)
        self.assertEqual(inspect.describe("AC", "AC")[0]["type"], "REF_EQUAL")
        self.assertTrue(inspect.describe("A", "AN")[0]["ambiguous"])

    def test_multiallelics_duplicates_and_posthoc_scope(self):
        result = inspect.summarize((f"chr1\t{p}\t{r}\t{a}\n" for p, r, a in records()), "chr1")
        self.assertEqual(result["site_counts"]["records"], 8)
        self.assertEqual(result["site_counts"]["multiallelic_records"], 2)
        self.assertEqual(result["site_counts"]["complex_net_lt50_but_residual_span_ge50"], 1)
        rows = {(r["scope"], r["type"]): r for r in result["allele_summary"]}
        dels = rows["net_ge50", "DEL"]
        self.assertEqual((dels["alleles"], dels["unique_alleles"], dels["exact_duplicate_alleles"]), (4, 3, 1))
        self.assertEqual(dels["complex_percent"], 75)
        self.assertEqual(dels["pos_ref_groups_with_multiple_unique_alts"], 1)
        self.assertEqual(rows["net_ge50", "INS"]["alleles"], 1)
        self.assertEqual(inspect.quantiles({1: 2, 5: 2})["median"], 3)
        with self.assertRaisesRegex(ValueError, "sorted"):
            inspect.summarize(["chr1\t2\tA\tAT\n", "chr1\t1\tA\tAT\n"], "chr1")

    @unittest.skipUnless(os.environ.get("SV_PILOT_INTEGRATION") == "1", "Requires Linux bcftools")
    def test_bcf_scan_and_cache_without_cloud(self):
        bcftools = os.environ["BCFTOOLS_BIN"]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            vcf = root / "fixture.vcf"
            vcf.write_text("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=10000>\n"
                           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
                           "".join(f"chr1\t{p}\t.\t{r}\t{a}\t.\tPASS\t.\n" for p, r, a in records()))
            bcf = root / "fixture.bcf"
            subprocess.run([bcftools, "view", "-Ob", "-o", str(bcf), str(vcf)], check=True, capture_output=True)
            subprocess.run([bcftools, "index", str(bcf)], check=True, capture_output=True)
            original_stamp = bcf.stat().st_mtime_ns
            cmd = [sys.executable, str(SCRIPTS / "inspect_sv_representation.py"), "--root", str(root / "run"),
                   "--input-bcf", str(bcf), "--compare-input", str(bcf), "--bcftools", bcftools]
            completed = subprocess.run(cmd, text=True, capture_output=True)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            report = root / "run/popped/summary.json"
            data = json.loads(report.read_text())
            self.assertEqual(data["site_counts"]["records"], 8)
            stamp = report.stat().st_mtime_ns
            subprocess.run(cmd, check=True, capture_output=True)
            self.assertEqual(report.stat().st_mtime_ns, stamp)
            self.assertEqual(bcf.stat().st_mtime_ns, original_stamp)


if __name__ == "__main__":
    unittest.main()
