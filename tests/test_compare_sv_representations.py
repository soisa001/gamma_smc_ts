import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
import compare_sv_representations as comparison


class RepresentationTests(unittest.TestCase):
    def test_trim_and_alignment(self):
        pos, ref, alt, _ = comparison.trim(199, "GA", "GATTT")
        self.assertEqual((pos, ref, alt), (201, "", "TTT"))
        self.assertEqual(comparison.alignment("CCCC", "", 100)["operation_bp"], {"D": 4})
        self.assertEqual(comparison.alignment("", "TTT", 100)["operation_bp"], {"I": 3})
        self.assertIn("skipped", comparison.alignment("A" * 101, "T", 100))

    def test_native_del_score_does_not_compare_alt(self):
        a = {"pos": 100, "ref": "A" * 80, "alt": "A" + "C" * 19, "svtype": ".", "svlen": "."}
        b = {**a, "alt": "A" + "G" * 19}
        for row in (a, b):
            pos, ref, alt, _ = comparison.trim(row["pos"], row["ref"], row["alt"])
            row["trimmed"] = {"pos": pos, "ref": ref, "alt": alt}
        scores = comparison.pair_scores(a, b, 1000)
        self.assertEqual(scores["native_seqsim"], 1)
        self.assertEqual(scores["trimmed_alt_direct_score"], .5)
        self.assertTrue(scores["passes_p90_sequence_size_only"])

    @unittest.skipUnless(os.environ.get("SV_PILOT_INTEGRATION") == "1", "Requires Linux bcftools")
    def test_three_files_row_packing_and_trimmed_equivalence(self):
        bcftools = os.environ["BCFTOOLS_BIN"]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            inputs = {
                "bubble": [(100, "ACCCC", "A,AT"), (200, "A", "ATTT")],
                "popped": [(100, "ACCCC", "A"), (100, "ACCCC", "AT"), (200, "A", "ATTT"), (200, "A", "ATTT")],
                "imputed": [(100, "ACCCC", "A"), (199, "GA", "GATTT"), (300, "ATTT", "A")],
            }
            cmd = [sys.executable, str(SCRIPTS / "compare_sv_representations.py"), "--root", str(root / "run"),
                   "--min-net-bp", "1", "--bcftools", bcftools, "--examples", "2"]
            for label, rows in inputs.items():
                raw, bcf = root / f"{label}.vcf", root / f"{label}.bcf"
                raw.write_text("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=10000>\n"
                               "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
                               "".join(f"chr1\t{p}\t{label}{i}\t{r}\t{a}\t.\tPASS\t.\n" for i, (p,r,a) in enumerate(rows)))
                subprocess.run([bcftools, "view", "-Ob", "-o", str(bcf), str(raw)], check=True, capture_output=True)
                subprocess.run([bcftools, "index", str(bcf)], check=True, capture_output=True)
                cmd += ["--" + label, str(bcf)]
            result = subprocess.run(cmd, text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            import csv
            with (root / "run/pairwise_summary.tsv").open() as handle:
                table = {(r["left"], r["right"], r["type"], r["comparison"]): r for r in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual(table["bubble", "popped", "DEL", "exact"]["sets_equal"], "True")
            self.assertEqual(table["bubble", "popped", "INS", "exact"]["sets_equal"], "True")
            self.assertEqual(table["bubble", "popped", "INS", "exact"]["multisets_equal"], "False")
            self.assertEqual(table["bubble", "imputed", "INS", "exact"]["shared_unique"], "0")
            self.assertEqual(table["bubble", "imputed", "INS", "canonical"]["shared_unique"], "1")
            self.assertEqual(table["bubble", "imputed", "DEL", "exact"]["left_only_unique"], "1")
            examples = [json.loads(line) for line in (root / "run/allele_examples.jsonl").read_text().splitlines()]
            self.assertTrue(any(r["original_alt_count"] == 2 for r in examples if r["dataset"] == "bubble"))
            self.assertTrue(all(r["original_alt_count"] == 1 for r in examples if r["dataset"] == "popped"))
            stamps = {p: p.stat().st_mtime_ns for p in (root / "run").glob("*/alleles.sqlite")}
            result = subprocess.run(cmd, text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue(all(p.stat().st_mtime_ns == stamp for p, stamp in stamps.items()))


if __name__ == "__main__":
    unittest.main()
