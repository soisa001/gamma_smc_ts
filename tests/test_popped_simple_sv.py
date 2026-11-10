"""Synthetic-only tests for the popped simple-indel trial."""
from __future__ import annotations

import csv
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import run_popped_simple_sv as simple


class ClassificationTests(unittest.TestCase):
    def test_shared_flanks_and_one_base_replacements(self):
        self.assertEqual(simple.classify(100, "GGAC", "GGA" + "T" * 60 + "C"),
                         ("INS", 60, 102, "A", "A" + "T" * 60, True))
        self.assertEqual(simple.classify(100, "GGA" + "T" * 60 + "C", "GGAC"),
                         ("DEL", 60, 102, "A" + "T" * 60, "A", True))
        self.assertTrue(simple.classify(100, "A", "T" * 60 + "A")[-1])
        self.assertFalse(simple.classify(100, "A", "T" * 61)[-1])
        self.assertFalse(simple.classify(100, "T" * 61, "A")[-1])
        self.assertFalse(simple.classify(100, "A", "A" + "N" * 60)[-1])
        self.assertFalse(simple.classify(100, "A", "A")[-1])
        self.assertFalse(simple.classify(100, "A", "G")[-1])
        with self.assertRaisesRegex(ValueError, "biallelic"):
            simple.classify(100, "A", "AT,AG")


@unittest.skipUnless(os.environ.get("SV_PILOT_INTEGRATION") == "1", "Requires Linux bcftools/Truvari")
class TrialTests(unittest.TestCase):
    def test_native_partitions_unchanged_replacements_counts_and_resume(self):
        import pysam
        import test_sv_phased_batch as old
        bcftools = os.environ.get("BCFTOOLS_BIN", "bcftools")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "source.vcf"
            old.fixture(source)
            header = "".join(line for line in source.read_text().splitlines(True) if line.startswith("#"))
            low, common = ("1|0:100:10", "0|0:100:20", "0|0:100:30"), ("0|0:100:10", "1|0:100:20", "0|1:100:30")
            definitions = [
                (100, "ins100", "A", "A" + "T" * 100, ".", low),
                (105, "ins80", "A", "A" + "T" * 80, ".", common),
                (2000, "del100", "A" + "C" * 100, "A", ".", low),
                (2005, "del80", "A" + "C" * 80, "A", ".", common),
                (3000, "replacement_ins", "A", "T" * 101, ".", low),
                (3000, "replacement_ins2", "A", "T" * 81, ".", common),
                (4000, "replacement_del", "C" * 101, "A", ".", low),
                (5000, "padded", "GGAC", "GGA" + "T" * 60 + "C", "END=5003;SVTYPE=INS;SVLEN=60", low),
                (6000, "right_anchor", "A", "T" * 60 + "A", ".", low),
                (7000, "small49", "A", "A" + "G" * 49, ".", low),
                (7000, "large50", "A", "A" + "G" * 50, ".", common),
                (8000, "ambiguous", "A", "A" + "N" * 60, ".", low),
                (9000, "symbolic", "N", "<INV>", "SVTYPE=INV;SVLEN=100;END=9100", low),
                (9500, "snp", "A", "G", ".", low),
                (9600, "small_replacement", "A" * 80, "C" * 70, ".", low),
            ]
            lines = {ident: "\t".join(map(str, ("chr1", pos, ident, ref, alt, ".", "PASS", info, "GT:PS:DP", *samples))) + "\n"
                     for pos, ident, ref, alt, info, samples in definitions}
            source.write_text(header + "".join(lines.values()))
            bcf = root / "source.bcf"
            subprocess.run([bcftools, "view", "-Ob", "-o", str(bcf), str(source)], check=True, capture_output=True)
            subprocess.run([bcftools, "index", str(bcf)], check=True, capture_output=True)
            run = root / "run"
            args = ["--root", str(run), "--input-bcf", str(bcf), "--bcftools", bcftools]
            try:
                simple.main(args)
            except BaseException:
                for log in (run / "logs").glob("*.log"):
                    print(log.name, log.read_text()[-4000:])
                raise
            with pysam.VariantFile(run / "prepare/INS.vcf.gz") as handle:
                prepared = {r.id: r for r in handle}
            self.assertIn("small49", prepared)
            self.assertNotIn("replacement_ins", prepared)
            self.assertEqual((prepared["padded"].pos, prepared["padded"].ref, prepared["padded"].stop), (5002, "A", 5002))
            for sample, fields in zip(prepared["padded"].samples.values(), low):
                self.assertTrue(sample.phased)
                self.assertEqual(sample["PS"], 100)
                self.assertEqual(sample["DP"], int(fields.split(":")[2]))
            for tag in ("p70", "p90"):
                with pysam.VariantFile(run / tag / "chr1.sv50.vcf.gz") as handle:
                    retained = {r.id: str(r) for r in handle}
                for name in ("replacement_ins", "replacement_ins2", "replacement_del", "ambiguous", "symbolic"):
                    self.assertEqual(retained[name], lines[name])
                for name in ("small49", "snp", "small_replacement"):
                    self.assertNotIn(name, retained)
                for kind in simple.KINDS:
                    # Every field, including genotypes, agrees with native collapse.
                    kept, gone = root / "direct.vcf.gz", root / "gone.vcf.gz"
                    threshold = "0.7" if tag == "p70" else "0.9"
                    subprocess.run(["truvari", "collapse", "-i", str(run / "prepare" / f"{kind}.vcf.gz"),
                        "-o", str(kept), "-c", str(gone), "-p", threshold, "-P", threshold,
                        "-r", "500", "-O", "0", "-s", "0", "-S", "-1", "--keep", "common", "--gt", "all"],
                        check=True, capture_output=True)
                    for native, name in ((kept, "retained"), (gone, "collapsed")):
                        with gzip.open(native, "rt") as h1, gzip.open(run / tag / kind / f"{name}.vcf.gz", "rt") as h2:
                            self.assertEqual(sorted(l for l in h1 if not l.startswith("#")), sorted(l for l in h2 if not l.startswith("#")))
            with open(run / "comparison.tsv") as handle:
                rows = {(r["similarity"], r["scope"], r["type"]): r for r in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual(int(rows["0.7", "all_INS_DEL", "TOTAL"]["before_sv50"]), 11)
            self.assertEqual(int(rows["0.7", "all_INS_DEL", "TOTAL"]["merged_out"]), 2)
            self.assertEqual(int(rows["0.9", "all_INS_DEL", "TOTAL"]["merged_out"]), 0)
            self.assertEqual(int(rows["0.7", "unchanged_passthrough", "TOTAL"]["merged_out"]), 0)
            saved = {p: p.stat().st_mtime_ns for p in run.rglob("*.vcf.gz")}
            simple.main(args)
            self.assertTrue(all(p.stat().st_mtime_ns == stamp for p, stamp in saved.items()))
            # Alter a valid cached VCF: checksums must trigger a repair, not reuse.
            victim = run / "p70/chr1.sv50.vcf.gz"
            victim.write_bytes(b"corrupt")
            simple.main(args)
            self.assertGreater(victim.stat().st_size, 7)

    def test_empty_del_partition(self):
        import test_sv_phased_batch as old
        bcftools = os.environ.get("BCFTOOLS_BIN", "bcftools")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "source.vcf"
            old.fixture(source)
            header = "".join(l for l in source.read_text().splitlines(True) if l.startswith("#"))
            source.write_text(header + "\t".join(["chr1", "100", "ins", "A", "A" + "T" * 60, ".", "PASS", ".", "GT", "1|0", "0|0", "0|0"]) + "\n")
            compressed = root / "source.vcf.gz"
            subprocess.run([bcftools, "view", "-Oz", "-o", str(compressed), str(source)], check=True, capture_output=True)
            subprocess.run([bcftools, "index", "-t", str(compressed)], check=True, capture_output=True)
            simple.main(["--root", str(root / "run"), "--input-bcf", str(compressed), "--bcftools", bcftools, "--similarities", "0.7"])
            summary = json.loads((root / "run/p70/DEL/summary.json").read_text())
            self.assertEqual(summary["all_sizes"], {"retained": 0, "collapsed": 0})


if __name__ == "__main__":
    unittest.main()
