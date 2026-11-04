"""Real Truvari regression and integration tests for GT/phase preservation."""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
import run_sv_phased_batch as batch
import run_chr1_sv_merge as pilot


def fixture(path, chrom="chr1"):
    header = ('##fileformat=VCFv4.2\n'
              f'##contig=<ID={chrom},length=1000000>\n'
              '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n'
              '##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="Length">\n'
              '##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n'
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
              '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\n'
              '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
              '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\n')
    low = ("1|0:100:10", "0|0:100:20", ".|.:.:.")
    common = ("0|0:200:15", "0|1:200:25", "1|0:200:30")
    same1 = ("1|0:4000:40", "0|0:4000:40", "0|0:4000:40")
    same2 = ("0|1:4000:41", "0|0:4000:41", "0|0:4000:41")
    rows = [
        (100, "low_ins", "AG", "A" + "T" * 101, ".", low),
        (105, "common_ins", "AG", "A" + "T" * 81, ".", common),
        (2000, "low_del", "A" + "C" * 101, "AG", ".", low),
        (2005, "common_del", "A" + "C" * 81, "AG", ".", common),
        (4000, "opposite1", "A", "A" + "G" * 100, ".", same1),
        (4005, "opposite2", "A", "A" + "G" * 100, ".", same2),
        (5000, "ambiguous", "A", "A" + "N" * 60, ".", same1),
        (7000, "small", "A", "A" + "T" * 49, ".", same1),
        (9000, "multi", "A", "A" + "G" * 60 + ",A" + "C" * 90, ".",
         ("1|2:9000:50", "0|0:9000:50", "0|0:9000:50")),
        (12000, "inv", "N", "<INV>", "SVTYPE=INV;SVLEN=100;END=12100", same1),
    ]
    with open(path, "w") as handle:
        handle.write(header)
        for pos, ident, ref, alt, info, samples in rows:
            handle.write("\t".join([chrom, str(pos), ident, ref, alt, ".", "PASS", info, "GT:PS:DP", *samples]) + "\n")


def key(record):
    return record.contig, record.pos, record.ref, record.alts


@unittest.skipUnless(os.environ.get("SV_PILOT_INTEGRATION") == "1", "Requires Linux bcftools/Truvari")
class PhasedTests(unittest.TestCase):
    def test_native_bug_and_phase_preserving_wrapper_through_actual_cli(self):
        import pysam
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raw = root / "source.vcf"
            fixture(raw)
            bcftools = os.environ.get("BCFTOOLS_BIN", "bcftools")
            source = root / "source.vcf.gz"
            subprocess.run([bcftools, "view", "-i", "POS<1000", "-Oz", "-o", str(source), str(raw)], check=True, capture_output=True)
            subprocess.run([bcftools, "index", "-t", str(source)], check=True, capture_output=True)
            for label, command in (("native", ["truvari", "collapse"]),
                                   ("phased", [sys.executable, str(SCRIPTS / "truvari_collapse_phased.py")])):
                subprocess.run(command + ["-i", str(source), "-o", str(root / f"{label}.vcf.gz"),
                    "-c", str(root / f"{label}.collapsed.vcf.gz"), "-p", ".7", "-P", ".7", "-r", "500",
                    "-s", "50", "-S", "-1", "--keep", "common", "--gt", "all"], check=True, capture_output=True)
                subprocess.run([bcftools, "index", "-t", str(root / f"{label}.vcf.gz")], check=True, capture_output=True)
            with pysam.VariantFile(root / "native.vcf.gz") as handle:
                native = list(handle)
            with pysam.VariantFile(root / "phased.vcf.gz") as handle:
                phased = list(handle)
            self.assertEqual(len(native), 1)
            self.assertEqual(key(native[0]), key(phased[0]))
            # Preserve this release's native --keep common choice, even though
            # 5.5.0 sorts its AC comparator ascending and chooses low_ins here.
            self.assertEqual(phased[0].id, native[0].id)
            # Truvari 5.5.0 does not write the summed GT back for this donor sample.
            self.assertEqual(native[0].samples["s2"]["GT"], (0, 0))
            for sample, gt, ps, dp in (("s1", (1, 0), 100, 10), ("s2", (0, 1), 200, 25), ("s3", (1, 0), 200, 30)):
                call = phased[0].samples[sample]
                self.assertEqual(call["GT"], gt)
                self.assertTrue(call.phased)
                self.assertEqual((call["PS"], call["DP"]), (ps, dp))
            self.assertEqual(phased[0].info["AC"], (3,))
            self.assertEqual(phased[0].info["AN"], 6)
            self.assertEqual(phased[0].info["AF"], (.5,))

    def test_two_chromosomes_thresholds_phase_membership_and_resume(self):
        import pysam
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bcftools = os.environ.get("BCFTOOLS_BIN", "bcftools")
            for chrom in ("chr1", "chr2"):
                raw = root / f"{chrom}.vcf"
                fixture(raw, chrom)
                subprocess.run([bcftools, "view", "-Ob", "-o", str(root / f"{chrom}.bcf"), str(raw)], check=True, capture_output=True)
                subprocess.run([bcftools, "index", str(root / f"{chrom}.bcf")], check=True, capture_output=True)
            args = ["--root", str(root / "run"), "--chromosomes", "chr1", "chr2", "--jobs", "2", "--phase", "all",
                    "--input-template", str(root / "{chrom}.bcf"), "--bcftools", bcftools]
            try:
                batch.main(args)
            except BaseException:
                for log in (root / "run").rglob("*.log"):
                    print(log, log.read_text()[-6000:])
                raise
            for chrom in ("chr1", "chr2"):
                with pysam.VariantFile(root / f"run/{chrom}/prepare/input.vcf.gz") as handle:
                    original = {key(record): record for record in handle}
                self.assertEqual(len(original), 10)
                for percent, expected in ((70, 8), (90, 10)):
                    directory = root / f"run/{chrom}/p{percent}"
                    with pysam.VariantFile(directory / f"{chrom}.vcf.gz") as handle:
                        kept = list(handle)
                    with pysam.VariantFile(directory / f"{chrom}.collapsed.vcf.gz") as handle:
                        removed = list(handle)
                    self.assertEqual(len(kept), expected)
                    self.assertEqual(len(kept) + len(removed), len(original))
                    self.assertEqual({key(r) for r in kept + removed}, set(original))
                    self.assertTrue({"opposite1", "opposite2", "inv"} <= {r.id for r in kept})
                    for rep in kept:
                        members = [original[key(rep)]] + [original[key(r)] for r in removed
                            if rep.info.get("CollapseId") in r.info.get("MatchId", ())]
                        for sample in rep.samples:
                            donors = [r.samples[sample] for r in members if 1 in r.samples[sample]["GT"]]
                            self.assertLessEqual(len(donors), 1)
                            donor = donors[0] if donors else original[key(rep)].samples[sample]
                            self.assertEqual(rep.samples[sample]["GT"], donor["GT"],
                                (rep.id, sample, dict(rep.info), [(r.id, dict(r.info)) for r in removed]))
                            self.assertEqual(rep.samples[sample].phased, donor.phased)
                            self.assertEqual(rep.samples[sample]["PS"], donor["PS"])
                            self.assertEqual(rep.samples[sample]["DP"], donor["DP"])
                    for gone in removed:
                        for sample in gone.samples:
                            self.assertEqual(dict(gone.samples[sample]), dict(original[key(gone)].samples[sample]))
                    self.assertFalse((directory / "new.unsorted.vcf").exists())
                self.assertTrue((root / f"{chrom}.bcf").exists())
            saved = {p: p.stat().st_mtime_ns for p in (root / "run").glob("chr*/p*/*.vcf.gz")}
            batch.main(args)
            self.assertTrue(all(p.stat().st_mtime_ns == mtime for p, mtime in saved.items()))

    def test_consolidator_preserves_partial_unphased_and_haploid_calls_and_rejects_overlap(self):
        import pysam
        from types import SimpleNamespace
        from truvari_collapse_phased import consolidate
        header = pysam.VariantHeader()
        header.contigs.add("chr1")
        header.formats.add("GT", 1, "String", "Genotype")
        header.formats.add("PS", 1, "Integer", "Phase set")
        header.add_sample("s1")
        for gt, phased, ps in (((1, None), True, 100), ((0, 1), False, None), ((1,), False, None)):
            a = header.new_record(contig="chr1", start=100, alleles=("A", "ATTT"))
            b = header.new_record(contig="chr1", start=105, alleles=("A", "ATTT"))
            a.samples[0]["GT"] = (0, 0)
            a.samples[0]["PS"] = 500
            b.samples[0]["GT"] = gt
            b.samples[0]["PS"] = ps
            b.samples[0].phased = phased
            expected_phase = b.samples[0].phased
            result, count = consolidate(a, [SimpleNamespace(comp=b)])
            self.assertEqual(count, 1)
            self.assertEqual(result.samples[0]["GT"], gt)
            self.assertEqual(result.samples[0].phased, expected_phase)
            self.assertEqual(result.samples[0]["PS"], ps)
            with self.assertRaisesRegex(ValueError, "Multiple ALT"):
                consolidate(result, [SimpleNamespace(comp=b)])


if __name__ == "__main__":
    unittest.main()
