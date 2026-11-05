"""Integration tests for unmodified Truvari batch outputs and resumability."""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
import run_sv_phased_batch as batch


def concurrent_stage(root, barrier):
    """Model gcloud's shared source-list file through the real stage call site."""
    import time
    root = Path(root)
    os.environ["XDG_CACHE_HOME"] = str(root.parent / "cache")
    args = batch.pilot.parse_args(["--root", str(root), "--scratch-reserve-gb", "0",
        "--source-uri", f"gs://fixture/{root.name}.bcf",
        "--confirm-source-uri", f"gs://fixture/{root.name}.bcf"])
    shared = root.parent / "shared-gcloud-list"
    def command(args, cmd, label, capture=False, **kwargs):
        if label == "rsync":
            shared.write_text(root.name)
            time.sleep(.2)
            if shared.read_text() != root.name:
                raise RuntimeError("Concurrent rsync replaced this chromosome's source listing")
        return "1"
    with mock.patch.object(batch.pilot.storage, "metadata_for", return_value={"bcf": {"size": 1}, "csi": None}), \
         mock.patch.object(batch.pilot.storage, "verify_object"), \
         mock.patch.object(batch.pilot, "manifest"), \
         mock.patch.object(batch.pilot, "command", side_effect=command):
        root.mkdir()
        barrier.wait(timeout=20)
        batch.pilot.stage(args)


@unittest.skipUnless(sys.platform == "linux", "VM staging uses Linux process locks")
class StagingTests(unittest.TestCase):
    def test_concurrent_staging_keeps_source_listings_separate(self):
        import multiprocessing
        with tempfile.TemporaryDirectory() as tmp:
            context = multiprocessing.get_context("fork")
            barrier = context.Barrier(4)
            workers = [context.Process(target=concurrent_stage, args=(str(Path(tmp) / f"chr{i}"), barrier))
                       for i in range(1, 5)]
            try:
                for worker in workers:
                    worker.start()
                for worker in workers:
                    worker.join(timeout=30)
                self.assertEqual([worker.exitcode for worker in workers], [0] * 4)
            finally:
                for worker in workers:
                    if worker.is_alive():
                        worker.terminate()
                        worker.join()


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
    def test_two_chromosomes_match_native_cli_and_resume(self):
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
                split = [r for r in original.values() if r.id == "multi"]
                self.assertEqual({r.samples["s1"]["GT"] for r in split}, {(1, 0), (0, 1)})
                self.assertTrue(all(r.samples["s1"].phased and r.samples["s1"]["PS"] == 9000 for r in split))
                self.assertFalse((root / f"run/{chrom}/p90").exists())
                for percent, expected in ((70, 8),):
                    directory = root / f"run/{chrom}/p{percent}"
                    with pysam.VariantFile(directory / f"{chrom}.vcf.gz") as handle:
                        kept = list(handle)
                    with pysam.VariantFile(directory / f"{chrom}.collapsed.vcf.gz") as handle:
                        removed = list(handle)
                    self.assertEqual(len(kept), expected)
                    self.assertEqual(len(kept) + len(removed), len(original))
                    self.assertEqual({key(r) for r in kept + removed}, set(original))
                    self.assertTrue({"opposite1", "opposite2", "inv"} <= {r.id for r in kept})
                    # Compare every emitted record, including GT and FORMAT, to a
                    # separate invocation of the unmodified pinned Truvari CLI.
                    direct = root / f"{chrom}.{percent}.direct.vcf.gz"
                    direct_removed = root / f"{chrom}.{percent}.removed.vcf.gz"
                    subprocess.run(["truvari", "collapse", "-i", str(root / f"run/{chrom}/prepare/input.vcf.gz"),
                        "-o", str(direct), "-c", str(direct_removed), "-p", str(percent / 100),
                        "-P", str(percent / 100), "-r", "500", "-s", "50", "-S", "-1",
                        "--keep", "common", "--gt", "all"], check=True, capture_output=True)
                    for path, actual in ((direct, kept), (direct_removed, removed)):
                        with pysam.VariantFile(path) as handle:
                            self.assertEqual(sorted(map(str, actual)), sorted(map(str, handle)))
                    for gone in removed:
                        for sample in gone.samples:
                            self.assertEqual(dict(gone.samples[sample]), dict(original[key(gone)].samples[sample]))
                    self.assertFalse((directory / "new.unsorted.vcf").exists())
                self.assertTrue((root / f"{chrom}.bcf").exists())
            saved = {p: p.stat().st_mtime_ns for p in (root / "run").glob("chr*/p*/*.vcf.gz")}
            batch.main(args)
            self.assertTrue(all(p.stat().st_mtime_ns == mtime for p, mtime in saved.items()))



if __name__ == "__main__":
    unittest.main()
