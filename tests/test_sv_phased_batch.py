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
    def test_sub50_alleles_participate_before_final_selection(self):
        import pysam
        import json
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bcftools = os.environ.get("BCFTOOLS_BIN", "bcftools")
            raw = root / "chr1.vcf"
            fixture(raw)
            header = "\n".join(line for line in raw.read_text().splitlines() if line.startswith("#")) + "\n"
            rows = [
                (100, "small49", "A", "A" + "T" * 49, "1|0", "0|0", "0|0"),
                (100, "large50", "A", "A" + "T" * 50, "0|0", "1|0", "1|0"),
                (2000, "snp", "A", "G", "1|0", "0|0", "0|0"),
                (4000, "complex10", "A" * 80, "C" * 70, "1|0", "0|0", "0|0"),
            ]
            raw.write_text(header + "".join("\t".join(["chr1", str(pos), ident, ref, alt,
                ".", "PASS", ".", "GT", s1, s2, s3]) + "\n"
                for pos, ident, ref, alt, s1, s2, s3 in rows))
            subprocess.run([bcftools, "view", "-Ob", "-o", str(root / "chr1.bcf"), str(raw)], check=True, capture_output=True)
            subprocess.run([bcftools, "index", str(root / "chr1.bcf")], check=True, capture_output=True)
            batch.main(["--root", str(root / "run"), "--chromosomes", "chr1", "--jobs", "1", "--phase", "all",
                        "--similarities", "0.7", "0.9", "--input-template", str(root / "{chrom}.bcf"), "--bcftools", bcftools])
            with pysam.VariantFile(root / "run/chr1/prepare/input.vcf.gz") as handle:
                self.assertEqual({r.id for r in handle}, {r[1] for r in rows})
            for tag in ("p70", "p90"):
                directory = root / "run/chr1" / tag
                with pysam.VariantFile(directory / "chr1.collapsed.vcf.gz") as handle:
                    self.assertEqual([r.id for r in handle], ["small49"])
                with pysam.VariantFile(directory / "chr1.vcf.gz") as handle:
                    self.assertEqual({r.id for r in handle}, {"large50", "snp", "complex10"})
                with pysam.VariantFile(directory / "chr1.sv50.vcf.gz") as handle:
                    self.assertEqual([r.id for r in handle], ["large50"])
                counts = json.loads((directory / "summary.json").read_text())["counts"]
                self.assertEqual(counts["final_sv50"], 1)
                self.assertEqual(counts["retained_excluded_by_final_filter"], 2)

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
                self.assertEqual(len(original), 11)
                self.assertIn("small", {r.id for r in original.values()})
                split = [r for r in original.values() if r.id == "multi"]
                self.assertEqual({r.samples["s1"]["GT"] for r in split}, {(1, 0), (0, 1)})
                self.assertTrue(all(r.samples["s1"].phased and r.samples["s1"]["PS"] == 9000 for r in split))
                self.assertFalse((root / f"run/{chrom}/p90").exists())
                for percent, expected in ((70, 9),):
                    directory = root / f"run/{chrom}/p{percent}"
                    with pysam.VariantFile(directory / f"{chrom}.vcf.gz") as handle:
                        kept = list(handle)
                    with pysam.VariantFile(directory / f"{chrom}.collapsed.vcf.gz") as handle:
                        removed = list(handle)
                    self.assertEqual(len(kept), expected)
                    with pysam.VariantFile(directory / f"{chrom}.sv50.vcf.gz") as handle:
                        final = list(handle)
                    self.assertEqual(len(final), 7)
                    self.assertTrue(all(abs(len(r.ref) - len(r.alts[0])) >= 50 for r in final))
                    self.assertNotIn("small", {r.id for r in final})
                    self.assertEqual(len(kept) + len(removed), len(original))
                    self.assertEqual({key(r) for r in kept + removed}, set(original))
                    self.assertTrue({"opposite1", "opposite2", "inv"} <= {r.id for r in kept})
                    # Compare every emitted record, including GT and FORMAT, to a
                    # separate invocation of the unmodified pinned Truvari CLI.
                    direct = root / f"{chrom}.{percent}.direct.vcf.gz"
                    direct_removed = root / f"{chrom}.{percent}.removed.vcf.gz"
                    subprocess.run(["truvari", "collapse", "-i", str(root / f"run/{chrom}/prepare/input.vcf.gz"),
                        "-o", str(direct), "-c", str(direct_removed), "-p", str(percent / 100),
                        "-P", str(percent / 100), "-r", "500", "-s", "0", "-S", "-1",
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

            import compare_sv_native
            import csv
            import json
            output = compare_sv_native.compare(root / "run", "chr1", bcftools, threads=4)
            with open(output / "changes.tsv") as handle:
                stats = {(r["type"], r["metric"]): r for r in csv.DictReader(handle, delimiter="\t")}
            for kind, before, after, bp_before, bp_after in (("INS", 7, 6, 590, 510), ("DEL", 2, 1, 180, 100)):
                self.assertEqual(float(stats[kind, "sv_records"]["before"]), before)
                self.assertEqual(float(stats[kind, "sv_records"]["after"]), after)
                self.assertAlmostEqual(float(stats[kind, "sv_records"]["percent_reduction"]), 100 * (before - after) / before)
                self.assertEqual(float(stats[kind, "raw_bp_sum_records"]["before"]), bp_before)
                self.assertEqual(float(stats[kind, "raw_bp_sum_records"]["after"]), bp_after)
            report = json.loads((output / "summary.json").read_text())
            self.assertEqual(report["all_record_counts"], {"before": 11, "after": 9, "removed": 2})
            self.assertEqual(report["posthoc_min_raw_bp"], 50)
            stamp = (output / "changes.tsv").stat().st_mtime_ns
            compare_sv_native.compare(root / "run", "chr1", bcftools, threads=4)
            self.assertEqual((output / "changes.tsv").stat().st_mtime_ns, stamp)
            import audit_sv_native
            from types import SimpleNamespace
            audit_root = root / "run/chr1/audit_native_p70"
            audit = audit_sv_native.run(SimpleNamespace(batch_root=root / "run", root=audit_root,
                                                       chrom="chr1", bcftools=bcftools))
            self.assertTrue(audit["exact_site_multiset_conserved"])
            self.assertTrue(audit["numcollapsed_matches_membership"])
            self.assertEqual(audit["record_counts"], {"before": 11, "after": 9, "removed": 2})
            for kind in ("INS", "DEL"):
                pairs = audit["removed_vs_representative"][kind]
                self.assertEqual(pairs["counts"]["pairs"], 1)
                self.assertAlmostEqual(pairs["distributions"]["raw_length_ratio"]["min"], .8)
                self.assertEqual(audit["merged_cluster_sizes"][kind]["max"], 2)
            # Adding chr1's 0.9 run must reuse 0.7 and leave other chromosomes alone.
            aggregate = root / "run/counts.tsv"
            aggregate.write_text("existing full-batch summary\n")
            batch.main(args + ["--chromosomes", "chr1", "--similarities", "0.7", "0.9"])
            self.assertTrue(all(p.stat().st_mtime_ns == mtime for p, mtime in saved.items()))
            self.assertFalse((root / "run/chr2/p90").exists())
            self.assertEqual(aggregate.read_text(), "existing full-batch summary\n")
            output90 = compare_sv_native.compare(root / "run", "chr1", bcftools, similarity=.9)
            joined = compare_sv_native.side_by_side(root / "run", "chr1", {.7: output, .9: output90})
            with open(joined) as handle:
                combined = {(r["type"], r["metric"]): r for r in csv.DictReader(handle, delimiter="\t")}
            for kind, p70, p90 in (("INS", 6, 7), ("DEL", 1, 2)):
                row = combined[kind, "sv_records"]
                self.assertEqual((int(row["after_p70"]), int(row["after_p90"])), (p70, p90))
                self.assertEqual(float(row["p90_minus_p70"]), 1)
                self.assertEqual(float(row["percent_reduction_p90"]), 0)
            audit90 = audit_sv_native.run(SimpleNamespace(batch_root=root / "run",
                root=root / "run/chr1/audit_native_p90", chrom="chr1", bcftools=bcftools, similarity=.9))
            self.assertEqual(audit90["record_counts"], {"before": 11, "after": 11, "removed": 0})
            # A structurally readable but altered retained VCF must fail membership validation.
            kept_path = root / "run/chr1/p70/chr1.vcf.gz"
            with pysam.VariantFile(kept_path) as handle:
                records, header = list(handle), handle.header.copy()
            records[0].pos += 1
            with pysam.VariantFile(kept_path, "wz", header=header) as handle:
                for record in records:
                    handle.write(record)
            with self.assertRaisesRegex(ValueError, "site multiset differs"):
                compare_sv_native.compare(root / "run", "chr1", bcftools, threads=4)



if __name__ == "__main__":
    unittest.main()
