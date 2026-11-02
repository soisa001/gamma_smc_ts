"""Portable logic tests plus opt-in real bcftools/Truvari integration test."""
from __future__ import annotations

import gzip
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
spec = importlib.util.spec_from_file_location("sv_pilot", SCRIPTS / "run_chr1_sv_merge.py")
sv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sv)


def fields(pos=100, ref="A", alt=None, ident="a", gt=("0/1", "0/0"), info="."):
    return ["chr1", str(pos), ident, ref, alt or "A" + "C" * 100, ".", "PASS", info, "GT", *gt]


def row(ident, pos, end=None, kind="INS", carrier_bits=1, ac=1, collapse=None, num=0):
    return {"id": ident, "pos": pos, "type": kind, "size": 100, "start0": pos,
            "end0": pos if end is None else end, "status": "eligible", "ac": ac,
            "an": 4, "af": ac / 4, "samples": 2, "missing_samples": 0,
            "carriers": carrier_bits.bit_count(), "carrier_bits": carrier_bits,
            "gt_sha256": str(carrier_bits), "collapse_id": collapse, "num_collapsed": num}


def write_vcf(path, entries):
    header = ('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=250000000>\n'
              '##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="Length">\n'
              '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n'
              '##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n'
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
              '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n')
    with sv.open_text(path, "wt") as handle:
        handle.write(header)
        handle.writelines("\t".join(e) + "\n" for e in entries)


class AlleleTests(unittest.TestCase):
    def test_trimmed_boundaries_and_exact_size_limits(self):
        # Insertion after anchor; padding at both ends must not become interval span.
        self.assertEqual(sv.allele("AC", "A" + "T" * 50 + "C", 100, {}),
                         ("INS", 50, 100, 100, "eligible"))
        self.assertEqual(sv.allele("A" + "T" * 50 + "C", "AC", 100, {}),
                         ("DEL", 50, 100, 150, "eligible"))
        self.assertIsNone(sv.allele("A", "A" + "T" * 49, 100, {}))
        self.assertEqual(sv.allele("A", "A" + "T" * 100, 100, {}, maximum=100)[-1], "eligible")
        self.assertEqual(sv.allele("A", "A" + "T" * 101, 100, {}, maximum=100)[-1], "above_sizemax")

    def test_no_silent_symbolic_or_ambiguous_sequence_merging(self):
        self.assertEqual(sv.allele("N", "<INS>", 100, {})[:2], ("INS", None))
        self.assertEqual(sv.allele("N", "<DEL>", 100, {"END": "250"}),
                         ("DEL", 150, 99, 250, "symbolic"))
        self.assertEqual(sv.allele("A", "A" + "N" * 80, 100, {})[-1], "ambiguous_sequence")
        self.assertEqual(sv.allele("AG", "A" + "T" * 80, 100, {})[-1], "complex_replacement")
        self.assertEqual(sv.allele("A", "A]chr2:20]", 100, {})[0], "BND")
        self.assertEqual(sv.allele("A", ".A", 100, {})[0], "BND")
        for alt in ("*", "<NON_REF>", "."):
            self.assertIsNone(sv.allele("A", alt, 100, {}))

    def test_genotype_counts_partial_missing_haploid_and_format_order(self):
        f = fields(gt=("9:1|.", "3:1", "8:0/0", "2:./."))
        f[8] = "DP:GT"
        result = sv.genotypes(f)
        self.assertEqual((result["ac"], result["an"], result["missing_samples"], result["carriers"]), (2, 4, 2, 2))
        self.assertEqual(result["carrier_bits"], 3)
        self.assertEqual(result["af"], .5)
        self.assertIsNone(sv.genotypes(fields(gt=("./.", ".")))["af"])
        with self.assertRaises(ValueError):
            sv.genotypes(fields(gt=("1/2",)))


class ComparisonTests(unittest.TestCase):
    def test_distances_include_ties_and_nested_interval_union(self):
        data = [row("a", 100, 1000), row("b", 100, 200), row("c", 300, 400), row("d", 1100, 1200)]
        gaps, nearest, intervals, overlaps = sv.distances(data)
        self.assertEqual(gaps, [0, 200, 800])
        self.assertEqual(nearest, [0, 0, 200, 800])
        self.assertEqual(intervals, [0, 0, 100])
        self.assertEqual(overlaps, 2)
        self.assertEqual(sv.distances([row("a", 1)]), ([], [], [], 0))
        self.assertEqual(sv.distances([]), ([], [], [], 0))

    def test_cluster_ids_are_namespaced_by_type_and_conserve_all_records(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            original = {r["id"]: r for r in [row("a", 100), row("b", 105),
                        row("c", 1000, kind="DEL", carrier_bits=2), row("d", 1005, kind="DEL")]}
            new = {"a": row("a", 100, collapse="1.0", num=1),
                   "c": row("c", 1000, kind="DEL", carrier_bits=3, ac=2, collapse="1.0", num=1)}
            paths = {}
            for kind, ident in (("INS", "b"), ("DEL", "d")):
                paths[kind] = root / f"{kind}.vcf"
                write_vcf(paths[kind], [fields(ident=ident, info="MatchId=1.0;PctSeqSimilarity=0.9;PctSizeSimilarity=0.8")])
            links, clusters = sv.audit_clusters(original, new, paths)
            self.assertEqual(len(links), 4)
            self.assertEqual([c["samples_carrying_multiple_members"] for c in clusters], [1, 0])
            self.assertEqual([c["lost_carriers"] for c in clusters], [0, 0])
            write_vcf(paths["INS"], [])
            with self.assertRaisesRegex(ValueError, "Conservation"):
                sv.audit_clusters(original, new, paths)

    def test_quantiles_empty_and_interpolation(self):
        self.assertIsNone(sv.quantiles([])["median"])
        self.assertEqual(sv.quantiles([0, 10, None])["median"], 5)

    def test_checksum_cache_invalidates_corruption_with_same_size_mtime(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            path = directory / "data.txt"
            path.write_text("good")
            old = path.stat()
            sv.manifest(directory, {"samples": "abc"}, [path])
            self.assertIsNotNone(sv.valid_manifest(directory, {"samples": "abc"}))
            self.assertIsNone(sv.valid_manifest(directory, {"samples": "changed"}))
            path.write_text("evil")
            os.utime(path, ns=(old.st_atime_ns, old.st_mtime_ns))
            self.assertIsNone(sv.valid_manifest(directory, {"samples": "abc"}))

    def test_stage_requires_exact_confirmation_before_cloud_access(self):
        with tempfile.TemporaryDirectory() as tmp:
            args = sv.parse_args(["--root", tmp, "--phase", "stage"])
            with mock.patch.object(sv.storage, "metadata_for") as remote:
                with self.assertRaisesRegex(ValueError, "confirm-source-uri"):
                    sv.stage(args)
                remote.assert_not_called()

    def test_rsync_only_the_exact_chr1_object_and_index(self):
        cmd = sv.storage.sync_command("gcloud", sv.SOURCE, Path("inputs"), "billing")
        pattern = cmd[cmd.index("--exclude") + 1]
        basename = sv.SOURCE.rsplit("/", 1)[1]
        import re
        self.assertIsNone(re.search(pattern, basename))
        self.assertIsNone(re.search(pattern, basename + ".csi"))
        self.assertIsNotNone(re.search(pattern, basename.replace("chr1.", "chr10.")))
        self.assertEqual(cmd[-2:], ["--billing-project", "billing"])

    def test_stage_does_not_require_or_infer_billing_project(self):
        with tempfile.TemporaryDirectory() as tmp:
            for environment in ({}, {"GOOGLE_PROJECT": "do-not-infer"}):
                with self.subTest(environment=environment), mock.patch.dict(os.environ, environment, clear=True):
                    args = sv.parse_args(["--root", tmp, "--phase", "stage", "--confirm-source-uri", sv.SOURCE])
                    self.assertIsNone(args.billing_project)
                    cmd = sv.storage.sync_command(args.gcloud, args.source_uri, Path(tmp), args.billing_project)
                    self.assertNotIn("--billing-project", cmd)
                    # Reaching metadata retrieval proves no local billing guard blocks staging.
                    with mock.patch.object(sv.storage, "metadata_for", side_effect=RuntimeError("metadata reached")) as remote:
                        with self.assertRaisesRegex(RuntimeError, "metadata reached"):
                            sv.stage(args)
                        remote.assert_called_once_with(args, sv.SOURCE)
            override = sv.parse_args(["--root", tmp, "--billing-project", "explicit-project"])
            self.assertEqual(override.billing_project, "explicit-project")

    def test_metadata_queries_omit_billing_flag_by_default(self):
        metadata = json.dumps({"size": "1", "generation": "123", "crc32c": "AAAAAA=="})
        with mock.patch.object(sv.storage, "command", return_value=metadata) as command:
            sv.storage.describe("gcloud", sv.SOURCE)
            self.assertNotIn("--billing-project", command.call_args.args[0])

    def test_real_comparison_outputs_and_cache(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for sub in ("prepare", "merge", "logs"):
                (root / sub).mkdir()
            baseline = [fields(ident="a", info="PM_STATUS=eligible"),
                        fields(pos=105, ident="b", info="PM_STATUS=eligible")]
            before = root / "prepare" / "original.vcf.gz"
            after = root / "merge" / "new.vcf.gz"
            write_vcf(before, baseline)
            write_vcf(after, [fields(ident="a", info="PM_STATUS=eligible;CollapseId=1.0;NumCollapsed=1")])
            removed = root / "merge" / "INS.removed.vcf.gz"
            write_vcf(removed, [fields(pos=105, ident="b", info="MatchId=1.0;PctSeqSimilarity=1;PctSizeSimilarity=1")])
            empty = root / "merge" / "DEL.removed.vcf.gz"
            write_vcf(empty, [])
            prepared = sv.manifest(root / "prepare", {"samples": "s1,s2"}, [before])
            sv.manifest(root / "merge", {"prepare": prepared}, [after, removed, empty])
            args = sv.parse_args(["--root", str(root), "--phase", "compare"])
            sv.compare(args)
            summary = json.loads((root / "compare" / "summary.json").read_text())
            self.assertEqual(summary["counts"][0]["original"], 2)
            self.assertEqual(summary["counts"][0]["new"], 1)
            self.assertEqual(summary["clusters_with_same_sample_cocarriage"], 1)
            output = root / "compare" / "counts.tsv"
            mtime = output.stat().st_mtime_ns
            sv.compare(args)
            self.assertEqual(output.stat().st_mtime_ns, mtime)
            self.assertIsNotNone(sv.valid_manifest(root / "compare"))
            # Optional synthetic-only visual QA copy outside the repository.
            if os.environ.get("SV_PILOT_QA_DIR"):
                shutil.copytree(root / "compare", os.environ["SV_PILOT_QA_DIR"], dirs_exist_ok=True)


@unittest.skipUnless(os.environ.get("SV_PILOT_INTEGRATION") == "1", "Set SV_PILOT_INTEGRATION=1 on Linux with bcftools/Truvari")
class IntegrationTests(unittest.TestCase):
    def test_real_tools_multiallelic_counts_and_idempotence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raw = root / "synthetic.vcf"
            write_vcf(raw, [
                fields(pos=10, alt="C", ident="snv"),
                fields(pos=20, alt="A" + "C" * 49, ident="small"),
                fields(pos=100, alt="A" * 101, ident="ins1", gt=("1|0", "0|0")),
                fields(pos=105, alt="A" * 111, ident="ins2", gt=("1|0", "0|1")),
                fields(pos=1000, ref="A" * 101, alt="A", ident="del1"),
                fields(pos=1005, ref="A" * 111, alt="A", ident="del2", gt=("0/0", "0/1")),
                fields(pos=10000, alt="A" + "G" * 100, ident="distant"),
                fields(pos=20000, alt="A" + "C" * 60 + ",A" + "T" * 90, ident="multi", gt=("1|2", "0|0")),
                fields(pos=30000, alt="A" + "N" * 80, ident="ambiguous"),
                fields(pos=40000, ref="N", alt="<INS>", ident="symbolic", info="SVLEN=100;SVTYPE=INS"),
                fields(pos=50000, ref="N", alt="<INV>", ident="inv", info="END=50200;SVTYPE=INV"),
                fields(pos=60000, ref="A" * 100002, alt="A", ident="large"),
            ])
            bcf = root / "synthetic.bcf"
            bcftools = os.environ.get("BCFTOOLS_BIN", "bcftools")
            subprocess.run([bcftools, "view", "-Ob", "-o", str(bcf), str(raw)], check=True, capture_output=True)
            subprocess.run([bcftools, "index", str(bcf)], check=True, capture_output=True)
            run = root / "run"
            argv = ["--phase", "all", "--root", str(run), "--input-bcf", str(bcf), "--bcftools", bcftools]
            sv.main(argv)
            summary = json.loads((run / "compare" / "summary.json").read_text())
            self.assertEqual([(r["type"], r["original"], r["new"]) for r in summary["counts"]],
                             [("ALL", 11, 9), ("INS", 7, 6), ("DEL", 3, 2), ("OTHER", 1, 1)])
            self.assertEqual(summary["clusters_with_same_sample_cocarriage"], 1)
            self.assertEqual(summary["clusters_with_lost_carriers"], 0)
            prep = json.loads((run / "prepare" / "summary.json").read_text())
            self.assertEqual(prep["source_rows_with_sv"], 10)
            self.assertEqual(prep["original_sv_alleles"], 11)
            original = sv.load_catalogue(run / "prepare" / "original.vcf.gz")
            self.assertEqual([r["ac"] for r in original.values() if r["pos"] == 20000], [1, 1])
            checksum = sv.common.digest(run / "merge" / "new.vcf.gz")
            mtime = (run / "merge" / "new.vcf.gz").stat().st_mtime_ns
            sv.main(argv)
            self.assertEqual(checksum, sv.common.digest(run / "merge" / "new.vcf.gz"))
            self.assertEqual(mtime, (run / "merge" / "new.vcf.gz").stat().st_mtime_ns)


if __name__ == "__main__":
    unittest.main()
