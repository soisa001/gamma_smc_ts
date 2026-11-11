"""Synthetic regional scores, genotype blocks and local native clustering."""
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
import audit_tlr_sv as audit


@unittest.skipUnless(os.environ.get("SV_PILOT_INTEGRATION") == "1", "Requires Linux native tools")
class LocusTests(unittest.TestCase):
    def test_scores_rules_clusters_and_cache(self):
        bcftools = os.environ.get("BCFTOOLS_BIN", "bcftools")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raw = root / "source.vcf"
            header = ('##fileformat=VCFv4.2\n##contig=<ID=chr4>\n'
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
                '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\tc\n')
            gt = [("1|0", "0|0", "0|0"), ("0|0", "0|1", "0|0"), ("0|0", "0|0", "1|0")]
            rows = [(1000, "del100", "A"+"C"*100, "A", 0),
                    (1000, "del80", "A"+"C"*80, "A", 1),
                    (1005, "shifted", "A"+"C"*100, "A", 2),
                    (2000, "ins_a", "A", "A"+"T"*60, 0),
                    (2000, "ins_shared", "A", "A"+"T"*60, 0),
                    (2100, "replacement", "G", "T"*61, 1),
                    (2200, "small", "A", "A"+"G"*49, 1)]
            raw.write_text(header + "".join("\t".join(map(str, ("chr4", pos, name, ref, alt, ".", "PASS", ".", "GT", *gt[sample]))) + "\n"
                                           for pos, name, ref, alt, sample in rows))
            source = root / "source.vcf.gz"
            subprocess.run([bcftools, "view", "-Oz", "-o", str(source), str(raw)], check=True, capture_output=True)
            subprocess.run([bcftools, "index", "-t", str(source)], check=True, capture_output=True)
            run = root / "run"
            args = ["--root", str(run), "--input-bcf", str(source), "--region", "chr4:900-2300", "--flank", "100", "--bcftools", bcftools]
            try:
                audit.main(args)
            except BaseException:
                for log in run.rglob("*.log"):
                    print(log.name, log.read_text()[-4000:])
                raise
            with gzip.open(run / "prepare/alleles.jsonl.gz", "rt") as handle:
                meta = [json.loads(line) for line in handle]
            labels = {r["original_id"]: r["label"] for r in meta}
            with open(run / "audit/pair_scores.tsv") as handle:
                scored = list(csv.DictReader(handle, delimiter="\t"))
            def pair(a, b, threshold=.7, policy="current_r500"):
                return next(r for r in scored if {r["allele_a"], r["allele_b"]} == {labels[a], labels[b]}
                            and float(r["threshold"]) == threshold and r["policy"] == policy)
            self.assertAlmostEqual(float(pair("del100", "del80")["size_similarity"]), .8)
            self.assertAlmostEqual(float(pair("del100", "del80")["payload_edit_similarity_max_length"]), .8)
            self.assertEqual(pair("del100", "del80")["pair_passes"], "True")
            self.assertEqual(pair("del100", "del80", .9)["pair_passes"], "False")
            self.assertIn("size_similarity", pair("del100", "del80", .9)["failed_criteria"])
            self.assertEqual(pair("del100", "shifted")["pair_passes"], "True")
            self.assertIn("different_trimmed_pos", pair("del100", "shifted", policy="same_start")["failed_criteria"])
            self.assertEqual(pair("del100", "shifted", policy="same_start")["same_local_cluster"], "False")
            shared = pair("ins_a", "ins_shared")
            self.assertEqual(float(shared["truvari_sequence_similarity"]), 1)
            self.assertEqual(shared["shared_carriers"], "1")
            self.assertEqual(shared["pair_passes"], "False")
            self.assertEqual(shared["same_local_cluster"], "False")
            self.assertIn("non_simple_passthrough", pair("ins_a", "replacement")["failed_criteria"])
            summary = json.loads((run / "audit/summary.json").read_text())
            self.assertEqual(summary["focus_sv_alleles"], 6)
            self.assertEqual(summary["distinct_pairs"], 15)
            self.assertEqual(len(meta), 7)  # Small indels participated in local collapse.
            saved = {p: p.stat().st_mtime_ns for p in (run / "audit").iterdir()}
            audit.main(args)
            self.assertTrue(all(p.stat().st_mtime_ns == stamp for p, stamp in saved.items()))


if __name__ == "__main__":
    unittest.main()
