"""Integration checks with real Truvari; Linux dependencies from the script."""

import csv
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pysam

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/compare_stringent_recollapse.py"
spec = importlib.util.spec_from_file_location("recollapse", SCRIPT)
recollapse = importlib.util.module_from_spec(spec)
spec.loader.exec_module(recollapse)


def make_fixture(path, same_start_cases=False):
    header = pysam.VariantHeader()
    header.add_line("##fileformat=VCFv4.2")
    # Deliberately missing contig length: the explicit BED must handle it.
    header.contigs.add("chr1")
    for key, number, typ in [("SVTYPE", 1, "String"), ("SVLEN", 1, "Integer"), ("END", 1, "Integer"),
                             ("N_DISCOVERY_SAMPLES", 1, "Integer"), ("NumCollapsed", 1, "Integer")]:
        header.add_line(f'##INFO=<ID={key},Number={number},Type={typ},Description="test">')
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    for name in ("s1", "s2", "s3", "s4"):
        header.add_sample(name)
    with pysam.VariantFile(path, "wb", header=header) as writer:
        def add(pos, kind, length, ref, alt, sample):
            rec = writer.new_record(contig="chr1", start=pos-1,
                                    stop=pos+length if kind != "INS" else pos,
                                    id=f"original_{pos}_{sample}", alleles=(ref, alt))
            rec.info["SVTYPE"] = kind
            rec.info["SVLEN"] = -length if kind == "DEL" else length
            rec.info["N_DISCOVERY_SAMPLES"] = 99
            rec.info["NumCollapsed"] = 777
            for name in header.samples:
                rec.samples[name]["GT"] = (0, 1) if name == sample else (None, None)
            writer.write(rec)
        # Truvari uses edit distance normalized by sum of lengths: 40% changed
        # bases gives ~80% similarity. Only 70% should collapse; both >50 kb.
        add(100, "INS", 60000, "A", "A" + "C"*60000, "s1")
        add(100, "INS", 60000, "A", "A" + "C"*36000 + "G"*24000, "s2")
        # Same carrier must block collapse even for identical alleles.
        add(100000, "INS", 60000, "A", "A" + "T"*60000, "s1")
        add(100000, "INS", 60000, "A", "A" + "T"*60000, "s1")
        # These would merge at 70% using size alone; now both must pass through.
        add(200000, "INS", 60000, "N", "<INS>", "s1")
        add(200000, "INS", 75000, "N", "<INS>", "s2")
        # Symbolic DELs within both start/end distance tolerances.
        add(300000, "DEL", 60000, "N", "<DEL>", "s1")
        add(300000, "DEL", 60300, "N", "<DEL>", "s2")
        # A DUP cannot merge with DEL at the same position/size.
        add(300000, "DUP", 60000, "N", "<DUP>", "s3")
        # Equal-length DELs with ~80% overlap: merge at 70%, not 90%.
        add(600000, "DEL", 60000, "N", "<DEL>", "s1")
        add(612000, "DEL", 60000, "N", "<DEL>", "s2")
        # Equal-length DELs with ~25% overlap: neither run should merge them.
        add(800000, "DEL", 60000, "N", "<DEL>", "s1")
        add(845000, "DEL", 60000, "N", "<DEL>", "s2")
        # Ambiguous literal sequence is passed through separately.
        add(1000000, "INS", 60000, "A", "A" + "N"*60000, "s4")
        if same_start_cases:
            # Identical insertions 1 bp apart must not merge in the new mode.
            add(1100000, "INS", 60000, "A", "A" + "C"*60000, "s1")
            add(1100001, "INS", 60000, "A", "A" + "C"*60000, "s2")
            # Different raw positions become identical after minimal trimming.
            add(1200000, "INS", 60000, "GGA", "GGA" + "C"*60000, "s1")
            add(1200002, "INS", 60000, "A", "A" + "C"*60000, "s2")
            # Identical raw positions become different after trimming.
            add(1300000, "INS", 60000, "GA", "GA" + "C"*60000, "s1")
            add(1300000, "INS", 60000, "A", "A" + "C"*60000, "s2")
            # Valid right anchors are retained without fetching reference bases.
            add(1400000, "INS", 60000, "A", "T"*60000 + "A", "s1")
            add(1400000, "INS", 60000, "A", "T"*60000 + "A", "s2")
            # Shared-start symbolic DEL size/overlap ratio .8: only 70% merges.
            add(1500000, "DEL", 60000, "N", "<DEL>", "s1")
            add(1500000, "DEL", 75000, "N", "<DEL>", "s2")
            # Even >99.99% overlap cannot merge DELs at different starts.
            add(1700000, "DEL", 60000, "N", "<DEL>", "s1")
            add(1700001, "DEL", 60000, "N", "<DEL>", "s2")


def test_actual_collapse_and_cache(tmp_path):
    source = tmp_path / "stringent.bcf"
    make_fixture(source)
    out = tmp_path / "out"
    command = [sys.executable, str(SCRIPT), "--bcf", str(source), "--out", str(out), "--workers", "4"]
    subprocess.run(command, check=True)
    assert json.loads((out / "run_manifest.json").read_text())["versions"]["truvari"] == "5.5.0"
    with (out / "comparison.tsv").open() as handle:
        rows = {(float(r["threshold"]), r["svtype"], r["mode"]): r for r in csv.DictReader(handle, delimiter="\t")}
    assert rows[0.9, "INS", "sequence_resolved"]["collapsed_records"] == "0"
    assert rows[0.7, "INS", "sequence_resolved"]["collapsed_records"] == "1"
    assert rows[0.7, "INS", "sequence_resolved"]["collapsed_gt50kb"] == "1"
    assert all(r["sizemax"] == "-1" for r in rows.values())
    for threshold in (0.9, 0.7):
        assert rows[threshold, "INS", "symbolic_passthrough"]["collapsed_records"] == "0"
        assert rows[threshold, "INS", "symbolic_passthrough"]["retained_records"] == "2"
        assert rows[threshold, "INS", "symbolic_passthrough"]["action"] == "passthrough"
        assert float(rows[threshold, "DEL", "symbolic_size_overlap"]["pctovl"]) == threshold
        assert float(rows[threshold, "DEL", "symbolic_size_overlap"]["pctseq"]) == 0
        folder = out / f"pct{int(threshold*100)}/INS.symbolic_passthrough"
        assert json.loads((folder / "complete.json").read_text())["command"] is None
        with pysam.VariantFile(out / "partitions/INS.symbolic_size_position.vcf.gz") as before:
            with pysam.VariantFile(folder / "retained.vcf.gz") as after:
                assert [str(r) for r in before] == [str(r) for r in after]
    assert rows[0.9, "DEL", "symbolic_size_overlap"]["collapsed_records"] == "1"
    assert rows[0.7, "DEL", "symbolic_size_overlap"]["collapsed_records"] == "2"
    assert rows[0.7, "DUP", "symbolic_size_position"]["collapsed_records"] == "0"
    assert rows[0.7, "INS", "unresolved_passthrough"]["retained_records"] == "1"
    final = out / "pct70/INS.sequence_resolved/retained.vcf.gz"
    with pysam.VariantFile(final) as vcf:
        recs = list(vcf)
        assert sorted(r.info["AC"][0] for r in recs) == [1, 1, 2]
        assert all("N_DISCOVERY_SAMPLES" not in r.info for r in recs)
        assert all(r.info.get("NumCollapsed", 0) != 777 for r in recs)
        assert all(r.info["RC_SOURCE_DISCOVERY"] == 99 for r in recs)
    first = recollapse.digest(final)
    # Old unchanged INS results can lack the new summary columns. They should
    # still be reused, with a consistent schema in the combined comparison TSV.
    receipt = final.parent / "complete.json"
    cached = json.loads(receipt.read_text())
    del cached["summary"]["pctovl"]
    del cached["summary"]["action"]
    receipt.write_text(json.dumps(cached))
    repeated = subprocess.run(command, check=True, text=True, capture_output=True)
    assert "Reusing validated partitions" in repeated.stdout
    assert "Reused pct70/INS.sequence_resolved" in repeated.stdout
    assert recollapse.digest(final) == first


def test_symbolic_span_required_and_ambiguous_sequence_passthrough(tmp_path):
    source = tmp_path / "input.bcf"
    make_fixture(source)
    with pysam.VariantFile(source) as vcf:
        records = list(vcf)
        duplication = next(r for r in records if r.info["SVTYPE"] == "DUP")
        assert recollapse.classify(duplication)[1] == "symbolic_size_position"
        duplication.stop = duplication.pos
        assert recollapse.classify(duplication)[1] == "unresolved_passthrough"
        ambiguous = records[-1]
        assert recollapse.classify(ambiguous)[1] == "unresolved_passthrough"


def test_same_start_native_adapter_and_reuse(tmp_path):
    source = tmp_path / "input.bcf"
    make_fixture(source, same_start_cases=True)
    baseline, out = tmp_path / "baseline", tmp_path / "same_start"
    command = [sys.executable, str(SCRIPT), "--bcf", str(source), "--workers", "4"]
    subprocess.run([*command, "--out", str(baseline)], check=True)
    stamps = {p: p.stat().st_mtime_ns for p in baseline.rglob("*") if p.is_file()}
    same = [*command, "--out", str(out), "--same-start", "--reuse-prepared", str(baseline)]
    subprocess.run(same, check=True)
    assert all(p.stat().st_mtime_ns == stamp for p, stamp in stamps.items())
    with (out / "comparison_vs_baseline.tsv").open() as handle:
        comparison = {(r["threshold"], r["svtype"]): r for r in csv.DictReader(handle, delimiter="\t")}
    assert comparison["0.7", "DEL"]["additional_retained"] == "2"
    for threshold in (90, 70):
        for kind, mode in (("DEL", "symbolic_size_overlap"), ("INS", "sequence_resolved")):
            folder = out / f"pct{threshold}/{kind}.{mode}"
            with pysam.VariantFile(folder / "retained.vcf.gz") as vcf:
                kept = list(vcf)
            with pysam.VariantFile(folder / "removed.vcf.gz") as vcf:
                removed = list(vcf)
            reps = {str(r.info["CollapseId"]): r for r in kept if "CollapseId" in r.info}
            for record in removed:
                match = record.info["MatchId"]
                cluster = match[0] if isinstance(match, tuple) else match.split(",")[0]
                assert (record.chrom, record.pos) == (reps[cluster].chrom, reps[cluster].pos)
            if kind == "DEL":
                assert {1700000, 1700001, 600000, 612000, 800000, 845000} <= {r.pos for r in kept}
                assert len(removed) == (1 if threshold == 90 else 2)
            else:
                assert {1100000, 1100001, 1300000, 1300001} <= {r.pos for r in kept}
                assert sum(r.pos == 1200002 for r in kept) == 1
                assert sum(r.pos == 1400000 for r in kept) == 1
            receipt = json.loads((folder / "complete.json").read_text())
            assert receipt["recipe"]["same_start"] is True
            assert "collapse_same_start.py" in receipt["command"][1]
        folder = out / f"pct{threshold}/INS.symbolic_passthrough"
        receipt = json.loads((folder / "complete.json").read_text())
        assert receipt["summary"]["retained_records"] == 2
        assert receipt["command"] is None
    stamps = {p: p.stat().st_mtime_ns for p in out.rglob("*.vcf.gz")}
    repeat = subprocess.run(same, check=True, capture_output=True, text=True)
    assert "Reusing validated same-start partitions" in repeat.stdout
    assert "Reused pct70/INS.sequence_resolved" in repeat.stdout
    assert all(p.stat().st_mtime_ns == stamp for p, stamp in stamps.items())
