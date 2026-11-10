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


def make_fixture(path):
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
        # Symbolic INS size ratio 0.8: merge only at 70%, sequence not claimed.
        add(200000, "INS", 60000, "N", "<INS>", "s1")
        add(200000, "INS", 75000, "N", "<INS>", "s2")
        # Symbolic DELs within both start/end distance tolerances.
        add(300000, "DEL", 60000, "N", "<DEL>", "s1")
        add(300000, "DEL", 60300, "N", "<DEL>", "s2")
        # A DUP cannot merge with DEL at the same position/size.
        add(300000, "DUP", 60000, "N", "<DUP>", "s3")
        # Ambiguous literal sequence is passed through separately.
        add(400000, "INS", 60000, "A", "A" + "N"*60000, "s4")


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
    assert rows[0.9, "INS", "symbolic_size_position"]["collapsed_records"] == "0"
    assert rows[0.7, "INS", "symbolic_size_position"]["collapsed_records"] == "1"
    assert rows[0.9, "DEL", "symbolic_size_position"]["collapsed_records"] == "1"
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
