"""Synthetic behavior tests; no controlled genomic input is used."""
import gzip
import importlib.util
import json
from pathlib import Path
import random
import subprocess

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "phase_singletons.py"
spec = importlib.util.spec_from_file_location("phase_singletons", SCRIPT)
ph = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ph)


def header(samples=("a", "b", "c", "d")):
    return "\n".join([
        "##fileformat=VCFv4.2", "##contig=<ID=chr1,length=1000000>",
        '##INFO=<ID=AC,Number=A,Type=Integer,Description="Original allele count">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Physical phase set">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples),
    ]) + "\n"


def row(pos, genotypes, alt="G", info=".", ref="A", ident=None):
    return "\t".join(["chr1", str(pos), ident or f"v{pos}", ref, alt, "39", "PASS", info, "GT:PS:DP"] + genotypes) + "\n"


def fixture(tmp_path):
    # a: flipped PS=10; b: consistent PS=20; c: conflicting PS=30.
    scaffold_rows = [
        row(100, ["0|1:1:10", "0|1:1:11", "0|1:1:12", "0|0:.:13"], info="AC=3"),
        row(200, ["0|1:1:10", "0|1:1:11", "0|1:1:12", "0|0:.:13"], info="AC=3"),
        row(230, ["1|0:1:10", "0|0:.:11", "0|0:.:12", "0|0:.:13"], alt="<DEL>", info="END=290", ident="original_sv"),
        # Missing homozygous GT allowed in authoritative original, not in SHAPEIT view.
        row(800, [".|.:.:10", "0|1:1:11", "0|1:1:12", "0|0:.:13"]),
    ]
    hip_rows = [
        row(100, ["1|0:10:10", "0|1:20:11", "0|1:30:12", "0|0:.:13"]),
        row(150, ["1|0:10:10", "0|0:.:11", "0|0:.:12", "0|0:.:13"], info="AC=1"),
        row(160, ["0|0:.:10", "0|1:20:11", "0|0:.:12", "0|0:.:13"]),
        row(170, ["0|0:.:10", "0|0:.:11", "0|1:30:12", "0|0:.:13"]),
        row(180, ["0|0:.:10", "0|0:.:11", "0|0:.:12", "0|1:40:13"]),
        row(190, ["0|0:.:10", "0|0:.:11", "0|0:.:12", "1/0:.:13"]),
        row(200, ["1|0:10:10", "0|1:20:11", "1|0:30:12", "0|0:.:13"]),
        # MAC=1 reference allele, not ALT=1; a is flipped to 1|0.
        row(210, ["0|1:10:10", "1|1:.:11", "1|1:.:12", "1|1:.:13"], info="AC=7"),
        row(220, ["0|1:10:10", "0|1:20:11", "0|0:.:12", "0|0:.:13"]),
        row(230, ["0|1:10:10", "0|0:.:11", "0|0:.:12", "0|0:.:13"], alt="<DEL>", info="END=290"),
        row(240, ["0|1:10:10", "0|0:.:11", "0|0:.:12", "./.:.:13"]),
        row(250, ["0|1:10:10", "0|0:.:11", "0|0:.:12", "0:.:13"]),
        row(260, ["0|2:10:10", "0|0:.:11", "0|0:.:12", "0|0:.:13"], alt="C,G"),
        row(270, ["0|1:10:10", "0|0:.:11", "0|0:.:12", "0|0:.:13"]).replace("\tPASS\t", "\tLowQual\t"),
    ]
    scaffold = tmp_path / "scaffold.vcf"
    hiphase = tmp_path / "hiphase.vcf"
    scaffold.write_text(header() + "".join(scaffold_rows))
    hiphase.write_text(header() + "".join(hip_rows))
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({"genome_build": "GRCh38", "scaffold_phase_scope": "chromosome", "samples": ["a", "b", "c", "d"], "jobs": [{"contig": "chr1", "scaffold": str(scaffold), "hiphase": str(hiphase)}]}))
    return scaffold, hiphase, manifest


def run(tmp_path, manifest, *extra):
    ph.main(["--manifest", str(manifest), "--out", str(tmp_path / "out"), "--workers", "1", *extra])
    directory = tmp_path / "out" / "chr1"
    with gzip.open(directory / "added_singletons.vcf.gz", "rt") as handle:
        records = {int(r[1]): r for r, _ in ph.VCF(handle, "chr1").records()}
    return directory, records


def test_transfer_both_orientations_and_reference_singleton_preserves_scaffold(tmp_path):
    scaffold, _, manifest = fixture(tmp_path)
    directory, records = run(tmp_path, manifest)
    assert set(records) == {150, 160, 170, 180, 190, 210}
    assert records[150][9] == "0|1:10"
    assert records[160][10] == "0|1:11"
    assert records[210][9] == "1|0:10"
    assert "AC=7" in records[210][7]
    assert "PS" not in records[150][8]
    assert "SP_ANCHORS=0,3" in records[150][7]
    assert "SP_SOURCE_PS=10" in records[150][7]
    assert "SP_SOURCE_GT=1|0" in records[150][7]
    assert records[170][11].startswith("0/1:")
    assert "CONFLICTING_ANCHORS" in records[170][7]
    assert "INSUFFICIENT_ANCHORS" in records[180][7]
    assert "NO_PHYSICAL_PS" in records[190][7]
    audit = ph.audit_preservation(scaffold, directory / "augmented.vcf.gz", "chr1")
    assert audit["preserved_records"] == 4
    assert audit["added_records"] == 6
    summary = json.loads((directory / "audit.json").read_text())
    assert summary["final_counts"] == {"physical": 3, "unresolved": 3}


def test_resume_checks_output_content_and_repairs_corruption(tmp_path):
    _, _, manifest = fixture(tmp_path)
    directory, _ = run(tmp_path, manifest)
    original = (directory / "augmented.vcf.gz").stat().st_mtime_ns
    run(tmp_path, manifest)
    assert (directory / "augmented.vcf.gz").stat().st_mtime_ns == original
    assert json.loads((tmp_path / "out" / "summary.json").read_text())["results"][0]["cached"]
    (directory / "augmented.vcf.gz").write_bytes(b"corrupt")
    run(tmp_path, manifest)
    assert not json.loads((tmp_path / "out" / "summary.json").read_text())["results"][0]["cached"]
    with gzip.open(directory / "augmented.vcf.gz", "rt") as handle:
        assert "original_sv" in handle.read()


@pytest.mark.parametrize("mutation,match", [
    (lambda x: x.replace("\ta\tb\tc\td\n", "\tb\ta\tc\td\n"), "Sample identity"),
    (lambda x: x + x.splitlines(keepends=True)[-1], "Duplicate"),
    (lambda x: x.replace("chr1\t200\t", "chr1\t99\t"), "not sorted"),
    (lambda x: x.replace("chr1\t", "chr2\t"), "shard or indexed"),
    (lambda x: x.replace("Type=Integer,Description=\"Depth", "Type=Float,Description=\"Depth"), "Incompatible header"),
    (lambda x: x.replace("1|0:10:10", "2|0:10:10"), "outside REF/ALT"),
])
def test_invalid_inputs_fail_closed(tmp_path, mutation, match):
    _, hip, manifest = fixture(tmp_path)
    hip.write_text(mutation(hip.read_text()))
    with pytest.raises(ValueError, match=match):
        run(tmp_path, manifest)
    assert not (tmp_path / "out" / "chr1" / "done.json").exists()


def test_unphased_scaffold_is_rejected(tmp_path):
    scaffold, _, manifest = fixture(tmp_path)
    scaffold.write_text(scaffold.read_text().replace("0|1:1:10", "0/1:1:10"))
    with pytest.raises(ValueError, match="Unphased scaffold"):
        run(tmp_path, manifest)


def test_audit_detects_original_sv_metadata_mutation(tmp_path):
    scaffold, _, manifest = fixture(tmp_path)
    directory, _ = run(tmp_path, manifest)
    output = directory / "augmented.vcf.gz"
    with gzip.open(output, "rt") as handle:
        content = handle.read()
    with gzip.open(output, "wt") as handle:
        handle.write(content.replace("END=290", "END=291"))
    with pytest.raises(ValueError, match="record changed"):
        ph.audit_preservation(scaffold, output, "chr1")


def test_existing_singleton_scaffold_record_not_replaced(tmp_path):
    scaffold, _, manifest = fixture(tmp_path)
    text = scaffold.read_text()
    line = row(150, ["1|0:999:99", "0|0:.:99", "0|0:.:99", "0|0:.:99"], info="AC=1", ident="keep_me")
    scaffold.write_text(text.replace("chr1\t200\t", line + "chr1\t200\t"))
    directory, records = run(tmp_path, manifest)
    assert 150 not in records
    with gzip.open(directory / "augmented.vcf.gz", "rt") as handle:
        assert line in handle.read()


def test_distinct_variant_at_same_position_preserves_order(tmp_path):
    scaffold, hip, manifest = fixture(tmp_path)
    hip.write_text(hip.read_text().replace("chr1\t150\t", "chr1\t100\t").replace("100\tv150\tA\tG", "100\tv150\tA\tC"))
    directory, records = run(tmp_path, manifest)
    assert 100 in records
    assert ph.audit_preservation(scaffold, directory / "augmented.vcf.gz", "chr1")["added_records"] == 6


def test_sparse_anchor_threshold(tmp_path):
    _, _, manifest = fixture(tmp_path)
    _, records = run(tmp_path, manifest, "--min-anchors", "4")
    assert all(ph.info_value(r, "SP_METHOD") == "UNRESOLVED" for r in records.values())


def test_global_scaffold_assertion_required(tmp_path):
    _, _, manifest = fixture(tmp_path)
    data = json.loads(manifest.read_text())
    del data["scaffold_phase_scope"]
    manifest.write_text(json.dumps(data))
    with pytest.raises(SystemExit):
        run(tmp_path, manifest)


def test_random_orientations_for_many_samples(tmp_path):
    rng = random.Random(ph.SEED)
    samples = [f"s{i}" for i in range(100)]
    flips = [rng.randrange(2) for _ in samples]
    truth = [rng.randrange(2) for _ in samples]
    scaffold_rows = [row(pos, ["0|1:1:30"] * len(samples)) for pos in (10, 20)]
    hip_rows = [row(pos, [f"{f}|{1-f}:10:30" for f in flips]) for pos in (10, 20)]
    for i, t in enumerate(truth):
        genotypes = ["0|0:.:30"] * len(samples)
        orientation = t ^ flips[i]
        genotypes[i] = f"{orientation}|{1-orientation}:10:30"
        hip_rows.append(row(30+i, genotypes))
    sf, hf = tmp_path / "scaffold.vcf", tmp_path / "hip.vcf"
    sf.write_text(header(samples) + "".join(scaffold_rows))
    hf.write_text(header(samples) + "".join(hip_rows))
    manifest = tmp_path / "m.json"
    manifest.write_text(json.dumps({"genome_build": "synthetic", "scaffold_phase_scope": "chromosome", "samples": samples, "jobs": [{"contig": "chr1", "scaffold": str(sf), "hiphase": str(hf)}]}))
    _, records = run(tmp_path, manifest)
    for i, t in enumerate(truth):
        assert records[30+i][9+i] == f"{t}|{1-t}:30"


def test_cli_subprocess(tmp_path):
    _, _, manifest = fixture(tmp_path)
    result = subprocess.run([__import__("sys").executable, str(SCRIPT), "--manifest", str(manifest), "--out", str(tmp_path / "cli"), "--workers", "1"], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "cli" / "chr1" / "done.json").exists()


def test_sv_only_anchor_lookup_and_incompatible_end(tmp_path):
    scaffold, hip, manifest = fixture(tmp_path)
    # Retain the common variants as homozygotes, leaving only the shared SV as
    # the heterozygous orientation anchor for sample a.
    text = scaffold.read_text().replace("0|1:1:10", "0|0:1:10")
    scaffold.write_text(text)
    directory, records = run(tmp_path, manifest, "--min-anchors", "1")
    assert records[150][9].startswith("0|1:")
    assert "SP_ANCHORS=0,1" in records[150][7]
    # An identically labelled symbolic DEL with a different end is not an anchor.
    hip.write_text(hip.read_text().replace("END=290", "END=295"))
    _, records = run(tmp_path, manifest, "--min-anchors", "1")
    assert records[150][9].startswith("0/1:")
    assert ph.info_value(records[150], "SP_METHOD") == "UNRESOLVED"


def test_sample_phase_set_isolation(tmp_path):
    scaffold, hip, manifest = fixture(tmp_path)
    # d's singleton gets PS=10, also used by a; a's anchors must not orient d.
    hip.write_text(hip.read_text().replace("0|1:40:13", "0|1:10:13"))
    _, records = run(tmp_path, manifest)
    assert records[180][12].startswith("0/1:")
    assert "SP_ANCHORS=0,0" in records[180][7]


def test_manifest_generator_checks_samples_without_emitting_ids(tmp_path):
    scaffold, hip, _ = fixture(tmp_path)
    output = tmp_path / "generated.json"
    proc = subprocess.run([__import__("sys").executable, str(SCRIPT.with_name("make_singleton_manifest.py")),
                           "--scaffold", str(scaffold), "--hiphase", str(hip), "--contigs", "chr1", "--genome-build", "GRCh38",
                           "--assert-chromosome-phase", "--output", str(output)], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    assert json.loads(output.read_text())["samples"] == ["a", "b", "c", "d"]
    assert "4 ordered samples" in proc.stdout


def test_truncated_gzip_rejected(tmp_path):
    _, hip, manifest = fixture(tmp_path)
    compressed = hip.with_suffix(".vcf.gz")
    compressed.write_bytes(gzip.compress(hip.read_bytes())[:-8])
    data = json.loads(manifest.read_text())
    data["jobs"][0]["hiphase"] = str(compressed)
    manifest.write_text(json.dumps(data))
    with pytest.raises(EOFError):
        run(tmp_path, manifest)


def test_new_code_commit_is_not_a_cache_key_and_worker_change_reuses_results(tmp_path):
    _, _, manifest = fixture(tmp_path)
    directory, _ = run(tmp_path, manifest)
    signature = json.loads((directory / "done.json").read_text())["signature"]
    assert "git_commit" not in signature
    ph.main(["--manifest", str(manifest), "--out", str(tmp_path / "out"), "--workers", "4"])
    assert json.loads((tmp_path / "out" / "summary.json").read_text())["results"][0]["cached"]


def test_failed_rerun_invalidates_stale_completion_marker(tmp_path):
    scaffold, _, manifest = fixture(tmp_path)
    directory, _ = run(tmp_path, manifest)
    scaffold.write_text(scaffold.read_text().replace("0|1:1:10", "0/1:1:10"))
    with pytest.raises(ValueError, match="Unphased scaffold"):
        run(tmp_path, manifest)
    assert not (directory / "done.json").exists()
    assert not (tmp_path / "out" / "summary.json").exists()
