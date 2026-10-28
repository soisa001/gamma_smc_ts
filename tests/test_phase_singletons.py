"""Synthetic behavior tests; no controlled genomic input is used."""
import csv
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
        row(300, ["0|1:1:10", "0|1:1:11", "0|1:1:12", "0|0:.:13"]),
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
        row(300, ["1|0:10:10", "0|1:20:11", "1|0:30:12", "0|0:.:13"]),
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


def exclusions(directory):
    with gzip.open(directory / "excluded_singletons.tsv.gz", "rt") as handle:
        return {int(r["pos"]): r for r in csv.DictReader(handle, delimiter="\t")}


def test_transfer_both_orientations_and_reference_singleton_preserves_scaffold(tmp_path):
    scaffold, _, manifest = fixture(tmp_path)
    directory, records = run(tmp_path, manifest)
    assert set(records) == {150, 160, 210}
    assert records[150][9] == "0|1:10"
    assert records[160][10] == "0|1:11"
    assert records[210][9] == "1|0:10"
    assert "AC=7" in records[210][7]
    assert "PS" not in records[150][8]
    assert "SP_ANCHORS=0,2" in records[150][7]
    assert "SP_LEFT_POS=100" in records[150][7]
    assert "SP_RIGHT_POS=200" in records[150][7]
    assert "SP_SOURCE_PS=10" in records[150][7]
    assert "SP_SOURCE_GT=1|0" in records[150][7]
    rejected = exclusions(directory)
    assert rejected[170]["reason"] == "CONFLICTING_NEAREST_ANCHORS"
    assert rejected[180]["reason"] == "NO_INFORMATIVE_ANCHOR"
    assert rejected[190]["reason"] == "NO_PHYSICAL_PS"
    audit = ph.audit_preservation(scaffold, directory / "augmented.vcf.gz", "chr1")
    assert audit["preserved_records"] == 5
    assert audit["added_records"] == 3
    summary = json.loads((directory / "audit.json").read_text())
    assert summary["final_counts"] == {"physical": 3, "excluded": 3}


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
    (lambda x: x.replace("Type=String,Description=\"Genotype", "Type=Integer,Description=\"Genotype"), "Incompatible header"),
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
    directory, records = run(tmp_path, manifest, "--min-anchors", "1")
    assert 100 in records
    assert ph.info_value(records[100], "SP_LEFT_POS") == "."
    assert ph.info_value(records[100], "SP_RIGHT_POS") == "200"
    assert ph.audit_preservation(scaffold, directory / "augmented.vcf.gz", "chr1")["added_records"] == 3


def test_two_left_anchors_do_not_substitute_for_right_anchor(tmp_path):
    _, hip, manifest = fixture(tmp_path)
    hip.write_text("\n".join(x for x in hip.read_text().splitlines() if not x.startswith("chr1\t300\t")) + "\n")
    directory, records = run(tmp_path, manifest)
    assert 210 not in records
    assert exclusions(directory)[210]["reason"] == "MISSING_RIGHT_ANCHOR"
    assert exclusions(directory)[210]["left_pos"] == "200"
    _, records = run(tmp_path, manifest, "--min-anchors", "1")
    assert records[210][9] == "1|0:10"


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
    scaffold_rows = [row(pos, ["0|1:1:30"] * len(samples)) for pos in (10, 20, 200)]
    hip_rows = [row(pos, [f"{f}|{1-f}:10:30" for f in flips]) for pos in (10, 20)]
    for i, t in enumerate(truth):
        genotypes = ["0|0:.:30"] * len(samples)
        orientation = t ^ flips[i]
        genotypes[i] = f"{orientation}|{1-orientation}:10:30"
        hip_rows.append(row(30+i, genotypes))
    hip_rows.append(row(200, [f"{f}|{1-f}:10:30" for f in flips]))
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


def test_sv_is_not_used_as_a_snv_anchor(tmp_path):
    scaffold, hip, manifest = fixture(tmp_path)
    # Retain the common variants as homozygotes, leaving only the shared SV as
    # the heterozygous orientation anchor for sample a.
    text = scaffold.read_text().replace("0|1:1:10", "0|0:1:10")
    scaffold.write_text(text)
    directory, records = run(tmp_path, manifest, "--min-anchors", "1")
    assert 150 not in records
    assert exclusions(directory)[150]["reason"] == "NO_INFORMATIVE_ANCHOR"
    assert ph.audit_preservation(scaffold, directory / "augmented.vcf.gz", "chr1")["preserved_records"] == 5


def test_sample_phase_set_isolation(tmp_path):
    scaffold, hip, manifest = fixture(tmp_path)
    # d's singleton gets PS=10, also used by a; a's anchors must not orient d.
    hip.write_text(hip.read_text().replace("0|1:40:13", "0|1:10:13"))
    directory, records = run(tmp_path, manifest)
    assert 180 not in records
    assert exclusions(directory)[180]["same"] == "0"
    assert exclusions(directory)[180]["flipped"] == "0"


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


def insert_record(path, text):
    lines = path.read_text().splitlines(keepends=True)
    meta = [x for x in lines if x.startswith("#")]
    data = sorted([x for x in lines if not x.startswith("#")] + [text], key=lambda x: int(x.split("\t")[1]))
    path.write_text("".join(meta + data))


def test_expand_past_homozygotes_and_do_not_cross_phase_sets(tmp_path):
    scaffold, hip, manifest = fixture(tmp_path)
    for pos in (140, 155):
        # Homozygotes in the same PS carry no orientation information.
        insert_record(scaffold, row(pos, ["1|1:1:30"] * 4))
        insert_record(hip, row(pos, ["1|1:10:30"] * 4))
    for pos in (145, 152):
        # These closer heterozygotes imply the opposite flip but belong to a
        # different physical block, so cannot orient the singleton at 150.
        insert_record(scaffold, row(pos, ["0|1:1:30"] * 4))
        insert_record(hip, row(pos, ["0|1:999:30"] * 4))
    _, records = run(tmp_path, manifest)
    assert records[150][9] == "0|1:10"
    assert ph.info_value(records[150], "SP_LEFT_POS") == "100"
    assert ph.info_value(records[150], "SP_RIGHT_POS") == "200"


def test_distant_disagreement_does_not_override_agreeing_nearest_flanks(tmp_path):
    scaffold, hip, manifest = fixture(tmp_path)
    insert_record(scaffold, row(50, ["0|1:1:30"] * 4))
    insert_record(hip, row(50, ["0|1:10:30"] * 4))
    _, records = run(tmp_path, manifest)
    # Same PS, but a distant conflicting anchor at 50 is outside the chosen pair.
    assert records[150][9] == "0|1:10"
    assert ph.info_value(records[150], "SP_LEFT_POS") == "100"


def test_cannot_skip_conflicting_nearest_marker_for_favorable_farther_marker(tmp_path):
    scaffold, hip, manifest = fixture(tmp_path)
    insert_record(scaffold, row(155, ["0|1:1:30"] * 4))
    insert_record(hip, row(155, ["0|1:10:30"] * 4))
    directory, records = run(tmp_path, manifest, "--min-anchors", "1")
    assert 150 not in records
    rejected = exclusions(directory)[150]
    assert rejected["reason"] == "CONFLICTING_NEAREST_ANCHORS"
    assert rejected["left_pos"] == "100"
    assert rejected["right_pos"] == "155"


def test_same_position_anchor_ties_cannot_hide_conflict(tmp_path):
    scaffold, hip, manifest = fixture(tmp_path)
    insert_record(scaffold, row(100, ["0|1:1:30"] * 4, alt="C"))
    insert_record(hip, row(100, ["0|1:10:30"] * 4, alt="C"))
    directory, records = run(tmp_path, manifest)
    assert 150 not in records
    assert exclusions(directory)[150]["reason"] == "CONFLICTING_NEAREST_ANCHORS"


def test_zero_shared_snvs_excludes_candidates_but_preserves_panel(tmp_path):
    scaffold, hip, manifest = fixture(tmp_path)
    hip.write_text("\n".join(x for x in hip.read_text().splitlines() if x.startswith("#") or int(x.split("\t")[1]) not in (100, 200, 300)) + "\n")
    directory, records = run(tmp_path, manifest)
    assert records == {}
    assert len(exclusions(directory)) == 6
    audit = ph.audit_preservation(scaffold, directory / "augmented.vcf.gz", "chr1")
    assert audit["added_records"] == 0 and audit["preserved_records"] == 5


def test_missing_left_anchor_is_excluded_by_default(tmp_path):
    _, hip, manifest = fixture(tmp_path)
    hip.write_text("\n".join(x for x in hip.read_text().splitlines() if not x.startswith("chr1\t100\t")) + "\n")
    directory, records = run(tmp_path, manifest)
    assert 150 not in records
    assert exclusions(directory)[150]["reason"] == "MISSING_LEFT_ANCHOR"
    _, records = run(tmp_path, manifest, "--min-anchors", "1")
    assert records[150][9] == "0|1:10"


def test_old_algorithm_cache_is_not_reused(tmp_path):
    _, _, manifest = fixture(tmp_path)
    directory, _ = run(tmp_path, manifest)
    for name in ("done.json", "prepare.done.json"):
        path = directory / name
        content = json.loads(path.read_text())
        content["signature"]["algorithm_version"] = 1
        path.write_text(json.dumps(content))
    run(tmp_path, manifest)
    assert not json.loads((tmp_path / "out" / "summary.json").read_text())["results"][0]["cached"]
