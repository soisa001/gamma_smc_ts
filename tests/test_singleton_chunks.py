"""Chunk boundaries must not change whole-chromosome phase decisions."""
import gzip
import io
import json
from pathlib import Path
import random
import subprocess
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import phase_singletons as ph
import singleton_chunks as chunks
from test_phase_singletons import fixture, header, row


def execute(manifest, out, count=4, overlap=10, workers=4):
    subprocess.run([sys.executable, str(Path(ph.__file__).resolve()), "--manifest", str(manifest), "--out", str(out),
                    "--workers", str(workers), "--chunks-per-chrom", str(count), "--overlap-bp", str(overlap)], check=True, capture_output=True, text=True)
    result = json.loads((out / "summary.json").read_text())["results"][0]
    with ph.text_input(result["output"]) as handle:
        records = list(ph.VCF(handle, "chr1").records())
    added = {ph.key(r): r for r, _ in records if ph.info_value(r, "SP_METHOD") == "PHYSICAL"}
    return result, added


def shrink_contig(*paths):
    for path in paths:
        path.write_text(path.read_text().replace("length=1000000", "length=1000"))


def test_sparse_reference_minor_boundary_lookup_and_exact_panel(tmp_path):
    scaffold, hiphase, manifest = fixture(tmp_path)
    shrink_contig(scaffold, hiphase)
    result, added = execute(manifest, tmp_path / "out")
    bypos = {k[1]: r for k, r in added.items()}
    assert set(bypos) == {150, 160, 210}
    assert bypos[150][8:] == ["GT", "0|1", "0|0", "0|0", "0|0"]
    assert bypos[210][8:] == ["GT", "1|0", "1|1", "1|1", "1|1"]
    assert "SP_RIGHT_POS=300" in bypos[210][7]  # Outside core 1..250 + halo 10.
    assert result["chunks"][0]["phase"]["cross_chunk_lookups"] >= 1
    assert result["preservation"]["preserved_records"] == 5
    assert result["preservation"]["exact_record_preservation"]
    directory = tmp_path / "out/chr1"
    assert not (directory / "added_singletons.vcf.gz").exists()
    with gzip.open(directory / "added_singletons.jsonl.gz", "rt") as handle:
        sparse = [json.loads(line) for line in handle]
    assert len(sparse) == 3
    assert all(set(r) == {"fields", "carrier", "gt", "background"} and len(r["fields"]) == 8 for r in sparse)


@pytest.mark.parametrize("overlap", [0, 10, 2000])
def test_boundary_ownership_and_distant_conflict_match_serial(tmp_path, overlap):
    scaffold, hiphase, manifest = fixture(tmp_path)
    # 1000-base chromosome; core ends are 250,500,750,1000.
    h = header().replace("length=1000000", "length=1000")
    anchors = [50, 450, 550, 950]
    sf, hf = {}, {}
    for pos in anchors:
        sf[pos] = row(pos, ["0|1:1:10"] * 4)
        hf[pos] = row(pos, [f"{'0|1' if pos == 950 else '1|0'}:42:10"] * 4)
    candidates = [249, 250, 251, 499, 500, 501, 749, 750, 751]
    for pos in candidates:
        hf[pos] = row(pos, ["1|0:42:10", "0|0:.:10", "0|0:.:10", "0|0:.:10"])
    scaffold.write_text(h + "".join(sf[p] for p in sorted(sf)))
    hiphase.write_text(h + "".join(hf[p] for p in sorted(hf)))
    serial, reference = execute(manifest, tmp_path / "serial", count=1)
    result, actual = execute(manifest, tmp_path / "chunked", overlap=overlap)
    assert set(actual) == set(reference)
    assert {k[1] for k in actual} == set(candidates[:6])
    assert result["final_counts"] == serial["final_counts"] == {"physical": 6, "excluded": 3}
    for k in actual:
        assert actual[k][:8] == reference[k][:8]
        assert ph.calls(actual[k])[0][0] == ph.calls(reference[k])[0][0]
    assert [c["extract"].get("singleton_candidates", 0) for c in result["chunks"]] == [2, 3, 3, 1]


def test_random_phase_blocks_homozygotes_and_conflicts_match_serial(tmp_path):
    scaffold, hiphase, manifest = fixture(tmp_path)
    shrink_contig(scaffold, hiphase)
    rng = random.Random(ph.SEED)
    h = header().replace("length=1000000", "length=1000")
    sf, hf = {}, {}
    for pos in range(10, 1000, 20):
        panel, physical = [], []
        for sample in range(4):
            gt = rng.choice(["0|1", "1|0", "0|0"])
            panel.append(f"{gt}:1:20")
            physical.append(f"{gt if rng.randrange(4) else gt[::-1]}:{rng.choice([10, 20])}:20")
        sf[pos] = row(pos, panel)
        hf[pos] = row(pos, physical)
    for pos in range(5, 1000, 10):
        cells = ["0|0:.:20"] * 4
        cells[rng.randrange(4)] = f"{rng.choice(['0|1', '1|0'])}:{rng.choice([10, 20, 30])}:20"
        hf[pos] = row(pos, cells)
    scaffold.write_text(h + "".join(sf[p] for p in sorted(sf)))
    hiphase.write_text(h + "".join(hf[p] for p in sorted(hf)))
    serial, expected = execute(manifest, tmp_path / "serial", count=1)
    result, actual = execute(manifest, tmp_path / "chunked", overlap=0)
    assert result["final_counts"] == serial["final_counts"]
    assert len(actual) > 0
    assert {k: v[:8] for k, v in actual.items()} == {k: v[:8] for k, v in expected.items()}
    for k in actual:
        assert [c[0] for c in ph.calls(actual[k])] == [c[0] for c in ph.calls(expected[k])]
    with gzip.open(tmp_path / "serial/chr1/excluded_singletons.tsv.gz", "rt") as a, gzip.open(tmp_path / "chunked/chr1/excluded_singletons.tsv.gz", "rt") as b:
        assert a.read() == b.read()


def test_chunk_restart_and_corruption_repair(tmp_path):
    scaffold, hiphase, manifest = fixture(tmp_path)
    shrink_contig(scaffold, hiphase)
    out = tmp_path / "out"
    original, expected = execute(manifest, out)
    cached, _ = execute(manifest, out, workers=2)
    assert cached["cached"]
    target = out / "chr1/chunk_000/added_singletons.jsonl.gz"
    target.write_bytes(b"corrupt")
    repaired, actual = execute(manifest, out)
    assert not repaired["cached"]
    assert repaired["chunks"][0]["extract_cached"]
    assert not repaired["chunks"][0]["phase_cached"]
    assert actual == expected
    (out / "chr1/chunk_000/anchors.sqlite").write_bytes(b"bad SQLite")
    rebuilt, actual = execute(manifest, out)
    assert not rebuilt["chunks"][0]["extract_cached"]
    assert actual == expected
    assert original["final_counts"] == rebuilt["final_counts"]


def test_indexed_region_command_uses_positions_not_sv_overlap(tmp_path, monkeypatch):
    path = tmp_path / "test.bcf"
    path.touch()
    commands = []
    class Process:
        def __init__(self, cmd, **kwargs):
            commands.append(cmd)
            self.stdout = io.StringIO(header())
        def __enter__(self):
            return self
        def __exit__(self, *args):
            pass
        def wait(self):
            return 0
    monkeypatch.setattr(ph.subprocess, "Popen", Process)
    with ph.text_input(path, "bcftools", "chr1", [251, 500]) as handle:
        assert ph.VCF(handle, "chr1").samples == ["a", "b", "c", "d"]
    assert commands == [["bcftools", "view", "--no-version", "-Ov", "-r", "chr1:251-500", "--regions-overlap", "0", str(path)]]


def test_missing_contig_length_fails_before_workers(tmp_path):
    scaffold, hiphase, manifest = fixture(tmp_path)
    scaffold.write_text(scaffold.read_text().replace(",length=1000000", ""))
    job = json.loads(manifest.read_text())["jobs"][0]
    from argparse import Namespace
    with pytest.raises(ValueError, match="contig length"):
        chunks.partition(job, Namespace(bcftools="bcftools", chunks_per_chrom=4, overlap_bp=10), ["a", "b", "c", "d"])


def test_same_coordinate_anchor_conflict_across_halo_is_not_double_counted(tmp_path):
    scaffold, hiphase, manifest = fixture(tmp_path)
    shrink_contig(scaffold, hiphase)
    # Two exact allele matches at 300 disagree for carrier a; no search past them.
    sr = row(300, ["0|1:1:10"] * 4, alt="T")
    hr = row(300, ["0|1:10:10", "0|1:20:10", "0|1:30:10", "0|1:40:10"], alt="T")
    lines = scaffold.read_text().splitlines(keepends=True)
    lines.insert(next(i for i, line in enumerate(lines) if line.startswith("chr1\t800\t")), sr)
    scaffold.write_text("".join(lines))
    with hiphase.open("a") as handle:
        handle.write(hr)
    serial, expected = execute(manifest, tmp_path / "serial", count=1)
    result, actual = execute(manifest, tmp_path / "chunked", overlap=60)
    assert {k: r[:8] for k, r in actual.items()} == {k: r[:8] for k, r in expected.items()}
    assert 210 not in {k[1] for k in actual}
    assert result["final_counts"] == serial["final_counts"]
