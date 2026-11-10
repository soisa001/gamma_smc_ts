"""Representation-aware, BND-free comparison of the four callsets."""

import importlib.util
import json
from pathlib import Path
import sqlite3

import pytest

spec = importlib.util.spec_from_file_location(
    "compare_v3", Path(__file__).resolve().parents[1] / "scripts/compare_v3_sv_callsets.py"
)
compare = importlib.util.module_from_spec(spec)
spec.loader.exec_module(compare)


def allele(alt, kind="INS", length="10001", end=None, ref="A"):
    return compare.describe_allele("chr1", 100, ref, alt, kind, length, end)


def test_multiallelic_and_split_literal_alleles_match():
    multi, flag, _ = compare.parse_line("chr1\t100\tA\tACCC,AGGGG\tINS\t3,4\t.\tPASS\n", {"SVLEN": "A", "SVTYPE": "1"})
    assert flag
    assert [a["length"] for a in multi] == [3, 4]
    for item, alt in zip(multi, ("ACCC", "AGGGG")):
        # Literal identity does not depend on a differing INFO length annotation.
        assert item["key"] == allele(alt)["key"]


def test_symbolic_deletion_end_and_length_matter():
    first = allele("<DEL>", "DEL", "10001", "10101", ref="N")
    other = allele("<DEL>", "DEL", "20000", "20100", ref="N")
    derived = allele("<DEL>", "DEL", None, "10101", ref="N")
    assert first["key"] != other["key"]
    assert first["key"] == derived["key"]
    assert first["evidence"] == "symbolic_metadata"


def test_symbolic_and_literal_insertions_are_not_exact_matches():
    assert allele("<INS>", length="3")["key"] != allele("ACCC", length="3")["key"]
    assert allele("<INS>", length=None)["key"] is None
    assert allele("ANNN", length="3")["evidence"] == "literal_with_N"


@pytest.mark.parametrize("alt,kind", [("N]chr2:500]", "INS"), ("<BND>", "INS"), (".A", ""), ("A.", ""), ("<TRA>", "TRA")])
def test_skip_breakends_even_without_correct_info(alt, kind):
    assert allele(alt, kind)["skip"] == "BND"


def test_ambiguous_multiallelic_scalar_length_not_reused():
    rows, _, _ = compare.parse_line("chr1\t100\tN\t<INS>,<INS>\tINS\t10000\t.\tPASS\n", {"SVLEN": "1"})
    assert all(a["length"] is None and a["key"] is None for a in rows)


def test_threshold_boundary():
    assert compare.length_bin(10000) == "1kb_to_10kb"
    assert compare.length_bin(10001) == "gt_10kb_to_100kb"


def test_overlap_counts_and_subset():
    # One allele in all four, one in main+all, one only in main.
    membership = [(15, "INS", "sequence_resolved", 1), (3, "DEL", "symbolic_metadata", 1), (1, "DEL", "symbolic_metadata", 1)]
    summaries = {label: {"unique_matchable_keys": n, "unmatchable_ALT_alleles": 0}
                 for label, n in zip(compare.LABELS, (3, 2, 1, 1))}
    rows = compare.overlap_rows(membership, summaries)
    assert rows[0]["shared"] == 2
    assert rows[0]["only_a"] == 1
    assert rows[0]["shared_sequence_resolved"] == 1
    assert rows[0]["shared_symbolic_metadata"] == 1
    assert rows[0]["b_subset_of_a"] == "yes"


def test_cache_validates_sample_manifest_and_database(tmp_path):
    db = tmp_path / "alleles.sqlite"
    metadata = tmp_path / "scan.json"
    fingerprint = {"sample_manifest_sha256": "original"}
    saved = {"fingerprint": fingerprint, "summary": {"unique_matchable_keys": 1}}
    with sqlite3.connect(db) as conn:
        conn.execute("CREATE TABLE alleles (key BLOB PRIMARY KEY)")
        conn.execute("INSERT INTO alleles VALUES (?)", (b"test",))
    metadata.write_text(json.dumps(saved))
    assert compare.valid_cache(db, metadata, fingerprint) == saved
    assert compare.valid_cache(db, metadata, {"sample_manifest_sha256": "changed"}) is None
    db.write_bytes(b"corrupt database")
    assert compare.valid_cache(db, metadata, fingerprint) is None


def test_empty_callset_overlap_is_reported():
    summaries = {label: {"unique_matchable_keys": 0, "unmatchable_ALT_alleles": 0} for label in compare.LABELS}
    rows = compare.overlap_rows([], summaries)
    assert len(rows) == 6
    assert rows[0]["shared"] == 0
    assert rows[0]["percent_a_shared"] is None
