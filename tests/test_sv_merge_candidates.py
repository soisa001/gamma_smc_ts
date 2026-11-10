"""Checks for conservative 70% same-position merge screening."""

import importlib.util
from pathlib import Path

import pytest

spec = importlib.util.spec_from_file_location(
    "merge_audit", Path(__file__).resolve().parents[1] / "scripts/audit_sv_merge_candidates.py"
)
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def record(seq, row=1, gt=("0/1", "./."), length=None):
    return audit.parse_record(
        f"chr1\t100\tr{row}\tA\tA{seq}\tINS\t{len(seq) if length is None else length}\t.\t"
        + "\t".join(gt) + "\n", row, 2
    )


def compare(a, b, limit=1000000):
    return audit.compare((a, b, 0.7, 0.7, limit))


def test_trim_common_prefix_suffix():
    assert audit.payload("AC", "AGGGGC", "INS") == ("GGGG", 1, "available")
    assert audit.payload("AGGGGC", "AC", "DEL") == ("GGGG", 1, "available")
    assert audit.payload("N", "NACGT", "INS") == ("ACGT", 1, "available")


def test_threshold_boundary_and_disjoint_carriers():
    a = record("C" * 10)
    b = record("C" * 7, 2, ("./.", "0/1"))
    result = compare(a, b)
    assert result["length_similarity"] == pytest.approx(0.7)
    assert result["sequence_similarity"] == pytest.approx(0.7)
    assert result["status"] == "candidate_disjoint_carriers"
    assert compare(a, record("C" * 6))["status"] == "length_below_threshold"


def test_equal_length_different_sequences_fail():
    result = compare(record("C" * 10), record("G" * 10))
    assert result["status"] == "sequence_below_threshold"
    assert result["sequence_similarity"] is None
    assert not result["similarity_pass"]


def test_symbolic_sequence_never_matches_itself():
    a = audit.parse_record("chr1\t100\tx\tN\t<INS>\tINS\t10001\t.\t0/1\t./.\n", 1, 2)
    result = compare(a, a)
    assert result["length_similarity"] == 1
    assert result["status"] == "sequence_unavailable"
    assert not result["similarity_pass"]


def test_shared_carriers_dosage_and_reference_conflicts():
    a = record("C" * 10)
    b = record("C" * 10, 2, ("1/1", "./."))
    result = compare(a, b)
    assert result["similarity_pass"]
    assert result["status"] == "candidate_genotype_review"
    assert result["shared_carriers"] == 1
    assert result["conflicting_called_genotypes"] == 1
    b = record("C" * 10, 2, ("0/0", "0/1"))
    assert compare(a, b)["status"] == "candidate_genotype_review"


def test_partial_unknown_sequence_and_length_conflicts():
    a = record("C" * 10)
    assert compare(a, record("C" * 10, gt=("./.", "1/.")))["status"] == "candidate_genotype_review"
    assert compare(a, record("C" * 9 + "N"))["status"] == "sequence_unavailable"
    assert compare(a, record("C" * 10, length=11))["status"] == "representation_review"
    assert compare(a, a, limit=9)["status"] == "sequence_length_limit"


def test_pairwise_similarity_is_not_transitive():
    a = record("C" * 10)
    b = record("C" * 7 + "G" * 3)
    c = record("C" * 4 + "G" * 6)
    assert compare(a, b)["similarity_pass"]
    assert compare(b, c)["similarity_pass"]
    assert not compare(a, c)["similarity_pass"]
