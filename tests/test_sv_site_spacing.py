"""Behavior checks for record-based spacing and ALT multiplicity."""

import importlib.util
from pathlib import Path
import sys

import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
spec = importlib.util.spec_from_file_location("sv_spacing", SCRIPTS / "plot_sv_site_spacing.py")
spacing = importlib.util.module_from_spec(spec)
spec.loader.exec_module(spacing)
sys.path.remove(str(SCRIPTS))


def test_unsorted_records_multiallelic_and_same_position(tmp_path):
    path = tmp_path / "sites.tsv"
    path.write_text(
        "chr1\t200\t<DEL>\tDEL\t-11000\t11200\t1\n"
        "chr2\t50\t<INS>\tINS\t11000\t.\t1\n"
        "chr1\t100\t<DEL>,<INS>\t.\t-11000,12000\t.\t1,1\n"
        "chr1\t100\t<INS>\tINS\t12000\t.\t1\n"
        "chr1\t400\t<DEL>\tDEL\t-11000\t11400\t1\n"
        "chr2\t150\t<INS>\tINS\t11000\t.\t1\n"
        "chr3\t10\t<DEL>\tDEL\t-11000\t11010\t1\n"
        "chr4\t2\t.\t.\t.\t.\t.\n",
        encoding="utf-8",
    )
    distances, summary = spacing.audit(path)
    assert summary["records"] == 8
    assert summary["multiallelic_records"] == 1
    assert summary["biallelic_records"] == 6
    assert summary["records_without_ALT"] == 1
    assert summary["multiallelic_records_by_group"]["Complex / other"] == 1
    assert sorted(distances["All SVs"]) == [100, 100, 200]
    assert distances["DEL"] == [200]
    assert distances["INS"] == [100]
    assert distances["Complex / other"] == []
    all_summary = summary["spacing"]["All SVs"]
    assert all_summary["adjacent_pairs"] == 4
    assert all_summary["zero_distance_pairs"] == 1
    assert all_summary["positions_with_multiple_records"] == 1
    assert all_summary["records_at_shared_positions"] == 2


def test_only_coincident_starts(tmp_path):
    path = tmp_path / "sites.tsv"
    path.write_text("chr1\t10\t<DEL>\tDEL\t-11000\t11010\t1\n" * 3)
    distances, summary = spacing.audit(path)
    assert distances["All SVs"] == []
    assert summary["multiallelic_records"] == 0
    assert summary["spacing"]["All SVs"]["zero_distance_pairs"] == 2


@pytest.mark.parametrize("text", ["", "chr1\t1\t<DEL>\n", "chr1\t0\t<DEL>\tDEL\t1\t1\t1\n"])
def test_reject_incomplete_or_invalid_input(tmp_path, text):
    path = tmp_path / "sites.tsv"
    path.write_text(text)
    with pytest.raises(ValueError):
        spacing.audit(path)
