from __future__ import annotations

import csv
import os
import subprocess
from pathlib import Path

import pytest


pytestmark = pytest.mark.skipif(
    "GAMMA_SMC_BIN" not in os.environ,
    reason="compiled Linux binary not available",
)


def test_samples_file_limits_within_pairs_to_requested_bcf_samples(tmp_path):
    fixture_dir = Path(__file__).parent / "data"
    summary = tmp_path / "subset.tsv"
    manifest = tmp_path / "subset.pairs.tsv"
    command = [
        os.environ["GAMMA_SMC_BIN"],
        "--input",
        str(fixture_dir / "sample_subset_panel.vcf"),
        "--input_format",
        "vcf",
        "--samples",
        str(fixture_dir / "sample_subset_samples.txt"),
        "--only_within",
        "--pairs_manifest",
        str(manifest),
        "--recent_summary",
        str(summary),
        "--scaled_mutation_rate",
        "0.00075",
        "--recombination_to_mutation_ratio",
        "0.8",
        "--unscaled_mutation_rate",
        "1.29e-8",
        "--generation_time",
        "25",
        "--recent_threshold_years",
        "4500",
        "--recent_call",
        "mean",
        "--output_at_hets=false",
        "--output_at_stride",
        "1000",
        "--cache_size",
        "10",
        "--threads",
        "1",
    ]
    completed = subprocess.run(command, check=True, capture_output=True, text=True)
    assert "Read 2 samples." in completed.stdout
    with summary.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows and {int(row["n_pairs"]) for row in rows} == {2}
    pair_rows = [
        line
        for line in manifest.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith("#")
    ]
    assert len(pair_rows) == 2
    assert pair_rows[0].split("\t")[2].startswith("B.")
    assert pair_rows[1].split("\t")[2].startswith("D.")


def test_exact_output_positions_file_preserves_requested_order(tmp_path):
    fixture_dir = Path(__file__).parent / "data"
    summary = tmp_path / "positions.tsv"
    positions = tmp_path / "positions.txt"
    positions.write_text("4567\n123\n4567\n", encoding="utf-8")
    command = [
        os.environ["GAMMA_SMC_BIN"],
        "--input",
        str(fixture_dir / "sample_subset_panel.vcf"),
        "--input_format",
        "vcf",
        "--samples",
        str(fixture_dir / "sample_subset_samples.txt"),
        "--only_within",
        "--recent_summary",
        str(summary),
        "--scaled_mutation_rate",
        "0.00075",
        "--recombination_to_mutation_ratio",
        "0.8",
        "--unscaled_mutation_rate",
        "1.29e-8",
        "--generation_time",
        "25",
        "--recent_threshold_years",
        "4500",
        "--recent_call",
        "mean",
        "--output_at_hets=false",
        "--output_at_stride",
        "-1",
        "--output_positions_file",
        str(positions),
        "--cache_size",
        "10",
        "--threads",
        "1",
    ]
    subprocess.run(command, check=True, capture_output=True, text=True)
    with summary.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [int(row["position_0based"]) for row in rows] == [123, 4567]
