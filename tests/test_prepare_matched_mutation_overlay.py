import json
import os
import subprocess
import sys
from pathlib import Path

import msprime
import pytest
import tskit


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "prepare_matched_mutation_overlay.py"


def mutation_free_tables(path: Path) -> tskit.TableCollection:
    tables = tskit.load(path).dump_tables()
    tables.sites.clear()
    tables.mutations.clear()
    tables.provenances.clear()
    return tables


@pytest.mark.parametrize(
    ("tree_name", "vcf_name", "record_file_names"),
    [
        ("selected_s0p05_af30.trees", "selected_s0p05_af30.vcf.gz", False),
        ("selected.trees", "selected.vcf.gz", True),
    ],
)
def test_matched_overlay_preserves_ancestry_and_records_provenance(
    tmp_path,
    tree_name,
    vcf_name,
    record_file_names,
):
    source = tmp_path / "source"
    output = tmp_path / "output"
    source.mkdir()

    ancestry = msprime.sim_ancestry(
        samples=2,
        sequence_length=10_000,
        recombination_rate=1e-6,
        random_seed=17,
    )
    template = msprime.sim_mutations(
        ancestry,
        rate=1e-4,
        model=msprime.BinaryMutationModel(),
        random_seed=23,
    )
    template.dump(source / tree_name)
    metrics = {"mutation_rate": 1.25e-8}
    if record_file_names:
        metrics["selected_tree_file"] = tree_name
        metrics["selected_vcf_file"] = vcf_name
    (source / "metrics.json").write_text(
        json.dumps(metrics) + "\n",
        encoding="utf-8",
    )
    (source / "copied_truth.tsv").write_text("position\tvalue\n0\t1\n")

    environment = os.environ.copy()
    environment["PYTHONPATH"] = os.pathsep.join(
        [str(ROOT / "python"), environment.get("PYTHONPATH", "")]
    )
    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--template-source",
            str(source),
            "--output-dir",
            str(output),
            "--mutation-rate",
            "1.29e-8",
            "--mutation-seed",
            "29",
        ],
        cwd=ROOT,
        env=environment,
        capture_output=True,
        text=True,
        check=True,
    )

    summary = json.loads(completed.stdout)
    metadata = json.loads(
        (output / "mutation_overlay_metadata.json").read_text()
    )
    assert summary["mutation_rate"] == 1.29e-8
    assert metadata["mutation_rate"] == 1.29e-8
    assert metadata["random_seed"] == 29
    assert metadata["sites_and_mutations_cleared"] is True
    assert metadata["mutation_free_tables_equal_ignoring_provenance"] is True
    assert mutation_free_tables(source / tree_name) == mutation_free_tables(
        output / tree_name
    )
    assert (output / vcf_name).is_file()
    assert (output / "copied_truth.tsv").read_text() == "position\tvalue\n0\t1\n"
