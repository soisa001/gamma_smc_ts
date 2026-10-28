"""Scientific contracts of the small phylogenetic benchmark."""

import importlib.util
import io
from pathlib import Path

import numpy as np
import pytest
import tskit
from Bio import Phylo

spec = importlib.util.spec_from_file_location(
    "phylo", Path(__file__).parents[1] / "scripts/sv_phylogeny_benchmark.py"
)
phylo = importlib.util.module_from_spec(spec)
spec.loader.exec_module(phylo)


def test_branch_placement_and_nonmonophyly():
    tree = Phylo.read(io.StringIO("(((a:10,b:10):90,c:100):900,chimp:1000);"), "newick")
    assert phylo.branch_bounds(tree, ["a", "b"], 1) == (10, 100)
    assert phylo.branch_bounds(tree, ["a", "c"], 1) is None
    result = phylo.summarize_bounds(([(10, 100)], 1), {"age_years": 75})
    assert result["age_contained"]
    assert result["estimate_years"] == 55


def test_chimp_preserves_human_genealogy():
    tables = tskit.TableCollection(100)
    tables.time_units = "generations"
    for _ in range(3):
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
    tables.nodes.add_row(time=10)
    tables.nodes.add_row(time=100)
    for p, c in [(3, 0), (3, 1), (4, 3), (4, 2)]:
        tables.edges.add_row(0, 100, p, c)
    tables.sort()
    human = tables.tree_sequence()
    cfg = dict(
        seed=20380101,
        chimp_split_years=6000000,
        generation_years=25,
        human_chimp_ancestor_diploid_ne=50000,
        recombination_rate=1e-8,
        mutation_rate=1.25e-8,
    )
    augmented = phylo.add_chimp(human, cfg, "test")
    for tree in augmented.trees():
        assert tree.tmrca(0, 1) == 10
        assert tree.tmrca(0, 2) == 100
        assert tree.tmrca(0, 3) >= 240000


def test_receipt_detects_corruption(tmp_path):
    (tmp_path / "file").write_text("complete")
    phylo.receipt(tmp_path, "phase", {"sample": "x"}, ["file"])
    assert phylo.valid_receipt(tmp_path, "phase", {"sample": "x"})
    (tmp_path / "file").write_text("corrupt")
    assert not phylo.valid_receipt(tmp_path, "phase", {"sample": "x"})


def test_nonclock_tree_rejected():
    tree = Phylo.read(io.StringIO("((a:10,b:10):10,c:100);"), "newick")
    with pytest.raises(AssertionError, match="time tree"):
        phylo.branch_bounds(tree, ["a", "b"])


def test_convergence_detects_stuck_or_shifted_chains():
    rng = np.random.default_rng(20380101)
    a, b = rng.normal(size=(2, 2000))
    assert phylo.ess(a) > 200
    assert phylo.split_rhat([a, b]) < 1.05
    assert phylo.split_rhat([a, b + 10]) > 1.05
    assert phylo.ess(np.ones(1000)) == 0


def test_lsd_dated_units_must_match_root_annotation(tmp_path):
    path = tmp_path / "dated.nex"
    path.write_text(
        '#NEXUS\nBegin trees;\ntree 1 = (a:100,b:100)[&date="-100"];\nEnd;\n'
    )
    assert len(phylo.read_lsd_time_tree(path).get_terminals()) == 2
    path.write_text(
        '#NEXUS\nBegin trees;\ntree 1 = (a:0.001,b:0.001)[&date="-100"];\nEnd;\n'
    )
    with pytest.raises(AssertionError, match="units disagree"):
        phylo.read_lsd_time_tree(path)


def test_completed_longer_chain_is_reused_only_for_same_model():
    requested = '<run chainLength="2000000" rate="1e-8"/>'
    assert phylo.at_least_requested_chain(
        requested.replace("2000000", "4000000"), requested
    )
    assert not phylo.at_least_requested_chain(
        requested.replace("2000000", "1000000"), requested
    )
    assert not phylo.at_least_requested_chain(
        requested.replace("1e-8", "2e-8"), requested
    )
