"""Check the biological meaning and strict cutoff of the lineage diagnostic."""
import importlib.util
from pathlib import Path

import numpy as np
import tskit

spec = importlib.util.spec_from_file_location(
    "evaluate_paused_array", Path(__file__).resolve().parents[1] / "scripts/evaluate_paused_array.py"
)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_lineages_match_pairwise_recency_including_equal_cutoff():
    tables = tskit.TableCollection(1)
    for _ in range(4):
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
    for age in (5, 10, 20):
        tables.nodes.add_row(time=age)
    for parent, child in ((4, 0), (4, 1), (5, 2), (5, 3), (6, 4), (6, 5)):
        tables.edges.add_row(0, 1, parent, child)
    tables.sort()
    tree = tables.tree_sequence().first()
    pairs = np.array([(a, b) for a in range(4) for b in range(a+1, 4)])
    ages = np.array([tree.tmrca(int(a), int(b)) for a, b in pairs])
    for cutoff in (5, 10, 10.001, 20, 20.001):
        groups = module.ancestral_groups(tree, np.arange(4), cutoff)
        np.testing.assert_array_equal(groups[pairs[:, 0]] == groups[pairs[:, 1]], ages < cutoff)
    groups = module.ancestral_groups(tree, np.arange(4), 10)
    assert len(np.unique(groups)) == 3
    assert groups[0] == groups[1]
    assert groups[2] != groups[3]


def test_wilson_includes_boundary_proportions():
    lo, hi = module.wilson(0, 100)
    assert abs(lo) < 1e-14 and 0 < hi < .05
    lo, hi = module.wilson(100, 100)
    assert .95 < lo < 1 and abs(hi-1) < 1e-14
