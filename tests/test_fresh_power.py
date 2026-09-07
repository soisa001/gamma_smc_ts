import io
import json

import msprime
import numpy as np
import pytest
import tskit
import pandas as pd

from gamma_smc_aou.fresh_power import (
    REPO,
    crop_at_site,
    empirical_p,
    focal_variant,
    model,
    ordered_nodes,
    plan,
    seed_for,
    selected_events,
    validate_config,
)
from gamma_smc_aou.tree_sequence import stream_tree_sequence_vcf
from gamma_smc_aou.fresh_power import ARTIFACTS, completed_record, digest


@pytest.fixture
def cfg():
    return json.loads((REPO / "configs/eas_q02_50k.json").read_text())


def test_grid_and_stable_unique_seeds(cfg):
    validate_config(cfg)
    tasks = plan(cfg)
    assert len(tasks) == 2000
    assert sum(t["mode"] == "neutral" for t in tasks) == 1000
    original = [seed_for(cfg, t) for t in tasks]
    assert len(set(original)) == len(original)
    extended = dict(
        cfg,
        selection_coefficients=[0.025] + list(reversed(cfg["selection_coefficients"])),
    )
    assert original == [seed_for(extended, t) for t in tasks]
    assert seed_for(cfg, tasks[0], 1) != seed_for(cfg, tasks[0], 0)


def test_empirical_upper_tail_includes_ties_and_never_zero():
    np.testing.assert_allclose(
        empirical_p([0, 0.2, 0.5, 1], [0.1, 0.2, 0.2, 0.5]), [1, 0.8, 0.4, 0.2]
    )
    assert empirical_p([1], np.zeros(1000))[0] == 1 / 1001
    with pytest.raises(ValueError):
        empirical_p([1], [])


def tiny_tree():
    tables = tskit.TableCollection(110)
    tables.reference_sequence.data = "A" * 110
    ind = tables.individuals.add_row()
    a = tables.nodes.add_row(flags=1, time=0, individual=ind)
    b = tables.nodes.add_row(flags=1, time=0, individual=ind)
    root = tables.nodes.add_row(time=100)
    tables.edges.add_row(0, 110, root, a)
    tables.edges.add_row(0, 110, root, b)
    # Fixed focal allele and nearby polymorphic allele.
    site = tables.sites.add_row(55, "0")
    tables.mutations.add_row(site, root, "1")
    site = tables.sites.add_row(56, "0")
    tables.mutations.add_row(site, a, "1")
    tables.sort()
    return tables.tree_sequence()


def test_fixation_retained_neutral_site_segregates_and_crop_preserves_coordinates():
    ts = tiny_tree()
    nodes = ordered_nodes(ts)
    assert focal_variant(ts, nodes, 55)["af"] == 1
    nearest = focal_variant(ts, nodes)
    assert nearest["position"] == 56 and nearest["af"] == 0.5
    cfg = dict(focal_position_bp=50, scored_length_bp=100)
    cropped, offset = crop_at_site(ts, 56, cfg)
    assert offset == 6 and cropped.sequence_length == 100
    assert cropped.reference_sequence.data == "A" * 100
    assert focal_variant(cropped, ordered_nodes(cropped), 50)["af"] == 0.5


def test_vcf_exact_site_coordinates(tmp_path):
    source = tmp_path / "test.trees"
    tiny_tree().dump(source)
    output = io.StringIO()
    stream_tree_sequence_vcf(source, output, position_transform="one_based")
    rows = [
        line.split("\t")
        for line in output.getvalue().splitlines()
        if not line.startswith("#")
    ]
    assert [int(row[1]) - 1 for row in rows] == [55, 56]


def test_matched_pulse_and_survival_only_events(cfg):
    built = model(cfg)
    pulse = [
        e
        for e in built.model.events
        if isinstance(e, msprime.MassMigration) and e.source == "EAS"
    ]
    assert len(pulse) == 1 and pulse[0].time == 2000 and pulse[0].proportion == 0.02
    events = selected_events(cfg, 0.001)
    conditions = [
        e
        for e in events
        if isinstance(e, __import__("stdpopsim").ConditionOnAlleleFrequency)
    ]
    assert len(conditions) == 1
    assert conditions[0].op == ">" and conditions[0].allele_frequency == 0
    assert events[1].start_time == 2000 and events[1].selection_coeff == 0.001


def test_resume_rejects_corrupted_completed_outputs(tmp_path):
    cfg = dict(
        scored_length_bp=20, stride_bp=10, tmrca_cutoffs_years=[5000], haplotype_pairs=2
    )
    for name in ARTIFACTS:
        (tmp_path / name).write_bytes(b"saved artifact")
    pd.DataFrame({"position_0based": [0, 10], "frac_recent_5000": [0, 0.5]}).to_csv(
        tmp_path / "frac_recent.tsv", sep="\t", index=False
    )
    record = dict(
        fingerprint="expected",
        artifacts={name: dict(sha256=digest(tmp_path / name)) for name in ARTIFACTS},
    )
    (tmp_path / "complete.json").write_text(json.dumps(record))
    assert completed_record(tmp_path, "expected", cfg) == record
    (tmp_path / "simulation.trees").write_bytes(b"truncated")
    with pytest.raises(ValueError, match="corrupt"):
        completed_record(tmp_path, "expected", cfg)
