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
    from gamma_smc_aou.fresh_power import exact_focal_variant

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
    for position in (55, 56):
        direct = exact_focal_variant(ts, nodes, position)
        scanned = focal_variant(ts, nodes, position)
        assert direct["af"] == scanned["af"]
        np.testing.assert_array_equal(direct["carriers"], scanned["carriers"])
    assert exact_focal_variant(ts, nodes, 57) is None
    assert exact_focal_variant(ts, nodes, 20) is None


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


def test_delayed_selection_keeps_pulse_and_survival_condition(cfg):
    import stdpopsim

    later = dict(cfg, selection_onset_years=10000)
    validate_config(later)
    assert model(later).model.asdict() == model(cfg).model.asdict()
    events = selected_events(later, 0.005)
    fitness = [e for e in events if isinstance(e, stdpopsim.ChangeMutationFitness)]
    condition = [
        e for e in events if isinstance(e, stdpopsim.ConditionOnAlleleFrequency)
    ]
    assert fitness[0].start_time == 400 and fitness[0].end_time == 0
    assert condition[0].start_time == 1995 and condition[0].end_time == 0
    assert condition[0].op == ">" and condition[0].allele_frequency == 0
    for invalid in (0, -1000, 50001):
        with pytest.raises(ValueError, match="Selection must"):
            validate_config(dict(cfg, selection_onset_years=invalid))


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


def test_simulation_only_cli_saves_resumes_and_audits_without_decoding(tmp_path, monkeypatch, cfg):
    from concurrent.futures import ThreadPoolExecutor
    from types import SimpleNamespace
    import gamma_smc_aou.fresh_power as fresh

    settings = dict(cfg, neutral_replicates=1, selected_replicates_per_coefficient=0,
                    selection_coefficients=[], sample_diploids=2, haplotype_pairs=3,
                    simulated_length_bp=11000, scored_length_bp=10000,
                    focal_position_bp=5000, stride_bp=1000, mutation_rate=1e-4,
                    storage_limit_bytes=10**16)
    config = tmp_path / "config.json"
    config.write_text(json.dumps(settings))
    binary = tmp_path / "binary-for-identity-only"
    binary.write_bytes(b"must not execute")
    out = tmp_path / "saved"
    out.mkdir()
    (out / "run_status.json").write_text('{"previous_decoding":"unchanged"}')
    demography = msprime.Demography()
    demography.add_population(name="EAS", initial_size=100)
    monkeypatch.setattr(fresh, "model", lambda _: SimpleNamespace(model=demography))
    monkeypatch.setattr(fresh, "ProcessPoolExecutor", ThreadPoolExecutor)

    def forbidden(*args, **kwargs):
        raise AssertionError("Simulation-only phase attempted decoding/truth profiling or resimulation")

    monkeypatch.setattr(fresh, "run_within_decoder", forbidden)
    monkeypatch.setattr(fresh, "truth_profile", forbidden)
    arguments = ["--config", str(config), "--out", str(out), "--slim", str(binary),
                 "--decoder", str(binary), "--workers", "1"]
    assert fresh.main(arguments + ["--phase", "simulate-only"]) == 0
    directory = out / "neutral/rep0000"
    before = {name: fresh.digest(directory / name) for name in fresh.SIMULATION_ARTIFACTS}
    assert not (directory / "complete.json").exists()
    assert not (directory / "frac_recent.tsv").exists()
    assert not (directory / "truth_frac_recent.tsv").exists()
    assert json.loads((out / "run_status.json").read_text()) == {"previous_decoding": "unchanged"}
    status = json.loads((out / "simulation_status.json").read_text())
    assert status["state"] == "complete" and status["new_simulations"] == 1
    assert status["decoding_requested"] is False
    inventory = pd.read_csv(out / "simulation_inventory.csv")
    assert len(inventory) == 1 and not inventory.already_decoded.any()
    monkeypatch.setattr(fresh, "model", forbidden)
    assert fresh.main(arguments + ["--phase", "simulate-only"]) == 0
    assert fresh.main(arguments + ["--phase", "audit-simulations"]) == 0
    assert before == {name: fresh.digest(directory / name) for name in fresh.SIMULATION_ARTIFACTS}
    # An internally plausible carrier file with an updated checksum must still
    # match the genotypes of the archived raw and cropped trees.
    path = directory / "focal_carriers.npy"
    np.save(path, ~np.load(path), allow_pickle=False)
    receipt = json.loads((directory / "simulated.json").read_text())
    receipt["simulation_hashes"][path.name] = fresh.digest(path)
    (directory / "simulated.json").write_text(json.dumps(receipt))
    with pytest.raises(ValueError, match="allele/carriers"):
        fresh.main(arguments + ["--phase", "audit-simulations"])


def test_simulation_archive_retains_fixed_allele_and_rejects_wrong_identity(tmp_path, cfg):
    import gamma_smc_aou.fresh_power as fresh

    settings = dict(cfg, simulated_length_bp=110, scored_length_bp=100,
                    focal_position_bp=50, sample_diploids=1)
    task = dict(mode="selected", s=0.005, replicate=0)
    raw = tiny_tree()
    cropped, offset = crop_at_site(raw, 55, settings)
    raw.dump(tmp_path / "simulation.trees")
    cropped.dump(tmp_path / "decoded_input.trees")
    np.save(tmp_path / "focal_carriers.npy", np.ones(2, dtype=bool), allow_pickle=False)
    seed = seed_for(settings, task)
    record = dict(task=task, fingerprint="expected", seed=seed,
                  attempts=[dict(accepted=True, attempt=0, seed=seed)],
                  sample_af=1.0, fixed=True, original_focal_position=55, crop_offset=offset,
                  simulation_hashes={name: fresh.digest(tmp_path / name)
                                     for name in fresh.SIMULATION_ARTIFACTS})
    fresh.atomic_json(tmp_path / "simulated.json", record)
    assert fresh.simulated_record(tmp_path, "expected", settings, task)["fixed"]
    with pytest.raises(ValueError, match="task identity"):
        fresh.simulated_record(tmp_path, "expected", settings, dict(task, s=0.001))
    (tmp_path / "simulation.trees").write_bytes(b"truncated")
    with pytest.raises(ValueError, match="corrupt"):
        fresh.simulated_record(tmp_path, "expected", settings, task)


def test_array_extension_preserves_original_manifest_bytes(tmp_path):
    from gamma_smc_aou.fresh_power import versioned_json

    path = tmp_path / "manifest.json"
    original = b'{"selection_coefficients": [0.005]}\n'
    path.write_bytes(original)
    value = dict(selection_coefficients=[0.001, 0.005])
    versioned_json(path, value)
    versions = list((tmp_path / "manifest_history").glob("manifest.*.json"))
    assert len(versions) == 1 and versions[0].read_bytes() == original
    versioned_json(path, value)
    assert list((tmp_path / "manifest_history").glob("manifest.*.json")) == versions
