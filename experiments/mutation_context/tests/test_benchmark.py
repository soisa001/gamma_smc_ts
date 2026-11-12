import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pytest

spec = importlib.util.spec_from_file_location("benchmark", Path(__file__).parents[1] / "benchmark.py")
b = importlib.util.module_from_spec(spec)
spec.loader.exec_module(b)
sys.modules["benchmark"] = b


def test_kernel_has_exact_requested_dependency():
    for scenario, k in (("constant", 1), ("three", 3), ("five_shift", 5)):
        kernel = b.make_kernel(scenario, .005)
        assert np.allclose(kernel.sum(axis=1), 1)
        assert (kernel > 0).all()
        assert np.allclose(kernel, kernel[b.rc_indices(7)][:, b.COMPLEMENT])
        groups = b.project(np.arange(4**7), k)
        for index in range(4**k):
            assert np.allclose(kernel[groups == index], kernel[groups == index][0])
    five = b.make_kernel("five_shift", .005)
    assert not np.allclose(five[b.project(np.arange(4**7), 3) == 0], five[0])


def test_global_alignment_recovers_snv_and_masks_indel():
    rng = np.random.default_rng(33)
    parent = b.sequence(rng, 500, .25)
    child = parent.copy()
    child[101] = (child[101] + 1) % 4
    calls, valid, distance = b.aligned_calls(parent, child)
    assert distance == 1
    assert np.array_equal(calls[valid], child[valid])
    child = np.insert(child, 250, [0, 1, 2])
    calls, valid, distance = b.aligned_calls(parent, child)
    assert distance == 4
    assert valid[101] and calls[101] != parent[101]
    assert not valid[250]
    assert np.array_equal(calls[valid & (np.arange(len(parent)) > 260)], parent[valid & (np.arange(len(parent)) > 260)])


def test_sequence_changes_outer_flanks_without_trimer_shift():
    rng = np.random.default_rng(7)
    seq = b.sequence(rng, 400_000, .8)
    assert abs((seq[:-4] == seq[4:]).mean() - .8) < .005
    tri = b.project(b.contexts(seq), 3)
    freq = np.bincount(tri, minlength=64) / len(tri)
    assert np.abs(freq - 1/64).max() < .002


def test_probability_scores_oracle_zero_and_worse_misspecification():
    truth = b.make_kernel("five_shift", .005)
    opp = np.ones(4**7)
    oracle = b.score_counts(opp, truth, truth, truth)
    assert abs(oracle["oracle_kl_per_mutation"]) < 1e-12
    assert abs(oracle["oracle_conditional_kl"]) < 1e-12
    assert abs(oracle["burden_relative_error"]) < 1e-12
    poor = b.score_counts(opp, truth, truth, b.make_kernel("three", .005))
    assert poor["oracle_kl_per_mutation"] > .01


def test_fitting_and_cache_integrity(tmp_path):
    truth = b.make_kernel("five_shift", .005)
    # Fractional exact expectation isolates fitting from Monte Carlo/alignment.
    train = truth * 10000
    validation = truth * 2500
    models, shrink, losses = b.fit_models(train, validation)
    assert losses["5mer"] < losses["3mer"]
    assert np.allclose(models["7mer"], models["5mer"][b.project(np.arange(4**7), 5)], atol=1e-5)
    for model in models.values():
        assert np.allclose(model.sum(axis=1), 1)
        assert (model > 0).all()
    path = tmp_path / "data.txt"
    path.write_text("correct")
    b.seal(tmp_path, "manifest.json", {"seed": 0}, [path])
    assert b.manifest_valid(tmp_path, "manifest.json", {"seed": 0})
    path.write_text("corrupt")
    assert not b.manifest_valid(tmp_path, "manifest.json")


def test_task_prompt_does_not_reveal_context():
    text = b.TASK_PROMPT.lower()
    for phrase in ("context", "5-mer", "5mer", "flank", "trinucleotide", "pentanucleotide"):
        assert phrase not in text


def test_public_only_solver_matches_evaluator_and_rejects_bad_submission(tmp_path):
    config = {"train_bases": 20000, "test_bases": 1000, "contig_bases": 1000,
              "mean_rate": .005, "target_outer_match": .8}
    job = (str(tmp_path), "five_shift", 0, config)
    b.prepare_one(job)
    assert b.prepare_one(job).startswith("cached")
    b.fit_one(job)
    b.evaluate_one(job)
    task = tmp_path / "five_shift" / "seed_000"
    modules = {}
    for name in ("solve", "grade", "release"):
        spec = importlib.util.spec_from_file_location(name, Path(__file__).parents[1] / f"{name}.py")
        modules[name] = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(modules[name])
    prediction = tmp_path / "prediction.tsv"
    modules["solve"].solve(task / "public", prediction, "5mer")
    scored = modules["grade"].grade(task, prediction, tmp_path / "grade.json")
    expected = next(row for row in json.loads((task / "scores.json").read_text()) if row["method"] == "5mer")
    for metric in ("oracle_kl_per_mutation", "observed_nll_per_site", "burden_relative_error"):
        assert scored[metric] == pytest.approx(expected[metric], abs=1e-8)
    lines = prediction.read_text().splitlines()
    bad = tmp_path / "bad.tsv"
    bad.write_text("\n".join(lines[:-1]) + "\n")
    with pytest.raises(ValueError, match="every query"):
        modules["grade"].grade(task, bad, tmp_path / "bad_grade.json")
    fields = lines[1].split("\t")
    fields[2:] = ["0.4"]*4
    lines[1] = "\t".join(fields)
    bad.write_text("\n".join(lines)+"\n")
    with pytest.raises(ValueError, match="sum to one"):
        modules["grade"].grade(task, bad, tmp_path / "bad_grade.json")
    modules["release"].release(task, tmp_path / "task.zip")
    released_hash = b.sha256(tmp_path / "task.zip")
    modules["release"].release(task, tmp_path / "task.zip")
    assert b.sha256(tmp_path / "task.zip") == released_hash
    # No hidden files are required by the standalone reference solver.
    (task / "private").rename(task / "withheld_elsewhere")
    modules["solve"].solve(task / "public", tmp_path / "blind.tsv", "5mer")
    assert b.sha256(tmp_path / "blind.tsv") == b.sha256(prediction)
    selected = modules["solve"].solve(task / "public", tmp_path / "auto.tsv", "auto")
    assert selected["selected_method"] in b.METHODS
