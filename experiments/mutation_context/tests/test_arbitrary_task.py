import csv
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).parents[1]))
import arbitrary_task as a
import benchmark as b
import grade
import directional_baselines as db


def test_arbitrary_kernel_has_exact_six_base_span():
    truth = a.kernel()
    assert np.allclose(truth.sum(axis=1), 1) and (truth > 0).all()
    assert np.isclose((1-truth[np.arange(4**7), b.digits(7)[:, 3]]).mean(), .005)
    for position in range(7):
        words = b.digits(7)
        words[:, position] = (words[:, position]+1)%4
        changed = truth[b.encode(words)]
        assert np.array_equal(changed, truth) == (position == 6)
    assert not np.allclose(truth, truth[b.rc_indices(7)][:, b.COMPLEMENT])
    assert np.array_equal(truth, a.kernel())


def test_circular_queries_and_oracle_grading(tmp_path):
    a.prepare(tmp_path, train_bases=2000, test_bases=200, contig_bases=100, workers=1)
    assert a.prepare(tmp_path, train_bases=2000, test_bases=200, contig_bases=100, workers=1)["cached"]
    parents = b.read_fasta(tmp_path/"public"/"new_parent.fasta")
    truth = np.load(tmp_path/"private"/"kernel.npy")
    with (tmp_path/"public"/"queries.tsv").open() as handle:
        queries = list(csv.DictReader(handle, delimiter="\t"))
    pairs = [(q["contig"], int(q["position"])) for q in queries]
    assert len(set(pairs)) == 200 and pairs != sorted(pairs)
    assert set(pairs) == {(name, i) for name in parents for i in range(1, 101)}
    prediction = tmp_path/"oracle.tsv"
    with prediction.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["contig", "position", "p_A", "p_C", "p_G", "p_T"])
        for name, pos in pairs:
            seq = parents[name]
            expected_word = np.array([seq[(pos-1+d)%len(seq)] for d in range(-3, 4)])
            code = int(b.encode(expected_word))
            assert a.circular_contexts(seq)[pos-1] == code
            writer.writerow([name, pos, *truth[code]])
    score = grade.grade(tmp_path, prediction, tmp_path/"grade.json")
    assert abs(score["oracle_kl_per_mutation"]) < 1e-12
    assert abs(score["burden_relative_error"]) < 1e-12
    prompt = (tmp_path/"public"/"TASK.md").read_text().lower()
    assert "circular" in prompt
    for word in ("6-mer", "six-base", "context", "flank", "window"):
        assert word not in prompt
    references = db.run(tmp_path, tmp_path/"references.json")
    assert references["alignment_audit"]["callable_sites"] == 2000
    assert references["calibration_margin"] == 0
    assert all(row["sites"] == 200 for row in references["scores"])
    assert all(np.isclose(row["expected_mutations"], score["expected_mutations"]) for row in references["scores"])
