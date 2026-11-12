from pathlib import Path
import sys

import numpy as np

root = Path(__file__).parents[1]
sys.path.insert(0, str(root))
import asymmetric_task as a
import benchmark as b
import directional_baselines as db


def test_exact_directional_six_base_dependence():
    p = a.kernel()
    assert np.allclose(p.sum(axis=1), 1)
    assert (p > 0).all()
    # Changing +3 cannot change the process; every other input base can.
    for position in range(7):
        words = b.digits(7)
        words[:, position] = (words[:, position]+1) % 4
        changed = p[b.encode(words)]
        if position == 6:
            assert np.array_equal(changed, p)
        else:
            assert np.max(np.abs(changed-p)) > 1e-4
    assert not np.allclose(p, p[b.rc_indices(7)][:, b.COMPLEMENT])


def test_six_base_target_preserves_shorter_composition_and_cache(tmp_path):
    rng = np.random.default_rng(0)
    seq = a.sequence(rng, 500_000, .8)
    assert abs((seq[:-5] == seq[5:]).mean()-.8) < .005
    words = b.project(b.contexts(seq), 5)
    freq = np.bincount(words, minlength=4**5)/len(words)
    assert np.abs(freq-1/1024).sum()/2 < .04
    a.prepare(tmp_path, train_bases=20000, test_bases=1000, contig_bases=1000)
    assert a.prepare(tmp_path, train_bases=20000, test_bases=1000, contig_bases=1000)["cached"]
    assert (tmp_path / "public" / "TASK.md").read_text() == b.TASK_PROMPT


def test_directional_projection_and_baseline_evaluation(tmp_path):
    words = b.digits(7)
    for left, right in db.WINDOWS:
        assert np.array_equal(db.project(np.arange(4**7), (left, right)), b.encode(words[:, 3-left:4+right]))
    a.prepare(tmp_path, train_bases=20000, test_bases=1000, contig_bases=1000)
    result = db.run(tmp_path, tmp_path/"reference.json")
    assert len(result["scores"]) == 6
    assert not result["strand_pooling"]
    assert all(np.isfinite(row["oracle_kl_per_mutation"]) for row in result["scores"])
