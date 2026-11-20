import csv
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).parents[1]))
import benchmark as b
import ten_task as e
import ten_learning_curve as lc
import grade


def test_ten_kernel_and_all_offsets():
    truth = e.kernel()
    words = b.digits(10)
    assert truth.shape == (1048576, 4)
    assert (truth > 0).all() and np.allclose(truth.sum(axis=1), 1)
    assert np.isclose((1-truth[np.arange(len(truth)), words[:, 5]]).mean(), .005)
    assert np.array_equal(truth, e.kernel())
    for pos in range(10):
        changed = words.copy()
        changed[:, pos] = (changed[:, pos]+1)%4
        assert not np.allclose(truth[b.encode(changed)], truth)


def test_ten_circular_oracle_and_integrity(tmp_path):
    options = dict(train_bases=2000, test_bases=200, contig_bases=100, workers=1)
    e.prepare(tmp_path, **options)
    assert e.prepare(tmp_path, **options)['cached']
    parents = b.read_fasta(tmp_path/'public/new_parent.fasta')
    truth = np.load(tmp_path/'private/kernel.npy')
    with (tmp_path/'public/queries.tsv').open() as f:
        queries = list(csv.DictReader(f, delimiter='\t'))
    pairs = [(q['contig'], int(q['position'])) for q in queries]
    assert set(pairs) == {(name, i) for name in parents for i in range(1, 101)}
    assert len(pairs) == len(set(pairs)) and pairs != sorted(pairs)
    output = tmp_path/'oracle.tsv'
    with output.open('w', newline='') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        writer.writerow(['contig', 'position', 'p_A', 'p_C', 'p_G', 'p_T'])
        for name, pos in pairs:
            seq = parents[name]
            code = int(b.encode(np.array([seq[(pos-1+d)%len(seq)] for d in range(-5, 5)])))
            assert e.circular_contexts(seq)[pos-1] == code
            writer.writerow([name, pos, *truth[code]])
    score = grade.grade(tmp_path, output, tmp_path/'score.json')
    assert abs(score['oracle_kl_per_mutation']) < 1e-12
    ids, outcomes, target_truth = lc.target_counts(tmp_path)
    counted_score = b.score_counts(outcomes.sum(axis=1), outcomes, target_truth, target_truth,
                                    ancestral=ids//4**5 % 4)
    assert abs(counted_score['oracle_kl_per_mutation']) < 1e-12
    assert np.isclose(counted_score['expected_mutations'], score['expected_mutations'])
    train, validation = lc.collect(tmp_path/'public', 2000)
    assert train.sum()+validation.sum() == 2000
    cal_parent = b.read_fasta(tmp_path/'public/calibration_parent.fasta')
    cal_child = b.read_fasta(tmp_path/'public/calibration_offspring.fasta')
    for window in lc.WINDOWS:
        counts = lc.aggregate(train+validation, window)
        expected = np.zeros_like(counts)
        left, right = window
        for name, seq in cal_parent.items():
            # Independent modular extraction, not the production rolling encoder.
            words = np.array([[seq[(i+d)%len(seq)] for d in range(-left, right+1)]
                              for i in range(len(seq))])
            np.add.at(expected, (b.encode(words), cal_child[name]), 1)
        assert np.array_equal(counts, expected)
    prompt = (tmp_path/'public/TASK.md').read_text().lower()
    assert not any(s in prompt for s in ('10-mer', 'ten-base', 'context', 'window', 'flank'))
    (tmp_path/'public/new_parent.fasta').write_text('corruption')
    assert not b.manifest_valid(tmp_path, 'samples.manifest.json')


def test_nested_data_and_parallel_generation(tmp_path):
    small, large = tmp_path/'small', tmp_path/'large'
    e.prepare(small, train_bases=200, test_bases=100, contig_bases=100, workers=1)
    e.prepare(large, train_bases=400, test_bases=100, contig_bases=100, workers=2)
    for filename in ('calibration_parent.fasta', 'calibration_offspring.fasta'):
        first = (small/'public'/filename).read_bytes()
        assert (large/'public'/filename).read_bytes().startswith(first)
    for relative in ('public/new_parent.fasta', 'public/queries.tsv', 'private/new_offspring.fasta',
                     'private/kernel.npy', 'private/kernel_encoding.json'):
        assert b.sha256(small/relative) == b.sha256(large/relative)
