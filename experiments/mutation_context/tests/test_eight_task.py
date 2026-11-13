import csv
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).parents[1]))
import benchmark as b
import eight_task as e
import eight_baselines as eb
import grade


def test_eight_kernel_and_all_offsets():
    truth = e.kernel()
    words = b.digits(8)
    assert truth.shape == (65536, 4)
    assert (truth > 0).all() and np.allclose(truth.sum(axis=1), 1)
    assert np.isclose((1-truth[np.arange(len(truth)), words[:, 4]]).mean(), .005)
    assert np.array_equal(truth, e.kernel())
    for pos in range(8):
        changed = words.copy()
        changed[:, pos] = (changed[:, pos]+1)%4
        assert not np.allclose(truth[b.encode(changed)], truth)


def test_eight_circular_oracle_and_integrity(tmp_path):
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
            code = int(b.encode(np.array([seq[(pos-1+d)%len(seq)] for d in range(-4, 4)])))
            assert e.circular_contexts(seq)[pos-1] == code
            writer.writerow([name, pos, *truth[code]])
    score = grade.grade(tmp_path, output, tmp_path/'score.json')
    assert abs(score['oracle_kl_per_mutation']) < 1e-12
    references = eb.run(tmp_path, tmp_path/'references.json')
    assert references['calibration_sites'] == 2000
    assert all(r['sites'] == 200 for r in references['scores'])
    assert all(np.isclose(r['expected_mutations'], score['expected_mutations']) for r in references['scores'])
    prompt = (tmp_path/'public/TASK.md').read_text().lower()
    assert not any(s in prompt for s in ('8-mer', 'eight-base', 'context', 'window', 'flank'))
    (tmp_path/'public/new_parent.fasta').write_text('corruption')
    assert not b.manifest_valid(tmp_path, 'samples.manifest.json')
