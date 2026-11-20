import csv
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parents[1]))
import platform_grade as g


def fixture(tmp_path, predictions):
    truth = np.array([[.97, .01, .01, .01], [.02, .95, .02, .01]])
    np.savez(tmp_path/'oracle.npz', probabilities=truth, parent=np.array([0, 1]), child=np.array([0, 2]))
    (tmp_path/'queries.tsv').write_text('contig\tposition\na\t2\na\t1\n')
    with (tmp_path/'predictions.tsv').open('w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(['contig', 'position', 'p_A', 'p_C', 'p_G', 'p_T'])
        for i, row in enumerate(predictions):
            w.writerow(['a', 2-i, *row])
    return truth


def test_oracle_and_independent_score(tmp_path):
    p = np.array([[.98, .005, .005, .01], [.01, .97, .01, .01]])
    truth = fixture(tmp_path, p)
    result = g.score(tmp_path/'predictions.tsv', tmp_path/'queries.tsv', tmp_path/'oracle.npz')
    expected = np.sum(truth*np.log(truth/p))/.08
    assert np.isclose(result['oracle_kl_per_mutation'], expected)
    assert np.isclose(result['reward'], np.exp(-expected))
    fixture(tmp_path, truth)
    assert g.score(tmp_path/'predictions.tsv', tmp_path/'queries.tsv', tmp_path/'oracle.npz')['reward'] == 1


@pytest.mark.parametrize('bad', [
    [[1, 0, 0, 0], [.01, .97, .01, .01]],
    [[.98, .005, .005, .02], [.01, .97, .01, .01]],
    [[float('nan'), .01, .01, .01], [.01, .97, .01, .01]],
])
def test_invalid_probabilities(tmp_path, bad):
    fixture(tmp_path, bad)
    with pytest.raises(ValueError):
        g.score(tmp_path/'predictions.tsv', tmp_path/'queries.tsv', tmp_path/'oracle.npz')


def test_query_integrity(tmp_path):
    fixture(tmp_path, [[.97,.01,.01,.01], [.02,.95,.02,.01]])
    p = tmp_path/'predictions.tsv'
    lines = p.read_text().splitlines()
    p.write_text('\n'.join([lines[0], lines[2], lines[1]])+'\n')
    with pytest.raises(ValueError, match='order'):
        g.score(p, tmp_path/'queries.tsv', tmp_path/'oracle.npz')
