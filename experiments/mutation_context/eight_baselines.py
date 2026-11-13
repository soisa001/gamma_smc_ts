"""Author-side window ablations for circular eight-base truth; private scoring last."""
import argparse
from pathlib import Path

import benchmark as b  # Set numeric-library thread limits before importing NumPy.
import numpy as np
from scipy.special import xlogy

from eight_task import circular_contexts


WINDOWS = ((0, 0), (1, 1), (2, 2), (3, 2), (2, 3), (3, 3), (4, 3), (3, 4), (4, 4))
STRENGTHS = np.r_[np.logspace(0, 5, 21), 1e12]


def project(ids, window):
    left, right = window
    return ids // 4**(4-right) % 4**(left+right+1)


def aggregate(counts, window):
    out = np.zeros((4**(sum(window)+1), 4))
    np.add.at(out, project(np.arange(4**9), window), counts)
    return out


def collect(public, group, child_directory=None):
    parents = b.read_fasta(public/f'{group}_parent.fasta')
    children = b.read_fasta((child_directory or public)/f'{group}_offspring.fasta')
    train, validation = np.zeros((4**9, 4)), np.zeros((4**9, 4))
    for index, (name, seq) in enumerate(parents.items()):
        if len(seq) != len(children[name]):
            raise ValueError('Unequal SNV-only sequence lengths')
        ids = circular_contexts(seq, 4, 4)
        counts = np.bincount(ids*4+children[name], minlength=4**9*4).reshape(-1, 4)
        (validation if index % 5 == 0 else train)[:] += counts
    return train, validation


def run(task, output):
    task = Path(task)
    if not b.manifest_valid(task, 'samples.manifest.json'):
        raise ValueError('Corrupt task inputs')
    train, validation = collect(task/'public', 'calibration')
    fits, partial, losses, strengths, parent_windows = {}, {}, {}, {}, {}
    for window in WINDOWS:
        tr, val = aggregate(train, window), aggregate(validation, window)
        full = tr+val
        if window == (0, 0):
            partial[window] = (tr+.5)/(tr.sum(axis=1, keepdims=True)+2)
            fits[window] = (full+.5)/(full.sum(axis=1, keepdims=True)+2)
            strengths[window], parent_windows[window] = 2., None
        else:
            contained = [w for w in fits if w[0] <= window[0] and w[1] <= window[1]]
            # Select the best already-fit contained prior using calibration only.
            parent = min(contained, key=losses.get)
            parent_windows[window] = list(parent)
            idx = np.arange(len(tr)) // 4**(window[1]-parent[1]) % len(fits[parent])
            best = None
            for strength in STRENGTHS:
                pred = (tr+strength*partial[parent][idx])/(tr.sum(axis=1, keepdims=True)+strength)
                loss = -float(xlogy(val, pred).sum())
                if best is None or loss < best[0]:
                    best = loss, strength, pred
            _, strength, partial[window] = best
            strengths[window] = float(strength)
            fits[window] = (full+strength*fits[parent][idx])/(full.sum(axis=1, keepdims=True)+strength)
        losses[window] = -float(xlogy(val, partial[window]).sum())

    # Freeze model choices above; only this section accesses hidden target data.
    target_parts = collect(task/'public', 'new', task/'private')
    target = target_parts[0]+target_parts[1]
    truth8 = np.load(task/'private/kernel.npy')
    truth9 = truth8[np.arange(4**9)//4]
    rows = []
    for window in WINDOWS:
        pred = fits[window][project(np.arange(4**9), window)]
        rows.append({'window': list(window), 'label': f'{sum(window)+1}-mer (-{window[0]}..+{window[1]})',
                     'validation_nll': losses[window], 'shrinkage': strengths[window],
                     'prior_window': parent_windows[window],
                     **b.score_counts(target.sum(axis=1), target, truth9, pred,
                                      ancestral=np.arange(4**9)//4**4 % 4)})
    counts8 = aggregate(train+validation, (4, 3))
    mask = np.ones_like(counts8, dtype=bool)
    mask[np.arange(4**8), np.arange(4**8)//4**3 % 4] = False
    alt = counts8[mask]
    result = {'scores': rows, 'calibration_sites': int(counts8.sum()),
              'calibration_mutations': int(alt.sum()), 'contexts': 4**8,
              'mean_sites_per_context': float(counts8.sum()/4**8),
              'minimum_sites_per_context': int(counts8.sum(axis=1).min()),
              'maximum_sites_per_context': int(counts8.sum(axis=1).max()),
              'mean_mutations_per_context': float(alt.sum()/4**8),
              'mean_events_per_channel': float(alt.mean()),
              'fraction_zero_observation_channels': float((alt == 0).mean()),
              'zero_observation_channels': int((alt == 0).sum()),
              'validation': 'contig index modulo 5 == 0; all circular positions; calibration only',
              'smoothing_grid': STRENGTHS.tolist(), 'strand_pooling': False,
              'sample_manifest_sha256': b.sha256(task/'samples.manifest.json'),
              'source_sha256': b.sha256(__file__)}
    b.write_json(output, result)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('task', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    print(run(args.task, args.output))
