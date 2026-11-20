"""Author-only nested-data window comparisons; no hidden data used to select fits."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import time

import benchmark as b
import numpy as np
from scipy.special import xlogy

from ten_task import circular_contexts

WINDOWS = ((0, 0), (2, 2), (3, 3), (4, 3), (4, 4), (5, 4), (4, 5), (5, 5))
STRENGTHS = np.r_[np.logspace(0, 5, 21), 1e12]
WIDTH = 11


def project(ids, window):
    left, right = window
    return ids // 4**(5-right) % 4**(left+right+1)


def aggregate(counts, window):
    left, right = window
    shape = (4**(5-left), 4**(left+right+1), 4**(5-right), 4)
    return counts.reshape(shape).sum(axis=(0, 2))


def fasta_records(path):
    name, pieces = None, []
    lookup = np.full(256, -1, dtype=np.int8)
    lookup[np.frombuffer(b'ACGT', dtype=np.uint8)] = np.arange(4)
    with Path(path).open('rb') as handle:
        for line in handle:
            if line.startswith(b'>'):
                if name is not None:
                    seq = lookup[np.frombuffer(b''.join(pieces), dtype=np.uint8)]
                    if (seq < 0).any():
                        raise ValueError('Non-ACGT sequence')
                    yield name, seq
                name, pieces = line[1:].strip().decode(), []
            else:
                pieces.append(line.strip())
        if name is not None:
            seq = lookup[np.frombuffer(b''.join(pieces), dtype=np.uint8)]
            if (seq < 0).any():
                raise ValueError('Non-ACGT sequence')
            yield name, seq


def collect(public, limit):
    train, val = np.zeros((4**WIDTH, 4), dtype=np.int32), np.zeros((4**WIDTH, 4), dtype=np.int32)
    total = 0
    pairs = zip(fasta_records(public/'calibration_parent.fasta'),
                fasta_records(public/'calibration_offspring.fasta'), strict=True)
    for index, ((name, seq), (child_name, child)) in enumerate(pairs):
        if name != child_name or len(seq) != len(child) or total+len(seq) > limit:
            raise ValueError('Homologous complete scaffolds required')
        ids = circular_contexts(seq, 5, 5)
        np.add.at((val if index % 5 == 0 else train).reshape(-1), ids*4+child, 1)
        total += len(seq)
        if total == limit:
            break
    if total != limit:
        raise ValueError('Calibration prefix shorter than requested')
    return train, val


def target_counts(task):
    ids, ys = [], []
    for (name, seq), (child_name, child) in zip(fasta_records(task/'public/new_parent.fasta'),
            fasta_records(task/'private/new_offspring.fasta'), strict=True):
        if name != child_name or len(seq) != len(child):
            raise ValueError('Target correspondence invalid')
        ids.append(circular_contexts(seq, 5, 5))
        ys.append(child)
    unique, inverse = np.unique(np.concatenate(ids), return_inverse=True)
    outcomes = np.bincount(inverse*4+np.concatenate(ys), minlength=len(unique)*4).reshape(-1, 4)
    truth = np.load(task/'private/kernel.npy')[unique//4]
    return unique, outcomes, truth


def fit_one(job):
    task, size, output = job
    task, output = Path(task), Path(output)
    started = time.monotonic()
    train, validation = collect(task/'public', size)
    fits, partial, losses, strengths, parents = {}, {}, {}, {}, {}
    for window in WINDOWS:
        tr, val = aggregate(train, window).astype(float), aggregate(validation, window).astype(float)
        full = tr+val
        if window == (0, 0):
            partial[window] = (tr+.5)/(tr.sum(axis=1, keepdims=True)+2)
            fits[window] = (full+.5)/(full.sum(axis=1, keepdims=True)+2)
            strengths[window], parents[window] = 2., None
        else:
            contained = [w for w in fits if w[0] <= window[0] and w[1] <= window[1]]
            parent = min(contained, key=losses.get)
            parents[window] = list(parent)
            idx = np.arange(len(tr)) // 4**(window[1]-parent[1]) % len(fits[parent])
            prior = partial[parent][idx]
            denom = tr.sum(axis=1, keepdims=True)
            best = None
            observed = val > 0
            weights = val[observed]
            for strength in STRENGTHS:
                pred = (tr+strength*prior)/(denom+strength)
                loss = -float(weights @ np.log(pred[observed]))
                if best is None or loss < best[0]:
                    best = loss, strength, pred
            _, strength, partial[window] = best
            strengths[window] = float(strength)
            fits[window] = (full+strength*fits[parent][idx])/(full.sum(axis=1, keepdims=True)+strength)
        losses[window] = -float(xlogy(val, partial[window]).sum())
        print(f'calibration={size} window={window} validation_nll={losses[window]:.6f}', flush=True)
    # Finish every model selection before accessing private target outcomes/truth.
    ids, outcomes, truth = target_counts(task)
    rows = []
    for window in WINDOWS:
        pred = fits[window][project(ids, window)]
        rows.append({'window': list(window), 'label': f'{sum(window)+1}-mer (-{window[0]}..+{window[1]})',
                     'validation_nll': losses[window], 'shrinkage': strengths[window],
                     'prior_window': parents[window],
                     **b.score_counts(outcomes.sum(axis=1), outcomes, truth, pred,
                                      ancestral=ids//4**5 % 4)})
    counts = aggregate(train+validation, (5, 4))
    exposure = counts.sum(axis=1)
    counts[np.arange(4**10), np.arange(4**10)//4**4 % 4] = 0
    observed_channels = int(np.count_nonzero(counts))
    result = {'calibration_sites': int(exposure.sum()), 'calibration_mutations': int(counts.sum()),
              'mean_sites_per_context': float(exposure.mean()),
              'mean_mutations_per_context': float(counts.sum()/4**10),
              'fraction_unseen_contexts': float((exposure == 0).mean()),
              'fraction_zero_event_channels': 1-observed_channels/(3*4**10),
              'zero_event_channels': 3*4**10-observed_channels,
              'target_fraction_unseen_contexts': float(outcomes.sum(axis=1)[exposure[ids//4] == 0].sum()/outcomes.sum()),
              'scores': rows, 'elapsed_seconds': time.monotonic()-started,
              'validation': 'Scaffold index modulo 5 == 0 in each nested prefix; whole scaffolds',
              'selected_window': list(min(losses, key=losses.get)),
              'source_sha256': b.sha256(__file__)}
    b.write_json(output, result)
    return {'calibration_sites': size, 'output': str(output), 'selected_window': result['selected_window']}


def run(task, output, sizes=(4_000_000, 16_000_000, 64_000_000), workers=4):
    task, output = Path(task), Path(output)
    if not b.manifest_valid(task, 'samples.manifest.json'):
        raise ValueError('Invalid sample manifest')
    if not 1 <= workers <= 4 or any(s <= 0 or s % 50_000 for s in sizes):
        raise ValueError('Use 1..4 workers and positive complete 50,000-base scaffolds')
    cfg = {'method': 'nested_hierarchical_tables_v1', 'sizes': list(sizes),
           'sample_manifest_sha256': b.sha256(task/'samples.manifest.json'),
           'windows': [list(w) for w in WINDOWS], 'smoothing_grid': STRENGTHS.tolist()}
    if b.manifest_valid(output, 'curve.manifest.json', cfg):
        return {'cached': True, 'directory': str(output)}
    output.mkdir(parents=True, exist_ok=True)
    jobs = [(task, size, output/f'calibration_{size}.json') for size in sizes]
    if workers == 1:
        results = list(map(fit_one, jobs))
    else:
        with ProcessPoolExecutor(max_workers=min(workers, len(jobs))) as pool:
            results = list(pool.map(fit_one, jobs))
    b.write_json(output/'provenance.json', {'config': cfg, 'workers': workers, 'memory_budget_gb': 26,
                 'source_sha256': b.sha256(__file__), 'results': results})
    b.seal(output, 'curve.manifest.json', cfg, [j[2] for j in jobs]+[output/'provenance.json'])
    return {'cached': False, 'results': results}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('task', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--sizes', type=int, nargs='+', default=[4_000_000, 16_000_000, 64_000_000])
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    print(run(args.task, args.output, tuple(args.sizes), args.workers))
