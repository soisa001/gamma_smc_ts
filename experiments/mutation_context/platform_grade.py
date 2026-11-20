"""Private Harbor verifier: score submitted probabilities against frozen truth."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np


def score(predictions, queries, oracle_file):
    predictions = Path(predictions)
    if not predictions.is_file() or predictions.is_symlink():
        raise ValueError('Missing regular predictions.tsv file')
    if predictions.stat().st_size > 300_000_000:
        raise ValueError('Prediction file exceeds 300 MB')
    with np.load(oracle_file, allow_pickle=False) as oracle:
        truth = oracle['probabilities']
        parent, child = oracle['parent'], oracle['child']
    n = len(parent)
    if truth.shape != (n, 4) or child.shape != (n,):
        raise ValueError('Invalid evaluator oracle')
    expected = float((1-truth[np.arange(n), parent]).sum())
    totals = dict(sites=0, oracle_kl=0., predicted_mutations=0., observed_nll=0.,
                  oracle_conditional_kl_total=0.)
    batch = []

    def flush():
        if not batch:
            return
        p = np.array(batch, dtype=float)
        if (not np.isfinite(p).all() or (p <= 0).any() or (p > 1).any()
                or not np.allclose(p.sum(axis=1), 1, rtol=0, atol=1e-9)):
            raise ValueError('Invalid probabilities or infinite log loss from zero probabilities')
        start, end = totals['sites'], totals['sites']+len(p)
        q, anc, obs = truth[start:end], parent[start:end], child[start:end]
        rows = np.arange(len(p))
        totals['oracle_kl'] += float((q*(np.log(q)-np.log(p))).sum())
        totals['predicted_mutations'] += float((1-p[rows, anc]).sum())
        totals['observed_nll'] -= float(np.log(p[rows, obs]).sum())
        qa, pa = q.copy(), p.copy()
        qa[rows, anc] = pa[rows, anc] = 0
        mu, pred_mu = qa.sum(axis=1), pa.sum(axis=1)
        mask = qa > 0
        logq, logp = np.zeros_like(qa), np.zeros_like(pa)
        np.log(qa, out=logq, where=mask)
        np.log(pa, out=logp, where=mask)
        difference = logq-np.log(mu[:, None])-logp+np.log(pred_mu[:, None])
        totals['oracle_conditional_kl_total'] += float((qa*difference).sum())
        totals['sites'] = end
        batch.clear()

    with Path(queries).open(newline='') as qh, predictions.open(newline='') as ph:
        expected_rows, actual = csv.reader(qh, delimiter='\t'), csv.reader(ph, delimiter='\t')
        if next(expected_rows, None) != ['contig', 'position']:
            raise ValueError('Invalid evaluator query header')
        if next(actual, None) != ['contig', 'position', 'p_A', 'p_C', 'p_G', 'p_T']:
            raise ValueError('Incorrect prediction header')
        count = 0
        for query in expected_rows:
            row = next(actual, None)
            if row is None or len(row) != 6 or row[:2] != query:
                raise ValueError(f'Query identity/order mismatch at row {count+1}')
            batch.append([float(x) for x in row[2:]])
            count += 1
            if count > n:
                raise ValueError('Evaluator query count exceeds oracle')
            if len(batch) == 50_000:
                flush()
        if next(actual, None) is not None or count != n:
            raise ValueError('Incorrect prediction row count')
        flush()
    kl = max(0., totals['oracle_kl']/expected)
    totals.update(expected_mutations=expected, observed_mutations=int((parent != child).sum()),
                  oracle_kl_per_mutation=kl,
                  oracle_conditional_kl=totals['oracle_conditional_kl_total']/expected,
                  observed_nll_per_site=totals['observed_nll']/n,
                  burden_relative_error=totals['predicted_mutations']/expected-1,
                  reward=math.exp(-kl), valid=True)
    h = hashlib.sha256()
    with predictions.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024*1024), b''):
            h.update(chunk)
    totals['prediction_sha256'] = h.hexdigest()
    return totals


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--predictions', type=Path, default=Path('/app/output/predictions.tsv'))
    parser.add_argument('--private', type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument('--logs', type=Path, default=Path('/logs/verifier'))
    args = parser.parse_args()
    args.logs.mkdir(parents=True, exist_ok=True)
    reward_file = args.logs/'reward.txt'
    reward_file.write_text('0\n')
    try:
        result = score(args.predictions, args.private/'queries.tsv', args.private/'oracle.npz')
    except (ValueError, OSError, TypeError, KeyError, OverflowError, csv.Error) as error:
        result = dict(reward=0., valid=False, error=str(error)[:2000])
    (args.logs/'score.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    reward_file.write_text(str(result['reward'])+'\n')
    print(json.dumps(result, allow_nan=False))


if __name__ == '__main__':
    main()
