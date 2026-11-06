#!/usr/bin/env python3
"""Describe mass/median agreement at the existing 420-kb example coordinate.

No significance tests: neighboring bins and shared haplotype pairs are dependent.
Reads checked caches only; no simulations, RNG draws, or dependency changes.
"""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import platform
import numpy as np
import pandas as pd
import scipy
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def corr(a, b):
    a, b = np.asarray(a).ravel(), np.asarray(b).ravel()
    assert a.shape == b.shape and np.isfinite(a).all() and np.isfinite(b).all()
    if np.ptp(a) == 0 or np.ptp(b) == 0:
        return dict(n=len(a), pearson=None, spearman=None)
    return dict(n=len(a), pearson=float(np.corrcoef(a, b)[0, 1]),
                spearman=float(spearmanr(a, b).statistic))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--out', type=Path, default=ROOT / 'docs/results/mass_median_correlation_10k_s001_20261106')
    args = parser.parse_args()
    assert 1 <= args.workers <= 24
    example_dir = ROOT / 'docs/results/selected_vs_100_neutrals_10k_s001_20261106'
    series_dir = ROOT / 'docs/results/selection_episode_10k_s001_20261106'
    paths = [example_dir / 'curves.npz', example_dir / 'manifest.json', series_dir / 'series_data.npz']
    hashes = {}
    for path in paths:
        inventory = json.loads((path.parent / 'output_hashes.json').read_text())
        hashes[str(path.relative_to(ROOT))] = sha(path)
        assert hashes[str(path.relative_to(ROOT))] == inventory[path.name]
    example = json.loads((example_dir / 'manifest.json').read_text())['example']
    with np.load(paths[0]) as data:
        coordinate = int(data['position_bp'])
        example_cdf = data['cdfs'][:, :1].copy()
        neutral_cdf = data['cdfs'][:, 1:].copy()
        neutral_ids = data['ids'][1:].copy()
    with np.load(paths[2]) as data:
        selected_cdf = data['cdfs'][1:, :100, coordinate // 10000].copy()
        selected_ids = data['ids'][:100].copy()
    assert coordinate == example['position_bp'] == 420000
    np.testing.assert_array_equal(selected_cdf[:, np.flatnonzero(selected_ids == example['id'])[0]], example_cdf[:, 0])
    assert len(set(neutral_ids)) == len(set(selected_ids)) == 100
    rows, per_replicate, across_replicates = [], [], []
    for label, cdf, ids in [('example', example_cdf, [example['id']]),
                             ('selected', selected_cdf, selected_ids), ('neutral', neutral_cdf, neutral_ids)]:
        for statistic, values, lower in [('one_kyr_fraction', np.diff(cdf, axis=-1), np.arange(1, 50)),
                                         ('cumulative_fraction', cdf, np.arange(1, 51))]:
            for interval, keep in [('1_50kya', np.ones(len(lower), dtype=bool)),
                                   ('1_10kya', lower < 10 if statistic == 'one_kyr_fraction' else lower <= 10)]:
                subset = values[..., keep]
                rows.append(dict(cohort=label, statistic=statistic, interval=interval,
                                 aggregation='pooled_replicates_and_time', **corr(subset[0], subset[1])))
                rows.append(dict(cohort=label, statistic=statistic, interval=interval,
                                 aggregation='cohort_mean_time_profile', **corr(subset[0].mean(axis=0), subset[1].mean(axis=0))))
                with ThreadPoolExecutor(max_workers=args.workers) as pool:
                    correlations = list(pool.map(lambda i: corr(subset[0, i], subset[1, i]), range(len(ids))))
                per_replicate.extend(dict(cohort=label, id=str(ident), statistic=statistic, interval=interval, **c)
                                     for ident, c in zip(ids, correlations))
            if label != 'example' and statistic == 'one_kyr_fraction':
                for index, age in enumerate(lower):
                    across_replicates.append(dict(cohort=label, lower_kya=int(age), upper_kya=int(age + 1),
                                                 **corr(values[0, :, index], values[1, :, index])))
    args.out.mkdir(parents=True, exist_ok=True)
    for name, records in [('correlations', rows), ('per_replicate', per_replicate), ('across_replicates_by_time', across_replicates)]:
        pd.DataFrame(records).to_csv(args.out / f'{name}.tsv', sep='\t', index=False, lineterminator='\n', na_rep='NA')
    frame = pd.DataFrame(per_replicate)
    distribution = []
    for keys, group in frame.groupby(['cohort', 'statistic', 'interval'], sort=False):
        record = dict(zip(['cohort', 'statistic', 'interval'], keys))
        for method in ('pearson', 'spearman'):
            valid = group[method].dropna()
            record[method] = dict(valid=len(valid), undefined=len(group) - len(valid),
                                 q25=float(valid.quantile(.25)), median=float(valid.median()), q75=float(valid.quantile(.75)))
        distribution.append(record)
    report = dict(example=example, position_bp=coordinate, correlations=rows, per_replicate_distribution=distribution,
                  interpretation='Primary statistic is the fraction in each 1-kyr bin, not cumulative fractions. Pooled, cohort-mean-profile, within-replicate, and between-replicate correlations answer different questions. Cumulative curves share monotonicity and can correlate highly despite differing temporal distributions. No correlation p-values or confidence intervals computed.',
                  design='100 selected survivors, s=0.01, selection 10 kya-present, EAS introgression; 100 neutrals (50 held-out + 50 calibration). Coordinate and example were chosen from a mass call, so these are descriptive results at an ascertained coordinate, not genome-wide estimates.',
                  provenance=dict(inputs=hashes, script_sha256=sha(Path(__file__)), workers=args.workers,
                                  python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__,
                                  randomness='No RNG draws; all fixed original simulations and pairs retained.'))
    (args.out / 'summary.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    (args.out / '.gitattributes').write_text('* -text\n', encoding='utf-8')
    (args.out / 'output_hashes.json').write_text(json.dumps({p.name: sha(p) for p in sorted(args.out.iterdir()) if p.is_file() and p.name != 'output_hashes.json'}, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'correlations': rows, 'per_replicate_distribution': distribution}, indent=2), flush=True)


if __name__ == '__main__':
    main()
