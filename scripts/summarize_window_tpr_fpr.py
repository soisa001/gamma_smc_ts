#!/usr/bin/env python3
"""Score selected TPR and independent-neutral FPR within the same time window."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import numpy as np
from summarize_selection_power_fpr import METHODS, RESULTS, ROOT, WINDOWS

BOUNDS = {'10k': (0, 10), '50k': (0, 50), '50_30k': (30, 50)}
CUTOFFS = np.arange(1, 51)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def window_mask(regime, statistic):
    low, high = BOUNDS[regime]
    if statistic == 'cumulative':
        # A called cutoff is a point on the time axis, not a call of all [0,t].
        return (CUTOFFS >= low) & (CUTOFFS <= high)
    assert statistic == 'interval'
    # Positive-duration overlap; touching the boundary alone does not count.
    return (CUTOFFS[:-1] < high) & (CUTOFFS[1:] > low)


def called(p, regime, statistic, alpha):
    return (p[:, window_mask(regime, statistic)] < alpha).any(axis=1)


def verify_ranks(data, null):
    """Independently recalculate each simulation's empirical rank, preserving ties."""
    assert null.shape == (100, 3, 50)
    for m in range(3):
        x, y = data['cdfs'][:, m], null[:, m]
        if m < 2:
            x, y = np.rint(x * 1000).astype(int), np.rint(y * 1000).astype(int)
        for statistic in ('cumulative', 'interval'):
            a, b = (x, y) if statistic == 'cumulative' else (np.diff(x), np.diff(y))
            expected = np.array([(1 + (b >= row).sum(axis=0)) / 101 for row in a]).reshape(len(a), a.shape[-1])
            np.testing.assert_array_equal(data['p_' + statistic][:, m], expected)


def boundary_checks():
    # An outside-window call does not count; endpoint-only bin contact is excluded.
    p = np.ones((3, 50)); p[0, 28] = .009; p[1, 29] = .009; p[2, 49] = .009
    np.testing.assert_array_equal(called(p, '50_30k', 'cumulative', .01), [False, True, True])
    q = np.ones((3, 49)); q[0, 28] = .009; q[1, 29] = .009; q[2, 48] = .009
    np.testing.assert_array_equal(called(q, '50_30k', 'interval', .01), [False, True, True])
    q = np.ones((3, 49)); q[0, 8] = .009; q[1, 9] = .009; q[2, 8] = .05
    np.testing.assert_array_equal(called(q, '10k', 'interval', .05), [True, False, False])
    assert not called(np.empty((0, 50)), '50k', 'cumulative', .05).size


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=RESULTS / 'window_tpr_fpr_20261116')
    parser.add_argument('--extension', default='selected_survivor_extension_20261112', help='Continuation report directory within docs/results.')
    parser.add_argument('--previous-summary', default='power_fpr_summary_20261116', help='Matching any-time summary directory within docs/results.')
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    assert 1 <= args.workers <= 20
    boundary_checks()
    provenance = {}

    def load(folder, name, expected=None):
        directory = RESULTS / folder; path = directory / name
        inventory = json.loads((directory / 'output_hashes.json').read_text())
        actual = digest(path)
        assert actual == inventory[name] and (expected is None or actual == expected), str(path)
        provenance[str(path.relative_to(ROOT))] = actual
        if path.suffix == '.json': return json.loads(path.read_text())
        with np.load(path, allow_pickle=False) as f: return {k: f[k].copy() for k in f.files}

    base, extension, neutral_folder = 'denovo_comparison_20261111', args.extension, 'segregating_neutral_fpr_100_20261110'
    meta = load(base, 'manifest.json'); ext = load(extension, 'manifest.json')
    old = load(args.previous_summary, 'summary.json')
    previous = {(r['origin'], r['regime'], r['selection'], r['method'], r['statistic'], r['alpha']): r for r in old['power']}
    cases = {}
    for c in meta['cases']:
        for origin in ('de_novo', 'introgression'):
            cases[(origin, c['regime'], c['selection'])] = dict(origin=origin, regime=c['regime'], selection=c['selection'],
                source=base, selected=meta['cohorts'][c['selected'][origin]], neutral=meta['cohorts'][c['neutral'][origin]])
    for c in ext['cases']:
        cases[(c['origin'], c['regime'], c['selection'])] = dict(origin=c['origin'], regime=c['regime'], selection=c['selection'],
            source=extension, selected=c, neutral=ext['nulls'][c['origin']])
    keys = sorted(cases, key=lambda k: (list(WINDOWS).index(k[1]), k[2], k[0]))
    assert len(keys) == 21

    def evaluate(key):
        c = cases[key]
        data = load(c['source'], c['selected']['cache'], c['selected']['cache_sha256'])
        ref = load(c['source'], c['neutral']['cache'], c['neutral']['cache_sha256'])
        assert len(set(data['ids'])) == len(data['ids'])
        assert set(data['ids']).isdisjoint(ref['ids'])
        for x in (data, ref): assert np.all((x['af'][:, -1] > 0) & (x['af'][:, -1] < 1))
        verify_ranks(data, ref['cdfs'])
        rows = []
        for m, method in enumerate(METHODS):
            for statistic in ('cumulative', 'interval'):
                for alpha in (.05, .01):
                    p = data['p_' + statistic][:, m]
                    flags = called(p, c['regime'], statistic, alpha); n = len(flags)
                    old_row = previous[key + (method, statistic, alpha)]
                    assert n == old_row['selected_n']
                    assert (int((p < alpha).any(1).sum()) if n else None) == old_row['any_time_count']
                    count = int(flags.sum()) if n else None
                    assert n == 0 or count <= old_row['any_time_count']
                    if c['regime'] == '50k': assert count == old_row['any_time_count']
                    rows.append(dict(origin=c['origin'], regime=c['regime'], selection=c['selection'], method=method,
                        statistic=statistic, alpha=alpha, selected_n=n, tp=count, tpr=count/n if n else None,
                        previous_any_time_count=old_row['any_time_count'], called_ids=data['ids'][flags].tolist(), source=c['source']))
        return rows

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        power = [r for rows in pool.map(evaluate, keys) for r in rows]
    neutral = load(neutral_folder, 'parameters.npz')
    np.testing.assert_array_equal(neutral['cutoffs_kya'], CUTOFFS)
    assert len(neutral['test_ids']) == len(neutral['reference_ids']) == 100
    assert set(neutral['test_ids']).isdisjoint(neutral['reference_ids'])
    intro = load(base, meta['cohorts']['introgression_neutral']['cache'])
    np.testing.assert_array_equal(intro['ids'], neutral['reference_ids'])
    np.testing.assert_array_equal(intro['cdfs'], neutral['reference_cdfs'])
    verify_ranks(dict(cdfs=neutral['test_cdfs'], p_cumulative=neutral['p_cumulative'], p_interval=neutral['p_interval']), neutral['reference_cdfs'])
    prior_fpr = {(r['method'], r['statistic'], r['threshold']): r for r in old['fpr']['results']}
    fpr = []
    for regime in BOUNDS:
        for m, method in enumerate(METHODS):
            for statistic in ('cumulative', 'interval'):
                for alpha in (.05, .01):
                    flags = called(neutral['p_' + statistic][:, m], regime, statistic, alpha)
                    count = int(flags.sum()); old_count = prior_fpr[(method, statistic, alpha)]['any_time_false_positive_count']
                    assert count <= old_count
                    if regime == '50k': assert count == old_count
                    fpr.append(dict(origin='introgression', regime=regime, method=method, statistic=statistic, alpha=alpha,
                        neutral_n=100, fp=count, fpr=count/100, previous_any_time_count=old_count,
                        called_ids=neutral['test_ids'][flags].tolist()))
    notes = [
        'TPR: fraction of selected simulations with >=1 call overlapping the true selection window.',
        'FPR: fraction of independent neutral simulations with >=1 call in that SAME window against the fixed null bank.',
        'Both endpoints use only the focal site (500 kb), strict pointwise p<alpha, and existing survivor-conditioned EAS simulations.',
        'Cumulative tests count called cutoff coordinates inside the closed window; a significant cutoff t is not treated as calling all ages below t.',
        'Interval tests require positive-duration overlap of the 1-kya bin with the window; boundary contact alone is excluded.',
        'Existing grid: cumulative cutoffs 1,...,50 kya; interval bins 1-2,...,49-50 kya. There is no 0-1 kya interval bin.',
        'FPR is available only for introgression (100 held-out neutrals vs 100 reference neutrals); de novo FPR remains unavailable.',
        'The same neutral cohort is evaluated in each window, so FPR does not depend on s. No new simulations or threshold changes.',
        'Calls outside the window do not contribute to either endpoint. Empty selected cohorts are NA, not zero TPR.',
        'Cumulative cutoff overlap is an operational detection definition; cumulative significance does not identify a coalescence interval.'
    ]
    text = ['WINDOW-MATCHED TPR AND FPR', '', *notes, '']
    index = {(r['origin'], r['regime'], r['selection'], r['method'], r['statistic'], r['alpha']): r for r in power}
    fi = {(r['regime'], r['method'], r['statistic'], r['alpha']): r for r in fpr}
    for method in METHODS:
        text += [method.upper(), 'TPR counts: model | window | s | cumulative p<.05 | interval p<.05 | cumulative p<.01 | interval p<.01']
        for key in keys:
            entries = [index[key + (method, statistic, alpha)] for alpha in (.05, .01) for statistic in ('cumulative', 'interval')]
            text.append(' | '.join([key[0], WINDOWS[key[1]], f'{key[2]:g}', *[f"{r['tp']}/{r['selected_n']}" if r['selected_n'] else 'NA (n=0)' for r in entries]]))
        text += ['', 'Introgression FPR counts: window | cumulative p<.05 | interval p<.05 | cumulative p<.01 | interval p<.01']
        for regime in BOUNDS:
            text.append(' | '.join([WINDOWS[regime], *[f"{fi[(regime, method, statistic, alpha)]['fp']}/100" for alpha in (.05, .01) for statistic in ('cumulative', 'interval')]]))
        text.append('')
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / '.gitattributes').write_text('* -text\n')
    (args.out / 'window_tpr_fpr_tables.txt').write_text('\n'.join(text).rstrip() + '\n')
    validation = dict(cohorts=21, tpr_rows=len(power), fpr_rows=len(fpr), workers=args.workers,
        checks='Source/cache hashes; independent tie-aware ranks; reference bank equality; unique survivor IDs; window boundary examples; subset of prior any-time counts; full-window equivalence.')
    assert len(power) == 252 and len(fpr) == 36
    payload = dict(notes=notes, bounds_kya=BOUNDS, tpr=power, fpr=fpr, source_sha256=dict(sorted(provenance.items())),
                   script_sha256=digest(Path(__file__)), validation=validation)
    (args.out / 'summary.json').write_text(json.dumps(payload, indent=2) + '\n')
    (args.out / 'output_hashes.json').write_text(json.dumps({p.name: digest(p) for p in sorted(args.out.iterdir()) if p.is_file() and p.name != 'output_hashes.json'}, indent=2) + '\n')
    print(json.dumps(validation))


if __name__ == '__main__':
    main()
