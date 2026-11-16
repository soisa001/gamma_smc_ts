#!/usr/bin/env python3
"""Consolidate validated survivor-cohort detection and independent neutral FPR.

No simulations or significance thresholds are changed. The November 12
continuation replaces earlier counts for its seven cohorts, rather than being
added as an independent cohort. Outputs are deterministic text and JSON.
"""
import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / 'docs/results'
METHODS = ['True TMRCA', 'Posterior median', 'Posterior-weighted mass']
WINDOWS = {'10k': '10 kya-present', '50k': '50 kya-present', '50_30k': '50-30 kya'}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=RESULTS / 'power_fpr_summary_20261116')
    args = parser.parse_args()
    provenance = {}

    def read(folder, filename):
        directory = RESULTS / folder
        inventory = json.loads((directory / 'output_hashes.json').read_text())
        path = directory / filename
        actual = digest(path)
        assert actual == inventory[filename], str(path)
        provenance[str(path.relative_to(ROOT))] = actual
        return json.loads(path.read_text())

    base = 'denovo_comparison_20261111'
    extension = 'selected_survivor_extension_20261112'
    cohort_rows = {}
    for row in read(base, 'summary.json')['results']:
        key = (row['origin'], row['regime'], row['selection'], row['method'], row['statistic'], row['alpha'])
        assert key not in cohort_rows
        cohort_rows[key] = dict(row, source=base)
    ext_meta = {c['key']: c for c in read(extension, 'manifest.json')['cases']}
    for row in read(extension, 'summary.json')['results']:
        case = ext_meta[row['case']]
        assert row['selected_n'] == case['n']
        key = (case['origin'], case['regime'], case['selection'], row['method'], row['statistic'], row['alpha'])
        cohort_rows[key] = dict(row, origin=case['origin'], regime=case['regime'],
                               selection=case['selection'], source=extension)
    rows = list(cohort_rows.values())
    assert len(rows) == 252
    cohorts = sorted({k[:3] for k in cohort_rows}, key=lambda k: (list(WINDOWS).index(k[1]), k[2], k[0]))
    assert len(cohorts) == 21
    ns = []
    for origin, regime, selection in cohorts:
        group = [r for r in rows if (r['origin'], r['regime'], r['selection']) == (origin, regime, selection)]
        assert len(group) == 12 and len({r['selected_n'] for r in group}) == 1
        ns.append(group[0]['selected_n'])
        for r in group:
            n, count = r['selected_n'], r['any_time_count']
            assert r['neutral_n'] == 100
            assert (count is None and r['any_time_fraction'] is None) if n == 0 else (
                0 <= count <= n and abs(r['any_time_fraction'] - count / n) < 1e-12)
    assert sorted(ns) == [0] * 5 + [2, 5] + [100] * 14
    fpr = read('segregating_neutral_fpr_100_20261110', 'summary.json')
    assert fpr['n_test'] == fpr['n_reference'] == 100
    focal = {(r['method'], r['statistic'], r['threshold']): r for r in fpr['results']}
    regional = {(r['method'], r['statistic'], r['threshold']): r for r in fpr['regional_any_position_and_time']}
    assert len(focal) == len(regional) == 12
    for key, r in focal.items():
        assert r['any_time_fpr'] == r['any_time_false_positive_count'] / 100
        assert regional[key]['fpr'] == regional[key]['false_positive_count'] / 100
        assert regional[key]['false_positive_count'] >= r['any_time_false_positive_count']

    notes = [
        'Current EAS three-window series, including the November 12 continuation; older pilot grids are excluded.',
        'Selected and reference cohorts require 0 < original focal allele population AF < 1 at present.',
        'Detection: >=1 strict pointwise p<alpha at the 500-kb focal site, over 1-50 kya cutoffs or 1-2,...,49-50 kya intervals.',
        'These are conditional any-time detection rates, not scan-wide calibrated power or selection-window localization accuracy.',
        'Rank p=(1 + number of reference values >= selected/test value)/101. 1,000 fixed pairwise comparisons per simulation.',
        'FPR: 100 independent introgression neutrals against the frozen 100-reference bank; no independent de novo FPR cohort.',
        'Regional FPR searches 100 positions (0-990 kb, 10-kb stride) and all times in each 1-Mb neutral simulation.',
        'No new independent FPR experiment, multiple-testing adjustment, interval-overlap TPR, or F1 calculation is supplied here.',
        'NA means no retained survivors; small n=2 and n=5 cohorts cannot support stable power estimates.',
        '0/100 observed false positives has a 95% Wilson upper bound of about 3.7%, conditional on this reference bank.'
    ]
    lines = ['SELECTION DETECTION AND FPR SUMMARY', '', *notes, '']
    for method in METHODS:
        lines += [method.upper(), 'Each entry is detections / retained simulations.',
                  'Model | Selection window | s | Cumulative p<.05 | Interval p<.05 | Cumulative p<.01 | Interval p<.01']
        for origin, regime, selection in cohorts:
            values = []
            for alpha in (.05, .01):
                for statistic in ('cumulative', 'interval'):
                    r = cohort_rows[(origin, regime, selection, method, statistic, alpha)]
                    values.append(f"{r['any_time_count']}/{r['selected_n']}" if r['selected_n'] else 'NA (n=0)')
            lines.append(' | '.join([origin, WINDOWS[regime], f'{selection:g}', *values]))
        lines.append('')
    lines += ['INDEPENDENT INTROGRESSION NULL FPR', 'Each count is out of 100 test neutrals (also the percentage).',
              'Method | Statistic | Focal any-time p<.05 | Focal any-time p<.01 | Region any-position/time p<.05 | Region any-position/time p<.01']
    for method in METHODS:
        for statistic in ('cumulative', 'interval'):
            values = [focal[(method, statistic, a)]['any_time_false_positive_count'] for a in (.05, .01)]
            values += [regional[(method, statistic, a)]['false_positive_count'] for a in (.05, .01)]
            lines.append(' | '.join([method, statistic, *[f'{v}/100' for v in values]]))
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / '.gitattributes').write_text('* -text\n')
    (args.out / 'power_fpr_tables.txt').write_text('\n'.join(lines) + '\n')
    compact_rows = [{k: v for k, v in r.items() if k != 'fraction_by_time'} for r in rows]
    payload = dict(notes=notes, power=compact_rows, fpr=fpr, source_sha256=provenance,
                   script_sha256=digest(Path(__file__)), validation=dict(cohorts=21, power_rows=252,
                   full_cohorts=14, small_cohorts=2, empty_cohorts=5, summary_hashes_and_counts_verified=True))
    (args.out / 'summary.json').write_text(json.dumps(payload, indent=2) + '\n')
    (args.out / 'output_hashes.json').write_text(json.dumps({p.name: digest(p) for p in sorted(args.out.iterdir())
        if p.is_file() and p.name != 'output_hashes.json'}, indent=2) + '\n')
    print(json.dumps(payload['validation']))


if __name__ == '__main__':
    main()
