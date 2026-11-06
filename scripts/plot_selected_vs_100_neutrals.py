#!/usr/bin/env python3
"""Illustrate one existing recent-selection scan call against 100 neutral regions.

Prepare a verified summary cache, then render without new simulations or RNG draws.
The example is deliberately chosen from existing mass calls, not an unbiased test.
"""
import os
os.environ['MPLBACKEND'] = 'Agg'
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import re
import shutil

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import PercentFormatter
import numpy as np
import pandas as pd

from plot_selection_episode_series import sha, save_json
from report_selection_arrays import local_path
from report_time_resolved_pilot import empirical_p

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / 'docs/results/selection_arrays_eas_20261105_4h'
SERIES = ROOT / 'docs/results/selection_episode_10k_s001_20261106'
METHODS = ('median', 'mass')
LABELS = ('Posterior median', 'Posterior mass')


def prepare(out, workers):
    inputs = {str(p.relative_to(ROOT)): sha(p) for p in (
        REPORT / 'completed_units.tsv', REPORT / 'calibration.npz',
        SERIES / 'series_data.npz', SERIES / 'input_manifest.json', SERIES / 'verified_calls.tsv')}
    for parent in (REPORT, SERIES):
        inventory = json.loads((parent / 'output_hashes.json').read_text())
        for name, digest in inputs.items():
            path = ROOT / name
            if path.parent == parent:
                assert digest == inventory[path.name], f'Changed source: {path}'
    contract = dict(version=1, inputs=inputs,
                    selection_rule='First replicate in frozen order with a mass joint_hit; maximum mass z within 400-600 kb and 1-10 kya; first flattened index breaks ties.',
                    null_rule='All 50 neutral_test plus first 50 calibration records by replicate; shared 10k introgression pool.')
    cache, manifest = out / 'curves.npz', out / 'manifest.json'
    if cache.exists() and manifest.exists():
        saved = json.loads(manifest.read_text())
        assert saved['contract'] == contract
        assert sha(cache) == saved['cache_sha256']
        print('Reusing verified selected/100-neutral curve cache.', flush=True)
        with np.load(cache) as data:
            return {k: data[k].copy() for k in data.files}, saved['example']
    records = pd.read_csv(REPORT / 'completed_units.tsv', sep='\t', low_memory=False)
    metrics = pd.read_csv(SERIES / 'verified_calls.tsv', sep='\t')
    hits = set(metrics[(metrics.method == 'mass') & (metrics.cohort == 'selected') & (metrics.joint_hit == 1)].id)
    with np.load(SERIES / 'series_data.npz') as data:
        selected_index = next(i for i, ident in enumerate(data['ids'][:100]) if str(ident) in hits)
        ident = str(data['ids'][selected_index])
        all_cdf = data['cdfs'][1:, selected_index].copy()
        mu, sd, maxima = data['mu'][1:].copy(), data['sd'][1:].copy(), data['maxima'][1:].copy()
    mass = np.diff(all_cdf, axis=-1)
    z = np.full(mass.shape, -np.inf)
    np.divide(mass - mu[:, None, :], sd[:, None, :], out=z, where=sd[:, None, :] > 0)
    eligible = np.full(z[1].shape, False)
    eligible[40:61, :9] = True
    pos_i, bin_i = np.unravel_index(np.argmax(np.where(eligible, z[1], -np.inf)), z[1].shape)
    assert empirical_p(z[1, pos_i, bin_i], maxima[1]) <= .05
    pool = records[(records.regime == '10k') & (records.origin == 'introgression')]
    null = pd.concat([pool[pool.cohort == 'neutral_test'].sort_values('replicate'),
                      pool[pool.cohort == 'calibration'].sort_values('replicate').head(50)])
    assert len(null) == 100 and null.id.nunique() == 100
    assert null.observed_sha256.nunique() == 100 and ident not in set(null.id)
    selected = records[records.id == ident]
    assert len(selected) == 1
    chosen = pd.concat([selected, null], ignore_index=True)

    def read(row):
        path = local_path(row.source_directory) / 'scores.npz'
        assert sha(path) == row.scores_sha256, f'Corrupt scores: {path}'
        with np.load(path) as data:
            np.testing.assert_array_equal(data['positions'], np.arange(100) * 10000)
            np.testing.assert_array_equal(data['cutoffs'], np.arange(1, 51) * 1000)
            np.testing.assert_array_equal(data['valid_pairs'], np.full(100, 1000))
            return data['cdfs'][1:, pos_i].copy()

    with ThreadPoolExecutor(max_workers=workers) as pool_executor:
        curves = np.stack(list(pool_executor.map(read, chosen.itertuples())), axis=1)
    np.testing.assert_array_equal(curves[:, 0], all_cdf[:, pos_i])
    assert np.isfinite(curves).all() and curves.min() >= 0 and curves.max() <= 1 + 1e-12
    assert np.diff(curves, axis=-1).min() >= -1e-12
    result = dict(cdfs=curves, mu=mu, sd=sd, maxima=maxima,
                  ids=chosen.id.to_numpy(str), bin_index=np.array(bin_i), position_bp=np.array(pos_i * 10000))
    example = dict(id=ident, selection=.01, selection_start_kya=10, selection_end_kya=0,
                   position_bp=int(pos_i * 10000), selected_coordinate_bp=500000,
                   bin_kya=[int(bin_i + 1), int(bin_i + 2)],
                   present_af=float(selected.iloc[0].final_frequency),
                   selection_onset_af=float(selected.iloc[0].selection_onset_frequency))
    chosen.to_csv(out / 'cohort_manifest.tsv', sep='\t', index=False, lineterminator='\n', na_rep='NA')
    np.savez_compressed(cache, **result)
    save_json(manifest, dict(contract=contract, example=example, cache_sha256=sha(cache)))
    print('Verified original score hashes and selected curve; original mass scan call confirmed.', flush=True)
    return result, example


def render(data, example, out):
    matplotlib.rcParams.update({'font.size': 15, 'axes.titlesize': 16, 'axes.labelsize': 16,
                               'xtick.labelsize': 13, 'ytick.labelsize': 13,
                               'pdf.fonttype': 42, 'ps.fonttype': 42})
    mass = np.diff(data['cdfs'], axis=-1)
    selected, neutral = mass[:, 0], mass[:, 1:]
    counts = (neutral >= selected[:, None]).sum(axis=1)
    point_p = (1 + counts) / 101
    threshold = np.sort(neutral, axis=1)[:, -5]
    np.testing.assert_array_equal(selected > threshold, point_p <= .05)
    z = np.full(selected.shape, -np.inf)
    np.divide(selected - data['mu'], data['sd'], out=z, where=data['sd'] > 0)
    scan_p = np.stack([empirical_p(z[j], data['maxima'][j]) for j in range(2)])
    k = int(data['bin_index'])
    assert scan_p[1, k] <= .05
    highlight_bins = [int(np.argmax(z[0, :9])), k]
    times = np.arange(1.5, 50)
    rows, summary = [], dict(example=example, methods={})
    for j, method in enumerate(METHODS):
        highlight = highlight_bins[j]
        summary['methods'][method] = dict(highlighted_bin_kya=[highlight + 1, highlight + 2],
            selected_bin_mass=float(selected[j, highlight]),
            neutrals_at_or_above=int(counts[j, highlight]), raw_tail_fraction=float(counts[j, highlight] / 100),
            pointwise_p=float(point_p[j, highlight]), original_scan_p=float(scan_p[j, highlight]),
            pointwise_significant_bins=(np.flatnonzero(point_p[j] <= .05) + 1).tolist(),
            scan_significant_bins=(np.flatnonzero(scan_p[j] <= .05) + 1).tolist())
        for t in range(49):
            rows.append(dict(method=method, lower_kya=t + 1, upper_kya=t + 2,
                             selected_mass=selected[j, t], neutral_threshold=threshold[j, t],
                             neutrals_at_or_above=int(counts[j, t]), pointwise_p=point_p[j, t], scan_p=scan_p[j, t]))
    pd.DataFrame(rows).to_csv(out / 'timewise_statistics.tsv', sep='\t', index=False, lineterminator='\n')
    save_json(out / 'summary.json', summary)
    with PdfPages(out / 'selected_vs_100_neutrals.pdf') as book:
        for kind in ('curves', 'pvalues'):
            fig, axes = plt.subplots(2, 1, figsize=(11, 8.5), sharex=True)
            fig.suptitle(f"s=0.01 · selection 10 kya–present · position {example['position_bp'] / 1000:g} kb", y=.99, fontsize=18)
            fig.text(.5, .945, f"{example['id'].split('/')[-1]} · present AF {example['present_af']:.1%} · chosen from mass scan calls", ha='center', fontsize=13)
            for j, ax in enumerate(axes):
                highlight = highlight_bins[j]
                info = summary['methods'][METHODS[j]]
                ax.axvspan(1, 10, color='#efefef', zorder=0)
                if kind == 'curves':
                    ax.plot(times, neutral[j].T, color='#7396b3', alpha=.23, lw=.7)
                    ax.plot(times, threshold[j], color='black', ls='--', lw=1.7)
                    ax.plot(times, selected[j], color='#D55E00', lw=2.7)
                    sig = point_p[j] <= .05
                    ax.scatter(times[sig], selected[j, sig], color='#D55E00', s=25, zorder=5)
                    ax.scatter(times[highlight], selected[j, highlight], facecolor='white', edgecolor='#D55E00', linewidth=2, s=80, zorder=6)
                    ax.yaxis.set_major_formatter(PercentFormatter(1))
                    ax.set_ylabel('Pairs / mass per 1-kyr bin')
                    ax.set_ylim(bottom=0)
                else:
                    ax.plot(times, point_p[j], color='#D55E00', lw=2.4)
                    ax.plot(times, scan_p[j], color='#0072B2', lw=2.4)
                    ax.axhline(.05, color='black', ls='--', lw=1.5)
                    ax.set(yscale='log', ylim=(.008, 1.2), ylabel='Empirical p-value')
                    ax.set_yticks([.01, .05, .1, 1], ['0.01', '0.05', '0.1', '1'])
                ax.set_title(f"{LABELS[j]} · {highlight+1}–{highlight+2} kya: {info['neutrals_at_or_above']}/100 ≥ selected; "
                             f"p={info['pointwise_p']:.3f}; scan p={info['original_scan_p']:.3f}", loc='left', fontsize=13)
                ax.set_xlim(1, 50)
                ax.set_xticks([1, 10, 20, 30, 40, 50])
                ax.grid(axis='y', alpha=.18)
            axes[-1].set_xlabel('Coalescence time (kya)')
            if kind == 'curves':
                handles = [Line2D([], [], color='#7396b3', label='100 neutral regions'),
                           Line2D([], [], color='#D55E00', lw=2.7, label='One selected region'),
                           Line2D([], [], color='black', ls='--', label='Pointwise 5% threshold')]
                footer = 'Dots: pointwise p ≤ 0.05; open circle: highlighted bin. Dashed threshold requires strict exceedance.'
            else:
                handles = [Line2D([], [], color='#D55E00', lw=2.4, label='Pointwise: 100 neutrals'),
                           Line2D([], [], color='#0072B2', lw=2.4, label='Scan: 99 neutral maxima'),
                           Line2D([], [], color='black', ls='--', label='p = 0.05')]
                footer = 'Scan p-values correct over 100 positions × 49 intervals; pointwise p-values do not.'
            fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.51, .91), ncol=3, frameon=False, fontsize=12)
            fig.text(.5, .055, footer, ha='center', fontsize=11)
            fig.text(.5, .025, 'Shading: evaluated part of selection episode. Neutral display: 50 held-out + 50 calibration; illustrative example.', ha='center', fontsize=10.5)
            fig.subplots_adjust(top=.79, bottom=.16, left=.12, right=.98, hspace=.42)
            name = '01_neutral_curves' if kind == 'curves' else '02_pointwise_and_scan_pvalues'
            fig.savefig(out / f'{name}.png', dpi=180)
            fig.savefig(out / f'{name}.pdf')
            book.savefig(fig)
            plt.close(fig)
    checks = {}
    for path in out.glob('*.pdf'):
        content = path.read_bytes()
        checks[path.name] = dict(embedded_truetype=b'/FontFile2' in content,
                                raster_images=len(re.findall(rb'/Subtype\s*/Image', content)))
        assert checks[path.name]['embedded_truetype'] and checks[path.name]['raster_images'] == 0
    save_json(out / 'pdf_validation.json', checks)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=ROOT / 'docs/results/selected_vs_100_neutrals_10k_s001_20261106')
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    assert 1 <= args.workers <= 24
    args.out.mkdir(parents=True, exist_ok=True)
    data, example = prepare(args.out, args.workers)
    summary = render(data, example, args.out)
    save_json(args.out / 'provenance.json', dict(script_sha256=sha(__file__), workers=args.workers,
        versions={n: __import__(n).__version__ for n in ('numpy', 'pandas', 'matplotlib')},
        parameters=dict(regime='10k', origin='introgression', selection=.01, selected_survival_condition='0<p<1 today',
                        neutral_condition='unconditioned; original fixed coordinates', pairs=1000,
                        stride_bp=10000, bin_kya=1, temporal_scan_kya=[1, 50], alpha=.05),
        randomness='No simulations or RNG draws. Fixed deterministic cohort/example rules in manifest.',
        interpretation='Example selected after inspecting scan calls; pointwise p-values are descriptive, not adjusted for choosing this example. Half of the displayed neutrals contributed to original calibration.'))
    text = f'''# One selected simulation versus 100 neutrals

Open [the two-page PDF](selected_vs_100_neutrals.pdf): posterior-median and posterior-mass bin curves, then pointwise versus scan-wide p-values.

Example: `{example['id']}`, position {example['position_bp']} bp (selected allele at 500000 bp), highlighted mass bin {example['bin_kya']} kya. Present AF {example['present_af']:.2%}. Selection s=0.01 from 10 kya to present; EAS demography, 1% Neanderthal pulse at 50 kya. Both methods use the same selected replicate, genomic position, and 1000 pairs. Each panel highlights its own maximum standardized score within 1–10 kya at that position; highlighted bins may differ.

The example is the first previously detected posterior-mass replicate in frozen order. The position/bin maximizes its standardized mass score within ±100 kb and 1–10 kya. This is an illustration chosen to show a call, not an unbiased performance estimate. The original scan correction covers positions and times within each region, not selection among 100 replicates for illustration.

The 100 displayed neutrals are all 50 held-out regions plus the first 50 calibration regions by replicate. They are distinct simulations, with input hashes verified; this combined set is not a new independent validation cohort. Both selected and neutral curves use the same genomic coordinate, not each neutral's nearest-segregating label.

For each time bin k is the number of neutral values **equal to or above** the selected value. Raw tail fraction k/100; conservative empirical p=(k+1)/101. Thus at most 4/100 gives p≤0.05; 5/100 gives p=0.0594. The dashed curve is the 96th ascending order statistic, which the selected curve must strictly exceed. It is a pointwise rank threshold, not a confidence band. Posterior-median curves count pairs with medians in the bin; posterior-mass curves sum pairwise posterior probabilities in the bin. Neither plot is a pairwise TMRCA credible interval.

The second page preserves the original scan-wide p-values from 99 independent calibration-region maxima, standardized using 20 development regions and searching 100 genomic positions × 49 time bins. No threshold was refit for the selected example. No 0–1 kya bin is evaluated. Detailed values are in `timewise_statistics.tsv` and `summary.json`.

## Reproduce from the verified summary cache

```bash
%%bash
set -euo pipefail
REPO=gamma_smc_ts_selection
if [[ -d "$REPO/.git" ]]; then
  git -C "$REPO" pull --ff-only
else
  git clone --branch AOU_run_opt git@github.com:soisa001/gamma_smc_ts.git "$REPO"
fi
cd "$REPO"
uv run --frozen python scripts/plot_selected_vs_100_neutrals.py --workers 4
```
'''
    (args.out / 'README.md').write_text(text, encoding='utf-8', newline='\n')
    shutil.copy2(__file__, args.out / Path(__file__).name)
    (args.out / '.gitattributes').write_text('* -text\n', encoding='utf-8')
    save_json(args.out / 'output_hashes.json', {p.name: sha(p) for p in sorted(args.out.iterdir()) if p.is_file() and p.name != 'output_hashes.json'})
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()
