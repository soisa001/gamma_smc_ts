#!/usr/bin/env python3
"""Plot a complete 100-survivor introgression selection cohort.

Phase 1 verifies and caches only the required existing outputs. Phase 2 validates
calls against the frozen report and writes figures. No simulations or RNG draws.
"""
import os
os.environ['MPLBACKEND'] = 'Agg'
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import re
import shutil

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import PercentFormatter
import numpy as np
import pandas as pd

from report_selection_arrays import evaluate, local_path
from report_time_resolved_pilot import empirical_p

ROOT = Path(__file__).resolve().parents[1]
METHODS = ('truth_all_pairs', 'median', 'mass')
LABELS = ('True TMRCAs · 19,900 pairs', 'Posterior median · 1,000 pairs',
          'Posterior mass · 1,000 pairs')
ORANGE, BLUE = '#D55E00', '#0072B2'
EPISODES = {'40_20k': (40, 20, '40–20 kya'), '10k': (10, 0, '10 kya–present')}


def milestones(regime):
    if regime == '10k':
        return [(50, 'Introgression'), (10, 'Selection starts'), (5, 'During selection'), (0, 'Present')]
    return [(50, 'Introgression'), (40, 'Selection starts'), (20, 'Selection ends'), (0, 'Present')]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save_json(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n', encoding='utf-8', newline='\n')


def prepare(report, out, workers, selection=.02, regime='40_20k'):
    start, end, _ = EPISODES[regime]
    names = ('completed_units.tsv', 'calibration.npz', 'replicate_metrics.tsv',
             'all_pairs_truth_cdfs.npz', 'all_pairs_truth_provenance.json',
             'all_pairs_truth_replicates.tsv', 'proposal_summary.tsv')
    inventory = json.loads((report / 'output_hashes.json').read_text(encoding='utf-8'))
    inputs = {name: sha(report / name) for name in names}
    assert all(inputs[name] == inventory[name] for name in names), 'Changed source report'
    records = pd.read_csv(report / 'completed_units.tsv', sep='\t', low_memory=False)
    selected = records[(records.regime == regime) & (records.origin == 'introgression') &
                       (records.cohort == 'selected') & (records.selection == selection)].sort_values('replicate')
    neutral = records[(records.regime == '10k') & (records.origin == 'introgression') &
                      (records.cohort == 'neutral_test')].sort_values('replicate')
    assert len(selected) == 100 and len(neutral) == 50
    chosen = pd.concat([selected, neutral], ignore_index=True)
    contract = dict(version=1, input_sha256=inputs, ids=chosen.id.tolist(),
                    selected=100, neutral=50, focal_position_bp=500000,
                    selection_start_kya=start, selection_end_kya=end, selection=selection)
    cache, manifest = out / 'series_data.npz', out / 'input_manifest.json'
    if cache.exists() and manifest.exists():
        saved = json.loads(manifest.read_text(encoding='utf-8'))
        assert saved['contract'] == contract, 'Cached input contract differs'
        assert saved['data_sha256'] == sha(cache), 'Corrupt plot cache'
        print('Reusing verified plot data; no raw simulation files needed.', flush=True)
    else:
        print('Checking trajectories and scores for 100 selected + 50 held-out null regions.', flush=True)

        def read(row):
            directory = local_path(row.source_directory)
            assert sha(directory / 'scores.npz') == row.scores_sha256
            with np.load(directory / 'scores.npz') as data:
                np.testing.assert_array_equal(data['positions'], np.arange(100) * 10000)
                np.testing.assert_array_equal(data['cutoffs'], np.arange(1, 51) * 1000)
                np.testing.assert_array_equal(data['valid_pairs'], np.full(100, 1000))
                scores = data['cdfs'].copy()
            trajectory = None
            if row.cohort == 'selected':
                assert sha(directory / 'trajectory.tsv') == row.trajectory_sha256
                tr = np.loadtxt(directory / 'trajectory.tsv', skiprows=1)
                np.testing.assert_array_equal(tr[:, 0], np.arange(50000, -1, -25))
                np.testing.assert_allclose(tr[:, 1], tr[:, 4] / (2 * tr[:, 2]), atol=1e-14)
                np.testing.assert_allclose(tr[:, 3], selection * ((tr[:, 0] < start * 1000) & (tr[:, 0] >= end * 1000)))
                assert np.isclose(tr[-1, 1], row.final_frequency) and 0 < tr[-1, 1] < 1
                trajectory = tr[:, 1]
            return scores, trajectory

        with ThreadPoolExecutor(max_workers=workers) as pool:
            loaded = list(pool.map(read, chosen.itertuples()))
        posterior = np.stack([item[0] for item in loaded])
        af = np.stack([item[1] for item in loaded[:100]])
        proof = json.loads((report / 'all_pairs_truth_provenance.json').read_text())
        assert proof['cdfs_sha256'] == inputs['all_pairs_truth_cdfs.npz']
        assert proof['completed_units_sha256'] == inputs['completed_units.tsv']
        with np.load(report / 'all_pairs_truth_cdfs.npz') as data:
            indices = {str(v): i for i, v in enumerate(data['ids'])}
            truth = data['cdfs']
            tree_hash = data['observed_sha256']
            for row in chosen.itertuples():
                assert tree_hash[indices[row.id]] == row.observed_sha256
            full = truth[[indices[i] for i in chosen.id]].copy()
            null_pool = (records.regime == '10k') & (records.origin == 'introgression')
            dev_ids = records[null_pool & (records.cohort == 'development')].id
            cal_ids = records[null_pool & (records.cohort == 'calibration')].id
            assert len(dev_ids) == 20 and len(cal_ids) == 99
            dev = np.diff(truth[[indices[i] for i in dev_ids]], axis=-1)
            cal = np.diff(truth[[indices[i] for i in cal_ids]], axis=-1)
            mu, sd = dev.mean(axis=(0, 1)), dev.std(axis=(0, 1), ddof=1)
            z = np.full(cal.shape, -np.inf)
            z[..., sd > 0] = (cal[..., sd > 0] - mu[sd > 0]) / sd[sd > 0]
            mus, sds, maxima = [mu], [sd], [z.max(axis=(1, 2))]
        with np.load(report / 'calibration.npz') as calibration:
            for method in METHODS[1:]:
                prefix = f'one_kyr_{regime}_introgression_' + method
                mus.append(calibration[prefix + '_mu'])
                sds.append(calibration[prefix + '_sd'])
                maxima.append(calibration[prefix + '_maxima'])
        np.savez_compressed(cache, ids=chosen.id.to_numpy(str), af=af,
                            years_ago=np.arange(50000, -1, -25),
                            cdfs=np.stack([full, posterior[:, 1], posterior[:, 2]]),
                            mu=np.stack(mus), sd=np.stack(sds), maxima=np.stack(maxima))
        save_json(manifest, dict(contract=contract, data_sha256=sha(cache)))
    with np.load(cache) as data:
        bundle = {name: data[name].copy() for name in data.files}
    assert bundle['cdfs'].shape == (3, 150, 100, 50)
    assert np.isfinite(bundle['cdfs']).all() and bundle['cdfs'].min() >= 0
    assert bundle['cdfs'].max() <= 1 + 1e-12 and np.diff(bundle['cdfs'], axis=-1).min() >= -1e-12
    assert bundle['af'].shape == (100, 2001) and np.isfinite(bundle['af']).all()
    assert bundle['af'].min() >= 0 and bundle['af'].max() <= 1
    np.testing.assert_allclose(bundle['af'][:, -1], selected.final_frequency, atol=1e-14)
    np.testing.assert_array_equal(bundle['ids'], chosen.id)
    chosen.to_csv(out / 'cohort_manifest.tsv', sep='\t', index=False, na_rep='NA', lineterminator='\n')
    return bundle, chosen


def calls_and_checks(bundle, report, out, regime='40_20k'):
    start, end, _ = EPISODES[regime]
    values = np.diff(bundle['cdfs'], axis=-1)
    calls, rows = [], []
    old = pd.read_csv(report / 'replicate_metrics.tsv', sep='\t')
    old_true = pd.read_csv(report / 'all_pairs_truth_replicates.tsv', sep='\t')
    for j, method in enumerate(METHODS):
        mu, sd = bundle['mu'][j], bundle['sd'][j]
        z = np.full(values[j].shape, -np.inf)
        z[..., sd > 0] = (values[j][..., sd > 0] - mu[sd > 0]) / sd[sd > 0]
        call = empirical_p(z, bundle['maxima'][j]) <= .05
        calls.append(call)
        expected = old_true if j == 0 else old[old.method == method]
        expected = expected[(expected.scan == 'one_kyr') & (expected.regime == regime) &
                            (expected.origin == 'introgression')].set_index('id')
        for i, ident in enumerate(bundle['ids']):
            metrics = evaluate(call[i], np.arange(100) * 10000, start, end, i < 100)
            for key, value in metrics.items():
                previous = expected.loc[ident, key]
                if value is None:
                    assert pd.isna(previous)
                elif isinstance(value, str):
                    assert value == previous, (ident, key)
                else:
                    assert np.isclose(value, previous, rtol=0, atol=1e-14), (ident, key)
            rows.append(dict(id=ident, method=method, cohort='selected' if i < 100 else 'neutral_test', **metrics))
    pd.DataFrame(rows).to_csv(out / 'verified_calls.tsv', sep='\t', index=False, na_rep='NA', lineterminator='\n')
    return np.stack(calls), pd.DataFrame(rows)


def figures(bundle, calls, metrics, out, selection=.02, regime='40_20k'):
    start, end, episode = EPISODES[regime]
    matplotlib.rcParams.update({'font.size': 14, 'axes.titlesize': 15, 'axes.labelsize': 15,
                               'xtick.labelsize': 12, 'ytick.labelsize': 12,
                               'pdf.fonttype': 42, 'ps.fonttype': 42})
    years = bundle['years_ago'] / 1000
    cdf, af = bundle['cdfs'], bundle['af']
    mass = np.diff(cdf, axis=-1)
    order = np.argsort(af[:, -1], kind='stable')
    summary = {'regime': regime, 'selection': selection, 'selected': 100, 'neutral_test': 50, 'methods': {}, 'af_milestones': {}}
    snapshots, profile_rows = [], []
    with PdfPages(out / 'selection_episode_series.pdf') as book:
        def save(fig, name):
            fig.savefig(out / (name + '.png'), dpi=180)
            fig.savefig(out / (name + '.pdf'))
            book.savefig(fig)
            plt.close(fig)

        fig, axes = plt.subplots(2, 1, figsize=(11, 8.5), gridspec_kw={'height_ratios': [1.4, 1]})
        fig.suptitle(f'AF histories · introgression, s={selection:g} · 100 survivors', y=.99, fontsize=17)
        axes[0].axvspan(end, start, color='.88', zorder=0)
        axes[0].plot(years, af.T, color=BLUE, alpha=.17, linewidth=.7)
        axes[0].plot(years, np.median(af, axis=0), color=ORANGE, linewidth=3)
        axes[0].set(xlim=(50, 0), ylim=(0, 1.02), xlabel='Time before present (kya)', ylabel='Focal allele frequency')
        axes[0].set_xticks([50, 40, 30, 20, 10, 0])
        axes[0].yaxis.set_major_formatter(PercentFormatter(1))
        fig.legend(handles=[Line2D([], [], color=BLUE, alpha=.5, label='Each of 100 replicates'),
                            Line2D([], [], color=ORANGE, lw=3, label='Median at each time'),
                            Patch(facecolor='.88', label=f'Selection: {episode}')],
                   loc='upper center', bbox_to_anchor=(.52, .95), ncol=3, frameon=False, fontsize=12)
        samples = []
        for age, _ in milestones(regime):
            freq = af[:, np.flatnonzero(years == age)[0]]
            samples.append(freq)
            quantiles = np.quantile(freq, [0, .1, .25, .5, .75, .9, 1])
            summary['af_milestones'][str(age)] = dict(zip(('min', 'q10', 'q25', 'median', 'q75', 'q90', 'max'), quantiles.tolist()))
            snapshots.extend(dict(id=i, time_kya=age, af=float(f)) for i, f in zip(bundle['ids'][:100], freq))
        axes[1].boxplot(samples, positions=np.arange(4), widths=.4, showfliers=False, whis=(10, 90),
                        medianprops={'color': ORANGE, 'linewidth': 2}, boxprops={'linewidth': 1.5})
        for k, freq in enumerate(samples):
            axes[1].scatter(k + np.linspace(-.2, .2, 100), freq[order], color=BLUE, alpha=.4, s=9)
        axes[1].set_xticks(np.arange(4), [f'{age} kya\n{label}' if age else label for age, label in milestones(regime)])
        axes[1].set(ylabel='Focal allele frequency', ylim=(0, 1.02))
        axes[1].yaxis.set_major_formatter(PercentFormatter(1))
        for ax in axes:
            ax.grid(axis='y', alpha=.2)
        fig.text(.5, .025, 'Whole-population AF; conditioned on segregation today. Boxes: IQR; whiskers: 10th–90th percentiles.', ha='center', fontsize=11)
        fig.subplots_adjust(top=.85, bottom=.13, left=.11, right=.98, hspace=.48)
        save(fig, '01_allele_frequency')

        fig, axes = plt.subplots(3, 2, figsize=(11, 8.5), sharex=True)
        fig.suptitle(f'Introgression, s={selection:g} · coalescence at the focal coordinate', y=.995, fontsize=17)
        for j, label in enumerate(LABELS):
            for col, (values, times) in enumerate(((cdf[j, :, 50], np.arange(1, 51)),
                                                   (mass[j, :, 50], np.arange(1.5, 50)))):
                ax = axes[j, col]
                ax.axvspan(end, start, color='.9', zorder=0)
                for cohort, sl, color in [('Neutral', slice(100, 150), BLUE), ('Selected', slice(0, 100), ORANGE)]:
                    v = values[sl]
                    lo, hi = np.quantile(v, [.1, .9], axis=0)
                    ax.fill_between(times, lo, hi, color=color, alpha=.14, linewidth=0)
                    ax.plot(times, v.mean(axis=0), color=color, lw=2)
                    profile_rows.extend(dict(method=METHODS[j], cohort=cohort, statistic='cdf' if col == 0 else 'one_kyr_mass',
                                             time_kya=float(t), mean=float(m), q10=float(l), q90=float(h))
                                        for t, m, l, h in zip(times, v.mean(axis=0), lo, hi))
                ax.set_title(label, loc='left', fontsize=13)
                ax.yaxis.set_major_formatter(PercentFormatter(1))
                ax.set(xlim=(1, 50), ylim=(0, None))
                if col == 0:
                    ax.set_ylim(0, 1.02)
                ax.set_xticks([1, 10, 20, 30, 40, 50])
                ax.grid(axis='y', alpha=.15)
                if j == 2:
                    ax.set_xlabel('Coalescence time (kya)')
                ax.set_ylabel('Pairs below cutoff' if col == 0 else 'Mass per 1-kyr bin', fontsize=13)
        fig.legend(handles=[Line2D([], [], color=ORANGE, lw=2, label='Selected: n=100'),
                            Line2D([], [], color=BLUE, lw=2, label='Neutral: n=50'),
                            Patch(facecolor='.9', label=f'Selection: {episode}')],
                   loc='upper center', bbox_to_anchor=(.52, .955), ncol=3, frameon=False, fontsize=12)
        fig.text(.5, .025, 'Lines: means. Bands: 10th–90th percentiles across regions, not confidence intervals.\nRight-hand y-scales vary. Mass is not renormalized to 1–50 kya; focal coordinate is 500 kb.', ha='center', fontsize=10.5)
        fig.subplots_adjust(top=.85, bottom=.15, left=.1, right=.98, wspace=.3, hspace=.55)
        save(fig, '02_selected_vs_neutral_tmrca')

        excess = mass[:, :100].mean(axis=1) - mass[:, 100:].mean(axis=1)
        limit = abs(excess).max()
        fig, axes = plt.subplots(3, 1, figsize=(11, 8.5), sharex=True, sharey=True)
        fig.suptitle(f'Introgression, s={selection:g} · selected-minus-neutral coalescence mass', y=.99, fontsize=16)
        for j, ax in enumerate(axes):
            im = ax.pcolormesh(np.arange(1, 51), np.arange(101) * 10 - 500, excess[j],
                               cmap='RdBu_r', vmin=-limit, vmax=limit, shading='flat', rasterized=False)
            for age in (end, start):
                ax.axvline(age, color='black', ls='--', lw=1.4)
            for pos in (-100, 100):
                ax.axhline(pos, color='black', ls=':', lw=1)
            ax.set_title(LABELS[j], loc='left', fontsize=14)
            ax.set_yticks([-500, 0, 500])
        axes[-1].set_xticks([1, 10, 20, 30, 40, 50])
        axes[-1].set_xlim(1, 50)
        axes[-1].set_xlabel('Coalescence time (kya)')
        fig.supylabel('Distance from selected coordinate (kb)', x=.025, fontsize=15)
        fig.text(.5, .935, f'Dashed: selection boundaries ({episode}). Dotted: ±100 kb. n=100 selected, 50 neutral.', ha='center', fontsize=11)
        fig.subplots_adjust(top=.87, bottom=.22, left=.13, right=.98, hspace=.43)
        cax = fig.add_axes([.28, .085, .57, .025])
        cb = fig.colorbar(im, cax=cax, orientation='horizontal', format=PercentFormatter(1))
        cb.set_label('Excess mass per 1-kyr bin', fontsize=13)
        cb.solids.set_rasterized(False)
        save(fig, '03_spatial_coalescence_excess')

        fig, axes = plt.subplots(1, 3, figsize=(11, 8.5), sharex=True, sharey=True)
        fig.suptitle(f'Introgression, s={selection:g} · temporal calls across 100 survivors', y=.99, fontsize=16)
        near = calls[:, :100, 40:61].any(axis=2)
        truth = (np.arange(1, 50) >= end) & (np.arange(1, 50) < start)
        for j, ax in enumerate(axes):
            codes = np.broadcast_to(truth, (100, 49)).astype(int).copy()
            codes[near[j, order]] = 2
            ax.pcolormesh(np.arange(1, 51), np.arange(101) + .5, codes,
                          cmap=ListedColormap(['white', '#e8e8e8', BLUE]), vmin=0, vmax=2,
                          shading='flat', rasterized=False)
            for age in (end, start):
                ax.axvline(age, color='black', ls='--', lw=1)
            selected = metrics[(metrics.method == METHODS[j]) & (metrics.cohort == 'selected')]
            null = metrics[(metrics.method == METHODS[j]) & (metrics.cohort == 'neutral_test')]
            info = dict(joint_hits=int(selected.joint_hit.sum()), only_mistimed=int(selected.near_mistimed.sum()),
                        any_outside_component=int(selected.near_off_target.sum()),
                        null_any_calls=int(null.any_call.sum()), center_hits=int(selected.center_hit.sum()),
                        mean_temporal_f1=float(selected.temporal_f1.mean()),
                        focal_mean_excess_peak_kya=float(1.5 + np.argmax(excess[j, 50])))
            summary['methods'][METHODS[j]] = info
            short = ('True TMRCAs\n19,900 pairs', 'Posterior median\n1,000 pairs', 'Posterior mass\n1,000 pairs')[j]
            ax.set_title(short + f"\n{info['joint_hits']} overlap; {info['only_mistimed']} mistimed", fontsize=13)
            ax.set(xlim=(1, 50), ylim=(.5, 100.5), xlabel='Coalescence time (kya)')
            ax.set_xticks([1, 10, 30, 50] if regime == '10k' else [1, 20, 40, 50])
            ax.set_yticks([1, 25, 50, 75, 100])
        axes[0].set_ylabel('Replicates ordered by present-day AF (low → high)', fontsize=14)
        fig.legend(handles=[Patch(facecolor=BLUE, label='Significant support'),
                            Patch(facecolor='#e8e8e8', label=f'True selection: {episode}')],
                   loc='upper center', bbox_to_anchor=(.53, .95), ncol=2, frameon=False, fontsize=12)
        fig.text(.5, .035, 'Support is the union over positions within ±100 kb; threshold is calibrated over the entire 1-Mb region.\n“Mistimed” means no near-site temporal overlap. Rows are identical across panels; every survivor is shown.',
                 ha='center', fontsize=10.5)
        fig.subplots_adjust(top=.79, bottom=.14, left=.11, right=.98, wspace=.16)
        save(fig, '04_called_intervals')
    pd.DataFrame(snapshots).to_csv(out / 'af_milestones.tsv', sep='\t', index=False, lineterminator='\n')
    pd.DataFrame(profile_rows).to_csv(out / 'tmrca_profiles.tsv', sep='\t', index=False, lineterminator='\n')
    pd.DataFrame({'rank_low_to_high': np.arange(1, 101), 'id': bundle['ids'][:100][order],
                  'present_af': af[order, -1]}).to_csv(out / 'plot_row_order.tsv', sep='\t', index=False, lineterminator='\n')
    save_json(out / 'summary.json', summary)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report', type=Path, default=ROOT / 'docs/results/selection_arrays_eas_20261105_4h')
    parser.add_argument('--out', type=Path)
    parser.add_argument('--selection', type=float, default=.02,
                        help='Selection coefficient of an existing complete 100-survivor cohort.')
    parser.add_argument('--regime', choices=tuple(EPISODES), default='40_20k')
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    start, end, episode = EPISODES[args.regime]
    assert 1 <= args.workers <= 24
    assert 0 < args.selection < 1
    if args.out is None:
        tag = format(args.selection, 'g').replace('.', '')
        args.out = ROOT / f'docs/results/selection_episode_{args.regime}_s{tag}_20261106'
    args.out.mkdir(parents=True, exist_ok=True)
    bundle, records = prepare(args.report, args.out, args.workers, args.selection, args.regime)
    calls, metrics = calls_and_checks(bundle, args.report, args.out, args.regime)
    print('All 450 region/method metric records match the frozen report.', flush=True)
    summary = figures(bundle, calls, metrics, args.out, args.selection, args.regime)
    counts = pd.read_csv(args.report / 'proposal_summary.tsv', sep='\t')
    count = counts[(counts.regime == args.regime) & (counts.origin == 'introgression') &
                   (counts.cohort == 'selected') & (counts.selection == args.selection)].iloc[0]
    summary['proposal_outcomes'] = {k: int(count[k]) for k in ('analyzed_prefix', 'accepted', 'lost', 'fixed')}
    save_json(args.out / 'summary.json', summary)
    pdf_checks = {}
    for path in args.out.glob('*.pdf'):
        data = path.read_bytes()
        check = dict(embedded_truetype=b'/FontFile2' in data,
                     raster_image_objects=len(re.findall(rb'/Subtype\s*/Image', data)),
                     pages=len(re.findall(rb'/Type\s*/Page\b', data)))
        assert check['embedded_truetype'] and check['raster_image_objects'] == 0
        assert check['pages'] == (4 if path.name == 'selection_episode_series.pdf' else 1)
        pdf_checks[path.name] = check
    save_json(args.out / 'pdf_validation.json', pdf_checks)
    lines = [f'# Selection at {episode}, s={args.selection:g}: complete introgression cohort', '',
        'Open [the four-page figure series](selection_episode_series.pdf). All pages are also saved individually as PNG and vector PDF; no inline plotting is used.', '',
        f'These are the first 100 present-day segregating survivors in the frozen proposal order, not a subset chosen for attractive trajectories or significant results. The 100 introgression survivors required {int(count.analyzed_prefix)} proposals: {int(count.lost)} lost the focal allele and {int(count.fixed)} fixed it. The null comparison uses the 50 held-out neutral introgression regions, without focal-allele survival or AF matching. Both cohorts have the same 1% Neanderthal pulse at 50 kya and EAS demographic history.', '',
        '## Figure guide', '',
        '1. [Allele-frequency trajectories](01_allele_frequency.pdf): all 100 whole-population trajectories, the pointwise median, and AF distributions at introgression, selection onset, selection cessation, and today. Selection is shaded; trajectory time runs from past to present. The selected cohort is conditioned on present-day segregation.',
        '2. [Selected versus neutral TMRCAs](02_selected_vs_neutral_tmrca.pdf): cumulative coalescence fractions and adjacent 1-kyr bin masses at the same fixed 500-kb coordinate in both cohorts. Lines are means; bands are pointwise 10th–90th percentiles across regions, not confidence or posterior credible intervals. Right-column y-scales vary by estimator. Mass is not conditioned on coalescence before 50 kya.',
        f'3. [Spatial and temporal excess](03_spatial_coalescence_excess.pdf): mean selected-minus-held-out-neutral coalescence mass across the 1-Mb region, sampled every 10 kb. All panels share the same color scale. Dashed lines mark selection boundaries at {end} and {start} kya when within the plotted range; dotted lines mark the ±100-kb detection tolerance. These are descriptive excesses, not significance maps.',
        '4. [Called intervals for every survivor](04_called_intervals.pdf): each row is one replicate, ordered by present-day AF from low to high and aligned across methods. Blue marks significant support at any scanned position within ±100 kb. Grey marks the true selected period. This union can combine support from different positions and does not preserve spatial connected-component identity.', '',
        '## Allele-frequency milestones', '',
        '| Time | Median AF | 10th–90th percentiles |', '|---|---:|---:|']
    for age, label in milestones(args.regime):
        values = summary['af_milestones'][str(age)]
        lines.append(f"| {age} kya: {label} | {values['median']:.2%} | {values['q10']:.2%}–{values['q90']:.2%} |")
    lines += ['', '## Timing and detection', '',
        '| Method | Overlapping near-site call / 100 | Only mistimed / 100 | Any wholly outside component / 100 | Center temporal F1 | Mean focal excess peak |',
        '|---|---:|---:|---:|---:|---:|']
    for method, label in zip(METHODS, LABELS):
        values = summary['methods'][method]
        lines.append(f"| {label} | {values['joint_hits']} | {values['only_mistimed']} | {values['any_outside_component']} | {values['mean_temporal_f1']:.3f} | {values['focal_mean_excess_peak_kya']:g} kya |")
    lines += ['', 'The table reports descriptive maxima of the mean selected-minus-neutral curve, not fitted onset/end dates or per-replicate significance calls. True TMRCAs use all 19,900 pairs; inference retains the original fixed 1,000 pairs, so the panels do not isolate inference error from pair subsampling. These are the original uncorrected Gamma-SMC estimates; no post-hoc temporal calibration or ASMC inference is applied.', '',
        'The primary scan covers [1,2), …, [49,50) kya; [0,1) is absent. Each method uses 20 independent development regions and 99 calibration regions, with conservative +1 empirical p-values against maxima over all 100 genomic positions and all 49 intervals. The 50 plotted neutral regions are held out from calibration. No thresholds were refitted to the selected cohort or the plotted null cohort. All 450 region/method metric records reproduce the frozen reports.', '',
        'A region can have both a correctly overlapping component and a wholly outside component. “Only mistimed” counts regions with near-site calls but no temporal overlap; “any wholly outside component” also includes errors alongside a correct call. Center temporal F1 refers to the exact focal coordinate, whereas the plotted call union uses ±100 kb.', '',
        '## Reproduce the figures without rerunning simulations', '',
        'The checked cache contains the plotted CDFs and all 100 AF trajectories. Cache validity depends on the source report, cohort identities, and file hashes, not Git commit identity. On a fresh clone, the saved cache supports rerendering without the large raw trees or posterior files. Initial cache construction requires the raw paths recorded in the cohort manifest. There are no RNG draws in this plotting workflow.', '',
        '```bash', '%%bash', 'set -euo pipefail', 'REPO=gamma_smc_ts_selection',
        'if [[ -d "$REPO/.git" ]]; then', '  git -C "$REPO" pull --ff-only', 'else',
        '  git clone --branch AOU_run_opt git@github.com:soisa001/gamma_smc_ts.git "$REPO"',
        'fi', 'cd "$REPO"',
        f'uv run --frozen python scripts/plot_selection_episode_series.py --regime {args.regime} --selection {args.selection:g} --workers 4', '```', '',
        'The script separates verified data extraction/cache reuse from plotting and metric validation. `series_data.npz` contains the full plotting arrays; `af_milestones.tsv` and `tmrca_profiles.tsv` contain tabular summaries; `plot_row_order.tsv` maps each displayed row to its simulation; `verified_calls.tsv` contains the recomputed metrics. `input_manifest.json`, `provenance.json`, `pdf_validation.json`, and `output_hashes.json` record provenance and integrity checks.']
    (args.out / 'README.md').write_text('\n'.join(lines) + '\n', encoding='utf-8', newline='\n')
    save_json(args.out / 'provenance.json', dict(script_sha256=sha(__file__), source_report=str(args.report.resolve()),
              workers=args.workers, parameters=vars(args) | {'report': str(args.report), 'out': str(args.out)},
              versions={n: __import__(n).__version__ for n in ('numpy', 'pandas', 'matplotlib')},
              calibration='20 development and 99 calibration neutral regions; maxima over 100 positions and 49 time bins; +1 empirical p<=.05.',
              verification='450 region/method metric records match frozen reports; original score/trajectory and true-tree cache hashes checked.',
              randomness='No simulations, bootstrap samples, or RNG draws. Plot offsets and cohort ordering are deterministic.',
              interpretation='100 selected survivors versus 50 unconditioned held-out neutral introgression regions. All curves are at fixed coordinates; no peak-based selection or alignment. Bands are between-region quantiles, not confidence intervals.',
              limits=f'Times 1–50 kya only; 0–1 kya absent from bin scan. All-pairs truth versus fixed-1000-pair inference; no exhaustive posterior decoding. Selection is {episode}, with a 1% Neanderthal pulse at 50 kya.'))
    shutil.copy2(__file__, args.out / 'plot_selection_episode_series.py')
    save_json(args.out / 'output_hashes.json', {p.name: sha(p) for p in sorted(args.out.iterdir()) if p.is_file() and p.name != 'output_hashes.json'})
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
