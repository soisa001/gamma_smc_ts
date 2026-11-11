#!/usr/bin/env python3
"""Cache and compare surviving de novo versus introgressed selected alleles."""
import os
os.environ['MPLBACKEND'] = 'Agg'
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import shutil
import sys
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import PercentFormatter
import numpy as np
from plot_selection_episode_series import sha, save_json
from report_selection_arrays import local_path
from pointwise_time_tests import exceedance_counts
from plot_cohort_significant_times import mass_order

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / 'docs/results'
REFERENCE = RESULTS / 'selection_calls_100_segregating_null_20261110'
METHODS = ['True TMRCA', 'Posterior median', 'Posterior-weighted mass']
METHOD_KEYS = ['true_tmrca', 'posterior_median', 'posterior_mass']
LABELS = {'10k': '10 kya–present', '50k': '50 kya–present', '50_30k': '50–30 kya'}
ORIGINS = ['de_novo', 'introgression']
ORIGIN_LABELS = ['De novo', '1% introgression']
COLORS = ['#0072B2', '#D55E00']


def load_npz(path):
    with np.load(path, allow_pickle=False) as src:
        return {k: src[k].copy() for k in src.files}


def prepare(args):
    raw, out = args.raw.resolve(), args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    contract = json.loads((raw / 'manifest.json').read_text())
    ledger = json.loads((raw / 'proposal_ledger.json').read_text())
    progress = json.loads((raw / ('progress.json' if args.live else 'finished.json')).read_text())
    pairs = np.loadtxt(raw / 'pairs.tsv', dtype=int)
    assert pairs.shape == (1000, 2)
    entries = {}
    source_refs = {}

    def cache(key, data, provenance):
        path = out / (key + '.npz')
        np.savez_compressed(path, **data)
        entries[key] = dict(cache=path.name, cache_sha256=sha(path), n=len(data['ids']), provenance=provenance)

    def read_record(record):
        directory = local_path(record['source_directory'])
        for name, field in [('trajectory.tsv', 'trajectory_sha256'), ('scores.npz', 'scores_sha256'),
                            ('observed.trees', 'observed_sha256'), ('posterior.zst', 'posterior_sha256')]:
            assert sha(directory / name) == record[field], str(directory / name)
        assert record['origin'] == 'de_novo' and 0 < record['final_frequency'] < 1
        trajectory = np.loadtxt(directory / 'trajectory.tsv', skiprows=1)
        years = np.arange(record['ancestry_start_years'], -1, -25)
        assert trajectory.shape == (len(years), 5)
        np.testing.assert_array_equal(trajectory[:, 0], years)
        np.testing.assert_allclose(trajectory[:, 1], trajectory[:, 4] / (2 * trajectory[:, 2]), atol=1e-14)
        assert trajectory[0, 4] == 1 and 0 < trajectory[-1, 1] < 1
        onset = json.loads((directory / 'onset_provenance.json').read_text())
        assert onset['donor_haplotypes'] == 0 and onset['selected_initial_copies'] == 1
        ticks = record['selected_reproduction_ticks']
        if isinstance(ticks, str): ticks = json.loads(ticks)
        gen = np.arange(1, len(years) + 1)
        first = 1 if ticks[0] == 2 else ticks[0]
        np.testing.assert_allclose(trajectory[:, 3], record['selection'] * ((gen >= first) & (gen <= ticks[1])), atol=1e-12)
        posterior_meta = json.loads((directory / 'posterior.zst.meta').read_text())
        np.testing.assert_array_equal(posterior_meta['pairs'], pairs)
        score = load_npz(directory / 'scores.npz')
        np.testing.assert_array_equal(score['positions'], np.arange(100) * 10000)
        np.testing.assert_array_equal(score['cutoffs'], np.arange(1000, 50001, 1000))
        np.testing.assert_array_equal(score['valid_pairs'], 1000)
        cdf = score['cdfs'][:, 50]
        assert np.isfinite(cdf).all() and np.all(np.diff(cdf) >= -1e-12)
        return trajectory[:, 1], cdf

    for stratum in contract['strata']:
        key = stratum['key']
        if args.live:
            rows = (contract['reused'][key] + [r for r in ledger[key] if r['accepted']])[:contract['target']]
        else:
            receipt = json.loads((raw / (key + '_cohort.json')).read_text())
            assert receipt['manifest_sha256'] == sha(raw / 'manifest.json')
            rows = receipt['accepted']
        assert len({r['observed_sha256'] for r in rows}) == len(rows)
        assert all(r['selection'] == stratum['selection'] for r in rows)
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            values = list(pool.map(read_record, rows))
        age = stratum['start']
        ids = [str(local_path(r['source_directory'])) for r in rows]
        data = dict(af=np.stack([v[0] for v in values]) if values else np.empty((0, age // 25 + 1)),
                    cdfs=np.stack([v[1] for v in values]) if values else np.empty((0, 3, 50)),
                    ids=np.array(ids, dtype=str), af_times_kya=np.arange(age, -1, -25) / 1000,
                    proposal_indices=np.arange(len(rows)), pairs=pairs)
        cache('de_novo_' + key.replace('.', 'p'), data, dict(stratum=stratum, records=rows))
        print(f'Cached de novo {key}: {len(rows)}', flush=True)

    ref_meta = json.loads((REFERENCE / 'manifest.json').read_text())
    ref_af_meta = json.loads((REFERENCE / 'neutral_af_manifest.json').read_text())
    assert sha(REFERENCE / 'neutral.npz') == ref_meta['neutral_cache_sha256']
    assert sha(REFERENCE / 'neutral_af.npz') == ref_af_meta['cache_sha256']
    ref = load_npz(REFERENCE / 'neutral.npz'); af = load_npz(REFERENCE / 'neutral_af.npz')
    np.testing.assert_array_equal(ref['ids'], af['ids']); np.testing.assert_array_equal(ref['pairs'], pairs)
    cache('introgression_neutral', dict(cdfs=ref['cdfs'], ids=ref['ids'], af=af['trajectories'][:, :, 1],
        af_times_kya=af['trajectories'][0, :, 0] / 1000, pairs=pairs, proposal_indices=np.arange(100)),
        dict(source=str(REFERENCE), cdf_sha256=sha(REFERENCE / 'neutral.npz'), af_sha256=sha(REFERENCE / 'neutral_af.npz')))
    cases = []
    for s in (r for r in contract['strata'] if r['role'] == 'selected'):
        key = s['key'].replace('.', 'p')
        if s['selection'] == .005:
            directory = RESULTS / 'selection_s0005_20261110'; oldkey = key
        elif s['selection'] == .01:
            directory = REFERENCE; oldkey = s['regime'] + '_s001'
        elif s['selection'] == .05:
            directory = RESULTS / 'strong_selection_20261110'; oldkey = key
        else: raise ValueError('No frozen introgression comparison for requested s')
        meta = json.loads((directory / 'manifest.json').read_text())
        entry = next((e for e in meta['cases'] if e['case']['key'] == oldkey), None)
        source_refs[str(directory)] = sha(directory / 'manifest.json')
        if entry is None:
            accounting = json.loads((directory / 'proposal_accounting.json').read_text())['strata']
            fate = next(r for r in accounting if r['key'] == s['key'])
            assert fate['retained'] == 0
            data = dict(cdfs=np.empty((0, 3, 50)), af=np.empty((0, 2001)), ids=np.array([], dtype=str), proposal_indices=np.array([], dtype=int))
            provenance = dict(source=str(directory), proposal_accounting=fate)
        else:
            path = directory / entry['cache']; assert sha(path) == entry['cache_sha256']
            old = load_npz(path)
            data = {k: old[k] for k in ('cdfs', 'af', 'ids', 'proposal_indices')}
            provenance = dict(source=str(path), sha256=sha(path))
        data.update(af_times_kya=np.arange(50000, -1, -25) / 1000, pairs=pairs)
        cache('introgression_' + key, data, provenance)
        neutral_key = contract['neutral_mapping'][s['regime']]
        cases.append(dict(key=key, regime=s['regime'], selection=s['selection'], window=[s['end']/1000, s['start']/1000],
            label=f"{LABELS[s['regime']]} · s = {s['selection']:g}",
            selected={'de_novo': 'de_novo_' + key, 'introgression': 'introgression_' + key},
            neutral={'de_novo': 'de_novo_' + neutral_key, 'introgression': 'introgression_neutral'}))
    for case in cases:
        for origin in ORIGINS:
            key = case['selected'][origin]; data = load_npz(out / entries[key]['cache'])
            null = load_npz(out / entries[case['neutral'][origin]]['cache'])
            assert set(data['ids']).isdisjoint(null['ids'])
            if len(null['ids']):
                counts = exceedance_counts(data['cdfs'], null['cdfs'])
                data.update(p_cumulative=(counts[0]+1)/(len(null['ids'])+1), p_interval=(counts[1]+1)/(len(null['ids'])+1))
            else:
                data.update(p_cumulative=np.full((len(data['ids']), 3, 50), np.nan), p_interval=np.full((len(data['ids']), 3, 49), np.nan))
            np.savez_compressed(out / entries[key]['cache'], **data)
            entries[key]['cache_sha256'] = sha(out / entries[key]['cache'])
    save_json(out / 'manifest.json', dict(cases=cases, cohorts=entries, simulation=contract, progress=progress,
        raw_manifest_sha256=sha(raw / 'manifest.json'), source_references=source_refs, live_snapshot=args.live,
        position_bp=500000, scan_stride_bp=10000, cutoffs_kya=list(range(1, 51)), pair_count=1000,
        p_definition='(1 + reference values >= query) / (n_reference + 1), strict p < alpha; own-origin, mutation-age-matched survivor null.',
        interpretation='Any-time detection at the focal site uses pointwise tests, not scan-wide calibrated power. No independent de novo FPR cohort.',
        conditioning=contract['conditioning'], script_sha256=sha(Path(__file__)), python=sys.version))
    save_json(out / 'proposal_accounting.json', dict(strata=progress['strata'],
        note='New-proposal accounting excludes reused survivors. Zero retained means not estimable, not zero power.'))
    for name in ('schedule.json', 'execution_environment.json'):
        shutil.copy2(raw / name, out / ('simulation_' + name))
    for name in ('founder_cache_validation.json', 'previous_phase_finished.json'):
        if (raw / name).exists(): shutil.copy2(raw / name, out / name)
    (out / '.gitattributes').write_text('* -text\n')


def report(args):
    out = args.out.resolve(); meta = json.loads((out / 'manifest.json').read_text())
    data = {}
    for key, entry in meta['cohorts'].items():
        path = out / entry['cache']; assert sha(path) == entry['cache_sha256']
        data[key] = load_npz(path)
    matplotlib.rcParams.update({'font.size': 15, 'axes.labelsize': 15, 'axes.titlesize': 16,
        'xtick.labelsize': 13, 'ytick.labelsize': 13, 'pdf.fonttype': 42, 'ps.fonttype': 42})
    summary = []; orders = {}; hashes = {}
    def save(fig, stem, book):
        for ext in ('png', 'pdf'):
            path = out / f'{stem}.{ext}'; fig.savefig(path, dpi=150); hashes[path.name] = sha(path)
        book.savefig(fig); plt.close(fig)
    def clean(ax):
        ax.grid(alpha=.15)
        ax.spines[['top', 'right']].set_visible(False)
    def header(fig, case, subtitle):
        fig.suptitle(case['label'], y=.985, fontsize=21)
        fig.text(.5, .93, subtitle, ha='center', fontsize=14)
    def footer(fig, text):
        fig.text(.5, .035, text, ha='center', fontsize=11)
    def boundary(values, alpha):
        n = len(values); max_exceed = np.flatnonzero((1 + np.arange(n+1))/(n+1) < alpha)
        if not len(max_exceed): return None
        return np.sort(values, axis=0)[-(int(max_exceed[-1])+1)]
    def format_count(mask, n, nref, alpha):
        if not n: return 'No segregating selected replicates'
        if 1/(nref+1) >= alpha: return f'Insufficient null: n = {nref}'
        return f'{mask.any(1).sum()}/{n} with ≥1 pointwise call'

    null_groups = [('de_novo_10k_neutral', 'De novo, 10 kya', '#0072B2'),
                   ('de_novo_50k_neutral', 'De novo, 50 kya', '#009E73'),
                   ('introgression_neutral', 'Introgressed, 50 kya', '#D55E00')]
    with PdfPages(out / 'neutral_reference_comparison.pdf') as book:
        fig, axes = plt.subplots(1, 3, figsize=(11, 8.5), sharey=True)
        fig.subplots_adjust(left=.10, right=.97, top=.77, bottom=.17, wspace=.16)
        fig.suptitle('Segregating neutral allele-frequency trajectories', y=.97, fontsize=20)
        fig.text(.5, .90, 'Each line is one neutral simulation; all original focal alleles survive to the present', ha='center', fontsize=13)
        for ax, (key, label, color) in zip(axes, null_groups):
            null = data[key]
            ax.plot(null['af_times_kya'], null['af'].T, color=color, alpha=.35, lw=.7)
            ax.set(title=f"{label}\nn = {len(null['ids'])}", xlim=(50, 0), ylim=(0, 1), xlabel='Time (kya)')
            ax.set_xticks([50, 25, 0]); ax.yaxis.set_major_formatter(PercentFormatter(1)); clean(ax)
        axes[0].set_ylabel('Focal allele frequency')
        footer(fig, 'De novo starts as one copy; introgression starts at 1%. EAS demography; whole-population segregation.')
        save(fig, 'neutral_af_comparison', book)
        for method, title in enumerate(METHODS):
            fig, axes = plt.subplots(1, 2, figsize=(11, 8.5))
            fig.subplots_adjust(left=.11, right=.97, top=.73, bottom=.16, wspace=.30)
            fig.suptitle('Neutral references: ' + title, y=.98, fontsize=20)
            handles = [Line2D([], [], color=color, lw=2, label=f"{label} (n={len(data[key]['ids'])})") for key, label, color in null_groups]
            handles += [Line2D([], [], color='black', lw=2, label='Solid: cohort median'),
                        Line2D([], [], color='black', ls='--', label='Dashed: p<0.05 boundary'),
                        Line2D([], [], color='black', ls=':', label='Dotted: p<0.01 boundary')]
            fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.5, .91), ncol=2, frameon=False, fontsize=12)
            for col, ax in enumerate(axes):
                x = np.arange(1, 51) if col==0 else np.arange(1.5, 50)
                for key, label, color in null_groups:
                    values = data[key]['cdfs'][:, method]
                    if col: values = np.diff(values)
                    ax.plot(x, values.T, color=color, alpha=.13, lw=.6)
                    if len(values): ax.plot(x, np.median(values, axis=0), color=color, lw=2)
                    cutoff = boundary(values, .05)
                    if cutoff is not None: ax.plot(x, cutoff, color=color, ls='--', lw=1.8)
                    cutoff = boundary(values, .01)
                    if cutoff is not None: ax.plot(x, cutoff, color=color, ls=':', lw=1.8)
                ax.set(title='Cumulative cutoff' if col==0 else '1-kya interval', xlim=(1, 50), ylim=(0, None),
                       xlabel='Cutoff / bin midpoint (kya)', ylabel='Pair fraction')
                ax.yaxis.set_major_formatter(PercentFormatter(1)); clean(ax)
            footer(fig, 'Focal site: 500 kb. Thin lines show every retained neutral. Pointwise thresholds, without scan correction.')
            save(fig, 'neutral_' + METHOD_KEYS[method] + '_comparison', book)

    with PdfPages(out / 'allele_frequency_comparison.pdf') as book:
        for case in meta['cases']:
            fig, axes = plt.subplots(1, 2, figsize=(11, 8.5), sharey=True)
            fig.subplots_adjust(top=.70, bottom=.17, left=.10, right=.97, wspace=.18)
            header(fig, case, 'Original focal allele segregating today in selected and neutral simulations')
            handles = [Line2D([], [], color='.5', label='Neutral trajectories'), Line2D([], [], color=COLORS[0], label='De novo selected'),
                       Line2D([], [], color=COLORS[1], label='Introgressed selected'), Patch(color='.92', label='Selection active')]
            fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.5, .885), ncol=2, frameon=False, fontsize=13)
            for ax, origin, label, color in zip(axes, ORIGINS, ORIGIN_LABELS, COLORS):
                sel, null = data[case['selected'][origin]], data[case['neutral'][origin]]
                ax.axvspan(*case['window'], color='.92')
                ax.plot(null['af_times_kya'], null['af'].T, color='.45', alpha=.35, lw=.65)
                ax.plot(sel['af_times_kya'], sel['af'].T, color=color, alpha=.4, lw=.8)
                ax.set(title=f"{label}\n{len(sel['ids'])} selected; {len(null['ids'])} neutrals", xlim=(50, 0), ylim=(0, 1), xlabel='Time before present (kya)')
                if not len(sel['ids']): ax.text(.5, .7, 'No surviving selected cohort', ha='center', transform=ax.transAxes, fontsize=13)
                ax.yaxis.set_major_formatter(PercentFormatter(1)); clean(ax)
            axes[0].set_ylabel('Focal allele frequency')
            footer(fig, 'De novo: one copy at selection onset. Introgression: 1% at 50 kya. Both use EAS demography.')
            save(fig, case['key'] + '_af_comparison', book)
            if case['regime']=='10k' and case['selection']<.05:
                sel, null = data[case['selected']['de_novo']], data[case['neutral']['de_novo']]
                fig, ax = plt.subplots(figsize=(11, 8.5))
                fig.subplots_adjust(top=.76, bottom=.16, left=.12, right=.97)
                header(fig, case, 'De novo AF zoom · one copy introduced at 10 kya')
                ax.plot(null['af_times_kya'], null['af'].T, color='.45', alpha=.4, lw=.7)
                ax.plot(sel['af_times_kya'], sel['af'].T, color=COLORS[0], alpha=.4, lw=.85)
                maximum = max([.01] + [float(d['af'].max()) for d in (sel, null) if d['af'].size])
                ax.set(xlim=(10, 0), ylim=(0, min(1, maximum*1.05)), xlabel='Time before present (kya)', ylabel='Focal allele frequency')
                ax.yaxis.set_major_formatter(PercentFormatter(1)); clean(ax)
                fig.legend(handles=[Line2D([], [], color=COLORS[0], label=f"{len(sel['ids'])} selected trajectories"),
                    Line2D([], [], color='.45', label=f"{len(null['ids'])} neutral trajectories")],
                    loc='upper center', bbox_to_anchor=(.5, .87), ncol=2, frameon=False, fontsize=14)
                footer(fig, 'Expanded AF scale; every retained trajectory is shown. Original focal alleles segregate today.')
                save(fig, case['key'] + '_denovo_af_zoom', book)

    with PdfPages(out / 'coalescence_comparison.pdf') as book, PdfPages(out / 'posterior_mass_comparison.pdf') as massbook:
        for case in meta['cases']:
            for method, title in enumerate(METHODS):
                fig, axes = plt.subplots(2, 2, figsize=(11, 8.5), sharex=True)
                fig.subplots_adjust(top=.76, bottom=.15, left=.11, right=.97, wspace=.35, hspace=.55)
                header(fig, case, title + ' · focal site 500 kb · 1,000 fixed pairs')
                fig.legend(handles=[Line2D([], [], color='.5', label='Each neutral'), Line2D([], [], color=COLORS[0], label='Each de novo selected'),
                    Line2D([], [], color=COLORS[1], label='Each introgressed selected'),
                    Line2D([], [], color='black', ls='--', label='Null p < 0.05 boundary'),
                    Line2D([], [], color='black', ls=':', label='Null p < 0.01 boundary')],
                    loc='upper center', bbox_to_anchor=(.5, .885), ncol=3, frameon=False, fontsize=11.5)
                for row, (origin, label, color) in enumerate(zip(ORIGINS, ORIGIN_LABELS, COLORS)):
                    sel, null = data[case['selected'][origin]], data[case['neutral'][origin]]
                    for col, kind in enumerate(['cumulative', 'interval']):
                        ax = axes[row, col]; x = np.arange(1, 51) if col == 0 else np.arange(1.5, 50)
                        y, n = sel['cdfs'][:, method], null['cdfs'][:, method]
                        if col: y, n = np.diff(y), np.diff(n)
                        ax.axvspan(*case['window'], color='.94')
                        ax.plot(x, n.T, color='.5', alpha=.25, lw=.6)
                        ax.plot(x, y.T, color=color, alpha=.3, lw=.7)
                        cutoff = boundary(n, .05)
                        if cutoff is not None: ax.plot(x, cutoff, color='black', ls='--', lw=1.5)
                        cutoff = boundary(n, .01)
                        if cutoff is not None: ax.plot(x, cutoff, color='black', ls=':', lw=1.5)
                        ax.set(title=f'{label}: {kind}', xlim=(1, 50), ylim=(0, None), ylabel='Pair fraction')
                        ax.yaxis.set_major_formatter(PercentFormatter(1)); clean(ax)
                        ax.text(.02, .96, f"Selected n={len(sel['ids'])}; null n={len(null['ids'])}", va='top', transform=ax.transAxes, fontsize=10)
                        if row == 1: ax.set_xlabel('Cutoff / bin midpoint (kya)')
                for col in range(2):
                    ymax = max(ax.get_ylim()[1] for ax in axes[:, col])
                    for ax in axes[:, col]: ax.set_ylim(0, ymax)
                footer(fig, 'Matched y scales by column. Thresholds are pointwise; shading shows the selection window.')
                if method == 2: massbook.savefig(fig)
                save(fig, case['key'] + '_' + METHOD_KEYS[method] + '_curves', book)

    for alpha, tag in [(.05, 'p005'), (.01, 'p001')]:
        with PdfPages(out / f'calls_comparison_{tag}.pdf') as book, PdfPages(out / f'posterior_mass_calls_{tag}.pdf') as massbook:
            for case in meta['cases']:
                for method, title in enumerate(METHODS):
                    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
                    fig.subplots_adjust(top=.79, bottom=.17, left=.12, right=.97, wspace=.24, hspace=.65)
                    header(fig, case, f'{title} · pointwise p < {alpha:g} · focal site 500 kb')
                    for row, (origin, label, color) in enumerate(zip(ORIGINS, ORIGIN_LABELS, COLORS)):
                        sel, null = data[case['selected'][origin]], data[case['neutral'][origin]]
                        n, nref = len(sel['ids']), len(null['ids'])
                        order = mass_order(sel) if n and method == 2 else np.arange(n)
                        orders[case['key'] + ':' + origin + ':' + str(method)] = sel['ids'][order].tolist()
                        for col, kind in enumerate(['cumulative', 'interval']):
                            ax = axes[row, col]; mask = sel['p_' + kind][:, method] < alpha
                            ax.axvspan(*case['window'], color='.92')
                            for display, i in enumerate(order):
                                bins = np.flatnonzero(mask[i])
                                if col:
                                    ax.broken_barh([(int(b)+1, 1) for b in bins], (display+.65, .7), facecolors=color)
                                else: ax.scatter(bins+1, np.full(len(bins), display+1), c=color, s=1)
                            ax.set(title=f'{label}: {kind}', xlim=(0, 50), ylim=(max(n, 1)+.8, .2), xlabel='Time before present (kya)')
                            ax.set_yticks(np.unique(np.linspace(1, n, min(n, 4)).astype(int)) if n else [])
                            if col == 0: ax.set_ylabel('Simulation (sorted)' if method == 2 else 'Simulation')
                            ax.text(.5, -.42, format_count(mask, n, nref, alpha) + (f' · null n={nref}' if n and nref and 1/(nref+1)<alpha else ''), ha='center', transform=ax.transAxes, fontsize=10)
                            clean(ax)
                            attainable = nref > 0 and 1/(nref+1) < alpha
                            summary.append(dict(case=case['key'], regime=case['regime'], selection=case['selection'], origin=origin,
                                method=title, statistic=kind, alpha=alpha, selected_n=n, neutral_n=nref,
                                threshold_attainable=attainable, any_time_count=int(mask.any(1).sum()) if n and attainable else None,
                                any_time_fraction=float(mask.any(1).mean()) if n and attainable else None,
                                median_final_af=float(np.median(sel['af'][:, -1])) if n else None,
                                fraction_by_time=mask.mean(0).tolist() if n and attainable else None))
                    footer(fig, 'Own-origin survivor null. Mass rows sorted at p<0.05. Pointwise calls do not estimate selection duration.')
                    if method == 2: massbook.savefig(fig)
                    save(fig, case['key'] + '_' + METHOD_KEYS[method] + '_calls_' + tag, book)
    with PdfPages(out / 'detection_comparison.pdf') as book:
        for alpha, tag in [(.05, 'p005'), (.01, 'p001')]:
            for method, title in enumerate(METHODS):
                fig, axes = plt.subplots(1, 2, figsize=(11, 8.5), sharey=True)
                fig.subplots_adjust(left=.23, right=.97, top=.80, bottom=.16, wspace=.18)
                fig.suptitle(f'{title}: pointwise p < {alpha:g}', y=.98, fontsize=21)
                fig.legend(handles=[Line2D([], [], marker='o', ls='', color=c, label=l)
                           for c, l in zip(COLORS, ORIGIN_LABELS)], loc='upper center', bbox_to_anchor=(.6, .91), ncol=2, frameon=False)
                for col, kind in enumerate(['cumulative', 'interval']):
                    ax = axes[col]
                    for i, case in enumerate(meta['cases']):
                        for offset, origin, color in zip([-.15, .15], ORIGINS, COLORS):
                            row = next(r for r in summary if r['case']==case['key'] and r['origin']==origin and r['method']==title and r['statistic']==kind and r['alpha']==alpha)
                            rate = row['any_time_fraction']
                            if rate is None:
                                ax.text(.015, i+offset, 'NA', color=color, va='center', fontsize=10)
                            else:
                                ax.plot(rate, i+offset, 'o', color=color, ms=6)
                    ax.set(title=kind.capitalize(), xlim=(-.03, 1.04), ylim=(len(meta['cases'])-.5, -.6), xlabel='Fraction with any call')
                    ax.xaxis.set_major_formatter(PercentFormatter(1)); ax.set_xticks([0, .5, 1])
                    ax.set_yticks(range(len(meta['cases'])))
                    if col == 0: ax.set_yticklabels([f"{LABELS[c['regime']]}\ns = {c['selection']:g}" for c in meta['cases']], fontsize=14)
                    clean(ax)
                footer(fig, 'Focal-site detection, without scan-wide correction. NA: no survivors or insufficient neutral resolution.')
                save(fig, METHOD_KEYS[method] + '_detection_' + tag, book)
    save_json(out / 'summary.json', dict(results=summary, interpretation=meta['interpretation'], report_script_sha256=sha(Path(__file__)),
        unattainable='No estimate when no selected survivors or insufficient neutral references for the strict p threshold.'))
    save_json(out / 'row_orders.json', orders)
    save_json(out / 'output_hashes.json', {p.name: sha(p) for p in sorted(out.iterdir())
        if p.is_file() and p.name != 'output_hashes.json' and p.suffix != '.log'})
    print(json.dumps([r for r in summary if r['origin']=='de_novo' and r['method']=='Posterior-weighted mass'], indent=2), flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('phase', choices=['prepare', 'report', 'all'])
    p.add_argument('--raw', type=Path, default=ROOT.parent / 'sim/denovo_comparison_cached_20261111')
    p.add_argument('--out', type=Path, default=RESULTS / 'denovo_comparison_20261111')
    p.add_argument('--workers', type=int, default=4)
    p.add_argument('--live', action='store_true', help='Explicit provisional snapshot of completed proposal prefixes.')
    args = p.parse_args()
    assert 1 <= args.workers <= 24
    if args.phase in ('prepare', 'all'): prepare(args)
    if args.phase in ('report', 'all'): report(args)


if __name__ == '__main__': main()
