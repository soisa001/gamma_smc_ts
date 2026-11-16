#!/usr/bin/env python3
"""AF trajectories and per-simulation call maps for ongoing introgression s=.005."""
import os
os.environ['MPLBACKEND'] = 'Agg'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import argparse
import hashlib
import json
from pathlib import Path
import re
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.ticker import PercentFormatter
import numpy as np
from PIL import Image
from plot_cohort_significant_times import mass_order
from summarize_window_tpr_fpr import verify_ranks

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT/'docs/results'
METHODS = [('Posterior-weighted mass', 2, 'posterior_mass'), ('Posterior median', 1, 'posterior_median'), ('True TMRCA', 0, 'true_tmrca')]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=ROOT/'output/pdf/ongoing_50k_s0005_detail_20261116')
    args = parser.parse_args()
    provenance = {}

    def load(folder, name):
        path = RESULTS/folder/name
        inventory = json.loads((path.parent/'output_hashes.json').read_text())
        signature = sha(path); assert signature == inventory[name], str(path)
        provenance[str(path.relative_to(ROOT))] = signature
        if path.suffix == '.json': return json.loads(path.read_text())
        with np.load(path, allow_pickle=False) as f: return {k:f[k].copy() for k in f.files}

    base = 'denovo_comparison_20261111'
    meta = load(base, 'manifest.json')
    case = next(c for c in meta['cases'] if c['regime'] == '50k' and c['selection'] == .005)
    selected = load(base, meta['cohorts'][case['selected']['introgression']]['cache'])
    reference = load(base, meta['cohorts'][case['neutral']['introgression']]['cache'])
    bank = load('segregating_neutral_fpr_100_20261110', 'parameters.npz')
    metrics = load('window_tpr_fpr_20261116', 'summary.json')
    np.testing.assert_array_equal(reference['ids'], bank['reference_ids'])
    np.testing.assert_array_equal(reference['cdfs'], bank['reference_cdfs'])
    np.testing.assert_array_equal(reference['af'], bank['reference_af'])
    np.testing.assert_array_equal(selected['af_times_kya'], reference['af_times_kya'])
    test = dict(ids=bank['test_ids'], af=bank['test_af'], cdfs=bank['test_cdfs'],
        proposal_indices=bank['proposal_indices'], p_cumulative=bank['p_cumulative'], p_interval=bank['p_interval'])
    for data in (selected, reference, test):
        assert len(data['ids']) == len(set(data['ids'])) == 100
        assert data['af'].shape == (100, 2001)
        assert np.all((data['af'][:, -1] > 0) & (data['af'][:, -1] < 1))
    for data in (selected, test):
        assert set(data['ids']).isdisjoint(reference['ids'])
        verify_ranks(data, reference['cdfs'])
    assert set(selected['ids']).isdisjoint(test['ids'])
    orders = {'selected': mass_order(selected), 'neutral': mass_order(test)}
    example_index = 0  # Fixed input order; never choose an example by its signal.
    example_id = str(selected['ids'][example_index])
    example_label = example_id.replace('\\', '/').split('/')[-1]
    example_row = int(np.flatnonzero(orders['selected'] == example_index)[0]) + 1
    args.out.mkdir(parents=True, exist_ok=True)
    matplotlib.rcParams.update({'font.size': 16, 'axes.titlesize': 18, 'axes.labelsize': 17,
        'xtick.labelsize': 14, 'ytick.labelsize': 14, 'pdf.fonttype': 42, 'ps.fonttype': 42})
    orange, grey = '#D55E00', '#737373'
    pdfmeta = {'Title': '50 kya-present, s=0.005: trajectories and called times', 'CreationDate': None, 'ModDate': None}
    numerical = []
    with PdfPages(args.out/'all_statistics_detail.pdf', metadata=pdfmeta) as all_book, PdfPages(args.out/'posterior_mass_detail.pdf', metadata=pdfmeta) as mass_book:
        def emit(fig, stem, include_mass):
            fig.savefig(args.out/(stem+'.png'), dpi=180)
            all_book.savefig(fig)
            if include_mass: mass_book.savefig(fig)
            plt.close(fig)

        for one in (True, False):
            fig, ax = plt.subplots(figsize=(11, 8.5))
            fig.subplots_adjust(left=.12, right=.97, bottom=.18, top=.74)
            fig.suptitle('1% introgression | 50 kya-present | s = 0.005', y=.98, fontsize=21, fontweight='bold')
            subtitle = f'One selected simulation: {example_label}' if one else 'All 100 selected simulations'
            fig.text(.55, .925, subtitle, ha='center', fontsize=17)
            handles = [Line2D([], [], color=grey, lw=2, label='100 neutral references'),
                Line2D([], [], color=orange, lw=2, label='Selected simulation' if one else '100 selected simulations')]
            fig.legend(handles=handles, ncol=2, loc='upper center', bbox_to_anchor=(.55, .875), frameon=False, fontsize=14)
            ax.set_facecolor('#fafafa')
            ax.plot(reference['af_times_kya'], reference['af'].T, color=grey, lw=.7, alpha=.35)
            if one:
                ax.plot(selected['af_times_kya'], selected['af'][example_index], color=orange, lw=2.1)
            else:
                ax.plot(selected['af_times_kya'], selected['af'].T, color=orange, lw=.8, alpha=.4)
            ax.set(xlim=(50, 0), ylim=(0, 1), xlabel='Time before present (kya)', ylabel='Focal allele frequency')
            ax.yaxis.set_major_formatter(PercentFormatter(1)); ax.grid(alpha=.18)
            ax.spines[['top', 'right']].set_visible(False)
            fig.text(.55, .09, 'Selection is active throughout 50 kya-present. All focal alleles segregate today.', ha='center', fontsize=12)
            extra = f'First retained simulation; AF today = {100*selected["af"][example_index,-1]:.1f}%; call-map row {example_row}.' if one else 'Thin lines show individual trajectories; grey references also calibrate the coalescence tests.'
            fig.text(.55, .045, extra, ha='center', fontsize=12)
            emit(fig, '01_one_selected_af' if one else '02_all_selected_af', True)

        for method, m, method_key in METHODS:
            for cohort, data in [('selected', selected), ('neutral', test)]:
                order = orders[cohort]
                color = orange if cohort == 'selected' else '#7B3294'
                for alpha, tag in ((.05, 'p005'), (.01, 'p001')):
                    fig, axes = plt.subplots(1, 2, figsize=(11, 8.5))
                    fig.subplots_adjust(left=.10, right=.97, bottom=.20, top=.70, wspace=.20)
                    heading = '100 selected simulations' if cohort == 'selected' else '100 held-out neutral simulations'
                    fig.suptitle(heading+f' | p < {alpha:g}', y=.98, fontsize=21, fontweight='bold')
                    fig.text(.535, .928, method+' | focal site 500 kb', ha='center', fontsize=17)
                    subtitle = '50 kya-present | s = 0.005 | 1% introgression' if cohort == 'selected' else 'False calls in 50 kya-present | s = 0 | 1% introgression'
                    fig.text(.535, .877, subtitle, ha='center', fontsize=14)
                    fig.legend(handles=[Line2D([], [], color=color, lw=7, label='Called'), Line2D([], [], color='#eeeeee', lw=7, label='Not called')],
                        ncol=2, loc='upper center', bbox_to_anchor=(.535, .845), frameon=False, fontsize=13)
                    for column, statistic in enumerate(('cumulative', 'interval')):
                        ax = axes[column]
                        flags = data['p_'+statistic][:, m] < alpha
                        count = int(flags.any(axis=1).sum())
                        records = metrics['tpr'] if cohort == 'selected' else metrics['fpr']
                        match = next(r for r in records if r['origin']=='introgression' and r['regime']=='50k' and
                            r['method']==method and r['statistic']==statistic and r['alpha']==alpha and
                            (cohort=='neutral' or r['selection']==.005))
                        assert count == match['tp' if cohort=='selected' else 'fp']
                        edges = np.arange(.5, 51) if statistic == 'cumulative' else np.arange(1, 51)
                        mesh = ax.pcolormesh(edges, np.arange(.5, 101), flags[order].astype(int),
                            cmap=ListedColormap(['#f0f0f0', color]), vmin=0, vmax=1, shading='flat', rasterized=False)
                        mesh.set_rasterized(False)
                        label = 'Cumulative cutoff' if statistic == 'cumulative' else '1-kya interval'
                        ax.set(title=label, xlim=(0, 50.5 if statistic=='cumulative' else 50), ylim=(100.5, .5), xlabel='Time before present (kya)')
                        ax.set_xticks([0, 10, 20, 30, 40, 50]); ax.set_yticks([1, 25, 50, 75, 100])
                        if column == 0: ax.set_ylabel('Simulation (sorted)')
                        for spine in ax.spines.values(): spine.set_visible(False)
                        ax.tick_params(length=0)
                        for t in (10, 20, 30, 40): ax.axvline(t, color='white', lw=.6, alpha=.5)
                        if cohort == 'selected':
                            ax.plot([1.015], [example_row], marker='<', ms=7, color='black', transform=ax.get_yaxis_transform(), clip_on=False)
                        endpoint = 'TPR' if cohort=='selected' else 'FPR'
                        ax.text(.5, -.18, f'{endpoint}: {count}/100 ({count}%)', ha='center', transform=ax.transAxes, fontsize=16, fontweight='bold')
                        numerical.append(dict(method=method, cohort=cohort, statistic=statistic, alpha=alpha,
                            count=count, n=100, order_ids=data['ids'][order].tolist(),
                            called_times_by_row=[(np.flatnonzero(row)+1).tolist() for row in flags[order]]))
                    note = f'Black marker: AF example at row {example_row}. ' if cohort == 'selected' else 'Test neutrals are independent of the 100-reference bank. '
                    fig.text(.535, .085, note+'Row order fixed across methods and thresholds.', ha='center', fontsize=11)
                    fig.text(.535, .05, 'Rows sorted using posterior-mass calls at p<0.05. All 100 simulations are shown.', ha='center', fontsize=11)
                    fig.text(.535, .017, 'Cumulative cells mark cutoff coordinates; interval cells cover 1-2,...,49-50 kya.', ha='center', fontsize=10.5)
                    emit(fig, f'{method_key}_{cohort}_{tag}', method_key=='posterior_mass')
    for path in args.out.glob('*.png'):
        with Image.open(path) as im: assert im.size == (1980, 1530); im.verify()
    for name, pages in [('all_statistics_detail.pdf', 14), ('posterior_mass_detail.pdf', 6)]:
        content = (args.out/name).read_bytes()
        assert content.startswith(b'%PDF-') and content.rstrip().endswith(b'%%EOF')
        assert b'/FontFile2' in content and not re.search(rb'/Subtype\s*/Image', content)
        assert len(re.findall(rb'/Type\s*/Page\b', content)) == pages
    assert len(numerical) == 24 and len(list(args.out.glob('*.png'))) == 14
    (args.out/'.gitattributes').write_text('* -text\n')
    manifest = dict(source_sha256=dict(sorted(provenance.items())), script_sha256=sha(Path(__file__)),
        origin='introgression', selection=.005, selection_window_kya=[0, 50], example_id=example_id,
        example_rule='First retained simulation in the existing manifest/cache order.', example_call_map_row=example_row,
        example_final_af=float(selected['af'][example_index,-1]), selected_n=100, reference_n=100, independent_neutral_test_n=100,
        notes=['AF plots use the fixed reference neutrals. Neutral call maps use independent test neutrals.',
               'Rows sorted by posterior-mass p<.05 calls, consistently across every method and threshold.',
               'Calls are pointwise; TPR and FPR count any overlapping call inside 0-50 kya at the focal site.'],
        call_maps=numerical, validation=dict(independent_ranks_verified=True, summary_counts_verified=True,
            pdf_vector_only=True, editable_text=True, pngs=14, pdfs=2, pdf_pages=20))
    (args.out/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    (args.out/'output_hashes.json').write_text(json.dumps({p.name:sha(p) for p in sorted(args.out.iterdir()) if p.is_file() and p.name!='output_hashes.json'}, indent=2)+'\n')
    print(json.dumps(dict(example=example_label, row=example_row, validation=manifest['validation'],
        posterior_mass=[{k:r[k] for k in ('cohort','statistic','alpha','count')} for r in numerical if r['method']=='Posterior-weighted mass'])))


if __name__ == '__main__':
    main()
