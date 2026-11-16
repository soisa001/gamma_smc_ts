#!/usr/bin/env python3
"""Render window-matched TPR/FPR as annotated, vector-safe letter-size heatmaps."""
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
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle
from matplotlib.cm import ScalarMappable
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
METHODS = [('Posterior-weighted mass', 'posterior_mass'), ('Posterior median', 'posterior_median'), ('True TMRCA', 'true_tmrca')]
REGIMES = ['10k', '50k', '50_30k']
S = [.005, .01, .05, .1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=ROOT/'docs/results/window_tpr_fpr_20261116')
    parser.add_argument('--out', type=Path, default=ROOT/'output/pdf/window_tpr_fpr_heatmaps_20261116')
    args = parser.parse_args()
    signature = sha(args.source/'summary.json')
    inventory = json.loads((args.source/'output_hashes.json').read_text())
    assert signature == inventory['summary.json']
    src = json.loads((args.source/'summary.json').read_text())
    tpr = {(r['origin'], r['regime'], r['selection'], r['method'], r['statistic'], r['alpha']): r for r in src['tpr']}
    fpr = {(r['regime'], r['method'], r['statistic'], r['alpha']): r for r in src['fpr']}
    assert len(tpr) == 252 and len(fpr) == 36
    args.out.mkdir(parents=True, exist_ok=True)
    matplotlib.rcParams.update({'font.size': 14, 'axes.titlesize': 16, 'xtick.labelsize': 13,
        'ytick.labelsize': 14, 'pdf.fonttype': 42, 'ps.fonttype': 42})
    cmaps = {'TPR': matplotlib.colormaps['Blues'], 'FPR': matplotlib.colormaps['Oranges']}
    norm = Normalize(0, 100)
    cell_records = []
    pdfmeta = {'Title': 'Window-matched selection TPR and neutral FPR', 'CreationDate': None, 'ModDate': None}
    with PdfPages(args.out/'all_statistics_heatmaps.pdf', metadata=pdfmeta) as all_book, PdfPages(args.out/'posterior_mass_heatmaps.pdf', metadata=pdfmeta) as mass_book:
        for method, stem in METHODS:
            for alpha, tag in ((.05, 'p005'), (.01, 'p001')):
                fig = plt.figure(figsize=(11, 8.5))
                fig.suptitle(f'{method} | p < {alpha:g}', x=.55, y=.98, fontsize=21, fontweight='bold')
                fig.text(.55, .927, 'Calls overlapping the same time window at the focal site', ha='center', fontsize=14)
                for label, left in [('TPR', .27), ('FPR', .65)]:
                    axbar = fig.add_axes([left, .835, .23, .018])
                    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmaps[label]), cax=axbar, orientation='horizontal', ticks=[0, 25, 50, 75, 100])
                    cb.solids.set_rasterized(False)
                    cb.ax.tick_params(labelsize=11, length=2)
                    cb.ax.set_title(label+' (%)', fontsize=13, pad=7)
                    cb.outline.set_visible(False)
                grid = fig.add_gridspec(2, 2, left=.17, right=.98, top=.735, bottom=.175, hspace=.58, wspace=.13)
                for oi, origin in enumerate(('de_novo', 'introgression')):
                    for si, statistic in enumerate(('cumulative', 'interval')):
                        ax = fig.add_subplot(grid[oi, si])
                        origin_label = 'De novo' if origin == 'de_novo' else '1% introgression'
                        stat_label = 'Cumulative' if statistic == 'cumulative' else '1-kya interval'
                        ax.set_title(origin_label+' | '+stat_label, pad=12)
                        centers = [.5, 1.5, 2.5, 3.5, 4.85]
                        for y, regime in enumerate(REGIMES):
                            for x, s in enumerate(S + [None]):
                                endpoint = 'FPR' if s is None else 'TPR'
                                left = x if s is not None else 4.35
                                row = fpr.get((regime, method, statistic, alpha)) if s is None and origin == 'introgression' else (
                                    tpr.get((origin, regime, s, method, statistic, alpha)) if s is not None else None)
                                n = row['neutral_n'] if row and s is None else row['selected_n'] if row else None
                                count = row['fp'] if row and s is None else row['tp'] if row else None
                                value = 100*count/n if n else None
                                bg = cmaps[endpoint](norm(value)) if value is not None else '#eeeeee'
                                ax.add_patch(Rectangle((left, y), 1, 1, facecolor=bg, edgecolor='white', linewidth=2, rasterized=False))
                                if value is None:
                                    label = 'NA\nn=0' if n == 0 else 'No run' if s is not None else 'NA'
                                    ax.text(left+.5, y+.5, label, ha='center', va='center', fontsize=12, color='#626262')
                                else:
                                    small = n < 100
                                    label = f'{value:.0f}%'+('*' if small else '')
                                    color = 'white' if value >= 65 else '#17252e'
                                    ax.text(left+.5, y+(.37 if small else .5), label, ha='center', va='center', fontsize=14 if small else 15, fontweight='bold', color=color)
                                    if small:
                                        ax.text(left+.5, y+.72, f'{count}/{n}', ha='center', va='center', fontsize=11, color=color)
                                        ax.add_patch(Rectangle((left+.025, y+.025), .95, .95, fill=False, edgecolor='#555555', linewidth=1, linestyle='--'))
                                cell_records.append(dict(method=method, alpha=alpha, origin=origin, regime=regime,
                                    statistic=statistic, endpoint=endpoint, selection=s, n=n, count=count, percent=value))
                        ax.set(xlim=(-.02, 5.37), ylim=(3.02, -.02))
                        ax.set_xticks(centers, ['.005', '.01', '.05', '.10', 'FPR'])
                        ax.set_yticks([.5, 1.5, 2.5], ['10 kya-\npresent', '50 kya-\npresent', '50-30 kya'] if si == 0 else ['', '', ''])
                        ax.tick_params(axis='both', length=0, pad=7)
                        for spine in ax.spines.values(): spine.set_visible(False)
                fig.text(.57, .115, 'TPR columns: selection coefficient s | FPR: held-out neutrals in the same window', ha='center', fontsize=12)
                fig.text(.55, .073, 'n=100 per cell unless marked. *Small cohort: count/n shown; NA n=0: no survivors.', ha='center', fontsize=11)
                fig.text(.55, .044, 'No run: treatment not simulated. Independent de novo FPR is unavailable.', ha='center', fontsize=11)
                fig.text(.55, .015, 'EAS | Present-day segregation | Cumulative: cutoff within window; interval: positive overlap.', ha='center', fontsize=10.5)
                filename = f'{stem}_{tag}.png'
                fig.savefig(args.out/filename, dpi=180)
                all_book.savefig(fig)
                if stem == 'posterior_mass': mass_book.savefig(fig)
                plt.close(fig)
    # Each source result appears exactly once; unavailable cells stay unavailable.
    observed_tpr = [r for r in cell_records if r['endpoint'] == 'TPR' and r['n'] is not None]
    observed_fpr = [r for r in cell_records if r['endpoint'] == 'FPR' and r['n'] is not None]
    assert len(observed_tpr) == len(tpr) and len(observed_fpr) == len(fpr)
    for r in observed_tpr:
        original = tpr[(r['origin'], r['regime'], r['selection'], r['method'], r['statistic'], r['alpha'])]
        assert r['count'] == original['tp'] and r['n'] == original['selected_n']
    for r in observed_fpr:
        original = fpr[(r['regime'], r['method'], r['statistic'], r['alpha'])]
        assert r['count'] == original['fp'] and r['n'] == 100
    for path in args.out.glob('*.png'):
        with Image.open(path) as im:
            assert im.size == (1980, 1530); im.verify()
    for name, pages in [('all_statistics_heatmaps.pdf', 6), ('posterior_mass_heatmaps.pdf', 2)]:
        content = (args.out/name).read_bytes()
        assert content.startswith(b'%PDF-') and content.rstrip().endswith(b'%%EOF')
        assert b'/FontFile2' in content and not re.search(rb'/Subtype\s*/Image', content)
        assert len(re.findall(rb'/Type\s*/Page\b', content)) == pages
    (args.out/'.gitattributes').write_text('* -text\n')
    manifest = dict(source=str((args.source/'summary.json').relative_to(ROOT)), source_sha256=signature,
        script_sha256=sha(Path(__file__)), color_scales=dict(TPR=[0, 100], FPR=[0, 100]),
        cells=cell_records, validation=dict(tpr_cells=len(observed_tpr), fpr_cells=len(observed_fpr),
        pngs=6, vector_pdfs=2, pdf_pages=8, editable_text=True, raster_images_in_pdf=False))
    (args.out/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    (args.out/'output_hashes.json').write_text(json.dumps({p.name:sha(p) for p in sorted(args.out.iterdir()) if p.is_file() and p.name!='output_hashes.json'}, indent=2)+'\n')
    print(json.dumps(manifest['validation']))


if __name__ == '__main__':
    main()
