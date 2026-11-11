#!/usr/bin/env python3
"""Plot every cached survivor-conditioned reference neutral; no simulations."""
import os
os.environ['MPLBACKEND'] = 'Agg'
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'

import argparse
import hashlib
import json
from pathlib import Path
import sys

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.ticker import PercentFormatter
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=ROOT / 'docs/results/selection_calls_100_segregating_null_20261110')
    parser.add_argument('--out', type=Path, default=ROOT / 'docs/results/neutral_reference_100_20261111')
    args = parser.parse_args()
    source, out = args.source.resolve(), args.out.resolve()
    meta = json.loads((source / 'manifest.json').read_text())
    af_meta = json.loads((source / 'neutral_af_manifest.json').read_text())
    assert meta['neutral_conditioned'] and meta['position_bp'] == 500000
    assert sha(source / 'neutral.npz') == meta['neutral_cache_sha256']
    assert sha(source / 'neutral_af.npz') == af_meta['cache_sha256']
    with np.load(source / 'neutral.npz', allow_pickle=False) as data:
        cdf, ids, pairs = (data[k].copy() for k in ('cdfs', 'ids', 'pairs'))
    with np.load(source / 'neutral_af.npz', allow_pickle=False) as data:
        trajectory = data['trajectories'].copy()
        np.testing.assert_array_equal(ids, data['ids'])
    assert cdf.shape == (100, 3, 50) and trajectory.shape == (100, 2001, 5)
    assert len(set(ids)) == 100 and len(pairs) == 1000
    assert np.isfinite(cdf).all() and np.isfinite(trajectory).all()
    assert ((cdf >= 0) & (cdf <= 1)).all() and (np.diff(cdf, axis=-1) >= -1e-12).all()
    np.testing.assert_array_equal(trajectory[:, :, 0], np.broadcast_to(np.arange(50000, -1, -25), (100, 2001)))
    np.testing.assert_array_equal(trajectory[:, :, 3], 0)
    np.testing.assert_allclose(trajectory[:, :, 1], trajectory[:, :, 4] / (2 * trajectory[:, :, 2]), atol=1e-14)
    final_af = trajectory[:, -1, 1]
    assert ((final_af > 0) & (final_af < 1)).all()
    interval = np.diff(cdf, axis=-1)
    counts = np.rint(cdf[:, :2] * len(pairs)).astype(int)
    np.testing.assert_allclose(cdf[:, :2] * len(pairs), counts, atol=1e-8)
    interval[:, :2] = np.diff(counts, axis=-1) / len(pairs)
    np.testing.assert_allclose(interval.sum(axis=-1), cdf[:, :, -1] - cdf[:, :, 0], atol=1e-12)
    out.mkdir(parents=True, exist_ok=True)
    matplotlib.rcParams.update({'font.size': 16, 'axes.titlesize': 17, 'axes.labelsize': 16,
                               'xtick.labelsize': 14, 'ytick.labelsize': 14,
                               'pdf.fonttype': 42, 'ps.fonttype': 42})
    norm = Normalize(0, .6)
    cmap = matplotlib.colormaps['viridis']
    colors = cmap(norm(final_af))
    order = np.argsort(final_af, kind='stable')
    outputs = []

    def clean(ax):
        ax.grid(alpha=.16)
        ax.yaxis.set_major_formatter(PercentFormatter(1))
        ax.spines[['top', 'right']].set_visible(False)

    def colorbar(fig, axes):
        bar = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axes,
                           orientation='horizontal', fraction=.065, pad=.15, aspect=45)
        bar.set_label('Present-day focal allele frequency (same color in every panel)')
        bar.ax.xaxis.set_major_formatter(PercentFormatter(1))
        bar.solids.set_rasterized(False)

    def save(fig, stem, book):
        for ext in ('png', 'pdf'):
            path = out / f'{stem}.{ext}'
            fig.savefig(path, dpi=160)
            outputs.append(path)
        book.savefig(fig)
        plt.close(fig)

    with PdfPages(out / 'all_100_neutrals.pdf') as book:
        fig, ax = plt.subplots(figsize=(11, 8.5))
        fig.subplots_adjust(left=.12, right=.96, top=.80, bottom=.12)
        fig.suptitle('100 neutral allele-frequency trajectories', y=.97, fontsize=21)
        fig.text(.5, .91, 'EAS demography | 1% Neanderthal introduction at 50 kya | s = 0', ha='center')
        fig.text(.5, .86, 'Original focal allele segregating in the present-day population', ha='center', fontsize=14)
        for i in order:
            ax.plot(trajectory[i, :, 0] / 1000, trajectory[i, :, 1], color=colors[i], lw=.85, alpha=.75)
        assert len(ax.lines) == 100
        ax.set(xlim=(50, 0), ylim=(0, 1), xlabel='Time before present (kya)', ylabel='Focal allele frequency')
        clean(ax)
        colorbar(fig, ax)
        save(fig, 'allele_frequency_100_neutrals', book)

        names = [('true_tmrca', 'True TMRCA'), ('posterior_median', 'Posterior median'), ('posterior_mass', 'Posterior-weighted mass')]
        handles = [Line2D([], [], color=cmap(.55), lw=1.5, label='100 individual neutrals'),
                   Line2D([], [], color='black', lw=2, label='Cohort median'),
                   Line2D([], [], color='black', lw=1.7, ls='--', label='Pointwise p < 0.05 boundary')]
        for method, (stem, title) in enumerate(names):
            fig, axes = plt.subplots(2, 1, figsize=(11, 8.5))
            fig.subplots_adjust(left=.13, right=.97, top=.78, bottom=.10, hspace=.70)
            fig.suptitle(f'100 neutrals: {title}', y=.98, fontsize=21)
            fig.text(.5, .925, 'Focal site: 500 kb | 1,000 fixed pairs | Present-day segregation required', ha='center', fontsize=14)
            fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.5, .90), ncol=3, fontsize=11.5, frameon=False)
            for ax, values, times, label in zip(axes, [cdf[:, method], interval[:, method]],
                    [np.arange(1, 51), np.arange(1.5, 50)],
                    ['Cumulative: TMRCA below cutoff', 'Interval: TMRCA within each 1-kya bin']):
                for i in order:
                    ax.plot(times, values[i], color=colors[i], lw=.85, alpha=.65)
                assert len(ax.lines) == 100
                ax.plot(times, np.median(values, axis=0), color='black', lw=2)
                ax.plot(times, np.sort(values, axis=0)[-5], color='black', lw=1.7, ls='--')
                ax.set(title=label, xlim=(1, 50), ylim=(0, None), ylabel='Pair fraction')
                ax.set_xlabel('TMRCA cutoff (kya)' if ax is axes[0] else 'TMRCA bin midpoint (kya)')
                clean(ax)
            colorbar(fig, axes)
            save(fig, f'{stem}_100_neutrals', book)
    outputs.append(out / 'all_100_neutrals.pdf')
    receipt = dict(source=str(source), output=str(out), script_sha256=sha(Path(__file__)),
        input_hashes={p: sha(source / p) for p in ('manifest.json', 'neutral.npz', 'neutral_af_manifest.json', 'neutral_af.npz')},
        ids=ids.tolist(), n=100, position_bp=500000, pair_count=len(pairs), selection=0,
        conditioning='Original introduced focal allele has 0 < whole-population AF < 1 at present.',
        cumulative_cutoffs_kya=list(range(1, 51)), interval_edges_kya=list(range(1, 51)),
        interval_note='49 bins from 1-2 through 49-50 kya; no inferred 0-1 kya bin.',
        threshold='Fifth-largest reference value; independent query must strictly exceed it for (1 + count(reference >= query))/101 < 0.05. Pointwise, not scan-wide; not self-testing neutrals.',
        plotting=dict(size_inches=[11, 8.5], dpi=160, color='viridis, whole-population final AF, fixed range 0 to 0.6',
                      order='ascending final AF, stable', random_draws=False),
        final_af_quantiles=dict(zip(['min', 'q25', 'median', 'q75', 'max'], np.quantile(final_af, [0, .25, .5, .75, 1]).tolist())),
        versions=dict(python=sys.version, numpy=np.__version__, matplotlib=matplotlib.__version__),
        outputs={p.name: sha(p) for p in outputs})
    (out / 'manifest.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(dict(out=str(out), figures=4, pdf_pages=4, validated_replicates=100, final_af=receipt['final_af_quantiles']), indent=2))


if __name__ == '__main__':
    main()
