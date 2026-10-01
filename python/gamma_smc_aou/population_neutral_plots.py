"""Large-text PNG/vector-PDF summaries of pointwise neutral cutoffs."""
import os
os.environ['MPLBACKEND'] = 'Agg'
import json
import zipfile
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'pdf.fonttype':42, 'ps.fonttype':42, 'font.size':17,
    'axes.titlesize':18, 'axes.labelsize':17, 'xtick.labelsize':14, 'ytick.labelsize':14})
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from . import origin_onset as oo
from .population_neutral import POPULATIONS


def plot(root, scores, cutoffs, counts):
    figures = root/'figures'
    figures.mkdir(exist_ok=True)
    pdf_path = root/'neutral_cutoffs.pdf'
    names = []
    times = sorted(cutoffs[cutoffs.cutoff_years > 0].cutoff_years.unique())
    provisional = bool(cutoffs.provisional.any()) or len(cutoffs.population.unique()) < len(POPULATIONS)
    scope = 'PROVISIONAL: incomplete null cohorts' if provisional else '1,000 null regions per population'
    with PdfPages(pdf_path) as pdf:
        def save(fig, name):
            for ax in fig.axes:
                for collection in ax.collections: collection.set_rasterized(False)
            fig.savefig(figures/(name+'.png'), dpi=160, bbox_inches='tight')
            fig.savefig(figures/(name+'.pdf'), bbox_inches='tight')
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)
            names.append(name)
        for source in ('truth', 'decoded'):
            fig, axes = plt.subplots(1, 2, figsize=(14,8), layout='constrained')
            for ax, method in zip(axes, ('all_pairs', 'alt_alt')):
                g = cutoffs[(cutoffs.source == source)&(cutoffs.method == method)]
                matrix = g.pivot(index='population', columns='cutoff_years', values='critical_score').reindex(index=POPULATIONS, columns=times).to_numpy()
                display = np.where(np.isposinf(matrix), 1., matrix)
                display = np.where(np.isneginf(display), 0., display)
                cmap = plt.get_cmap('viridis').copy()
                cmap.set_bad('#ececec')
                mesh = ax.imshow(display, vmin=0, vmax=1, cmap=cmap, aspect='auto')
                ax.set_xticks(range(len(times)), [t//1000 for t in times])
                totals = counts.set_index('population').completed.to_dict()
                ax.set_yticks(range(len(POPULATIONS)), [f'{p} (n={totals[p]})' for p in POPULATIONS])
                ax.set(xlabel='TMRCA cutoff (kya)', title=method.replace('_', ' ').upper())
                for i in range(len(POPULATIONS)):
                    for j in range(len(times)):
                        v = matrix[i,j]
                        label = 'pending' if np.isnan(v) else ('no cutoff' if v >= 1 else f'>{v:.4f}')
                        ax.text(j, i, label, ha='center', va='center', fontsize=12,
                                color='black' if np.isnan(v) or v > .55 else 'white')
            bar = fig.colorbar(mesh, ax=axes, fraction=.025, pad=.025, label='Raw recent-pair fraction boundary')
            bar.solids.set_rasterized(False)
            fig.suptitle(f'{source.capitalize()} TMRCA | site-level p ≤ 0.05\n{scope}', fontsize=21)
            save(fig, 'cutoffs_'+source)
        for cutoff in times:
            fig, axes = plt.subplots(2, 2, figsize=(14,11), layout='constrained')
            for i, source in enumerate(('truth','decoded')):
                for j, method in enumerate(('all_pairs','alt_alt')):
                    ax = axes[i,j]
                    g = scores[(scores.source==source)&(scores.method==method)&(scores.cutoff_years==cutoff)]
                    ax.boxplot([g[g.population==p].score.dropna().to_numpy() for p in POPULATIONS],
                        tick_labels=POPULATIONS, showfliers=True)
                    for x, p in enumerate(POPULATIONS, 1):
                        row = cutoffs[(cutoffs.population==p)&(cutoffs.source==source)&
                                      (cutoffs.method==method)&(cutoffs.cutoff_years==cutoff)]
                        if len(row) and np.isfinite(row.critical_score.iloc[0]):
                            ax.hlines(row.critical_score.iloc[0], x-.4, x+.4, color='#b62828', ls='--', lw=2)
                    ax.set(ylim=(-.03,1.03), title=f'{source} | {method.replace("_", " ")}',
                           ylabel='Raw P(TMRCA < T)')
            fig.suptitle(f'T = {cutoff//1000} kya | dashed red: p ≤ 0.05 boundary (must exceed)\n{scope}', fontsize=19)
            save(fig, f'null_distributions_T{cutoff}')
        fig, ax = plt.subplots(figsize=(11,8.5), layout='constrained')
        af = scores[scores.method=='AF']
        ax.boxplot([af[af.population==p].score.to_numpy() for p in POPULATIONS], tick_labels=POPULATIONS)
        for x,p in enumerate(POPULATIONS,1):
            row = cutoffs[(cutoffs.population==p)&(cutoffs.method=='AF')]
            if len(row) and np.isfinite(row.critical_score.iloc[0]):
                ax.hlines(row.critical_score.iloc[0], x-.4, x+.4, color='#b62828', ls='--', lw=2)
        ax.set(ylim=(-.03,1.03), ylabel='Present-day focal ALT frequency',
            title='Neutral allele frequencies | dashed red: upper 5% boundary\n'+scope)
        save(fig, 'allele_frequency')
    outputs = [root/name for name in ('pointwise_cutoffs.csv','null_scores.csv','cohort_counts.csv','manifest.json')]
    outputs += [pdf_path]+[figures/(name+ext) for name in names for ext in ('.png','.pdf')]
    oo.legacy.atomic_json(root/'report_complete.json', dict(provisional=provisional, figures=len(names),
        calibrated_regions=int(counts.completed.sum()), planned_regions=int(counts.planned.sum()),
        artifacts=oo.artifact_specs(root, [str(p.relative_to(root)) for p in outputs])))
    with zipfile.ZipFile(root/'neutral_cutoffs_bundle.zip', 'w', compression=zipfile.ZIP_DEFLATED) as archive:
        for path in outputs+[root/'report_complete.json']:
            archive.write(path, path.relative_to(root))
