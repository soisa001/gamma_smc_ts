"""Empirical focal-site CDFs with equal-weight, complete demographic-history groups."""
from __future__ import annotations
import os
os.environ['MPLBACKEND']='Agg'
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42,'font.size':15})
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from . import origin_onset as oo
from . import population_neutral as pn


def history_cdfs(groups, grid):
    """No-call nulls remain in the denominator, below the bounded score support."""
    return np.stack([np.searchsorted(np.sort(np.where(np.isnan(v),-np.inf,v)),grid,side='right')/len(v)
                     for v in groups])


def report(root):
    root=Path(root)
    manifest=json.loads((root/'manifest.json').read_text())
    expected=manifest['regions_per_history']
    out=root/'cdf'
    out.mkdir(exist_ok=True)
    records=[]
    inputs={}
    for pop in pn.POPULATIONS:
        for path in sorted((root/pop/'regions').glob('I50/null/rep*/calibration.json')):
            r=json.loads(path.read_text())
            cfg=r['identity']['simulation']['parameters']
            history=cfg['demographic_draw']['history_id']
            for source,sha in r['scores_sha256'].items():
                assert oo.legacy.digest(root/pop/'scores/regions'/r['task']['id']/source)==sha
            inputs[str(path.relative_to(root))]=oo.legacy.digest(path)
            for score in r['scores']:
                records.append(dict(population=pop,history_id=history,replicate=r['task']['replicate'],**score))
    if not records:
        return
    data=pd.DataFrame(records)
    data.to_csv(out/'focal_scores.csv',index=False)
    cuts=[]
    bandrows=[]
    for source in ('truth','decoded'):
        for method in ('all_pairs','alt_alt'):
            figures={t:plt.subplots(2,3,figsize=(17,10),layout='constrained') for t in (10000,50000)}
            for pi,pop in enumerate(pn.POPULATIONS):
                individual,iaxes=plt.subplots(1,2,figsize=(13,6),layout='constrained')
                for ti,t in enumerate((10000,50000)):
                    subset=data[(data.population==pop)&(data.source==source)&(data.method==method)&(data.cutoff_years==t)]
                    counts=subset.groupby('history_id').size()
                    complete=counts[counts==expected].index
                    balanced=subset[subset.history_id.isin(complete)]
                    axes=(figures[t][1].flat[pi],iaxes[ti])
                    if balanced.empty:
                        for ax in axes:
                            ax.text(.5,.5,'Awaiting complete history groups',ha='center',va='center',transform=ax.transAxes,fontsize=11)
                            ax.set(title=f'{pop}: T = {t:,} years',xlim=(0,1),ylim=(0,1))
                        continue
                    groups=[g.score.to_numpy(dtype=float) for _,g in balanced.groupby('history_id')]
                    values=balanced.score.to_numpy(dtype=float)
                    # Use every observed score as a knot; quantile steps are exact.
                    grid=np.unique(np.r_[0.,values[np.isfinite(values)],1.])
                    cdfs=history_cdfs(groups,grid)
                    lo,med,hi=np.quantile(cdfs,[.025,.5,.975],axis=0)
                    mean=cdfs.mean(axis=0)
                    sorted_values=np.sort(np.where(np.isnan(values),-np.inf,values))
                    q95,q99=np.quantile(sorted_values,[.95,.99],method='inverted_cdf')
                    boundaries={a:pn.critical_summary(values,a)['critical_score'] for a in (.05,.01)}
                    cuts.append(dict(population=pop,source=source,method=method,cutoff_years=t,n=len(values),
                        complete_histories=len(groups),regions_per_history=expected,available=int(np.isfinite(values).sum()),
                        q95=q95,q99=q99,p05_strict_boundary=boundaries[.05],p01_strict_boundary=boundaries[.01],
                        decision='target score > boundary; missing target is no call',provisional=len(values)<manifest['replicates_per_population']))
                    bandrows.extend(dict(population=pop,source=source,method=method,cutoff_years=t,score=x,
                        mean_cdf=a,median_cdf=b,q025_cdf=c,q975_cdf=d,complete_histories=len(groups))
                        for x,a,b,c,d in zip(grid,mean,med,lo,hi))
                    for ax in axes:
                        ax.fill_between(grid,lo,hi,step='post',color='#1683cd',alpha=.18,label='95% between-history range')
                        ax.step(grid,mean,where='post',color='#1476c9',lw=2,label='Mean CDF')
                        ax.step(grid,med,where='post',color='#17446e',ls='--',lw=1.6,label='Median CDF')
                        ax.axhline(.95,color='.55',ls='--',lw=.8)
                        ax.axhline(.99,color='.55',ls=':',lw=.8)
                        for q in (q95,q99):
                            if np.isfinite(q):
                                y=np.mean(sorted_values<=q)
                                ax.scatter([q],[y],color='#ef713b',edgecolor='white',s=55,zorder=5)
                        ax.set(title=f'{pop} | {len(groups)} histories, n = {len(values)}\nTop 5%: {q95:.1%}   Top 1%: {q99:.1%}',
                               xlim=(0,1),ylim=(0,1.015),xlabel=f'Fraction of pairs with TMRCA < {t:,} years',ylabel='Fraction of null focal sites at or below')
                        ax.xaxis.set_major_formatter(PercentFormatter(1))
                        ax.yaxis.set_major_formatter(PercentFormatter(1))
                        ax.grid(alpha=.18)
                handles,labels=iaxes[0].get_legend_handles_labels()
                if handles:
                    individual.legend(handles,labels,loc='outside upper center',ncol=3,frameon=False,fontsize=12)
                individual.savefig(out/f'{pop}_{source}_{method}.png',dpi=180)
                individual.savefig(out/f'{pop}_{source}_{method}.pdf')
                plt.close(individual)
            for t,(fig,axes) in figures.items():
                handles,labels=axes.flat[0].get_legend_handles_labels()
                if handles:
                    fig.legend(handles,labels,loc='outside upper center',ncol=3,frameon=False)
                fig.savefig(out/f'all_{source}_{method}_T{t}.png',dpi=150)
                fig.savefig(out/f'all_{source}_{method}_T{t}.pdf')
                plt.close(fig)
    pd.DataFrame(cuts).to_csv(out/'four_cutoffs_long.csv',index=False)
    pd.DataFrame(bandrows).to_csv(out/'cdf_curves_and_history_bands.csv',index=False)
    oo.legacy.atomic_json(out/'manifest.json',dict(inputs=inputs,
        quantiles='empirical inverted CDF, no interpolation; exact finite-null p boundaries also supplied',
        balancing='Only complete history groups are used for curves and cutoffs; all completed scores saved separately',
        bands='pointwise 2.5/97.5 percentiles of history-specific CDFs; include finite-region simulation noise; not simultaneous or mean-CDF confidence intervals',
        missing_scores='below valid score support in ranking and CDF denominator; missing target is no call'))
