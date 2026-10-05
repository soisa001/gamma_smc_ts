"""Read verified completed calibration receipts; export ECDFs without running sims."""
import os
os.environ['MPLBACKEND'] = 'Agg'
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[name] = '1'
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import zipfile
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42,'font.size':17,
    'axes.labelsize':18,'axes.titlesize':19,'xtick.labelsize':15,'ytick.labelsize':15})
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import PercentFormatter
from gamma_smc_aou.population_neutral import POPULATIONS, critical_summary


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_record(path):
    record = json.loads(path.read_text())
    pop = record['population']
    task = record['task']
    assert task['role']=='null' and task['s']==0 and task['ascertainment_years']==50000
    root = path.parents[5]
    frames = []
    for source in ('truth','decoded'):
        score_path = root/pop/'scores'/'regions'/task['id']/(source+'.csv')
        assert digest(score_path)==record['scores_sha256'][source+'.csv'], score_path
        frame = pd.read_csv(score_path,float_precision='round_trip')
        expected = pd.DataFrame([s for s in record['scores'] if s['source']==source])
        keys=['source','method','cutoff_years','position_0based','n_pairs']
        a=frame.sort_values(keys).reset_index(drop=True)
        b=expected.sort_values(keys).reset_index(drop=True)
        pd.testing.assert_frame_equal(a[keys],b[keys])
        np.testing.assert_allclose(a.score,b.score.astype(float),rtol=0,atol=1e-14,equal_nan=True)
        frame['population']=pop
        frame['task_id']=task['id']
        frames.append(frame[frame.cutoff_years.isin([10000,50000])])
    return pd.concat(frames),dict(path=str(path),sha256=digest(path),population=pop,task_id=task['id'])


def ecdf(values):
    x=np.sort(np.asarray(values,dtype=float))
    x=x[np.isfinite(x)]
    return np.r_[0,x,1]*100,np.r_[0,np.arange(1,len(x)+1)/len(x),1]*100


def summaries(scores):
    rows=[]
    for key,g in scores.groupby(['population','source','method','cutoff_years']):
        available=g.score.dropna().to_numpy()
        for alpha in (.05,.01):
            q=float(np.quantile(available,1-alpha,method='inverted_cdf')) if len(available) else np.nan
            exact=critical_summary(g.score,alpha)
            boundary=exact['critical_score']
            possible=exact['rejection_possible_for_bounded_score']
            reason=('available' if possible else 'insufficient_null_count' if len(g)+1 < 1/alpha else 'bounded_score_or_ties')
            rows.append(dict(zip(['population','source','method','tmrca_years'],key),
                alpha=alpha,null_n=len(g),available_n=len(available),missing_n=len(g)-len(available),
                empirical_percentile_pct=q*100,
                exact_p_boundary_pct=boundary*100 if np.isfinite(boundary) else np.nan,
                exact_p_attainable=possible,exact_p_status=reason,
                minimum_p=1/(len(g)+1),decision='finite score strictly greater than boundary',provisional=True))
    return pd.DataFrame(rows)


def panel(ax, scores, cuts, pop, source, method, time):
    selected=scores[(scores.source==source)&(scores.method==method)&(scores.cutoff_years==time)]
    for other in POPULATIONS:
        values=selected[selected.population==other].score.dropna()
        if len(values):
            ax.step(*ecdf(values),where='post',color='#c5c4ba',lw=1.2,zorder=1)
    group=selected[selected.population==pop]
    available=group.score.dropna()
    if len(available): ax.step(*ecdf(available),where='post',color='#2074c9',lw=2.8,zorder=3)
    info=cuts[(cuts.population==pop)&(cuts.source==source)&(cuts.method==method)&(cuts.tmrca_years==time)]
    labels=[]
    for alpha in (.05,.01):
        row=info[info.alpha==alpha].iloc[0]
        q=row.empirical_percentile_pct
        level=100*(1-alpha)
        ax.axhline(level,color='#97978f',lw=1,ls='--',zorder=0)
        if np.isfinite(q):
            ax.plot(q,level,'o',color='#ed733c',mec='white',ms=8,zorder=4)
            ax.vlines(q,0,level,color='#ed733c',lw=.8,ls=':',alpha=.6)
        labels.append(f'Top {alpha*100:.0f}%: {q:.2f}%')
    ax.set_title(f'{pop} | n = {len(group)}\n'+ '    '.join(labels),fontsize=16,loc='left',pad=10)
    if len(available)!=len(group):
        ax.text(.97,.04,f'{len(available)} available',transform=ax.transAxes,ha='right',fontsize=13)
    ax.set(xlim=(0,100),ylim=(0,102),xticks=range(0,101,20),yticks=range(0,101,20))
    ax.xaxis.set_major_formatter(PercentFormatter(100,decimals=0))
    ax.yaxis.set_major_formatter(PercentFormatter(100,decimals=0))
    ax.grid(alpha=.18)
    ax.spines[['top','right']].set_visible(False)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    paths=sorted(args.root.glob('*/regions/I50/null/rep*/calibration.json'))
    with ThreadPoolExecutor(max_workers=4) as pool:
        data=list(pool.map(read_record,paths))
    scores=pd.concat([d[0] for d in data],ignore_index=True)
    assert not scores.duplicated(['population','task_id','source','method','cutoff_years']).any()
    cuts=summaries(scores)
    scores.to_csv(args.out/'null_scores_10k_50k.csv',index=False)
    cuts.to_csv(args.out/'cutoffs_all_methods.csv',index=False)
    for source in ('decoded','truth'):
        for method in ('all_pairs','alt_alt'):
            subset=cuts[(cuts.source==source)&(cuts.method==method)].copy()
            subset['column']=subset.apply(lambda r:f'T{r.tmrca_years//1000}k_top{round(100*r.alpha)}_pct',axis=1)
            wide=subset.pivot(index='population',columns='column',values='empirical_percentile_pct').reindex(POPULATIONS)
            wide=wide.reindex(columns=['T10k_top5_pct','T10k_top1_pct','T50k_top5_pct','T50k_top1_pct'])
            wide.insert(0,'null_n',subset.groupby('population').null_n.first())
            wide.insert(1,'available_n',subset.groupby('population').available_n.first())
            wide.to_csv(args.out/f'four_cutoffs_{source}_{method}.csv')
    notes=('PROVISIONAL completed subset | step CDF; descriptive empirical percentiles\n'
           'Blue: focal population; gray: other populations. Percentiles are not finite-null p-value boundaries.')
    names=[]
    with PdfPages(args.out/'null_cdf_all.pdf') as combined:
        def save(fig,name):
            fig.savefig(args.out/(name+'.png'),dpi=150,bbox_inches='tight')
            fig.savefig(args.out/(name+'.pdf'),bbox_inches='tight')
            combined.savefig(fig,bbox_inches='tight')
            plt.close(fig)
            names.append(name)
        for source in ('decoded','truth'):
            for method in ('all_pairs','alt_alt'):
                title=f'{source.capitalize()} TMRCA | '+('All sampled pairs' if method=='all_pairs' else 'ALT/ALT pairs')
                for time in (10000,50000):
                    fig,axes=plt.subplots(2,3,figsize=(18,12),layout='constrained')
                    for ax,pop in zip(axes.flat,POPULATIONS): panel(ax,scores,cuts,pop,source,method,time)
                    fig.suptitle(title+f' | T < {time:,} years',fontsize=24)
                    fig.supxlabel(f'Pairs with TMRCA < {time:,} years (%)\n'+notes,fontsize=15)
                    fig.supylabel('Null focal sites at or below (%)',fontsize=19)
                    save(fig,f'cdf_{source}_{method}_T{time}')
                for pop in POPULATIONS:
                    fig,axes=plt.subplots(1,2,figsize=(14,8.5),layout='constrained')
                    for ax,time in zip(axes,(10000,50000)):
                        panel(ax,scores,cuts,pop,source,method,time)
                        ax.set_xlabel(f'Pairs with TMRCA < {time:,} years (%)')
                    axes[0].set_ylabel('Null focal sites at or below (%)')
                    fig.suptitle(title,fontsize=23)
                    fig.supxlabel(notes,fontsize=13)
                    save(fig,f'cdf_{pop}_{source}_{method}')
    explanation='''# Provisional neutral CDFs: completed simulations only

Each observation is one prespecified focal site from an independent retained
neutral region, not every position within the 10-Mb region. Focal alleles are
archaic-derived, segregating immediately after the 2% pulse at 50 kya, and
observed today; subsequent fixation is retained. No new simulations or decoding
were run. Both simulation campaigns remain paused.

The four main numbers are percentages of pairs with TMRCA below 10,000 or
50,000 years, not TMRCA ages. The blue step is the focal population ECDF; gray
steps are the other populations. Orange markers give descriptive 95th and 99th
percentiles using the inverse empirical CDF (nearest-rank order statistic).
They are placed on the nominal 95%/99% guides; a finite ECDF can jump past these
levels. No smoothing or fitted tail model is used. Missing ALT/ALT pair classes
are excluded from the displayed CDF and their counts are reported.

The descriptive percentiles are NOT the exact significance thresholds.
cutoffs_all_methods.csv separately reports conservative p=(1+#null>=score)/(n+1)
boundaries for p<=0.05 and p<=0.01; a finite target score must strictly exceed
the boundary. Missing null pair classes remain in that test's denominator,
matching the campaign convention. Fewer than 19 nulls cannot attain p<=0.05;
fewer than 99 cannot attain p<=0.01. Saturation at score=1 can also prevent
rejection. Unsupported thresholds are flagged, never replaced with percentiles.

All populations are below the planned 1,000 retained nulls. Small-sample tail
percentiles often equal the observed maximum. Completed trajectories may be
biased toward faster runs. These plots are provisional descriptions, not a
finished null calibration. The all-pairs decoded tables/figures are the primary
view; truth and ALT/ALT are supplied separately.

All source calibration receipts and truth/decoded score-file hashes were checked.
PNG and editable vector PDF outputs share the same figure code.
'''
    (args.out/'README.md').write_text(explanation)
    counts=scores.groupby('population').task_id.nunique().reindex(POPULATIONS)
    manifest=dict(schema='population-null-ecdf/v1',regions=int(counts.sum()),counts=counts.to_dict(),
        percentiles='numpy quantile method=inverted_cdf; available scores',provisional=True,
        source_receipts=[d[1] for d in data],script_sha256=digest(Path(__file__)),figure_count=len(names),
        simulations_started=False)
    files=[p for p in args.out.iterdir() if p.suffix in ('.pdf','.png','.csv','.md')]
    manifest['outputs']={p.name:dict(bytes=p.stat().st_size,sha256=digest(p)) for p in files}
    (args.out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    with zipfile.ZipFile(args.out/'null_cdf_bundle.zip','w',compression=zipfile.ZIP_DEFLATED) as archive:
        for p in files+[args.out/'manifest.json']: archive.write(p,p.name)
    print(json.dumps({k:v for k,v in manifest.items() if k not in ('source_receipts','outputs')}))


if __name__=='__main__': main()
