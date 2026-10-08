"""Read verified completed calibration receipts; export ECDFs without running sims."""
import os
os.environ['MPLBACKEND'] = 'Agg'
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[name] = '1'
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
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
from gamma_smc_aou.population_neutral import POPULATIONS


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_record(path):
    record = json.loads(path.read_text())
    pop = record['population']
    task = record['task']
    assert task['role']=='null' and task['s']==0 and task['ascertainment_years']==50000
    root = path.parents[5]
    assert path == root/pop/'regions'/task['id']/'calibration.json'
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
    parameters=record['identity']['simulation']['parameters']
    return pd.concat(frames),dict(path=str(path),sha256=digest(path),population=pop,task_id=task['id'],
        sample_af=record['sample_af'],onset_af=record['onset_af'],attempts=record['attempts'],
        history_id=parameters.get('demographic_draw',{}).get('history_id'))


def ecdf(values):
    x=np.sort(np.asarray(values,dtype=float))
    n=len(x)
    assert n and not np.isinf(x).any()
    x=x[np.isfinite(x)]
    missing=n-len(x)
    return np.r_[0,x,1]*100,np.r_[missing/n,(missing+np.arange(1,len(x)+1))/n,1]*100


def empirical_p(values, target):
    """Inclusive empirical tail; null no-calls stay in n, target no-call is NaN."""
    values=np.asarray(values,dtype=float)
    assert len(values) and not np.isinf(values).any()
    if not np.isfinite(target):
        return np.nan
    return float(np.count_nonzero(values>=target)/len(values))


def empirical_summary(values, alpha):
    values=np.asarray(values,dtype=float)
    assert len(values) and not np.isinf(values).any() and 0<alpha<1
    finite=values[np.isfinite(values)]
    reference=np.sort(np.where(np.isnan(values),-np.inf,values))
    q=float(reference[math.ceil((1-alpha)*len(values))-1])
    allowed=int(np.flatnonzero(np.arange(len(values)+1)/len(values)<=alpha)[-1])
    boundary=float(reference[-allowed-1])
    return dict(empirical_percentile_pct=q*100,empirical_p_boundary_pct=boundary*100,
        available_only_percentile_pct=float(np.quantile(finite,1-alpha,method='inverted_cdf'))*100 if len(finite) else np.nan,
        p_at_cutoff=empirical_p(values,q),tied_n=int(np.count_nonzero(values==q)),
        tied_fraction=float(np.count_nonzero(values==q)/len(values)),
        p_at_score_one=empirical_p(values,1.),p_grid_step=1/len(values),
        empirical_p_attainable=empirical_p(values,1.)<=alpha)


def summaries(scores):
    rows=[]
    for key,g in scores.groupby(['population','source','method','cutoff_years']):
        available=g.score.dropna().to_numpy()
        for alpha in (.05,.01):
            summary=empirical_summary(g.score,alpha)
            rows.append(dict(zip(['population','source','method','tmrca_years'],key),
                alpha=alpha,null_n=len(g),available_n=len(available),missing_n=len(g)-len(available),
                **summary,p_definition='count(null >= observed) / n; no add-one adjustment',
                decision='finite score strictly greater than boundary; missing target is no call',provisional=True))
    return pd.DataFrame(rows)


def panel(ax, scores, cuts, pop, source, method, time):
    selected=scores[(scores.source==source)&(scores.method==method)&(scores.cutoff_years==time)]
    for other in POPULATIONS:
        values=selected[selected.population==other].score
        if len(values):
            ax.step(*ecdf(values),where='post',color='#c5c4ba',lw=1.2,zorder=1)
    group=selected[selected.population==pop]
    available=group.score.dropna()
    if len(group): ax.step(*ecdf(group.score),where='post',color='#2074c9',lw=2.8,zorder=3)
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
        labels.append(f'{100*(1-alpha):.0f}th: {q:.2f}%')
    ax.set_title(f'{pop} | n = {len(group)}\n'+ '    '.join(labels),fontsize=16,loc='left',pad=10)
    if len(available)!=len(group):
        ax.text(.97,.04,f'{len(available)} available',transform=ax.transAxes,ha='right',fontsize=13)
    ax.set(xlim=(0,100),ylim=(0,102),xticks=range(0,101,20),yticks=range(0,101,20))
    ax.xaxis.set_major_formatter(PercentFormatter(100,decimals=0))
    if method=='all_pairs' and time==10000:
        upper=max(1.,float(np.ceil(115*selected.score.max())))
        ax.set(xlim=(0,upper),xticks=np.linspace(0,upper,5))
        ax.xaxis.set_major_formatter(PercentFormatter(100,decimals=1 if upper % 4 else 0))
    ax.yaxis.set_major_formatter(PercentFormatter(100,decimals=0))
    ax.grid(alpha=.18)
    ax.spines[['top','right']].set_visible(False)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--reuse-snapshot',action='store_true',help='Redraw the existing frozen receipt list, not newly completed simulations')
    parser.add_argument('--first-per-population',type=int,help='Require and use exactly replicates 0 through N-1 in every population')
    args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    if args.reuse_snapshot:
        old=json.loads((args.out/'manifest.json').read_text())
        paths=[Path(r['path']) for r in old['source_receipts']]
        for path,record in zip(paths,old['source_receipts']):
            assert digest(path)==record['sha256'],f'Changed frozen receipt: {path}'
    else:
        paths=sorted(args.root.glob('*/regions/I50/null/rep*/calibration.json'))
    if args.first_per_population is not None:
        assert args.first_per_population>0
        expected={args.root/pop/'regions/I50/null'/f'rep{rep:04d}'/'calibration.json'
                  for pop in POPULATIONS for rep in range(args.first_per_population)}
        assert expected.issubset(set(paths)),f'Missing required receipts: {sorted(expected-set(paths))}'
        if args.reuse_snapshot:
            assert set(paths)==expected,'Frozen snapshot differs from requested cohort'
        paths=sorted(expected)
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
    notes=('Empirical p = count(null >= observed) / n; ties count; no +1 adjustment\n'
           'Blue: focal population; gray: others. Null no-calls remain in n, below finite support.')
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
    explanation='''# Empirical neutral CDFs: frozen completed cohort

Each observation is one prespecified focal site from an independent retained
neutral region, not every position within the 10-Mb region. Focal alleles are
archaic-derived, segregating immediately after the 2% pulse at 50 kya, and
observed today; subsequent fixation is retained. No new simulations or decoding
were run by this reporting script. The simulation runner is not modified.

The four main numbers are percentages of pairs with TMRCA below 10,000 or
50,000 years, not TMRCA ages. The blue step is the focal population ECDF; gray
steps are the other populations. Orange markers give empirical 95th and 99th
percentiles using the inverse empirical CDF (nearest-rank order statistic).
They are placed on the nominal 95%/99% guides; a finite ECDF can jump past these
levels. No smoothing or fitted tail model is used. Missing ALT/ALT pair classes
remain in the denominator below finite score support; their CDF mass appears
at the left edge. These are no-calls, not biological zero scores.

The user's specified p-value is count(null >= observed) / n, without any +1
adjustment. Ties count. Undefined targets are no-calls. With n=100, p<=0.05
requires strictly exceeding sorted score 95; p<=0.01 requires strictly
exceeding sorted score 99. Equality at those boundaries does not reject.
The primary percentiles use all n nulls, with null no-calls below finite
support. Available-only percentiles are separately labeled in the CSV.
The CSV also exports tied counts/fractions, p at the cutoff and p at score 1.
When all null scores equal the observation, p=1. At score 1, p is the fraction
of null scores equal to 1. Saturation can prevent rejection at either level.

All populations are below the planned 1,000 retained nulls. With 100 nulls,
the p grid has 0.01 steps and the 99th percentile is the second-largest score;
the 1% tail is imprecisely estimated. A zero empirical tail is possible and
does not establish zero population probability. These plots are provisional
calibration results. The all-pairs decoded tables/figures are the primary
view; truth and ALT/ALT are supplied separately.

All source calibration receipts and truth/decoded score-file hashes were checked.
PNG and editable vector PDF outputs share the same figure code.
'''
    (args.out/'README.md').write_text(explanation)
    counts=scores.groupby('population').task_id.nunique().reindex(POPULATIONS)
    history_rows=[]
    for pop in POPULATIONS:
        histories=[d[1]['history_id'] for d in data if d[1]['population']==pop and d[1]['history_id'] is not None]
        sizes=pd.Series(histories,dtype=int).value_counts()
        history_rows.append(dict(population=pop,completed_regions=int(counts[pop]),distinct_histories=len(sizes),
            minimum_regions_per_history=int(sizes.min()) if len(sizes) else 0,
            maximum_regions_per_history=int(sizes.max()) if len(sizes) else 0,
            complete_ten_region_histories=int((sizes==10).sum())))
    pd.DataFrame(history_rows).to_csv(args.out/'history_coverage.csv',index=False)
    if any(row['distinct_histories'] for row in history_rows):
        explanation+='''
This snapshot uses the mapped-autosomal-region MVN-demography campaign. The
CDFs pool completed regions. Demographic-history CDF bands are not inferred
from singleton histories; history_coverage.csv reports replication. Completion
order can favor faster histories, so incomplete-cohort tails remain provisional.
'''
        (args.out/'README.md').write_text(explanation)
    manifest=dict(schema='population-null-ecdf/v2',regions=int(counts.sum()),counts=counts.to_dict(),
        first_per_population=args.first_per_population,
        p_definition='count(null >= observed) / n; no add-one adjustment',
        percentiles='nearest rank; all nulls; missing nulls below finite support',provisional=True,
        source_receipts=[d[1] for d in data],script_sha256=digest(Path(__file__)),figure_count=len(names),
        simulations_started=False)
    base_names={f'{name}{suffix}' for name in names for suffix in ('.pdf','.png')}
    base_names.update(('null_cdf_all.pdf','null_scores_10k_50k.csv','cutoffs_all_methods.csv','history_coverage.csv','README.md'))
    base_names.update(f'four_cutoffs_{source}_{method}.csv' for source in ('decoded','truth') for method in ('all_pairs','alt_alt'))
    files=sorted(p for p in args.out.iterdir() if p.name in base_names)
    manifest['outputs']={p.name:dict(bytes=p.stat().st_size,sha256=digest(p)) for p in files}
    (args.out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    with zipfile.ZipFile(args.out/'null_cdf_bundle.zip','w',compression=zipfile.ZIP_DEFLATED) as archive:
        for p in files+[args.out/'manifest.json']: archive.write(p,p.name)
    print(json.dumps({k:v for k,v in manifest.items() if k not in ('source_receipts','outputs')}))


if __name__=='__main__': main()
