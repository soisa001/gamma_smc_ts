#!/usr/bin/env python3
"""Matched-null performance for the unconditional ongoing-selection pilot."""
import os
os.environ['MPLBACKEND']='Agg'
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.ndimage import label
from report_time_resolved_pilot import empirical_p, interval_values, expand_support, wilson

METHODS=('truth','median','mass')


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1<<20),b''):h.update(block)
    return h.hexdigest()


def evaluate(calls,positions,selected):
    """All 49 scanned support cells are inside the ongoing selection period."""
    mid=len(positions)//2;near=abs(positions-positions[mid])<=100000
    hit=bool(calls[near].any());n=int(calls[mid].sum())
    _,components=label(calls)
    tp=int(selected and hit)
    return dict(any_call=int(calls.any()),spatial_hit=int(hit),center_hit=int(n>0),
        called_kyr=n,temporal_f1=2*n/(49+n),temporal_recall=n/49,
        event_tp=tp,event_fp=int(components-tp),event_fn=int(selected)-tp)


def read(root):
    manifest=json.loads((root/'manifest.json').read_text());signature=sha(root/'manifest.json')
    receipts=sorted(root.glob('**/complete.json'))
    def unit(path):
        rec=json.loads(path.read_text());assert rec['config_sha256']==signature
        for key,name in [('scores_sha256','scores.npz'),('observed_sha256','observed.trees'),('posterior_sha256','posterior.zst'),('onset_sha256','onset.trees')]:
            assert sha(path.parent/name)==rec[key],f'Corrupt {path.parent/name}'
        trajectory=np.loadtxt(path.parent/'trajectory.tsv',skiprows=1)
        np.testing.assert_array_equal(trajectory[:,0],np.arange(50000,-1,-25))
        with np.load(path.parent/'scores.npz') as a:
            np.testing.assert_array_equal(a['positions'],np.arange(0,1000000,10000))
            np.testing.assert_array_equal(a['cutoffs'],np.arange(1000,50001,1000))
            np.testing.assert_array_equal(a['valid_pairs'],np.full(100,1000))
            cdf=a['cdfs'].copy()
            assert cdf.shape==(3,100,50) and np.isfinite(cdf).all()
            assert cdf.min()>=0 and cdf.max()<=1+1e-12 and np.diff(cdf,axis=-1).min()>=-1e-12
        return rec,cdf
    with ThreadPoolExecutor(max_workers=4) as pool:units=list(pool.map(unit,receipts))
    records=pd.DataFrame([u[0] for u in units]);assert records.seed.is_unique
    expected=json.loads((root/'tasks.json').read_text())
    completed=set(records.unit);missing=[t['unit'] for t in expected if t['unit'] not in completed]
    return manifest,records,np.stack([u[1] for u in units]),missing


def rate(indicator):
    n=len(indicator);k=int(np.sum(indicator));lo,hi=wilson(k,n)
    return dict(k=k,n=n,rate=k/n if n else None,ci95_low=lo,ci95_high=hi)


def report(root,out):
    out.mkdir(parents=True,exist_ok=True)
    manifest,records,cdfs,missing=read(root)
    positions=np.arange(0,1000000,10000)
    rows=[];summaries=[];calibration={};profiles=[]
    for origin in ('de_novo','introgression'):
        originmask=(records.origin==origin).to_numpy()
        development=originmask&(records.cohort=='development').to_numpy()
        cal=originmask&(records.cohort=='calibration').to_numpy()
        test=originmask&records.cohort.isin(['neutral_test','selected']).to_numpy()
        assert development.sum()>=10 and cal.sum()>=19,'Insufficient matched nulls'
        for scan,widths in [('one_kyr',[1]),('multiscale',[1,2,5,10]),('fixed_50kya',None)]:
            values,intervals=interval_values(cdfs,widths) if widths else (cdfs[...,-1:],None)
            for j,method in enumerate(METHODS):
                mu=values[development,j].mean(axis=(0,1));sd=values[development,j].std(axis=(0,1),ddof=1)
                valid=sd>0;z=np.full(values[:,j].shape,-np.inf)
                z[...,valid]=(values[:,j][...,valid]-mu[valid])/sd[valid]
                maxima=z[cal].max(axis=(1,2))
                calibration[f'{origin}_{scan}_{method}_maxima']=maxima
                calibration[f'{origin}_{scan}_{method}_mu']=mu
                calibration[f'{origin}_{scan}_{method}_sd']=sd
                for i in np.flatnonzero(test):
                    rec=records.iloc[i];called=empirical_p(z[i],maxima)<=.05
                    support=expand_support(called,intervals) if widths else np.repeat(called,49,axis=1)
                    metrics=evaluate(support,positions,rec.cohort=='selected')
                    if widths is None:
                        # A cumulative cutoff does not estimate an interval.
                        for name in ('called_kyr','temporal_f1','temporal_recall'):metrics[name]=None
                    rows.append(dict(unit=rec.unit,origin=origin,scan=scan,method=method,cohort=rec.cohort,
                        selection=rec.selection,final_frequency=rec.final_frequency,survived=rec.survived,
                        reached_10_percent=rec.reached_10_percent,global_p=float(empirical_p(z[i].max(),maxima)),**metrics))
                frame=pd.DataFrame([r for r in rows if r['origin']==origin and r['scan']==scan and r['method']==method])
                null=frame[frame.cohort=='neutral_test']
                for s,sel in frame[frame.cohort=='selected'].groupby('selection'):
                    joint=pd.concat([sel,null]);tp=int(sel.spatial_hit.sum());fp=int(null.any_call.sum());fn=len(sel)-tp
                    etp=int(joint.event_tp.sum());efp=int(joint.event_fp.sum());efn=int(joint.event_fn.sum())
                    summary=dict(origin=origin,scan=scan,method=method,selection=float(s),
                        region_f1=2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else None,
                        event_f1=2*etp/(2*etp+efp+efn) if 2*etp+efp+efn else None,
                        mean_temporal_f1=float(sel.temporal_f1.mean()) if widths else None,
                        mean_called_kyr=float(sel.called_kyr.mean()) if widths else None)
                    for metric,data in [('neutral_fpr',null.any_call),('any_call_power',sel.any_call),
                        ('spatial_tpr',sel.spatial_hit),('center_tpr',sel.center_hit),
                        ('survivor_tpr',sel[sel.survived].spatial_hit),('established_tpr',sel[sel.reached_10_percent].spatial_hit)]:
                        summaries.append(dict(summary,metric=metric,**rate(data)))
        for s,indices in records[originmask&(records.cohort=='selected').to_numpy()].groupby('selection').groups.items():
            for j,method in enumerate(METHODS):
                excess=np.diff(cdfs[list(indices),j,50].mean(axis=0)-cdfs[development,j,50].mean(axis=0))
                for k,v in enumerate(excess):profiles.append(dict(origin=origin,selection=s,method=method,time_kya=k+1.5,excess_mass=v))
    frame=pd.DataFrame(rows);summary=pd.DataFrame(summaries)
    frame.to_csv(out/'replicate_metrics.tsv',sep='\t',index=False,na_rep='NA')
    summary.to_csv(out/'performance.tsv',sep='\t',index=False,na_rep='NA')
    records.to_csv(out/'completed_units.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(profiles).to_csv(out/'center_excess_profiles.tsv',sep='\t',index=False)
    np.savez_compressed(out/'calibration.npz',**calibration)
    outcomes=[]
    for (origin,s),g in records[records.cohort=='selected'].groupby(['origin','selection']):
        outcomes.append(dict(origin=origin,selection=float(s),n=len(g),survived=int(g.survived.sum()),
            reached_10_percent=int(g.reached_10_percent.sum()),fixed=int((g.final_frequency==1).sum()),
            mean_final_frequency=float(g.final_frequency.mean()),median_final_frequency=float(g.final_frequency.median())))
    pd.DataFrame(outcomes).to_csv(out/'allele_outcomes.tsv',sep='\t',index=False)
    result=dict(manifest=manifest,root=str(root),missing_units=missing,report_source_sha256=sha(__file__),
        cohort_counts=records.groupby(['origin','cohort']).size().to_dict())
    result['cohort_counts']={':'.join(k):int(v) for k,v in result['cohort_counts'].items()}
    result['integrity']='All completed tree, onset, posterior and score hashes verified; all trajectories have 2,001 time points.'
    (out/'provenance.json').write_text(json.dumps(result,indent=2)+'\n')
    plot(out,summary)
    lines=['# Ongoing EAS selection pilot','',
        'Selection coefficient active from 50 kya to present; 25-year generations; Q=1; h=0.5. Each s is a homozygous fitness advantage.',
        'De novo: one copy at 50 kya. Introgression: nearest integer to 1% of 2N donor haplotypes, random pairing into EAS diploids; allele fixed in donors. Human/Neanderthal split 500 kya; Neanderthal Ne=3,600.',
        'EAS PHLASH median history; 1 Mb; 100 diploids; same 1,000 unfiltered random haplotype pairs; 10 kb stride; cutoffs 1–50 kya. This samples the all-pairs statistic, rather than enumerating every pair.',
        'Separate neutral cohorts for each origin. Interval scores use development-null means/SDs and calibration maxima across all positions and windows. Empirical upper-tail p includes ties and +1 correction; alpha=0.05.',
        'No retries or survival conditioning. Conditional power is descriptive and can have very small denominators. Each scan/method is calibrated separately; taking their union would need another correction.',
        'Spatial TPR requires a call within 100 kb of the known center. Neutral FPR means any call in a neutral region. Region F1 uses spatial TP, selected misses and any-call neutral FP, so it depends on the test mixture. Event F1 additionally penalizes unmatched connected components.',
        'All scanned times lie in the selection interval. Temporal overlap adds no timing specificity; outside-interval FPR is not estimable. Temporal F1 measures coverage of [1,50) kya and is not evidence of accurate onset/end dating. The 0–1 kya interval is outside this grid.',
        'Selection remains enabled after fixation, but there is no segregating allele left to respond. Loss/fixation times are retained. The old 30–20 kya comparison used standing variation near 5%, so this is not an isolated duration experiment.',
        '',f'Incomplete units: {len(missing)}.','',
        '| Origin | s | Scan | Posterior-mass spatial TPR | Neutral FPR | Survivor TPR | Region F1 |',
        '|---|---:|---|---:|---:|---:|---:|']
    mass=summary[summary.method=='mass']
    for (origin,s,scan),g in mass.groupby(['origin','selection','scan']):
        lookup=g.set_index('metric');fmt=lambda name:f"{int(lookup.loc[name,'k'])}/{int(lookup.loc[name,'n'])}"
        lines.append(f"| {origin} | {s:.2f} | {scan} | {fmt('spatial_tpr')} | {fmt('neutral_fpr')} | {fmt('survivor_tpr')} | {g.iloc[0].region_f1:.3f} |")
    lines.extend(['','Wilson 95% intervals and all truth/median/mass results are in performance.tsv. Raw simulation artifacts remain outside Git. PNG and editable vector PDF figures are saved to disk only.'])
    (out/'README.md').write_text('\n'.join(lines)+'\n')


def plot(out,summary):
    import matplotlib
    matplotlib.rcParams.update({'font.size':16,'axes.titlesize':17,'axes.labelsize':16,'pdf.fonttype':42,'ps.fonttype':42})
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(11,8.5),sharex=True,sharey=True)
    colors={'truth':'#333333','median':'#e69f00','mass':'#0072b2'}
    for i,origin in enumerate(('de_novo','introgression')):
        for j,scan in enumerate(('one_kyr','multiscale')):
            ax=axes[i,j]
            for method in METHODS:
                g=summary[(summary.origin==origin)&(summary.scan==scan)&(summary.method==method)&(summary.metric=='spatial_tpr')].sort_values('selection')
                x=g.selection.to_numpy()+(METHODS.index(method)-1)*.001
                ax.errorbar(x,g.rate,yerr=np.array([g.rate-g.ci95_low,g.ci95_high-g.rate]),fmt='o-',capsize=3,color=colors[method],label=method)
            ax.set(title=f"{origin.replace('_',' ')} · {'1 kyr' if j==0 else '1, 2, 5, 10 kyr'}",ylim=(-.03,1.03),xticks=[.01,.02,.03,.04,.05])
            ax.grid(axis='y',alpha=.2)
            if j==0:ax.set_ylabel('Spatial detection rate')
            if i==1:ax.set_xlabel('Selection coefficient s')
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,1),ncols=3,frameon=False)
    fig.tight_layout(rect=(0,0,1,.92))
    for ext in ('png','pdf'):fig.savefig(out/f'power.{ext}',dpi=160)
    plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',required=True,type=Path);p.add_argument('--out',required=True,type=Path)
    a=p.parse_args();report(a.root,a.out)
