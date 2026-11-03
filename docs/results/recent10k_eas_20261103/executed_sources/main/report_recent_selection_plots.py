#!/usr/bin/env python3
"""Independent whole-region calibration and temporal recovery for recent sweeps."""
import os
os.environ['MPLBACKEND']='Agg'
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.ndimage import label
from report_time_resolved_pilot import empirical_p, interval_values, expand_support
from report_ongoing_selection_pilot import sha, rate, METHODS


def evaluate(calls,positions,onset_kya=10,selected=True):
    mid=len(positions)//2;near=np.abs(positions-positions[mid])<=100000
    truth=np.arange(1,50)<onset_kya
    joint=bool(np.any(calls&near[:,None]&truth[None,:]))
    _,components=label(calls);tp=int(selected and joint)
    center=calls[mid];overlap=int(np.sum(center&truth));n=int(center.sum());total=int(truth.sum())
    center_labels,center_components=label(center)
    off=sum(not np.any(truth[center_labels==i]) for i in range(1,center_components+1))
    return dict(any_call=int(calls.any()),spatial_hit=int(calls[near].any()),joint_hit=int(joint),
        center_hit=int(overlap>0),center_mistimed=int(n>0 and overlap==0),center_off_target=int(off>0),
        temporal_precision=overlap/n if n else None,temporal_recall=overlap/total,
        temporal_f1=2*overlap/(n+total),temporal_iou=overlap/(n+total-overlap),
        called_kyr=n,overlap_kyr=overlap,event_tp=tp,event_fp=int(components-tp),event_fn=int(selected)-tp)


def select_prefix(tasks,completed,target):
    """First target accepted proposals, in proposal order, without skipping gaps."""
    used=[];naccepted=0
    for task in sorted(tasks,key=lambda t:t['replicate']):
        record=completed.get(task['unit'])
        if record is None:break
        used.append(task['unit']);naccepted+=record['accepted']
        if naccepted>=target:break
    return used


def read(root):
    manifest=json.loads((root/'manifest.json').read_text());cfg=manifest['config'];signature=sha(root/'manifest.json')
    schedule=json.loads((root/'schedule.json').read_text());tasks=json.loads((root/'tasks.json').read_text())
    def check(task):
        directory=root/task['unit'];receipt=directory/'complete.json'
        if not receipt.exists():receipt=directory/'rejected.json'
        if not receipt.exists():return None
        record=json.loads(receipt.read_text());assert record['config_sha256']==signature
        if 'source_receipt' in record:assert sha(Path(record['source_receipt']))==record['source_receipt_sha256']
        assert sha(directory/'trajectory.tsv')==record['trajectory_sha256']
        assert sha(directory/'onset.trees')==record['onset_sha256']
        tr=np.loadtxt(directory/'trajectory.tsv',skiprows=1,ndmin=2)
        age=record['ancestry_start_years'];np.testing.assert_array_equal(tr[:,0],age-np.arange(len(tr))*25)
        np.testing.assert_allclose(tr[:,1],tr[:,4]/(2*tr[:,2]),atol=1e-14)
        if task['origin']=='introgression':assert np.all(tr[tr[:,0]>=cfg['selection_start_years'],3]==0)
        if task['cohort']=='selected':
            alive=tr[:,4]>0;after=tr[:,0]<cfg['selection_start_years']
            np.testing.assert_allclose(tr[alive&after,3],task['selection'],atol=1e-12)
            assert record['accepted']==(0<record['final_frequency']<1)
        else:assert record['accepted'] and np.all(tr[:,3]==0)
        cdf=None
        if record['accepted']:
            assert tr[-1,0]==0
            for key,name in [('scores_sha256','scores.npz'),('posterior_sha256','posterior.zst'),('observed_sha256','observed.trees')]:assert sha(directory/name)==record[key]
            with np.load(directory/'scores.npz') as a:
                np.testing.assert_array_equal(a['positions'],np.arange(0,cfg['length'],cfg['stride']))
                np.testing.assert_array_equal(a['cutoffs'],np.arange(1000,50001,1000))
                np.testing.assert_array_equal(a['valid_pairs'],np.full(len(a['positions']),cfg['pairs']))
                cdf=a['cdfs'].copy();assert np.isfinite(cdf).all() and cdf.min()>=0 and cdf.max()<=1+1e-12
                assert np.diff(cdf,axis=-1).min()>=-1e-12
        record['receipt_sha256']=sha(receipt)
        return record,cdf
    with ThreadPoolExecutor(max_workers=4) as pool:checked=list(pool.map(check,tasks))
    completed={r[0]['unit']:r[0] for r in checked if r is not None}
    arrays={r[0]['unit']:r[1] for r in checked if r is not None and r[1] is not None}
    assert len({r['seed'] for r in completed.values()})==len(completed)
    used=set();counts=[]
    for stratum in schedule['strata']:
        group=[t for t in tasks if all(t[k]==stratum[k] for k in ('origin','cohort','selection'))]
        chosen=select_prefix(group,completed,stratum['target']);used.update(chosen)
        rows=[completed[u] for u in chosen]
        counts.append(dict(**stratum,submitted=len(group),completed=sum(t['unit'] in completed for t in group),
            analyzed_prefix=len(chosen),accepted=sum(r['accepted'] for r in rows),
            lost=sum(r['rejection']=='loss' for r in rows),fixed=sum(r['rejection']=='fixation' for r in rows)))
    proposals=pd.DataFrame(completed.values());proposals['included_prefix']=proposals.unit.isin(used)
    records=proposals[proposals.included_prefix&proposals.accepted].reset_index(drop=True)
    return manifest,schedule,proposals,records,np.stack([arrays[u] for u in records.unit]),pd.DataFrame(counts)


def report(root,out):
    out.mkdir(parents=True,exist_ok=True)
    manifest,schedule,proposals,records,cdfs,counts=read(root);cfg=manifest['config']
    positions=np.arange(0,cfg['length'],cfg['stride']);onset=cfg['selection_start_years']/1000
    rows=[];summaries=[];calibration={};profiles=[];unavailable=[]
    for origin in ('de_novo','introgression'):
        originmask=(records.origin==origin).to_numpy();dev=originmask&(records.cohort=='development').to_numpy()
        cal=originmask&(records.cohort=='calibration').to_numpy();test=originmask&(records.cohort=='neutral_test').to_numpy()
        selected=originmask&(records.cohort=='selected').to_numpy()
        if dev.sum()<5 or cal.sum()<19 or test.sum()<1:
            unavailable.append(dict(origin=origin,development=int(dev.sum()),calibration=int(cal.sum()),test=int(test.sum())));continue
        for scan,widths in [('one_kyr',[1]),('multiscale',[1,2,5,10]),('fixed_10kya',None),('fixed_50kya',None)]:
            if widths:values,intervals=interval_values(cdfs,widths)
            else:values,intervals=cdfs[...,[9 if scan=='fixed_10kya' else 49]],None
            for j,method in enumerate(METHODS):
                mu=values[dev,j].mean(axis=(0,1));sd=values[dev,j].std(axis=(0,1),ddof=1);eligible=sd>0
                z=np.full(values[:,j].shape,-np.inf);z[...,eligible]=(values[:,j][...,eligible]-mu[eligible])/sd[eligible]
                maxima=z[cal].max(axis=(1,2));prefix=f'{origin}_{scan}_{method}'
                calibration[prefix+'_mu']=mu;calibration[prefix+'_sd']=sd;calibration[prefix+'_maxima']=maxima
                tested=[]
                for i in np.flatnonzero(test|selected):
                    called=empirical_p(z[i],maxima)<=.05
                    support=expand_support(called,intervals) if widths else np.repeat(called,49,axis=1)
                    metrics=evaluate(support,positions,onset,bool(selected[i]))
                    mid=len(positions)//2
                    transitions=np.diff(np.r_[False,support[mid],False].astype(int))
                    starts=np.flatnonzero(transitions==1)+1;ends=np.flatnonzero(transitions==-1)+1
                    metrics['center_intervals_kya']=json.dumps(list(zip(starts.tolist(),ends.tolist()))) if widths else None
                    metrics['center_peak_kya']=float(np.argmax(z[i,mid])+1.5) if scan=='one_kyr' else None
                    if not widths:
                        for key in ('joint_hit','center_hit','center_mistimed','center_off_target','temporal_precision','temporal_recall','temporal_f1','temporal_iou','called_kyr','overlap_kyr'):metrics[key]=None
                    rec=records.iloc[i]
                    row=dict(origin=origin,scan=scan,method=method,unit=rec.unit,cohort=rec.cohort,selection=rec.selection,
                        final_frequency=rec.final_frequency,selection_onset_frequency=rec.selection_onset_frequency,
                        global_p=float(empirical_p(z[i].max(),maxima)),**metrics)
                    rows.append(row);tested.append(row)
                frame=pd.DataFrame(tested);null=frame[frame.cohort=='neutral_test']
                for s in schedule['parameters']['selection']:
                    sel=frame[(frame.cohort=='selected')&(frame.selection==s)]
                    tp=int(sel.joint_hit.sum()) if widths else int(sel.spatial_hit.sum());fp=int(null.any_call.sum());fn=len(sel)-tp
                    alltest=pd.concat([sel,null]);etp=int(alltest.event_tp.sum());efp=int(alltest.event_fp.sum());efn=int(alltest.event_fn.sum())
                    summary=dict(origin=origin,scan=scan,method=method,selection=s,
                        region_f1=2*tp/(2*tp+fp+fn) if len(sel) and 2*tp+fp+fn else None,
                        event_f1=2*etp/(2*etp+efp+efn) if len(sel) and 2*etp+efp+efn else None,
                        mean_temporal_f1=float(sel.temporal_f1.mean()) if widths and len(sel) else None,
                        mean_called_kyr=float(sel.called_kyr.mean()) if widths and len(sel) else None)
                    metrics=[('neutral_fpr',null.any_call),('any_call_power',sel.any_call),('spatial_tpr',sel.spatial_hit)]
                    if widths:metrics.extend([('joint_tpr',sel.joint_hit),('center_overlap',sel.center_hit),
                        ('center_mistimed',sel.center_mistimed),('center_off_target',sel.center_off_target)])
                    for metric,data in metrics:summaries.append(dict(summary,metric=metric,**rate(data)))
        for s,g in records[selected].groupby('selection'):
            for j,method in enumerate(METHODS):
                excess=np.diff(cdfs[g.index,j,len(positions)//2].mean(axis=0)-cdfs[dev,j,len(positions)//2].mean(axis=0))
                for k,value in enumerate(excess):profiles.append(dict(origin=origin,selection=s,method=method,time_kya=k+1.5,excess_mass=float(value)))
    summary=pd.DataFrame(summaries);frame=pd.DataFrame(rows)
    for name,data in [('performance',summary),('replicate_metrics',frame),('proposal_summary',counts),('proposals',proposals),('completed_units',records),('center_excess_profiles',pd.DataFrame(profiles))]:
        data.to_csv(out/f'{name}.tsv',sep='\t',index=False,na_rep='NA')
    records[records.cohort!='selected'][['unit','origin','cohort','neutral_focal_position_bp','neutral_focal_offset_bp','neutral_focal_minor_count']].to_csv(out/'neutral_focals.tsv',sep='\t',index=False)
    frequency=records[records.cohort=='selected'].groupby(['origin','selection']).agg(
        n=('unit','size'),onset_af_median=('selection_onset_frequency','median'),
        final_af_min=('final_frequency','min'),final_af_median=('final_frequency','median'),
        final_af_max=('final_frequency','max')).reset_index()
    frequency.to_csv(out/'selected_frequencies.tsv',sep='\t',index=False)
    np.savez_compressed(out/'calibration.npz',**calibration)
    status=json.loads((root/'status.json').read_text())
    result=dict(manifest=manifest,schedule=schedule,unavailable=unavailable,root=str(root),report_sha256=sha(__file__),
        example_decoder_command=json.loads((root/records.iloc[0].unit/'decoder.json').read_text())['command'],
        completed_proposals=len(proposals),analyzed_regions=len(records),analyzed_selected=int((records.cohort=='selected').sum()),
        failed_units=[r for r in status if r['state']=='failed'],
        integrity='Accepted tree/posterior/score hashes, all receipt trajectory/onset hashes, exact allele counts, activation date, pair/grid alignment and CDF monotonicity verified.',
        notes=['Fixed target is first 50 accepted replicates in proposal order; actual target is in schedule.',
            'Deadline-limited strata use a complete proposal prefix; Wilson intervals are exploratory.',
            'All methods and scan variants have separate null calibration; their union is not a calibrated discovery rule.',
            'Cumulative 10-kya comparator uses the known simulated onset and is an oracle comparator, not an onset inference method.',
            'Selection truth [0,10) kya is evaluated on observable [1,10) support: nine one-kyr cells; 0-1 kya is not scanned.',
            'Region F1 counts joint hits (spatial hits for cumulative baselines), selected misses, and neutral any-call false positives; it depends on cohort mixture.',
            'Event F1 also penalizes unmatched connected components; each selected region has at most one true event.'])
    (out/'provenance.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    lines=['# EAS selection beginning 10 kya','',
        'Selection coefficients 0.005, 0.01, 0.02, 0.03; h=0.5; 25 years/generation; unscaled WF forward simulation. Selected focal allele must be strictly segregating in the whole population at present.',
        'De novo: introduce one copy at 10 kya. Introgression: a 1% Neanderthal pulse at 50 kya, with allele fixed in donors, human–Neanderthal split 500 kya and Neanderthal Ne=3,600. The introgressed allele is neutral until 10 kya. Pre-onset losses/fixations are rejected because they cannot pass the final segregation condition.',
        'Neutral regions are unconditioned, with focal label at the nearest present-day sample-polymorphic variant. Nulls match the origin-specific demographic simulation and forward duration. EAS PHLASH median history; 1 Mb; 100 diploids; 1,000 unfiltered pairs; 10 kb stride; cutoffs 1–50 kya.',
        'Independent development and calibration nulls; conservative empirical max calibration over all positions and time windows at alpha=0.05. Primary intervals are one kyr. Secondary widths 1/2/5/10 kyr and cumulative 10/50-kya comparators are separately calibrated.',
        'Spatial TPR allows any time within 100 kb of the known center. Joint TPR additionally requires overlap with [1,10) kya. Temporal F1 is computed at the exact selected coordinate against those nine cells. Neutral FPR is any call anywhere in a held-out neutral region. Center mistiming means a center call with no temporal overlap; this differs from neutral FPR.',
        '',f'Completed proposals: {len(proposals)}; analyzed selected regions: {int((records.cohort=="selected").sum())}.','',
        '| Origin | s | Accepted | Proposals in prefix | Lost | Fixed |','|---|---:|---:|---:|---:|---:|']
    for r in counts[counts.cohort=='selected'].itertuples():lines.append(f'| {r.origin} | {r.selection:g} | {r.accepted} | {r.analyzed_prefix} | {r.lost} | {r.fixed} |')
    lines.extend(['','| Origin | s | Mass spatial TPR | Mass joint TPR | Neutral FPR | Temporal F1 |','|---|---:|---:|---:|---:|---:|'])
    for origin in ('de_novo','introgression'):
        for s in schedule['parameters']['selection']:
            g=summary[(summary.origin==origin)&(summary.selection==s)&(summary.scan=='one_kyr')&(summary.method=='mass')].set_index('metric')
            if g.empty:continue
            fmt=lambda metric:f"{int(g.loc[metric,'k'])}/{int(g.loc[metric,'n'])}" if g.loc[metric,'n'] else 'Not estimable'
            f1=g.loc['joint_tpr','mean_temporal_f1'];f1='NA' if pd.isna(f1) else f'{f1:.3f}'
            lines.append(f"| {origin} | {s:g} | {fmt('spatial_tpr')} | {fmt('joint_tpr')} | {fmt('neutral_fpr')} | {f1} |")
    lines.extend(['','## One-kyr estimator comparison','',
        '| Origin | s | True-tree joint TPR | Median joint TPR | Mass joint TPR | Mass region F1 | Mass event F1 |',
        '|---|---:|---:|---:|---:|---:|---:|'])
    for origin in ('de_novo','introgression'):
        for s in schedule['parameters']['selection']:
            g=summary[(summary.origin==origin)&(summary.selection==s)&(summary.scan=='one_kyr')&(summary.metric=='joint_tpr')].set_index('method')
            if g.empty:continue
            fmt=lambda method:f"{int(g.loc[method,'k'])}/{int(g.loc[method,'n'])}" if g.loc[method,'n'] else 'Not estimable'
            f1=lambda key:'NA' if pd.isna(g.loc['mass',key]) else f"{g.loc['mass',key]:.3f}"
            lines.append(f"| {origin} | {s:g} | {fmt('truth')} | {fmt('median')} | {fmt('mass')} | {f1('region_f1')} | {f1('event_f1')} |")
    lines.extend(['','## Held-out neutral error','',
        '| Origin | Estimator | False-positive regions | FPR | Wilson 95% interval |',
        '|---|---|---:|---:|---:|'])
    nulls=summary[(summary.scan=='one_kyr')&(summary.metric=='neutral_fpr')].drop_duplicates(['origin','method'])
    for r in nulls.itertuples():lines.append(f'| {r.origin} | {r.method} | {r.k}/{r.n} | {r.rate:.1%} | {r.ci95_low:.1%}–{r.ci95_high:.1%} |')
    lines.extend(['','Wilson intervals, region/event F1 and alternative scans are in performance.tsv. All simulation artifacts remain outside Git. Figures are saved only to PNG and editable vector PDF.']+['- '+n for n in result['notes']])
    (out/'README.md').write_text('\n'.join(lines)+'\n')
    plot(out,summary,pd.DataFrame(profiles),onset)


def plot(out,summary,profiles,onset):
    import matplotlib
    matplotlib.rcParams.update({'font.size':16,'axes.titlesize':17,'pdf.fonttype':42,'ps.fonttype':42})
    import matplotlib.pyplot as plt
    colors={'truth':'#444444','median':'#e69f00','mass':'#0072b2'}
    fig,axes=plt.subplots(2,2,figsize=(11,8.5),sharey=True)
    for i,origin in enumerate(('de_novo','introgression')):
        for j,metric in enumerate(('spatial_tpr','joint_tpr')):
            ax=axes[i,j]
            for method in METHODS:
                g=summary[(summary.origin==origin)&(summary.scan=='one_kyr')&(summary.method==method)&(summary.metric==metric)].sort_values('selection')
                ax.errorbar(np.arange(len(g))+(METHODS.index(method)-1)*.08,g.rate,
                    yerr=np.maximum(0,np.array([g.rate-g.ci95_low,g.ci95_high-g.rate])),fmt='o',capsize=3,color=colors[method],label=method)
            denominators=summary[(summary.origin==origin)&(summary.scan=='one_kyr')&(summary.method=='mass')&(summary.metric==metric)].set_index('selection').n
            ticks=[.005,.01,.02,.03]
            ax.set(title=f"{origin.replace('_',' ')} · {'location' if j==0 else 'location + time'}",ylim=(-.03,1.03),xticks=np.arange(len(ticks)),xlim=(-.4,3.4),
                xticklabels=[f"{s:g}\n(n={int(denominators.get(s,0))})" for s in ticks])
            ax.grid(axis='y',alpha=.2)
            if j==0:ax.set_ylabel('Detection rate')
            if i==1:ax.set_xlabel('Selection coefficient s')
    handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,1),ncols=3,frameon=False)
    fig.tight_layout(rect=(0,0,1,.92))
    for ext in ('png','pdf'):fig.savefig(out/f'performance.{ext}',dpi=160)
    plt.close(fig)
    fig,axes=plt.subplots(2,1,figsize=(11,8.5),sharex=True)
    selection_colors=dict(zip(sorted(profiles.selection.unique()),plt.get_cmap('tab10').colors))
    for ax,origin in zip(axes,('de_novo','introgression')):
        ax.axvspan(1,onset,color='#cccccc',alpha=.35);ax.axhline(0,color='#777777',lw=1)
        for s,g in profiles[(profiles.origin==origin)&(profiles.method=='mass')].groupby('selection'):
            ax.plot(g.time_kya,g.excess_mass,label=f's={s:g}',lw=2,color=selection_colors[s])
        ax.set(title=origin.replace('_',' '),ylabel='Excess posterior mass')
    axes[-1].set_xlabel('Coalescence time (kya)');handles,labels=axes[-1].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,1),ncols=4,frameon=False)
    fig.tight_layout(rect=(0,0,1,.92))
    for ext in ('png','pdf'):fig.savefig(out/f'center_excess.{ext}',dpi=160)
    plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--plots-only',action='store_true',help='Render saved metric tables without recalculating simulations or statistics')
    a=p.parse_args()
    if a.plots_only:
        provenance=json.loads((a.out/'provenance.json').read_text())
        assert Path(provenance['root']).resolve()==a.root.resolve()
        plot(a.out,pd.read_csv(a.out/'performance.tsv',sep='\t'),pd.read_csv(a.out/'center_excess_profiles.tsv',sep='\t'),
             provenance['manifest']['config']['selection_start_years']/1000)
    else:report(a.root,a.out)
    (a.out/'plot_provenance.json').write_text(json.dumps(dict(source_sha256=sha(__file__),plots_only=a.plots_only,
        files={name:sha(a.out/name) for name in ('performance.tsv','center_excess_profiles.tsv','performance.png','performance.pdf','center_excess.png','center_excess.pdf')}),indent=2)+'\n')
