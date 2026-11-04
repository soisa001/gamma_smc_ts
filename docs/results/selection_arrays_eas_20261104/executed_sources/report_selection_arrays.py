#!/usr/bin/env python3
"""Validate and compare three selection episodes using the entire temporal scan."""
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
from report_time_resolved_pilot import empirical_p,interval_values,expand_support
from report_ongoing_selection_pilot import sha,rate,METHODS

REGIMES=('10k','50k','40_20k')
LABELS={'10k':'10 kya → present','50k':'50 kya → present','40_20k':'40–20 kya'}


def local_path(value):
    """Resolve saved Windows/WSL paths on either side of the shared filesystem."""
    value=str(value).replace('\\','/')
    if os.name=='nt' and value.startswith('/mnt/') and len(value)>7 and value[6]=='/':
        value=value[5].upper()+':/'+value[7:]
    elif os.name!='nt' and len(value)>2 and value[1:3]==':/':
        value='/mnt/'+value[0].lower()+'/'+value[3:]
    return Path(value)


def evaluate(calls,positions,start,end,selected):
    mid=len(positions)//2;near=abs(positions-500000)<=100000
    lower=np.arange(1,50);truth=(lower>=end)&(lower<start)
    joint=bool(np.any(calls&near[:,None]&truth[None,:]))
    _,components=label(calls);tp=int(selected and joint)
    center=calls[mid];overlap=int(np.sum(center&truth));n=int(center.sum());total=int(truth.sum())
    labels,ncomponents=label(center)
    outside=sum(not np.any(truth[labels==i]) for i in range(1,ncomponents+1))
    near_labels,near_components=label(calls[near])
    near_outside=sum(not np.any((near_labels==i)&truth[None,:]) for i in range(1,near_components+1))
    changes=np.diff(np.r_[False,center,False].astype(int))
    intervals=list(zip((np.flatnonzero(changes==1)+1).tolist(),(np.flatnonzero(changes==-1)+1).tolist()))
    return dict(any_call=int(calls.any()),spatial_hit=int(calls[near].any()),joint_hit=int(joint),
        center_hit=int(overlap>0),center_mistimed=int(n>0 and overlap==0),center_off_target=int(outside>0),
        center_outside_components=outside,center_components=ncomponents,
        near_mistimed=int(calls[near].any() and not joint),near_off_target=int(near_outside>0),near_outside_components=near_outside,
        temporal_precision=overlap/n if n else None,temporal_recall=overlap/total,
        temporal_f1=2*overlap/(n+total),temporal_iou=overlap/(n+total-overlap),
        called_kyr=n,overlap_kyr=overlap,center_intervals_kya=json.dumps(intervals),
        event_tp=tp,event_fp=int(components-tp),event_fn=int(selected)-tp)


def read(root):
    manifest=json.loads((root/'manifest.json').read_text());contract=manifest['contract'];signature=sha(root/'manifest.json')
    schedule=json.loads((root/'schedule.json').read_text());tasks=json.loads((root/'tasks.json').read_text())
    def check(task):
        directory=local_path(task.get('source_directory',root/task['id']))
        receipt=next((directory/n for n in ('complete.json','rejected.json') if (directory/n).exists()),None)
        if receipt is None:return None
        r=json.loads(receipt.read_text());cfg=contract['configs'][task['regime']]
        if 'source_directory' in task:assert sha(receipt)==task['source_receipt_sha256']
        else:assert r['config_sha256']==signature
        if 'source_receipt' in r:assert sha(local_path(r['source_receipt']))==r['source_receipt_sha256']
        for filename,key in [('trajectory.tsv','trajectory_sha256'),('onset.trees','onset_sha256')]:assert sha(directory/filename)==r[key]
        tr=np.loadtxt(directory/'trajectory.tsv',skiprows=1,ndmin=2)
        age=r['ancestry_start_years'];np.testing.assert_array_equal(tr[:,0],age-np.arange(len(tr))*25)
        np.testing.assert_allclose(tr[:,1],tr[:,4]/(2*tr[:,2]),atol=1e-14)
        assert r['final_frequency']==tr[-1,1]
        if task['cohort']=='selected':
            active=(tr[:,0]<cfg['selection_start_years'])&(tr[:,0]>=cfg['selection_end_years'])
            # A marker introduced at onset already has s in its parental generation,
            # including introgression when the pulse and onset both occur at 50 kya.
            if age==cfg['selection_start_years']:active[0]=True
            alive=tr[:,4]>0
            np.testing.assert_allclose(tr[alive,3],np.where(active,task['selection'],0.)[alive],atol=1e-12)
            assert r['accepted']==(0<r['final_frequency']<1)
        else:assert r['accepted'] and np.all(tr[:,3]==0)
        cdf=None
        if r['accepted']:
            assert tr[-1,0]==0
            for filename,key in [('scores.npz','scores_sha256'),('observed.trees','observed_sha256'),('posterior.zst','posterior_sha256')]:assert sha(directory/filename)==r[key]
            with np.load(directory/'scores.npz') as a:
                np.testing.assert_array_equal(a['positions'],np.arange(0,1000000,10000))
                np.testing.assert_array_equal(a['cutoffs'],np.arange(1000,50001,1000))
                np.testing.assert_array_equal(a['valid_pairs'],np.full(100,1000))
                cdf=a['cdfs'].copy();assert cdf.shape==(3,100,50) and np.isfinite(cdf).all()
                assert cdf.min()>=0 and cdf.max()<=1+1e-12 and np.diff(cdf,axis=-1).min()>=-1e-12
        r.update(task,source_directory=str(directory),receipt_sha256=sha(receipt))
        return r,cdf
    with ThreadPoolExecutor(max_workers=4) as pool:
        checked=[]
        for i,result in enumerate(pool.map(check,tasks),1):
            checked.append(result)
            if i%500==0:print(f'Validated {i}/{len(tasks)} proposals',flush=True)
    completed={r[0]['id']:r[0] for r in checked if r is not None}
    arrays={r[0]['id']:r[1] for r in checked if r is not None and r[1] is not None}
    included=set();counts=[]
    for s in schedule['strata']:
        group=sorted([t for t in tasks if all(t[k]==s[k] for k in ('regime','origin','cohort','selection'))],key=lambda t:t['replicate'])
        chosen=[];accepted=0
        for t in group:
            if t['id'] not in completed or accepted>=s['target']:break
            r=completed[t['id']];chosen.append(r);included.add(t['id']);accepted+=r['accepted']
        counts.append(dict(**s,submitted=len(group),completed=sum(t['id'] in completed for t in group),analyzed_prefix=len(chosen),
            accepted=accepted,lost=sum(r['rejection']=='loss' for r in chosen),fixed=sum(r['rejection']=='fixation' for r in chosen),
            reused=sum('source_receipt_sha256' in r for r in chosen)))
    proposals=pd.DataFrame(completed.values());proposals['included_prefix']=proposals.id.isin(included)
    assert proposals.id.is_unique and proposals.seed.is_unique,'Duplicate proposal identity or seed'
    records=proposals[proposals.included_prefix&proposals.accepted].reset_index(drop=True)
    assert records.id.is_unique
    cdfs=np.stack([arrays[i] for i in records.id])
    return manifest,schedule,proposals,records,cdfs,pd.DataFrame(counts)


def report(root,out):
    out.mkdir(parents=True,exist_ok=True)
    manifest,schedule,proposals,records,cdfs,counts=read(root)
    rows=[];summaries=[];profiles=[];calibration={};unavailable=[]
    for scan,widths in [('one_kyr',[1]),('multiscale',[1,2,5,10])]:
        values,intervals=interval_values(cdfs,widths)
        for regime in REGIMES:
            cfg=manifest['contract']['configs'][regime]
            for origin in ('de_novo','introgression'):
                pool=schedule['null_pool_mapping'][regime][origin]
                mask=(records.origin==origin)&(records.regime==pool)
                dev=(mask&(records.cohort=='development')).to_numpy();cal=(mask&(records.cohort=='calibration')).to_numpy();test=(mask&(records.cohort=='neutral_test')).to_numpy()
                selected=((records.regime==regime)&(records.origin==origin)&(records.cohort=='selected')).to_numpy()
                if dev.sum()<5 or cal.sum()<19 or test.sum()<1:
                    unavailable.append(dict(regime=regime,origin=origin,development=int(dev.sum()),calibration=int(cal.sum()),test=int(test.sum())));continue
                for j,method in enumerate(METHODS):
                    mu=values[dev,j].mean(axis=(0,1));sd=values[dev,j].std(axis=(0,1),ddof=1);eligible=sd>0
                    z=np.full(values[:,j].shape,-np.inf);z[...,eligible]=(values[:,j][...,eligible]-mu[eligible])/sd[eligible]
                    maxima=z[cal].max(axis=(1,2));key=f'{scan}_{regime}_{origin}_{method}'
                    calibration[key+'_mu']=mu;calibration[key+'_sd']=sd;calibration[key+'_maxima']=maxima
                    tested=[]
                    for i in np.flatnonzero(test|selected):
                        called=expand_support(empirical_p(z[i],maxima)<=.05,intervals)
                        metrics=evaluate(called,np.arange(100)*10000,cfg['selection_start_years']/1000,cfg['selection_end_years']/1000,bool(selected[i]))
                        rec=records.iloc[i]
                        row=dict(scan=scan,regime=regime,origin=origin,method=method,id=rec.id,cohort=rec.cohort,selection=rec.selection,
                            final_frequency=rec.final_frequency,selection_onset_frequency=rec.selection_onset_frequency,
                            global_p=float(empirical_p(z[i].max(),maxima)),center_peak_kya=float(np.mean(intervals[np.argmax(z[i,50])])),**metrics)
                        rows.append(row);tested.append(row)
                    frame=pd.DataFrame(tested);null=frame[frame.cohort=='neutral_test']
                    for s in manifest['contract']['selection']:
                        sel=frame[(frame.cohort=='selected')&(frame.selection==s)]
                        tp=int(sel.joint_hit.sum());fp=int(null.any_call.sum());fn=len(sel)-tp
                        joined=pd.concat([sel,null]);etp=int(joined.event_tp.sum());efp=int(joined.event_fp.sum());efn=int(joined.event_fn.sum())
                        summary=dict(scan=scan,regime=regime,origin=origin,method=method,selection=s,n_selected=len(sel),n_development=int(dev.sum()),n_calibration=int(cal.sum()),n_null=len(null),
                            region_f1=2*tp/(2*tp+fp+fn) if len(sel) else None,
                            balanced_region_f1=2*(tp/len(sel))/(1+tp/len(sel)+fp/len(null)) if len(sel) else None,
                            event_f1=2*etp/(2*etp+efp+efn) if len(sel) else None,
                            mean_temporal_f1=float(sel.temporal_f1.mean()) if len(sel) else None,
                            mean_called_kyr=float(sel.called_kyr.mean()) if len(sel) else None)
                        for metric,data in [('joint_tpr',sel.joint_hit),('spatial_tpr',sel.spatial_hit),('neutral_fpr',null.any_call),
                            ('center_mistimed',sel.center_mistimed),('center_off_target',sel.center_off_target),('center_overlap',sel.center_hit),
                            ('near_mistimed',sel.near_mistimed),('near_off_target',sel.near_off_target)]:
                            summaries.append(dict(summary,metric=metric,**rate(data)))
                        if len(sel) and scan=='one_kyr':
                            indices=records.index[selected&(records.selection==s).to_numpy()]
                            excess=values[indices,j,50].mean(axis=0)-values[dev,j,50].mean(axis=0)
                            for k,v in enumerate(excess):profiles.append(dict(regime=regime,origin=origin,selection=s,method=method,time_kya=k+1.5,excess_mass=float(v)))
    summary=pd.DataFrame(summaries);frame=pd.DataFrame(rows)
    for name,data in [('performance',summary),('replicate_metrics',frame),('proposal_summary',counts),('proposals',proposals),('completed_units',records),('center_excess_profiles',pd.DataFrame(profiles))]:
        data.to_csv(out/f'{name}.tsv',sep='\t',index=False,na_rep='NA')
    records[records.cohort=='selected'].groupby(['regime','origin','selection']).agg(n=('id','size'),onset_af_median=('selection_onset_frequency','median'),
        final_af_min=('final_frequency','min'),final_af_median=('final_frequency','median'),final_af_max=('final_frequency','max')).reset_index().to_csv(out/'selected_frequencies.tsv',sep='\t',index=False)
    np.savez_compressed(out/'calibration.npz',**calibration)
    status=json.loads((root/'status.json').read_text())
    status_by_id={r['id']:r for r in status}
    pd.DataFrame([dict(t,**{k:v for k,v in status_by_id.get(t['id'],{'state':'incomplete_at_deadline'}).items() if k!='id'})
        for t in json.loads((root/'tasks.json').read_text())]).to_csv(out/'task_status.tsv',sep='\t',index=False,na_rep='NA')
    provenance=dict(manifest=manifest,schedule=schedule,root=str(root),unavailable=unavailable,report_sha256=sha(__file__),
        example_decoder_command=json.loads((local_path(records.iloc[0].source_directory)/'decoder.json').read_text())['command'],
        finished=json.loads((root/'finished.json').read_text()) if (root/'finished.json').exists() else None,
        finish_submitted=json.loads((root/'finish_submitted.json').read_text()) if (root/'finish_submitted.json').exists() else None,
        failed_units=[r for r in status if r['state']=='failed'],completed_proposals=len(proposals),analyzed_regions=len(records),
        analyzed_selected=int((records.cohort=='selected').sum()),
        integrity='Validated trajectory/onset/tree/posterior/score hashes, exact allele counts, onset and cessation, strict present-day segregation, grids and CDF bounds/monotonicity.',
        notes=['First 100 acceptances in uninterrupted proposal order per treatment; report actual deadline-limited counts.',
        'Null introgression pool shared across time windows because s=0 throughout; neutral de novo pools match forward duration.',
        'Primary: 49 adjacent intervals [1,2),...,[49,50) kya. Secondary: widths 1/2/5/10 kyr with 1-kyr starts, jointly calibrated over widths and positions. No data for [0,1) kya. Truth used only for evaluation.',
        'Conservative +1 empirical p-values from independent calibration-region maxima over all 100 positions and all candidate intervals; alpha .05 separately per estimator and scan. The union of the two scan variants is not calibrated.',
        'Fixed 1,000 unfiltered pairs out of 19,900; neither an exhaustive all-pairs posterior analysis nor independent pair tests.',
        'True-tree benchmarks apply the same interval statistic to true TMRCAs; they are not an optimal selection detector or a theoretical upper bound on power.',
        'Near-site mistiming is a call within 100 kb with zero temporal overlap; near-site off-target flags any wholly outside connected component, even alongside a correct component. Center versions restrict evaluation to the exact selected coordinate.',
        'Ongoing 50-kya selection occupies the full observed range, so outside-time false calls cannot be assessed there.',
        'center_peak_kya is the midpoint of the highest standardized window even when it is not significant; it is not by itself an inferred selection date.',
        'Region F1 depends on selected/null cohort sizes; balanced_region_f1 standardizes to equal selected and neutral prevalence. Temporal F1 compares center support to the observable selection interval.',
        'Reuse is restricted to prior 10-kya proposals and their nulls; previous simulation compute is outside the current budget.'])
    (out/'provenance.json').write_text(json.dumps(provenance,indent=2,allow_nan=False)+'\n')
    primary=summary[summary.scan=='one_kyr'] if len(summary) else summary
    lines=['# Three EAS selection episodes','',
        'Selection 10 kya to present, 50 kya to present, and 40–20 kya; s = 0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05; target 100 surviving segregating replicates per treatment (42 treatments).',
        'EAS PHLASH history, unscaled WF, h=0.5, 25 years/generation, 1 Mb, 100 diploids. De novo starts at one copy at onset. Introgression injects 1% donor haplotypes at 50 kya, with a fixed donor allele, human–Neanderthal split 500 kya and donor Ne 3,600. Selected loss/fixation rejected; neutral regions unconditioned with nearest segregating focal label.',
        'Detection scans all 49 adjacent one-kyr intervals from 1 to 50 kya, every 10 kb. Development nulls standardize interval coalescence fractions. Independent calibration maxima control the whole-region scan at nominal 5% per estimator; held-out nulls estimate FPR. All three estimators use the same 1,000 unfiltered pairs (of 19,900).',
        'Joint TPR requires a call within 100 kb of the selected site and overlapping the true selection interval. Spatial TPR ignores time. Temporal F1 is the Dice overlap at the exact focal coordinate. A wholly outside temporal call is reported separately from held-out neutral FPR. Region F1 depends on the selected/null mixture.',
        '',f'Completed proposals: {len(proposals)}; selected regions analyzed: {int((records.cohort=="selected").sum())}.','',
        '| Window | Origin | s | Accepted / target | Prefix proposals | Lost | Fixed |','|---|---|---:|---:|---:|---:|---:|']
    for r in counts[counts.cohort=='selected'].itertuples():lines.append(f'| {LABELS[r.regime]} | {r.origin} | {r.selection:g} | {r.accepted}/{r.target} | {r.analyzed_prefix} | {r.lost} | {r.fixed} |')
    lines+=['','| Window | Origin | s | Truth joint TPR | Median joint TPR | Mass joint TPR | Mass region F1 | Mass temporal F1 |','|---|---|---:|---:|---:|---:|---:|---:|']
    if len(summary):
        for (regime,origin,s),g in primary[primary.metric=='joint_tpr'].groupby(['regime','origin','selection'],sort=False):
            g=g.set_index('method');fmt=lambda m:f"{int(g.loc[m,'k'])}/{int(g.loc[m,'n'])}" if g.loc[m,'n'] else 'NA'
            val=lambda k:'NA' if pd.isna(g.loc['mass',k]) else f"{g.loc['mass',k]:.3f}"
            lines.append(f'| {LABELS[regime]} | {origin} | {s:g} | {fmt("truth")} | {fmt("median")} | {fmt("mass")} | {val("region_f1")} | {val("mean_temporal_f1")} |')
        lines+=['','| Window | Origin | Method | Neutral false-positive regions | 95% Wilson interval |','|---|---|---|---:|---:|']
        for r in primary[primary.metric=='neutral_fpr'].drop_duplicates(['regime','origin','method']).itertuples():lines.append(f'| {LABELS[r.regime]} | {r.origin} | {r.method} | {r.k}/{r.n} | {r.ci95_low:.3f}–{r.ci95_high:.3f} |')
        lines+=['','## Posterior mass: temporal errors near the selected site','',
            'Denominators are surviving selected regions. “Only outside” means a call within 100 kb but no temporal overlap. “Any outside component” also counts an outside component alongside a correct call. Outside-time error is unobservable for the 50-kya ongoing episode because it spans the entire scan.',
            '', '| Window | Origin | s | Only outside time | Any outside component |','|---|---|---:|---:|---:|']
        for (regime,origin,s),g in primary[(primary.method=='mass')&(primary.n>0)].groupby(['regime','origin','selection'],sort=False):
            g=g.set_index('metric')
            if 'near_mistimed' not in g.index:continue
            fmt=lambda metric:'Unobservable' if regime=='50k' else f"{int(g.loc[metric,'k'])}/{int(g.loc[metric,'n'])}"
            lines.append(f'| {LABELS[regime]} | {origin} | {s:g} | {fmt("near_mistimed")} | {fmt("near_off_target")} |')
        lines+=['','## Secondary scan: widths 1, 2, 5, 10 kyr','',
            'All windows start on the same 1-kyr grid. Each estimator is calibrated against maxima over all widths and positions; the primary/secondary union is not calibrated.',
            '', '| Window | Origin | s | Mass joint TPR | Mass temporal F1 | Neutral false-positive regions |','|---|---|---:|---:|---:|---:|']
        for (regime,origin,s),g in summary[(summary.scan=='multiscale')&(summary.method=='mass')].groupby(['regime','origin','selection'],sort=False):
            g=g.set_index('metric');n=int(g.loc['joint_tpr','n'])
            hit=f"{int(g.loc['joint_tpr','k'])}/{n}" if n else 'NA'
            f1=f"{g.loc['joint_tpr','mean_temporal_f1']:.3f}" if n else 'NA'
            null=f"{int(g.loc['neutral_fpr','k'])}/{int(g.loc['neutral_fpr','n'])}"
            lines.append(f'| {LABELS[regime]} | {origin} | {s:g} | {hit} | {f1} | {null} |')
    lines+=['','## Limits','']+['- '+s for s in provenance['notes']]
    if unavailable:lines+=['','Calibration unavailable for: '+json.dumps(unavailable)]
    (out/'README.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    if len(summary):
        plots(out)
        case_study(out,records,cdfs,frame,calibration)
    for path in out.iterdir():
        if path.suffix in ('.tsv','.json','.md'):path.write_bytes(path.read_bytes().replace(b'\r\n',b'\n'))
    (out/'output_hashes.json').write_text(json.dumps({p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name!='output_hashes.json'},indent=2)+'\n',encoding='utf-8',newline='\n')


def case_study(out,records,cdfs,frame,calibration):
    """Illustrative first qualifying case, not an independently chosen test set."""
    primary=frame[(frame.scan=='one_kyr')&(frame.regime=='40_20k')&(frame.cohort=='selected')]
    correct=set(primary[(primary.method=='truth')&(primary.joint_hit==1)].id)
    candidates=primary[(primary.method=='mass')&(primary.center_mistimed==1)&primary.id.isin(correct)]
    if candidates.empty:return
    chosen=candidates.iloc[0];i=records.index[records.id==chosen.id][0];record=records.loc[i]
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    fig,axes=plt.subplots(3,1,figsize=(11,8.5),sharex=True)
    times=np.arange(1.5,50,1);rows=[]
    for j,(ax,method) in enumerate(zip(axes,METHODS)):
        key=f'one_kyr_40_20k_{record.origin}_{method}';mu=calibration[key+'_mu'];sd=calibration[key+'_sd'];maxima=calibration[key+'_maxima']
        z=np.full(49,-np.inf);valid=sd>0;values=np.diff(cdfs[i,j,50]);z[valid]=(values[valid]-mu[valid])/sd[valid]
        called=empirical_p(z,maxima)<=.05
        allowed=int(np.floor(.05*(len(maxima)+1)-1+1e-9));threshold=np.sort(maxima)[len(maxima)-allowed-1]
        color={'truth':'#555555','median':'#0072B2','mass':'#D55E00'}[method]
        ax.axvspan(20,40,color='#dddddd',zorder=0);ax.axhline(threshold,color='black',linestyle='--',linewidth=1.5)
        ax.plot(times,np.where(valid,z,np.nan),color=color,linewidth=2)
        ax.scatter(times[called],z[called],color=color,s=24,zorder=3)
        ax.set_title({'truth':'True trees (1,000 pairs)','median':'Posterior median','mass':'Posterior mass'}[method],loc='left',pad=5)
        ax.grid(axis='y',alpha=.15);ax.set_xlim(1,50)
        rows.extend(dict(id=record.id,method=method,time_kya=t,z=float(v) if np.isfinite(v) else None,called=bool(c),threshold=float(threshold)) for t,v,c in zip(times,z,called))
    axes[-1].set_xticks([1,10,20,30,40,50]);axes[-1].set_xlabel('Coalescence time (kya)');fig.supylabel('Standardized interval excess')
    fig.suptitle(f'Example: {record.origin.replace("_"," ")}, s = {record.selection:g}, present AF = {record.final_frequency:.2f}',fontsize=17,y=.995)
    fig.legend(handles=[Patch(facecolor='#dddddd',label='True selection period'),Line2D([0],[0],color='black',linestyle='--',label='5% region threshold')],loc='upper center',ncol=2,bbox_to_anchor=(.54,.95),frameon=False,fontsize=14)
    fig.subplots_adjust(top=.84,bottom=.11,left=.12,right=.98,hspace=.5)
    for ext in ('png','pdf'):fig.savefig(out/f'older_episode_example.{ext}',dpi=180)
    plt.close(fig)
    pd.DataFrame(rows).to_csv(out/'older_episode_example.tsv',sep='\t',index=False,na_rep='NA')
    (out/'older_episode_example.json').write_text(json.dumps(dict(id=record.id,origin=record.origin,selection=float(record.selection),
        final_frequency=float(record.final_frequency),scores_sha256=record.scores_sha256,
        choice='First analyzed 40–20 kya region with true-tree joint detection and a posterior-mass center call wholly outside the selection period; illustrative, selected after examining outcomes.'),indent=2)+'\n')
    with (out/'README.md').open('a',encoding='utf-8') as f:
        f.write('\n## Illustrative dating error\n\n'
            f'`{record.id}` is the first analyzed 40–20 kya region with true-tree joint detection and a posterior-mass center call wholly outside the selected period. '
            'This example was selected after inspecting outcomes and is illustrative; aggregate performance above supplies the denominators. '
            'The three panels in `older_episode_example.pdf` use the same 1,000 pairs, with independently calibrated thresholds for each estimator.\n')


def plots(out):
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle,Patch
    matplotlib.rcParams.update({'font.size':15,'axes.titlesize':17,'axes.labelsize':16,'pdf.fonttype':42,'ps.fonttype':42})
    perf=pd.read_csv(out/'performance.tsv',sep='\t');counts=pd.read_csv(out/'proposal_summary.tsv',sep='\t')
    perf=perf[perf.scan=='one_kyr']
    ss=[.001,.005,.01,.02,.03,.04,.05];colors={'truth':'#555555','median':'#0072B2','mass':'#D55E00'}
    # Separate letter-size pages for each origin keep text and seven coefficients readable.
    for origin in ('de_novo','introgression'):
        fig,axes=plt.subplots(3,1,figsize=(11,8.5),sharex=True,sharey=True)
        for ax,regime in zip(axes,REGIMES):
            for j,method in enumerate(METHODS):
                g=perf[(perf.origin==origin)&(perf.regime==regime)&(perf.method==method)&(perf.metric=='joint_tpr')].set_index('selection').reindex(ss)
                x=np.arange(7)+(j-1)*.13;y=g['rate'].to_numpy()
                ax.errorbar(x,y,yerr=np.maximum(0,np.array([y-g.ci95_low,g.ci95_high-y])),fmt='o',capsize=3,color=colors[method],label={'truth':'True trees','median':'Posterior median','mass':'Posterior mass'}[method])
            ax.set_title(LABELS[regime],loc='left',pad=5);ax.set_ylim(-.07,1.07);ax.set_yticks([0,.5,1]);ax.grid(axis='y',alpha=.2)
            n=counts[(counts.regime==regime)&(counts.origin==origin)&(counts.cohort=='selected')].set_index('selection').reindex(ss).accepted
            ax.text(.99,1.04,'n = '+', '.join(str(int(v)) for v in n),transform=ax.transAxes,ha='right',va='bottom',fontsize=12)
        axes[-1].set_xticks(np.arange(7),[f'{s:g}' for s in ss]);axes[-1].set_xlabel('Selection coefficient s')
        fig.supylabel('Joint detection rate');handles,labels=axes[0].get_legend_handles_labels()
        fig.legend(handles,labels,loc='upper center',ncol=3,bbox_to_anchor=(.53,1),frameon=False)
        fig.subplots_adjust(top=.91,bottom=.1,left=.1,right=.98,hspace=.45)
        for ext in ('png','pdf'):fig.savefig(out/f'power_{origin}.{ext}',dpi=180)
        plt.close(fig)
    profiles=pd.read_csv(out/'center_excess_profiles.tsv',sep='\t')
    mass=profiles[profiles.method=='mass'];limit=max(abs(mass.excess_mass).max(),1e-4)
    cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#dddddd')
    for origin in ('de_novo','introgression'):
        fig,axes=plt.subplots(3,1,figsize=(11,8.5),sharex=True)
        for ax,regime in zip(axes,REGIMES):
            g=mass[(mass.origin==origin)&(mass.regime==regime)].pivot(index='selection',columns='time_kya',values='excess_mass').reindex(index=ss,columns=np.arange(1.5,50,1))
            im=ax.pcolormesh(np.arange(1,51),np.arange(8)-.5,g.to_numpy(),cmap=cmap,vmin=-limit,vmax=limit,rasterized=False,shading='flat')
            start,end={'10k':(10,1),'50k':(50,1),'40_20k':(40,20)}[regime]
            ax.add_patch(Rectangle((end,-.49),start-end,6.98,fill=False,edgecolor='black',linewidth=2))
            ax.set_yticks(np.arange(7),[f'{s:g}' for s in ss]);ax.tick_params(axis='y',labelsize=12);ax.set_title(LABELS[regime],loc='left',pad=5)
        axes[-1].set_xticks([1,10,20,30,40,50]);axes[-1].set_xlabel('Coalescence time (kya)');fig.supylabel('Selection coefficient s')
        fig.legend(handles=[Patch(facecolor='none',edgecolor='black',label='Simulated selection period'),
            Patch(facecolor='#dddddd',label='No surviving replicates')],loc='upper center',ncol=2,bbox_to_anchor=(.55,1),frameon=False,fontsize=13)
        fig.subplots_adjust(top=.90,bottom=.19,left=.13,right=.98,hspace=.45)
        cax=fig.add_axes([.26,.07,.6,.025]);cb=fig.colorbar(im,cax=cax,orientation='horizontal');cb.set_label('Posterior mass excess over neutral');cb.solids.set_rasterized(False)
        for ext in ('png','pdf'):fig.savefig(out/f'time_profiles_{origin}.{ext}',dpi=180)
        plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',required=True);p.add_argument('--out',required=True);p.add_argument('--plots-only',action='store_true');a=p.parse_args()
    if a.plots_only:plots(Path(a.out))
    else:report(Path(a.root),Path(a.out))
