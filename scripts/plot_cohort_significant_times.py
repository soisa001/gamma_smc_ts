#!/usr/bin/env python3
"""Show each selected simulation's pointwise calls against a 100-neutral bank.

prepare verifies and caches the selected-site scores. --available-only permits
caching finished cohorts while the 50-30-kya simulations run. report needs only
the portable cache and renders vector PDFs and PNGs at p<.05 and p<.01.
The report function also accepts explicitly labelled partial selected cohorts.
"""
import os
os.environ['MPLBACKEND']='Agg'
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import re
import time
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from plot_pointwise_significant_times import runs
from plot_selection_episode_series import sha,save_json
from report_selection_arrays import local_path
from pointwise_time_tests import exceedance_counts

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'docs/results/selection_arrays_eas_20261105_4h'
NEUTRAL=ROOT/'docs/results/selected_site_posteriors_20261106'
RAW=ROOT.parent/'sim/selection_50_30k_demo_20261106'
CASES=[dict(key='50_30k_s001',regime='50_30k',selection=.01,window=[30,50],label='50–30 kya · s=0.01'),
       dict(key='10k_s001',regime='10k',selection=.01,window=[0,10],label='10 kya–present · s=0.01'),
       dict(key='50k_s001',regime='50k',selection=.01,window=[0,50],label='50 kya–present · s=0.01'),
       dict(key='50k_s0005',regime='50k',selection=.005,window=[0,50],label='50 kya–present · s=0.005')]
METHODS=['True TMRCA','Posterior median','Posterior-weighted mass']
METHOD_KEYS=['true_tmrca','posterior_median','posterior_mass']
ROW_ORDER='Posterior mass: earliest significant time ascending, summed significant spans across panels descending, significant-test count descending, proposal index ascending; based on p<0.05 and shared across panels/thresholds. No-call simulations last. Other methods retain proposal order.'


def mass_order(data):
    masks=[data['p_cumulative'][:,2]<.05,data['p_interval'][:,2]<.05]
    def key(i):
        hits=[np.flatnonzero(mask[i]) for mask in masks]
        first=min((h[0] for h in hits if len(h)),default=1000)
        span=sum(int(h[-1]-h[0]+1) for h in hits if len(h))
        count=sum(len(h) for h in hits)
        return first,-span,-count,int(data['proposal_indices'][i])
    order=np.array(sorted(range(len(masks[0])),key=key))
    np.testing.assert_array_equal(np.sort(order),np.arange(len(masks[0])))
    return order


def prepare(out,workers,available_only):
    inventory=json.loads((BASE/'output_hashes.json').read_text())
    assert sha(BASE/'completed_units.tsv')==inventory['completed_units.tsv']
    base_configs=json.loads((BASE/'manifest.json').read_text())['contract']['configs']
    episode_config=json.loads((RAW/'manifest.json').read_text())['contract']['config']
    all_rows=pd.read_csv(BASE/'completed_units.tsv',sep='\t')
    neutral_meta=json.loads((NEUTRAL/'manifest.json').read_text())
    assert sha(NEUTRAL/'parameters.npz')==neutral_meta['cache_sha256']
    with np.load(NEUTRAL/'parameters.npz') as data:
        neutral=data['cdfs'][3:].copy();pairs=data['pairs'].copy();neutral_ids=data['ids'][3:].copy()
    assert neutral.shape==(100,3,50) and len(set(neutral_ids))==100
    if (out/'neutral.npz').exists():
        with np.load(out/'neutral.npz') as data:
            np.testing.assert_array_equal(data['cdfs'],neutral)
            np.testing.assert_array_equal(data['ids'],neutral_ids)
            np.testing.assert_array_equal(data['pairs'],pairs)
    else:np.savez_compressed(out/'neutral.npz',cdfs=neutral,ids=neutral_ids,pairs=pairs)
    completed,pending=[],[]
    for case in CASES:
        if case['regime']=='50_30k':
            cohort=RAW/'cohort_00100.json'
            if not cohort.exists():
                pending.append(case['key'])
                if available_only:
                    print('Pending 50–30-kya cohort; caching other treatments.',flush=True)
                    continue
                raise FileNotFoundError(f'Await completed 100-survivor cohort: {cohort}')
            receipt=json.loads(cohort.read_text())
            assert receipt['target']==100 and receipt['start_replicate']==0
            rows=[dict(r,source_directory=str(RAW/r['id'])) for r in receipt['selected_cohort']]
            source=dict(path=str(cohort),sha256=sha(cohort))
        else:
            subset=all_rows[(all_rows.regime==case['regime'])&(all_rows.origin=='introgression')&
                            (all_rows.cohort=='selected')&(all_rows.selection==case['selection'])].sort_values('replicate')
            assert len(subset)>=100
            rows=subset.iloc[:100].to_dict('records')
            source=dict(path=str(BASE/'completed_units.tsv'),sha256=inventory['completed_units.tsv'])
        assert len(rows)==100 and len({r['id'] for r in rows})==100
        assert len({r['observed_sha256'] for r in rows})==100, 'Duplicate selected simulation trees.'
        assert [r['replicate'] for r in rows]==sorted(r['replicate'] for r in rows)
        assert all(r['accepted'] and 0<r['final_frequency']<1 for r in rows)
        assert all(r['origin']=='introgression' and r['cohort']=='selected' and
                   r['regime']==case['regime'] and r['selection']==case['selection'] for r in rows)
        contract=dict(case=case,position_bp=500000,ids=[r['id'] for r in rows],
            score_hashes=[r['scores_sha256'] for r in rows],trajectory_hashes=[r['trajectory_sha256'] for r in rows],
            observed_hashes=[r['observed_sha256'] for r in rows],posterior_hashes=[r['posterior_sha256'] for r in rows],
            neutral_cache_sha256=sha(out/'neutral.npz'),choice='First 100 accepted simulations in proposal order; no ascertainment by inferred signal.',
            conditioning='Selected: 0 < whole-population focal AF < 1 today. Neutral: unconditioned.')
        cache=out/f"{case['key']}.npz";manifest=out/f"{case['key']}_manifest.json"
        if manifest.exists():
            saved=json.loads(manifest.read_text())
            assert saved['contract']==contract and sha(cache)==saved['cache_sha256']
            if saved.get('p_value_version',1)<2:
                with np.load(cache) as previous:content={k:previous[k].copy() for k in previous.files}
                count_c,count_i=exceedance_counts(content['cdfs'],neutral)
                content.update(p_cumulative=(count_c+1)/101,p_interval=(count_i+1)/101,
                    neutral_ge_cumulative=count_c,neutral_ge_interval=count_i)
                np.savez_compressed(cache,**content)
                saved.update(cache_sha256=sha(cache),p_value_version=2,
                    comparison='Exact integer pair counts for truth/median; original continuous values for posterior mass.')
                save_json(manifest,saved)
            completed.append(dict(case=case,cache=cache.name,manifest=manifest.name,cache_sha256=sha(cache)))
            print(f"Verified existing cache: {case['key']}",flush=True)
            continue
        def read(row):
            directory=local_path(row['source_directory'])
            for name,key in [('scores.npz','scores_sha256'),('trajectory.tsv','trajectory_sha256'),
                             ('observed.trees','observed_sha256'),('posterior.zst','posterior_sha256')]:
                assert sha(directory/name)==row[key],str(directory/name)
            meta=json.loads((directory/'posterior.zst.meta').read_text())
            np.testing.assert_array_equal(meta['pairs'],pairs)
            np.testing.assert_array_equal(meta['output_positions'],np.arange(100)*10000)
            with np.load(directory/'scores.npz') as data:
                np.testing.assert_array_equal(data['positions'],np.arange(100)*10000)
                np.testing.assert_array_equal(data['cutoffs'],np.arange(1000,50001,1000))
                np.testing.assert_array_equal(data['valid_pairs'],np.full(100,1000))
                cdfs=data['cdfs'][:,50,:].copy()
            assert cdfs.shape==(3,50) and np.isfinite(cdfs).all()
            assert cdfs.min()>=0 and cdfs.max()<=1+1e-12 and np.diff(cdfs,axis=-1).min()>=-1e-12
            trajectory=np.loadtxt(directory/'trajectory.tsv',skiprows=1)
            np.testing.assert_array_equal(trajectory[:,0],np.arange(50000,-1,-25))
            ticks=row['selected_reproduction_ticks']
            if isinstance(ticks,str):ticks=json.loads(ticks)
            # Tick 1 is the initial state; when selection starts at introduction
            # it already has s. Later onset activates at its first reproduction tick.
            first_tick=1 if ticks[0]==2 else ticks[0]
            generation_ticks=np.arange(1,2002)
            expected=case['selection']*((generation_ticks>=first_tick)&(generation_ticks<=ticks[1]))
            np.testing.assert_allclose(trajectory[:,3],expected,atol=1e-12)
            np.testing.assert_allclose(trajectory[:,1],trajectory[:,4]/(2*trajectory[:,2]),atol=1e-14)
            assert 0<trajectory[-1,1]<1
            np.testing.assert_allclose(trajectory[-1,1],row['final_frequency'],atol=1e-14)
            return cdfs,trajectory[:,1],dict(id=row['id'],replicate=int(row['replicate']),seed=int(row['seed']),
                directory=str(directory),metadata_sha256=sha(directory/'posterior.zst.meta'))
        with ThreadPoolExecutor(max_workers=workers) as pool:results=list(pool.map(read,rows))
        cdfs=np.stack([r[0] for r in results])
        count_c,count_i=exceedance_counts(cdfs,neutral)
        p_c=(count_c+1)/101;p_i=(count_i+1)/101
        for alpha,max_count in [(.05,4),(.01,0)]:
            np.testing.assert_array_equal(p_c<alpha,count_c<=max_count)
            np.testing.assert_array_equal(p_i<alpha,count_i<=max_count)
        np.savez_compressed(cache,cdfs=cdfs,af=np.stack([r[1] for r in results]),ids=np.array(contract['ids']),
            proposal_indices=np.array([r['replicate'] for r in rows]),p_cumulative=p_c,p_interval=p_i,
            neutral_ge_cumulative=count_c,neutral_ge_interval=count_i,cutoffs_kya=np.arange(1,51),
            interval_left_kya=np.arange(1,50),interval_right_kya=np.arange(2,51))
        save_json(manifest,dict(contract=contract,source=source,inputs=[r[2] for r in results],
            cache_sha256=sha(cache),p_value_version=2,
            comparison='Exact integer pair counts for truth/median; original continuous values for posterior mass.',
            validation='100 unique accepted simulations; verified score/tree/posterior/trajectory hashes, pair panel, 500-kb coordinate, selection schedule, exact AF and all thresholds.'))
        completed.append(dict(case=case,cache=cache.name,manifest=manifest.name,cache_sha256=sha(cache)))
        print(f"Prepared {case['key']}: 100 selected simulations.",flush=True)
    save_json(out/'manifest.json',dict(cases=completed,pending=pending,neutral_cache_sha256=sha(out/'neutral.npz'),
        neutral_source_parameters_sha256=neutral_meta['cache_sha256'],position_bp=500000,pair_count=1000,
        neutral_count=100,neutral_pool='50 held-out neutral_test plus 50 calibration; same pool across all cases.',
        p_definition='One-sided empirical p=(1+# neutrals >= selected)/101; strict p<alpha.',
        thresholds=[.05,.01],multiple_testing='Pointwise only; no adjustment across times or loci.',
        comparison='Exact integer pair counts for truth/median; original continuous values for posterior mass.',
        rows=ROW_ORDER,
        selected_ascertainment='Present-day segregation only; no selection by inferred signal.',workers=workers,
        simulation_configurations={case['key']:(episode_config if case['regime']=='50_30k' else base_configs[case['regime']]) for case in CASES},
        source_simulation_manifests=[dict(path=str(path),sha256=sha(path)) for path in (BASE/'manifest.json',RAW/'manifest.json')]))


def report(out,available_only=False,calls_only=False):
    meta=json.loads((out/'manifest.json').read_text())
    assert meta['cases'], 'No completed cohorts to render.'
    if not available_only:
        assert not meta['pending'] and len(meta['cases'])==4, 'Complete all four cohorts before rendering.'
    assert sha(out/'neutral.npz')==meta['neutral_cache_sha256']
    datasets=[]
    for entry in meta['cases']:
        assert sha(out/entry['cache'])==entry['cache_sha256']
        with np.load(out/entry['cache']) as data:
            datasets.append({k:data[k].copy() for k in ('p_cumulative','p_interval','ids','proposal_indices','cdfs')})
    with np.load(out/'neutral.npz') as data:neutral=data['cdfs'].copy()
    orders={entry['case']['key']:mass_order(data) for entry,data in zip(meta['cases'],datasets)}
    meta['rows']=ROW_ORDER
    save_json(out/'manifest.json',meta)
    save_json(out/'posterior_mass_row_order.json',dict(rule=ROW_ORDER,cases={
        entry['case']['key']:[dict(display_row=row+1,accepted_index=int(i)+1,
            proposal_index=int(data['proposal_indices'][i]),id=str(data['ids'][i]))
            for row,i in enumerate(orders[entry['case']['key']])]
        for entry,data in zip(meta['cases'],datasets)}))
    matplotlib.rcParams.update({'font.size':16,'axes.labelsize':17,'axes.titlesize':17,
        'xtick.labelsize':13,'ytick.labelsize':13,'pdf.fonttype':42,'ps.fonttype':42})
    blue,purple='#0072B2','#AA3377'
    neutral_label='100 segregating neutrals' if meta.get('neutral_conditioned',False) else '100 neutrals'
    summary=[];replicate_records=[]
    # Each underlying statistic curve is retained; bold segments mark its own
    # empirical call. The separate row plots preserve individual call identity.
    for threshold,tag in ([] if calls_only else [(.05,'p005'),(.01,'p001')]):
        with PdfPages(out/f'all_statistic_curves_{tag}.pdf') as book, PdfPages(out/f'posterior_mass_curves_{tag}.pdf') as mass_book:
            for entry,data in zip(meta['cases'],datasets):
                case=entry['case']
                n_selected=len(data['ids'])
                for method,name in enumerate(METHODS):
                    selected=data['cdfs'][:,method]
                    null=neutral[:,method]
                    values=[selected,np.diff(selected,axis=-1)]
                    null_values=[null,np.diff(null,axis=-1)]
                    masks=[data['p_cumulative'][:,method]<threshold,data['p_interval'][:,method]<threshold]
                    fig,axes=plt.subplots(1,2,figsize=(11,8.5))
                    for col,(ax,y,n,mask,color,title) in enumerate(zip(axes,values,null_values,masks,[blue,purple],['Cumulative cutoff','1-kyr interval'])):
                        x=np.arange(1,51) if col==0 else np.arange(1.5,50)
                        ax.axvspan(*case['window'],color='.93',zorder=0)
                        ax.plot(x,n.T,color='.55',alpha=.18,lw=.65,zorder=1)
                        ax.plot(x,y.T,color=color,alpha=.19,lw=.7,zorder=2)
                        for row in range(n_selected):
                            ax.plot(x,np.where(mask[row],y[row],np.nan),color=color,alpha=.65,lw=1.2,marker='o',markersize=1.8,zorder=3)
                        ax.set(title=title,xlim=(0,50.5),xlabel='Time before present (kya)',ylim=(0,None))
                        ax.set_ylabel('Pairwise coalescence fraction' if col==0 else 'Pairwise fraction in 1-kyr bin')
                        ax.set_xticks(np.arange(0,51,10));ax.ticklabel_format(axis='y',style='plain',useOffset=False)
                        ax.grid(alpha=.15)
                        for spine in ('top','right'):ax.spines[spine].set_visible(False)
                    fig.suptitle(case['label'],y=.98,fontsize=21)
                    fig.text(.5,.925,f'{name} · all {n_selected} selected simulations · selected site: 500 kb',ha='center',va='center',fontsize=15)
                    handles=[Line2D([],[],color='.55',lw=1,label=neutral_label),Line2D([],[],color=blue,alpha=.4,lw=1,label=f'{n_selected} selected curves'),
                        Line2D([],[],color=blue,lw=2,marker='o',markersize=4,label=f'Selected points: p < {threshold:g}'),Patch(facecolor='.93',label='Selection active')]
                    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.87),ncol=2,frameon=False,fontsize=12)
                    fig.subplots_adjust(top=.74,bottom=.16,left=.10,right=.97,wspace=.40)
                    fig.text(.5,.06,f'Bold points use each simulation’s own p value against the same {neutral_label}.',ha='center',fontsize=11.5)
                    fig.text(.5,.025,'Pointwise tests; no correction across times/loci. Interval values are plotted at bin midpoints.',ha='center',fontsize=11)
                    stem=f"{case['key']}_{METHOD_KEYS[method]}_curves_{tag}"
                    for extension in ('png','pdf'):fig.savefig(out/f'{stem}.{extension}',dpi=180)
                    book.savefig(fig)
                    if method==2:mass_book.savefig(fig)
                    plt.close(fig)
                print(f"Rendered curves: {case['key']} at p<{threshold:g}.",flush=True)
    for threshold,tag in [(.05,'p005'),(.01,'p001')]:
        with PdfPages(out/f'all_statistics_{tag}.pdf') as book, PdfPages(out/f'posterior_mass_{tag}.pdf') as mass_book:
            for case_entry,data in zip(meta['cases'],datasets):
                case=case_entry['case']
                n_selected=len(data['ids'])
                for method,name in enumerate(METHODS):
                    masks=[data['p_cumulative'][:,method]<threshold,data['p_interval'][:,method]<threshold]
                    order=orders[case['key']] if method==2 else np.arange(n_selected)
                    display_rows=np.argsort(order)+1
                    fig,axes=plt.subplots(1,2,figsize=(11,8.5),sharey=True)
                    for col,(ax,mask,color,title) in enumerate(zip(axes,masks,[blue,purple],['Cumulative cutoff','1-kyr interval'])):
                        lower,upper=case['window']
                        ax.axvspan(lower,upper,color='.91',zorder=0)
                        for row,i in enumerate(order):
                            for start,end in runs(mask[i]):
                                if col==0:
                                    x=np.arange(start+1,end+2)
                                    ax.plot(x,np.full(len(x),row+1),color=color,lw=1.25,marker='o',markersize=1.5)
                                else:
                                    ax.broken_barh([(start+1,end-start+1)],(row+.66,.68),facecolors=color,edgecolors='none')
                        ax.set(title=title,xlim=(0,50.5),ylim=(n_selected+.8,.2),xlabel='Time before present (kya)')
                        ax.set_xticks(np.arange(0,51,10));ax.set_yticks([1,25,50,75,100] if n_selected==100 else np.unique(np.linspace(1,n_selected,min(5,n_selected)).astype(int)))
                        ax.grid(axis='x',alpha=.15)
                        for spine in ('top','right'):ax.spines[spine].set_visible(False)
                        count=int(mask.any(1).sum())
                        ax.text(.5,-.19,f'{count}/{n_selected} with ≥1 significant test',transform=ax.transAxes,ha='center',fontsize=12)
                    axes[0].set_ylabel('Selected simulation (sorted)' if method==2 else 'Selected simulation (proposal order)')
                    fig.suptitle(case['label'],y=.98,fontsize=21)
                    fig.text(.5,.925,f'{name} · pointwise p < {threshold:g} · selected site: 500 kb',ha='center',va='center',fontsize=15)
                    handles=[Line2D([],[],color=blue,marker='o',lw=2,markersize=4,label='Significant cutoff'),
                        Line2D([],[],color=purple,lw=6,label='Significant bin'),Patch(facecolor='.91',label='Selection active')]
                    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.87),ncol=3,frameon=False,fontsize=12)
                    if method==2:
                        fig.text(.5,.075,'Order: first significant time, then longer spans · based on p < 0.05 in both panels',ha='center',fontsize=11)
                    fig.text(.5,.05,f'Same {neutral_label} · 1,000 pairs · p = (1 + neutral exceedances) / 101',ha='center',fontsize=11.5)
                    fig.text(.5,.025,'Blank rows have no pointwise call. No correction across times/loci; bars do not estimate selection duration.',ha='center',fontsize=10.5)
                    fig.subplots_adjust(top=.8,bottom=.21,left=.11,right=.97,wspace=.16)
                    stem=f"{case['key']}_{METHOD_KEYS[method]}_{tag}"
                    for extension in ('png','pdf'):fig.savefig(out/f'{stem}.{extension}',dpi=180)
                    book.savefig(fig)
                    if method==2:mass_book.savefig(fig)
                    plt.close(fig)
                    summary.append(dict(case=case['key'],label=case['label'],method=name,threshold=threshold,n_selected=n_selected,
                        cumulative_any_count=int(masks[0].any(1).sum()),interval_any_count=int(masks[1].any(1).sum()),
                        cumulative_fraction_by_cutoff=masks[0].mean(0).tolist(),interval_fraction_by_bin=masks[1].mean(0).tolist()))
                    for i in range(n_selected):
                        replicate_records.append(dict(case=case['key'],accepted_index=i+1,id=str(data['ids'][i]),
                            display_row=int(display_rows[i]),
                            proposal_index=int(data['proposal_indices'][i]),method=name,threshold=threshold,
                            cumulative_cutoff_runs_kya=[[a+1,b+1] for a,b in runs(masks[0][i])],
                            interval_runs_kya=[[a+1,b+2] for a,b in runs(masks[1][i])]))
                print(f"Rendered {case['key']} at p<{threshold:g}.",flush=True)
    for path in out.glob('*.pdf'):
        contents=path.read_bytes();assert b'/FontFile2' in contents and not re.search(rb'/Subtype\s*/Image',contents),path.name
    save_json(out/'summary.json',dict(results=summary,interpretation='Counts are fractions with any pointwise call, not calibrated power or FPR.',
        script_sha256=sha(Path(__file__)),comparison_source_sha256=sha(Path(__file__).with_name('pointwise_time_tests.py'))))
    save_json(out/'replicate_calls.json',replicate_records)
    (out/'.gitattributes').write_text('* -text\n')
    save_json(out/'output_hashes.json',{p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file() and p.name!='output_hashes.json'})
    print(json.dumps([r for r in summary if r['method']=='Posterior-weighted mass'],indent=2),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase',choices=['prepare','report','all'])
    parser.add_argument('--out',type=Path,default=ROOT/'docs/results/selection_calls_100_20261110')
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--available-only',action='store_true')
    parser.add_argument('--calls-only',action='store_true',help='Regenerate call plots and row mappings without rendering statistic curves.')
    args=parser.parse_args();assert 1<=args.workers<=24
    args.out.mkdir(parents=True,exist_ok=True)
    if args.phase in ('prepare','all'):prepare(args.out,args.workers,args.available_only)
    if args.phase in ('report','all'):report(args.out,args.available_only,args.calls_only)


if __name__=='__main__':main()
