#!/usr/bin/env python3
"""Prepare and plot bounded strong-selection cohorts against the conditioned null.

Zero-survivor treatments retain proposal accounting and have no power estimate.
Partial cohorts are explicitly labelled by their actual number of simulations.
"""
import os
os.environ['MPLBACKEND']='Agg'
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import shutil
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import PercentFormatter
import numpy as np
from plot_selection_episode_series import sha,save_json
from report_selection_arrays import local_path
from pointwise_time_tests import exceedance_counts
from plot_cohort_significant_times import report as report_calls,ROW_ORDER

ROOT=Path(__file__).resolve().parents[1]
REFERENCE=ROOT/'docs/results/selection_calls_100_segregating_null_20261110'


def prepare(args):
    raw=args.raw;out=args.out
    manifest=json.loads((raw/'manifest.json').read_text())
    finished=json.loads((raw/'finished.json').read_text())
    ledger=json.loads((raw/'proposal_ledger.json').read_text())
    ref_meta=json.loads((args.reference/'manifest.json').read_text())
    assert ref_meta['neutral_conditioned'] and sha(args.reference/'neutral.npz')==ref_meta['neutral_cache_sha256']
    historical_path=ROOT/'docs/results/selection_arrays_eas_20261105_4h/manifest.json'
    historical_configs=json.loads(historical_path.read_text())['contract']['configs']
    for stratum in manifest['strata']:
        if manifest['reused'][stratum['key']]:
            previous=historical_configs[stratum['regime']]
            current=manifest['configs'][stratum['regime']]
            assert {k:v for k,v in previous.items() if k!='seed_namespace'}=={k:v for k,v in current.items() if k!='seed_namespace'}
    out.mkdir(parents=True,exist_ok=True)
    for name in ('neutral.npz','neutral_af.npz','neutral_af_manifest.json'):
        dest=out/name
        if dest.exists():assert sha(dest)==sha(args.reference/name)
        else:shutil.copy2(args.reference/name,dest)
    with np.load(out/'neutral.npz') as data:neutral=data['cdfs'].copy();pairs=data['pairs'].copy()
    entries=[];missing=[];accounting=[]
    labels={'10k':'10 kya–present','50k':'50 kya–present','50_30k':'50–30 kya'}
    for s in manifest['strata']:
        receipt_path=raw/(s['key']+'_cohort.json')
        receipt=json.loads(receipt_path.read_text())
        assert receipt['manifest_sha256']==sha(raw/'manifest.json')
        rows=receipt['accepted'];n=len(rows)
        # Any incomplete treatment remains visible in accounting. Do not infer
        # a zero detection rate when no segregating simulation was observed.
        status=next(r for r in finished['strata'] if r['key']==s['key'])
        proposals=ledger[s['key']]
        assert len(proposals)==status['new_completed']
        assert [r['replicate'] for r in proposals]==list(range(len(proposals)))
        for proposal in proposals:
            assert proposal['selection']==s['selection'] and proposal['regime']==s['regime']
            af=proposal['final_frequency']
            assert proposal['accepted']==(0<af<1)
            assert proposal['rejection']==(None if 0<af<1 else 'loss' if af==0 else 'fixation')
        fixation_times=[r['stop_years_ago']/1000 for r in proposals if r['rejection']=='fixation']
        new_survivors=sum(r['accepted'] for r in proposals)
        accounting.append(dict(status,retained=n,complete=receipt['complete'],cohort_sha256=sha(receipt_path),
            detection_estimable=n>0,selection_condition=manifest['conditioning'],
            lost_before_selection_onset=sum(r['rejection']=='loss' and r['stop_years_ago']>=s['start'] for r in proposals),
            lost_after_selection_onset=sum(r['rejection']=='loss' and r['stop_years_ago']<s['start'] for r in proposals),
            fixation_time_kya_quantiles_05_50_95=np.quantile(fixation_times,[.05,.5,.95]).tolist() if fixation_times else None,
            zero_survivor_probability_upper95=(float(1-.05**(1/len(proposals))) if len(proposals) and new_survivors==0 else None)))
        if n==0:
            missing.append(s['key'])
            continue
        assert len({r['observed_sha256'] for r in rows})==n
        def read(r):
            directory=local_path(r['source_directory'])
            for name,key in [('trajectory.tsv','trajectory_sha256'),('scores.npz','scores_sha256'),('posterior.zst','posterior_sha256'),('observed.trees','observed_sha256')]:
                assert sha(directory/name)==r[key],str(directory/name)
            trajectory=np.loadtxt(directory/'trajectory.tsv',skiprows=1)
            assert trajectory.shape==(2001,5) and 0<trajectory[-1,1]<1
            np.testing.assert_array_equal(trajectory[:,0],np.arange(50000,-1,-25))
            np.testing.assert_allclose(trajectory[:,1],trajectory[:,4]/(2*trajectory[:,2]),atol=1e-14)
            ticks=r['selected_reproduction_ticks']
            if isinstance(ticks,str):ticks=json.loads(ticks)
            first=1 if ticks[0]==2 else ticks[0]
            generation=np.arange(1,2002)
            np.testing.assert_allclose(trajectory[:,3],s['selection']*((generation>=first)&(generation<=ticks[1])),atol=1e-12)
            metadata=json.loads((directory/'posterior.zst.meta').read_text())
            np.testing.assert_array_equal(metadata['pairs'],pairs)
            with np.load(directory/'scores.npz') as data:
                np.testing.assert_array_equal(data['positions'],np.arange(100)*10000)
                np.testing.assert_array_equal(data['cutoffs'],np.arange(1000,50001,1000))
                np.testing.assert_array_equal(data['valid_pairs'],1000)
                return trajectory[:,1].copy(),data['cdfs'][:,50].copy()
        with ThreadPoolExecutor(max_workers=args.workers) as pool:parts=list(pool.map(read,rows))
        af=np.stack([r[0] for r in parts]);cdfs=np.stack([r[1] for r in parts])
        ids=[('strong_v1/'+r['id'] if local_path(r['source_directory']).is_relative_to(raw.resolve()) else r['id']) for r in rows]
        assert len(set(ids))==n
        counts_c,counts_i=exceedance_counts(cdfs,neutral)
        key=s['key'].replace('.','p')
        case=dict(key=key,regime=s['regime'],selection=s['selection'],window=[s['end']/1000,s['start']/1000],label=f"{labels[s['regime']]} · s={s['selection']:g}")
        cache=out/(key+'.npz')
        np.savez_compressed(cache,af=af,cdfs=cdfs,ids=np.array(ids),source_ids=np.array([r['id'] for r in rows]),
            proposal_indices=np.arange(n),original_proposal_indices=np.array([r['replicate'] for r in rows]),
            cutoffs_kya=np.arange(1,51),interval_left_kya=np.arange(1,50),interval_right_kya=np.arange(2,51),
            p_cumulative=(1+counts_c)/101,p_interval=(1+counts_i)/101,
            neutral_ge_cumulative=counts_c,neutral_ge_interval=counts_i)
        cache_manifest=out/(key+'_manifest.json')
        save_json(cache_manifest,dict(case=case,contract=dict(ids=ids,position_bp=500000,
            choice=manifest['selection_rule'],conditioning=manifest['conditioning']),inputs=rows,
            cache_sha256=sha(cache),source=dict(path=str(receipt_path),sha256=sha(receipt_path)),
            row_identity='proposal_indices is the cohort ordinal (historical survivors then new proposals); original_proposal_indices retains each source proposal index. IDs disambiguate the two seed namespaces.',
            validation='All score/tree/posterior/AF hashes, selection schedules, AF copy counts, pair panel and focal coordinate checked.'))
        entries.append(dict(case=case,cache=cache.name,manifest=cache_manifest.name,cache_sha256=sha(cache),
            final_population_af=dict(quantile_levels=[0,.25,.5,.75,1],
                quantiles=np.quantile(af[:,-1],[0,.25,.5,.75,1]).tolist(),
                above_99_percent=int(np.count_nonzero(af[:,-1]>.99)),
                conditioning='Whole-population segregation; does not additionally require polymorphism in the 200 sampled haplotypes.')))
    save_json(out/'proposal_accounting.json',dict(strata=accounting,finished=finished,
        interpretation='Rates are conditioned on present-day segregation. Zero surviving replicates means detection cannot be estimated, not 0% power. Rejected fixed trajectories are not used as segregating examples.'))
    save_json(out/'manifest.json',dict(cases=entries,pending=missing,neutral_conditioned=True,
        neutral_cache_sha256=sha(out/'neutral.npz'),neutral_pool='Same 100 original-focal-allele segregation-conditioned neutrals as the weak-selection report.',
        position_bp=500000,pair_count=1000,neutral_count=100,rows=ROW_ORDER,
        simulation_configuration=manifest,raw_manifest_sha256=sha(raw/'manifest.json'),
        reference_manifest_sha256=sha(args.reference/'manifest.json'),
        historical_configuration_sha256=sha(historical_path),
        p_definition='One-sided empirical p=(1+# reference neutrals >= selected)/101; strict thresholds.',
        thresholds=[.05,.01],multiple_testing='Pointwise only; no correction over time or genomic positions.',
        source_script_sha256=sha(Path(__file__))))
    for name in ('schedule.json','finished.json','execution_environment.json'):
        shutil.copy2(raw/name,out/('simulation_'+name))
    schedule=json.loads((raw/'schedule.json').read_text())
    driver=raw/'executed_sources'/f"run_strong_selection_cohorts_{schedule['script_sha256'][:12]}.py"
    assert sha(driver)==schedule['script_sha256']
    shutil.copy2(driver,out/'executed_simulation_driver.py')
    (out/'.gitattributes').write_text('* -text\n')


def report(args):
    out=args.out;meta=json.loads((out/'manifest.json').read_text())
    if meta['cases']:report_calls(out,available_only=True)
    if args.neutral_fpr.exists() and meta['cases']:
        null_manifest=json.loads((args.neutral_fpr/'manifest.json').read_text())
        assert null_manifest['contract']['reference_cdfs_sha256']==meta['neutral_cache_sha256']
        fpr=json.loads((args.neutral_fpr/'summary.json').read_text())['results']
        selected=json.loads((out/'summary.json').read_text())['results']
        benchmark=[]
        for row in selected:
            for kind in ('cumulative','interval'):
                control=next(r for r in fpr if r['method']==row['method'] and r['statistic']==kind and r['threshold']==row['threshold'])
                tp=row[kind+'_any_count'];fn=row['n_selected']-tp
                fp=control['any_time_false_positive_count'];tn=100-fp
                record=dict(case=row['case'],method=row['method'],statistic=kind,threshold=row['threshold'],
                    selected_n=row['n_selected'],neutral_test_n=100,true_positives=tp,false_negatives=fn,
                    false_positives=fp,true_negatives=tn,detection_recall=tp/row['n_selected'],false_positive_rate=fp/100)
                if row['n_selected']==100:
                    record.update(precision=tp/(tp+fp) if tp+fp else None,f1=2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else None)
                else:
                    record.update(precision=None,f1=None,limitation='Incomplete positive cohort; display individual outcomes rather than a balanced-cohort F1.')
                benchmark.append(record)
        save_json(out/'focal_site_detection_benchmark.json',dict(results=benchmark,
            decision_rule='At least one unadjusted pointwise p < threshold anywhere in 1–50 kya at the fixed 500-kb focal site, separately for cumulative cutoffs and intervals.',
            interpretation='Replicate-level selected-versus-neutral detection. Not a temporal-overlap endpoint or a genome-wide F1. F1 is reported only for complete, balanced 100-selected versus 100-neutral-test cohorts.',
            neutral_test_summary_sha256=sha(args.neutral_fpr/'summary.json'),neutral_test_manifest_sha256=sha(args.neutral_fpr/'manifest.json')))
    with np.load(out/'neutral_af.npz') as data:neutral=data['trajectories'][:,:,1].copy()
    assert np.all((neutral[:,-1]>0)&(neutral[:,-1]<1))
    matplotlib.rcParams.update({'font.size':17,'axes.labelsize':19,'xtick.labelsize':16,'ytick.labelsize':16,'pdf.fonttype':42,'ps.fonttype':42})
    accounting=json.loads((out/'proposal_accounting.json').read_text())['strata']
    fig,ax=plt.subplots(figsize=(11,8.5))
    ypos=np.arange(len(accounting));left=np.zeros(len(accounting))
    totals=np.array([r['new_completed'] for r in accounting])
    labels={'10k':'10 kya–present','50k':'50 kya–present','50_30k':'50–30 kya'}
    counts=[np.array([r['lost'] for r in accounting]),np.array([r['fixed'] for r in accounting]),
        np.array([r['accepted']-r['reused'] for r in accounting])]
    np.testing.assert_array_equal(sum(counts),totals)
    for values,color,label in zip(counts,['#999999','#D55E00','#0072B2'],['Lost','Fixed','Segregating']):
        rates=np.divide(values,totals,out=np.zeros(len(totals)),where=totals>0)
        ax.barh(ypos,rates,left=left,color=color,height=.62,label=label)
        left+=rates
    for i,row in enumerate(accounting):
        survivors=row['accepted']-row['reused']
        label=f"{survivors:,} / {row['new_completed']:,}" if row['new_completed'] else '—'
        ax.text(1.025,i,label,va='center',fontsize=14)
        if not row['new_completed']:
            ax.text(.03,i,f"{row['reused']} historical simulations reused; no new proposals",va='center',fontsize=13,color='.35')
    ax.set(yticks=ypos,yticklabels=[f"{labels[r['regime']]}\ns = {r['selection']:g}" for r in accounting],
        xlim=(0,1.27),ylim=(len(accounting)-.45,-.8),xlabel='Fraction of new simulation proposals')
    ax.set_xticks([0,.25,.5,.75,1]);ax.xaxis.set_major_formatter(PercentFormatter(1))
    for spine in ('top','right'):ax.spines[spine].set_visible(False)
    ax.text(1.025,-.60,'Segregating / total',fontsize=11)
    fig.suptitle('Present-day allele fate under selection',fontsize=21,y=.97)
    fig.legend(*ax.get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.5,.90),ncol=3,frameon=False,fontsize=14)
    fig.subplots_adjust(left=.25,right=.98,top=.79,bottom=.19)
    reused=sum(r['reused'] for r in accounting)
    fig.text(.5,.075,f'New proposals only; the selected cohorts also reuse {reused} historical survivors.',ha='center',fontsize=11.5)
    fig.text(.5,.035,'Only segregating alleles enter the coalescence comparison. Fixed and lost proposals are excluded.',ha='center',fontsize=11)
    for ext in ('png','pdf'):fig.savefig(out/f'new_proposal_allele_fates.{ext}',dpi=180)
    plt.close(fig)
    times=np.arange(50000,-1,-25)/1000
    if meta['cases']:
        with PdfPages(out/'allele_frequency_selected_vs_neutral.pdf') as book:
            for entry in meta['cases']:
                with np.load(out/entry['cache']) as data:af=data['af'].copy()
                n=len(af);case=entry['case']
                fig,ax=plt.subplots(figsize=(11,8.5))
                ax.axvspan(*case['window'],color='#E9E3D7',alpha=.55,zorder=0)
                ax.plot(times,neutral.T,color='.4',lw=.9,alpha=.50)
                ax.plot(times,af.T,color='#D55E00',lw=1.5 if n<10 else .85,alpha=.85 if n<10 else .37)
                ax.set(xlim=(50,0),ylim=(-.015,1.015),xlabel='Time before present (kya)',ylabel='Focal-allele frequency')
                ax.yaxis.set_major_formatter(PercentFormatter(1));ax.grid(axis='y',alpha=.18)
                fig.suptitle(case['label'],fontsize=23,y=.97)
                fig.legend(handles=[Line2D([],[],color='#D55E00',lw=2,label=f'{n} selected simulations'),
                    Line2D([],[],color='.4',lw=2,label='100 segregating neutrals'),Patch(facecolor='#E9E3D7',label='Selection active')],
                    loc='upper center',bbox_to_anchor=(.5,.90),ncol=3,frameon=False,fontsize=12)
                fig.subplots_adjust(left=.12,right=.97,top=.79,bottom=.17)
                fig.text(.5,.065,'EAS + 1% introgression at 50 kya · focal site: 500 kb · original focal allele',ha='center',fontsize=12)
                fig.text(.5,.03,'Selected and neutral alleles required to be segregating today; fixed/lost proposals excluded.',ha='center',fontsize=11.5)
                for ext in ('png','pdf'):fig.savefig(out/f"{case['key']}_allele_frequency.{ext}",dpi=180)
                book.savefig(fig);plt.close(fig)
    save_json(out/'output_hashes.json',{p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file() and p.name!='output_hashes.json'})


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('phase',choices=['prepare','report','all'])
    p.add_argument('--raw',type=Path,default=ROOT.parent/'sim/strong_selection_20261110')
    p.add_argument('--reference',type=Path,default=REFERENCE)
    p.add_argument('--out',type=Path,default=ROOT/'docs/results/strong_selection_20261110')
    p.add_argument('--neutral-fpr',type=Path,default=ROOT/'docs/results/segregating_neutral_fpr_100_20261110')
    p.add_argument('--workers',type=int,default=4)
    args=p.parse_args();assert 1<=args.workers<=24
    if args.phase in ('prepare','all'):prepare(args)
    if args.phase in ('report','all'):report(args)


if __name__=='__main__':main()
