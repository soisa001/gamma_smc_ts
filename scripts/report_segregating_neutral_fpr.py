#!/usr/bin/env python3
"""Independent 100-vs-100 conditioned-neutral false-positive evaluation.

prepare verifies raw test scores and caches the exact focal-site comparison.
report requires only the committed cache and writes PNG/vector PDF figures.
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
import shutil
import sys
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.ticker import PercentFormatter
import numpy as np
from pointwise_time_tests import exceedance_counts
from plot_selection_episode_series import sha,save_json
from plot_pointwise_significant_times import runs
from plot_cohort_significant_times import mass_order
from report_selection_arrays import local_path

ROOT=Path(__file__).resolve().parents[1]
REFERENCE=ROOT/'docs/results/selection_calls_100_segregating_null_20261110'
RAW=ROOT.parent/'sim/segregating_neutral_fpr_100_20261110'
OUT=ROOT/'docs/results/segregating_neutral_fpr_100_20261110'
METHODS=['True TMRCA','Posterior median','Posterior-weighted mass']
COLORS=['#009E73','#D55E00','#0072B2']


def wilson(count,n=100):
    """95% Wilson binomial interval, conditional on the fixed reference bank."""
    z=1.959963984540054
    p=np.asarray(count)/n;denom=1+z*z/n
    center=(p+z*z/(2*n))/denom
    radius=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/denom
    # Exact boundaries avoid tiny negative error-bar lengths from roundoff.
    lower=np.where(p==0,0,np.maximum(0,center-radius))
    upper=np.where(p==1,1,np.minimum(1,center+radius))
    return lower,upper


def prepare(args):
    receipt_path=args.raw/'cohort.json'
    receipt=json.loads(receipt_path.read_text())
    contract=json.loads((args.raw/'manifest.json').read_text())
    assert contract['independent_test'] and not receipt['reused']
    assert receipt['manifest_sha256']==sha(args.raw/'manifest.json')
    records=receipt['accepted']
    assert len(records)==receipt['target']==100
    ref_meta=json.loads((args.reference/'manifest.json').read_text())
    ref_receipt=json.loads((args.reference/'neutral_cohort_receipt.json').read_text())
    assert ref_meta['neutral_conditioned'] and sha(args.reference/'neutral.npz')==ref_meta['neutral_cache_sha256']
    ref_af_meta=json.loads((args.reference/'neutral_af_manifest.json').read_text())
    assert sha(args.reference/'neutral_af.npz')==ref_af_meta['cache_sha256']
    with np.load(args.reference/'neutral.npz') as data:
        ref_cdfs=data['cdfs'].copy();ref_ids=data['ids'].copy();pairs=data['pairs'].copy()
    with np.load(args.reference/'neutral_af.npz') as data:
        np.testing.assert_array_equal(ref_ids,data['ids'])
        ref_af=data['trajectories'][:,:,1].copy()
    np.testing.assert_array_equal(pairs,np.loadtxt(args.raw/'pairs.tsv',dtype=int))
    assert not set(r['id'] for r in records)&set(ref_ids)
    assert not set(r['observed_sha256'] for r in records)&set(r['observed_sha256'] for r in ref_receipt['accepted'])
    assert not set(r['seed'] for r in records)&set(r['seed'] for r in ref_receipt['accepted'])
    # Assert independent deterministic phase seeds, beyond top-level receipt IDs.
    from run_time_resolved_pilot import seed
    for phase in ('ancestry','pairing','slim','sample','mutations'):
        test_seeds={seed(r['seed_namespace']+r['unit']+':'+phase) for r in records}
        reference_seeds={seed(('recent-v1:' if not r.get('seed_namespace') else r['seed_namespace'])+r['unit']+':'+phase) for r in ref_receipt['accepted']}
        assert len(test_seeds)==100 and not test_seeds&reference_seeds,phase
    def read(r):
        directory=local_path(r['source_directory'])
        for name,key in [('trajectory.tsv','trajectory_sha256'),('observed.trees','observed_sha256'),('posterior.zst','posterior_sha256'),('scores.npz','scores_sha256')]:
            assert sha(directory/name)==r[key],str(directory/name)
        trajectory=np.loadtxt(directory/'trajectory.tsv',skiprows=1)
        assert trajectory.shape==(2001,5)
        np.testing.assert_array_equal(trajectory[:,0],np.arange(50000,-1,-25))
        np.testing.assert_array_equal(trajectory[:,3],0)
        np.testing.assert_allclose(trajectory[:,1],trajectory[:,4]/(2*trajectory[:,2]),atol=1e-14)
        assert 0<trajectory[-1,1]<1 and r['selection']==0 and r['accepted']
        posterior_meta=json.loads((directory/'posterior.zst.meta').read_text())
        np.testing.assert_array_equal(posterior_meta['pairs'],pairs)
        with np.load(directory/'scores.npz') as data:
            np.testing.assert_array_equal(data['positions'],np.arange(100)*10000)
            np.testing.assert_array_equal(data['cutoffs'],np.arange(1000,50001,1000))
            np.testing.assert_array_equal(data['valid_pairs'],1000)
            return data['cdfs'].copy(),trajectory[:,1].copy()
    cache_contract=dict(test_cohort_sha256=sha(receipt_path),reference_cdfs_sha256=sha(args.reference/'neutral.npz'),
        reference_af_sha256=sha(args.reference/'neutral_af.npz'),test_ids=[r['id'] for r in records],
        position_bp=500000,region_length_bp=1000000,scanned_positions_bp=list(range(0,1000000,10000)),
        conditioning='Original introduced focal allele segregating today in both test and reference cohorts.',
        test_ascertainment='First 100 accepted new proposals in fixed index order, independent seeds, no reuse or signal ascertainment.',
        p_definition='(1 + number of reference neutrals >= test neutral) / 101; exact integer comparisons for true TMRCA and median.',
        comparison_source_sha256=sha(Path(__file__).with_name('pointwise_time_tests.py')))
    args.out.mkdir(parents=True,exist_ok=True)
    cache=args.out/'parameters.npz';manifest=args.out/'manifest.json'
    if manifest.exists():
        saved=json.loads(manifest.read_text())
        assert saved['contract']==cache_contract and sha(cache)==saved['cache_sha256']
        print('Verified existing independent FPR cache.',flush=True)
        return
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        parts=list(pool.map(read,records));ref_parts=list(pool.map(read,ref_receipt['accepted']))
    test_all=np.stack([p[0] for p in parts]);reference_all=np.stack([p[0] for p in ref_parts])
    test=test_all[:,:,50,:];test_af=np.stack([p[1] for p in parts])
    assert test_all.shape==reference_all.shape==(100,3,100,50)
    np.testing.assert_array_equal(ref_cdfs,reference_all[:,:,50,:])
    np.testing.assert_array_equal(ref_af,np.stack([p[1] for p in ref_parts]))
    assert test.shape==ref_cdfs.shape==(100,3,50) and len(set(r['id'] for r in records))==100
    assert len({r['observed_sha256'] for r in records})==100
    counts_c,counts_i=exceedance_counts(test,ref_cdfs)
    np.savez_compressed(cache,test_cdfs=test,reference_cdfs=ref_cdfs,test_af=test_af,reference_af=ref_af,
        test_ids=np.array([r['id'] for r in records]),reference_ids=ref_ids,pairs=pairs,
        proposal_indices=np.array([r['replicate'] for r in records]),cutoffs_kya=np.arange(1,51),
        test_all_cdfs=test_all,reference_all_cdfs=reference_all,positions_bp=np.arange(100)*10000,
        p_cumulative=(1+counts_c)/101,p_interval=(1+counts_i)/101,
        reference_ge_cumulative=counts_c,reference_ge_interval=counts_i)
    save_json(manifest,dict(contract=cache_contract,cache_sha256=sha(cache),workers=args.workers,
        simulation_configuration=contract,inputs=records,reference_receipt_sha256=sha(args.reference/'neutral_cohort_receipt.json'),
        validation='Independent cohort IDs, observed-tree hashes, top-level and phase seeds; original focal AF strictly between 0 and 1; s=0; all input hashes, pair panels, cutoffs, coordinates and exact copy-count AF checked.',
        script_sha256=sha(Path(__file__)),versions=dict(python=sys.version,numpy=np.__version__,matplotlib=matplotlib.__version__)))
    for name in ('cohort.json','progress.json','schedule.json','execution_environment.json'):
        shutil.copy2(args.raw/name,args.out/('test_'+name))
    (args.out/'.gitattributes').write_text('* -text\n')


def report(args):
    out=args.out
    meta=json.loads((out/'manifest.json').read_text())
    assert sha(out/'parameters.npz')==meta['cache_sha256']
    with np.load(out/'parameters.npz') as cached:data={k:cached[k].copy() for k in cached.files}
    counts=exceedance_counts(data['test_cdfs'],data['reference_cdfs'])
    for kind,count in zip(('cumulative','interval'),counts):
        np.testing.assert_array_equal(data['reference_ge_'+kind],count)
        np.testing.assert_array_equal(data['p_'+kind],(count+1)/101)
    def position_counts(i):
        return exceedance_counts(data['test_all_cdfs'][:,:,i,:],data['reference_all_cdfs'][:,:,i,:])
    with ThreadPoolExecutor(max_workers=args.workers) as pool:regional_counts=list(pool.map(position_counts,range(100)))
    regional={kind:(1+np.stack([r[i] for r in regional_counts],axis=2))/101 for i,kind in enumerate(('cumulative','interval'))}
    for kind in ('cumulative','interval'):np.testing.assert_array_equal(regional[kind][:,:,50,:],data['p_'+kind])
    np.savez_compressed(out/'regional_p_values.npz',p_cumulative=regional['cumulative'],p_interval=regional['interval'],positions_bp=data['positions_bp'])
    matplotlib.rcParams.update({'font.size':16,'axes.labelsize':17,'axes.titlesize':17,
        'xtick.labelsize':13,'ytick.labelsize':13,'pdf.fonttype':42,'ps.fonttype':42})
    rows=[];replicates=[]
    def save(fig,stem,book=None):
        for ext in ('png','pdf'):fig.savefig(out/f'{stem}.{ext}',dpi=180)
        if book is not None:book.savefig(fig)
        plt.close(fig)
    def clean(ax):
        ax.grid(axis='y',alpha=.18)
        for spine in ('top','right'):ax.spines[spine].set_visible(False)
    with PdfPages(out/'fpr_pointwise.pdf') as book:
        for alpha,tag in [(.05,'p005'),(.01,'p001')]:
            fig,axes=plt.subplots(1,2,figsize=(11,8.5),sharey=True)
            for col,(ax,kind,title) in enumerate(zip(axes,('cumulative','interval'),('Cumulative cutoff','1-kyr interval'))):
                x=np.arange(1,51) if col==0 else np.arange(1.5,50)
                for method,(label,color) in enumerate(zip(METHODS,COLORS)):
                    mask=data['p_'+kind][:,method]<alpha
                    counts=mask.sum(0);rate=counts/100
                    ax.plot(x,rate,color=color,lw=2,label=label)
                    any_count=int(mask.any(1).sum());lo,hi=wilson(any_count)
                    rows.append(dict(method=label,statistic=kind,threshold=alpha,n_test=100,
                        any_time_false_positive_count=any_count,any_time_fpr=any_count/100,
                        any_time_wilson95=[float(lo),float(hi)],false_positive_count_by_time=counts.tolist(),
                        fpr_by_time=rate.tolist(),mean_pointwise_fpr=float(rate.mean()),max_pointwise_fpr=float(rate.max())))
                ax.axhline(alpha,color='.4',ls='--',lw=1.4)
                ax.set(title=title,xlim=(0,50),xlabel='Time before present (kya)')
                ax.yaxis.set_major_formatter(PercentFormatter(1,decimals=0));clean(ax)
            axes[0].set_ylabel('Test neutrals called significant')
            ymax=max(.15,max(np.max(line.get_ydata()) for ax in axes for line in ax.lines)*1.18)
            axes[0].set_ylim(0,ymax)
            fig.suptitle(f'Pointwise false-positive rate · p < {alpha:g}',fontsize=21,y=.97)
            handles=[Line2D([],[],color=c,lw=2,label=m) for c,m in zip(COLORS,METHODS)]
            fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.90),ncol=3,frameon=False,fontsize=12)
            fig.subplots_adjust(left=.11,right=.97,top=.78,bottom=.19,wspace=.20)
            fig.text(.5,.08,'100 independent test neutrals vs 100 fixed reference neutrals · focal site: 500 kb',ha='center',fontsize=12)
            fig.text(.5,.04,'Both cohorts segregating today · dashed line: nominal threshold · 1,000 pairs',ha='center',fontsize=12)
            save(fig,f'fpr_pointwise_{tag}',book)
    fig,axes=plt.subplots(1,2,figsize=(11,8.5),sharey=True)
    for ax,alpha in zip(axes,[.05,.01]):
        for offset,kind,color,label in [(-.18,'cumulative','#0072B2','Cumulative: any cutoff'),(.18,'interval','#AA3377','Interval: any bin')]:
            selected=[next(r for r in rows if r['method']==m and r['statistic']==kind and r['threshold']==alpha) for m in METHODS]
            heights=np.array([r['any_time_fpr'] for r in selected]);cis=np.array([r['any_time_wilson95'] for r in selected])
            x=np.arange(3)+offset
            ax.bar(x,heights,width=.34,color=color,label=label)
            ax.errorbar(x,heights,yerr=np.array([heights-cis[:,0],cis[:,1]-heights]),fmt='none',ecolor='.2',capsize=4)
            for xx,h,ci in zip(x,heights,cis):ax.text(xx,ci[1]+.025,f'{h:.0%}',ha='center',fontsize=13)
        ax.set(title=f'p < {alpha:g}',ylim=(0,1.09),xticks=np.arange(3),xticklabels=['True\nTMRCA','Posterior\nmedian','Posterior\nmass'])
        ax.yaxis.set_major_formatter(PercentFormatter(1,decimals=0));clean(ax)
    axes[0].set_ylabel('Test neutrals with ≥1 false-positive call')
    fig.suptitle('False positives when searching across 1–50 kya',fontsize=21,y=.97)
    fig.legend(*axes[0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.5,.90),ncol=2,frameon=False,fontsize=13)
    fig.subplots_adjust(left=.11,right=.97,top=.79,bottom=.22,wspace=.18)
    fig.text(.5,.095,'100 independent test neutrals · bars: observed FPR · whiskers: 95% Wilson interval',ha='center',fontsize=12)
    fig.text(.5,.055,'Unadjusted pointwise p values; any-time FPR need not equal 5% or 1%.',ha='center',fontsize=12)
    fig.text(.5,.02,'One focal site. Intervals describe test-cohort uncertainty conditional on this fixed reference pool.',ha='center',fontsize=10.5)
    save(fig,'fpr_any_time')
    regional_rows=[]
    fig,axes=plt.subplots(1,2,figsize=(11,8.5),sharey=True)
    for ax,alpha in zip(axes,[.05,.01]):
        for offset,kind,color,label in [(-.18,'cumulative','#0072B2','Cumulative cutoffs'),(.18,'interval','#AA3377','1-kyr intervals')]:
            rates=[];intervals=[]
            for method,name in enumerate(METHODS):
                mask=regional[kind][:,method]<alpha
                count=int(mask.any(axis=(1,2)).sum());lo,hi=wilson(count)
                rates.append(count/100);intervals.append([float(lo),float(hi)])
                regional_rows.append(dict(method=name,statistic=kind,threshold=alpha,false_positive_count=count,
                    fpr=count/100,wilson95=[float(lo),float(hi)],any_time_fpr_by_position=mask.any(2).mean(0).tolist()))
            rates=np.array(rates);intervals=np.array(intervals);x=np.arange(3)+offset
            ax.bar(x,rates,width=.34,color=color,label=label)
            ax.errorbar(x,rates,yerr=np.array([rates-intervals[:,0],intervals[:,1]-rates]),fmt='none',ecolor='.2',capsize=4)
            label_pad=.025 if kind=='cumulative' else .075
            for xx,h,ci in zip(x,rates,intervals):ax.text(xx,ci[1]+label_pad,f'{h:.0%}',ha='center',fontsize=13)
        ax.set(title=f'p < {alpha:g}',ylim=(0,1.15),xticks=np.arange(3),xticklabels=['True\nTMRCA','Posterior\nmedian','Posterior\nmass'])
        ax.yaxis.set_major_formatter(PercentFormatter(1,decimals=0));clean(ax)
    axes[0].set_ylabel('Neutral regions with ≥1 false-positive call')
    fig.suptitle('False positives across the 1-Mb region and 1–50 kya',fontsize=20,y=.97)
    fig.legend(*axes[0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.5,.90),ncol=2,frameon=False,fontsize=13)
    fig.subplots_adjust(left=.11,right=.97,top=.79,bottom=.22,wspace=.18)
    fig.text(.5,.095,'100 independent test regions · 100 positions at 10-kb strides · 95% Wilson intervals',ha='center',fontsize=12)
    fig.text(.5,.055,'At least one call at any tested position/time; no correction over positions or times.',ha='center',fontsize=12)
    fig.text(.5,.02,'Both cohorts conditioned on the original 500-kb focal allele segregating today.',ha='center',fontsize=11)
    save(fig,'fpr_any_position_and_time')
    order=mass_order(data);inverse=np.argsort(order)+1
    with PdfPages(out/'neutral_false_calls.pdf') as book:
        for alpha,tag in [(.05,'p005'),(.01,'p001')]:
            for method,label in enumerate(METHODS):
                ordering=order if method==2 else np.arange(100)
                fig,axes=plt.subplots(1,2,figsize=(11,8.5),sharey=True)
                for col,(ax,kind,color,title) in enumerate(zip(axes,('cumulative','interval'),('#0072B2','#AA3377'),('Cumulative cutoff','1-kyr interval'))):
                    mask=data['p_'+kind][:,method]<alpha
                    for row,i in enumerate(ordering):
                        for start,end in runs(mask[i]):
                            if col==0:
                                x=np.arange(start+1,end+2);ax.plot(x,np.full(len(x),row+1),color=color,lw=1.2,marker='o',markersize=1.5)
                            else:ax.broken_barh([(start+1,end-start+1)],(row+.66,.68),facecolors=color,edgecolors='none')
                    ax.set(title=title,xlim=(0,50.5),ylim=(100.8,.2),xlabel='Time before present (kya)')
                    ax.set_yticks([1,25,50,75,100]);ax.grid(axis='x',alpha=.15)
                    ax.text(.5,-.19,f'{int(mask.any(1).sum())}/100 with ≥1 false call',transform=ax.transAxes,ha='center',fontsize=12)
                    for i in range(100):
                        replicates.append(dict(id=str(data['test_ids'][i]),proposal_index=int(data['proposal_indices'][i]),
                            method=label,statistic=kind,threshold=alpha,display_row=int(inverse[i] if method==2 else i+1),
                            called=bool(mask[i].any()),runs_kya=[[a+1,b+(1 if col==0 else 2)] for a,b in runs(mask[i])]))
                axes[0].set_ylabel('Neutral test simulation (sorted)' if method==2 else 'Neutral test simulation (proposal order)')
                fig.suptitle(f'{label} · false calls at p < {alpha:g}',fontsize=20,y=.97)
                fig.text(.5,.90,'100 independent test neutrals vs 100 reference neutrals · both segregating today',ha='center',fontsize=13)
                fig.subplots_adjust(top=.80,bottom=.21,left=.12,right=.97,wspace=.16)
                fig.text(.5,.055,'Every colored mark is a false positive: all test simulations have s = 0.',ha='center',fontsize=12)
                fig.text(.5,.025,'Posterior mass sorted using p < 0.05 calls, same rows at both thresholds. No multiple-testing correction.',ha='center',fontsize=10.5)
                save(fig,f"neutral_false_calls_{['truth','median','mass'][method]}_{tag}",book)
    fig,ax=plt.subplots(figsize=(11,8.5))
    times=np.arange(50000,-1,-25)/1000
    ax.plot(times,data['reference_af'].T,color='.5',alpha=.35,lw=.85)
    ax.plot(times,data['test_af'].T,color='#0072B2',alpha=.35,lw=.85)
    ax.set(xlim=(50,0),ylim=(-.015,1.015),xlabel='Time before present (kya)',ylabel='Original focal-allele frequency')
    ax.yaxis.set_major_formatter(PercentFormatter(1));clean(ax)
    fig.suptitle('Independent neutral cohorts · allele-frequency trajectories',fontsize=19,y=.97)
    fig.legend(handles=[Line2D([],[],color='.5',lw=2,label='100 reference neutrals'),Line2D([],[],color='#0072B2',lw=2,label='100 test neutrals')],
        loc='upper center',bbox_to_anchor=(.5,.90),ncol=2,frameon=False,fontsize=14)
    fig.subplots_adjust(left=.12,right=.97,top=.80,bottom=.16)
    fig.text(.5,.055,'EAS + 1% introgression at 50 kya · s = 0 · all 200 focal alleles segregating today',ha='center',fontsize=12)
    save(fig,'neutral_test_vs_reference_af')
    combined=[]
    for alpha in (.05,.01):
        for method,label in enumerate(METHODS):
            mask=(data['p_cumulative'][:,method]<alpha).any(1)|(data['p_interval'][:,method]<alpha).any(1)
            count=int(mask.sum());lo,hi=wilson(count)
            combined.append(dict(method=label,threshold=alpha,false_positive_count=count,fpr=count/100,wilson95=[float(lo),float(hi)]))
    for path in out.glob('*.pdf'):
        content=path.read_bytes()
        assert content.rstrip().endswith(b'%%EOF') and b'/FontFile2' in content and not re.search(rb'/Subtype\s*/Image',content),path.name
    save_json(out/'summary.json',dict(results=rows,either_statistic_any_time=combined,regional_any_position_and_time=regional_rows,
        interpretation='Independent neutral-test FPR: pointwise and unadjusted any-time calls at the fixed focal site, plus a separate any-position/any-time regional endpoint. Wilson intervals condition on this fixed reference bank and do not include variability from resampling that bank.',
        n_test=100,n_reference=100,position_bp=500000,condition=meta['contract']['conditioning'],script_sha256=sha(Path(__file__))))
    save_json(out/'replicate_calls.json',replicates)
    save_json(out/'output_hashes.json',{f.name:sha(f) for f in sorted(out.iterdir()) if f.is_file() and f.name!='output_hashes.json'})
    print(json.dumps([{k:v for k,v in r.items() if not isinstance(v,list)} for r in rows],indent=2),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('phase',choices=['prepare','report','all'])
    p.add_argument('--raw',type=Path,default=RAW)
    p.add_argument('--reference',type=Path,default=REFERENCE)
    p.add_argument('--out',type=Path,default=OUT)
    p.add_argument('--workers',type=int,default=4)
    args=p.parse_args();assert 1<=args.workers<=24
    if args.phase in ('prepare','all'):prepare(args)
    if args.phase in ('report','all'):report(args)


if __name__=='__main__':main()
