#!/usr/bin/env python3
"""Compare significant cumulative cutoffs and 1-kyr intervals at 500 kb.

Use the previously displayed selected examples and exactly the same 100 neutral
simulations. Empirical one-sided p=(1+# neutral >= selected)/101; threshold p<.05.
These are pointwise tests, without correction for searching across times or loci.
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
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
from plot_selection_episode_series import sha,save_json

ROOT=Path(__file__).resolve().parents[1]
SOURCES=[ROOT/'docs/results/selected_site_posteriors_20261106',
         ROOT/'docs/results/ongoing_50k_examples_20261109']
METHODS=['True TMRCA','Posterior median','Posterior-weighted mass']
STEMS=['true_tmrca','posterior_median','posterior_mass']
LABELS=['50–30 kya\ns=0.01','10 kya–present\ns=0.01','50 kya–present\ns=0.01','50 kya–present\ns=0.005']
WINDOWS=[(30,50),(0,10),(0,50),(0,50)]


def load_source(path):
    meta=json.loads((path/'manifest.json').read_text())
    assert meta['contract']['position_bp']==500000
    assert sha(path/'parameters.npz')==meta['cache_sha256']
    with np.load(path/'parameters.npz') as data:
        arrays={k:data[k].copy() for k in ('cdfs','ids','pairs')}
    return arrays,dict(path=str(path),manifest_sha256=sha(path/'manifest.json'),parameters_sha256=meta['cache_sha256'])


def runs(mask):
    """Inclusive index ranges of consecutive True cells, without bridging gaps."""
    indices=np.flatnonzero(mask)
    return [(int(group[0]),int(group[-1])) for group in np.split(indices,np.flatnonzero(np.diff(indices)>1)+1) if len(group)]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=ROOT/'docs/results/pointwise_significant_times_20261110')
    parser.add_argument('--workers',type=int,default=4)
    args=parser.parse_args();assert 1<=args.workers<=24
    args.out.mkdir(parents=True,exist_ok=True)
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        loaded=list(pool.map(load_source,SOURCES))
    first,second=[x[0] for x in loaded]
    np.testing.assert_array_equal(first['pairs'],second['pairs'])
    np.testing.assert_array_equal(first['ids'][3:],second['ids'][2:])
    np.testing.assert_array_equal(first['cdfs'][3:],second['cdfs'][2:])
    neutral=first['cdfs'][3:]
    assert neutral.shape==(100,3,50)
    indices=[(first,0),(first,2),(second,1),(second,0)]
    selected=np.stack([data['cdfs'][i] for data,i in indices])
    ids=[str(data['ids'][i]) for data,i in indices]
    assert ids==['50_30k/introgression/s0.01/selected_00000','10k/introgression/s0.01/selected_00020',
                 '50k/introgression/s0.01/selected_00001','50k/introgression/s0.005/selected_00001']
    assert np.isfinite(selected).all() and np.isfinite(neutral).all()
    assert np.diff(selected,axis=-1).min()>=-1e-12 and np.diff(neutral,axis=-1).min()>=-1e-12
    counts_c=(neutral[None,:,:,:]>=selected[:,None,:,:]).sum(axis=1)
    counts_i=(np.diff(neutral,axis=-1)[None,:,:,:]>=np.diff(selected,axis=-1)[:,None,:,:]).sum(axis=1)
    p_c=(1+counts_c)/101
    p_i=(1+counts_i)/101
    sig_c=p_c<.05;sig_i=p_i<.05
    np.testing.assert_array_equal(sig_c,counts_c<=4)
    np.testing.assert_array_equal(sig_i,counts_i<=4)
    np.savez_compressed(args.out/'p_values.npz',cumulative_p=p_c,interval_p=p_i,
        cumulative_neutral_ge=counts_c,interval_neutral_ge=counts_i,
        cutoffs_kya=np.arange(1,51),interval_left_kya=np.arange(1,50),interval_right_kya=np.arange(2,51),
        selected_ids=np.array(ids),neutral_ids=first['ids'][3:])
    summary=[]
    for i,label in enumerate(LABELS):
        entry=dict(label=label.replace('\n','; '),id=ids[i],selection_window_kya=WINDOWS[i],methods={})
        for j,method in enumerate(METHODS):
            entry['methods'][method]=dict(
                significant_cumulative_cutoff_runs_kya=[[a+1,b+1] for a,b in runs(sig_c[i,j])],
                significant_1kyr_interval_runs_kya=[[a+1,b+2] for a,b in runs(sig_i[i,j])],
                min_cumulative_p=float(p_c[i,j].min()),min_interval_p=float(p_i[i,j].min()))
        summary.append(entry)
    matplotlib.rcParams.update({'font.size':16,'axes.labelsize':17,'axes.titlesize':19,
        'xtick.labelsize':14,'ytick.labelsize':15,'pdf.fonttype':42,'ps.fonttype':42})
    blue,purple='#0072B2','#AA3377'
    handles=[Line2D([],[],color=blue,marker='o',lw=2,markersize=5,label='Cumulative cutoff'),
        Line2D([],[],color=purple,lw=7,label='1-kyr interval'),Patch(facecolor='.9',label='Selection active')]
    with PdfPages(args.out/'all_statistics.pdf') as book:
        for j,(method,stem) in enumerate(zip(METHODS,STEMS)):
            fig,ax=plt.subplots(figsize=(11,8.5))
            for i,(lower,upper) in enumerate(WINDOWS):
                y=3-i
                ax.broken_barh([(lower,upper-lower)],(y-.38,.76),facecolors='.9',edgecolors='none',zorder=0)
                for start,end in runs(sig_c[i,j]):
                    x=np.arange(start+1,end+2)
                    ax.plot(x,np.full(len(x),y+.15),color=blue,marker='o',lw=2,markersize=4)
                for start,end in runs(sig_i[i,j]):
                    ax.broken_barh([(start+1,end-start+1)],(y-.24,.18),facecolors=purple,edgecolors='none')
                if not sig_c[i,j].any():ax.text(1,y+.15,'None',va='center',color=blue,fontsize=13)
                if not sig_i[i,j].any():ax.text(1,y-.15,'None',va='center',color=purple,fontsize=13)
            ax.set(xlim=(0,50.5),ylim=(-.6,3.6),xlabel='Time before present (kya)')
            ax.set_yticks([3,2,1,0],LABELS)
            ax.tick_params(axis='y',length=0,pad=12)
            ax.set_xticks(np.arange(0,51,5));ax.grid(axis='x',alpha=.18)
            for spine in ('top','right','left'):ax.spines[spine].set_visible(False)
            fig.suptitle(f'Pointwise p < 0.05: {method.lower()}',y=.97,fontsize=20)
            fig.text(.5,.92,'Selected site: 500 kb · 100 neutrals · same 1,000 pairs for all methods',ha='center',fontsize=13)
            fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.89),ncol=3,frameon=False,fontsize=13)
            fig.text(.5,.085,'Dots: tested cumulative cutoffs. Bars: significant 1-kyr bins (1–2 through 49–50 kya).',ha='center',fontsize=11)
            fig.text(.5,.055,'p = (1 + number of neutrals ≥ selected) / 101; p < 0.05 requires at most 4 neutrals.',ha='center',fontsize=11)
            fig.text(.5,.025,'Pointwise comparisons; no correction across times or loci. These are not inferred selection durations.',ha='center',fontsize=11)
            fig.subplots_adjust(left=.25,right=.97,top=.8,bottom=.19)
            for extension in ('png','pdf'):fig.savefig(args.out/f'{stem}.{extension}',dpi=180)
            book.savefig(fig);plt.close(fig)
    for path in args.out.glob('*.pdf'):
        contents=path.read_bytes();assert b'/FontFile2' in contents and not re.search(rb'/Subtype\s*/Image',contents)
    save_json(args.out/'summary.json',summary)
    save_json(args.out/'manifest.json',dict(inputs=[x[1] for x in loaded],selected_ids=ids,
        position_bp=500000,methods=METHODS,selection_windows=WINDOWS,
        grid='Cumulative cutoffs 1..50 kya; intervals [1,2),...,[49,50) kya.',
        p_definition='(1 + number of neutral simulations with statistic >= selected) / 101',threshold=.05,comparison='strictly less',
        multiple_testing='None; pointwise descriptive comparison, not a scan-adjusted call or timing confidence interval.',
        neutral_pool='Same 50 held-out plus 50 calibration simulations across all scenarios.',
        ascertainment='50–30 and ongoing examples first surviving proposals; 10-kya example previously chosen for a mass call at 420 kb.',
        validation='Verified both caches, identical neutral IDs/CDFs and identical pair panel; threshold masks agree exactly with neutral exceedance counts <=4.',
        script_sha256=sha(Path(__file__)),workers=args.workers,versions=dict(numpy=np.__version__,matplotlib=matplotlib.__version__)))
    (args.out/'.gitattributes').write_text('* -text\n')
    save_json(args.out/'output_hashes.json',{p.name:sha(p) for p in sorted(args.out.iterdir()) if p.is_file() and p.name!='output_hashes.json'})
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':main()
