#!/usr/bin/env python3
"""Plot the existing 100-selected cohorts against their same 100 neutral replicates.

Uses verified portable caches only; no simulations or random draws. AF tracks
the original allele introduced at 500 kb, including neutral losses/fixation.
"""
import os
os.environ['MPLBACKEND']='Agg'
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import re
import sys
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import PercentFormatter
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
DEFAULT_OUT=ROOT/'docs/results/selection_calls_100_20261110'
DEFAULT_NEUTRAL=ROOT/'docs/results/selected_site_posteriors_20261106'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save_json(path,value):
    path.write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8',newline='\n')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=DEFAULT_OUT)
    parser.add_argument('--neutral-source',type=Path,default=DEFAULT_NEUTRAL)
    parser.add_argument('--workers',type=int,default=4)
    args=parser.parse_args()
    assert 1<=args.workers<=24
    out=args.out
    meta=json.loads((out/'manifest.json').read_text())
    assert len(meta['cases'])==4 and not meta['pending']
    null_meta=json.loads((args.neutral_source/'neutral_af_manifest.json').read_text())
    null_path=args.neutral_source/'neutral_af.npz'
    assert sha(null_path)==null_meta['cache_sha256']
    assert sha(out/'neutral.npz')==meta['neutral_cache_sha256']
    with np.load(out/'neutral.npz') as data:expected_ids=data['ids'].copy()
    with np.load(null_path) as data:
        trajectories=data['trajectories'].copy()
        neutral_ids=data['ids'].copy()
    np.testing.assert_array_equal(neutral_ids,expected_ids)
    assert len(set(neutral_ids))==100 and trajectories.shape==(100,2001,5)
    years=np.arange(50000,-1,-25)
    times=years/1000
    np.testing.assert_array_equal(trajectories[:,:,0],np.broadcast_to(years,(100,2001)))
    np.testing.assert_array_equal(trajectories[:,:,3],0)
    np.testing.assert_allclose(trajectories[:,:,1],trajectories[:,:,4]/(2*trajectories[:,:,2]),atol=1e-14)
    neutral=trajectories[:,:,1]
    assert np.isfinite(neutral).all() and np.all((neutral>=0)&(neutral<=1))
    conditioned=meta.get('neutral_conditioned',False)
    if conditioned:assert np.all((neutral[:,-1]>0)&(neutral[:,-1]<1))

    def read(entry):
        path=out/entry['cache']
        manifest=json.loads((out/entry['manifest']).read_text())
        assert sha(path)==entry['cache_sha256']==manifest['cache_sha256']
        with np.load(path) as data:
            af=data['af'].copy();ids=data['ids'].copy()
        np.testing.assert_array_equal(ids,manifest['contract']['ids'])
        assert len(set(ids))==100 and af.shape==(100,2001)
        assert np.isfinite(af).all() and np.all((af>=0)&(af<=1))
        assert np.all((af[:,-1]>0)&(af[:,-1]<1))
        return entry,af,ids
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        selected=list(pool.map(read,meta['cases']))

    matplotlib.rcParams.update({'font.size':17,'axes.labelsize':19,'axes.titlesize':20,
        'xtick.labelsize':16,'ytick.labelsize':16,'pdf.fonttype':42,'ps.fonttype':42})
    orange='#D55E00';gray='#666666';shade='#E9E3D7'
    handles=[Line2D([],[],color=orange,lw=2,label='100 selected simulations'),
        Line2D([],[],color=gray,lw=2,label='100 segregating neutrals' if conditioned else '100 neutral simulations'),
        Patch(facecolor=shade,label='Selection active')]
    def draw(ax,case,af):
        ax.axvspan(*case['window'],facecolor=shade,alpha=.55,zorder=0)
        selected_lines=ax.plot(times,af.T,color=orange,lw=.85,alpha=.37,zorder=2)
        neutral_lines=ax.plot(times,neutral.T,color=gray,lw=.9,alpha=.50,zorder=3)
        assert len(selected_lines)==len(neutral_lines)==100
        ax.set(xlim=(50,0),ylim=(-.015,1.015),xlabel='Time before present (kya)',ylabel='Focal-allele frequency')
        ax.set_xticks(np.arange(0,51,10));ax.set_yticks(np.arange(0,1.01,.2))
        ax.yaxis.set_major_formatter(PercentFormatter(1,decimals=0))
        ax.grid(axis='y',alpha=.18)
        for spine in ('top','right'):ax.spines[spine].set_visible(False)

    outputs=[]
    with PdfPages(out/'allele_frequency_100_selected_vs_neutral.pdf') as book:
        for entry,af,ids in selected:
            case=entry['case']
            fig,ax=plt.subplots(figsize=(11,8.5))
            draw(ax,case,af)
            fig.suptitle(case['label'],fontsize=23,y=.97)
            fig.text(.5,.915,'Allele-frequency trajectories at the focal site (500 kb)',ha='center',fontsize=17)
            fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.875),ncol=3,frameon=False,fontsize=12)
            fig.subplots_adjust(left=.12,right=.97,bottom=.18,top=.79)
            fig.text(.5,.07,'EAS demography · 1% introgression at 50 kya · same cohorts as the coalescence plots',ha='center',fontsize=12)
            conditioning_note=('Selected and neutral: original introduced focal allele required to be segregating today.' if conditioned
                else 'Selected: segregating today. Neutral: original introduced allele, including losses and fixation.')
            fig.text(.5,.035,conditioning_note,ha='center',fontsize=11.5)
            stem=f"{case['key']}_allele_frequency_100"
            for extension in ('png','pdf'):
                path=out/f'{stem}.{extension}'
                fig.savefig(path,dpi=180)
                outputs.append(path)
            book.savefig(fig)
            plt.close(fig)
            print(f"Rendered {case['key']}: 100 selected + 100 neutral trajectories.",flush=True)
    outputs.append(out/'allele_frequency_100_selected_vs_neutral.pdf')
    for path in outputs:
        if path.suffix=='.pdf':
            contents=path.read_bytes()
            assert b'/FontFile2' in contents and not re.search(rb'/Subtype\s*/Image',contents)
            assert contents.rstrip().endswith(b'%%EOF')
            expected=4 if path.name=='allele_frequency_100_selected_vs_neutral.pdf' else 1
            assert len(re.findall(rb'/Type\s*/Page\b',contents))==expected
    final=neutral[:,-1]
    provenance=dict(
        source_manifest=dict(path=str(out/'manifest.json'),sha256=sha(out/'manifest.json')),
        neutral_source=dict(path=str(null_path),sha256=sha(null_path),manifest_sha256=sha(args.neutral_source/'neutral_af_manifest.json')),
        selected_sources=[dict(case=entry['case'],cache=entry['cache'],sha256=entry['cache_sha256'],ids=ids.tolist(),
            present_af_quantiles=np.quantile(af[:,-1],[0,.25,.5,.75,1]).tolist()) for entry,af,ids in selected],
        neutral_ids=neutral_ids.tolist(),neutral_present_status=dict(lost=int((final==0).sum()),fixed=int((final==1).sum()),
            segregating=int(((final>0)&(final<1)).sum())),
        tracking=null_meta['contract']['tracked_allele'],conditioning=null_meta['contract']['conditioning'],
        times_kya=dict(start=50,end=0,step=.025,count=2001),position_bp=500000,
        plot=dict(curves_per_group=100,selected_color=orange,neutral_color=gray,figsize_inches=[11,8.5],dpi=180,full_af_range=[0,1]),
        workers=args.workers,randomness='None; existing fixed-seed simulations, unchanged.',
        validation='Verified source hashes, identical 100 neutral IDs and selected cohort IDs, full time grid, neutral zero s and copy-count AF, selected present-day segregation, 100 line artists per group and vector PDF integrity.',
        outputs={path.name:sha(path) for path in outputs},script_sha256=sha(Path(__file__)),
        versions=dict(python=sys.version,numpy=np.__version__,matplotlib=matplotlib.__version__))
    save_json(out/'allele_frequency_provenance.json',provenance)
    save_json(out/'output_hashes.json',{p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file() and p.name!='output_hashes.json'})
    print(json.dumps(provenance['neutral_present_status']),flush=True)


if __name__=='__main__':main()
