#!/usr/bin/env python3
"""Prepare portable survivor caches, render updates, and audit their ranks."""
import os
os.environ['MPLBACKEND']='Agg'
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import PercentFormatter
import numpy as np
from PIL import Image
from pointwise_time_tests import exceedance_counts
from plot_cohort_significant_times import mass_order

ROOT=Path(__file__).resolve().parents[1]
PRIOR=ROOT/'docs/results/denovo_comparison_20261111'
LABELS={'10k':'10 kya–present','50k':'50 kya–present','50_30k':'50–30 kya'}
METHODS=['True TMRCA','Posterior median','Posterior-weighted mass']
METHOD_KEYS=['true_tmrca','posterior_median','posterior_mass']


def digest(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path): return json.loads(Path(path).read_text())
def save(path,value): Path(path).write_text(json.dumps(value,indent=2)+'\n')
def npz(path):
    with np.load(path,allow_pickle=False) as f: return {k:f[k].copy() for k in f.files}


def prepare(args):
    raw,out=args.raw.resolve(),args.out.resolve(); out.mkdir(parents=True,exist_ok=True)
    simulation=read(raw/'manifest.json'); prior=read(PRIOR/'manifest.json')
    if not args.preview:
        status=read(raw/'finished.json')
        assert status['manifest_sha256']==digest(raw/'manifest.json')
    else: status=read(raw/'progress.json')
    nulls={}
    for origin,key in [('de_novo','de_novo_50k_neutral'),('introgression','introgression_neutral')]:
        entry=prior['cohorts'][key]; path=PRIOR/entry['cache']
        assert digest(path)==entry['cache_sha256']
        d=npz(path); assert len(d['ids'])==100
        destination=out/(origin+'_neutral.npz'); shutil.copy2(path,destination)
        nulls[origin]=dict(cache=destination.name,cache_sha256=digest(destination),source=str(path))
    cases=[]
    for s in simulation['strata']:
        rows=s['reused'] if args.preview else read(raw/(s['key']+'_cohort.json'))['accepted']
        af=[]; cdfs=[]; ids=[]
        for r in rows:
            directory=Path(r['source_directory'])
            for name,field in [('scores.npz','scores_sha256'),('trajectory.tsv','trajectory_sha256'),
                               ('observed.trees','observed_sha256'),('posterior.zst','posterior_sha256')]:
                assert digest(directory/name)==r[field]
            assert 0<r['final_frequency']<1 and r['origin']==s['origin'] and r['selection']==s['selection']
            trajectory=np.loadtxt(directory/'trajectory.tsv',skiprows=1)
            np.testing.assert_array_equal(trajectory[:,0],np.arange(50000,-1,-25))
            np.testing.assert_allclose(trajectory[:,1],trajectory[:,4]/(2*trajectory[:,2]),atol=1e-14)
            af.append(trajectory[:,1]); ids.append(str(directory))
            scores=npz(directory/'scores.npz')
            np.testing.assert_array_equal(scores['cutoffs'],np.arange(1000,50001,1000))
            np.testing.assert_array_equal(scores['positions'],np.arange(0,1000000,10000))
            np.testing.assert_array_equal(scores['valid_pairs'],1000)
            cdfs.append(scores['cdfs'][:,50])
        n=len(rows)
        data=dict(ids=np.array(ids,dtype=str),proposal_indices=np.arange(n),af=np.stack(af) if n else np.empty((0,2001)),
            af_times_kya=np.arange(50000,-1,-25)/1000,cdfs=np.stack(cdfs) if n else np.empty((0,3,50)))
        neutral=npz(out/nulls[s['origin']]['cache'])
        counts=exceedance_counts(data['cdfs'],neutral['cdfs'])
        data.update(p_cumulative=(counts[0]+1)/101,p_interval=(counts[1]+1)/101)
        key=s['key'].replace('.','p'); cache=out/(key+'.npz'); np.savez_compressed(cache,**data)
        row=next(r for r in status['strata'] if r['key']==s['key'])
        cases.append(dict(key=key,origin=s['origin'],regime=s['regime'],selection=s['selection'],
            window=[s['end']/1000,s['start']/1000],cache=cache.name,cache_sha256=digest(cache),
            n=n,reused_n=len(s['reused']),records=rows,progress=row))
    save(out/'manifest.json',dict(cases=cases,nulls=nulls,simulation_manifest_sha256=digest(raw/'manifest.json'),
        preview=args.preview,position_bp=500000,pairs=1000,cutoffs_kya=list(range(1,51)),scan_stride_bp=10000,
        conditioning=simulation['conditioning'],p_definition='(1 + number of null values >= query) / 101; strict pointwise p < alpha.',
        interpretation='Focal-site any-time pointwise detection; no multiple-testing correction or independent new FPR cohort.',
        script_sha256=digest(__file__)))
    shutil.copy2(raw/'manifest.json',out/'simulation_manifest.json')
    shutil.copy2(raw/'schedule.json',out/'simulation_schedule.json')
    save(out/'simulation_summary.json',status)
    (out/'.gitattributes').write_text('* -text\n')
    (out/'.gitignore').write_text('# Figures are regenerated from the portable NPZ caches.\n*.png\n*.pdf\n')
    print(json.dumps([{k:c[k] for k in ('key','n','reused_n')} for c in cases]),flush=True)


def render(args):
    out=args.out.resolve(); meta=read(out/'manifest.json')
    data={}
    for c in meta['cases']:
        assert digest(out/c['cache'])==c['cache_sha256']; data[c['key']]=npz(out/c['cache'])
    nulls={}
    for origin,entry in meta['nulls'].items():
        assert digest(out/entry['cache'])==entry['cache_sha256']; nulls[origin]=npz(out/entry['cache'])
    matplotlib.rcParams.update({'font.size':15,'axes.labelsize':15,'axes.titlesize':16,
        'xtick.labelsize':13,'ytick.labelsize':13,'pdf.fonttype':42,'ps.fonttype':42})
    def clean(ax):
        ax.grid(alpha=.15); ax.spines[['top','right']].set_visible(False)
    def title(c):
        return ('De novo' if c['origin']=='de_novo' else '1% introgression')+' · '+LABELS[c['regime']]+f" · s = {c['selection']:g}"
    def emit(fig,stem,book):
        fig.savefig(out/(stem+'.png'),dpi=150)
        fig.savefig(out/(stem+'.pdf'))
        book.savefig(fig); plt.close(fig)
    rows=[]; orders={}
    for c in meta['cases']:
        d=data[c['key']]; n=c['n']
        for method,name in enumerate(METHODS):
            for kind in ('cumulative','interval'):
                for alpha in (.05,.01):
                    mask=d['p_'+kind][:,method]<alpha
                    rows.append(dict(case=c['key'],method=name,statistic=kind,alpha=alpha,selected_n=n,neutral_n=100,
                        any_time_count=int(mask.any(1).sum()) if n else None,any_time_fraction=float(mask.any(1).mean()) if n else None,
                        median_final_af=float(np.median(d['af'][:,-1])) if n else None))
    positive=[c for c in meta['cases'] if c['n']]
    with PdfPages(out/'allele_frequency.pdf') as book:
        for c in positive:
            d=data[c['key']]; neutral=nulls[c['origin']]; color='#0072B2' if c['origin']=='de_novo' else '#D55E00'
            fig,ax=plt.subplots(figsize=(11,8.5)); fig.subplots_adjust(top=.74,bottom=.17,left=.11,right=.97)
            fig.suptitle(title(c),y=.985,fontsize=20)
            fig.text(.5,.92,f"{c['n']} selected survivors ({c['reused_n']} previous + {c['n']-c['reused_n']} added)",ha='center',fontsize=15)
            fig.legend(handles=[Line2D([],[],color='.5',label='100 neutral survivors'),
                Line2D([],[],color=color,label='Selected survivors'),Patch(color='.92',label='Selection active')],
                loc='upper center',bbox_to_anchor=(.5,.86),ncol=3,frameon=False,fontsize=13)
            ax.axvspan(*c['window'],color='.92')
            ax.plot(neutral['af_times_kya'],neutral['af'].T,color='.45',alpha=.35,lw=.7)
            ax.plot(d['af_times_kya'],d['af'].T,color=color,alpha=.55,lw=1)
            ax.set(xlim=(50,0),ylim=(0,1),xlabel='Time before present (kya)',ylabel='Focal allele frequency')
            ax.yaxis.set_major_formatter(PercentFormatter(1)); clean(ax)
            fig.text(.5,.045,'EAS demography. All original focal alleles remain segregating in the present-day population.',ha='center',fontsize=11)
            emit(fig,c['key']+'_af',book)
    for alpha,tag in ((.05,'p005'),(.01,'p001')):
        with PdfPages(out/('coalescence_and_calls_'+tag+'.pdf')) as book:
            for c in positive:
                d=data[c['key']]; neutral=nulls[c['origin']]; color='#0072B2' if c['origin']=='de_novo' else '#D55E00'
                for method,name in enumerate(METHODS):
                    fig,axes=plt.subplots(2,2,figsize=(11,8.5))
                    fig.subplots_adjust(top=.75,bottom=.17,left=.12,right=.97,wspace=.28,hspace=.65)
                    fig.suptitle(title(c),y=.985,fontsize=20)
                    fig.text(.5,.93,f'{name} · focal site 500 kb · pointwise p < {alpha:g}',ha='center',fontsize=14)
                    fig.legend(handles=[Line2D([],[],color='.5',label='Each neutral'),Line2D([],[],color=color,label='Each selected'),
                        Line2D([],[],color='black',ls='--',label=f'Null p < {alpha:g} boundary')],loc='upper center',
                        bbox_to_anchor=(.5,.87),ncol=3,frameon=False,fontsize=13)
                    order=mass_order(d) if method==2 else np.arange(c['n'])
                    orders[c['key']+':'+str(method)]=d['ids'][order].tolist()
                    for col,kind in enumerate(('cumulative','interval')):
                        x=np.arange(1,51) if col==0 else np.arange(1.5,50)
                        a=d['cdfs'][:,method]; b=neutral['cdfs'][:,method]
                        if col: a,b=np.diff(a),np.diff(b)
                        ax=axes[0,col]; ax.axvspan(*c['window'],color='.92')
                        ax.plot(x,b.T,color='.5',alpha=.2,lw=.6); ax.plot(x,a.T,color=color,alpha=.5,lw=.8)
                        ax.plot(x,np.sort(b,axis=0)[-5 if alpha==.05 else -1],color='black',ls='--',lw=2)
                        ax.set(title='Cumulative cutoff' if col==0 else '1-kya interval',xlim=(1,50),ylim=(0,None),ylabel='Pair fraction')
                        ax.yaxis.set_major_formatter(PercentFormatter(1)); clean(ax)
                        mask=d['p_'+kind][:,method]<alpha
                        ax=axes[1,col]; ax.axvspan(*c['window'],color='.92')
                        for display,index in enumerate(order):
                            bins=np.flatnonzero(mask[index])
                            if col: ax.broken_barh([(int(v)+1,1) for v in bins],(display+.65,.7),facecolors=color)
                            else: ax.scatter(bins+1,np.full(len(bins),display+1),s=2,c=color)
                        ax.set(xlim=(0,50),ylim=(c['n']+.8,.2),xlabel='Time before present (kya)',ylabel='Simulation (sorted)' if method==2 else 'Simulation')
                        ax.set_yticks(np.unique(np.linspace(1,c['n'],min(c['n'],4)).astype(int))); clean(ax)
                        ax.text(.5,-.42,f"{mask.any(1).sum()}/{c['n']} with ≥1 pointwise call",ha='center',transform=ax.transAxes,fontsize=11)
                    fig.text(.5,.035,'100 frozen survivor neutrals. Mass rows sorted at p<0.05; pointwise calls do not estimate selection duration.',ha='center',fontsize=10.5)
                    emit(fig,c['key']+'_'+METHOD_KEYS[method]+'_'+tag,book)
    save(out/'summary.json',dict(results=rows,interpretation=meta['interpretation']))
    save(out/'row_orders.json',orders)
    save(out/'output_hashes.json',{p.name:digest(p) for p in sorted(out.iterdir()) if p.is_file() and p.name!='output_hashes.json' and p.suffix!='.log'})
    print(json.dumps({'rendered_cohorts':len(positive),'numeric_rows':len(rows)}))


def check(args):
    out=args.out.resolve(); meta=read(out/'manifest.json'); rows=read(out/'summary.json')['results']
    for name,value in read(out/'output_hashes.json').items(): assert digest(out/name)==value,name
    assert len(rows)==len(meta['cases'])*12
    for c in meta['cases']:
        assert digest(out/c['cache'])==c['cache_sha256']
        d=npz(out/c['cache']); entry=meta['nulls'][c['origin']]
        assert digest(out/entry['cache'])==entry['cache_sha256']; ref=npz(out/entry['cache'])
        assert len(d['ids'])==c['n'] and len(set(d['ids']))==c['n']
        assert set(d['ids']).isdisjoint(ref['ids']) and len(ref['ids'])==100
        assert d['af'].shape==(c['n'],2001) and d['cdfs'].shape==(c['n'],3,50)
        for x in (d,ref):
            assert np.all((x['af'][:,-1]>0)&(x['af'][:,-1]<1))
            assert np.isfinite(x['cdfs']).all() and np.all(np.diff(x['cdfs'])>=-1e-12)
        for method,name in enumerate(METHODS):
            x,y=d['cdfs'][:,method],ref['cdfs'][:,method]
            if method<2: x,y=np.rint(x*1000).astype(int),np.rint(y*1000).astype(int)
            for kind in ('cumulative','interval'):
                xx,yy=(x,y) if kind=='cumulative' else (np.diff(x),np.diff(y))
                expected=np.array([(1+(yy>=r).sum(0))/101 for r in xx]).reshape(c['n'],xx.shape[-1])
                np.testing.assert_array_equal(d['p_'+kind][:,method],expected)
                for alpha in (.05,.01):
                    r=next(r for r in rows if r['case']==c['key'] and r['method']==name and r['statistic']==kind and r['alpha']==alpha)
                    assert r['any_time_count']==(int((expected<alpha).any(1).sum()) if c['n'] else None)
    positive=sum(c['n']>0 for c in meta['cases'])
    pngs,pdfs=list(out.glob('*.png')),list(out.glob('*.pdf'))
    assert len(pngs)==7*positive and len(pdfs)==7*positive+3
    pages=0
    for path in pngs:
        with Image.open(path) as im: assert im.size==(1650,1275); im.verify()
    for path in pdfs:
        content=path.read_bytes()
        assert content.startswith(b'%PDF-') and content.rstrip().endswith(b'%%EOF')
        assert b'/FontFile2' in content and not re.search(rb'/Subtype\s*/Image',content)
        pages+=len(re.findall(rb'/Type\s*/Page\b',content))
    assert pages==14*positive
    print(json.dumps(dict(validated=True,cohorts=len(meta['cases']),numeric_rows=len(rows),pngs=len(pngs),vector_pdfs=len(pdfs),pdf_pages=pages,
        checks='Input/output hashes, unique simulations, segregation, independent tie-aware ranks, summary counts, image integrity, vector PDFs and editable text.'),indent=2))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('phase',choices=['prepare','report','check','all'])
    p.add_argument('--raw',type=Path,default=ROOT.parent/'sim/selected_survivor_extension_20261112')
    p.add_argument('--out',type=Path,default=ROOT/'docs/results/selected_survivor_extension_20261112')
    p.add_argument('--preview',action='store_true',help='Prepare only frozen pre-extension survivors for layout QA, not scientific results.')
    args=p.parse_args()
    if args.phase in ('prepare','all'): prepare(args)
    if args.phase in ('report','all'): render(args)
    if args.phase in ('check','all'): check(args)


if __name__=='__main__': main()
