"""Reproduce cumulative within-window power and neutral false-positive rates."""
import argparse
import csv
import json
import os
from pathlib import Path
import shutil

os.environ['MPLBACKEND']='Agg'
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
import report_cumulative_grid as base

OUT=base.ROOT/'docs/results/eas_strength_series_20261118'
RAW=Path.home()/'sim/eas_strength_series_20261118'
SIM=base.ROOT.parent/'sim'
WINDOWS=base.WINDOWS
WINDOW_NAMES=dict(zip(WINDOWS,base.LABELS))


def prepare(raw):
    manifest=raw/'manifest.json'
    contract=base.read(manifest)
    signature=base.sha(manifest)
    curves, sources={},[]
    archives={}
    for name,digest in contract['base_banks'].items():
        path=base.ROOT/'docs/results'/name/'mass.npz'
        assert base.sha(path)==digest==base.read(path.parent/'inputs.json')['cache_sha256']
        archives[name]=np.load(path)
    neutral=archives[base.OUT.name]
    for field in neutral.files:
        if field.endswith(('_neutral','_test','_neutral_ids','_test_ids')):
            curves[field]=neutral[field]
    for c in contract['cases']:
        k=c['key']
        if c['frozen']:
            bank=archives[c['frozen']]
            curves[k]=bank[k]
            ids=bank[k+'_ids']
            af=bank[k+'_af']
            assert len(set(ids))==100 and np.all((af>0)&(af<1))
            for slot,source in enumerate(ids):
                sources.append(dict(case=k,slot=slot,source=str(source),bank=c['frozen'],tag=None))
            curves[k+'_fixed']=np.zeros(100,dtype=bool)
        else:
            records=[base.read(raw/'slots'/k/f'{i:03d}.json') for i in range(100)]
            assert len({r['onset'] for r in records})==len(set(r['source'] for r in records))==100
            group=[]
            for i,r in enumerate(records):
                assert r['manifest']==signature
                path=Path(r['scores'])
                assert base.sha(path)==r['scores_sha256']
                with np.load(path) as d:
                    positions=d['positions']
                    ix,=np.flatnonzero(positions==r['position'])
                    np.testing.assert_array_equal(d['cutoffs'],np.arange(1000,50001,1000))
                    np.testing.assert_array_equal(d['valid_pairs'],1000)
                    group.append(d['cdfs'][2,ix])
                if r['tag'] is not None:
                    assert r['final_frequency']==1 and 0<r['tag']['modern_tag_af']<1 and r['tag']['position']==r['position']
                else:
                    assert 0<r['final_frequency']<1 and r['position']==500000
                sources.append(dict(case=k,slot=i,**r))
            curves[k]=np.array(group)
            curves[k+'_fixed']=np.array([r['tag'] is not None for r in records])
    for a in archives.values(): a.close()
    # Audit the complete held-out neutral proposal sets for fixed alleles.
    audit={}
    for age in ('10k','50k'):
        folder=SIM/'cumulative_grid_20261117'/f'heldout_de_novo_{age}'
        status=[base.read(p) for p in sorted(folder.glob('slot_*/status.json'))]
        assert len(status)==100 and all(v['accepted'] for v in status)
        audit['de_novo_'+age]=dict(fixed=sum(v['fixed'] for v in status),lost=sum(v['lost'] for v in status),tested=100)
    progress=base.read(SIM/'segregating_neutral_fpr_100_20261110/progress.json')
    audit['introgression']=dict(fixed=progress['fixed'],lost=progress['lost'],tested=100,
                                note='First 100 retained of the completed proposals')
    assert not any(x['fixed'] for x in audit.values())
    assert len(contract['cases'])==72 and len(sources)==7200
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'.gitattributes').write_text('* -text\n')
    shutil.copyfile(manifest, OUT/'simulation_manifest.json')
    np.savez_compressed(OUT/'mass.npz',**curves)
    base.save(OUT/'inputs.json',dict(cache_sha256=base.sha(OUT/'mass.npz'),records=sources,neutral_audit=audit))


def write_table(path, rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def rank_p(ref,obs):
    for a in (ref,obs):
        assert a.shape==(100,50) and np.isfinite(a).all() and np.diff(a,axis=1).min()>=-1e-12
        assert a.min()>=0 and a.max()<=1+1e-12
    return (1+(ref[:,None,:]>=obs[None,:,:]).sum(axis=0))/101


def report():
    contract=base.read(OUT/'simulation_manifest.json')
    assert base.sha(OUT/'mass.npz')==base.read(OUT/'inputs.json')['cache_sha256']
    rows, nulls, pvals=[],{},{}
    with np.load(OUT/'mass.npz') as data:
        for c in contract['cases']:
            k,origin,win=c['key'],c['origin'],c['regime']
            low,high=WINDOWS[win]
            name='introgression' if origin=='introgression' else ('de_novo_10k' if win=='10k' else 'de_novo_50k')
            assert set(data[name+'_test_ids']).isdisjoint(data[name+'_neutral_ids'])
            ref=data[name+'_neutral']
            p=rank_p(ref,data[k])
            pvals[k]=p
            if name+'_neutral_test' not in pvals:
                pvals[name+'_neutral_test']=rank_p(ref,data[name+'_test'])
            fixed=data[k+'_fixed']
            for a in (.05,.01):
                selected=(p[:,low-1:high]<a).any(axis=1)
                neutral=(pvals[name+'_neutral_test'][:,low-1:high]<a).any(axis=1)
                nulls[origin,win,a]=dict(origin=origin,window=win,alpha=a,fp=int(neutral.sum()),n=100,fpr=neutral.mean())
                rows.append(dict(origin=origin,window=win,s=c['selection'],alpha=a,tp=int(selected.sum()),n=100,tpr=selected.mean(),
                            tagged=int(fixed.sum()),tagged_calls=int(selected[fixed].sum()),segregating=int((~fixed).sum()),
                            segregating_calls=int(selected[~fixed].sum()),fp=int(neutral.sum()),neutral_n=100,fpr=neutral.mean()))
    assert len(rows)==144 and len(nulls)==12
    with (base.OUT/'rates.tsv').open() as f:
        for old in csv.DictReader(f,delimiter='\t'):
            if old['endpoint']=='FPR':
                assert nulls[old['origin'],old['window'],float(old['alpha'])]['fp']==int(old['detections'])
    write_table(OUT/'power.tsv',rows)
    write_table(OUT/'fpr.tsv',list(nulls.values()))
    np.savez_compressed(OUT/'pointwise_p.npz',**pvals)
    plot(rows)


def plot(rows):
    plt.rcParams.update({'font.size':14,'pdf.fonttype':42,'ps.fonttype':42})
    colors=('#145488','#a94621')
    strengths=sorted({r['s'] for r in rows})
    for a in (.05,.01):
        fig, axes=plt.subplots(1,2,figsize=(10.8,8.5),sharey=True)
        fig.subplots_adjust(top=.80,left=.10,right=.99,bottom=.10,wspace=.10)
        fig.suptitle(f'Within-window detection at pointwise p < {a:g}',fontsize=22,y=.975)
        for col,(ax,origin) in enumerate(zip(axes,base.ORIGINS)):
            subset={(r['s'],r['window']):r for r in rows if r['origin']==origin and r['alpha']==a}
            ax.set_title('De novo' if col==0 else '1% introgression',fontsize=17,pad=8)
            for y,s in enumerate([*strengths,0.]):
                for x,win in enumerate(WINDOWS):
                    row=subset[s or strengths[0],win]
                    pct=row['tp'] if s else row['fp']
                    cmap=plt.get_cmap('Blues' if s else 'Oranges')
                    ax.add_patch(Rectangle((x,y),1,1,facecolor=cmap(.08+.89*pct/100),edgecolor='white',linewidth=1.5))
                    ax.text(x+.5,y+.5,f'{pct}%',ha='center',va='center',fontsize=12.5,color='white' if pct>=60 else '#17212b')
            ax.set(xlim=(0,3),ylim=(13,0),xticks=np.arange(3)+.5,xticklabels=['10–0 kya','50–0 kya','50–30 kya'],
                   yticks=np.arange(13)+.5,yticklabels=[f'{s:.3f}' for s in strengths]+['FPR'])
            ax.xaxis.tick_top()
            ax.tick_params(length=0,pad=6)
            for spine in ax.spines.values():spine.set_visible(False)
        axes[0].set_ylabel('Selection coefficient')
        fig.text(.52,.045,'Each cell uses 100 replicates. Orange shows the held-out neutral FPR for each window.',ha='center',fontsize=13)
        save_figure(fig,f'heatmap_{a:g}')
    fig, axes=plt.subplots(2,3,figsize=(11,8.5),sharex=True,sharey=True)
    fig.subplots_adjust(top=.78,left=.09,right=.98,bottom=.13,wspace=.14,hspace=.26)
    fig.suptitle('Cumulative posterior-mass detection',fontsize=22,y=.975)
    for i,a in enumerate((.05,.01)):
        for j,win in enumerate(WINDOWS):
            ax=axes[i,j]
            for color,origin in zip(colors,base.ORIGINS):
                r=sorted((r for r in rows if (r['alpha'],r['window'],r['origin'])==(a,win,origin)),key=lambda r:r['s'])
                ax.plot(strengths,[x['tp'] for x in r],color=color,marker='o',markersize=4,linewidth=1.9)
                ax.axhline(r[0]['fp'],color=color,linestyle='--',linewidth=1.4)
            ax.grid(axis='y',alpha=.22)
            ax.spines[['top','right']].set_visible(False)
            ax.set(ylim=(-3,103),yticks=[0,25,50,75,100],xticks=[.001,.01,.02,.03],xticklabels=['.001','.01','.02','.03'])
            ax.tick_params(labelsize=11)
            if not i:ax.set_title(WINDOW_NAMES[win],fontsize=16)
            if not j:ax.set_ylabel(f'p < {a:g}\nCall rate (%)',fontsize=14)
    fig.supxlabel('Selection coefficient',fontsize=16,y=.04)
    legend=[Line2D([0],[0],color=c,linewidth=3) for c in colors]+[Line2D([0],[0],color='#333333',marker='o'),Line2D([0],[0],color='#333333',linestyle='--')]
    fig.legend(legend,['De novo','1% introgression','TPR','Neutral FPR'],loc='upper center',bbox_to_anchor=(.5,.92),ncol=4,frameon=False)
    save_figure(fig,'power_by_strength')


def save_figure(fig,name):
    for ext in ('png','pdf'):
        fig.savefig(OUT/f'{name}.{ext}',dpi=170)
    plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('phase',choices=('prepare','report'))
    p.add_argument('--raw',type=Path,default=RAW)
    args=p.parse_args()
    if args.phase=='prepare':prepare(args.raw.resolve())
    else:report()
