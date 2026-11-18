"""Reproduce the reported 93-replicate and 20-tag interim snapshot."""
import argparse
import csv
from pathlib import Path
import report_cumulative_grid as base
import numpy as np

OUT = base.ROOT/'docs/results/cumulative_snapshot_20261118'


def prepare():
    raw = Path.home()/'sim'
    weak = raw/'cumulative_grid_20261118'
    strong = raw/'fixed_tag_gaps_20261118'
    rows = []
    with np.load(base.OUT/'mass.npz') as bank:
        def add(kind, origin, window, s, curves):
            refkey = 'introgression' if origin=='introgression' else 'de_novo_'+('10k' if window=='10k' else '50k')
            reference = bank[refkey+'_neutral']
            left,right = base.WINDOWS[window]
            for population, sample in [('selected',np.stack(curves)),('neutral',bank[refkey+'_test'])]:
                p = (1+(reference[:,None,:]>=sample[None,:,:]).sum(axis=0))/101
                for a in (.05,.01):
                    rows.append(dict(cohort=kind,origin=origin,window=window,s=s,population=population,alpha=a,
                                     n=len(sample),called=int((p[:,left-1:right]<a).any(axis=1).sum())))
        for g in base.read(weak/'manifest.json')['strata']:
            curves = []
            for i in range(93):
                rec = base.read(weak/'slots'/g['key']/f'{i:03d}.json')
                with np.load(Path(rec['source_directory'])/'scores.npz') as z:
                    curves.append(z['cdfs'][2,50])
            add('segregating',g['origin'],g['regime'],g['selection'],curves)
        for g in base.read(strong/'manifest.json')['cases']:
            curves = []
            for i in range(20):
                folder = strong/g['key']/f'{i:03d}'
                rec = base.read(folder/'tag.json')
                with np.load(folder/'scores.npz') as z:
                    np.testing.assert_array_equal(z['positions'],[rec['tag']['position']])
                    curves.append(z['cdfs'][2,0])
            add('tag',g['origin'],g['regime'],g['s'],curves)
    OUT.mkdir(parents=True,exist_ok=True)
    with (OUT/'counts.tsv').open('w',newline='') as f:
        writer = csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def report():
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    with (OUT/'counts.tsv').open() as f:
        rows = list(csv.DictReader(f,delimiter='\t'))
    lookup = {(r['cohort'],r['origin'],r['window'],float(r['s']),r['population'],float(r['alpha'])):
              100*int(r['called'])/int(r['n']) for r in rows}
    windows = dict(zip(base.WINDOWS,base.LABELS))
    origins = dict(de_novo='De novo', introgression='Introgressed')
    plt.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42,'font.size':13})
    for kind in ('segregating','tag'):
        is_tag = kind=='tag'
        groups = sorted({(r['origin'],r['window'],float(r['s']) if is_tag else 0.) for r in rows if r['cohort']==kind},
                        key=lambda g:(g[0],list(base.WINDOWS).index(g[1]),g[2]))
        cols = ['Tag calls','Neutral FPR'] if is_tag else ['.002','.003','.004','FPR']
        labels = [f'{origins[o]}   {windows[w]}' + (f'   s = {s:g}' if is_tag else '') for o,w,s in groups]
        fig, axes = plt.subplots(1,2,figsize=(11,7.1),sharey=True)
        fig.subplots_adjust(left=.39 if is_tag else .28,right=.985,top=.77,bottom=.17,wspace=.12)
        fig.suptitle('Detection at surviving tags after fixation' if is_tag else 'Detection while the focal allele segregates',
                     x=.5,y=.97,fontsize=22)
        fig.text(.5,.895,'First 20 fixed histories per treatment' if is_tag else 'First 93 selected replicates per treatment',
                 ha='center',fontsize=16)
        for ax,alpha in zip(axes,(.05,.01)):
            ax.set_title(f'Pointwise p < {alpha:g}',fontsize=16,pad=12)
            for i,(origin,win,s) in enumerate(groups):
                selected = [s] if is_tag else [.002,.003,.004]
                for j,value in enumerate([*selected,selected[0]]):
                    neutral = j==len(selected)
                    rate = lookup[kind,origin,win,value,'neutral' if neutral else 'selected',alpha]
                    color = plt.get_cmap('Oranges' if neutral else 'Blues')(.10+.87*rate/100)
                    ax.add_patch(Rectangle((j,i),1,1,facecolor=color,edgecolor='white',linewidth=2.5))
                    label = f'{rate:.1f}'.removesuffix('.0')+'%'
                    ax.text(j+.5,i+.5,label,ha='center',va='center',color='white' if rate>=58 else '#152331',fontsize=14)
            ax.set(xlim=(0,len(cols)),ylim=(len(groups),0),xticks=np.arange(len(cols))+.5,xticklabels=cols,
                   yticks=np.arange(len(groups))+.5,yticklabels=labels)
            ax.xaxis.tick_top()
            ax.tick_params(length=0,pad=6)
            for spine in ax.spines.values(): spine.set_visible(False)
        note = ('Orange gives the FPR in 100 segregating neutrals. The FPR for tags of neutral fixations is unmeasured.'
                if is_tag else 'Blue gives TPR at each selection coefficient. Orange gives the window-matched FPR in 100 held-out neutrals.')
        fig.text(.5,.075,note,ha='center',fontsize=12)
        for ext in ('pdf','png'):
            fig.savefig(OUT/f'{kind}_detection.{ext}',dpi=180)
        plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase',choices=('prepare','report'))
    if parser.parse_args().phase=='prepare':
        prepare()
    else:
        report()
