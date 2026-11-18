"""Build a portable focal-site cumulative-mass dataset and report window-matched detection."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path

os.environ['MPLBACKEND'] = 'Agg'
import numpy as np
import matplotlib.pyplot as plt
import shutil
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'docs/results/denovo_comparison_20261111'
OUT = ROOT/'docs/results/cumulative_grid_20261117'
WINDOWS = {'10k':(1,10),'50k':(1,50),'50_30k':(30,50)}
LABELS = ['10 kya–present','50 kya–present','50–30 kya']
ORIGINS = ('de_novo','introgression')
STRENGTHS = (.001,.005,.01)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def save(path,data):
    path.write_text(json.dumps(data,indent=2)+'\n')


def prepare():
    from run_cumulative_grid import RAW, ARCHIVE, archive_rows, local_path
    contract = read(RAW/'manifest.json')
    sha_manifest=sha(RAW/'manifest.json')
    old = archive_rows()
    portable = {}
    provenance = {}
    def store(key,curves,ids,af):
        assert len(set(ids)) == len(curves) == len(af) == 100, key
        assert np.all((af>0)&(af<1)), key
        portable[key] = curves
        portable[key+'_ids'] = np.asarray(ids,dtype=str)
        portable[key+'_af'] = af
    meta=read(BASE/'manifest.json')
    for key, entry in meta['cohorts'].items():
        if not (key.endswith(('_s0p005','_s0p01','neutral'))) or key.endswith('_50k_s0p01'):
            continue
        path=BASE/entry['cache']
        assert sha(path)==entry['cache_sha256']
        with np.load(path) as d: store(key,d['cdfs'][:,2],d['ids'].tolist(),d['af'][:,-1])
        provenance[key]=entry['cache_sha256']
    intro_file=ROOT/'docs/results/segregating_neutral_fpr_100_20261110/parameters.npz'
    with np.load(intro_file) as d:
        np.testing.assert_array_equal(d['reference_cdfs'][:,2],portable['introgression_neutral'])
        store('introgression_test',d['test_cdfs'][:,2],d['test_ids'].tolist(),d['test_af'][:,-1])
    provenance['introgression_test']=sha(intro_file)
    sources=[]
    for s in contract['strata']:
        archived=old.get((s['regime'],s['origin']),[]) if s['role']=='selected' else []
        assert len(archived)==s['reused']
        records=[]
        for row in archived:
            directory=local_path(row['source_directory'])
            assert sha(directory/'complete.json')==row['receipt_sha256']
            records.append((directory,read(directory/'complete.json')))
        for slot in range(s['reused'],100):
            record=read(RAW/'slots'/s['key']/f'{slot:03d}.json')
            directory=Path(record['source_directory'])
            assert record['config_sha256']==sha_manifest and record['manifest_sha256']==sha_manifest
            records.append((directory,record))
        curves=[]; ids=[]; af=[]
        for directory,record in records:
            assert record['origin']==s['origin'] and record['selection']==s['selection'] and record['stop_years_ago']==0
            assert record.get('selection_start_years',s['start'])==s['start'] and record.get('selection_end_years',s['end'])==s['end']
            assert sha(directory/'scores.npz')==record['scores_sha256'] and sha(directory/'trajectory.tsv')==record['trajectory_sha256']
            with np.load(directory/'scores.npz') as d:
                np.testing.assert_array_equal(d['positions'],np.arange(0,1000000,10000))
                np.testing.assert_array_equal(d['cutoffs'],np.arange(1000,50001,1000))
                np.testing.assert_array_equal(d['valid_pairs'],1000)
                curves.append(d['cdfs'][2,50])
            ids.append(str(directory)); af.append(record['final_frequency'])
            sources.append(dict(id=str(directory),group=s['key'],onset_sha256=record['onset_sha256'],scores_sha256=record['scores_sha256'],seed=record['seed']))
        key=f"{s['origin']}_{s['regime']}_s0p001" if s['role']=='selected' else f"de_novo_{s['regime']}_test"
        store(key,np.asarray(curves),ids,np.asarray(af))
    for age in ('10k','50k'):
        refs=meta['cohorts'][f'de_novo_{age}_neutral']['provenance']['records']
        held={r['onset_sha256'] for r in sources if r['group']==f'heldout_de_novo_{age}'}
        assert len(held)==100 and not held.intersection(r['onset_sha256'] for r in refs)
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'.gitattributes').write_text('* -text\n')
    shutil.copyfile(RAW/'manifest.json', OUT/'simulation_manifest.json')
    (OUT/'README.md').write_text('Cumulative posterior-weighted mass at 500 kb, over 1,000 fixed pairs. Each replicate contributes one detection if any 1-kya cutoff inside the closed selection window has strict pointwise p < alpha. The empirical rank is (1 + number of reference values at least as large) / 101. Selected, test, and reference focal alleles segregate in the present-day population. The independent 50-kya de novo neutral bank serves both windows beginning at 50 kya. The ongoing 50-kya s=0.01 groups are excluded for both origins because retained alleles were usually near fixation. Rates are conditional on segregation.\n')
    np.savez_compressed(OUT/'mass.npz',**portable)
    save(OUT/'inputs.json',dict(raw_manifest_sha256=sha_manifest,frozen_sources=provenance,
        archived_index_sha256=sha(ARCHIVE/'completed_units.tsv'),sources=sources,cache_sha256=sha(OUT/'mass.npz')))


def report():
    assert sha(OUT/'mass.npz')==read(OUT/'inputs.json')['cache_sha256']
    data=np.load(OUT/'mass.npz')
    rows=[]; detail={}
    for origin in ORIGINS:
        for regime,(low,high) in WINDOWS.items():
            age='10k' if regime=='10k' else '50k'
            ref_key='introgression_neutral' if origin=='introgression' else f'de_novo_{age}_neutral'
            test_key='introgression_test' if origin=='introgression' else f'de_novo_{age}_test'
            for s in (*STRENGTHS,0.):
                if regime=='50k' and s==.01: continue
                query=f"{origin}_{regime}_s{str(s).replace('.','p')}" if s else test_key
                ref, obs=data[ref_key], data[query]
                for arr in (ref,obs):
                    assert arr.shape==(100,50) and np.isfinite(arr).all() and np.diff(arr,axis=1).min()>=-1e-12
                assert set(data[ref_key+'_ids']).isdisjoint(data[query+'_ids'])
                p=(1+(ref[:,None,:]>=obs[None,:,:]).sum(axis=0))/101
                detail[f'{origin}_{regime}_{s}']=p
                for alpha in (.05,.01):
                    flags=(p[:,low-1:high]<alpha).any(axis=1)
                    rows.append(dict(origin=origin,window=regime,s=s,endpoint='TPR' if s else 'FPR',alpha=alpha,
                                     detections=int(flags.sum()),n=len(flags),rate=float(flags.mean())))
    np.savez_compressed(OUT/'pointwise_p.npz',**detail)
    with (OUT/'rates.tsv').open('w', newline='') as f:
        writer = csv.DictWriter(f, delimiter='\t', lineterminator='\n', fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    historical=read(ROOT/'docs/results/window_tpr_fpr_20261117/summary.json')
    for r in rows:
        if r['s'] in (.005,.01) or (r['origin']=='introgression' and not r['s']):
            hits=[x for x in historical[r['endpoint'].lower()] if x['regime']==r['window'] and x['origin']==r['origin']
                  and x.get('selection',0)==r['s'] and x['statistic']=='cumulative' and x['alpha']==r['alpha'] and x['method']=='Posterior-weighted mass']
            assert len(hits)==1 and r['detections']==hits[0].get('tp',hits[0].get('fp'))
    plt.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42,'font.size':15})
    for alpha in (.05,.01):
        fig,axes=plt.subplots(2,1,figsize=(8.5,11),layout='constrained')
        fig.suptitle(f'Cumulative posterior mass • pointwise p < {alpha:g}',fontsize=19)
        for ax,origin in zip(axes,ORIGINS):
            ax.set_title('De novo' if origin=='de_novo' else '1% introgression')
            for y,regime in enumerate(WINDOWS):
                for x,s in enumerate((*STRENGTHS,0.)):
                    row=next((r for r in rows if (r['origin'],r['window'],r['s'],r['alpha'])==(origin,regime,s,alpha)),None)
                    cmap=plt.get_cmap('Blues' if s else 'Oranges')
                    color=cmap(row['rate']) if row else '#eeeeee'
                    ax.add_patch(Rectangle((x,y),1,1,facecolor=color,edgecolor='white',linewidth=3))
                    ax.text(x+.5,y+.5,f"{row['detections']}%" if row else 'Excluded',ha='center',va='center',
                            color='white' if row and row['rate']>.6 else 'black')
            ax.set(xlim=(0,4),ylim=(3,0),xticks=np.arange(4)+.5,xticklabels=['s = .001','s = .005','s = .01','Neutral FPR'],
                   yticks=np.arange(3)+.5,yticklabels=LABELS)
            ax.tick_params(length=0)
            for spine in ax.spines.values(): spine.set_visible(False)
        fig.supxlabel('Each cell uses 100 segregating replicates at the focal site.',fontsize=13)
        for ext in ('png','pdf'): fig.savefig(OUT/f'cumulative_{alpha:g}.{ext}',dpi=160)
        plt.close(fig)
    print(json.dumps(rows),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase',choices=['prepare','report'])
    args=parser.parse_args()
    if args.phase=='prepare':
        prepare()
    else: report()
