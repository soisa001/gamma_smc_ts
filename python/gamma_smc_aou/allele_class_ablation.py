"""Fixed ALT/ALT, REF/REF and no-ALT/REF-penalty region scans at onset 50 kya."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import time

os.environ['MPLBACKEND']='Agg'
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[name]='1'
import numpy as np
import pandas as pd

from .fresh_power import atomic_json, digest
from .joint_scan_evaluation import load_one, metrics, rank_p

SCHEMES=('grid10kb','archaic_sites')
SOURCES=('truth','decoded')
FAMILIES=('aa','rr','contrast','weighted_contrast','joint','mass')


def methods(cutoffs):
    result=[dict(name='af_r1',family='af',gate=0,cutoff=-1,years=0)]
    for j,t in enumerate(cutoffs):
        result.append(dict(name=f'all_r1_T{t}',family='all',gate=0,cutoff=j,years=t))
        for family in FAMILIES:
            for gate in (0,50):
                result.append(dict(name=f'{family}_g{gate}_r1_T{t}',family=family,gate=gate,cutoff=j,years=t))
    return result


def local_scores(counts,n,valid):
    fractions=np.divide(counts,n[None,:,:,None],out=np.zeros_like(counts,dtype=float),where=n[None,:,:,None]>0)
    rr,ar,aa,all_pairs=[fractions[:,:,j] for j in range(4)]
    aa=aa.copy();rr=rr.copy()
    aa[:,~valid]=0;rr[:,~valid]=0
    contrast=np.maximum(aa-rr,0)
    weight=(n[:,2]/n[:,3])[None,:,None]
    weighted=weight*contrast
    return dict(all=all_pairs,aa=aa,rr=rr,contrast=contrast,
                weighted_contrast=weighted,joint=weighted*(1-ar),mass=weight*aa)


def score_one(payload):
    item,features,gate_values,method_list,cutoffs=payload
    result=np.zeros((5,2,2,len(method_list)))
    diagnostics=[]
    for scheme_index,scheme in enumerate(SCHEMES):
        prefix='grid' if scheme=='grid10kb' else 'site'
        n=features[prefix+'_n'];c=features[prefix+'_counts']
        valid=features['grid_markers']>=0 if prefix=='grid' else np.ones(len(n),dtype=bool)
        local=local_scores(c,n,valid)
        if prefix=='site':
            af=features['site_af']
            assert not np.any(c[0,:,1]),'Exact-site true mixed pairs must predate the mutation'
            np.testing.assert_array_equal(local['joint'][0],local['weighted_contrast'][0])
        else:
            af=np.zeros(len(n));idx=features['grid_markers']
            af[valid]=features['site_af'][idx[valid]]
        cached={}
        for family,values in local.items():
            cached[(family,0)]=np.max(values,axis=1,initial=0)
        for source_index,source in enumerate(SOURCES):
            for j,t in enumerate(cutoffs):
                for family,klass in (('aa',2),('rr',0)):
                    values=local[family][source_index,:,j]
                    maximum=float(np.max(values,initial=0))
                    support=valid & (n[:,klass]>0) & (values==maximum)
                    diagnostics.append(dict(key=item['key'],onset_years=item['onset_years'],
                        scheme=scheme,source=source,years=t,family=family,maximum=maximum,
                        saturated_sites=int(np.sum(valid & (n[:,klass]>0) & (values==1))),
                        minimum_pairs_at_max=int(np.min(n[support,klass])) if support.any() else 0,
                        eligible_sites=int(valid.sum())))
        for fold in range(5):
            for family in FAMILIES:
                values=np.where(local['all']>gate_values[fold,:,None,:],local[family],0)
                cached[(family,50)]=np.max(values,axis=1,initial=0)
            for k,m in enumerate(method_list):
                result[fold,scheme_index,:,k]=float(np.max(af,initial=0)) if m['family']=='af' else cached[(m['family'],m['gate'])][:,m['cutoff']]
    return result,diagnostics


def calibrate(scores,items,folds,method_list):
    onset=np.array([i['onset_years'] for i in items]);neutral=onset==0
    predictions=[];targets=[]
    for fold in range(5):
        test=folds==fold
        calibration=neutral & np.isin(folds,[(fold+3)%5,(fold+4)%5])
        training=(~neutral)&(~test)
        if (int(calibration.sum()),int((neutral&test).sum()),int(training.sum()),int((~neutral&test).sum()))!=(400,200,80,20):
            raise ValueError('Unexpected fold sizes')
        for si,scheme in enumerate(SCHEMES):
            for src,source in enumerate(SOURCES):
                matrix=scores[fold,si,src]
                for k,m in enumerate(method_list):
                    ids=np.flatnonzero(test)
                    pvalues=rank_p(matrix[calibration,k],matrix[test,k])
                    for alpha in (.01,.05,.1):
                        for i,pv in zip(ids,pvalues):
                            predictions.append(dict(key=items[i]['key'],onset_years=int(onset[i]),fold=fold,
                                scheme=scheme,source=source,method=m['name'],alpha=alpha,
                                score=float(matrix[i,k]),p=float(pv),called=bool(pv<=alpha)))
                    if m['years'] in (0,50000):
                        threshold=float(np.quantile(matrix[training,k],.3,method='lower'))
                        for i in ids:
                            targets.append(dict(key=items[i]['key'],onset_years=int(onset[i]),fold=fold,
                                scheme=scheme,source=source,method=m['name'],alpha=.7,
                                threshold=threshold,score=float(matrix[i,k]),
                                called=bool(matrix[i,k]>=threshold and matrix[i,k]>0)))
    return pd.DataFrame(predictions),pd.DataFrame(targets)


def figures(table,out):
    import matplotlib
    matplotlib.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42,'font.size':17,'axes.labelsize':18,'axes.titlesize':18,'legend.fontsize':14})
    import matplotlib.pyplot as plt
    checks=[]
    for scheme in SCHEMES:
        fig,axes=plt.subplots(1,2,figsize=(11,8.5),sharey=True)
        fig.subplots_adjust(left=.1,right=.975,bottom=.15,top=.74,wspace=.14)
        names=['all_r1_T50000','aa_g0_r1_T50000','rr_g0_r1_T50000','contrast_g0_r1_T50000',
               'weighted_contrast_g50_r1_T50000','joint_g50_r1_T50000']
        labels=['All pairs','ALT/ALT','REF/REF','AA minus RR','Weighted, no AR penalty','Original joint']
        colors=['#777777','#228833','#AA3377','#EEAA33','#0077BB','#EE6677']
        for ax,source in zip(axes,SOURCES):
            for method,label,color in zip(names,labels,colors):
                sub=table[(table.scheme==scheme)&(table.source==source)&(table.method==method)].sort_values('neutral_call_fraction')
                ax.plot(sub.neutral_call_fraction*100,sub.power*100,'o-',label=label,color=color,lw=2)
            ax.set(xlim=(-.5,12),ylim=(-2,102),xlabel='Neutral regions called (%)',title=source.title())
            ax.axhline(70,ls=':',color='#999999',lw=1);ax.grid(alpha=.15)
        axes[0].set_ylabel('Selected regions called (%)')
        handles,labels=axes[0].get_legend_handles_labels()
        legend=fig.legend(handles,labels,ncol=3,loc='upper center',bbox_to_anchor=(.52,.92),frameon=False)
        fig.suptitle('50-kya onset: '+('10 kb grid' if scheme=='grid10kb' else 'every archaic SNP'),fontsize=21,y=.98)
        fig.canvas.draw();renderer=fig.canvas.get_renderer();box=legend.get_window_extent(renderer)
        if any(box.overlaps(ax.get_tightbbox(renderer)) for ax in axes):raise ValueError('Legend overlaps panel or title')
        for box in [legend.get_window_extent(renderer)]+[ax.get_tightbbox(renderer) for ax in axes]:
            if box.x0<0 or box.y0<0 or box.x1>fig.bbox.x1+1 or box.y1>fig.bbox.y1+1:raise ValueError('Figure content outside canvas')
        name=f'allele_class_{scheme}'
        fig.savefig(out/f'{name}.png',dpi=180)
        fig.savefig(out/f'{name}.pdf',metadata={'CreationDate':None,'ModDate':None})
        if b'/Subtype /Image' in (out/f'{name}.pdf').read_bytes():raise ValueError('Rasterized PDF')
        checks.append(dict(name=name,legend_clear=True,labels_within_canvas=True,vector_pdf=True));plt.close(fig)
    atomic_json(out/'figure_checks.json',checks)


def evaluate(input_dir,out,workers):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=True)
    old=input_dir/'analysis'
    old_artifacts=json.loads((old/'artifact_manifest.json').read_text())
    for name in ('region_inventory.csv','gate_fits.json','methods.json','fold_region_scores.npz'):
        if digest(old/name)!=old_artifacts[name]['sha256']:raise ValueError(f'Corrupt reference {name}')
    manifest=json.loads((input_dir/'profile_manifest.json').read_text())
    cfg=manifest['config'];items=[i for i in manifest['items'] if i['onset_years'] in (0,50000)]
    if len(items)!=1100 or sum(i['onset_years']==50000 for i in items)!=100:raise ValueError('Expected 1000 neutral and 100 immediate-onset regions')
    identity=dict(source=digest(Path(__file__)),profile_manifest=digest(input_dir/'profile_manifest.json'),reference_manifest=digest(old/'artifact_manifest.json'))
    if (out/'artifact_manifest.json').exists():
        previous=json.loads((out/'provenance.json').read_text())
        if previous['identity']!=identity:raise ValueError('Analysis identity changed; use a new output directory')
        for name,spec in json.loads((out/'artifact_manifest.json').read_text()).items():
            if (out/name).stat().st_size!=spec['bytes'] or digest(out/name)!=spec['sha256']:raise ValueError(f'Corrupt cached output {name}')
        print('Verified cached analysis',flush=True);return
    inventory=pd.read_csv(old/'region_inventory.csv').set_index('key')
    folds=np.array([inventory.loc[i['key'],'fold'] for i in items],dtype=int)
    with ThreadPoolExecutor(max_workers=workers) as pool:loaded=list(pool.map(load_one,[(i,input_dir) for i in items]))
    print(f'Loaded and hash-verified {len(loaded)} profiles',flush=True)
    features=[x[0] for x in loaded]
    all_grid=np.array([f['grid_counts'][:,:,3,:]/f['grid_n'][None,:,3,None] for f in features])
    neutral=np.array([i['onset_years']==0 for i in items]);gates=[];fold_records=[]
    old_gates=json.loads((old/'gate_fits.json').read_text())
    for fold in range(5):
        fit=neutral&np.isin(folds,[(fold+1)%5,(fold+2)%5]);cal=neutral&np.isin(folds,[(fold+3)%5,(fold+4)%5])
        assert fit.sum()==cal.sum()==400
        gate=np.quantile(all_grid[fit],.5,axis=(0,2));np.testing.assert_array_equal(gate,old_gates[fold]['values']['0.5'])
        gates.append(gate)
        fold_records.append(dict(fold=fold,fit_ids=[items[i]['key'] for i in np.flatnonzero(fit)],
            calibration_ids=[items[i]['key'] for i in np.flatnonzero(cal)],median=gate.tolist()))
    method_list=methods(cfg['tmrca_cutoffs_years'])
    with ThreadPoolExecutor(max_workers=workers) as pool:
        scored=list(pool.map(score_one,[(i,f,np.array(gates),method_list,cfg['tmrca_cutoffs_years']) for i,f in zip(items,features)]))
    scores=np.stack([x[0] for x in scored],axis=3)
    print(f'Calculated fixed region scores: {scores.shape}',flush=True)
    # Verify every unchanged fixed rule against the original full scan.
    old_names={m['name']:k for k,m in enumerate(json.loads((old/'methods.json').read_text()))}
    original_index={i['key']:k for k,i in enumerate(manifest['items'])};ids=[original_index[i['key']] for i in items]
    comparisons=0
    with np.load(old/'fold_region_scores.npz',allow_pickle=False) as saved:
        for fold in range(5):
            for si,scheme in enumerate(SCHEMES):
                for src,source in enumerate(SOURCES):
                    reference=saved[f'fold{fold}_{scheme}_{source}'][ids]
                    for k,m in enumerate(method_list):
                        if m['name'] in old_names:
                            np.testing.assert_allclose(scores[fold,si,src,:,k],reference[:,old_names[m['name']]],rtol=0,atol=1e-14)
                            comparisons+=len(items)
    print(f'Verified {comparisons} unchanged reference score values',flush=True)
    pred,target=calibrate(scores,items,folds,method_list)
    compression=dict(method='gzip',mtime=0)
    pred.to_csv(out/'predictions.csv.gz',index=False,compression=compression)
    target.to_csv(out/'target70_predictions.csv.gz',index=False,compression=compression)
    table=metrics(pred);target_table=metrics(target)
    table.to_csv(out/'metrics.csv',index=False);target_table.to_csv(out/'target70_metrics.csv',index=False)
    diagnostics=pd.DataFrame([r for x in scored for r in x[1]])
    diagnostics.to_csv(out/'class_maxima.csv.gz',index=False,compression=compression)
    saturation=diagnostics.assign(saturated=diagnostics.maximum==1,empty=diagnostics.eligible_sites==0)
    saturation.groupby(['scheme','source','years','family','onset_years']).agg(regions=('key','size'),
        saturated_regions=('saturated','sum'),marker_free_regions=('empty','sum'),median_minimum_pairs_at_max=('minimum_pairs_at_max','median')).reset_index().to_csv(out/'saturation.csv',index=False)
    np.savez_compressed(out/'scores.npz',scores=scores)
    pd.DataFrame([dict(key=i['key'],onset_years=i['onset_years'],fold=int(folds[k])) for k,i in enumerate(items)]).to_csv(out/'regions.csv',index=False)
    atomic_json(out/'methods.json',method_list);atomic_json(out/'folds.json',fold_records)
    atomic_json(out/'provenance.json',dict(identity=identity,input_dir=str(input_dir),config=cfg,workers=workers,
        profile_hashes={i['key']:loaded[k][2] for k,i in enumerate(items)},regions=1100,selected_regions=100,neutral_regions=1000,
        included_onsets=[0,50000],excluded_onset=10000,unchanged_score_values_verified=comparisons,
        statistic='frac_recent_T',cutoff_selection='none',frequency_or_pair_count_filter='none',
        validation_scope='Exploratory existing-cohort cross-validation; prior outcomes informed design',
        primary_comparison='Exact archaic sites; ungated AA, RR and AA-RR at T=50000',
        missing_class='No stand-alone evidence; RR comparison baseline zero at fixation',
        minimum_calibration_p=1/401,seconds=time.monotonic()-start))
    figures(table,out)
    names=['predictions.csv.gz','target70_predictions.csv.gz','metrics.csv','target70_metrics.csv','class_maxima.csv.gz',
           'saturation.csv','scores.npz','regions.csv','methods.json','folds.json','provenance.json','figure_checks.json']
    names += [f'allele_class_{scheme}.{ext}' for scheme in SCHEMES for ext in ('png','pdf')]
    atomic_json(out/'artifact_manifest.json',{name:dict(bytes=(out/name).stat().st_size,sha256=digest(out/name)) for name in names})
    print(table[(table.scheme=='archaic_sites')&(table.alpha==.05)&table.method.str.endswith('T50000')][['source','method','power','neutral_call_fraction']].to_string(index=False),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--workers',type=int,default=20);args=p.parse_args()
    evaluate(args.input,args.out,args.workers)


if __name__=='__main__':main()
