#!/usr/bin/env python3
"""Exact all-pair truth audit for completed three-episode arrays; no new simulation."""
import os
os.environ['MPLBACKEND']='Agg'
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import numpy as np
import pandas as pd
import tskit
from audit_recent_all_pairs_truth import all_pairs_cdf
from report_selection_arrays import evaluate,REGIMES,local_path
from report_ongoing_selection_pilot import sha,rate
from report_time_resolved_pilot import empirical_p,interval_values,expand_support


def calculate(record):
    path=local_path(record['source_directory'])/'observed.trees'
    assert sha(path)==record['observed_sha256']
    ts=tskit.load(path);assert ts.num_samples==200
    return all_pairs_cdf(ts,np.arange(100)*10000,np.arange(1000,50001,1000))


def audit(root,out,workers):
    records=pd.read_csv(out/'completed_units.tsv',sep='\t')
    provenance=json.loads((out/'provenance.json').read_text());assert local_path(provenance['root']).resolve()==root.resolve()
    cache=out/'all_pairs_truth_cdfs.npz';old={}
    proof=out/'all_pairs_truth_provenance.json'
    algorithm_sha=sha(Path(__file__).with_name('audit_recent_all_pairs_truth.py'))
    if cache.exists() and proof.exists():
        previous=json.loads(proof.read_text());assert sha(cache)==previous['cdfs_sha256'],'Corrupt all-pair cache'
        if previous['algorithm_source_sha256']==algorithm_sha and previous['pair_count']==19900:
            with np.load(cache) as a:
                assert a['cdfs'].shape[1:]==(100,50) and np.isfinite(a['cdfs']).all()
                assert a['cdfs'].min()>=0 and a['cdfs'].max()<=1 and np.diff(a['cdfs'],axis=-1).min()>=0
                old={str(i):(str(h),v.copy()) for i,h,v in zip(a['ids'],a['observed_sha256'],a['cdfs'])}
    missing=[r for r in records.to_dict('records') if r['id'] not in old or old[r['id']][0]!=r['observed_sha256']]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for r,cdf in zip(missing,pool.map(calculate,missing)):old[r['id']]=(r['observed_sha256'],cdf)
    cdfs=np.stack([old[i][1] for i in records.id])
    np.savez_compressed(cache,ids=records.id.to_numpy(str),observed_sha256=records.observed_sha256.to_numpy(str),cdfs=cdfs)
    rows=[];summary=[]
    for scan,widths in [('one_kyr',[1]),('multiscale',[1,2,5,10])]:
        values,intervals=interval_values(cdfs,widths)
        for regime in REGIMES:
            cfg=provenance['manifest']['contract']['configs'][regime]
            for origin in ('de_novo','introgression'):
                pool=provenance['schedule']['null_pool_mapping'][regime][origin]
                mask=(records.regime==pool)&(records.origin==origin)
                dev=mask&(records.cohort=='development');cal=mask&(records.cohort=='calibration');test=mask&(records.cohort=='neutral_test')
                selected=(records.regime==regime)&(records.origin==origin)&(records.cohort=='selected')
                if dev.sum()<5 or cal.sum()<19 or test.sum()<1:continue
                mu=values[dev].mean(axis=(0,1));sd=values[dev].std(axis=(0,1),ddof=1);valid=sd>0
                z=np.full(values.shape,-np.inf);z[...,valid]=(values[...,valid]-mu[valid])/sd[valid]
                maxima=z[cal].max(axis=(1,2));tested=[]
                for i in np.flatnonzero(test|selected):
                    rec=records.iloc[i];calls=expand_support(empirical_p(z[i],maxima)<=.05,intervals)
                    metrics=evaluate(calls,np.arange(100)*10000,cfg['selection_start_years']/1000,cfg['selection_end_years']/1000,bool(selected.iloc[i]))
                    row=dict(scan=scan,regime=regime,origin=origin,selection=rec.selection,id=rec.id,cohort=rec.cohort,method='truth_all_pairs',**metrics)
                    rows.append(row);tested.append(row)
                frame=pd.DataFrame(tested);null=frame[frame.cohort=='neutral_test']
                for s in provenance['manifest']['contract']['selection']:
                    sel=frame[(frame.cohort=='selected')&(frame.selection==s)];tp=int(sel.joint_hit.sum());fp=int(null.any_call.sum());fn=len(sel)-tp
                    for metric,data in [('joint_tpr',sel.joint_hit),('spatial_tpr',sel.spatial_hit),('neutral_fpr',null.any_call),('near_mistimed',sel.near_mistimed),('near_off_target',sel.near_off_target)]:
                        summary.append(dict(scan=scan,regime=regime,origin=origin,selection=s,method='truth_all_pairs',metric=metric,
                            region_f1=2*tp/(2*tp+fp+fn) if len(sel) else None,
                            mean_temporal_f1=float(sel.temporal_f1.mean()) if len(sel) else None,**rate(data)))
    pd.DataFrame(rows).to_csv(out/'all_pairs_truth_replicates.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(summary).to_csv(out/'all_pairs_truth_performance.tsv',sep='\t',index=False,na_rep='NA')
    (out/'all_pairs_truth_provenance.json').write_text(json.dumps(dict(source_sha256=sha(__file__),pair_count=19900,
        algorithm_source_sha256=algorithm_sha,
        completed_units_sha256=sha(out/'completed_units.tsv'),report_provenance_sha256=sha(out/'provenance.json'),
        cdfs_sha256=sha(cache),n_regions=len(records),newly_computed=len(missing),workers=workers),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',required=True);p.add_argument('--out',required=True);p.add_argument('--workers',type=int,default=4);a=p.parse_args()
    audit(Path(a.root),Path(a.out),a.workers)
