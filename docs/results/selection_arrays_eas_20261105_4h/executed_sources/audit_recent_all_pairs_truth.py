#!/usr/bin/env python3
"""Audit pair subsampling using exact genealogies, without new simulations."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import numpy as np
import pandas as pd
import tskit
from report_recent_selection_pilot import evaluate
from report_time_resolved_pilot import empirical_p
from report_ongoing_selection_pilot import sha,rate


def all_pairs_cdf(ts,positions,cutoffs,generation_time=25):
    """Count pairs by their unique MRCA; include all unordered sample pairs."""
    n=ts.num_samples;denominator=n*(n-1)//2
    result=np.zeros((len(positions),len(cutoffs)))
    for i,position in enumerate(positions):
        tree=ts.at(position);times=[];weights=[]
        for u in tree.nodes():
            size=tree.num_samples(u);weight=size*(size-1)//2
            for child in tree.children(u):
                size=tree.num_samples(child);weight-=size*(size-1)//2
            if weight:
                times.append(tree.time(u)*generation_time);weights.append(weight)
        assert sum(weights)==denominator,'Incomplete genealogy'
        order=np.argsort(times);times=np.asarray(times)[order]
        cumulative=np.r_[0,np.cumsum(np.asarray(weights)[order])]
        result[i]=cumulative[np.searchsorted(times,cutoffs,side='left')]/denominator
    return result


def calculate(payload):
    root,record,positions,cutoffs=payload
    path=Path(root)/record['unit']/'observed.trees'
    assert sha(path)==record['observed_sha256']
    ts=tskit.load(path)
    assert ts.num_samples==200
    carrier_record=None
    if record['cohort']=='selected':
        import run_time_resolved_pilot as base
        forward=path.parent/'forward.trees'
        sampled=base._sample_diploids(tskit.load(forward),100,base.seed('recent-v1:'+record['unit']+':sample'))
        sampled.tables.nodes.assert_equals(ts.tables.nodes)
        sampled.tables.edges.assert_equals(ts.tables.edges)
        carriers=0;states=np.zeros(sampled.num_nodes,dtype=bool)
        for variant in sampled.variants():
            assert variant.site.position==ts.sequence_length//2
            assert np.all((variant.genotypes==0)|(variant.genotypes==1))
            carriers=int(np.sum(variant.genotypes==1))
            states[sampled.samples()]=variant.genotypes==1
        pairs=np.loadtxt(Path(root)/'pairs.tsv',dtype=int,ndmin=2)
        meta=json.loads((path.parent/'posterior.zst.meta').read_text())
        np.testing.assert_array_equal(pairs,meta['pairs'])
        ordered=base.ordered_nodes(sampled)
        sampled_carrier_pairs=int(np.sum(states[ordered[pairs[:,0]]]&states[ordered[pairs[:,1]]]))
        carrier_record=dict(unit=record['unit'],origin=record['origin'],selection=record['selection'],
            population_frequency=record['final_frequency'],sample_carriers=carriers,sample_haplotypes=200,
            carrier_pairs=carriers*(carriers-1)//2,all_pairs=19900,
            sampled_carrier_pairs=sampled_carrier_pairs,sampled_pairs=len(pairs),forward_sha256=sha(forward))
    return all_pairs_cdf(ts,positions,cutoffs),carrier_record


def audit(root,report,workers):
    records=pd.read_csv(report/'completed_units.tsv',sep='\t')
    provenance=json.loads((report/'provenance.json').read_text());cfg=provenance['manifest']['config']
    assert sha(root/'manifest.json')==records.config_sha256.iloc[0]
    positions=np.arange(0,cfg['length'],cfg['stride']);cutoffs=np.arange(1000,50001,1000)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        results=list(pool.map(calculate,[(str(root),r,positions,cutoffs) for r in records.to_dict('records')]))
    cdfs=np.stack([r[0] for r in results])
    pd.DataFrame([r[1] for r in results if r[1] is not None]).to_csv(report/'selected_sample_carriers.tsv',sep='\t',index=False)
    values=np.diff(cdfs,axis=-1);rows=[];summaries=[]
    for origin in ('de_novo','introgression'):
        same=records.origin==origin;dev=same&(records.cohort=='development');cal=same&(records.cohort=='calibration')
        assert dev.sum()>=5 and cal.sum()>=19
        mu=values[dev].mean(axis=(0,1));sd=values[dev].std(axis=(0,1),ddof=1);valid=sd>0
        z=np.full(values.shape,-np.inf);z[...,valid]=(values[...,valid]-mu[valid])/sd[valid]
        maxima=z[cal].max(axis=(1,2));tested=[]
        for i in np.flatnonzero(same&records.cohort.isin(['selected','neutral_test'])):
            rec=records.iloc[i];called=empirical_p(z[i],maxima)<=.05
            row=dict(unit=rec.unit,origin=origin,selection=rec.selection,cohort=rec.cohort,
                **evaluate(called,positions,cfg['selection_start_years']/1000,rec.cohort=='selected'))
            rows.append(row);tested.append(row)
        frame=pd.DataFrame(tested);null=frame[frame.cohort=='neutral_test']
        for s in provenance['schedule']['parameters']['selection']:
            selected=frame[(frame.cohort=='selected')&(frame.selection==s)]
            for metric,indicator in [('joint_tpr',selected.joint_hit),('spatial_tpr',selected.spatial_hit),('neutral_fpr',null.any_call)]:
                summaries.append(dict(origin=origin,selection=s,method='truth_all_19900_pairs',metric=metric,**rate(indicator)))
    pd.DataFrame(rows).to_csv(report/'all_pairs_truth_replicates.tsv',sep='\t',index=False)
    pd.DataFrame(summaries).to_csv(report/'all_pairs_truth_performance.tsv',sep='\t',index=False)
    np.savez_compressed(report/'all_pairs_truth_cdfs.npz',cdfs=cdfs,units=records.unit.to_numpy(dtype=str),positions=positions,cutoffs=cutoffs)
    (report/'all_pairs_truth_provenance.json').write_text(json.dumps(dict(source_sha256=sha(__file__),
        report_provenance_sha256=sha(report/'provenance.json'),completed_units_sha256=sha(report/'completed_units.tsv'),
        pairs=19900,workers=workers,method='Exact pair counts by MRCA; separate whole-region max calibration; no additional simulations'),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--report',type=Path,required=True);p.add_argument('--workers',type=int,default=4)
    a=p.parse_args();audit(a.root,a.report,a.workers)
