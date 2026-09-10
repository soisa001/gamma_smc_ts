"""Independently audit region predictions, calibration ranks and reported rates."""
import argparse
from collections import defaultdict
import csv
import gzip
import hashlib
import json
from pathlib import Path
import numpy as np


def sha256(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def audit(out):
    analysis=out/'analysis'
    artifact=json.loads((analysis/'artifact_manifest.json').read_text())
    for name,spec in artifact.items():
        path=analysis/name
        assert path.stat().st_size==spec['bytes'] and sha256(path)==spec['sha256'],name
    manifest=json.loads((out/'profile_manifest.json').read_text())
    items=manifest['items'];index={x['key']:i for i,x in enumerate(items)}
    assert len(index)==1200
    with (analysis/'region_inventory.csv').open() as stream:
        inventory={row['key']:row for row in csv.DictReader(stream)}
    assert set(inventory)==set(index)
    gates=json.loads((analysis/'gate_fits.json').read_text())
    for gate in gates:
        test={k for k,v in inventory.items() if int(v['fold'])==gate['fold']}
        fit,cal=set(gate['fit_ids']),set(gate['calibration_ids'])
        assert len(fit)==len(cal)==400 and len(test)==240
        assert not fit&cal and not fit&test and not cal&test
        assert all(inventory[k]['onset_years']=='0' for k in fit|cal)
    methods=json.loads((analysis/'methods.json').read_text())
    mi={x['name']:i for i,x in enumerate(methods)}
    with (analysis/'training_choices.csv').open() as stream:
        choices={(int(r['fold']),r['scheme'],r['source']):r['method'] for r in csv.DictReader(stream)}
    scores=np.load(analysis/'fold_region_scores.npz',allow_pickle=False)
    references={}
    for fold in range(5):
        cal=[index[k] for k in gates[fold]['calibration_ids']]
        for scheme in ('grid10kb','archaic_sites'):
            for source in ('truth','decoded'):
                references[(fold,scheme,source)]=np.sort(scores[f'fold{fold}_{scheme}_{source}'][cal],axis=0)
    counts=defaultdict(lambda:[0,0])
    seen=set();rows=0
    with gzip.open(analysis/'heldout_predictions.csv.gz','rt',newline='') as stream:
        for row in csv.DictReader(stream):
            fold=int(row['fold']);onset=int(row['onset_years']);pv=float(row['p']);value=float(row['score'])
            scheme,source,method=row['scheme'],row['source'],row['method']
            alpha=float(row['alpha']);called=row['called']=='True'
            assert int(inventory[row['key']]['fold'])==fold
            assert int(inventory[row['key']]['onset_years'])==onset
            resolved=choices[(fold,scheme,source)] if method=='tuned_joint' else method
            k=mi[resolved]
            expected_score=scores[f'fold{fold}_{scheme}_{source}'][index[row['key']],k]
            assert abs(expected_score-value)<1e-14
            reference=references[(fold,scheme,source)][:,k]
            expected_p=(1+400-np.searchsorted(reference,value,side='left'))/401
            assert abs(expected_p-pv)<1e-14 and called==(pv<=alpha)
            key=(scheme,source,method,alpha,onset)
            unique=key+(row['key'],)
            assert unique not in seen
            seen.add(unique);counts[key][0]+=int(called);counts[key][1]+=1;rows+=1
    metrics_rows=0
    with (analysis/'heldout_metrics.csv').open() as stream:
        for row in csv.DictReader(stream):
            base=(row['scheme'],row['source'],row['method'],float(row['alpha']))
            tp,n=counts[base+(int(row['onset_years']),)]
            fp,nn=counts[base+(0,)]
            assert (n,nn)==(100,1000)
            assert (tp,fp)==(int(row['selected_called']),int(row['neutral_called']))
            power,fpr=tp/n,fp/nn
            assert abs(power-float(row['power']))<1e-14
            assert abs(fpr-float(row['neutral_call_fraction']))<1e-14
            for prevalence in (.5,.1,.01):
                den=prevalence*power+(1-prevalence)*fpr
                if den:
                    fdr=(1-prevalence)*fpr/den
                    assert abs(fdr-float(row[f'fdr_prevalence_{prevalence}']))<1e-14
            metrics_rows+=1
    target_counts=defaultdict(lambda:[0,0])
    target_reference={}
    with (analysis/'target70_predictions.csv').open() as stream:
        for row in csv.DictReader(stream):
            fold=int(row['fold']);target=int(row['target_onset']);onset=int(row['onset_years'])
            scheme,source,method=row['scheme'],row['source'],row['method']
            resolved=choices[(fold,scheme,source)] if method=='tuned_joint' else method
            k=mi[resolved];matrix=scores[f'fold{fold}_{scheme}_{source}']
            reference_key=(fold,scheme,source,method,target)
            if reference_key not in target_reference:
                train=[index[key] for key,v in inventory.items() if int(v['onset_years'])==target and int(v['fold'])!=fold]
                assert len(train)==80
                target_reference[reference_key]=(float(np.quantile(matrix[train,k],.3,method='lower')),{items[i]['key'] for i in train})
            threshold,training_ids=target_reference[reference_key]
            assert row['key'] not in training_ids
            value=matrix[index[row['key']],k]
            assert abs(threshold-float(row['threshold']))<1e-14
            assert abs(value-float(row['score']))<1e-14
            called=row['called']=='True'
            assert called==(value>=threshold and value>0)
            key=(scheme,source,method,target,onset)
            target_counts[key][0]+=int(called);target_counts[key][1]+=1
    with (analysis/'target70_metrics.csv').open() as stream:
        for row in csv.DictReader(stream):
            base=(row['scheme'],row['source'],row['method'],int(row['target_onset']))
            tp,n=target_counts[base+(int(row['onset_years']),)];fp,nn=target_counts[base+(0,)]
            assert (n,nn)==(100,1000)
            assert (tp,fp)==(int(row['selected_called']),int(row['neutral_called']))
            assert abs(tp/n-float(row['power']))<1e-14 and abs(fp/nn-float(row['neutral_call_fraction']))<1e-14
    ablation=out/'af_ablation/metrics.csv'
    if ablation.exists():
        with ablation.open() as stream:
            for row in csv.DictReader(stream):
                for source in ('truth','decoded'):
                    base=(row['scheme'],source,row['method'],float(row['alpha']))
                    assert counts[base+(int(row['onset_years']),)]==[int(row['selected_called']),int(row['selected_n'])]
                    assert counts[base+(0,)]==[int(row['neutral_called']),int(row['neutral_n'])]
    # Every derived class count stays within its denominator and partitions all pairs.
    positions=0
    for item in items:
        dest=out/'profiles'/item['key']
        receipt=json.loads((dest/'complete.json').read_text())
        spec=receipt['outputs']['features.npz']
        assert sha256(dest/'features.npz')==spec['sha256']
        with np.load(dest/'features.npz',allow_pickle=False) as f:
            for prefix in ('grid','site'):
                c,n=f[prefix+'_counts'],f[prefix+'_n']
                assert np.array_equal(c[:,:,:3].sum(axis=2),c[:,:,3])
                assert np.all(c<=n[None,:,:,None])
                assert np.all(np.diff(c.astype(int),axis=-1)>=0)
                if prefix=='site':
                    # A single archaic mutation older than every cutoff forces
                    # mixed-allele pairs to coalesce before that mutation.
                    assert not np.any(c[0,:,1]),item['key']
                positions+=c.shape[1]
    result=dict(status='passed',regions=len(items),prediction_rows=rows,metric_rows=metrics_rows,
                profiled_positions=positions,calibration_ranks_recomputed=True,
                all_regions_retained=True,disjoint_fit_calibration_test=True,
                class_partitions_verified=True,artifact_hashes_verified=True,
                target70_training_thresholds_recomputed=True,af_ablation_matches=ablation.exists())
    (analysis/'audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    audit(parser.parse_args().out)
