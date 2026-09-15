"""Audit saved allele-class ranks/counts and optionally export compact artifacts."""
import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import shutil

import numpy as np
import pandas as pd


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(out,export=None):
    manifest=json.loads((out/'artifact_manifest.json').read_text())
    for name,spec in manifest.items():
        assert (out/name).stat().st_size==spec['bytes'] and digest(out/name)==spec['sha256'],name
    inventory=pd.read_csv(out/'regions.csv');keys=inventory.key.tolist()
    index={key:k for k,key in enumerate(keys)};folds=inventory.fold.to_numpy();onset=inventory.onset_years.to_numpy()
    assert len(index)==1100 and sum(onset==0)==1000 and sum(onset==50000)==100 and set(onset)=={0,50000}
    method_list=json.loads((out/'methods.json').read_text());mi={m['name']:k for k,m in enumerate(method_list)}
    schemes={'grid10kb':0,'archaic_sites':1};sources={'truth':0,'decoded':1}
    with np.load(out/'scores.npz',allow_pickle=False) as saved:scores=saved['scores']
    assert scores.shape==(5,2,2,1100,len(mi)) and np.isfinite(scores).all() and np.all((scores>=0)&(scores<=1))
    fold_records=json.loads((out/'folds.json').read_text())
    refs={};thresholds={}
    for f in range(5):
        test=set(inventory[inventory.fold==f].key)
        fit=set(fold_records[f]['fit_ids']);cal=set(fold_records[f]['calibration_ids'])
        assert len(test)==220 and len(fit)==len(cal)==400
        assert not fit&cal and not fit&test and not cal&test
        assert all(onset[index[k]]==0 for k in fit|cal)
        cal_indices=[index[k] for k in cal]
        train=(onset==50000)&(folds!=f);assert train.sum()==80
        for si in range(2):
            for src in range(2):
                refs[f,si,src]=np.sort(scores[f,si,src,cal_indices],axis=0)
                thresholds[f,si,src]=np.quantile(scores[f,si,src,train],.3,axis=0,method='lower')
    prediction_rows=0
    for filename,metric_name,target in [('predictions.csv.gz','metrics.csv',False),('target70_predictions.csv.gz','target70_metrics.csv',True)]:
        pred=pd.read_csv(out/filename)
        columns=['key','scheme','source','method','alpha']
        assert not pred.duplicated(columns).any()
        for (f,scheme,source,method),group in pred.groupby(['fold','scheme','source','method']):
            si,src,k=schemes[scheme],sources[source],mi[method]
            ids=np.array([index[key] for key in group.key])
            assert np.all(folds[ids]==f) and np.array_equal(onset[ids],group.onset_years.to_numpy())
            assert set(group.key)==set(inventory[inventory.fold==f].key)
            value=scores[f,si,src,ids,k]
            np.testing.assert_allclose(value,group.score,atol=1e-14,rtol=0)
            if target:
                threshold=thresholds[f,si,src][k]
                np.testing.assert_allclose(group.threshold,threshold,atol=1e-14,rtol=0)
                expected=(value>=threshold)&(value>0)
            else:
                ref=refs[f,si,src][:,k]
                expected_p=(401-np.searchsorted(ref,value,side='left'))/401
                np.testing.assert_allclose(expected_p,group.p,atol=1e-14,rtol=0)
                expected=expected_p<=group.alpha.to_numpy()
            assert np.array_equal(expected,group.called.to_numpy())
        for row in pd.read_csv(out/metric_name).to_dict('records'):
            subset=pred[(pred.scheme==row['scheme'])&(pred.source==row['source'])&(pred.method==row['method'])&(pred.alpha==row['alpha'])]
            selected=subset[subset.onset_years==50000];neutral=subset[subset.onset_years==0]
            assert len(selected)==100 and len(neutral)==1000
            tp,fp=int(selected.called.sum()),int(neutral.called.sum())
            assert tp==row['selected_called'] and fp==row['neutral_called']
            power,fpr=tp/100,fp/1000
            assert abs(power-row['power'])<1e-14 and abs(fpr-row['neutral_call_fraction'])<1e-14
            for pi in (.5,.1,.01):
                den=pi*power+(1-pi)*fpr
                if den:assert abs((1-pi)*fpr/den-row[f'fdr_prevalence_{pi}'])<1e-14
                expected_f1=2*pi*power/(pi*(1+power)+(1-pi)*fpr)
                assert abs(expected_f1-row[f'f1_prevalence_{pi}'])<1e-14
        prediction_rows+=len(pred)
    # Direct ablation property across every region/fold and all six cutoffs.
    equality=0
    for name,k in mi.items():
        if name.startswith('joint_'):
            other=mi[name.replace('joint_','weighted_contrast_',1)]
            np.testing.assert_array_equal(scores[:,1,0,:,k],scores[:,1,0,:,other]);equality+=5500
    result=dict(status='passed',regions=1100,selected=100,neutral=1000,prediction_rows=prediction_rows,
                no_later_onset=True,calibration_ranks_recomputed=True,target70_thresholds_recomputed=True,
                truth_exact_site_penalty_equality_checks=equality,artifact_hashes_verified=True,
                script_sha256=digest(Path(__file__)))
    (out/'audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    if export:
        export.mkdir(parents=True,exist_ok=True)
        omitted={'predictions.csv.gz','target70_predictions.csv.gz','scores.npz','class_maxima.csv.gz'}
        names=[name for name in manifest if name not in omitted]+['audit.json']
        for name in names:shutil.copyfile(out/name,export/name)
        shutil.copyfile(out/'artifact_manifest.json',export/'full_analysis_manifest.json');names.append('full_analysis_manifest.json')
        target_path=export/'primary_predictions.csv.gz'
        with gzip.open(out/'predictions.csv.gz','rt',newline='') as source,target_path.open('wb') as raw:
            with gzip.GzipFile(fileobj=raw,mode='wb',mtime=0) as zipped,io.TextIOWrapper(zipped,newline='') as dest:
                reader=csv.DictReader(source);writer=csv.DictWriter(dest,fieldnames=reader.fieldnames);writer.writeheader()
                for row in reader:
                    if row['method'].endswith('T50000') or row['method']=='af_r1':writer.writerow(row)
        names.append(target_path.name)
        (export/'artifact_manifest.json').write_text(json.dumps(dict(source=str(out),files={name:dict(bytes=(export/name).stat().st_size,sha256=digest(export/name)) for name in names}),indent=2)+'\n')
        print(f'Exported {len(names)+1} files to {export}',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--export',type=Path);args=p.parse_args();audit(args.out,args.export)
