"""Oracle focal-allele diagnostic, separate from blind whole-region detection."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import shutil

import numpy as np
import pandas as pd

from gamma_smc_aou.allele_class_ablation import local_scores,methods
from gamma_smc_aou.fresh_power import atomic_json,digest
from gamma_smc_aou.joint_scan_evaluation import load_one,metrics,rank_p


def choose_marker(sites,focal,selected):
    if not len(sites):
        if selected:raise ValueError('Selected focal allele is missing')
        return None
    k=int(np.argmin(np.abs(sites-focal)))
    if selected and sites[k]!=focal:raise ValueError('Selected focal SNP not found exactly')
    return k


def evaluate(input_dir,out,workers,export):
    out.mkdir(parents=True,exist_ok=True)
    manifest=json.loads((input_dir/'profile_manifest.json').read_text())
    items=[i for i in manifest['items'] if i['onset_years'] in (0,50000)]
    if len(items)!=1100:raise ValueError('Expected 1100 regions')
    cfg=manifest['config'];method_list=[m for m in methods(cfg['tmrca_cutoffs_years']) if m['gate']==0]
    inventory=pd.read_csv(input_dir/'analysis/region_inventory.csv').set_index('key')
    folds=np.array([inventory.loc[i['key'],'fold'] for i in items],dtype=int)
    with ThreadPoolExecutor(max_workers=workers) as pool:loaded=list(pool.map(load_one,[(i,input_dir) for i in items]))
    scores=np.zeros((2,1100,len(method_list)));records=[]
    for i,(item,(f,receipt,_)) in enumerate(zip(items,loaded)):
        focal=int(item['cfg']['focal_position_bp'])
        k=choose_marker(f['sites'],focal,item['onset_years']==50000)
        record=dict(key=item['key'],onset_years=item['onset_years'],fold=int(folds[i]),
                    marker_position=None if k is None else int(f['sites'][k]),
                    distance_to_anchor=None if k is None else int(abs(f['sites'][k]-focal)),
                    af=0 if k is None else float(f['site_af'][k]),
                    nAA=0 if k is None else int(f['site_n'][k,2]),nRR=0 if k is None else int(f['site_n'][k,0]))
        if k is not None:
            values=local_scores(f['site_counts'][:,k:k+1],f['site_n'][k:k+1],np.array([True]))
            for j,m in enumerate(method_list):
                scores[:,i,j]=f['site_af'][k] if m['family']=='af' else values[m['family']][:,0,m['cutoff']]
        records.append(record)
    regions=pd.DataFrame(records);neutral=regions.onset_years.to_numpy()==0;predictions=[]
    assert neutral.sum()==1000 and (~neutral).sum()==100
    for fold in range(5):
        test=folds==fold;cal=neutral&np.isin(folds,[(fold+3)%5,(fold+4)%5])
        assert cal.sum()==400 and (test&neutral).sum()==200 and (test&~neutral).sum()==20
        for src,source in enumerate(('truth','decoded')):
            for j,m in enumerate(method_list):
                ids=np.flatnonzero(test);p=rank_p(scores[src,cal,j],scores[src,test,j])
                for alpha in (.01,.05,.1):
                    for idx,pv in zip(ids,p):
                        predictions.append(dict(key=items[idx]['key'],onset_years=items[idx]['onset_years'],fold=fold,
                            scheme='oracle_focal_allele',source=source,method=m['name'],alpha=alpha,
                            score=float(scores[src,idx,j]),p=float(pv),called=bool(pv<=alpha)))
    pred=pd.DataFrame(predictions);table=metrics(pred)
    pred.to_csv(out/'predictions.csv.gz',index=False,compression=dict(method='gzip',mtime=0))
    table.to_csv(out/'metrics.csv',index=False);regions.to_csv(out/'regions.csv',index=False)
    np.savez_compressed(out/'scores.npz',scores=scores)
    atomic_json(out/'methods.json',method_list)
    atomic_json(out/'provenance.json',dict(source_sha256=digest(Path(__file__)),profile_manifest_sha256=digest(input_dir/'profile_manifest.json'),
        class_score_code_sha256=digest(Path(__file__).resolve().parents[1]/'python/gamma_smc_aou/allele_class_ablation.py'),
        input_dir=str(input_dir),workers=workers,config=cfg,
        selected_coordinate='Known selected focal SNP; exact eligible-site match required',
        neutral_coordinate='Nearest eligible archaic SNP to the prespecified focal coordinate; ties left',
        endpoint='Oracle single-candidate diagnostic, not blind region scan power',
        introduced_after_regional_ablation_results=True,later_onset_excluded=True,
        profile_hashes={item['key']:loaded[i][2] for i,item in enumerate(items)}))
    # Re-read persisted values and audit by direct counting, not rank_p.
    checked=pd.read_csv(out/'predictions.csv.gz');index={item['key']:i for i,item in enumerate(items)}
    method_index={m['name']:j for j,m in enumerate(method_list)}
    with np.load(out/'scores.npz',allow_pickle=False) as saved:stored=saved['scores']
    for (fold,source,method),group in checked.groupby(['fold','source','method']):
        src=0 if source=='truth' else 1;j=method_index[method]
        ids=np.array([index[k] for k in group.key]);assert np.all(folds[ids]==fold)
        cal=neutral&np.isin(folds,[(fold+3)%5,(fold+4)%5])
        values=stored[src,ids,j];reference=stored[src,cal,j]
        pv=(1+(reference[:,None]>=values).sum(axis=0))/401
        np.testing.assert_allclose(values,group.score,rtol=0,atol=1e-14)
        np.testing.assert_allclose(pv,group.p,rtol=0,atol=1e-14)
        np.testing.assert_array_equal(pv<=group.alpha,group.called)
    counts=checked.groupby(['source','method','alpha','onset_years']).called.agg(['size','sum'])
    for row in pd.read_csv(out/'metrics.csv').to_dict('records'):
        base=(row['source'],row['method'],row['alpha'])
        sn,nn=counts.loc[base+(50000,)],counts.loc[base+(0,)]
        assert sn['size']==100 and nn['size']==1000
        assert sn['sum']==row['selected_called'] and nn['sum']==row['neutral_called']
        assert abs(sn['sum']/100-row['power'])<1e-14 and abs(nn['sum']/1000-row['neutral_call_fraction'])<1e-14
    atomic_json(out/'audit.json',dict(status='passed',prediction_rows=len(checked),selected_focal_sites_matched=100,
        neutral_marker_free=int(((regions.onset_years==0)&regions.marker_position.isna()).sum()),calibration_ranks_directly_counted=True))
    names=['metrics.csv','regions.csv','methods.json','provenance.json','audit.json','predictions.csv.gz','scores.npz']
    atomic_json(out/'artifact_manifest.json',{name:dict(bytes=(out/name).stat().st_size,sha256=digest(out/name)) for name in names})
    if export:
        export.mkdir(parents=True,exist_ok=True)
        for name in names+['artifact_manifest.json']:shutil.copyfile(out/name,export/name)
    print(table[(table.alpha==.05)&table.method.str.endswith('T50000')][['source','method','power','neutral_call_fraction']].to_string(index=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--workers',type=int,default=20);p.add_argument('--export',type=Path);a=p.parse_args()
    evaluate(a.input,a.out,a.workers,a.export)
