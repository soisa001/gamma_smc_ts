"""Compute the prespecified AF-only ablation directly from recovered marker labels.

This does not use TMRCA, so it can run after labels finish while decoding runs.
The final joint evaluation must reproduce these values from its feature files.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import numpy as np
import pandas as pd
from gamma_smc_aou.joint_scan_evaluation import make_folds,binned_region_score,rank_p,metrics
from gamma_smc_aou.joint_scan_profiles import nearest_markers
from gamma_smc_aou.fresh_power import digest,atomic_json


def one(payload):
    item,out,cfg=payload
    label=out/'labels'/item['key']
    receipt=json.loads((label/'complete.json').read_text())
    assert digest(label/'markers.npz')==receipt['outputs']['markers.npz']['sha256']
    with np.load(label/'markers.npz',allow_pickle=False) as f:
        positions,af=f['positions'],f['af']
    grid=np.arange(0,cfg['scored_length_bp'],cfg['stride_bp'])
    index=nearest_markers(grid,positions,cfg['stride_bp']//2)
    grid_af=np.zeros(len(grid));grid_af[index>=0]=af[index[index>=0]]
    return np.array([[binned_region_score(values,pos,len(grid),cfg['stride_bp'],run) for run in (1,2,3)]
                     for values,pos in ((grid_af,grid),(af,positions))])


def main(out):
    status=json.loads((out/'label_status.json').read_text())
    if status['state']!='complete' or status['completed']!=1200:raise RuntimeError('Wait for all label checks')
    manifest=json.loads((out/'label_manifest.json').read_text());cfg=manifest['config']
    items=[]
    for key in manifest['samples']:
        arm,rep=key.split('/')
        items.append(dict(key=key,onset_years=0 if arm=='neutral' else int(arm.removeprefix('onset')),replicate=int(rep.removeprefix('rep'))))
    folds=make_folds(items,cfg['seed_base']);onsets=np.array([i['onset_years'] for i in items])
    with ThreadPoolExecutor(max_workers=4) as pool:scores=np.array(list(pool.map(one,[(i,out,cfg) for i in items])))
    rows=[]
    for fold in range(5):
        cal=(onsets==0)&np.isin(folds,[(fold+3)%5,(fold+4)%5]);test=folds==fold
        for scheme_index,scheme in enumerate(('grid10kb','archaic_sites')):
            for k,run in enumerate((1,2,3)):
                p=rank_p(scores[cal,scheme_index,k],scores[test,scheme_index,k])
                for alpha in (.01,.05,.1):
                    for i,pv in zip(np.where(test)[0],p):
                        rows.append(dict(key=items[i]['key'],onset_years=int(onsets[i]),fold=fold,scheme=scheme,
                                         source='label_AF',method=f'af_r{run}',alpha=alpha,p=float(pv),called=bool(pv<=alpha)))
    dest=out/'af_ablation';dest.mkdir(exist_ok=True)
    predictions=pd.DataFrame(rows);table=metrics(predictions)
    predictions.to_csv(dest/'predictions.csv.gz',index=False,compression=dict(method='gzip',mtime=0))
    table.to_csv(dest/'metrics.csv',index=False)
    atomic_json(dest/'provenance.json',dict(label_manifest_sha256=digest(out/'label_manifest.json'),seed=cfg['seed_base'],source_sha256=digest(Path(__file__))))
    print(table[(table.alpha==.05)&(table.method=='af_r1')][['scheme','onset_years','selected_called','selected_n','neutral_called','neutral_n','power','neutral_call_fraction']].to_string(index=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True)
    main(p.parse_args().out)
