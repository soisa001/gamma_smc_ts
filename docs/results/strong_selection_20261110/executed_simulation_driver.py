#!/usr/bin/env python3
"""Bounded six-treatment strong-selection extension with segregation conditioning.

Keep the original EAS/introgression physics. Reuse the frozen 91 accepted recent
s=.05 simulations; generate independent additional proposals in fixed order.
Completed and rejected receipts are resumable and all parameters are recorded.
"""
import os
os.environ['MPLBACKEND']='Agg'
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import argparse
from concurrent.futures import ProcessPoolExecutor,wait,FIRST_COMPLETED
import json
from pathlib import Path
import shutil
import time
import numpy as np
import pandas as pd
import run_selection_arrays as arrays
from report_selection_arrays import local_path

base=arrays.base
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'docs/results/selection_arrays_eas_20261105_4h'
REFERENCE=ROOT/'docs/results/selection_calls_100_segregating_null_20261110'
STRATA=[dict(key=f'{regime}_s{s:g}',regime=regime,selection=s,start=start,end=end)
        for regime,start,end in [('10k',10000,0),('50k',50000,0),('50_30k',50000,30000)] for s in (.05,.1)]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--raw',type=Path,default=ROOT.parent/'sim/strong_selection_20261110')
    p.add_argument('--target',type=int,default=100)
    p.add_argument('--workers',type=int,default=4)
    p.add_argument('--seconds',type=int,default=14400)
    p.add_argument('--deadline-epoch',type=float)
    p.add_argument('--attempts',type=int,default=5000,help='Maximum new proposals per treatment, also bounded by wall time.')
    args=p.parse_args();assert 1<=args.workers<=24 and args.target>0 and args.seconds>0
    raw=args.raw;raw.mkdir(parents=True,exist_ok=True)
    reference=json.loads((REFERENCE/'manifest.json').read_text())
    template=reference['simulation_configurations']['50k_s001']
    configs={s['regime']:dict(template,selection_start_years=s['start'],selection_end_years=s['end'],
        seed_namespace=f"strong-selection-v1:{s['regime']}:") for s in STRATA}
    for key,path in [('ancestry_source_sha256',arrays.ancestor.__file__),('score_source_sha256',base.__file__),
                     ('slim_sha256',template['slim']),('decoder_sha256',template['decoder']),
                     ('phlash_sha256',ROOT/'no_introgression_EAS_sim/resources/EAS.npz')]:
        assert base.digest(path)==template[key],str(path)
    lib=Path(template['slim']).parent.parent/'lib'
    os.environ['LD_LIBRARY_PATH']=str(lib)+(':'+os.environ['LD_LIBRARY_PATH'] if os.environ.get('LD_LIBRARY_PATH') else '')
    inventory=json.loads((BASE/'output_hashes.json').read_text())
    assert base.digest(BASE/'completed_units.tsv')==inventory['completed_units.tsv']
    table=pd.read_csv(BASE/'completed_units.tsv',sep='\t')
    reused={s['key']:[] for s in STRATA}
    for s in STRATA:
        rows=table[(table.origin=='introgression')&(table.cohort=='selected')&(table.regime==s['regime'])&(table.selection==s['selection'])].sort_values('replicate').iloc[:args.target]
        if s['regime']=='50_30k':assert rows.empty
        for r in json.loads(rows.to_json(orient='records',double_precision=15)):
            directory=local_path(r['source_directory'])
            arrays.ancestor.verify(directory,r)
            assert base.digest(directory/'trajectory.tsv')==r['trajectory_sha256']
            assert r['accepted'] and 0<r['final_frequency']<1
            reused[s['key']].append(r)
    contract=dict(configs=configs,strata=STRATA,target=args.target,reused=reused,
        conditioning='Original focal allele remains segregating today: 0 < whole-population AF < 1; fixed/lost proposals excluded.',
        selection_rule='Frozen surviving historical cohort first, followed by new surviving proposals in ascending index order; no AF or inferred-signal ascertainment beyond segregation.',
        source_table_sha256=inventory['completed_units.tsv'],reference_manifest_sha256=base.digest(REFERENCE/'manifest.json'),
        source_hashes={Path(path).name:base.digest(path) for path in (arrays.__file__,arrays.ancestor.__file__,base.__file__)})
    manifest=raw/'manifest.json'
    if manifest.exists():assert json.loads(manifest.read_text())==contract
    else:base.save(manifest,contract)
    with np.load(REFERENCE/'neutral.npz') as data:pairs=data['pairs'].copy()
    for name,values in [('pairs.tsv',pairs),('positions.txt',np.arange(100)*10000)]:
        path=raw/name
        if path.exists():np.testing.assert_array_equal(np.loadtxt(path,dtype=int),values)
        else:np.savetxt(path,values,fmt='%d',delimiter='\t')
    snapshot=raw/'executed_sources';snapshot.mkdir(exist_ok=True)
    for path in map(Path,(__file__,arrays.__file__,arrays.ancestor.__file__,base.__file__)):
        target=snapshot/f'{path.stem}_{base.digest(path)[:12]}.py'
        if not target.exists():shutil.copy2(path,target)
    start=time.time();deadline=min(start+args.seconds,args.deadline_epoch or float('inf'))
    if (raw/'schedule.json').exists():
        previous=json.loads((raw/'schedule.json').read_text())
        shutil.copy2(raw/'schedule.json',raw/f"schedule_{int(previous['started_epoch'])}.json")
    base.save(raw/'schedule.json',dict(started_epoch=start,deadline_epoch=deadline,seconds=args.seconds,workers=args.workers,
        attempts_per_treatment=args.attempts,target=args.target,script_sha256=base.digest(__file__),native_library_path=str(lib)))
    completed={s['key']:[] for s in STRATA};next_index={s['key']:0 for s in STRATA};pending={};cursor=0
    last_ledger=0.;last_log=0.
    def accepted(s):return len(reused[s['key']])+sum(r['accepted'] for r in completed[s['key']])
    def progress():
        return dict(elapsed_seconds=time.time()-start,deadline_epoch=deadline,pending=len(pending),strata=[dict(s,
            reused=len(reused[s['key']]),new_completed=len(completed[s['key']]),accepted=accepted(s),target=args.target,
            lost=sum(r['rejection']=='loss' for r in completed[s['key']]),fixed=sum(r['rejection']=='fixation' for r in completed[s['key']]),
            submitted=next_index[s['key']]) for s in STRATA])
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        while True:
            while len(pending)<args.workers and time.time()<deadline-120:
                chosen=None
                for _ in STRATA:
                    s=STRATA[cursor%len(STRATA)];cursor+=1
                    if accepted(s)<args.target and next_index[s['key']]<args.attempts:chosen=s;break
                if chosen is None:break
                s=chosen;i=next_index[s['key']];next_index[s['key']]+=1
                unit=f"introgression/s{s['selection']:g}/selected_{i:05d}"
                task=dict(id=f"{s['regime']}/{unit}",unit=unit,regime=s['regime'],origin='introgression',cohort='selected',selection=s['selection'],replicate=i)
                f=pool.submit(arrays.work,(str(raw),configs[s['regime']],task,deadline));pending[f]=(s,task)
            if not pending:break
            done,_=wait(pending,return_when=FIRST_COMPLETED)
            for f in done:
                s,task=pending.pop(f);result=f.result()
                assert result['state'] in ('accepted','rejected'),result
                path=raw/task['id']/('complete.json' if result['state']=='accepted' else 'rejected.json')
                record=json.loads(path.read_text());record['source_directory']=str(path.parent)
                completed[s['key']].append(record);completed[s['key']].sort(key=lambda r:r['replicate'])
            # Each simulation already writes its own complete/rejected receipt.
            # Avoid rewriting an ever-growing combined ledger for every result.
            if time.time()-last_ledger>=30:
                arrays.save_ledger(raw/'proposal_ledger.json',completed);last_ledger=time.time()
            value=progress();arrays.save_ledger(raw/'progress.json',value)
            if time.time()-last_log>=10:
                print(json.dumps(value),flush=True);last_log=time.time()
    arrays.save_ledger(raw/'proposal_ledger.json',completed)
    for s in STRATA:
        rows=completed[s['key']]
        assert [r['replicate'] for r in rows]==list(range(next_index[s['key']]))
        kept=(reused[s['key']]+[r for r in rows if r['accepted']])[:args.target]
        base.save(raw/(s['key']+'_cohort.json'),dict(stratum=s,target=args.target,complete=len(kept)==args.target,
            accepted=kept,reused_count=len(reused[s['key']]),new_proposals=len(rows),
            new_lost=sum(r['rejection']=='loss' for r in rows),new_fixed=sum(r['rejection']=='fixation' for r in rows),
            manifest_sha256=base.digest(manifest)))
    base.save(raw/'finished.json',progress())
    print('Finished bounded strong-selection arrays.',flush=True)


if __name__=='__main__':main()
