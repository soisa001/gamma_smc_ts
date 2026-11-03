#!/usr/bin/env python3
"""Exact rejection on present-day segregation, with absorbing-state early exits."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
from pathlib import Path
import shutil
import signal
import sys
import time
import run_ongoing_selection_pilot as model
base=model.base


def accepted(frequency):
    return 0 < frequency < 1


def work(payload):
    root,source,cfg,task,deadline=payload
    if time.time()>deadline:return dict(unit=task['unit'],state='not_started_deadline')
    root=Path(root);source=Path(source);directory=root/task['unit'];directory.mkdir(parents=True,exist_ok=True)
    record=dict(task,seed=base.seed('ongoing-v1:'+task['unit']),config_sha256=base.digest(root/'manifest.json'))
    for name in ('complete.json','rejected.json'):
        if (directory/name).exists():
            previous=json.loads((directory/name).read_text());assert previous['config_sha256']==record['config_sha256']
            if name=='complete.json':model.verify(directory,previous)
            return dict(unit=task['unit'],state='accepted' if name=='complete.json' else 'rejected')
    old=source/task['unit']/'simulation.json'
    if not old.exists() and (source/task['unit']/'rejected.json').exists():old=source/task['unit']/'rejected.json'
    if old.exists():
        prior=json.loads(old.read_text());record.update({k:v for k,v in prior.items() if k!='config_sha256'})
        record['source_receipt']=str(old);record['source_receipt_sha256']=base.digest(old)
        record['accepted']=accepted(record['final_frequency'])
        if record['accepted']:
            for name in ('onset.trees','observed.trees','trajectory.tsv','onset_provenance.json'):
                shutil.copy2(old.parent/name,directory/name)
            assert base.digest(directory/'observed.trees')==record['observed_sha256']
        else:record['rejection']='loss' if record['final_frequency']==0 else 'fixation'
    else:
        record.update(model.simulate(cfg,task,directory));record['accepted']=accepted(record['final_frequency'])
    if not record['accepted']:
        base.save(directory/'rejected.json',record)
        return dict(unit=task['unit'],state='rejected')
    base.save(directory/'simulation.json',record)
    record.update(base.score(cfg,directory,root));base.save(directory/'complete.json',record)
    return dict(unit=task['unit'],state='accepted')


def launch(a):
    root=Path(a.out).resolve();source=Path(a.source).resolve();root.mkdir(parents=True,exist_ok=True)
    old=json.loads((source/'manifest.json').read_text())
    cfg=dict(old['config'],schema='ongoing-eas-segregating/v1',condition_segregating=True,
        conditioning='0 < whole-population focal allele frequency at present < 1 for selected simulations',
        condition_neutral=a.condition_neutral,
        source_sha256=base.digest(__file__),model_source_sha256=base.digest(model.__file__))
    if (root/'manifest.json').exists():assert json.loads((root/'manifest.json').read_text())['config']==cfg
    else:base.save(root/'manifest.json',dict(config=cfg,versions=old['versions'],source_root=str(source),created_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())))
    for name in ('positions.txt','pairs.tsv'):shutil.copy2(source/name,root/name)
    tasks=[]
    # Preassigned cohort membership is retained even after rejection. No accepted
    # replicate is moved between development, calibration, and test sets.
    for i in range(max(a.attempts,a.calibration_attempts)):
        for origin in ('de_novo','introgression'):
            for cohort in (('development','calibration','neutral_test') if a.condition_neutral else ()):
                if i<(a.calibration_attempts if cohort=='calibration' else a.attempts):
                    tasks.append(dict(unit=f'{origin}/{cohort}_{i:04d}',origin=origin,cohort=cohort,selection=0.,replicate=i))
            if i<a.attempts:
                for s in a.selection:
                    tasks.append(dict(unit=f'{origin}/s{s:g}/selected_{i:04d}',origin=origin,cohort='selected',selection=s,replicate=i))
    base.save(root/'tasks.json',tasks)
    deadline=time.time()+a.seconds
    base.save(root/'schedule.json',dict(workers=a.workers,deadline_epoch=deadline,attempts_per_stratum=a.attempts,
        calibration_attempts=a.calibration_attempts,selection_coefficients=a.selection,source_root=str(source)))
    records=[]
    # Each child establishes its own process group, so the deadline terminates
    # its SLiM/decoder descendants too. Persistent workers avoid import overhead.
    pool=ProcessPoolExecutor(max_workers=a.workers,initializer=os.setsid)
    futures=[pool.submit(work,(str(root),str(source),cfg,t,deadline)) for t in tasks]
    try:
        for f in as_completed(futures,timeout=a.seconds):
            try:result=f.result()
            except Exception as exc:result=dict(unit=tasks[futures.index(f)]['unit'],state='failed',error=repr(exc))
            records.append(result);base.save(root/'status.json',records)
    except TimeoutError:
        for process in pool._processes.values():
            try:os.killpg(process.pid,signal.SIGKILL)
            except ProcessLookupError:pass
    finally:pool.shutdown(wait=True,cancel_futures=True)
    base.save(root/'status.json',records)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',required=True);p.add_argument('--source',required=True)
    p.add_argument('--workers',type=int,default=4);p.add_argument('--seconds',type=int,default=600);p.add_argument('--attempts',type=int,default=100)
    p.add_argument('--calibration-attempts',type=int,default=200)
    p.add_argument('--selection',type=float,nargs='+',default=[.005,.01,.02,.03])
    p.add_argument('--condition-neutral',action='store_true',help='Experimental alternative; default uses a separate unconditioned neutral reference')
    launch(p.parse_args())
