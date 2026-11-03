#!/usr/bin/env python3
"""Finish only previously submitted recent-selection units within a hard cap."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse
from concurrent.futures import ProcessPoolExecutor,wait,FIRST_COMPLETED
import importlib.util
import json
from pathlib import Path
import signal
import time
import run_recent_selection_pilot as current


def finish(root,seconds,workers,absolute_deadline=None):
    root=root.resolve();manifest=json.loads((root/'manifest.json').read_text());cfg=manifest['config']
    source=root/'executed_sources/run_recent_selection_pilot.py'
    assert current.base.digest(source)==cfg['source_sha256'],'Executed source differs from manifest'
    assert current.base.digest(current.ancestor.__file__)==cfg['ancestry_source_sha256']
    assert current.base.digest(current.base.__file__)==cfg['score_source_sha256']
    for name in ('slim','decoder'):assert current.base.digest(cfg[name])==cfg[name+'_sha256']
    spec=importlib.util.spec_from_file_location('executed_recent_runner',source)
    module=importlib.util.module_from_spec(spec)
    import sys
    sys.modules[spec.name]=module;spec.loader.exec_module(module)
    tasks=json.loads((root/'tasks.json').read_text());schedule=json.loads((root/'schedule.json').read_text())
    assert json.loads((root/'finished.json').read_text())['finished_epoch']>=schedule['started_epoch'],'Main runner has not finished'
    pending_tasks=[t for t in tasks if not any((root/t['unit']/name).exists() for name in ('complete.json','rejected.json'))]
    if not pending_tasks and (root/'finish_submitted.json').exists():
        print((root/'finish_submitted.json').read_text());return
    started=time.time();deadline=min(started+seconds,absolute_deadline) if absolute_deadline is not None else started+seconds;results=[]
    previous=json.loads((root/'status.json').read_text());byunit={r['unit']:r for r in previous}
    if pending_tasks:
        pool=ProcessPoolExecutor(max_workers=min(workers,len(pending_tasks)),initializer=os.setsid)
        pending={pool.submit(module.work,(str(root),cfg,t,deadline,schedule['parameters'].get('reuse_root'))):t for t in pending_tasks}
        try:
            while pending and time.time()<deadline:
                ready,_=wait(pending,timeout=min(1,max(0,deadline-time.time())),return_when=FIRST_COMPLETED)
                for future in ready:
                    task=pending.pop(future)
                    try:result=future.result()
                    except Exception as exc:result=dict(unit=task['unit'],state='failed',error=repr(exc))
                    results.append(result);byunit[task['unit']]=result
                current.base.save(root/'status.json',list(byunit.values()))
                current.base.save(root/'finish_status.json',results)
            if pending:
                for process in pool._processes.values():
                    try:os.killpg(process.pid,signal.SIGKILL)
                    except ProcessLookupError:pass
        finally:pool.shutdown(wait=True,cancel_futures=True)
    result=dict(started_epoch=started,finished_epoch=time.time(),deadline_epoch=deadline,seconds=seconds,workers=workers,
        units=[t['unit'] for t in pending_tasks],returned=len(results),source_sha256=current.base.digest(__file__),
        executed_runner_sha256=current.base.digest(source),new_proposals=0)
    current.base.save(root/'finish_submitted.json',result)
    print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',required=True,type=Path)
    p.add_argument('--seconds',type=int,default=150);p.add_argument('--workers',type=int,default=4)
    p.add_argument('--deadline-epoch',type=float,help='Optional absolute UTC Unix deadline, which overrides a longer relative cap')
    a=p.parse_args();finish(a.root,a.seconds,a.workers,a.deadline_epoch)
