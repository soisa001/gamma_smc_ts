#!/usr/bin/env python3
"""Exact de novo rejection sampling with one independent founder per slot.

Pre-mutation ancestry is independent of focal-allele survival: there are no
other selected mutations and every initial genotype configuration has one
heterozygote. Retrying forward histories against a slot's fixed ancestry thus
leaves its founder distribution unchanged. Different retained slots never share
founders. Mutation, sampling, forward and ancestry seeds are separately recorded.
"""
import os
os.environ['MPLBACKEND'] = 'Agg'
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'
import argparse
import ast
import hashlib
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
import json
from pathlib import Path
import shutil
import signal
import sys
import time
import numpy as np
import run_selection_arrays as arrays

ROOT = Path(__file__).resolve().parents[1]
base = arrays.base
ORIGINAL_ONSET = arrays.ancestor.onset


def install_founder(founder, cfg, origin, label, directory):
    assert origin == 'de_novo'
    meta = json.loads((founder / 'founder.json').read_text())
    assert cfg['ancestry_start_generations'] == meta['age_generations']
    for name in ('onset.trees', 'onset_provenance.json'):
        destination = directory / name
        if destination.exists():
            assert base.digest(destination) == base.digest(founder / name)
        else:
            os.symlink(os.path.relpath(founder / name, directory), destination)
            assert base.digest(destination) == base.digest(founder / name)
    return meta['carrier_ids'], meta['sizes']


def worker(raw, cfg, stratum, slot, attempts, deadline):
    raw = Path(raw); key = stratum['key']
    slotdir = raw / key / f'slot_{slot:03d}'; slotdir.mkdir(parents=True, exist_ok=True)
    signature = base.digest(raw / 'manifest.json')
    founder = slotdir / 'founder'; founder.mkdir(exist_ok=True)
    founder_label = f'denovo-independent-founder-v1:{key}:slot_{slot:03d}'
    founder_meta = founder / 'founder.json'
    if founder_meta.exists():
        meta = json.loads(founder_meta.read_text())
        assert meta['config_sha256'] == signature
        assert base.digest(founder / 'onset.trees') == meta['onset_sha256']
    else:
        age = stratum['start'] // 25
        carriers, sizes = ORIGINAL_ONSET(dict(cfg, ancestry_start_generations=age), 'de_novo', founder_label, founder)
        assert len(carriers) == 1
        meta = dict(carrier_ids=carriers, sizes=sizes, age_generations=age, seed_label=founder_label,
                    ancestry_seed=base.seed(founder_label+':ancestry'), pairing_seed=base.seed(founder_label+':pairing'),
                    config_sha256=signature, onset_sha256=base.digest(founder / 'onset.trees'))
        base.save(founder_meta, meta)
    records = []
    arrays.ancestor.onset = lambda c, o, l, d: install_founder(founder, c, o, l, d)
    try:
        for attempt in range(attempts):
            if time.time() >= deadline-90: break
            unit = f"de_novo/{key}/slot_{slot:03d}/attempt_{attempt:05d}"
            task = dict(id=unit, unit=unit, regime=stratum['regime'], origin='de_novo', cohort='selected',
                scientific_role=stratum['role'], selection=stratum['selection'], replicate=slot*attempts+attempt,
                slot=slot, attempt=attempt, founder_seed_label=founder_label,
                founder_ancestry_seed=meta['ancestry_seed'], founder_pairing_seed=meta['pairing_seed'])
            result = arrays.work((str(raw), cfg, task, deadline))
            assert result['state'] in ('accepted', 'rejected'), result
            path = raw / unit / ('complete.json' if result['state']=='accepted' else 'rejected.json')
            record = json.loads(path.read_text()); record['source_directory'] = str(path.parent)
            records.append(record)
            status = dict(key=key, slot=slot, proposals=len(records), lost=sum(r['rejection']=='loss' for r in records),
                fixed=sum(r['rejection']=='fixation' for r in records), accepted=record if record['accepted'] else None,
                complete=record['accepted'] or len(records)==attempts, founder_sha256=meta['onset_sha256'])
            arrays.save_ledger(slotdir / 'status.json', status)
            if record['accepted']: break
        return dict(key=key, slot=slot, records=records, complete=bool(records and (records[-1]['accepted'] or len(records)==attempts)))
    finally:
        arrays.ancestor.onset = ORIGINAL_ONSET


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--previous', type=Path, default=ROOT.parent/'sim/denovo_comparison_20261111')
    parser.add_argument('--raw', type=Path, default=ROOT.parent/'sim/denovo_comparison_cached_20261111')
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--seconds', type=int, default=14400)
    parser.add_argument('--deadline-epoch', type=float)
    parser.add_argument('--attempts-per-slot', type=int, default=10000)
    parser.add_argument('--older-strong-attempts-per-slot', type=int, default=50)
    args = parser.parse_args(); assert 1 <= args.workers <= 24
    previous, raw = args.previous.resolve(), args.raw.resolve(); raw.mkdir(parents=True, exist_ok=True)
    old = json.loads((previous/'manifest.json').read_text())
    strata, target = old['strata'], old['target']
    reused = {s['key']: json.loads((previous/(s['key']+'_cohort.json')).read_text())['accepted'] for s in strata}
    configs = {r: dict(c, seed_namespace=f'denovo-cached-forward-v1:{r}:') for r, c in old['configs'].items()}
    for key, digest in old['scientific_source_hashes'].items():
        assert base.digest(Path(arrays.__file__).parent/key) == digest
    contract = dict(old, configs=configs, reused=reused, previous_manifest_sha256=base.digest(previous/'manifest.json'),
        founder_design='Each missing replicate slot receives an independent msprime founder genealogy. Only forward history is retried within that slot. At most one retained simulation per founder.',
        independence_basis='One initial focal heterozygote; no other fitness mutations; fixed demographic schedule. Pre-mutation genealogy does not enter focal fitness or survival probability.',
        wrapper_sha256=base.digest(__file__),
        attempt_limits=dict(default_per_slot=args.attempts_per_slot, older_s005_per_slot=args.older_strong_attempts_per_slot),
        selection_rule='Previous verified survivors, then completed slots in fixed stratum/slot order. An incomplete slot stops the reportable prefix; exhausted slots count as completed without a retained simulation.')
    manifest = raw/'manifest.json'
    if manifest.exists():
        saved = json.loads(manifest.read_text())
        if saved != contract:
            assert {k:v for k,v in saved.items() if k!='wrapper_sha256'} == {k:v for k,v in contract.items() if k!='wrapper_sha256'}
            initial = raw/'executed_sources'/f"run_denovo_cached_founders_{saved['wrapper_sha256'][:12]}.py"
            assert base.digest(initial) == saved['wrapper_sha256']
            def physics(path):
                return {n.name:ast.dump(n,include_attributes=False) for n in ast.parse(path.read_text()).body
                        if isinstance(n,ast.FunctionDef) and n.name in ('worker','install_founder')}
            old_physics, new_physics = physics(initial), physics(Path(__file__))
            assert old_physics['worker'] == new_physics['worker'], 'Simulation worker changed'
            # Audited storage-only change: same founder bytes and returned
            # identities/sizes; symlinks avoid NTFS's 1,024-hard-link limit.
            allowed_installation_hashes = {
                '4878cced8263752c7e855cfdaaea05e1d37382572fd31e781e0f0447f4ee4dfa',
                '9a0b9d512494a4ead8bcd4da303b30cd2157f7fcbc3a5299c0759fd25d4826fd',
            }
            old_install = hashlib.sha256(old_physics['install_founder'].encode()).hexdigest()
            new_install = hashlib.sha256(new_physics['install_founder'].encode()).hexdigest()
            assert old_install == new_install or {old_install, new_install} <= allowed_installation_hashes, 'Founder installation changed'
        contract = saved
    else: base.save(manifest, contract)
    for name in ('pairs.tsv', 'positions.txt'):
        if (raw/name).exists(): assert base.digest(raw/name) == base.digest(previous/name)
        else: shutil.copy2(previous/name, raw/name)
    native = str(Path(configs['10k']['slim']).parent.parent/'lib')
    os.environ['LD_LIBRARY_PATH'] = native + (':'+os.environ['LD_LIBRARY_PATH'] if os.environ.get('LD_LIBRARY_PATH') else '')
    start = time.time(); deadline = min(start+args.seconds, args.deadline_epoch or float('inf'))
    base.save(raw/'schedule.json', dict(started_epoch=start, deadline_epoch=deadline, workers=args.workers,
        seconds=args.seconds, original_run_started=json.loads((previous/'schedule.json').read_text())['started_epoch'],
        script_sha256=base.digest(__file__), previous=str(previous)))
    shutil.copy2(previous/'execution_environment.json', raw/'execution_environment.json')
    snapshots=raw/'executed_sources'; snapshots.mkdir(exist_ok=True)
    for path in map(Path, (__file__, arrays.__file__, arrays.ancestor.__file__, base.__file__)):
        dest=snapshots/f'{path.stem}_{base.digest(path)[:12]}.py'
        if not dest.exists(): shutil.copy2(path,dest)
    tasks=[]
    for ordinal in range(target):
        for s in strata:
            slot=len(reused[s['key']])+ordinal
            if slot<target:
                limit=args.older_strong_attempts_per_slot if s['selection']==.05 and s['start']==50000 else args.attempts_per_slot
                tasks.append((s,slot,limit))
    completed={s['key']:{} for s in strata}; pending={}; cursor=0; last=0

    def checkpoint(final=False):
        status=dict(elapsed_seconds=time.time()-start, deadline_epoch=deadline, pending=len(pending), strata=[])
        ledger={}
        for s in strata:
            key=s['key']; rows=[]; slots=completed[key]; first=len(reused[key])
            for slot in range(first,target):
                if slot not in slots or not slots[slot]['complete']: break
                rows.extend(slots[slot]['records'])
            ledger[key]=rows
            kept=reused[key]+[r for r in rows if r['accepted']]
            # Include all in-flight receipts in monitoring, without using them
            # to select the reportable cohort by completion speed.
            active_status=[]
            for path in (raw/key).glob('slot_*/status.json') if (raw/key).exists() else []:
                try: active_status.append(json.loads(path.read_text()))
                except (FileNotFoundError, json.JSONDecodeError):
                    # WSL/Windows can briefly hide a concurrently replaced
                    # monitoring file. Durable simulation receipts are separate.
                    continue
            status['strata'].append(dict(s,reused=len(reused[key]),new_completed=len(rows),accepted=len(kept),target=target,
                lost=sum(r['rejection']=='loss' for r in rows),fixed=sum(r['rejection']=='fixation' for r in rows),
                monitored_new_proposals=sum(r['proposals'] for r in active_status),
                monitored_new_survivors=sum(r['accepted'] is not None for r in active_status),
                slots_complete=sum(r['complete'] for r in active_status)))
            if final:
                assert len({r['onset_sha256'] for r in kept})==len(kept)
                base.save(raw/(key+'_cohort.json'),dict(stratum=s,accepted=kept,target=target,complete=len(kept)==target,
                    reused_count=len(reused[key]),new_proposals=len(rows),manifest_sha256=base.digest(manifest)))
        arrays.save_ledger(raw/'proposal_ledger.json',ledger)
        arrays.save_ledger(raw/('finished.json' if final else 'progress.json'),status)
        print(json.dumps(status),flush=True)

    pool=ProcessPoolExecutor(max_workers=args.workers,initializer=os.setsid)
    try:
        while time.time()<deadline:
            while len(pending)<args.workers and cursor<len(tasks) and time.time()<deadline-180:
                s,slot,limit=tasks[cursor];cursor+=1
                future=pool.submit(worker,str(raw),configs[s['regime']],s,slot,limit,deadline)
                pending[future]=(s,slot)
            if not pending: break
            done,_=wait(pending,timeout=min(10,max(0,deadline-time.time())),return_when=FIRST_COMPLETED)
            for future in done:
                s,slot=pending.pop(future);completed[s['key']][slot]=future.result()
            if time.time()-last>=60:checkpoint();last=time.time()
    finally:
        if pending:
            for process in pool._processes.values():
                try:os.killpg(process.pid,signal.SIGKILL)
                except ProcessLookupError:pass
            pool.shutdown(wait=True,cancel_futures=True)
        else:pool.shutdown(wait=True)
        checkpoint(final=True)


if __name__=='__main__':main()
