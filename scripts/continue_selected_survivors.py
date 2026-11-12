#!/usr/bin/env python3
"""Continue incomplete selected cohorts with unchanged seeds and exact rejection.

The immutable manifest and dispatch journal identify every proposal. A restart
replays durable receipts and preserves the original deadline. De novo proposals
reuse the independently drawn founders from the previous run, with at most one
retained history per founder. Introgression draws independent ancestry each time.
"""
import os
os.environ['MPLBACKEND'] = 'Agg'
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
import hashlib
import json
import math
from pathlib import Path
import shutil
import signal
import sys
import time
import numpy as np
import run_selection_arrays as arrays
import run_denovo_cached_founders as cached

ROOT = Path(__file__).resolve().parents[1]
BASE = arrays.base


def job(raw, cfg, task, deadline, founder):
    original = arrays.ancestor.onset
    if founder:
        arrays.ancestor.onset = lambda c, o, l, d: cached.install_founder(Path(founder), c, o, l, d)
    try:
        return arrays.work((raw, cfg, task, deadline))
    finally:
        arrays.ancestor.onset = original


def contract(raw, workers, memory_gb):
    sources = {'de_novo': ROOT.parent/'sim/denovo_comparison_cached_20261111',
               'introgression': ROOT.parent/'sim/strong_selection_20261110'}
    strata = []
    for origin, source in sources.items():
        old = json.loads((source/'manifest.json').read_text())
        finished = json.loads((source/'finished.json').read_text())
        hashes = old.get('scientific_source_hashes', old.get('source_hashes'))
        for name, expected in hashes.items():
            assert BASE.digest(ROOT/'scripts'/name) == expected, name
        for entry in finished['strata']:
            if entry.get('role', 'selected') != 'selected' or entry['accepted'] >= 100:
                continue
            key = entry['key']; receipt = source/(key+'_cohort.json')
            cohort = json.loads(receipt.read_text())
            assert len(cohort['accepted']) == entry['accepted']
            cfg = old['configs'][entry['regime']]
            for name, path in [('slim_sha256', cfg['slim']), ('decoder_sha256', cfg['decoder']),
                               ('phlash_sha256', ROOT/'no_introgression_EAS_sim/resources/EAS.npz')]:
                assert BASE.digest(path) == cfg[name]
            for r in cohort['accepted']:
                assert 0 < r['final_frequency'] < 1
                arrays.ancestor.verify(Path(r['source_directory']), r)
            s = dict(key=origin+'__'+key, old_key=key, origin=origin, regime=entry['regime'],
                selection=entry['selection'], start=entry['start'], end=entry['end'], config=cfg,
                previous=str(source), previous_manifest_sha256=BASE.digest(source/'manifest.json'),
                previous_cohort_sha256=BASE.digest(receipt), reused=cohort['accepted'],
                previous_proposals=entry['new_completed'], previous_lost=entry['lost'], previous_fixed=entry['fixed'])
            if origin == 'de_novo':
                assert not s['reused'] and entry['slots_complete'] == 100
                s['slot_attempt_starts'] = []
                s['founders'] = []
                for slot in range(100):
                    directory = source/key/f'slot_{slot:03d}'
                    status = json.loads((directory/'status.json').read_text())
                    assert status['complete'] and status['accepted'] is None
                    founder = directory/'founder'
                    fm = json.loads((founder/'founder.json').read_text())
                    assert fm['config_sha256'] == s['previous_manifest_sha256']
                    assert BASE.digest(founder/'onset.trees') == fm['onset_sha256']
                    s['slot_attempt_starts'].append(status['proposals'])
                    s['founders'].append(dict(path=str(founder), **fm))
            else:
                assert entry['submitted'] == entry['new_completed']
                s['proposal_start'] = entry['submitted']
            strata.append(s)
    value = dict(schema='selected-survivor-continuation/v1', strata=strata, target=100,
        workers=workers, memory_budget_gb=memory_gb, threads_per_worker=1,
        conditioning='Original focal allele: 0 < whole-population AF < 1 today; fixed and lost alleles excluded.',
        selection_rule='Existing survivors followed by successes in the completed dispatch prefix within each stratum. At most one de novo survivor per independent founder.',
        deadline_rule='Five-hour wall clock including inference; stop dispatch 180 seconds before deadline, terminate remaining worker process groups at deadline.',
        sources={Path(p).name: BASE.digest(p) for p in (__file__, arrays.__file__, arrays.ancestor.__file__, BASE.__file__, cached.__file__)},
        note='Only incomplete cohorts in the current 10k, 50k, and 50-30k series; completed selected cohorts and neutral references are frozen.')
    path = raw/'manifest.json'
    if path.exists(): assert json.loads(path.read_text()) == value, 'Manifest mismatch'
    else: BASE.save(path, value)
    source = sources['de_novo']
    for name in ('pairs.tsv', 'positions.txt'):
        target = raw/name
        if target.exists(): assert BASE.digest(target) == BASE.digest(source/name)
        else: shutil.copy2(source/name, target)
    snapshots = raw/'executed_sources'; snapshots.mkdir(exist_ok=True)
    for filename, digest in value['sources'].items():
        target = snapshots/(Path(filename).stem+'_'+digest[:12]+'.py')
        if not target.exists(): shutil.copy2(ROOT/'scripts'/filename, target)
    return value


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--raw', type=Path, default=ROOT.parent/'sim/selected_survivor_extension_20261112')
    p.add_argument('--workers', type=int, default=4)
    p.add_argument('--seconds', type=int, default=18000)
    p.add_argument('--memory-gb', type=float, default=26)
    p.add_argument('--max-new-proposals', type=int, default=1000000, help='Per treatment, also bounded by time and target.')
    p.add_argument('--prepare-only', action='store_true')
    args = p.parse_args()
    assert 1 <= args.workers <= 20 and 0 < args.seconds <= 18000 and 0 < args.memory_gb <= 200
    raw = args.raw.resolve(); raw.mkdir(parents=True, exist_ok=True)
    spec = contract(raw, args.workers, args.memory_gb)
    signature = BASE.digest(raw/'manifest.json')
    states = {}
    for s in spec['strata']:
        states[s['key']] = dict(spec=s, next_ordinal=0, prefix=0, waiting={}, kept=list(s['reused']),
            lost=0, fixed=0, active=0, completed=0, all_survivors=0, busy_slots=set(), successful_slots=set(),
            attempts=list(s.get('slot_attempt_starts', [])), slot_cursor=0,
            proposal=s.get('proposal_start', 0))
    if args.prepare_only:
        print(json.dumps({'prepared': True, 'strata': [{k:s[k] for k in ('key','selection','start','end')} for s in spec['strata']]}))
        return
    lib = str(Path(spec['strata'][0]['config']['slim']).parent.parent/'lib')
    os.environ['LD_LIBRARY_PATH'] = lib + (':'+os.environ['LD_LIBRARY_PATH'] if os.environ.get('LD_LIBRARY_PATH') else '')
    schedule_path = raw/'schedule.json'
    if schedule_path.exists():
        schedule = json.loads(schedule_path.read_text()); deadline = schedule['deadline_epoch']
    else:
        start = time.time(); deadline = start+args.seconds
        schedule = dict(started_epoch=start, deadline_epoch=deadline, workers=args.workers, memory_budget_gb=args.memory_gb,
            seconds=args.seconds, max_new_proposals_per_treatment=args.max_new_proposals,
            native_library_path=lib, python=sys.version, manifest_sha256=signature)
        BASE.save(schedule_path, schedule)
    assert schedule['manifest_sha256'] == signature
    history_path = raw/'dispatch.jsonl'
    replay = []
    verified_onsets = {}

    def consume(event, verify=False):
        st = states[event['key']]; task = event['task']; directory = raw/task['id']
        receipt = next((directory/n for n in ('complete.json','rejected.json') if (directory/n).exists()), None)
        if receipt is None: return False
        record = json.loads(receipt.read_text())
        assert record['config_sha256'] == signature and record['seed'] == BASE.seed(st['spec']['config']['seed_namespace']+task['unit'])
        assert record['selection'] == st['spec']['selection'] and record['origin'] == st['spec']['origin']
        if verify:
            assert BASE.digest(directory/'trajectory.tsv') == record['trajectory_sha256']
            onset = (directory/'onset.trees').resolve()
            if onset not in verified_onsets: verified_onsets[onset] = BASE.digest(onset)
            assert verified_onsets[onset] == record['onset_sha256']
            if record['accepted']: arrays.ancestor.verify(directory, record)
        if record['accepted']:
            assert 0 < record['final_frequency'] < 1 and record['stop_years_ago'] == 0
            record['source_directory'] = str(directory)
            st['all_survivors'] += 1
            if event['founder']: st['successful_slots'].add(task['slot'])
        else:
            assert record['rejection'] in ('loss','fixation')
            assert record['final_frequency'] == (0 if record['rejection']=='loss' else 1)
        st['completed'] += 1
        st['waiting'][event['ordinal']] = record
        while st['prefix'] in st['waiting']:
            r = st['waiting'].pop(st['prefix']); st['prefix'] += 1
            if r['accepted']:
                if len(st['kept']) < spec['target']: st['kept'].append(r)
            elif r['rejection']=='loss': st['lost'] += 1
            else: st['fixed'] += 1
        return True

    if history_path.exists():
        for line in history_path.read_text().splitlines():
            event = json.loads(line); st = states[event['key']]; task = event['task']
            assert event['ordinal'] == st['next_ordinal']; st['next_ordinal'] += 1
            if event['founder']:
                st['attempts'][task['slot']] = max(st['attempts'][task['slot']], task['attempt']+1)
            else: st['proposal'] = max(st['proposal'], task['replicate']+1)
            if not consume(event, verify=True):
                replay.append(event)
                if event['founder']: st['busy_slots'].add(task['slot'])

    def new_event(st):
        s = st['spec']; founder = None
        if len(st['kept']) >= spec['target'] or st['next_ordinal'] >= args.max_new_proposals: return None
        if s['origin']=='de_novo':
            slot = None
            for _ in range(100):
                candidate = st['slot_cursor'] % 100; st['slot_cursor'] += 1
                if candidate not in st['successful_slots'] and candidate not in st['busy_slots']:
                    slot = candidate; break
            if slot is None: return None
            attempt = st['attempts'][slot]; st['attempts'][slot] += 1
            unit = f"de_novo/{s['old_key']}/slot_{slot:03d}/attempt_{attempt:05d}"
            fm = s['founders'][slot]; founder = fm['path']
            task = dict(id=unit, unit=unit, regime=s['regime'], origin='de_novo', cohort='selected', scientific_role='selected',
                selection=s['selection'], replicate=slot*1000000+attempt, slot=slot, attempt=attempt,
                founder_seed_label=fm['seed_label'], founder_ancestry_seed=fm['ancestry_seed'], founder_pairing_seed=fm['pairing_seed'])
            st['busy_slots'].add(slot)
        else:
            i = st['proposal']; st['proposal'] += 1
            unit = f"introgression/s{s['selection']:g}/selected_{i:05d}"
            task = dict(id=f"{s['regime']}/{unit}", unit=unit, regime=s['regime'], origin='introgression', cohort='selected', selection=s['selection'], replicate=i)
        event = dict(key=s['key'], ordinal=st['next_ordinal'], task=task, founder=founder)
        st['next_ordinal'] += 1
        return event

    pending = {}; last = 0; cursor = 0; replay_cursor = 0
    keys = list(states); share = math.ceil(args.workers/len(keys))
    def checkpoint(final=False):
        status = dict(elapsed_seconds=time.time()-schedule['started_epoch'], deadline_epoch=deadline,
            pending=len(pending), manifest_sha256=signature, strata=[])
        for key, st in states.items():
            s = st['spec']
            status['strata'].append(dict(key=key, origin=s['origin'], regime=s['regime'], selection=s['selection'], start=s['start'], end=s['end'],
                previous_accepted=len(s['reused']), accepted=len(st['kept']), target=spec['target'],
                new_submitted=st['next_ordinal'], new_completed=st['completed'], completed_prefix=st['prefix'],
                prefix_lost=st['lost'], prefix_fixed=st['fixed'], all_completed_survivors=st['all_survivors'],
                previous_proposals=s['previous_proposals'], previous_lost=s['previous_lost'], previous_fixed=s['previous_fixed']))
            if final:
                assert len({r['observed_sha256'] for r in st['kept']}) == len(st['kept'])
                assert len({r['onset_sha256'] for r in st['kept']}) == len(st['kept'])
                BASE.save(raw/(key+'_cohort.json'),dict(stratum={k:s[k] for k in ('key','origin','regime','selection','start','end')},
                    accepted=st['kept'], target=spec['target'], manifest_sha256=signature))
        arrays.save_ledger(raw/('finished.json' if final else 'progress.json'), status)
        print(json.dumps(status), flush=True)
    pool = ProcessPoolExecutor(max_workers=args.workers, initializer=os.setsid)
    try:
        with history_path.open('a', buffering=1) as journal:
            while time.time() < deadline:
                while len(pending) < args.workers and time.time() < deadline-180:
                    event = None
                    if replay_cursor < len(replay):
                        event = replay[replay_cursor]; replay_cursor += 1
                    else:
                        for _ in keys:
                            st = states[keys[cursor % len(keys)]]; cursor += 1
                            if st['active'] < share:
                                event = new_event(st)
                                if event is not None:
                                    journal.write(json.dumps(event)+'\n'); break
                    if event is None: break
                    st = states[event['key']]; st['active'] += 1
                    future = pool.submit(job, str(raw), st['spec']['config'], event['task'], deadline, event['founder'])
                    pending[future] = event
                if not pending: break
                done, _ = wait(pending, timeout=min(5,max(0,deadline-time.time())), return_when=FIRST_COMPLETED)
                for future in done:
                    event = pending.pop(future); st = states[event['key']]; st['active'] -= 1
                    result = future.result()
                    assert result['state'] in ('accepted','rejected'), result
                    assert consume(event)
                    if event['founder']: st['busy_slots'].remove(event['task']['slot'])
                if time.time()-last >= 60: checkpoint(); last = time.time()
    finally:
        if pending:
            for process in pool._processes.values():
                try: os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError: pass
            pool.shutdown(wait=True, cancel_futures=True)
            for event in pending.values(): consume(event, verify=True)
        else: pool.shutdown(wait=True)
        checkpoint(final=True)


if __name__ == '__main__': main()
