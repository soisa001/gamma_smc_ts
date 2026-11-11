#!/usr/bin/env python3
"""Bounded de novo selected and segregation-conditioned neutral cohorts.

Uses the existing EAS simulation and decoder unchanged. A de novo allele starts
as one copy at selection onset. Neutral references match mutation age (10/50 kya).
"""
import os
os.environ['MPLBACKEND'] = 'Agg'
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'
import argparse
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
import importlib.metadata
import json
from pathlib import Path
import shutil
import signal
import sys
import time
import numpy as np
import pandas as pd
import run_selection_arrays as arrays
from report_selection_arrays import local_path

base = arrays.base
ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / 'docs/results/selection_calls_100_segregating_null_20261110'
HISTORY = ROOT / 'docs/results/selection_arrays_eas_20261105_4h'


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--raw', type=Path, default=ROOT.parent / 'sim/denovo_comparison_20261111')
    p.add_argument('--selection', type=float, nargs='+', default=[.005, .01, .05])
    p.add_argument('--target', type=int, default=100)
    p.add_argument('--workers', type=int, default=4)
    p.add_argument('--seconds', type=int, default=14400)
    p.add_argument('--attempts', type=int, default=50000)
    p.add_argument('--drain-seconds', type=int, default=180)
    args = p.parse_args()
    assert 1 <= args.workers <= 24 and args.target > 0 and args.attempts > 0
    assert len(set(args.selection)) == len(args.selection) and all(s > 0 and np.isfinite(s) for s in args.selection)
    raw = args.raw.resolve(); raw.mkdir(parents=True, exist_ok=True)
    reference = json.loads((REFERENCE / 'manifest.json').read_text())
    template = reference['simulation_configurations']['50k_s001']
    regimes = [('10k', 10000, 0), ('50k', 50000, 0), ('50_30k', 50000, 30000)]
    configs = {r: dict(template, selection_start_years=start, selection_end_years=end,
                      seed_namespace=f'denovo-comparison-v1:{r}:',
                      null_condition='Original single-copy de novo focal allele segregating in the present-day population')
               for r, start, end in regimes}
    strata = [dict(key=f'{r}_s{s:g}', regime=r, selection=s, start=start, end=end, role='selected')
              for r, start, end in regimes for s in args.selection]
    strata += [dict(key=f'{r}_neutral', regime=r, selection=0., start=start, end=0, role='neutral')
               for r, start in [('10k', 10000), ('50k', 50000)]]
    for key, path in [('ancestry_source_sha256', arrays.ancestor.__file__), ('score_source_sha256', base.__file__),
                      ('slim_sha256', template['slim']), ('decoder_sha256', template['decoder']),
                      ('phlash_sha256', ROOT / 'no_introgression_EAS_sim/resources/EAS.npz')]:
        assert base.digest(path) == template[key], str(path)
    native = str(Path(template['slim']).parent.parent / 'lib')
    os.environ['LD_LIBRARY_PATH'] = native + (':' + os.environ['LD_LIBRARY_PATH'] if os.environ.get('LD_LIBRARY_PATH') else '')
    inventory = json.loads((HISTORY / 'output_hashes.json').read_text())
    assert base.digest(HISTORY / 'completed_units.tsv') == inventory['completed_units.tsv']
    table = pd.read_csv(HISTORY / 'completed_units.tsv', sep='\t')
    old_configs = json.loads((HISTORY / 'manifest.json').read_text())['contract']['configs']
    reused = {s['key']: [] for s in strata}
    for s in strata:
        rows = table[(table.origin == 'de_novo') & (table.regime == s['regime']) &
                     (table.selection == s['selection']) & (table.final_frequency > 0) & (table.final_frequency < 1)]
        if s['role'] == 'selected': rows = rows[rows.cohort == 'selected']
        else: rows = rows[rows.cohort != 'selected']
        rows = rows.sort_values(['cohort', 'replicate']).iloc[:args.target]
        if len(rows):
            # Reuse only the same physics, with unchanged ancestry/score code.
            for key in ('length', 'diploids', 'pairs', 'stride', 'generation_time', 'dominance',
                        'scaling', 'mutation_rate', 'recombination_rate', 'base_seed', 'pair_seed',
                        'slim_sha256', 'decoder_sha256', 'phlash_sha256', 'ancestry_source_sha256',
                        'score_source_sha256', 'selection_start_years', 'selection_end_years'):
                assert old_configs[s['regime']][key] == configs[s['regime']][key], key
        for record in json.loads(rows.to_json(orient='records', double_precision=15)):
            directory = local_path(record['source_directory'])
            arrays.ancestor.verify(directory, record)
            assert base.digest(directory / 'trajectory.tsv') == record['trajectory_sha256']
            trajectory = np.loadtxt(directory / 'trajectory.tsv', skiprows=1)
            assert trajectory[0, 4] == 1 and 0 < trajectory[-1, 1] < 1
            reused[s['key']].append(record)
    contract = dict(configs=configs, strata=strata, target=args.target, reused=reused,
        origin='de_novo', de_novo_initial_copies=1, introgression=False,
        conditioning='0 < whole-population AF of the original focal allele < 1 at present, for selected and neutral.',
        neutral_mapping={'10k': '10k_neutral', '50k': '50k_neutral', '50_30k': '50k_neutral'},
        selection_rule='Frozen compatible survivors first, then new survivors in proposal-index order; no signal/AF ascertainment beyond segregation.',
        neutral_gate='Existing selected-cohort absorbing-state rejection gate with s=0; no neutral allele relabeling.',
        source_table_sha256=inventory['completed_units.tsv'],
        scientific_source_hashes={Path(f).name: base.digest(f) for f in (arrays.__file__, arrays.ancestor.__file__, base.__file__)})
    manifest = raw / 'manifest.json'
    if manifest.exists(): assert json.loads(manifest.read_text()) == contract
    else: base.save(manifest, contract)
    with np.load(REFERENCE / 'neutral.npz') as data: pairs = data['pairs'].copy()
    for name, values in [('pairs.tsv', pairs), ('positions.txt', np.arange(100) * 10000)]:
        path = raw / name
        if path.exists(): np.testing.assert_array_equal(np.loadtxt(path, dtype=int), values)
        else: np.savetxt(path, values, fmt='%d', delimiter='\t')
    snapshots = raw / 'executed_sources'; snapshots.mkdir(exist_ok=True)
    for f in map(Path, (__file__, arrays.__file__, arrays.ancestor.__file__, base.__file__)):
        dest = snapshots / f'{f.stem}_{base.digest(f)[:12]}.py'
        if not dest.exists(): shutil.copy2(f, dest)
    start = time.time(); deadline = start + args.seconds
    if (raw / 'schedule.json').exists():
        previous = json.loads((raw / 'schedule.json').read_text())
        shutil.copy2(raw / 'schedule.json', raw / f"schedule_{int(previous['started_epoch'])}.json")
    base.save(raw / 'schedule.json', dict(started_epoch=start, deadline_epoch=deadline, workers=args.workers,
        seconds=args.seconds, attempts_per_treatment=args.attempts, drain_seconds=args.drain_seconds,
        script_sha256=base.digest(__file__), raw=str(raw)))
    base.save(raw / 'execution_environment.json', dict(python=sys.version, workers=args.workers, threads_per_worker=1,
        native_library_path=native, versions={n: importlib.metadata.version(n) for n in
            ('numpy', 'scipy', 'msprime', 'pyslim', 'tskit', 'matplotlib', 'pandas', 'zstandard')}))
    completed = {s['key']: {} for s in strata}
    next_index = {s['key']: 0 for s in strata}
    pending = {}; cursor = 0; last_write = 0

    def prefix(s):
        rows = completed[s['key']]; n = 0
        while n in rows: n += 1
        return [rows[i] for i in range(n)]

    def accepted(s):
        return len(reused[s['key']]) + sum(r['accepted'] for r in completed[s['key']].values())

    def checkpoint(final=False):
        status = dict(elapsed_seconds=time.time()-start, deadline_epoch=deadline, pending=len(pending), strata=[])
        ledger = {}
        for s in strata:
            rows = prefix(s)
            ledger[s['key']] = rows
            kept = (reused[s['key']] + [r for r in rows if r['accepted']])[:args.target]
            status['strata'].append(dict(s, reused=len(reused[s['key']]), new_completed=len(rows),
                accepted=len(kept), target=args.target, submitted=next_index[s['key']],
                lost=sum(r['rejection']=='loss' for r in rows), fixed=sum(r['rejection']=='fixation' for r in rows)))
            if final:
                base.save(raw / (s['key'] + '_cohort.json'), dict(stratum=s, accepted=kept, target=args.target,
                    complete=len(kept)==args.target, reused_count=len(reused[s['key']]), new_proposals=len(rows),
                    manifest_sha256=base.digest(manifest)))
        arrays.save_ledger(raw / 'proposal_ledger.json', ledger)
        arrays.save_ledger(raw / ('finished.json' if final else 'progress.json'), status)
        print(json.dumps(status), flush=True)

    pool = ProcessPoolExecutor(max_workers=args.workers, initializer=os.setsid)
    try:
        while time.time() < deadline:
            while len(pending) < args.workers and time.time() < deadline - args.drain_seconds:
                chosen = None
                for _ in strata:
                    s = strata[cursor % len(strata)]; cursor += 1
                    if accepted(s) < args.target and next_index[s['key']] < args.attempts:
                        chosen = s; break
                if chosen is None: break
                s = chosen; i = next_index[s['key']]; next_index[s['key']] += 1
                suffix = f"s{s['selection']:g}/selected" if s['role']=='selected' else 'neutral_conditioned'
                unit = f'de_novo/{suffix}_{i:05d}'
                # cohort='selected' activates exact rejection even for neutral s=0.
                task = dict(id=f"{s['regime']}/{unit}", unit=unit, regime=s['regime'], origin='de_novo',
                            cohort='selected', scientific_role=s['role'], selection=s['selection'], replicate=i)
                future = pool.submit(arrays.work, (str(raw), configs[s['regime']], task, deadline))
                pending[future] = (s, task)
            if not pending: break
            done, _ = wait(pending, timeout=min(10, max(0, deadline-time.time())), return_when=FIRST_COMPLETED)
            for future in done:
                s, task = pending.pop(future); result = future.result()
                assert result['state'] in ('accepted', 'rejected'), result
                path = raw / task['id'] / ('complete.json' if result['state']=='accepted' else 'rejected.json')
                record = json.loads(path.read_text()); record['source_directory'] = str(path.parent)
                completed[s['key']][task['replicate']] = record
            if time.time()-last_write >= 60:
                checkpoint(); last_write = time.time()
    finally:
        if pending:
            for process in pool._processes.values():
                try: os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError: pass
            pool.shutdown(wait=True, cancel_futures=True)
        else: pool.shutdown(wait=True)
        checkpoint(final=True)


if __name__ == '__main__': main()
