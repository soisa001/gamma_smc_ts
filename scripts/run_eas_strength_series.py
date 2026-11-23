"""Assemble and simulate the fixation-inclusive EAS strength series."""
import argparse
import csv
import json
import os
import shutil
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import run_cumulative_grid as grid
import replay_fixed_tag
import report_cumulative_grid as old

RAW = Path.home()/'sim/eas_strength_series_20261118'
OUT = old.ROOT/'docs/results/eas_strength_series_20261118'
SIM = old.ROOT.parent/'sim'
MID = old.ROOT/'docs/results/cumulative_grid_20261118'
STRENGTHS = tuple(i/1000 for i in range(1,11)) + (.02,.03)
arrays = grid.arrays


def key(origin, regime, s):
    return f"{origin}_{regime}_s{str(s).replace('.', 'p')}"


def make_plan(raw):
    raw.mkdir(parents=True, exist_ok=True)
    manifest = raw/'manifest.json'
    if manifest.exists():
        return old.read(manifest)
    baseline = old.read(old.OUT/'simulation_manifest.json')
    cfgs = baseline['configs']
    for c in cfgs.values():
        c['seed_namespace']='eas-strength-series-v1:'
        for field, path in [('slim_sha256', c['slim']), ('decoder_sha256', c['decoder']),
                            ('phlash_sha256',old.ROOT/'no_introgression_EAS_sim/resources/EAS.npz')]:
            assert old.sha(Path(path))==c[field]
    sources = defaultdict(list)
    seen = defaultdict(set)
    def register(k, directory, rec):
        if len(sources[k])>=100 or rec['onset_sha256'] in seen[k]:
            return
        if rec['accepted'] or rec['rejection']=='fixation':
            seen[k].add(rec['onset_sha256'])
            sources[k].append(dict(path=str(directory),fixed=not rec['accepted'],onset=rec['onset_sha256']))

    table = grid.ARCHIVE/'completed_units.tsv'
    with table.open() as f:
        recent = [r for r in csv.DictReader(f, delimiter='\t') if r['regime']=='10k' and r['cohort']=='selected'
                  and float(r['selection']) in (.02,.03) and r['included_prefix']=='True']
    for r in sorted(recent,key=lambda r:int(r['replicate'])):
        path=grid.local_path(r['source_directory'])
        assert old.sha(path/'complete.json')==r['receipt_sha256']
        register(key(r['origin'],'10k',float(r['selection'])),path,old.read(path/'complete.json'))
    # Take BOTH possible surviving outcomes, in original proposal order.
    def scan(bank, regime, origin, s, limit):
        k=key(origin,regime,s)
        folder=bank/regime/origin/f's{s:g}'
        for i in range(limit):
            if len(sources[k])==100: break
            path=folder/f'selected_{i:05d}'
            receipt=path/'complete.json'
            if not receipt.is_file(): receipt=path/'rejected.json'
            register(k,path,old.read(receipt))
    for s in (.005,.01,.02,.03):
        scan(SIM/'selection_arrays_20261104','50k','de_novo',s,1225)
        if s!=.005:
            scan(SIM/'selection_arrays_20261104','50k','introgression',s,1225 if s!=.01 else 325)
    scan(SIM/'selection_50_30k_demo_20261106','50_30k','introgression',.01,201)
    for name in ('denovo_comparison_20261111','denovo_comparison_cached_20261111'):
        ledger=old.read(SIM/name/'proposal_ledger.json')
        for regime,s in [('50k',.005),('50k',.01),('50_30k',.01)]:
            for r in ledger[f'{regime}_s{s:g}']:
                register(key('de_novo',regime,s),r['source_directory'],r)
    cases=[]
    for s in STRENGTHS:
        for regime,(start,end) in grid.WINDOWS.items():
            for origin in old.ORIGINS:
                k=key(origin,regime,s)
                frozen=None
                if s in (.002,.003,.004):
                    frozen=MID.name
                elif s in (.001,.005,.01) and k not in sources:
                    frozen=old.OUT.name
                cases.append(dict(key=k,origin=origin,regime=regime,selection=s,start=start,end=end,
                                  frozen=frozen,archive=sources.get(k,[])))
    provenance = {p.name:old.sha(p) for p in map(Path,[__file__,replay_fixed_tag.__file__,arrays.__file__,
                      grid.base.__file__,grid.cached.__file__,arrays.ancestor.__file__])}
    contract=dict(target=100,configs=cfgs,cases=cases,script_sha256=provenance,input_sha256=baseline['input_hashes'],
                  base_banks={p.name:old.sha(p/'mass.npz') for p in (old.OUT,MID)},
                  archive_list_sha256=old.sha(table),conditioning='First nonlost and scoreable replicate per slot, include fixes at the best surviving historical LD tag.')
    arrays.save_ledger(manifest,contract)
    return contract


def founder(raw, cfg, case, slot):
    cached=grid.cached
    directory=raw/'founders'/case['key']/f'{slot:03d}'
    directory.mkdir(parents=True,exist_ok=True)
    receipt=directory/'founder.json'
    if not receipt.exists():
        label=f"eas-full-series-founder:{case['key']}:{slot}"
        age=case['start']//cfg['generation_time']
        carriers,sizes=cached.ORIGINAL_ONSET(dict(cfg,ancestry_start_generations=age),'de_novo',label,directory)
        assert len(carriers)==1
        arrays.save_ledger(receipt,dict(carrier_ids=carriers,sizes=sizes,age_generations=age,label=label,
                                       onset_sha256=old.sha(directory/'onset.trees')))
    assert old.sha(directory/'onset.trees')==old.read(receipt)['onset_sha256']
    return directory


def slot_worker(raw, case, slot, cfg, signature):
    receipt=raw/'slots'/case['key']/f'{slot:03d}.json'
    if receipt.exists():
        rec=old.read(receipt)
        assert rec['manifest']==signature and old.sha(Path(rec['scores']))==rec['scores_sha256']
        return
    receipt.parent.mkdir(parents=True,exist_ok=True)
    founder_dir=None
    archived=case['archive'][slot] if slot<len(case['archive']) else None
    try:
        for attempt in range(100000):
            if archived and attempt==0:
                directory=Path(archived['path'])
                fixed=archived['fixed']
            else:
                if case['origin']=='de_novo' and founder_dir is None:
                    founder_dir=founder(raw,cfg,case,slot)
                    arrays.ancestor.onset=lambda *args:grid.cached.install_founder(founder_dir,*args)
                unit=f"{case['key']}/slot_{slot:03d}/attempt_{attempt:05d}"
                directory=raw/unit
                task=dict(id=unit,unit=unit,regime=case['regime'],origin=case['origin'],cohort='selected',
                          selection=case['selection'],replicate=slot*100000+attempt)
                state=arrays.work((raw,cfg,task,float('inf')))['state']
                fixed=state=='rejected'
            rec=old.read(directory/('rejected.json' if fixed else 'complete.json'))
            if archived and attempt==0: assert rec['onset_sha256']==archived['onset']
            assert (rec['origin'],rec['selection'])==(case['origin'],case['selection'])
            if fixed:
                if rec['rejection']=='loss': continue
                assert rec['rejection']=='fixation'
                tagged=raw/'tags'/case['key']/f'{slot:03d}'/f'{attempt:05d}'
                result=replay_fixed_tag.replay(directory,tagged,cfg,raw/'pairs.tsv',signature)
                tag=result['tag']
                if tag is None: continue
                scores=tagged/'scores.npz'
                digest=result['scores_sha256']
                position=tag['position']
            else:
                assert rec['accepted'] and rec['stop_years_ago']==0 and 0<rec['final_frequency']<1
                scores=directory/'scores.npz'
                digest=rec['scores_sha256']
                tag=None
                position=cfg['length']//2
            assert old.sha(scores)==digest
            arrays.save_ledger(receipt,dict(source=str(directory),scores=str(scores),scores_sha256=digest,
                        seed=rec['seed'],onset=rec['onset_sha256'],position=position,tag=tag,attempt=attempt,
                        final_frequency=rec['final_frequency'],manifest=signature))
            return
        raise RuntimeError(f"No usable history for {case['key']} slot {slot}")
    finally:
        arrays.ancestor.onset=grid.cached.ORIGINAL_ONSET


def simulate(raw, workers):
    contract=make_plan(raw)
    sig=old.sha(raw/'manifest.json')
    for name,digest in contract['input_sha256'].items():
        if not (raw/name).exists(): shutil.copyfile(SIM/'cumulative_grid_20261117'/name,raw/name)
        assert old.sha(raw/name)==digest
    cfgs=contract['configs']
    os.environ['LD_LIBRARY_PATH']=str(Path(cfgs['50k']['slim']).parent.parent/'lib')+':'+os.environ.get('LD_LIBRARY_PATH','')
    cases=[c for c in contract['cases'] if not c['frozen']]
    total=len(cases)*100
    print(json.dumps(dict(workers=workers,raw=str(raw),total=7200,frozen=7200-total,
          archived=sum(len(c['archive']) for c in cases),remaining_slots=total,started=time.time())),flush=True)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        jobs=[pool.submit(slot_worker,raw,c,i,cfgs[c['regime']],sig) for i in range(100) for c in cases]
        for i,j in enumerate(as_completed(jobs),1):
            j.result()
            if i%100==0: print(f'{i}/{total} completed at {time.time()}',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--workers',type=int,default=4,choices=range(1,21))
    p.add_argument('--raw',type=Path,default=RAW)
    args=p.parse_args()
    simulate(args.raw.resolve(),args.workers)
