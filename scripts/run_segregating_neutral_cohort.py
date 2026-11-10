#!/usr/bin/env python3
"""Fill a 100-replicate EAS introgression null conditioned on focal segregation.

run: reuse surviving members of the old neutral pool, then exact rejection at
s=0 with fixed proposal seeds. prepare: verify scores/trajectories and recalibrate
the unchanged selected cohorts. Raw simulations stay outside the repository.
"""
import os
os.environ['MPLBACKEND']='Agg'
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import argparse
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, wait, FIRST_COMPLETED
import json
import importlib.metadata
from pathlib import Path
import shutil
import sys
import time
import numpy as np
import pandas as pd
import run_selection_arrays as arrays
from report_selection_arrays import local_path
from pointwise_time_tests import exceedance_counts
from plot_cohort_significant_times import ROW_ORDER

base=arrays.base
ROOT=Path(__file__).resolve().parents[1]
OLD=ROOT/'docs/results/selection_calls_100_20261110'
NEUTRAL=ROOT/'docs/results/selected_site_posteriors_20261106'
TABLE=ROOT/'docs/results/selected_vs_100_neutrals_10k_s001_20261106'
CONDITION='0 < whole-population frequency of the original introduced focal allele at present < 1'


def verify(record):
    directory=local_path(record['source_directory'])
    for name,key in [('trajectory.tsv','trajectory_sha256'),('onset.trees','onset_sha256')]:
        assert base.digest(directory/name)==record[key],str(directory/name)
    if record['accepted']:
        arrays.ancestor.verify(directory,record)
        assert 0<record['final_frequency']<1 and record['stop_years_ago']==0


def work(payload):
    raw,cfg,signature,i=payload
    raw=Path(raw)
    unit=f'introgression/neutral_conditioned_{i:05d}'
    directory=raw/unit;directory.mkdir(parents=True,exist_ok=True)
    task=dict(id=cfg.get('cohort_id_prefix','segregating_null')+'/'+unit,unit=unit,origin='introgression',cohort='neutral_conditioned',
        selection=0.,replicate=i,regime='10k',source_directory=str(directory))
    for name in ('complete.json','rejected.json'):
        path=directory/name
        if path.exists():
            record=json.loads(path.read_text());assert record['config_sha256']==signature
            verify(record)
            return record
    sim=directory/'simulation.json'
    if sim.exists():
        record=json.loads(sim.read_text());assert record['config_sha256']==signature
        assert base.digest(directory/'observed.trees')==record['observed_sha256']
    else:
        record=dict(task,seed=base.seed(cfg['seed_namespace']+unit),config_sha256=signature,
            conditioning=CONDITION,simulation_gate='Use existing selected-cohort rejection gate with selection coefficient exactly zero; neutral physics, original focal allele.')
        # The model's selected flag only enables the absorbing-state rejection
        # gate and suppresses neutral-site relabelling. s stays zero throughout.
        record.update(arrays.simulate(cfg,dict(task,cohort='selected'),directory))
        trajectory=np.loadtxt(directory/'trajectory.tsv',skiprows=1,ndmin=2)
        np.testing.assert_array_equal(trajectory[:,3],0)
        assert record['accepted']==bool(0<trajectory[-1,1]<1)
        if not record['accepted']:
            base.save(directory/'rejected.json',record)
            return record
        base.save(sim,record)
    record.update(base.score(cfg,directory,raw))
    base.save(directory/'complete.json',record)
    return record


def sources():
    old=json.loads((OLD/'manifest.json').read_text())
    inventory=json.loads((TABLE/'output_hashes.json').read_text())
    assert base.digest(TABLE/'cohort_manifest.tsv')==inventory['cohort_manifest.tsv']
    rows=pd.read_csv(TABLE/'cohort_manifest.tsv',sep='\t').iloc[1:]
    with np.load(OLD/'neutral.npz') as data:
        np.testing.assert_array_equal(rows.id.to_numpy(str),data['ids'])
        pairs=data['pairs'].copy()
    assert len(rows)==100 and rows.id.nunique()==100 and np.all(rows.selection==0)
    survivors=rows[(rows.final_frequency>0)&(rows.final_frequency<1)]
    reused=json.loads(survivors.to_json(orient='records',double_precision=15))
    for record in reused:verify(record)
    return old,pairs,reused


def run(args):
    old,pairs,reused=sources()
    if args.independent_test:reused=[]
    raw=args.raw;raw.mkdir(parents=True,exist_ok=True)
    cfg=dict(old['simulation_configurations']['10k_s001'],seed_namespace='neutral-segregating-v1:',
        null_condition=CONDITION)
    if args.independent_test:
        cfg.update(seed_namespace='neutral-segregating-independent-test-v1:',cohort_id_prefix='segregating_null_test')
    for key,path in [('ancestry_source_sha256',arrays.ancestor.__file__),('score_source_sha256',base.__file__),
                     ('slim_sha256',cfg['slim']),('decoder_sha256',cfg['decoder']),
                     ('phlash_sha256',ROOT/'no_introgression_EAS_sim/resources/EAS.npz')]:
        assert base.digest(path)==cfg[key],str(path)
    native_lib=Path(cfg['slim']).parent.parent/'lib'
    os.environ['LD_LIBRARY_PATH']=str(native_lib)+(':'+os.environ['LD_LIBRARY_PATH'] if os.environ.get('LD_LIBRARY_PATH') else '')
    contract=dict(config=cfg,condition=CONDITION,target=args.target,reused_ids=[r['id'] for r in reused],
        parent_manifest_sha256=base.digest(OLD/'manifest.json'),source_table_sha256=base.digest(TABLE/'cohort_manifest.tsv'),
        scientific_sources={Path(p).name:base.digest(p) for p in (arrays.__file__,arrays.ancestor.__file__,base.__file__)},
        selection='All surviving members of the original 100-neutral pool, then first accepted new proposals in index order; no selection by AF beyond segregation or inferred signal.',
        seed_scheme='Existing base seed 20380101 and deterministic phase hashes; unique namespace neutral-segregating-v1:.')
    if args.independent_test:
        contract.update(selection='First 100 accepted new proposals in index order; no reuse from reference cohort and no ascertainment by AF beyond segregation or by inferred signal.',
            seed_scheme='Existing base seed 20380101 and deterministic phase hashes; independent namespace neutral-segregating-independent-test-v1:.',
            independent_test=True)
    manifest=raw/'manifest.json'
    if manifest.exists():assert json.loads(manifest.read_text())==contract
    else:base.save(manifest,contract)
    for name,values in [('pairs.tsv',pairs),('positions.txt',np.arange(100)*10000)]:
        path=raw/name
        if path.exists():np.testing.assert_array_equal(np.loadtxt(path,dtype=int),values)
        else:np.savetxt(path,values,fmt='%d',delimiter='\t')
    snapshots=raw/'executed_sources';snapshots.mkdir(exist_ok=True)
    for p in (Path(__file__),Path(arrays.__file__),Path(arrays.ancestor.__file__),Path(base.__file__)):
        dest=snapshots/f'{p.stem}_{base.digest(p)[:12]}.py'
        if not dest.exists():shutil.copy2(p,dest)
    started=time.time();deadline=min(started+args.seconds,args.deadline_epoch or float('inf'))
    if (raw/'schedule.json').exists():
        prior=json.loads((raw/'schedule.json').read_text())
        shutil.copy2(raw/'schedule.json',raw/f"schedule_{int(prior['started_epoch'])}.json")
    base.save(raw/'schedule.json',dict(workers=args.workers,seconds=args.seconds,attempts=args.attempts,
        target=args.target,started_epoch=started,deadline_epoch=deadline,raw=str(raw),out=str(args.out),
        script_sha256=base.digest(__file__),native_library_path=str(native_lib)))
    base.save(raw/'execution_environment.json',dict(python=sys.version,threads_per_worker=1,workers=args.workers,
        native_library_path=str(native_lib),packages={name:importlib.metadata.version(name) for name in ('numpy','scipy','msprime','pyslim','tskit','matplotlib','pandas','zstandard')}))
    records=[];accepted=list(reused);next_index=0;pending={}
    print(f'Reused {len(reused)} segregating neutrals; target {args.target}.',flush=True)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        while True:
            while len(pending)<args.workers and next_index<args.attempts and len(accepted)<args.target and time.time()<deadline:
                future=pool.submit(work,(str(raw),cfg,base.digest(manifest),next_index))
                pending[future]=next_index;next_index+=1
            if not pending:break
            done,_=wait(pending,return_when=FIRST_COMPLETED)
            for future in done:
                records.append(future.result());pending.pop(future)
            records.sort(key=lambda r:r['replicate'])
            accepted=list(reused)+[r for r in records if r['accepted']]
            base.save(raw/'proposal_ledger.json',records)
            progress=dict(target=args.target,reused=len(reused),new_completed=len(records),accepted=len(accepted),
                new_accepted=sum(r['accepted'] for r in records),lost=sum(r['rejection']=='loss' for r in records),
                fixed=sum(r['rejection']=='fixation' for r in records),elapsed_seconds=time.time()-started,
                deadline_epoch=deadline,last_proposal=records[-1]['replicate'],submitted=next_index,pending=len(pending))
            base.save(raw/'progress.json',progress)
            print(json.dumps(progress),flush=True)
    assert len(accepted)>=args.target,'Neutral target not reached within budget; resume using verified caches.'
    # Finish every submitted proposal before retaining the first 100 accepted
    # IDs. Thus completion speed cannot change which simulations are retained.
    assert [r['replicate'] for r in records]==list(range(next_index))
    assert len({r['id'] for r in accepted[:args.target]})==args.target
    base.save(raw/'cohort.json',dict(target=args.target,accepted=accepted[:args.target],reused=reused,
        new_proposals=len(records),new_accepted=sum(r['accepted'] for r in records),manifest_sha256=base.digest(manifest),
        script_sha256=base.digest(__file__)))


def prepare(args):
    old,pairs,reused=sources()
    cohort=json.loads((args.raw/'cohort.json').read_text())
    assert cohort['target']==100 and len(cohort['accepted'])==100
    assert cohort['manifest_sha256']==base.digest(args.raw/'manifest.json')
    def read(record):
        verify(record)
        directory=local_path(record['source_directory'])
        trajectory=np.loadtxt(directory/'trajectory.tsv',skiprows=1)
        assert trajectory.shape==(2001,5)
        np.testing.assert_array_equal(trajectory[:,0],np.arange(50000,-1,-25))
        np.testing.assert_array_equal(trajectory[:,3],0)
        np.testing.assert_allclose(trajectory[:,1],trajectory[:,4]/(2*trajectory[:,2]),atol=1e-14)
        assert 0<trajectory[-1,1]<1
        meta=json.loads((directory/'posterior.zst.meta').read_text())
        np.testing.assert_array_equal(meta['pairs'],pairs)
        with np.load(directory/'scores.npz') as data:
            np.testing.assert_array_equal(data['positions'],np.arange(100)*10000)
            np.testing.assert_array_equal(data['cutoffs'],np.arange(1000,50001,1000))
            np.testing.assert_array_equal(data['valid_pairs'],1000)
            return trajectory,data['cdfs'][:,50,:].copy()
    with ThreadPoolExecutor(max_workers=args.workers) as pool:parts=list(pool.map(read,cohort['accepted']))
    trajectories=np.stack([p[0] for p in parts]);cdfs=np.stack([p[1] for p in parts])
    ids=np.array([r['id'] for r in cohort['accepted']])
    assert len(set(ids))==100 and len({r['observed_sha256'] for r in cohort['accepted']})==100
    out=args.out;out.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(out/'neutral.npz',cdfs=cdfs,ids=ids,pairs=pairs)
    np.savez_compressed(out/'neutral_af.npz',trajectories=trajectories,ids=ids)
    base.save(out/'neutral_af_manifest.json',dict(cache_sha256=base.digest(out/'neutral_af.npz'),
        contract=dict(ids=ids.tolist(),position_bp=500000,tracked_allele='Original focal allele introduced on Neanderthal donor haplotypes at 50 kya; no nearest-variant relabelling.',
            conditioning='Selected and neutral: '+CONDITION,neutral_pool='100 independent segregation-conditioned neutrals, shared across all selected cases.'),
        inputs=cohort['accepted'],raw_manifest_sha256=base.digest(args.raw/'manifest.json')))
    entries=[]
    for entry in old['cases']:
        assert base.digest(OLD/entry['cache'])==entry['cache_sha256']
        with np.load(OLD/entry['cache']) as data:content={k:data[k].copy() for k in data.files}
        count_c,count_i=exceedance_counts(content['cdfs'],cdfs)
        content.update(p_cumulative=(count_c+1)/101,p_interval=(count_i+1)/101,
            neutral_ge_cumulative=count_c,neutral_ge_interval=count_i)
        np.savez_compressed(out/entry['cache'],**content)
        previous=json.loads((OLD/entry['manifest']).read_text())
        previous['contract'].update(neutral_cache_sha256=base.digest(out/'neutral.npz'),conditioning='Selected and neutral: '+CONDITION)
        previous.update(cache_sha256=base.digest(out/entry['cache']),p_value_version=2,
            previous_selected_cache=dict(path=str(OLD/entry['cache']),sha256=entry['cache_sha256']),
            recalibration='Selected CDFs, AF trajectories, IDs and proposal indices unchanged; only neutral comparison and p values replaced.')
        base.save(out/entry['manifest'],previous)
        entries.append(dict(entry,cache_sha256=base.digest(out/entry['cache'])))
    meta=dict(old,cases=entries,neutral_cache_sha256=base.digest(out/'neutral.npz'),neutral_conditioned=True,
        neutral_source_parameters_sha256=base.digest(args.raw/'cohort.json'),neutral_pool='100 independent neutrals conditioned on original focal allele segregation today; shared across all cases.',
        neutral_condition=CONDITION,neutral_simulation_configuration=json.loads((args.raw/'manifest.json').read_text()),
        rows=ROW_ORDER,parent_report_manifest=dict(path=str(OLD/'manifest.json'),sha256=base.digest(OLD/'manifest.json')))
    base.save(out/'manifest.json',meta)
    shutil.copy2(args.raw/'cohort.json',out/'neutral_cohort_receipt.json')
    (out/'.gitattributes').write_text('* -text\n')
    base.save(out/'output_hashes.json',{p.name:base.digest(p) for p in sorted(out.iterdir()) if p.is_file() and p.name!='output_hashes.json'})
    print(f'Prepared conditioned analysis: {out}',flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('phase',choices=['run','prepare','all'])
    p.add_argument('--raw',type=Path,default=ROOT.parent/'sim/segregating_neutral_100_20261110')
    p.add_argument('--out',type=Path,default=ROOT/'docs/results/selection_calls_100_segregating_null_20261110')
    p.add_argument('--workers',type=int,default=4)
    p.add_argument('--target',type=int,default=100)
    p.add_argument('--seconds',type=int,default=14400)
    p.add_argument('--deadline-epoch',type=float)
    p.add_argument('--attempts',type=int,default=5000)
    p.add_argument('--independent-test',action='store_true',help='No reused simulations; separate seeds and IDs for a held-out neutral FPR cohort.')
    args=p.parse_args()
    assert 1<=args.workers<=24 and args.target>=8 and args.seconds>0 and args.attempts>0
    assert not args.independent_test or args.phase=='run','Use the independent FPR report script to prepare the test bank.'
    if args.phase in ('run','all'):run(args)
    if args.phase in ('prepare','all'):prepare(args)


if __name__=='__main__':main()
