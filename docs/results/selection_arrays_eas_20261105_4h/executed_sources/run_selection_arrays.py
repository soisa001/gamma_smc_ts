#!/usr/bin/env python3
"""Bounded recent-onset EAS sweeps with independent nulls and exact rejection."""
import os
os.environ['MPLBACKEND']='Agg'
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse
import ast
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
import json
from pathlib import Path
import signal
import shutil
import subprocess
import sys
import time
import numpy as np
import msprime
import tskit
import run_ongoing_selection_pilot as ancestor
base=ancestor.base


def timing(origin,onset_years=10000):
    assert onset_years%25==0 and 0<onset_years<=50000
    selection_generations=onset_years//25
    age=2000 if origin=='introgression' else selection_generations
    first_selected_tick=age-selection_generations+2
    return age,first_selected_tick


def simulate(cfg,task,directory):
    started=time.monotonic();origin=task['origin'];label=cfg['seed_namespace']+task['unit']
    age,first=timing(origin,cfg['selection_start_years']);last=age+1
    end=cfg['selection_end_years'];assert end%25==0 and 0<=end<cfg['selection_start_years']
    stop=age-end//25+2
    carrier_ids,sizes=ancestor.onset(dict(cfg,ancestry_start_generations=age),origin,label,directory)
    s=task['selection'];selected=task['cohort']=='selected';initial_s=s if first==2 else 0.
    schedules=[f'{i+1} early() {{ p0.setSubpopulationSize({n}); }}' for i,n in enumerate(sizes) if i and n!=sizes[i-1]]
    quote=lambda name:json.dumps(str(directory/name))
    reject='if (af==0.0 | af==1.0) sim.simulationFinished();' if selected else ''
    activate=f'''{first} early() {{
 muts=sim.mutationsOfType(m2); if (size(muts)) muts.setSelectionCoeff({s});
 sim.recalculateFitness();
}}''' if first>2 else ''
    deactivate=f'{stop} early() {{ muts=sim.mutationsOfType(m2); if (size(muts)) muts.setSelectionCoeff(0.0); sim.recalculateFitness(); }}' if end else ''
    script=f'''initialize() {{
 initializeSLiMOptions(keepPedigrees=T); initializeTreeSeq(); initializeMutationRate(0);
 initializeMutationType("m0",0.5,"f",0.0);
 initializeMutationType("m2",0.5,"f",{initial_s}); m2.convertToSubstitution=F;
 initializeGenomicElementType("g1",m0,1.0); initializeGenomicElement(g1,0,{cfg['length']-1});
 initializeRecombinationRate(1e-8);
}}
1 late() {{
 sim.readFromPopulationFile({quote('onset.trees')});
 carriers=p0.haplosomes[match(p0.haplosomes.haplosomePedigreeID,c({','.join(map(str,carrier_ids))}))>=0];
 if (size(carriers)!={len(carrier_ids)}) stop("carrier identity mismatch");
 mut=carriers[0].addNewDrawnMutation(m2,{cfg['length']//2});
 if (size(carriers)>1) carriers[1:(size(carriers)-1)].addMutations(mut);
 sim.recalculateFitness();
}}
{chr(10).join(schedules)}
{activate}
{deactivate}
1:{last} late() {{
 muts=sim.mutationsOfType(m2); af=0.0; coefficient=0.0;
 if (size(muts)) {{ af=sim.mutationFrequencies(p0,muts); coefficient=muts.selectionCoeff; }}
 copies=0; if (size(muts)) copies=sim.mutationCounts(p0,muts);
 catn("TRAJECTORY "+community.tick+" "+af+" "+p0.individualCount+" "+coefficient+" "+copies);
 {reject}
}}
{last} late() {{ sim.treeSeqOutput({quote('forward.trees')}); sim.simulationFinished(); }}
'''
    (directory/'model.slim').write_text(script)
    r=subprocess.run([cfg['slim'],'-s',str(base.seed(label+':slim')),str(directory/'model.slim')],capture_output=True,text=True)
    (directory/'slim.stdout.log').write_text(r.stdout);(directory/'slim.stderr.log').write_text(r.stderr);r.check_returncode()
    trajectory=np.array([list(map(float,line.split()[1:])) for line in r.stdout.splitlines() if line.startswith('TRAJECTORY ')])
    np.testing.assert_array_equal(trajectory[:,0],np.arange(1,len(trajectory)+1))
    np.testing.assert_array_equal(trajectory[:,2],sizes[:len(trajectory)])
    exact=trajectory[:,4]/(2*trajectory[:,2])
    np.testing.assert_allclose(exact,trajectory[:,1],atol=5.1e-7)
    np.testing.assert_equal(int(trajectory[0,4]),len(carrier_ids))
    # Initial mutation already has s for a de novo onset; introgressed alleles
    # must remain neutral through the last pre-selection parental generation.
    expected=np.where(trajectory[:,0]>=stop,0.,np.where(trajectory[:,0]>=first,s,initial_s))
    if selected:np.testing.assert_allclose(trajectory[:,3][trajectory[:,4]>0],expected[trajectory[:,4]>0],atol=1e-10)
    else:assert np.all(trajectory[:,3]==0)
    np.savetxt(directory/'trajectory.tsv',np.c_[(last-trajectory[:,0])*25,exact,trajectory[:,2:]],delimiter='\t',
        header='years_ago\taf\tdiploids\tselection_coefficient\tfocal_copies',comments='')
    accepted=bool(not selected or 0<exact[-1]<1)
    at_onset=age-cfg['selection_start_years']//25
    record=dict(accepted=accepted,initial_frequency=float(exact[0]),final_frequency=float(exact[-1]),
        selection_onset_frequency=float(exact[at_onset]) if at_onset<len(exact) else None,
        rejection=None if accepted else ('loss' if exact[-1]==0 else 'fixation'),
        stop_years_ago=float((last-trajectory[-1,0])*25),selected_reproduction_ticks=[first,min(last,stop-1)],
        selection_start_years=cfg['selection_start_years'],selection_end_years=end,seed_namespace=cfg['seed_namespace'],
        ancestry_start_years=age*25,trajectory_rows=len(trajectory),simulation_seconds=time.monotonic()-started,
        onset_sha256=base.digest(directory/'onset.trees'),trajectory_sha256=base.digest(directory/'trajectory.tsv'))
    if not accepted:return record
    assert len(trajectory)==last
    sampled=base._sample_diploids(tskit.load(directory/'forward.trees'),cfg['diploids'],base.seed(label+':sample'))
    tables=sampled.dump_tables();tables.sites.clear();tables.mutations.clear();tables.time_units='generations'
    observed=msprime.sim_mutations(tables.tree_sequence(),rate=1.29e-8,random_seed=base.seed(label+':mutations'),model=msprime.JC69())
    observed.dump(directory/'observed.trees')
    record.update(observed_sha256=base.digest(directory/'observed.trees'),sites=observed.num_sites,trees=observed.num_trees)
    if not selected:
        variants=[]
        for variant in observed.variants():
            alleles,counts=np.unique(variant.genotypes,return_counts=True)
            if len(alleles)>1:variants.append((abs(variant.site.position-cfg['length']/2),variant.site.position,int(counts.min())))
        if not variants:raise ValueError('No present-day segregating neutral focal variant')
        distance,position,minor=min(variants)
        record.update(neutral_focal_position_bp=position,neutral_focal_minor_count=minor,
            neutral_focal_offset_bp=position-cfg['length']/2)
    return record


def work(payload):
    root,cfg,task,deadline=payload;root=Path(root)
    if time.time()>=deadline:return dict(id=task['id'],state='not_started_deadline')
    directory=root/task['id'];directory.mkdir(parents=True,exist_ok=True)
    signature=base.digest(root/'manifest.json')
    for filename in ('complete.json','rejected.json'):
        if (directory/filename).exists():
            record=json.loads((directory/filename).read_text());assert record['config_sha256']==signature
            for name,key in [('trajectory.tsv','trajectory_sha256'),('onset.trees','onset_sha256')]:
                assert base.digest(directory/name)==record[key]
            if filename=='complete.json':ancestor.verify(directory,record)
            return dict(id=task['id'],state='accepted' if record['accepted'] else 'rejected')
    receipt=directory/'simulation.json'
    if receipt.exists():
        record=json.loads(receipt.read_text());assert record['config_sha256']==signature
        assert base.digest(directory/'observed.trees')==record['observed_sha256']
    else:
        record=dict(task,seed=base.seed(cfg['seed_namespace']+task['unit']),config_sha256=signature)
        record.update(simulate(cfg,task,directory))
        if not record['accepted']:
            base.save(directory/'rejected.json',record);return dict(id=task['id'],state='rejected')
        base.save(receipt,record)
    record.update(base.score(cfg,directory,root));base.save(directory/'complete.json',record)
    return dict(id=task['id'],state='accepted')


def launch(a):
    root=Path(a.out).resolve();root.mkdir(parents=True,exist_ok=True)
    common=dict(length=1000000,diploids=100,pairs=1000,stride=10000,generation_time=25,dominance=.5,
        pulse_years=50000,pulse_fraction=.01,archaic_split_years=500000,archaic_ne=3600,scaling=1,
        mutation_rate=1.29e-8,recombination_rate=1e-8,base_seed=base.BASE_SEED,pair_seed=1729,
        slim=str(Path(a.slim).resolve()),decoder=str(Path(a.decoder).resolve()),
        selected_condition='0 < whole-population focal AF < 1 at present',
        null_condition='unconditioned; focal label nearest sample-polymorphic site',
        source_sha256=base.digest(__file__),ancestry_source_sha256=base.digest(ancestor.__file__),
        score_source_sha256=base.digest(base.__file__),slim_sha256=base.digest(a.slim),decoder_sha256=base.digest(a.decoder),
        phlash_sha256=base.digest(base.REPO/'no_introgression_EAS_sim/resources/EAS.npz'))
    configs={name:dict(common,selection_start_years=start,selection_end_years=end,
        seed_namespace='recent-v1:' if name=='10k' else f'episode-v1:{name}:')
        for name,start,end in [('10k',10000,0),('50k',50000,0),('40_20k',40000,20000)]}
    contract=dict(schema='selection-arrays/v1',configs=configs,selection=a.selection,target=a.accepted,
        development=a.development,calibration=a.calibration,neutral_test=a.test,reuse_recent=str(Path(a.reuse_recent).resolve()))
    if (root/'manifest.json').exists():assert json.loads((root/'manifest.json').read_text())['contract']==contract
    else:base.save(root/'manifest.json',dict(contract=contract,versions={n:__import__(n).__version__ for n in ('numpy','scipy','msprime','pyslim','tskit')},python=sys.version))
    snapshots=root/'executed_sources';snapshots.mkdir(exist_ok=True)
    for filename in (__file__,ancestor.__file__,base.__file__):
        dest=snapshots/Path(filename).name
        if dest.exists():assert base.digest(dest)==base.digest(filename)
        else:shutil.copy2(filename,dest)
    possible=np.column_stack(np.triu_indices(200,1))
    pairs=possible[np.sort(np.random.default_rng(1729).choice(len(possible),1000,replace=False))]
    np.savetxt(root/'pairs.tsv',pairs,fmt='%d',delimiter='\t');np.savetxt(root/'positions.txt',np.arange(0,1000000,10000),fmt='%d')
    old=Path(a.reuse_recent);oldcfg=json.loads((old/'manifest.json').read_text())['config']
    # Reuse is explicitly limited to unchanged 10-kya physics and neutral introgression.
    for key in ('length','diploids','pairs','stride','generation_time','dominance','pulse_years','pulse_fraction','archaic_split_years','archaic_ne','mutation_rate','recombination_rate','base_seed','pair_seed','slim_sha256','decoder_sha256','phlash_sha256','ancestry_source_sha256','score_source_sha256'):
        assert common[key]==oldcfg[key],f'Reuse contract differs: {key}'
    assert oldcfg['selection_start_years']==10000 and oldcfg['selection_end_years']==0
    np.testing.assert_array_equal(pairs,np.loadtxt(old/'pairs.tsv',dtype=int))
    strata=[]
    for regime in configs:
        for origin in ('de_novo','introgression'):
            if regime=='10k' or origin=='de_novo':
                for cohort,n in [('development',a.development),('calibration',a.calibration),('neutral_test',a.test)]:
                    strata.append(dict(regime=regime,origin=origin,cohort=cohort,selection=0.,target=n,limit=n))
            for s in a.selection:strata.append(dict(regime=regime,origin=origin,cohort='selected',selection=s,target=a.accepted,limit=a.attempts))
    started=time.time();deadline=min(started+a.seconds,a.deadline_epoch or float('inf'))
    base.save(root/'schedule.json',dict(parameters=vars(a),strata=strata,started_epoch=started,deadline_epoch=deadline,
        null_pool_mapping={r:{o:('10k' if o=='introgression' else r) for o in ('de_novo','introgression')} for r in configs}))
    tasks=[];results=[];next_index=[0]*len(strata);done=[{} for _ in strata];pending={};cursor=0
    def make_task(k,i):
        s=strata[k];suffix=f"s{s['selection']:g}/selected" if s['cohort']=='selected' else s['cohort']
        unit=f"{s['origin']}/{suffix}_{i:05d}"
        return dict(id=f"{s['regime']}/{unit}",unit=unit,regime=s['regime'],origin=s['origin'],cohort=s['cohort'],selection=s['selection'],replicate=i)
    # Register contiguous historical proposals without copying or changing receipts.
    for k,s in enumerate(strata):
        if s['regime']!='10k':continue
        accepted=0
        for i in range(s['limit']):
            if accepted>=s['target']:break
            task=make_task(k,i);directory=old/task['unit']
            receipt=next((directory/n for n in ('complete.json','rejected.json') if (directory/n).exists()),None)
            if receipt is None:break
            rec=json.loads(receipt.read_text());assert rec['config_sha256']==base.digest(old/'manifest.json')
            assert rec['selection']==task['selection'] and rec['origin']==task['origin'] and rec['cohort']==task['cohort']
            task.update(source_directory=str(directory),source_receipt_sha256=base.digest(receipt))
            state='accepted' if rec['accepted'] else 'rejected';accepted+=rec['accepted']
            tasks.append(task);results.append(dict(id=task['id'],state=state,reused=True));done[k][i]=state;next_index[k]=i+1
    def prefix_accepts(k):
        n=0
        for i in range(next_index[k]):
            if i not in done[k]:break
            n+=done[k][i]=='accepted'
        return n
    def next_task():
        nonlocal cursor
        for _ in strata:
            k=cursor%len(strata);cursor+=1;s=strata[k]
            if next_index[k]>=s['limit'] or prefix_accepts(k)>=s['target']:continue
            i=next_index[k];next_index[k]+=1;task=make_task(k,i);tasks.append(task);return k,task
        return None
    pool=ProcessPoolExecutor(max_workers=a.workers,initializer=os.setsid)
    try:
        while time.time()<deadline:
            while len(pending)<a.workers and time.time()<deadline-a.drain_seconds:
                item=next_task()
                if item is None:break
                k,task=item;f=pool.submit(work,(str(root),configs[task['regime']],task,deadline));pending[f]=(k,task)
            base.save(root/'tasks.json',tasks)
            if not pending:break
            finished,_=wait(pending,timeout=min(1,max(0,deadline-time.time())),return_when=FIRST_COMPLETED)
            for f in finished:
                k,task=pending.pop(f)
                try:result=f.result()
                except Exception as exc:result=dict(id=task['id'],state='failed',error=repr(exc))
                done[k][task['replicate']]=result['state'];results.append(result)
                base.save(root/'status.json',results)
        if pending:
            for process in pool._processes.values():
                try:os.killpg(process.pid,signal.SIGKILL)
                except ProcessLookupError:pass
    finally:pool.shutdown(wait=True,cancel_futures=True)
    base.save(root/'tasks.json',tasks);base.save(root/'status.json',results)
    base.save(root/'finished.json',dict(finished_epoch=time.time(),submitted=len(tasks),returned=len(results),
        prefix_accepts=[prefix_accepts(k) for k in range(len(strata))],strata=strata,pending_at_deadline=len(pending)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('out','slim','decoder','reuse-recent'):p.add_argument('--'+key,required=True)
    p.add_argument('--selection',type=float,nargs='+',default=[.001,.005,.01,.02,.03,.04,.05])
    p.add_argument('--workers',type=int,default=4);p.add_argument('--seconds',type=int,default=3300)
    p.add_argument('--deadline-epoch',type=float);p.add_argument('--drain-seconds',type=int,default=120)
    p.add_argument('--development',type=int,default=20);p.add_argument('--calibration',type=int,default=99);p.add_argument('--test',type=int,default=50)
    p.add_argument('--accepted',type=int,default=100);p.add_argument('--attempts',type=int,default=10000)
    launch(p.parse_args())

