#!/usr/bin/env python3
"""Bounded EAS ongoing-selection pilot; unconditional de novo and archaic pulse."""
import os
os.environ['MPLBACKEND']='Agg'
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
from pathlib import Path
import signal
import subprocess
import sys
import time
import msprime
import numpy as np
import pyslim
import tskit
import run_time_resolved_pilot as base


def onset(cfg, origin, label, directory):
    times, ne=base.history()
    lookup=lambda t: float(ne[np.searchsorted(times,t,side='right')-1])
    age=2000
    size=round(lookup(age))
    donors=round(2*size*.01) if origin=='introgression' else 0
    dem=msprime.Demography()
    dem.add_population(name='EAS',initial_size=size)
    if donors:
        dem.add_population(name='Neanderthal',initial_size=3600)
        dem.add_mass_migration(time=18000,source='Neanderthal',dest='EAS',proportion=1)
    for t,n in zip(times[times>age],ne[times>age]):
        dem.add_population_parameters_change(time=float(t)-age,initial_size=float(n),population='EAS')
    dem.sort_events()
    samples=[msprime.SampleSet(2*size-donors,population='EAS',ploidy=1)]
    if donors: samples.append(msprime.SampleSet(donors,population='Neanderthal',ploidy=1))
    initial=msprime.sim_ancestry(samples=samples,demography=dem,sequence_length=cfg['length'],
        recombination_rate=1e-8,random_seed=base.seed(label+':ancestry'),
        model=[msprime.DiscreteTimeWrightFisher(duration=100),msprime.StandardCoalescent()])
    nodes=initial.samples()
    if donors:
        tree=initial.at(cfg['length']//2)
        assert all(tree.tmrca(int(a),int(b))>=18000 for a in nodes[:10] for b in nodes[-donors:][:10]), 'Cross-population coalescence predates split'
    carrier_nodes=nodes[-donors:] if donors else nodes[:1]
    tables=initial.dump_tables()
    pop=tables.nodes.population.copy();pop[nodes]=0;tables.nodes.population=pop
    tables.individuals.clear()
    individual=np.full(initial.num_nodes,tskit.NULL,dtype=np.int32)
    pairs=np.random.default_rng(base.seed(label+':pairing')).permutation(nodes).reshape(-1,2)
    for pair in pairs: individual[pair]=tables.individuals.add_row()
    tables.nodes.individual=individual
    order=np.r_[pairs.ravel(),np.setdiff1d(np.arange(initial.num_nodes),nodes)]
    inverse=np.empty(initial.num_nodes,dtype=np.int32);inverse[order]=np.arange(initial.num_nodes)
    rows=list(tables.nodes);tables.nodes.clear()
    for i in order: tables.nodes.append(rows[int(i)])
    tables.edges.parent=inverse[tables.edges.parent];tables.edges.child=inverse[tables.edges.child]
    tables.sort()
    annotated=pyslim.annotate(tables.tree_sequence(),model_type='WF',tick=1,stage='late')
    annotated.dump(directory/'onset.trees')
    # These are the exact SLiM haplosome identities of donor chromosomes.
    carrier_ids=[annotated.node(int(inverse[n])).metadata['slim_id'] for n in carrier_nodes]
    base.save(directory/'onset_provenance.json',dict(population_size=size,donor_haplotypes=donors,
        selected_initial_copies=len(carrier_ids),carrier_slim_ids=carrier_ids,
        pulse_fraction_realized=donors/(2*size),split_years=500000 if donors else None,
        onset_sha256=base.digest(directory/'onset.trees')))
    return carrier_ids,[max(2,round(lookup(age-i))) for i in range(age+1)]


def simulate(cfg,task,directory):
    start=time.monotonic();origin=task['origin'];label='ongoing-v1:'+task['unit']
    carrier_ids,sizes=onset(cfg,origin,label,directory)
    s=task['selection'];length=cfg['length']
    schedules=[f'{i+1} early() {{ p0.setSubpopulationSize({n}); }}' for i,n in enumerate(sizes) if i and n!=sizes[i-1]]
    quote=lambda name:json.dumps(str(directory/name))
    reject_early='if (af==0.0 | af==1.0) { sim.simulationFinished(); }' if cfg.get('condition_segregating') else ''
    script=f'''initialize() {{
 initializeSLiMOptions(keepPedigrees=T); initializeTreeSeq(); initializeMutationRate(0);
 initializeMutationType("m0",0.5,"f",0.0);
 initializeMutationType("m2",0.5,"f",{s}); m2.convertToSubstitution=F;
 initializeGenomicElementType("g1",m0,1.0); initializeGenomicElement(g1,0,{length-1});
 initializeRecombinationRate(1e-8);
}}
1 late() {{
 sim.readFromPopulationFile({quote('onset.trees')});
 carriers=p0.haplosomes[match(p0.haplosomes.haplosomePedigreeID,c({','.join(map(str,carrier_ids))}))>=0];
 if (size(carriers)!={len(carrier_ids)}) stop("carrier identity mismatch");
 mut=carriers[0].addNewDrawnMutation(m2,{length//2});
 if (size(carriers)>1) carriers[1:(size(carriers)-1)].addMutations(mut);
 sim.recalculateFitness();
}}
{chr(10).join(schedules)}
1:2001 late() {{
 muts=sim.mutationsOfType(m2); af=0.0;
 if (size(muts)) af=sim.mutationFrequencies(p0,muts);
 catn("TRAJECTORY "+community.tick+" "+af+" "+p0.individualCount);
 {reject_early}
}}
2001 late() {{ sim.treeSeqOutput({quote('forward.trees')}); sim.simulationFinished(); }}
'''
    (directory/'model.slim').write_text(script)
    r=subprocess.run([cfg['slim'],'-s',str(base.seed(label+':slim')),str(directory/'model.slim')],capture_output=True,text=True)
    (directory/'slim.stdout.log').write_text(r.stdout);(directory/'slim.stderr.log').write_text(r.stderr)
    r.check_returncode()
    trajectory=np.array([list(map(float,line.split()[1:])) for line in r.stdout.splitlines() if line.startswith('TRAJECTORY ')])
    if not cfg.get('condition_segregating'):assert len(trajectory)==2001
    np.testing.assert_array_equal(trajectory[:,0],np.arange(1,len(trajectory)+1))
    np.testing.assert_array_equal(trajectory[:,2],sizes[:len(trajectory)])
    np.testing.assert_allclose(trajectory[0,1],len(carrier_ids)/(2*sizes[0]),rtol=1e-5)
    np.savetxt(directory/'trajectory.tsv',np.c_[(2001-trajectory[:,0])*25,trajectory[:,1:]],delimiter='\t',header='years_ago\taf\tdiploids',comments='')
    if cfg.get('condition_segregating') and not 0<trajectory[-1,1]<1:
        return dict(accepted=False,initial_frequency=len(carrier_ids)/(2*sizes[0]),
            final_frequency=float(trajectory[-1,1]),rejection='loss' if trajectory[-1,1]==0 else 'fixation',
            stop_years_ago=float((2001-trajectory[-1,0])*25),simulation_seconds=time.monotonic()-start,
            onset_sha256=base.digest(directory/'onset.trees'))
    assert len(trajectory)==2001
    sampled=base._sample_diploids(tskit.load(directory/'forward.trees'),cfg['diploids'],base.seed(label+':sample'))
    tables=sampled.dump_tables();tables.sites.clear();tables.mutations.clear();tables.time_units='generations'
    observed=msprime.sim_mutations(tables.tree_sequence(),rate=1.29e-8,random_seed=base.seed(label+':mutations'),model=msprime.JC69())
    observed.dump(directory/'observed.trees')
    return dict(initial_frequency=len(carrier_ids)/(2*sizes[0]),final_frequency=float(trajectory[-1,1]),
        maximum_frequency=float(trajectory[:,1].max()),survived=bool(trajectory[-1,1]>0),
        loss_years_ago=float((2001-trajectory[np.flatnonzero(trajectory[:,1]==0)[0],0])*25) if np.any(trajectory[:,1]==0) else None,
        fixation_years_ago=float((2001-trajectory[np.flatnonzero(trajectory[:,1]==1)[0],0])*25) if np.any(trajectory[:,1]==1) else None,
        reached_10_percent=bool(trajectory[:,1].max()>=.1),simulation_seconds=time.monotonic()-start,
        observed_sha256=base.digest(directory/'observed.trees'),onset_sha256=base.digest(directory/'onset.trees'),
        sites=observed.num_sites,trees=observed.num_trees)


def verify(directory,record):
    for key,name in [('scores_sha256','scores.npz'),('observed_sha256','observed.trees'),('posterior_sha256','posterior.zst')]:
        assert base.digest(directory/name)==record[key],f'Corrupt {directory/name}'


def worker(args):
    root=Path(args.out).resolve();manifest=json.loads((root/'manifest.json').read_text())
    task=next(t for t in json.loads((root/'tasks.json').read_text()) if t['unit']==args.unit)
    directory=root/args.unit;directory.mkdir(parents=True,exist_ok=True)
    record=dict(task,seed=base.seed('ongoing-v1:'+args.unit),config_sha256=base.digest(root/'manifest.json'))
    if (directory/'complete.json').exists():
        old=json.loads((directory/'complete.json').read_text());assert old['config_sha256']==record['config_sha256'];verify(directory,old);return
    if (directory/'simulation.json').exists():
        record=json.loads((directory/'simulation.json').read_text())
        assert record['config_sha256']==base.digest(root/'manifest.json')
        assert record['observed_sha256']==base.digest(directory/'observed.trees')
    else:
        record.update(simulate(manifest['config'],task,directory));base.save(directory/'simulation.json',record)
    record.update(base.score(manifest['config'],directory,root));base.save(directory/'complete.json',record)


def tasks(args):
    result=[]
    for i in range(max(args.development,args.calibration,args.test,args.selected)):
        for origin in ('de_novo','introgression'):
            for cohort,n in [('development',args.development),('calibration',args.calibration),('neutral_test',args.test)]:
                if i<n: result.append(dict(unit=f'{origin}/{cohort}_{i:04d}',origin=origin,cohort=cohort,selection=0.,replicate=i))
            if i<args.selected:
                for s in args.selection:
                    result.append(dict(unit=f'{origin}/s{s:.2f}/selected_{i:04d}',origin=origin,cohort='selected',selection=s,replicate=i))
    return result


def launch(args):
    root=Path(args.out).resolve();root.mkdir(parents=True,exist_ok=True)
    cfg=dict(schema='ongoing-eas-pilot/v1',length=1000000,diploids=100,pairs=1000,
        slim=str(Path(args.slim).resolve()),decoder=str(Path(args.decoder).resolve()),
        source_sha256=base.digest(__file__),score_source_sha256=base.digest(base.__file__),
        slim_sha256=base.digest(args.slim),decoder_sha256=base.digest(args.decoder),
        base_seed=base.BASE_SEED,pair_seed=1729,generation_time=25,scaling=1,dominance=.5,
        selection_start_years=50000,selection_end_years=0,pulse_fraction=.01,archaic_split_years=500000,
        archaic_ne=3600,mutation_rate=1.29e-8,recombination_rate=1e-8,stride=10000,
        cutoffs=list(range(1000,50001,1000)),prehistory='DTWF first 100 generations, then Hudson',
        donor_selected_allele='fixed in Neanderthal donors; absent in EAS',de_novo_initial_copies=1,
        phlash_sha256=base.digest(base.REPO/'no_introgression_EAS_sim/resources/EAS.npz'))
    manifest=root/'manifest.json'
    if manifest.exists(): assert json.loads(manifest.read_text())['config']==cfg,'Scientific contract changed'
    else: base.save(manifest,dict(config=cfg,python=sys.version,versions={n:__import__(n).__version__ for n in ('numpy','scipy','msprime','pyslim','tskit')},created_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())))
    possible=np.column_stack(np.triu_indices(200,1))
    pairs=possible[np.sort(np.random.default_rng(1729).choice(len(possible),1000,replace=False))]
    np.savetxt(root/'pairs.tsv',pairs,fmt='%d',delimiter='\t');np.savetxt(root/'positions.txt',np.arange(0,1000000,10000),fmt='%d')
    schedule=tasks(args);base.save(root/'tasks.json',schedule)
    deadline=time.time()+args.seconds
    base.save(root/'schedule.json',dict(workers=args.workers,deadline_epoch=deadline,timeout_per_unit=args.unit_timeout,targets=vars(args)))
    def execute(task):
        unit=task['unit'];directory=root/unit
        remaining=deadline-time.time()
        if remaining<5:return dict(unit=unit,state='not_started_deadline')
        directory.mkdir(parents=True,exist_ok=True)
        with (directory/'worker.stdout.log').open('w') as out,(directory/'worker.stderr.log').open('w') as err:
            proc=subprocess.Popen([sys.executable,__file__,'worker','--out',str(root),'--unit',unit],stdout=out,stderr=err,start_new_session=True)
            try:
                code=proc.wait(timeout=min(args.unit_timeout,remaining));state='complete' if code==0 else 'failed'
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid,signal.SIGKILL);proc.wait();state='timeout'
        result=dict(unit=unit,state=state);print(json.dumps(result),flush=True);return result
    records=[]
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for f in as_completed([pool.submit(execute,t) for t in schedule]):
            records.append(f.result());base.save(root/'status.json',records)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('phase',choices=['run','worker']);p.add_argument('--out',required=True);p.add_argument('--unit')
    p.add_argument('--slim');p.add_argument('--decoder');p.add_argument('--workers',type=int,default=4)
    p.add_argument('--seconds',type=int,default=1200);p.add_argument('--unit-timeout',type=int,default=240)
    p.add_argument('--selection',type=float,nargs='+',default=[.01,.02,.03,.04,.05])
    p.add_argument('--development',type=int,default=20);p.add_argument('--calibration',type=int,default=39)
    p.add_argument('--test',type=int,default=20);p.add_argument('--selected',type=int,default=20)
    a=p.parse_args();worker(a) if a.phase=='worker' else launch(a)
