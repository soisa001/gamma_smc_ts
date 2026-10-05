"""Mapped genomic-region nulls, grouped by correlated PHLASH demographic draw."""
from __future__ import annotations
import os
os.environ['MPLBACKEND']='Agg'
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
import argparse
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
import importlib.metadata
import json
from pathlib import Path
import subprocess
import time
import traceback
import numpy as np
import pandas as pd
from . import population_neutral as pn
from . import origin_onset as oo
from . import segregating_introgression as si

VERSION='mapped-regions-phlash-mvn/v1'


def initialise(root, resources, histories, n_histories, replicates, workers, slim):
    assert 1<=workers<=20 and replicates % n_histories == 0
    study=json.loads((oo.REPO/'configs/all_population_neutral.json').read_text())
    study['parameters'].update(null_replicates=replicates,workers=workers,mutation_rate=1.29e-8,
        decoder_recombination_model='cropped_region_mean', decoder_recombination_to_mutation_ratio=1e-8/1.29e-8)
    study.update(reuse_eas_root=str(root/'no_legacy_reuse'),reuse_eas_report=str(root/'no_legacy_reuse'))
    configs=pn.configs(study)
    regions=pd.read_csv(resources/'regions.csv').iloc[:replicates]
    assert len(regions)==replicates
    region_receipt=json.loads((resources/'manifest.json').read_text())
    for filename,sha in region_receipt['files'].items():
        assert oo.legacy.digest(resources/filename)==sha
    result=subprocess.run([slim,'-v'],capture_output=True,text=True,check=True)
    assert 'SLiM version 5.' in result.stdout+result.stderr
    runtime=dict(versions={p:importlib.metadata.version(p) for p in ('msprime','numpy','pyslim','stdpopsim','tskit')},slim_sha256=oo.legacy.digest(slim))
    decoder_hash=oo.legacy.digest(oo.REPO/'bin/gamma_smc')
    identity=dict(schema=VERSION,parameters={k:v for k,v in study['parameters'].items() if k not in ('workers','memory_limit_gb','worker_memory_budget_gb')},
        n_histories=n_histories,regions=oo.legacy.digest(resources/'manifest.json'),history_sources={},runtime=runtime,decoder_sha256=decoder_hash)
    choices=np.column_stack(np.triu_indices(2*study['parameters']['sample_diploids'],k=1))
    chosen=np.random.default_rng(study['parameters']['pairs_seed']).choice(len(choices),study['parameters']['haplotype_pairs'],replace=False)
    pairs=choices[np.sort(chosen)]
    draw_specs={}
    for pop in pn.POPULATIONS:
        pop_root=root/pop
        pop_root.mkdir(parents=True,exist_ok=True)
        pairpath=pop_root/'pairs.tsv'
        if pairpath.exists():
            np.testing.assert_array_equal(np.loadtxt(pairpath,dtype=int),pairs)
        else:
            np.savetxt(pairpath,pairs,fmt='%d',delimiter='\t')
        source=histories/f'{pop}_draws.npz'
        receipt=json.loads(source.with_suffix('.json').read_text())
        assert oo.legacy.digest(source)==receipt['sha256']
        assert receipt['identity']['source_sha256']==configs[pop]['phlash_sha256']
        identity['history_sources'][pop]=receipt
        drawdir=pop_root/'demographic_draws'
        drawdir.mkdir(exist_ok=True)
        with np.load(source,allow_pickle=False) as d:
            time_grid, ne = d['time_generations'],d['ne']
        assert n_histories<=len(ne)
        draw_specs[pop]=[]
        for i in range(n_histories):
            target=drawdir/f'history{i:04d}.npz'
            if target.exists():
                with np.load(target,allow_pickle=False) as d:
                    np.testing.assert_array_equal(d['ne'],ne[i])
                    np.testing.assert_array_equal(d['time_generations'],time_grid)
            else:
                np.savez_compressed(target,ne=ne[i],time_generations=time_grid)
            draw_specs[pop].append(dict(path=str(target),sha256=oo.legacy.digest(target),history_id=i,source_sha256=receipt['sha256']))
    manifest=root/'manifest.json'
    if manifest.exists():
        assert json.loads(manifest.read_text())['identity']==identity,'Incompatible campaign manifest'
    else:
        oo.legacy.atomic_json(manifest,dict(identity=identity,study=study,created=time.time(),
            replicates_per_population=replicates,histories_per_population=n_histories,regions_per_history=replicates//n_histories,
            conditioning='nearest archaic-derived allele segregating at 50 kya; observed in present sample; fixation retained',
            uncertainty='Between-history empirical CDF bands also contain finite-region Monte Carlo variation; not confidence bands for a mean CDF',
            mask_policy=region_receipt['masking'],recombination_model='Full local map in ancestry and SLiM; regional mean in native decoder',
            resources=str(resources),histories=str(histories)))
    record=json.loads(manifest.read_text())
    if 'implementation_sha256' not in record:
        record['implementation_sha256']={name:oo.legacy.digest(Path(__file__).with_name(name)) for name in
            ('genomic_neutral.py','genomic_neutral_report.py','population_neutral.py','segregating_introgression.py',
             'origin_onset.py','stopped_campaign_report.py')}
        oo.legacy.atomic_json(manifest,record)
    plan=[]
    for pop,task in pn.jobs(configs):
        rep=task['replicate']
        region=regions.iloc[rep]
        cfg=dict(configs[pop], demographic_draw=draw_specs[pop][rep % n_histories],
            genomic_region=dict(path=str(resources/region.map_file),sha256=region.sha256,
                chromosome=region.chromosome,start_0based=int(region.start_0based),end_exclusive=int(region.end_exclusive),
                callable_fraction=float(region.callable_fraction),map_coverage=float(region.map_coverage)),
            recombination_rate=float(region.mean_recombination_rate))
        plan.append((pop,task,cfg))
    return study,configs,runtime,decoder_hash,plan


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--regions',type=Path,required=True)
    p.add_argument('--histories',type=Path,required=True)
    p.add_argument('--n-histories',type=int,default=100)
    p.add_argument('--replicates',type=int,default=1000)
    p.add_argument('--workers',type=int,default=20)
    p.add_argument('--slim',required=True)
    p.add_argument('--phase',choices=('plan','simulate','decode','run','report'),default='plan')
    p.add_argument('--limit',type=int,help='Execute only the first N planned tasks; leave the full manifest intact')
    args=p.parse_args()
    args.out=args.out.resolve()
    args.regions=args.regions.resolve()
    args.histories=args.histories.resolve()
    args.out.mkdir(parents=True,exist_ok=True)
    import fcntl
    with (args.out/'study.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        study,configs,runtime,decoder,plan=initialise(args.out,args.regions,args.histories,args.n_histories,args.replicates,args.workers,args.slim)
        if args.phase=='plan':
            print(json.dumps(dict(planned=len(plan),workers=args.workers,output=str(args.out))),flush=True)
            return
        if args.phase=='report':
            pn.summarize(args.out,configs,figures=True)
            from .genomic_neutral_report import report
            report(args.out)
            return
        planned=plan[:args.limit] if args.limit else plan
        remaining=iter(planned)
        active={}
        results=[]
        failure=False
        def status(state):
            oo.legacy.atomic_json(args.out/'run_status.json',dict(state=state,pid=os.getpid(),workers=args.workers,
                memory_limit_gb=200,worker_memory_budget_gb=180,total_planned=len(plan),scheduled=len(planned),
                completed=len(results),failed=[r for r in results if r['state']=='failed'],active=list(active.values()),updated=time.time()))
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            def submit():
                item=next(remaining,None)
                if item is None:
                    return False
                pop,task,cfg=item
                payload=(str(args.out),study,cfg,task,runtime,args.slim,decoder,args.phase)
                active[pool.submit(pn.run_job,payload)]=f'{pop}/{task["id"]}'
                return True
            for _ in range(args.workers):
                if not submit(): break
            status('running')
            previous_summary=0
            while active:
                done,_=wait(active,timeout=30,return_when=FIRST_COMPLETED)
                for future in done:
                    name=active.pop(future)
                    try:
                        result=future.result()
                    except Exception:
                        failure=True
                        result=dict(task_id=name,state='failed',error=traceback.format_exc())
                    results.append(result)
                    print(json.dumps(result),flush=True)
                if not failure:
                    while len(active)<args.workers and submit(): pass
                status('finishing_after_failure' if failure else 'running')
                if len(results)-previous_summary>=20:
                    pn.summarize(args.out,configs)
                    previous_summary=len(results)
        pn.summarize(args.out,configs,figures=not failure)
        if not failure:
            from .genomic_neutral_report import report
            report(args.out)
        status('failed' if failure else 'complete')
        if failure:
            raise RuntimeError('Task failed; see run_status.json')


if __name__=='__main__':
    main()
