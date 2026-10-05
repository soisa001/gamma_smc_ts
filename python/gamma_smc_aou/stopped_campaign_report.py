"""Read-only analysis of a frozen, incomplete campaign; never creates simulations."""
from __future__ import annotations
import os
os.environ['MPLBACKEND'] = 'Agg'
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
from pathlib import Path
import time
import tempfile
import numpy as np
import pandas as pd
import tskit
from . import origin_onset as oo
from .origin_onset_analysis import rank_p, wilson
from .pair_class_profiles import pair_classes, posterior_arrays, summarize_decoded, build_frame
from .decoder import run_within_decoder
from .tree_sequence import diploid_individuals

VERSION = 'stopped-campaign-report/v1'

def checked(directory, record, name):
    path = directory/name
    spec = record['artifacts'][name]
    if path.stat().st_size != spec['bytes'] or oo.legacy.digest(path) != spec['sha256']:
        raise ValueError(f'Corrupt report input: {path}')
    return path

def collect(root, out):
    manifest = json.loads((root/'manifest.json').read_text())
    cfg = manifest['config']
    rows, missing = [], []
    for task in manifest['tasks']:
        directory = root/'regions'/task['id']
        receipt = directory/'simulation.json'
        if not receipt.exists():
            missing.append(task)
            continue
        sim = json.loads(receipt.read_text())
        assert sim['identity'] == oo.simulation_identity(cfg,task,manifest['runtime'])
        assert sim['task'] == task
        assert sim['seed'] == oo.seed_for(cfg,task,sim['attempts'][-1]['attempt'])
        rows.append(dict(task_id=task['id'],family=task['family'],arm=task['family'],
            role=task['role'],s=task['s'],replicate=task['replicate'],seed=sim['seed'],
            sample_af=sim['sample_af'],final_census_af=sim['final_census_af'],
            onset_af=sim['focal_choice']['onset_af'],attempts=len(sim['attempts']),
            rejected_attempts=len(sim['attempts'])-1,
            attempt_seconds=sum(a['seconds'] for a in sim['attempts']),
            focal_original=sim['focal_position_original'],
            focal_distance_bp=abs(sim['focal_position_original']-cfg['simulated_length_bp']/2),
            mutation_age_years=sim['focal_choice']['mutation_time_generations']*cfg['generation_time_years'],
            simulation_sha256=oo.legacy.digest(receipt)))
    inventory = pd.DataFrame(rows)
    assert not inventory.task_id.duplicated().any() and not inventory.seed.duplicated().any()
    inventory.to_csv(out/'inventory.csv',index=False)
    pd.DataFrame(missing).to_csv(out/'unfinished_tasks.csv',index=False)
    plan = pd.DataFrame(manifest['tasks'])
    counts = plan.groupby(['family','role','s']).size().rename('planned').reset_index()
    observed = inventory.groupby(['family','role','s']).size().rename('completed').reset_index()
    counts = counts.merge(observed,on=['family','role','s'],how='left').fillna({'completed':0})
    counts['completed'] = counts.completed.astype(int)
    counts['missing'] = counts.planned-counts.completed
    counts.to_csv(out/'cohort_counts.csv',index=False)
    oo.legacy.atomic_json(out/'snapshot.json',dict(schema=VERSION,created=time.time(),
        root=str(root),config=cfg,completed=len(inventory),planned=len(plan),
        source_manifest_sha256=oo.legacy.digest(root/'manifest.json'),
        source_pairs_sha256=oo.legacy.digest(root/'pairs.tsv'),
        source_stop_status=json.loads((root/'run_status.json').read_text()),
        analysis_source_sha256=oo.legacy.digest(Path(__file__)),
        task_ids=inventory.task_id.tolist(),simulation_generation=False))
    return cfg, inventory

def truth_at_positions(ts, pairs, classes, positions, cfg):
    nodes = oo.legacy.ordered_nodes(ts)[pairs]
    thresholds = np.asarray(cfg['tmrca_cutoffs_years'])/cfg['generation_time_years']
    rows = []
    for pos in positions:
        tree = ts.at(float(pos))
        times = np.fromiter((tree.tmrca(int(a),int(b)) for a,b in nodes),float,count=len(nodes))
        assert np.isfinite(times).all() and (times>=0).all()
        for method, mask in [('all_pairs',np.ones(len(pairs),bool)),('alt_alt',classes==2)]:
            n = int(mask.sum())
            for cutoff, generations in zip(cfg['tmrca_cutoffs_years'],thresholds):
                rows.append(dict(source='truth',method=method,position_0based=int(pos),
                    cutoff_years=cutoff,n_pairs=n,
                    score=float(np.mean(times[mask]<generations)) if n else np.nan))
    return pd.DataFrame(rows)

def compute_one(payload):
    root_s,out_s,row,cfg,pairs_hash,mode,example = payload
    root,out = Path(root_s),Path(out_s)
    directory = root/'regions'/row['task_id']
    dest = out/'regions'/row['task_id']
    dest.mkdir(parents=True,exist_ok=True)
    receipt_path = directory/'simulation.json'
    assert oo.legacy.digest(receipt_path)==row['simulation_sha256']
    sim = json.loads(receipt_path.read_text())
    positions = (np.arange(cfg['focal_position_bp']-500000,cfg['focal_position_bp']+500001,10000)
                 if example else np.array([cfg['focal_position_bp']]))
    identity = dict(schema=VERSION,simulation=row['simulation_sha256'],pairs=pairs_hash,
        positions=positions.tolist(),mode=mode,cutoffs=cfg['tmrca_cutoffs_years'])
    if mode=='decoded':
        identity.update(decoder=oo.legacy.digest(oo.REPO/'bin/gamma_smc'),
            parameters={k:v for k,v in cfg.items() if k.startswith('decoder_') or k in
                ('generation_time_years','mutation_rate','recent_call')})
    cache_path = dest/(mode+'.json')
    score_path = dest/(mode+'.csv')
    if cache_path.exists():
        cache = json.loads(cache_path.read_text())
        if cache['identity']==identity and oo.legacy.digest(score_path)==cache['scores_sha256']:
            return dict(task_id=row['task_id'],mode=mode,state='cached')
        # A revised representative requests more/fewer positions from the same input.
        # Recompute only this derived profile; do not reuse incompatible values.
        old_identity={k:v for k,v in cache['identity'].items() if k!='positions'}
        new_identity={k:v for k,v in identity.items() if k!='positions'}
        if old_identity!=new_identity or cache['identity']['positions']==identity['positions']:
            raise ValueError(f'Incompatible or corrupt report cache {cache_path}')
    started = time.monotonic()
    tree_path = checked(directory,sim,'decoded_input.trees')
    carriers = np.load(checked(directory,sim,'focal_carriers.npy'),allow_pickle=False)
    pairs = np.loadtxt(root/'pairs.tsv',dtype=int,ndmin=2)
    ts = tskit.load(tree_path)
    focal = oo.legacy.exact_focal_variant(ts,oo.legacy.ordered_nodes(ts),cfg['focal_position_bp'])
    assert focal is not None
    np.testing.assert_array_equal(focal['carriers'],carriers)
    assert np.isclose(float(carriers.mean()),row['sample_af'],rtol=0,atol=1e-12)
    assert ts.sequence_length==cfg['scored_length_bp']
    classes = pair_classes(carriers,pairs)
    if mode=='truth':
        checked(directory,sim,'trajectory.csv')
        scores = truth_at_positions(ts,pairs,classes,positions,cfg)
    else:
        position_path = dest/'positions.txt'
        np.savetxt(position_path,positions,fmt='%d')
        raw = dest/'posteriors.zst'
        # The worker already loaded tskit. Avoid a fresh Python import subprocess
        # per region; an in-memory temporary VCF preserves the exact same input.
        temporary=tempfile.TemporaryDirectory(prefix='gamma-report-',dir='/dev/shm')
        vcf=Path(temporary.name)/'input.vcf'
        site_positions=ts.tables.sites.position
        assert np.array_equal(site_positions,np.floor(site_positions))
        with vcf.open('w') as stream:
            ts.write_vcf(stream,individuals=diploid_individuals(ts),
                position_transform=lambda x: np.asarray(x,dtype=np.int64)+1)
        vcf_sha256=oo.legacy.digest(vcf)
        decoder_mask = None
        recombination_ratio = cfg['decoder_recombination_to_mutation_ratio']
        if cfg.get('genomic_region'):
            from .segregating_introgression import rate_map
            offset = sim['crop_offset']
            end = offset + cfg['scored_length_bp']
            mutation_map = rate_map(cfg, 'mutation').slice(left=offset, right=end, trim=True)
            decoder_mask = dest/'simulation_callable.bed'
            with decoder_mask.open('w') as stream:
                for left, right, rate in zip(mutation_map.position[:-1], mutation_map.position[1:], mutation_map.rate):
                    if rate > 0:
                        stream.write(f'1\t{int(left)}\t{int(right)}\n')
            # The native decoder currently accepts a scalar recombination rate.
            # Record/use the cropped region's mean, while simulation retains hotspots.
            recombination_ratio = float(rate_map(cfg, 'recombination').slice(left=offset, right=end, trim=True).mean_rate / cfg['mutation_rate'])
        result = run_within_decoder(oo.REPO/'bin/gamma_smc',vcf,dest/'native_mean.tsv',input_format='vcf',
            raw_output=raw,threshold_years=cfg['tmrca_cutoffs_years'],generation_time=cfg['generation_time_years'],
            mutation_rate=cfg['mutation_rate'],scaled_mutation_rate=cfg['decoder_scaled_mutation_rate'],
            recombination_to_mutation_ratio=recombination_ratio, mask=decoder_mask,
            output_at_stride=-1,output_at_hets=False,only_within=False,
            output_positions_file=position_path,pairs_file=root/'pairs.tsv',
            vcf_position_transform='one_based',recent_call='mean',threads=1,
            cache_size=cfg['decoder_cache_size'],pair_block=cfg['decoder_pair_block'],
            exp10=cfg['decoder_exp10'],backward_alignment=cfg['decoder_backward_alignment'],
            extra_args=['--no_recent_probability'])
        result.update(source_tree=str(tree_path),vcf_sha256=vcf_sha256,vcf_position_transform='one_based')
        if decoder_mask is not None:
            result.update(simulation_mask_sha256=oo.legacy.digest(decoder_mask),
                recombination_model='cropped-region arithmetic mean; simulation uses full map',
                recombination_to_mutation_ratio=recombination_ratio)
        temporary.cleanup()
        for channel in ('stdout','stderr'):
            (dest/f'decoder.{channel}.log').write_text(result.pop(channel))
        meta = json.loads(Path(str(raw)+'.meta').read_text())
        np.testing.assert_array_equal(meta['pairs'],pairs)
        np.testing.assert_array_equal(meta['output_positions'],positions)
        alpha,beta = posterior_arrays(raw,len(pairs),len(positions))
        counts,sums = summarize_decoded(alpha,beta,classes,cfg)
        frame = build_frame(counts,sums,classes,positions,cfg,'decoded')
        cols=[f'frac_recent_{t}' for t in cfg['tmrca_cutoffs_years']]
        native=pd.read_csv(dest/'native_mean.tsv',sep='\t')
        np.testing.assert_allclose(frame[frame.pair_class=='all pairs'][cols],native[cols],atol=1e-7,rtol=0)
        rows=[]
        for r in frame[frame.pair_class.isin(['all pairs','alt/alt'])].itertuples():
            for cutoff in cfg['tmrca_cutoffs_years']:
                rows.append(dict(source='decoded',method={'all pairs':'all_pairs','alt/alt':'alt_alt'}[r.pair_class],
                    position_0based=r.position_0based,cutoff_years=cutoff,n_pairs=r.n_pairs,
                    score=getattr(r,f'frac_recent_{cutoff}')))
        scores=pd.DataFrame(rows)
        oo.legacy.atomic_json(dest/'decoder_run.json',result)
    scores.to_csv(score_path,index=False)
    oo.legacy.atomic_json(cache_path,dict(identity=identity,scores_sha256=oo.legacy.digest(score_path),
        seconds=time.monotonic()-started,verified_inputs=['decoded_input.trees','focal_carriers.npy']+
        (['trajectory.csv'] if mode=='truth' else [])))
    return dict(task_id=row['task_id'],mode=mode,state='computed',seconds=time.monotonic()-started)

def examples(inventory):
    rows=[]
    for _,g in inventory.groupby(['family','role','s']):
        # AF is a count / 400. Rounding removes floating point noise in exact ties.
        chosen=g.assign(distance=(g.sample_af-g.sample_af.median()).abs().round(12)).sort_values(['distance','task_id']).iloc[0]
        rows.append(chosen.drop(labels='distance').to_dict())
    return pd.DataFrame(rows)

def process(root,out,cfg,inventory,mode,workers):
    chosen=examples(inventory)
    chosen.to_csv(out/'examples.csv',index=False)
    example_ids=set(chosen.task_id)
    pair_hash=oo.legacy.digest(root/'pairs.tsv')
    # Nulls first, then neutral tests, then selected. This also makes partial decode coverage explicit.
    ordered=inventory.assign(priority=inventory.role.map({'null':0,'neutral_target':1,'selected':2})).sort_values(['priority','replicate','family','s'])
    payloads=[(str(root),str(out),r,cfg,pair_hash,mode,r['task_id'] in example_ids) for r in ordered.to_dict('records')]
    results=[]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures=[pool.submit(compute_one,p) for p in payloads]
        for future in as_completed(futures):
            results.append(future.result())
            if len(results)%20==0 or len(results)==len(payloads):
                state=dict(phase=mode,completed=len(results),target=len(payloads),workers=workers,updated=time.time())
                oo.legacy.atomic_json(out/f'{mode}_status.json',state)
                print(json.dumps(state),flush=True)
    oo.legacy.atomic_json(out/f'{mode}_complete.json',dict(results=results,completed=len(results),updated=time.time()))

def evaluate_available(focal,cfg):
    predictions=[]
    keys=['family','source','method','cutoff_years']
    for key,g in focal.groupby(keys):
        null=g[g.role=='null']
        target=g[g.role!='null'].copy()
        if null.empty or target.empty:
            continue
        assert not null.task_id.duplicated().any()
        assert not set(null.task_id)&set(target.task_id)
        target['p']=rank_p(null.score,target.score)
        target['called']=target.p<=cfg['alpha']
        target['evaluable']=target.score.notna()
        target['null_n']=len(null)
        target['null_evaluable']=int(null.score.notna().sum())
        target['p_min']=1/(len(null)+1)
        predictions.append(target)
    if not predictions:
        return pd.DataFrame(),pd.DataFrame()
    predictions=pd.concat(predictions,ignore_index=True)
    rows=[]
    keys=['family','role','s','source','method','cutoff_years']
    for key,g in predictions.groupby(keys):
        n=len(g); calls=int(g.called.sum()); available=int(g.evaluable.sum())
        low,high=wilson(calls,n)
        rows.append(dict(zip(keys,key),n=n,calls=calls,rate=calls/n,ci_low=low,ci_high=high,
            evaluable=available,coverage=available/n,evaluable_rate=calls/available if available else np.nan,
            null_n=int(g.null_n.iloc[0]),null_evaluable=int(g.null_evaluable.iloc[0]),
            p_min=float(g.p_min.iloc[0]),alpha=cfg['alpha'],planned_target_n=cfg['target_replicates'],
            planned_null_n=cfg['null_replicates'],provisional=True,
            endpoint='focal power' if key[1]=='selected' else 'neutral FPR'))
    return predictions,pd.DataFrame(rows)

def assemble(out,cfg,inventory):
    rows=[]
    # Never calibrate an interim report on whichever decoder jobs finish first.
    # Include the decoded source only once the frozen cohort is fully processed.
    sources=['truth']
    complete=out/'decoded_complete.json'
    if complete.exists() and json.loads(complete.read_text())['completed']==len(inventory):
        sources.append('decoded')
    for r in inventory.to_dict('records'):
        common={k:r[k] for k in ['task_id','family','arm','role','s']}
        for source in sources:
            path=out/'regions'/r['task_id']/(source+'.csv')
            receipt=path.with_suffix('.json')
            if not receipt.exists():
                continue
            assert oo.legacy.digest(path)==json.loads(receipt.read_text())['scores_sha256']
            scores=pd.read_csv(path)
            for score in scores[scores.position_0based==cfg['focal_position_bp']].to_dict('records'):
                rows.append(dict(common,**score))
        rows.append(dict(common,source='genotypes',method='AF',position_0based=cfg['focal_position_bp'],
            cutoff_years=0,n_pairs=np.nan,score=r['sample_af']))
    focal=pd.DataFrame(rows)
    predictions,metrics=evaluate_available(focal,cfg)
    focal.to_csv(out/'focal_scores.csv',index=False)
    predictions.to_csv(out/'predictions.csv',index=False)
    metrics.to_csv(out/'results_table.csv',index=False)
    summary=inventory.groupby(['family','role','s']).agg(n=('task_id','size'),
        af_median=('sample_af','median'),af_mean=('sample_af','mean'),af_min=('sample_af','min'),
        af_max=('sample_af','max'),onset_af_median=('onset_af','median'),
        attempts_median=('attempts','median'),attempts_max=('attempts','max'),
        rejected_attempts=('rejected_attempts','sum')).reset_index()
    summary.to_csv(out/'allele_frequency_and_attempts.csv',index=False)
    raw=focal.groupby(['family','role','s','source','method','cutoff_years']).score.agg(
        n='size',available='count',mean='mean',median='median',minimum='min',maximum='max').reset_index()
    raw.to_csv(out/'raw_fraction_summary.csv',index=False)
    diagnostic=[]
    for key,g in focal[focal.role=='null'].groupby(['family','source','method','cutoff_years']):
        values=g.score.dropna()
        n=len(g)
        diagnostic.append(dict(zip(['family','source','method','cutoff_years'],key),
            null_n=n,available=len(values),unavailable=n-len(values),
            score_one_count=int((values==1).sum()),score_one_fraction=float((values==1).sum()/n),
            score_zero_count=int((values==0).sum()),median=float(values.median()),
            q95=float(values.quantile(.95)),p_at_max_possible_score=(1+int((values>=1).sum()))/(n+1),
            p_grid_step=1/(n+1),largest_attainable_significant_p=np.floor(cfg['alpha']*(n+1))/(n+1),
            any_bounded_score_can_reject=(1+int((values>=1).sum()))/(n+1)<=cfg['alpha']))
    pd.DataFrame(diagnostic).to_csv(out/'calibration_diagnostics.csv',index=False)
    return focal,metrics

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--phase',choices=['truth','decoded','report'],required=True)
    parser.add_argument('--workers',type=int,default=4)
    args=parser.parse_args()
    assert 1<=args.workers<=24
    args.out.mkdir(parents=True,exist_ok=True)
    if (args.out/'snapshot.json').exists():
        snapshot=json.loads((args.out/'snapshot.json').read_text())
        cfg=snapshot['config']
        # "null" is a cohort label, not a missing value (pandas' default NA token).
        inventory=pd.read_csv(args.out/'inventory.csv',float_precision='round_trip',keep_default_na=False)
        assert inventory.task_id.tolist()==snapshot['task_ids']
    else:
        cfg,inventory=collect(args.root,args.out)
    if args.phase in ['truth','decoded']:
        process(args.root,args.out,cfg,inventory,args.phase,args.workers)
    else:
        focal,metrics=assemble(args.out,cfg,inventory)
        from .stopped_campaign_plots import make_report
        make_report(args.root,args.out,cfg,inventory,focal,metrics)

if __name__=='__main__':
    main()
