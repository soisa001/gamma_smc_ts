"""Restartable, site-level archaic-allele null calibration across populations."""
from __future__ import annotations
import os
os.environ['MPLBACKEND'] = 'Agg'
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'
import argparse
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from contextlib import redirect_stdout, redirect_stderr
import importlib.metadata
import json
from pathlib import Path
import shutil
import time
import traceback
import numpy as np
import pandas as pd
from . import origin_onset as oo
from . import segregating_introgression as si
from . import stopped_campaign_report as scoring
from .origin_onset_analysis import rank_p

VERSION = 'all-population-neutral-site/v1'
POPULATIONS = ('AFR', 'AMR', 'EAS', 'EUR', 'MID', 'SAS')


def configs(study):
    assert set(study['populations']) == set(POPULATIONS)
    result = {}
    for pop in POPULATIONS:
        cfg = dict(study['parameters'], **study['populations'][pop], population=pop)
        assert cfg['selection_onsets_years'] == [50000] and cfg['pulse_years'] == 50000
        assert cfg['selection_coefficients'] == [] and cfg['target_replicates'] == 0
        assert cfg['introgression_proportion'] == .02
        assert cfg['focal_ascertainment'] == si.SPEC
        assert set(cfg['scaling_by_origin'].values()) == {1}
        assert cfg['null_replicates'] > 0
        result[pop] = cfg
    return result


def jobs(pop_configs):
    by_pop = {pop: si.tasks(cfg) for pop, cfg in pop_configs.items()}
    assert len({len(tasks) for tasks in by_pop.values()}) == 1
    return [(pop, by_pop[pop][rep]) for rep in range(len(next(iter(by_pop.values()))))
            for pop in POPULATIONS]


def analysis_identity(cfg, runtime, task, decoder_hash):
    return dict(schema=VERSION, simulation=oo.simulation_identity(cfg, task, runtime),
        decoder=decoder_hash, analysis={k: v for k, v in cfg.items() if k.startswith('decoder_') or k in
        ('pairs_seed', 'haplotype_pairs', 'tmrca_cutoffs_years', 'recent_call')})


def score_identity(cfg, row, pairs_hash, source, decoder_hash):
    result = dict(schema=scoring.VERSION, simulation=row['simulation_sha256'], pairs=pairs_hash,
        positions=[cfg['focal_position_bp']], mode=source, cutoffs=cfg['tmrca_cutoffs_years'])
    if source == 'decoded':
        result.update(decoder=decoder_hash, parameters={k: v for k, v in cfg.items()
            if k.startswith('decoder_') or k in ('generation_time_years', 'mutation_rate', 'recent_call')})
    return result


def reuse_scores(archive, dest, cfg, row, pairs_hash, source, decoder_hash):
    old = archive/'regions'/row['task_id']/(source+'.json')
    if not old.exists():
        return False
    receipt = json.loads(old.read_text())
    expected = score_identity(cfg, row, pairs_hash, source, decoder_hash)
    if {k:v for k,v in receipt['identity'].items() if k != 'positions'} != {
            k:v for k,v in expected.items() if k != 'positions'}:
        return False
    scores = old.with_suffix('.csv')
    if oo.legacy.digest(scores) != receipt['scores_sha256']:
        raise ValueError(f'Corrupt cached score file: {scores}')
    frame = pd.read_csv(scores, float_precision='round_trip')
    frame = frame[frame.position_0based == cfg['focal_position_bp']]
    if len(frame) != 2*len(cfg['tmrca_cutoffs_years']):
        raise ValueError('Archive lacks the complete focal statistic set')
    dest.mkdir(parents=True, exist_ok=True)
    target = dest/(source+'.csv')
    frame.to_csv(target, index=False)
    oo.legacy.atomic_json(target.with_suffix('.json'), dict(identity=expected,
        scores_sha256=oo.legacy.digest(target), reused_from=str(scores),
        source_receipt_sha256=oo.legacy.digest(old)))
    return True


def run_job(payload):
    root_s, study, cfg, task, runtime, slim, decoder_hash, phase = payload
    root = Path(root_s)
    pop = cfg['population']
    pop_root = root/pop
    directory = pop_root/'regions'/task['id']
    directory.mkdir(parents=True, exist_ok=True)
    final = directory/'calibration.json'
    identity = analysis_identity(cfg, runtime, task, decoder_hash)
    started = time.time()
    if final.exists() and phase != 'simulate':
        record = json.loads(final.read_text())
        if record['identity'] != identity:
            raise ValueError(f'Incompatible completed calibration: {final}')
        for name, digest in record['scores_sha256'].items():
            assert oo.legacy.digest(pop_root/'scores'/'regions'/task['id']/name) == digest
        return dict(population=pop, task_id=task['id'], state='cached', record=str(final))
    import resource
    # New workers plus the paused older campaign stay within the 200-GB budget.
    budget = cfg.get('worker_memory_budget_gb', cfg['memory_limit_gb'])
    limit = int(budget*1e9/cfg['workers'])
    resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
    free = shutil.disk_usage(root)
    if free.free < 20_000_000_000 or free.used >= cfg['storage_limit_bytes']:
        raise RuntimeError('Storage reserve/budget reached')
    with (directory/'calibration.stdout.log').open('a') as stdout, (directory/'calibration.stderr.log').open('a') as stderr:
        with redirect_stdout(stdout), redirect_stderr(stderr):
            sim_identity = identity['simulation']
            source_root = pop_root
            old_root = Path(study['reuse_eas_root'])
            if pop == 'EAS' and (old_root/'regions'/task['id']/'simulation.json').exists():
                old_manifest = json.loads((old_root/'manifest.json').read_text())
                if oo.simulation_identity(old_manifest['config'], task, old_manifest['runtime']) == sim_identity:
                    source_root = old_root
            source_dir = source_root/'regions'/task['id']
            receipt = source_dir/'simulation.json'
            if phase != 'decode' and source_root == pop_root:
                si.simulate(cfg, task, directory, sim_identity, slim)
            sim = oo.read_receipt(source_dir, 'simulation.json', sim_identity)
            if sim is None:
                raise ValueError(f'Missing simulation for decode: {source_dir}')
            oo.validate_trees(source_dir, cfg, sim)
            assert sim['seed'] == oo.seed_for(cfg, task, sim['attempts'][-1]['attempt'])
            pairs_hash = oo.legacy.digest(pop_root/'pairs.tsv')
            assert oo.legacy.digest(source_root/'pairs.tsv') == pairs_hash
            oo.legacy.atomic_json(directory/'simulation_source.json', dict(source_root=str(source_root),
                receipt_sha256=oo.legacy.digest(receipt), population=pop, task=task,
                reused=source_root != pop_root))
            if phase == 'simulate':
                return dict(population=pop, task_id=task['id'], state='simulated_or_reused')
            row = dict(task_id=task['id'], sample_af=sim['sample_af'], simulation_sha256=oo.legacy.digest(receipt))
            dest = pop_root/'scores'/'regions'/task['id']
            for source in ('truth', 'decoded'):
                if not (dest/(source+'.json')).exists() and source_root != pop_root:
                    reuse_scores(Path(study['reuse_eas_report']), dest, cfg, row, pairs_hash, source, decoder_hash)
                scoring.compute_one((str(source_root), str(pop_root/'scores'), row, cfg, pairs_hash, source, False))
            scores = pd.concat([pd.read_csv(dest/(source+'.csv'), float_precision='round_trip')
                                for source in ('truth', 'decoded')], ignore_index=True)
            assert len(scores) == 4*len(cfg['tmrca_cutoffs_years'])
            values = [{k: (None if pd.isna(v) else v) for k,v in r.items()} for r in scores.to_dict('records')]
            oo.legacy.atomic_json(final, dict(identity=identity, population=pop, task=task,
                source_root=str(source_root), simulation_sha256=row['simulation_sha256'],
                sample_af=sim['sample_af'], onset_af=sim['focal_choice']['onset_af'],
                attempts=len(sim['attempts']), scores=values,
                scores_sha256={source+'.csv': oo.legacy.digest(dest/(source+'.csv')) for source in ('truth','decoded')},
                completed=time.time(), seconds=time.time()-started))
    return dict(population=pop, task_id=task['id'], state='complete', record=str(final))


def critical_summary(values, alpha):
    """Exact conservative upper-tail boundary, retaining null no-calls in n."""
    values = np.asarray(values, dtype=float)
    assert len(values) and not np.isinf(values).any()
    reference = np.sort(np.where(np.isnan(values), -np.inf, values))
    allowed = np.flatnonzero((1+np.arange(len(values)+1))/(len(values)+1) <= alpha)
    boundary = float(reference[-(int(allowed[-1])+1)]) if len(allowed) else np.inf
    finite = values[np.isfinite(values)]
    return dict(null_n=len(values), available=len(finite), unavailable=len(values)-len(finite),
        critical_score=boundary, decision='score > critical_score; missing score is no call',
        q95_among_available=float(np.quantile(finite, .95)) if len(finite) else np.nan,
        p_at_score_one=float(rank_p(values, [1])[0]),
        rejection_possible_for_bounded_score=bool(rank_p(values, [1])[0] <= alpha),
        p_grid_step=1/(len(values)+1), alpha=alpha)


def summarize(root, pop_configs, figures=False):
    rows = []
    counts = []
    for pop, cfg in pop_configs.items():
        n = 0
        for task in si.tasks(cfg):
            path = root/pop/'regions'/task['id']/'calibration.json'
            if not path.exists():
                continue
            r = json.loads(path.read_text())
            n += 1
            common = dict(population=pop, replicate=task['replicate'], task_id=task['id'])
            rows.extend(dict(common, **score) for score in r['scores'])
            rows.append(dict(common, source='genotypes', method='AF', cutoff_years=0,
                position_0based=cfg['focal_position_bp'], n_pairs=0, score=r['sample_af']))
        counts.append(dict(population=pop, completed=n, planned=cfg['null_replicates']))
    table = pd.DataFrame(counts)
    table.to_csv(root/'cohort_counts.csv', index=False)
    if not rows:
        return table
    scores = pd.DataFrame(rows)
    scores.to_csv(root/'null_scores.csv', index=False)
    cuts = []
    for key, g in scores.groupby(['population', 'source', 'method', 'cutoff_years']):
        cfg = pop_configs[key[0]]
        cuts.append(dict(zip(['population', 'source', 'method', 'cutoff_years'], key),
            **critical_summary(g.score, cfg['alpha']), planned_null_n=cfg['null_replicates'],
            provisional=len(g) < cfg['null_replicates']))
    cutoffs = pd.DataFrame(cuts)
    cutoffs.to_csv(root/'pointwise_cutoffs.csv', index=False)
    if figures:
        from .population_neutral_plots import plot
        plot(root, scores, cutoffs, table)
    return table


def initialise(root, study, slim):
    import subprocess
    result = subprocess.run([slim, '-v'], capture_output=True, text=True, check=True)
    assert 'SLiM version 5.' in result.stdout+result.stderr
    pop_configs = configs(study)
    cfg = study['parameters']
    assert 1 <= cfg['workers'] <= 20 and 0 < cfg['memory_limit_gb'] <= 200
    runtime = dict(versions={p:importlib.metadata.version(p) for p in ('msprime','numpy','pyslim','stdpopsim','tskit')},
                   slim_sha256=oo.legacy.digest(slim))
    decoder_hash = oo.legacy.digest(oo.REPO/'bin/gamma_smc')
    scientific = {pop: {k:v for k,v in c.items() if k not in ('workers','memory_limit_gb','worker_memory_budget_gb','storage_limit_bytes')}
                  for pop,c in pop_configs.items()}
    scientific.update(runtime=runtime, decoder_sha256=decoder_hash, schema=VERSION)
    manifest_path = root/'manifest.json'
    if manifest_path.exists():
        assert json.loads(manifest_path.read_text())['scientific_identity'] == scientific
    possible = np.column_stack(np.triu_indices(2*cfg['sample_diploids'], k=1))
    choices = np.random.default_rng(cfg['pairs_seed']).choice(len(possible), cfg['haplotype_pairs'], replace=False)
    pairs = possible[np.sort(choices)]
    for pop, c in pop_configs.items():
        artifact = si.load_phlash_npz(oo.REPO/c['phlash_resource'], expected_sha256=c['phlash_sha256'], expected_population=pop)
        pop_root = root/pop
        pop_root.mkdir(exist_ok=True)
        if (pop_root/'pairs.tsv').exists():
            np.testing.assert_array_equal(np.loadtxt(pop_root/'pairs.tsv', dtype=int, ndmin=2), pairs)
        else:
            np.savetxt(pop_root/'pairs.tsv', pairs, fmt='%d', delimiter='\t')
        oo.legacy.atomic_json(pop_root/'demography.json', dict(population=pop, hash=c['phlash_sha256'],
            time_generations=[0.]+artifact.time_generations.tolist(),
            median_ne=[float(artifact.median_ne[0])]+artifact.median_ne.tolist(),
            mutation_rate=c['mutation_rate'], extrapolation='constant endpoint Ne; original generation scale retained'))
    plan = jobs(pop_configs)
    assert len({oo.seed_for(pop_configs[p], t) for p,t in plan}) == len(plan)
    oo.legacy.atomic_json(manifest_path, dict(scientific_identity=scientific, study=study,
        created=time.time(), jobs=len(plan), conditional_on='archaic-derived allele segregating at 50 kya and observed in present sample; fixation retained',
        testing='prespecified focal site; separate p<=0.05 tests for each TMRCA cutoff and pair class',
        implementation_sha256={name:oo.legacy.digest(Path(__file__).with_name(name)) for name in
            ('population_neutral.py','segregating_introgression.py','eas_sweep_models.py','origin_onset.py','stopped_campaign_report.py')}))
    return pop_configs, runtime, decoder_hash


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--slim', required=True)
    parser.add_argument('--phase', choices=['plan','simulate','decode','report','run'], default='plan')
    parser.add_argument('--replicates', type=int, help='Bounded validation in a separate output directory')
    args = parser.parse_args()
    study = json.loads(args.config.read_text(encoding='utf-8-sig'))
    if args.replicates is not None:
        assert args.replicates > 0
        study['parameters']['null_replicates'] = args.replicates
    args.out.mkdir(parents=True, exist_ok=True)
    import fcntl
    with (args.out/'study.lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX|fcntl.LOCK_NB)
        pop_configs, runtime, decoder_hash = initialise(args.out, study, args.slim)
        plan = jobs(pop_configs)
        if args.phase == 'plan':
            print(json.dumps(dict(jobs=len(plan), populations=list(pop_configs), output=str(args.out))))
            return
        if args.phase == 'report':
            summarize(args.out, pop_configs, figures=True)
            return
        # Validate/reuse already completed EAS nulls before scheduling new work.
        # Task seeds depend on their identities, never on completion/submission order.
        archive = Path(study['reuse_eas_root'])
        plan.sort(key=lambda item: not (item[0] == 'EAS' and
            (archive/'regions'/item[1]['id']/'simulation.json').exists()))
        results, active, failure = [], {}, False
        remaining = iter(plan)
        workers = study['parameters']['workers']
        def status(state):
            oo.legacy.atomic_json(args.out/'run_status.json', dict(state=state, phase=args.phase,
                pid=os.getpid(), workers=workers, memory_limit_gb=study['parameters']['memory_limit_gb'],
                worker_memory_budget_gb=study['parameters'].get('worker_memory_budget_gb', study['parameters']['memory_limit_gb']),
                target=len(plan), completed=len(results), failed=[r for r in results if r['state']=='failed'],
                active=list(active.values()), updated=time.time()))
        with ProcessPoolExecutor(max_workers=workers) as pool:
            def submit():
                item = next(remaining, None)
                if item is not None:
                    pop, task = item
                    payload = (str(args.out), study, pop_configs[pop], task, runtime, args.slim, decoder_hash, args.phase)
                    active[pool.submit(run_job, payload)] = f'{pop}/{task["id"]}'
            for _ in range(workers): submit()
            status('running')
            last_summary = 0
            while active:
                done, _ = wait(active, timeout=30, return_when=FIRST_COMPLETED)
                for future in done:
                    name = active.pop(future)
                    try:
                        result = future.result()
                    except Exception:
                        result = dict(task_id=name, state='failed', error=traceback.format_exc())
                        failure = True
                    results.append(result)
                    print(json.dumps(result), flush=True)
                if not failure:
                    while len(active) < workers:
                        n = len(active)
                        submit()
                        if len(active) == n: break
                status('finishing_after_failure' if failure else 'running')
                if len(results)-last_summary >= 20:
                    summarize(args.out, pop_configs)
                    last_summary = len(results)
        summarize(args.out, pop_configs, figures=True)
        status('failed' if failure else 'complete')
        if failure:
            raise RuntimeError('Calibration task failed; see run_status.json and per-region logs')


if __name__ == '__main__':
    main()
