#!/usr/bin/env python3
"""Simulate survival-conditioned 50-30-kya selection examples.

Phases: run (bounded four-worker batches until the survivor target) and report (cache
only, no raw simulations needed). Reuses the established EAS SLiM simulator and
Gamma-SMC scorer. Fixed original RNG scheme; select by proposal order, not signal.
"""
import os
os.environ['MPLBACKEND'] = 'Agg'
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import ast
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
from pathlib import Path
import re
import shutil
import subprocess
import time
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np
from scipy.stats import gamma
from scipy.special import gammainc
import zstandard
import run_selection_arrays as arrays
from plot_fitted_gamma_examples import read_exact
from plot_selection_episode_series import sha, save_json

ROOT = Path(__file__).resolve().parents[1]
PRIOR = ROOT / 'docs/results/fitted_gamma_examples_20261106'
BASE_REPORT = ROOT / 'docs/results/selection_arrays_eas_20261105_4h'


def run(args):
    args.raw.mkdir(parents=True, exist_ok=True)
    original = json.loads((BASE_REPORT / 'manifest.json').read_text())['contract']['configs']['50k']
    cfg = dict(original, selection_start_years=50000, selection_end_years=30000,
               seed_namespace='episode-v1:50_30k:')
    native_lib = Path(cfg['slim']).parent.parent / 'lib'
    os.environ['LD_LIBRARY_PATH'] = str(native_lib) + (':' + os.environ['LD_LIBRARY_PATH'] if os.environ.get('LD_LIBRARY_PATH') else '')
    linked = subprocess.run(['ldd', cfg['decoder']], capture_output=True, text=True, check=True)
    assert 'not found' not in linked.stdout, linked.stdout
    save_json(args.raw / 'execution_environment.json', dict(native_library_directory=str(native_lib),
        libhts_sha256=sha(native_lib / 'libhts.so.3'), decoder_linkage=linked.stdout))
    for key, path in [('source_sha256', arrays.__file__), ('ancestry_source_sha256', arrays.ancestor.__file__),
                      ('score_source_sha256', arrays.base.__file__), ('slim_sha256', cfg['slim']),
                      ('decoder_sha256', cfg['decoder']), ('phlash_sha256', ROOT / 'no_introgression_EAS_sim/resources/EAS.npz')]:
        current_hash = sha(Path(path))
        if key == 'source_sha256' and current_hash != original[key]:
            snapshot = ROOT.parent / 'sim/selection_arrays_20261104/executed_sources/run_selection_arrays.py'
            assert sha(snapshot) == original[key]
            def scientific_functions(source):
                return {node.name: ast.dump(node, include_attributes=False) for node in ast.parse(Path(source).read_text()).body
                        if isinstance(node, ast.FunctionDef) and node.name in ('timing', 'simulate', 'work')}
            assert scientific_functions(path) == scientific_functions(snapshot), 'Simulation implementation changed'
            cfg['historical_source_sha256'] = original[key]
            cfg[key] = current_hash
            cfg['source_compatibility'] = 'Simulation/timing/worker AST exactly matches historical executed source; scheduler-only changes allowed.'
        else:
            assert current_hash == original[key], f'Historical simulation dependency changed: {path}'
    contract = dict(config=cfg, selection=.01, selection_rule='First accepted selected replicate in proposal order; no ascertainment by AF or inferred signal.',
                    null='First existing held-out introgression neutral; no focal survival conditioning.',
                    position_bp=420000, schema='single-episode-demo/v1')
    manifest = args.raw / 'manifest.json'
    if manifest.exists():
        assert json.loads(manifest.read_text())['contract'] == contract
    else:
        save_json(manifest, dict(contract=contract, python=__import__('sys').version,
                                versions={n: __import__(n).__version__ for n in ('numpy', 'scipy', 'msprime', 'pyslim', 'tskit')}))
    possible = np.column_stack(np.triu_indices(200, 1))
    pairs = possible[np.sort(np.random.default_rng(1729).choice(len(possible), 1000, replace=False))]
    for name, values in [('pairs.tsv', pairs), ('positions.txt', np.arange(0, 1000000, 10000))]:
        path = args.raw / name
        if path.exists():
            np.testing.assert_array_equal(np.loadtxt(path, dtype=int), values)
        else:
            np.savetxt(path, values, fmt='%d', delimiter='\t')
    snapshots = args.raw / 'executed_sources'
    snapshots.mkdir(exist_ok=True)
    for filename in (arrays.__file__, arrays.ancestor.__file__, arrays.base.__file__, __file__):
        dest = snapshots / Path(filename).name
        if dest.exists() and sha(dest) != sha(Path(filename)):
            shutil.copy2(filename, snapshots / f'{Path(filename).stem}_{sha(Path(filename))[:12]}.py')
        elif not dest.exists():
            shutil.copy2(filename, dest)
    deadline = time.time() + args.seconds
    receipts, accepted = [], None
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for first in range(args.start_replicate, args.start_replicate + args.attempts, args.workers):
            if time.time() >= deadline:
                break
            tasks = []
            for i in range(first, min(first + args.workers, args.start_replicate + args.attempts)):
                unit = f'introgression/s0.01/selected_{i:05d}'
                tasks.append(dict(id='50_30k/' + unit, unit=unit, origin='introgression', cohort='selected',
                                  selection=.01, replicate=i, regime='50_30k'))
            futures = {pool.submit(arrays.work, (str(args.raw), cfg, task, deadline)): task for task in tasks}
            for future in as_completed(futures):
                result = future.result()
                task = futures[future]
                path = args.raw / task['id'] / ('complete.json' if result['state'] == 'accepted' else 'rejected.json')
                record = json.loads(path.read_text()) if path.exists() else dict(task, state=result['state'])
                receipts.append(record)
                print(f"{task['id']}: {result['state']}; AF={record.get('final_frequency')}", flush=True)
            receipts.sort(key=lambda r: r['replicate'])
            ledger = args.raw / 'proposal_ledger.json'
            previous = json.loads(ledger.read_text()) if ledger.exists() else []
            combined = {r['id']: r for r in previous + receipts}
            save_json(ledger, sorted(combined.values(), key=lambda r: r['replicate']))
            hits = [r for r in receipts if r.get('accepted')]
            save_json(args.raw / 'cohort_progress.json', dict(target=args.target, accepted=len(hits),
                completed=len(receipts), lost=sum(r.get('rejection')=='loss' for r in receipts),
                fixed=sum(r.get('rejection')=='fixation' for r in receipts),
                last_proposal=max(r['replicate'] for r in receipts), updated_epoch=time.time(),
                deadline_epoch=deadline, workers=args.workers))
            if len(hits) >= args.target:
                accepted = hits[0]
                break
    if accepted is None:
        raise RuntimeError('Survivor target not reached within bounded proposal budget; proposal ledger retained for resume.')
    if args.simulation_only:
        filename = f'rerun_{args.start_replicate:05d}.json' if args.target==1 else f'cohort_{args.target:05d}.json'
        save_json(args.raw / filename, dict(selected=accepted, selected_cohort=hits[:args.target],
            proposals=receipts, start_replicate=args.start_replicate, workers=args.workers,
            seconds=args.seconds, attempts=args.attempts, target=args.target, script_sha256=sha(Path(__file__))))
        return
    directory = args.raw / accepted['id']
    trajectory = np.loadtxt(directory / 'trajectory.tsv', skiprows=1)
    np.testing.assert_array_equal(trajectory[:, 0], np.arange(50000, -1, -25))
    np.testing.assert_allclose(trajectory[:, 3], .01 * (trajectory[:, 0] >= 30000))
    assert 0 < trajectory[-1, 1] < 1
    raw = directory / 'posterior.zst'
    assert sha(raw) == accepted['posterior_sha256']
    assert sha(directory / 'scores.npz') == accepted['scores_sha256']
    meta = json.loads(Path(str(raw) + '.meta').read_text())
    np.testing.assert_array_equal(meta['pairs'], pairs)
    np.testing.assert_array_equal(meta['output_positions'], np.arange(100) * 10000)
    chunk = meta['chunk_size']
    size = 100 * chunk * 4
    parts = []
    with raw.open('rb') as handle, zstandard.ZstdDecompressor().stream_reader(handle) as reader:
        for first in range(0, 1000, chunk):
            n = min(chunk, 1000 - first)
            alpha = np.frombuffer(read_exact(reader, size), dtype=np.float32).reshape(100, chunk)[42, :n].astype(float)
            beta = np.frombuffer(read_exact(reader, size), dtype=np.float32).reshape(100, chunk)[42, :n].astype(float)
            assert np.isfinite(alpha).all() and np.isfinite(beta).all() and alpha.min() > 0 and beta.min() > 0
            parts.append((alpha, beta))
        assert not reader.read(1)
    alpha = np.concatenate([p[0] for p in parts])
    scale = meta['scaled_mutation_rate'] / (2 * 1.29e-8) * 25 / 1000 / np.concatenate([p[1] for p in parts])
    with np.load(directory / 'scores.npz') as scores:
        reconstructed = gammainc(alpha[:, None], np.arange(1, 51) / scale[:, None]).mean(axis=0)
        np.testing.assert_allclose(reconstructed, scores['cdfs'][2, 42], rtol=1e-12, atol=1e-14)
    old_hashes = json.loads((PRIOR / 'output_hashes.json').read_text())
    assert sha(PRIOR / 'parameters.npz') == old_hashes['parameters.npz']
    with np.load(PRIOR / 'parameters.npz') as old:
        np.testing.assert_array_equal(old['pairs'][1], pairs)
        alpha = np.stack([alpha, old['alpha'][1]])
        scale = np.stack([scale, old['scale_kyr'][1]])
        neutral_id = str(old['ids'][1])
    args.out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.out / 'parameters.npz', alpha=alpha, scale_kyr=scale, trajectory=trajectory)
    through_first = [r for r in receipts if r['replicate'] <= accepted['replicate']]
    save_json(args.out / 'manifest.json', dict(contract=contract, raw_root=str(args.raw), selected=accepted,
        neutral_id=neutral_id, neutral_parameters_sha256=old_hashes['parameters.npz'],
        selected_metadata_sha256=sha(Path(str(raw) + '.meta')), parameters_sha256=sha(args.out / 'parameters.npz'),
        proposals_through_first_survivor=len(through_first), total_completed_in_batch=len(receipts),
        rejected_before_first=[r.get('rejection') for r in through_first if not r.get('accepted')],
        computation_budget=dict(seconds=args.seconds, attempts=args.attempts, workers=args.workers,
                                note='No new batch after deadline; an active batch is allowed to finish.'),
        validation='Exact selection schedule and whole-population segregation checked. Original RNG, pairs, and demographic/binary hashes retained. Decoded gamma parameters reproduce all 50 saved posterior-mass cutoffs at 420 kb. Neutral unchanged.'))
    save_json(args.out / 'proposal_ledger.json', receipts)
    shutil.copy2(directory / 'trajectory.tsv', args.out / 'trajectory.tsv')


def report(args):
    meta = json.loads((args.out / 'manifest.json').read_text())
    assert sha(args.out / 'parameters.npz') == meta['parameters_sha256']
    with np.load(args.out / 'parameters.npz') as data:
        alpha, scales, trajectory = data['alpha'], data['scale_kyr'], data['trajectory']
    matplotlib.rcParams.update({'font.size': 15, 'axes.labelsize': 16, 'axes.titlesize': 16,
                               'xtick.labelsize': 13, 'ytick.labelsize': 13, 'pdf.fonttype': 42, 'ps.fonttype': 42})
    colors = ['#D55E00', '#0072B2']
    times = np.linspace(.1, 200, 2400)
    fig, ax = plt.subplots(figsize=(11, 8.5))
    for j, label in enumerate(('Selected: 50–30 kya', 'Neutral')):
        ax.plot(times, gamma.pdf(times[:, None], a=alpha[j], scale=scales[j]).mean(axis=1), color=colors[j], lw=3, label=label)
    ax.axvspan(30, 50, color='.9', zorder=0)
    ax.set(xlim=(0, 200), ylim=(0, None), xlabel='Pairwise TMRCA (kya)', ylabel='Mean posterior density (% per kyr)')
    ax.yaxis.set_major_formatter(PercentFormatter(1))
    ax.grid(axis='y', alpha=.2)
    fig.suptitle('Selected versus neutral: selection at 50–30 kya', y=.97, fontsize=18)
    fig.text(.5, .92, 's=0.01 · EAS + 1% introgression at 50 kya · 1,000 pairwise posteriors · position 420 kb', ha='center', fontsize=12)
    fig.legend(loc='upper center', bbox_to_anchor=(.5, .89), ncol=2, frameon=False)
    tails = [float(gamma.sf(200, a=alpha[j], scale=scales[j]).mean()) for j in range(2)]
    fig.text(.5, .075, f"First surviving proposal: {meta['selected']['replicate']}; present AF {trajectory[-1,1]:.1%}. Same neutral as earlier figure.", ha='center', fontsize=12)
    fig.text(.5, .04, f'Probability beyond 200 kya: selected {tails[0]:.1%}; neutral {tails[1]:.1%}. No renormalization.', ha='center', fontsize=12)
    fig.subplots_adjust(top=.81, bottom=.18, left=.12, right=.97)
    for extension in ('png', 'pdf'):
        fig.savefig(args.out / f'01_selected_vs_neutral.{extension}', dpi=180)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(11, 8.5))
    ax.axvspan(30, 50, color='.9', label='Selection active: 50–30 kya')
    ax.plot(trajectory[:, 0]/1000, trajectory[:, 1], color=colors[0], lw=3, label='Selected allele')
    ax.set(xlim=(50, 0), ylim=(0, 1), xlabel='Time before present (kya)', ylabel='Whole-population allele frequency')
    ax.yaxis.set_major_formatter(PercentFormatter(1))
    ax.grid(axis='y', alpha=.2)
    fig.suptitle('Allele-frequency trajectory of the new selected example', y=.96, fontsize=17)
    fig.legend(loc='upper center', bbox_to_anchor=(.5, .91), ncol=2, frameon=False, fontsize=13)
    fig.text(.5, .045, 'First surviving proposal in fixed seed order; conditioned on present-day segregation, not inferred signal.', ha='center', fontsize=11)
    fig.subplots_adjust(top=.82, bottom=.15, left=.12, right=.97)
    for extension in ('png', 'pdf'):
        fig.savefig(args.out / f'02_allele_frequency.{extension}', dpi=180)
    plt.close(fig)
    for path in args.out.glob('*.pdf'):
        content = path.read_bytes()
        assert b'/FontFile2' in content and not re.search(rb'/Subtype\s*/Image', content)
    ages = (50, 30, 0)
    summary = dict(selected_id=meta['selected']['id'], neutral_id=meta['neutral_id'],
                   af={str(t): float(trajectory[trajectory[:,0] == t*1000,1][0]) for t in ages},
                   posterior_mass_30_50_kya={label: float((gamma.cdf(50,a=alpha[j],scale=scales[j])-gamma.cdf(30,a=alpha[j],scale=scales[j])).mean()) for j,label in enumerate(('selected','neutral'))},
                   tail_above_200_kya=dict(zip(('selected','neutral'), tails)),
                   interpretation='Single illustrative survivor versus unchanged neutral. Does not estimate power or selection dates. Averaged fitted posteriors are mixtures, not necessarily gamma.',
                   script_sha256=sha(Path(__file__)), versions=dict(numpy=np.__version__, matplotlib=matplotlib.__version__))
    save_json(args.out / 'summary.json', summary)
    (args.out / '.gitattributes').write_text('* -text\n', encoding='utf-8')
    save_json(args.out / 'output_hashes.json', {p.name: sha(p) for p in sorted(args.out.iterdir()) if p.is_file() and p.name != 'output_hashes.json'})
    print(json.dumps(summary, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=('run', 'report'))
    parser.add_argument('--raw', type=Path, default=ROOT.parent / 'sim/selection_50_30k_demo_20261106')
    parser.add_argument('--out', type=Path, default=ROOT / 'docs/results/selection_50_30k_demo_20261106')
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--seconds', type=int, default=900)
    parser.add_argument('--attempts', type=int, default=64)
    parser.add_argument('--start-replicate', type=int, default=0)
    parser.add_argument('--simulation-only', action='store_true', help='Skip historical 420-kb report; use plot_selected_site_posteriors.py for 500-kb figures.')
    parser.add_argument('--target', type=int, default=1, help='Number of accepted simulations in proposal order; use --simulation-only for cohorts.')
    args = parser.parse_args()
    assert 1 <= args.workers <= 24 and args.seconds > 0 and args.attempts > 0
    assert args.target > 0 and args.start_replicate >= 0
    assert args.target == 1 or args.simulation_only
    if args.phase == 'run':
        run(args)
    if not args.simulation_only:
        report(args)


if __name__ == '__main__':
    main()
