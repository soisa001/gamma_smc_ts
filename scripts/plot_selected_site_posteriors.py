#!/usr/bin/env python3
"""Actual fitted posteriors at the selected site, versus 100 neutral simulations.

prepare reads and verifies raw data; report regenerates figures from portable cache.
The rerun is the first survivor after proposal 3, not selected for inferred signal.
"""
import os
os.environ['MPLBACKEND'] = 'Agg'
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import re
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import PercentFormatter
import numpy as np
import pandas as pd
from scipy.special import gammainc
from scipy.stats import gamma
import tskit
import zstandard
from gamma_smc_aou.fresh_power import ordered_nodes
from plot_fitted_gamma_examples import read_exact
from plot_selection_episode_series import sha, save_json
from report_selection_arrays import local_path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'docs/results/selected_vs_100_neutrals_10k_s001_20261106'
RAW = ROOT.parent / 'sim/selection_50_30k_demo_20261106'
POSITION = 500000
LABELS = ['50–30 kya: original', '50–30 kya: fresh rerun', '10 kya–present: existing example']
WINDOWS = [(30, 50), (30, 50), (0, 10)]
STEMS = ['50_30k_original', '50_30k_rerun', '10k_existing']


def prepare(out, workers):
    inventory = json.loads((SOURCE / 'output_hashes.json').read_text())
    assert sha(SOURCE / 'cohort_manifest.tsv') == inventory['cohort_manifest.tsv']
    old = pd.read_csv(SOURCE / 'cohort_manifest.tsv', sep='\t').to_dict('records')
    assert len(old) == 101 and old[0]['cohort'] == 'selected'
    rerun = json.loads((RAW / 'rerun_00004.json').read_text())
    selected = []
    for record in [json.loads((RAW / '50_30k/introgression/s0.01/selected_00000/complete.json').read_text()), rerun['selected']]:
        selected.append(dict(record, source_directory=str(RAW / record['id'])))
    rows = selected + [old[0]] + old[1:]
    assert len({r['observed_sha256'] for r in rows[3:]}) == 100
    contract = dict(position_bp=POSITION, pair_count=1000, mutation_rate=1.29e-8,
        generation_years=25, source_manifest_sha256=inventory['cohort_manifest.tsv'],
        rerun_receipt_sha256=sha(RAW / 'rerun_00004.json'),
        ids=[r['id'] for r in rows], expected_scores=[r['scores_sha256'] for r in rows],
        neutral_pool='50 held-out neutral_test plus first 50 calibration, same 100 as previous comparison',
        choice='Original 50–30 first survivor; rerun first survivor from proposal 4; existing 10k example previously chosen for a mass call at 420 kb. No new signal ascertainment.',
        schema='selected-site-posteriors/v1')
    if (out / 'manifest.json').exists():
        saved = json.loads((out / 'manifest.json').read_text())
        assert saved['contract'] == contract
        assert sha(out / 'parameters.npz') == saved['cache_sha256']
        return

    def extract(row):
        directory = local_path(row['source_directory'])
        files = {name: directory / name for name in ('posterior.zst', 'posterior.zst.meta', 'observed.trees', 'scores.npz')}
        hashes = {name: sha(path) for name, path in files.items()}
        for name, key in [('posterior.zst', 'posterior_sha256'), ('observed.trees', 'observed_sha256'), ('scores.npz', 'scores_sha256')]:
            assert hashes[name] == row[key], str(files[name])
        meta = json.loads(files['posterior.zst.meta'].read_text())
        positions = np.asarray(meta['output_positions'])
        np.testing.assert_array_equal(positions, np.arange(100) * 10000)
        index = int(np.flatnonzero(positions == POSITION)[0])
        assert meta['num_pairs'] == 1000 and meta['sequence_length'] == 100
        pairs = np.array(meta['pairs'])
        chunk, parts = meta['chunk_size'], []
        size = 100 * chunk * 4
        with files['posterior.zst'].open('rb') as handle, zstandard.ZstdDecompressor().stream_reader(handle) as reader:
            for first in range(0, 1000, chunk):
                n = min(chunk, 1000-first)
                a = np.frombuffer(read_exact(reader, size), np.float32).reshape(100, chunk)[index, :n].astype(float)
                b = np.frombuffer(read_exact(reader, size), np.float32).reshape(100, chunk)[index, :n].astype(float)
                assert np.isfinite(a).all() and np.isfinite(b).all() and a.min() > 0 and b.min() > 0
                parts.append((a, b))
            assert not reader.read(1)
        a = np.concatenate([p[0] for p in parts])
        scale = meta['scaled_mutation_rate'] / (2*1.29e-8) * 25/1000 / np.concatenate([p[1] for p in parts])
        ts = tskit.load(files['observed.trees'])
        nodes = ordered_nodes(ts)
        assert len(nodes) == 200
        for i in range(200):
            assert meta['sample_names'][str(i)] == f'tsk_{i//2}.{i%2}'
        tree = ts.at(POSITION)
        truth = np.array([tree.tmrca(int(nodes[x]), int(nodes[y])) * .025 for x, y in pairs])
        probs = gammainc(a[:, None], np.arange(1, 51) / scale[:, None])
        with np.load(files['scores.npz']) as scored:
            np.testing.assert_allclose(probs.mean(0), scored['cdfs'][2, index], rtol=1e-12, atol=1e-14)
            np.testing.assert_allclose((probs > .5).mean(0), scored['cdfs'][1, index], atol=1e-14)
            np.testing.assert_allclose((truth[:, None] < np.arange(1, 51)).mean(0), scored['cdfs'][0, index], atol=1e-14)
            cdfs = scored['cdfs'][:, index].copy()
        return a, scale, truth, pairs, cdfs, dict(id=row['id'], directory=str(directory), hashes=hashes)

    with ThreadPoolExecutor(max_workers=workers) as pool:
        result = list(pool.map(extract, rows))
    for r in result[1:]:
        np.testing.assert_array_equal(r[3], result[0][3])
    trajectories = []
    for row in rows[:3]:
        path = local_path(row['source_directory']) / 'trajectory.tsv'
        assert sha(path) == row['trajectory_sha256']
        trajectory = np.loadtxt(path, skiprows=1)
        assert 0 < trajectory[-1, 1] < 1
        trajectories.append(trajectory)
    np.savez_compressed(out / 'parameters.npz', alpha=np.stack([r[0] for r in result]),
        scale_kyr=np.stack([r[1] for r in result]), true_kya=np.stack([r[2] for r in result]),
        pairs=result[0][3], cdfs=np.stack([r[4] for r in result]),
        trajectories=np.stack(trajectories), ids=np.array(contract['ids']))
    save_json(out / 'manifest.json', dict(contract=contract, inputs=[r[5] for r in result],
        cache_sha256=sha(out / 'parameters.npz'), selected_records=rows[:3], rerun=rerun,
        validation='Verified source hashes, all 103 posterior streams, identical 1000 pairs, positive parameters, all 50 cutoffs for truth/median/mass at exactly 500000 bp.'))
    print('Verified and cached 103 simulations at 500 kb.', flush=True)


def report(out, workers):
    manifest = json.loads((out / 'manifest.json').read_text())
    assert sha(out / 'parameters.npz') == manifest['cache_sha256']
    with np.load(out / 'parameters.npz') as data:
        alpha, scales, truth, cdfs, trajectories = [data[k].copy() for k in ('alpha', 'scale_kyr', 'true_kya', 'cdfs', 'trajectories')]
    matplotlib.rcParams.update({'font.size': 15, 'axes.labelsize': 16, 'axes.titlesize': 16,
        'xtick.labelsize': 13, 'ytick.labelsize': 13, 'pdf.fonttype': 42, 'ps.fonttype': 42})
    times = np.linspace(.05, 200, 1000)
    def density(i):
        return gamma.pdf(times[:, None], a=alpha[i], scale=scales[i]).mean(1)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        pdfs = np.stack(list(pool.map(density, range(len(alpha)))))
    edges = np.arange(51)
    cumulative = np.concatenate([np.zeros((103, 3, 1)), cdfs], axis=2)
    bins = np.diff(cumulative, axis=2)
    handles = [Line2D([], [], color='.55', lw=1, label='100 neutral simulations'),
               Line2D([], [], color='#0072B2', lw=2.5, label='Neutral mean'),
               Line2D([], [], color='#D55E00', lw=3, label='Selected simulation')]
    def background(ax, x, neutral, selected, window):
        ax.axvspan(*window, color='#F0D8AD', alpha=.5, zorder=0)
        for y in neutral:
            ax.plot(x, y, color='.6', lw=.7, alpha=.38)
        ax.plot(x, neutral.mean(0), color='#0072B2', lw=2.5)
        ax.plot(x, selected, color='#D55E00', lw=3)
        ax.yaxis.set_major_formatter(PercentFormatter(1))
        ax.grid(axis='y', alpha=.2)
        ax.set_ylim(bottom=0)
    def save(fig, name):
        for extension in ('png', 'pdf'):
            fig.savefig(out / f'{name}.{extension}', dpi=180)
        plt.close(fig)
    summaries = []
    for i, (label, window, stem) in enumerate(zip(LABELS, WINDOWS, STEMS)):
        fig, axes = plt.subplots(2, 1, figsize=(11, 8.5))
        background(axes[0], times, pdfs[3:], pdfs[i], window)
        axes[0].set(xlim=(0, 200), xlabel='Pairwise TMRCA (kya)', ylabel='Posterior density\n(% per kyr)')
        background(axes[1], edges, cumulative[3:, 2], cumulative[i, 2], window)
        axes[1].set(xlim=(0, 50), xlabel='TMRCA cutoff (kya)', ylabel='Cumulative posterior\ncoalescence fraction')
        fig.suptitle(f'{label} · s=0.01 · selected site (500 kb)', y=.98, fontsize=18)
        fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.5, .94), ncol=3, frameon=False, fontsize=12)
        fig.text(.5, .025, 'EAS + 1% introgression · 1,000 pairs/simulation · shaded: selection active · no renormalization', ha='center', fontsize=11)
        fig.subplots_adjust(top=.84, bottom=.12, left=.16, right=.97, hspace=.6)
        save(fig, f'{stem}_posterior')
        fig, axes = plt.subplots(3, 1, figsize=(11, 8.5), sharex=True)
        for method, ax in enumerate(axes):
            background(ax, edges[:-1]+.5, bins[3:, method], bins[i, method], window)
            ax.set_ylabel('Pair fraction')
            ax.set_title(['True TMRCA', 'Posterior median', 'Posterior-weighted mass'][method], loc='left', fontsize=14)
            ax.set_xlim(0, 50)
        axes[-1].set_xlabel('Pairwise TMRCA: 1-kyr intervals (kya)')
        fig.suptitle(f'{label} · s=0.01 · selected site (500 kb)', y=.98, fontsize=18)
        fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.5, .94), ncol=3, frameon=False, fontsize=12)
        fig.text(.5, .025, 'Same 1,000 pairs for all methods · 50 held-out + 50 calibration neutrals · descriptive comparison', ha='center', fontsize=11)
        fig.subplots_adjust(top=.82, bottom=.12, left=.13, right=.97, hspace=.57)
        save(fig, f'{stem}_intervals')
        intervals = {}
        for lower, upper in [window, (0, 50)]:
            values = cumulative[:, :, upper] - cumulative[:, :, lower]
            intervals[f'{lower}_{upper}_kya'] = {method: dict(selected=float(values[i,j]),
                neutral_mean=float(values[3:,j].mean()), neutral_95_quantile=float(np.quantile(values[3:,j],.95)),
                neutral_ge_selected=int(np.sum(values[3:,j]>=values[i,j])),
                descriptive_rank_p=float((1+np.sum(values[3:,j]>=values[i,j]))/101))
                for j,method in enumerate(('truth','median','mass'))}
        summaries.append(dict(label=label, id=manifest['contract']['ids'][i], intervals=intervals,
            af={str(t):float(trajectories[i, trajectories[i,:,0]==t*1000,1][0]) for t in [50,30,10,0]},
            posterior_tail_above_200=float(gamma.sf(200,a=alpha[i],scale=scales[i]).mean())))
    fig, ax = plt.subplots(figsize=(11, 8.5))
    for i, color in enumerate(('#D55E00','#009E73','#CC79A7')):
        ax.plot(trajectories[i,:,0]/1000, trajectories[i,:,1], color=color, lw=2.5, label=LABELS[i])
    ax.set(xlim=(50,0), ylim=(0,1), xlabel='Time before present (kya)', ylabel='Selected-allele frequency')
    ax.yaxis.set_major_formatter(PercentFormatter(1))
    ax.grid(axis='y',alpha=.2)
    fig.suptitle('Allele-frequency trajectories · all selected-site examples', y=.97, fontsize=18)
    fig.legend(loc='upper center',bbox_to_anchor=(.5,.92),frameon=False,ncol=1,fontsize=13)
    fig.text(.5,.045,'s=0.01 · whole-population frequencies · conditioned on present-day segregation',ha='center',fontsize=12)
    fig.subplots_adjust(top=.73,bottom=.14,left=.13,right=.97)
    save(fig,'allele_frequencies')
    for path in out.glob('*.pdf'):
        content=path.read_bytes()
        assert b'/FontFile2' in content and not re.search(rb'/Subtype\s*/Image',content)
    save_json(out/'summary.json',dict(examples=summaries,position_bp=POSITION,neutral_count=100,
        caveats='Descriptive fixed-site comparisons, not genome-wide adjusted selection calls or estimates of power. 10k example previously ascertained at 420 kb. 50 held-out plus 50 calibration neutrals.',
        script_sha256=sha(Path(__file__)),versions=dict(numpy=np.__version__,matplotlib=matplotlib.__version__)))
    (out/'.gitattributes').write_text('* -text\n',encoding='utf-8')
    save_json(out/'output_hashes.json',{p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file() and p.name!='output_hashes.json'})
    print(json.dumps(summaries,indent=2),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase',choices=['prepare','report','all'])
    parser.add_argument('--out',type=Path,default=ROOT/'docs/results/selected_site_posteriors_20261106')
    parser.add_argument('--workers',type=int,default=4)
    args=parser.parse_args()
    assert 1<=args.workers<=24
    args.out.mkdir(parents=True,exist_ok=True)
    if args.phase in ('prepare','all'):prepare(args.out,args.workers)
    if args.phase in ('report','all'):report(args.out,args.workers)


if __name__=='__main__':main()
