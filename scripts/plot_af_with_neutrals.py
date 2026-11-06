#!/usr/bin/env python3
"""Plot selected and matched neutral focal-allele trajectories at 500 kb.

prepare verifies and caches all 100 neutral trajectories; report needs only the
committed caches. Neutral trajectories track the original introduced allele,
including loss, not the nearest present-day segregating neutral variant.
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
from plot_selection_episode_series import sha, save_json
from report_selection_arrays import local_path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'docs/results/selected_vs_100_neutrals_10k_s001_20261106'


def prepare(out, workers):
    source_hashes = json.loads((SOURCE / 'output_hashes.json').read_text())
    assert sha(SOURCE / 'cohort_manifest.tsv') == source_hashes['cohort_manifest.tsv']
    rows = pd.read_csv(SOURCE / 'cohort_manifest.tsv', sep='\t').iloc[1:]
    assert len(rows) == 100 and rows.id.nunique() == 100
    parent = json.loads((out / 'manifest.json').read_text())
    assert rows.id.tolist() == parent['contract']['ids'][3:]
    contract = dict(position_bp=500000, source_manifest_sha256=source_hashes['cohort_manifest.tsv'],
        ids=rows.id.tolist(), trajectory_sha256=rows.trajectory_sha256.tolist(),
        tracked_allele='Original focal allele introduced on Neanderthal donor haplotypes at 50 kya; not the relabelled nearest segregating neutral site.',
        conditioning='Neutral unconditioned, including losses/fixation; selected examples conditioned on present-day segregation.',
        neutral_pool='50 held-out plus 50 calibration; same simulations and order as posterior comparison.')
    manifest_path, cache = out / 'neutral_af_manifest.json', out / 'neutral_af.npz'
    if manifest_path.exists():
        saved = json.loads(manifest_path.read_text())
        assert saved['contract'] == contract and sha(cache) == saved['cache_sha256']
        return

    def read(row):
        path = local_path(row.source_directory) / 'trajectory.tsv'
        assert sha(path) == row.trajectory_sha256, str(path)
        data = np.loadtxt(path, skiprows=1)
        assert data.shape == (2001, 5) and np.isfinite(data).all()
        np.testing.assert_array_equal(data[:, 0], np.arange(50000, -1, -25))
        np.testing.assert_allclose(data[:, 1], data[:, 4] / (2 * data[:, 2]), atol=1e-14)
        assert np.all((data[:, 1] >= 0) & (data[:, 1] <= 1))
        np.testing.assert_array_equal(data[:, 3], 0)
        np.testing.assert_allclose(data[[0,-1], 1], [row.initial_frequency, row.final_frequency], atol=1e-14)
        return data, dict(id=row.id, path=str(path), sha256=row.trajectory_sha256)

    with ThreadPoolExecutor(max_workers=workers) as pool:
        results = list(pool.map(read, rows.itertuples()))
    trajectories = np.stack([r[0] for r in results])
    np.savez_compressed(cache, trajectories=trajectories, ids=rows.id.to_numpy(str))
    save_json(manifest_path, dict(contract=contract, inputs=[r[1] for r in results],
        cache_sha256=sha(cache), validation='All 100 hashes, time grids, zero selection coefficients, AF bounds, exact copy-count frequencies and receipt endpoints verified.'))


def report(out):
    parent = json.loads((out / 'manifest.json').read_text())
    neutral_manifest = json.loads((out / 'neutral_af_manifest.json').read_text())
    assert sha(out / 'parameters.npz') == parent['cache_sha256']
    assert sha(out / 'neutral_af.npz') == neutral_manifest['cache_sha256']
    with np.load(out / 'parameters.npz') as data:
        selected = data['trajectories'].copy()
    with np.load(out / 'neutral_af.npz') as data:
        neutral = data['trajectories'].copy()
        np.testing.assert_array_equal(data['ids'], parent['contract']['ids'][3:])
    np.testing.assert_array_equal(neutral[:, :, 0], np.broadcast_to(selected[0, :, 0], (100, 2001)))
    np.testing.assert_array_equal(neutral[:, :, 2], np.broadcast_to(selected[0, :, 2], (100, 2001)))
    matplotlib.rcParams.update({'font.size':16, 'axes.labelsize':16, 'axes.titlesize':16,
        'xtick.labelsize':14, 'ytick.labelsize':14, 'pdf.fonttype':42, 'ps.fonttype':42})
    times = neutral[0, :, 0] / 1000
    fig, axes = plt.subplots(2, 1, figsize=(11, 8.5), sharex=True, sharey=True)
    colors = ['#D55E00', '#009E73', '#CC79A7']
    for ax, window, indices, title in zip(axes, [(30,50), (0,10)], [(0,1), (2,)],
            ['Selection from 50 to 30 kya', 'Selection from 10 kya to present']):
        ax.axvspan(*window, color='#F0D8AD', alpha=.5, zorder=0)
        for trajectory in neutral:
            ax.plot(times, trajectory[:,1], color='.6', lw=.8, alpha=.45)
        ax.plot(times, neutral[:,:,1].mean(0), color='#0072B2', lw=2.5)
        for i in indices:
            ax.plot(times, selected[i,:,1], color=colors[i], lw=2.7)
        ax.set(xlim=(50,0), ylim=(0,1), ylabel='Focal-allele frequency', title=title)
        ax.yaxis.set_major_formatter(PercentFormatter(1))
        ax.grid(axis='y',alpha=.2)
    axes[-1].set_xlabel('Time before present (kya)')
    handles = [Line2D([],[],color='.6',lw=1,label='100 neutrals'),
        Line2D([],[],color='#0072B2',lw=2.5,label='Neutral mean')]
    handles += [Line2D([],[],color=c,lw=2.7,label=l) for c,l in zip(colors,
        ['50–30 kya: original','50–30 kya: rerun','<10 kya: existing'])]
    fig.suptitle('Allele-frequency trajectories at the focal site (500 kb)',y=.98,fontsize=19)
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.935),ncol=3,frameon=False,fontsize=12)
    fig.text(.5,.055,'EAS + 1% introgression at 50 kya · selected s=0.01 · shaded: selection active',ha='center',fontsize=12)
    fig.text(.5,.025,'Neutral: original introduced allele, including losses · selected: present-day segregation required',ha='center',fontsize=11)
    fig.subplots_adjust(top=.79,bottom=.15,left=.13,right=.97,hspace=.36)
    for extension in ('png','pdf'):
        fig.savefig(out / f'allele_frequencies_with_neutrals.{extension}',dpi=180)
    plt.close(fig)
    pdf = (out / 'allele_frequencies_with_neutrals.pdf').read_bytes()
    assert b'/FontFile2' in pdf and not re.search(rb'/Subtype\s*/Image',pdf)
    final = neutral[:,-1,1]
    summary = dict(neutral_count=100,lost=int(np.sum(final==0)),fixed=int(np.sum(final==1)),
        segregating=int(np.sum((final>0)&(final<1))),mean_present_frequency=float(final.mean()),
        max_present_frequency=float(final.max()), selected_present_frequencies=selected[:,-1,1].tolist(),
        script_sha256=sha(Path(__file__)), tracking=neutral_manifest['contract']['tracked_allele'])
    save_json(out / 'neutral_af_summary.json',summary)
    save_json(out / 'output_hashes.json',{p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file() and p.name!='output_hashes.json'})
    print(json.dumps(summary,indent=2),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase',choices=['prepare','report','all'])
    parser.add_argument('--out',type=Path,default=ROOT/'docs/results/selected_site_posteriors_20261106')
    parser.add_argument('--workers',type=int,default=4)
    args=parser.parse_args()
    assert 1<=args.workers<=24
    if args.phase in ('prepare','all'):prepare(args.out,args.workers)
    if args.phase in ('report','all'):report(args.out)


if __name__=='__main__':main()
