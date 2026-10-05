"""Seeded, unclipped draws from the archived low-rank log-Ne MVN, with diagnostics."""
from __future__ import annotations

import os
os.environ['MPLBACKEND'] = 'Agg'
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '4'
import argparse
import hashlib
import json
from pathlib import Path
import zipfile
import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams.update({'pdf.fonttype': 42, 'ps.fonttype': 42, 'font.size': 15,
                           'axes.spines.top': False, 'axes.spines.right': False})
import matplotlib.pyplot as plt
from gamma_smc_aou.eas_sweep_models import load_phlash_npz

REPO = Path(__file__).resolve().parents[1]
POPS = ('AFR', 'AMR', 'EAS', 'EUR', 'MID', 'SAS')


def digest(path):
    return hashlib.file_digest(Path(path).open('rb'), 'sha256').hexdigest()


def draw(artifact, n, seed):
    rng = np.random.default_rng(seed)
    values = artifact.mean_log_ne + rng.normal(size=(n, artifact.covariance_factor.shape[0])) @ artifact.covariance_factor
    if artifact.jitter:
        values += rng.normal(scale=artifact.jitter, size=values.shape)
    values = np.exp(values)
    if not np.all(np.isfinite(values) & (values > 0)):
        raise ValueError('MVN generated invalid Ne; do not clip or silently redraw')
    return values


def panel(ax, years, original, generated, title):
    for curves, color, label, linestyle in ((original, '#555555', 'Archived fits', '--'),
                                          (generated, '#1476c9', 'MVN draws', '-')):
        q = np.quantile(curves, [.025, .5, .975], axis=0)
        ax.fill_between(years, q[0], q[2], color=color, alpha=.13)
        ax.plot(years, q[1], color=color, lw=2.2, ls=linestyle, label=label + ': median + 95% band')
        ax.plot(years, q[0], color=color, lw=.9, ls=linestyle)
        ax.plot(years, q[2], color=color, lw=.9, ls=linestyle)
    ax.set(xscale='log', yscale='log', title=title, xlabel='Years ago', ylabel='Diploid effective size')
    ax.grid(alpha=.18)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True, type=Path)
    parser.add_argument('--draws', type=int, default=1000)
    parser.add_argument('--seed', type=int, default=20380101)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    config = json.loads((REPO/'configs/all_population_neutral.json').read_text())
    fig, axes = plt.subplots(3, 2, figsize=(15, 15), layout='constrained')
    metrics, manifest = [], {'seed_base': args.seed, 'draws_per_population': args.draws,
        'generation_time_years': 25, 'model': 'log_Ne = mean_log_ne + z @ covariance_factor + jitter * epsilon',
        'clipping': False, 'uncertainty': 'pointwise 2.5/97.5 percentiles across fit-level histories',
        'sources': {}, 'files': {}}
    for index, pop in enumerate(POPS):
        spec = config['populations'][pop]
        artifact = load_phlash_npz(REPO/spec['phlash_resource'], expected_sha256=spec['phlash_sha256'], expected_population=pop)
        seed = args.seed + index
        path = args.out/f'{pop}_draws.npz'
        identity = dict(source_sha256=artifact.actual_sha256, seed=seed, draws=args.draws, algorithm='numpy.default_rng/low_rank_log_mvn/v1')
        receipt = path.with_suffix('.json')
        if receipt.exists():
            cached = json.loads(receipt.read_text())
            assert cached['identity'] == identity and digest(path) == cached['sha256'], 'Incompatible or corrupt draws'
            with np.load(path, allow_pickle=False) as d:
                generated = d['ne']
        else:
            generated = draw(artifact, args.draws, seed)
            temp = path.with_suffix('.tmp')
            with temp.open('wb') as stream:
                np.savez_compressed(stream, time_generations=artifact.time_generations, ne=generated, seed=seed, population=pop)
            assert zipfile.ZipFile(temp).testzip() is None
            temp.replace(path)
            receipt.write_text(json.dumps(dict(identity=identity, sha256=digest(path)), indent=2)+'\n')
        # Prefix stability and RNG reproducibility are verified independently of the cache.
        np.testing.assert_allclose(generated[:3], draw(artifact, args.draws if artifact.jitter else 3, seed)[:3], rtol=0, atol=1e-9)
        orig = artifact.bootstrap_ne
        oq = np.quantile(orig, [.025, .5, .975], axis=0)
        dq = np.quantile(generated, [.025, .5, .975], axis=0)
        analytic = np.exp(artifact.mean_log_ne[None, :] + np.array([-1.959963984540054, 0, 1.959963984540054])[:, None] *
                          np.sqrt(np.sum(artifact.covariance_factor**2, axis=0) + artifact.jitter**2)[None, :])
        table = {'time_generations': artifact.time_generations, 'years_ago': artifact.time_generations*25}
        for i, quantile in enumerate(('q025','median','q975')):
            table.update({f'original_{quantile}': oq[i], f'draw_{quantile}': dq[i], f'analytic_mvn_{quantile}': analytic[i]})
            error = np.abs(dq[i]/oq[i]-1)
            fit_error = np.abs(analytic[i]/oq[i]-1)
            metrics.append(dict(population=pop, quantile=quantile, median_relative_difference=float(np.median(error)),
                p95_relative_difference=float(np.quantile(error,.95)), maximum_relative_difference=float(error.max()),
                analytic_mvn_median_relative_difference=float(np.median(fit_error)),
                min_draw_ne=float(generated.min()), max_draw_ne=float(generated.max())))
        table['original_mean'] = orig.mean(axis=0)
        table['draw_mean'] = generated.mean(axis=0)
        pd.DataFrame(table).to_csv(args.out/f'{pop}_quantiles.csv', index=False)
        panel(axes.flat[index], artifact.time_generations*25, orig, generated, pop)
        single, ax = plt.subplots(figsize=(11, 8.5), layout='constrained')
        panel(ax, artifact.time_generations*25, orig, generated, f'{pop}: 100 archived fits vs {args.draws:,} MVN histories')
        ax.legend(loc='lower center', bbox_to_anchor=(.5,1.02), ncol=2, frameon=False)
        single.savefig(args.out/f'{pop}_history_validation.png', dpi=180)
        single.savefig(args.out/f'{pop}_history_validation.pdf')
        plt.close(single)
        manifest['sources'][pop] = identity
        print(json.dumps(dict(population=pop, state='draws_and_diagnostics_complete', **identity)), flush=True)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='outside upper center', ncol=2, frameon=False)
    fig.savefig(args.out/'history_validation_all.png', dpi=150)
    fig.savefig(args.out/'history_validation_all.pdf')
    plt.close(fig)
    pd.DataFrame(metrics).to_csv(args.out/'history_validation_metrics.csv', index=False)
    manifest['files'] = {p.name: digest(p) for p in sorted(args.out.iterdir()) if p.is_file() and p.name != 'manifest.json'}
    (args.out/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')


if __name__ == '__main__':
    main()
