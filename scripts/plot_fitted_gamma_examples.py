#!/usr/bin/env python3
"""Plot actual Gamma-SMC pair posteriors from selected and neutral simulations.

Phase 1 verifies raw decodes and caches parameters and tree truth at 420 kb.
Phase 2 plots paired examples near 10/100 kya and the average over all 1000 pairs.
No posterior parameters are invented, refit, or rescaled to match target ages.
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
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.ticker import PercentFormatter
import numpy as np
import pandas as pd
from scipy.stats import gamma
from scipy.special import gammainc, gammaincinv
import tskit
import zstandard
from gamma_smc_aou.fresh_power import ordered_nodes
from plot_selection_episode_series import sha, save_json
from report_selection_arrays import local_path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'docs/results/selected_vs_100_neutrals_10k_s001_20261106'
COLORS = ('#D55E00', '#0072B2')
COHORTS = ('Selected', 'Neutral')


def read_exact(reader, size):
    result = bytearray()
    while len(result) < size:
        block = reader.read(size - len(result))
        if not block:
            raise EOFError('Truncated posterior file')
        result.extend(block)
    return result


def prepare(out, workers):
    contract = dict(version=1, source_manifest_sha256=sha(SOURCE / 'cohort_manifest.tsv'),
                    source_curves_sha256=sha(SOURCE / 'curves.npz'), position_bp=420000,
                    mutation_rate=1.29e-8, generation_years=25,
                    choice='Selected example and first held-out neutral; nearest posterior median or true TMRCA to each target, ties by pair order.')
    inventory = json.loads((SOURCE / 'output_hashes.json').read_text())
    assert contract['source_manifest_sha256'] == inventory['cohort_manifest.tsv']
    assert contract['source_curves_sha256'] == inventory['curves.npz']
    manifest, cache = out / 'manifest.json', out / 'parameters.npz'
    if cache.exists() and manifest.exists():
        saved = json.loads(manifest.read_text())
        assert saved['contract'] == contract and saved['cache_sha256'] == sha(cache)
        with np.load(cache) as data:
            return {k: data[k].copy() for k in data.files}, saved
    rows = pd.read_csv(SOURCE / 'cohort_manifest.tsv', sep='\t').iloc[:2]
    assert rows.cohort.tolist() == ['selected', 'neutral_test']

    def extract(row):
        directory = local_path(row.source_directory)
        paths = {name: directory / name for name in ('posterior.zst', 'posterior.zst.meta', 'observed.trees', 'scores.npz')}
        hashes = {name: sha(path) for name, path in paths.items()}
        for name, expected in [('posterior.zst', row.posterior_sha256), ('observed.trees', row.observed_sha256), ('scores.npz', row.scores_sha256)]:
            assert hashes[name] == expected, f'Changed input: {paths[name]}'
        meta = json.loads(paths['posterior.zst.meta'].read_text())
        positions = np.asarray(meta['output_positions'])
        np.testing.assert_array_equal(positions, np.arange(100) * 10000)
        index = int(np.flatnonzero(positions == 420000)[0])
        assert meta['num_pairs'] == 1000 and meta['sequence_length'] == 100
        pairs = np.asarray(meta['pairs'])
        assert pairs.shape == (1000, 2)
        for i in range(200):
            assert meta['sample_names'][str(i)] == f'tsk_{i//2}.{i%2}'
        chunk = meta['chunk_size']
        size = len(positions) * chunk * 4
        parts = []
        with paths['posterior.zst'].open('rb') as handle, zstandard.ZstdDecompressor().stream_reader(handle) as reader:
            for first in range(0, 1000, chunk):
                n = min(chunk, 1000 - first)
                a = np.frombuffer(read_exact(reader, size), dtype=np.float32).reshape(100, chunk)[index, :n].astype(float)
                b = np.frombuffer(read_exact(reader, size), dtype=np.float32).reshape(100, chunk)[index, :n].astype(float)
                assert np.isfinite(a).all() and np.isfinite(b).all() and a.min() > 0 and b.min() > 0
                parts.append((a, b))
            assert not reader.read(1), 'Trailing posterior bytes'
        alpha = np.concatenate([x[0] for x in parts])
        beta = np.concatenate([x[1] for x in parts])
        scale_kyr = meta['scaled_mutation_rate'] / (2 * 1.29e-8) * 25 / 1000 / beta
        ts = tskit.load(paths['observed.trees'])
        nodes = ordered_nodes(ts)
        assert len(nodes) == 200
        tree = ts.at(420000)
        truth = np.array([tree.tmrca(int(nodes[a]), int(nodes[b])) * 25 / 1000 for a, b in pairs])
        cutoffs = np.arange(1, 51)
        probability = gammainc(alpha[:, None], cutoffs[None, :] / scale_kyr[:, None])
        with np.load(paths['scores.npz']) as scores:
            np.testing.assert_allclose(probability.mean(axis=0), scores['cdfs'][2, index], rtol=1e-12, atol=1e-14)
            np.testing.assert_allclose((probability > .5).mean(axis=0), scores['cdfs'][1, index], atol=1e-14)
            np.testing.assert_allclose((truth[:, None] < cutoffs).mean(axis=0), scores['cdfs'][0, index], atol=1e-14)
        return alpha, scale_kyr, truth, pairs, dict(id=row.id, directory=str(directory), hashes=hashes,
                                                  theta=meta['scaled_mutation_rate'], rho=meta['scaled_recombination_rate'])

    with ThreadPoolExecutor(max_workers=workers) as pool:
        result = list(pool.map(extract, rows.itertuples()))
    arrays = dict(alpha=np.stack([r[0] for r in result]), scale_kyr=np.stack([r[1] for r in result]),
                  true_kya=np.stack([r[2] for r in result]), pairs=np.stack([r[3] for r in result]), ids=rows.id.to_numpy(str))
    with np.load(SOURCE / 'curves.npz') as source:
        for j in range(2):
            reconstructed = gammainc(arrays['alpha'][j, :, None], np.arange(1, 51) / arrays['scale_kyr'][j, :, None]).mean(axis=0)
            np.testing.assert_allclose(reconstructed, source['cdfs'][1, j], rtol=1e-12, atol=1e-14)
    np.savez_compressed(cache, **arrays)
    saved = dict(contract=contract, inputs=[r[4] for r in result], cache_sha256=sha(cache),
                 validation='Reconstructed all 50 stored mass and median cutoffs and sampled-pair truth at 420 kb for both simulations. Confirmed pair ordering, stream length, positive parameters, and original file hashes.')
    save_json(manifest, saved)
    return arrays, saved


def render(data, out):
    matplotlib.rcParams.update({'font.size': 15, 'axes.labelsize': 16, 'axes.titlesize': 16,
                               'xtick.labelsize': 13, 'ytick.labelsize': 13, 'pdf.fonttype': 42, 'ps.fonttype': 42})
    alpha, scales = data['alpha'], data['scale_kyr']
    medians = gammaincinv(alpha, .5) * scales
    np.testing.assert_allclose(gammainc(alpha, medians / scales), .5, atol=1e-12)
    records = []
    for basis, ages in [('posterior_median', medians), ('true_tmrca', data['true_kya'])]:
        for target in (10, 100):
            for j, cohort in enumerate(COHORTS):
                i = int(np.argmin(abs(ages[j] - target)))
                a, scale = float(alpha[j, i]), float(scales[j, i])
                low, high = gamma.ppf([.025, .975], a=a, scale=scale)
                records.append(dict(basis=basis, target_kya=target, cohort=cohort, simulation=str(data['ids'][j]),
                                    pair_index=i, haplotype_a=int(data['pairs'][j, i, 0]), haplotype_b=int(data['pairs'][j, i, 1]),
                                    alpha=a, scale_kyr=scale, mean_kya=a * scale, median_kya=float(medians[j, i]),
                                    mode_kya=max(a - 1, 0) * scale, true_kya=float(data['true_kya'][j, i]),
                                    credible_low_kya=float(low), credible_high_kya=float(high)))
    pd.DataFrame(records).to_csv(out / 'pair_examples.tsv', sep='\t', index=False, lineterminator='\n')
    with PdfPages(out / 'fitted_gamma_curves.pdf') as book:
        def save(fig, name):
            fig.savefig(out / f'{name}.png', dpi=180)
            fig.savefig(out / f'{name}.pdf')
            book.savefig(fig)
            plt.close(fig)

        fig, ax = plt.subplots(figsize=(11, 8.5))
        times = np.linspace(.1, 200, 2400)
        for j, cohort in enumerate(COHORTS):
            densities = gamma.pdf(times[:, None], a=alpha[j], scale=scales[j]).mean(axis=1)
            ax.plot(times, densities, lw=3, color=COLORS[j], label=cohort)
        ax.axvspan(0, 10, color='.9', zorder=0)
        ax.set(xlim=(0, 200), ylim=(0, None), xlabel='Pairwise TMRCA (kya)', ylabel='Mean posterior density (% per kyr)')
        ax.yaxis.set_major_formatter(PercentFormatter(1))
        ax.grid(axis='y', alpha=.2)
        fig.suptitle('Selected versus neutral: actual Gamma-SMC posteriors', y=.97, fontsize=18)
        fig.text(.5, .92, 'One simulation per curve · 1,000 pairwise posteriors averaged · position 420 kb', ha='center', fontsize=13)
        fig.legend(loc='upper center', bbox_to_anchor=(.5, .89), ncol=2, frameon=False)
        tails = [float(gamma.sf(200, a=alpha[j], scale=scales[j]).mean()) for j in range(2)]
        fig.text(.5, .075, 'Selected: s=0.01 from 10 kya to present (shaded). Neutral: first held-out replicate.', ha='center', fontsize=12)
        fig.text(.5, .04, f'Probability beyond 200 kya: selected {tails[0]:.1%}; neutral {tails[1]:.1%}. Curves are not renormalized.', ha='center', fontsize=12)
        fig.subplots_adjust(top=.81, bottom=.18, left=.12, right=.97)
        save(fig, '01_selected_vs_neutral')
        for number, basis in enumerate(('posterior_median', 'true_tmrca'), start=2):
            fig, axes = plt.subplots(2, 1, figsize=(11, 8.5))
            basis_label = 'inferred median' if basis == 'posterior_median' else 'true TMRCA'
            fig.suptitle(f'Individual gamma posteriors: pairs chosen by {basis_label}', y=.98, fontsize=17)
            fig.text(.5, .935, 'Actual fitted parameters · same two simulations and position · nearest pair to each target', ha='center', fontsize=12)
            handles = [Line2D([], [], color=COLORS[j], lw=3, label=COHORTS[j]) for j in range(2)]
            handles += [Line2D([], [], color='.3', ls='--', label='Inferred median'), Line2D([], [], color='.3', ls=':', label='True TMRCA')]
            fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.51, .91), ncol=4, frameon=False, fontsize=12)
            for target, ax in zip((10, 100), axes):
                group = [r for r in records if r['basis'] == basis and r['target_kya'] == target]
                xmax = np.ceil(max(max(r['credible_high_kya'], r['true_kya']) for r in group) / 10) * 10
                times = np.linspace(.05, xmax, 3000)
                for j, record in enumerate(group):
                    ax.plot(times, gamma.pdf(times, a=record['alpha'], scale=record['scale_kyr']), color=COLORS[j], lw=3)
                    ax.axvline(record['median_kya'], color=COLORS[j], ls='--', lw=1.5)
                    ax.axvline(record['true_kya'], color=COLORS[j], ls=':', lw=1.5)
                ax.set_title(f"Target {basis_label}: {target} kya", loc='left')
                ax.set(xlim=(0, xmax), ylim=(0, None), xlabel='Pairwise TMRCA (kya)', ylabel='Posterior density\n(% per kyr)')
                ax.yaxis.set_major_formatter(PercentFormatter(1))
                ax.grid(axis='y', alpha=.18)
                labels = [f"{r['cohort']}: median {r['median_kya']:.1f}; true {r['true_kya']:.1f} kya" for r in group]
                ax.text(.99, .94, '\n'.join(labels), transform=ax.transAxes, ha='right', va='top', fontsize=11.5,
                        bbox=dict(facecolor='white', edgecolor='none', alpha=.85))
            fig.text(.5, .04, 'Pairs were chosen by age proximity, not curve shape. X/y scales differ between panels; older tails retained.', ha='center', fontsize=11)
            fig.subplots_adjust(top=.79, bottom=.13, left=.12, right=.98, hspace=.6)
            save(fig, f'0{number}_pairs_by_{basis}')
    pdf_checks = {}
    for path in out.glob('*.pdf'):
        content = path.read_bytes()
        pdf_checks[path.name] = dict(embedded_truetype=b'/FontFile2' in content,
                                   raster_objects=len(re.findall(rb'/Subtype\s*/Image', content)))
        assert pdf_checks[path.name]['embedded_truetype'] and pdf_checks[path.name]['raster_objects'] == 0
    save_json(out / 'pdf_validation.json', pdf_checks)
    save_json(out / 'summary.json', dict(examples=records, average_posterior_tail_above_200_kya=dict(zip(COHORTS, tails)),
                                       script_sha256=sha(Path(__file__)), versions=dict(numpy=np.__version__, matplotlib=matplotlib.__version__)))
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=ROOT / 'docs/results/fitted_gamma_examples_20261106')
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    assert 1 <= args.workers <= 24
    args.out.mkdir(parents=True, exist_ok=True)
    data, proof = prepare(args.out, args.workers)
    records = render(data, args.out)
    (args.out / '.gitattributes').write_text('* -text\n', encoding='utf-8')
    save_json(args.out / 'provenance.json', dict(workers=args.workers, randomness='No RNG draws or new simulations; deterministic nearest-age examples.',
        interpretation='Actual original Gamma-SMC gamma parameters, transformed from coalescent units into kyr using theta/(2*mu)*25/1000. Selected region ascertained from earlier mass call; neutral first held-out replicate. Individual pairs need not share genotype class or true age when chosen by inferred median. Average curve is a mixture, not itself necessarily gamma.',
        density_grid='Positive lower endpoint 0.05 kya (individual) or 0.1 kya (mixture), avoiding possible gamma singularity at zero; no truncation renormalization.',
        input_validation=proof['validation']))
    save_json(args.out / 'output_hashes.json', {p.name: sha(p) for p in sorted(args.out.iterdir()) if p.is_file() and p.name != 'output_hashes.json'})
    print(pd.DataFrame(records)[['basis', 'target_kya', 'cohort', 'alpha', 'scale_kyr', 'median_kya', 'true_kya', 'credible_low_kya', 'credible_high_kya']].to_string(index=False), flush=True)


if __name__ == '__main__':
    main()
