#!/usr/bin/env python3
"""Small matched-pair ASMC/Gamma-SMC comparison; reuse existing simulated data.

Run prepare, decode, and report as separate, restartable phases. ASMC runs in
sequence mode with CSFS disabled. No new simulations or random pair draws.
"""
import os
os.environ['MPLBACKEND'] = 'Agg'
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
import hashlib
import importlib.metadata
import json
from pathlib import Path
import time
import shutil

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MU, REC, GEN = 1.29e-8, 1e-8, 25
POSITIONS = np.arange(0, 1000000, 10000)
SOURCE = ROOT / 'docs/results/selection_episode_40_20k_s002_20261106/cohort_manifest.tsv'
METHODS = ['truth', 'gamma_median', 'gamma_mass', 'asmc_constant', 'asmc_eas']


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')


def local(path):
    s = str(path).replace('\\', '/')
    if os.name != 'nt' and len(s) > 2 and s[1:3] == ':/':
        s = '/mnt/' + s[0].lower() + s[2:]
    return Path(s)


def prepare(a):
    import tskit
    work = a.work.resolve(); work.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(SOURCE.open(), delimiter='\t'))
    rows = ([r for r in rows if r['cohort'] == 'selected'][:a.replicates] +
            [r for r in rows if r['cohort'] == 'neutral_test'][:a.replicates])
    assert len(rows) == 2 * a.replicates
    contracts = []
    for i, row in enumerate(rows):
        directory = local(row['source_directory'])
        for name, key in [('observed.trees', 'observed_sha256'), ('scores.npz', 'scores_sha256')]:
            assert sha(directory / name) == row[key]
        meta = json.loads((directory / 'posterior.zst.meta').read_text())
        assert np.isclose(meta['scaled_mutation_rate'], .00075, rtol=1e-8)
        assert np.isclose(meta['scaled_recombination_rate'], .00075*REC/MU, rtol=1e-6)
        np.testing.assert_array_equal(meta['output_positions'], POSITIONS)
        pairs = np.array(meta['pairs'], dtype=int)
        assert pairs.shape == (1000, 2)
        if i == 0:
            np.savetxt(work / 'pairs.tsv', pairs, fmt='%d')
        else:
            np.testing.assert_array_equal(pairs, np.loadtxt(work / 'pairs.tsv', dtype=int))
        d = work / f'region_{i:02d}'; d.mkdir(exist_ok=True)
        ts = tskit.load(directory / 'observed.trees')
        sample_set = set(map(int, ts.samples()))
        groups = [[int(n) for n in ind.nodes if n in sample_set] for ind in ts.individuals()]
        groups = [g for g in groups if g]
        assert len(groups) == 100 and all(len(g) == 2 for g in groups)
        nodes = np.array([n for g in groups for n in g])
        assert set(nodes) == sample_set
        records = {}; skipped = 0
        for v in ts.variants(samples=nodes):
            # Match Gamma-SMC's VCF reader: it skips non-biallelic sites.
            if len(v.alleles) != 2:
                skipped += 1
                continue
            assert v.genotypes.min() >= 0
            assert v.site.position == int(v.site.position)
            records[int(v.site.position)] = v.genotypes.copy()
        for pos in POSITIONS:
            records.setdefault(int(pos), np.zeros(200, dtype=np.int32))
        records.setdefault(999999, np.zeros(200, dtype=np.int32))
        with (d / 'input.hap').open('w') as h, (d / 'input.map').open('w') as m:
            for pos, g in sorted(records.items()):
                h.write(f'1 p{pos} {pos+1} A C ' + ' '.join(map(str, g)) + '\n')
                m.write(f'1\tp{pos}\t{pos*REC*100:.12g}\t{pos+1}\n')
        (d / 'input.samples').write_text('ID_1 ID_2 missing\n0 0 0\n' +
                                       ''.join(f'1 s{i} 0\n' for i in range(100)))
        np.save(d / 'marker_positions.npy', np.array(sorted(records)))
        with np.load(directory / 'scores.npz') as old:
            assert old['cdfs'].shape == (3, 100, 50)
            # Independently check sample ordering and focal truth against tree.
            ages = np.array([ts.at(500000).tmrca(nodes[x], nodes[y]) * GEN for x, y in pairs])
            expected = (ages[:, None] < np.arange(1, 51) * 1000).mean(0)
            np.testing.assert_allclose(expected, old['cdfs'][0, 50], atol=1e-14)
            np.savez_compressed(d / 'baseline.npz', cdfs=old['cdfs'], true_ages=ages)
        contracts.append(dict(id=row['id'], cohort=row['cohort'], source=str(directory), skipped_non_biallelic=skipped,
                              observed_sha256=row['observed_sha256'], scores_sha256=row['scores_sha256'],
                              input_sha256={name: sha(d / name) for name in
                                            ['input.hap', 'input.map', 'input.samples', 'baseline.npz']}))
    # Neutral allele frequencies only; CSFS/array ascertainment are unused in decoding.
    with (work / 'neutral.frq').open('w') as f:
        f.write('CHR SNP A1 A2 MAF NCHROBS\n')
        for i in range(a.replicates, 2*a.replicates):
            for line in (work / f'region_{i:02d}/input.hap').read_text().splitlines():
                v = line.split(); p = sum(map(int, v[5:])) / 200
                if 0 < p < 1:
                    f.write(f'1 r{i}_{v[1]} A C {min(p,1-p):.12g} 200\n')
    # Preserve cumulative coalescent hazard at each demographic boundary.
    with np.load(ROOT / 'no_introgression_EAS_sim/resources/EAS.npz') as z:
        n = np.median(z['bootstrap_ne'], axis=0).astype(float)
        t = np.r_[0., z['time'].astype(float)]; n = np.r_[n[0], n]
    H = np.r_[0., np.cumsum(np.diff(t) / (2*n[:-1]))]
    edges = np.unique(np.r_[0., np.arange(40, 2001, 40), np.geomspace(100, 40000, 256)])
    idx = np.searchsorted(t, edges, side='right') - 1
    h = H[idx] + (edges-t[idx]) / (2*n[idx])
    haploid = np.r_[np.diff(edges)/np.diff(h), 2*n[-1]]
    np.savetxt(work / 'eas.demo', np.c_[edges, haploid], fmt='%.15g')
    np.savetxt(work / 'constant.demo', [[0., .00075 / (2*MU)]], fmt='%.15g')
    # Identical time states in both arms, including 0-1 kya and an old tail.
    bounds = np.unique(np.r_[np.arange(0, 2001, 40), np.geomspace(2000, 200000, 25)])
    np.savetxt(work / 'times.disc', bounds, fmt='%.15g')
    probe = np.unique(np.r_[t, edges])
    ix = np.searchsorted(t, probe, side='right')-1
    exact_h = H[ix] + (probe-t[ix])/(2*n[ix])
    iy = np.searchsorted(edges, probe, side='right')-1
    approx_h = h[iy] + (probe-edges[iy])/haploid[iy]
    prior_error = float(np.max(np.abs(np.exp(-exact_h)-np.exp(-approx_h))))
    assert prior_error < .001
    contract = dict(schema=1, regions=contracts, mu=MU, recombination_rate=REC, generation_time=GEN,
                    source_manifest_sha256=sha(SOURCE),
                    eas_source_sha256=sha(ROOT/'no_introgression_EAS_sim/resources/EAS.npz'),
                    pairs_sha256=sha(work/'pairs.tsv'), positions=POSITIONS.tolist(),
                    demography_max_absolute_cdf_error=prior_error, asmc_csfs=False,
                    selected_episode_kya=[20,40], selection=.02, origin='introgression',
                    files={s: sha(work/s) for s in ['eas.demo','constant.demo','times.disc','neutral.frq']})
    save(work / 'manifest.json', contract)
    print(f'Prepared {len(rows)} regions; demographic CDF approximation error {prior_error:.3g}', flush=True)


def quantities(a):
    from asmc.preparedecoding import prepare_decoding
    for model in ['constant', 'eas']:
        prefix = a.work / model
        receipt = a.work / f'{model}.receipt.json'
        contract = dict(demo=sha(a.work/f'{model}.demo'), times=sha(a.work/'times.disc'),
                        frequencies=sha(a.work/'neutral.frq'), mu=MU, samples=200,
                        version=importlib.metadata.version('asmc-preparedecoding'))
        if receipt.exists():
            old = json.loads(receipt.read_text())
            assert old['contract'] == contract
            assert sha(str(prefix)+'.decodingQuantities.gz') == old['sha256']
            continue
        start = time.monotonic()
        dq = prepare_decoding(demography=str(a.work/f'{model}.demo'),
                              discretization=str(a.work/'times.disc'), frequencies=str(a.work/'neutral.frq'),
                              samples=200, mutation_rate=MU)
        dq.save_decoding_quantities(str(prefix)); dq.save_intervals(str(prefix))
        save(receipt, dict(contract=contract, seconds=time.monotonic()-start,
                           sha256=sha(str(prefix)+'.decodingQuantities.gz')))
        print(f'{model}: decoding quantities ready', flush=True)


def decode_one(args):
    work, i, model = args
    from asmc.asmc import ASMC, DecodingParams
    d = work / f'region_{i:02d}'
    out, receipt = d/f'{model}.npz', d/f'{model}.json'
    contract = dict(dq=sha(work/f'{model}.decodingQuantities.gz'), pairs=sha(work/'pairs.tsv'),
                    hap=sha(d/'input.hap'), map=sha(d/'input.map'), samples=sha(d/'input.samples'),
                    version=importlib.metadata.version('asmc-asmc'), sequence=True, csfs=False)
    if receipt.exists():
        old = json.loads(receipt.read_text())
        assert old['contract'] == contract and sha(out) == old['sha256']
        return dict(region=i, model=model, cached=True)
    # Native C++ output needs OS-level redirection, not redirect_stdout.
    with (d/f'{model}.log').open('w') as log:
        oldout, olderr = os.dup(1), os.dup(2)
        try:
            os.dup2(log.fileno(), 1); os.dup2(log.fileno(), 2)
            started = time.monotonic()
            p = DecodingParams(in_file_root=str(d/'input'), dq_file=str(work/f'{model}.decodingQuantities.gz'),
                               decoding_mode_string='sequence', decoding_sequence=True,
                               using_CSFS=False, compress=False, skip_CSFS_distance=float('inf'),
                               do_posterior_sums=True)
            engine = ASMC(p)
            actual = engine.get_decoding_params()
            assert actual.decodingSequence and not actual.usingCSFS
            # ASMC 1.4 accumulates posterior sums inside the per-pair-mean path.
            engine.set_store_per_pair_posterior_mean(True)
            engine.set_store_per_pair_map(True)
            engine.set_store_sum_of_posterior(True)
            pairs = np.loadtxt(work/'pairs.tsv', dtype=int)
            engine.decode_pairs(pairs[:,0].tolist(), pairs[:,1].tolist())
            res = engine.get_ref_of_results()
            returned = np.array([[x[0], x[2]] for x in res.per_pair_indices], dtype=int)
            np.testing.assert_array_equal(returned, pairs)
            pos = np.asarray(engine.get_physical_positions(), dtype=int)-1
            np.testing.assert_array_equal(pos, np.load(d/'marker_positions.npy'))
            np.testing.assert_allclose(engine.get_genetic_positions(), pos*REC, rtol=2e-6, atol=1e-10)
            take = np.searchsorted(pos, POSITIONS)
            np.testing.assert_array_equal(pos[take], POSITIONS)
            mass = np.array(res.sum_of_posteriors[:,take], dtype=float) / len(pairs)
            assert np.isfinite(mass).all() and mass.min() >= -1e-8
            np.testing.assert_allclose(mass.sum(0), 1., atol=2e-5)
            # First 50 states exactly cover [0,1),...,[49,50) kya.
            cdfs = np.cumsum(mass, axis=0)[:50].T
            np.savez_compressed(out, cdfs=cdfs, state_mass=mass, positions=POSITIONS,
                                expected_times=np.asarray(engine.get_expected_times()))
            seconds = time.monotonic()-started
            save(receipt, dict(contract=contract, sha256=sha(out), seconds=seconds,
                               normalization_max_error=float(abs(mass.sum(0)-1).max())))
        finally:
            os.dup2(oldout, 1); os.dup2(olderr, 2); os.close(oldout); os.close(olderr)
    return dict(region=i, model=model, seconds=seconds)


def decode(a):
    manifest = json.loads((a.work/'manifest.json').read_text())
    tasks = [(a.work, i, model) for i in range(len(manifest['regions'])) for model in a.models]
    for model in a.models:
        ready = json.loads((a.work/f'{model}.receipt.json').read_text())
        assert sha(a.work/f'{model}.decodingQuantities.gz') == ready['sha256']
    if a.limit:
        tasks = tasks[:a.limit]
    with ProcessPoolExecutor(a.workers) as pool:
        for result in pool.map(decode_one, tasks):
            print(json.dumps(result), flush=True)


def report(a):
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.ticker import PercentFormatter
    out = a.out.resolve(); out.mkdir(parents=True, exist_ok=True)
    cache = out/'comparison.npz'
    if a.from_cache:
        saved = json.loads((out/'provenance.json').read_text())
        assert sha(cache) == saved['cache_sha256']
        manifest = json.loads((out/'manifest.json').read_text())
    else:
        manifest = json.loads((a.work/'manifest.json').read_text())
        regions, receipts = [], []
        for i, row in enumerate(manifest['regions']):
            d = a.work/f'region_{i:02d}'
            with np.load(d/'baseline.npz') as z:
                values = list(z['cdfs'])
            for model in ['constant', 'eas']:
                receipt = json.loads((d/f'{model}.json').read_text())
                assert receipt['sha256'] == sha(d/f'{model}.npz')
                receipts.append(dict(region=i, model=model, **receipt))
                with np.load(d/f'{model}.npz') as z:
                    values.append(z['cdfs'])
            regions.append(np.stack(values))
        cdfs = np.stack(regions, axis=1)
        assert cdfs.shape == (5, len(regions), 100, 50)
        assert np.isfinite(cdfs).all() and cdfs.min() >= -1e-7 and cdfs.max() <= 1+2e-5
        assert np.diff(cdfs, axis=-1).min() >= -1e-7
        np.savez_compressed(cache, cdfs=cdfs, positions=POSITIONS,
                            methods=np.array(METHODS), ids=np.array([r['id'] for r in manifest['regions']]),
                            selected=np.array([r['cohort']=='selected' for r in manifest['regions']]))
        shutil.copy2(a.work/'manifest.json', out/'manifest.json')
        for name in ['eas.demo','constant.demo','times.disc','pairs.tsv','constant.receipt.json','eas.receipt.json']:
            shutil.copy2(a.work/name, out/name)
        for name in ['audit.json','audit.log','quantities.log']:
            if (a.work/name).exists():
                if name.endswith('.log'):
                    (out/name).write_text('\n'.join(line.rstrip() for line in (a.work/name).read_text().splitlines())+'\n')
                else:
                    shutil.copy2(a.work/name,out/name)
        for model in ['constant','eas']:
            shutil.copy2(a.work/f'region_00/{model}.log',out/f'example_{model}.log')
        save(out/'decoding_receipts.json', receipts)
        save(out/'provenance.json', dict(cache_sha256=sha(cache), script_sha256=sha(__file__),
             packages={x: importlib.metadata.version(x) for x in ['asmc-asmc','asmc-preparedecoding','numpy','tskit']},
             workers=a.workers, RNG='No new random draws; original fixed pairs and simulation seeds reused.',
             evaluation='Same 1000 pairs and exact 10-kb coordinates; invariant anchors supply no new variation.',
             scope=f'{int(sum(r["cohort"]=="selected" for r in manifest["regions"]))} selected plus '
                   f'{int(sum(r["cohort"]=="neutral_test" for r in manifest["regions"]))} neutral regions; '
                   'descriptive inference accuracy, not a calibrated power study.'))
    with np.load(cache) as z:
        cdfs, selected = z['cdfs'], z['selected']
    rows = []
    for j, method in enumerate(METHODS):
        for i, record in enumerate(manifest['regions']):
            focal = cdfs[j,i,50]
            rows.append(dict(id=record['id'], cohort=record['cohort'], method=method,
                fraction_lt20=float(focal[19]), fraction_20to40=float(focal[39]-focal[19]),
                fraction_lt50=float(focal[49]),
                focal_cdf_mae=float(abs(focal-cdfs[0,i,50]).mean()),
                region_cdf_mae=float(abs(cdfs[j,i]-cdfs[0,i]).mean())))
    with (out/'replicate_metrics.tsv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n'); writer.writeheader(); writer.writerows(rows)
    summary = {}
    for j, method in enumerate(METHODS):
        excess = cdfs[j,selected,50].mean(0)-cdfs[j,~selected,50].mean(0)
        info = dict(excess_peak_1kyr_kya=float(np.diff(excess).argmax()+1.5))
        for yes, name in [(True,'selected'),(False,'neutral')]:
            group = [r for r in rows if r['method']==method and (r['cohort']=='selected')==yes]
            info[name] = {k:float(np.mean([r[k] for r in group])) for k in
                         ['fraction_lt20','fraction_20to40','fraction_lt50','focal_cdf_mae','region_cdf_mae']}
        summary[method] = info
    save(out/'summary.json', summary)
    matplotlib.rcParams.update({'font.size':14, 'axes.labelsize':16, 'axes.titlesize':17,
                                'pdf.fonttype':42,'ps.fonttype':42, 'axes.spines.top':False,
                                'axes.spines.right':False})
    labels = ['Truth', 'Gamma median', 'Gamma mass', 'ASMC constant', 'ASMC EAS']
    colors = ['#222222','#CC79A7','#D55E00','#0072B2','#009E73']
    years = np.arange(1,51)
    fig, ax = plt.subplots(2,2,figsize=(11,8.5))
    for row, mask in enumerate([selected,~selected]):
        for j in range(5):
            mean = cdfs[j,mask,50].mean(0)
            style = ':' if j==1 else '-'
            ax[row,0].plot(years,mean,color=colors[j],ls=style,lw=2,label=labels[j])
            ax[row,1].plot(years[:-1]+.5,np.diff(mean),color=colors[j],ls=style,lw=2)
        for col in range(2):
            ax[row,col].axvspan(20,40,color='gray',alpha=.12,zorder=-1)
            ax[row,col].set_xlim(1,50); ax[row,col].yaxis.set_major_formatter(PercentFormatter(1))
            ax[row,col].set_xlabel('TMRCA (kya)')
        ax[row,0].set_ylabel('Cumulative pair fraction')
        ax[row,1].set_ylabel('Pair mass per 1-kyr bin')
        title = 'Selected: 40–20 kya, s = 0.02' if row==0 else 'Neutral introgression'
        ax[row,0].set_title(title); ax[row,1].set_title('Coalescence-time distribution')
    fig.legend(*ax[0,0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.5,.985),
               ncol=5,frameon=False,fontsize=12,handlelength=1.5,columnspacing=1)
    fig.subplots_adjust(top=.88,bottom=.09,left=.11,right=.98,wspace=.36,hspace=.52)
    for ext in ['png','pdf']: fig.savefig(out/f'01_focal_tmrca.{ext}',dpi=180)
    plt.close(fig)
    fig, ax = plt.subplots(1,2,figsize=(11,6))
    for j in [0,2,3,4]:
        delta=cdfs[j,selected,50].mean(0)-cdfs[j,~selected,50].mean(0)
        ax[0].plot(years[:-1]+.5,np.diff(delta),lw=2,color=colors[j],label=labels[j])
        # Five-kyr masses use boundaries 0,5,...,50, preserving full normalization.
        edges = np.r_[0.,delta[4::5]]
        ax[1].plot(np.arange(2.5,50,5),np.diff(edges),lw=2,marker='o',color=colors[j])
    for k in range(2):
        ax[k].axvspan(20,40,color='gray',alpha=.12,zorder=-1)
        ax[k].axhline(0,color='gray',lw=.7); ax[k].set_xlim(0,50)
        ax[k].set_xlabel('TMRCA (kya)'); ax[k].yaxis.set_major_formatter(PercentFormatter(1))
    ax[0].set_ylabel('Selected − neutral pair mass'); ax[0].set_title('1-kyr bins')
    ax[1].set_title('5-kyr bins')
    fig.legend(*ax[0].get_legend_handles_labels(),loc='upper center',ncol=4,frameon=False,fontsize=13)
    fig.subplots_adjust(top=.82,bottom=.15,left=.1,right=.98,wspace=.3)
    for ext in ['png','pdf']: fig.savefig(out/f'02_selected_excess.{ext}',dpi=180)
    plt.close(fig)
    print(json.dumps(summary,indent=2),flush=True)


def audit(a):
    """Check ASMC units and whether invariant output anchors perturb inference."""
    import tskit
    from asmc.asmc import ASMC, DecodingParams
    manifest = json.loads((a.work/'manifest.json').read_text())
    directory = a.work/'region_00'
    check = a.work/'anchor_audit'; check.mkdir(exist_ok=True)
    ts = tskit.load(local(manifest['regions'][0]['source'])/'observed.trees')
    sites = {int(v.site.position) for v in ts.variants() if len(v.alleles)==2}
    for suffix in ['hap','map']:
        lines = (directory/f'input.{suffix}').read_text().splitlines()
        lines = [line for line in lines if int(line.split()[2 if suffix=='hap' else 3])-1 in sites]
        (check/f'input.{suffix}').write_text('\n'.join(lines)+'\n')
    shutil.copy2(directory/'input.samples',check/'input.samples')
    pairs = np.loadtxt(a.work/'pairs.tsv',dtype=int)
    results = {}
    for model in ['constant','eas']:
        comparisons = []
        for base in [directory,check]:
            p = DecodingParams(in_file_root=str(base/'input'),dq_file=str(a.work/f'{model}.decodingQuantities.gz'),
                 decoding_mode_string='sequence',decoding_sequence=True,using_CSFS=False,
                 skip_CSFS_distance=float('inf'),do_posterior_sums=True)
            engine = ASMC(p)
            assert engine.get_decoding_params().decodingSequence
            assert not engine.get_decoding_params().usingCSFS
            engine.set_store_per_pair_posterior_mean(True); engine.set_store_per_pair_map(True)
            engine.set_store_sum_of_posterior(True)
            engine.decode_pairs(pairs[:,0].tolist(),pairs[:,1].tolist())
            pos = np.asarray(engine.get_physical_positions())-1
            np.testing.assert_allclose(engine.get_genetic_positions(),pos*REC,rtol=2e-6,atol=1e-10)
            vals = engine.get_ref_of_results()
            np.testing.assert_array_equal([[v[0],v[2]] for v in vals.per_pair_indices],pairs)
            mass = np.asarray(vals.sum_of_posteriors).copy()/len(pairs)
            np.testing.assert_allclose(mass.sum(0),1,atol=2e-5)
            take = np.isin(pos,list(sites))
            comparisons.append((pos[take], np.cumsum(mass,axis=0)[:50,take]))
        np.testing.assert_array_equal(comparisons[0][0],comparisons[1][0])
        error = abs(comparisons[0][1]-comparisons[1][1])
        center=(comparisons[0][0]>=400000)&(comparisons[0][0]<=600000)
        nearest=int(np.argmin(abs(comparisons[0][0]-500000)))
        results[model] = dict(anchor_cdf_mae=float(error.mean()),anchor_cdf_max=float(error.max()),
                              central_anchor_cdf_mae=float(error[:,center].mean()),
                              central_anchor_cdf_max=float(error[:,center].max()),
                              nearest_focal_marker=int(comparisons[0][0][nearest]),
                              nearest_focal_cdf_max=float(error[:,nearest].max()),
                              max_error_position=int(comparisons[0][0][np.unravel_index(error.argmax(),error.shape)[1]]))
        print(json.dumps(results[model]),flush=True)
        # Report the complete sensitivity, including boundary outliers. This is
        # not an invariance test: ASMC approximates transitions across gaps.
        intervals = np.loadtxt(a.work/f'{model}.intervalsInfo')
        np.testing.assert_allclose(intervals[:50,0],np.arange(0,2000,40))
        np.testing.assert_allclose(intervals[:50,2],np.arange(40,2001,40))
        if model=='constant':
            scale=.00075/(2*MU)
            expected = 40/np.expm1(40/scale)
            expected = scale-expected+np.arange(0,2000,40)
            np.testing.assert_allclose(intervals[:50,1],expected,atol=1e-6)
    save(a.work/'audit.json',results)
    print(json.dumps(results,indent=2),flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('phase', choices=['prepare','quantities','decode','report','audit'])
    p.add_argument('--work', type=Path, default=ROOT.parent/'sim/asmc_gamma_demo_20261106')
    p.add_argument('--replicates', type=int, default=10)
    p.add_argument('--workers', type=int, default=4)
    p.add_argument('--limit', type=int, default=0)
    p.add_argument('--models', nargs='+', choices=['constant','eas'], default=['constant','eas'])
    p.add_argument('--out', type=Path, default=ROOT/'docs/results/asmc_gamma_demo_20261106')
    p.add_argument('--from-cache', action='store_true')
    a = p.parse_args(); a.work = a.work.resolve()
    assert 1 <= a.workers <= 24 and 1 <= a.replicates <= 50
    globals()[a.phase](a)


if __name__ == '__main__':
    main()
