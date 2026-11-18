"""Simulate intermediate selection strengths and score them against the frozen EAS nulls."""
import argparse
import csv
import json
import os
import shutil
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import report_cumulative_grid as previous
import numpy as np

STRENGTHS = (.002, .003, .004)
RAW = Path.home() / 'sim/cumulative_grid_20261118'
OUT = previous.ROOT / 'docs/results/cumulative_grid_20261118'


def simulate(raw, workers, slim=None, decoder=None):
    import run_cumulative_grid as grid

    manifest_source = previous.OUT / 'simulation_manifest.json'
    baseline = previous.read(manifest_source)
    configs = baseline['configs']
    for cfg in configs.values():
        if slim:
            cfg['slim'] = str(slim.resolve(strict=True))
        if decoder:
            cfg['decoder'] = str(decoder.resolve(strict=True))
        cfg['seed_namespace'] = 'cumulative-eas-intermediate-20261118:'
        for path, key in ((cfg['slim'], 'slim_sha256'), (cfg['decoder'], 'decoder_sha256'),
                          (previous.ROOT / 'no_introgression_EAS_sim/resources/EAS.npz', 'phlash_sha256'),
                          (grid.arrays.ancestor.__file__, 'ancestry_source_sha256'),
                          (grid.base.__file__, 'score_source_sha256')):
            assert previous.sha(Path(path)) == cfg[key], path
    assert previous.sha(Path(grid.arrays.__file__)) == baseline['script_hashes']['run_selection_arrays.py']
    strata = [
        dict(key=f"{origin}_{regime}_s{str(s).replace('.', 'p')}", origin=origin, regime=regime,
             start=start, end=end, selection=s, role='selected')
        for s in STRENGTHS for regime, (start, end) in grid.WINDOWS.items() for origin in previous.ORIGINS
    ]
    raw.mkdir(parents=True, exist_ok=True)
    contract = dict(configs=configs, strata=strata, target=100, input_hashes=baseline['input_hashes'],
                    parent_manifest_sha256=previous.sha(manifest_source),
                    scripts={Path(p).name: previous.sha(Path(p)) for p in
                             (grid.__file__, grid.cached.__file__, grid.arrays.__file__, grid.arrays.ancestor.__file__, grid.base.__file__)})
    manifest = raw / 'manifest.json'
    if manifest.exists():
        assert previous.read(manifest) == contract
    else:
        grid.base.save(manifest, contract)
    cfg = configs['50k']
    possible_pairs = np.column_stack(np.triu_indices(2 * cfg['diploids'], 1))
    sample = np.random.default_rng(cfg['pair_seed']).choice(len(possible_pairs), cfg['pairs'], replace=False)
    inputs = {'pairs.tsv': possible_pairs[np.sort(sample)], 'positions.txt': np.arange(0, cfg['length'], cfg['stride'])}
    for name, digest in contract['input_hashes'].items():
        if not (raw / name).exists():
            np.savetxt(raw / name, inputs[name], fmt='%d', delimiter='\t')
        assert previous.sha(raw / name) == digest
    lib = Path(configs['50k']['slim']).parent.parent / 'lib'
    os.environ['LD_LIBRARY_PATH'] = str(lib) + ':' + os.environ.get('LD_LIBRARY_PATH', '')
    print(json.dumps(dict(started=time.time(), raw=str(raw), workers=workers, manifest=str(manifest), slots=1800)), flush=True)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(grid.simulate_slot, (raw, configs[s['regime']], s, slot))
                   for slot in range(contract['target']) for s in strata]
        for count, future in enumerate(as_completed(futures), 1):
            future.result()
            if count % 50 == 0:
                print(f'{count}/1800 finished at {time.time()}', flush=True)


def prepare(raw):
    manifest = raw / 'manifest.json'
    contract = previous.read(manifest)
    signature = previous.sha(manifest)
    source = previous.OUT / 'mass.npz'
    assert previous.sha(source) == previous.read(previous.OUT / 'inputs.json')['cache_sha256']
    with np.load(source) as data:
        curves = {k: data[k] for k in data.files if '_neutral' in k or '_test' in k}
    sources = []
    for s in contract['strata']:
        group, af, ids = [], [], []
        for slot in range(contract['target']):
            rec = previous.read(raw / 'slots' / s['key'] / f'{slot:03d}.json')
            directory = Path(rec['source_directory'])
            assert rec['config_sha256'] == rec['manifest_sha256'] == signature
            assert (rec['origin'], rec['selection'], rec['selection_start_years'], rec['selection_end_years']) == (
                s['origin'], s['selection'], s['start'], s['end'])
            assert rec['accepted'] and rec['stop_years_ago'] == 0 and 0 < rec['final_frequency'] < 1
            for name, field in (('scores.npz', 'scores_sha256'), ('trajectory.tsv', 'trajectory_sha256')):
                assert previous.sha(directory / name) == rec[field]
            with np.load(directory / 'scores.npz') as score:
                np.testing.assert_array_equal(score['cutoffs'], np.arange(1000, 50001, 1000))
                np.testing.assert_array_equal(score['positions'], np.arange(0, 1000000, 10000))
                np.testing.assert_array_equal(score['valid_pairs'], 1000)
                group.append(score['cdfs'][2, 50])
            ids.append(str(directory))
            af.append(rec['final_frequency'])
            sources.append(dict(group=s['key'], slot=slot, id=str(directory), seed=rec['seed'],
                                onset_sha256=rec['onset_sha256'], scores_sha256=rec['scores_sha256']))
        key = s['key']
        curves[key], curves[key+'_ids'], curves[key+'_af'] = np.array(group), np.array(ids), np.array(af)
    assert len(sources) == len({r['onset_sha256'] for r in sources}) == 1800
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / '.gitattributes').write_text('* -text\n')
    shutil.copyfile(manifest, OUT / 'simulation_manifest.json')
    np.savez_compressed(OUT / 'mass.npz', **curves)
    previous.save(OUT / 'inputs.json', dict(raw_manifest_sha256=signature, neutral_parent_sha256=previous.sha(source),
                  sources=sources, cache_sha256=previous.sha(OUT / 'mass.npz')))


def report():
    previous.report(OUT, STRENGTHS)
    def neutral_counts(folder):
        with (folder / 'rates.tsv').open() as stream:
            return {(r['origin'], r['window'], r['alpha']): (int(r['detections']), int(r['n']))
                    for r in csv.DictReader(stream, delimiter='\t') if r['endpoint'] == 'FPR'}
    reference = neutral_counts(previous.OUT)
    assert len(reference) == 12 and neutral_counts(OUT) == reference


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=('simulate', 'prepare', 'report'))
    parser.add_argument('--raw', type=Path, default=RAW)
    parser.add_argument('--workers', type=int, choices=range(1, 21), default=4)
    parser.add_argument('--slim', type=Path, help='SLiM binary matching the baseline checksum')
    parser.add_argument('--decoder', type=Path, help='Gamma-SMC binary matching the baseline checksum')
    args = parser.parse_args()
    if args.phase == 'simulate':
        simulate(args.raw.resolve(), args.workers, args.slim, args.decoder)
    elif args.phase == 'prepare':
        prepare(args.raw.resolve())
    else:
        report()
