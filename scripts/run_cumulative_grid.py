"""Run missing cumulative-mass cohorts, with independent de novo neutral test founders."""
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import argparse
import csv
import json
import shutil
import time
import run_denovo_cached_founders as cached

arrays = cached.arrays
base = arrays.base
ROOT = Path(__file__).resolve().parents[1]
BASELINE = ROOT/'docs/results/denovo_comparison_20261111'
ARCHIVE = ROOT/'docs/results/selection_arrays_eas_20261105_4h'
RAW = ROOT.parent/'sim/cumulative_grid_20261117'
WINDOWS = {'10k': (10000, 0), '50k': (50000, 0), '50_30k': (50000, 30000)}


def archive_rows():
    selected = {}
    with (ARCHIVE/'completed_units.tsv').open() as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            if row['cohort'] != 'selected' or float(row['selection']) != .001 or row['regime'] not in ('10k','50k') or row['included_prefix'] != 'True':
                continue
            selected.setdefault((row['regime'], row['origin']), []).append(row)
    return {key:sorted(rows, key=lambda r:int(r['replicate']))[:100] for key,rows in selected.items()}


def local_path(value):
    return Path(value.replace('D:', '/mnt/d').replace('\\','/'))


def simulate_slot(item):
    raw, cfg, spec, slot = item
    receipt = raw/'slots'/spec['key']/f'{slot:03d}.json'
    if receipt.exists():
        entry = json.loads(receipt.read_text())
        assert entry['manifest_sha256'] == base.digest(raw/'manifest.json')
        arrays.ancestor.verify(Path(entry['source_directory']),entry)
        return
    # Each slot has a deterministic sequence of forward histories, independent of job order.
    if spec['origin'] == 'de_novo':
        result = cached.worker(raw, cfg, spec, slot, 100000, float('inf'))
        record = result['records'][-1]
    else:
        for attempt in range(100000):
            unit = f"{spec['key']}/slot_{slot:03d}/attempt_{attempt:05d}"
            task = dict(id=unit, unit=unit, regime=spec['regime'], origin=spec['origin'],
                        cohort='selected', selection=spec['selection'], replicate=slot*100000+attempt)
            result = arrays.work((raw, cfg, task, float('inf')))
            if result['state'] == 'accepted':
                directory = raw/unit
                record = json.loads((directory/'complete.json').read_text())
                record['source_directory'] = str(directory)
                break
        else:
            raise RuntimeError('Introgression slot exhausted')
    assert record['accepted'] and 0 < record['final_frequency'] < 1
    receipt.parent.mkdir(parents=True, exist_ok=True)
    arrays.save_ledger(receipt,dict(record, manifest_sha256=base.digest(raw/'manifest.json')))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers',type=int,default=4)
    args = parser.parse_args()
    assert 1 <= args.workers <= 20
    previous = ROOT.parent/'sim/denovo_comparison_20261111'
    configs = json.loads((previous/'manifest.json').read_text())['configs']
    RAW.mkdir(parents=True,exist_ok=True)
    old = archive_rows()
    strata = []
    for regime,(start,end) in WINDOWS.items():
        for origin in ('de_novo','introgression'):
            strata.append(dict(key=f'weak_{origin}_{regime}',origin=origin,regime=regime,
                start=start,end=end,selection=.001,role='selected',reused=len(old.get((regime,origin),[]))))
        if end == 0:
            strata.append(dict(key=f'heldout_de_novo_{regime}',origin='de_novo',regime=regime,
                start=start,end=end,selection=0.,role='neutral_test',reused=0))
    for config in configs.values():
        config['seed_namespace']='cumulative-eas-20261117:'
        for field, path in [('slim_sha256',config['slim']),('decoder_sha256',config['decoder']),('phlash_sha256', ROOT/'no_introgression_EAS_sim/resources/EAS.npz')]:
            assert base.digest(path)==config[field]
    provenance = dict(configs=configs,strata=strata,target=100,archive_table_sha256=base.digest(ARCHIVE/'completed_units.tsv'),
        script_hashes={Path(p).name:base.digest(p) for p in [__file__, cached.__file__,arrays.__file__,arrays.ancestor.__file__,base.__file__]},
        input_hashes={name:base.digest(previous/name) for name in ('pairs.tsv','positions.txt')})
    manifest = RAW/'manifest.json'
    if manifest.exists(): assert json.loads(manifest.read_text())==provenance
    else: base.save(manifest,provenance)
    for name, digest in provenance['input_hashes'].items():
        if not (RAW/name).exists(): shutil.copy2(previous/name, RAW/name)
        assert base.digest(RAW/name)==digest
    os.environ['LD_LIBRARY_PATH']=str(Path(configs['50k']['slim']).parent.parent/'lib')+':'+os.environ.get('LD_LIBRARY_PATH','')
    tasks=[(RAW,configs[s['regime']],s,i) for i in range(100) for s in strata if i>=s['reused']]
    print(json.dumps(dict(workers=args.workers,slots=len(tasks),start=time.time(),strata=strata)),flush=True)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for i,_ in enumerate(pool.map(simulate_slot,tasks),1):
            if i%10==0: print(f'{i}/{len(tasks)} slots',flush=True)

if __name__=='__main__':
    main()
