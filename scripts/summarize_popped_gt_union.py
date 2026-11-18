#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = ["pysam==0.24.1", "google-crc32c==1.8.0"]
# ///
"""Summarize completed popped panel GT repairs using the original and final allele strings."""
import argparse
import csv
import fcntl
import json
import os
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pysam
from inspect_sv_representation import describe
from phase_singletons import digest, identity, json_write


def read(path):
    return json.loads(path.read_text())


def scan(path, chrom, expected_records):
    counts = Counter()
    with pysam.VariantFile(str(path), drop_samples=True) as source:
        for record in source.fetch(chrom):
            counts['records'] += 1
            if len(record.alts or ()) != 1:
                counts['excluded_non_biallelic_records'] += 1
                continue
            ref, alt = record.ref, record.alts[0]
            if max(len(ref), len(alt)) < 50:
                continue
            details, _ = describe(ref, alt)
            if details is None:
                continue
            if details['complex']:
                if details['net_bp'] >= 50:
                    counts['complex_net50'] += 1
                elif max(details['residual_ref_bp'], details['residual_alt_bp']) >= 50:
                    counts['complex_span50_net_lt50'] += 1
            elif details['net_bp'] >= 50 and not details['ambiguous']:
                counts[details['type']] += 1
    if counts['records'] != expected_records:
        raise ValueError(f'{chrom}: record count differs from completed run at {path}')
    return counts


def summarize(task):
    chrom, original_root, repaired_root, out = task
    merged = original_root/chrom/'final'
    repaired = repaired_root/chrom
    old_manifest, new_manifest = (read(d/'done.json') for d in (merged, repaired))
    if old_manifest['artifacts'] != new_manifest['key']['original_artifacts']:
        raise ValueError(f'{chrom}: repair refers to a different original merge')
    summaries = []
    for path, manifest in ((merged, old_manifest), (repaired, new_manifest)):
        summary = path/'summary.json'
        if digest(summary) != manifest['artifacts'][summary.name]:
            raise ValueError(f'{chrom}: incorrect summary checksum at {path}')
        summaries.append(read(summary))
    merge, fix = summaries
    if merge['chrom'] != chrom or fix['chrom'] != chrom or merge['after_records'] != fix['records']:
        raise ValueError(f'{chrom}: inconsistent chromosome summaries')
    source_files = old_manifest['key']['source']
    for snapshot in source_files:
        if identity(snapshot['path']) != snapshot:
            raise ValueError(f'{chrom}: original popped input or index changed')
    source = Path(source_files[0]['path'])
    final = repaired/f'{chrom}.popped.samepos.p70.gt_union.vcf.gz'
    index = Path(str(final)+'.csi')
    key = dict(schema=1, source=source_files, final=[identity(final),identity(index)],
               committed=new_manifest['artifacts'], original_summary=old_manifest['artifacts']['summary.json'])
    cache = out/f'{chrom}.json'
    if cache.exists():
        saved = read(cache)
        if saved['key'] == key:
            return saved['row']
    print(f'Scanning {chrom}', flush=True)
    before = scan(source, chrom, merge['before_records'])
    after = scan(final, chrom, fix['records'])
    row = dict(chrom=chrom)
    for kind in ('INS', 'DEL'):
        start, end = before[kind], after[kind]
        if start != merge['eligible_sv_before'].get(kind,0) or start-end != merge.get('removed_'+kind,0):
            raise ValueError(f'{chrom}: independently counted simple {kind} disagrees with saved merge')
        row.update({f'{kind}_before':start, f'{kind}_after':end, f'{kind}_removed':start-end})
    for group in ('complex_net50','complex_span50_net_lt50'):
        row[group+'_before'] = before[group]
        row[group+'_after'] = after[group]
    row['GT_changed'] = fix.get('changed_GT_or_phase', 0)
    json_write(cache, dict(key=key, row=row))
    return row


def write_table(out, rows):
    total = {'chrom':'ALL', **{k:sum(r[k] for r in rows) for k in rows[0] if k!='chrom'}}
    result = []
    for row in [*rows, total]:
        line = dict(chrom=row['chrom'])
        for kind in ('INS','DEL'):
            start = row[kind+'_before']
            for k in ('before','after','removed'):
                line[kind+'_'+k] = row[kind+'_'+k]
            line[kind+'_pct_removed'] = f"{100*row[kind+'_removed']/start:.2f}" if start else 'NA'
        line.update({k:v for k,v in row.items() if k.startswith('complex') or k=='GT_changed'})
        result.append(line)
    path = out/'chr1_21_summary.tsv'
    tmp = path.with_suffix('.tmp')
    with tmp.open('w', newline='') as f:
        writer = csv.DictWriter(f,fieldnames=list(result[0]),delimiter='\t')
        writer.writeheader()
        writer.writerows(result)
    tmp.replace(path)
    print(path.read_text(),end='')
    print(f'Saved {path}')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--original',type=Path,default=Path.home()/'popped_full_same_original_pos_p70')
    p.add_argument('--repaired',type=Path,default=Path.home()/'popped_full_same_original_pos_p70_gt_union')
    p.add_argument('--jobs',type=int,choices=range(1,5),default=4)
    args = p.parse_args()
    os.umask(0o077)
    chroms = [f'chr{i}' for i in range(1,22)]
    missing = [c for c in chroms if not (args.repaired/c/'done.json').is_file()]
    if missing:
        p.exit(1, 'Waiting for completed repairs: '+', '.join(missing)+'\n')
    out = args.repaired/'sv_gt_summary'
    out.mkdir(exist_ok=True)
    with (out/'.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX | fcntl.LOCK_NB)
        with ProcessPoolExecutor(max_workers=args.jobs) as pool:
            rows = list(pool.map(summarize, [(c,args.original,args.repaired,out) for c in chroms]))
        write_table(out,rows)
