#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = ["pysam==0.24.1", "google-crc32c==1.8.0"]
# ///
"""Rebuild GTs in completed same-position p70 VCFs using the saved clusters."""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from contextlib import ExitStack
import fcntl
from functools import lru_cache
import gzip
import json
import logging
import os
from pathlib import Path
import re
import sys

import pysam

import run_chr1_sv_merge as pilot
from phased_sv_union import union_calls


def checked(path, manifest):
    if pilot.common.digest(path) != manifest['artifacts'][path.name]:
        raise ValueError(f'Checksum mismatch: {path}')
    return path


@lru_cache(maxsize=8192)
def decode(gt, ps):
    phased = b'|' in gt
    alleles = tuple(None if a == b'.' else int(a) for a in re.split(rb'[|/]', gt))
    return dict(gt=alleles, phased=phased, ps=None if ps == b'.' else ps)


def combine(calls):
    result = union_calls(calls)
    sep = b'|' if result['phased'] else b'/'
    gt = sep.join(b'.' if a is None else str(a).encode() for a in result['gt'])
    return gt, result['ps'] or b'.', result['gt'].count(1), sum(a is not None for a in result['gt']), result['status']


def read_cluster(members, inputs, offsets, samples):
    calls = []
    for m in members:
        handle = inputs[m['kind']]
        handle.seek(offsets[m['working_id']])
        fields = handle.readline().split(b'\t', 9)
        if fields[2].decode() != m['working_id']:
            raise ValueError('Invalid member offset')
        fmt = fields[8].split(b':')
        gt_idx = fmt.index(b'GT')
        ps_idx = fmt.index(b'PS') if b'PS' in fmt else None
        column = []
        for value in fields[9].split(b'\t'):
            parts = value.split(b':')
            gt = parts[gt_idx] if gt_idx < len(parts) else b'.'
            ps = parts[ps_idx] if ps_idx is not None and ps_idx < len(parts) else b'.'
            column.append(decode(gt, ps))
        if len(column) != samples:
            raise ValueError('Wrong number of prepared samples')
        calls.append(column)
    return calls


def repaired_line(fields, calls, use_ps, stats):
    fmt, raw = fields[8].split(b'\t', 1)
    gt_idx = fmt.split(b':').index(b'GT')
    old_samples = raw.split(b'\t')
    if len(old_samples) != len(calls[0]):
        raise ValueError('Wrong number of final samples')
    genotypes, ac, an = [], 0, 0
    # Most samples share one of a few GT/PS patterns. The cache is per cluster.
    cache = {}
    for idx, members in enumerate(zip(*calls, strict=True)):
        key = tuple((v['gt'], v['phased'], v['ps']) for v in members)
        if key not in cache:
            value = combine(members)
            if len(cache) < 1024:
                cache[key] = value
        else:
            value = cache[key]
        gt, ps, dosage, called, status = value
        old = old_samples[idx].split(b':')[gt_idx]
        stats['changed_GT_or_phase'] += old != gt
        stats['known_ALT_copies_added_net'] += dosage - old.replace(b'|', b'/').split(b'/').count(b'1')
        stats[status] += 1
        genotypes.append(gt + b':' + ps if use_ps else gt)
        ac += dosage
        an += called
    # Restrict repaired records to structural annotations and cluster provenance.
    keep = {b'SVTYPE', b'SVLEN', b'END', b'CollapseId', b'NumCollapsed',
            b'CollapseStart', b'CollapseEnd', b'CollapseSize'}
    info = [v for v in fields[7].split(b';') if v.split(b'=', 1)[0] in keep]
    af = format(ac / an, '.10g') if an else '.'
    info += [f'AC={ac}'.encode(), f'AN={an}'.encode(), f'AF={af}'.encode()]
    return b'\t'.join([*fields[:7], b';'.join(info), b'GT:PS' if use_ps else b'GT', *genotypes])


def repair(task):
    chrom, batch, out = task
    base = batch/chrom
    dest = out/chrom
    dest.mkdir(parents=True, exist_ok=True)
    with ExitStack() as stack:
        lock = stack.enter_context((dest/'.lock').open('a'))
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        logging.basicConfig(filename=dest/'repair.log', level=logging.INFO, force=True,
                            format='%(asctime)s %(levelname)s %(message)s')
        final_dir = base/'final'
        manifest = json.loads((final_dir/'done.json').read_text())
        policy = manifest['key']['policy']
        if (policy['gt'], policy['pctseq'], policy['pctsize'], policy['grouping']) != ('all', .7, .7, 'original_CHROM_POS_and_type'):
            raise ValueError(f'{chrom}: unexpected source merge policy')
        prepared = manifest['key']['prepared']
        if prepared['key']['settings']['chrom'] != chrom:
            raise ValueError('Chromosome manifest mismatch')
        key = dict(repair_schema=1, original_artifacts=manifest['artifacts'],
                   prepared_artifacts=prepared['artifacts'], policy=policy)
        if pilot.valid_manifest(dest, key):
            return json.loads((dest/'summary.json').read_text())
        (dest/'done.json').unlink(missing_ok=True)
        logging.info('PARAMETERS chrom=%s batch=%s out=%s seed=%s python=%s pysam=%s policy=%s',
                     chrom,batch,out,pilot.SEED,sys.version,pysam.__version__,policy)
        print(f'Starting {chrom}', flush=True)
        source = checked(final_dir/f'{chrom}.popped.samepos.p70.vcf.gz', manifest)
        mapping = checked(final_dir/'old_to_new.jsonl.gz', manifest)
        previous = json.loads(checked(final_dir/'summary.json', manifest).read_text())
        samples = checked(base/'prepare/samples.txt', prepared).read_bytes().splitlines()
        groups = defaultdict(list)
        needed = {}
        with gzip.open(mapping, 'rt') as rows:
            for row in map(json.loads, rows):
                if row['members'] <= 1:
                    continue
                m = row['old']
                k = tuple(str(x).encode() for x in (chrom, row['new_pos'], row['new_id'], row['new_ref'], row['new_alt']))
                groups[k].append(m)
                if m['working_id'] in needed:
                    raise ValueError('Duplicate working ID in saved mapping')
                needed[m['working_id']] = m
        if len(groups) != previous.get('merged_representatives', 0) or len(needed)-len(groups) != previous.get('removed_records', 0):
            raise ValueError('Saved membership count mismatch')
        offsets, inputs = {}, {}
        for kind in ('INS','DEL'):
            path = checked(base/'prepare'/f'{kind}.vcf.gz', prepared)
            handle = stack.enter_context(pysam.BGZFile(str(path), 'r'))
            inputs[kind] = handle
            while True:
                offset = handle.tell()
                line = handle.readline()
                if not line:
                    break
                if line.startswith(b'#'):
                    if line.startswith(b'#CHROM') and line.split(b'\t')[9:] != samples:
                        raise ValueError('Prepared sample order mismatch')
                    continue
                site = line.split(b'\t', 5)[:5]
                ident = site[2].decode()
                if ident in needed:
                    m = needed[ident]
                    if (kind, chrom, int(site[1]), site[3].decode(), site[4].decode()) != (
                        m['kind'], site[0].decode(), m['trimmed_pos'], m['trimmed_ref'], m['trimmed_alt']
                    ) or ident in offsets:
                        raise ValueError('Prepared allele identity mismatch')
                    offsets[ident] = offset
        if offsets.keys() != needed.keys():
            raise ValueError('Missing original member genotypes')
        stats = Counter()
        target = dest/f'{chrom}.popped.samepos.p70.gt_union.vcf.gz'
        partial = dest/'repair.partial.vcf.gz'
        with pysam.BGZFile(str(source), 'r') as src, pysam.BGZFile(str(partial),'w') as dst:
            use_ps = False
            declared = set()
            for line in src:
                if line.startswith(b'#'):
                    for tag in (b'AC', b'AN', b'AF'):
                        if line.startswith(b'##INFO=<ID=' + tag + b','):
                            declared.add(tag)
                    use_ps |= line.startswith(b'##FORMAT=<ID=PS,')
                    if line.startswith(b'#CHROM'):
                        if line.split(b'\t')[9:] != samples:
                            raise ValueError('Final sample order mismatch')
                        for tag, number, dtype in [(b'AC',b'A',b'Integer'),(b'AN',b'1',b'Integer'),(b'AF',b'A',b'Float')]:
                            if tag not in declared:
                                dst.write(b'##INFO=<ID='+tag+b',Number='+number+b',Type='+dtype+b',Description="Derived from GT on repaired records">\n')
                else:
                    stats['records'] += 1
                    site = tuple(line.split(b'\t', 5)[:5])
                    if site in groups:
                        members = groups.pop(site)
                        calls = read_cluster(members, inputs, offsets, len(samples))
                        line = repaired_line(line.split(b'\t',8), calls, use_ps, stats)
                        stats['repaired_'+members[0]['kind']] += 1
                        del calls
                    if stats['records'] % 100000 == 0:
                        logging.info('PROGRESS %s', dict(stats))
                dst.write(line + b'\n')
        if groups or stats['records'] != previous['after_records']:
            raise ValueError('Missing representatives or final record count mismatch')
        pysam.tabix_index(str(partial), preset='vcf', csi=True, force=True)
        partial_idx = Path(str(partial)+'.csi')
        with pysam.VariantFile(str(partial), index_filename=str(partial_idx)) as check:
            if list(check.header.samples) != [s.decode() for s in samples] or (stats['records'] and next(check.fetch(chrom), None) is None):
                raise ValueError('Indexed output validation failed')
        partial.replace(target)
        index = Path(str(target)+'.csi')
        partial_idx.replace(index)
        report = dict(chrom=chrom, samples=len(samples), **stats)
        summary = dest/'summary.json'
        pilot.common.json_write(summary,report)
        pilot.manifest(dest,key,[target,index,summary])
        return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch-root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--jobs', type=int, choices=range(1,5), default=4)
    autosomes = [f'chr{i}' for i in range(1,22)]
    parser.add_argument('--chromosomes', nargs='+', choices=autosomes, default=autosomes)
    args = parser.parse_args()
    if len(set(args.chromosomes)) != len(args.chromosomes):
        parser.error('Duplicate chromosomes')
    os.umask(0o077)
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        for report in pool.map(repair, [(c,args.batch_root.resolve(),args.out.resolve()) for c in args.chromosomes]):
            print(json.dumps(report), flush=True)
