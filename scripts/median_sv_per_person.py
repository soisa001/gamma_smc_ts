# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy>=2,<3", "pysam==0.24.1", "google-crc32c==1.8.0"]
# ///
"""Count observed SV allele records per person across the unified panel."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re

import numpy as np
import pysam

from build_unified_sv_panel import info_dict, profile, write_tsv
from phase_singletons import digest, identity


def genotype_counts(fmt, tail, n):
    if fmt.split(b':', 1)[0] != b'GT':
        raise ValueError('Expected GT first on unified SV record')
    raw = re.sub(rb':[^\t]*', b'', tail) if b':' in fmt else tail
    data = np.frombuffer(raw, dtype=np.uint8)
    # Fast path for uniform diploid GT. Validate field boundaries before reshaping.
    if len(raw) == 4*n-1 and np.all(data[3::4] == 9) and np.all((data[1::4] == 47) | (data[1::4] == 124)):
        left, right = data[::4], data[2::4]
        copies = (left == 49).astype(np.uint8) + (right == 49)
        called = (left != 46).astype(np.uint8) + (right != 46)
        return copies, called
    starts = np.r_[0, np.flatnonzero(data == 9)+1]
    if len(starts) != n or not len(data) or starts[-1] >= len(data):
        raise ValueError('GT column count mismatch')
    return (np.add.reduceat(data == 49, starts, dtype=np.int64),
            np.add.reduceat((data == 48) | (data == 49), starts, dtype=np.int64))


def sample_header(path):
    with pysam.VariantFile(path) as vcf:
        names = list(vcf.header.samples)
        if not names or len(names) != len(set(names)):
            raise ValueError('Invalid sample list')
        return names


def count_chr(task):
    root, out, chrom, sample_hash = task
    folder = root/chrom
    path = folder/f'{chrom}.unified.vcf.gz'
    names = sample_header(path)
    if hashlib.sha256('\n'.join(names).encode()).hexdigest() != sample_hash:
        raise ValueError(f'{chrom}: sample identity/order differs')
    stats = folder/'stats.json'
    artifacts = json.loads((folder/'done.json').read_text())['artifacts']
    if digest(stats) != artifacts[stats.name]:
        raise ValueError(f'{chrom}: panel statistics checksum mismatch')
    expected = {(row['source'],row['sv_class']):row for row in json.loads(stats.read_text())['rows']
                if row['source'] != 'combined' and row['sv_class'] != 'ALL'
                and not row['sv_class'].endswith('_literal_and_symbolic')}
    keys = sorted(expected)
    signature = json.dumps([1,identity(path),artifacts[path.name],artifacts[stats.name],sample_hash],sort_keys=True)
    cache = out/f'{chrom}.npz'
    if cache.exists():
        with np.load(cache,allow_pickle=False) as saved:
            if saved['signature'].item() == signature:
                return keys,saved['counts'],saved['records']
    print(f'Scanning {chrom}',flush=True)
    n = len(names)
    counts = np.zeros((len(keys),3,n),dtype=np.int64)
    records = np.zeros(len(keys),dtype=np.int64)
    index = {key:i for i,key in enumerate(keys)}
    with pysam.BGZFile(str(path),'r') as vcf:
        for line in vcf:
            if line.startswith(b'#') or b'UP_KEY=' not in line:
                continue
            fields = line.split(b'\t',8)
            tags = info_dict(fields[7].decode())
            key = tags.get('UP_KEY')
            if key is None:
                continue
            tier = 'popped' if key.startswith('P.') else 'lenient' if key.startswith('L.') else None
            kind = profile(int(fields[1]),fields[3].decode(),fields[4].decode(),tags,already_sv=True)[0]
            i = index[tier,kind]
            fmt,tail = fields[8].split(b'\t',1)
            copies,called = genotype_counts(fmt,tail,n)
            counts[i,0] += copies > 0
            counts[i,1] += copies
            counts[i,2] += called
            records[i] += 1
    for i,key in enumerate(keys):
        row = expected[key]
        if (records[i],counts[i,1].sum(),counts[i,2].sum()) != (row['allele_records'],row['AC_sum'],row['AN_sum']):
            raise ValueError(f'{chrom}: record count, AC or AN disagrees with completed panel for {key}')
    tmp = cache.with_suffix('.tmp.npz')
    np.savez_compressed(tmp,signature=signature,counts=counts,records=records)
    tmp.replace(cache)
    return keys,counts,records


def summarize(results,n):
    totals = {}
    denom = Counter()
    for keys,counts,records in results:
        for i,(tier,kind) in enumerate(keys):
            classes = [kind,'ALL']
            if kind in {'simple_DEL','symbolic_DEL','simple_INS','symbolic_INS'}:
                classes.append(kind.rsplit('_',1)[1]+'_literal_and_symbolic')
            for source in (tier,'combined'):
                for label in classes:
                    key = source,label
                    if key not in totals:
                        totals[key] = np.zeros((3,n),dtype=np.int64)
                    totals[key] += counts[i]
                    denom[key] += int(records[i])
    for key,(carried,copies,called) in sorted(totals.items()):
        p25,median,p75 = np.percentile(carried,[25,50,75])
        yield dict(source=key[0],sv_class=key[1],people=n,panel_records=denom[key],
                   median_SVs_carried=median,p25=p25,p75=p75,median_ALT_copies=np.median(copies),
                   people_with_observed_ALT=int(np.count_nonzero(carried)),
                   median_called_alleles_pct=round(100*np.median(called)/(2*denom[key]),3))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--panel',type=Path,default=Path.home()/'unified_popped_lenient_p70')
    parser.add_argument('--jobs',type=int,choices=range(1,5),default=4)
    args = parser.parse_args()
    os.umask(0o077)
    chroms = [f'chr{i}' for i in range(1,22)]
    missing = [c for c in chroms if not (args.panel/c/'done.json').is_file()]
    if missing:
        parser.error('Missing completed panel: '+', '.join(missing))
    out = args.panel/'per_person_summary'; out.mkdir(exist_ok=True)
    names = sample_header(args.panel/'chr1/chr1.unified.vcf.gz')
    sha = hashlib.sha256('\n'.join(names).encode()).hexdigest()
    with (out/'.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        with ProcessPoolExecutor(max_workers=args.jobs) as pool:
            rows = list(summarize(pool.map(count_chr,[(args.panel,out,c,sha) for c in chroms]),len(names)))
        target = out/'chr1_21_medians.tsv'
        write_tsv(target,rows)
        print('sv_class\tmedian_SVs_carried\tp25\tp75\tmedian_ALT_copies')
        for row in rows:
            if row['source']=='combined':
                print('\t'.join(str(row[k]) for k in ('sv_class','median_SVs_carried','p25','p75','median_ALT_copies')))
        print(f'Full aggregate results: {target}')
