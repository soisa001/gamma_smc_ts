"""Sample length-weighted autosomal windows; preserve the deCODE map intervals."""
from __future__ import annotations
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
import pyBigWig


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def merge(intervals):
    result = []
    for left, right in sorted(intervals):
        if right <= left:
            continue
        if result and left <= result[-1][1]:
            result[-1][1] = max(right, result[-1][1])
        else:
            result.append([left, right])
    return result


def region_arrays(intervals, excluded, start, length, mutation_rate):
    end = start + length
    pieces = [(max(start,l)-start, min(end,r)-start, float(v)*1e-8)
              for l,r,v in intervals if r > start and l < end and np.isfinite(v)]
    positions = np.unique([0, length] + [x for l,r,v in pieces for x in (l,r)]).astype(np.int64)
    rates = np.full(len(positions)-1, 1e-8)
    covered = np.zeros(len(rates), dtype=bool)
    for l,r,v in pieces:
        assert v >= 0
        lo, hi = np.searchsorted(positions, [l,r])
        rates[lo:hi] = v
        covered[lo:hi] = True
    missing = [(int(l),int(r)) for l,r,ok in zip(positions[:-1],positions[1:],covered) if not ok]
    bad = merge(missing + [(max(start,l)-start,min(end,r)-start) for l,r in excluded if r>start and l<end])
    mpos = np.unique([0,length] + [x for l,r in bad for x in (l,r)]).astype(np.int64)
    mrate = np.full(len(mpos)-1, mutation_rate)
    for l,r in bad:
        lo,hi = np.searchsorted(mpos,[l,r])
        mrate[lo:hi] = 0
    width = np.diff(positions)
    return dict(recombination_position=positions, recombination_rate=rates,
                mutation_position=mpos, mutation_rate=mrate), dict(
        map_coverage=float(width[covered].sum()/length), callable_fraction=float(np.diff(mpos)[mrate>0].sum()/length),
        mean_recombination_rate=float(np.sum(width*rates)/length),
        map_min_rate=float(rates.min()), map_max_rate=float(rates.max()))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--maps', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--replicates', type=int, default=1000)
    p.add_argument('--length', type=int, default=11000000)
    p.add_argument('--seed', type=int, default=20380101)
    p.add_argument('--minimum-callable', type=float, default=.9)
    p.add_argument('--mutation-rate', type=float, default=1.29e-8)
    args=p.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    sources={name:digest(args.maps/name) for name in ('recombAvg.bw','gap.txt.gz','hg38.chrom.sizes')}
    identity=dict(sources=sources, seed=args.seed, length=args.length, replicates=args.replicates,
        minimum_callable=args.minimum_callable, mutation_rate=args.mutation_rate,
        algorithm='uniform-autosomal-start/reject-low-coverage/v2', rng_domain=0x52454749)
    receipt=args.out/'manifest.json'
    if receipt.exists():
        r=json.loads(receipt.read_text())
        assert r['identity']==identity
        for name,sha in r['files'].items():
            assert digest(args.out/name)==sha, f'Corrupt region resource: {name}'
        print('Validated cached genomic regions',flush=True)
        return
    chroms=[f'chr{i}' for i in range(1,23)]
    sizes={x.split()[0]:int(x.split()[1]) for x in (args.maps/'hg38.chrom.sizes').read_text().splitlines()}
    gaps={c:[] for c in chroms}
    with gzip.open(args.maps/'gap.txt.gz','rt') as f:
        for line in f:
            x=line.split()
            if x[1] in gaps:
                gaps[x[1]].append((int(x[2]),int(x[3])))
    bw=pyBigWig.open(str(args.maps/'recombAvg.bw'))
    assert all(bw.chroms()[c]==sizes[c] for c in chroms)
    nstarts=np.array([sizes[c]-args.length+1 for c in chroms],dtype=np.int64)
    cumulative=np.cumsum(nstarts)
    # Independent of demographic normals even though both retain the study seed.
    rng=np.random.default_rng(np.random.SeedSequence([args.seed, 0x52454749]))
    rows=[]
    proposals=0
    for rep in range(args.replicates):
        while True:
            proposals+=1
            draw=int(rng.integers(cumulative[-1]))
            ci=int(np.searchsorted(cumulative,draw,side='right'))
            chrom=chroms[ci]
            start=draw-(int(cumulative[ci-1]) if ci else 0)
            end=start+args.length
            intervals=bw.intervals(chrom,start,end) or ()
            arrays,metrics=region_arrays(intervals,gaps[chrom],start,args.length,args.mutation_rate)
            if metrics['callable_fraction'] >= args.minimum_callable and metrics['map_coverage'] >= args.minimum_callable:
                break
        target=args.out/f'region{rep:04d}.npz'
        with target.with_suffix('.tmp').open('wb') as f:
            np.savez_compressed(f,**arrays)
        target.with_suffix('.tmp').replace(target)
        rows.append(dict(replicate=rep, chromosome=chrom, start_0based=start,end_exclusive=end,
                         map_file=target.name,sha256=digest(target),**metrics))
    bw.close()
    pd.DataFrame(rows).to_csv(args.out/'regions.csv',index=False)
    record=dict(identity=identity, proposals=proposals, accepted=len(rows),
        chromosome_counts=pd.DataFrame(rows).chromosome.value_counts().to_dict(),
        genome_build='GRCh38', recombination_source='https://hgdownload.soe.ucsc.edu/gbdb/hg38/recombRate/recombAvg.bw',
        rate_conversion='cM/Mb multiplied by 1e-8 = per-base per-generation',
        masking='assembly gaps and missing recombination-map coverage; HMMix strict mask is NOT used',
        missing_recombination_rate=1e-8, mutation_exposure='1.29e-8 on included bases; excluded bases have no observed mutations and are masked in decoding',
        focal_window='11 Mb simulated; focal chosen near middle at onset; crop centered 10 Mb after simulation',
        files={x.name:digest(x) for x in args.out.iterdir() if x.is_file() and x.name!='manifest.json'})
    receipt.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({k:v for k,v in record.items() if k!='files'}),flush=True)


if __name__=='__main__':
    main()
