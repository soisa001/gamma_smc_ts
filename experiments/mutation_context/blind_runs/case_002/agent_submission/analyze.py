import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='4'
import resource
resource.setrlimit(resource.RLIMIT_AS,(26*1024**3,26*1024**3))
import json, hashlib, sys
from pathlib import Path
import numpy as np
from scipy.special import xlogy
np.random.seed(0)
ROOT=Path(__file__).resolve().parent
BASES='ACGT'
LUT=np.full(256,255,dtype=np.uint8)
for i,b in enumerate(BASES): LUT[ord(b)]=i

def fasta(name):
    seq={}; current=None; chunks=[]
    for line in (ROOT/name).read_text().splitlines():
        if line.startswith('>'):
            if current is not None: seq[current]=LUT[np.frombuffer(''.join(chunks).encode(),dtype=np.uint8)]
            current=line[1:].split()[0]; assert current not in seq; chunks=[]
        else: chunks.append(line.strip())
    if current is not None: seq[current]=LUT[np.frombuffer(''.join(chunks).encode(),dtype=np.uint8)]
    assert all(np.all(v<4) for v in seq.values())
    return seq

def context(seq,r):
    n=len(seq); c=np.zeros(n-2*r,dtype=np.int64)
    for offset in range(2*r+1): c=c*4+seq[offset:n-2*r+offset]
    return c

def main():
    p=fasta('calibration_parent.fasta'); y=fasta('calibration_offspring.fasta'); new=fasta('new_parent.fasta')
    assert p.keys()==y.keys()
    assert all(len(p[k])==len(y[k]) for k in p)
    rows=[]; mat=np.zeros((4,4),dtype=np.int64); blocks=[]
    for k in p:
        x=p[k]; z=y[k]; changes=x!=z
        mat+=np.bincount(x.astype(int)*4+z,minlength=16).reshape(4,4)
        rows.append({'contig':k,'length':len(x),'changes':int(changes.sum()),'gc':float(np.mean((x==1)|(x==2)))})
        for a in range(0,len(x),10000): blocks.append([k,a,int(changes[a:a+10000].sum()),len(x[a:a+10000])])
    inputs=['TASK.md','metadata.json','queries.tsv','calibration_parent.fasta','calibration_offspring.fasta','new_parent.fasta']
    report={'input_sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in inputs},'calibration_contigs':rows,'substitution_counts':mat.tolist(),'new_contigs':{k:len(v) for k,v in new.items()},'blocks_10000':blocks}
    (ROOT/'input_diagnostics.json').write_text(json.dumps(report,indent=2))
    print('calibration contigs',len(p),'lengths',np.unique([len(v) for v in p.values()],return_counts=True))
    print('new contigs',len(new),'lengths',np.unique([len(v) for v in new.values()],return_counts=True))
    print('matrix',mat,'rates',mat/mat.sum(1)[:,None],sep='\n')
    print('contig mutation min/median/max',np.quantile([r['changes']/r['length'] for r in rows],[0,.5,1]))
    print('blocks min median max',np.quantile([b[2]/b[3] for b in blocks],[0,.5,1]))
    qn=0; qmin=10**10; qmax=0; seen=set(); last=None; duplicates=0
    with (ROOT/'queries.tsv').open() as f:
        assert next(f).strip()=='contig\tposition'
        for line in f:
            k,pos=line.rstrip().split('\t'); pos=int(pos)
            assert k in new and 1<=pos<=len(new[k]); qn+=1; qmin=min(qmin,pos);qmax=max(qmax,pos)
            if (k,pos)==last: duplicates+=1
            last=(k,pos);seen.add(k)
    print('queries',qn,'contigs',len(seen),'coordinate bounds',qmin,qmax,'adjacent duplicates',duplicates)
    report['queries']={'count':qn,'contigs':len(seen),'min':qmin,'max':qmax,'adjacent_duplicates':duplicates}
    (ROOT/'input_diagnostics.json').write_text(json.dumps(report,indent=2))
    for r in [1,2,3]:
        counts=np.zeros((4**(2*r+1),4),dtype=np.int64)
        for k in p:
            c=context(p[k],r); z=y[k][r:-r]
            counts+=np.bincount(c*4+z,minlength=counts.size).reshape(counts.shape)
        np.save(ROOT/f'context_counts_r{r}.npy',counts)
        total=counts.sum(1); parent=(np.arange(len(total))//4**r)%4
        mutations=total-counts[np.arange(len(total)),parent]
        top=np.argsort(mutations/np.maximum(total,1))[-15:][::-1]
        print('radius',r,'top mutation contexts')
        for i in top:
            ctx=''.join(BASES[(i//4**j)%4] for j in range(2*r,-1,-1))
            print(ctx,int(total[i]),int(mutations[i]),np.round(counts[i]/max(total[i],1),4))

if __name__=='__main__': main()
