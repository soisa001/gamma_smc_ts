import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[k]='4'
import json, hashlib
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
BASES='ACGT'
def fasta(path):
    result={}; name=None; chunks=[]
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            if name is not None: result[name]=''.join(chunks)
            name=line[1:].split()[0]; chunks=[]
        else: chunks.append(line.strip())
    if name is not None: result[name]=''.join(chunks)
    return result
def encode(seq):
    lut=np.full(256,255,dtype=np.uint8)
    lut[np.frombuffer(b'ACGT',dtype=np.uint8)]=np.arange(4,dtype=np.uint8)
    x=lut[np.frombuffer(seq.encode(),dtype=np.uint8)]
    assert (x<4).all()
    return x
def inputs():
    a=fasta(ROOT/'calibration_parent.fasta'); b=fasta(ROOT/'calibration_offspring.fasta')
    assert a.keys()==b.keys()
    names=list(a)
    xs=[encode(a[n]) for n in names]; ys=[encode(b[n]) for n in names]
    assert all(len(x)==len(y) for x,y in zip(xs,ys))
    return names,xs,ys
if __name__=='__main__':
    names,xs,ys=inputs(); x=np.concatenate(xs); y=np.concatenate(ys); mut=x!=y
    mat=np.bincount(x.astype(int)*4+y,minlength=16).reshape(4,4)
    out={'scaffolds':len(names),'lengths':[len(s) for s in xs],'n':len(x),'changes':int(mut.sum()),'rate':float(mut.mean()),'counts':mat.tolist(),'matrix':(mat/mat.sum(1,keepdims=True)).tolist(),'per_scaffold':[{'name':n,'n':len(a),'mutations':int((a!=b).sum()),'GC':float(np.isin(a,[1,2]).mean())} for n,a,b in zip(names,xs,ys)]}
    new=fasta(ROOT/'new_parent.fasta');out['new_lengths']={n:len(s) for n,s in new.items()}
    out['sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.iterdir() if p.name in ['TASK.md','metadata.json','calibration_parent.fasta','calibration_offspring.fasta','new_parent.fasta','queries.tsv']}
    (ROOT/'input_audit.json').write_text(json.dumps(out,indent=2))
    print(json.dumps(out,indent=2))
