import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[k]='4'
import json, hashlib
from pathlib import Path
import numpy as np
from scipy.stats import chi2_contingency
ROOT=Path(__file__).resolve().parent
def fasta(name):
    out={}; name0=None; parts=[]
    for line in (ROOT/name).read_text().splitlines():
        if line.startswith('>'):
            if name0 is not None: out[name0]=''.join(parts)
            name0=line[1:].split()[0]; parts=[]
        else: parts.append(line.strip())
    if name0 is not None: out[name0]=''.join(parts)
    return out
def encode(s):
    lut=np.full(256,255,dtype=np.uint8)
    lut[np.frombuffer(b'ACGT',dtype=np.uint8)]=np.arange(4)
    a=lut[np.frombuffer(s.encode(),dtype=np.uint8)]
    assert np.all(a<4)
    return a
def contexts(a,r):
    z=np.zeros(len(a)-2*r,dtype=np.int32)
    for d in range(2*r+1):z=z*4+a[d:len(a)-2*r+d]
    return z
if __name__=='__main__':
    p=fasta('calibration_parent.fasta'); o=fasta('calibration_offspring.fasta'); n=fasta('new_parent.fasta')
    assert p.keys()==o.keys()
    print('contigs',len(p),len(n),'lengths',sorted(set(map(len,p.values()))), sorted(set(map(len,n.values()))))
    mats=[]; groups=[]; arrays=[]; yy=[]; ci=[]; pos=[]
    for i,(key,s) in enumerate(p.items()):
        a=encode(s); b=encode(o[key]); assert len(a)==len(b)
        m=np.bincount(a*4+b,minlength=16).reshape(4,4); mats.append(m)
        groups.append([key,len(a),int(np.sum(a!=b))])
        arrays.append(a); yy.append(b); ci.append(np.full(len(a),i));pos.append(np.arange(len(a)))
    print('substitution_counts',np.sum(mats,axis=0).tolist())
    print('contig_changes',groups)
    x=np.concatenate(arrays); y=np.concatenate(yy); c=np.concatenate(ci); positions=np.concatenate(pos); changed=x!=y
    print('total',len(x),'changed',int(changed.sum()),'rate',changed.mean())
    print('by_position_decile',[float(changed[(positions//(positions.max()+1)*10)==i].mean()) for i in []])
    for d in range(-12,13):
        if d==0:continue
        valid=(positions>=12)&(positions<positions.max()-12)
        ix=np.where(valid)[0]; feature=x[ix+d]
        counts=np.bincount(x[ix]*4+feature,minlength=16).reshape(4,4)
        changes=np.bincount(x[ix]*4+feature,weights=changed[ix],minlength=16).reshape(4,4)
        print('offset',d,'rates',np.round(changes/counts,5).tolist())
    for r in range(1,5):
        ctx=np.concatenate([contexts(a,r) for a in arrays]); tar=np.concatenate([b[r:-r] for b in yy]); cen=np.concatenate([a[r:-r] for a in arrays]); mut=cen!=tar
        counts=np.bincount(ctx,minlength=4**(2*r+1)); changes=np.bincount(ctx,weights=mut,minlength=len(counts))
        ix=np.argsort(changes/np.maximum(counts,1))[-20:][::-1]
        alphabet='ACGT'
        def decode(z):return ''.join(alphabet[(int(z)//4**k)%4] for k in range(2*r,-1,-1))
        print('context',r,'top',[(decode(z),int(counts[z]),int(changes[z]),round(changes[z]/counts[z],5)) for z in ix])
    for lag in [1,2,3,4,5,10,20,100]:
        valid=c[:-lag]==c[lag:]
        print('mutation_lag',lag,'conditional',float(changed[lag:][valid&changed[:-lag]].mean()),'baseline',changed.mean())
    print('new_composition',np.bincount(np.concatenate([encode(s) for s in n.values()]),minlength=4).tolist())
    with (ROOT/'provenance.json').open('w') as f:
        json.dump({'inputs':{k:hashlib.sha256((ROOT/k).read_bytes()).hexdigest() for k in ['TASK.md','metadata.json','calibration_parent.fasta','calibration_offspring.fasta','new_parent.fasta','queries.tsv']},'seed':0,'threads':4},f,indent=2)
