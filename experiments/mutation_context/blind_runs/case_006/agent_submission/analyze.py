import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[key]='4'
import json, hashlib, time
from pathlib import Path
import numpy as np
from scipy.stats import chi2_contingency
ROOT=Path(__file__).resolve().parent
ALPH='ACGT'
LUT=np.full(256,255,dtype=np.uint8)
LUT[np.frombuffer(b'ACGT',dtype=np.uint8)]=np.arange(4)
def fasta(path):
    out={}; name=None; seq=[]
    with open(path) as f:
        for line in f:
            line=line.strip()
            if line.startswith('>'):
                if name is not None: out[name]=LUT[np.frombuffer(''.join(seq).encode(),dtype=np.uint8)]
                name=line[1:].split()[0]; seq=[]
            else: seq.append(line)
    if name is not None: out[name]=LUT[np.frombuffer(''.join(seq).encode(),dtype=np.uint8)]
    assert all(np.max(x)<4 for x in out.values())
    return out
def contexts(x,k):
    code=np.zeros(len(x),dtype=np.int32)
    for offset in range(-(k//2),k//2+1): code=4*code+np.roll(x,-offset)
    return code
def main():
    p=fasta(ROOT/'calibration_parent.fasta'); y=fasta(ROOT/'calibration_offspring.fasta'); new=fasta(ROOT/'new_parent.fasta')
    assert p.keys()==y.keys()
    rng=np.random.default_rng(0); names=list(p); perm=rng.permutation(len(names)); split=np.zeros(len(names),int);split[perm[:len(names)//5]]=1;split[perm[len(names)//5:2*(len(names)//5)]]=2
    ks=[1,3,5,7,9]; counts={k:np.zeros((3,4**k,4),dtype=np.int64) for k in ks}
    offsets=np.arange(-30,31); assoc=np.zeros((len(offsets),4,4,2),dtype=np.int64)
    sc=[];poscounts=np.zeros((100,2),dtype=np.int64);periods={k:np.zeros((k,2),dtype=np.int64) for k in range(2,31)}
    for i,(name,x) in enumerate(p.items()):
        z=y[name];assert len(x)==len(z);m=x!=z
        sc.append({'name':name,'n':len(x),'changes':int(m.sum()),'gc':float(np.mean((x==1)|(x==2))),'split':int(split[i])})
        for k in ks: counts[k][split[i]]+=np.bincount(contexts(x,k)*4+z,minlength=4**(k+1)).reshape(-1,4)
        if split[i]==0:
            for j,offset in enumerate(offsets):assoc[j]+=np.bincount((x.astype(int)*4+np.roll(x,-int(offset)))*2+m,minlength=32).reshape(4,4,2)
        bins=np.minimum(np.arange(len(x))*100//len(x),99)
        poscounts+=np.bincount(bins*2+m,minlength=200).reshape(100,2)
        for period in periods:periods[period]+=np.bincount(np.arange(len(x))%period*2+m,minlength=period*2).reshape(period,2)
    np.savez_compressed(ROOT/'context_counts.npz',**{f'k{k}':v for k,v in counts.items()},assoc=assoc,poscounts=poscounts)
    report={'n_scaffolds':len(p),'new_scaffolds':{n:len(x) for n,x in new.items()},'total_sites':sum(len(x) for x in p.values()),'total_changes':sum(s['changes'] for s in sc),'scaffolds':sc,'spectrum':counts[1].sum(axis=0).tolist(),'offsets':[],'periods':{}}
    for j,offset in enumerate(offsets):
        a=assoc[j];stat=sum(chi2_contingency(b)[0] for b in a if np.all(b.sum(axis=1)>0)) if offset!=0 else 0
        report['offsets'].append({'offset':int(offset),'chi2':float(stat),'rates':(a[:,:,1]/np.maximum(a.sum(axis=2),1)).tolist()})
    for k,a in periods.items():report['periods'][k]={'chi2':float(chi2_contingency(a)[0]),'rates':(a[:,1]/a.sum(axis=1)).tolist()}
    (ROOT/'exploration.json').write_text(json.dumps(report,indent=2))
    print(json.dumps({k:v for k,v in report.items() if k not in ['scaffolds','offsets','periods','new_scaffolds']},indent=2))
    print('Lengths:',sorted(set(len(x) for x in p.values())), 'new:',report['new_scaffolds'])
    print('Scaffold rate quantiles',np.quantile([s['changes']/s['n'] for s in sc],[0,.01,.1,.5,.9,.99,1]))
    print('Offsets:',[(a['offset'],round(a['chi2'],1)) for a in report['offsets']])
    print('Periods:',[(k,round(a['chi2'],1)) for k,a in report['periods'].items()])
if __name__=='__main__':main()
