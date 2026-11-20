import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'): os.environ[k]='4'
import json, hashlib
from pathlib import Path
import numpy as np
from scipy.stats import chi2_contingency

ROOT=Path(__file__).parent
BASES='ACGT'
def fasta(name):
    out={}
    for line in (ROOT/name).read_text().splitlines():
        if line.startswith('>'): key=line[1:].split()[0]; out[key]=[]
        else: out[key].append(line.strip())
    return {k:np.frombuffer(''.join(v).encode(),dtype=np.uint8) for k,v in out.items()}
def encoded(d):
    lut=np.full(256,255,np.uint8)
    lut[np.frombuffer(b'ACGT',dtype=np.uint8)]=np.arange(4)
    return {k:lut[v] for k,v in d.items()}

if __name__=='__main__':
    p=encoded(fasta('calibration_parent.fasta')); o=encoded(fasta('calibration_offspring.fasta')); new=encoded(fasta('new_parent.fasta'))
    assert p.keys()==o.keys()
    result={'inputs':{},'scaffolds':[]}
    for f in ['TASK.md','metadata.json','calibration_parent.fasta','calibration_offspring.fasta','new_parent.fasta','queries.tsv']:
        result['inputs'][f]={'sha256':hashlib.sha256((ROOT/f).read_bytes()).hexdigest(),'bytes':(ROOT/f).stat().st_size}
    x=np.concatenate(list(p.values())); y=np.concatenate(list(o.values())); changed=x!=y
    assert x.max()<4 and y.max()<4
    result['n']=len(x); result['changed']=int(changed.sum()); result['rate']=float(changed.mean())
    result['spectrum']=np.bincount(x.astype(int)*4+y,minlength=16).reshape(4,4).tolist()
    for k,v in p.items():
        assert len(v)==len(o[k])
        result['scaffolds'].append({'id':k,'length':len(v),'changed':int((v!=o[k]).sum()),'gc':float(np.isin(v,[1,2]).mean())})
    result['new_scaffolds']={k:len(v) for k,v in new.items()}
    result['offset_associations']=[]
    # Stratify by source base, compare mutations to sequence at each circular offset.
    for offset in range(-30,31):
        if offset==0: continue
        z=np.concatenate([np.roll(v,-offset) for v in p.values()])
        table=np.bincount((x.astype(int)*4+z)*2+changed,minlength=32).reshape(4,4,2)
        chisq=sum(chi2_contingency(t,correction=False)[0] for t in table)
        result['offset_associations'].append({'offset':offset,'chi2':chisq,'counts':table.tolist()})
    result['position_moduli']=[]
    pos=np.concatenate([np.arange(len(v)) for v in p.values()])
    for mod in range(2,41):
        table=np.bincount((pos%mod)*2+changed,minlength=mod*2).reshape(mod,2)
        result['position_moduli'].append({'mod':mod,'chi2':chi2_contingency(table,correction=False)[0], 'rates':(table[:,1]/table.sum(1)).tolist()})
    (ROOT/'exploration.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k not in ['inputs','offset_associations','position_moduli']},indent=2))
    print('Offsets', [(r['offset'],round(r['chi2'],1)) for r in result['offset_associations']])
    print('Moduli', [(r['mod'],round(r['chi2'],1)) for r in result['position_moduli']])
