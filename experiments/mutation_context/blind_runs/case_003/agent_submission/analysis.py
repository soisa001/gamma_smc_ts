import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']: os.environ[k]='4'
import numpy as np, json
from pathlib import Path
BASE=Path(__file__).parent
ABC=np.array(list('ACGT'))
def fasta(name):
 d={}; key=None
 for line in (BASE/name).read_text().splitlines():
  if line.startswith('>'): key=line[1:].split()[0]; d[key]=''
  else: d[key]+=line.strip()
 trans=np.full(256,255,dtype=np.uint8)
 for i,c in enumerate(b'ACGT'): trans[c]=i
 return {k:trans[np.frombuffer(v.encode(),dtype=np.uint8)] for k,v in d.items()}
def load():
 p=fasta('calibration_parent.fasta'); o=fasta('calibration_offspring.fasta'); n=fasta('new_parent.fasta')
 assert p.keys()==o.keys()
 assert all(len(v)==len(o[k]) and np.all(v<4) and np.all(o[k]<4) for k,v in p.items())
 return p,o,n
if __name__=='__main__':
 p,o,n=load(); out={}
 out['train_lengths']={k:len(v) for k,v in p.items()}; out['new_lengths']={k:len(v) for k,v in n.items()}
 out['train_scaffolds']=len(p); out['new_scaffolds']=len(n)
 x=np.concatenate(list(p.values())); y=np.concatenate(list(o.values())); m=x!=y
 out['substitution_counts']=np.bincount(x.astype(int)*4+y,minlength=16).reshape(4,4).tolist()
 out['total_changed']=int(m.sum()); out['total_bases']=len(x)
 out['scaffold_rates']={k:float((v!=o[k]).mean()) for k,v in p.items()}
 out['offset_assoc']={}
 for offset in range(-20,21):
  xx=np.concatenate([np.roll(v,-offset) for v in p.values()]); counts=np.bincount(xx,minlength=4); changes=np.bincount(xx,weights=m,minlength=4)
  out['offset_assoc'][offset]=(changes/counts).tolist()
 print(json.dumps(out,indent=2)); (BASE/'initial_statistics.json').write_text(json.dumps(out,indent=2))
