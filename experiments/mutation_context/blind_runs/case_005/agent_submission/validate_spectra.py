import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'): os.environ[k]='4'
import json
import numpy as np
from scipy.stats import chi2, chi2_contingency
from explore import fasta, encoded, ROOT
from validate import features, metrics

p=encoded(fasta('calibration_parent.fasta')); o=encoded(fasta('calibration_offspring.fasta'))
seqs=list(p.values()); x=np.concatenate(seqs); y=np.concatenate(list(o.values())); changed=x!=y; xor=x^y
group=np.concatenate([np.full(len(v),i,np.int16) for i,v in enumerate(seqs)])
previous=json.loads((ROOT/'validation.json').read_text()); groupfold=np.array(previous['fold_assignment_scaffolds'])
f=features(seqs)
records=[]
for scheme, fold in [('scaffold',groupfold[group]),('contiguous_region',np.concatenate([np.arange(len(v))*5//len(v) for v in seqs]))]:
    baseline=next(r for r in previous['models'] if r['scheme']==scheme and r['model']=='pooled')
    for name in ['dinuc_-1','dinuc_1','context_3','context_rc_3','context_5','context_rc_5','context_rc_7']:
        code=f[name]; ncat=int(code.max())+1
        for alpha in ([0,1000,10000] if ncat<=64 else [10000,100000]):
            pro=np.empty((len(x),4))
            for k in range(5):
                train=fold!=k; test=~train; r=changed[train].mean()
                counts=np.bincount(code[train].astype(int)*4+xor[train],minlength=ncat*4).reshape(ncat,4)
                prior=np.array([1-r,r/3,r/3,r/3])
                probs=(counts+alpha*prior)/(counts.sum(1,keepdims=True)+alpha)
                probs=np.maximum(probs,1e-12); probs/=probs.sum(1,keepdims=True)
                px=probs[code[test]]
                pro[test]=np.take_along_axis(px,x[test,None]^np.arange(4)[None,:],axis=1)
            rec={'scheme':scheme,'model':name+'_full_spectrum','prior_strength':alpha,**metrics(pro,x,y,group)}
            rec['nll_delta_vs_pooled']=rec['nll_total']-baseline['nll_total']; records.append(rec)
            print(scheme,name,alpha,'delta NLL',round(rec['nll_delta_vs_pooled'],4),flush=True)
checks=[]
for off in range(-30,31):
    if not off: continue
    z=np.concatenate([np.roll(v,-off) for v in seqs])[changed]
    tab=np.bincount((x[changed].astype(int)*4+z)*4+y[changed],minlength=64).reshape(4,4,4)
    stat=sum(chi2_contingency(t[:,np.arange(4)!=b],correction=False)[0] for b,t in enumerate(tab))
    checks.append({'offset':off,'chi2':stat,'df':24,'p':float(chi2.sf(stat,24))})
(ROOT/'validation_spectra.json').write_text(json.dumps({'models':records,'offset_spectrum_associations':checks},indent=2))
print('Smallest offset p-values',sorted(checks,key=lambda v:v['p'])[:6])
