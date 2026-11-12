import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='4'
import numpy as np
from pathlib import Path
import json
from explore import fasta,encode,contexts,ROOT

def rc_codes(z,k):
    out=np.zeros_like(z)
    for j in range(k):out=out*4+(3-(z//(4**j))%4)
    return out

def load():
    p=fasta('calibration_parent.fasta');o=fasta('calibration_offspring.fasta')
    names=list(p);a=[encode(p[k]) for k in names];b=[encode(o[k]) for k in names]
    ctx=np.concatenate([contexts(v,3) for v in a]); y=np.concatenate([v[3:-3] for v in b]);c=np.repeat(np.arange(len(a)),[len(v)-6 for v in a])
    return ctx,y,c

def features(ctx,y=None):
    x=(ctx//64)%4
    reverse=x>=2
    # Orient A and C as reference alleles, complementing G/T and reversing context.
    z=np.where(reverse,rc_codes(ctx,7),ctx)
    target=None if y is None else np.where(reverse,3-y,y)
    center=(z//64)%4
    cpg=(center==1)&((z//16)%4==2)
    f={'base':x,'cpg_unpooled':x*2+(((x==1)&((ctx//16)%4==2))|((x==2)&((ctx//256)%4==1))),
       'symmetric_base':center,'cpg':center+cpg.astype(int)}
    for r in [1,2,3]:
        f[f'context{2*r+1}']=((z//(4**(3-r)))%(4**(2*r+1))).astype(int)
        f[f'raw_context{2*r+1}']=((ctx//(4**(3-r)))%(4**(2*r+1))).astype(int)
    return f,target,reverse

def counts(f,y,n):return np.bincount(f*4+y,minlength=n*4).reshape(n,4).astype(float)
def smooth(count,prior,strength):return (count+strength*prior)/(count.sum(1)[:,None]+strength)

def main():
    ctx,y,c=load();f,t,reverse=features(ctx,y)
    rng=np.random.default_rng(0); perm=rng.permutation(80); fold=np.empty(80,int);fold[perm]=np.arange(80)%5
    # All model selection is by held-out whole contigs.
    configs=[('base',0),('symmetric_base',0),('cpg_unpooled',0),('cpg',0)]
    for name in ['context3','context5','context7','raw_context3','raw_context5','raw_context7']:
        for strength in [10,100,500,2000,10000,50000,1e6]:configs.append((name,strength))
    scores={str(z):[] for z in configs}; llcont={str(z):np.zeros(80) for z in configs}
    sizes={k:int(v.max()+1) for k,v in f.items()}
    for fi in range(5):
        train=fold[c]!=fi;test=~train
        cc=counts(f['cpg'][train],t[train],3); prior=(cc+.5)/(cc.sum(1)[:,None]+2)
        for name,strength in configs:
            raw=(name in ['base','cpg_unpooled']) or name.startswith('raw_')
            target=y if raw else t
            ct=counts(f[name][train],target[train],sizes[name])
            if strength==0:p=(ct+.5)/(ct.sum(1)[:,None]+2)
            else:
                # Map each context class to parent CpG model probabilities via observed membership.
                ix=np.unique(f[name],return_index=True)[1]
                pp=np.empty((sizes[name],4));pp[:]=.25
                pp[f[name][ix]]=prior[f['cpg'][ix]]
                if raw:
                    pp[f[name][ix]]=np.where(reverse[ix,None],pp[f[name][ix]][:,::-1],pp[f[name][ix]])
                p=smooth(ct,pp,strength)
            pred=p[f[name][test]]
            losses=-np.log(pred[np.arange(test.sum()),target[test]])
            scores[str((name,strength))].append(float(losses.mean()))
            llcont[str((name,strength))]+=np.bincount(c[test],weights=losses,minlength=80)
        print('finished fold',fi,flush=True)
    results=[{'model':name,'strength':strength,'mean_log_loss':float(np.mean(scores[str((name,strength))])),'fold_log_losses':scores[str((name,strength))]} for name,strength in configs]
    results.sort(key=lambda z:z['mean_log_loss'])
    (ROOT/'cv_results.json').write_text(json.dumps(results,indent=2))
    np.savez_compressed(ROOT/'cv_contig_losses.npz',**llcont)
    print(json.dumps(results[:20],indent=2))
    # Test residual dispersion at every resolution after accounting for central base/CpG.
    cc=counts(f['cpg'],t,3); prior=cc/cc.sum(1)[:,None]
    print('cpg_counts',cc.tolist(),'cpg_probabilities',prior.tolist())
    for name in ['context3','context5','context7']:
        ct=counts(f[name],t,sizes[name]); ix=np.unique(f[name],return_index=True)[1]
        n=ct.sum(1);p0=np.empty_like(ct);p0[f[name][ix]]=prior[f['cpg'][ix]]
        valid=n>0
        stat=np.sum(((ct[valid]-n[valid,None]*p0[valid])**2)/(n[valid,None]*p0[valid]))
        print('dispersion',name,stat,'df',(valid.sum()-3)*3,'ratio',stat/((valid.sum()-3)*3))

if __name__=='__main__':main()
