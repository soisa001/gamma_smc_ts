"""Deterministic, full-population validation of candidate probability models."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'): os.environ[k]='4'
import json
from pathlib import Path
import numpy as np
from scipy.stats import chi2, chi2_contingency, binomtest
from explore import fasta, encoded, ROOT

def contexts(seqs, radius):
    codes=[]
    for v in seqs:
        code=np.zeros(len(v),np.int64)
        for off in range(-radius,radius+1): code=code*4+np.roll(v,-off)
        codes.append(code)
    return np.concatenate(codes)

def features(seqs):
    x=np.concatenate(seqs)
    f={'pooled':np.zeros(len(x),np.int8),'source':x, 'source_gc':np.isin(x,[1,2]).astype(np.int8)}
    for r in [1,2,3]:
        code=contexts(seqs,r)
        f[f'context_{2*r+1}']=code
        rc=np.zeros(4**(2*r+1),np.int64); tmp=np.arange(len(rc))
        for _ in range(2*r+1): rc=rc*4+3-tmp%4; tmp//=4
        canon=np.minimum(np.arange(len(rc)),rc)
        _, remap=np.unique(canon,return_inverse=True)
        f[f'context_rc_{2*r+1}']=remap[code]
    for off in [-1,1]: f[f'dinuc_{off}']=x.astype(np.int16)*4+np.concatenate([np.roll(v,-off) for v in seqs])
    pos=np.concatenate([np.arange(len(v)) for v in seqs])
    for bins in [10,50,100]: f[f'position_bins_{bins}']=np.concatenate([np.arange(len(v))*bins//len(v) for v in seqs])
    for mod in [2,3,7,13,39]: f[f'position_mod_{mod}']=pos%mod
    for radius in [5,25,100,500]:
        gc=[]
        for v in seqs:
            g=np.isin(v,[1,2]).astype(np.int32)
            padded=np.concatenate([g[-radius:],g,g[:radius]])
            cs=np.r_[0,np.cumsum(padded)]
            vals=cs[2*radius+1:]-cs[:-(2*radius+1)]
            gc.append(np.minimum(9,(vals*10//(2*radius+1))).astype(np.int8))
        f[f'gc_window_{2*radius+1}']=np.concatenate(gc)
    f['homopolymer_neighbor']=np.concatenate([(np.roll(v,1)==v).astype(np.int8)+(np.roll(v,-1)==v).astype(np.int8) for v in seqs])
    return f

def metrics(probs,x,y,groups):
    n=len(y); change=x!=y
    ll=-np.log(probs[np.arange(n),y]); br=(probs**2).sum(1)-2*probs[np.arange(n),y]+1
    predchange=1-probs[np.arange(n),x]
    return {'log_loss':float(ll.mean()),'nll_total':float(ll.sum()),'brier':float(br.mean()),'expected_changes':float(predchange.sum()),'observed_changes':int(change.sum()),'group_nll':np.bincount(groups,weights=ll).tolist()}

def evaluate(f,x,y,fold,groups,scheme):
    records=[]; base=None
    for name,code in f.items():
        for alpha in ([0] if int(code.max())<100 else [1000,10000,100000]):
            ncat=int(code.max())+1
            pro=np.empty((len(x),4))
            for k in range(5):
                train=fold!=k; test=~train
                rate=float((x[train]!=y[train]).mean())
                nt=np.bincount(code[train],minlength=ncat)
                nm=np.bincount(code[train],weights=x[train]!=y[train],minlength=ncat)
                rates=np.divide(nm+alpha*rate,nt+alpha,out=np.full(ncat,rate),where=(nt+alpha)>0)
                r=np.clip(rates[code[test]],1e-10,1-1e-10)
                pt=np.repeat((r/3)[:,None],4,axis=1); pt[np.arange(test.sum()),x[test]]=1-r
                pro[test]=pt
            rec={'scheme':scheme,'model':name,'prior_strength':alpha,**metrics(pro,x,y,groups)}
            if base is None: base=rec
            rec['nll_delta_vs_pooled']=rec['nll_total']-base['nll_total']
            records.append(rec)
            print(scheme,name,alpha,'delta NLL',round(rec['nll_delta_vs_pooled'],4),flush=True)
    # Richer substitution spectra; XOR yields three universally non-reference classes.
    for name in ['source_spectrum','xor_spectrum','transition_spectrum','source_rate_xor_spectrum']:
        pro=np.empty((len(x),4))
        for k in range(5):
            train=fold!=k; test=~train
            c=x[train]!=y[train]; rate=c.mean(); xt=x[test]
            if name=='source_spectrum':
                tab=np.bincount(x[train].astype(int)*4+y[train],minlength=16).reshape(4,4)+0.5
                pt=(tab/tab.sum(1,keepdims=True))[xt]
            else:
                bins=np.bincount(x[train][c]^y[train][c],minlength=4)[1:]+0.5
                spec=bins/bins.sum()
                if name=='transition_spectrum': spec[[0,2]]=spec[[0,2]].mean()
                if name=='source_rate_xor_spectrum':
                    rates=np.bincount(x[train],weights=c,minlength=4)/np.bincount(x[train],minlength=4)
                    r=rates[xt]
                else: r=np.full(test.sum(),rate)
                pt=np.empty((test.sum(),4))
                for a in range(4):
                    dest=xt^a
                    pt[np.arange(len(xt)),dest]=(1-r) if a==0 else r*spec[a-1]
            pro[test]=pt
        rec={'scheme':scheme,'model':name,'prior_strength':None,**metrics(pro,x,y,groups)}
        rec['nll_delta_vs_pooled']=rec['nll_total']-base['nll_total']; records.append(rec)
        print(scheme,name,'delta NLL',round(rec['nll_delta_vs_pooled'],4),flush=True)
    return records

if __name__=='__main__':
    p=encoded(fasta('calibration_parent.fasta')); o=encoded(fasta('calibration_offspring.fasta'))
    seqs=list(p.values()); x=np.concatenate(seqs); y=np.concatenate(list(o.values())); change=x!=y
    group=np.concatenate([np.full(len(v),i,np.int16) for i,v in enumerate(seqs)])
    rng=np.random.default_rng(0); perm=rng.permutation(len(seqs)); groupfold=np.zeros(len(seqs),int); groupfold[perm]=np.arange(len(seqs))%5
    folds={'scaffold':groupfold[group],'contiguous_region':np.concatenate([np.arange(len(v))*5//len(v) for v in seqs])}
    f=features(seqs)
    records=[]
    for scheme,fold in folds.items(): records+=evaluate(f,x,y,fold,group,scheme)
    (ROOT/'validation.json').write_text(json.dumps({'seed':0,'fold_assignment_scaffolds':groupfold.tolist(),'models':records},indent=2))
    diagnostics={}
    for name in ['source','source_gc','position_bins_10','gc_window_1001','homopolymer_neighbor']:
        code=f[name]; ncat=int(code.max())+1
        tab=np.bincount(code.astype(int)*2+change,minlength=ncat*2).reshape(ncat,2); tab=tab[tab.sum(1)>0]
        stat,pvalue,dof,_=chi2_contingency(tab,correction=False)
        diagnostics[name]={'chi2':float(stat),'p':float(pvalue),'dof':dof,'counts':tab.tolist()}
    nc=np.bincount(group,weights=change); n=np.bincount(group); rate=change.mean()
    stat=np.sum((nc-n*rate)**2/(n*rate*(1-rate)))
    diagnostics['scaffold_heterogeneity']={'chi2':stat,'df':len(n)-1,'p':float(chi2.sf(stat,len(n)-1)),'mean_changes':nc.mean(),'variance_changes':nc.var(ddof=1)}
    spec=np.bincount((x[change]^y[change]),minlength=4)[1:]; stat=np.sum((spec-spec.mean())**2/spec.mean())
    diagnostics['xor_spectrum']={'counts':spec.tolist(),'chi2':stat,'p':float(chi2.sf(stat,2))}
    tab=np.bincount(x.astype(int)*4+y,minlength=16).reshape(4,4); arr=tab[~np.eye(4,dtype=bool)]
    exposure=np.repeat(np.bincount(x,minlength=4),3)*rate/3
    stat=np.sum((arr-exposure)**2/exposure)
    diagnostics['all_12_substitution_rates']={'chi2':stat,'df':11,'p':float(chi2.sf(stat,11))}
    pairs=[]
    for lag in [1,2,3,5,10,50,100,1000]:
        k=sum(int(np.sum((v!=o[key]) & np.roll(v!=o[key],-lag))) for key,v in p.items())
        pairs.append({'lag':lag,'observed_pairs':k,'iid_expected_pairs':len(x)*rate**2})
    diagnostics['mutation_pair_counts']=pairs
    ci=binomtest(int(change.sum()),len(change)).proportion_ci(confidence_level=.95,method='exact')
    diagnostics['rate_ci95']=[ci.low,ci.high]
    (ROOT/'diagnostics.json').write_text(json.dumps(diagnostics,indent=2))
    print(json.dumps(diagnostics,indent=2))
