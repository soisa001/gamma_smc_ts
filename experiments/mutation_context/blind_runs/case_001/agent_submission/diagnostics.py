import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='4'
import numpy as np,json
from scipy.stats import chi2, chi2_contingency
from model import load,features,counts,rc_codes
from explore import ROOT
ctx,y,c=load();f,t,rev=features(ctx,y);x=(ctx//64)%4;mut=x!=y
cc=counts(f['cpg'],t,3);pr=cc/cc.sum(1)[:,None]
res={}
for name in ['context3','context5','context7','raw_context3','raw_context5','raw_context7']:
    raw=name.startswith('raw');ct=counts(f[name],y if raw else t,int(f[name].max()+1)); ix=np.unique(f[name],return_index=True)[1]
    p0=np.empty_like(ct);p0[f[name][ix]]=pr[f['cpg'][ix]]
    if raw:p0[f[name][ix]]=np.where(rev[ix,None],p0[f[name][ix]][:,::-1],p0[f[name][ix]])
    n=ct.sum(1);valid=n>0;e=n[valid,None]*p0[valid];obs=ct[valid]
    stat=np.sum((obs-e)**2/e);df=(valid.sum()-3)*3
    res[name]={'pearson_stat':float(stat),'df':int(df),'dispersion':float(stat/df)}
    if name=='context5':
        idx=np.argsort(np.max((obs-e)/np.sqrt(e),axis=1))[-30:][::-1]; valid_codes=np.where(valid)[0]
        for j in idx:
            z=valid_codes[j];s=''.join('ACGT'[(z//4**k)%4] for k in range(4,-1,-1))
            print(s,obs[j].astype(int).tolist(),np.round(e[j],2).tolist())
expected=1-pr[f['cpg'],(np.where(rev,3-x,x))]
for name,g in [('contig',c),('position_decile',np.tile(np.arange(49994)//5000,80)),('position_mod3',np.tile(np.arange(3,49997)%3,80))]:
    obs=np.bincount(g,weights=mut);exp=np.bincount(g,weights=expected);stat=np.sum((obs-exp)**2/exp)
    res[name]={'observed':obs.tolist(),'expected':exp.tolist(),'pearson_stat':float(stat),'df':len(obs)-1,'p':float(chi2.sf(stat,len(obs)-1))}
# Test strand symmetry in each CpG class.
res['strand_tests']=[]
for k in range(3):
    tab=np.array([np.bincount(t[(f['cpg']==k)&(rev==r)],minlength=4) for r in [0,1]])
    st,p,df,e=chi2_contingency(tab)
    res['strand_tests'].append({'class':k,'counts':tab.tolist(),'p':float(p),'stat':float(st)})
(ROOT/'diagnostics.json').write_text(json.dumps(res,indent=2))
print(json.dumps(res,indent=2))
