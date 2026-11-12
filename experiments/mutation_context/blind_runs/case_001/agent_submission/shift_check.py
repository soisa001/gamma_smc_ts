import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='4'
import json,numpy as np
from scipy.stats import chi2
from explore import ROOT
from model import load,counts,rc_codes
from advanced import motif_features,tied_fit
ctx,y,c=load();g,eq,f,rev=motif_features(ctx);x=(ctx//64)%4;out=y^x;z=np.where(rev,rc_codes(ctx,7),ctx)
p=tied_fit(counts(g,out,4));q=json.loads((ROOT/'summary.json').read_text())['query_group_exposures'];weights=np.array(q)/np.bincount(g,minlength=4)
features={'shared_base_minus2':g*4+(z//1024)%4,'left1':g*4+(z//256)%4,'right1':g*4+(z//16)%4,'left3':g*4+(z//4096)%4,'right3':g*4+z%4}
rng=np.random.default_rng(0);perm=rng.permutation(80);fold=np.empty(80,int);fold[perm]=np.arange(80)%5
res={}
for name,feat in features.items():
    ct=counts(feat,out,16);pr=ct.sum(1)[:,None]*p[np.arange(16)//4];valid=ct.sum(1)>0;df=int(valid.sum()*3-5);stat=float(np.sum((ct[valid]-pr[valid])**2/pr[valid]));res[name]={'counts':ct.tolist(),'pearson_stat':stat,'df':df,'p':float(chi2.sf(stat,df))}
configs=['tied5','tied6','motif4',*features]
loss={name:np.zeros(80) for name in configs}
for fi in range(5):
    tr=fold[c]!=fi;te=~tr;ct=counts(g[tr],out[tr],4)
    for name in configs:
        if name.startswith('tied'):pred=tied_fit(ct,name)[g[te]]
        elif name=='motif4':pred=((ct+.5)/(ct.sum(1)[:,None]+2))[g[te]]
        else:
            cc=counts(features[name][tr],out[tr],16);pred=((cc+.5)/(cc.sum(1)[:,None]+2))[features[name][te]]
        ls=-np.log(pred[np.arange(te.sum()),out[te]])*weights[g[te]]
        loss[name]+=np.bincount(c[te],weights=ls,minlength=80)
for name in configs:
    delta=loss[name]-loss['tied5'];res[name]=dict(res.get(name,{}),weighted_log_loss=float(loss[name].sum()/sum(q)),weighted_delta_vs_tied5=float(delta.sum()/sum(q)),weighted_delta_se=float(delta.std(ddof=1)*np.sqrt(80)/sum(q)))
(ROOT/'shift_diagnostics.json').write_text(json.dumps(res,indent=2));print(json.dumps(res,indent=2))
