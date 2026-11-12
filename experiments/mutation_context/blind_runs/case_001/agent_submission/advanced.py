import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='4'
import numpy as np,json
from scipy.optimize import minimize
from scipy.stats import chi2
from model import load,features,counts
from explore import ROOT

def motif_features(ctx):
    f,_,rev=features(ctx)
    eq=((ctx//1024)%4)==((ctx//4)%4)
    group=np.where(f['cpg']==0,eq.astype(int),f['cpg']+1)
    return group.astype(int),eq,f,rev

def tied_fit(ct,variant='tied5'):
    # Outcome columns in XOR order: unchanged, non-complement transversion,
    # transition, complement transversion.
    def make(v):
        q=np.exp(v)
        if variant=='tied5': tv,ts,av,ats,cpts=q; rates=np.array([[tv,ts,tv],[av,ats,tv],[tv,ts,tv],[tv,cpts,tv]])
        elif variant=='tied6': tv1,ts,tv2,av,ats,cpts=q;rates=np.array([[tv1,ts,tv2],[av,ats,tv2],[tv1,ts,tv2],[tv1,cpts,tv2]])
        else: raise ValueError(variant)
        return np.column_stack([1-rates.sum(1),rates])
    init=np.log([.0008,.0016,.0032,.0064,.0096] if variant=='tied5' else [.0008,.0016,.0008,.0032,.0064,.0096])
    ids=np.array([[0,1,0],[2,3,0],[0,1,0],[0,4,0]] if variant=='tied5' else [[0,1,2],[3,4,2],[0,1,2],[0,5,2]])
    def fun(v):
        p=make(v)
        if np.any(p<=0):return 1e20,np.zeros_like(v)
        derivative=(ct[:,0]/p[:,0])[:,None]-ct[:,1:]/p[:,1:]
        grad=np.bincount(ids.ravel(),weights=(derivative*p[:,1:]).ravel(),minlength=len(v))/ct.sum()
        return -np.sum(ct*np.log(p))/ct.sum(),grad
    fit=minimize(fun,init,jac=True,method='L-BFGS-B',options={'ftol':1e-14,'gtol':1e-10,'maxiter':300})
    if not fit.success:
        # Near a solution, quasi-Newton line searches can hit floating-point
        # precision. Solve the concave multinomial likelihood directly in
        # scaled probability parameters with an exact Hessian.
        design=np.zeros((4,4,len(init)))
        for gi in range(4):
            for oi in range(3):
                design[gi,oi+1,ids[gi,oi]]+=.001
                design[gi,0,ids[gi,oi]]-=.001
        q=np.exp(fit.x)/.001;offset=np.zeros((4,4));offset[:,0]=1
        for iteration in range(100):
            pp=offset+np.einsum('gok,k->go',design,q)
            gradient=-np.einsum('go,gok->k',ct/pp,design)/ct.sum()
            if np.max(np.abs(gradient))<1e-10:return pp
            hessian=np.einsum('go,gok,gol->kl',ct/pp**2,design,design)/ct.sum()
            step=np.linalg.solve(hessian,gradient);factor=1.
            old=-np.sum(ct*np.log(pp))/ct.sum()
            for _ in range(40):
                candidate=q-factor*step
                pp1=offset+np.einsum('gok,k->go',design,candidate)
                if np.all(pp1>0) and -np.sum(ct*np.log(pp1))/ct.sum()<=old+1e-15:
                    q=candidate;break
                factor*=.5
            else:raise RuntimeError('Exact Newton line search failed')
        raise RuntimeError('Exact Newton score did not converge')
    return make(fit.x)

def main():
    ctx,y,c=load();group,eq,f,rev=motif_features(ctx);center=(ctx//64)%4;out=y^center
    f['motif4']=group;f['motif6']=f['cpg']*2+eq;f['motif8']=group*2+rev
    # Raw context models already tested; pooled ones are the parsimonious challengers.
    for d in [1,3]:
        equal=((ctx//4**(3+d))%4)==((ctx//4**(3-d))%4)
        f[f'equal{d}']=group*2+equal
    z=np.where(rev,0,0) # All other features are reverse-complement invariant below.
    for name in ['context3','context5','context7']:
        f[name]=f[name]
    configs=[('motif4',0),('motif6',0),('motif8',0),('equal1',0),('equal3',0),('tied5',0),('tied6',0)]
    for name in ['context5','context7']:
        for strength in [500,2000,10000,50000,200000,1000000]:configs.append((name,strength))
    uniq={k:np.unique(v,return_index=True)[1] for k,v in f.items() if k in ['context5','context7']}
    rng=np.random.default_rng(0);perm=rng.permutation(80);fold=np.empty(80,int);fold[perm]=np.arange(80)%5
    score={str(z):np.zeros(80) for z in configs};brier={str(z):0. for z in configs};mutation={str(z):[0.,0.] for z in configs}
    group_validation=[]
    for fi in range(5):
        train=fold[c]!=fi;test=~train
        ct=counts(group[train],out[train],4);pr=(ct+.5)/(ct.sum(1)[:,None]+2)
        for name,strength in configs:
            if name.startswith('tied'): p=tied_fit(ct,name);pred=p[group[test]]
            else:
                sz=int(f[name].max()+1);cc=counts(f[name][train],out[train],sz)
                if strength:
                    pp=np.empty_like(cc);pp[:]=.25;ix=uniq[name];pp[f[name][ix]]=pr[group[ix]]
                    p=(cc+strength*pp)/(cc.sum(1)[:,None]+strength)
                else:p=(cc+.5)/(cc.sum(1)[:,None]+2)
                pred=p[f[name][test]]
            loss=-np.log(pred[np.arange(test.sum()),out[test]])
            key=str((name,strength));score[key]+=np.bincount(c[test],weights=loss,minlength=80)
            brier[key]+=float(np.sum(np.sum(pred**2,axis=1)-2*pred[np.arange(test.sum()),out[test]]+1))
            mutation[key][0]+=float(np.sum(1-pred[:,0]));mutation[key][1]+=int(np.sum(out[test]!=0))
            if name=='tied5':
                for g in range(4):
                    gg=group[test]==g
                    group_validation.append({'fold':fi,'group':g,'n':int(gg.sum()),'observed_changes':int(np.sum(out[test][gg]!=0)),'expected_changes':float(np.sum(1-pred[gg,0]))})
        print('completed advanced fold',fi,flush=True)
    results=[];ref=score[str(('tied5',0))]
    for name,strength in configs:
        key=str((name,strength));delta=(score[key]-ref)/49994
        results.append({'model':name,'strength':strength,'mean_log_loss':float(score[key].sum()/len(y)),'brier':brier[key]/len(y),'delta_log_loss_vs_tied5':float(delta.mean()),'delta_se_contig':float(delta.std(ddof=1)/np.sqrt(80)),'expected_and_observed_changes':mutation[key]})
    results.sort(key=lambda z:z['mean_log_loss'])
    ct=counts(group,out,4);fit5=tied_fit(ct,'tied5');fit6=tied_fit(ct,'tied6');full=ct/ct.sum(1)[:,None]
    lr=2*np.sum(ct*np.log(full/fit5)); lr6=2*np.sum(ct*np.log(fit6/fit5))
    output={'cv':results,'counts_xor':ct.tolist(),'fit_tied5_xor':fit5.tolist(),'fit_tied6_xor':fit6.tolist(),'fit_motif4_xor':full.tolist(),'tied5_vs_full_LR':float(lr),'tied5_vs_full_p':float(chi2.sf(lr,7)),'tied5_vs_tied6_LR':float(lr6),'tied5_vs_tied6_p':float(chi2.sf(lr6,1)),'oof_group_validation':group_validation}
    (ROOT/'advanced_results.json').write_text(json.dumps(output,indent=2));print(json.dumps(output,indent=2))

if __name__=='__main__':main()
