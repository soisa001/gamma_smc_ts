from explore import *
from window_scan import window,get_counts,fit_mix,mix_rates
from context_models import score
from scipy.special import betaln, logsumexp
from scipy.optimize import minimize

def sufficient(n,m):
    nm,freq=np.unique(np.column_stack([n,m]),axis=0,return_counts=True)
    return *nm.T,freq
def beta_fit(n,m):
    nn,mm,freq=sufficient(n,m)
    def obj(t):
        a,b=np.exp(t);return -float((betaln(mm+a,nn-mm+b)-betaln(a,b))@freq)
    fit=minimize(obj,np.log([1,199]),method='Nelder-Mead',options={'maxiter':2000,'xatol':1e-8,'fatol':1e-7})
    return {'a':float(np.exp(fit.x[0])),'b':float(np.exp(fit.x[1])),'loglik':-float(fit.fun),'success':bool(fit.success)}
def npmix_fit(n,m):
    nn,mm,freq=sufficient(n,m);p=np.r_[1e-8,np.geomspace(1e-5,.3,90)];w=np.full(len(p),1/len(p))
    kernel=mm[:,None]*np.log(p)+(nn-mm)[:,None]*np.log1p(-p)
    for i in range(5000):
        lp=kernel+np.log(w+1e-300);norm=logsumexp(lp,axis=1);resp=np.exp(lp-norm[:,None])*freq[:,None];w=resp.sum(0)/freq.sum()
    return {'probabilities':p.tolist(),'weights':w.tolist(),'loglik':float(norm@freq),'iterations':i+1}
def dest_counts(xs,ys):
    d=np.zeros((4**8,4))
    for x,y in zip(xs,ys):
        z=window(x,-4,8);m=x!=y
        np.add.at(d,(z[m],y[m]),1)
    return d
def run():
    _,xs,ys=inputs();s=json.loads((ROOT/'splits.json').read_text());ix=s['discovery'];vi=s['validation'];xx=[xs[i]for i in ix];yy=[ys[i]for i in ix];vx=[xs[i]for i in vi];vy=[ys[i]for i in vi]
    n,m=get_counts(xx,yy,-4,8);vn,vm=get_counts(vx,vy,-4,8);p=m.sum()/n.sum();base=score(vn,vm,np.full(len(n),p))
    beta=beta_fit(n,m);mix=fit_mix(n,m);npmix=npmix_fit(n,m)
    fits={'beta':beta,'two_point_mixture':mix,'nonparametric_mixture':npmix}
    pred={'uniform':np.full(len(n),p),'beta':(m+beta['a'])/(n+beta['a']+beta['b']),'two_point_mixture':mix_rates(n,m,mix),'nonparametric_mixture':mix_rates(n,m,npmix)}
    out={'fits':fits,'rate_models':{k:{'validation_binary_nll':-score(vn,vm,v)/vn.sum(),'validation_gain':score(vn,vm,v)-base}for k,v in pred.items()}}
    d=dest_counts(xx,yy);vd=dest_counts(vx,vy);central=np.arange(4**8)//4**3%4
    dest_rows=[]
    for tau in [.1,.3,1,3,10,30,100,300,1000,1e12]:
        q=(d+tau/3)/(m[:,None]+tau);q[np.arange(len(q)),central]=0
        dest_rows.append({'tau':tau,'validation_gain_vs_uniform':float(np.sum(vd*np.log(np.maximum(q,1e-300)))+vd.sum()*np.log(3))})
    out['destination_models']=dest_rows
    # Corroborating destination sharing: among pairs of observed changes with the same 8-mer.
    pairs=float(np.sum(m*(m-1)/2));same=float(np.sum(d*(d-1)/2));out['same_destination_mutation_pairs']={'total':pairs,'same':same,'fraction':same/pairs}
    (ROOT/'refinement_results.json').write_text(json.dumps(out,indent=2));np.savez(ROOT/'discovery_candidate_rates.npz',**pred)
    display={k:v for k,v in out.items()if k!='fits'};display['beta_fit']=beta;display['mixture']=mix;display['np_mixture_top_components']=sorted(zip(npmix['probabilities'],npmix['weights']),key=lambda a:a[1],reverse=True)[:10];print(json.dumps(display,indent=2))
if __name__=='__main__':run()
