from refine_model import *
from scipy.special import gammaln
from scipy.optimize import minimize_scalar
from itertools import permutations

def triples(d):
    central=np.arange(len(d))//4**3%4
    return d[np.arange(len(d))[:,None],np.array([[j for j in range(4) if j!=b] for b in central])]
def dirichlet_fit(t):
    tt,freq=np.unique(np.sort(t,axis=1),axis=0,return_counts=True);m=tt.sum(1)
    def obj(logtau):
        tau=np.exp(logtau);return -float((gammaln(tau)-gammaln(m+tau)+(gammaln(tt+tau/3)-gammaln(tau/3)).sum(1))@freq)
    opt=minimize_scalar(obj,bounds=(-5,5),method='bounded')
    return {'tau':float(np.exp(opt.x)),'loglik':-float(opt.fun)}
def preferred_fit(t):
    tt,freq=np.unique(np.sort(t,axis=1),axis=0,return_counts=True)
    def obj(h):
        probs=np.full((3,3),(1-h)/2);np.fill_diagonal(probs,h)
        return -float((logsumexp(tt@np.log(probs).T,axis=1)-np.log(3))@freq)
    opt=minimize_scalar(obj,bounds=(1/3+.00001,.99999),method='bounded');h=float(opt.x)
    probs=np.full((3,3),(1-h)/2);np.fill_diagonal(probs,h)
    return {'probabilities':probs.tolist(),'weights':[1/3]*3,'loglik':-float(opt.fun)}
def qpredict(t,fit):
    p=np.array(fit['probabilities']);w=np.array(fit['weights']);lp=t@np.log(p).T+np.log(w+1e-300)
    return np.exp(lp-logsumexp(lp,axis=1)[:,None])@p
def grid_fit(t):
    # Symmetric distribution over a coarse simplex, then learned weights by shape.
    grid=[]
    for a in range(1,31):
        for b in range(a,31):
            c=33-a-b
            if c>=b:grid.append(sorted(set(permutations(np.array([a,b,c],float)/33))))
    tt,freq=np.unique(np.sort(t,axis=1),axis=0,return_counts=True)
    kernels=np.column_stack([logsumexp(tt@np.log(np.array(g)).T,axis=1)-np.log(len(g))for g in grid]);w=np.full(len(grid),1/len(grid))
    for i in range(2000):
        lp=kernels+np.log(w+1e-300);norm=logsumexp(lp,axis=1);resp=np.exp(lp-norm[:,None])*freq[:,None];w=resp.sum(0)/freq.sum()
    return {'probabilities':[list(p)for g in grid for p in g],'weights':[float(wi/len(g))for wi,g in zip(w,grid)for p in g],'loglik':float(norm@freq)}
def run():
    _,xs,ys=inputs();s=json.loads((ROOT/'splits.json').read_text());a=s['discovery'];b=s['validation'];d=dest_counts([xs[i]for i in a],[ys[i]for i in a]);vd=dest_counts([xs[i]for i in b],[ys[i]for i in b]);t=triples(d);vt=triples(vd)
    fits={'dirichlet':dirichlet_fit(t),'one_preferred':preferred_fit(t),'simplex_grid':grid_fit(t)};res={}
    for key,fit in fits.items():
        if key=='dirichlet':q=(t+fit['tau']/3)/(t.sum(1)[:,None]+fit['tau'])
        else:q=qpredict(t,fit)
        res[key]={'validation_gain':float(np.sum(vt*np.log(q))+vt.sum()*np.log(3)),'loglik':fit['loglik']}
    out={'fits':fits,'scores':res};(ROOT/'destination_results.json').write_text(json.dumps(out,indent=2));print(json.dumps({'dirichlet':fits['dirichlet'],'one_preferred':fits['one_preferred'],'scores':res},indent=2))
if __name__=='__main__':run()
