from explore import *
from context_models import score
from scipy.special import expit, logsumexp

def window(x,start,L):
    z=np.zeros(len(x),dtype=np.int64)
    for k in range(start,start+L):z=z*4+np.roll(x,-k)
    return z
def get_counts(xs,ys,start,L):
    K=4**L;n=np.zeros(K);m=np.zeros(K)
    for x,y in zip(xs,ys):
        z=window(x,start,L);n+=np.bincount(z,minlength=K);m+=np.bincount(z,weights=x!=y,minlength=K)
    return n,m
def fit_mix(n,m):
    active=n>0;nn=n[active];mm=m[active]
    # Group sufficient statistics for efficient EM iterations.
    nm,freq=np.unique(np.column_stack([nn,mm]),axis=0,return_counts=True);nn,mm=nm.T
    p=m.sum()/n.sum();best=None
    for initial_hot in [.02,.08,.3]:
        probs=np.array([p*.6,initial_hot]);weights=np.array([.97,.03]);prev=-np.inf
        for it in range(400):
            lp=mm[:,None]*np.log(probs)+(nn-mm)[:,None]*np.log1p(-probs)+np.log(weights)
            norm=logsumexp(lp,axis=1);ll=float(norm@freq);resp=np.exp(lp-norm[:,None])*freq[:,None]
            weights=np.maximum(resp.sum(0)/freq.sum(),1e-300);weights/=weights.sum()
            numerator=(resp*mm[:,None]).sum(0);denominator=(resp*nn[:,None]).sum(0)
            probs=np.clip(np.divide(numerator,denominator,out=probs.copy(),where=denominator>0),1e-8,1-1e-8)
            if ll-prev<1e-7:break
            prev=ll
        if best is None or ll>best['loglik']:best={'probabilities':probs.tolist(),'weights':weights.tolist(),'loglik':ll,'iterations':it+1}
    return best
def mix_rates(n,m,fit):
    p=np.array(fit['probabilities']);w=np.array(fit['weights']);lp=m[:,None]*np.log(p)+(n-m)[:,None]*np.log1p(-p)+np.log(w)
    return np.exp(lp-logsumexp(lp,axis=1)[:,None])@p
def run():
    _,xs,ys=inputs();s=json.loads((ROOT/'splits.json').read_text());a=s['discovery'];b=s['validation'];xx=[xs[i]for i in a];yy=[ys[i]for i in a];vx=[xs[i]for i in b];vy=[ys[i]for i in b]
    p=sum((x!=y).sum()for x,y in zip(xx,yy))/sum(map(len,xx));res=[]
    for L in range(3,11):
        for start in range(-L+1,1):
            n,m=get_counts(xx,yy,start,L);vn,vm=get_counts(vx,vy,start,L);base=score(vn,vm,np.full(len(n),p))
            fit=fit_mix(n,m);rate=mix_rates(n,m,fit)
            row={'length':L,'start':start,'mixture':fit,'mixture_validation_gain':score(vn,vm,rate)-base}
            candidates=[{'alpha':alpha,'gain':score(vn,vm,(m+alpha*p)/(n+alpha))-base} for alpha in [50,100,200,500,1000,2000,5000,20000,1e12]]
            row['best_beta']=max(candidates,key=lambda z:z['gain']);res.append(row)
        print('L',L,'best',max([z for z in res if z['length']==L],key=lambda z:z['mixture_validation_gain']),flush=True)
        (ROOT/'window_results.json').write_text(json.dumps(res,indent=2))
if __name__=='__main__':run()
