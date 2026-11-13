from destination_models import *
from scipy.special import roots_hermitenorm, expit
from scipy.optimize import minimize
def logitnormal_fit(n,m):
    nn,mm,freq=sufficient(n,m);nodes,weights=roots_hermitenorm(80);weights/=np.sqrt(2*np.pi)
    def obj(v):
        p=np.clip(expit(v[0]+np.exp(v[1])*nodes),1e-12,1-1e-12)
        return -float(logsumexp(mm[:,None]*np.log(p)+(nn-mm)[:,None]*np.log1p(-p)+np.log(weights),axis=1)@freq)
    opt=minimize(obj,[-6,0],method='Nelder-Mead',options={'xatol':1e-9,'maxiter':2000});p=expit(opt.x[0]+np.exp(opt.x[1])*nodes)
    return {'probabilities':p.tolist(),'weights':weights.tolist(),'mu':float(opt.x[0]),'sigma':float(np.exp(opt.x[1])),'loglik':-float(opt.fun)}
def run():
    _,xs,ys=inputs();s=json.loads((ROOT/'splits.json').read_text());a=s['discovery'];b=s['validation'];n,m=get_counts([xs[i]for i in a],[ys[i]for i in a],-4,8);vn,vm=get_counts([xs[i]for i in b],[ys[i]for i in b],-4,8);p=m.sum()/n.sum();base=score(vn,vm,np.full(len(n),p));fit=logitnormal_fit(n,m);gain=score(vn,vm,mix_rates(n,m,fit))-base
    ids=np.arange(4**8);res=[]
    for removed in range(8):
        z=np.zeros(len(ids),dtype=int)
        for j in range(8):
            if j!=removed:z=z*4+(ids//4**(7-j)%4)
        nn=np.bincount(z,weights=n);mm=np.bincount(z,weights=m);vnn=np.bincount(z,weights=vn);vmm=np.bincount(z,weights=vm);f=beta_fit(nn,mm);rate=(mm+f['a'])/(nn+f['a']+f['b']);res.append({'removed_offset':removed-4,'validation_gain':score(vnn,vmm,rate)-base})
    out={'logitnormal_fit':fit,'logitnormal_validation_gain':gain,'leave_one_offset_out':res};(ROOT/'structure_results.json').write_text(json.dumps(out,indent=2));print(json.dumps({'logitnormal_validation_gain':gain,'mu':fit['mu'],'sigma':fit['sigma'],'leave_one_offset_out':res},indent=2))
if __name__=='__main__':run()
