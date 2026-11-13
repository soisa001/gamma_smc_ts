from analysis import *
from scipy.special import logsumexp, gammaln, expit
from scipy.optimize import minimize

def context_codes(v,left=-3,right=2):
 return sum(np.roll(v,-j).astype(np.int32)*4**(right-j) for j in range(left,right+1))
def sufficient(p,o,keys):
 ct=np.zeros(4096);alt=np.zeros((4096,4))
 for k in keys:
  c=context_codes(p[k]); y=o[k];m=p[k]!=y
  ct+=np.bincount(c,minlength=4096);alt+=np.bincount(c[m]*4+y[m],minlength=16384).reshape(4096,4)
 return ct,alt
PARENTS=(np.arange(4096)//16)%4
OTHER=np.array([[a for a in range(4) if a!=p] for p in PARENTS])
def unpack(alt):return alt[np.arange(4096)[:,None],OTHER]
def score(rate,dist,ct,alt):
 e=alt.sum(1); return float(-(np.sum((ct-e)*np.log1p(-rate))+np.sum(e*np.log(rate))+np.sum(unpack(alt)*np.log(dist)))/ct.sum())
def fit_beta(ct,e):
 def f(z):
  a,b=np.exp(z);return -np.sum(gammaln(e+a)+gammaln(ct-e+b)-gammaln(ct+a+b)-gammaln(a)-gammaln(b)+gammaln(a+b))
 opt=minimize(f,np.log([1,200]),method='Nelder-Mead',options={'maxiter':500});a,b=np.exp(opt.x)
 return (e+a)/(ct+a+b),dict(a=a,b=b,prior_mean=a/(a+b),success=bool(opt.success))
def fit_grid(ct,e,kind='logitnormal'):
 z=np.linspace(-15,0,500);rate=expit(z);ll=e[:,None]*np.log(rate)+(ct-e)[:,None]*np.log1p(-rate)
 def f(par):
  mu,ls=par; lp=-.5*((z-mu)/np.exp(ls))**2-ls;lp-=logsumexp(lp)
  return -logsumexp(ll+lp,axis=1).sum()
 opt=minimize(f,[-5.8,.1],method='Nelder-Mead',options={'maxiter':500});mu,ls=opt.x
 w=np.exp(-.5*((z-mu)/np.exp(ls))**2);w/=w.sum()
 if kind=='npmle':
  for j in range(300):
   post=np.exp(ll+np.log(np.maximum(w,1e-300))-logsumexp(ll+np.log(np.maximum(w,1e-300)),axis=1)[:,None]);wn=post.mean(0)
   if abs(wn-w).sum()<1e-7:break
   w=wn
 post=np.exp(ll+np.log(np.maximum(w,1e-300))-logsumexp(ll+np.log(np.maximum(w,1e-300)),axis=1)[:,None]);pred=post@rate
 return pred,dict(mu=float(mu),sd=float(np.exp(ls)),prior_mean=float(w@rate),success=bool(opt.success)),(z,w)
def fit_alt(alt,kind='dirichlet'):
 a=unpack(alt);e=a.sum(1)
 def marginal(v):return gammaln(3*v)-gammaln(e+3*v)+np.sum(gammaln(a+v)-gammaln(v),axis=1)
 opt=minimize(lambda z:-marginal(np.exp(z[0])).sum(),[0.],method='Nelder-Mead');v=float(np.exp(opt.x[0]))
 if kind=='dirichlet':return (a+v)/(e[:,None]+3*v),dict(alpha=v)
 vals=np.logspace(-2,2,35);ll=np.stack([marginal(vv) for vv in vals],axis=1);w=np.ones(len(vals))/len(vals)
 for _ in range(300):
  post=np.exp(ll+np.log(w)-logsumexp(ll+np.log(w),axis=1)[:,None]);w=post.mean(0)
 dist=np.sum(post[:,:,None]*(a[:,None,:]+vals[None,:,None])/(e[:,None,None]+3*vals[None,:,None]),axis=1)
 return dist,dict(alpha_grid=vals.tolist(),weight=w.tolist())
if __name__=='__main__':
 p,o,n=load();sp=json.loads((BASE/'split.json').read_text());ct,alt=sufficient(p,o,[k for k in p if sp[k]==0]);vc,va=sufficient(p,o,[k for k in p if sp[k]==1]);e=alt.sum(1)
 rates={'fixed300':((e+300*e.sum()/ct.sum())/(ct+300),{})}
 rates['beta']=fit_beta(ct,e);rates['logitnormal']=fit_grid(ct,e)[:2];rates['npmle']=fit_grid(ct,e,'npmle')[:2]
 dists={'fixed3':((unpack(alt)+1)/(e[:,None]+3),{})}
 for kind in ['dirichlet','mixture']:dists[kind]=fit_alt(alt,kind)
 res=[]
 for rn,(r,rp) in rates.items():
  for dn,(d,dp) in dists.items():res.append(dict(rate=rn,destination=dn,dev_loss=score(r,d,vc,va),rate_parameters=rp,destination_parameters=dp))
 res.sort(key=lambda a:a['dev_loss']);(BASE/'bayes_validation.json').write_text(json.dumps(res,indent=2));print(json.dumps(res,indent=2))
 np.savez(BASE/'development_counts.npz',count=ct,alt=alt,validation_count=vc,validation_alt=va)
