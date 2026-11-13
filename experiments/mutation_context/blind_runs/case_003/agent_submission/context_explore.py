from analysis import *
from scipy.stats import chi2
p,o,n=load(); keys=list(p); order=np.random.default_rng(0).permutation(len(keys)); split=np.empty(len(keys),dtype=int); split[order[:48]]=0;split[order[48:64]]=1;split[order[64:]]=2
(BASE/'split.json').write_text(json.dumps({k:int(s) for k,s in zip(keys,split)},indent=2))
x=np.concatenate(list(p.values())); y=np.concatenate(list(o.values())); change=x!=y; ss=np.repeat(split,50000)
results=[]
def codes(left,right):
 return np.concatenate([sum(np.roll(v,-j).astype(np.int32)*4**(right-j) for j in range(left,right+1)) for v in p.values()])
def desc(code,k):
 s=''
 for _ in range(k):s='ACGT'[code%4]+s;code//=4
 return s
for radius in range(5):
 c=codes(-radius,radius);k=2*radius+1;K=4**k
 count=np.bincount(c,minlength=K); events=np.bincount(c,weights=change,minlength=K)
 top=np.argsort((events-count*change.mean())/np.sqrt(count*change.mean()+1))[-20:][::-1]
 print('RADIUS',radius,'top',[(desc(int(i),k),int(count[i]),int(events[i]),round(events[i]/max(count[i],1),4)) for i in top])
 ct=np.bincount(c[ss==0],minlength=K); ev=np.bincount(c[ss==0],weights=change[ss==0],minlength=K);base=change[ss==0].mean()
 for alpha in [1,10,100,1000,10000,100000,1000000]:
  rate=(ev+alpha*base)/(ct+alpha)
  v=ss==1; pr=rate[c[v]]; truth=change[v]
  ll=-np.mean(truth*np.log(pr)+(1-truth)*np.log1p(-pr))
  results.append(dict(radius=radius,alpha=alpha,dev_logloss=ll,predicted_changes=float(pr.sum()),observed_changes=int(truth.sum())))
print(json.dumps(results,indent=2));(BASE/'context_validation.json').write_text(json.dumps(results,indent=2))
# spatial and phase checks
for mod in [2,3,4,5,10,100]:
 idx=np.tile(np.arange(50000)%mod,80);ct=np.bincount(idx);ev=np.bincount(idx,weights=change);exp=ct*change.mean();stat=np.sum((ev-exp)**2/exp)
 print('phase',mod,'p',chi2.sf(stat,mod-1),'rates range',(ev/ct).min(),(ev/ct).max())
print('scaffold chi2',chi2.sf(sum(((v!=o[k]).sum()-change.mean()*len(v))**2/(change.mean()*len(v)) for k,v in p.items()),79))
