from analysis import *
p,o,n=load(); sp=json.loads((BASE/'split.json').read_text()); ss=np.repeat([sp[k] for k in p],50000)
x=np.concatenate(list(p.values())).astype(int);y=np.concatenate(list(o.values())).astype(int);m=x!=y
base=m[ss==0].mean(); results=[]
for k in range(1,8):
 for left in range(-k+1,1):
  right=left+k-1; K=4**k
  c=np.concatenate([sum(np.roll(v,-j).astype(np.int32)*4**(right-j) for j in range(left,right+1)) for v in p.values()]); a=ss==0;b=ss==1
  ct=np.bincount(c[a],minlength=K);ev=np.bincount(c[a],weights=m[a],minlength=K)
  alt=np.bincount(c[a&m]*4+y[a&m],minlength=K*4).reshape(K,4)
  for alpha in [10,30,100,300,1000,3000,10000,30000,100000]:
   rate=(ev+alpha*base)/(ct+alpha); pr=rate[c[b]]; truth=m[b]
   loss=-np.mean(truth*np.log(pr)+(1-truth)*np.log1p(-pr))
   # Conditional alt prior uniform over nonparent possibilities.
   parents=(np.arange(K)//4**right)%4
   prior=np.ones((K,4))/3; prior[np.arange(K),parents]=0
   for beta in [0.3,3,10,30,100,1000,1e9]:
    dist=(alt+beta*prior)/(ev[:,None]+beta)
    full=loss-np.log(dist[c[b&m],y[b&m]]).sum()/b.sum()
    results.append(dict(left=left,right=right,alpha=alpha,beta=beta,dev_binary_loss=loss,dev_full_loss=full))
  if left==-3 and right==2:
   print('sixmer event quantiles',np.quantile(ev,[0,.1,.25,.5,.75,.9,.95,.99,1])); print('sixmer hist',np.bincount(ev.astype(int)).tolist())
   inds=np.argsort(ev)[-25:][::-1]
   print('hot sixmers',[dict(code=int(i),count=int(ct[i]),changes=int(ev[i]),alt=alt[i].tolist()) for i in inds])
results.sort(key=lambda z:z['dev_full_loss']);print(json.dumps(results[:30],indent=2));(BASE/'window_validation.json').write_text(json.dumps(results,indent=2))
print('BEST WINDOWS')
for left,right in sorted(set((r['left'],r['right']) for r in results)):
 print(next(r for r in results if (r['left'],r['right'])==(left,right)))
