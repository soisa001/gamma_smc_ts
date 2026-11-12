from models import *
from scipy.stats import poisson
from itertools import combinations
raw=np.load(ROOT/'all_context_outcomes.npy'); c=np.arange(16384);b=c//64%4;n=raw.sum(1)
pr=fit_tables(raw,False)['cpg'];mu=n*(1-pr[c,b]);obs=n-raw[c,b]; residual=obs-mu
print('mutations',int(obs.sum()),'expected',mu.sum(),'Pearson',np.sum(residual**2/np.maximum(mu,1e-20)))
for k in range(4,13): print('count at least',k,'observed',int((obs>=k).sum()),'expected',float(poisson.sf(k-1,mu).sum()))
for i in np.argsort(poisson.sf(obs-1,mu))[:60]:
    motif=''.join(BASES[i//4**j%4] for j in range(6,-1,-1))
    print(motif,int(n[i]),int(obs[i]),round(mu[i],3),raw[i].tolist(),round(float(-np.log10(poisson.sf(obs[i]-1,mu[i]))),3))
print('aggregate context feature associations conditioned on CpG baseline')
digits=np.array([c//4**j%4 for j in range(6,-1,-1)]).T
res=[]
for k in [1,2,3]:
    for inds in combinations([0,1,2,4,5,6],k):
        ids=b.copy()
        for j in inds: ids=ids*4+digits[:,j]
        nn=4**(k+1);oo=np.bincount(ids,weights=obs,minlength=nn);ee=np.bincount(ids,weights=mu,minlength=nn)
        stat=np.sum((oo-ee)**2/ee);p=chi2.sf(stat,nn-4)
        res.append((float(p),inds,float(stat),nn-4))
for row in sorted(res)[:25]: print(row)
