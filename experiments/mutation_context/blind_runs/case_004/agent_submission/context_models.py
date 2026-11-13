from diagnose import context
from explore import *
from scipy.stats import chi2, binomtest
def counts(xs,ys,r):
    K=4**(2*r+1);ns=np.zeros(K);ms=np.zeros(K);dest=np.zeros((K,4))
    for x,y in zip(xs,ys):
        z=context(x,r);m=x!=y; ns+=np.bincount(z,minlength=K);ms+=np.bincount(z,weights=m,minlength=K)
        for b in range(4):dest[:,b]+=np.bincount(z,weights=m&(y==b),minlength=K)
    return ns,ms,dest
def score(ns,ms,rates):return float(np.sum(ms*np.log(rates)+(ns-ms)*np.log1p(-rates)))
def text_context(z,r):
    s=''
    for _ in range(2*r+1):s=BASES[z%4]+s;z//=4
    return s
def run():
    names,xs,ys=inputs();s=json.loads((ROOT/'splits.json').read_text());idx=s['discovery'];vidx=s['validation'];p=sum((xs[i]!=ys[i]).sum() for i in idx)/sum(len(xs[i])for i in idx)
    res=[];top={}
    for r in range(5):
        n,m,d=counts([xs[i]for i in idx],[ys[i]for i in idx],r);vn,vm,vd=counts([xs[i]for i in vidx],[ys[i]for i in vidx],r)
        for a in [50,200,1000,5000,20000,100000,1e6,1e12]:
            rate=(m+a*p)/(n+a)
            res.append({'radius':r,'prior_strength':a,'validation_nll':-score(vn,vm,rate)/vn.sum(),'validation_gain_vs_uniform':score(vn,vm,rate)-score(vn,vm,np.full(len(n),p))})
        rank=np.argsort((m-n*p)/np.sqrt(n*p+1e-10))[-30:][::-1]
        top[str(r)]=[{'context':text_context(int(k),r),'n':int(n[k]),'mutations':int(m[k]),'validation_n':int(vn[k]),'validation_mutations':int(vm[k]),'p_raw':float(binomtest(int(m[k]),int(n[k]),p,alternative='greater').pvalue)}for k in rank]
    out={'discovery_rate':p,'validation_models':res,'top_contexts':top}
    (ROOT/'context_results.json').write_text(json.dumps(out,indent=2));print(json.dumps({'best':sorted(res,key=lambda z:z['validation_nll'])[:15],'top_7mer':top['3'],'top_9mer':top['4'][:10]},indent=2))
if __name__=='__main__':run()
