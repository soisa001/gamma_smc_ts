from analyze import *
from scipy.stats import chi2, pearsonr, norm
def model_uniform(c,k=1):
    mid=(np.arange(4**k)//(4**(k//2)))%4
    r=1-np.trace(c)/c.sum()
    probs=np.full((4**k,4),r/3);probs[np.arange(4**k),mid]=1-r
    return probs
def ll(c,p):return float(np.sum(c*np.log(p)))
def main():
    data=np.load(ROOT/'context_counts.npz'); counts={k:data[f'k{k}'] for k in [1,3,5,7,9]}
    c=counts[1][0];parent=np.arange(4)
    uniform=model_uniform(c)
    models={'pooled_uniform':uniform,'parent_specific_rate':None,'transition_transversion':None,'full_substitution':(c+.5)/(c.sum(axis=1)[:,None]+2)}
    r=(c.sum(axis=1)-np.diag(c))/c.sum(axis=1)
    p=np.tile(r[:,None]/3,(1,4));p[parent,parent]=1-r;models['parent_specific_rate']=p
    rate=1-np.trace(c)/c.sum();transitions=np.array([2,3,0,1]);ts=c[parent,transitions].sum()/(c.sum()-np.trace(c))
    p=np.full((4,4),rate*(1-ts)/2);p[parent,parent]=1-rate;p[parent,transitions]=rate*ts;models['transition_transversion']=p
    results=[]
    for name,p in models.items():
        results.append({'model':name,'tuning_log_likelihood':ll(counts[1][1],p),'test_log_likelihood':ll(counts[1][2],p)})
    for k in [3,5,7,9]:
        prior=model_uniform(c,k);n=counts[k][0].sum(axis=1)[:,None]
        for a in [10,100,1000,10000,100000,1000000]:
            p=(counts[k][0]+a*prior)/(n+a)
            results.append({'model':f'context_{k}_alpha_{a}','tuning_log_likelihood':ll(counts[k][1],p),'test_log_likelihood':ll(counts[k][2],p)})
    base=results[0]
    for x in results:
        x['tuning_improvement']=x['tuning_log_likelihood']-base['tuning_log_likelihood'];x['test_improvement']=x['test_log_likelihood']-base['test_log_likelihood']
    exp=json.loads((ROOT/'exploration.json').read_text()); sc=exp['scaffolds'];nr=np.array([s['changes'] for s in sc]);r=nr.sum()/exp['total_sites'];mu=50000*r
    hetero=float(np.sum((nr-mu)**2/(mu*(1-r))))
    positional=data['poscounts'];st,pv,_,_=chi2_contingency(positional)
    offp=[chi2.sf(x['chi2'],12) for x in exp['offsets'] if x['offset']!=0]
    summary={'comparisons':results,'split_sites':counts[1].sum(axis=(1,2)).tolist(),'split_changes':[(a.sum()-np.trace(a)).item() for a in counts[1]],'scaffold_heterogeneity':{'chi2':hetero,'df':len(sc)-1,'p':float(chi2.sf(hetero,len(sc)-1))},'scaffold_gc_rate_correlation':list(map(float,pearsonr([s['gc'] for s in sc],nr))),'relative_position_100bins':{'chi2':float(st),'p':float(pv)},'minimum_offset_p':float(min(offp)),'bonferroni_offset_p':float(min(1,60*min(offp))),'rate':r,'rate_wald_95_interval':[float(r-1.96*np.sqrt(r*(1-r)/exp['total_sites'])),float(r+1.96*np.sqrt(r*(1-r)/exp['total_sites']))]}
    print(json.dumps(summary,indent=2));(ROOT/'validation.json').write_text(json.dumps(summary,indent=2))
if __name__=='__main__':main()
