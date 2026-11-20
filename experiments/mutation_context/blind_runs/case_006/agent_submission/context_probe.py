from analyze import *
from scipy.stats import pearsonr
def main():
    d=np.load(ROOT/'context_counts.npz');report={}
    for k in [3,5,7,9]:
        a=d[f'k{k}'];mid=(np.arange(4**k)//4**(k//2))%4
        n=a.sum(axis=2);mut=n-a[:,np.arange(4**k),mid]
        x=mut[0]/n[0];z=(mut[1]+mut[2])/(n[1]+n[2])
        print(k,'rate corr:',pearsonr(x,z),'mutcounts',np.quantile(mut.sum(axis=0),[0,.5,.9,.99,.999,1]))
        total=a.sum(axis=0);num=total.sum(axis=1);m=num-total[np.arange(4**k),mid]
        top=np.argsort(m/num)[-30:][::-1]
        print('TOP', [(''.join(ALPH[(int(i)//4**j)%4] for j in range(k-1,-1,-1)),int(num[i]),int(m[i]),total[i].tolist()) for i in top[:15]])
        tr=a[0];val=a[1];rate=(tr.sum()-np.trace(d['k1'][0]))/tr.sum() if False else 192005/38400000
        p=np.full(tr.shape,rate/3);p[np.arange(4**k),mid]=1-rate
        for alpha in [100,300,1000,3000,10000]:
            q=(tr+alpha*p)/(tr.sum(axis=1)[:,None]+alpha)
            rateq=1-q[np.arange(4**k),mid]
            change_ll=float(np.sum((val.sum(axis=1)-val[np.arange(4**k),mid])*np.log(rateq/rate)+val[np.arange(4**k),mid]*np.log((1-rateq)/(1-rate))))
            total_ll=float(np.sum(val*np.log(q/p)))
            print(k,alpha,'gain rate',change_ll,'gain allele',total_ll-change_ll)
    # Base-independent and reverse-complement context pooling, generalizing variants.
if __name__=='__main__':main()
