from analyze import *
from scipy.special import logsumexp
def relative(c,k):
    mid=np.arange(4**k)//4**(k//2)%4
    return np.take_along_axis(c,(mid[:,None]+np.arange(4))%4,axis=1)
def absolute(p,k):
    mid=np.arange(4**k)//4**(k//2)%4
    return np.take_along_axis(p,(np.arange(4)-mid[:,None])%4,axis=1)
def fit(c,h=1,steps=100):
    c,multiplicity=np.unique(c,axis=0,return_counts=True)
    weights=np.r_[.98,np.full(3*h,.02/(3*h))];rates=np.geomspace(.025,.08,h)
    ps=np.empty((1+3*h,4));ps[0]=[.996,.004/3,.004/3,.004/3]
    for j,r in enumerate(rates):
        for a in range(3):ps[1+3*j+a]=[1-r-.002,.001,.001,.001];ps[1+3*j+a,a+1]=r
    old=-np.inf
    for i in range(steps):
        weights=np.maximum(weights,1e-100);weights/=weights.sum()
        score=c@np.log(ps).T+np.log(weights);lse=logsumexp(score,axis=1);resp=np.exp(score-lse[:,None])*multiplicity[:,None];total=lse@multiplicity
        wc=resp.T@c+.01
        weights[0]=resp[:,0].sum()/multiplicity.sum()
        ps[0,0]=wc[0,0]/wc[0].sum();ps[0,1:]=(1-ps[0,0])/3
        for j in range(h):
            sl=slice(1+3*j,4+3*j);weights[sl]=resp[:,sl].sum()/multiplicity.sum()/3
            nc=wc[sl];den=nc.sum();ref=nc[:,0].sum()/den;target=sum(nc[a,a+1] for a in range(3))/den;other=(1-ref-target)/2
            for a in range(3):ps[1+3*j+a]=[ref,other,other,other];ps[1+3*j+a,a+1]=target
        if total-old<1e-6:break
        old=total
    weights=np.maximum(weights,1e-100);weights/=weights.sum()
    return weights,ps,{'iterations':i+1,'weights':weights.tolist(),'components':ps.tolist(),'train_marginal_log_likelihood':float(total)}
def predict(c,w,ps):
    score=c@np.log(ps).T+np.log(w);resp=np.exp(score-logsumexp(score,axis=1)[:,None]);return resp@ps
def rc_map(k):
    idx=np.arange(4**k);v=idx.copy();rc=np.zeros_like(idx)
    for _ in range(k):rc=4*rc+3-v%4;v//=4
    return rc
def main():
    d=np.load(ROOT/'context_counts.npz');out=[]
    for k in [5,7,9]:
        cc=d[f'k{k}'].astype(float);train=relative(cc[0],k);val=relative(cc[1],k);test=relative(cc[2],k)
        for rc in [False,True]:
            if rc:
                mapping=rc_map(k);train2=train+train[mapping][:,[0,3,2,1]]
                canonical=np.arange(len(train))<mapping;fitting=train2[canonical]
            else:train2=train;fitting=train
            for h in [1,2,3]:
                w,ps,info=fit(fitting,h);pred=predict(train2,w,ps)
                row={'k':k,'rc':rc,'hot_classes':h,'tuning_ll':float(np.sum(val*np.log(pred))),'assessment_ll':float(np.sum(test*np.log(pred))),**info};out.append(row);print(json.dumps(row),flush=True)
    (ROOT/'mixture_validation.json').write_text(json.dumps(out,indent=2))
if __name__=='__main__':main()
