from beta_model import *
def fit_np(c,steps=3000):
    n=c.sum(axis=1)
    pairs=np.concatenate([np.c_[n,c[:,j]] for j in [1,2,3]])
    u,mult=np.unique(pairs,axis=0,return_counts=True)
    grid=np.r_[0,np.geomspace(.00001,.3,99)]
    score=u[:,1,None]*np.log(np.maximum(grid,1e-200))+(u[:,0]-u[:,1])[:,None]*np.log1p(-grid)
    weights=np.ones(len(grid))/len(grid);old=-np.inf
    for i in range(steps):
        z=score+np.log(np.maximum(weights,1e-200));l=logsumexp(z,axis=1);resp=np.exp(z-l[:,None]);weights=(resp.T@mult)/mult.sum();total=float(l@mult)
        if total-old<1e-7:break
        old=total
    return grid,weights,{'iterations':i+1,'grid':grid.tolist(),'weights':weights.tolist(),'log_likelihood':total}
def predict_np(c,grid,w):
    n=c.sum(axis=1);p=np.empty(c.shape,float)
    for j in [1,2,3]:
        u,inv=np.unique(np.c_[n,c[:,j]],axis=0,return_inverse=True)
        score=u[:,1,None]*np.log(np.maximum(grid,1e-200))+(u[:,0]-u[:,1])[:,None]*np.log1p(-grid)+np.log(np.maximum(w,1e-200))
        post=np.exp(score-logsumexp(score,axis=1)[:,None]);p[:,j]=(post@grid)[inv]
    p[:,0]=1-p[:,1:].sum(axis=1);assert np.all(p>0)
    return p
def get10():
    cc=np.load(ROOT/'counts11.npy',mmap_mode='r').reshape(3,4**10,4,4).sum(axis=2)
    mid=np.arange(4**10)//4**4%4
    return np.take_along_axis(cc,(mid[:,None]+np.arange(4))[None,:,:]%4,axis=2)
def main():
    c=get10();out=[]
    g,w,info=fit_np(c[0]);q=predict_np(c[0],g,w)
    row={'model':'independent_nonparametric_binomial','tuning_ll':float(np.sum(c[1]*np.log(q))),'assessment_ll':float(np.sum(c[2]*np.log(q))),**info};out.append(row);print(json.dumps(row),flush=True)
    for h in [1,2,3,5]:
        w,ps,info=fit(c[0],h,steps=2000);q=predict(c[0],w,ps)
        row={'model':f'mixture_{h}','tuning_ll':float(np.sum(c[1]*np.log(q))),'assessment_ll':float(np.sum(c[2]*np.log(q))),**info};out.append(row);print(json.dumps(row),flush=True)
    (ROOT/'refined_validation.json').write_text(json.dumps(out,indent=2))
if __name__=='__main__':main()
