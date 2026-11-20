from mixture import *
from scipy.special import betaln,gammaln
from scipy.optimize import minimize
def fit_beta(c):
    c,mult=np.unique(c,axis=0,return_counts=True);n=c.sum(axis=1);m=n-c[:,0]
    def objective(params):
        a,b,e=np.exp(params)
        score=betaln(m+a,n-m+b)-betaln(a,b)
        score+=gammaln(3*e)-gammaln(m+3*e)+np.sum(gammaln(c[:,1:]+e)-gammaln(e),axis=1)
        return -float(score@mult)
    opt=minimize(objective,np.log([.2,40,.1]),method='L-BFGS-B',bounds=[(-12,15)]*3,options={'maxiter':500,'ftol':1e-12})
    return np.exp(opt.x),{'success':bool(opt.success),'message':str(opt.message),'parameters':np.exp(opt.x).tolist(),'objective':float(opt.fun)}
def pred_beta(c,params):
    a,b,e=params;n=c.sum(axis=1);m=n-c[:,0];r=(m+a)/(n+a+b)
    return np.c_[1-r,r[:,None]*(c[:,1:]+e)/(m[:,None]+3*e)]
def main():
    c11=np.load(ROOT/'counts11.npy',mmap_mode='r');out=[]
    for left,right in [(4,4),(5,3),(5,4),(5,5)]:
        k=left+right+1;cc=c11.reshape(3,4**(5-left),4**k,4**(5-right),4).sum(axis=(1,3));mid=np.arange(4**k)//4**right%4
        c=np.take_along_axis(cc,(mid[:,None]+np.arange(4))[None,:,:]%4,axis=2)
        par,info=fit_beta(c[0]);q=pred_beta(c[0],par)
        row={'left':left,'right':right,'tuning_ll':float(np.sum(c[1]*np.log(q))),'assessment_ll':float(np.sum(c[2]*np.log(q))),**info};out.append(row);print(json.dumps(row),flush=True)
    (ROOT/'beta_validation.json').write_text(json.dumps(out,indent=2))
if __name__=='__main__':main()
