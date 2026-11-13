from explore import *
from scipy.stats import chi2_contingency, chi2, binomtest
from scipy.fft import rfft, irfft

def table_test(z,mut,n):
    count=np.bincount(z,minlength=n); changed=np.bincount(z,weights=mut,minlength=n)
    active=count>0
    t=np.column_stack([changed[active], count[active]-changed[active]])
    stat,p,dof,_=chi2_contingency(t,correction=False)
    return {'chi2':float(stat),'p':float(p),'df':int(dof),'counts':count.tolist(),'changed':changed.tolist()}
def context(x,r):
    z=np.zeros(len(x),dtype=np.int64)
    for k in range(-r,r+1):z=z*4+np.roll(x,-k)
    return z
def run():
    names,xs,ys=inputs(); order=np.random.default_rng(0).permutation(len(xs)); splits={'discovery':order[:40].tolist(),'validation':order[40:60].tolist(),'test':order[60:].tolist()}
    (ROOT/'splits.json').write_text(json.dumps(splits,indent=2))
    ix=splits['discovery']; xx=[xs[i] for i in ix]; yy=[ys[i] for i in ix]
    x=np.concatenate(xx);y=np.concatenate(yy);m=x!=y
    out={'n':len(x),'changed':int(m.sum()),'tests':{}}
    for r in range(4):out['tests'][f'context_radius_{r}']=table_test(np.concatenate([context(a,r) for a in xx]),m,4**(2*r+1))
    pos=np.tile(np.arange(50000),len(ix))
    for mod in [2,3,4,5,7,10,16,20,50,100]:out['tests'][f'position_mod_{mod}']=table_test(pos%mod,m,mod)
    out['tests']['position_1000bp_bins']=table_test(pos//1000,m,50)
    out['tests']['scaffold']=table_test(np.repeat(np.arange(len(ix)),50000),m,len(ix))
    for w in [5,21,101,501,2001]:
        gc=[]
        for a in xx:
            g=np.isin(a,[1,2]).astype(float);g=np.concatenate([g[-w//2:],g,g[:w//2+1]])
            c=np.r_[0,np.cumsum(g)];gc.append((c[w:]-c[:-w])[:len(a)])
        g=np.concatenate(gc);cuts=np.unique(np.quantile(g,np.linspace(0,1,11)))
        out['tests'][f'GC_window_{w}']=table_test(np.searchsorted(cuts[1:-1],g),m,len(cuts)-1)
    # Circular cross-correlations. Record all offsets without cherry-picking their validation.
    L=len(xx[0]); match=np.zeros(L); source=np.zeros(L); comp=np.zeros(L); corr=np.zeros(L)
    for a,b in zip(xx,yy):
        mut=a!=b;fm=rfft(mut.astype(float),workers=4)
        corr+=irfft(np.conj(fm)*fm,n=L,workers=4)
        for base in range(4):
            f=rfft((a==base).astype(float),workers=4)
            match+=irfft(np.conj(rfft((mut&(b==base)).astype(float),workers=4))*f,n=L,workers=4)
            source+=irfft(np.conj(rfft((mut&(a==base)).astype(float),workers=4))*f,n=L,workers=4)
            comp+=irfft(np.conj(rfft((mut&(b==3-base)).astype(float),workers=4))*f,n=L,workers=4)
    def peaks(v,expect,sd):
        z=(v-expect)/sd;order=np.argsort(np.abs(z[1:]))[-15:][::-1]+1
        return [{'offset':int(k if k<L//2 else k-L),'count':float(v[k]),'z':float(z[k])}for k in order]
    out['offset_match_peaks']=peaks(match,m.sum()/4,np.sqrt(m.sum()*.25*.75))
    out['offset_source_peaks']=peaks(source,m.sum()/4,np.sqrt(m.sum()*.25*.75))
    out['offset_complement_match_peaks']=peaks(comp,m.sum()/4,np.sqrt(m.sum()*.25*.75))
    out['mutation_pair_peaks']=peaks(corr,len(m)*m.mean()**2,np.sqrt(len(m)*m.mean()**2))
    np.savez(ROOT/'discovery_offset_scan.npz',match=match,source=source,complement=comp,mutation_pairs=corr)
    (ROOT/'diagnostics.json').write_text(json.dumps(out,indent=2))
    small={k:v for k,v in out.items() if k!='tests'};small['tests']={k:{j:v[j] for j in ['chi2','p','df']} for k,v in out['tests'].items()};print(json.dumps(small,indent=2))
if __name__=='__main__':run()
