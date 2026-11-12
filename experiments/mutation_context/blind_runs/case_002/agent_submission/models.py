from analyze import *
from scipy.stats import chi2

def rc_codes(c,k):
    out=np.zeros_like(c)
    for j in range(k): out=out*4+(3-(c//4**j)%4)
    return out

def load_data():
    p=fasta('calibration_parent.fasta'); y=fasta('calibration_offspring.fasta')
    codes=[]; targets=[]; groups=[]
    for i,k in enumerate(p):
        codes.append(context(p[k],3)); targets.append(y[k][3:-3]); groups.append(np.full(len(p[k])-6,i,dtype=np.uint8))
    return np.concatenate(codes),np.concatenate(targets),np.concatenate(groups)

def features(c,symmetric):
    b=(c//64)%4
    if symmetric: c=np.where(b>=2,rc_codes(c,7),c);b=(c//64)%4
    tri=(c//16)%64;five=(c//4)%1024
    cpg=(((b==1)&((tri%4)==2))|((b==2)&((tri//16)==1))).astype(int)
    return {'base':b,'cpg':b*2+cpg,'tri':tri,'five':five,'seven':c}

SIZES={'base':4,'cpg':8,'tri':64,'five':1024,'seven':16384}

def counts(c,y): return np.bincount(c*4+y,minlength=16384*4).reshape(16384,4)

def fit_tables(raw,sym):
    allc=np.arange(16384); feat=features(allc,sym); vals=raw.copy()
    if sym:
        flip=(allc//64)%4>=2;vals[flip]=vals[flip,::-1]
    cts={}
    for key,f in feat.items():
        tab=np.zeros((SIZES[key],4))
        for j in range(4): np.add.at(tab[:,j],f,vals[:,j])
        cts[key]=tab
    tab=cts['base']; base=(tab+.5)/(tab.sum(1)[:,None]+2)
    b8=np.arange(8)//2
    tab=cts['cpg'];cpg=(tab+100*base[b8])/(tab.sum(1)[:,None]+100)
    b64=(np.arange(64)//4)%4
    cpg64=b64*2+(((b64==1)&((np.arange(64)%4)==2))|((b64==2)&((np.arange(64)//16)==1)))
    outputs={'base':base[feat['base']],'cpg':cpg[feat['cpg']]}
    for a in [100,1000,10000,100000,1000000]:
        tab=cts['tri'];tri=(tab+a*cpg[cpg64])/(tab.sum(1)[:,None]+a)
        outputs[f'tri_{a}']=tri[feat['tri']]
    # Higher contexts shrink directly to the low-dimensional CpG model.
    for key in ['five','seven']:
        n=SIZES[key]; k=5 if key=='five' else 7
        inner=(np.arange(n)//4**((k-3)//2))%64
        prior=cpg[cpg64[inner]]
        for a in [100,1000,10000,100000,1000000]:
            tab=cts[key]; pred=(tab+a*prior)/(tab.sum(1)[:,None]+a)
            outputs[f'{key}_{a}']=pred[feat[key]]
    if sym:
        flip=(allc//64)%4>=2
        for pred in outputs.values(): pred[flip]=pred[flip,::-1]
    return outputs

def main():
    c,y,g=load_data(); allraw=counts(c,y)
    rng=np.random.default_rng(0); order=rng.permutation(80); folds=np.zeros(80,dtype=int);folds[order]=np.arange(80)%5
    scores={}; per_contig={}; calibration={}
    for fold in range(5):
        test=folds[g]==fold; train=allraw-counts(c[test],y[test])
        for sym in [False,True]:
            tables=fit_tables(train,sym)
            for name,tab in tables.items():
                name=('symmetric_' if sym else 'oriented_')+name
                pc=tab[c[test]];yt=y[test];gt=g[test]
                losses=-np.log(pc[np.arange(len(pc)),yt])
                b=(c[test]//64)%4
                expected=np.sum(1-pc[np.arange(len(pc)),b])
                score=scores.setdefault(name,{'nll':0.,'expected_changes':0.,'observed_changes':0,'n':0})
                score['nll']+=float(losses.sum());score['expected_changes']+=float(expected);score['observed_changes']+=int(np.sum(yt!=b));score['n']+=len(pc)
                con=per_contig.setdefault(name,{})
                for gi in np.unique(gt): con[int(gi)]=float(losses[gt==gi].sum())
        print('finished fold',fold,flush=True)
    best=min(scores,key=lambda x:scores[x]['nll'])
    for name,v in scores.items():
        v['mean_nll']=v['nll']/v['n'];v['excess_nll_vs_best']=v['nll']-scores[best]['nll']
        delta=np.array([per_contig[name][i]-per_contig[best][i] for i in range(80)])
        v['excess_nll_se']=float(delta.std(ddof=1)*np.sqrt(80))
    result={'seed':0,'split':'5 folds of whole contigs','fold_by_contig_index':folds.tolist(),'scores':scores,'contig_nll':per_contig,'best':best}
    (ROOT/'validation.json').write_text(json.dumps(result,indent=2))
    np.save(ROOT/'all_context_outcomes.npy',allraw)
    print('best',best)
    for name in sorted(scores,key=lambda n:scores[n]['nll']):
        v=scores[name]; print(name,round(v['nll'],3),round(v['excess_nll_vs_best'],3),round(v['excess_nll_se'],3),round(v['expected_changes'],2))

if __name__=='__main__': main()
