from models import *
from itertools import combinations

def rule_features(c,pair=(0,5)):
    b=c//64%4; cp=((b==1)&(c//16%4==2))|((b==2)&(c//256%4==1))
    left=c//(4**(6-pair[0]))%4;right=c//(4**(6-pair[1]))%4
    eq=left==right;at=(b==0)|(b==3);hot=np.where(at,eq,cp)
    return b,cp,eq,at,hot,left,right

def fit_rule(raw,name,pair=(0,5)):
    c=np.arange(16384);b,cp,eq,at,hot,left,right=rule_features(c,pair)
    sym=name.startswith('pooled');flip=b>=2; canonical=np.minimum(b,3-b) if sym else b
    vals=raw.copy()
    if sym: vals[flip]=vals[flip,::-1]
    if 'identity' in name:
        feature=np.where(at,np.where(eq,left+1,0),cp.astype(int));ids=canonical*5+feature; size=20
    elif 'pair' in name:
        feature=np.where(at,left*4+right,cp.astype(int));ids=canonical*16+feature;size=64
    elif 'all_eq' in name:
        ids=canonical*4+cp.astype(int)*2+eq;size=16
    else: ids=canonical*2+hot;size=8
    tab=np.zeros((size,4))
    for j in range(4): np.add.at(tab[:,j],ids,vals[:,j])
    probs=(tab+.5)/(tab.sum(1)[:,None]+2)
    pred=probs[ids]
    if sym: pred[flip]=pred[flip,::-1]
    if name in ['pooled_reduced','pooled_three_rules']:
        # Share complementary event rates and rates uninfluenced by the trigger.
        # Estimate each substitution from its complete applicable exposure set.
        pred=np.zeros((16384,4)); n=raw.sum(1)
        categories=['at_other','at_transition_cold','at_transition_hot','at_transversion_cold','at_transversion_hot','cg_transversions','cg_transition_cold','cg_transition_hot']
        if name=='pooled_three_rules':categories=['at_other_cold','at_other_hot']+categories[1:]
        for category in categories:
            if category.startswith('at_other'):
                mask=at;targets=3-b
                if name=='pooled_three_rules':mask=mask&((c//1024%4==c//256%4)==category.endswith('_hot'))
            elif category.startswith('at_'):
                mask=at&(eq==category.endswith('_hot'))
                targets=np.where(b==0,2,1) if '_transition_' in category else np.where(b==0,1,2)
            elif category=='cg_transversions':
                mask=~at
                # Two exchangeable alternate alleles per C/G site.
                numer=raw[c[mask],np.where(b[mask]==1,0,3)].sum()+raw[c[mask],3-b[mask]].sum()
                rate=(numer+.5)/(2*n[mask].sum()+1)
                pred[c[mask],np.where(b[mask]==1,0,3)]=rate;pred[c[mask],3-b[mask]]=rate
                continue
            else: mask=(~at)&(cp==category.endswith('_hot'));targets=np.where(b==1,3,0)
            rate=(raw[c[mask],targets[mask]].sum()+.5)/(n[mask].sum()+1)
            pred[c[mask],targets[mask]]=rate
        pred[c,b]=1-pred.sum(1)
    return pred

def discover_pair(raw):
    c=np.arange(16384);n=raw.sum(1);b=c//64%4;obs=n-raw[c,b];pr=fit_tables(raw,False)['cpg'];mu=n*(1-pr[c,b]);out=[]
    for pair in combinations([0,1,2,4,5,6],2):
        _,_,eq,at,_,_,_=rule_features(c,pair);ids=b*2+eq
        oo=np.bincount(ids,weights=obs,minlength=8);ee=np.bincount(ids,weights=mu,minlength=8)
        out.append((float(np.sum((oo-ee)**2/ee)),pair))
    return sorted(out,reverse=True)

def main():
    c,y,g=load_data();allraw=counts(c,y);folds=np.array(json.loads((ROOT/'validation.json').read_text())['fold_by_contig_index'])
    names=['oriented_rule','pooled_rule','oriented_identity','pooled_identity','oriented_pair','pooled_pair','oriented_all_eq','pooled_all_eq','pooled_reduced']
    results={name:{'nll':0,'expected_changes':0,'contig_nll':{}} for name in names};selected=[];oof=np.zeros((len(c),4))
    for fold in range(5):
        test=folds[g]==fold;train=allraw-counts(c[test],y[test]);disc=discover_pair(train);pair=disc[0][1]
        selected.append({'fold':fold,'ranked_equalities':disc}); print('fold',fold,'best equalities',disc[:3],flush=True)
        for name in names:
            tab=fit_rule(train,name,pair);pc=tab[c[test]];b=c[test]//64%4;gt=g[test]
            losses=-np.log(pc[np.arange(len(pc)),y[test]])
            results[name]['nll']+=float(losses.sum());results[name]['expected_changes']+=float((1-pc[np.arange(len(pc)),b]).sum())
            for gi in np.unique(gt):results[name]['contig_nll'][int(gi)]=float(losses[gt==gi].sum())
            if name=='pooled_reduced': oof[test]=pc
    best=min(results,key=lambda k:results[k]['nll'])
    for name,v in results.items():
        delta=np.array([v['contig_nll'][i]-results[best]['contig_nll'][i] for i in range(80)])
        v['excess_nll']=float(delta.sum());v['excess_nll_se']=float(delta.std(ddof=1)*np.sqrt(80))
    result={'best':best,'fold_discovery':selected,'scores':results}
    (ROOT/'rule_validation.json').write_text(json.dumps(result,indent=2));np.save(ROOT/'rule_oof_probabilities.npy',oof)
    for name in sorted(results,key=lambda k:results[k]['nll']):print(name,{k:v for k,v in results[name].items() if k!='contig_nll'})

if __name__=='__main__':main()
