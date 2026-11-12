from rule_models import *

def residual_search(raw,pr):
    c=np.arange(16384);b=c//64%4;n=raw.sum(1);flip=b>=2;can=np.minimum(b,3-b);items=[]
    for pair in combinations([0,1,2,4,5,6],2):
        _,_,eq,at,_,_,_=rule_features(c,pair)
        for center in [0,1]:
            for dest in [j for j in range(4) if j!=center]:
                mask=can==center;target=np.where(flip,3-dest,dest)
                obs=np.bincount(eq[mask],weights=raw[c[mask],target[mask]],minlength=2)
                expected=np.bincount(eq[mask],weights=n[mask]*pr[c[mask],target[mask]],minlength=2)
                stat=np.sum((obs-expected)**2/expected)
                items.append({'pair_offsets':[pair[0]-3,pair[1]-3],'canonical_center':BASES[center],'canonical_target':BASES[dest],'statistic':float(stat),'observed':obs.tolist(),'expected':expected.tolist()})
    return sorted(items,key=lambda x:-x['statistic'])

def main():
    c,y,g=load_data();raw=counts(c,y);folds=np.array(json.loads((ROOT/'validation.json').read_text())['fold_by_contig_index']);scores={};searches=[];oof=np.zeros((len(c),4))
    for fold in range(5):
        test=folds[g]==fold;train=raw-counts(c[test],y[test]);first=fit_rule(train,'pooled_reduced');searches.append(residual_search(train,first)[:3])
        for name in ['pooled_reduced','pooled_three_rules']:
            pc=fit_rule(train,name)[c[test]];b=c[test]//64%4;yt=y[test];gt=g[test];ll=-np.log(pc[np.arange(len(pc)),yt])
            rec=scores.setdefault(name,{'nll':0,'expected_changes':0,'brier_sum':0,'contig_nll':{}})
            rec['nll']+=float(ll.sum());rec['expected_changes']+=float((1-pc[np.arange(len(pc)),b]).sum());rec['brier_sum']+=float((np.sum(pc**2,axis=1)-2*pc[np.arange(len(pc)),yt]+1).sum())
            for gi in np.unique(gt):rec['contig_nll'][int(gi)]=float(ll[gt==gi].sum())
            if name=='pooled_three_rules':oof[test]=pc
    final=fit_rule(raw,'pooled_three_rules');res=residual_search(raw,final)
    b=c//64%4;mutation=y!=b;pchange=1-oof[np.arange(len(c)),b]; spatial=[]
    for gi in range(80):
        mask=g==gi;obs=int(mutation[mask].sum());exp=float(pchange[mask].sum());var=float(np.sum(pchange[mask]*(1-pchange[mask])))
        spatial.append({'contig':gi,'observed':obs,'expected':exp,'variance':var,'z':(obs-exp)/np.sqrt(var)})
    chi=sum((r['observed']-r['expected'])**2/r['variance'] for r in spatial)
    calibration=[];bb,cp,eq,at,hot,_,_=rule_features(c);e2=c//1024%4==c//256%4;ids=np.minimum(bb,3-bb)*4+hot*2+e2
    for group in range(8):
        mask=ids==group;calibration.append({'group':group,'sites':int(mask.sum()),'observed':int(mutation[mask].sum()),'expected':float(pchange[mask].sum())})
    delta=np.array([scores['pooled_reduced']['contig_nll'][i]-scores['pooled_three_rules']['contig_nll'][i] for i in range(80)])
    result={'scores':scores,'secondary_fold_discovery':searches,'residual_checks_top10':res[:10],'spatial_contigs':spatial,'spatial_chi_square':chi,'spatial_approx_p':float(chi2.sf(chi,79)),'group_calibration':calibration,'n':len(c),'observed_changes':int(mutation.sum()),'three_rule_improvement':float(delta.sum()),'three_rule_improvement_se':float(delta.std(ddof=1)*np.sqrt(80))}
    (ROOT/'final_validation.json').write_text(json.dumps(result,indent=2));np.save(ROOT/'final_context_probabilities.npy',final);np.save(ROOT/'final_oof_probabilities.npy',oof)
    print(json.dumps({k:v for k,v in result.items() if k not in ['scores','spatial_contigs','secondary_fold_discovery']},indent=2))
    print('scores', {name:{k:v for k,v in rec.items() if k!='contig_nll'} for name,rec in scores.items()})
    print('secondary discovery',[s[0] for s in searches])

if __name__=='__main__':main()
