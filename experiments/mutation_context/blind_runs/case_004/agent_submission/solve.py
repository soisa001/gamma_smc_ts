from check_structure import *
from scipy.stats import norm
from sklearn.metrics import roc_auc_score
import sys,platform,time,importlib.metadata

def fit_model(xs,ys):
    n,m=get_counts(xs,ys,-4,8);d=dest_counts(xs,ys);rate_fit=logitnormal_fit(n,m);destination_fit=dirichlet_fit(triples(d))
    rates=mix_rates(n,m,rate_fit);tau=destination_fit['tau'];q=(d+tau/3)/(m[:,None]+tau);central=np.arange(len(n))//4**3%4;q[np.arange(len(q)),central]=0
    probs=rates[:,None]*q;probs[np.arange(len(q)),central]=1-rates
    prior_p=np.array(rate_fit['probabilities']);prior_w=np.array(rate_fit['weights']);lp=m[:,None]*np.log(prior_p)+(n-m)[:,None]*np.log1p(-prior_p)+np.log(prior_w);post=np.exp(lp-logsumexp(lp,axis=1)[:,None]);var=np.maximum(post@(prior_p**2)-rates**2,0)
    return probs,rates,var,{'rate_prior':rate_fit,'destination_prior':destination_fit,'n':int(n.sum()),'changes':int(m.sum()),'context_counts_min_max':[int(n.min()),int(n.max())]},n,m
def per_scaffold_metrics(x,y,pr):
    m=x!=y;idx=np.arange(len(x));rate=1-pr[idx,x];return {'n':len(x),'changes':int(m.sum()),'expected_changes':float(rate.sum()),'nll':float(-np.log(pr[idx,y]).sum()),'binary_nll':float(-(m*np.log(rate)+(~m)*np.log1p(-rate)).sum()),'brier_sum':float(np.sum(pr**2)-2*pr[idx,y].sum()+len(x))}
def evaluation(xs,ys,splits):
    train=splits['discovery']+splits['validation'];test=splits['test'];xx=[xs[i]for i in train];yy=[ys[i]for i in train]
    pr,rate,var,fit,n,m=fit_model(xx,yy);p=m.sum()/n.sum();central=np.arange(len(n))//4**3%4
    uniform=np.full((4,4),p/3);np.fill_diagonal(uniform,1-p)
    source=np.zeros((4,4))
    for x,y in zip(xx,yy):source+=np.bincount(x.astype(int)*4+y,minlength=16).reshape(4,4)
    source=(source+.5)/(source.sum(1,keepdims=True)+2)
    no_dest=np.repeat((rate/3)[:,None],4,axis=1);no_dest[np.arange(len(rate)),central]=1-rate
    models={'uniform':[],'source_matrix':[],'context_rate_uniform_destination':[],'selected':[]}
    true=[];rs=[];gain=[]
    for i in test:
        x,y=xs[i],ys[i];z=window(x,-4,8)
        for name,preds in [('uniform',uniform[x]),('source_matrix',source[x]),('context_rate_uniform_destination',no_dest[z]),('selected',pr[z])]:models[name].append(per_scaffold_metrics(x,y,preds))
        true.append(x!=y);rs.append(rate[z]);gain.append(models['uniform'][-1]['nll']-models['selected'][-1]['nll'])
    out={'fit':fit,'train_scaffolds':train,'test_scaffolds':test,'models':{}}
    for name,rows in models.items():
        total={k:sum(r[k]for r in rows)for k in rows[0]};total['mean_nll']=total['nll']/total['n'];total['mean_binary_nll']=total['binary_nll']/total['n'];total['mean_brier']=total['brier_sum']/total['n'];out['models'][name]={'total':total,'per_scaffold':rows}
    rng=np.random.default_rng(0);gain=np.array(gain);bs=gain[rng.integers(0,len(gain),size=(5000,len(gain)))].sum(1);out['selected_gain_nats_vs_uniform']={'total':float(gain.sum()),'scaffold_bootstrap_95_percentile_interval':np.quantile(bs,[.025,.975]).tolist(),'resamples':5000,'seed':0}
    true=np.concatenate(true);rs=np.concatenate(rs);out['mutation_auc']=float(roc_auc_score(true,rs));cuts=np.unique(np.quantile(rs,np.linspace(0,1,11)));groups=np.searchsorted(cuts[1:-1],rs);out['calibration_bins']=[{'n':int((groups==i).sum()),'expected_changes':float(rs[groups==i].sum()),'observed_changes':int(true[groups==i].sum()),'mean_predicted_rate':float(rs[groups==i].mean())}for i in range(len(cuts)-1)]
    (ROOT/'validation.json').write_text(json.dumps(out,indent=2));print('FINAL TEST',json.dumps({k:v for k,v in out.items()if k not in ['fit','models','calibration_bins','train_scaffolds','test_scaffolds']},indent=2),flush=True)
    print('TEST MODELS',json.dumps({k:v['total']for k,v in out['models'].items()},indent=2),flush=True)
    return out
def run():
    started=time.time();names,xs,ys=inputs();splits=json.loads((ROOT/'splits.json').read_text());validation=evaluation(xs,ys,splits)
    probs,rates,var,fit,n,m=fit_model(xs,ys);fit.update({'model':'circular_8mer_logistic_normal_rate_symmetric_dirichlet_destination','offsets':list(range(-4,4)),'seed':0,'threads':4,'memory_limit_gb':26,'calibration_contigs':names});(ROOT/'model.json').write_text(json.dumps(fit,indent=2));np.savez(ROOT/'model_context_table.npz',probabilities=probs,change_rate=rates,posterior_rate_variance=var,opportunities=n,changes=m)
    import pandas as pd
    new=fasta(ROOT/'new_parent.fasta');newx={k:encode(v)for k,v in new.items()};newz={k:window(v,-4,8)for k,v in newx.items()}
    queries=pd.read_csv(ROOT/'queries.tsv',sep='\t',dtype={'contig':str,'position':np.int64});assert list(queries.columns)==['contig','position'];assert not queries.isna().any().any()
    N=len(queries);z=np.empty(N,dtype=np.int64);parent=np.empty(N,dtype=np.uint8)
    for name,idx in queries.groupby('contig',sort=False).groups.items():
        assert name in newx;pos=queries.loc[idx,'position'].to_numpy()-1;assert np.all((pos>=0)&(pos<len(newx[name])));z[idx]=newz[name][pos];parent[idx]=newx[name][pos]
    pred=probs[z];assert pred.shape==(N,4) and np.isfinite(pred).all() and (pred>=0).all();assert np.max(np.abs(pred.sum(1)-1))<1e-12
    expected=float((1-pred[np.arange(N),parent]).sum());counts=np.bincount(z,minlength=len(rates));parameter_var=float(np.sum(counts**2*var));predictive_var=float(np.sum(counts*rates*(1-rates)+counts*(counts-1)*var));duplicate_rows=int(queries.duplicated().sum())
    assert duplicate_rows==0,'Predictive-variance calculation assumes distinct queried sites.'
    for i,b in enumerate(BASES):queries['p_'+b]=pred[:,i]
    tmp=ROOT/'predictions.tsv.tmp';queries.to_csv(tmp,sep='\t',index=False,float_format='%.12g');tmp.replace(ROOT/'predictions.tsv')
    # Re-read final on-disk output to validate rounding and row integrity.
    checked=pd.read_csv(ROOT/'predictions.tsv',sep='\t');orig=pd.read_csv(ROOT/'queries.tsv',sep='\t');assert checked[['contig','position']].equals(orig);pp=checked[['p_'+b for b in BASES]].to_numpy();assert np.isfinite(pp).all() and (pp>=0).all();err=float(np.max(np.abs(pp.sum(1)-1)));assert err<1e-10
    actual_expected=float((1-pp[np.arange(N),parent]).sum());assert abs(actual_expected-expected)<1e-6
    audit=json.loads((ROOT/'input_audit.json').read_text());hashes={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()for name in audit['sha256']};assert hashes==audit['sha256']
    summary={'expected_number_of_changed_alleles':actual_expected,'expected_total_changed_alleles':actual_expected,'number_of_queried_sites':N,'unique_queried_sites':N-duplicate_rows,'mean_predicted_change_probability':actual_expected/N,'model':fit['model'],'context_offsets':fit['offsets'],'calibration_sites':int(n.sum()),'calibration_changes':int(m.sum()),'fitted_logit_rate_mean':fit['rate_prior']['mu'],'fitted_logit_rate_sd':fit['rate_prior']['sigma'],'fitted_destination_dirichlet_total_concentration':fit['destination_prior']['tau'],'conditional_posterior_sd_expected_change_count':parameter_var**.5,'conditional_posterior_predictive_sd_change_count':predictive_var**.5,'uncertainty_caveat':'Conditional on the selected model and fitted hyperparameters; excludes hyperparameter and model uncertainty.','max_probability_row_sum_error':err,'input_hashes_verified':True,'seed':0,'threads':4,'runtime_seconds':time.time()-started}
    (ROOT/'summary.json').write_text(json.dumps(summary,indent=2));provenance={'python':sys.version,'platform':platform.platform(),'packages':{k:importlib.metadata.version(k)for k in ['numpy','scipy','pandas','scikit-learn']},'input_sha256':hashes,'output_sha256':{'predictions.tsv':hashlib.sha256((ROOT/'predictions.tsv').read_bytes()).hexdigest()},'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in ROOT.glob('*.py')},'seed':0,'threads':4};(ROOT/'provenance.json').write_text(json.dumps(provenance,indent=2));print('SUMMARY',json.dumps(summary,indent=2),flush=True)
if __name__=='__main__':run()
