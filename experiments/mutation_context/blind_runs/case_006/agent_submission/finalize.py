from nonparam import *
import pandas as pd
import platform,sys
def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def main():
    c=get10();grid,w,info=fit_np(c[0]);q=predict_np(c[0],grid,w)
    rate=192005/38400000;base=np.array([1-rate,rate/3,rate/3,rate/3]);val=c[2];n=val.sum(axis=1);m=n-val[:,0]
    bins=np.array([0,.002,.004,.006,.01,.02,.05,.1,1.]);ix=np.digitize(1-q[:,0],bins)-1;cal=[]
    for i in range(len(bins)-1):
        take=ix==i;den=n[take].sum()
        cal.append({'lower':float(bins[i]),'upper':float(bins[i+1]),'sites':int(den),'expected_changes':float(n[take]@(1-q[take,0])),'observed_changes':int(m[take].sum())})
    exp=json.loads((ROOT/'exploration.json').read_text());splits={s['name']:s['split'] for s in exp['scaffolds']};parents=fasta(ROOT/'calibration_parent.fasta');children=fasta(ROOT/'calibration_offspring.fasta');diff=[]
    for name,x in parents.items():
        if splits[name]!=2:continue
        code=contexts(x,11)//4;z=(children[name].astype(int)-x)%4;llgain=np.log(q[code,z]/base[z]);diff.append(float(llgain.sum()))
    diff=np.array(diff);validation={'selected_model':'independent_nonparametric_binomial_prior_on_oriented_10mer','window':[-5,4],'assessment_log_likelihood':float(np.sum(val*np.log(q))),'pooled_log_likelihood':float(np.sum(val*np.log(base))),'assessment_log_likelihood_gain':float(diff.sum()),'gain_standard_error_scaffold_clustered':float(diff.std(ddof=1)*np.sqrt(len(diff))),'assessment_scaffolds':len(diff),'calibration_bins':cal,'assessment_expected_changes':float(n@(1-q[:,0])),'assessment_observed_changes':int(m.sum()),'assessment_multiclass_brier':float(np.sum(n*np.sum(q*q,axis=1)-2*np.sum(val*q,axis=1)+n)/n.sum())}
    (ROOT/'selected_validation.json').write_text(json.dumps(validation,indent=2))
    del parents,children,q
    total=c.sum(axis=0);grid,w,info=fit_np(total);pred=predict_np(total,grid,w)
    np.savez_compressed(ROOT/'final_model.npz',relative_probabilities=pred,grid=grid,weights=w)
    (ROOT/'final_model_parameters.json').write_text(json.dumps(info,indent=2))
    new=fasta(ROOT/'new_parent.fasta');queries=pd.read_csv(ROOT/'queries.tsv',sep='\t');assert list(queries.columns)==['contig','position']
    output=np.empty((len(queries),4));parent=np.empty(len(queries),np.uint8)
    for name,group in queries.groupby('contig',sort=False):
        x=new[name];pos=group.position.to_numpy(dtype=int)-1;assert np.all((pos>=0)&(pos<len(x)));code=contexts(x,11)//4;rel=pred[code[pos]];src=x[pos];absolute_prob=np.take_along_axis(rel,(np.arange(4)[None,:]-src.astype(int)[:,None])%4,axis=1);output[group.index]=absolute_prob;parent[group.index]=src
    assert np.isfinite(output).all() and np.all(output>=0) and np.max(np.abs(output.sum(axis=1)-1))<1e-12
    out=queries.copy();out[['p_A','p_C','p_G','p_T']]=output;out.to_csv(ROOT/'predictions.tsv',sep='\t',index=False,float_format='%.15g')
    reread=pd.read_csv(ROOT/'predictions.tsv',sep='\t');assert reread[['contig','position']].equals(queries)
    pp=reread[['p_A','p_C','p_G','p_T']].to_numpy();assert np.isfinite(pp).all() and np.all(pp>=0) and np.max(np.abs(pp.sum(axis=1)-1))<1e-12
    expected=float(np.sum(1-pp[np.arange(len(pp)),parent]));unique=int(queries.drop_duplicates().shape[0])
    inputfiles=['TASK.md','metadata.json','calibration_parent.fasta','calibration_offspring.fasta','new_parent.fasta','queries.tsv']
    summary={'expected_number_of_changed_alleles':expected,'expected_total_changed_alleles':expected,'number_of_queries':len(queries),'number_of_unique_queried_sites':unique,'calibration_sites':64000000,'calibration_changed_alleles':320288,'calibration_pooled_change_rate':.0050045,'model':'Empirical Bayes independent binomial mixture for each alternate allele, conditioned on oriented 10-mer offsets -5 through +4; posterior mean allele probabilities.','seed':0,'threads':4,'input_sha256':{f:sha256(ROOT/f) for f in inputfiles},'output_sha256':sha256(ROOT/'predictions.tsv'),'maximum_row_sum_error_after_serialization':float(np.max(np.abs(pp.sum(axis=1)-1))),'min_probability':float(pp.min()),'max_probability':float(pp.max()),'runtime':{'python':sys.version,'numpy':np.__version__,'platform':platform.platform()},'validation':validation}
    (ROOT/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps({k:v for k,v in summary.items() if k not in ['validation','input_sha256']},indent=2),flush=True)
if __name__=='__main__':main()
