from model import *
from scipy.stats import chi2
import hashlib, platform, importlib.metadata
INPUTS=['TASK.md','metadata.json','calibration_parent.fasta','calibration_offspring.fasta','new_parent.fasta','queries.tsv']
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def fulltable(rate,dist):
 t=np.zeros((4096,4));t[np.arange(4096),PARENTS]=1-rate;t[np.arange(4096)[:,None],OTHER]=rate[:,None]*dist;return t

def main():
 manifest={f:{'sha256':digest(BASE/f),'bytes':(BASE/f).stat().st_size} for f in INPUTS}
 (BASE/'input_manifest.json').write_text(json.dumps(manifest,indent=2))
 p,o,n=load();sp=json.loads((BASE/'split.json').read_text());fitkeys=[k for k in p if sp[k]!=2];auditkeys=[k for k in p if sp[k]==2]
 ct,alt=sufficient(p,o,fitkeys);ac,aa=sufficient(p,o,auditkeys)
 r,rpar,grid=fit_grid(ct,alt.sum(1));d,dpar=fit_alt(alt)
 simple=(alt.sum(1)+300*alt.sum()/ct.sum())/(ct+300);simpled=(unpack(alt)+1)/(alt.sum(1)[:,None]+3)
 baseline_rate=alt.sum()/ct.sum(); br=np.full(4096,baseline_rate);bd=np.full((4096,3),1/3)
 tables={'selected':fulltable(r,d),'simple_sixmer':fulltable(simple,simpled),'constant':fulltable(br,bd)}
 byscaffold=[];quantiles=np.quantile(r,np.linspace(0,1,11));bins=np.searchsorted(quantiles[1:-1],r,side='right');calibration=[]
 for b in range(10):
  mask=bins==b;calibration.append(dict(bin=b,sites=int(ac[mask].sum()),predicted_changes=float((ac[mask]*r[mask]).sum()),observed_changes=int(aa[mask].sum())))
 for k in auditkeys:
  c=context_codes(p[k]);y=o[k];m=y!=p[k];row={'contig':k,'changed':int(m.sum()),'predicted_changes':float(r[c].sum())}
  for name,tab in tables.items():
   pr=tab[c];row[name+'_logloss']=float(-np.log(pr[np.arange(len(y)),y]).mean());row[name+'_brier']=float(np.mean(np.sum(pr**2,axis=1)-2*pr[np.arange(len(y)),y]+1))
  byscaffold.append(row)
 arr=np.array([a['constant_logloss']-a['selected_logloss'] for a in byscaffold]);rng=np.random.default_rng(0);boot=arr[rng.integers(0,len(arr),size=(20000,len(arr)))].mean(1)
 audit={'fit_scaffolds':len(fitkeys),'audit_scaffolds':len(auditkeys),'audit_sites':int(ac.sum()),'observed_changes':int(aa.sum()),'predicted_changes':float(ac@r),'rate_parameters':rpar,'destination_parameters':dpar,'calibration_deciles':calibration,'scaffolds':byscaffold,'baseline_minus_selected_logloss':float(arr.mean()),'scaffold_bootstrap_95_interval':np.quantile(boot,[.025,.975]).tolist()}
 for name in tables:
  audit[name+'_logloss']=float(np.mean([a[name+'_logloss'] for a in byscaffold]));audit[name+'_brier']=float(np.mean([a[name+'_brier'] for a in byscaffold]))
 (BASE/'audit_validation.json').write_text(json.dumps(audit,indent=2))
 print('AUDIT',json.dumps({k:v for k,v in audit.items() if k not in ['scaffolds','calibration_deciles']},indent=2),flush=True)
 # Fit frozen model on all calibration data.
 ct,alt=sufficient(p,o,list(p));r,rpar,(z,w)=fit_grid(ct,alt.sum(1));d,dpar=fit_alt(alt);table=fulltable(r,d)
 assert np.all(np.isfinite(table)) and np.all(table>0) and np.max(abs(table.sum(1)-1))<1e-12
 np.savez(BASE/'fitted_model.npz',context_counts=ct,destination_counts=alt,probabilities=table,logit_grid=z,prior_weights=w)
 import pandas as pd
 q=pd.read_csv(BASE/'queries.tsv',sep='\t',dtype={'contig':str,'position':np.int64});assert list(q)==['contig','position']
 codes={k:context_codes(v) for k,v in n.items()}; probs=np.zeros((len(q),4));parents=np.zeros(len(q),dtype=np.uint8)
 for k,ix in q.groupby('contig',sort=False).groups.items():
  assert k in n;ix=np.array(ix);pos=q.loc[ix,'position'].to_numpy()-1;assert np.all((pos>=0)&(pos<len(n[k])))
  probs[ix]=table[codes[k][pos]];parents[ix]=n[k][pos]
 exp=float(np.sum(1-probs[np.arange(len(q)),parents]))
 result=q.copy()
 for j,a in enumerate('ACGT'):result['p_'+a]=probs[:,j]
 result.to_csv(BASE/'predictions.tsv',sep='\t',index=False,float_format='%.12g')
 check=pd.read_csv(BASE/'predictions.tsv',sep='\t');vv=check[['p_A','p_C','p_G','p_T']].to_numpy()
 assert check[['contig','position']].equals(q);assert np.isfinite(vv).all() and (vv>=0).all();err=float(abs(vv.sum(1)-1).max());assert err<2e-11
 rounded_exp=float(np.sum(1-vv[np.arange(len(q)),parents]));assert abs(rounded_exp-exp)<1e-6
 for f in INPUTS:assert digest(BASE/f)==manifest[f]['sha256']
 query_unique=int(len(q.drop_duplicates()));newhash={hashlib.sha256(v.tobytes()).hexdigest() for v in n.values()};trainhash={hashlib.sha256(v.tobytes()).hexdigest() for v in p.values()}
 summary={'expected_changed_alleles':exp,'expected_total_number_of_changed_alleles':exp,'query_rows':len(q),'unique_queried_sites':query_unique,'mean_change_probability':exp/len(q),'calibration_sites':int(ct.sum()),'calibration_changes':int(alt.sum()),'model':'circular sixmer offsets -3..+2; empirical-Bayes logit-normal rate and symmetric Dirichlet destination','rate_parameters':rpar,'destination_parameters':dpar,'seed':0,'threads':4,'audit_logloss':audit['selected_logloss'],'audit_baseline_logloss':audit['constant_logloss'],'output_sha256':digest(BASE/'predictions.tsv'),'maximum_probability_sum_error_after_serialization':err,'expected_changes_from_serialized_file':rounded_exp,'query_order_verified':True,'input_hashes_unchanged':True,'identical_train_target_sequences':len(newhash&trainhash),'python':platform.python_version(),'packages':{a:importlib.metadata.version(a) for a in ['numpy','scipy','pandas','scikit-learn']}}
 (BASE/'summary.json').write_text(json.dumps(summary,indent=2));print('SUMMARY',json.dumps(summary,indent=2),flush=True)
if __name__=='__main__':main()
