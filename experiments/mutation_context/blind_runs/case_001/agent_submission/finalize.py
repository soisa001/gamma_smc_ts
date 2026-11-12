import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='4'
import json,hashlib,sys,platform
from concurrent.futures import ThreadPoolExecutor
import numpy as np,pandas as pd,scipy,sklearn
from scipy.stats import chi2
from explore import ROOT,fasta,encode,contexts
from model import load,counts
from advanced import motif_features,tied_fit

ctx,y,c=load();g,eq,f,rev=motif_features(ctx);center=(ctx//64)%4;out=y^center
ct=counts(g,out,4);p=tied_fit(ct,'tied5')
queries=pd.read_csv(ROOT/'queries.tsv',sep='\t',dtype={'contig':str,'position':np.int64})
assert list(queries.columns)==['contig','position']
new=fasta('new_parent.fasta');ans=np.full((len(queries),4),np.nan);qg=np.full(len(queries),-1,int);qref=np.full(len(queries),-1,int)
for name,rows in queries.groupby('contig',sort=False):
    assert name in new
    a=encode(new[name]);pos=rows.position.to_numpy()-1
    assert np.all((pos>=3)&(pos<len(a)-3))
    z=contexts(a,3)[pos-3];groups,_,_,_=motif_features(z)
    ref=a[pos]
    prob=np.take_along_axis(p[groups],np.arange(4)[None,:]^ref[:,None],axis=1)
    ans[rows.index]=prob;qg[rows.index]=groups;qref[rows.index]=ref
assert np.all(np.isfinite(ans)) and np.all(ans>=0) and np.allclose(ans.sum(1),1,atol=1e-14)
assert np.all(qg>=0)
change=1-ans[np.arange(len(ans)),qref];expected=float(change.sum());exposure=np.bincount(qg,minlength=4)
output=queries.copy()
for j,base in enumerate('ACGT'):output['p_'+base]=ans[:,j]
tmp=ROOT/'predictions.tsv.tmp';output.to_csv(tmp,sep='\t',index=False,float_format='%.12g');tmp.replace(ROOT/'predictions.tsv')
# Independent round-trip verification against original query order and reference bases.
read=pd.read_csv(ROOT/'predictions.tsv',sep='\t');assert list(read.columns)==['contig','position','p_A','p_C','p_G','p_T']
assert read[['contig','position']].equals(queries)
arr=read.iloc[:,2:].to_numpy();assert np.all(np.isfinite(arr)) and np.all(arr>=0) and np.max(np.abs(arr.sum(1)-1))<1e-10
expected_roundtrip=float(np.sum(1-arr[np.arange(len(arr)),qref]));assert abs(expected_roundtrip-expected)<1e-5

# Seeded contig bootstrap estimates calibration-parameter uncertainty, conditional on model selection.
bycontig=np.array([counts(g[c==j],out[c==j],4) for j in range(80)])
rng=np.random.default_rng(0);sample=rng.integers(0,80,size=(500,80))
def boot(ix):return float(np.sum(exposure*(1-tied_fit(bycontig[ix].sum(0),'tied5')[:,0])))
with ThreadPoolExecutor(max_workers=4) as pool: boots=np.array(list(pool.map(boot,sample)))
np.savetxt(ROOT/'bootstrap_expected_changes.tsv',boots,header='expected_changes',comments='')
residual={}
for name,feat in [('context3',f['context3']),('context5',f['context5']),('context7',f['context7']),('contig',c),('position_decile',np.tile(np.minimum(np.arange(49994)//5000,9),80))]:
    observed=counts(feat,out,int(feat.max()+1));expected_table=np.zeros_like(observed)
    for k in range(4):
        expected_table[:,k]=np.bincount(feat,weights=p[g,k],minlength=len(observed))
    valid=expected_table.sum(1)>0;stat=float(np.sum((observed[valid]-expected_table[valid])**2/expected_table[valid]));df=int(valid.sum()*3-5)
    residual[name]={'pearson_stat':stat,'df':df,'dispersion':stat/df,'p_approx':float(chi2.sf(stat,df))}
prov=json.loads((ROOT/'provenance.json').read_text())
for name,digest in prov['inputs'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest
summary={'expected_total_changed_alleles':expected,'expected_changed_alleles':expected,'query_count':len(queries),'expected_unchanged_alleles':len(queries)-expected,'model':'five-parameter reverse-complement symmetric CpG and +/-2 equality motif model','model_selection':'parsimonious model within one paired-contig standard error of best cross-validated log loss','calibration_sites_used':len(y),'calibration_changes_used':int(np.sum(out!=0)),'seed':0,'threads':4,'bootstrap_replicates':500,'parameter_uncertainty_bootstrap_se':float(boots.std(ddof=1)),'parameter_uncertainty_bootstrap_95_percent_interval':np.quantile(boots,[.025,.975]).tolist(),'conditional_independent_site_count_sd':float(np.sqrt(np.sum(change*(1-change)))),'query_group_exposures':exposure.tolist(),'xor_outcome_probabilities':p.tolist(),'max_roundtrip_row_sum_error':float(np.max(np.abs(arr.sum(1)-1))),'expected_changes_from_serialized_predictions':expected_roundtrip,'predictions_sha256':hashlib.sha256((ROOT/'predictions.tsv').read_bytes()).hexdigest(),'residual_diagnostics':residual,'software':{'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'scikit_learn':sklearn.__version__,'platform':platform.platform()}}
(ROOT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
