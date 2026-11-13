from analysis import *
from model import context_codes
import pandas as pd, hashlib
p,o,n=load();sp=json.loads((BASE/'split.json').read_text());rng=np.random.default_rng(0)
assert sorted(np.bincount(list(sp.values())).tolist())==[16,16,48]
checks=0
for k,v in n.items():
 c=context_codes(v)
 for pos in np.unique(np.r_[0,1,2,len(v)-2,len(v)-1,rng.integers(0,len(v),size=1000)]):
  manual=0
  for offset in [-3,-2,-1,0,1,2]:manual=manual*4+int(v[(int(pos)+offset)%len(v)])
  assert c[pos]==manual;checks+=1
q=pd.read_csv(BASE/'queries.tsv',sep='\t');out=pd.read_csv(BASE/'predictions.tsv',sep='\t');assert out.iloc[:,:2].equals(q)
pr=out.iloc[:,2:].to_numpy();assert pr.shape==(1000000,4);assert np.isfinite(pr).all();assert ((pr>=0)&(pr<=1)).all();assert np.max(abs(pr.sum(1)-1))<1e-11
summary=json.loads((BASE/'summary.json').read_text());assert hashlib.sha256((BASE/'predictions.tsv').read_bytes()).hexdigest()==summary['output_sha256']
result={'scalar_circular_context_checks':checks,'split_counts':[list(sp.values()).count(i) for i in range(3)],'prediction_rows':len(out),'unique_query_sites':len(q.drop_duplicates()),'all_probabilities_valid':True,'query_order_and_identity_verified':True,'sha256_verified':True}
(BASE/'verification.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
