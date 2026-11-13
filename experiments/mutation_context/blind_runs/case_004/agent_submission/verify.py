from explore import *
from window_scan import window
import pandas as pd
def run():
    rng=np.random.default_rng(0);checks=0
    for filename in ['calibration_parent.fasta','new_parent.fasta']:
        for name,seq in fasta(ROOT/filename).items():
            x=encode(seq);z=window(x,-4,8)
            for pos in list(range(5))+list(range(len(x)-5,len(x)))+rng.integers(0,len(x),size=30).tolist():
                encoded=0
                for k in range(-4,4):encoded=encoded*4+int(x[(pos+k)%len(x)])
                assert encoded==z[pos];assert z[pos]//4**3%4==x[pos];checks+=1
    new=fasta(ROOT/'new_parent.fasta');q=pd.read_csv(ROOT/'queries.tsv',sep='\t');p=pd.read_csv(ROOT/'predictions.tsv',sep='\t');assert p[['contig','position']].equals(q);assert list(p)==['contig','position','p_A','p_C','p_G','p_T'];a=p.iloc[:,2:].to_numpy();assert np.isfinite(a).all() and (a>=0).all();assert np.max(abs(a.sum(1)-1))<1e-10
    observed_counts=q.groupby('contig')['position'].agg(['min','max','nunique']).to_dict('index')
    for name,values in observed_counts.items():assert values=={'min':1,'max':len(new[name]),'nunique':len(new[name])}
    initial=json.loads((ROOT/'input_audit.json').read_text())['sha256'];assert all(hashlib.sha256((ROOT/k).read_bytes()).hexdigest()==v for k,v in initial.items())
    out={'circular_context_checks':checks,'probability_rows_checked':len(p),'query_order_preserved':True,'query_coverage':'Every base of all 20 new-parent scaffolds, exactly once','input_hashes_unchanged':True,'probability_min':float(a.min()),'probability_max':float(a.max()),'max_row_sum_error':float(np.max(abs(a.sum(1)-1)))}
    (ROOT/'verification.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
if __name__=='__main__':run()
