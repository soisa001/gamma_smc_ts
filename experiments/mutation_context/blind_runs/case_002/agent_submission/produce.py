from rule_models import *
import csv, itertools, platform
from datetime import datetime, timezone

def parameter_summary(raw, query_counts):
    c=np.arange(16384);b,cp,eq,at,_,_,_=rule_features(c);n=raw.sum(1);e2=c//1024%4==c//256%4
    rows=[]
    for label in ['at_other_cold','at_other_hot','at_transition_cold','at_transition_hot','at_transversion_cold','at_transversion_hot','cg_transversions','cg_transition_cold','cg_transition_hot']:
        multiplicity=1
        if label.startswith('at_other'):
            mask=at&(e2==label.endswith('_hot'));target=3-b
        elif label.startswith('at_'):
            mask=at&(eq==label.endswith('_hot'));target=np.where(b==0,2,1) if '_transition_' in label else np.where(b==0,1,2)
        elif label=='cg_transversions':
            mask=~at;multiplicity=2;target=np.where(b==1,0,3)
        else:mask=(~at)&(cp==label.endswith('_hot'));target=np.where(b==1,3,0)
        events=int(raw[c[mask],target[mask]].sum())
        if multiplicity==2:events+=int(raw[c[mask],3-b[mask]].sum())
        exposure=int(multiplicity*n[mask].sum());query_exposure=int(multiplicity*query_counts[mask].sum())
        rate=(events+.5)/(exposure+1)
        variance=(events+.5)*(exposure-events+.5)/((exposure+1)**2*(exposure+2))
        rows.append({'parameter':label,'events':events,'calibration_exposure':exposure,'query_exposure':query_exposure,'rate':rate,'approximate_rate_sd':float(np.sqrt(variance)),'expected_query_changes':query_exposure*rate,'expected_count_parameter_variance_contribution':query_exposure**2*variance})
    return rows

def main():
    diagnostics=json.loads((ROOT/'input_diagnostics.json').read_text())
    for name,h in diagnostics['input_sha256'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h,name
    raw=np.load(ROOT/'all_context_outcomes.npy');tab=fit_rule(raw,'pooled_three_rules');new=fasta('new_parent.fasta')
    codes={k:context(v,3) for k,v in new.items()};nq=0;expect=0.;process_var=0.;qcounts=np.zeros(16384,dtype=int)
    minimum=1.;maximum=0.
    with (ROOT/'queries.tsv').open(newline='') as f,(ROOT/'predictions.tsv.tmp').open('w',newline='') as out:
        reader=csv.DictReader(f,delimiter='\t');out.write('contig\tposition\tp_A\tp_C\tp_G\tp_T\n')
        for row in reader:
            k=row['contig'];pos=int(row['position']);assert 4<=pos<=len(new[k])-3
            code=int(codes[k][pos-4]);p=tab[code];parent=int(new[k][pos-1]);change=1-p[parent]
            out.write(k+'\t'+str(pos)+'\t'+'\t'.join(format(float(v),'.15g') for v in p)+'\n')
            nq+=1;expect+=change;process_var+=change*(1-change);qcounts[code]+=1;minimum=min(minimum,float(p.min()));maximum=max(maximum,float(p.max()))
    (ROOT/'predictions.tsv.tmp').replace(ROOT/'predictions.tsv')
    reread_count=0;reread_expected=0.;max_error=0.
    with (ROOT/'queries.tsv').open(newline='') as qf,(ROOT/'predictions.tsv').open(newline='') as pf:
        qreader=csv.DictReader(qf,delimiter='\t');preader=csv.DictReader(pf,delimiter='\t')
        assert preader.fieldnames==['contig','position','p_A','p_C','p_G','p_T']
        for q,p in itertools.zip_longest(qreader,preader):
            assert q is not None and p is not None
            assert q['contig']==p['contig'] and q['position']==p['position']
            vals=np.array([float(p['p_'+base]) for base in BASES]);assert np.all(np.isfinite(vals)) and np.all(vals>=0) and np.all(vals<=1)
            error=abs(vals.sum()-1);assert error<1e-12;max_error=max(max_error,float(error))
            reread_count+=1;reread_expected+=1-vals[int(new[q['contig']][int(q['position'])-1])]
    assert reread_count==nq==999880 and abs(reread_expected-expect)<1e-6
    parameters=parameter_summary(raw,qcounts);parameter_var=sum(p['expected_count_parameter_variance_contribution'] for p in parameters)
    assert abs(sum(p['expected_query_changes'] for p in parameters)-expect)<1e-5
    summary={'expected_number_of_changed_alleles':expect,'number_of_queried_sites':nq,'expected_changed_fraction':expect/nq,'model':'nine-rate pooled substitution model with CpG, oriented (-3,+2) equality, and oriented (-2,-1) equality','calibration_sites_used':int(raw.sum()),'calibration_changes':int(raw.sum()-raw[np.arange(16384),np.arange(16384)//64%4].sum()),'seed':0,'threads':4,'memory_limit_gib':26,'parameters':parameters,'uncertainty':{'approximate_parameter_sd_of_expected_count':float(np.sqrt(parameter_var)),'conditional_independent_site_sd_of_realized_count':float(np.sqrt(process_var)),'approximate_combined_predictive_sd':float(np.sqrt(parameter_var+process_var)),'note':'Approximate: ignores cross-parameter covariance, adaptive model-selection uncertainty, and any residual dependence across sites.'},'verification':{'rows_revalidated':reread_count,'max_absolute_row_sum_error':max_error,'minimum_probability':minimum,'maximum_probability':maximum,'reread_expected_changes':reread_expected,'six_input_hashes_unchanged':True,'predictions_sha256':hashlib.sha256((ROOT/'predictions.tsv').read_bytes()).hexdigest()},'provenance':{'created_utc':datetime.now(timezone.utc).isoformat(),'task_directory':str(ROOT),'python':platform.python_version(),'numpy':np.__version__,'input_sha256':diagnostics['input_sha256']}}
    (ROOT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ['parameters','provenance']},indent=2));print('parameters',json.dumps(parameters,indent=2))

if __name__=='__main__':main()
