import numpy as np
import pandas as pd
import tskit
from io import StringIO
from gamma_smc_aou.stopped_campaign_report import evaluate_available, truth_at_positions, examples
from gamma_smc_aou.origin_onset_analysis import upper_tail_boundary, rank_p


def test_af_critical_boundary_matches_empirical_calls_with_ties():
    for n in [18, 19, 97, 98, 1000]:
        null = np.round(np.linspace(0, 1, n), 1)
        boundary = upper_tail_boundary(null)
        scores = np.unique(np.concatenate([null, np.nextafter(null, np.inf), [0, 1]]))
        np.testing.assert_array_equal(scores > boundary, rank_p(null, scores) <= .05)

def test_actual_null_denominators_and_missing_pair_no_call():
    rows=[]
    for family,n in [('I50',97),('I10',98)]:
        for i in range(n):
            rows.append(dict(task_id=f'{family}/null/{i}',family=family,role='null',s=0,
                source='truth',method='alt_alt',cutoff_years=10000,score=0.5))
        for i,score in enumerate([1.0,0.5,np.nan]):
            rows.append(dict(task_id=f'{family}/selected/{i}',family=family,role='selected',s=.01,
                source='truth',method='alt_alt',cutoff_years=10000,score=score))
    pred,metrics=evaluate_available(pd.DataFrame(rows),{'alpha':.05,'null_replicates':1000,'target_replicates':100})
    for family,n in [('I50',97),('I10',98)]:
        g=pred[pred.family==family]
        np.testing.assert_allclose(g.p,[1/(n+1),1,1])
        row=metrics[metrics.family==family].iloc[0]
        assert row.n==3 and row.calls==1 and row.null_n==n
        assert row.rate==1/3 and row.coverage==2/3 and row.evaluable_rate==1/2
        assert row.ci_low < row.rate < row.ci_high

def test_truth_strict_threshold_and_unavailable_class():
    tables=tskit.TableCollection(10)
    individual=tables.individuals.add_row()
    tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE,time=0,individual=individual)
    tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE,time=0,individual=individual)
    ancestor=tables.nodes.add_row(time=400)
    for node in [0,1]: tables.edges.add_row(0,10,ancestor,node)
    tables.sort()
    scores=truth_at_positions(tables.tree_sequence(),np.array([[0,1]]),np.array([1]),[5],
        {'tmrca_cutoffs_years':[10000,20000],'generation_time_years':25})
    assert scores[scores.method=='all_pairs'].score.tolist()==[0,1]
    assert scores[scores.method=='alt_alt'].score.isna().all()
    assert (scores[scores.method=='alt_alt'].n_pairs==0).all()

def test_representative_is_nearest_median_with_task_id_tie_break():
    data=pd.DataFrame([dict(task_id=f't{i}',family='I10',role='selected',s=.01,sample_af=af)
        for i,af in enumerate([.1,.3,.5,.9])])
    assert examples(data).task_id.tolist()==['t1']

def test_inventory_roundtrip_preserves_calibration_null_role():
    frame=pd.DataFrame([dict(task_id='I50/null/rep0000',family='I50',role='null',s=0,sample_af=.1)])
    restored=pd.read_csv(StringIO(frame.to_csv(index=False)),keep_default_na=False,float_precision='round_trip')
    assert restored.role.tolist()==['null']
    assert examples(restored).role.tolist()==['null']
