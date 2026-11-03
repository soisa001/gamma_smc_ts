"""Recent-onset timing and interval scoring contracts."""
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from report_recent_selection_pilot import evaluate,select_prefix


def test_recent_truth_stops_at_ten_kya_and_counts_mistiming():
    positions=np.arange(100)*10000;calls=np.zeros((100,49),bool)
    calls[50,9]=True  # [10,11) only touches the onset boundary.
    r=evaluate(calls,positions)
    assert r['spatial_hit']==1 and r['joint_hit']==0 and r['center_mistimed']==1
    calls[50,:9]=True
    r=evaluate(calls,positions)
    assert r['joint_hit']==1 and r['overlap_kyr']==9
    assert r['temporal_f1']==18/19


def test_recent_joint_hit_requires_same_time_and_place():
    positions=np.arange(100)*10000;calls=np.zeros((100,49),bool)
    calls[0,2]=True;calls[50,30]=True
    r=evaluate(calls,positions)
    assert r['any_call']==1 and r['spatial_hit']==1 and r['joint_hit']==0
    assert r['event_fp']==2 and r['event_fn']==1


def test_prefix_does_not_skip_incomplete_proposals_or_pick_faster_acceptances():
    tasks=[dict(unit=str(i),replicate=i) for i in range(5)]
    complete={'0':dict(accepted=False),'1':dict(accepted=True),'3':dict(accepted=True)}
    assert select_prefix(tasks,complete,2)==['0','1']
    complete['2']=dict(accepted=True)
    assert select_prefix(tasks,complete,2)==['0','1','2']


def test_selection_is_active_for_exactly_four_hundred_reproductions():
    from run_recent_selection_pilot import timing
    age,first=timing('de_novo',10000)
    assert (age,first)==(400,2) and age+1-first+1==400
    age,first=timing('introgression',10000)
    assert (age,first)==(2000,1602) and age+1-first+1==400
    assert (age+1-(first-1))*25==10000


def test_reuse_accepts_count_optimization_but_rejects_changed_scientific_code(tmp_path):
    import json
    import pytest
    import run_recent_selection_pilot as runner
    snapshot=tmp_path/'executed_sources/run_recent_selection_pilot.py'
    snapshot.parent.mkdir()
    source=Path(runner.__file__).read_text()
    snapshot.write_text(source.replace('copies=sim.mutationCounts(p0,muts);','copies=sum(p0.haplosomes.containsMutations(muts));'))
    old=dict(source_sha256=runner.base.digest(snapshot),ancestry_source_sha256='ancestry',score_source_sha256='score')
    manifest=tmp_path/'manifest.json';manifest.write_text(json.dumps(dict(config=old)))
    proof=tmp_path/'count_method_validation/results.json';proof.parent.mkdir()
    proof.write_text(json.dumps([dict(origin=o,trajectory_exact_match=True,genealogical_tables_exact_match=True)
                                for o in ('de_novo','introgression')]))
    cfg=dict(old,source_sha256=runner.base.digest(runner.__file__))
    runner.validate_reuse(tmp_path,cfg)
    with pytest.raises(AssertionError,match='ancestry/scoring'):
        runner.validate_reuse(tmp_path,dict(cfg,score_source_sha256='changed'))
    snapshot.write_text(snapshot.read_text().replace('initializeRecombinationRate(1e-8)','initializeRecombinationRate(2e-8)'))
    old['source_sha256']=runner.base.digest(snapshot);manifest.write_text(json.dumps(dict(config=old)))
    with pytest.raises(AssertionError,match='implementation differs'):
        runner.validate_reuse(tmp_path,cfg)


def test_all_pair_counts_match_explicit_pairs_including_polytomies_and_cutoff_ties():
    import itertools
    import tskit
    from audit_recent_all_pairs_truth import all_pairs_cdf
    tables=tskit.TableCollection(100)
    for _ in range(5):tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE,time=0)
    a=tables.nodes.add_row(time=40);b=tables.nodes.add_row(time=80)
    for child in (0,1,2):tables.edges.add_row(0,100,a,child)
    for child in (a,3,4):tables.edges.add_row(0,100,b,child)
    tables.sort();ts=tables.tree_sequence();cutoffs=np.array([1000,1001,2000,2001])
    ages=np.array([ts.first().tmrca(x,y)*25 for x,y in itertools.combinations(ts.samples(),2)])
    expected=np.mean(ages[:,None]<cutoffs,axis=0)
    np.testing.assert_array_equal(all_pairs_cdf(ts,[0,50],cutoffs),np.tile(expected,(2,1)))
