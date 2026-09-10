"""Scientific contracts for ancestry labeling and regional calibration."""
import numpy as np
import tskit
from gamma_smc_aou.joint_scan_labels import ancestry_mask
from gamma_smc_aou.joint_scan_profiles import nearest_markers,count_classes
from gamma_smc_aou.joint_scan_evaluation import components,extra_components,binned_region_score,rank_p,make_folds
from gamma_smc_aou.joint_truth import pair_ages


def test_migration_intervals_and_haplotype_classes():
    tables=tskit.TableCollection(100)
    for _ in range(4):tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE,time=0)
    ancestor=tables.nodes.add_row(time=10)
    root=tables.nodes.add_row(time=20)
    for u in (0,2):tables.edges.add_row(0,100,ancestor,u)
    for u in (ancestor,1,3):tables.edges.add_row(0,100,root,u)
    tables.sort();tree=tables.tree_sequence().first()
    mapping=np.array([0,1,2,3,-1,-1])
    migrations=np.array([[10,30,ancestor]],dtype=float)
    assert ancestry_mask(tree,migrations,mapping,10).tolist()==[True,False,True,False]
    assert not ancestry_mask(tree,migrations,mapping,30).any()
    # ALT haplotypes 0 and 2 belong to different heterozygous diploid individuals.
    pairs=np.array([[0,2],[1,3],[0,1],[2,3]])
    classes=np.array([1,0,1,0])[pairs].sum(axis=1)
    recent=np.array([[1,1],[0,1],[0,0],[0,0]],dtype=bool)
    counts=count_classes(recent,classes)
    assert counts[:,0].tolist()==[0,0,1,1]
    assert np.array_equal(counts[:3].sum(axis=0),counts[3])


def test_joint_fraction_fixation_and_missing_pairs():
    # RR=2, AR=3, AA=5, all=10. Recent fractions .5,0,1,.6.
    n=np.array([[2,3,5,10],[0,0,10,10],[10,0,0,10]])
    counts=np.array([[[[1],[0],[5],[6]],[[0],[0],[8],[8]],[[3],[0],[0],[3]]]])
    all_pairs,joint=components(counts,n)
    np.testing.assert_allclose(joint[0,:,0],[.25,.8,0])
    np.testing.assert_allclose(all_pairs[0,:,0],[.6,.8,.3])
    mass,excess=extra_components(counts,n)
    np.testing.assert_allclose(mass[0,:,0],[.5,.8,0])
    np.testing.assert_allclose(excess[0,:,0],[.4,.8,0])


def test_spatial_runs_require_distinct_adjacent_physical_bins():
    assert binned_region_score([.8,.9],[1,2],5,10,2)==0
    assert binned_region_score([.8,.9],[1,11],5,10,2)==.8
    assert binned_region_score([.8,.9],[1,21],5,10,2)==0
    assert nearest_markers(np.array([0,10,20]),np.array([5,16]),5).tolist()==[0,0,1]
    assert nearest_markers(np.array([0]),np.array([],dtype=int),5).tolist()==[-1]


def test_conservative_ranks_ties_and_whole_region_folds():
    np.testing.assert_allclose(rank_p([0,1,1,2],[0,1,3]),[1,.8,.2])
    # Zero evidence never acquires significance from empty pair classes.
    assert rank_p(np.zeros(400),[0])[0]==1
    items=[dict(onset_years=o,replicate=r) for o,n in [(0,1000),(10000,100),(50000,100)] for r in range(n)]
    fold=make_folds(items,20380101)
    assert np.array_equal(fold,make_folds(items,20380101))
    onset=np.array([i['onset_years'] for i in items])
    for f in range(5):
        assert [np.sum((fold==f)&(onset==o)) for o in (0,10000,50000)]==[200,20,20]
        fit=(onset==0)&np.isin(fold,[(f+1)%5,(f+2)%5])
        calibration=(onset==0)&np.isin(fold,[(f+3)%5,(f+4)%5])
        assert not np.any(fit&calibration) and not np.any((fit|calibration)&(fold==f))


def test_native_pair_ages_against_tskit_every_marginal_tree():
    import msprime
    ts=msprime.sim_ancestry(samples=10,population_size=100,sequence_length=1000,
                           recombination_rate=.0001,random_seed=20380101)
    pairs=np.array([(a,b) for a in ts.samples() for b in ts.samples() if a<b],dtype=np.int32)
    for tree in ts.trees():
        expected=np.array([tree.tmrca(int(a),int(b)) for a,b in pairs])
        np.testing.assert_array_equal(pair_ages(tree,pairs,ts.nodes_time),expected)


def test_end_to_end_region_evaluation_and_figures(tmp_path,monkeypatch):
    """A known separable signal validates fold aggregation and persisted outputs."""
    import json
    import pandas as pd
    import gamma_smc_aou.joint_scan_evaluation as evaluation
    items=[dict(key=f'{o}/{r}',onset_years=o,replicate=r) for o,n in [(0,1000),(10000,100),(50000,100)] for r in range(n)]
    config=dict(seed_base=20380101,tmrca_cutoffs_years=[5000,10000,20000,30000,40000,50000],scored_length_bp=30,stride_bp=10)
    (tmp_path/'profile_manifest.json').write_text(json.dumps(dict(items=items,config=config)))
    selected_methods=[m for m in evaluation.methods(config['tmrca_cutoffs_years']) if m['name'] in [evaluation.PRIMARY,'all_r1_T50000','af_r1','mass_g50_r1_T50000','excess_g50_r1_T50000']]
    monkeypatch.setattr(evaluation,'methods',lambda cutoffs:selected_methods)
    def fake_load(payload):
        item,_=payload
        recent=[1,0,5,6] if item['onset_years'] else [1,0,1,2]
        counts=np.broadcast_to(np.array(recent)[None,None,:,None],(2,3,4,6)).copy()
        feature=dict(grid_counts=counts,site_counts=counts,grid_n=np.tile([2,3,5,10],(3,1)),site_n=np.tile([2,3,5,10],(3,1)),
                     grid=np.array([0,10,20]),sites=np.array([0,10,20]),site_af=np.array([.5,.5,.5]),grid_markers=np.array([0,1,2]))
        receipt=dict(markers=3,valid_grid_positions=3,verification={})
        return feature,receipt,'synthetic'
    monkeypatch.setattr(evaluation,'load_one',fake_load)
    evaluation.evaluate(tmp_path)
    result=pd.read_csv(tmp_path/'analysis/heldout_metrics.csv')
    primary=result[(result.method==evaluation.PRIMARY)&(result.alpha==.05)]
    assert (primary.power==1).all() and (primary.neutral_call_fraction==0).all()
    assert (primary.selected_n==100).all() and (primary.neutral_n==1000).all()
    extra=result[result.method.isin(['mass_g50_r1_T50000','excess_g50_r1_T50000'])]
    assert len(extra)>0 and (extra.power==1).all() and (extra.neutral_call_fraction==0).all()
    af=result[result.method=='af_r1']
    assert (af.power==0).all()
    target=pd.read_csv(tmp_path/'analysis/target70_metrics.csv')
    assert (target[target.method=='af_r1'].neutral_call_fraction==1).all()
    checks=json.loads((tmp_path/'analysis/figure_checks.json').read_text())
    assert len(checks)==2 and all(c['vector_pdf'] and c['legend_clear'] for c in checks)
