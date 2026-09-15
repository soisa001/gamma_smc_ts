import numpy as np

from gamma_smc_aou.allele_class_ablation import local_scores, score_one, calibrate


def test_remove_only_mixed_pair_penalty_and_retain_fixation():
    n=np.array([[2500,5000,2500,10000],[0,0,10000,10000],[10000,0,0,10000]])
    counts=np.array([[[[500],[100],[2000],[2600]],[[0],[0],[8000],[8000]],[[2000],[0],[0],[2000]]]])
    values=local_scores(counts,n,np.ones(3,dtype=bool))
    np.testing.assert_allclose(values['weighted_contrast'][0,:,0],[.15,.8,0])
    np.testing.assert_allclose(values['joint'][0,:,0],[.147,.8,0])
    np.testing.assert_allclose(values['aa'][0,:,0],[.8,.8,0])
    np.testing.assert_allclose(values['rr'][0,:,0],[.2,0,.2])
    changed=counts.copy();changed[0,0,1,0]=2500;changed[0,0,3,0]=5000
    updated=local_scores(changed,n,np.ones(3,dtype=bool))
    for family in ('aa','rr','contrast','weighted_contrast','mass'):
        np.testing.assert_array_equal(values[family],updated[family])
    np.testing.assert_allclose(updated['joint'][0,0,0],.075)


def test_unassigned_grid_is_not_ref_evidence_and_empty_sites_are_retained():
    counts=np.array([[[[10],[0],[0],[10]],[[1],[0],[5],[6]]]]*2)
    n=np.array([[10,0,0,10],[2,3,5,10]])
    features=dict(grid_counts=counts,grid_n=n,grid_markers=np.array([-1,0]),
                  site_counts=np.zeros((2,0,4,1)),site_n=np.zeros((0,4)),site_af=np.zeros(0))
    # Use no assigned markers for this fixture; all-pair can still see the grid.
    features['grid_markers'][:]=-1
    method_list=[dict(name='rr',family='rr',gate=0,cutoff=0,years=50000),
                 dict(name='all',family='all',gate=0,cutoff=0,years=50000)]
    result,_=score_one((dict(key='fixture',onset_years=0),features,np.zeros((5,2,1)),method_list,[50000]))
    assert not result[:,:,:,0].any()
    assert np.all(result[:,0,:,1]==1) and not result[:,1].any()


def test_saturated_fraction_ties_do_not_become_significant():
    # One recent pair can produce fAA=1 under both hypotheses; no false power
    # from tie-breaking. A separate known signal verifies whole-fold counts.
    items=[dict(key=f'{onset}/{r}',onset_years=onset) for onset,n in [(0,1000),(50000,100)] for r in range(n)]
    folds=np.array([r%5 for n in (1000,100) for r in range(n)])
    scores=np.ones((5,2,2,1100,2))
    scores[:,:,:,:1000,1]=.01;scores[:,:,:,1000:,1]=.5
    method_list=[dict(name='saturated',years=50000),dict(name='separable',years=50000)]
    pred,target=calibrate(scores,items,folds,method_list)
    assert not pred[pred.method=='saturated'].called.any()
    assert pred[(pred.method=='separable')&(pred.onset_years==50000)].called.all()
    assert not pred[(pred.method=='separable')&(pred.onset_years==0)].called.any()
    assert target[target.method=='saturated'].called.all()
    assert not target[(target.method=='separable')&(target.onset_years==0)].called.any()
    assert set(pred[pred.fold==0].key)=={i['key'] for j,i in enumerate(items) if folds[j]==0}


def test_focal_choice_requires_the_selected_snp_and_breaks_neutral_ties_left():
    from pathlib import Path
    import runpy
    import pytest
    choose=runpy.run_path(str(Path(__file__).resolve().parents[1]/'scripts/allele_focal_comparison.py'))['choose_marker']
    assert choose(np.array([100,200]),150,False)==0
    assert choose(np.array([100,200]),200,True)==1
    assert choose(np.array([]),150,False) is None
    with pytest.raises(ValueError):choose(np.array([100,200]),150,True)
    with pytest.raises(ValueError):choose(np.array([]),150,True)
