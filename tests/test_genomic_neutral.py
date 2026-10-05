import importlib.util
from pathlib import Path
import numpy as np
from gamma_smc_aou import segregating_introgression as si
from gamma_smc_aou import origin_onset as oo


def load_script(name):
    spec=importlib.util.spec_from_file_location(name,oo.REPO/'scripts'/f'{name}.py')
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_mapped_region_preserves_hotspot_and_separate_mutation_exposure(tmp_path):
    prep=load_script('prepare_genomic_null_regions')
    maps,metrics=prep.region_arrays([(100,140,1),(140,150,12),(160,200,2)],[(125,130)],100,100,1.29e-8)
    np.testing.assert_array_equal(maps['recombination_position'],[0,40,50,60,100])
    np.testing.assert_allclose(maps['recombination_rate'],[1e-8,12e-8,1e-8,2e-8])
    assert metrics['map_coverage']==.9
    assert metrics['callable_fraction']==.85
    target=tmp_path/'map.npz'
    np.savez(target,**maps)
    cfg={'genomic_region':{'path':str(target),'sha256':oo.legacy.digest(target)}}
    mu=si.rate_map(cfg,'mutation',protected_site=42)
    for position in (26,42,55):
        assert mu.get_rate(position)==0
    assert mu.get_rate(43)==1.29e-8
    np.testing.assert_allclose(si.rate_map(cfg,'recombination').get_rate(43),12e-8,rtol=1e-15,atol=0)
    cropped=mu.slice(left=20,right=80,trim=True)
    assert cropped.sequence_length==60 and cropped.get_rate(6)==0


def test_draws_keep_correlations_and_do_not_clip():
    from types import SimpleNamespace
    module=load_script('draw_population_histories')
    artifact=SimpleNamespace(mean_log_ne=np.log([.5,2.,2e8]),covariance_factor=np.array([[.01,.01,-.01]]),jitter=0)
    a=module.draw(artifact,100,1729)
    b=module.draw(artifact,100,1729)
    np.testing.assert_array_equal(a,b)
    np.testing.assert_allclose(a[:,1]/a[:,0],4)
    np.testing.assert_allclose(a[:,0]*a[:,2],1e8)
    assert a[:,0].max()<1 and a[:,2].min()>1e8


def test_history_uses_original_time_grid_and_constant_recent_endpoint(tmp_path):
    target=tmp_path/'draw.npz'
    np.savez(target,time_generations=[100.,200.,1000.],ne=[700.,900.,1100.])
    cfg={'demographic_draw':dict(path=str(target),sha256=oo.legacy.digest(target))}
    np.testing.assert_array_equal(si.population_sizes(cfg,[0,99,100,199,200,1500]),[700,700,700,700,900,1100])


def test_history_cdfs_equal_weight_groups_and_include_null_no_calls():
    from gamma_smc_aou.genomic_neutral_report import history_cdfs
    grid=np.array([0,.2,.5,1])
    result=history_cdfs([np.array([.2,.5]),np.array([np.nan,1])],grid)
    np.testing.assert_array_equal(result,[[0,.5,1,1],[.5,.5,.5,1]])
    np.testing.assert_array_equal(result.mean(axis=0),[.25,.5,.75,1])
