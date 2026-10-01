import json
import numpy as np
import pytest
from gamma_smc_aou import population_neutral as pn
from gamma_smc_aou import origin_onset as oo
from gamma_smc_aou import segregating_introgression as si
from gamma_smc_aou.origin_onset_analysis import rank_p
from gamma_smc_aou.eas_sweep_models import load_phlash_eas_npz


def study():
    return json.loads((oo.REPO/'configs/all_population_neutral.json').read_text(encoding='utf-8-sig'))


def test_grid_is_only_neutral_50k_and_seeds_preserve_eas():
    configs = pn.configs(study())
    jobs = pn.jobs(configs)
    assert len(jobs)==6000
    assert all(t['s']==0 and t['role']=='null' and t['ascertainment_years']==50000 for _,t in jobs)
    seeds=[oo.seed_for(configs[p],t) for p,t in jobs]
    assert len(set(seeds))==6000
    old=json.loads((oo.REPO/'configs/segregating_introgression.json').read_text())
    for p,t in jobs:
        if p=='EAS':
            assert oo.seed_for(configs[p],t)==oo.seed_for(old,t)
            assert oo.simulation_identity(configs[p],t,{})==oo.simulation_identity(old,t,{})


@pytest.mark.parametrize('population',pn.POPULATIONS)
def test_population_demographies_use_correct_history_and_pulse(population):
    cfg=pn.configs(study())[population]
    artifact=si.load_phlash_npz(oo.REPO/cfg['phlash_resource'],expected_sha256=cfg['phlash_sha256'],expected_population=population)
    dem=si.demography(cfg)
    assert dem.populations[0].name==population
    assert dem[population].initial_size==artifact.median_ne[0]
    assert dem['Neanderthal'].initial_size==3600
    pulses=[e for e in dem.events if hasattr(e,'proportion')]
    assert [(e.time,e.source,e.dest,e.proportion) for e in pulses]==[
        (2000,population,'Neanderthal',.02),(20000,'Neanderthal',population,1)]
    if population!='EAS':
        with pytest.raises(ValueError,match='artifact population'):
            load_phlash_eas_npz(oo.REPO/cfg['phlash_resource'],expected_sha256=cfg['phlash_sha256'])


@pytest.mark.parametrize('n',[18,19,97,1000])
def test_critical_cutoff_matches_rank_p_including_ties_and_no_calls(n):
    null=np.round(np.linspace(0,1,n),1)
    null[0]=np.nan
    result=pn.critical_summary(null,.05)
    scores=np.r_[np.linspace(0,1,401), np.nextafter(null[1:],np.inf)]
    assert result['null_n']==n and result['unavailable']==1
    np.testing.assert_array_equal(scores>result['critical_score'], rank_p(null,scores)<=.05)
    saturated=pn.critical_summary(np.ones(n),.05)
    assert not saturated['rejection_possible_for_bounded_score']


def test_cached_focal_scores_can_be_reused_from_larger_profile(tmp_path):
    cfg=pn.configs(study())['EAS']
    row=dict(task_id='I50/null/rep0000',simulation_sha256='sim')
    expected=pn.score_identity(cfg,row,'pairs','truth','decoder')
    archive=tmp_path/'archive'
    old=archive/'regions'/row['task_id']
    old.mkdir(parents=True)
    import pandas as pd
    frame=pd.DataFrame([dict(source='truth',method=m,position_0based=p,cutoff_years=t,n_pairs=1,score=.5)
        for m in ('all_pairs','alt_alt') for p in [0,5000000] for t in cfg['tmrca_cutoffs_years']])
    frame.to_csv(old/'truth.csv',index=False)
    oo.legacy.atomic_json(old/'truth.json',dict(identity=dict(expected,positions=[0,5000000]),scores_sha256=oo.legacy.digest(old/'truth.csv')))
    dest=tmp_path/'dest'
    assert pn.reuse_scores(archive,dest,cfg,row,'pairs','truth','decoder')
    kept=pd.read_csv(dest/'truth.csv')
    assert len(kept)==10 and set(kept.position_0based)=={5000000}
    (old/'truth.csv').write_text('corrupted')
    with pytest.raises(ValueError,match='Corrupt'):
        pn.reuse_scores(archive,dest,cfg,row,'pairs','truth','decoder')
