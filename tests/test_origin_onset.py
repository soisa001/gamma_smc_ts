import io
import json
from pathlib import Path

import numpy as np
import pytest
import stdpopsim
from stdpopsim import slim_engine

from gamma_smc_aou import origin_onset as oo


@pytest.fixture
def cfg():
    return json.loads((oo.REPO / "configs/origin_onset_pilot.json").read_text())


def test_cohort_sizes_shared_nulls_and_independent_reproducible_seeds(cfg):
    plan = oo.tasks(cfg)
    assert len(plan) == 3510
    assert sum(t['role'] == 'selected' for t in plan) == 480
    assert sum(t['role'] == 'null' for t in plan) == 3000
    assert sum(t['role'] == 'neutral_target' for t in plan) == 30
    assert {t['s'] for t in plan if t['role'] == 'selected'} == {i/1000 for i in range(1,11)} | {.02,.03}
    seeds = {t['id']: oo.seed_for(cfg,t) for t in plan}
    assert len(set(seeds.values())) == len(plan)
    expanded = dict(cfg, target_replicates=100, selection_coefficients=list(reversed(cfg['selection_coefficients'])))
    assert all(oo.seed_for(expanded,t)==seeds[t['id']] for t in oo.tasks(expanded) if t['id'] in seeds)
    assert all(oo.seed_for(cfg,t,1)!=seeds[t['id']] for t in plan)


def test_neutral_targets_and_calibration_share_events_but_not_seeds(cfg):
    plan = oo.tasks(cfg)
    for family in oo.FAMILIES:
        null = next(t for t in plan if t['family']==family and t['role']=='null')
        target = next(t for t in plan if t['family']==family and t['role']=='neutral_target')
        assert [vars(e) for e in oo.events(cfg,null)] == [vars(e) for e in oo.events(cfg,target)]
        assert oo.seed_for(cfg,null) != oo.seed_for(cfg,target)
        assert not any(isinstance(e,stdpopsim.ChangeMutationFitness) for e in oo.events(cfg,null))
        assert any(isinstance(e,stdpopsim.ConditionOnAlleleFrequency) for e in oo.events(cfg,null))


def test_delayed_introgression_keeps_old_pulse_and_denovo_is_born_at_onset(cfg):
    plan = oo.tasks(cfg)
    delayed = next(t for t in plan if t['family']=='I50' and t['onset_years']==10000)
    ev = oo.events(cfg,delayed)
    assert ev[0].population == 'Neanderthal'
    assert ev[1].start_time == 400
    assert ev[-1].start_time == 1995
    denovo = next(t for t in plan if t['family']=='D10' and t['role']=='selected')
    ev = oo.events(cfg,denovo)
    assert ev[0].population == 'EAS' and ev[0].time == 400
    assert ev[1].start_time == 400 and denovo['scaling_factor']==1


def test_script_patch_restores_global_functions_on_error(cfg):
    task = next(t for t in oo.tasks(cfg) if t['family']=='D10')
    old_functions, old_makescript = slim_engine._slim_functions, slim_engine.slim_makescript
    with pytest.raises(RuntimeError):
        with oo.simulation_patch(cfg,task):
            assert 'targets = sample(pop.genomes, 1)' in slim_engine._slim_functions
            assert 'origin_placement_alt_copies' in slim_engine._slim_functions
            assert 'origin_record();' in slim_engine._slim_functions
            raise RuntimeError('sentinel')
    assert slim_engine._slim_functions == old_functions
    assert slim_engine.slim_makescript is old_makescript


def test_cache_checks_identity_and_corruption(tmp_path):
    data=tmp_path/'data.csv'
    data.write_text('value\n1\n')
    identity={'source':'test'}
    oo.legacy.atomic_json(tmp_path/'receipt.json',dict(identity=identity,artifacts=oo.artifact_specs(tmp_path,['data.csv'])))
    assert oo.read_receipt(tmp_path,'receipt.json',identity)
    with pytest.raises(ValueError,match='identity'):
        oo.read_receipt(tmp_path,'receipt.json',{'source':'different'})
    data.write_text('value\n2\n')
    with pytest.raises(ValueError,match='Corrupt'):
        oo.read_receipt(tmp_path,'receipt.json',identity)


def test_focal_calibration_keeps_ties_and_unavailable_pair_scores():
    from gamma_smc_aou.origin_onset_analysis import rank_p
    np.testing.assert_allclose(rank_p([.1,.2,.2,.5],[0,.2,.5,1,np.nan]),[1,.8,.4,.2,1])
    np.testing.assert_allclose(rank_p([np.nan,np.nan,1],[0,1,np.nan]),[.5,.5,1])
    assert rank_p(np.zeros(1000),[1])[0]==1/1001


def test_independent_neutral_targets_are_the_fpr_denominator():
    import pandas as pd
    from gamma_smc_aou.origin_onset_analysis import evaluate
    rows=[]
    for i in range(1000):
        rows.append(dict(task_id=f'null{i}',family='I50',arm='I50',role='null',s=0,
                         source='truth',method='alt_alt',cutoff_years=10000,score=i/1000))
    for role,s in [('selected',.003),('neutral_target',0)]:
        for i in range(10):
            # Exactly one positive, one unavailable, eight negative targets.
            value=1 if i==0 else np.nan if i==1 else 0
            rows.append(dict(task_id=f'{role}{i}',family='I50',arm='I50',role=role,s=s,
                             source='truth',method='alt_alt',cutoff_years=10000,score=value))
    predictions,metrics=evaluate(pd.DataFrame(rows),dict(null_replicates=1000,target_replicates=10,alpha=.05))
    assert len(predictions)==20 and set(predictions.role)=={'selected','neutral_target'}
    assert (metrics.n==10).all() and (metrics.calls==1).all()
    assert (metrics.rate==.1).all() and (metrics.evaluable==9).all()
    bad=pd.DataFrame(rows)
    bad.loc[bad.task_id=='neutral_target0','task_id']='null0'
    with pytest.raises(ValueError,match='calibration'):
        evaluate(bad,dict(null_replicates=1000,target_replicates=10,alpha=.05))
