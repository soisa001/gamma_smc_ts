import json

import msprime
import numpy as np
import pytest
import tskit

from gamma_smc_aou import origin_onset as oo
from gamma_smc_aou import segregating_introgression as si


@pytest.fixture
def cfg():
    return json.loads((oo.REPO/'configs/segregating_introgression.json').read_text())


def test_cohort_and_ascertainment_nulls_are_separate(cfg):
    plan = oo.tasks(cfg)
    assert len(plan) == 4600
    assert sum(t['role']=='selected' for t in plan)==2400
    assert sum(t['role']=='neutral_target' for t in plan)==200
    assert sum(t['role']=='null' for t in plan)==2000
    assert {t['origin'] for t in plan}=={'introgressed'}
    assert {t['family'] for t in plan}=={'I50','I10'}
    assert {t['scaling_factor'] for t in plan}=={1}
    seeds = {t['id']:oo.seed_for(cfg,t) for t in plan}
    assert len(set(seeds.values())) == len(plan)
    for t in si.tasks(dict(cfg,target_replicates=101)):
        if t['id'] in seeds:
            assert seeds[t['id']]==oo.seed_for(cfg,t)
    for family,onset in [('I50',50000),('I10',10000)]:
        assert {t['ascertainment_years'] for t in plan if t['family']==family}=={onset}


def test_demography_has_constant_donor_without_hidden_bottleneck(cfg):
    dem=si.demography(cfg)
    donor_changes=[e for e in dem.events if isinstance(e,msprime.PopulationParametersChange) and e.population=='Neanderthal']
    assert donor_changes==[]
    assert dem['Neanderthal'].initial_size==3600
    pulses=[e for e in dem.events if isinstance(e,msprime.MassMigration)]
    assert [(e.time,e.proportion) for e in pulses]==[(2000,.02),(20000,1)]


def miniature():
    t=tskit.TableCollection(110)
    t.populations.add_row()
    t.populations.add_row()
    for _ in range(4):
        t.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=400,population=0)
    t.nodes.add_row(time=3000,population=1)
    t.nodes.add_row(time=25000,population=0)
    for child in [0,1]:
        t.edges.add_row(0,110,4,child)
        t.migrations.add_row(left=0,right=110,node=child,source=0,dest=1,time=2000)
    for child in [2,3,4]:
        t.edges.add_row(0,110,5,child)
    # Nearest variant is modern-derived. Candidates equidistant from center
    # are resolved by coordinate; their mutations occur above the pulse crossing.
    for pos,age in [(49,2500),(52,2500),(55,1000),(58,2500),(61,2500)]:
        site=t.sites.add_row(pos,'0')
        t.mutations.add_row(site=site,node=0,time=age,derived_state='1')
    t.sort()
    return t.tree_sequence()


def test_nearest_choice_uses_mutation_origin_not_lower_node_population(cfg):
    cfg=dict(cfg,scored_length_bp=100,focal_position_bp=50)
    chosen=si.choose_focal(miniature(),cfg)
    assert chosen['position']==52
    assert chosen['onset_af']==.25
    assert chosen['eligible_window']==[50,60]
    assert chosen['mutation_population']=='Neanderthal'
    assert si.choose_focal(miniature(),dict(cfg,pulse_years=70000)) is None


def test_cache_identity_includes_new_ascertainment(cfg):
    task=oo.tasks(cfg)[0]
    identity=oo.simulation_identity(cfg,task,{})
    assert identity['parameters']['focal_ascertainment']==si.SPEC
    changed=oo.simulation_identity(dict(cfg,focal_ascertainment='other'),task,{})
    assert identity!=changed


def test_crop_moves_observed_focal_to_center(cfg):
    cfg=dict(cfg,scored_length_bp=100,focal_position_bp=50)
    tables=miniature().dump_tables()
    tables.migrations.clear()
    cropped,offset=oo.legacy.crop_at_site(tables.tree_sequence(),52,cfg)
    assert offset==2 and cropped.sequence_length==100
    assert any(site.position==50 for site in cropped.sites())


def test_analysis_uses_onset_specific_control_families(cfg):
    from gamma_smc_aou.origin_onset_analysis import arms,family_for_arm
    assert arms(cfg)==('I50','I10')
    assert family_for_arm('I10',cfg)=='I10'
    assert family_for_arm('I50',cfg)=='I50'


def test_immediate_onset_has_post_pulse_ancestry_at_exact_date(cfg, tmp_path):
    cfg=dict(cfg,simulated_length_bp=110000,scored_length_bp=100000,focal_position_bp=50000)
    task=next(t for t in si.tasks(cfg) if t['family']=='I50' and t['role']=='null')
    focal=si.initial_state(cfg,task,oo.seed_for(cfg,task),tmp_path)
    assert focal is not None
    original=tskit.load(tmp_path/'ascertainment.trees')
    assert np.all(original.nodes_time[original.samples()]==2000)
    origins=original.nodes_population[original.samples()]
    assert 0 < np.sum(origins==1) < len(origins)
    assert si.choose_focal(original,cfg)=={k:v for k,v in focal.items() if k!='slim_mutation_id'}
    onset=tskit.load(tmp_path/'onset.trees')
    assert np.all(onset.nodes_time[onset.samples()]==0)
    assert np.all(onset.nodes_population[onset.samples()]==0)
    for individual in onset.individuals():
        assert len(individual.nodes)==2
        assert [onset.node(n).metadata['slim_id'] for n in individual.nodes]==[2*individual.metadata['pedigree_id']+j for j in range(2)]
    mutation=onset.mutation(focal['mutation_id'])
    assert mutation.metadata['mutation_list'][0]['mutation_type']==2
    assert onset.at(focal['position']).num_samples(mutation.node)==focal['onset_alt_copies']


def test_retry_resume_keeps_completed_attempt_seeds(cfg,tmp_path,monkeypatch):
    cfg=dict(cfg,maximum_sample_attempts=2)
    task=si.tasks(cfg)[0]
    identity=oo.simulation_identity(cfg,task,{})
    oo.legacy.atomic_json(tmp_path/'attempt_identity.json',identity)
    oo.legacy.atomic_json(tmp_path/'attempts.json',[dict(attempt=0,seed=oo.seed_for(cfg,task,0),
        outcome='no_eligible_onset_allele',seconds=1)])
    observed=[]
    def no_candidate(cfg,task,seed,directory):
        observed.append(seed)
        return None
    monkeypatch.setattr(si,'initial_state',no_candidate)
    with pytest.raises(RuntimeError,match='attempt limit'):
        si.simulate(cfg,task,tmp_path,identity,'unused')
    assert observed==[oo.seed_for(cfg,task,1)]
    attempts=json.loads((tmp_path/'attempts.json').read_text())
    assert [r['attempt'] for r in attempts]==[0,1]
