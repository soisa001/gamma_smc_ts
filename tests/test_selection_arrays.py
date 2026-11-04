"""Selection episode timing and temporal recovery contracts."""
from pathlib import Path
from types import SimpleNamespace
import sys
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))


@pytest.mark.parametrize('origin,age,first,stop',[('de_novo',1600,2,802),('introgression',2000,402,1202)])
def test_episode_switches_off_before_first_post_episode_reproduction(tmp_path,monkeypatch,origin,age,first,stop):
    import run_selection_arrays as runner
    def onset(cfg,origin,label,directory):
        (directory/'onset.trees').write_bytes(b'test fixture')
        return [1],[10]*(age+1)
    monkeypatch.setattr(runner.ancestor,'onset',onset)
    coefficient=.02 if first==2 else 0.
    stdout=f'TRAJECTORY 1 0.05 10 {coefficient} 1\nTRAJECTORY 2 0 10 0 0\n'
    monkeypatch.setattr(runner.subprocess,'run',lambda *a,**k:SimpleNamespace(stdout=stdout,stderr='',check_returncode=lambda:None))
    cfg=dict(seed_namespace='test:',selection_start_years=40000,selection_end_years=20000,length=1000000,slim='slim')
    rec=runner.simulate(cfg,dict(origin=origin,unit='test',selection=.02,cohort='selected'),tmp_path)
    script=(tmp_path/'model.slim').read_text()
    assert f'{stop} early()' in script and 'muts.setSelectionCoeff(0.0)' in script
    if first>2:assert f'{first} early()' in script
    assert rec['selected_reproduction_ticks']==[first,stop-1]
    assert stop-first==800  # Exactly 20,000 years, then 800 neutral generations.
    assert (age+1-(stop-1))*25==20000
    assert not rec['accepted'] and rec['rejection']=='loss'


def test_temporal_scoring_does_not_credit_boundary_touch_or_remote_signal():
    from report_selection_arrays import evaluate
    calls=np.zeros((100,49),bool);positions=np.arange(100)*10000
    calls[50,18]=True;calls[50,39]=True  # [19,20), [40,41): outside.
    calls[0,25]=True  # Correct time, wrong place.
    r=evaluate(calls,positions,40,20,True)
    assert r['joint_hit']==0 and r['center_mistimed']==1 and r['center_off_target']==1
    calls[50,19:39]=True
    r=evaluate(calls,positions,40,20,True)
    assert r['joint_hit']==1 and r['overlap_kyr']==20 and r['temporal_f1']==40/42
    assert r['center_mistimed']==0  # Outside components joined into one broad call.


def test_ongoing_fifty_kya_has_no_observable_outside_time_cells():
    from report_selection_arrays import evaluate
    calls=np.ones((100,49),bool)
    r=evaluate(calls,np.arange(100)*10000,50,0,True)
    assert r['temporal_f1']==1 and r['center_off_target']==0


def test_ledger_survives_transient_windows_reader_lock(tmp_path,monkeypatch):
    import run_selection_arrays as runner
    import json
    original=runner.base.save;attempts=[]
    def locked_save(path,value):
        attempts.append(path)
        if len(attempts)<3:raise PermissionError('reader holds rename lock')
        original(path,value)
    monkeypatch.setattr(runner.base,'save',locked_save)
    monkeypatch.setattr(runner.time,'sleep',lambda _:None)
    path=tmp_path/'status.json'
    runner.save_ledger(path,{'completed':17})
    assert json.loads(path.read_text())=={'completed':17} and len(attempts)==3
    assert not path.with_suffix('.json.tmp').exists()


@pytest.mark.parametrize('error,expected_attempts',[(PermissionError,51),(ValueError,1)])
def test_ledger_failures_are_bounded_and_not_silenced(tmp_path,monkeypatch,error,expected_attempts):
    import run_selection_arrays as runner
    attempts=[]
    def fail(*args):
        attempts.append(1)
        raise error('persistent failure')
    monkeypatch.setattr(runner.base,'save',fail)
    monkeypatch.setattr(runner.time,'sleep',lambda _:None)
    with pytest.raises(error,match='persistent failure'):
        runner.save_ledger(tmp_path/'status.json',{})
    assert len(attempts)==expected_attempts


@pytest.mark.parametrize('origin',['de_novo','introgression'])
def test_zero_end_preserves_recent_simulation_and_rng(tmp_path,monkeypatch,origin):
    import run_selection_arrays as new
    import run_recent_selection_pilot as old
    age,first=old.timing(origin,10000)
    def onset(cfg,origin,label,directory):
        assert label=='recent-v1:example'
        (directory/'onset.trees').write_bytes(b'test fixture')
        return [1],[10]*(age+1)
    monkeypatch.setattr(new.ancestor,'onset',onset)
    coefficient=.02 if first==2 else 0.
    stdout=f'TRAJECTORY 1 0.05 10 {coefficient} 1\nTRAJECTORY 2 0 10 0 0\n'
    commands=[]
    def run(command,**kwargs):
        commands.append(command)
        return SimpleNamespace(stdout=stdout,stderr='',check_returncode=lambda:None)
    monkeypatch.setattr(new.subprocess,'run',run)
    cfg=dict(seed_namespace='recent-v1:',selection_start_years=10000,selection_end_years=0,length=1000000,slim='slim')
    task=dict(origin=origin,unit='example',selection=.02,cohort='selected')
    old.simulate(cfg,task,tmp_path);script=(tmp_path/'model.slim').read_text();trajectory=(tmp_path/'trajectory.tsv').read_bytes()
    new.simulate(cfg,task,tmp_path)
    assert script.split()==(tmp_path/'model.slim').read_text().split()
    assert trajectory==(tmp_path/'trajectory.tsv').read_bytes() and commands[0]==commands[1]
