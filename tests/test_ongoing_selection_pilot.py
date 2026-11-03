"""Check ongoing-period scoring without claiming identifiable dates."""
import importlib.util
from pathlib import Path
import sys
import numpy as np

scripts=Path(__file__).resolve().parents[1]/'scripts'
sys.path.insert(0,str(scripts))
spec=importlib.util.spec_from_file_location('ongoing_report',scripts/'report_ongoing_selection_pilot.py')
report=importlib.util.module_from_spec(spec)
spec.loader.exec_module(report)


def test_ongoing_temporal_coverage_is_not_short_episode_overlap():
    calls=np.zeros((100,49),bool);positions=np.arange(100)*10000
    assert report.evaluate(calls,positions,True)['temporal_f1']==0
    calls[50,0]=True
    result=report.evaluate(calls,positions,True)
    assert result['spatial_hit']==1 and result['temporal_f1']==2/50
    calls[50,:]=True
    assert report.evaluate(calls,positions,True)['temporal_f1']==1


def test_distant_signal_does_not_count_as_selected_locus_recovery():
    calls=np.zeros((100,49),bool);positions=np.arange(100)*10000
    calls[0,0]=True
    result=report.evaluate(calls,positions,True)
    assert result['any_call']==1 and result['spatial_hit']==0
    assert result['event_fp']==1 and result['event_fn']==1
    calls[50,30]=True
    result=report.evaluate(calls,positions,True)
    assert result['event_tp']==1 and result['event_fp']==1


def test_no_survivors_gives_undefined_conditional_power():
    assert report.rate([])==dict(k=0,n=0,rate=None,ci95_low=None,ci95_high=None)
    result=report.rate([False]*20)
    assert result['rate']==0 and result['ci95_high']>.1


def test_strict_segregation_excludes_both_absorbing_states():
    from report_conditioned_ongoing_pilot import is_segregating
    assert not is_segregating(0)
    assert not is_segregating(1)
    assert is_segregating(1/7264)
    assert is_segregating(1-1/7264)
