"""Observable contracts for scan calibration and interval recovery metrics."""
import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd

path=Path(__file__).resolve().parents[1]/"scripts/report_time_resolved_pilot.py"
spec=importlib.util.spec_from_file_location("time_pilot_report",path)
report=importlib.util.module_from_spec(spec)
spec.loader.exec_module(report)


def test_empirical_upper_tail_counts_ties_and_has_finite_resolution():
    cal=np.arange(19.)
    np.testing.assert_allclose(report.empirical_p(np.array([-1,18,19]),cal),[1,.1,.05])
    assert report.empirical_p(1.,np.ones(99))==1.


def test_probability_mass_differences_preserve_outside_mass():
    cdf=np.linspace(.1,.6,50)[None,:]
    values,intervals=report.interval_values(cdf,[1])
    assert values.shape==(1,49)
    np.testing.assert_allclose(values.sum(),.5)
    np.testing.assert_array_equal(intervals[[0,-1]],[[1,2],[49,50]])
    assert report.interval_values(cdf,[1,2,5,10])[0].shape==(1,182)


def test_no_call_exact_truth_and_overly_broad_intervals():
    calls=np.zeros(49,dtype=bool)
    assert report.temporal_stats(calls)["temporal_f1"]==0
    calls[19:29]=True
    result=report.temporal_stats(calls)
    assert result["temporal_f1"]==1 and result["center_overlap"]==1
    calls[:]=True
    assert np.isclose(report.temporal_stats(calls)["temporal_f1"],20/59)


def test_joint_match_cannot_mix_distant_spatial_and_temporal_hits():
    positions=np.arange(0,1_000_000,10000)
    calls=np.zeros((100,49),dtype=bool)
    calls[0,20]=True
    calls[50,0]=True
    result=report.evaluate_calls(calls,positions,True)
    assert result["any_call"]==1 and result["spatial_hit"]==1
    assert result["joint_hit"]==0 and result["event_fp"]==2 and result["event_fn"]==1
    calls[50,20]=True
    result=report.evaluate_calls(calls,positions,True)
    assert result["joint_hit"]==1 and result["event_tp"]==1 and result["event_fp"]==2


def test_multiscale_windows_expand_without_counting_overlap_twice():
    calls=np.array([[True,True,False],[False,False,True]])
    intervals=np.array([[20,25],[23,30],[1,3]])
    support=report.expand_support(calls,intervals)
    assert support[0].sum()==10 and support[1].sum()==2
    assert report.temporal_stats(support[0])["temporal_f1"]==1


def test_cumulative_baselines_use_independent_calibration_and_spatial_hits():
    records=pd.DataFrame({"cohort":["development"]*10+["calibration"]*19+["neutral_test","selected"]})
    levels=np.r_[np.linspace(.1,.4,10),np.linspace(.1,.4,19),.2,.9]
    cdfs=np.broadcast_to(levels[:,None,None,None],(31,3,3,50))
    result=report.cumulative_baselines(records,cdfs,np.arange(3)*10000)
    assert (result.loc[result.metric=="neutral_fpr","rate"]==0).all()
    assert (result.loc[result.metric=="selected_any_call","rate"]==1).all()
    assert (result.loc[result.metric=="spatial_tpr","rate"]==1).all()
