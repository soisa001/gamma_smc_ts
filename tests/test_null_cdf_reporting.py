"""Behavioral checks for the user-specified inclusive, unadjusted empirical tail."""
import importlib.util
from pathlib import Path
import numpy as np
import pytest

spec=importlib.util.spec_from_file_location('null_cdf_reporting',Path(__file__).parents[1]/'scripts/plot_population_null_cdfs.py')
report=importlib.util.module_from_spec(spec)
spec.loader.exec_module(report)


def test_equal_nulls_give_one_and_above_all_gives_zero():
    values=np.full(100,.4)
    assert report.empirical_p(values,.4)==1
    assert report.empirical_p(values,.5)==0
    assert report.empirical_p(np.ones(100),1)==1
    assert np.isnan(report.empirical_p(values,np.nan))


@pytest.mark.parametrize('alpha,rank',[(.05,95),(.01,99)])
@pytest.mark.parametrize('values',[np.arange(1,101)/100,np.r_[np.zeros(95),np.ones(5)],
                                 np.ones(100),np.r_[np.full(25,np.nan),np.linspace(0,1,75)]])
def test_percentile_boundary_matches_direct_tail_with_ties_and_no_calls(alpha,rank,values):
    expected=np.sort(np.where(np.isnan(values),-np.inf,values))[rank-1]
    result=report.empirical_summary(values,alpha)
    assert result['empirical_percentile_pct']==expected*100
    assert result['empirical_p_boundary_pct']==expected*100
    finite=values[np.isfinite(values)]
    candidates=np.unique(np.r_[finite,np.nextafter(finite,np.inf),0,1])
    for target in candidates:
        direct=sum(v>=target for v in values)/100
        assert report.empirical_p(values,target)==direct
        assert (direct<=alpha)==(target>expected)
    assert result['tied_fraction']==sum(v==expected for v in values)/100
    assert result['p_at_score_one']==sum(v>=1 for v in values)/100


def test_null_no_calls_stay_in_both_cdf_and_tail_denominator():
    values=np.r_[np.full(25,np.nan),np.full(75,.5)]
    x,y=report.ecdf(values)
    assert x[0]==0 and y[0]==25
    assert y[-1]==100
    assert report.empirical_p(values,.5)==.75
    assert report.empirical_summary(values,.05)['p_at_cutoff']==.75
