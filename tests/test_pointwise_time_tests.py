"""Tie handling and finite-null resolution for selected-site comparisons."""
import importlib.util
from pathlib import Path
import numpy as np
import pytest

path=Path(__file__).resolve().parents[1]/'scripts/pointwise_time_tests.py'
spec=importlib.util.spec_from_file_location('pointwise_time_tests',path)
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_discrete_bin_ties_survive_different_cumulative_baselines():
    selected=np.array([[[.19,.20],[.19,.20],[.19,.20]]])
    neutral=np.array([[[.01,.02],[.01,.02],[.01,.02]]])
    cumulative,interval=module.exceedance_counts(selected,neutral)
    np.testing.assert_array_equal(cumulative,0)
    np.testing.assert_array_equal(interval[0,:2,0],[1,1])
    # Continuous mass keeps its original double-precision ordering.
    assert interval[0,2,0]==0


def test_identical_simulations_are_ties_not_excess():
    selected=np.array([[[.01,.02,.05],[.02,.02,.03],[.0137,.0245,.0567]]])
    neutral=np.repeat(selected,100,axis=0)
    for counts in module.exceedance_counts(selected,neutral):
        np.testing.assert_array_equal(counts,100)
        np.testing.assert_array_equal((1+counts)/101,1)


def test_strict_thresholds_with_100_neutrals():
    selected=np.tile([0,.2],(1,3,1))
    for exceeds,called005,called001 in [(0,True,True),(4,True,False),(5,False,False)]:
        neutral=np.tile([0,.1],(100,3,1))
        neutral[:exceeds,:,1]=.3
        cumulative,interval=module.exceedance_counts(selected,neutral)
        np.testing.assert_array_equal(cumulative[:,:,1],exceeds)
        np.testing.assert_array_equal(interval,exceeds)
        p=(1+exceeds)/101
        assert (p<.05)==called005 and (p<.01)==called001


def test_discrete_methods_reject_noninteger_pair_fractions():
    selected=np.array([[[0,.0105],[0,.01],[0,.0105]]])
    with pytest.raises(AssertionError):
        module.exceedance_counts(selected,selected)
