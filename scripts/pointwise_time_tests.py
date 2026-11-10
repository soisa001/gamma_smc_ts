"""Empirical comparisons that preserve exact ties in discrete pair fractions."""
import numpy as np


def exceedance_counts(selected, neutral, pair_count=1000):
    """Return neutral >= selected counts for CDFs and adjacent-bin fractions.

    Arrays are simulations x (truth, posterior median, posterior mass) x times.
    Truth and median are integer pair counts / pair_count. Subtract their integer
    counts before comparing bins, avoiding floating subtraction breaking ties.
    Posterior-weighted mass retains its original continuous values.
    """
    selected=np.asarray(selected)
    neutral=np.asarray(neutral)
    assert selected.ndim==neutral.ndim==3 and selected.shape[1:]==neutral.shape[1:]
    assert selected.shape[1]==3 and pair_count>0
    assert np.isfinite(selected).all() and np.isfinite(neutral).all()
    cumulative=(neutral[None]>=selected[:,None]).sum(1)
    interval=(np.diff(neutral,axis=-1)[None]>=np.diff(selected,axis=-1)[:,None]).sum(1)
    selected_counts=np.rint(selected[:,:2]*pair_count).astype(np.int64)
    neutral_counts=np.rint(neutral[:,:2]*pair_count).astype(np.int64)
    np.testing.assert_allclose(selected[:,:2]*pair_count,selected_counts,rtol=0,atol=1e-8)
    np.testing.assert_allclose(neutral[:,:2]*pair_count,neutral_counts,rtol=0,atol=1e-8)
    cumulative[:,:2]=(neutral_counts[None]>=selected_counts[:,None]).sum(1)
    interval[:,:2]=(np.diff(neutral_counts,axis=-1)[None]>=np.diff(selected_counts,axis=-1)[:,None]).sum(1)
    return cumulative,interval
