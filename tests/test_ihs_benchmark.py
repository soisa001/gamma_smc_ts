"""Independent EHH integration and calibration-summary checks for iHS."""
import importlib.util
from pathlib import Path

import numpy as np

spec = importlib.util.spec_from_file_location("benchmark_ihs", Path(__file__).resolve().parents[1] / "scripts/benchmark_ihs.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def manual_ihh(h, pos, core, allele, cutoff=.05):
    panel = h[:, h[core] == allele]
    denominator = panel.shape[1]*(panel.shape[1]-1)
    total = 0.
    for direction in (-1, 1):
        previous = 1.
        k = core+direction
        while 0 <= k < len(h):
            interval = panel[min(core, k):max(core, k)+1].T
            counts = np.unique(interval, axis=0, return_counts=True)[1]
            ehh = np.sum(counts*(counts-1))/denominator
            total += abs(pos[k]-pos[k-direction])*(previous+ehh)/2
            if ehh <= cutoff:
                break
            previous = ehh
            k += direction
        else:
            return np.nan
    return total


def test_library_matches_direct_haplotype_partition_integration():
    rng = np.random.default_rng(1729)
    h = rng.integers(0, 2, size=(61, 40), dtype=np.int8)
    core = 30
    h[core, :20] = 0
    h[core, 20:] = 1
    h[core-4:core+5, 20:] = 1  # Longer shared ALT background.
    pos = np.arange(len(h))*1000+1
    cfg = dict(recombination_rate=1e-8, min_ehh=.05, min_maf=.05,
               include_edges=False, gap_scale=None, max_gap=None)
    result = module.raw_ihs(h, pos, cfg)
    expected = np.log(manual_ihh(h, pos, core, 1)/manual_ihh(h, pos, core, 0))
    np.testing.assert_allclose(result[core], expected, atol=1e-12, rtol=0)
    assert result[core] > 0
    scaled = module.raw_ihs(h, pos, dict(cfg, recombination_rate=1e-7))
    np.testing.assert_allclose(result, scaled, atol=1e-12, rtol=0, equal_nan=True)
    assert np.isnan(result[0]) and np.isnan(result[-1])


def test_frequency_moments_and_exact_count_bin_boundaries():
    counts = np.array([20, 21, 200, 380])
    values = np.array([1., 3., 8., 2.])
    result = module.moments(counts, values, 20)
    np.testing.assert_array_equal(module.bin_ids(counts, 20), [1, 1, 10, 19])
    np.testing.assert_array_equal(result[:, 1], [2, 4, 10])
    assert result[0].sum() == 4 and result[1].sum() == 14


def test_region_window_score_retains_empty_or_undefined_focal():
    cfg = dict(scored_length_bp=300, window_bp=100, window_stride_bp=50,
               focal_position_bp=150, minimum_window_sites=2, extreme_abs_ihs=2, core_radius_bp=5)
    pos = np.array([10, 20, 110, 120, 160, 170, 210, 220])
    z = np.array([2.5, 2.5, 0, 0, 3, 3, 0, 0])
    region, focal = module.score_region(pos, z, cfg)
    np.testing.assert_array_equal(region, [3, 1])
    np.testing.assert_array_equal(focal, [0, .5])
    assert np.all(focal <= region)
    _, nearby = module.score_region(pos, z, dict(cfg, core_radius_bp=10))
    np.testing.assert_array_equal(nearby, [3, .5])
    region, focal = module.score_region(np.array([]), np.array([]), cfg)
    np.testing.assert_array_equal(region, [0, 0])
    np.testing.assert_array_equal(focal, [0, 0])
