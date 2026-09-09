import json

import numpy as np
import pandas as pd
import pytest

from gamma_smc_aou.fresh_power import digest
from gamma_smc_aou.recent_call_comparison import (
    cached_decode,
    grid_summary,
    longest_run,
)


def test_run_counts_match_explicit_connected_hits():
    rng = np.random.default_rng(20380101)
    masks = rng.random((100, 37)) < 0.4
    masks[0] = False
    masks[1] = True
    for gap in (0, 1):
        expected = []
        for mask in masks:
            hits = np.flatnonzero(mask)
            groups = np.split(hits, np.flatnonzero(np.diff(hits) > gap + 1) + 1)
            expected.append(max(map(len, groups)))
        np.testing.assert_array_equal(longest_run(masks, gap), expected)


def test_grid_calls_anywhere_and_retains_fixed_selected():
    meta = pd.DataFrame(
        dict(mode=["neutral"] * 40 + ["selected"], fixed=[False] * 40 + [True])
    )
    values = np.zeros((41, 20, 1))
    values[-1, -1, 0] = 1
    cfg = dict(haplotype_pairs=10000, stride_bp=10000, tmrca_cutoffs_years=[10000])
    rows = pd.DataFrame(grid_summary(values, meta, cfg, "mean", 10000))
    rule = rows[
        (rows.alpha == 0.05)
        & (rows.stride_bp == 10000)
        & (rows.max_gap_strides == 0)
        & (rows.min_significant_strides == 1)
    ]
    assert rule.set_index("mode").loc["selected", "called"] == 1
    assert rule.set_index("mode").loc["neutral", "called"] == 0
    assert set(rule.n) == {40, 1}


def test_decode_cache_checks_bytes_and_scientific_identity(tmp_path):
    cfg = dict(
        scored_length_bp=20,
        stride_bp=10,
        tmrca_cutoffs_years=[10000],
        haplotype_pairs=2,
    )
    assert not cached_decode(tmp_path, "a", cfg)
    file = tmp_path / "frac_recent.tsv"
    file.write_text("position_0based\tfrac_recent_10000\n0\t0.5\n10\t1.0\n")
    receipt = dict(
        fingerprint="a",
        outputs={file.name: dict(bytes=file.stat().st_size, sha256=digest(file))},
    )
    (tmp_path / "complete.json").write_text(json.dumps(receipt))
    assert cached_decode(tmp_path, "a", cfg)
    with pytest.raises(ValueError, match="Different decoding inputs"):
        cached_decode(tmp_path, "b", cfg)
    file.write_text(file.read_text().replace("0.5", "0.0"))
    with pytest.raises(ValueError, match="Corrupt saved decoding"):
        cached_decode(tmp_path, "a", cfg)
