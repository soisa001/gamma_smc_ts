import os
import subprocess

import pandas as pd
import pytest

from test_simulation import diploid_ts


@pytest.mark.skipif("GAMMA_SMC_BIN" not in os.environ, reason="compiled Linux binary not available")
def test_cpp_tree_input_streaming_recent_summary(tmp_path):
    ts = diploid_ts(length=20_000)
    source = tmp_path / "input.trees"
    output = tmp_path / "recent.tsv"
    ts.dump(source)
    subprocess.run(
        [
            os.environ["GAMMA_SMC_BIN"], "--input", str(source),
            "--input_format", "trees", "--only_within",
            "--scaled_mutation_rate", "0.0008",
            "--recombination_to_mutation_ratio", "0.05",
            "--unscaled_mutation_rate", "2e-7",
            "--recent_threshold_years", "4500", "--generation_time", "30",
            "--recent_summary", str(output),
            "--cache_size", "10",
            # Per-site output is no longer the default; ask for it explicitly.
            "--output_at_hets", "--output_at_stride", "-1",
        ],
        check=True,
    )
    result = pd.read_csv(output, sep="\t")
    assert len(result) == ts.num_sites
    assert result["n_pairs"].eq(4).all()
    assert result["mean_p_tmrca_lt_threshold"].between(0, 1).all()
