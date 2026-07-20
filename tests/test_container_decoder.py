import json
from pathlib import Path

import numpy as np
import zstandard
from scipy.special import gammainc

from gamma_smc_aou.container_decoder import (
    build_container_command,
    summarize_posteriors,
)


def test_docker_command_uses_official_stride_interface(tmp_path):
    command = build_container_command(
        runtime="docker",
        image="docker.io/regevsch/gamma_smc:v0.2",
        work_dir=tmp_path,
        input_name="input.vcf.gz",
        output_name="output.zst",
        scaled_mutation_rate=0.0005,
        recombination_to_mutation_ratio=0.8,
        stride=1000,
    )
    assert command[:3] == ["docker", "run", "--rm"]
    assert "--only_within" in command
    assert command[-2:] == ["--output_at_stride", "1000"]
    assert "--output_at_hets=false" in command


def test_singularity_command_adds_docker_transport(tmp_path):
    command = build_container_command(
        runtime="singularity",
        image="docker.io/regevsch/gamma_smc:v0.2",
        work_dir=tmp_path,
        input_name="in.vcf.gz",
        output_name="out.zst",
        scaled_mutation_rate=0.0005,
        recombination_to_mutation_ratio=0.8,
        stride=1000,
    )
    assert command[0:2] == ["singularity", "run"]
    assert "docker://docker.io/regevsch/gamma_smc:v0.2" in command


def test_streaming_posterior_summary_handles_partial_last_chunk(tmp_path):
    path = tmp_path / "posterior.zst"
    positions = [0, 1000]
    chunk_size = 4
    n_pairs = 5
    chunks = []
    expected_probabilities = []
    expected_tmrca = []
    for chunk_index, real_pairs in enumerate((4, 1)):
        alpha = np.full((2, chunk_size), 2.0 + chunk_index, dtype=np.float32)
        beta = np.full((2, chunk_size), 4.0, dtype=np.float32)
        chunks.extend((alpha.tobytes(), beta.tobytes()))
        expected_probabilities.extend(
            [gammainc(2.0 + chunk_index, 4.0 * 0.009)] * real_pairs
        )
        expected_tmrca.extend(
            [(2.0 + chunk_index) / 4.0 * 20_000] * real_pairs
        )
    path.write_bytes(zstandard.ZstdCompressor().compress(b"".join(chunks)))
    path.with_name(path.name + ".meta").write_text(
        json.dumps({
            "num_pairs": n_pairs,
            "chunk_size": chunk_size,
            "sequence_length": len(positions),
            "output_positions": positions,
            "scaled_mutation_rate": 0.0005,
            "pairs": [[2 * i, 2 * i + 1] for i in range(n_pairs)],
        }),
        encoding="utf-8",
    )
    profile = summarize_posteriors(
        path,
        mutation_rate=1.25e-8,
        threshold_years=4500,
        generation_time=25,
    )
    assert profile["n_pairs"].tolist() == [5, 5]
    assert np.allclose(
        profile["mean_p_tmrca_lt_threshold"], np.mean(expected_probabilities)
    )
    assert np.allclose(
        profile["mean_tmrca_generations"], np.mean(expected_tmrca)
    )
