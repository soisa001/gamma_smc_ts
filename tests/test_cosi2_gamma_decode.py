from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from gamma_smc_aou.cosi2_gamma_decode import (
    DEFAULT_MUTATION_RATE,
    DEFAULT_PAIR_COUNT,
    DEFAULT_PAIR_PANEL_HAPLOTYPES,
    DEFAULT_PAIR_PANEL_SEED,
    TMRCA_THRESHOLDS_YEARS,
    decode_cosi_unit,
    parse_cosi_ms,
    prepare_cosi_gamma_input,
    reduce_gamma_summary,
)


def _haplotypes(n_sites: int = 3) -> list[str]:
    rows = []
    for index in range(200):
        values = [
            "1" if index % 2 else "0",
            "1" if index % 5 == 0 else "0",
            "1" if index < 40 else "0",
        ]
        rows.append("".join(values[:n_sites]))
    return rows


def _write_ms(
    path: Path,
    *,
    mutation_times: bool = True,
    arg_edges: bool = True,
    trailing: str = "",
) -> Path:
    lines = ["ms 200 1", "cosi_rand 12345", "", "// seed=12345"]
    if arg_edges:
        lines.append("E R 0 200 0 1 0 0 0.25 1")
    lines.extend(
        [
            "segsites: 3",
            "positions: 0.1 0.10000001 0.5",
        ]
    )
    if mutation_times:
        lines.append("muttimes: 1 2 3")
    lines.extend(_haplotypes())
    if trailing:
        lines.append(trailing)
    path.write_text("\n".join(lines) + "\n", encoding="ascii", newline="\n")
    return path


def _summary_frame(
    *,
    sequence_length: int,
    stride: int,
    n_pairs: int = DEFAULT_PAIR_COUNT,
) -> pd.DataFrame:
    positions = np.arange(0, sequence_length, stride, dtype=np.int64)
    frame = pd.DataFrame(
        {
            "position_0based": positions,
            "position_1based": positions + 1,
            "n_pairs": n_pairs,
            "mean_tmrca_generations": 500.0,
        }
    )
    for threshold_index, threshold in enumerate(TMRCA_THRESHOLDS_YEARS):
        counts = np.full(len(frame), 5 + 10 * threshold_index, dtype=np.int64)
        hard = counts / n_pairs
        soft = hard + 0.01 + (positions / sequence_length) * 0.001
        frame[f"n_recent_{threshold}"] = counts
        frame[f"frac_recent_{threshold}"] = hard
        frame[f"mean_p_lt_{threshold}"] = soft
    frame["mean_p_tmrca_lt_threshold"] = frame[f"mean_p_lt_{TMRCA_THRESHOLDS_YEARS[0]}"]
    return frame


def test_streaming_parser_accepts_optional_edges_and_mutation_times(tmp_path):
    parsed = parse_cosi_ms(_write_ms(tmp_path / "with_edges.ms"))
    assert parsed.sample_haplotypes == 200
    assert parsed.header_seed == 12345
    assert parsed.replicate_seed == 12345
    assert parsed.arg_edge_count == 1
    assert parsed.positions == (0.1, 0.10000001, 0.5)
    assert parsed.mutation_times_generations == (1.0, 2.0, 3.0)
    assert len(parsed.haplotypes) == 200

    no_optional = parse_cosi_ms(
        _write_ms(
            tmp_path / "without_optional.ms",
            mutation_times=False,
            arg_edges=False,
        )
    )
    assert no_optional.arg_edge_count == 0
    assert no_optional.mutation_times_generations is None


@pytest.mark.parametrize(
    "transform,match",
    [
        (
            lambda text: text.replace(
                "E R 0 200 0 1 0 0 0.25 1",
                "E X 0 200 0 1 0 0 0.25 1",
            ),
            "malformed ARG edge",
        ),
        (
            lambda text: text.replace(_haplotypes()[0], "01", 1),
            "haplotype 0",
        ),
        (lambda text: text + "// second replicate\n", "unexpected trailing record"),
    ],
)
def test_streaming_parser_rejects_malformed_inputs(tmp_path, transform, match):
    path = _write_ms(tmp_path / "bad.ms")
    path.write_text(transform(path.read_text(encoding="ascii")), encoding="ascii")
    with pytest.raises(ValueError, match=match):
        parse_cosi_ms(path)


def test_conversion_retains_duplicate_positions_and_writes_explicit_pairs(tmp_path):
    source = _write_ms(tmp_path / "input.ms")
    artifacts = prepare_cosi_gamma_input(
        source,
        tmp_path / "converted",
        sequence_length_bp=100,
        focal_position_0based=50,
    )
    assert not artifacts.cache_hit
    records = [
        line.split("\t")
        for line in artifacts.vcf_path.read_text(encoding="ascii").splitlines()
        if line and not line.startswith("#")
    ]
    assert [record[1] for record in records] == ["11", "11", "51"]
    assert len({record[2] for record in records}) == 3
    assert all(len(record) == 109 for record in records)
    assert all("|" in genotype for record in records for genotype in record[9:])
    pairs = [
        line
        for line in artifacts.pairs_path.read_text(encoding="ascii").splitlines()
        if line and not line.startswith("#")
    ]
    assert len(pairs) == DEFAULT_PAIR_COUNT
    metadata = json.loads(artifacts.metadata_path.read_text(encoding="utf-8"))
    panel = metadata["pair_panel_indices"]
    expected_pairs = [
        f"{panel[left]}\t{panel[right]}"
        for left in range(DEFAULT_PAIR_PANEL_HAPLOTYPES)
        for right in range(left + 1, DEFAULT_PAIR_PANEL_HAPLOTYPES)
    ]
    assert pairs == expected_pairs
    assert metadata["pair_panel_seed"] == DEFAULT_PAIR_PANEL_SEED
    assert metadata["pair_panel_haplotypes"] == DEFAULT_PAIR_PANEL_HAPLOTYPES
    assert metadata["pair_count"] == DEFAULT_PAIR_COUNT
    assert metadata["discrete_position_collisions"] == 1
    assert metadata["exact_focal_site_present"] is True
    assert metadata["exact_focal_sample_alt_count"] == 40
    assert metadata["exact_focal_pair_panel_alt_count"] == sum(
        index < 40 for index in panel
    )

    cached = prepare_cosi_gamma_input(
        source,
        tmp_path / "converted",
        sequence_length_bp=100,
        focal_position_0based=50,
    )
    assert cached.cache_hit
    artifacts.pairs_path.write_text("corrupt\n", encoding="ascii")
    with pytest.raises(ValueError, match="checksum failed"):
        prepare_cosi_gamma_input(
            source,
            tmp_path / "converted",
            sequence_length_bp=100,
            focal_position_0based=50,
        )


def test_conversion_maps_closed_right_endpoint_to_last_vcf_base(tmp_path):
    source = _write_ms(tmp_path / "terminal.ms")
    source.write_text(
        source.read_text(encoding="ascii").replace(
            "positions: 0.1 0.10000001 0.5",
            "positions: 0.1 0.10000001 1",
        ),
        encoding="ascii",
    )
    artifacts = prepare_cosi_gamma_input(
        source,
        tmp_path / "terminal_converted",
        sequence_length_bp=100,
        focal_position_0based=50,
    )
    positions = [
        int(line.split("\t")[1])
        for line in artifacts.vcf_path.read_text(encoding="ascii").splitlines()
        if line and not line.startswith("#")
    ]
    assert positions == [11, 11, 100]
    metadata = json.loads(artifacts.metadata_path.read_text(encoding="utf-8"))
    assert metadata["normalized_terminal_position_count"] == 1


def test_conversion_retains_ambiguous_duplicate_focal_records(tmp_path):
    source = _write_ms(tmp_path / "focal_collision.ms")
    source.write_text(
        source.read_text(encoding="ascii").replace(
            "positions: 0.1 0.10000001 0.5",
            "positions: 0.1 0.5 0.5",
        ),
        encoding="ascii",
    )
    artifacts = prepare_cosi_gamma_input(
        source,
        tmp_path / "focal_collision_converted",
        sequence_length_bp=100,
        focal_position_0based=50,
    )
    metadata = json.loads(artifacts.metadata_path.read_text(encoding="utf-8"))
    assert metadata["exact_focal_site_count"] == 2
    assert metadata["exact_focal_sample_alt_count"] is None
    assert metadata["exact_focal_pair_panel_alt_count"] is None


def test_summary_reducer_emits_hard_and_soft_focal_and_local_scores():
    summary = _summary_frame(sequence_length=1_000, stride=10)
    scores = reduce_gamma_summary(
        summary,
        sequence_length_bp=1_000,
        focal_position_0based=500,
        output_stride_bp=10,
        local_half_width_bp=100,
    )
    assert len(scores) == 2 * len(TMRCA_THRESHOLDS_YEARS)
    assert set(scores["window"]) == {"focal", "local_100kb"}
    assert set(scores.loc[scores["window"] == "focal", "n_output_positions"]) == {1}
    assert set(scores.loc[scores["window"] == "local_100kb", "n_output_positions"]) == {
        21
    }
    focal = scores[
        (scores["window"] == "focal")
        & (scores["threshold_years"] == TMRCA_THRESHOLDS_YEARS[0])
    ].iloc[0]
    assert focal["hard_mean_call_fraction"] == pytest.approx(5 / DEFAULT_PAIR_COUNT)
    assert focal["soft_mean_cdf_score"] == pytest.approx(
        5 / DEFAULT_PAIR_COUNT + 0.0105
    )


def test_summary_reducer_accepts_gamma_ten_significant_digit_output(tmp_path):
    summary = _summary_frame(sequence_length=1_000, stride=10)
    path = tmp_path / "gamma.summary.tsv"
    summary.to_csv(path, sep="\t", index=False, float_format="%.10g")
    reread = pd.read_csv(path, sep="\t")
    scores = reduce_gamma_summary(
        reread,
        sequence_length_bp=1_000,
        focal_position_0based=500,
        output_stride_bp=10,
        local_half_width_bp=100,
    )
    assert len(scores) == 2 * len(TMRCA_THRESHOLDS_YEARS)


def test_summary_reducer_rejects_grid_and_probability_contract_changes():
    summary = _summary_frame(sequence_length=1_000, stride=10)
    with pytest.raises(ValueError, match="exact ordered stride grid"):
        reduce_gamma_summary(
            summary.iloc[:-1],
            sequence_length_bp=1_000,
            focal_position_0based=500,
            output_stride_bp=10,
            local_half_width_bp=100,
        )
    bad = summary.copy()
    bad.loc[0, "frac_recent_4500"] = 0.0
    bad.loc[0, "n_recent_4500"] = 0
    with pytest.raises(ValueError, match="not monotone"):
        reduce_gamma_summary(
            bad,
            sequence_length_bp=1_000,
            focal_position_0based=500,
            output_stride_bp=10,
            local_half_width_bp=100,
        )


def test_decode_helper_is_checksum_bound_and_uses_mean_calls(tmp_path):
    source = _write_ms(tmp_path / "input.ms")
    decoder = tmp_path / "gamma_smc"
    decoder.write_bytes(b"fake decoder binary")
    calls = []

    def fake_runner(executable, input_path, output_summary, **kwargs):
        calls.append((executable, input_path, kwargs))
        assert kwargs["input_format"] == "vcf"
        assert kwargs["recent_call"] == "mean"
        assert kwargs["only_within"] is False
        assert kwargs["generation_time"] == 29.0
        pair_rows = [
            line
            for line in Path(kwargs["pairs_file"])
            .read_text(encoding="ascii")
            .splitlines()
            if line and not line.startswith("#")
        ]
        assert len(pair_rows) == DEFAULT_PAIR_COUNT
        assert kwargs["scaled_mutation_rate"] == pytest.approx(
            4 * 65_835 * DEFAULT_MUTATION_RATE
        )
        _summary_frame(sequence_length=10_000_000, stride=10_000).to_csv(
            output_summary, sep="\t", index=False
        )
        raw = Path(kwargs["raw_output"])
        raw.write_bytes(b"posterior")
        Path(str(raw) + ".meta").write_text("{}\n", encoding="utf-8")
        Path(kwargs["pairs_manifest"]).write_text(
            Path(kwargs["pairs_file"]).read_text(encoding="ascii"),
            encoding="ascii",
        )
        return {"command": ["fake"], "decode_seconds": 1.0}

    artifacts = decode_cosi_unit(
        source,
        tmp_path / "unit",
        decoder,
        unit_id="han_selected_0000",
        demography_id="out_of_africa_archaic_admixture_5r19_chb",
        present_ne=65_835,
        generation_time_years=29.0,
        runner=fake_runner,
    )
    assert not artifacts.cache_hit
    assert len(calls) == 1
    scores = pd.read_csv(artifacts.scores_path, sep="\t")
    assert len(scores) == 14
    assert set(scores["generation_time_years"]) == {29.0}
    completion = json.loads(artifacts.completion_path.read_text(encoding="utf-8"))
    assert completion["contract"]["settings"]["recent_call"] == "mean"
    assert completion["contract"]["settings"]["n_pairs"] == DEFAULT_PAIR_COUNT
    assert completion["contract"]["settings"]["pair_panel_haplotypes"] == 100
    assert completion["contract"]["settings"]["pair_panel_seed"] == 1729
    decoded_pairs = [
        line
        for line in artifacts.pair_manifest_path.read_text(
            encoding="ascii"
        ).splitlines()
        if line and not line.startswith("#")
    ]
    assert len(decoded_pairs) == DEFAULT_PAIR_COUNT
    run_manifest = json.loads(artifacts.run_manifest_path.read_text(encoding="utf-8"))
    assert run_manifest["execution_paths_use_ephemeral_atomic_stage"] is True
    assert run_manifest["atomic_stage_promoted_to_decoded_directory"] is True
    assert "result_paths_are_relative_to_decoded_directory" not in run_manifest

    def unexpected_runner(*_args, **_kwargs):
        raise AssertionError("valid decode cache should not rerun")

    cached = decode_cosi_unit(
        source,
        tmp_path / "unit",
        decoder,
        unit_id="han_selected_0000",
        demography_id="out_of_africa_archaic_admixture_5r19_chb",
        present_ne=65_835,
        generation_time_years=29.0,
        runner=unexpected_runner,
    )
    assert cached.cache_hit
