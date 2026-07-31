from __future__ import annotations

import io
import json
from types import SimpleNamespace
from pathlib import Path

import numpy as np
import pandas as pd
import zstandard

from gamma_smc_aou import workbench_candidates


def _summary(path: Path) -> None:
    positions = np.arange(0, 90_000, 10_000)
    fractions = [0.06, 0.0, 0.0, 0.08, 0.0, 0.0, 0.0, 0.07, 0.0]
    pd.DataFrame({
        "position_0based": positions,
        "position_1based": positions + 1,
        "n_pairs": 100,
        "mean_p_tmrca_lt_threshold": np.asarray(fractions) / 2,
        "mean_tmrca_generations": np.linspace(1000, 2000, len(positions)),
        "frac_recent_4500": fractions,
    }).to_csv(path, sep="\t", index=False)


def test_candidate_regions_merge_only_gaps_up_to_20kb(tmp_path):
    summary = tmp_path / "summary.tsv"
    regions_path = tmp_path / "regions.tsv"
    positions_path = tmp_path / "positions.txt"
    _summary(summary)
    regions = workbench_candidates.build_candidate_regions(
        summary,
        population="afr",
        chromosome=1,
        sequence_length=100_000,
        output_path=regions_path,
        positions_path=positions_path,
        minimum_fraction=0.05,
        merge_gap=20_000,
        stride=10_000,
        profile_half_width=10_000,
    )
    assert len(regions) == 2
    assert regions.iloc[0]["start_0based"] == 0
    assert regions.iloc[0]["end_0based_exclusive"] == 40_000
    assert regions.iloc[0]["peak_position_0based"] == 30_000
    assert regions.iloc[1]["peak_position_0based"] == 70_000
    assert positions_path.read_text(encoding="utf-8").splitlines() == [
        "20000",
        "30000",
        "40000",
        "60000",
        "70000",
        "80000",
    ]


def test_raw_tmrca_array_preserves_literal_pair_manifest_order(tmp_path):
    raw_path = tmp_path / "candidate.zst"
    manifest = tmp_path / "pairs.tsv"
    ordered = tmp_path / "candidate_tmrca.f32.zst"
    ordered_meta = tmp_path / "candidate_tmrca.meta.json"
    working = tmp_path / "working.f32"
    pairs = [[0, 2], [1, 3], [0, 3]]
    manifest.write_text(
        "# manifest\n"
        "0\t2\tA.0\tB.0\n"
        "1\t3\tA.1\tB.1\n"
        "0\t3\tA.0\tB.1\n",
        encoding="utf-8",
    )
    metadata = {
        "scaled_mutation_rate": 0.00075,
        "sequence_length": 2,
        "chunk_size": 8,
        "num_pairs": 3,
        "sample_names": {"0": "A.0", "1": "A.1", "2": "B.0", "3": "B.1"},
        "output_positions": [100, 200],
        "pairs": pairs,
    }
    raw_path.with_name(raw_path.name + ".meta").write_text(
        json.dumps(metadata), encoding="utf-8"
    )
    alpha = np.zeros((2, 8), dtype=np.float32)
    beta = np.ones((2, 8), dtype=np.float32)
    alpha[:, :3] = [[1, 2, 3], [4, 5, 6]]
    payload = np.stack([alpha, beta]).tobytes()
    raw_path.write_bytes(zstandard.ZstdCompressor().compress(payload))

    tmrca, result = workbench_candidates.convert_raw_to_ordered_tmrca(
        raw_path,
        pairs_manifest_path=manifest,
        mutation_rate=1.29e-8,
        output_path=ordered,
        metadata_path=ordered_meta,
        working_path=working,
    )
    two_ne = 0.00075 / (2 * 1.29e-8)
    np.testing.assert_allclose(
        np.asarray(tmrca), np.asarray([[1, 2, 3], [4, 5, 6]]) * two_ne
    )
    assert result["pair_order"] == "column k is data row k of pairs_manifest"
    assert result["shape"] == [2, 3]
    mapping = tmrca._mmap
    del tmrca
    mapping.close()


def test_variant_scoring_uses_haplotypes_from_decoded_pairs(monkeypatch, tmp_path):
    class FakeProcess:
        def __init__(self, *args, **kwargs):
            self.stdout = io.StringIO(
                "chr1\t101\trs1\tA\tG\t50\tPASS\t0|0\t0|0\t1|1\t1|1\n"
            )
            self.stderr = io.StringIO("")

        def wait(self):
            return 0

        def kill(self):
            return None

    monkeypatch.setattr(workbench_candidates.subprocess, "Popen", FakeProcess)
    pairs = np.asarray([[0, 2], [4, 6], [0, 4], [2, 6]], dtype=np.int64)
    tmrca = np.asarray([[100.0, 10.0, 60.0, 70.0]], dtype=np.float32)
    scores = workbench_candidates._score_region_variants(
        bcftools="bcftools",
        bcf_path=tmp_path / "panel.bcf",
        sample_list_path=tmp_path / "samples.txt",
        contig="chr1",
        start_0based=0,
        end_0based_exclusive=200,
        population="AFR",
        chromosome=1,
        region_id="chr1_0001",
        tmrca=tmrca,
        tmrca_positions=np.asarray([100]),
        pairs=pairs,
        n_samples=4,
        minimum_genotype_pairs=1,
    )
    assert len(scores) == 1
    assert scores.iloc[0]["n_hom_ref_pairs"] == 1
    assert scores.iloc[0]["n_hom_alt_pairs"] == 1
    assert scores.iloc[0]["mean_hom_ref_tmrca_generations"] == 100
    assert scores.iloc[0]["mean_hom_alt_tmrca_generations"] == 10
    assert scores.iloc[0]["r2_tmrca"] == 1


def test_candidate_analysis_writes_ranked_profiles_plots_and_ordered_array(
    monkeypatch, tmp_path
):
    summary = tmp_path / "summary.tsv"
    pd.DataFrame({
        "position_0based": [0, 100, 200],
        "position_1based": [1, 101, 201],
        "n_pairs": [5, 5, 5],
        "mean_p_tmrca_lt_threshold": [0.01, 0.08, 0.01],
        "mean_tmrca_generations": [1000.0, 500.0, 1000.0],
        "frac_recent_4500": [0.0, 0.08, 0.0],
    }).to_csv(summary, sep="\t", index=False)
    regions = tmp_path / "regions.tsv"
    positions = tmp_path / "positions.txt"
    workbench_candidates.build_candidate_regions(
        summary,
        population="AFR",
        chromosome=1,
        sequence_length=300,
        output_path=regions,
        positions_path=positions,
        minimum_fraction=0.05,
        merge_gap=20,
        stride=100,
        profile_half_width=100,
    )

    samples = tmp_path / "samples.txt"
    samples.write_text("A\nB\nC\nD\n", encoding="utf-8")
    pairs = [[0, 2], [1, 3], [4, 6], [5, 7], [0, 4]]
    manifest = tmp_path / "pairs.tsv"
    manifest.write_text(
        "\n".join(
            f"{first}\t{second}\t{['A.0', 'A.1', 'B.0', 'B.1', 'C.0', 'C.1', 'D.0', 'D.1'][first]}\t"
            f"{['A.0', 'A.1', 'B.0', 'B.1', 'C.0', 'C.1', 'D.0', 'D.1'][second]}"
            for first, second in pairs
        )
        + "\n",
        encoding="utf-8",
    )
    raw = tmp_path / "candidate.zst"
    raw.with_name(raw.name + ".meta").write_text(
        json.dumps({
            "scaled_mutation_rate": 0.00075,
            "sequence_length": 3,
            "chunk_size": 8,
            "num_pairs": 5,
            "sample_names": {
                str(index): label
                for index, label in enumerate(
                    ["A.0", "A.1", "B.0", "B.1", "C.0", "C.1", "D.0", "D.1"]
                )
            },
            "output_positions": [0, 100, 200],
            "pairs": pairs,
        }),
        encoding="utf-8",
    )
    alpha = np.ones((3, 8), dtype=np.float32)
    beta = np.ones((3, 8), dtype=np.float32)
    alpha[:, :5] = [
        [5.0, 6.0, 1.0, 1.2, 3.0],
        [8.0, 9.0, 0.2, 0.3, 4.0],
        [5.0, 6.0, 1.0, 1.2, 3.0],
    ]
    raw.write_bytes(
        zstandard.ZstdCompressor().compress(np.stack([alpha, beta]).tobytes())
    )

    query_output = "chr1\t101\trs1\tA\tG\t50\tPASS\t0|0\t0|0\t1|1\t1|1\n"

    class FakeProcess:
        def __init__(self, *args, **kwargs):
            self.stdout = io.StringIO(query_output)

        def wait(self):
            return 0

        def kill(self):
            return None

    monkeypatch.setattr(workbench_candidates.subprocess, "Popen", FakeProcess)
    monkeypatch.setattr(
        workbench_candidates.subprocess,
        "run",
        lambda *args, **kwargs: SimpleNamespace(
            stdout="101\tA\tG\t0|0\t0|0\t1|1\t1|1\n",
            stderr="",
            returncode=0,
        ),
    )
    output_dir = tmp_path / "analysis"
    result = workbench_candidates.analyze_candidate_regions(
        summary_path=summary,
        regions_path=regions,
        raw_path=raw,
        pairs_manifest_path=manifest,
        bcf_path=tmp_path / "panel.bcf",
        sample_list_path=samples,
        bcftools="bcftools",
        contig="chr1",
        population="AFR",
        chromosome=1,
        output_dir=output_dir,
        mutation_rate=1.29e-8,
        threshold_years=4500,
        profile_half_width=100,
        variant_half_width=100,
        minimum_genotype_pairs=2,
        minimum_fraction=0.05,
        merge_gap=20,
    )

    assert result["n_regions"] == 1
    assert result["n_regions_with_representative_variant"] == 1
    representative = pd.read_csv(
        output_dir / "representative_variants.tsv", sep="\t"
    ).iloc[0]
    assert representative["position_1based"] == 101
    assert representative["n_hom_ref_pairs"] == 2
    assert representative["n_hom_alt_pairs"] == 2
    profile = pd.read_csv(output_dir / "pair_tmrca_profiles.tsv", sep="\t")
    assert set(profile["pair_class"]) == {"hom_ref", "hom_alt"}
    assert (output_dir / "candidate_tmrca.f32.zst").stat().st_size > 0
    assert (output_dir / "candidate_tmrca.meta.json").stat().st_size > 0
    assert len(list((output_dir / "figures").glob("*.png"))) == 2
    assert len(list((output_dir / "figures").glob("*.pdf"))) == 2
