from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest
import zstandard
from gamma_smc_aou import focused_selection_decode as decode
from scipy.special import gammainc


def _raw_fixture(
    tmp_path,
    *,
    trailing: bool = False,
    mismatched_pairs: bool = False,
    missing_focal: bool = False,
):
    raw = tmp_path / "posterior.zst"
    pairs_path = tmp_path / "pair_table.tsv"
    pair_table = pd.DataFrame(
        {
            "gamma_smc_haplotype_0": [0, 2, 4],
            "gamma_smc_haplotype_1": [1, 3, 5],
            "genotype_class": ["hom_ref", "heterozygous", "hom_alt"],
        }
    )
    pair_table.to_csv(pairs_path, sep="\t", index=False)
    pairs = [[0, 1], [2, 3], [4, 5]]
    if mismatched_pairs:
        pairs[1] = [2, 5]
    metadata = {
        "scaled_mutation_rate": 1.0,
        "sequence_length": 3,
        "chunk_size": 4,
        "num_pairs": 3,
        "output_positions": [0, 6, 10] if missing_focal else [0, 5, 10],
        "pairs": pairs,
    }
    raw.with_name(raw.name + ".meta").write_text(json.dumps(metadata), encoding="utf-8")
    alpha = np.zeros((3, 4), dtype=np.float32)
    beta = np.zeros((3, 4), dtype=np.float32)
    alpha[:, :3] = np.asarray([[1, 2, 3], [2, 3, 4], [3, 4, 5]])
    beta[:, :3] = 1
    payload = alpha.tobytes() + beta.tobytes()
    if trailing:
        payload += b"x"
    raw.write_bytes(zstandard.ZstdCompressor().compress(payload))
    return raw, pairs_path


def _decode(tmp_path, raw, pairs_path, output_name="decoded", **kwargs):
    return decode.postprocess_raw_posteriors(
        raw,
        pairs_path,
        tmp_path / output_name,
        unit_record={
            "unit_id": "eas__af10__selected__rep001",
            "demography_id": "eas_phlash_median",
            "simulation_class": "selected",
            "target_allele_frequency": 0.1,
            "replicate_index": 1,
            "sequence_length_bp": 11,
        },
        mutation_rate=0.5,
        generation_time_years=1,
        thresholds_years=(1, 2),
        focal_position_0based=5,
        **kwargs,
    )


def test_one_pass_decode_preserves_pair_order_and_exact_gamma_probabilities(tmp_path):
    raw, pairs_path = _raw_fixture(tmp_path)
    artifacts = _decode(tmp_path, raw, pairs_path)

    assert not artifacts.cache_hit
    pair_summary = pd.read_csv(artifacts.pair_summaries, sep="\t")
    focal = pair_summary[
        (pair_summary["pair_index"] == 0) & (pair_summary["threshold_years"] == 1)
    ].iloc[0]
    assert focal["focal_p_tmrca_lt_threshold"] == pytest.approx(gammainc(2, 1))
    assert focal["focal_tmrca_generations"] == pytest.approx(2)
    assert focal["region_mean_tmrca_generations"] == pytest.approx(2)

    classes = pd.read_csv(artifacts.class_summaries, sep="\t")
    hom_alt = classes[
        (classes["genotype_class"] == "hom_alt") & (classes["threshold_years"] == 2)
    ].iloc[0]
    assert hom_alt["n_pairs"] == 1
    assert hom_alt["focal_mean_p_tmrca_lt_threshold"] == pytest.approx(gammainc(4, 2))
    overall = classes[
        (classes["genotype_class"] == "overall") & (classes["threshold_years"] == 1)
    ].iloc[0]
    assert overall["n_pairs"] == 3

    spatial = pd.read_csv(artifacts.spatial_profiles, sep="\t")
    assert len(spatial) == 4 * 2 * 3
    assert set(spatial["position_0based"]) == {0, 5, 10}


def test_decode_cache_validates_after_raw_compaction_and_fails_on_corruption(tmp_path):
    raw, pairs_path = _raw_fixture(tmp_path)
    artifacts = _decode(tmp_path, raw, pairs_path)
    cached = _decode(
        tmp_path,
        raw,
        pairs_path,
        remove_raw_after_success=True,
    )
    assert cached.cache_hit
    assert not raw.exists()
    assert not raw.with_name(raw.name + ".meta").exists()
    cached_again = _decode(tmp_path, raw, pairs_path)
    assert cached_again.cache_hit

    artifacts.class_summaries.write_text("corrupt\n", encoding="utf-8")
    with pytest.raises(ValueError, match="checksum"):
        _decode(tmp_path, raw, pairs_path)


def test_decode_rejects_pair_order_mismatch_and_trailing_stream_bytes(tmp_path):
    mismatch_dir = tmp_path / "mismatch"
    mismatch_dir.mkdir()
    raw, pairs_path = _raw_fixture(mismatch_dir, mismatched_pairs=True)
    with pytest.raises(ValueError, match="pair order"):
        _decode(mismatch_dir, raw, pairs_path)

    trailing_dir = tmp_path / "trailing"
    trailing_dir.mkdir()
    raw, pairs_path = _raw_fixture(trailing_dir, trailing=True)
    with pytest.raises(ValueError, match="trailing"):
        _decode(trailing_dir, raw, pairs_path)

    missing_focal_dir = tmp_path / "missing-focal"
    missing_focal_dir.mkdir()
    raw, pairs_path = _raw_fixture(missing_focal_dir, missing_focal=True)
    with pytest.raises(ValueError, match="exact focal position"):
        _decode(missing_focal_dir, raw, pairs_path)


def test_decode_recovers_partial_outputs_from_validated_raw_stream(tmp_path):
    raw, pairs_path = _raw_fixture(tmp_path)
    output = tmp_path / "decoded"
    output.mkdir()
    (output / "class_summaries.tsv").write_text("partial\n", encoding="utf-8")
    artifacts = _decode(tmp_path, raw, pairs_path)

    assert artifacts.completion.is_file()
    classes = pd.read_csv(artifacts.class_summaries, sep="\t")
    assert set(classes["genotype_class"]) == {
        "overall",
        "hom_ref",
        "heterozygous",
        "hom_alt",
    }


def test_decode_manifests_auxiliary_outputs_and_postprocessor_source(tmp_path):
    raw, pairs_path = _raw_fixture(tmp_path)
    auxiliary = {
        "overall_summary": tmp_path / "overall.summary.tsv",
        "run_manifest": tmp_path / "overall.summary.tsv.run.json",
    }
    auxiliary["overall_summary"].write_text("position\n", encoding="utf-8")
    auxiliary["run_manifest"].write_text("{}\n", encoding="utf-8")

    artifacts = _decode(
        tmp_path,
        raw,
        pairs_path,
        auxiliary_output_paths=auxiliary,
    )
    completion = json.loads(artifacts.completion.read_text(encoding="utf-8"))

    assert set(auxiliary).issubset(completion["outputs"])
    assert completion["contract"]["static"]["implementation"][
        "postprocessor_sha256"
    ] == decode.sha256_file(decode.POSTPROCESSOR_SOURCE)
    assert _decode(
        tmp_path,
        raw,
        pairs_path,
        auxiliary_output_paths=auxiliary,
    ).cache_hit
