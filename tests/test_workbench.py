from __future__ import annotations

import copy
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from gamma_smc_aou import workbench


def _write_summary(path: Path, chromosome: int = 1) -> None:
    offset = chromosome * 100
    pd.DataFrame(
        {
            "position_0based": [offset, offset + 10_000, offset + 20_000],
            "position_1based": [offset + 1, offset + 10_001, offset + 20_001],
            "n_pairs": [10, 10, 9],
            "mean_p_tmrca_lt_threshold": [0.1, 0.3, 0.2],
            "mean_tmrca_generations": [900.0, 400.0, 700.0],
            "frac_recent_4500": [0.0, 0.4, 0.2],
        }
    ).to_csv(path, sep="\t", index=False)


def test_zoom_plot_labels_only_peaks_above_two_percent_horizontally(
    tmp_path: Path, monkeypatch
) -> None:
    called_column = "frac_recent_4500"
    genome = pd.DataFrame(
        {
            "population": ["AFR", "AFR", "AFR"],
            "chromosome": [1, 1, 1],
            "genome_position_0based": [0, 10_000, 20_000],
            called_column: [0.005, 0.021, 0.035],
        }
    )
    labels = pd.DataFrame(
        {
            "population": ["AFR", "AFR", "AFR"],
            "chromosome": [1, 1, 1],
            "peak_genome_position_0based": [10_000, 20_000, 29_900],
            "peak_ranking_value": [0.019, 0.021, 0.03],
            "probable_gene": ["BELOW", "ABOVE", "RIGHT"],
        }
    )
    captured = []
    monkeypatch.setattr(workbench, "_atomic_figure", lambda *args, **kwargs: None)
    monkeypatch.setattr(workbench.plt, "close", captured.append)
    workbench._plot_combined_recent_genome(
        {"AFR": genome},
        chromosomes=[1],
        called_column=called_column,
        threshold_years=4500,
        ticks=[10_000],
        boundaries=[0, 30_000],
        signal_fraction=0.05,
        output_stem=tmp_path / "zoom",
        fixed_ymax=0.04,
        hit_labels=labels,
        hit_label_min_fraction=0.02,
    )
    figure = captured[0]
    axis = figure.axes[0]
    annotations = {text.get_text(): text for text in axis.texts}
    assert "ABOVE" in annotations
    assert "BELOW" not in annotations
    label = annotations["ABOVE"]
    assert label.get_rotation() == 0
    assert label.get_fontweight() == "bold"
    assert label.get_fontsize() == pytest.approx(workbench.FIGURE_ANNOTATION_SIZE)
    assert label.get_ha() == "left"
    assert label.get_position() == (8, 0)
    right_label = annotations["RIGHT"]
    assert right_label.get_ha() == "right"
    assert right_label.get_position() == (-8, 0)
    assert axis.get_ylim() == pytest.approx((0.0, 0.04))


def _contract(
    tmp_path: Path,
    *,
    with_mask: bool = False,
    mask_source_semantics: str = "excluded_intervals",
) -> tuple[dict, Path, Path, Path]:
    summary = tmp_path / "chr1.gamma_smc.tsv"
    pairs = tmp_path / "chr1.pairs.tsv"
    input_path = tmp_path / "AFR.chr1.phased.bcf"
    index_path = tmp_path / "AFR.chr1.phased.bcf.csi"
    hardmask_path = tmp_path / "hardmask.bed"
    mask_path = tmp_path / "chr1.callable.bed"
    mask_audit = tmp_path / "chr1.callable.audit.json"
    ancestry_path = tmp_path / "ancestry_preds.tsv"
    qc_path = tmp_path / "flagged_samples.tsv"
    relatedness_path = tmp_path / "relatedness_flagged_samples.tsv"
    bcf_samples_path = tmp_path / "bcf_samples.txt"
    sample_list = tmp_path / "chr1.samples.txt"
    sample_audit = tmp_path / "chr1.samples.audit.json"
    _write_summary(summary)
    pairs.write_text(
        "# manifest\n"
        + "".join(
            f"{2 * index}\t{2 * index + 1}\t{1_000_000 + index}.0\t"
            f"{1_000_000 + index}.1\n"
            for index in range(10)
        ),
        encoding="utf-8",
    )
    sample_ids = [str(1_000_000 + index) for index in range(10)]
    pd.DataFrame(
        {
            "research_id": sample_ids,
            "ancestry_pred_other": ["afr"] * len(sample_ids),
        }
    ).to_csv(ancestry_path, sep="\t", index=False)
    pd.DataFrame({"s": []}).to_csv(qc_path, sep="\t", index=False)
    pd.DataFrame({"sample_id.s": []}).to_csv(
        relatedness_path, sep="\t", index=False
    )
    bcf_samples_path.write_text("\n".join(sample_ids) + "\n", encoding="utf-8")
    workbench.build_workbench_sample_list(
        ancestry_path=ancestry_path,
        qc_exclusions_path=qc_path,
        relatedness_exclusions_path=relatedness_path,
        bcf_samples_path=bcf_samples_path,
        population="AFR",
        output_path=sample_list,
        audit_path=sample_audit,
    )
    if with_mask:
        hardmask_path.write_text("chr1\t100\t200\n", encoding="utf-8")
        workbench.build_workbench_callable_mask(
            hardmask_path=hardmask_path,
            contig="chr1",
            sequence_length=1_000_000,
            output_path=mask_path,
            audit_path=mask_audit,
            source_semantics=mask_source_semantics,
        )
    contract = workbench.build_workbench_contract(
        population="afr",
        chromosome=1,
        input_uri="gs://bucket/gamma_smc/inputs/AFR/AFR.chr1.phased.bcf",
        input_fingerprint="123\t456\tcrc",
        local_input=input_path,
        index_uri="gs://bucket/gamma_smc/inputs/AFR/AFR.chr1.phased.bcf.csi",
        index_fingerprint="124\t45\tindex-crc",
        local_index=index_path,
        output_summary=summary,
        pairs_manifest=pairs,
        sample_list=sample_list,
        sample_audit=sample_audit,
        ancestry_uri="gs://bucket/ancestry_preds.tsv",
        ancestry_fingerprint="200\t300\tancestry-crc",
        local_ancestry=ancestry_path,
        qc_exclusions_uri="gs://bucket/flagged_samples.tsv",
        qc_exclusions_fingerprint="201\t301\tqc-crc",
        local_qc_exclusions=qc_path,
        relatedness_exclusions_uri="gs://bucket/relatedness_flagged_samples.tsv",
        relatedness_exclusions_fingerprint="202\t302\trelated-crc",
        local_relatedness_exclusions=relatedness_path,
        mask_uri=("gs://bucket/gamma_smc/inputs/masks/chr1.callable.bed" if with_mask else None),
        mask_fingerprint=("321\t654\tmask-crc" if with_mask else None),
        local_mask_source=(hardmask_path if with_mask else None),
        local_mask=(mask_path if with_mask else None),
        mask_audit=(mask_audit if with_mask else None),
        mask_source_semantics=(mask_source_semantics if with_mask else None),
        theta=0.00075,
        rho_over_theta=0.8,
        mutation_rate=1.29e-8,
        generation_time=25,
        threshold_years=4500,
        recent_call="mean",
        stride=10_000,
        cache_size=1_000,
        threads=12,
        pair_block=256,
        exp10="accurate",
        backward_alignment="fixed",
        code_commit="abc123",
    )
    command = [
        "gamma_smc",
        "--input",
        str(input_path),
        "--recent_summary",
        str(summary),
        "--pairs_manifest",
        str(pairs),
        "--samples",
        str(sample_list),
        "--scaled_mutation_rate",
        "0.00075",
        "--recombination_to_mutation_ratio",
        "0.8",
        "--unscaled_mutation_rate",
        "1.29e-8",
        "--generation_time",
        "25",
        "--recent_threshold_years",
        "4500",
        "--recent_call",
        "mean",
        "--output_at_stride",
        "10000",
        "--cache_size",
        "1000",
        "--threads",
        "12",
        "--pair_block",
        "256",
        "--exp10",
        "accurate",
        "--backward_alignment",
        "fixed",
        "--output_at_hets=false",
        "--only_within",
    ]
    if with_mask:
        command.extend(["--mask", str(mask_path)])
    run_json = summary.with_suffix(summary.suffix + ".run.json")
    run_json.write_text(
        json.dumps(
            {
                "command": command,
                "decode_seconds": 12.5,
                "stride_bp": 10_000,
                "cache_size_bp": 1_000,
                "n_output_positions": 3,
                "pairs_manifest": str(pairs),
                "n_pairs_recorded": 10,
            }
        ),
        encoding="utf-8",
    )
    return contract, summary, run_json, pairs


def test_completion_contract_detects_changed_output(tmp_path):
    contract, summary, run_json, _ = _contract(tmp_path, with_mask=True)
    completion_path = tmp_path / "chr1.complete.json"
    written = workbench.write_workbench_completion(
        summary_path=summary,
        run_json_path=run_json,
        completion_path=completion_path,
        contract=contract,
    )
    assert written["summary"]["n_output_positions"] == 3
    assert written["contract"]["decoder"]["threads"] == 12
    workbench.validate_workbench_completion(
        summary_path=summary,
        run_json_path=run_json,
        completion_path=completion_path,
        contract=contract,
    )

    summary.write_text(summary.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="summary hash"):
        workbench.validate_workbench_completion(
            summary_path=summary,
            run_json_path=run_json,
            completion_path=completion_path,
            contract=contract,
        )


def test_completion_cache_ignores_git_commit_but_keeps_manifest_hashes(tmp_path):
    contract, summary, run_json, pairs = _contract(tmp_path, with_mask=True)
    completion_path = tmp_path / "chr1.complete.json"
    workbench.write_workbench_completion(
        summary_path=summary,
        run_json_path=run_json,
        completion_path=completion_path,
        contract=contract,
    )
    legacy = json.loads(completion_path.read_text(encoding="utf-8"))
    legacy.pop("cache_contract_sha256")
    completion_path.write_text(json.dumps(legacy), encoding="utf-8")

    requested = copy.deepcopy(contract)
    requested["code_commit"] = "plot-only-change"
    workbench.validate_workbench_completion(
        summary_path=summary,
        run_json_path=run_json,
        completion_path=completion_path,
        contract=requested,
    )

    pairs.write_text(pairs.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="pair-manifest hash"):
        workbench.validate_workbench_completion(
            summary_path=summary,
            run_json_path=run_json,
            completion_path=completion_path,
            contract=requested,
        )


def test_strict_mask_semantics_are_part_of_completion_contract(tmp_path):
    contract, summary, run_json, _ = _contract(
        tmp_path,
        with_mask=True,
        mask_source_semantics="included_intervals",
    )
    completion_path = tmp_path / "chr1.complete.json"
    written = workbench.write_workbench_completion(
        summary_path=summary,
        run_json_path=run_json,
        completion_path=completion_path,
        contract=contract,
    )
    assert written["contract"]["mask"]["source_semantics"] == "included_intervals"
    workbench.validate_workbench_completion(
        summary_path=summary,
        run_json_path=run_json,
        completion_path=completion_path,
        contract=contract,
    )

    incompatible = copy.deepcopy(contract)
    incompatible["mask"]["source_semantics"] = "excluded_intervals_complemented"
    with pytest.raises(ValueError, match="completion contract"):
        workbench.validate_workbench_completion(
            summary_path=summary,
            run_json_path=run_json,
            completion_path=completion_path,
            contract=incompatible,
        )


def test_bitmatrix_contract_accepts_native_float32_theta_rounding(
    tmp_path, monkeypatch
):
    bits = tmp_path / "chr1.recent.bits"
    bits.write_bytes(b"packed")
    bits.with_name(bits.name + ".meta").write_text("{}", encoding="utf-8")
    pairs = tmp_path / "chr1.pairs.tsv"
    pairs.write_text("0\t1\n", encoding="utf-8")
    samples = tmp_path / "chr1.samples.txt"
    samples.write_text("sample1\n", encoding="utf-8")
    theta = 0.00075
    mutation_rate = 1.29e-8
    native_two_ne = float(np.float32(theta)) / (2.0 * mutation_rate)
    metadata = {
        "pairs": [[0, 1]],
        "n_pairs": 1,
        "recent_call": "mean",
        "thresholds_years": [4500.0],
        "two_ne_generations": native_two_ne,
        "sample_names": {"0": "sample1.0", "1": "sample1.1"},
    }
    monkeypatch.setattr(workbench.bitmatrix_reader, "read_meta", lambda _: metadata)
    monkeypatch.setattr(
        workbench.bitmatrix_reader,
        "output_positions",
        lambda _: np.asarray([0], dtype=np.int64),
    )
    contract = {
        "analysis_outputs": {"bitmatrix": str(bits)},
        "pairs_manifest": str(pairs),
        "sample_selection": {"sample_list": str(samples)},
        "decoder": {
            "theta": theta,
            "mutation_rate": mutation_rate,
            "threshold_years": 4500.0,
            "recent_call": "mean",
        },
    }
    summary = pd.DataFrame({"position_0based": [0]})

    metrics = workbench.validate_workbench_bitmatrix(
        contract, summary_frame=summary, required=True
    )
    assert metrics["size_bytes"] == len(b"packed")

    metadata["two_ne_generations"] = native_two_ne * 1.001
    with pytest.raises(ValueError, match="time scale"):
        workbench.validate_workbench_bitmatrix(
            contract, summary_frame=summary, required=True
        )


def test_completion_contract_requires_all_mask_fields(tmp_path):
    with pytest.raises(ValueError, match="supplied together"):
        workbench.build_workbench_contract(
            population="AFR",
            chromosome=1,
            input_uri="gs://bucket/input.bcf",
            input_fingerprint="fingerprint",
            local_input=tmp_path / "input.bcf",
            index_uri="gs://bucket/input.bcf.csi",
            index_fingerprint="index-fingerprint",
            local_index=tmp_path / "input.bcf.csi",
            output_summary=tmp_path / "output.tsv",
            pairs_manifest=tmp_path / "pairs.tsv",
            sample_list=tmp_path / "samples.txt",
            sample_audit=tmp_path / "samples.audit.json",
            ancestry_uri="gs://bucket/ancestry.tsv",
            ancestry_fingerprint="ancestry-fingerprint",
            local_ancestry=tmp_path / "ancestry.tsv",
            qc_exclusions_uri="gs://bucket/qc.tsv",
            qc_exclusions_fingerprint="qc-fingerprint",
            local_qc_exclusions=tmp_path / "qc.tsv",
            relatedness_exclusions_uri="gs://bucket/related.tsv",
            relatedness_exclusions_fingerprint="related-fingerprint",
            local_relatedness_exclusions=tmp_path / "related.tsv",
            mask_uri="gs://bucket/mask.bed",
            mask_fingerprint=None,
            local_mask_source=None,
            local_mask=None,
            mask_audit=None,
            theta=0.00075,
            rho_over_theta=0.8,
            mutation_rate=1.29e-8,
            generation_time=25,
            threshold_years=4500,
            recent_call="mean",
            stride=10_000,
            cache_size=1_000,
            threads=12,
            pair_block=256,
            exp10="accurate",
            backward_alignment="fixed",
            code_commit="abc123",
        )


def test_sample_builder_uses_other_label_exclusions_and_bcf_order(tmp_path):
    ancestry = tmp_path / "ancestry.tsv"
    qc = tmp_path / "qc.tsv"
    relatedness = tmp_path / "relatedness.tsv"
    bcf_samples = tmp_path / "bcf_samples.txt"
    output = tmp_path / "AFR.samples.txt"
    audit_path = tmp_path / "AFR.samples.audit.json"
    pd.DataFrame(
        {
            "research_id": ["1", "2", "3", "4", "5", "6"],
            "ancestry_pred_other": ["afr", "AFR", "oth", "eur", "afr", "afr"],
        }
    ).to_csv(ancestry, sep="\t", index=False)
    pd.DataFrame({"s": ["2"]}).to_csv(qc, sep="\t", index=False)
    pd.DataFrame({"sample_id.s": ["5"]}).to_csv(
        relatedness, sep="\t", index=False
    )
    bcf_samples.write_text("6\n5\n4\n3\n2\n1\n7\n", encoding="utf-8")

    audit = workbench.build_workbench_sample_list(
        ancestry_path=ancestry,
        qc_exclusions_path=qc,
        relatedness_exclusions_path=relatedness,
        bcf_samples_path=bcf_samples,
        population="afr",
        output_path=output,
        audit_path=audit_path,
    )
    assert output.read_text(encoding="utf-8").splitlines() == ["6", "1"]
    assert audit["ancestry_label_column"] == "ancestry_pred_other"
    assert audit["id_columns"] == {
        "ancestry": "research_id",
        "qc_exclusions": "s",
        "relatedness_exclusions": "sample_id.s",
    }
    assert audit["counts"]["selected_samples"] == 2
    assert audit["counts"]["bcf_samples_absent_from_ancestry"] == 1
    assert audit["counts"]["target_bcf_samples_excluded_by_qc"] == 1
    assert audit["counts"]["target_bcf_samples_excluded_by_relatedness"] == 1


def test_sample_builder_rejects_conflicting_ancestry_labels(tmp_path):
    ancestry = tmp_path / "ancestry.tsv"
    qc = tmp_path / "qc.tsv"
    relatedness = tmp_path / "relatedness.tsv"
    bcf_samples = tmp_path / "bcf_samples.txt"
    pd.DataFrame(
        {
            "research_id": ["1", "1"],
            "ancestry_pred_other": ["afr", "eur"],
        }
    ).to_csv(ancestry, sep="\t", index=False)
    pd.DataFrame({"s": []}).to_csv(qc, sep="\t", index=False)
    pd.DataFrame({"sample_id.s": []}).to_csv(
        relatedness, sep="\t", index=False
    )
    bcf_samples.write_text("1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="conflicting labels"):
        workbench.build_workbench_sample_list(
            ancestry_path=ancestry,
            qc_exclusions_path=qc,
            relatedness_exclusions_path=relatedness,
            bcf_samples_path=bcf_samples,
            population="AFR",
            output_path=tmp_path / "samples.txt",
            audit_path=tmp_path / "audit.json",
        )


def test_hardmask_is_complemented_into_positive_callable_intervals(tmp_path):
    hardmask = tmp_path / "hardmask.bed"
    callable_mask = tmp_path / "chr1.callable.bed"
    audit_path = tmp_path / "chr1.callable.audit.json"
    hardmask.write_text(
        "# excluded intervals\n"
        "1\t10\t20\n"
        "chr1\t15\t30\n"
        "chr1\t90\t120\n"
        "chr2\t0\t100\n",
        encoding="utf-8",
    )
    audit = workbench.build_workbench_callable_mask(
        hardmask_path=hardmask,
        contig="chr1",
        sequence_length=100,
        output_path=callable_mask,
        audit_path=audit_path,
    )
    assert callable_mask.read_text(encoding="utf-8").splitlines() == [
        "chr1\t0\t10",
        "chr1\t30\t90",
    ]
    assert audit["source_semantics"] == "excluded_intervals"
    assert audit["decoder_mask_semantics"] == "included_intervals"
    assert audit["counts"]["excluded_bases"] == 30
    assert audit["counts"]["callable_bases"] == 70


def test_strict_callable_mask_is_merged_without_complementing(tmp_path):
    strict = tmp_path / "strict.bed"
    callable_mask = tmp_path / "chr1.callable.bed"
    audit_path = tmp_path / "chr1.callable.audit.json"
    strict.write_text(
        "1\t10\t20\nchr1\t15\t30\nchr1\t90\t120\nchr2\t0\t100\n",
        encoding="utf-8",
    )
    audit = workbench.build_workbench_callable_mask(
        hardmask_path=strict,
        contig="chr1",
        sequence_length=100,
        output_path=callable_mask,
        audit_path=audit_path,
        source_semantics="included_intervals",
    )
    assert callable_mask.read_text(encoding="utf-8").splitlines() == [
        "chr1\t10\t30",
        "chr1\t90\t100",
    ]
    assert audit["source_semantics"] == "included_intervals"
    assert audit["counts"]["callable_bases"] == 30
    assert audit["counts"]["excluded_bases"] == 70


def test_population_plots_are_separate_and_whole_genome_is_complete(tmp_path, monkeypatch):
    one_summary = tmp_path / "one.tsv"
    _write_summary(one_summary)
    one_output = tmp_path / "one_plot"
    result = workbench.plot_workbench_population(
        {1: one_summary},
        population="afr",
        output_dir=one_output,
    )
    assert result["population"] == "AFR"
    assert (one_output / "chromosomes" / "chr1.gamma_smc.png").is_file()
    assert (one_output / "chromosome_scan_summary.tsv").is_file()
    assert not (one_output / "AFR.whole_genome.gamma_smc.png").exists()

    summaries = {}
    for chromosome in workbench.AUTOSOMES:
        path = tmp_path / f"chr{chromosome}.tsv"
        _write_summary(path, chromosome)
        summaries[chromosome] = path
    def lightweight_chromosome_plot(*args, output_stem, **kwargs):
        Path(f"{output_stem}.png").write_bytes(b"png")
        Path(f"{output_stem}.pdf").write_bytes(b"pdf")

    monkeypatch.setattr(workbench, "_plot_chromosome", lightweight_chromosome_plot)
    genome_output = tmp_path / "whole_genome"
    result = workbench.plot_workbench_population(
        summaries,
        population="AFR",
        output_dir=genome_output,
        whole_genome=True,
        top_n=2,
    )
    assert result["whole_genome_complete"] is True
    assert result["whole_genome_rows"] == 66
    assert (genome_output / "AFR.whole_genome.gamma_smc.png").is_file()
    assert (genome_output / "AFR.whole_genome.gamma_smc.pdf").is_file()
    assert (genome_output / "AFR.whole_genome.diagnostics.png").is_file()
    top = pd.read_csv(genome_output / "whole_genome_top_windows.tsv", sep="\t")
    assert len(top) == 4
    reused = workbench.plot_workbench_population(
        summaries,
        population="AFR",
        output_dir=genome_output,
        whole_genome=True,
        top_n=2,
    )
    assert reused["reused"] is True


def test_whole_genome_plot_refuses_partial_autosomes(tmp_path):
    summary = tmp_path / "chr1.tsv"
    _write_summary(summary)
    with pytest.raises(ValueError, match="requires all autosomes"):
        workbench.plot_workbench_population(
            {1: summary},
            population="AFR",
            output_dir=tmp_path / "plots",
            whole_genome=True,
        )


def test_workbench_run_report_counts_regions_and_is_idempotent(tmp_path, monkeypatch):
    results_root = tmp_path / "results"
    for population in ("AFR", "EUR"):
        chromosome_root = results_root / population / "chromosomes"
        chromosome_root.mkdir(parents=True)
        for chromosome in workbench.AUTOSOMES:
            _write_summary(
                chromosome_root / f"chr{chromosome}.gamma_smc.tsv", chromosome
            )

    gene_annotation = tmp_path / "gencode.test.gtf"
    gene_annotation.write_text(
        "##description: test GRCh38 annotation\n"
        + "".join(
            f'chr{chromosome}\ttest\tgene\t1\t50000\t.\t+\t.\t'
            f'gene_id "ENSG{chromosome:011d}.1"; gene_type "protein_coding"; '
            f'gene_name "GENE{chromosome}";\n'
            for chromosome in workbench.AUTOSOMES
        ),
        encoding="utf-8",
    )

    combined_calls = []

    def lightweight_combined_plot(*args, output_stem, **kwargs):
        combined_calls.append(output_stem)
        outputs = [Path(f"{output_stem}.png"), Path(f"{output_stem}.pdf")]
        outputs[0].write_bytes(b"png")
        outputs[1].write_bytes(b"pdf")
        return outputs

    monkeypatch.setattr(
        workbench, "_plot_combined_recent_genome", lightweight_combined_plot
    )
    output_dir = tmp_path / "report"
    result = workbench.summarize_workbench_run(
        results_root,
        populations=["AFR", "EUR"],
        chromosomes=list(workbench.AUTOSOMES),
        output_dir=output_dir,
        whole_genome=True,
        top_n=2,
        gene_annotation=gene_annotation,
    )
    assert result["signal_fraction"] == pytest.approx(0.02)
    assert result["population_region_counts"] == {"AFR": 22, "EUR": 22}
    assert result["total_regions"] == 44
    summary = pd.read_csv(output_dir / "regions_by_population.tsv", sep="\t")
    assert summary.set_index("population")["regions_found"].to_dict() == {
        "AFR": 22,
        "EUR": 22,
    }
    assert (output_dir / "all_candidate_regions.tsv").is_file()
    candidate_loci = pd.read_csv(output_dir / "candidate_loci.tsv", sep="\t")
    assert candidate_loci["candidate_merge_max_gap_bp"].eq(20_000).all()
    assert candidate_loci["signal_fraction_strictly_greater_than"].eq(0.02).all()
    raw_scan = pd.read_csv(output_dir / "raw_scan_windows.tsv.gz", sep="\t")
    assert len(raw_scan) == 2 * 22 * 3
    layers = pd.read_csv(output_dir / "report_data_layers.tsv", sep="\t")
    assert layers.set_index("layer")["maximum_gap_bp"].to_dict()[
        "candidate_loci"
    ] == 20_000
    assert layers.set_index("layer")["maximum_gap_bp"].to_dict()[
        "plot_loci"
    ] == 1_000_000
    assert (output_dir / "combined.whole_genome.gamma_smc.png").is_file()
    assert (output_dir / "combined.whole_genome.gamma_smc.zoom4pct.png").is_file()
    gene_list = pd.read_csv(output_dir / "gene_list.tsv", sep="\t")
    plot_loci = pd.read_csv(output_dir / "plot_loci.tsv", sep="\t")
    pd.testing.assert_frame_equal(gene_list, plot_loci)
    assert set(gene_list["probable_gene"]) == {"GENE1", "GENE2"}
    assert len(gene_list) == 4
    assert gene_list["plot_merge_max_gap_bp"].eq(1_000_000).all()
    assert len(combined_calls) == 2

    reused = workbench.summarize_workbench_run(
        results_root,
        populations=["AFR", "EUR"],
        chromosomes=list(workbench.AUTOSOMES),
        output_dir=output_dir,
        whole_genome=True,
        top_n=2,
        gene_annotation=gene_annotation,
    )
    assert reused["reused"] is True
    assert len(combined_calls) == 2

    changed_threshold = workbench.summarize_workbench_run(
        results_root,
        populations=["AFR", "EUR"],
        chromosomes=list(workbench.AUTOSOMES),
        output_dir=output_dir,
        whole_genome=True,
        top_n=2,
        gene_annotation=gene_annotation,
        signal_fraction=0.3,
    )
    assert changed_threshold["reused"] is False
    assert changed_threshold["signal_fraction"] == pytest.approx(0.3)
    recomputed = pd.read_csv(output_dir / "all_candidate_regions.tsv", sep="\t")
    assert len(recomputed) == 44
    assert recomputed["n_signal_windows"].eq(1).all()
    assert len(combined_calls) == 4


def test_ranked_gene_list_merges_top_windows_with_at_most_one_megabase_gap(tmp_path):
    annotation = tmp_path / "genes.gtf"
    annotation.write_text(
        "chr1\ttest\tgene\t1500001\t1600000\t.\t+\t.\t"
        'gene_id "ENSG1.1"; gene_type "protein_coding"; gene_name "GENE_A";\n'
        "chr1\ttest\tgene\t5400001\t5500000\t.\t+\t.\t"
        'gene_id "ENSG2.1"; gene_type "protein_coding"; gene_name "GENE_B";\n',
        encoding="utf-8",
    )
    genes = workbench._load_protein_coding_genes(annotation)
    genome = pd.DataFrame(
        {
            "population": ["AFR"] * 4,
            "chromosome": [1] * 4,
            "position_0based": [1_500_000, 2_100_000, 5_450_000, 8_000_000],
            "position_1based": [1_500_001, 2_100_001, 5_450_001, 8_000_001],
            "frac_recent_4500": [0.04, 0.03, 0.02, 0.001],
            "mean_p_tmrca_lt_threshold": [0.1, 0.08, 0.06, 0.01],
            "mean_tmrca_generations": [100.0, 200.0, 300.0, 20_000.0],
            "genome_position_0based": [
                1_500_000,
                2_100_000,
                5_450_000,
                8_000_000,
            ],
        }
    )
    hits = workbench._build_ranked_gene_list(
        {"AFR": genome},
        called_column="frac_recent_4500",
        genes=genes,
        top_n=3,
        plot_merge_gap=1_000_000,
        context_flank=500_000,
    )
    assert len(hits) == 2
    first = hits.iloc[0]
    assert first["merged_start_0based"] == 1_500_000
    assert first["n_top_windows"] == 2
    assert first["plot_merge_max_gap_bp"] == 1_000_000
    assert first["maximum_joined_gap_bp"] == 600_000
    assert first["probable_gene"] == "GENE_A"
    assert hits.iloc[1]["probable_gene"] == "GENE_B"
