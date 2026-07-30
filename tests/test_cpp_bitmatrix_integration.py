"""End-to-end tests of the compiled decoder's new outputs.

These need GAMMA_SMC_BIN and are skipped without it. What they check that the
pure-Python tests cannot: that the bit matrix the C++ writes agrees with the
summary TSV the same run produced, that threading does not change the answer,
and that random pair sampling is reproducible.
"""

import json
import os
import re
import struct
import subprocess
from pathlib import Path

import msprime
import numpy as np
import pandas as pd
import pytest

from gamma_smc_aou import bitmatrix


pytestmark = pytest.mark.skipif(
    "GAMMA_SMC_BIN" not in os.environ, reason="compiled Linux binary not available"
)

THRESHOLDS = (4500.0, 10000.0)
STRIDE = 1000


def panel_ts(n_individuals=24, length=60_000):
    """A panel big enough that C(2N, 2) comfortably exceeds the sampled pairs.

    24 diploids give 48 haplotypes and 1,128 distinct pairs, so a 64-pair sample
    spans several 8-pair chunks and more than one default block.
    """
    ancestry = msprime.sim_ancestry(
        samples=[msprime.SampleSet(n_individuals, ploidy=2)],
        population_size=1_000,
        sequence_length=length,
        recombination_rate=1e-8,
        model=msprime.StandardCoalescent(),
        random_seed=42,
    )
    return msprime.sim_mutations(
        ancestry, rate=2e-7, model=msprime.BinaryMutationModel(), random_seed=43
    )


def decode(tmp_path, ts, *, name, extra=(), n_random_pairs=0, threads=1, bits=True,
           only_within=None):
    source = tmp_path / f"{name}.trees"
    ts.dump(source)
    summary = tmp_path / f"{name}.tsv"
    command = [
        os.environ["GAMMA_SMC_BIN"],
        "--input", str(source), "--input_format", "trees",
        "--scaled_mutation_rate", "0.0008",
        "--recombination_to_mutation_ratio", "0.05",
        "--unscaled_mutation_rate", "2e-7",
        "--recent_threshold_years", ",".join(str(value) for value in THRESHOLDS),
        "--generation_time", "25",
        "--recent_summary", str(summary),
        "--output_at_hets=false",
        "--output_at_stride", str(STRIDE),
        "--cache_size", "10",
        "--threads", str(threads),
    ]
    if only_within is None:
        only_within = (n_random_pairs == 0) and not any("--pairs_file" in str(a) for a in extra)
    if n_random_pairs:
        command += ["--n_random_pairs", str(n_random_pairs), "--pairs_seed", "20240727"]
    elif only_within:
        command += ["--only_within"]
    if bits:
        command += ["--recent_bitmatrix", str(tmp_path / f"{name}.bits")]
    command += list(extra)

    completed = subprocess.run(command, check=True, text=True, capture_output=True)
    return summary, tmp_path / f"{name}.bits", completed


def test_bitmatrix_counts_match_the_summary_tsv(tmp_path):
    ts = panel_ts()
    summary, bits_path, _ = decode(tmp_path, ts, name="within", n_random_pairs=0)

    frame = pd.read_csv(summary, sep="\t")
    counts = bitmatrix.position_counts(bits_path)
    meta = bitmatrix.read_meta(bits_path)

    assert list(meta["thresholds_years"]) == list(THRESHOLDS)
    assert len(frame) == meta["n_positions"]
    np.testing.assert_array_equal(
        bitmatrix.output_positions(meta), frame["position_0based"].to_numpy()
    )
    np.testing.assert_array_equal(counts[0], frame["n_recent_4500"].to_numpy())
    np.testing.assert_array_equal(counts[1], frame["n_recent_10000"].to_numpy())


def test_summary_keeps_the_legacy_schema(tmp_path):
    ts = panel_ts()
    summary, _, _ = decode(tmp_path, ts, name="legacy")
    frame = pd.read_csv(summary, sep="\t")

    for column in (
        "position_0based", "position_1based", "n_pairs",
        "mean_p_tmrca_lt_threshold", "mean_tmrca_generations",
    ):
        assert column in frame.columns

    # The legacy alias is the first threshold's mean probability.
    np.testing.assert_allclose(
        frame["mean_p_tmrca_lt_threshold"].to_numpy(),
        frame["mean_p_lt_4500"].to_numpy(),
        rtol=0, atol=0,
    )
    assert frame["mean_p_tmrca_lt_threshold"].between(0, 1).all()
    assert frame["frac_recent_4500"].between(0, 1).all()
    # A longer horizon can only catch more pairs.
    assert (frame["n_recent_10000"] >= frame["n_recent_4500"]).all()
    assert (frame["mean_p_lt_10000"] >= frame["mean_p_lt_4500"] - 1e-6).all()


def test_thread_count_does_not_change_the_result(tmp_path):
    ts = panel_ts()
    one, bits_one, _ = decode(tmp_path, ts, name="t1", n_random_pairs=64, threads=1)
    many, bits_many, _ = decode(tmp_path, ts, name="t8", n_random_pairs=64, threads=8)

    frame_one = pd.read_csv(one, sep="\t")
    frame_many = pd.read_csv(many, sep="\t")
    pd.testing.assert_frame_equal(frame_one, frame_many)

    np.testing.assert_array_equal(
        bitmatrix.position_counts(bits_one), bitmatrix.position_counts(bits_many)
    )


def test_pair_block_does_not_change_the_result(tmp_path):
    ts = panel_ts()
    small, bits_small, _ = decode(
        tmp_path, ts, name="b8", n_random_pairs=50, threads=4, extra=("--pair_block", "8")
    )
    large, bits_large, _ = decode(
        tmp_path, ts, name="b512", n_random_pairs=50, threads=4, extra=("--pair_block", "512")
    )
    pd.testing.assert_frame_equal(pd.read_csv(small, sep="\t"), pd.read_csv(large, sep="\t"))
    np.testing.assert_array_equal(
        bitmatrix.position_counts(bits_small), bitmatrix.position_counts(bits_large)
    )


def test_random_pairs_are_reproducible_and_well_formed(tmp_path):
    ts = panel_ts()
    _, bits_a, _ = decode(tmp_path, ts, name="ra", n_random_pairs=20)
    _, bits_b, _ = decode(tmp_path, ts, name="rb", n_random_pairs=20)

    meta_a = bitmatrix.read_meta(bits_a)
    meta_b = bitmatrix.read_meta(bits_b)
    assert meta_a["pairs"] == meta_b["pairs"]
    assert meta_a["n_pairs"] == 20

    pairs = bitmatrix.pairs(meta_a)
    n_haplotypes = 2 * len(ts.individuals())
    assert pairs.min() >= 0 and pairs.max() < n_haplotypes
    assert (pairs[:, 0] < pairs[:, 1]).all()
    assert len({tuple(pair) for pair in pairs}) == 20


def test_exclude_within_drops_same_individual_pairs(tmp_path):
    ts = panel_ts()
    _, bits_path, _ = decode(
        tmp_path, ts, name="xw", n_random_pairs=20, extra=("--exclude_within",)
    )
    pairs = bitmatrix.pairs(bitmatrix.read_meta(bits_path))
    assert ((pairs[:, 0] >> 1) != (pairs[:, 1] >> 1)).all()


def test_pairs_file_is_honoured(tmp_path):
    ts = panel_ts()
    wanted = [(0, 3), (1, 6), (2, 7)]
    pairs_file = tmp_path / "pairs.txt"
    pairs_file.write_text(
        "# haplotype pairs\n" + "\n".join(f"{i}\t{j}" for i, j in wanted) + "\n"
    )
    _, bits_path, _ = decode(
        tmp_path, ts, name="pf", n_random_pairs=0, extra=("--pairs_file", str(pairs_file))
    )
    meta = bitmatrix.read_meta(bits_path)
    assert [tuple(pair) for pair in meta["pairs"]] == wanted


def test_asking_for_more_pairs_than_exist_fails_cleanly(tmp_path):
    ts = panel_ts()
    with pytest.raises(subprocess.CalledProcessError) as excinfo:
        decode(tmp_path, ts, name="toomany", n_random_pairs=10_000)
    assert "exceeds" in excinfo.value.stdout


def test_lookup_tables_self_check_within_tolerance(tmp_path):
    ts = panel_ts()
    _, _, completed = decode(tmp_path, ts, name="selfcheck")
    log = completed.stdout

    error = float(log.split("max abs error vs Boost")[1].split()[0])
    disagreement = float(log.split("disagreement with Boost")[1].split()[0])
    print(f"\nP(T<t) table max abs error vs Boost: {error:.3e}")
    print(f"hard-call disagreement vs Boost:     {disagreement:.3e}")
    assert error < 1e-3, f"P(T<t) table error {error:g} is worse than expected"
    assert disagreement < 1e-4, f"call table disagreement {disagreement:g} is too high"


def test_mean_call_rule_is_a_threshold_on_the_posterior_mean(tmp_path):
    ts = panel_ts()
    summary, _, _ = decode(
        tmp_path, ts, name="meanrule", extra=("--recent_call", "mean"), bits=False
    )
    frame = pd.read_csv(summary, sep="\t")
    # 10000 years at 25 years/generation is 400 generations.
    assert frame["frac_recent_10000"].between(0, 1).all()
    assert (frame["n_recent_10000"] >= frame["n_recent_4500"]).all()


def test_raw_posteriors_still_round_trip(tmp_path):
    ts = panel_ts()
    raw = tmp_path / "posteriors.zst"
    summary, _, _ = decode(
        tmp_path, ts, name="raw", n_random_pairs=20, threads=4,
        extra=("--output", str(raw)), bits=False,
    )
    meta = json.loads((tmp_path / "posteriors.zst.meta").read_text())
    assert meta["num_pairs"] == 20
    assert meta["chunk_size"] == 8
    assert meta["sequence_length"] == len(pd.read_csv(summary, sep="\t"))

    import zstandard

    raw_floats = np.frombuffer(
        zstandard.ZstdDecompressor().stream_reader(raw.open("rb")).read(), dtype=np.float32
    )
    n_chunks = -(-meta["num_pairs"] // meta["chunk_size"])
    assert raw_floats.size == n_chunks * 2 * meta["sequence_length"] * meta["chunk_size"]


def test_pair_subset_counts_agree_with_the_whole(tmp_path):
    ts = panel_ts()
    _, bits_path, _ = decode(tmp_path, ts, name="subset", n_random_pairs=40, threads=4)
    meta = bitmatrix.read_meta(bits_path)
    total = bitmatrix.position_counts(bits_path, meta)
    first = bitmatrix.position_counts(bits_path, meta, pair_indices=range(0, 20))
    second = bitmatrix.position_counts(bits_path, meta, pair_indices=range(20, 40))
    np.testing.assert_array_equal(total, first + second)


def manifest_path(completed):
    """The path the binary reports, rather than re-deriving its anchor rule."""
    match = re.search(r"Pair manifest: (\S+)", strip_ansi(completed.stdout))
    assert match, strip_ansi(completed.stdout)[-2000:]
    return Path(match.group(1))


def strip_ansi(text):
    return re.sub(r"\[[0-9;]*m", "", text)


def read_manifest(path):
    """Split a pair manifest into its '#' header fields and its pair rows."""
    header, pairs = {}, []
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            body = line[1:].strip()
            parts = body.split("\t")
            if len(parts) == 2:
                header[parts[0]] = parts[1]
        elif line.strip():
            fields = line.split("\t")
            pairs.append((int(fields[0]), int(fields[1])))
    return header, pairs


def test_random_draw_writes_a_manifest_without_being_asked(tmp_path):
    # A seed alone does not identify a draw once the panel can change, so the
    # manifest is derived from the output path rather than left to a flag.
    summary, _, completed = decode(tmp_path, panel_ts(), name="auto", n_random_pairs=64)
    manifest = manifest_path(completed)
    assert manifest.exists()
    assert manifest == summary.with_name(summary.name + ".pairs.tsv")

    header, pairs = read_manifest(manifest)
    assert header["mode"] == "random"
    assert header["n_pairs"] == "64"
    assert header["n_haplotypes"] == "48"
    assert header["pairs_seed"] == "20240727"
    assert header["panel_digest"].startswith("0x")
    assert len(pairs) == 64
    assert len(set(pairs)) == 64
    assert all(0 <= i < j < 48 for i, j in pairs)


def test_manifest_replays_as_a_pairs_file(tmp_path):
    ts = panel_ts()
    drawn, _, completed = decode(tmp_path, ts, name="drawn", n_random_pairs=64)
    manifest = manifest_path(completed)

    replayed, _, _ = decode(
        tmp_path, ts, name="replayed", extra=["--pairs_file", str(manifest)]
    )
    pd.testing.assert_frame_equal(
        pd.read_csv(drawn, sep="\t"), pd.read_csv(replayed, sep="\t")
    )


def test_truncated_manifest_is_refused(tmp_path):
    ts = panel_ts()
    _, _, completed = decode(tmp_path, ts, name="complete", n_random_pairs=16)
    manifest = manifest_path(completed)
    lines = manifest.read_text().splitlines()
    pair_rows = [index for index, line in enumerate(lines) if line and not line.startswith("#")]
    assert pair_rows
    lines.pop(pair_rows[-1])
    truncated = tmp_path / "truncated.pairs.tsv"
    truncated.write_text("\n".join(lines) + "\n")

    with pytest.raises(subprocess.CalledProcessError) as caught:
        decode(
            tmp_path,
            ts,
            name="truncated",
            extra=["--pairs_file", str(truncated)],
        )
    assert "does not match its manifest header" in strip_ansi(caught.value.stdout)


def test_manifest_records_haplotype_names_alongside_indices(tmp_path):
    _, _, completed = decode(tmp_path, panel_ts(), name="named", n_random_pairs=16)
    manifest = manifest_path(completed)
    rows = [
        line.split("\t")
        for line in manifest.read_text().splitlines()
        if line and not line.startswith("#")
    ]
    for hap_i, hap_j, label_i, label_j in rows:
        # Haplotype h is haplotype (h & 1) of diploid (h >> 1).
        assert label_i.endswith(f".{int(hap_i) & 1}")
        assert label_j.endswith(f".{int(hap_j) & 1}")


def reduced_manifest(tmp_path, completed, n_haplotypes, destination):
    """A manifest keeping only pairs that are valid indices in a smaller panel.

    This is what makes the digest check meaningful: every index is in range, so
    range validation cannot object, and only the recorded panel identity
    distinguishes the two panels.
    """
    source = manifest_path(completed)
    header = [line for line in source.read_text().splitlines() if line.startswith("#")]
    _, pairs = read_manifest(source)
    usable = [pair for pair in pairs if max(pair) < n_haplotypes]
    assert usable, "no drawn pair fits the smaller panel"
    digest = 1469598103934665603
    for pair in usable:
        for value in struct.pack("=ii", *pair):
            digest ^= value
            digest = (digest * 1099511628211) & ((1 << 64) - 1)
    header = [
        f"# n_pairs\t{len(usable)}" if line.startswith("# n_pairs\t")
        else f"# pairs_digest\t0x{digest:016x}" if line.startswith("# pairs_digest\t")
        else line
        for line in header
    ]
    rows = [f"{i}\t{j}" for i, j in usable]
    destination.write_text("\n".join(header + rows) + "\n")
    return destination, usable


def test_reusing_a_manifest_against_another_panel_is_refused(tmp_path):
    _, _, completed = decode(
        tmp_path, panel_ts(n_individuals=24), name="big", n_random_pairs=64
    )
    reduced, _ = reduced_manifest(tmp_path, completed, 40, tmp_path / "big.reduced.tsv")

    with pytest.raises(subprocess.CalledProcessError) as caught:
        decode(
            tmp_path, panel_ts(n_individuals=20), name="small",
            extra=["--pairs_file", str(reduced)],
        )
    stdout = strip_ansi(caught.value.stdout)
    assert "different panel" in stdout
    # Every index is valid for the smaller panel, so the range check in
    # read_pairs_file would have passed this file without complaint.
    assert "outside" not in stdout


def test_panel_mismatch_can_be_overridden(tmp_path):
    _, _, completed = decode(
        tmp_path, panel_ts(n_individuals=24), name="ovbig", n_random_pairs=32
    )
    reduced, usable = reduced_manifest(tmp_path, completed, 40, tmp_path / "ov.reduced.tsv")

    summary, _, completed = decode(
        tmp_path, panel_ts(n_individuals=20), name="ovsmall",
        extra=["--pairs_file", str(reduced), "--allow_panel_mismatch"],
    )
    assert "Warning" in strip_ansi(completed.stdout)
    assert pd.read_csv(summary, sep="\t")["n_pairs"].max() == len(usable)


def test_explicit_manifest_path_is_honoured_for_any_mode(tmp_path):
    target = tmp_path / "nested" / "within.pairs.tsv"
    decode(
        tmp_path, panel_ts(), name="within", only_within=True,
        extra=["--pairs_manifest", str(target)],
    )
    header, pairs = read_manifest(target)
    assert header["mode"] == "only_within"
    assert "pairs_seed" not in header      # meaningless for a deterministic list
    assert pairs == [(2 * i, 2 * i + 1) for i in range(24)]


def test_runs_with_no_rate_flags_at_all(tmp_path):
    # The three rates now have reference defaults, so a decode needs only an
    # input, an output and a pair selection.
    ts = panel_ts()
    source = tmp_path / "defaults.trees"
    ts.dump(source)
    summary = tmp_path / "defaults.tsv"
    completed = subprocess.run(
        [
            os.environ["GAMMA_SMC_BIN"], "--input", str(source),
            "--input_format", "trees", "--only_within",
            "--recent_summary", str(summary),
            "--output_at_hets=false", "--output_at_stride", str(STRIDE),
            "--cache_size", "10",
            "--threads", "1",
        ],
        check=True, text=True, capture_output=True,
    )
    stdout = strip_ansi(completed.stdout)
    assert "Scaled mutation rate: 0.000750" in stdout
    assert "Scaled recombination rate: 0.000600" in stdout
    frame = pd.read_csv(summary, sep="\t")
    assert frame["mean_p_tmrca_lt_threshold"].between(0, 1).all()


def test_stride_output_covers_declared_invariant_tail(tmp_path):
    ts = panel_ts(length=60_000)
    summary, _, _ = decode(tmp_path, ts, name="tail", bits=False)
    positions = pd.read_csv(summary, sep="\t")["position_0based"].tolist()
    assert positions == list(range(0, 60_000, STRIDE))


def test_multi_contig_vcf_is_refused(tmp_path):
    vcf = tmp_path / "two_contigs.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=2000>\n"
        "##contig=<ID=chr2,length=2000>\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\n"
        "chr1\t101\t.\tA\tG\t.\tPASS\t.\tGT\t0|0\t1|1\n"
        "chr2\t101\t.\tA\tG\t.\tPASS\t.\tGT\t0|0\t1|1\n"
    )
    completed = subprocess.run(
        [
            os.environ["GAMMA_SMC_BIN"],
            "--input", str(vcf),
            "--input_format", "vcf",
            "--only_within",
            "--recent_summary", str(tmp_path / "multi.tsv"),
            "--output_at_stride", "1000",
            "--cache_size", "10",
        ],
        text=True,
        capture_output=True,
    )
    assert completed.returncode != 0
    assert "more than one contig" in strip_ansi(completed.stdout)


def test_mask_must_cover_the_input_contig(tmp_path):
    vcf = tmp_path / "chr1.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=2000>\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\n"
        "chr1\t101\t.\tA\tG\t.\tPASS\t.\tGT\t0|0\t1|1\n"
    )
    mask = tmp_path / "wrong.bed"
    mask.write_text("chr2\t0\t2000\n")
    completed = subprocess.run(
        [
            os.environ["GAMMA_SMC_BIN"],
            "--input", str(vcf),
            "--input_format", "vcf",
            "--only_within",
            "--recent_summary", str(tmp_path / "wrong.tsv"),
            "--mask", str(mask),
            "--output_at_stride", "1000",
            "--cache_size", "10",
        ],
        text=True,
        capture_output=True,
    )
    assert completed.returncode != 0
    assert "has no records for input contig chr1" in strip_ansi(completed.stdout)


def test_estimating_theta_is_opt_in(tmp_path):
    # Omitting --scaled_mutation_rate used to estimate it from heterozygosity,
    # which silently rescaled the time axis per file. That is now explicit.
    ts = panel_ts()
    source = tmp_path / "est.trees"
    ts.dump(source)
    common = [
        os.environ["GAMMA_SMC_BIN"], "--input", str(source),
        "--input_format", "trees", "--only_within",
        "--output_at_hets=false", "--output_at_stride", str(STRIDE),
        "--cache_size", "10",
        "--threads", "1", "--recombination_to_mutation_ratio", "0.8",
    ]
    estimated = subprocess.run(
        common + ["--recent_summary", str(tmp_path / "est.tsv"),
                  "--estimate_mutation_rate"],
        check=True, text=True, capture_output=True,
    )
    assert "Estimating mutation rate from heterozygosity" in strip_ansi(estimated.stdout)
    assert "Scaled mutation rate: 0.000750" not in strip_ansi(estimated.stdout)

    # And it conflicts with pinning theta explicitly.
    with pytest.raises(subprocess.CalledProcessError) as caught:
        subprocess.run(
            common + ["--recent_summary", str(tmp_path / "conflict.tsv"),
                      "--estimate_mutation_rate", "--scaled_mutation_rate", "0.0005"],
            check=True, text=True, capture_output=True,
        )
    assert "mutually exclusive" in strip_ansi(caught.value.stdout)
