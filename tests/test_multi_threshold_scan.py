from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
import pytest

from gamma_smc_aou import bitmatrix, workbench
from gamma_smc_aou.defaults import DEFAULT_THRESHOLD_YEARS, MAX_THRESHOLDS

ROOT = Path(__file__).resolve().parents[1]


def _write_multi_threshold_summary(path: Path) -> None:
    """A scan summary decoded at 10,000 and 50,000 years.

    ``mean_p_tmrca_lt_threshold`` aliases the *first* threshold, exactly as the
    decoder writes it.
    """
    pd.DataFrame(
        {
            "position_0based": [0, 10_000, 20_000],
            "position_1based": [1, 10_001, 20_001],
            "n_pairs": [100, 100, 98],
            "mean_p_tmrca_lt_threshold": [0.011, 0.012, 0.013],
            "mean_tmrca_generations": [30_000.0, 29_000.0, 28_000.0],
            "n_recent_10000": [1, 2, 3],
            "frac_recent_10000": [0.01, 0.02, 0.03],
            "mean_p_lt_10000": [0.011, 0.012, 0.013],
            "n_recent_50000": [6, 7, 8],
            "frac_recent_50000": [0.06, 0.07, 0.08],
            "mean_p_lt_50000": [0.061, 0.071, 0.081],
        }
    ).to_csv(path, sep="\t", index=False)


# --------------------------------------------------------------------------
# threshold_suffix / format_threshold agreement
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "years,expected",
    [
        (4500, "4500"),
        (10_000, "10000"),
        (50_000, "50000"),
        (100_000.0, "100000"),
        (4500.5, "4500.5"),
        (2.25, "2.25"),
        # Six *significant* digits would have rendered these in scientific
        # notation and collapsed them onto one another.
        (1_234_567.5, "1234567.5"),
        (1_234_567.1, "1234567.1"),
        (1_234_567.2, "1234567.2"),
    ],
)
def test_threshold_suffix_is_decimal_and_lossless(years, expected):
    assert workbench.threshold_suffix(years) == expected
    assert bitmatrix._format_threshold(years) == expected


def test_every_python_formatter_shares_one_rule():
    """There must be exactly one implementation of the column-suffix rule.

    container_study carried a fourth copy using str(float), which agreed with
    the others on integral thresholds and silently diverged otherwise.
    """
    from gamma_smc_aou import container_study

    assert container_study._threshold_suffix is workbench.threshold_suffix
    for value in (4500, 50_000, 4500.5, 1_234_567.5):
        assert (
            container_study._threshold_suffix(value)
            == workbench.threshold_suffix(value)
            == bitmatrix._format_threshold(value)
        )


def test_threshold_suffix_matches_the_cpp_rule():
    """The C++ writes the column names the Python readers look up.

    Guards the divergence this replaced: the C++ used to fall through to the
    stream's default %g, which is six significant digits, while Python used
    the full decimal value.
    """
    source = (ROOT / "src" / "gamma_smc.h").read_text()
    body = source[source.index("static string format_threshold(double years)") :]
    body = body[: body.index("\n    }")]
    assert "std::fixed" in body
    assert "std::setprecision(6)" in body
    assert "find_last_not_of('0')" in body
    # The bare `stream << years` fallback is what produced scientific notation.
    assert not re.search(r"stream << years\s*;", body)


# --------------------------------------------------------------------------
# parse_threshold_years
# --------------------------------------------------------------------------


def test_parse_threshold_years_preserves_decode_order():
    assert workbench.parse_threshold_years([50_000, 10_000]) == [50_000.0, 10_000.0]
    assert workbench.parse_threshold_years(4500) == [4500.0]


@pytest.mark.parametrize(
    "values",
    [
        [],
        [0],
        [-1],
        [float("inf")],
        [float("nan")],
        ["not-a-number"],
        list(range(1, MAX_THRESHOLDS + 2)),
        [100, 100.0],
        # Distinct floats that collide once formatted to six decimals.
        [1.0000001, 1.0000002],
    ],
)
def test_parse_threshold_years_rejects_bad_input(values):
    with pytest.raises(ValueError):
        workbench.parse_threshold_years(values)


def test_default_thresholds_are_within_the_cap():
    assert 1 <= len(DEFAULT_THRESHOLD_YEARS) <= MAX_THRESHOLDS
    assert list(DEFAULT_THRESHOLD_YEARS) == [10_000.0, 50_000.0]


# --------------------------------------------------------------------------
# split_summary_by_threshold
# --------------------------------------------------------------------------


@pytest.mark.parametrize("threshold", [10_000, 50_000])
def test_split_rewrites_the_soft_alias_for_its_own_threshold(tmp_path, threshold):
    combined = tmp_path / "chr1.gamma_smc.tsv"
    _write_multi_threshold_summary(combined)
    suffix = workbench.threshold_suffix(threshold)

    destination = workbench.split_summary_by_threshold(
        combined, threshold, tmp_path / f"t{suffix}" / "chr1.gamma_smc.tsv"
    )
    frame = pd.read_csv(destination, sep="\t")

    # The alias must follow the requested threshold, not the decoded first one.
    assert frame[workbench.SOFT_COLUMN].tolist() == frame[f"mean_p_lt_{suffix}"].tolist()
    assert workbench.called_fraction_column(threshold) in frame.columns
    # The other threshold's block is gone.
    other = 50_000 if threshold == 10_000 else 10_000
    assert workbench.called_fraction_column(other) not in frame.columns
    # Historical five-column prefix is preserved for existing consumers.
    assert list(frame.columns)[:5] == [
        "position_0based",
        "position_1based",
        "n_pairs",
        workbench.SOFT_COLUMN,
        workbench.TMRCA_COLUMN,
    ]
    # And the result is a valid single-threshold summary.
    workbench.validate_summary(destination, threshold_years=threshold)


def test_split_refuses_a_threshold_that_was_not_decoded(tmp_path):
    combined = tmp_path / "chr1.gamma_smc.tsv"
    _write_multi_threshold_summary(combined)
    with pytest.raises(ValueError, match="no block for threshold 4500"):
        workbench.split_summary_by_threshold(combined, 4500, tmp_path / "out.tsv")


# --------------------------------------------------------------------------
# completion contract
# --------------------------------------------------------------------------


def _contract(**overrides):
    base = dict(
        population="AFR",
        chromosome=1,
        input_uri="gs://bucket/chr1.bcf",
        input_fingerprint="fingerprint",
        local_input="chr1.bcf",
        index_uri="gs://bucket/chr1.bcf.csi",
        index_fingerprint="fingerprint",
        local_index="chr1.bcf.csi",
        output_summary="summary.tsv",
        pairs_manifest="pairs.tsv",
        sample_list="samples.txt",
        sample_audit="samples.json",
        ancestry_uri="gs://bucket/ancestry.tsv",
        ancestry_fingerprint="fingerprint",
        local_ancestry="ancestry.tsv",
        qc_exclusions_uri="gs://bucket/qc.tsv",
        qc_exclusions_fingerprint="fingerprint",
        local_qc_exclusions="qc.tsv",
        relatedness_exclusions_uri="gs://bucket/rel.tsv",
        relatedness_exclusions_fingerprint="fingerprint",
        local_relatedness_exclusions="rel.tsv",
        mask_uri=None,
        mask_fingerprint=None,
        local_mask_source=None,
        local_mask=None,
        mask_audit=None,
        theta=0.00075,
        rho_over_theta=0.8,
        mutation_rate=1.25e-8,
        generation_time=25,
        recent_call="mean",
        stride=10_000,
        cache_size=1_000,
        threads=12,
        pair_block=256,
        exp10="accurate",
        backward_alignment="fixed",
        code_commit="deadbeef",
    )
    base.update(overrides)
    return workbench.build_workbench_contract(**base)


def test_contract_records_every_decoded_threshold():
    contract = _contract(threshold_years=10_000, threshold_years_all=[10_000, 50_000])
    decoder = contract["decoder"]
    assert decoder["threshold_years"] == 10_000.0
    assert decoder["threshold_years_all"] == [10_000.0, 50_000.0]
    assert workbench._contract_thresholds(decoder) == [10_000.0, 50_000.0]


def test_contract_defaults_the_list_to_the_single_threshold():
    decoder = _contract(threshold_years=4500)["decoder"]
    assert decoder["threshold_years_all"] == [4500.0]


def test_legacy_contract_without_the_list_still_reads():
    # Completion records written before multi-threshold support.
    assert workbench._contract_thresholds({"threshold_years": 4500.0}) == [4500.0]


def test_contract_rejects_a_candidate_threshold_outside_the_decoded_set():
    with pytest.raises(ValueError, match="not among the decoded thresholds"):
        _contract(threshold_years=4500, threshold_years_all=[10_000, 50_000])


def test_contract_rejects_more_than_the_cap():
    with pytest.raises(ValueError, match=f"at most {MAX_THRESHOLDS}"):
        _contract(
            threshold_years=1,
            threshold_years_all=list(range(1, MAX_THRESHOLDS + 2)),
        )
