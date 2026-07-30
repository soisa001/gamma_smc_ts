from importlib import resources
from pathlib import Path

from gamma_smc_aou.cli import parser


def test_uv_portability_files_are_present():
    root = Path(__file__).resolve().parents[1]
    for relative in (
        "uv.lock",
        ".python-version",
        "scripts/bootstrap_uv.sh",
        "scripts/bootstrap_uv.ps1",
        "scripts/aou.sh",
        "scripts/aou.ps1",
        "scripts/verify_install.py",
    ):
        assert (root / relative).is_file(), relative
    assert resources.files("gamma_smc_aou").joinpath(
        "slim/recent_sweep.slim"
    ).is_file()


def test_decode_defaults_to_wrapper_decoder_environment(monkeypatch):
    monkeypatch.setenv("GAMMA_SMC_BIN", "/portable/bin/gamma_smc")
    arguments = parser().parse_args(
        [
            "decode",
            "--input",
            "input.trees",
            "--output",
            "output.tsv",
            "--theta",
            "0.0005",
            "--rho-over-theta",
            "0.8",
            "--mutation-rate",
            "1.25e-8",
        ]
    )
    assert arguments.executable == "/portable/bin/gamma_smc"


def test_workbench_defaults_are_shared_by_the_cli():
    arguments = parser().parse_args(
        ["decode", "--input", "input.vcf.gz", "--output", "output.tsv"]
    )
    assert arguments.mutation_rate == 1.29e-9
    assert arguments.output_at_stride == 10_000
    assert arguments.cache_size == 1_000

    study = parser().parse_args(
        [
            "run-native-study",
            "--source-dir", "source",
            "--output-dir", "output",
        ]
    )
    assert study.output_at_stride == 10_000
    assert study.cache_size == 1_000
    assert study.workers == 12
    assert study.recent_call == "median"
    assert study.compare_recent_call is None
    assert study.calibration_statistic == "mean-posterior-probability"

    comparison = parser().parse_args(
        [
            "run-native-study",
            "--source-dir", "source",
            "--output-dir", "output",
            "--recent-call", "mean",
            "--compare-recent-call", "median",
            "--calibration-statistic", "called-fraction",
        ]
    )
    assert comparison.recent_call == "mean"
    assert comparison.compare_recent_call == "median"
    assert comparison.calibration_statistic == "called-fraction"

    high_af = parser().parse_args(
        [
            "prepare-high-af-selected",
            "--null-truth-dir", "null",
            "--output-dir", "output",
            "--selected-attempt", "2514",
        ]
    )
    assert high_af.selected_attempt == 2514
