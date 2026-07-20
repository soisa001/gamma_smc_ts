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
