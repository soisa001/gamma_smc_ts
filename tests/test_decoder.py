from pathlib import Path
from types import SimpleNamespace

from gamma_smc_aou.decoder import run_within_decoder


def test_output_at_hets_false_is_attached_to_boolean_option(tmp_path, monkeypatch):
    captured = {}

    def fake_run(command, **kwargs):
        captured["command"] = command
        output = Path(command[command.index("--recent_summary") + 1])
        output.write_text(
            "position_0based\tmean_p_tmrca_lt_threshold\n0\t0.0\n",
            encoding="utf-8",
        )
        return SimpleNamespace(stdout="", stderr="")

    monkeypatch.setattr("gamma_smc_aou.decoder.subprocess.run", fake_run)
    run_within_decoder(
        "gamma_smc",
        "input.vcf.gz",
        tmp_path / "summary.tsv",
        output_at_hets=False,
    )

    assert "--output_at_hets=false" in captured["command"]
    assert "--output_at_hets" not in captured["command"]


def test_output_at_hets_true_is_attached_to_boolean_option(tmp_path, monkeypatch):
    captured = {}

    def fake_run(command, **kwargs):
        captured["command"] = command
        output = Path(command[command.index("--recent_summary") + 1])
        output.write_text(
            "position_0based\tmean_p_tmrca_lt_threshold\n0\t0.0\n",
            encoding="utf-8",
        )
        return SimpleNamespace(stdout="", stderr="")

    monkeypatch.setattr("gamma_smc_aou.decoder.subprocess.run", fake_run)
    run_within_decoder(
        "gamma_smc",
        "input.vcf.gz",
        tmp_path / "summary.tsv",
        output_at_hets=True,
    )

    assert "--output_at_hets=true" in captured["command"]
    assert "--output_at_hets" not in captured["command"]
