import io
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
        return SimpleNamespace(stdout="", stderr="", returncode=0)

    monkeypatch.setattr("gamma_smc_aou.decoder.subprocess.run", fake_run)
    run_within_decoder(
        "gamma_smc",
        "input.vcf.gz",
        tmp_path / "summary.tsv",
        output_at_hets=False,
        samples=tmp_path / "samples.txt",
    )

    assert "--output_at_hets=false" in captured["command"]
    assert "--output_at_hets" not in captured["command"]
    assert captured["command"][captured["command"].index("--samples") + 1] == str(
        tmp_path / "samples.txt"
    )


def test_output_at_hets_true_is_attached_to_boolean_option(tmp_path, monkeypatch):
    captured = {}

    def fake_run(command, **kwargs):
        captured["command"] = command
        output = Path(command[command.index("--recent_summary") + 1])
        output.write_text(
            "position_0based\tmean_p_tmrca_lt_threshold\n0\t0.0\n",
            encoding="utf-8",
        )
        return SimpleNamespace(stdout="", stderr="", returncode=0)

    monkeypatch.setattr("gamma_smc_aou.decoder.subprocess.run", fake_run)
    run_within_decoder(
        "gamma_smc",
        "input.vcf.gz",
        tmp_path / "summary.tsv",
        output_at_hets=True,
    )

    assert "--output_at_hets=true" in captured["command"]
    assert "--output_at_hets" not in captured["command"]


def test_tsz_input_streams_vcf_to_decoder_stdin(tmp_path, monkeypatch):
    captured = {}

    class FakeProducer:
        def __init__(self, command, **kwargs):
            captured["producer_command"] = command
            self.stdout = io.BytesIO(b"vcf stream")
            self.stderr = io.BytesIO()

        def wait(self):
            return 0

    def fake_run(command, **kwargs):
        captured["command"] = command
        captured["stdin"] = kwargs["stdin"]
        output = Path(command[command.index("--recent_summary") + 1])
        output.write_text(
            "position_0based\tmean_p_tmrca_lt_threshold\n0\t0.0\n",
            encoding="utf-8",
        )
        return SimpleNamespace(stdout="", stderr="", returncode=0)

    monkeypatch.setattr("gamma_smc_aou.decoder.subprocess.Popen", FakeProducer)
    monkeypatch.setattr("gamma_smc_aou.decoder.subprocess.run", fake_run)
    run = run_within_decoder(
        "gamma_smc",
        "input.tsz",
        tmp_path / "summary.tsv",
        input_format="tsz",
    )

    command = captured["command"]
    assert command[command.index("--input") + 1] == "/dev/stdin"
    assert command[command.index("--input_format") + 1] == "vcf"
    assert captured["stdin"] is not None
    assert run["input_path"] == "input.tsz"
    assert run["input_format"] == "tsz"
    assert run["tree_sequence_vcf_producer_command"] == captured["producer_command"]
