"""Local staging-control tests with mocked cloud/bcftools; no cloud data access."""
import argparse
import importlib.util
import json
from pathlib import Path
import re
import subprocess
import sys

import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
spec = importlib.util.spec_from_file_location("run_singleton_chr", SCRIPTS / "run_singleton_chr.py")
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


def args(tmp_path):
    return argparse.Namespace(root=tmp_path, chrom="chr1", gcloud="gcloud", bcftools="bcftools",
                              billing_project="test-billing", io_threads=4, scratch_reserve_gb=0)


def setup_stage(tmp_path, monkeypatch, have_csi=True):
    cfg = args(tmp_path)
    uri = runner.source_paths("chr1")["hiphase"]
    source = tmp_path / "mock_source"
    source.write_bytes(b"synthetic bytes, not a real BCF")
    checksum = runner.crc32c(source)
    metadata = {"bcf": {"uri": uri, "size": source.stat().st_size, "generation": "1", "crc32c": checksum},
                "csi": {"uri": uri + ".csi", "size": source.stat().st_size, "generation": "2", "crc32c": checksum} if have_csi else None}
    invocations = []
    def fake_command(cmd, capture=False, env=None):
        cmd = list(map(str, cmd))
        invocations.append((cmd, env))
        if cmd[1:3] == ["storage", "rsync"]:
            dest = Path(cmd[4])
            (dest / uri.rsplit("/", 1)[1]).write_bytes(source.read_bytes())
            if have_csi:
                (dest / (uri.rsplit("/", 1)[1] + ".csi")).write_bytes(source.read_bytes())
        elif cmd[1:3] == ["index", "--force"]:
            Path(cmd[-1] + ".csi").write_bytes(b"synthetic index")
        elif cmd[1:3] == ["index", "--nrecords"]:
            return "100\n"
        else:
            raise AssertionError(cmd)
    monkeypatch.setattr(runner, "command", fake_command)
    monkeypatch.setattr(runner, "metadata_for", lambda *a: metadata)
    return cfg, metadata, invocations


def test_exact_chr1_filter_and_local_destination(tmp_path):
    for name, uri in runner.source_paths("chr1").items():
        cmd = runner.sync_command("gcloud", uri, tmp_path / name, "billing")
        regex = cmd[cmd.index("--exclude") + 1]
        basename = uri.rsplit("/", 1)[1]
        assert not re.search(regex, basename)
        assert not re.search(regex, basename + ".csi")
        assert re.search(regex, basename.replace("chr1.", "chr10."))
        assert re.search(regex, basename.replace("chr1.", "chr2."))
        assert re.search(regex, "folder/" + basename)
        assert not cmd[4].startswith("gs://")
        assert "--delete-unmatched-destination-objects" not in cmd
        assert cmd[-2:] == ["--billing-project", "billing"]
    assert "/panel/panel/panel_bubble_split_vcf/" in runner.source_paths("chr1")["scaffold"]


def test_dry_run_calls_no_tools_and_creates_no_files(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(runner, "command", lambda *a, **k: pytest.fail("dry run invoked tool"))
    root = tmp_path / "not_created"
    runner.main(["--root", str(root), "--dry-run"])
    plan = json.loads(capsys.readouterr().out)
    assert plan["chromosome"] == "chr1"
    assert plan["active_chromosome_workers"] == 1
    assert plan["workers"] == 4
    assert plan["uploads"] is False
    assert not root.exists()


@pytest.mark.parametrize("error_text,optional_missing", [("HTTPError 404 object not found", True), ("403 Permission denied", False)])
def test_optional_index_missing_does_not_mask_access_errors(monkeypatch, error_text, optional_missing):
    def fail(cmd, **kwargs):
        raise subprocess.CalledProcessError(1, cmd, stderr=error_text)
    monkeypatch.setattr(runner, "command", fail)
    if optional_missing:
        assert runner.describe("gcloud", "gs://synthetic/index.csi", "billing", optional=True) is None
    else:
        with pytest.raises(subprocess.CalledProcessError):
            runner.describe("gcloud", "gs://synthetic/index.csi", "billing", optional=True)


def test_crc_rejects_same_length_corruption(tmp_path):
    path = tmp_path / "payload"
    path.write_bytes(b"original")
    metadata = {"size": 8, "crc32c": runner.crc32c(path)}
    runner.verify_object(path, metadata)
    path.write_bytes(b"corrupt!")
    with pytest.raises(ValueError, match="CRC32C"):
        runner.verify_object(path, metadata)


def test_valid_cache_skips_download_and_changed_file_is_repaired(tmp_path, monkeypatch):
    cfg, metadata, calls = setup_stage(tmp_path, monkeypatch)
    bcf = runner.stage_one(cfg, "hiphase", metadata)
    first_calls = len(calls)
    assert runner.stage_one(cfg, "hiphase", metadata) == bcf
    assert len(calls) == first_calls
    bcf.write_bytes(b"bad")
    runner.stage_one(cfg, "hiphase", metadata)
    assert len(calls) > first_calls
    runner.verify_object(bcf, metadata["bcf"])
    sync_env = next(env for cmd, env in calls if cmd[1:3] == ["storage", "rsync"])
    assert sync_env == {"CLOUDSDK_STORAGE_PROCESS_COUNT": "1", "CLOUDSDK_STORAGE_THREAD_COUNT": "4"}


def test_missing_csi_is_built_locally(tmp_path, monkeypatch):
    cfg, metadata, calls = setup_stage(tmp_path, monkeypatch, have_csi=False)
    bcf = runner.stage_one(cfg, "hiphase", metadata)
    assert Path(str(bcf) + ".csi").is_file()
    assert any(cmd[1:3] == ["index", "--force"] and "--threads" in cmd for cmd, _ in calls)
    marker = json.loads((bcf.parent / "stage.done.json").read_text())
    assert marker["files"]["csi"]["sha256"] == runner.ph.digest(Path(str(bcf) + ".csi"))


def test_source_generation_change_prevents_completion(tmp_path, monkeypatch):
    cfg, metadata, calls = setup_stage(tmp_path, monkeypatch)
    changed = {**metadata, "bcf": {**metadata["bcf"], "generation": "999"}}
    monkeypatch.setattr(runner, "metadata_for", lambda *a: changed)
    with pytest.raises(ValueError, match="changed during staging"):
        runner.stage_one(cfg, "hiphase", metadata)
    assert not (tmp_path / "inputs/hiphase/stage.done.json").exists()


def test_disk_preflight_prevents_download(tmp_path, monkeypatch):
    cfg, metadata, calls = setup_stage(tmp_path, monkeypatch)
    monkeypatch.setattr(runner.shutil, "disk_usage", lambda *a: argparse.Namespace(free=0))
    with pytest.raises(RuntimeError, match="GiB free"):
        runner.stage(cfg)
    assert calls == []


def test_prepare_and_run_use_local_staging_without_cloud_calls(tmp_path, monkeypatch):
    paths = {"scaffold": tmp_path / "s.bcf", "hiphase": tmp_path / "h.bcf"}
    calls = []
    monkeypatch.setattr(runner.shutil, "which", lambda executable: executable)
    monkeypatch.setattr(runner, "staged_paths", lambda cfg: paths)
    monkeypatch.setattr(runner, "stage", lambda *a: pytest.fail("unexpected cloud staging"))
    monkeypatch.setattr(runner, "command", lambda cmd, **kwargs: calls.append(list(map(str, cmd))) or "test version")
    runner.main(["--root", str(tmp_path), "--phase", "run"])
    assert not any("gcloud" in cmd[0] for cmd in calls)
    assert any("make_singleton_manifest.py" in " ".join(cmd) for cmd in calls)
    assert any("phase_singletons.py" in " ".join(cmd) for cmd in calls)
    assert not (tmp_path / ".runner.lock").exists()
