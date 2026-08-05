from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest


def _runner_command(repo: Path, arguments: list[str]) -> tuple[list[str], Path]:
    if os.name == "nt":
        git_bash = (
            Path(os.environ.get("ProgramFiles", r"C:\Program Files"))
            / "Git/usr/bin/bash.exe"
        )
        if not git_bash.is_file():
            pytest.skip("Git Bash is not available")
        command = "bash scripts/run_aou_workbench.sh " + " ".join(arguments)
        return [str(git_bash), "-lc", command], repo
    bash = shutil.which("bash")
    if bash is None:
        pytest.skip("bash is not available")
    return [bash, str(repo / "scripts" / "run_aou_workbench.sh"), *arguments], repo


def test_runner_dry_run_resolves_case_insensitive_defaults():
    repo = Path(__file__).resolve().parents[1]
    environment = os.environ.copy()
    environment["WORKSPACE_BUCKET"] = "gs://test-workspace"
    environment["GOOGLE_PROJECT"] = "test-billing-project"
    command, cwd = _runner_command(repo, ["-chr", "1", "-pops", "afr", "--dry-run"])
    completed = subprocess.run(
        command,
        check=True,
        text=True,
        capture_output=True,
        env=environment,
        cwd=cwd,
    )
    output = completed.stdout
    assert "populations: AFR" in output
    assert "requester-pays billing project: test-billing-project" in output
    assert "chromosomes: 1" in output
    assert "threads=12" in output
    assert "recent_call=mean" in output
    assert "stride=10000 bp, cache=1000 bp" in output
    assert "100000 random haplotype pairs/pop, seed=1729, exclude_within=0" in output
    assert "candidates: fraction>0.02, merge_gap=20000 bp" in output
    assert "aou_lr_phase2_v1.chr1.bubble.split.bcf" in output
    assert "ancestry_preds.tsv (column ancestry_pred_other)" in output
    assert "flagged_samples.tsv" in output
    assert "relatedness_flagged_samples.tsv" in output
    assert "hardmask.hg38.v4.over99.bed" in output
    assert "mask mode: default (excluded_intervals)" in output
    assert "gencode.v50.basic.annotation.gtf.gz" in output
    assert (
        "merge_gap=1000000 bp (display only), gene_flank=+/-500000 bp, zoom_ymax=0.04, "
        "label_min=0.02"
    ) in output


def test_runner_has_a_workbench_output_bucket_default():
    repo = Path(__file__).resolve().parents[1]
    environment = os.environ.copy()
    environment.pop("WORKSPACE_BUCKET", None)
    environment.pop("AOU_GAMMA_OUTPUT_PREFIX", None)
    command, cwd = _runner_command(repo, ["-chr", "1", "-pops", "AFR", "--dry-run"])
    completed = subprocess.run(
        command,
        check=True,
        text=True,
        capture_output=True,
        env=environment,
        cwd=cwd,
    )
    assert (
        "gs://rw-migration-aou-rw-fa99430f/gamma_smc/results/AFR/{chromosomes,plots}/"
    ) in completed.stdout


def test_strict_hardmask_uses_callable_semantics_and_separate_roots():
    repo = Path(__file__).resolve().parents[1]
    environment = os.environ.copy()
    environment["WORKSPACE_BUCKET"] = "gs://test-workspace"
    environment["GOOGLE_PROJECT"] = "test-billing-project"
    command, cwd = _runner_command(
        repo,
        ["-chr", "1", "-pops", "AFR", "--strict-hardmask", "--dry-run"],
    )
    completed = subprocess.run(
        command,
        check=True,
        text=True,
        capture_output=True,
        env=environment,
        cwd=cwd,
    )
    output = completed.stdout
    assert "mask mode: strict (included_intervals)" in output
    assert "/home/jupyter/gamma_smc_workbench_strict_hardmask" in output
    assert "gs://test-workspace/gamma_smc/results_strict_hardmask/AFR/" in output
    assert "gs://test-workspace/hmmix-static/hg38_strick_callability_mask.bed" in output


def test_both_hardmask_runs_resolve_default_and_strict_namespaces():
    repo = Path(__file__).resolve().parents[1]
    environment = os.environ.copy()
    environment["WORKSPACE_BUCKET"] = "gs://test-workspace"
    environment["GOOGLE_PROJECT"] = "test-billing-project"
    command, cwd = _runner_command(
        repo,
        ["-chr", "1", "-pops", "AFR", "--both-hardmask", "--dry-run"],
    )
    completed = subprocess.run(
        command,
        check=True,
        text=True,
        capture_output=True,
        env=environment,
        cwd=cwd,
    )
    output = completed.stdout
    assert output.count("Gamma-SMC Workbench plan") == 2
    assert "mask mode: default (excluded_intervals)" in output
    assert "mask mode: strict (included_intervals)" in output
    assert "/home/jupyter/gamma_smc_workbench\n" in output
    assert "/home/jupyter/gamma_smc_workbench_strict_hardmask" in output


def test_runner_requires_explicit_scope():
    repo = Path(__file__).resolve().parents[1]
    environment = os.environ.copy()
    environment["WORKSPACE_BUCKET"] = "gs://test-workspace"
    command, cwd = _runner_command(repo, ["--dry-run"])
    completed = subprocess.run(
        command,
        text=True,
        capture_output=True,
        env=environment,
        cwd=cwd,
    )
    assert completed.returncode == 2
    assert "-chr is required" in completed.stderr


def test_all_gcloud_storage_calls_include_requester_pays_billing():
    repo = Path(__file__).resolve().parents[1]
    runner = (repo / "scripts/run_aou_workbench.sh").read_text()

    assert "AOU_GAMMA_BILLING_PROJECT:-${GOOGLE_PROJECT:-}" in runner
    assert runner.count("gcloud storage ") == runner.count(
        '--billing-project "$BILLING_PROJECT"'
    )


def test_runner_recovers_local_outputs_before_destructive_restart():
    repo = Path(__file__).resolve().parents[1]
    runner = (repo / "scripts/run_aou_workbench.sh").read_text()

    recovery = 'echo "Recovered complete local outputs; skipping decode."'
    destructive_restart = 'safe_clear_run "$population_summary_dir"'
    assert recovery in runner
    assert runner.index(recovery) < runner.index(destructive_restart)


def test_resume_helpers_do_not_expand_locals_during_their_declaration():
    repo = Path(__file__).resolve().parents[1]
    runner = (repo / "scripts/run_aou_workbench.sh").read_text()

    assert 'destination="$2" temporary="${destination}' not in runner
    assert 'local_file="$2" temporary="${local_file}' not in runner
    assert 'local temporary="${destination}.partial.$$"' in runner
    assert 'local temporary="${local_file}.remote.$$"' in runner


def test_runner_locks_reports_and_checksum_syncs_outputs():
    repo = Path(__file__).resolve().parents[1]
    runner = (repo / "scripts/run_aou_workbench.sh").read_text()

    assert 'flock -n 9 || die "another Workbench runner' in runner
    assert "workbench-report" in runner
    assert '--gene-annotation "$local_gene_annotation"' in runner
    assert '--gene-label-overrides "$GENE_LABEL_OVERRIDES"' in runner
    assert '--plot-merge-gap "$PLOT_MERGE_GAP"' in runner
    assert '--hit-label-min-fraction "$HIT_LABEL_MIN_FRACTION"' in runner
    assert "gene_list.tsv" in (repo / "python/gamma_smc_aou/workbench.py").read_text()
    assert (
        "raw_scan_windows.tsv.gz"
        in (repo / "python/gamma_smc_aou/workbench.py").read_text()
    )
    assert 'gcloud storage rsync "$source_directory" "$remote_directory"' in runner
    assert "--recursive --checksums-only" in runner
    assert "--delete-unmatched-destination-objects" not in runner
    assert 'rsync_single_file "$population_plot_dir/plot_manifest.json"' in runner
    assert 'rsync_single_file "$report_dir/run_report_manifest.json"' in runner
    assert '"$AOU" --sync-only' in runner
    assert "export AOU_UV_NO_SYNC=1" in runner
    aou = (repo / "scripts/aou.sh").read_text()
    assert "run_args+=(--no-sync)" in aou
    assert (
        "regions_by_population.tsv"
        in (repo / "python/gamma_smc_aou/workbench.py").read_text()
    )


def test_runner_preserves_chromosome_candidate_directory_in_cloud():
    repo = Path(__file__).resolve().parents[1]
    runner = (repo / "scripts/run_aou_workbench.sh").read_text()

    expected = (
        '"$remote_directory/$(basename "$candidate_directory")/'
        '$(basename "$candidate_manifest")"'
    )
    assert runner.count(expected) == 2
    assert '"$remote_directory/candidates/' not in runner
    assert 'cp -al -- "$candidate_directory" "$stage/"' in runner


def test_runner_syncs_completion_markers_after_payloads():
    repo = Path(__file__).resolve().parents[1]
    runner = (repo / "scripts/run_aou_workbench.sh").read_text()

    function = runner.split("upload_completed_chromosome() {", 1)[1].split(
        "\nsafe_clear_run()", 1
    )[0]
    assert function.index(
        'rsync_directory "$stage" "$remote_directory"'
    ) < function.index('rsync_single_file "$completion" "$remote_directory"')
