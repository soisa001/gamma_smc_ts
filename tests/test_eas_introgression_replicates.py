from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import gamma_smc_aou.eas_introgression_replicates as replicates
import gamma_smc_aou.eas_sweep_study as study
import pandas as pd
import pytest


def _runtime(tmp_path, **overrides):
    values = {
        "repo_root": str(tmp_path),
        "sample_diploids": 10,
        "campaign_subdirectory": replicates.DEFAULT_CAMPAIGN,
    }
    values.update(overrides)
    return study.StudyRuntime(**values)


def _base_specifications(tmp_path):
    return study.build_specifications(_runtime(tmp_path), ("introgression",))[
        "introgression"
    ]


def test_expansion_is_replicate_major_unique_and_stable_when_extended(tmp_path):
    base = _base_specifications(tmp_path)

    two = replicates.expand_introgression_specifications(base, 2)
    three = replicates.expand_introgression_specifications(base, 3)

    assert len(two) == 54
    assert two == three[:54]
    assert [row["replicate_index"] for row in two[:27]] == [1] * 27
    assert [row["replicate_index"] for row in two[27:]] == [2] * 27
    assert two[0]["base_specification_id"] == base[0]["specification_id"]
    assert two[0]["replicate_id"] == f"{base[0]['specification_id']}__rep001"
    assert two[0]["specification_id"] == two[0]["replicate_id"]
    assert two[27]["specification_id"] == f"{base[0]['specification_id']}__rep002"
    assert len({row["specification_id"] for row in two}) == 54
    assert (
        study._stable_seed(study.BASE_SEED, two[0]["replicate_id"], 0)  # noqa: SLF001
        != study._stable_seed(  # noqa: SLF001
            study.BASE_SEED, two[27]["replicate_id"], 0
        )
    )
    assert base[0]["specification_id"].endswith("af10")


@pytest.mark.parametrize("value", [0, -1, 1.5, True])
def test_expansion_rejects_invalid_replicate_count(tmp_path, value):
    with pytest.raises(ValueError, match="positive integer"):
        replicates.expand_introgression_specifications(
            _base_specifications(tmp_path), value
        )


def test_plan_records_fixed_contract_hashes_and_is_idempotent(tmp_path):
    runtime = _runtime(tmp_path)

    base, expanded = replicates.write_introgression_replicate_plan(runtime, 10)
    campaign_dir = study._runtime_study_dir(  # noqa: SLF001
        runtime, "introgression"
    )
    table_path = campaign_dir / "specifications.tsv"
    design_path = campaign_dir / "study_design.json"
    first_table = table_path.read_bytes()
    first_design = design_path.read_bytes()

    base_again, expanded_again = replicates.write_introgression_replicate_plan(
        runtime, 10
    )

    assert base_again == base
    assert expanded_again == expanded
    assert table_path.read_bytes() == first_table
    assert design_path.read_bytes() == first_design
    table = pd.read_csv(table_path, sep="\t")
    design = json.loads(design_path.read_text(encoding="utf-8"))
    contract = design["campaign_contract"]
    assert len(base) == 27
    assert len(expanded) == len(table) == 270
    assert contract["replicates_per_biological_specification"] == 10
    assert contract["target_execution_slots"] == 270
    assert contract["han_final_frequency_floors"] == pytest.approx(
        [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    )
    assert (
        contract["acceptance"]["archaic_specific_no_ils"]["han_entry_route"]
        == "introgression_pulse_only"
    )
    assert contract["acceptance"]["genotype_classes"] == {
        "hom_ref_minimum_diploids": 2,
        "heterozygous_minimum_diploids": 1,
        "hom_alt_minimum_diploids": 2,
    }
    assert design["rng"]["derivation"] == (
        "sha256(base_seed:replicate_id:external_draw)"
    )
    assert design["slim_compatibility"]["required_version"] == "4.2.2"
    table_sha256 = study._sha256(table_path)  # noqa: SLF001
    assert design["hashes"]["specifications_tsv_sha256"] == table_sha256
    assert "replicate_count" not in table.columns


def test_immutable_plan_rejects_changed_target_count_and_preserves_work(tmp_path):
    runtime = _runtime(tmp_path, campaign_subdirectory="extensible_campaign")
    _, first = replicates.write_introgression_replicate_plan(runtime, 2)
    campaign_dir = study._runtime_study_dir(  # noqa: SLF001
        runtime, "introgression"
    )
    sentinel = campaign_dir / "work" / first[0]["specification_id"] / "sentinel"
    sentinel.parent.mkdir(parents=True)
    sentinel.write_text("keep", encoding="utf-8")
    plan_bytes = (campaign_dir / "study_design.json").read_bytes()
    table_bytes = (campaign_dir / "specifications.tsv").read_bytes()

    with pytest.raises(ValueError, match="incompatible"):
        replicates.write_introgression_replicate_plan(runtime, 3)

    assert sentinel.read_text(encoding="utf-8") == "keep"
    assert (campaign_dir / "study_design.json").read_bytes() == plan_bytes
    assert (campaign_dir / "specifications.tsv").read_bytes() == table_bytes
    assert (
        (campaign_dir / "work")
        .as_posix()
        .endswith("introgression_EAS_sim/extensible_campaign/work")
    )


def test_existing_plan_allows_execution_knob_changes_without_rewrite(tmp_path):
    initial = _runtime(tmp_path, campaign_subdirectory="execution_knobs")
    replicates.write_introgression_replicate_plan(initial, 10)
    campaign_dir = study._runtime_study_dir(initial, "introgression")  # noqa: SLF001
    plan_path = campaign_dir / "study_design.json"
    table_path = campaign_dir / "specifications.tsv"
    before = (plan_path.read_bytes(), table_path.read_bytes())

    resumed = _runtime(
        tmp_path,
        campaign_subdirectory="execution_knobs",
        simulation_workers=12,
        threads=8,
        max_external_draws=99,
        max_launched_draws=120,
        max_specification_seconds=600,
        max_cumulative_specification_seconds=3_600,
    )
    replicates.write_introgression_replicate_plan(resumed, 10)

    assert (plan_path.read_bytes(), table_path.read_bytes()) == before

    incompatible = _runtime(
        tmp_path,
        campaign_subdirectory="execution_knobs",
        recombination_rate=2e-8,
    )
    with pytest.raises(ValueError, match="recombination_rate"):
        replicates.write_introgression_replicate_plan(incompatible, 10)
    assert (plan_path.read_bytes(), table_path.read_bytes()) == before


def test_existing_plan_fails_closed_on_table_corruption(tmp_path):
    runtime = _runtime(tmp_path, campaign_subdirectory="corrupt_table")
    replicates.write_introgression_replicate_plan(runtime, 10)
    campaign_dir = study._runtime_study_dir(runtime, "introgression")  # noqa: SLF001
    table_path = campaign_dir / "specifications.tsv"
    table_path.write_bytes(table_path.read_bytes() + b"# corruption\n")

    with pytest.raises(ValueError, match="checksum"):
        replicates.write_introgression_replicate_plan(runtime, 10)


def test_existing_plan_fails_closed_when_only_one_plan_file_exists(tmp_path):
    runtime = _runtime(tmp_path, campaign_subdirectory="partial_plan")
    campaign_dir = study._runtime_study_dir(runtime, "introgression")  # noqa: SLF001
    campaign_dir.mkdir(parents=True)
    (campaign_dir / "specifications.tsv").write_text("incomplete\n", encoding="utf-8")

    with pytest.raises(ValueError, match="plan is incomplete"):
        replicates.write_introgression_replicate_plan(runtime, 10)


@pytest.mark.parametrize(
    ("output", "accepted"),
    [
        ("SLiM 4.2.2, built for x86_64", True),
        ("SLiM version 4.2.2", True),
        ("SLiM 5.2.0", False),
        ("SLiM 4.2.1", False),
        ("unknown", False),
    ],
)
def test_slim_version_is_strictly_4_2_2(tmp_path, monkeypatch, output, accepted):
    executable = tmp_path / "slim"
    executable.write_bytes(b"test executable")
    monkeypatch.setattr(
        replicates.subprocess,
        "run",
        lambda *args, **kwargs: SimpleNamespace(stdout=output, stderr="", returncode=0),
    )

    if accepted:
        observed = replicates.verify_slim_version(executable)
        assert observed["version"] == "4.2.2"
        assert observed["sha256"] == study._sha256(executable)  # noqa: SLF001
    else:
        with pytest.raises(RuntimeError, match="requires SLiM 4.2.2"):
            replicates.verify_slim_version(executable)


def test_cli_defaults_are_bounded_and_parallel():
    args = replicates.build_parser().parse_args(["plan"])

    assert args.campaign == "replicate_study_n10"
    assert args.replicates == 10
    assert args.workers == 4
    assert args.threads == 4
    assert args.max_external_draws == 20
    assert args.max_launched_draws == 40
    assert args.spec_timeout_minutes == 5.0
    assert args.cumulative_spec_timeout_minutes == 30.0


@pytest.mark.parametrize("phase", ["plot", "all"])
def test_cli_forbids_filtered_campaign_wide_aggregation(tmp_path, phase):
    with pytest.raises(SystemExit) as error:
        replicates.main(
            [
                phase,
                "--repo-root",
                str(tmp_path),
                "--campaign",
                "filtered_plot",
                "--spec",
                "af10",
            ]
        )

    assert error.value.code == 2
    assert not (tmp_path / study.INTROGRESSION_DIR / "filtered_plot").exists()


def test_cli_phase_invocation_provenance_is_append_only(tmp_path):
    arguments = [
        "plan",
        "--repo-root",
        str(tmp_path),
        "--campaign",
        "logged_plan",
    ]
    assert replicates.main(arguments) == 0
    runtime = _runtime(tmp_path, campaign_subdirectory="logged_plan")
    campaign_dir = study._runtime_study_dir(runtime, "introgression")  # noqa: SLF001
    invocation_dir = campaign_dir / replicates.PHASE_INVOCATION_DIRECTORY
    first_files = sorted(invocation_dir.glob("*.json"))
    first_bytes = {path.name: path.read_bytes() for path in first_files}

    assert replicates.main([*arguments, "--workers", "12"]) == 0

    all_files = sorted(invocation_dir.glob("*.json"))
    assert len(first_files) == 2
    assert len(all_files) == 4
    assert all(
        (invocation_dir / name).read_bytes() == payload
        for name, payload in first_bytes.items()
    )
    events = [json.loads(path.read_text(encoding="utf-8")) for path in all_files]
    assert {record["event"] for record in events} == {"started", "completed"}
    assert all(
        record["schema"] == replicates.PHASE_INVOCATION_SCHEMA for record in events
    )
    invocation_ids = {record["invocation_id"] for record in events}
    assert len(invocation_ids) == 2
    assert all(
        sum(record["invocation_id"] == identifier for record in events) == 2
        for identifier in invocation_ids
    )
    assert all(record["plan"]["specifications_sha256"] for record in events)


def test_cli_records_failed_phase_without_rewriting_plan(tmp_path, monkeypatch):
    campaign = "logged_failure"
    fake_slim = tmp_path / "slim"
    fake_slim.write_bytes(b"not used")
    monkeypatch.setattr(
        replicates,
        "verify_slim_version",
        lambda path: (_ for _ in ()).throw(RuntimeError("version rejected")),
    )

    with pytest.raises(RuntimeError, match="version rejected"):
        replicates.main(
            [
                "simulate",
                "--repo-root",
                str(tmp_path),
                "--campaign",
                campaign,
                "--slim-bin",
                str(fake_slim),
            ]
        )

    runtime = _runtime(tmp_path, campaign_subdirectory=campaign)
    campaign_dir = study._runtime_study_dir(runtime, "introgression")  # noqa: SLF001
    events = [
        json.loads(path.read_text(encoding="utf-8"))
        for path in sorted(
            (campaign_dir / replicates.PHASE_INVOCATION_DIRECTORY).glob("*.json")
        )
    ]
    assert [record["event"] for record in events] == ["started", "failed"]
    assert events[-1]["error"] == {
        "type": "RuntimeError",
        "message": "version rejected",
    }
    assert (campaign_dir / "study_design.json").is_file()
    assert (campaign_dir / "specifications.tsv").is_file()


def test_recover_requires_an_aware_stop_time_before_creating_campaign(tmp_path):
    with pytest.raises(SystemExit) as missing:
        replicates.main(
            [
                "recover",
                "--repo-root",
                str(tmp_path),
                "--campaign",
                "missing_stop",
            ]
        )
    assert missing.value.code == 2
    assert not (tmp_path / study.INTROGRESSION_DIR / "missing_stop").exists()

    with pytest.raises(ValueError, match="must include a UTC offset"):
        replicates.main(
            [
                "recover",
                "--repo-root",
                str(tmp_path),
                "--campaign",
                "naive_stop",
                "--interruption-stop-utc",
                "2026-08-12T12:00:00",
            ]
        )
    assert not (tmp_path / study.INTROGRESSION_DIR / "naive_stop").exists()


def test_recover_closes_external_interruptions_and_is_idempotent(tmp_path, monkeypatch):
    campaign = "recovery_campaign"
    runtime = study.StudyRuntime(
        repo_root=str(tmp_path), campaign_subdirectory=campaign
    )
    _, expanded = replicates.write_introgression_replicate_plan(runtime, 10)
    campaign_dir = study._runtime_study_dir(runtime, "introgression")  # noqa: SLF001
    invocation_dir = campaign_dir / replicates.PHASE_INVOCATION_DIRECTORY
    invocation_dir.mkdir(parents=True)
    prior_id = "prior-simulate-invocation"
    prior_context = {
        "phase": "simulate",
        "phase_actions": ["simulate"],
        "campaign": campaign,
        "operator_context": "must be preserved",
    }
    prior_started = {
        "schema": replicates.PHASE_INVOCATION_SCHEMA,
        "invocation_id": prior_id,
        "event": "started",
        "recorded_at_utc": "2026-08-12T12:00:00+00:00",
        "process_id": 12345,
        **prior_context,
    }
    (invocation_dir / f"20260812T120000.000000Z_{prior_id}_started.json").write_text(
        json.dumps(prior_started), encoding="utf-8"
    )

    specification = expanded[0]
    specification_dir = campaign_dir / "work" / specification["specification_id"]
    specification_dir.mkdir(parents=True)
    lock_path = specification_dir / ".simulation.lock"
    lock_path.write_text("{}", encoding="utf-8")
    recorded_slim = tmp_path / "recorded" / "slim"
    (specification_dir / "simulation_contract.json").write_text(
        json.dumps({"software": {"slim_path": str(recorded_slim)}}),
        encoding="utf-8",
    )
    recovered_tasks = []

    def fake_recover(task, *, stopped_at_utc):
        recovered_tasks.append((task, stopped_at_utc))
        Path(task["specification_dir"], ".simulation.lock").unlink()
        return {
            "status": "failed",
            "failure": {
                "elapsed_seconds": 300.0,
                "error_type": "ExternalProcessInterruption",
            },
        }

    monkeypatch.setattr(study, "record_external_simulation_interruption", fake_recover)
    monkeypatch.setattr(
        replicates,
        "verify_slim_version",
        lambda path: pytest.fail("recover must not execute or version-check SLiM"),
    )
    arguments = [
        "recover",
        "--repo-root",
        str(tmp_path),
        "--campaign",
        campaign,
        "--interruption-stop-utc",
        "2026-08-12T12:05:00Z",
        "--interruption-note",
        "desktop task was interrupted",
    ]
    assert replicates.main(arguments) == 0
    assert len(recovered_tasks) == 1
    recovered_task, recovered_stop = recovered_tasks[0]
    assert (
        recovered_task["specification"]["specification_id"]
        == specification["specification_id"]
    )
    assert Path(recovered_task["slim_path"]) == recorded_slim.resolve()
    assert recovered_stop == "2026-08-12T12:05:00+00:00"
    assert not lock_path.exists()

    events = [
        json.loads(path.read_text(encoding="utf-8"))
        for path in sorted(invocation_dir.glob("*.json"))
    ]
    prior_failed = next(
        event
        for event in events
        if event["invocation_id"] == prior_id and event["event"] == "failed"
    )
    assert prior_failed["elapsed_seconds"] == 300.0
    assert prior_failed["error"]["type"] == "ExternalProcessInterruption"
    assert {key: prior_failed[key] for key in prior_context} == prior_context
    recovery_completed = [
        event
        for event in events
        if event["phase"] == "recover" and event["event"] == "completed"
    ][-1]
    summary = recovery_completed["tool_provenance"]["recovery"]
    assert summary["recovered_simulation_lock_count"] == 1
    assert summary["phase_invocations"]["closed_invocation_ids"] == [prior_id]
    assert summary["no_op"] is False

    assert replicates.main(arguments) == 0
    assert len(recovered_tasks) == 1
    events = [
        json.loads(path.read_text(encoding="utf-8"))
        for path in sorted(invocation_dir.glob("*.json"))
    ]
    recovery_completed = [
        event
        for event in events
        if event["phase"] == "recover" and event["event"] == "completed"
    ]
    assert len(recovery_completed) == 2
    assert recovery_completed[-1]["tool_provenance"]["recovery"]["no_op"] is True
    assert (
        sum(
            event["invocation_id"] == prior_id and event["event"] == "failed"
            for event in events
        )
        == 1
    )


def test_phase_event_context_cannot_override_reserved_fields(tmp_path):
    with pytest.raises(ValueError, match="reserved event fields"):
        replicates._append_phase_invocation_event(  # noqa: SLF001
            tmp_path,
            invocation_id="invocation",
            event="started",
            context={"event": "completed"},
        )
