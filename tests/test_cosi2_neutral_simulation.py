from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from gamma_smc_aou import cosi2_neutral_simulation as neutral

REPO_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def plan_bundle(tmp_path_factory: pytest.TempPathFactory) -> Path:
    output = tmp_path_factory.mktemp("cosi2-neutral") / "plan"
    result = neutral.generate_neutral_plan(REPO_ROOT, output_dir=output)
    assert result["status"] == "complete"
    return output


def _unit(plan_bundle: Path, model_id: str, replicate_index: int = 0) -> dict:
    frame = pd.read_csv(plan_bundle / neutral.PLAN_FILENAME, sep="\t")
    row = frame[
        (frame["model_id"] == model_id) & (frame["replicate_index"] == replicate_index)
    ]
    assert len(row) == 1
    return row.iloc[0].to_dict()


def _trajectory_frame(unit: dict, *, endpoint: float = 0.03) -> pd.DataFrame:
    model_id = str(unit["model_id"])
    birth = int(unit["mutation_birth_generations_ago"])
    generations = np.arange(birth, -1, -1, dtype=int)
    if model_id == neutral.EAS_MODEL_ID:
        size = np.full(len(generations), 3_197.0)
        frequency = np.linspace(1 / (2 * size[0]), endpoint, len(generations))
        return pd.DataFrame(
            {
                "sim": 1,
                "gen": generations,
                "selfreq_1": frequency,
                "popsize_1": size,
            }
        )
    sizes = {
        identifier: np.full(len(generations), 3_600.0) for identifier in range(1, 6)
    }
    frequencies = {
        identifier: np.zeros(len(generations), dtype=float)
        for identifier in range(1, 6)
    }
    frequencies[4] = np.linspace(1 / 7_200, 0.2, len(generations))
    frequencies[3] = np.linspace(0.0, endpoint, len(generations))
    boundary_index = int(np.flatnonzero(generations == 645)[0])
    frequencies[3][boundary_index] = 0.40  # deliberately outside the old 1.1%-1.3% gate
    data: dict[str, np.ndarray | int] = {"sim": 1, "gen": generations}
    for identifier in range(1, 6):
        data[f"selfreq_{identifier}"] = frequencies[identifier]
    for identifier in range(1, 6):
        data[f"popsize_{identifier}"] = sizes[identifier]
    return pd.DataFrame(data)


def _write_trajectory(path: Path, unit: dict, *, endpoint: float = 0.03) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    _trajectory_frame(unit, endpoint=endpoint).to_csv(
        path, sep="\t", index=False, lineterminator="\n"
    )


def _ms_text(
    unit: dict,
    *,
    positions: tuple[float, ...] = (0.2, 0.8),
    haplotypes: list[str] | None = None,
) -> str:
    count = len(positions)
    if haplotypes is None:
        haplotypes = ["0" * count] * int(unit["sample_haploids"])
    return (
        f"ms {int(unit['sample_haploids'])} 1\n"
        f"cosi_rand {int(unit['seed'])}\n\n"
        f"// seed={int(unit['seed'])}\n"
        f"segsites: {count}\n"
        + "positions: "
        + " ".join(f"{value:.12g}" for value in positions)
        + "\n"
        + "\n".join(haplotypes)
        + "\n\n"
    )


def test_configs_are_s0_survival_only_and_origin_age_matched():
    cells = neutral.build_neutral_model_cells(REPO_ROOT).set_index("model_id")
    assert (
        int(cells.loc[neutral.EAS_MODEL_ID, "mutation_birth_generations_ago"]) == 1475
    )
    assert (
        int(cells.loc[neutral.HAN_MODEL_ID, "mutation_birth_generations_ago"]) == 2400
    )
    assert (
        int(cells.loc[neutral.HAN_MODEL_ID, "selection_onset_generations_ago"]) == 645
    )
    texts = neutral.render_neutral_parameter_files(REPO_ROOT)
    expected = {
        neutral.EAS_MODEL_ID: (
            'pop_event sweep_mult_standing "eas_neutral_survival" 1 1475 0 '
            "0.5 0.000000000001-1 1 1475"
        ),
        neutral.HAN_MODEL_ID: (
            'pop_event sweep_mult_standing "han_introgressed_neutral_survival" '
            "4 2400 0 0.5 0.000000000001-1 3 645"
        ),
    }
    for model_id, text in texts.items():
        assert expected[model_id] in text
        assert "1e-12-1" not in text  # invalid under CoSi2 v2.4 ValRange parsing
        assert "0.195-0.205" not in text
        assert "0.011-0.013" not in text
        assert "applied post hoc" not in text
        assert "random_seed 0" in text
        audit = neutral.lint_neutral_parameter_text(
            text, expected_cell=cells.loc[model_id].to_dict()
        )
        assert audit["survival_only"] is True
        assert audit["han_intermediate_af_conditioned"] is False


def test_plan_predeclares_exact_disjoint_100_plus_100_seed_streams(plan_bundle: Path):
    frame = pd.read_csv(plan_bundle / neutral.PLAN_FILENAME, sep="\t")
    assert len(frame) == 200
    assert frame.groupby("model_id").size().to_dict() == {
        neutral.EAS_MODEL_ID: 100,
        neutral.HAN_MODEL_ID: 100,
    }
    eas = frame[frame["model_id"] == neutral.EAS_MODEL_ID]
    han = frame[frame["model_id"] == neutral.HAN_MODEL_ID]
    np.testing.assert_array_equal(eas["seed"], neutral.EAS_SEED_START + np.arange(100))
    np.testing.assert_array_equal(han["seed"], neutral.HAN_SEED_START + np.arange(100))
    assert frame["seed"].is_unique
    assert not set(frame["seed"]).intersection(neutral.SELECTED_SEEDS_EXCLUDED)
    assert not frame["adaptive_stopping"].any()
    assert not frame["sample_af_conditioned"].any()
    assert not frame["han_intermediate_af_conditioned"].any()
    manifest = json.loads((plan_bundle / neutral.MANIFEST_FILENAME).read_text())
    assert manifest["seed_contract"]["formulas"] == neutral.SEED_FORMULAS
    assert manifest["all_predeclared_units_are_retained"] is True
    assert manifest["adaptive_stopping"] is False


def test_plan_generation_is_idempotent_and_checksum_bound(plan_bundle: Path):
    first = neutral.verify_neutral_plan(REPO_ROOT, plan_dir=plan_bundle)
    second = neutral.generate_neutral_plan(REPO_ROOT, output_dir=plan_bundle)
    assert second["status"] == "verified_cached"
    assert second["plan_tsv_sha256"] == first["plan_tsv_sha256"]
    assert second["plan_manifest_sha256"] == first["plan_manifest_sha256"]
    assert second["plan_completion_sha256"] == first["plan_completion_sha256"]


def test_trajectory_survival_does_not_gate_han_boundary_or_sample_carriers(
    tmp_path: Path, plan_bundle: Path
):
    eas = _unit(plan_bundle, neutral.EAS_MODEL_ID)
    eas_path = tmp_path / "eas.tsv"
    _write_trajectory(eas_path, eas, endpoint=1e-6)
    eas_audit = neutral.validate_neutral_trajectory(eas_path, unit=eas)
    assert eas_audit["present_population_frequency"] == pytest.approx(1e-6)
    assert eas_audit["sample_carrier_conditioned"] is False

    han = _unit(plan_bundle, neutral.HAN_MODEL_ID)
    han_path = tmp_path / "han.tsv"
    _write_trajectory(han_path, han, endpoint=2e-5)
    han_audit = neutral.validate_neutral_trajectory(han_path, unit=han)
    assert han_audit["present_population_frequency"] == pytest.approx(2e-5)
    assert han_audit["han_generation_645_chb_af_descriptive"] == pytest.approx(0.4)
    assert han_audit["han_intermediate_af_conditioned"] is False

    _write_trajectory(eas_path, eas, endpoint=0.0)
    with pytest.raises(ValueError, match="does not survive"):
        neutral.validate_neutral_trajectory(eas_path, unit=eas)


def test_ms_validation_does_not_require_a_sampled_focal_carrier_and_is_exact(
    tmp_path: Path, plan_bundle: Path
):
    unit = _unit(plan_bundle, neutral.EAS_MODEL_ID)
    no_focal = tmp_path / "no_focal.ms"
    no_focal.write_text(_ms_text(unit), encoding="utf-8")
    audit = neutral.validate_cosi_ms(no_focal, unit=unit)
    assert audit["reported_focal_position_count"] == 0
    assert audit["focal_alt_count_if_unambiguous"] is None
    assert audit["sample_carrier_conditioned"] is False

    near_focal = tmp_path / "near_focal.ms"
    near_focal.write_text(_ms_text(unit, positions=(0.2, 0.500004)), encoding="utf-8")
    near_audit = neutral.validate_cosi_ms(near_focal, unit=unit)
    assert near_audit["reported_focal_position_count"] == 0

    exact_focal = tmp_path / "exact_focal.ms"
    exact_focal.write_text(_ms_text(unit, positions=(0.2, 0.5)), encoding="utf-8")
    exact_audit = neutral.validate_cosi_ms(exact_focal, unit=unit)
    assert exact_audit["reported_focal_position_count"] == 1

    terminal = tmp_path / "terminal.ms"
    terminal.write_text(_ms_text(unit, positions=(0.2, 1.0)), encoding="utf-8")
    terminal_audit = neutral.validate_cosi_ms(terminal, unit=unit)
    assert terminal_audit["status"] == "valid"
    assert exact_audit["focal_alt_count_if_unambiguous"] == 0

    forbidden = tmp_path / "forbidden.ms"
    forbidden.write_text(
        _ms_text(unit).replace("positions:", "muttimes: 1 2\npositions:"),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="forbidden"):
        neutral.validate_cosi_ms(forbidden, unit=unit)


def test_dead_lock_is_archived_but_live_lock_is_never_stolen(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    unit_id = "eas_neutral_survival_r000"
    lock = tmp_path / "locks" / f"{unit_id}.lock"
    lock.parent.mkdir(parents=True)
    lock.write_text(
        json.dumps(
            {
                "unit_id": unit_id,
                "pid": 999_999_999,
                "runtime": neutral._runtime_identity(),
                "started_utc": "2000-01-01T00:00:00Z",
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr(neutral, "_pid_is_live", lambda _pid: False)
    descriptor, acquired = neutral._acquire_lock(tmp_path, unit_id)
    os.close(descriptor)
    acquired.unlink()
    archived = list((tmp_path / "lock_provenance" / unit_id).glob("*.lock"))
    audits = list((tmp_path / "lock_provenance" / unit_id).glob("*.json"))
    assert len(archived) == len(audits) == 1
    assert json.loads(audits[0].read_text())["status"] == (
        "reclaimed_verified_dead_lock"
    )

    lock.write_text(
        json.dumps(
            {
                "unit_id": unit_id,
                "pid": os.getpid(),
                "runtime": neutral._runtime_identity(),
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr(neutral, "_pid_is_live", lambda _pid: True)
    with pytest.raises(RuntimeError, match="live or unverifiable"):
        neutral._acquire_lock(tmp_path, unit_id)
    assert lock.is_file()


def test_plan_and_output_validation_reject_symlinks(tmp_path: Path, plan_bundle: Path):
    if not hasattr(os, "symlink"):
        pytest.skip("symlinks unavailable")
    copied = tmp_path / "plan"
    shutil.copytree(plan_bundle, copied)
    config = copied / "configs" / f"{neutral.EAS_MODEL_ID}.par"
    target = copied / "configs" / "real.par"
    config.replace(target)
    try:
        config.symlink_to(target.name)
    except OSError as error:
        pytest.skip(f"symlink creation unavailable: {error}")
    with pytest.raises(ValueError, match="symlink"):
        neutral.verify_neutral_plan(REPO_ROOT, plan_dir=copied)

    unit = _unit(plan_bundle, neutral.EAS_MODEL_ID)
    target_ms = tmp_path / "target.ms"
    target_ms.write_text(_ms_text(unit), encoding="utf-8")
    linked_ms = tmp_path / "linked.ms"
    linked_ms.symlink_to(target_ms.name)
    with pytest.raises(ValueError, match="symlink"):
        neutral.validate_cosi_ms(linked_ms, unit=unit)


def test_atomic_runner_cache_status_collect_and_exact_completion_contract(
    tmp_path: Path,
    plan_bundle: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    unit = _unit(plan_bundle, neutral.EAS_MODEL_ID)
    binary = tmp_path / "coalescent"
    binary.write_bytes(b"fake-cosi2-binary\n")
    run_root = tmp_path / "runs"
    calls: list[list[str]] = []

    def fake_run(command, *, cwd, env, stdout, stderr, check, timeout):
        del cwd, check, timeout
        calls.append(list(command))
        assert command[-1] == "-m"
        assert "-e" not in command and "-M" not in command
        assert env["COSI_MAXATTEMPTS"] == "1234"
        trajectory = Path(env["COSI_SAVE_TRAJ"])
        _write_trajectory(trajectory, unit, endpoint=0.017)
        stdout.write(_ms_text(unit).encode("utf-8"))
        stderr.write(b"fake deterministic CoSi2 run\n")
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr(neutral.subprocess, "run", fake_run)
    first = neutral.run_neutral_simulations(
        REPO_ROOT,
        cosi2_binary=binary,
        plan_dir=plan_bundle,
        run_dir=run_root,
        max_workers=1,
        max_attempts=1234,
        unit_ids=[str(unit["unit_id"])],
    )
    assert first.iloc[0]["status"] == "complete"
    assert len(calls) == 1
    second = neutral.run_neutral_simulations(
        REPO_ROOT,
        cosi2_binary=binary,
        plan_dir=plan_bundle,
        run_dir=run_root,
        max_workers=1,
        max_attempts=1234,
        unit_ids=[str(unit["unit_id"])],
    )
    assert second.iloc[0]["status"] == "cached"
    assert len(calls) == 1

    statuses = neutral.neutral_simulation_status(
        REPO_ROOT,
        plan_dir=plan_bundle,
        run_dir=run_root,
        cosi2_binary=binary,
    )
    assert len(statuses) == 200
    assert (statuses["status"] == "complete").sum() == 1
    metadata = neutral.collect_neutral_metadata(
        REPO_ROOT,
        plan_dir=plan_bundle,
        run_dir=run_root,
        cosi2_binary=binary,
    )
    assert len(metadata) == 1
    assert metadata.iloc[0]["present_population_frequency"] == pytest.approx(0.017)
    overall = neutral.verify_neutral_simulations(
        REPO_ROOT,
        plan_dir=plan_bundle,
        run_dir=run_root,
        cosi2_binary=binary,
        require_all=False,
    )
    assert overall["status"] == "partial"
    assert overall["validated_complete_units"] == 1
    assert overall["all_predeclared_units_retained"] is True

    completion_path = (
        run_root / "units" / str(unit["unit_id"]) / neutral.UNIT_COMPLETION_FILENAME
    )
    original = json.loads(completion_path.read_text())
    tampered = dict(original)
    tampered["command"] = [*original["command"], "--unexpected"]
    completion_path.write_text(json.dumps(tampered) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="exact run contract"):
        neutral.verify_neutral_unit(
            plan_bundle,
            run_root,
            unit=unit,
            cosi2_binary_sha256=neutral.sha256_file(binary),
        )

    with pytest.raises(ValueError, match="max_workers"):
        neutral.run_neutral_simulations(
            REPO_ROOT,
            cosi2_binary=binary,
            plan_dir=plan_bundle,
            run_dir=run_root,
            max_workers=21,
            unit_ids=[str(unit["unit_id"])],
        )
