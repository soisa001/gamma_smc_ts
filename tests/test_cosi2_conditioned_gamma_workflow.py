from concurrent.futures import Future
import os
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from gamma_smc_aou import cosi2_conditioned_gamma_workflow as workflow
from gamma_smc_aou import cosi2_gamma_analysis as analysis


def _synthetic_base_manifest() -> dict:
    return {
        "schema": "synthetic-base-manifest",
        "manifest_payload_sha256": "1" * 64,
        "base_failure_diagnostics": {
            "pilot_failed_quality_gates": ["global_weight_sum_over_max"],
            "v2_failed_quality_gates": ["global_weight_sum_over_max"],
            "proposal_bank_binding_sha256": "2" * 64,
        },
    }


def _synthetic_extension_plan() -> dict:
    return {
        "schema": "synthetic-extension-plan",
        "plan_payload_sha256": "3" * 64,
    }


def _stub_plan_inputs(monkeypatch) -> None:
    monkeypatch.setattr(
        workflow.conditioned,
        "build_plan_tables",
        lambda **_kwargs: {
            "eas_units": pd.DataFrame({"unit_id": ["synthetic_eas"], "seed": [101]}),
            "han_candidates": pd.DataFrame(
                {"candidate_id": ["synthetic_han"], "seed": [202]}
            ),
        },
    )
    monkeypatch.setattr(
        workflow.conditioned,
        "render_han_parameter_text",
        lambda _root, *, screen: f"synthetic_screen={str(screen).lower()}\n",
    )
    monkeypatch.setattr(
        workflow.conditioned,
        "validate_eas_base_bank_manifest",
        lambda *_args, **_kwargs: {"status": "valid"},
    )
    monkeypatch.setattr(
        workflow.conditioned,
        "validate_eas_extension_plan",
        lambda *_args, **_kwargs: {"status": "valid"},
    )


def _write_synthetic_legacy_v2_plan(
    monkeypatch,
    *,
    repo_root: Path,
    work_root: Path,
) -> tuple[Path, dict, dict]:
    _stub_plan_inputs(monkeypatch)
    manifest = _synthetic_base_manifest()
    extension = _synthetic_extension_plan()
    legacy = work_root / workflow.LEGACY_V2_PLAN_DIRNAME
    legacy.mkdir(parents=True)
    payloads = workflow._plan_payloads(
        repo_root,
        eas_base_bank_manifest=manifest,
        eas_extension_plan=extension,
    )
    row_counts = {}
    for label, filename in workflow.LEGACY_V2_PLAN_OUTPUTS.items():
        payload = payloads[label]
        path = legacy / filename
        if isinstance(payload, pd.DataFrame):
            workflow._atomic_frame(path, payload)
            row_counts[label] = len(payload)
        elif isinstance(payload, dict):
            workflow._atomic_json(path, payload)
        else:
            workflow._atomic_text(path, payload)
    outputs = {
        label: workflow._output_record(
            legacy / filename,
            legacy,
            rows=row_counts.get(label),
        )
        for label, filename in workflow.LEGACY_V2_PLAN_OUTPUTS.items()
    }
    contract = workflow._eas_bank_plan_contract(
        manifest,
        extension,
        source_binding=workflow._legacy_v2_source_binding(repo_root),
    )
    workflow._atomic_json(
        legacy / workflow.PLAN_COMPLETION_FILENAME,
        {
            "schema": workflow.LEGACY_V2_PLAN_SCHEMA,
            "status": "complete",
            "contract": contract,
            "contract_sha256": workflow._canonical_sha256(contract),
            "outputs": outputs,
            "output_count": len(outputs),
        },
    )
    monkeypatch.setattr(
        workflow,
        "LEGACY_V2_PLAN_COMPLETION_SHA256",
        workflow.sha256_file(legacy / workflow.PLAN_COMPLETION_FILENAME),
    )
    return legacy, manifest, extension


def _write_synthetic_legacy_v3_plan(
    monkeypatch,
    *,
    repo_root: Path,
    work_root: Path,
) -> tuple[Path, dict, dict, dict]:
    legacy_v2, manifest, extension = _write_synthetic_legacy_v2_plan(
        monkeypatch,
        repo_root=repo_root,
        work_root=work_root,
    )
    proposal_root = work_root / workflow.EAS_PROPOSAL_DIRNAME
    v2 = workflow._verify_legacy_v2_plan(
        repo_root,
        legacy_plan_dir=legacy_v2,
        proposal_dir=proposal_root,
        current_replay_workflow_sha256=workflow.LEGACY_V3_WORKFLOW_SHA256,
    )
    replay_plan = {
        "schema": workflow.LEGACY_V3_PARALLEL_REPLAY_PLAN_SCHEMA,
        "status": "complete",
        "execution_method": workflow.EAS_PARALLEL_REPLAY_METHOD,
        "source_binding": {
            "workflow_sha256": workflow.LEGACY_V3_WORKFLOW_SHA256,
            "core_sha256": workflow.sha256_file(workflow.conditioned.MODULE_PATH),
        },
        "legacy_v2_plan_attestation": v2["attestation"],
    }
    replay_plan["plan_payload_sha256"] = workflow._canonical_sha256(replay_plan)
    legacy_v3 = work_root / workflow.LEGACY_V3_PLAN_DIRNAME
    legacy_v3.mkdir(parents=True)
    for label, filename in workflow.LEGACY_V2_PLAN_OUTPUTS.items():
        (legacy_v3 / filename).write_bytes((legacy_v2 / filename).read_bytes())
    workflow._atomic_json(
        legacy_v3 / workflow.LEGACY_V3_PLAN_OUTPUTS["eas_parallel_replay_plan"],
        replay_plan,
    )
    outputs = {
        label: workflow._output_record(legacy_v3 / filename, legacy_v3)
        for label, filename in workflow.LEGACY_V3_PLAN_OUTPUTS.items()
    }
    for label in ("eas_units", "eas_proposal_blocks", "han_candidates"):
        outputs[label]["rows"] = len(
            pd.read_csv(legacy_v3 / workflow.LEGACY_V3_PLAN_OUTPUTS[label], sep="\t")
        )
    contract = workflow._eas_bank_plan_contract(
        manifest,
        extension,
        source_binding=workflow._plan_source_binding(
            repo_root,
            workflow_sha256=workflow.LEGACY_V3_WORKFLOW_SHA256,
        ),
    )
    contract.update(
        {
            "legacy_v2_plan_attestation": v2["attestation"],
            "eas_parallel_replay_plan_payload_sha256": replay_plan[
                "plan_payload_sha256"
            ],
        }
    )
    workflow._atomic_json(
        legacy_v3 / workflow.PLAN_COMPLETION_FILENAME,
        {
            "schema": workflow.LEGACY_V3_PLAN_SCHEMA,
            "status": "complete",
            "contract": contract,
            "contract_sha256": workflow._canonical_sha256(contract),
            "outputs": outputs,
            "output_count": len(outputs),
        },
    )
    legacy_replay_root = work_root / workflow.LEGACY_V3_PARALLEL_REPLAY_DIRNAME
    legacy_replay_root.mkdir()
    monkeypatch.setattr(workflow, "LEGACY_V3_REPLAY_BLOCK_INDEXES", ())
    monkeypatch.setattr(
        workflow,
        "LEGACY_V3_REPLAY_MANIFEST_SHA256",
        workflow._canonical_sha256([]),
    )
    monkeypatch.setattr(
        workflow,
        "LEGACY_V3_REPLAY_PLAN_FILE_SHA256",
        workflow.sha256_file(
            legacy_v3 / workflow.LEGACY_V3_PLAN_OUTPUTS["eas_parallel_replay_plan"]
        ),
    )
    monkeypatch.setattr(
        workflow,
        "LEGACY_V3_REPLAY_PLAN_PAYLOAD_SHA256",
        replay_plan["plan_payload_sha256"],
    )
    monkeypatch.setattr(
        workflow,
        "LEGACY_V3_PLAN_COMPLETION_SHA256",
        workflow.sha256_file(legacy_v3 / workflow.PLAN_COMPLETION_FILENAME),
    )
    return legacy_v3, manifest, extension, replay_plan


def _stub_synthetic_legacy_v4_migration(
    monkeypatch,
    *,
    repo_root: Path,
    work_root: Path,
) -> tuple[Path, dict, dict, dict]:
    legacy_v3, manifest, extension, _v3_replay_plan = _write_synthetic_legacy_v3_plan(
        monkeypatch,
        repo_root=repo_root,
        work_root=work_root,
    )
    legacy_v4 = work_root / workflow.LEGACY_V4_PLAN_DIRNAME
    legacy_v4.mkdir()
    for label, filename in workflow.LEGACY_V2_PLAN_OUTPUTS.items():
        (legacy_v4 / filename).write_bytes((legacy_v3 / filename).read_bytes())
    replay_plan = {
        "schema": workflow.LEGACY_V4_PARALLEL_REPLAY_PLAN_SCHEMA,
        "status": "complete",
        "plan_payload_sha256": "4" * 64,
    }
    workflow._atomic_json(
        legacy_v4 / workflow.LEGACY_V4_PLAN_OUTPUTS["eas_parallel_replay_plan"],
        replay_plan,
    )
    outputs = {
        label: workflow._output_record(legacy_v4 / filename, legacy_v4)
        for label, filename in workflow.LEGACY_V4_PLAN_OUTPUTS.items()
    }
    workflow._atomic_json(
        legacy_v4 / workflow.PLAN_COMPLETION_FILENAME,
        {
            "schema": workflow.LEGACY_V4_PLAN_SCHEMA,
            "status": "complete",
            "outputs": outputs,
        },
    )
    attestation = {
        "schema": "synthetic-v4-attestation",
        "legacy_plan_completion_sha256": workflow.sha256_file(
            legacy_v4 / workflow.PLAN_COMPLETION_FILENAME
        ),
        "legacy_plan_producer_workflow_sha256": workflow.LEGACY_V4_WORKFLOW_SHA256,
        "legacy_output_bindings": {
            label: {
                "size_bytes": int(record["size_bytes"]),
                "sha256": str(record["sha256"]),
            }
            for label, record in outputs.items()
        },
        "legacy_replay_root_manifest": {
            "legacy_replay_root_manifest_sha256": "5" * 64,
            "legacy_file_count": 244,
            "legacy_block_count": 81,
        },
    }
    migration = {
        "attestation": attestation,
        "base_manifest": manifest,
        "extension_plan": extension,
        "legacy_replay_plan": replay_plan,
        "legacy_replay_manifest": attestation["legacy_replay_root_manifest"],
    }
    monkeypatch.setattr(
        workflow,
        "_verify_legacy_v4_plan",
        lambda *_args, **_kwargs: migration,
    )
    return legacy_v4, manifest, extension, migration


def _stub_synthetic_legacy_v5_migration(
    monkeypatch,
    *,
    repo_root: Path,
    work_root: Path,
) -> tuple[Path, dict, dict, dict]:
    legacy_v4, manifest, extension, v4 = _stub_synthetic_legacy_v4_migration(
        monkeypatch,
        repo_root=repo_root,
        work_root=work_root,
    )
    legacy_v5 = work_root / workflow.LEGACY_V5_PLAN_DIRNAME
    legacy_v5.mkdir()
    for label, filename in workflow.LEGACY_V5_PLAN_OUTPUTS.items():
        (legacy_v5 / filename).write_bytes((legacy_v4 / filename).read_bytes())
    outputs = {
        label: workflow._output_record(legacy_v5 / filename, legacy_v5)
        for label, filename in workflow.LEGACY_V5_PLAN_OUTPUTS.items()
    }
    workflow._atomic_json(
        legacy_v5 / workflow.PLAN_COMPLETION_FILENAME,
        {
            "schema": workflow.LEGACY_V5_PLAN_SCHEMA,
            "status": "complete",
            "outputs": outputs,
        },
    )
    attestation = {
        "schema": "synthetic-v5-attestation",
        "legacy_plan_completion_sha256": workflow.sha256_file(
            legacy_v5 / workflow.PLAN_COMPLETION_FILENAME
        ),
        "legacy_plan_producer_workflow_sha256": workflow.LEGACY_V5_WORKFLOW_SHA256,
        "legacy_output_bindings": {
            label: {
                "size_bytes": int(record["size_bytes"]),
                "sha256": str(record["sha256"]),
            }
            for label, record in outputs.items()
        },
        "legacy_replay_root_manifest_sha256": (
            workflow.LEGACY_V4_REPLAY_ROOT_MANIFEST_SHA256
        ),
        "legacy_selection_root_manifest_sha256": (
            workflow.LEGACY_V5_SELECTION_ROOT_MANIFEST_SHA256
        ),
    }
    migration = {
        "attestation": attestation,
        "base_manifest": manifest,
        "extension_plan": extension,
        "legacy_replay_plan": v4["legacy_replay_plan"],
        "selected": pd.DataFrame(),
        "diagnostics": {},
        "selection_completion": {"status": "synthetic"},
    }
    monkeypatch.setattr(
        workflow,
        "_verify_legacy_v5_plan",
        lambda *_args, **_kwargs: migration,
    )
    return legacy_v5, manifest, extension, migration


def _stub_synthetic_legacy_v6_migration(
    monkeypatch,
    *,
    repo_root: Path,
    work_root: Path,
) -> tuple[Path, dict, dict, dict]:
    legacy_v5, manifest, extension, _v5 = _stub_synthetic_legacy_v5_migration(
        monkeypatch,
        repo_root=repo_root,
        work_root=work_root,
    )
    legacy_v6 = work_root / workflow.LEGACY_V6_PLAN_DIRNAME
    legacy_v6.mkdir()
    for label, filename in workflow.LEGACY_V5_PLAN_OUTPUTS.items():
        (legacy_v6 / filename).write_bytes((legacy_v5 / filename).read_bytes())
    (legacy_v6 / workflow.LEGACY_V6_PLAN_OUTPUTS["recombination_map"]).write_bytes(
        (repo_root / workflow.RECOMBINATION_MAP_RELATIVE).read_bytes()
    )
    outputs = {
        label: workflow._output_record(legacy_v6 / filename, legacy_v6)
        for label, filename in workflow.LEGACY_V6_PLAN_OUTPUTS.items()
    }
    workflow._atomic_json(
        legacy_v6 / workflow.PLAN_COMPLETION_FILENAME,
        {
            "schema": workflow.LEGACY_V6_PLAN_SCHEMA,
            "status": "complete",
            "outputs": outputs,
        },
    )
    attestation = {
        "schema": "synthetic-v6-attestation",
        "legacy_plan_completion_sha256": workflow.sha256_file(
            legacy_v6 / workflow.PLAN_COMPLETION_FILENAME
        ),
        "legacy_plan_producer_workflow_sha256": workflow.LEGACY_V6_WORKFLOW_SHA256,
        "legacy_output_bindings": {
            label: {
                "size_bytes": int(record["size_bytes"]),
                "sha256": str(record["sha256"]),
            }
            for label, record in outputs.items()
        },
    }
    migration = {
        "completion": {"status": "synthetic"},
        "attestation": attestation,
        "base_manifest": manifest,
        "extension_plan": extension,
        "selected": pd.DataFrame(),
        "diagnostics": {},
    }
    monkeypatch.setattr(
        workflow,
        "_verify_legacy_v6_plan",
        lambda *_args, **_kwargs: migration,
    )
    monkeypatch.setattr(
        workflow,
        "_verify_legacy_v6_execution_state",
        lambda *_args, **_kwargs: {
            "schema": "synthetic-v6-execution-state",
            "status": "exact_pre_v7_boundary",
        },
    )
    return legacy_v6, manifest, extension, migration


def _synthetic_selected_and_blocks() -> tuple[pd.DataFrame, list[dict]]:
    proposal_count = workflow.conditioned.EAS_PROPOSALS_PER_BLOCK
    selected_blocks = (1, 196)
    block_rows: dict[int, list[tuple[int, int]]] = {
        index: [] for index in selected_blocks
    }
    rows = []
    for accepted_index in range(workflow.EXPECTED_NEUTRAL_PER_DEMOGRAPHY):
        block_index = 1 if accepted_index % 2 == 0 else 196
        proposal_index = 7 if accepted_index in {0, 2} else accepted_index + 20
        block_rows[block_index].append((accepted_index, proposal_index))
        duplicate = accepted_index in {0, 2}
        rows.append(
            {
                "accepted_index": accepted_index,
                "unit_id": (
                    f"{workflow.conditioned.EAS_MODEL_ID}_r{accepted_index:03d}"
                ),
                "cosi2_seed": workflow.conditioned.EAS_SEED_START + accepted_index,
                "block_index": block_index,
                "block_seed": (
                    workflow.conditioned.EAS_PROPOSAL_BLOCK_SEED_START + block_index
                ),
                "proposal_index": proposal_index,
                "global_proposal_index": block_index * proposal_count + proposal_index,
                "birth_generation": 1,
                "log_forward": -2.0,
                "log_reverse": -4.0,
                "log_mutation": 1.0,
                "log_weight": 3.0,
                "source_proposal_block_sha256": f"{block_index + 10:064x}",
                "source_proposal_completion_sha256": f"{block_index + 20:064x}",
                "normalized_global_weight": 0.001,
                "source_particle_occurrence_index": (
                    accepted_index // 2 if duplicate else 0
                ),
                "source_particle_draw_multiplicity": 2 if duplicate else 1,
                "selection_method": (
                    "seeded_iid_categorical_with_replacement_from_global_weighted_bank"
                ),
            }
        )
    blocks = []
    for block_index in range(workflow.conditioned.EAS_PROPOSAL_BLOCK_COUNT):
        block = {
            "block_index": block_index,
            "block_seed": (
                workflow.conditioned.EAS_PROPOSAL_BLOCK_SEED_START + block_index
            ),
            "proposal_count": proposal_count,
            "proposal_block_sha256": f"{block_index + 10:064x}",
            "completion_sha256": f"{block_index + 20:064x}",
        }
        if block_index in selected_blocks:
            outcomes = np.full(proposal_count, "invalid_loss", dtype="U12")
            birth_generations = np.full(proposal_count, -1, dtype=np.int64)
            log_forward = np.zeros(proposal_count, dtype=float)
            log_reverse = np.zeros(proposal_count, dtype=float)
            log_mutation = np.zeros(proposal_count, dtype=float)
            log_weights = np.full(proposal_count, -np.inf, dtype=float)
            for _accepted_index, proposal_index in block_rows[block_index]:
                outcomes[proposal_index] = "valid"
                birth_generations[proposal_index] = 1
                log_forward[proposal_index] = -2.0
                log_reverse[proposal_index] = -4.0
                log_mutation[proposal_index] = 1.0
                log_weights[proposal_index] = 3.0
            block.update(
                {
                    "outcomes": outcomes,
                    "birth_generations": birth_generations,
                    "log_forward": log_forward,
                    "log_reverse": log_reverse,
                    "log_mutation": log_mutation,
                    "log_weights": log_weights,
                }
            )
        blocks.append(block)
    return pd.DataFrame(rows), blocks


def _write_synthetic_legacy_replay_block(
    *,
    bundle_dir: Path,
    selected_group: pd.DataFrame,
    block: dict,
    base_manifest: dict,
    extension_plan: dict,
    legacy_replay_plan: dict,
) -> dict:
    contract = workflow._legacy_v3_parallel_replay_block_contract(
        selected_group=selected_group,
        block=block,
        legacy_parallel_replay_plan=legacy_replay_plan,
        base_manifest=base_manifest,
        extension_plan=extension_plan,
    )
    paths = {unit_id: (14_704, 1) for unit_id in selected_group["unit_id"].astype(str)}
    bundle_dir.mkdir(parents=True)
    workflow._atomic_frame(bundle_dir / "selected_rows.tsv", selected_group)
    np.savez_compressed(
        bundle_dir / "paths.npz",
        **workflow._pack_selected_paths(selected_group, paths),
    )
    outputs = {
        "selected_rows": workflow._output_record(
            bundle_dir / "selected_rows.tsv",
            bundle_dir,
            rows=len(selected_group),
        ),
        "paths": workflow._output_record(bundle_dir / "paths.npz", bundle_dir),
    }
    workflow._atomic_json(
        bundle_dir / "completion.json",
        {
            "schema": workflow.LEGACY_V3_PARALLEL_REPLAY_BLOCK_SCHEMA,
            "status": "complete",
            "contract": contract,
            "contract_sha256": workflow._canonical_sha256(contract),
            "outputs": outputs,
        },
    )
    files = sorted(bundle_dir.iterdir(), key=lambda path: path.name)
    return {
        "block_dir": bundle_dir.name,
        "files": {
            path.name: {
                "size_bytes": path.stat().st_size,
                "sha256": workflow.sha256_file(path),
            }
            for path in files
        },
    }


def _units() -> pd.DataFrame:
    rows = []
    for demography, settings in workflow.DEMOGRAPHY_SETTINGS.items():
        for index in range(workflow.EXPECTED_NEUTRAL_PER_DEMOGRAPHY):
            count = index % 25
            eas = demography == "EAS"
            birth = 2_000 + index if eas else workflow.conditioned.HAN_BIRTH_GENERATION
            rows.append(
                {
                    "unit_id": f"{settings['neutral_model_id']}_r{index:03d}",
                    "demography": demography,
                    "demography_id": settings["demography_id"],
                    "model_id": settings["neutral_model_id"],
                    "simulation_class": "neutral",
                    "selection_coefficient": 0.0,
                    "seed": (
                        workflow.conditioned.EAS_SEED_START + index
                        if eas
                        else workflow.conditioned.HAN_NEW_SEED_START + index
                    ),
                    "generation_time_years": settings["generation_time_years"],
                    "present_ne": settings["present_ne"],
                    "final_population_af": (
                        workflow.conditioned.EAS_SELECTED_PRESENT_AF
                        if eas
                        else 0.001 + index / 1_000
                    ),
                    "sample_alt_count": count,
                    "sample_af": count / analysis.PANEL_HAPLOTYPES,
                    "ms_path": f"/read-only/{demography}/{index}.ms",
                    "ms_sha256": f"{index + 1000 * (not eas):064x}"[-64:],
                    "source_completion_path": f"/read-only/{demography}/{index}.json",
                    "source_completion_sha256": f"{index + 5000:064x}"[-64:],
                    "conditioning": "synthetic_conditioned_null",
                    "birth_generation": birth,
                    "allele_age_years": birth * settings["generation_time_years"],
                    "eas_selected_log_importance_weight": -100.0 + index
                    if eas
                    else float("nan"),
                    "eas_global_effective_sample_size": 20_000.0
                    if eas
                    else float("nan"),
                    "eas_global_weight_sum_over_max": 5_000.0 if eas else float("nan"),
                    "eas_global_max_normalized_weight": 0.001 if eas else float("nan"),
                    "eas_time_reversal_is_approximation": eas,
                    "eas_source_particle_id": f"particle_{index}" if eas else "",
                    "eas_draw_occurrence": index if eas else -1,
                    "eas_particle_multiplicity": 1 if eas else -1,
                    "han_candidate_id": f"han_candidate_{index:04d}" if not eas else "",
                    "han_candidate_order": index if not eas else -1,
                    "han_candidate_source_kind": "new_screen_then_replay"
                    if not eas
                    else "",
                    "han_generation_645_chb_af": 0.012 if not eas else float("nan"),
                }
            )
        selected_count = 40 if demography == "EAS" else 32
        selected_birth = (
            1_475 if demography == "EAS" else workflow.conditioned.HAN_BIRTH_GENERATION
        )
        rows.append(
            {
                "unit_id": settings["selected_unit_id"],
                "demography": demography,
                "demography_id": settings["demography_id"],
                "model_id": settings["selected_cell_id"],
                "simulation_class": "selected",
                "selection_coefficient": 0.01,
                "seed": settings["selected_seed"],
                "generation_time_years": settings["generation_time_years"],
                "present_ne": settings["present_ne"],
                "final_population_af": 0.1952,
                "sample_alt_count": selected_count,
                "sample_af": selected_count / analysis.PANEL_HAPLOTYPES,
                "ms_path": f"/read-only/{demography}/selected.ms",
                "ms_sha256": "a" * 64,
                "source_completion_path": f"/read-only/{demography}/selected.json",
                "source_completion_sha256": "b" * 64,
                "conditioning": "selected_endpoint",
                "birth_generation": selected_birth,
                "allele_age_years": selected_birth * settings["generation_time_years"],
                "eas_selected_log_importance_weight": float("nan"),
                "eas_global_effective_sample_size": float("nan"),
                "eas_global_weight_sum_over_max": float("nan"),
                "eas_global_max_normalized_weight": float("nan"),
                "eas_time_reversal_is_approximation": False,
                "eas_source_particle_id": "",
                "eas_draw_occurrence": -1,
                "eas_particle_multiplicity": -1,
                "han_candidate_id": "",
                "han_candidate_order": -1,
                "han_candidate_source_kind": "",
                "han_generation_645_chb_af": float("nan"),
            }
        )
    return workflow.validate_unit_inventory(pd.DataFrame(rows))


def _scores(units: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for unit in units.to_dict(orient="records"):
        for statistic_index, statistic in enumerate(analysis.SUPPORTED_STATISTICS):
            for window_index, window in enumerate(analysis.SUPPORTED_WINDOWS):
                for threshold_index, threshold in enumerate(
                    analysis.TMRCA_THRESHOLDS_YEARS
                ):
                    rows.append(
                        {
                            **{
                                column: unit[column]
                                for column in analysis.UNIT_METADATA_COLUMNS
                            },
                            "statistic": statistic,
                            "window": window,
                            "threshold_years": threshold,
                            "score": min(
                                0.99,
                                0.04
                                + threshold_index * 0.08
                                + statistic_index * 0.01
                                + window_index * 0.005
                                + 0.1 * (unit["simulation_class"] == "selected"),
                            ),
                        }
                    )
    return analysis.validate_scores(pd.DataFrame(rows))


def _decodes(units: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "unit_id": units["unit_id"],
            "status": "complete",
            "decode_completion_sha256": [
                f"{index + 8000:064x}" for index in range(len(units))
            ],
            "decode_full_sample_alt_count": units["sample_alt_count"].astype(int),
            "decode_pair_panel_alt_count": (
                units["sample_alt_count"].astype(int) // 2
            ).clip(upper=100),
            "decode_focal_count_status": "unique_exact_focal_site",
        }
    )


def _fake_plot(_data, output_dir, *, stem=workflow.PROFILE_FIGURE_STEM):
    destination = Path(output_dir)
    png = destination / f"{stem}.png"
    pdf = destination / f"{stem}.pdf"
    png.write_bytes(b"\x89PNG\r\n\x1a\nsynthetic")
    pdf.write_bytes(b"%PDF-1.4\nsynthetic\n%%EOF\n")
    return {"png": png, "pdf": pdf}


def _write_minimal_v6_han_plan(repo_root: Path, plan_root: Path) -> Path:
    plan_root.mkdir(parents=True)
    config = plan_root / workflow.PLAN_OUTPUTS["han_screen_config"]
    config.write_text(
        workflow.conditioned.render_han_parameter_text(repo_root, screen=True),
        encoding="utf-8",
    )
    recombination_map = plan_root / workflow.PLAN_OUTPUTS["recombination_map"]
    recombination_map.write_bytes(
        (repo_root / workflow.RECOMBINATION_MAP_RELATIVE).read_bytes()
    )
    outputs = {
        "han_screen_config": workflow._output_record(config, plan_root),
        "recombination_map": workflow._output_record(recombination_map, plan_root),
    }
    workflow._atomic_json(
        plan_root / workflow.PLAN_COMPLETION_FILENAME,
        {
            "schema": workflow.PLAN_SCHEMA,
            "status": "complete",
            "outputs": outputs,
        },
    )
    return plan_root


def _write_minimal_v7_han_replay_plan(repo_root: Path, plan_root: Path) -> Path:
    plan_root.mkdir(parents=True)
    config = plan_root / workflow.PLAN_OUTPUTS["han_replay_config"]
    config.write_text(
        workflow.conditioned.render_han_parameter_text(repo_root, screen=False),
        encoding="utf-8",
    )
    recombination_map = plan_root / workflow.PLAN_OUTPUTS["recombination_map"]
    recombination_map.write_bytes(
        (repo_root / workflow.RECOMBINATION_MAP_RELATIVE).read_bytes()
    )
    outputs = {
        "han_replay_config": workflow._output_record(config, plan_root),
        "recombination_map": workflow._output_record(recombination_map, plan_root),
    }
    workflow._atomic_json(
        plan_root / workflow.PLAN_COMPLETION_FILENAME,
        {
            "schema": workflow.PLAN_SCHEMA,
            "status": "complete",
            "outputs": outputs,
        },
    )
    return plan_root


def test_inventory_and_conditioning_contrasts():
    units = _units()
    assert units.groupby(["demography", "simulation_class"]).size().to_dict() == {
        ("EAS", "neutral"): 100,
        ("EAS", "selected"): 1,
        ("Han", "neutral"): 100,
        ("Han", "selected"): 1,
    }
    contrasts = workflow.conditioning_contrasts(units).set_index("demography")
    assert contrasts.loc["EAS", "neutral_extreme_count"] == 0
    assert contrasts.loc["EAS", "plus_one_simulation_p"] == pytest.approx(1 / 101)
    assert contrasts.loc["Han", "tail"] == "upper"


def test_plan_bundle_is_atomic_and_reusable(monkeypatch, tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    destination = tmp_path / workflow.PLAN_DIRNAME
    legacy, _manifest, _extension, _migration = _stub_synthetic_legacy_v6_migration(
        monkeypatch,
        repo_root=repo_root,
        work_root=tmp_path,
    )
    created = workflow.generate_plan(
        repo_root,
        output_dir=destination,
        proposal_dir=tmp_path / "proposal_blocks",
    )
    assert created["cache_hit"] is False
    assert {path.name for path in destination.iterdir()} == {
        workflow.PLAN_COMPLETION_FILENAME,
        *workflow.PLAN_OUTPUTS.values(),
    }
    ledger = pd.read_csv(
        destination / workflow.PLAN_OUTPUTS["eas_proposal_blocks"], sep="\t"
    )
    assert len(ledger) == 200
    assert ledger.loc[0, "source_contract"] == "frozen_v1_pilot_hash_attested"
    assert ledger.loc[20, "source_contract"] == "post_pilot_v2_extension_current_source"
    completion = workflow._read_json(
        destination / workflow.PLAN_COMPLETION_FILENAME, "test plan completion"
    )
    bank_contract = completion["contract"]["eas_proposal_bank_contract"]
    assert bank_contract["target_block_count"] == 200
    assert bank_contract["total_proposals"] == 10_000_000
    assert bank_contract["adaptive_stopping_allowed"] is False
    assert (
        workflow.sha256_file(destination / workflow.PLAN_OUTPUTS["recombination_map"])
        == workflow.RECOMBINATION_MAP_SHA256
    )
    assert completion["contract"]["legacy_v6_plan_attestation"][
        "legacy_plan_completion_sha256"
    ] == workflow.sha256_file(legacy / workflow.PLAN_COMPLETION_FILENAME)
    assert completion["contract"]["v7_revision_scope"]["loader_projection_schema"] == (
        workflow.HAN_LOADER_TRAJECTORY_SCHEMA
    )
    cached = workflow.generate_plan(
        repo_root,
        output_dir=destination,
        proposal_dir=tmp_path / "proposal_blocks",
    )
    assert cached["cache_hit"] is True

    (tmp_path / workflow.EAS_SELECTION_DIRNAME).mkdir()
    assert (
        workflow.verify_plan(
            repo_root,
            plan_dir=destination,
            proposal_dir=tmp_path / "proposal_blocks",
        )["cache_hit"]
        is True
    )

    inherited = destination / workflow.PLAN_OUTPUTS["han_screen_config"]
    inherited.write_bytes(inherited.read_bytes().replace(b"\n", b"\r\n"))
    completion = workflow._read_json(
        destination / workflow.PLAN_COMPLETION_FILENAME, "test plan completion"
    )
    completion["outputs"]["han_screen_config"] = workflow._output_record(
        inherited, destination
    )
    workflow._atomic_json(destination / workflow.PLAN_COMPLETION_FILENAME, completion)
    with pytest.raises(ValueError, match="inherited v6 bytes changed"):
        workflow.verify_plan(
            repo_root,
            plan_dir=destination,
            proposal_dir=tmp_path / "proposal_blocks",
        )


def test_v7_plan_rejects_missing_or_tampered_recombination_map(monkeypatch, tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    _stub_synthetic_legacy_v6_migration(
        monkeypatch,
        repo_root=repo_root,
        work_root=tmp_path,
    )
    destination = tmp_path / workflow.PLAN_DIRNAME
    workflow.generate_plan(
        repo_root,
        output_dir=destination,
        proposal_dir=tmp_path / "proposal_blocks",
    )
    map_path = destination / workflow.PLAN_OUTPUTS["recombination_map"]
    original = map_path.read_bytes()
    map_path.unlink()
    with pytest.raises(ValueError, match="inventory mismatch"):
        workflow.verify_plan(
            repo_root,
            plan_dir=destination,
            proposal_dir=tmp_path / "proposal_blocks",
        )
    map_path.write_bytes(original + b"tamper")
    with pytest.raises(ValueError, match="checksum failed"):
        workflow.verify_plan(
            repo_root,
            plan_dir=destination,
            proposal_dir=tmp_path / "proposal_blocks",
        )


def test_legacy_v4_absence_gate_allows_only_v5_selection(tmp_path):
    current = tmp_path / workflow.EAS_SELECTION_DIRNAME
    current.mkdir()
    state = workflow._verify_legacy_v4_unpublished_state(tmp_path)
    assert state["current_v5_selection_allowed"] is True
    assert workflow._resolve_paths(tmp_path, tmp_path, None)["eas_selection"] == current

    legacy = tmp_path / workflow.LEGACY_V4_SELECTION_DIRNAME
    legacy.mkdir()
    with pytest.raises(ValueError, match="selection or replay staging"):
        workflow._verify_legacy_v4_unpublished_state(tmp_path)
    legacy.rmdir()

    legacy_stage = tmp_path / f".{workflow.LEGACY_V4_SELECTION_DIRNAME}.staging.test"
    legacy_stage.mkdir()
    with pytest.raises(ValueError, match="selection or replay staging"):
        workflow._verify_legacy_v4_unpublished_state(tmp_path)
    legacy_stage.rmdir()

    replay_stage = (
        tmp_path / f".{workflow.LEGACY_V4_PARALLEL_REPLAY_DIRNAME}.staging.test"
    )
    replay_stage.mkdir()
    with pytest.raises(ValueError, match="selection or replay staging"):
        workflow._verify_legacy_v4_unpublished_state(tmp_path)


def test_eas_extension_uses_process_pool_and_fixed_indexes(
    monkeypatch, tmp_path, capsys
):
    used = {"entered": 0}

    class ImmediateProcessPool:
        def __init__(self, *, max_workers):
            assert max_workers == 2

        def __enter__(self):
            used["entered"] += 1
            return self

        def __exit__(self, *_args):
            return False

        def submit(self, function, *args):
            future = Future()
            future.set_result(function(*args))
            return future

    monkeypatch.setattr(workflow, "ProcessPoolExecutor", ImmediateProcessPool)
    monkeypatch.setattr(
        workflow.conditioned,
        "validate_eas_base_bank_manifest",
        lambda *_args, **_kwargs: {"status": "valid"},
    )
    monkeypatch.setattr(
        workflow.conditioned,
        "validate_eas_extension_plan",
        lambda *_args, **_kwargs: {"status": "valid"},
    )
    monkeypatch.setattr(
        workflow,
        "_run_eas_extension_block_process",
        lambda _root, _proposal, index, _manifest, _plan: {
            "block_index": index,
            "status": "complete",
            "cache_hit": False,
            "worker_process_id": 123,
            "elapsed_seconds": 1.0,
            "valid_proposals": 10,
            "proposal_block_sha256": f"{index:064x}",
            "completion_sha256": f"{index + 1:064x}",
        },
    )
    result = workflow.run_eas_proposal_blocks(
        tmp_path,
        proposal_dir=tmp_path / "blocks",
        base_manifest=_synthetic_base_manifest(),
        extension_plan=_synthetic_extension_plan(),
        max_workers=2,
        block_indexes=[20, 21],
    )
    assert used["entered"] == 1
    assert result["block_index"].tolist() == [20, 21]
    progress = capsys.readouterr().err
    assert "eas_extension_process_pool_start" in progress
    assert progress.count("eas_extension_block_complete") == 2
    assert "eas_extension_process_pool_finished" in progress
    with pytest.raises(ValueError, match="20..199"):
        workflow.run_eas_proposal_blocks(
            tmp_path,
            proposal_dir=tmp_path / "blocks",
            base_manifest=_synthetic_base_manifest(),
            extension_plan=_synthetic_extension_plan(),
            block_indexes=[19],
        )


def test_eas_extension_empty_request_does_not_start_pool(monkeypatch, tmp_path):
    monkeypatch.setattr(
        workflow.conditioned,
        "validate_eas_base_bank_manifest",
        lambda *_args, **_kwargs: {"status": "valid"},
    )
    result = workflow.run_eas_proposal_blocks(
        tmp_path,
        proposal_dir=tmp_path / "blocks",
        base_manifest=_synthetic_base_manifest(),
        extension_plan=_synthetic_extension_plan(),
        block_indexes=[],
    )
    assert result.empty


def test_eas_loader_rejects_extra_inventory(monkeypatch, tmp_path):
    proposal_root = tmp_path / "blocks"
    for index in range(workflow.conditioned.EAS_PROPOSAL_BLOCK_COUNT):
        (proposal_root / f"block_{index:02d}").mkdir(parents=True)
    monkeypatch.setattr(
        workflow.conditioned,
        "validate_eas_base_bank_manifest",
        lambda *_args, **_kwargs: {"status": "valid"},
    )
    monkeypatch.setattr(
        workflow.conditioned,
        "load_eas_extended_proposal_bank",
        lambda *_args, **_kwargs: [{"block_index": index} for index in range(200)],
    )
    blocks = workflow.load_eas_proposal_blocks(
        tmp_path,
        proposal_dir=proposal_root,
        base_manifest=_synthetic_base_manifest(),
        extension_plan=_synthetic_extension_plan(),
    )
    assert [block["block_index"] for block in blocks] == list(range(200))
    (proposal_root / "unexpected.txt").write_text("tamper\n", encoding="utf-8")
    with pytest.raises(ValueError, match="fixed-200"):
        workflow.load_eas_proposal_blocks(
            tmp_path,
            proposal_dir=proposal_root,
            base_manifest=_synthetic_base_manifest(),
            extension_plan=_synthetic_extension_plan(),
        )


def test_selected_tsv_reader_round_trips_production_magnitude_floats(tmp_path):
    source_values = np.asarray(
        [-44163.69483891431, -44159.01385987655],
        dtype=np.float64,
    )
    frame = pd.DataFrame(
        {
            "unit_id": ["eas_0", "eas_1"],
            "selection_method": ["synthetic", "synthetic"],
            "source_proposal_block_sha256": ["a" * 64, "b" * 64],
            "source_proposal_completion_sha256": ["c" * 64, "d" * 64],
            "log_forward": source_values,
            "log_reverse": [-44160.0, -44155.0],
            "log_mutation": [1.25, 1.5],
            "log_weight": [-2.444838914307649, -2.513859876547474],
        }
    )
    path = tmp_path / "selected_rows.tsv"
    workflow._atomic_frame(path, frame)

    observed = workflow._read_eas_selected_frame(path)
    observed_values = observed["log_forward"].to_numpy(dtype=np.float64)
    assert (
        observed_values.view(np.uint64).tolist()
        == source_values.view(np.uint64).tolist()
    )
    assert [value.hex() for value in observed_values] == [
        value.hex() for value in source_values
    ]

    default_values = pd.read_csv(path, sep="\t")["log_forward"].to_numpy(float)
    assert (
        default_values.view(np.uint64).tolist()
        != source_values.view(np.uint64).tolist()
    )


def test_selected_particle_index_adapter_calls_real_core_writer(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    selected, _blocks = _synthetic_selected_and_blocks()
    selected = selected.iloc[[0]].copy()
    counts = (workflow.conditioned.EAS_SELECTED_PRESENT_COUNT, 1)
    components = workflow.conditioned._reverse_path_log_components(
        workflow.conditioned.eas_integer_ne_by_generation(repo_root), counts
    )
    for label, value in components.items():
        selected.loc[selected.index[0], label] = value
    selected_index = selected.set_index("unit_id", drop=False)
    unit_id = str(selected.iloc[0]["unit_id"])
    row = workflow._selected_particle_row(selected_index, unit_id)
    assert row["unit_id"] == unit_id

    diagnostics = {"test_contract": "real_core_writer_unit_id_adapter"}
    bundle = tmp_path / unit_id
    written = workflow.conditioned.write_eas_selected_path_bundle(
        repo_root,
        bundle,
        selected_row=row,
        derived_counts_ascending=counts,
        global_diagnostics=diagnostics,
    )
    assert written["unit_id"] == unit_id
    assert {path.name for path in bundle.iterdir()} == {
        "bridge_counts.tsv",
        "trajectory.tsv",
        "parameters.par",
        "selection.json",
        "completion.json",
    }
    verified = workflow.conditioned.verify_eas_selected_path_bundle(
        repo_root,
        bundle,
        expected_selected_row=row,
        expected_counts=counts,
        expected_global_diagnostics=diagnostics,
    )
    assert verified["completion_sha256"] == written["completion_sha256"]

    dropped = selected.set_index("unit_id")
    with pytest.raises(ValueError, match="drop=True"):
        workflow._selected_particle_row(dropped, unit_id)


def test_legacy_v3_replay_block_migrates_byte_exact_and_atomically(
    monkeypatch, tmp_path
):
    selected, blocks = _synthetic_selected_and_blocks()
    selected_group = workflow._selected_block_frame(selected, 1)
    block = blocks[1]
    manifest = _synthetic_base_manifest()
    extension = _synthetic_extension_plan()
    legacy_plan = {"plan_payload_sha256": workflow.LEGACY_V3_REPLAY_PLAN_PAYLOAD_SHA256}
    legacy_bundle = tmp_path / workflow.LEGACY_V3_PARALLEL_REPLAY_DIRNAME / "block_01"
    binding = _write_synthetic_legacy_replay_block(
        bundle_dir=legacy_bundle,
        selected_group=selected_group,
        block=block,
        base_manifest=manifest,
        extension_plan=extension,
        legacy_replay_plan=legacy_plan,
    )
    monkeypatch.setattr(workflow, "LEGACY_V3_REPLAY_BLOCK_INDEXES", (1,))
    v4_plan = {
        "plan_payload_sha256": "4" * 64,
        "legacy_v3_plan_attestation": {
            "legacy_replay_manifest": {
                "legacy_block_bindings": [binding],
            }
        },
    }
    monkeypatch.setattr(
        workflow.conditioned,
        "verify_eas_proposal_block",
        lambda *_args, **_kwargs: block,
    )
    monkeypatch.setattr(
        workflow.conditioned,
        "eas_integer_ne_by_generation",
        lambda _root: np.asarray([10_000, 10_000, 10_000], dtype=np.int64),
    )
    monkeypatch.setattr(
        workflow.conditioned,
        "eas_selected_path_frames",
        lambda _sizes, _path: (pd.DataFrame(), pd.DataFrame()),
    )
    monkeypatch.setattr(
        workflow.conditioned,
        "_reverse_path_log_components",
        lambda _sizes, _path: {
            "log_forward": -2.0,
            "log_reverse": -4.0,
            "log_mutation": 1.0,
            "log_weight": 3.0,
        },
    )
    destination = tmp_path / workflow.EAS_PARALLEL_REPLAY_DIRNAME / "block_01"
    call_args = (
        str(tmp_path),
        str(tmp_path / workflow.EAS_PROPOSAL_DIRNAME),
        str(destination),
        [
            workflow._portable_record(record)
            for record in selected_group.to_dict(orient="records")
        ],
        manifest,
        extension,
        v4_plan,
        "5" * 64,
        str(legacy_bundle),
        legacy_plan,
    )

    original_copyfile = workflow.shutil.copyfile
    copy_calls = 0

    def fail_second_copy(source, target):
        nonlocal copy_calls
        copy_calls += 1
        if copy_calls == 2:
            raise OSError("synthetic interrupted migration")
        return original_copyfile(source, target)

    with monkeypatch.context() as context:
        context.setattr(workflow.shutil, "copyfile", fail_second_copy)
        with pytest.raises(OSError, match="interrupted migration"):
            workflow._write_eas_parallel_replay_block(*call_args)
    assert not destination.exists()
    assert not list(destination.parent.glob(".block_01.legacy-v3-staging.*"))

    migrated = workflow._write_eas_parallel_replay_block(*call_args)
    assert migrated["artifact_provenance_kind"] == "exact_legacy_v3_migration"
    assert migrated["cache_hit"] is True
    assert migrated["legacy_v3_migration_cache_hit"] is True
    assert {path.name for path in destination.iterdir()} == set(binding["files"])
    for filename, record in binding["files"].items():
        assert workflow.sha256_file(destination / filename) == record["sha256"]
        assert workflow.sha256_file(legacy_bundle / filename) == record["sha256"]

    tampered = legacy_bundle / "selected_rows.tsv"
    original = tampered.read_bytes()
    tampered.write_bytes(original + b"\n")
    with pytest.raises(ValueError, match="frozen manifest"):
        workflow._write_eas_parallel_replay_block(
            *(call_args[:2] + (str(destination.parent / "block_02"),) + call_args[3:])
        )
    tampered.write_bytes(original)


def test_v5_entire_replay_root_migration_is_atomic_and_byte_identical(
    monkeypatch, tmp_path
):
    def forbidden_replay(*_args, **_kwargs):
        raise AssertionError("v5 migration must not replay or resample")

    monkeypatch.setattr(
        workflow.conditioned,
        "replay_eas_global_selection",
        forbidden_replay,
    )
    source = tmp_path / workflow.LEGACY_V4_PARALLEL_REPLAY_DIRNAME
    (source / "block_01").mkdir(parents=True)
    (source / "block_02").mkdir()
    files = {
        workflow.EAS_PARALLEL_REPLAY_COMPLETION_FILENAME: b"root-completion\n",
        "block_01/completion.json": b"block-one\n",
        "block_01/selected_rows.tsv": b"selected-one\n",
        "block_01/paths.npz": b"paths-one\n",
        "block_02/completion.json": b"block-two\n",
        "block_02/selected_rows.tsv": b"selected-two\n",
        "block_02/paths.npz": b"paths-two\n",
    }
    for relative, content in files.items():
        (source / relative).write_bytes(content)

    def synthetic_manifest(directory):
        root = Path(directory)
        bindings = [
            {
                "path": path.relative_to(root).as_posix(),
                "size_bytes": path.stat().st_size,
                "sha256": workflow.sha256_file(path),
            }
            for path in sorted(
                root.rglob("*"), key=lambda path: path.relative_to(root).as_posix()
            )
            if path.is_file()
        ]
        return {
            "legacy_file_bindings": bindings,
            "legacy_replay_root_manifest_sha256": workflow._canonical_sha256(bindings),
            "legacy_file_count": len(bindings),
            "legacy_block_indexes": [1, 2],
        }

    monkeypatch.setattr(workflow, "_legacy_v4_replay_root_manifest", synthetic_manifest)
    selected, _blocks = _synthetic_selected_and_blocks()
    paths = {unit_id: (14_704, 1) for unit_id in selected["unit_id"].astype(str)}

    def synthetic_semantics(*_args, replay_dir, **_kwargs):
        assert synthetic_manifest(replay_dir)["legacy_file_count"] == len(files)
        return {
            "status": "complete",
            "cache_hit": True,
            "completion_sha256": workflow.sha256_file(
                Path(replay_dir) / workflow.EAS_PARALLEL_REPLAY_COMPLETION_FILENAME
            ),
            "unique_selected_block_count": 2,
            "block_bindings": [],
            "paths": paths,
        }

    monkeypatch.setattr(
        workflow,
        "_verify_exact_legacy_v4_replay_root_semantics",
        synthetic_semantics,
    )
    source_manifest = synthetic_manifest(source)
    replay_plan = {
        "legacy_v4_plan_attestation": {
            "legacy_replay_root_manifest": source_manifest,
        }
    }
    destination = tmp_path / workflow.EAS_PARALLEL_REPLAY_DIRNAME
    kwargs = {
        "proposal_dir": tmp_path / workflow.EAS_PROPOSAL_DIRNAME,
        "source_replay_dir": source,
        "destination_replay_dir": destination,
        "selected": selected,
        "base_manifest": _synthetic_base_manifest(),
        "extension_plan": _synthetic_extension_plan(),
        "current_replay_plan": replay_plan,
        "legacy_v4_parallel_replay_plan": {},
        "legacy_v3_parallel_replay_plan": {},
    }
    original_copyfile = workflow.shutil.copyfile
    copy_count = 0

    def interrupted_copy(source_path, destination_path):
        nonlocal copy_count
        copy_count += 1
        if copy_count == 3:
            raise OSError("synthetic whole-root copy interruption")
        return original_copyfile(source_path, destination_path)

    with monkeypatch.context() as context:
        context.setattr(workflow.shutil, "copyfile", interrupted_copy)
        with pytest.raises(OSError, match="whole-root copy interruption"):
            workflow._migrate_exact_legacy_v4_replay_root(tmp_path, **kwargs)
    assert not destination.exists()
    assert not list(tmp_path.glob(f".{destination.name}.legacy-v4-staging.*"))

    migrated = workflow._migrate_exact_legacy_v4_replay_root(tmp_path, **kwargs)
    assert migrated["cache_hit"] is False
    assert {
        path.relative_to(destination).as_posix(): workflow.sha256_file(path)
        for path in destination.rglob("*")
        if path.is_file()
    } == {
        path.relative_to(source).as_posix(): workflow.sha256_file(path)
        for path in source.rglob("*")
        if path.is_file()
    }
    cached = workflow._migrate_exact_legacy_v4_replay_root(tmp_path, **kwargs)
    assert cached["cache_hit"] is True

    (destination / "block_02" / "paths.npz").write_bytes(b"tampered\n")
    with pytest.raises(ValueError, match="destination bytes changed"):
        workflow._migrate_exact_legacy_v4_replay_root(tmp_path, **kwargs)


def test_v5_replay_verifier_accepts_legacy_completion_only_via_hard_attestation(
    monkeypatch, tmp_path
):
    manifest = {
        "legacy_replay_root_manifest_sha256": (
            workflow.LEGACY_V4_REPLAY_ROOT_MANIFEST_SHA256
        )
    }
    attestation = {
        "legacy_plan_completion_sha256": workflow.LEGACY_V4_PLAN_COMPLETION_SHA256,
        "legacy_plan_producer_workflow_sha256": workflow.LEGACY_V4_WORKFLOW_SHA256,
        "legacy_replay_root_manifest": manifest,
    }
    replay_plan = workflow._build_eas_parallel_replay_plan(attestation)
    monkeypatch.setattr(
        workflow,
        "_legacy_v4_replay_root_manifest",
        lambda _root: manifest,
    )
    sentinel = {
        "status": "complete",
        "completion_interpretation": "hard_bound_legacy_v4_only",
    }
    monkeypatch.setattr(
        workflow,
        "_verify_exact_legacy_v4_replay_root_semantics",
        lambda *_args, **_kwargs: sentinel,
    )
    selected, _blocks = _synthetic_selected_and_blocks()
    verified = workflow.verify_eas_parallel_replay_completion(
        tmp_path,
        proposal_dir=tmp_path / workflow.EAS_PROPOSAL_DIRNAME,
        replay_dir=tmp_path / workflow.EAS_PARALLEL_REPLAY_DIRNAME,
        selected=selected,
        base_manifest=_synthetic_base_manifest(),
        extension_plan=_synthetic_extension_plan(),
        parallel_replay_plan=replay_plan,
        legacy_v3_parallel_replay_plan={},
        legacy_v4_parallel_replay_plan={
            "plan_payload_sha256": workflow.LEGACY_V4_REPLAY_PLAN_PAYLOAD_SHA256
        },
        plan_completion_sha256="f" * 64,
    )
    assert verified is sentinel

    tampered = dict(replay_plan)
    tampered["legacy_v4_plan_attestation"] = {
        **attestation,
        "legacy_plan_producer_workflow_sha256": "0" * 64,
    }
    with pytest.raises(ValueError, match="not the exact attested v4 predecessor"):
        workflow.verify_eas_parallel_replay_completion(
            tmp_path,
            proposal_dir=tmp_path / workflow.EAS_PROPOSAL_DIRNAME,
            replay_dir=tmp_path / workflow.EAS_PARALLEL_REPLAY_DIRNAME,
            selected=selected,
            base_manifest=_synthetic_base_manifest(),
            extension_plan=_synthetic_extension_plan(),
            parallel_replay_plan=tampered,
            legacy_v3_parallel_replay_plan={},
            legacy_v4_parallel_replay_plan={
                "plan_payload_sha256": (workflow.LEGACY_V4_REPLAY_PLAN_PAYLOAD_SHA256)
            },
            plan_completion_sha256="f" * 64,
        )


def test_parallel_eas_replay_is_block_grouped_restartable_and_ordered(
    monkeypatch, tmp_path, capsys
):
    selected, blocks = _synthetic_selected_and_blocks()
    indexed_blocks = {int(block["block_index"]): block for block in blocks}
    manifest = _synthetic_base_manifest()
    extension = _synthetic_extension_plan()
    legacy_replay_plan = {
        "plan_payload_sha256": workflow.LEGACY_V3_REPLAY_PLAN_PAYLOAD_SHA256
    }
    legacy_binding = _write_synthetic_legacy_replay_block(
        bundle_dir=(tmp_path / workflow.LEGACY_V3_PARALLEL_REPLAY_DIRNAME / "block_01"),
        selected_group=workflow._selected_block_frame(selected, 1),
        block=blocks[1],
        base_manifest=manifest,
        extension_plan=extension,
        legacy_replay_plan=legacy_replay_plan,
    )
    monkeypatch.setattr(workflow, "LEGACY_V3_REPLAY_BLOCK_INDEXES", (1,))
    plan_root = tmp_path / workflow.PLAN_DIRNAME
    plan_root.mkdir()
    replay_plan = {
        "plan_payload_sha256": "4" * 64,
        "legacy_v3_plan_attestation": {
            "legacy_replay_manifest": {
                "legacy_block_bindings": [legacy_binding],
            }
        },
    }
    workflow._atomic_json(
        plan_root / workflow.PLAN_OUTPUTS["eas_parallel_replay_plan"], replay_plan
    )
    legacy_plan_root = tmp_path / workflow.LEGACY_V3_PLAN_DIRNAME
    legacy_plan_root.mkdir()
    workflow._atomic_json(
        legacy_plan_root / workflow.LEGACY_V3_PLAN_OUTPUTS["eas_parallel_replay_plan"],
        legacy_replay_plan,
    )
    workflow._atomic_json(
        plan_root / workflow.PLAN_COMPLETION_FILENAME, {"status": "synthetic"}
    )
    monkeypatch.setattr(workflow, "verify_plan", lambda *_args, **_kwargs: {})
    monkeypatch.setattr(
        workflow.conditioned,
        "verify_eas_proposal_block",
        lambda _root, _bundle, *, expected_block_index, **_kwargs: indexed_blocks[
            int(expected_block_index)
        ],
    )
    monkeypatch.setattr(
        workflow.conditioned,
        "eas_integer_ne_by_generation",
        lambda _root: np.asarray([10_000, 10_000, 10_000], dtype=np.int64),
    )
    monkeypatch.setattr(
        workflow.conditioned,
        "eas_selected_path_frames",
        lambda _sizes, _path: (pd.DataFrame(), pd.DataFrame()),
    )
    monkeypatch.setattr(
        workflow.conditioned,
        "_reverse_path_log_components",
        lambda _sizes, _path: {
            "log_forward": -2.0,
            "log_reverse": -4.0,
            "log_mutation": 1.0,
            "log_weight": 3.0,
        },
    )
    replay_calls = []
    fail_blocks = {196}

    def fake_core_replay(_sizes, block_subset, group):
        block_index = int(block_subset[0]["block_index"])
        replay_calls.append(
            (
                block_index,
                tuple(group["unit_id"].astype(str)),
            )
        )
        if block_index in fail_blocks:
            raise RuntimeError("synthetic selected-block replay failure")
        return {unit_id: (14_704, 1) for unit_id in group["unit_id"].astype(str)}

    monkeypatch.setattr(
        workflow.conditioned, "replay_eas_global_selection", fake_core_replay
    )
    pool_workers = []
    submitted = []

    class ImmediateProcessPool:
        def __init__(self, *, max_workers):
            pool_workers.append(max_workers)

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        def submit(self, function, *args):
            submitted.append(int(args[2].rsplit("_", 1)[1]))
            future = Future()
            try:
                future.set_result(function(*args))
            except Exception as error:  # noqa: BLE001 - emulate ProcessPool failure
                future.set_exception(error)
            return future

    monkeypatch.setattr(workflow, "ProcessPoolExecutor", ImmediateProcessPool)
    monkeypatch.setattr(workflow, "as_completed", lambda futures: reversed(futures))
    replay_root = tmp_path / workflow.EAS_PARALLEL_REPLAY_DIRNAME
    with pytest.raises(RuntimeError, match="1 parallel EAS replays failed"):
        workflow._run_eas_parallel_selected_replays_v4_legacy_implementation(
            tmp_path,
            plan_dir=plan_root,
            proposal_dir=tmp_path / workflow.EAS_PROPOSAL_DIRNAME,
            replay_dir=replay_root,
            selected=selected,
            blocks=blocks,
            base_manifest=manifest,
            extension_plan=extension,
            max_workers=20,
        )
    assert (replay_root / "block_01" / "completion.json").is_file()
    assert not (replay_root / "block_196" / "completion.json").exists()
    assert not (replay_root / workflow.EAS_PARALLEL_REPLAY_COMPLETION_FILENAME).exists()
    assert not (tmp_path / workflow.EAS_SELECTION_DIRNAME).exists()
    fail_blocks.clear()
    result = workflow._run_eas_parallel_selected_replays_v4_legacy_implementation(
        tmp_path,
        plan_dir=plan_root,
        proposal_dir=tmp_path / workflow.EAS_PROPOSAL_DIRNAME,
        replay_dir=replay_root,
        selected=selected,
        blocks=blocks,
        base_manifest=manifest,
        extension_plan=extension,
        max_workers=20,
    )
    assert pool_workers == [2, 2]
    assert submitted == [1, 196, 1, 196]
    assert [call[0] for call in replay_calls] == [196, 196]
    assert list(result["paths"]) == selected["unit_id"].tolist()
    assert (
        result["paths"][selected.loc[0, "unit_id"]]
        == result["paths"][selected.loc[2, "unit_id"]]
    )
    assert result["worker_results"]["block_index"].tolist() == [1, 196]
    assert result["worker_results"]["cache_hit"].tolist() == [True, False]
    assert result["block_bindings"][0]["artifact_provenance_kind"] == (
        "exact_legacy_v3_migration"
    )
    assert result["block_bindings"][1]["artifact_provenance_kind"] == "current_v4"
    progress = capsys.readouterr().err
    assert progress.count("eas_parallel_replay_process_pool_start") == 2
    assert progress.count("eas_parallel_replay_block_complete") == 3
    assert progress.count("eas_parallel_replay_block_failed") == 1
    assert progress.count("eas_parallel_replay_process_pool_finished") == 2

    cached = workflow._run_eas_parallel_selected_replays_v4_legacy_implementation(
        tmp_path,
        plan_dir=plan_root,
        proposal_dir=tmp_path / workflow.EAS_PROPOSAL_DIRNAME,
        replay_dir=replay_root,
        selected=selected,
        blocks=blocks,
        base_manifest=manifest,
        extension_plan=extension,
        max_workers=20,
    )
    assert cached["cache_hit"] is True
    assert cached["worker_results"].empty
    assert submitted == [1, 196, 1, 196]

    (replay_root / "unexpected.txt").write_text("tamper\n", encoding="utf-8")
    with pytest.raises(ValueError, match="unexpected state"):
        workflow._run_eas_parallel_selected_replays_v4_legacy_implementation(
            tmp_path,
            plan_dir=plan_root,
            proposal_dir=tmp_path / workflow.EAS_PROPOSAL_DIRNAME,
            replay_dir=replay_root,
            selected=selected,
            blocks=blocks,
            base_manifest=manifest,
            extension_plan=extension,
        )
    (replay_root / "unexpected.txt").unlink()
    monkeypatch.setattr(
        workflow.conditioned,
        "_reverse_path_log_components",
        lambda _sizes, _path: {
            "log_forward": -2.0,
            "log_reverse": -4.0,
            "log_mutation": 1.0,
            "log_weight": 4.0,
        },
    )
    with pytest.raises(ValueError, match="cached replay path changed log_weight"):
        workflow._run_eas_parallel_selected_replays_v4_legacy_implementation(
            tmp_path,
            plan_dir=plan_root,
            proposal_dir=tmp_path / workflow.EAS_PROPOSAL_DIRNAME,
            replay_dir=replay_root,
            selected=selected,
            blocks=blocks,
            base_manifest=manifest,
            extension_plan=extension,
        )


def test_parallel_eas_replay_rejects_selected_provenance_tamper():
    selected, blocks = _synthetic_selected_and_blocks()
    group = workflow._selected_block_frame(selected, 1)
    workflow._validate_selected_group_against_block(group, blocks[1])
    for column, value in (
        ("unit_id", "wrong_unit"),
        ("cosi2_seed", workflow.conditioned.EAS_SEED_START + 99),
        ("global_proposal_index", -1),
        ("source_proposal_block_sha256", "f" * 64),
    ):
        tampered = group.copy()
        tampered.loc[0, column] = value
        with pytest.raises(ValueError):
            workflow._validate_selected_group_against_block(tampered, blocks[1])


def test_global_selection_verifier_binds_units_to_parallel_replay_paths(
    monkeypatch, tmp_path
):
    selected, blocks = _synthetic_selected_and_blocks()
    manifest = _synthetic_base_manifest()
    extension = _synthetic_extension_plan()
    replay_plan = {"plan_payload_sha256": "4" * 64}
    diagnostics = {"proposal_bank_binding_sha256": "5" * 64}
    replay_paths = {unit_id: (14_704, 1) for unit_id in selected["unit_id"].astype(str)}
    plan_root = tmp_path / "synthetic_v5_writer_plan"
    plan_root.mkdir()
    workflow._atomic_json(
        plan_root / workflow.PLAN_OUTPUTS["eas_base_bank_manifest"], manifest
    )
    workflow._atomic_json(
        plan_root / workflow.PLAN_OUTPUTS["eas_extension_plan"], extension
    )
    workflow._atomic_json(
        plan_root / workflow.PLAN_OUTPUTS["eas_parallel_replay_plan"], replay_plan
    )
    legacy_plan_root = tmp_path / workflow.LEGACY_V3_PLAN_DIRNAME
    legacy_plan_root.mkdir()
    workflow._atomic_json(
        legacy_plan_root / workflow.LEGACY_V3_PLAN_OUTPUTS["eas_parallel_replay_plan"],
        {"plan_payload_sha256": workflow.LEGACY_V3_REPLAY_PLAN_PAYLOAD_SHA256},
    )
    legacy_v4_plan_root = tmp_path / workflow.LEGACY_V4_PLAN_DIRNAME
    legacy_v4_plan_root.mkdir()
    workflow._atomic_json(
        legacy_v4_plan_root
        / workflow.LEGACY_V4_PLAN_OUTPUTS["eas_parallel_replay_plan"],
        {"plan_payload_sha256": workflow.LEGACY_V4_REPLAY_PLAN_PAYLOAD_SHA256},
    )
    workflow._atomic_json(
        plan_root / workflow.PLAN_COMPLETION_FILENAME, {"status": "synthetic"}
    )
    plan_completion_sha256 = workflow.sha256_file(
        plan_root / workflow.PLAN_COMPLETION_FILENAME
    )
    selection_root = tmp_path / workflow.EAS_SELECTION_DIRNAME
    selection_root.mkdir()
    workflow._atomic_frame(selection_root / "selected_occurrences.tsv", selected)
    workflow._atomic_json(selection_root / "global_diagnostics.json", diagnostics)
    unit_outputs = {}
    for unit_id in selected["unit_id"].astype(str):
        unit_dir = selection_root / "units" / unit_id
        unit_dir.mkdir(parents=True)
        (unit_dir / "completion.json").write_text("{}\n", encoding="utf-8")
        (unit_dir / "trajectory.tsv").write_text("synthetic\n", encoding="utf-8")
        unit_outputs[unit_id] = {
            "completion_sha256": workflow.sha256_file(unit_dir / "completion.json"),
            "trajectory_sha256": workflow.sha256_file(unit_dir / "trajectory.tsv"),
        }
    replay_completion_sha256 = "6" * 64
    contract = {
        "schema": workflow.EAS_BANK_WORKFLOW_SCHEMA,
        "proposal_bank_contract": (
            workflow.conditioned.EasProposalBankContract().to_record()
        ),
        "plan_completion_sha256": plan_completion_sha256,
        "legacy_v2_plan_completion_sha256": (workflow.LEGACY_V2_PLAN_COMPLETION_SHA256),
        "legacy_v3_plan_completion_sha256": (workflow.LEGACY_V3_PLAN_COMPLETION_SHA256),
        "legacy_v4_plan_completion_sha256": (workflow.LEGACY_V4_PLAN_COMPLETION_SHA256),
        "legacy_v4_replay_completion_sha256": (
            workflow.LEGACY_V4_REPLAY_COMPLETION_SHA256
        ),
        "legacy_v4_replay_root_manifest_sha256": (
            workflow.LEGACY_V4_REPLAY_ROOT_MANIFEST_SHA256
        ),
        "base_manifest_payload_sha256": manifest["manifest_payload_sha256"],
        "extension_plan_payload_sha256": extension["plan_payload_sha256"],
        "parallel_replay_plan_payload_sha256": replay_plan["plan_payload_sha256"],
        "parallel_replay_execution_method": workflow.EAS_PARALLEL_REPLAY_METHOD,
        "parallel_replay_completion_interpretation": (
            "hard_bound_legacy_v4_bytes_not_current_v5_completion"
        ),
        "parallel_replay_completion_sha256": replay_completion_sha256,
        "parallel_replay_unique_selected_block_count": 2,
        "workflow_sha256": workflow.sha256_file(workflow.MODULE_PATH),
        "proposal_bank_binding_sha256": diagnostics["proposal_bank_binding_sha256"],
        "proposal_completion_sha256": {
            str(block["block_index"]): str(block["completion_sha256"])
            for block in blocks
        },
        "conditioned_null_sha256": workflow.sha256_file(
            workflow.conditioned.MODULE_PATH
        ),
        "selection_method": "fixed_seed_categorical_multinomial_with_replacement",
        "resample_seed": workflow.conditioned.EAS_CATEGORICAL_RESAMPLE_SEED,
        "selected_count": workflow.EXPECTED_NEUTRAL_PER_DEMOGRAPHY,
        "selected_row_adapter": "set_index_drop_false_preserve_unit_id",
    }
    outputs = {
        "selected_occurrences": workflow._output_record(
            selection_root / "selected_occurrences.tsv", selection_root
        ),
        "global_diagnostics": workflow._output_record(
            selection_root / "global_diagnostics.json", selection_root
        ),
    }
    workflow._atomic_json(
        selection_root / workflow.EAS_SELECTION_COMPLETION,
        {
            "schema": workflow.EAS_SELECTION_SCHEMA,
            "status": "complete",
            "contract": contract,
            "contract_sha256": workflow._canonical_sha256(contract),
            "outputs": outputs,
            "unit_outputs": unit_outputs,
        },
    )
    monkeypatch.setattr(workflow, "verify_plan", lambda *_args, **_kwargs: {})
    monkeypatch.setattr(
        workflow, "load_eas_proposal_blocks", lambda *_args, **_kwargs: blocks
    )
    monkeypatch.setattr(
        workflow.conditioned,
        "select_eas_global_proposals",
        lambda *_args, **_kwargs: (selected, diagnostics),
    )
    monkeypatch.setattr(
        workflow,
        "verify_eas_parallel_replay_completion",
        lambda *_args, **_kwargs: {
            "completion_sha256": replay_completion_sha256,
            "unique_selected_block_count": 2,
            "paths": replay_paths,
        },
    )
    received_counts = {}

    def fake_verify_path(
        _root,
        bundle_dir,
        *,
        expected_selected_row,
        expected_counts,
        expected_global_diagnostics,
    ):
        unit_id = str(expected_selected_row["unit_id"])
        assert expected_global_diagnostics == diagnostics
        received_counts[unit_id] = tuple(expected_counts)
        return {"birth_generation": int(expected_selected_row["birth_generation"])}

    monkeypatch.setattr(
        workflow.conditioned,
        "verify_eas_selected_path_bundle",
        fake_verify_path,
    )
    verified = workflow.verify_eas_global_selection(
        tmp_path,
        plan_dir=plan_root,
        proposal_dir=tmp_path / workflow.EAS_PROPOSAL_DIRNAME,
        selection_dir=selection_root,
        replay_dir=tmp_path / workflow.EAS_PARALLEL_REPLAY_DIRNAME,
    )
    assert verified["cache_hit"] is True
    assert received_counts == replay_paths
    assert selection_root.name == "eas_global_selection_v5"
    assert not (tmp_path / workflow.LEGACY_V4_SELECTION_DIRNAME).exists()
    cached = workflow.verify_eas_global_selection(
        tmp_path,
        plan_dir=plan_root,
        proposal_dir=tmp_path / workflow.EAS_PROPOSAL_DIRNAME,
        selection_dir=selection_root,
        replay_dir=tmp_path / workflow.EAS_PARALLEL_REPLAY_DIRNAME,
    )
    assert cached["cache_hit"] is True


def test_v6_eas_selection_is_verify_only_exact_v5_reuse(monkeypatch, tmp_path):
    selected, _blocks = _synthetic_selected_and_blocks()
    diagnostics = {"proposal_bank_binding_sha256": "5" * 64}
    work = tmp_path / "work"
    plan_root = work / workflow.PLAN_DIRNAME
    replay_root = work / workflow.LEGACY_V5_PARALLEL_REPLAY_DIRNAME
    selection_root = work / workflow.LEGACY_V5_SELECTION_DIRNAME
    for path in (plan_root, replay_root, selection_root):
        path.mkdir(parents=True)
    state = {
        "selection_completion": {"schema": workflow.EAS_SELECTION_SCHEMA},
        "selected": selected,
        "diagnostics": diagnostics,
    }
    monkeypatch.setattr(
        workflow,
        "verify_plan",
        lambda *_args, **_kwargs: {
            "legacy_v5_selected": selected,
            "legacy_v5_diagnostics": diagnostics,
        },
    )
    monkeypatch.setattr(
        workflow,
        "_verify_legacy_v5_eas_state",
        lambda *_args, **_kwargs: state,
    )
    monkeypatch.setattr(
        workflow.conditioned,
        "select_eas_global_proposals",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("v6 must not resample the frozen v5 EAS bank")
        ),
    )
    verified = workflow.verify_eas_global_selection(
        tmp_path,
        plan_dir=plan_root,
        proposal_dir=work / workflow.EAS_PROPOSAL_DIRNAME,
        selection_dir=selection_root,
        replay_dir=replay_root,
    )
    assert verified["legacy_v5_reuse"] is True
    assert verified["selected"].equals(selected)


def test_inherited_han_cache_revalidates_source(monkeypatch, tmp_path):
    source = tmp_path / "natural" / "units" / "han_source"
    source.mkdir(parents=True)
    trajectory = source / "trajectory.tsv"
    trajectory.write_text("synthetic trajectory\n", encoding="utf-8")
    workflow._atomic_json(
        source / workflow.natural_neutral.UNIT_COMPLETION_FILENAME,
        {
            "status": "complete",
            "unit_id": "han_source",
            "outputs": {"trajectory.tsv": {"sha256": workflow.sha256_file(trajectory)}},
        },
    )
    candidate = {
        "candidate_order": 0,
        "candidate_id": "han_candidate_0000",
        "seed": 101,
        "source_kind": "existing_complete_natural_neutral_candidate",
        "source_unit_id": "han_source",
        "source_unit_dir": str(source),
    }
    validation = {
        "present_chb_af": 0.02,
        "han_generation_645_chb_af": 0.012,
        "han_generation_645_gate_passed": True,
        "present_strictly_segregating": True,
    }
    monkeypatch.setattr(
        workflow.conditioned,
        "validate_han_trajectory",
        lambda *_args, **_kwargs: validation,
    )
    binary = tmp_path / "coalescent"
    binary.write_bytes(b"synthetic")
    screen = tmp_path / "screen"
    workflow._run_han_screen_candidate(
        tmp_path,
        candidate,
        plan_dir=tmp_path / "unused-plan",
        screen_root=screen,
        binary=binary,
        max_attempts=1,
        timeout_seconds=1,
    )
    trajectory.write_text("tampered source\n", encoding="utf-8")
    with pytest.raises(ValueError, match="source binding"):
        workflow.verify_han_screen_candidate(
            candidate,
            plan_dir=tmp_path / "unused-plan",
            screen_dir=screen,
            cosi2_binary=binary,
        )


def test_real_cosi2_han_screen_resolves_plan_cwd_map_and_cache_binding(tmp_path):
    if os.name != "posix":
        pytest.skip("real CoSi2 smoke runs in the exact WSL production runtime")
    repo_root = Path(__file__).resolve().parents[1]
    binary = repo_root / workflow.DEFAULT_COSI2_BINARY_RELATIVE
    if not binary.is_file() or not os.access(binary, os.X_OK):
        pytest.skip("CoSi2 binary is unavailable")
    plan_root = _write_minimal_v6_han_plan(repo_root, tmp_path / "plan")
    candidate = workflow.conditioned.build_han_candidate_ledger().iloc[100].to_dict()
    screen_root = tmp_path / "screen"
    result = workflow._run_han_screen_candidate(
        repo_root,
        candidate,
        plan_dir=plan_root,
        screen_root=screen_root,
        binary=binary,
        max_attempts=10_000,
        timeout_seconds=180,
    )
    assert result["status"] == "complete"
    assert Path(result["trajectory_path"]).is_file()
    completion = workflow._read_json(
        Path(result["screen_completion_path"]), "real Han smoke completion"
    )
    assert completion["source_binding"]["recombination_map_sha256"] == (
        workflow.RECOMBINATION_MAP_SHA256
    )
    assert completion["source_binding"]["working_directory"] == str(plan_root)
    assert completion["command"][1] == "-p"
    plan_completion = plan_root / workflow.PLAN_COMPLETION_FILENAME
    plan_payload = workflow._read_json(plan_completion, "minimal v6 plan")
    plan_payload["test_stale_plan_marker"] = True
    workflow._atomic_json(plan_completion, plan_payload)
    with pytest.raises(ValueError, match="completion changed"):
        workflow.verify_han_screen_candidate(
            candidate,
            plan_dir=plan_root,
            screen_dir=screen_root,
            cosi2_binary=binary,
        )


def _accepted_han_unit(repo_root: Path, candidate_id: str) -> dict:
    accepted = pd.read_csv(
        repo_root
        / workflow.DEFAULT_WORK_RELATIVE
        / workflow.HAN_SCREEN_DIRNAME
        / workflow.HAN_ACCEPTED_FILENAME,
        sep="\t",
    )
    matches = accepted[accepted["candidate_id"].astype(str) == candidate_id]
    assert len(matches) == 1
    unit = matches.iloc[0].to_dict()
    unit["demography"] = "Han"
    unit["model_id"] = workflow.conditioned.HAN_MODEL_ID
    return unit


def test_han_loader_projection_is_exact_lexical_12_to_7_and_tamper_closed(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    unit = _accepted_han_unit(repo_root, "han_candidate_0005")
    source = Path(str(unit["trajectory_path"]))
    loader = tmp_path / "loader.tsv"
    metadata = workflow._write_han_loader_trajectory(source, loader)
    assert metadata["source_sha256"] == unit["trajectory_sha256"]
    assert metadata["loader_sha256"] == (
        "4a93374129a40bd1af4e57e06898b7980365db81709ea6ba6243d8100fcf8d0c"
    )
    source_lines = source.read_bytes().splitlines()
    loader_lines = loader.read_bytes().splitlines()
    assert len(source_lines) == len(loader_lines) == 2402
    assert (
        loader_lines[0]
        == b"sim\tgen\tselfreq_1\tselfreq_2\tselfreq_3\tselfreq_4\tselfreq_5"
    )
    assert all(
        loader_row.split(b"\t") == source_row.split(b"\t")[:7]
        for source_row, loader_row in zip(
            source_lines[1:], loader_lines[1:], strict=True
        )
    )
    loader.write_bytes(
        loader.read_bytes().replace(b"0.01138721555", b"0.01138721556", 1)
    )
    with pytest.raises(ValueError, match="exact lexical source projection"):
        workflow._validate_han_loader_trajectory(source, loader)
    malformed = tmp_path / "malformed.tsv"
    malformed.write_bytes(
        source.read_bytes().replace(b"\tpopsize_1", b"\r\tpopsize_1", 1)
    )
    with pytest.raises(ValueError, match="canonical UTF-8/LF"):
        workflow._han_loader_trajectory_projection(malformed)


def test_v7_failed_han_replay_retains_atomic_diagnostics(monkeypatch, tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    binary = repo_root / workflow.DEFAULT_COSI2_BINARY_RELATIVE
    plan_root = _write_minimal_v7_han_replay_plan(repo_root, tmp_path / "plan")
    unit = _accepted_han_unit(repo_root, "han_candidate_0005")
    run_root = tmp_path / workflow.SIMULATION_DIRNAME
    failure_root = tmp_path / workflow.SIMULATION_FAILURE_DIRNAME

    def fail_run(*_args, stdout, stderr, **_kwargs):
        stdout.write(b"partial-ms\n")
        stderr.write(b"bad traj file header\n")
        return SimpleNamespace(returncode=1)

    monkeypatch.setattr(workflow.subprocess, "run", fail_run)
    with pytest.raises(RuntimeError, match="conditioned replay exited 1"):
        workflow._run_conditioned_unit(
            repo_root,
            unit,
            plan_dir=plan_root,
            run_root=run_root,
            binary=binary,
            failure_root=failure_root,
            max_attempts=1_000_000,
            timeout_seconds=60,
        )
    attempt = failure_root / unit["unit_id"] / "attempt_0001"
    assert not (run_root / "units" / unit["unit_id"]).exists()
    assert not list((run_root / "units").glob(".*.stage.*"))
    assert (attempt / "simulation.stderr.log").read_bytes() == b"bad traj file header\n"
    assert (attempt / "loader_trajectory.tsv").is_file()
    failure = workflow._read_json(attempt / "failure.json", "retained failure")
    assert failure["schema"] == workflow.SIMULATION_FAILURE_SCHEMA
    assert failure["outputs"]["simulation.stderr.log"]["sha256"] == (
        workflow.sha256_file(attempt / "simulation.stderr.log")
    )


def test_v7_simulation_inventory_and_status_use_legacy_eas_native_han_union(tmp_path):
    work = tmp_path / "work"
    legacy_units = work / workflow.LEGACY_V6_SIMULATION_DIRNAME / "units"
    native_root = work / workflow.SIMULATION_DIRNAME
    native_units = native_root / "units"
    legacy_units.mkdir(parents=True)
    native_units.mkdir(parents=True)
    for index in range(workflow.EXPECTED_NEUTRAL_PER_DEMOGRAPHY):
        unit = legacy_units / f"{workflow.conditioned.EAS_MODEL_ID}_r{index:03d}"
        unit.mkdir()
        (unit / "completion.json").write_text("{}\n", encoding="utf-8")
    for index in range(3):
        unit = native_units / f"{workflow.conditioned.HAN_MODEL_ID}_r{index:03d}"
        unit.mkdir()
        (unit / "completion.json").write_text("{}\n", encoding="utf-8")
    status = workflow.workflow_status(tmp_path, work_dir=work)
    assert status["conditioned_simulations"] == {
        "complete": 103,
        "planned": 200,
        "legacy_v6_eas_complete": 100,
        "native_v7_han_complete": 3,
        "native_v7_root_completion": "absent",
    }
    expected_han = {
        f"{workflow.conditioned.HAN_MODEL_ID}_r{index:03d}"
        for index in range(workflow.EXPECTED_NEUTRAL_PER_DEMOGRAPHY)
    }
    for index in range(3, workflow.EXPECTED_NEUTRAL_PER_DEMOGRAPHY):
        (native_units / f"{workflow.conditioned.HAN_MODEL_ID}_r{index:03d}").mkdir()
    (native_root / workflow.SIMULATION_ROOT_COMPLETION_FILENAME).write_text(
        "{}\n", encoding="utf-8"
    )
    workflow._validate_v7_simulation_root_inventory(
        native_root, expected_han_ids=expected_han, require_complete=True
    )
    unexpected_eas = native_units / f"{workflow.conditioned.EAS_MODEL_ID}_r000"
    unexpected_eas.mkdir()
    with pytest.raises(ValueError, match="unit inventory changed"):
        workflow._validate_v7_simulation_root_inventory(
            native_root, expected_han_ids=expected_han, require_complete=True
        )


@pytest.mark.parametrize("candidate_id", ["han_candidate_0005", "han_candidate_0140"])
def test_real_cosi2_v7_han_full_replay_uses_lexical_loader_view(tmp_path, candidate_id):
    if os.name != "posix":
        pytest.skip("real CoSi2 smoke runs in the exact WSL production runtime")
    repo_root = Path(__file__).resolve().parents[1]
    binary = repo_root / workflow.DEFAULT_COSI2_BINARY_RELATIVE
    if not binary.is_file() or not os.access(binary, os.X_OK):
        pytest.skip("CoSi2 binary is unavailable")
    plan_root = _write_minimal_v7_han_replay_plan(repo_root, tmp_path / "plan")
    unit = _accepted_han_unit(repo_root, candidate_id)
    run_root = tmp_path / workflow.SIMULATION_DIRNAME
    result = workflow._run_conditioned_unit(
        repo_root,
        unit,
        plan_dir=plan_root,
        run_root=run_root,
        binary=binary,
        failure_root=tmp_path / workflow.SIMULATION_FAILURE_DIRNAME,
        max_attempts=1_000_000,
        timeout_seconds=300,
    )
    assert result["status"] == "complete"
    assert result["schema"] == workflow.SIMULATION_UNIT_SCHEMA
    unit_dir = run_root / "units" / unit["unit_id"]
    assert {path.name for path in unit_dir.iterdir()} == {
        "completion.json",
        "loader_trajectory.tsv",
        "model.params",
        "simulation.ms",
        "simulation.stderr.log",
        "trajectory.tsv",
    }
    completion = workflow._read_json(unit_dir / "completion.json", "v7 Han smoke")
    assert completion["environment_contract"]["COSI_LOAD_TRAJ"] == (
        "checksum_bound_lexical_loader_trajectory"
    )
    assert completion["loader_trajectory_validation"] == (
        workflow._validate_han_loader_trajectory(
            unit_dir / "trajectory.tsv", unit_dir / "loader_trajectory.tsv"
        )
    )
    assert Path(result["simulation_ms_path"]).stat().st_size > 1_000_000


def test_boolean_string_tamper_is_rejected():
    units = _units()
    units["eas_time_reversal_is_approximation"] = units[
        "eas_time_reversal_is_approximation"
    ].astype(object)
    mask = (units["demography"] == "EAS") & (units["simulation_class"] == "neutral")
    units.loc[mask, "eas_time_reversal_is_approximation"] = "False"
    with pytest.raises(ValueError, match="global-bank approximation"):
        workflow.validate_unit_inventory(units)


def test_atomic_results_cache_and_tamper(monkeypatch, tmp_path):
    units = _units()
    scores = _scores(units)
    decodes = _decodes(units)
    cosi2 = tmp_path / "coalescent"
    decoder = tmp_path / "gamma_smc"
    cosi2.write_bytes(b"cosi2")
    decoder.write_bytes(b"gamma")
    monkeypatch.setattr(workflow, "plot_conditioned_score_profiles", _fake_plot)
    monkeypatch.setattr(
        workflow,
        "plot_conditioning_diagnostics",
        lambda data, output_dir, *, stem=workflow.CONDITIONING_FIGURE_STEM: _fake_plot(
            data, output_dir, stem=stem
        ),
    )
    destination = tmp_path / "results"
    result = workflow.write_results_bundle(
        tmp_path,
        results_dir=destination,
        units=units,
        scores=scores,
        decode_records=decodes,
        decoder_bin=decoder,
        cosi2_binary=cosi2,
    )
    assert result["cache_hit"] is False
    assert {path.name for path in destination.iterdir()} == {
        workflow.RESULTS_COMPLETION_FILENAME,
        *workflow.RESULT_OUTPUTS.values(),
    }
    assert len(list(destination.glob("*.pdf"))) == 2
    contrasts = pd.read_csv(
        destination / workflow.RESULT_OUTPUTS["conditioning_contrasts"], sep="\t"
    )
    assert set(contrasts["demography"]) == {"EAS", "Han"}
    cached = workflow.write_results_bundle(
        tmp_path,
        results_dir=destination,
        units=units,
        scores=scores,
        decode_records=decodes,
        decoder_bin=decoder,
        cosi2_binary=cosi2,
    )
    assert cached["cache_hit"] is True
    target = destination / workflow.RESULT_OUTPUTS["scores"]
    target.write_bytes(target.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="checksum failed"):
        workflow.verify_results_bundle(tmp_path, results_dir=destination)


def test_cli_exposes_restartable_conditioned_phases():
    action = next(
        action for action in workflow.build_parser()._actions if action.dest == "phase"
    )
    assert set(action.choices) == {
        "plan",
        "eas-bridge",
        "han-screen",
        "simulate",
        "decode",
        "analyze",
        "status",
        "verify",
        "all",
    }


def test_default_han_screen_fills_through_later_partial_completion(
    monkeypatch, tmp_path
):
    work = tmp_path / "work"
    plan = work / workflow.PLAN_DIRNAME
    plan.mkdir(parents=True)
    ledger = workflow.conditioned.build_han_candidate_ledger()
    ledger.to_csv(plan / workflow.PLAN_OUTPUTS["han_candidates"], sep="\t", index=False)
    late = work / workflow.HAN_SCREEN_DIRNAME / "candidates" / "han_candidate_0400"
    late.mkdir(parents=True)
    (late / "completion.json").write_text("{}\n", encoding="utf-8")
    binary = tmp_path / "coalescent"
    binary.write_bytes(b"cosi")
    calls = []

    def fake_screens(*_args, candidate_ids, **_kwargs):
        calls.append(tuple(candidate_ids))
        return pd.DataFrame(
            {
                "candidate_id": candidate_ids,
                "status": "complete",
                "gate_passed": True,
            }
        )

    monkeypatch.setattr(workflow, "run_han_screens", fake_screens)
    monkeypatch.setattr(workflow, "verify_plan", lambda *_args, **_kwargs: {})
    monkeypatch.setattr(
        workflow,
        "freeze_han_acceptances",
        lambda *_args, **_kwargs: pd.DataFrame({"unit_id": range(100)}),
    )
    result = workflow.han_screen_study(
        tmp_path,
        work_dir=work,
        cosi2_binary=binary,
    )
    assert result["accepted"] == 100
    assert len(calls) == 5
    assert calls[-1][0] == "han_candidate_0400"
    assert calls[-1][-1] == "han_candidate_0499"
