"""Audit low-AF selected runs against saved trees, counts, seeds, and inputs."""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
import json
import os
from pathlib import Path

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "1"
import numpy as np
import pandas as pd
import tskit

from gamma_smc_aou.fresh_power import digest, ordered_nodes, seed_for, task_id

ROOT = Path("/mnt/d/phase2simselection/sim/eas_q02_h400")
OUT = ROOT.parent / "eas_lab_meeting_20260921/analysis/trajectory_audit"


def inspect(payload):
    row, cfg, fingerprint = payload
    path = ROOT/row["task_id"]
    record = json.loads((path/"simulated.json").read_text())
    task = record["task"]
    assert task_id(task) == row["task_id"]
    assert task["s"] == row["s"] and task["replicate"] == row["replicate"]
    assert record["fingerprint"] == fingerprint
    accepted = record["attempts"][-1]
    assert accepted["accepted"] and accepted["attempt"] == row["accepted_attempt"]
    assert seed_for(cfg, task, accepted["attempt"]) == row["seed"] == record["seed"]
    for name, column in (("simulation.trees", "simulation_sha256"),
                         ("decoded_input.trees", "decoder_input_sha256"),
                         ("focal_carriers.npy", "carriers_sha256")):
        assert digest(path/name) == row[column] == record["simulation_hashes"][name], (path, name)
    raw = tskit.load(path/"simulation.trees")
    cropped = tskit.load(path/"decoded_input.trees")
    carriers = np.load(path/"focal_carriers.npy", allow_pickle=False)
    assert len(carriers) == 400 and raw.num_samples == cropped.num_samples == 400
    assert np.all(raw.tables.nodes.time[raw.samples()] == 0)
    for ts, coordinate in ((raw, cfg["simulated_length_bp"]//2), (cropped, cfg["focal_position_bp"])):
        index = int(np.searchsorted(ts.tables.sites.position, coordinate))
        assert ts.site(index).position == coordinate
        variant = tskit.Variant(ts, samples=ordered_nodes(ts))
        variant.decode(index)
        assert len(variant.alleles) == 2
        np.testing.assert_array_equal(variant.genotypes != 0, carriers)
    af = float(carriers.mean())
    assert af == row["sample_af"] == record["sample_af"]
    meta = raw.metadata["SLiM"]["user_metadata"]
    census_af = meta["run7_final_census_af"][0]
    alt, total = meta["run7_final_census_alt_count"][0], meta["run7_final_census_total_count"][0]
    assert 0 <= alt <= total and alt/total == census_af
    assert meta["run7_placement_frequency"] == [1.0]
    assert meta["Q"] == [cfg["slim_scaling_factor"]]
    historical = [k for k in meta if any(x in k.lower() for x in ("trajectory", "frequency_history", "af_history"))]
    return dict(task_id=row["task_id"], s=row["s"], replicate=int(row["replicate"]),
        seed=int(row["seed"]), accepted_attempt=int(row["accepted_attempt"]),
        sample_af=af, census_af=census_af, census_alt=alt, census_total=total,
        allele_placement_tick=meta["run7_placement_tick"][0], final_tick=raw.metadata["SLiM"]["tick"],
        final_cycle=raw.metadata["SLiM"]["cycle"], historical_af_fields=";".join(historical),
        accepted_run_seconds=accepted["seconds"], archive_checks="passed")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((ROOT/"manifest.json").read_text())
    repo = Path(__file__).resolve().parents[1]
    for name in ("python/gamma_smc_aou/fresh_power.py", "python/gamma_smc_aou/run7_models.py", "python/gamma_smc_aou/run2_models.py"):
        assert digest(repo/name) == manifest["source_hashes"][name], name
    inventory = pd.read_csv(ROOT/"simulation_inventory.csv", float_precision="round_trip")
    selected = inventory[inventory["mode"] == "selected"].sort_values(["s", "replicate"])
    assert len(selected) == 1000 and (selected.groupby("s").size() == 100).all()
    payloads = [(r, manifest["config"], manifest["fingerprint"]) for r in selected.to_dict("records")]
    with ProcessPoolExecutor(max_workers=20) as pool:
        records = []
        for result in pool.map(inspect, payloads):
            records.append(result)
            if len(records) % 100 == 0:
                print(json.dumps(dict(audited=len(records), target=len(payloads))), flush=True)
    frame = pd.DataFrame(records)
    frame.to_csv(OUT/"selected_endpoint_audit.csv", index=False)
    tail = frame[(frame.s == .006) & (frame.sample_af < .5)].sort_values("sample_af")
    tail.to_csv(OUT/"s006_lower_tail.csv", index=False)
    summary = []
    for s, group in frame.groupby("s"):
        summary.append(dict(s=s, n=len(group), mean=group.sample_af.mean(), median=group.sample_af.median(),
            q05=group.sample_af.quantile(.05), q95=group.sample_af.quantile(.95),
            below_50=int((group.sample_af < .5).sum()), census_below_50=int((group.census_af < .5).sum())))
    pd.DataFrame(summary).to_csv(OUT/"af_tail_by_s.csv", index=False)
    report = dict(status="passed", selected_audited=len(frame), source_sha256=digest(Path(__file__)),
        input_inventory_sha256=digest(ROOT/"simulation_inventory.csv"), workers=20,
        s006_below_50=len(tail), s006_min_sample_af=float(tail.sample_af.min()),
        s006_min_sample_corresponding_census_af=float(tail.iloc[0].census_af),
        historical_af_fields_found=int(frame.historical_af_fields.str.len().gt(0).sum()),
        accepted_replay_seconds_sum=float(frame.accepted_run_seconds.sum()),
        interpretation="Endpoint integrity passed; stochastic historical paths cannot be inferred from terminal samples. Replay is needed to check the actual trajectory.",
        new_simulations=False, replacements=False)
    (OUT/"audit.json").write_text(json.dumps(report, indent=2)+"\n")
    print(tail.to_string(index=False), flush=True)
    print(json.dumps(report), flush=True)


if __name__ == "__main__":
    main()
