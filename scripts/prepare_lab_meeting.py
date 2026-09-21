"""Refresh focal truth from the completed archive, reusing verified cached summaries."""
from concurrent.futures import ProcessPoolExecutor
import json
import os
from pathlib import Path

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "1"
import numpy as np
import pandas as pd

from evaluate_paused_array import focal_one, digest, atomic_json

ROOT = Path("/mnt/d/phase2simselection/sim")
OUT = ROOT / "eas_lab_meeting_20260921/analysis"


def extract(payload):
    row, cfg, pairs, destination = payload
    # New records are checked directly against the audited archive inventory.
    if destination == OUT:
        for name, column in (("decoder_input_path", "decoder_input_sha256"), ("carriers_path", "carriers_sha256")):
            if digest(Path(row[name])) != row[column]:
                raise ValueError(f"Corrupt input: {row['task_id']}")
    return focal_one(payload)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    study = ROOT / "eas_q02_h400"
    status = json.loads((study / "simulation_status.json").read_text())
    assert status["state"] == "complete" and status["completed"] == 2000 and not status["failed"]
    inventory = pd.read_csv(study / "simulation_inventory.csv", float_precision="round_trip")
    selected = inventory[inventory["mode"] == "selected"].sort_values(["s", "replicate"])
    assert len(inventory) == 2000 and len(selected) == 1000
    assert (selected.groupby("s").size() == 100).all()
    cfg = json.loads((study / "manifest.json").read_text())["config"]
    pairs = np.loadtxt(study / "pairs.tsv", dtype=int)
    old = ROOT / "eas_h400_pause_evaluation_20260917"
    payloads, reused = [], 0
    for row in selected.to_dict("records"):
        present = (old / "cache" / (row["task_id"].replace("/", "_")+".json")).exists()
        reused += int(present)
        payloads.append((row, cfg, pairs, old if present else OUT))
    with ProcessPoolExecutor(max_workers=20) as pool:
        results = list(pool.map(extract, payloads))
    frame = pd.DataFrame(results)
    frame.to_csv(OUT / "focal_truth.csv", index=False)
    atomic_json(OUT / "audit.json", dict(status="passed", selected_regions=len(frame),
        selected_per_coefficient=100, cached_focal_records_verified=reused,
        new_focal_records_verified=len(frame)-reused, new_simulations=False, new_decoding=False))
    atomic_json(OUT / "provenance.json", dict(config=cfg, source_sha256=digest(Path(__file__)),
        extraction_source_sha256=digest(Path(__file__).with_name("evaluate_paused_array.py")),
        inventory_sha256=digest(study / "simulation_inventory.csv"), original_cache=str(old / "cache")))
    names = ("focal_truth.csv", "audit.json", "provenance.json")
    atomic_json(OUT / "artifact_manifest.json", {n: dict(sha256=digest(OUT/n), bytes=(OUT/n).stat().st_size) for n in names})
    print(json.dumps(dict(selected_regions=len(frame), reused=reused, newly_extracted=len(frame)-reused)), flush=True)


if __name__ == "__main__":
    main()
