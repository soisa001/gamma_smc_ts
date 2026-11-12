"""Freeze a completed blind submission and grade it without modifying the solver files."""
import argparse
import json
from pathlib import Path
import shutil

import benchmark as b
import grade


def capture(agent_directory, task_directory, release_manifest, output):
    agent_directory, task_directory, output = map(Path, (agent_directory, task_directory, output))
    manifest = json.loads(Path(release_manifest).read_text())
    for name, expected in manifest["files"].items():
        if b.sha256(agent_directory/name) != expected:
            raise ValueError(f"Agent input changed: {name}")
        if b.sha256(task_directory/"public"/name) != expected:
            raise ValueError(f"Evaluator input mismatch: {name}")
    for required in ("REPORT.md", "summary.json", "decision_log.md", "predictions.tsv"):
        if not (agent_directory/required).is_file():
            raise ValueError(f"Incomplete submission: {required}")
    output.mkdir(parents=True, exist_ok=True)
    copied = {}
    source_dir = output/"agent_submission"
    source_dir.mkdir(exist_ok=True)
    for path in sorted(agent_directory.iterdir()):
        if not path.is_file() or path.name in manifest["files"] or path.name == "predictions.tsv":
            continue
        if path.suffix not in (".py", ".md", ".json", ".txt", ".toml", ".lock", ".sh", ".tsv", ".csv", ".log", ".npz", ".npy"):
            continue
        if path.stat().st_size > 10*1024*1024:
            if path.suffix == ".npy":
                continue  # Regenerable per-site probability caches remain local.
            raise ValueError(f"Unexpectedly large source/report file: {path}")
        shutil.copyfile(path, source_dir/path.name)
        copied[path.name] = b.sha256(path)
    prediction = agent_directory/"predictions.tsv"
    frozen_hash = b.sha256(prediction)
    b.write_json(output/"submission_manifest.json", {"agent_directory": str(agent_directory.resolve()),
                 "input_manifest_sha256": b.sha256(release_manifest), "inputs": manifest["files"],
                 "prediction_sha256": frozen_hash, "prediction_bytes": prediction.stat().st_size,
                 "source_and_report_files": copied, "freeze_policy": "No hidden feedback or revisions before grading"})
    score = grade.grade(task_directory, prediction, output/"grade.json")
    if b.sha256(prediction) != frozen_hash:
        raise ValueError("Prediction file changed during grading")
    with prediction.open() as handle, (output/"predictions_head.tsv").open("w", newline="") as sample:
        for _ in range(21):
            line = handle.readline()
            if not line:
                break
            sample.write(line)
    return score


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("agent_directory", type=Path)
    parser.add_argument("task_directory", type=Path)
    parser.add_argument("release_manifest", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(capture(args.agent_directory, args.task_directory, args.release_manifest, args.output))
