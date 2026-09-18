"""Export verified summary tables; simulations and per-site iHS remain on D."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def export(source, destination, names, supplement=()):
    audit = json.loads((source / "audit.json").read_text())
    if audit["status"] != "passed":
        raise ValueError(f"Analysis has not passed its audit: {source}")
    manifest = json.loads((source / "artifact_manifest.json").read_text())
    if supplement:
        manifest.update(json.loads((source / "supplement_manifest.json").read_text()))
    destination.mkdir(parents=True, exist_ok=True)
    copied = []
    for name in (*names, *supplement):
        entry = manifest[name]
        if (source/name).stat().st_size != entry["bytes"] or digest(source/name) != entry["sha256"]:
            raise ValueError(f"Corrupt summary: {source/name}")
        shutil.copyfile(source/name, destination/name)
        copied.append(name)
    shutil.copyfile(source / "artifact_manifest.json", destination / "full_analysis_manifest.json")
    copied.append("full_analysis_manifest.json")
    if supplement:
        shutil.copyfile(source / "supplement_manifest.json", destination / "full_supplement_manifest.json")
        copied.append("full_supplement_manifest.json")
    exported = dict(source=str(source), source_manifest_sha256=digest(source/"artifact_manifest.json"),
        exporter_sha256=digest(Path(__file__)),
        files={name: dict(bytes=(destination/name).stat().st_size, sha256=digest(destination/name)) for name in copied})
    (destination / "artifact_manifest.json").write_text(json.dumps(exported, indent=2)+"\n")
    print(json.dumps(dict(destination=str(destination), files=len(copied)+1,
        bytes=sum(item["bytes"] for item in exported["files"].values()))))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sim-root", type=Path, default=Path("/mnt/d/phase2simselection/sim"))
    parser.add_argument("--destination", type=Path, default=Path(__file__).resolve().parents[1]/"docs/results")
    parser.add_argument("--analysis", choices=("positional", "distance", "ihs", "all"), default="all")
    args = parser.parse_args()
    if args.analysis in ("positional", "all"):
        export(args.sim_root / "eas_positional_h400", args.destination / "positional_eas_evaluation",
            ("metrics.csv", "audit.json", "provenance.json"),
            ("calibration_thresholds.csv", "paired_af_comparison.csv", "pooled_rank_metrics.csv", "supplement_audit.json"))
    if args.analysis in ("ihs", "all"):
        export(args.sim_root / "eas_ihs_h400", args.destination / "ihs_eas_h400",
            ("metrics.csv", "audit.json", "analysis_provenance.json", "regions.csv"))
    if args.analysis in ("distance", "all"):
        export(args.sim_root / "eas_positional_distance_h400", args.destination / "positional_distance_h400",
            ("selected_distances.csv", "audit.json", "provenance.json"))
