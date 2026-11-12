"""Package only the solver-visible files; fixed ZIP metadata makes releases repeatable."""
import argparse
from pathlib import Path
import zipfile

import benchmark as b

PUBLIC_FILES = ("TASK.md", "metadata.json", "calibration_parent.fasta",
                "calibration_offspring.fasta", "new_parent.fasta", "queries.tsv")


def release(task, output):
    task, output = Path(task), Path(output)
    if not b.manifest_valid(task, "samples.manifest.json"):
        raise ValueError("Missing or corrupt sample files")
    output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for name in PUBLIC_FILES:
            info = zipfile.ZipInfo(name, date_time=(2026, 11, 12, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            archive.writestr(info, (task / "public" / name).read_bytes())
    with zipfile.ZipFile(output) as archive:
        if archive.testzip() is not None or tuple(archive.namelist()) != PUBLIC_FILES:
            raise ValueError("Release archive failed integrity validation")
    summary = {"archive": output.name, "archive_sha256": b.sha256(output),
               "files": {name: b.sha256(task / "public" / name) for name in PUBLIC_FILES}}
    b.write_json(output.with_suffix(".manifest.json"), summary)
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("task", type=Path)
    parser.add_argument("--output", type=Path, default=Path("releases/task_001.zip"))
    args = parser.parse_args()
    print(release(args.task, args.output))
