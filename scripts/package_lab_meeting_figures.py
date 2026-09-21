"""Check figure files, create a portable ZIP, and optionally export to the repo."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import zipfile

from PIL import Image, ImageStat
from pypdf import PdfReader


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2)+"\n", encoding="utf-8", newline="\n")


def verify(root, visual_review):
    book = root/"EAS_lab_meeting_figures.pdf"
    pdfs = [book, *sorted((root/"figures").glob("*.pdf"))]
    pngs = sorted((root/"figures").glob("*.png"))
    if len(pdfs) != 16 or len(pngs) != 15:
        raise ValueError("Expected 15 figure pairs and one combined PDF")
    pdf_records, png_records = [], []
    for path in pdfs:
        reader = PdfReader(path, strict=True)
        assert len(reader.pages) == (15 if path == book else 1), path
        fonts = set()
        for page in reader.pages:
            assert len(page.images) == 0, (path, "raster image")
            assert len(page.extract_text()) > 140, (path, "missing text")
            assert abs(float(page.mediabox.width)/float(page.mediabox.height)-16/9) < .0001
            for obj in page["/Resources"].get("/Font", {}).values():
                font = obj.get_object()
                subtype = str(font["/Subtype"])
                fonts.add(subtype)
                assert subtype != "/Type3", path
                for face in font.get("/DescendantFonts", [font]):
                    desc = face.get_object().get("/FontDescriptor")
                    assert desc is not None, (path, "missing font descriptor")
                    assert any(k in desc.get_object() for k in ("/FontFile", "/FontFile2", "/FontFile3"))
        pdf_records.append(dict(file=path.relative_to(root).as_posix(), pages=len(reader.pages),
            sha256=digest(path), raster_images=0, font_types=sorted(fonts), fonts_embedded=True))
    for path in pngs:
        with Image.open(path) as img:
            img.load()
            assert img.size == (3199, 1800), path
            assert min(ImageStat.Stat(img.convert("RGB")).stddev) > 10, (path, "blank image")
            png_records.append(dict(file=path.relative_to(root).as_posix(), width=img.width, height=img.height))
    assert all(not row["text_outside_canvas"] for row in json.loads((root/"layout_checks.json").read_text()))
    if not visual_review and (root/"quality_checks.json").exists():
        old = json.loads((root/"quality_checks.json").read_text())
        # A manual review remains valid only for the exact same PDF files.
        if old.get("pdfs") == pdf_records:
            visual_review = old.get("visual_review")
    report = dict(status="passed", pdfs=pdf_records, pngs=png_records,
        visual_review=visual_review or "Not recorded for these PDF hashes; render and inspect before presentation.")
    write_json(root/"quality_checks.json", report)


def main(root, export_to, visual_review):
    verify(root, visual_review)
    audits = root/"analysis_audit"
    audits.mkdir(exist_ok=True)
    for source, prefix in ((root/"analysis", "focal"), (root/"analysis/positional", "positional"),
                           (root.parent/"eas_ihs_h400", "ihs")):
        for name in ("audit.json", "provenance.json", "analysis_provenance.json", "supplement_audit.json"):
            if (source/name).exists():
                shutil.copyfile(source/name, audits/f"{prefix}_{name}")
    tail_audit = root/"analysis/trajectory_audit"
    if (tail_audit/"audit.json").exists():
        assert json.loads((tail_audit/"audit.json").read_text())["status"] == "passed"
        for name in ("audit.json", "selected_endpoint_audit.csv", "s006_lower_tail.csv", "af_tail_by_s.csv"):
            shutil.copyfile(tail_audit/name, audits/f"af_tail_{name}")
    shutil.copyfile(root/"analysis/positional/calibration_thresholds.csv", root/"figure_data/calibration_thresholds.csv")
    (root/"REPRODUCE.md").write_text("""# Reproduce this figure collection

Repository: git@github.com:soisa001/gamma_smc_ts.git, branch AOU_run_opt.
The scripts run in WSL Ubuntu and use the saved study on D: at
`/mnt/d/phase2simselection/sim`. Source data are not included in Git.

```bash
git clone --branch AOU_run_opt git@github.com:soisa001/gamma_smc_ts.git
cd gamma_smc_ts
bash scripts/launch_lab_meeting_figures.sh build
bash scripts/launch_lab_meeting_figures.sh package
```

`build` renders the saved analyses; `package` checks all output PDFs/PNGs,
writes SHA-256 inventories, and verifies the ZIP. Both are idempotent.
Use `refresh` instead of `build` only when focal summaries need to be rebuilt:
it checks the completed archive, reuses verified focal/iHS caches, refreshes
positional calibration, and builds figures. It never starts SLiM or Gamma-SMC.
Focal extraction and iHS use 20 workers; plotting uses four profile readers
and one thread per numerical-library process. Seeds are preserved from the
simulation manifest and no figure uses random jitter or resampling.

Required saved inputs include `eas_q02_h400`, `eas_h400_pause_evaluation_20260917`,
`eas_allele_class_ablation`, `eas_joint_scan`, `eas_ihs_h400`,
`eas_positional_distance_h400`, and this collection's `analysis` directory.
The exact files and hashes used to construct figures are recorded in
`figure_provenance.json`. Original parameter and analysis audits accompany
the package. ZIP contents are only figures, summarized simulation data, captions,
and provenance; tree sequences and per-site raw decoding files stay on D:.

After changing figure code, render the combined PDF with Poppler and inspect
every page. The PDF/PNG checks cannot replace visual inspection. A manual
review recorded in `quality_checks.json` applies only to its exact PDF hashes.
""", encoding="utf-8", newline="\n")
    names = ["README.md", "REPRODUCE.md", "figure_index.json", "figure_provenance.json",
             "layout_checks.json", "quality_checks.json", "EAS_lab_meeting_figures.pdf"]
    if (root/"AF_TAIL_AUDIT.md").exists():
        names.append("AF_TAIL_AUDIT.md")
    paths = [root/name for name in names]
    for folder in ("figures", "figure_data", "analysis_audit"):
        paths.extend(sorted((root/folder).glob("*")))
    assert all(path.is_file() for path in paths)
    manifest = dict(package="EAS lab meeting figures, 2026-09-21", source_sha256=digest(Path(__file__)),
        files={p.relative_to(root).as_posix(): dict(sha256=digest(p), bytes=p.stat().st_size) for p in paths})
    write_json(root/"artifact_manifest.json", manifest)
    paths.append(root/"artifact_manifest.json")
    archive = root/"EAS_lab_meeting_figure_pack.zip"
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zipped:
        for path in paths:
            info = zipfile.ZipInfo(path.relative_to(root).as_posix(), date_time=(2026, 9, 21, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            zipped.writestr(info, path.read_bytes())
    with zipfile.ZipFile(archive) as zipped:
        assert zipped.testzip() is None
        for name, entry in manifest["files"].items():
            assert hashlib.sha256(zipped.read(name)).hexdigest() == entry["sha256"], name
    (root/(archive.name+".sha256")).write_text(f"{digest(archive)}  {archive.name}\n", encoding="ascii", newline="\n")
    if export_to:
        for path in paths:
            target = export_to/path.relative_to(root)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(path, target)
            assert digest(target) == digest(path), target
    print(json.dumps(dict(status="passed", files=len(paths), zip=str(archive),
        zip_bytes=archive.stat().st_size, exported_to=str(export_to) if export_to else None)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("/mnt/d/phase2simselection/sim/eas_lab_meeting_20260921"))
    parser.add_argument("--export-to", type=Path)
    parser.add_argument("--visual-review-note")
    args = parser.parse_args()
    main(args.root, args.export_to, args.visual_review_note)
