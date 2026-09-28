"""Verify this manuscript package locally without accessing cohort data or networks."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path


def main() -> None:
    root = Path(__file__).resolve().parent
    manifest = json.loads((root / "package_manifest.json").read_text(encoding="utf-8"))
    for name, expected in manifest["sha256"].items():
        path = root / name
        if path.parent != root or not path.is_file():
            raise ValueError(f"Invalid or missing package file: {name}")
        content = path.read_bytes()
        if hashlib.sha256(content).hexdigest() != expected:
            raise ValueError(f"Checksum mismatch: {name}")
        text = content.decode("utf-8")
        if "\ufffd" in text or "\x00" in text:
            raise ValueError(f"Invalid text content: {name}")
    summary = json.loads((root / "verified_summary.json").read_text(encoding="utf-8"))
    compact = summary["compact_sv"]
    assert compact["rows"] == compact["unique_event_ids"] == 38219
    assert sum(compact["sv_types"].values()) == compact["rows"]
    exclusive = summary["exclusive_sv"]
    assert sum(exclusive["functional_classes"].values()) == exclusive["unique_events"]
    assert exclusive["annotated_table_matches_existing_sha256"]
    assert summary["gamma_50k"]["threshold_generations"] == 2000
    assert abs(summary["gamma_50k"]["threshold_coalescent_units"] - 0.0688) < 1e-12
    assert not summary["gamma_50k"]["validated_empirical_results_available"]
    assert sum(summary["lct"]["n_by_population"].values()) == 12238
    demographic = json.loads((root / "phlash_summary.json").read_text(encoding="utf-8"))
    populations = demographic["populations"]
    assert set(populations) == {"AFR", "AMR", "EAS", "EUR", "MID", "SAS"}
    assert sum(p["num_fits"] for p in populations.values()) == 600
    assert not demographic["calendar_time_conversion_applied"]
    draft = (root / "introgression_selection_draft.md").read_text(encoding="utf-8")
    for name, population in populations.items():
        assert population["checks_passed"]
        assert population["num_time_points"] == 10000
        assert population["grid_min_generations"] == 100
        assert population["grid_max_generations"] == 40000
        for point in population["at_target_generations"].values():
            assert all(math.isfinite(value) for value in point.values())
            assert 0 < point["q025_ne"] <= point["median_ne"] <= point["q975_ne"]
        point = population["at_target_generations"]["2000"]
        assert abs(point["time_generations"] - 2000) < 1
        expected_row = (f"| {name} | {point['median_ne']:,.0f} | "
                        f"{point['q025_ne']:,.0f}--{point['q975_ne']:,.0f} |")
        assert expected_row in draft, f"PHLASH table differs from verified summary: {name}"
    print("Manuscript package checks passed; unresolved production fields remain explicitly marked.")


if __name__ == "__main__":
    main()
