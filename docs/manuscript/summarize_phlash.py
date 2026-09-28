"""Validate published aggregate PHLASH artifacts and summarize resampling fits.

Requires NumPy; never opens individual-level data or pickle files. No sampling,
model fitting, calendar-time rescaling, or endpoint extrapolation is performed.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import zipfile

import numpy as np

POPULATIONS = ("AFR", "AMR", "EAS", "EUR", "MID", "SAS")
TARGETS = (100, 500, 1000, 2000, 10000, 40000)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def summarize(root, population):
    path = root / f"{population}.npz"
    sidecar = path.with_suffix(".json")
    meta = json.loads(sidecar.read_text(encoding="utf-8"))
    with zipfile.ZipFile(path) as archive:
        assert archive.testzip() is None, path
    with np.load(path, allow_pickle=False) as data:
        assert str(data["schema"]) == "phlash.aou.log-ne-mvn/v1"
        assert str(data["population"]) == population
        time = np.asarray(data["time"], dtype=float)
        curves = np.asarray(data["bootstrap_ne"], dtype=float)
        mean = np.asarray(data["mean_log_ne"], dtype=float)
        factor = np.asarray(data["covariance_factor"], dtype=float)
        jitter = float(data["jitter"])
    assert np.array_equal(time, np.geomspace(100, 40000, 10000).astype(np.float32))
    assert curves.shape == factor.shape == (100, 10000)
    assert mean.shape == (10000,)
    assert all(np.isfinite(a).all() for a in (time, curves, mean, factor))
    assert (curves > 0).all() and jitter == 0
    log_ne = np.log(curves)
    mean_error = float(np.max(np.abs(mean - log_ne.mean(axis=0))))
    factor_error = float(np.max(np.abs(factor - (log_ne - log_ne.mean(axis=0)) / np.sqrt(99))))
    assert mean_error < 2e-6 and factor_error < 2e-7
    assert meta["artifact_sha256"] == sha256(path)
    assert meta["artifact_size_bytes"] == path.stat().st_size
    assert meta["population"] == population
    assert meta["schema"] == "phlash.aou.log-ne-mvn/v1"
    assert meta["num_fits"] == 100 and meta["num_time_points"] == 10000
    assert meta["min_time_generations"] == 100 and meta["max_time_generations"] == 40000
    assert meta["covariance_rank_at_most"] == 99 and meta["jitter_log_ne_sd"] == 0
    assert sorted(Path(s["path"]).name for s in meta["sources"]) == [f"fit{i:03d}.pkl" for i in range(100)]
    quantiles = np.quantile(curves, [0.025, 0.5, 0.975], axis=0, method="linear")

    def at(index):
        return dict(time_generations=float(time[index]),
                    q025_ne=float(quantiles[0, index]),
                    median_ne=float(quantiles[1, index]),
                    q975_ne=float(quantiles[2, index]))

    result = {
        "artifact_sha256": sha256(path), "sidecar_sha256": sha256(sidecar),
        "artifact_size_bytes": path.stat().st_size,
        "num_fits": 100, "num_time_points": 10000,
        "grid_min_generations": 100, "grid_max_generations": 40000,
        "source_fit_directory": str(Path(meta["sources"][0]["path"]).parent).replace("\\", "/"),
        "mean_reconstruction_max_abs_error": mean_error,
        "factor_reconstruction_max_abs_error": factor_error,
        "minimum_of_pointwise_median": at(int(np.argmin(quantiles[1]))),
        "maximum_of_pointwise_median": at(int(np.argmax(quantiles[1]))),
        "at_target_generations": {str(target): at(int(np.argmin(np.abs(time - target)))) for target in TARGETS},
        "checks_passed": True,
    }
    return population, result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("--threads must be positive")
    print(f"Reading aggregate artifacts: {args.artifact_dir.resolve()}; threads={args.threads}")
    with ThreadPoolExecutor(max_workers=args.threads) as pool:
        populations = dict(pool.map(lambda p: summarize(args.artifact_dir, p), POPULATIONS))
    result = {
        "schema_version": 1,
        "source_repository": "https://github.com/soisa001/phlash_ld",
        "source_revision_reviewed": "70e177ea0e77fe38d408eb8042e13a44c94d15c9",
        "artifact_relative_directory": "artifacts/aou-mvn",
        "summary_definition": "Pointwise median and linear 2.5th/97.5th percentiles of 100 resampling-fit posterior-median curves; nearest stored grid point to each target generation.",
        "uncertainty": "Between-fit resampling spread; not pooled within-fit posterior uncertainty or simultaneous confidence bands.",
        "calendar_time_conversion_applied": False,
        "original_fit_metadata_available": False,
        "numpy_version": np.__version__,
        "populations": populations,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes((json.dumps(result, indent=2, sort_keys=True) + "\n").encode("utf-8"))
    print(f"Validated six populations / 600 fit curves; wrote {args.output.resolve()}")
    for population, values in populations.items():
        print(population, json.dumps(values["at_target_generations"]["2000"]))


if __name__ == "__main__":
    main()
