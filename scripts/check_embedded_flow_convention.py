#!/usr/bin/env python3
"""Compare five embedded Gamma-SMC flow velocities with their source projection.

This deterministic numerical audit reads source files and uses SciPy to reproduce
the existing constant-size generator and its weighted least-squares projection.
It does not decode data, run simulations, or modify the production flow field.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import re
import warnings

import numpy as np
import scipy
from scipy.special import digamma, gammaincc, gammainccinv, gammaln, hyp1f1


PROBES = ((25, 49), (30, 42), (35, 35), (40, 49), (40, 42))
REPO = Path(__file__).resolve().parents[1]
SOURCE_PATHS = (
    "src/io.h",
    "src/generate_canonical_flow_field.cpp",
    "src/flow_field.h",
    "src/gamma_smc.cpp",
    "src/data_processor.h",
    "python/gamma_smc_aou/decoder.py",
)


def file_identity(path: Path) -> dict:
    content = path.read_bytes()
    return {"bytes": len(content), "sha256": hashlib.sha256(content).hexdigest()}


def read_embedded_table() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    text = (REPO / "src/io.h").read_text(encoding="utf-8")
    section = text.split("void read_flow_field_default(", 1)[1]
    section = section.split("void read_flow_field_raw(", 1)[0]
    grids = []
    for name in ("mean_grid_def", "cv_grid_def"):
        values = re.findall(rf"{name}\.push_back\(([^)]+)\)", section)
        grid = np.array([float(value) for value in values])
        if grid.size != 3 or grid[2] != int(grid[2]):
            raise ValueError(f"Unexpected embedded {name}: {values}")
        grids.append(grid)
    match = re.search(r"flow_field_unravelled\s*=\s*\{([^}]+)\}", section)
    if match is None:
        raise ValueError("Embedded flow field not found")
    field = np.array([float(value) for value in match.group(1).split(",")])
    expected = 2 * int(grids[0][2]) * int(grids[1][2])
    if field.size != expected or not np.isfinite(field).all():
        raise ValueError("Embedded flow-field size or values are invalid")
    return grids[0], grids[1], field


def projected_velocity(mean: float, cv: float, settings: dict) -> np.ndarray:
    alpha = 1.0 / cv**2
    beta = alpha / mean
    lower = settings["integration_min_time"]
    steps = settings["integration_n_steps"]
    gamma_end = gammainccinv(alpha, settings["gamma_tail_probability"]) / beta
    first_step = (gamma_end - lower) / steps
    exp_step = -math.log(settings["exponential_tail_probability"]) / (steps - 1)
    # The C++ grid advances once beyond the last gamma-grid point, then adds
    # the first exponential step. Keep that convention for this source audit.
    times = np.r_[
        lower + np.arange(steps + 1) * first_step,
        gamma_end + first_step + np.arange(1, steps + 1) * exp_step,
    ]
    weights = np.sqrt(np.r_[np.diff(times), 0.0])
    density = np.exp(
        alpha * np.log(beta)
        - gammaln(alpha)
        + (alpha - 1) * np.log(times)
        - beta * times
    )
    # This is the existing distribution_difference_pdf expression. It uses
    # 2*s*q(t|s) - 2*t*f(t), not the standard rho*s recombination convention.
    prefactor = np.exp(-times + alpha * np.log(times * beta) - gammaln(alpha + 1))
    change = (
        prefactor
        * (
            hyp1f1(alpha, alpha + 1, -(beta - 1) * times)
            - hyp1f1(alpha, alpha + 1, -(beta + 1) * times)
        )
        + (np.exp(-2 * times) / 2 - 0.5 - times) * density
        + (1 - np.exp(-2 * times)) * gammaincc(alpha, beta * times)
    )
    # Partial derivatives with respect to log10(alpha), log10(beta), including
    # their chain-rule factors. sqrt(delta_t) implements the L2 quadrature.
    columns = np.c_[
        density * (-digamma(alpha) + np.log(beta) + np.log(times)) * alpha * np.log(10),
        density * (alpha / beta - times) * beta * np.log(10),
    ]
    weighted_columns = columns * weights[:, None]
    weighted_change = change * weights
    if (
        not np.isfinite(weighted_columns).all()
        or not np.isfinite(weighted_change).all()
    ):
        raise ValueError("Nonfinite projection input")
    solution, _, rank, _ = np.linalg.lstsq(
        weighted_columns, weighted_change, rcond=None
    )
    if rank != 2:
        raise ValueError("Projection tangent matrix is rank deficient")
    return np.array([solution[0] - solution[1], -0.5 * solution[0]])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    mean_grid, cv_grid, field = read_embedded_table()
    settings = {
        "integration_min_time": 1e-10,
        "integration_n_steps": 1000,
        "gamma_tail_probability": 1e-3,
        "exponential_tail_probability": 1e-3,
        "log_coordinates": True,
        "least_squares_weight": "sqrt(delta_t); final point has zero weight",
        "special_function_backend": "scipy.special.hyp1f1; production generator uses Arb",
        "relative_ratio_tolerance": 1e-4,
    }
    means = np.linspace(
        np.log10(mean_grid[0]), np.log10(mean_grid[1]), int(mean_grid[2])
    )
    cvs = np.linspace(np.log10(cv_grid[0]), np.log10(cv_grid[1]), int(cv_grid[2]))
    cv_count = int(cv_grid[2])
    grid_size = int(mean_grid[2]) * cv_count
    probes = []
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        for row, col in PROBES:
            mean = float(10.0 ** means[row])
            cv = float(10.0 ** cvs[col])
            projected = projected_velocity(mean, cv, settings)
            embedded = field[[row * cv_count + col, grid_size + row * cv_count + col]]
            if np.any(np.abs(projected) < 1e-8):
                raise ValueError(
                    "Probe velocity too small for a stable ratio comparison"
                )
            ratio = embedded / projected
            probes.append(
                {
                    "mean_index_zero_based": row,
                    "cv_index_zero_based": col,
                    "mean_scaled_time": mean,
                    "cv": cv,
                    "alpha": 1.0 / cv**2,
                    "beta_rate": 1.0 / (cv**2 * mean),
                    "velocity_order": ["log10_mean", "log10_cv"],
                    "source_projection": projected.tolist(),
                    "embedded_velocity": embedded.tolist(),
                    "embedded_over_source": ratio.tolist(),
                    "max_absolute_ratio_error_from_one": float(
                        np.max(np.abs(ratio - 1))
                    ),
                }
            )
    errors = [probe["max_absolute_ratio_error_from_one"] for probe in probes]
    result = {
        "purpose": "Five-point embedded/source flow normalization audit; not decoder validation",
        "rng": "No random numbers used; probe indices are fixed",
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "source_files": {name: file_identity(REPO / name) for name in SOURCE_PATHS},
        "audit_script": file_identity(Path(__file__)),
        "scientific_assumptions": {
            "population_model": "Constant diploid effective population size N0",
            "scaled_tmrca": "t = physical generations / (2*N0)",
            "scaled_rates": "theta = 4*N0*mu; rho = 4*N0*r",
            "conditional_transition": "Pairwise SMC prime, including invisible recombinations",
            "standard_first_order_recombination_probability": "rho*s per physical base",
            "existing_source_generator": "Uses 2*s times conditional kernel, minus 2*t*f(t)",
            "comparison": "Embedded table versus this existing source formula, before runtime rho scaling",
            "runtime_operations": "No interpolation, cache building, clipping, or sequence decoding in this script",
            "scope": "Ratios near one exclude a global hidden 1/2 at these entries; not full-grid accuracy or power",
            "quadrature_origin": "Default settings in generate_canonical_flow_field.cpp; historical table build metadata unavailable",
        },
        "settings": settings,
        "embedded_mean_grid": mean_grid.tolist(),
        "embedded_cv_grid": cv_grid.tolist(),
        "embedded_value_count": int(field.size),
        "probes": probes,
        "max_absolute_ratio_error_from_one": max(errors),
        "passed": all(error < settings["relative_ratio_tolerance"] for error in errors),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    temporary.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    temporary.replace(args.output)
    print(
        json.dumps(
            {
                "output": str(args.output.resolve()),
                "probes": len(probes),
                "max_absolute_ratio_error_from_one": max(errors),
                "passed": result["passed"],
            }
        )
    )
    if not result["passed"]:
        raise SystemExit("Embedded/source convention check failed")


if __name__ == "__main__":
    main()
