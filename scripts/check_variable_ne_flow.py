#!/usr/bin/env python3
"""Deterministic numerical checks for the variable-Ne research derivation.

This is a quadrature calculation, not a decoder or simulation workflow. It does
not load simulation data or alter the production constant-Ne flow field.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import platform

import numpy as np
import scipy
from scipy.integrate import quad
from scipy.stats import gamma


class History:
    def __init__(self, starts, rates):
        self.starts = tuple(starts)
        self.rates = tuple(rates)
        assert self.starts[0] == 0 and len(self.starts) == len(self.rates)
        assert all(r > 0 for r in self.rates)
        assert all(a < b for a, b in zip(self.starts, self.starts[1:]))

    def quantities(self, t):
        hazard = 0.0
        j = 0.0
        for i, (start, rate) in enumerate(zip(self.starts, self.rates)):
            end = self.starts[i + 1] if i + 1 < len(self.starts) else math.inf
            width = max(0.0, min(t, end) - start)
            hazard += rate * width
            decay = math.exp(-2.0 * rate * width)
            j = j * decay - math.expm1(-2.0 * rate * width) / (2.0 * rate)
            if t < end:
                return rate, hazard, j
        raise AssertionError("history lacks a tail")

    def prior(self, t):
        rate, hazard, _ = self.quantities(t)
        return rate * math.exp(-hazard)

    def integral(self, function, lower=0.0, upper=math.inf, extra=()):
        bounds = sorted(
            {
                lower,
                upper,
                *(x for x in self.starts if lower < x < upper),
                *(x for x in extra if lower < x < upper),
            }
        )
        return sum(
            quad(function, a, b, epsabs=2e-11, epsrel=2e-11, limit=200)[0]
            for a, b in zip(bounds, bounds[1:])
        )

    def kernel(self, s, t):
        # Twice the standard generator: exactly the local canonical generator's
        # convention. A prospective decoder must resolve rho scaling separately.
        rate_t, ht, _ = self.quantities(t)
        _, hs, _ = self.quantities(s)
        _, _, jmin = self.quantities(min(s, t))
        return 2.0 * rate_t * jmin * math.exp(-max(ht - hs, 0.0))

    def exit_rate(self, s):
        return s + self.quantities(s)[2]

    def apply(self, density, t):
        inflow = self.integral(lambda s: self.kernel(s, t) * density(s), extra=(t,))
        return inflow - self.exit_rate(t) * density(t)


def max_abs(values):
    return max(abs(float(value)) for value in values)


def check_history(history):
    positions = np.geomspace(0.002, 8.0, 14)
    row_errors = []
    balance_errors = []
    j_errors = []
    for s in positions:
        row_errors.append(
            history.integral(lambda t: history.kernel(s, t), extra=(s,))
            - history.exit_rate(s)
        )
        _, hs, js = history.quantities(s)
        integrated_j = history.integral(
            lambda u: math.exp(-2 * (hs - history.quantities(u)[1])), upper=s
        )
        j_errors.append(js - integrated_j)
        for t in positions:
            lhs = history.prior(s) * history.kernel(s, t)
            rhs = history.prior(t) * history.kernel(t, s)
            balance_errors.append(lhs - rhs)
    stationary = [history.apply(history.prior, t) for t in positions]
    mean = history.integral(lambda t: t * history.prior(t))
    second = history.integral(lambda t: t * t * history.prior(t))
    return {
        "epoch_starts_scaled": history.starts,
        "coalescent_rates": history.rates,
        "prior_integral_error": abs(history.integral(history.prior) - 1),
        "j_direct_integral_max_abs_error": max_abs(j_errors),
        "exit_rate_integral_max_abs_error": max_abs(row_errors),
        "detailed_balance_max_abs_error": max_abs(balance_errors),
        "stationary_generator_max_abs_error": max_abs(stationary),
        "prior_mean_scaled": mean,
        "prior_cv": math.sqrt(second - mean * mean) / mean,
    }


def tilted_density(history, a, b):
    def unnormalized(t):
        if t == 0:
            return 0.0 if a > 0 else history.prior(0)
        rate, hazard, _ = history.quantities(t)
        return math.exp(math.log(rate) - hazard + a * math.log(t) - b * t)

    z = history.integral(unnormalized)
    return lambda t: unnormalized(t) / z


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    constant = History([0], [1])
    variable = History([0, 0.05, 0.15, 0.4], [0.25, 4.0, 0.5, 1.0])
    positions = np.geomspace(0.002, 6, 16)
    kernel_errors = []
    generator_errors = []
    gamma_errors = []
    for s in positions:
        for t in positions:
            expected = -math.expm1(-2 * min(s, t)) * math.exp(-max(t - s, 0))
            kernel_errors.append(constant.kernel(s, t) - expected)
    for alpha, beta in [(1, 1), (2, 3), (6, 1.2)]:
        density = lambda t: gamma.pdf(t, a=alpha, scale=1 / beta)
        tilt = tilted_density(constant, alpha - 1, beta - 1)
        for t in positions:
            integral = quad(
                lambda s: (math.exp(s - t) - math.exp(-s - t)) * density(s),
                0,
                t,
                epsabs=2e-11,
                epsrel=2e-11,
            )[0]
            expected = (
                integral
                + (math.exp(-2 * t) / 2 - 0.5 - t) * density(t)
                + (-math.expm1(-2 * t)) * gamma.sf(t, a=alpha, scale=1 / beta)
            )
            generator_errors.append(constant.apply(density, t) - expected)
            gamma_errors.append(tilt(t) - density(t))

    # Independently normalize products, verifying the proposed family is closed
    # under emissions and the prior-corrected forward/backward combination.
    q1 = tilted_density(variable, 0.7, 0.4)
    q2 = tilted_density(variable, 1.3, 0.8)
    combined = tilted_density(variable, 2.0, 1.2)
    emission_updated = tilted_density(variable, 1.7, 0.43)
    product = lambda t: q1(t) * q2(t) / variable.prior(t) if t < 500 else 0.0
    emission_product = lambda t: q1(t) * t * math.exp(-0.03 * t)
    z_product = variable.integral(product)
    z_emission = variable.integral(emission_product)
    merge_error = max_abs(product(t) / z_product - combined(t) for t in positions)
    emission_error = max_abs(
        emission_product(t) / z_emission - emission_updated(t) for t in positions
    )
    result = {
        "purpose": "Deterministic algebra/quadrature checks only; not a variable-Ne decoder",
        "seed": "No RNG used",
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "generator_convention": "Twice standard rho-scaled SMC prime, matching existing local generator",
        "absolute_tolerance": 1e-8,
        "histories": {
            "constant": check_history(constant),
            "illustrative_variable": check_history(variable),
        },
        "constant_kernel_max_abs_error": max_abs(kernel_errors),
        "constant_generator_max_abs_error": max_abs(generator_errors),
        "constant_tilt_gamma_max_abs_error": max_abs(gamma_errors),
        "variable_tilt_emission_max_abs_error": emission_error,
        "variable_tilt_merge_max_abs_error": merge_error,
        "entropy_counterexample": {
            "theta": 0.00075,
            "prior_entropy": 1.0,
            "one_heterozygote_posterior_entropy": float(
                gamma.entropy(a=2, scale=1 / 1.00075)
            ),
        },
    }
    errors = [v for k, v in result.items() if k.endswith("_error")]
    errors.extend(
        v
        for h in result["histories"].values()
        for k, v in h.items()
        if k.endswith("_error")
    )
    result["passed"] = all(error < result["absolute_tolerance"] for error in errors)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    temporary.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    temporary.replace(args.output)
    print(json.dumps(result, indent=2))
    if not result["passed"]:
        raise SystemExit("Numerical feasibility check failed")


if __name__ == "__main__":
    main()
