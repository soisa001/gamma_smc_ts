#!/usr/bin/env python
"""Evaluate run8's coalescent predictions and check them against the simulations.

Predictions are computed from :mod:`run8_theory` independently of the simulated
data; whatever has finished simulating is then compared to them. Anything still
running simply prints as unavailable rather than blocking the rest.
"""

from __future__ import annotations

import argparse
import glob
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "python"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from gamma_smc_aou.run7_config import (  # noqa: E402
    ADMIXTURE_PROPORTION,
    GENERATION_TIME_YEARS,
    PULSE_GENERATIONS,
    SPLIT_GENERATIONS,
)
from gamma_smc_aou.run8_config import SCENARIOS, THRESHOLDS_YEARS  # noqa: E402
from gamma_smc_aou.run8_theory import (  # noqa: E402
    establishment_frequency,
    load_trajectory,
    neutral_p_tmrca_below,
    predicted_final_frequency,
    sweep_coalescence_probability,
    sweep_coalescence_probability_trajectory,
    sweep_duration_generations,
)

STUDY = REPO / "sim_results_run8"


def observed_neutral(thresholds_years) -> np.ndarray | None:
    paths = sorted(glob.glob(str(STUDY / "neutral" / "replicates" / "*" / "tmrca_matrix.npz")))
    if not paths:
        return None
    values = []
    for path in paths:
        with np.load(path) as data:
            tmrca = np.asarray(data["tmrca_generations"], dtype=float)
        # Every position, not just the focal base: the neutral profile is
        # stationary along the contig, so this is simply a larger sample.
        values.append(
            [
                float(np.mean(tmrca < float(y) / GENERATION_TIME_YEARS))
                for y in thresholds_years
            ]
        )
    return np.asarray(values).mean(axis=0)


def observed_scenario(scenario_id: str) -> pd.DataFrame | None:
    rows = []
    for path in sorted(
        glob.glob(str(STUDY / scenario_id / "selected" / "replicates" / "*" / "endpoint.json"))
    ):
        record = json.loads(Path(path).read_text())
        if record.get("status") != "completed":
            continue
        frequency = record["final_allele_frequency"]
        rows.append(
            {
                "sample_af": frequency["sample_af"],
                "het": frequency.get("sample_het_diploids", 0),
                "hom": frequency.get("sample_hom_carrier_diploids", 0),
            }
        )
    return pd.DataFrame(rows) if rows else None


def class_floor_check(scenario_id: str, onset_generations: float, origin: str):
    """Heterozygotes cannot coalesce below the split (or below the mutation)."""
    floor = SPLIT_GENERATIONS if origin == "introgressed_pulse" else onset_generations
    matrices = sorted(
        glob.glob(str(STUDY / scenario_id / "selected" / "replicates" / "*" / "tmrca_matrix.npz"))
    )
    if not matrices:
        return None
    worst = 0.0
    n = 0
    for path in matrices:
        directory = Path(path).parent
        genotypes_path = directory / "focal_genotypes.npy"
        if not genotypes_path.is_file():
            continue
        genotypes = np.load(genotypes_path)
        with np.load(path) as data:
            tmrca = np.asarray(data["tmrca_generations"], dtype=float)
            focal = int(data["focal_index"])
        het = genotypes == 1
        if het.sum() < 1:
            continue
        n += 1
        # Fraction of het pairs violating the floor; theory says exactly zero.
        worst = max(worst, float(np.mean(tmrca[focal, het] < floor)))
    return {"floor_generations": floor, "replicates": n, "max_violating_fraction": worst}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=str(REPO))
    args = parser.parse_args(argv)

    trajectory = load_trajectory(args.repo_root)
    pd.set_option("display.width", 240)

    print("=" * 78)
    print("1. NEUTRAL BASELINE   P(T<t) = 1 - exp(-int_0^t dtau / 2N(tau))")
    print("=" * 78)
    predicted = neutral_p_tmrca_below(trajectory, THRESHOLDS_YEARS)
    observed = observed_neutral(THRESHOLDS_YEARS)
    table = pd.DataFrame(
        {
            "years": list(THRESHOLDS_YEARS),
            "generations": [y / GENERATION_TIME_YEARS for y in THRESHOLDS_YEARS],
            "harmonic_Ne": [
                trajectory.harmonic_mean_to(y / GENERATION_TIME_YEARS)
                for y in THRESHOLDS_YEARS
            ],
            "predicted": predicted,
        }
    )
    if observed is not None:
        table["observed"] = observed
        table["ratio"] = table["observed"] / table["predicted"].replace(0, np.nan)
    print(table.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    print()

    print("=" * 78)
    print("2. COALESCENCE DURING THE SWEEP   P ~ 1 - exp(-1/(N s p0))")
    print("=" * 78)
    rows = []
    for scenario in SCENARIOS:
        ne_at_onset = float(trajectory.size_at(scenario.onset_generations))
        if scenario.origin == "introgressed_pulse":
            p0 = ADMIXTURE_PROPORTION
            p0_label = "m = 0.025"
        else:
            p0 = 1.0 / (2.0 * ne_at_onset)
            p0_label = "1/(2N)"
        approx = 1.0 - np.exp(-1.0 / (ne_at_onset * scenario.selection_coefficient * p0))
        exact = sweep_coalescence_probability(
            ne_at_onset, scenario.selection_coefficient, p0
        )
        # The closed forms hold N fixed at its onset value; this one integrates
        # over the real trajectory, which grows tenfold across the sweep.
        with_growth = sweep_coalescence_probability_trajectory(
            trajectory, scenario.selection_coefficient, p0, scenario.onset_generations
        )
        rows.append(
            {
                "scenario": scenario.scenario_id,
                "s": scenario.selection_coefficient,
                "Ne_onset": ne_at_onset,
                "Ne_now": float(trajectory.size_at(0.0)),
                "p0": p0,
                "N*s*p0": ne_at_onset * scenario.selection_coefficient * p0,
                "P_constN_approx": approx,
                "P_constN_exact": exact,
                "P_with_growth": with_growth,
            }
        )
    print(pd.DataFrame(rows).to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    print()
    print("  Softness: P_sweep falls as s rises, so stronger selection retains more")
    print("  founding haplotypes and pushes carrier TMRCA deeper. N*s is invariant")
    print("  under Q-rescaling, so Q = 5 does not affect this prediction.")
    print()
    print("  The constant-N columns saturate at 1 and are useless as quantitative")
    print("  predictions here: they evaluate N at the onset, where EAS is at its")
    print("  smallest. P_with_growth integrates 1/(2 N(t) p(t)) over the real")
    print("  trajectory and is the column comparable to the measured fraction of")
    print("  carrier pairs coalescing after the onset.")
    print()

    # run7 measured this across a four-point s gradient on the same model, which
    # is a sharper test of the *ordering* than run8's two introgressed arms.
    softness = REPO / "sim_results_run7" / "tables" / "sweep_softness_summary.tsv"
    if softness.is_file():
        print("  Cross-check against run7's four-point gradient (same model):")
        measured = pd.read_csv(softness, sep="	")
        measured = measured[measured["origin"] == "introgressed"].sort_values(
            "selection_coefficient"
        )
        rows = []
        for _, row in measured.iterrows():
            coefficient = float(row["selection_coefficient"])
            rows.append(
                {
                    "s": coefficient,
                    "predicted_P_sweep": sweep_coalescence_probability_trajectory(
                        trajectory, coefficient, ADMIXTURE_PROPORTION, PULSE_GENERATIONS
                    ),
                    "measured_post_pulse": float(row["post_pulse_fraction"]),
                    "measured_median_TMRCA_ky": float(row["median_tmrca_ky"]),
                }
            )
        frame = pd.DataFrame(rows)
        print(frame.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
        pred_mono = bool(np.all(np.diff(frame["predicted_P_sweep"]) < 0))
        obs_mono = bool(np.all(np.diff(frame["measured_post_pulse"]) < 0))
        tail = frame[frame["s"] >= 0.003]
        tail_mono = bool(np.all(np.diff(tail["measured_post_pulse"]) < 0))
        rho = float(
            np.corrcoef(
                np.argsort(np.argsort(frame["s"])),
                np.argsort(np.argsort(frame["measured_post_pulse"])),
            )[0, 1]
        )
        print(f"  predicted strictly decreasing in s : {pred_mono}")
        print(f"  measured  strictly decreasing in s : {obs_mono}"
              f"   (rank correlation with s: {rho:+.2f})")
        print(f"  measured decreasing for s >= 0.003 : {tail_mono}")
        print()
        print("  So the prediction is only partly borne out. The trend holds from")
        print("  s = 0.003 upward, but s = 0.002 sits out of order, below s = 0.003")
        print("  rather than above it. That is the low-s subset bias documented in")
        print("  run7: at s = 0.002 only 35 of 80 replicates carry the >=10 hom")
        print("  carriers needed to be scored, and the ones that qualify are those")
        print("  that drifted up fastest -- exactly the trajectories that escape low")
        print("  frequency early and therefore coalesce least. The measurement at")
        print("  that point is conditioned on the very thing being measured.")
        print()
        print("  The absolute level is overpredicted throughout (0.99 vs 0.89 at")
        print("  s = 0.002). Theory assumes a deterministic logistic; conditioning")
        print("  on survival keeps only trajectories that left low frequency fast,")
        print("  which spend less time where carriers coalesce.")
        print()

    print("=" * 78)
    print("3. SWEEP DURATION AND FINAL FREQUENCY")
    print("=" * 78)
    rows = []
    for scenario in SCENARIOS:
        ne_at_onset = float(trajectory.size_at(scenario.onset_generations))
        if scenario.origin == "introgressed_pulse":
            p0 = ADMIXTURE_PROPORTION
            p_eff = p0
        else:
            p0 = 1.0 / (2.0 * ne_at_onset)
            p_eff = establishment_frequency(ne_at_onset, scenario.selection_coefficient)
        row = {
            "scenario": scenario.scenario_id,
            "s": scenario.selection_coefficient,
            "onset_gen": scenario.onset_generations,
            "p_eff": p_eff,
            "t_to_reach_0.5": sweep_duration_generations(
                scenario.selection_coefficient, min(p_eff, 0.49), 0.5
            ),
            "AF_pred": predicted_final_frequency(
                scenario.selection_coefficient, min(p_eff, 0.99), scenario.onset_generations
            ),
        }
        seen = observed_scenario(scenario.scenario_id)
        if seen is not None:
            row["AF_obs"] = float(seen["sample_af"].mean())
            row["n"] = int(len(seen))
            row["hom_obs"] = float(seen["hom"].mean())
        rows.append(row)
    print(pd.DataFrame(rows).to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    print()
    print("  p_eff uses the establishment size 1/(2 N alpha) for de novo origins:")
    print("  a lineage that survives drift does so by reaching ~1/alpha copies, so")
    print("  a single-copy start is not the right initial condition once we have")
    print("  conditioned on survival.")
    print()

    print("=" * 78)
    print("4. HARD FLOORS   (theory says exactly zero below these)")
    print("=" * 78)
    print(f"  human-Neanderthal split : {SPLIT_GENERATIONS:,.0f} gen "
          f"({SPLIT_GENERATIONS * GENERATION_TIME_YEARS / 1000:,.0f} kya)")
    print(f"  introgression pulse     : {PULSE_GENERATIONS:,.0f} gen "
          f"({PULSE_GENERATIONS * GENERATION_TIME_YEARS / 1000:,.0f} kya)")
    print()
    for scenario in SCENARIOS:
        check = class_floor_check(
            scenario.scenario_id, scenario.onset_generations, scenario.origin
        )
        if check is None:
            print(f"  {scenario.scenario_id:24s} (not simulated yet)")
            continue
        verdict = "OK" if check["max_violating_fraction"] == 0.0 else "VIOLATED"
        print(
            f"  {scenario.scenario_id:24s} floor {check['floor_generations']:>9,.0f} gen  "
            f"n={check['replicates']:3d}  max violating fraction "
            f"{check['max_violating_fraction']:.4f}  [{verdict}]"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
