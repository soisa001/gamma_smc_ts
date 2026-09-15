"""Paired regional comparison of carrier mass and the REF/REF contrast.

Uses existing audited predictions only; no simulations, decoding, or refitting.
The conditional McNemar p-values are exploratory because cross-validation
shares training/calibration regions and this cohort informed method development.
"""

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path

from scipy.stats import binomtest


NEW = "mass_g0_r1_T50000"
OLD = "weighted_contrast_g0_r1_T50000"
INPUTS = (
    "predictions.csv.gz", "target70_predictions.csv.gz", "metrics.csv",
    "target70_metrics.csv", "regions.csv", "provenance.json",
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as stream:
        yield from csv.DictReader(stream)


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def paired_counts(pairs):
    both = sum(new and old for new, old in pairs)
    new_only = sum(new and not old for new, old in pairs)
    old_only = sum(old and not new for new, old in pairs)
    discordant = new_only + old_only
    p = float(binomtest(new_only, discordant).pvalue) if discordant else 1.0
    n = len(pairs)
    return dict(
        n=n, both=both, mass_only=new_only, contrast_only=old_only,
        neither=n - both - discordant, mass_called=both + new_only,
        contrast_called=both + old_only,
        mass_rate=(both + new_only) / n,
        contrast_rate=(both + old_only) / n,
        mass_minus_contrast_pp=100 * (new_only - old_only) / n,
        mcnemar_p_nominal=p,
    )


def compare(source, out):
    manifest = json.loads((source / "artifact_manifest.json").read_text())
    for name in INPUTS:
        assert (source / name).stat().st_size == manifest[name]["bytes"], name
        assert digest(source / name) == manifest[name]["sha256"], name
    audit = json.loads((source / "audit.json").read_text())
    assert audit["status"] == "passed" and audit["no_later_onset"]
    config = json.loads((source / "provenance.json").read_text())["config"]
    assert config["statistic"] == "frac_recent_T"
    assert config["selection_onset_years"] == config["pulse_years"] == 50000
    assert config["selection_coefficients"] == [0.005]
    regions = {r["key"]: (int(r["onset_years"]), int(r["fold"]))
               for r in read_csv(source / "regions.csv")}
    assert len(regions) == 1100
    results, fold_results, matched_rows, gate_checks = [], [], [], []
    for label, prediction_file, metric_file, alpha in (
        ("regional_p_0.05", "predictions.csv.gz", "metrics.csv", 0.05),
        ("training_power_target_0.70", "target70_predictions.csv.gz", "target70_metrics.csv", 0.7),
    ):
        methods = (NEW, OLD, NEW.replace("_g0_", "_g50_"), OLD.replace("_g0_", "_g50_"))
        calls = {}
        for r in read_csv(source / prediction_file):
            if r["method"] not in methods or float(r["alpha"]) != alpha:
                continue
            assert r["called"] in ("True", "False")
            assert regions[r["key"]] == (int(r["onset_years"]), int(r["fold"]))
            key = (r["scheme"], r["source"], r["key"], r["method"])
            assert key not in calls
            calls[key] = r["called"] == "True"
        assert len(calls) == 17600
        metrics = {(r["scheme"], r["source"], r["method"]): r
                   for r in read_csv(source / metric_file)
                   if r["method"] in (NEW, OLD) and float(r["alpha"]) == alpha}
        assert len(metrics) == 8
        for scheme in ("grid10kb", "archaic_sites"):
            for tmrca in ("decoded", "truth"):
                for method in (NEW, OLD):
                    changes = sum(calls[scheme, tmrca, key, method] !=
                                  calls[scheme, tmrca, key, method.replace("_g0_", "_g50_")]
                                  for key in regions)
                    gate_checks.append(dict(operating_point=label, scheme=scheme,
                                            source=tmrca, method=method, changed_calls=changes))
                for onset, cohort in ((50000, "selected"), (0, "neutral")):
                    keys = [key for key, (year, _) in regions.items() if year == onset]
                    base = dict(operating_point=label, scheme=scheme, source=tmrca, cohort=cohort)
                    pairs = [(calls[scheme, tmrca, key, NEW], calls[scheme, tmrca, key, OLD])
                             for key in keys]
                    row = dict(**base, **paired_counts(pairs))
                    assert row["n"] == (100 if onset else 1000)
                    for method, column in ((NEW, "mass_called"), (OLD, "contrast_called")):
                        expected = metrics[scheme, tmrca, method]
                        assert row[column] == int(expected[f"{cohort}_called"])
                        assert row["n"] == int(expected[f"{cohort}_n"])
                    results.append(row)
                    for fold in range(5):
                        fp = [pair for key, pair in zip(keys, pairs) if regions[key][1] == fold]
                        fold_results.append(dict(**base, fold=fold, **paired_counts(fp)))
                    for key, (new, old) in zip(keys, pairs):
                        matched_rows.append(dict(**base, key=key, fold=regions[key][1],
                                                 mass_called=new, contrast_called=old))
    assert len(results) == 16 and len(matched_rows) == 8800
    # Holm family = all 16 comparisons displayed here, not all historical analyses.
    running = 0.0
    for rank, row in enumerate(sorted(results, key=lambda r: r["mcnemar_p_nominal"])):
        running = max(running, min(1.0, (len(results) - rank) * row["mcnemar_p_nominal"]))
        row["mcnemar_p_holm_16"] = running
    out.mkdir(parents=True, exist_ok=True)
    write_csv(out / "paired_comparisons.csv", results)
    write_csv(out / "fold_comparisons.csv", fold_results)
    write_csv(out / "paired_calls.csv", matched_rows)
    write_csv(out / "gate_sensitivity.csv", gate_checks)
    lines = ["# Carrier mass versus REF/REF subtraction", "",
             "Mass: w_AA * frac_recent_T_AA. Contrast: w_AA * max(frac_recent_T_AA - frac_recent_T_RR, 0).",
             "Both omit ALT/REF and the all-pair gate. T=50 kya; immediate-onset EAS s=0.005.",
             "The unit of comparison is an independently simulated 10 Mb region, not a haplotype pair.", "",
             "| Operating point | Positions | TMRCA | Cohort | Mass calls | Contrast calls | Mass only | Contrast only | Nominal paired p | Holm p (16) |",
             "|---|---|---|---|---:|---:|---:|---:|---:|---:|"]
    for r in results:
        lines.append(f"| {r['operating_point']} | {r['scheme']} | {r['source']} | {r['cohort']} | "
                     f"{r['mass_called']}/{r['n']} | {r['contrast_called']}/{r['n']} | "
                     f"{r['mass_only']} | {r['contrast_only']} | {r['mcnemar_p_nominal']:.6g} | "
                     f"{r['mcnemar_p_holm_16']:.6g} |")
    lines.extend(["", "Two-sided exact McNemar calculation: binomial test of mass-only versus contrast-only calls.",
                  "These are conditional, exploratory p-values on the saved fold decisions. Shared training/calibration",
                  "regions create dependence among cross-validated predictions; the calculation does not propagate",
                  "threshold-fitting uncertainty. Holm adjustment covers only these 16 reported comparisons, not",
                  "earlier score/cutoff development. Independent new-seed validation is needed for confirmation.", "",
                  "The 70% operating point trains each score's threshold separately toward 70% training power;",
                  "held-out power is not forced to match. It is not a p=0.70 significance threshold.",
                  "Neutral call fraction is regional false-positive rate, not prevalence-dependent FDR.", ""])
    (out / "report.md").write_text("\n".join(lines))
    provenance = dict(source_directory=str(source), config=config, mass_method=NEW, contrast_method=OLD,
                      script_sha256=digest(Path(__file__)), source_files={n: manifest[n] for n in INPUTS},
                      source_audit_sha256=digest(source / "audit.json"),
                      no_simulation_decoding_or_refitting=True, rng_used=False,
                      test="two-sided conditional exact McNemar; exploratory", holm_family_size=16,
                      primary_comparison="regional_p_0.05 / grid10kb / decoded")
    (out / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    names = ("paired_comparisons.csv", "fold_comparisons.csv", "paired_calls.csv", "gate_sensitivity.csv",
             "report.md", "provenance.json")
    (out / "artifact_manifest.json").write_text(json.dumps(
        {n: dict(bytes=(out / n).stat().st_size, sha256=digest(out / n)) for n in names}, indent=2) + "\n")
    print(json.dumps(dict(status="verified", comparisons=len(results), matched_region_calls=len(matched_rows),
                          historical_gate_changed_calls=sum(r["changed_calls"] for r in gate_checks),
                          out=str(out)), indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=Path("/mnt/d/phase2simselection/sim/eas_allele_class_ablation"))
    parser.add_argument("--out", type=Path, default=Path("/mnt/d/phase2simselection/sim/eas_alt_alt_ref_comparison"))
    args = parser.parse_args()
    compare(args.source, args.out)
