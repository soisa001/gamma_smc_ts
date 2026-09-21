"""Render the s=.005 all-pair-only pack from checked saved decoding and analyses."""
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path

import s002_lab_meeting_figures as plots
import lab_meeting_figures as v
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

ROOT = v.ROOT
OUT = ROOT/"eas_lab_meeting_s005_all_pairs_20260921"
POS = ROOT/"eas_lab_meeting_20260921/analysis/positional"
S = .005


def main():
    v.OUT, v.FIG, v.DATA = OUT, OUT/"figures", OUT/"figure_data"
    v.FIG.mkdir(parents=True, exist_ok=True)
    v.DATA.mkdir(parents=True, exist_ok=True)
    v.FIGURES.clear()
    v.LAYOUT.clear()
    metrics = v.table(POS, "metrics.csv")
    metrics = metrics[(metrics.s == S) & metrics.method.str.startswith("all_")].copy()
    pooled = v.table(POS, "pooled_rank_metrics.csv", manifest="supplement_manifest.json")
    pooled = pooled[(pooled.s == S) & pooled.method.str.startswith("all_")].copy()
    distance = v.table(ROOT/"eas_positional_distance_h400", "metrics.csv")
    distance = distance[distance.method.str.startswith("all_")].copy()
    focal = v.table(ROOT/"eas_lab_meeting_20260921/analysis", "focal_truth.csv")
    regions = v.table(ROOT/"eas_allele_class_ablation", "regions.csv")
    with ThreadPoolExecutor(max_workers=4) as pool:
        profiles = list(pool.map(v.profile, regions.key.tolist()))
    assert len(profiles) == 1100 and len(metrics) == 24 and len(pooled) == 12
    assert (metrics.selected_regions == 100).all() and (metrics.neutral_positions == 1000).all()
    source_manifest = ROOT/"eas_joint_scan/profile_manifest.json"
    source = json.loads(source_manifest.read_text())
    v.INPUTS[str(source_manifest)] = v.digest(source_manifest)
    items = {item["key"]: item for item in source["items"]}
    selected = regions.key.str.startswith("onset50000/").to_numpy()
    assert selected.sum() == 100 and (~selected).sum() == 1000
    for key in regions.key:
        item = items[key]
        if key.startswith("onset50000/"):
            assert item["s"] == S and item["onset_years"] == 50000
        receipt = json.loads((ROOT/"eas_joint_scan/profiles"/key/"complete.json").read_text())
        assert receipt["identity"]["simulation"] == item["receipt_sha256"]
    # Independently reproduce focal metrics from raw integer profile counts.
    folds = regions.fold.to_numpy()
    scores = np.array([p["counts"][:, 500, 3, :]/10000 for p in profiles])
    for src, source_name in enumerate(("truth", "decoded")):
        for t_index, t in enumerate(v.T):
            pvalues = np.full(len(profiles), np.nan)
            for fold in range(5):
                calibration = ~selected & np.isin(folds, [(fold+3)%5, (fold+4)%5])
                test = folds == fold
                assert calibration.sum() == 400
                pvalues[test] = (1+(scores[calibration, src, t_index, None] >= scores[test, src, t_index]).sum(axis=0))/401
            for alpha in (.01, .05):
                row = metrics[(metrics.source == source_name) & (metrics.method == f"all_{t}") & (metrics.alpha == alpha)].iloc[0]
                assert int((pvalues[selected] <= alpha).sum()) == row.selected_called
                assert int((pvalues[~selected] <= alpha).sum()) == row.neutral_called
                at_center = distance[(distance.source == source_name) & (distance.method == f"all_{t}") & (distance.alpha == alpha) & (distance.distance_bp == 0)].iloc[0]
                assert at_center.selected_called == row.selected_called and at_center.neutral_called == row.neutral_called
    for name, frame in (("positional_metrics", metrics), ("pooled_rank_metrics", pooled), ("distance_metrics", distance)):
        frame.to_csv(v.DATA/f"{name}.csv", index=False)
    cfg = source["config"]
    with PdfPages(OUT/"EAS_lab_meeting_figures.pdf", metadata={"Title": "EAS s=0.005 - all pairs only", "CreationDate": v.DATE, "ModDate": v.DATE}) as book:
        plots.design(book, all_only=True, selection_s=S, newly_decoded=False)
        plots.all_null(book, profiles, selection_s=S)
        plots.power(book, metrics, None, all_only=True, selection_s=S)
        plots.heatmaps(book, metrics, all_only=True, selection_s=S)
        plots.heatmaps(book, metrics, all_only=True, fpr=True, selection_s=S)
        plots.validation(book, metrics, all_only=True, selection_s=S)
        plots.spatial(book, distance, selection_s=S)
        plots.examples(book, profiles, focal, selection_s=S)
        plots.stringency(book, metrics, pooled, all_only=True, selection_s=S)
        v.demography(book, cfg)
    assert len(v.FIGURES) == 10
    policy = dict(primary_tmrca="Gamma-SMC decoded posterior-mean hard calls; frac_recent_T",
        truth_comparison_figures=["07_truth_vs_decoding"], decoded_selection_coefficients=[S],
        decoded_selected_regions=100, decoded_neutral_regions=1000, all_pairs_only=True,
        simulations_started=False, decoding_started=False, simulation_replays=False)
    lines = ["# EAS s=0.005: all pairs only", "", "10 figures in the same layout as the s=0.002 all-pair pack. Each has a PNG and a vector PDF with editable text; the combined PDF has one figure per page.", "",
        "All primary TMRCA figures use saved Gamma-SMC decoded frac_recent_T across all 10,000 sampled pairs, without carrier conditioning. Figure 07 alone plots true TMRCA as an explicit validation reference. The cohort contains 100 selected and 1,000 neutral regions, with 400 haplotypes per region. No new simulations, decoding or ancestry replays were run.", "",
        "Selection starts at the 2% pulse, 50 kya. The primary null is matched to the same coordinate. Five-fold evaluation uses 400 calibration and 200 held-out neutrals per fold; all fixed selected replicates remain included. Power is conditional on original allele survival/observation. Positional FPR is not discovery FDR or genome-wide family-wise error.", ""]
    for item in v.FIGURES:
        lines += [f"## {item['number']:02d}. {item['title']}", "", f"[PNG](figures/{item['stem']}.png) | [Vector PDF](figures/{item['stem']}.pdf)", "", "**Takeaway:** "+item["takeaway"], "", "**Caption:** "+item["caption"], ""]
    lines += ["## Interpretation limits", "", "Uniform mutation/recombination, isolated 10-Mb regions, a shared strong archaic bottleneck and a fixed-Ne decoder. Cutoffs are compared separately, without per-replicate tuning. Wilson intervals omit uncertainty in the estimated neutral null. Pooled p=0.001 is a sensitivity analysis with only one extreme neutral rank. The data do not establish empirical genome-wide FDR control or causal localization.", ""]
    (OUT/"README.md").write_text("\n".join(lines))
    for name, data in (("figure_index.json", v.FIGURES), ("layout_checks.json", v.LAYOUT), ("plot_sources.json", policy)):
        (OUT/name).write_text(json.dumps(data, indent=2)+"\n")
    (OUT/"figure_provenance.json").write_text(json.dumps(dict(source_sha256=v.digest(Path(__file__)),
        shared_plot_source_sha256=v.digest(Path(plots.__file__)), shared_style_source_sha256=v.digest(Path(v.__file__)),
        inputs=v.INPUTS, config=cfg, plot_source_policy=policy,
        seed_policy="No new random draws; original simulation and pair seeds retained.",
        matplotlib=v.matplotlib.__version__, numpy=np.__version__), indent=2)+"\n")
    (OUT/"analysis_audit").mkdir(exist_ok=True)
    audit = dict(status="passed", selection_coefficient=S, selected_regions=100, neutral_regions=1000,
        verified_profile_count=1100, focal_metrics_recomputed_from_integer_counts=True,
        focal_spatial_center_matches=True, only_all_pair_methods=True, new_simulations=False, new_decoding=False)
    (OUT/"analysis_audit/source_validation.json").write_text(json.dumps(audit, indent=2)+"\n")
    print(json.dumps(dict(figures=len(v.FIGURES), destination=str(OUT), audit=audit)), flush=True)


if __name__ == "__main__":
    main()
