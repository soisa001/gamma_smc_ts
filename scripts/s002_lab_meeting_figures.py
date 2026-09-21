"""Two s=.002 packs: focal carrier comparison, and all-pairs-only decoding.

Spatial carrier inference is deliberately absent: only the saved selected-allele
labels are known. No missing spatial labels are interpreted as zero ancestry.
"""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path

import lab_meeting_figures as v
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

ROOT = v.ROOT
OUT = ROOT / "eas_lab_meeting_s002_20260921"
ALL_OUT = ROOT / "eas_lab_meeting_s002_all_pairs_20260921"
POS = OUT / "analysis/positional"
T = v.T


def value(metrics, method, alpha=.05, source="decoded", field="power"):
    row = metrics[(metrics.method == method) & (metrics.alpha == alpha) & (metrics.source == source)]
    assert len(row) == 1
    return row.iloc[0][field]


def profile(key):
    base = "eas_joint_scan" if key.startswith("neutral/") else "eas_s002_saved_decoding"
    directory = ROOT/base/"profiles"/key
    path = directory/"features.npz"
    receipt = json.loads((directory/"complete.json").read_text())
    assert v.digest(path) == receipt["outputs"]["features.npz"]["sha256"]
    v.INPUTS[str(path)] = v.digest(path)
    with np.load(path, allow_pickle=False) as f:
        np.testing.assert_array_equal(f["cutoffs"], T)
        np.testing.assert_array_equal(f["sources"], ["truth", "decoded"])
        assert (f["grid_n"][:, 3] == 10000).all()
        if base == "eas_s002_saved_decoding":
            assert str(f["archaic_labels_scope"]) == "focal_only"
            np.testing.assert_array_equal(np.flatnonzero(f["carrier_labels_known"]), [500])
        markers = f["grid_markers"]
        af = np.full(len(markers), np.nan)
        af[markers >= 0] = f["site_af"][markers[markers >= 0]]
        return dict(key=key, position=f["grid"], af=af, counts=f["grid_counts"], n=f["grid_n"], valid=markers >= 0)


def design(book, all_only):
    fig, axes = v.canvas("All-pair coalescence scan" if all_only else "Focal archaic-carrier coalescence", "s = 0.002 | saved simulations, newly decoded with Gamma-SMC")
    ax = axes[0, 0]
    ax.axis("off")
    ax.text(.5, .90, "2% introgression pulse + selection starting at 50 kya\n100 selected regions + 1,000 neutral regions\n400 haplotypes; 10,000 fixed sampled pairs; isolated 10 Mb", ha="center", va="top", transform=ax.transAxes, fontsize=22, linespacing=1.65)
    formula = (r"$S_T=\frac{n_{\mathrm{pairs\ with\ decoded\ TMRCA}<T}}{10{,}000}$" if all_only else
        r"$S_T=\mathrm{frac\_recent}_T(\mathrm{ALT/ALT})\times\frac{n_{\mathrm{ALT/ALT}}}{10{,}000}$")
    ax.text(.5, .29, formula, transform=ax.transAxes, ha="center", fontsize=28)
    ax.text(.5, .02, "Posterior-mean TMRCA hard calls | fixed-coordinate neutral calibration", ha="center", transform=ax.transAxes, fontsize=17)
    v.finish(book, fig, "01_design_and_score", "Design and score", "Same saved 100 s=0.002 selected simulations; no new simulations or ancestry replays. Selection begins at the 2% pulse, 50 kya. Fixed mutation/recombination rates, 25-year generations, h=0.5. Selected power is conditional on survival and observation; all fixed replicates are retained. " +
        ("All-pair scoring uses every sampled pair regardless of archaic allele status." if all_only else "Carrier information is available at the selected allele only. The carrier score has no ALT/REF or REF/REF penalty. Spatial plots show all-pair recency only."),
        "Compare the same decoded statistic at the same coordinate in selected and neutral regions.")


def power(book, metrics, ihs, all_only):
    methods = ["all_50000"] if all_only else ["af", "mass_50000", "all_50000", "nearest_core_abs_ihs"]
    labels = {"af": "Archaic\nAF", "mass_50000": "Carrier\nmass", "all_50000": "All-pair\nrecency", "nearest_core_abs_ihs": "Positional\niHS"}
    colors = {"af": v.C["af"], "mass_50000": v.C["mass"], "all_50000": v.C["all"], "nearest_core_abs_ihs": v.C["ihs"]}
    fig, axes = v.canvas("Decoded all-pair power at s = 0.002" if all_only else "Gamma-SMC power at s = 0.002", "T = 50 kya | 100 selected replicates | error bars: 95% Wilson intervals", cols=2)
    records = []
    for ax, alpha in zip(axes[0], (.05, .01)):
        for j, method in enumerate(methods):
            df = ihs[ihs.s == .002] if method == "nearest_core_abs_ihs" else metrics[metrics.source == "decoded"]
            row = df[(df.method == method) & (df.alpha == alpha)].iloc[0]
            lo, hi = v.wilson(row.selected_called, row.selected_regions)
            ax.bar(j, row.power, color=colors[method], width=.5, zorder=3)
            ax.errorbar(j, row.power, yerr=[[row.power-lo], [hi-row.power]], fmt="none", color="#172E41", capsize=5)
            ax.text(j, hi+.025, f"{row.power:.0%}", ha="center", fontsize=22, weight="bold")
            records.append(dict(method=method, alpha=alpha, s=.002, power=row.power, selected_called=int(row.selected_called), selected_regions=int(row.selected_regions)))
        ax.set_xticks(range(len(methods)), [labels[m] for m in methods], fontsize=14)
        ax.set_title(f"Nominal p <= {alpha:g}")
        v.percent(ax, upper=1.13)
        ax.set_yticks([0, .25, .5, .75, 1])
        if all_only:
            ax.set_xlim(-1, 1)
    pd.DataFrame(records).to_csv(v.DATA/"figure04_decoded_power.csv", index=False)
    v.finish(book, fig, "04_positional_power", "Decoded focal power", "Each selected replicate is ranked against 400 neutral values at the same focal coordinate in its assigned fold; 200 held-out neutral positions per fold measure FPR. Five folds cover 1,000 neutrals. Wilson intervals describe binomial uncertainty conditional on the estimated null, not its Monte Carlo uncertainty. " +
        ("The all-pair score uses all 10,000 pairs." if all_only else "Carrier and all-pair scores use Gamma-SMC frac_recent_T; AF and iHS use genotypes. iHS is a separately calibrated endpoint."),
        f"Decoded all-pair power at T=50 kya is {value(metrics, 'all_50000'):.0%} at nominal 5%.")


def heatmaps(book, metrics, all_only, fpr=False):
    families = ["all"] if all_only else ["all", "mass"]
    fig, axes = v.canvas("Decoded neutral false-positive rates" if fpr else "TMRCA cutoffs change decoded power", "s = 0.002 | 1,000 neutral positions" if fpr else "s = 0.002 | 100 selected replicates per cell | percentages", cols=2)
    for ax, alpha in zip(axes[0], (.05, .01)):
        data = np.array([[100*value(metrics, f"{f}_{t}", alpha, field="positional_fpr" if fpr else "power") for t in T] for f in families])
        v.heatmap(ax, data, [t//1000 for t in T], ["All pairs"] if all_only else ["All pairs", "Carrier mass"], maximum=(7 if alpha == .05 else 2) if fpr else 100, fmt=".1f" if fpr else ".0f", cmap="Blues" if fpr else "YlGnBu")
        ax.set(title=f"Nominal p <= {alpha:g}", xlabel="TMRCA cutoff (kya)")
        ax.tick_params(axis="y", labelsize=13)
    fig.subplots_adjust(left=.15, wspace=.55)
    v.finish(book, fig, "06_positional_false_positive_rate" if fpr else "05_time_cutoff_power", "Positional neutral FPR" if fpr else "Decoded power across time cutoffs",
        "All scores are Gamma-SMC decoded and separately calibrated at each cutoff (5, 10, 20, 30, 40, 50 kya). No maximum across cutoffs enters the test. " +
        ("FPR is the fraction of neutral positions called, not discovery FDR. Each neutral is held out from its own calibration set; all 1,000 positions remain in the denominator. Color scales differ between nominal thresholds." if fpr else "Cells use all 100 selected replicates. Cutoff comparisons are exploratory, not per-replicate tuning."),
        "Calibration and power use the same decoded score and genomic coordinate.")


def validation(book, metrics, all_only):
    families = ["all"] if all_only else ["all", "mass"]
    fig, axes = v.canvas("Validation: true versus decoded TMRCA", "s = 0.002 | same saved genealogies, selected allele, and 10,000 pairs", cols=len(families))
    for ax, family in zip(axes[0], families):
        for src in ("truth", "decoded"):
            for alpha, style in ((.05, "-"), (.01, "--")):
                y = [value(metrics, f"{family}_{t}", alpha, src) for t in T]
                ax.plot(np.array(T)/1000, y, marker="o", color=v.C[src], ls=style, label=f"{src.capitalize()}, p <= {alpha:g}")
        v.percent(ax)
        ax.set(xlabel="TMRCA cutoff (kya)", xticks=np.array(T)/1000, title="All-pair recency" if family == "all" else "ALT/ALT carrier mass")
    v.legend(fig, ncol=4)
    v.finish(book, fig, "07_truth_vs_decoding", "Truth versus decoded validation", "This is the only figure that uses true TMRCA as a plotted reference. Each source has its own matching neutral calibration. The decoded rule thresholds posterior-mean TMRCA per pair; it does not average posterior probability mass. Matching the null calibrates errors but need not restore discrimination lost in decoding.",
        f"At 50 kya and nominal 5%, all-pair power is {value(metrics, 'all_50000', source='truth'):.0%} from truth and {value(metrics, 'all_50000'):.0%} after decoding.", category="validation")


def gain(book, metrics, paired):
    fig, axes = v.canvas("Does decoded coalescence add to archaic AF?", "s = 0.002 | nominal p <= 0.01 | paired comparison of the same 100 replicates", cols=2)
    ax = axes[0, 0]
    af = value(metrics, "af", .01)
    mass = [value(metrics, f"mass_{t}", .01) for t in T]
    ax.plot(np.array(T)/1000, [af]*len(T), color=v.C["af"], ls=":", marker="s", label="Archaic AF")
    ax.plot(np.array(T)/1000, mass, color=v.C["mass"], marker="o", label="Decoded carrier mass")
    v.percent(ax)
    ax.set(xlabel="TMRCA cutoff (kya)", xticks=np.array(T)/1000)
    p = paired[(paired.source == "decoded") & (paired.alpha == .01)].set_index("method").loc[[f"mass_{t}" for t in T]]
    ax = axes[0, 1]
    x = np.arange(len(T))
    ax.bar(x-.18, p.mass_only, width=.36, color=v.C["mass"], label="Carrier only (+)")
    ax.bar(x+.18, -p.af_only, width=.36, color=v.C["af"], label="AF only (-)")
    span = max(2, int(max(p.mass_only.max(), p.af_only.max())))+1
    ax.set(xticks=x, xticklabels=np.array(T)//1000, xlabel="TMRCA cutoff (kya)", ylabel="Discordant detections / 100", ylim=(-span, span))
    ax.axhline(0, color="#46596A", lw=1)
    h1, l1 = axes[0, 0].get_legend_handles_labels()
    h2, l2 = ax.get_legend_handles_labels()
    v.legend(fig, h1+h2, l1+l2, ncol=2)
    fig.subplots_adjust(top=.66)
    v.finish(book, fig, "08_gain_beyond_allele_frequency", "Carrier mass versus AF", "Both panels use decoded carrier mass at the focal selected allele. Negative bars are AF-only detections, positive bars carrier-only detections. Exact paired tests in the accompanying table are exploratory and unadjusted across cutoffs. No tuning of T per replicate is used.",
        f"At nominal 1%, AF detects {af:.0%}; decoded carrier mass detects {min(mass):.0%}-{max(mass):.0%} across cutoffs.")


def spatial(book, distance):
    fig, axes = v.canvas("Spatial extent of the decoded all-pair signal", "s = 0.002 | T = 50 kya | each stride uses its own matched neutral coordinate", rows=2, cols=2, gridspec_kw=dict(height_ratios=[2.3, 1]))
    fig.subplots_adjust(top=.73, bottom=.21, hspace=.28)
    shown = distance[(distance.source == "decoded") & (distance.method == "all_50000") & (distance.distance_bp.abs() <= 500000)]
    power_upper = max(.20, np.ceil(shown.power.max()*1.2*10)/10)
    for col, alpha in enumerate((.05, .01)):
        z = distance[(distance.source == "decoded") & (distance.method == "all_50000") & (distance.alpha == alpha) & (distance.distance_bp.abs() <= 500000)].sort_values("distance_bp")
        axes[0, col].plot(z.distance_bp/1000, z.power, color=v.C["all"])
        axes[1, col].plot(z.distance_bp/1000, z.positional_fpr, color=v.C["all"])
        axes[0, col].set_title(f"Nominal p <= {alpha:g}")
        v.percent(axes[0, col], "Detection (%)", upper=power_upper)
        v.percent(axes[1, col], "FPR (%)", upper=alpha*2)
        axes[1, col].set_yticks([0, alpha, alpha*2])
        axes[1, col].axhline(alpha, color="#172E41", ls="--", lw=1)
        axes[1, col].set_xlabel("Distance from selected allele (kb)")
        axes[0, col].tick_params(labelbottom=False)
        for row in range(2):
            axes[row, col].axvline(0, color="#9FAAB4", ls=":", lw=1)
            axes[row, col].set_xlim(-500, 500)
    v.finish(book, fig, "10_spatial_decay", "All-pair detection across the region", "Only all-pair scores are used. Each point is the fraction of 100 selected or 1,000 neutral regions called at that stride. Full 10-Mb results are saved; +/-500 kb is displayed without smoothing. Linked signal is not counted as a new causal target or localization error. Spatial carrier scores are unavailable from the saved labels and are omitted.", "The positional null is matched separately at every stride; unrelated peaks never set the cutoff.")


def examples(book, profiles, focal):
    selected = focal[focal.s == .002].copy()
    selected["delta"] = abs(selected.sample_af-selected.sample_af.median())
    rep = int(selected.sort_values(["delta", "replicate"]).iloc[0].replicate)
    neutrals = [p for p in profiles if p["key"].startswith("neutral/")]
    focal_values = np.array([p["counts"][1, 500, 3, -1]/10000 for p in neutrals])
    neutral_key = sorted(zip(abs(focal_values-np.median(focal_values)), [p["key"] for p in neutrals]))[0][1]
    keys = [neutral_key, f"onset50000/rep{rep:04d}"]
    by_key = {p["key"]: p for p in profiles}
    fig, axes = v.canvas("Example decoded all-pair profiles", "Examples chosen by median scores or allele frequency, independently of significance", rows=2)
    fig.subplots_adjust(top=.73, hspace=.60)
    records = []
    for ax, key in zip(axes[:, 0], keys):
        p = by_key[key]
        y = p["counts"][1, :, 3, -1]/10000
        ax.plot(p["position"]/1e6, y, color=v.C["all"], lw=1.6)
        ax.axvline(5, color="#AD7C14", ls="--", lw=1.2)
        ax.axvspan(4.95, 5.05, color="#D8A82F", alpha=.25)
        ax.set_title(("Neutral: " if key.startswith("neutral/") else "Selected s=0.002: ")+key.split("/")[-1], loc="left", fontsize=16)
        v.percent(ax, "frac_recent_T", upper=1)
        ax.set_xlim(0, 10)
        records.append(pd.DataFrame(dict(key=key, position_bp=p["position"], decoded_all_T50000=y)))
    axes[1, 0].set_xlabel("Position in scored region (Mb)")
    pd.concat(records).to_csv(v.DATA/"example_all_pair_profiles.csv", index=False)
    v.finish(book, fig, "11_example_regions", "Example all-pair profiles", f"T=50 kya. Neutral {neutral_key} has focal decoded all-pair score closest to the neutral median. Selected replicate {rep:04d} has selected AF closest to the s=0.002 median; ties use IDs. Gold shading marks a fixed 100-kb interval around 5 Mb. These examples were not selected for significant peaks and do not estimate power.", "All-pair profiles are available across the entire saved 10-Mb region.", category="backup")


def stringency(book, metrics, pooled, all_only):
    fig, axes = v.canvas("Power at stricter positional thresholds", "s = 0.002 | decoded TMRCA < 50 kya | p=0.001 uses pooled ranks")
    ax = axes[0, 0]
    methods = ["all_50000"] if all_only else ["af", "mass_50000", "all_50000"]
    x = np.arange(len(methods))
    for j, (alpha, color) in enumerate(zip((.05, .01, .001), (v.C["mass"], v.C["all"], v.C["ihs"]))):
        y = [value(pooled if alpha == .001 else metrics, m, alpha) for m in methods]
        bx = x+(j-1)*.24
        ax.bar(bx, y, width=.22, color=color, label=f"p <= {alpha:g}"+(" (pooled)" if alpha == .001 else " (5-fold)"))
        for xx, yy in zip(bx, y):
            ax.text(xx, yy+.025, f"{yy:.0%}", ha="center", fontsize=18)
    ax.set_xticks(x, ["Decoded all-pair recency"] if all_only else ["Archaic AF", "Decoded carrier mass", "Decoded all-pair recency"])
    v.percent(ax, upper=1.10)
    v.legend(fig, ncol=3)
    v.finish(book, fig, "12_threshold_stringency", "Stringency sensitivity", "Nominal 5% and 1% use 400 calibration neutrals per fold; 0.1% uses all 1,000 neutrals to rank selected regions and leave-one-out ranks for neutrals. This is a different calibration size, and only one extreme neutral rank is available. It is not precise independent validation of a 0.1% tail or genome-wide error control.", f"All-pair power at pooled p<=0.001 is {value(pooled, 'all_50000', .001):.0%}.", category="backup")


def components(book, profiles):
    selected = [p for p in profiles if p["key"].startswith("onset50000/")]
    records = []
    for p in selected:
        n = int(p["n"][500, 2])
        for t in (20000, 50000):
            count = int(p["counts"][1, 500, 2, T.index(t)])
            records.append(dict(key=p["key"], s=.002, source="decoded", T=t, af=p["af"][500], aa_pairs=n, mass=count/10000, aa_pair_fraction=n/10000, aa_frac_recent=count/n if n else np.nan))
    frame = pd.DataFrame(records)
    frame.to_csv(v.DATA/"decoded_carrier_score_components.csv", index=False)
    finite = frame.aa_pairs > 0
    np.testing.assert_allclose(frame.loc[finite, "mass"], frame.loc[finite, "aa_pair_fraction"]*frame.loc[finite, "aa_frac_recent"], atol=1e-15, rtol=0)
    fig, axes = v.canvas("Focal carrier mass combines frequency and recency", "s = 0.002 | Gamma-SMC decoded | one dot per saved selected replicate", cols=2)
    for t, color in ((20000, v.C["mass"]), (50000, v.C["all"])):
        z = frame[frame["T"] == t]
        for ax, column in zip(axes[0], ("mass", "aa_frac_recent")):
            ax.scatter(z.af, z[column], s=28, color=color, alpha=.65, linewidths=0, label=f"T={t//1000} kya")
    x = np.linspace(0, 1, 200)
    axes[0, 0].plot(x, x*x, color="#172E41", ls="--", lw=1.5, label="AF squared (approximation)")
    for ax, label in zip(axes[0], ("ALT/ALT carrier-mass score", "ALT/ALT frac_recent_T")):
        v.percent(ax, label)
        ax.set(xlabel="Selected allele frequency", xlim=(0, 1))
        ax.xaxis.set_major_formatter(v.PercentFormatter(1, decimals=0))
    v.legend(fig, ncol=3)
    v.finish(book, fig, "13_carrier_genealogies", "Focal decoded carrier-score components", "Carrier mass is exactly within-ALT/ALT frac_recent_T times the ALT/ALT fraction of the 10,000-pair panel. AF squared only approximates this sampled-pair fraction. Replicates with no sampled ALT/ALT pair have mass zero and undefined within-class recency; they remain in the power denominator. All labels are from the saved selected allele; no spatial archaic labels are inferred.", "The same pairs must both carry the allele and have recent decoded TMRCA.", category="backup")


def all_null(book, profiles):
    fig, axes = v.canvas("The all-pair score at the fixed focal coordinate", "Gamma-SMC | T = 50 kya | 100 selected versus 1,000 neutral regions")
    ax = axes[0, 0]
    records = []
    for selected, color, label in ((False, v.C["af"], "Neutral"), (True, v.C["all"], "Selected s=0.002")):
        group = [p for p in profiles if p["key"].startswith("onset50000/") == selected]
        values = np.sort([p["counts"][1, 500, 3, -1]/10000 for p in group])
        ax.step(values, np.arange(1, len(values)+1)/len(values), where="post", color=color, label=label)
        records.extend(dict(cohort=label, score=float(x)) for x in values)
    v.percent(ax, "Cumulative fraction of regions")
    ax.set(xlabel="Decoded all-pair frac_recent_T", xlim=(0, 1))
    v.legend(fig, ncol=2)
    pd.DataFrame(records).to_csv(v.DATA/"all_pair_focal_score_distribution.csv", index=False)
    v.finish(book, fig, "02_matched_position_null", "Focal all-pair score distributions", "The score uses all 10,000 sampled pairs at the same 5-Mb coordinate in each region. There is no requirement for an archaic marker and no maximum across the region. Selected simulations were conditioned on allele survival and observation when generated; no further selected or neutral filtering is performed here.", "Separation of the score distributions determines positional power.")


def write_index(all_only):
    title = "EAS s=0.002: all pairs only" if all_only else "EAS s=0.002: focal carrier comparison"
    lines = [f"# {title}", "", f"{len(v.FIGURES)} figures: PNG and vector PDF with editable text, plus one combined PDF.", "",
        "Primary TMRCA figures use Gamma-SMC decoded frac_recent_T; only figure 07 plots truth as an explicit validation reference. 100 saved selected regions, 1,000 existing decoded neutrals, 400 haplotypes and 10,000 identical sampled pairs. No simulations, ancestry recovery, or trajectory replays were run.", "",
        "All-pair spatial scores are available across the 10 Mb region. Archaic carrier classes are known only at the focal selected allele in this arm. Spatial carrier plots are omitted, and missing labels are never treated as observed absence of archaic ancestry.", "",
        "Error rates are matched-position FPR, not discovery FDR or genome-wide family-wise error. Selected power is conditional on survival and observation; fixed alleles stay in the denominator. 5-fold evaluation uses 400 neutral calibration and 200 held-out positions per fold. The p=0.001 pooled-rank sensitivity has only one extreme neutral rank. Cutoff comparisons are exploratory. Wilson intervals omit uncertainty in the estimated null.", ""]
    for f in v.FIGURES:
        lines += [f"## {f['number']:02d}. {f['title']}", "", f"[PNG](figures/{f['stem']}.png) | [Vector PDF](figures/{f['stem']}.pdf)", "", "**Takeaway:** "+f["takeaway"], "", "**Caption:** "+f["caption"], ""]
    lines += ["## Model limits", "", "Uniform mutation/recombination rates; isolated 10-Mb regions; shared strong archaic bottleneck; fixed-Ne decoder matched in neutral and selected sims. These comparisons do not test variable-Ne flow fields, empirical maps, genome-wide FDR, or independent multiple-haplotype introgression. Carrier labels use an oracle, not a finite archaic-reference/outgroup experiment. Trajectory figures 03_a/03_b remain deferred.", ""]
    (v.OUT/"README.md").write_text("\n".join(lines))
    (v.OUT/"figure_index.json").write_text(json.dumps(v.FIGURES, indent=2)+"\n")


def main():
    metrics = v.table(POS, "metrics.csv")
    assert sorted(metrics.s.unique()) == [.002]
    neutral = v.table(POS, "neutral_positions.csv")
    paired = v.table(POS, "paired_af_comparison.csv", manifest="supplement_manifest.json")
    pooled = v.table(POS, "pooled_rank_metrics.csv", manifest="supplement_manifest.json")
    focal = v.table(ROOT/"eas_lab_meeting_20260921/analysis", "focal_truth.csv")
    ihs = v.table(ROOT/"eas_ihs_h400", "metrics.csv")
    ihs_regions = v.table(ROOT/"eas_ihs_h400", "regions.csv")
    distance = v.table(OUT/"analysis/distance", "metrics.csv")
    regions = v.table(ROOT/"eas_allele_class_ablation", "regions.csv")
    methods = json.loads(v.checked(ROOT/"eas_allele_class_ablation", "methods.json").read_text())
    with np.load(v.checked(ROOT/"eas_allele_class_ablation", "scores.npz"), allow_pickle=False) as f:
        region_max = f["scores"][0, 0, 0, regions.key.str.startswith("neutral/"), next(i for i, m in enumerate(methods) if m["name"] == "af_r1")]
    with ThreadPoolExecutor(max_workers=4) as pool:
        profiles = list(pool.map(profile, regions.key.tolist()))
    cfg = json.loads((ROOT/"eas_q02_h400/manifest.json").read_text())["config"]
    for all_only in (False, True):
        v.OUT = ALL_OUT if all_only else OUT
        v.FIG = v.OUT/"figures"
        v.DATA = v.OUT/"figure_data"
        v.FIG.mkdir(parents=True, exist_ok=True)
        v.DATA.mkdir(parents=True, exist_ok=True)
        v.FIGURES.clear()
        v.LAYOUT.clear()
        for name, df in (("positional_metrics", metrics), ("pooled_rank_metrics", pooled), ("distance_metrics", distance)):
            data = df[df.method.str.startswith("all_")] if all_only else df
            data.to_csv(v.DATA/f"{name}.csv", index=False)
        if not all_only:
            for name, df in (("focal_truth", focal), ("ihs_metrics", ihs), ("ihs_regions", ihs_regions), ("paired_af_comparison", paired), ("neutral_focal_af", neutral)):
                df.to_csv(v.DATA/f"{name}.csv", index=False)
        with PdfPages(v.OUT/"EAS_lab_meeting_figures.pdf", metadata={"Title": "EAS s=0.002 - "+("all pairs only" if all_only else "focal carrier comparison"), "CreationDate": v.DATE, "ModDate": v.DATE}) as book:
            design(book, all_only)
            if all_only:
                all_null(book, profiles)
            else:
                v.local_null(book, neutral, region_max)
                v.allele_frequency(book, focal, neutral, focus_s=.002)
            power(book, metrics, ihs, all_only)
            heatmaps(book, metrics, all_only)
            heatmaps(book, metrics, all_only, fpr=True)
            validation(book, metrics, all_only)
            if not all_only:
                gain(book, metrics, paired)
                v.ihs_comparison(book, ihs, ihs_regions, focus_s=.002)
            spatial(book, distance)
            examples(book, profiles, focal)
            stringency(book, metrics, pooled, all_only)
            if not all_only:
                components(book, profiles)
            v.demography(book, cfg)
            if not all_only:
                v.pair_classes(book, profiles, selection_s=.002)
        write_index(all_only)
        policy = dict(primary_tmrca="Gamma-SMC decoded posterior-mean hard calls; frac_recent_T", truth_comparison_figures=["07_truth_vs_decoding"],
            decoded_selection_coefficients=[.002], decoded_selected_regions=100, decoded_neutral_regions=1000,
            all_pairs_only=all_only, carrier_label_scope="focal selected allele only", spatial_carrier_plots=False,
            simulation_replays=False, new_simulations=False, existing_selected_trees_newly_decoded=True)
        (v.OUT/"plot_sources.json").write_text(json.dumps(policy, indent=2)+"\n")
        (v.OUT/"layout_checks.json").write_text(json.dumps(v.LAYOUT, indent=2)+"\n")
        (v.OUT/"figure_provenance.json").write_text(json.dumps(dict(source_sha256=v.digest(Path(__file__)), shared_figure_source_sha256=v.digest(Path(v.__file__)), inputs=v.INPUTS, config=cfg, plot_source_policy=policy, seed_policy="Original simulation and pair seeds retained; no resampling or jitter.", matplotlib=v.matplotlib.__version__, numpy=np.__version__), indent=2)+"\n")
        print(json.dumps(dict(figures=len(v.FIGURES), destination=str(v.OUT))), flush=True)


if __name__ == "__main__":
    main()
