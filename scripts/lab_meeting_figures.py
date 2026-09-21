"""Build a verified, slide-sized figure collection from the completed EAS array."""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path

os.environ["MPLBACKEND"] = "Agg"
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "1"
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import FancyBboxPatch
from matplotlib.ticker import PercentFormatter
import numpy as np
import pandas as pd

from gamma_smc_aou.fresh_power import model

ROOT = Path("/mnt/d/phase2simselection/sim")
OUT = ROOT / "eas_lab_meeting_20260921"
POS = OUT / "analysis/positional"
FIG = OUT / "figures"
DATA = OUT / "figure_data"
S = np.arange(1, 11)/1000
T = (5000, 10000, 20000, 30000, 40000, 50000)
C = dict(af="#556575", mass="#087F8C", all="#CA6424", ihs="#8256A7", truth="#087F8C", decoded="#C55B42")
INPUTS, FIGURES, LAYOUT = {}, [], []
DATE = datetime(2026, 9, 21, tzinfo=timezone.utc)
matplotlib.rcParams.update({"font.family":"DejaVu Sans", "font.size":17, "axes.titlesize":20,
    "axes.labelsize":18, "xtick.labelsize":15, "ytick.labelsize":15,
    "legend.fontsize":15, "axes.spines.top":False, "axes.spines.right":False,
    "pdf.fonttype":42, "ps.fonttype":42, "axes.unicode_minus":False,
    "savefig.facecolor":"white", "figure.facecolor":"white", "lines.linewidth":2.5})


def digest(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def checked(directory, name, manifest="artifact_manifest.json"):
    path = directory/name
    expected = json.loads((directory/manifest).read_text())[name]["sha256"]
    actual = digest(path)
    if expected != actual:
        raise ValueError(f"Corrupt figure input: {path}")
    INPUTS[str(path)] = actual
    return path


def table(directory, name, **kwargs):
    return pd.read_csv(checked(directory, name, **kwargs), float_precision="round_trip")


def canvas(title, subtitle, rows=1, cols=1, **kwargs):
    fig, axes = plt.subplots(rows, cols, figsize=(13.3333, 7.5), squeeze=False, **kwargs)
    fig.subplots_adjust(left=.085, right=.965, bottom=.21, top=.70, hspace=.60, wspace=.33)
    fig.text(.04, .948, title, fontsize=25, weight="bold", ha="left", va="top", color="#142E42")
    fig.text(.04, .884, subtitle, fontsize=15, ha="left", va="top", color="#46596A")
    return fig, axes


def legend(fig, handles=None, labels=None, ncol=3):
    if handles is None:
        handles, labels = fig.axes[0].get_legend_handles_labels()
    item=fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(.52, .837),
               ncol=ncol, frameon=False, columnspacing=1.7, handlelength=2.5)
    for handle in item.legend_handles:
        if isinstance(handle, matplotlib.collections.PathCollection):
            handle.set_sizes([45]);handle.set_alpha(1)


def percent(ax, label="Power (%)", upper=1.04):
    ax.set_ylim(0, upper)
    ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
    ax.set_ylabel(label)
    ax.grid(axis="y", color="#E5E9EE", zorder=0)


def scoeff(ax):
    ax.set_xticks(S[::2], [f"{s:.3f}" for s in S[::2]])
    ax.set_xlim(.0007, .0103)
    ax.set_xlabel("Selection coefficient (s)")


def rows(data, source="truth", method="mass_50000", alpha=.05):
    return data[(data.source == source) & (data.method == method) & (data.alpha == alpha)].sort_values("s")


def wilson(k, n):
    k, n = np.asarray(k), np.asarray(n)
    z = 1.959963984540054
    center = (k/n+z*z/(2*n))/(1+z*z/n)
    radius = z*np.sqrt(k/n*(1-k/n)/n+z*z/(4*n*n))/(1+z*z/n)
    return np.maximum(center-radius, 0), np.minimum(center+radius, 1)


def curve(ax, frame, label, color, style="-", band=False):
    is_af=color==C["af"]
    ax.plot(frame.s, frame.power, marker="s" if is_af else "o", ms=5,
        ls=":" if is_af and style=="-" else style, color=color, label=label,
        zorder=5 if is_af else 3, markerfacecolor="white" if is_af else color)
    if band:
        lo, hi = wilson(frame.selected_called, frame.selected_regions)
        ax.fill_between(frame.s, lo, hi, color=color, alpha=.09, linewidth=0)


def finish(book, fig, stem, title, caption, takeaway, category="main"):
    number = len(FIGURES)+1
    fig.text(.04, .040, "EAS | 2% pulse + selection at 50 kya | isolated 10 Mb | 400 haplotypes", fontsize=11, color="#46596A")
    fig.text(.96, .040, f"{number:02d}  /  21 Sep 2026", fontsize=11, color="#46596A", ha="right")
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    hidden_ticks = set()
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            lo, hi = sorted(axis.get_view_interval())
            for tick in axis.get_major_ticks()+axis.get_minor_ticks():
                if not lo-1e-10 <= tick.get_loc() <= hi+1e-10:
                    hidden_ticks.update((id(tick.label1), id(tick.label2)))
    outside = []
    for text in fig.findobj(matplotlib.text.Text):
        if not text.get_visible() or not text.get_text() or id(text) in hidden_ticks:
            continue
        box = text.get_window_extent(renderer)
        if box.width > 0 and box.height > 0 and (box.x0 < -2 or box.y0 < -2 or box.x1 > fig.bbox.width+2 or box.y1 > fig.bbox.height+2):
            outside.append(text.get_text())
    LAYOUT.append(dict(figure=stem, text_outside_canvas=outside))
    if outside:
        raise ValueError(f"Text outside canvas in {stem}: {outside}")
    # Matplotlib paths/markers remain vector; embedded PDF text uses TrueType.
    for artist in fig.findobj():
        if artist.get_rasterized():
            artist.set_rasterized(False)
    fig.savefig(FIG/f"{stem}.png", dpi=240)
    fig.savefig(FIG/f"{stem}.pdf",metadata={"Title":title,"CreationDate":DATE,"ModDate":DATE})
    book.attach_note(caption)
    book.savefig(fig)
    plt.close(fig)
    FIGURES.append(dict(number=number, stem=stem, title=title, category=category, caption=caption, takeaway=takeaway))


def design(book):
    fig, aa = canvas("An archaic-carrier coalescence scan", "Completed simulation design and the primary score")
    ax = aa[0,0]
    ax.axis("off")
    boxes = [(.02,.58,.27,.40,"2% introgression pulse\nSelection starts at 50 kya\ns = 0.001 to 0.010"),
             (.365,.58,.27,.40,"1,000 neutral regions\n100 selected per s\nAllele survival required"),
             (.71,.58,.27,.40,"400 haplotypes\n10,000 fixed sampled pairs\n11 Mb simulated; 10 Mb scored")]
    for x,y,w,h,label in boxes:
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=.014", transform=ax.transAxes,
            facecolor="#EAF3F6", edgecolor="#BFD7DE"))
        ax.text(x+w/2,y+h/2,label,transform=ax.transAxes,ha="center",va="center",fontsize=16,linespacing=1.5)
    for x in (.29,.65):
        ax.annotate("", xy=(x+.06,.78), xytext=(x,.78), xycoords="axes fraction", arrowprops=dict(arrowstyle="->",color="#46596A",lw=2))
    ax.text(.5,.36,r"$S_T=\mathrm{frac\_recent}_T(\mathrm{ALT/ALT})\times\frac{n_{\mathrm{ALT/ALT}}}{n_{\mathrm{all\ pairs}}}=\frac{n_{\mathrm{recent\ ALT/ALT},T}}{n_{\mathrm{all\ pairs}}}$",ha="center",va="center",transform=ax.transAxes,fontsize=23)
    ax.text(.5,.12,"Eligible archaic ALT: oracle-specific origin after the archaic-modern split, before the pulse.\nTruth uses local TMRCA; decoding uses posterior-mean TMRCA hard calls.",transform=ax.transAxes,ha="center",va="center",fontsize=17,linespacing=1.6)
    finish(book,fig,"01_design_and_score","Design and carrier score",
        "EAS; generation time 25 years; mutation rate 1.25e-8 and recombination rate 1e-8 per base per generation; dominance 0.5. Selected alleles must survive and be observed; fixation is retained. ALT/ALT denotes two haplotypes, not a diploid homozygote. The score uses frac_recent_T and no ALT/REF penalty or REF/REF subtraction.",
        "Require the same haplotype pairs to carry the archaic allele and coalesce recently.")


def local_null(book, neutral, region_max):
    fig, axes = canvas("The null must match the position being tested", "Same 1,000 neutral simulations; different spatial questions", cols=2)
    ax=axes[0,0]
    for values,label,color in ((neutral.af,"Focal position",C["mass"]),(region_max,"Maximum in 10 Mb",C["af"])):
        values=np.sort(values)
        ax.step(values,np.arange(1,len(values)+1)/len(values),where="post",label=label,color=color)
    ax.axvline(.18,color="#B84945",ls="--",lw=1.7)
    ax.text(.20,.25,"18% AF",color="#B84945",fontsize=16)
    percent(ax,"Cumulative fraction of neutrals")
    ax.set(xlabel="Archaic marker frequency",xlim=(-.015,.91))
    ax.xaxis.set_major_formatter(PercentFormatter(1,decimals=0))
    ax=axes[0,1]
    values=[float((neutral.af >= .18).mean()),float((region_max >= .18).mean())]
    bars=ax.bar([0,1],values,color=[C["mass"],C["af"]],width=.58,zorder=3)
    for b,v in zip(bars,values):
        ax.text(b.get_x()+b.get_width()/2,v+.03,f"{v:.1%}",ha="center",weight="bold",fontsize=24)
    ax.set_xticks([0,1],["At the focal\nposition","Somewhere\nin 10 Mb"])
    percent(ax,"Neutral fraction with AF >=18%",upper=1.05)
    legend(fig,ncol=2)
    finish(book,fig,"02_matched_position_null","Positional versus regional null",
        "The focal score assigns the nearest eligible archaic marker within 5 kb of the fixed focal stride. All 1,000 neutral positions are retained; 169 have an assigned marker. Marker-free positions score zero. The whole-10-Mb maximum is shown only to explain the earlier calibration mismatch. A pre-specified gene requires calibration of the same summary in a matching neutral gene window.",
        "Only 4.3% of neutral focal positions reach 18% AF; unrelated peaks elsewhere do not enter the local test.")


def allele_frequency(book, focal, neutral):
    fig, axes=canvas("Selection shifts local archaic-allele frequency", "Completed cohort: 100 selected replicates per coefficient; 1,000 neutral focal positions")
    ax=axes[0,0]
    fig.subplots_adjust(bottom=.26)
    groups=[neutral.af.to_numpy()]+[focal.loc[focal.s == s,"sample_af"].to_numpy() for s in S]
    bp=ax.boxplot(groups,positions=np.arange(11),widths=.60,whis=(5,95),showfliers=False,patch_artist=True,
        medianprops=dict(color="white",lw=2.2))
    for j,box in enumerate(bp["boxes"]):box.set(facecolor=C["af"] if j==0 else C["mass"],alpha=.92)
    ax.plot(np.arange(11),[g.mean() for g in groups],"D",color="#172E41",ms=5,label="Mean")
    ax.axhline(.02,color="#B84945",ls="--",label="Nominal 2% pulse",lw=1.6)
    ax.set_xticks(np.arange(11),["Neutral"]+[f"{s:.3f}" for s in S],rotation=25,ha="right")
    percent(ax,"Archaic allele frequency")
    ax.set_xlabel("Selection coefficient (s)")
    legend(fig,ncol=2)
    finish(book,fig,"03_allele_frequency_distribution","Allele frequency by selection strength",
        "Boxes show interquartile range and median; whiskers show the 5th and 95th percentiles; diamonds show means. Selected data are conditional on survival and sample observation. Neutral data are unconditional pre-specified positions, with absent nearby markers scored zero. These are not frequency-matched cohorts.",
        "Weak selection produces a broad AF distribution; its mean does not imply every replicate is detectable.")


def power(book, metrics, ihs):
    fig, axes=canvas("Local calibration recovers partial-sweep power", "T = 50 kya for TMRCA scores | shaded bands: 95% Wilson intervals conditional on the fitted null",cols=2)
    for ax,alpha in zip(axes[0],(.05,.01)):
        for method,label,key in (("af","Archaic AF","af"),("mass_50000","ALT/ALT carrier mass","mass"),("all_50000","All-pair recency","all")):
            curve(ax,rows(metrics,method=method,alpha=alpha),label,C[key],band=True)
        z=ihs[(ihs.method == "nearest_core_abs_ihs") & (ihs.alpha == alpha)].sort_values("s")
        curve(ax,z,"Positional iHS",C["ihs"],style="--",band=True)
        scoeff(ax);percent(ax);ax.set_title(f"Nominal p <= {alpha:g}",pad=14)
    legend(fig,ncol=4)
    finish(book,fig,"04_positional_power","Power across selection coefficients",
        "All TMRCA curves use true TMRCA. Each position is compared with 400 neutral calibration positions in one of five folds; 200 neutral positions per fold are held out. All 1,000 selected and 1,000 neutrals are included. iHS uses the nearest scorable core within 5 kb. Wilson bands do not include uncertainty from estimating the neutral reference or cross-fold dependence.",
        "Carrier-based local tests are much more sensitive than all-pair recency at weak selection.")


def heatmap(ax, values, xticklabels, yticklabels, maximum=100, fmt=".0f", cmap="YlGnBu"):
    # Explicit vector cells: imshow embeds a bitmap even with rasterized=False.
    nr, nc = values.shape
    im=ax.pcolormesh(np.arange(nc+1)-.5,np.arange(nr+1)-.5,values,
        vmin=0,vmax=maximum,cmap=cmap,shading="flat",rasterized=False)
    ax.set_xlim(-.5,nc-.5);ax.set_ylim(nr-.5,-.5)
    ax.set_xticks(range(len(xticklabels)),xticklabels)
    ax.set_yticks(range(len(yticklabels)),yticklabels)
    for (i,j),v in np.ndenumerate(values):
        ax.text(j,i,format(v,fmt),ha="center",va="center",fontsize=14.5,color="white" if v>maximum*.62 else "#142E42")
    return im


def time_power(book, metrics):
    fig,axes=canvas("Recency cutoffs change power", "True TMRCA | p <= 0.05 | every cell is a separately calibrated test",cols=2)
    for ax,family,title in zip(axes[0],("all","mass"),("All-pair recency","ALT/ALT carrier mass")):
        values=np.array([rows(metrics,method=f"{family}_{t}").power.to_numpy()*100 for t in T])
        heatmap(ax,values,list(range(1,11)),[t//1000 for t in T])
        ax.set(title=title,xlabel="s (x 0.001)",ylabel="TMRCA cutoff (kya)")
    finish(book,fig,"05_time_cutoff_power","Power across TMRCA cutoffs",
        "Cutoffs are 5, 10, 20, 30, 40 and 50 kya; there is no per-replicate maximum over cutoffs. Each cell uses 100 selected replicates and the matching neutral statistic. Comparing cutoffs is exploratory, not an independently validated tuning procedure.",
        "Shorter cutoffs can add discrimination beyond frequency, particularly at s=0.001.")


def fpr(book, metrics):
    fig,axes=canvas("Held-out positional false-positive rates", "1,000 neutral positions, tested once each | these are FPRs, not FDRs among discoveries",cols=2)
    labels=["Truth: all pairs","Truth: carrier mass","Decoded: all pairs","Decoded: carrier mass"]
    for ax,alpha in zip(axes[0],(.05,.01)):
        data=np.array([[rows(metrics,source=src,method=f"{fam}_{t}",alpha=alpha).positional_fpr.iloc[0]*100 for t in T]
            for src,fam in (("truth","all"),("truth","mass"),("decoded","all"),("decoded","mass"))])
        heatmap(ax,data,[t//1000 for t in T],labels,maximum=7 if alpha==.05 else 2,fmt=".1f",cmap="Blues")
        ax.set(title=f"Nominal {alpha:.0%}",xlabel="TMRCA cutoff (kya)")
        ax.tick_params(axis="y",labelsize=13)
    fig.subplots_adjust(left=.17,wspace=.6)
    finish(book,fig,"06_positional_false_positive_rate","Neutral FPR calibration",
        "Cell labels are percentages, with separate color scales for the 5% and 1% panels. Neutral marker-free positions contribute zero carrier mass. The false discovery rate among selected discoveries also depends on the prevalence of selection and is not estimated by this figure. No chromosome-maximum error rate is substituted for the positional rate.",
        "Matched nulls control local false calls close to the nominal level for both truth and decoding.")


def decoding(book, metrics):
    fig,axes=canvas("Carrier mass retains power after decoding", "s = 0.005 only | matched truth and decoded nulls | 100 selected + 1,000 neutral regions",cols=2)
    for ax,family,title in zip(axes[0],("all","mass"),("All-pair recency","ALT/ALT carrier mass")):
        for source in ("truth","decoded"):
            for alpha,style in ((.05,"-"),(.01,"--")):
                values=[rows(metrics,source=source,method=f"{family}_{t}",alpha=alpha).query("s == .005").power.iloc[0] for t in T]
                ax.plot(np.array(T)/1000,values,marker="o",color=C[source],ls=style,label=f"{source.capitalize()}, p <= {alpha:g}")
        ax.set(title=title,xlabel="TMRCA cutoff (kya)",xticks=[5,10,20,30,40,50]);percent(ax)
    legend(fig,ncol=4)
    finish(book,fig,"07_truth_vs_decoding","Truth versus decoded TMRCA",
        "This comparison is available only for the already decoded s=0.005 cohort. Both sources use the same 10,000 sampled pairs and local marker assignments. Decoded frac_recent_T thresholds posterior-mean TMRCA per pair; it does not average posterior mass. Weak-selection arms have not been decoded. Matching nulls calibrates the test but does not remove decoding-related loss of discrimination.",
        "At T=50 kya and p<=0.05, all-pair power falls from 86% to 61%; carrier mass remains at 100%.")


def gain(book, metrics, paired):
    fig,axes=canvas("Coalescence can add information beyond AF", "True TMRCA | p <= 0.05 | comparisons across T and s are exploratory",cols=2)
    ax=axes[0,0]
    for method,label,color in (("af","Archaic AF",C["af"]),("mass_20000","Carrier mass, T=20 kya",C["mass"]),("mass_50000","Carrier mass, T=50 kya",C["all"])):
        curve(ax,rows(metrics,method=method).query("s <= .005"),label,color)
    percent(ax);ax.set_xticks(S[:5],[f"{s:.3f}" for s in S[:5]],rotation=20)
    ax.set_xlabel("Selection coefficient (s)")
    ax=axes[0,1]
    p=paired[(paired.source=="truth") & (paired.alpha==.05) & (paired.s==.001)]
    p=p.set_index("method").loc[[f"mass_{t}" for t in T]]
    x=np.arange(len(T))
    ax.bar(x-.18,p.mass_only,width=.36,color=C["mass"],label="Carrier only (+)")
    ax.bar(x+.18,-p.af_only,width=.36,color=C["af"],label="AF only (-)")
    ax.axhline(0,color="#46596A",lw=1)
    ax.set(xticks=x,xticklabels=[t//1000 for t in T],xlabel="TMRCA cutoff (kya)",ylabel="Discordant detections / 100",title="Paired comparison at s=0.001")
    ax.grid(axis="y",color="#E5E9EE")
    h1,l1=axes[0,0].get_legend_handles_labels();h2,l2=ax.get_legend_handles_labels()
    legend(fig,h1+h2,l1+l2,ncol=3)
    fig.subplots_adjust(top=.64,bottom=.24)
    ax.set_title("Paired comparison at s=0.001",fontsize=18)
    finish(book,fig,"08_gain_beyond_allele_frequency","Incremental gain from coalescence",
        "Left: power for pre-specified AF and carrier-mass scores. Right: within-replicate gains and losses relative to AF at s=0.001; negative bars are losses, not negative probabilities. Null FPR is approximately 5% for each method but not identical. Paired p-values and all cutoffs are in paired_af_comparison.csv; they are unadjusted for examining several cutoffs and arms.",
        "The added coalescence signal is clearest at a shorter cutoff in the weakest selection arm.")


def ihs_comparison(book, ihs, regions):
    fig,axes=canvas("iHS is strongest before the target approaches fixation", "Conventional frequency-standardized iHS | all selected replicates retained",cols=2)
    ax=axes[0,0]
    for method,label,color in (("nearest_core_abs_ihs","Positional iHS",C["ihs"]),("focal_100kb_extreme_fraction","Fixed 100-kb window",C["mass"])):
        z=ihs[(ihs.method==method)&(ihs.alpha==.05)].sort_values("s")
        curve(ax,z,label,color,band=True)
    percent(ax);scoeff(ax);ax.set_title("Power at nominal p <= 0.05")
    ax=axes[0,1]
    selected=regions[regions["mode"]=="selected"]
    g=selected.groupby("s")
    ax.plot(S,g.focal_valid.mean().reindex(S),marker="o",color=C["ihs"],label="Focal SNP has finite iHS")
    fixed=selected.assign(fixed=selected.focal_af==1).groupby("s").fixed.mean().reindex(S)
    ax.plot(S,fixed,marker="s",color=C["all"],label="Selected allele fixed")
    scoeff(ax);percent(ax,"Selected replicates (%)");ax.set_title("Target-site availability")
    h1,l1=axes[0,0].get_legend_handles_labels();h2,l2=ax.get_legend_handles_labels()
    legend(fig,h1+h2,l1+l2,ncol=2)
    fig.subplots_adjust(top=.68)
    finish(book,fig,"09_ihs_and_sweep_completion","iHS power and target-site eligibility",
        "The two iHS curves represent separately calibrated endpoints. Positional iHS uses the nearest finite-scoring core with MAF>=5% within 5 kb; the fixed 100-kb endpoint uses the fraction of finite-scoring cores with |standardized iHS|>2, requiring at least 20 cores. Right: exact focal-SNP eligibility, not eligibility of the fallback positional statistic. Fixation and missing/undefined scores never remove replicates from the denominator.",
        "Near fixation, the selected SNP is often unscorable by iHS; nearby SNPs may still carry signal.")


def spatial(book, distance):
    fig,axes=canvas("The selected signal decays with distance", "s = 0.005; T = 50 kya; p <= 0.05 | every stride uses its own matched neutral position",rows=2,cols=2,
        gridspec_kw=dict(height_ratios=[2.3,1]))
    fig.subplots_adjust(top=.70,bottom=.21,hspace=.26)
    for col,source in enumerate(("truth","decoded")):
        for method,label,key in (("af","Archaic AF","af"),("mass_50000","Carrier mass","mass"),("all_50000","All-pair recency","all")):
            z=distance[(distance.source==source)&(distance.method==method)&(distance.alpha==.05)&(distance.distance_bp.abs()<=500000)].sort_values("distance_bp")
            axes[0,col].plot(z.distance_bp/1000,z.power,color=C[key],label=label)
            axes[1,col].plot(z.distance_bp/1000,z.positional_fpr,color=C[key],lw=1.6)
        axes[0,col].set_title(source.capitalize(),pad=12)
        percent(axes[0,col],"Detection (%)")
        percent(axes[1,col],"FPR (%)",upper=.10)
        axes[1,col].set_yticks([0,.05,.10])
        axes[1,col].axhline(.05,color="#172E41",ls="--",lw=1)
        for row in range(2):
            axes[row,col].axvline(0,color="#A1ABB4",lw=1,ls=":")
            axes[row,col].set_xlim(-500,500)
        axes[0,col].tick_params(labelbottom=False)
        axes[1,col].set_xlabel("Distance from selected allele (kb)")
    legend(fig,ncol=3)
    finish(book,fig,"10_spatial_decay","Detection and FPR across the region",
        "Each point is the fraction of 100 selected or 1,000 neutral replicates called at that coordinate. No smoothing or maximum over positions is used. Detection at linked positions describes the extent of the signal, not a separate causal target or a localization false discovery. Curves are restricted to +/-500 kb for readability; full 10-Mb pointwise results are in the source table.",
        "Strong local carrier signal extends into linked sequence while neutral per-position FPR stays near 5%.")


def profile(key):
    directory=ROOT/"eas_joint_scan/profiles"/key
    receipt=json.loads((directory/"complete.json").read_text())
    path=directory/"features.npz"
    if digest(path)!=receipt["outputs"]["features.npz"]["sha256"]:raise ValueError(f"Corrupt profile {key}")
    INPUTS[str(path)]=digest(path)
    with np.load(path,allow_pickle=False) as f:
        markers=f["grid_markers"];af=np.zeros(len(markers));valid=markers>=0
        af[valid]=f["site_af"][markers[valid]]
        return dict(key=key,position=f["grid"],af=af,counts=f["grid_counts"],n=f["grid_n"],valid=valid)


def examples(book, neutral, region_max, region_keys, focal):
    candidates=pd.DataFrame(dict(key=region_keys,maximum=region_max)).merge(neutral[["key","af"]],on="key")
    candidates=candidates[(candidates.af<.05)&(candidates.maximum>=.4)]
    med=candidates.maximum.median()
    selected_key=candidates.assign(delta=(candidates.maximum-med).abs()).sort_values(["delta","key"]).iloc[0].key
    selected=focal[focal.s==.005].copy()
    selected["distance"]=(selected.sample_af-selected.sample_af.median()).abs()
    rep=int(selected.sort_values(["distance","replicate"]).iloc[0].replicate)
    examples=[profile(selected_key),profile(f"onset50000/rep{rep:04d}")]
    fig,axes=canvas("A peak elsewhere does not invalidate the focal region", "Illustrative examples selected by fixed rules; not an estimate of power",rows=2)
    fig.subplots_adjust(top=.70,bottom=.21,hspace=.50)
    export=[]
    for ax,data,title in zip(axes[:,0],examples,(f"Neutral: {selected_key}",f"Selected s=0.005: replicate {rep:04d}")):
        x=data["position"]/1e6
        for y,label,color in ((data["af"],"Archaic AF",C["af"]),(data["counts"][1,:,2,-1]/10000,"Decoded carrier mass",C["mass"])):
            ax.plot(x,y,label=label,color=color,lw=1.7)
        ax.axvspan(4.95,5.05,color="#D8A82F",alpha=.3)
        ax.axvline(5,color="#AD7C14",ls="--",lw=1)
        ax.set_title(title,loc="left",fontsize=16,pad=7)
        percent(ax,"Score",upper=1.02);ax.set_xlim(0,10)
        export.append(pd.DataFrame(dict(key=data["key"],position_bp=data["position"],af=data["af"],decoded_mass_T50000=data["counts"][1,:,2,-1]/10000)))
    axes[1,0].set_xlabel("Position in scored region (Mb)")
    legend(fig,ncol=2)
    pd.concat(export).to_csv(DATA/"example_profiles.csv",index=False)
    finish(book,fig,"11_example_regions","Example neutral and selected regions",
        f"Neutral example {selected_key}: choose focal AF<5% and whole-region maximum AF>=40%, then the maximum nearest the median of eligible examples (ties by ID). Selected example {rep:04d}: AF nearest the s=0.005 median (ties by replicate). Gold shading marks a fixed 100-kb interval centered at 5 Mb. AF and carrier mass have different scales of interpretation despite both lying in [0,1].",
        "The gene or position defines the scope of the null; peaks outside it are irrelevant.",category="backup")


def stringency(book, metrics, pooled):
    fig,axes=canvas("Greater stringency trades power for fewer false calls", "Pointwise tests | p=0.001 uses pooled ranks and has very limited neutral-tail resolution",cols=2)
    colors=(C["mass"],C["all"],C["ihs"])
    for ax,method,title in zip(axes[0],("af","mass_50000"),("Archaic AF","Carrier mass, T=50 kya")):
        for alpha,color in zip((.05,.01,.001),colors):
            z=rows(pooled if alpha==.001 else metrics,method=method,alpha=alpha)
            ax.plot(z.s,z.power,marker="o",color=color,label=f"p <= {alpha:g}"+(" (pooled)" if alpha==.001 else " (5-fold)"))
        percent(ax);scoeff(ax);ax.set_title(title)
    legend(fig,ncol=3)
    finish(book,fig,"12_threshold_stringency","Power at stricter thresholds",
        "p<=0.05 and 0.01 use 400 neutral calibration positions per fold. The separate p<=0.001 sensitivity uses all 1,000 neutrals for selected ranks and leave-one-out ranks for neutrals. One extreme neutral rank gives 1/1,000 calls for the shown scores; this is not precise independent validation of a 0.1% tail. Comparisons remain positional, not family-wise over a genome.",
        "The s=0.005 signal remains strong at stringent ranks; very weak selection loses substantial power.",category="backup")


def genealogies(book, focal):
    fig,axes=canvas("Carrier mass and allele frequency are related", "True focal genealogies | 1,000 selected regions | ancestral branches are not a count of founding haplotypes",cols=2)
    ax=axes[0,0]
    for cutoff,color,label in ((20000,C["mass"],"T=20 kya"),(50000,C["all"],"T=50 kya")):
        ax.scatter(focal.sample_af,focal[f"mass_{cutoff}"],s=12,color=color,alpha=.32,linewidths=0,label=label)
    x=np.linspace(0,1,200)
    ax.plot(x,x*x,color="#172E41",ls="--",label="AF squared (approximation)",lw=1.5)
    percent(ax,"ALT/ALT carrier-mass score");ax.set_xlabel("Selected allele frequency")
    ax.xaxis.set_major_formatter(PercentFormatter(1,decimals=0));ax.set_xlim(0,1)
    ax=axes[0,1]
    lineage=focal.assign(category=np.minimum(focal.carrier_lineages_at_50k,3))
    values=pd.crosstab(lineage.s,lineage.category,normalize="index").reindex(index=S,columns=[1,2,3],fill_value=0)
    bottom=np.zeros(10)
    for cat,color,label in ((1,C["mass"],"One branch"),(2,C["all"],"Two branches"),(3,C["ihs"],"Three or more")):
        ax.bar(S,values[cat],bottom=bottom,width=.00068,color=color,label=label)
        bottom+=values[cat].to_numpy()
    percent(ax,"Selected regions (%)");scoeff(ax);ax.set_title("Carrier branches just before 50 kya")
    h1,l1=axes[0,0].get_legend_handles_labels();h2,l2=ax.get_legend_handles_labels()
    legend(fig,h1+h2,l1+l2,ncol=3)
    fig.subplots_adjust(top=.68)
    finish(book,fig,"13_carrier_genealogies","Why AF and carrier mass can be similar",
        "Carrier mass equals the sampled ALT/ALT pair fraction times within-carrier frac_recent_T. AF squared approximates the first factor for random distinct pairs, with finite-panel and pair-sampling differences; it is not an exact upper bound. Branches are the sampled-carrier ancestral groups immediately below the 50-kya cutoff, verified against pairwise TMRCA, not a census of introgressing founders. The model has a strong archaic bottleneck.",
        "When most carriers share recent ancestry by T=50 kya, weighting by coalescence adds less beyond AF.",category="backup")


def demography(book, cfg):
    dem=model(cfg).model
    times=[0.];sizes=[dem.populations[0].initial_size]
    for event in dem.events:
        if type(event).__name__=="PopulationParametersChange" and event.population=="EAS" and event.initial_size is not None:
            times.append(event.time*cfg["generation_time_years"]);sizes.append(event.initial_size)
    history=pd.DataFrame(dict(years=times,Ne=sizes)).sort_values("years")
    history.to_csv(DATA/"eas_demography.csv",index=False)
    fig,axes=canvas("Demography is shared by neutral and selected simulations", "The generative EAS history varies with time; the current decoder uses a fixed Ne",cols=2)
    ax=axes[0,0]
    z=history[history.years<=50000]
    x=np.r_[z.years,50000]/1000;y=np.r_[z.Ne,z.Ne.iloc[-1]]
    ax.step(x,y,where="post",color=C["truth"],label="Simulated EAS Ne(t)")
    decoder_ne=cfg["decoder_scaled_mutation_rate"]/(4*cfg["mutation_rate"])
    ax.axhline(decoder_ne,color=C["decoded"],ls="--",label=f"Decoder Ne = {decoder_ne:,.0f}")
    ax.set(xlabel="Time before present (kya)",ylabel="Diploid effective population size",yscale="log",xlim=(0,50))
    ax.grid(axis="y",color="#E5E9EE")
    ax=axes[0,1];ax.axis("off")
    text="Archaic branch\n\nSplit: 700 kya\nBaseline Ne: 3,600\nPulse bottleneck: Ne=10\nfor 100 generations\n\nMutation rate: 1.25 x 10^-8\nRecombination rate: 1 x 10^-8\nGeneration time: 25 years"
    ax.text(.02,.98,text,transform=ax.transAxes,ha="left",va="top",fontsize=16.5,linespacing=1.35)
    legend(fig,ncol=2)
    finish(book,fig,"14_demography_and_assumptions","Demography and decoder assumptions",
        "The EAS trajectory is read from the same PHLASH-derived model used by the archived simulations. The Ne history is displayed only to 50 kya; simulation genealogies also contain the older archaic branch. Decoder Ne is theta/(4*mu)=15,000. Neutral and selected simulations share the archaic bottleneck and rate assumptions. No variable-Ne Gamma-SMC flow field or empirical recombination-map benchmark is claimed here.",
        "Matched nulls address calibration under the model; they do not establish robustness to every model misspecification.",category="backup")


def pair_classes(book, all_profiles):
    fig,axes=canvas("Recent coalescence among archaic ALT/ALT pairs", "Focal position | pairs pooled within each defined class | s=0.005 selected versus neutral",cols=2)
    records=[]
    for ax,source,src in zip(axes[0],("truth","decoded"),range(2)):
        for selected,style in ((True,"-"),(False,"--")):
            group=[p for p in all_profiles if p["key"].startswith("onset50000/")==selected and p["valid"][500]]
            counts=np.array([p["counts"][src,500] for p in group]).sum(axis=0)
            totals=np.array([p["n"][500] for p in group]).sum(axis=0)
            for cls,label,color in ((2,"ALT/ALT",C["mass"]),(0,"REF/REF",C["af"])):
                values=counts[cls]/totals[cls]
                cohort="Selected" if selected else "Neutral"
                ax.plot(np.array(T)/1000,values,ls=style,color=color,marker="o",label=f"{cohort}: {label}")
                records.extend(dict(source=source,cohort=cohort,pair_class=label,T=t,frac_recent=float(v),regions=len(group),pairs=int(totals[cls])) for t,v in zip(T,values))
        percent(ax,"Within-class frac_recent_T")
        ax.set(title=source.capitalize(),xlabel="TMRCA cutoff (kya)",xticks=[5,10,20,30,40,50])
    legend(fig,ncol=4)
    pd.DataFrame(records).to_csv(DATA/"pair_class_recency.csv",index=False)
    finish(book,fig,"15_pair_class_coalescence","Pair-class coalescence diagnostic",
        "ALT/ALT and REF/REF are haplotype-pair classes at the focal archaic marker. Only regions with an assigned marker contribute to this within-class diagnostic: 100 selected and 169 neutral. Counts are pooled across pairs, so regions with larger classes contribute more weight; this is not a replicate-average or the unconditional neutral denominator used for FPR. REF/REF is a diagnostic and is not subtracted from the carrier score; ALT/REF is not scored.",
        "The carrier score combines this within-ALT/ALT recency with the fraction of the complete pair panel that is ALT/ALT.",category="backup")


def write_index():
    lines=["# EAS lab meeting figures - 21 September 2026", "",
        "15 figures in slide-sized 16:9 layout. Each has a 3199 x 1800 PNG and a vector PDF with editable text. The combined PDF has one figure per page; full captions and suggested speaking points are below.", "",
        "The figures use the completed 1,000 neutral + 1,000 selected array (100 selected per s). Gamma-SMC decoding is available only for s=0.005 and the neutral cohort. Simulation generation and new Gamma-SMC decoding were not started for this figure request.", "",
        "Suggested main narrative: 01 -> 02 -> 03 -> 04 -> 05 -> 07 -> 08 -> 09 -> 10. Figure 06 verifies FPR calibration; 11-15 are backup figures.", "",
        "All primary error rates are positional FPR, not FDR among discoveries. Tests at a pre-specified gene require the same statistic in a matching neutral interval. No maximum across unrelated positions or T cutoffs enters the primary p-values. The 100-kb iHS window is a separate endpoint.", ""]
    for item in FIGURES:
        lines.extend([f"## {item['number']:02d}. {item['title']} ({item['category']})", "",
            f"[PNG](figures/{item['stem']}.png) | [Vector PDF](figures/{item['stem']}.pdf)", "",
            "**Takeaway:** "+item["takeaway"], "", "**Caption:** "+item["caption"], ""])
    lines.extend(["## Interpretation limits", "",
        "Selected power is conditional on allele survival and observation. The unconditional positional null includes marker-free positions; it is not a null conditional on first selecting an observed archaic SNP. The archaic-marker oracle is more restrictive than simply lying on an introgressed segment and is not a literal finite-reference Neanderthal/Denisovan/African ascertainment experiment.", "",
        "This array uses a strong shared archaic bottleneck, uniform mutation/recombination rates and isolated 10-Mb regions. It does not demonstrate genome-wide FDR control, causal localization, or a gain from multiple independent introgressing haplotypes. Cutoff comparisons are exploratory. Wilson intervals are conditional on the estimated null and do not include its Monte Carlo uncertainty.", "",
        "## Methods references", "",
        "Voight et al. (2006), iHS: https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0040072", "",
        "scikit-allel implementation: https://scikit-allel.readthedocs.io/en/stable/_modules/allel/stats/selection.html", ""])
    (OUT/"README.md").write_text("\n".join(lines))
    (OUT/"figure_index.json").write_text(json.dumps(FIGURES,indent=2)+"\n")


def main():
    FIG.mkdir(parents=True,exist_ok=True);DATA.mkdir(parents=True,exist_ok=True)
    metrics=table(POS,"metrics.csv")
    focal=table(OUT/"analysis","focal_truth.csv")
    neutral=table(POS,"neutral_positions.csv")
    paired=table(POS,"paired_af_comparison.csv",manifest="supplement_manifest.json")
    pooled=table(POS,"pooled_rank_metrics.csv",manifest="supplement_manifest.json")
    ihs=table(ROOT/"eas_ihs_h400","metrics.csv")
    ihs_regions=table(ROOT/"eas_ihs_h400","regions.csv")
    distance=table(ROOT/"eas_positional_distance_h400","metrics.csv")
    assert len(focal)==1000 and (focal.groupby("s").size()==100).all()
    assert (ihs.selected_regions==100).all() and (metrics.selected_regions==100).all()
    ref=ROOT/"eas_allele_class_ablation"
    regions=table(ref,"regions.csv")
    methods=json.loads(checked(ref,"methods.json").read_text())
    with np.load(checked(ref,"scores.npz"),allow_pickle=False) as saved:
        region_max=saved["scores"][0,0,0,regions.key.str.startswith("neutral/"),next(i for i,m in enumerate(methods) if m["name"]=="af_r1")]
    region_keys=regions.loc[regions.key.str.startswith("neutral/"),"key"].to_numpy()
    with ThreadPoolExecutor(max_workers=4) as pool:
        all_profiles=list(pool.map(profile,regions.key.tolist()))
    cfg=json.loads((ROOT/"eas_q02_h400/manifest.json").read_text())["config"]
    for name,frame in (("positional_metrics",metrics),("focal_truth",focal),("ihs_metrics",ihs),("ihs_regions",ihs_regions),
                       ("paired_af_comparison",paired),("pooled_rank_metrics",pooled),("distance_metrics",distance),
                       ("neutral_focal_af",neutral)):
        frame.to_csv(DATA/f"{name}.csv",index=False)
    pd.DataFrame(dict(key=region_keys,maximum_af=region_max)).to_csv(DATA/"neutral_region_max_af.csv",index=False)
    with PdfPages(OUT/"EAS_lab_meeting_figures.pdf",metadata={"Title":"EAS introgression selection - lab meeting figures", "CreationDate":DATE,"ModDate":DATE}) as book:
        design(book)
        local_null(book,neutral,region_max)
        allele_frequency(book,focal,neutral)
        power(book,metrics,ihs)
        time_power(book,metrics)
        fpr(book,metrics)
        decoding(book,metrics)
        gain(book,metrics,paired)
        ihs_comparison(book,ihs,ihs_regions)
        spatial(book,distance)
        examples(book,neutral,region_max,region_keys,focal)
        stringency(book,metrics,pooled)
        genealogies(book,focal)
        demography(book,cfg)
        pair_classes(book,all_profiles)
    write_index()
    (OUT/"layout_checks.json").write_text(json.dumps(LAYOUT,indent=2)+"\n")
    (OUT/"figure_provenance.json").write_text(json.dumps(dict(source_sha256=digest(Path(__file__)),
        inputs=INPUTS,config=cfg,matplotlib=matplotlib.__version__,numpy=np.__version__,
        seed_policy="No random sampling/jitter; original simulation seeds and pair seed retained.",
        simulations_started=False,gamma_smc_decoding_started=False,selected_regions=1000,neutral_regions=1000,
        decoded_selected_regions=100,format="16:9; PNG 240 dpi; vector PDF, TrueType editable text"),indent=2)+"\n")
    print(json.dumps(dict(figures=len(FIGURES),destination=str(OUT),layout="passed")),flush=True)


if __name__=="__main__":
    main()
