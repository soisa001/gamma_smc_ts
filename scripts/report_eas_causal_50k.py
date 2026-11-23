"""Summarize causal-site simulations and render power and frequency plots."""
import argparse
import json
import os
from pathlib import Path
import shutil

os.environ["MPLBACKEND"] = "Agg"
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import PercentFormatter
import numpy as np
import report_eas_strength_series as previous

base = previous.base
OUT = previous.OUT.parent/"eas_causal_50k_20261123"
RAW = Path.home()/"sim/eas_causal_50k_20261123"


def neutral_banks():
    frozen = OUT/"neutral_reference.npz"
    receipt = OUT/"neutral_reference.json"
    if receipt.exists():
        assert base.sha(frozen) == base.read(receipt)["sha256"]
        with np.load(frozen) as data:
            return {key:data[key] for key in data.files}
    path = previous.OUT/"mass.npz"
    audit = base.read(previous.OUT/"inputs.json")
    assert base.sha(path) == audit["cache_sha256"]
    for pop in ("de_novo_50k", "introgression"):
        assert audit["neutral_audit"][pop]["fixed"] == 0
    for run in ("denovo_comparison_20261111", "denovo_comparison_cached_20261111"):
        report = base.read(previous.SIM/run/"progress.json")
        neutral = next(r for r in report["strata"] if r["key"] == "50k_neutral")
        assert neutral["fixed"] == 0
    assert base.read(previous.SIM/"segregating_neutral_100_20261110/progress.json")["fixed"] == 0
    result = {}
    with np.load(path) as data:
        for origin in base.ORIGINS:
            name = "de_novo_50k" if origin == "de_novo" else origin
            for role in ("neutral","test","neutral_ids","test_ids"):
                result[origin+"_"+role] = data[name+"_"+role]
            assert set(result[origin+"_neutral_ids"]).isdisjoint(result[origin+"_test_ids"])
    for origin in base.ORIGINS:
        if origin == "introgression":
            path = previous.OUT.parent/"segregating_neutral_fpr_100_20261110/parameters.npz"
            assert base.sha(path) == base.read(path.parent/"output_hashes.json")[path.name]
            with np.load(path) as archive:
                np.testing.assert_array_equal(archive["test_ids"], result[origin+"_test_ids"])
                np.testing.assert_array_equal(archive["test_cdfs"][:,2], result[origin+"_test"])
                result[origin+"_test_af"] = archive["test_af"]
        else:
            traces = []
            for source in result[origin+"_test_ids"]:
                path = Path(str(source))
                rec = base.read(path/"complete.json")
                trace, coeff = trajectory(path/"trajectory.tsv", rec["trajectory_sha256"])
                assert np.all(coeff == 0)
                traces.append(trace)
            result[origin+"_test_af"] = np.array(traces)
        assert result[origin+"_test_af"].shape == (100,2001)
    np.savez_compressed(frozen, **result)
    base.save(receipt, dict(sha256=base.sha(frozen), source_bank=audit["cache_sha256"],
              focal_bp=500000, age_years=50000, calibration=100, heldout=100,
              fixed_in_audited_neutral_proposals=0))
    return result


def null_rows(data):
    rows = []
    for origin in base.ORIGINS:
        p = previous.rank_p(data[origin+"_neutral"], data[origin+"_test"])
        for window, (low, high) in previous.WINDOWS.items():
            for alpha in (.05,.01):
                calls = (p[:,low-1:high] < alpha).any(axis=1)
                rows.append(dict(origin=origin, window=window, alpha=alpha, fp=int(calls.sum()), n=len(calls), fpr=calls.mean()))
    return rows


def trajectory(path, expected_digest):
    assert base.sha(path) == expected_digest
    years, af, diploids, coeff, copies = np.loadtxt(path, skiprows=1).T
    np.testing.assert_array_equal(years, np.arange(50000,-1,-25))
    np.testing.assert_array_equal(af, copies/(2*diploids))
    assert np.all(af > 0) and np.all(af <= 1)
    assert not np.any(np.diff((af==1).astype(int)) < 0)
    return af, coeff


def prepare(raw):
    raw = raw.resolve()
    OUT.mkdir(exist_ok=True, parents=True)
    contract = base.read(raw/"manifest.json")
    data = neutral_banks()
    af = {}
    index = []
    for origin in base.ORIGINS:
        af[origin+"_neutral"] = data.pop(origin+"_test_af")
        for window in previous.WINDOWS:
            cfg = contract["configs"][window]
            for s in contract["strengths"]:
                key = f"{origin}_{window}_s{s:g}"
                curves, fixed, recent, founders = [], [], [], set()
                for slot in range(contract["target"]):
                    directory = raw/"groups"/origin/window/f"{slot:03d}"/f"s{s:g}"
                    rec = base.read(directory/"complete.json")
                    assert rec["slot"]==slot and rec["selection"]==s and rec["origin"]==origin
                    assert rec["position"]==500000 and rec["ancestry_start_years"]==50000 and rec["stop_years_ago"]==0
                    assert rec["onset_sha256"] not in founders
                    founders.add(rec["onset_sha256"])
                    assert base.sha(directory/"scores.npz")==rec["scores_sha256"]
                    with np.load(directory/"scores.npz") as scores:
                        np.testing.assert_array_equal(scores["positions"],[500000])
                        np.testing.assert_array_equal(scores["valid_pairs"],[1000])
                        np.testing.assert_array_equal(scores["cutoffs"],np.arange(1000,50001,1000))
                        curves.append(scores["cdfs"][2,0])
                    trace, coeff = trajectory(directory/"trajectory.tsv", rec["trajectory_sha256"])
                    years = np.arange(50000,-1,-25)
                    active = (years <= cfg["selection_start_years"])
                    if cfg["selection_end_years"]:
                        active &= (years > cfg["selection_end_years"])
                    np.testing.assert_allclose(coeff, np.where(active, s, 0))
                    assert trace[-1] == rec["final_frequency"]
                    fixed.append(trace[-1] == 1)
                    if window=="10k" and s in (.01,.03):
                        recent.append(trace)
                    index.append(dict(case=key,slot=slot,attempt=rec["attempt"],trajectory=rec["trajectory_sha256"],
                                      scores=rec["scores_sha256"],founder=rec["onset_sha256"],unit=rec["unit"]))
                data[key] = np.array(curves)
                data[key+"_fixed"] = np.array(fixed)
                if recent:
                    af[key] = np.array(recent)
    np.savez_compressed(OUT/"mass.npz", **data)
    np.savez_compressed(OUT/"recent_af.npz", **af)
    shutil.copyfile(raw/"manifest.json", OUT/"simulation_manifest.json")
    base.save(OUT/"inputs.json", dict(records=index, hashes={p:base.sha(OUT/p) for p in ("mass.npz","recent_af.npz","simulation_manifest.json")}))


def af_plot(af, origin, s):
    years = np.arange(50000,-1,-25)/1000
    fig, axes = plt.subplots(2,1,figsize=(11,8.5),sharex=True,sharey=True)
    fig.subplots_adjust(top=.80,left=.11,right=.97,bottom=.12,hspace=.36)
    name = "De novo" if origin=="de_novo" else "1% Neanderthal introgression"
    fig.suptitle(f"{name} at 50 kya", y=.97, fontsize=23)
    for ax, key, title, color in zip(axes, (origin+"_neutral",f"{origin}_10k_s{s:g}"),
                                     ("Neutral",f"Selection from 10 kya to present, s = {s:g}"),("#506679","#994226")):
        a = af[key]
        ax.plot(years,a.T, color=color, lw=.6, alpha=.28)
        ax.plot(years,np.median(a, axis=0),color=color,lw=2.6)
        if key != origin+"_neutral":
            ax.axvline(10,color="#222222",ls="--",lw=1.3)
        ax.set_title(title, loc="left", fontsize=16)
        ax.set(xlim=(50,0),ylim=(-.015,1.02),ylabel="Allele frequency")
        ax.yaxis.set_major_formatter(PercentFormatter(1))
        ax.spines[["top","right"]].set_visible(False)
        ax.grid(axis="y",alpha=.16)
    fig.legend([Line2D([],[],color="#627381",lw=.7),Line2D([],[],color="#384853",lw=2.6),Line2D([],[],color="#222222",ls="--")],
               ["100 individual trajectories","Median","Selection onset"],ncol=3,loc="upper center",bbox_to_anchor=(.51,.92),frameon=False)
    axes[-1].set_xlabel("Time before present (kya)")
    for ext in ("pdf","png"):
        fig.savefig(OUT/f"af_{origin}_10k_s{s:g}.{ext}",dpi=175)
    plt.close(fig)


def timing_examples(raw):
    """Pair two identifiable selected histories with their 100 held-out controls."""
    specs = [("de_novo","10k",3,.01), ("introgression","50_30k",1,.03)]
    controls = neutral_banks()
    years = np.arange(50000,-1,-25)/1000
    plt.rcParams.update({"font.size":15,"pdf.fonttype":42,"ps.fonttype":42})
    fig,axes = plt.subplots(2,1,figsize=(11,8.5),sharex=True)
    fig.subplots_adjust(top=.81,left=.12,right=.96,bottom=.11,hspace=.38)
    fig.suptitle("Focal allele frequency from 50 kya to present",fontsize=22,y=.97)
    for ax,(origin,window,slot,s) in zip(axes,specs):
        directory = raw/"groups"/origin/window/f"{slot:03d}"/f"s{s:g}"
        r = base.read(directory/"complete.json")
        af, coeff = trajectory(directory/"trajectory.tsv",r["trajectory_sha256"])
        assert r["position"]==500000 and af[-1]==r["final_frequency"]
        lo, hi = r["selection_end_years"]/1000, r["selection_start_years"]/1000
        ax.axvspan(lo,hi,color="#f6dfba",alpha=.6,zorder=0)
        ax.plot(years,controls[origin+"_test_af"].T,color="#586f83",lw=.6,alpha=.25)
        ax.plot(years,af,color="#993f23",lw=2.2)
        if origin=="de_novo":
            assert np.all(coeff[years>10]==0) and np.all(coeff[years<=10]==s)
            ax.scatter([10],[af[1600]],s=32,color="#993f23",zorder=4)
            ax.set_title(f"De novo, s = {s:g}   ·   Selection from 10 kya to present",loc="left",fontsize=16)
        else:
            first_fixed = np.flatnonzero(af==1)[0]
            ax.scatter([years[first_fixed]],[1],s=30,color="#993f23",zorder=4)
            ax.annotate(f"Fixed at {years[first_fixed]:.1f} kya",(years[first_fixed],1),
                        xytext=(12,-20),textcoords="offset points",ha="left",fontsize=13)
            ax.set_title(f"1% introgression, s = {s:g}   ·   Selection from 50 to 30 kya",loc="left",fontsize=16)
        ax.set(xlim=(50,0),ylim=(-.02,1.045),ylabel="Allele frequency")
        ax.yaxis.set_major_formatter(PercentFormatter(1))
        ax.grid(axis="y",alpha=.15)
        ax.spines[["top","right"]].set_visible(False)
    from matplotlib.patches import Patch
    fig.legend([Line2D([],[],color="#586f83",lw=1),Line2D([],[],color="#993f23",lw=2.2),Patch(color="#f6dfba")],
               ["100 neutral trajectories","One selected trajectory","Selection window"],ncol=3,
               loc="upper center",bbox_to_anchor=(.51,.914),frameon=False,fontsize=14)
    axes[-1].set_xlabel("Time before present (kya)")
    for ext in ("png","pdf"):
        fig.savefig(OUT/f"timing_examples.{ext}",dpi=175)
    plt.close(fig)


def report():
    for name, digest in base.read(OUT/"inputs.json")["hashes"].items():
        assert base.sha(OUT/name)==digest
    with np.load(OUT/"mass.npz") as data:
        fpr = null_rows(data)
        rows = []
        for r in fpr:
            o, w, a = r["origin"],r["window"],r["alpha"]
            low, high = previous.WINDOWS[w]
            for s in base.read(OUT/"simulation_manifest.json")["strengths"]:
                k = f"{o}_{w}_s{s:g}"
                p = previous.rank_p(data[o+"_neutral"], data[k])
                calls = (p[:,low-1:high]<a).any(axis=1)
                fixed = data[k+"_fixed"]
                rows.append(dict(origin=o,window=w,s=s,alpha=a,tp=int(calls.sum()),n=len(calls),tpr=calls.mean(),
                                 fixed=int(fixed.sum()),fixed_calls=int(calls[fixed].sum()),fp=r["fp"],neutral_n=r["n"],fpr=r["fpr"]))
    previous.write_table(OUT/"power.tsv",rows)
    previous.write_table(OUT/"fpr.tsv",fpr)
    previous.plot(rows, out=OUT)
    with np.load(OUT/"recent_af.npz") as data:
        for origin in base.ORIGINS:
            for s in (.01,.03):
                af_plot(data,origin,s)


if __name__=="__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("phase",choices=("null","examples","prepare","report"))
    p.add_argument("--raw",type=Path,default=RAW)
    a = p.parse_args()
    if a.phase=="null":
        OUT.mkdir(exist_ok=True,parents=True)
        previous.write_table(OUT/"fpr.tsv",null_rows(neutral_banks()))
    elif a.phase=="examples":
        timing_examples(a.raw)
    elif a.phase=="prepare":
        prepare(a.raw)
    else:
        report()
