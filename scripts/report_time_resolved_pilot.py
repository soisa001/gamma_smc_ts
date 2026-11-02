#!/usr/bin/env python3
"""Held-out, whole-region calibration of a completed bounded EAS pilot."""
import os
os.environ["MPLBACKEND"]="Agg"
for key in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"):
    os.environ[key]="1"
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import label
from scipy.stats import norm

METHODS=("truth","median","mass")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def wilson(k,n):
    if not n: return [None,None]
    z=norm.ppf(.975); p=k/n; den=1+z*z/n
    mid=(p+z*z/(2*n))/den
    half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return [float(mid-half),float(mid+half)]


def empirical_p(scores,calibration_maxima):
    """Conservative +1 upper-tail probability, including ties."""
    ordered=np.sort(calibration_maxima)
    return (1+len(ordered)-np.searchsorted(ordered,scores,side="left"))/(1+len(ordered))


def interval_values(cdf,widths):
    parts=[]; intervals=[]
    for width in widths:
        parts.append(cdf[...,width:]-cdf[...,:-width])
        intervals.extend((a,a+width) for a in range(1,51-width))
    return np.concatenate(parts,axis=-1),np.asarray(intervals)


def expand_support(calls,intervals):
    """Union called windows into the original 49 one-kyr support cells."""
    result=np.zeros((calls.shape[0],49),dtype=bool)
    for k,(a,b) in enumerate(intervals):
        result[:,a-1:b-1] |= calls[:,k,None]
    return result


def temporal_stats(calls):
    """One-kyr support at the selected coordinate; truth [20,30)."""
    truth=np.zeros(49,dtype=bool);truth[19:29]=True
    overlap=int(np.sum(calls&truth));n=int(calls.sum())
    components,ncomp=label(calls)
    hit=overlap>0
    off=sum(not np.any(truth[components==i]) for i in range(1,ncomp+1))
    return dict(center_overlap=int(hit),temporal_precision=overlap/n if n else None,
        temporal_recall=overlap/10,temporal_f1=2*overlap/(n+10),iou=overlap/(n+10-overlap),
        called_kyr=n,overlap_kyr=overlap,center_off_target=int(off>0),
        center_mistimed=int(n>0 and not hit),center_components=ncomp)


def evaluate_calls(calls,positions,selected):
    center=len(positions)//2
    near=np.abs(positions-positions[center])<=100000
    truth=np.zeros(49,dtype=bool);truth[19:29]=True
    components,ncomp=label(calls)
    joint=bool(np.any(calls & near[:,None] & truth[None,:]))
    tp=int(selected and joint)
    result=dict(any_call=int(calls.any()),spatial_hit=int(np.any(calls[near])),joint_hit=int(joint),
        event_tp=tp,event_fp=int(ncomp-tp),event_fn=int(selected)-tp,components=int(ncomp))
    if selected: result.update(temporal_stats(calls[center]))
    return result


def read_units(root):
    manifest=json.loads((root/"manifest.json").read_text())
    signature=sha(root/"manifest.json")
    records=[];arrays=[]
    positions=np.loadtxt(root/"positions.txt",dtype=int)
    seen=set()
    for receipt in sorted(root.glob("*/complete.json")):
        rec=json.loads(receipt.read_text())
        assert rec["config_sha256"]==signature
        assert rec["seed"] not in seen
        seen.add(rec["seed"])
        path=receipt.parent/"scores.npz"
        assert sha(path)==rec["scores_sha256"]
        assert sha(receipt.parent/"observed.trees")==rec["observed_sha256"]
        with np.load(path) as a:
            np.testing.assert_array_equal(a["positions"],positions)
            np.testing.assert_array_equal(a["cutoffs"],np.arange(1000,50001,1000))
            np.testing.assert_array_equal(a["valid_pairs"],np.full(len(positions),manifest["config"]["pairs"]))
            arrays.append(a["cdfs"])
        rec["cohort"]=rec["unit"].rsplit("_",1)[0]
        records.append(rec)
    return manifest,pd.DataFrame(records),np.stack(arrays),positions


def cumulative_baselines(records,cdfs,positions):
    """Prespecified cumulative comparators; their cutoffs are not episode dates."""
    cohort=records.cohort.to_numpy()
    train=cohort=="development";cal=cohort=="calibration"
    selected=cohort=="selected";neutral=cohort=="neutral_test"
    near=abs(positions-positions[len(positions)//2])<=100000
    rows=[]
    for name,indices in (("fixed_50kya",[49]),("cumulative_scan_1_to_50kya",list(range(50)))):
        for j,method in enumerate(METHODS):
            values=cdfs[:,j][:,:,indices]
            center=values[train].mean(axis=(0,1));scale=values[train].std(axis=(0,1),ddof=1)
            valid=scale>0
            z=np.full_like(values,-np.inf)
            z[...,valid]=(values[...,valid]-center[valid])/scale[valid]
            maxima=z[cal].max(axis=(1,2))
            called=empirical_p(z,maxima)<=.05
            any_call=called.any(axis=(1,2))
            spatial=called[:,near].any(axis=(1,2))
            for metric,mask,indicator in (("neutral_fpr",neutral,any_call),("selected_any_call",selected,any_call),("spatial_tpr",selected,spatial)):
                k=int(indicator[mask].sum());n=int(mask.sum())
                lo,hi=wilson(k,n)
                rows.append(dict(baseline=name,method=method,metric=metric,k=k,n=n,rate=k/n,ci95_low=lo,ci95_high=hi))
    return pd.DataFrame(rows)


def report(root,out,widths=(1,),reference=None):
    out.mkdir(parents=True,exist_ok=True)
    manifest,records,cdfs,positions=read_units(root)
    if reference is not None:
        ref_manifest,ref_records,ref_cdfs,ref_positions=read_units(reference)
        a={k:v for k,v in manifest["config"].items() if k!="selection"}
        b={k:v for k,v in ref_manifest["config"].items() if k!="selection"}
        assert a==b, "Reference nulls differ in more than the unused selected-arm coefficient"
        np.testing.assert_array_equal(positions,ref_positions)
        keep=ref_records.cohort.to_numpy()!="selected"
        target=records.cohort.to_numpy()=="selected"
        records=pd.concat([ref_records[keep],records[target]],ignore_index=True)
        cdfs=np.concatenate([ref_cdfs[keep],cdfs[target]],axis=0)
        assert records.seed.is_unique
    cohort=records.cohort.to_numpy()
    development=cohort=="development"; calibration=cohort=="calibration"
    assert development.sum()>=10 and calibration.sum()>=19, "Insufficient independent nulls to report p<=.05"
    mass,intervals=interval_values(cdfs,widths)
    rows=[]; summaries=[]; centers=[]; scales=[]; maxima=[]
    test=(cohort=="neutral_test")|(cohort=="selected")
    for j,method in enumerate(METHODS):
        # Pooled positions define nuisance scaling, not independent replicates.
        # Held-out calibration still uses one full-grid maximum per region.
        center=mass[development,j].mean(axis=(0,1))
        scale=mass[development,j].std(axis=(0,1),ddof=1)
        eligible=scale>0
        z=np.full(mass[:,j].shape,-np.inf)
        z[...,eligible]=(mass[:,j][...,eligible]-center[eligible])/scale[eligible]
        cal_max=z[calibration].max(axis=(1,2))
        centers.append(center);scales.append(scale);maxima.append(cal_max)
        for i in np.flatnonzero(test):
            p=empirical_p(z[i],cal_max)
            calls=expand_support(p<=.05,intervals)
            metrics=evaluate_calls(calls,positions,cohort[i]=="selected")
            best=np.unravel_index(np.argmax(z[i]),z[i].shape)
            center_peak=float(intervals[int(np.argmax(z[i,len(positions)//2]))].mean())
            rows.append(dict(unit=records.iloc[i].unit,cohort=cohort[i],method=method,
                global_p=float(empirical_p(z[i].max(),cal_max)),peak_position_bp=int(positions[best[0]]),
                peak_time_kya=float(intervals[best[1]].mean()),center_peak_time_kya=center_peak,
                center_peak_error_kya=abs(center_peak-25),final_frequency=float(records.iloc[i].final_frequency),**metrics))
        frame=pd.DataFrame([r for r in rows if r["method"]==method])
        sel=frame[frame.cohort=="selected"]; null=frame[frame.cohort=="neutral_test"]
        counts={"neutral_fpr":(int(null.any_call.sum()),len(null)),
                "selected_any_call":(int(sel.any_call.sum()),len(sel)),
                "spatial_tpr":(int(sel.spatial_hit.sum()),len(sel)),
                "joint_tpr":(int(sel.joint_hit.sum()),len(sel)),
                "center_overlap_tpr":(int(sel.center_overlap.sum()),len(sel)),
                "center_off_target_rate":(int(sel.center_off_target.sum()),len(sel)),
                "center_mistimed_rate":(int(sel.center_mistimed.sum()),len(sel))}
        rates={name:dict(k=k,n=n,rate=k/n if n else None,ci95=wilson(k,n)) for name,(k,n) in counts.items()}
        tp=int(frame.event_tp.sum());fp=int(frame.event_fp.sum());fn=int(frame.event_fn.sum())
        summaries.append(dict(method=method,eligible_bins=int(eligible.sum()),rates=rates,
            mean_temporal_f1=float(sel.temporal_f1.mean()),mean_temporal_recall=float(sel.temporal_recall.mean()),
            mean_called_kyr=float(sel.called_kyr.mean()),mean_iou=float(sel.iou.mean()),
            center_peak_mae_kya=float(sel.center_peak_error_kya.mean()),
            event_tp=tp,event_fp=fp,event_fn=fn,event_f1=2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else None))
    frame=pd.DataFrame(rows);frame.to_csv(out/"replicate_metrics.tsv",sep="\t",index=False,na_rep="NA")
    cumulative_baselines(records,cdfs,positions).to_csv(out/"cumulative_baselines.tsv",sep="\t",index=False)
    records.to_csv(out/"completed_units.tsv",sep="\t",index=False)
    descriptions=[]
    mid=len(positions)//2
    for j,method in enumerate(METHODS):
        for group in ("development","calibration","neutral_test","selected"):
            subset=cdfs[cohort==group,j,mid]
            for kya in (1,10,20,30,40,50):
                descriptions.append(dict(method=method,cohort=group,summary=f"below_{kya}_kya",mean_fraction=float(subset[:,kya-1].mean())))
            descriptions.append(dict(method=method,cohort=group,summary="within_20_30_kya",mean_fraction=float((subset[:,29]-subset[:,19]).mean())))
    pd.DataFrame(descriptions).to_csv(out/"center_coalescence_summaries.tsv",sep="\t",index=False)
    np.savez_compressed(out/"calibration.npz",center=np.stack(centers),scale=np.stack(scales),maxima=np.stack(maxima))
    # Paired bootstrap over test regions, never pairs or grid cells.
    rng=np.random.default_rng(20380101)
    comparison={}
    for metric in ("any_call","joint_hit","center_overlap","temporal_f1"):
        table=frame[frame.cohort=="selected"].pivot(index="unit",columns="method",values=metric)
        delta=(table["mass"]-table["median"]).to_numpy()
        boot=delta[rng.integers(0,len(delta),size=(2000,len(delta)))].mean(axis=1)
        comparison[metric]=dict(mass_minus_median=float(delta.mean()),ci95=np.quantile(boot,[.025,.975]).tolist())
    status=json.loads((root/"status.json").read_text()) if (root/"status.json").exists() else []
    results=dict(config=manifest["config"],window_widths_kyr=list(widths),cohort_counts=records.cohort.value_counts().to_dict(),
        scientific_parameters=dict(generation_time_years=25,mutation_rate=1.29e-8,recombination_rate=1e-8,
            selection_start_years=30000,selection_end_years=20000,dominance=.5,standing_frequency_target=.05,
            standing_frequency_bounds=[.02,.10],no_introgression=True,spatial_tolerance_bp=100000,
            decoder_scaled_mutation_rate=.00075,decoder_threads_per_worker=1,decoder_cache_size=1000,
            decoder_exp10="accurate",decoder_backward_alignment="fixed",alpha=.05,bootstrap_draws=2000,
            phlash_resource_sha256="a5bdfd843629d48050cd56a84da4d70e356d0339017d3c1f67be646d727db928"),
        statuses={state:sum(r["state"]==state for r in status) for state in sorted(set(r["state"] for r in status))},
        summaries=summaries,paired_comparison=comparison,root=str(root),reference_null_root=str(reference) if reference else None,
        selection_frequency_summary=records[records.selected].final_frequency.describe().to_dict(),
        notes=["Independent calibration and neutral test regions; no time/pair pseudoreplication.",
        f"Test is maximum over all genomic positions and {len(intervals)} time intervals; widths={list(widths)} kyr.",
        "Event matching requires a cell within 100 kb of the selected center and in [20,30) kya.",
        "Temporal F1 is computed at the selected coordinate, with no calls scored zero.",
        "This is a 1-Mb bounded pilot, not the planned 10-Mb/10,000-pair benchmark.",
        "Standing allele chosen on a valid onset genealogy near 5% frequency; no final-frequency conditioning.",
        "Analytic Gamma mass is a marginal posterior expectation; full-region calibration handles dependence."])
    (out/"performance.json").write_text(json.dumps(results,indent=2,allow_nan=False)+"\n")
    plot(out,records,cdfs,positions,summaries,frame,widths)
    lines=["# Bounded no-introgression EAS time-scan pilot", "", f"Completed 2026-11-02. Independent whole-region calibration at alpha=0.05; time-window widths {list(widths)} kyr.","",f"Cohort counts: {results['cohort_counts']}","", "| Method | Neutral FPR | Any-call power | Joint space/time TPR | Center overlap TPR | Mean temporal F1 | Event F1 |", "|---|---:|---:|---:|---:|---:|---:|"]
    for s in summaries:
        rate=lambda key:f"{s['rates'][key]['k']}/{s['rates'][key]['n']} ({s['rates'][key]['rate']:.1%})"
        lines.append(f"| {s['method']} | {rate('neutral_fpr')} | {rate('selected_any_call')} | {rate('joint_tpr')} | {rate('center_overlap_tpr')} | {s['mean_temporal_f1']:.3f} | {s['event_f1']:.3f} |")
    lines.extend(["","## Design and limitations",""]+["- "+n for n in results["notes"]])
    lines.extend(["",f"Settings: {manifest['config']}","","Rates, Wilson 95% intervals, paired differences and uncertainty are in `performance.json`. Raw simulation trees/posteriors remain in the local campaign directory; only aggregate simulation results and provenance are versioned."])
    (out/"README.md").write_text("\n".join(lines)+"\n")
    artifacts={p.name:dict(bytes=p.stat().st_size,sha256=sha(p)) for p in out.iterdir() if p.is_file() and p.name!="manifest.json"}
    (out/"manifest.json").write_text(json.dumps(dict(source_manifest=sha(root/"manifest.json"),reference_manifest=sha(reference/"manifest.json") if reference else None,report_source=sha(__file__),artifacts=artifacts),indent=2)+"\n")
    print(json.dumps(results,indent=2))


def plot(out,records,cdfs,positions,summaries,frame,widths):
    import matplotlib
    matplotlib.use("Agg")
    matplotlib.rcParams.update({"font.size":16,"axes.labelsize":17,"axes.titlesize":18,"pdf.fonttype":42,"ps.fonttype":42})
    import matplotlib.pyplot as plt
    colors=["#444444","#d97706","#0072b2"]
    fig,axes=plt.subplots(2,2,figsize=(14,10),layout="constrained")
    fig.suptitle(f"EAS no-introgression pilot | window widths: {', '.join(map(str,widths))} kyr",fontsize=20)
    sel=records.selected.to_numpy(dtype=bool);null=~sel
    mid=len(positions)//2
    for j,name in enumerate(METHODS):
        difference=np.diff(cdfs[sel,j,mid].mean(axis=0)-cdfs[null,j,mid].mean(axis=0))
        axes[0,0].plot(np.arange(1.5,50),100*difference,label=name,color=colors[j],lw=2)
    axes[0,0].axvspan(20,30,color="gray",alpha=.15)
    axes[0,0].axhline(0,color="black",lw=.5)
    axes[0,0].set(xlabel="TMRCA (kya)",ylabel="Excess mass (percentage points / kyr)",title="Mean excess at selected coordinate")
    axes[0,0].legend(loc="lower center",bbox_to_anchor=(.5,1.10),ncols=3,frameon=False,fontsize=13)
    x=np.arange(3)
    for k,(metric,label_) in enumerate((("neutral_fpr","Neutral FPR"),("joint_tpr","Joint TPR"))):
        vals=np.array([s["rates"][metric]["rate"] for s in summaries])
        bounds=np.array([s["rates"][metric]["ci95"] for s in summaries])
        axes[0,1].errorbar(x+(k-.5)*.18,vals,yerr=[vals-bounds[:,0],bounds[:,1]-vals],fmt="o",capsize=4,label=label_)
    axes[0,1].set(xticks=x,xticklabels=METHODS,ylim=(-.03,1.03),ylabel="Fraction of test regions",title="Held-out scan performance")
    axes[0,1].legend(loc="lower center",bbox_to_anchor=(.5,1.10),ncols=2,frameon=False,fontsize=13)
    axes[1,0].bar(x,[s["mean_temporal_f1"] for s in summaries],color=colors)
    for i,s in enumerate(summaries):
        axes[1,0].text(i,s["mean_temporal_f1"]+.025,f"{s['mean_temporal_f1']:.3f}",ha="center",fontsize=15)
    axes[1,0].set(xticks=x,xticklabels=METHODS,ylim=(0,1),ylabel="Mean temporal F1",title="Interval recovery at selected coordinate")
    axes[1,1].hist(records[sel].final_frequency,bins=np.linspace(0,1,11),color="#777777",edgecolor="white")
    axes[1,1].set(xlabel="Final selected-allele frequency",ylabel="Independent selected regions",title="All selected outcomes retained")
    for ext in ("png","pdf"): fig.savefig(out/f"performance.{ext}",dpi=160)
    plt.close(fig)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root",type=Path,required=True);parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--widths",type=int,nargs="+",default=[1])
    parser.add_argument("--reference-nulls",type=Path)
    args=parser.parse_args();report(args.root,args.out,args.widths,args.reference_nulls)
