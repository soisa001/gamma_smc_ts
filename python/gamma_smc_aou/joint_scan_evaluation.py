"""Whole-region held-out evaluation of archaic pair-class confirmation.

Five outer folds: 400 neutral fit, 400 neutral calibration, 200 neutral test;
80 selected training and 20 selected test per onset. The whole site/time/run
search is reduced to one region score before calibration. No independent-pair
assumption, no causal-coordinate lookup, and no dropping undetectable regions.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path

os.environ["MPLBACKEND"]="Agg"
for _n in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"):
    os.environ[_n]="1"
import numpy as np
import pandas as pd
from scipy.stats import norm
from .fresh_power import atomic_json, digest

GATES=(0.,.5,.8,.95)
PRIMARY="joint_g50_r1_T50000"


def methods(cutoffs):
    result=[]
    for run in (1,2,3):
        result.append(dict(name=f"af_r{run}",family="af",gate=0.,run=run,cutoff=-1,years=0))
        for j,t in enumerate(cutoffs):
            result.append(dict(name=f"all_r{run}_T{t}",family="all",gate=0.,run=run,cutoff=j,years=int(t)))
            for gate in GATES:
                result.append(dict(name=f"joint_g{int(gate*100)}_r{run}_T{t}",family="joint",gate=gate,run=run,cutoff=j,years=int(t)))
        # Union across all six cutoffs is itself calibrated as one region score.
        for gate in GATES:
            result.append(dict(name=f"joint_g{int(gate*100)}_r{run}_Tall",family="joint",gate=gate,run=run,cutoff=-1,years=-1))
    return result


def components(counts,n):
    fractions=np.divide(counts,n[None,:,:,None],out=np.zeros_like(counts,dtype=float),where=n[None,:,:,None]>0)
    rr,ar,aa,all_pairs=[fractions[:,:,j] for j in range(4)]
    weight=n[:,2]/n[:,3]
    # Empty RR/AR classes at fixation use zero comparison baselines. This
    # reduces to the all-pair recent fraction, retaining fixed sweeps.
    joint=weight[None,:,None]*np.maximum(aa-rr,0)*(1-ar)
    return all_pairs,joint


def binned_region_score(local,positions,nbins,stride,run):
    local=np.asarray(local)
    binned=np.zeros(nbins,dtype=float)
    if len(local):
        np.maximum.at(binned,np.asarray(positions,dtype=int)//stride,local)
    if run==1:return float(binned.max())
    return float(np.lib.stride_tricks.sliding_window_view(binned,run).min(axis=1).max())


def score_region(features,scheme,source,method,gate_values,stride,nbins):
    prefix="grid" if scheme=="grid10kb" else "site"
    counts,n=features[prefix+"_counts"],features[prefix+"_n"]
    positions=features["grid"] if prefix=="grid" else features["sites"]
    all_pairs,joint=components(counts,n)
    if method["family"]=="af":
        if prefix=="site": local=features["site_af"]
        else:
            idx=features["grid_markers"]
            local=np.zeros(len(idx))
            local[idx>=0]=features["site_af"][idx[idx>=0]]
    elif method["family"]=="all":
        local=all_pairs[source,:,method["cutoff"]]
    else:
        values=joint[source].copy()
        if method["gate"]:
            # Gate at the SAME time cutoff and coordinate as the carrier score.
            values[all_pairs[source]<=gate_values[method["gate"]][source]]=0
        local=values.max(axis=1) if method["cutoff"]<0 else values[:,method["cutoff"]]
    return binned_region_score(local,positions,nbins,stride,method["run"])


def rank_p(calibration,scores):
    calibration=np.sort(calibration)
    return (1+len(calibration)-np.searchsorted(calibration,scores,side="left"))/(len(calibration)+1)


def make_folds(items,seed):
    rng=np.random.default_rng(seed)
    result=np.zeros(len(items),dtype=int)
    for onset in (0,10000,50000):
        group=sorted([i for i,x in enumerate(items) if x["onset_years"]==onset],key=lambda i:items[i]["replicate"])
        result[rng.permutation(group)]=np.arange(len(group))%5
    return result


def wilson(k,n):
    z=norm.ppf(.975)
    p=k/n
    center=(p+z*z/(2*n))/(1+z*z/n)
    half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/(1+z*z/n)
    return max(0,center-half),min(1,center+half)


def metrics(predictions):
    rows=[]
    keys=["scheme","source","method","alpha"]
    for key,group in predictions.groupby(keys,sort=True):
        neutral=group[group.onset_years==0]
        for onset in (10000,50000):
            selected=group[group.onset_years==onset]
            if not len(selected) or not len(neutral):continue
            tp,fp=int(selected.called.sum()),int(neutral.called.sum())
            sensitivity,fpr=tp/len(selected),fp/len(neutral)
            row=dict(zip(keys,key),onset_years=onset,selected_called=tp,selected_n=len(selected),
                     neutral_called=fp,neutral_n=len(neutral),power=sensitivity,neutral_call_fraction=fpr)
            row["power_low"],row["power_high"]=wilson(tp,len(selected))
            row["neutral_low"],row["neutral_high"]=wilson(fp,len(neutral))
            for prevalence in (.5,.1,.01):
                den=prevalence*sensitivity+(1-prevalence)*fpr
                row[f"fdr_prevalence_{prevalence}"]=(1-prevalence)*fpr/den if den else np.nan
                row[f"f1_prevalence_{prevalence}"]=2*prevalence*sensitivity/(prevalence*(1+sensitivity)+(1-prevalence)*fpr) if (sensitivity+fpr)>0 else 0.
            rows.append(row)
    return pd.DataFrame(rows)


def load_one(payload):
    item,out=payload
    dest=out/"profiles"/item["key"]
    receipt=json.loads((dest/"complete.json").read_text())
    if receipt["identity"]["simulation"]!=item["receipt_sha256"]:raise ValueError("Profile simulation mismatch")
    spec=receipt["outputs"]["features.npz"]
    if (dest/"features.npz").stat().st_size!=spec["bytes"] or digest(dest/"features.npz")!=spec["sha256"]:
        raise ValueError("Corrupt scan features")
    with np.load(dest/"features.npz",allow_pickle=False) as saved:
        features={k:saved[k] for k in saved.files}
    return features,receipt,spec["sha256"]


def evaluate(out):
    manifest=json.loads((out/"profile_manifest.json").read_text())
    items,cfg=manifest["items"],manifest["config"]
    if len(items)!=1200:raise ValueError("Final evaluation requires all 1200 independent regions")
    with ThreadPoolExecutor(max_workers=20) as pool:
        loaded=list(pool.map(load_one,[(i,out) for i in items]))
    features=[v[0] for v in loaded]
    methods_list=methods(cfg["tmrca_cutoffs_years"])
    folds=make_folds(items,cfg["seed_base"])
    onset=np.array([i["onset_years"] for i in items])
    neutral=onset==0
    output=out/"analysis";output.mkdir(exist_ok=True)
    all_grid=np.array([components(f["grid_counts"],f["grid_n"])[0] for f in features])
    # shape regions x sources x 1000 strides x 6 cutoffs
    predictions,chosen,targets,gate_records=[],[],[],[]
    score_archive={}
    for fold in range(5):
        test=folds==fold
        fit=neutral & np.isin(folds,[(fold+1)%5,(fold+2)%5])
        calibration=neutral & np.isin(folds,[(fold+3)%5,(fold+4)%5])
        if (fit.sum(),calibration.sum(),(test&neutral).sum())!=(400,400,200):
            raise ValueError("Fold sizes violate calibration contract")
        # Pooled neutral strides define a translation-invariant prescreen.
        gate_values={g:np.quantile(all_grid[fit],g,axis=(0,2)) for g in GATES if g}
        gate_records.append(dict(fold=fold,fit_ids=[items[i]["key"] for i in np.where(fit)[0]],
                                 calibration_ids=[items[i]["key"] for i in np.where(calibration)[0]],
                                 values={str(g):v.tolist() for g,v in gate_values.items()}))
        for scheme in ("grid10kb","archaic_sites"):
            for source_index,source in enumerate(("truth","decoded")):
                scores=np.empty((len(items),len(methods_list)))
                # Components are computed once per region below for efficiency.
                for i,f in enumerate(features):
                    prefix="grid" if scheme=="grid10kb" else "site"
                    positions=f["grid"] if prefix=="grid" else f["sites"]
                    all_pairs,joint=components(f[prefix+"_counts"],f[prefix+"_n"])
                    af=f["site_af"] if prefix=="site" else np.zeros(len(positions))
                    if prefix=="grid":
                        idx=f["grid_markers"];af[idx>=0]=f["site_af"][idx[idx>=0]]
                    gated={0.:joint[source_index]}
                    for g in GATES[1:]:
                        gated[g]=np.where(all_pairs[source_index]>gate_values[g][source_index],joint[source_index],0.)
                    nbins=cfg["scored_length_bp"]//cfg["stride_bp"]
                    # Bin once per local feature, then evaluate run lengths.
                    local_features={"af":af}
                    for j,t in enumerate(cfg["tmrca_cutoffs_years"]):
                        local_features[f"all_T{t}"]=all_pairs[source_index,:,j]
                        for g in GATES:local_features[f"joint_g{int(g*100)}_T{t}"]=gated[g][:,j]
                    for g in GATES:local_features[f"joint_g{int(g*100)}_Tall"]=gated[g].max(axis=1)
                    cached={}
                    for name,local in local_features.items():
                        binned=np.zeros(nbins)
                        np.maximum.at(binned,positions//cfg["stride_bp"],local)
                        cached[(name,1)]=float(binned.max())
                        for run in (2,3):cached[(name,run)]=float(np.lib.stride_tricks.sliding_window_view(binned,run).min(axis=1).max())
                    for k,m in enumerate(methods_list):
                        name=m["name"].replace(f"_r{m['run']}","")
                        scores[i,k]=cached[(name,m["run"])]
                score_archive[f"fold{fold}_{scheme}_{source}"]=scores
                selected_train=(~neutral)&(~test)
                joint_index=np.array([k for k,m in enumerate(methods_list) if m["family"]=="joint"])
                # Family/cutoff/run selection uses FIT regions only. Calibration
                # and test regions never participate in choosing the method.
                training_p=np.column_stack([rank_p(scores[fit,k],scores[selected_train,k]) for k in joint_index])
                train_power=(training_p<=.05).mean(axis=0)
                best=int(joint_index[np.argmax(train_power)])
                chosen.append(dict(fold=fold,scheme=scheme,source=source,method=methods_list[best]["name"],training_power=float(train_power.max())))
                selected_methods=list(enumerate(methods_list))+[(best,dict(name="tuned_joint"))]
                for k,m in selected_methods:
                    ids=np.where(test)[0]
                    pvalues=rank_p(scores[calibration,k],scores[test,k])
                    for alpha in (.01,.05,.1):
                        for i,pv in zip(ids,pvalues):
                            predictions.append(dict(key=items[i]["key"],onset_years=int(onset[i]),fold=fold,
                                scheme=scheme,source=source,method=m["name"],alpha=alpha,score=float(scores[i,k]),p=float(pv),called=bool(pv<=alpha)))
                    if m["name"] in (PRIMARY,"tuned_joint"):
                        for target_onset in (10000,50000):
                            train=(onset==target_onset)&(~test)
                            threshold=float(np.quantile(scores[train,k],.3,method="lower"))
                            for i in np.where(test & ((onset==target_onset)|neutral))[0]:
                                targets.append(dict(key=items[i]["key"],onset_years=int(onset[i]),target_onset=target_onset,
                                    fold=fold,scheme=scheme,source=source,method=m["name"],threshold=threshold,
                                    score=float(scores[i,k]),called=bool(scores[i,k]>=threshold and scores[i,k]>0)))
                print(json.dumps(dict(fold=fold,scheme=scheme,source=source,status="scored")),flush=True)
    pred=pd.DataFrame(predictions)
    pred.to_csv(output/"heldout_predictions.csv.gz",index=False,compression=dict(method="gzip",mtime=0))
    table=metrics(pred);table.to_csv(output/"heldout_metrics.csv",index=False)
    pd.DataFrame(chosen).to_csv(output/"training_choices.csv",index=False)
    target=pd.DataFrame(targets)
    target.to_csv(output/"target70_predictions.csv",index=False)
    target_tables=[]
    for target_onset,group in target.groupby("target_onset"):
        group=group.assign(alpha=.7)
        part=metrics(group)
        target_tables.append(part.assign(target_onset=target_onset))
    target_metrics=pd.concat(target_tables,ignore_index=True)
    target_metrics.to_csv(output/"target70_metrics.csv",index=False)
    np.savez_compressed(output/"fold_region_scores.npz",**score_archive)
    atomic_json(output/"gate_fits.json",gate_records)
    atomic_json(output/"methods.json",methods_list)
    pd.DataFrame([dict(key=i["key"],onset_years=i["onset_years"],fold=int(folds[k]),
                       markers=loaded[k][1]["markers"],valid_grid_positions=loaded[k][1]["valid_grid_positions"],
                       **loaded[k][1]["verification"]) for k,i in enumerate(items)]).to_csv(output/"region_inventory.csv",index=False)
    # Fixed, ungated scores need no fitting. At p=.001 the complete 1000-null
    # reference can resolve a selected score exceeding all neutral scores.
    # Neutral LOO rank counts are descriptive, NOT an independent validation.
    extreme=[]
    for scheme in ("grid10kb","archaic_sites"):
        for source in ("truth","decoded"):
            scores=score_archive[f"fold0_{scheme}_{source}"]
            for k,m in enumerate(methods_list):
                if m["gate"]!=0:continue
                null=scores[neutral,k]
                for target_onset in (10000,50000):
                    selected=scores[onset==target_onset,k]
                    called=int(np.sum(rank_p(null,selected)<=.001))
                    loo=(np.sum(null[None,:]>=null[:,None],axis=1))/len(null)
                    extreme.append(dict(scheme=scheme,source=source,method=m["name"],onset_years=target_onset,
                        selected_called=called,selected_n=len(selected),power=called/len(selected),
                        neutral_loo_called=int(np.sum(loo<=.001)),neutral_n=len(null),
                        independent_neutral_validation=False))
    pd.DataFrame(extreme).to_csv(output/"p001_full_null_exploratory.csv",index=False)
    atomic_json(output/"provenance.json",dict(config=cfg,profile_manifest_sha256=digest(out/"profile_manifest.json"),
        source_sha256=digest(Path(__file__)),profile_hashes={i["key"]:loaded[k][2] for k,i in enumerate(items)},
        primary_method=PRIMARY,seed=cfg["seed_base"],fold_sizes=dict(fit_neutral=400,calibration_neutral=400,test_neutral=200,train_selected_per_onset=80,test_selected_per_onset=20),
        minimum_crossfit_p=1/401,p001_note="Exploratory fixed-score calibration against all 1000 neutrals; neutral LOO ranks are constrained by construction and not independent validation",
        operational_fdr="Fraction of whole neutral 10 Mb regions with any call (a regional false-positive rate)",
        selected_endpoint="Any call anywhere in region; all 100 surviving-allele regions per onset retained",
        primary_frozen_before_outcomes=True,methods_are_exploratory_except_primary=True,
        confidence_interval_note="Wilson intervals summarize held-out binomial counts; shared fitted thresholds induce dependence not included in these intervals"))
    figures(table,target_metrics,output)
    atomic_json(output/"artifact_manifest.json",{p.name:dict(bytes=p.stat().st_size,sha256=digest(p)) for p in sorted(output.iterdir()) if p.is_file() and p.name!="artifact_manifest.json"})
    print(table[(table.method.isin([PRIMARY,"tuned_joint","af_r1","all_r1_T50000"]))&(table.alpha==.05)][["scheme","source","method","onset_years","power","neutral_call_fraction"]].to_string(index=False),flush=True)


def figures(table,target,output):
    import matplotlib
    matplotlib.rcParams.update({"pdf.fonttype":42,"ps.fonttype":42,"font.size":17,"axes.titlesize":18,"axes.labelsize":18,"legend.fontsize":12})
    import matplotlib.pyplot as plt
    selected_methods=["all_r1_T50000","af_r1",PRIMARY,"tuned_joint"]
    names=["All pairs, 50 ky","Archaic AF","Joint, 50 ky","Trained joint"]
    colors=["#666666","#CC8833","#228833","#4477AA"]
    checks=[]
    for scheme in ("grid10kb","archaic_sites"):
        fig,axes=plt.subplots(2,2,figsize=(11,8.5),sharex=True,sharey=True)
        fig.subplots_adjust(left=.105,right=.97,bottom=.12,top=.79,hspace=.30,wspace=.20)
        for row,source in enumerate(("truth","decoded")):
            for col,onset in enumerate((50000,10000)):
                ax=axes[row,col]
                for m,label,color in zip(selected_methods,names,colors):
                    sub=table[(table.scheme==scheme)&(table.source==source)&(table.onset_years==onset)&(table.method==m)].sort_values("neutral_call_fraction")
                    ax.plot(sub.neutral_call_fraction*100,sub.power*100,"o-",color=color,label=label,lw=2)
                ax.axhline(70,color="#999999",ls=":",lw=1)
                ax.axvline(5,color="#999999",ls=":",lw=1)
                ax.set_title(f"{source.title()}, onset {onset//1000} ky")
                ax.set(xlim=(0,15),ylim=(0,102))
                ax.grid(alpha=.15)
                if col==0:ax.set_ylabel("Selected regions called (%)")
                if row==1:ax.set_xlabel("Neutral regions called (%)")
        handles,labels=axes[0,0].get_legend_handles_labels()
        legend=fig.legend(handles,labels,loc="upper center",bbox_to_anchor=(.52,.915),ncol=4,frameon=False)
        fig.suptitle("10 kb grid" if scheme=="grid10kb" else "Every observed archaic SNP",y=.98,fontsize=21)
        fig.canvas.draw()
        renderer=fig.canvas.get_renderer()
        box=legend.get_window_extent(renderer)
        if any(box.overlaps(ax.get_window_extent(renderer)) for ax in axes.flat):
            raise ValueError("Figure legend overlaps a data panel")
        for ax in axes.flat:
            tight=ax.get_tightbbox(renderer)
            if tight.x0<0 or tight.y0<0 or tight.x1>fig.bbox.x1+1 or tight.y1>fig.bbox.y1+1:
                raise ValueError("Figure labels exceed canvas")
        name=f"power_neutral_calls_{scheme}"
        fig.savefig(output/f"{name}.png",dpi=180)
        fig.savefig(output/f"{name}.pdf",metadata={"CreationDate":None,"ModDate":None})
        if b"/Subtype /Image" in (output/f"{name}.pdf").read_bytes():
            raise ValueError("Scientific PDF was rasterized")
        checks.append(dict(name=name,legend_clear=True,labels_within_canvas=True,vector_pdf=True))
        plt.close(fig)
    atomic_json(output/"figure_checks.json",checks)


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args(argv)
    evaluate(args.out)


if __name__=="__main__":main()
