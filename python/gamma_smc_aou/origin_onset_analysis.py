"""Focal power and independent neutral-target FPR for the origin/onset pilot."""
from __future__ import annotations

import os
os.environ["MPLBACKEND"] = "Agg"
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 17,
    "axes.titlesize": 18, "axes.labelsize": 17, "xtick.labelsize": 13, "ytick.labelsize": 14})
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import origin_onset as oo

PAIR_LABELS = {"all pairs": "all_pairs", "alt/alt": "alt_alt"}
ARMS = ("I50", "I10", "D50", "D10")


def arms(cfg):
    return tuple(f'I{t//1000}' for t in cfg['selection_onsets_years']) if cfg.get('focal_ascertainment') else ARMS


def family_for_arm(arm, cfg):
    return arm if cfg.get('focal_ascertainment') else ('I50' if arm.startswith('I') else arm)


def rank_p(null, scores):
    null, scores = np.asarray(null, dtype=float), np.asarray(scores, dtype=float)
    if len(null) == 0 or np.isinf(null).any() or np.isinf(scores).any():
        raise ValueError("Invalid empirical-rank inputs")
    # Missing pair classes are unavailable/no-call outcomes, not raw zero scores.
    reference = np.sort(np.where(np.isnan(null), -np.inf, null))
    p = (1 + len(null) - np.searchsorted(reference, scores, side="left")) / (len(null) + 1)
    return np.where(np.isnan(scores), 1.0, p)


def wilson(k, n):
    if n == 0:
        return np.nan, np.nan
    z = 1.959963984540054
    p = k/n
    center = (p + z*z/(2*n))/(1+z*z/n)
    half = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n))/(1+z*z/n)
    return center-half, center+half


def evaluate(focal, cfg):
    """Same 1,000 reference IDs for both target roles, separately by source/score/T."""
    predictions = []
    grouping = ["family", "source", "method", "cutoff_years"]
    for key, group in focal.groupby(grouping, dropna=False):
        null = group[group.role == "null"]
        targets = group[group.role != "null"].copy()
        if len(null) != cfg["null_replicates"] or null.task_id.duplicated().any():
            raise ValueError(f"Wrong null cohort: {key}")
        if set(null.task_id) & set(targets.task_id):
            raise ValueError("A test target entered the calibration cohort")
        targets["p"] = rank_p(null.score, targets.score)
        targets["called"] = targets.p <= cfg["alpha"]
        targets["evaluable"] = targets.score.notna()
        targets["null_n"] = len(null)
        targets["null_evaluable"] = int(null.score.notna().sum())
        predictions.append(targets)
    predictions = pd.concat(predictions, ignore_index=True)
    rows = []
    grouping = ["family", "arm", "role", "s", "source", "method", "cutoff_years"]
    for key, group in predictions.groupby(grouping, dropna=False):
        n, calls = len(group), int(group.called.sum())
        if n != cfg["target_replicates"]:
            raise ValueError(f"Incomplete target cell {key}: {n}")
        lo, hi = wilson(calls, n)
        evaluable = int(group.evaluable.sum())
        rows.append(dict(zip(grouping,key), n=n, calls=calls, rate=calls/n,
            ci_low=lo, ci_high=hi, evaluable=evaluable, coverage=evaluable/n,
            evaluable_rate=calls/evaluable if evaluable else np.nan,
            null_n=int(group.null_n.iloc[0]), null_evaluable=int(group.null_evaluable.iloc[0]),
            alpha=cfg["alpha"], endpoint="focal power" if key[2]=="selected" else "neutral detection (FPR)"))
    return predictions, pd.DataFrame(rows)


def save_figure(fig, directory, name):
    directory.mkdir(parents=True, exist_ok=True)
    # Reserve distinct title and horizontal-legend bands above every data panel.
    extras=[]
    if fig.legends:
        has_title=fig._suptitle is not None
        fig.set_layout_engine("constrained",rect=(0,0,1,.83 if has_title else .91))
        if has_title:
            fig._suptitle.set_y(.99)
            fig._suptitle.set_in_layout(False)
            extras.append(fig._suptitle)
        for legend in fig.legends:
            legend.set_loc("upper center")
            legend.set_bbox_to_anchor((.5,.90 if has_title else .995))
            legend.set_in_layout(False)
            extras.append(legend)
    for ax in fig.axes:
        for collection in ax.collections:
            collection.set_rasterized(False)
    fig.savefig(directory / (name + ".png"), dpi=180, bbox_inches="tight",bbox_extra_artists=extras or None)
    fig.savefig(directory / (name + ".pdf"), bbox_inches="tight",bbox_extra_artists=extras or None)
    plt.close(fig)


def heatmap(ax, values, xlabels, ylabels, maximum, title):
    mesh = ax.pcolormesh(np.ma.masked_invalid(values), vmin=0, vmax=maximum,
                         cmap="viridis", edgecolors="white", linewidth=.8, rasterized=False)
    ax.set_xticks(np.arange(len(xlabels))+.5, xlabels)
    ax.set_yticks(np.arange(len(ylabels))+.5, ylabels)
    ax.invert_yaxis()
    ax.set_title(title)
    ax.set_xlabel("TMRCA cutoff (kya)")
    for r,c in np.ndindex(values.shape):
        value=values[r,c]
        ax.text(c+.5,r+.5,"NA" if np.isnan(value) else f"{value:.0f}",ha="center",va="center",
                fontsize=13,color="white" if np.isnan(value) or value<maximum*.55 else "black")
    return mesh


def plot_results(out, cfg, focal, metrics):
    directory=out/"figures"
    cutoffs=cfg["tmrca_cutoffs_years"]
    coefficients=cfg["selection_coefficients"]
    s_labels=[f"{s:g}" for s in coefficients]
    selected_arms = arms(cfg)
    families = list(dict.fromkeys(family_for_arm(arm,cfg) for arm in selected_arms))
    layout = (1,2) if len(selected_arms)==2 else (2,2)
    figsize = (14,7) if len(selected_arms)==2 else (14,11)
    neutral=metrics[(metrics.role=="neutral_target") & (metrics.method!="AF")]
    fpr_max=max(10.,float(np.ceil(neutral.rate.max()*100/10)*10))
    for source in ("truth","decoded"):
        for method in PAIR_LABELS.values():
            selected=metrics[(metrics.source==source)&(metrics.method==method)&(metrics.role=="selected")]
            fig,axes=plt.subplots(*layout,figsize=figsize,layout="constrained")
            for ax,arm in zip(axes.flat,selected_arms):
                z=selected[selected.arm==arm].pivot(index="s",columns="cutoff_years",values="rate")
                values=z.reindex(index=coefficients,columns=cutoffs).to_numpy()*100
                mesh=heatmap(ax,values,[t//1000 for t in cutoffs],s_labels,100,arm)
                ax.set_ylabel("Selection coefficient s")
            bar=fig.colorbar(mesh,ax=axes,fraction=.025,pad=.025,label="Focal detection (%)")
            bar.solids.set_rasterized(False)
            fig.suptitle(f"{source.capitalize()} TMRCA | {method.replace('_',' ')} | p <= 0.05\n{cfg['target_replicates']} targets per cell; {cfg['null_replicates']} matched nulls",fontsize=20)
            save_figure(fig,directory,f"power_{source}_{method}")

            z=neutral[(neutral.source==source)&(neutral.method==method)].pivot(index="family",columns="cutoff_years",values="rate")
            values=z.reindex(index=families,columns=cutoffs).to_numpy()*100
            fig,ax=plt.subplots(figsize=(11,7),layout="constrained")
            mesh=heatmap(ax,values,[t//1000 for t in cutoffs],families,fpr_max,"Independent neutral targets")
            bar=fig.colorbar(mesh,ax=ax,fraction=.025,pad=.025,label="Neutral detection / FPR (%)")
            bar.solids.set_rasterized(False)
            fig.suptitle(f"{source.capitalize()} TMRCA | {method.replace('_',' ')}\n{cfg['target_replicates']} neutral targets; {cfg['null_replicates']} separate calibration nulls",fontsize=19)
            save_figure(fig,directory,f"fpr_{source}_{method}")

    af=focal[focal.method=="AF"]
    fig,axes=plt.subplots(*layout,figsize=figsize,layout="constrained")
    for ax,arm in zip(axes.flat,selected_arms):
        family=family_for_arm(arm,cfg)
        datasets=[af[(af.family==family)&(af.role=="neutral_target")].score.to_numpy()]
        datasets += [af[(af.arm==arm)&(af.role=="selected")&(af.s==s)].score.to_numpy() for s in coefficients]
        ax.boxplot(datasets,tick_labels=["0"]+s_labels,showfliers=True)
        null=af[(af.family==family)&(af.role=="null")].score
        lo,median,hi=null.quantile([.25,.5,.75])
        ax.axhspan(lo,hi,color="grey",alpha=.18,label="Null IQR")
        ax.axhline(median,color="grey",ls="--",label="Null median")
        ax.set(ylim=(0,1.03),title=arm,xlabel="Selection coefficient (0 = neutral targets)",ylabel="Focal ALT frequency")
        ax.tick_params(axis="x",rotation=60)
    handles,labels=axes.flat[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc="outside upper center",ncol=2,frameon=False)
    save_figure(fig,directory,"allele_frequency")

    # One cutoff per page keeps all 12 coefficients legible at letter-page scale.
    for arm in selected_arms:
        family=family_for_arm(arm,cfg)
        for cutoff in cutoffs:
            fig,axes=plt.subplots(2,2,figsize=(14,11),layout="constrained")
            for row,source in enumerate(("truth","decoded")):
                for col,method in enumerate(PAIR_LABELS.values()):
                    ax=axes[row,col]
                    z=focal[(focal.source==source)&(focal.method==method)&(focal.cutoff_years==cutoff)]
                    datasets=[z[(z.family==family)&(z.role=="neutral_target")].score.dropna().to_numpy()]
                    datasets += [z[(z.arm==arm)&(z.role=="selected")&(z.s==s)].score.dropna().to_numpy() for s in coefficients]
                    ax.boxplot(datasets,tick_labels=["0"]+s_labels)
                    null=z[(z.family==family)&(z.role=="null")].score.dropna()
                    if len(null):
                        lo,median,hi=null.quantile([.25,.5,.75])
                        ax.axhspan(lo,hi,color="grey",alpha=.18,label="Null IQR")
                        ax.axhline(median,color="grey",ls="--",label="Null median")
                    ax.set(ylim=(-.03,1.03),title=f"{source} | {method.replace('_',' ')}",xlabel="s (0 = neutral targets)",ylabel="Raw frac_recent_T")
                    ax.tick_params(axis="x",rotation=60)
            fig.suptitle(f"{arm}: focal raw fractions at T = {cutoff//1000} kya\nUnavailable pair classes omitted here; coverage is reported in tables",fontsize=19)
            handles,labels=axes.flat[0].get_legend_handles_labels()
            if handles:
                fig.legend(handles,labels,loc="outside upper center",ncol=2,frameon=False)
            save_figure(fig,directory,f"raw_{arm}_T{cutoff}")

    af_metrics=metrics[(metrics.method=="AF")&(metrics.role=="selected")]
    fig,ax=plt.subplots(figsize=(11,8.5),layout="constrained")
    for arm in selected_arms:
        z=af_metrics[af_metrics.arm==arm].sort_values("s")
        ax.errorbar(z.s,z.rate*100,yerr=np.array([z.rate-z.ci_low,z.ci_high-z.rate])*100,marker="o",capsize=3,label=arm)
    ax.set(xlabel="Selection coefficient s",ylabel="AF-only focal detection (%)",ylim=(-3,103),title="AF-only power at p <= 0.05")
    fig.legend(loc="outside upper center",ncol=4,frameon=False)
    save_figure(fig,directory,"af_only_power")
    oo.legacy.atomic_json(out/"plot_scales.json",dict(power_percent=[0,100],neutral_fpr_percent=[0,fpr_max],raw_fraction=[0,1],af=[0,1]))


def example_profiles(out, cfg, targets):
    records=[]
    group_keys=["family","arm","role","s"]
    for key,group in targets.groupby(group_keys):
        median=group.sample_af.median()
        chosen=group.assign(distance=(group.sample_af-median).abs()).sort_values(["distance","task_id"]).iloc[0]
        records.append(dict(zip(group_keys,key), task_id=chosen.task_id, sample_af=chosen.sample_af,
                            selection_rule="closest to cohort median final AF; tie by task ID"))
        path=out/"regions"/chosen.task_id
        profile=pd.read_csv(path/"profiles.csv")
        local=profile[np.abs(profile.position_0based-cfg["focal_position_bp"])<=500000]
        fig,axes=plt.subplots(2,2,figsize=(14,11),layout="constrained")
        for r,source in enumerate(("truth","decoded")):
            for c,pair in enumerate(PAIR_LABELS):
                ax=axes[r,c]
                z=local[(local.source==source)&(local.pair_class==pair)]
                if int(z.n_pairs.iloc[0])==0:
                    ax.text(.5,.5,"No ALT/ALT pair available",ha="center",va="center",transform=ax.transAxes)
                else:
                    for cutoff in cfg["tmrca_cutoffs_years"]:
                        ax.plot((z.position_0based-cfg["focal_position_bp"])/1000,z[f"frac_recent_{cutoff}"],label=f"{cutoff//1000} kya",lw=1.8)
                ax.axvline(0,color="black",ls=":")
                ax.set(ylim=(-.03,1.03),xlabel="Distance from focal allele (kb)",ylabel="Raw frac_recent_T",title=f"{source} | {pair} | n={int(z.n_pairs.iloc[0])}")
        fig.suptitle(f"{chosen.task_id} | AF={chosen.sample_af:.3f}",fontsize=17)
        handles,labels=axes[0,0].get_legend_handles_labels()
        if handles:
            fig.legend(handles,labels,loc="outside upper center",ncol=5,frameon=False)
        save_figure(fig,out/"figures"/"examples",chosen.task_id.replace("/","_"))
    pd.DataFrame(records).to_csv(out/"examples.csv",index=False)


def analyze(root, cfg, plan):
    out=root/"analysis"
    out.mkdir(exist_ok=True)
    inventory=[]
    rows=[]
    source_receipts=[]
    for task in plan:
        directory=root/"regions"/task["id"]
        sim=json.loads((directory/"simulation.json").read_text())
        receipt=json.loads((directory/"decoded.json").read_text())
        oo.read_receipt(directory,"decoded.json",receipt["identity"])
        assert receipt['identity']['parameters']['tmrca_cutoffs_years']==cfg['tmrca_cutoffs_years']
        source_receipts.append(dict(task_id=task['id'],simulation=oo.legacy.digest(directory/'simulation.json'),decoded=oo.legacy.digest(directory/'decoded.json')))
        arm=("I50" if task["onset_years"]==50000 else "I10") if task["origin"]=="introgressed" and task["role"]=="selected" else task["family"]
        common=dict(task_id=task["id"],family=task["family"],arm=arm,role=task["role"],s=task["s"])
        inventory.append(dict(common,sample_af=sim["sample_af"],seed=sim["seed"],trajectory_available=sim["trajectory_available"]))
        frame=pd.read_csv(directory/"profiles.csv")
        focal=frame[(frame.position_0based==cfg["focal_position_bp"]) & frame.pair_class.isin(PAIR_LABELS)]
        if len(focal)!=4:
            raise ValueError(f"Missing or duplicated focal pair/source records: {directory}")
        for row in focal.itertuples():
            for cutoff in cfg["tmrca_cutoffs_years"]:
                rows.append(dict(common,source=row.source,method=PAIR_LABELS[row.pair_class],cutoff_years=cutoff,
                    score=getattr(row,f"frac_recent_{cutoff}"),n_pairs=row.n_pairs))
        rows.append(dict(common,source="genotypes",method="AF",cutoff_years=0,score=sim["sample_af"],n_pairs=np.nan))
    inventory=pd.DataFrame(inventory)
    if inventory.seed.duplicated().any():
        raise ValueError("Accepted simulation seeds are not disjoint")
    focal=pd.DataFrame(rows)
    predictions,metrics=evaluate(focal,cfg)
    inventory.to_csv(out/"inventory.csv",index=False)
    focal.to_csv(out/"focal_scores.csv",index=False)
    predictions.to_csv(out/"predictions.csv",index=False)
    metrics.to_csv(out/"metrics.csv",index=False)
    oo.legacy.atomic_json(out/"source_receipts.json",source_receipts)
    plot_results(out,cfg,focal,metrics)
    # Profiles live under root, while report products live under analysis.
    examples_root=out/"regions"
    if not examples_root.exists():
        examples_root.symlink_to(root/"regions",target_is_directory=True)
    example_profiles(out,cfg,inventory[inventory.role!="null"])
    # Cohort trajectories: accepted population trajectories only; no lost-origin inference.
    for arm in arms(cfg):
        group=inventory[(inventory.arm==arm)&(inventory.role=="selected")]
        fig,axes=plt.subplots(3,4,figsize=(18,13),layout="constrained")
        for ax,s in zip(axes.flat,cfg["selection_coefficients"]):
            available=0
            for row in group[group.s==s].itertuples():
                if row.trajectory_available:
                    trajectory=pd.read_csv(root/"regions"/row.task_id/"trajectory.csv")
                    ax.plot(trajectory.years_ago/1000,trajectory.population_af,lw=1,alpha=.55)
                    available+=1
            ax.set(title=f"s={s:g}; n={available}",ylim=(0,1.02),xlabel="kya",ylabel="Population AF")
            ax.invert_xaxis()
            if not available:
                ax.text(.5,.5,"Not saved in legacy trees",ha="center",va="center",transform=ax.transAxes,fontsize=12)
        fig.suptitle(f"{arm}: accepted-survivor trajectories (not establishment probabilities)")
        save_figure(fig,out/"figures",f"trajectories_{arm}")
    report=("# Origin/onset focal pilot\n\n"
        f"{cfg['target_replicates']} targets per cell; {cfg['null_replicates']} separate calibration nulls per family. "
        "All scores use the focal allele position and p<=0.05. Truth and decoded TMRCA are calibrated separately. "
        "ALT/ALT uses the raw within-class fraction, without carrier-mass weighting.\n\n"
        "Neutral detection rate is FPR among independent unselected targets, not mixed-prevalence discovery FDR. "
        "The calibration nulls are not test targets. "
        + ("I50 and I10 use distinct onset-ascertained controls. " if cfg.get('focal_ascertainment') else "I50 and I10 share introgressed controls. ") +
        "Missing ALT/ALT scores are NA and yield no discovery; tables include availability. "
        "Wilson intervals are conditional on the fitted null and omit null-estimation uncertainty. "
        f"The {cfg['target_replicates']}-target estimate changes in {100/cfg['target_replicates']:g}-percentage-point steps.\n\n"
        "Power is conditional on survival and sample observation; fixation is retained. "
        "All cutoffs are separate prespecified tests, not a best-cutoff or any-cutoff rule. "
        f"Population scaling factors by origin: {cfg['scaling_by_origin']}.\n")
    (out/"README.md").write_text(report)
    files=[p for p in out.rglob('*') if p.is_file() and 'regions' not in p.relative_to(out).parts]
    oo.legacy.atomic_json(out/"complete.json",dict(schema="origin-analysis/v1",config=cfg,
        task_count=len(plan),metric_cells=len(metrics),artifacts=oo.artifact_specs(out,[str(p.relative_to(out)) for p in files if p.name!='complete.json'])))
