"""Publication-style figures and a linked report for the stopped campaign."""
import os
os.environ['MPLBACKEND']='Agg'
from html import escape
import json
from pathlib import Path
import textwrap
import zipfile
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42,'font.size':17,
    'axes.titlesize':19,'axes.labelsize':17,'xtick.labelsize':13,'ytick.labelsize':14})
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
from . import origin_onset as oo
from .origin_onset_analysis import heatmap

LABELS={'all_pairs':'All sampled pairs','alt_alt':'ALT/ALT pairs'}

def make_report(root,out,cfg,inventory,focal,metrics):
    figures=out/'figures'; figures.mkdir(exist_ok=True)
    families=['I50','I10']; cutoffs=cfg['tmrca_cutoffs_years']; coeff=cfg['selection_coefficients']
    slabels=[f'{s:g}' for s in coeff]
    sources=[s for s in ['truth','decoded'] if s in set(focal.source)]
    cohort=pd.read_csv(out/'cohort_counts.csv',keep_default_na=False)
    null_counts=inventory[inventory.role=='null'].groupby('family').size().to_dict()
    neutral_counts=inventory[inventory.role=='neutral_target'].groupby('family').size().to_dict()
    calibration=pd.read_csv(out/'calibration_diagnostics.csv',keep_default_na=False)
    decoded_n=focal[focal.source=='decoded'].task_id.nunique()
    truth_n=focal[focal.source=='truth'].task_id.nunique()
    notes=[
        f'Stopped snapshot: {len(inventory):,} of 4,600 planned regions. No new simulations were generated.',
        f'Calibration nulls: I50 n={null_counts.get("I50",0)}, I10 n={null_counts.get("I10",0)}; planned n=1,000 each. Independent neutral targets: I50 n={neutral_counts.get("I50",0)}, I10 n={neutral_counts.get("I10",0)}.',
        f'True TMRCA scored for {truth_n:,} regions. Gamma-SMC scored for {decoded_n:,} regions. Each source is separately calibrated against its available nulls.',
        'All power and FPR estimates are provisional. Completion time depends on allele loss and retry behavior; completed runs may overrepresent faster-finishing trajectories. They are not a random subsample of the planned grid.',
        'Power is conditional on a focal archaic-derived allele segregating at onset and observed in the present-day sample. Lost or unobserved focal alleles were retried; fixation after onset was retained.',
        'Neutral detection among independent unselected introgressed targets estimates the false-positive rate (FPR), the neutral-only endpoint requested as FDR. Discovery FDR requires a mixture/prevalence and is not estimated here.',
        'At each focal site: p=(1 + number of calibration null scores >= target score)/(1 + actual null count); alpha=0.05. Ties are conservative. Each TMRCA cutoff is a separate test; there is no best-cutoff or any-cutoff claim.',
        'Raw scores are P(TMRCA < T) within the fixed sampled pair set, with no carrier-mass weighting. ALT/ALT uses carrier pairs among the same 10,000 sampled pairs, with membership fixed by the focal allele across the window. Zero available ALT/ALT pairs give an unavailable score and no call; coverage is tabulated.',
        'Truth uses strict TMRCA < T. The established Gamma-SMC mean-call implementation counts posterior mean TMRCA <= T using its native float32 threshold rule. Decoded scores are fractions of pairs, not posterior P(TMRCA < T) averages.',
        'Wilson 95% intervals condition on the fitted null distribution and omit null-estimation uncertainty. Only 97/98 nulls also mean coarser p-values and greater calibration uncertainty than the planned 1,000.',
        'Representative profiles are chosen closest to the cohort median present-day sample AF, with task-ID tie breaks. They illustrate individual replicates and are not averages. Profiles cover +/-500 kb around the focal site; decoding uses the full 10-Mb region.',
        'Demography: constant archaic Ne=3,600 diploids; split 500 kya; 2% EAS pulse at 50 kya; selection/focal ascertainment at 50 or 10 kya; generation time 25 years; Q=1. Focal choice is the nearest eligible archaic-derived EAS polymorphism in an 11-Mb region, cropped to 10 Mb.',
    ]
    artifacts=[]
    report_pdf=out/'comprehensive_report.pdf'
    pdf=PdfPages(report_pdf,metadata={'Title':'Introgressed allele selection: stopped campaign report',
        'Author':'Simulation analysis','Subject':'Provisional focal power and independent neutral FPR'})
    def save(fig,name,caption):
        if fig.legends:
            fig.set_layout_engine('constrained',rect=(0,0,1,.83))
            if fig._suptitle:
                fig._suptitle.set_y(.99); fig._suptitle.set_in_layout(False)
            for legend in fig.legends:
                legend.set_loc('upper center'); legend.set_bbox_to_anchor((.5,.91)); legend.set_in_layout(False)
        for ax in fig.axes:
            for coll in ax.collections:
                coll.set_rasterized(False)
        extra=list(fig.legends)+([fig._suptitle] if fig._suptitle else [])
        fig.savefig(figures/(name+'.png'),dpi=170,bbox_inches='tight',bbox_extra_artists=extra)
        fig.savefig(figures/(name+'.pdf'),bbox_inches='tight',bbox_extra_artists=extra)
        pdf.savefig(fig,bbox_inches='tight',bbox_extra_artists=extra)
        artifacts.append(dict(name=name,caption=caption))
        plt.close(fig)
    findings=[]
    saturated=calibration[(calibration.source=='truth')&(calibration.method=='alt_alt')&(calibration.cutoff_years>=20000)]
    if len(saturated)==8 and not saturated.any_bounded_score_can_reject.any():
        findings.append('For true TMRCA, ALT/ALT cannot reject at T=20–50 kya in this snapshot: '
            f'even a target score of 1 has p={saturated.p_at_max_possible_score.min():.3f}–{saturated.p_at_max_possible_score.max():.3f}. '
            'Many matched unselected-introgressed nulls already have a raw fraction of 1. This is score saturation with conservative ties, not evidence that selected genealogies lack recent coalescence.')
    for family,cutoff in [('I50',50000),('I10',10000)]:
        g=metrics[(metrics.family==family)&(metrics.source=='truth')&(metrics.method=='all_pairs')&
            (metrics.role=='selected')&(metrics.cutoff_years==cutoff)&(metrics.s.isin([.005,.01,.03]))].sort_values('s')
        pieces=[f's={r.s:g}: {r.calls}/{r.n} ({100*r.rate:.1f}%)' for r in g.itertuples()]
        findings.append(f'{family}, true all-pairs detection at T={cutoff//1000} kya: '+ '; '.join(pieces)+'. These fixed-cutoff examples are not an optimized-cutoff test.')
    truth_fpr=metrics[(metrics.source=='truth')&(metrics.role=='neutral_target')]
    if len(truth_fpr):
        findings.append(f'True-TMRCA neutral FPR ranges from {100*truth_fpr.rate.min():.1f}% to {100*truth_fpr.rate.max():.2f}% across the tested pair classes and cutoffs. '
            'The 96/98 neutral test targets and 97/98 calibration nulls provide limited precision; uncertainty from calibration is not included in the Wilson intervals.')
    fig=plt.figure(figsize=(11,8.5))
    fig.text(.06,.94,'Findings in the completed snapshot',fontsize=23,weight='bold')
    fig.text(.06,.88,'Provisional estimates; full denominators and intervals in results_table.csv',fontsize=14)
    y=.80
    for paragraph in findings:
        wrapped=textwrap.wrap(paragraph,94)
        fig.text(.06,y,'\n'.join(wrapped),va='top',fontsize=13.5,linespacing=1.3)
        y-=.032*len(wrapped)+.035
    save(fig,'findings','Descriptive findings from the completed snapshot; not a replacement for the planned calibration cohort.')
    for page in range(3):
        fig=plt.figure(figsize=(11,8.5))
        fig.text(.06,.94,'Introgressed allele selection',fontsize=25,weight='bold')
        fig.text(.06,.88,f'Stopped campaign | provisional report | methods {page+1}/3',fontsize=18)
        y=.80
        for paragraph in notes[page*4:(page+1)*4]:
            lines=textwrap.wrap(paragraph,90)
            fig.text(.06,y,'\n'.join(lines),va='top',fontsize=14,linespacing=1.35)
            y-=.035*len(lines)+.035
        fig.text(.06,.04,'Detailed denominators, scores, p-values, intervals and provenance accompany this PDF.',fontsize=12)
        save(fig,f'methods_{page+1}','Scope, definitions and limits of this stopped snapshot.')
    # Denominators for every selected treatment are visible alongside the controls.
    fig,axes=plt.subplots(1,2,figsize=(14,8.5),layout='constrained')
    for ax,family in zip(axes,families):
        z=cohort[cohort.family==family]
        selected=z[z.role=='selected'].set_index('s').reindex(coeff)
        ax.barh(np.arange(len(coeff)),selected.completed,color='#4477AA')
        ax.set_yticks(np.arange(len(coeff)),slabels); ax.invert_yaxis()
        ax.set(xlim=(0,110),xlabel='Completed of 100 planned targets',ylabel='Selection coefficient s',
            title=f'{family}: null n={null_counts[family]}, neutral n={neutral_counts[family]}')
        for i,n in enumerate(selected.completed): ax.text(n+1,i,str(int(n)),va='center',fontsize=13)
        ax.axvline(100,color='grey',ls='--')
    fig.suptitle('Completed cohort sizes | calibration nulls planned: 1,000 per onset',fontsize=20)
    save(fig,'cohort_completion','Actual denominators; all missing tasks are listed in unfinished_tasks.csv.')
    # A fixed FPR scale is shared by truth, decoded and both pair classes.
    fpr_max=100
    for source in sources:
        for method in LABELS:
            z=metrics[(metrics.source==source)&(metrics.method==method)]
            if z.empty: continue
            fig,axes=plt.subplots(1,2,figsize=(14,8.5),layout='constrained')
            for ax,family in zip(axes,families):
                selected=z[(z.family==family)&(z.role=='selected')]
                values=selected.pivot(index='s',columns='cutoff_years',values='rate').reindex(index=coeff,columns=cutoffs).to_numpy()*100
                mesh=heatmap(ax,values,[t//1000 for t in cutoffs],slabels,100,f'{family} | null n={null_counts[family]}')
                ax.set_ylabel('Selection coefficient s')
            bar=fig.colorbar(mesh,ax=axes,fraction=.025,pad=.025,label='Focal power (%)'); bar.solids.set_rasterized(False)
            fig.suptitle(f'{source.capitalize()} | {LABELS[method]} | p <= 0.05\nProvisional; actual target counts in results table',fontsize=20)
            save(fig,f'power_{source}_{method}','Conditional focal power; no-call targets remain in the denominator.')
            neutral=z[z.role=='neutral_target']
            values=neutral.pivot(index='family',columns='cutoff_years',values='rate').reindex(index=families,columns=cutoffs).to_numpy()*100
            fig,ax=plt.subplots(figsize=(11,6.5),layout='constrained')
            mesh=heatmap(ax,values,[t//1000 for t in cutoffs],families,fpr_max,'Independent unselected introgressed targets')
            bar=fig.colorbar(mesh,ax=ax,fraction=.025,pad=.025,label='False-positive rate (%)'); bar.solids.set_rasterized(False)
            fig.suptitle(f'{source.capitalize()} | {LABELS[method]} | provisional FPR\nCommon color scale: 0–100%',fontsize=20)
            save(fig,f'fpr_{source}_{method}','The same 0–100% scale is used for every FPR heatmap.')
    # Intervals make the limited neutral target precision visible without changing heatmap scales.
    fig,axes=plt.subplots(1,2,figsize=(14,7.5),layout='constrained')
    for ax,family in zip(axes,families):
        for j,(source,method) in enumerate((s,m) for s in sources for m in LABELS):
            g=metrics[(metrics.family==family)&(metrics.role=='neutral_target')&(metrics.source==source)&(metrics.method==method)].sort_values('cutoff_years')
            if g.empty: continue
            ax.errorbar(g.cutoff_years/1000+(j-1.5)*.6,g.rate*100,
                yerr=np.maximum(0,np.array([g.rate-g.ci_low,g.ci_high-g.rate])*100),marker='o',capsize=3,
                label=f'{source} / {LABELS[method]}')
        ax.axhline(5,color='grey',ls='--')
        ax.set(ylim=(0,100),xlabel='TMRCA cutoff (kya)',ylabel='FPR (%) with 95% Wilson interval',title=family)
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,ncol=2,frameon=False,fontsize=13)
    fig.suptitle('Independent neutral FPR | fitted-null uncertainty not included',fontsize=20)
    save(fig,'fpr_intervals','Neutral FPR with conditional Wilson intervals and the nominal 5% reference.')
    fig,axes=plt.subplots(1,2,figsize=(14,8.5),layout='constrained')
    for ax,family in zip(axes,families):
        g=inventory[inventory.family==family]
        data=[g[g.role=='neutral_target'].sample_af.to_numpy()]+[g[(g.role=='selected')&(g.s==s)].sample_af.to_numpy() for s in coeff]
        ax.boxplot(data,tick_labels=['0']+slabels,showfliers=True)
        null=g[g.role=='null'].sample_af
        median=null.median()
        ax.axhline(median,color='grey',ls='--',label='Calibration-null median')
        ax.set(ylim=(-.02,1.03),title=family,xlabel='s (0 = independent neutral targets)',ylabel='Present-day sample ALT frequency')
        ax.tick_params(axis='x',rotation=60)
    h,l=axes[0].get_legend_handles_labels(); fig.legend(h,l,ncol=2,frameon=False,fontsize=14)
    fig.suptitle('Allele-frequency distributions | retained focal alleles',fontsize=22)
    save(fig,'allele_frequency','Final sample AF, with matched calibration null median.')
    # Frequency-only benchmark, calibrated exactly like the pair statistics.
    fig,axes=plt.subplots(1,2,figsize=(14,7.5),layout='constrained')
    for ax,family in zip(axes,families):
        g=metrics[(metrics.family==family)&(metrics.method=='AF')&(metrics.role=='selected')].sort_values('s')
        ax.errorbar(np.arange(len(g)),g.rate*100,yerr=np.maximum(0,np.array([g.rate-g.ci_low,g.ci_high-g.rate])*100),marker='o',capsize=4)
        ax.set_xticks(np.arange(len(g)),[f'{s:g}' for s in g.s],rotation=60)
        neutral=metrics[(metrics.family==family)&(metrics.method=='AF')&(metrics.role=='neutral_target')]
        ax.set(ylim=(-2,102),xlabel='Selection coefficient s',ylabel='AF-only detection (%)',title=f'{family}; neutral FPR={100*neutral.rate.iloc[0]:.1f}%')
    fig.suptitle('Allele-frequency-only benchmark | p <= 0.05',fontsize=22)
    save(fig,'af_only_power','AF-only focal power with Wilson intervals; independent neutral FPR in panel titles.')
    for family in families:
        for cutoff in cutoffs:
            fig,axes=plt.subplots(len(sources),2,figsize=(14,5.5*len(sources)),layout='constrained',squeeze=False)
            for i,source in enumerate(sources):
                for j,method in enumerate(LABELS):
                    ax=axes[i,j]
                    g=focal[(focal.family==family)&(focal.source==source)&(focal.method==method)&(focal.cutoff_years==cutoff)]
                    data=[g[g.role=='neutral_target'].score.dropna().to_numpy()]+[g[(g.role=='selected')&(g.s==s)].score.dropna().to_numpy() for s in coeff]
                    ax.boxplot(data,tick_labels=['0']+slabels)
                    null=g[g.role=='null'].score.dropna()
                    if len(null):
                        ax.axhline(null.median(),color='grey',ls='--',label='Null median')
                    ax.set(ylim=(-.03,1.03),title=f'{source} | {LABELS[method]}',xlabel='s (0 = neutral targets)',ylabel=f'Raw P(TMRCA < {cutoff//1000} kya)')
                    ax.tick_params(axis='x',rotation=60)
            h,l=axes[0,0].get_legend_handles_labels(); fig.legend(h,l,ncol=2,frameon=False,fontsize=14)
            fig.suptitle(f'{family} | raw focal fractions | T = {cutoff//1000} kya',fontsize=22)
            save(fig,f'raw_{family}_T{cutoff}','Unweighted within-class fractions; unavailable ALT/ALT classes are omitted from boxplots and tabulated as no-calls.')
    for source in sources:
        fig,axes=plt.subplots(1,2,figsize=(14,7.5),layout='constrained')
        for ax,family in zip(axes,families):
            g=focal[(focal.family==family)&(focal.source==source)&(focal.method=='alt_alt')&(focal.cutoff_years==cutoffs[0])]
            groups=[g[g.role=='neutral_target']]+[g[(g.role=='selected')&(g.s==s)] for s in coeff]
            values=[100*d.score.notna().mean() if len(d) else np.nan for d in groups]
            ax.bar(np.arange(len(values)),values,color='#228833'); ax.set_xticks(np.arange(len(values)),['0']+slabels,rotation=60)
            ax.set(ylim=(0,105),title=family,xlabel='s (0 = neutral targets)',ylabel='ALT/ALT score availability (%)')
        fig.suptitle(f'{source.capitalize()} | at least one sampled ALT/ALT pair',fontsize=22)
        save(fig,f'coverage_{source}','Pair availability is separate from power; no-call targets count as nondetections.')
    for source in sources:
        fig,axes=plt.subplots(1,2,figsize=(14,7.5),layout='constrained')
        for ax,method in zip(axes,LABELS):
            for family in families:
                g=calibration[(calibration.family==family)&(calibration.source==source)&(calibration.method==method)].sort_values('cutoff_years')
                ax.plot(g.cutoff_years/1000,g.score_one_fraction*100,marker='o',lw=2,label=family)
            ax.set(ylim=(-2,102),xlabel='TMRCA cutoff (kya)',ylabel='Calibration nulls with raw score = 1 (%)',title=LABELS[method])
        h,l=axes[0].get_legend_handles_labels(); fig.legend(h,l,ncol=2,frameon=False)
        fig.suptitle(f'{source.capitalize()}: null-score saturation at the upper bound',fontsize=22)
        save(fig,f'null_saturation_{source}','When too many nulls equal 1, even a target score of 1 cannot achieve p <= 0.05. See calibration_diagnostics.csv.')
    # Representative local raw profiles for every selected and neutral treatment.
    selected_examples=pd.read_csv(out/'examples.csv',keep_default_na=False)
    for r in selected_examples.itertuples():
        fig,axes=plt.subplots(len(sources),2,figsize=(14,5*len(sources)),layout='constrained',squeeze=False)
        for i,source in enumerate(sources):
            path=out/'regions'/r.task_id/(source+'.csv')
            profile=pd.read_csv(path) if path.exists() else pd.DataFrame()
            for j,method in enumerate(LABELS):
                ax=axes[i,j]
                z=profile[profile.method==method] if len(profile) else profile
                if z.empty or not z.score.notna().any():
                    ax.text(.5,.5,'Unavailable pair class / source',ha='center',va='center',transform=ax.transAxes,fontsize=14)
                else:
                    for cutoff in cutoffs:
                        local=z[z.cutoff_years==cutoff].sort_values('position_0based')
                        ax.plot((local.position_0based-cfg['focal_position_bp'])/1000,local.score,label=f'{cutoff//1000} kya',lw=1.7)
                n=int(z.n_pairs.iloc[0]) if len(z) else 0
                ax.axvline(0,color='black',ls=':')
                ax.set(ylim=(-.03,1.03),xlabel='Distance from focal allele (kb)',ylabel='Raw frac_recent_T',title=f'{source} | {LABELS[method]} | pairs={n:,}')
        h,l=axes[0,0].get_legend_handles_labels()
        if h: fig.legend(h,l,ncol=5,frameon=False,fontsize=14)
        fig.suptitle(f'{r.task_id}\nSample AF={r.sample_af:.3f}; representative nearest median AF',fontsize=18)
        save(fig,'example_'+r.task_id.replace('/','_'),'One prespecified median-AF representative for this cohort; raw spatial profile, not a treatment mean.')
    # Accepted trajectories: individual lines plus medians and IQR, one family per page.
    for family in families:
        fig,axes=plt.subplots(3,4,figsize=(18,13),layout='constrained')
        for ax,s in zip(axes.flat,coeff):
            g=inventory[(inventory.family==family)&(inventory.role=='selected')&(inventory.s==s)]
            traces=[]; years=None
            for r in g.itertuples():
                trajectory=pd.read_csv(root/'regions'/r.task_id/'trajectory.csv')
                years=trajectory.years_ago.to_numpy()/1000
                traces.append(trajectory.population_af.to_numpy())
            if traces:
                a=np.asarray(traces); q1,median,q3=np.quantile(a,[.25,.5,.75],axis=0)
                for values in a: ax.plot(years,values,color='#4477AA',alpha=.10,lw=.7)
                ax.fill_between(years,q1,q3,color='#4477AA',alpha=.25)
                ax.plot(years,median,color='black',lw=2)
                ax.invert_xaxis()
            ax.set(ylim=(0,1.02),title=f's={s:g}; n={len(g)}',xlabel='Time (kya)',ylabel='Population AF')
        fig.suptitle(f'{family}: accepted-survivor trajectories | black: median; band: IQR',fontsize=22)
        save(fig,'trajectories_'+family,'Only accepted and completed trajectories; these are not establishment probabilities.')
    fig,axes=plt.subplots(1,2,figsize=(14,7.5),layout='constrained')
    for ax,family in zip(axes,families):
        g=inventory[inventory.family==family]
        data=[g[g.role=='null'].attempts.to_numpy(),g[g.role=='neutral_target'].attempts.to_numpy()]+[g[(g.role=='selected')&(g.s==s)].attempts.to_numpy() for s in coeff]
        ax.boxplot(data,tick_labels=['Null','0']+slabels); ax.set_yscale('log'); ax.tick_params(axis='x',rotation=60)
        ax.set(title=family,xlabel='Cohort / selection coefficient s',ylabel='Attempts per retained region (log scale)')
    fig.suptitle('Retry burden among completed regions',fontsize=22)
    save(fig,'retry_burden','Only completed regions; interrupted attempts and unfinished tasks are listed separately.')
    pdf.close()
    oo.legacy.atomic_json(out/'plot_scales.json',dict(power_percent=[0,100],neutral_fpr_percent=[0,100],raw_fraction=[0,1],af=[0,1]))
    tables=['cohort_counts.csv','results_table.csv','allele_frequency_and_attempts.csv','raw_fraction_summary.csv','calibration_diagnostics.csv',
        'inventory.csv','focal_scores.csv','predictions.csv','examples.csv','unfinished_tasks.csv']
    # A concise comparison table joins truth and decoded on exactly the same treatment endpoints.
    compare=metrics[metrics.source.isin(['truth','decoded'])].pivot(index=['family','role','s','method','cutoff_years'],columns='source',values='rate').reset_index()
    if all(s in compare for s in ['truth','decoded']): compare['decoded_minus_truth_pp']=100*(compare.decoded-compare.truth)
    compare.to_csv(out/'truth_decoded_comparison.csv',index=False); tables.append('truth_decoded_comparison.csv')
    lines=['# Introgressed allele selection: stopped campaign report','','## Findings','',*sum(([n,''] for n in findings),[]),
        '## Scope and interpretation','',*sum(([n,''] for n in notes),[]),
        '## Files','', '[Full figure and methods PDF](comprehensive_report.pdf)', '[Browsable report](index.html)','']
    lines += [f'- [{t}]({t})' for t in tables]
    lines += ['','## Provenance','', 'snapshot.json freezes the included task IDs and original configuration. Per-region truth/decoded receipts hash the source simulation receipt and pair manifest. Tree, carrier-mask and trajectory inputs are hash-checked; focal genotypes are independently compared with the carrier mask. This is a report-input audit, not a new full ascertainment re-audit.','',
        'The original run remains stopped. Report processing never invokes SLiM or generates replacement nulls.']
    (out/'README.md').write_text('\n'.join(lines)+'\n')
    parts=['<!doctype html><html><head><meta charset="utf-8"><title>Stopped simulation report</title>',
        '<style>body{max-width:1200px;margin:36px auto;padding:0 24px;font:18px/1.5 system-ui;color:#18212d}h1,h2{line-height:1.2}img{width:100%;height:auto}section{margin:3em 0}table{border-collapse:collapse;font-size:15px}td,th{padding:7px 12px;border:1px solid #ddd}a{color:#145ca1}.notice{padding:20px;background:#fff4d7;border-left:5px solid #d69625}</style></head><body>',
        '<h1>Introgressed allele selection</h1><p>Stopped campaign · provisional report · 28 September 2026</p>',
        f'<div class="notice"><strong>{len(inventory):,}/4,600 completed regions.</strong> Truth: {truth_n:,}; Gamma-SMC: {decoded_n:,}. Calibration nulls: I50 {null_counts["I50"]}, I10 {null_counts["I10"]} of 1,000 planned each.</div>',
        '<p><a href="comprehensive_report.pdf">Download comprehensive PDF</a> · <a href="report_bundle.zip">Download figures and tables ZIP</a></p>',
        '<h2>Findings</h2>',*['<p>'+escape(n)+'</p>' for n in findings],
        '<h2>Interpretation and limitations</h2>',*['<p>'+escape(n)+'</p>' for n in notes],
        '<h2>Results tables</h2><ul>',*['<li><a href="'+t+'">'+t+'</a></li>' for t in tables],'</ul>',
        cohort.to_html(index=False),'<h2>Figures</h2>']
    for item in artifacts:
        name=item['name']; caption=item['caption']
        parts += [f'<section><h2>{escape(name.replace("_"," "))}</h2><p>{escape(caption)}</p>',
            f'<a href="figures/{name}.pdf">Vector PDF</a> · <a href="figures/{name}.png">PNG</a>',
            f'<img loading="lazy" src="figures/{name}.png" alt="{escape(caption)}"></section>']
    parts += ['</body></html>']; (out/'index.html').write_text('\n'.join(parts))
    deliverables=[out/t for t in tables]+[out/'README.md',out/'index.html',report_pdf,out/'snapshot.json',out/'plot_scales.json']
    deliverables+=list(figures.glob('*.png'))+list(figures.glob('*.pdf'))
    manifest=dict(schema='stopped-report-deliverables/v1',regions=len(inventory),truth_n=truth_n,decoded_n=decoded_n,
        metric_rows=len(metrics),figure_count=len(artifacts),provisional=True,
        report_code_sha256=oo.legacy.digest(Path(__file__)),
        analysis_code_sha256=oo.legacy.digest(Path(__file__).with_name('stopped_campaign_report.py')),
        artifacts=oo.artifact_specs(out,[str(p.relative_to(out)) for p in deliverables]))
    oo.legacy.atomic_json(out/'report_complete.json',manifest)
    with zipfile.ZipFile(out/'report_bundle.zip','w',compression=zipfile.ZIP_DEFLATED) as archive:
        for path in deliverables+[out/'report_complete.json']:
            archive.write(path,path.relative_to(out))
    print(json.dumps({k:v for k,v in manifest.items() if k!='artifacts'}),flush=True)
