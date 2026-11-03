#!/usr/bin/env python3
"""Report only strictly segregating proposals with independently conditioned nulls."""
import os
os.environ['MPLBACKEND']='Agg'
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from report_ongoing_selection_pilot import sha, rate, evaluate, METHODS
from report_time_resolved_pilot import empirical_p, interval_values, expand_support


def is_segregating(f):
    return 0 < f < 1


def report(root,out,neutral_reference=None):
    out.mkdir(parents=True,exist_ok=True)
    manifest=json.loads((root/'manifest.json').read_text());signature=sha(root/'manifest.json')
    tasks=pd.DataFrame(json.loads((root/'tasks.json').read_text()))
    proposals=[];arrays={}
    for task in tasks.to_dict('records'):
        directory=root/task['unit']
        path=directory/'complete.json'
        if not path.exists():path=directory/'rejected.json'
        if not path.exists():continue
        r=json.loads(path.read_text());assert r['config_sha256']==signature
        keep=path.name=='complete.json';assert keep==is_segregating(r['final_frequency'])
        if 'source_receipt' in r:
            # Linux path recorded by the runner; report can also run on Windows.
            source=r['source_receipt']
            if os.name=='nt' and source.startswith('/mnt/d/'):source='D:/'+source[7:]
            assert sha(source)==r['source_receipt_sha256']
        if keep:
            for key,name in [('scores_sha256','scores.npz'),('observed_sha256','observed.trees'),('posterior_sha256','posterior.zst'),('onset_sha256','onset.trees')]:assert sha(directory/name)==r[key]
            trajectory=np.loadtxt(directory/'trajectory.tsv',skiprows=1)
            np.testing.assert_array_equal(trajectory[:,0],np.arange(50000,-1,-25))
            assert is_segregating(trajectory[-1,1])
            with np.load(directory/'scores.npz') as a:
                np.testing.assert_array_equal(a['positions'],np.arange(0,1000000,10000))
                np.testing.assert_array_equal(a['cutoffs'],np.arange(1000,50001,1000))
                np.testing.assert_array_equal(a['valid_pairs'],np.full(100,1000))
                arrays[r['unit']]=a['cdfs'].copy()
        proposals.append(dict(task,accepted=keep,final_frequency=r['final_frequency'],
            rejected_loss=not keep and r['final_frequency']==0,rejected_fixation=not keep and r['final_frequency']==1,
            seed=r['seed'],receipt_sha256=sha(path)))
    records=pd.DataFrame(proposals);assert records.seed.is_unique
    # Deadline censoring favors fast rejections. Use only fully completed
    # deterministic proposal prefixes, never a completion-speed-selected set.
    usable=set();counts=[]
    for key,g in tasks.groupby(['origin','cohort','selection']):
        g=g.sort_values('replicate');done=records[records.unit.isin(g.unit)]
        known=set(done.unit);prefix=[]
        for unit in g.unit:
            if unit not in known:break
            prefix.append(unit)
        usable.update(prefix);d=done[done.unit.isin(prefix)]
        counts.append(dict(origin=key[0],cohort=key[1],selection=key[2],scheduled=len(g),
            completed=len(done),completed_prefix=len(d),accepted=int(d.accepted.sum()),
            rejected_loss=int(d.rejected_loss.sum()),rejected_fixation=int(d.rejected_fixation.sum()),
            accepted_outside_prefix=int(done.accepted.sum()-d.accepted.sum())))
    counts=pd.DataFrame(counts);counts.to_csv(out/'proposal_summary.tsv',sep='\t',index=False)
    records['included_prefix']=records.unit.isin(usable)
    records.to_csv(out/'proposals.tsv',sep='\t',index=False)
    accepted=records[records.accepted&records.included_prefix].reset_index(drop=True)
    accepted_proposal_count=len(accepted)
    neutral_focals=[]
    if neutral_reference is not None:
        from report_ongoing_selection_pilot import read
        import tskit
        reference_manifest,reference,reference_cdfs,_=read(neutral_reference)
        for key in ('length','diploids','pairs','generation_time','mutation_rate','recombination_rate','stride',
                    'phlash_sha256','pulse_fraction','archaic_split_years','archaic_ne','decoder_sha256'):
            assert reference_manifest['config'][key]==manifest['config'][key],f'Unmatched neutral parameter: {key}'
        # Prespecified complete blocks from the two original scheduling lanes:
        # first ten main calibration units and first ten reverse-ordered tail
        # units (98..89); first ten development and held-out test units.
        first_twenty=all(set(range(20)).issubset(set(g.replicate)) for _,g in reference[reference.cohort=='calibration'].groupby('origin'))
        calibration_mask=(reference.replicate<20) if first_twenty else ((reference.replicate<10)|(reference.replicate>=89))
        mask=(reference.cohort!='selected') & np.where(reference.cohort=='calibration',calibration_mask,reference.replicate<10)
        null=reference[mask].copy()
        assert all(len(g)==(20 if cohort=='calibration' else 10) for (_,cohort),g in null.groupby(['origin','cohort']))
        for i,row in null.iterrows():
            arrays[row.unit]=reference_cdfs[i]
            ts=tskit.load(neutral_reference/row.unit/'observed.trees')
            center=ts.sequence_length/2
            found=[]
            for variant in ts.variants():
                genotypes=variant.genotypes
                alleles,counts_=np.unique(genotypes[genotypes>=0],return_counts=True)
                if len(alleles)>1:
                    found.append((abs(variant.site.position-center),variant.site.position,int(counts_.min()),len(genotypes)))
            if not found:raise ValueError(f'No segregating neutral focal variant: {row.unit}')
            distance,position,minor,total=min(found)
            neutral_focals.append(dict(unit=row.unit,origin=row.origin,cohort=row.cohort,position_bp=position,
                offset_from_center_bp=position-center,minor_allele_count=minor,sampled_haplotypes=total,
                definition='Nearest variant polymorphic in the present-day sample; ties choose lower coordinate'))
        accepted=pd.concat([accepted[accepted.cohort=='selected'],null],ignore_index=True)
        pd.DataFrame(neutral_focals).to_csv(out/'neutral_focal_variants.tsv',sep='\t',index=False)
    metrics=[];rows=[];unavailable=[]
    positions=np.arange(0,1000000,10000)
    for origin in ('de_novo','introgression'):
        r=accepted[accepted.origin==origin].reset_index(drop=True)
        development=(r.cohort=='development').to_numpy();cal=(r.cohort=='calibration').to_numpy()
        neutral=(r.cohort=='neutral_test').to_numpy()
        if development.sum()<5 or cal.sum()<19 or neutral.sum()<1:
            unavailable.append(dict(origin=origin,development=int(development.sum()),calibration=int(cal.sum()),neutral_test=int(neutral.sum()),
                reason='Need at least 5 independent development nulls, 19 calibration nulls for p<=0.05, and one held-out neutral region.'))
            continue
        cdfs=np.stack([arrays[u] for u in r.unit])
        for scan,widths in [('one_kyr',[1]),('multiscale',[1,2,5,10]),('fixed_50kya',None)]:
            values,intervals=interval_values(cdfs,widths) if widths else (cdfs[...,-1:],None)
            for j,method in enumerate(METHODS):
                mu=values[development,j].mean(axis=(0,1));sd=values[development,j].std(axis=(0,1),ddof=1)
                valid=sd>0;z=np.full(values[:,j].shape,-np.inf);z[...,valid]=(values[:,j][...,valid]-mu[valid])/sd[valid]
                maxima=z[cal].max(axis=(1,2));tested=[]
                for i in np.flatnonzero(neutral|(r.cohort=='selected').to_numpy()):
                    calls=empirical_p(z[i],maxima)<=.05
                    support=expand_support(calls,intervals) if widths else np.repeat(calls,49,axis=1)
                    m=evaluate(support,positions,r.iloc[i].cohort=='selected')
                    if widths is None:m['temporal_f1']=None;m['called_kyr']=None
                    rec=dict(origin=origin,scan=scan,method=method,unit=r.iloc[i].unit,
                        cohort=r.iloc[i].cohort,selection=r.iloc[i].selection,**m)
                    tested.append(rec);rows.append(rec)
                frame=pd.DataFrame(tested);null=frame[frame.cohort=='neutral_test']
                for s in sorted(tasks[tasks.cohort=='selected'].selection.unique()):
                    sel=frame[(frame.cohort=='selected')&(frame.selection==s)]
                    tp=int(sel.spatial_hit.sum());fp=int(null.any_call.sum());fn=len(sel)-tp
                    f1=2*tp/(2*tp+fp+fn) if len(sel) and 2*tp+fp+fn else None
                    for metric,data in [('spatial_tpr',sel.spatial_hit),('neutral_fpr',null.any_call)]:
                        metrics.append(dict(origin=origin,scan=scan,method=method,selection=s,metric=metric,
                            region_f1=f1,mean_temporal_f1=float(sel.temporal_f1.mean()) if len(sel) and widths else None,**rate(data)))
    pd.DataFrame(rows).to_csv(out/'replicate_metrics.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(metrics).to_csv(out/'performance.tsv',sep='\t',index=False,na_rep='NA')
    result=dict(manifest=manifest,unavailable=unavailable,completed_proposals=len(records),
        unavailable_selected_conditions=[dict(origin=origin,selection=float(s),reason='No accepted segregating selected replicate')
            for origin in ('de_novo','introgression') for s in sorted(tasks[tasks.cohort=='selected'].selection.unique())
            if not np.any((accepted.origin==origin)&(accepted.cohort=='selected')&(accepted.selection==s))],
        usable_prefix_proposals=int(records.included_prefix.sum()),accepted_prefix_proposals=accepted_proposal_count,
        selected_accepted_prefix_proposals=int((accepted.cohort=='selected').sum()),analysis_regions=len(accepted),
        rejected_loss=int(records.rejected_loss.sum()),rejected_fixation=int(records.rejected_fixation.sum()),
        scheduled=len(tasks),report_sha256=sha(__file__),root=str(root),neutral_reference=str(neutral_reference) if neutral_reference else None,
        null_policy='Unconditioned region; focal label moved to nearest present-day segregating variant' if neutral_reference else 'Focal survival conditioned',
        integrity='All accepted score, tree and posterior hashes checked; exact trajectory grid and strict segregation gate verified.')
    (out/'provenance.json').write_text(json.dumps(result,indent=2)+'\n')
    lines=['# Present-day-segregating EAS selection pilot','',
        'The primary condition is strictly 0 < whole-population focal allele frequency at present < 1. Losses and fixations are both rejected. The same condition applies independently to development, calibration, neutral-test and selected cohorts.',
        'Selection starts at 50 kya and remains enabled to present. The corrected coefficient grid is 0.005, 0.01, 0.02, 0.03. Stronger recent-onset comparisons (around 5 kya) are reserved for a separate experiment; no recent-onset result is claimed here.',
        'The introgression pulse is 1% at 50 kya, with 500 kya human–Neanderthal divergence, Neanderthal Ne=3,600, and an allele fixed in donors. De novo starts from one copy. EAS PHLASH median demography; 25 years/generation; h=0.5; Q=1; 1 Mb; 100 diploids; 1,000 unfiltered pairs; 10 kb stride; cutoffs 1–50 kya.',
        'Rejected trajectories stop at absorbing loss/fixation. This preserves exact rejection conditioning because no recurrent focal mutation is introduced. Previously completed proposals are reused with receipt and artifact hashes. Cohorts keep their originally assigned seeds; accepted simulations are never moved between cohorts.',
        'Analysis uses only the fully completed proposal prefix within each origin/cohort/coefficient stratum; later completions are recorded but excluded. This avoids choosing individual replicates by completion order, but deadline-dependent stopping and small cohorts remain limitations. Wilson intervals are exploratory.',
        'Calibration uses a maximum across the genomic/time grid, independent null development, and conservative +1 upper-tail p-values. At least 19 independent calibration nulls are needed to resolve p<=0.05. Missing selected cohorts are not reported as zero power.',
        'All scanned times lie inside the selection interval, so temporal overlap does not demonstrate dating accuracy. Temporal F1 measures [1,50) kya coverage. Whole-population segregation does not require both alleles to appear in the sampled 100 diploids.',
        '',f"Completed proposals: {len(records)}; usable prefix proposals: {int(records.included_prefix.sum())}; accepted selected proposals: {int((accepted.cohort=='selected').sum())}.",'',
        '| Origin | Cohort | s | Completed prefix | Accepted | Lost | Fixed |','|---|---|---:|---:|---:|---:|---:|']
    displayed=counts[counts.cohort=='selected'] if neutral_reference is not None else counts
    for c in displayed.to_dict('records'):lines.append(f"| {c['origin']} | {c['cohort']} | {c['selection']:g} | {c['completed_prefix']} | {c['accepted']} | {c['rejected_loss']} | {c['rejected_fixation']} |")
    lines.extend(['','## Performance availability','']+[f"- {u['origin']}: development={u['development']}, calibration={u['calibration']}, held-out neutral={u['neutral_test']}. {u['reason']}" for u in unavailable])
    if metrics:lines.append('Conditional TPR, neutral FPR, region F1 and Wilson 95% intervals are in performance.tsv.')
    if neutral_reference is not None:
        lines[2]='Selected simulations require strictly 0 < focal allele frequency at present < 1. Neutral genealogies are unconditioned; their focal label is the nearest present-day segregating variant. This supersedes the exploratory rejection of neutral focal alleles.'
        lines.insert(3,'Calibration uses 10 independent development, 20 calibration and 10 held-out neutral regions per origin, from complete blocks of the original main and supplemental schedules. The 10 kb grid is unchanged; whole-region maxima and FPR are unaffected by the neutral focal label. The variant is polymorphic in the sampled 200 haplotypes, which establishes present-day segregation.')
        lines.append('The exploratory neutral-rejection outcomes remain in proposal_summary.tsv for provenance; they are not used for calibration or FPR in this final analysis.')
        lines.extend(['','| s | Posterior-mass 1-kyr TPR | Neutral FPR | Region F1 | Temporal coverage F1 |','|---:|---:|---:|---:|---:|'])
        for m in metrics:
            if m['origin']=='introgression' and m['scan']=='one_kyr' and m['method']=='mass' and m['metric']=='spatial_tpr':
                fmt=lambda value:'NA' if value is None else f'{value:.3f}'
                fpr=next(x for x in metrics if x['origin']==m['origin'] and x['scan']==m['scan'] and x['method']==m['method'] and x['selection']==m['selection'] and x['metric']=='neutral_fpr')
                tpr=f"{m['k']}/{m['n']}" if m['n'] else 'Not estimable (0 accepted)'
                lines.append(f"| {m['selection']:g} | {tpr} | {fpr['k']}/{fpr['n']} | {fmt(m['region_f1'])} | {fmt(m['mean_temporal_f1'])} |")
    (out/'README.md').write_text('\n'.join(lines)+'\n')
    plot(out,counts)


def plot(out,counts):
    import matplotlib
    matplotlib.rcParams.update({'font.size':16,'axes.titlesize':18,'pdf.fonttype':42,'ps.fonttype':42})
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,6),sharey=True)
    for ax,origin in zip(axes,('de_novo','introgression')):
        g=counts[(counts.origin==origin)&(counts.cohort=='selected')].sort_values('selection')
        x=np.arange(len(g));values=np.divide(g.accepted.to_numpy(),g.completed_prefix.to_numpy(),out=np.zeros(len(g)),where=g.completed_prefix.to_numpy()>0)
        ax.bar(x,values,color='#0072b2',width=.65)
        for i,row in enumerate(g.itertuples()):ax.text(i,values[i]+.025,f'{row.accepted}/{row.completed_prefix}',ha='center',fontsize=14)
        ax.set(xticks=x,xticklabels=[f'{s:g}' for s in g.selection],xlabel='Selection coefficient s',title=origin.replace('_',' '),ylim=(0,1.05))
        ax.grid(axis='y',alpha=.2)
    axes[0].set_ylabel('Fraction segregating at present')
    fig.suptitle('Acceptance under selection from 50 kya',y=.98,fontsize=19)
    fig.tight_layout(rect=(0,0,1,.92))
    for ext in ('png','pdf'):fig.savefig(out/f'acceptance.{ext}',dpi=160)
    plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',required=True,type=Path);p.add_argument('--out',required=True,type=Path)
    p.add_argument('--neutral-reference',type=Path)
    a=p.parse_args();report(a.root,a.out,a.neutral_reference)
