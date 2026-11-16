"""Private, read-only genotype plots for the completed popped TLR comparison.

The original plot's lane table supplies sample identity and called introgressed
tracts. Genotypes always come from the two real VCFs, never the identity map.
"""
import os
os.environ['MPLBACKEND'] = 'Agg'

from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime, timezone
import gzip
import logging

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import pysam

import compare_tlr_popped_merged as report


def table(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def lane_inputs(folder, window):
    path = folder/'popped_selected_lanes.tsv.gz'
    if not path.is_file():
        raise ValueError('Haplotype plotting requires popped_selected_lanes.tsv.gz in --plot-dir')
    rows = table(path)
    required = {'sample_id','hap','chrom','start','end'}
    if not rows or not required.issubset(rows[0]):
        raise ValueError('Existing popped lane table lacks sample/haplotype/tract coordinates')
    tracts = defaultdict(list)
    samples = set()
    for row in rows:
        if row['chrom'] != window[0]:
            continue
        sample, hap = row['sample_id'], int(row['hap'])
        start, end = int(row['start']), int(row['end'])
        if not sample or hap not in (1, 2) or start < 0 or end <= start:
            raise ValueError('Invalid sample/haplotype/tract in popped lane table')
        samples.add(sample)
        start, end = max(start, window[1]-1), min(end, window[2])
        if start < end:
            tracts[sample, hap].append((start,end))
    if not samples:
        raise ValueError('No chr4 samples in existing popped lane table')
    # Both phase slots are candidates; only explicitly annotated tracts receive
    # shading. Absence in this selected-lane table does not prove no introgression.
    candidates = [(sample,hap) for sample in sorted(samples) for hap in (1,2)]
    return path, candidates, tracts


def gt_state(gt, phased, allele, hap):
    """ALT/ref/unknown for a phase slot; never place unphased heterozygotes."""
    if len(gt)!=2:
        return 'unknown'
    if not phased and not (gt[0] is not None and gt[0]==gt[1]):
        return 'unknown'
    value = gt[hap-1]
    return 'unknown' if value is None else 'alt' if value==allele else 'other'


def read_calls(path, rows, samples, window, threads):
    wanted = {report.key(r) for r in rows}
    calls = {}
    with pysam.VariantFile(str(path), threads=threads) as vcf:
        if not set(samples).issubset(vcf.header.samples):
            raise ValueError('Selected plot samples are absent from a VCF; check plot input provenance')
        vcf.subset_samples(samples)
        for rec in vcf.fetch(window[0], window[1]-1, window[2]):
            for allele,alt in enumerate(rec.alts or (), 1):
                key = (rec.chrom,rec.pos,rec.id or '.',rec.ref,alt)
                if key not in wanted:
                    continue
                if key in calls:
                    raise ValueError('Duplicated exact allele identity prevents unambiguous haplotype plotting')
                calls[key] = {}
                for sample in samples:
                    value = rec.samples[sample]
                    gt = tuple(value.get('GT') or ())
                    for hap in (1,2):
                        calls[key][sample,hap] = dict(
                            state=gt_state(gt,value.phased,allele,hap),
                            gt=('|' if value.phased else '/').join('.' if x is None else str(x) for x in gt),
                            ps=str(value.get('PS', '.')),phased=value.phased)
    if set(calls)!=wanted:
        raise ValueError('Genotype extraction did not exactly recover the requested regional alleles')
    return calls


def clusters_for(before, after, mappings):
    bkeys, akeys = set(map(report.key,before)), set(map(report.key,after))
    grouped = defaultdict(list)
    for m in mappings:
        if m['members']>1:
            grouped[m['representative_working_id']].append(m)
    clusters = []
    for members in grouped.values():
        old = sorted({report.old_key(m) for m in members})
        dest = report.key(report.target(members[0]))
        if dest in akeys and any(k in bkeys for k in old):
            if not set(old).issubset(bkeys):
                raise ValueError('A merged cluster crosses the TLR boundary; widen the plot before comparing phase slots')
            clusters.append(dict(old=old, new=dest))
    clusters.sort(key=lambda c:(c['old'][0][1],c['old'][0]))
    for i,c in enumerate(clusters,1):
        c['label'] = f'C{i}'
    return clusters


def select_haps(candidates, tracts, calls, clusters, limit):
    features = {h:{k for c in clusters for k in c['old'] if calls[k][h]['state']=='alt'}
                for h in candidates}
    total = {h:sum(r[h]['state']=='alt' for r in calls.values()) for h in candidates}
    # Prioritize introgressed carriers of merged alleles, then other carriers,
    # then introgressed context lanes. Stable ties; no stochastic sampling.
    useful = [h for h in candidates if total[h] or tracts.get(h)]
    selected, covered = [], set()
    while useful and len(selected)<limit:
        h = min(useful, key=lambda h:(-bool(tracts.get(h) and features[h]),
                -bool(features[h]),-bool(tracts.get(h)),
                -len(features[h]-covered),-len(features[h]),-total[h],h))
        selected.append(h); covered |= features[h]; useful.remove(h)
    return selected, features


def source_union(cluster, calls, hap):
    states = [calls[k][hap]['state'] for k in cluster['old']]
    return 'alt' if 'alt' in states else 'unknown' if 'unknown' in states else 'other'


def include_phase_partners(ranked, before, after, clusters, limit):
    """Show the other phase slot when merging changes a sample's ALT placement."""
    changed_samples = set()
    for h in ranked:
        for c in clusters:
            old, new = source_union(c,before,h), after[c['new']][h]['state']
            if 'unknown' not in (old,new) and old!=new:
                changed_samples.add(h[0])
    selected = []
    for h in ranked:
        if h in selected:
            continue
        unit = [h]
        partner = (h[0],3-h[1])
        if h[0] in changed_samples and partner not in selected:
            unit.append(partner)
        if len(selected)+len(unit)<=limit:
            selected.extend(unit)
    return selected


def write_table(path, rows, fields):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path,'wt',newline='') as handle:
        writer = csv.DictWriter(handle,fieldnames=fields,delimiter='\t')
        writer.writeheader(); writer.writerows(rows)


def render(out, window, selected, tracts, stages, clusters, changed=0):
    plt.rcParams.update({'font.size':18,'axes.titlesize':23,'axes.labelsize':20,
                         'pdf.fonttype':42,'ps.fonttype':42})
    colors = {'INS':'#D55E00','DEL':'#0072B2','complex_replacement':'#8E4585',
              'ambiguous_sequence':'#777777'}
    shapes = {'INS':'^','DEL':'v','complex_replacement':'D','ambiguous_sequence':'s'}
    cluster_lookup = {}
    for c in clusters:
        for key in [*c['old'],c['new']]:
            cluster_lookup[key] = c['label']
    sample_ids = {sample:i+1 for i,sample in enumerate(dict.fromkeys(h[0] for h in selected))}
    labels = [f'S{sample_ids[h[0]]:02d} / h{h[1]}' + (' *' if tracts.get(h) else '') for h in selected]
    fig,axes = plt.subplots(2,1,figsize=(15,max(10,4+.95*len(selected))),sharex=True)
    fig.subplots_adjust(left=.18,right=.98,bottom=.20,top=.78,hspace=.34)
    fig.suptitle('TLR locus · same haplotypes before and after merging',y=.985,fontsize=24)
    handles = [Patch(color='#D5DED6',label='Called introgressed tract')]
    for name in ('INS','DEL','complex_replacement','ambiguous_sequence'):
        if name in ('INS','DEL') or any(report.allele_class(r)==name for _,rows,_ in stages for r in rows):
            handles.append(Line2D([],[],ls='',marker=shapes[name],color=colors[name],markersize=11,
                                  label=name.replace('_',' ').capitalize() if '_' in name else name))
    handles.append(Line2D([],[],ls='',marker='x',color='#555',markersize=9,label='GT placement unknown'))
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.56,.947),ncol=3,
               frameon=False,fontsize=16,columnspacing=1.5)
    xmin,xmax = (window[1]-1)/1e6,window[2]/1e6
    span = xmax-xmin
    for ax,(name,rows,calls) in zip(axes,stages):
        ax.set_title(name,loc='left',pad=18)
        ax.set_xlim(xmin,xmax); ax.set_ylim(len(selected)-.55,-.85)
        ax.set_yticks(range(len(selected)),labels,fontsize=17)
        ax.spines[['top','right','left']].set_visible(False)
        ax.tick_params(axis='y',length=0)
        for i,h in enumerate(selected):
            ax.hlines(i,xmin,xmax,color='#E7E7E7',lw=1,zorder=0)
            for start,end in tracts.get(h,[]):
                ax.plot([start/1e6,end/1e6],[i,i],color='#D5DED6',lw=13,
                        solid_capstyle='butt',zorder=1)
            annotations = []
            at_position = defaultdict(int)
            for row in sorted(rows,key=lambda r:(r['pos'],r['net_bp'],r['id'])):
                key = report.key(row); state = calls[key][h]['state']
                if state=='other':
                    continue
                desc,_ = report.rep.describe(row['ref'],row['alt'])
                x = (row['pos']-1+desc['shared_prefix_bp'])/1e6
                x = min(xmax,max(xmin,x))
                category = report.allele_class(row)
                style = row['kind'] if category=='simple' else category
                if state=='unknown':
                    ax.scatter(x,i,marker='x',s=30,color='#777',linewidths=1,zorder=2)
                    continue
                offset = .16*at_position[x]; at_position[x] += 1
                y = i+offset
                if row['kind']=='DEL':
                    end = (row['pos']-1+desc['shared_prefix_bp']+desc['residual_ref_bp'])/1e6
                    ax.plot([x,min(xmax,end)],[y,y],color=colors[style],lw=3,zorder=3)
                ax.scatter(x,y,marker=shapes[style],s=105,color=colors[style],
                           edgecolors='black' if key in cluster_lookup else 'white',
                           linewidths=1,zorder=4)
                if key in cluster_lookup:
                    annotations.append((x,y,f"{cluster_lookup[key]}:{row['net_bp']}"))
            # Labels show cluster and bp, not sequence. Space them horizontally,
            # retaining leader lines to exact genomic positions.
            annotations.sort()
            gap = min(span*.105,span*.9/max(1,len(annotations)-1))
            positions = [min(xmax-span*.045,max(xmin+span*.045,a[0])) for a in annotations]
            for j in range(1,len(positions)):
                positions[j] = max(positions[j],positions[j-1]+gap)
            if positions:
                positions[-1] = min(positions[-1],xmax-span*.045)
                for j in range(len(positions)-2,-1,-1):
                    positions[j] = min(positions[j],positions[j+1]-gap)
            for (x,y,label),label_x in zip(annotations,positions):
                ax.annotate(label,(x,y),xytext=(label_x,i-.31),ha='center',va='bottom',fontsize=13,
                            arrowprops=dict(arrowstyle='-',lw=.6,color='#555'),zorder=5)
        ax.grid(axis='x',alpha=.12)
    axes[-1].set_xlabel('Chromosome 4 position (Mb)',labelpad=12)
    axes[-1].ticklabel_format(axis='x',style='plain',useOffset=False)
    fig.text(.18,.025,'* Introgression documented in saved lanes. C#:bp identifies merged-cluster alleles.\n'
             'GT phase slots are read independently; blank tracts do not establish ancestry.\n'
             f'{changed} cluster × haplotype ALT placements changed; see the placement audit.',
             fontsize=14,ha='left')
    for suffix in ('png','pdf'):
        path = out/f'tlr_haplotypes_before_after.{suffix}'
        temp = path.with_name(path.stem+'.partial.'+suffix)
        fig.savefig(temp,dpi=200,facecolor='white')
        os.replace(temp,path)
    plt.close(fig)


def build(source, merged, before, after, mappings, folder, out, window, limit=8, threads=4):
    started = datetime.now(timezone.utc).isoformat()
    if not 2<=limit<=20 or threads<1:
        raise ValueError('Use 2–20 haplotypes and a positive thread count')
    lane_path,candidates,tracts = lane_inputs(folder,window)
    lane_identity = report.pilot.common.identity(lane_path)
    lane_sha = report.pilot.common.digest(lane_path)
    b,a = [[r for r in rows if report.inside(r,window)] for rows in (before,after)]
    samples = sorted({h[0] for h in candidates})
    logging.info('Haplotype plot: window=%s samples=%d limit=%d threads=%d lanes=%s',
                 window,len(samples),limit,threads,lane_path)
    with ThreadPoolExecutor(max_workers=min(2,threads)) as pool:
        futures = [pool.submit(read_calls,path,rows,samples,window,max(1,threads//2))
                   for path,rows in ((source,b),(merged,a))]
        bcalls,acalls = [f.result() for f in futures]
    clusters = clusters_for(b,a,mappings)
    cluster_labels = {key:c['label'] for c in clusters for key in [*c['old'],c['new']]}
    ranked,features = select_haps(candidates,tracts,bcalls,clusters,len(candidates))
    selected = include_phase_partners(ranked,bcalls,acalls,clusters,limit)
    if not selected:
        raise ValueError('No informative TLR haplotypes among saved plot samples')
    sample_ids = {sample:i+1 for i,sample in enumerate(dict.fromkeys(h[0] for h in selected))}
    ids = {h:f'S{sample_ids[h[0]]:02d}.h{h[1]}' for h in selected}
    selection_rows = [dict(label=ids[h],sample_id=h[0],hap=h[1],
        introgression_documented=bool(tracts.get(h)),merged_source_alleles=len(features[h])) for h in selected]
    write_table(out/'selected_haplotypes.tsv',selection_rows,list(selection_rows[0]))
    states,phase = [],[]
    stages = [('Unmerged',b,bcalls),('Merged · 70% same original position',a,acalls)]
    for stage,rows,calls in stages:
        for row in rows:
            for h in selected:
                call = calls[report.key(row)][h]
                states.append(dict(stage=stage,label=ids[h],sample_id=h[0],hap=h[1],
                    chrom=row['chrom'],pos=row['pos'],id=row['id'],ref=row['ref'],alt=row['alt'],
                    cluster=cluster_labels.get(report.key(row),''),
                    allele_class=report.allele_class(row),type=row['kind'],length_bp=row['net_bp'],**call))
    for c in clusters:
        for h in selected:
            expected = source_union(c,bcalls,h)
            actual = acalls[c['new']][h]['state']
            phase.append(dict(cluster=c['label'],label=ids[h],sample_id=h[0],hap=h[1],
                merged_id=c['new'][2],merged_pos=c['new'][1],
                source_union=expected,merged_state=actual,
                comparison='unknown' if 'unknown' in (expected,actual) else 'match' if expected==actual else 'changed'))
    write_table(out/'haplotype_genotypes.tsv.gz',states,
                list(states[0]) if states else ['stage','label','state'])
    write_table(out/'haplotype_placement_audit.tsv',phase,
                list(phase[0]) if phase else ['cluster','label','comparison'])
    changed = sum(r['comparison']=='changed' for r in phase)
    render(out,window,selected,tracts,stages,clusters,changed)
    if report.pilot.common.identity(lane_path)!=lane_identity:
        raise ValueError('Introgression lane table changed during plotting')
    paths = ['selected_haplotypes.tsv','haplotype_genotypes.tsv.gz','haplotype_placement_audit.tsv',
             'tlr_haplotypes_before_after.png','tlr_haplotypes_before_after.pdf']
    result = dict(status='completed',started_utc=started,finished_utc=datetime.now(timezone.utc).isoformat(),
        window_1based_inclusive=window,lanes=lane_identity,lanes_sha256=lane_sha,
        selection='Saved popped-plot samples; introgressed merged-allele carriers, other merged-allele carriers, then context lanes; favor distinct source alleles, deterministic ties. Include partner phase slot when ALT placement changes.',
        max_haplotypes=limit,threads=threads,seed=report.pilot.SEED,
        phase_assumption='VCF GT hap1/hap2 matches original decoder orientation. No rephasing or proof of cross-block phase consistency. PS is saved per call.',
        unknown='Unphased heterozygotes, missing phase slots and non-diploid GTs are unassigned; homozygotes are phase-independent.',
        aggregate=dict(selected_haplotypes=len(selected),
            selected_with_documented_introgression=sum(bool(tracts.get(h)) for h in selected),
            merged_clusters_in_TLR=len(clusters),
            merged_clusters_shown=sum(any(bcalls[k][h]['state']=='alt' for k in c['old'] for h in selected) for c in clusters),
            changed_cluster_haplotype_placements=sum(r['comparison']=='changed' for r in phase),
            unknown_cluster_haplotype_comparisons=sum(r['comparison']=='unknown' for r in phase)),
        artifacts={p:report.pilot.common.digest(out/p) for p in paths},
        versions=dict(pysam=pysam.__version__,matplotlib=matplotlib.__version__),
        script_sha256=report.pilot.common.digest(__file__))
    report.pilot.common.json_write(out/'haplotype_plot_summary.json',result)
    return result
