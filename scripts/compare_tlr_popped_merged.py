#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = ["google-crc32c==1.8.0"]
# ///
"""Read-only TLR comparison of original popped panel and completed full merge.

Only aggregate counts go to stdout. Individual identities/crosswalks stay local.
"""
import argparse
from collections import Counter, defaultdict, deque
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
import json
import logging
from pathlib import Path
import sys

import inspect_sv_representation as rep
import run_chr1_sv_merge as pilot

WINDOWS = {'plot_window':('chr4',38660000,38960000),
           'TLR_interval':('chr4',38780001,38890001)}


def key(row):
    return row['chrom'], int(row['pos']), row['id'], row['ref'], row['alt']


def inside(row, window):
    chrom, start, end = window
    return row['chrom']==chrom and row['pos']<=end and row['pos']+len(row['ref'])-1>=start


def read_sites(args, path, label):
    rows = []
    cmd = [args.bcftools,'query','-r','chr4:38660000-38960000','--regions-overlap','1',
           '-i','ILEN>=50 || ILEN<=-50','-f','%CHROM\t%POS\t%ID\t%REF\t%ALT\n',path]
    with pilot.stream(args,cmd,label) as lines:
        for line in lines:
            chrom, pos, ident, ref, alts = line.rstrip('\n').split('\t')
            for alt in alts.split(','):
                desc, _ = rep.describe(ref,alt)
                if not desc or desc['type'] not in ('INS','DEL') or desc['net_bp']<50:
                    continue
                row = dict(chrom=chrom,pos=int(pos),id=ident,ref=ref,alt=alt,
                           kind=desc['type'],net_bp=desc['net_bp'],
                           biallelic=',' not in alts,
                           simple=not desc['complex'] and not desc['ambiguous'])
                if inside(row,WINDOWS['plot_window']):
                    rows.append(row)
    return rows


def target(mapping):
    return dict(chrom=mapping['old']['chrom'],pos=int(mapping['new_pos']),
                id=mapping['new_id'],ref=mapping['new_ref'],alt=mapping['new_alt'])


def old_key(mapping):
    old = mapping['old']
    return old['chrom'],old['original_pos'],old['original_id'],old['original_ref'],old['original_alt']


def check_mapping(before, after, mappings):
    """Verify regional after alleles against source + completed old-to-new map.

    Count distinct representatives once, while preserving unrelated duplicate
    source records. Trimming across a window boundary is explicitly accounted for.
    """
    lookup = defaultdict(deque)
    for item in mappings:
        lookup[old_key(item)].append(item)
    expected, details = {}, []
    for number,row in enumerate(before):
        found = lookup.get(key(row))
        match = found.popleft() if found else None
        if row['simple'] and row.get('biallelic',True) and match is None:
            raise ValueError('A simple source SV is missing from the completed old-to-new map')
        dest = target(match) if match else row
        token = ('mapped',match['representative_working_id']) if match else ('unchanged',number)
        expected[token] = dest
        details.append(dict(old_id=row['id'],old_pos=row['pos'],type=row['kind'],old_length=row['net_bp'],
            new_id=dest['id'],new_pos=dest['pos'],new_length=abs(len(dest['alt'])-len(dest['ref'])),
            status=match['status'] if match else 'unchanged_passthrough',
            cluster_members=match['members'] if match else 1))
    expected_counts = Counter(key(r) for r in expected.values() if inside(r,WINDOWS['plot_window']))
    if expected_counts != Counter(map(key,after)):
        raise ValueError('Regional merged VCF disagrees with the completed old-to-new mapping; inspect inputs')
    return details


def deleted_union(rows, window):
    intervals = []
    for row in rows:
        if row['kind']!='DEL' or not row['simple']:
            continue
        desc, _ = rep.describe(row['ref'],row['alt'])
        start = row['pos']-1+desc['shared_prefix_bp']
        end = start+desc['residual_ref_bp']
        start,end = max(start,window[1]-1),min(end,window[2])
        if start<end:
            intervals.append((start,end))
    total, previous_end = 0, -1
    for start,end in sorted(intervals):
        total += max(0,end-max(start,previous_end))
        previous_end = max(previous_end,end)
    return total


def summarize(before, after):
    counts, bp = [], []
    for scope,window in WINDOWS.items():
        b,a = ([r for r in rows if inside(r,window)] for rows in (before,after))
        for kind in ('INS','DEL','TOTAL'):
            left,right = ([r for r in rows if kind=='TOTAL' or r['kind']==kind] for rows in (b,a))
            n0,n1 = len(left),len(right)
            counts.append(dict(scope=scope,type=kind,before=n0,after=n1,removed=n0-n1,
                percent_removed=round(100*(n0-n1)/n0,3) if n0 else 0,
                simple_before=sum(r['simple'] for r in left),simple_after=sum(r['simple'] for r in right),
                unique_before=len({(r['chrom'],r['pos'],r['ref'],r['alt']) for r in left}),
                unique_after=len({(r['chrom'],r['pos'],r['ref'],r['alt']) for r in right})))
            if kind!='TOTAL':
                measures = {'raw_net_bp_sum':(sum(r['net_bp'] for r in left),sum(r['net_bp'] for r in right))}
                if kind=='DEL':
                    measures['unique_deleted_reference_bp_simple_only'] = (deleted_union(left,window),deleted_union(right,window))
                for metric,(v0,v1) in measures.items():
                    bp.append(dict(scope=scope,type=kind,metric=metric,before=v0,after=v1,
                                   reduction=v0-v1,percent_reduction=100*(v0-v1)/v0 if v0 else 0))
    return counts,bp


def plot_crosswalk(folder, before, mappings, after):
    """Map displayed popped identities; does not infer post-merge phased carriers."""
    allele_file = folder/'popped_alleles.tsv.gz'
    color_files = sorted(folder.glob('custom_TLR10_TLR1_TLR6_sv_alleles_*.tsv'))
    if not allele_file.is_file() or not color_files:
        return [],dict(status='plot_tables_not_found')
    displayed = {}
    for path in color_files:
        with open(path) as handle:
            for row in csv.DictReader(handle,delimiter='\t'):
                displayed[row['sv_event_id']] = row
    with gzip.open(allele_file,'rt') as handle:
        alleles = [r for r in csv.DictReader(handle,delimiter='\t') if r['sv_event_id'] in displayed]
    mapping_lookup = defaultdict(list)
    for m in mappings:
        mapping_lookup[old_key(m)].append(m)
    before_keys,after_keys = set(map(key,before)),set(map(key,after))
    results,new_keys = [],set()
    for row in alleles:
        src = dict(chrom=row['chrom'],pos=int(row['original_pos']),id=row['original_id'],ref=row['ref'],alt=row['alt'])
        if key(src) not in before_keys:
            raise ValueError('A plotted allele does not match the original popped input; check plot provenance')
        matches = mapping_lookup.get(key(src),[])
        if len(matches)>1:
            raise ValueError('Ambiguous duplicated plotted source identity')
        dest = target(matches[0]) if matches else src
        if key(dest) not in after_keys:
            raise ValueError('A plotted allele maps outside the queried merged window or is absent')
        new_keys.add(key(dest))
        color = displayed[row['sv_event_id']]
        results.append(dict(plot_label=color['allele_label'],marker_color=color['marker_color'],
            original_id=src['id'],original_pos=src['pos'],type=row['svtype'],original_length=row['sv_length_bp'],
            new_id=dest['id'],new_pos=dest['pos'],new_length=abs(len(dest['alt'])-len(dest['ref'])),
            status=matches[0]['status'] if matches else 'unchanged_passthrough',
            cluster_members=matches[0]['members'] if matches else 1))
    if len(results)!=len(displayed):
        raise ValueError('Plot allele table and marker/color table disagree')
    return results,dict(status='verified_identity_crosswalk',plotted_popped_alleles=len(results),
        distinct_mapped_merged_alleles=len(new_keys),
        note='Identity mapping only; post-merge GT carrier/phase placement is not assessed')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch-root',type=Path,default=Path.home()/'popped_full_same_original_pos_p70')
    parser.add_argument('--plot-dir',type=Path,default=Path.home()/'ashg_revised/23_shortlist_followup/custom_loci/TLR_popped/internal_only')
    parser.add_argument('--out',type=Path)
    parser.add_argument('--bcftools',default='bcftools')
    opts = parser.parse_args()
    base = opts.batch_root.resolve()/'chr4'; final = base/'final'
    prepared = json.loads((base/'prepare/done.json').read_text())
    completed = json.loads((final/'done.json').read_text())
    if completed['key']['prepared']!=prepared:
        raise ValueError('Completed chr4 final VCF does not match the preparation manifest')
    policy = completed['key']['policy']
    if policy['grouping']!='original_CHROM_POS_and_type' or policy['pctseq']!=.7 or policy['pctsize']!=.7:
        raise ValueError('Expected the original-POS/type p70 full merge')
    source = Path(prepared['key']['settings']['source'][0]['path'])
    merged = final/'chr4.popped.samepos.p70.vcf.gz'
    mapping = final/'old_to_new.jsonl.gz'
    for entry in prepared['key']['settings']['source']:
        if pilot.common.identity(Path(entry['path']))!=entry:
            raise ValueError('Original source/index changed since the completed merge')
    # Validate the small lineage artifact fully; scan indexed regional VCF data
    # rather than rehashing a potentially very large chromosome VCF.
    if pilot.common.digest(mapping)!=completed['artifacts'][mapping.name]:
        raise ValueError('Corrupt old-to-new mapping')
    out = (opts.out or opts.batch_root/'tlr_before_after').resolve()
    (out/'logs').mkdir(parents=True,exist_ok=True)
    args = argparse.Namespace(root=out,bcftools=opts.bcftools)
    logging.basicConfig(filename=out/'logs/comparison.log',level=logging.INFO,force=True)
    bcftools_version = pilot.command(args,[opts.bcftools,'--version'],'version',True)
    snapshots = [pilot.common.identity(p) for p in
        (*[Path(r['path']) for r in prepared['key']['settings']['source']],
         merged,Path(str(merged)+'.tbi'),mapping,final/'done.json')]
    import fcntl
    with open(out/'.compare.lock','a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        with ThreadPoolExecutor(max_workers=2) as pool:
            results = list(pool.map(lambda item:read_sites(args,*item),[(source,'popped'),(merged,'merged')]))
        before,after = results
        with gzip.open(mapping,'rt') as handle:
            mappings = []
            for line in handle:
                item = json.loads(line)
                old = item['old']
                row = dict(chrom=old['chrom'],pos=old['original_pos'],ref=old['original_ref'])
                if inside(row,WINDOWS['plot_window']):
                    mappings.append(item)
        details = check_mapping(before,after,mappings)
        counts,bp = summarize(before,after)
        plot_rows,plot_summary = plot_crosswalk(opts.plot_dir,before,mappings,after)
        for name,rows in [('counts.tsv',counts),('base_pairs.tsv',bp),('allele_changes.tsv',details),('plotted_allele_changes.tsv',plot_rows)]:
            if rows:
                pilot.tsv(out/name,rows,list(rows[0]))
            else:
                (out/name).unlink(missing_ok=True)
        if [pilot.common.identity(Path(r['path'])) for r in snapshots]!=snapshots:
            raise ValueError('Input changed during comparison')
        pilot.common.json_write(out/'summary.json',dict(source=str(source),merged=str(merged),policy=policy,
            windows_1based_inclusive=WINDOWS,counts=counts,plot=plot_summary,
            selection='Literal INS/DEL with abs(len(ALT)-len(REF))>=50; REF-span overlap, including anchor',
            validation='Source stat provenance, mapping SHA256, indexed regional read and exact lineage-to-VCF reconciliation',
            snapshots=snapshots,bcftools=bcftools_version,script_sha256=pilot.common.digest(__file__),seed=pilot.SEED))
        print((out/'counts.tsv').read_text(),end='')
        print(json.dumps(plot_summary,indent=2))
        print(f'All individual allele details remain on the VM: {out}')


if __name__=='__main__':
    main()
