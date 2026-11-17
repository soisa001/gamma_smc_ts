#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = ["truvari==5.5.0", "google-crc32c==1.8.0"]
# ///
"""Audit original -> native retained -> final GT, and write separate locus trials.

No existing callset is changed. Trial VCFs contain only the selected site's
simple-SV representatives, GT/PS and recomputed AC/AN/AF (not pooled likelihoods).
"""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import importlib.metadata
import inspect
import json
import logging
from pathlib import Path
import sys
from types import SimpleNamespace

import pysam
import truvari
import truvari.collapse as native

import run_chr1_sv_merge as pilot
from phased_sv_union import snapshot, union_calls, apply_union


def text_gt(call):
    return ('|' if call['phased'] else '/').join('.' if a is None else str(a) for a in call['gt'])


def old_key(m):
    o=m['old']; return o['chrom'],o['original_pos'],o['original_id'],o['original_ref'],o['original_alt']


def new_key(m):
    return m['old']['chrom'],m['new_pos'],m['new_id'],m['new_ref'],m['new_alt']


def fetch_keys(path, keys):
    result = {}
    chrom = next(iter(keys))[0]
    with pysam.VariantFile(str(path)) as vcf:
        samples = list(vcf.header.samples)
        header = vcf.header.copy()
        lo,hi = min(k[1] for k in keys),max(k[1]+len(k[3]) for k in keys)
        for record in vcf.fetch(chrom,lo-1,hi):
            if len(record.alts or ())!=1:
                continue
            key = record.chrom,record.pos,record.id or '.',record.ref,record.alts[0]
            if key in keys:
                if key in result:
                    raise ValueError('Duplicate exact record identity')
                result[key] = record.copy()
    if set(result)!=set(keys):
        raise ValueError(f'Could not recover the requested alleles from {path}')
    return header,samples,result


def trial_header(source_header):
    header = pysam.VariantHeader()
    for name,contig in source_header.contigs.items():
        if contig.length:
            header.contigs.add(name,length=contig.length)
        else:
            header.contigs.add(name)
    header.formats.add('GT',1,'String','Haplotype union across cluster members; not dosage addition')
    ps_type = source_header.formats['PS'].type if 'PS' in source_header.formats else 'Integer'
    header.formats.add('PS',1,ps_type,'Shared source phase set when supported')
    for name,number,kind,desc in [
        ('AC','A','Integer','ALT allele count recomputed from GT'),
        ('AN',1,'Integer','Number of called alleles'),('AF','A','Float','AC/AN from GT'),
        ('SVTYPE',1,'String','SV class'),('SVLEN',1,'Integer','Signed net length'),
        ('PM_MEMBERS','.', 'String','Original working IDs in this cluster'),
        ('PM_UNRESOLVED',1,'Integer','Samples whose haplotype union is ambiguous or unsupported')]:
        header.info.add(name,number,kind,desc)
    for sample in source_header.samples:
        header.add_sample(sample)
    return header


def write_trial(path, header, groups, metadata, originals):
    stats, mapping = Counter(), []
    with pysam.VariantFile(str(path),'wz',header=header) as out:
        records = []
        for representative,members in groups:
            m = metadata[representative]['old']
            merged = len(members)>1
            pos,ref,alt = ((m['trimmed_pos'],m['trimmed_ref'],m['trimmed_alt']) if merged else
                           (m['original_pos'],m['original_ref'],m['original_alt']))
            record = header.new_record(contig=m['chrom'],start=pos-1,stop=pos-1+len(ref),
                                       id=m['original_id'],alleles=(ref,alt))
            record.info['SVTYPE'] = m['kind']; record.info['SVLEN'] = len(alt)-len(ref)
            record.info['PM_MEMBERS'] = tuple(members)
            source_records = [originals[ident] for ident in members]
            if merged:
                statuses = apply_union(record,source_records)
            else:
                statuses = {}
                for sample in record.samples:
                    call=snapshot(source_records[0].samples[sample])
                    record.samples[sample]['GT']=call['gt']; record.samples[sample].phased=call['phased']
                    record.samples[sample]['PS']=call['ps']
                    statuses[sample]='unchanged'
            unresolved = sum(v in ('phase_set_conflict','unphased_union_ambiguous','unsupported_gt') for v in statuses.values())
            record.info['PM_UNRESOLVED']=unresolved
            stats.update(statuses.values())
            ac = sum(tuple(s['GT']).count(1) for s in record.samples.values())
            an = sum(sum(a is not None for a in s['GT']) for s in record.samples.values())
            record.info['AC']=(ac,); record.info['AN']=an; record.info['AF']=(ac/an if an else None,)
            records.append(record)
            for ident in members:
                mapping.append(dict(old_working_id=ident,old_id=metadata[ident]['old']['original_id'],
                    representative_working_id=representative,new_id=record.id,new_pos=record.pos,
                    new_length=abs(len(alt)-len(ref)),members=len(members)))
        for record in sorted(records,key=lambda r:(r.chrom,r.pos,r.id)):
            out.write(record)
    pysam.tabix_index(str(path),preset='vcf',force=True)
    pilot.tsv(Path(str(path)+'.mapping.tsv'),mapping,list(mapping[0]))
    with pysam.VariantFile(str(path)) as reread:
        if sum(1 for _ in reread)!=len(groups):
            raise ValueError('Trial round-trip count mismatch')
    return dict(records=len(groups),sample_cluster_statuses=dict(stats),
                format_policy='GT/PS only; original per-allele quality/depth/likelihood fields are not pooled')


def verify_existing_trial(path, groups, metadata, originals, final_records, aliases):
    _, samples, repaired = fetch_keys(path, {new_key(metadata[r]) for r in groups})
    mismatches = comparisons = changed = 0
    focus = 'unverified: focused plot alias or target cluster unavailable'
    for representative, members in groups.items():
        key = new_key(metadata[representative])
        fixed_record = repaired[key]
        if set(fixed_record.info['PM_MEMBERS']) != set(members):
            raise ValueError('Trial VCF membership differs from saved mapping')
        donor = next((originals[i] for i in members if originals[i].id == '4_92312'), None)
        for sample in samples:
            calls = [snapshot(originals[i].samples[sample]) for i in members]
            expected = union_calls(calls) if len(members) > 1 else calls[0]
            actual = snapshot(fixed_record.samples[sample])
            old = snapshot(final_records[key].samples[sample])
            comparisons += 1
            mismatches += any(actual[k] != expected[k] for k in ('gt', 'phased', 'ps'))
            changed += (actual['gt'], actual['phased']) != (old['gt'], old['phased'])
            if key[2] == '4_92311' and donor is not None and 'S02.h1' in aliases.get(sample, ()):
                source = snapshot(donor.samples[sample])
                if source['phased'] and source['gt'][0] == 1 and old['gt'] == (0, 0):
                    focus = 'restored' if actual['phased'] and actual['gt'][0] == 1 else 'failed'
                else:
                    focus = 'unverified: source or prior GT differs from reported case'
    return dict(sample_cluster_comparisons=comparisons, gt_phase_ps_mismatches=mismatches,
                changed_gt_or_phase=changed, focused_h1=focus)


def off_groups(originals, metadata, out):
    args = native.parse_args(['-i',str(out/'unused.gz'),'-p','.7','-P','.7','-r','500','-O','0',
                              '-s','50','-S','-1','--gt','off','--keep','common'])
    params=truvari.VariantParams(args=args,short_circuit=True,skip_gt=True,sizefilt=50,chunksize=500)
    params.keep,params.hap,params.gt,params.chain='common',False,'off',False
    params.sorter,params.consolidate=native.SORTS['common'],False
    entries=[]
    for ident,record in originals.items():
        m=metadata[ident]['old']; r=record.copy()
        r.id=ident; r.pos=m['trimmed_pos']; r.alleles=(m['trimmed_ref'],m['trimmed_alt'])
        r.stop=r.start+len(r.ref)
        entries.append(truvari.VariantRecord(r,params))
    calls=native.collapse_chunk(({'base':entries,'__filtered':[]},1),params)
    return [(c.entry.id,[c.entry.id,*[m.comp.id for m in c.matches]]) for c in calls]


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch-root',type=Path,default=Path.home()/'popped_full_same_original_pos_p70')
    parser.add_argument('--pos',type=int,default=38754377)
    parser.add_argument('--kind',choices=['INS','DEL'],default='INS')
    parser.add_argument('--out',type=Path)
    parser.add_argument('--plot-audit-dir',type=Path)
    parser.add_argument('--existing-only', action='store_true', help='Use only saved cluster membership')
    opts=parser.parse_args(argv)
    if importlib.metadata.version('truvari')!='5.5.0':
        raise ValueError('Requires Truvari 5.5.0 for native replay')
    base=opts.batch_root.resolve()/'chr4'
    manifest=json.loads((base/'final/done.json').read_text())
    if manifest['key']['policy']['gt']!='all':
        raise ValueError('Expected completed --gt all run for this audit')
    source=Path(manifest['key']['source'][0]['path'])
    for identity in manifest['key']['source']:
        if pilot.common.identity(Path(identity['path']))!=identity:
            raise ValueError('Original input changed since the merge')
    mapping_path=base/'final/old_to_new.jsonl.gz'
    if pilot.common.digest(mapping_path)!=manifest['artifacts'][mapping_path.name]:
        raise ValueError('Corrupt old-to-new mapping')
    metadata={}
    with gzip.open(mapping_path,'rt') as handle:
        for m in map(json.loads,handle):
            if m['old']['original_pos']==opts.pos and m['old']['kind']==opts.kind:
                metadata[m['old']['working_id']]=m
    if not metadata:
        raise ValueError('No candidate alleles at the requested original position/type')
    header,samples,source_records=fetch_keys(source,{old_key(m) for m in metadata.values()})
    originals={ident:source_records[old_key(m)] for ident,m in metadata.items()}
    groups=defaultdict(list)
    for ident,m in metadata.items():
        groups[m['representative_working_id']].append(ident)
    final=base/'final/chr4.popped.samepos.p70.vcf.gz'
    _,final_samples,final_records=fetch_keys(final,{new_key(metadata[r]) for r in groups})
    native_path=base/'merge'/opts.kind/'retained.vcf.gz'
    native_keys={(m['old']['chrom'],m['old']['trimmed_pos'],ident,m['old']['trimmed_ref'],m['old']['trimmed_alt'])
                 for ident,m in metadata.items() if ident in groups}
    _,native_samples,native_records=fetch_keys(native_path,native_keys)
    if samples!=final_samples or samples!=native_samples:
        raise ValueError('Sample identity/order mismatch')
    native_by_id={r.id:r for r in native_records.values()}
    out=(opts.out or opts.batch_root/f'tlr_gt_transfer_audit_{opts.pos}').resolve()
    out.mkdir(parents=True,exist_ok=True)
    logging.basicConfig(filename=out/'audit.log',level=logging.INFO,force=True)
    logging.info('PARAMETERS %s', json.dumps(vars(opts), default=str, sort_keys=True))
    logging.info('TRIAL_POLICY p=0.7 P=0.7 refdist=500 overlap=0 sizemin=50 sizemax=-1 keep=common gt=off chain=false seed=%s', pilot.SEED)
    import fcntl
    with open(out/'.audit.lock','a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        (out/'summary.json').unlink(missing_ok=True)
        tracked=[source,final,native_path,mapping_path,base/'final/done.json']
        identities=[pilot.common.identity(p) for p in tracked]
        logging.info('INPUTS %s', json.dumps(identities, sort_keys=True))
        plot_dir=opts.plot_audit_dir or opts.batch_root/'tlr_merge_carrier_examples'
        aliases=defaultdict(list)
        selection=plot_dir/'selected_haplotypes.tsv'
        if selection.is_file():
            with selection.open() as handle:
                for row in csv.DictReader(handle,delimiter='\t'):
                    aliases[row['sample_id']].append(row['label'])
        totals,detail=Counter(),[]
        for representative,members in groups.items():
            first=originals[representative]
            donors=[originals[i] for i in members if i!=representative]
            replay,_=native.gt_aware_consolidate(first.copy(),[SimpleNamespace(comp=r.copy()) for r in donors])
            stored=native_by_id[representative]; assembled=final_records[new_key(metadata[representative])]
            for sample in samples:
                original_calls=[snapshot(originals[i].samples[sample]) for i in members]
                expected=union_calls(original_calls)
                n,f,r=[snapshot(record.samples[sample]) for record in (stored,assembled,replay)]
                totals['sample_cluster_comparisons']+=1
                multi=sum(1 in c['gt'] for c in original_calls)>1
                totals['existing_cluster_samples_with_multiple_carried_members']+=multi
                totals['native_final_gt_phase_ps_mismatches']+=n!=f
                totals['native_replay_gt_phase_ps_mismatches']+=n!=r
                comparable=expected['status'] not in ('phase_set_conflict','unphased_union_ambiguous','unsupported_gt')
                lost=comparable and None not in f['gt'] and expected['gt'].count(1)>f['gt'].count(1)
                totals['samples_losing_known_ALT_dosage_in_native']+=comparable and None not in n['gt'] and expected['gt'].count(1)>n['gt'].count(1)
                totals['samples_losing_known_ALT_dosage']+=lost
                totals['samples_with_missing_final_GT']+=None in f['gt']
                phase_changed=comparable and expected['phased'] and f['phased'] and expected['gt'].count(1)==f['gt'].count(1) and expected['gt']!=f['gt']
                totals['samples_with_changed_phase_slots_same_dosage']+=phase_changed
                if any(1 in c['gt'] for c in original_calls) or n!=f or n!=r:
                    detail.append(dict(sample_id=sample,plot_labels=';'.join(aliases[sample]),representative_id=first.id,cluster_members=len(members),
                        source_GTs=';'.join(f'{originals[i].id}={text_gt(snapshot(originals[i].samples[sample]))}' for i in members),
                        expected_GT=text_gt(expected),native_retained_GT=text_gt(n),final_GT=text_gt(f),replayed_native_GT=text_gt(r),
                        expected_PS=expected['ps'],native_PS=n['ps'],final_PS=f['ps'],
                        expected_status=expected['status'],lost_ALT_dosage=lost,phase_slots_changed=phase_changed,
                        native_equals_final=n==f,native_equals_replay=n==r))
        if detail:
            with gzip.open(out/'sample_gt_audit.tsv.gz','wt',newline='') as handle:
                writer=csv.DictWriter(handle,fieldnames=list(detail[0]),delimiter='\t');writer.writeheader();writer.writerows(detail)
        # Check the exact plotted GT strings against source/final independently.
        plot_dir=opts.plot_audit_dir or opts.batch_root/'tlr_merge_carrier_examples'
        plot_file=plot_dir/'haplotype_genotypes.tsv.gz'
        plot_counts=Counter()
        if plot_file.is_file():
            with gzip.open(plot_file,'rt') as handle:
                for row in csv.DictReader(handle,delimiter='\t'):
                    key=(row['chrom'],int(row['pos']),row['id'],row['ref'],row['alt'])
                    records=source_records if row['stage']=='Unmerged' else final_records
                    if key not in records:
                        continue
                    actual=snapshot(records[key].samples[row['sample_id']])
                    plot_counts['checked_genotype_rows']+=1
                    plot_counts['GT_string_mismatches']+=row['gt']!=text_gt(actual)
        legacy_groups=[(r,[r,*[i for i in members if i!=r]]) for r,members in groups.items()]
        th=trial_header(header)
        same=write_trial(out/'existing_clusters.gt_union.vcf.gz',th,legacy_groups,metadata,originals)
        checked=verify_existing_trial(out/'existing_clusters.gt_union.vcf.gz',groups,metadata,originals,final_records,aliases)
        if checked['gt_phase_ps_mismatches']:
            raise ValueError(f"Repaired GT/phase/PS mismatches: {checked['gt_phase_ps_mismatches']}")
        intra=None
        if not opts.existing_only:
            intra_groups=off_groups(originals,metadata,out)
            intra=write_trial(out/'intrasample.gt_union.vcf.gz',th,intra_groups,metadata,originals)
        if identities!=[pilot.common.identity(p) for p in tracked]:
            raise ValueError('Inputs changed during audit')
        summary=dict(position=opts.pos,type=opts.kind,samples=len(samples),original_alleles=len(metadata),
            existing_clusters=len(groups),existing_run_audit=dict(totals),plot_audit=dict(plot_counts),
            existing_cluster_union_trial=same,intrasample_union_trial=intra,existing_cluster_verification=checked,
            input_snapshots=identities,seed=pilot.SEED,
            phase_assumption='Same original-position phased GT slots share panel orientation; conflicting explicit PS values yield missing trial GT and a flag',
            native_version=importlib.metadata.version('truvari'),
            native_function_sha256=hashlib.sha256(inspect.getsource(native.gt_aware_consolidate).encode()).hexdigest(),
            script_sha256=pilot.common.digest(__file__),originals_modified=False)
        pilot.common.json_write(out/'summary.json',summary)
        print(json.dumps({k:summary[k] for k in ('position','type','samples','original_alleles','existing_clusters',
            'existing_run_audit','plot_audit','existing_cluster_union_trial','existing_cluster_verification','intrasample_union_trial')},indent=2))
        print(f'Private genotype audit and separate locus trials: {out}')


if __name__=='__main__':
    main()
