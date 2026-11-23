#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = ["pysam==0.24.1", "google-crc32c==1.8.0"]
# ///
"""Combine repaired popped and same-start p70 lenient panels on the Workbench VM."""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from contextlib import ExitStack
import csv
import fcntl
import gzip
import heapq
import json
import logging
import os
from pathlib import Path
import re
import subprocess
import sys

import pysam
import run_chr1_sv_merge as pilot
from inspect_sv_representation import describe
from phase_singletons import digest, json_write

CUT = 10000
VERSION = 1


def load(path):
    return json.loads(path.read_text())


def verify(directory, name, hashes):
    path = directory/name
    if digest(path) != hashes[name]:
        raise ValueError(f'Input checksum mismatch: {path}')
    return path


def long_inputs(run):
    folder = run/'partitions'
    manifest = load(folder/'manifest.json')
    sig = manifest['input']['source']
    if Path(sig['path']).name != 'v3_ultralong_lenient.bcf' or manifest['input'].get('trim_version') != 1:
        raise ValueError('Expected the lenient same-start run')
    staged=Path(sig['path'])
    if staged.exists() and (staged.stat().st_size,staged.stat().st_mtime_ns)!=(sig['size'],sig['mtime_ns']):
        raise ValueError('Staged lenient BCF changed since recollapse')
    if any(manifest['skipped'].values()):
        raise ValueError('Long run skipped source variants; cannot represent this as a complete lenient panel')
    verify(folder,'source_records.tsv',manifest['hashes'])
    with (run/'comparison.tsv').open() as f:
        rows = [r for r in csv.DictReader(f,delimiter='\t') if float(r['threshold']) == .7]
    parts = []
    total = 0
    for row in rows:
        part = run/'pct70'/(re.sub(r'[^A-Za-z0-9_.-]', '_', row['svtype'])+'.'+row['mode'])
        cert = load(part/'complete.json')
        recipe = cert['recipe']
        if (recipe.get('same_start'),recipe['pctsize'],recipe['gt']) != (True,.7,'all'):
            raise ValueError(f'Wrong long-merge recipe: {part}')
        if any(int(row[k]) != cert['summary'][k] for k in ('input_records','retained_records','collapsed_records')):
            raise ValueError(f'Long completion counts disagree: {part}')
        for name in ('retained.vcf.gz','retained.vcf.gz.tbi','removed.vcf.gz'):
            verify(part,name,cert['hashes'])
        if recipe['input_sha256'] not in manifest['hashes'].values():
            raise ValueError(f'Long partition provenance mismatch: {part}')
        total += int(row['input_records'])
        parts.append(part)
    if (not rows or len(parts)!=len(set(parts)) or len(parts)!=len(manifest['partitions'])
            or total!=sum(p['records'] for p in manifest['partitions'].values())):
        raise ValueError('Incomplete p70 results')
    return parts


def info_dict(raw):
    return dict((v.split('=',1)+[None])[:2] for v in raw.split(';') if v and v!='.')


def profile(pos, ref, alt, tags, long_size=None):
    """Breakpoint is 0-based; ref intervals are half-open. size governs tier ownership."""
    if ',' in alt:
        raise ValueError('Expected biallelic panel records; a multiallelic row requires allele-specific routing')
    if alt in ('.','*','<NON_REF>','<NONREF>'):
        return None
    if alt.startswith('<') or any(c in alt for c in '[]') or alt.startswith('.') or alt.endswith('.'):
        typ = (tags.get('SVTYPE') or (alt.strip('<>') if alt.startswith('<') else 'BND')).split(':')[0]
        span = int(tags.get('END') or pos)-pos
        svlen = tags.get('SVLEN')
        size = long_size or (abs(int(svlen)) if svlen and svlen!='.' else (span if span>0 else None))
        interval = (pos,pos+span) if typ in ('DEL','DUP','INV') and span>0 and (not size or size==span) else None
        return ('symbolic_'+typ,size,interval,None,(pos,ref,alt,tags.get('END'),size))
    if max(len(ref),len(alt)) < 50 and long_size is None:
        return None
    d,_ = describe(ref,alt)
    if d is None or max(d['residual_ref_bp'],d['residual_alt_bp'])<50 and long_size is None:
        return None
    prefix = d['shared_prefix_bp']
    a = alt[prefix:prefix+d['residual_alt_bp']].upper()
    r = ref[prefix:prefix+d['residual_ref_bp']].upper()
    start = pos-1+prefix
    kind = 'complex' if r and a else 'simple_INS' if a else 'simple_DEL'
    if d['ambiguous']:
        kind = 'ambiguous_'+kind
    return (kind,long_size if long_size is not None else abs(len(a)-len(r)),
            (start,start+len(r)) if r else None, a or None, (start,r,a))


def union_bp(intervals):
    total = end = 0
    for a,b in sorted(intervals):
        total += max(0,b-max(a,end))
        end = max(end,b)
    return total


class Metrics:
    def __init__(self):
        self.pos=set(); self.keys=set(); self.intervals=[]; self.alt=set()
        self.counts=Counter(); self.bins=Counter()

    def add(self,pos,ref,alt,p,ac,an,n):
        kind,size,interval,payload,key = p
        self.pos.add(pos); self.keys.add(key)
        if interval:
            self.intervals.append(interval)
        if payload:
            self.alt.add((key[0],payload))
        c=self.counts
        c['allele_records']+=1
        c['AC_sum']+=ac; c['AN_sum']+=an; c['diploid_slots']+=2*n
        c['singletons_AC1']+=ac==1
        if kind.startswith('symbolic_'):
            c['symbolic_records']+=1
            c['symbolic_reported_length_sum']+=size or 0
            c['symbolic_DEL_without_verified_span']+=kind=='symbolic_DEL' and interval is None
            c['symbolic_INS_without_payload']+=kind=='symbolic_INS'
        af=ac/an if an else None
        label=('AN=0' if not an else 'AF=0' if not ac else '0<AF<0.001' if af<.001
               else '0.001<=AF<0.01' if af<.01 else '0.01<=AF<0.05' if af<.05
               else '0.05<=AF<0.5' if af<.5 else '0.5<=AF<1' if af<1 else 'AF=1')
        self.bins['called_AF',label]+=1
        ac_bin=('AN=0' if not an else '0' if not ac else '1' if ac==1 else '2' if ac==2
                else '3-5' if ac<=5 else '6-10' if ac<=10 else '11-100' if ac<=100 else '>100')
        self.bins['observed_AC',ac_bin]+=1

    def row(self):
        return dict(self.counts,positions=len(self.pos),unique_allele_keys=len(self.keys),
                    unique_reference_span_bp=union_bp(self.intervals),
                    unique_alt_payload_bp=sum(len(a) for _,a in self.alt))


def gt_counts(fmt, tail,n):
    if fmt.split(b':',1)[0]!=b'GT':
        raise ValueError('Expected leading GT on SV record')
    gts=re.sub(rb':[^\t]*',b'',tail) if b':' in fmt else tail
    if gts.count(b'\t')+1!=n or re.search(rb'[^01./|\t]',gts):
        raise ValueError('Invalid biallelic GT/sample columns')
    ac=gts.count(b'1')
    return ac, ac+gts.count(b'0')


def short_stream(handle,chrom):
    for i,line in enumerate(handle,1):
        if line.startswith(b'#'):
            continue
        c,pos,_=line.split(b'\t',2)
        if c.decode()!=chrom:
            raise ValueError('Unexpected chromosome in repaired panel')
        yield int(pos), 'popped', i, line


def long_stream(vcf,contig,chrom,order):
    for rec in vcf.fetch(contig):
        fields=str(rec).rstrip('\n').encode().split(b'\t',9)
        fields[0]=chrom.encode()
        if order is not None:
            raw=fields[9].split(b'\t')
            fields[9]=b'\t'.join(raw[j] for j in order)
        yield rec.pos,'lenient',rec.id,b'\t'.join(fields)


def write_tsv(path,rows):
    rows=list(rows)
    fields=list(dict.fromkeys(k for row in rows for k in row))
    tmp=path.with_suffix('.tmp')
    with tmp.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t',restval=0)
        writer.writeheader(); writer.writerows(rows)
    tmp.replace(path)


def short_receipts(args,chrom):
    """Only the completed repair and the original crosswalk are used in assembly."""
    repaired=args.repaired/chrom
    fix=load(repaired/'done.json')
    base=args.original/chrom/'final'
    old=load(base/'done.json')
    if fix['key']['original_artifacts']!=old['artifacts']:
        raise ValueError(f'{chrom}: repair and original merge differ')
    panel=repaired/f'{chrom}.popped.samepos.p70.gt_union.vcf.gz'
    verify(repaired,panel.name,fix['artifacts'])
    verify(repaired,'summary.json',fix['artifacts'])
    mapping=verify(base,'old_to_new.jsonl.gz',old['artifacts'])
    return fix,old,panel,mapping


def prepare_headers(short, vcfs, chrom):
    header=short.copy()
    samples=list(header.samples)
    if not samples or 'UP_KEY' in header.info:
        raise ValueError('Unexpected popped header')
    header.add_line('##INFO=<ID=UP_KEY,Number=1,Type=String,Description="Stable SV crosswalk key within this chromosome">')
    plans=[]
    for vcf in vcfs:
        h=vcf.header
        if set(h.samples)!=set(samples):
            raise ValueError('Popped and lenient sample sets differ')
        order=None
        if list(h.samples)!=samples:
            idx={s:i for i,s in enumerate(h.samples)}
            order=[idx[s] for s in samples]
        contig=chrom if chrom in h.contigs else chrom.removeprefix('chr')
        if contig in h.contigs:
            a,b=header.contigs[chrom].length,h.contigs[contig].length
            if a and b and a!=b:
                raise ValueError('Reference contig lengths disagree')
        renames={}
        for rec in h.records:
            cat=rec.type
            if cat not in ('INFO','FORMAT','FILTER','ALT'):
                continue
            name=rec['ID']
            existing={'INFO':header.info,'FORMAT':header.formats,'FILTER':header.filters,'ALT':header.alts}[cat]
            if name not in existing:
                header.add_line(str(rec).strip())
            elif cat in ('INFO','FORMAT') and (existing[name].type, str(existing[name].number)) != (rec['Type'],str(rec['Number'])):
                if name in ('AC','AF','SVLEN') and existing[name].type==rec['Type'] and {str(existing[name].number),str(rec['Number'])}<={'1','A','.'}:
                    continue
                if name in ('GT','END'):
                    raise ValueError(f'Incompatible reserved {cat}/{name}')
                target='UL_'+name
                if target not in existing:
                    header.add_line(str(rec).strip().replace('<ID='+name+',','<ID='+target+',',1))
                renames[cat,name]=target
        plans.append((contig,order,renames))
    return header,plans


def run_chrom(task):
    args,chrom,parts=task
    dest=args.out/chrom; dest.mkdir(parents=True,exist_ok=True)
    with ExitStack() as stack:
        lock=stack.enter_context((dest/'.lock').open('a'))
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        repaired=args.repaired/chrom
        try:
            fix,old,panel,mapping=short_receipts(args,chrom)
        except (OSError,ValueError,KeyError):
            if not args.recover_short:
                raise
            # The upstream batch runner takes a global lock and may bootstrap tools.
            scripts=Path(__file__).resolve().parent
            env=dict(os.environ,POPPED_FULL_ROOT=str(args.original),POPPED_FULL_JOBS='1',
                     POPPED_FULL_MEMORY_GB=str(max(4,26//args.jobs)))
            with (args.out/'.recovery.lock').open('a') as recovery:
                fcntl.flock(recovery,fcntl.LOCK_EX)
                subprocess.run(['bash',str(scripts/'run_popped_full_batch_workbench.sh'),'--chromosomes',chrom],
                               env=env,check=True)
                subprocess.run([sys.executable,str(scripts/'repair_popped_cluster_gts.py'),
                    '--batch-root',str(args.original),'--out',str(args.repaired),'--jobs','1','--chromosomes',chrom],check=True)
            fix,old,panel,mapping=short_receipts(args,chrom)
        key=dict(schema=VERSION,short=fix['artifacts'],map=old['artifacts']['old_to_new.jsonl.gz'],
                 long=[load(p/'complete.json')['hashes'] for p in parts],
                 long_source_map=digest(args.long/'partitions/source_records.tsv'),cut=CUT)
        if pilot.valid_manifest(dest,key):
            return load(dest/'stats.json')
        (dest/'done.json').unlink(missing_ok=True)
        print(f'Building {chrom}',flush=True)
        with (args.long/'partitions/source_records.tsv').open() as f:
            long_ids={r['audit_ID']:r for r in csv.DictReader(f,delimiter='\t')}
        removed=defaultdict(list)
        vcfs=[stack.enter_context(pysam.VariantFile(p/'retained.vcf.gz')) for p in parts]
        with pysam.VariantFile(panel) as vcf:
            hdr,plans=prepare_headers(vcf.header,vcfs,chrom)
        n=len(hdr.samples)
        if not any(contig in v.header.contigs for v,(contig,_,_) in zip(vcfs,plans)):
            raise ValueError(f'{chrom}: no compatible chromosome in lenient headers')
        streams=[short_stream(stack.enter_context(pysam.BGZFile(str(panel),'r')),chrom)]
        rules={}
        for p,vcf,(contig,order,names) in zip(parts,vcfs,plans):
            if contig not in vcf.header.contigs:
                continue
            streams.append(long_stream(vcf,contig,chrom,order))
            rules[p.name]=names
            if load(p/'complete.json')['summary']['collapsed_records']:
                with pysam.VariantFile(p/'removed.vcf.gz') as discarded:
                    for rec in discarded:
                        if rec.chrom==contig:
                            removed[p.name,rec.info['CollapseId']].append(rec.id)
        # Audit IDs are globally unique. Look up partition by source type/mode.
        by_type={}
        prepared=load(args.long/'partitions/manifest.json')
        for p in parts:
            cert=load(p/'complete.json')
            summary=cert['summary']
            for m in prepared['partitions'].values():
                if m['svtype']==summary['svtype'] and cert['recipe']['input_sha256']==prepared['hashes'][m['filename']]:
                    by_type[m['svtype'],m['mode']]=p.name
        group_rows=defaultdict(list)
        with gzip.open(mapping,'rt') as f:
            for row in map(json.loads,f):
                group_rows[(row['new_pos'],row['new_id'],row['new_ref'],row['new_alt'])].append(row)
        target=dest/f'{chrom}.unified.vcf.gz'
        tmp=dest/'panel.partial.vcf.gz'; cross=dest/'sv_mapping.jsonl.gz'
        metrics=defaultdict(Metrics); totals=Counter(); previous=0
        with pysam.BGZFile(str(tmp),'w') as out, gzip.open(cross,'wt') as maps:
            out.write(str(hdr).encode())
            for pos,source,ident,line in heapq.merge(*streams,key=lambda r:r[0]):
                if pos<previous:
                    raise ValueError('Input not coordinate sorted')
                previous=pos; totals[source+'_read']+=1
                if source=='popped' and totals['popped_read']%1000000==0:
                    print(f"{chrom}: read {totals['popped_read']:,} popped records",flush=True)
                fields=line.split(b'\t',8)
                ref,alt=fields[3].decode(),fields[4].decode()
                # Ordinary small alleles dominate I/O. Preserve entire original line.
                if source=='popped' and max(len(ref),len(alt))<50 and re.fullmatch('[ACGTNacgtn]+',alt):
                    out.write(line+b'\n'); totals['written']+=1; continue
                tags=info_dict(fields[7].decode())
                size=int(long_ids[ident]['length_bp']) if source=='lenient' else None
                prof=profile(pos,ref,alt,tags,size)
                if prof is None:
                    if source=='lenient':
                        raise ValueError('Unclassified lenient allele')
                    out.write(line+b'\n'); totals['written']+=1; continue
                kind,length,_,_,_=prof
                retained=source=='lenient' or length is None or length<=CUT
                pk=f'P.{chrom}.{ident}' if source=='popped' else f'L.{ident}'
                fmt,tail=fields[8].split(b'\t',1)
                old_id=fields[2].decode()
                members=group_rows.pop((pos,old_id,ref,alt),[]) if source=='popped' else []
                if members:
                    if len(members)!=members[0]['members']:
                        raise ValueError('Popped cluster cardinality mismatch')
                    pk='P.'+chrom+'.'+members[0]['representative_working_id']
                    if any((m['old']['net_bp']<=CUT)!=retained for m in members):
                        totals['boundary_crossing_clusters']+=1
                    for m in members:
                        maps.write(json.dumps(dict(source=source,panel_key=pk if retained else None,
                            disposition='included' if retained else 'excluded_by_tier',**m))+'\n')
                elif source=='lenient':
                    src=long_ids[ident]; part=by_type[src['SVTYPE'],src['mode']]
                    donors=removed.pop((part,tags.get('CollapseId')),[])
                    for member in [ident,*donors]:
                        maps.write(json.dumps(dict(source=source,panel_key=pk,old=long_ids[member],
                            new_id=old_id,new_pos=pos,new_ref=ref,new_alt=alt,members=1+len(donors),
                            status='collapsed' if member!=ident else 'merged_representative' if donors else 'unchanged'))+'\n')
                    names=rules[part]
                    tags={names.get(('INFO',k),k):v for k,v in tags.items()}
                    fmt=':'.join(names.get(('FORMAT',k),k) for k in fmt.decode().split(':')).encode()
                    if length<=CUT:
                        raise ValueError('Lenient input violates >10 kb policy')
                else:
                    maps.write(json.dumps(dict(source=source,panel_key=pk if retained else None,
                        old=dict(chrom=chrom,pos=pos,id=old_id,ref=ref,alt=alt),status='unchanged',
                        disposition='included' if retained else 'excluded_by_tier',members=1))+'\n')
                if not retained:
                    totals['popped_excluded_above_cut']+=1; continue
                ac,an=gt_counts(fmt,tail,n)
                classes=[kind,'ALL']
                if kind in ('simple_DEL','symbolic_DEL','simple_INS','symbolic_INS'):
                    classes.append(kind.rsplit('_',1)[1]+'_literal_and_symbolic')
                for label in classes:
                    for tier in (source,'combined'):
                        metrics[tier,label].add(pos,ref,alt,prof,ac,an,n)
                tags['UP_KEY']=pk
                fields[7]=';'.join(k if v is None else k+'='+v for k,v in tags.items()).encode()
                fields[8]=fmt+b'\t'+tail
                out.write(b'\t'.join(fields)+b'\n'); totals['written']+=1
        if group_rows or removed:
            raise ValueError('Unmatched saved cluster mapping')
        expected=load(repaired/'summary.json')['records']
        if totals['popped_read']!=expected or totals['written']!=expected+totals['lenient_read']-totals['popped_excluded_above_cut']:
            raise ValueError('Record conservation failed')
        pysam.tabix_index(str(tmp),preset='vcf',csi=True,force=True)
        with pysam.VariantFile(tmp,index_filename=str(tmp)+'.csi') as check:
            if list(check.header.samples)!=list(hdr.samples) or (totals['written'] and next(check.fetch(chrom),None) is None):
                raise ValueError('Indexed unified panel validation failed')
        tmp.replace(target); Path(str(tmp)+'.csi').replace(Path(str(target)+'.csi'))
        rows=[]; af=[]
        for (source,kind),m in sorted(metrics.items()):
            rows.append(dict(chrom=chrom,source=source,sv_class=kind,**m.row()))
            af.extend(dict(chrom=chrom,source=source,sv_class=kind,measure=measure,bin=b,alleles=v)
                      for (measure,b),v in sorted(m.bins.items()))
        result=dict(rows=rows,afs=af,counts=totals)
        json_write(dest/'stats.json',result)
        pilot.manifest(dest,key,[target,Path(str(target)+'.csi'),cross,dest/'stats.json'])
        return result


def aggregate(rows, group_cols, value_cols):
    totals=defaultdict(Counter)
    for r in rows:
        totals[tuple(r[k] for k in group_cols)].update({k:r.get(k,0) for k in value_cols})
    return rows+[dict(chrom='ALL',**dict(zip(group_cols,key)),**values) for key,values in sorted(totals.items())]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--original',type=Path,default=Path.home()/'popped_full_same_original_pos_p70')
    p.add_argument('--repaired',type=Path,default=Path.home()/'popped_full_same_original_pos_p70_gt_union')
    p.add_argument('--long',type=Path,default=Path.home()/'ashg_revised/v3_lenient_recollapse_same_start_90_70')
    p.add_argument('--out',type=Path,default=Path.home()/'unified_popped_lenient_p70')
    p.add_argument('--chromosomes',nargs='+',default=[f'chr{i}' for i in range(1,22)])
    p.add_argument('--jobs',type=int,choices=range(1,5),default=4)
    p.add_argument('--recover-short',action='store_true',help='Resume original merge/repair if its required files fail validation')
    args=p.parse_args()
    if len(set(args.chromosomes))!=len(args.chromosomes) or any(c not in [f'chr{i}' for i in range(1,22)] for c in args.chromosomes):
        p.error('Select distinct chromosomes from chr1 through chr21')
    os.umask(0o077); args.out.mkdir(parents=True,exist_ok=True)
    with (args.out/'.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        logging.basicConfig(filename=args.out/'run.log',level=logging.INFO)
        logging.info('PARAMETERS %s VERSION %s SEED %s',args,VERSION,pilot.SEED)
        parts=long_inputs(args.long)
        with ProcessPoolExecutor(max_workers=args.jobs) as pool:
            results=list(pool.map(run_chrom,[(args,c,parts) for c in args.chromosomes]))
        rows=[r for result in results for r in result['rows']]
        values=sorted({k for r in rows for k in r}-{'chrom','source','sv_class'})
        rows=aggregate(rows,['source','sv_class'],values)
        for r in rows:
            r['called_allele_fraction']=r['AN_sum']/r['diploid_slots'] if r['diploid_slots'] else 'NA'
        write_tsv(args.out/'sv_summary.tsv',rows)
        write_tsv(args.out/'allele_frequency_spectrum.tsv',aggregate([r for result in results for r in result['afs']],['source','sv_class','measure','bin'],['alleles']))
        write_tsv(args.out/'record_accounting.tsv',[dict(chrom=c,**r['counts']) for c,r in zip(args.chromosomes,results)])
        print('Overall SV counts, chr set:',','.join(args.chromosomes))
        print('class\tpositions\tallele_records\tunique_allele_keys\tunique_reference_span_bp\tunique_alt_payload_bp')
        for r in rows:
            if r['chrom']=='ALL' and r['source']=='combined':
                print('\t'.join(str(r[k]) for k in ('sv_class','positions','allele_records',
                    'unique_allele_keys','unique_reference_span_bp','unique_alt_payload_bp')))
        print('Completed:',args.out,flush=True)


if __name__=='__main__':
    main()
