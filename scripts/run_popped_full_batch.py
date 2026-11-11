#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = ["truvari==5.5.0", "google-crc32c==1.8.0"]
# ///
"""Full popped-panel VCFs: native simple-SV collapse within original POS/type.

Only actual merged representatives are replaced with trimmed native outputs.
Every other source record is copied unchanged. Controlled data stay on the VM.
"""
import argparse
from collections import Counter, defaultdict, deque
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import ExitStack
import csv
import gzip
import importlib.metadata
import json
import logging
import os
from pathlib import Path
import subprocess
import sys
import time

import pysam

import run_chr1_sv_merge as pilot
import run_popped_simple_sv as simple

VERSION = 1
TEMPLATE = simple.representation.SOURCE.replace('.chr1.', '.{chrom}.')
POLICY = dict(version=VERSION, eligibility='trim_then_pure_INS_DEL_net_ge50',
    grouping='original_CHROM_POS_and_type', pctseq=.7, pctsize=.7, refdist=500,
    keep='common', gt='all', chain=False, hap=False, left_align=False,
    output='trim_actual_merged_representatives_only; other_source_records_unchanged',
    truvari='5.5.0', genotype_consolidation='native')


def source_paths(path):
    indexes = [p for ext in ('.csi','.tbi') if (p:=Path(str(path)+ext)).is_file()]
    if not path.is_file() or not indexes:
        raise ValueError('Source must exist with its index')
    return [path, *indexes]


def source_identity(source):
    return [pilot.common.identity(p) for p in source_paths(source)]


def choose_source(options, root, chrom, uri):
    if options.input_template:
        return Path(options.input_template.format(chrom=chrom)).resolve(), False
    staged = root/'inputs'/uri.rsplit('/',1)[1]
    if options.source_template != TEMPLATE:
        return staged, True
    # Reuse immutable known popped inputs, including the earlier chr1/TLR trials.
    # Never select bubble or imputed files or modify the existing trial roots.
    candidates = []
    marker = root/'prepare/done.json'
    if marker.is_file():
        try:
            candidates.append(Path(json.loads(marker.read_text())['key']['settings']['source'][0]['path']))
        except (KeyError, ValueError, IndexError):
            pass
    candidates.append(staged)
    home = Path.home()
    if chrom == 'chr1':
        metadata = home/'chr1_sv_representation_comparison/comparison_metadata.json'
        if metadata.is_file():
            try:
                candidates.append(Path(json.loads(metadata.read_text())['inputs']['popped']))
            except (KeyError, ValueError, TypeError):
                pass
    folders = [f'{chrom}_popped_simple_merge', f'{chrom}_confirmed_popped_representation', f'{chrom}_popped_representation']
    if chrom == 'chr4':
        folders.append('tlr_popped_sv_audit')
    candidates.extend(home/folder/'inputs'/staged.name for folder in folders)
    candidates.append(home/'panel_popped_vcf'/staged.name)
    for path in candidates:
        if path.name == staged.name and path.is_file() and Path(str(path)+'.csi').is_file():
            return path.resolve(), path.resolve()==staged.resolve()
    return staged, True


def prepare(args, source):
    directory = args.root/'prepare'; directory.mkdir(exist_ok=True)
    snapshot = source_identity(source)
    settings = dict(policy=POLICY, source=snapshot, chrom=args.chrom)
    old = pilot.valid_manifest(directory)
    if old and old['key']['settings'] == settings:
        return old
    (directory/'done.json').unlink(missing_ok=True)
    samples = pilot.command(args, [args.bcftools,'query','-l',source], 'samples', True)
    if not samples.strip():
        raise ValueError('Full sample genotypes are required')
    (directory/'samples.txt').write_text(samples)
    counts, candidates, identifier = Counter(), Counter(), 0
    with ExitStack() as stack:
        handles = {kind:stack.enter_context(pysam.BGZFile(str(directory/f'{kind}.vcf.gz'), 'w')) for kind in simple.KINDS}
        mapping = stack.enter_context(gzip.open(directory/'source_alleles.jsonl.gz','wt'))
        lines = stack.enter_context(pilot.stream(args, [args.bcftools,'view','--no-version','-r',args.chrom,
            '-i','ILEN>=50 || ILEN<=-50','-Ov',source], 'prepare-candidates'))
        for line in lines:
            if line.startswith('##INFO=<ID=PM_ORIG_POS,'):
                raise ValueError('Reserved PM_ORIG_POS already exists in source header')
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    for handle in handles.values():
                        handle.write(b'##INFO=<ID=PM_ORIG_POS,Number=1,Type=Integer,Description="Original POS used for merge grouping">\n')
                for handle in handles.values():
                    handle.write(line.encode())
                continue
            fields = line.rstrip('\n').split('\t',8)
            if ',' in fields[4]:
                counts['multiallelic_records_unchanged'] += 1
                continue
            kind, size, pos, ref, alt, eligible = simple.classify(int(fields[1]), fields[3], fields[4])
            simple.bump(counts,kind,size)
            if not eligible or size < 50:
                continue
            identifier += 1
            working = f'PM{identifier:09d}'
            meta = dict(working_id=working, chrom=fields[0], original_pos=int(fields[1]),
                original_id=fields[2], original_ref=fields[3], original_alt=fields[4],
                trimmed_pos=pos, trimmed_ref=ref, trimmed_alt=alt, kind=kind, net_bp=size)
            mapping.write(json.dumps(meta)+'\n')
            info = pilot.info_parse(fields[7])
            if ('SVTYPE' in info and info['SVTYPE'] != kind) or ('SVLEN' in info and abs(int(info['SVLEN'])) != size):
                raise ValueError('Inconsistent simple-SV type/length annotations in source')
            info['PM_ORIG_POS'] = fields[1]
            if 'END' in info:
                info['END'] = str(pos+len(ref)-1)
            fields[1:5] = [str(pos),working,ref,alt]
            fields[7] = pilot.info_string(info)
            handles[kind].write(('\t'.join(fields)+'\n').encode())
            candidates[kind] += 1
            if identifier % 10000 == 0:
                logging.info('PREPARE_PROGRESS eligible_sv_records=%d',identifier)
    if source_identity(source) != snapshot:
        raise ValueError('Source changed during prepare')
    report = dict(settings=settings, before_sv50=counts, candidates=candidates,
                  samples=len(samples.splitlines()))
    pilot.common.json_write(directory/'summary.json',report)
    return pilot.manifest(directory, report, [directory/n for n in
        ('INS.vcf.gz','DEL.vcf.gz','source_alleles.jsonl.gz','samples.txt','summary.json')])


def merge(args, prepared):
    manifests = {}
    for kind in simple.KINDS:
        directory = args.root/'merge'/kind; directory.mkdir(parents=True,exist_ok=True)
        key = dict(policy=POLICY, prepared=prepared, kind=kind)
        cached = pilot.valid_manifest(directory,key)
        if cached:
            manifests[kind] = cached
            continue
        (directory/'done.json').unlink(missing_ok=True)
        retained, removed = [directory/f'{name}.partial.vcf.gz' for name in ('retained','collapsed')]
        cmd = [sys.executable,Path(__file__).with_name('collapse_original_start.py'),
            '-i',args.root/'prepare'/f'{kind}.vcf.gz','-o',retained,'-c',removed,
            '-p','.7','-P','.7','-r','500','-O','0','-s','50','-S','-1','--keep','common','--gt','all']
        pilot.command(args,cmd,kind+'-collapse')
        files, numbers = [], []
        for partial, name in ((retained,'retained'),(removed,'collapsed')):
            final = directory/f'{name}.vcf.gz'
            numbers.append(simple.finish(args,partial,final,kind+'-'+name,sort=True))
            files.extend(simple.artifacts([final]))
        if sum(numbers) != prepared['key']['candidates'].get(kind,0):
            raise ValueError('Retained + collapsed count differs from input')
        mapping = directory/'mapping.tsv'
        Path(str(retained)+'.mapping.tsv').replace(mapping)
        with open(mapping) as handle:
            rows = list(csv.DictReader(handle,delimiter='\t'))
        if len(rows) != sum(numbers) or len({r['working_id'] for r in rows}) != len(rows):
            raise ValueError('Incomplete/duplicate membership mapping')
        pilot.common.json_write(directory/'summary.json',dict(retained=numbers[0],collapsed=numbers[1]))
        manifests[kind] = pilot.manifest(directory,key,files+[mapping,directory/'summary.json'])
    return manifests


def original_key(meta):
    return (meta['chrom'],str(meta['original_pos']),meta['original_id'],meta['original_ref'],meta['original_alt'])


def finalize(args, source, prepared, merged):
    directory = args.root/'final'; directory.mkdir(exist_ok=True)
    key = dict(policy=POLICY,prepared=prepared,merged=merged,source=source_identity(source))
    if key['source'] != prepared['key']['settings']['source']:
        raise ValueError('Source changed since prepare')
    if pilot.valid_manifest(directory,key):
        return json.loads((directory/'summary.json').read_text())
    (directory/'done.json').unlink(missing_ok=True)
    with gzip.open(args.root/'prepare/source_alleles.jsonl.gz','rt') as handle:
        meta = {r['working_id']:r for r in map(json.loads,handle)}
    membership = {}
    for kind in simple.KINDS:
        with open(args.root/'merge'/kind/'mapping.tsv') as handle:
            for row in csv.DictReader(handle,delimiter='\t'):
                row['members'] = int(row['members'])
                membership[row['working_id']] = row
    if set(meta) != set(membership):
        raise ValueError('Source-to-cluster mapping incomplete')
    for ident, row in membership.items():
        a, b = meta[ident], meta[row['representative']]
        if (a['chrom'],a['original_pos'],a['kind']) != (b['chrom'],b['original_pos'],b['kind']):
            raise ValueError('Mapping crossed original position/type boundary')
    wanted = {r['representative'] for r in membership.values() if r['members'] > 1}
    offsets, extra_headers = {}, {}
    # BGZF virtual offsets let final assembly load one merged GT record at a
    # time rather than retaining all cohort genotypes in Python memory.
    for kind in simple.KINDS:
        path = args.root/'merge'/kind/'retained.vcf.gz'
        with pysam.BGZFile(str(path),'r') as handle:
            while True:
                offset = handle.tell(); line = handle.readline()
                if not line:
                    break
                text = line.decode()
                if text.startswith('##INFO=') and 'ID=PM_ORIG_POS,' not in text:
                    extra_headers[text.split(',',1)[0]] = text+'\n'
                if text.startswith('#'):
                    continue
                ident = text.split('\t',3)[2]
                if ident in wanted:
                    offsets[ident] = (kind,offset)
    if set(offsets) != wanted:
        raise ValueError('Missing merged representative')
    lookup = defaultdict(deque)
    for ident,m in meta.items():
        lookup[original_key(m)].append(ident)
    positions = {m['original_pos'] for m in meta.values()}
    partial = directory/f'{args.chrom}.popped.samepos.p70.partial.vcf.gz'
    final = directory/f'{args.chrom}.popped.samepos.p70.vcf.gz'
    counts = Counter()
    sort_cmd = [args.bcftools,'sort','-m',f'{max(1,args.memory_gb//2)}G',
                '-T',directory/'sort-tmp','-Oz','-o',partial,'-']
    logging.info('FINAL_STREAM %s',list(map(str,sort_cmd)))
    with ExitStack() as stack:
        retained = {k:stack.enter_context(pysam.BGZFile(str(args.root/'merge'/k/'retained.vcf.gz'),'r')) for k in simple.KINDS}
        log = stack.enter_context(open(args.root/'logs/final-sort.log','a'))
        lines = stack.enter_context(pilot.stream(args,[args.bcftools,'view','--no-version','-r',args.chrom,'-Ov',source],'final-source'))
        proc = subprocess.Popen(list(map(str,sort_cmd)),stdin=subprocess.PIPE,stdout=log,stderr=log)
        try:
            header_keys = set()
            for line in lines:
                if line.startswith('##INFO='):
                    header_keys.add(line.split(',',1)[0])
                if line.startswith('#CHROM'):
                    for name, value in extra_headers.items():
                        if name not in header_keys:
                            proc.stdin.write(value.encode())
                if line.startswith('#'):
                    proc.stdin.write(line.encode()); continue
                counts['before_records'] += 1
                fields = line.rstrip('\n').split('\t',8)
                ident = None
                if int(fields[1]) in positions:
                    queue = lookup.get(tuple(fields[:5]))
                    if queue:
                        ident = queue.popleft()
                if ident is not None:
                    row = membership[ident]
                    if row['representative'] != ident:
                        counts['removed_records'] += 1
                        counts['removed_'+meta[ident]['kind']] += 1
                        continue
                    if row['members'] > 1:
                        kind,offset = offsets[ident]
                        retained[kind].seek(offset)
                        fields = retained[kind].readline().decode().split('\t',8)
                        if fields[2] != ident:
                            raise ValueError('Invalid representative BGZF offset')
                        fields[2] = meta[ident]['original_id']
                        info = pilot.info_parse(fields[7]); info.pop('PM_ORIG_POS',None)
                        # Native cluster counters restart per type; namespace the final annotation.
                        if 'CollapseId' in info:
                            info['CollapseId'] = kind+'_'+info['CollapseId']
                        fields[7] = pilot.info_string(info)
                        line = '\t'.join(fields).rstrip('\n')+'\n'
                        counts['merged_representatives'] += 1
                    else:
                        counts['unchanged_records'] += 1
                else:
                    counts['unchanged_records'] += 1
                proc.stdin.write(line.encode())
                counts['after_records'] += 1
                if counts['before_records'] % 100000 == 0:
                    logging.info('FINAL_PROGRESS %s',dict(counts))
            proc.stdin.close()
            if proc.wait():
                raise RuntimeError('Final sort failed; inspect final-sort.log')
        finally:
            if proc.poll() is None:
                proc.terminate(); proc.wait()
    if any(lookup.values()):
        raise ValueError('A prepared source allele was not encountered in final assembly')
    if source_identity(source) != key['source']:
        raise ValueError('Source changed during final assembly')
    total = simple.finish(args,partial,final,'final',sort=False)
    if total != counts['after_records'] or counts['before_records']-total != counts['removed_records']:
        raise ValueError('Final count conservation failed')
    samples = pilot.command(args,[args.bcftools,'query','-l',final],'final-samples',True)
    if samples != (args.root/'prepare/samples.txt').read_text():
        raise ValueError('Final sample order differs from source')
    mapping_path = directory/'old_to_new.jsonl.gz'
    with gzip.open(mapping_path,'wt') as handle:
        for ident,m in meta.items():
            member = membership[ident]; rep = meta[member['representative']]
            change = member['members'] > 1
            handle.write(json.dumps(dict(old=m, representative_working_id=rep['working_id'],
                new_id=rep['original_id'], new_pos=rep['trimmed_pos'] if change else rep['original_pos'],
                new_ref=rep['trimmed_ref'] if change else rep['original_ref'],
                new_alt=rep['trimmed_alt'] if change else rep['original_alt'],
                members=member['members'], status='collapsed' if ident!=rep['working_id'] else 'merged_representative' if change else 'unchanged'))+'\n')
    summary = dict(chrom=args.chrom,**counts,policy=POLICY,
                   eligible_sv_before=prepared['key']['candidates'],before_sv50=prepared['key']['before_sv50'])
    summary['after_sv50'] = {k:prepared['key']['before_sv50'].get(k,0)-counts['removed_'+k] for k in simple.KINDS}
    pilot.common.json_write(directory/'summary.json',summary)
    pilot.manifest(directory,key,simple.artifacts([final])+[mapping_path,directory/'summary.json'])
    return summary


def worker(options, chrom):
    root = options.root/chrom
    (root/'logs').mkdir(parents=True,exist_ok=True)
    uri = options.source_template.format(chrom=chrom)
    args = pilot.parse_args(['--root',str(root),'--chrom',chrom,'--threads','1',
        '--memory-gb',str(max(4,options.memory_gb//options.jobs)),
        '--source-uri',uri,'--confirm-source-uri',uri,'--bcftools',options.bcftools,'--gcloud',options.gcloud])
    args.billing_project = options.billing_project
    logging.basicConfig(filename=root/'logs/workflow.log',level=logging.INFO,
                        format='%(asctime)s %(process)d %(levelname)s %(message)s',force=True)
    import fcntl
    with open(root/'.full_merge.lock','a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        provenance = dict(parameters=json.loads(json.dumps(vars(options),default=str)),chrom=chrom,source_uri=uri,seed=pilot.SEED,policy=POLICY,
            python=sys.version,script_sha256=pilot.common.digest(__file__),
            bcftools=pilot.command(args,[args.bcftools,'--version'],'version',True),
            per_chromosome=dict(threads=1,memory_gb=args.memory_gb))
        pilot.common.json_write(root/'logs'/f'run-{time.time_ns()}.json',provenance)
        logging.info('PROVENANCE %s',json.dumps(provenance,default=str))
        print(f'Starting {chrom}; phase={options.phase}',flush=True)
        try:
            source, own_staging = choose_source(options,root,chrom,uri)
            logging.info('SOURCE local_path=%s owned_staging=%s',source,own_staging)
            if own_staging and options.phase in ('stage','prepare','all'):
                pilot.stage(args)
            if options.phase == 'stage':
                return dict(chrom=chrom,stage='complete')
            prepared = prepare(args,source) if options.phase in ('prepare','all') else pilot.require_manifest(root/'prepare')
            if options.phase == 'prepare':
                return dict(chrom=chrom,prepare='complete')
            merged = merge(args,prepared) if options.phase in ('merge','all') else {
                k:pilot.require_manifest(root/'merge'/k) for k in simple.KINDS}
            if options.phase == 'merge':
                return dict(chrom=chrom,merge='complete')
            result = finalize(args,source,prepared,merged)
            print(f"Completed {chrom}: removed {result.get('removed_records',0):,} records",flush=True)
            return result
        except Exception:
            logging.exception('CHROMOSOME_FAILED')
            raise RuntimeError(f'{chrom} failed; inspect {root}/logs/workflow.log') from None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--chromosomes',nargs='+',default=[f'chr{i}' for i in range(1,22)])
    parser.add_argument('--jobs',type=int,default=4)
    parser.add_argument('--memory-gb',type=int,default=26)
    parser.add_argument('--source-template',default=TEMPLATE)
    parser.add_argument('--input-template',help='Optional local indexed popped file path containing {chrom}')
    parser.add_argument('--phase',choices=['stage','prepare','merge','finalize','all'],default='all')
    parser.add_argument('--bcftools',default=os.environ.get('BCFTOOLS_BIN','bcftools'))
    parser.add_argument('--gcloud',default='gcloud')
    parser.add_argument('--billing-project',default=None)
    options = parser.parse_args()
    if not 1<=options.jobs<=4 or options.memory_gb<4*options.jobs:
        parser.error('Use 1–4 jobs and at least 4 GB budget per job')
    if len(set(options.chromosomes))!=len(options.chromosomes) or not set(options.chromosomes)<={f'chr{i}' for i in range(1,22)}:
        parser.error('Select unique chromosomes from chr1–chr21')
    if importlib.metadata.version('truvari')!='5.5.0':
        raise ValueError('Requires Truvari 5.5.0')
    options.root = options.root.resolve(); options.root.mkdir(parents=True,exist_ok=True)
    import fcntl
    with open(options.root/'.batch.lock','a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        rows, failed = [], []
        with ProcessPoolExecutor(max_workers=options.jobs) as pool:
            futures = {pool.submit(worker,options,ch):ch for ch in options.chromosomes}
            for future in as_completed(futures):
                try:
                    rows.append(future.result())
                except Exception as error:
                    failed.append(futures[future]); print(str(error),file=sys.stderr,flush=True)
                pilot.common.json_write(options.root/'batch_status.json',dict(completed=rows,failed=failed))
        if failed:
            raise SystemExit('Failed chromosomes: '+', '.join(failed)+'; completed phases are reusable')
        if options.phase in ('all','finalize'):
            comparison = []
            for row in sorted(rows,key=lambda x:int(x['chrom'][3:])):
                for kind in (*simple.KINDS,'TOTAL'):
                    kinds = simple.KINDS if kind=='TOTAL' else (kind,)
                    before = sum(row['before_sv50'].get(k,0) for k in kinds)
                    after = sum(row['after_sv50'].get(k,0) for k in kinds)
                    comparison.append(dict(chrom=row['chrom'],type=kind,before_sv50=before,after_sv50=after,
                        removed=before-after,percent_removed=100*(before-after)/before if before else 0))
            pilot.tsv(options.root/'comparison.tsv',comparison,list(comparison[0]))
        print(f'Completed batch; outputs stay in {options.root}',flush=True)


if __name__=='__main__':
    main()
