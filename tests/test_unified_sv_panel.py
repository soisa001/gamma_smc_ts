"""Integration with completed synthetic short and long runs, including missing long GTs."""
import csv
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
import pysam

scripts=Path(__file__).resolve().parents[1]/'scripts'


def test_unified(tmp_path):
    work=tmp_path
    def run(*cmd):
        subprocess.run(list(map(str,cmd)),check=True,capture_output=True,text=True)
    header=('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=248956422>\n'
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n')
    sh=work/'short.vcf'
    sh.write_text(header+'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\n'+''.join(
        f'chr1\t{pos}\t{id_}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{a}\t{b}\n'
        for pos,id_,ref,alt,a,b in [(1,'snv','A','T','0|0','1|0'),(2,'indel','A','AT','0|0','1|0'),
        (20,'a','A','A'+'T'*60,'1|0','0|0'),(20,'b','A','A'+'T'*55,'0|0','1|0'),
        (30,'boundary','A','A'+'C'*10000,'0|1','0|0'),(50,'out','C','C'+'A'*10001,'1|0','0|0'),
        (60,'replacement','T'*61,'C'*62,'0|0','0|1')]))
    long=work/'long.vcf'
    long.write_text(header+'##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n'
        +'##INFO=<ID=SVTYPE,Number=1,Type=String,Description="type">\n'
        +'##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="length">\n'
        +'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tb\ta\n'
        +'chr1\t25\torigdel\tA\t<DEL>\t.\tPASS\tEND=12025;SVTYPE=DEL;SVLEN=-12000\tGT\t0/1\t./.\n')
    bcf=work/'v3_ultralong_lenient.bcf'
    for src,target in [(sh,work/'short.bcf'),(long,bcf)]:
        run('bcftools','view','-Ob','-o',target,src);run('bcftools','index',target)
    batch=work/'original'; repaired=work/'repaired'; lr=work/'longrun'; out=work/'unified'
    run(sys.executable,scripts/'run_popped_full_batch.py','--root',batch,'--input-template',work/'short.bcf','--chromosomes','chr1','--jobs',1)
    run(sys.executable,scripts/'repair_popped_cluster_gts.py','--batch-root',batch,'--out',repaired,'--chromosomes','chr1','--jobs',1)
    run(sys.executable,scripts/'compare_stringent_recollapse.py','--bcf',bcf,'--out',lr,'--same-start','--workers',1)
    command=(sys.executable,scripts/'build_unified_sv_panel.py','--original',batch,'--repaired',repaired,'--long',lr,'--out',out,'--chromosomes','chr1','--jobs',1)
    run(*command)
    target=out/'chr1/chr1.unified.vcf.gz'
    with pysam.VariantFile(target) as vcf:
        recs=list(vcf.fetch('chr1'))
        assert [r.pos for r in recs]==[1,2,20,25,30,60]
        assert recs[3].samples['a']['GT']==(None,None) and recs[3].samples['b']['GT']==(0,1)
        assert recs[2].samples['b']['GT']==(1,0)
    with gzip.open(out/'chr1/sv_mapping.jsonl.gz','rt') as f:
        rows=[json.loads(line) for line in f]
    assert sum(r['panel_key'] is None for r in rows)==1
    assert len({r['panel_key'] for r in rows if r['panel_key']})==4
    with (out/'sv_summary.tsv').open() as f:
        totals=next(r for r in csv.DictReader(f,delimiter='\t') if (r['chrom'],r['source'],r['sv_class'])==('ALL','combined','ALL'))
    assert (int(totals['positions']),int(totals['allele_records']))==(4,4)
    assert int(totals['unique_reference_span_bp'])==12000
    assert int(totals['AN_sum'])==14
    with (out/'allele_frequency_spectrum.tsv').open() as f:
        spectrum=list(csv.DictReader(f,delimiter='\t'))
    assert [(r['bin'],r['alleles']) for r in spectrum if
            (r['chrom'],r['source'],r['sv_class'],r['measure'])==('ALL','lenient','symbolic_DEL','called_AF')]==[('0.5<=AF<1','1')]
    before=target.stat().st_mtime_ns;run(*command);assert target.stat().st_mtime_ns==before
