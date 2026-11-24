import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import pysam


def test_counts_observed_carriers_before_taking_cohort_median(tmp_path):
    script = Path(__file__).resolve().parents[1]/'scripts/median_sv_per_person.py'
    for i in range(1,22):
        chrom = f'chr{i}'
        folder = tmp_path/chrom
        folder.mkdir()
        path = folder/f'{chrom}.unified.vcf.gz'
        with pysam.BGZFile(str(path),'w') as f:
            header = (f'##fileformat=VCFv4.2\n##contig=<ID={chrom},length=100000>\n'
                      '##INFO=<ID=UP_KEY,Number=1,Type=String,Description="key">\n'
                      '##FORMAT=<ID=GT,Number=1,Type=String,Description="genotype">\n'
                      '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="depth">\n'
                      '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\tc\n')
            f.write(header.encode())
            # Rotating carriers exercise summation by person before taking the median.
            calls = ['1|1:9','0|1:4','./.:.']
            k = i%3
            calls = calls[k:]+calls[:k]
            f.write(('\t'.join([chrom,'1','.','A','A'+'C'*50,'.','PASS','UP_KEY=P.ins','GT:DP',*calls])+'\n').encode())
            f.write(f'{chrom}\t2\t.\tA\t<DEL>\t.\tPASS\tUP_KEY=L.del\tGT\t1/.\t.\t0/0\n'.encode())
            f.write(f'{chrom}\t3\t.\tA\tT\t.\tPASS\t.\tGT\t1/1\t1/1\t1/1\n'.encode())
        pysam.tabix_index(str(path),preset='vcf',csi=True)
        stats = folder/'stats.json'
        stats.write_text(json.dumps({'rows':[
            dict(source='popped',sv_class='simple_INS',allele_records=1,AC_sum=3,AN_sum=4),
            dict(source='lenient',sv_class='symbolic_DEL',allele_records=1,AC_sum=1,AN_sum=3)]}))
        (folder/'done.json').write_text(json.dumps({'artifacts':{
            p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (path,stats)}}))
    command=[sys.executable,str(script),'--jobs','2','--panel',str(tmp_path)]
    subprocess.run(command,check=True,capture_output=True)
    with (tmp_path/'per_person_summary/chr1_21_medians.tsv').open() as f:
        rows={r['sv_class']:r for r in csv.DictReader(f,delimiter='\t') if r['source']=='combined'}
    assert float(rows['simple_INS']['median_SVs_carried'])==14
    assert float(rows['simple_INS']['median_ALT_copies'])==21
    assert float(rows['symbolic_DEL']['median_SVs_carried'])==0
    assert int(rows['symbolic_DEL']['people_with_observed_ALT'])==1
    assert float(rows['ALL']['median_SVs_carried'])==14
    assert int(rows['ALL']['panel_records'])==42
