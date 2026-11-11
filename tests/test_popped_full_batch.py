"""Native integration: original-position grouping and full-callset preservation."""
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


@unittest.skipUnless(os.environ.get('SV_PILOT_INTEGRATION')=='1','Requires native Linux tools')
class FullBatchTests(unittest.TestCase):
    def test_two_chromosomes_original_position_trim_only_merged_and_resume(self):
        bcftools = os.environ['BCFTOOLS_BIN']
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            inputs = root/'inputs'; inputs.mkdir()
            gt = [('1|0','0|0','0|0'),('0|0','0|1','0|0'),('0|0','0|0','1|0')]
            rows = [
                (100,'same_orig_a','AA','AA'+'C'*60,0),
                (100,'same_orig_b','A','A'+'C'*60,1),
                (199,'different_orig_a','AA','AA'+'G'*60,0),
                (200,'different_orig_b','A','A'+'G'*60,1),
                (300,'small','A','A'+'T'*49,0),
                (300,'sv50','A','A'+'T'*50,1),
                (400,'unmerged_padded','AA','AA'+'T'*100,0),
                (500,'replacement','G','T'*61,0),
                (600,'snv','A','C',0),
                (700,'del50','A'+'C'*50,'A',0),
                (700,'del54','A'+'C'*54,'A',1),
                (700,'ins50','A','A'+'C'*50,2),
                (800,'multiallelic','A','C,G',0),
            ]
            for chrom in ('chr1','chr2'):
                header = ('##fileformat=VCFv4.2\n'+f'##contig=<ID={chrom}>\n'+
                    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'+
                    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\tc\n')
                chosen = rows if chrom=='chr1' else [(100,'only_snv','A','G',0)]
                raw = inputs/f'{chrom}.vcf'
                raw.write_text(header+''.join('\t'.join(map(str,(chrom,pos,name,ref,alt,'.','PASS','.','GT',*gt[g])))+'\n'
                                             for pos,name,ref,alt,g in chosen))
                subprocess.run([bcftools,'view','-Ob','-o',str(inputs/f'{chrom}.bcf'),str(raw)],check=True,capture_output=True)
                subprocess.run([bcftools,'index',str(inputs/f'{chrom}.bcf')],check=True,capture_output=True)
            out = root/'run'
            cmd = [sys.executable,str(ROOT/'scripts/run_popped_full_batch.py'),'--root',str(out),
                '--input-template',str(inputs/'{chrom}.bcf'),'--chromosomes','chr1','chr2',
                '--jobs','2','--memory-gb','8','--bcftools',bcftools]
            def run():
                result = subprocess.run(cmd,text=True,capture_output=True)
                if result.returncode:
                    logs = '\n'.join(p.read_text()[-4000:] for p in out.rglob('*.log'))
                    self.fail(result.stdout+result.stderr+logs)
            run()
            def records(path):
                result = subprocess.check_output([bcftools,'view','--no-version','-H',str(path)],text=True)
                return {line.split('\t')[2]:line for line in result.splitlines()}
            original = records(inputs/'chr1.bcf')
            final = records(out/'chr1/final/chr1.popped.samepos.p70.vcf.gz')
            self.assertEqual(len(final),len(original)-2)
            with gzip.open(out/'chr1/final/old_to_new.jsonl.gz','rt') as handle:
                mapping = [json.loads(line) for line in handle]
            by_id = {r['old']['original_id']:r for r in mapping}
            native = {}
            for kind in ('INS','DEL'):
                native.update(records(out/'chr1/merge'/kind/'retained.vcf.gz'))
            self.assertEqual(by_id['same_orig_a']['representative_working_id'],by_id['same_orig_b']['representative_working_id'])
            self.assertNotEqual(by_id['different_orig_a']['representative_working_id'],by_id['different_orig_b']['representative_working_id'])
            self.assertEqual(by_id['del50']['representative_working_id'],by_id['del54']['representative_working_id'])
            for name in ['different_orig_a','different_orig_b','small','sv50','unmerged_padded','replacement','snv','ins50','multiallelic']:
                self.assertEqual(final[name],original[name],name)
            for row in mapping:
                if row['status']=='merged_representative':
                    fields = final[row['new_id']].split('\t')
                    self.assertEqual((int(fields[1]),fields[3],fields[4]),(row['new_pos'],row['new_ref'],row['new_alt']))
                    self.assertEqual(min(len(fields[3]),len(fields[4])),1)
                    self.assertNotIn('PM_ORIG_POS',fields[7])
                    self.assertEqual(fields[8:],native[row['representative_working_id']].split('\t')[8:])
            self.assertEqual(records(inputs/'chr2.bcf'),records(out/'chr2/final/chr2.popped.samepos.p70.vcf.gz'))
            marker = out/'chr1/final/done.json'
            stamp = marker.stat().st_mtime_ns
            run()
            self.assertEqual(stamp,marker.stat().st_mtime_ns)
            # Corrupt only final output: prepare/collapse must remain reusable.
            prep_stamp = (out/'chr1/prepare/done.json').stat().st_mtime_ns
            (out/'chr1/final/chr1.popped.samepos.p70.vcf.gz').write_bytes(b'broken')
            run()
            self.assertEqual(prep_stamp,(out/'chr1/prepare/done.json').stat().st_mtime_ns)
            self.assertEqual(final,records(out/'chr1/final/chr1.popped.samepos.p70.vcf.gz'))


if __name__=='__main__':
    unittest.main()
