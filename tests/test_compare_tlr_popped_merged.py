import csv
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import compare_tlr_popped_merged as report


class MetricTests(unittest.TestCase):
    def test_union_excludes_anchor_clips_and_counts_overlap_once(self):
        rows = [dict(chrom='chr4',pos=100,ref='A'+'C'*50,alt='A',kind='DEL',simple=True),
                dict(chrom='chr4',pos=100,ref='A'+'C'*54,alt='A',kind='DEL',simple=True)]
        self.assertEqual(report.deleted_union(rows,('chr4',1,1000)),54)
        self.assertEqual(report.deleted_union(rows,('chr4',110,120)),11)
        self.assertEqual(report.deleted_union(rows,('chr4',100,100)),0)


@unittest.skipUnless(os.environ.get('SV_PILOT_INTEGRATION')=='1','Requires native tools')
class IntegrationTests(unittest.TestCase):
    def test_counts_mapping_plot_identity_and_boundary_overlap(self):
        bcftools = os.environ['BCFTOOLS_BIN']
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            header = ('##fileformat=VCFv4.2\n##contig=<ID=chr4>\n'
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\n')
            records = [(38659980,'boundary_del','A'+'C'*50,'A',0),
                       (38754377,'ins60','A','A'+'T'*60,0),
                       (38754377,'ins65','A','A'+'T'*65,1),
                       (38800000,'replacement','G','T'*61,0),
                       (38807500,'del50','A'+'C'*50,'A',0),
                       (38807500,'del54','A'+'C'*54,'A',1),
                       (38810000,'snp','G','C',0)]
            raw = root/'source.vcf'
            raw.write_text(header+''.join('\t'.join(map(str,('chr4',pos,name,ref,alt,'.','PASS','.','GT',
                *(['1|0','0|0'] if sample==0 else ['0|0','0|1']))))+'\n' for pos,name,ref,alt,sample in records))
            source = root/'source.bcf'
            subprocess.run([bcftools,'view','-Ob','-o',str(source),str(raw)],check=True,capture_output=True)
            subprocess.run([bcftools,'index',str(source)],check=True,capture_output=True)
            batch = root/'batch'
            cmd = [sys.executable,str(ROOT/'scripts/run_popped_full_batch.py'),'--root',str(batch),
                '--input-template',str(source),'--chromosomes','chr4','--jobs','1','--bcftools',bcftools]
            done = subprocess.run(cmd,capture_output=True,text=True)
            if done.returncode:
                self.fail(done.stderr+'\n'.join(p.read_text()[-3000:] for p in batch.rglob('*.log')))
            plot = root/'plot'; plot.mkdir()
            with gzip.open(plot/'popped_alleles.tsv.gz','wt') as handle:
                writer = csv.writer(handle,delimiter='\t')
                writer.writerow(['sv_event_id','chrom','original_pos','original_id','ref','alt','svtype','sv_length_bp'])
                writer.writerow(['plot1','chr4',38754377,'ins60','A','A'+'T'*60,'INS',60])
                writer.writerow(['plot2','chr4',38807500,'del50','A'+'C'*50,'A','DEL',50])
            (plot/'custom_TLR10_TLR1_TLR6_sv_alleles_001.tsv').write_text(
                'sv_event_id\tallele_label\tmarker_color\nplot1\tSV1\tblue\nplot2\tSV2\tgreen\n')
            cmd = [sys.executable,str(ROOT/'scripts/compare_tlr_popped_merged.py'),'--batch-root',str(batch),
                   '--plot-dir',str(plot),'--bcftools',bcftools]
            done = subprocess.run(cmd,capture_output=True,text=True)
            self.assertEqual(done.returncode,0,done.stdout+done.stderr)
            output = batch/'tlr_before_after'
            with open(output/'counts.tsv') as handle:
                counts = {(r['scope'],r['type']):r for r in csv.DictReader(handle,delimiter='\t')}
            self.assertEqual((counts['plot_window','TOTAL']['before'],counts['plot_window','TOTAL']['after']),('6','4'))
            self.assertEqual((counts['plot_window','INS']['before'],counts['plot_window','INS']['after']),('3','2'))
            self.assertEqual((counts['TLR_interval','TOTAL']['before'],counts['TLR_interval','TOTAL']['after']),('3','2'))
            summary = json.loads((output/'summary.json').read_text())
            self.assertEqual(summary['plot']['plotted_popped_alleles'],2)
            self.assertEqual(summary['plot']['distinct_mapped_merged_alleles'],2)
            self.assertNotIn('ins60',done.stdout)
            with open(output/'plotted_allele_changes.tsv') as handle:
                plotted = list(csv.DictReader(handle,delimiter='\t'))
            self.assertEqual([r['cluster_members'] for r in plotted],['2','2'])


if __name__=='__main__':
    unittest.main()
