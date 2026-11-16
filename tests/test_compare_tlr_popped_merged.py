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
import plot_tlr_before_after as haplotypes


class MetricTests(unittest.TestCase):
    def test_complex_and_ambiguous_alleles_are_separate_and_cannot_enter_merge(self):
        def row(ref,alt):
            desc,_ = report.rep.describe(ref,alt)
            return dict(chrom='chr4',pos=38800000,id='x',ref=ref,alt=alt,
                        kind=desc['type'],net_bp=desc['net_bp'],
                        simple=not desc['complex'] and not desc['ambiguous'])
        pure = row('AC','A'+'T'*50+'C')
        complex_allele = row('AC','A'+'T'*51)
        ambiguous = row('A','A'+'N'*50)
        rows = [pure,complex_allele,ambiguous]
        self.assertEqual([report.allele_class(r) for r in rows],
                         ['simple','complex_replacement','ambiguous_sequence'])
        counts = report.summarize_classes(rows,rows)
        self.assertEqual([r['before'] for r in counts if r['scope']=='TLR_interval' and r['type']=='INS'],[1,1,1])
        mapping = dict(old=dict(chrom='chr4',original_pos=38800000,original_id='x',
                               original_ref=complex_allele['ref'],original_alt=complex_allele['alt']))
        with self.assertRaisesRegex(ValueError,'complex or ambiguous'):
            report.check_mapping([complex_allele],[complex_allele],[mapping])

    def test_phase_slots_never_inferred_from_unphased_heterozygotes(self):
        state = haplotypes.gt_state
        self.assertEqual([state((0,1),False,1,h) for h in (1,2)],['unknown','unknown'])
        self.assertEqual([state((1,1),False,1,h) for h in (1,2)],['alt','alt'])
        self.assertEqual([state((None,1),True,1,h) for h in (1,2)],['unknown','alt'])
        self.assertEqual([state((1,2),True,2,h) for h in (1,2)],['other','alt'])
        self.assertEqual(state((1,),True,1,1),'unknown')

    def test_selection_prioritizes_introgressed_carriers_and_distinct_alleles(self):
        haps = [('a',1),('b',1),('c',1),('d',1)]
        calls = {'x':{h:dict(state='alt' if h in haps[:2] else 'other') for h in haps},
                 'y':{h:dict(state='alt' if h==haps[2] else 'other') for h in haps}}
        tracts = {h:[(0,10)] for h in haps}
        selected,_ = haplotypes.select_haps(haps,tracts,calls,[dict(old=['x','y'])],3)
        self.assertEqual(selected,[haps[0],haps[2],haps[1]])

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
            with gzip.open(plot/'popped_selected_lanes.tsv.gz','wt') as handle:
                handle.write('sample_id\thap\tchrom\tstart\tend\n'
                             'a\t1\tchr4\t38780000\t38850000\n'
                             'b\t2\tchr4\t38800000\t38890001\n')
            cmd = [sys.executable,str(ROOT/'scripts/compare_tlr_popped_merged.py'),'--batch-root',str(batch),
                   '--plot-dir',str(plot),'--bcftools',bcftools,'--plot-haplotypes','--max-haplotypes','4']
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
            with open(output/'counts_by_class.tsv') as handle:
                stratified = {(r['scope'],r['allele_class'],r['type']):r for r in csv.DictReader(handle,delimiter='\t')}
            self.assertEqual(stratified['TLR_interval','complex_replacement','INS']['removed'],'0')
            selected = haplotypes.table(output/'selected_haplotypes.tsv')
            selected_keys = {(r['sample_id'],r['hap']) for r in selected}
            self.assertTrue({('a','1'),('b','2')}.issubset(selected_keys))
            # Native --gt all can change a donor's phase slot. Its mate must be
            # visible instead of making the ALT appear to vanish from the plot.
            audit = haplotypes.table(output/'haplotype_placement_audit.tsv')
            for row in audit:
                if row['comparison']=='changed':
                    self.assertIn((row['sample_id'],str(3-int(row['hap']))),selected_keys)
            states = haplotypes.table(output/'haplotype_genotypes.tsv.gz')
            # Source deletion alleles differ across the selected haps. Read the
            # actual output GTs independently to verify their plotted placement.
            self.assertEqual({r['id'] for r in states if r['stage']=='Unmerged' and r['state']=='alt'},
                             {'replacement','del50','del54'})
            import pysam
            final = batch/'chr4/final/chr4.popped.samepos.p70.vcf.gz'
            with pysam.VariantFile(str(final)) as vcf:
                actual = {rec.id:rec.copy() for rec in vcf.fetch('chr4',38780000,38890001)}
            for r in states:
                if r['stage']!='Unmerged':
                    sample = actual[r['id']].samples[r['sample_id']]
                    expected = haplotypes.gt_state(sample['GT'],sample.phased,1,int(r['hap']))
                    self.assertEqual(r['state'],expected)
            png = output/'tlr_haplotypes_before_after.png'
            pdf = output/'tlr_haplotypes_before_after.pdf'
            self.assertTrue(png.read_bytes().startswith(b'\x89PNG'))
            self.assertTrue(pdf.read_bytes().startswith(b'%PDF'))
            self.assertNotIn(b'/Subtype /Image',pdf.read_bytes())
            self.assertEqual(summary['haplotypes']['aggregate']['merged_clusters_shown'],1)
            # Re-run after corruption: report-only plotting repairs the artifact
            # without changing the completed merge or its input manifest.
            final_hash = report.pilot.common.digest(final)
            png.write_bytes(b'corrupt')
            rerun = subprocess.run(cmd,capture_output=True,text=True)
            self.assertEqual(rerun.returncode,0,rerun.stderr)
            self.assertEqual(report.pilot.common.digest(final),final_hash)
            self.assertTrue(png.read_bytes().startswith(b'\x89PNG'))
            if os.environ.get('TLR_SYNTHETIC_QA_DIR'):
                import shutil
                dest = Path(os.environ['TLR_SYNTHETIC_QA_DIR']); dest.mkdir(parents=True,exist_ok=True)
                for path in (png,pdf,output/'haplotype_plot_summary.json'):
                    shutil.copyfile(path,dest/path.name)


if __name__=='__main__':
    unittest.main()
