"""Small phased genotype consolidation regression cases."""
from pathlib import Path
import sys
import tempfile
import unittest
import json
import os
import subprocess
from types import SimpleNamespace

import pysam
import truvari.collapse as native

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from phased_sv_union import snapshot, union_calls, apply_union


def records():
    header = pysam.VariantHeader()
    header.contigs.add('chr4',length=190214555)
    header.formats.add('GT',1,'String','Genotype')
    header.formats.add('PS',1,'Integer','Phase set')
    header.add_sample('synthetic_donor')
    kept = header.new_record(contig='chr4',start=38754376,alleles=('A','A'+'T'*652),id='representative')
    removed = header.new_record(contig='chr4',start=38754376,alleles=('A','A'+'T'*574),id='donor')
    kept.samples[0]['GT'] = (0,0); kept.samples[0].phased = True
    removed.samples[0]['GT'] = (1,0); removed.samples[0].phased = True
    removed.samples[0]['PS'] = 38750000
    return kept,removed


class TransferTests(unittest.TestCase):
    def test_reference_representative_receives_donor_gt(self):
        kept,removed = records()
        original = kept.copy()
        result,_ = native.gt_aware_consolidate(kept,[SimpleNamespace(comp=removed)])
        # Pin the actual native failure, then verify the proposed union fix.
        self.assertEqual(result.samples[0]['GT'],(0,0))
        apply_union(result,[original,removed])
        self.assertEqual(result.samples[0]['GT'],(1,0))
        self.assertTrue(result.samples[0].phased)
        self.assertEqual(result.samples[0]['PS'],38750000)

    def test_intrasample_union_not_dosage_addition(self):
        def call(gt,phased=True,ps=10):
            return dict(gt=gt,phased=phased,ps=ps)
        self.assertEqual(union_calls([call((1,0)),call((1,0))])['gt'],(1,0))
        self.assertEqual(union_calls([call((1,0)),call((0,1))])['gt'],(1,1))
        self.assertEqual(union_calls([call((0,0)),call((0,1))])['gt'],(0,1))
        self.assertEqual(union_calls([call((None,0)),call((0,1))])['gt'],(None,1))
        self.assertEqual(union_calls([call((1,0),ps=10),call((0,1),ps=20)])['status'],'phase_set_conflict')
        value = union_calls([call((0,0)),call((0,1),phased=False)])
        self.assertEqual(value['gt'],(0,1)); self.assertFalse(value['phased'])
        self.assertEqual(union_calls([call((0,1),False),call((0,1),False)])['status'],'unphased_union_ambiguous')

    def test_gt_off_native_does_not_union_two_het_haplotypes(self):
        kept,removed = records()
        kept.samples[0]['GT'] = (1,0); kept.samples[0].phased=True
        removed.samples[0]['GT'] = (0,1); removed.samples[0].phased=True
        originals = [kept.copy(),removed.copy()]
        # With a representative already ALT, the gt-off consolidation branch
        # leaves it alone; switching the flag alone is insufficient.
        result,_ = native.collapse_into_entry(kept,[SimpleNamespace(comp=removed)],False)
        self.assertEqual(result.samples[0]['GT'],(1,0))
        apply_union(result,originals)
        self.assertEqual(result.samples[0]['GT'],(1,1))

    def test_gt_all_rejects_shared_carriers_before_consolidation(self):
        kept,removed=records()
        kept.samples[0]['GT']=(0,1);kept.samples[0].phased=True
        cluster=native.CollapsedCalls(kept,'test')
        cluster.genotype_mask=cluster.make_genotype_mask(kept,'all')
        self.assertTrue(cluster.gt_conflict(removed,'all'))


@unittest.skipUnless(os.environ.get('SV_PILOT_INTEGRATION')=='1','Requires native tools')
class AuditTests(unittest.TestCase):
    def test_original_native_final_audit_and_intrasample_trial(self):
        root_code=Path(__file__).resolve().parents[1]
        bcftools=os.environ['BCFTOOLS_BIN']
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);raw=root/'source.vcf';source=root/'source.bcf'
            header=('##fileformat=VCFv4.2\n##contig=<ID=chr4>\n'
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
                '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="PS">\n'
                '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\n')
            rows=[('4_92311',652,'1|0:100','0|0:100'),('4_92312',574,'0|0:100','1|0:100'),
                  ('c520',520,'0|0:100','0|1:100')]
            raw.write_text(header+''.join('\t'.join(['chr4','38754377',name,'A','A'+'T'*length,
                '.','PASS','.','GT:PS',a,b])+'\n' for name,length,a,b in rows)
                +'chr4\t38754380\tsnv\tA\tG\t.\tPASS\t.\tGT:PS\t1|0:100\t0|1:100\n')
            subprocess.run([bcftools,'view','-Ob','-o',str(source),str(raw)],check=True,capture_output=True)
            subprocess.run([bcftools,'index',str(source)],check=True,capture_output=True)
            batch=root/'batch'
            command=[sys.executable,str(root_code/'scripts/run_popped_full_batch.py'),'--root',str(batch),
                '--input-template',str(source),'--chromosomes','chr4','--jobs','1','--bcftools',bcftools]
            run=subprocess.run(command,capture_output=True,text=True)
            self.assertEqual(run.returncode,0,run.stdout+run.stderr)
            from run_chr1_sv_merge import common
            final=batch/'chr4/final/chr4.popped.samepos.p70.vcf.gz'
            before_hash=common.digest(final)
            command=[sys.executable,str(root_code/'scripts/audit_tlr_gt_transfer.py'),'--batch-root',str(batch)]
            run=subprocess.run(command,capture_output=True,text=True)
            self.assertEqual(run.returncode,0,run.stdout+run.stderr)
            output=batch/'tlr_gt_transfer_audit_38754377'
            result=json.loads((output/'summary.json').read_text())
            self.assertEqual(result['original_alleles'],3)
            self.assertEqual(result['existing_clusters'],2)
            self.assertEqual(result['existing_run_audit']['samples_losing_known_ALT_dosage'],1)
            self.assertEqual(result['existing_run_audit']['native_final_gt_phase_ps_mismatches'],0)
            self.assertEqual(result['existing_run_audit']['native_replay_gt_phase_ps_mismatches'],0)
            self.assertEqual(result['intrasample_union_trial']['records'],1)
            with pysam.VariantFile(str(output/'intrasample.gt_union.vcf.gz')) as vcf:
                records=list(vcf)
                self.assertEqual(len(records),1)
                self.assertEqual(records[0].samples['b']['GT'],(1,1))
                self.assertTrue(records[0].samples['b'].phased)
                self.assertEqual(records[0].samples['b']['PS'],100)
                self.assertEqual(records[0].info['AC'],(3,))
                self.assertEqual(records[0].info['AN'],4)
            self.assertEqual(common.digest(final),before_hash)
            plot_dir=batch/'tlr_merge_carrier_examples'
            plot_dir.mkdir()
            (plot_dir/'selected_haplotypes.tsv').write_text('label\tsample_id\thap\nS02.h1\tb\t1\n')
            isolated=root/'existing_only'
            run=subprocess.run(command+['--existing-only','--out',str(isolated)],capture_output=True,text=True)
            self.assertEqual(run.returncode,0,run.stdout+run.stderr)
            checked=json.loads((isolated/'summary.json').read_text())['existing_cluster_verification']
            self.assertEqual(checked['gt_phase_ps_mismatches'],0)
            self.assertEqual(checked['focused_h1'],'restored')
            self.assertFalse((isolated/'intrasample.gt_union.vcf.gz').exists())
            repaired_root=root/'repaired'
            command=[sys.executable,str(root_code/'scripts/repair_popped_cluster_gts.py'),
                     '--batch-root',str(batch),'--out',str(repaired_root),'--chromosomes','chr4','--jobs','1']
            run=subprocess.run(command,capture_output=True,text=True)
            self.assertEqual(run.returncode,0,run.stdout+run.stderr)
            target=repaired_root/'chr4/chr4.popped.samepos.p70.gt_union.vcf.gz'
            with pysam.VariantFile(str(target)) as vcf, pysam.VariantFile(str(isolated/'existing_clusters.gt_union.vcf.gz')) as trial:
                fixed={r.id:r for r in vcf.fetch('chr4')}
                for r in trial:
                    for sample in r.samples:
                        self.assertEqual(snapshot(fixed[r.id].samples[sample]),snapshot(r.samples[sample]))
                    if r.id=='4_92311':
                        self.assertEqual(fixed[r.id].info['AC'],r.info['AC'])
                with pysam.VariantFile(str(final)) as old:
                    for r in old:
                        if r.id!='4_92311':
                            self.assertEqual(str(fixed[r.id]),str(r))
            modified=target.stat().st_mtime_ns
            subprocess.run(command,check=True,capture_output=True)
            self.assertEqual(target.stat().st_mtime_ns,modified)
            from summarize_popped_gt_union import summarize, write_table, scan
            reports=root/'summaries'
            reports.mkdir()
            row=summarize(('chr4',batch,repaired_root,reports))
            self.assertEqual((row['INS_before'],row['INS_after'],row['GT_changed']),(3,2,1))
            write_table(reports,[row])
            self.assertIn('33.33',(reports/'chr1_21_summary.tsv').read_text())
            complex_fixture=root/'replacements.vcf'
            complex_fixture.write_text(header+''.join(
                f'chr4\t{i}\tx{i}\t{ref}\t{alt}\t.\tPASS\t.\tGT:PS\t0|0:1\t1|0:1\n'
                for i,(ref,alt) in enumerate([('A','T'*60),('C'*80,'G'*75),
                                               ('A','A'+'C'*60),('A','A'+'N'*60)],1)))
            bgz=pysam.tabix_index(str(complex_fixture),preset='vcf')
            counts=scan(Path(bgz),'chr4',4)
            self.assertEqual((counts['INS'],counts['complex_net50'],counts['complex_span50_net_lt50']),(1,1,1))


if __name__=='__main__':
    unittest.main()
