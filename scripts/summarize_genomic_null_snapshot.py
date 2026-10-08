"""Add a concise PDF report and allele-frequency plots to a verified CDF snapshot."""
import os
os.environ['MPLBACKEND']='Agg'
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[k]='1'
import argparse
import hashlib
import json
from pathlib import Path
import zipfile
import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42,'font.size':16})
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from pypdf import PdfReader, PdfWriter
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import letter, landscape
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak

POPS=('AFR','AMR','EAS','EUR','MID','SAS')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--snapshot',type=Path,required=True)
    p.add_argument('--campaign',type=Path,required=True)
    p.add_argument('--report-date',default='2026-10-07',help='ISO date for the report heading')
    args=p.parse_args()
    out=args.snapshot
    manifest=json.loads((out/'manifest.json').read_text())
    assert manifest['schema']=='population-null-ecdf/v2','Use the empirical inclusive-tail CDF snapshot'
    # Verify the existing snapshot before extending it; do not read new scores.
    for name,spec in manifest['outputs'].items():
        assert sha(out/name)==spec['sha256'],name
    source=pd.DataFrame(manifest['source_receipts'])
    scores=pd.read_csv(out/'null_scores_10k_50k.csv',float_precision='round_trip')
    cuts=pd.read_csv(out/'cutoffs_all_methods.csv',float_precision='round_trip')
    coverage=pd.read_csv(out/'history_coverage.csv')
    status=json.loads((args.campaign/'run_status.json').read_text())
    source[['population','task_id','history_id','sample_af','onset_af','attempts']].to_csv(out/'focal_allele_frequencies.csv',index=False)
    afrows=[]
    fig,axes=plt.subplots(2,3,figsize=(17,10),layout='constrained')
    for ax,pop in zip(axes.flat,POPS):
        group=source[source.population==pop]
        af=group.sample_af.to_numpy()
        q95=float(np.quantile(af,.95,method='inverted_cdf'))
        ax.hist(af,bins=np.linspace(0,1,26),weights=np.ones(len(af))/len(af),color='#267bc1',edgecolor='white')
        ax.axvline(q95,color='#e77735',lw=2,ls='--',label='Empirical 95th percentile')
        ax.set(title=f'{pop} | n = {len(af)}\n95th percentile: {q95:.1%}',xlim=(0,1),xlabel='Present-day focal ALT frequency',ylabel='Fraction of retained nulls')
        ax.xaxis.set_major_formatter(PercentFormatter(1))
        ax.yaxis.set_major_formatter(PercentFormatter(1))
        ax.grid(axis='y',alpha=.2)
        afrows.append(dict(population=pop,n=len(af),median_af=float(np.median(af)),q95_af=q95,
            mean_af=float(np.mean(af)),fixed_n=int(np.sum(af==1)),median_attempts=float(np.median(group.attempts)),
            maximum_attempts=int(group.attempts.max())))
    fig.legend(*axes.flat[0].get_legend_handles_labels(),loc='outside upper center',ncol=1,frameon=False)
    fig.savefig(out/'allele_frequency_distributions.png',dpi=160)
    fig.savefig(out/'allele_frequency_distributions.pdf')
    plt.close(fig)
    afstats=pd.DataFrame(afrows)
    afstats.to_csv(out/'allele_frequency_summary.csv',index=False)
    comparisons=[]
    for (pop,method,t),g in scores.groupby(['population','method','cutoff_years']):
        pivot=g.pivot(index='task_id',columns='source',values='score').dropna()
        delta=pivot.decoded-pivot.truth
        comparisons.append(dict(population=pop,method=method,cutoff_years=t,paired_n=len(pivot),
            median_decoded_minus_truth=float(delta.median()),mean_absolute_difference=float(delta.abs().mean())))
    pd.DataFrame(comparisons).to_csv(out/'truth_decoded_comparison.csv',index=False)
    styles=getSampleStyleSheet()
    styles.add(ParagraphStyle('TitleLarge',fontName='Helvetica-Bold',fontSize=24,leading=28,spaceAfter=14))
    styles.add(ParagraphStyle('BodyLarge',fontName='Helvetica',fontSize=11.5,leading=14,spaceAfter=7))
    styles.add(ParagraphStyle('SubLarge',fontName='Helvetica-Bold',fontSize=15,leading=17,spaceBefore=6,spaceAfter=6))
    story=[]
    def text(value,kind='BodyLarge'):
        story.append(Paragraph(value,styles[kind]))
    def table(rows,widths):
        t=Table(rows,colWidths=widths,repeatRows=1,hAlign='LEFT')
        t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e4eef7')),
            ('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),('FONTNAME',(0,1),(-1,-1),'Helvetica'),
            ('FONTSIZE',(0,0),(-1,-1),12),('LEADING',(0,0),(-1,-1),15),('BOTTOMPADDING',(0,0),(-1,-1),6),
            ('TOPPADDING',(0,0),(-1,-1),6),('LINEBELOW',(0,0),(-1,0),.7,colors.HexColor('#607d95')),
            ('LINEBELOW',(0,1),(-1,-1),.3,colors.HexColor('#ccd5dd')),('ALIGN',(1,0),(-1,-1),'RIGHT')]))
        story.append(t)
    n=manifest['regions']
    text('Genomic null calibration: first 100 per population' if manifest.get('first_per_population')==100
         else 'Genomic null calibration: interim results','TitleLarge')
    text(f'{args.report_date} | Frozen snapshot: <b>{n:,} / 6,000</b> completed simulations ({n/6000:.1%}). '
         f'Target: 1,000 per population. Live runner at report generation: {status["state"]}; '
         f'{status["workers"]} workers; {len(status["failed"])} reported failures.')
    text('Decoded TMRCA: all sampled haplotype pairs','SubLarge')
    text('Four requested numbers per population. Entries are percentages of pairs with TMRCA below the specified age; '
         'they are descriptive empirical 95th/99th percentiles across retained null focal sites.')
    primary=pd.read_csv(out/'four_cutoffs_decoded_all_pairs.csv').set_index('population')
    rows=[['Population','Nulls','10 kya: top 5%','10 kya: top 1%','50 kya: top 5%','50 kya: top 1%']]
    for pop in POPS:
        r=primary.loc[pop]
        rows.append([pop,str(int(r.null_n))]+[f'{r[k]:.2f}%' for k in ('T10k_top5_pct','T10k_top1_pct','T50k_top5_pct','T50k_top1_pct')])
    table(rows,[86,47,136,136,136,136])
    text('Empirical p-values: ties included, no +1 adjustment','SubLarge')
    text('<b>p = number of null scores greater than or equal to the observation / n.</b> '
         'With n = 100, p &lt;= 0.05 requires strictly exceeding the 95th-percentile boundary; '
         'p &lt;= 0.01 requires strictly exceeding the 99th-percentile boundary. Equality does not reject.')
    text('The p-value grid has 0.01 steps. The 99th percentile is the second-largest null score, '
         'so the 1% tail is still imprecisely estimated. An observation above every null has empirical p = 0, '
         'which does not establish zero population tail probability.')
    history_note=('Every represented history currently has one completed region; ten are planned per history. '
                  if (coverage.maximum_regions_per_history==1).all() else 'History replication counts are supplied in history_coverage.csv. ')
    text('Pooled CDFs are shown without demographic confidence bands. '+history_note+
         'This report freezes the prespecified cohort; later completed replicates are excluded.')
    story.append(PageBreak())
    text('True TMRCA and allele-pair availability','TitleLarge')
    text('For the decoded all-pairs statistic, reject at p &lt;= 0.05 only when the target score is <b>strictly greater</b> '
         'than the top-5% value on the previous page. The true-TMRCA boundaries below use the same rule.')
    text('True TMRCA: all sampled haplotype pairs','SubLarge')
    truth=pd.read_csv(out/'four_cutoffs_truth_all_pairs.csv').set_index('population')
    rows=[['Population','Nulls','10 kya: top 5%','10 kya: top 1%','50 kya: top 5%','50 kya: top 1%']]
    for pop in POPS:
        r=truth.loc[pop]
        rows.append([pop,str(int(r.null_n))]+[f'{r[k]:.2f}%' for k in ('T10k_top5_pct','T10k_top1_pct','T50k_top5_pct','T50k_top1_pct')])
    table(rows,[86,47,136,136,136,136])
    text('ALT/ALT availability and saturation at 50 kya','SubLarge')
    rows=[['Population','Available / n','Missing','p at 100%: true','p at 100%: decoded']]
    for pop in POPS:
        g=cuts[(cuts.population==pop)&(cuts.method=='alt_alt')&(cuts.tmrca_years==50000)&(cuts.alpha==.05)].set_index('source')
        r=g.loc['truth']
        rows.append([pop,f'{int(r.available_n)} / {int(r.null_n)}',str(int(r.missing_n)),
                     f'{g.loc["truth","p_at_score_one"]:.2f}',f'{g.loc["decoded","p_at_score_one"]:.2f}'])
    table(rows,[110,120,85,162,200])
    text('At a target score of 100%, p is the fraction of null scores also equal to 100%. '
         'This table counts ties in the numerator. If every null equals the target, p = 1.')
    story.append(PageBreak())
    text('Raw ALT/ALT fractions: empirical cutoffs','TitleLarge')
    text('Percentages of carrier pairs; no carrier-mass weighting. The denominator remains all 100 null simulations. '
         'Undefined null carrier-pair scores are retained below finite support; undefined targets are no-calls. '
         'Available-only percentiles and tied counts are separately labeled in the full cutoff CSV.')
    for mode in ('decoded','truth'):
        text(f'{mode.capitalize()} TMRCA: ALT/ALT pairs','SubLarge')
        wide=pd.read_csv(out/f'four_cutoffs_{mode}_alt_alt.csv').set_index('population')
        rows=[['Population','Nulls','10 kya: top 5%','10 kya: top 1%','50 kya: top 5%','50 kya: top 1%']]
        for pop in POPS:
            r=wide.loc[pop]
            rows.append([pop,str(int(r.null_n))]+[f'{r[k]:.2f}%' for k in ('T10k_top5_pct','T10k_top1_pct','T50k_top5_pct','T50k_top1_pct')])
        table(rows,[86,47,136,136,136,136])
    story.append(PageBreak())
    text('Model, allele frequencies and limitations','TitleLarge')
    text('Neutral archaic-derived alleles, segregating immediately after 2% introgression at 50 kya and observed in the present sample '
         '(fixation retained). Neanderthal Ne = 3,600; split = 500 kya; 25 years/generation. '
         '200 sampled diploids; 10,000 fixed sampled haplotype pairs; no rescaling.')
    text('Random autosomal regions use the local deCODE GRCh38 map in simulation. Gamma-SMC uses the cropped region mean '
         'recombination rate, accurate exp10, and fixed backward alignment. Mutation rate = 1.29e-8. '
         'Assembly gaps and missing map coverage are masked; the HMMix strict mask is not used.')
    text('True versus decoded and ALT/ALT results','SubLarge')
    text('The following pages include allele-frequency distributions, separate true/decoded CDFs, and all-pairs/ALT-ALT comparisons. '
         'ALT/ALT scores use raw carrier-pair fractions. Null no-calls remain below finite score support in the CDF and percentile '
         'denominator. See the CSVs for available counts, tied fractions and p at each cutoff.')
    text('This campaign contains calibration nulls only. It does not supply new selected-treatment power or independent neutral-target FPR estimates.')
    saturated=cuts[(cuts.method=='alt_alt')&(cuts.tmrca_years==50000)&(cuts.alpha==.05)]
    if (~saturated.empirical_p_attainable).all():
        text('<b>ALT/ALT at 50 kya is saturated:</b> every population has a 100% upper-tail boundary for both true and decoded TMRCA. '
             'Even a target score of 100% cannot reject at p &lt;= 0.05 with these nulls.')
    text('Focal ALT frequencies among retained nulls','SubLarge')
    rows=[['Population','Median AF','95th-percentile AF','Fixed alleles','Median attempts']]
    for r in afstats.itertuples():
        rows.append([r.population,f'{r.median_af:.1%}',f'{r.q95_af:.1%}',str(r.fixed_n),f'{r.median_attempts:g}'])
    table(rows,[97,120,175,140,145])
    SimpleDocTemplate(str(out/'summary.pdf'),pagesize=landscape(letter),leftMargin=36,rightMargin=36,topMargin=30,bottomMargin=30).build(story)
    writer=PdfWriter()
    for name in ('summary.pdf','allele_frequency_distributions.pdf','null_cdf_all.pdf'):
        writer.append(out/name)
    with (out/'genomic_null_report.pdf').open('wb') as stream:
        writer.write(stream)
    writer.close()
    notes=dict(snapshot_regions=n,counts=manifest['counts'],live_run_status=status,
        history_coverage=coverage.to_dict('records'),simulations_changed=False,
        p_definition=manifest['p_definition'],
        contents='Inclusive empirical tail tables, tied fractions, AF distributions, 32 true/decoded all-pairs/ALT-ALT CDF pages',
        source_snapshot_manifest_sha256=sha(out/'manifest.json'),script_sha256=sha(Path(__file__)),
        report_pages=len(PdfReader(out/'genomic_null_report.pdf').pages))
    files=sorted(f for f in out.iterdir() if f.suffix in ('.pdf','.png','.csv','.md','.json') and f.name!='report_manifest.json')
    notes['outputs']={f.name:sha(f) for f in files}
    (out/'report_manifest.json').write_text(json.dumps(notes,indent=2)+'\n')
    with zipfile.ZipFile(out/'genomic_null_report_bundle.zip','w',compression=zipfile.ZIP_DEFLATED) as z:
        for f in files+[out/'report_manifest.json']:
            z.write(f,f.name)
    print(json.dumps({k:v for k,v in notes.items() if k not in ('outputs','live_run_status','history_coverage')}),flush=True)


if __name__=='__main__':
    main()
