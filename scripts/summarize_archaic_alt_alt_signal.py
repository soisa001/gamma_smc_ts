"""Present the adopted ALT/ALT signal from audited existing-cohort results."""
import argparse
import csv
import hashlib
import json
from pathlib import Path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path):
    with path.open(newline='') as stream:return list(csv.DictReader(stream))


def write_csv(path,rows):
    with path.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)


def summarize(source,config_path,out):
    cfg=json.loads(config_path.read_text())
    assert cfg['statistic']=='frac_recent_T' and cfg['source_family']=='mass'
    assert cfg['all_pair_gate_quantile']==0 and cfg['consecutive_bins']==1
    assert not cfg['ref_ref_subtraction'] and not cfg['alt_ref_penalty']
    audit=json.loads((source/'audit.json').read_text());assert audit['status']=='passed' and audit['no_later_onset']
    manifest=json.loads((source/'artifact_manifest.json').read_text())
    for name in ('metrics.csv','target70_metrics.csv','provenance.json'):
        assert digest(source/name)==manifest[name]['sha256'],name
        assert (source/name).stat().st_size==manifest[name]['bytes'],name
    provenance=json.loads((source/'provenance.json').read_text())
    assert provenance['config']['selection_onset_years']==cfg['selection_onset_years']==cfg['pulse_years']
    prefix=f"{cfg['source_family']}_g0_r1_T"
    rows=[r for r in read_csv(source/'metrics.csv') if r['method'].startswith(prefix)]
    target=[r for r in read_csv(source/'target70_metrics.csv') if r['method']==cfg['primary_method']]
    assert len(rows)==72 and len(target)==4
    for r in rows+target:
        assert int(r['onset_years'])==cfg['selection_onset_years']
        assert (int(r['selected_n']),int(r['neutral_n']))==(100,1000)
        assert abs(float(r['power'])-int(r['selected_called'])/100)<1e-14
        assert abs(float(r['neutral_call_fraction'])-int(r['neutral_called'])/1000)<1e-14
    out.mkdir(parents=True,exist_ok=True)
    write_csv(out/'metrics.csv',rows);write_csv(out/'target70_metrics.csv',target)
    names={'grid10kb':'10 kb grid','archaic_sites':'Every archaic site'}
    lines=['# Archaic ALT/ALT signal: existing results','',
           'Score: `frac_recent_T_among_ALT_ALT * n_ALT_ALT / n_total_pairs`.',
           'The saved method identifiers retain the historical `mass_g0_r1_` prefix.',
           'Only immediate selection at the 50-kya pulse is included.','',
           '## Regional p <= 0.05, T = 50 kya','',
           '| Positions | TMRCA source | Selected called / 100 | Neutral called / 1,000 |',
           '|---|---|---:|---:|']
    primary=[r for r in rows if r['method']==cfg['primary_method'] and float(r['alpha'])==cfg['primary_region_alpha']]
    assert len(primary)==4
    for r in sorted(primary,key=lambda r:(r['scheme']!='grid10kb',r['source']!='truth')):
        lines.append(f"| {names[r['scheme']]} | {r['source'].title()} | {r['selected_called']} ({float(r['power']):.0%}) | {r['neutral_called']} ({float(r['neutral_call_fraction']):.1%}) |")
    lines.extend(['','## Thresholds trained toward 70% power','',
                  '| Positions | TMRCA source | Selected called / 100 | Neutral called / 1,000 |',
                  '|---|---|---:|---:|'])
    for r in sorted(target,key=lambda r:(r['scheme']!='grid10kb',r['source']!='truth')):
        lines.append(f"| {names[r['scheme']]} | {r['source'].title()} | {r['selected_called']} ({float(r['power']):.0%}) | {r['neutral_called']} ({float(r['neutral_call_fraction']):.1%}) |")
    lines.extend(['','Thresholds use the selected training regions; reported outcomes use held-out regions.',
                  'These are exploratory existing-cohort estimates, not new-seed validation.',''])
    (out/'tables.md').write_text('\n'.join(lines))
    record=dict(display_name=cfg['display_name'],config=cfg,config_sha256=digest(config_path),
                source_directory=str(source),source_audit_sha256=digest(source/'audit.json'),
                source_files={n:dict(sha256=digest(source/n)) for n in ('metrics.csv','target70_metrics.csv','provenance.json')},
                script_sha256=digest(Path(__file__)),selected_metric_rows=len(rows),target70_rows=len(target),
                no_simulation_or_decoding_rerun=True)
    (out/'provenance.json').write_text(json.dumps(record,indent=2)+'\n')
    files=['metrics.csv','target70_metrics.csv','tables.md','provenance.json']
    (out/'artifact_manifest.json').write_text(json.dumps({n:dict(bytes=(out/n).stat().st_size,sha256=digest(out/n)) for n in files},indent=2)+'\n')
    print(json.dumps(dict(status='verified',rows=len(rows),target70_rows=len(target),out=str(out)),indent=2))


if __name__=='__main__':
    repo=Path(__file__).resolve().parents[1]
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,default=Path('/mnt/d/phase2simselection/sim/eas_allele_class_ablation'))
    p.add_argument('--config',type=Path,default=repo/'configs/archaic_alt_alt_signal.json')
    p.add_argument('--out',type=Path,default=repo/'docs/results/archaic_alt_alt_signal')
    a=p.parse_args();summarize(a.source,a.config,a.out)
