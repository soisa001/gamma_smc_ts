"""Export checked summary tables and figures; large simulation data stay on D."""
import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import shutil


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def main(out,destination):
    source=out/'analysis'
    audit=json.loads((source/'audit.json').read_text())
    if audit['status']!='passed':raise ValueError('Export requires a completed independent audit')
    manifest=json.loads((source/'artifact_manifest.json').read_text())
    excluded={'heldout_predictions.csv.gz','fold_region_scores.npz','target70_predictions.csv'}
    destination.mkdir(parents=True,exist_ok=True)
    copied=[]
    for name,spec in manifest.items():
        if name in excluded:continue
        if digest(source/name)!=spec['sha256'] or (source/name).stat().st_size!=spec['bytes']:
            raise ValueError(f'Corrupt summary input: {name}')
        shutil.copyfile(source/name,destination/name);copied.append(name)
    shutil.copyfile(source/'audit.json',destination/'audit.json');copied.append('audit.json')
    shutil.copyfile(source/'artifact_manifest.json',destination/'full_analysis_manifest.json');copied.append('full_analysis_manifest.json')
    focus={'af_r1','joint_g50_r1_T50000','mass_g50_r1_T50000','excess_g50_r1_T50000','tuned_all','tuned_joint'}
    compact=destination/'primary_predictions.csv.gz';count=0
    with gzip.open(source/'heldout_predictions.csv.gz','rt',newline='') as stream, compact.open('wb') as raw:
        reader=csv.DictReader(stream)
        with gzip.GzipFile(fileobj=raw,mode='wb',mtime=0) as zipped, io.TextIOWrapper(zipped,newline='') as text:
            writer=csv.DictWriter(text,fieldnames=reader.fieldnames);writer.writeheader()
            for row in reader:
                if row['method'] in focus:writer.writerow(row);count+=1
    copied.append(compact.name)
    for name in ('calibration_thresholds.csv','region_scores.csv','provenance.json'):
        path=out/'af_ablation'/name
        if path.exists():
            target='af_'+name;shutil.copyfile(path,destination/target);copied.append(target)
    record=dict(source=str(source),source_manifest_sha256=digest(source/'artifact_manifest.json'),
                compact_prediction_rows=count,large_intermediates_remain_at=str(out),
                files={name:dict(bytes=(destination/name).stat().st_size,sha256=digest(destination/name)) for name in sorted(copied)})
    (destination/'artifact_manifest.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(dict(files=len(copied)+1,bytes=sum(v['bytes'] for v in record['files'].values()),destination=str(destination)),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--destination',type=Path,required=True)
    args=parser.parse_args();main(args.out,args.destination)
