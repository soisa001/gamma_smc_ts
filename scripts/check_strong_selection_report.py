#!/usr/bin/env python3
"""Audit strong-selection cohorts, pointwise calls, benchmark counts and figures."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import argparse
import hashlib
import json
from pathlib import Path
import re
import numpy as np
from PIL import Image


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=Path(__file__).resolve().parents[1]/'docs/results/strong_selection_20261110')
    args=parser.parse_args();out=args.out
    digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    inventory=json.loads((out/'output_hashes.json').read_text())
    assert all(digest(out/name)==value for name,value in inventory.items())
    meta=json.loads((out/'manifest.json').read_text())
    accounting=json.loads((out/'proposal_accounting.json').read_text())['strata']
    summary=json.loads((out/'summary.json').read_text())['results']
    benchmark=json.loads((out/'focal_site_detection_benchmark.json').read_text())['results']
    n_cases=len(meta['cases']);assert len(accounting)==6
    assert len(summary)==6*n_cases and len(benchmark)==12*n_cases
    with np.load(out/'neutral.npz') as data:neutral=data['cdfs'].copy()
    with np.load(out/'neutral_af.npz') as data:
        assert np.all((data['trajectories'][:,-1,1]>0)&(data['trajectories'][:,-1,1]<1))
    for row in accounting:
        assert row['accepted']+row['lost']+row['fixed']==row['new_completed']+row['reused']
        assert row['lost_before_selection_onset']+row['lost_after_selection_onset']==row['lost']
        assert row['retained']==min(row['accepted'],row['target'])
        assert row['detection_estimable']==(row['retained']>0)
    names=['True TMRCA','Posterior median','Posterior-weighted mass']
    for entry in meta['cases']:
        key=entry['case']['key']
        with np.load(out/entry['cache']) as src:data={k:src[k].copy() for k in src.files}
        n=len(data['ids']);assert len(set(data['ids']))==n
        assert np.all((data['af'][:,-1]>0)&(data['af'][:,-1]<1))
        for method,name in enumerate(names):
            selected=data['cdfs'][:,method].copy();ref=neutral[:,method].copy()
            if method<2:selected=np.rint(selected*1000).astype(int);ref=np.rint(ref*1000).astype(int)
            for kind in ('cumulative','interval'):
                x,y=(selected,ref) if kind=='cumulative' else (np.diff(selected),np.diff(ref))
                p=np.stack([(1+np.count_nonzero(y>=values,axis=0))/101 for values in x])
                np.testing.assert_array_equal(p,data['p_'+kind][:,method])
                for alpha in (.05,.01):
                    mask=p<alpha;count=int(mask.any(1).sum())
                    row=next(r for r in summary if r['case']==key and r['method']==name and r['threshold']==alpha)
                    assert row['n_selected']==n and row[kind+'_any_count']==count
                    b=next(r for r in benchmark if r['case']==key and r['method']==name and r['threshold']==alpha and r['statistic']==kind)
                    assert b['true_positives']==count and b['false_negatives']==n-count
                    assert b['false_positives']+b['true_negatives']==100
                    if n==100:
                        denominator=2*count+b['false_positives']+(n-count)
                        expected=2*count/denominator if denominator else None
                        assert b['f1']==expected
                    else:assert b['f1'] is None
    pngs=list(out.glob('*.png'));pdfs=list(out.glob('*.pdf'))
    assert len(pngs)==13*n_cases+1 and len(pdfs)==len(pngs)+9
    for path in pngs:
        with Image.open(path) as im:assert im.size==(1980,1530);im.verify()
    pages=0
    for path in pdfs:
        content=path.read_bytes()
        assert content.startswith(b'%PDF-') and content.rstrip().endswith(b'%%EOF')
        assert b'/FontFile2' in content and not re.search(rb'/Subtype\s*/Image',content)
        expected=3*n_cases if path.name.startswith('all_statistic') else n_cases if path.name.startswith('posterior_mass_') or path.stem=='allele_frequency_selected_vs_neutral' else 1
        count=len(re.findall(rb'/Type\s*/Page\b',content));assert count==expected,(path.name,count,expected)
        pages+=count
    assert pages==30*n_cases+1
    print(f'PASS: six-treatment accounting, segregation, independent empirical ranks, benchmark counts, hashes, {len(pngs)} PNGs and {len(pdfs)} vector PDFs ({pages} pages).')


if __name__=='__main__':main()
