#!/usr/bin/env python3
"""Independently audit the cached 100-versus-100 neutral FPR calculation."""
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
    parser.add_argument('--out',type=Path,default=Path(__file__).resolve().parents[1]/'docs/results/segregating_neutral_fpr_100_20261110')
    args=parser.parse_args();out=args.out
    digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    inventory=json.loads((out/'output_hashes.json').read_text())
    assert all(digest(out/name)==value for name,value in inventory.items())
    summary=json.loads((out/'summary.json').read_text())
    with np.load(out/'parameters.npz') as src:data={k:src[k].copy() for k in src.files}
    with np.load(out/'regional_p_values.npz') as src:regional={k:src[k].copy() for k in src.files}
    assert not set(data['test_ids']) & set(data['reference_ids'])
    for key in ('test_af','reference_af'):
        assert data[key].shape==(100,2001)
        assert np.all((data[key][:,-1]>0)&(data[key][:,-1]<1))
    names=['True TMRCA','Posterior median','Posterior-weighted mass']
    # Separate scalar-method implementation: convert discrete fractions before
    # differencing, then count ranks independently at every position and time.
    for method,name in enumerate(names):
        test=data['test_all_cdfs'][:,method].copy()
        ref=data['reference_all_cdfs'][:,method].copy()
        if method<2:
            test=np.rint(test*1000).astype(int);ref=np.rint(ref*1000).astype(int)
        for kind in ('cumulative','interval'):
            x,y=(test,ref) if kind=='cumulative' else (np.diff(test),np.diff(ref))
            p=np.empty(x.shape)
            for i in range(100):p[i]=(1+np.count_nonzero(y>=x[i],axis=0))/101
            np.testing.assert_array_equal(p,regional['p_'+kind][:,method])
            np.testing.assert_array_equal(p[:,50],data['p_'+kind][:,method])
            for alpha in (.05,.01):
                mask=p<alpha;count=int(mask[:,50].any(1).sum())
                row=next(r for r in summary['results'] if r['method']==name and r['statistic']==kind and r['threshold']==alpha)
                assert row['any_time_false_positive_count']==count
                np.testing.assert_array_equal(row['false_positive_count_by_time'],mask[:,50].sum(0))
                assert abs(row['mean_pointwise_fpr']-mask[:,50].mean())<1e-14
                lo,hi=row['any_time_wilson95'];assert 0<=lo<=count/100<=hi<=1
                regional_row=next(r for r in summary['regional_any_position_and_time'] if r['method']==name and r['statistic']==kind and r['threshold']==alpha)
                assert regional_row['false_positive_count']==int(mask.any(axis=(1,2)).sum())
                np.testing.assert_array_equal(regional_row['any_time_fpr_by_position'],mask.any(2).mean(0))
    assert len(summary['results'])==len(summary['regional_any_position_and_time'])==12
    calls=json.loads((out/'replicate_calls.json').read_text());assert len(calls)==1200
    for row in calls:
        i=list(data['test_ids']).index(row['id']);m=names.index(row['method'])
        assert row['called']==bool((data['p_'+row['statistic']][i,m]<row['threshold']).any())
    pngs=list(out.glob('*.png'));pdfs=list(out.glob('*.pdf'))
    assert len(pngs)==11 and len(pdfs)==13
    for path in pngs:
        with Image.open(path) as im:assert im.size==(1980,1530);im.verify()
    pages=0
    for path in pdfs:
        content=path.read_bytes()
        assert content.startswith(b'%PDF-') and content.rstrip().endswith(b'%%EOF')
        assert b'/FontFile2' in content and not re.search(rb'/Subtype\s*/Image',content)
        expected=6 if path.stem=='neutral_false_calls' else 2 if path.stem=='fpr_pointwise' else 1
        count=len(re.findall(rb'/Type\s*/Page\b',content));assert count==expected
        pages+=count
    assert pages==19
    print('PASS: independent regional/focal ranks, FPR summaries, 1,200 replicate calls, AF segregation, hashes, 11 PNGs and 13 vector PDFs (19 pages).')


if __name__=='__main__':main()
