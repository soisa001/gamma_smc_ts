#!/usr/bin/env python3
"""Independently audit cached ranks, cohort counts, and comparison artifacts."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import argparse
import hashlib
import json
from pathlib import Path
import re
import numpy as np
from PIL import Image


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, default=Path(__file__).resolve().parents[1] / 'docs/results/denovo_comparison_20261111')
    args = p.parse_args(); out = args.out
    digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    meta = json.loads((out / 'manifest.json').read_text())
    inventory = json.loads((out / 'output_hashes.json').read_text())
    assert all(digest(out / name) == value for name, value in inventory.items())
    summary = json.loads((out / 'summary.json').read_text())['results']
    assert len(summary) == len(meta['cases']) * 24
    data = {}
    for key, entry in meta['cohorts'].items():
        assert digest(out / entry['cache']) == entry['cache_sha256']
        with np.load(out / entry['cache'], allow_pickle=False) as src:
            d = {k: src[k].copy() for k in src.files}
        data[key] = d
        n = len(d['ids']); assert n == entry['n'] and len(set(d['ids'])) == n
        assert d['cdfs'].shape == (n, 3, 50) and d['af'].shape == (n, len(d['af_times_kya']))
        assert np.isfinite(d['cdfs']).all() and np.isfinite(d['af']).all()
        assert np.all((d['af'][:, -1] > 0) & (d['af'][:, -1] < 1))
        assert np.all((d['cdfs'] >= 0) & (d['cdfs'] <= 1+1e-12))
        assert np.all(np.diff(d['cdfs']) >= -1e-12)
    founders, observed = [], []
    for key, entry in meta['cohorts'].items():
        if not key.startswith('de_novo_'): continue
        records = entry['provenance']['records']
        assert len(records) == entry['n']
        assert all(r['origin']=='de_novo' and 0<r['final_frequency']<1 for r in records)
        founders.extend(r['onset_sha256'] for r in records)
        observed.extend(r['observed_sha256'] for r in records)
    assert len(set(founders)) == len(founders), 'Retained simulations share a founder genealogy'
    assert len(set(observed)) == len(observed), 'Duplicate retained simulation'
    if not meta['live_snapshot']:
        for row in meta['progress']['strata']:
            key = 'de_novo_' + row['key'].replace('.', 'p')
            assert row['accepted'] == meta['cohorts'][key]['n']
            assert row['reused'] == len(meta['simulation']['reused'][row['key']])
            assert row['accepted']-row['reused']+row['lost']+row['fixed'] == row['new_completed']
    names = ['True TMRCA', 'Posterior median', 'Posterior-weighted mass']
    for case in meta['cases']:
        for origin in ('de_novo', 'introgression'):
            selected = data[case['selected'][origin]]; null = data[case['neutral'][origin]]
            n, nref = len(selected['ids']), len(null['ids'])
            assert set(selected['ids']).isdisjoint(null['ids'])
            if origin == 'de_novo':
                expected_age = 10 if case['regime']=='10k' else 50
                assert selected['af_times_kya'][0] == null['af_times_kya'][0] == expected_age
            for method, name in enumerate(names):
                x, y = selected['cdfs'][:, method], null['cdfs'][:, method]
                if method < 2: x, y = np.rint(x*1000).astype(int), np.rint(y*1000).astype(int)
                for kind in ('cumulative', 'interval'):
                    xx, yy = (x, y) if kind=='cumulative' else (np.diff(x), np.diff(y))
                    if nref:
                        expected = np.array([(1 + (yy >= row).sum(0))/(nref+1) for row in xx]).reshape(n, xx.shape[-1])
                    else: expected = np.full((n, xx.shape[-1]), np.nan)
                    np.testing.assert_array_equal(selected['p_' + kind][:, method], expected)
                    for alpha in (.05, .01):
                        row = next(r for r in summary if r['case']==case['key'] and r['origin']==origin and r['method']==name and r['statistic']==kind and r['alpha']==alpha)
                        attainable = nref > 0 and 1/(nref+1) < alpha
                        assert row['threshold_attainable'] == attainable
                        assert row['selected_n'] == n and row['neutral_n'] == nref
                        count = int((expected < alpha).any(1).sum()) if n and attainable else None
                        assert row['any_time_count'] == count
                        assert row['any_time_fraction'] == (count/n if count is not None else None)
    pngs, pdfs = list(out.glob('*.png')), list(out.glob('*.pdf'))
    zooms = sum(c['regime']=='10k' and c['selection']<.05 for c in meta['cases'])
    assert len(pngs) == 10*len(meta['cases'])+10+zooms and len(pdfs) == len(pngs)+9
    for path in pngs:
        with Image.open(path) as im:
            assert im.size == (1650, 1275)
            im.verify()
    pages = 0
    for path in pdfs:
        content = path.read_bytes()
        assert content.startswith(b'%PDF-') and content.rstrip().endswith(b'%%EOF')
        assert b'/FontFile2' in content and not re.search(rb'/Subtype\s*/Image', content)
        count = len(re.findall(rb'/Type\s*/Page\b', content)); assert count >= 1
        pages += count
    assert pages == 23*len(meta['cases'])+20+2*zooms
    result = dict(validated=True, cases=len(meta['cases']), p_value_rows=len(summary), pngs=len(pngs), vector_pdfs=len(pdfs), pdf_pages=pages,
                  checks='Input/output hashes; unique IDs; final population segregation; age-matched nulls; independent integer-tie-aware ranks; summary counts; image integrity; editable PDF text and no raster images.')
    print(json.dumps(result, indent=2))


if __name__ == '__main__': main()
