#!/usr/bin/env python3
"""Audit a conditioned-null report against the unchanged selected simulation data."""
import argparse
import json
from pathlib import Path
import numpy as np
from pointwise_time_tests import exceedance_counts
from plot_selection_episode_series import sha,save_json

ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--old',type=Path,default=ROOT/'docs/results/selection_calls_100_20261110')
    p.add_argument('--new',type=Path,default=ROOT/'docs/results/selection_calls_100_segregating_null_20261110')
    args=p.parse_args()
    before=json.loads((args.old/'summary.json').read_text())['results']
    after=json.loads((args.new/'summary.json').read_text())['results']
    new_meta=json.loads((args.new/'manifest.json').read_text())
    assert new_meta['neutral_conditioned'] is True
    with np.load(args.new/'neutral.npz') as data:
        null=data['cdfs'].copy();ids=data['ids'].copy()
        pairs=data['pairs'].copy()
    with np.load(args.old/'neutral.npz') as data:np.testing.assert_array_equal(pairs,data['pairs'])
    with np.load(args.new/'neutral_af.npz') as data:
        np.testing.assert_array_equal(ids,data['ids'])
        assert np.all((data['trajectories'][:,-1,1]>0)&(data['trajectories'][:,-1,1]<1))
        np.testing.assert_array_equal(data['trajectories'][:,:,3],0)
    assert null.shape==(100,3,50) and len(set(ids))==100
    receipt=json.loads((args.new/'neutral_cohort_receipt.json').read_text())
    assert len({r['seed'] for r in receipt['accepted']})==100
    for entry in new_meta['cases']:
        with np.load(args.old/entry['cache']) as old,np.load(args.new/entry['cache']) as new:
            for name in ('af','cdfs','ids','proposal_indices','cutoffs_kya','interval_left_kya','interval_right_kya'):
                np.testing.assert_array_equal(old[name],new[name])
            counts=exceedance_counts(new['cdfs'],null)
            for kind,count in zip(('cumulative','interval'),counts):
                np.testing.assert_array_equal(new['neutral_ge_'+kind],count)
                np.testing.assert_array_equal(new['p_'+kind],(count+1)/101)
            for record in (r for r in after if r['case']==entry['case']['key']):
                method=['True TMRCA','Posterior median','Posterior-weighted mass'].index(record['method'])
                for kind in ('cumulative','interval'):
                    mask=new['p_'+kind][:,method]<record['threshold']
                    assert int(mask.any(1).sum())==record[kind+'_any_count']
    index={(r['case'],r['method'],r['threshold']):r for r in before}
    comparison=[]
    for result in after:
        key=(result['case'],result['method'],result['threshold'])
        original=index[key]
        comparison.append(dict(case=key[0],method=key[1],threshold=key[2],
            old_cumulative=original['cumulative_any_count'],new_cumulative=result['cumulative_any_count'],
            old_interval=original['interval_any_count'],new_interval=result['interval_any_count']))
    save_json(args.new/'conditioning_comparison.json',dict(results=comparison,
        validation='All selected AF/CDF values, cohort IDs, proposal indices and time grids are exactly unchanged. Same pair panel. All 100 neutral focal alleles segregate today at s=0, with unique IDs/seeds. P values and summary call counts independently recomputed.',
        old_summary_sha256=sha(args.old/'summary.json'),new_summary_sha256=sha(args.new/'summary.json'),
        script_sha256=sha(Path(__file__)),interpretation='Fractions with any pointwise call, not scan-adjusted power or FPR.'))
    save_json(args.new/'output_hashes.json',{f.name:sha(f) for f in sorted(args.new.iterdir()) if f.is_file() and f.name!='output_hashes.json'})
    print(json.dumps([r for r in comparison if r['method']=='Posterior-weighted mass'],indent=2))


if __name__=='__main__':main()
