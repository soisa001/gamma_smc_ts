#!/usr/bin/env python3
"""Audit continuation seed ranges, retained simulations, and proposal accounting."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import numpy as np
import run_selection_arrays as arrays

ROOT = Path(__file__).resolve().parents[1]
BASE = arrays.base


def audit(raw, sample=0):
    spec = json.loads((raw/'manifest.json').read_text())
    signature = BASE.digest(raw/'manifest.json')
    strata = {s['key']:s for s in spec['strata']}
    counts = {k:dict(submitted=0, checked=0, completed=0, lost=0, fixed=0, accepted=0, prefix=0, prefix_lost=0, prefix_fixed=0, prefix_accepted=[]) for k in strata}
    ids = set(); onsets = {}; founder_accepts = set(); gap = set()
    for line in (raw/'dispatch.jsonl').read_text().splitlines():
        event = json.loads(line); key=event['key']; s=strata[key]; c=counts[key]; task=event['task']
        assert task['id'] not in ids; ids.add(task['id'])
        assert event['ordinal']==c['submitted']; c['submitted']+=1
        assert task['selection']==s['selection'] and task['origin']==s['origin']
        assert task['cohort']=='selected' and task['regime']==s['regime']
        if s['origin']=='de_novo':
            slot, attempt=task['slot'],task['attempt']
            assert attempt>=s['slot_attempt_starts'][slot]
            assert event['founder']==s['founders'][slot]['path']
            assert task['unit']==f"de_novo/{s['old_key']}/slot_{slot:03d}/attempt_{attempt:05d}"
        else:
            assert task['replicate']>=s['proposal_start']
            assert task['unit']==f"introgression/s{s['selection']:g}/selected_{task['replicate']:05d}"
        directory=raw/task['id']
        receipt=next((directory/n for n in ('complete.json','rejected.json') if (directory/n).exists()),None)
        if receipt is None: gap.add(key); continue
        r=json.loads(receipt.read_text())
        assert r['config_sha256']==signature
        assert r['seed']==BASE.seed(s['config']['seed_namespace']+task['unit'])
        assert r['seed_namespace']==s['config']['seed_namespace']
        if s['origin']=='de_novo': assert r['onset_sha256']==s['founders'][slot]['onset_sha256']
        c['completed']+=1
        if not sample or c['checked']<sample or r['accepted']:
            assert BASE.digest(directory/'trajectory.tsv')==r['trajectory_sha256']
            onset=(directory/'onset.trees').resolve()
            if onset not in onsets: onsets[onset]=BASE.digest(onset)
            assert onsets[onset]==r['onset_sha256']
            trajectory=np.loadtxt(directory/'trajectory.tsv',skiprows=1)
            np.testing.assert_allclose(trajectory[:,1],trajectory[:,4]/(2*trajectory[:,2]),atol=1e-14)
            assert trajectory[-1,1]==r['final_frequency']
            if s['origin']=='de_novo': assert trajectory[0,4]==1
            c['checked']+=1
        if r['accepted']:
            assert 0<r['final_frequency']<1 and r['stop_years_ago']==0
            arrays.ancestor.verify(directory,r)
            if s['origin']=='de_novo':
                assert (key,slot) not in founder_accepts; founder_accepts.add((key,slot))
            c['accepted']+=1
        else:
            assert r['rejection'] in ('loss','fixation')
            assert r['final_frequency']==(0 if r['rejection']=='loss' else 1)
            c['lost' if r['rejection']=='loss' else 'fixed']+=1
        if key not in gap:
            c['prefix']+=1
            if r['accepted']: c['prefix_accepted'].append(str(directory))
            else: c['prefix_lost' if r['rejection']=='loss' else 'prefix_fixed']+=1
    final=raw/'finished.json'
    if final.exists():
        status=json.loads(final.read_text())
        for row in status['strata']:
            key=row['key']; c=counts[key]; s=strata[key]
            for left,right in [('new_submitted','submitted'),('new_completed','completed'),('completed_prefix','prefix'),
                               ('prefix_lost','prefix_lost'),('prefix_fixed','prefix_fixed'),('all_completed_survivors','accepted')]:
                assert row[left]==c[right], (key,left,row[left],c[right])
            kept=json.loads((raw/(key+'_cohort.json')).read_text())['accepted']
            expected=([r['source_directory'] for r in s['reused']]+c['prefix_accepted'])[:100]
            assert [r['source_directory'] for r in kept]==expected
            assert len(kept)==row['accepted']
            for r in kept:
                assert 0<r['final_frequency']<1
                arrays.ancestor.verify(Path(r['source_directory']),r)
    return dict(validated=True,final=final.exists(),sample_per_stratum=sample,
        manifest_sha256=signature,strata=counts,
        checks='Unique proposal identities; unused continuation seeds; unchanged configurations; source-founder hashes; trajectory copy counts; strict population segregation; all retained output hashes; completed-prefix accounting.')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--raw',type=Path,default=ROOT.parent/'sim/selected_survivor_extension_20261112')
    p.add_argument('--sample',type=int,default=0,help='Zero checks every completed proposal; positive checks first N plus all survivors.')
    args=p.parse_args(); assert args.sample>=0
    print(json.dumps(audit(args.raw.resolve(),args.sample),indent=2))


if __name__=='__main__': main()
