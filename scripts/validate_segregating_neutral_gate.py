#!/usr/bin/env python3
"""Check that s=0 rejection changes cohort membership, not neutral trajectories."""
import argparse
import json
from pathlib import Path
import numpy as np
import tskit
import run_selection_arrays as arrays


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--raw',type=Path,required=True)
    args=p.parse_args()
    raw=args.raw
    manifest=json.loads((raw/'manifest.json').read_text())
    assert arrays.base.digest(arrays.__file__)==manifest['scientific_sources']['run_selection_arrays.py']
    cfg=manifest['config']
    checks=[]
    for name in ('complete.json','rejected.json'):
        receipt=sorted((raw/'introgression').glob('*/'+name))[0]
        record=json.loads(receipt.read_text())
        dest=raw/'validation'/receipt.parent.name
        dest.mkdir(parents=True,exist_ok=True)
        cached=dest/'replay.json'
        contract=dict(source_receipt_sha256=arrays.base.digest(receipt),manifest_sha256=arrays.base.digest(raw/'manifest.json'))
        if cached.exists():
            result=json.loads(cached.read_text());assert result['contract']==contract
        else:
            result=dict(contract=contract,replay=arrays.simulate(cfg,dict(record,cohort='neutral_test'),dest))
            arrays.base.save(cached,result)
        original=np.loadtxt(receipt.parent/'trajectory.tsv',skiprows=1,ndmin=2)
        replay=np.loadtxt(dest/'trajectory.tsv',skiprows=1)
        assert replay.shape==(2001,5)
        np.testing.assert_array_equal(original,replay[:len(original)])
        np.testing.assert_array_equal(replay[:,3],0)
        if record['accepted']:
            assert 0<replay[-1,1]<1
            first=tskit.load(receipt.parent/'observed.trees')
            second=tskit.load(dest/'observed.trees')
            tables=[first.dump_tables(),second.dump_tables()]
            # SLiM's nonspatial individual coordinates can contain subnormal
            # values that differ across processes; they do not enter ancestry,
            # mutation simulation or the decoder. Compare all other columns.
            for table in tables:
                assert table.metadata['SLiM']['spatial_dimensionality']==''
                assert np.isfinite(table.individuals.location).all()
                assert np.max(np.abs(table.individuals.location))<1e-100
                table.individuals.location=np.zeros_like(table.individuals.location)
            assert tables[0].equals(tables[1],ignore_provenance=True)
            outcome='Exact full trajectory, genealogy, mutation tables and sample identities match the unconditioned s=0 replay. Only timestamp provenance and unused nonspatial individual coordinates are excluded.'
        else:
            assert replay[-1,1] in (0,1)
            np.testing.assert_array_equal(replay[len(original)-1:,1],original[-1,1])
            outcome='Exact trajectory prefix matches; after the rejection point, the full unconditioned replay stays in the same absorbing state.'
        checks.append(dict(id=record['id'],source_receipt=str(receipt),source_receipt_sha256=arrays.base.digest(receipt),
            replay_directory=str(dest),outcome=outcome,trajectory_rows=len(original)))
    arrays.base.save(raw/'neutral_gate_validation.json',dict(checks=checks,passed=True,
        manifest_sha256=arrays.base.digest(raw/'manifest.json'),script_sha256=arrays.base.digest(__file__)))
    print(json.dumps(checks,indent=2),flush=True)


if __name__=='__main__':main()
