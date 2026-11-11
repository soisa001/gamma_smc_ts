"""Replay a retained simulation using the cached founder installation path."""
import argparse
import os,json,shutil,time
from pathlib import Path
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import numpy as np,tskit
import run_denovo_cached_founders as cached
from report_selection_arrays import local_path
root=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--raw',type=Path,default=root.parent/'sim/denovo_comparison_20261111')
parser.add_argument('--out',type=Path,default=root.parent/'sim/denovo_cached_validation_20261111')
args=parser.parse_args()
raw=args.raw.resolve()
cfg=json.loads((raw/'manifest.json').read_text())['configs']['10k']
rows=json.loads((raw/'10k_s0.01_cohort.json').read_text())['accepted']
r=next(r for r in rows if local_path(r['source_directory']).is_relative_to(raw))
source=local_path(r['source_directory'])
out=args.out.resolve();out.mkdir(parents=True,exist_ok=True)
founder=out/'founder';founder.mkdir(exist_ok=True)
for name in ['onset.trees','onset_provenance.json']:shutil.copy2(source/name,founder/name)
p=json.loads((founder/'onset_provenance.json').read_text())
t,ne=cached.base.history();age=400
sizes=[max(2,round(float(ne[np.searchsorted(t,age-i,side='right')-1]))) for i in range(age+1)]
meta=dict(age_generations=age,carrier_ids=p['carrier_slim_ids'],sizes=sizes)
(founder/'founder.json').write_text(json.dumps(meta))
os.environ['LD_LIBRARY_PATH']=str(Path(cfg['slim']).parent.parent/'lib')
cached.arrays.ancestor.onset=lambda c,o,l,d:cached.install_founder(founder,c,o,l,d)
start=time.time()
result=cached.arrays.simulate(cfg,r,out)
np.testing.assert_array_equal(np.loadtxt(out/'trajectory.tsv',skiprows=1),np.loadtxt(source/'trajectory.tsv',skiprows=1))
a=tskit.load(out/'observed.trees');b=tskit.load(source/'observed.trees')
for name in ('nodes','edges','mutations','sites','individuals','populations','migrations'):
    left=getattr(a.tables,name).asdict();right=getattr(b.tables,name).asdict()
    assert left.keys()==right.keys()
    for key in left:
        if name=='individuals' and key=='location':continue  # Unused, uninitialized coordinates in the nonspatial SLiM model.
        if isinstance(left[key],np.ndarray):np.testing.assert_array_equal(left[key],right[key],err_msg=name+':'+key)
        else:assert left[key]==right[key],name+':'+key
assert a.metadata==b.metadata and a.sequence_length==b.sequence_length and a.time_units==b.time_units
assert result['accepted'] and result['onset_sha256']==r['onset_sha256']
receipt=dict(passed=True,source_id=r['id'],source_directory=str(source),trajectory_exact=True,genealogy_and_mutation_tables_exact=True,excluded_fields=['provenance','unused nonspatial individual.location'],onset_hash_exact=True,seconds=time.time()-start,wrapper_sha256=cached.base.digest(cached.__file__))
(out/'validation.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
