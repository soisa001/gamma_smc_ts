"""Continue the active phased run after verified label recovery finishes."""
import argparse
import fcntl
import json
import os
from pathlib import Path
import subprocess
import sys

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--out',type=Path,required=True)
args=p.parse_args()
repo=Path(__file__).resolve().parents[1]
with (args.out/'labels.lock').open('r') as lock:
    fcntl.flock(lock,fcntl.LOCK_SH)
    state=json.loads((args.out/'label_status.json').read_text())
    if state['state']!='complete' or state['completed']!=1200 or state['failed']:
        raise RuntimeError('Label recovery is incomplete or failed; no downstream run launched')
for phase in ('profiles','analyse'):
    print(f'Starting {phase}',flush=True)
    subprocess.run(['bash',str(repo/'scripts/launch_joint_scan.sh'),phase],check=True,
                   env={**os.environ,'JOINT_SCAN_OUTPUT_DIR':str(args.out)})
subprocess.run([sys.executable,str(repo/'scripts/audit_joint_scan.py'),'--out',str(args.out)],check=True)
print('Joint scan, held-out evaluation and audit complete',flush=True)
