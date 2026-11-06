"""Read-only one-minute ledger monitor; saves compact continuation progress."""
import collections
import datetime
import json
from pathlib import Path
import time

root = Path(__file__).resolve().parent
baseline = root / 'run_history' / 'before_4h_20261104T2123Z'
old_ids = {t['id'] for t in json.loads((baseline / 'tasks.json').read_text())}
schedule = json.loads((root / 'schedule.json').read_text())
log = root / 'continuation_progress.jsonl'
while True:
    try:
        states = json.loads((root / 'status.json').read_text())
        tasks = json.loads((root / 'tasks.json').read_text())
        finished = json.loads((root / 'finished.json').read_text())
    except (json.JSONDecodeError, OSError):
        time.sleep(1)
        continue
    state_by_id = {r['id']: r for r in states}
    counts = collections.Counter(r['state'] for r in states)
    selected = collections.Counter()
    proposals = collections.Counter()
    for t in tasks:
        if t['cohort'] == 'selected' and t['id'] in state_by_id:
            key = f"{t['regime']}/{t['origin']}/s{t['selection']:g}"
            proposals[key] += 1
            selected[key] += state_by_id[t['id']]['state'] == 'accepted'
    done = finished['finished_epoch'] >= schedule['started_epoch']
    added = [r for r in states if r['id'] not in old_ids]
    row = dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        elapsed_minutes=round((time.time()-schedule['started_epoch'])/60, 2),
        submitted=len(tasks), returned=len(states), states=dict(counts),
        new_completed=len(added), new_accepted=sum(r['state']=='accepted' for r in added),
        selected_accepted=sum(selected.values()), at_target=sum(n>=100 for n in selected.values()),
        selected_counts=dict(selected), proposal_counts=dict(proposals), finished=done)
    with log.open('a', encoding='utf-8', newline='\n') as stream:
        stream.write(json.dumps(row)+'\n')
    print(json.dumps({k:v for k,v in row.items() if k not in ('selected_counts','proposal_counts')}), flush=True)
    if done or time.time()>schedule['deadline_epoch']+90:
        break
    time.sleep(60)
