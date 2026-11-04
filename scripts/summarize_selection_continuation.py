#!/usr/bin/env python3
"""Validate a resumed selection array against its earlier report and summarize gains."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def summarize(before, after, started_epoch=None):
    before, after = Path(before), Path(after)
    read = lambda root, name: pd.read_csv(root / f'{name}.tsv', sep='\t')
    old_proposals, new_proposals = (read(p, 'proposals') for p in (before, after))
    old_ids, new_ids = (set(p.id) for p in (old_proposals, new_proposals))
    assert old_ids <= new_ids, 'Continuation omitted previously completed proposals'
    assert new_proposals.id.is_unique
    old_records, new_records = (read(p, 'completed_units').set_index('id') for p in (before, after))
    assert set(old_records.index) <= set(new_records.index)
    unchanged = ['seed', 'observed_sha256', 'posterior_sha256', 'scores_sha256', 'receipt_sha256']
    pd.testing.assert_frame_equal(old_records[unchanged], new_records.loc[old_records.index, unchanged])
    with np.load(before / 'calibration.npz') as old, np.load(after / 'calibration.npz') as new:
        assert set(old.files) == set(new.files)
        for key in old.files:
            np.testing.assert_array_equal(old[key], new[key])
    # Reused observations must retain exactly the same calls and metric values.
    keys = ['scan', 'regime', 'origin', 'method', 'id', 'cohort', 'selection']
    old_metrics, new_metrics = (read(p, 'replicate_metrics').set_index(keys) for p in (before, after))
    pd.testing.assert_frame_equal(old_metrics, new_metrics.loc[old_metrics.index], check_exact=True)
    counts = read(before, 'proposal_summary').merge(
        read(after, 'proposal_summary'), on=['regime', 'origin', 'cohort', 'selection'],
        suffixes=('_before', '_after'), validate='one_to_one')
    for field in ['accepted', 'analyzed_prefix', 'lost', 'fixed']:
        counts[f'additional_{field}'] = counts[f'{field}_after'] - counts[f'{field}_before']
        assert (counts[f'additional_{field}'] >= 0).all()
    counts['reached_target'] = counts.accepted_after >= counts.target_after
    counts['segregating_fraction'] = counts.accepted_after / counts.analyzed_prefix_after.replace(0, np.nan)
    counts.to_csv(after / 'continuation_counts.tsv', sep='\t', index=False, na_rep='NA', lineterminator='\n')
    keys = ['scan', 'regime', 'origin', 'method', 'selection', 'metric']
    power = read(before, 'performance').merge(read(after, 'performance'), on=keys,
        suffixes=('_before', '_after'), validate='one_to_one')
    power['rate_change'] = power.rate_after - power.rate_before
    power.to_csv(after / 'continuation_performance.tsv', sep='\t', index=False, na_rep='NA', lineterminator='\n')
    provenance = json.loads((after / 'provenance.json').read_text())
    schedule, finished = provenance['schedule'], provenance['finished']
    continuation_start = schedule['started_epoch'] if started_epoch is None else started_epoch
    assert continuation_start <= schedule['started_epoch']
    assert finished['finished_epoch'] >= schedule['started_epoch']
    assert finished['finished_epoch'] <= schedule['deadline_epoch'] + 10, 'Run exceeded deadline'
    added = new_proposals[~new_proposals.id.isin(old_ids)]
    selected = counts[counts.cohort == 'selected']
    proof = dict(
        before=str(before.resolve()), after=str(after.resolve()), source_sha256=digest(__file__),
        before_provenance_sha256=digest(before / 'provenance.json'),
        after_provenance_sha256=digest(after / 'provenance.json'),
        previous_completed_proposals=len(old_ids), total_completed_proposals=len(new_ids),
        additional_completed_proposals=len(added), additional_selected_proposals=int((added.cohort == 'selected').sum()),
        additional_selected_survivors=int(selected.additional_accepted.sum()),
        total_selected_survivors=int(selected.accepted_after.sum()),
        selected_treatments_at_target=int(selected.reached_target.sum()), selected_treatments=len(selected),
        elapsed_seconds=finished['finished_epoch'] - continuation_start,
        final_invocation_elapsed_seconds=finished['finished_epoch'] - schedule['started_epoch'],
        started_epoch=continuation_start, final_invocation_started_epoch=schedule['started_epoch'],
        finished_epoch=finished['finished_epoch'],
        deadline_epoch=schedule['deadline_epoch'], workers=schedule['parameters']['workers'],
        unchanged_checks=['All previous completed proposal IDs retained',
            'Previously analyzed region seeds, tree/posterior/score/receipt hashes identical',
            'All calibration arrays exactly identical', 'All earlier replicate metrics exactly identical'],
        interpretation='Before/after results are nested cohorts, not independent experiments. '
            'Power is conditional on present-day segregation. Zero surviving replicates means power is unestimable.')
    (after / 'continuation_provenance.json').write_text(json.dumps(proof, indent=2) + '\n', encoding='utf-8', newline='\n')
    lines = ['# Four-hour continuation', '',
        f"Added {len(added):,} completed proposals and {proof['additional_selected_survivors']:,} selected survivors in {proof['elapsed_seconds']/3600:.3f} hours using {proof['workers']} workers.",
        f"Total: {proof['total_selected_survivors']:,} selected survivors; {proof['selected_treatments_at_target']}/{len(selected)} treatments reached 100.",
        '', 'Earlier proposals, analyzed output hashes, calibration arrays, and replicate metrics were checked against the previous report. Seeds and simulation settings were unchanged.',
        '', '| Selection period | Origin | s | Before | After / target | New proposals in analyzed prefix | Lost total | Fixed total |',
        '|---|---|---:|---:|---:|---:|---:|---:|']
    for r in selected.itertuples():
        lines.append(f'| {r.regime} | {r.origin} | {r.selection:g} | {r.accepted_before} | {r.accepted_after}/{r.target_after} | {r.additional_analyzed_prefix} | {r.lost_after} | {r.fixed_after} |')
    lines += ['', 'Power is conditional on present-day segregation. Zero survivors does not establish zero detection power. Before/after estimates share the earlier observations, so their changes are not independent comparisons.',
        '', 'The full method comparison and confidence intervals are in `performance.tsv`; `continuation_performance.tsv` contains the nested before/after comparison. `continuation_counts.tsv` reports all selected and neutral strata.']
    (after / 'CONTINUATION.md').write_text('\n'.join(lines) + '\n', encoding='utf-8', newline='\n')
    print(json.dumps(proof, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--before', required=True)
    parser.add_argument('--after', required=True)
    parser.add_argument('--started-epoch', type=float, help='First continuation invocation, including any restart overhead')
    args = parser.parse_args()
    summarize(args.before, args.after, args.started_epoch)
