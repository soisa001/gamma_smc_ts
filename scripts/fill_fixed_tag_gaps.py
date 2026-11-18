"""Fill the seven historical EAS survivor-cohort gaps with independent archived fixations."""
import argparse
import csv
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import shutil

import report_cumulative_grid as prior
import numpy as np

ROOT = prior.ROOT
HISTORY = ROOT/'docs/results/selected_survivor_extension_20261117'
RAW = Path.home()/'sim/fixed_tag_gaps_20261118'
OUT = ROOT/'docs/results/fixed_tag_gaps_20261118'


def plan(raw):
    import replay_fixed_tag
    raw.mkdir(parents=True, exist_ok=True)
    cases = prior.read(HISTORY/'manifest.json')['cases']
    covered = {c['key']: {r['onset_sha256'] for r in c['records']} for c in cases}
    pools = {c['key']: [] for c in cases}
    cfg_path = prior.OUT/'simulation_manifest.json'
    contract = dict(configs=prior.read(cfg_path)['configs'], baseline=prior.sha(cfg_path),
                    previous=prior.sha(HISTORY/'manifest.json'), replay_sha256=prior.sha(Path(replay_fixed_tag.__file__)),
                    rule='Last generation with focal AF <= 0.8 before fixation. Rank historical biallelic SNPs across 1 Mb by r2, distance, then position. Use first segregating tag among 100 present-day diploids.',
                    historical_diploids=1000, ledgers={}, cases=[])
    for name, origin in [('strong_selection_20261110','introgression'), ('denovo_comparison_20261111','de_novo'),
                         ('denovo_comparison_cached_20261111','de_novo')]:
        path = ROOT.parent/'sim'/name/'proposal_ledger.json'
        ledger = prior.read(path)
        contract['ledgers'][str(path)] = prior.sha(path)
        for c in cases:
            if c['origin'] != origin:
                continue
            key = c['key']
            rows = ledger[f"{c['regime']}_s{c['selection']:g}"]
            for r in rows:
                if r['rejection'] != 'fixation' or r['onset_sha256'] in covered[key]:
                    continue
                assert (r['origin'],r['regime'],r['selection']) == (origin,c['regime'],c['selection'])
                pools[key].append(r['source_directory'])
                covered[key].add(r['onset_sha256'])
                if len(pools[key]) >= 130:
                    break
    for c in cases:
        key = c['key']
        assert len(pools[key]) >= 100-c['n'], key
        contract['cases'].append(dict(key=key, origin=c['origin'], regime=c['regime'], s=c['selection'],
                                      kept=c['n'], cache=c['cache'], cache_sha256=c['cache_sha256'], sources=pools[key]))
    manifest = raw/'manifest.json'
    if manifest.exists():
        assert prior.read(manifest) == contract
    else:
        prior.save(manifest, contract)
    return contract


def run(raw, workers):
    from replay_fixed_tag import replay
    manifest = raw/'manifest.json'
    contract = plan(raw)
    signature = prior.sha(manifest)
    pairs = ROOT.parent/'sim/cumulative_grid_20261117/pairs.tsv'
    assert prior.sha(pairs) == prior.read(prior.OUT/'simulation_manifest.json')['input_hashes']['pairs.tsv']
    for cfg in contract['configs'].values():
        for src, key in ((cfg['slim'],'slim_sha256'),(cfg['decoder'],'decoder_sha256')):
            assert prior.sha(Path(src)) == cfg[key]
    cases = contract['cases']
    next_id = {c['key']:0 for c in cases}
    counts = {c['key']:c['kept'] for c in cases}
    with ProcessPoolExecutor(max_workers=workers) as pool:
        while any(n < 100 for n in counts.values()):
            futures = {}
            for offset in range(100):
                for c in cases:
                    key = c['key']
                    if offset >= 100-counts[key]:
                        continue
                    idx = next_id[key]
                    if idx >= len(c['sources']):
                        raise ValueError(f'Insufficient distinct tagged histories for {key}')
                    next_id[key] += 1
                    future = pool.submit(replay, c['sources'][idx], raw/key/f'{idx:03d}', contract['configs'][c['regime']], pairs, signature)
                    futures[future] = key
            for completed, future in enumerate(as_completed(futures),1):
                if future.result()['tag'] is not None:
                    counts[futures[future]] += 1
                if completed % 25 == 0:
                    print(counts, flush=True)
    print(counts, flush=True)


def prepare(raw):
    manifest = raw/'manifest.json'
    contract = prior.read(manifest)
    signature = prior.sha(manifest)
    curves, records = {}, []
    for c in contract['cases']:
        key = c['key']
        cache = HISTORY/c['cache']
        assert prior.sha(cache) == c['cache_sha256']
        with np.load(cache) as data:
            masses = list(data['cdfs'][:,2])
            ids = list(map(str,data['ids']))
        assert len(masses) == c['kept']
        for i, src in enumerate(c['sources']):
            if len(masses) == 100:
                break
            folder = raw/key/f'{i:03d}'
            rec = prior.read(folder/'tag.json')
            assert rec['source']==src and rec['contract']==signature
            records.append(dict(case=key, **rec))
            if rec['tag'] is None:
                continue
            assert prior.sha(folder/'scores.npz') == rec['scores_sha256']
            with np.load(folder/'scores.npz') as score:
                np.testing.assert_array_equal(score['positions'],[rec['tag']['position']])
                np.testing.assert_array_equal(score['cutoffs'],np.arange(1000,50001,1000))
                masses.append(score['cdfs'][2,0])
            ids.append(src)
        assert len(masses) == len(set(ids)) == 100
        curves[key] = np.array(masses)
        curves[key+'_ids'] = np.array(ids)
    already_complete = []
    old = prior.read(prior.BASE/'manifest.json')['cohorts']
    for origin in prior.ORIGINS:
        key = origin+'_50k_s0p01'
        entry = old[key]
        path = prior.BASE/entry['cache']
        assert prior.sha(path)==entry['cache_sha256']
        with np.load(path) as bank:
            curves[key] = bank['cdfs'][:,2]
            curves[key+'_ids'] = bank['ids']
        already_complete.append(dict(key=key, origin=origin, regime='50k', s=.01, kept=100,
                                     cache=entry['cache'], cache_sha256=entry['cache_sha256']))
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT/'.gitattributes').write_text('* -text\n')
    shutil.copyfile(manifest, OUT/'simulation_manifest.json')
    np.savez_compressed(OUT/'mass.npz', **curves)
    prior.save(OUT/'inputs.json', dict(records=records, already_complete=already_complete, cache_sha256=prior.sha(OUT/'mass.npz'),
                                      baseline_mass_sha256=prior.sha(prior.OUT/'mass.npz')))


def report():
    inputs = prior.read(OUT/'inputs.json')
    assert prior.sha(OUT/'mass.npz')==inputs['cache_sha256'] and prior.sha(prior.OUT/'mass.npz')==inputs['baseline_mass_sha256']
    rows, detail = [], {}
    contract = prior.read(OUT/'simulation_manifest.json')
    with np.load(OUT/'mass.npz') as queries, np.load(prior.OUT/'mass.npz') as banks:
        for c in contract['cases'] + inputs['already_complete']:
            origin, regime, s = c['origin'], c['regime'], c['s']
            refkey = 'introgression' if origin=='introgression' else ('de_novo_10k' if regime=='10k' else 'de_novo_50k')
            ref = banks[refkey+'_neutral']
            null_ranks = 1+(ref[:,None,:]>=banks[refkey+'_test'][None,:,:]).sum(axis=0)
            cur = queries[c['key']]
            assert cur.shape==ref.shape==(100,50) and np.isfinite(cur).all() and np.diff(cur,axis=1).min()>=-1e-12
            p = (1+(ref[:,None,:]>=cur[None,:,:]).sum(axis=0))/101
            detail[c['key']] = p
            lo, hi = prior.WINDOWS[regime]
            for threshold in (.05,.01):
                call = (p[:,lo-1:hi]<threshold).any(axis=1)
                kept = c['kept']
                rows.append(dict(origin=origin,window=regime,s=s,alpha=threshold,original_n=kept,original_calls=int(call[:kept].sum()),
                                 tagged_n=100-kept,tagged_calls=int(call[kept:].sum()),total_calls=int(call.sum()),total_n=100,
                                 segregating_neutral_calls=int((null_ranks[:,lo-1:hi]/101<threshold).any(axis=1).sum()),
                                 segregating_neutral_n=100))
    np.savez_compressed(OUT/'pointwise_p.npz', **detail)
    with (OUT/'rates.tsv').open('w',newline='') as f:
        w = csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    import matplotlib.pyplot as plt
    plt.rcParams.update({'pdf.fonttype':42, 'ps.fonttype':42, 'font.size':14})
    fig, ax = plt.subplots(figsize=(11,8.5), layout='constrained')
    treatments = sorted({(r['origin'], r['window'], r['s']) for r in rows},
                        key=lambda t: (t[0], list(prior.WINDOWS).index(t[1]),t[2]))
    labels = []
    for y,(origin, regime, s) in enumerate(treatments):
        entries = [r for r in rows if (r['origin'], r['window'], r['s'])==(origin, regime, s)]
        model = 'De novo' if origin=='de_novo' else 'Introgression'
        window = dict(zip(prior.WINDOWS,prior.LABELS))[regime]
        labels.append(f'{model}    {window}    s = {s:g}')
        for j,r in enumerate(entries):
            for offset, (field, palette) in enumerate([('total_calls','Blues'),('segregating_neutral_calls','Oranges')]):
                x = 2*j+offset
                n = r[field]
                ax.add_patch(prior.Rectangle((x,y),1,1,facecolor=plt.get_cmap(palette)(.10+.87*n/100),
                                             edgecolor='white',linewidth=2))
                ax.text(x+.5,y+.5,f'{n}%',ha='center',va='center',color='white' if n>60 else 'black')
        ax.text(4.5,y+.5,str(entries[0]['tagged_n']),ha='center',va='center')
    ax.set(xlim=(0,5),ylim=(len(treatments),0),yticks=np.arange(len(labels))+.5,yticklabels=labels,
           xticks=np.arange(5)+.5,xticklabels=['Calls\np < .05','FPR\np < .05','Calls\np < .01','FPR\np < .01','Tagged\nper 100'])
    ax.xaxis.tick_top()
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.suptitle('Within-window calls in stronger selection treatments',fontsize=19)
    fig.supxlabel('FPR uses 100 held-out segregating neutrals. The FPR for tags of neutral fixations is unmeasured.',fontsize=13)
    for ext in ('png','pdf'):
        fig.savefig(OUT/f'tagged_detection.{ext}',dpi=160)
    plt.close(fig)
    print(rows,flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=('plan','simulate','prepare','report'))
    parser.add_argument('--raw', type=Path, default=RAW)
    parser.add_argument('--workers', type=int, default=4, choices=range(1,21))
    args = parser.parse_args()
    if args.phase == 'plan':
        plan(args.raw)
    elif args.phase == 'simulate':
        run(args.raw,args.workers)
    elif args.phase == 'prepare':
        prepare(args.raw)
    else:
        report()
