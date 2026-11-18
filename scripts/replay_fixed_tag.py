"""Replay one fixed simulation through the present and score its best observable historical tag."""
import json
import os
from pathlib import Path
import subprocess

import run_selection_arrays as arrays
import msprime
import numpy as np
import pyslim
import tskit

base = arrays.base


def choose_tag(ts, historical, modern, focal_carriers, center):
    """Rank by historical haplotype r², then distance, then coordinate."""
    n = len(historical)
    n_carriers = int(focal_carriers.sum())
    assert 0 < n_carriers < n
    ranked = []
    for var in ts.variants(samples=np.r_[historical, modern], isolated_as_missing=False):
        pos = int(var.site.position)
        if pos == center:
            continue
        past = var.genotypes[:n]
        alleles = np.unique(past)
        if len(alleles) != 2:
            continue
        carries = past == alleles[1]
        count = int(carries.sum())
        overlap = int(np.count_nonzero(carries & focal_carriers))
        numerator = n * overlap - count * n_carriers
        if numerator == 0:
            continue
        r2 = numerator**2 / (count * (n-count) * n_carriers * (n-n_carriers))
        allele = int(alleles[1] if numerator > 0 else alleles[0])
        modern_count = int(np.count_nonzero(var.genotypes[n:] == allele))
        ranked.append((round(r2, 14), abs(pos-center), pos, var.alleles[allele],
                       count if numerator > 0 else n-count, modern_count))
    ranked.sort(key=lambda r: (-r[0], r[1], r[2]))
    for rank, (r2, distance, pos, allele, old_count, new_count) in enumerate(ranked, 1):
        if 0 < new_count < len(modern):
            return dict(position=pos, distance_bp=distance, allele=allele, r2=r2, rank=rank,
                        candidates=len(ranked), historical_tag_af=old_count/n, modern_tag_af=new_count/len(modern))
    return None


def replay(source, output, cfg, pairs, contract):
    source, output = Path(source), Path(output)
    output.mkdir(parents=True, exist_ok=True)
    receipt = output / 'tag.json'
    if receipt.exists():
        record = json.loads(receipt.read_text())
        assert record['contract'] == contract and record['source'] == str(source)
        if record['tag'] is not None:
            assert base.digest(output/'scores.npz') == record['scores_sha256']
        return record

    archived = json.loads((source/'rejected.json').read_text())
    assert archived['rejection'] == 'fixation' and archived['final_frequency'] == 1
    for name, field in [('onset.trees', 'onset_sha256'), ('trajectory.tsv', 'trajectory_sha256')]:
        assert base.digest(source/name) == archived[field]
    trajectory = np.loadtxt(source/'trajectory.tsv', skiprows=1)
    past_index = int(np.flatnonzero((trajectory[:,1] > 0) & (trajectory[:,1] <= .8))[-1])
    tick = past_index + 1
    label = archived['seed_namespace'] + archived['unit']
    assert base.seed(label) == archived['seed']
    size = int(trajectory[past_index,2])
    n_past = min(size, 1000)
    old_indices = np.random.default_rng(base.seed(label+':tag-history')).choice(size, n_past, replace=False)
    remember = f"{tick} late() {{ sim.treeSeqRememberIndividuals(p0.individuals[c({','.join(map(str,sorted(old_indices)))} )]); }}"
    original = (source/'model.slim').read_text()
    stop = 'if (af==0.0 | af==1.0) sim.simulationFinished();'
    assert original.count(stop) == 1 and original.count(str(source/'forward.trees')) == 1
    model = original.replace(stop, '').replace(str(source/'forward.trees'), str(output/'forward.trees'))
    (output/'replay.slim').write_text(model+'\n'+remember+'\n')
    native_lib = str(Path(cfg['slim']).parent.parent/'lib')
    lib_paths = os.environ.get('LD_LIBRARY_PATH', '').split(':')
    env = dict(os.environ, LD_LIBRARY_PATH=':'.join(dict.fromkeys([native_lib, *filter(None,lib_paths)])))
    with (output/'slim.log').open('w') as stdout, (output/'slim.err').open('w') as stderr:
        subprocess.run([cfg['slim'], '-s', str(base.seed(label+':slim')), str(output/'replay.slim')],
                       stdout=stdout, stderr=stderr, env=env, check=True)
    rows = np.array([list(map(float,l.split()[1:])) for l in (output/'slim.log').read_text().splitlines()
                     if l.startswith('TRAJECTORY ')])
    np.testing.assert_array_equal(rows[:len(trajectory), 4], trajectory[:,4])
    np.testing.assert_array_equal(rows[:len(trajectory), 2], trajectory[:,2])
    assert np.all(rows[len(trajectory)-1:,4] == 2*rows[len(trajectory)-1:,2])
    assert len(rows) == 1+archived['ancestry_start_years']//cfg['generation_time']
    forward = tskit.load(output/'forward.trees')
    old = forward.samples(time=len(rows)-tick)
    assert len(old) == 2*n_past
    alive = pyslim.individuals_alive_at(forward, 0)
    individuals = np.random.default_rng(base.seed(label+':sample')).choice(alive, cfg['diploids'], replace=False)
    new = np.concatenate([forward.individual(int(i)).nodes for i in individuals])
    joint = forward.simplify(np.r_[old,new])
    old, new = joint.samples()[:len(old)], joint.samples()[len(old):]
    center = cfg['length']//2
    focal, = joint.variants()
    assert focal.site.position == center
    focal_carriers = (focal.genotypes[:len(old)] > 0)
    assert np.all(focal.genotypes[len(old):] > 0)
    tables = joint.dump_tables()
    tables.sites.clear()
    tables.mutations.clear()
    tables.time_units = 'generations'
    observed = msprime.sim_mutations(tables.tree_sequence(), rate=cfg['mutation_rate'],
                                     random_seed=base.seed(label+':mutations'), model=msprime.JC69())
    tag = choose_tag(observed, old, new, focal_carriers, center)
    record = dict(contract=contract, source=str(source), onset_sha256=archived['onset_sha256'],
                  source_trajectory_sha256=archived['trajectory_sha256'], selection=archived['selection'],
                  origin=archived['origin'], regime=archived['regime'], fixation_kya=archived['stop_years_ago']/1000,
                  snapshot_kya=trajectory[past_index,0]/1000, historical_population_af=trajectory[past_index,1],
                  historical_sample_af=float(focal_carriers.mean()), historical_diploids=n_past, modern_diploids=len(individuals),
                  seed=archived['seed'], tag=tag)
    if tag is not None:
        observed.simplify(new).dump(output/'observed.trees')
        pairfile = output/'pairs.tsv'
        if not pairfile.exists():
            pairfile.symlink_to(pairs)
        assert base.digest(pairfile)==base.digest(pairs)
        np.savetxt(output/'positions.txt', [tag['position']], fmt='%d')
        os.environ['LD_LIBRARY_PATH'] = env['LD_LIBRARY_PATH']
        record.update(base.score(cfg, output, output))
        with np.load(output/'scores.npz') as d:
            np.testing.assert_array_equal(d['valid_pairs'], 1000)
            np.testing.assert_array_equal(d['positions'], [tag['position']])
    arrays.save_ledger(receipt, record)
    return record
