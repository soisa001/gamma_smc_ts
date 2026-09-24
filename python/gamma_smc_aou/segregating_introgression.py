"""Select a natural archaic-derived EAS polymorphism at the treatment onset.

Neutral prehistory uses unscaled DTWF through the archaic split, then Hudson
ancestry deeper in the shared human branch. SLiM runs the subsequent unscaled
WF generations. Mutations present at onset are retained, not redrawn afterwards.
"""
from __future__ import annotations

import json
from pathlib import Path
import subprocess
import time
import warnings

import msprime
import numpy as np
import pandas as pd
import pyslim
import tskit

from . import fresh_power as fp
from .selection import _sample_diploids

SPEC = "nearest-archaic-segregating-at-onset/v1"


def tasks(cfg):
    result = []
    for rep in range(max(cfg['null_replicates'], cfg['target_replicates'])):
        for onset in cfg['selection_onsets_years']:
            family = f'I{onset//1000}'
            common = dict(family=family, origin='introgressed', introduction_years=cfg['pulse_years'],
                          ascertainment_years=onset, replicate=rep, scaling_factor=1)
            if rep < cfg['null_replicates']:
                result.append(dict(common, id=f'{family}/null/rep{rep:04d}', role='null', s=0., onset_years=0))
            if rep < cfg['target_replicates']:
                result.append(dict(common, id=f'{family}/neutral_target/rep{rep:04d}', role='neutral_target', s=0., onset_years=0))
                for s in cfg['selection_coefficients']:
                    arm = f'onset{onset}_s{s:.3f}'.replace('.', 'p')
                    result.append(dict(common, id=f'{family}/{arm}/rep{rep:04d}', role='selected', s=s, onset_years=onset))
    return result


def stage_seed(seed, stage):
    return int(fp.canonical_hash([seed, stage])[:16], 16) % (2**31 - 2) + 1


def history(cfg):
    artifact = fp.load_phlash_eas_npz(fp.REPO / cfg['phlash_resource'], expected_sha256=cfg['phlash_sha256'])
    return fp.build_eas_demography_models(artifact)['median']


def population_sizes(cfg, times):
    h = history(cfg)
    indices = np.searchsorted(h.time_generations, times, side='right') - 1
    return np.maximum(2, np.rint(h.ne[indices]).astype(int))


def demography(cfg):
    h = history(cfg)
    dem = msprime.Demography()
    dem.add_population(name='EAS', initial_size=float(h.ne[0]))
    dem.add_population(name='Neanderthal', initial_size=cfg['archaic_effective_size'])
    for t, ne in zip(h.time_generations[1:], h.ne[1:]):
        dem.add_population_parameters_change(time=float(t), initial_size=float(ne), population='EAS')
    dem.add_mass_migration(time=cfg['pulse_years']/cfg['generation_time_years'], source='EAS',
                           dest='Neanderthal', proportion=cfg['introgression_proportion'])
    dem.add_mass_migration(time=cfg['archaic_split_years']/cfg['generation_time_years'], source='Neanderthal',
                           dest='EAS', proportion=1)
    dem.sort_events()
    dem.validate()
    return dem


def mutation_population(ts, mutation, position, migrations_by_node):
    """Population on a mutation's branch at its age, including pulse crossings.

    A branch's lower node can be in EAS even when its older portion is archaic.
    Migration records must be read before simplification removes their nodes.
    """
    population = ts.node(mutation.node).population
    for migration in migrations_by_node.get(mutation.node, ()):
        if migration.left <= position < migration.right and migration.time <= mutation.time:
            population = migration.dest
    return population


def choose_focal(ts, cfg):
    """Nearest unique, biallelic, archaic-origin mutation segregating at onset."""
    low = cfg['focal_position_bp']
    high = ts.sequence_length - (cfg['scored_length_bp'] - cfg['focal_position_bp'])
    center = ts.sequence_length / 2
    pulse = cfg['pulse_years']/cfg['generation_time_years']
    split = cfg['archaic_split_years']/cfg['generation_time_years']
    migrations = {}
    for migration in sorted(ts.migrations(), key=lambda m: m.time):
        migrations.setdefault(migration.node, []).append(migration)
    choices = sorted((site for site in ts.sites() if low <= site.position <= high),
                     key=lambda site: (abs(site.position-center), site.position))
    inspected = 0
    for site in choices:
        inspected += 1
        if len(site.mutations) != 1:
            continue
        mutation = site.mutations[0]
        if not pulse <= mutation.time < split:
            continue
        if mutation_population(ts, mutation, site.position, migrations) != 1:
            continue
        tree = ts.at(site.position)
        carriers = tree.num_samples(mutation.node)
        if not 0 < carriers < ts.num_samples:
            continue
        return dict(position=int(site.position), site_id=site.id, mutation_id=mutation.id,
                    mutation_time_generations=float(mutation.time), mutation_population='Neanderthal',
                    onset_alt_copies=int(carriers), onset_total_copies=ts.num_samples,
                    onset_af=carriers/ts.num_samples, distance_bp=abs(site.position-center),
                    eligible_window=[low, high], sites_inspected=inspected,
                    rule='nearest archaic-origin biallelic EAS polymorphism; ties use lower coordinate')
    return None


def onset_ancestry(cfg, task, seed):
    age = int(task['ascertainment_years']/cfg['generation_time_years'])
    size = int(population_sizes(cfg, [age])[0])
    split = cfg['archaic_split_years']/cfg['generation_time_years']
    dem = demography(cfg)
    samples = [msprime.SampleSet(size, population='EAS', time=age, ploidy=2)]
    if task['ascertainment_years'] == cfg['pulse_years']:
        # Ancient samples activate after same-time demographic events in msprime.
        # Represent the post-pulse cohort explicitly, without shifting its date:
        # each of the 2N haplotypes draws its source with probability q.
        count = np.random.default_rng(stage_seed(seed,'pulse_allocation')).binomial(2*size, cfg['introgression_proportion'])
        samples = [msprime.SampleSet(2*size-count, population='EAS', time=age, ploidy=1)]
        if count:
            samples.append(msprime.SampleSet(count, population='Neanderthal', time=age, ploidy=1))
        dem.events = [e for e in dem.events if not (isinstance(e,msprime.MassMigration) and e.source=='EAS')]
    # All diploid genomes alive at onset are needed for a forward WF restart.
    return msprime.sim_ancestry(samples=samples,
        demography=dem, sequence_length=cfg['simulated_length_bp'], recombination_rate=cfg['recombination_rate'],
        model=[msprime.DiscreteTimeWrightFisher(duration=split-age+1), msprime.StandardCoalescent()],
        record_migrations=True, random_seed=stage_seed(seed, 'ancestry'))


def initial_state(cfg, task, seed, directory):
    age = int(task['ascertainment_years']/cfg['generation_time_years'])
    initial = onset_ancestry(cfg, task, seed)
    initial = msprime.sim_mutations(initial, rate=cfg['mutation_rate'],
        model=msprime.SLiMMutationModel(type=0, slim_generation=age+1), random_seed=stage_seed(seed, 'onset_mutations'))
    focal = choose_focal(initial, cfg)
    if focal is None:
        return None
    initial.dump(directory/'ascertainment.trees')
    fp.atomic_json(directory/'focal_choice.json', dict(focal, seed=seed, ascertainment_years=task['ascertainment_years']))
    tables = initial.dump_tables()
    # Selection eligibility was resolved using absolute ages and migration records.
    tables.migrations.clear()
    # All sampled haplotypes now belong to the admixed EAS census. Randomly
    # pair them into diploids; source labels above were only for prehistory.
    nodes = initial.samples()
    population = tables.nodes.population.copy()
    population[nodes] = 0
    tables.nodes.population = population
    tables.individuals.clear()
    individual = np.full(initial.num_nodes, tskit.NULL, dtype=np.int32)
    pairs = np.random.default_rng(stage_seed(seed,'onset_pairing')).permutation(nodes).reshape(-1,2)
    for pair in pairs:
        individual[pair] = tables.individuals.add_row()
    tables.nodes.individual = individual
    # SLiM also requires the two nodes of every individual to be adjacent.
    order = np.concatenate((pairs.ravel(), np.setdiff1d(np.arange(initial.num_nodes), nodes)))
    inverse = np.empty(initial.num_nodes, dtype=np.int32)
    inverse[order] = np.arange(initial.num_nodes)
    node_rows = list(tables.nodes)
    tables.nodes.clear()
    for old_id in order:
        tables.nodes.append(node_rows[int(old_id)])
    tables.edges.parent = inverse[tables.edges.parent]
    tables.edges.child = inverse[tables.edges.child]
    tables.mutations.node = inverse[tables.mutations.node]
    tables.sort()
    tables.nodes.time -= age
    tables.mutations.time -= age
    annotated = pyslim.annotate(tables.tree_sequence(), model_type='WF', tick=1, stage='late', annotate_mutations=False)
    tables = annotated.dump_tables()
    assert annotated.site(annotated.mutation(focal['mutation_id']).site).position == focal['position']
    rows = list(tables.mutations)
    tables.mutations.clear()
    for i, row in enumerate(rows):
        metadata = row.metadata
        if i == focal['mutation_id']:
            metadata['mutation_list'][0]['mutation_type'] = 2
            metadata['mutation_list'][0]['selection_coeff'] = task['s']
        tables.mutations.append(row.replace(metadata=metadata))
    # pyslim 1.1 writes SLiM 5 metadata; the neutral mutation type is m0.
    annotated = tables.tree_sequence()
    focal['slim_mutation_id'] = annotated.mutation(focal['mutation_id']).derived_state
    annotated.dump(directory/'onset.trees')
    fp.atomic_json(directory/'focal_choice.json', dict(focal, seed=seed, ascertainment_years=task['ascertainment_years']))
    return focal


def slim_script(cfg, task, directory, focal):
    age = int(task['ascertainment_years']/cfg['generation_time_years'])
    sizes = population_sizes(cfg, np.arange(age, -1, -1))
    schedule = []
    for index in range(1, len(sizes)):
        if sizes[index] != sizes[index-1]:
            schedule.append(f'{index+1} early() {{ p0.setSubpopulationSize({sizes[index]}); }}')
    quoted = lambda name: json.dumps(str((directory/name).resolve()))
    script = f'''// Unscaled WF: Q=1. Natural focal polymorphism chosen before selection.
initialize() {{
    initializeSLiMOptions(keepPedigrees=T);
    initializeTreeSeq();
    initializeMutationRate(0);
    initializeMutationType("m0", 0.5, "f", 0.0);
    initializeMutationType("m2", {cfg['dominance']}, "f", {task['s']});
    m2.convertToSubstitution = F;
    initializeGenomicElementType("g1", m0, 1.0);
    initializeGenomicElement(g1, 0, {cfg['simulated_length_bp']-1});
    initializeRecombinationRate({cfg['recombination_rate']});
}}
1 late() {{
    sim.readFromPopulationFile({quoted('onset.trees')});
    muts = sim.mutationsOfType(m2);
    if (size(muts) != 1) stop("Expected exactly one chosen focal mutation");
    if (muts.position != {focal['position']}) stop("Focal position mismatch");
    if (sum(p0.haplosomes.containsMutations(muts)) != {focal['onset_alt_copies']}) stop("Focal carrier mismatch");
    if (p0.individualCount != {sizes[0]}) stop("Initial population size mismatch");
    sim.recalculateFitness();
}}
{chr(10).join(schedule)}
1:{age+1} late() {{
    muts = sim.mutationsOfType(m2);
    af = 0.0;
    if (size(muts) == 1) af = sim.mutationFrequencies(p0, muts);
    catn("TRAJECTORY " + community.tick + " " + af + " " + p0.individualCount);
    if (af == 0.0) {{
        catn("FOCAL_LOST");
        sim.simulationFinished();
    }}
}}
{age+1} late() {{
    sim.treeSeqOutput({quoted('forward.trees')});
    catn("FORWARD_COMPLETE");
    sim.simulationFinished();
}}
'''
    (directory/'model.slim').write_text(script)
    return age


def simulate(cfg, task, directory, identity, slim):
    from . import origin_onset as oo
    if oo.read_receipt(directory, 'simulation.json', identity):
        record = oo.read_receipt(directory, 'simulation.json', identity)
        oo.validate_trees(directory, cfg, record)
        return 'cached'
    attempt_identity = directory/'attempt_identity.json'
    attempt_path = directory/'attempts.json'
    if attempt_identity.exists():
        if json.loads(attempt_identity.read_text()) != identity:
            raise ValueError('Incompatible incomplete simulation attempts')
    else:
        if attempt_path.exists():
            raise ValueError('Incomplete attempts lack an identity receipt')
        fp.atomic_json(attempt_identity, identity)
    attempts = json.loads(attempt_path.read_text()) if attempt_path.exists() else []
    for index, row in enumerate(attempts):
        assert row['attempt']==index and row['seed']==oo.seed_for(cfg,task,index)
    # A crash during postprocessing replays only the last accepted attempt.
    if attempts and attempts[-1]['outcome']=='accepted':
        attempts.pop()
    for attempt in range(len(attempts),cfg['maximum_sample_attempts']):
        seed = oo.seed_for(cfg, task, attempt)
        started = time.monotonic()
        print(json.dumps(dict(stage='neutral_prehistory', attempt=attempt, seed=seed)), flush=True)
        focal = initial_state(cfg, task, seed, directory)
        outcome = 'no_eligible_onset_allele'
        row = dict(attempt=attempt, seed=seed)
        if focal is not None:
            row['focal'] = focal
            print(json.dumps(dict(stage='forward_selection', attempt=attempt, focal=focal)), flush=True)
            age = slim_script(cfg, task, directory, focal)
            command = [slim, '-s', str(stage_seed(seed, 'forward')), str(directory/'model.slim')]
            with (directory/f'attempt{attempt:03d}.stdout.log').open('w') as stdout, (directory/f'attempt{attempt:03d}.stderr.log').open('w') as stderr:
                completed = subprocess.run(command, stdout=stdout, stderr=stderr)
            if completed.returncode:
                raise RuntimeError(f'SLiM failed: {command}; see attempt{attempt:03d}.stderr.log')
            lines = (directory/f'attempt{attempt:03d}.stdout.log').read_text().splitlines()
            if 'FOCAL_LOST' not in lines and 'FORWARD_COMPLETE' not in lines:
                raise RuntimeError('SLiM exited without a recognized biological outcome')
            trajectory = [line.split()[1:] for line in lines if line.startswith('TRAJECTORY ')]
            outcome = 'lost_after_onset'
            if 'FORWARD_COMPLETE' in lines:
                raw = tskit.load(directory/'forward.trees')
                sampled = _sample_diploids(raw, cfg['sample_diploids'], stage_seed(seed, 'sample'))
                observed = fp.exact_focal_variant(sampled, fp.ordered_nodes(sampled), focal['position'])
                outcome = 'not_observed_in_sample'
                if observed is not None:
                    outcome = 'accepted'
        row.update(outcome=outcome, seconds=time.monotonic()-started)
        attempts.append(row)
        fp.atomic_json(directory/'attempts.json', attempts)
        print(json.dumps(row), flush=True)
        if outcome == 'accepted':
            break
    else:
        raise RuntimeError('No observed focal allele within the recorded attempt limit')
    finalize(cfg,task,directory,identity,seed,attempts,focal)
    return 'simulated'


def finalize(cfg,task,directory,identity,seed,attempts,focal):
    """Rebuild derived samples/overlays deterministically from the accepted run."""
    from . import origin_onset as oo
    age = int(task['ascertainment_years']/cfg['generation_time_years'])
    attempt = attempts[-1]['attempt']
    lines = (directory/f'attempt{attempt:03d}.stdout.log').read_text().splitlines()
    assert 'FORWARD_COMPLETE' in lines
    trajectory = [line.split()[1:] for line in lines if line.startswith('TRAJECTORY ')]
    raw = tskit.load(directory/'forward.trees')
    sampled = _sample_diploids(raw,cfg['sample_diploids'],stage_seed(seed,'sample'))
    observed = fp.exact_focal_variant(sampled,fp.ordered_nodes(sampled),focal['position'])
    assert observed is not None
    trajectory = np.asarray(trajectory, dtype=float)
    assert len(trajectory) == age+1 and np.array_equal(trajectory[:,0], np.arange(1,age+2))
    assert np.allclose(trajectory[:,2], population_sizes(cfg, np.arange(age,-1,-1)))
    pd.DataFrame(dict(years_ago=(age+1-trajectory[:,0])*cfg['generation_time_years'],
                      population_af=trajectory[:,1], population_diploids=trajectory[:,2])).to_csv(directory/'trajectory.csv', index=False)
    # Existing mutations are one realized neutral history. Add mutations only
    # during the forward interval, and protect the chosen site's allele identity.
    pos = focal['position']
    positions = [0, pos, pos+1, cfg['simulated_length_bp']]
    rates = [cfg['mutation_rate'], 0, cfg['mutation_rate']]
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=msprime.TimeUnitsMismatchWarning)
        sampled = msprime.sim_mutations(sampled, rate=msprime.RateMap(position=positions, rate=rates),
            end_time=age, model=msprime.SLiMMutationModel(type=0, next_id=pyslim.next_slim_mutation_id(sampled), slim_generation=age+1),
            keep=True, random_seed=stage_seed(seed, 'forward_mutations'))
    sampled = pyslim.convert_alleles(pyslim.generate_nucleotides(sampled, seed=stage_seed(seed, 'nucleotides')))
    observed_final = fp.exact_focal_variant(sampled, fp.ordered_nodes(sampled), pos)
    assert observed_final is not None
    np.testing.assert_array_equal(observed['carriers'], observed_final['carriers'])
    cropped, offset = fp.crop_at_site(sampled, pos, cfg)
    for name, tree in [('simulation.trees',sampled), ('decoded_input.trees',cropped)]:
        temporary = directory/(name+'.tmp')
        tree.dump(temporary)
        temporary.replace(directory/name)
    np.save(directory/'focal_carriers.npy', observed_final['carriers'])
    record = dict(identity=identity, task=task, seed=seed, sample_af=observed_final['af'], attempts=attempts,
                  focal_position_original=pos, focal_choice=focal, crop_offset=offset, trajectory_available=True,
                  final_census_af=float(trajectory[-1,1]), conditional_on='segregating at onset; focal observation in present-day sample; fixation retained',
                  artifacts=oo.artifact_specs(directory, (*fp.SIMULATION_ARTIFACTS, 'trajectory.csv', 'model.slim', 'focal_choice.json', 'onset.trees', 'ascertainment.trees')))
    oo.validate_trees(directory, cfg, record)
    fp.atomic_json(directory/'simulation.json', record)
    return record
