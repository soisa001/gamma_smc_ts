#!/usr/bin/env python3
"""Bounded, restartable no-introgression EAS truth/median/mass pilot (Linux)."""
from __future__ import annotations

import os
os.environ["MPLBACKEND"] = "Agg"
for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[key] = "1"
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import signal
import subprocess
import sys
import time

import msprime
import numpy as np
import pyslim
from scipy.special import gammainc
import tskit
import zstandard

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "python"))
os.environ["PYTHONPATH"] = str(REPO / "python") + os.pathsep + os.environ.get("PYTHONPATH", "")
from gamma_smc_aou.decoder import run_within_decoder
from gamma_smc_aou.selection import _sample_diploids
from gamma_smc_aou.fresh_power import ordered_nodes

VERSION = "time-resolved-pilot/v1"
BASE_SEED = 20380101
METHODS = ("truth", "median", "mass")


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def save(path, value):
    path = Path(path)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temp.replace(path)


def seed(label):
    return int(hashlib.sha256(f"{BASE_SEED}:time-resolved-pilot:{label}".encode()).hexdigest()[:12], 16) % (2**31-2) + 1


def history():
    path = REPO / "no_introgression_EAS_sim/resources/EAS.npz"
    assert digest(path) == "a5bdfd843629d48050cd56a84da4d70e356d0339017d3c1f67be646d727db928"
    with np.load(path) as a:
        return np.r_[0., a["time"].astype(float)], np.r_[np.median(a["bootstrap_ne"], axis=0)[0], np.median(a["bootstrap_ne"], axis=0)]


def simulate(cfg, label, selected, directory):
    start = time.monotonic()
    q, length = cfg["scaling"], cfg["length"]
    age, end = 1200 // q, 800 // q
    assert 1200 % q == 0 and 800 % q == 0
    times, ne = history()
    lookup = lambda t: float(ne[np.searchsorted(times, t, side="right")-1])
    size = round(lookup(1200) / q)
    dem = msprime.Demography()
    dem.add_population(name="EAS", initial_size=size)
    for t, n in zip(times[times > 1200], ne[times > 1200]):
        dem.add_population_parameters_change(time=(float(t)-1200)/q, initial_size=float(n)/q)
    initial = msprime.sim_ancestry(samples=size, demography=dem, sequence_length=length,
        recombination_rate=1e-8*q, random_seed=seed(label+":ancestry"),
        model=[msprime.DiscreteTimeWrightFisher(duration=100/q), msprime.StandardCoalescent()])
    center = length // 2
    tree = initial.at(center)
    candidates = [n for n in tree.nodes() if tree.parent(n) != tskit.NULL and .02 <= tree.num_samples(n)/initial.num_samples <= .10]
    if not candidates:
        raise ValueError("No eligible 2-10% standing branch at prespecified center; not silently retried")
    focal = min(candidates, key=lambda n: (abs(tree.num_samples(n)/initial.num_samples-.05), n))
    af0 = tree.num_samples(focal)/initial.num_samples
    tables = initial.dump_tables()
    site = tables.sites.add_row(center, ancestral_state="0")
    mutation_age = (tree.time(focal)+tree.time(tree.parent(focal)))/2
    tables.mutations.add_row(site=site, node=focal, derived_state="1", time=mutation_age)
    initial = pyslim.annotate(tables.tree_sequence(), model_type="WF", tick=1, stage="late")
    tables = initial.dump_tables()
    m = tables.mutations[0]
    metadata = m.metadata
    metadata["mutation_list"][0]["mutation_type"] = 2
    metadata["mutation_list"][0]["selection_coeff"] = cfg["selection"]*q if selected else 0.
    tables.mutations.clear()
    tables.mutations.append(m.replace(metadata=metadata))
    initial = tables.tree_sequence()
    initial.dump(directory / "onset.trees")
    sizes = [max(2,round(lookup((age-i)*q)/q)) for i in range(age+1)]
    schedules = [f"{i+1} early() {{ p0.setSubpopulationSize({n}); }}" for i,n in enumerate(sizes) if i and n != sizes[i-1]]
    # Reproduction at ticks 2 through 1+(1200-800)/Q experiences selection.
    stop_tick = age-end+2
    selected_s = cfg["selection"]*q if selected else 0.
    quote = lambda name: json.dumps(str(directory / name))
    script = f'''initialize() {{
 initializeSLiMOptions(keepPedigrees=T);
 initializeTreeSeq(); initializeMutationRate(0);
 initializeMutationType("m0", 0.5, "f", 0.0);
 initializeMutationType("m2", 0.5, "f", {selected_s});
 m2.convertToSubstitution=F;
 initializeGenomicElementType("g1",m0,1.0);
 initializeGenomicElement(g1,0,{length-1});
 initializeRecombinationRate({1e-8*q});
}}
1 late() {{
 sim.readFromPopulationFile({quote('onset.trees')});
 muts=sim.mutationsOfType(m2);
 if (size(muts)!=1) stop("missing focal mutation");
 if (sum(p0.haplosomes.containsMutations(muts))!={tree.num_samples(focal)}) stop("carrier mismatch");
 sim.recalculateFitness();
}}
{chr(10).join(schedules)}
{stop_tick} early() {{
 muts=sim.mutationsOfType(m2);
 if (size(muts)) muts.setSelectionCoeff(0.0);
 sim.recalculateFitness();
}}
1:{age+1} late() {{
 muts=sim.mutationsOfType(m2); af=0.0;
 if (size(muts)) af=sim.mutationFrequencies(p0,muts);
 catn("TRAJECTORY " + community.tick + " " + af + " " + p0.individualCount);
}}
{age+1} late() {{ sim.treeSeqOutput({quote('forward.trees')}); sim.simulationFinished(); }}
'''
    (directory / "model.slim").write_text(script)
    result = subprocess.run([cfg["slim"], "-s", str(seed(label+":slim")), str(directory/"model.slim")], capture_output=True, text=True, check=True)
    (directory/"slim.stdout.log").write_text(result.stdout)
    (directory/"slim.stderr.log").write_text(result.stderr)
    trajectory = np.array([list(map(float,line.split()[1:])) for line in result.stdout.splitlines() if line.startswith("TRAJECTORY ")])
    assert len(trajectory)==age+1
    np.testing.assert_array_equal(trajectory[:,0],np.arange(1,age+2))
    np.testing.assert_array_equal(trajectory[:,2],sizes)
    np.testing.assert_allclose(trajectory[0,1],af0,atol=1e-5)
    np.savetxt(directory/"trajectory.tsv",np.column_stack(((age+1-trajectory[:,0])*q*25,trajectory[:,1:])),delimiter="\t",header="years_ago\taf\tdiploids",comments="")
    raw = tskit.load(directory/"forward.trees")
    sampled = _sample_diploids(raw,cfg["diploids"],seed(label+":sample"))
    # Preserve genealogy but remove the deliberately ascertained focal marker;
    # independent neutral sequence variation is the decoder's observation.
    tables=sampled.dump_tables()
    tables.sites.clear(); tables.mutations.clear()
    tables.nodes.time = tables.nodes.time*q
    tables.time_units="generations"
    sampled=tables.tree_sequence()
    observed=msprime.sim_mutations(sampled,rate=1.29e-8,random_seed=seed(label+":mutations"),model=msprime.JC69())
    observed.dump(directory/"observed.trees")
    return dict(initial_frequency=af0,final_frequency=float(trajectory[-1,1]),
        mutation_age_years=(mutation_age*q+1200)*25,simulation_seconds=time.monotonic()-start,
        selection_start_years=30000,selection_end_years=20000,selected_reproduction_ticks=[2,stop_tick-1],
        sites=observed.num_sites,trees=observed.num_trees)


def read_exact(reader, size):
    data=bytearray()
    while len(data)<size:
        block=reader.read(size-len(data))
        if not block: raise EOFError("Truncated posterior stream")
        data.extend(block)
    return data


def score(cfg,directory,root):
    start=time.monotonic()
    positions=np.loadtxt(root/"positions.txt",dtype=int,ndmin=1)
    pairs=np.loadtxt(root/"pairs.tsv",dtype=int,ndmin=2)
    ts=tskit.load(directory/"observed.trees")
    nodes=ordered_nodes(ts)
    cutoffs=np.arange(1000,50001,1000,dtype=float)
    truth=np.zeros((len(positions),len(cutoffs)))
    for i,x in enumerate(positions):
        tree=ts.at(x)
        ages=np.array([tree.tmrca(int(nodes[a]),int(nodes[b]))*25 for a,b in pairs])
        truth[i]=np.mean(ages[:,None]<cutoffs,axis=0)
    raw=directory/"posterior.zst"
    result=run_within_decoder(cfg["decoder"],directory/"observed.trees",directory/"native.tsv",
        raw_output=raw,mutation_rate=1.29e-8,recombination_to_mutation_ratio=1e-8/1.29e-8,
        threshold_years=[10000,50000],generation_time=25,only_within=False,pairs_file=root/"pairs.tsv",
        output_positions_file=root/"positions.txt",output_at_stride=-1,threads=1,cache_size=1000,
        vcf_position_transform="one_based",extra_args=["--no_recent_probability"])
    save(directory/"decoder.json",result)
    meta=json.loads(Path(str(raw)+".meta").read_text())
    np.testing.assert_array_equal(meta["pairs"],pairs)
    np.testing.assert_array_equal(meta["output_positions"],positions)
    assert meta["num_pairs"]==len(pairs) and meta["sequence_length"]==len(positions)
    mass=np.zeros_like(truth); median=np.zeros_like(truth)
    counts=np.zeros(len(positions),dtype=int)
    chunk=meta["chunk_size"]
    nbytes=len(positions)*chunk*4
    scaled=cutoffs/25/(meta["scaled_mutation_rate"]/(2*1.29e-8))
    with raw.open("rb") as handle,zstandard.ZstdDecompressor().stream_reader(handle) as reader:
        for first in range(0,len(pairs),chunk):
            n=min(chunk,len(pairs)-first)
            alpha=np.frombuffer(read_exact(reader,nbytes),dtype=np.float32).reshape(len(positions),chunk)[:,:n].astype(float)
            beta=np.frombuffer(read_exact(reader,nbytes),dtype=np.float32).reshape(len(positions),chunk)[:,:n].astype(float)
            valid=np.isfinite(alpha)&np.isfinite(beta)&(alpha>0)&(beta>0)
            if not valid.all(): raise ValueError("Invalid posterior; no changing pair denominator allowed")
            counts+=valid.sum(axis=1)
            for k,t in enumerate(scaled):
                probabilities=gammainc(alpha,beta*t)
                mass[:,k]+=probabilities.sum(axis=1)
                median[:,k]+=(probabilities>.5).sum(axis=1)
        if reader.read(1): raise ValueError("Trailing posterior bytes")
    np.testing.assert_array_equal(counts,np.full(len(positions),len(pairs)))
    cdfs=np.stack((truth,median/counts[:,None],mass/counts[:,None]))
    assert np.isfinite(cdfs).all() and cdfs.min()>=0 and cdfs.max()<=1+1e-12
    assert np.diff(cdfs,axis=-1).min()>=-1e-12
    np.savez_compressed(directory/"scores.npz",cdfs=cdfs,positions=positions,cutoffs=cutoffs,valid_pairs=counts)
    return dict(score_seconds=time.monotonic()-start,posterior_sha256=digest(raw),scores_sha256=digest(directory/"scores.npz"))


def worker(args):
    root=Path(args.out).resolve(); cfg=json.loads((root/"manifest.json").read_text())["config"]
    directory=root/args.unit;directory.mkdir(exist_ok=True)
    receipt=directory/"complete.json"
    if receipt.exists():
        record=json.loads(receipt.read_text())
        assert record["scores_sha256"]==digest(directory/"scores.npz")
        return
    selected=args.unit.startswith("selected")
    label=f"q{cfg['scaling']}:{args.unit}"
    record=dict(unit=args.unit,selected=selected,seed=seed(label),config_sha256=digest(root/"manifest.json"))
    sim_receipt=directory/"simulation.json"
    if sim_receipt.exists():
        previous=json.loads(sim_receipt.read_text())
        assert previous["config_sha256"]==record["config_sha256"]
        assert previous["observed_sha256"]==digest(directory/"observed.trees")
        record=previous
    else:
        record.update(simulate(cfg,label,selected,directory))
        record["observed_sha256"]=digest(directory/"observed.trees")
        save(sim_receipt,record)
    record.update(score(cfg,directory,root))
    save(receipt,record)
    print(json.dumps(record),flush=True)


def launch(args):
    root=Path(args.out).resolve();root.mkdir(parents=True,exist_ok=True)
    cfg=dict(schema=VERSION,scaling=args.scaling,length=args.length,diploids=args.diploids,pairs=args.pairs,
        selection=args.selection,slim=str(Path(args.slim).resolve()),decoder=str(Path(args.decoder).resolve()),
        base_seed=BASE_SEED,pair_seed=1729,stride=10000,generation_time=25,cutoffs=list(range(1000,50001,1000)),
        source_sha256=digest(__file__),slim_sha256=digest(args.slim),decoder_sha256=digest(args.decoder))
    manifest=root/"manifest.json"
    if manifest.exists():
        old=json.loads(manifest.read_text())["config"]
        assert old==cfg, "Different parameters or implementation; use a new output directory"
    else:
        save(manifest,dict(config=cfg,python=sys.version,versions={name:__import__(name).__version__ for name in ("numpy","scipy","msprime","pyslim","tskit")},created_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ",time.gmtime())))
    possible=np.column_stack(np.triu_indices(2*args.diploids,1))
    pairs=possible[np.sort(np.random.default_rng(1729).choice(len(possible),args.pairs,replace=False))]
    np.savetxt(root/"pairs.tsv",pairs,fmt="%d",delimiter="\t")
    np.savetxt(root/"positions.txt",np.arange(0,args.length,10000),fmt="%d")
    # Round-robin cohorts; never select successes by effect size or completion speed.
    targets={"development":args.development,"calibration":args.calibration,"neutral_test":args.test,"selected":args.selected}
    tasks=[f"{cohort}_{i:04d}" for i in range(max(targets.values())) for cohort,n in targets.items() if i<n]
    deadline=time.time()+args.seconds
    save(root/"schedule.json",dict(tasks=tasks,workers=args.workers,deadline_epoch=deadline,timeout_per_unit=args.unit_timeout))
    def execute(unit):
        if (root/unit/"complete.json").exists(): return dict(unit=unit,state="cached")
        remaining=deadline-time.time()
        if remaining<10: return dict(unit=unit,state="not_started_deadline")
        directory=root/unit;directory.mkdir(exist_ok=True)
        with (directory/"worker.stdout.log").open("w") as out,(directory/"worker.stderr.log").open("w") as err:
            proc=subprocess.Popen([sys.executable,__file__,"worker","--out",str(root),"--unit",unit],stdout=out,stderr=err,start_new_session=True)
            try:
                code=proc.wait(timeout=min(args.unit_timeout,remaining))
                result=dict(unit=unit,state="complete" if code==0 else "failed",returncode=code)
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid,signal.SIGKILL);proc.wait()
                result=dict(unit=unit,state="timeout")
        print(json.dumps(result),flush=True)
        return result
    results=[]
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for future in as_completed([pool.submit(execute,t) for t in tasks]):
            results.append(future.result());save(root/"status.json",results)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("phase",choices=["run","worker"]);p.add_argument("--out",required=True);p.add_argument("--unit")
    p.add_argument("--slim");p.add_argument("--decoder");p.add_argument("--scaling",type=int,default=1)
    p.add_argument("--length",type=int,default=1_000_000);p.add_argument("--diploids",type=int,default=100)
    p.add_argument("--pairs",type=int,default=1000);p.add_argument("--selection",type=float,default=.01)
    p.add_argument("--development",type=int,default=40);p.add_argument("--calibration",type=int,default=99)
    p.add_argument("--test",type=int,default=40);p.add_argument("--selected",type=int,default=40)
    p.add_argument("--workers",type=int,default=4);p.add_argument("--seconds",type=int,default=1200)
    p.add_argument("--unit-timeout",type=int,default=180)
    args=p.parse_args()
    worker(args) if args.phase=="worker" else launch(args)


if __name__=="__main__": main()
