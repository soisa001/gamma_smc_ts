"""Recover archaic-origin SNP truth without changing the completed simulations.

Neutral replay records migrations without extra unary nodes. Selected replay remembers every
archaic individual at the pulse, consuming no additional random draws. Original
trees, mutations, pairs and decoding inputs stay authoritative. Equality checks
must pass before recovered ancestry can label any marker.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import redirect_stdout, redirect_stderr
import inspect
import json
import os
from pathlib import Path
import time
import traceback

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_name] = "1"
import msprime
import numpy as np
import stdpopsim
from stdpopsim import slim_engine
import tskit

from .fresh_power import atomic_json, canonical_hash, digest, model, ordered_nodes, selected_events, SITE_ID
from .posterior_replay import checked_source
from .recent_call_comparison import inventory
from .run7_models import scoped_slim_patch


def young_edges(ts, scale, ceiling):
    """Canonical SLiM pedigree edges below the archaic split, over all bases."""
    ids = np.array([(n.metadata or {}).get("slim_id", -1) for n in ts.nodes()], dtype=np.int64)
    e = ts.tables.edges
    keep = ts.nodes_time[e.parent] * scale < ceiling
    if np.any(ids[e.parent[keep]] < 0) or np.any(ids[e.child[keep]] < 0):
        raise ValueError("Missing pedigree identity in the forward simulation")
    return sorted(zip(e.left[keep], e.right[keep], ids[e.parent[keep]], ids[e.child[keep]]))


def replay_ancestry(item, slim):
    cfg = item["cfg"]
    original = tskit.load(checked_source(item, "simulation.trees"))
    receipt = json.loads((Path(item["directory"]) / "simulated.json").read_text())
    pulse = cfg["pulse_years"] / cfg["generation_time_years"]
    split = cfg["archaic_split_years"] / cfg["generation_time_years"]
    if item["mode"] == "neutral":
        ts = msprime.sim_ancestry(
            samples={"EAS": cfg["sample_diploids"]}, demography=model(cfg).model,
            sequence_length=cfg["simulated_length_bp"], recombination_rate=cfg["recombination_rate"],
            ploidy=2, random_seed=receipt["seed"], record_migrations=True,
        )
        tables = ts.dump_tables()
        tables.migrations.clear()
        collapsed = tables.tree_sequence().simplify()
        reference = original.simplify()
        for key in ("nodes", "edges", "individuals", "populations"):
            if getattr(collapsed.tables, key) != getattr(reference.tables, key):
                raise ValueError(f"Migration replay changed neutral {key}")
        samples = ordered_nodes(ts)
        ancestors = np.array([(m.left,m.right,m.node) for m in ts.migrations()
                              if m.time == pulse and m.dest == 1],dtype=float).reshape(-1,3)
        # Migration records can name nodes absent in a marginal tree. Never
        # silently interpret missing descendants as absence of introgression.
        for left,right,node in ancestors:
            tree = ts.at(left)
            while tree.interval.left < right:
                if tree.num_samples(int(node)) == 0:
                    raise ValueError("Orphan migration record: ancestry not recoverable from this replay")
                if not tree.next(): break
        evidence = "All simplified nodes, edges, individuals and populations equal original; every pulse migration has represented descendants over its entire interval"
    else:
        contig = stdpopsim.get_species("HomSap").get_contig(
            length=cfg["simulated_length_bp"], mutation_rate=cfg["mutation_rate"],
            recombination_rate=cfg["recombination_rate"],
        )
        contig.add_single_site(id=SITE_ID, coordinate=cfg["simulated_length_bp"] // 2)
        main = slim_engine._slim_main
        anchor = "    // Sample individuals."
        if main.count(anchor) != 1:
            raise ValueError("SLiM callback insertion contract changed")
        inserted = (
            f"    arc_g = time_to_tick({cfg['pulse_years']});\n"
            '    community.registerLateEvent(NULL, "{sim.treeSeqRememberIndividuals(p1.individuals);}", arc_g, arc_g);\n'
        )
        slim_engine._slim_main = main.replace(anchor, inserted + anchor)
        try:
            with scoped_slim_patch("selected", 0, patch_placement=True):
                ts = stdpopsim.get_engine("slim").simulate(
                    model(cfg), contig, {"EAS": cfg["sample_diploids"]},
                    seed=receipt["seed"], extended_events=selected_events(cfg, item["s"]),
                    slim_path=slim, slim_scaling_factor=cfg["slim_scaling_factor"],
                    slim_burn_in=cfg["slim_burn_in"], _recap_and_rescale=False,
                )
        finally:
            slim_engine._slim_main = main
        lookup = {(n.metadata or {}).get("slim_id"): n.id for n in ts.nodes()}
        samples = np.array([lookup[original.node(int(u)).metadata["slim_id"]]
                            for u in ordered_nodes(original)], dtype=int)
        collapsed = ts.simplify(samples=samples, filter_populations=False)
        if young_edges(collapsed, cfg["slim_scaling_factor"], split) != young_edges(original, 1, split):
            raise ValueError("Remembering archaic parents changed selected pedigree edges")
        old_nodes = {(n.metadata or {}).get("slim_id"): (n.time, n.population)
                     for n in original.nodes() if n.time < split}
        for n in collapsed.nodes():
            t = n.time * cfg["slim_scaling_factor"]
            if t < split and old_nodes.get(n.metadata.get("slim_id")) != (t, n.population):
                raise ValueError("Selected replay changed a node time or population")
        ancestors = np.where((ts.nodes_time * cfg["slim_scaling_factor"] == pulse)
                             & (ts.nodes_population == 1))[0]
        if not len(ancestors):
            raise ValueError("Archaic census is absent")
        evidence = "All pedigree edges, node times and populations below split equal original; sample pedigree order identical"
    return ts, samples, ancestors, receipt, evidence


def ancestry_mask(tree, ancestors, sample_lookup, position):
    mask = np.zeros(int(np.max(sample_lookup)) + 1, dtype=bool)
    if ancestors.ndim == 2:
        ancestors = ancestors[(ancestors[:,0] <= position) & (position < ancestors[:,1]),2]
    for ancestor in ancestors:
        # A retained ancestor can be isolated outside its transmitted interval.
        for node in tree.samples(int(ancestor)):
            index = sample_lookup[node]
            if index >= 0:
                mask[index] = True
    return mask


def label_sites(cropped, ancestry, samples, ancestors, offset, cfg):
    lookup = np.full(ancestry.num_nodes, -1, dtype=int)
    lookup[samples] = np.arange(len(samples))
    pulse = cfg["pulse_years"] / cfg["generation_time_years"]
    split = cfg["archaic_split_years"] / cfg["generation_time_years"]
    tree = ancestry.first()
    output_positions, site_ids, genotypes = [], [], []
    checked, mixed, nonbinary = 0, 0, 0
    last_index, mask = -1, None
    for variant in cropped.variants(samples=ordered_nodes(cropped)):
        site = variant.site
        if len(site.mutations) != 1 or len(variant.alleles) != 2:
            nonbinary += 1
            continue
        mutation = site.mutations[0]
        if not pulse <= mutation.time < split:
            continue
        alt = variant.genotypes == 1
        if not np.any(alt):
            continue
        position = site.position + offset
        while position >= tree.interval.right:
            if not tree.next():
                raise ValueError("Ancestry replay does not cover marker")
        if tree.index != last_index or ancestors.ndim == 2:
            mask = ancestry_mask(tree, ancestors, lookup, position)
            last_index = tree.index
        checked += 1
        if not np.all(mask[alt]):
            mixed += 1
            continue
        output_positions.append(int(site.position))
        site_ids.append(site.id)
        genotypes.append(alt)
    g = np.array(genotypes, dtype=bool).reshape(-1, len(samples))
    return np.asarray(output_positions, dtype=np.int64), np.asarray(site_ids, dtype=np.int32), g, dict(
        eligible_age_single_mutation=checked, rejected_nonarchaic=mixed,
        excluded_nonbinary_or_multimutation=nonbinary,
    )


def label_one(payload):
    item, out, slim, code_hash = payload
    dest = Path(out) / "labels" / item["key"]
    dest.mkdir(parents=True, exist_ok=True)
    identity = dict(simulation=item["receipt_sha256"], computation=code_hash, slim=digest(slim), schema="joint-archaic-labels/v1")
    fingerprint = canonical_hash(identity)
    start = time.monotonic()
    try:
        receipt_path = dest / "complete.json"
        if receipt_path.exists():
            receipt = json.loads(receipt_path.read_text())
            if receipt["fingerprint"] != fingerprint:
                raise ValueError("Label inputs or computation changed")
            for name, spec in receipt["outputs"].items():
                if (dest / name).stat().st_size != spec["bytes"] or digest(dest / name) != spec["sha256"]:
                    raise ValueError(f"Corrupt label artifact {name}")
            return dict(key=item["key"], status="cached", markers=receipt["markers"])
        with (dest / "execution.log").open("a", buffering=1) as log, redirect_stdout(log), redirect_stderr(log):
            ts, samples, ancestors, original, evidence = replay_ancestry(item, slim)
            cropped = tskit.load(checked_source(item, "decoded_input.trees"))
            positions, site_ids, genotypes, counts = label_sites(cropped, ts, samples, ancestors, original["crop_offset"], item["cfg"])
            if item["mode"] == "selected":
                i = np.where(positions == item["cfg"]["focal_position_bp"])[0]
                if len(i) != 1 or not np.array_equal(genotypes[i[0]], np.load(checked_source(item, "focal_carriers.npy"))):
                    raise ValueError("Recovered archaic SNP labels miss/change the selected allele")
            if not np.all(np.diff(positions) > 0):
                raise ValueError("Marker positions are not unique and sorted")
            ts.dump(dest / "ancestry_replay.trees")
            np.savez_compressed(dest / "markers.npz", positions=positions, site_ids=site_ids,
                                carriers=np.packbits(genotypes, axis=1), n_haplotypes=len(samples),
                                af=genotypes.mean(axis=1))
            np.savez_compressed(dest / "ancestry_mapping.npz", samples=samples, ancestors=ancestors,
                                crop_offset=original["crop_offset"])
            atomic_json(receipt_path, dict(identity=identity, fingerprint=fingerprint,
                markers=len(positions), fixed_markers=int(np.sum(genotypes.all(axis=1))),
                ascertainment="Single mutation, biallelic, observed ALT; pulse <= mutation age < archaic split; all ALT haplotypes inherit archaic pulse ancestry. Same rule in both cohorts; fixation retained.",
                verification=evidence, seed=original["seed"], counts=counts,
                seconds=time.monotonic()-start,
                outputs={name: dict(bytes=(dest/name).stat().st_size,sha256=digest(dest/name))
                         for name in ("markers.npz","ancestry_replay.trees","ancestry_mapping.npz")}))
        return dict(key=item["key"], status="complete", markers=len(positions), seconds=time.monotonic()-start)
    except Exception:
        error = traceback.format_exc()
        atomic_json(dest / "failed.json", dict(error=error))
        return dict(key=item["key"], status="failed", error=error)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("baseline", "later", "out", "decoder", "slim"):
        p.add_argument("--"+name, type=Path, required=True)
    p.add_argument("--workers", type=int, default=20)
    p.add_argument("--limit-per-arm", type=int)
    a = p.parse_args(argv)
    a.out.mkdir(parents=True, exist_ok=True)
    import fcntl
    with (a.out/"labels.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        items,cfg,studies = inventory(a.baseline,a.later,a.decoder)
        if a.limit_per_arm is not None:
            items = [i for i in items if i["replicate"] < a.limit_per_arm]
        # Do expensive selected replays first while neutral jobs fill remaining slots.
        items.sort(key=lambda i: (i["mode"] == "neutral", i["replicate"], i["onset_years"]))
        code_hash = canonical_hash([inspect.getsource(f) for f in (young_edges,replay_ancestry,ancestry_mask,label_sites)])
        atomic_json(a.out/"label_manifest.json",dict(config=cfg,studies=studies,samples=[i["key"] for i in items],
                    computation=code_hash,workers=a.workers,source_sha256=digest(Path(__file__))))
        done,failures = 0,[]
        with ProcessPoolExecutor(max_workers=a.workers) as pool:
            jobs=[pool.submit(label_one,(i,str(a.out),str(a.slim),code_hash)) for i in items]
            for future in as_completed(jobs):
                result=future.result()
                if result["status"]=="failed": failures.append(result)
                else: done+=1
                print(json.dumps(result),flush=True)
                atomic_json(a.out/"label_status.json",dict(completed=done,target=len(items),failed=failures,state="running"))
        atomic_json(a.out/"label_status.json",dict(completed=done,target=len(items),failed=failures,state="failed" if failures else "complete"))
        if failures: raise RuntimeError("Label recovery failed; inspect label_status.json")


if __name__ == "__main__":
    main()
