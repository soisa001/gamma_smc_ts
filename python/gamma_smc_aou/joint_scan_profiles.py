"""Truth and Gamma-SMC frac_recent_T at grid and observed archaic SNP positions."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import inspect
import json
import os
from pathlib import Path
import shutil
import time
import traceback

for _n in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_n] = "1"
import numpy as np
import pandas as pd
import tskit

from .decoder import run_within_decoder
from .fresh_power import atomic_json, canonical_hash, digest, ordered_nodes, read_profile
from .pair_class_profiles import posterior_arrays
from .posterior_replay import checked_source
from .recent_call_comparison import inventory
from .joint_truth import pair_ages


def nearest_markers(grid, sites, maximum_distance):
    if len(sites) == 0:
        return np.full(len(grid), -1, dtype=int)
    right = np.clip(np.searchsorted(sites, grid), 0, len(sites)-1)
    left = np.maximum(right-1, 0)
    index = np.where(abs(sites[left]-grid) <= abs(sites[right]-grid), left, right)
    return np.where(abs(sites[index]-grid) <= maximum_distance, index, -1)


def count_classes(recent, classes):
    """Recent is pairs x thresholds; all classes partition the SAME fixed pairs."""
    result = np.zeros((4, recent.shape[1]), dtype=np.uint16)
    for c in range(3):
        result[c] = recent[classes == c].sum(axis=0)
    result[3] = result[:3].sum(axis=0)
    return result


def extract_profiles(ts, pairs, grid, sites, carriers, alpha, beta, output_positions, cfg):
    if len(pairs) > np.iinfo(np.uint16).max:
        raise ValueError("Count storage dtype cannot hold the requested panel")
    grid_markers = nearest_markers(grid, sites, cfg["stride_bp"] // 2)
    classes_site = carriers[:, pairs].sum(axis=2).astype(np.int8)
    classes_grid = np.zeros((len(grid), len(pairs)), dtype=np.int8)
    valid = grid_markers >= 0
    classes_grid[valid] = classes_site[grid_markers[valid]]
    grid_out = np.searchsorted(output_positions, grid)
    site_out = np.searchsorted(output_positions, sites)
    outputs = {}
    denominators = {}
    for scheme, classes in (("grid",classes_grid),("sites",classes_site)):
        counts = np.array([(classes == c).sum(axis=1) for c in range(3)]).T
        denominators[scheme] = np.column_stack((counts,np.full(len(classes),len(pairs)))).astype(np.uint16)
        outputs[scheme] = np.zeros((2,len(classes),4,len(cfg["tmrca_cutoffs_years"])),dtype=np.uint16)
    two_ne = float(np.float32(cfg["decoder_scaled_mutation_rate"]))/(2*cfg["mutation_rate"])
    native_threshold = np.array([(y/cfg["generation_time_years"])/two_ne for y in cfg["tmrca_cutoffs_years"]],dtype=np.float32)
    truth_threshold = np.array(cfg["tmrca_cutoffs_years"])/cfg["generation_time_years"]
    node_pairs = ordered_nodes(ts)[pairs]
    node_times = ts.nodes_time
    tree = ts.first()
    last_tree, times = -1, None
    grid_lookup = {int(j):i for i,j in enumerate(grid_out)}
    site_lookup = {int(j):i for i,j in enumerate(site_out)}
    for j,pos in enumerate(output_positions):
        while pos >= tree.interval.right:
            if not tree.next(): raise ValueError("Truth sequence has incomplete coverage")
        if tree.index != last_tree:
            times = pair_ages(tree,node_pairs,node_times)
            last_tree = tree.index
        if not np.isfinite(times).all() or np.any(times<0):
            raise ValueError("Invalid true pairwise TMRCA")
        recent_truth = times[:,None] < truth_threshold
        recent_decoded = beta[j,:,None] * native_threshold >= np.clip(alpha[j,:,None],np.float32(2**-16),np.float32(2**16))
        for scheme,lookup,classes in (("grid",grid_lookup,classes_grid),("sites",site_lookup,classes_site)):
            if j in lookup:
                i=lookup[j]
                outputs[scheme][0,i] = count_classes(recent_truth,classes[i])
                outputs[scheme][1,i] = count_classes(recent_decoded,classes[i])
    for scheme in outputs:
        x = outputs[scheme]
        if not np.array_equal(x[:,:,:3].sum(axis=2),x[:,:,3]):
            raise ValueError("Class partition does not recover all-pair counts")
        if np.any(np.diff(x.astype(int),axis=-1)<0):
            raise ValueError("Nonmonotone recent counts")
        if np.any(x>denominators[scheme][None,:,:,None]):
            raise ValueError("Recent count exceeds its class denominator")
    return outputs, denominators, grid_markers


def profile_one(payload):
    item,out,decoder,code_hash = payload
    out=Path(out)
    dest=out/"profiles"/item["key"]
    dest.mkdir(parents=True,exist_ok=True)
    label=out/"labels"/item["key"]
    start=time.monotonic()
    try:
        label_receipt=json.loads((label/"complete.json").read_text())
        for name in ("markers.npz",):
            spec=label_receipt["outputs"][name]
            if (label/name).stat().st_size != spec["bytes"] or digest(label/name)!=spec["sha256"]:
                raise ValueError("Corrupt marker labels")
        cfg=item["cfg"]
        pairs_path=Path(item["root"])/"pairs.tsv"
        grid=np.loadtxt(Path(item["root"])/"positions.txt",dtype=int)
        pairs=np.loadtxt(pairs_path,dtype=int)
        markers=np.load(label/"markers.npz",allow_pickle=False)
        sites=markers["positions"]
        carriers=np.unpackbits(markers["carriers"],axis=1,count=int(markers["n_haplotypes"])).astype(bool)
        positions=np.union1d(grid,sites)
        identity=dict(simulation=item["receipt_sha256"],labels=digest(label/"complete.json"),
                      pairs=digest(pairs_path),positions=canonical_hash(positions.tolist()),
                      decoder=digest(decoder),computation=code_hash,schema="joint-scan-profiles/v1")
        fingerprint=canonical_hash(identity)
        receipt_path=dest/"complete.json"
        if receipt_path.exists():
            old=json.loads(receipt_path.read_text())
            if old["fingerprint"]!=fingerprint: raise ValueError("Profile identity changed")
            for name,spec in old["outputs"].items():
                if (dest/name).stat().st_size!=spec["bytes"] or digest(dest/name)!=spec["sha256"]:
                    raise ValueError("Corrupt profile output")
            return dict(key=item["key"],status="cached")
        tree_path=checked_source(item,"decoded_input.trees")
        np.savetxt(dest/"positions.txt",positions,fmt="%d")
        decode_id=canonical_hash({k:v for k,v in identity.items() if k not in ("computation","labels")})
        raw_receipt=dest/"posterior_complete.json"
        if raw_receipt.exists():
            raw=json.loads(raw_receipt.read_text())
            if raw["fingerprint"]!=decode_id: raise ValueError("Posterior identity changed")
            for name,spec in raw["outputs"].items():
                if (dest/name).stat().st_size!=spec["bytes"] or digest(dest/name)!=spec["sha256"]:
                    raise ValueError("Corrupt saved posterior")
        else:
            result=run_within_decoder(
                decoder,tree_path,dest/"native_mean.tsv",raw_output=dest/"posteriors.zst",
                scaled_mutation_rate=cfg["decoder_scaled_mutation_rate"],
                recombination_to_mutation_ratio=cfg["decoder_recombination_to_mutation_ratio"],
                mutation_rate=cfg["mutation_rate"],threshold_years=cfg["tmrca_cutoffs_years"],
                generation_time=cfg["generation_time_years"],output_at_stride=-1,output_at_hets=False,
                only_within=False,output_positions_file=dest/"positions.txt",pairs_file=pairs_path,
                recent_call="mean",threads=1,cache_size=cfg["decoder_cache_size"],
                pair_block=cfg["decoder_pair_block"],exp10=cfg["decoder_exp10"],
                backward_alignment=cfg["decoder_backward_alignment"],vcf_position_transform="one_based",
                extra_args=["--no_recent_probability"],
            )
            for channel in ("stdout","stderr"):
                (dest/f"decoder.{channel}.log").write_text(result.pop(channel))
            atomic_json(raw_receipt,dict(fingerprint=decode_id,decode=result,outputs={
                name:dict(bytes=(dest/name).stat().st_size,sha256=digest(dest/name))
                for name in ("posteriors.zst","posteriors.zst.meta","native_mean.tsv")}))
        meta=json.loads((dest/"posteriors.zst.meta").read_text())
        if not np.array_equal(meta["pairs"],pairs) or not np.array_equal(meta["output_positions"],positions) or meta["chunk_size"]!=8:
            raise ValueError("Posterior pair/position order mismatch")
        alpha,beta=posterior_arrays(dest/"posteriors.zst",len(pairs),len(positions))
        outputs,denominators,grid_markers=extract_profiles(tskit.load(tree_path),pairs,grid,sites,carriers,alpha,beta,positions,cfg)
        native=pd.read_csv(dest/"native_mean.tsv",sep="\t")
        native_counts=np.rint(native[[f"frac_recent_{t}" for t in cfg["tmrca_cutoffs_years"]]].to_numpy()*len(pairs)).astype(int)
        for scheme,pos in (("grid",grid),("sites",sites)):
            if not np.array_equal(outputs[scheme][1,:,3],native_counts[np.searchsorted(positions,pos)]):
                raise ValueError("Raw posterior mean rule differs from native integer counts")
        old_truth=read_profile(checked_source(item,"truth_frac_recent.tsv"),cfg).iloc[:,1:].to_numpy()
        if not np.array_equal(outputs["grid"][0,:,3],np.rint(old_truth*len(pairs))):
            raise ValueError("True all-pair grid differs from original")
        old_decoded=read_profile(checked_source(item,"frac_recent.tsv"),cfg).iloc[:,1:].to_numpy()
        difference=outputs["grid"][1,:,3].astype(int)-np.rint(old_decoded*len(pairs)).astype(int)
        reference_path=dest/"python_reference_features.npz"
        if reference_path.exists():
            with np.load(reference_path,allow_pickle=False) as reference:
                for scheme in ("grid","site"):
                    actual=outputs["sites" if scheme=="site" else "grid"]
                    if not np.array_equal(actual,reference[scheme+"_counts"]):
                        raise ValueError("Native exact-MRCA extraction differs from Python reference")
        np.savez_compressed(dest/"features.npz",grid_counts=outputs["grid"],site_counts=outputs["sites"],
                            grid_n=denominators["grid"],site_n=denominators["sites"],grid=grid,sites=sites,
                            grid_markers=grid_markers,site_af=markers["af"],
                            sources=np.array(["truth","decoded"]),classes=np.array(["ref/ref","alt/ref","alt/alt","all"]),
                            cutoffs=np.array(cfg["tmrca_cutoffs_years"]))
        atomic_json(receipt_path,dict(identity=identity,fingerprint=fingerprint,markers=len(sites),
            valid_grid_positions=int(np.sum(grid_markers>=0)),seconds=time.monotonic()-start,
            verification=dict(class_partition=True,original_truth_exact=True,native_decoded_counts_exact=True,
                              original_decoded_grid_max_abs_count_difference=int(abs(difference).max()),
                              original_decoded_grid_changed_cells=int(np.count_nonzero(difference))),
            outputs={"features.npz":dict(bytes=(dest/"features.npz").stat().st_size,sha256=digest(dest/"features.npz"))}))
        return dict(key=item["key"],status="complete",markers=len(sites),seconds=time.monotonic()-start,
                    grid_max_difference=int(abs(difference).max()))
    except Exception:
        error=traceback.format_exc()
        atomic_json(dest/"failed.json",dict(error=error))
        return dict(key=item["key"],status="failed",error=error)


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("baseline","later","out","decoder"):
        p.add_argument("--"+name,type=Path,required=True)
    p.add_argument("--workers",type=int,default=20)
    p.add_argument("--limit-per-arm",type=int)
    a=p.parse_args(argv)
    a.out.mkdir(parents=True,exist_ok=True)
    import fcntl
    with (a.out/"profiles.lock").open("w") as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        items,cfg,studies=inventory(a.baseline,a.later,a.decoder)
        if a.limit_per_arm is not None:items=[i for i in items if i["replicate"]<a.limit_per_arm]
        if shutil.disk_usage(a.out).used+500_000_000_000>cfg["storage_limit_bytes"]:
            raise ValueError("Insufficient storage within 3 TB allowance")
        code_hash=canonical_hash(dict(functions=[inspect.getsource(f) for f in (nearest_markers,count_classes,extract_profiles,posterior_arrays,pair_ages)],
                                     native_source=digest(Path(__file__).resolve().parents[2]/"cpp/joint_pair_tmrca.cpp"),
                                     native_binary=digest(Path(__file__).resolve().parents[2]/"bin/joint_pair_tmrca.so")))
        atomic_json(a.out/"profile_manifest.json",dict(config=cfg,studies=studies,items=items,
                    computation=code_hash,workers=a.workers,source_sha256=digest(Path(__file__))))
        done,failures=0,[]
        with ProcessPoolExecutor(max_workers=a.workers) as pool:
            jobs=[pool.submit(profile_one,(i,str(a.out),str(a.decoder),code_hash)) for i in items]
            for future in as_completed(jobs):
                result=future.result()
                if result["status"]=="failed":failures.append(result)
                else:done+=1
                print(json.dumps(result),flush=True)
                atomic_json(a.out/"profile_status.json",dict(completed=done,target=len(items),failed=failures,state="running"))
        atomic_json(a.out/"profile_status.json",dict(completed=done,target=len(items),failed=failures,state="failed" if failures else "complete"))
        if failures:raise RuntimeError("Joint profiles failed; inspect profile_status.json")


if __name__=="__main__":main()
