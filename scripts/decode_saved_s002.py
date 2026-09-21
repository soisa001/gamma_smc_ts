"""Decode the 100 saved s=.002 trees; never simulate or recover ancestry by replay.

Only the selected allele has an observed archaic carrier label in this archive.
Other grid coordinates support all-pair recency, NOT spatial carrier inference.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import inspect
import json
import os
from pathlib import Path
import time
import traceback

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "1"
import numpy as np
import pandas as pd
import tskit

from gamma_smc_aou.decoder import run_within_decoder
from gamma_smc_aou.fresh_power import atomic_json, canonical_hash, digest, ordered_nodes
from gamma_smc_aou.joint_scan_profiles import extract_profiles
from gamma_smc_aou.pair_class_profiles import posterior_arrays

ROOT = Path("/mnt/d/phase2simselection/sim")
STUDY = ROOT / "eas_q02_h400"
OUT = ROOT / "eas_s002_saved_decoding"
REPO = Path(__file__).resolve().parents[1]


def checked(path, sha):
    if digest(path) != sha:
        raise ValueError(f"Corrupt source: {path}")
    return path


def valid_receipt(path, fingerprint):
    if not path.exists():
        return False
    old = json.loads(path.read_text())
    if old["fingerprint"] != fingerprint:
        raise ValueError(f"Changed computation or inputs: {path}")
    for name, spec in old["outputs"].items():
        output = path.parent / name
        if output.stat().st_size != spec["bytes"] or digest(output) != spec["sha256"]:
            raise ValueError(f"Corrupt cache: {output}")
    return True


def artifacts(directory, names):
    return {name: dict(bytes=(directory/name).stat().st_size, sha256=digest(directory/name)) for name in names}


def one(payload):
    row, cfg, identity_common = payload
    dest = OUT / "profiles" / f"onset50000/rep{int(row['replicate']):04d}"
    dest.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    try:
        tree_path = checked(Path(row["decoder_input_path"]), row["decoder_input_sha256"])
        carrier_path = checked(Path(row["carriers_path"]), row["carriers_sha256"])
        identity = dict(identity_common, tree=digest(tree_path), carriers=digest(carrier_path))
        fingerprint = canonical_hash(identity)
        if valid_receipt(dest/"complete.json", fingerprint):
            return dict(task_id=row["task_id"], status="cached")
        raw_identity = {k: v for k, v in identity.items() if k not in ("extraction", "carriers")}
        raw_fingerprint = canonical_hash(raw_identity)
        if not valid_receipt(dest/"posterior_complete.json", raw_fingerprint):
            result = run_within_decoder(
                REPO/"bin/gamma_smc", tree_path, dest/"native_mean.tsv", raw_output=dest/"posteriors.zst",
                scaled_mutation_rate=cfg["decoder_scaled_mutation_rate"],
                recombination_to_mutation_ratio=cfg["decoder_recombination_to_mutation_ratio"],
                mutation_rate=cfg["mutation_rate"], threshold_years=cfg["tmrca_cutoffs_years"],
                generation_time=cfg["generation_time_years"], output_at_stride=-1, output_at_hets=False,
                only_within=False, output_positions_file=STUDY/"positions.txt", pairs_file=STUDY/"pairs.tsv",
                recent_call="mean", threads=1, cache_size=cfg["decoder_cache_size"],
                pair_block=cfg["decoder_pair_block"], exp10=cfg["decoder_exp10"],
                backward_alignment=cfg["decoder_backward_alignment"], vcf_position_transform="one_based",
                extra_args=["--no_recent_probability"],
            )
            for channel in ("stdout", "stderr"):
                (dest/f"decoder.{channel}.log").write_text(result.pop(channel))
            atomic_json(dest/"posterior_complete.json", dict(fingerprint=raw_fingerprint, identity=raw_identity,
                decode=result, outputs=artifacts(dest, ("posteriors.zst", "posteriors.zst.meta", "native_mean.tsv"))))
        pairs = np.loadtxt(STUDY/"pairs.tsv", dtype=int)
        grid = np.loadtxt(STUDY/"positions.txt", dtype=int)
        carriers = np.load(carrier_path).astype(bool)[None, :]
        meta = json.loads((dest/"posteriors.zst.meta").read_text())
        np.testing.assert_array_equal(meta["pairs"], pairs)
        np.testing.assert_array_equal(meta["output_positions"], grid)
        assert meta["chunk_size"] == 8 and carriers.shape == (1, 400)
        ts = tskit.load(tree_path)
        focal = cfg["focal_position_bp"]
        variants = list(ts.variants(samples=ordered_nodes(ts), left=focal, right=focal+1))
        assert len(variants) == 1 and variants[0].site.position == focal
        np.testing.assert_array_equal(variants[0].genotypes == 1, carriers[0])
        alpha, beta = posterior_arrays(dest/"posteriors.zst", len(pairs), len(grid))
        counts, denominators, markers = extract_profiles(ts, pairs, grid, np.array([focal]), carriers, alpha, beta, grid, cfg)
        native = pd.read_csv(dest/"native_mean.tsv", sep="\t")
        native_counts = np.rint(native[[f"frac_recent_{t}" for t in cfg["tmrca_cutoffs_years"]]].to_numpy()*len(pairs)).astype(int)
        np.testing.assert_array_equal(counts["grid"][1, :, 3], native_counts)
        assert carriers.mean() == row["sample_af"]
        np.savez_compressed(dest/"features.npz", grid_counts=counts["grid"], grid_n=denominators["grid"],
            grid=grid, sites=np.array([focal]), grid_markers=markers, site_af=carriers.mean(axis=1),
            sources=np.array(["truth", "decoded"]), classes=np.array(["ref/ref", "alt/ref", "alt/alt", "all"]),
            cutoffs=np.array(cfg["tmrca_cutoffs_years"]), archaic_labels_scope=np.array("focal_only"),
            carrier_labels_known=grid == focal)
        # Original archive is immutable; check it again after decoding/extraction.
        checked(tree_path, row["decoder_input_sha256"])
        checked(carrier_path, row["carriers_sha256"])
        atomic_json(dest/"complete.json", dict(identity=identity, fingerprint=fingerprint,
            task_id=row["task_id"], s=.002, seconds=time.monotonic()-started,
            scope="All-pair recency across the 10 Mb grid; archaic carrier classes known only at the selected allele",
            verification=dict(class_partition=True, native_decoded_counts_exact=True,
                saved_allele_genotypes_exact=True, source_hashes_unchanged=True),
            outputs=artifacts(dest, ("features.npz",))))
        return dict(task_id=row["task_id"], status="complete", seconds=time.monotonic()-started)
    except Exception:
        error = traceback.format_exc()
        atomic_json(dest/"failed.json", dict(error=error))
        return dict(task_id=row["task_id"], status="failed", error=error)


def main(workers, limit):
    OUT.mkdir(parents=True, exist_ok=True)
    import fcntl
    with (OUT/"run.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        cfg = json.loads((STUDY/"manifest.json").read_text())["config"]
        frame = pd.read_csv(STUDY/"simulation_inventory.csv", float_precision="round_trip")
        frame = frame[(frame["mode"] == "selected") & (frame.s == .002)].sort_values("replicate")
        assert len(frame) == 100
        if limit:
            frame = frame.head(limit)
        identity = dict(schema="saved-focal-only-profiles/v1", decoder=digest(REPO/"bin/gamma_smc"),
            pairs=digest(STUDY/"pairs.tsv"), grid=digest(STUDY/"positions.txt"),
            parameters={k: v for k, v in cfg.items() if k.startswith("decoder_") or k in
                ("mutation_rate", "generation_time_years", "tmrca_cutoffs_years", "recent_call")},
            extraction=canonical_hash([inspect.getsource(extract_profiles), inspect.getsource(posterior_arrays), inspect.getsource(one)]))
        atomic_json(OUT/"manifest.json", dict(config=cfg, selected_s=.002, workers=workers,
            identity=identity, source_sha256=digest(Path(__file__)), inventory_sha256=digest(STUDY/"simulation_inventory.csv"),
            simulation_replays=False, new_simulations=False, spatial_carrier_scores_available=False,
            items=frame.to_dict("records")))
        results = []
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(one, (row, cfg, identity)) for row in frame.to_dict("records")]
            for future in as_completed(futures):
                result = future.result()
                results.append(result)
                print(json.dumps(result), flush=True)
                atomic_json(OUT/"status.json", dict(state="running", target=len(frame),
                    completed=sum(r["status"] != "failed" for r in results),
                    failures=[r for r in results if r["status"] == "failed"]))
        failures = [r for r in results if r["status"] == "failed"]
        atomic_json(OUT/"status.json", dict(state="failed" if failures else "complete", target=len(frame),
            completed=len(results)-len(failures), failures=failures))
        if failures:
            raise RuntimeError("Saved decoding failed; see status.json")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--workers", type=int, default=20)
    p.add_argument("--limit", type=int)
    a = p.parse_args()
    main(a.workers, a.limit)
