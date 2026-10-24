"""Neutral SV-age pilot with independent SNP observations and region-held-out evaluation.

Focal SVs are single-origin binary event surrogates, NOT structural read simulations.
Candidate events and SNP observations are independent overlays of the same genealogy.
Truth is used only for ascertainment, oracle comparisons, training and evaluation.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import os
import platform
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import reduce
from pathlib import Path

os.environ["MPLBACKEND"] = "Agg"
for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_name] = "1"

import msprime
import numpy as np
import pandas as pd
import scipy
import tskit
import zstandard
from scipy.spatial.distance import cdist
from scipy.special import gammainc, gammaincc, logsumexp

from .decoder import run_within_decoder
from .segregating_introgression import demography, mutation_population

REPO = Path(__file__).resolve().parents[2]
HISTORIES = ("post_introgression", "imported_old", "ancestral_retained")


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def canonical(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def write_json(path, value):
    path = Path(path)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    tmp.replace(path)


def seed(cfg, rep, stage):
    return int(canonical([cfg["seed_base"], rep, stage])[:16], 16) % (2**31 - 2) + 1


def contract(cfg, phase):
    # Commit/source hashes are provenance; explicit scientific schema controls semantics.
    keys = [k for k in cfg if k not in ("replicates", "workers", "memory_limit_gb")]
    return {"phase": phase, "config": {k: cfg[k] for k in keys}}


def cached(directory, name, expected):
    receipt = directory / (name + ".json")
    if not receipt.exists():
        return False
    record = json.loads(receipt.read_text())
    if record["contract"] != expected:
        raise ValueError(f"Incompatible cache: {receipt}; use a new output directory")
    for filename, checksum in record["outputs"].items():
        if (
            not (directory / filename).exists()
            or digest(directory / filename) != checksum
        ):
            raise ValueError(
                f"Corrupt or missing cached output: {directory / filename}"
            )
    return True


def complete(directory, name, expected, filenames, **details):
    write_json(
        directory / (name + ".json"),
        {
            "contract": expected,
            "outputs": {f: digest(directory / f) for f in filenames},
            "details": details,
        },
    )


def recombination_map(cfg):
    length = cfg["sequence_length"]
    boundaries = {0, length}
    for start in range(0, length, cfg["hotspot_spacing"]):
        boundaries.update((start, min(length, start + cfg["hotspot_width"])))
    positions = np.array(sorted(boundaries))
    multipliers = np.where(
        positions[:-1] % cfg["hotspot_spacing"] < cfg["hotspot_width"],
        cfg["hotspot_multiplier"],
        1.0,
    )
    rates = (
        multipliers
        * cfg["recombination_rate"]
        / np.average(multipliers, weights=np.diff(positions))
    )
    return msprime.RateMap(position=positions, rate=rates)


def introgressed_samples(tree, samples, migrations, position, pulse):
    nodes = [
        m.node
        for m in migrations
        if m.left <= position < m.right
        and m.source == 0
        and m.dest == 1
        and m.time == pulse
    ]
    return np.array(
        [any(tree.is_descendant(int(s), u) for u in nodes) for s in samples]
    )


def exact_bounds(tree, samples, carriers):
    derived, ancestral = np.asarray(samples)[carriers], np.asarray(samples)[~carriers]
    crown = reduce(tree.mrca, (int(x) for x in derived))
    lower = tree.time(crown)
    upper = min(tree.tmrca(int(derived[0]), int(x)) for x in ancestral)
    return lower, upper


def pair_classes(carriers, ancestry, cap, rng):
    c, n = np.flatnonzero(carriers), np.flatnonzero(~carriers)
    ci, ni = np.flatnonzero(carriers & ancestry), np.flatnonzero(~carriers & ancestry)
    candidates = {
        "CC": list(itertools.combinations(c, 2)),
        "CN": list(itertools.product(c, n)),
        "ICN": list(itertools.product(ci, ni)),
    }
    result = {}
    for name, pairs in candidates.items():
        pairs = sorted({tuple(sorted((int(a), int(b)))) for a, b in pairs})
        if len(pairs) > cap:
            pairs = [
                pairs[i] for i in sorted(rng.choice(len(pairs), cap, replace=False))
            ]
        result[name] = pairs
    return result


def simulate_one(cfg, output, rep):
    started = time.monotonic()
    directory = Path(output) / f"rep{rep:04d}"
    directory.mkdir(parents=True, exist_ok=True)
    spec = contract(cfg, "simulate/v1")
    if cached(directory, "simulation", spec):
        return {"replicate": rep, "phase": "simulate", "cached": True}
    gen = cfg["generation_time_years"]
    pulse, split = cfg["pulse_years"] / gen, cfg["archaic_split_years"] / gen
    n = 2 * cfg["sample_diploids"]
    ts = msprime.sim_ancestry(
        samples=[
            msprime.SampleSet(n, population="EAS", time=0, ploidy=1),
            msprime.SampleSet(
                2,
                population="Neanderthal",
                time=cfg["archaic_sample_years"] / gen,
                ploidy=1,
            ),
        ],
        demography=demography(cfg),
        sequence_length=cfg["sequence_length"],
        recombination_rate=recombination_map(cfg),
        record_migrations=True,
        model=[
            msprime.DiscreteTimeWrightFisher(duration=100),
            msprime.StandardCoalescent(),
        ],
        random_seed=seed(cfg, rep, "ancestry"),
    )
    modern = ts.samples()[:n]
    archaic = ts.samples()[n:]
    assert np.all(ts.tables.nodes.time[modern] == 0)
    migrations = list(ts.migrations())
    by_node = {}
    for m in sorted(migrations, key=lambda m: m.time):
        by_node.setdefault(m.node, []).append(m)
    # Continuous sites avoid recurrent-event conditioning in candidate generation.
    candidates = msprime.sim_mutations(
        ts,
        rate=cfg["candidate_mutation_rate"],
        random_seed=seed(cfg, rep, "candidate_events"),
        model=msprime.BinaryMutationModel(),
        discrete_genome=False,
    )
    rng = np.random.default_rng(seed(cfg, rep, "ascertainment"))
    sites = list(candidates.sites())
    rng.shuffle(sites)
    truth, observations, found = [], [], set()
    for site in sites:
        position = int(site.position)
        if not cfg["sequence_length"] * 0.2 <= position <= cfg["sequence_length"] * 0.8:
            continue
        mutation = site.mutations[0]
        tree = ts.at(position)
        carrier = np.array([tree.is_descendant(int(s), mutation.node) for s in modern])
        count = int(carrier.sum())
        if not cfg["minimum_carriers"] <= count <= n - cfg["minimum_carriers"]:
            continue
        age = mutation.time
        pop = mutation_population(ts, mutation, position, by_node)
        # Test population/time first to avoid unnecessarily mapping ancestry.
        possible = (
            (age < pulse and pop == 0)
            or (pulse <= age < split and pop == 1)
            or age >= split
        )
        if not possible:
            continue
        ancestry = introgressed_samples(tree, modern, migrations, position, pulse)
        archaic_carriers = sum(
            tree.is_descendant(int(s), mutation.node) for s in archaic
        )
        if age < pulse and np.all(ancestry[carrier]):
            history = "post_introgression"
        elif pulse <= age < split and pop == 1 and np.all(ancestry[carrier]):
            history = "imported_old"
        elif age >= split and archaic_carriers > 0 and np.sum(carrier & ~ancestry) >= 2:
            history = "ancestral_retained"
        else:
            continue
        if history in found or position in {o["position"] for o in observations}:
            continue
        found.add(history)
        low, high = exact_bounds(tree, modern, carrier)
        assert low <= age <= high
        label = f"r{rep:04d}_{len(observations)}"  # no origin encoded in observation ID
        classes = pair_classes(
            carrier,
            ancestry,
            cfg["pairs_per_class"],
            np.random.default_rng(seed(cfg, rep, f"pairs_{cfg['pairs_seed']}_{label}")),
        )
        observations.append(
            {
                "id": label,
                "replicate": rep,
                "position": position,
                "carriers": carrier.tolist(),
                "oracle_local_ancestry": ancestry.tolist(),
                "classes": classes,
            }
        )
        truth.append(
            {
                "id": label,
                "replicate": rep,
                "history": history,
                "age_years": age * gen,
                "lower_years": low * gen,
                "upper_years": high * gen,
                "carrier_count": count,
                "introgressed_count": int(ancestry.sum()),
                "introgressed_noncarriers": int(np.sum(ancestry & ~carrier)),
                "archaic_carriers": int(archaic_carriers),
                "candidate_position": site.position,
                "mutation_node": mutation.node,
                "mutation_population": pop,
            }
        )
        if len(found) == len(HISTORIES):
            break
    # Separate observation overlay; no candidate mutation is inserted in these SNP data.
    snps = msprime.sim_mutations(
        ts,
        rate=cfg["mutation_rate"],
        model=msprime.JC69(),
        random_seed=seed(cfg, rep, "observed_snps"),
    )
    ts.dump(directory / "ancestry.trees")
    snps.dump(directory / "snps.trees")
    # Pair modern haplotypes into diploids for the decoder, retaining exact node order.
    tables = snps.dump_tables()
    tables.migrations.clear()
    modern_ts = tables.tree_sequence().simplify(samples=modern)
    tables = modern_ts.dump_tables()
    tables.individuals.clear()
    individual = np.full(tables.nodes.num_rows, -1, dtype=np.int32)
    for i in range(n // 2):
        ind = tables.individuals.add_row()
        individual[2 * i : 2 * i + 2] = ind
    tables.nodes.individual = individual
    modern_ts = tables.tree_sequence()
    with (directory / "observed.vcf").open("w") as f:
        modern_ts.write_vcf(
            f, position_transform=lambda x: np.asarray(x, dtype=int) + 1
        )
    write_json(directory / "observations.json", observations)
    write_json(directory / "truth.json", truth)
    complete(
        directory,
        "simulation",
        spec,
        [
            "ancestry.trees",
            "snps.trees",
            "observed.vcf",
            "observations.json",
            "truth.json",
        ],
        replicate=rep,
        seeds={
            s: seed(cfg, rep, s)
            for s in ("ancestry", "candidate_events", "ascertainment", "observed_snps")
        },
        candidates=len(sites),
        found=sorted(found),
        missing=sorted(set(HISTORIES) - found),
        num_trees=ts.num_trees,
        snps=snps.num_sites,
        seconds=time.monotonic() - started,
    )
    return {
        "replicate": rep,
        "phase": "simulate",
        "found": sorted(found),
        "seconds": time.monotonic() - started,
    }


def read_raw(path):
    meta = json.loads(Path(str(path) + ".meta").read_text())
    n, p, chunk = meta["sequence_length"], meta["num_pairs"], meta["chunk_size"]
    arrays = [np.empty((n, p)) for _ in range(2)]
    block = n * chunk * 4
    with (
        Path(path).open("rb") as f,
        zstandard.ZstdDecompressor().stream_reader(f) as reader,
    ):
        for start in range(0, p, chunk):
            size = min(chunk, p - start)
            for array in arrays:
                content = reader.read(block)
                if len(content) != block:
                    raise ValueError("Truncated posterior stream")
                array[:, start : start + size] = np.frombuffer(
                    content, dtype=np.float32
                ).reshape(n, chunk)[:, :size]
        if reader.read(1):
            raise ValueError("Trailing posterior bytes")
    if not all(np.all(np.isfinite(a) & (a > 0)) for a in arrays):
        raise ValueError("Invalid posterior parameters")
    return meta, arrays[0], arrays[1]


def mask_intervals(length, positions, radius):
    intervals = sorted(
        (max(0, p - radius), min(length, p + radius + 1)) for p in positions
    )
    result, start = [], 0
    for left, right in intervals:
        if left > start:
            result.append((start, left))
        start = max(start, right)
    if start < length:
        result.append((start, length))
    return result


def decode_one(cfg, output, rep, decoder):
    directory = Path(output) / f"rep{rep:04d}"
    if not cached(directory, "simulation", contract(cfg, "simulate/v1")):
        raise ValueError(f"Simulation not complete: {directory}")
    observations = json.loads((directory / "observations.json").read_text())
    if not observations:
        return {"replicate": rep, "phase": "decode", "empty": True}
    pairs = sorted(
        {tuple(p) for o in observations for cls in o["classes"].values() for p in cls}
    )
    positions = sorted(o["position"] for o in observations)
    np.savetxt(directory / "pairs.txt", pairs, fmt="%d", delimiter="\t")
    np.savetxt(directory / "positions.txt", positions, fmt="%d")
    for radius in cfg["mask_radii"]:
        name = f"decode_mask{radius}"
        spec = {
            **contract(cfg, name + "/v1"),
            "vcf": digest(directory / "observed.vcf"),
            "observations": digest(directory / "observations.json"),
            "pairs": pairs,
            "positions": positions,
            "decoder": digest(decoder),
        }
        # JSON canonicalize tuples for exact cache comparisons.
        spec = json.loads(json.dumps(spec))
        if cached(directory, name, spec):
            continue
        mask = directory / f"mask{radius}.bed"
        mask.write_text(
            "".join(
                f"1\t{a}\t{b}\n"
                for a, b in mask_intervals(cfg["sequence_length"], positions, radius)
            )
        )
        raw = directory / f"posterior_mask{radius}.bin"
        receipt = run_within_decoder(
            decoder,
            directory / "observed.vcf",
            directory / f"summary_mask{radius}.tsv",
            input_format="vcf",
            only_within=False,
            pairs_file=directory / "pairs.txt",
            output_positions_file=directory / "positions.txt",
            output_at_stride=-1,
            raw_output=raw,
            mask=mask,
            threads=1,
            cache_size=cfg["decoder_cache_size"],
            scaled_mutation_rate=cfg["decoder_scaled_mutation_rate"],
            recombination_to_mutation_ratio=cfg[
                "decoder_recombination_to_mutation_ratio"
            ],
            mutation_rate=cfg["mutation_rate"],
            generation_time=cfg["generation_time_years"],
            threshold_years=cfg["pulse_years"],
        )
        write_json(directory / f"log_mask{radius}.json", receipt)
        meta, alpha, beta = read_raw(raw)
        assert meta["output_positions"] == positions
        assert [tuple(p) for p in meta["pairs"]] == pairs
        np.savez_compressed(
            directory / f"posterior_mask{radius}.npz",
            alpha=alpha,
            beta=beta,
            positions=positions,
            pairs=pairs,
        )
        complete(
            directory,
            name,
            spec,
            [
                raw.name,
                raw.name + ".meta",
                f"posterior_mask{radius}.npz",
                mask.name,
                f"summary_mask{radius}.tsv",
                f"log_mask{radius}.json",
            ],
        )
    return {
        "replicate": rep,
        "phase": "decode",
        "pairs": len(pairs),
        "positions": positions,
    }


def audit_one(cfg, output, rep):
    directory = Path(output) / f"rep{rep:04d}"
    if not cached(directory, "simulation", contract(cfg, "simulate/v1")):
        raise ValueError(f"Missing simulation {rep}")
    ts = tskit.load(directory / "ancestry.trees")
    snps = tskit.load(directory / "snps.trees")
    assert ts.num_samples == snps.num_samples == 2 * cfg["sample_diploids"] + 2
    assert (
        ts.num_trees == snps.num_trees and ts.sequence_length == cfg["sequence_length"]
    )
    observations = json.loads((directory / "observations.json").read_text())
    truth = {r["id"]: r for r in json.loads((directory / "truth.json").read_text())}
    modern = ts.samples()[: 2 * cfg["sample_diploids"]]
    migrations = list(ts.migrations())
    by_node = {}
    for m in sorted(migrations, key=lambda m: m.time):
        by_node.setdefault(m.node, []).append(m)
    for o in observations:
        t = truth[o["id"]]
        tree = ts.at(o["position"])
        carriers = np.array(
            [tree.is_descendant(int(s), t["mutation_node"]) for s in modern]
        )
        assert carriers.tolist() == o["carriers"]
        low, high = exact_bounds(tree, modern, carriers)
        age = t["age_years"] / cfg["generation_time_years"]
        assert low <= age <= high
        assert (
            tree.time(t["mutation_node"])
            <= age
            <= tree.time(tree.parent(t["mutation_node"]))
        )
        ancestry = introgressed_samples(
            tree,
            modern,
            migrations,
            o["position"],
            cfg["pulse_years"] / cfg["generation_time_years"],
        )
        assert ancestry.tolist() == o["oracle_local_ancestry"]
        for kind, pairs in o["classes"].items():
            for a, b in pairs:
                assert a < b < len(modern)
                assert (
                    carriers[a] and carriers[b]
                    if kind == "CC"
                    else carriers[a] != carriers[b]
                )
                if kind == "ICN":
                    assert ancestry[a] and ancestry[b]
    result = {
        "replicate": rep,
        "phase": "audit",
        "events": len(observations),
        "verified": True,
    }
    write_json(directory / "audit.json", result)
    return result


def age_grid(cfg):
    edges = np.geomspace(
        cfg["age_grid_min_years"], cfg["age_grid_max_years"], cfg["age_grid_size"] + 1
    )
    return np.sqrt(edges[:-1] * edges[1:]), np.diff(edges)


def distribution_summary(grid, log_mass, pulse):
    weights = np.exp(log_mass - logsumexp(log_mass))
    cumulative = np.cumsum(weights)
    low, median, high = np.interp([0.05, 0.5, 0.95], cumulative, grid)
    return dict(
        estimate=median, low=low, high=high, p_old=float(weights[grid > pulse].sum())
    )


def composite_curves(alpha, beta, classes, grid, scale):
    terms = {}
    for name in ("CC", "CN", "ICN"):
        idx = classes[name]
        if name == "ICN":
            idx = [i for i in idx if i not in set(classes["CN"])]
        if not idx:
            terms[name] = np.zeros((0, len(grid)))
            continue
        function = gammainc if name == "CC" else gammaincc
        probabilities = function(
            alpha[idx, None], beta[idx, None] * grid[None, :] / scale
        )
        terms[name] = np.log(np.clip(probabilities, 1e-250, 1))
    return terms


def observation_features(alpha, beta, cc, cn, count, n, scale, pulse):
    # No truth ages, migration labels, history categories or coordinates enter kNN.
    features = [math.log(count), math.log(n - count)]
    for indices in (cc, cn):
        means = alpha[indices] / beta[indices] * scale
        features.extend(np.log(np.quantile(means, [0, 0.1, 0.25, 0.5, 0.75, 0.9, 1])))
        features.append(
            float(np.mean(gammainc(alpha[indices], beta[indices] * pulse / scale)))
        )
    return features


def split_name(rep):
    return "train" if rep % 4 < 2 else ("calibration" if rep % 4 == 2 else "test")


def conformal_padding(y, low, high, alpha):
    scores = np.maximum.reduce((low - y, y - high, np.zeros_like(y)))
    rank = math.ceil((len(scores) + 1) * (1 - alpha))
    if rank > len(scores):
        raise ValueError("Too few calibration regions for requested coverage")
    return float(np.sort(scores)[rank - 1])


def analyze(cfg, output, replicates):
    output = Path(output)
    records, features, truths = [], [], []
    grid, widths = age_grid(cfg)
    scale = (
        cfg["decoder_scaled_mutation_rate"]
        / (2 * cfg["mutation_rate"])
        * cfg["generation_time_years"]
    )
    for rep in range(replicates):
        directory = output / f"rep{rep:04d}"
        if not cached(directory, "simulation", contract(cfg, "simulate/v1")):
            raise ValueError(f"Missing simulation {rep}")
        observations = json.loads((directory / "observations.json").read_text())
        truths.extend(json.loads((directory / "truth.json").read_text()))
        for radius in cfg["mask_radii"]:
            if not observations:
                continue
            receipt = json.loads((directory / f"decode_mask{radius}.json").read_text())
            cached(directory, f"decode_mask{radius}", receipt["contract"])
            with np.load(
                directory / f"posterior_mask{radius}.npz", allow_pickle=False
            ) as data:
                pairs = {tuple(p): i for i, p in enumerate(data["pairs"])}
                positions = list(data["positions"])
                for o in observations:
                    row = positions.index(o["position"])
                    a, b = data["alpha"][row], data["beta"][row]
                    classes = {
                        k: [pairs[tuple(p)] for p in v] for k, v in o["classes"].items()
                    }
                    terms = composite_curves(a, b, classes, grid, scale)
                    base = {
                        "id": o["id"],
                        "replicate": rep,
                        "mask_radius": radius,
                        "split": split_name(rep),
                    }
                    means = a[classes["CC"]] / b[classes["CC"]] * scale
                    records.append(
                        {
                            **base,
                            "method": "carrier_mean",
                            "estimate": float(np.mean(means)),
                            "low": np.nan,
                            "high": np.nan,
                            "p_old": np.nan,
                        }
                    )
                    for method, curve in (
                        ("composite_raw", terms["CC"].sum(0) + terms["CN"].sum(0)),
                        (
                            "composite_balanced",
                            terms["CC"].mean(0) + terms["CN"].mean(0),
                        ),
                        (
                            "composite_with_ancestry",
                            terms["CC"].sum(0)
                            + terms["CN"].sum(0)
                            + terms["ICN"].sum(0),
                        ),
                    ):
                        records.append(
                            {
                                **base,
                                "method": method,
                                **distribution_summary(
                                    grid, curve + np.log(widths), cfg["pulse_years"]
                                ),
                            }
                        )
                    features.append(
                        {
                            **base,
                            "log_compatibility": terms["CC"].sum(0)
                            + terms["CN"].sum(0),
                            "x": observation_features(
                                a,
                                b,
                                classes["CC"],
                                classes["CN"],
                                sum(o["carriers"]),
                                len(o["carriers"]),
                                scale,
                                cfg["pulse_years"],
                            ),
                        }
                    )
    truth = pd.DataFrame(truths)
    write_json(
        output / "analysis_config.json",
        {
            "config": cfg,
            "replicates": replicates,
            "split": "region modulo 4: 0/1 train, 2 calibration, 3 test",
            "uniform_age_prior_years": [
                cfg["age_grid_min_years"],
                cfg["age_grid_max_years"],
            ],
        },
    )
    feature_frame = pd.DataFrame(features).merge(
        truth[["id", "age_years"]], on="id", validate="many_to_one"
    )
    calibration_details = []
    for radius, frame in feature_frame.groupby("mask_radius"):
        x_full = np.stack(frame.x)
        y = np.log(frame.age_years.to_numpy())
        train, cal, test = [
            np.flatnonzero(frame.split.to_numpy() == s)
            for s in ("train", "calibration", "test")
        ]
        if len(cal) < 10 or len(train) < cfg["knn_neighbors"]:
            continue
        tempered, tempered_details = fit_tempered_composite(
            frame, grid, widths, cfg, train, cal, test
        )
        records.extend(tempered)
        calibration_details.append(tempered_details)
        for feature_mode in ("full", "frequency_only"):
            x = x_full if feature_mode == "full" else x_full[:, :2]
            location, spread = np.median(x[train], axis=0), np.std(x[train], axis=0)
            spread = np.maximum(spread, 0.1)
            distance = cdist((x - location) / spread, (x[train] - location) / spread)
            nearest = np.argsort(distance, axis=1)[
                :, : min(cfg["knn_neighbors"], len(train))
            ]
            ages = y[train][nearest]
            low, median, high = np.quantile(ages, [0.05, 0.5, 0.95], axis=1)
            # Region-cluster calibration: maximum error across that region's selected SVs.
            cal_scores = (
                pd.DataFrame(
                    {
                        "rep": frame.iloc[cal].replicate.to_numpy(),
                        "score": np.maximum.reduce(
                            (low[cal] - y[cal], y[cal] - high[cal], np.zeros(len(cal)))
                        ),
                    }
                )
                .groupby("rep")
                .score.max()
                .to_numpy()
            )
            rank = math.ceil((len(cal_scores) + 1) * (1 - cfg["interval_alpha"]))
            if rank > len(cal_scores):
                continue
            padding = float(np.sort(cal_scores)[rank - 1])
            calibration_details.append(
                {
                    "mask_radius": int(radius),
                    "training_events": len(train),
                    "feature_mode": feature_mode,
                    "calibration_regions": len(cal_scores),
                    "log_padding": padding,
                    "neighbors": cfg["knn_neighbors"],
                }
            )
            suffix = "" if feature_mode == "full" else "_frequency"
            for j in np.r_[cal, test]:
                row = frame.iloc[j]
                for method, pad in (
                    ("knn" + suffix, 0.0),
                    ("knn" + suffix + "_conformal", padding),
                ):
                    records.append(
                        {
                            "id": row.id,
                            "replicate": int(row.replicate),
                            "mask_radius": int(radius),
                            "split": row.split,
                            "method": method,
                            "estimate": float(np.exp(median[j])),
                            "low": float(np.exp(low[j] - pad)),
                            "high": float(np.exp(high[j] + pad)),
                            "p_old": float(
                                np.mean(np.exp(ages[j]) > cfg["pulse_years"])
                            ),
                        }
                    )
    predictions = pd.DataFrame(records).merge(
        truth, on=["id", "replicate"], validate="many_to_one"
    )
    predictions["abs_log2_error"] = abs(
        np.log2(predictions.estimate / predictions.age_years)
    )
    predictions["covered"] = (predictions.low <= predictions.age_years) & (
        predictions.age_years <= predictions.high
    )
    predictions["correct_side"] = (predictions.estimate > cfg["pulse_years"]) == (
        predictions.age_years > cfg["pulse_years"]
    )
    predictions["confident_call"] = (predictions.low > cfg["pulse_years"]) | (
        predictions.high < cfg["pulse_years"]
    )
    predictions["wrong_confident"] = predictions.confident_call & ~np.where(
        predictions.low > cfg["pulse_years"],
        predictions.age_years > cfg["pulse_years"],
        predictions.age_years < cfg["pulse_years"],
    )
    truth.to_csv(output / "truth_summary.tsv", sep="\t", index=False)
    predictions.to_csv(output / "predictions.tsv", sep="\t", index=False)
    feature_frame.drop(columns=["x", "log_compatibility"]).join(
        pd.DataFrame(feature_frame.x.tolist()).add_prefix("feature_")
    ).to_csv(output / "features.tsv", sep="\t", index=False)
    write_json(output / "calibration.json", calibration_details)
    metrics = []
    test = predictions[predictions.split == "test"]
    for (radius, method), group in test.groupby(["mask_radius", "method"]):
        for history in ("ALL", "POST_OR_IMPORTED", *HISTORIES):
            g = (
                group
                if history == "ALL"
                else group[group.history != "ancestral_retained"]
                if history == "POST_OR_IMPORTED"
                else group[group.history == history]
            )
            if not len(g):
                continue
            metrics.append(
                {
                    "mask_radius": radius,
                    "method": method,
                    "history": history,
                    "n": len(g),
                    "median_fold_error": 2 ** g.abs_log2_error.median(),
                    "within_2fold": float((g.abs_log2_error <= 1).mean()),
                    "median_bias_ratio": float((g.estimate / g.age_years).median()),
                    "side_accuracy": g.correct_side.mean(),
                    "interval_coverage": g.covered.mean()
                    if g.low.notna().all()
                    else np.nan,
                    "median_interval_fold_width": (g.high / g.low).median(),
                    "call_fraction": g.confident_call.mean(),
                    "wrong_confident_fraction": g.wrong_confident.mean(),
                    "error_among_calls": g.loc[
                        g.confident_call, "wrong_confident"
                    ].mean(),
                }
            )
    pd.DataFrame(metrics).to_csv(output / "metrics.tsv", sep="\t", index=False)
    oracle = truth.copy()
    oracle["split"] = oracle.replicate.map(split_name)
    oracle["spans_pulse"] = (oracle.lower_years < cfg["pulse_years"]) & (
        oracle.upper_years > cfg["pulse_years"]
    )
    oracle["crown_younger_than_pulse"] = oracle.lower_years < cfg["pulse_years"]
    oracle["bound_fold_width"] = oracle.upper_years / oracle.lower_years
    oracle["bound_contains_age"] = (oracle.lower_years <= oracle.age_years) & (
        oracle.age_years <= oracle.upper_years
    )
    assert oracle.bound_contains_age.all()
    oracle.to_csv(output / "oracle_bounds.tsv", sep="\t", index=False)
    oracle[oracle.split == "test"].groupby("history").agg(
        n=("id", "count"),
        ambiguous_fraction=("spans_pulse", "mean"),
        young_crown_fraction=("crown_younger_than_pulse", "mean"),
        median_bound_fold_width=("bound_fold_width", "median"),
        median_age_years=("age_years", "median"),
        median_carrier_count=("carrier_count", "median"),
        fraction_with_introgressed_noncarriers=(
            "introgressed_noncarriers",
            lambda s: (s > 0).mean(),
        ),
    ).to_csv(output / "oracle_summary.tsv", sep="\t")
    bootstrap_metrics(test, output, cfg)
    plot_results(output, predictions)
    return {
        "phase": "analyze",
        "events": len(truth),
        "held_out_events": test.id.nunique(),
        "calibration": calibration_details,
    }


def fit_tempered_composite(frame, grid, widths, cfg, train, cal, test):
    """Choose dependence attenuation and age prior using training loci only."""
    curves = np.stack(frame.log_compatibility)
    true_ages = frame.age_years.to_numpy()
    candidate_fits = []
    for temperature in (1.0, 3.0, 10.0, 30.0, 100.0):
        for prior in ("uniform", "log_uniform"):
            log_density = curves / temperature - (
                np.log(grid) if prior == "log_uniform" else 0
            )
            log_density -= logsumexp(log_density + np.log(widths), axis=1)[:, None]
            loss = -np.mean(
                [np.interp(true_ages[i], grid, log_density[i]) for i in train]
            )
            candidate_fits.append((loss, temperature, prior))
    _, temperature, prior = min(candidate_fits)
    log_mass = (
        curves / temperature
        + np.log(widths)
        - (np.log(grid) if prior == "log_uniform" else 0)
    )
    summaries = [distribution_summary(grid, m, cfg["pulse_years"]) for m in log_mass]
    low = np.log([s["low"] for s in summaries])
    high = np.log([s["high"] for s in summaries])
    y = np.log(true_ages)
    scores = (
        pd.DataFrame(
            {
                "rep": frame.iloc[cal].replicate.to_numpy(),
                "score": np.maximum.reduce(
                    (low[cal] - y[cal], y[cal] - high[cal], np.zeros(len(cal)))
                ),
            }
        )
        .groupby("rep")
        .score.max()
        .to_numpy()
    )
    rank = math.ceil((len(scores) + 1) * (1 - cfg["interval_alpha"]))
    if rank > len(scores):
        return [], {"method": "composite_tempered", "insufficient_calibration": True}
    padding = float(np.sort(scores)[rank - 1])
    records = []
    for i in np.r_[cal, test]:
        row = frame.iloc[i]
        for name, pad in (
            ("composite_tempered", 0.0),
            ("composite_tempered_conformal", padding),
        ):
            records.append(
                {
                    "id": row.id,
                    "replicate": int(row.replicate),
                    "mask_radius": int(row.mask_radius),
                    "split": row.split,
                    "method": name,
                    **summaries[i],
                    "low": float(np.exp(low[i] - pad)),
                    "high": float(np.exp(high[i] + pad)),
                }
            )
    return records, {
        "method": "composite_tempered",
        "mask_radius": int(frame.iloc[0].mask_radius),
        "temperature": temperature,
        "prior": prior,
        "training_events": len(train),
        "calibration_regions": len(scores),
        "log_padding": padding,
        "training_grid": [
            {"loss": float(v), "temperature": t, "prior": p}
            for v, t, p in candidate_fits
        ],
    }


def bootstrap_metrics(test, output, cfg):
    """Region bootstrap, preserving dependence among events within a region."""
    rows = []
    for (radius, method), group in test.groupby(["mask_radius", "method"]):
        for scope in ("ALL", "POST_OR_IMPORTED"):
            g = (
                group
                if scope == "ALL"
                else group[group.history != "ancestral_retained"]
            )
            reps = np.unique(g.replicate)
            rng = np.random.default_rng(cfg["pairs_seed"])
            if len(reps) < 2:
                continue
            by_rep = {r: g[g.replicate == r] for r in reps}
            samples = []
            for _ in range(500):
                draw = pd.concat(
                    [by_rep[r] for r in rng.choice(reps, len(reps), replace=True)]
                )
                samples.append(
                    [
                        2 ** draw.abs_log2_error.median(),
                        (draw.abs_log2_error <= 1).mean(),
                        draw.covered.mean() if draw.low.notna().all() else np.nan,
                        draw.correct_side.mean(),
                    ]
                )
            samples = np.array(samples)
            for j, name in enumerate(
                (
                    "median_fold_error",
                    "within_2fold",
                    "interval_coverage",
                    "side_accuracy",
                )
            ):
                if np.isnan(samples[:, j]).all():
                    continue
                low, high = np.quantile(samples[:, j], [0.025, 0.975])
                rows.append(
                    {
                        "mask_radius": radius,
                        "method": method,
                        "scope": scope,
                        "metric": name,
                        "bootstrap_low": low,
                        "bootstrap_high": high,
                        "independent_regions": len(reps),
                    }
                )
    pd.DataFrame(rows).to_csv(output / "metric_bootstrap.tsv", sep="\t", index=False)


def plot_results(output, predictions):
    import matplotlib
    import matplotlib.pyplot as plt

    matplotlib.rcParams.update(
        {
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.size": 17,
            "axes.titlesize": 19,
            "axes.labelsize": 18,
            "xtick.labelsize": 14,
            "ytick.labelsize": 14,
        }
    )
    colors = dict(zip(HISTORIES, ["#2780a5", "#d56a26", "#8064a2"]))
    test = predictions[(predictions.split == "test") & (predictions.mask_radius == 0)]
    methods = [
        m
        for m in (
            "carrier_mean",
            "composite_raw",
            "composite_balanced",
            "knn_conformal",
        )
        if m in set(test.method)
    ]
    fig, axes = plt.subplots(2, 2, figsize=(11, 10), layout="constrained")
    for ax, method in zip(axes.flat, methods):
        for history, color in colors.items():
            g = test[(test.method == method) & (test.history == history)]
            ax.scatter(
                g.age_years / 1000,
                g.estimate / 1000,
                s=28,
                alpha=0.75,
                color=color,
                label=history.replace("_", " "),
            )
        ax.plot([0.1, 5000], [0.1, 5000], color=".5", linestyle="--", linewidth=1)
        ax.axvline(50, color=".6", linestyle=":")
        ax.axhline(50, color=".6", linestyle=":")
        ax.set(
            xscale="log",
            yscale="log",
            xlim=(0.1, 5000),
            ylim=(0.1, 5000),
            title=method.replace("_", " "),
            xlabel="True mutation age (kya)",
            ylabel="Estimated age (kya)",
        )
    handles, _ = axes.flat[0].get_legend_handles_labels()
    fig.legend(
        handles,
        ["New on introgressed", "Old imported", "Ancestral retained"],
        loc="outside upper center",
        ncol=3,
        fontsize=14,
        frameon=False,
    )
    fig.savefig(output / "age_estimates.png", dpi=170)
    fig.savefig(output / "age_estimates.pdf")
    plt.close(fig)
    metrics = pd.read_csv(output / "metrics.tsv", sep="\t")
    g = metrics[(metrics.history == "ALL") & (metrics.mask_radius == 0)].dropna(
        subset=["interval_coverage"]
    )
    fig, ax = plt.subplots(figsize=(11, 7.5), layout="constrained")
    x = np.arange(len(g))
    ax.barh(
        x - 0.18,
        g.interval_coverage,
        height=0.36,
        label="90% interval coverage",
        color="#2780a5",
    )
    ax.barh(
        x + 0.18,
        g.call_fraction,
        height=0.36,
        label="Pre/post call fraction",
        color="#d56a26",
    )
    ax.axvline(0.9, color=".4", linestyle="--")
    method_labels = {
        "composite_raw": "CDF: all pairs",
        "composite_balanced": "CDF: balanced",
        "composite_with_ancestry": "CDF: + oracle ancestry",
        "composite_tempered": "CDF: trained attenuation",
        "composite_tempered_conformal": "CDF: trained + calibration",
        "knn": "kNN",
        "knn_conformal": "kNN + calibration",
        "knn_frequency": "Frequency only",
        "knn_frequency_conformal": "Frequency + calibration",
    }
    ax.set(
        xlim=(0, 1.03),
        xlabel="Fraction of held-out SVs",
        yticks=x,
        yticklabels=[method_labels.get(s, s) for s in g.method],
    )
    ax.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 1.01),
        ncol=2,
        frameon=False,
        fontsize=15,
    )
    fig.savefig(output / "interval_calibration.png", dpi=170)
    fig.savefig(output / "interval_calibration.pdf")
    plt.close(fig)
    # Prespecified first twelve held-out events per focal history; no selection by error.
    fig, axes = plt.subplots(1, 2, figsize=(11, 8.5), layout="constrained", sharex=True)
    examples = test.drop_duplicates("id").sort_values("replicate")
    for ax, history, title in zip(
        axes,
        ["imported_old", "post_introgression"],
        ["Old imported SVs", "New mutations"],
    ):
        g = examples[examples.history == history].head(12)
        y = np.arange(len(g))
        ax.hlines(
            y,
            g.lower_years / 1000,
            g.upper_years / 1000,
            color=colors[history],
            linewidth=3,
            label="Exact tree bounds",
        )
        ax.scatter(
            g.age_years / 1000,
            y,
            color="black",
            s=45,
            zorder=3,
            label="True mutation age",
        )
        ax.axvline(
            50, color=".4", linestyle="--", linewidth=2, label="Introgression (50 kya)"
        )
        ax.set(
            xscale="log",
            xlabel="Age (kya)",
            title=title,
            yticks=y,
            yticklabels=[f"Region {int(r)}" for r in g.replicate],
        )
        ax.invert_yaxis()
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles, labels, loc="outside upper center", ncol=3, frameon=False, fontsize=13
    )
    fig.savefig(output / "oracle_age_bounds.png", dpi=170)
    fig.savefig(output / "oracle_age_bounds.pdf")
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path, default=REPO / "configs/sv_mutation_age.json"
    )
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--decoder", type=Path, default=REPO / "bin/gamma_smc")
    parser.add_argument(
        "--phase",
        choices=["simulate", "audit", "decode", "analyze", "run"],
        default="run",
    )
    parser.add_argument("--replicates", type=int)
    args = parser.parse_args(argv)
    cfg = json.loads(args.config.read_text())
    reps = args.replicates or cfg["replicates"]
    args.out.mkdir(parents=True, exist_ok=True)
    if cfg["workers"] > 24 or cfg["workers"] < 1:
        raise ValueError("workers must lie in 1..24")
    # Linux pilot workers inherit a total address-space budget divided among them.
    if platform.system() == "Linux":
        import resource

        limit = int(cfg["memory_limit_gb"] * 1024**3 / cfg["workers"])
        resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
    provenance = {
        "config": cfg,
        "replicates": reps,
        "output": str(args.out.resolve()),
        "phase": args.phase,
        "versions": {
            "python": platform.python_version(),
            "msprime": msprime.__version__,
            "tskit": tskit.__version__,
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "source_sha256": digest(__file__),
        "git_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=REPO, text=True
        ).strip(),
        "decoder_path": str(args.decoder.resolve()),
        "decoder_sha256": digest(args.decoder),
        "time_unix": time.time(),
    }
    # Prevent concurrent writers without leaving a stale lock after failure.
    import fcntl

    with (args.out / "study.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        write_json(args.out / "provenance.json", provenance)
        write_json(args.out / f"provenance_{time.time_ns()}.json", provenance)
        for phase in (
            ["simulate", "audit", "decode", "analyze"]
            if args.phase == "run"
            else [args.phase]
        ):
            if phase == "analyze":
                result = analyze(cfg, args.out, reps)
                print(json.dumps(result), flush=True)
                write_json(args.out / "status.json", result)
                continue
            function = {
                "simulate": simulate_one,
                "audit": audit_one,
                "decode": decode_one,
            }[phase]
            with ProcessPoolExecutor(max_workers=cfg["workers"]) as pool:
                futures = [
                    pool.submit(
                        function,
                        cfg,
                        str(args.out),
                        rep,
                        *([str(args.decoder)] if phase == "decode" else []),
                    )
                    for rep in range(reps)
                ]
                for done, future in enumerate(as_completed(futures), 1):
                    result = future.result()
                    print(json.dumps(result), flush=True)
                    write_json(
                        args.out / "status.json",
                        {
                            "phase": phase,
                            "completed": done,
                            "total": reps,
                            "latest": result,
                        },
                    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
