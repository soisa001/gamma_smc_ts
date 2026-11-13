"""Assembly-to-offspring probability benchmark; evaluator truth never enters fitting."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
import hashlib
import importlib.metadata
import json
import logging
import os
from pathlib import Path
import platform
import re
import subprocess
import sys
import time

# Set BEFORE importing NumPy/Matplotlib, including on spawned Windows workers.
for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_key] = "1"
os.environ["MPLBACKEND"] = "Agg"
import edlib
import numpy as np
from scipy.special import xlogy

BASES = np.array(list("ACGT"))
COMPLEMENT = np.array([3, 2, 1, 0])
METHODS = ("1mer", "3mer", "5mer", "7mer", "7mer_flat")
SCENARIOS = ("constant", "three", "five_iid", "five_shift", "seven_shift", "empirical_five")
SCHEMA = 1
FIT_VERSION = 2
EVALUATION_VERSION = 2
SHRINKAGE = (10.0, 100.0, 1000.0, 10000.0, 100000.0, 1e12)
LOG = logging.getLogger("benchmark")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(obj, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    temp.replace(path)


def digits(k):
    return ((np.arange(4**k)[:, None] // 4 ** np.arange(k - 1, -1, -1)) % 4).astype(np.int8)


def encode(array):
    return (array * 4 ** np.arange(array.shape[-1] - 1, -1, -1)).sum(axis=-1)


def rc_indices(k):
    return encode(3 - digits(k)[:, ::-1])


def project(indices, k):
    return (indices // 4 ** ((7 - k) // 2)) % (4**k)


def sequence(rng, length, outer_match):
    """Four independent stationary chains: uniform 4mers, P(X_i=X_i-4)=q."""
    steps = rng.integers(1, 4, size=((length + 3) // 4, 4))
    steps[rng.random(steps.shape) < outer_match] = 0
    steps[0] = rng.integers(0, 4, size=4)
    return (np.cumsum(steps, axis=0) % 4).astype(np.int8).reshape(-1)[:length]


def contexts(seq):
    return encode(np.lib.stride_tricks.sliding_window_view(seq, 7))


def make_kernel(scenario, mean_rate):
    d = digits(7)
    ancestral = d[:, 3]
    rates = np.ones((4**7, 4), dtype=float)
    rates[np.arange(4**7), ancestral] = 0
    rates[np.arange(4**7), ancestral ^ 2] *= 2  # A<->G and C<->T
    if scenario != "constant":
        cpg_c = (ancestral == 1) & (d[:, 4] == 2)
        cpg_g = (ancestral == 2) & (d[:, 2] == 1)
        rates[cpg_c, 3] *= 6
        rates[cpg_g, 0] *= 6
    if scenario in ("five_iid", "five_shift", "seven_shift"):
        outer = d[:, 1] == d[:, 5]
        weak = (ancestral == 0) | (ancestral == 3)
        rates[outer & weak, 1] *= 4
        rates[outer & weak, 2] *= 4
    if scenario == "seven_shift":
        outer7 = d[:, 0] == d[:, 6]
        rates[outer7 & (ancestral == 1), 0] *= 12
        rates[outer7 & (ancestral == 2), 3] *= 12
    if scenario == "empirical_five":
        rates[:] = 0
        path = Path(__file__).parent / "resources" / "ERV_7bp_rates.txt"
        provenance = json.loads(path.with_suffix(".provenance.json").read_text())
        if sha256(path) != provenance["sha256"]:
            raise ValueError("Empirical source table failed its recorded checksum")
        with path.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                motif = row["Sequence"].split("(")[0]
                b = np.array(["ACGT".index(x) for x in motif])
                i = int(encode(b))
                j = int(encode(3 - b[::-1]))
                mapping = {0: {1: "AT_CG", 2: "AT_GC", 3: "AT_TA"},
                           1: {0: "GC_TA", 2: "GC_CG", 3: "GC_AT"}}[int(b[3])]
                for alt, column in mapping.items():
                    rates[i, alt] = float(row[column])
                    rates[j, 3 - alt] = float(row[column])
        reduced = np.zeros((4**5, 4))
        np.add.at(reduced, project(np.arange(4**7), 5), rates / 16)
        rates = reduced[project(np.arange(4**7), 5)]
    rates *= mean_rate / rates.sum(axis=1).mean()
    if rates.sum(axis=1).max() >= 0.5:
        raise ValueError("Choose a lower mean-rate; per-site mutation probability must be <0.5")
    rates[np.arange(4**7), ancestral] = 1 - rates.sum(axis=1)
    if not np.allclose(rates, rates[rc_indices(7)][:, COMPLEMENT]):
        raise AssertionError("Kernel is not reverse-complement symmetric")
    return rates


def mutate(rng, parent, kernel):
    probs = kernel[contexts(parent)]
    child = parent.copy()
    child[3:-3] = (rng.random(len(probs))[:, None] > np.cumsum(probs, axis=1)).sum(axis=1)
    return child


def write_fasta(path, records):
    with Path(path).open("w", encoding="ascii", newline="\n") as handle:
        for name, seq in records:
            handle.write(f">{name}\n")
            letters = "".join(BASES[seq])
            for start in range(0, len(letters), 80):
                handle.write(letters[start:start + 80] + "\n")


def read_fasta(path):
    records, name, parts = {}, None, []
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    records[name] = np.array(["ACGT".index(x) for x in "".join(parts)], dtype=np.int8)
                name, parts = line[1:].strip(), []
            else:
                parts.append(line.strip())
    if name is not None:
        records[name] = np.array(["ACGT".index(x) for x in "".join(parts)], dtype=np.int8)
    if not records:
        raise ValueError(f"Empty FASTA: {path}")
    return records


def aligned_calls(parent, child):
    """Global edit-distance alignment, then parent-coordinate calls and gap masking."""
    ref, query = "".join(BASES[parent]), "".join(BASES[child])
    result = edlib.align(query, ref, mode="NW", task="path")
    called = np.full(len(parent), -1, dtype=np.int8)
    bad = np.zeros(len(parent), dtype=bool)
    r = q = 0
    for number, op in re.findall(r"(\d+)([=XIDM])", result["cigar"]):
        n = int(number)
        if op in "=XM":
            called[r:r+n] = child[q:q+n]
            r, q = r+n, q+n
        elif op == "I":
            bad[max(0, r-3):min(len(parent), r+4)] = True
            q += n
        else:
            bad[max(0, r-3):min(len(parent), r+n+3)] = True
            r += n
    assert r == len(parent) and q == len(child)
    valid = (called >= 0) & ~bad
    valid[:3] = valid[-3:] = False
    return called, valid, result["editDistance"]


TASK_PROMPT = """# Predict offspring alleles

You are given independently assembled, haploid parent and offspring DNA sequences
from a synthetic inheritance process. Each calibration pair is separated by one
transmission. Assemblies are error-free, and matching contig identifiers identify
homologous sequence. The same process applies to the new parent assembly. Only
single-nucleotide changes are generated. Rates in this synthetic system are higher
than typical human germline rates. Predict the offspring's allele on the transmitted
haplotype; diploid segregation is outside this task.

Use calibration_parent.fasta and calibration_offspring.fasta to learn the process.
For each coordinate in queries.tsv, predict the allele in one offspring of
new_parent.fasta. Coordinates are one-based parent-assembly coordinates. Return
predictions.tsv with columns contig, position, p_A, p_C, p_G, p_T, in query order.
Each probability must be finite and nonnegative, and each row must sum to one.
Include the possibility that the offspring allele equals the parent's allele.
Also return the expected total number of changed alleles across the queried sites.
Explain your estimation and validation approach. The offspring outcomes for the
new parent are withheld.
"""


def manifest_valid(directory, filename, expected=None):
    path = directory / filename
    if not path.exists():
        return False
    manifest = json.loads(path.read_text())
    if expected is not None and manifest["config"] != expected:
        return False
    return all((directory / name).is_file() and sha256(directory / name) == value
               for name, value in manifest["files"].items())


def seal(directory, filename, config, files, extra=None):
    result = {"config": config, "files": {str(p.relative_to(directory)): sha256(p) for p in files}}
    if extra:
        result.update(extra)
    write_json(directory / filename, result)


def prepare_one(job):
    root, scenario, seed, config = job
    directory = Path(root) / scenario / f"seed_{seed:03d}"
    cfg = {**config, "scenario": scenario, "seed": seed, "schema": SCHEMA}
    if scenario == "empirical_five":
        cfg["empirical_table_sha256"] = sha256(Path(__file__).parent / "resources" / "ERV_7bp_rates.txt")
    if manifest_valid(directory, "samples.manifest.json", cfg):
        return f"cached samples {scenario} seed={seed}"
    public, private = directory / "public", directory / "private"
    public.mkdir(parents=True, exist_ok=True)
    private.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    kernel = make_kernel(scenario, config["mean_rate"])
    target_q = 0.25 if scenario in ("constant", "three", "five_iid") else config["target_outer_match"]
    records = {}
    for group, total, q in (("calibration", config["train_bases"], 0.25),
                            ("new", config["test_bases"], target_q)):
        parents, children = [], []
        for index, start in enumerate(range(0, total, config["contig_bases"])):
            parent = sequence(rng, min(config["contig_bases"], total-start), q)
            name = f"scaffold_{index:04d}"
            parents.append((name, parent))
            children.append((name, mutate(rng, parent, kernel)))
        write_fasta(public / f"{group}_parent.fasta", parents)
        write_fasta((public if group == "calibration" else private) / f"{group}_offspring.fasta", children)
        records[group] = parents
    with (public / "queries.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["contig", "position"])
        for name, seq in records["new"]:
            for pos in range(3, len(seq)-3):
                writer.writerow([name, pos+1])
    (public / "TASK.md").write_text(TASK_PROMPT, encoding="utf-8")
    write_json(public / "metadata.json", {"transmissions_per_pair": 1, "coordinates": "1-based",
                                          "ploidy": "haploid transmitted allele", "assembly_error_rate": 0})
    np.save(private / "kernel.npy", kernel)
    seal(directory, "samples.manifest.json", cfg, sorted(public.glob("*")) + sorted(private.glob("*")))
    return f"prepared {scenario} seed={seed}"


def collect(public, return_blocks=False):
    parents = read_fasta(public / "calibration_parent.fasta")
    children = read_fasta(public / "calibration_offspring.fasta")
    train, validation = np.zeros((4**7, 4)), np.zeros((4**7, 4))
    validation_blocks = []
    audit = {"true_mutations": 0, "called_mutations": 0, "correct_mutations": 0,
             "callable_sites": 0, "masked_sites": 0, "alignment_edit_distance": 0}
    for index, (name, parent) in enumerate(parents.items()):
        child = children[name]
        called, valid, distance = aligned_calls(parent, child)
        ids = contexts(parent)[valid[3:-3]]
        outcome = called[valid]
        matrix = validation if index % 5 == 0 else train
        chunk = np.bincount(ids*4 + outcome, minlength=4**7*4).reshape(-1, 4)
        matrix += chunk
        if return_blocks and index % 5 == 0:
            validation_blocks.append(chunk)
        predicted = valid & (called != parent)
        truth = child != parent
        audit["true_mutations"] += int(truth.sum())
        audit["called_mutations"] += int(predicted.sum())
        audit["correct_mutations"] += int((predicted & truth & (called == child)).sum())
        audit["callable_sites"] += int(valid.sum())
        audit["masked_sites"] += int(len(parent)-6-valid.sum())
        audit["alignment_edit_distance"] += int(distance)
    result = (train, validation, audit)
    return (*result, validation_blocks) if return_blocks else result


def reduce_counts(counts, k):
    out = np.zeros((4**k, 4))
    np.add.at(out, project(np.arange(4**7), k), counts)
    return out


def pooled_counts(counts, k):
    return counts + counts[rc_indices(k)][:, COMPLEMENT]


def fit_models(train, validation, return_training_models=False):
    fitted, chosen, validation_nll = {}, {}, {}
    training_models = {}
    parent_train = parent_full = None
    for k in (1, 3, 5, 7):
        raw_train = reduce_counts(train, k)
        raw_full = reduce_counts(train + validation, k)
        cv = reduce_counts(validation, k)
        tr, full = pooled_counts(raw_train, k), pooled_counts(raw_full, k)
        if k == 1:
            # Fixed weak Dirichlet prior prevents zero probabilities.
            pred_train = (tr + 0.5) / (tr.sum(axis=1, keepdims=True) + 2)
            pred_full = (full + 0.5) / (full.sum(axis=1, keepdims=True) + 2)
            base_train = pred_train
            chosen["1mer"] = 2.0
        else:
            parent_map = (np.arange(4**k) // 4) % (4**(k-2))
            prior = parent_train[parent_map]
            candidates = [(tr + s * prior) / (tr.sum(axis=1, keepdims=True) + s) for s in SHRINKAGE]
            losses = [-float(xlogy(cv, p).sum()) for p in candidates]
            best = int(np.argmin(losses))
            s = SHRINKAGE[best]
            pred_train = candidates[best]
            pred_full = (full + s * parent_full[parent_map]) / (full.sum(axis=1, keepdims=True) + s)
            chosen[f"{k}mer"] = s
        validation_nll[f"{k}mer"] = -float(xlogy(cv, pred_train).sum())
        fitted[f"{k}mer"] = pred_full
        training_models[f"{k}mer"] = pred_train
        parent_train, parent_full = pred_train, pred_full
    # Independently regularized 7mer model: the same CV grid, but its prior is
    # the single-base model. The main 7mer comparator can also pool through 5mer.
    centers = digits(7)[:, 3]
    tr = pooled_counts(train, 7)
    full = pooled_counts(train + validation, 7)
    candidates = [(tr + s*base_train[centers]) / (tr.sum(axis=1, keepdims=True)+s) for s in SHRINKAGE]
    losses = [-float(xlogy(validation, p).sum()) for p in candidates]
    best = int(np.argmin(losses))
    s = SHRINKAGE[best]
    fitted["7mer_flat"] = (full + s * fitted["1mer"][centers]) / (full.sum(axis=1, keepdims=True) + s)
    chosen["7mer_flat"] = s
    validation_nll["7mer_flat"] = losses[best]
    training_models["7mer_flat"] = candidates[best]
    result = (fitted, chosen, validation_nll)
    return (*result, training_models) if return_training_models else result


def fit_one(job):
    root, scenario, seed, _config = job
    directory = Path(root) / scenario / f"seed_{seed:03d}"
    if not manifest_valid(directory, "samples.manifest.json"):
        raise ValueError(f"Corrupt or missing samples: {directory}; rerun prepare")
    sample_manifest = json.loads((directory / "samples.manifest.json").read_text())
    # Only public inputs determine fit cache. Neither hidden outcomes nor kernel enters fitting.
    cfg = {"schema": SCHEMA, "fitting_version": FIT_VERSION,
           "inputs": {k: v for k, v in sample_manifest["files"].items() if k.startswith("public/")},
           "shrinkage_grid": list(SHRINKAGE), "validation": "contig index modulo 5 == 0"}
    if manifest_valid(directory, "fit.manifest.json", cfg):
        return f"cached fit {scenario} seed={seed}"
    if manifest_valid(directory, "fit.manifest.json") and json.loads((directory / "fit.manifest.json").read_text())["config"]["inputs"] == cfg["inputs"]:
        with np.load(directory / "fit" / "counts.npz") as counts:
            train, validation = counts["train"], counts["validation"]
        audit = json.loads((directory / "fit" / "audit.json").read_text())["alignment"]
    else:
        train, validation, audit = collect(directory / "public")
    fitted, chosen, validation_nll = fit_models(train, validation)
    fit_dir = directory / "fit"
    fit_dir.mkdir(exist_ok=True)
    np.savez_compressed(fit_dir / "models.npz", **fitted)
    np.savez_compressed(fit_dir / "counts.npz", train=train, validation=validation)
    audit["precision"] = audit["correct_mutations"] / max(1, audit["called_mutations"])
    audit["recall"] = audit["correct_mutations"] / max(1, audit["true_mutations"])
    write_json(fit_dir / "audit.json", {"alignment": audit, "shrinkage": chosen,
                                        "validation_nll": validation_nll,
                                        "training_sites": int(train.sum()), "validation_sites": int(validation.sum())})
    seal(directory, "fit.manifest.json", cfg, sorted(fit_dir.glob("*")))
    return f"fitted {scenario} seed={seed} precision={audit['precision']:.6f} recall={audit['recall']:.6f} shrink={chosen}"


def score_counts(opportunities, outcomes, truth, prediction, ancestral=None):
    anc = digits(7)[:, 3] if ancestral is None else np.asarray(ancestral)
    if truth.shape != (len(anc), 4) or prediction.shape != truth.shape or outcomes.shape != truth.shape:
        raise ValueError("Scoring arrays must agree with the ancestral-base encoding")
    off = np.ones_like(truth, dtype=bool)
    off[np.arange(len(anc)), anc] = False
    true_alt, pred_alt = np.where(off, truth, 0), np.where(off, prediction, 0)
    true_mu, pred_mu = true_alt.sum(axis=1), pred_alt.sum(axis=1)
    sites = opportunities.sum()
    expected = float(opportunities @ true_mu)
    observed = float(np.where(off, outcomes, 0).sum())
    tc, pc = true_alt / true_mu[:, None], pred_alt / pred_mu[:, None]
    kl = (xlogy(truth, truth) - xlogy(truth, prediction)).sum(axis=1)
    ck = (xlogy(tc, tc) - xlogy(tc, pc)).sum(axis=1)
    weak = (anc == 0) | (anc == 3)
    ws_mask = weak[:, None] & np.isin(np.arange(4), [1, 2])[None, :]
    sw_mask = ~weak[:, None] & np.isin(np.arange(4), [0, 3])[None, :]
    true_ws = float((opportunities[:, None] * truth * ws_mask).sum())
    true_sw = float((opportunities[:, None] * truth * sw_mask).sum())
    pred_ws = float((opportunities[:, None] * prediction * ws_mask).sum())
    pred_sw = float((opportunities[:, None] * prediction * sw_mask).sum())
    conditional_nll = -float(xlogy(np.where(off, outcomes, 0), pc).sum() / max(observed, 1))
    full_nll = -float(xlogy(outcomes, prediction).sum() / sites)
    return {"sites": int(sites), "observed_mutations": int(observed), "expected_mutations": expected,
            "predicted_mutations": float(opportunities @ pred_mu),
            "burden_relative_error": float(opportunities @ pred_mu / expected - 1),
            "gc_odds_log_error": float(np.log((pred_ws / pred_sw) / (true_ws / true_sw))),
            "oracle_kl_per_mutation": float(opportunities @ kl / expected),
            "oracle_conditional_kl": float((opportunities * true_mu) @ ck / expected),
            "oracle_binary_kl_per_mutation": float((opportunities @ kl - (opportunities*true_mu) @ ck) / expected),
            "observed_nll_per_site": full_nll,
            "observed_binary_nll_per_site": full_nll - conditional_nll*observed/sites,
            "observed_conditional_nll": conditional_nll,
            "oracle_brier_excess": float(opportunities @ ((prediction-truth)**2).sum(axis=1) / sites)}


def evaluate_one(job):
    root, scenario, seed, _config = job
    directory = Path(root) / scenario / f"seed_{seed:03d}"
    for filename in ("samples.manifest.json", "fit.manifest.json"):
        if not manifest_valid(directory, filename):
            raise ValueError(f"Missing/corrupt {filename}: {directory}")
    cfg = {"schema": SCHEMA, "evaluation_version": EVALUATION_VERSION, "samples": sha256(directory / "samples.manifest.json"),
           "fit": sha256(directory / "fit.manifest.json")}
    if manifest_valid(directory, "evaluation.manifest.json", cfg):
        return f"cached evaluation {scenario} seed={seed}"
    parents = read_fasta(directory / "public" / "new_parent.fasta")
    children = read_fasta(directory / "private" / "new_offspring.fasta")
    outcomes = np.zeros((4**7, 4))
    for name, parent in parents.items():
        # Evaluator knows synthetic correspondence; this does not enter fitting.
        outcomes += np.bincount(contexts(parent)*4 + children[name][3:-3], minlength=4**7*4).reshape(-1, 4)
    opportunities = outcomes.sum(axis=1)
    truth = np.load(directory / "private" / "kernel.npy")
    audit = json.loads((directory / "fit" / "audit.json").read_text())
    training_counts = np.load(directory / "fit" / "counts.npz")
    train_opp3 = reduce_counts(training_counts["train"] + training_counts["validation"], 3).sum(axis=1)
    test_opp3 = reduce_counts(outcomes, 3).sum(axis=1)
    tri_tv = float(np.abs(train_opp3/train_opp3.sum()-test_opp3/test_opp3.sum()).sum()/2)
    records = []
    with np.load(directory / "fit" / "models.npz") as fitted:
        for method in METHODS:
            k = int(method[0])
            prediction = fitted[method][project(np.arange(4**7), k)]
            records.append({"scenario": scenario, "seed": seed, "method": method,
                            "trinucleotide_tv": tri_tv, "alignment_precision": audit["alignment"]["precision"],
                            "alignment_recall": audit["alignment"]["recall"],
                            "shrinkage": audit["shrinkage"].get(method, 2),
                            **score_counts(opportunities, outcomes, truth, prediction)})
    write_json(directory / "scores.json", records)
    seal(directory, "evaluation.manifest.json", cfg, [directory / "scores.json"])
    return f"evaluated {scenario} seed={seed}"


def export_predictions(directory, method, output):
    directory, output = Path(directory), Path(output)
    with np.load(directory / "fit" / "models.npz") as models:
        model = models[method]
    k = int(method[0])
    parents = read_fasta(directory / "public" / "new_parent.fasta")
    lookup = {name: model[project(contexts(seq), k)] for name, seq in parents.items()}
    expected, rows = 0.0, 0
    output.parent.mkdir(parents=True, exist_ok=True)
    with (directory / "public" / "queries.tsv").open() as source, output.open("w", newline="") as target:
        reader = csv.DictReader(source, delimiter="\t")
        writer = csv.writer(target, delimiter="\t", lineterminator="\n")
        writer.writerow(["contig", "position", "p_A", "p_C", "p_G", "p_T"])
        for row in reader:
            name, pos = row["contig"], int(row["position"])
            p = lookup[name][pos-4]
            assert np.isfinite(p).all() and (p >= 0).all() and np.isclose(p.sum(), 1)
            expected += 1 - p[parents[name][pos-1]]
            writer.writerow([name, pos, *[f"{value:.12g}" for value in p]])
            rows += 1
    write_json(output.with_suffix(".summary.json"), {"method": method, "query_count": rows,
                                                     "expected_changed_alleles": expected,
                                                     "prediction_sha256": sha256(output)})


def report(root, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    matplotlib.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 16,
                                "axes.labelsize": 17, "axes.titlesize": 18, "legend.fontsize": 13})
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    rows = []
    for path in sorted(Path(root).glob("*/seed_*/scores.json")):
        if not manifest_valid(path.parent, "evaluation.manifest.json"):
            raise ValueError(f"Corrupt evaluation: {path}")
        rows.extend(json.loads(path.read_text()))
    if not rows:
        raise ValueError("No scored runs")
    with (output / "replicate_scores.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    summaries, paired = [], []
    scenarios = [s for s in SCENARIOS if any(r["scenario"] == s for r in rows)]
    metrics = ("oracle_kl_per_mutation", "oracle_conditional_kl", "oracle_binary_kl_per_mutation", "observed_conditional_nll",
               "burden_relative_error", "gc_odds_log_error", "observed_nll_per_site", "observed_binary_nll_per_site", "oracle_brier_excess")
    rng = np.random.default_rng(20261112)
    for scenario in scenarios:
        for method in METHODS:
            selected = [r for r in rows if r["scenario"] == scenario and r["method"] == method]
            summaries.append({"scenario": scenario, "method": method, "replicates": len(selected),
                              **{m: float(np.mean([r[m] for r in selected])) for m in metrics},
                              "burden_mean_absolute_relative_error": float(np.mean([abs(r["burden_relative_error"]) for r in selected]))})
        for a, b in (("3mer", "5mer"), ("5mer", "7mer"), ("5mer", "7mer_flat")):
            first = {r["seed"]: r for r in rows if r["scenario"] == scenario and r["method"] == a}
            second = {r["seed"]: r for r in rows if r["scenario"] == scenario and r["method"] == b}
            for metric in ("oracle_kl_per_mutation", "observed_conditional_nll"):
                delta = np.array([second[s][metric]-first[s][metric] for s in sorted(first)])
                boot = rng.choice(delta, size=(10000, len(delta))).mean(axis=1)
                paired.append({"scenario": scenario, "contrast": f"{b} minus {a}", "metric": metric,
                               "mean_difference": float(delta.mean()), "ci_low": float(np.quantile(boot, .025)),
                               "ci_high": float(np.quantile(boot, .975)), "second_better_replicates": int((delta < 0).sum()),
                               "replicates": len(delta)})
    write_json(output / "summary.json", summaries)
    write_json(output / "paired_comparisons.json", paired)
    labels = {"constant": "1-base\ntruth", "three": "3-base\ntruth", "five_iid": "5-base\nsame mix",
              "five_shift": "5-base\nshifted mix", "seven_shift": "7-base\ntruth", "empirical_five": "Published\n5-base kernel"}
    fig, axes = plt.subplots(2, 1, figsize=(11, 8.5), sharex=True, layout="constrained")
    colors = ["#777777", "#CC6677", "#228833", "#4477AA", "#AA3377"]
    for j, method in enumerate(METHODS):
        x = np.arange(len(scenarios)) + (j-2)*.11
        for ax, metric in zip(axes, ("oracle_kl_per_mutation", "burden_mean_absolute_relative_error")):
            y = [next(s[metric] for s in summaries if s["scenario"] == scenario and s["method"] == method) for scenario in scenarios]
            ax.plot(x, y, "o", color=colors[j], markersize=7, label=method.replace("mer_flat", "mer flat"))
    axes[0].set_yscale("symlog", linthresh=.002)
    axes[0].set_ylabel("Excess log loss\n(nats / expected mutation)")
    axes[0].set_title("Allele probability prediction from assembled genomes", pad=47)
    axes[0].legend(loc="lower center", bbox_to_anchor=(.5, 1.02), ncol=5, frameon=False)
    axes[1].set_ylabel("Mutation burden\nmean absolute relative error")
    axes[1].set_xticks(np.arange(len(scenarios)), [labels[s] for s in scenarios], fontsize=13)
    for ax in axes:
        ax.grid(axis="y", alpha=.2)
        ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(output / "prediction_comparison.png", dpi=180)
    fig.savefig(output / "prediction_comparison.pdf")
    plt.close(fig)
    write_json(output / "report_manifest.json", {"source_scores": {str(p.resolve()): sha256(p) for p in sorted(Path(root).glob("*/seed_*/scores.json"))},
                                                 "bootstrap_seed": 20261112, "bootstrap_replicates": 10000,
                                                 "files": {p.name: sha256(p) for p in output.iterdir() if p.is_file() and p.name != "report_manifest.json"}})
    print(json.dumps(summaries, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=("prepare", "fit", "evaluate", "report", "export", "all"))
    parser.add_argument("--root", type=Path, default=Path("work/main"))
    parser.add_argument("--output", type=Path, default=Path("results"))
    parser.add_argument("--replicates", type=int, default=12)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--memory-gb", type=float, default=26, help="Logged resource budget; this small-contig pipeline is bounded by workers")
    parser.add_argument("--train-bases", type=int, default=4_000_000)
    parser.add_argument("--test-bases", type=int, default=1_000_000)
    parser.add_argument("--contig-bases", type=int, default=50_000)
    parser.add_argument("--mean-rate", type=float, default=.005)
    parser.add_argument("--target-outer-match", type=float, default=.8)
    parser.add_argument("--scenarios", nargs="+", choices=SCENARIOS, default=list(SCENARIOS[:-1]))
    parser.add_argument("--task-dir", type=Path)
    parser.add_argument("--method", choices=METHODS, default="5mer")
    args = parser.parse_args()
    if not 1 <= args.workers <= 24:
        parser.error("Local workers must be 1..24")
    if args.train_bases < 10*args.contig_bases or args.contig_bases < 20 or args.test_bases < args.contig_bases:
        parser.error("Use at least 10 training contigs, >=1 target contig, and >=20 bases/contig")
    if args.train_bases % args.contig_bases or args.test_bases % args.contig_bases:
        parser.error("Base counts must be divisible by contig size")
    if not 0 < args.mean_rate < .1 or not 0 < args.target_outer_match < 1 or args.replicates < 1:
        parser.error("Invalid probability or replicate count")
    args.root.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                        handlers=[logging.StreamHandler(), logging.FileHandler(args.root / "run.log", encoding="utf-8")])
    logging.getLogger("fontTools").setLevel(logging.WARNING)
    revision = subprocess.run(["git", "rev-parse", "HEAD"], text=True, capture_output=True).stdout.strip()
    LOG.info("parameters=%s cwd=%s python=%s platform=%s git=%s script_sha256=%s", vars(args), Path.cwd(), sys.version, platform.platform(), revision, sha256(__file__))
    config = {k: getattr(args, k) for k in ("train_bases", "test_bases", "contig_bases", "mean_rate", "target_outer_match")}
    jobs = [(str(args.root.resolve()), scenario, seed, config) for scenario in args.scenarios for seed in range(args.replicates)]
    write_json(args.root / f"invocation_{time.time_ns()}.json", {"arguments": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
                                                               "config": config, "git": revision, "script_sha256": sha256(__file__),
                                                               "python": sys.version, "seeds": list(range(args.replicates)),
                                                               "dependency_lock_sha256": sha256(Path(__file__).parent / "uv.lock"),
                                                               "packages": {name: importlib.metadata.version(name) for name in ("numpy", "scipy", "edlib", "matplotlib")}})
    phases = ("prepare", "fit", "evaluate", "report") if args.phase == "all" else (args.phase,)
    for phase in phases:
        if phase in ("prepare", "fit", "evaluate"):
            function = {"prepare": prepare_one, "fit": fit_one, "evaluate": evaluate_one}[phase]
            with ProcessPoolExecutor(max_workers=args.workers) as pool:
                for result in pool.map(function, jobs):
                    LOG.info(result)
        elif phase == "report":
            report(args.root, args.output)
        elif phase == "export":
            if args.task_dir is None:
                parser.error("export needs --task-dir")
            export_predictions(args.task_dir, args.method, args.output)


if __name__ == "__main__":
    main()
