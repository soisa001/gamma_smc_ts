"""Author-only arbitrary directed six-base spectrum with clue-free query boundaries."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
from pathlib import Path

import numpy as np

import benchmark as b


def rng_for(seed, *key):
    return np.random.default_rng(np.random.SeedSequence(seed, spawn_key=key))


def kernel(seed=0, mean_rate=.005, log_sd=1.25):
    words = b.digits(6)
    ancestral = words[:, 3]
    # Each applicable channel is an independent draw before one common scaling.
    rates = rng_for(seed, 0).lognormal(-log_sd**2/2, log_sd, size=(4**6, 4))
    rates[np.arange(4**6), ancestral] = 0
    rates *= mean_rate/rates.sum(axis=1).mean()
    if rates.sum(axis=1).max() >= .5:
        raise ValueError("Mutation probabilities too large; lower the declared mean rate")
    rates[np.arange(4**6), ancestral] = 1-rates.sum(axis=1)
    return np.repeat(rates, 4, axis=0)  # +3 is irrelevant in the grader's 7-base encoding.


def circular_contexts(sequence):
    return b.contexts(np.pad(sequence, (3, 3), mode="wrap"))


def generate_contig(job):
    seed, group, index, length, truth = job
    rng = rng_for(seed, 1, group, index)
    parent = rng.integers(0, 4, size=length, dtype=np.int8)
    probs = truth[circular_contexts(parent)]
    child = (rng.random(length)[:, None] > np.cumsum(probs, axis=1)).sum(axis=1).astype(np.int8)
    return f"scaffold_{index:04d}", parent, child


def prepare(directory, seed=0, train_bases=4_000_000, test_bases=1_000_000,
            contig_bases=50_000, mean_rate=.005, log_sd=1.25, workers=4):
    directory = Path(directory)
    cfg = {"generator": "arbitrary_directed_six_v1", "seed": seed, "train_bases": train_bases,
           "test_bases": test_bases, "contig_bases": contig_bases, "mean_rate": mean_rate,
           "log_sd": log_sd, "topology": "circular", "query_policy": "all positions, independently shuffled",
           "offsets": [-3, -2, -1, 0, 1, 2], "strand_symmetry": False,
           "streams": "SeedSequence(seed,spawn_key): kernel=(0,), contig=(1,group,index), query=(2,)"}
    if b.manifest_valid(directory, "samples.manifest.json", cfg):
        return {"cached": True, "directory": str(directory)}
    if train_bases % contig_bases or test_bases % contig_bases or contig_bases < 20:
        raise ValueError("Use complete contigs of at least 20 bases")
    if not 1 <= workers <= 24:
        raise ValueError("Use 1..24 workers")
    public, private = directory/"public", directory/"private"
    public.mkdir(parents=True, exist_ok=True)
    private.mkdir(parents=True, exist_ok=True)
    truth = kernel(seed, mean_rate, log_sd)
    for group, (label, total) in enumerate((("calibration", train_bases), ("new", test_bases))):
        jobs = [(seed, group, i, contig_bases, truth) for i in range(total//contig_bases)]
        if workers == 1:
            records = list(map(generate_contig, jobs))
        else:
            with ProcessPoolExecutor(max_workers=workers) as pool:
                records = list(pool.map(generate_contig, jobs))
        b.write_fasta(public/f"{label}_parent.fasta", [(n, p) for n, p, _ in records])
        b.write_fasta((public if group == 0 else private)/f"{label}_offspring.fasta", [(n, c) for n, _, c in records])
    order = rng_for(seed, 2).permutation(test_bases)
    with (public/"queries.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["contig", "position"])
        writer.writerows((f"scaffold_{int(i)//contig_bases:04d}", int(i)%contig_bases+1) for i in order)
    prompt = b.TASK_PROMPT.replace("Only\nsingle-nucleotide changes are generated.",
        "Each scaffold represents a circular DNA molecule; coordinates wrap from its\nfinal base to its first base. Only single-nucleotide changes are generated.")
    assert prompt != b.TASK_PROMPT
    (public/"TASK.md").write_text(prompt, encoding="utf-8")
    b.write_json(public/"metadata.json", {"transmissions_per_pair": 1, "coordinates": "1-based",
                 "ploidy": "haploid transmitted allele", "assembly_error_rate": 0, "topology": "circular"})
    np.save(private/"kernel.npy", truth)
    b.seal(directory, "samples.manifest.json", cfg, sorted(public.glob("*"))+sorted(private.glob("*")))
    b.write_json(directory/"generation_provenance.json", {"config": cfg, "workers": workers,
                 "memory_budget_gb": 26, "source_sha256": b.sha256(__file__),
                 "max_mutation_probability": float(1-truth[np.arange(4**7), b.digits(7)[:, 3]].min()),
                 "numpy_version": np.__version__})
    return {"cached": False, "directory": str(directory), "config": cfg}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    print(prepare(args.directory, seed=args.seed, workers=args.workers))
