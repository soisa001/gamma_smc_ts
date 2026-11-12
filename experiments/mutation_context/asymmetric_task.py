"""Author-only generator for a directed six-base mock process (offsets -3 through +2)."""
import argparse
import csv
from pathlib import Path

import numpy as np

import benchmark as b


def kernel(mean_rate=.005):
    d = b.digits(7)
    ancestral = d[:, 3]
    rates = np.ones((4**7, 4))
    rates[np.arange(4**7), ancestral] = 0
    rates[np.arange(4**7), ancestral ^ 2] *= 2
    rates[(ancestral == 1) & (d[:, 4] == 2), 3] *= 6
    rates[(ancestral == 2) & (d[:, 2] == 1), 0] *= 6
    # Long-range interaction needs position -3, outside a centered five-base window.
    outer = d[:, 0] == d[:, 5]
    weak = np.isin(ancestral, [0, 3])
    rates[outer & weak, 1] *= 4
    rates[outer & weak, 2] *= 4
    # Give the remaining leading base (-2) a genuine, separately learnable role.
    adjacent = d[:, 1] == d[:, 2]
    rates[adjacent & (ancestral == 0), 3] *= 2
    rates[adjacent & (ancestral == 3), 0] *= 2
    rates *= mean_rate/rates.sum(axis=1).mean()
    assert rates.sum(axis=1).max() < .5
    rates[np.arange(4**7), ancestral] = 1-rates.sum(axis=1)
    assert np.allclose(rates.reshape(4**6, 4, 4), rates.reshape(4**6, 4, 4)[:, :1])
    return rates


def sequence(rng, length, match_probability):
    """Uniform consecutive five-base words, tunable dependence at lag five."""
    increments = rng.integers(1, 4, size=((length+4)//5, 5))
    increments[rng.random(increments.shape) < match_probability] = 0
    increments[0] = rng.integers(0, 4, size=5)
    return (np.cumsum(increments, axis=0) % 4).astype(np.int8).reshape(-1)[:length]


def prepare(directory, seed=0, train_bases=4_000_000, test_bases=1_000_000, contig_bases=50_000, mean_rate=.005):
    directory = Path(directory)
    cfg = {"generator": "directed_six_base_v1", "seed": seed, "train_bases": train_bases,
           "test_bases": test_bases, "contig_bases": contig_bases, "mean_rate": mean_rate,
           "training_match_probability": .25, "target_match_probability": .8,
           "transmissions": 1, "offsets": [-3, -2, -1, 0, 1, 2], "strand_symmetry": False}
    if b.manifest_valid(directory, "samples.manifest.json", cfg):
        return {"cached": True, "directory": str(directory)}
    if train_bases % contig_bases or test_bases % contig_bases or contig_bases < 20:
        raise ValueError("Use complete contigs of at least 20 bases")
    public, private = directory / "public", directory / "private"
    public.mkdir(parents=True, exist_ok=True)
    private.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    truth = kernel(mean_rate)
    for group, total, q in (("calibration", train_bases, .25), ("new", test_bases, .8)):
        parents, children = [], []
        for index in range(total // contig_bases):
            parent = sequence(rng, contig_bases, q)
            name = f"scaffold_{index:04d}"
            parents.append((name, parent))
            children.append((name, b.mutate(rng, parent, truth)))
        b.write_fasta(public / f"{group}_parent.fasta", parents)
        b.write_fasta((public if group == "calibration" else private) / f"{group}_offspring.fasta", children)
        if group == "new":
            with (public / "queries.tsv").open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["contig", "position"])
                for name, parent in parents:
                    writer.writerows((name, pos+1) for pos in range(3, len(parent)-3))
    (public / "TASK.md").write_text(b.TASK_PROMPT, encoding="utf-8")
    b.write_json(public / "metadata.json", {"transmissions_per_pair": 1, "coordinates": "1-based",
                                            "ploidy": "haploid transmitted allele", "assembly_error_rate": 0})
    np.save(private / "kernel.npy", truth)
    b.seal(directory, "samples.manifest.json", cfg, sorted(public.glob("*"))+sorted(private.glob("*")))
    return {"cached": False, "directory": str(directory), "config": cfg}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    print(prepare(args.directory, seed=args.seed))
