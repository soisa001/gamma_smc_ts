"""Author-side directional window comparators; no strand symmetry is assumed."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import xlogy

import benchmark as b


WINDOWS = ((0, 0), (1, 1), (2, 2), (3, 2), (2, 3), (3, 3))


def project(ids, window):
    left, right = window
    return (ids // 4**(3-right)) % 4**(left+right+1)


def aggregate(counts, window):
    out = np.zeros((4**(sum(window)+1), 4))
    np.add.at(out, project(np.arange(4**7), window), counts)
    return out


def run(task, output):
    task = Path(task)
    if not b.manifest_valid(task, "samples.manifest.json"):
        raise ValueError("Corrupt task inputs")
    circular = json.loads((task/"public"/"metadata.json").read_text()).get("topology") == "circular"
    if circular:
        # The public task guarantees homologous SNV-only sequences. Use the
        # corresponding ungapped alignment at every site, including wraparound.
        parents = b.read_fasta(task/"public"/"calibration_parent.fasta")
        children = b.read_fasta(task/"public"/"calibration_offspring.fasta")
        train, validation = np.zeros((4**7, 4)), np.zeros((4**7, 4))
        changed = 0
        for i, (name, seq) in enumerate(parents.items()):
            if len(seq) != len(children[name]):
                raise ValueError("Unequal lengths violate SNV-only correspondence")
            ids = b.contexts(np.pad(seq, (3, 3), mode="wrap"))
            counts = np.bincount(ids*4+children[name], minlength=4**7*4).reshape(-1, 4)
            (validation if i % 5 == 0 else train)[:] += counts
            changed += int((seq != children[name]).sum())
        audit = {"method": "ungapped SNV-only alignment with circular contexts",
                 "callable_sites": int(train.sum()+validation.sum()), "called_mutations": changed,
                 "masked_sites": 0}
    else:
        train, validation, audit = b.collect(task / "public")
    fitted, partial, cv_losses, smoothing = {}, {}, {}, {}
    for window in WINDOWS:
        tr, full, val = [aggregate(counts, window) for counts in (train, train+validation, validation)]
        if window == (0, 0):
            partial[window] = (tr+.5)/(tr.sum(axis=1, keepdims=True)+2)
            fitted[window] = (full+.5)/(full.sum(axis=1, keepdims=True)+2)
            smoothing[window] = 2.
        else:
            if window == (3, 3):
                parent = min(((3, 2), (2, 3)), key=cv_losses.get)
            elif sum(window) == 5:
                parent = (2, 2)
            else:
                parent = (window[0]-1, window[1]-1)
            parent_index = (np.arange(len(tr)) // 4**(window[1]-parent[1])) % len(partial[parent])
            options = [(tr+s*partial[parent][parent_index])/(tr.sum(axis=1, keepdims=True)+s) for s in b.SHRINKAGE]
            best = int(np.argmin([-xlogy(val, p).sum() for p in options]))
            strength = b.SHRINKAGE[best]
            partial[window] = options[best]
            fitted[window] = (full+strength*fitted[parent][parent_index])/(full.sum(axis=1, keepdims=True)+strength)
            smoothing[window] = strength
        cv_losses[window] = -float(xlogy(val, partial[window]).sum())
    # Only the following evaluator section accesses withheld target outcomes/truth.
    parents = b.read_fasta(task / "public" / "new_parent.fasta")
    children = b.read_fasta(task / "private" / "new_offspring.fasta")
    counts = np.zeros((4**7, 4))
    for name, seq in parents.items():
        ids = b.contexts(np.pad(seq, (3, 3), mode="wrap") if circular else seq)
        outcomes = children[name] if circular else children[name][3:-3]
        counts += np.bincount(ids*4+outcomes, minlength=4**7*4).reshape(-1, 4)
    truth = np.load(task / "private" / "kernel.npy")
    rows = []
    for window in WINDOWS:
        pred = fitted[window][project(np.arange(4**7), window)]
        rows.append({"window": list(window), "label": f"{sum(window)+1}-mer (-{window[0]}..+{window[1]})",
                     "validation_nll": cv_losses[window], "shrinkage": smoothing[window],
                     **b.score_counts(counts.sum(axis=1), counts, truth, pred)})
    result = {"alignment_audit": audit, "strand_pooling": False, "target_topology": "circular" if circular else "linear",
              "calibration_margin": 0 if circular else 3, "scores": rows,
              "samples_manifest_sha256": b.sha256(task / "samples.manifest.json"),
              "source_sha256": b.sha256(__file__)}
    b.write_json(output, result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("task", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(run(args.task, args.output))
