"""Evaluator-only grading of a probability TSV against withheld offspring and truth."""
import argparse
import csv
import json
from itertools import zip_longest
from pathlib import Path

import numpy as np
from scipy.special import xlogy

import benchmark as b


def grade(task, predictions, output):
    task = Path(task)
    parents = b.read_fasta(task / "public" / "new_parent.fasta")
    children = b.read_fasta(task / "private" / "new_offspring.fasta")
    kernel = np.load(task / "private" / "kernel.npy")
    circular = json.loads((task/"public"/"metadata.json").read_text()).get("topology") == "circular"
    ids = {name: b.contexts(np.pad(seq, (3, 3), mode="wrap") if circular else seq) for name, seq in parents.items()}
    totals = {"sites": 0, "observed_mutations": 0, "expected_mutations": 0., "predicted_mutations": 0.,
              "observed_nll": 0., "observed_conditional_nll_total": 0., "oracle_kl": 0., "oracle_conditional_kl_total": 0.}
    # Bound memory independently of query count.
    batch = []

    def flush():
        if not batch:
            return
        ancestor = np.array([x[0] for x in batch])
        descendant = np.array([x[1] for x in batch])
        truth = kernel[np.array([x[2] for x in batch])]
        p = np.array([x[3] for x in batch])
        if not np.isfinite(p).all() or (p < 0).any() or (p > 1).any() or not np.allclose(p.sum(axis=1), 1, atol=1e-9, rtol=0):
            raise ValueError("Probabilities must be finite, in [0,1], and sum to one within 1e-9")
        n = len(batch)
        i = np.arange(n)
        changed = ancestor != descendant
        true_mu, pred_mu = 1-truth[i, ancestor], 1-p[i, ancestor]
        # Zero probabilities at outcomes with positive mass imply infinite log loss.
        # Emit an explicit failed score instead of silently clipping the submission.
        if (p <= 0).any():
            raise ValueError("Submission has infinite oracle log loss: all four outcomes have positive probability in this task")
        pc, tc = p.copy(), truth.copy()
        pc[i, ancestor] = tc[i, ancestor] = 0
        pc /= pred_mu[:, None]
        tc /= true_mu[:, None]
        totals["sites"] += n
        totals["observed_mutations"] += int(changed.sum())
        totals["expected_mutations"] += float(true_mu.sum())
        totals["predicted_mutations"] += float(pred_mu.sum())
        totals["observed_nll"] -= float(np.log(p[i, descendant]).sum())
        totals["observed_conditional_nll_total"] -= float(np.log(pc[i[changed], descendant[changed]]).sum())
        totals["oracle_kl"] += float((xlogy(truth, truth)-xlogy(truth, p)).sum())
        totals["oracle_conditional_kl_total"] += float((true_mu * (xlogy(tc, tc)-xlogy(tc, pc)).sum(axis=1)).sum())
        batch.clear()

    with (task / "public" / "queries.tsv").open() as queries, Path(predictions).open() as submitted:
        expected = csv.DictReader(queries, delimiter="\t")
        actual = csv.DictReader(submitted, delimiter="\t")
        if actual.fieldnames != ["contig", "position", "p_A", "p_C", "p_G", "p_T"]:
            raise ValueError("Incorrect prediction columns")
        for query, row in zip_longest(expected, actual):
            if query is None or row is None or row["contig"] != query["contig"] or row["position"] != query["position"]:
                raise ValueError("Predictions must match every query exactly, in order")
            name, pos = query["contig"], int(query["position"])-1
            batch.append((parents[name][pos], children[name][pos], ids[name][pos if circular else pos-3],
                          [float(row[f"p_{base}"]) for base in "ACGT"]))
            if len(batch) == 50000:
                flush()
        flush()
    totals.update({"oracle_kl_per_mutation": totals["oracle_kl"] / totals["expected_mutations"],
                   "oracle_conditional_kl": totals["oracle_conditional_kl_total"] / totals["expected_mutations"],
                   "observed_nll_per_site": totals["observed_nll"] / totals["sites"],
                   "observed_conditional_nll": totals["observed_conditional_nll_total"] / max(1, totals["observed_mutations"]),
                   "burden_relative_error": totals["predicted_mutations"] / totals["expected_mutations"] - 1,
                   "prediction_sha256": b.sha256(predictions)})
    b.write_json(output, totals)
    return totals


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("task", type=Path)
    parser.add_argument("predictions", type=Path)
    parser.add_argument("--output", type=Path, default=Path("grade.json"))
    args = parser.parse_args()
    print(grade(args.task, args.predictions, args.output))
