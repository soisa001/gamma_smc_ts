"""Reference solver that reads only the public task directory."""
import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.special import xlogy

import benchmark as b


def solve(public, output, method="auto"):
    public, output = Path(public), Path(output)
    train, validation, audit, blocks = b.collect(public, return_blocks=True)
    models, shrinkage, losses, training_models = b.fit_models(train, validation, return_training_models=True)
    selection = {}
    if method == "auto":
        # Paired block one-standard-error rule: prefer a smaller window whose
        # validation disadvantage is within its contig-level sampling error.
        best = min(losses, key=losses.get)
        p_best = training_models[best][b.project(np.arange(4**7), int(best[0]))]
        for candidate in b.METHODS:
            p = training_models[candidate][b.project(np.arange(4**7), int(candidate[0]))]
            differences = np.array([float((xlogy(block, p_best)-xlogy(block, p)).sum()) for block in blocks])
            se = float(differences.std(ddof=1)/np.sqrt(len(blocks))) if len(blocks) > 1 else 0.0
            eligible = bool(float(differences.mean()) <= se + 1e-8*validation.sum()/len(blocks))
            selection[candidate] = {"mean_extra_nll_per_contig": float(differences.mean()),
                                    "paired_standard_error": se, "eligible": eligible}
        method = next(name for name in b.METHODS if selection[name]["eligible"])
    model = models[method]
    parents = b.read_fasta(public / "new_parent.fasta")
    lookup = {name: model[b.project(b.contexts(seq), int(method[0]))] for name, seq in parents.items()}
    expected, rows = 0.0, 0
    output.parent.mkdir(parents=True, exist_ok=True)
    with (public / "queries.tsv").open() as source, output.open("w", newline="") as target:
        reader = csv.DictReader(source, delimiter="\t")
        writer = csv.writer(target, delimiter="\t", lineterminator="\n")
        writer.writerow(["contig", "position", "p_A", "p_C", "p_G", "p_T"])
        for row in reader:
            name, position = row["contig"], int(row["position"])
            p = lookup[name][position-4]
            if not np.isfinite(p).all() or (p < 0).any() or not np.isclose(p.sum(), 1):
                raise ValueError("Invalid fitted probability vector")
            expected += 1-p[parents[name][position-1]]
            writer.writerow([name, position, *[f"{x:.12g}" for x in p]])
            rows += 1
    summary = {"selected_method": method, "selection_rule": "paired contig one-standard-error, smallest window first",
               "selection_diagnostics": selection, "shrinkage": shrinkage, "validation_nll": losses,
               "alignment": audit, "expected_changed_alleles": expected, "rows": rows,
               "prediction_sha256": b.sha256(output)}
    b.write_json(output.with_suffix(".summary.json"), summary)
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("public", type=Path)
    parser.add_argument("--output", type=Path, default=Path("predictions.tsv"))
    parser.add_argument("--method", choices=("auto", *b.METHODS), default="auto")
    args = parser.parse_args()
    print(solve(args.public, args.output, args.method))
