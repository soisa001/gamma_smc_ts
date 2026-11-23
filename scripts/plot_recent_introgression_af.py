"""Plot the selected and held-out neutral AF trajectories for the recent EAS sweep."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path

os.environ["MPLBACKEND"] = "Agg"
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import PercentFormatter
import numpy as np

RESULTS = Path(__file__).resolve().parents[1] / "docs/results"
OUT = RESULTS / "eas_strength_series_20261118"
CASE = "introgression_10k_s0p01"
EXTRA = ("de_novo_10k_s0p01", "de_novo_10k_s0p03", "introgression_10k_s0p03")
CACHE = OUT / "recent_af_trajectory_inputs.npz"


def read(path):
    return json.loads(path.read_text())


def original():
    selected = RESULTS / "denovo_comparison_20261111" / f"{CASE}.npz"
    neutral = RESULTS / "segregating_neutral_fpr_100_20261110/parameters.npz"
    for path, expected in (
        (selected, read(selected.parent / "manifest.json")["cohorts"][CASE]["cache_sha256"]),
        (neutral, read(neutral.parent / "output_hashes.json")[neutral.name]),
    ):
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected
    with np.load(selected) as s, np.load(neutral) as n, np.load(OUT / "mass.npz") as bank:
        selected_ids = [r["source"] for r in read(OUT / "inputs.json")["records"] if r["case"] == CASE]
        np.testing.assert_array_equal(s["ids"], selected_ids)
        np.testing.assert_array_equal(n["test_ids"], bank["introgression_test_ids"])
        assert not bank[CASE + "_fixed"].any()
        time = s["af_times_kya"]
        np.testing.assert_array_equal(time, np.arange(50000, -1, -25) / 1000)
        cohorts = [n["test_af"], s["af"]]
    plot(time, cohorts, "introgression", .01)


def load_history(item):
    source, expected_score, strength, oldest = item
    directory = Path(source)
    receipt = read(directory / "complete.json")
    trajectory = directory / "trajectory.tsv"
    assert receipt["scores_sha256"] == expected_score and receipt["selection"] == strength
    assert hashlib.sha256(trajectory.read_bytes()).hexdigest() == receipt["trajectory_sha256"]
    years, af, size, _, copies = np.loadtxt(trajectory, skiprows=1).T
    np.testing.assert_array_equal(years, np.arange(oldest * 1000, -1, -25))
    np.testing.assert_allclose(af, copies / (2 * size), rtol=1e-12)
    assert af[-1] == receipt["final_frequency"] and 0 < af[-1] < 1
    return af


def prepare():
    previous = RESULTS / "denovo_comparison_20261111"
    records = read(OUT / "inputs.json")["records"]
    selected = {k: [r for r in records if r["case"] == k] for k in EXTRA}
    entries = {}
    with np.load(OUT / "mass.npz") as bank, ThreadPoolExecutor(max_workers=4) as pool:
        for k, rows in selected.items():
            assert len(rows) == 100 and not bank[k + "_fixed"].any()
            np.testing.assert_array_equal([r["slot"] for r in rows], np.arange(100))
            entries[k + "_ids"] = np.array([r["source"] for r in rows])
            if k.endswith("s0p01"):
                path = previous / (k + ".npz")
                assert hashlib.sha256(path.read_bytes()).hexdigest() == read(previous / "manifest.json")["cohorts"][k]["cache_sha256"]
                with np.load(path) as data:
                    np.testing.assert_array_equal(data["ids"], entries[k + "_ids"])
                    np.testing.assert_array_equal(data["af_times_kya"], np.arange(10000, -1, -25) / 1000)
                    entries[k] = data["af"]
            else:
                oldest = 10 if k.startswith("de_novo") else 50
                entries[k] = np.stack(list(pool.map(load_history, [(r["source"], r["scores_sha256"], .03, oldest) for r in rows])))
        ids = bank["de_novo_10k_test_ids"]
        previous_sources = read(RESULTS / "cumulative_grid_20261117/inputs.json")["sources"]
        hashes = {r["id"]: r["scores_sha256"] for r in previous_sources if r["group"] == "heldout_de_novo_10k"}
        assert len(set(ids)) == len(hashes) == 100
        entries["de_novo_neutral_ids"] = ids
        entries["de_novo_neutral"] = np.stack(list(pool.map(load_history, [(p, hashes[p], 0., 10) for p in ids])))
        path = RESULTS / "segregating_neutral_fpr_100_20261110/parameters.npz"
        assert hashlib.sha256(path.read_bytes()).hexdigest() == read(path.parent / "output_hashes.json")[path.name]
        with np.load(path) as controls:
            np.testing.assert_array_equal(controls["test_ids"], bank["introgression_test_ids"])
            entries["introgression_neutral_ids"] = controls["test_ids"]
            entries["introgression_neutral"] = controls["test_af"]
    np.savez_compressed(CACHE, **entries)


def extended():
    with np.load(CACHE) as data:
        for k in EXTRA:
            origin, strength = k.split("_10k_s")
            oldest = 10 if origin == "de_novo" else 50
            time = np.arange(oldest * 1000, -1, -25) / 1000
            plot(time, [data[origin + "_neutral"], data[k]], origin, float(strength.replace("p", ".")))


def plot(time, cohorts, origin, strength):
    introgressed = origin == "introgression"
    colors = ["#495C70", "#A7461F"]
    names = ["Neutral", "Selected"]
    plt.rcParams.update({"font.size": 14, "pdf.fonttype": 42, "ps.fonttype": 42})
    fig = plt.figure(figsize=(11, 8.5))
    grid = fig.add_gridspec(2, 2, height_ratios=[1.15, 1], left=.095, right=.975,
                            top=.79, bottom=.095, wspace=.16, hspace=.54)
    full = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])]
    recent = fig.add_subplot(grid[1, :])
    title = "Allele-frequency trajectories after 1% introgression" if introgressed else "Allele-frequency trajectories from a de novo mutation"
    fig.suptitle(title, y=.975, fontsize=23)
    fig.text(.5, .902, "EAS simulations  ·  100 replicates per cohort  ·  Focal allele segregating today",
             ha="center", fontsize=15)
    recent.set_title("Final 10,000 years" if introgressed else "Median trajectories", fontsize=17, loc="left")
    top = 1.0 if introgressed or strength == .03 else min(1., .01 * np.ceil(max(a.max() for a in cohorts) * 100))
    peak = max(np.percentile(a[:, time <= 10], 75, axis=0).max() for a in cohorts)
    recent_top = .65 if introgressed and strength == .01 else min(1., peak * 1.18)
    for ax, af, name, color in zip(full, cohorts, names, colors):
        assert af.shape == (100, len(time)) and np.isfinite(af).all()
        assert np.all((af >= 0) & (af <= 1)) and np.all((af[:, -1] > 0) & (af[:, -1] < 1))
        q1, median, q3 = np.percentile(af, [25, 50, 75], axis=0)
        ax.plot(time, af.T, color=color, lw=.55, alpha=.25)
        ax.plot(time, median, color=color, lw=2.6)
        ax.set(xlim=(time[0], 0), ylim=(0, top), xticks=np.linspace(time[0], 0, 6))
        ax.set_title("Neutral" if name == "Neutral" else f"Selected from 10 kya to present (s = {strength:g})", fontsize=15)
        if name == "Selected" and introgressed:
            ax.axvline(10, color="#222222", lw=1.1, ls="--")
        keep = time <= 10
        recent.fill_between(time[keep], q1[keep], q3[keep], color=color, alpha=.18, linewidth=0)
        recent.plot(time[keep], median[keep], color=color, lw=2.6)
        label = f"{median[-1]:.2%}" if median[-1] < .01 else f"{median[-1]:.1%}"
        label_offset = 5 if median[-1] < .02 * recent_top else 0
        recent.annotate(label, (0, median[-1]), xytext=(6, label_offset),
                        textcoords="offset points", va="center", color=color, fontsize=13)
    full[0].set_ylabel("Allele frequency")
    full[1].tick_params(labelleft=False)
    recent.set(xlim=(10, -.8), ylim=(0, recent_top), xticks=[10, 8, 6, 4, 2, 0],
               xlabel="Time before present (kya)", ylabel="Allele frequency")
    for ax in [*full, recent]:
        ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=2 if ax.get_ylim()[1] < .05 else 0))
        ax.set_axisbelow(True)
        ax.grid(axis="y", color="#E4E7E9", lw=.7)
        ax.spines[["top", "right"]].set_visible(False)
    handles = [Line2D([], [], color="#637282", lw=.8), Line2D([], [], color="#333333", lw=2.6)]
    labels = ["Individual replicate", "Pointwise median"]
    if introgressed:
        handles.append(Line2D([], [], color="#222222", lw=1.1, ls="--"))
        labels.append("Selection onset")
    fig.legend(handles, labels,
               ncol=3, frameon=False, loc="center", bbox_to_anchor=(.50, .85), fontsize=13)
    recent.legend([Line2D([], [], color=c, lw=2.6) for c in colors] + [Patch(color="#7D8993", alpha=.22)],
                  ["Neutral median", "Selected median", "Middle 50%"], ncol=3, frameon=False,
                  loc="lower right", bbox_to_anchor=(1, 1.0), fontsize=12)
    for ext in ("png", "pdf"):
        value = str(strength).replace(".", "p")
        fig.savefig(OUT / f"af_10k_s{value}_{origin}.{ext}", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", nargs="?", default="original", choices=("original", "prepare", "extended"))
    phase = parser.parse_args().phase
    {"original": original, "prepare": prepare, "extended": extended}[phase]()
