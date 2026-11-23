"""Plot the selected and held-out neutral AF trajectories for the recent EAS sweep."""
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


def read(path):
    return json.loads(path.read_text())


def main():
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

    colors = ["#495C70", "#A7461F"]
    names = ["Neutral", "Selected"]
    plt.rcParams.update({"font.size": 14, "pdf.fonttype": 42, "ps.fonttype": 42})
    fig = plt.figure(figsize=(11, 8.5))
    grid = fig.add_gridspec(2, 2, height_ratios=[1.15, 1], left=.095, right=.975,
                            top=.79, bottom=.095, wspace=.16, hspace=.54)
    full = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])]
    recent = fig.add_subplot(grid[1, :])
    fig.suptitle("Allele-frequency trajectories after 1% introgression", y=.975, fontsize=23)
    fig.text(.5, .902, "EAS simulations  ·  100 replicates per cohort  ·  Focal allele segregating today",
             ha="center", fontsize=15)
    recent.set_title("Final 10,000 years", fontsize=17, loc="left")
    for ax, af, name, color in zip(full, cohorts, names, colors):
        assert af.shape == (100, len(time)) and np.isfinite(af).all()
        assert np.all((af >= 0) & (af <= 1)) and np.all((af[:, -1] > 0) & (af[:, -1] < 1))
        q1, median, q3 = np.percentile(af, [25, 50, 75], axis=0)
        ax.plot(time, af.T, color=color, lw=.55, alpha=.25)
        ax.plot(time, median, color=color, lw=2.6)
        ax.set(xlim=(50, 0), ylim=(0, 1), xticks=[50, 40, 30, 20, 10, 0])
        ax.set_title("Neutral" if name == "Neutral" else "Selected from 10 kya to present (s = 0.01)", fontsize=15)
        if name == "Selected":
            ax.axvline(10, color="#222222", lw=1.1, ls="--")
        keep = time <= 10
        recent.fill_between(time[keep], q1[keep], q3[keep], color=color, alpha=.18, linewidth=0)
        recent.plot(time[keep], median[keep], color=color, lw=2.6)
        recent.annotate(f"{median[-1]:.1%}", (0, median[-1]), xytext=(6, 0),
                        textcoords="offset points", va="center", color=color, fontsize=13)
    full[0].set_ylabel("Allele frequency")
    full[1].tick_params(labelleft=False)
    recent.set(xlim=(10, -.7), ylim=(0, .65), xticks=[10, 8, 6, 4, 2, 0],
               yticks=[0, .2, .4, .6], xlabel="Time before present (kya)", ylabel="Allele frequency")
    for ax in [*full, recent]:
        ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
        ax.set_axisbelow(True)
        ax.grid(axis="y", color="#E4E7E9", lw=.7)
        ax.spines[["top", "right"]].set_visible(False)
    fig.legend([Line2D([], [], color="#637282", lw=.8), Line2D([], [], color="#333333", lw=2.6),
                Line2D([], [], color="#222222", lw=1.1, ls="--")],
               ["Individual replicate", "Pointwise median", "Selection onset"],
               ncol=3, frameon=False, loc="center", bbox_to_anchor=(.50, .85), fontsize=13)
    recent.legend([Line2D([], [], color=c, lw=2.6) for c in colors] + [Patch(color="#7D8993", alpha=.22)],
                  ["Neutral median", "Selected median", "Middle 50%"], ncol=3, frameon=False,
                  loc="lower right", bbox_to_anchor=(1, 1.0), fontsize=12)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"af_10k_s0p01_introgression.{ext}", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
