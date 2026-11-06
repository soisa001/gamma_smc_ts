#!/usr/bin/env python3
"""Draw a conceptual gamma posterior and its mass/median bin summaries.

Analytic example only: Gamma(shape=1, scale=20 kyr), not fitted simulation output.
No RNG draws. Outputs PNG and vector PDF with editable text.
"""
import os
os.environ['MPLBACKEND'] = 'Agg'
import argparse
import hashlib
import json
from pathlib import Path
import re
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import PercentFormatter
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=ROOT / 'docs/results/posterior_mass_intuition_20261106')
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    scale = 20.
    median = scale * np.log(2)
    times = np.linspace(0, 80, 1601)
    density = np.exp(-times / scale) / scale
    edges = np.arange(51)
    mass = np.exp(-edges[:-1] / scale) - np.exp(-edges[1:] / scale)
    median_weights = ((edges[:-1] <= median) & (median < edges[1:])).astype(float)
    assert np.isclose(mass.sum() + np.exp(-50 / scale), 1)
    assert median_weights.sum() == 1 and np.isclose(1 - np.exp(-median / scale), .5)
    matplotlib.rcParams.update({'font.size': 15, 'axes.labelsize': 16, 'axes.titlesize': 16,
                               'xtick.labelsize': 13, 'ytick.labelsize': 13,
                               'pdf.fonttype': 42, 'ps.fonttype': 42})
    fig = plt.figure(figsize=(11, 8.5))
    grid = fig.add_gridspec(2, 2, height_ratios=[1.25, 1])
    top = fig.add_subplot(grid[0, :])
    left, right = fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])
    fig.suptitle('One broad posterior: mass and median summarize it differently', y=.98, fontsize=18)
    fig.text(.5, .935, 'Illustrative Gamma(shape = 1, scale = 20 kyr) = exponential posterior; not fitted data', ha='center', fontsize=12.5)
    top.plot(times, density, color='#0072B2', lw=3)
    shade_x = np.linspace(0, median, 300)
    top.fill_between(shade_x, np.exp(-shade_x / scale) / scale, color='#0072B2', alpha=.15)
    top.axvline(median, color='#D55E00', lw=2.3)
    top.axvline(scale, color='#555555', ls='--', lw=2)
    top.text(6.8, .011, '50% of\nprobability', ha='center', fontsize=13)
    top.set(xlim=(0, 80), ylim=(0, .055), xlabel='Pairwise TMRCA (kya)', ylabel='Posterior density (% per kyr)')
    top.yaxis.set_major_formatter(PercentFormatter(1))
    fig.legend(handles=[Line2D([], [], color='#0072B2', lw=3, label='Posterior curve; mode = 0'),
                        Line2D([], [], color='#D55E00', lw=2.3, label=f'Median = {median:.2f} kya'),
                        Line2D([], [], color='#555555', ls='--', lw=2, label='Mean = 20 kya')],
               loc='upper center', bbox_to_anchor=(.52, .895), ncol=3, frameon=False, fontsize=12)
    left.bar(edges[:-1] + .5, mass, width=.9, color='#0072B2')
    left.set_title('Posterior-mass contribution', loc='left')
    left.set_ylim(0, .055)
    left.text(.98, .95, f'1–2 kya: {mass[1]:.2%}\n13–14 kya: {mass[13]:.2%}',
              transform=left.transAxes, ha='right', va='top', fontsize=13)
    right.bar(edges[:-1] + .5, median_weights, width=.9, color='#D55E00')
    right.set_title('Posterior-median contribution', loc='left')
    right.set_ylim(0, 1.08)
    right.text(.98, .95, '1–2 kya: 0%\n13–14 kya: 100%', transform=right.transAxes,
               ha='right', va='top', fontsize=13)
    for ax in (left, right):
        ax.set(xlim=(0, 50), xlabel='Pairwise TMRCA bin (kya)', ylabel='Pair weight per 1-kyr bin')
        ax.yaxis.set_major_formatter(PercentFormatter(1))
        ax.set_xticks([0, 10, 20, 30, 40, 50])
    for ax in (top, left, right):
        ax.grid(axis='y', alpha=.18)
        ax.set_axisbelow(True)
    fig.text(.5, .045, 'Bottom panels have different y-scales. Mass retains uncertainty; the median assigns one entire pair to one bin.',
             ha='center', fontsize=10.5)
    fig.text(.5, .018, f'Posterior probability older than 50 kya: {np.exp(-50 / scale):.2%}; retained, not renormalized away.',
             ha='center', fontsize=11)
    fig.subplots_adjust(left=.11, right=.98, top=.82, bottom=.14, hspace=.62, wspace=.32)
    fig.savefig(args.out / 'posterior_curve.png', dpi=180)
    fig.savefig(args.out / 'posterior_curve.pdf')
    plt.close(fig)
    pdf = (args.out / 'posterior_curve.pdf').read_bytes()
    assert b'/FontFile2' in pdf and not re.search(rb'/Subtype\s*/Image', pdf)
    summary = dict(example='Conceptual distribution, not fitted to simulation data', shape=1, scale_kyr=scale,
                   mode_kya=0, median_kya=median, mean_kya=scale, mass_1_2_kya=float(mass[1]),
                   mass_13_14_kya=float(mass[13]), mass_above_50_kya=float(np.exp(-50 / scale)),
                   validation='Analytic CDF at median equals 0.5; bin masses plus older tail sum to one; median weights sum to one. Embedded TrueType fonts; no raster PDF objects.',
                   randomness='None; analytic deterministic functions.',
                   script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   versions=dict(numpy=np.__version__, matplotlib=matplotlib.__version__))
    (args.out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    (args.out / '.gitattributes').write_text('* -text\n', encoding='utf-8')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
