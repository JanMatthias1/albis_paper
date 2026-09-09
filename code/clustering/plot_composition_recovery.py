#!/usr/bin/env python
"""
Figure 3B (reframed): the proposed paper panel for "how well does structure
survive bin/spot/cell aggregation and re-clustering", scoring against the SOFT
ground truth instead of hard-label ARI.

Reads the summary written by composition_recovery.py:
    data/figure_3/composition_recovery/composition_recovery_summary.csv

Writes to the same directory:
    composition_recovery_fig3b.png          2 rows (cell_type / domain) x 3 cols
    composition_recovery_fig3b_celltype.png  the cell_type row alone (talk slide)
    composition_recovery_fig3b_domain.png    the domain row alone

Panel order per row, per PI feedback 2026-09-06 (lead with the fair metrics,
demote hard-label ARI):
    1. Local kNN label concordance vs the sum(p_k^2) chance line
       -- does the embedding keep same-dominant-category units together?
    2. Composition variance explained: pipeline vs the KMeans-on-truth oracle
       ceiling -- of the composition structure that is recoverable at all, how
       much does the pipeline get?
    3. Hard-label ARI (+ V-measure), greyed -- the previous metric, shown only
       for reference; it treats a mixture-of-types aggregate like a labelled
       cell, which is the thing this reframe is fixing.

The cell_type row x-labels also carry the effective number of cell-types per
observation (exp Shannon entropy of the true fraction vector) -- 1.0 for cell,
up to ~6 for spot -- which is the aggregation confound the whole panel is about.

Expected environment:
    source code/count_distribution/_env.sh
Example:
    python sim_paper/code/clustering/plot_composition_recovery.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent
SIM_PAPER_DIR = SCRIPT_DIR.parents[1]
CR_DIR = SIM_PAPER_DIR / "data" / "figure_3" / "composition_recovery"

MODS = ["cell", "bin", "bin16um", "spot"]
MOD_LABEL = {"cell": "cell", "bin": "bin 8µm", "bin16um": "bin 16µm", "spot": "spot"}

# two-value encoding: the pipeline result vs a reference (chance / ceiling /
# deprecated). Not a categorical ramp -- one accent, one recessive grey.
C_METHOD = "#1f6fb2"      # pipeline / observed
C_CEILING = "#c7ccd1"     # oracle ceiling / reference fill
C_REF_INK = "#6b7280"     # chance line, reference text
C_DEPRECATED = "#b8bec6"  # greyed hard-ARI bars
INK = "#1f2933"
MUTED = "#6b7280"
GRID = "#d9dde1"

plt.rcParams.update({
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "xtick.labelsize": 9, "ytick.labelsize": 8, "legend.fontsize": 8,
    "axes.edgecolor": "#b0b6bc", "axes.linewidth": 0.8,
    "text.color": INK, "axes.labelcolor": INK,
    "xtick.color": "#5b6169", "ytick.color": "#5b6169",
    "figure.facecolor": "white", "axes.facecolor": "white",
})


def _bars(ax, values, color, width=0.62, x=None):
    x = np.arange(len(values)) if x is None else x
    b = ax.bar(x, values, width, color=color, edgecolor="white", linewidth=1.4, zorder=3)
    return b


def _label_bars(ax, bars, values, fmt="{:.2f}", dy=0.012, color=INK):
    for rect, v in zip(bars, values):
        if v is None or (isinstance(v, float) and np.isnan(v)):
            continue
        ax.text(rect.get_x() + rect.get_width() / 2, rect.get_height() + dy,
                fmt.format(v), ha="center", va="bottom", fontsize=8, color=color)


def _tidy(ax, ymax=1.0, ylabel=""):
    ax.set_ylim(0, ymax)
    ax.set_ylabel(ylabel)
    ax.grid(axis="y", color=GRID, linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def panel_knn(ax, d, xlabels):
    x = np.arange(len(MODS))
    obs = d["knn_concordance"].to_numpy()
    chance = float(np.nanmean(d["knn_chance"]))
    bars = _bars(ax, obs, C_METHOD)
    _label_bars(ax, bars, obs)
    ax.axhline(chance, color=C_REF_INK, linestyle=(0, (4, 3)), linewidth=1.4, zorder=4,
              label=f"chance (Σpₖ²) ≈ {chance:.2f}")
    ax.set_xticks(x); ax.set_xticklabels(xlabels)
    _tidy(ax, 1.0, "fraction of 15 NN sharing dominant label")
    ax.set_title("Local kNN label concordance", loc="left")
    ax.legend(frameon=False, loc="upper right", handlelength=1.6)


def panel_r2(ax, d, xlabels):
    x = np.arange(len(MODS))
    w = 0.38
    ceil = d["comp_r2_oracle"].to_numpy()
    meth = d["comp_r2"].to_numpy()
    b1 = ax.bar(x - w / 2, ceil, w, color=C_CEILING, edgecolor="white", linewidth=1.2,
                zorder=3, label="oracle ceiling (k-means on true composition)")
    b2 = ax.bar(x + w / 2, meth, w, color=C_METHOD, edgecolor="white", linewidth=1.2,
                zorder=3, label="pipeline (Leiden on Harmony PCA)")
    _label_bars(ax, b1, ceil)
    _label_bars(ax, b2, meth)
    for xi, m, c in zip(x, meth, ceil):
        if c > 0.01:
            ax.text(xi + w / 2, m + 0.055, f"{m / c:.0%}\nof ceiling",
                    ha="center", va="bottom", fontsize=7, color=MUTED, linespacing=0.95)
    ax.set_xticks(x); ax.set_xticklabels(xlabels)
    _tidy(ax, 1.32, "composition variance explained  (R²)")
    ax.set_title("Composition recovery vs oracle ceiling", loc="left")
    ax.legend(frameon=False, loc="upper right", handlelength=1.3)


def panel_ari(ax, d, xlabels):
    x = np.arange(len(MODS))
    w = 0.38
    ari = d["ari"].to_numpy()
    vme = d["v_measure"].to_numpy()
    ax.set_facecolor("#faf6f6")
    b1 = ax.bar(x - w / 2, ari, w, color=C_DEPRECATED, edgecolor="white", linewidth=1.2,
                zorder=3, label="ARI")
    b2 = ax.bar(x + w / 2, vme, w, color="#d7d2c4", edgecolor="white", linewidth=1.2,
                zorder=3, label="V-measure")
    _label_bars(ax, b1, ari, color=MUTED)
    _label_bars(ax, b2, vme, color=MUTED)
    ax.set_xticks(x); ax.set_xticklabels(xlabels)
    _tidy(ax, 1.0, "score")
    ax.set_title("Hard-label ARI  —  deprecated for aggregates", loc="left", color=MUTED)
    ax.legend(frameon=False, loc="upper right", handlelength=1.3)


def _xlabels(df_level, with_eff):
    out = []
    for m in MODS:
        row = df_level[df_level.modality == m]
        lab = MOD_LABEL[m]
        if with_eff and len(row):
            lab += f"\n≈{float(row['eff_cats_per_obs_median'].iloc[0]):.1f} types/obs"
        out.append(lab)
    return out


def make_row(axes, df, level, with_eff):
    d = df[df.level == level].set_index("modality").reindex(MODS).reset_index()
    xl = _xlabels(d, with_eff)
    panel_knn(axes[0], d, xl)
    panel_r2(axes[1], d, xl)
    panel_ari(axes[2], d, xl)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", type=Path, default=CR_DIR / "composition_recovery_summary.csv")
    ap.add_argument("--out-dir", type=Path, default=CR_DIR)
    args = ap.parse_args()
    df = pd.read_csv(args.csv)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    # combined 2 x 3
    fig, axes = plt.subplots(2, 3, figsize=(14.5, 9.2))
    make_row(axes[0], df, "cell_type", with_eff=True)
    make_row(axes[1], df, "domain", with_eff=False)
    fig.text(0.006, 0.735, "cell-type structure", rotation=90, va="center",
             ha="center", fontsize=11, weight="bold", color=INK)
    fig.text(0.006, 0.265, "spatial-domain structure", rotation=90, va="center",
             ha="center", fontsize=11, weight="bold", color=INK)
    fig.suptitle("Figure 3B  ·  recovery of aggregated structure after re-clustering",
                 x=0.5, y=0.995, fontsize=12.5, weight="bold")
    fig.text(0.5, 0.945,
             "cell keeps hard-label ARI (a cell has one type); bin / spot are mixtures "
             "— read the kNN and composition panels, not ARI",
             ha="center", fontsize=8.5, color=MUTED)
    fig.tight_layout(rect=[0.014, 0.0, 1, 0.935])
    p = args.out_dir / "composition_recovery_fig3b.png"
    fig.savefig(p, dpi=200, bbox_inches="tight")
    print(f"[plot] {p}")

    # single-row cut-outs
    for level, tag, eff in [("cell_type", "celltype", True), ("domain", "domain", False)]:
        f2, ax2 = plt.subplots(1, 3, figsize=(14.5, 4.6))
        make_row(ax2, df, level, with_eff=eff)
        f2.suptitle(f"Figure 3B · {level.replace('_', '-')} structure recovery",
                    x=0.5, y=1.02, fontsize=12, weight="bold")
        f2.tight_layout()
        pp = args.out_dir / f"composition_recovery_fig3b_{tag}.png"
        f2.savefig(pp, dpi=200, bbox_inches="tight")
        print(f"[plot] {pp}")


if __name__ == "__main__":
    main()
