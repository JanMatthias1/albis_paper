#!/usr/bin/env python3
"""
Figure 4A (step 1, layout option B) -- "error-vector" infographic of
pairwise translation error, cell resolution, all 4 algorithms.

For every algorithm x axis x slice-pair, evaluation_moving_metrics.py
already anchors the reference slice onto its own ground truth (a similarity
transform of the reference true->aligned frame) and reports the *moving*
slice's residual translation in that anchored frame:
    moving_relative_translation_x/y   (um, in the truth-anchored frame)
    moving_relative_translation_norm  (um, == |vector|)
This is a direct measurement of "if the algorithm's own frame is correctly
oriented via its better-constrained reference slice, how far off is the
other slice's placement" -- i.e. translation error. Reads the same
already-computed CSVs as figure4a_metric_comparison.py; nothing is rerun.

One panel per (axis, algorithm): 45 arrows from the origin to each pair's
residual (moving_relative_translation_x, moving_relative_translation_y),
plus a dashed circle at the median |translation|.

Produces:
    figure4a_error_vectors.png
"""

import argparse
import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

METRIC_ROOT = (
    "/dcs04/hicks/data/multi-sample-alignment-benchmark/evaluation_result/"
    "simulated_pairwise_relative_metrics"
)
OUT = "/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_4/alignment_translation_error"

ALGO_DIRNAME = {"STAIR": "stair", "SLAT": "slat", "SPACEL": "spacel", "STAligner": "staligner"}
ALGO_COLOR = {"STAIR": "#4878d0", "SLAT": "#ee854a", "SPACEL": "#6acc64", "STAligner": "#d65f5f"}
ALGOS = ["STAIR", "SLAT", "SPACEL", "STAligner"]
AXES = ["X", "Y", "Z"]
DATASET = "cell_obs_sectioned"


def load(dataset, axis, algo):
    p = os.path.join(
        METRIC_ROOT, ALGO_DIRNAME[algo],
        f"{dataset}_{axis}_relative_alignment_metrics_all_pairs.csv",
    )
    df = pd.read_csv(p)
    return df[df["status"] == "success"].copy()


def plot(dataset, axes, outpath):
    fig, grid = plt.subplots(len(axes), len(ALGOS),
                              figsize=(3.4 * len(ALGOS), 3.4 * len(axes)),
                              squeeze=False)

    for row, axis in enumerate(axes):
        panels = {algo: load(dataset, axis, algo) for algo in ALGOS}
        # shared, symmetric axis limits across algorithms within this axis row
        lim = max(
            np.abs(np.concatenate([df["moving_relative_translation_x"].values,
                                    df["moving_relative_translation_y"].values])).max()
            for df in panels.values()
        ) * 1.08

        for col, algo in enumerate(ALGOS):
            ax = grid[row][col]
            df = panels[algo]
            tx, ty = df["moving_relative_translation_x"].values, df["moving_relative_translation_y"].values
            norm = df["moving_relative_translation_norm"].values
            median_norm = np.median(norm)

            ax.quiver(np.zeros(len(tx)), np.zeros(len(ty)), tx, ty,
                      angles="xy", scale_units="xy", scale=1,
                      color=ALGO_COLOR[algo], alpha=0.45, width=0.006,
                      headwidth=3, headlength=4)
            theta = np.linspace(0, 2 * np.pi, 100)
            ax.plot(median_norm * np.cos(theta), median_norm * np.sin(theta),
                    color="black", ls="--", lw=1.1, alpha=0.8)
            ax.scatter([0], [0], color="black", s=14, zorder=5)

            ax.set_xlim(-lim, lim)
            ax.set_ylim(-lim, lim)
            ax.set_aspect("equal")
            ax.axhline(0, color="0.85", lw=0.6, zorder=0)
            ax.axvline(0, color="0.85", lw=0.6, zorder=0)
            ax.set_title(f"{algo}\nmedian |t| = {median_norm:,.0f} µm", fontsize=10.5)
            if col == 0:
                ax.set_ylabel(f"axis {axis}\nΔy (µm)", fontsize=10)
            else:
                ax.set_yticklabels([])
            if row == len(axes) - 1:
                ax.set_xlabel("Δx (µm)", fontsize=10)
            else:
                ax.set_xticklabels([])

    fig.suptitle("Figure 4A -- residual translation error per slice pair, cell resolution\n"
                  "(45 pairs/algorithm/axis; arrow = moving slice's position error in the "
                  "truth-anchored frame; dashed circle = median magnitude)", fontsize=11.5)
    fig.patch.set_facecolor("white")
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    fig.savefig(outpath, dpi=200, facecolor="white")
    plt.close(fig)
    print("wrote", outpath)


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default=DATASET)
    ap.add_argument("--axes", nargs="+", default=AXES, choices=AXES)
    ap.add_argument("--outdir", default=OUT)
    return ap.parse_args()


if __name__ == "__main__":
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    plot(args.dataset, args.axes, os.path.join(args.outdir, "figure4a_error_vectors.png"))
