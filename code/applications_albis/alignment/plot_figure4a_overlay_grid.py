#!/usr/bin/env python3
"""
Figure 4A (step 1, layout option A) -- overlay grid of all 10 slices, before
vs. after alignment vs. ground truth, extending the existing Fig4C style
(plot_figure4c.py) to multiple algorithms.

Reads the joint (all-10-slices-at-once, not pairwise) alignment results from
the multi-sample-alignment-benchmark repo:
    results/sim_data/3D/<ALGO>/cell_obs_sectioned_<axis>/adata_results/
        Sim_3D_<ALGO>_dataset_cell_obs_sectioned_<axis>.h5ad
Nothing is rerun -- these already exist for STAIR, SPACEL, STAligner.

SLAT is NOT included: its joint 10-slice alignment (3D_slat.py, run_SLAT_multi)
OOM'd on GPU at cell resolution for every axis (confirmed from
code/Jan/sim_data/3D_alignment/logs_slat_3D/*cell_obs_sectioned*.out -- "CUDA
out of memory, tried to allocate 20.60 GiB"), so no such file exists to read.
SLAT's cell-resolution translation error is shown instead in the companion
figure4a_error_vectors.png, built from its (successful) pairwise runs.

Produces one file per axis:
    figure4a_overlay_2d_<axis>.png   rows = algorithm, cols = [Original |
                                      Aligned | Ground truth], all 10 slices
                                      overlaid and colored by slice_id
"""

import argparse
import os

import numpy as np
import pandas as pd
import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

JOINT_3D_ROOT = "/dcs04/hicks/data/multi-sample-alignment-benchmark/results/sim_data/3D"
OUT = "/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_4/alignment_translation_error"

DATASET = "cell_obs_sectioned"
AXES = ["X", "Y", "Z"]
COL_LABELS = ["Original", "Aligned", "Ground truth"]

# algorithms with a joint (all-10-slice) result at cell resolution
ALGOS = ["STAIR", "SPACEL", "STAligner"]
ALGO_KEYS = {
    "STAIR": dict(original="spatial", aligned="transform_fine", true="spatial_3d_true"),
    "SPACEL": dict(original="spatial", aligned="spatial_aligned", true="spatial_3d_true"),
    "STAligner": dict(original="spatial_orig", aligned="spatial", true="spatial_3d_true"),
}
SLAT_NOTE = ("SLAT: joint 10-slice alignment OOM'd on GPU at cell resolution\n"
             "(run_SLAT_multi, 20.6 GiB alloc failed on every axis).\n"
             "See figure4a_error_vectors.png for SLAT's pairwise-based result.")


def procrustes_fit(X, Y):
    """Best rigid map X -> Y (rotation + translation, no scale/reflection)."""
    X, Y = np.asarray(X, float), np.asarray(Y, float)
    mx, my = X.mean(0), Y.mean(0)
    U, _, Vt = np.linalg.svd((X - mx).T @ (Y - my))
    R = U @ Vt
    if np.linalg.det(R) < 0:
        U[:, -1] *= -1
        R = U @ Vt
    return (X - mx) @ R + my


def load(algo, axis):
    p = os.path.join(
        JOINT_3D_ROOT, algo, f"{DATASET}_{axis}", "adata_results",
        f"Sim_3D_{algo}_dataset_{DATASET}_{axis}.h5ad",
    )
    return ad.read_h5ad(p)


def subsample(n, max_points, seed=0):
    if n <= max_points:
        return np.arange(n)
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(n, max_points, replace=False))


def slice_colors(adata):
    cats = sorted(adata.obs["slice_id"].unique(), key=float)
    cmap = plt.get_cmap("tab10")
    lut = {c: cmap(i % 10) for i, c in enumerate(cats)}
    colors = adata.obs["slice_id"].astype(str).map({str(k): v for k, v in lut.items()}).values
    return np.vstack(colors), cats, lut


def panels_for(adata, keys):
    orig_xy = np.asarray(adata.obsm[keys["original"]], float)[:, :2]
    align_xy = np.asarray(adata.obsm[keys["aligned"]], float)[:, :2]
    true_xy = np.asarray(adata.obsm[keys["true"]], float)[:, :2]
    align_xy = procrustes_fit(align_xy, true_xy)  # display-only anchor
    return [orig_xy, align_xy, true_xy]


def plot_axis(axis, outdir):
    data = {}
    for algo in ALGOS:
        try:
            data[algo] = load(algo, axis)
        except FileNotFoundError as e:
            print(f"skip {algo}/{axis}: {e}")

    n_rows = len(ALGOS) + 1  # + SLAT placeholder row
    fig, grid = plt.subplots(n_rows, 3, figsize=(11.5, 3.6 * n_rows), squeeze=False)

    cats_all, lut_all = None, None
    for row, algo in enumerate(ALGOS):
        adata = data.get(algo)
        if adata is None:
            for c in range(3):
                grid[row][c].axis("off")
            continue
        col_vec, cats, lut = slice_colors(adata)
        cats_all, lut_all = cats, lut
        idx = subsample(adata.n_obs, 60000)
        cols3 = panels_for(adata, ALGO_KEYS[algo])
        s = 0.8 if adata.n_obs > 100000 else 3.0
        for c in range(3):
            ax = grid[row][c]
            xy = cols3[c][idx]
            ax.scatter(xy[:, 0], xy[:, 1], s=s, c=col_vec[idx], linewidths=0, alpha=0.8)
            ax.set_aspect("equal")
            ax.axis("off")
            if row == 0:
                ax.set_title(COL_LABELS[c], fontsize=13)
        grid[row][0].text(-0.06, 0.5, algo, transform=grid[row][0].transAxes,
                           fontsize=13, fontweight="bold", ha="right", va="center",
                           rotation=90)

    # SLAT placeholder row
    slat_row = n_rows - 1
    for c in range(3):
        ax = grid[slat_row][c]
        ax.axis("off")
        if c == 1:
            ax.text(0.5, 0.5, SLAT_NOTE, transform=ax.transAxes, fontsize=9.5,
                     ha="center", va="center", style="italic", color="0.35",
                     bbox=dict(boxstyle="round", fc="0.96", ec="0.8"))
    grid[slat_row][0].text(-0.06, 0.5, "SLAT", transform=grid[slat_row][0].transAxes,
                            fontsize=13, fontweight="bold", ha="right", va="center",
                            rotation=90, color="0.5")

    if cats_all is not None:
        handles = [Line2D([0], [0], marker="o", ls="", mfc=lut_all[c], mec="none", label=str(c))
                   for c in cats_all]
        fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.0),
                   borderaxespad=0.2, ncol=min(len(cats_all), 10), frameon=False,
                   title="slice_id")

    fig.suptitle(f"Figure 4A -- all-10-slice alignment overlay, cell resolution, axis {axis}",
                 fontsize=13)
    fig.patch.set_facecolor("white")
    fig.tight_layout(rect=[0.02, 0.035, 1, 0.95])
    outpath = os.path.join(outdir, f"figure4a_overlay_2d_{axis}.png")
    fig.savefig(outpath, dpi=200, facecolor="white")
    plt.close(fig)
    print("wrote", outpath)


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--axes", nargs="+", default=AXES, choices=AXES)
    ap.add_argument("--outdir", default=OUT)
    return ap.parse_args()


if __name__ == "__main__":
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    for axis in args.axes:
        plot_axis(axis, args.outdir)
