"""
Figure 4E -- before/after STAIR cross-technology alignment (bin16um + spot +
cell of the same physical ALBIS slice_id).

Reads the output of cross_tech_stair.py
    sim_paper/data/figure_4/cross_modality_alignment/STAIR/cross_tech/slice_<n>/adata_results/
        Sim_CrossTech_STAIR_slice_<n>.h5ad
        metrics.json

Produces, under sim_paper/data/figure_4/cross_modality_alignment/plots/ :
    figure4e_overlay_2d_<colorby>_slice_<n>.png   cols = [unaligned | STAIR | truth],
                                                   all 3 technologies overlaid, 2D
    figure4e_rmse_slice_<n>.png                   joint-Procrustes RMSE bar, unaligned
                                                   vs stair_init vs stair_fine

Mirrors plot_figure4a.py's conventions (same Procrustes-for-display helper,
same panel layout) but colors by `technology` instead of `slice_id`, and has
no 3D sphere panel (cross-tech is a single 2D slice, not a z-stack).
"""

import argparse
import json
import os

import numpy as np
import pandas as pd
import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D

BASE = "/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_4/cross_modality_alignment"
COL_LABELS = ["Unaligned", "STAIR aligned", "Ground truth"]
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from manuscript_style import MODALITY_LOOKUP, category_color, save_figure
TECH_COLORS = {key: MODALITY_LOOKUP[key] for key in ('bin16um', 'spot', 'cell')}


def procrustes_fit(X, Y):
    """Best rigid map X -> Y (rotation + translation, no scale/reflection)."""
    X = np.asarray(X, float)
    Y = np.asarray(Y, float)
    mx, my = X.mean(0), Y.mean(0)
    U, _, Vt = np.linalg.svd((X - mx).T @ (Y - my))
    R = U @ Vt
    if np.linalg.det(R) < 0:
        U[:, -1] *= -1
        R = U @ Vt
    return (X - mx) @ R + my


def load(slice_id):
    outdir = os.path.join(BASE, "STAIR", "cross_tech", f"slice_{slice_id}")
    p = os.path.join(outdir, "adata_results", f"Sim_CrossTech_STAIR_slice_{slice_id}.h5ad")
    a = ad.read_h5ad(p)
    mp = os.path.join(outdir, "adata_results", "metrics.json")
    m = json.load(open(mp)) if os.path.exists(mp) else {}
    return a, m


def color_vec(adata, color_by):
    if color_by == "technology":
        cats = ["bin16um", "spot", "cell"]
        lut = {c: mcolors.to_rgba(TECH_COLORS[c]) for c in cats}
    else:
        cats = list(pd.Categorical(adata.obs[color_by]).categories)
        lut = {c: mcolors.to_rgba(category_color(c, color_by)) for c in cats}
    colors = adata.obs[color_by].astype(str).map({str(k): v for k, v in lut.items()}).values
    return np.vstack(colors), cats, lut


def subsample(n, max_points, seed=0):
    if n <= max_points:
        return np.arange(n)
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(n, max_points, replace=False))


def panels_for(adata, do_procrustes):
    true_xy = np.asarray(adata.obsm["spatial_true"], float)
    unal_xy = np.asarray(adata.obsm["spatial"], float)
    fine_xy = np.asarray(adata.obsm["transform_fine"], float)
    if do_procrustes:
        fine_xy = procrustes_fit(fine_xy, true_xy)
    return [unal_xy, fine_xy, true_xy]


def plot_overlay_2d(adata, color_by, do_procrustes, outdir, slice_id):
    col, cats, lut = color_vec(adata, color_by)
    cols3 = panels_for(adata, do_procrustes)
    tech = adata.obs["technology"].astype(str).values
    # per-technology subsample cap (NOT a pooled cap across all 3 technologies --
    # a pooled subsample.min(80000, n_obs) applied then intersected per-tech keeps
    # the same ~fraction of every technology's rows, which silently drops most of
    # spot's already-small point count while being invisible for bin16um/cell)
    per_tech_idx = {t: np.where(tech == t)[0][subsample(int((tech == t).sum()), 60000)]
                    for t in TECH_COLORS}
    # draw largest-n technology first, sparse spot points last so they aren't buried
    plot_order = ["cell", "bin16um", "spot"]
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.3), squeeze=False)
    for c in range(3):
        ax = axes[0][c]
        xy = cols3[c]
        for t in plot_order:
            sel = per_tech_idx[t]
            colors = [lut[t]] * len(sel) if color_by == "technology" else col[sel]
            if t == "spot":
                # spot is ~40-100x sparser per unit area than bin16um/cell (Visium-pitch
                # grid vs. bin/cell resolution) -- small semi-transparent markers blend
                # into the dense point clouds underneath and read as "missing" even
                # though every point is there. Bigger, opaque, outlined markers fix that.
                ax.scatter(xy[sel, 0], xy[sel, 1], s=14, c=colors, linewidths=0.3,
                           edgecolors="black", alpha=1.0)
            else:
                s = 1.0 if (tech == t).sum() > 20000 else 8.0
                ax.scatter(xy[sel, 0], xy[sel, 1], s=s, c=colors, linewidths=0, alpha=0.85)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(COL_LABELS[c], fontsize=13)
    handles = [Line2D([0], [0], marker="o", ls="", mfc=lut[c], mec="none",
                      label=str(c)) for c in cats]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.0),
               borderaxespad=0.2, ncol=min(len(cats), 10), frameon=False,
               title="Domain True" if color_by == "domain_true" else color_by)
    fig.patch.set_facecolor("white")
    fig.tight_layout(rect=[0, 0.035, 1, 0.94])
    outpath = os.path.join(outdir, f"figure4e_overlay_2d_{color_by}_slice_{slice_id}.png")
    save_figure(fig, outpath, dpi=500, facecolor="white")
    plt.close(fig)
    print("wrote", outpath)


def plot_rmse(metrics, outdir, slice_id):
    if not metrics:
        print("no metrics.json found, skipping RMSE bar")
        return
    row = dict(unaligned=metrics.get("unaligned_rmse_um", np.nan),
               stair_init=metrics.get("stair_init_rmse_um", np.nan),
               stair_fine=metrics.get("stair_fine_rmse_um", np.nan))
    df = pd.DataFrame([row], index=[f"slice_{slice_id}"])
    ax = df.plot.bar(figsize=(5, 4), color=["#B9C0C7", "#7393B3", "#4C8FD5"])
    ax.set_ylabel("joint-Procrustes RMSE to truth (um)")
    ax.set_xlabel("")
    ax.set_title("Figure 4E -- cross-tech alignment error\nbefore vs after STAIR")
    plt.xticks(rotation=0)
    plt.tight_layout()
    p = os.path.join(outdir, f"figure4e_rmse_slice_{slice_id}.png")
    save_figure(plt.gcf(), p, dpi=500)
    plt.close()
    print("wrote", p)
    print(df)


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--slice", type=int, default=5)
    ap.add_argument("--color-by", nargs="+", default=["technology", "domain_true"],
                    choices=["technology", "domain_true", "cell_type_true"])
    ap.add_argument("--no-procrustes", action="store_true")
    ap.add_argument("--outdir", default=os.path.join(BASE, "plots"))
    return ap.parse_args()


if __name__ == "__main__":
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    do_pro = not args.no_procrustes

    adata, metrics = load(args.slice)
    print(f"loaded slice_{args.slice}: {adata.shape}")

    for cb in args.color_by:
        plot_overlay_2d(adata, cb, do_pro, args.outdir, args.slice)
    plot_rmse(metrics, args.outdir, args.slice)
