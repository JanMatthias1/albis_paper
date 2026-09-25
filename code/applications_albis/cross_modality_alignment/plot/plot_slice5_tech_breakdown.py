"""
Figure 4E, "convince me" diagnostic -- shows each technology's own unaligned
disc standalone (proving each technology's full point count really does form
a complete, non-truncated circle) before showing how the three technologies'
independently-perturbed discs relate to each other once overlaid.

Top row: bin16um / spot / cell, each alone, unaligned coords, own axes.
Bottom row: the existing Figure 4E overlay -- Unaligned / STAIR aligned /
Ground truth, all 3 technologies overlaid, colored by technology. Answers
"why does spot look like it lost half its points in the overlay" -- it
hasn't; two same-radius circles offset by roughly 3/4 of their radius only
partially overlap, which is what the bottom-left panel shows once you've
seen each technology is a complete circle on its own in the top row.

Reads the same cross_tech_stair.py output as plot_cross_tech_stair.py:
    sim_paper/data/figure_4/cross_modality_alignment/STAIR/cross_tech/slice_<n>/
        adata_results/Sim_CrossTech_STAIR_slice_<n>.h5ad
"""

import argparse
import os

import numpy as np
import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.patches import Circle

BASE = "/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_4/cross_modality_alignment"
TECHS = ["bin16um", "spot", "cell"]
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from manuscript_style import MODALITY_LOOKUP, category_color
TECH_COLORS = {key: MODALITY_LOOKUP[key] for key in ('bin16um', 'spot', 'cell')}
COL_LABELS = ["Unaligned", "STAIR aligned", "Ground truth"]


def procrustes_fit(X, Y):
    X = np.asarray(X, float)
    Y = np.asarray(Y, float)
    mx, my = X.mean(0), Y.mean(0)
    U, _, Vt = np.linalg.svd((X - mx).T @ (Y - my))
    R = U @ Vt
    if np.linalg.det(R) < 0:
        U[:, -1] *= -1
        R = U @ Vt
    return (X - mx) @ R + my


def subsample(n, max_points, seed=0):
    if n <= max_points:
        return np.arange(n)
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(n, max_points, replace=False))


def load(slice_id):
    outdir = os.path.join(BASE, "STAIR", "cross_tech", f"slice_{slice_id}")
    p = os.path.join(outdir, "adata_results", f"Sim_CrossTech_STAIR_slice_{slice_id}.h5ad")
    return ad.read_h5ad(p)


def panels_for(adata, do_procrustes):
    true_xy = np.asarray(adata.obsm["spatial_true"], float)
    unal_xy = np.asarray(adata.obsm["spatial"], float)
    fine_xy = np.asarray(adata.obsm["transform_fine"], float)
    if do_procrustes:
        fine_xy = procrustes_fit(fine_xy, true_xy)
    return [unal_xy, fine_xy, true_xy]


def plot(adata, outdir, slice_id, do_procrustes=True):
    tech = adata.obs["technology"].astype(str).values
    unal_xy = np.asarray(adata.obsm["spatial"], float)
    cols3 = panels_for(adata, do_procrustes)
    # per-technology subsample cap for the bottom row (NOT a pooled cap across all
    # 3 technologies -- that keeps the same ~fraction of every technology's rows,
    # which silently drops most of spot's already-small point count while being
    # invisible for bin16um/cell)
    per_tech_idx = {t: np.where(tech == t)[0][subsample(int((tech == t).sum()), 60000)]
                    for t in TECHS}

    fig, axes = plt.subplots(2, 3, figsize=(12, 8.6), squeeze=False)

    # --- top row: each technology alone, unaligned, own axes -------------- #
    for j, t in enumerate(TECHS):
        ax = axes[0][j]
        m = tech == t
        pts = unal_xy[m]
        idx = subsample(pts.shape[0], 60000)
        if t == "spot":
            ax.scatter(pts[idx, 0], pts[idx, 1], s=14, c=TECH_COLORS[t], linewidths=0.3,
                       edgecolors="black", alpha=1.0)
        else:
            s = 1.0 if m.sum() > 20000 else 8.0
            ax.scatter(pts[idx, 0], pts[idx, 1], s=s, c=TECH_COLORS[t], linewidths=0, alpha=0.85)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(f"{t} alone, unaligned\n(n={m.sum()})", fontsize=12)

    # --- bottom row: the 3-stage overlay, colored by technology ----------- #
    # plus a fitted circle outline per technology (centroid + max radius,
    # computed on the SAME points shown), so each technology's disc boundary
    # is explicit rather than implied by scatter density.
    plot_order = ["cell", "bin16um", "spot"]
    for c in range(3):
        ax = axes[1][c]
        xy = cols3[c]
        for t in plot_order:
            m = tech == t
            sel = per_tech_idx[t]
            if t == "spot":
                ax.scatter(xy[sel, 0], xy[sel, 1], s=14, c=TECH_COLORS[t], linewidths=0.3,
                           edgecolors="black", alpha=1.0)
            else:
                s = 1.0 if m.sum() > 20000 else 8.0
                ax.scatter(xy[sel, 0], xy[sel, 1], s=s, c=TECH_COLORS[t], linewidths=0, alpha=0.85)
        for t in TECHS:
            m = tech == t
            pts = xy[m]
            centroid = pts.mean(0)
            radius = np.linalg.norm(pts - centroid, axis=1).max()
            ax.add_patch(Circle(centroid, radius, fill=False, edgecolor=TECH_COLORS[t],
                                linewidth=2.0, linestyle="--", alpha=0.9))
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(COL_LABELS[c], fontsize=12)

    handles = [Line2D([0], [0], marker="o", ls="", mfc=TECH_COLORS[t], mec="none", label=t) for t in TECHS]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.0),
               borderaxespad=0.2, ncol=3, frameon=False, title="technology")
    fig.patch.set_facecolor("white")
    fig.tight_layout(rect=[0, 0.035, 1, 0.96])

    outpath = os.path.join(outdir, f"slice{slice_id}_tech_breakdown.png")
    fig.savefig(outpath, dpi=200, facecolor="white")
    plt.close(fig)
    print("wrote", outpath)


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--slice", type=int, default=5)
    ap.add_argument("--outdir", default=os.path.join(BASE, "plots"))
    return ap.parse_args()


if __name__ == "__main__":
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    adata = load(args.slice)
    print(f"loaded slice_{args.slice}: {adata.shape}")
    plot(adata, args.outdir, args.slice)
