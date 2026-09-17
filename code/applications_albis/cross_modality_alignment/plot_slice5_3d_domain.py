"""
3D view of the Figure 4E cross-tech alignment at slice 5, colored by
domain_true: Unaligned / STAIR aligned / Ground truth, all 3 technologies
overlaid.

z is a clean artificial per-technology stack offset (bin16um=0, spot=+900,
cell=+1800 um), NOT each technology's true z. Earlier attempt used real z
(bin16um/spot flat at z~205; cell genuinely spans ~0-1200um since cell bins
by z-range rather than cutting a thin optical plane like the imaging
technologies do) -- that put cell's sparse, spread-out points in the same
z-band as bin16um/spot's dense flat wafer, reading as diffuse noise around a
plate rather than 3 coherent objects, even after fixing the aspect-ratio/
tilt issue (matplotlib's default per-axis autoscale was stretching a true
~15:1 xy:z ratio into looking near-1:1, which also read as a spurious
diagonal tilt -- confirmed not real, corr(z,x/y) ~ 0 in the true data).
Since the actual 2D alignment (what STAIR/the RMSE numbers evaluate) is
XY-only and fully shown in the 2D panels already, z here is purely an
organizing device -- like Fig4C's sphere_3d plots, which get their shape
from stacking many z-slices of ONE technology at real depths; we only have
one slice of 3 technologies, so there's no true depth to show, and forcing
one produced clutter instead of clarity. Distinct flat layers, uniform small
markers (no more per-technology size/alpha special-casing -- that was for
the 2D density-occlusion problem, which stacked layers don't have).
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

BASE = "/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_4/cross_modality_alignment"
TECH_COLORS = {"bin16um": "#4878d0", "spot": "#ee854a", "cell": "#6acc64"}
COL_LABELS = ["Unaligned", "STAIR aligned", "Ground truth"]
PASTEL = ["#4878d0", "#ee854a", "#6acc64", "#d65f5f", "#956cb4",
          "#8c613c", "#dc7ec0", "#797979", "#d5bb67", "#82c6e2"]
LAYER_Z = {"bin16um": 0.0, "spot": 900.0, "cell": 1800.0}
LAYER_ORDER = ["bin16um", "spot", "cell"]


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


def plot(adata, slice_id, outdir, do_procrustes=True):
    tech = adata.obs["technology"].astype(str).values
    domain = adata.obs["domain_true"].astype(str).values
    cats = sorted(set(domain))
    lut = {c: mcolors.to_rgba(PASTEL[i % len(PASTEL)]) for i, c in enumerate(cats)}

    true_xy = np.asarray(adata.obsm["spatial_true"], float)
    unal_xy = np.asarray(adata.obsm["spatial"], float)
    fine_xy = np.asarray(adata.obsm["transform_fine"], float)
    if do_procrustes:
        fine_xy = procrustes_fit(fine_xy, true_xy)
    z = np.array([LAYER_Z[t] for t in tech])
    cols3 = [np.column_stack([unal_xy, z]), np.column_stack([fine_xy, z]),
             np.column_stack([true_xy, z])]

    per_tech_idx = {t: np.where(tech == t)[0][subsample(int((tech == t).sum()), 40000)]
                    for t in TECH_COLORS}

    xy_diam = np.ptp(true_xy, axis=0).max()
    z_span = LAYER_Z["cell"] - LAYER_Z["bin16um"]
    box_aspect = (xy_diam, xy_diam, z_span * 1.3)

    fig = plt.figure(figsize=(13, 4.6))
    for c in range(3):
        ax = fig.add_subplot(1, 3, c + 1, projection="3d")
        for t in LAYER_ORDER:
            sel = per_tech_idx[t]
            colors = [lut[d] for d in domain[sel]]
            xyz = cols3[c][sel]
            ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], s=2.0, c=colors,
                       linewidths=0, alpha=0.85)
        ax.set_axis_off()
        ax.set_box_aspect(box_aspect, zoom=1.35)
        ax.view_init(elev=18, azim=-60)
        ax.set_title(COL_LABELS[c], fontsize=13, y=0.95)

    handles = [Line2D([0], [0], marker="o", ls="", mfc=lut[c], mec="none", label=c) for c in cats]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.02),
               borderaxespad=0.2, ncol=min(len(cats), 10), frameon=False, title="domain_true")
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(left=0.01, right=0.99, top=0.92, bottom=0.14, wspace=0.0)

    outpath = os.path.join(outdir, f"slice{slice_id}_sphere_3d_domain_true.png")
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
    plot(adata, args.slice, args.outdir)
