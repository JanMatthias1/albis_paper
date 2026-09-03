"""
Figure 4C -- before/after STAIR 3D alignment of the ALBIS sphere z-stack.

Reads the STAIR outputs written by 3D_stair.py
    sim_paper/data/figure_4/alignment/STAIR/<dataset>/adata_results/
        Sim_3D_STAIR_<dataset>.h5ad
        metrics.json

Produces, under sim_paper/data/figure_4/alignment/plots/ :
    figure4c_overlay_2d_<colorby>.png   rows = modality, cols = [unaligned | STAIR | truth]
                                        all 10 slices overlaid, 2D
    figure4c_sphere_3d_<colorby>.png    same columns, 3D scatter (the reconstructed sphere)
    figure4c_rmse.png                   joint-Procrustes RMSE, unaligned vs STAIR, per modality
    figure4c_metrics.csv                the numbers behind the bar chart

The STAIR result is anchored on slice 0, so for a shared frame with the
ground-truth column we rigidly Procrustes-rotate the STAIR stack onto the truth
for *display only* (disable with --no-procrustes). RMSE numbers already use a
joint Procrustes fit and are unaffected.
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
from matplotlib.lines import Line2D

BASE = "/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_4/alignment"
ALL_DATASETS = ["bin16um", "spot", "cell"]
FIG4 = "/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_4"
DATASET_CHOICES = ALL_DATASETS + ["bin16um_realwindow", "spot_realwindow"]
COL_LABELS = ["Unaligned (input)", "STAIR aligned", "Ground truth"]


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


def load(dataset):
    p = os.path.join(BASE, "STAIR", dataset, "adata_results", f"Sim_3D_STAIR_{dataset}.h5ad")
    a = ad.read_h5ad(p)
    m = {}
    mp = os.path.join(BASE, "STAIR", dataset, "adata_results", "metrics.json")
    if os.path.exists(mp):
        m = json.load(open(mp))
    return a, m


def color_vec(adata, color_by):
    if color_by == "slice_id":
        cats = sorted(adata.obs["slice_id"].unique(), key=float)
        cmap = plt.get_cmap("tab10")
        lut = {c: cmap(i % 10) for i, c in enumerate(cats)}
    else:
        cats = list(pd.Categorical(adata.obs[color_by]).categories)
        cmap = plt.get_cmap("tab10" if len(cats) <= 10 else "tab20")
        lut = {c: cmap(i % cmap.N) for i, c in enumerate(cats)}
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
    z = np.asarray(adata.obsm["spatial_3d_true"])[:, 2]
    return {
        "2d": [unal_xy, fine_xy, true_xy],
        "3d": [np.column_stack([unal_xy, z]),
               np.column_stack([fine_xy, z]),
               np.column_stack([true_xy, z])],
    }


def plot_overlay_2d(data, color_by, do_procrustes, outpath):
    nrow = len(data)
    fig, axes = plt.subplots(nrow, 3, figsize=(12, 4 * nrow), squeeze=False)
    for r, (ds, (adata, _)) in enumerate(data.items()):
        col, cats, lut = color_vec(adata, color_by)
        idx = subsample(adata.n_obs, 60000)
        cols3 = panels_for(adata, do_procrustes)["2d"]
        s = 1.0 if adata.n_obs > 20000 else 6.0
        for c in range(3):
            ax = axes[r][c]
            xy = cols3[c][idx]
            ax.scatter(xy[:, 0], xy[:, 1], s=s, c=col[idx], linewidths=0, alpha=0.6)
            ax.set_aspect("equal")
            ax.set_xticks([]); ax.set_yticks([])
            if r == 0:
                ax.set_title(COL_LABELS[c], fontsize=13)
            if c == 0:
                ax.set_ylabel(ds, fontsize=13)
    handles = [Line2D([0], [0], marker="o", ls="", mfc=lut[c], mec="none",
                      label=str(c)) for c in cats]
    fig.legend(handles=handles, loc="lower center", ncol=min(len(cats), 10),
               frameon=False, title=color_by)
    fig.suptitle("Figure 4C -- ALBIS sphere z-stack, before / after STAIR alignment", fontsize=15)
    fig.tight_layout(rect=[0, 0.05, 1, 0.96])
    fig.savefig(outpath, dpi=200)
    plt.close(fig)
    print("wrote", outpath)


def plot_sphere_3d(data, color_by, do_procrustes, outpath):
    nrow = len(data)
    fig = plt.figure(figsize=(12, 4 * nrow))
    for r, (ds, (adata, _)) in enumerate(data.items()):
        col, cats, lut = color_vec(adata, color_by)
        idx = subsample(adata.n_obs, 40000)
        cols3 = panels_for(adata, do_procrustes)["3d"]
        for c in range(3):
            ax = fig.add_subplot(nrow, 3, r * 3 + c + 1, projection="3d")
            xyz = cols3[c][idx]
            ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], s=1.5, c=col[idx],
                       linewidths=0, alpha=0.5)
            ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
            ax.view_init(elev=12, azim=-60)
            if r == 0:
                ax.set_title(COL_LABELS[c], fontsize=13)
            if c == 0:
                ax.text2D(-0.1, 0.5, ds, transform=ax.transAxes, fontsize=13,
                          rotation=90, va="center")
    fig.suptitle("Figure 4C -- reconstructed 10-slice sphere", fontsize=15)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(outpath, dpi=200)
    plt.close(fig)
    print("wrote", outpath)


def plot_rmse(data, outdir):
    rows = []
    for ds, (_, m) in data.items():
        if not m:
            continue
        rows.append(dict(dataset=ds,
                         unaligned=m.get("unaligned_rmse_um", np.nan),
                         stair_init=m.get("stair_init_rmse_um", np.nan),
                         stair_fine=m.get("stair_fine_rmse_um", np.nan),
                         n_obs=m.get("n_obs", np.nan)))
    if not rows:
        print("no metrics.json found, skipping RMSE bar")
        return
    df = pd.DataFrame(rows).set_index("dataset")
    df.to_csv(os.path.join(outdir, "figure4c_metrics.csv"))
    ax = df[["unaligned", "stair_init", "stair_fine"]].plot.bar(figsize=(7, 4))
    ax.set_ylabel("joint-Procrustes RMSE to truth (um)")
    ax.set_xlabel("")
    ax.set_title("Figure 4C -- alignment error before vs after STAIR")
    plt.xticks(rotation=0)
    plt.tight_layout()
    p = os.path.join(outdir, "figure4c_rmse.png")
    plt.savefig(p, dpi=200)
    plt.close()
    print("wrote", p)
    print(df)


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="+", default=ALL_DATASETS, choices=DATASET_CHOICES)
    ap.add_argument("--color-by", default="slice_id",
                    choices=["slice_id", "domain_true", "cell_type_true"])
    ap.add_argument("--no-procrustes", action="store_true",
                    help="show the STAIR result in its own (slice-0) frame")
    ap.add_argument("--base", default=BASE,
                    help="dir holding STAIR/<dataset>/adata_results/ "
                         "(e.g. .../figure_4/alignment_window_sizing); "
                         "--outdir defaults to <base>/plots")
    ap.add_argument("--outdir", default=None)
    return ap.parse_args()


if __name__ == "__main__":
    args = parse_args()
    BASE = args.base
    outdir = args.outdir or os.path.join(BASE, "plots")
    args.outdir = outdir
    os.makedirs(outdir, exist_ok=True)
    do_pro = not args.no_procrustes

    data = {}
    for ds in args.datasets:
        try:
            data[ds] = load(ds)
        except FileNotFoundError as e:
            print("skip", ds, "-", e)
    if not data:
        raise SystemExit("no STAIR outputs found")

    tag = args.color_by
    plot_overlay_2d(data, args.color_by, do_pro,
                    os.path.join(args.outdir, f"figure4c_overlay_2d_{tag}.png"))
    plot_sphere_3d(data, args.color_by, do_pro,
                   os.path.join(args.outdir, f"figure4c_sphere_3d_{tag}.png"))
    plot_rmse(data, args.outdir)
