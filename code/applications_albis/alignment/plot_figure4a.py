"""
Figure 4C -- before/after STAIR 3D alignment of the ALBIS sphere z-stack.

Reads the STAIR outputs written by 3D_stair.py
    sim_paper/data/figure_4/alignment/STAIR/<dataset>/adata_results/
        Sim_3D_STAIR_<dataset>.h5ad
        metrics.json

Produces, under sim_paper/data/figure_4/alignment/plots/ :
    figure4c_overlay_2d_<colorby>_<dataset>.png  one modality per image, cols =
                                        [unaligned | STAIR | truth], all 10 slices overlaid, 2D
    figure4c_sphere_3d_<colorby>_<dataset>.png   same columns, 3D scatter (the reconstructed sphere)
    figure4c_rmse.png                   joint-Procrustes RMSE, unaligned vs STAIR, per modality
    figure4c_metrics.csv                the numbers behind the bar chart
    figure4c_mnn_domain_agreement.csv    domain-recovery quantification, unaligned vs STAIR: for
                                        every pair of the 10 slices, cross-slice mutual nearest
                                        neighbors (in the 2D spatial embedding) and what fraction
                                        share the same domain_true label. Skip with
                                        --no-mnn-metrics. Plotted by plot_stair_mnn_accuracy.py
                                        (stair_mnn_accuracy.png) -- not plotted from this script.

The STAIR result is anchored on slice 0, so for a shared frame with the
ground-truth column we rigidly Procrustes-rotate the STAIR stack onto the truth
for *display only* (disable with --no-procrustes). RMSE numbers already use a
joint Procrustes fit and are unaffected.
"""

import argparse
import itertools
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
from sklearn.neighbors import NearestNeighbors

BASE = "/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_4/alignment"
ALL_DATASETS = ["bin16um", "spot", "cell"]
FIG4 = "/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_4"
DATASET_CHOICES = ALL_DATASETS + ["bin16um_realwindow", "spot_realwindow", "cell_r6000"]
COL_LABELS = ["Unaligned", "STAIR aligned", "Ground truth"]

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from manuscript_style import category_color, category_order


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
    cats = category_order(adata.obs[color_by].astype(str))
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
    z = np.asarray(adata.obsm["spatial_3d_true"])[:, 2]
    return {
        "2d": [unal_xy, fine_xy, true_xy],
        "3d": [np.column_stack([unal_xy, z]),
               np.column_stack([fine_xy, z]),
               np.column_stack([true_xy, z])],
    }


def plot_overlay_2d(data, color_by, do_procrustes, outdir, tag):
    for ds, (adata, _) in data.items():
        col, cats, lut = color_vec(adata, color_by)
        idx = subsample(adata.n_obs, 60000)
        cols3 = panels_for(adata, do_procrustes)["2d"]
        s = 1.0 if adata.n_obs > 20000 else 6.0
        fig, axes = plt.subplots(1, 3, figsize=(12, 4.3), squeeze=False)
        for c in range(3):
            ax = axes[0][c]
            xy = cols3[c][idx]
            ax.scatter(xy[:, 0], xy[:, 1], s=s, c=col[idx], linewidths=0, alpha=0.85)
            ax.set_aspect("equal")
            ax.axis("off")
            ax.set_title(COL_LABELS[c], fontsize=13)
        handles = [Line2D([0], [0], marker="o", ls="", mfc=lut[c], mec="none",
                          label=str(c)) for c in cats]
        fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.0),
                   borderaxespad=0.2, ncol=min(len(cats), 10), frameon=False,
                   title={"domain_true": "Domain True", "cell_type_true": "Cell Type",
                          "slice_id": "Slice ID"}.get(color_by, color_by))
        fig.patch.set_facecolor("white")
        fig.tight_layout(rect=[0, 0.035, 1, 0.94])
        outpath = os.path.join(outdir, f"figure4c_overlay_2d_{tag}_{ds}.png")
        fig.savefig(outpath, dpi=200, facecolor="white")
        plt.close(fig)
        print("wrote", outpath)


def plot_sphere_3d(data, color_by, do_procrustes, outdir, tag):
    for ds, (adata, _) in data.items():
        col, cats, lut = color_vec(adata, color_by)
        idx = subsample(adata.n_obs, 40000)
        cols3 = panels_for(adata, do_procrustes)["3d"]
        fig = plt.figure(figsize=(12, 4.3))
        for c in range(3):
            ax = fig.add_subplot(1, 3, c + 1, projection="3d")
            xyz = cols3[c][idx]
            ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], s=1.5, c=col[idx],
                       linewidths=0, alpha=0.75)
            ax.set_axis_off()
            ax.view_init(elev=12, azim=-60)
            ax.set_box_aspect(None, zoom=1.3)  # crop the wide default 3D margins
            ax.set_title(COL_LABELS[c], fontsize=13, y=0.95)
        handles = [Line2D([0], [0], marker="o", ls="", mfc=lut[c], mec="none",
                          label=str(c)) for c in cats]
        fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.0),
                   borderaxespad=0.2, ncol=min(len(cats), 10), frameon=False,
                   title={"domain_true": "Domain True", "cell_type_true": "Cell Type",
                          "slice_id": "Slice ID"}.get(color_by, color_by))
        fig.patch.set_facecolor("white")
        fig.subplots_adjust(left=0.01, right=0.99, top=0.92, bottom=0.08, wspace=0.0)
        outpath = os.path.join(outdir, f"figure4c_sphere_3d_{tag}_{ds}.png")
        fig.savefig(outpath, dpi=200, facecolor="white")
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


def _mutual_nn_domain_agreement(xa, da, xb, db):
    """Fraction of cross-set mutual-nearest-neighbor pairs (point i in xa whose
    nearest neighbor in xb is j, and vice versa) whose domain labels agree."""
    jb = NearestNeighbors(n_neighbors=1).fit(xb).kneighbors(xa, return_distance=False)[:, 0]
    ia = NearestNeighbors(n_neighbors=1).fit(xa).kneighbors(xb, return_distance=False)[:, 0]
    mutual = ia[jb] == np.arange(len(xa))
    if not mutual.any():
        return np.nan, 0
    agree = (da[mutual] == db[jb[mutual]]).mean()
    return agree, int(mutual.sum())


def cross_slice_mnn_domain_agreement(data, max_points_per_slice=None, seed=0):
    """For every pair of slices, find cross-slice mutual nearest neighbors in
    the 2D spatial embedding (points that are each other's closest match
    across the two slices) and report what fraction share a domain_true label
    -- unaligned vs STAIR-aligned. Plain percent agreement, no clustering step
    and no chance-adjustment. One row per (dataset, method, slice pair).

    max_points_per_slice: optionally cap each section to this many observations
    (fixed seed) before the nearest-neighbour search. None (default) uses every
    observation -- the metric is a fraction and is insensitive to the cap."""
    rows = []
    for ds, (adata, _) in data.items():
        domain = adata.obs["domain_true"].astype(str).values
        slice_id = adata.obs["slice_id"].astype(str).values
        slices = sorted(adata.obs["slice_id"].unique(), key=float)
        per_slice_idx = {}
        for sl in slices:
            full_idx = np.where(slice_id == str(sl))[0]
            if max_points_per_slice is None or len(full_idx) <= max_points_per_slice:
                per_slice_idx[sl] = full_idx
            else:
                per_slice_idx[sl] = full_idx[subsample(len(full_idx), max_points_per_slice, seed=seed)]
        panels = panels_for(adata, do_procrustes=True)["2d"]
        for method, xy in zip(["unaligned", "stair"], panels[:2]):
            for a, b in itertools.combinations(slices, 2):
                ia, ib = per_slice_idx[a], per_slice_idx[b]
                agree, n_mnn = _mutual_nn_domain_agreement(xy[ia], domain[ia], xy[ib], domain[ib])
                if n_mnn == 0:
                    continue
                rows.append(dict(dataset=ds, method=method, slice_a=a, slice_b=b,
                                  agreement=agree, n_mnn=n_mnn))
    return pd.DataFrame(rows)


def write_mnn_metrics(data, outdir):
    """Compute cross_slice_mnn_domain_agreement and write it to CSV. The
    actual plot (with chance-level reference line) lives in the sibling
    script plot_stair_mnn_accuracy.py, which reads this CSV -- kept separate
    so there's a single source of truth for the final figure."""
    mnn_df = cross_slice_mnn_domain_agreement(data)
    csv_path = os.path.join(outdir, "figure4c_mnn_domain_agreement.csv")
    mnn_df.to_csv(csv_path, index=False)
    print("wrote", csv_path)
    print((100 * mnn_df.groupby(["dataset", "method"])["agreement"].mean()).round(1))


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
    ap.add_argument("--no-mnn-metrics", action="store_true",
                    help="skip the cross-slice MNN domain-agreement CSV (nearest-neighbor "
                         "search over all slice pairs, slower than the other plots)")
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
    plot_overlay_2d(data, args.color_by, do_pro, args.outdir, tag)
    plot_sphere_3d(data, args.color_by, do_pro, args.outdir, tag)
    plot_rmse(data, args.outdir)
    if not args.no_mnn_metrics:
        write_mnn_metrics(data, args.outdir)
