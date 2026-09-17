"""
Figure 4C -- per-slice translation error, before vs after STAIR.

A rigid perturbation is rotation + translation. The existing figure4c_rmse.png
(joint-Procrustes RMSE) mixes both together; this isolates just the
translation component: for each slice_id, the distance between that slice's
centroid and the corresponding ground-truth centroid.

  unaligned translation error = || centroid(spatial[slice]) - centroid(spatial_true[slice]) ||
  STAIR translation error     = || centroid(transform_fine_procrustes[slice]) - centroid(spatial_true[slice]) ||

STAIR's own reconstruction is anchored on slice 0 in an arbitrary frame (no
fixed relationship to the ground-truth frame), so transform_fine is first
Procrustes-fit (rotation + translation, no scale) as ONE RIGID BODY onto
spatial_true across all slices together -- same "for display" fit
plot_figure4c.py uses -- before computing per-slice centroids. That removes
the one global degree of freedom STAIR's anchor choice introduces and leaves
the genuine per-slice residual. unaligned needs no such fit: each slice's
'spatial'/'spatial_true' already share the same base frame by construction
(only that slice's own local rigid_perturb_inplane offsets it).

Reads the STAIR outputs written by 3D_stair.py
    sim_paper/data/figure_4/alignment/STAIR/<dataset>/adata_results/
        Sim_3D_STAIR_<dataset>.h5ad

Produces, under sim_paper/data/figure_4/alignment/plots/ :
    figure4c_translation_error.csv     dataset, slice_id, unaligned_um, stair_um
    figure4c_translation_error.png     per-modality before/after strip+bar plot
"""

import argparse
import os

import numpy as np
import pandas as pd
import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_4/alignment"
ALL_DATASETS = ["bin16um", "spot", "cell"]
MOD_COLORS = {"bin16um": "#4878d0", "spot": "#ee854a", "cell": "#6acc64"}


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
    return ad.read_h5ad(p)


def per_slice_translation_error(adata):
    true_xy = np.asarray(adata.obsm["spatial_true"], float)
    unal_xy = np.asarray(adata.obsm["spatial"], float)
    fine_xy = procrustes_fit(np.asarray(adata.obsm["transform_fine"], float), true_xy)
    slice_id = adata.obs["slice_id"].astype(str).values

    rows = []
    for sl in sorted(set(slice_id), key=float):
        m = slice_id == sl
        true_c = true_xy[m].mean(0)
        unal_err = np.linalg.norm(unal_xy[m].mean(0) - true_c)
        stair_err = np.linalg.norm(fine_xy[m].mean(0) - true_c)
        rows.append(dict(slice_id=sl, n_obs=int(m.sum()),
                          unaligned_um=unal_err, stair_um=stair_err))
    return pd.DataFrame(rows)


def plot(df, outdir):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    datasets = df["dataset"].unique()
    x = np.arange(len(datasets))
    width = 0.32
    for i, stage in enumerate(["unaligned_um", "stair_um"]):
        means = [df.loc[df["dataset"] == d, stage].mean() for d in datasets]
        offset = (i - 0.5) * width
        color = "#797979" if stage == "unaligned_um" else "#4878d0"
        ax.bar(x + offset, means, width=width, color=color,
               label="Unaligned" if stage == "unaligned_um" else "STAIR aligned")
    for i, d in enumerate(datasets):
        sub = df[df["dataset"] == d]
        for stage, offset in [("unaligned_um", -0.5 * width), ("stair_um", 0.5 * width)]:
            jitter = (np.random.default_rng(0).random(len(sub)) - 0.5) * width * 0.5
            ax.scatter(np.full(len(sub), i + offset) + jitter, sub[stage],
                       s=14, c="black", alpha=0.6, zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels(datasets)
    ax.set_ylabel("per-slice centroid translation error (um)")
    ax.set_title("Figure 4C -- translation error, unaligned vs STAIR\n(dots = individual slices)")
    ax.legend(frameon=False)
    fig.tight_layout()
    p = os.path.join(outdir, "figure4c_translation_error.png")
    fig.savefig(p, dpi=200)
    plt.close(fig)
    print("wrote", p)


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="+", default=ALL_DATASETS, choices=ALL_DATASETS)
    ap.add_argument("--outdir", default=os.path.join(BASE, "plots"))
    return ap.parse_args()


if __name__ == "__main__":
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    all_rows = []
    for ds in args.datasets:
        adata = load(ds)
        df = per_slice_translation_error(adata)
        df["dataset"] = ds
        all_rows.append(df)
        print(f"{ds}: unaligned mean={df['unaligned_um'].mean():.1f}um, "
              f"STAIR mean={df['stair_um'].mean():.1f}um "
              f"({100 * (1 - df['stair_um'].mean() / df['unaligned_um'].mean()):.1f}% reduction)")

    full = pd.concat(all_rows, ignore_index=True)
    csv_path = os.path.join(args.outdir, "figure4c_translation_error.csv")
    full.to_csv(csv_path, index=False)
    print("wrote", csv_path)

    plot(full, args.outdir)
