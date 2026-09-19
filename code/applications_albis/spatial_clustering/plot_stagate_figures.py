#!/usr/bin/env python
"""
Figure 4B plots -- domain scatter + UMAP -- from an already-written
3D_stagate.py output. Pure CPU: reads the saved h5ad (obsm['STAGATE'],
obs['mclust_3d'], obs['domain_true'], obs['slice_id']) and only needs
scanpy/matplotlib/sklearn, not torch/STAGATE_pyG/R -- so replotting (e.g.
after a title/format tweak) doesn't need a GPU allocation.

3D_stagate.py imports these same functions to plot right after training;
this script lets you re-run just the plotting step against an existing
adata_results/Sim_3D_STAGATE_<dataset>.h5ad.

Usage:
    python plot_stagate_figures.py --dataset cell_strongmix_bs0_seed2025

Output, under sim_paper/data/figure_4/spatial_clustering/STAGATE/<run_tag>/plots/ :
    domains_3d_true_vs_stagate.png
    umap_stagate3d.png            domain_true + STAGATE 3D Domains
    umap_stagate3d_by_slice.png   + Slice ID
"""
import argparse
import os

import numpy as np
import scanpy as sc
import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

BASE_OUTDIR = "/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_4/spatial_clustering/STAGATE"


# --------------------------------------------------------------------------- #
# plotting                                                                   #
# --------------------------------------------------------------------------- #
def plot_3d_panels(adata, z, outdir):
    # Each column gets its own color LUT, built from its own label set --
    # domain_true is "D0".."D5"/"unassigned" but mclust_3d/mclust_2d are R
    # mclust's own cluster IDs ("1".."6", no relation to domain_true's
    # naming), so a single shared LUT keyed on domain_true's labels left
    # every mclust point unmatched and grey.
    cols = [("domain_true", "True Domains"),
            ("mclust_3d", "STAGATE-3D"),
            ("mclust_2d", "STAGATE-2D")]
    fig = plt.figure(figsize=(13, 4.5))
    domain_lut = None
    for i, (key, title) in enumerate(cols):
        ax = fig.add_subplot(1, 3, i + 1, projection="3d")
        vals = adata.obs[key].astype(str).values
        labels = sorted(np.unique(vals))
        lut = {lab: plt.cm.tab10(j % 10) for j, lab in enumerate(labels)}
        if key == "domain_true":
            domain_lut = lut
        for lab in labels:
            m = vals == lab
            ax.scatter(adata.obsm["spatial"][m, 0], adata.obsm["spatial"][m, 1],
                       z[m], s=0.5, marker="o", label=lab,
                       color=lut[lab])
        ax.set_title(title)
        ax.set_xticklabels([]); ax.set_yticklabels([]); ax.set_zticklabels([])
        ax.elev = 15; ax.azim = -60
        ax.grid(False)
        ax.set_axis_off()
        if key != "domain_true":
            # mclust cluster IDs aren't aligned to domain_true's identity, so
            # each STAGATE panel needs its own small legend rather than
            # sharing the bottom domain_true one.
            handles = [Line2D([0], [0], marker="o", ls="", mfc=lut[lab], mec="none", label=lab)
                       for lab in labels]
            ax.legend(handles=handles, loc="upper right", fontsize=6, frameon=False,
                      title=key, title_fontsize=6, markerscale=1.5)
    handles = [Line2D([0], [0], marker="o", ls="", mfc=domain_lut[lab], mec="none", label=lab)
               for lab in sorted(domain_lut)]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.0),
               borderaxespad=0.2, ncol=min(len(handles), 10), frameon=False,
               title="domain_true")
    fig.subplots_adjust(left=0.02, right=0.98, wspace=0.05)
    fig.tight_layout(rect=[0, 0.08, 1, 1])
    fig.savefig(os.path.join(outdir, "domains_3d_true_vs_stagate.png"), dpi=200)
    plt.close(fig)


def plot_umap(adata, outdir):
    try:
        sc.pp.neighbors(adata, use_rep="STAGATE")
        sc.tl.umap(adata)

        fig = sc.pl.umap(adata, color=["domain_true", "mclust_3d"],
                         title=["Domain True", "STAGATE 3D Domains"],
                         show=False, return_fig=True)
        fig.savefig(os.path.join(outdir, "umap_stagate3d.png"), dpi=200, bbox_inches="tight")
        plt.close(fig)

        fig = sc.pl.umap(adata, color=["domain_true", "mclust_3d", "slice_id"],
                         title=["Domain True", "STAGATE 3D Domains", "Slice ID"],
                         show=False, return_fig=True)
        fig.savefig(os.path.join(outdir, "umap_stagate3d_by_slice.png"), dpi=200, bbox_inches="tight")
        plt.close(fig)
    except Exception as e:  # noqa: BLE001
        print("[warn] UMAP plot skipped:", e)


# --------------------------------------------------------------------------- #
# main                                                                       #
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True,
                    help="base dataset name, e.g. cell_strongmix (matches 3D_stagate.py's --dataset)")
    args = ap.parse_args()

    run_tag = args.dataset
    outdir = os.path.join(BASE_OUTDIR, run_tag)
    out_adata = os.path.join(outdir, "adata_results")
    out_plots = os.path.join(outdir, "plots")
    os.makedirs(out_plots, exist_ok=True)

    h5ad_path = os.path.join(out_adata, f"Sim_3D_STAGATE_{args.dataset}.h5ad")
    print(f"[load] {h5ad_path}")
    adata = ad.read_h5ad(h5ad_path)
    print(f"[load] {adata.shape}")

    z = (np.asarray(adata.obsm["spatial_3d"])[:, 2]
         if "spatial_3d" in adata.obsm else adata.obs["slice_id"].astype(float).values)
    plot_3d_panels(adata, z, out_plots)
    plot_umap(adata, out_plots)
    print("[done]", out_plots)


if __name__ == "__main__":
    main()
