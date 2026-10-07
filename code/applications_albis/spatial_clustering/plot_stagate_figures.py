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

Output, under albis_paper/data/figure_4/spatial_clustering/STAGATE/<run_tag>/plots/ :
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

BASE_OUTDIR = "/dcs04/hicks/data/Jan/sim_project/albis_paper/data/figure_4/spatial_clustering/STAGATE"


import sys
from pathlib import Path
from matplotlib.transforms import Bbox
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from manuscript_style import (apply_style, save_figure, scatter_colors, matched_labels,
    PANEL_FIGSIZE, PANEL_MARGINS, PANEL_EXPORT_BOTTOM)
import json
apply_style()


def panel_figure(n, projection=None):
    # Match Figure 3 axes in physical inches, including the three-panel view.
    width, height = PANEL_FIGSIZE
    panel_width = width * (PANEL_MARGINS['right'] - PANEL_MARGINS['left']) / (2 + PANEL_MARGINS['wspace'])
    gap = panel_width * PANEL_MARGINS['wspace']
    total_width = width + (n - 2) * (panel_width + gap)
    fig = plt.figure(figsize=(total_width, height))
    axes = [fig.add_axes([(width * PANEL_MARGINS['left'] + i * (panel_width + gap)) / total_width,
                         PANEL_MARGINS['bottom'], panel_width / total_width,
                         PANEL_MARGINS['top'] - PANEL_MARGINS['bottom']], projection=projection)
            for i in range(n)]
    return fig, axes


def draw_panel(ax, adata, key, title, coords):
    args, labels, colors = scatter_colors(adata.obs, key,
        'domain_true' if key.startswith('mclust_') else None)
    # Clusters keep their matched domain colours, but the legend shows only the
    # cluster number (no "1 → D0" correspondence).
    labels = [label.split(' → ')[0] for label in labels]
    ax.scatter(*coords.T, **args, s=4, linewidths=0, alpha=.85)
    ax.set_title(title)
    handles = [Line2D([0], [0], marker='o', linestyle='', color=color,
                     markersize=12, label=label) for label, color in zip(labels, colors)]
    heading = 'Cluster label' if key.startswith('mclust_') else ('Slice ID' if key == 'slice_id' else 'Domain True')
    ax.legend(handles=handles, title=heading, loc='upper center',
              bbox_to_anchor=(.5, -.17), ncol=max(1, (len(labels)+1)//2),
              frameon=False, fontsize=12, title_fontsize=13, columnspacing=.7,
              handlelength=1, handletextpad=.3)


def save_panel(fig, outdir, name):
    crop = Bbox.from_extents(0, PANEL_EXPORT_BOTTOM, fig.get_figwidth(), fig.get_figheight())
    save_figure(fig, os.path.join(outdir, name), dpi=500, bbox_inches=crop)
    plt.close(fig)


def plot_3d_panels(adata, z, outdir):
    fig, axes = panel_figure(3, '3d')
    coords = np.column_stack([adata.obsm['spatial'][:, :2], z])
    for ax, key, title in zip(axes, ['domain_true','mclust_3d','mclust_2d'],
                             ['Domain True','STAGATE-3D','STAGATE-2D']):
        draw_panel(ax, adata, key, title, coords)
        ax.view_init(elev=12, azim=-60)
        ax.set_box_aspect(None, zoom=1.65)
        ax.set_axis_off()
    save_panel(fig, outdir, 'domains_3d_true_vs_stagate.png')
    Path(outdir, 'domain_colors.json').write_text(json.dumps({
        key: matched_labels(adata.obs, key, 'domain_true')
        for key in ['mclust_3d','mclust_2d']}, indent=2))


def plot_umap(adata, outdir):
    # Preserve the saved layout; computing a new embedding is unnecessary for styling.
    if 'X_umap' not in adata.obsm:
        sc.pp.neighbors(adata, use_rep='STAGATE', random_state=0)
        sc.tl.umap(adata, random_state=0)
    coords = np.asarray(adata.obsm['X_umap'])
    for keys, titles, name in [
        (['domain_true','mclust_3d'], ['Domain True','STAGATE-3D'], 'umap_stagate3d.png'),
        (['domain_true','mclust_3d','slice_id'], ['Domain True','STAGATE-3D','Slice ID'], 'umap_stagate3d_by_slice.png')]:
        fig, axes = panel_figure(len(keys))
        for ax, key, title in zip(axes, keys, titles):
            draw_panel(ax, adata, key, title, coords)
            ax.set_xlabel('UMAP 1'); ax.set_ylabel('UMAP 2')
            ax.set_xlim(coords[:,0].min()-1, coords[:,0].max()+1)
            ax.set_ylim(coords[:,1].min()-1, coords[:,1].max()+1)
        save_panel(fig, outdir, name)


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
    import h5py
    from anndata._io.specs import read_elem
    with h5py.File(h5ad_path) as f:
        adata = ad.AnnData(obs=read_elem(f['obs']),
                          obsm={key: read_elem(f['obsm'][key]) for key in
                                ('spatial', 'spatial_3d', 'STAGATE', 'X_umap') if key in f['obsm']})
    print(f"[load] {adata.shape}")

    z = (np.asarray(adata.obsm["spatial_3d"])[:, 2]
         if "spatial_3d" in adata.obsm else adata.obs["slice_id"].astype(float).values)
    plot_3d_panels(adata, z, out_plots)
    plot_umap(adata, out_plots)
    print("[done]", out_plots)


if __name__ == "__main__":
    main()
