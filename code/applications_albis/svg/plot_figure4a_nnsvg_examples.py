#!/usr/bin/env python
"""
Figure 4A (illustrative panel) -- what an SVG analysis with this simulator
looks like in practice.

The AUROC/ROC panel (plot_figure4a_nnsvg.py) is the quantitative benchmark;
this one is qualitative: on a single tissue slice of the strong-domain-mix
cell dataset it shows

  * the ground-truth spatial layout (domains / cell types), then
  * a few genes nnSVG ranks as most spatially variable -- each a true
    `marker_typeN` gene whose expression tracks the domain its cell type
    concentrates in, and
  * a few true `noise` genes (no planted spatial signal), which show only
    uniform speckle regardless of nnSVG rank.

So a reader sees that the simulator plants genes with genuine, known spatial
structure and that a standard off-the-shelf SVG method recovers them --
"here is how you would use the package for an SVG study", not a realism
claim.

Expected environment:
    conda activate /dcs04/hicks/data/Jan/sim_project/sim_paper/env/albis-tutorial

Example:
    python sim_paper/code/applications_albis/svg/plot_figure4a_nnsvg_examples.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402

SCRIPT_DIR = Path(__file__).resolve().parent
SIM_PAPER_DIR = SCRIPT_DIR.parents[2]

STRONGMIX_DIR = SIM_PAPER_DIR / "data" / "figure_4" / "spatial_clustering" / "sim_data"
NNSVG_DIR = SIM_PAPER_DIR / "data" / "figure_4" / "svg" / "nnsvg"

# per-dataset: (h5ad, nnsvg gene_results, signal genes, noise genes). Gene
# picks are top-ranked markers of distinct cell types (-> distinct domains)
# plus true noise genes spanning the nnSVG rank range -- read off each run's
# gene_results_*.csv. bin16um is the default: aggregation over ~4 cells per
# 16um bin makes the domain-level trend visible where per-cell NB noise
# buries it (cell strongmix AUROC 0.65 vs bin16um 0.74).
DATASETS = {
    "bin16um": (
        STRONGMIX_DIR / "packing_pf0p04_bin16um_log_mu_-2.5_bsigma07_strongmix" / "simulation_bin_z_qc.h5ad",
        NNSVG_DIR / "bin16um_strongmix" / "gene_results_bin16um_strongmix.csv",
        ["G44", "G358", "G427"],    # marker_type1 (r1), _type6 (r6), _type7 (r10)
        ["G510", "G501", "G529"],   # noise, nnSVG rank 61 / 78 / 140
    ),
    "cell": (
        STRONGMIX_DIR / "log_mu_-2.3_theta_0.40_jitter0.15_bsigma15_strongmix" / "simulation_cell_z_qc.h5ad",
        NNSVG_DIR / "cell_strongmix" / "gene_results_cell_strongmix.csv",
        ["G180", "G241", "G361"],   # marker_type3 (r1), _type4 (r9), _type6 (r10)
        ["G501", "G524", "G556"],   # noise, nnSVG rank 20 / 253 / 556
    ),
}
DEFAULT_OUT = NNSVG_DIR / "figure_4a_nnsvg_examples.png"


def _slice_expr(adata: ad.AnnData, gene: str) -> np.ndarray:
    x = adata[:, gene].X
    return np.asarray(x.todense()).ravel() if hasattr(x, "todense") else np.asarray(x).ravel()


def _cat_scatter(ax, xy, c, *, title="", s=2):
    cats = pd.Categorical(c)
    colors = plt.get_cmap("tab10")(np.linspace(0, 1, 10))
    for i, cat in enumerate(cats.categories):
        m = cats.codes == i
        ax.scatter(xy[m, 0], xy[m, 1], s=s, color=colors[i % 10], label=str(cat), linewidths=0)
    ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=6, markerscale=3,
              frameon=False, handletextpad=0.2)
    ax.set_title(title, fontsize=8)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_aspect("equal")


def _expr_hexbin(ax, xy, c, *, title="", cmap="viridis", gridsize=55, mincnt=4):
    """Mean log-expression over a spatial hex grid -- averages out cell-level
    NB noise so the domain-level spatial trend the SVG test picks up on is
    actually visible (per-cell scatter buries it at ~49% domain purity).
    vmax clipped to the 97th percentile of bin means so a few sparse edge
    bins don't flatten the colour scale."""
    hb = ax.hexbin(xy[:, 0], xy[:, 1], C=c, gridsize=gridsize, reduce_C_function=np.mean,
                   cmap=cmap, mincnt=mincnt, linewidths=0)
    vals = hb.get_array()
    if vals is not None and len(vals):
        hb.set_clim(float(np.nanmin(vals)), float(np.nanpercentile(vals, 97)))
    ax.set_title(title, fontsize=8)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_aspect("equal")
    return hb


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dataset", choices=sorted(DATASETS), default="bin16um",
                   help="which strong-mix modality to illustrate (default: bin16um)")
    p.add_argument("--h5ad", type=Path, default=None, help="override the dataset's h5ad")
    p.add_argument("--nnsvg-results", type=Path, default=None, help="override the dataset's nnSVG CSV")
    p.add_argument("--slice-id", type=int, default=5)
    p.add_argument("--signal-genes", nargs="+", default=None)
    p.add_argument("--noise-genes", nargs="+", default=None)
    p.add_argument("--max-points", type=int, default=80000, help="subsample the slice for plotting speed")
    p.add_argument("--output", type=Path, default=None)
    args = p.parse_args()

    d_h5ad, d_nnsvg, d_sig, d_noise = DATASETS[args.dataset]
    args.h5ad = args.h5ad or d_h5ad
    args.nnsvg_results = args.nnsvg_results or d_nnsvg
    args.signal_genes = args.signal_genes or d_sig
    args.noise_genes = args.noise_genes or d_noise
    args.output = args.output or (DEFAULT_OUT if args.dataset == "bin16um"
                                 else DEFAULT_OUT.with_name(f"figure_4a_nnsvg_examples_{args.dataset}.png"))

    adata = sc.read_h5ad(args.h5ad)
    adata = adata[adata.obs["slice_id"] == args.slice_id].copy()
    print(f"[load] {args.h5ad.name}  slice {args.slice_id}: {adata.n_obs} cells x {adata.n_vars} genes")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    if adata.n_obs > args.max_points:
        idx = np.random.default_rng(0).choice(adata.n_obs, args.max_points, replace=False)
        adata = adata[np.sort(idx)].copy()
        print(f"[subsample] {adata.n_obs} cells for plotting")

    nn = pd.read_csv(args.nnsvg_results).set_index("gene_names")
    xy = np.asarray(adata.obsm["spatial"], dtype=float)

    n_ex = max(len(args.signal_genes), len(args.noise_genes))
    fig, axes = plt.subplots(2, n_ex + 1, figsize=(3.0 * (n_ex + 1), 6.2))

    _cat_scatter(axes[0, 0], xy, adata.obs["domain_true"].astype(str).to_numpy(),
                 title="ground truth: spatial domains")
    _cat_scatter(axes[1, 0], xy, adata.obs["cell_type_true"].astype(str).to_numpy(),
                 title="ground truth: cell types")

    def _panel(ax, gene):
        gclass = str(adata.var.loc[gene, "gene_class"]) if gene in adata.var_names else "?"
        rank = int(nn.loc[gene, "rank"]) if gene in nn.index else -1
        pval = float(nn.loc[gene, "pval"]) if gene in nn.index else float("nan")
        hb = _expr_hexbin(ax, xy, _slice_expr(adata, gene),
                          title=f"{gene} · {gclass}\nnnSVG rank {rank}/{len(nn)}  (p={pval:.1e})")
        cb = fig.colorbar(hb, ax=ax, fraction=0.046, pad=0.02)
        cb.ax.tick_params(labelsize=6)
        cb.set_label("mean log-expr", fontsize=6)

    for j, g in enumerate(args.signal_genes):
        _panel(axes[0, j + 1], g)
    for j, g in enumerate(args.noise_genes):
        _panel(axes[1, j + 1], g)
    for j in range(len(args.signal_genes), n_ex):
        axes[0, j + 1].axis("off")
    for j in range(len(args.noise_genes), n_ex):
        axes[1, j + 1].axis("off")

    axes[0, 1].text(-0.12, 0.5, "recovered SVGs\n(true marker genes)", transform=axes[0, 1].transAxes,
                    rotation=90, va="center", ha="center", fontsize=9, fontweight="bold")
    axes[1, 1].text(-0.12, 0.5, "true noise genes\n(no planted signal)", transform=axes[1, 1].transAxes,
                    rotation=90, va="center", ha="center", fontsize=9, fontweight="bold")

    fig.suptitle(f"Figure 4A — SVG identification in practice: nnSVG on one slice of the "
                 f"strong-domain-mix cell simulation (n={adata.n_obs:,} cells)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[save] {args.output}")


if __name__ == "__main__":
    main()
