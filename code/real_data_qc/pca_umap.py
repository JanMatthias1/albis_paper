#!/usr/bin/env python
"""
PCA + UMAP for a QC'd real slice/sample, as a rough "how crazy is our data"
reference point next to the sim pipeline's step01_pca_harmony.py plots.

Each real dataset here is a single sample (no batch/slice_id column), so
there is nothing for Harmony to correct -- this only runs normalize/log/scale
-> PCA (--n-pcs, default 30, matching step01_pca_harmony.py's default) -> neighbors
-> Leiden (for an unsupervised coloring, since there's no cell_type_true/
domain_true ground truth) -> UMAP. low_qc-flagged cells/bins are dropped by
default (use --keep-low-qc to keep them). Visium HD's whole-transcriptome
panel (18,085 genes) is subset to the top --n-hvg highly variable genes
before scaling/PCA to keep the dense step tractable; Xenium's ~400-gene
panel is small enough to use as-is, matching how step01_pca_harmony.py uses all
genes for the sim data.

Expected environment:
    conda activate /dcs04/hicks/data/Jan/sim_project/albis/env/albis-tutorial
    pip install scikit-misc   # one-time; not yet in requirements -- needed by
                              # sc.pp.highly_variable_genes(flavor="seurat_v3")

Example:
    python sim_paper/code/real_data_qc/pca_umap.py --dataset non_diseased_lung
    python sim_paper/code/real_data_qc/pca_umap.py --dataset breast_cancer_visium_hd
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

SCRIPT_DIR = Path(__file__).resolve().parent
SIM_PAPER_DIR = SCRIPT_DIR.parents[1]
REAL_DATA_QC_ROOT = SIM_PAPER_DIR / "data" / "real_data_qc"

VALID_DATASETS = (
    "human_pancreas_visium_hd",
    "breast_cancer_visium_hd",
    "non_diseased_lung",
    "lung_cancer",
)
EXTRA_COLOR_CANDIDATES = ["pct_counts_mt", "count_density"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run PCA, Leiden, and UMAP on a QC'd real dataset (no Harmony: single sample)."
    )
    parser.add_argument("--dataset", choices=VALID_DATASETS, required=True)
    parser.add_argument("--input", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--resolution", type=float, default=0.5)
    parser.add_argument("--target-sum", type=float, default=1e4)
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument(
        "--hvg-threshold",
        type=int,
        default=3000,
        help="Subset to --n-hvg genes if the dataset has more genes than this (whole-transcriptome Visium HD panels).",
    )
    parser.add_argument("--n-hvg", type=int, default=2000)
    parser.add_argument("--keep-low-qc", action="store_true", help="Skip dropping obs['low_qc'] cells/bins.")
    parser.add_argument("--point-size", type=float, default=2.0)
    parser.add_argument("--alpha", type=float, default=0.5)
    parser.add_argument(
        "--plot-colors",
        nargs="+",
        default=None,
        help="obs columns to color by. Defaults to leiden + total_counts + n_genes_by_counts, "
        "plus pct_counts_mt/count_density when present.",
    )
    args = parser.parse_args()

    if args.input is None:
        args.input = REAL_DATA_QC_ROOT / args.dataset / f"{args.dataset}_qc.h5ad"
    if args.output_dir is None:
        args.output_dir = REAL_DATA_QC_ROOT / args.dataset
    return args


def color_values(obs: pd.DataFrame, key: str):
    values = obs[key]
    if pd.api.types.is_numeric_dtype(values):
        return values.to_numpy(), None

    categories = pd.Categorical(values.astype(str))
    return categories.codes, list(categories.categories)


def plot_two_dims(adata, embedding_key: str, color_key: str, output_path: Path, point_size: float, alpha: float) -> None:
    coords = np.asarray(adata.obsm[embedding_key])
    colors, categories = color_values(adata.obs, color_key)

    fig, ax = plt.subplots(figsize=(6, 5))
    scatter = ax.scatter(
        coords[:, 0],
        coords[:, 1],
        c=colors,
        s=point_size,
        linewidths=0,
        cmap="tab20" if categories is not None else "viridis",
        alpha=alpha,
    )
    ax.set_xlabel(f"{embedding_key} 1")
    ax.set_ylabel(f"{embedding_key} 2")
    ax.set_title(f"{embedding_key}: {color_key}")

    if categories is None:
        fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04, label=color_key)
    else:
        handles = [
            plt.Line2D(
                [0], [0], marker="o", color="w",
                markerfacecolor=scatter.cmap(scatter.norm(i)), markersize=6, label=label,
            )
            for i, label in enumerate(categories)
        ]
        ax.legend(handles=handles, title=color_key, bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)

    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def resolve_plot_colors(adata, requested: list[str] | None) -> list[str]:
    if requested is not None:
        return requested
    colors = ["leiden", "total_counts", "n_genes_by_counts"]
    colors += [c for c in EXTRA_COLOR_CANDIDATES if c in adata.obs]
    return colors


def main() -> None:
    args = parse_args()
    if not args.input.is_file():
        raise SystemExit(f"Input file not found: {args.input}\nRun xenium_qc.py / visium_hd_qc.py first.")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f"[load] {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"[load] AnnData shape: {adata.n_obs} obs x {adata.n_vars} genes")

    if not args.keep_low_qc and "low_qc" in adata.obs:
        n_before = adata.n_obs
        adata = adata[~adata.obs["low_qc"]].copy()
        print(f"[qc] Dropped {n_before - adata.n_obs} low_qc obs, {adata.n_obs} remaining")

    adata.var_names_make_unique()
    adata.layers["counts"] = adata.X.copy()

    hvg_selected = adata.n_vars > args.hvg_threshold
    if hvg_selected:
        print(f"[hvg] {adata.n_vars} genes > {args.hvg_threshold}; selecting top {args.n_hvg} highly variable genes")
        sc.pp.highly_variable_genes(adata, n_top_genes=args.n_hvg, flavor="seurat_v3", layer="counts")
        adata = adata[:, adata.var["highly_variable"]].copy()
        print(f"[hvg] Subset to {adata.n_vars} genes")

    print(f"[pca] Normalizing, scaling, and running PCA with {args.n_pcs} components")
    sc.pp.normalize_total(adata, target_sum=args.target_sum)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(
        adata,
        n_comps=args.n_pcs,
        use_highly_variable=False,
        svd_solver="arpack",
        random_state=args.random_state,
    )

    print(f"[neighbors] n_neighbors={args.n_neighbors} on X_pca")
    sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, use_rep="X_pca", random_state=args.random_state)

    print(f"[leiden] resolution={args.resolution} (unsupervised coloring; no ground-truth labels here)")
    sc.tl.leiden(
        adata,
        resolution=args.resolution,
        key_added="leiden",
        random_state=args.random_state,
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )

    print("[umap] Computing UMAP")
    sc.tl.umap(adata, random_state=args.random_state)

    plot_colors = resolve_plot_colors(adata, args.plot_colors)
    valid_colors = [key for key in plot_colors if key in adata.obs]
    missing = sorted(set(plot_colors) - set(valid_colors))
    if missing:
        print(f"[plot] Skipping missing obs columns: {', '.join(missing)}")

    for color_key in valid_colors:
        plot_two_dims(adata, "X_pca", color_key, args.output_dir / f"pca_by_{color_key}.png", args.point_size, args.alpha)
        plot_two_dims(adata, "X_umap", color_key, args.output_dir / f"umap_by_{color_key}.png", args.point_size, args.alpha)

    summary = {
        "dataset": args.dataset,
        "input": str(args.input),
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "n_pcs": int(args.n_pcs),
        "n_neighbors": int(args.n_neighbors),
        "resolution": float(args.resolution),
        "hvg_selected": hvg_selected,
        "dropped_low_qc": not args.keep_low_qc,
        "n_leiden_clusters": int(adata.obs["leiden"].nunique()),
        "plot_colors": valid_colors,
    }
    with open(args.output_dir / f"{args.dataset}_pca_umap_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[save] Summary written to {args.output_dir / f'{args.dataset}_pca_umap_summary.json'}")

    adata.uns["pca_umap"] = summary
    out_h5ad = args.output_dir / f"{args.dataset}_pca_umap.h5ad"
    adata.write_h5ad(out_h5ad)
    print(f"[save] {out_h5ad}")


if __name__ == "__main__":
    main()
