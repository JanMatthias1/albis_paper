#!/usr/bin/env python
"""
Run Leiden or Louvain clustering and UMAP from a Harmony-corrected representation.

This is separate from pca_harmony.py/gene_harmony_umap.py so clustering
algorithm/resolution can be iterated without recomputing normalization, PCA,
or Harmony. Two independent Harmony pipelines exist upstream, and either can
be clustered on via --pipeline:
  - pca_harmony:  Harmony corrects the PCA embedding (pca_harmony.py output).
                  Uses the first --n-pcs dimensions (variance-ordered).
  - gene_harmony: Harmony corrects the full gene matrix directly
                  (gene_harmony_umap.py output). Uses all dimensions as-is,
                  since raw genes have no natural "first N" ordering.

Example:
    python sim_paper/code/clustering/clustering_leiden_louvain.py --resolution 0.5
    python sim_paper/code/clustering/clustering_leiden_louvain.py --algorithm louvain --resolution 1.0
    python sim_paper/code/clustering/clustering_leiden_louvain.py --pipeline gene_harmony --resolution 0.5
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


SCRIPT_DIR = Path(__file__).resolve().parent
SIM_PAPER_DIR = SCRIPT_DIR.parents[1]
CLUSTERING_ROOT = SIM_PAPER_DIR / "data" / "clustering"
DEFAULT_INPUT = CLUSTERING_ROOT / "spot" / "pca_harmony" / "simulation_spot_z_pca_harmony.h5ad"
DEFAULT_OUTPUT_DIR = CLUSTERING_ROOT
VALID_MODALITIES = ("spot", "bin", "cell")
VALID_ALGORITHMS = ("leiden", "louvain")
VALID_PIPELINES = ("pca_harmony", "gene_harmony")
PIPELINE_TAGS = {"pca_harmony": "pca", "gene_harmony": "gene"}
PIPELINE_REP_KEYS = {"pca_harmony": "X_pca_harmony", "gene_harmony": "X_gene_harmony"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run neighbors, Leiden clustering, UMAP, and UMAP plots."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--modality", choices=VALID_MODALITIES, default="spot")
    parser.add_argument("--algorithm", choices=VALID_ALGORITHMS, default="leiden")
    parser.add_argument("--pipeline", choices=VALID_PIPELINES, default="pca_harmony")
    parser.add_argument("--resolution", type=float, default=0.5)
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument(
        "--rep-key",
        default=None,
        help="obsm key to cluster on. Defaults to X_pca_harmony or X_gene_harmony based on --pipeline.",
    )
    parser.add_argument("--cluster-key", default="cluster_label")
    parser.add_argument(
        "--plot-colors",
        nargs="+",
        default=["domain_true", "cell_type_true", "cluster_label"],
        help="obs columns to use for UMAP plots.",
    )
    return parser.parse_args()


def resolution_tag(resolution: float) -> str:
    return str(resolution).replace(".", "p").replace("-", "m")


def color_values(obs: pd.DataFrame, key: str):
    values = obs[key]
    if pd.api.types.is_numeric_dtype(values):
        return values.to_numpy(), None

    categories = pd.Categorical(values.astype(str))
    return categories.codes, list(categories.categories)


def plot_umap(adata, color_key: str, output_path: Path) -> None:
    coords = np.asarray(adata.obsm["X_umap"])
    colors, categories = color_values(adata.obs, color_key)

    fig, ax = plt.subplots(figsize=(6, 5))
    scatter = ax.scatter(
        coords[:, 0],
        coords[:, 1],
        c=colors,
        s=5,
        linewidths=0,
        cmap="tab20" if categories is not None else "viridis",
        alpha=0.85,
    )
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title(f"UMAP: {color_key}")

    if categories is None:
        fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04, label=color_key)
    else:
        handles = [
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor=scatter.cmap(scatter.norm(i)),
                markersize=6,
                label=label,
            )
            for i, label in enumerate(categories)
        ]
        ax.legend(
            handles=handles,
            title=color_key,
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            frameon=False,
        )

    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    pipeline_dirname = "pca_harmony" if args.pipeline == "pca_harmony" else "gene_harmony_umap"
    if args.input == DEFAULT_INPUT:
        args.input = (
            CLUSTERING_ROOT / args.modality / pipeline_dirname / f"simulation_{args.modality}_z_{pipeline_dirname}.h5ad"
        )
    if args.rep_key is None:
        args.rep_key = PIPELINE_REP_KEYS[args.pipeline]

    if not args.input.is_file():
        raise SystemExit(
            f"Input file not found: {args.input}\n"
            f"Run sim_paper/code/clustering/{'pca_harmony.py' if args.pipeline == 'pca_harmony' else 'gene_harmony_umap.py'} first."
        )

    tag = resolution_tag(args.resolution)
    pipeline_tag = PIPELINE_TAGS[args.pipeline]
    if args.output_dir == DEFAULT_OUTPUT_DIR:
        args.output_dir = CLUSTERING_ROOT / args.modality / f"{args.algorithm}_{pipeline_tag}_res{tag}"
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output = args.output_dir / f"simulation_{args.modality}_z_{args.algorithm}_{pipeline_tag}_res{tag}.h5ad"
    plot_dir = args.output_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    print(f"[load] {args.input}")
    adata = sc.read_h5ad(args.input)

    if args.rep_key not in adata.obsm:
        available = ", ".join(adata.obsm.keys())
        raise KeyError(
            f"Representation {args.rep_key!r} was not found in adata.obsm. "
            f"Available obsm keys: {available}"
        )

    rep = np.asarray(adata.obsm[args.rep_key])
    if args.pipeline == "pca_harmony":
        # PCA components are variance-ordered, so truncating to the top n_pcs is meaningful.
        n_pcs = min(args.n_pcs, rep.shape[1])
        active_rep_key = f"{args.rep_key}_{n_pcs}"
        adata.obsm[active_rep_key] = rep[:, :n_pcs].copy()
        print(
            f"[neighbors] Using first {n_pcs} dimensions of {args.rep_key} "
            f"with n_neighbors={args.n_neighbors}"
        )
    else:
        # Raw gene-space Harmony has no natural "first N" ordering; use all dimensions.
        n_pcs = rep.shape[1]
        active_rep_key = args.rep_key
        print(
            f"[neighbors] Using all {n_pcs} dimensions of {args.rep_key} "
            f"with n_neighbors={args.n_neighbors}"
        )
    sc.pp.neighbors(
        adata,
        n_neighbors=args.n_neighbors,
        use_rep=active_rep_key,
        random_state=args.random_state,
    )

    print(f"[{args.algorithm}] resolution={args.resolution}, key={args.cluster_key}")
    cluster_fn = sc.tl.leiden if args.algorithm == "leiden" else sc.tl.louvain
    cluster_fn(
        adata,
        resolution=args.resolution,
        key_added=args.cluster_key,
        random_state=args.random_state,
    )

    print(f"[umap] Computing UMAP from the {args.pipeline} neighbor graph")
    sc.tl.umap(adata, random_state=args.random_state)

    valid_colors = [key for key in args.plot_colors if key in adata.obs]
    missing = sorted(set(args.plot_colors) - set(valid_colors))
    if missing:
        print(f"[plot] Skipping missing obs columns: {', '.join(missing)}")

    for color_key in valid_colors:
        plot_umap(adata, color_key, plot_dir / f"umap_by_{color_key}.png")

    adata.uns[f"{args.algorithm}_umap"] = {
        "algorithm": args.algorithm,
        "pipeline": args.pipeline,
        "resolution": float(args.resolution),
        "cluster_key": args.cluster_key,
        "rep_key": args.rep_key,
        "n_pcs": int(n_pcs),
        "n_neighbors": int(args.n_neighbors),
    }
    adata.write_h5ad(output)
    print(f"[save] {output}")


if __name__ == "__main__":
    main()
