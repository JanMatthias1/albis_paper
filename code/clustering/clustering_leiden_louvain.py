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

# Shared with pca_harmony.py's plot_umap_before_after so the two side-by-side UMAP
# comparison panels in Figure 3 come out at the same aspect ratio (see that
# function for why bbox_inches="tight" isn't used here).
PANEL_FIGSIZE = (12, 6)
PANEL_DPI = 180
PANEL_RECT = (0.0, 0.10, 1.0, 0.93)
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
    parser.add_argument(
        "--packing-tag", default=None,
        help="If set, use data/clustering_<packing-tag>/ as the root instead of data/clustering/ "
        "(matches ari_vs_ground_truth.py's --packing-tag). Ignored if --input/--output-dir given.",
    )
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


def plot_umap_true_vs_predicted(adata, true_key: str, pred_key: str, output_path: Path) -> None:
    """Side-by-side UMAP: ground truth on the left, predicted clusters on the right,
    same coordinates/axis limits so shapes are directly comparable without flipping
    between two separate images."""
    coords = np.asarray(adata.obsm["X_umap"])
    xlim = (coords[:, 0].min() - 1, coords[:, 0].max() + 1)
    ylim = (coords[:, 1].min() - 1, coords[:, 1].max() + 1)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
    for ax, key, title in ((axes[0], true_key, f"Ground truth: {true_key}"), (axes[1], pred_key, f"Predicted: {pred_key}")):
        colors, categories = color_values(adata.obs, key)
        scatter = ax.scatter(coords[:, 0], coords[:, 1], c=colors, s=5, linewidths=0, cmap="tab20", alpha=0.85)
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_title(title)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        if categories is not None:
            handles = [
                plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=scatter.cmap(scatter.norm(i)), markersize=6, label=label)
                for i, label in enumerate(categories)
            ]
            ax.legend(handles=handles, title=key, bbox_to_anchor=(0.5, -0.15), loc="upper center", ncol=min(len(categories), 8), frameon=False, fontsize=7)

    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_contingency_heatmap(adata, true_key: str, pred_key: str, output_path: Path) -> None:
    """Row-normalized contingency table (true label x predicted cluster): for each
    true category, what fraction of its cells landed in each predicted cluster.
    A clean recovery shows up as (close to) one bright cell per row; a smeared row
    means that true category isn't separable from the others in this embedding."""
    true_labels = adata.obs[true_key].astype(str)
    pred_labels = adata.obs[pred_key].astype(str)
    table = pd.crosstab(true_labels, pred_labels)
    table = table.reindex(sorted(table.index, key=str), axis=0)
    table = table.reindex(sorted(table.columns, key=lambda x: int(x) if x.lstrip("-").isdigit() else x), axis=1)
    frac = table.div(table.sum(axis=1), axis=0)

    fig, ax = plt.subplots(figsize=(max(5, 0.6 * frac.shape[1] + 2), max(4, 0.5 * frac.shape[0] + 1.5)))
    im = ax.imshow(frac.to_numpy(), cmap="viridis", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(frac.shape[1]))
    ax.set_xticklabels(frac.columns, rotation=90)
    ax.set_yticks(range(frac.shape[0]))
    ax.set_yticklabels(frac.index)
    ax.set_xlabel(f"Predicted cluster ({pred_key})")
    ax.set_ylabel(f"Ground truth ({true_key})")
    ax.set_title(f"Fraction of each {true_key} category per predicted cluster")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="fraction of row")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    clustering_root = (
        SIM_PAPER_DIR / "data" / f"clustering_{args.packing_tag}" if args.packing_tag else CLUSTERING_ROOT
    )
    pipeline_dirname = "pca_harmony" if args.pipeline == "pca_harmony" else "gene_harmony_umap"
    if args.input == DEFAULT_INPUT:
        args.input = (
            clustering_root / args.modality / pipeline_dirname / f"simulation_{args.modality}_z_{pipeline_dirname}.h5ad"
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
        args.output_dir = clustering_root / args.modality / f"{args.algorithm}_{pipeline_tag}_res{tag}"
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
    extra_kwargs = (
        {"flavor": "igraph", "n_iterations": 2, "directed": False}
        if args.algorithm == "leiden" else {}
    )
    cluster_fn(
        adata,
        resolution=args.resolution,
        key_added=args.cluster_key,
        random_state=args.random_state,
        **extra_kwargs,
    )

    print(f"[umap] Computing UMAP from the {args.pipeline} neighbor graph")
    sc.tl.umap(adata, random_state=args.random_state)

    valid_colors = [key for key in args.plot_colors if key in adata.obs]
    missing = sorted(set(args.plot_colors) - set(valid_colors))
    if missing:
        print(f"[plot] Skipping missing obs columns: {', '.join(missing)}")

    for color_key in valid_colors:
        plot_umap(adata, color_key, plot_dir / f"umap_by_{color_key}.png")

    # Ground-truth-vs-predicted comparisons: side-by-side UMAP + contingency heatmap,
    # for each true-label column that's actually present alongside cluster_label.
    for true_key in ("cell_type_true", "domain_true"):
        if true_key not in adata.obs or args.cluster_key not in adata.obs:
            continue
        plot_umap_true_vs_predicted(adata, true_key, args.cluster_key, plot_dir / f"umap_true_vs_predicted_{true_key}.png")
        plot_contingency_heatmap(adata, true_key, args.cluster_key, plot_dir / f"contingency_{true_key}.png")

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
