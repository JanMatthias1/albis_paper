#!/usr/bin/env python
"""
Run Harmony directly on the normalized all-gene feature matrix and compare UMAPs.

This intentionally skips the PCA projection used by clustering.py. The feature
matrix passed to Harmony is the scaled log-normalized expression matrix with one
dimension per gene.

Example:
    python sim_paper/code/clustering/gene_harmony_umap.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import harmonypy as hm
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


SCRIPT_DIR = Path(__file__).resolve().parent
SIM_PAPER_DIR = SCRIPT_DIR.parents[1]
DEFAULT_INPUT = SIM_PAPER_DIR / "data" / "simulation_spot_z.h5ad"
DEFAULT_OUTPUT_DIR = SIM_PAPER_DIR / "data" / "clustering"
DEFAULT_OUTPUT = DEFAULT_OUTPUT_DIR / "simulation_spot_z_gene_harmony_umap.h5ad"
VALID_MODALITIES = ("spot", "bin", "cell")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run all-gene Harmony and plot UMAPs before/after correction."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--modality", choices=VALID_MODALITIES, default="spot")
    parser.add_argument("--batch-key", default="slice_id")
    parser.add_argument("--target-sum", type=float, default=1e4)
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument("--point-size", type=float, default=4.0)
    parser.add_argument("--alpha", type=float, default=0.85)
    parser.add_argument(
        "--plot-colors",
        nargs="+",
        default=["slice_id", "domain_true", "cell_type_true"],
        help="obs columns to use for before/after UMAP plots.",
    )
    return parser.parse_args()


def dense_float32(matrix) -> np.ndarray:
    if sparse.issparse(matrix):
        matrix = matrix.toarray()
    return np.asarray(matrix, dtype=np.float32)


def preprocess_all_genes(adata, target_sum: float) -> np.ndarray:
    adata.var_names_make_unique()
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, max_value=10)
    return dense_float32(adata.X)


def run_harmony_matrix(features: np.ndarray, obs: pd.DataFrame, batch_key: str) -> np.ndarray:
    harmony_out = hm.run_harmony(features.astype(np.float64, copy=False), obs, batch_key)
    corrected = np.asarray(harmony_out.Z_corr)

    if corrected.shape == features.shape:
        return corrected.astype(np.float32, copy=False)

    if corrected.T.shape == features.shape:
        return corrected.T.astype(np.float32, copy=False)

    raise ValueError(
        f"Harmony returned shape {corrected.shape}; expected {features.shape} "
        "or its transpose."
    )


def compute_umap(adata, rep_key: str, umap_key: str, n_neighbors: int, random_state: int) -> None:
    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        use_rep=rep_key,
        random_state=random_state,
    )
    sc.tl.umap(adata, random_state=random_state)
    adata.obsm[umap_key] = np.asarray(adata.obsm["X_umap"]).copy()


def color_values(obs: pd.DataFrame, key: str):
    values = obs[key]
    if pd.api.types.is_numeric_dtype(values):
        return values.to_numpy(), None

    categories = pd.Categorical(values.astype(str))
    return categories.codes, list(categories.categories)


def plot_umap_before_after(
    adata,
    color_key: str,
    output_path: Path,
    point_size: float,
    alpha: float,
) -> None:
    colors, categories = color_values(adata.obs, color_key)
    pre = np.asarray(adata.obsm["X_umap_gene_pre_harmony"])
    post = np.asarray(adata.obsm["X_umap_gene_harmony"])

    fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharex=False, sharey=False)
    scatters = []
    for ax, coords, title in [
        (axes[0], pre, "Before Harmony"),
        (axes[1], post, "After Harmony"),
    ]:
        scatter = ax.scatter(
            coords[:, 0],
            coords[:, 1],
            c=colors,
            s=point_size,
            linewidths=0,
            cmap="tab20" if categories is not None else "viridis",
            alpha=alpha,
        )
        scatters.append(scatter)
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_title(title)

    if categories is None:
        fig.colorbar(scatters[-1], ax=axes, fraction=0.046, pad=0.04, label=color_key)
    else:
        handles = [
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor=scatters[-1].cmap(scatters[-1].norm(i)),
                markersize=6,
                label=label,
            )
            for i, label in enumerate(categories)
        ]
        axes[-1].legend(
            handles=handles,
            title=color_key,
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            frameon=False,
        )

    fig.suptitle(f"All-gene UMAP before vs after Harmony by {color_key}")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def save_umap_plots(adata, plot_dir: Path, color_keys: list[str], point_size: float, alpha: float) -> None:
    plot_dir.mkdir(parents=True, exist_ok=True)
    valid_colors = [key for key in color_keys if key in adata.obs]
    missing = sorted(set(color_keys) - set(valid_colors))
    if missing:
        print(f"[plot] Skipping missing obs columns: {', '.join(missing)}")

    for color_key in valid_colors:
        plot_umap_before_after(
            adata,
            color_key,
            plot_dir / f"umap_gene_harmony_before_after_by_{color_key}.png",
            point_size,
            alpha,
        )


def main() -> None:
    args = parse_args()
    if args.input == DEFAULT_INPUT:
        args.input = SIM_PAPER_DIR / "data" / f"simulation_{args.modality}_z.h5ad"
    if args.output == DEFAULT_OUTPUT:
        args.output = DEFAULT_OUTPUT_DIR / f"simulation_{args.modality}_z_gene_harmony_umap.h5ad"

    if not args.input.is_file():
        raise SystemExit(
            f"Input file not found: {args.input}\n"
            "Generate it first with sim_paper/code/data/generate_simulation.py."
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    plot_dir = args.output.parent / "plots" / f"{args.modality}_gene_harmony_umap"

    print(f"[load] {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"[load] AnnData shape: {adata.n_obs} observations x {adata.n_vars} genes")

    if args.batch_key not in adata.obs:
        raise KeyError(f"Batch key {args.batch_key!r} was not found in adata.obs.")

    print("[features] Normalizing, log-transforming, and scaling all genes")
    features = preprocess_all_genes(adata, args.target_sum)
    adata.obsm["X_gene_pre_harmony"] = features

    print(f"[harmony] Correcting X_gene_pre_harmony by obs['{args.batch_key}']")
    adata.obsm["X_gene_harmony"] = run_harmony_matrix(features, adata.obs, args.batch_key)

    print("[umap] Computing UMAP before Harmony from all-gene features")
    compute_umap(
        adata,
        "X_gene_pre_harmony",
        "X_umap_gene_pre_harmony",
        args.n_neighbors,
        args.random_state,
    )

    print("[umap] Computing UMAP after Harmony from all-gene features")
    compute_umap(
        adata,
        "X_gene_harmony",
        "X_umap_gene_harmony",
        args.n_neighbors,
        args.random_state,
    )

    print(f"[plot] Saving before/after UMAP plots under {plot_dir}")
    save_umap_plots(adata, plot_dir, args.plot_colors, args.point_size, args.alpha)

    adata.uns["gene_harmony_umap"] = {
        "batch_key": args.batch_key,
        "target_sum": float(args.target_sum),
        "n_neighbors": int(args.n_neighbors),
        "used_all_genes_without_pca": True,
        "pre_harmony_obsm": "X_gene_pre_harmony",
        "post_harmony_obsm": "X_gene_harmony",
        "pre_harmony_umap_obsm": "X_umap_gene_pre_harmony",
        "post_harmony_umap_obsm": "X_umap_gene_harmony",
        "input": str(args.input),
    }
    adata.write_h5ad(args.output)
    print(f"[save] {args.output}")


if __name__ == "__main__":
    main()
