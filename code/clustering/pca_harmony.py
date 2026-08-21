#!/usr/bin/env python
"""
Run PCA and Harmony correction for the simulated sim_app dataset.

This script uses all genes by default because the simulated data has only
about 550 genes. Harmony is run over ``obs['slice_id']``, which is the same
field used by sim_app when adding slice-specific batch effects.

Expected environment:
    conda activate /dcs04/hicks/data/Jan/sim_project/sim_app/env/sim-app-tutorial
    python -m pip install -r sim_paper/code/clustering/requirements.txt

Example:
    python sim_paper/code/clustering/pca_harmony.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import harmonypy as hm
import numpy as np
import pandas as pd
import scanpy as sc

from pc_pairs import plot_pc_pairs, sampled_indices


SCRIPT_DIR = Path(__file__).resolve().parent
SIM_PAPER_DIR = SCRIPT_DIR.parents[1]
DEFAULT_INPUT = SIM_PAPER_DIR / "data" / "simulation_spot_z.h5ad"
CLUSTERING_ROOT = SIM_PAPER_DIR / "data" / "clustering"
DEFAULT_OUTPUT = CLUSTERING_ROOT / "spot" / "pca_harmony" / "simulation_spot_z_pca_harmony.h5ad"
VALID_MODALITIES = ("spot", "bin", "cell")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run PCA, Harmony correction, and PCA diagnostic plots."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--modality", choices=VALID_MODALITIES, default="spot")
    parser.add_argument("--batch-key", default="slice_id")
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--target-sum", type=float, default=1e4)
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--point-size", type=float, default=4.0)
    parser.add_argument("--alpha", type=float, default=0.85)
    parser.add_argument(
        "--umap-max-obs",
        type=int,
        default=50000,
        help="Subsample to at most this many observations before computing the before/after-Harmony "
        "UMAPs (neighbors+UMAP twice at full scale is the slow step for million-cell datasets; "
        "downstream clustering reads X_pca_harmony from the full, non-subsampled data, not this UMAP, "
        "so subsampling only affects the diagnostic plots). Use --no-umap-sample to disable.",
    )
    parser.add_argument(
        "--no-umap-sample",
        action="store_true",
        help="Compute the before/after UMAPs on all observations instead of subsampling (slow for large datasets).",
    )
    parser.add_argument(
        "--plot-colors",
        nargs="+",
        default=["slice_id", "domain_true", "cell_type_true"],
        help="obs columns to use for PCA before/after Harmony plots.",
    )
    return parser.parse_args()


def preprocess_for_pca(adata, n_pcs: int, random_state: int, target_sum: float) -> None:
    adata.var_names_make_unique()
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(
        adata,
        n_comps=n_pcs,
        use_highly_variable=False,
        svd_solver="arpack",
        random_state=random_state,
    )
    adata.obsm["X_pca_pre_harmony"] = adata.obsm["X_pca"].copy()


def run_harmony(adata, batch_key: str, basis: str, adjusted_basis: str) -> None:
    pca = np.asarray(adata.obsm[basis], dtype=np.float64)
    harmony_out = hm.run_harmony(pca, adata.obs, batch_key)
    corrected = np.asarray(harmony_out.Z_corr)

    if corrected.shape == pca.shape:
        adata.obsm[adjusted_basis] = corrected
        return

    if corrected.T.shape == pca.shape:
        adata.obsm[adjusted_basis] = corrected.T
        return

    raise ValueError(
        f"Harmony returned shape {corrected.shape}; expected {pca.shape} "
        f"or its transpose for adata.obsm[{adjusted_basis!r}]."
    )


def color_values(obs: pd.DataFrame, key: str):
    values = obs[key]
    if pd.api.types.is_numeric_dtype(values):
        return values.to_numpy(), None

    categories = pd.Categorical(values.astype(str))
    return categories.codes, list(categories.categories)


def plot_two_dims(adata, embedding_key: str, color_key: str, output_path: Path) -> None:
    coords = np.asarray(adata.obsm[embedding_key])
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
    ax.set_xlabel(f"{embedding_key} 1")
    ax.set_ylabel(f"{embedding_key} 2")
    ax.set_title(f"{embedding_key}: {color_key}")

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


def save_pca_plots(adata, plot_dir: Path, color_keys: list[str]) -> None:
    plot_dir.mkdir(parents=True, exist_ok=True)
    valid_colors = [key for key in color_keys if key in adata.obs]
    missing = sorted(set(color_keys) - set(valid_colors))
    if missing:
        print(f"[plot] Skipping missing obs columns: {', '.join(missing)}")

    for color_key in valid_colors:
        plot_two_dims(
            adata,
            "X_pca_pre_harmony",
            color_key,
            plot_dir / f"pca_before_harmony_by_{color_key}.png",
        )
        plot_two_dims(
            adata,
            "X_pca_post_harmony",
            color_key,
            plot_dir / f"pca_after_harmony_by_{color_key}.png",
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


def plot_umap_before_after(
    adata,
    color_key: str,
    output_path: Path,
    point_size: float,
    alpha: float,
) -> None:
    colors, categories = color_values(adata.obs, color_key)
    pre = np.asarray(adata.obsm["X_umap_pca_pre_harmony"])
    post = np.asarray(adata.obsm["X_umap_pca_post_harmony"])

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

    fig.suptitle(f"PCA UMAP before vs after Harmony by {color_key}")
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
            plot_dir / f"umap_pca_harmony_before_after_by_{color_key}.png",
            point_size,
            alpha,
        )


def main() -> None:
    args = parse_args()
    if args.input == DEFAULT_INPUT:
        args.input = SIM_PAPER_DIR / "data" / f"simulation_{args.modality}_z.h5ad"
    if args.output == DEFAULT_OUTPUT:
        args.output = (
            CLUSTERING_ROOT
            / args.modality
            / "pca_harmony"
            / f"simulation_{args.modality}_z_pca_harmony.h5ad"
        )

    if not args.input.is_file():
        raise SystemExit(
            f"Input file not found: {args.input}\n"
            "Generate it first with sim_paper/code/data/generate_simulation.py."
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    plot_dir = args.output.parent / "plots"

    print(f"[load] {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"[load] AnnData shape: {adata.n_obs} observations x {adata.n_vars} genes")

    if args.batch_key not in adata.obs:
        raise KeyError(f"Batch key {args.batch_key!r} was not found in adata.obs.")

    print(f"[pca] Running PCA with {args.n_pcs} components over all genes")
    preprocess_for_pca(adata, args.n_pcs, args.random_state, args.target_sum)

    print(f"[harmony] Correcting X_pca_pre_harmony by obs['{args.batch_key}']")
    run_harmony(adata, args.batch_key, "X_pca_pre_harmony", "X_pca_post_harmony")
    adata.obsm["X_pca"] = adata.obsm["X_pca_pre_harmony"].copy()
    adata.obsm["X_pca_harmony"] = adata.obsm["X_pca_post_harmony"].copy()

    print(f"[plot] Saving PCA before/after plots under {plot_dir}")
    save_pca_plots(adata, plot_dir, args.plot_colors)

    obs_idx = sampled_indices(adata.n_obs, args.umap_max_obs, args.no_umap_sample, args.random_state)
    if len(obs_idx) < adata.n_obs:
        print(
            f"[umap] Subsampling {len(obs_idx):,} / {adata.n_obs:,} observations for the before/after-Harmony "
            "UMAP diagnostic plots (full neighbors+UMAP twice is the slow step; not used downstream)"
        )
    umap_adata = adata[obs_idx].copy() if len(obs_idx) < adata.n_obs else adata

    print("[umap] Computing UMAP before Harmony from PCA")
    compute_umap(
        umap_adata,
        "X_pca_pre_harmony",
        "X_umap_pca_pre_harmony",
        args.n_neighbors,
        args.random_state,
    )

    print("[umap] Computing UMAP after Harmony from PCA")
    compute_umap(
        umap_adata,
        "X_pca_post_harmony",
        "X_umap_pca_post_harmony",
        args.n_neighbors,
        args.random_state,
    )

    print(f"[plot] Saving before/after UMAP plots under {plot_dir}")
    save_umap_plots(umap_adata, plot_dir, args.plot_colors, args.point_size, args.alpha)

    print(f"[pc_pairs] Plotting PC pairs by {args.batch_key} under {plot_dir / 'pc_pairs'}")
    plot_pc_pairs(
        adata,
        pre_key="X_pca_pre_harmony",
        post_key="X_pca_post_harmony",
        output_dir=plot_dir / "pc_pairs",
        color_key=args.batch_key,
    )

    adata.uns["clustering_pca_harmony"] = {
        "batch_key": args.batch_key,
        "n_pcs": int(args.n_pcs),
        "used_all_genes": True,
        "pre_harmony_obsm": "X_pca_pre_harmony",
        "post_harmony_obsm": "X_pca_post_harmony",
        "n_neighbors": int(args.n_neighbors),
        "umap_obsm_note": "before/after UMAP is diagnostic-plot-only, computed on a subsample "
        "(see umap_n_obs) and not stored in this file; downstream clustering uses pre_harmony_obsm/"
        "post_harmony_obsm from the full data.",
        "umap_n_obs": int(umap_adata.n_obs),
        "input": str(args.input),
    }
    adata.write_h5ad(args.output)
    print(f"[save] {args.output}")


if __name__ == "__main__":
    main()
