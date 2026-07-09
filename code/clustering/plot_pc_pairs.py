#!/usr/bin/env python
"""
Plot pairwise PCA component scatter plots before and after Harmony.

Example:
    python sim_paper/code/clustering/plot_pc_pairs.py
    python sim_paper/code/clustering/plot_pc_pairs.py --color-key cell_type_true
    python sim_paper/code/clustering/plot_pc_pairs.py --n-pcs 12 --no-sample
"""

from __future__ import annotations

import argparse
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


SCRIPT_DIR = Path(__file__).resolve().parent
SIM_PAPER_DIR = SCRIPT_DIR.parents[1]
DEFAULT_INPUT = SIM_PAPER_DIR / "data" / "clustering" / "simulation_spot_z_pca_harmony.h5ad"
DEFAULT_OUTPUT_DIR = SIM_PAPER_DIR / "data" / "clustering" / "plots" / "pc_pairs"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot all requested PC-vs-PC pairs before and after Harmony."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--pre-key", default="X_pca_pre_harmony")
    parser.add_argument("--post-key", default="X_pca_post_harmony")
    parser.add_argument("--color-key", default="slice_id")
    parser.add_argument(
        "--n-pcs",
        type=int,
        default=None,
        help="Number of leading PCs to plot. Defaults to every PC in the input.",
    )
    parser.add_argument(
        "--max-points",
        type=int,
        default=20000,
        help="Randomly sample this many observations per plot. Use --no-sample for all points.",
    )
    parser.add_argument("--no-sample", action="store_true")
    parser.add_argument("--point-size", type=float, default=3.0)
    parser.add_argument("--alpha", type=float, default=0.75)
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument("--dpi", type=int, default=160)
    return parser.parse_args()


def color_values(obs: pd.DataFrame, key: str):
    values = obs[key]
    if pd.api.types.is_numeric_dtype(values):
        return values.to_numpy(), None

    categories = pd.Categorical(values.astype(str))
    return categories.codes, list(categories.categories)


def sampled_indices(n_obs: int, max_points: int, no_sample: bool, random_state: int) -> np.ndarray:
    if no_sample or max_points <= 0 or max_points >= n_obs:
        return np.arange(n_obs)

    rng = np.random.default_rng(random_state)
    return np.sort(rng.choice(n_obs, size=max_points, replace=False))


def plot_pair(
    coords: np.ndarray,
    colors: np.ndarray,
    categories: list[str] | None,
    pc_x: int,
    pc_y: int,
    label: str,
    color_key: str,
    output_path: Path,
    point_size: float,
    alpha: float,
    dpi: int,
) -> None:
    fig, ax = plt.subplots(figsize=(5.5, 5))
    scatter = ax.scatter(
        coords[:, pc_x],
        coords[:, pc_y],
        c=colors,
        s=point_size,
        linewidths=0,
        cmap="tab20" if categories is not None else "viridis",
        alpha=alpha,
    )
    ax.set_xlabel(f"PC{pc_x + 1}")
    ax.set_ylabel(f"PC{pc_y + 1}")
    ax.set_title(f"{label}: PC{pc_x + 1} vs PC{pc_y + 1} by {color_key}")

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
                label=category,
            )
            for i, category in enumerate(categories)
        ]
        ax.legend(
            handles=handles,
            title=color_key,
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            frameon=False,
        )

    fig.tight_layout()
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_all_pairs(
    adata,
    embedding_key: str,
    label: str,
    output_dir: Path,
    color_key: str,
    n_pcs: int,
    obs_idx: np.ndarray,
    point_size: float,
    alpha: float,
    dpi: int,
) -> int:
    coords = np.asarray(adata.obsm[embedding_key])[:, :n_pcs][obs_idx]
    all_colors, categories = color_values(adata.obs, color_key)
    colors = all_colors[obs_idx]

    output_dir.mkdir(parents=True, exist_ok=True)
    n_plots = 0
    for pc_x, pc_y in combinations(range(n_pcs), 2):
        output_path = output_dir / f"pc{pc_x + 1:02d}_vs_pc{pc_y + 1:02d}_by_{color_key}.png"
        plot_pair(
            coords,
            colors,
            categories,
            pc_x,
            pc_y,
            label,
            color_key,
            output_path,
            point_size,
            alpha,
            dpi,
        )
        n_plots += 1

    return n_plots


def main() -> None:
    args = parse_args()

    if not args.input.is_file():
        raise SystemExit(
            f"Input file not found: {args.input}\n"
            "Run sim_paper/code/clustering/clustering.py first."
        )

    print(f"[load] {args.input}")
    adata = sc.read_h5ad(args.input)

    for key in [args.pre_key, args.post_key]:
        if key not in adata.obsm:
            available = ", ".join(adata.obsm.keys())
            raise KeyError(f"Embedding {key!r} not found in adata.obsm. Available: {available}")

    if args.color_key not in adata.obs:
        available = ", ".join(adata.obs.columns)
        raise KeyError(f"Color key {args.color_key!r} not found in adata.obs. Available: {available}")

    pre = np.asarray(adata.obsm[args.pre_key])
    post = np.asarray(adata.obsm[args.post_key])
    if pre.shape != post.shape:
        raise ValueError(f"Pre/post embeddings have different shapes: {pre.shape} vs {post.shape}")

    available_pcs = pre.shape[1]
    n_pcs = available_pcs if args.n_pcs is None else min(args.n_pcs, available_pcs)
    if n_pcs < 2:
        raise ValueError(f"Need at least 2 PCs to plot pairs; got n_pcs={n_pcs}.")

    obs_idx = sampled_indices(adata.n_obs, args.max_points, args.no_sample, args.random_state)
    sample_msg = "all observations" if len(obs_idx) == adata.n_obs else f"{len(obs_idx)} sampled observations"
    print(f"[plot] Plotting {n_pcs} PCs ({n_pcs * (n_pcs - 1) // 2} pairs) using {sample_msg}")

    color_dir = args.output_dir / f"by_{args.color_key}"
    pre_count = plot_all_pairs(
        adata,
        args.pre_key,
        "pre-Harmony",
        color_dir / "pre_harmony",
        args.color_key,
        n_pcs,
        obs_idx,
        args.point_size,
        args.alpha,
        args.dpi,
    )
    post_count = plot_all_pairs(
        adata,
        args.post_key,
        "post-Harmony",
        color_dir / "post_harmony",
        args.color_key,
        n_pcs,
        obs_idx,
        args.point_size,
        args.alpha,
        args.dpi,
    )

    print(f"[save] {pre_count} pre-Harmony plots -> {color_dir / 'pre_harmony'}")
    print(f"[save] {post_count} post-Harmony plots -> {color_dir / 'post_harmony'}")


if __name__ == "__main__":
    main()
