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
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import scanpy as sc


SCRIPT_DIR = Path(__file__).resolve().parent
SIM_PAPER_DIR = SCRIPT_DIR.parents[1]
DEFAULT_INPUT = SIM_PAPER_DIR / "data" / "clustering" / "simulation_spot_z_pca_harmony.h5ad"
DEFAULT_OUTPUT_DIR = SIM_PAPER_DIR / "data" / "clustering" / "plots" / "pc_pairs"
VALID_MODALITIES = ("spot", "bin", "cell")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot adjacent PC pairs before and after Harmony as two grid figures."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--modality", choices=VALID_MODALITIES, default="spot")
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
    parser.add_argument("--ncols", type=int, default=5)
    parser.add_argument("--legend-max-levels", type=int, default=30)
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument("--dpi", type=int, default=160)
    return parser.parse_args()


def make_color_lookup(labels: pd.Series) -> dict[str, tuple[float, float, float]]:
    label_list = sorted(pd.unique(labels.astype(str)))
    n_labels = len(label_list)
    if n_labels <= 10:
        colors = list(plt.get_cmap("tab10").colors)[:n_labels]
    elif n_labels <= 20:
        colors = list(plt.get_cmap("tab20").colors)[:n_labels]
    elif n_labels <= 60:
        colors = (
            list(plt.get_cmap("tab20").colors)
            + list(plt.get_cmap("tab20b").colors)
            + list(plt.get_cmap("tab20c").colors)
        )[:n_labels]
    else:
        colors = [mcolors.hsv_to_rgb([i / n_labels, 0.65, 0.95]) for i in range(n_labels)]

    return dict(zip(label_list, colors))


def sampled_indices(n_obs: int, max_points: int, no_sample: bool, random_state: int) -> np.ndarray:
    if no_sample or max_points <= 0 or max_points >= n_obs:
        return np.arange(n_obs)

    rng = np.random.default_rng(random_state)
    return np.sort(rng.choice(n_obs, size=max_points, replace=False))


def pc_variance_fractions(coords: np.ndarray) -> np.ndarray:
    variances = coords.var(axis=0)
    total = variances.sum()
    if total == 0:
        return np.zeros_like(variances)
    return variances / total


def plot_pc_pair_grid(
    coords: np.ndarray,
    color_labels: pd.Series,
    label: str,
    color_key: str,
    n_pcs: int,
    ncols: int,
    output_path: Path,
    point_size: float,
    alpha: float,
    dpi: int,
    legend_max_levels: int,
) -> None:
    pairs = [(i, i + 1) for i in range(0, n_pcs - 1, 2)]
    n_panels = len(pairs)
    nrows = math.ceil(n_panels / ncols)
    var_frac = pc_variance_fractions(coords[:, :n_pcs])

    labels = color_labels.astype(str)
    lookup = make_color_lookup(labels)
    point_colors = labels.map(lookup).to_numpy()

    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3.2, nrows * 3.2))
    axes = np.atleast_1d(axes).ravel()

    for ax, (pc_x, pc_y) in zip(axes, pairs):
        ax.scatter(
            coords[:, pc_x],
            coords[:, pc_y],
            c=point_colors,
            s=point_size,
            linewidths=0,
            alpha=alpha,
            rasterized=True,
        )
        ax.set_xlabel(f"PC{pc_x + 1} ({var_frac[pc_x] * 100:.1f}%)", fontsize=7)
        ax.set_ylabel(f"PC{pc_y + 1} ({var_frac[pc_y] * 100:.1f}%)", fontsize=7)
        ax.tick_params(labelsize=6)

    for ax in axes[n_panels:]:
        ax.axis("off")

    legend_note = ""
    if len(lookup) <= legend_max_levels:
        handles = [mpatches.Patch(color=color, label=category) for category, color in lookup.items()]
        fig.legend(
            handles=handles,
            fontsize=6,
            loc="center left",
            bbox_to_anchor=(1.0, 0.5),
            frameon=False,
            ncol=max(1, math.ceil(len(lookup) / 25)),
        )
    else:
        legend_note = f" | legend suppressed: {len(lookup)} levels"

    if n_pcs % 2:
        legend_note += f" | PC{n_pcs} unpaired"

    fig.suptitle(
        f"{label} PC pairs by {color_key} | N={coords.shape[0]:,}{legend_note}",
        fontsize=11,
    )
    fig.tight_layout(rect=[0, 0, 0.92, 0.97])
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_embedding(
    adata,
    embedding_key: str,
    label: str,
    output_path: Path,
    color_key: str,
    n_pcs: int,
    obs_idx: np.ndarray,
    point_size: float,
    alpha: float,
    dpi: int,
    ncols: int,
    legend_max_levels: int,
) -> None:
    coords = np.asarray(adata.obsm[embedding_key])[:, :n_pcs][obs_idx]
    color_labels = adata.obs.iloc[obs_idx][color_key]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plot_pc_pair_grid(
        coords,
        color_labels,
        label,
        color_key,
        n_pcs,
        ncols,
        output_path,
        point_size,
        alpha,
        dpi,
        legend_max_levels,
    )


def main() -> None:
    args = parse_args()
    if args.input == DEFAULT_INPUT:
        args.input = SIM_PAPER_DIR / "data" / "clustering" / f"simulation_{args.modality}_z_pca_harmony.h5ad"
    if args.output_dir == DEFAULT_OUTPUT_DIR:
        args.output_dir = DEFAULT_OUTPUT_DIR / args.modality

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
    n_pairs = n_pcs // 2
    print(f"[plot] Plotting {n_pcs} PCs ({n_pairs} adjacent pairs) using {sample_msg}")

    color_dir = args.output_dir / f"by_{args.color_key}"
    pre_output = color_dir / f"pc_pairs_pre_harmony_by_{args.color_key}.png"
    post_output = color_dir / f"pc_pairs_post_harmony_by_{args.color_key}.png"

    plot_embedding(
        adata,
        args.pre_key,
        "pre-Harmony",
        pre_output,
        args.color_key,
        n_pcs,
        obs_idx,
        args.point_size,
        args.alpha,
        args.dpi,
        args.ncols,
        args.legend_max_levels,
    )
    plot_embedding(
        adata,
        args.post_key,
        "post-Harmony",
        post_output,
        args.color_key,
        n_pcs,
        obs_idx,
        args.point_size,
        args.alpha,
        args.dpi,
        args.ncols,
        args.legend_max_levels,
    )

    print(f"[save] pre-Harmony plot -> {pre_output}")
    print(f"[save] post-Harmony plot -> {post_output}")


if __name__ == "__main__":
    main()
