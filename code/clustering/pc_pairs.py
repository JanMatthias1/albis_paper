"""Shared PC-pair grid plotting, used by step01_pca_harmony.py.

Plots adjacent PC pairs (PC1 vs PC2, PC3 vs PC4, ...) before and after Harmony
as two grid figures, colored by a single obs column.
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd


def make_color_lookup(labels: pd.Series):
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from manuscript_style import category_color, category_order
    key = labels.name or "slice_id"
    return {label: category_color(label, key) for label in category_order(labels)}


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
            fontsize=12,
            title="Slice ID" if color_key == "slice_id" else color_key,
            title_fontsize=13,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.0),
            frameon=False,
            ncol=min(len(lookup), 10),
        )
    else:
        legend_note = f" | legend suppressed: {len(lookup)} levels"

    if n_pcs % 2:
        legend_note += f" | PC{n_pcs} unpaired"

    fig.suptitle(
        f"{label} PC pairs by {color_key} | N={coords.shape[0]:,}{legend_note}",
        fontsize=17,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.97])
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


def plot_pc_pairs(
    adata,
    pre_key: str,
    post_key: str,
    output_dir: Path,
    color_key: str = "slice_id",
    n_pcs: int | None = None,
    max_points: int = 20000,
    no_sample: bool = False,
    point_size: float = 3.0,
    alpha: float = 0.75,
    ncols: int = 5,
    legend_max_levels: int = 30,
    random_state: int = 0,
    dpi: int = 160,
) -> None:
    """Plot pre/post-Harmony PC-pair grids for adata.obsm[pre_key]/[post_key].

    Silently no-ops (with a printed note) if either key is missing, the color
    key is missing, or fewer than 2 PCs are available, since this is a
    diagnostic plot that should not fail the calling pipeline.
    """
    for key in (pre_key, post_key):
        if key not in adata.obsm:
            print(f"[pc_pairs] Skipping: {key!r} not found in adata.obsm")
            return

    if color_key not in adata.obs:
        print(f"[pc_pairs] Skipping: color key {color_key!r} not found in adata.obs")
        return

    pre = np.asarray(adata.obsm[pre_key])
    post = np.asarray(adata.obsm[post_key])
    if pre.shape != post.shape:
        print(f"[pc_pairs] Skipping: {pre_key} shape {pre.shape} != {post_key} shape {post.shape}")
        return

    available_pcs = pre.shape[1]
    n_pcs_eff = available_pcs if n_pcs is None else min(n_pcs, available_pcs)
    if n_pcs_eff < 2:
        print(f"[pc_pairs] Skipping: need at least 2 PCs, got {n_pcs_eff}")
        return

    obs_idx = sampled_indices(adata.n_obs, max_points, no_sample, random_state)
    sample_msg = "all observations" if len(obs_idx) == adata.n_obs else f"{len(obs_idx)} sampled observations"
    n_pairs = n_pcs_eff // 2
    print(f"[pc_pairs] Plotting {n_pcs_eff} PCs ({n_pairs} adjacent pairs) using {sample_msg}")

    color_dir = output_dir / f"by_{color_key}"
    pre_output = color_dir / f"pc_pairs_pre_harmony_by_{color_key}.png"
    post_output = color_dir / f"pc_pairs_post_harmony_by_{color_key}.png"

    plot_embedding(
        adata, pre_key, "pre-Harmony", pre_output, color_key, n_pcs_eff, obs_idx,
        point_size, alpha, dpi, ncols, legend_max_levels,
    )
    plot_embedding(
        adata, post_key, "post-Harmony", post_output, color_key, n_pcs_eff, obs_idx,
        point_size, alpha, dpi, ncols, legend_max_levels,
    )

    print(f"[pc_pairs] Saved -> {pre_output}")
    print(f"[pc_pairs] Saved -> {post_output}")
