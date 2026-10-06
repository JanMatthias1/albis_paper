#!/usr/bin/env python
"""
Run Leiden or Louvain clustering and UMAP from a Harmony-corrected representation.

This is separate from step01_pca_harmony.py so clustering algorithm/resolution can be
iterated without recomputing normalization, PCA, or Harmony. Clusters on
step01_pca_harmony.py's output (Harmony-corrected PCA embedding), using the first
--n-pcs dimensions (variance-ordered).

--pipeline is kept as a single-choice flag (every caller passes
"pca_harmony" explicitly) so another embedding can be added later.

Example:
    python sim_paper/code/clustering/step03_cluster_and_plot.py --resolution 0.5
    python sim_paper/code/clustering/step03_cluster_and_plot.py --algorithm louvain --resolution 1.0
"""

from __future__ import annotations

import argparse
import sys
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
import numpy as np
import pandas as pd
import scanpy as sc


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from manuscript_style import (
    save_figure, apply_style, scatter_colors, legend_handles, pretty_label,
    LEGEND_MARKERSIZE, LEGEND_TITLE_SIZE, PANEL_FIGSIZE, PANEL_MARGINS, PANEL_EXPORT_BOTTOM,
    matched_labels, category_order,
)
apply_style()

SCRIPT_DIR = Path(__file__).resolve().parent
SIM_PAPER_DIR = SCRIPT_DIR.parents[1]
CLUSTERING_ROOT = SIM_PAPER_DIR / "data" / "clustering"
DEFAULT_INPUT = CLUSTERING_ROOT / "spot" / "pca_harmony" / "simulation_spot_z_pca_harmony.h5ad"
DEFAULT_OUTPUT_DIR = CLUSTERING_ROOT

PANEL_DPI = 500
VALID_MODALITIES = ("spot", "bin", "cell")
VALID_ALGORITHMS = ("leiden", "louvain")
VALID_PIPELINES = ("pca_harmony",)
PIPELINE_TAG = "pca"
PIPELINE_REP_KEY = "X_pca_harmony"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run neighbors, Leiden clustering, UMAP, and UMAP plots."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--packing-tag", default=None,
        help="If set, use data/clustering_<packing-tag>/ as the root instead of data/clustering/ "
        "(matches step02_leiden_resolution_sweep.py's --packing-tag). Ignored if --input/--output-dir given.",
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
        help="obsm key to cluster on. Defaults to X_pca_harmony.",
    )
    parser.add_argument("--cluster-key", default="cluster_label")
    parser.add_argument(
        "--plot-colors",
        nargs="+",
        default=["domain_true", "cell_type_true", "cluster_label"],
        help="obs columns to use for UMAP plots.",
    )
    parser.add_argument(
        "--plots-only",
        action="store_true",
        help="Skip neighbors/clustering/UMAP and just regenerate plots from the existing "
        "output h5ad for this --modality/--algorithm/--pipeline/--resolution (must already "
        "have --cluster-key and X_umap) -- for re-plotting after a plot-formatting change "
        "without rerunning clustering.",
    )
    return parser.parse_args()


def resolution_tag(resolution: float) -> str:
    return str(resolution).replace(".", "p").replace("-", "m")


def fixed_margin_legend_ncol(categories, max_wide: int = 10) -> int:
    """Column count for plot_umap_true_vs_predicted's bottom-anchored legend,
    which sits inside PANEL_MARGINS' fixed strip with no tight_layout/
    bbox_inches to rescue a horizontal overflow. Short labels (cluster_label
    "0".."7") fit a full wide row even at 8 categories -- verified by
    rendering, same fix as step01_pca_harmony.py's identically-
    named helper. Longer labels (cell_type_true "type1".."type8") overflow
    past the axis edge at the same count, so those still wrap to a narrower
    4-column row."""
    max_label_len = max(len(str(c)) for c in categories)
    return min(len(categories), max_wide if max_label_len <= 2 else 4)


def color_values(obs: pd.DataFrame, key: str):
    values = obs[key]
    if pd.api.types.is_numeric_dtype(values):
        return values.to_numpy(), None

    categories = pd.Categorical(values.astype(str))
    return categories.codes, list(categories.categories)


def plot_umap(adata, color_key: str, output_path: Path) -> None:
    coords = np.asarray(adata.obsm["X_umap"])
    color_args, categories, swatches = scatter_colors(adata.obs, color_key,
        "cell_type_true" if color_key == "cluster_label" and "cell_type_true" in adata.obs else None)

    fig, ax = plt.subplots(figsize=(6, 5))
    scatter = ax.scatter(
        coords[:, 0],
        coords[:, 1],
        **color_args,
        s=5,
        linewidths=0,
        alpha=0.85,
    )
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title(f"UMAP: {pretty_label(color_key)}")

    if categories is None:
        fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04, label=pretty_label(color_key))
    else:
        handles = legend_handles(categories, swatches)
        ax.legend(
            handles=handles,
            title=pretty_label(color_key),
            title_fontsize=LEGEND_TITLE_SIZE,
            bbox_to_anchor=(0.5, -0.18),
            loc="upper center",
            ncol=2 if max(map(len, categories)) > 5 else min(len(categories), 4),
            frameon=False,
        )

    fig.tight_layout()
    save_figure(fig, output_path, dpi=500, bbox_inches="tight")
    plt.close(fig)


def plot_umap_true_vs_predicted(adata, true_key: str, pred_key: str, output_path: Path) -> None:
    """Side-by-side UMAP: ground truth on the left, predicted clusters on the right,
    same coordinates/axis limits so shapes are directly comparable without flipping
    between two separate images."""
    coords = np.asarray(adata.obsm["X_umap"])
    xlim = (coords[:, 0].min() - 1, coords[:, 0].max() + 1)
    ylim = (coords[:, 1].min() - 1, coords[:, 1].max() + 1)

    fig, axes = plt.subplots(1, 2, figsize=PANEL_FIGSIZE)
    legends = []
    for ax, key, title in (
        (axes[0], true_key, f"Ground truth: {pretty_label(true_key)}"),
        (axes[1], pred_key, f"Predicted: {pretty_label(pred_key)}"),
    ):
        color_args, categories, swatches = scatter_colors(
            adata.obs, key, true_key if key == pred_key else None)
        ax.scatter(coords[:, 0], coords[:, 1], **color_args, s=5, linewidths=0, alpha=0.85)
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_title(title)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        if categories is not None:
            if key == pred_key:
                categories = category_order(adata.obs[key].astype(str))
            handles = legend_handles(categories, swatches)
            ncol = len(categories) // 2 if len(categories) in (6, 8) else (4 if len(categories) <= 6 else 2)
            if len(categories) in (6, 8):
                # Matplotlib fills columns first; keep the displayed rows sequential.
                handles = [handles[i] for col in range(ncol) for i in (col, col + ncol)]
            legend = ax.legend(handles=handles, title=pretty_label(key),
                      title_fontsize=LEGEND_TITLE_SIZE, bbox_to_anchor=(0.5, -0.25),
                      loc="upper center", ncol=ncol, frameon=False, fontsize=12,
                      columnspacing=1.0, handlelength=2.0, handletextpad=0.4)
            legends.append(legend)
    # Record the display-only correspondence; cluster IDs and scores are unchanged.
    output_path.with_suffix(".colors.json").write_text(json.dumps({
        "ground_truth": true_key, "predicted": pred_key,
        "matching": matched_labels(adata.obs, pred_key, true_key),
        "unmatched_color": "#B9C0C7",
    }, indent=2))

    # Preserve the full canvas for the two-row legends.
    fig.subplots_adjust(**PANEL_MARGINS)
    crop = Bbox.from_extents(0, PANEL_EXPORT_BOTTOM,
                             fig.get_figwidth(), fig.get_figheight())
    save_figure(fig, output_path, dpi=PANEL_DPI, bbox_inches=crop)
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
    ax.set_xlabel(f"Predicted cluster ({pretty_label(pred_key)})")
    ax.set_ylabel(f"Ground truth ({pretty_label(true_key)})")
    ax.set_title(f"{pretty_label(true_key)} vs. predicted cluster")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="fraction of row")
    fig.tight_layout()
    save_figure(fig, output_path, dpi=500, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    clustering_root = (
        SIM_PAPER_DIR / "data" / f"clustering_{args.packing_tag}" if args.packing_tag else CLUSTERING_ROOT
    )
    pipeline_dirname = "pca_harmony"
    if args.input == DEFAULT_INPUT:
        args.input = (
            clustering_root / args.modality / pipeline_dirname / f"simulation_{args.modality}_z_{pipeline_dirname}.h5ad"
        )
    if args.rep_key is None:
        args.rep_key = PIPELINE_REP_KEY

    if not args.input.is_file():
        raise SystemExit(
            f"Input file not found: {args.input}\n"
            f"Run sim_paper/code/clustering/step01_pca_harmony.py first."
        )

    tag = resolution_tag(args.resolution)
    pipeline_tag = PIPELINE_TAG
    if args.output_dir == DEFAULT_OUTPUT_DIR:
        args.output_dir = clustering_root / args.modality / f"{args.algorithm}_{pipeline_tag}_res{tag}"
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output = args.output_dir / f"simulation_{args.modality}_z_{args.algorithm}_{pipeline_tag}_res{tag}.h5ad"
    plot_dir = args.output_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    if args.plots_only:
        if not output.is_file():
            raise SystemExit(f"--plots-only requires an existing output file: {output}")
        print(f"[plots-only] Loading existing output {output}")
        adata = sc.read_h5ad(output)
        if "X_umap" not in adata.obsm or args.cluster_key not in adata.obs:
            raise KeyError(
                f"--plots-only requires 'X_umap' in obsm and {args.cluster_key!r} in obs; "
                f"{output} looks like it wasn't produced by this script."
            )
        valid_colors = [key for key in args.plot_colors if key in adata.obs]
        for color_key in valid_colors:
            plot_umap(adata, color_key, plot_dir / f"umap_by_{color_key}.png")
        for true_key in ("cell_type_true", "domain_true"):
            if true_key not in adata.obs or args.cluster_key not in adata.obs:
                continue
            plot_umap_true_vs_predicted(adata, true_key, args.cluster_key, plot_dir / f"umap_true_vs_predicted_{true_key}.png")
            plot_contingency_heatmap(adata, true_key, args.cluster_key, plot_dir / f"contingency_{true_key}.png")
        print("[plots-only] Done, not touching the existing output file")
        return

    print(f"[load] {args.input}")
    adata = sc.read_h5ad(args.input)

    if args.rep_key not in adata.obsm:
        available = ", ".join(adata.obsm.keys())
        raise KeyError(
            f"Representation {args.rep_key!r} was not found in adata.obsm. "
            f"Available obsm keys: {available}"
        )

    rep = np.asarray(adata.obsm[args.rep_key])
    # PCA components are variance-ordered, so truncating to the top n_pcs is meaningful.
    n_pcs = min(args.n_pcs, rep.shape[1])
    active_rep_key = f"{args.rep_key}_{n_pcs}"
    adata.obsm[active_rep_key] = rep[:, :n_pcs].copy()
    print(
        f"[neighbors] Using first {n_pcs} dimensions of {args.rep_key} "
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

    shared_umap_key = "X_umap_pca_post_harmony"
    if shared_umap_key in adata.obsm:
        print(
            f"[umap] Reusing {shared_umap_key} from the step01_pca_harmony.py output -- this figure's "
            "UMAP is then the exact same embedding as the before/after-Harmony panel "
            "(clustering/Leiden still uses its own neighbor graph above)."
        )
        adata.obsm["X_umap"] = np.asarray(adata.obsm[shared_umap_key]).copy()
    else:
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
