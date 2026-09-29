#!/usr/bin/env python
"""
Plot count-distribution diagnostics for an albis (or real) dataset.

Produces six PNGs:
  - mean_variance.png: per-gene mean vs variance (log-log), with the Poisson
    line (var = mean) and a single method-of-moments NB line
    (var = mean + mean^2/theta) overlaid, computed from the raw counts.
  - mean_dropout.png: per-gene mean vs observed zero fraction, with the
    Poisson- and NB-predicted zero-probability curves overlaid (same theta
    as mean_variance.png), computed from the raw counts. Tests whether the
    NB fit from the first two moments also explains the zero rate, or
    whether the data needs zero-inflation on top of NB.
  - total_counts.png: per-cell library-size (total counts) histogram,
    computed from the raw counts.
  - raw_norm_log.png: histogram of nonzero matrix entries at three pipeline
    stages -- raw counts, target-sum normalized, and log1p(normalized) --
    matching the normalize_total/log1p steps in clustering.py. Zeros are
    excluded at every stage, so this says nothing about sparsity -- see
    sparsity_summary.png for that.
  - genes_per_cell.png: per-cell number of distinct genes detected
    (nonzero entries per row), computed from the raw counts. A cell/bin
    could match total_counts.png well while still concentrating counts
    into fewer genes than real, which this catches and total_counts
    doesn't.
  - sparsity_summary.png: two headline sparsity numbers computed from the
    raw counts -- percent of all matrix entries that are zero, and percent
    of cells/bins that are fully empty (zero total counts). mean_dropout.png
    is conditional on each gene's mean, so it can't tell you the overall
    sparsity if gene means themselves differ between datasets; this does.

"Raw counts" means adata.X (post-batch-effect counts) by default, so the sim
side carries the same baked-in technical noise the real reference does --
apples-to-apples. Pass --no-batch-effect to instead use
adata.layers["counts_pre_batch"] when present (counts before the synthetic
batch-effect multiplier). The raw/norm/log panel always uses adata.X
regardless, since that is what the clustering pipeline actually consumes.

By default --input is also restricted to a single slice (--slice-id 4, the
5th of 10 slices counted from the bottom; slice_id 0 is the lowest z) rather
than pooling all 10 z-planes, matching the single-section nature of every
real reference. Changed from 5 (5th from the top) on 2026-09-29 so Figure 2
and Figure 5 use the same central slice. Pass --all-slices to pool.
--compare-input data is never filtered by slice (real data has no slice_id);
pre-subset a simulated --compare-input yourself.

Pass --compare-input (plus --compare-label) to overlay a second dataset --
e.g. real Xenium data -- on the same four diagnostics instead of plotting
--input alone. Both datasets are density-normalized so they're comparable
regardless of how many cells/genes each one has.

Expected environment:
    conda activate /dcs04/hicks/data/Jan/sim_project/albis/env/albis-tutorial

Example:
    python sim_paper/code/count_distribution/count_distribution.py --modality cell
    python sim_paper/code/count_distribution/count_distribution.py --modality cell \\
        --compare-input sim_paper/data/real_data_qc/non_diseased_lung/non_diseased_lung_qc.h5ad \\
        --compare-label non_diseased_lung
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import anndata as ad
import numpy as np
import scanpy as sc
from scipy import sparse
from scipy.spatial.distance import jensenshannon

SCRIPT_DIR = Path(__file__).resolve().parent
SIM_PAPER_DIR = SCRIPT_DIR.parents[1]

# Fixed categorical assignment (dataviz palette slots 1/2) -- color follows the
# dataset role (primary vs compare), never plot order.
PRIMARY_COLOR = "#2a78d6"  # blue
COMPARE_COLOR = "#eb6834"  # orange
# Comparison-plot title suffix; set from --title-context (default = Figure 2's).
TITLE_CONTEXT = "simulated vs real"

# Restored from typography_refresh_20260919/count_distribution_before.py,
# with the final refresh axis/base font size (16 pt), verified against saved plots.
# Sized for legibility once these PNGs are shrunk into a multi-panel print
# figure -- default matplotlib sizes (title ~12, legend ~10) read fine full-size
# on screen but wash out at print scale.
TITLE_SIZE = 17
LABEL_SIZE = 16
TICK_SIZE = 12
LEGEND_SIZE = 12
ANNOT_SIZE = 13
plt.rcParams.update({
    "font.size": LABEL_SIZE,
    "axes.titlesize": TITLE_SIZE,
    "axes.titleweight": "bold",
    "axes.labelsize": LABEL_SIZE,
    "xtick.labelsize": TICK_SIZE,
    "ytick.labelsize": TICK_SIZE,
    "legend.fontsize": LEGEND_SIZE,
    "figure.titlesize": TITLE_SIZE,
    "figure.titleweight": "bold",
})


def compute_jsd(x1: np.ndarray, x2: np.ndarray, bins: np.ndarray) -> float:
    """Jensen-Shannon divergence between the empirical distributions of x1 and x2,
    estimated from histograms over shared bins and normalized to probability mass
    functions. Symmetric, log base 2, bounded in [0, 1] (0 = identical, 1 =
    non-overlapping support). scipy's jensenshannon returns the JS *distance*
    (sqrt of the divergence), hence the squaring."""
    p, _ = np.histogram(x1, bins=bins)
    q, _ = np.histogram(x2, bins=bins)
    p = p / p.sum()
    q = q / q.sum()
    return float(jensenshannon(p, q, base=2) ** 2)


def display_label(raw: str) -> str:
    """Human-readable form of a --modality/--compare-label value for plot text, e.g.
    "breast_cancer_visium_hd" -> "Breast Cancer Visium HD". The raw string itself is left
    untouched everywhere else (JSON summary keys, prints, file paths)."""
    acronyms = {"hd"}
    words = [w.upper() if w.lower() in acronyms else w[0].upper() + w[1:] for w in raw.split("_") if w]
    return " ".join(words)


def annotate_jsd(ax, jsd: float, loc: str = "upper left") -> None:
    x, ha = (0.03, "left") if "left" in loc else (0.97, "right")
    y, va = (0.97, "top") if "upper" in loc else (0.03, "bottom")
    ax.text(
        x, y, f"JSD = {jsd:.3f}", transform=ax.transAxes, ha=ha, va=va,
        fontsize=ANNOT_SIZE, fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="0.6", alpha=0.85),
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot per-gene mean-variance, total-counts, and raw/norm/log1p diagnostics."
    )
    parser.add_argument("--input", type=Path, default=None)
    parser.add_argument(
        "--modality",
        default="spot",
        help="Label used for default input/output paths and plot titles "
        "(spot/bin/cell for albis data; any label when --input is a real dataset).",
    )
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--target-sum", type=float, default=1e4)
    parser.add_argument(
        "--sample-size",
        type=int,
        default=2_000_000,
        help="Max nonzero matrix entries sampled for the raw/norm/log1p histograms.",
    )
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument(
        "--compare-input",
        type=Path,
        default=None,
        help="Second dataset (e.g. real Xenium data) to overlay against --input instead of "
        "plotting --input alone.",
    )
    parser.add_argument(
        "--compare-label",
        default="real",
        help="Label used for --compare-input in plot legends/titles and default output path.",
    )
    parser.add_argument(
        "--slice-id",
        type=int,
        default=4,
        help="Restrict --input (the sim/primary dataset) to this single obs['slice_id'] value "
        "before computing any stats or plots, instead of pooling all slices. Defaults to 4 (the "
        "5th of 10 slices from the bottom, a central section); pass --all-slices to pool instead. "
        "--compare-input data is never filtered by this flag.",
    )
    parser.add_argument(
        "--title-context",
        default="simulated vs real",
        help="Text after the colon in comparison-plot titles, e.g. 'ALBIS vs SPIDER' for a "
        "method-vs-method comparison. Default keeps Figure 2's titles.",
    )
    parser.add_argument(
        "--all-slices",
        action="store_true",
        help="Pool every slice of --input instead of restricting to --slice-id. Overrides --slice-id.",
    )
    parser.add_argument(
        "--no-batch-effect",
        dest="use_batch_effect",
        action="store_false",
        default=True,
        help="Use --input's 'counts_pre_batch' layer (counts before the synthetic batch-effect "
        "multiplier) when present, instead of the default X (post-batch-effect counts). By default "
        "X is used uniformly for every stat/plot so the sim side matches real --compare-input data, "
        "which has no pre-batch version -- whatever technical noise it carries is just baked into "
        "its counts -- making it apples-to-apples rather than artificially cleaner than real can "
        "ever be. (The raw_norm_log_compare panel always used X regardless.)",
    )
    parser.add_argument(
        "--match-panel-size",
        action="store_true",
        help="HVG-subset whichever of --input/--compare-input has the larger gene panel down to "
        "the other's gene count, before computing any comparison metric. Usually shrinks the real "
        "reference (bin/spot vs. Visium(HD): 18,085-36,601 real genes vs. 556 sim genes), but is "
        "symmetric -- for cell vs. Xenium, sim's panel (556) is actually larger than Xenium's "
        "(392), so this shrinks sim instead. Deliberately NOT a random subsample -- HVG selection "
        "keeps the marker/signal-heavy genes on the larger side, giving the closer like-for-like "
        "comparison, and it directly removes apparent total-counts/sparsity gaps that are really "
        "just panel-size artifacts rather than simulator fidelity issues.",
    )
    args = parser.parse_args()

    if args.all_slices:
        args.slice_id = None

    if args.input is None:
        args.input = SIM_PAPER_DIR / "data" / f"simulation_{args.modality}_z.h5ad"
    if args.output_dir is None:
        if args.compare_input is None:
            args.output_dir = SIM_PAPER_DIR / "data" / "count_distribution" / args.modality
        else:
            args.output_dir = (
                SIM_PAPER_DIR / "data" / "count_distribution" / f"{args.modality}_vs_{args.compare_label}"
            )
    return args


def gene_mean_var(X) -> tuple[np.ndarray, np.ndarray]:
    if sparse.issparse(X):
        mean = np.asarray(X.mean(axis=0)).ravel()
        mean_sq = np.asarray(X.multiply(X).mean(axis=0)).ravel()
    else:
        mean = X.mean(axis=0)
        mean_sq = (X**2).mean(axis=0)
    var = mean_sq - mean**2
    return mean, var


def fit_common_dispersion(mean: np.ndarray, var: np.ndarray) -> float:
    overdispersed = var > mean
    theta_hat = mean[overdispersed] ** 2 / (var[overdispersed] - mean[overdispersed])
    theta_hat = theta_hat[np.isfinite(theta_hat) & (theta_hat > 0)]
    return float(np.median(theta_hat)) if theta_hat.size else float("nan")


def plot_mean_variance(mean: np.ndarray, var: np.ndarray, theta: float, modality: str, output_path: Path) -> None:
    keep = (mean > 0) & (var > 0)
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    ax.scatter(mean[keep], var[keep], s=8, alpha=0.5, linewidths=0, label="genes")

    x = np.logspace(np.log10(mean[keep].min()), np.log10(mean[keep].max()), 200)
    ax.plot(x, x, "k--", label="Poisson (var = mean)")
    if np.isfinite(theta):
        ax.plot(x, x + x**2 / theta, "r-", label=rf"NB fit ($\hat\theta$={theta:.1f})")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Mean count per gene")
    ax.set_ylabel("Variance per gene")
    ax.set_title(f"{modality}: gene mean-variance")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def gene_zero_fraction(X, n_obs: int) -> np.ndarray:
    nnz = X.getnnz(axis=0) if sparse.issparse(X) else np.count_nonzero(X, axis=0)
    return 1.0 - np.asarray(nnz).ravel() / n_obs


def plot_mean_dropout(mean: np.ndarray, zero_frac: np.ndarray, theta: float, modality: str, output_path: Path) -> None:
    keep = mean > 0
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    ax.scatter(mean[keep], zero_frac[keep], s=8, alpha=0.5, linewidths=0, label="genes")

    x = np.logspace(np.log10(mean[keep].min()), np.log10(mean[keep].max()), 200)
    ax.plot(x, np.exp(-x), "k--", label="Poisson-predicted")
    if np.isfinite(theta):
        ax.plot(x, (theta / (theta + x)) ** theta, "r-", label=rf"NB-predicted ($\hat\theta$={theta:.1f})")

    ax.set_xscale("log")
    ax.set_xlabel("Mean count per gene")
    ax.set_ylabel("Fraction of zero cells")
    ax.set_title(f"{modality}: gene mean-dropout")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_total_counts(total_counts: np.ndarray, modality: str, output_path: Path) -> None:
    positive = total_counts[total_counts > 0]
    bins = np.logspace(np.log10(positive.min()), np.log10(positive.max()), 60)

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.hist(positive, bins=bins, color="steelblue", edgecolor="none")
    ax.set_xscale("log")
    ax.set_xlabel("Total counts per cell")
    ax.set_ylabel("Number of cells")
    ax.set_title(f"{modality}: total counts per cell")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def genes_per_cell(X) -> np.ndarray:
    return np.asarray(X.getnnz(axis=1) if sparse.issparse(X) else np.count_nonzero(X, axis=1)).ravel()


def plot_genes_per_cell(n_genes: np.ndarray, modality: str, output_path: Path) -> None:
    positive = n_genes[n_genes > 0]
    bins = np.logspace(np.log10(positive.min()), np.log10(positive.max()), 60)

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.hist(positive, bins=bins, color="steelblue", edgecolor="none")
    ax.set_xscale("log")
    ax.set_xlabel("Genes detected per cell")
    ax.set_ylabel("Number of cells")
    ax.set_title(f"{modality}: genes detected per cell")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def compute_sparsity_stats(X, n_genes: np.ndarray) -> dict[str, float]:
    """n_genes: per-cell nonzero-gene counts from genes_per_cell(X), reused rather than
    recomputed. matrix_zero_frac is the overall fraction of zero entries in X; empty_row_frac
    is the fraction of cells/bins with zero genes detected."""
    matrix_zero_frac = 1.0 - n_genes.sum() / (X.shape[0] * X.shape[1])
    empty_row_frac = float(np.mean(n_genes == 0))
    return {"matrix_zero_frac": float(matrix_zero_frac), "empty_row_frac": empty_row_frac}


def plot_sparsity_summary(stats: dict[str, float], modality: str, output_path: Path) -> None:
    labels = ["Zero matrix\nentries", "Fully-empty\ncells"]
    values = [stats["matrix_zero_frac"] * 100, stats["empty_row_frac"] * 100]

    fig, ax = plt.subplots(figsize=(5, 5.5))
    bars = ax.bar(labels, values, color="steelblue", width=0.5)
    for bar, v in zip(bars, values):
        ax.annotate(f"{v:.1f}%", (bar.get_x() + bar.get_width() / 2, v), ha="center", va="bottom",
                    fontsize=ANNOT_SIZE, fontweight="bold")
    ax.set_ylabel("Percent (%)")
    ax.set_ylim(0, 105)
    ax.set_title(f"{modality}: sparsity summary")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def sample_values(data: np.ndarray, sample_size: int, rng: np.random.Generator) -> np.ndarray:
    if data.size <= sample_size:
        return data
    idx = rng.choice(data.size, size=sample_size, replace=False)
    return data[idx]


def plot_raw_norm_log(
    X, target_sum: float, sample_size: int, rng: np.random.Generator, modality: str, output_path: Path
) -> None:
    tmp = ad.AnnData(X=X.copy())
    raw_vals = sample_values(tmp.X.data.copy(), sample_size, rng)

    sc.pp.normalize_total(tmp, target_sum=target_sum)
    norm_vals = sample_values(tmp.X.data.copy(), sample_size, rng)

    sc.pp.log1p(tmp)
    log_vals = sample_values(tmp.X.data.copy(), sample_size, rng)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    stages = [
        (axes[0], raw_vals, "Raw counts", True),
        (axes[1], norm_vals, f"Normalized (target_sum={target_sum:g})", True),
        (axes[2], log_vals, "log1p(normalized)", False),
    ]
    for ax, vals, title, log_x in stages:
        vals = vals[vals > 0]
        if log_x:
            bins = np.logspace(np.log10(vals.min()), np.log10(vals.max()), 60)
            ax.set_xscale("log")
        else:
            bins = 60
        ax.hist(vals, bins=bins, color="darkorange", edgecolor="none")
        ax.set_title(title)
        ax.set_xlabel("Value (nonzero matrix entries)")
    axes[0].set_ylabel("Count")
    fig.suptitle(f"{modality}: raw -> normalized -> log1p")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_mean_variance_compare(datasets: list[dict], output_path: Path) -> None:
    """datasets: [{"label", "display_label", "color", "mean", "var", "theta"}, ...] (primary first)."""
    fig, ax = plt.subplots(figsize=(8, 6.5))

    all_mean = np.concatenate([d["mean"][(d["mean"] > 0) & (d["var"] > 0)] for d in datasets])
    x = np.logspace(np.log10(all_mean.min()), np.log10(all_mean.max()), 200)
    ax.plot(x, x, "k--", label="Poisson (var = mean)", zorder=1)

    for d in datasets:
        keep = (d["mean"] > 0) & (d["var"] > 0)
        ax.scatter(d["mean"][keep], d["var"][keep], s=8, alpha=0.4, linewidths=0, color=d["color"], label=d["display_label"])
        if np.isfinite(d["theta"]):
            ax.plot(x, x + x**2 / d["theta"], "-", color=d["color"], linewidth=1.5,
                     label=rf"{d['display_label']} NB fit ($\hat\theta$={d['theta']:.1f})")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Mean count per gene")
    ax.set_ylabel("Variance per gene")
    ax.set_title(f"Gene mean-variance: {TITLE_CONTEXT}")
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_mean_dropout_compare(datasets: list[dict], output_path: Path) -> None:
    """datasets: [{"label", "display_label", "color", "mean", "zero_frac", "theta"}, ...] (primary first)."""
    fig, ax = plt.subplots(figsize=(8, 6.5))

    all_mean = np.concatenate([d["mean"][d["mean"] > 0] for d in datasets])
    x = np.logspace(np.log10(all_mean.min()), np.log10(all_mean.max()), 200)
    ax.plot(x, np.exp(-x), "k--", label="Poisson-predicted", zorder=1)

    for d in datasets:
        keep = d["mean"] > 0
        ax.scatter(d["mean"][keep], d["zero_frac"][keep], s=8, alpha=0.4, linewidths=0, color=d["color"], label=d["display_label"])
        if np.isfinite(d["theta"]):
            ax.plot(x, (d["theta"] / (d["theta"] + x)) ** d["theta"], "-", color=d["color"], linewidth=1.5,
                     label=rf"{d['display_label']} NB-predicted ($\hat\theta$={d['theta']:.1f})")

    ax.set_xscale("log")
    ax.set_xlabel("Mean count per gene")
    ax.set_ylabel("Fraction of zero cells")
    ax.set_title(f"Gene mean-dropout: {TITLE_CONTEXT}")
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_total_counts_compare(datasets: list[dict], output_path: Path) -> None:
    """datasets: [{"label", "display_label", "color", "total_counts"}, ...] (primary first).
    Density-normalized so datasets with very different cell counts are still comparable by
    shape. No legend is drawn here -- see plot_dataset_legend for a standalone legend image
    to insert by hand instead."""
    positive = [d["total_counts"][d["total_counts"] > 0] for d in datasets]
    all_positive = np.concatenate(positive)
    bins = np.logspace(np.log10(all_positive.min()), np.log10(all_positive.max()), 60)

    fig, ax = plt.subplots(figsize=(7, 6.5))
    for d, pos in zip(datasets, positive):
        ax.hist(pos, bins=bins, density=True, color=d["color"], alpha=0.5, label=d["display_label"])
    ax.set_xscale("log")
    ax.set_xlabel("Total counts per cell")
    ax.set_ylabel("Density (fraction of cells)")
    ax.set_title(f"Total counts per cell: {TITLE_CONTEXT}")
    annotate_jsd(ax, compute_jsd(positive[0], positive[1], bins), loc="upper left")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_genes_per_cell_compare(datasets: list[dict], output_path: Path) -> None:
    """datasets: [{"label", "display_label", "color", "n_genes"}, ...] (primary first)."""
    positive = [d["n_genes"][d["n_genes"] > 0] for d in datasets]
    all_positive = np.concatenate(positive)
    bins = np.logspace(np.log10(all_positive.min()), np.log10(all_positive.max()), 60)

    fig, ax = plt.subplots(figsize=(7, 6.5))
    for d, pos in zip(datasets, positive):
        ax.hist(pos, bins=bins, density=True, color=d["color"], alpha=0.5, label=d["display_label"])
    ax.set_xscale("log")
    ax.set_xlabel("Genes detected per cell")
    ax.set_ylabel("Density (fraction of cells)")
    ax.set_title(f"Genes detected per cell: {TITLE_CONTEXT}")
    ax.legend(frameon=False, loc="upper right")
    annotate_jsd(ax, compute_jsd(positive[0], positive[1], bins), loc="upper left")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_sparsity_summary_compare(datasets: list[dict], output_path: Path) -> None:
    """datasets: [{"label", "display_label", "color", "sparsity"}, ...] (primary first), sparsity
    is the dict returned by compute_sparsity_stats."""
    metrics = [("matrix_zero_frac", "Zero matrix\nentries"), ("empty_row_frac", "Fully-empty\ncells")]
    x = np.arange(len(metrics))
    width = 0.8 / len(datasets)

    fig, ax = plt.subplots(figsize=(6.5, 6))
    for i, d in enumerate(datasets):
        values = [d["sparsity"][key] * 100 for key, _ in metrics]
        offset = (i - (len(datasets) - 1) / 2) * width
        bars = ax.bar(x + offset, values, width=width, color=d["color"], label=d["display_label"])
        for bar, v in zip(bars, values):
            ax.annotate(f"{v:.1f}%", (bar.get_x() + bar.get_width() / 2, v), ha="center", va="bottom",
                        fontsize=ANNOT_SIZE, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels([label for _, label in metrics])
    ax.set_ylabel("Percent (%)")
    ax.set_ylim(0, 105)
    ax.set_title(f"Sparsity summary: {TITLE_CONTEXT}")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_dataset_legend(datasets: list[dict], output_path: Path) -> None:
    """Standalone legend image (dataset display-name -> color patch), for hand-insertion into
    total_counts_compare.png and raw_norm_log_compare.png, whose inline legends are omitted to
    avoid crowding those panels. datasets: [{"display_label", "color"}, ...] (primary first)."""
    handles = [Patch(facecolor=d["color"], alpha=0.5, label=d["display_label"]) for d in datasets]
    fig, ax = plt.subplots(figsize=(3.5, 0.6 + 0.45 * len(handles)))
    ax.axis("off")
    ax.legend(handles=handles, loc="center", frameon=False)
    fig.savefig(output_path, dpi=180, bbox_inches="tight", pad_inches=0.02, transparent=True)
    plt.close(fig)


def plot_raw_norm_log_compare(
    datasets: list[dict], target_sum: float, sample_size: int, rng: np.random.Generator, output_path: Path
) -> None:
    """datasets: [{"label", "display_label", "color", "X"}, ...] (primary first). No legend is
    drawn here -- see plot_dataset_legend for a standalone legend image to insert by hand instead."""
    staged = []
    for d in datasets:
        tmp = ad.AnnData(X=d["X"].copy())
        raw_vals = sample_values(tmp.X.data.copy(), sample_size, rng)
        sc.pp.normalize_total(tmp, target_sum=target_sum)
        norm_vals = sample_values(tmp.X.data.copy(), sample_size, rng)
        sc.pp.log1p(tmp)
        log_vals = sample_values(tmp.X.data.copy(), sample_size, rng)
        staged.append({"display_label": d["display_label"], "color": d["color"], "raw": raw_vals, "norm": norm_vals, "log": log_vals})

    fig, axes = plt.subplots(1, 3, figsize=(17, 6))
    stage_specs = [
        (axes[0], "raw", "Raw counts", True),
        (axes[1], "norm", "Normalized", True),
        (axes[2], "log", "log1p(normalized)", False),
    ]
    for ax, key, title, log_x in stage_specs:
        vals_by_dataset = [s[key][s[key] > 0] for s in staged]
        all_vals = np.concatenate(vals_by_dataset)
        if log_x:
            bins = np.logspace(np.log10(all_vals.min()), np.log10(all_vals.max()), 60)
            ax.set_xscale("log")
        else:
            bins = np.linspace(all_vals.min(), all_vals.max(), 60)
        for s, vals in zip(staged, vals_by_dataset):
            ax.hist(vals, bins=bins, density=True, color=s["color"], alpha=0.5, label=s["display_label"])
        ax.set_title(title)
        ax.set_xlabel("Value (nonzero matrix entries)")
        annotate_jsd(ax, compute_jsd(vals_by_dataset[0], vals_by_dataset[1], bins), loc="upper right")
    axes[0].set_ylabel("Density")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def select_counts(adata: ad.AnnData, use_batch_effect: bool):
    if use_batch_effect:
        return adata.X
    return adata.layers["counts_pre_batch"] if "counts_pre_batch" in adata.layers else adata.X


def restrict_to_slice(adata: ad.AnnData, slice_id: int) -> ad.AnnData:
    if "slice_id" not in adata.obs.columns:
        raise SystemExit(f"--slice-id {slice_id} given but 'slice_id' not found in obs columns")
    keep = adata.obs["slice_id"] == slice_id
    if not keep.any():
        raise SystemExit(f"--slice-id {slice_id} matches no observations (available: "
                          f"{sorted(adata.obs['slice_id'].unique().tolist())})")
    print(f"[slice] restricting to slice_id={slice_id}: {adata.n_obs} -> {int(keep.sum())} observations")
    return adata[keep].copy()


def run_compare(args: argparse.Namespace, rng: np.random.Generator) -> None:
    global TITLE_CONTEXT
    TITLE_CONTEXT = args.title_context
    print(f"[load] primary: {args.input}")
    primary = sc.read_h5ad(args.input)
    print(f"[load] primary shape: {primary.n_obs} observations x {primary.n_vars} genes")
    if args.slice_id is not None:
        primary = restrict_to_slice(primary, args.slice_id)

    print(f"[load] compare: {args.compare_input}")
    compare = sc.read_h5ad(args.compare_input)
    print(f"[load] compare shape: {compare.n_obs} observations x {compare.n_vars} genes")

    if args.match_panel_size and compare.n_vars != primary.n_vars:
        # Whichever side has the larger panel gets HVG-subsetted down to the
        # smaller side's gene count -- e.g. compare (real) is much larger for
        # bin/spot vs. Visium(HD), but primary (sim) is actually the larger
        # panel for cell vs. Xenium (556 vs. 392 genes), so this must be
        # symmetric rather than always shrinking `compare`.
        if compare.n_vars > primary.n_vars:
            larger_name, n_target = "compare", primary.n_vars
        else:
            larger_name, n_target = "primary", compare.n_vars
        larger = compare if larger_name == "compare" else primary
        print(f"[panel-match] subsetting {larger_name} from {larger.n_vars} to top-{n_target} "
              f"highly-variable genes (matching the smaller panel)")
        larger_counts = select_counts(larger, args.use_batch_effect)
        hvg_source = ad.AnnData(X=larger_counts.copy(), var=larger.var.copy())
        sc.pp.highly_variable_genes(hvg_source, n_top_genes=n_target, flavor="seurat_v3")
        larger = larger[:, hvg_source.var["highly_variable"]].copy()
        if larger_name == "compare":
            compare = larger
        else:
            primary = larger
        print(f"[panel-match] {larger_name} shape after HVG subset: {larger.n_obs} observations x {larger.n_vars} genes")
    elif args.match_panel_size:
        print(f"[panel-match] skipped -- compare ({compare.n_vars} genes) and {args.modality} "
              f"({primary.n_vars} genes) are already the same size")

    specs = [
        {"label": args.modality, "display_label": display_label(args.modality), "color": PRIMARY_COLOR, "adata": primary},
        {"label": args.compare_label, "display_label": display_label(args.compare_label), "color": COMPARE_COLOR, "adata": compare},
    ]
    for spec in specs:
        adata = spec["adata"]
        raw_counts = select_counts(adata, args.use_batch_effect)
        mean, var = gene_mean_var(raw_counts)
        spec["mean"] = mean
        spec["var"] = var
        spec["theta"] = fit_common_dispersion(mean, var)
        spec["zero_frac"] = gene_zero_fraction(raw_counts, adata.n_obs)
        spec["total_counts"] = np.asarray(raw_counts.sum(axis=1)).ravel()
        spec["n_genes"] = genes_per_cell(raw_counts)
        spec["sparsity"] = compute_sparsity_stats(raw_counts, spec["n_genes"])
        used = "X" if (args.use_batch_effect or "counts_pre_batch" not in adata.layers) else "counts_pre_batch"
        print(f"[nb] {spec['label']}: theta_hat = {spec['theta']:.2f}, {used} used for mean-variance/total-counts")
        print(f"[sparsity] {spec['label']}: {spec['sparsity']['matrix_zero_frac'] * 100:.2f}% zero entries, "
              f"{spec['sparsity']['empty_row_frac'] * 100:.2f}% fully-empty cells")

    plot_mean_variance_compare(specs, args.output_dir / "mean_variance_compare.png")
    plot_mean_dropout_compare(specs, args.output_dir / "mean_dropout_compare.png")
    plot_total_counts_compare(specs, args.output_dir / "total_counts_compare.png")
    plot_genes_per_cell_compare(specs, args.output_dir / "genes_per_cell_compare.png")
    plot_sparsity_summary_compare(specs, args.output_dir / "sparsity_summary_compare.png")
    plot_raw_norm_log_compare(
        [{"display_label": s["display_label"], "color": s["color"], "X": s["adata"].X} for s in specs],
        args.target_sum, args.sample_size, rng, args.output_dir / "raw_norm_log_compare.png",
    )
    plot_dataset_legend(specs, args.output_dir / "dataset_legend.png")
    print(f"[save] Comparison plots written to {args.output_dir}")

    # Numeric summary alongside the plots -- lets a sweep across many configs be scored
    # programmatically (theta_hat, total_counts, empty-bin fraction, ...) instead of having
    # to eyeball every PNG, which is how prior sweeps (see figure.md/DATA_VERSIONS.md bin
    # packing-fraction x log_mu history) were tracked by hand.
    summary = {
        "modality": args.modality,
        "compare_label": args.compare_label,
        "match_panel_size": bool(args.match_panel_size),
        "use_batch_effect": bool(args.use_batch_effect),
        "slice_id": args.slice_id,
        "n_obs": {s["label"]: int(s["adata"].n_obs) for s in specs},
        "n_vars": {s["label"]: int(s["adata"].n_vars) for s in specs},
        "theta_hat": {s["label"]: s["theta"] for s in specs},
        "total_counts_median": {s["label"]: float(np.median(s["total_counts"])) for s in specs},
        "total_counts_mean": {s["label"]: float(np.mean(s["total_counts"])) for s in specs},
        "genes_per_cell_median": {s["label"]: float(np.median(s["n_genes"])) for s in specs},
        "matrix_zero_frac": {s["label"]: s["sparsity"]["matrix_zero_frac"] for s in specs},
        "empty_row_frac": {s["label"]: s["sparsity"]["empty_row_frac"] for s in specs},
    }
    primary_label, compare_label = args.modality, args.compare_label
    summary["ratio_primary_over_compare"] = {
        "theta_hat": summary["theta_hat"][primary_label] / summary["theta_hat"][compare_label]
        if summary["theta_hat"][compare_label] not in (0, None) and np.isfinite(summary["theta_hat"][compare_label])
        else None,
        "total_counts_median": summary["total_counts_median"][primary_label] / summary["total_counts_median"][compare_label]
        if summary["total_counts_median"][compare_label] else None,
    }
    summary_path = args.output_dir / "comparison_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[save] Numeric summary written to {summary_path}")


def main() -> None:
    args = parse_args()
    if not args.input.is_file():
        raise SystemExit(f"Input file not found: {args.input}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.random_state)

    if args.compare_input is not None:
        if not args.compare_input.is_file():
            raise SystemExit(f"Compare input file not found: {args.compare_input}")
        run_compare(args, rng)
        return

    print(f"[load] {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"[load] AnnData shape: {adata.n_obs} observations x {adata.n_vars} genes")
    if args.slice_id is not None:
        adata = restrict_to_slice(adata, args.slice_id)

    raw_counts = select_counts(adata, args.use_batch_effect)
    used = "X" if (args.use_batch_effect or "counts_pre_batch" not in adata.layers) else "counts_pre_batch"
    print(f"[nb] Using {used} for mean-variance/total-counts")

    print("[mean_variance] Computing per-gene mean/variance and fitting a common NB dispersion")
    mean, var = gene_mean_var(raw_counts)
    theta = fit_common_dispersion(mean, var)
    print(f"[mean_variance] theta_hat (median, method of moments) = {theta:.2f}")
    plot_mean_variance(mean, var, theta, args.modality, args.output_dir / "mean_variance.png")

    print("[mean_dropout] Computing per-gene zero fraction")
    zero_frac = gene_zero_fraction(raw_counts, adata.n_obs)
    plot_mean_dropout(mean, zero_frac, theta, args.modality, args.output_dir / "mean_dropout.png")

    print("[total_counts] Computing per-cell total counts")
    total_counts = np.asarray(raw_counts.sum(axis=1)).ravel()
    plot_total_counts(total_counts, args.modality, args.output_dir / "total_counts.png")

    print("[genes_per_cell] Computing per-cell genes-detected")
    n_genes = genes_per_cell(raw_counts)
    plot_genes_per_cell(n_genes, args.modality, args.output_dir / "genes_per_cell.png")

    print("[sparsity] Computing global sparsity stats")
    sparsity = compute_sparsity_stats(raw_counts, n_genes)
    print(f"[sparsity] {sparsity['matrix_zero_frac'] * 100:.2f}% zero entries, "
          f"{sparsity['empty_row_frac'] * 100:.2f}% fully-empty cells")
    plot_sparsity_summary(sparsity, args.modality, args.output_dir / "sparsity_summary.png")

    print("[raw_norm_log] Sampling nonzero entries through normalize_total -> log1p")
    plot_raw_norm_log(
        adata.X, args.target_sum, args.sample_size, rng, args.modality, args.output_dir / "raw_norm_log.png"
    )

    print(f"[save] Plots written to {args.output_dir}")


if __name__ == "__main__":
    main()
