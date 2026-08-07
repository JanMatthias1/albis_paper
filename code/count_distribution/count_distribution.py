#!/usr/bin/env python
"""
Plot count-distribution diagnostics for a sim_app (or real) dataset.

Produces three PNGs:
  - mean_variance.png: per-gene mean vs variance (log-log), with the Poisson
    line (var = mean) and a single method-of-moments NB line
    (var = mean + mean^2/theta) overlaid, computed from the raw counts.
  - total_counts.png: per-cell library-size (total counts) histogram,
    computed from the raw counts.
  - raw_norm_log.png: histogram of nonzero matrix entries at three pipeline
    stages -- raw counts, target-sum normalized, and log1p(normalized) --
    matching the normalize_total/log1p steps in clustering.py.

"Raw counts" means adata.layers["counts_pre_batch"] when present (the counts
before the synthetic batch-effect multiplier is applied), otherwise adata.X.
The raw/norm/log panel always uses adata.X, since that is what the
clustering pipeline actually consumes.

Expected environment:
    conda activate /dcs04/hicks/data/Jan/sim_project/sim_app/env/sim-app-tutorial

Example:
    python sim_paper/code/count_distribution/count_distribution.py --modality cell
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import anndata as ad
import numpy as np
import scanpy as sc
from scipy import sparse

SCRIPT_DIR = Path(__file__).resolve().parent
SIM_PAPER_DIR = SCRIPT_DIR.parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot per-gene mean-variance, total-counts, and raw/norm/log1p diagnostics."
    )
    parser.add_argument("--input", type=Path, default=None)
    parser.add_argument(
        "--modality",
        default="spot",
        help="Label used for default input/output paths and plot titles "
        "(spot/bin/cell for sim_app data; any label when --input is a real dataset).",
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
    args = parser.parse_args()

    if args.input is None:
        args.input = SIM_PAPER_DIR / "data" / f"simulation_{args.modality}_z.h5ad"
    if args.output_dir is None:
        args.output_dir = SIM_PAPER_DIR / "data" / "count_distribution" / args.modality
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
    fig, ax = plt.subplots(figsize=(6, 5))
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
    fig.savefig(output_path, dpi=180)
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
    fig.savefig(output_path, dpi=180)
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
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    if not args.input.is_file():
        raise SystemExit(f"Input file not found: {args.input}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.random_state)

    print(f"[load] {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"[load] AnnData shape: {adata.n_obs} observations x {adata.n_vars} genes")

    raw_counts = adata.layers["counts_pre_batch"] if "counts_pre_batch" in adata.layers else adata.X
    print(f"[nb] Using {'counts_pre_batch' if 'counts_pre_batch' in adata.layers else 'X'} for mean-variance/total-counts")

    print("[mean_variance] Computing per-gene mean/variance and fitting a common NB dispersion")
    mean, var = gene_mean_var(raw_counts)
    theta = fit_common_dispersion(mean, var)
    print(f"[mean_variance] theta_hat (median, method of moments) = {theta:.2f}")
    plot_mean_variance(mean, var, theta, args.modality, args.output_dir / "mean_variance.png")

    print("[total_counts] Computing per-cell total counts")
    total_counts = np.asarray(raw_counts.sum(axis=1)).ravel()
    plot_total_counts(total_counts, args.modality, args.output_dir / "total_counts.png")

    print("[raw_norm_log] Sampling nonzero entries through normalize_total -> log1p")
    plot_raw_norm_log(
        adata.X, args.target_sum, args.sample_size, rng, args.modality, args.output_dir / "raw_norm_log.png"
    )

    print(f"[save] Plots written to {args.output_dir}")


if __name__ == "__main__":
    main()
