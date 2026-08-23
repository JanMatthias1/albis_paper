#!/usr/bin/env python
"""
Plot count-distribution diagnostics for a sim_app (or real) dataset.

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

"Raw counts" means adata.layers["counts_pre_batch"] when present (the counts
before the synthetic batch-effect multiplier is applied), otherwise adata.X.
The raw/norm/log panel always uses adata.X, since that is what the
clustering pipeline actually consumes.

Pass --compare-input (plus --compare-label) to overlay a second dataset --
e.g. real Xenium data -- on the same four diagnostics instead of plotting
--input alone. Both datasets are density-normalized so they're comparable
regardless of how many cells/genes each one has.

Expected environment:
    conda activate /dcs04/hicks/data/Jan/sim_project/sim_app/env/sim-app-tutorial

Example:
    python sim_paper/code/count_distribution/count_distribution.py --modality cell
    python sim_paper/code/count_distribution/count_distribution.py --modality cell \\
        --compare-input sim_paper/data/real_data_qc/non_diseased_lung/non_diseased_lung_qc.h5ad \\
        --compare-label non_diseased_lung
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

# Fixed categorical assignment (dataviz palette slots 1/2) -- color follows the
# dataset role (primary vs compare), never plot order.
PRIMARY_COLOR = "#2a78d6"  # blue
COMPARE_COLOR = "#eb6834"  # orange


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
    args = parser.parse_args()

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


def gene_zero_fraction(X, n_obs: int) -> np.ndarray:
    nnz = X.getnnz(axis=0) if sparse.issparse(X) else np.count_nonzero(X, axis=0)
    return 1.0 - np.asarray(nnz).ravel() / n_obs


def plot_mean_dropout(mean: np.ndarray, zero_frac: np.ndarray, theta: float, modality: str, output_path: Path) -> None:
    keep = mean > 0
    fig, ax = plt.subplots(figsize=(6, 5))
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
    fig.savefig(output_path, dpi=180)
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

    fig, ax = plt.subplots(figsize=(4.5, 5))
    bars = ax.bar(labels, values, color="steelblue", width=0.5)
    for bar, v in zip(bars, values):
        ax.annotate(f"{v:.1f}%", (bar.get_x() + bar.get_width() / 2, v), ha="center", va="bottom", fontsize=9)
    ax.set_ylabel("Percent (%)")
    ax.set_ylim(0, 105)
    ax.set_title(f"{modality}: sparsity summary")
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


def plot_mean_variance_compare(datasets: list[dict], output_path: Path) -> None:
    """datasets: [{"label", "color", "mean", "var", "theta"}, ...] (primary first)."""
    fig, ax = plt.subplots(figsize=(6.5, 5.5))

    all_mean = np.concatenate([d["mean"][(d["mean"] > 0) & (d["var"] > 0)] for d in datasets])
    x = np.logspace(np.log10(all_mean.min()), np.log10(all_mean.max()), 200)
    ax.plot(x, x, "--", color="#898781", label="Poisson (var = mean)", zorder=1)

    for d in datasets:
        keep = (d["mean"] > 0) & (d["var"] > 0)
        ax.scatter(d["mean"][keep], d["var"][keep], s=8, alpha=0.4, linewidths=0, color=d["color"], label=d["label"])
        if np.isfinite(d["theta"]):
            ax.plot(x, x + x**2 / d["theta"], "-", color=d["color"], linewidth=1.5,
                     label=rf"{d['label']} NB fit ($\hat\theta$={d['theta']:.1f})")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Mean count per gene")
    ax.set_ylabel("Variance per gene")
    ax.set_title("Gene mean-variance: simulated vs real")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_mean_dropout_compare(datasets: list[dict], output_path: Path) -> None:
    """datasets: [{"label", "color", "mean", "zero_frac", "theta"}, ...] (primary first)."""
    fig, ax = plt.subplots(figsize=(6.5, 5.5))

    all_mean = np.concatenate([d["mean"][d["mean"] > 0] for d in datasets])
    x = np.logspace(np.log10(all_mean.min()), np.log10(all_mean.max()), 200)
    ax.plot(x, np.exp(-x), "--", color="#898781", label="Poisson-predicted", zorder=1)

    for d in datasets:
        keep = d["mean"] > 0
        ax.scatter(d["mean"][keep], d["zero_frac"][keep], s=8, alpha=0.4, linewidths=0, color=d["color"], label=d["label"])
        if np.isfinite(d["theta"]):
            ax.plot(x, (d["theta"] / (d["theta"] + x)) ** d["theta"], "-", color=d["color"], linewidth=1.5,
                     label=rf"{d['label']} NB-predicted ($\hat\theta$={d['theta']:.1f})")

    ax.set_xscale("log")
    ax.set_xlabel("Mean count per gene")
    ax.set_ylabel("Fraction of zero cells")
    ax.set_title("Gene mean-dropout: simulated vs real")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_total_counts_compare(datasets: list[dict], output_path: Path) -> None:
    """datasets: [{"label", "color", "total_counts"}, ...] (primary first). Density-normalized
    so datasets with very different cell counts are still comparable by shape."""
    positive = [d["total_counts"][d["total_counts"] > 0] for d in datasets]
    all_positive = np.concatenate(positive)
    bins = np.logspace(np.log10(all_positive.min()), np.log10(all_positive.max()), 60)

    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    for d, pos in zip(datasets, positive):
        ax.hist(pos, bins=bins, density=True, color=d["color"], alpha=0.5, label=d["label"])
    ax.set_xscale("log")
    ax.set_xlabel("Total counts per cell")
    ax.set_ylabel("Density (fraction of cells)")
    ax.set_title("Total counts per cell: simulated vs real")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_genes_per_cell_compare(datasets: list[dict], output_path: Path) -> None:
    """datasets: [{"label", "color", "n_genes"}, ...] (primary first)."""
    positive = [d["n_genes"][d["n_genes"] > 0] for d in datasets]
    all_positive = np.concatenate(positive)
    bins = np.logspace(np.log10(all_positive.min()), np.log10(all_positive.max()), 60)

    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    for d, pos in zip(datasets, positive):
        ax.hist(pos, bins=bins, density=True, color=d["color"], alpha=0.5, label=d["label"])
    ax.set_xscale("log")
    ax.set_xlabel("Genes detected per cell")
    ax.set_ylabel("Density (fraction of cells)")
    ax.set_title("Genes detected per cell: simulated vs real")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_sparsity_summary_compare(datasets: list[dict], output_path: Path) -> None:
    """datasets: [{"label", "color", "sparsity"}, ...] (primary first), sparsity is the dict
    returned by compute_sparsity_stats."""
    metrics = [("matrix_zero_frac", "Zero matrix\nentries"), ("empty_row_frac", "Fully-empty\ncells")]
    x = np.arange(len(metrics))
    width = 0.8 / len(datasets)

    fig, ax = plt.subplots(figsize=(6, 5.5))
    for i, d in enumerate(datasets):
        values = [d["sparsity"][key] * 100 for key, _ in metrics]
        offset = (i - (len(datasets) - 1) / 2) * width
        bars = ax.bar(x + offset, values, width=width, color=d["color"], label=d["label"])
        for bar, v in zip(bars, values):
            ax.annotate(f"{v:.1f}%", (bar.get_x() + bar.get_width() / 2, v), ha="center", va="bottom", fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels([label for _, label in metrics])
    ax.set_ylabel("Percent (%)")
    ax.set_ylim(0, 105)
    ax.set_title("Sparsity summary: simulated vs real")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_raw_norm_log_compare(
    datasets: list[dict], target_sum: float, sample_size: int, rng: np.random.Generator, output_path: Path
) -> None:
    """datasets: [{"label", "color", "X"}, ...] (primary first)."""
    staged = []
    for d in datasets:
        tmp = ad.AnnData(X=d["X"].copy())
        raw_vals = sample_values(tmp.X.data.copy(), sample_size, rng)
        sc.pp.normalize_total(tmp, target_sum=target_sum)
        norm_vals = sample_values(tmp.X.data.copy(), sample_size, rng)
        sc.pp.log1p(tmp)
        log_vals = sample_values(tmp.X.data.copy(), sample_size, rng)
        staged.append({"label": d["label"], "color": d["color"], "raw": raw_vals, "norm": norm_vals, "log": log_vals})

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    stage_specs = [
        (axes[0], "raw", "Raw counts", True),
        (axes[1], "norm", f"Normalized (target_sum={target_sum:g})", True),
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
            ax.hist(vals, bins=bins, density=True, color=s["color"], alpha=0.5, label=s["label"])
        ax.set_title(title)
        ax.set_xlabel("Value (nonzero matrix entries)")
    axes[0].set_ylabel("Density")
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle("Raw -> normalized -> log1p: simulated vs real")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def run_compare(args: argparse.Namespace, rng: np.random.Generator) -> None:
    print(f"[load] primary: {args.input}")
    primary = sc.read_h5ad(args.input)
    print(f"[load] primary shape: {primary.n_obs} observations x {primary.n_vars} genes")

    print(f"[load] compare: {args.compare_input}")
    compare = sc.read_h5ad(args.compare_input)
    print(f"[load] compare shape: {compare.n_obs} observations x {compare.n_vars} genes")

    specs = [
        {"label": args.modality, "color": PRIMARY_COLOR, "adata": primary},
        {"label": args.compare_label, "color": COMPARE_COLOR, "adata": compare},
    ]
    for spec in specs:
        adata = spec["adata"]
        raw_counts = adata.layers["counts_pre_batch"] if "counts_pre_batch" in adata.layers else adata.X
        mean, var = gene_mean_var(raw_counts)
        spec["mean"] = mean
        spec["var"] = var
        spec["theta"] = fit_common_dispersion(mean, var)
        spec["zero_frac"] = gene_zero_fraction(raw_counts, adata.n_obs)
        spec["total_counts"] = np.asarray(raw_counts.sum(axis=1)).ravel()
        spec["n_genes"] = genes_per_cell(raw_counts)
        spec["sparsity"] = compute_sparsity_stats(raw_counts, spec["n_genes"])
        print(f"[nb] {spec['label']}: theta_hat = {spec['theta']:.2f}, "
              f"{'counts_pre_batch' if 'counts_pre_batch' in adata.layers else 'X'} used for mean-variance/total-counts")
        print(f"[sparsity] {spec['label']}: {spec['sparsity']['matrix_zero_frac'] * 100:.2f}% zero entries, "
              f"{spec['sparsity']['empty_row_frac'] * 100:.2f}% fully-empty cells")

    plot_mean_variance_compare(specs, args.output_dir / "mean_variance_compare.png")
    plot_mean_dropout_compare(specs, args.output_dir / "mean_dropout_compare.png")
    plot_total_counts_compare(specs, args.output_dir / "total_counts_compare.png")
    plot_genes_per_cell_compare(specs, args.output_dir / "genes_per_cell_compare.png")
    plot_sparsity_summary_compare(specs, args.output_dir / "sparsity_summary_compare.png")
    plot_raw_norm_log_compare(
        [{"label": s["label"], "color": s["color"], "X": s["adata"].X} for s in specs],
        args.target_sum, args.sample_size, rng, args.output_dir / "raw_norm_log_compare.png",
    )
    print(f"[save] Comparison plots written to {args.output_dir}")


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

    raw_counts = adata.layers["counts_pre_batch"] if "counts_pre_batch" in adata.layers else adata.X
    print(f"[nb] Using {'counts_pre_batch' if 'counts_pre_batch' in adata.layers else 'X'} for mean-variance/total-counts")

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
