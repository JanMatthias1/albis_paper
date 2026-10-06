"""Diagnostic plots for a simulated dataset (01-06 under plots/<stem>/).

generate_simulation_noisy.py calls save_diagnostic_plots() right before it
saves the h5ad. Run this file directly to redraw the plots from an existing
h5ad without regenerating the data:

    python simulation_diagnostic_plots.py --input <run>/<stem>.h5ad --modality bin

Plots go to <h5ad dir>/plots/<stem>/ (same place the generator writes them),
as PNG plus PDF/SVG via manuscript_style.save_figure.
"""

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless-safe backend, must be set before pyplot/albis plotting is used
import matplotlib.pyplot as plt
import numpy as np

import albis as ab

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from manuscript_style import save_figure

MAX_PLOT_OBS = 100_000


def save_diagnostic_plots(adata, plots_dir, output_modality):
    plots_dir = Path(plots_dir)
    plots_dir.mkdir(parents=True, exist_ok=True)

    def savefig(fig, name):
        path = plots_dir / name
        save_figure(fig, path, dpi=500, bbox_inches="tight")
        plt.close(fig)
        print(f"  saved {path}")

    print("\nGenerating diagnostic plots...")

    if adata.n_obs > MAX_PLOT_OBS:
        rng = np.random.default_rng(2025)
        plot_idx = np.sort(rng.choice(adata.n_obs, size=MAX_PLOT_OBS, replace=False))
        plot_adata = adata[plot_idx].copy()
        print(f"Plotting a deterministic subset of {MAX_PLOT_OBS:,} / {adata.n_obs:,} observations.")
    else:
        plot_adata = adata

    # 1. Aligned spatial coordinates colored by ground-truth domain
    fig = ab.plot(plot_adata, view="2d", coordinates="aligned", color="domain_true", point_size=4)
    savefig(fig, "01_aligned_domain_true.png")

    # 2. Aligned spatial coordinates colored by ground-truth cell type
    fig = ab.plot(plot_adata, view="2d", coordinates="aligned", color="cell_type_true", point_size=4)
    savefig(fig, "02_aligned_cell_type_true.png")

    # 3. Unaligned (per-slice, pre-registration) coordinates colored by slice id
    fig = ab.plot(plot_adata, view="2d", coordinates="unaligned", color="slice_id", point_size=4)
    savefig(fig, "03_unaligned_slice_id.png")

    # 4. QC: total counts / genes detected per observation
    X = adata.X
    total_counts = np.asarray(X.sum(axis=1)).ravel()
    genes_detected = np.asarray((X > 0).sum(axis=1)).ravel()

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].hist(total_counts, bins=50, color="steelblue")
    axes[0].set_xlabel(f"Total counts per {output_modality}")
    axes[0].set_ylabel(f"Number of {output_modality}s")
    axes[0].set_title("Total counts distribution")

    axes[1].hist(genes_detected, bins=50, color="darkorange")
    axes[1].set_xlabel(f"Genes detected per {output_modality}")
    axes[1].set_ylabel(f"Number of {output_modality}s")
    axes[1].set_title("Genes detected distribution")

    fig.tight_layout()
    savefig(fig, "04_qc_count_distributions.png")

    # 5. Number of observations per slice
    slice_counts = adata.obs["slice_id"].value_counts().sort_index()
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(slice_counts.index.astype(str), slice_counts.values, color="slategray")
    ax.set_xlabel("Slice ID")
    ax.set_ylabel(f"Number of {output_modality}s")
    ax.set_title(f"{output_modality.capitalize()}s per slice")
    fig.tight_layout()
    savefig(fig, f"05_{output_modality}s_per_slice.png")

    # 6. Domain composition per slice (stacked bar)
    comp = adata.obs.groupby(["slice_id", "domain_true"]).size().unstack(fill_value=0)
    comp_frac = comp.div(comp.sum(axis=1), axis=0)
    fig, ax = plt.subplots(figsize=(8, 5))
    bottom = np.zeros(len(comp_frac))
    for domain in comp_frac.columns:
        ax.bar(comp_frac.index.astype(str), comp_frac[domain], bottom=bottom, label=str(domain))
        bottom += comp_frac[domain].values
    ax.set_xlabel("Slice ID")
    ax.set_ylabel(f"Fraction of {output_modality}s")
    ax.set_title("Domain composition per slice")
    ax.legend(title="domain_true", bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    savefig(fig, "06_domain_composition_per_slice.png")

    print(f"\nAll diagnostic plots saved under {plots_dir}")


def main():
    parser = argparse.ArgumentParser(description="Redraw the 01-06 diagnostic plots from a saved simulation h5ad.")
    parser.add_argument("--input", type=Path, required=True, help="raw (pre-QC) simulation h5ad")
    parser.add_argument("--modality", required=True, help="output modality word used in labels/file names (cell, bin, spot)")
    parser.add_argument("--plots-dir", type=Path, default=None, help="default: <input dir>/plots/<input stem>")
    args = parser.parse_args()
    import anndata as ad
    adata = ad.read_h5ad(args.input)
    plots_dir = args.plots_dir or args.input.parent / "plots" / args.input.stem
    save_diagnostic_plots(adata, plots_dir, args.modality)


if __name__ == "__main__":
    main()
