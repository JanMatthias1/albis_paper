#!/usr/bin/env python
"""
QC diagnostics for a real Xenium slice, matched to the albis cell-level
data for comparison.

Loads a raw Xenium `outs/` bundle (cell_feature_matrix.h5 + cells.csv.gz),
restricts the count matrix to Gene Expression features, computes MAD-based
QC thresholds on total_counts / n_genes_by_counts / count_density, and flags
low-quality cells and control-probe/codeword contamination. Flagged cells are
dropped before saving; no normalization or log-transform is applied.

Expected environment:
    conda activate /dcs04/hicks/data/Jan/sim_project/albis/env/albis-tutorial

Example:
    python sim_paper/code/real_data_qc/xenium_qc.py --slice non_diseased_lung
    python sim_paper/code/real_data_qc/xenium_qc.py --slice lung_cancer
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import median_abs_deviation

SCRIPT_DIR = Path(__file__).resolve().parent
SIM_PAPER_DIR = SCRIPT_DIR.parents[1]

SLICE_RAW_DIRS = {
    "non_diseased_lung": "Xenium_Preview_Human_Non_diseased_Lung_With_Add_on_FFPE",
    "lung_cancer": "Xenium_Preview_Human_Lung_Cancer_With_Add_on_2_FFPE",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="MAD-based QC diagnostics for a raw Xenium slice.")
    parser.add_argument("--slice", choices=sorted(SLICE_RAW_DIRS), required=True)
    parser.add_argument("--input-dir", type=Path, default=None, help="Xenium outs/ directory (overrides --slice default).")
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--nmads", type=float, default=4)
    args = parser.parse_args()

    if args.input_dir is None:
        args.input_dir = SIM_PAPER_DIR / "data" / "real_data" / SLICE_RAW_DIRS[args.slice] / "outs"
    if args.output_dir is None:
        args.output_dir = SIM_PAPER_DIR / "data" / "real_data_qc" / args.slice
    return args


def load_xenium_slice(input_dir: Path) -> sc.AnnData:
    adata = sc.read_10x_h5(input_dir / "cell_feature_matrix.h5", gex_only=False)
    adata.var_names_make_unique()
    adata = adata[:, adata.var["feature_types"] == "Gene Expression"].copy()

    cells = pd.read_csv(input_dir / "cells.csv.gz", index_col="cell_id")
    adata.obs = cells.loc[adata.obs_names]

    sc.pp.calculate_qc_metrics(adata, percent_top=None, inplace=True)
    adata.obs["count_density"] = adata.obs["total_counts"] / adata.obs["cell_area"]
    return adata


def apply_qc_filters(adata, nmads=4):
    # tresholding functions
    # We are currently filtering, and removing cells that are below 3 mad: median_abs_deviation(vals)
    # MAD= median of absolute(vals-median) --> spread around the median, but without being skewed by the extreme outliers

    # we filter based on total_counts, number of detected features and if control probes are present

    thresholds = {}

    # --- Count density ---> we don't filter based on this
    # --> count_density: total amount of transcripts per unit area of the cell
    vals = np.log1p(adata.obs["count_density"])
    med = np.median(vals)
    mad = median_abs_deviation(vals)
    thresholds["count_density"] = med - nmads * mad
    adata.obs["low_count_density"] = vals < thresholds["count_density"]

    # --- Total counts ---
    # --> total_counts: total amount of transcripts detected in the cell
    vals = np.log1p(adata.obs["total_counts"])
    med = np.median(vals)
    mad = median_abs_deviation(vals)
    thresholds["total_counts"] = med - nmads * mad
    adata.obs["low_total_count"] = vals < thresholds["total_counts"]

    # --- Detected features ---
    # --> n_genes_by_counts: number of genes detected in the cell, low number of detected genes may indicate low quality cells
    vals = np.log1p(adata.obs["n_genes_by_counts"])
    med = np.median(vals)
    mad = median_abs_deviation(vals)
    thresholds["detected_features"] = med - nmads * mad
    adata.obs["low_detected_features"] = vals < thresholds["detected_features"]

    # --- Control probes ---
    # --> control_probe_counts: number of counts from control probes, should be zero in real cells, presence indicates that the cell is likely low quality or background
    adata.obs["qc_controls"] = (
        (adata.obs["control_probe_counts"] > 0) |
        (adata.obs["control_codeword_counts"] > 0)
    )

    # --- Combined flag ---
    adata.obs["low_qc"] = (
        #adata.obs["low_count_density"]
        # low count density is to stringent, we are only excluding cells with 0 count density
        (adata.obs["count_density"] == 0)|
        adata.obs["low_total_count"] |
        adata.obs["low_detected_features"] |
        adata.obs["qc_controls"]
    )

    return adata, thresholds


def plot_qc_exclusions(adata, x_col="x_centroid", y_col="y_centroid", s=1, sample_name=None, output_path=None):
    """
    Plot Xenium QC exclusions by metric without modifying adata.obs.
    Expects boolean QC flags already stored in adata.obs:
    - low_count_density
    - low_total_count
    - low_detected_features
    - qc_controls
    """

    # Assign QC reason (first failure wins, else 'kept')
    qc_labels = []
    for i in range(adata.n_obs):
        obs_i = adata.obs.iloc[i]
        if obs_i["count_density"] == 0:
            qc_labels.append("zero_count_density")
        elif obs_i["low_total_count"]:
            qc_labels.append("low_total_count")
        elif obs_i["low_detected_features"]:
            qc_labels.append("low_detected_features")
        elif obs_i["qc_controls"]:
            qc_labels.append("control_probe")
        else:
            qc_labels.append("kept")

    # Color mapping
    color_map = {
        "kept": "lightgrey",
        "zero_count_density": "blue",
        "low_total_count": "green",
        "low_detected_features": "orange",
        "control_probe": "red",
    }
    colors = [color_map[label] for label in qc_labels]

    # Plot scatter
    plt.figure(figsize=(7, 7))
    plt.scatter(
        adata.obs[x_col], adata.obs[y_col],
        c=colors, s=s
    )
    plt.gca().invert_yaxis()
    plt.title(f"QC-excluded cells by metric {sample_name}")

    # Legend
    legend_elements = [
        mpatches.Patch(color=color, label=label)
        for label, color in color_map.items()
    ]
    plt.legend(handles=legend_elements, markerscale=6, loc="upper right")
    plt.savefig(output_path, dpi=180)
    plt.close()


def plot_qc_metrics(adata, thresholds, sample_name=None, output_path=None):
    """
    Plot QC histograms for one AnnData object.
    """
    fig, axs = plt.subplots(2, 3, figsize=(15, 7.5))

    # total counts
    axs[0,0].hist(np.log1p(adata.obs["total_counts"]), bins=100, color="black")
    axs[0,0].set_title("log(total_counts)")
    axs[0,0].axvline(thresholds["total_counts"], color="red", linestyle="--")

    # cell area
    axs[0,1].hist(np.log1p(adata.obs["cell_area"]), bins=100, color="black")
    axs[0,1].set_title("log(cell_area)")

    # counts per area
    axs[0, 2].hist(np.log1p(adata.obs["count_density"]), bins=100, color="black")
    axs[0, 2].axvline(thresholds["count_density"], color="red", linestyle="--", label="Low QC cutoff")
    axs[0, 2].axvline(0, color="blue", linestyle=":", label="count_density = 0")
    axs[0, 2].set_title("log(total_counts / cell_area)")
    axs[0, 2].legend()

    # n genes by counts
    axs[1,0].hist(np.log1p(adata.obs["n_genes_by_counts"]), bins=100, color="black")
    axs[1,0].set_title("log(n_genes_by_counts)")
    axs[1,0].axvline(thresholds["detected_features"], color="red", linestyle="--")

    # total counts, raw scale
    axs[1,1].hist(adata.obs["total_counts"], bins=100, color="black")
    axs[1,1].set_title("total_counts (raw)")
    axs[1,1].axvline(np.expm1(thresholds["total_counts"]), color="red", linestyle="--")

    # n genes by counts, raw scale
    axs[1,2].hist(adata.obs["n_genes_by_counts"], bins=100, color="black")
    axs[1,2].set_title("n_genes_by_counts (raw)")
    axs[1,2].axvline(np.expm1(thresholds["detected_features"]), color="red", linestyle="--")

    # optional title
    if sample_name:
        fig.suptitle(f"QC metrics for {sample_name}", fontsize=16)

    plt.tight_layout()
    plt.savefig(output_path, dpi=180)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    if not args.input_dir.is_dir():
        raise SystemExit(f"Input dir not found: {args.input_dir}")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f"[load] {args.input_dir}")
    adata = load_xenium_slice(args.input_dir)
    print(f"[load] AnnData shape: {adata.n_obs} cells x {adata.n_vars} genes")

    print(f"[qc] Applying MAD-based QC filters (nmads={args.nmads})")
    adata, thresholds = apply_qc_filters(adata, nmads=args.nmads)
    print(f"[qc] thresholds (log1p scale): {thresholds}")
    print(f"[qc] low_qc cells: {int(adata.obs['low_qc'].sum())} / {adata.n_obs}")

    plot_qc_metrics(adata, thresholds, sample_name=args.slice, output_path=args.output_dir / "qc_metrics.png")
    plot_qc_exclusions(adata, sample_name=args.slice, output_path=args.output_dir / "qc_exclusions.png")

    n_cells_before = adata.n_obs
    summary = {
        "slice": args.slice,
        "input_dir": str(args.input_dir),
        "nmads": args.nmads,
        "n_cells_before_qc": int(n_cells_before),
        "n_genes": int(adata.n_vars),
        "thresholds_log1p": {k: float(v) for k, v in thresholds.items()},
        "n_low_qc": int(adata.obs["low_qc"].sum()),
        "n_zero_count_density": int((adata.obs["count_density"] == 0).sum()),
        "n_low_total_count": int(adata.obs["low_total_count"].sum()),
        "n_low_detected_features": int(adata.obs["low_detected_features"].sum()),
        "n_qc_controls": int(adata.obs["qc_controls"].sum()),
    }

    adata = adata[~adata.obs["low_qc"]].copy()
    print(f"[qc] Dropped {n_cells_before - adata.n_obs} low_qc cells; {adata.n_obs} remain")
    summary["n_cells_after_qc"] = int(adata.n_obs)

    with open(args.output_dir / "qc_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[save] Summary written to {args.output_dir / 'qc_summary.json'}")

    adata.write_h5ad(args.output_dir / f"{args.slice}_qc.h5ad")
    print(f"[save] QC-filtered AnnData written to {args.output_dir / f'{args.slice}_qc.h5ad'}")

    print(f"[save] Plots written to {args.output_dir}")


if __name__ == "__main__":
    main()
