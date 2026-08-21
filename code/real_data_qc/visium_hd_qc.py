#!/usr/bin/env python
"""
QC diagnostics for a real Visium HD sample using SpotSweeper's local-outlier
detection, matched to the sim_app bin-level data for comparison.

Loads a raw Visium HD `binned_outputs/square_XXXum/` bundle
(filtered_feature_bc_matrix.h5 + spatial/tissue_positions_list.csv), computes
total_counts / n_genes_by_counts / pct_counts_mt, and flags spatially-local
outliers on each metric via spotsweeper.local_outliers -- without
normalizing, log-transforming, or dropping any bins.

Raw downloads must already be extracted via extract_breast_visium_hd.py.

Expected environment:
    conda activate /dcs04/hicks/data/Jan/sim_project/sim_app/env/sim-app-tutorial
    pip install spotsweeper   # one-time; not yet in requirements

Example:
    python sim_paper/code/real_data_qc/visium_hd_qc.py --sample breast_cancer
    python sim_paper/code/real_data_qc/visium_hd_qc.py --sample human_pancreas
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import spotsweeper.local_outliers as lo
import spotsweeper.plot_QC as plot_QC

SCRIPT_DIR = Path(__file__).resolve().parent
SIM_PAPER_DIR = SCRIPT_DIR.parents[1]

SAMPLE_RAW_DIRS = {
    "breast_cancer": "breast_cancer/Visium_HD_11mm_Human_Breast_Cancer",
    "human_pancreas": "human_pancreas_HD/Visium_HD_11mm_Human_Pancreas",
}

POSITION_COLS = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]

# metric -> (direction, log, renamed outlier column)
METRICS = {
    "total_counts": ("lower", True, "total_counts_outliers"),
    "n_genes_by_counts": ("lower", True, "n_genes_outliers"),
    "pct_counts_mt": ("higher", False, "mt_outliers"),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="SpotSweeper local-outlier QC for a raw Visium HD sample.")
    parser.add_argument("--sample", choices=sorted(SAMPLE_RAW_DIRS), required=True)
    parser.add_argument("--bin-size-um", type=int, default=8, help="Matches sim_app's bin_size_um.")
    parser.add_argument("--input-dir", type=Path, default=None, help="binned_outputs/square_XXXum dir (overrides --sample/--bin-size-um default).")
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--n-neighbors", type=int, default=36)
    parser.add_argument("--cutoff", type=float, default=3.0, help="Robust z-score cutoff for local_outliers.")
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()

    if args.input_dir is None:
        args.input_dir = (
            SIM_PAPER_DIR / "data" / "real_data" / SAMPLE_RAW_DIRS[args.sample]
            / "binned_outputs" / f"square_{args.bin_size_um:03d}um"
        )
    if args.output_dir is None:
        args.output_dir = SIM_PAPER_DIR / "data" / "real_data_qc" / f"{args.sample}_visium_hd"
    return args


def load_tissue_positions(spatial_dir: Path) -> pd.DataFrame:
    """Handles both the older, headerless tissue_positions_list.csv and the newer
    tissue_positions.parquet (10x has shipped either depending on download vintage)."""
    csv_path = spatial_dir / "tissue_positions_list.csv"
    if csv_path.exists():
        return pd.read_csv(csv_path, header=None, names=POSITION_COLS, index_col="barcode")
    parquet_path = spatial_dir / "tissue_positions.parquet"
    return pd.read_parquet(parquet_path).set_index("barcode")


def load_visium_hd_sample(input_dir: Path) -> sc.AnnData:
    adata = sc.read_10x_h5(input_dir / "filtered_feature_bc_matrix.h5")
    adata.var_names_make_unique()

    positions = load_tissue_positions(input_dir / "spatial")
    adata.obs = positions.loc[adata.obs_names]
    adata = adata[adata.obs["in_tissue"] == 1].copy()
    adata.obsm["spatial"] = adata.obs[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy()

    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, inplace=True)
    return adata


def apply_qc_filters(adata, n_neighbors=36, cutoff=3.0, workers=4):
    adata.obs["region"] = "sample1"  # single-sample dataset; spotsweeper groups outlier detection by this key

    for metric, (direction, log, renamed) in METRICS.items():
        lo.local_outliers(
            adata, metric=metric, sample_key="region", direction=direction, log=log,
            n_neighbors=n_neighbors, cutoff=cutoff, workers=workers,
        )
        adata.obs.rename(columns={f"{metric}_outliers": renamed}, inplace=True)

    adata.obs["low_qc"] = (
        adata.obs["total_counts_outliers"] |
        adata.obs["n_genes_outliers"] |
        adata.obs["mt_outliers"]
    )
    return adata


def plot_qc_metrics_spatial(adata, sample_name, output_dir: Path) -> None:
    """One SpotSweeper spatial diagnostic per metric (value + local-outlier ring), ring_overlay=False
    per spotsweeper's own guidance for Visium HD's bin density."""
    for metric, (_direction, _log, renamed) in METRICS.items():
        fig = plot_QC.plot_qc_metrics(
            adata, sample_id="region", sample="sample1", metric=metric, outliers=renamed,
            coord_key="spatial", ring_overlay=False, legend=True,
            title=f"{sample_name}: {metric} (SpotSweeper local outliers)",
        )
        fig.savefig(output_dir / f"qc_{metric}.png", dpi=180)
        fig.close()


def plot_qc_exclusions(adata, sample_name, output_path: Path) -> None:
    """Combined exclusion map: which criterion (if any) flagged each bin. First failure wins."""
    qc_labels = []
    for i in range(adata.n_obs):
        obs_i = adata.obs.iloc[i]
        if obs_i["total_counts_outliers"]:
            qc_labels.append("low_total_count")
        elif obs_i["n_genes_outliers"]:
            qc_labels.append("low_detected_features")
        elif obs_i["mt_outliers"]:
            qc_labels.append("high_mt")
        else:
            qc_labels.append("kept")

    color_map = {
        "kept": "lightgrey",
        "low_total_count": "green",
        "low_detected_features": "orange",
        "high_mt": "red",
    }
    colors = [color_map[label] for label in qc_labels]

    plt.figure(figsize=(7, 7))
    plt.scatter(adata.obs["array_col"], adata.obs["array_row"], c=colors, s=1)
    plt.gca().invert_yaxis()
    plt.title(f"QC-excluded bins by metric: {sample_name}")

    legend_elements = [mpatches.Patch(color=color, label=label) for label, color in color_map.items()]
    plt.legend(handles=legend_elements, markerscale=6, loc="upper right")
    plt.savefig(output_path, dpi=180)
    plt.close()


def main() -> None:
    args = parse_args()
    if not args.input_dir.is_dir():
        raise SystemExit(f"Input dir not found: {args.input_dir}")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f"[load] {args.input_dir}")
    adata = load_visium_hd_sample(args.input_dir)
    print(f"[load] AnnData shape: {adata.n_obs} bins x {adata.n_vars} genes")

    print(f"[qc] Running SpotSweeper local_outliers (n_neighbors={args.n_neighbors}, cutoff={args.cutoff})")
    adata = apply_qc_filters(adata, n_neighbors=args.n_neighbors, cutoff=args.cutoff, workers=args.workers)
    print(f"[qc] low_qc bins: {int(adata.obs['low_qc'].sum())} / {adata.n_obs}")

    print("[plot] Per-metric spatial diagnostics")
    plot_qc_metrics_spatial(adata, args.sample, args.output_dir)
    plot_qc_exclusions(adata, args.sample, args.output_dir / "qc_exclusions.png")

    summary = {
        "sample": args.sample,
        "input_dir": str(args.input_dir),
        "bin_size_um": args.bin_size_um,
        "method": "spotsweeper.local_outliers",
        "n_neighbors": args.n_neighbors,
        "cutoff": args.cutoff,
        "n_bins": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_low_qc": int(adata.obs["low_qc"].sum()),
        "n_low_total_count": int(adata.obs["total_counts_outliers"].sum()),
        "n_low_detected_features": int(adata.obs["n_genes_outliers"].sum()),
        "n_high_mt": int(adata.obs["mt_outliers"].sum()),
    }
    with open(args.output_dir / "qc_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[save] Summary written to {args.output_dir / 'qc_summary.json'}")

    out_h5ad = args.output_dir / f"{args.sample}_visium_hd_qc.h5ad"
    adata.write_h5ad(out_h5ad)
    print(f"[save] AnnData with QC flags written to {out_h5ad}")

    print(f"[save] Plots written to {args.output_dir}")


if __name__ == "__main__":
    main()
