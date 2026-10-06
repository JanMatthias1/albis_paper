#!/usr/bin/env python
"""
QC diagnostics for a real (non-HD) Visium sample using SpotSweeper's
local-outlier detection, matched to the albis spot-level data for
comparison.

Loads a raw Visium sample bundle -- this download ships spatial/ as a
separate <prefix>_spatial.tar.gz alongside a top-level
<prefix>_filtered_feature_bc_matrix.h5, rather than a single self-contained
sample directory -- computes total_counts / n_genes_by_counts /
pct_counts_mt, and flags spatially-local outliers on each metric via
spotsweeper.local_outliers. Flagged spots are dropped before saving; no
normalization or log-transform is applied.

Expected environment:
    conda activate /dcs04/hicks/data/Jan/sim_project/albis/env/albis-tutorial
    pip install spotsweeper   # one-time; not yet in requirements

Example:
    python sim_paper/code/real_data_qc/visium_qc.py --sample breast_cancer
"""

from __future__ import annotations

import argparse
import json
import tarfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import scanpy as sc
import spotsweeper.local_outliers as lo
import spotsweeper.plot_QC as plot_QC

SCRIPT_DIR = Path(__file__).resolve().parent
SIM_PAPER_DIR = SCRIPT_DIR.parents[1]

SAMPLE_RAW_DIRS = {
    "breast_cancer": "breast_cancer_visium",
    "lymph_node": "lymph_node_visium",
    "tonsil": "tonsil_visium",
}

# metric -> (direction, log, renamed outlier column)
METRICS = {
    "total_counts": ("lower", True, "total_counts_outliers"),
    "n_genes_by_counts": ("lower", True, "n_genes_outliers"),
    "pct_counts_mt": ("higher", False, "mt_outliers"),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="SpotSweeper local-outlier QC for a raw Visium sample.")
    parser.add_argument("--sample", choices=sorted(SAMPLE_RAW_DIRS), required=True)
    parser.add_argument("--input-dir", type=Path, default=None, help="Raw 10x Visium download dir (overrides --sample default).")
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--n-neighbors", type=int, default=36)
    parser.add_argument("--cutoff", type=float, default=3.0, help="Robust z-score cutoff for local_outliers.")
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()

    if args.input_dir is None:
        args.input_dir = SIM_PAPER_DIR / "data" / "real_data" / SAMPLE_RAW_DIRS[args.sample]
    if args.output_dir is None:
        args.output_dir = SIM_PAPER_DIR / "data" / "real_data_qc" / SAMPLE_RAW_DIRS[args.sample]
    return args


def ensure_spatial_extracted(input_dir: Path) -> None:
    """This download ships spatial/ as a separate <prefix>_spatial.tar.gz rather than an
    already-extracted spatial/ folder; sc.read_visium needs the latter, so extract it in
    place (idempotent -- skipped once spatial/ exists)."""
    if (input_dir / "spatial").is_dir():
        return
    spatial_tars = list(input_dir.glob("*_spatial.tar.gz"))
    if not spatial_tars:
        raise SystemExit(f"No spatial/ dir and no *_spatial.tar.gz found in {input_dir}")
    print(f"[extract] {spatial_tars[0].name} -> {input_dir / 'spatial'}")
    with tarfile.open(spatial_tars[0], "r:gz") as tar:
        tar.extractall(input_dir)


def load_visium_sample(input_dir: Path) -> sc.AnnData:
    ensure_spatial_extracted(input_dir)

    count_files = list(input_dir.glob("*filtered_feature_bc_matrix.h5"))
    if len(count_files) != 1:
        raise SystemExit(f"Expected exactly one *filtered_feature_bc_matrix.h5 in {input_dir}, found {len(count_files)}")

    adata = sc.read_visium(input_dir, count_file=count_files[0].name)
    adata.var_names_make_unique()
    adata = adata[adata.obs["in_tissue"] == 1].copy()

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
    """One SpotSweeper spatial diagnostic per metric (value + local-outlier ring)."""
    for metric, (_direction, _log, renamed) in METRICS.items():
        # plot_qc_metrics returns the plt module (operating on its current figure), not a
        # Figure instance -- save/close via plt, not fig.savefig()/fig.close().
        plot_QC.plot_qc_metrics(
            adata, sample_id="region", sample="sample1", metric=metric, outliers=renamed,
            coord_key="spatial", legend=True,
            title=f"{sample_name}: {metric} (SpotSweeper local outliers)",
        )
        plt.savefig(output_dir / f"qc_{metric}.png", dpi=500)
        plt.close()


def plot_qc_exclusions(adata, sample_name, output_path: Path) -> None:
    """Combined exclusion map: which criterion (if any) flagged each spot. First failure wins."""
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

    spatial = adata.obsm["spatial"]
    plt.figure(figsize=(7, 7))
    plt.scatter(spatial[:, 0], spatial[:, 1], c=colors, s=4)
    plt.gca().invert_yaxis()
    plt.gca().set_aspect("equal")
    plt.title(f"QC-excluded spots by metric: {sample_name}")

    legend_elements = [mpatches.Patch(color=color, label=label) for label, color in color_map.items()]
    plt.legend(handles=legend_elements, markerscale=6, loc="upper right")
    plt.savefig(output_path, dpi=500)
    plt.close()


def main() -> None:
    args = parse_args()
    if not args.input_dir.is_dir():
        raise SystemExit(f"Input dir not found: {args.input_dir}")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f"[load] {args.input_dir}")
    adata = load_visium_sample(args.input_dir)
    print(f"[load] AnnData shape: {adata.n_obs} spots x {adata.n_vars} genes")

    print(f"[qc] Running SpotSweeper local_outliers (n_neighbors={args.n_neighbors}, cutoff={args.cutoff})")
    adata = apply_qc_filters(adata, n_neighbors=args.n_neighbors, cutoff=args.cutoff, workers=args.workers)
    print(f"[qc] low_qc spots: {int(adata.obs['low_qc'].sum())} / {adata.n_obs}")

    print("[plot] Per-metric spatial diagnostics")
    plot_qc_metrics_spatial(adata, args.sample, args.output_dir)
    plot_qc_exclusions(adata, args.sample, args.output_dir / "qc_exclusions.png")

    n_spots_before = adata.n_obs
    summary = {
        "sample": args.sample,
        "input_dir": str(args.input_dir),
        "method": "spotsweeper.local_outliers",
        "n_neighbors": args.n_neighbors,
        "cutoff": args.cutoff,
        "n_spots_before_qc": int(n_spots_before),
        "n_genes": int(adata.n_vars),
        "n_low_qc": int(adata.obs["low_qc"].sum()),
        "n_low_total_count": int(adata.obs["total_counts_outliers"].sum()),
        "n_low_detected_features": int(adata.obs["n_genes_outliers"].sum()),
        "n_high_mt": int(adata.obs["mt_outliers"].sum()),
    }

    adata = adata[~adata.obs["low_qc"]].copy()
    print(f"[qc] Dropped {n_spots_before - adata.n_obs} low_qc spots; {adata.n_obs} remain")
    summary["n_spots_after_qc"] = int(adata.n_obs)

    with open(args.output_dir / "qc_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[save] Summary written to {args.output_dir / 'qc_summary.json'}")

    out_h5ad = args.output_dir / f"{args.sample}_visium_qc.h5ad"
    adata.write_h5ad(out_h5ad)
    print(f"[save] QC-filtered AnnData written to {out_h5ad}")

    print(f"[save] Plots written to {args.output_dir}")


if __name__ == "__main__":
    main()
