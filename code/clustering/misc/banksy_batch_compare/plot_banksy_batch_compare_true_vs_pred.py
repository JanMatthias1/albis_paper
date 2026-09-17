#!/usr/bin/env python
"""
After-Harmony "ground truth vs. predicted cluster" UMAP + contingency heatmap
for the BANKSY batch-compare job (35608415) -- the qualitative panel type
that's already generated for the plain PCA+Harmony cell-type pipeline
(leiden_pca_qc_celltype_matched/plots/umap_true_vs_predicted_*.png,
contingency_*.png) but was never run for this job: banksy_batch_compare.sh
only ran 01_build_banksy_matrix.py -> 02_leiden_resolution_sweep.py ->
composition_recovery.py, no step03_cluster_and_plot.py step.

02_leiden_resolution_sweep.py's output h5ad already carries the resolution-matched
predicted labels for BOTH ground truths (leiden_domain_true, matched to
domain_true's k; leiden_cell_type_true, matched to cell_type_true's k) plus
X_pca_harmony -- so this reuses that instead of re-clustering, and reuses
step03_cluster_and_plot.py's own plotting functions for a pixel-matched
panel (same PANEL_FIGSIZE/MARGINS/DPI, same true-vs-predicted + contingency
layout).

Usage:
    conda activate /dcs04/hicks/data/Jan/sim_project/sim_paper/env/albis-tutorial
    python plot_banksy_batch_compare_true_vs_pred.py --modality cell --batch tuned
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import scanpy as sc

CLUSTERING_DIR = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(CLUSTERING_DIR))
# 01.2_pca_harmony.py's numeric-leading name isn't a valid Python identifier,
# so it can't be `from ... import`-ed directly -- load it by file path instead.
_pca_harmony_spec = importlib.util.spec_from_file_location(
    "pca_harmony", CLUSTERING_DIR / "01.2_pca_harmony.py"
)
_pca_harmony = importlib.util.module_from_spec(_pca_harmony_spec)
_pca_harmony_spec.loader.exec_module(_pca_harmony)
compute_umap = _pca_harmony.compute_umap
from step03_cluster_and_plot import plot_umap_true_vs_predicted, plot_contingency_heatmap  # noqa: E402
from pc_pairs import sampled_indices  # noqa: E402

SIM_PAPER_DIR = CLUSTERING_DIR.parents[1]
ROOT = SIM_PAPER_DIR / "data" / "figure_3" / "banksy_batch_compare"

# (ground-truth column, resolution-matched predicted-cluster column)
GT_PRED_PAIRS = [("domain_true", "leiden_domain_true"), ("cell_type_true", "leiden_cell_type_true")]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--modality", required=True, choices=["cell", "bin", "spot"])
    p.add_argument("--batch", default="prebatch", choices=["prebatch", "tuned"], help="Ignored if --input is given (only used for the default path/label).")
    p.add_argument("--input", type=Path, default=None, help="Override the default banksy_batch_compare ari_recovery.h5ad path (e.g. for a batch_sigma-slide run).")
    p.add_argument("--label", default=None, help="Label used in printouts when --input is given (e.g. 'bs0.5').")
    p.add_argument("--n-neighbors", type=int, default=15)
    p.add_argument("--random-state", type=int, default=0)
    p.add_argument("--umap-max-obs", type=int, default=50000)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    if args.input is not None:
        h5ad_path = args.input
        plot_dir = args.input.parent / "plots"
    else:
        run_dir = ROOT / args.modality / args.batch
        h5ad_path = run_dir / "ari" / f"simulation_{args.modality}_z_ari_recovery.h5ad"
        plot_dir = run_dir / "ari" / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    print(f"[load] {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    print(f"[load] {adata.n_obs:,} x {adata.n_vars:,}")

    obs_idx = sampled_indices(adata.n_obs, args.umap_max_obs, False, args.random_state)
    if len(obs_idx) < adata.n_obs:
        print(f"[sample] {len(obs_idx):,} / {adata.n_obs:,} observations")
    adata = adata[obs_idx].copy() if len(obs_idx) < adata.n_obs else adata

    print("[umap] after Harmony (X_pca_harmony)")
    compute_umap(adata, "X_pca_harmony", "X_umap", args.n_neighbors, args.random_state)

    for true_key, pred_key in GT_PRED_PAIRS:
        if true_key not in adata.obs or pred_key not in adata.obs:
            print(f"[skip] missing {true_key} or {pred_key}")
            continue
        out1 = plot_dir / f"umap_true_vs_predicted_{true_key}.png"
        out2 = plot_dir / f"contingency_{true_key}.png"
        plot_umap_true_vs_predicted(adata, true_key, pred_key, out1)
        plot_contingency_heatmap(adata, true_key, pred_key, out2)
        print(f"[save] {out1}")
        print(f"[save] {out2}")


if __name__ == "__main__":
    main()
