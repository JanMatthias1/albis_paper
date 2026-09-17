#!/usr/bin/env python
"""
Sanity check for the plain cell-type pipeline: does Harmony's large (~37-44%)
embedding shift actually recover real structure, or is the good post-Harmony
ARI an artifact of how far Harmony moves things? Runs the exact same
resolution-matched Leiden/ARI procedure as ari_vs_ground_truth.py, but on
X_pca_pre_harmony (same batch-corrupted PCA, before Harmony) instead of
X_pca_harmony -- direct pre/post comparison on identical data.

Usage:
    conda activate /dcs04/hicks/data/Jan/sim_project/sim_paper/env/albis-tutorial
    python check_pre_harmony_ari.py --modality cell
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import scanpy as sc
from sklearn.metrics import adjusted_rand_score

CLUSTERING_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CLUSTERING_DIR))
from ari_vs_ground_truth import find_resolution_for_k  # noqa: E402

SIM_PAPER_DIR = CLUSTERING_DIR.parents[1]

INPUTS = {
    "cell": SIM_PAPER_DIR / "data/figure_3/pca_harmony_single_cell/cell/pca_harmony_qc/simulation_cell_z_pca_harmony_qc.h5ad",
    "bin": SIM_PAPER_DIR / "data/figure_3/pca_harmony_single_cell/bin/pca_harmony_qc/simulation_bin_z_pca_harmony_qc.h5ad",
    "bin16um": SIM_PAPER_DIR / "data/figure_3/pca_harmony_single_cell/bin16um/pca_harmony_qc/simulation_bin_z_pca_harmony_qc.h5ad",
    "spot": SIM_PAPER_DIR / "data/figure_3/pca_harmony_single_cell/spot/pca_harmony_qc/simulation_spot_z_pca_harmony_qc.h5ad",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--modality", required=True, choices=list(INPUTS))
    p.add_argument("--n-pcs", type=int, default=30)
    p.add_argument("--n-neighbors", type=int, default=15)
    p.add_argument("--random-state", type=int, default=0)
    p.add_argument("--res-lo", type=float, default=0.01)
    p.add_argument("--res-hi", type=float, default=3.0)
    p.add_argument("--max-iter", type=int, default=15)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    path = INPUTS[args.modality]
    print(f"[load] {path}")
    adata = sc.read_h5ad(path)
    print(f"[load] {adata.n_obs:,} x {adata.n_vars:,}")

    n_pcs = min(args.n_pcs, adata.obsm["X_pca_pre_harmony"].shape[1])
    rep_key = "X_pca_pre_harmony_ari"
    adata.obsm[rep_key] = np.asarray(adata.obsm["X_pca_pre_harmony"])[:, :n_pcs].copy()
    print(f"[neighbors] use_rep=X_pca_pre_harmony (NO Harmony) n_pcs={n_pcs} n_neighbors={args.n_neighbors}")
    sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, use_rep=rep_key, random_state=args.random_state)

    for truth_col in ("cell_type_true", "domain_true"):
        if truth_col not in adata.obs:
            continue
        target_k = int(adata.obs[truth_col].nunique())
        cluster_key = f"leiden_pre_harmony_{truth_col}"
        best_res, best_k, tried = find_resolution_for_k(
            adata, target_k, args.res_lo, args.res_hi, args.max_iter, args.random_state, cluster_key,
        )
        ari = adjusted_rand_score(adata.obs[truth_col], adata.obs[cluster_key])
        leak = adjusted_rand_score(adata.obs["slice_id"], adata.obs[cluster_key]) if "slice_id" in adata.obs else float("nan")
        print(
            f"[RESULT] {args.modality} / {truth_col} (PRE-Harmony): "
            f"resolution={best_res:.4f} achieved_k={best_k} (target {target_k})  "
            f"ARI={ari:.4f}  slice_leak_ARI={leak:.4f}"
        )


if __name__ == "__main__":
    main()
