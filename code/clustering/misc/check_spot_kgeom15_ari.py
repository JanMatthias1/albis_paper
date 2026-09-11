#!/usr/bin/env python
"""
Spot, prebatch, at spot's OWN k_geom=15 (not the shared cell-tuned k_geom=200
used by banksy_batch_compare.sh) -- pre- vs post-Harmony ARI, to check
whether the earlier "Harmony rescues spot from 0.021 to 0.678" result was
actually a k_geom-oversmoothing artifact rather than a real batch/Harmony
effect (k_geom=200 was ~20% of an average spot slice's population, i.e. not
local; k_geom=15 is ~1.5%, matching cell/bin16um's locality).
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

N_PCS = 30
N_NEIGHBORS = 15
RANDOM_STATE = 0


def run(adata, rep_key: str, label: str) -> None:
    ari_rep = f"{rep_key}_ari"
    adata.obsm[ari_rep] = np.asarray(adata.obsm[rep_key])[:, :N_PCS].copy()
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, use_rep=ari_rep, random_state=RANDOM_STATE)
    for truth_col in ("domain_true", "cell_type_true"):
        target_k = int(adata.obs[truth_col].nunique())
        cluster_key = f"leiden_{label}_{truth_col}"
        best_res, best_k, _ = find_resolution_for_k(
            adata, target_k, 0.01, 3.0, 15, RANDOM_STATE, cluster_key,
        )
        ari = adjusted_rand_score(adata.obs[truth_col], adata.obs[cluster_key])
        leak = adjusted_rand_score(adata.obs["slice_id"], adata.obs[cluster_key])
        print(f"[RESULT] spot k_geom=15 {label:14s} / {truth_col:15s}: k={best_k} (target {target_k})  ARI={ari:.4f}  slice_leak_ARI={leak:.4f}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--condition", required=True, choices=["prebatch", "tuned"])
    args = parser.parse_args()
    path = SIM_PAPER_DIR / f"data/figure_3/banksy_kgeom_check/spot_{args.condition}_kgeom15/banksy_matrix/simulation_spot_z_banksy_pca_harmony_qc.h5ad"

    print(f"[load] {path}")
    adata = sc.read_h5ad(path)
    print(f"[load] {adata.n_obs:,} x {adata.n_vars:,}")

    pre = np.asarray(adata.obsm["X_pca_pre_harmony"])
    post = np.asarray(adata.obsm["X_pca_harmony"])
    shift = np.linalg.norm(post - pre, axis=1)
    scale = np.linalg.norm(pre, axis=1)
    print(f"[shift] mean|pre|={scale.mean():.3f}  mean|post-pre|={shift.mean():.4f}  ratio={shift.mean()/scale.mean():.4%}")

    run(adata, "X_pca_pre_harmony", "no_harmony")
    run(adata, "X_pca_harmony", "with_harmony")


if __name__ == "__main__":
    main()
