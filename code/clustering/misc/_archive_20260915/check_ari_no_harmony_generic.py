#!/usr/bin/env python
"""
BANKSY -> PCA -> Leiden, NO Harmony: resolution-matched Leiden/ARI on
X_pca_pre_harmony (BANKSY's own PCA, before Harmony) for an arbitrary
banksy_matrix h5ad already built by 01_build_banksy_matrix.py. No rebuild
needed -- X_pca_pre_harmony is always stored regardless of whether Harmony
was run downstream. Writes a summary JSON in the same shape as
ari_vs_ground_truth.py's output (ground_truth/target_k/resolution/
achieved_k/ari) plus slice_leak_ari, so existing tabulation code can read
either file with minimal changes.

Usage:
    python check_ari_no_harmony_generic.py --input <banksy_matrix.h5ad> --output-dir <dir> --modality spot
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import scanpy as sc
from sklearn.metrics import adjusted_rand_score

CLUSTERING_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CLUSTERING_DIR))
from ari_vs_ground_truth import find_resolution_for_k  # noqa: E402


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--modality", required=True)
    p.add_argument("--n-pcs", type=int, default=30)
    p.add_argument("--n-neighbors", type=int, default=15)
    p.add_argument("--random-state", type=int, default=0)
    p.add_argument("--res-lo", type=float, default=0.01)
    p.add_argument("--res-hi", type=float, default=3.0)
    p.add_argument("--max-iter", type=int, default=15)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    print(f"[load] {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"[load] {adata.n_obs:,} x {adata.n_vars:,}")

    n_pcs = min(args.n_pcs, adata.obsm["X_pca_pre_harmony"].shape[1])
    rep_key = "X_pca_pre_harmony_ari"
    adata.obsm[rep_key] = np.asarray(adata.obsm["X_pca_pre_harmony"])[:, :n_pcs].copy()
    print(f"[neighbors] BANKSY -> PCA -> Leiden, NO Harmony. n_pcs={n_pcs} n_neighbors={args.n_neighbors}")
    sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, use_rep=rep_key, random_state=args.random_state)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    results = []
    for truth_col in ("domain_true", "cell_type_true"):
        if truth_col not in adata.obs:
            continue
        target_k = int(adata.obs[truth_col].nunique())
        cluster_key = f"leiden_no_harmony_{truth_col}"
        best_res, best_k, tried = find_resolution_for_k(
            adata, target_k, args.res_lo, args.res_hi, args.max_iter, args.random_state, cluster_key,
        )
        ari = adjusted_rand_score(adata.obs[truth_col], adata.obs[cluster_key])
        leak = adjusted_rand_score(adata.obs["slice_id"], adata.obs[cluster_key]) if "slice_id" in adata.obs else float("nan")
        print(f"[RESULT] {args.modality} / {truth_col} (NO Harmony): k={best_k} (target {target_k})  ARI={ari:.4f}  slice_leak_ARI={leak:.4f}")
        results.append({
            "modality": args.modality, "ground_truth": truth_col, "target_k": target_k,
            "resolution": float(best_res), "achieved_k": int(best_k), "ari": float(ari),
            "slice_leak_ari": float(leak), "n_trials": len(tried), "trials": tried,
        })

    out_path = args.output_dir / f"ari_summary_no_harmony_{args.modality}.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"[save] {out_path}")


if __name__ == "__main__":
    main()
