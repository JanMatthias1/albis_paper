#!/usr/bin/env python
"""
ARI between re-clustered data and ground-truth domain/cell-type labels, for
Figure 3B: how well does spatial-domain and cell-type structure survive
bin/spot/cell aggregation and re-clustering?

Builds a neighbor graph from the pca_harmony-corrected PCA embedding
(matching 03_clustering_plots.py's pca_harmony defaults: first --n-pcs
dimensions, --n-neighbors), then separately for domain_true and
cell_type_true, binary-searches the Leiden resolution until the number of
clusters found matches the true number of categories -- so ARI isn't
confounded by over/under-clustering -- and reports
sklearn.metrics.adjusted_rand_score against each.

Expected environment:
    conda activate /dcs04/hicks/data/Jan/sim_project/albis/env/albis-tutorial

Example:
    python sim_paper/code/clustering/02_leiden_resolution_sweep.py --modality bin --packing-tag packing_pf0p04
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import scanpy as sc
from sklearn.metrics import adjusted_rand_score

SCRIPT_DIR = Path(__file__).resolve().parent
SIM_PAPER_DIR = SCRIPT_DIR.parents[1]

GROUND_TRUTH_COLS = ("domain_true", "cell_type_true")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="ARI between Leiden clusters (resolution matched to true category count) and ground truth."
    )
    parser.add_argument("--modality", choices=("bin", "spot", "cell"), required=True)
    parser.add_argument(
        "--packing-tag", default="packing_pf0p04",
        help="Reads data/clustering_<packing-tag>/<modality>/pca_harmony/ as input, unless --input is given.",
    )
    parser.add_argument("--input", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument("--res-lo", type=float, default=0.01)
    parser.add_argument("--res-hi", type=float, default=3.0)
    parser.add_argument("--max-iter", type=int, default=15)
    args = parser.parse_args()

    root = SIM_PAPER_DIR / "data" / f"clustering_{args.packing_tag}"
    if args.input is None:
        args.input = root / args.modality / "pca_harmony" / f"simulation_{args.modality}_z_pca_harmony.h5ad"
    if args.output_dir is None:
        args.output_dir = root / args.modality / "ari_recovery"
    return args


def find_resolution_for_k(adata, target_k, res_lo, res_hi, max_iter, random_state, cluster_key):
    """Bisection search on Leiden resolution to hit target_k clusters. Leiden's cluster
    count is not strictly monotonic in resolution, so this keeps the closest-to-target
    trial seen rather than assuming the final bisection step is best."""
    lo, hi = res_lo, res_hi
    best_res, best_k, best_diff = None, None, None
    tried = []
    for _ in range(max_iter):
        mid = (lo + hi) / 2
        sc.tl.leiden(
            adata, resolution=mid, key_added=cluster_key, random_state=random_state,
            flavor="igraph", n_iterations=2, directed=False,
        )
        k = adata.obs[cluster_key].nunique()
        tried.append({"resolution": mid, "n_clusters": int(k)})
        diff = abs(k - target_k)
        if best_diff is None or diff < best_diff:
            best_res, best_k, best_diff = mid, k, diff
        if k == target_k:
            break
        elif k < target_k:
            lo = mid
        else:
            hi = mid

    # Re-run at best_res so adata.obs[cluster_key] reflects the winning trial, not
    # whichever resolution the loop happened to try last.
    sc.tl.leiden(
        adata, resolution=best_res, key_added=cluster_key, random_state=random_state,
        flavor="igraph", n_iterations=2, directed=False,
    )
    return best_res, best_k, tried


def main() -> None:
    args = parse_args()
    if not args.input.is_file():
        raise SystemExit(f"Input not found: {args.input}\nRun pca_harmony.py first.")

    print(f"[load] {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"[load] AnnData shape: {adata.n_obs} x {adata.n_vars}")

    n_pcs = min(args.n_pcs, adata.obsm["X_pca_harmony"].shape[1])
    rep_key = "X_pca_harmony_ari"
    adata.obsm[rep_key] = np.asarray(adata.obsm["X_pca_harmony"])[:, :n_pcs].copy()
    print(f"[neighbors] n_pcs={n_pcs} n_neighbors={args.n_neighbors}")
    sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, use_rep=rep_key, random_state=args.random_state)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    results = []
    for truth_col in GROUND_TRUTH_COLS:
        target_k = int(adata.obs[truth_col].nunique())
        cluster_key = f"leiden_{truth_col}"
        print(f"[ari] {truth_col}: target_k={target_k}")
        best_res, best_k, tried = find_resolution_for_k(
            adata, target_k, args.res_lo, args.res_hi, args.max_iter, args.random_state, cluster_key,
        )
        ari = adjusted_rand_score(adata.obs[truth_col], adata.obs[cluster_key])
        print(f"[ari] {truth_col}: resolution={best_res:.4f} achieved_k={best_k} (target {target_k})  ARI={ari:.4f}")
        results.append({
            "modality": args.modality,
            "ground_truth": truth_col,
            "target_k": target_k,
            "resolution": float(best_res),
            "achieved_k": int(best_k),
            "ari": float(ari),
            "n_trials": len(tried),
            "trials": tried,
        })

    summary_path = args.output_dir / f"ari_summary_{args.modality}.json"
    with open(summary_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"[save] {summary_path}")

    out_h5ad = args.output_dir / f"simulation_{args.modality}_z_ari_recovery.h5ad"
    adata.write_h5ad(out_h5ad)
    print(f"[save] {out_h5ad}")


if __name__ == "__main__":
    main()
