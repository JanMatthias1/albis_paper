#!/usr/bin/env python
"""
Draft Figure 3B (spatial-domain) qualitative panel for a BANKSY run whose
matrix was built with --skip-umap (e.g. the banksy_lambda_kgeom_sweep_16um
outputs). Loads one *_ari_recovery.h5ad, computes UMAP on the pre- and
post-Harmony BANKSY-PCA embeddings, and writes the same plot set the
*_domain_panel.sh scripts produce via step03_cluster_and_plot.py:

  pre/post-Harmony UMAP, before/after side-by-side, by slice_id / domain_true
  UMAP by domain_true / predicted-cluster / slice_id
  UMAP true (domain_true) vs predicted, same coords
  contingency heatmap  predicted x domain_true

Reuses the plotting functions from 01.2_pca_harmony.py + step03_cluster_and_plot.py
so the draft matches the real panels exactly.

Usage:
  env/albis-tutorial/bin/python code/clustering/misc/plot_banksy_domain_panel.py \
      --h5ad data/figure_3/banksy_lambda_kgeom_sweep_16um/cell/lam0.5_kg200/ari/simulation_cell_z_ari_recovery.h5ad \
      --modality cell --out-dir <dir>
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np
import scanpy as sc
from sklearn.metrics import adjusted_rand_score

CLUSTERING_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CLUSTERING_DIR))
# 01.2_pca_harmony.py's numeric-leading name isn't a valid Python identifier,
# so it can't be `from ... import`-ed directly -- load it by file path instead.
_pca_harmony_spec = importlib.util.spec_from_file_location(
    "pca_harmony", CLUSTERING_DIR / "01.2_pca_harmony.py"
)
_pca_harmony = importlib.util.module_from_spec(_pca_harmony_spec)
_pca_harmony_spec.loader.exec_module(_pca_harmony)
plot_umap_before_after = _pca_harmony.plot_umap_before_after
from step03_cluster_and_plot import (  # noqa: E402
    plot_umap, plot_umap_true_vs_predicted, plot_contingency_heatmap,
)


def _umap(adata, rep_key, out_obsm, seed):
    sc.pp.neighbors(adata, use_rep=rep_key, n_neighbors=15, random_state=seed)
    sc.tl.umap(adata, random_state=seed)
    adata.obsm[out_obsm] = adata.obsm["X_umap"].copy()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", type=Path, required=True)
    ap.add_argument("--modality", required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--pred-key", default="leiden_domain_true",
                    help="obs column with the domain-matched Leiden partition")
    ap.add_argument("--pre-key", default="X_pca_pre_harmony")
    ap.add_argument("--post-key", default="X_pca_harmony")
    ap.add_argument("--n", type=int, default=60000, help="subsample for plotting/UMAP speed")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--point-size", type=float, default=4.0)
    ap.add_argument("--alpha", type=float, default=0.85)
    args = ap.parse_args()

    adata = sc.read_h5ad(args.h5ad)
    print(f"[load] {args.h5ad.name}  {adata.n_obs} obs")
    for k in (args.pre_key, args.post_key):
        if k not in adata.obsm:
            sys.exit(f"obsm['{k}'] missing -- have {list(adata.obsm)}")
    if args.pred_key not in adata.obs:
        sys.exit(f"obs['{args.pred_key}'] missing -- have {list(adata.obs)}")

    adata.obs[args.pred_key] = adata.obs[args.pred_key].astype(str).astype("category")
    for k in ("domain_true", "cell_type_true", "slice_id"):
        if k in adata.obs:
            adata.obs[k] = adata.obs[k].astype(str).astype("category")

    ari = adjusted_rand_score(adata.obs["domain_true"], adata.obs[args.pred_key])
    print(f"[ari] domain_true vs {args.pred_key}: {ari:.4f}  "
          f"(k_pred={adata.obs[args.pred_key].nunique()}, k_true={adata.obs['domain_true'].nunique()})")

    if adata.n_obs > args.n:
        idx = np.sort(np.random.default_rng(args.seed).choice(adata.n_obs, args.n, replace=False))
        adata = adata[idx].copy()
        print(f"[subsample] {adata.n_obs} obs")

    _umap(adata, args.pre_key, "X_umap_pca_pre_harmony", args.seed)
    _umap(adata, args.post_key, "X_umap_pca_post_harmony", args.seed)
    adata.obsm["X_umap"] = adata.obsm["X_umap_pca_post_harmony"].copy()   # post-Harmony for the panels

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    ba, lp = out / "pre_post_harmony", out / "leiden_domain_matched"
    ba.mkdir(exist_ok=True); lp.mkdir(exist_ok=True)

    for key in ("slice_id", "domain_true", "cell_type_true"):
        if key in adata.obs:
            plot_umap_before_after(adata, key, ba / f"umap_before_after_by_{key}.png",
                                   args.point_size, args.alpha)
    for key in ("domain_true", args.pred_key, "slice_id"):
        if key in adata.obs:
            plot_umap(adata, key, lp / f"umap_by_{key}.png")
    plot_umap_true_vs_predicted(adata, "domain_true", args.pred_key,
                                lp / "umap_true_vs_predicted_domain_true.png")
    plot_contingency_heatmap(adata, "domain_true", args.pred_key,
                             lp / "contingency_domain_true.png")

    (out / "ari.txt").write_text(f"domain_true vs {args.pred_key}: ARI = {ari:.4f}\n")
    print(f"[done] {out}")


if __name__ == "__main__":
    main()
