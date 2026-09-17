#!/usr/bin/env python
"""
Before/after-Harmony UMAPs for the BANKSY batch-compare job (35608415),
matching 01.2_pca_harmony.py's plain-pipeline demo panels
(umap_pca_harmony_before_after_by_<color>.png) pixel-for-pixel -- same
subsampling, neighbors/UMAP params, panel geometry, and color/legend
handling, reused directly from 01.2_pca_harmony.py rather than reimplemented.

Reads data/figure_3/banksy_batch_compare/<modality>/<batch>/banksy_matrix/
simulation_<modality>_z_banksy_pca_harmony_qc.h5ad (already has
obsm["X_pca_pre_harmony"] / obsm["X_pca_harmony"] from 01_build_banksy_matrix.py,
--skip-umap so no UMAP was computed there -- this script adds it as a
lightweight post-hoc step on the already-built PCA, subsampled the same way
01.2_pca_harmony.py subsamples for its own diagnostic UMAPs).

Usage:
    conda activate /dcs04/hicks/data/Jan/sim_project/sim_paper/env/albis-tutorial
    python plot_banksy_batch_compare_umap.py --modality cell --batch tuned
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
plot_umap_before_after = _pca_harmony.plot_umap_before_after
PANEL_FIGSIZE = _pca_harmony.PANEL_FIGSIZE
from pc_pairs import sampled_indices  # noqa: E402

SIM_PAPER_DIR = CLUSTERING_DIR.parents[1]
ROOT = SIM_PAPER_DIR / "data" / "figure_3" / "banksy_batch_compare"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--modality", required=True, choices=["cell", "bin", "spot"])
    p.add_argument("--batch", default="prebatch", choices=["prebatch", "tuned"], help="Ignored if --input is given (only used for the default path/label).")
    p.add_argument("--input", type=Path, default=None, help="Override the default banksy_batch_compare path (e.g. for a batch_sigma-slide run).")
    p.add_argument("--label", default=None, help="Label used in printouts when --input is given (e.g. 'bs0.5').")
    p.add_argument("--n-neighbors", type=int, default=15)
    p.add_argument("--random-state", type=int, default=0)
    p.add_argument("--umap-max-obs", type=int, default=50000)
    p.add_argument("--point-size", type=float, default=4.0)
    p.add_argument("--alpha", type=float, default=0.85)
    p.add_argument("--plot-colors", nargs="+", default=["domain_true", "slice_id", "cell_type_true"])
    return p.parse_args()


def main() -> None:
    args = parse_args()
    if args.input is not None:
        h5ad_path = args.input
        plot_dir = args.input.parent / "plots"
    else:
        run_dir = ROOT / args.modality / args.batch
        h5ad_path = run_dir / "banksy_matrix" / f"simulation_{args.modality}_z_banksy_pca_harmony_qc.h5ad"
        plot_dir = run_dir / "banksy_matrix" / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    print(f"[load] {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    print(f"[load] {adata.n_obs:,} x {adata.n_vars:,}")

    obs_idx = sampled_indices(adata.n_obs, args.umap_max_obs, False, args.random_state)
    if len(obs_idx) < adata.n_obs:
        print(f"[sample] {len(obs_idx):,} / {adata.n_obs:,} observations (matches 01.2_pca_harmony.py's default)")
    umap_adata = adata[obs_idx].copy() if len(obs_idx) < adata.n_obs else adata

    # 01_build_banksy_matrix.py writes X_pca_harmony (post-Harmony); 01.2_pca_harmony.py's
    # shared plotting code expects X_pca_post_harmony -- alias rather than rename on disk.
    umap_adata.obsm["X_pca_post_harmony"] = umap_adata.obsm["X_pca_harmony"]

    print("[umap] before Harmony")
    compute_umap(umap_adata, "X_pca_pre_harmony", "X_umap_pca_pre_harmony", args.n_neighbors, args.random_state)
    print("[umap] after Harmony")
    compute_umap(umap_adata, "X_pca_post_harmony", "X_umap_pca_post_harmony", args.n_neighbors, args.random_state)

    for color_key in args.plot_colors:
        if color_key not in umap_adata.obs:
            print(f"[skip] {color_key} not in obs")
            continue
        out_path = plot_dir / f"umap_pca_harmony_before_after_by_{color_key}.png"
        plot_umap_before_after(umap_adata, color_key, out_path, args.point_size, args.alpha)
        print(f"[save] {out_path}")


if __name__ == "__main__":
    main()
