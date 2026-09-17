#!/usr/bin/env python
"""Lightweight subsampled re-plot of the Figure 3 celltype-panel UMAPs.

The full-data panels in data/figure_3/pca_harmony_single_cell/<mod>/ plot every
observation (600k-707k points) -- correct, but heavy PNGs. This regenerates the
same two UMAP figure sets from a fixed random SUBSAMPLE of the *already-computed*
embedding (no re-clustering, no re-UMAP -- it reuses X_umap_pca_pre_harmony /
X_umap_pca_post_harmony / X_umap and cluster_label straight out of the leiden
output h5ad), writing to a sibling tree:

    data/figure_3/pca_harmony_single_cell_sub_sample_umap/<mod>/
        pca_harmony_qc/plots/umap_pca_harmony_before_after_by_{slice_id,cell_type_true,domain_true}.png
        leiden_pca_qc_celltype_matched/plots/umap_true_vs_predicted_{cell_type_true,domain_true}.png
        leiden_pca_qc_celltype_matched/plots/umap_by_{cluster_label,cell_type_true,domain_true,slice_id}.png
        leiden_pca_qc_celltype_matched/plots/contingency_{cell_type_true,domain_true}.png  (copied, full-data)
        ari_recovery_qc/ari_summary_<mod>.json                                             (copied)

Because it subsamples one shared coordinate set, the before/after panel and the
true-vs-predicted panel stay on the exact same embedding, same as the full run.
"""
import argparse
import importlib.util
import shutil
import sys
from pathlib import Path

import numpy as np
import scanpy as sc

_here = Path(__file__).resolve().parent
CLUSTERING_DIR = _here if (_here / "01.2_pca_harmony.py").exists() else _here.parent
sys.path.insert(0, str(CLUSTERING_DIR))
# 01.2_pca_harmony.py's numeric-leading name isn't a valid Python identifier,
# so it can't be `from ... import`-ed directly -- load it by file path instead.
_pca_harmony_spec = importlib.util.spec_from_file_location(
    "pca_harmony", CLUSTERING_DIR / "01.2_pca_harmony.py"
)
_pca_harmony = importlib.util.module_from_spec(_pca_harmony_spec)
_pca_harmony_spec.loader.exec_module(_pca_harmony)
plot_umap_before_after = _pca_harmony.plot_umap_before_after
from step03_cluster_and_plot import plot_umap, plot_umap_true_vs_predicted  # noqa: E402

FIG3 = Path("/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_3")
SRC_ROOT = FIG3 / "pca_harmony_single_cell"
DST_ROOT = FIG3 / "pca_harmony_single_cell_sub_sample_umap"


def replot(mod: str, n: int, seed: int, point_size: float, alpha: float) -> None:
    mm = "bin" if mod == "bin16um" else mod
    leiden_dir = SRC_ROOT / mod / "leiden_pca_qc_celltype_matched"
    hits = sorted(leiden_dir.glob(f"simulation_{mm}_z_leiden_pca_res*.h5ad"))
    if not hits:
        print(f"[{mod}] SKIP -- no leiden output under {leiden_dir} yet (full run not finished?)")
        return
    src_h5ad = hits[-1]
    print(f"[{mod}] reading {src_h5ad.name}")
    adata = sc.read_h5ad(src_h5ad)

    if adata.n_obs > n:
        rng = np.random.default_rng(seed)
        idx = np.sort(rng.choice(adata.n_obs, size=n, replace=False))
        sub = adata[idx].copy()
    else:
        sub = adata
    print(f"[{mod}] subsample {sub.n_obs:,} / {adata.n_obs:,}")

    ba_dir = DST_ROOT / mod / "pca_harmony_qc" / "plots"
    ba_dir.mkdir(parents=True, exist_ok=True)
    for ck in ("slice_id", "cell_type_true", "domain_true"):
        if ck in sub.obs:
            plot_umap_before_after(
                sub, ck, ba_dir / f"umap_pca_harmony_before_after_by_{ck}.png", point_size, alpha
            )

    tp_dir = DST_ROOT / mod / "leiden_pca_qc_celltype_matched" / "plots"
    tp_dir.mkdir(parents=True, exist_ok=True)
    for ck in ("cluster_label", "cell_type_true", "domain_true", "slice_id"):
        if ck in sub.obs:
            plot_umap(sub, ck, tp_dir / f"umap_by_{ck}.png")
    for tk in ("cell_type_true", "domain_true"):
        if tk in sub.obs and "cluster_label" in sub.obs:
            plot_umap_true_vs_predicted(
                sub, tk, "cluster_label", tp_dir / f"umap_true_vs_predicted_{tk}.png"
            )

    # contingency heatmaps + ARI json are cheap / full-data -- copy, don't subsample.
    for tk in ("cell_type_true", "domain_true"):
        src_png = leiden_dir / "plots" / f"contingency_{tk}.png"
        if src_png.exists():
            shutil.copy2(src_png, tp_dir / f"contingency_{tk}.png")
    ari_src = SRC_ROOT / mod / "ari_recovery_qc" / f"ari_summary_{mm}.json"
    if ari_src.exists():
        (DST_ROOT / mod / "ari_recovery_qc").mkdir(parents=True, exist_ok=True)
        shutil.copy2(ari_src, DST_ROOT / mod / "ari_recovery_qc" / ari_src.name)

    print(f"[{mod}] done -> {DST_ROOT / mod}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--modalities", nargs="+", default=["cell", "bin", "bin16um", "spot"])
    ap.add_argument("--n", type=int, default=50_000, help="Subsample size (plot all if fewer obs).")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--point-size", type=float, default=4.0)
    ap.add_argument("--alpha", type=float, default=0.85)
    args = ap.parse_args()
    for mod in args.modalities:
        replot(mod, args.n, args.seed, args.point_size, args.alpha)


if __name__ == "__main__":
    main()
