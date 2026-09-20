#!/usr/bin/env python
"""
Build a BANKSY-augmented embedding for a simulated albis dataset: BANKSY ->
PCA -> Harmony -> (this script stops here; step03_cluster_and_plot.py does
Leiden/Louvain from the output). This is the spatially-aware counterpart to
01.2_pca_harmony.py, which runs plain PCA -> Harmony with no spatial information
at all.

Spatial coordinates are staggered PER SLICE before the BANKSY neighbor graph
is built: every slice reuses the same local x/y coordinate range (they're all
the same capture window, just originally at different Z depths), so an
unstaggered neighbor graph would treat cells from different slices as spatial
neighbors of each other, which is meaningless.

IMPORTANT: this script requires the isolated `sim-app-banksy` conda env
(banksy_py pins a much older scanpy/numpy/anndata/pandas/scikit-learn/scipy
stack than the rest of this project -- see
sim_paper/env/create_banksy_env.sh). Everything downstream of this script's
output (step03_cluster_and_plot.py, etc.) runs in the normal
albis-tutorial env, since nothing after this step needs banksy_py itself.

Expected environment:
    conda activate /dcs04/hicks/data/Jan/sim_project/albis/env/sim-app-banksy

Example:
    python sim_paper/code/clustering/01_build_banksy_matrix.py \\
        --modality cell --packing-tag packing_pf0p04
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import umap as umap_lib
import harmonypy as hm

from banksy.initialize_banksy import initialize_banksy
from banksy.embed_banksy import generate_banksy_matrix
from banksy_utils.umap_pca import pca_umap


SCRIPT_DIR = Path(__file__).resolve().parent
SIM_PAPER_DIR = SCRIPT_DIR.parents[1]
VALID_MODALITIES = ("spot", "bin", "cell")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a BANKSY -> PCA -> Harmony embedding for a simulated dataset."
    )
    parser.add_argument("--modality", choices=VALID_MODALITIES, default="cell")
    parser.add_argument(
        "--packing-tag", default="packing_pf0p04",
        help="Reads data/noisy/<packing-tag>/simulation_<modality>_z.h5ad as input, unless --input is given.",
    )
    parser.add_argument("--input", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--batch-key", default="slice_id")
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument(
        "--n-hvg", type=int, default=None,
        help="Number of highly-variable genes to select. Default (None) uses every gene, matching "
        "01.2_pca_harmony.py's convention -- the simulated gene panel is only ~550 genes already, all "
        "either markers, shared markers, or deliberately-uninformative noise genes.",
    )
    parser.add_argument("--k-geom", type=int, default=15, help="Number of spatial neighbors per cell for BANKSY.")
    parser.add_argument("--max-m", type=int, default=1, help="Highest azimuthal-gradient-feature order BANKSY computes.")
    parser.add_argument(
        "--lambda", dest="lambda_value", type=float, default=0.8,
        help="BANKSY's neighborhood-weighting parameter (0 = pure cell-typing signal, ~0.8 = domain-segmentation signal).",
    )
    parser.add_argument("--nbr-weight-decay", default="scaled_gaussian", choices=("scaled_gaussian", "reciprocal", "ranked"))
    parser.add_argument(
        "--skip-umap", action="store_true",
        help="Skip the pre/post-Harmony UMAP computation and the UMAP diagnostic "
        "plots. Only X_pca_harmony is needed downstream (ARI / composition "
        "scoring); UMAP on 600k-700k obs dominates runtime in a parameter sweep.",
    )
    parser.add_argument(
        "--stagger-scale", type=float, default=1.5,
        help="Multiple of the widest slice's x-range used as the per-slice x-offset when staggering.",
    )
    parser.add_argument(
        "--plot-colors", nargs="+", default=["slice_id", "domain_true", "cell_type_true"],
        help="obs columns to color the diagnostic PCA/UMAP plots by.",
    )
    args = parser.parse_args()

    root = SIM_PAPER_DIR / "data" / "noisy" / args.packing_tag
    if args.input is None:
        args.input = root / f"simulation_{args.modality}_z.h5ad"
    if args.output_dir is None:
        # "_qc" in the dirname since this script is meant to run on QC-filtered input
        # (see 00_qc_filter.py / figure.md 2026-08-24) -- named explicitly so it's
        # never ambiguous with a hypothetical non-QC BANKSY run.
        args.output_dir = (
            SIM_PAPER_DIR / "data" / f"clustering_{args.packing_tag}" / args.modality / "banksy_pca_harmony_qc"
        )
    return args


# ── plotting helpers ──────────────────────────────────────────────────────

def make_color_lookup(labels, cmap=None):
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from manuscript_style import category_color, category_order
    key = labels.name or "slice_id"
    return {label: category_color(label, key) for label in category_order(labels)}


def scatter_by_label(ax, coords, labels_series, title, cmap=None, axis_labels=("dim 1", "dim 2")):
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from manuscript_style import apply_style, scatter_colors, legend_handles, pretty_label
    apply_style()
    key = labels_series.name or "slice_id"
    kwargs, labels, colors = scatter_colors(labels_series.to_frame(name=key), key)
    ax.scatter(coords[:, 0], coords[:, 1], **kwargs, s=0.3, linewidths=0,
               alpha=0.4, rasterized=True)
    ax.legend(handles=legend_handles(labels, colors), title=pretty_label(key),
              fontsize=12, title_fontsize=13, loc="upper center",
              bbox_to_anchor=(0.5, -0.18), frameon=False, ncol=min(len(labels), 4))
    ax.set_title(title, fontsize=17, fontweight="bold")
    ax.set_xlabel(axis_labels[0], fontsize=16)
    ax.set_ylabel(axis_labels[1], fontsize=16)
    ax.tick_params(labelsize=12)


def save_embedding_plot(embedding, obs, col, title, path, cmap="tab20", axis_labels=("dim 1", "dim 2")):
    if col not in obs.columns:
        print(f"  [plot] skipping {col} (not in obs)")
        return
    fig, ax = plt.subplots(figsize=(8, 6))
    scatter_by_label(ax, embedding, obs[col].astype(str), title, cmap=cmap, axis_labels=axis_labels)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  [plot] saved {path.name}")


def main() -> None:
    args = parse_args()
    if not args.input.is_file():
        raise SystemExit(f"Input not found: {args.input}\nRun the packing-fraction generator script first.")

    plot_dir = args.output_dir / "plots"
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    print(f"[load] {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"[load] AnnData shape: {adata.n_obs} x {adata.n_vars}")

    # Use post-resampling X at every batch sigma, including zero.

    # ── stagger spatial coordinates PER SLICE ────────────────────────────────
    # Every slice reuses the same local x/y coordinate range (they're all the
    # same capture window, just at different Z depths originally), so BANKSY's
    # neighbor graph -- which only looks at raw x/y -- would otherwise treat
    # cells from different slices as spatial neighbors of each other.
    print("[stagger] offsetting each slice's x coordinate so slices don't overlap")
    coord_keys = ("x_stagger", "y_stagger", "spatial")
    slice_ids = sorted(adata.obs[args.batch_key].unique())

    xy = np.asarray(adata.obsm["spatial"], dtype=float).copy()
    x_ranges = []
    for s in slice_ids:
        mask = (adata.obs[args.batch_key] == s).to_numpy()
        x_ranges.append(float(xy[mask, 0].max() - xy[mask, 0].min()))
    step = max(x_ranges) * args.stagger_scale

    for i, s in enumerate(slice_ids):
        mask = (adata.obs[args.batch_key] == s).to_numpy()
        x_min = xy[mask, 0].min()
        xy[mask, 0] = xy[mask, 0] - x_min + i * step
        print(f"  slice {s}: offset x by {i * step:.1f}")

    adata.obsm["spatial"] = xy
    adata.obs[coord_keys[0]] = xy[:, 0]
    adata.obs[coord_keys[1]] = xy[:, 1]

    # ── normalize + (optional) HVG selection ─────────────────────────────────
    print("\n[normalize] normalize_total (no log1p -- BANKSY expects normalized, non-log counts here)")
    adata.layers["counts"] = adata.X.copy()
    if args.n_hvg is not None:
        sc.pp.highly_variable_genes(
            adata, n_top_genes=args.n_hvg, flavor="seurat_v3", layer="counts", batch_key=args.batch_key,
        )
        hvg_genes = adata.var_names[adata.var["highly_variable"]].tolist()
        print(f"[hvg] selected {len(hvg_genes)} of {adata.n_vars} genes")
        adata = adata[:, hvg_genes].copy()
    else:
        print(f"[hvg] skipped -- using all {adata.n_vars} genes (matches 01.2_pca_harmony.py's convention)")
    sc.pp.normalize_total(adata, target_sum=1e4)

    # ── BANKSY ────────────────────────────────────────────────────────────────
    print(f"\n[banksy] k_geom={args.k_geom} max_m={args.max_m} lambda={args.lambda_value} decay={args.nbr_weight_decay}")
    banksy_dict = initialize_banksy(
        adata,
        coord_keys,
        args.k_geom,
        nbr_weight_decay=args.nbr_weight_decay,
        max_m=args.max_m,
        plt_edge_hist=False,
        plt_nbr_weights=False,
        plt_agf_angles=False,
        plt_theta=False,
    )
    banksy_dict, _ = generate_banksy_matrix(adata, banksy_dict, [args.lambda_value], args.max_m)

    bdata = banksy_dict[args.nbr_weight_decay][args.lambda_value]["adata"]
    bdata.obsm["spatial"] = adata.obsm["spatial"].copy()
    bdata.obs = adata.obs.copy()
    # carry over the soft ground-truth composition vectors (BANKSY doesn't
    # subset obs, so they stay row-aligned) -- composition_recovery.py needs
    # them to score bin/spot as mixtures rather than argmax labels.
    for k in list(adata.obsm):
        if k.endswith("_frac_true") and k not in bdata.obsm:
            bdata.obsm[k] = np.asarray(adata.obsm[k]).copy()
    print(f"[banksy] BANKSY matrix shape: {bdata.shape[0]} x {bdata.shape[1]}")
    del adata

    banksy_dict = {args.nbr_weight_decay: {args.lambda_value: {"adata": bdata}}}

    # ── PCA + UMAP (pre-Harmony) ─────────────────────────────────────────────
    print(f"\n[pca] pre-Harmony PCA{'' if args.skip_umap else ' + UMAP'}")
    pca_umap(banksy_dict, pca_dims=[args.n_pcs], add_umap=not args.skip_umap, plt_remaining_var=True)

    pc_key = f"reduced_pc_{args.n_pcs}"
    umap_key = f"{pc_key}_umap"
    bdata.obsm["X_pca_pre_harmony"] = bdata.obsm[pc_key].copy()
    if not args.skip_umap:
        bdata.obsm["X_umap_pre_harmony"] = bdata.obsm[umap_key].copy()

    for col in args.plot_colors:
        save_embedding_plot(
            bdata.obsm[pc_key], bdata.obs, col, f"BANKSY PCA pre-Harmony | {col}",
            plot_dir / f"pca_preharmony_by_{col}.png", axis_labels=("PC1", "PC2"),
        )
        if not args.skip_umap:
            save_embedding_plot(
                bdata.obsm[umap_key], bdata.obs, col, f"BANKSY UMAP pre-Harmony | {col}",
                plot_dir / f"umap_preharmony_by_{col}.png", axis_labels=("UMAP1", "UMAP2"),
            )

    # ── Harmony batch correction (harmonypy -- matches 01.2_pca_harmony.py) ───────
    print(f"\n[harmony] batch_key={args.batch_key}")
    pca = np.asarray(bdata.obsm[pc_key], dtype=np.float64)
    harmony_out = hm.run_harmony(pca, bdata.obs, args.batch_key)
    corrected = np.asarray(harmony_out.Z_corr)
    corrected = corrected.T if corrected.T.shape == pca.shape else corrected
    if corrected.shape != pca.shape:
        raise ValueError(f"Harmony returned shape {corrected.shape}; expected {pca.shape} or its transpose.")
    bdata.obsm["X_pca_harmony"] = corrected

    # ── UMAP (post-Harmony) ───────────────────────────────────────────────────
    if not args.skip_umap:
        print("\n[umap] post-Harmony UMAP")
        reducer = umap_lib.UMAP(random_state=args.random_state)
        bdata.obsm["X_umap_post_harmony"] = reducer.fit_transform(corrected)

    for col in args.plot_colors:
        save_embedding_plot(
            bdata.obsm["X_pca_harmony"], bdata.obs, col, f"BANKSY PCA post-Harmony | {col}",
            plot_dir / f"pca_postharmony_by_{col}.png", axis_labels=("PC1", "PC2"),
        )
        if not args.skip_umap:
            save_embedding_plot(
                bdata.obsm["X_umap_post_harmony"], bdata.obs, col, f"BANKSY UMAP post-Harmony | {col}",
                plot_dir / f"umap_postharmony_by_{col}.png", axis_labels=("UMAP1", "UMAP2"),
            )

    # ── write ──────────────────────────────────────────────────────────────
    # obsm["X_pca_harmony"] is the key step03_cluster_and_plot.py's
    # --pipeline pca_harmony expects -- point --input at this file and it
    # needs no changes to run Leiden/Louvain on the BANKSY-derived embedding.
    out_path = args.output_dir / f"simulation_{args.modality}_z_banksy_pca_harmony_qc.h5ad"
    print(f"\n[save] {out_path}")
    bdata.write_h5ad(out_path)
    print("[save] done")


if __name__ == "__main__":
    main()
