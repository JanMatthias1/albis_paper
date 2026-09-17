"""
Figure 4C -- 3D spatial alignment of an ALBIS sphere z-stack with STAIR.

Adapted from the multi-sample-alignment-benchmark script
    code/Jan/sim_data/3D_alignment/3D_stair.py

Differences from the benchmark version:
  * Inputs are ALBIS `simulation_<modality>_z.h5ad` files (Figure 2 tags). These
    already carry a full alignment problem + ground truth:
        obsm['spatial']              true, aligned in-plane coords (ground truth)
        obsm['spatial_unaligned']    per-slice rigid perturbation (rot + shift)
        obsm['spatial_3d'] / '..._unaligned'  the same, z stacked
        uns['rigid_perturb_inplane'] exact per-slice R, t, seed
        obs['slice_id'] 0..9, obs['domain_true'] (6), obs['cell_type_true'] (8)
  * Before running STAIR we stash the ground truth and feed STAIR the
    *misaligned* coords (spatial_unaligned -> spatial).
  * Known slice order only: batch_order stays sorted by slice_id, the STAIR
    z-reconstruction path (sort_slices / loc_predict_z) is left disabled.
  * Three modalities driven from a dataset table (bin16um / spot / cell) rather
    than a filename-prefix switch.

Output (per dataset) under
    sim_paper/data/figure_4/alignment/STAIR/<dataset>/
        adata_results/Sim_3D_STAIR_<dataset>.h5ad   full result, obsm:
            spatial              -> STAIR input (misaligned)
            spatial_true         -> ground-truth aligned XY
            spatial_3d_true      -> ground-truth aligned XYZ
            transform_init       -> STAIR initial (MNN) alignment
            transform_fine       -> STAIR fine (ICP) alignment   [the result]
            spatial_3d_stair     -> [transform_fine | true z]    reconstructed stack
            spatial_3d_unaligned -> [spatial_unaligned | true z]  disorganized stack
            STAIR                -> integrated spatial embedding
        embeddings/attention_<dataset>.csv
        align/...                STAIR alignment logs / edge plots
        metrics.json             joint-Procrustes RMSE before vs after
"""

import argparse
import json
import os
import random
import sys

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402  (kept for parity with benchmark script)

random.seed(42)
np.random.seed(42)

print(sys.executable)

import STAIR  # noqa: E402
from STAIR.emb_alignment import Emb_Align  # noqa: E402
from STAIR.loc_alignment import Loc_Align  # noqa: E402
from STAIR.loc_prediction import sort_slices, loc_predict_z  # noqa: E402,F401
from STAIR.utils import *  # noqa: E402,F401,F403  (provides cluster_func)

print("R in PATH:", os.popen("which R").read())
import rpy2.robjects as ro  # noqa: E402
from rpy2.robjects import r, numpy2ri  # noqa: E402

numpy2ri.activate()
r('.libPaths("/users/jmatthia/R/4.3")')
r("library(mclust)")

import torch  # noqa: E402

torch.manual_seed(42)
torch.cuda.manual_seed_all(42)
print("Torch:", torch.__version__, "| CUDA:", torch.cuda.is_available(), torch.version.cuda)


# --------------------------------------------------------------------------- #
# dataset table                                                              #
# --------------------------------------------------------------------------- #
FIG2 = "/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_2/smaller_sphere/data"
CDIST = "/dcs04/hicks/data/Jan/sim_project/sim_paper/data/count_distribution"
BASE_OUTDIR = "/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_4/alignment/STAIR"

# 2026-09-01: realwindow is now the baked-in default for the Figure 2 tags, so
# bin16um / spot carry a real 6.5 mm platform window with ~77% off-tissue
# bins/spots labelled domain_true="unassigned" / is_empty=True pre-QC. Point
# STAIR at the QC-filtered h5ad (00_qc_filter.py drops the empty border, leaving
# the clean 6-domain disc) -- the raw simulation_<mod>_z.h5ad would feed STAIR a
# mostly-empty stack and a 7th "unassigned" domain. cell has no off-tissue
# observations (12x24 mm window fully covered) so it stays on the raw h5ad.
# The former *_realwindow DATASETS entries + run_stair_3D_realwindow.sh were
# dropped here -- they pointed at data/count_distribution/..._realwindow/ paths
# that no longer exist, and are now redundant with these defaults.
#
# 2026-09-05: spot's Figure 2 config was retuned onto the two CytAssist probe
# references (log_mu -2.5 -> -2.0, theta 2.0 -> 0.25; see
# code/count_distribution/FIGURE2_METHODOLOGY.md). spot tag bumped here to
# match; pre-retune STAIR outputs archived at
# data/figure_4/alignment/STAIR/spot_pre_probe_retune_20260905/. bin16um / cell
# tags unchanged.
# 2026-09-06: spot tag bumped again -- ..._theta_0.25_bsigma03 ->
# ..._theta_0.25_jitter0.10_bsigma03. The shared spot config had no
# --theta-jitter so it fell back to the generator default 1.0, flooring ~40%
# of per-gene NB dispersions to 1e-3 (the disjoint upper cloud in Figure 2's
# mean_variance panel). --theta-jitter 0.10 collapses it. Same data fix that
# was applied to `cell` on 2026-08-25. Pre-jitter-fix STAIR outputs archived
# at data/figure_4/alignment/STAIR/spot_pre_jitter_fix_20260906/.
# 2026-09-17: figure_2 was split into parallel smaller_sphere (r=2050) /
# larger_sphere (r=6000) tracks (see project_figure2_smaller_larger_sphere
# memory); the same 4 tags below moved from bare `data/figure_2/<tag>/` to
# `data/figure_2/smaller_sphere/data/<tag>/` (FIG2 updated above), which had
# silently broken this script (FileNotFoundError) since that move -- never
# run end-to-end after it, same failure pattern as the clustering/ path bugs.
# Also added `bin8um` here for the first time: no prior STAIR run of this
# modality anywhere in the codebase, so n_neigh_hom/c_neigh_het below are an
# untuned first pass (reusing bin16um's values, the nearest analog) rather
# than a tuned choice. All 4 pre-reorg STAIR outputs (bin16um/spot/cell were
# stale vs. the reorg regardless of bin8um) archived to
# data/figure_4/alignment/STAIR/<dataset>_pre_reorg_20260917/.
DATASETS = {
    # 8 um Visium HD bin z-stack. Untuned first pass -- n_neigh_hom/c_neigh_het
    # reused from bin16um (nearest analog), not empirically tuned for this
    # modality. Also by far the largest of the 4 datasets (6.6M obs, ~4x
    # bin16um's 1.65M) -- resourced separately in run_stair_3D.sh; may need a
    # resubmit at higher mem/time if it OOMs or times out on the first try.
    "bin8um": dict(
        h5ad=f"{FIG2}/packing_pf0p04_log_mu_0.0_bsigma08/simulation_bin_z_qc.h5ad",
        n_neigh_hom=4,
        c_neigh_het=0.97,
    ),
    # 16 um Visium HD bin z-stack
    "bin16um": dict(
        h5ad=f"{FIG2}/packing_pf0p04_bin16um_log_mu_-2.5_bsigma07/simulation_bin_z_qc.h5ad",
        n_neigh_hom=4,
        c_neigh_het=0.97,
    ),
    # Visium spot z-stack
    "spot": dict(
        h5ad=f"{FIG2}/packing_pf0p04_log_mu_-2.0_theta_0.25_jitter0.10_bsigma03/simulation_spot_z_qc.h5ad",
        n_neigh_hom=8,
        c_neigh_het=0.90,
    ),
    # single-cell z-stack
    "cell": dict(
        h5ad=f"{FIG2}/log_mu_-2.3_theta_0.40_jitter0.15_bsigma15/simulation_cell_z.h5ad",
        n_neigh_hom=8,
        c_neigh_het=0.90,
    ),
}

N_DOMAINS = 6  # obs['domain_true'] has D0..D5


# --------------------------------------------------------------------------- #
# helpers                                                                    #
# --------------------------------------------------------------------------- #
def prepare_input(adata):
    """Stash ALBIS ground truth, hand STAIR the misaligned coords."""
    adata.obs["slice_id"] = adata.obs["slice_id"].astype(str)

    sp_true = np.asarray(adata.obsm["spatial"], dtype=float)[:, :2].copy()
    sp_unaligned = np.asarray(adata.obsm["spatial_unaligned"], dtype=float)[:, :2].copy()

    adata.obsm["spatial_true"] = sp_true
    if "spatial_3d" in adata.obsm:
        adata.obsm["spatial_3d_true"] = np.asarray(adata.obsm["spatial_3d"], dtype=float).copy()
    # STAIR reads obsm['spatial'] -> give it the disorganized stack
    adata.obsm["spatial"] = sp_unaligned
    return adata


def procrustes_rmse(X, Y, allow_scale=False, allow_reflection=True):
    """RMSE of the best similarity/rigid map X -> Y (per-row euclidean)."""
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)
    Xc = X - X.mean(0)
    Yc = Y - Y.mean(0)
    U, S, Vt = np.linalg.svd(Xc.T @ Yc)
    R = U @ Vt
    if not allow_reflection and np.linalg.det(R) < 0:
        U2 = U.copy()
        U2[:, -1] *= -1
        R = U2 @ Vt
    s = (S.sum() / (Xc ** 2).sum()) if allow_scale else 1.0
    Xhat = s * Xc @ R + Y.mean(0)
    d = np.linalg.norm(Xhat - Y, axis=1)
    return float(np.sqrt((d ** 2).mean())), float(np.median(d))


def alignment_metrics(adata):
    """Joint-Procrustes RMSE of unaligned vs STAIR-fine vs ground truth."""
    true = adata.obsm["spatial_true"]
    out = {"n_obs": int(adata.n_obs)}
    for name, key in [("unaligned", "spatial"), ("stair_init", "transform_init"),
                      ("stair_fine", "transform_fine")]:
        if key not in adata.obsm:
            continue
        rmse, med = procrustes_rmse(adata.obsm[key], true)
        out[f"{name}_rmse_um"] = rmse
        out[f"{name}_median_um"] = med
    # per-slice centroid error after the joint fit (fine)
    per_slice = {}
    for sid in sorted(adata.obs["slice_id"].unique(), key=float):
        m = (adata.obs["slice_id"] == sid).values
        rmse, _ = procrustes_rmse(adata.obsm["transform_fine"][m], true[m])
        per_slice[str(sid)] = rmse
    out["stair_fine_rmse_per_slice_um"] = per_slice
    return out


# --------------------------------------------------------------------------- #
# main routine                                                               #
# --------------------------------------------------------------------------- #
def run_3d_alignment(adata, dataset, n_neigh_hom, c_neigh_het,
                     output_adata, output_embeddings, output_align, used_device):
    adata = prepare_input(adata)

    key_use = sorted(adata.obs["slice_id"].unique(), key=lambda x: float(x))
    key_use = [str(x) for x in key_use]
    print("-----------------------", key_use, "-----------------------")

    slices = [adata[adata.obs["slice_id"] == sid].copy() for sid in key_use]
    for i, sl in enumerate(slices):
        sl.var_names_make_unique()
        sl.obs_names = sl.obs_names.astype(str) + f"_slice{i}"
        sc.pp.filter_cells(sl, min_genes=10)
        sc.pp.filter_genes(sl, min_cells=3)
        print(f"Slice {key_use[i]} after filtering: {sl.shape}")

    adata = ad.concat(slices, join="inner", label="batch", keys=key_use, uns_merge="unique")
    adata.obs["slice_id"] = adata.obs["batch"].astype(str)
    print(f"POST-CONCAT shape: {adata.shape}")

    # ---- embedding integration (HGAT) -------------------------------------- #
    emb_align = Emb_Align(adata, batch_key="slice_id", result_path=output_embeddings, device=used_device)
    emb_align.prepare()
    emb_align.preprocess()
    emb_align.latent()
    emb_align.prepare_hgat(
        spatial_key="spatial",
        slice_order=key_use,
        n_neigh_hom=n_neigh_hom,
        c_neigh_het=c_neigh_het,
    )
    emb_align.train_hgat(gamma=0.8, epoch_hgat=150)
    adata, atte = emb_align.predict_hgat()
    atte.to_csv(os.path.join(output_embeddings, f"attention_{dataset}.csv"))

    # spatial clustering on the integrated embedding (matches benchmark)
    adata = cluster_func(adata, clustering="mclust", use_rep="STAIR",
                         cluster_num=N_DOMAINS, key_add="STAIR")

    # ---- location alignment (known order, no z-reconstruction) ------------ #
    loc_align = Loc_Align(adata, batch_key="slice_id", batch_order=key_use, result_path=output_align)
    loc_align.init_align(emb_key="STAIR", spatial_key="spatial", num_mnn=5)
    loc_align.detect_fine_points(
        domain_key="STAIR",
        slice_boundary=True,
        domain_boundary=False,
        num_domains=1,
        alpha=500,
        return_result=False,
    )
    adata = loc_align.fine_align()

    # ---- reconstructed 3D stacks ----------------------------------------- #
    if "spatial_3d_true" in adata.obsm:
        z = np.asarray(adata.obsm["spatial_3d_true"])[:, 2]
    else:
        z = adata.obs["slice_id"].astype(float).values
    adata.obsm["spatial_3d_stair"] = np.column_stack([adata.obsm["transform_fine"], z])
    adata.obsm["spatial_3d_unaligned"] = np.column_stack([adata.obsm["spatial"], z])

    # ---- metrics -------------------------------------------------------- #
    metrics = alignment_metrics(adata)
    metrics["dataset"] = dataset
    with open(os.path.join(output_adata, "metrics.json"), "w") as fh:
        json.dump(metrics, fh, indent=2)
    print("METRICS:", json.dumps(metrics, indent=2))

    out_h5ad = os.path.join(output_adata, f"Sim_3D_STAIR_{dataset}.h5ad")
    adata.write(out_h5ad)
    print("wrote", out_h5ad)


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, choices=sorted(DATASETS),
                    help="which ALBIS z-stack to align")
    return ap.parse_args()


if __name__ == "__main__":
    args = parse_args()
    spec = DATASETS[args.dataset]

    data_path = spec["h5ad"]
    if not os.path.exists(data_path):
        raise FileNotFoundError(data_path)

    outdir = os.path.join(BASE_OUTDIR, args.dataset)
    output_adata = os.path.join(outdir, "adata_results")
    output_embeddings = os.path.join(outdir, "embeddings")
    output_align = os.path.join(outdir, "align")
    for d in (output_adata, output_embeddings, output_align):
        os.makedirs(d, exist_ok=True)

    print(f"Running dataset: {args.dataset}")
    print(f"Input:  {data_path}")
    print(f"Output: {outdir}")

    adata = sc.read_h5ad(data_path)
    used_device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print("device:", used_device)

    run_3d_alignment(
        adata=adata,
        dataset=args.dataset,
        n_neigh_hom=spec["n_neigh_hom"],
        c_neigh_het=spec["c_neigh_het"],
        output_adata=output_adata,
        output_embeddings=output_embeddings,
        output_align=output_align,
        used_device=used_device,
    )
