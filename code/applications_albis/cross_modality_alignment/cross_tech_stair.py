"""
Figure 4E -- cross-technology spatial alignment with STAIR.

Takes the SAME physical tissue slice (one slice_id from the ALBIS 10-slice
z-stack design) captured at three different technologies -- bin16um, spot,
cell -- and asks STAIR to align them into one shared spatial frame using
expression alone, exactly as `3D_stair.py` (Figure 4C) aligns different
z-sections of one technology. Here the three "slices" STAIR sees are the
three technologies of one physical slice, not ten z-planes of one technology.

Adapted from a multi-sample-alignment-benchmark cross-tech script (user-
supplied, itself adapted from that project's `3D_alignment/3D_stair.py`).
One change vs. that version, forced by a property of *this* project's
Figure 2 simulations -- checked empirically before writing this script:

  * 2026-09-15: cell used to be simulated at a much larger ~5985 um sphere
    radius (600k cells, ~0.16% 3D packing -- needed at that n_cells to avoid
    unrealistically crowding the sphere) while bin16um/spot shared a ~2050 um
    radius (packing_pf0p04, ~4% packing -- a compromise tuned to fix a
    DIFFERENT problem, bin/spot's empty-grid-cell artifact, that doesn't
    apply to cell). rescale_to_ref() used to isotropically rescale cell's
    coords onto the bin/spot disc by the resulting ~0.343 radius ratio
    before anything else. cell is now regenerated at the SAME sphere_r_um
    (2050) directly, with n_cells scaled down to 24207 (600000 *
    (2050/6000)**3) to preserve its original ~0.16% packing fraction rather
    than being forced up to bin/spot's 4% -- see DATA_VERSIONS.md and the
    2026-09-15 conversation this was built from. rescale_to_ref() is now a
    no-op for all three modalities (kept below for robustness/documentation,
    not because it does anything at these tags).
  * Reuses the *existing* obsm['spatial_unaligned'] per (modality,
    slice_id) as the "already misaligned" cross-tech input, rather than
    obsm['spatial'] directly. Each modality's slice_id==N unaligned coords
    carry a rigid (rotation + shift, NO scale -- verified: pairwise
    distances match obsm['spatial'] exactly) perturbation from
    uns['rigid_perturb_inplane']. Reusing this (rather than handing STAIR
    three already-co-registered point clouds) makes the task non-trivial
    and mirrors 3D_stair.py's prepare_input() exactly: obsm['spatial'] is
    stashed as ground truth, STAIR gets obsm['spatial_unaligned'] instead.
  * 2026-09-16: the three modalities' per-slice perturbations are now
    SYNCED -- Figure 2's generate_simulation_noisy.py calls pass
    --sync-unaligned-seed, so cell/bin16um/spot draw the identical
    rotation+translation for a given slice_id (previously each modality
    used its own independent draw, offset by a fixed per-modality constant;
    see sync_unaligned_seed in simulation_sphere.py). All three technologies
    for one slice now start from one shared ground-truth misalignment
    rather than three unrelated ones.

batch_key/slice_order becomes the three technology labels (bin16um / spot /
cell) instead of z-plane ids 0..9; n_neigh_hom / c_neigh_het take one shared
value across the whole run (8 / 0.90, CLI-overridable) rather than the
per-modality-tuned values 3D_stair.py uses for same-tech z-stacks, since here
all three modalities sit in one HGAT call together.

Output, per slice, under
    sim_paper/data/figure_4/cross_modality_alignment/STAIR/cross_tech/slice_<n>/
        adata_results/Sim_CrossTech_STAIR_slice_<n>.h5ad   obsm:
            spatial              -> STAIR input (per-technology unaligned coords)
            spatial_true         -> ground truth (cell rescaled onto bin/spot disc)
            transform_init       -> STAIR initial (MNN) alignment
            transform_fine       -> STAIR fine (ICP) alignment   [the result]
            STAIR                -> integrated cross-tech spatial embedding
        embeddings/attention_slice_<n>.csv
        align/...                STAIR alignment logs / edge plots
        metrics.json             joint-Procrustes RMSE before vs after, per technology
"""

import argparse
import json
import os
import random
import sys

import numpy as np
import anndata as ad
import scanpy as sc
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402,F401  (kept for parity with 3D_stair.py)

random.seed(42)
np.random.seed(42)

print(sys.executable)

import STAIR  # noqa: E402
from STAIR.emb_alignment import Emb_Align  # noqa: E402
from STAIR.loc_alignment import Loc_Align  # noqa: E402
from STAIR.utils import *  # noqa: E402,F401,F403  (provides cluster_func)

print("R in PATH:", os.popen("which R").read())
import rpy2.robjects as ro  # noqa: E402,F401
from rpy2.robjects import r, numpy2ri  # noqa: E402

numpy2ri.activate()
r('.libPaths("/users/jmatthia/R/4.3")')
r("library(mclust)")

import torch  # noqa: E402

torch.manual_seed(42)
torch.cuda.manual_seed_all(42)
print("Torch:", torch.__version__, "| CUDA:", torch.cuda.is_available(), torch.version.cuda)


# --------------------------------------------------------------------------- #
# dataset table -- same Figure 2 tags 3D_stair.py / Figure 4C use            #
# --------------------------------------------------------------------------- #
# 2026-09-17: figure_2 was split into parallel smaller_sphere (r=2050) /
# larger_sphere (r=6000) tracks; the tags below moved from bare
# `data/figure_2/<tag>/` to `data/figure_2/smaller_sphere/data/<tag>/`, which
# had silently broken this script (FileNotFoundError) since that move --
# same fix as 3D_stair.py's (Figure 4C), see that script's dataset-table
# comment for the full reorg context. Prior slice_5 output (Sep 16, predates
# the reorg) archived to STAIR/cross_tech/slice_5_pre_reorg_20260917/.
FIG2 = "/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_2/smaller_sphere/data"
BASE_OUTDIR = "/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_4/cross_modality_alignment/STAIR/cross_tech"

DATASETS = {
    "bin16um": f"{FIG2}/packing_pf0p04_bin16um_log_mu_-2.5_bsigma07/simulation_bin_z_qc.h5ad",
    "spot":    f"{FIG2}/packing_pf0p04_log_mu_-2.0_theta_0.25_jitter0.10_bsigma03/simulation_spot_z_qc.h5ad",
    "cell":    f"{FIG2}/log_mu_-2.3_theta_0.40_jitter0.15_bsigma15/simulation_cell_z.h5ad",
}
TECHS = ["bin16um", "spot", "cell"]   # fixed slice_order STAIR sees
REF_TECH = "bin16um"                  # rescale every other modality onto this disc
N_DOMAINS = 6                         # obs['domain_true'] has D0..D5


# --------------------------------------------------------------------------- #
# helpers                                                                     #
# --------------------------------------------------------------------------- #
def load_slice(tech, path, slice_id):
    a = sc.read_h5ad(path)
    a.obs["orig_slice_id"] = a.obs["slice_id"].astype(str)
    m = (a.obs["orig_slice_id"] == str(slice_id)).values
    if m.sum() == 0:
        raise ValueError(f"{tech}: no obs at slice_id={slice_id}")
    a = a[m].copy()
    a.obs["technology"] = tech
    return a


def rescale_to_ref(tech_adatas, ref_tech):
    """Isotropic rescale of every non-ref modality's spatial coords onto the
    ref modality's disc radius (max |xy| over obsm['spatial']). bin16um and
    spot already share a radius (~2050 um) -- only cell (~5985 um) needs
    this. Applies the same scale factor to 'spatial' and 'spatial_unaligned'
    since the latter is a rigid (rotation + shift, NO scale) perturbation of
    the former."""
    ref_radius = np.abs(np.asarray(tech_adatas[ref_tech].obsm["spatial"])[:, :2]).max()
    for tech, a in tech_adatas.items():
        if tech == ref_tech:
            a.uns["cross_tech_scale"] = 1.0
            continue
        radius = np.abs(np.asarray(a.obsm["spatial"])[:, :2]).max()
        scale = float(ref_radius / radius)
        for key in ("spatial", "spatial_unaligned"):
            xy = np.asarray(a.obsm[key], dtype=float).copy()
            xy[:, :2] *= scale
            a.obsm[key] = xy
        a.uns["cross_tech_scale"] = scale
        print(f"[rescale] {tech}: radius {radius:.1f} -> {ref_radius:.1f} um (x{scale:.4f})")
    return tech_adatas


def prepare_input(a):
    """Stash ground truth ('spatial' -> 'spatial_true'), hand STAIR the
    per-technology misaligned coords ('spatial_unaligned' -> 'spatial')."""
    sp_true = np.asarray(a.obsm["spatial"], dtype=float)[:, :2].copy()
    sp_unaligned = np.asarray(a.obsm["spatial_unaligned"], dtype=float)[:, :2].copy()
    a.obsm["spatial_true"] = sp_true
    a.obsm["spatial"] = sp_unaligned
    return a


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
    """Joint-Procrustes RMSE of unaligned vs STAIR-fine vs ground truth,
    overall and per technology (mirrors 3D_stair.py's per-slice scoring)."""
    true = adata.obsm["spatial_true"]
    out = {"n_obs": int(adata.n_obs)}
    for name, key in [("unaligned", "spatial"), ("stair_init", "transform_init"),
                      ("stair_fine", "transform_fine")]:
        if key not in adata.obsm:
            continue
        rmse, med = procrustes_rmse(adata.obsm[key], true)
        out[f"{name}_rmse_um"] = rmse
        out[f"{name}_median_um"] = med
    per_tech = {}
    for tech in sorted(adata.obs["slice_id"].unique()):
        m = (adata.obs["slice_id"] == tech).values
        rmse, _ = procrustes_rmse(adata.obsm["transform_fine"][m], true[m])
        per_tech[str(tech)] = rmse
    out["stair_fine_rmse_per_technology_um"] = per_tech
    return out


# --------------------------------------------------------------------------- #
# main routine                                                                #
# --------------------------------------------------------------------------- #
def run_cross_tech_alignment(tech_adatas, slice_id, n_neigh_hom, c_neigh_het,
                             output_adata, output_embeddings, output_align, used_device):
    tech_adatas = rescale_to_ref(tech_adatas, REF_TECH)
    tech_adatas = {t: prepare_input(a) for t, a in tech_adatas.items()}

    key_use = list(TECHS)
    slices = []
    for tech in key_use:
        sl = tech_adatas[tech]
        sl.var_names_make_unique()
        sl.obs_names = sl.obs_names.astype(str) + f"_{tech}"
        sc.pp.filter_cells(sl, min_genes=10)
        sc.pp.filter_genes(sl, min_cells=3)
        print(f"{tech} after filtering: {sl.shape}")
        slices.append(sl)

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
    atte.to_csv(os.path.join(output_embeddings, f"attention_slice_{slice_id}.csv"))

    # spatial clustering on the integrated embedding (matches 3D_stair.py)
    adata = cluster_func(adata, clustering="mclust", use_rep="STAIR",
                         cluster_num=N_DOMAINS, key_add="STAIR")

    # ---- location alignment (known technology order) ----------------------- #
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

    metrics = alignment_metrics(adata)
    metrics["slice_id"] = slice_id
    metrics["technologies"] = key_use
    metrics["n_neigh_hom"] = n_neigh_hom
    metrics["c_neigh_het"] = c_neigh_het
    with open(os.path.join(output_adata, "metrics.json"), "w") as fh:
        json.dump(metrics, fh, indent=2)
    print("METRICS:", json.dumps(metrics, indent=2))

    out_h5ad = os.path.join(output_adata, f"Sim_CrossTech_STAIR_slice_{slice_id}.h5ad")
    adata.write(out_h5ad)
    print("wrote", out_h5ad)


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--slice", type=int, default=5,
                    help="ALBIS slice_id (0-9) to take from each technology")
    ap.add_argument("--n-neigh-hom", type=int, default=8,
                    help="within-technology spatial neighbour count for the HGAT (shared across all 3 modalities)")
    ap.add_argument("--c-neigh-het", type=float, default=0.90,
                    help="cross-technology (heterogeneous) MNN edge cutoff for the HGAT")
    return ap.parse_args()


if __name__ == "__main__":
    args = parse_args()

    tech_adatas = {}
    for tech, path in DATASETS.items():
        if not os.path.exists(path):
            raise FileNotFoundError(path)
        tech_adatas[tech] = load_slice(tech, path, args.slice)
        print(f"{tech}: {tech_adatas[tech].shape} at slice_id={args.slice}")

    outdir = os.path.join(BASE_OUTDIR, f"slice_{args.slice}")
    output_adata = os.path.join(outdir, "adata_results")
    output_embeddings = os.path.join(outdir, "embeddings")
    output_align = os.path.join(outdir, "align")
    for d in (output_adata, output_embeddings, output_align):
        os.makedirs(d, exist_ok=True)

    print(f"Running cross-tech alignment at slice_id={args.slice}")
    print(f"Output: {outdir}")

    used_device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print("device:", used_device)

    run_cross_tech_alignment(
        tech_adatas=tech_adatas,
        slice_id=args.slice,
        n_neigh_hom=args.n_neigh_hom,
        c_neigh_het=args.c_neigh_het,
        output_adata=output_adata,
        output_embeddings=output_embeddings,
        output_align=output_align,
        used_device=used_device,
    )
