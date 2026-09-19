"""
Figure 4B -- 3D spatial-domain identification on an ALBIS sphere z-stack with
STAGATE (PyG port).

Adapted from STAGATE Tutorial 5 ("3D spatial domain identification"):
    https://stagate.readthedocs.io/en/latest/T5_3D.html
The tutorial builds a 3D spatial network (2D SNN within each section + edges
between adjacent sections) so the autoencoder can share signal across
consecutive sections and smooth out section-specific technical noise.

Differences from the tutorial:
  * Input is an ALBIS `simulation_<modality>_z_qc.h5ad` file (a locked Figure 2
    tag, post-batch QC). It already carries:
        obsm['spatial']        true aligned in-plane XY (ground truth)
        obsm['spatial_3d']     the same, with the real z stacked on
        obs['slice_id']        0..9  -> the "sections"
        obs['domain_true']     D0..D5 (6)  -> the domain label STAGATE targets
        obs['cell_type_true']  8 types      -> secondary reference
    Unlike Figure 4C (alignment), this task assumes the stack is ALREADY
    aligned -- we feed STAGATE the true coords and ask it to recover domains.
  * Sections are obs['slice_id'] (not a Puck id); section_order is the sorted
    slice_id list.
  * rad_cutoff_2D / rad_cutoff_Zaxis default to a multiple of the median
    nearest-neighbour spacing measured per section, so one script works across
    bin16um / spot / cell without hand-tuned radii (override with --rad-2d /
    --rad-z / --rad-mult).
  * Runs STAGATE twice -- with the 3D network and with the 2D-only network --
    and reports ARI/NMI for both, which is the comparison Figure 4B makes.

Output (per dataset) under
    sim_paper/data/figure_4/spatial_clustering/STAGATE/<dataset>/
        adata_results/Sim_3D_STAGATE_<dataset>.h5ad
            obsm['STAGATE'] / obsm['STAGATE_2D']   the two embeddings
            obs['mclust_3d'] / obs['mclust_2d']    mclust labels (k = 6)
        metrics.json    ARI / NMI vs domain_true (primary) and cell_type_true,
                        3D vs 2D, plus per-slice ARI for the 3D run
        plots/          3D domain scatter (true | STAGATE-3D | STAGATE-2D),
                        UMAP of the 3D embedding coloured by domain / slice
"""

import argparse
import json
import os
import random
import sys

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

from sklearn.neighbors import NearestNeighbors  # noqa: E402
from sklearn.metrics import (  # noqa: E402
    adjusted_rand_score,
    normalized_mutual_info_score,
)

print(sys.executable)

import torch  # noqa: E402
import STAGATE_pyG as ST  # noqa: E402

from plot_stagate_figures import plot_3d_panels, plot_umap  # noqa: E402

SEED = 0
random.seed(SEED)
np.random.seed(SEED)
torch.manual_seed(SEED)
torch.cuda.manual_seed_all(SEED)
print("Torch:", torch.__version__, "| CUDA:", torch.cuda.is_available(), torch.version.cuda)


# --------------------------------------------------------------------------- #
# dataset table                                                              #
# --------------------------------------------------------------------------- #
from stagate_inputs import DATASETS

BASE_OUTDIR = "/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_4/spatial_clustering/STAGATE"

N_DOMAINS = 6  # obs['domain_true'] has D0..D5
N_TOP_GENES = 3000  # HVG cap; the ALBIS panel is 556 genes so this keeps them all


# --------------------------------------------------------------------------- #
# helpers                                                                    #
# --------------------------------------------------------------------------- #
def prepare_input(adata, n_top_genes=N_TOP_GENES):
    """Tutorial-style normalisation; section key + 2D coords for STAGATE."""
    adata.var_names_make_unique()
    adata.obs["slice_id"] = adata.obs["slice_id"].astype(str)

    # true aligned in-plane coords -> obsm['spatial'] (this is the domain-ID
    # task: the stack is assumed already aligned)
    adata.obsm["spatial"] = np.asarray(adata.obsm["spatial"], dtype=float)[:, :2].copy()

    sc.pp.highly_variable_genes(adata, flavor="seurat_v3",
                                n_top_genes=min(n_top_genes, adata.n_vars))
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


def auto_radius(adata, key_section="slice_id", mult=2.5):
    """Median nearest-neighbour spacing across sections, times `mult`."""
    dists = []
    for sec in adata.obs[key_section].unique():
        xy = adata.obsm["spatial"][(adata.obs[key_section] == sec).values]
        if xy.shape[0] < 2:
            continue
        nn = NearestNeighbors(n_neighbors=2).fit(xy)
        d, _ = nn.kneighbors(xy)
        dists.append(np.median(d[:, 1]))
    nn_med = float(np.median(dists))
    return nn_med, nn_med * mult


def _knn_edges(coor_a, idx_a, coor_b, idx_b, k):
    """k nearest points of `b` for every point of `a`; (Cell1, Cell2, Distance) rows."""
    k = int(min(k, len(idx_b)))
    if k < 1:
        return []
    nn = NearestNeighbors(n_neighbors=k).fit(coor_b)
    dist, ind = nn.kneighbors(coor_a)
    rows = []
    for i in range(ind.shape[0]):
        for j in range(ind.shape[1]):
            rows.append((idx_a[i], idx_b[ind[i, j]], float(dist[i, j])))
    return rows


def cal_spatial_net_3d_knn(adata, k_2d, k_z, key_section="slice_id",
                           section_order=None, verbose=True):
    """KNN analogue of ST.Cal_Spatial_Net_3D.

    A fixed radius starves the irregular cell-resolution point cloud (see the
    'cell' note in DATASETS); KNN gives every cell exactly k neighbours no
    matter the local density. Populates the same uns keys the library builds:
    Spatial_Net_2D (within-section, k_2d NN), Spatial_Net_Zaxis (between
    adjacent sections, k_z NN queried BOTH directions so cross-section edges
    are guaranteed), and their concat Spatial_Net.
    """
    obs_names = np.asarray(adata.obs_names)
    xy = np.asarray(adata.obsm["spatial"], dtype=float)
    sec = adata.obs[key_section].astype(str).values
    if section_order is None:
        section_order = sorted(np.unique(sec), key=float)

    net_2d = []
    for s in section_order:
        m = np.where(sec == s)[0]
        df = pd.DataFrame(_knn_edges(xy[m], obs_names[m], xy[m], obs_names[m], k_2d + 1),
                          columns=["Cell1", "Cell2", "Distance"])
        df = df[df["Distance"] > 0]
        df["SNN"] = s
        net_2d.append(df)
        if verbose:
            print(f"------2D KNN section {s}: {len(df)} edges, {len(m)} cells "
                  f"({len(df) / max(len(m), 1):.2f}/cell)")
    net_2d = pd.concat(net_2d, ignore_index=True)

    net_z = []
    for s1, s2 in zip(section_order[:-1], section_order[1:]):
        m1, m2 = np.where(sec == s1)[0], np.where(sec == s2)[0]
        rows = (_knn_edges(xy[m1], obs_names[m1], xy[m2], obs_names[m2], k_z)
                + _knn_edges(xy[m2], obs_names[m2], xy[m1], obs_names[m1], k_z))
        df = pd.DataFrame(rows, columns=["Cell1", "Cell2", "Distance"])
        df["SNN"] = f"{s1}-{s2}"
        net_z.append(df)
        if verbose:
            print(f"------Z KNN {s1}-{s2}: {len(df)} cross-section edges")
    net_z = pd.concat(net_z, ignore_index=True)

    adata.uns["Spatial_Net_2D"] = net_2d
    adata.uns["Spatial_Net_Zaxis"] = net_z
    adata.uns["Spatial_Net"] = pd.concat([net_2d, net_z], ignore_index=True)


def run_stagate(adata, spatial_net, n_epochs, key_added):
    """One STAGATE fit against a given (2D or 3D) spatial network."""
    a = adata.copy()
    a.uns["Spatial_Net"] = spatial_net.copy()
    a = ST.train_STAGATE(a, n_epochs=n_epochs, key_added=key_added,
                         device=torch.device("cuda" if torch.cuda.is_available() else "cpu"))
    return a.obsm[key_added]


def cluster_scores(labels_pred, adata):
    out = {}
    for ref in ("domain_true", "cell_type_true"):
        if ref not in adata.obs:
            continue
        y = adata.obs[ref].astype(str).values
        out[f"ari_{ref}"] = float(adjusted_rand_score(y, labels_pred))
        out[f"nmi_{ref}"] = float(normalized_mutual_info_score(y, labels_pred))
    return out


def per_slice_ari(labels_pred, adata, ref="domain_true"):
    out = {}
    lab = pd.Series(labels_pred, index=adata.obs_names)
    for sid in sorted(adata.obs["slice_id"].unique(), key=float):
        m = (adata.obs["slice_id"] == sid).values
        out[str(sid)] = float(adjusted_rand_score(
            adata.obs.loc[m, ref].astype(str).values, lab[m].values))
    return out


# --------------------------------------------------------------------------- #
# plotting -- plot_3d_panels / plot_umap now live in plot_stagate_figures.py #
# (CPU-only) so plots can be regenerated without a GPU allocation           #
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, choices=sorted(DATASETS))
    ap.add_argument("--n-epochs", type=int, default=500)
    ap.add_argument("--rad-mult", type=float, default=1.5,
                    help="rad_cutoff = rad_mult * median NN spacing (per-section); "
                         "~1.5 -> first neighbour ring, ~2.5 -> ~3 rings")
    ap.add_argument("--rad-2d", type=float, default=None, help="override rad_cutoff_2D (um)")
    ap.add_argument("--rad-z", type=float, default=None, help="override rad_cutoff_Zaxis (um)")
    ap.add_argument("--graph-model", choices=("radius", "knn"), default=None,
                    help="spatial graph: 'radius' (default) or 'knn' (fixed degree; "
                         "overrides DATASETS[...]['graph_model'])")
    ap.add_argument("--k-2d", type=int, default=None, help="KNN: within-section neighbours (default 6)")
    ap.add_argument("--k-z", type=int, default=None, help="KNN: between-section neighbours each way (default 3)")
    ap.add_argument("--subsample", type=int, default=0,
                    help="randomly keep this many obs before building the graph (0 = all)")
    args = ap.parse_args()

    spec = DATASETS[args.dataset]
    if not os.path.exists(spec["h5ad"]):
        raise FileNotFoundError(spec["h5ad"])

    run_tag = args.dataset
    outdir = os.path.join(BASE_OUTDIR, run_tag)
    out_adata = os.path.join(outdir, "adata_results")
    out_plots = os.path.join(outdir, "plots")
    for d in (out_adata, out_plots):
        os.makedirs(d, exist_ok=True)

    print(f"[load] {spec['h5ad']}")
    adata = sc.read_h5ad(spec["h5ad"])
    print(f"[load] {adata.shape}")

    # Always use post-resampling X; sigma=0 retains Poisson sampling.

    if args.subsample and args.subsample < adata.n_obs:
        rng = np.random.default_rng(SEED)
        keep = rng.choice(adata.n_obs, size=args.subsample, replace=False)
        adata = adata[np.sort(keep)].copy()
        print(f"[subsample] -> {adata.shape}")

    adata = prepare_input(adata)

    # section_order MUST hold the same string values as obs['slice_id']
    # (Cal_Spatial_Net_3D matches them with .isin) -- just sorted numerically.
    section_order = sorted(adata.obs["slice_id"].unique(), key=float)
    print("[sections]", section_order)

    nn_med, rad_auto = auto_radius(adata, mult=args.rad_mult)

    # precedence: CLI flag > per-dataset spec > default. 'knn' builds a fixed-degree
    # graph (needed for the sparse/irregular cell point cloud); 'radius' is the
    # tutorial default and stays the choice for bin16um / spot.
    graph_model = (args.graph_model if args.graph_model is not None
                   else spec.get("graph_model", "radius"))
    rad_2d = rad_z = float("nan")
    k_2d = k_z = None

    if graph_model == "knn":
        k_2d = args.k_2d if args.k_2d is not None else spec.get("k_2d", 6)
        k_z = args.k_z if args.k_z is not None else spec.get("k_z", 3)
        print(f"[graph] KNN 3D net: k_2d = {k_2d}, k_z = {k_z} "
              f"(median NN spacing = {nn_med:.2f} um)")
        cal_spatial_net_3d_knn(adata, k_2d=k_2d, k_z=k_z, key_section="slice_id",
                               section_order=section_order, verbose=True)
    else:
        rad_2d = (args.rad_2d if args.rad_2d is not None
                  else spec.get("rad_2d", rad_auto))
        rad_z = (args.rad_z if args.rad_z is not None
                 else spec.get("rad_z", rad_auto))
        print(f"[radius] median NN spacing = {nn_med:.2f} um -> "
              f"rad_2d = {rad_2d:.2f}, rad_z = {rad_z:.2f}")
        ST.Cal_Spatial_Net_3D(adata, rad_cutoff_2D=rad_2d, rad_cutoff_Zaxis=rad_z,
                              key_section="slice_id", section_order=section_order,
                              verbose=True)
    net_3d = adata.uns["Spatial_Net"].copy()
    net_2d = adata.uns["Spatial_Net_2D"].copy()
    deg_3d = len(net_3d) / adata.n_obs
    deg_2d = len(net_2d) / adata.n_obs
    print(f"[net] 3D edges = {len(net_3d)} ({deg_3d:.1f}/cell)  | "
          f"2D-only edges = {len(net_2d)} ({deg_2d:.1f}/cell)")
    # STAGATE needs a real spatial graph; below ~3 nbrs/cell the GAT has nothing
    # to propagate and mclust on the embedding collapses to ~chance ARI.
    if deg_2d < 3.0:
        knob = (f"k_2d={k_2d} (--k-2d)" if graph_model == "knn"
                else f"rad_2d={rad_2d:.1f} um (--rad-2d / --rad-mult)")
        raise RuntimeError(
            f"spatial graph too sparse: {deg_2d:.2f} 2D nbrs/cell [{knob}]. "
            "Raise it (or set it in DATASETS) so degree is ~6-10.")

    # ---- STAGATE with the 3D network ----------------------------------------
    print("[stagate] 3D network")
    adata.obsm["STAGATE"] = run_stagate(adata, net_3d, args.n_epochs, "STAGATE")
    adata = ST.mclust_R(adata, N_DOMAINS, used_obsm="STAGATE")
    adata.obs["mclust_3d"] = adata.obs["mclust"].astype(str)

    # ---- STAGATE with the 2D-only network (tutorial comparison) ------------
    print("[stagate] 2D-only network")
    adata.obsm["STAGATE_2D"] = run_stagate(adata, net_2d, args.n_epochs, "STAGATE_2D")
    adata = ST.mclust_R(adata, N_DOMAINS, used_obsm="STAGATE_2D")
    adata.obs["mclust_2d"] = adata.obs["mclust"].astype(str)

    # ---- metrics ---------------------------------------------------------
    metrics = {
        "dataset": args.dataset,
        "input_h5ad": os.path.realpath(spec["h5ad"]),
        "input_mtime_ns": os.stat(spec["h5ad"]).st_mtime_ns,
        "training_seed": SEED,
        "simulation_seed": spec["simulation_seed"],
        "technology": spec["technology"],
        "batch_sigma": spec["batch_sigma"],
        "use_pre_batch": False,
        "input_layer": "X",
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "n_epochs": args.n_epochs,
        "median_nn_spacing_um": nn_med,
        "graph_model": graph_model,
        "rad_cutoff_2D": float(rad_2d),
        "rad_cutoff_Zaxis": float(rad_z),
        "k_2d": k_2d,
        "k_z": k_z,
        "n_edges_3d": int(len(net_3d)),
        "n_edges_2d": int(len(net_2d)),
        "deg_3d": deg_3d,
        "deg_2d": deg_2d,
        "stagate_3d": cluster_scores(adata.obs["mclust_3d"].values, adata),
        "stagate_2d": cluster_scores(adata.obs["mclust_2d"].values, adata),
        "stagate_3d_ari_domain_true_per_slice": per_slice_ari(
            adata.obs["mclust_3d"].values, adata),
    }
    with open(os.path.join(out_adata, "metrics.json"), "w") as fh:
        json.dump(metrics, fh, indent=2)
    print("METRICS:", json.dumps(metrics, indent=2))

    # ---- plots + save --------------------------------------------------
    z = (np.asarray(adata.obsm["spatial_3d"])[:, 2]
         if "spatial_3d" in adata.obsm else adata.obs["slice_id"].astype(float).values)
    plot_3d_panels(adata, z, out_plots)
    plot_umap(adata, out_plots)

    out_h5ad = os.path.join(out_adata, f"Sim_3D_STAGATE_{args.dataset}.h5ad")
    adata.write(out_h5ad)
    print("wrote", out_h5ad)


if __name__ == "__main__":
    main()
