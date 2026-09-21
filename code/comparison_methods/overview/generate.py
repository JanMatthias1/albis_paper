"""One native 3D tissue per method, three outputs; use --help for invocation.

The parent uses only stdlib and dispatches to the existing method environments.
No empirical inputs, pretrained atlas, or custom spatial/aggregation algorithm.
"""
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import random
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
OLD = ROOT / "comparison_methods"
ENVS = OLD / "env"
TYPES = [f"type{i}" for i in range(1, 9)]
SPIDER_SOURCE = ROOT / "sim_paper/env/spider-overview-src"
SPIDER_COMMIT = "6ccd4da77257f2807c430f8f42fbe2dc175991de"


def run_command(command, log, env):
    with log.open("w") as handle:
        subprocess.run(list(map(str, command)), stdout=handle, stderr=subprocess.STDOUT,
                       env=env, check=True)


def dispatch(args, cfg):
    if args.spider_spots == "circle" and not (SPIDER_SOURCE / "COMMIT").exists():
        raise FileNotFoundError("Run setup_spider.sh before requesting upstream circular Spider spots")
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    if any(out.iterdir()):
        raise FileExistsError(f"Use a new output directory: {out}")
    config = out / "settings.json"
    config.write_text(json.dumps(cfg, indent=2) + "\n")
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", NUMBA_NUM_THREADS="1",
               MPLCONFIGDIR=str(out / "cache/matplotlib"),
               NUMBA_CACHE_DIR=str(out / "cache/numba"),
               XDG_CACHE_HOME=str(out / "cache"), PYTHONHASHSEED=str(cfg["seed"]))
    for p in ["cache/matplotlib", "cache/numba"]:
        (out / p).mkdir(parents=True)
    contract = out / "splatter_contract.json"
    contract.write_text(json.dumps(dict(n_cells=cfg["n_cells"], n_genes=cfg["n_genes"],
                                       n_cell_types=len(TYPES), cell_type_proportions=cfg["proportions"])))
    print(f"Generating synthetic Splatter input; logs and outputs: {out}", flush=True)
    run_command([ENVS / "splatter/bin/Rscript", OLD / "code/workflows/splatter_expression.R",
                 "--contract", contract, "--seed", cfg["seed"], "--out-dir", out / "splatter"],
                out / "splatter.log", env)
    manifest = json.loads((out / "splatter/manifest.json").read_text())
    if manifest["status"] != "ok":
        raise RuntimeError(manifest["failure_reason"])
    for method, environment in [("albis", "our_method"), ("sccube", "sccube"), ("spider", "spider")]:
        print(f"Generating {method}: one tissue, then bin/spot/cell outputs", flush=True)
        run_command([ENVS / environment / "bin/python", __file__, "--method", method,
                     "--settings", config, "--out", out, "--spider-spots", args.spider_spots],
                    out / f"{method}.log", env)
    run_command([ENVS / "analysis/bin/python", HERE / "plot.py", "--input", out],
                out / "plot.log", env)
    print(f"Finished: {out / 'overview.pdf'}", flush=True)


def cell_data(X, genes, labels, xyz, cfg):
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy.sparse import csr_matrix
    a = ad.AnnData(csr_matrix(X), obs=pd.DataFrame({"cell_type_true": list(labels)},
                   index=[f"cell{i}" for i in range(len(labels))]),
                   var=pd.DataFrame(index=list(genes)))
    a.obsm["spatial_3d"] = np.asarray(xyz, dtype=float)
    a.obsm["spatial"] = a.obsm["spatial_3d"][:, :2].copy()
    # Known volume boundaries, not per-slice min/max; include the upper edge.
    a.obs["slice_id"] = np.clip(np.floor(a.obsm["spatial_3d"][:, 2] /
                               cfg["extent_um"] * cfg["n_slices"]), 0, cfg["n_slices"] - 1).astype(int)
    return a


def capture_data(X, genes, xy, counts, slice_id, cfg):
    import numpy as np
    z = (slice_id + .5) * cfg["extent_um"] / cfg["n_slices"]
    n = counts.sum(axis=1)
    labels = np.array(TYPES)[counts.argmax(axis=1)]
    labels = np.where(n > 0, labels, "unassigned")
    a = cell_data(X, genes, labels, np.column_stack([xy, np.full(len(xy), z)]), cfg)
    a.obs["slice_id"] = slice_id
    a.obs["n_source_cells"] = n
    a.obs["is_empty"] = n == 0
    a.obsm["cell_type_counts"] = counts
    return a


def albis_outputs(cfg):
    import numpy as np
    sys.path.insert(0, str(ROOT / "albis"))
    from albis.simulation_sphere import simulate_3d_molecule_sphere_multires
    # Preserve ALBIS's native gene model: 80 markers/type => 556 genes for eight types.
    mix = np.full((len(TYPES), len(TYPES)), .30 / (len(TYPES) - 1))
    np.fill_diagonal(mix, .70)
    result = simulate_3d_molecule_sphere_multires(
        n_cells=cfg["n_cells"], n_cell_types=len(TYPES), n_domains=len(TYPES), domain_type_mix=mix,
        marker_genes_per_type=80, sphere_R_um=cfg["extent_um"] / 2,
        n_slices=cfg["n_slices"], batch_sigma=0, max_deg=0, max_shift=0,
        capture_window_um=False, xenium_capture_window_um=False,
        bin_size_um=cfg["bin_width_um"], spot_spacing_um=cfg["spot_spacing_um"],
        spot_radius_um=cfg["spot_radius_um"], seed=cfg["seed"],
        output_modalities=("cell", "bin", "spot"), slice_axes=("Z",))
    outputs = {"cell": result["adata_cell_sectioned"]["Z"],
               "bin": result["bin_adatas"]["Z"], "spot": result["spot_adatas"]["Z"]}
    for a in outputs.values():
        a.obsm["spatial_3d_native"] = a.obsm["spatial_3d"].copy()
        a.obsm["spatial_3d"] = a.obsm["spatial_3d"] + cfg["extent_um"] / 2
        a.obsm["spatial"] = a.obsm["spatial_3d"][:, :2].copy()
    return outputs


def competitor_outputs(method, cfg, out, spider_spots):
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy.io import mmread
    pool = out / "splatter"
    X = mmread(pool / "counts.mtx").T.tocsr()
    genes = (pool / "genes.tsv").read_text().splitlines()
    meta = pd.read_csv(pool / "cells.tsv", sep="\t")
    if set(meta.cell_type_true) != set(TYPES):
        raise ValueError("Splatter pool must contain all eight types; increase pilot cell count")
    if method == "spider":
        if spider_spots == "circle":
            if (SPIDER_SOURCE / "COMMIT").read_text().strip() != SPIDER_COMMIT:
                raise ValueError("Unexpected Spider source revision; run setup_spider.sh")
            sys.path.insert(0, str(SPIDER_SOURCE))
        import spider
        # Top-level alias still points to legacy sim_naive; use the corrected
        # native sim_expr implementation explicitly, without modifying its code.
        from spider.sim_expr import get_sim_spot_level_expr
        prior = np.array(cfg["proportions"])
        trans = np.tile(prior, (len(TYPES), 1))
        np.fill_diagonal(trans, 0)
        trans *= (1 - cfg["spider_self_probability"]) / trans.sum(axis=1, keepdims=True)
        np.fill_diagonal(trans, cfg["spider_self_probability"])
        labels, xyz = spider.simulate_10X_3d(
            cell_num=cfg["n_cells"], Num_celltype=len(TYPES), prior=prior, target_trans=trans,
            image_width=cfg["extent_um"], image_height=cfg["extent_um"], image_depth=cfg["extent_um"],
            smallsample_max_iter=cfg["spider_max_iterations"])
        reference = ad.AnnData(X, obs=meta.set_index("cell_id"), var=pd.DataFrame(index=genes))
        expression = spider.get_sim_cell_level_expr(
            celltype_assignment=labels, adata=reference, Num_celltype=len(TYPES),
            Num_ct_sample=np.bincount(labels, minlength=len(TYPES)), match_list=TYPES, ct_key="cell_type_true")
        cells = cell_data(expression.X, genes, np.array(TYPES)[labels], xyz, cfg)
    else:
        import torch
        from scCube.sccube import scCube
        torch.manual_seed(cfg["seed"])
        torch.set_num_threads(1)
        model = scCube()
        metadata = pd.DataFrame({"Cell": meta.cell_id, "Cell_type": meta.cell_type_true})
        reference = model.pre_process(pd.DataFrame(X.T.toarray(), index=genes, columns=meta.cell_id),
                                      metadata, is_normalized=False)
        target_num = meta.cell_type_true.value_counts().reindex(TYPES).astype(int).to_dict()
        generated_meta, data = model.train_vae_and_generate_cell(
            reference, celltype_key="Cell_type", cell_key="Cell", target_num=target_num,
            epoch_num=cfg["sccube_epochs"], used_device="cpu", save_model=True,
            save_path=str(out / "sccube_model"), project_name="synthetic_only")
        data, generated_meta = model.generate_pattern_random(
            data, generated_meta, spatial_dim=3, spatial_size=cfg["sccube_grid_size"],
            delta=cfg["sccube_delta"], lamda=cfg["sccube_lamda"], is_split=False,
            set_seed=True, seed=cfg["seed"])
        generated_meta = generated_meta.loc[data.columns]
        xyz = generated_meta[["point_x", "point_y", "point_z"]].to_numpy()
        scale = cfg["extent_um"] / cfg["sccube_grid_size"]
        cells = cell_data(data.T.to_numpy(), data.index, generated_meta.Cell_type, xyz * scale, cfg)
        cells.obsm["spatial_3d_native"] = xyz
    outputs = {"cell": cells}
    for modality in ("bin", "spot"):
        sections = []
        for sid in range(cfg["n_slices"]):
            sub = cells[cells.obs.slice_id == sid].copy()
            if sub.n_obs == 0:
                raise ValueError(f"Empty native cell section {sid}; increase n_cells or reduce n_slices")
            if method == "spider":
                codes = pd.Categorical(sub.obs.cell_type_true, categories=TYPES).codes
                aggregate = get_sim_spot_level_expr if spider_spots == "circle" else spider.get_sim_spot_level_expr
                result = aggregate(
                    Num_sample=sub.n_obs,
                    spot_diameter=(cfg["bin_width_um"] if modality == "bin" else
                                   2 * cfg["spot_radius_um"] if spider_spots == "circle" else cfg["spot_spacing_um"]),
                    image_width=cfg["extent_um"], image_height=cfg["extent_um"],
                    celltype_assignment=codes, cell_spatial=sub.obsm["spatial"], sim_cell_expr=sub,
                    gap=cfg["spot_spacing_um"], coord_type="generic",
                    spot_generate_type="square" if modality == "bin" else spider_spots, cell_coord_type="generic")
                expr, xy, membership = result[:3]
                counts = np.asarray(membership @ np.eye(len(TYPES), dtype=int)[codes])
                if modality == "spot" and spider_spots == "circle":
                    distance2 = ((xy[:, None, :] - sub.obsm["spatial"][None, :, :]) ** 2).sum(axis=2)
                    expected = distance2 <= cfg["spot_radius_um"] ** 2
                    if not np.array_equal(membership.toarray(), expected):
                        raise ValueError("Native circular capture does not match the requested radius")
                    if not np.allclose(expr.toarray(), expected @ sub.X.toarray()):
                        raise ValueError("Native circular capture expression differs from contributing cells")
                # Check native aggregation against its own membership matrix.
                if modality == "bin" and not np.allclose(np.asarray(expr.sum(axis=0)), np.asarray(sub.X.sum(axis=0))):
                    raise ValueError("Spider square aggregation did not conserve expression")
                if modality == "bin" and not np.all(np.asarray(membership.sum(axis=0)).ravel() == 1):
                    raise ValueError("Spider did not assign every cell exactly once")
            else:
                ids = sub.obs_names
                # Native routine relies on columns 2/3 being x/y. Drop other metadata.
                native_meta = pd.DataFrame({"Cell": ids, "Cell_type": sub.obs.cell_type_true.to_numpy(),
                                           "point_x": sub.obsm["spatial"][:, 0] / scale,
                                           "point_y": sub.obsm["spatial"][:, 1] / scale})
                source = pd.DataFrame(sub.X.toarray().T, index=sub.var_names, columns=ids)
                expr_df, loc, membership = model.generate_spot_data_random(
                    source, native_meta, platform="ST" if modality == "bin" else "Visium",
                    gene_type="whole", n_cell=cfg[f"sccube_cells_per_{modality}"])
                if len(membership) != sub.n_obs or membership.Cell.nunique() != sub.n_obs:
                    raise ValueError("scCube lost or duplicated cell memberships")
                loc = loc.set_index("spot").loc[expr_df.columns]
                counts = pd.crosstab(membership.spot, membership.Cell_type).reindex(
                    index=expr_df.columns, columns=TYPES, fill_value=0).to_numpy()
                if counts.sum() != sub.n_obs:
                    raise ValueError("scCube membership does not cover this section")
                expr, xy = expr_df.T.to_numpy(), loc[["spot_x", "spot_y"]].to_numpy(dtype=float) * scale
            sections.append(capture_data(expr, cells.var_names, xy, counts, sid, cfg))
        outputs[modality] = ad.concat(sections, index_unique="-slice", merge="same")
    return outputs


# albis_outputs() never reads n_genes/proportions/sccube_*/spider_*; those only
# feed the shared Splatter pool (dispatch()) or the competitor methods. Reporting
# the full shared cfg under albis's manifest wrongly implies it consumed them.
ALBIS_SETTINGS_KEYS = ["seed", "n_cells", "n_slices", "extent_um",
                       "bin_width_um", "spot_spacing_um", "spot_radius_um"]


def worker(args, cfg):
    import numpy as np
    random.seed(cfg["seed"])
    np.random.seed(cfg["seed"])
    destination = args.out / args.method
    destination.mkdir(parents=True, exist_ok=False)
    outputs = albis_outputs(cfg) if args.method == "albis" else competitor_outputs(args.method, cfg, args.out, args.spider_spots)
    reported_settings = {k: cfg[k] for k in ALBIS_SETTINGS_KEYS} if args.method == "albis" else cfg
    report = dict(method=args.method, settings=reported_settings, empirical_input=False, pretrained_atlas=False,
                  expression_input="internal" if args.method == "albis" else "Splatter synthetic pool",
                  coordinate_units="um; scCube units scaled uniformly to requested extent",
                  slice_ids="zero-based; capture z is section midpoint; cell z is continuous",
                  spider_spots=args.spider_spots, outputs={})
    if args.method == "spider":
        report["spider_upstream_commit"] = SPIDER_COMMIT if args.spider_spots == "circle" else "installed"
        import spider
        report["spider_source_path"] = spider.__file__
        report["spider_runtime_version"] = getattr(spider, "__version__", "unknown")
        report["aggregation_function"] = "spider.sim_expr.get_sim_spot_level_expr" if args.spider_spots == "circle" else "spider.get_sim_spot_level_expr"
    report["versions"] = {}
    for package in ["numpy", "anndata", "scCube", "spider", "albis", "torch"]:
        try:
            report["versions"][package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            pass
    report["script_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    for modality, a in outputs.items():
        xyz = np.asarray(a.obsm["spatial_3d"])
        if xyz.shape != (a.n_obs, 3) or not np.isfinite(xyz).all():
            raise ValueError(f"Invalid coordinates: {modality}")
        if set(a.obs.slice_id.astype(int)) != set(range(cfg["n_slices"])):
            raise ValueError(f"Missing sections: {modality}")
        values = a.X.data if hasattr(a.X, "tocsr") else np.asarray(a.X)
        if not np.isfinite(values).all() or (values < 0).any():
            raise ValueError(f"Invalid expression: {modality}")
        if not set(a.obs.cell_type_true.astype(str)).issubset(set(TYPES) | {"unassigned"}):
            raise ValueError(f"Unexpected cell labels: {modality}")
        if not a.obs_names.is_unique:
            raise ValueError(f"Duplicate observation identifiers: {modality}")
        a.uns["overview_method"] = args.method
        if args.method == "spider":
            a.uns["spider_spot_geometry"] = args.spider_spots
        a.uns["overview_settings_json"] = json.dumps(cfg)
        a.write_h5ad(destination / f"{modality}.h5ad")
        report["outputs"][modality] = dict(n_observations=a.n_obs, n_genes=a.n_vars,
                                          slices=a.obs.slice_id.value_counts().sort_index().to_dict())
    (destination / "manifest.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--settings", type=Path, default=HERE / "settings.json")
    parser.add_argument("--out", type=Path, required=True, help="New directory; never overwrite existing runs")
    parser.add_argument("--method", choices=["all", "albis", "sccube", "spider"], default="all")
    parser.add_argument("--spider-spots", choices=["circle", "square"], default="circle",
                        help="circle: pinned upstream; square: existing installed Spider")
    args = parser.parse_args()
    cfg = json.loads(args.settings.read_text())
    if len(cfg["proportions"]) != len(TYPES) or abs(sum(cfg["proportions"]) - 1) > 1e-8:
        raise ValueError("This overview expects eight types with proportions summing to 1")
    dispatch(args, cfg) if args.method == "all" else worker(args, cfg)
