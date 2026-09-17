#!/usr/bin/env python
"""
Figure 4A -- spatially variable gene (SVG) identification with nnSVG
(Bioconductor, Weber lab) on a single representative ALBIS slice.

Ground truth: see 3D_scbsp.py's docstring -- same `is_noise` binary-recovery
setup (56 `noise` genes = true negative, remaining 500 = true positive,
`marker_shared` counted as signal per 2026-09-04 decision).

Unlike 3D_scbsp.py / SPARK-X, nnSVG does NOT support 3D coordinates: its
BRISC::BRISC_order dependency hard-requires exactly 2 columns (see project
memory, 2026-09-05 -- no 3D mode exists in nnSVG's vignette/NEWS/help either,
despite that being a reasonable expectation). So this script restricts to a
single slice (`--slice-id`, default 5, matching Figure 2's own slice-5
convention) rather than pooling the z-stack.

nnSVG also expects log-normalized counts (`sc.pp.normalize_total` + `log1p`),
not raw counts like scBSP/SPARK-X take -- a real methodological difference,
not just a CLI default.

By default `--input` resolves to the canonical weak-mix Figure 2 h5ad. Pass
`--strongmix` to point instead at the strong-domain-mix sibling dataset
(`data/figure_4/spatial_clustering/sim_data/<tag>_strongmix/`, the same data
Figure 4B's STAGATE work uses) -- outputs are tagged `<dataset>_strongmix`
so a weak-mix and a strongmix run never collide. `--use-pre-batch` appends
`_prebatch` on top of that (oracle run, no batch noise). Use `--run-tag` to
override the output stem entirely if you need a non-standard combination.

Per-gene cost is much higher than scBSP/SPARK-X (~2s total): nnSVG fits a
per-gene nearest-neighbor Gaussian process via BRISC, roughly linear in both
n_genes and n_obs. Empirically (2026-09-05, single slice): spot (~1.3k obs)
~150s, bin16um (~50k obs) ~50min, cell (~88k obs) ~90min projected, all with
n_threads=4 (BiocParallel across genes; each gene's own BRISC fit is
single-threaded).

Expected environment:
    conda activate /dcs04/hicks/data/Jan/sim_project/sim_paper/env/stagate-pyg
    (nnSVG installed 2026-09-05 -- also needs its Rscript for run_nnsvg.R)

Example:
    python sim_paper/code/applications_albis/svg/3D_nnsvg.py --dataset spot
    python sim_paper/code/applications_albis/svg/3D_nnsvg.py --dataset spot --strongmix
    python sim_paper/code/applications_albis/svg/3D_nnsvg.py --dataset spot --strongmix --use-pre-batch
"""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io as sio
from sklearn.metrics import average_precision_score, roc_auc_score

SCRIPT_DIR = Path(__file__).resolve().parent
SIM_PAPER_DIR = SCRIPT_DIR.parents[2]
R_BIN_PREFIX = Path("/dcs04/hicks/data/Jan/sim_project/sim_paper/env/stagate-pyg")
RSCRIPT = R_BIN_PREFIX / "lib" / "R" / "bin" / "Rscript"
RUN_NNSVG_R = SCRIPT_DIR / "run_nnsvg.R"

# Same DATASETS mapping as 3D_scbsp.py.
#
# NOTE (2026-09-06): the "spot" entry here is DELIBERATELY NOT the spot dataset
# used by Figure 2 / Figure 3 / Figure 4C. Those moved to
# `packing_pf0p04_log_mu_-2.0_theta_0.25_jitter0.10_bsigma03` (2026-09-05
# CytAssist probe-reference retune + 2026-09-06 `--theta-jitter 0.10`
# two-cloud fix). Figure 4A (and 4B/STAGATE) stay on the older
# `packing_pf0p04_log_mu_-2.5_bsigma03` spot config on purpose: the 4A/4B
# conclusions are about spatial-signal strength (domain-mix coupling), not
# real-data count-distribution fidelity, so the spot count tuning doesn't
# change them, and re-running would also need a matching `..._jitter0.10_..._strongmix`
# sibling that doesn't exist. bin16um and cell already match the current
# shared tags. If 4A ever needs to be consistent with the rest of Figure 4,
# bump the "spot" tag here + generate the strongmix sibling.
DATASETS = {
    "bin16um": ("packing_pf0p04_bin16um_log_mu_-2.5_bsigma07", "bin"),
    "spot": ("packing_pf0p04_log_mu_-2.5_bsigma03", "spot"),
    "cell": ("log_mu_-2.3_theta_0.40_jitter0.15_bsigma15", "cell"),
}

# Canonical weak-mix h5ads live under data/figure_2/<tag>/; the strong-domain-mix
# siblings (--strongmix) live here, generated for the Figure 4B STAGATE work.
STRONGMIX_DIR = "data/figure_4/spatial_clustering/sim_data"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dataset", choices=sorted(DATASETS), required=True)
    parser.add_argument("--input", type=Path, default=None, help="Override the auto-resolved h5ad path.")
    parser.add_argument("--strongmix", action="store_true",
                        help="Resolve --input to the strong-domain-mix sibling dataset "
                             f"({STRONGMIX_DIR}/<tag>_strongmix/) instead of the canonical weak-mix "
                             "Figure 2 h5ad, and tag outputs <dataset>_strongmix. Ignored if --input "
                             "is passed explicitly.")
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--run-tag", default=None,
                        help="Override the output subdir / filename stem (default: "
                             "<dataset>[_strongmix][_prebatch]). Set this to keep a weak-mix and a "
                             "strongmix --use-pre-batch run from colliding on <dataset>_prebatch when "
                             "--input is given by hand.")
    parser.add_argument("--slice-id", type=int, default=5, help="nnSVG has no 3D mode -- pick one slice.")
    parser.add_argument("--n-threads", type=int, default=4)
    parser.add_argument("--seed", type=int, default=42,
                        help="passed to run_nnsvg.R: set.seed() + MulticoreParam(RNGseed=) "
                             "so nnSVG/BRISC stochasticity is reproducible.")
    parser.add_argument("--use-pre-batch", action="store_true",
                         help="oracle run: swap X <- layers['counts_pre_batch'] before normalizing. "
                              "Mirrors 3D_scbsp.py's --use-pre-batch.")
    args = parser.parse_args()

    tag, file_modality = DATASETS[args.dataset]
    if args.input is None:
        if args.strongmix:
            args.input = (SIM_PAPER_DIR / STRONGMIX_DIR / f"{tag}_strongmix"
                          / f"simulation_{file_modality}_z_qc.h5ad")
        else:
            args.input = SIM_PAPER_DIR / "data" / "figure_2" / tag / f"simulation_{file_modality}_z_qc.h5ad"
    args.resolved_tag = tag + ("_strongmix" if args.strongmix else "")
    if args.run_tag is None:
        args.run_tag = (args.dataset + ("_strongmix" if args.strongmix else "")
                        + ("_prebatch" if args.use_pre_batch else ""))
    if args.output_dir is None:
        args.output_dir = SIM_PAPER_DIR / "data" / "figure_4" / "svg" / "nnsvg" / args.run_tag
    return args


def main() -> None:
    args = parse_args()
    if not args.input.is_file():
        raise SystemExit(f"Input not found: {args.input}\nRun the Figure 2 generate+QC pipeline first.")

    print(f"[load] {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"[load] AnnData shape: {adata.n_obs} x {adata.n_vars}")

    if args.use_pre_batch:
        if "counts_pre_batch" not in adata.layers:
            raise KeyError("--use-pre-batch: layers['counts_pre_batch'] not in this h5ad")
        adata.X = adata.layers["counts_pre_batch"].copy()
        print("[oracle] X <- counts_pre_batch (no batch effect)")

    adata = adata[adata.obs["slice_id"] == args.slice_id].copy()
    print(f"[slice] restricted to slice_id={args.slice_id}: {adata.n_obs} obs")

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    workdir = args.output_dir / "nnsvg_input"
    workdir.mkdir(parents=True, exist_ok=True)

    sio.mmwrite(str(workdir / "logcounts.mtx"), adata.X.T.tocoo())  # genes x cells
    (workdir / "cell_names.txt").write_text("\n".join(adata.obs_names))
    (workdir / "gene_names.txt").write_text("\n".join(adata.var_names))
    adata.var[["is_noise"]].to_csv(workdir / "gene_meta.csv")
    coords = adata.obsm["spatial"]
    pd.DataFrame(coords, columns=["x", "y"], index=adata.obs_names).to_csv(workdir / "locations.csv")
    print(f"[save] {workdir}")

    print(f"[nnsvg] running {RUN_NNSVG_R} n_threads={args.n_threads} seed={args.seed}")
    t0 = time.time()
    subprocess.run(
        [str(RSCRIPT), str(RUN_NNSVG_R), str(workdir), str(args.n_threads), str(args.seed)],
        check=True,
    )
    elapsed = time.time() - t0
    print(f"[nnsvg] done in {elapsed:.1f}s")

    res = pd.read_csv(workdir / "nnsvg_results.csv")
    gmeta = pd.read_csv(workdir / "gene_meta.csv")
    gmeta.columns = ["gene_names", "is_noise"]
    df = res.merge(gmeta, on="gene_names", how="left", validate="one_to_one")

    y_true = (~df["is_noise"]).astype(int)
    eps = np.finfo(np.float64).tiny
    score = -np.log10(df["pval"].to_numpy(dtype=float).clip(min=eps))
    auroc = roc_auc_score(y_true, score)
    auprc = average_precision_score(y_true, score)
    n_signal, n_noise = int(y_true.sum()), int((~y_true.astype(bool)).sum())
    print(f"[eval] n_signal={n_signal} n_noise={n_noise} AUROC={auroc:.4f} AUPRC={auprc:.4f}")

    genes_path = args.output_dir / f"gene_results_{args.run_tag}.csv"
    df.to_csv(genes_path, index=False)
    print(f"[save] {genes_path}")

    metrics = {
        "dataset": args.dataset,
        "run_tag": args.run_tag,
        "tag": args.resolved_tag,
        "strongmix": bool(args.strongmix),
        "input": str(args.input),
        "use_pre_batch": bool(args.use_pre_batch),
        "slice_id": args.slice_id,
        "n_obs": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_threads": args.n_threads,
        "seed": args.seed,
        "elapsed_sec": elapsed,
        "n_signal": n_signal,
        "n_noise": n_noise,
        "auroc": float(auroc),
        "auprc": float(auprc),
    }
    metrics_path = args.output_dir / f"metrics_{args.run_tag}.json"
    with open(metrics_path, "w") as f:
        json.dump(metrics, f, indent=2)
    print(f"[save] {metrics_path}")


if __name__ == "__main__":
    main()
