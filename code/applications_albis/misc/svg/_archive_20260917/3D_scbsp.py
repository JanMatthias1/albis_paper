#!/usr/bin/env python
"""
Figure 4A -- spatially variable gene (SVG) identification with scBSP
(https://github.com/YQ-Wang/scBSP) on an ALBIS z-stack.

Ground truth is already carried in every simulated h5ad's `.var` (see
SVG_PLAN.md): 56 `noise` genes have no true spatial signal by construction,
the remaining 500 (`marker_typeN` x8 + `marker_shared`) do, since cell types
are spatially organized into domains. This is exactly `is_noise` (False =
signal, True = noise) -- `marker_shared` counts as a true positive here (user
decision 2026-09-04: shared markers still track domain structure, just not
domain-specifically).

Runs scBSP's `granp()` on the FULL aligned 3D z-stack (`obsm['spatial_3d']`,
all 10 slices), not a single representative slice -- unlike Figure 2's
slice-5 convention, this uses the sim's true 3D coordinates directly (no
per-slice neighbor-graph staggering is involved, so there's no analogue of
the BANKSY cross-slice leakage bug), so pooling is the more informative
choice for genuinely testing spatial-signal recovery in 3D.

`granp()` returns a p-value per gene; lower = more spatially variable. Scored
as a binary-recovery problem (same shape as the ARI panel and the parked SVG
plan): AUROC/AUPRC of -log10(p) separating is_noise==False (signal, positive)
from is_noise==True (noise, negative).

Input: `input_exp_mat_raw` is passed the raw (integer, post-batch) counts
directly from `adata.X` as a CSR sparse matrix -- matches the parameter's own
naming and avoids a ~2-3GB dense conversion for the largest (cell, 600k obs)
modality. Gene names are reattached by position afterward (verified against a
DataFrame-input smoke test: identical p-values, positions preserved).

Expected environment:
    conda activate /dcs04/hicks/data/Jan/sim_project/sim_paper/env/albis-tutorial
    (scbsp installed 2026-09-04, no new env needed -- pure numpy/pandas/scipy/sklearn)

Example:
    python sim_paper/code/applications_albis/svg/3D_scbsp.py --dataset bin16um
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scbsp
from sklearn.metrics import average_precision_score, roc_auc_score

SCRIPT_DIR = Path(__file__).resolve().parent
SIM_PAPER_DIR = SCRIPT_DIR.parents[2]

# modality -> (Figure 2 canonical tag, modality name in the h5ad filename).
# bin16um / cell match Figure 2/3/4B/4C. The "spot" tag does NOT (see note).
#
# NOTE (2026-09-06): the "spot" entry is DELIBERATELY NOT the spot dataset
# used by Figure 2 / Figure 3 / Figure 4C. Those moved to
# `packing_pf0p04_log_mu_-2.0_theta_0.25_jitter0.10_bsigma03` (2026-09-05
# CytAssist probe-reference retune + 2026-09-06 `--theta-jitter 0.10`
# two-cloud fix). Figure 4A (and 4B/STAGATE) stay on the older
# `packing_pf0p04_log_mu_-2.5_bsigma03` spot config on purpose: the 4A/4B
# conclusions are about spatial-signal strength (domain-mix coupling), not
# real-data count-distribution fidelity. If 4A ever needs to match the rest
# of Figure 4, bump this "spot" tag + generate the strongmix sibling.
DATASETS = {
    "bin16um": ("packing_pf0p04_bin16um_log_mu_-2.5_bsigma07", "bin"),
    "spot": ("packing_pf0p04_log_mu_-2.5_bsigma03", "spot"),
    "cell": ("log_mu_-2.3_theta_0.40_jitter0.15_bsigma15", "cell"),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dataset", choices=sorted(DATASETS), required=True)
    parser.add_argument("--input", type=Path, default=None, help="Override the auto-resolved Figure 2 h5ad path.")
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--d1", type=float, default=1.0)
    parser.add_argument("--d2", type=float, default=3.0)
    parser.add_argument("--leaf-size", type=int, default=80)
    parser.add_argument("--use-gpu", action="store_true")
    parser.add_argument("--use-pre-batch", action="store_true",
                         help="oracle run: swap X <- layers['counts_pre_batch'] (the clean "
                              "pre-batch counts) instead of X (post-batch) -- bounds what scBSP "
                              "can recover with zero batch noise; writes to <dataset>_prebatch/. "
                              "Mirrors 3D_stagate.py's --use-pre-batch.")
    args = parser.parse_args()

    tag, file_modality = DATASETS[args.dataset]
    if args.input is None:
        args.input = SIM_PAPER_DIR / "data" / "figure_2" / tag / f"simulation_{file_modality}_z_qc.h5ad"
    args.run_tag = args.dataset + ("_prebatch" if args.use_pre_batch else "")
    if args.output_dir is None:
        args.output_dir = SIM_PAPER_DIR / "data" / "figure_4" / "svg" / "scbsp" / args.run_tag
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

    sp_mat = np.asarray(adata.obsm["spatial_3d"], dtype=np.float64)
    exp_mat = adata.X  # raw counts (post-batch by default, or pre-batch oracle above)

    print(f"[scbsp] d1={args.d1} d2={args.d2} leaf_size={args.leaf_size} use_gpu={args.use_gpu}")
    t0 = time.time()
    result = scbsp.granp(sp_mat, exp_mat, d1=args.d1, d2=args.d2, leaf_size=args.leaf_size, use_gpu=args.use_gpu)
    elapsed = time.time() - t0
    print(f"[scbsp] done in {elapsed:.1f}s, {len(result)} genes")

    assert len(result) == adata.n_vars, f"granp returned {len(result)} rows, expected {adata.n_vars}"
    result = result.reset_index(drop=True)
    result["gene_names"] = adata.var_names.to_numpy()  # positional reattach; sparse input drops real names

    var = adata.var.reset_index(names="gene_names")
    df = result.merge(var, on="gene_names", how="left", validate="one_to_one")

    # lower p -> more spatially variable -> higher score. Clip to avoid -inf on p=0.
    eps = np.finfo(np.float64).tiny
    df["neg_log10_p"] = -np.log10(df["p_values"].clip(lower=eps))

    y_true = (~df["is_noise"]).astype(int)  # 1 = true signal (incl. marker_shared), 0 = true noise
    score = df["neg_log10_p"].to_numpy()
    auroc = roc_auc_score(y_true, score)
    auprc = average_precision_score(y_true, score)
    n_signal, n_noise = int(y_true.sum()), int((~y_true.astype(bool)).sum())
    print(f"[eval] n_signal={n_signal} n_noise={n_noise} AUROC={auroc:.4f} AUPRC={auprc:.4f}")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    genes_path = args.output_dir / f"gene_pvalues_{args.run_tag}.csv"
    df.to_csv(genes_path, index=False)
    print(f"[save] {genes_path}")

    metrics = {
        "dataset": args.dataset,
        "run_tag": args.run_tag,
        "tag": DATASETS[args.dataset][0],
        "input": str(args.input),
        "use_pre_batch": bool(args.use_pre_batch),
        "n_obs": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "d1": args.d1,
        "d2": args.d2,
        "leaf_size": args.leaf_size,
        "use_gpu": args.use_gpu,
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
