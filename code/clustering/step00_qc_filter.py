#!/usr/bin/env python
"""
Drop degenerate observations (empty or near-empty cells/bins/spots) from a
simulated dataset before clustering.

This isn't about matching real data's full QC pipeline (see
code/real_data_qc/visium_hd_qc.py's SpotSweeper local-outlier QC for that) --
it's a minimal, targeted fix for a specific clustering failure: at
packing_pf0p04, bin-level Leiden never collapsed below ~190 clusters even at
the resolution-search floor, far short of the target 6-8. A likely cause is
graph fragmentation from empty/near-empty bins (9.1% of bins are still fully
empty at this packing fraction) -- their near-identical, degenerate
embeddings can form many tiny disconnected components that no resolution
merges into the main structure.

Filters (computed on the POST-batch count matrix -- see note at the filter
site; this is the matrix every downstream step consumes):
  - total_counts > 0 (drop fully-empty observations)
  - n_genes_detected >= --min-genes (drop near-empty observations with too
    few nonzero genes to carry real signal)

Expected environment:
    conda activate /dcs04/hicks/data/Jan/sim_project/albis/env/albis-tutorial

Example:
    python albis_paper/code/clustering/step00_qc_filter.py --modality bin --packing-tag packing_pf0p04
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import scanpy as sc
from scipy import sparse

SCRIPT_DIR = Path(__file__).resolve().parent
SIM_PAPER_DIR = SCRIPT_DIR.parents[1]
VALID_MODALITIES = ("spot", "bin", "cell")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Drop empty/near-empty observations from a simulated dataset before clustering."
    )
    parser.add_argument("--modality", choices=VALID_MODALITIES, required=True)
    parser.add_argument(
        "--packing-tag", default="packing_pf0p04",
        help="Reads data/noisy/<packing-tag>/simulation_<modality>_z.h5ad as input, unless --input is given.",
    )
    parser.add_argument("--input", type=Path, default=None)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument(
        "--min-genes", type=int, default=3,
        help="Minimum number of distinct genes detected (nonzero) to keep an observation.",
    )
    args = parser.parse_args()

    root = SIM_PAPER_DIR / "data" / "noisy" / args.packing_tag
    if args.input is None:
        args.input = root / f"simulation_{args.modality}_z.h5ad"
    if args.output is None:
        args.output = root / f"simulation_{args.modality}_z_qc.h5ad"
    return args


def main() -> None:
    args = parse_args()
    if not args.input.is_file():
        raise SystemExit(f"Input not found: {args.input}")

    print(f"[load] {args.input}")
    adata = sc.read_h5ad(args.input)
    n_before = adata.n_obs
    print(f"[load] AnnData shape: {n_before} x {adata.n_vars}")

    # Filter on the POST-batch counts -- the matrix every downstream step
    # Filter on the POST-batch counts -- the matrix every downstream step
    # actually consumes (count_distribution.py, step01_pca_harmony.py, STAIR). Using
    # counts_pre_batch here let near-empty bins pass the >=min-genes floor on
    # their pre-batch counts and then get zeroed by the batch noise, entering
    # PCA as empty rows (the detached-island artifact in the Fig 3 UMAPs).
    # In the raw noisy h5ad, adata.X is the post-batch integer counts; the
    # counts_pre_batch layer is the only pre-batch store.
    X = adata.layers["counts"] if "counts" in adata.layers else adata.X
    n_genes = np.asarray(X.getnnz(axis=1) if sparse.issparse(X) else np.count_nonzero(X, axis=1)).ravel()
    total_counts = np.asarray(X.sum(axis=1)).ravel()

    keep = (total_counts > 0) & (n_genes >= args.min_genes)
    n_empty = int((total_counts == 0).sum())
    n_near_empty = int(((total_counts > 0) & (n_genes < args.min_genes)).sum())

    print(f"[qc] fully-empty (total_counts == 0): {n_empty} ({n_empty / n_before * 100:.2f}%)")
    print(f"[qc] near-empty (0 < total_counts, n_genes < {args.min_genes}): "
          f"{n_near_empty} ({n_near_empty / n_before * 100:.2f}%)")

    adata = adata[keep].copy()
    n_after = adata.n_obs
    print(f"[qc] dropped {n_before - n_after} observations ({(n_before - n_after) / n_before * 100:.2f}%); "
          f"{n_after} remain")

    summary = {
        "modality": args.modality,
        "input": str(args.input),
        "min_genes": args.min_genes,
        "n_obs_before_qc": n_before,
        "n_fully_empty": n_empty,
        "n_near_empty": n_near_empty,
        "n_obs_after_qc": n_after,
    }
    summary_path = args.output.with_name(args.output.stem + "_summary.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[save] {summary_path}")

    adata.write_h5ad(args.output)
    print(f"[save] {args.output}")


if __name__ == "__main__":
    main()
