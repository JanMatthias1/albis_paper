"""One real Xenium sample vs ALBIS, SPIDER and scCube (Figure 5B cell inputs), HVG-matched.

Panels and statistics are Figure 5B's (plot_three_methods.py, i.e. Figure 2's
count_distribution.py functions); the real sample takes ALBIS's place as the first dataset, so
every JSD is real vs each method.

Panel matching follows Figure 2's qc_and_hvg_matched mode: each simulated dataset is cut down
to the real sample's gene count by keeping its top highly variable genes. ALBIS and SPIDER are
counts and use seurat_v3, as in Figure 2. scCube returns log1p(normalized) expression, so it
uses flavor="seurat" (made for log data); its values are not changed. The real sample keeps its
whole panel.
"""
import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc

import plot_three_methods as ptm
from plot_three_methods import cd

REAL_COLOR = "#4d4d4d"


def hvg_subset(a, n_top, flavor):
    hvg = ad.AnnData(X=a.X.copy(), var=a.var[[]].copy())
    sc.pp.highly_variable_genes(hvg, n_top_genes=n_top, flavor=flavor)
    return a[:, hvg.var["highly_variable"].to_numpy()].copy()


def as_dataset(a, label, color):
    d = dict(label=label, display_label=label, color=color, adata=a, X=a.X)
    d["n_genes"] = cd.genes_per_cell(a.X)
    d["sparsity"] = cd.compute_sparsity_stats(a.X, d["n_genes"])
    return d


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--real", type=Path, required=True, help="QC'd real Xenium h5ad (raw counts)")
    parser.add_argument("--real-label", required=True)
    for method, _, _ in ptm.METHODS:
        parser.add_argument(f"--{method}", type=Path, required=True, help=f"{method} cell slice h5ad")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)

    real = ad.read_h5ad(args.real)
    n_genes = real.n_vars
    sims = []
    for method, label, color in ptm.METHODS:
        a = ad.read_h5ad(getattr(args, method))
        flavor = "seurat" if method == "sccube" else "seurat_v3"
        print(f"[hvg] {label}: {a.n_vars} -> top {n_genes} genes ({flavor})")
        sims.append(as_dataset(hvg_subset(a, n_genes, flavor), label, color))
    albis, spider, sccube = sims
    real = as_dataset(real, args.real_label, REAL_COLOR)
    real["jsd_label"] = "Xenium"  # short, so the JSD box clears the legend

    counts, everyone = [real, albis, spider], [real, albis, spider, sccube]
    for d in counts:
        ptm.add_count_stats(d)
    ptm.mean_variance(counts, out)
    ptm.mean_dropout(counts, out)
    jsd = dict(total_counts=ptm.total_counts(counts, out), genes_per_cell=ptm.genes_detected(everyone, out))
    ptm.sparsity(everyone, out)
    jsd.update(ptm.raw_norm_log(counts, sccube, np.random.default_rng(ptm.RANDOM_STATE), out))
    ptm.legend(everyone, out)

    summary = dict(
        real=args.real_label, inputs={"real": str(args.real), **{m: str(getattr(args, m)) for m, _, _ in ptm.METHODS}},
        hvg_matched_to=n_genes,
        n_obs={d["label"]: int(d["adata"].n_obs) for d in everyone},
        n_vars={d["label"]: int(d["adata"].n_vars) for d in everyone},
        theta_hat={d["label"]: d["theta"] for d in counts},
        total_counts_median={d["label"]: float(np.median(d["total_counts"])) for d in counts},
        genes_per_cell_median={d["label"]: float(np.median(d["n_genes"])) for d in everyone},
        matrix_zero_frac={d["label"]: d["sparsity"]["matrix_zero_frac"] for d in everyone},
        empty_row_frac={d["label"]: d["sparsity"]["empty_row_frac"] for d in everyone},
        jsd_real_vs=jsd)
    (out / "comparison_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
