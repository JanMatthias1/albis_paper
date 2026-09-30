"""ALBIS vs scCube with Figure 2's functions, on the scales where scCube is comparable.

scCube outputs log-normalized expression (log1p of normalize_total 1e4), not
counts, so Figure 2's count panels (mean-variance/dropout NB fits, total counts,
raw/normalized histograms) do not apply. This reuses count_distribution.py's
own functions for what does apply:
  - genes detected per observation (JSD) and sparsity: all modalities
  - log1p(normalized) histogram (JSD): cell only. ALBIS goes through the same
    normalize_total(1e4) -> log1p sampling sequence as Figure 2's raw_norm_log
    panel; scCube values are used unchanged.
  - bin/spot: scCube values are sums of cell log-expression, so no log-scale
    JSD; a descriptive native-expression panel instead.
"""
import argparse
import json
import sys
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "count_distribution"))
import count_distribution as cd  # noqa: E402

TARGET_SUM, SAMPLE_SIZE, RANDOM_STATE = 1e4, 2_000_000, 0  # count_distribution.py defaults


def load(path, label, display, color):
    a = ad.read_h5ad(path)
    n_genes = cd.genes_per_cell(a.X)
    return dict(label=label, display_label=display, color=color, adata=a, X=a.X,
                n_genes=n_genes, sparsity=cd.compute_sparsity_stats(a.X, n_genes))


def log_panel(albis, sccube, rng, output_path):
    # Same draw order as cd.plot_raw_norm_log_compare: ALBIS raw, norm, log; then the next dataset.
    tmp = ad.AnnData(X=albis["X"].copy())
    cd.sample_values(tmp.X.data.copy(), SAMPLE_SIZE, rng)
    sc.pp.normalize_total(tmp, target_sum=TARGET_SUM)
    cd.sample_values(tmp.X.data.copy(), SAMPLE_SIZE, rng)
    sc.pp.log1p(tmp)
    albis_log = cd.sample_values(tmp.X.data.copy(), SAMPLE_SIZE, rng)
    X = sccube["X"]
    sccube_log = cd.sample_values((X.data if hasattr(X, "data") else np.asarray(X).ravel()).copy(), SAMPLE_SIZE, rng)
    vals = [albis_log[albis_log > 0], sccube_log[sccube_log > 0]]
    all_vals = np.concatenate(vals)
    bins = np.linspace(all_vals.min(), all_vals.max(), 60)
    fig, ax = plt.subplots(figsize=(6, 6))
    for d, v in zip([albis, sccube], vals):
        ax.hist(v, bins=bins, density=True, color=d["color"], alpha=0.5, label=d["display_label"])
    jsd = cd.compute_jsd(vals[0], vals[1], bins)
    cd.annotate_jsd(ax, jsd, loc="upper right")
    ax.set_title("log1p(normalized)")
    ax.set_xlabel("Value (nonzero matrix entries)")
    ax.set_ylabel("Density")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return jsd


def native_panel(d, output_path):
    mean, var = cd.gene_mean_var(d["X"])
    total = np.asarray(d["X"].sum(axis=1)).ravel()
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    axes[0].scatter(mean, var, s=9, alpha=.5, color=d["color"])
    axes[0].set(xlabel="Mean native expression", ylabel="Variance of native expression", title="scCube native mean–variance")
    axes[1].hist(total, bins=60, color=d["color"], alpha=.7)
    axes[1].set(xlabel="Total native expression per observation", ylabel="Number of observations",
                title="scCube native expression totals")
    fig.suptitle("Sums of reconstructed cell log-expression")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--modality", choices=["cell", "bin", "spot"], required=True)
    parser.add_argument("--albis", type=Path, required=True)
    parser.add_argument("--sccube", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    cd.TITLE_CONTEXT = "ALBIS vs scCube"
    albis = load(args.albis, "ALBIS", "ALBIS", cd.PRIMARY_COLOR)
    sccube = load(args.sccube, "scCube", "scCube", cd.COMPARE_COLOR)
    specs = [albis, sccube]
    cd.plot_genes_per_cell_compare(specs, args.output_dir / "genes_per_cell_compare.png")
    cd.plot_sparsity_summary_compare(specs, args.output_dir / "sparsity_summary_compare.png")
    cd.plot_dataset_legend(specs, args.output_dir / "dataset_legend.png")
    positive = [d["n_genes"][d["n_genes"] > 0] for d in specs]
    all_pos = np.concatenate(positive)
    gbins = np.logspace(np.log10(all_pos.min()), np.log10(all_pos.max()), 60)  # as in cd.plot_genes_per_cell_compare
    jsd = {"genes_per_cell": cd.compute_jsd(positive[0], positive[1], gbins)}
    if args.modality == "cell":
        jsd["log1p_normalized"] = log_panel(albis, sccube, np.random.default_rng(RANDOM_STATE),
                                            args.output_dir / "log_normalized_compare.png")
    else:
        native_panel(sccube, args.output_dir / "sccube_native_expression.png")
    summary = dict(modality=args.modality, inputs={"ALBIS": str(args.albis), "scCube": str(args.sccube)},
                   n_obs={d["label"]: int(d["adata"].n_obs) for d in specs},
                   n_vars={d["label"]: int(d["adata"].n_vars) for d in specs},
                   genes_per_cell_median={d["label"]: float(np.median(d["n_genes"])) for d in specs},
                   matrix_zero_frac={d["label"]: d["sparsity"]["matrix_zero_frac"] for d in specs},
                   empty_row_frac={d["label"]: d["sparsity"]["empty_row_frac"] for d in specs},
                   jsd_albis_vs_sccube=jsd)
    (args.output_dir / "comparison_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
