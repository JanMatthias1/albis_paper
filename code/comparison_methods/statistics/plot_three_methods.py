"""Figure 5B: ALBIS, SPIDER and scCube on the same panels, with Figure 2's statistics.

All statistics come from code/count_distribution/count_distribution.py (gene mean/variance,
NB dispersion fit, dropout, totals, genes detected, sparsity, value sampling, JSD on shared
bins) with its defaults: .X counts, 2,000,000 sampled values, random state 0, target sum 1e4.

scCube returns log-normalized expression (log1p of normalize_total 1e4), not counts, so it
is shown only where it is comparable:
  - genes detected per observation and sparsity: all three methods, every modality
  - log1p(normalized) histogram: cell only (scCube values unchanged; ALBIS and SPIDER go
    through Figure 2's normalize_total -> log1p sequence)
  - bin/spot: scCube values are sums of cell log-expression, so they get their own fourth
    panel in raw_norm_log_compare.png (unchanged values, no JSD) plus a separate
    native-expression figure.
Count panels (mean-variance, mean-dropout, total counts, raw and normalized values) show
ALBIS and SPIDER. JSDs are reported for ALBIS vs each other method.
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
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "count_distribution"))
import count_distribution as cd  # noqa: E402  (also applies Figure 2's font sizes)

TARGET_SUM, SAMPLE_SIZE, RANDOM_STATE = 1e4, 2_000_000, 0  # count_distribution.py defaults
METHODS = [("albis", "ALBIS", "#2c7fb8"), ("spider", "SPIDER", "#d95f02"), ("sccube", "scCube", "#7570b3")]


def load(path, label, color):
    a = ad.read_h5ad(path)
    X = a.X
    d = dict(label=label, display_label=label, color=color, adata=a, X=X)
    d["n_genes"] = cd.genes_per_cell(X)
    d["sparsity"] = cd.compute_sparsity_stats(X, d["n_genes"])
    return d


def add_count_stats(d):
    d["mean"], d["var"] = cd.gene_mean_var(d["X"])
    d["theta"] = cd.fit_common_dispersion(d["mean"], d["var"])
    d["zero_frac"] = cd.gene_zero_fraction(d["X"], d["adata"].n_obs)
    d["total_counts"] = np.asarray(d["X"].sum(axis=1)).ravel()


def ref_label(d):
    """Name of the JSD reference in the annotation box (optional short "jsd_label")."""
    return d.get("jsd_label", d["display_label"])


def annotate(ax, jsds, ref="ALBIS", loc="upper left"):
    x, ha = (0.03, "left") if "left" in loc else (0.97, "right")
    ax.text(x, 0.97, "\n".join(f"JSD({ref} vs {k}) = {v:.3f}" for k, v in jsds.items()),
            transform=ax.transAxes, ha=ha, va="top", fontsize=cd.ANNOT_SIZE, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="0.6", alpha=0.85))


def grid(values, log_x):
    """Figure 2's 60 shared bins over the pooled values."""
    pooled = np.concatenate(values)
    return (np.logspace(np.log10(pooled.min()), np.log10(pooled.max()), 60) if log_x
            else np.linspace(pooled.min(), pooled.max(), 60))


def hist_panel(ax, values, datasets, log_x):
    """Overlaid histograms; returns JSD of the first dataset (ALBIS) vs each other one. Each JSD
    uses bins pooled over that pair only, exactly as count_distribution.py compares two datasets,
    so ALBIS vs SPIDER matches Figure 2's definition regardless of scCube being drawn.
    Log-x panels show the fraction per bin: density divides by the linear bin width, which on
    log bins shrinks high-value distributions by orders of magnitude (e.g. SPIDER totals)."""
    bins = grid(values, log_x)
    for d, v in zip(datasets, values):
        weights = np.full(len(v), 1 / len(v)) if log_x else None
        ax.hist(v, bins=bins, density=not log_x, weights=weights, color=d["color"], alpha=0.5,
                label=d["display_label"])
    ax.set_ylabel("Fraction per bin" if log_x else "Density")
    if log_x:
        ax.set_xscale("log")
    ax.set_ylim(top=ax.get_ylim()[1] * 1.22)  # headroom so the JSD box sits above the bars
    return {d["display_label"]: cd.compute_jsd(values[0], v, grid([values[0], v], log_x))
            for d, v in zip(datasets[1:], values[1:])}


def save(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def mean_variance(counts, out):
    fig, ax = plt.subplots(figsize=(8, 6.5))
    pooled = np.concatenate([d["mean"][(d["mean"] > 0) & (d["var"] > 0)] for d in counts])
    x = np.logspace(np.log10(pooled.min()), np.log10(pooled.max()), 200)
    ax.plot(x, x, "k--", label="Poisson (var = mean)", zorder=1)
    for d in counts:
        keep = (d["mean"] > 0) & (d["var"] > 0)
        ax.scatter(d["mean"][keep], d["var"][keep], s=8, alpha=0.4, linewidths=0, color=d["color"], label=d["display_label"])
        if np.isfinite(d["theta"]):
            ax.plot(x, x + x ** 2 / d["theta"], "-", color=d["color"], linewidth=1.5,
                    label=rf"{d['display_label']} NB fit ($\hat\theta$={d['theta']:.1f})")
    ax.set(xscale="log", yscale="log", xlabel="Mean count per gene", ylabel="Variance per gene",
           title="Gene mean-variance")
    ax.legend(frameon=False, loc="upper left")
    save(fig, out / "mean_variance_compare.png")


def mean_dropout(counts, out):
    fig, ax = plt.subplots(figsize=(8, 6.5))
    pooled = np.concatenate([d["mean"][d["mean"] > 0] for d in counts])
    x = np.logspace(np.log10(pooled.min()), np.log10(pooled.max()), 200)
    ax.plot(x, np.exp(-x), "k--", label="Poisson-predicted", zorder=1)
    for d in counts:
        keep = d["mean"] > 0
        ax.scatter(d["mean"][keep], d["zero_frac"][keep], s=8, alpha=0.4, linewidths=0, color=d["color"], label=d["display_label"])
        if np.isfinite(d["theta"]):
            ax.plot(x, (d["theta"] / (d["theta"] + x)) ** d["theta"], "-", color=d["color"], linewidth=1.5,
                    label=rf"{d['display_label']} NB-predicted ($\hat\theta$={d['theta']:.1f})")
    ax.set(xscale="log", xlabel="Mean count per gene", ylabel="Fraction of zero observations",
           title="Gene mean-dropout")
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0)
    save(fig, out / "mean_dropout_compare.png")


def total_counts(counts, out):
    fig, ax = plt.subplots(figsize=(7, 6.5))
    jsd = hist_panel(ax, [d["total_counts"][d["total_counts"] > 0] for d in counts], counts, log_x=True)
    ax.set(xlabel="Total counts per observation", title="Total counts per observation")
    ax.legend(frameon=False, loc="upper right")
    annotate(ax, jsd, ref_label(counts[0]))
    save(fig, out / "total_counts_compare.png")
    return jsd


def genes_detected(datasets, out):
    fig, ax = plt.subplots(figsize=(7, 6.5))
    jsd = hist_panel(ax, [d["n_genes"][d["n_genes"] > 0] for d in datasets], datasets, log_x=True)
    ax.set(xlabel="Genes detected per observation", title="Genes detected per observation")
    ax.legend(frameon=False, loc="upper right")
    annotate(ax, jsd, ref_label(datasets[0]))
    save(fig, out / "genes_per_cell_compare.png")
    return jsd


def sparsity(datasets, out):
    metrics = [("matrix_zero_frac", "Zero matrix\nentries"), ("empty_row_frac", "Fully-empty\nobservations")]
    x, width = np.arange(len(metrics)), 0.8 / len(datasets)
    fig, ax = plt.subplots(figsize=(7, 6))
    for i, d in enumerate(datasets):
        values = [d["sparsity"][k] * 100 for k, _ in metrics]
        bars = ax.bar(x + (i - (len(datasets) - 1) / 2) * width, values, width=width, color=d["color"], label=d["display_label"])
        for bar, v in zip(bars, values):
            ax.annotate(f"{v:.1f}%", (bar.get_x() + bar.get_width() / 2, v), ha="center", va="bottom",
                        fontsize=cd.ANNOT_SIZE - 2, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels([label for _, label in metrics])
    ax.set(ylabel="Percent (%)", ylim=(0, 105), title="Sparsity summary")
    ax.legend(frameon=False)
    save(fig, out / "sparsity_summary_compare.png")


def raw_norm_log(counts, sccube, rng, out, native_panel=False):
    """Figure 2's raw -> normalized -> log1p sampling for ALBIS then SPIDER (same RNG order as
    count_distribution.py), then scCube's own values: on the log1p panel (cell), or on a fourth
    panel of its own (bin/spot, native_panel=True) since there they are sums of cell log1p values."""
    staged = []
    for d in counts:
        tmp = ad.AnnData(X=d["X"].copy())
        raw = cd.sample_values(tmp.X.data.copy(), SAMPLE_SIZE, rng)
        sc.pp.normalize_total(tmp, target_sum=TARGET_SUM)
        norm = cd.sample_values(tmp.X.data.copy(), SAMPLE_SIZE, rng)
        sc.pp.log1p(tmp)
        staged.append(dict(raw=raw, norm=norm, log=cd.sample_values(tmp.X.data.copy(), SAMPLE_SIZE, rng)))
    log_sets, log_vals = list(counts), [s["log"] for s in staged]
    if sccube is not None:
        X = sccube["X"]
        native = cd.sample_values((X.data if hasattr(X, "data") else np.asarray(X).ravel()).copy(), SAMPLE_SIZE, rng)
        if not native_panel:
            log_sets.append(sccube)
            log_vals.append(native)
    fig, axes = plt.subplots(1, 4 if native_panel else 3, figsize=(22.5 if native_panel else 17, 6))
    jsds = {}
    for ax, key, title, log_x in [(axes[0], "raw", "Raw counts", True), (axes[1], "norm", "Normalized", True)]:
        jsds[key] = hist_panel(ax, [s[key][s[key] > 0] for s in staged], counts, log_x)
        ax.set(title=title, xlabel="Value (nonzero matrix entries)")
        annotate(ax, jsds[key], ref_label(counts[0]), loc="upper right")
    jsds["log1p_normalized"] = hist_panel(axes[2], [v[v > 0] for v in log_vals], log_sets, log_x=False)
    axes[2].set(title="log1p(normalized)", xlabel="Value (nonzero matrix entries)")
    annotate(axes[2], jsds["log1p_normalized"], ref_label(counts[0]), loc="upper right")
    axes[2].legend(frameon=False, loc="center right")
    if native_panel:
        axes[3].hist(native[native > 0], bins=60, density=True, color=sccube["color"], alpha=0.5,
                     label=sccube["display_label"])
        axes[3].set(title="scCube native\n(sum of cell log1p(normalized))", xlabel="Value (nonzero matrix entries)",
                    ylabel="Density")
        axes[3].legend(frameon=False, loc="upper right")
    save(fig, out / "raw_norm_log_compare.png")
    return jsds


def sccube_native(d, out):
    mean, var = cd.gene_mean_var(d["X"])
    total = np.asarray(d["X"].sum(axis=1)).ravel()
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    axes[0].scatter(mean, var, s=9, alpha=.5, color=d["color"])
    axes[0].set(xlabel="Mean native expression", ylabel="Variance of native expression", title="scCube native mean–variance")
    axes[1].hist(total, bins=60, color=d["color"], alpha=.7)
    axes[1].set(xlabel="Total native expression per observation", ylabel="Number of observations",
                title="scCube native expression totals")
    fig.suptitle("Sums of cell log-expression")
    save(fig, out / "sccube_native_expression.png")


def legend(datasets, out):
    fig, ax = plt.subplots(figsize=(3.5, 0.6 + 0.45 * len(datasets)))
    ax.axis("off")
    ax.legend(handles=[Patch(facecolor=d["color"], alpha=0.5, label=d["display_label"]) for d in datasets],
              loc="center", frameon=False)
    fig.savefig(out / "dataset_legend.png", dpi=180, bbox_inches="tight", pad_inches=0.02, transparent=True)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--modality", choices=["cell", "bin", "spot"], required=True)
    for method, _, _ in METHODS:
        parser.add_argument(f"--{method}", type=Path, required=True, help=f"{method} slice h5ad")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    albis, spider, sccube = [load(getattr(args, m), label, color) for m, label, color in METHODS]
    counts, everyone = [albis, spider], [albis, spider, sccube]
    for d in counts:
        add_count_stats(d)
    mean_variance(counts, out)
    mean_dropout(counts, out)
    jsd = dict(total_counts=total_counts(counts, out), genes_per_cell=genes_detected(everyone, out))
    sparsity(everyone, out)
    jsd.update(raw_norm_log(counts, sccube, np.random.default_rng(RANDOM_STATE), out,
                            native_panel=args.modality != "cell"))
    if args.modality != "cell":
        sccube_native(sccube, out)
    legend(everyone, out)
    summary = dict(
        modality=args.modality, inputs={d["label"]: str(getattr(args, m)) for (m, _, _), d in zip(METHODS, everyone)},
        n_obs={d["label"]: int(d["adata"].n_obs) for d in everyone},
        n_vars={d["label"]: int(d["adata"].n_vars) for d in everyone},
        theta_hat={d["label"]: d["theta"] for d in counts},
        total_counts_median={d["label"]: float(np.median(d["total_counts"])) for d in counts},
        genes_per_cell_median={d["label"]: float(np.median(d["n_genes"])) for d in everyone},
        matrix_zero_frac={d["label"]: d["sparsity"]["matrix_zero_frac"] for d in everyone},
        empty_row_frac={d["label"]: d["sparsity"]["empty_row_frac"] for d in everyone},
        jsd_albis_vs=jsd)
    (out / "comparison_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
