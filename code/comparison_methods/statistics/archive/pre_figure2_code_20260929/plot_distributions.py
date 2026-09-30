#!/usr/bin/env python
"""Figure 5B: shared-tissue Figure 4 ALBIS and current Figure 5A competitors.

Uses displayed Slice 5 (zero-based slice_id=4) for cell/bin/spot. Both versions
remove explicitly empty/unassigned captures. --hvg-match additionally selects
556 genes by variance on normalized/log expression; raw X is retained for stats.
ALBIS has 556 genes; SPIDER/scCube have 2000. scCube expression is continuous.
This is a descriptive comparison, not a matched technical-noise benchmark.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import scanpy as sc

SIM_PAPER_DIR = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "count_distribution"))
import count_distribution as cd  # noqa: E402

SLICE_ID = 4
OVERVIEW_RUN = SIM_PAPER_DIR / "data/comparison_methods/figure_5A_600k"
FIGURE_5B_DIR = SIM_PAPER_DIR / "data/comparison_methods/figure_5B"
_COLORS = {"our_method": "#2c7fb8", "splatter_spider": "#d95f02", "splatter_sccube": "#7570b3"}
_DISPLAY_LABELS = {"our_method": "ALBIS", "splatter_spider": "SPIDER", "splatter_sccube": "scCube"}

def sources_for(modality: str) -> list[dict]:
    methods = {"our_method": "albis", "splatter_spider": "spider", "splatter_sccube": "sccube"}
    return [{"label": label, "display_label": _DISPLAY_LABELS[label], "color": _COLORS[label],
             "path": OVERVIEW_RUN / method / f"{modality}.h5ad"} for label, method in methods.items()]

TARGET_SUM = 1e4
SAMPLE_SIZE = 2_000_000
RANDOM_STATE = 0


def plot_mean_variance_compare_n(datasets: list[dict], output_path: Path) -> None:
    """Local copy of count_distribution.plot_mean_variance_compare with a title
    that fits a 3-way method comparison instead of the original's "vs real"."""
    fig, ax = plt.subplots(figsize=(8, 6.5))

    all_mean = np.concatenate([d["mean"][(d["mean"] > 0) & (d["var"] > 0)] for d in datasets])
    x = np.logspace(np.log10(all_mean.min()), np.log10(all_mean.max()), 200)
    ax.plot(x, x, "k--", label="Poisson (var = mean)", zorder=1)

    for d in datasets:
        keep = (d["mean"] > 0) & (d["var"] > 0)
        ax.scatter(d["mean"][keep], d["var"][keep], s=8, alpha=0.4, linewidths=0, color=d["color"], label=d["display_label"])
        if np.isfinite(d["theta"]):
            ax.plot(x, x + x**2 / d["theta"], "-", color=d["color"], linewidth=1.5,
                     label=rf"{d['display_label']} NB fit ($\hat\theta$={d['theta']:.1f})")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Mean count per gene")
    ax.set_ylabel("Variance per gene")
    ax.set_title("Gene mean-variance")
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_mean_dropout_compare_n(datasets: list[dict], output_path: Path) -> None:
    """Local copy of count_distribution.plot_mean_dropout_compare, same rationale."""
    fig, ax = plt.subplots(figsize=(8, 6.5))

    all_mean = np.concatenate([d["mean"][d["mean"] > 0] for d in datasets])
    x = np.logspace(np.log10(all_mean.min()), np.log10(all_mean.max()), 200)
    ax.plot(x, np.exp(-x), "k--", label="Poisson-predicted", zorder=1)

    for d in datasets:
        keep = d["mean"] > 0
        ax.scatter(d["mean"][keep], d["zero_frac"][keep], s=8, alpha=0.4, linewidths=0, color=d["color"], label=d["display_label"])
        if np.isfinite(d["theta"]):
            ax.plot(x, (d["theta"] / (d["theta"] + x)) ** d["theta"], "-", color=d["color"], linewidth=1.5,
                     label=rf"{d['display_label']} NB-predicted ($\hat\theta$={d['theta']:.1f})")

    ax.set_xscale("log")
    ax.set_xlabel("Mean count per gene")
    ax.set_ylabel("Fraction of zero observations")
    ax.set_title("Gene mean-dropout")
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_sparsity_summary_compare_n(datasets: list[dict], output_path: Path) -> None:
    """Local copy of count_distribution.plot_sparsity_summary_compare, same rationale."""
    metrics = [("matrix_zero_frac", "Zero matrix\nentries"), ("empty_row_frac", "Fully-empty\nobservations")]
    x = np.arange(len(metrics))
    width = 0.8 / len(datasets)

    fig, ax = plt.subplots(figsize=(6.5, 6))
    for i, d in enumerate(datasets):
        values = [d["sparsity"][key] * 100 for key, _ in metrics]
        offset = (i - (len(datasets) - 1) / 2) * width
        bars = ax.bar(x + offset, values, width=width, color=d["color"], label=d["display_label"])
        for bar, v in zip(bars, values):
            ax.annotate(f"{v:.1f}%", (bar.get_x() + bar.get_width() / 2, v), ha="center", va="bottom",
                        fontsize=cd.ANNOT_SIZE, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels([label for _, label in metrics])
    ax.set_ylabel("Percent (%)")
    ax.set_ylim(0, 105)
    ax.set_title("Sparsity summary")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_dataset_legend(datasets: list[dict], output_path: Path) -> None:
    """Standalone legend image (dataset display-name -> color patch), for
    hand-insertion into total_counts_compare.png / raw_norm_log_compare.png,
    whose inline legends are kept compact. Local copy: count_distribution.py's
    on-disk state doesn't currently have this helper."""
    handles = [Patch(facecolor=d["color"], alpha=0.5, label=d["display_label"]) for d in datasets]
    fig, ax = plt.subplots(figsize=(3.5, 0.6 + 0.45 * len(handles)))
    ax.axis("off")
    ax.legend(handles=handles, loc="center", frameon=False)
    fig.savefig(output_path, dpi=180, bbox_inches="tight", pad_inches=0.02, transparent=True)
    plt.close(fig)


def annotate_jsd_pairs(ax, pairs: list[tuple[str, float]], loc: str = "upper left") -> None:
    """pairs: [(competitor_display_label, jsd_vs_our_method), ...]. Same visual
    style as count_distribution.annotate_jsd, extended to multiple lines."""
    x, ha = (0.03, "left") if "left" in loc else (0.97, "right")
    y, va = (0.97, "top") if "upper" in loc else (0.03, "bottom")
    text = "\n".join(f"JSD(ALBIS vs {label}) = {jsd:.3f}" for label, jsd in pairs)
    ax.text(x, y, text, transform=ax.transAxes, ha=ha, va=va,
            fontsize=cd.ANNOT_SIZE, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="0.6", alpha=0.85))


def plot_total_counts_compare_n(datasets: list[dict], output_path: Path) -> None:
    positive = [d["total_counts"][d["total_counts"] > 0] for d in datasets]
    all_positive = np.concatenate(positive)
    bins = np.logspace(np.log10(all_positive.min()), np.log10(all_positive.max()), 60)

    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    for d, pos in zip(datasets, positive):
        ax.hist(pos, bins=bins, density=True, color=d["color"], alpha=0.5, label=d["display_label"])
    ax.set_xscale("log")
    ax.set_xlabel("Total counts per observation")
    ax.set_ylabel("Density")
    ax.set_title("Total counts per observation")
    ax.legend(frameon=False, loc="upper right")
    pairs = [(d["display_label"], cd.compute_jsd(positive[0], pos, bins))
             for d, pos in zip(datasets[1:], positive[1:])]
    annotate_jsd_pairs(ax, pairs, loc="upper left")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_genes_per_cell_compare_n(datasets: list[dict], output_path: Path) -> None:
    positive = [d["n_genes"][d["n_genes"] > 0] for d in datasets]
    all_positive = np.concatenate(positive)
    bins = np.logspace(np.log10(all_positive.min()), np.log10(all_positive.max()), 60)

    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    for d, pos in zip(datasets, positive):
        ax.hist(pos, bins=bins, density=True, color=d["color"], alpha=0.5, label=d["display_label"])
    ax.set_xscale("log")
    ax.set_xlabel("Genes with nonzero expression per observation")
    ax.set_ylabel("Density")
    ax.set_title("Genes with nonzero expression per observation")
    ax.legend(frameon=False, loc="upper right")
    pairs = [(d["display_label"], cd.compute_jsd(positive[0], pos, bins))
             for d, pos in zip(datasets[1:], positive[1:])]
    annotate_jsd_pairs(ax, pairs, loc="upper left")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_raw_norm_log_compare_n(datasets: list[dict], output_path: Path, rng: np.random.Generator, cell_scale: bool = False) -> None:
    staged = []
    for d in datasets:
        if d["label"] == "splatter_sccube":
            values = d["X"].data if hasattr(d["X"], "tocsr") else np.asarray(d["X"]).ravel()
            staged.append({"display_label": d["display_label"], "color": d["color"],
                           ("log" if cell_scale else "native"): cd.sample_values(values.copy(), SAMPLE_SIZE, rng)})
            continue
        tmp = sc.AnnData(X=d["X"].copy())
        raw_vals = cd.sample_values(tmp.X.data.copy() if hasattr(tmp.X, "data") else np.asarray(tmp.X).ravel(),
                                     SAMPLE_SIZE, rng)
        sc.pp.normalize_total(tmp, target_sum=TARGET_SUM)
        norm_vals = cd.sample_values(tmp.X.data.copy() if hasattr(tmp.X, "data") else np.asarray(tmp.X).ravel(),
                                      SAMPLE_SIZE, rng)
        sc.pp.log1p(tmp)
        log_vals = cd.sample_values(tmp.X.data.copy() if hasattr(tmp.X, "data") else np.asarray(tmp.X).ravel(),
                                     SAMPLE_SIZE, rng)
        staged.append({"display_label": d["display_label"], "color": d["color"],
                        "raw": raw_vals, "norm": norm_vals, "log": log_vals})

    fig, axes = plt.subplots(1, 3 if cell_scale else 4, figsize=(18 if cell_scale else 24, 6))
    stage_specs = [(axes[0], "raw", "Raw counts", True),
                   (axes[1], "norm", "Normalized", True),
                   (axes[2], "log", "log1p(normalized)", False)]
    if not cell_scale:
        stage_specs.append((axes[3], "native", "scCube: summed cell log-expression", False))
    for ax, key, title, log_x in stage_specs:
        active = [s for s in staged if key in s]
        vals_by_dataset = [s[key][s[key] > 0] for s in active]
        all_vals = np.concatenate(vals_by_dataset)
        if log_x:
            bins = np.logspace(np.log10(all_vals.min()), np.log10(all_vals.max()), 60)
            ax.set_xscale("log")
        else:
            bins = np.linspace(all_vals.min(), all_vals.max(), 60)
        for s, vals in zip(active, vals_by_dataset):
            ax.hist(vals, bins=bins, density=True, color=s["color"], alpha=0.5, label=s["display_label"])
        ax.set_title(title)
        ax.set_xlabel("Value (nonzero matrix entries)")
        pairs = [(s["display_label"], cd.compute_jsd(vals_by_dataset[0], vals, bins))
                 for s, vals in zip(active[1:], vals_by_dataset[1:])]
        if pairs:
            annotate_jsd_pairs(ax, pairs, loc="upper right")
    axes[0].set_ylabel("Density")
    axes[0].legend(frameon=False, loc="upper left")
    if cell_scale:
        axes[2].legend(frameon=False, loc="center right")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def apply_qc(adata, label: str):
    """Drop rows explicitly marked off-tissue/empty in the Figure 5A inputs. None of the three cell-level
    sources here carry is_empty or cell_type_true=="unassigned" -- those are
    only set by the bin/spot aggregation path -- so this is expected to be a
    documented no-op at cell resolution, not a silent assumption."""
    keep = np.ones(adata.n_obs, dtype=bool)
    if "is_empty" in adata.obs.columns:
        keep &= ~adata.obs["is_empty"].to_numpy()
    if "cell_type_true" in adata.obs.columns:
        keep &= (adata.obs["cell_type_true"].astype(str) != "unassigned").to_numpy()
    n_dropped = int((~keep).sum())
    print(f"[qc] {label}: dropped {n_dropped}/{adata.n_obs} off-tissue/empty rows")
    return adata[keep].copy() if n_dropped else adata


def hvg_match_panels(adatas: dict[str, "sc.AnnData"], cell_scale: bool = False) -> dict[str, "sc.AnnData"]:
    """HVG-subset every source with more genes than the smallest panel among
    them, down to that count -- controls for panel-size as a confound
    (Splatter's 2000-gene pool vs ALBIS's 556 native marker genes) rather
    than comparing distributions computed over very differently-sized gene
    sets. Ranks genes by variance on a normalize_total+log1p working copy and
    keeps the top `target`, computed directly (not scanpy's
    highly_variable_genes): "seurat_v3" needs the scikit-misc package this
    env doesn't have, and "seurat"/"cell_ranger" both bin genes by mean
    expression via pandas.cut, which errors on scCube's output (many
    near-identical near-zero means produce duplicate bin edges)."""
    target = min(a.n_vars for a in adatas.values())
    out = {}
    for label, adata in adatas.items():
        if adata.n_vars <= target:
            print(f"[hvg-match] {label}: already at/below target ({adata.n_vars} <= {target}), unchanged")
            out[label] = adata
            continue
        tmp = adata.copy()
        if label != "splatter_sccube":
            sc.pp.normalize_total(tmp, target_sum=TARGET_SUM)
            sc.pp.log1p(tmp)
        _, var = cd.gene_mean_var(tmp.X)
        top = np.argsort(var)[::-1][:target]
        out[label] = adata[:, np.sort(top)].copy()
        print(f"[hvg-match] {label}: {adata.n_vars} -> {out[label].n_vars} genes (target {target})")
    return out


def load_and_qc(spec: dict):
    print(f"[load] {spec['label']}: {spec['path']}")
    adata = sc.read_h5ad(spec["path"], backed="r")
    print(f"[load] {spec['label']} shape before slice restriction: {adata.n_obs} x {adata.n_vars}")
    if "slice_id" not in adata.obs.columns:
        raise SystemExit(f"{spec['label']}: 'slice_id' not found in obs columns")
    keep = adata.obs["slice_id"].astype(int) == SLICE_ID
    if not keep.any():
        raise SystemExit(f"{spec['label']}: slice_id={SLICE_ID} matches no observations "
                          f"(available: {sorted(adata.obs['slice_id'].unique().tolist())})")
    full = adata
    adata = full[keep].to_memory()
    full.file.close()
    print(f"[load] {spec['label']} restricted to slice_id={SLICE_ID}: {adata.n_obs} x {adata.n_vars}")
    return apply_qc(adata, spec["label"])


def compute_stats(spec: dict, adata) -> dict:
    X = adata.X
    mean, var = cd.gene_mean_var(X)
    theta = cd.fit_common_dispersion(mean, var) if spec["label"] != "splatter_sccube" else float("nan")
    n_genes = cd.genes_per_cell(X)
    out = dict(spec)
    out.update({
        "adata": adata, "X": X, "mean": mean, "var": var, "theta": theta,
        "zero_frac": cd.gene_zero_fraction(X, adata.n_obs),
        "total_counts": np.asarray(X.sum(axis=1)).ravel(),
        "n_genes": n_genes,
        "sparsity": cd.compute_sparsity_stats(X, n_genes),
    })
    print(f"[nb] {spec['label']}: theta_hat = {theta:.2f}")
    print(f"[sparsity] {spec['label']}: {out['sparsity']['matrix_zero_frac'] * 100:.2f}% zero entries, "
          f"{out['sparsity']['empty_row_frac'] * 100:.2f}% fully-empty cells")
    return out


def plot_sccube_native(d, output_dir, modality):
    """Describe native scCube output without a count-model interpretation."""
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    axes[0].scatter(d['mean'], d['var'], s=9, alpha=.5, color=d['color'])
    axes[0].set(xlabel='Mean native expression', ylabel='Variance of native expression', title='scCube native mean–variance')
    axes[1].hist(d['total_counts'], bins=60, color=d['color'], alpha=.7)
    axes[1].set(xlabel='Total native expression per observation', ylabel='Number of observations', title='scCube native expression totals')
    fig.suptitle('Reconstructed cell log-expression' if modality == 'cell' else 'Sums of reconstructed cell log-expression')
    fig.tight_layout()
    fig.savefig(output_dir/'sccube_native_expression.png',dpi=180,bbox_inches='tight')
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--modality", choices=["cell", "bin", "spot"], default="cell")
    parser.add_argument("--hvg-match", action="store_true",
                         help="Also apply QC-row-drop + HVG panel-size matching; "
                              "writes to distribution_comparison_qc_hvg/ instead.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    sources = sources_for(args.modality)
    for spec in sources:
        if not spec["path"].is_file():
            raise SystemExit(f"Missing input for {spec['label']}: {spec['path']}")
    output_dir = (FIGURE_5B_DIR / args.modality /
                  ("distribution_comparison_qc_hvg" if args.hvg_match else "distribution_comparison"))
    output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RANDOM_STATE)

    adatas = {spec["label"]: load_and_qc(spec) for spec in sources}
    if args.hvg_match:
        adatas = hvg_match_panels(adatas, cell_scale=args.modality == "cell")
    datasets = [compute_stats(spec, adatas[spec["label"]]) for spec in sources]

    counts = [d for d in datasets if d["label"] != "splatter_sccube"]
    plot_sccube_native(next(d for d in datasets if d["label"] == "splatter_sccube"), output_dir, args.modality)
    plot_mean_variance_compare_n(counts, output_dir / "mean_variance_compare.png")
    plot_mean_dropout_compare_n(counts, output_dir / "mean_dropout_compare.png")
    plot_sparsity_summary_compare_n(datasets, output_dir / "sparsity_summary_compare.png")
    plot_total_counts_compare_n(counts, output_dir / "total_counts_compare.png")
    plot_genes_per_cell_compare_n(datasets, output_dir / "genes_per_cell_compare.png")
    plot_raw_norm_log_compare_n(datasets, output_dir / "raw_norm_log_compare.png", rng, cell_scale=args.modality == "cell")
    plot_dataset_legend(datasets, output_dir / "dataset_legend.png")
    print(f"[save] Comparison plots written to {output_dir}")

    summary = {
        "modality": args.modality,
        "slice_id": SLICE_ID,
        "hvg_match": bool(args.hvg_match),
        "sources": {d["label"]: str(d["path"].relative_to(SIM_PAPER_DIR)) for d in datasets},
        "figure": "5B",
        "cell_scale_handling": "scCube unchanged, shown only in log panel; native variance for scCube gene selection" if args.modality == "cell" else "scCube unchanged sums of cell log-expression; separate native panel; no count-based NB fit or total-count comparison",
        "slice_number": SLICE_ID + 1,
        "resolved_sources": {d["label"]: str(d["path"].resolve()) for d in datasets},
        "caveat": "All methods use Figure 5A inputs with 600k source cells each. ALBIS uses the shared Figure 4 strong_domain_mix_shift3x tissue across cell/bin16um/spot. Expression/noise models and capture counts differ. Both versions remove explicitly empty/unassigned rows. scCube X is continuous. " + ("Panels matched by normalized-log variance selection; genes are not homologous across methods." if args.hvg_match else "Native gene panels: ALBIS 556; SPIDER/scCube 2000."),
        "n_obs": {d["label"]: int(d["adata"].n_obs) for d in datasets},
        "n_vars": {d["label"]: int(d["adata"].n_vars) for d in datasets},
        "theta_hat": {d["label"]: d["theta"] if np.isfinite(d["theta"]) else None for d in datasets},
        "total_counts_median": {d["label"]: float(np.median(d["total_counts"])) for d in counts},
        "native_total_expression_median": {d["label"]: float(np.median(d["total_counts"])) for d in datasets},
        "total_counts_mean": {d["label"]: float(np.mean(d["total_counts"])) for d in counts},
        "genes_per_cell_median": {d["label"]: float(np.median(d["n_genes"])) for d in datasets},
        "matrix_zero_frac": {d["label"]: d["sparsity"]["matrix_zero_frac"] for d in datasets},
        "empty_row_frac": {d["label"]: d["sparsity"]["empty_row_frac"] for d in datasets},
    }
    summary_path = output_dir / "comparison_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    print(f"[save] Numeric summary written to {summary_path}")


if __name__ == "__main__":
    main()
