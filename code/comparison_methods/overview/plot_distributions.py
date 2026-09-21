#!/usr/bin/env python
"""Statistical-distribution comparison: our_method vs Splatter->Spider vs
Splatter->scCube, at slice_id=5, for cell/bin/spot (--modality).

Unlike code/comparison_methods/figure5's controllability/adherence benchmark
(which holds all three workflows to one shared high-level contract), this
script deliberately mixes sources:

  our_method (ALBIS)  -- data/figure_2/smaller_sphere/data/<canonical config>/
                          simulation_<modality>_z_qc.h5ad
                          The canonical, manuscript-tuned Figure 2 config for
                          this modality (dispersion calibrated against real
                          reference data -- see MODALITY_SOURCES for which
                          config per modality). Real batch effects are baked
                          into X.
  splatter_spider,
  splatter_sccube    -- data/comparison_methods/overview_full_35831924/
                          {spider,sccube}/<modality>.h5ad
                          The existing qualitative-overview run (generic
                          contract: n_cells=10000, 2000-gene Splatter pool,
                          8 equal-proportion types, no batch-effect model).

Rationale (see conversation, 2026-09-21): our_method's Figure 2 parameters
are already optimized/validated against real reference data, so reuse that
output directly rather than regenerating a fresh generic-contract ALBIS run
-- and reuse the overview run's existing Spider/scCube outputs rather than
regenerating those too. This is an illustrative comparison, not a
contract-controlled benchmark: n_cells, n_genes, and technical-noise model
differ across the three inputs, so treat differences in scale (not just
shape) with that in mind. scCube's X is continuous VAE output, not integer
counts -- treated as-is, matching code/common/targets.py's convention of
never rounding it.

Reuses the mean-variance/mean-dropout/sparsity plotting + NB-fit machinery
from code/count_distribution/count_distribution.py (imported, not
reimplemented) via sys.path. Its total_counts/genes_per_cell/raw_norm_log
"_compare" variants are hardcoded to exactly two datasets (for one JSD
annotation); this script's local N-way versions instead annotate JSD of
each competitor against our_method (the anchor).

Run with the comparison_methods analysis env's python (has scanpy):
    comparison_methods/env/analysis/bin/python \\
        sim_paper/code/comparison_methods/overview/plot_distributions.py \\
        --modality {cell,bin,spot}

Pass --hvg-match to additionally: (1) explicitly drop off-tissue/empty rows
(is_empty=True or cell_type_true="unassigned", mirroring Figure 2's QC
convention -- a no-op at cell resolution for all three sources (none of them
carry that concept at native-cell granularity), but a REAL filter at bin/spot
resolution -- e.g. splatter_spider/bin.h5ad is ~50% is_empty at slice 5);
and (2) HVG-subset (variance-ranked on a normalize_total+log1p working copy;
see hvg_match_panels docstring for why not scanpy's highly_variable_genes)
whichever source(s) have more genes than the smallest panel (556, our_method's
native marker-gene count) down to that count -- controlling for the
panel-size confound behind the raw comparison's total-counts/sparsity gaps
(Splatter's 2000-gene pool vs ALBIS's 556). Written to
distribution_comparison_qc_hvg/ instead, so the unmatched raw comparison in
distribution_comparison/ is left untouched.

Output: data/comparison_methods/figure_4A/<modality>/distribution_comparison[_qc_hvg]/
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

SLICE_ID = 5
FIGURE2_DATA = SIM_PAPER_DIR / "data/figure_2/smaller_sphere/data"
OVERVIEW_RUN = SIM_PAPER_DIR / "data/comparison_methods/overview_full_35831924"
FIGURE_4A_DIR = SIM_PAPER_DIR / "data/comparison_methods/figure_4A"

# Match code/figures/figure5{a,b}.py's workflow palette for visual consistency
# across the comparison_methods project.
_COLORS = {"our_method": "#2c7fb8", "splatter_spider": "#d95f02", "splatter_sccube": "#7570b3"}
_DISPLAY_LABELS = {"our_method": "Our method", "splatter_spider": "Splatter → Spider",
                    "splatter_sccube": "Splatter → scCube"}

# our_method path per modality: Figure 2's canonical config for that modality
# (see data/figure_2/smaller_sphere/data/README.md) -- bin uses the 16um
# config to match the overview run's bin_width_um=16.
MODALITY_OUR_METHOD_PATH = {
    "cell": FIGURE2_DATA / "log_mu_-2.5_theta_0.40_jitter0.15_bsigma15/simulation_cell_z_qc.h5ad",
    "bin": FIGURE2_DATA / "packing_pf0p04_bin16um_log_mu_-2.5_bsigma07/simulation_bin_z_qc.h5ad",
    "spot": FIGURE2_DATA / "packing_pf0p04_log_mu_-2.0_theta_0.25_jitter0.10_bsigma03/simulation_spot_z_qc.h5ad",
}


def sources_for(modality: str) -> list[dict]:
    paths = {
        "our_method": MODALITY_OUR_METHOD_PATH[modality],
        "splatter_spider": OVERVIEW_RUN / f"spider/{modality}.h5ad",
        "splatter_sccube": OVERVIEW_RUN / f"sccube/{modality}.h5ad",
    }
    return [{"label": label, "display_label": _DISPLAY_LABELS[label], "color": _COLORS[label], "path": path}
            for label, path in paths.items()]

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
    ax.set_ylabel("Fraction of zero cells")
    ax.set_title("Gene mean-dropout")
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_sparsity_summary_compare_n(datasets: list[dict], output_path: Path) -> None:
    """Local copy of count_distribution.plot_sparsity_summary_compare, same rationale."""
    metrics = [("matrix_zero_frac", "Zero matrix\nentries"), ("empty_row_frac", "Fully-empty\ncells")]
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
    text = "\n".join(f"JSD(ours, {label}) = {jsd:.3f}" for label, jsd in pairs)
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
    ax.set_xlabel("Total counts per cell")
    ax.set_ylabel("Density (fraction of cells)")
    ax.set_title("Total counts per cell")
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
    ax.set_xlabel("Genes detected per cell")
    ax.set_ylabel("Density (fraction of cells)")
    ax.set_title("Genes detected per cell")
    ax.legend(frameon=False, loc="upper right")
    pairs = [(d["display_label"], cd.compute_jsd(positive[0], pos, bins))
             for d, pos in zip(datasets[1:], positive[1:])]
    annotate_jsd_pairs(ax, pairs, loc="upper left")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_raw_norm_log_compare_n(datasets: list[dict], output_path: Path, rng: np.random.Generator) -> None:
    staged = []
    for d in datasets:
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

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    stage_specs = [(axes[0], "raw", "Raw counts", True),
                   (axes[1], "norm", "Normalized", True),
                   (axes[2], "log", "log1p(normalized)", False)]
    for ax, key, title, log_x in stage_specs:
        vals_by_dataset = [s[key][s[key] > 0] for s in staged]
        all_vals = np.concatenate(vals_by_dataset)
        if log_x:
            bins = np.logspace(np.log10(all_vals.min()), np.log10(all_vals.max()), 60)
            ax.set_xscale("log")
        else:
            bins = np.linspace(all_vals.min(), all_vals.max(), 60)
        for s, vals in zip(staged, vals_by_dataset):
            ax.hist(vals, bins=bins, density=True, color=s["color"], alpha=0.5, label=s["display_label"])
        ax.set_title(title)
        ax.set_xlabel("Value (nonzero matrix entries)")
        pairs = [(s["display_label"], cd.compute_jsd(vals_by_dataset[0], vals, bins))
                 for s, vals in zip(staged[1:], vals_by_dataset[1:])]
        annotate_jsd_pairs(ax, pairs, loc="upper right")
    axes[0].set_ylabel("Density")
    axes[0].legend(frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def apply_qc(adata, label: str):
    """Drop off-tissue/empty rows, mirroring Figure 2's QC convention (see
    data/figure_2/smaller_sphere/data/README.md). None of the three cell-level
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


def hvg_match_panels(adatas: dict[str, "sc.AnnData"]) -> dict[str, "sc.AnnData"]:
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
        sc.pp.normalize_total(tmp, target_sum=TARGET_SUM)
        sc.pp.log1p(tmp)
        _, var = cd.gene_mean_var(tmp.X)
        top = np.argsort(var)[::-1][:target]
        out[label] = adata[:, np.sort(top)].copy()
        print(f"[hvg-match] {label}: {adata.n_vars} -> {out[label].n_vars} genes (target {target})")
    return out


def load_and_qc(spec: dict):
    print(f"[load] {spec['label']}: {spec['path']}")
    adata = sc.read_h5ad(spec["path"])
    print(f"[load] {spec['label']} shape before slice restriction: {adata.n_obs} x {adata.n_vars}")
    if "slice_id" not in adata.obs.columns:
        raise SystemExit(f"{spec['label']}: 'slice_id' not found in obs columns")
    keep = adata.obs["slice_id"] == SLICE_ID
    if not keep.any():
        raise SystemExit(f"{spec['label']}: slice_id={SLICE_ID} matches no observations "
                          f"(available: {sorted(adata.obs['slice_id'].unique().tolist())})")
    adata = adata[keep].copy()
    print(f"[load] {spec['label']} restricted to slice_id={SLICE_ID}: {adata.n_obs} x {adata.n_vars}")
    return apply_qc(adata, spec["label"])


def compute_stats(spec: dict, adata) -> dict:
    X = adata.X
    mean, var = cd.gene_mean_var(X)
    theta = cd.fit_common_dispersion(mean, var)
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
    output_dir = (FIGURE_4A_DIR / args.modality /
                  ("distribution_comparison_qc_hvg" if args.hvg_match else "distribution_comparison"))
    output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RANDOM_STATE)

    adatas = {spec["label"]: load_and_qc(spec) for spec in sources}
    if args.hvg_match:
        adatas = hvg_match_panels(adatas)
    datasets = [compute_stats(spec, adatas[spec["label"]]) for spec in sources]

    plot_mean_variance_compare_n(datasets, output_dir / "mean_variance_compare.png")
    plot_mean_dropout_compare_n(datasets, output_dir / "mean_dropout_compare.png")
    plot_sparsity_summary_compare_n(datasets, output_dir / "sparsity_summary_compare.png")
    plot_total_counts_compare_n(datasets, output_dir / "total_counts_compare.png")
    plot_genes_per_cell_compare_n(datasets, output_dir / "genes_per_cell_compare.png")
    plot_raw_norm_log_compare_n(datasets, output_dir / "raw_norm_log_compare.png", rng)
    plot_dataset_legend(datasets, output_dir / "dataset_legend.png")
    print(f"[save] Comparison plots written to {output_dir}")

    summary = {
        "modality": args.modality,
        "slice_id": SLICE_ID,
        "hvg_match": bool(args.hvg_match),
        "sources": {d["label"]: str(d["path"].relative_to(SIM_PAPER_DIR)) for d in datasets},
        "caveat": "our_method uses Figure 2's manuscript-tuned dispersion parameters (real batch "
                  "effects included); splatter_spider/splatter_sccube use the generic overview "
                  "qualitative-demo contract. n_cells and technical-noise model are NOT matched "
                  "across sources" + (" (n_genes IS matched via HVG selection)"
                  if args.hvg_match else "; n_genes is not matched either") + " -- see script docstring.",
        "n_obs": {d["label"]: int(d["adata"].n_obs) for d in datasets},
        "n_vars": {d["label"]: int(d["adata"].n_vars) for d in datasets},
        "theta_hat": {d["label"]: d["theta"] for d in datasets},
        "total_counts_median": {d["label"]: float(np.median(d["total_counts"])) for d in datasets},
        "total_counts_mean": {d["label"]: float(np.mean(d["total_counts"])) for d in datasets},
        "genes_per_cell_median": {d["label"]: float(np.median(d["n_genes"])) for d in datasets},
        "matrix_zero_frac": {d["label"]: d["sparsity"]["matrix_zero_frac"] for d in datasets},
        "empty_row_frac": {d["label"]: d["sparsity"]["empty_row_frac"] for d in datasets},
    }
    summary_path = output_dir / "comparison_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    print(f"[save] Numeric summary written to {summary_path}")


if __name__ == "__main__":
    main()
