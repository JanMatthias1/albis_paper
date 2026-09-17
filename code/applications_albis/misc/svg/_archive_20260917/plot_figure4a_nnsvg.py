#!/usr/bin/env python
"""
Figure 4A -- SVG (spatially variable gene) recovery with nnSVG.

nnSVG counterpart to plot_figure4a.py (which draws the scBSP version). Reads
the per-run ``metrics_<run_tag>.json`` / ``gene_results_<run_tag>.csv`` written
by 3D_nnsvg.py under data/figure_4/svg/nnsvg/<run_tag>/ and draws:

  top row     AUROC of nnSVG's -log10(p) separating true signal genes
              (is_noise==False, incl. marker_shared -- 500 genes) from true
              noise genes (is_noise==True, 56 genes), one grouped bar per
              modality: canonical weak domain mix vs. the strengthened
              domain mix vs. the strengthened mix with the per-slice batch
              effect removed (oracle). Chance = 0.5.
  bottom row  per-modality ROC curve for the STRONG-MIX run (the curve whose
              area is the matching bar above). ROC rather than raw -log10(p)
              distributions because nnSVG p-values underflow to ~1e-320 for a
              chunk of genes at bin/spot, which flattens any distribution plot.

The story this panel is meant to carry (see figure.md / project memory):
SVG detectability tracks the strength of the simulated spatial organization
-- near chance on the canonical weak-mix data for every method tried
(nnSVG/scBSP/SPARK-X), recovered by nnSVG once the domain/cell-type coupling
is strengthened. It is a controllability result, not a method failure.

Expected environment:
    conda activate /dcs04/hicks/data/Jan/sim_project/sim_paper/env/albis-tutorial

Example:
    python sim_paper/code/applications_albis/svg/plot_figure4a_nnsvg.py
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from sklearn.metrics import roc_curve  # noqa: E402

SCRIPT_DIR = Path(__file__).resolve().parent
SIM_PAPER_DIR = SCRIPT_DIR.parents[2]
DEFAULT_ROOT = SIM_PAPER_DIR / "data" / "figure_4" / "svg" / "nnsvg"

MODALITY_ORDER = ["cell", "bin16um", "spot"]
MODALITY_DISPLAY = {"cell": "Cell", "bin16um": "Bin (16µm)", "spot": "Spot"}

# (suffix on the run_tag, legend label, bar colour)
CONDITIONS = [
    ("", "canonical (weak domain mix)", "#bdbdbd"),
    ("_strongmix", "strong domain mix", "#4393c3"),
    ("_strongmix_prebatch", "strong mix, no batch (oracle)", "#2166ac"),
]
SIGNAL_COLOR, NOISE_COLOR = "#2166ac", "#b2182b"
EPS = float(np.finfo(np.float64).tiny)


def _neg_log10_p(df: pd.DataFrame) -> np.ndarray:
    return -np.log10(df["pval"].to_numpy(dtype=float).clip(min=EPS))


def load_all(root: Path) -> tuple[dict, dict]:
    """metrics[(modality, cond_suffix)] -> dict ; genes[(modality, cond_suffix)] -> DataFrame"""
    metrics, genes = {}, {}
    for modality in MODALITY_ORDER:
        for suffix, _label, _colour in CONDITIONS:
            run_tag = f"{modality}{suffix}"
            mpath = root / run_tag / f"metrics_{run_tag}.json"
            gpath = root / run_tag / f"gene_results_{run_tag}.csv"
            if not mpath.is_file() or not gpath.is_file():
                print(f"[skip] {run_tag}: no results on disk")
                continue
            metrics[(modality, suffix)] = json.loads(mpath.read_text())
            genes[(modality, suffix)] = pd.read_csv(gpath)
            m = metrics[(modality, suffix)]
            print(f"[load] {run_tag:28s} AUROC={m['auroc']:.4f} "
                  f"n_obs={m['n_obs']} slice={m.get('slice_id')}")
    return metrics, genes


def plot_auroc_bars(ax, metrics: dict) -> None:
    present_mod = [m for m in MODALITY_ORDER if any((m, s) in metrics for s, _l, _c in CONDITIONS)]
    x = np.arange(len(present_mod))
    n = len(CONDITIONS)
    width = 0.8 / n
    for j, (suffix, label, colour) in enumerate(CONDITIONS):
        vals = [metrics.get((m, suffix), {}).get("auroc", np.nan) for m in present_mod]
        off = (j - (n - 1) / 2) * width
        bars = ax.bar(x + off, vals, width, label=label, color=colour)
        ax.bar_label(bars, fmt="%.2f", fontsize=8, padding=2)
    ax.axhline(0.5, color="black", lw=0.9, ls=":", zorder=0)
    ax.text(len(present_mod) - 0.5, 0.505, "chance", fontsize=7, va="bottom", ha="right", color="0.3")
    ax.set_xticks(x)
    ax.set_xticklabels([MODALITY_DISPLAY[m] for m in present_mod])
    ax.set_ylabel("AUROC  (nnSVG $-\\log_{10}p$ vs. is_noise)")
    ax.set_ylim(0, 1.02)
    ax.set_title("SVG recovery by nnSVG: detectability tracks spatial-organization strength")
    ax.legend(loc="upper left", fontsize=8, frameon=False)
    ax.grid(axis="y", ls=":", alpha=0.4)
    ax.set_axisbelow(True)


def plot_roc_curves(axes, genes: dict, metrics: dict, feature_suffix: str) -> None:
    """One ROC curve per modality for the featured run -- the curve whose area
    is the matching bar in the top panel. Avoids the -log10(p) underflow-spike
    scaling problem the raw distributions have at bin/spot."""
    for ax, modality in zip(axes, MODALITY_ORDER):
        ax.plot([0, 1], [0, 1], ls=":", color="0.5", lw=1)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_aspect("equal")
        ax.set_xlabel("false positive rate")
        ax.grid(ls=":", alpha=0.4)
        ax.set_axisbelow(True)
        key = (modality, feature_suffix)
        if key not in genes:
            ax.set_title(f"{MODALITY_DISPLAY[modality]}  (no run)", fontsize=10)
            continue
        df = genes[key]
        y_true = (~df["is_noise"].astype(bool)).astype(int).to_numpy()
        fpr, tpr, _ = roc_curve(y_true, _neg_log10_p(df))
        auroc = metrics.get(key, {}).get("auroc")
        ax.plot(fpr, tpr, color="#2166ac", lw=2)
        ax.fill_between(fpr, tpr, alpha=0.12, color="#2166ac")
        ax.set_title(f"{MODALITY_DISPLAY[modality]}   AUROC {auroc:.2f}" if auroc is not None
                     else MODALITY_DISPLAY[modality], fontsize=10)
    axes[0].set_ylabel("true positive rate\n(strong domain mix)")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--feature-condition", default="_strongmix",
                        help="which run's -log10(p) distributions to show in the bottom row "
                             "(default: _strongmix)")
    args = parser.parse_args()
    out_dir = args.output_dir or args.root
    out_dir.mkdir(parents=True, exist_ok=True)

    metrics, genes = load_all(args.root)
    if not metrics:
        raise SystemExit(f"No nnSVG results under {args.root} -- run 3D_nnsvg.py first.")

    fig = plt.figure(figsize=(11, 8))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.15, 1])
    ax_bars = fig.add_subplot(gs[0, :])
    axes_dist = [fig.add_subplot(gs[1, i]) for i in range(3)]

    plot_auroc_bars(ax_bars, metrics)
    plot_roc_curves(axes_dist, genes, metrics, args.feature_condition)

    fig.suptitle("Figure 4A — spatially variable gene identification (nnSVG, single slice)", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out_path = out_dir / "figure_4a_nnsvg_recovery.png"
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[save] {out_path}")

    rows = []
    for modality in MODALITY_ORDER:
        for suffix, label, _c in CONDITIONS:
            m = metrics.get((modality, suffix))
            if m is None:
                continue
            rows.append({
                "modality": modality, "condition": label,
                "run_tag": m["run_tag"], "tag": m["tag"],
                "n_obs": m["n_obs"], "slice_id": m.get("slice_id"),
                "use_pre_batch": m.get("use_pre_batch"),
                "auroc": m["auroc"], "auprc": m["auprc"],
            })
    csv_path = out_dir / "figure_4a_nnsvg_summary.csv"
    pd.DataFrame(rows).to_csv(csv_path, index=False)
    print(f"[save] {csv_path}")


if __name__ == "__main__":
    main()
