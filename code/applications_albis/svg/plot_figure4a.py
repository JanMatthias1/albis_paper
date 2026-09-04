#!/usr/bin/env python
"""
Figure 4A -- SVG (spatially variable gene) recovery with scBSP.

Reads the per-modality `metrics_<dataset>.json` / `gene_pvalues_<dataset>.csv`
written by 3D_scbsp.py (data/figure_4/svg/scbsp/<dataset>/) and draws a two-row
summary:

  top row:    AUROC and AUPRC of -log10(p) separating true signal genes
              (is_noise==False, incl. marker_shared) from true noise genes
              (is_noise==True, 56 genes), one grouped bar per modality --
              same visual language as plot_ari_recovery.py / Figure 3.
  bottom row: per-modality distribution of -log10(p) split by ground-truth
              signal/noise, so a reader can see *how* separated the two
              populations are, not just the summary AUROC/AUPRC number.

A modality is silently skipped (with a printed note) if its output isn't on
disk yet -- same convention as plot_ari_recovery.py / plot_stagate_summary.py.

Expected environment:
    conda activate /dcs04/hicks/data/Jan/sim_project/sim_paper/env/albis-tutorial

Example:
    python sim_paper/code/applications_albis/svg/plot_figure4a.py
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

SCRIPT_DIR = Path(__file__).resolve().parent
SIM_PAPER_DIR = SCRIPT_DIR.parents[2]
DEFAULT_ROOT = SIM_PAPER_DIR / "data" / "figure_4" / "svg" / "scbsp"

MODALITY_ORDER = ["cell", "bin16um", "spot"]
MODALITY_DISPLAY = {"cell": "Cell", "bin16um": "Bin (16µm)", "spot": "Spot"}
SIGNAL_COLOR, NOISE_COLOR = "#2166ac", "#b2182b"


def load_results(root: Path) -> tuple[dict[str, dict], dict[str, pd.DataFrame]]:
    metrics, genes = {}, {}
    for modality in MODALITY_ORDER:
        mpath = root / modality / f"metrics_{modality}.json"
        gpath = root / modality / f"gene_pvalues_{modality}.csv"
        if not mpath.is_file() or not gpath.is_file():
            print(f"[skip] not found: {mpath.relative_to(SIM_PAPER_DIR) if mpath.exists() else mpath}")
            continue
        metrics[modality] = json.loads(mpath.read_text())
        genes[modality] = pd.read_csv(gpath)
        m = metrics[modality]
        print(f"[load] {modality}: AUROC={m['auroc']:.4f} AUPRC={m['auprc']:.4f} "
              f"n_signal={m['n_signal']} n_noise={m['n_noise']}")
    return metrics, genes


def plot_metrics_bar(ax, metrics: dict[str, dict]) -> None:
    present = [m for m in MODALITY_ORDER if m in metrics]
    x = np.arange(len(present))
    width = 0.35
    auroc = [metrics[m]["auroc"] for m in present]
    auprc = [metrics[m]["auprc"] for m in present]

    b1 = ax.bar(x - width / 2, auroc, width, label="AUROC", color="#4393c3")
    b2 = ax.bar(x + width / 2, auprc, width, label="AUPRC", color="#f4a582")
    ax.bar_label(b1, fmt="%.3f", fontsize=8, padding=2)
    ax.bar_label(b2, fmt="%.3f", fontsize=8, padding=2)

    ax.set_xticks(x)
    ax.set_xticklabels([MODALITY_DISPLAY[m] for m in present])
    ax.set_ylabel("Score")
    ax.set_title("SVG recovery (scBSP p-value vs. is_noise ground truth)")
    ax.axhline(0.5, color="black", linewidth=0.8, linestyle=":", label="chance (AUROC)")
    ax.set_ylim(0, 1.05)
    ax.legend(loc="lower right", fontsize=8, frameon=False)
    ax.grid(axis="y", ls=":", alpha=0.4)
    ax.set_axisbelow(True)


def plot_pvalue_distributions(axes, genes: dict[str, pd.DataFrame]) -> None:
    for ax, modality in zip(axes, MODALITY_ORDER):
        if modality not in genes:
            ax.axis("off")
            continue
        df = genes[modality]
        signal = df.loc[~df["is_noise"], "neg_log10_p"].to_numpy()
        noise = df.loc[df["is_noise"], "neg_log10_p"].to_numpy()

        parts = ax.violinplot([noise, signal], positions=[0, 1], showmedians=True, widths=0.8)
        for body, color in zip(parts["bodies"], [NOISE_COLOR, SIGNAL_COLOR]):
            body.set_facecolor(color)
            body.set_alpha(0.5)
        for key in ("cmedians", "cmins", "cmaxes", "cbars"):
            parts[key].set_color("0.3")

        rng = np.random.default_rng(0)
        for pos, vals, color in [(0, noise, NOISE_COLOR), (1, signal, SIGNAL_COLOR)]:
            jitter = rng.uniform(-0.08, 0.08, size=len(vals))
            ax.scatter(pos + jitter, vals, s=4, color=color, alpha=0.4, edgecolors="none")

        ax.set_xticks([0, 1])
        ax.set_xticklabels([f"noise\n(n={len(noise)})", f"signal\n(n={len(signal)})"])
        ax.set_title(MODALITY_DISPLAY[modality], fontsize=10)
        ax.grid(axis="y", ls=":", alpha=0.4)
        ax.set_axisbelow(True)
    axes[0].set_ylabel(r"$-\log_{10}(p)$  (scBSP)")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()
    output_dir = args.output_dir or args.root
    output_dir.mkdir(parents=True, exist_ok=True)

    metrics, genes = load_results(args.root)
    if not metrics:
        raise SystemExit(f"No scBSP results found under {args.root} -- run 3D_scbsp.py first.")

    fig = plt.figure(figsize=(11, 8))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.1, 1])
    ax_metrics = fig.add_subplot(gs[0, :])
    axes_dist = [fig.add_subplot(gs[1, i]) for i in range(3)]

    plot_metrics_bar(ax_metrics, metrics)
    plot_pvalue_distributions(axes_dist, genes)

    fig.suptitle("Figure 4A -- SVG identification (scBSP, full 3D z-stack)", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    out_path = output_dir / "figure_4a_svg_recovery.png"
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[save] {out_path}")

    summary_rows = [
        {"modality": m, **{k: v for k, v in metrics[m].items() if k not in ("input",)}}
        for m in MODALITY_ORDER if m in metrics
    ]
    csv_path = output_dir / "figure_4a_svg_summary.csv"
    pd.DataFrame(summary_rows).to_csv(csv_path, index=False)
    print(f"[save] {csv_path}")


if __name__ == "__main__":
    main()
