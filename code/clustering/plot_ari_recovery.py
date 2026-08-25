#!/usr/bin/env python
"""
Figure 3 summary: grouped bar charts of ARI recovery (Leiden vs. ground
truth, resolution-matched to true category count -- see
ari_vs_ground_truth.py) across modality (cell/bin/spot) and pipeline
(plain PCA+Harmony vs. BANKSY+PCA+Harmony), one chart for domain_true and
one for cell_type_true.

Auto-discovers every ari_summary_<modality>.json under
data/clustering_<packing-tag>/*/*ari_recovery*/ -- so pipelines that
haven't finished yet (e.g. BANKSY ARI jobs still pending) are simply
absent from the plot rather than needing this script to be told about
them by name. Rerun once more jobs land to pick up new bars.

Expected environment:
    conda activate /dcs04/hicks/data/Jan/sim_project/sim_app/env/sim-app-tutorial

Example:
    python sim_paper/code/clustering/plot_ari_recovery.py --packing-tag packing_pf0p04
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
SIM_PAPER_DIR = SCRIPT_DIR.parents[1]

MODALITY_ORDER = ["cell", "bin", "spot"]
GROUND_TRUTH_LABELS = {"domain_true": "Spatial domain recovery", "cell_type_true": "Cell type recovery"}


def pipeline_label(dirname: str) -> str:
    is_banksy = "banksy" in dirname
    is_qc = "_qc" in dirname
    base = "BANKSY + PCA + Harmony" if is_banksy else "PCA + Harmony"
    return f"{base} (QC)" if is_qc else base


def load_results(packing_tag: str) -> dict[str, dict[str, dict[str, float]]]:
    """Returns {ground_truth_col: {pipeline_label: {modality: ari}}}."""
    root = SIM_PAPER_DIR / "data" / f"clustering_{packing_tag}"
    out: dict[str, dict[str, dict[str, float]]] = {gt: {} for gt in GROUND_TRUTH_LABELS}

    for summary_path in sorted(root.glob("*/*ari_recovery*/ari_summary_*.json")):
        modality = summary_path.parent.parent.name
        pipeline = pipeline_label(summary_path.parent.name)
        results = json.loads(summary_path.read_text())
        for r in results:
            gt = r["ground_truth"]
            out.setdefault(gt, {}).setdefault(pipeline, {})[modality] = r["ari"]
            print(f"[load] {summary_path.relative_to(SIM_PAPER_DIR)}: {pipeline} / {modality} / {gt} = {r['ari']:.4f}")

    return out


def plot_one(ax, data_by_pipeline: dict[str, dict[str, float]], title: str) -> None:
    pipelines = sorted(data_by_pipeline.keys())
    n_pipelines = len(pipelines)
    x = np.arange(len(MODALITY_ORDER))
    width = 0.8 / max(n_pipelines, 1)

    for i, pipeline in enumerate(pipelines):
        values = [data_by_pipeline[pipeline].get(m, np.nan) for m in MODALITY_ORDER]
        offset = (i - (n_pipelines - 1) / 2) * width
        bars = ax.bar(x + offset, values, width, label=pipeline)
        ax.bar_label(bars, fmt="%.3f", fontsize=8, padding=2)

    ax.set_xticks(x)
    ax.set_xticklabels([m.capitalize() for m in MODALITY_ORDER])
    ax.set_ylabel("Adjusted Rand Index")
    ax.set_title(title)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylim(min(-0.05, ax.get_ylim()[0]), 1.05)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--packing-tag", default="packing_pf0p04")
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()

    if args.output_dir is None:
        args.output_dir = SIM_PAPER_DIR / "data" / f"clustering_{args.packing_tag}" / "ari_recovery_summary"
    args.output_dir.mkdir(parents=True, exist_ok=True)

    data = load_results(args.packing_tag)

    for gt, label in GROUND_TRUTH_LABELS.items():
        if not data.get(gt):
            print(f"[skip] no results found for {gt}")
            continue
        fig, ax = plt.subplots(figsize=(6, 4.5))
        plot_one(ax, data[gt], label)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2, frameon=False)
        fig.tight_layout()
        out_path = args.output_dir / f"ari_recovery_{gt}.png"
        fig.savefig(out_path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"[save] {out_path}")

    # Combined 1x2 figure for the paper.
    present_gts = [gt for gt in GROUND_TRUTH_LABELS if data.get(gt)]
    if len(present_gts) == 2:
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
        for ax, gt in zip(axes, present_gts):
            plot_one(ax, data[gt], GROUND_TRUTH_LABELS[gt])
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.02), ncol=2, frameon=False)
        fig.tight_layout()
        out_path = args.output_dir / "ari_recovery_combined.png"
        fig.savefig(out_path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"[save] {out_path}")


if __name__ == "__main__":
    main()
