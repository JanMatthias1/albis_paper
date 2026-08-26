#!/usr/bin/env python
"""
Figure 3B: grouped bar charts of ARI recovery (Leiden vs. ground truth,
resolution-matched to true category count -- see ari_vs_ground_truth.py)
across modality (cell/bin/spot) and pipeline (plain PCA+Harmony vs.
BANKSY+PCA+Harmony), one chart for domain_true and one for cell_type_true.

Reads the two known Figure 3 output locations directly -- NOT an
auto-discovery glob under one shared --packing-tag, because the cell-type and
domain panels don't live under one shared root/tag (see code/clustering/README.md,
"Output layout", flagged 2026-08-26): cell-type-panel results
(plain PCA+Harmony) always live at
data/figure_3/pca_harmony_single_cell/<modality>/ari_recovery_qc/, while
domain-panel results (BANKSY) live at data/clustering_<tag>/<modality>/banksy_ari_recovery/
with a DIFFERENT tag per modality (batch_sigma/domain_type_mix differ by
modality -- see --bin-domain-tag/--cell-domain-tag/--spot-domain-tag below).
A modality/pipeline combo simply doesn't appear in the chart if its
ari_summary_<modality>.json isn't on disk yet -- rerun once more jobs land.

Every summary file has both ground-truth columns computed regardless of which
pipeline it's from (ari_vs_ground_truth.py always does both), so this
necessarily also plots the "wrong" pairings (e.g. BANKSY's cell_type_true
collapse, plain's near-zero domain_true) -- keep them in rather than filter,
since they're exactly what motivated the pipeline split documented in
figure.md's 2026-08-25 "task/pipeline split" decision.

Expected environment:
    conda activate /dcs04/hicks/data/Jan/sim_project/sim_app/env/sim-app-tutorial

Example:
    python sim_paper/code/clustering/plot_ari_recovery.py
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

# Domain-panel (BANKSY) tag per modality -- see each modality's
# *_domain_panel.sh header for the current SIM_TAG; these get retuned as work
# continues, so override via --<modality>-domain-tag if a script's tag moves on
# without this default being updated.
DEFAULT_DOMAIN_TAGS = {
    "bin": "packing_pf0p04_bsigma05_strongdomainmix",
    "cell": "log_mu_-2.5_theta_0.25_strongdomainmix",
    "spot": "packing_pf0p04_bsigma015_strongdomainmix",
}


def pipeline_label(is_banksy: bool) -> str:
    base = "BANKSY + PCA + Harmony" if is_banksy else "PCA + Harmony"
    return f"{base} (QC)"


def load_results(domain_tags: dict[str, str]) -> dict[str, dict[str, dict[str, float]]]:
    """Returns {ground_truth_col: {pipeline_label: {modality: ari}}}."""
    out: dict[str, dict[str, dict[str, float]]] = {gt: {} for gt in GROUND_TRUTH_LABELS}

    def load_one(summary_path: Path, modality: str, is_banksy: bool) -> None:
        if not summary_path.is_file():
            print(f"[skip] not found: {summary_path.relative_to(SIM_PAPER_DIR)}")
            return
        pipeline = pipeline_label(is_banksy)
        results = json.loads(summary_path.read_text())
        for r in results:
            gt = r["ground_truth"]
            out.setdefault(gt, {}).setdefault(pipeline, {})[modality] = r["ari"]
            print(f"[load] {summary_path.relative_to(SIM_PAPER_DIR)}: {pipeline} / {modality} / {gt} = {r['ari']:.4f}")

    for modality in MODALITY_ORDER:
        load_one(
            SIM_PAPER_DIR / "data" / "figure_3" / "pca_harmony_single_cell" / modality
            / "ari_recovery_qc" / f"ari_summary_{modality}.json",
            modality, is_banksy=False,
        )
        load_one(
            SIM_PAPER_DIR / "data" / f"clustering_{domain_tags[modality]}" / modality
            / "banksy_ari_recovery" / f"ari_summary_{modality}.json",
            modality, is_banksy=True,
        )

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
    parser.add_argument("--bin-domain-tag", default=DEFAULT_DOMAIN_TAGS["bin"])
    parser.add_argument("--cell-domain-tag", default=DEFAULT_DOMAIN_TAGS["cell"])
    parser.add_argument("--spot-domain-tag", default=DEFAULT_DOMAIN_TAGS["spot"])
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()
    domain_tags = {"bin": args.bin_domain_tag, "cell": args.cell_domain_tag, "spot": args.spot_domain_tag}

    if args.output_dir is None:
        args.output_dir = SIM_PAPER_DIR / "data" / "figure_3" / "ari_recovery_summary"
    args.output_dir.mkdir(parents=True, exist_ok=True)

    data = load_results(domain_tags)

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
