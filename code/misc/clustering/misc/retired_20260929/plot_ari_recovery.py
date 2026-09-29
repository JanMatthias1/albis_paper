#!/usr/bin/env python
"""
Figure 3B: grouped bar chart of cell-type-recovery ARI (Leiden vs.
cell_type_true, resolution-matched to true category count -- see
02_leiden_resolution_sweep.py) across modality (cell/bin/bin16um/spot), plain
PCA+Harmony only.

Reads data/figure_3/weak_domain_mix/pca_harmony/<modality>/bs<sigma>/ari_recovery_qc/
ari_summary_<modality>.json. A modality simply doesn't appear in the chart if
that file isn't on disk yet -- rerun once more jobs land.

2026-09-16: dropped the domain_true / BANKSY half (and the --bin/cell/spot-
domain-tag flags, --ground-truths, and the combined 1x2 figure that paired
them) -- that side only ever read data/clustering_<tag>/<modality>/
banksy_ari_recovery/, which only legacy/*_domain_panel.sh populates (moved to
legacy/ as superseded, last built 2026-08-25). The current domain-vs-celltype
comparison is misc/plot_domain_vs_celltype_recovery.py, reading fresher
domain data from banksy_batch_compare/ + cellbin_batch_sigma_slide/bin/.

Expected environment:
    conda activate /dcs04/hicks/data/Jan/sim_project/albis/env/albis-tutorial

Example:
    python sim_paper/code/clustering/ari_recovery_summary/plot_ari_recovery.py
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from manuscript_style import MODALITY_LOOKUP, apply_style
apply_style()

SCRIPT_DIR = Path(__file__).resolve().parent
SIM_PAPER_DIR = SCRIPT_DIR.parents[2]

MODALITY_ORDER = ["cell", "bin", "bin16um", "spot"]
MODALITY_DISPLAY = {"cell": "Cell", "bin": "Bin (8 µm)", "bin16um": "Bin (16 µm)", "spot": "Spot"}
# label -> (dir under data/figure_3/weak_domain_mix/pca_harmony/, modality name in the
# summary filename). bin16um's cell-type panel lives in its own bin16um/ dir but
# 02_leiden_resolution_sweep.py was run with --modality bin, so the file is ari_summary_bin.json.
CELLTYPE_PANEL = {
    "cell": ("cell/bs1.5", "cell"),
    "bin": ("bin8um/bs0.7", "bin"),
    "bin16um": ("bin16um/bs0.7", "bin"),
    "spot": ("spot/bs0.3", "spot"),
}
GROUND_TRUTH = "cell_type_true"
GROUND_TRUTH_LABEL = "Cell type recovery"
PIPELINE_LABEL = "PCA + Harmony (QC)"


def load_results() -> dict[str, dict[str, float]]:
    """Returns {pipeline_label: {modality: ari}} for cell_type_true."""
    out: dict[str, dict[str, float]] = {PIPELINE_LABEL: {}}

    for modality in MODALITY_ORDER:
        ct_dir, ct_file_mod = CELLTYPE_PANEL[modality]
        summary_path = (
            SIM_PAPER_DIR / "data" / "figure_3" / "weak_domain_mix" / "pca_harmony" / ct_dir
            / "ari_recovery_qc" / f"ari_summary_{ct_file_mod}.json"
        )
        if not summary_path.is_file():
            print(f"[skip] not found: {summary_path.relative_to(SIM_PAPER_DIR)}")
            continue
        results = json.loads(summary_path.read_text())
        for r in results:
            if r["ground_truth"] != GROUND_TRUTH:
                continue
            out[PIPELINE_LABEL][modality] = r["ari"]
            print(f"[load] {summary_path.relative_to(SIM_PAPER_DIR)}: {modality} / {GROUND_TRUTH} = {r['ari']:.4f}")

    return out


def plot_one(ax, data_by_pipeline: dict[str, dict[str, float]], title: str) -> None:
    pipelines = sorted(data_by_pipeline.keys())
    n_pipelines = len(pipelines)
    x = np.arange(len(MODALITY_ORDER))
    width = 0.8 / max(n_pipelines, 1)

    for i, pipeline in enumerate(pipelines):
        values = [data_by_pipeline[pipeline].get(m, np.nan) for m in MODALITY_ORDER]
        offset = (i - (n_pipelines - 1) / 2) * width
        bars = ax.bar(x + offset, values, width, label=pipeline, color=[MODALITY_LOOKUP[m] for m in MODALITY_ORDER])
        ax.bar_label(bars, fmt="%.3f", fontsize=12, padding=2)

    ax.set_xticks(x)
    ax.set_xticklabels([MODALITY_DISPLAY[m] for m in MODALITY_ORDER])
    ax.set_ylabel("Adjusted Rand Index")
    ax.set_title(title)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylim(min(-0.05, ax.get_ylim()[0]), 1.05)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()

    if args.output_dir is None:
        args.output_dir = SIM_PAPER_DIR / "data" / "figure_3" / "ari_recovery_summary"
    args.output_dir.mkdir(parents=True, exist_ok=True)

    data = load_results()
    if not data[PIPELINE_LABEL]:
        print(f"[skip] no results found for {GROUND_TRUTH}")
        return

    fig, ax = plt.subplots(figsize=(6, 4.5))
    plot_one(ax, data, GROUND_TRUTH_LABEL)
    ax.text(0.5, -0.18, PIPELINE_LABEL, transform=ax.transAxes, ha="center", fontsize=12)
    fig.tight_layout()
    out_path = args.output_dir / f"ari_recovery_{GROUND_TRUTH}.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[save] {out_path}")


if __name__ == "__main__":
    main()
