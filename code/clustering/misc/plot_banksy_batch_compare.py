#!/usr/bin/env python
"""
Spatial-domain recovery (ARI) vs. per-slice batch effect: BANKSY at a fixed
(lambda=0.5, k_geom=200 -- the lambda_kgeom_sweep winner on near-batch-free
data) run at prebatch vs. tuned batch_sigma, across cell / bin16um / spot.

Reads data/figure_3/banksy_batch_compare/<mod>/<batch>/ari/ari_summary_<mod>.json
(hard ARI) and .../scores/composition_recovery_<mod>_<batch>.json (slice_id
leakage ARI, annotated on the bars -- an ARI that's mostly slice-leakage is
not real domain recovery, see banksy_batch_compare.sh header).

Usage:
    conda activate /dcs04/hicks/data/Jan/sim_project/sim_paper/env/albis-tutorial
    python sim_paper/code/clustering/misc/plot_banksy_batch_compare.py
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
SIM_PAPER_DIR = SCRIPT_DIR.parents[2]
ROOT = SIM_PAPER_DIR / "data" / "figure_3" / "banksy_batch_compare"

ALL_MODALITIES = ["cell", "bin", "spot"]
# bin's official pair now lives under cellbin_batch_sigma_slide/bin/{bs0,bs0.7}
# (k_geom=100, see figure3_banksy_domain_sweep memory 2026-09-14 declutter) --
# skip any modality not present here rather than hardcode bin's absence, so
# this still works unmodified if a modality's dir structure changes again.
MODALITIES = [m for m in ALL_MODALITIES if (ROOT / m).is_dir()]
MODALITY_DISPLAY = {"cell": "Cell", "bin": "Bin (16µm)", "spot": "Spot"}
BATCHES = ["prebatch", "tuned"]
BATCH_DISPLAY = {"prebatch": "No batch effect", "tuned": "Tuned batch (canonical)"}
LEVELS = ["domain", "cell_type"]
LEVEL_DISPLAY = {"domain": "Spatial domain recovery (domain_true)", "cell_type": "Cell type recovery (cell_type_true)"}
LEVEL_TO_ARI_GT = {"domain": "domain_true", "cell_type": "cell_type_true"}


def load() -> dict[str, dict[str, dict[str, dict]]]:
    """{level: {modality: {batch: {"ari":..., "leak":...}}}}"""
    out = {lvl: {m: {} for m in MODALITIES} for lvl in LEVELS}
    for m in MODALITIES:
        for b in BATCHES:
            ari_path = ROOT / m / b / "ari" / f"ari_summary_{m}.json"
            score_path = ROOT / "scores" / f"composition_recovery_{m}_{b}.json"
            ari_by_gt = {r["ground_truth"]: r["ari"] for r in json.loads(ari_path.read_text())}
            scores = json.loads(score_path.read_text())["levels"]
            for lvl in LEVELS:
                gt = LEVEL_TO_ARI_GT[lvl]
                out[lvl][m][b] = {
                    "ari": ari_by_gt[gt],
                    "leak": scores[lvl]["slice_id_leakage_ari"],
                }
    return out


LEVEL_TITLE = {"domain": "Spatial domain recovery", "cell_type": "Cell type recovery"}


def plot_level(ax, data_by_mod: dict[str, dict[str, dict]], title: str) -> None:
    """Same visual convention as plot_ari_recovery.py's plot_one: grouped bars,
    one group per batch condition, bar_label values, legend below, axhline(0)."""
    x = np.arange(len(MODALITIES))
    n = len(BATCHES)
    width = 0.8 / n
    for i, b in enumerate(BATCHES):
        values = [data_by_mod[m][b]["ari"] for m in MODALITIES]
        leaks = [data_by_mod[m][b]["leak"] for m in MODALITIES]
        offset = (i - (n - 1) / 2) * width
        bars = ax.bar(x + offset, values, width, label=BATCH_DISPLAY[b])
        bar_labels = [f"{v:.3f}" + (" †" if leak > 0.3 else "") for v, leak in zip(values, leaks)]
        ax.bar_label(bars, labels=bar_labels, fontsize=8, padding=2)
    ax.set_xticks(x)
    ax.set_xticklabels([MODALITY_DISPLAY[m] for m in MODALITIES])
    ax.set_ylabel("Adjusted Rand Index")
    ax.set_title(title)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylim(min(-0.05, ax.get_ylim()[0]), 1.05)


def main() -> None:
    data = load()
    out_dir = ROOT / "plots"
    out_dir.mkdir(parents=True, exist_ok=True)

    for lvl in LEVELS:
        fig, ax = plt.subplots(figsize=(6, 4.5))
        plot_level(ax, data[lvl], LEVEL_TITLE[lvl])
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2, frameon=False)
        fig.text(0.5, -0.02, "† = slice_id leakage ARI > 0.3 (recovered \"clusters\" track slice, not signal)",
                  ha="center", fontsize=7.5, color="firebrick")
        fig.tight_layout()
        out_path = out_dir / f"ari_recovery_{lvl}_true_prebatch_vs_tuned.png"
        fig.savefig(out_path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"[save] {out_path}")

    # Combined 1x2 figure, same layout convention as ari_recovery_combined.png.
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for ax, lvl in zip(axes, LEVELS):
        plot_level(ax, data[lvl], LEVEL_TITLE[lvl])
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.02), ncol=2, frameon=False)
    fig.text(0.5, -0.04, "† = slice_id leakage ARI > 0.3 (recovered \"clusters\" track slice, not signal)",
              ha="center", fontsize=7.5, color="firebrick")
    fig.tight_layout()
    out_path = out_dir / "ari_recovery_combined_prebatch_vs_tuned.png"
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[save] {out_path}")

    print()
    print("Summary (ARI, slice-leakage ARI):")
    for lvl in LEVELS:
        print(f"  {LEVEL_DISPLAY[lvl]}")
        for m in MODALITIES:
            row = data[lvl][m]
            print(
                f"    {MODALITY_DISPLAY[m]:12s} "
                f"prebatch ARI={row['prebatch']['ari']:.4f} (leak {row['prebatch']['leak']:.3f})   "
                f"tuned ARI={row['tuned']['ari']:.4f} (leak {row['tuned']['leak']:.3f})"
            )


if __name__ == "__main__":
    main()
