#!/usr/bin/env python
"""
Same combined ARI bar chart as plot_banksy_batch_compare.py, but with bin
added back in. Bin's official prebatch/tuned pair doesn't live under
banksy_batch_compare/ (see that script's header) -- it's at
cellbin_batch_sigma_slide/bin/{bs0,bs0.7}/ (k_geom=100, the current config).
This script reads cell/spot from banksy_batch_compare/ and bin from
cellbin_batch_sigma_slide/, then plots all three together.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
SIM_PAPER_DIR = SCRIPT_DIR.parents[3]
FIG3_DIR = SIM_PAPER_DIR / "data" / "figure_3"
COMPARE_ROOT = FIG3_DIR / "banksy_batch_compare"
SLIDE_ROOT = FIG3_DIR / "cellbin_batch_sigma_slide"

MODALITIES = ["cell", "bin", "spot"]
MODALITY_DISPLAY = {"cell": "Cell", "bin": "Bin (16µm)", "spot": "Spot"}
BATCHES = ["prebatch", "tuned"]
BATCH_DISPLAY = {"prebatch": "No batch effect", "tuned": "Tuned batch (canonical)"}
LEVELS = ["domain", "cell_type"]
LEVEL_TITLE = {"domain": "Spatial domain recovery", "cell_type": "Cell type recovery"}
LEVEL_TO_ARI_GT = {"domain": "domain_true", "cell_type": "cell_type_true"}

# (ari_summary path, composition_recovery score path) per modality/batch
PATHS = {
    "cell": {
        "prebatch": (COMPARE_ROOT / "cell/prebatch/ari/ari_summary_cell.json", COMPARE_ROOT / "scores/composition_recovery_cell_prebatch.json"),
        "tuned": (COMPARE_ROOT / "cell/tuned/ari/ari_summary_cell.json", COMPARE_ROOT / "scores/composition_recovery_cell_tuned.json"),
    },
    "spot": {
        "prebatch": (COMPARE_ROOT / "spot/prebatch/ari/ari_summary_spot.json", COMPARE_ROOT / "scores/composition_recovery_spot_prebatch.json"),
        "tuned": (COMPARE_ROOT / "spot/tuned/ari/ari_summary_spot.json", COMPARE_ROOT / "scores/composition_recovery_spot_tuned.json"),
    },
    "bin": {
        "prebatch": (SLIDE_ROOT / "bin/bs0/ari/ari_summary_bin.json", SLIDE_ROOT / "scores/composition_recovery_bin_bs0.json"),
        "tuned": (SLIDE_ROOT / "bin/bs0.7/ari/ari_summary_bin.json", SLIDE_ROOT / "scores/composition_recovery_bin_bs0.7.json"),
    },
}


def load() -> dict[str, dict[str, dict[str, dict]]]:
    out = {lvl: {m: {} for m in MODALITIES} for lvl in LEVELS}
    for m in MODALITIES:
        for b in BATCHES:
            ari_path, score_path = PATHS[m][b]
            ari_by_gt = {r["ground_truth"]: r["ari"] for r in json.loads(ari_path.read_text())}
            scores = json.loads(score_path.read_text())["levels"]
            for lvl in LEVELS:
                gt = LEVEL_TO_ARI_GT[lvl]
                out[lvl][m][b] = {
                    "ari": ari_by_gt[gt],
                    "leak": scores[lvl]["slice_id_leakage_ari"],
                }
    return out


def plot_level(ax, data_by_mod: dict[str, dict[str, dict]], title: str) -> None:
    x = np.arange(len(MODALITIES))
    n = len(BATCHES)
    width = 0.8 / n
    for i, b in enumerate(BATCHES):
        values = [data_by_mod[m][b]["ari"] for m in MODALITIES]
        offset = (i - (n - 1) / 2) * width
        bars = ax.bar(x + offset, values, width, label=BATCH_DISPLAY[b])
        ax.bar_label(bars, fmt="%.3f", fontsize=8, padding=2)
    ax.set_xticks(x)
    ax.set_xticklabels([MODALITY_DISPLAY[m] for m in MODALITIES])
    ax.set_ylabel("Adjusted Rand Index")
    ax.set_title(title)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylim(min(-0.05, ax.get_ylim()[0]), 1.05)


def main() -> None:
    data = load()
    out_dir = COMPARE_ROOT / "plots"
    out_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    for ax, lvl in zip(axes, LEVELS):
        plot_level(ax, data[lvl], LEVEL_TITLE[lvl])
    # Reserve a dedicated bottom strip for the legend so tight_layout doesn't
    # let it collide with the axes (fig.legend sits outside the axes it
    # optimizes around).
    fig.tight_layout(rect=[0, 0.08, 1, 1])
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.05), ncol=2, frameon=False)
    out_path = out_dir / "ari_recovery_combined_prebatch_vs_tuned_3mod.png"
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[save] {out_path}")

    print()
    print("Summary (ARI, slice-leakage ARI):")
    for lvl in LEVELS:
        print(f"  {LEVEL_TITLE[lvl]}")
        for m in MODALITIES:
            row = data[lvl][m]
            print(
                f"    {MODALITY_DISPLAY[m]:12s} "
                f"prebatch ARI={row['prebatch']['ari']:.4f} (leak {row['prebatch']['leak']:.3f})   "
                f"tuned ARI={row['tuned']['ari']:.4f} (leak {row['tuned']['leak']:.3f})"
            )


if __name__ == "__main__":
    main()
