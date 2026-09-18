#!/usr/bin/env python
"""
Manuscript panel: BANKSY spatial-domain recovery (left) vs. plain
Harmony -> Leiden cell-type recovery, no BANKSY (right), across cell /
bin16um / spot.

Left = banksy_batch_compare/<mod>/{prebatch,tuned} (+ bin from
cellbin_batch_sigma_slide/bin/{bs0,bs0.7}), same data as
plot_banksy_batch_compare_3mod.py's domain panel.

Right = pca_harmony_single_cell/<mod>/ari_recovery_qc/ari_summary_<mod>.json
(cell_type_true), the plain PCA -> Harmony -> Leiden pipeline used for
Figure 2's cell-type panels -- no BANKSY spatial smoothing. Only one
condition exists here (the canonical/tuned batch QC data, matching Figure 2;
no prebatch companion run), so this panel is single-bar-per-modality, not
grouped -- and no slice_id-leakage check was computed for this pipeline, so
no dagger annotation on this side.

Usage:
    conda activate /dcs04/hicks/data/Jan/sim_project/sim_paper/env/albis-tutorial
    python sim_paper/code/clustering/misc/plot_domain_vs_celltype_recovery.py
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
SIM_PAPER_DIR = SCRIPT_DIR.parents[2]
FIG3_DIR = SIM_PAPER_DIR / "data" / "figure_3"
COMPARE_ROOT = FIG3_DIR / "banksy_batch_compare"
SLIDE_ROOT = FIG3_DIR / "cellbin_batch_sigma_slide"
PLAIN_ROOT = FIG3_DIR / "pca_harmony_single_cell"

MODALITIES = ["cell", "bin", "spot"]
MODALITY_DISPLAY = {"cell": "Cell", "bin": "Bin (16µm)", "spot": "Spot"}
BATCHES = ["prebatch", "tuned"]
BATCH_DISPLAY = {"prebatch": "No batch effect", "tuned": "Tuned batch"}
# each modality's own tuned per-slice batch effect magnitude -- "canonical" means
# this batch_sigma, not a shared value across modalities (see batch_sigma_slide
# cliff analysis, cellbin_batch_sigma_slide/batch_sigma_slide_domain_ari_final.png)
BATCH_SIGMA = {"cell": 1.5, "bin": 0.7, "spot": 0.3}

# validated categorical palette (dataviz skill, references/palette.md), slots 1-2
COLOR_NO_BATCH = "#2a78d6"   # slot 1, blue
COLOR_TUNED = "#eb6834"      # slot 2, orange
COLOR_PLAIN = "#1baf7a"      # slot 3, aqua -- distinct method, own color
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRIDLINE = "#e1e0d9"
BASELINE = "#c3c2b7"

DOMAIN_PATHS = {
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
PLAIN_PATHS = {
    "cell": PLAIN_ROOT / "cell/ari_recovery_qc/ari_summary_cell.json",
    "bin": PLAIN_ROOT / "bin16um/ari_recovery_qc/ari_summary_bin.json",
    "spot": PLAIN_ROOT / "spot/ari_recovery_qc/ari_summary_spot.json",
}


def load_domain() -> dict[str, dict[str, float]]:
    out: dict[str, dict[str, float]] = {m: {} for m in MODALITIES}
    for m in MODALITIES:
        for b in BATCHES:
            ari_path, _ = DOMAIN_PATHS[m][b]
            ari_by_gt = {r["ground_truth"]: r["ari"] for r in json.loads(ari_path.read_text())}
            out[m][b] = ari_by_gt["domain_true"]
    return out


def load_plain_celltype() -> dict[str, float]:
    out = {}
    for m in MODALITIES:
        d = json.loads(PLAIN_PATHS[m].read_text())
        out[m] = next(r["ari"] for r in d if r["ground_truth"] == "cell_type_true")
    return out


def style_axis(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_color(BASELINE)
    ax.yaxis.grid(True, color=GRIDLINE, linewidth=0.9, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(axis="both", colors=INK_MUTED, labelsize=10)
    ax.axhline(0, color=BASELINE, linewidth=1.0, zorder=1)


def plot_domain_panel(ax, data: dict[str, dict[str, float]]) -> None:
    x = np.arange(len(MODALITIES))
    width = 0.34
    for i, b in enumerate(BATCHES):
        color = COLOR_NO_BATCH if b == "prebatch" else COLOR_TUNED
        values = [data[m][b] for m in MODALITIES]
        offset = (i - 0.5) * width
        bars = ax.bar(x + offset, values, width, color=color, label=BATCH_DISPLAY[b],
                       edgecolor="white", linewidth=0.6, zorder=2)
        if b == "tuned":
            labels = [f"{v:.2f}\nσ={BATCH_SIGMA[m]}" for v, m in zip(values, MODALITIES)]
        else:
            labels = [f"{v:.2f}" for v in values]
        ax.bar_label(bars, labels=labels, fontsize=9, padding=3, color=INK_SECONDARY, linespacing=1.6)
    ax.set_xticks(x)
    ax.set_xticklabels([MODALITY_DISPLAY[m] for m in MODALITIES], fontsize=11, color=INK_PRIMARY)
    ax.set_ylabel("Adjusted Rand Index", fontsize=10.5, color=INK_SECONDARY)
    ax.set_title("Spatial Domain Recovery", fontsize=13, color=INK_PRIMARY, pad=12)
    ax.set_ylim(-0.05, 1.08)
    style_axis(ax)


def plot_plain_panel(ax, data: dict[str, float]) -> None:
    x = np.arange(len(MODALITIES))
    values = [data[m] for m in MODALITIES]
    bars = ax.bar(x, values, 0.5, color=COLOR_PLAIN, edgecolor="white", linewidth=0.6, zorder=2)
    labels = [f"{v:.2f}\nσ={BATCH_SIGMA[m]}" for v, m in zip(values, MODALITIES)]
    ax.bar_label(bars, labels=labels, fontsize=9, padding=3, color=INK_SECONDARY, linespacing=1.6)
    ax.set_xticks(x)
    ax.set_xticklabels([MODALITY_DISPLAY[m] for m in MODALITIES], fontsize=11, color=INK_PRIMARY)
    ax.set_ylabel("Adjusted Rand Index", fontsize=10.5, color=INK_SECONDARY)
    ax.set_title("Cell-type Recovery", fontsize=13, color=INK_PRIMARY, pad=12)
    ax.set_ylim(-0.05, 1.08)
    style_axis(ax)


def main() -> None:
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial", "Helvetica"]

    domain_data = load_domain()
    plain_data = load_plain_celltype()

    fig, axes = plt.subplots(1, 2, figsize=(11, 3.4))
    plot_domain_panel(axes[0], domain_data)
    plot_plain_panel(axes[1], plain_data)

    fig.tight_layout(rect=[0, 0.09, 1, 1])
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.28, 0.02), ncol=2,
               frameon=False, fontsize=10, labelcolor=INK_SECONDARY)

    out_dir = FIG3_DIR / "banksy_batch_compare" / "plots"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "domain_banksy_vs_celltype_plain.png"
    fig.savefig(out_path, dpi=300, facecolor="white")
    plt.close(fig)
    print(f"[save] {out_path}")

    print("\nDomain recovery (BANKSY):")
    for m in MODALITIES:
        row = domain_data[m]
        print(f"  {MODALITY_DISPLAY[m]:12s} prebatch={row['prebatch']:.4f}  tuned={row['tuned']:.4f} (σ={BATCH_SIGMA[m]})")
    print("\nCell-type recovery (plain Harmony->Leiden, canonical batch):")
    for m in MODALITIES:
        print(f"  {MODALITY_DISPLAY[m]:12s} ari={plain_data[m]:.4f}")


if __name__ == "__main__":
    main()
