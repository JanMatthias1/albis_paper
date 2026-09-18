#!/usr/bin/env python
"""
Manuscript panel: spatial-domain recovery ceiling (left, batch_sigma=0, the
strong-domain-mix BANKSY+Harmony testbed) vs. real cell-type recovery (right,
each modality's own canonical/tuned batch_sigma, plain PCA+Harmony -> Leiden,
weak/manuscript-baseline domain-mix) across cell / bin16um / spot.

The two panels intentionally use DIFFERENT datasets and batch conditions --
this is the point of the comparison, not an inconsistency: the left panel
asks "how much domain structure exists in principle, with no batch effect at
all to fight," the right panel asks "how well does the actual real-calibration
pipeline recover cell types once the real per-modality batch effect is
baked in." Domain recovery under the SAME canonical batch_sigma used on the
right is a separate, much lower number (see
cellbin_batch_sigma_slide/batch_sigma_slide_domain_ari_final.png for that
full cliff curve) -- deliberately not shown here.

bin8um is excluded: no strong-domain-mix data exists for it (out of scope
per the 2026-09-18 Figure 3B closure decision), so it has no left-panel value.

Left = data/figure_3/cellbin_batch_sigma_slide/<folder>/bs0/ari/ari_summary_<true_mod>.json
Right = data/figure_3/pca_harmony_single_cell/<dir>/ari_recovery_qc/ari_summary_<mod>.json

Supersedes code/clustering/misc/plot_domain_vs_celltype_recovery.py (archived
to misc/legacy/ 2026-09-18 -- it read domain data from banksy_batch_compare/,
which no longer exists, and paired prebatch+tuned bars rather than the single
batch_sigma=0 bar asked for here).

Usage:
    conda activate /dcs04/hicks/data/Jan/sim_project/sim_paper/env/albis-tutorial
    python sim_paper/code/clustering/ari_recovery_summary/plot_domain_vs_celltype.py
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
SIM_PAPER_DIR = SCRIPT_DIR.parents[2]
FIG3_DIR = SIM_PAPER_DIR / "data" / "figure_3"
SLIDE_ROOT = FIG3_DIR / "cellbin_batch_sigma_slide"
PLAIN_ROOT = FIG3_DIR / "pca_harmony_single_cell"

MODALITIES = ["cell", "bin16um", "spot"]
MODALITY_DISPLAY = {"cell": "Cell", "bin16um": "Bin (16µm)", "spot": "Spot"}
# folder under cellbin_batch_sigma_slide/, and the modality name baked into
# that folder's ari_summary_<...>.json filename (bin16um's pipeline was run
# with --modality bin, so its file is ari_summary_bin.json, not _bin16um.json)
DOMAIN_TRUE_MOD = {"cell": "cell", "bin16um": "bin", "spot": "spot"}
# the no-batch point's folder name isn't spelled consistently across
# modalities (cell/bin16um's generator formats 0 as "0", spot's as "0.0")
DOMAIN_BS0_DIR = {"cell": "bs0", "bin16um": "bs0", "spot": "bs0.0"}
# folder under pca_harmony_single_cell/, and that panel's summary-file modality
PLAIN_DIR = {"cell": "cell", "bin16um": "bin16um", "spot": "spot"}
PLAIN_TRUE_MOD = {"cell": "cell", "bin16um": "bin", "spot": "spot"}
# each modality's own tuned per-slice batch_sigma (see batch_sigma_slide_domain_ari_final.png)
BATCH_SIGMA = {"cell": 1.5, "bin16um": 0.7, "spot": 0.3}

COLOR_DOMAIN = "#2a78d6"    # dataviz palette slot 1, blue
COLOR_CELLTYPE = "#eb6834"  # dataviz palette slot 2, orange
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRIDLINE = "#e1e0d9"
BASELINE = "#c3c2b7"


def load_domain_bs0() -> dict[str, float]:
    out = {}
    for m in MODALITIES:
        true_mod = DOMAIN_TRUE_MOD[m]
        ari_path = SLIDE_ROOT / m / DOMAIN_BS0_DIR[m] / "ari" / f"ari_summary_{true_mod}.json"
        data = json.loads(ari_path.read_text())
        out[m] = next(r["ari"] for r in data if r["ground_truth"] == "domain_true")
    return out


def load_celltype_tuned() -> dict[str, float]:
    out = {}
    for m in MODALITIES:
        ari_path = PLAIN_ROOT / PLAIN_DIR[m] / "ari_recovery_qc" / f"ari_summary_{PLAIN_TRUE_MOD[m]}.json"
        data = json.loads(ari_path.read_text())
        out[m] = next(r["ari"] for r in data if r["ground_truth"] == "cell_type_true")
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


def plot_panel(ax, data: dict[str, float], color: str, title: str, value_suffix: dict[str, str] | None = None) -> None:
    x = np.arange(len(MODALITIES))
    values = [data[m] for m in MODALITIES]
    bars = ax.bar(x, values, 0.5, color=color, edgecolor="white", linewidth=0.6, zorder=2)
    if value_suffix:
        labels = [f"{v:.3f}\n{value_suffix[m]}" for v, m in zip(values, MODALITIES)]
    else:
        labels = [f"{v:.3f}" for v in values]
    ax.bar_label(bars, labels=labels, fontsize=9, padding=3, color=INK_SECONDARY, linespacing=1.6)
    ax.set_xticks(x)
    ax.set_xticklabels([MODALITY_DISPLAY[m] for m in MODALITIES], fontsize=11, color=INK_PRIMARY)
    ax.set_ylabel("Adjusted Rand Index", fontsize=10.5, color=INK_SECONDARY)
    ax.set_title(title, fontsize=13, color=INK_PRIMARY, pad=12)
    ax.set_ylim(-0.05, 1.08)
    style_axis(ax)


def main() -> None:
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial", "Helvetica"]

    domain_data = load_domain_bs0()
    celltype_data = load_celltype_tuned()
    domain_sigma_labels = {m: "σ=0" for m in MODALITIES}
    celltype_sigma_labels = {m: f"σ={BATCH_SIGMA[m]}" for m in MODALITIES}

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.2))
    plot_panel(axes[0], domain_data, COLOR_DOMAIN, "Spatial Domain Recovery", value_suffix=domain_sigma_labels)
    plot_panel(axes[1], celltype_data, COLOR_CELLTYPE, "Cell-type Recovery", value_suffix=celltype_sigma_labels)
    fig.tight_layout()

    out_dir = FIG3_DIR / "ari_recovery_summary"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "domain_vs_celltype_recovery.png"
    fig.savefig(out_path, dpi=300, facecolor="white")
    plt.close(fig)
    print(f"[save] {out_path}")

    print("\nDomain recovery (BANKSY+Harmony, batch_sigma=0):")
    for m in MODALITIES:
        print(f"  {MODALITY_DISPLAY[m]:12s} ari={domain_data[m]:.4f}")
    print("\nCell-type recovery (plain PCA+Harmony->Leiden, canonical/tuned batch):")
    for m in MODALITIES:
        print(f"  {MODALITY_DISPLAY[m]:12s} ari={celltype_data[m]:.4f}  (batch_sigma={BATCH_SIGMA[m]})")


if __name__ == "__main__":
    main()
