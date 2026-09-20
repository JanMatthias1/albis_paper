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

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from manuscript_style import MODALITY_LOOKUP, apply_style
apply_style()

SCRIPT_DIR = Path(__file__).resolve().parent
SIM_PAPER_DIR = SCRIPT_DIR.parents[2]
FIG3_DIR = SIM_PAPER_DIR / "data" / "figure_3"
SLIDE_ROOT = FIG3_DIR / "cellbin_batch_sigma_slide"
PLAIN_ROOT = FIG3_DIR / "pca_harmony_single_cell"

MODALITIES = ["cell", "bin16um", "spot"]
MODALITY_DISPLAY = {"cell": "Cell", "bin16um": "Bin (16 µm)", "spot": "Spot"}
# folder under cellbin_batch_sigma_slide/, and the modality name baked into
# that folder's ari_summary_<...>.json filename (bin16um's pipeline was run
# with --modality bin, so its file is ari_summary_bin.json, not _bin16um.json)
DOMAIN_TRUE_MOD = {"cell": "cell", "bin16um": "bin", "spot": "spot"}
# the no-batch point's folder name isn't spelled consistently across
# modalities (cell/bin16um's generator formats 0 as "0", spot's as "0.0")
DOMAIN_BS0_DIR = {"cell": "bs0", "bin16um": "bs0", "spot": "bs0.0"}
# bin16um/bs0/ was mid-regen ("poisson baseline" fix) as of 2026-09-18 --
# temporarily pointed at bs0_pre_poisson_baseline_20260918/ (the archived
# pre-regen ari=0.799) while job 35789683_0 (generate_strong_mix_bin16um.sh)
# was still queued. Reverted back to "bs0" here since this script is only
# meant to run once that job has landed fresh ari/ output (see
# replot_bin16um_bs0.sh, chained on that job via --dependency=afterok).
# folder under pca_harmony_single_cell/, and that panel's summary-file modality
PLAIN_DIR = {"cell": "cell", "bin16um": "bin16um", "spot": "spot"}
PLAIN_TRUE_MOD = {"cell": "cell", "bin16um": "bin", "spot": "spot"}
# each modality's own tuned per-slice batch_sigma (see batch_sigma_slide_domain_ari_final.png)
BATCH_SIGMA = {"cell": 1.5, "bin16um": 0.7, "spot": 0.3}

INK_PRIMARY = "#000000"
INK_SECONDARY = "#000000"
INK_MUTED = "#3d3d3d"
GRIDLINE = "#DDDDDD"
BASELINE = "#c3c2b7"

# Same type scale as code/count_distribution/count_distribution.py (Figure 2's
# plotting script) -- sized for legibility once shrunk into a multi-panel
# print figure, not just on-screen. Reused verbatim rather than re-derived so
# every figure in the paper reads at the same visual weight.
TITLE_SIZE = 17
LABEL_SIZE = 16
TICK_SIZE = 12
ANNOT_SIZE = 13


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
    ax.tick_params(axis="both", colors=INK_MUTED, labelsize=TICK_SIZE)
    ax.axhline(0, color=BASELINE, linewidth=1.0, zorder=1)


def plot_panel(ax, data: dict[str, float], title: str,
                value_suffix: dict[str, str] | None = None) -> None:
    x = np.arange(len(MODALITIES))
    values = [data[m] for m in MODALITIES]
    bars = ax.bar(x, values, 0.5, color=[MODALITY_LOOKUP[m] for m in MODALITIES], edgecolor="white", linewidth=0.6, zorder=2)
    if value_suffix:
        labels = [f"{v:.3f}\n{value_suffix[m]}" for v, m in zip(values, MODALITIES)]
    else:
        labels = [f"{v:.3f}" for v in values]
    ax.bar_label(bars, labels=labels, fontsize=ANNOT_SIZE, padding=3, color=INK_SECONDARY,
                 fontweight="bold", linespacing=1.6)
    ax.set_xticks(x)
    ax.set_xticklabels([MODALITY_DISPLAY[m] for m in MODALITIES], fontsize=LABEL_SIZE, color=INK_PRIMARY)
    ax.set_ylabel("Adjusted Rand Index", fontsize=LABEL_SIZE, color=INK_SECONDARY)
    ax.set_title(title, fontsize=TITLE_SIZE, fontweight="bold", color=INK_PRIMARY, pad=14)
    ax.set_ylim(-0.05, 1.08)
    style_axis(ax)


def main() -> None:
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial", "Helvetica"]
    plt.rcParams.update({
        "font.size": LABEL_SIZE,
        "axes.titlesize": TITLE_SIZE,
        "axes.titleweight": "bold",
        "axes.labelsize": LABEL_SIZE,
        "xtick.labelsize": TICK_SIZE,
        "ytick.labelsize": TICK_SIZE,
    })

    domain_data = load_domain_bs0()
    celltype_data = load_celltype_tuned()
    domain_sigma_labels = {m: "σ=0" for m in MODALITIES}
    celltype_sigma_labels = {m: f"σ={BATCH_SIGMA[m]}" for m in MODALITIES}

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 5.4))
    plot_panel(axes[0], domain_data, "Spatial Domain Recovery", value_suffix=domain_sigma_labels)
    plot_panel(axes[1], celltype_data, "Cell-type Recovery", value_suffix=celltype_sigma_labels)
    fig.tight_layout()

    out_dir = FIG3_DIR / "ari_recovery_summary"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "domain_vs_celltype_recovery.png"
    fig.savefig(out_path, dpi=300, facecolor="white", bbox_inches="tight")
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
