#!/usr/bin/env python
"""
Domain recovery (domain_true ARI) on the REAL (non-strengthened,
manuscript-baseline domain_type_mix) data, prebatch (batch_sigma=0) only --
across cell / bin16um / spot. Isolates whether there is any BANKSY-findable
domain signal at all before batch effect enters the picture.

Sources (each modality's own tuned lambda/k_geom, prebatch):
  cell:    data/figure_3/pca_harmony_domain_figure2_data/cell/prebatch/ari/ari_summary_cell.json
  bin16um: data/figure_3/pca_harmony_domain_figure2_data/bin16um/prebatch/ari/ari_summary_bin.json
  spot:    data/figure_3/cellbin_batch_sigma_slide/spot_weakmix/bs0.0/ari/ari_summary_spot.json
           (spot's pca_harmony_domain_check/spot/prebatch point, once that
           task lands, is a --use-pre-batch cross-check of this same value --
           see figure3_banksy_domain_sweep memory, 2026-09-14.)

Usage:
    conda activate /dcs04/hicks/data/Jan/sim_project/sim_paper/env/albis-tutorial
    python plot_domain_recovery_realmix_prebatch.py
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent
SIM_PAPER_DIR = SCRIPT_DIR.parents[2]

SOURCES = {
    "Cell\n(λ=0.5, k_geom=200)": SIM_PAPER_DIR / "data/figure_3/pca_harmony_domain_figure2_data/cell/prebatch/ari/ari_summary_cell.json",
    "Bin 16µm\n(λ=0.5, k_geom=100)": SIM_PAPER_DIR / "data/figure_3/pca_harmony_domain_figure2_data/bin16um/prebatch/ari/ari_summary_bin.json",
    "Spot\n(λ=0.1, k_geom=8)": SIM_PAPER_DIR / "data/figure_3/cellbin_batch_sigma_slide/spot_weakmix/bs0.0/ari/ari_summary_spot.json",
}


def domain_ari(path: Path) -> float:
    data = json.loads(path.read_text())
    return next(r["ari"] for r in data if r["ground_truth"] == "domain_true")


def main() -> None:
    labels = list(SOURCES.keys())
    values = [domain_ari(p) for p in SOURCES.values()]

    fig, ax = plt.subplots(figsize=(6, 4.5))
    bars = ax.bar(labels, values, color=["C0", "C1", "C4"], width=0.55)
    ax.bar_label(bars, labels=[f"{v:.4f}" for v in values], fontsize=10, padding=3)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("Domain ARI (BANKSY + Harmony)")
    ax.set_title("Domain recovery on REAL domain-mix data, prebatch\n(no batch effect, no strong-domain-mix strengthening)")
    ax.set_ylim(-0.05, 1.0)
    fig.tight_layout()

    out_path = SIM_PAPER_DIR / "data/figure_3/pca_harmony_domain_figure2_data/domain_recovery_realmix_prebatch.png"
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[save] {out_path}")

    print("\nSummary:")
    for label, v in zip(labels, values):
        print(f"  {label.splitlines()[0]:10s} domain_true ARI = {v:.4f}")


if __name__ == "__main__":
    main()
