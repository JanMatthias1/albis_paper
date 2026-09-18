#!/usr/bin/env python
"""
Compact summary table of the Figure 3B batch_sigma-sweep configuration and
where each modality's domain recovery breaks down -- the companion piece to
batch_sigma_slide_domain_ari.png's curves (plot_batch_sigma_slide.py).
Values transcribed from that script's printed summary + the
figure3_banksy_domain_sweep memory (2026-09-13/14 sections); not
recomputed from json here since "where it breaks" is an interval read off
the curve, not a single stored number.

Usage:
    conda activate /dcs04/hicks/data/Jan/sim_project/sim_paper/env/albis-tutorial
    python plot_sweep_parameters_table.py
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent
SIM_PAPER_DIR = SCRIPT_DIR.parents[2]
OUT_PATH = SIM_PAPER_DIR / "data" / "figure_3" / "cellbin_batch_sigma_slide" / "batch_sigma_sweep_parameters_table.png"

COLUMNS = ["Modality", "BANKSY\n(λ, k_geom)", "Canonical\nbatch_σ", "Clean-recovery\nrange", "Break behavior", "Canonical vs. break"]
COL_WIDTHS = [0.12, 0.11, 0.09, 0.15, 0.32, 0.21]
WRAP_CHARS = [16, 14, 10, 18, 42, 26]

ROWS = [
    ["Cell", "λ=0.5\nk_geom=60", "1.5", "σ ≤ 0.1 (ARI 0.41–0.61)",
     "Sharp cliff at σ≈0.1–0.2 → ARI ~0.02, pure slice leakage (k_geom re-tuned "
     "2026-09-18; only 9 coarse points, cliff not finely bracketed)", "~15x past clean recovery"],
    ["Bin 16µm", "λ=0.5\nk_geom=100", "0.7", "σ ≤ 0.25 (ARI 0.63–0.80)",
     "Sharp cliff at σ≈0.25–0.3 → ARI 0.015, pure slice leakage", "Just past the cliff"],
    ["Spot\n(strongmix)", "λ=0.1\nk_geom=8", "0.3", "σ ≤ 0.12 (ARI 0.42–0.49)",
     "Non-monotonic notch: collapses to ARI≈0 (leak≈0, not leakage) for σ≈0.13–0.28, then sharply recovers to ARI 0.26–0.38 by σ≈0.29–0.3",
     "Sits inside the recovery band, not past a simple cliff"],
    ["Spot\n(real domain-mix)", "λ=0.1\nk_geom=8", "0.3", "none",
     "ARI ≈ 0 at every σ tested, including σ=0 (no batch effect at all)",
     "N/A — no domain signal to lose"],
]

ROW_COLORS = ["#eef3fb", "#eef3fb", "#fdf1e3", "#f6f6f6"]


def wrapped(rows: list[list[str]]) -> list[list[str]]:
    out = []
    for row in rows:
        out.append([
            "\n".join(line for part in cell.split("\n") for line in (textwrap.wrap(part, w) or [""]))
            for cell, w in zip(row, WRAP_CHARS)
        ])
    return out


def main() -> None:
    fig, ax = plt.subplots(figsize=(16, 5.2))
    ax.axis("off")

    table = ax.table(
        cellText=wrapped(ROWS),
        colLabels=COLUMNS,
        colWidths=COL_WIDTHS,
        cellLoc="left",
        colLoc="left",
        loc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)

    cells = table.get_celld()
    n_rows = max(row for row, _ in cells) + 1
    row_lines = {
        row: max(cells[(row, col)].get_text().get_text().count("\n") + 1 for col in range(len(COLUMNS)))
        for row in range(n_rows)
    }
    row_height = {row: 0.075 + 0.048 * (n - 1) for row, n in row_lines.items()}

    for (row, col), cell in cells.items():
        cell.set_edgecolor("#cccccc")
        cell.PAD = 0.015
        cell.get_text().set_verticalalignment("center")
        cell.set_height(row_height[row])
        if row == 0:
            cell.set_text_props(weight="bold", color="white", verticalalignment="center")
            cell.set_facecolor("#4a5a70")
        else:
            cell.set_facecolor(ROW_COLORS[row - 1])

    ax.set_title(
        "Figure 3B batch_sigma sweep: config and where domain recovery breaks\n"
        "(BANKSY + Harmony + Leiden, on the strong-domain-mix testbed unless noted)",
        fontsize=11, pad=14,
    )
    fig.tight_layout()
    fig.savefig(OUT_PATH, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[save] {OUT_PATH}")


if __name__ == "__main__":
    main()
