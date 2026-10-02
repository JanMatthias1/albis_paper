#!/usr/bin/env python
"""
STAGATE domain recovery on the strong-mix data with NO batch effect
(batch_sigma = 0), 3D vs 2D spatial graph, per resolution.

Bars = mean over simulation seeds 101 / 202 / 2025 error bars = sample SD, labels = mean.
Same runs and colours as Figure 4B
(plot_stagate_fig4b.py), restricted to its sigma = 0 column.

Writes <STAGATE>/stagate_ari_strongmix_prebatch_seeds.{png,pdf,csv}.

Env: sim_paper/env/albis-tutorial.
"""
import os

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

from plot_stagate_fig4b import (C_2D, C_3D, GRID, MOD_DISPLAY, MODALITIES,
                                STAGATE_DIR, load, summarize)

METHODS = [("ari_domain_3d", C_3D, "STAGATE-3D"), ("ari_domain_2d", C_2D, "STAGATE-2D")]
OUT = os.path.join(STAGATE_DIR, "stagate_ari_strongmix_prebatch_seeds")


def main():
    df = load()
    strong = df[df.mix == "strong"].copy()
    summary = summarize(strong)  # validates exactly seeds 101/202/2025 per condition
    stats = summary[summary.batch_sigma == 0]
    stats.to_csv(OUT + ".csv", index=False)
    print(stats.to_string(index=False))

    fig, ax = plt.subplots(figsize=(8, 5.6))
    x = np.arange(len(MODALITIES))
    width = 0.36
    for i, (col, color, _) in enumerate(METHODS):
        offset = (i - 0.5) * width
        s = stats[stats.method == col].set_index("modality").loc[MODALITIES]
        ax.bar(x + offset, s["mean"], width, color=color, edgecolor="white", zorder=2)
        ax.errorbar(x + offset, s["mean"], yerr=s["sd"], fmt="none", ecolor="#2F2F2F",
                    elinewidth=1.3, capsize=5, zorder=4)
        for j, m in enumerate(MODALITIES):
            ax.annotate(f"{s.loc[m, 'mean']:.2f}", (x[j] + offset, s.loc[m, "mean"] + s.loc[m, "sd"]),
                        xytext=(0, 8), textcoords="offset points", ha="center",
                        fontsize=10, fontweight="bold", color="#0b0b0b")
    ax.set_xticks(x)
    ax.set_xticklabels([MOD_DISPLAY[m] for m in MODALITIES], fontsize=13)
    ax.set_ylabel("Domain ARI", fontsize=16)
    ax.set_ylim(0, 1)
    ax.tick_params(labelsize=12)
    ax.grid(axis="y", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("STAGATE domain recovery", fontsize=16, fontweight="bold", pad=12)
    ax.legend(handles=[Patch(facecolor=c, label=label) for _, c, label in METHODS],
              frameon=False, fontsize=11, loc="upper right")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{OUT}.{ext}", dpi=300, facecolor="white")
    plt.close(fig)
    print("wrote", OUT + ".png")


if __name__ == "__main__":
    main()
