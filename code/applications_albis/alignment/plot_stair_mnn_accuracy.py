#!/usr/bin/env python
"""
Figure 4C -- STAIR cross-section alignment accuracy.

Reads figure4c_mnn_domain_agreement.csv (written by plot_figure4c.py's
cross_slice_mnn_domain_agreement). Exact quantity computed there:

  * for each of the C(10, 2) = 45 pairs of tissue sections, each section is
    subsampled to 4000 observations;
  * a cross-section mutual-nearest-neighbour (MNN) pair is an observation i in
    section A and an observation j in section B that are each other's single
    nearest neighbour (Euclidean) in the 2D coordinates -- the *unaligned*
    coordinates for "unaligned", the *STAIR-aligned* coordinates for "stair";
  * the value for that section pair is the fraction of its MNN pairs with
    domain_true[i] == domain_true[j].

Each box below is the distribution of that value over the 45 section pairs.
No clustering, no chance-adjustment (unlike ARI) -- a plain percent-agreement
metric.

Manuscript wording (chosen 2026-09-09): the axis / Methods / caption call an
MNN pair a "matched spot" and the metric "% of matched spots sharing a true
domain". "Matched spot" is defined once as a cross-section mutual nearest
neighbour; note it is loose for the bin16um and cell resolutions (bins /
cells, not spots).

Writes sim_paper/data/figure_4/alignment/plots/stair_mnn_accuracy.png.
Env: sim_paper/env/albis-tutorial (or any with pandas + matplotlib).
"""
import os

import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

PLOTS_DIR = "/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_4/alignment/plots"
CSV = os.path.join(PLOTS_DIR, "figure4c_mnn_domain_agreement.csv")

MODALITIES = ["bin16um", "spot", "cell"]
MOD_DISPLAY = {"bin16um": "Bin (16 µm)", "spot": "Spot", "cell": "Cell"}
METHOD_LABELS = {"unaligned": "Unaligned", "stair": "STAIR aligned"}
# baseline (unaligned) = muted grey-blue, method (STAIR) = deeper teal-blue.
# Shared with plot_stagate_fig4b.py so the Figure 4 panels read as one.
METHOD_COLORS = {"unaligned": "#B7C0C8", "stair": "#2F6F8F"}
EDGE = "#3a3f44"


def _grouped_boxplot(ax, df):
    width = 0.35
    for mi, method in enumerate(("unaligned", "stair")):
        boxes = [100 * df.loc[(df.dataset == ds) & (df.method == method), "agreement"].values
                 for ds in MODALITIES]
        bp = ax.boxplot(boxes, positions=[p + (mi - 0.5) * width for p in range(len(MODALITIES))],
                        widths=width * 0.9, patch_artist=True, manage_ticks=False)
        for patch in bp["boxes"]:
            patch.set(facecolor=METHOD_COLORS[method], edgecolor=EDGE, linewidth=0.8)
        for line in bp["whiskers"] + bp["caps"]:
            line.set_color(EDGE)
        for med in bp["medians"]:
            med.set(color=EDGE, linewidth=1.4)
        for fly in bp["fliers"]:
            fly.set(marker="o", markersize=3, markerfacecolor="none", markeredgecolor=EDGE)
    ax.set_xticks(range(len(MODALITIES)))
    ax.set_xticklabels([MOD_DISPLAY[m] for m in MODALITIES])
    ax.spines[["top", "right"]].set_visible(False)


def plot(df, out_png):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    _grouped_boxplot(ax, df)
    ax.set_ylim(0, 100)
    ax.set_ylabel("% of matched spots sharing a true domain")
    ax.grid(axis="y", ls=":", alpha=0.5)
    ax.set_axisbelow(True)

    handles = [plt.Rectangle((0, 0), 1, 1, fc=METHOD_COLORS[m], ec=EDGE, lw=0.8,
                             label=METHOD_LABELS[m]) for m in ("unaligned", "stair")]
    ax.legend(handles=handles, frameon=False, loc="upper left",
              bbox_to_anchor=(1.0, 1.0), borderaxespad=0)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200, facecolor="white")
    plt.close(fig)
    print("wrote", out_png)


def main():
    df = pd.read_csv(CSV)
    summary = df.groupby(["dataset", "method"])["agreement"].agg(
        n="count", median="median",
        q1=lambda s: s.quantile(.25), q3=lambda s: s.quantile(.75))
    for c in ("median", "q1", "q3"):
        summary[c] = (100 * summary[c]).round(1)          # -> percent
    print(summary.to_string())
    plot(df, os.path.join(PLOTS_DIR, "stair_mnn_accuracy.png"))


if __name__ == "__main__":
    main()
