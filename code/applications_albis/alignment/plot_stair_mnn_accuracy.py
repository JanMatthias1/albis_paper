#!/usr/bin/env python
"""
STAIR domain-recovery accuracy, unaligned vs STAIR-aligned -- cross-slice
mutual-nearest-neighbor (MNN) percentage agreement.

Reads figure4c_mnn_domain_agreement.csv (written by plot_figure4c.py's
cross_slice_mnn_domain_agreement): for every pair of the 10 sections, find
mutual nearest neighbors in the 2D spatial embedding (point i in section A
whose nearest neighbor in section B is point j, and vice versa), and report
what fraction of those matched pairs share the same domain_true label. One row
per (dataset, method, slice pair) -- so each box below is a real distribution
over the 45 slice pairs, not synthetic clustering-seed noise.

No clustering step, no chance-adjustment (unlike ARI): this is a plain percent
-agreement metric, easy to state as "of the spots that align across sections,
X% land in the same true domain."

Output, under sim_paper/data/figure_4/alignment/plots/ :
    stair_mnn_accuracy.png
"""
import os

import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

PLOTS_DIR = "/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_4/alignment/plots"
CSV = os.path.join(PLOTS_DIR, "figure4c_mnn_domain_agreement.csv")
MODALITIES = ["bin16um", "spot", "cell"]
METHOD_LABELS = {"unaligned": "Unaligned", "stair": "STAIR aligned"}
METHOD_COLORS = {"unaligned": "#d65f5f", "stair": "#4878d0"}  # match plot_figure4c.py's METHOD_COLORS


def _grouped_boxplot(ax, df, datasets, methods):
    width = 0.35
    positions = range(len(datasets))
    for mi, method in enumerate(methods):
        boxes = [100 * df.loc[(df.dataset == ds) & (df.method == method), "agreement"].values
                 for ds in datasets]
        bp = ax.boxplot(boxes, positions=[p + (mi - 0.5) * width for p in positions],
                         widths=width * 0.9, patch_artist=True, manage_ticks=False)
        for patch in bp["boxes"]:
            patch.set_facecolor(METHOD_COLORS[method])
            patch.set_alpha(0.85)
        for med in bp["medians"]:
            med.set_color("black")
    ax.set_xticks(list(positions))
    ax.set_xticklabels(datasets)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot(df, out_png):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    _grouped_boxplot(ax, df, MODALITIES, ["unaligned", "stair"])

    ax.set_ylim(0, 100)
    ax.set_ylabel("% of matched spots sharing a true domain")
    ax.set_title("STAIR Alignment", fontsize=12)
    ax.grid(axis="y", ls=":", alpha=0.5)
    ax.set_axisbelow(True)

    handles = [plt.Rectangle((0, 0), 1, 1, fc=METHOD_COLORS[m], alpha=0.85,
                             label=METHOD_LABELS[m]) for m in ["unaligned", "stair"]]
    ax.legend(handles=handles, frameon=False, loc="upper left",
             bbox_to_anchor=(1.0, 1.0), borderaxespad=0)
    fig.patch.set_facecolor("white")
    fig.tight_layout()
    fig.savefig(out_png, dpi=200, facecolor="white")
    plt.close(fig)
    print("wrote", out_png)


def main():
    df = pd.read_csv(CSV)
    print((100 * df.groupby(["dataset", "method"])["agreement"].mean()).round(1))
    plot(df, os.path.join(PLOTS_DIR, "stair_mnn_accuracy.png"))


if __name__ == "__main__":
    main()
