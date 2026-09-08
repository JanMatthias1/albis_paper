#!/usr/bin/env python
"""
Figure 4B -- STAGATE spatial-domain recovery, 2D vs 3D spatial graph, as a
function of the per-section batch effect.

Trimmed cut of plot_stagate_summary.py's grid: the strong domain-composition
row only (one column per resolution), x-axis = batch condition

    no batch effect   STAGATE run on the pre-batch counts
                      (3D_stagate.py --use-pre-batch, layers['counts_pre_batch'])
    sigma = 0.05      per-section batch shift at SD 0.05
    tuned             per-resolution batch SD used elsewhere in the study
                      (bin16um 0.7 / spot 0.3 / cell 1.5)

The weak domain-composition condition (domain ARI <= 0.03 for every
resolution / batch level) is stated in text, not drawn.

Reads <STAGATE>/<dataset>/adata_results/metrics.json (written by 3D_stagate.py).
Writes <STAGATE>/stagate_fig4b_batch.png and stagate_fig4b_strongmix.csv.

Env: sim_paper/env/albis-tutorial (or any with pandas + matplotlib).
"""
import glob
import json
import os

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

STAGATE_DIR = ("/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_4/"
               "spatial_clustering/STAGATE")
TUNED_SIGMA = {"bin16um": 0.7, "spot": 0.3, "cell": 1.5}
MODALITIES = ["bin16um", "spot", "cell"]
MOD_DISPLAY = {"bin16um": "Bin (16 µm)", "spot": "Spot", "cell": "Cell"}
BATCH_ORDER = {"none": 0, "0.05": 1, "tuned": 2}

# baseline (2D graph) = muted grey-blue, method (3D graph) = deeper teal-blue.
# Shared with plot_stair_mnn_accuracy.py so the Figure 4 panels read as one.
C_3D, C_2D, EDGE = "#2F6F8F", "#B7C0C8", "#3a3f44"


def parse_dataset(name):
    """<modality>[_strongmix][_lowbatch|_prebatch] -> (modality, mix, batch, sigma)."""
    modality = next((m for m in MODALITIES if name.startswith(m)), None)
    if modality is None:
        return None
    mix = "strong" if "strongmix" in name else "weak"
    if name.endswith("_prebatch"):
        batch, sigma = "none", 0.0
    elif "lowbatch" in name:
        batch, sigma = "0.05", 0.05
    else:
        batch, sigma = "tuned", TUNED_SIGMA[modality]
    return modality, mix, batch, sigma


def load():
    rows = []
    for mj in sorted(glob.glob(os.path.join(STAGATE_DIR, "*/adata_results/metrics.json"))):
        parsed = parse_dataset(mj.split(os.sep)[-3])
        if parsed is None:
            continue
        modality, mix, batch, sigma = parsed
        m = json.load(open(mj))
        rows.append(dict(
            dataset=mj.split(os.sep)[-3], modality=modality, mix=mix,
            batch=batch, batch_sigma=sigma,
            ari_domain_3d=m.get("stagate_3d", {}).get("ari_domain_true"),
            ari_domain_2d=m.get("stagate_2d", {}).get("ari_domain_true"),
            n_obs=m.get("n_obs"),
        ))
    df = pd.DataFrame(rows)
    df["_ord"] = df.batch.map(BATCH_ORDER)
    return df.sort_values(["modality", "mix", "_ord"]).reset_index(drop=True)


def _xtick(row):
    if row.batch == "none":
        return "no batch\neffect"
    if row.batch == "0.05":
        return r"$\sigma$=0.05"
    return f"tuned\n$\\sigma$={row.batch_sigma:g}"


def plot(strong, out_png):
    """strong: the mix == 'strong' rows only, one per (modality, batch)."""
    fig, axes = plt.subplots(1, len(MODALITIES), figsize=(11, 3.9), sharey=True)
    for ax, modality in zip(axes, MODALITIES):
        sub = strong[strong.modality == modality].sort_values("_ord")
        x = np.arange(len(sub))
        w = 0.38
        for offset, col, colour, label in [
            (-w / 2, "ari_domain_3d", C_3D, "STAGATE-3D"),
            (+w / 2, "ari_domain_2d", C_2D, "STAGATE-2D"),
        ]:
            bars = ax.bar(x + offset, sub[col], w, label=label, color=colour,
                          edgecolor=EDGE, linewidth=0.8)
            ax.bar_label(bars, fmt="%.2f", fontsize=8, padding=2)
        ax.set_xticks(x)
        ax.set_xticklabels([_xtick(r) for _, r in sub.iterrows()], fontsize=9)
        ax.set_ylim(0, 0.8)
        ax.set_title(MOD_DISPLAY[modality], fontsize=11)
        ax.grid(axis="y", ls=":", alpha=0.5)
        ax.set_axisbelow(True)
        ax.spines[["top", "right"]].set_visible(False)

    axes[0].set_ylabel("spatial-domain ARI (mclust vs. true domains)", fontsize=9)
    # Cell bars top out ~0.26, so the legend fits in its own empty upper area.
    axes[-1].legend(fontsize=8, frameon=False, loc="upper right", borderaxespad=0.6)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print("wrote", out_png)


def main():
    df = load()
    strong = df[df.mix == "strong"].copy()

    csv = os.path.join(STAGATE_DIR, "stagate_fig4b_strongmix.csv")
    strong.drop(columns="_ord").to_csv(csv, index=False)
    print("wrote", csv)
    print(strong[["modality", "batch", "batch_sigma",
                  "ari_domain_3d", "ari_domain_2d"]].to_string(index=False))

    weak_max = df.loc[df.mix == "weak", ["ari_domain_3d", "ari_domain_2d"]].max().max()
    print(f"\nweak domain-composition: max domain ARI across all conditions = {weak_max:.3f}")

    plot(strong, os.path.join(STAGATE_DIR, "stagate_fig4b_batch.png"))


if __name__ == "__main__":
    main()
