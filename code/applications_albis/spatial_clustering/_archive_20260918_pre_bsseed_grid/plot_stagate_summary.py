#!/usr/bin/env python
"""
Figure 4B summary -- STAGATE 3D vs 2D spatial-domain recovery across the batch grid.

Scans every  <STAGATE>/<dataset>/adata_results/metrics.json  written by
3D_stagate.py and draws domain-ARI (STAGATE-3D vs STAGATE-2D) as a function of
batch strength, one column per modality, split by domain-mix strength.

    batch conditions, strong mix:  none (oracle: counts_pre_batch, --use-pre-batch)
                                   -> sigma 0.05 (very low)  ->  tuned (Figure 2 sigma)
    batch conditions, weak mix:    sigma 0.05  ->  tuned

Output (both under <STAGATE>/):
    stagate_batch_comparison.png
    stagate_results_summary.csv
"""
import glob
import json
import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

STAGATE_DIR = ("/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_4/"
               "spatial_clustering/STAGATE")

# tuned batch_sigma per modality (Figure 2 finalized values; see 3D_stagate.py)
TUNED_SIGMA = {"bin16um": 0.7, "spot": 0.3, "cell": 1.5}
MODALITIES = ["bin16um", "spot", "cell"]
BATCH_ORDER = {"none": 0, "0.05": 1, "tuned": 2}

C3D, C2D = "#2166ac", "#b2182b"


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
        ds = mj.split(os.sep)[-3]
        parsed = parse_dataset(ds)
        if parsed is None:
            continue
        modality, mix, batch, sigma = parsed
        m = json.load(open(mj))
        s3, s2 = m.get("stagate_3d", {}), m.get("stagate_2d", {})
        rows.append(dict(
            dataset=ds, modality=modality, mix=mix, batch=batch, batch_sigma=sigma,
            ari_domain_3d=s3.get("ari_domain_true"), ari_domain_2d=s2.get("ari_domain_true"),
            nmi_domain_3d=s3.get("nmi_domain_true"), nmi_domain_2d=s2.get("nmi_domain_true"),
            ari_celltype_3d=s3.get("ari_cell_type_true"), ari_celltype_2d=s2.get("ari_cell_type_true"),
            n_obs=m.get("n_obs"), graph_model=m.get("graph_model"),
            deg_2d=m.get("deg_2d"), deg_3d=m.get("deg_3d"),
        ))
    df = pd.DataFrame(rows)
    df["_ord"] = df.batch.map(BATCH_ORDER)
    return df.sort_values(["modality", "mix", "_ord"]).reset_index(drop=True)


def _xlabel(row):
    if row.batch == "none":
        return "none\n(oracle)"
    if row.batch == "0.05":
        return r"$\sigma$=0.05"
    return f"tuned\n$\\sigma$={row.batch_sigma:g}"


def plot(df, out_png):
    fig, axes = plt.subplots(2, 3, figsize=(13, 8), sharey=True)
    for j, modality in enumerate(MODALITIES):
        for i, mix in enumerate(["strong", "weak"]):
            ax = axes[i, j]
            sub = df[(df.modality == modality) & (df.mix == mix)].sort_values("_ord")
            x = np.arange(len(sub))
            w = 0.38
            b3 = ax.bar(x - w / 2, sub.ari_domain_3d, w, label="STAGATE-3D", color=C3D)
            b2 = ax.bar(x + w / 2, sub.ari_domain_2d, w, label="STAGATE-2D", color=C2D)
            ax.bar_label(b3, fmt="%.2f", fontsize=7, padding=2)
            ax.bar_label(b2, fmt="%.2f", fontsize=7, padding=2)
            ax.set_xticks(x)
            ax.set_xticklabels([_xlabel(r) for _, r in sub.iterrows()], fontsize=8)
            ax.set_ylim(0, 1)
            ax.set_title(f"{modality}  |  {mix} domain mix", fontsize=10)
            ax.grid(axis="y", ls=":", alpha=0.5)
            ax.set_axisbelow(True)
            if j == 0:
                ax.set_ylabel("domain ARI  (mclust vs domain_true)")
            if i == 0 and j == 0:
                ax.legend(fontsize=8, loc="upper left", framealpha=0.9)
    fig.suptitle("STAGATE 3D vs 2D spatial-domain recovery  —  batch strength x domain-mix strength",
                 fontsize=12)
    fig.text(0.5, 0.005,
             "batch axis: none = counts_pre_batch (oracle ceiling)  ->  sigma 0.05  ->  "
             "tuned Figure-2 sigma (bin16um 0.7 / spot 0.3 / cell 1.5)",
             ha="center", fontsize=8, color="0.35")
    fig.tight_layout(rect=[0, 0.03, 1, 0.96])
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    print("wrote", out_png)


def main():
    df = load()
    csv = os.path.join(STAGATE_DIR, "stagate_results_summary.csv")
    df.drop(columns="_ord").to_csv(csv, index=False)
    print("wrote", csv)
    print(df[["dataset", "modality", "mix", "batch", "batch_sigma",
              "ari_domain_3d", "ari_domain_2d"]].to_string(index=False))
    plot(df, os.path.join(STAGATE_DIR, "stagate_batch_comparison.png"))


if __name__ == "__main__":
    main()
