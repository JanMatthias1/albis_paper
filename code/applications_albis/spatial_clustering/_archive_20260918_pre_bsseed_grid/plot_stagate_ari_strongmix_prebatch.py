#!/usr/bin/env python
"""
STAGATE domain-recovery ARI, strongmix_prebatch condition only.

Reads the existing metrics.json already written by 3D_stagate.py for the
<modality>_strongmix_prebatch runs (strong domain/expression mix, batch effect
removed via --use-pre-batch) -- no new STAGATE training, just plots what's on
disk under
    sim_paper/data/figure_4/spatial_clustering/STAGATE/<modality>_strongmix_prebatch/adata_results/metrics.json

Output, under sim_paper/data/figure_4/spatial_clustering/STAGATE/ :
    stagate_ari_strongmix_prebatch.png
    stagate_ari_strongmix_prebatch.csv
"""
import json
import os

import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

STAGATE_DIR = ("/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_4/"
               "spatial_clustering/STAGATE")
MODALITIES = ["bin16um", "spot", "cell"]
C3D, C2D = "#2166ac", "#b2182b"


def load():
    rows = []
    for modality in MODALITIES:
        p = os.path.join(STAGATE_DIR, f"{modality}_strongmix_prebatch",
                         "adata_results", "metrics.json")
        if not os.path.exists(p):
            print("[warn] missing", p)
            continue
        m = json.load(open(p))
        s3, s2 = m.get("stagate_3d", {}), m.get("stagate_2d", {})
        rows.append(dict(modality=modality, n_obs=m.get("n_obs"),
                         ari_domain_3d=s3.get("ari_domain_true"),
                         ari_domain_2d=s2.get("ari_domain_true")))
    return pd.DataFrame(rows)


def plot(df, out_png):
    x = range(len(df))
    w = 0.35
    fig, ax = plt.subplots(figsize=(6, 4.3))
    b3 = ax.bar([i - w / 2 for i in x], df.ari_domain_3d, w, label="STAGATE-3D", color=C3D)
    b2 = ax.bar([i + w / 2 for i in x], df.ari_domain_2d, w, label="STAGATE-2D", color=C2D)
    ax.bar_label(b3, fmt="%.2f", fontsize=9, padding=2)
    ax.bar_label(b2, fmt="%.2f", fontsize=9, padding=2)
    ax.set_xticks(list(x))
    ax.set_xticklabels(df.modality)
    ax.set_ylim(0, 1)
    ax.set_ylabel("domain ARI")
    ax.set_title("STAGATE Domain Recovery", fontsize=12)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", ls=":", alpha=0.5)
    ax.set_axisbelow(True)
    ax.legend(frameon=False)
    fig.patch.set_facecolor("white")
    fig.tight_layout()
    fig.savefig(out_png, dpi=200, facecolor="white")
    plt.close(fig)
    print("wrote", out_png)


def main():
    df = load()
    if df.empty:
        raise SystemExit("no strongmix_prebatch STAGATE metrics.json found")
    csv = os.path.join(STAGATE_DIR, "stagate_ari_strongmix_prebatch.csv")
    df.to_csv(csv, index=False)
    print("wrote", csv)
    print(df.to_string(index=False))
    plot(df, os.path.join(STAGATE_DIR, "stagate_ari_strongmix_prebatch.png"))


if __name__ == "__main__":
    main()
