"""
Per-slice spatial view of a Figure 2 modality's QC'd data -- one panel per
slice_id (0-9), so the off-tissue exclusion / pole-thinning pattern
(see DATA_VERSIONS.md, "realwindow is now the Figure 2 default" +
"Known issue: bin/spot zero-inflation from sparse tissue packing") is
visible slice by slice rather than pooled into one overlaid scatter.

Reads simulation_<modality>_z_qc.h5ad for the given tag and plots
obsm['spatial'], colored by domain_true, one subplot per slice_id in a
2x5 grid.

Usage:
    python plot_qc_slices.py --modality spot \
        --tag packing_pf0p04_log_mu_-2.0_theta_0.25_jitter0.10_bsigma03
"""

import argparse
import os

import numpy as np
import pandas as pd
import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D

FIG2 = "/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_2"
PASTEL = ["#4878d0", "#ee854a", "#6acc64", "#d65f5f", "#956cb4",
          "#8c613c", "#dc7ec0", "#797979", "#d5bb67", "#82c6e2"]

DEFAULT_TAGS = {
    "bin16um": "packing_pf0p04_bin16um_log_mu_-2.5_bsigma07",
    "spot": "packing_pf0p04_log_mu_-2.0_theta_0.25_jitter0.10_bsigma03",
    "cell": "log_mu_-2.3_theta_0.40_jitter0.15_bsigma15",
}
MODALITY_FILE_PREFIX = {"bin16um": "bin", "spot": "spot", "cell": "cell"}


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--modality", required=True, choices=list(DEFAULT_TAGS))
    ap.add_argument("--tag", default=None, help="defaults to the current canonical Fig2 tag")
    ap.add_argument("--color-by", default="domain_true", choices=["domain_true", "cell_type_true"])
    ap.add_argument("--outdir", default=None, help="defaults to <tag_dir>/plots")
    return ap.parse_args()


if __name__ == "__main__":
    args = parse_args()
    tag = args.tag or DEFAULT_TAGS[args.modality]
    tag_dir = os.path.join(FIG2, tag)
    prefix = MODALITY_FILE_PREFIX[args.modality]
    h5ad_path = os.path.join(tag_dir, f"simulation_{prefix}_z_qc.h5ad")
    outdir = args.outdir or os.path.join(tag_dir, "plots")
    os.makedirs(outdir, exist_ok=True)

    a = ad.read_h5ad(h5ad_path)
    print(f"loaded {h5ad_path}: {a.shape}")

    slices = sorted(a.obs["slice_id"].unique(), key=lambda v: int(v))
    cats = list(pd.Categorical(a.obs[args.color_by].astype(str)).categories)
    lut = {c: mcolors.to_rgba(PASTEL[i % len(PASTEL)]) for i, c in enumerate(cats)}

    n = len(slices)
    ncols = 5
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.2 * ncols, 3.2 * nrows), squeeze=False)

    for i, sl in enumerate(slices):
        ax = axes[i // ncols][i % ncols]
        m = (a.obs["slice_id"].astype(str) == str(sl)).values
        xy = np.asarray(a.obsm["spatial"])[m, :2]
        colors = a.obs[args.color_by].astype(str).values[m]
        colors = np.vstack([lut[c] for c in colors])
        s = 2.0 if m.sum() > 5000 else 10.0
        ax.scatter(xy[:, 0], xy[:, 1], s=s, c=colors, linewidths=0, alpha=0.85)
        ax.set_aspect("equal")
        ax.set_title(f"slice_id={sl}  (n={m.sum()})", fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
    for j in range(n, nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")

    handles = [Line2D([0], [0], marker="o", ls="", mfc=lut[c], mec="none", label=c) for c in cats]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.0),
               borderaxespad=0.2, ncol=min(len(cats), 10), frameon=False, title=args.color_by)
    fig.suptitle(f"{args.modality} (QC'd) -- {tag}", fontsize=12)
    fig.patch.set_facecolor("white")
    fig.tight_layout(rect=[0, 0.04, 1, 0.95])

    outpath = os.path.join(outdir, f"qc_slices_{args.modality}_{args.color_by}.png")
    fig.savefig(outpath, dpi=180, facecolor="white")
    plt.close(fig)
    print("wrote", outpath)
