"""
Cross-modality check -- does slice_id=5 look the same across technologies
when EACH technology is aligned entirely on its own (the existing
per-technology 3D_stair.py z-stack runs used for Figure 4C), with no
cross-tech training at all? This is the complement to cross_tech_stair.py
(Figure 4E), which jointly trains STAIR across all three technologies
together -- here each technology never sees the other two during alignment.

Reads the three existing per-technology 3D_stair.py outputs
    sim_paper/data/figure_4/alignment/STAIR/<tech>/adata_results/Sim_3D_STAIR_<tech>.h5ad
subsets each to obs['slice_id'] == <slice>, isotropically rescales
non-reference technologies onto the reference technology's disc radius
(same rescale_to_ref logic as cross_tech_stair.py -- bin16um and spot
already share a radius, only cell needs it), then for display Procrustes-
fits each technology's OWN transform_fine onto its OWN (rescaled)
spatial_true. Because spatial_true for all three technologies already sits
in the same absolute simulated frame (differing only by that isotropic
scale -- verified in cross_tech_stair.py's docstring, ~87% domain-NN match
after rescale), overlaying the three per-tech display fits afterwards is a
meaningful comparison, not just three independently-rotated point clouds.

Produces, under sim_paper/data/figure_4/cross_modality_alignment/plots/ :
    slice5_independent_crossmodality_<colorby>.png   cols = [unaligned |
        independently STAIR-aligned | ground truth], all 3 technologies
        overlaid, 2D

--order controls the technology draw order (last = on top); default keeps
the sparse spot points visible above cell/bin16um.
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

FIG4C_STAIR = "/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_4/alignment/STAIR"
OUTDIR_DEFAULT = "/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_4/cross_modality_alignment/plots"
COL_LABELS = ["Unaligned", "Independently STAIR-aligned\n(no cross-tech training)", "Ground truth"]
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from manuscript_style import MODALITY_LOOKUP, category_color
TECH_COLORS = {key: MODALITY_LOOKUP[key] for key in ('bin16um', 'spot', 'cell')}
TECHS = ["bin16um", "spot", "cell"]
REF_TECH = "bin16um"


def procrustes_fit(X, Y):
    """Best rigid map X -> Y (rotation + translation, no scale/reflection)."""
    X = np.asarray(X, float)
    Y = np.asarray(Y, float)
    mx, my = X.mean(0), Y.mean(0)
    U, _, Vt = np.linalg.svd((X - mx).T @ (Y - my))
    R = U @ Vt
    if np.linalg.det(R) < 0:
        U[:, -1] *= -1
        R = U @ Vt
    return (X - mx) @ R + my


def load_tech_slice(tech, slice_id):
    p = os.path.join(FIG4C_STAIR, tech, "adata_results", f"Sim_3D_STAIR_{tech}.h5ad")
    a = ad.read_h5ad(p)
    m = (a.obs["slice_id"].astype(str) == str(slice_id)).values
    if m.sum() == 0:
        raise ValueError(f"{tech}: no obs at slice_id={slice_id}")
    return a[m].copy()


def rescale_to_ref(tech_adatas, ref_tech):
    """Isotropic rescale of every non-ref technology's spatial coords onto the
    ref technology's disc radius (max |xy| over this slice's spatial_true).
    Mirrors cross_tech_stair.py's rescale_to_ref, applied post-hoc to the
    independently-aligned 4C outputs instead of pre-training inputs."""
    ref_radius = np.abs(np.asarray(tech_adatas[ref_tech].obsm["spatial_true"])[:, :2]).max()
    for tech, a in tech_adatas.items():
        if tech == ref_tech:
            continue
        radius = np.abs(np.asarray(a.obsm["spatial_true"])[:, :2]).max()
        scale = float(ref_radius / radius)
        for key in ("spatial_true", "spatial_unaligned", "transform_fine"):
            xy = np.asarray(a.obsm[key], dtype=float).copy()
            xy[:, :2] *= scale
            a.obsm[key] = xy
        print(f"[rescale] {tech}: radius {radius:.1f} -> {ref_radius:.1f} um (x{scale:.4f})")
    return tech_adatas


def panels_for(tech_adatas):
    """unaligned / independently-STAIR-aligned (display-Procrustes'd onto own
    truth) / ground truth, each as a dict tech -> (n_tech, 2) array."""
    unal, fine, true = {}, {}, {}
    for tech, a in tech_adatas.items():
        t = np.asarray(a.obsm["spatial_true"], float)[:, :2]
        u = np.asarray(a.obsm["spatial_unaligned"], float)[:, :2]
        f = np.asarray(a.obsm["transform_fine"], float)[:, :2]
        unal[tech], true[tech] = u, t
        fine[tech] = procrustes_fit(f, t)
    return [unal, fine, true]


def plot_overlay(tech_adatas, color_by, order, outdir, slice_id):
    cols3 = panels_for(tech_adatas)
    if color_by == "technology":
        cats = TECHS
        lut = {c: mcolors.to_rgba(TECH_COLORS[c]) for c in cats}
    else:
        all_vals = pd.concat([tech_adatas[t].obs[color_by].astype(str) for t in TECHS])
        cats = sorted(all_vals.unique())
        lut = {c: mcolors.to_rgba(category_color(c, color_by)) for c in cats}

    fig, axes = plt.subplots(1, 3, figsize=(12, 4.3), squeeze=False)
    for c in range(3):
        ax = axes[0][c]
        for tech in order:
            xy = cols3[c][tech]
            n = xy.shape[0]
            idx = np.arange(n) if n <= 60000 else np.sort(
                np.random.default_rng(0).choice(n, 60000, replace=False))
            if color_by == "technology":
                colors = [lut[tech]] * len(idx)
            else:
                vals = tech_adatas[tech].obs[color_by].astype(str).values[idx]
                colors = [lut[v] for v in vals]
            if tech == "spot":
                # spot is far sparser per unit area than bin16um/cell -- small
                # semi-transparent markers blend into the denser layers underneath
                # and read as "missing" even though every point is there.
                ax.scatter(xy[idx, 0], xy[idx, 1], s=14, c=colors, linewidths=0.3,
                           edgecolors="black", alpha=1.0)
            else:
                s = 1.0 if n > 20000 else 8.0
                ax.scatter(xy[idx, 0], xy[idx, 1], s=s, c=colors, linewidths=0, alpha=0.85)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(COL_LABELS[c], fontsize=12)
    handles = [Line2D([0], [0], marker="o", ls="", mfc=lut[c], mec="none",
                      label=str(c)) for c in cats]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.0),
               borderaxespad=0.2, ncol=min(len(cats), 10), frameon=False,
               title=color_by)
    fig.patch.set_facecolor("white")
    fig.tight_layout(rect=[0, 0.05, 1, 0.92])
    outpath = os.path.join(outdir, f"slice{slice_id}_independent_crossmodality_{color_by}.png")
    fig.savefig(outpath, dpi=200, facecolor="white")
    plt.close(fig)
    print("wrote", outpath)


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--slice", type=int, default=5)
    ap.add_argument("--color-by", nargs="+", default=["technology", "domain_true"],
                    choices=["technology", "domain_true", "cell_type_true"])
    ap.add_argument("--order", nargs="+", default=["cell", "bin16um", "spot"],
                    choices=TECHS, help="draw order, last = on top (default keeps "
                                        "sparse spot points visible above cell/bin16um)")
    ap.add_argument("--outdir", default=OUTDIR_DEFAULT)
    return ap.parse_args()


if __name__ == "__main__":
    args = parse_args()
    assert set(args.order) == set(TECHS), "--order must list all three technologies"
    os.makedirs(args.outdir, exist_ok=True)

    tech_adatas = {t: load_tech_slice(t, args.slice) for t in TECHS}
    for t, a in tech_adatas.items():
        print(f"{t}: {a.shape} at slice_id={args.slice}")
    tech_adatas = rescale_to_ref(tech_adatas, REF_TECH)

    for cb in args.color_by:
        plot_overlay(tech_adatas, cb, args.order, args.outdir, args.slice)
