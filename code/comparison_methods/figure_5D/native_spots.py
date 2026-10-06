"""Figure 5D: each method's own spots on slice 5 (stored slice_id 4).

Central 290 x 290 um field of each tissue. Everything drawn is native output:
- ALBIS: native molecules (cached by actual_capture.py) and ALBIS's own spots from
  figure_5A_600k/albis/spot.h5ad; the label is the number of molecules within the
  spot radius, asserted equal to the spot's saved pre-batch total.
- SPIDER: cells and circular spots from the seed-20260922 run (--spider-base).
  Membership = cells within the spot radius, which generate.py asserted equals
  SPIDER's own membership matrix; the label is asserted equal to n_source_cells.
- scCube: Figure 5A cells and spots; memberships from sccube_slice_membership.py
  (scCube's own generate_spot_data_random, reproduces the saved spots). Each spot is
  the convex hull of its member cells; the label is the spot's cell count.
Run with sim_paper/env/albis-tutorial/bin/python.
"""
import argparse
import json
from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.colors import to_rgb
from matplotlib.patches import Circle, Patch, Polygon
from scipy.spatial import ConvexHull

PAPER = Path(__file__).resolve().parents[3]
BASE = PAPER / "data/figure_5/figure_5A_600k"
ALBIS_DIR = PAPER / "data/figure_5/figure_5D_actual_coordinates"  # shared ALBIS molecule cache
SLICE_ID, HALF, OFFSET = 4, 145.0, 2050.0  # field half-width; scCube/SPIDER cube centre

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--spider-base", type=Path, required=True, help="run folder with spider/{cell,spot}.h5ad")
parser.add_argument("--out", type=Path, required=True, help="folder with sccube_membership.csv.gz; figure goes here")
args = parser.parse_args()
cfg = json.loads((BASE / "settings.json").read_text())
palette, radius = cfg["celltype_colors"], cfg["spot_radius_um"]


def in_field(xy, margin=0.0):
    return (np.abs(xy[:, 0]) <= HALF + margin) & (np.abs(xy[:, 1]) <= HALF + margin)


def slice_obs(path, offset):
    a = ad.read_h5ad(path, backed="r")
    keep = a.obs.slice_id.to_numpy().astype(int) == SLICE_ID
    xy = np.asarray(a.obsm["spatial"])[keep] - offset
    obs = a.obs[keep].copy()
    layer = a.layers["counts_pre_batch"][np.where(keep)[0]] if "counts_pre_batch" in a.layers else None
    a.file.close()
    return xy, obs, layer


panels, reports = {}, []

# ALBIS: molecules within the radius of ALBIS's own spots.
mol = np.load(ALBIS_DIR / "albis_molecules_roi.npz")
mol_xy, mol_type = mol["xyz"][:, :2], np.array([f"type{i + 1}" for i in mol["source_type"]])
spot_xy, spot_obs, pre = slice_obs(BASE / "albis/spot.h5ad", 0.0)
# ALBIS's spot grid is anchored to its capture window, so no spot sits at the tissue
# centre; centre ALBIS's field on its spot nearest the centre so it shows the same
# 3 x 3 block as SPIDER (the fields then sit at different tissue positions).
albis_centre = spot_xy[np.argmin(np.abs(spot_xy).sum(1))].copy()
spot_xy, mol_xy = spot_xy - albis_centre, mol_xy - albis_centre
f = in_field(spot_xy, -radius)  # whole capture inside the field (the molecule cache covers +-200 um)
spot_xy, pre_totals = spot_xy[f], np.asarray(pre[np.where(f)[0]].sum(axis=1)).ravel()
d2 = ((mol_xy[:, None, :] - spot_xy[None, :, :]) ** 2).sum(2)
inside = d2.min(1) <= radius ** 2
counts = (d2 <= radius ** 2).sum(0)
assert np.array_equal(counts, pre_totals.astype(int)), (counts, pre_totals)
panels["ALBIS"] = dict(points=mol_xy, types=mol_type, inside=inside, spots=spot_xy, counts=counts)
reports.append(dict(method="ALBIS", field_centre_um=albis_centre.tolist(), spots=spot_xy.tolist(), counts=counts.tolist(),
                    check="molecules within radius == saved counts_pre_batch totals"))

# SPIDER: cells within the radius of SPIDER's own circular spots.
cell_xy, cell_obs, _ = slice_obs(args.spider_base / "spider/cell.h5ad", OFFSET)
spot_xy, spot_obs, _ = slice_obs(args.spider_base / "spider/spot.h5ad", OFFSET)
f = in_field(spot_xy)
spot_xy, saved_n = spot_xy[f], spot_obs.n_source_cells.to_numpy()[f]
d2 = ((cell_xy[:, None, :] - spot_xy[None, :, :]) ** 2).sum(2)
counts = (d2 <= radius ** 2).sum(0)
assert np.array_equal(counts, saved_n), (counts, saved_n)
c = in_field(cell_xy)
panels["SPIDER"] = dict(points=cell_xy[c], types=cell_obs.cell_type_true.astype(str).to_numpy()[c],
                        inside=(d2.min(1) <= radius ** 2)[c], spots=spot_xy, counts=counts)
reports.append(dict(method="SPIDER", spots=spot_xy.tolist(), counts=counts.tolist(),
                    check="cells within radius == saved n_source_cells"))

# scCube: native memberships (every cell belongs to one spot). scCube assigns cells to
# square grid tiles; platform="Visium" then shifts every second column's reported centre
# by half a spacing without changing membership. Each spot is drawn as the convex hull
# of its own member cells (native membership only; scCube's grid is not re-derived).
cell_xy, cell_obs, _ = slice_obs(BASE / "sccube/cell.h5ad", OFFSET)
spots = pd.read_csv(args.out / "sccube_spots.csv").set_index("spot")
member = pd.read_csv(args.out / "sccube_membership.csv.gz").set_index("Cell").spot
spot_of = member.reindex(cell_obs.index).to_numpy()
spot_xy = spots[["x_um", "y_um"]].to_numpy() - OFFSET
f = in_field(spot_xy)
c = in_field(cell_xy)
tiles, tile_counts = [], []
for name in pd.unique(spot_of[c]):  # every spot with a member cell in the field
    pts = cell_xy[spot_of == name]
    tiles.append(pts[ConvexHull(pts).vertices] if len(pts) >= 3 else pts)
    tile_counts.append((pts.mean(0), int(spots.loc[name, "n_cells"])))
panels["scCube"] = dict(points=cell_xy[c], types=cell_obs.cell_type_true.astype(str).to_numpy()[c],
                        inside=np.ones(c.sum(), bool), spots=spot_xy[f], counts=spots.n_cells.to_numpy()[f],
                        tiles=tiles, tile_counts=tile_counts)
reports.append(dict(method="scCube", n_spots_in_field=int(f.sum()), counts=spots.n_cells.to_numpy()[f].tolist(),
                    check="memberships from scCube generate_spot_data_random reproduce saved Figure 5A spots"))

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "pdf.fonttype": 42, "svg.fonttype": "none"})
fig, axes = plt.subplots(1, 3, figsize=(13, 5.8))
subtitles = {"ALBIS": "Native mRNA instances and spots",
             "scCube": f"Native cells and spots (target {cfg['sccube_cells_per_spot']} cells per spot)",
             "SPIDER": "Native cells and spots"}
for ax, name in zip(axes, ["ALBIS", "scCube", "SPIDER"]):
    p = panels[name]
    albis = name == "ALBIS"
    if "tiles" in p:
        ax.add_collection(PatchCollection([Polygon(t, closed=True) for t in p["tiles"]], facecolor="#F4F4F4",
                                          edgecolor="#9E9E9E", linewidths=.8, zorder=1))
    out = ~p["inside"]
    ax.scatter(*p["points"][out].T, s=3 if albis else 10, c="#BDBDBD", alpha=.35, marker="x" if albis else "o",
               linewidths=.35 if albis else 0, rasterized=True, zorder=2)
    ax.scatter(*p["points"][~out].T, s=5 if albis else 20, c=[palette[t] for t in p["types"][~out]],
               marker="x" if albis else "o", linewidths=.45 if albis else 0, rasterized=True, zorder=3)
    for centre, n in p.get("tile_counts", []):  # scCube: cell count inside each tile within the field
        if np.all(np.abs(centre) <= HALF - 8):
            ax.text(*centre, str(n), ha="center", va="center", fontsize=8, fontweight="bold", zorder=5, clip_on=True,
                    bbox=dict(boxstyle="round,pad=.15", facecolor="white", edgecolor="none", alpha=.8))
    for centre, n in zip(p["spots"], p["counts"]):
        if name == "scCube":  # scCube's reported spot centre (every second column shifted for "Visium")
            ax.plot(*centre, marker="+", color="#333333", ms=6, mew=1, zorder=4)
        else:
            ax.add_patch(Circle(centre, radius, fill=False, edgecolor="#333333", lw=1.2, zorder=4))
            ax.text(centre[0], centre[1] - radius - 4, str(n), ha="center", va="top", fontsize=9, zorder=5)
    ax.text(.5, 1.10, name, transform=ax.transAxes, ha="center", va="bottom", fontweight="bold", fontsize=15)
    ax.set_title(subtitles[name], fontweight="bold", fontsize=11, pad=12)
    ax.set(xlim=(-HALF, HALF), ylim=(-HALF, HALF), aspect="equal")
    ax.set_axis_off()
    # scale bar just below the field, clear of the spot-count labels
    ax.plot([-HALF, -HALF + 50], [-HALF - 10, -HALF - 10], color="black", lw=2, clip_on=False)
    ax.text(-HALF + 25, -HALF - 14, "50 µm", ha="center", va="top", fontsize=9, clip_on=False)
fig.legend(handles=[Patch(color=palette[f"type{i}"], label=f"Type {i}") for i in range(1, 9)],
           loc="lower center", bbox_to_anchor=(.5, .03), ncol=8, frameon=False, fontsize=10)
fig.subplots_adjust(left=.025, right=.975, top=.80, bottom=.20, wspace=.15)
for ext in ["png", "pdf", "svg"]:
    fig.savefig(args.out / f"figure5d_native_spots.{ext}", dpi=500, bbox_inches="tight", facecolor="white")
(args.out / "provenance.json").write_text(json.dumps(dict(
    albis_and_sccube=str(BASE), spider=str(args.spider_base / "spider"), slice_id=SLICE_ID,
    field_um=[-HALF, HALF], spot_radius_um=radius,
    labels="ALBIS: pre-batch molecules per spot; SPIDER/scCube: cells per spot",
    reports=reports), indent=2) + "\n")
print(json.dumps(reports), flush=True)
