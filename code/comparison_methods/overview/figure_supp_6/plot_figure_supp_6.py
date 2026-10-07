"""Supplementary Figure 6: scCube spot aggregation vs cells per spot, Figure 5A slice 5.

Native scCube only: generate_spot_data_random(platform="Visium", n_cell=target)
on the Figure 5A scCube cells of slice 5 (stored slice_id 4), then
calculate_spot_prop for each spot's cell-type proportions. Spots are coloured by
their most abundant type (ties -> lowest type number, as in Figure 5A).
Run with comparison_methods/env/sccube/bin/python.
"""
import argparse
import hashlib
import json
from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parents[5]
BASE = ROOT / "albis_paper/data/figure_5/figure_5A_600k"
OUT = ROOT / "albis_paper/data/figure_5/figure_supp_6"
TARGETS = [3, 100, 1000, 5000, 10000]
TYPES = [f"type{i}" for i in range(1, 9)]
SLICE_ID = 4


def native_spots(model, calculate_spot_prop, source, meta, target, scale):
    expr, loc, membership = model.generate_spot_data_random(
        source, meta.copy(), platform="Visium", gene_type="whole", n_cell=target)
    loc = loc.set_index("spot").loc[expr.columns]
    prop = calculate_spot_prop(membership).reindex(index=expr.columns, columns=TYPES, fill_value=0)
    n = membership.spot.value_counts().reindex(expr.columns, fill_value=0).to_numpy()
    assert n.sum() == len(meta) and membership.Cell.nunique() == len(meta)
    xy = loc[["spot_x", "spot_y"]].to_numpy(dtype=float) * scale
    return expr.columns.to_numpy(), xy, prop.to_numpy(), n, membership[["Cell", "spot"]]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--redraw-only", action="store_true", help="Redraw from saved spot tables.")
    args = parser.parse_args()
    OUT.mkdir(exist_ok=True)
    cfg = json.loads((BASE / "settings.json").read_text())
    scale = cfg["extent_um"] / cfg["sccube_grid_size"]
    colors = cfg["celltype_colors"]
    a = ad.read_h5ad(BASE / "sccube/cell.h5ad")
    a = a[a.obs.slice_id.to_numpy().astype(int) == SLICE_ID].copy()
    ids, labels = a.obs_names.to_numpy(), a.obs.cell_type_true.astype(str).to_numpy()
    native_xy = np.asarray(a.obsm["spatial_3d_native"])[:, :2]
    N = a.n_obs

    tables = {}
    if args.redraw_only:
        for t in TARGETS:
            tables[t] = pd.read_csv(OUT / f"spots_target{t}.csv")
    else:
        from scCube.sccube import scCube
        from scCube.utils import calculate_spot_prop
        model = scCube()
        # Same inputs generate.py gives scCube: its own grid coordinates, columns 2/3 = x/y.
        meta = pd.DataFrame({"Cell": ids, "Cell_type": labels,
                             "point_x": native_xy[:, 0], "point_y": native_xy[:, 1]})
        X = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)
        source = pd.DataFrame(X.T, index=a.var_names, columns=ids)
        for t in TARGETS:
            names, xy, prop, n, membership = native_spots(model, calculate_spot_prop, source, meta, t, scale)
            table = pd.DataFrame(dict(spot=names, x_um=xy[:, 0], y_um=xy[:, 1], n_cells=n))
            for k, name in enumerate(TYPES):
                table[name] = prop[:, k]
            tables[t] = table
            table.to_csv(OUT / f"spots_target{t}.csv", index=False)
            membership.to_csv(OUT / f"membership_target{t}.csv.gz", index=False)
            print(f"target {t}: {len(table)} spots", flush=True)
        # Consistency: target 10 must reproduce the saved Figure 5A slice-5 spots.
        _, xy10, prop10, _, _ = native_spots(model, calculate_spot_prop, source, meta, 10, scale)
        b = ad.read_h5ad(BASE / "sccube/spot.h5ad", backed="r")
        sel = b.obs.slice_id.to_numpy().astype(int) == SLICE_ID
        saved_xy, saved_prop = np.asarray(b.obsm["spatial"])[sel], np.asarray(b.obsm["sccube_spot_prop"])[sel]
        order = lambda p: np.lexsort((p[:, 1], p[:, 0]))
        i, j = order(xy10), order(saved_xy)
        assert len(xy10) == len(saved_xy) and np.allclose(xy10[i], saved_xy[j]) and np.allclose(prop10[i], saved_prop[j])
        print("target 10 matches the saved Figure 5A spots", flush=True)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig, axes = plt.subplots(1, 6, figsize=(19, 4.7))
    cells_um = native_xy * scale
    axes[0].scatter(*cells_um.T, c=[colors[t] for t in labels], s=.15, linewidths=0, rasterized=True)
    axes[0].set_title(f"Cell input\n{N:,} cells", fontweight="bold", fontsize=12)
    rows, centers = [], [(cells_um.min(0) + cells_um.max(0)) / 2]
    for ax, t in zip(axes[1:], TARGETS):
        table = tables[t]
        xy, prop, n = table[["x_um", "y_um"]].to_numpy(), table[TYPES].to_numpy(), table.n_cells.to_numpy()
        dominant, purity = np.array(TYPES)[prop.argmax(1)], prop.max(1)
        step = float(np.median(np.diff(np.unique(np.round(xy[:, 0], 6))))) if len(np.unique(xy[:, 0])) > 1 else 4100.
        centers.append((xy.min(0) + xy.max(0)) / 2)
        # Display glyph diameter scales with grid spacing, not a physical capture radius.
        ax.scatter(*xy.T, c=[colors[d] for d in dominant], s=max(1.2, (step / 4100 * 165 * .68) ** 2),
                   linewidths=0, rasterized=True)
        ax.set_title(f"Target: {t:,} cells/spot\n{len(xy):,} occupied spots", fontsize=12, fontweight="bold")
        ax.text(.5, -.08, f"Mean occupancy: {n.mean():.1f}\nMean dominant fraction: {purity.mean():.2f}",
                transform=ax.transAxes, ha="center", va="top", fontsize=10)
        rows.append(dict(target_cells_per_spot=t, n_spots=len(xy), n_cells=int(n.sum()), mean_occupancy=n.mean(),
                         median_occupancy=float(np.median(n)), min_occupancy=int(n.min()), max_occupancy=int(n.max()),
                         mean_dominant_fraction=purity.mean(), grid_step_um=step))
    for ax, c in zip(axes, centers):
        ax.set(xlim=(c[0] - 3450, c[0] + 3450), ylim=(-400, 6500), aspect="equal")
        ax.set_xticks([]); ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.plot([100, 1100], [100, 100], color="black", lw=2)
    axes[0].text(600, 190, "1 mm", ha="center", fontsize=9)
    fig.suptitle("scCube spot aggregation · Slice 5", fontsize=17, fontweight="bold", y=1.05)
    fig.legend(handles=[Patch(color=colors[t], label=f"Type {i + 1}") for i, t in enumerate(TYPES)],
               loc="lower center", bbox_to_anchor=(.5, .035), ncol=8, frameon=False, fontsize=10)
    fig.subplots_adjust(left=.01, right=.99, top=.84, bottom=.28, wspace=.12)
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(OUT / f"slice5_spot_aggregation.{ext}", dpi=500, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    summary = pd.DataFrame(rows)
    print(summary.to_string(index=False))
    if args.redraw_only:
        return
    summary.to_csv(OUT / "aggregation_summary.csv", index=False)
    (OUT / "provenance.json").write_text(json.dumps(dict(
        input=str((BASE / "sccube/cell.h5ad").resolve()), slice_id=SLICE_ID, slice_number=SLICE_ID + 1, n_cells=N,
        targets=TARGETS, platform="Visium",
        native_calls=["scCube.generate_spot_data_random(platform='Visium', gene_type='whole', n_cell=target)",
                      "scCube.utils.calculate_spot_prop"],
        validation="target 10 reproduces the saved Figure 5A slice-5 spot centres and proportions",
        our_conventions="dominant type = argmax of native proportions (ties -> lowest type); glyph size ~ grid spacing",
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        sccube_utils_sha256=hashlib.sha256((ROOT / "comparison_methods/env/scCube_src/scCube/utils.py").read_bytes()).hexdigest()),
        indent=2) + "\n")


if __name__ == "__main__":
    main()
