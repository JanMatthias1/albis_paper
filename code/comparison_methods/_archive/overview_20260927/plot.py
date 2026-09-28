"""Illustrate section stacks: rows ALBIS/scCube/Spider; columns cell/bin/spot.

Display only: project each observation to its section plane and enlarge the
distance between planes. Native coordinates in the input files are untouched.
Axis-free styling follows applications_albis/alignment/plot_figure4a.py.
"""
import argparse
import json
from pathlib import Path
import sys

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
from matplotlib.collections import PolyCollection
from matplotlib.colors import to_rgb
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from manuscript_style import CELLTYPE_COLORS


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--label", default="", help="Optional heading, e.g. Smoke test")
    parser.add_argument("--section-spacing", type=float, default=1.0,
                        help="Display plane spacing as a multiple of physical section thickness")
    parser.add_argument("--output-prefix", default="overview", help="Output filenames within --input")
    parser.add_argument("--trim-side-margins", action="store_true",
                        help="Remove the unused left/right canvas margins")
    parser.add_argument("--slice-number", type=int,
                        help="Show one slice in 2D, numbered from 1 (5 selects stored slice_id=4)")
    args = parser.parse_args()
    cfg = json.loads((args.input / "settings.json").read_text())
    if args.slice_number is not None and not 1 <= args.slice_number <= cfg["n_slices"]:
        parser.error("--slice-number must be between 1 and n_slices")
    colors = {f"type{i}": CELLTYPE_COLORS[f"Cell Type {i}"] for i in range(1, 9)}
    if not np.isfinite(args.section_spacing) or args.section_spacing < 1:
        parser.error("--section-spacing must be finite and at least 1")
    # Match the 12 × 4.3 in canvas used for the Figure 4a reference panel.
    fig = plt.figure(figsize=(4.3 if args.slice_number is not None else 12, 4.3),
                     facecolor="white")
    plt.rcParams.update({"pdf.fonttype": 42, "svg.fonttype": "none", "font.size": 9})
    methods = ["albis", "sccube", "spider"]
    row_labels = ["ALBIS", "scCube", "Spider"]
    modalities = ["cell", "bin", "spot"]
    # Only coordinates and observation metadata are needed for this display.
    # Keep the expression matrices on disk instead of loading all nine of them.
    datasets = {}
    for m in methods:
        for t in modalities:
            print(f"Loading plot metadata: {m}/{t}", flush=True)
            datasets[m, t] = ad.read_h5ad(args.input / m / f"{t}.h5ad", backed="r")
    # Native scCube staggered centers can extend beyond the source tissue.
    # Keep all outputs visible in one shared frame, without clipping or rescaling.
    lo = min(0., min(float(a.obsm["spatial_3d"].min()) for a in datasets.values()))
    hi = max(cfg["extent_um"], max(float(a.obsm["spatial_3d"].max()) for a in datasets.values()))
    # Leave a real margin around the outermost markers; otherwise 3D axes clip
    # the rim of spheres and grids at their coordinate limits.
    padding = .10 * (hi - lo)
    thickness = cfg["extent_um"] / cfg["n_slices"]
    plane_spacing = thickness * args.section_spacing
    planes = (np.arange(cfg["n_slices"]) + .5) * plane_spacing
    summaries = []
    for ri, method in enumerate(methods):
        for ci, modality in enumerate(modalities):
            a = datasets[method, modality]
            xyz = np.asarray(a.obsm["spatial_3d"]).copy()
            sections = a.obs.slice_id.to_numpy(dtype=int)
            if not np.isin(sections, np.arange(cfg["n_slices"])).all():
                raise ValueError(f"Invalid slice IDs: {method}/{modality}")
            # Exploded section view, including cells: do not imply that their
            # native continuous z coordinates were generated on flat planes.
            xyz[:, 2] = planes[sections]
            labels = a.obs.cell_type_true.astype(str).to_numpy()
            keep = labels != "unassigned"
            if "is_empty" in a.obs:
                keep &= ~a.obs.is_empty.to_numpy(dtype=bool)
            if args.slice_number is not None:
                keep &= sections == args.slice_number - 1
            xyz, labels = xyz[keep], labels[keep]
            order = np.random.default_rng(cfg["seed"]).permutation(len(xyz))
            if args.slice_number is not None:
                # Keep the original square panel size and row positions while
                # reducing the horizontal distance between panel centers.
                ax = fig.add_axes(((.5 + 1.25 * ci) / 4.3, .63 - .235 * ri, .22, .22))
                if modality == "bin":
                    # Derive each method's native grid pitch before omitting
                    # empty bins. Draw tiles in data units, not point-sized
                    # scatter markers, and keep vector boundaries in PDF/SVG.
                    grid_xy = np.asarray(a.obsm["spatial_3d"])[
                        sections == args.slice_number - 1, :2]
                    pitch = np.array([np.median(np.diff(np.unique(grid_xy[:, j])))
                                      for j in range(2)])
                    if not np.all(np.isfinite(pitch) & (pitch > 0)):
                        raise ValueError(f"Cannot determine bin grid pitch: {method}")
                    corners = np.array([[-.5, -.5], [.5, -.5], [.5, .5], [-.5, .5]])
                    vertices = xyz[:, None, :2] + corners[None, :, :] * pitch
                    # Blend with white to retain the existing .75-alpha palette
                    # without transparent edges overlapping at tile boundaries.
                    faces = [.75 * np.array(to_rgb(colors[t])) + .25 for t in labels]
                    ax.add_collection(PolyCollection(vertices, facecolors=faces,
                                      edgecolors="white", linewidths=.18,
                                      antialiaseds=True, rasterized=False))
                else:
                    ax.scatter(*xyz[order, :2].T, c=[colors[t] for t in labels[order]],
                               s=9 if modality == "spot" else 1.5, alpha=.75,
                               linewidths=0, rasterized=True)
                ax.set(xlim=(lo - padding, hi + padding), ylim=(lo - padding, hi + padding))
                ax.set_aspect("equal")
                ax.set_axis_off()
                if ri == 0:
                    ax.set_title(modality.capitalize(), fontsize=13, y=1.08)
                summaries.append({"method": method, "modality": modality,
                                  "slice_number": args.slice_number, "slice_id": args.slice_number - 1,
                                  "n_observations": int((sections == args.slice_number - 1).sum()),
                                  "n_displayed": len(xyz), "n_genes": a.n_vars})
                continue
            # Fixed placements avoid Matplotlib shrinking each 3D view into the
            # centre of a wide subplot-grid cell.
            ax = fig.add_axes((.19 + .18 * ci, .63 - .235 * ri, .32, .22), projection="3d")
            ax.scatter(*xyz[order].T, c=[colors[t] for t in labels[order]],
                       s={"bin": 1.5, "spot": 9, "cell": 1.5}[modality],
                       marker="s" if modality == "bin" else "o", alpha=.75,
                       # Axes3D makes a square clipping patch inside the wide
                       # panel. The zoomed projection extends beyond that patch;
                       # allow the full silhouette into the surrounding margin.
                       linewidths=0, depthshade=False, rasterized=True, clip_on=False)
            ax.set_proj_type("ortho")
            ax.view_init(elev=12, azim=-60)
            # Shared x/y frame and display section positions for all panels.
            ax.set(xlim=(lo - padding, hi + padding), ylim=(lo - padding, hi + padding),
                   zlim=(0, cfg["n_slices"] * plane_spacing))
            # A mildly compressed display z scale lets the three native section
            # stacks use the short Figure 4a-style landscape panel effectively.
            ax.set_box_aspect((hi - lo + 2 * padding, hi - lo + 2 * padding,
                               .65 * cfg["n_slices"] * plane_spacing), zoom=1.9)
            ax.set_axis_off()
            if ri == 0:
                ax.set_title(modality.capitalize(), fontsize=13, y=1.08)
            summaries.append({"method": method, "modality": modality, "n_observations": a.n_obs,
                              "n_displayed": len(xyz), "n_genes": a.n_vars})
    # Keep method labels outside the 3D clipping region and away from the data.
    for y, label in zip((.74, .505, .27), row_labels):
        fig.text(.045 if args.slice_number is not None else .275, y, label,
                 rotation=90, ha="center", va="center", fontsize=12)
    handles = [plt.Line2D([], [], markerfacecolor=c, markeredgecolor="none", marker="o", ls="", label=f"type{i}")
               for i, c in enumerate(colors.values(), 1)]
    # Matplotlib fills legend columns first; reorder for type1–4 above type5–8.
    handles = [handles[i] for i in (0, 4, 1, 5, 2, 6, 3, 7)]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.52, .005), ncol=4,
               frameon=False, title="Cell Type", borderaxespad=.2)
    if args.label:
        fig.suptitle(args.label, fontsize=10, y=.99)
    for ext in ("png", "pdf", "svg"):
        print(f"Saving {args.output_prefix}.{ext}", flush=True)
        # Explicit canvas bounds include the unclipped 3D scatter silhouettes;
        # automatic tight bounds include the much wider, invisible axes boxes.
        bounds = (Bbox.from_extents(3.12, 0, 9.36, 4.3)
                  if args.trim_side_margins and args.slice_number is None else None)
        fig.savefig(args.input / f"{args.output_prefix}.{ext}", dpi=300,
                    facecolor="white", bbox_inches=bounds)
    plt.close(fig)
    for a in datasets.values():
        a.file.close()
    summary_name = "plot_summary.json" if args.output_prefix == "overview" else f"{args.output_prefix}_summary.json"
    (args.input / summary_name).write_text(json.dumps(summaries, indent=2) + "\n")


if __name__ == "__main__":
    main()
