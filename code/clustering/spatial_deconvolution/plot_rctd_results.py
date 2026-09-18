#!/usr/bin/env python
"""
Figure 3 -- RCTD spot deconvolution, presentation plots.

Reads run_rctd_spot.R/.sh's output (data/figure_3/spatial_deconvolution/RCTD/spot/)
and produces two figures:

  rctd_spatial_zoom.png   -- one representative slice (slice_id=5, the
      project's usual "representative near-equatorial interior section"),
      cropped to a sub-region so individual spots are visible as actual
      spot-radius circles (not a dense blur). Rows = cell type, columns =
      [true fraction, RCTD estimated fraction]. Each row is colored on its
      own 0-1 sequential ramp built from that cell type's own identity hue
      (CELLTYPE_COLORS below), matching the tab20-derived palette already
      used for cell_type_true across every Figure 3 UMAP (verified directly
      off data/figure_3/pca_harmony_single_cell/cell/leiden_pca_qc_celltype_matched/
      plots/umap_true_vs_predicted_cell_type_true.png's legend swatches,
      2026-09-18) -- so a reader who already knows "type3 = light green"
      from those UMAPs recognizes the same identity here, rather than
      learning a second, unrelated color code for the same 8 categories.
      Since hue now differs per row, each row gets its own compact
      colorbar (a single shared bar can't represent 8 different hues).

  rctd_error_boxplot.png -- per-cell-type quantification: box plot of
      per-spot absolute error |estimated - true| (the standard
      distributional error metric), one box per cell type, with that
      type's overall Pearson r (from per_celltype_metrics.csv, already
      computed by run_rctd_spot.R) annotated above each box. RMSE is the
      root-mean-square of the same per-spot errors the box shows; PCC is
      an aggregate across all spots and has no natural per-spot
      distribution, hence the annotation rather than a second box.

Color: sequential blue (one hue, light->dark) for magnitude (fraction),
per the project's dataviz convention -- never a rainbow for a continuous
quantity. The 8 error boxes share one neutral color (cell type is already
encoded by x-position + label; coloring each box a different hue would add
an 8-way categorical rainbow the palette's own all-pairs cap doesn't
support and the data doesn't need).
"""

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize, LinearSegmentedColormap, to_rgb
import numpy as np
import pandas as pd
import anndata as ad

SPOT_H5AD = "/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_2/smaller_sphere/data/packing_pf0p04_log_mu_-2.0_theta_0.25_jitter0.10_bsigma03/simulation_spot_z_qc.h5ad"
OUT_DIR = "/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_3/spatial_deconvolution/RCTD/spot"
CELL_TYPES = [f"type{i}" for i in range(1, 9)]

# This project's established cell_type_true identity colors -- tab20 applied
# to the 8 alphabetically-sorted category codes, exactly as
# code/clustering/01.2_pca_harmony.py's color_values()/plot_umap_before_after()
# (and step03_cluster_and_plot.py's plot_umap_true_vs_predicted) already do
# for every Figure 3 cell-type UMAP. Confirmed against the actual rendered
# legend, not just re-derived from the tab20 formula.
CELLTYPE_COLORS = {
    "type1": "#1f77b4",  # blue
    "type2": "#ff7f0e",  # orange
    "type3": "#98df8a",  # light green
    "type4": "#9467bd",  # purple
    "type5": "#c49c94",  # tan
    "type6": "#7f7f7f",  # gray
    "type7": "#dbdb8d",  # olive
    "type8": "#9edae5",  # light cyan
}


def sequential_cmap_for(hex_color, light_frac=0.85, dark_frac=0.35):
    """One-hue sequential ramp (light tint -> base hue -> darkened base),
    same light->dark structure as the generic SEQ_CMAP below, anchored on
    a specific cell type's own identity color instead of blue."""
    base = np.array(to_rgb(hex_color))
    light = np.array([1.0, 1.0, 1.0]) * light_frac + base * (1 - light_frac)
    dark = base * (1 - dark_frac)
    return LinearSegmentedColormap.from_list(f"seq_{hex_color}", [light, base, dark])


# dataviz skill palette: sequential blue ramp (magnitude), one hue light->dark
# -- still used for the box plot below, which isn't per-cell-type-colored.
SEQ_CMAP = matplotlib.colors.LinearSegmentedColormap.from_list(
    "seq_blue", ["#cde2fb", "#6da7ec", "#2a78d6", "#0d366b"]
)
BOX_COLOR = "#6da7ec"   # sequential step 300 -- neutral, not a per-type hue
BOX_EDGE = "#184f95"    # sequential step 600
ANNOT_COLOR = "#eb6834"  # categorical slot 2 (orange) -- visually distinct "this is a different metric" cue
GRID_COLOR = "#e1e0d9"
MUTED = "#898781"

SLICE_ID = "5"
ZOOM_XMIN, ZOOM_XMAX = 0, 1200
ZOOM_YMIN, ZOOM_YMAX = 0, 1200


def load_data():
    spot = ad.read_h5ad(SPOT_H5AD)
    true_frac = pd.DataFrame(
        np.asarray(spot.obsm["cell_type_frac_true"]), columns=CELL_TYPES, index=spot.obs_names
    )
    spatial = pd.DataFrame(
        np.asarray(spot.obsm["spatial"])[:, :2], columns=["x", "y"], index=spot.obs_names
    )
    meta = spot.obs[["slice_id"]].copy()
    spot_radius_um = float(spot.obs["spot_radius_um"].iloc[0])

    est_wide = pd.read_csv(os.path.join(OUT_DIR, "estimated_fractions_wide.csv"))
    # spot_id was assigned "spot_<1-based positional index>" in run_rctd_spot.R,
    # in the same row order the h5ad was read in -- map back to obs_names.
    positional_idx = est_wide["spot_id"].str.replace("spot_", "", regex=False).astype(int) - 1
    est_wide.index = spot.obs_names[positional_idx.to_numpy()]
    est_frac = est_wide[CELL_TYPES]

    common = est_frac.index
    return (
        true_frac.loc[common],
        est_frac.loc[common],
        spatial.loc[common],
        meta.loc[common],
        spot_radius_um,
    )


def plot_spatial_zoom(true_frac, est_frac, spatial, meta, spot_radius_um,
                       cell_types=CELL_TYPES, pairs_per_row=1, out_name="rctd_spatial_zoom.png"):
    """cell_types: which types to include (in order). pairs_per_row: how many
    [True, Estimated] column-pairs sit side by side per row -- 1 gives the
    original one-type-per-row strip, 2 gives a 4-row x 4-col grid for 8
    types (2 types per row), etc. Each type keeps its own identity color +
    its own colorbar regardless of layout."""
    m = (meta["slice_id"].astype(str) == SLICE_ID).to_numpy()
    m &= (spatial["x"].between(ZOOM_XMIN, ZOOM_XMAX) & spatial["y"].between(ZOOM_YMIN, ZOOM_YMAX)).to_numpy()
    xy = spatial.loc[m]
    n = m.sum()
    print(f"[spatial zoom -> {out_name}] slice {SLICE_ID}, cropped to x/y in "
          f"[{ZOOM_XMIN},{ZOOM_XMAX}]x[{ZOOM_YMIN},{ZOOM_YMAX}]: {n} spots, "
          f"{len(cell_types)} cell types, {pairs_per_row} pair(s)/row")

    # Draw circles a few times their true radius -- the real footprint
    # (27.5um) is tiny relative to spot spacing (~100um), so at true size
    # the spots would read as barely-visible dots. A fixed visual multiplier
    # is standard practice in these figures; the crop/zoom already conveys
    # true spacing, this only exaggerates each spot's own disc.
    draw_radius = spot_radius_um * 1.9
    n_rows = -(-len(cell_types) // pairs_per_row)  # ceil
    n_cols = pairs_per_row * 2

    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(2.5 * n_cols, 2.05 * n_rows),
        constrained_layout=True, squeeze=False,
    )
    norm = Normalize(vmin=0, vmax=1)

    for idx, ct in enumerate(cell_types):
        row, pair = divmod(idx, pairs_per_row)
        col_true, col_est = pair * 2, pair * 2 + 1
        row_cmap = sequential_cmap_for(CELLTYPE_COLORS[ct])
        for col, label, frac_df in [(col_true, "True", true_frac), (col_est, "RCTD estimated", est_frac)]:
            ax = axes[row, col]
            vals = frac_df.loc[xy.index, ct].to_numpy()
            circles = [Circle((x, y), draw_radius) for x, y in zip(xy["x"], xy["y"])]
            pc = PatchCollection(circles, array=vals, cmap=row_cmap, norm=norm,
                                  edgecolor="white", linewidth=0.3)
            ax.add_collection(pc)
            ax.set_xlim(ZOOM_XMIN, ZOOM_XMAX)
            ax.set_ylim(ZOOM_YMIN, ZOOM_YMAX)
            ax.set_aspect("equal")
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_color(GRID_COLOR)
            if row == 0:
                ax.set_title(label, fontsize=11, color="#0b0b0b")
            if col == col_true:
                ct_label = f"Cell Type\n{ct.replace('type', '')}"
                ax.text(-0.18, 0.5, ct_label, transform=ax.transAxes, fontsize=11,
                        color=CELLTYPE_COLORS[ct], fontweight="bold",
                        ha="right", va="center", linespacing=1.3)

        # Per-type colorbar -- hue differs by type, so one shared bar can't
        # represent every identity color in the figure.
        sm = ScalarMappable(norm=norm, cmap=row_cmap)
        fig.colorbar(sm, ax=axes[row, col_est], fraction=0.12, pad=0.03, aspect=8)

    fig.suptitle(f"RCTD deconvolution, slice {SLICE_ID} ({n} spots)", fontsize=11, fontweight="bold")
    out = os.path.join(OUT_DIR, out_name)
    fig.savefig(out, dpi=180)
    plt.close(fig)
    print("wrote", out)


def plot_error_boxplot(true_frac, est_frac):
    """Violin (distribution shape) + narrow box (median/IQR) + jittered
    per-spot points, one per cell type, colored with that type's own
    identity color (CELLTYPE_COLORS) -- same palette as the spatial zoom
    figures, so this plot reads as the same 8 categories rather than an
    unrelated color code. r (Pearson correlation, an aggregate that has no
    per-spot distribution) is still annotated above each violin."""
    per_type = pd.read_csv(os.path.join(OUT_DIR, "per_celltype_metrics.csv"))
    pcc = dict(zip(per_type["cell_type"], per_type["pearson_r"]))

    abs_err = (est_frac[CELL_TYPES] - true_frac[CELL_TYPES]).abs()
    data = [abs_err[ct].to_numpy() for ct in CELL_TYPES]
    positions = np.arange(len(CELL_TYPES))

    # Violin shape is a KDE spanning [min(data), max(data)] by default -- a
    # a handful of true outlier spots per type (up to ~0.25) stretch that
    # into a long thin spike reaching almost to the top of the axis, right
    # where the r-annotation sits. Cap the values FED TO THE VIOLIN ONLY
    # (not the box/whiskers/median, which still reflect the real data) at a
    # shared 99th-percentile cutoff across all 8 types, so every violin's
    # visual extent is compressed to the same, much shorter, typical range.
    pctile_cap = float(np.percentile(np.concatenate(data), 99))
    data_for_violin = [np.clip(vals, None, pctile_cap) for vals in data]

    fig, ax = plt.subplots(figsize=(10, 5.5))

    vp = ax.violinplot(data_for_violin, positions=positions, widths=0.8,
                        showmedians=False, showextrema=False)
    for i, body in enumerate(vp["bodies"]):
        ct = CELL_TYPES[i]
        body.set_facecolor(CELLTYPE_COLORS[ct])
        body.set_edgecolor(CELLTYPE_COLORS[ct])
        body.set_alpha(0.45)
        body.set_zorder(3)

    bp = ax.boxplot(
        data, positions=positions, widths=0.12,
        patch_artist=True, showfliers=False, zorder=4,
        medianprops=dict(color="#0b0b0b", linewidth=1.6),
        boxprops=dict(facecolor="white", edgecolor="#0b0b0b", linewidth=1.0),
        whiskerprops=dict(color="#0b0b0b", linewidth=1.0),
        capprops=dict(color="#0b0b0b", linewidth=1.0),
    )
    # Headroom now scales off the (compressed) violin cap, not the raw data
    # max, so the annotation row sits in clean whitespace above every violin.
    annot_y = pctile_cap * 1.1
    for i, ct in enumerate(CELL_TYPES):
        ax.annotate(f"r = {pcc[ct]:.2f}", xy=(i, annot_y), ha="center",
                    fontsize=9, color="#0b0b0b", fontweight="bold")

    ax.set_xticks(positions)
    ax.set_xticklabels(
        [f"Cell Type\n{ct.replace('type', '')}" for ct in CELL_TYPES]
    )
    for tick, ct in zip(ax.get_xticklabels(), CELL_TYPES):
        tick.set_color(CELLTYPE_COLORS[ct])
        tick.set_fontweight("bold")
    ax.set_ylabel("per-spot absolute error  |estimated - true fraction|")
    ax.set_title("RCTD deconvolution accuracy by cell type", fontweight="bold")
    ax.set_ylim(0, pctile_cap * 1.25)
    ax.grid(axis="y", color=GRID_COLOR, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(MUTED)

    fig.tight_layout()
    out = os.path.join(OUT_DIR, "rctd_error_boxplot.png")
    fig.savefig(out, dpi=180)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    true_frac, est_frac, spatial, meta, spot_radius_um = load_data()
    print(f"[load] {len(true_frac)} spots with an RCTD estimate")

    # Full panel, all 8 types, one-type-per-row strip.
    plot_spatial_zoom(true_frac, est_frac, spatial, meta, spot_radius_um,
                       cell_types=CELL_TYPES, pairs_per_row=1, out_name="rctd_spatial_zoom.png")
    # Full panel, all 8 types, 4x4 grid (2 types per row) -- more compact.
    plot_spatial_zoom(true_frac, est_frac, spatial, meta, spot_radius_um,
                       cell_types=CELL_TYPES, pairs_per_row=2, out_name="rctd_spatial_zoom_4x4.png")
    # Main-figure version: cell types 1 and 4 only.
    plot_spatial_zoom(true_frac, est_frac, spatial, meta, spot_radius_um,
                       cell_types=["type1", "type4"], pairs_per_row=1, out_name="rctd_spatial_zoom_main.png")
    # 4-type subset, same 2-types-per-row grid style as the 8-type 4x4 (comes
    # out 2 rows x 4 cols for 4 types).
    plot_spatial_zoom(true_frac, est_frac, spatial, meta, spot_radius_um,
                       cell_types=["type1", "type2", "type3", "type4"], pairs_per_row=2,
                       out_name="rctd_spatial_zoom_4types.png")

    plot_error_boxplot(true_frac, est_frac)
