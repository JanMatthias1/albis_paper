"""
Plot capture-window selection and actual Visium HD-like binned sections.

Figure 1B: cell centroids in the intact sphere, with cells outside the
planned capture window shown in grey and captured cells colored by domain.

Figure 1C: the same base sphere after sectioning into Z slices and aggregating
transcripts into Visium HD-like bins.

Figure 1D: the sectioned bins again, but using each slice's unaligned
(randomly rotated and translated) coordinates, showing the naive stack
before any registration/coordinate-shift correction is applied.
"""

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def find_repo_root(start):
    for candidate in [start, *start.parents]:
        nested = candidate / "albis"
        if (nested / "pyproject.toml").is_file() and (nested / "albis" / "__init__.py").is_file():
            return nested
        if (candidate / "pyproject.toml").is_file() and (candidate / "albis" / "__init__.py").is_file():
            return candidate
    raise RuntimeError("Could not find the albis repository root.")


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = find_repo_root(SCRIPT_DIR)
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

sys.modules.pop("albis", None)
import albis as ab


DOMAIN_COLORS = {
    "D0": "#8E63C7",
    "D1": "#4C8FD5",
    "D2": "#43B7A5",
    "D3": "#72C69C",
    "D4": "#F2A65A",
    "D5": "#E65F5C",
}
OUTSIDE_COLOR = "#C9CDD3"


def manuscript_simulation_kwargs():
    return dict(
        sphere_R_um=6000.0,
        capture_window_um=(6500.0, 6500.0),
        n_domains=6,
        core_frac=0.55,
        core_bump_amp=0.25,
        wedge_angle_amp_deg=25.0,
        noise_terms=16,
        noise_freq_range=(3.0, 6.0),
        boundary_fuzz_width_deg=6.0,
        boundary_fuzz_flip_prob=0.15,
        core_fuzz_width_um=300.0,
        core_fuzz_flip_prob=0.25,
        n_cells=600_000,
        cell_radius_kwargs=dict(
            radius_dist="lognormal",
            r_mean=7.5,
            r_sigma=0.28,
            r_min=4.0,
            r_max=14.0,
        ),
        n_cell_types=8,
        domain_type_mix=np.array(
            [
                [0.18, 0.18, 0.13, 0.12, 0.11, 0.10, 0.09, 0.09],
                [0.11, 0.12, 0.18, 0.18, 0.13, 0.10, 0.09, 0.09],
                [0.10, 0.11, 0.12, 0.13, 0.18, 0.18, 0.09, 0.09],
                [0.10, 0.10, 0.11, 0.12, 0.13, 0.13, 0.16, 0.15],
                [0.14, 0.13, 0.12, 0.11, 0.12, 0.13, 0.13, 0.12],
                [0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125],
            ],
            dtype=float,
        ),
        seed=2025,
        output_modalities=("bin",),
    )


def parse_args():
    parser = argparse.ArgumentParser(description="Plot capture-window and binned-section panels.")
    parser.add_argument("--max-cells", type=int, default=180_000)
    parser.add_argument("--max-bins", type=int, default=220_000)
    parser.add_argument("--outdir", type=Path, default=SCRIPT_DIR / "outputs")
    parser.add_argument("--dpi", type=int, default=600)
    parser.add_argument("--max-shift", type=float, default=1500.0,
                        help="Maximum per-slice translation in µm (use 200 for the original schematic).")
    return parser.parse_args()


def style_3d_axis(ax, radius=6000, zlim=None):
    ax.set_axis_off()
    ax.set_xlim(-radius, radius)
    ax.set_ylim(-radius, radius)
    ax.set_zlim(*(zlim if zlim is not None else (-radius, radius)))
    ax.set_box_aspect((1, 1, 0.9))
    ax.view_init(elev=21, azim=-47)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.fill = False
        axis.pane.set_edgecolor("none")


def save_figure(fig, outdir, stem, dpi):
    for suffix in ("png", "pdf"):
        out = outdir / f"{stem}.{suffix}"
        fig.savefig(out, dpi=dpi, bbox_inches="tight", pad_inches=0.02)
        print(f"Saved: {out}")
    plt.close(fig)


def plot_capture_window(base, outdir, dpi, max_cells):
    adata = base["adata_cell_true"]
    coords = np.asarray(adata.obsm["spatial"])
    domains = adata.obs["domain_true"].astype(str).to_numpy()

    capture_size = np.asarray(base["meta"]["captures"]["capture_window_um"], dtype=float)
    capture_center = np.asarray(base["meta"]["captures"]["capture_window_center_um"], dtype=float)
    half_size = capture_size / 2.0
    in_capture = (
        (coords[:, 0] >= capture_center[0] - half_size[0])
        & (coords[:, 0] <= capture_center[0] + half_size[0])
        & (coords[:, 1] >= capture_center[1] - half_size[1])
        & (coords[:, 1] <= capture_center[1] + half_size[1])
    )

    rng = np.random.default_rng(2026)
    if coords.shape[0] > max_cells:
        idx = np.sort(rng.choice(coords.shape[0], size=max_cells, replace=False))
        coords = coords[idx]
        domains = domains[idx]
        in_capture = in_capture[idx]

    inside_colors = np.array([DOMAIN_COLORS.get(domain, "#8A8F98") for domain in domains[in_capture]])

    width_in = 92 / 25.4
    fig = plt.figure(figsize=(width_in, width_in), constrained_layout=False)
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(
        coords[~in_capture, 0],
        coords[~in_capture, 1],
        coords[~in_capture, 2],
        c=OUTSIDE_COLOR,
        s=0.18,
        alpha=0.13,
        linewidths=0,
        depthshade=False,
        rasterized=True,
    )
    ax.scatter(
        coords[in_capture, 0],
        coords[in_capture, 1],
        coords[in_capture, 2],
        c=inside_colors,
        s=0.2,
        alpha=0.86,
        linewidths=0,
        depthshade=False,
        rasterized=True,
    )

    x0, x1 = capture_center[0] - half_size[0], capture_center[0] + half_size[0]
    y0, y1 = capture_center[1] - half_size[1], capture_center[1] + half_size[1]
    z = -6100
    ax.plot([x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0], [z] * 5, color="#222222", lw=0.8)

    style_3d_axis(ax)
    ax.set_title("Planned Visium HD capture", fontsize=8, fontweight="bold", y=0.93)
    save_figure(fig, outdir, "figure_1b_capture_window", dpi)


def run_sectioning(base, max_shift=1500.0):
    # Larger translations make slice displacement visible in the schematic.
    # Retain the original 180-degree rotation range.
    return ab.section_3d_molecule_sphere(
        base,
        n_slices=10,
        batch_sigma=0.22,
        capture_window_um=(6500.0, 6500.0),
        capture_window_center_um=(0.0, 0.0),
        bin_size_um=8.0,
        output_modalities=("bin",),
        slice_axes=("Z",),
        max_deg=180.0,
        max_shift=max_shift,
    )


def _plot_stacked_bins(adata, outdir, dpi, max_bins, spatial_key, title, stem):
    coords = np.asarray(adata.obsm[spatial_key])
    domains = adata.obs["domain_true"].astype(str).to_numpy()
    slice_ids = adata.obs["slice_id"].astype(int).to_numpy()

    rng = np.random.default_rng(2027)
    nonzero = np.asarray(adata.X.sum(axis=1)).ravel() > 0
    keep_idx = np.flatnonzero(nonzero)
    if keep_idx.size > max_bins:
        keep_idx = np.sort(rng.choice(keep_idx, size=max_bins, replace=False))

    coords = coords[keep_idx]
    domains = domains[keep_idx]
    slice_ids = slice_ids[keep_idx]
    colors = np.array([DOMAIN_COLORS.get(domain, "#8A8F98") for domain in domains])

    z_offset = 14000.0
    z_stack = slice_ids.astype(float) * z_offset
    order = rng.permutation(coords.shape[0])

    width_in = 92 / 25.4
    fig = plt.figure(figsize=(width_in, width_in * 1.05), constrained_layout=False)
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(
        coords[order, 0],
        coords[order, 1],
        z_stack[order],
        c=colors[order],
        s=0.08,
        alpha=0.74,
        linewidths=0,
        depthshade=False,
        rasterized=True,
    )
    # Keep displaced slices within view if a larger shift is requested.
    radius = max(6000.0, float(np.max(np.abs(coords[:, :2]))) * 1.05)
    style_3d_axis(ax, radius=radius, zlim=(-0.5 * z_offset, 9.5 * z_offset))
    ax.set_box_aspect((1, 1, 1.25))
    ax.view_init(elev=22, azim=-60)
    ax.set_title(title, fontsize=8, fontweight="bold", y=0.94)
    save_figure(fig, outdir, stem, dpi)


def plot_binned_sections(adata, outdir, dpi, max_bins):
    _plot_stacked_bins(
        adata,
        outdir,
        dpi,
        max_bins,
        spatial_key="spatial_3d",
        title="Sectioned Visium HD bins",
        stem="figure_1c_sectioned_bins",
    )


def plot_unaligned_sections(adata, outdir, dpi, max_bins):
    _plot_stacked_bins(
        adata,
        outdir,
        dpi,
        max_bins,
        spatial_key="spatial_3d_unaligned",
        title="Reconstructed bins, coordinate-shifted",
        stem="figure_1d_unaligned_bins",
    )


def main():
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    print("Generating base sphere with full molecule stream for binning...")
    base = ab.simulate_3d_molecule_sphere_base(**manuscript_simulation_kwargs())
    plot_capture_window(base, args.outdir, args.dpi, args.max_cells)

    print("Sectioning sphere into Visium HD-like bins...")
    sim = run_sectioning(base, max_shift=args.max_shift)
    adata = sim["bin_adatas"]["Z"]
    plot_binned_sections(adata, args.outdir, args.dpi, args.max_bins)
    plot_unaligned_sections(adata, args.outdir, args.dpi, args.max_bins)


if __name__ == "__main__":
    main()
