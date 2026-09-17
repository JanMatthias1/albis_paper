"""
Plot the intact simulated 3D sphere before slicing, capture, or batch effects.

This script uses the same manuscript-scale simulation parameters as
sim_paper/code/data/generate_simulation.py, but stops at the new base-sphere
checkpoint and plots cell centroids colored by ground-truth domain.
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
PROJECT_ROOT = SCRIPT_DIR.parents[2]
REPO_ROOT = find_repo_root(SCRIPT_DIR)
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

sys.modules.pop("albis", None)
import albis as ab


DOMAIN_COLORS = {
    "D0": "#7B2CBF",
    "D1": "#1E88E5",
    "D2": "#00A6A6",
    "D3": "#66BB6A",
    "D4": "#F6A21A",
    "D5": "#D7263D",
}


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
        output_modalities=("cell",),
    )


def parse_args():
    parser = argparse.ArgumentParser(description="Plot the intact simulated sphere.")
    parser.add_argument("--max-cells", type=int, default=180_000)
    parser.add_argument("--outdir", type=Path, default=SCRIPT_DIR / "outputs")
    parser.add_argument("--dpi", type=int, default=600)
    return parser.parse_args()


def style_3d_axis(ax, radius):
    ax.set_axis_off()
    ax.set_xlim(-radius, radius)
    ax.set_ylim(-radius, radius)
    ax.set_zlim(-radius, radius)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=21, azim=-47)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.fill = False
        axis.pane.set_edgecolor("none")


def main():
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    print("Generating intact base sphere...")
    base = ab.simulate_3d_molecule_sphere_base(**manuscript_simulation_kwargs())
    adata = base["adata_cell_true"]
    coords = np.asarray(adata.obsm["spatial"])
    domains = adata.obs["domain_true"].astype(str).to_numpy()

    rng = np.random.default_rng(2025)
    if coords.shape[0] > args.max_cells:
        idx = np.sort(rng.choice(coords.shape[0], size=args.max_cells, replace=False))
        coords = coords[idx]
        domains = domains[idx]

    draw_order = rng.permutation(coords.shape[0])
    colors = np.array([DOMAIN_COLORS.get(domain, "#8A8F98") for domain in domains])

    width_in = 92 / 25.4
    fig = plt.figure(figsize=(width_in, width_in), constrained_layout=False)
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(
        coords[draw_order, 0],
        coords[draw_order, 1],
        coords[draw_order, 2],
        c=colors[draw_order],
        s=0.18,
        alpha=0.82,
        linewidths=0,
        depthshade=False,
        rasterized=True,
    )
    style_3d_axis(ax, radius=6000)
    ax.set_title("Intact 3D Tissue Sphere", fontsize=8, fontweight="bold", y=0.93)

    handles = [
        plt.Line2D([0], [0], marker="o", color="none", label=domain, markerfacecolor=color, markersize=4)
        for domain, color in DOMAIN_COLORS.items()
    ]
    ax.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.07),
        ncol=6,
        frameon=False,
        fontsize=5,
        handletextpad=0.2,
        columnspacing=0.7,
    )

    for suffix in ("png", "pdf"):
        out = args.outdir / f"figure_1a_intact_sphere.{suffix}"
        fig.savefig(out, dpi=args.dpi, bbox_inches="tight", pad_inches=0.02)
        print(f"Saved: {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
