"""
Generate a sim_app synthetic spatial-transcriptomics dataset and save it,
along with diagnostic plots, under sim_paper/data.

Run from anywhere inside the sim_project tree, e.g.:
    cd /dcs04/hicks/data/Jan/sim_project/sim_paper
    python generate_simulation.py
"""

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # headless-safe backend, must be set before pyplot/sim_app plotting is used
import matplotlib.pyplot as plt
import numpy as np


# ---------------------------------------------------------------------------
# Locate and import sim_app (same logic as the tutorial notebook)
# ---------------------------------------------------------------------------
def find_repo_root(start):
    for candidate in [start, *start.parents]:
        if (candidate / "pyproject.toml").is_file() and (candidate / "sim_app" / "__init__.py").is_file():
            return candidate
        nested = candidate / "sim_app"
        if (nested / "pyproject.toml").is_file() and (nested / "sim_app" / "__init__.py").is_file():
            return nested
    raise RuntimeError("Could not find the sim_app repository root.")


repo_root = find_repo_root(Path.cwd().resolve())
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))

sys.modules.pop("sim_app", None)
import sim_app
from sim_app.simulation_sphere import _to_serializable, simulate_3d_molecule_sphere_multires

print("Python:", sys.executable)
print("sim_app module:", getattr(sim_app, "__file__", "<no __file__>"))
print("sim_app version:", getattr(sim_app, "__version__", "<no __version__>"))


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
DATA_DIR = Path("/dcs04/hicks/data/Jan/sim_project/sim_paper/data")
PLOTS_DIR = DATA_DIR / "plots"
DATA_DIR.mkdir(parents=True, exist_ok=True)
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

H5AD_PATH = DATA_DIR / "simulation_bin_z.h5ad"
SUMMARY_PATH = DATA_DIR / "simulation_bin_z_summary.json"


# ---------------------------------------------------------------------------
# Generate data
# ---------------------------------------------------------------------------
# The previous 200 um sphere could not physically fit 20,000 non-overlapping
# cells with the requested cell radii. This larger tissue matches the intended
# 12 mm-diameter simulation and avoids that placement failure.
cell_radius_kwargs = dict(
    radius_dist="lognormal",
    r_mean=7.5,
    r_sigma=0.28,
    r_min=4.0,
    r_max=14.0,
)

domain_type_mix = np.array(
    [
        [0.18, 0.18, 0.13, 0.12, 0.11, 0.10, 0.09, 0.09],
        [0.11, 0.12, 0.18, 0.18, 0.13, 0.10, 0.09, 0.09],
        [0.10, 0.11, 0.12, 0.13, 0.18, 0.18, 0.09, 0.09],
        [0.10, 0.10, 0.11, 0.12, 0.13, 0.13, 0.16, 0.15],
        [0.14, 0.13, 0.12, 0.11, 0.12, 0.13, 0.13, 0.12],
        [0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125],
    ],
    dtype=float,
)

SIM_PARAMS = dict(
    # tissue: 12 mm diameter
    sphere_R_um=6000.0,

    # capture windows
    xenium_capture_size_um=(12000.0, 24000.0),
    visium_capture_size_um=(6500.0, 6500.0),
    visium_capture_center_um=(0.0, 0.0),

    # domains: noisier boundaries
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

    # cells
    n_cells=600_000,
    cell_radius_kwargs=cell_radius_kwargs,
    allow_cell_overlap=False,

    # cell types / genes
    n_cell_types=8,
    domain_type_mix=domain_type_mix,
    marker_genes_per_type=80,
    noise_gene_frac=0.10,
    shared_marker_frac=0.25,

    # slices + batch effect
    n_slices=10,
    batch_sigma=0.22,

    # VisiumHD bins / Visium spots
    bin_size_um=8.0,
    spot_spacing_um=100.0,
    spot_radius_um=27.5,

    # unaligned perturbation
    max_deg=270.0,
    max_shift=3000.0,
    base_seed_unaligned=12345,

    # Generate only the output needed by the clustering pipeline.
    output_modalities=("bin",),
    slice_axes=("Z",),

    seed=2025,
)

print("\nGenerating data with parameters:")
for k, v in SIM_PARAMS.items():
    if k == "domain_type_mix":
        print(f"  {k}: shape={v.shape}")
    else:
        print(f"  {k}: {v}")

sim = simulate_3d_molecule_sphere_multires(**SIM_PARAMS)
adata = sim["bin_adatas"]["Z"]
adata.uns["sim_params"] = _to_serializable(sim["meta"])
adata.uns["output"] = {
    "platform": "bin",
    "slice_axis": "Z",
}
print("\nGenerated AnnData:")
print(adata)


# ---------------------------------------------------------------------------
# Inspect / summarize
# ---------------------------------------------------------------------------
summary = sim_app.describe(adata)
print("\nSummary:")
print(json.dumps(summary, indent=2, default=str))

with open(SUMMARY_PATH, "w") as f:
    json.dump(summary, f, indent=2, default=str)
print(f"\nSaved summary -> {SUMMARY_PATH}")


# ---------------------------------------------------------------------------
# Diagnostic plots
# ---------------------------------------------------------------------------
def savefig(fig, name):
    path = PLOTS_DIR / name
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {path}")


print("\nGenerating diagnostic plots...")

# 1. Aligned spatial coordinates colored by ground-truth domain
fig = sim_app.plot(adata, view="2d", coordinates="aligned", color="domain_true", point_size=4)
savefig(fig, "01_aligned_domain_true.png")

# 2. Aligned spatial coordinates colored by ground-truth cell type
fig = sim_app.plot(adata, view="2d", coordinates="aligned", color="cell_type_true", point_size=4)
savefig(fig, "02_aligned_cell_type_true.png")

# 3. Unaligned (per-slice, pre-registration) coordinates colored by slice id
fig = sim_app.plot(adata, view="2d", coordinates="unaligned", color="slice_id", point_size=4)
savefig(fig, "03_unaligned_slice_id.png")

# 4. QC: total counts / genes detected per bin
X = adata.X
if hasattr(X, "toarray"):
    total_counts = np.asarray(X.sum(axis=1)).ravel()
    genes_detected = np.asarray((X > 0).sum(axis=1)).ravel()
else:
    total_counts = np.asarray(X.sum(axis=1)).ravel()
    genes_detected = np.asarray((X > 0).sum(axis=1)).ravel()

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
axes[0].hist(total_counts, bins=50, color="steelblue")
axes[0].set_xlabel("Total counts per bin")
axes[0].set_ylabel("Number of bins")
axes[0].set_title("Total counts distribution")

axes[1].hist(genes_detected, bins=50, color="darkorange")
axes[1].set_xlabel("Genes detected per bin")
axes[1].set_ylabel("Number of bins")
axes[1].set_title("Genes detected distribution")

fig.tight_layout()
savefig(fig, "04_qc_count_distributions.png")

# 5. Number of bins per slice
slice_counts = adata.obs["slice_id"].value_counts().sort_index()
fig, ax = plt.subplots(figsize=(6, 4))
ax.bar(slice_counts.index.astype(str), slice_counts.values, color="slategray")
ax.set_xlabel("Slice ID")
ax.set_ylabel("Number of bins")
ax.set_title("Bins per slice")
fig.tight_layout()
savefig(fig, "05_bins_per_slice.png")

# 6. Domain composition per slice (stacked bar)
comp = adata.obs.groupby(["slice_id", "domain_true"]).size().unstack(fill_value=0)
comp_frac = comp.div(comp.sum(axis=1), axis=0)
fig, ax = plt.subplots(figsize=(8, 5))
bottom = np.zeros(len(comp_frac))
for domain in comp_frac.columns:
    ax.bar(comp_frac.index.astype(str), comp_frac[domain], bottom=bottom, label=str(domain))
    bottom += comp_frac[domain].values
ax.set_xlabel("Slice ID")
ax.set_ylabel("Fraction of bins")
ax.set_title("Domain composition per slice")
ax.legend(title="domain_true", bbox_to_anchor=(1.02, 1), loc="upper left")
fig.tight_layout()
savefig(fig, "06_domain_composition_per_slice.png")

print(f"\nAll diagnostic plots saved under {PLOTS_DIR}")


# ---------------------------------------------------------------------------
# Save the AnnData object
# ---------------------------------------------------------------------------
output_path = sim_app.save(adata, H5AD_PATH)
print(f"\nSaved dataset -> {output_path}")
