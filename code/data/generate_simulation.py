"""
Generate an albis synthetic spatial-transcriptomics dataset and save it,
along with diagnostic plots, under sim_paper/data.

Run from anywhere inside the sim_project tree, e.g.:
    cd /dcs04/hicks/data/Jan/sim_project/sim_paper
    python generate_simulation.py
"""

import json
import sys
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # headless-safe backend, must be set before pyplot/albis plotting is used
import matplotlib.pyplot as plt
import numpy as np


# ---------------------------------------------------------------------------
# Locate and import albis as ab (same logic as the tutorial notebook)
# ---------------------------------------------------------------------------
def find_repo_root(start):
    for candidate in [start, *start.parents]:
        if (candidate / "pyproject.toml").is_file() and (candidate / "albis" / "__init__.py").is_file():
            return candidate
        nested = candidate / "albis"
        if (nested / "pyproject.toml").is_file() and (nested / "albis" / "__init__.py").is_file():
            return nested
    raise RuntimeError("Could not find the albis repository root.")


repo_root = find_repo_root(Path.cwd().resolve())
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))

sys.modules.pop("albis", None)
import albis as ab

print("Python:", sys.executable)
print("albis module:", getattr(ab, "__file__", "<no __file__>"))
print("albis version:", getattr(ab, "__version__", "<no __version__>"))


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
DATA_DIR = Path("/dcs04/hicks/data/Jan/sim_project/sim_paper/data")
PLOTS_DIR = DATA_DIR / "plots"
DATA_DIR.mkdir(parents=True, exist_ok=True)
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

SLICE_AXIS = "Z"
VALID_MODALITIES = ("spot", "bin", "cell")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate one albis modality for the paper dataset."
    )
    parser.add_argument(
        "--modality",
        choices=VALID_MODALITIES,
        default="spot",
        help="Output modality to generate. Use SLURM arrays to run these in parallel.",
    )
    parser.add_argument("--slice-axis", choices=("X", "Y", "Z"), default=SLICE_AXIS)
    return parser.parse_args()


args = parse_args()
OUTPUT_MODALITY = args.modality
SLICE_AXIS = args.slice_axis
OUTPUT_STEM = f"simulation_{OUTPUT_MODALITY}_{SLICE_AXIS.lower()}"
H5AD_PATH = DATA_DIR / f"{OUTPUT_STEM}.h5ad"
SUMMARY_PATH = DATA_DIR / f"{OUTPUT_STEM}_summary.json"
PLOTS_DIR = PLOTS_DIR / OUTPUT_STEM
PLOTS_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Generate data
# ---------------------------------------------------------------------------
def generate_modality(output_modality, slice_axis):
    return ab.generate_data(
        output=output_modality,
        slice_axis=slice_axis,
        sphere_radius_um=6000.0,
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
        n_slices=10,
        bin_size_um=8.0,
        batch_sigma=0.22,
        max_deg=270.0,
        max_shift=3000.0,
    )


print(f"\nGenerating {OUTPUT_MODALITY} data for slice axis {SLICE_AXIS}")
adata = generate_modality(OUTPUT_MODALITY, SLICE_AXIS)
print("\nGenerated AnnData:")
print(adata)


# ---------------------------------------------------------------------------
# Inspect / summarize
# ---------------------------------------------------------------------------
summary = ab.describe(adata)
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

MAX_PLOT_OBS = 100_000
if adata.n_obs > MAX_PLOT_OBS:
    rng = np.random.default_rng(2025)
    plot_idx = np.sort(rng.choice(adata.n_obs, size=MAX_PLOT_OBS, replace=False))
    plot_adata = adata[plot_idx].copy()
    print(f"Plotting a deterministic subset of {MAX_PLOT_OBS:,} / {adata.n_obs:,} observations.")
else:
    plot_adata = adata

# 1. Aligned spatial coordinates colored by ground-truth domain
fig = ab.plot(plot_adata, view="2d", coordinates="aligned", color="domain_true", point_size=4)
savefig(fig, "01_aligned_domain_true.png")

# 2. Aligned spatial coordinates colored by ground-truth cell type
fig = ab.plot(plot_adata, view="2d", coordinates="aligned", color="cell_type_true", point_size=4)
savefig(fig, "02_aligned_cell_type_true.png")

# 3. Unaligned (per-slice, pre-registration) coordinates colored by slice id
fig = ab.plot(plot_adata, view="2d", coordinates="unaligned", color="slice_id", point_size=4)
savefig(fig, "03_unaligned_slice_id.png")

# 4. QC: total counts / genes detected per observation
X = adata.X
if hasattr(X, "toarray"):
    total_counts = np.asarray(X.sum(axis=1)).ravel()
    genes_detected = np.asarray((X > 0).sum(axis=1)).ravel()
else:
    total_counts = np.asarray(X.sum(axis=1)).ravel()
    genes_detected = np.asarray((X > 0).sum(axis=1)).ravel()

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
axes[0].hist(total_counts, bins=50, color="steelblue")
axes[0].set_xlabel(f"Total counts per {OUTPUT_MODALITY}")
axes[0].set_ylabel(f"Number of {OUTPUT_MODALITY}s")
axes[0].set_title("Total counts distribution")

axes[1].hist(genes_detected, bins=50, color="darkorange")
axes[1].set_xlabel(f"Genes detected per {OUTPUT_MODALITY}")
axes[1].set_ylabel(f"Number of {OUTPUT_MODALITY}s")
axes[1].set_title("Genes detected distribution")

fig.tight_layout()
savefig(fig, "04_qc_count_distributions.png")

# 5. Number of observations per slice
slice_counts = adata.obs["slice_id"].value_counts().sort_index()
fig, ax = plt.subplots(figsize=(6, 4))
ax.bar(slice_counts.index.astype(str), slice_counts.values, color="slategray")
ax.set_xlabel("Slice ID")
ax.set_ylabel(f"Number of {OUTPUT_MODALITY}s")
ax.set_title(f"{OUTPUT_MODALITY.capitalize()}s per slice")
fig.tight_layout()
savefig(fig, f"05_{OUTPUT_MODALITY}s_per_slice.png")

# 6. Domain composition per slice (stacked bar)
comp = adata.obs.groupby(["slice_id", "domain_true"]).size().unstack(fill_value=0)
comp_frac = comp.div(comp.sum(axis=1), axis=0)
fig, ax = plt.subplots(figsize=(8, 5))
bottom = np.zeros(len(comp_frac))
for domain in comp_frac.columns:
    ax.bar(comp_frac.index.astype(str), comp_frac[domain], bottom=bottom, label=str(domain))
    bottom += comp_frac[domain].values
ax.set_xlabel("Slice ID")
ax.set_ylabel(f"Fraction of {OUTPUT_MODALITY}s")
ax.set_title("Domain composition per slice")
ax.legend(title="domain_true", bbox_to_anchor=(1.02, 1), loc="upper left")
fig.tight_layout()
savefig(fig, "06_domain_composition_per_slice.png")

print(f"\nAll diagnostic plots saved under {PLOTS_DIR}")


# ---------------------------------------------------------------------------
# Save the AnnData object
# ---------------------------------------------------------------------------
output_path = ab.save(adata, H5AD_PATH)
print(f"\nSaved dataset -> {output_path}")
