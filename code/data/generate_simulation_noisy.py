"""
Generate a sim_app synthetic spatial-transcriptomics dataset with higher
per-gene count noise (lower NB dispersion), along with diagnostic plots,
under sim_paper/data/noisy.

Identical to generate_simulation.py except this bypasses sim_app.generate_data()
and calls sim_app.simulate_3d_molecule_sphere_multires() directly, since the
NB dispersion knobs (theta, theta_jitter, noise_scale) are not forwarded by
generate_data() -- see sim_app/README.md, "Low-level simulator".

Defaults were picked by comparing against real Xenium data with
sim_paper/code/count_distribution/count_distribution.py --compare-input:
the manuscript-baseline simulation (theta=25 input) fits an empirical
theta_hat ~= 1.56, while real slices (non_diseased_lung, lung_cancer) fit
theta_hat ~= 0.15-0.16 -- i.e. real data is far more overdispersed. These
defaults (theta=2.0, theta_jitter=1.0, noise_scale=1.3) are a first pass at
closing that gap; rerun count_distribution.py --compare-input against this
output to check the resulting theta_hat and retune --theta if needed.

Run from anywhere inside the sim_project tree, e.g.:
    cd /dcs04/hicks/data/Jan/sim_project/sim_paper
    python generate_simulation_noisy.py
"""

import json
import sys
import argparse
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

print("Python:", sys.executable)
print("sim_app module:", getattr(sim_app, "__file__", "<no __file__>"))
print("sim_app version:", getattr(sim_app, "__version__", "<no __version__>"))


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE_DATA_DIR = Path("/dcs04/hicks/data/Jan/sim_project/sim_paper/data/noisy")

SLICE_AXIS = "Z"
VALID_MODALITIES = ("spot", "bin", "cell")

# Manuscript baseline uses theta=25 (default), theta_jitter=2.0, noise_scale=0.9.
# These push dispersion down / background noise up to better match real data
# (see module docstring for the theta_hat comparison that motivated them).
NOISY_THETA = 2.0
NOISY_THETA_JITTER = 1.0
NOISY_NOISE_SCALE = 1.3
MANUSCRIPT_BATCH_SIGMA = 0.22
# Manuscript baseline (and simulate_3d_molecule_sphere_multires's own default)
# is base_gene_lognormal=(0.7, 0.7) -- median per-gene baseline expression
# exp(0.7) ~= 2.0. Lower the first value (log_mu) to bring down average
# per-cell total counts / genes detected without touching dispersion/noise.
NOISY_BASE_GENE_LOGNORMAL = (0.7, 0.7)

# Manuscript-baseline geometry (sphere_R_um=6000, n_cells=600_000) packs cells
# at only ~0.16% of sphere volume -- real tissue is essentially fully packed.
# A fixed 8um bin grid over that mostly lands in empty interstitial space:
# 66.3% of bins came back with zero genes detected vs. 0.0% in real
# breast_cancer_visium_hd (see DATA_VERSIONS.md, "Known issue" section,
# 2026-08-20). Raising theta does not fix this -- tested up to theta=200 with
# almost no effect -- because it's geometric (empty grid cells), not
# per-molecule count noise. Local validation swept target 3D packing fraction
# via --sphere-r-um (holding --n-cells fixed): packing fully solves emptiness
# by ~10% but overshoots real counts/bin by 5-9x, because it shrinks cell
# spacing enough that neighboring cells' molecule clouds start overlapping
# into the same bin. There's no packing fraction that hits both targets
# exactly; ~0.8-1.5% packing was the best compromise found (empty bins cut
# from 66% to ~7-18%, counts/bin within ~1.1-1.6x of real). See
# packing_fraction_to_sphere_R() below for the conversion.
NOISY_N_CELLS = 600_000
NOISY_SPHERE_R_UM = 6000.0
WINDOW_TO_R = 6500.0 / 6000.0
CORE_FUZZ_TO_R = 300.0 / 6000.0
MAX_SHIFT_TO_R = 3000.0 / 6000.0


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate one sim_app modality for the paper dataset, with higher count noise."
    )
    parser.add_argument(
        "--modality",
        choices=VALID_MODALITIES,
        default="spot",
        help="Output modality to generate. Use SLURM arrays to run these in parallel.",
    )
    parser.add_argument("--slice-axis", choices=("X", "Y", "Z"), default=SLICE_AXIS)
    parser.add_argument(
        "--theta",
        type=float,
        default=NOISY_THETA,
        help="NB dispersion (variance = mean + mean^2/theta); lower = noisier.",
    )
    parser.add_argument(
        "--theta-jitter",
        type=float,
        default=NOISY_THETA_JITTER,
        help="Per-gene spread of theta around --theta.",
    )
    parser.add_argument(
        "--noise-scale",
        type=float,
        default=NOISY_NOISE_SCALE,
        help="Multiplier applied to non-marker 'noise gene' expression.",
    )
    parser.add_argument(
        "--marker-foldchange",
        type=float,
        default=3.5,
        help="Expression multiplier for a cell type's own unique marker genes. Pooling variance "
        "across all cell types inflates a marker gene's variance far above the NB fit (between-type "
        "mean differences dominate over within-type NB variance) -- higher values make that upper "
        "branch in the mean-variance plot more pronounced; see the 2026-08-24 count_distribution "
        "mean-variance investigation.",
    )
    parser.add_argument(
        "--shared-marker-foldchange",
        type=float,
        default=2.5,
        help="Expression multiplier for markers shared across all cell types (same variance-inflation "
        "caveat as --marker-foldchange, but shared markers don't differ between types so the effect "
        "is much smaller in practice).",
    )
    parser.add_argument(
        "--base-gene-lognormal",
        type=float,
        nargs=2,
        metavar=("LOG_MU", "LOG_SIGMA"),
        default=NOISY_BASE_GENE_LOGNORMAL,
        help="Lognormal (log_mu, log_sigma) for per-gene baseline expression "
        "(base_gene = exp(Normal(log_mu, log_sigma))); lower log_mu = lower "
        "average per-cell total counts / genes detected.",
    )
    parser.add_argument(
        "--batch-sigma",
        type=float,
        default=MANUSCRIPT_BATCH_SIGMA,
        help="Standard deviation of the per-slice, per-gene log-fold-change batch effect.",
    )
    parser.add_argument(
        "--domain-size-factors",
        type=float,
        nargs=6,
        default=None,
        metavar=("D0", "D1", "D2", "D3", "D4", "D5"),
        help="Per-domain expression scale factor (6 values, one per domain -- this script's "
        "n_domains=6 is hardcoded below). Multiplies every gene's NB mean uniformly for cells "
        "in that domain (sim_app's existing domain_size_factors kwarg -- a per-domain library-size "
        "shift, independent of cell type; NOT a per-gene profile the way cell-type markers are). "
        "Defaults to None (all domains =1.0, no effect, current manuscript behavior).",
    )
    parser.add_argument(
        "--sphere-r-um",
        type=float,
        default=NOISY_SPHERE_R_UM,
        help="Sphere radius; with --n-cells fixed, smaller = higher 3D cell-packing "
        "fraction (see count_distribution bin/spot zero-inflation investigation, "
        "2026-08-20 -- default 6000um gives ~0.16%% packing, which leaves the "
        "majority of 8um bins empty).",
    )
    parser.add_argument(
        "--n-cells",
        type=int,
        default=NOISY_N_CELLS,
    )
    parser.add_argument(
        "--capture-window-um",
        type=float,
        default=None,
        help="Square capture window side length. Defaults to "
        f"sphere_r_um * {WINDOW_TO_R:.4f} (same ratio as the manuscript baseline: "
        "6500/6000).",
    )
    parser.add_argument(
        "--core-fuzz-width-um",
        type=float,
        default=None,
        help=f"Defaults to sphere_r_um * {CORE_FUZZ_TO_R:.4f} (manuscript ratio: 300/6000).",
    )
    parser.add_argument(
        "--max-shift",
        type=float,
        default=None,
        help=f"Defaults to sphere_r_um * {MAX_SHIFT_TO_R:.4f} (manuscript ratio: 3000/6000).",
    )
    parser.add_argument(
        "--out-tag",
        default="",
        help="If set, write output under data/noisy/<out-tag>/ instead of data/noisy/ directly "
        "(e.g. a parameter-sweep label like 'log_mu_-2.5'), so this run doesn't overwrite the "
        "existing simulation_<modality>_<axis>.h5ad.",
    )
    return parser.parse_args()


args = parse_args()
OUTPUT_MODALITY = args.modality
SLICE_AXIS = args.slice_axis
THETA = args.theta
THETA_JITTER = args.theta_jitter
NOISE_SCALE = args.noise_scale
MARKER_FOLDCHANGE = args.marker_foldchange
SHARED_MARKER_FOLDCHANGE = args.shared_marker_foldchange
BATCH_SIGMA = args.batch_sigma
DOMAIN_SIZE_FACTORS = args.domain_size_factors
BASE_GENE_LOGNORMAL = tuple(args.base_gene_lognormal)
N_CELLS = args.n_cells
SPHERE_R_UM = args.sphere_r_um
CAPTURE_WINDOW_UM = args.capture_window_um if args.capture_window_um is not None else SPHERE_R_UM * WINDOW_TO_R
CORE_FUZZ_WIDTH_UM = args.core_fuzz_width_um if args.core_fuzz_width_um is not None else SPHERE_R_UM * CORE_FUZZ_TO_R
MAX_SHIFT = args.max_shift if args.max_shift is not None else SPHERE_R_UM * MAX_SHIFT_TO_R
OUT_TAG = args.out_tag
DATA_DIR = (BASE_DATA_DIR / OUT_TAG) if OUT_TAG else BASE_DATA_DIR
DATA_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_STEM = f"simulation_{OUTPUT_MODALITY}_{SLICE_AXIS.lower()}"
H5AD_PATH = DATA_DIR / f"{OUTPUT_STEM}.h5ad"
SUMMARY_PATH = DATA_DIR / f"{OUTPUT_STEM}_summary.json"
PLOTS_DIR = DATA_DIR / "plots" / OUTPUT_STEM
PLOTS_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Generate data
# ---------------------------------------------------------------------------
def generate_modality(
    output_modality, slice_axis, theta, theta_jitter, noise_scale, batch_sigma, base_gene_lognormal,
    n_cells, sphere_r_um, capture_window_um, core_fuzz_width_um, max_shift,
    marker_foldchange, shared_marker_foldchange, domain_size_factors,
):
    # sim_app.generate_data() doesn't forward theta/theta_jitter/noise_scale/
    # base_gene_lognormal, so call the lower-level simulator directly and pull
    # out the requested modality/axis ourselves (mirrors what generate_data()
    # does internally).
    simulation = sim_app.simulate_3d_molecule_sphere_multires(
        sphere_R_um=sphere_r_um,
        capture_window_um=(capture_window_um, capture_window_um),
        n_domains=6,
        core_frac=0.55,
        core_bump_amp=0.25,
        wedge_angle_amp_deg=25.0,

        noise_terms=16,
        noise_freq_range=(3.0, 6.0),
        boundary_fuzz_width_deg=6.0,
        boundary_fuzz_flip_prob=0.15,
        core_fuzz_width_um=core_fuzz_width_um,
        core_fuzz_flip_prob=0.25,

        n_cells=n_cells,
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
        theta=theta,
        theta_jitter=theta_jitter,
        noise_scale=noise_scale,
        marker_foldchange=marker_foldchange,
        shared_marker_foldchange=shared_marker_foldchange,
        base_gene_lognormal=base_gene_lognormal,
        batch_sigma=batch_sigma,
        domain_size_factors=domain_size_factors,
        max_deg=270.0,
        max_shift=max_shift,
        output_modalities=(output_modality,),
        slice_axes=(slice_axis,),
    )

    output_key = {"cell": "adata_cell_sectioned", "bin": "bin_adatas", "spot": "spot_adatas"}[output_modality]
    adata = simulation[output_key][slice_axis]
    adata.uns["output"] = {"platform": output_modality, "slice_axis": slice_axis}
    adata.uns["sim_params"] = {
        "theta": theta,
        "theta_jitter": theta_jitter,
        "n_cells": n_cells,
        "sphere_r_um": sphere_r_um,
        "capture_window_um": capture_window_um,
        "core_fuzz_width_um": core_fuzz_width_um,
        "max_shift": max_shift,
        "noise_scale": noise_scale,
        "marker_foldchange": marker_foldchange,
        "shared_marker_foldchange": shared_marker_foldchange,
        "batch_sigma": batch_sigma,
        "base_gene_lognormal": list(base_gene_lognormal),
        "domain_size_factors": list(domain_size_factors) if domain_size_factors is not None else None,
    }
    return adata


print(
    f"\nGenerating {OUTPUT_MODALITY} data for slice axis {SLICE_AXIS} with "
    f"theta={THETA}, theta_jitter={THETA_JITTER}, noise_scale={NOISE_SCALE}, "
    f"batch_sigma={BATCH_SIGMA}, base_gene_lognormal={BASE_GENE_LOGNORMAL}, "
    f"n_cells={N_CELLS}, sphere_r_um={SPHERE_R_UM}, capture_window_um={CAPTURE_WINDOW_UM}, "
    f"core_fuzz_width_um={CORE_FUZZ_WIDTH_UM}, max_shift={MAX_SHIFT}"
)
adata = generate_modality(
    OUTPUT_MODALITY, SLICE_AXIS, THETA, THETA_JITTER, NOISE_SCALE, BATCH_SIGMA, BASE_GENE_LOGNORMAL,
    N_CELLS, SPHERE_R_UM, CAPTURE_WINDOW_UM, CORE_FUZZ_WIDTH_UM, MAX_SHIFT,
    MARKER_FOLDCHANGE, SHARED_MARKER_FOLDCHANGE, DOMAIN_SIZE_FACTORS,
)
print("\nGenerated AnnData:")
print(adata)


# ---------------------------------------------------------------------------
# Quick empirical dispersion check (same method as count_distribution.py)
# ---------------------------------------------------------------------------
def fit_empirical_theta(X):
    if hasattr(X, "multiply"):
        mean = np.asarray(X.mean(axis=0)).ravel()
        mean_sq = np.asarray(X.multiply(X).mean(axis=0)).ravel()
    else:
        mean = X.mean(axis=0)
        mean_sq = (X**2).mean(axis=0)
    var = mean_sq - mean**2
    overdispersed = var > mean
    theta_hat = mean[overdispersed] ** 2 / (var[overdispersed] - mean[overdispersed])
    theta_hat = theta_hat[np.isfinite(theta_hat) & (theta_hat > 0)]
    return float(np.median(theta_hat)) if theta_hat.size else float("nan")


raw_counts = adata.layers["counts_pre_batch"] if "counts_pre_batch" in adata.layers else adata.X
empirical_theta = fit_empirical_theta(raw_counts)
print(
    f"\n[nb] empirical theta_hat = {empirical_theta:.2f} "
    "(real data fits ~0.15-0.16 via count_distribution.py --compare-input; "
    "lower --theta further if this is still much higher)"
)


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

MAX_PLOT_OBS = 100_000
if adata.n_obs > MAX_PLOT_OBS:
    rng = np.random.default_rng(2025)
    plot_idx = np.sort(rng.choice(adata.n_obs, size=MAX_PLOT_OBS, replace=False))
    plot_adata = adata[plot_idx].copy()
    print(f"Plotting a deterministic subset of {MAX_PLOT_OBS:,} / {adata.n_obs:,} observations.")
else:
    plot_adata = adata

# 1. Aligned spatial coordinates colored by ground-truth domain
fig = sim_app.plot(plot_adata, view="2d", coordinates="aligned", color="domain_true", point_size=4)
savefig(fig, "01_aligned_domain_true.png")

# 2. Aligned spatial coordinates colored by ground-truth cell type
fig = sim_app.plot(plot_adata, view="2d", coordinates="aligned", color="cell_type_true", point_size=4)
savefig(fig, "02_aligned_cell_type_true.png")

# 3. Unaligned (per-slice, pre-registration) coordinates colored by slice id
fig = sim_app.plot(plot_adata, view="2d", coordinates="unaligned", color="slice_id", point_size=4)
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
output_path = sim_app.save(adata, H5AD_PATH)
print(f"\nSaved dataset -> {output_path}")
