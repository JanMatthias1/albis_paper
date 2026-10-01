"""
Generate one ALBIS simulation (one modality: cell, bin or spot) for the
manuscript, with its h5ad, a _summary.json and diagnostic plots.

Calls albis.simulate_3d_molecule_sphere_multires() directly rather than
albis.generate_data(), because the count-noise settings (theta, theta_jitter,
noise_scale) are only exposed by the low-level simulator.

The defaults below are a starting point. The settings used in each figure
(gene mean, dispersion, batch effect, tissue size, domain mix, ...) are passed
on the command line by the figure scripts, e.g.
code/count_distribution/smaller_sphere/*.sh.

Output: data/noisy/[<--out-tag>/], or --output-dir if given.

Example:
    python generate_simulation_noisy.py --modality cell --out-tag my_run
    sbatch run_generate_simulation_noisy.sh [extra flags]   # all 3 modalities
"""

import json
import sys
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # headless-safe backend, must be set before pyplot/albis plotting is used
import matplotlib.pyplot as plt
import numpy as np


import albis as ab

REQUIRED_ALBIS_VERSION = "0.1.2"
if ab.__version__ != REQUIRED_ALBIS_VERSION:
    raise RuntimeError(
        f"albis {ab.__version__} found at {ab.__file__}; this script requires albis "
        f"{REQUIRED_ALBIS_VERSION} (pip install albis=={REQUIRED_ALBIS_VERSION})."
    )

print("Python:", sys.executable)
print("albis module:", getattr(ab, "__file__", "<no __file__>"))
print("albis version:", getattr(ab, "__version__", "<no __version__>"))


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

# Manuscript-baseline domain_type_mix (6 domains x 8 cell types) makes several
# domains only weakly distinguishable by composition: domain 5 (core) is
# exactly uniform, domain 4 deviates only 12% from uniform, domain 3 only
# 28% -- pairwise L1 distance between domains ranges just 0.06-0.36. No
# clustering method (plain, BANKSY, or a domain expression shift) can recover
# domains whose true compositions barely differ. This alternative gives every
# domain two strongly-enriched "signature" cell types (0.30 each vs. a 0.125
# uniform baseline) and leaves no domain uniform -- pairwise L1 distance
# becomes 0.467-0.933, i.e. more than 3x the old matrix's *maximum* at its
# *minimum*. Opt-in only (--strong-domain-mix) -- the default stays the
# original matrix so this never silently affects Figure 2 count-distribution
# comparisons, only dedicated Figure 3 domain-recovery runs.
MANUSCRIPT_DOMAIN_TYPE_MIX = np.array(
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
STRONG_DOMAIN_TYPE_MIX = np.array(
    [
        [0.30, 0.30, 0.0667, 0.0667, 0.0667, 0.0667, 0.0667, 0.0667],
        [0.0667, 0.0667, 0.30, 0.30, 0.0667, 0.0667, 0.0667, 0.0667],
        [0.0667, 0.0667, 0.0667, 0.0667, 0.30, 0.30, 0.0667, 0.0667],
        [0.0667, 0.0667, 0.0667, 0.0667, 0.0667, 0.0667, 0.30, 0.30],
        [0.30, 0.0667, 0.0667, 0.0667, 0.30, 0.0667, 0.0667, 0.0667],
        [0.0667, 0.0667, 0.30, 0.0667, 0.0667, 0.0667, 0.0667, 0.30],
    ],
    dtype=float,
)

# Manuscript-baseline geometry (sphere_R_um=6000, n_cells=600_000) packs cells
# at only ~0.16% of sphere volume -- real tissue is essentially fully packed.
# A fixed 8um bin grid over that mostly lands in empty interstitial space:
# 66.3% of bins came back with zero genes detected vs. 0.0% in real
# breast_cancer_visium_hd. Raising theta does not fix this -- tested up to theta=200 with
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
# Capture area = the real instrument's, per modality, independent of how much
# tissue sits in it -- it does NOT scale with --sphere-r-um. bin/spot get the
# 6.5 mm Visium / Visium HD square; cell gets the 12 x 24 mm Xenium window
# (applied via xenium_capture_window_um, which crops cell-level output).
PLATFORM_CAPTURE_WINDOW_UM = {
    "bin": (6500.0, 6500.0),
    "spot": (6500.0, 6500.0),
    "cell": (12000.0, 24000.0),
}
CORE_FUZZ_TO_R = 300.0 / 6000.0
MAX_SHIFT_TO_R = 3000.0 / 6000.0


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate one albis modality for the paper dataset, with higher count noise."
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
        "branch in the mean-variance plot more pronounced.",
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
        "--seed",
        type=int,
        default=2025,
        help="Master seed for simulate_3d_molecule_sphere_multires -- drives the base tissue "
        "(cell placement, domain boundaries, cell types, baseline gene expression) directly, "
        "and the per-slice batch-effect factors indirectly (seed + a fixed per-modality offset). "
        "Does NOT affect obsm['spatial_unaligned'] (see --sync-unaligned-seed / base_seed_unaligned, "
        "seeded independently). Default (2025) matches the library default and all prior runs; "
        "override to build an independent replicate of the same config for variability/error bars.",
    )
    parser.add_argument(
        "--strong-domain-mix",
        action="store_true",
        help="Use STRONG_DOMAIN_TYPE_MIX instead of the manuscript-baseline domain_type_mix -- "
        "opt-in only, for dedicated Figure 3 domain-recovery test runs. Never affects Figure 2 "
        "count-distribution comparisons since it defaults off.",
    )
    parser.add_argument(
        "--domain-size-factors",
        type=float,
        nargs=6,
        default=None,
        metavar=("D0", "D1", "D2", "D3", "D4", "D5"),
        help="Per-domain expression scale factor (6 values, one per domain -- this script's "
        "n_domains=6 is hardcoded below). Multiplies every gene's NB mean uniformly for cells "
        "in that domain (albis's existing domain_size_factors kwarg -- a per-domain library-size "
        "shift, independent of cell type; NOT a per-gene profile the way cell-type markers are). "
        "Defaults to None (all domains =1.0, no effect, current manuscript behavior).",
    )
    parser.add_argument(
        "--sphere-r-um",
        type=float,
        default=NOISY_SPHERE_R_UM,
        help="Sphere radius; with --n-cells fixed, smaller = higher 3D cell-packing "
        "fraction (default 6000um gives ~0.16%% packing, which leaves the "
        "majority of 8um bins empty).",
    )
    parser.add_argument(
        "--n-cells",
        type=int,
        default=NOISY_N_CELLS,
    )
    parser.add_argument(
        "--cell-r-mean",
        type=float,
        default=7.5,
        help="Lognormal mean cell radius (um), which also sets the mean molecule-spillover "
        "radius per cell (sample_molecule_coords_for_cell scales molecule spread to this). "
        "Affects how many genes each bin detects.",
    )
    parser.add_argument(
        "--allow-cell-overlap",
        action="store_true",
        help="Passed through to simulate_3d_molecule_sphere_multires (default False). Needed to "
        "test --cell-r-mean values large enough that non-overlapping placement at the given "
        "--n-cells/--sphere-r-um becomes geometrically infeasible. Already validated elsewhere "
        "in this project as not changing bin-level count statistics at the manuscript radius.",
    )
    parser.add_argument(
        "--bin-size-um",
        type=float,
        default=8.0,
        help="Bin modality grid spacing (um). Manuscript baseline is 8um (Visium-HD-like). "
        "A larger bin aggregates more of each cell's molecule cloud per observation, which "
        "directly addresses two issues diagnosed at 8um: the near-binary empty/non-empty bin "
        "sampling caused by bin size being comparable to --cell-r-mean, and weak "
        "per-bin domain-compositional signal for BANKSY domain recovery (same "
        "aggregation-reveals-composition mechanism as for spot vs. bin/cell). Real Visium HD ships 8um/16um (and "
        "2um) bins, so this is a legitimate alternate real configuration, not just a knob.",
    )
    parser.add_argument(
        "--capture-window-um",
        type=float,
        default=None,
        help="Square capture window side length (um). Default (unset) uses the "
        "real platform window for the modality: 6500 for bin/spot (Visium / "
        "Visium HD), 12000x24000 for cell (Xenium). NOT scaled by --sphere-r-um.",
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
        "--sync-unaligned-seed",
        action="store_true",
        help="Pass sync_unaligned_seed=True through to simulate_3d_molecule_sphere_multires so "
        "obsm['spatial_unaligned']'s per-slice rigid (rotation+translation) perturbation is seeded "
        "identically across modalities (same base_seed_unaligned + slice_id, no per-modality offset) "
        "instead of each modality drawing its own independent perturbation (default). Needed so "
        "cross_tech_stair.py's cross-technology STAIR alignment starts all three modalities from the "
        "SAME misalignment for a given slice_id.",
    )
    parser.add_argument(
        "--out-tag",
        default="",
        help="If set, write output under data/noisy/<out-tag>/ instead of data/noisy/ directly "
        "(e.g. a parameter-sweep label like 'log_mu_-2.5'), so this run doesn't overwrite the "
        "existing simulation_<modality>_<axis>.h5ad.",
    )
    parser.add_argument("--output-dir", type=Path, default=None,
                        help="Explicit output directory; overrides data/noisy and --out-tag.")
    parser.add_argument("--output-stem", default=None,
                        help="File stem for the h5ad, its _summary.json and plots/<stem>/ (default: "
                        "simulation_<modality>_<axis>). Figure 3's strong-mix scripts pass a "
                        "parameter-derived name so the file itself records how it was generated.")
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
DOMAIN_TYPE_MIX = STRONG_DOMAIN_TYPE_MIX if args.strong_domain_mix else None
BASE_GENE_LOGNORMAL = tuple(args.base_gene_lognormal)
N_CELLS = args.n_cells
CELL_R_MEAN = args.cell_r_mean
ALLOW_CELL_OVERLAP = args.allow_cell_overlap
BIN_SIZE_UM = args.bin_size_um
SPHERE_R_UM = args.sphere_r_um
CAPTURE_WINDOW_UM = (
    (args.capture_window_um, args.capture_window_um)
    if args.capture_window_um is not None
    else PLATFORM_CAPTURE_WINDOW_UM[OUTPUT_MODALITY]
)
CORE_FUZZ_WIDTH_UM = args.core_fuzz_width_um if args.core_fuzz_width_um is not None else SPHERE_R_UM * CORE_FUZZ_TO_R
MAX_SHIFT = args.max_shift if args.max_shift is not None else SPHERE_R_UM * MAX_SHIFT_TO_R
SYNC_UNALIGNED_SEED = args.sync_unaligned_seed
SEED = args.seed
OUT_TAG = args.out_tag
DATA_DIR = args.output_dir if args.output_dir is not None else ((BASE_DATA_DIR / OUT_TAG) if OUT_TAG else BASE_DATA_DIR)
DATA_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_STEM = args.output_stem or f"simulation_{OUTPUT_MODALITY}_{SLICE_AXIS.lower()}"
H5AD_PATH = DATA_DIR / f"{OUTPUT_STEM}.h5ad"
SUMMARY_PATH = DATA_DIR / f"{OUTPUT_STEM}_summary.json"
PLOTS_DIR = DATA_DIR / "plots" / OUTPUT_STEM
PLOTS_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Generate data
# ---------------------------------------------------------------------------
def generate_modality(
    output_modality, slice_axis, theta, theta_jitter, noise_scale, batch_sigma, base_gene_lognormal,
    n_cells, cell_r_mean, allow_cell_overlap, sphere_r_um, capture_window_um, core_fuzz_width_um, max_shift,
    marker_foldchange, shared_marker_foldchange, domain_size_factors, domain_type_mix, bin_size_um,
    sync_unaligned_seed, seed,
):
    # ab.generate_data() doesn't forward theta/theta_jitter/noise_scale/
    # base_gene_lognormal, so call the lower-level simulator directly and pull
    # out the requested modality/axis ourselves (mirrors what generate_data()
    # does internally).
    # Platform capture window: bin/spot crop via capture_window_um; cell crops
    # via xenium_capture_window_um.
    capture_kw = (
        {"xenium_capture_window_um": tuple(capture_window_um)}
        if output_modality == "cell"
        else {"capture_window_um": tuple(capture_window_um)}
    )
    simulation = ab.simulate_3d_molecule_sphere_multires(
        sphere_R_um=sphere_r_um,
        **capture_kw,
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
            r_mean=cell_r_mean,
            r_sigma=0.28,
            # r_min/r_max scaled proportionally to cell_r_mean (manuscript ratio: 4.0/7.5,
            # 14.0/7.5) -- sample_cell_radii() clips to [r_min, r_max] regardless of r_mean,
            # so leaving these fixed at the manuscript's absolute values would silently
            # clamp any --cell-r-mean override back down near the old default.
            r_min=cell_r_mean * (4.0 / 7.5),
            r_max=cell_r_mean * (14.0 / 7.5),
        ),
        allow_cell_overlap=allow_cell_overlap,
        n_cell_types=8,
        domain_type_mix=domain_type_mix if domain_type_mix is not None else MANUSCRIPT_DOMAIN_TYPE_MIX,
        n_slices=10,
        bin_size_um=bin_size_um,
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
        sync_unaligned_seed=sync_unaligned_seed,
        seed=seed,
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
        "cell_r_mean": cell_r_mean,
        "allow_cell_overlap": allow_cell_overlap,
        "bin_size_um": bin_size_um,
        "sphere_r_um": sphere_r_um,
        "capture_window_um": list(capture_window_um),
        "core_fuzz_width_um": core_fuzz_width_um,
        "max_shift": max_shift,
        "noise_scale": noise_scale,
        "marker_foldchange": marker_foldchange,
        "shared_marker_foldchange": shared_marker_foldchange,
        "batch_sigma": batch_sigma,
        "base_gene_lognormal": list(base_gene_lognormal),
        "domain_size_factors": list(domain_size_factors) if domain_size_factors is not None else None,
        "strong_domain_mix": domain_type_mix is not None,
        "sync_unaligned_seed": sync_unaligned_seed,
        "seed": seed,
    }
    return adata


print(
    f"\nGenerating {OUTPUT_MODALITY} data for slice axis {SLICE_AXIS} with "
    f"theta={THETA}, theta_jitter={THETA_JITTER}, noise_scale={NOISE_SCALE}, "
    f"batch_sigma={BATCH_SIGMA}, base_gene_lognormal={BASE_GENE_LOGNORMAL}, "
    f"n_cells={N_CELLS}, cell_r_mean={CELL_R_MEAN}, allow_cell_overlap={ALLOW_CELL_OVERLAP}, "
    f"sphere_r_um={SPHERE_R_UM}, capture_window_um={CAPTURE_WINDOW_UM}, "
    f"core_fuzz_width_um={CORE_FUZZ_WIDTH_UM}, max_shift={MAX_SHIFT}, bin_size_um={BIN_SIZE_UM}, "
    f"sync_unaligned_seed={SYNC_UNALIGNED_SEED}, seed={SEED}"
)
adata = generate_modality(
    OUTPUT_MODALITY, SLICE_AXIS, THETA, THETA_JITTER, NOISE_SCALE, BATCH_SIGMA, BASE_GENE_LOGNORMAL,
    N_CELLS, CELL_R_MEAN, ALLOW_CELL_OVERLAP, SPHERE_R_UM, CAPTURE_WINDOW_UM, CORE_FUZZ_WIDTH_UM, MAX_SHIFT,
    MARKER_FOLDCHANGE, SHARED_MARKER_FOLDCHANGE, DOMAIN_SIZE_FACTORS, DOMAIN_TYPE_MIX, BIN_SIZE_UM,
    SYNC_UNALIGNED_SEED, SEED,
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
