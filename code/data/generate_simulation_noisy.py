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
BASE_DATA_DIR = Path("/dcs04/hicks/data/Jan/sim_project/albis_paper/data/noisy")

SLICE_AXIS = "Z"
VALID_MODALITIES = ("spot", "bin", "cell")

# Defaults only; the figure scripts pass their own values.
NOISY_THETA = 2.0
NOISY_THETA_JITTER = 1.0
NOISY_NOISE_SCALE = 1.3
MANUSCRIPT_BATCH_SIGMA = 0.22
NOISY_BASE_GENE_LOGNORMAL = (0.7, 0.7)

# Domain x cell-type composition (6 domains x 8 cell types): weak baseline mix,
# and the strong mix (two enriched cell types per domain) used with --strong-domain-mix.
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

NOISY_N_CELLS = 600_000
NOISY_SPHERE_R_UM = 6000.0
# Real instrument capture area per modality (Visium / Visium HD, Xenium); not scaled with sphere radius.
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
        help="Expression multiplier for a cell type's own marker genes.",
    )
    parser.add_argument(
        "--shared-marker-foldchange",
        type=float,
        default=2.5,
        help="Expression multiplier for markers shared across all cell types.",
    )
    parser.add_argument(
        "--base-gene-lognormal",
        type=float,
        nargs=2,
        metavar=("LOG_MU", "LOG_SIGMA"),
        default=NOISY_BASE_GENE_LOGNORMAL,
        help="Lognormal (log_mu, log_sigma) of per-gene baseline expression; lower log_mu = fewer counts.",
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
        help="Master seed (tissue, cell types, expression, batch factors); change for replicates.",
    )
    parser.add_argument(
        "--strong-domain-mix",
        action="store_true",
        help="Use STRONG_DOMAIN_TYPE_MIX instead of the baseline domain x cell-type mix.",
    )
    parser.add_argument(
        "--domain-size-factors",
        type=float,
        nargs=6,
        default=None,
        metavar=("D0", "D1", "D2", "D3", "D4", "D5"),
        help="Per-domain library-size factor (one per domain); default: no effect.",
    )
    parser.add_argument(
        "--sphere-r-um",
        type=float,
        default=NOISY_SPHERE_R_UM,
        help="Sphere radius (um); with --n-cells fixed, smaller = denser tissue.",
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
        help="Lognormal mean cell radius (um); also sets each cell's molecule spread.",
    )
    parser.add_argument(
        "--allow-cell-overlap",
        action="store_true",
        help="Allow overlapping cells during placement (default off; off for all manuscript data).",
    )
    parser.add_argument(
        "--bin-size-um",
        type=float,
        default=8.0,
        help="Bin grid spacing (um), e.g. 8 or 16 as in Visium HD.",
    )
    parser.add_argument(
        "--capture-window-um",
        type=float,
        default=None,
        help="Square capture window side (um). Default: platform window "
        "(6500 for bin/spot, 12000x24000 for cell).",
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
        help="Use the same per-slice misalignment (obsm['spatial_unaligned']) across modalities.",
    )
    parser.add_argument(
        "--out-tag",
        default="",
        help="Write output under data/noisy/<out-tag>/ instead of data/noisy/.",
    )
    parser.add_argument("--output-dir", type=Path, default=None,
                        help="Explicit output directory; overrides data/noisy and --out-tag.")
    parser.add_argument("--output-stem", default=None,
                        help="File stem for the h5ad, _summary.json and plots/<stem>/ "
                        "(default: simulation_<modality>_<axis>).")
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
    # bin/spot crop via capture_window_um; cell via xenium_capture_window_um.
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
            # clip bounds scale with cell_r_mean (ratios of the 7.5 um default)
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
# Diagnostic plots (simulation_diagnostic_plots.py can redraw them from the saved h5ad)
# ---------------------------------------------------------------------------
from simulation_diagnostic_plots import save_diagnostic_plots

save_diagnostic_plots(adata, PLOTS_DIR, OUTPUT_MODALITY)


# ---------------------------------------------------------------------------
# Save the AnnData object
# ---------------------------------------------------------------------------
output_path = ab.save(adata, H5AD_PATH)
print(f"\nSaved dataset -> {output_path}")
