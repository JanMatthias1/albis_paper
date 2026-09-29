#!/bin/bash
#SBATCH --job-name=rctd_ref_seed
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/spatial_deconvolution/logs/rctd_ref_seed_%A_%a.out
#SBATCH --array=0-1
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# 2-way seed-replication array (indices 0-1 -> SEEDS below). With the
# canonical seed-2025 reference that makes 3 draws, the project's "3 total
# draws" convention, and mirrors the strong-mix design (2025 + seeds 101/202,
# ../strong_domain/). Seed 271828 was a third independent seed until
# 2026-09-23; dropped by user decision (last in the list, not by result),
# data archived to data/figure_2/misc/archive/rctdref_seed271828_dropped_20260923/.
#
# Independently-seeded single-cell reference for the RCTD spot-deconvolution
# benchmark (run_rctd_spot.R). Not a replacement for Figure 2/3/4's shared
# canonical cell tag (log_mu_-2.5_theta_0.40_jitter0.15_bsigma15, seed=2025,
# the library default) -- this is a SEPARATE, RCTD-only dataset.
#
# Why: RCTD's reference (cell) and query (spot) are already independently
# placed (n_cells 24,207 vs 600,000) and independently dispersion-tuned
# (theta/jitter/batch_sigma/base_gene_lognormal all differ, matching
# Xenium-like vs Visium-like calibration) -- but both nominally use
# --seed 2025 (the CLI default; neither generation script overrides it).
# Traced through albis/albis/simulation_sphere.py
# (simulate_3d_molecule_sphere_multires): a single rng = default_rng(seed)
# is consumed by cell placement and cell-type assignment (both n_cells-
# dependent) BEFORE the per-gene/per-cell-type baseline expression template
# is drawn -- so the ~25x n_cells gap between cell (24,207) and spot
# (600,000) already decorrelates the two templates in practice. That
# decorrelation is an RNG-architecture side effect, not a guarantee: if a
# future change ever brought the two n_cells values closer together, the
# templates could start overlapping again. This script removes that
# dependency entirely by using an explicit, clearly-different --seed, so the
# methods text can state independence outright rather than relying on an
# incidental consequence of draw order.
#
# Every OTHER parameter is copied verbatim from the canonical cell tag
# (log_mu_-2.5_theta_0.40_jitter0.15_bsigma15) -- sphere_r_um, n_cells,
# base_gene_lognormal, theta, theta_jitter, batch_sigma all unchanged.
# Changing those too would confound "is this seed-independent" with "does a
# differently-calibrated reference perform differently" -- matching the
# Xenium-like calibration is the realistic part of this benchmark (a
# reference is supposed to match the tissue's true cell-type vocabulary and
# technology-appropriate noise), not something to perturb here.
#
# Output tag distinguishes this from the shared Figure 2/3/4 cell tag so it
# can never be silently picked up by any other consumer.
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/spatial_deconvolution/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

SEEDS=(999999 314159)
SEED="${SEEDS[${SLURM_ARRAY_TASK_ID:-0}]}"
TAG="log_mu_-2.5_theta_0.40_jitter0.15_bsigma15_rctdref_seed${SEED}"
SIM_RAW="sim_paper/data/figure_2/smaller_sphere/data/${TAG}/simulation_cell_z.h5ad"
SIM_QC="sim_paper/data/figure_2/smaller_sphere/data/${TAG}/simulation_cell_z_qc.h5ad"
NOISY_DIR="sim_paper/data/noisy/${TAG}"

if [[ ! -f "${SIM_RAW}" ]]; then
    echo "[generate] ${TAG} (seed=${SEED})"
    "${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py \
        --modality cell \
        --sphere-r-um 2050 \
        --n-cells 24207 \
        --base-gene-lognormal -2.5 0.7 \
        --theta 0.40 \
        --theta-jitter 0.15 \
        --batch-sigma 1.5 \
        --seed "${SEED}" \
        --sync-unaligned-seed \
        --out-tag "${TAG}"
    mkdir -p sim_paper/data/figure_2/smaller_sphere/data
    mv "${NOISY_DIR}" "sim_paper/data/figure_2/smaller_sphere/data/${TAG}"
fi

if [[ ! -f "${SIM_QC}" ]]; then
    echo "[qc] ${TAG}"
    "${PYTHON_BIN}" sim_paper/code/clustering/step00_qc_filter.py \
        --modality cell --input "${SIM_RAW}" --output "${SIM_QC}"
fi

echo "[done] ${SIM_QC}"
