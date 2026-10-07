#!/bin/bash
#SBATCH --job-name=rctd_ref_seed
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/albis_paper/code/clustering/spatial_deconvolution/logs/rctd_ref_seed_%A_%a.out
#SBATCH --array=0-1
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Independently seeded single-cell references for the RCTD spot benchmark
# (run_rctd_spot.R): array tasks 0-1 = the seeds below; with the canonical
# seed-2025 reference this gives 3 references, mirroring the strong-mix design.
# RCTD-only data: every parameter equals the canonical Figure 2 cell settings
# (log_mu_-2.5_theta_0.40_jitter0.15_bsigma15) except --seed. The canonical
# reference and the spot query both use seed 2025; their expression templates
# differ only because their cell counts differ, so an explicit different seed
# makes reference and query independent by construction. The output tag keeps
# these files apart from the shared Figure 2/3/4 cell data.
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/albis_paper/code/clustering/spatial_deconvolution/logs
source /dcs04/hicks/data/Jan/sim_project/albis_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

SEEDS=(999999 314159)
SEED="${SEEDS[${SLURM_ARRAY_TASK_ID:-0}]}"
TAG="log_mu_-2.5_theta_0.40_jitter0.15_bsigma15_rctdref_seed${SEED}"
SIM_RAW="albis_paper/data/figure_2/smaller_sphere/data/${TAG}/simulation_cell_z.h5ad"
SIM_QC="albis_paper/data/figure_2/smaller_sphere/data/${TAG}/simulation_cell_z_qc.h5ad"
NOISY_DIR="albis_paper/data/noisy/${TAG}"

if [[ ! -f "${SIM_RAW}" ]]; then
    echo "[generate] ${TAG} (seed=${SEED})"
    "${PYTHON_BIN}" albis_paper/code/data/generate_simulation_noisy.py \
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
    mkdir -p albis_paper/data/figure_2/smaller_sphere/data
    mv "${NOISY_DIR}" "albis_paper/data/figure_2/smaller_sphere/data/${TAG}"
fi

if [[ ! -f "${SIM_QC}" ]]; then
    echo "[qc] ${TAG}"
    "${PYTHON_BIN}" albis_paper/code/clustering/step00_qc_filter.py \
        --modality cell --input "${SIM_RAW}" --output "${SIM_QC}"
fi

echo "[done] ${SIM_QC}"
