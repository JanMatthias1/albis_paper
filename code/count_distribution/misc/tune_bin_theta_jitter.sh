#!/bin/bash
#SBATCH --job-name=tune_bin_theta_jitter
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/logs/tune_bin_theta_jitter_%A_%a.out
#SBATCH --time=08:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#SBATCH --array=0-1
#
# One-off theta/theta_jitter sweep for BIN vs. breast_cancer_visium_hd.
# bin has used the script default (theta=2.0, theta_jitter=1.0) all
# session -- never tuned, unlike cell (theta=0.25, theta_jitter=0.15).
# Testing whether lowering theta (more overdispersion -> higher NB
# zero-probability at fixed mean) closes the panel-matched detection-rate
# gap (sim: 65.67% zeros vs. real: 95.10% zeros, ~7x over-detection) --
# see figure.md, 2026-08-25, "raw_norm_log / genes_per_cell root cause".
# Tradeoff to watch: bin's current empirical theta_hat (~0.4-0.5) is
# already MORE overdispersed than real's target (~1.2), so lowering
# generative theta further may pull dispersion match the wrong way while
# fixing sparsity -- check both qc_filtered (mean_variance) and
# qc_and_hvg_matched (genes_per_cell/sparsity/raw_norm_log) outputs.
# 2 array tasks: 0=theta=1.0/jitter=0.3, 1=theta=0.5/jitter=0.15.
# Not a permanent script -- provisional, output under
# data/count_distribution/sweeps/bin_theta_jitter_sweep/.

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

case "${SLURM_ARRAY_TASK_ID}" in
  0) THETA=1.0; JITTER=0.3 ;;
  1) THETA=0.5; JITTER=0.15 ;;
esac

TAG="packing_pf0p04_log_mu_0.0_theta${THETA}_jitter${JITTER}"
MODALITY="bin"
REAL_LABEL="breast_cancer_visium_hd"
SIM_RAW="sim_paper/data/noisy/${TAG}/simulation_${MODALITY}_z.h5ad"
SIM_QC="sim_paper/data/noisy/${TAG}/simulation_${MODALITY}_z_qc.h5ad"
REAL_INPUT="sim_paper/data/real_data_qc/${REAL_LABEL}/${REAL_LABEL}_qc.h5ad"
OUT_ROOT="sim_paper/data/count_distribution/sweeps/bin_theta_jitter_sweep/${TAG}"

echo "[generate] theta=${THETA} theta_jitter=${JITTER}"
"${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py \
    --modality "${MODALITY}" --sphere-r-um 2050 \
    --base-gene-lognormal 0.0 0.7 \
    --theta "${THETA}" --theta-jitter "${JITTER}" \
    --out-tag "${TAG}"

echo "[qc] filtering"
"${PYTHON_BIN}" sim_paper/code/clustering/00_qc_filter.py --modality "${MODALITY}" --packing-tag "${TAG}"

echo "[compare] qc_filtered (mean_variance, unbiased)"
"${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
    --modality "${MODALITY}" --input "${SIM_QC}" \
    --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
    --output-dir "${OUT_ROOT}_qc_filtered"

echo "[compare] qc_and_hvg_matched (genes_per_cell/sparsity/raw_norm_log)"
"${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
    --modality "${MODALITY}" --input "${SIM_QC}" \
    --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
    --match-panel-size \
    --output-dir "${OUT_ROOT}_qc_hvg_matched"

echo "[done] ${TAG}"
