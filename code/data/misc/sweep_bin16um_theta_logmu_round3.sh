#!/bin/bash
#SBATCH --job-name=sweep_bin16um_r3
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs/sweep_bin16um_r3_%A_%a.out
#SBATCH --array=0-5
#SBATCH --time=04:00:00
#SBATCH --mem=320G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Round 3 of the bin_size_um=16 calibration sweep: 2D grid over generative
# theta x log_mu, following two corrections from rounds 1-2:
#
# 1. theta_hat was essentially FLAT across all round-1/2 log_mu values
#    (1.11-1.15, from the trustworthy qc_filtered/full-panel mode -- HVG-matched
#    theta_hat is NOT trustworthy, HVG selection on the real side biases it down,
#    same caveat as the 8um work) -- confirms theta_hat is controlled by
#    generative --theta/--theta-jitter, not --base-gene-lognormal. Real is
#    1.67-1.72 vs sim's ~1.13 at the default --theta 2.0 (empirical/generative
#    ratio ~0.57) -- extrapolating that ratio, --theta 3.0 should land empirical
#    theta_hat near real. Testing 3.0 and 4.0 to bracket in case the relationship
#    isn't linear (it often isn't -- see the cell theta/theta_hat tuning history).
# 2. total_counts must be judged on the QC+HVG-matched comparison (556 sim genes
#    vs. real's full 18,085-gene whole-transcriptome panel makes the unmatched
#    comparison unfair) -- under that fair comparison, round 1/2's best point was
#    actually log_mu=-1.5 (0.86x pancreas, 1.97x breast_cancer), NOT the higher
#    log_mu values that looked better only under the unmatched comparison.
#    Narrowing to -2.5/-2.0/-1.5 to bracket log_mu=-1.5 more finely.
#
# Runs BOTH qc_filtered (theta_hat/mean_variance) and qc_and_hvg_matched
# (total_counts/genes_per_cell) against both real 16um references for every
# point, per the QC+HVG standard going forward.

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

THETAS=(3.0 3.0 3.0 4.0 4.0 4.0)
LOG_MUS=(-2.5 -2.0 -1.5 -2.5 -2.0 -1.5)
IDX="${SLURM_ARRAY_TASK_ID:-0}"
THETA="${THETAS[$IDX]}"
LOG_MU="${LOG_MUS[$IDX]}"
TAG="packing_pf0p04_bin16um_theta${THETA}_log_mu_${LOG_MU}"
MODALITY="bin"

SIM_RAW="sim_paper/data/noisy/${TAG}/simulation_${MODALITY}_z.h5ad"
SIM_QC="sim_paper/data/noisy/${TAG}/simulation_${MODALITY}_z_qc.h5ad"
OUT_ROOT="sim_paper/data/count_distribution/sweeps/bin16um_logmu_sweep"

echo "[config] theta=${THETA} log_mu=${LOG_MU} tag=${TAG}"

if [[ ! -f "${SIM_RAW}" ]]; then
    echo "[generate] ${TAG} not found, generating"
    "${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py \
        --modality "${MODALITY}" \
        --sphere-r-um 2050 \
        --bin-size-um 16 \
        --theta "${THETA}" \
        --base-gene-lognormal "${LOG_MU}" 0.7 \
        --out-tag "${TAG}"
fi

if [[ ! -f "${SIM_QC}" ]]; then
    echo "[qc] ${SIM_QC} not found, running 00_qc_filter.py"
    "${PYTHON_BIN}" sim_paper/code/clustering/00_qc_filter.py --modality "${MODALITY}" --packing-tag "${TAG}"
fi

for real in breast_cancer_visium_hd human_pancreas_visium_hd; do
    echo "[compare] ${TAG} vs ${real}_16um (qc_filtered)"
    "${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
        --modality "${MODALITY}" --input "${SIM_QC}" \
        --compare-input "sim_paper/data/real_data_qc/${real}_16um/${real}_qc.h5ad" \
        --compare-label "${real}_16um" \
        --output-dir "${OUT_ROOT}/${TAG}_qc_filtered_vs_${real}_16um"

    echo "[compare] ${TAG} vs ${real}_16um (qc_and_hvg_matched)"
    "${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
        --modality "${MODALITY}" --input "${SIM_QC}" \
        --compare-input "sim_paper/data/real_data_qc/${real}_16um/${real}_qc.h5ad" \
        --compare-label "${real}_16um" \
        --match-panel-size \
        --output-dir "${OUT_ROOT}/${TAG}_qc_and_hvg_matched_vs_${real}_16um"
done

echo "[done] ${TAG}"
