#!/bin/bash
#SBATCH --job-name=sweep_bin16um_logmu_r2
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs/sweep_bin16um_logmu_r2_%A_%a.out
#SBATCH --array=0-1
#SBATCH --time=08:00:00
#SBATCH --mem=320G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Round 2 of the bin_size_um=16 log_mu sweep (see sweep_bin16um_logmu.sh for
# full rationale). Round 1 (log_mu -2.5/-1.5/-0.5) landed with theta_hat
# already much closer to real than 8um ever achieved (ratio ~0.66-0.69 vs.
# 8um's best ~0.3-0.4) and empty-bin fraction resolved (0% post-QC, ~4.5%
# pre-QC at the sparsest point, vs. 8um's 9.1%) -- confirms the 16um
# fragmentation hypothesis. But total_counts is still far under real at every
# round-1 point (12-59x low), and scales ~exp(log_mu) as expected (measured:
# 131 -> 352 from -2.5 -> -1.5, matching e^1). Extrapolating that trend to
# breast_cancer_visium_hd_16um's real median (1572.5) lands log_mu near 0.0
# -- notably close to 8um's own eventual winner (log_mu=0.0), not the shifted
# -down value originally hypothesized from the 4x bin-area difference alone.
# This round brackets that: log_mu 0.0 and 0.5.

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

LOG_MUS=(0.0 0.5)
LOG_MU="${LOG_MUS[${SLURM_ARRAY_TASK_ID:-0}]}"
TAG="packing_pf0p04_bin16um_log_mu_${LOG_MU}"
MODALITY="bin"

SIM_RAW="sim_paper/data/noisy/${TAG}/simulation_${MODALITY}_z.h5ad"
SIM_QC="sim_paper/data/noisy/${TAG}/simulation_${MODALITY}_z_qc.h5ad"
OUT_ROOT="sim_paper/data/count_distribution/sweeps/bin16um_logmu_sweep"

echo "[config] log_mu=${LOG_MU} tag=${TAG}"

if [[ ! -f "${SIM_RAW}" ]]; then
    echo "[generate] ${TAG} not found, generating"
    "${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py \
        --modality "${MODALITY}" \
        --sphere-r-um 2050 \
        --bin-size-um 16 \
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
done

echo "[done] ${TAG}"
