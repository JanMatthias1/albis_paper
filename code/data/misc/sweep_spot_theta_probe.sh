#!/bin/bash
#SBATCH --job-name=sweep_spot_theta_probe
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs/sweep_spot_theta_probe_%A_%a.out
#SBATCH --array=0-5
#SBATCH --time=08:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# spot theta (NB dispersion) bracketing sweep vs. the two 18k CytAssist
# probe-based real refs (lymph_node_visium, tonsil_visium). Follow-up to
# sweep_spot_logmu_probe.sh: once real refs are HVG-matched down to sim's
# 556-gene panel (--match-panel-size), log_mu=-2.5 (the pre-existing
# incumbent) turned out to already be the best log_mu match on
# total_counts_median/matrix_zero_frac -- but theta_hat stayed off by 2-4x
# (sim 1.43 vs real 0.34 tonsil / 0.76 lymph_node) across the ENTIRE log_mu
# grid, since log_mu shifts the mean, not dispersion. This sweep isolates
# that gap by holding log_mu fixed and moving theta down from the 2.0
# default (lower theta = more overdispersion = lower empirical theta_hat).
#
# base_gene_lognormal fixed at (-2.5, 0.7) and batch_sigma fixed at 0.3
# (both the winning spot log_mu config) -- only --theta varies. theta_jitter
# left at its default (1.0), same "don't touch jitter" precedent as the
# bin16um log_mu sweep.
#
# Uses --match-panel-size from the start this time (HVG-subset the real refs
# down to top-556 genes before computing stats) -- the log_mu sweep showed
# this is required for an apples-to-apples theta_hat/zero_frac read; the
# sim's 556-gene panel vs real's full ~18k-gene panel comparison is not
# apples-to-apples for dispersion.
#
# Prereq: real_data_qc/run_visium_tonsil_qc.sh + run_visium_lymph_qc.sh
# already run (data/real_data_qc/{tonsil,lymph_node}_visium/*_qc.h5ad).
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

THETAS=(1.5 1.0 0.75 0.5 0.35 0.25)
THETA="${THETAS[${SLURM_ARRAY_TASK_ID:-0}]}"
LOG_MU=-2.5
TAG="spot_theta_sweep_${THETA}"
MODALITY="spot"

SIM_RAW="sim_paper/data/noisy/${TAG}/simulation_${MODALITY}_z.h5ad"
SIM_QC="sim_paper/data/noisy/${TAG}/simulation_${MODALITY}_z_qc.h5ad"
OUT_ROOT="sim_paper/data/count_distribution/sweeps/spot_theta_sweep_probe"

echo "[config] theta=${THETA} log_mu=${LOG_MU} tag=${TAG}"

if [[ ! -f "${SIM_RAW}" ]]; then
    echo "[generate] ${TAG} not found, generating"
    "${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py \
        --modality "${MODALITY}" \
        --sphere-r-um 2050 \
        --base-gene-lognormal "${LOG_MU}" 0.7 \
        --theta "${THETA}" \
        --batch-sigma 0.3 \
        --out-tag "${TAG}"
fi

if [[ ! -f "${SIM_QC}" ]]; then
    echo "[qc] ${SIM_QC} not found, running 00_qc_filter.py"
    "${PYTHON_BIN}" sim_paper/code/clustering/00_qc_filter.py --modality "${MODALITY}" --packing-tag "${TAG}"
fi

for real in tonsil_visium lymph_node_visium; do
    echo "[compare] ${TAG} vs ${real} (hvg-matched)"
    "${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
        --modality "${MODALITY}" --input "${SIM_QC}" \
        --compare-input "sim_paper/data/real_data_qc/${real}/${real}_qc.h5ad" \
        --compare-label "${real}" \
        --match-panel-size \
        --output-dir "${OUT_ROOT}/${TAG}_hvg_matched_vs_${real}"
done

echo "[done] ${TAG}"
