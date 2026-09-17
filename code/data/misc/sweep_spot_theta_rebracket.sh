#!/bin/bash
#SBATCH --job-name=spot_theta_rebracket
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs/spot_theta_rebracket_%A_%a.out
#SBATCH --array=0-4
#SBATCH --time=08:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Closes the open item flagged in FIGURE2_METHODOLOGY.md / figure.md
# (2026-09-06): the 2026-09-05 joint log_mu x theta sweep that picked
# theta=0.25 ran with the BUGGY theta_jitter=1.0 (fixed to 0.10 the next
# day) -- the floored-gene subpopulation from that bug drags median
# theta_hat down, so theta=0.25 may not be the true optimum now that jitter
# is sane. This re-brackets theta in {0.15,0.20,0.25,0.30,0.35} with
# theta_jitter=0.10 FIXED (the corrected value) and log_mu=-2.0 FIXED (the
# already-settled winner from the log_mu sweep), against both probe refs
# (lymph_node_visium, tonsil_visium), hvg-matched -- same composite scoring
# as summary_spot_logmu_theta_joint.sh (sum of |ln(sim/real ratio)| over
# theta_hat/total_counts_median/matrix_zero_frac, summed across both refs).
#
# If the winner isn't 0.25, promote it into the shared SIM_TAG used by
# spot_vs_{lymph_node,tonsil}_visium.sh, spot_celltype_panel.sh, and
# 3D_stair.py (spot) -- same "byte-identical shared dataset" convention as
# the jitter fix.
#
# Prereq: real_data_qc/run_visium_tonsil_qc.sh + run_visium_lymph_qc.sh
# already run (data/real_data_qc/{tonsil,lymph_node}_visium/*_qc.h5ad).
#
# Tabulate with: sim_paper/code/data/misc/summary_spot_theta_rebracket.sh
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

THETAS=(0.15 0.20 0.25 0.30 0.35)
THETA="${THETAS[${SLURM_ARRAY_TASK_ID:-0}]}"
LOG_MU=-2.0
JITTER=0.10
BATCH_SIGMA=0.3
SPHERE_R_UM=2050
TAG="spot_theta_rebracket_${THETA}"
MODALITY="spot"

SIM_RAW="sim_paper/data/noisy/${TAG}/simulation_${MODALITY}_z.h5ad"
SIM_QC="sim_paper/data/noisy/${TAG}/simulation_${MODALITY}_z_qc.h5ad"
OUT_ROOT="sim_paper/data/count_distribution/sweeps/spot_theta_rebracket_probe"

echo "[config] theta=${THETA} jitter=${JITTER} log_mu=${LOG_MU} batch_sigma=${BATCH_SIGMA} tag=${TAG}"

if [[ ! -f "${SIM_RAW}" ]]; then
    echo "[generate] ${TAG} not found, generating"
    "${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py \
        --modality "${MODALITY}" \
        --sphere-r-um "${SPHERE_R_UM}" \
        --base-gene-lognormal "${LOG_MU}" 0.7 \
        --theta "${THETA}" \
        --theta-jitter "${JITTER}" \
        --batch-sigma "${BATCH_SIGMA}" \
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
