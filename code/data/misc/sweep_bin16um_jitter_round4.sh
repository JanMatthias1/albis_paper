#!/bin/bash
#SBATCH --job-name=sweep_bin16um_r4
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs/sweep_bin16um_r4_%A_%a.out
#SBATCH --array=0-2
#SBATCH --time=04:00:00
#SBATCH --mem=320G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Round 4: theta_jitter sweep, testing whether tightening the per-gene theta
# spread (still at its default 1.0 everywhere so far) closes the theta_hat gap
# more efficiently than continuing to push the baseline --theta higher (round
# 3: --theta 2.0->3.0 only moved empirical theta_hat ratio 0.67->0.75, not
# proportionally -- diminishing returns already visible). Same lever that
# unblocked cell's dispersion match at 8um (jitter=1.0 relative to a small
# baseline theta=0.25 was hugely unstable there; here the baseline is much
# larger, but jitter=1.0 around --theta 3.0 is still a wide ~33% relative
# per-gene spread that could be dragging genes down to low effective theta,
# pulling the empirical median down regardless of the baseline).
#
# Fixed at this round's best-so-far combo: --theta 3.0 (round 3, confirmed
# moving theta_hat the right direction), --base-gene-lognormal -1.5 0.7 (best
# QC+HVG-matched total_counts balance from rounds 1-2, see figure.md/session
# notes). Sweeping --theta-jitter 0.15/0.3/0.5 (well below the 1.0 default).

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

JITTERS=(0.15 0.3 0.5)
JITTER="${JITTERS[${SLURM_ARRAY_TASK_ID:-0}]}"
THETA=3.0
LOG_MU=-1.5
TAG="packing_pf0p04_bin16um_theta${THETA}_jitter${JITTER}_log_mu_${LOG_MU}"
MODALITY="bin"

SIM_RAW="sim_paper/data/noisy/${TAG}/simulation_${MODALITY}_z.h5ad"
SIM_QC="sim_paper/data/noisy/${TAG}/simulation_${MODALITY}_z_qc.h5ad"
OUT_ROOT="sim_paper/data/count_distribution/sweeps/bin16um_logmu_sweep"

echo "[config] theta=${THETA} theta_jitter=${JITTER} log_mu=${LOG_MU} tag=${TAG}"

if [[ ! -f "${SIM_RAW}" ]]; then
    echo "[generate] ${TAG} not found, generating"
    "${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py \
        --modality "${MODALITY}" \
        --sphere-r-um 2050 \
        --bin-size-um 16 \
        --theta "${THETA}" \
        --theta-jitter "${JITTER}" \
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
