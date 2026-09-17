#!/bin/bash
#SBATCH --job-name=sweep_spot_logmu_theta_joint_probe
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs/sweep_spot_logmu_theta_joint_probe_%A_%a.out
#SBATCH --array=0-2
#SBATCH --time=08:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# spot joint log_mu x theta mini-grid vs. the two 18k CytAssist probe-based
# real refs (lymph_node_visium, tonsil_visium). Follow-up to
# sweep_spot_logmu_probe.sh + sweep_spot_theta_probe.sh, which between them
# established: (1) log_mu=-2.5 is the best log_mu on total_counts_median /
# matrix_zero_frac at theta=default; (2) dropping theta to ~0.35-0.5 (log_mu
# held at -2.5) matches lymph_node theta_hat but leaves lymph_node
# total_counts_median stuck ~0.55x (sim too shallow) -- a gap theta cannot
# move. This grid tests whether log_mu=-2.0 (which alone overshoots theta_hat)
# combined with a lower theta closes that total_counts gap while keeping
# theta_hat in range -- i.e. whether the real joint optimum is off the two
# 1-D axes swept so far.
#
# Grid: log_mu in {-2.5, -2.0} x theta in {0.25, 0.35, 0.5}. The log_mu=-2.5
# row (theta 0.25/0.35/0.5) is ALREADY COMPUTED by sweep_spot_theta_probe.sh
# -- byte-identical config (same --base-gene-lognormal -2.5 0.7, same
# --batch-sigma 0.3, default --theta-jitter, HVG-matched, slice 5, both refs)
# -- its outputs live at data/count_distribution/sweeps/spot_theta_sweep_probe/
# spot_theta_sweep_{0.25,0.35,0.5}_hvg_matched_vs_{tonsil,lymph_node}_visium/.
# So this script only runs the 3 NEW log_mu=-2.0 points (array 0-2). Tabulate
# all 6 with the companion summary_spot_logmu_theta_joint.sh once these land.
#
# base_gene_lognormal sigma fixed at 0.7, batch_sigma fixed at 0.3, theta_jitter
# left at default (1.0) -- same "only move the axis under test" discipline as the
# two prior probe sweeps. --match-panel-size from the start (HVG-subset the real
# refs to sim's top-556 genes) -- required for an apples-to-apples
# theta_hat/zero_frac read, per the log_mu sweep.
#
# Prereq: real_data_qc/run_visium_tonsil_qc.sh + run_visium_lymph_qc.sh already
# run (data/real_data_qc/{tonsil,lymph_node}_visium/*_qc.h5ad).
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

LOG_MU=-2.0
THETAS=(0.25 0.35 0.5)
THETA="${THETAS[${SLURM_ARRAY_TASK_ID:-0}]}"
TAG="spot_joint_logmu_${LOG_MU}_theta_${THETA}"
MODALITY="spot"

SIM_RAW="sim_paper/data/noisy/${TAG}/simulation_${MODALITY}_z.h5ad"
SIM_QC="sim_paper/data/noisy/${TAG}/simulation_${MODALITY}_z_qc.h5ad"
OUT_ROOT="sim_paper/data/count_distribution/sweeps/spot_joint_logmu_theta_probe"

echo "[config] log_mu=${LOG_MU} theta=${THETA} tag=${TAG}"

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
