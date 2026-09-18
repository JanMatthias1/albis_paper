#!/bin/bash
#SBATCH --job-name=sweep_cell_logmu_theta_joint_xenium
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs/sweep_cell_logmu_theta_joint_xenium_%A_%a.out
#SBATCH --array=0-8
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# 2026-09-17: cell log_mu x theta joint grid vs. the two Xenium real refs
# (lung_cancer, non_diseased_lung), mirroring sweep_spot_logmu_theta_joint_probe.sh's
# methodology exactly (composite = sum of |ln(sim/real ratio)| over theta_hat,
# total_counts_median, matrix_zero_frac; see summary_cell_logmu_theta_joint.sh).
#
# Why this exists: FIGURE2_METHODOLOGY.md's own "Known open gaps" section
# documents that cell's current config (log_mu=-2.3, theta=0.40, picked
# 2026-08-27) was a single manual adjustment targeting theta_hat/ARI, NOT a
# joint grid search like spot got (spot has 4 dedicated sweep scripts + a
# documented composite-loss decision). Cell has exactly one prior sweep
# script anywhere in the codebase (run_sweep_cell_log_mu-2.5_theta0.5.sh, a
# single point that doesn't even match the final config) and an explicitly
# deferred total-counts gap ("Accepted for the shared dataset; revisit if
# Figure 2 cell total-counts fidelity becomes priority"). Now that cell's
# smaller_sphere config is about to become Figure 3's strong-mix reference
# too, this sweep answers whether log_mu=-2.3/theta=0.40 is actually
# near-best or whether a real search finds something better.
#
# Grid: log_mu in {-2.5, -2.3, -2.0} x theta in {0.30, 0.40, 0.50} -- brackets
# the current config on both axes. theta_jitter held at 0.15 (already
# validated 2026-08-25 as artifact-free for cell, "only move the axis under
# test" discipline matching spot's sweeps), batch_sigma held at 1.5 (cell's
# canonical Figure 2/3 value -- FIGURE2_METHODOLOGY.md's cited theta_hat~0.083
# is the WITH-batch number, so the sweep must match that to be comparable).
# sphere_r_um=2050 / n_cells=24207 fixed throughout -- the smaller_sphere disc
# (see project_figure2_smaller_larger_sphere memory), NOT cell's old native
# 6000/600000 -- this sweep's whole point is to pick dispersion params for
# the diameter-matched config, so it must already run at that geometry.
#
# Real refs: sim_paper/data/real_data_qc/{lung_cancer,non_diseased_lung}/
# <label>_qc.h5ad (Xenium, 392-gene targeted panel) -- same refs
# code/count_distribution/smaller_sphere/cell_vs_*.sh already use.
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

LOG_MUS=(-2.5 -2.5 -2.5 -2.3 -2.3 -2.3 -2.0 -2.0 -2.0)
THETAS=(0.30  0.40  0.50 0.30  0.40  0.50 0.30  0.40  0.50)
LOG_MU="${LOG_MUS[${SLURM_ARRAY_TASK_ID:-0}]}"
THETA="${THETAS[${SLURM_ARRAY_TASK_ID:-0}]}"
MODALITY="cell"
TAG="cell_joint_logmu_${LOG_MU}_theta_${THETA}"

SIM_RAW="sim_paper/data/noisy/${TAG}/simulation_${MODALITY}_z.h5ad"
SIM_QC="sim_paper/data/noisy/${TAG}/simulation_${MODALITY}_z_qc.h5ad"
OUT_ROOT="sim_paper/data/count_distribution/sweeps/cell_joint_logmu_theta_xenium"

echo "[config] log_mu=${LOG_MU} theta=${THETA} tag=${TAG}"

if [[ ! -f "${SIM_RAW}" ]]; then
    echo "[generate] ${TAG} not found, generating"
    "${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py \
        --modality "${MODALITY}" \
        --sphere-r-um 2050 --n-cells 24207 \
        --base-gene-lognormal "${LOG_MU}" 0.7 \
        --theta "${THETA}" --theta-jitter 0.15 \
        --batch-sigma 1.5 \
        --out-tag "${TAG}"
fi

if [[ ! -f "${SIM_QC}" ]]; then
    echo "[qc] ${SIM_QC} not found, running 00_qc_filter.py"
    "${PYTHON_BIN}" sim_paper/code/clustering/00_qc_filter.py --modality "${MODALITY}" --packing-tag "${TAG}"
fi

for real in lung_cancer non_diseased_lung; do
    echo "[compare] ${TAG} vs ${real} (hvg-matched)"
    "${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
        --modality "${MODALITY}" --input "${SIM_QC}" \
        --compare-input "sim_paper/data/real_data_qc/${real}/${real}_qc.h5ad" \
        --compare-label "${real}" \
        --match-panel-size \
        --output-dir "${OUT_ROOT}/${TAG}_hvg_matched_vs_${real}"
done

echo "[done] ${TAG}"
