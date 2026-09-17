#!/bin/bash
#SBATCH --job-name=sweep_spot_logmu_probe
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs/sweep_spot_logmu_probe_%A_%a.out
#SBATCH --array=0-5
#SBATCH --time=08:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# spot log_mu bracketing sweep vs. the two 18k CytAssist probe-based real
# refs (lymph_node_visium, tonsil_visium), replacing breast_cancer_visium
# (36k, non-probe WTA -- dropped per user decision 2026-09-04) as spot's
# tuning target. Motivated by spot_vs_lymph_node_visium.sh's mismatch at the
# incumbent log_mu=-2.5 (tuned against breast_cancer_visium): qc_and_hvg_matched
# theta_hat ratio flipped to 1.87 (sim more dispersed than real, opposite of
# the breast_cancer_visium fit), qc_filtered total_counts ratio 0.05.
#
# batch_sigma is held FIXED at 0.3 throughout (explicit user constraint --
# don't disturb the Figure 3 Harmony before/after UMAP demo, which was tuned
# at this value). Only base_gene_lognormal's log_mu is swept -- per user
# decision (AskUserQuestion, 2026-09-04): log_mu only, not a log_mu x theta
# grid. theta/theta_jitter/marker_foldchange/sphere_r_um (packing_pf0p04)
# all stay at spot's current defaults.
#
# Uses qc_filtered mode only (bracketing pass, not a finalized Figure 2
# panel) -- same rationale as sweep_bin16um_logmu.sh's first pass. Winner
# picked via composite log-distance score on theta_hat/total_counts_median/
# matrix_zero_frac, SUMMED ACROSS BOTH real refs (user: "optimize for both
# simultaneously"), not picked independently per-ref.
#
# Prereq: real_data_qc/run_visium_tonsil_qc.sh (writes
# data/real_data_qc/tonsil_visium/tonsil_visium_qc.h5ad). Submit with:
#   sbatch --dependency=afterok:<tonsil_qc_job_id> sweep_spot_logmu_probe.sh
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

LOG_MUS=(-2.5 -2.0 -1.5 -1.0 -0.5 0.0)
LOG_MU="${LOG_MUS[${SLURM_ARRAY_TASK_ID:-0}]}"
TAG="spot_logmu_sweep_${LOG_MU}"
MODALITY="spot"

SIM_RAW="sim_paper/data/noisy/${TAG}/simulation_${MODALITY}_z.h5ad"
SIM_QC="sim_paper/data/noisy/${TAG}/simulation_${MODALITY}_z_qc.h5ad"
OUT_ROOT="sim_paper/data/count_distribution/sweeps/spot_logmu_sweep_probe"

echo "[config] log_mu=${LOG_MU} tag=${TAG}"

if [[ ! -f "${SIM_RAW}" ]]; then
    echo "[generate] ${TAG} not found, generating"
    "${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py \
        --modality "${MODALITY}" \
        --sphere-r-um 2050 \
        --base-gene-lognormal "${LOG_MU}" 0.7 \
        --batch-sigma 0.3 \
        --out-tag "${TAG}"
fi

if [[ ! -f "${SIM_QC}" ]]; then
    echo "[qc] ${SIM_QC} not found, running 00_qc_filter.py"
    "${PYTHON_BIN}" sim_paper/code/clustering/00_qc_filter.py --modality "${MODALITY}" --packing-tag "${TAG}"
fi

for real in tonsil_visium lymph_node_visium; do
    echo "[compare] ${TAG} vs ${real} (qc_filtered)"
    "${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
        --modality "${MODALITY}" --input "${SIM_QC}" \
        --compare-input "sim_paper/data/real_data_qc/${real}/${real}_qc.h5ad" \
        --compare-label "${real}" \
        --output-dir "${OUT_ROOT}/${TAG}_qc_filtered_vs_${real}"
done

echo "[done] ${TAG}"
