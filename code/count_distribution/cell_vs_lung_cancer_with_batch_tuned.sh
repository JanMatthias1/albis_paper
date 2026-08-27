#!/bin/bash
#SBATCH --job-name=fig2_cell_vs_lung_cancer_with_batch_tuned
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/logs/fig2_cell_vs_lung_cancer_with_batch_tuned_%j.out
#SBATCH --time=06:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Same pairing as cell_vs_lung_cancer_with_batch.sh, but reading from
# log_mu_-2.5_theta_0.25_jitter0.15_bsigma022 instead of the plain tag -- see
# cell_vs_non_diseased_lung_with_batch_tuned.sh for the full rationale (same
# sim config, different real reference). Writes to
# figure_2_with_batch_tuned/.

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

SIM_TAG="log_mu_-2.5_theta_0.25_jitter0.15_bsigma022"
MODALITY="cell"
REAL_LABEL="lung_cancer"
SIM_RAW="sim_paper/data/figure_2/${SIM_TAG}/simulation_${MODALITY}_z.h5ad"
SIM_QC="sim_paper/data/figure_2/${SIM_TAG}/simulation_${MODALITY}_z_qc.h5ad"
REAL_INPUT="sim_paper/data/real_data_qc/${REAL_LABEL}/${REAL_LABEL}_qc.h5ad"
OUT_ROOT="sim_paper/data/count_distribution/figure_2_with_batch_tuned/${MODALITY}_vs_${REAL_LABEL}"

if [[ ! -f "${SIM_QC}" ]]; then
    if [[ ! -f "${SIM_RAW}" ]]; then
        echo "[generate] ${SIM_TAG} not found, generating"
        "${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py \
            --modality "${MODALITY}" \
            --base-gene-lognormal -2.5 0.7 \
            --theta 0.25 \
            --theta-jitter 0.15 \
            --batch-sigma 0.22 \
            --out-tag "${SIM_TAG}"
        echo "[qc] ${SIM_TAG}"
        "${PYTHON_BIN}" sim_paper/code/clustering/00_qc_filter.py --modality "${MODALITY}" --packing-tag "${SIM_TAG}"
        echo "[move] data/noisy/${SIM_TAG} -> data/figure_2/${SIM_TAG}"
        mkdir -p sim_paper/data/figure_2
        mv "sim_paper/data/noisy/${SIM_TAG}" "sim_paper/data/figure_2/${SIM_TAG}"
    else
        echo "[qc] ${SIM_QC} not found but raw exists, running 00_qc_filter.py directly against figure_2/"
        "${PYTHON_BIN}" sim_paper/code/clustering/00_qc_filter.py --modality "${MODALITY}" \
            --input "${SIM_RAW}" --output "${SIM_QC}"
    fi
fi

"${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
    --modality "${MODALITY}" --input "${SIM_RAW}" --slice-id 5 --use-batch-effect \
    --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
    --output-dir "${OUT_ROOT}/full_panel"

"${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
    --modality "${MODALITY}" --input "${SIM_RAW}" --slice-id 5 --use-batch-effect \
    --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
    --match-panel-size \
    --output-dir "${OUT_ROOT}/hvg_matched"

"${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
    --modality "${MODALITY}" --input "${SIM_QC}" --slice-id 5 --use-batch-effect \
    --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
    --output-dir "${OUT_ROOT}/qc_filtered"

"${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
    --modality "${MODALITY}" --input "${SIM_QC}" --slice-id 5 --use-batch-effect \
    --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
    --match-panel-size \
    --output-dir "${OUT_ROOT}/qc_and_hvg_matched"

echo "[done] All 4 modes written under ${OUT_ROOT}"
