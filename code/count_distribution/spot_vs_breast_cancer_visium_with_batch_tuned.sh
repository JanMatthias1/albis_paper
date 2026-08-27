#!/bin/bash
#SBATCH --job-name=fig2_spot_vs_breast_cancer_visium_with_batch_tuned
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/logs/fig2_spot_vs_breast_cancer_visium_with_batch_tuned_%j.out
#SBATCH --time=08:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Same pairing as spot_vs_breast_cancer_visium_with_batch.sh, but reading
# from packing_pf0p04_log_mu_-2.5_bsigma015 (spot's Figure-3-decided
# batch_sigma=0.15) instead of the flat-0.22 tag -- see
# bin_vs_breast_cancer_visium_hd_with_batch_tuned.sh for the full rationale.
# Fully self-contained: generates the sim data (if not already present),
# QC-filters it, then runs all 4 count_distribution.py modes. Writes to
# figure_2_with_batch_tuned/.

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

SIM_TAG="packing_pf0p04_log_mu_-2.5_bsigma015"
MODALITY="spot"
REAL_LABEL="breast_cancer_visium"
SIM_RAW="sim_paper/data/figure_2/${SIM_TAG}/simulation_${MODALITY}_z.h5ad"
SIM_QC="sim_paper/data/figure_2/${SIM_TAG}/simulation_${MODALITY}_z_qc.h5ad"
REAL_INPUT="sim_paper/data/real_data_qc/${REAL_LABEL}/${REAL_LABEL}_qc.h5ad"
OUT_ROOT="sim_paper/data/count_distribution/figure_2_with_batch_tuned/${MODALITY}_vs_${REAL_LABEL}"

if [[ ! -f "${SIM_QC}" ]]; then
    if [[ ! -f "${SIM_RAW}" ]]; then
        echo "[generate] ${SIM_TAG} not found, generating"
        "${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py \
            --modality "${MODALITY}" \
            --sphere-r-um 2050 \
            --base-gene-lognormal -2.5 0.7 \
            --batch-sigma 0.15 \
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
