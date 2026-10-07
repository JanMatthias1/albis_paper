#!/bin/bash
#SBATCH --job-name=fig2_larger_bin16um_vs_human_pancreas_visium_hd
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/albis_paper/code/count_distribution/larger_sphere/logs/fig2_larger_bin16um_vs_human_pancreas_visium_hd_%j.out
#SBATCH --time=24:00:00
#SBATCH --mem=400G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Larger-sphere version of smaller_sphere/bin16um_vs_human_pancreas_visium_hd.sh (supplementary, not part of
# Figure 2): the same comparison on tissue at r = 6000 µm with the smaller-sphere 3D packing (~15.0 M cells); generation needs several hundred GB of memory.
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/albis_paper/code/count_distribution/larger_sphere/logs
source /dcs04/hicks/data/Jan/sim_project/albis_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

SIM_TAG="packing_pf0p04_bin16um_log_mu_-2.5_bsigma07_r6000"
MODALITY="bin"
REAL_LABEL="human_pancreas_visium_hd_16um"
SIM_RAW="albis_paper/data/figure_2/larger_sphere/data/${SIM_TAG}/simulation_${MODALITY}_z.h5ad"
SIM_QC="albis_paper/data/figure_2/larger_sphere/data/${SIM_TAG}/simulation_${MODALITY}_z_qc.h5ad"
NOISY_DIR="albis_paper/data/noisy/${SIM_TAG}"
REAL_INPUT="albis_paper/data/real_data_qc/${REAL_LABEL}/human_pancreas_visium_hd_qc.h5ad"
OUT_ROOT="albis_paper/data/figure_2/larger_sphere/plots/${MODALITY}_vs_${REAL_LABEL}"

if [[ ! -f "${SIM_RAW}" ]]; then
    echo "[generate] ${SIM_TAG} not present under figure_2/larger_sphere/data/, building"
    if [[ ! -f "${NOISY_DIR}/simulation_${MODALITY}_z.h5ad" ]]; then
        "${PYTHON_BIN}" albis_paper/code/data/generate_simulation_noisy.py \
            --modality "${MODALITY}" \
            --sphere-r-um 6000 \
            --n-cells 15043310 \
            --bin-size-um 16 \
            --base-gene-lognormal -2.5 0.7 \
            --batch-sigma 0.7 \
            --sync-unaligned-seed \
            --out-tag "${SIM_TAG}"
    fi
    if [[ ! -f "${NOISY_DIR}/simulation_${MODALITY}_z_qc.h5ad" ]]; then
        "${PYTHON_BIN}" albis_paper/code/clustering/step00_qc_filter.py \
            --modality "${MODALITY}" --packing-tag "${SIM_TAG}"
    fi
    mkdir -p albis_paper/data/figure_2/larger_sphere/data
    mv "${NOISY_DIR}" "albis_paper/data/figure_2/larger_sphere/data/${SIM_TAG}"
fi

if [[ ! -f "${SIM_QC}" ]]; then
    echo "[qc] ${SIM_QC} not found, running step00_qc_filter.py"
    "${PYTHON_BIN}" albis_paper/code/clustering/step00_qc_filter.py \
        --modality "${MODALITY}" --input "${SIM_RAW}" --output "${SIM_QC}"
fi

"${PYTHON_BIN}" albis_paper/code/count_distribution/count_distribution.py \
    --modality "${MODALITY}" --input "${SIM_RAW}" \
    --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
    --output-dir "${OUT_ROOT}/full_panel"

"${PYTHON_BIN}" albis_paper/code/count_distribution/count_distribution.py \
    --modality "${MODALITY}" --input "${SIM_RAW}" \
    --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
    --match-panel-size \
    --output-dir "${OUT_ROOT}/hvg_matched"

"${PYTHON_BIN}" albis_paper/code/count_distribution/count_distribution.py \
    --modality "${MODALITY}" --input "${SIM_QC}" \
    --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
    --output-dir "${OUT_ROOT}/qc_filtered"

"${PYTHON_BIN}" albis_paper/code/count_distribution/count_distribution.py \
    --modality "${MODALITY}" --input "${SIM_QC}" \
    --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
    --match-panel-size \
    --output-dir "${OUT_ROOT}/qc_and_hvg_matched"

echo "[done] All 4 modes written under ${OUT_ROOT}"
