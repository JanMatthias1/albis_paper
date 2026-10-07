#!/bin/bash
#SBATCH --job-name=fig2_bin8um_same_tissue_as_16um
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/albis_paper/code/count_distribution/smaller_sphere/logs/fig2_bin8um_same_tissue_as_16um_%j.out
#SBATCH --time=08:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Figure 2, bin 8 µm: simulated 8 µm bins vs Visium HD human breast cancer and
# human pancreas (8 µm; QC by real_data_qc/), one run per reference.
# The simulation is the bin16um tissue at 8 µm: the bin16um settings with only
# --bin-size-um 8. The simulator bins the same molecules and draws batch factors
# independently of bin size, so these bins are an exact 2 × 2 subdivision of the
# 16 µm bins. QC modes only (most bins are off-tissue before QC).
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/albis_paper/code/count_distribution/smaller_sphere/logs
source /dcs04/hicks/data/Jan/sim_project/albis_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

SIM_TAG="packing_pf0p04_bin8um_from16umcfg_log_mu_-2.5_jitter0.6_bsigma07"
MODALITY="bin"
SIM_RAW="albis_paper/data/figure_2/smaller_sphere/data/${SIM_TAG}/simulation_${MODALITY}_z.h5ad"
SIM_QC="albis_paper/data/figure_2/smaller_sphere/data/${SIM_TAG}/simulation_${MODALITY}_z_qc.h5ad"
NOISY_DIR="albis_paper/data/noisy/${SIM_TAG}"

if [[ ! -f "${SIM_RAW}" ]]; then
    echo "[generate] ${SIM_TAG} not present under figure_2/, building"
    if [[ ! -f "${NOISY_DIR}/simulation_${MODALITY}_z.h5ad" ]]; then
        "${PYTHON_BIN}" albis_paper/code/data/generate_simulation_noisy.py \
            --modality "${MODALITY}" \
            --sphere-r-um 2050 \
            --bin-size-um 8 \
            --base-gene-lognormal -2.5 0.7 \
            --theta 2.0 \
            --theta-jitter 0.6 \
            --batch-sigma 0.7 \
            --sync-unaligned-seed \
            --out-tag "${SIM_TAG}"
    fi
    if [[ ! -f "${NOISY_DIR}/simulation_${MODALITY}_z_qc.h5ad" ]]; then
        "${PYTHON_BIN}" albis_paper/code/clustering/step00_qc_filter.py \
            --modality "${MODALITY}" --packing-tag "${SIM_TAG}"
    fi
    mkdir -p albis_paper/data/figure_2/smaller_sphere/data
    mv "${NOISY_DIR}" "albis_paper/data/figure_2/smaller_sphere/data/${SIM_TAG}"
fi

if [[ ! -f "${SIM_QC}" ]]; then
    echo "[qc] ${SIM_QC} not found, running step00_qc_filter.py"
    "${PYTHON_BIN}" albis_paper/code/clustering/step00_qc_filter.py \
        --modality "${MODALITY}" --input "${SIM_RAW}" --output "${SIM_QC}"
fi

for REAL_LABEL in breast_cancer_visium_hd human_pancreas_visium_hd; do
    REAL_INPUT="albis_paper/data/real_data_qc/${REAL_LABEL}/${REAL_LABEL}_qc.h5ad"
    OUT_ROOT="albis_paper/data/figure_2/smaller_sphere/plots/bin8um_same_tissue_as_16um_vs_${REAL_LABEL}"

    "${PYTHON_BIN}" albis_paper/code/count_distribution/count_distribution.py \
        --modality "${MODALITY}" --input "${SIM_QC}" \
        --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
        --output-dir "${OUT_ROOT}/qc_filtered"

    "${PYTHON_BIN}" albis_paper/code/count_distribution/count_distribution.py \
        --modality "${MODALITY}" --input "${SIM_QC}" \
        --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
        --match-panel-size \
        --output-dir "${OUT_ROOT}/qc_and_hvg_matched"

    echo "[done] ${OUT_ROOT}"
done
