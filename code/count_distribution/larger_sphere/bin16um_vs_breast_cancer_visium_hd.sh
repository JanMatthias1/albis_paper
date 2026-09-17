#!/bin/bash
#SBATCH --job-name=fig2_larger_bin16um_vs_breast_cancer_visium_hd
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/larger_sphere/logs/fig2_larger_bin16um_vs_breast_cancer_visium_hd_%j.out
#SBATCH --time=24:00:00
#SBATCH --mem=400G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# larger_sphere counterpart of
# smaller_sphere/bin16um_vs_breast_cancer_visium_hd.sh -- resubmit of
# sweep_bin16um_sphere_r6000.sh's n_cells=15043310 point (the one that
# TIMED OUT at 64G/8h on 2026-09-16/17, MaxRSS still only 3.1G at
# cancellation -- i.e. it never even got close to finishing generation, not
# a near-miss). This version uses a fresh SIM_TAG (not the old
# bin16um_sphere_r6000_ncells15043310 tag, which left an empty directory
# from the failed attempt) and a much larger 400G/24h budget, sized off the
# lower two sweep points that DID complete: 3.0M cells -> ~60GB/38min,
# 5.6M cells -> ~97GB/1h43m (both exceeding their own 64G request without
# being killed, and scaling worse than linearly with n_cells).
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/larger_sphere/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

SIM_TAG="packing_pf0p04_bin16um_log_mu_-2.5_bsigma07_r6000"
MODALITY="bin"
REAL_LABEL="breast_cancer_visium_hd_16um"
SIM_RAW="sim_paper/data/figure_2/larger_sphere/data/${SIM_TAG}/simulation_${MODALITY}_z.h5ad"
SIM_QC="sim_paper/data/figure_2/larger_sphere/data/${SIM_TAG}/simulation_${MODALITY}_z_qc.h5ad"
NOISY_DIR="sim_paper/data/noisy/${SIM_TAG}"
REAL_INPUT="sim_paper/data/real_data_qc/${REAL_LABEL}/breast_cancer_visium_hd_qc.h5ad"
OUT_ROOT="sim_paper/data/figure_2/larger_sphere/plots/${MODALITY}_vs_${REAL_LABEL}"

if [[ ! -f "${SIM_RAW}" ]]; then
    echo "[generate] ${SIM_TAG} not present under figure_2/larger_sphere/data/, building"
    if [[ ! -f "${NOISY_DIR}/simulation_${MODALITY}_z.h5ad" ]]; then
        "${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py \
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
        "${PYTHON_BIN}" sim_paper/code/clustering/00_qc_filter.py \
            --modality "${MODALITY}" --packing-tag "${SIM_TAG}"
    fi
    mkdir -p sim_paper/data/figure_2/larger_sphere/data
    mv "${NOISY_DIR}" "sim_paper/data/figure_2/larger_sphere/data/${SIM_TAG}"
fi

if [[ ! -f "${SIM_QC}" ]]; then
    echo "[qc] ${SIM_QC} not found, running 00_qc_filter.py"
    "${PYTHON_BIN}" sim_paper/code/clustering/00_qc_filter.py \
        --modality "${MODALITY}" --input "${SIM_RAW}" --output "${SIM_QC}"
fi

"${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
    --modality "${MODALITY}" --input "${SIM_RAW}" \
    --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
    --output-dir "${OUT_ROOT}/full_panel"

"${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
    --modality "${MODALITY}" --input "${SIM_RAW}" \
    --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
    --match-panel-size \
    --output-dir "${OUT_ROOT}/hvg_matched"

"${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
    --modality "${MODALITY}" --input "${SIM_QC}" \
    --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
    --output-dir "${OUT_ROOT}/qc_filtered"

"${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
    --modality "${MODALITY}" --input "${SIM_QC}" \
    --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
    --match-panel-size \
    --output-dir "${OUT_ROOT}/qc_and_hvg_matched"

echo "[done] All 4 modes written under ${OUT_ROOT}"
